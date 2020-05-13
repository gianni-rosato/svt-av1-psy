/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decode MT process related functions

/**************************************
 * Includes
 **************************************/

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecMemInit.h"

#include "EbObuParse.h"
#include "EbDecParseFrame.h"
#include "EbDecProcessFrame.h"
#include "EbDecLF.h"
#include "EbDecCdef.h"
#include "EbDecRestoration.h"

#include "EbDecBitstream.h"
#include "EbTime.h"

#include "EbDecInverseQuantize.h"
#include "EbLog.h"

#include "EbUtility.h"

#include <stdlib.h>

#ifdef _WIN32
extern uint8_t        num_groups;
extern GROUP_AFFINITY group_affinity;
extern EbBool         alternate_groups;
#elif defined(__linux__)
extern cpu_set_t      group_affinity;
#endif
void *dec_all_stage_kernel(void *input_ptr);
/*ToDo : Remove all these replications */
void eb_av1_loop_filter_frame_init(FrameHeader *frm_hdr, LoopFilterInfoN *lfi, int32_t plane_start,
                                   int32_t plane_end);
void dec_loop_filter_row(EbDecHandle *dec_handle_ptr,
                         EbPictureBufferDesc *recon_picture_buf,
                         LfCtxt *lf_ctxt,
                         uint32_t y_sb_index,
                         int32_t plane_start, int32_t plane_end);
void save_deblock_boundary_lines(uint8_t *src_buf, int32_t src_stride, int32_t src_width,
                                 int32_t src_height, const Av1Common *cm, int32_t plane,
                                 int32_t row, int32_t stripe, int32_t use_highbd, int32_t is_above,
                                 RestorationStripeBoundaries *boundaries);
void save_cdef_boundary_lines(uint8_t *src_buf, int32_t src_stride, int32_t src_width,
                              const Av1Common *cm, int32_t plane, int32_t row, int32_t stripe,
                              int32_t use_highbd, int32_t is_above,
                              RestorationStripeBoundaries *boundaries);

EbErrorType dec_dummy_ctor(DecMtNode *context_ptr, EbPtr object_init_data_ptr) {
    context_ptr->node_index = *(uint32_t *)object_init_data_ptr;

    return EB_ErrorNone;
}

#if MT_WAIT_PROFILE
void dec_timer_start(struct EbDecTimer *t) {
#if defined(_WIN32)
    QueryPerformanceCounter(&t->begin);
#else
    gettimeofday(&t->begin, NULL);
#endif
}

void dec_timer_mark(struct EbDecTimer *t) {
#if defined(_WIN32)
    QueryPerformanceCounter(&t->end);
#else
    gettimeofday(&t->end, NULL);
#endif
}

int64_t dec_timer_elapsed(struct EbDecTimer *t) {
#if defined(_WIN32)
    LARGE_INTEGER freq, diff;

    diff.QuadPart = t->end.QuadPart - t->begin.QuadPart;

    QueryPerformanceFrequency(&freq);
    return diff.QuadPart * 1000000 / freq.QuadPart;
#else
    struct timeval diff;

    timersub(&t->end, &t->begin, &diff);
    return ((int64_t)diff.tv_sec) * 1000000 + diff.tv_usec;
#endif
}

void dec_display_timer(char s[], struct EbDecTimer *timer,
    int th_cnt, FILE *fp)
{
    dec_timer_mark(timer);
    uint64_t dx_time = dec_timer_elapsed(timer);
    fprintf(fp, "\n%s Tid %d %d", s, th_cnt, dx_time);
    fflush(fp);
}
#endif

/* Row MT Simple Q functions */
/* Return the sb_row_to_process */
int32_t get_sb_row_to_process(DecMtRowInfo *sb_row_info) {

    int32_t     sb_row_to_process = -1;
    //lock mutex
    eb_block_on_mutex(sb_row_info->sbrow_mutex);

    //pick up a row and increment the sb row counter
    if (sb_row_info->sb_row_to_process !=
        sb_row_info->num_sb_rows) {
        sb_row_to_process = sb_row_info->sb_row_to_process;
        sb_row_info->sb_row_to_process++;
        }

    //unlock mutex
    eb_release_mutex(sb_row_info->sbrow_mutex);

    return sb_row_to_process;
}

EbErrorType dec_dummy_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    DecMtNode *obj;
    *object_dbl_ptr = NULL;
    EB_NEW(obj, dec_dummy_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

/************************************
* System Resource Managers & Fifos
************************************/
EbErrorType dec_system_resource_init(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info) {
    EbErrorType     return_error = EB_ErrorNone;
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    memory_map_start_address = svt_dec_memory_map;
    memset(&dec_mt_frame_data->prev_frame_info, 0, sizeof(PrevFrameMtCheck));

    int32_t num_tiles = tiles_info->tile_cols * tiles_info->tile_rows;

    assert(dec_handle_ptr->dec_config.threads > 1);
#if MT_WAIT_PROFILE
    dec_mt_frame_data->fp = fopen("profile.txt", "w"); // stdout;
#endif
    /************************************
    * System Resource Managers & Fifos
    ************************************/

    /* Motion Filed Projection*/
    dec_mt_frame_data->motion_proj_info.num_motion_proj_rows = -1;
    EB_CREATE_MUTEX(dec_mt_frame_data->motion_proj_info.motion_proj_mutex);

    int32_t  sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    uint32_t picture_height_in_sb = (dec_handle_ptr->frame_header.
        frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    /************************************
    * Contexts
    ************************************/
    /* Parse */
    DecMtRowInfo *parse_tile_info = &dec_mt_frame_data->parse_tile_info;

    EB_CREATE_MUTEX(parse_tile_info->sbrow_mutex);
    parse_tile_info->num_sb_rows        = num_tiles;
    parse_tile_info->sb_row_to_process  = 0;

    /* Recon */
    EB_MALLOC_DEC(uint32_t *,
                  dec_mt_frame_data->sb_recon_row_map,
                  picture_height_in_sb * tiles_info->tile_cols * sizeof(uint32_t),
                  EB_N_PTR);

    DecMtRowInfo *recon_tile_info = &dec_mt_frame_data->recon_tile_info;
    EB_CREATE_MUTEX(recon_tile_info->sbrow_mutex);

    recon_tile_info->num_sb_rows        = num_tiles;
    recon_tile_info->sb_row_to_process  = 0;
    /* recon top right sync */
    {
        int32_t tiles_ctr;

        EB_CREATE_MUTEX(dec_mt_frame_data->tile_switch_mutex);

        EB_MALLOC_DEC(DecMtParseReconTileInfo *,
                      dec_mt_frame_data->parse_recon_tile_info_array,
                      num_tiles * sizeof(DecMtParseReconTileInfo),
                      EB_N_PTR);

        for (tiles_ctr = 0; tiles_ctr < num_tiles; tiles_ctr++) {
            int32_t   tile_row = tiles_ctr / tiles_info->tile_cols;
            int32_t   tile_col = tiles_ctr % tiles_info->tile_cols;
            int32_t   tile_num_sb_rows;
            TileInfo *tile_info;

            tile_info = &dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].tile_info;

            /* init tile info */
            svt_tile_init(tile_info, &dec_handle_ptr->frame_header, tile_row, tile_col);

            tile_num_sb_rows = ((((tile_info->mi_row_end - 1) << MI_SIZE_LOG2) >>
                                 dec_handle_ptr->seq_header.sb_size_log2) -
                                ((tile_info->mi_row_start << MI_SIZE_LOG2) >>
                                 dec_handle_ptr->seq_header.sb_size_log2) +
                                1);

            dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].tile_num_sb_rows =
                tile_num_sb_rows;

            EB_MALLOC_DEC(
                uint32_t *,
                dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_recon_row_parsed,
                tile_num_sb_rows * sizeof(uint32_t),
                EB_N_PTR);

            EB_MALLOC_DEC(
                uint32_t *,
                dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_recon_completed_in_row,
                tile_num_sb_rows * sizeof(uint32_t),
                EB_N_PTR);

            EB_MALLOC_DEC(
                uint32_t *,
                dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_recon_row_started,
                tile_num_sb_rows * sizeof(uint32_t),
                EB_N_PTR);

            EB_CREATE_MUTEX(
                dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].tile_sbrow_mutex);

            /* SB row queue */
            //EB_NEW(dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].recon_tile_sbrow_resource_ptr,
            //    eb_system_resource_ctor,
            //    tile_num_sb_rows, /* object_total_count */
            //    1, /* producer procs cnt : 1 Q per cnt is created inside, so kept 1*/
            //    1, /* consumer prcos cnt : 1 Q per cnt is created inside, so kept 1*/
            //    &dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].recon_tile_sbrow_producer_fifo_ptr, /* producer_fifo */
            //    &dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].recon_tile_sbrow_consumer_fifo_ptr, /* consumer_fifo */
            //    EB_TRUE, /* Full Queue*/
            //    dec_dummy_creator,
            //    &node_idx,
            //    NULL);
        }
    }

    /* LF */
    EB_MALLOC_DEC(int32_t *,
                  dec_mt_frame_data->lf_frame_info.sb_lf_completed_in_row,
                  picture_height_in_sb * sizeof(int32_t),
                  EB_N_PTR);

    EB_MALLOC_DEC(uint32_t *,
                  dec_mt_frame_data->lf_row_map,
                  picture_height_in_sb * sizeof(uint32_t),
                  EB_N_PTR);

    DecMtRowInfo *lf_sb_row_info = &dec_mt_frame_data->lf_frame_info.lf_sb_row_info;

    EB_CREATE_MUTEX(lf_sb_row_info->sbrow_mutex);
    lf_sb_row_info->num_sb_rows = picture_height_in_sb;
    lf_sb_row_info->sb_row_to_process = 0;

    /* CDEF */
    const int32_t num_planes = av1_num_planes(&dec_handle_ptr->seq_header.color_config);
    uint32_t      mi_cols    = 2 * ((dec_handle_ptr->seq_header.max_frame_width + 7) >> 3);

    const int32_t nhfb = (mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nvfb =
        (dec_handle_ptr->seq_header.max_frame_height + (MI_SIZE_64X64 << MI_SIZE_LOG2) - 1) /
        (MI_SIZE_64X64 << MI_SIZE_LOG2);

    const int32_t stride                   = (mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
    dec_mt_frame_data->cdef_linebuf_stride = stride;

    /*ToDo: Linebuff memory we can allocate min(sb_rows , threads)*/
    /*Currently we r allocating for every (64x64 +1 )rows*/
    EB_MALLOC_DEC(
        uint16_t ***, dec_mt_frame_data->cdef_linebuf, (nvfb + 1) * sizeof(uint16_t **), EB_N_PTR);
    for (int32_t sb_row = 0; sb_row < (nvfb + 1); sb_row++) {
        uint16_t **p_linebuf;
        EB_MALLOC_DEC(uint16_t **,
                      dec_mt_frame_data->cdef_linebuf[sb_row],
                      num_planes * sizeof(uint16_t **),
                      EB_N_PTR);
        p_linebuf = dec_mt_frame_data->cdef_linebuf[sb_row];
        for (int32_t pli = 0; pli < num_planes; pli++) {
            EB_MALLOC_DEC(
                uint16_t *, p_linebuf[pli], sizeof(uint16_t) * CDEF_VBORDER * stride, EB_N_PTR);
        }
    }

    dec_mt_frame_data->cdef_map_stride = nhfb + 2;
    /*For fbr=0, previous row cdef points some junk memory, if we allocate memory only for nvfb 64x64 blocks,
    to avoid to pointing junck memory, we allocate nvfb+1 64x64 blocks*/
    EB_MALLOC_DEC(uint8_t *,
                  dec_mt_frame_data->row_cdef_map,
                  (nvfb + 1) * dec_mt_frame_data->cdef_map_stride * sizeof(uint8_t),
                  EB_N_PTR);
    memset(dec_mt_frame_data->row_cdef_map,
           1,
           (nvfb + 1) * dec_mt_frame_data->cdef_map_stride * sizeof(uint8_t));

    EB_MALLOC_DEC(uint32_t *,
                  dec_mt_frame_data->cdef_completed_in_row,
                  (nvfb + 2) * sizeof(uint32_t),
                  EB_N_PTR);

    memset(dec_mt_frame_data->cdef_completed_in_row,
           0,
           (nvfb + 2) * //Rem here nhbf+2 u replaced with nvfb + 2
               sizeof(uint32_t));

    EB_MALLOC_DEC(uint32_t *,
                  dec_mt_frame_data->cdef_completed_for_row_map,
                  picture_height_in_sb * sizeof(uint32_t),
                  EB_N_PTR);

    DecMtRowInfo *cdef_sb_row_info = &dec_mt_frame_data->cdef_sb_row_info;

    EB_CREATE_MUTEX(cdef_sb_row_info->sbrow_mutex);

    cdef_sb_row_info->num_sb_rows       = picture_height_in_sb;
    cdef_sb_row_info->sb_row_to_process = 0;
    /* LR */
    EB_MALLOC_DEC(int32_t *,
                  dec_mt_frame_data->sb_lr_completed_in_row,
                  picture_height_in_sb * sizeof(int32_t),
                  EB_N_PTR);

    EB_MALLOC_DEC(uint32_t *,
                  dec_mt_frame_data->lr_row_map,
                  picture_height_in_sb * sizeof(uint32_t),
                  EB_N_PTR);

    DecMtRowInfo *lr_sb_row_info = &dec_mt_frame_data->lr_sb_row_info;

    EB_CREATE_MUTEX(lr_sb_row_info->sbrow_mutex);

    lr_sb_row_info->num_sb_rows         = picture_height_in_sb;
    lr_sb_row_info->sb_row_to_process   = 0;

    dec_mt_frame_data->temp_mutex = eb_create_mutex();

    dec_mt_frame_data->start_motion_proj  = EB_FALSE;
    dec_mt_frame_data->start_parse_frame  = EB_FALSE;
    dec_mt_frame_data->start_decode_frame = EB_FALSE;
    dec_mt_frame_data->start_lf_frame     = EB_FALSE;
    dec_mt_frame_data->start_cdef_frame   = EB_FALSE;
    dec_mt_frame_data->start_lr_frame     = EB_FALSE;
    dec_mt_frame_data->num_threads_cdefed = 0;
    dec_mt_frame_data->num_threads_lred   = 0;

    /************************************
    * Thread Handles
    ************************************/

    /* Decode Library Threads */
    uint32_t num_lib_threads = (int32_t)dec_handle_ptr->dec_config.threads - 1;

    /* Use a scratch memory so that the memory allocated within
       init_dec_mod_ctxt reallocated when required */
    void **dec_mod_ctxt_arr = NULL;
    dec_mod_ctxt_arr = (void **)malloc(num_lib_threads * sizeof(DecModCtxt*));

    for (uint32_t i = 0; i < num_lib_threads; i++) {
        init_dec_mod_ctxt(dec_handle_ptr,
            &dec_mod_ctxt_arr[i]);
    }

    memory_map_end_address = svt_dec_memory_map;

    if (EB_FALSE == dec_handle_ptr->start_thread_process) {
        dec_mt_frame_data->end_flag           = EB_FALSE;
        dec_mt_frame_data->num_threads_exited = 0;

        if (num_lib_threads > 0) {
            DecThreadCtxt *thread_ctxt_pa;
            EB_MALLOC_DEC(
                DecThreadCtxt *, thread_ctxt_pa, num_lib_threads * sizeof(DecThreadCtxt), EB_N_PTR);
            dec_handle_ptr->thread_ctxt_pa = thread_ctxt_pa;
            EB_CREATE_SEMAPHORE(dec_handle_ptr->thread_semaphore, 0, 100000);

            for (uint32_t i = 0; i < num_lib_threads; i++) {
                thread_ctxt_pa[i].thread_cnt     = i + 1;
                thread_ctxt_pa[i].dec_handle_ptr = dec_handle_ptr;
                thread_ctxt_pa[i].dec_mod_ctxt = dec_mod_ctxt_arr[i];
                EB_CREATE_SEMAPHORE(thread_ctxt_pa[i].thread_semaphore,
                    0, 100000);
                int use_highbd = (dec_handle_ptr->seq_header.color_config.bit_depth > EB_8BIT ||
                    dec_handle_ptr->is_16bit_pipeline);
                EB_MALLOC_DEC(uint8_t *,
                              thread_ctxt_pa[i].dst,
                              (MAX_SB_SIZE + 8) * RESTORATION_PROC_UNIT_SIZE *
                                    sizeof(uint8_t) << use_highbd,
                              EB_N_PTR);
            }
            EB_CREATE_THREAD_ARRAY(dec_handle_ptr->decode_thread_handle_array,
                                   num_lib_threads,
                                   dec_all_stage_kernel,
                                   (void **)&thread_ctxt_pa);
        }
    } else {
        for (uint32_t i = 0; i < num_lib_threads; i++) {
            dec_handle_ptr->thread_ctxt_pa[i].dec_mod_ctxt =
                dec_mod_ctxt_arr[i];
        }
    }
    free(dec_mod_ctxt_arr);
    return return_error;
}

/* Scan through the Tiles to find Bitstream offsets */
void svt_av1_scan_tiles(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info, ObuHeader *obu_header,
                        Bitstrm *bs, uint32_t tg_start, uint32_t tg_end) {
    size_t           tile_size;
    MasterParseCtxt *master_parse_ctxt = (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;
    ParseTileData *  parse_tile_data   = master_parse_ctxt->parse_tile_data;

    for (uint32_t tile_num = tg_start; tile_num <= tg_end; tile_num++) {
        if (tile_num == tg_end)
            tile_size = obu_header->payload_size;
        else {
            tile_size = dec_get_bits_le(bs, tiles_info->tile_size_bytes) + 1;
            obu_header->payload_size -= (tiles_info->tile_size_bytes + tile_size);
        }
        PRINT_FRAME("tile_size", (tile_size));

        // Assign to ParseCtxt
        parse_tile_data[tile_num].data      = get_bitsteam_buf(bs);
        parse_tile_data[tile_num].data_end  = bs->buf_max;
        parse_tile_data[tile_num].tile_size = tile_size;

        dec_bits_init(bs, (get_bitsteam_buf(bs) + tile_size), obu_header->payload_size);
    }
}

void svt_av1_queue_parse_jobs(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    int32_t sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    uint32_t picture_height_in_sb =
        (dec_handle_ptr->frame_header.frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    EB_MEMSET(dec_mt_frame_data->sb_recon_row_map,
              0,
              picture_height_in_sb * tiles_info->tile_cols * sizeof(uint32_t));

    dec_mt_frame_data->parse_tile_info.sb_row_to_process = 0;
    dec_mt_frame_data->recon_tile_info.sb_row_to_process = 0;
    //dec_handle_ptr->start_thread_process = EB_TRUE;
}

EbErrorType parse_tile_job(EbDecHandle *dec_handle_ptr, int32_t tile_num) {
    EbErrorType status = EB_ErrorNone;

    TilesInfo *      tiles_info        = &dec_handle_ptr->frame_header.tiles_info;
    MasterParseCtxt *master_parse_ctxt = (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;
    ParseCtxt *      parse_ctxt        = &master_parse_ctxt->tile_parse_ctxt[tile_num];

    parse_ctxt->seq_header   = &dec_handle_ptr->seq_header;
    parse_ctxt->frame_header = &dec_handle_ptr->frame_header;

    parse_ctxt->parse_above_nbr4x4_ctxt = &master_parse_ctxt->parse_above_nbr4x4_ctxt[tile_num];
    parse_ctxt->parse_left_nbr4x4_ctxt  = &master_parse_ctxt->parse_left_nbr4x4_ctxt[tile_num];

    start_parse_tile(dec_handle_ptr, parse_ctxt, tiles_info, tile_num, 1);

    return status;
}

void parse_frame_tiles(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile EbBool *start_parse_frame = &dec_mt_frame_data->start_parse_frame;
#if MT_WAIT_PROFILE
    FILE *            fp = dec_mt_frame_data->fp;
    struct EbDecTimer timer;
    int               th_cnt = NULL == thread_ctxt ? 0 : thread_ctxt->thread_cnt;
    dec_timer_start(&timer);
#endif
    while (*start_parse_frame != EB_TRUE)
        eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle_ptr->thread_semaphore
                                                  : thread_ctxt->thread_semaphore);

#if MT_WAIT_PROFILE
    dec_display_timer("SPF", &timer, th_cnt, fp);
#endif
    int32_t tile_num;
    while (1) {
#if MT_WAIT_PROFILE
        dec_timer_start(&timer);
#endif
        tile_num = get_sb_row_to_process(&dec_mt_frame_data->parse_tile_info);
        if (-1 != tile_num) {
            dec_mt_frame_data->start_decode_frame = EB_TRUE;
            if (EB_ErrorNone != parse_tile_job(dec_handle_ptr, tile_num)) {
                SVT_LOG("\nParse Issue for Tile %d", tile_num);
                break;
            }
            eb_post_semaphore(dec_handle_ptr->thread_semaphore);
            for (uint32_t lib_thrd = 0;
                lib_thrd < dec_handle_ptr->dec_config.threads - 1;
                 lib_thrd++)
            {
                eb_post_semaphore(dec_handle_ptr->
                    thread_ctxt_pa[lib_thrd].thread_semaphore);
            }
        } else
            break;
    }
}

EbErrorType decode_tile_job(EbDecHandle *dec_handle_ptr, int32_t tile_num,
                            DecModCtxt *dec_mod_ctxt) {
    EbErrorType status     = EB_ErrorNone;
    TilesInfo * tiles_info = &dec_handle_ptr->frame_header.tiles_info;
    status                 = start_decode_tile(dec_handle_ptr, dec_mod_ctxt, tiles_info, tile_num);
    return status;
}

void decode_frame_tiles(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile EbBool *start_decode_frame = &dec_mt_frame_data->start_decode_frame;
#if MT_WAIT_PROFILE
    FILE *            fp = dec_mt_frame_data->fp;
    struct EbDecTimer timer;
    int               th_cnt = NULL == thread_ctxt ? 0 : thread_ctxt->thread_cnt;
    dec_timer_start(&timer);
#endif
    while (*start_decode_frame != EB_TRUE)
        eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle_ptr->thread_semaphore
                                                  : thread_ctxt->thread_semaphore);

#if MT_WAIT_PROFILE
    dec_display_timer("SDF", &timer, th_cnt, fp);
#endif

    int32_t tile_num;
    while (1) {
        DecModCtxt *dec_mod_ctxt = (DecModCtxt *)dec_handle_ptr->pv_dec_mod_ctxt;

#if MT_WAIT_PROFILE
        dec_timer_start(&timer);
#endif
        tile_num = get_sb_row_to_process(&dec_mt_frame_data->recon_tile_info);
#if MT_WAIT_PROFILE
        dec_display_timer("GFDT", &timer, th_cnt, fp);
#endif
        if (thread_ctxt != NULL) {
            dec_mod_ctxt = thread_ctxt->dec_mod_ctxt;

            /* TODO : Calling this function at a tile level is
                   excessive. Move this call to operate at a frame level.*/
            setup_segmentation_dequant(thread_ctxt->dec_mod_ctxt);
        }

        if (-1 != tile_num) {
            if (EB_ErrorNone != decode_tile_job(dec_handle_ptr,
                tile_num, dec_mod_ctxt))
            {
                SVT_LOG("\nDecode Issue for Tile %d", tile_num);
                break;
            }
        } else {
            int32_t max_rows_pend = -1;
            int32_t tile_idx, next_tile_idx = -1;
            int32_t num_tiles = dec_handle_ptr->frame_header.tiles_info.tile_cols *
                                dec_handle_ptr->frame_header.tiles_info.tile_rows;

            eb_block_on_mutex(dec_mt_frame_data->tile_switch_mutex);

            //logic for switching acroos tile in SB decode
            for (tile_idx = 0; tile_idx < num_tiles; tile_idx++) {
                DecMtParseReconTileInfo *parse_recon_tile_info_array =
                    &dec_mt_frame_data->parse_recon_tile_info_array[tile_idx];
                int32_t tile_num_sb_rows = parse_recon_tile_info_array->tile_num_sb_rows;

                //check for incompleted tile with maximum rows to be processed
                if (parse_recon_tile_info_array->sb_row_to_process !=
                    parse_recon_tile_info_array->tile_num_sb_rows) {
                    int32_t ctr;
                    int32_t curr_tile_rows_pend;

                    for (ctr = 0; ctr < tile_num_sb_rows; ctr++) {
                        if (0 == parse_recon_tile_info_array->sb_recon_row_started[ctr]) { break; }
                    }

                    curr_tile_rows_pend = tile_num_sb_rows - ctr;

                    //check for min completed rows in a tile
                    if (max_rows_pend < curr_tile_rows_pend) {
                        max_rows_pend = curr_tile_rows_pend;
                        next_tile_idx = tile_idx;
                    }
                }
            }

            eb_release_mutex(dec_mt_frame_data->tile_switch_mutex);

            //break from while 1 loop if no tile to be processed
            if (-1 == next_tile_idx) {
                break;
            } else {
                if (EB_ErrorNone != decode_tile_job(dec_handle_ptr, next_tile_idx, dec_mod_ctxt)) {
                    SVT_LOG("\nDecode Issue for Tile %d", next_tile_idx);
                    break;
                }
            }
        }
    }
    return;
}

void svt_av1_queue_lf_jobs(EbDecHandle *dec_handle_ptr) {
    DecMtlfFrameInfo *lf_frame_info =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data.lf_frame_info;
    int32_t sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    /* ToDo : picture_height_in_sb used many places. Reuse! */
    uint32_t picture_height_in_sb =
        (dec_handle_ptr->frame_header.frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    EB_MEMSET(dec_mt_frame_data->lf_row_map, 0, picture_height_in_sb * sizeof(uint32_t));

    memset(lf_frame_info->sb_lf_completed_in_row, -1, picture_height_in_sb * sizeof(int32_t));
    lf_frame_info->lf_info_init_done = EB_FALSE;
    lf_frame_info->lf_sb_row_info.sb_row_to_process = 0;
}

/* Store LF_boundary_line req for LR */
static INLINE void dec_save_lf_boundary_lines_sb_row(EbDecHandle *  dec_handle,
                                                     Av1PixelRect **tile_rect, int32_t sb_row,
                                                     uint8_t **src, int32_t *stride,
                                                     int32_t num_planes) {
    Av1Common *cm         = &dec_handle->cm;
    FrameSize *frame_size = &dec_handle->frame_header.frame_size;
    EbBool     sb_128     = dec_handle->seq_header.sb_size == BLOCK_128X128;
    int32_t    num64s     = sb_128 ? 1 : 0;
    const int use_highbd = (dec_handle->seq_header.color_config.bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);
    LrCtxt *   lr_ctxt    = (LrCtxt *)dec_handle->pv_lr_ctxt;
    int32_t frame_stripe /* 64 strip */, plane_height;
    for (int32_t p = 0; p < num_planes; ++p) {
        int32_t                      ss_x          = p ? cm->subsampling_x : 0;
        int32_t                      ss_y          = p ? cm->subsampling_y : 0;
        const int32_t                stripe_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;
        const int32_t                stripe_off    = RESTORATION_UNIT_OFFSET >> ss_y;
        RestorationStripeBoundaries *boundaries    = &lr_ctxt->boundaries[p];

        plane_height = ROUND_POWER_OF_TWO(cm->frm_size.frame_height, ss_y);

        int32_t src_width  = frame_size->frame_width >> ss_x;
        int32_t src_height = frame_size->frame_height >> ss_y;

        for (int32_t row_cnt = 0; row_cnt <= num64s; row_cnt++) {
            frame_stripe         = (sb_row << num64s) + row_cnt;
            const int32_t rel_y0 = AOMMAX(0, frame_stripe * stripe_height - stripe_off);
            const int32_t y0     = tile_rect[p]->top + rel_y0;
            if (y0 >= tile_rect[p]->bottom) break;

            const int32_t rel_y1 = (frame_stripe + 1) * stripe_height - stripe_off;
            const int32_t y1     = AOMMIN(tile_rect[p]->top + rel_y1, tile_rect[p]->bottom);

            int32_t use_deblock_above, use_deblock_below;
            // In this case, we should only use CDEF pixels at the top
            // and bottom of the frame as a whole; internal tile boundaries
            // can use deblocked pixels from adjacent tiles for context.
            use_deblock_above = (frame_stripe > 0);
            use_deblock_below = (y1 < plane_height);

            // Save deblocked context where needed.
            if (use_deblock_above) {
                save_deblock_boundary_lines(src[p],
                                            stride[p],
                                            src_width,
                                            src_height,
                                            cm,
                                            p,
                                            y0 - RESTORATION_CTX_VERT,
                                            frame_stripe,
                                            use_highbd,
                                            1,
                                            boundaries);
            }
            if (use_deblock_below) {
                save_deblock_boundary_lines(src[p],
                                            stride[p],
                                            src_width,
                                            src_height,
                                            cm,
                                            p,
                                            y1,
                                            frame_stripe,
                                            use_highbd,
                                            0,
                                            boundaries);
            }
        }
    }
}

/* Store CDEF_boundary_line req for LR */
static INLINE void dec_save_CDEF_boundary_lines_SB_row(
    EbDecHandle *  dec_handle, Av1PixelRect **tile_rect,
    int32_t sb_row, uint8_t **src, int32_t *stride, int32_t num_planes)
{
    Av1Common *     cm         = &dec_handle->cm;
    FrameSize *     frame_size = &dec_handle->frame_header.frame_size;
    const int use_highbd = (dec_handle->seq_header.color_config.bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);
    LrCtxt *        lr_ctxt    = (LrCtxt *)dec_handle->pv_lr_ctxt;
    int32_t         frame_stripe /* 64 strip */;
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    for (int32_t p = 0; p < num_planes; ++p) {
        int32_t         ss_x          = p ? cm->subsampling_x : 0;
        int32_t         ss_y          = p ? cm->subsampling_y : 0;
        const int32_t   stripe_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;
        const int32_t   stripe_off    = RESTORATION_UNIT_OFFSET >> ss_y;

        RestorationStripeBoundaries *boundaries = &lr_ctxt->boundaries[p];
        int32_t                      src_width  = frame_size->frame_width >> ss_x;

        frame_stripe = 0;
        if (sb_row == dec_mt_frame_data->sb_rows - 1)
            frame_stripe = frame_size->frame_height >> MIN_SB_SIZE_LOG2;

        const int32_t rel_y0 = AOMMAX(0, frame_stripe * stripe_height - stripe_off);
        const int32_t y0     = tile_rect[p]->top + rel_y0;

        const int32_t rel_y1 = (frame_stripe + 1) * stripe_height - stripe_off;
        int32_t y1  = AOMMIN(tile_rect[p]->top + rel_y1, tile_rect[p]->bottom);

        int32_t plane_height =
            ROUND_POWER_OF_TWO(cm->frm_size.frame_height, ss_y);

        // Save CDEF context where needed.
        if (frame_stripe == 0) {
            save_cdef_boundary_lines(
                src[p], stride[p], src_width, cm, p, y0,
                frame_stripe, use_highbd, 1, boundaries);
        }
        if (y1 >= plane_height) {
            save_cdef_boundary_lines(
                src[p], stride[p], src_width, cm, p, y1 - 1,
                frame_stripe, use_highbd, 0, boundaries);
        }
    }
}

/*Frame level function to trigger loop filter for each superblock*/
void dec_av1_loop_filter_frame_mt(EbDecHandle *dec_handle,
                                  EbPictureBufferDesc *recon_picture_buf,
                                  LfCtxt *lf_ctxt,
                                  int32_t plane_start,
                                  int32_t plane_end,
                                  DecThreadCtxt *thread_ctxt) {
    int32_t         sb_row;
    DecMtFrameData *dec_mt_frame_data1 =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    volatile EbBool *start_lf_frame = &dec_mt_frame_data1->start_lf_frame;
#if MT_WAIT_PROFILE
    FILE *            fp = dec_mt_frame_data1->fp;
    struct EbDecTimer timer;
    int               th_cnt = NULL == thread_ctxt ? 0 : thread_ctxt->thread_cnt;
    dec_timer_start(&timer);
#endif
    while (*start_lf_frame != EB_TRUE)
        eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle->thread_semaphore
                                                  : thread_ctxt->thread_semaphore);
#if MT_WAIT_PROFILE
    dec_display_timer("SLF", &timer, th_cnt, fp);
#endif

    FrameHeader *frm_hdr = &dec_handle->frame_header;

    lf_ctxt->delta_lf_stride = dec_handle->master_frame_buf.sb_cols * FRAME_LF_COUNT;
    frm_hdr->loop_filter_params.combine_vert_horz_lf = 1;

    DecMtlfFrameInfo *dec_mt_lf_frame_info = &dec_mt_frame_data1->lf_frame_info;
    //lock mutex
    eb_block_on_mutex(dec_mt_lf_frame_info->lf_sb_row_info.sbrow_mutex);
    //Check lf_info init done or not
    if (dec_mt_lf_frame_info->lf_info_init_done == EB_FALSE) {
        /*init hev threshold const vectors*/
        for (int lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
             memset(lf_ctxt->lf_info.lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);

        eb_av1_loop_filter_frame_init(frm_hdr,
                                      &lf_ctxt->lf_info,
                                      plane_start,
                                      plane_end);

        dec_mt_lf_frame_info->lf_info_init_done = EB_TRUE;
    }
    //unlock mutex
    eb_release_mutex(dec_mt_lf_frame_info->lf_sb_row_info.sbrow_mutex);

    set_lbd_lf_filter_tap_functions();
    set_hbd_lf_filter_tap_functions();

    /* For LF boundary store : Can optimize based on CDEF & LR flag */
    // Get the tile rectangle, with height rounded up to the next multiple of 8
    // luma pixels (only relevant for the bottom tile of the frame)
    Av1PixelRect         tile_rect[MAX_MB_PLANE];
    Av1PixelRect *       tile_rect_p[MAX_MB_PLANE];
    const int            num_planes  = av1_num_planes(&dec_handle->seq_header.color_config);
    EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    uint8_t *            src[MAX_MB_PLANE];
    int32_t              stride[MAX_MB_PLANE];
    for (int p = 0; p < num_planes; ++p) {
        int32_t is_uv  = p ? 1 : 0;
        int32_t sx     = is_uv ? dec_handle->cm.subsampling_x : 0;
        int32_t sy     = is_uv ? dec_handle->cm.subsampling_y : 0;
        tile_rect[p]   = whole_frame_rect(&dec_handle->cm.frm_size, sx, sy, is_uv);
        tile_rect_p[p] = &tile_rect[p];
        derive_blk_pointers(cur_pic_buf, p, 0, 0, (void *)&src[p], &stride[p], sx, sy);
    }

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    while (1) {
#if MT_WAIT_PROFILE
        dec_timer_start(&timer);
#endif
        sb_row = get_sb_row_to_process(&dec_mt_frame_data->lf_frame_info.
            lf_sb_row_info);

#if MT_WAIT_PROFILE
        dec_display_timer("GFLF", &timer, th_cnt, fp);
#endif
        if (-1 != sb_row) {
            TilesInfo *tiles_info = &dec_handle->frame_header.tiles_info;

            /* Ensure all TileRecon jobs are over for (row, row-1, row+1)      */
            /* row-1 : To ensure line buf copy with TopR sync if LF skips row  */
            /* This prevent issues across Tiles where recon sync is not ensured*/
            /* row+1 : This is for CDEF actually, should be moved to CDEF stage*/
            int32_t start_lf[3] = {0};
            int32_t row_index[3];
            row_index[0] = (sb_row)*tiles_info->tile_cols;
            row_index[1] = (sb_row - (sb_row == 0 ? 0 : 1)) * tiles_info->tile_cols;
            row_index[2] = (sb_row + (sb_row == (dec_mt_frame_data->sb_rows - 1) ? 0 : 1)) *
                           tiles_info->tile_cols;
#if MT_WAIT_PROFILE
            dec_timer_start(&timer);
#endif
            while ((!start_lf[0]) || (!start_lf[1]) || (!start_lf[2])) {
                start_lf[0] = 1;
                start_lf[1] = 1;
                start_lf[2] = 1;
                for (int i = 0; i < tiles_info->tile_cols; i++) {
                    start_lf[0] &= dec_mt_frame_data->sb_recon_row_map[row_index[0] + i];
                    start_lf[1] &= dec_mt_frame_data->sb_recon_row_map[row_index[1] + i];
                    start_lf[2] &= dec_mt_frame_data->sb_recon_row_map[row_index[2] + i];
                }
            }
#if MT_WAIT_PROFILE
            dec_display_timer("LFWR", &timer, th_cnt, fp);
#endif
            if (!dec_handle->frame_header.allow_intrabc) {
                if (dec_handle->frame_header.loop_filter_params.filter_level[0] ||
                    dec_handle->frame_header.loop_filter_params.filter_level[1]) {
                    dec_loop_filter_row(dec_handle,
                                        recon_picture_buf,
                                        lf_ctxt,
                                        sb_row,
                                        plane_start,
                                        plane_end);
                }
            }

            /* Store LR_save_boundary_lines at 64 lines : After LF         */
            /* Store Above 64 line always, for SB 128 store Middle 64 also */
            /* Bottom 64 won't be ready yet as next SB row Lf can modify it*/
            /* TO DO: Should be based on LR flag! */
            if (sb_row != 0) {
                dec_save_lf_boundary_lines_sb_row(
                    dec_handle, tile_rect_p, sb_row - 1, src, stride, num_planes);

                /* Update LF done map */
                dec_mt_frame_data1->lf_row_map[sb_row - 1] = 1;
            }
            if (sb_row == dec_mt_frame_data->sb_rows - 1) {
                dec_save_lf_boundary_lines_sb_row(
                    dec_handle, tile_rect_p, sb_row, src, stride, num_planes);

                /* Update LF done map */
                dec_mt_frame_data1->lf_row_map[sb_row] = 1;
            }
        } else
            break;
    }
}

void svt_av1_queue_cdef_jobs(EbDecHandle *dec_handle_ptr) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    int32_t  sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    uint32_t picture_height_in_sb =
        (dec_handle_ptr->frame_header.frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    const int32_t nvfb = (dec_handle_ptr->frame_header.mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;

    memset(dec_mt_frame_data->cdef_completed_in_row, 0, nvfb * sizeof(uint32_t));

    EB_MEMSET(
        dec_mt_frame_data->cdef_completed_for_row_map, 0, picture_height_in_sb * sizeof(uint32_t));

    dec_mt_frame_data->cdef_sb_row_info.sb_row_to_process = 0;
}

void svt_cdef_frame_mt(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt) {
    uint8_t *       curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t         curr_recon_stride[MAX_MB_PLANE];
    DecMtFrameData *dec_mt_frame_data1 =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile EbBool *start_cdef_frame = &dec_mt_frame_data1->start_cdef_frame;
#if MT_WAIT_PROFILE
    FILE *            fp = dec_mt_frame_data1->fp;
    struct EbDecTimer timer;
    int               th_cnt = NULL == thread_ctxt ? 0 : thread_ctxt->thread_cnt;
    dec_timer_start(&timer);
#endif
    while (*start_cdef_frame != EB_TRUE)
        eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle_ptr->thread_semaphore
                                                  : thread_ctxt->thread_semaphore);

#if MT_WAIT_PROFILE
    dec_display_timer("SCF", &timer, th_cnt, fp);
#endif
    EbPictureBufferDesc *recon_picture_ptr = dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf;
    const int32_t        num_planes = av1_num_planes(&dec_handle_ptr->seq_header.color_config);

    DECLARE_ALIGNED(16, uint16_t, src[CDEF_INBUF_SIZE]);
    uint16_t *colbuf[2 * 3];
    int32_t   mi_wide_l2[3];
    int32_t   mi_high_l2[3];

    Av1PixelRect  tile_rect[MAX_MB_PLANE];
    Av1PixelRect *tile_rect_p[MAX_MB_PLANE];

    EbBool no_ibc      = !dec_handle_ptr->frame_header.allow_intrabc;
    EbBool do_upscale  = no_ibc &&
        !av1_superres_unscaled(&dec_handle_ptr->frame_header.frame_size);
    LrParams *lr_param = dec_handle_ptr->frame_header.lr_params;
    EbBool    do_lr    = no_ibc &&
        (lr_param[AOM_PLANE_Y].frame_restoration_type != RESTORE_NONE ||
         lr_param[AOM_PLANE_U].frame_restoration_type != RESTORE_NONE ||
         lr_param[AOM_PLANE_V].frame_restoration_type != RESTORE_NONE);

    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t is_uv   = pli ? 1 : 0;
        int32_t sub_x   = !is_uv ?
            0 : dec_handle_ptr->seq_header.color_config.subsampling_x;
        int32_t sub_y   = !is_uv ?
            0 : dec_handle_ptr->seq_header.color_config.subsampling_y;
        mi_wide_l2[pli] = MI_SIZE_LOG2 - sub_x;
        mi_high_l2[pli] = MI_SIZE_LOG2 - sub_y;

        tile_rect[pli]   =
            whole_frame_rect(&dec_handle_ptr->cm.frm_size, sub_x, sub_y, is_uv);
        tile_rect_p[pli] = &tile_rect[pli];

        /*Deriveing  recon pict buffer ptr's*/
        derive_blk_pointers(recon_picture_ptr,
                            pli,
                            0,
                            0,
                            (void *)&curr_blk_recon_buf[pli],
                            &curr_recon_stride[pli],
                            sub_x,
                            sub_y);

        if (dec_handle_ptr->seq_header.sb_size == BLOCK_128X128) {
            /*For SB SIZE 128x128, we need two colbuf because, because we do cdef for
            each 64x64 in SB block in raster scan order,
            i.e for transversing across 0 - 3 64x64s in SB block*/
            for (int32_t i = 0; i < 4; i += 3)
                colbuf[pli + i] = (uint16_t *)eb_aom_malloc(
                    sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
                    CDEF_HBORDER);
        } else {
            colbuf[pli] = (uint16_t *)eb_aom_malloc(
                sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
                CDEF_HBORDER);
        }
    }

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    int32_t sb_row;

    while (1) {
#if MT_WAIT_PROFILE
        dec_timer_start(&timer);
#endif
        sb_row = get_sb_row_to_process(&dec_mt_frame_data->cdef_sb_row_info);
#if MT_WAIT_PROFILE
        dec_display_timer("GFCF", &timer, th_cnt, fp);
#endif
        if (-1 != sb_row) {
            /* Ensure all LF jobs are over for row_index (row / row+1) */
            int32_t offset = sb_row == dec_mt_frame_data->sb_rows - 1 ? 0 : 1;
#if MT_WAIT_PROFILE
            dec_timer_start(&timer);
#endif
            volatile int32_t *start_cdef = (volatile int32_t *)&dec_mt_frame_data
                                                ->lf_row_map[sb_row + offset];
            while (!*start_cdef)
                ;
            assert(*start_cdef == 1);
#if MT_WAIT_PROFILE
            dec_display_timer("CWLF", &timer, th_cnt, fp);
#endif
            FrameHeader *frame_header = &dec_handle_ptr->frame_header;
            if (!frame_header->allow_intrabc) {
                const int32_t do_cdef = !frame_header->coded_lossless &&
                                        (frame_header->cdef_params.cdef_bits ||
                                         frame_header->cdef_params.cdef_y_strength[0] ||
                                         frame_header->cdef_params.cdef_uv_strength[0]);
                if (do_cdef) {
                    /* SB Row Index */
                    int32_t sb_fbr = (int32_t)sb_row;

                    svt_cdef_sb_row_mt(dec_handle_ptr,
                                       mi_wide_l2,
                                       mi_high_l2,
                                       &colbuf[0],
                                       sb_fbr,
                                       &src[0],
                                       &curr_recon_stride[0],
                                       &curr_blk_recon_buf[0]);
                }
            }

            if (do_lr && !do_upscale) {
                // In this case, we should only use CDEF pixels at the top
                // and bottom of the frame as a whole; internal tile boundaries
                // can use deblocked pixels from adjacent tiles for context.
                if (sb_row == 0 || sb_row == dec_mt_frame_data->sb_rows - 1) {
                    dec_save_CDEF_boundary_lines_SB_row(dec_handle_ptr,
                                                        tile_rect_p,
                                                        sb_row,
                                                        curr_blk_recon_buf,
                                                        curr_recon_stride,
                                                        num_planes);
                }
            }
            /* Update CDEF done map */
            dec_mt_frame_data1->cdef_completed_for_row_map[sb_row] = 1;

        } else
            break;
    }
    if (dec_handle_ptr->seq_header.sb_size == BLOCK_128X128) {
        for (int32_t i = 0; i < 4; i += 3) {
            for (int32_t pli = 0; pli < num_planes; pli++) { eb_aom_free(colbuf[pli + i]); }
        }
    } else
        for (int32_t pli = 0; pli < num_planes; pli++) { eb_aom_free(colbuf[pli]); }

    eb_block_on_mutex(dec_mt_frame_data->temp_mutex);
    dec_mt_frame_data->num_threads_cdefed++;
    eb_release_mutex(dec_mt_frame_data->temp_mutex);
    if (do_upscale) {
        volatile uint32_t *num_threads_cdefed = &dec_mt_frame_data->num_threads_cdefed;
        while (*num_threads_cdefed != dec_handle_ptr->dec_config.threads)
            ;
    }
}

void svt_av1_queue_lr_jobs(EbDecHandle *dec_handle_ptr) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    int32_t  sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    uint32_t picture_height_in_sb =
        (dec_handle_ptr->frame_header.frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    EB_MEMSET(dec_mt_frame_data->lr_row_map, 0, picture_height_in_sb * sizeof(uint32_t));

    memset(dec_mt_frame_data->sb_lr_completed_in_row, -1, picture_height_in_sb * sizeof(int32_t));
    dec_mt_frame_data->lr_sb_row_info.sb_row_to_process = 0;
}

void pad_pre_lr(EbPictureBufferDesc *recon_picture_buf, int32_t sb_row, int32_t sb_size,
                int32_t num_rows, uint8_t **curr_blk_recon_buf, int *rec_stride,
                uint32_t frame_width, uint32_t frame_height, int sx, int sy) {
    /* Current SB row */
    uint32_t row = sb_row * sb_size;

    /* Plane buffers pointing to the current row */
    EbByte src_y  = curr_blk_recon_buf[AOM_PLANE_Y] + (row * rec_stride[AOM_PLANE_Y]);
    EbByte src_cb = curr_blk_recon_buf[AOM_PLANE_U] + ((row >> sy) * rec_stride[AOM_PLANE_U]);
    EbByte src_cr = curr_blk_recon_buf[AOM_PLANE_V] + ((row >> sy) * rec_stride[AOM_PLANE_V]);

    /* Directions in which padding needs to be done */
    PadDir flags = get_neighbour_flags(sb_row, 0, num_rows, 1);

    int32_t row_height = MIN(sb_size, (int32_t)(frame_height - row));
    pad_row(recon_picture_buf,
            src_y,
            src_cb,
            src_cr,
            frame_width,
            row_height,
            LR_PAD_SIDE,
            LR_PAD_SIDE,
            sx,
            sy,
            flags);
}

void pad_post_lr(EbPictureBufferDesc *recon_picture_buf, int32_t sb_row, int32_t sb_size,
                 int32_t num_rows, int *rec_stride, uint32_t pad_width, uint32_t pad_height,
                 uint32_t shift, uint32_t frame_width, uint32_t frame_height, int sx, int sy) {
    /* Do padding of Nth row after the completion of
       LR of (N+1) row to avoid possible corruption of
       recon buffers. */
    if (sb_row != 0 || num_rows == 1) {
        int processing_row = sb_row;
        if (num_rows != 1) processing_row--;
        int loop = 1;

        /* Pad last 2 rows when the last row's LR
           is being processed*/
        if (processing_row == num_rows - 2 && num_rows > 1) loop = 2;
        for (int i = 0; i < loop; i++, processing_row++) {
            uint32_t row   = processing_row * sb_size;
            EbByte   src_y = recon_picture_buf->buffer_y + (pad_width << shift) +
                           ((pad_height + row) * rec_stride[AOM_PLANE_Y]);
            EbByte src_cb = recon_picture_buf->buffer_cb + (pad_width >> sx << shift) +
                            (((pad_height + row) >> sy) * rec_stride[AOM_PLANE_U]);
            EbByte src_cr = recon_picture_buf->buffer_cr + (pad_width >> sx << shift) +
                            (((pad_height + row) >> sy) * rec_stride[AOM_PLANE_V]);

            int32_t row_height = MIN(sb_size, (int32_t)(frame_height - row));

            PadDir flags = get_neighbour_flags(processing_row, 0, num_rows, 1);

            pad_row(recon_picture_buf,
                    src_y,
                    src_cb,
                    src_cr,
                    frame_width,
                    row_height,
                    pad_width,
                    pad_height,
                    sx,
                    sy,
                    flags);
        }
    }
}

void dec_av1_loop_restoration_filter_frame_mt(
    EbDecHandle *dec_handle, DecThreadCtxt *thread_ctxt)
{
    uint8_t *    curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t      curr_recon_stride[MAX_MB_PLANE];

    Av1PixelRect tile_rect[MAX_MB_PLANE];
    Av1PixelRect *tile_rect_p[MAX_MB_PLANE];

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile EbBool *start_lr_frame = &dec_mt_frame_data->start_lr_frame;
    while (*start_lr_frame != EB_TRUE)
        eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle->thread_semaphore
                                                  : thread_ctxt->thread_semaphore);

    EbPictureBufferDesc *recon_picture_ptr = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    const int32_t        num_planes        = av1_num_planes(&dec_handle->seq_header.color_config);

    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t sub_x = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_x;
        int32_t sub_y = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_y;

        /*Deriveing  recon pict buffer ptr's*/
        derive_blk_pointers(recon_picture_ptr,
                            pli,
                            0,
                            0,
                            (void *)&curr_blk_recon_buf[pli],
                            &curr_recon_stride[pli],
                            sub_x,
                            sub_y);

        tile_rect[pli] =
            whole_frame_rect(&dec_handle->frame_header.frame_size, sub_x, sub_y, pli > 0);
        tile_rect_p[pli] = &tile_rect[pli];
    }
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;

    uint32_t frame_width  = dec_handle->frame_header.frame_size.superres_upscaled_width;
    uint32_t frame_height = dec_handle->frame_header.frame_size.frame_height;

    int sx = dec_handle->seq_header.color_config.subsampling_x;
    int sy = dec_handle->seq_header.color_config.subsampling_y;

    int32_t sb_size           = dec_handle->seq_header.use_128x128_superblock ? 128 : 64;
    int32_t sb_size_log2      = dec_handle->seq_header.sb_size_log2;
    int32_t sb_aligned_height = ALIGN_POWER_OF_TWO(frame_height, sb_size_log2);
    int32_t num_rows          = sb_aligned_height >> sb_size_log2;

    uint32_t pad_width  = recon_picture_buf->origin_x;
    uint32_t pad_height = recon_picture_buf->origin_y;

    int32_t shift = 0;
    if ((recon_picture_buf->bit_depth != EB_8BIT) ||
        recon_picture_buf->is_16bit_pipeline) shift = 1;

    int32_t recon_stride[MAX_MB_PLANE];
    recon_stride[AOM_PLANE_Y] = recon_picture_buf->stride_y << shift;
    recon_stride[AOM_PLANE_U] = recon_picture_buf->stride_cb << shift;
    recon_stride[AOM_PLANE_V] = recon_picture_buf->stride_cr << shift;

    int32_t sb_row;

    FrameHeader *frame_header = &dec_handle->frame_header;

    EbBool    no_ibc   = !frame_header->allow_intrabc;
    LrParams *lr_param = frame_header->lr_params;
    EbBool    do_lr    = no_ibc &&
        (lr_param[AOM_PLANE_Y].frame_restoration_type != RESTORE_NONE ||
         lr_param[AOM_PLANE_U].frame_restoration_type != RESTORE_NONE ||
         lr_param[AOM_PLANE_V].frame_restoration_type != RESTORE_NONE);
    EbBool    do_upscale = no_ibc &&
        !av1_superres_unscaled(&dec_handle->frame_header.frame_size);
    int       th_cnt     = NULL == thread_ctxt ? 0 : thread_ctxt->thread_cnt;
    while (1) {
        sb_row = get_sb_row_to_process(&dec_mt_frame_data->lr_sb_row_info);
        if (-1 != sb_row) {
            /* Ensure all CDEF jobs are over for row_index row  */
            volatile int32_t *start_lr =
                (volatile int32_t *)&dec_mt_frame_data->
                cdef_completed_for_row_map[sb_row];
            while (!*start_lr)
                ;

            LrCtxt * lr_ctxt = (LrCtxt *)dec_handle->pv_lr_ctxt;

            uint8_t *dst     = NULL == thread_ctxt ? lr_ctxt->dst : thread_ctxt->dst;

            if (do_lr && !do_upscale) {
                // In this case, we should only use CDEF pixels at the top
                // and bottom of the frame as a whole; internal tile boundaries
                // can use deblocked pixels from adjacent tiles for context.
                if (sb_row == 0 || sb_row == dec_mt_frame_data->sb_rows - 1) {
                    dec_save_CDEF_boundary_lines_SB_row(dec_handle,
                                                        tile_rect_p,
                                                        sb_row,
                                                        curr_blk_recon_buf,
                                                        curr_recon_stride,
                                                        num_planes);
                }
            }

            /* Pad LR_PAD_SIDE pixels for each row before the
               LR process starts for the current row. */
            pad_pre_lr(recon_picture_buf,
                       sb_row,
                       sb_size,
                       num_rows,
                       &curr_blk_recon_buf[AOM_PLANE_Y],
                       &recon_stride[AOM_PLANE_Y],
                       frame_width,
                       frame_height,
                       sx,
                       sy);

            /* Row level LR */
            if (do_lr)
                dec_av1_loop_restoration_filter_row(dec_handle,
                                                    sb_row,
                                                    &curr_blk_recon_buf[AOM_PLANE_Y],
                                                    &curr_recon_stride[AOM_PLANE_Y],
                                                    tile_rect,
                                                    0 /*opt_lr*/,
                                                    dst,
                                                    th_cnt);

            /* Pad pixels for the previous row to avoid recon buffer */
            pad_post_lr(recon_picture_buf,
                        sb_row,
                        sb_size,
                        num_rows,
                        &recon_stride[AOM_PLANE_Y],
                        pad_width,
                        pad_height,
                        shift,
                        frame_width,
                        frame_height,
                        sx,
                        sy);

            /* Update LR done map */
            dec_mt_frame_data->lr_row_map[sb_row] = 1;
        } else
            break;
    }

    eb_block_on_mutex(dec_mt_frame_data->temp_mutex);
    dec_mt_frame_data->num_threads_lred++;
    if (dec_handle->dec_config.threads == dec_mt_frame_data->num_threads_lred) {
        dec_mt_frame_data->start_motion_proj  = EB_FALSE;
        dec_mt_frame_data->start_parse_frame  = EB_FALSE;
        dec_mt_frame_data->start_decode_frame = EB_FALSE;
        dec_mt_frame_data->start_lf_frame     = EB_FALSE;
        dec_mt_frame_data->start_cdef_frame   = EB_FALSE;
        dec_mt_frame_data->start_lr_frame     = EB_FALSE;
    }
    eb_release_mutex(dec_mt_frame_data->temp_mutex);

    volatile uint32_t *num_threads_lred = &dec_mt_frame_data->num_threads_lred;
    while (*num_threads_lred != dec_handle->dec_config.threads &&
            EB_FALSE == dec_mt_frame_data->end_flag);
}

void *dec_all_stage_kernel(void *input_ptr) {
    // Context
    DecThreadCtxt * thread_ctxt    = (DecThreadCtxt *)input_ptr;
    EbDecHandle *   dec_handle_ptr = thread_ctxt->dec_handle_ptr;
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile EbBool *start_thread = (volatile EbBool *)&dec_handle_ptr->start_thread_process;
    while (*start_thread == EB_FALSE)
        ;

    while (1) {
        /* Motion Field Projection */
        svt_setup_motion_field(dec_handle_ptr, thread_ctxt);
        /* Parse Tiles */
        parse_frame_tiles(dec_handle_ptr, thread_ctxt);
        /* Decode Tiles */
        decode_frame_tiles(dec_handle_ptr, thread_ctxt);
        dec_av1_loop_filter_frame_mt(dec_handle_ptr,
                                     dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf,
                                     dec_handle_ptr->pv_lf_ctxt,
                                     AOM_PLANE_Y,
                                     MAX_MB_PLANE,
                                     thread_ctxt);
        /*Frame CDEF*/
        svt_cdef_frame_mt(dec_handle_ptr, thread_ctxt);

        /*Frame LR */
        dec_av1_loop_restoration_filter_frame_mt(dec_handle_ptr, thread_ctxt);

        if (EB_TRUE == dec_mt_frame_data->end_flag) {
            eb_block_on_mutex(dec_mt_frame_data->temp_mutex);
            dec_mt_frame_data->num_threads_exited++;
            eb_release_mutex(dec_mt_frame_data->temp_mutex);
            break;
        }
    }
    return NULL;
}

void dec_sync_all_threads(EbDecHandle *dec_handle_ptr) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    dec_mt_frame_data->end_flag = EB_TRUE;

    /* To make all worker exit except main thread! */
    dec_mt_frame_data->num_threads_cdefed = 1;
    dec_mt_frame_data->num_threads_lred   = 1;

    /* To make all worker exit except main thread! */
    dec_mt_frame_data->num_threads_header          = 1;
    dec_handle_ptr->frame_header.use_ref_frame_mvs = 0;
    dec_mt_frame_data->start_motion_proj           = EB_TRUE;

    dec_mt_frame_data->start_parse_frame = EB_TRUE;
    eb_post_semaphore(dec_handle_ptr->thread_semaphore);
    for (uint32_t lib_thrd = 0; lib_thrd < dec_handle_ptr->dec_config.threads - 1; lib_thrd++)
        eb_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
    dec_mt_frame_data->start_decode_frame = EB_TRUE;
    eb_post_semaphore(dec_handle_ptr->thread_semaphore);
    for (uint32_t lib_thrd = 0; lib_thrd < dec_handle_ptr->dec_config.threads - 1; lib_thrd++)
        eb_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
    dec_mt_frame_data->start_lf_frame = EB_TRUE;
    eb_post_semaphore(dec_handle_ptr->thread_semaphore);
    for (uint32_t lib_thrd = 0; lib_thrd < dec_handle_ptr->dec_config.threads - 1; lib_thrd++)
        eb_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
    dec_mt_frame_data->start_cdef_frame = EB_TRUE;
    eb_post_semaphore(dec_handle_ptr->thread_semaphore);
    for (uint32_t lib_thrd = 0; lib_thrd < dec_handle_ptr->dec_config.threads - 1; lib_thrd++)
        eb_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
    dec_mt_frame_data->start_lr_frame = EB_TRUE;
    eb_post_semaphore(dec_handle_ptr->thread_semaphore);
    for (uint32_t lib_thrd = 0; lib_thrd < dec_handle_ptr->dec_config.threads - 1; lib_thrd++)
        eb_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);

    while (dec_mt_frame_data->num_threads_exited != dec_handle_ptr->dec_config.threads - 1)
        eb_sleep_ms(5);

    /*Destroying lib created thread's*/
    EB_DESTROY_THREAD_ARRAY(dec_handle_ptr->decode_thread_handle_array,
                            dec_handle_ptr->dec_config.threads - 1);
}
