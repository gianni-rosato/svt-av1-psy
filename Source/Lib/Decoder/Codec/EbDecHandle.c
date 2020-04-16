/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the API component functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecMemInit.h"
#include "EbDecPicMgr.h"
#include "grainSynthesis.h"

#ifndef _WIN32
#include <pthread.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#include "common_dsp_rtcd.h"

#include "EbLog.h"

/**************************************
* Globals
**************************************/
uint8_t num_groups = 0;
#ifdef _WIN32
GROUP_AFFINITY group_affinity;
EbBool         alternate_groups = 0;
#elif defined(__linux__)
cpu_set_t group_affinity;
typedef struct logicalProcessorGroup {
    uint32_t num;
    uint32_t group[1024];
} processorGroup;
#define INITIAL_PROCESSOR_GROUP 16
processorGroup *lp_group = NULL;
#endif

EbMemoryMapEntry *svt_dec_memory_map;
uint32_t *        svt_dec_memory_map_index;
uint64_t *        svt_dec_total_lib_memory;

uint32_t svt_dec_lib_malloc_count = 0;

//TODO: Should be removed! Check
EbMemoryMapEntry *memory_map;
uint32_t *        memory_map_index;
uint64_t *        total_lib_memory;

uint32_t lib_malloc_count    = 0;
uint32_t lib_semaphore_count = 0;
uint32_t lib_mutex_count     = 0;

/* Track address for memory during change in resoultion*/
EbMemoryMapEntry    *memory_map_start_address;
EbMemoryMapEntry    *memory_map_end_address;

void        asm_set_convolve_asm_table(void);
void        init_intra_dc_predictors_c_internal(void);
void        asm_set_convolve_hbd_asm_table(void);
void        init_intra_predictors_internal(void);
extern void av1_init_wedge_masks(void);
void        dec_sync_all_threads(EbDecHandle *dec_handle_ptr);

EbErrorType decode_multiple_obu(EbDecHandle *dec_handle_ptr, uint8_t **data, size_t data_size,
                                uint32_t is_annexb);

void switch_to_real_time() {
#ifndef _WIN32

    struct sched_param schedParam = {.sched_priority = 99};

    int32_t retValue = pthread_setschedparam(pthread_self(), SCHED_FIFO, &schedParam);
    UNUSED(retValue);
#endif
}

/***********************************
* Decoder Library Handle Constructor
************************************/
/*TODO : Add more features*/
static EbErrorType eb_dec_handle_ctor(EbDecHandle **decHandleDblPtr, EbComponentType *ebHandlePtr) {
    (void)ebHandlePtr;
    EbErrorType return_error = EB_ErrorNone;

    // Allocate Memory
    EbDecHandle *dec_handle_ptr = (EbDecHandle *)malloc(sizeof(EbDecHandle));
    *decHandleDblPtr            = dec_handle_ptr;
    if (dec_handle_ptr == (EbDecHandle *)EB_NULL) return EB_ErrorInsufficientResources;
    dec_handle_ptr->memory_map       = (EbMemoryMapEntry *)malloc(sizeof(EbMemoryMapEntry));
    dec_handle_ptr->memory_map_index = 0;
    dec_handle_ptr->total_lib_memory =
        sizeof(EbComponentType) + sizeof(EbDecHandle) + sizeof(EbMemoryMapEntry);
    dec_handle_ptr->memory_map_init_address = dec_handle_ptr->memory_map;
    // Save Memory Map Pointers
    svt_dec_total_lib_memory = &dec_handle_ptr->total_lib_memory;
    svt_dec_memory_map       = dec_handle_ptr->memory_map;
    svt_dec_memory_map_index = &dec_handle_ptr->memory_map_index;
    svt_dec_lib_malloc_count = 0;

    dec_handle_ptr->start_thread_process = EB_FALSE;
    memory_map_start_address = NULL;
    memory_map_end_address = NULL;

    return return_error;
}

/* FilmGrain module req. even dim. for internal operation
   For odd luma, chroma uses data from even pixel in luma*/
static void copy_even(uint8_t *luma, uint32_t wd,
                 uint32_t ht, uint32_t stride, int32_t use_hbd) {
    if ((wd & 1) == 0 && (ht & 1) == 0) return;
    if (wd & 1) {
        for (uint32_t i = 0; i < ht; ++i)
            luma[i * (stride << use_hbd) + (wd << use_hbd)] =
                luma[i *(stride << use_hbd) + ((wd - 1) << use_hbd)];
        wd = wd + 1;
    }
    if (ht & 1) {
        memcpy(&luma[ht *  (stride << use_hbd)], &luma[(ht - 1) * (stride << use_hbd)],
            sizeof(*luma) * (wd << use_hbd));
    }
}
/* Copy from recon buffer to out buffer! */
int svt_dec_out_buf(EbDecHandle *dec_handle_ptr, EbBufferHeaderType *p_buffer) {
    EbPictureBufferDesc *recon_picture_buf = dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf;
    EbSvtIOFormat *      out_img           = (EbSvtIOFormat *)p_buffer->p_buffer;

    uint8_t *luma = NULL;
    uint8_t *cb   = NULL;
    uint8_t *cr   = NULL;

    /* TODO: Should add logic for show_existing_frame */
    if (0 == dec_handle_ptr->show_frame) {
        assert(0 == dec_handle_ptr->show_existing_frame);
        return 0;
    }

    uint32_t wd = dec_handle_ptr->frame_header.frame_size.superres_upscaled_width;
    uint32_t ht = dec_handle_ptr->frame_header.frame_size.frame_height;
    uint32_t i, sx = 0, sy = 0;
    /* FilmGrain module req. even dim. for internal operation */
    int even_w = (wd & 1) ? (wd + 1) : wd;
    int even_h = (ht & 1) ? (ht + 1) : ht;

    if (out_img->height != ht || out_img->width != wd ||
        out_img->color_fmt != recon_picture_buf->color_format ||
        out_img->bit_depth != (EbBitDepth)recon_picture_buf->bit_depth) {
        int size = (dec_handle_ptr->seq_header.color_config.bit_depth == EB_EIGHT_BIT)
                       ? sizeof(uint8_t)
                       : sizeof(uint16_t);

        int luma_size = size * even_w * even_h;
        int chroma_size = -1;
        out_img->color_fmt = recon_picture_buf->color_format;
        switch (recon_picture_buf->color_format) {
        case EB_YUV400:
            out_img->cb_stride = INT32_MAX;
            out_img->cr_stride = INT32_MAX;
            break;
        case EB_YUV420:
            out_img->cb_stride = (wd + 1) >> 1;
            out_img->cr_stride = (wd + 1) >> 1;
            chroma_size        = size * (((wd + 1) >> 1) * ((ht + 1) >> 1));
            break;
        case EB_YUV422:
            out_img->cb_stride = (wd + 1) >> 1;
            out_img->cr_stride = (wd + 1) >> 1;
            chroma_size        = size * (((wd + 1) >> 1) * ht);
            break;
        case EB_YUV444:
            out_img->cb_stride = wd;
            out_img->cr_stride = wd;
            chroma_size        = size * ht * wd;
            break;
        default: SVT_LOG("Unsupported colour format. \n"); return 0;
        }

        /* FilmGrain module req. even dim. for internal operation */
        out_img->y_stride = even_w;
        out_img->width    = wd;
        out_img->height   = ht;
        if (out_img->bit_depth != (EbBitDepth)recon_picture_buf->bit_depth) {
            SVT_LOG(
                "Warning : Output bit depth conversion not supported."
                " Output depth set to %d. ",
                recon_picture_buf->bit_depth);
            out_img->bit_depth = (EbBitDepth)recon_picture_buf->bit_depth;
        }

        free(out_img->luma);
        if (recon_picture_buf->color_format != EB_YUV400) {
            free(out_img->cb);
            free(out_img->cr);
        }
        out_img->luma = (uint8_t *)malloc(luma_size);
        if (recon_picture_buf->color_format != EB_YUV400) {
            out_img->cb = (uint8_t *)malloc(chroma_size);
            out_img->cr = (uint8_t *)malloc(chroma_size);
        }
    }

    switch (recon_picture_buf->color_format) {
    case EB_YUV400:
        sx = -1;
        sy = -1;
        break;
    case EB_YUV420:
        sx = 1;
        sy = 1;
        break;
    case EB_YUV422:
        sx = 1;
        sy = 0;
        break;
    case EB_YUV444:
        sx = 0;
        sy = 0;
        break;
    default: assert(0);
    }

    int32_t use_high_bit_depth = recon_picture_buf->bit_depth == EB_8BIT ? 0 : 1;

    luma = out_img->luma +
           ((out_img->origin_y * out_img->y_stride + out_img->origin_x) << use_high_bit_depth);
    if (recon_picture_buf->color_format != EB_YUV400) {
        cb = out_img->cb +
             ((out_img->cb_stride * (out_img->origin_y >> sy) + (out_img->origin_x >> sx))
              << use_high_bit_depth);
        cr = out_img->cr +
             ((out_img->cr_stride * (out_img->origin_y >> sy) + (out_img->origin_x >> sx))
              << use_high_bit_depth);
    }

    /* Memcpy to dst buffer */
    {
        if (recon_picture_buf->bit_depth == EB_8BIT) {

            if (dec_handle_ptr->is_16bit_pipeline) {
                uint8_t *dst;
                uint16_t *pu2_src;
                uint32_t  j;
                ASSERT(recon_picture_buf->is_16bit_pipeline);

                /* Luma */
                dst = luma;
                pu2_src = (uint16_t *)recon_picture_buf->buffer_y + recon_picture_buf->origin_x +
                    (recon_picture_buf->origin_y * recon_picture_buf->stride_y);

                for (i = 0; i < ht; i++) {
                    for (j = 0; j < wd; j++)
                        dst[j] = (uint8_t)pu2_src[j];
                    dst += out_img->y_stride;
                    pu2_src += recon_picture_buf->stride_y;
                }

                if (recon_picture_buf->color_format != EB_YUV400) {
                    /* Cb */
                    dst = cb;
                    pu2_src = (uint16_t *)recon_picture_buf->buffer_cb +
                        (recon_picture_buf->origin_x >> sx) +
                        ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cb);

                    for (i = 0; i < ((ht + sy) >> sy); i++) {
                        for (j = 0; j < ((wd + sx) >> sx); j++)
                            dst[j] = (uint8_t)pu2_src[j];
                        dst += out_img->cb_stride;
                        pu2_src += recon_picture_buf->stride_cb;
                    }

                    /* Cr */
                    dst = cr;
                    pu2_src = (uint16_t *)recon_picture_buf->buffer_cr +
                        (recon_picture_buf->origin_x >> sx) +
                        ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cr);

                    for (i = 0; i < ((ht + sy) >> sy); i++) {
                        for (j = 0; j < ((wd + sx) >> sx); j++)
                            dst[j] = (uint8_t)pu2_src[j];
                        dst += out_img->cr_stride;
                        pu2_src += recon_picture_buf->stride_cr;
                    }
                }
            }
            else {
            uint8_t *src, *dst;
            dst = luma;
            src = recon_picture_buf->buffer_y + recon_picture_buf->origin_x +
                (recon_picture_buf->origin_y * recon_picture_buf->stride_y);

            for (i = 0; i < ht; i++) {
                memcpy(dst, src, wd);
                dst += out_img->y_stride;
                src += recon_picture_buf->stride_y;
            }

            if (recon_picture_buf->color_format != EB_YUV400) {
                /* Cb */
                dst = cb;
                src = recon_picture_buf->buffer_cb + (recon_picture_buf->origin_x >> sx) +
                    ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cb);

                for (i = 0; i < ((ht + sy) >> sy); i++) {
                    memcpy(dst, src, ((wd + sx) >> sx));
                    dst += out_img->cb_stride;
                    src += recon_picture_buf->stride_cb;
                }

                /* Cr */
                dst = cr;
                src = recon_picture_buf->buffer_cr + (recon_picture_buf->origin_x >> sx) +
                    ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cr);

                for (i = 0; i < ((ht + sy) >> sy); i++) {
                    memcpy(dst, src, ((wd + sx) >> sx));
                    dst += out_img->cr_stride;
                    src += recon_picture_buf->stride_cr;
                }
            }
        }
        } else {
            uint16_t *pu2_dst;
            uint16_t *pu2_src;

            /* Luma */
            pu2_dst = (uint16_t *)luma;
            pu2_src = (uint16_t *)recon_picture_buf->buffer_y + recon_picture_buf->origin_x +
                      (recon_picture_buf->origin_y * recon_picture_buf->stride_y);

            for (i = 0; i < ht; i++) {
                memcpy(pu2_dst, pu2_src, sizeof(uint16_t) * wd);
                pu2_dst += out_img->y_stride;
                pu2_src += recon_picture_buf->stride_y;
            }

            if (recon_picture_buf->color_format != EB_YUV400) {
                /* Cb */
                pu2_dst = (uint16_t *)cb;
                pu2_src = (uint16_t *)recon_picture_buf->buffer_cb +
                          (recon_picture_buf->origin_x >> sx) +
                          ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cb);

                for (i = 0; i < ((ht + sy) >> sy); i++) {
                    memcpy(pu2_dst, pu2_src, sizeof(uint16_t) * ((wd + sx) >> sx));
                    pu2_dst += out_img->cb_stride;
                    pu2_src += recon_picture_buf->stride_cb;
                }

                /* Cr */
                pu2_dst = (uint16_t *)cr;
                pu2_src = (uint16_t *)recon_picture_buf->buffer_cr +
                          (recon_picture_buf->origin_x >> sx) +
                          ((recon_picture_buf->origin_y >> sy) * recon_picture_buf->stride_cr);

                for (i = 0; i < ((ht + sy) >> sy); i++) {
                    memcpy(pu2_dst, pu2_src, sizeof(uint16_t) * ((wd + sx) >> sx));
                    pu2_dst += out_img->cr_stride;
                    pu2_src += recon_picture_buf->stride_cr;
                }
            }
        }
    }

    if (!dec_handle_ptr->dec_config.skip_film_grain) {
        /* Need to fill the dst buf with recon data before calling film_grain */
        AomFilmGrain *film_grain_ptr = &dec_handle_ptr->cur_pic_buf[0]->film_grain_params;
        if (film_grain_ptr->apply_grain) {
            switch (recon_picture_buf->bit_depth) {
            case EB_8BIT: film_grain_ptr->bit_depth = 8; break;
            case EB_10BIT: film_grain_ptr->bit_depth = 10; break;
            default: assert(0);
            }
            copy_even(luma, wd, ht, out_img->y_stride, use_high_bit_depth);
            eb_av1_add_film_grain_run(film_grain_ptr,
                                      luma,
                                      cb,
                                      cr,
                                      even_h,/*(ht & 1 ? ht + 1 : ht),*/
                                      even_w,/*(wd & 1 ? wd + 1 : ht),*/
                                      out_img->y_stride,
                                      out_img->cb_stride,
                                      use_high_bit_depth,
                                      sy,
                                      sx);
        }
    }

    return 1;
}

/**********************************
Set Default Library Params
**********************************/
EbErrorType eb_svt_dec_set_default_parameter(EbSvtAv1DecConfiguration *config_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    if (config_ptr == NULL) return EB_ErrorBadParameter;

    config_ptr->operating_point           = -1;
    config_ptr->output_all_layers         = 0;
    config_ptr->skip_film_grain           = 0;
    config_ptr->skip_frames               = 0;
    config_ptr->frames_to_be_decoded      = 0;
    config_ptr->compressed_ten_bit_format = 0;
    config_ptr->eight_bit_output          = 0;

    /* Picture parameters */
    config_ptr->max_picture_width  = 0;
    config_ptr->max_picture_height = 0;
    config_ptr->max_bit_depth      = EB_EIGHT_BIT;
    config_ptr->is_16bit_pipeline = 0;
    config_ptr->max_color_format   = EB_YUV420;
    config_ptr->threads            = 1;

    // Application Specific parameters
    config_ptr->channel_id           = 0;
    config_ptr->active_channel_count = 1;
    config_ptr->stat_report          = 0;

    /* Multi-thread parameters */
    config_ptr->threads      = 1;
    config_ptr->num_p_frames = 1;

    return return_error;
}

/**********************************
* Decoder Handle Initialization
**********************************/
static EbErrorType init_svt_av1_decoder_handle(EbComponentType *hComponent) {
    EbErrorType      return_error      = EB_ErrorNone;
    EbComponentType *svt_dec_component = (EbComponentType *)hComponent;

    SVT_LOG("SVT [version]:\tSVT-AV1 Decoder Lib v%d.%d.%d\n",
            SVT_VERSION_MAJOR,
            SVT_VERSION_MINOR,
            SVT_VERSION_PATCHLEVEL);
#if ( defined( _MSC_VER ) && (_MSC_VER >= 1920) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2019");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1910) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2017");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1900) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2015");
#elif ( defined( _MSC_VER ) && (_MSC_VER < 1900) )
    SVT_LOG("SVT [build]  :\tVisual Studio (old)");
#elif defined(__GNUC__)
    SVT_LOG("SVT [build]  :\tGCC %d.%d.%d\t", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    SVT_LOG("SVT [build]  :\tunknown compiler");
#endif
    SVT_LOG(" %u bit\n", (unsigned)sizeof(void *) * 8);
    SVT_LOG("LIB Build date: %s %s\n", __DATE__, __TIME__);
    SVT_LOG("-------------------------------------------\n");

    switch_to_real_time();

    // Set Component Size & Version
    svt_dec_component->size = sizeof(EbComponentType);

    // Decoder Private Handle Ctor
    return_error = (EbErrorType)eb_dec_handle_ctor(
        (EbDecHandle **)&(svt_dec_component->p_component_private), svt_dec_component);

    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_init_handle(EbComponentType **p_handle, void *p_app_data,
                   EbSvtAv1DecConfiguration *config_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    if (p_handle == NULL) return EB_ErrorBadParameter;

    svt_log_init();

    *p_handle = (EbComponentType *)malloc(sizeof(EbComponentType));

    if (*p_handle != (EbComponentType *)NULL) {
        // Init Component OS objects (threads, semaphores, etc.)
        // also links the various Component control functions
        return_error = init_svt_av1_decoder_handle(*p_handle);

        if (return_error == EB_ErrorNone)
            ((EbComponentType *)(*p_handle))->p_application_private = p_app_data;
        else if (return_error == EB_ErrorInsufficientResources) {
            svt_av1_dec_deinit((EbComponentType *)NULL);
            *p_handle = (EbComponentType *)NULL;
        } else
            return_error = EB_ErrorInvalidComponent;
    } else {
        //SVT_LOG("Error: Component Struct Malloc Failed\n");
        return_error = EB_ErrorInsufficientResources;
    }

    if (return_error == EB_ErrorNone) return_error = eb_svt_dec_set_default_parameter(config_ptr);

    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_set_parameter(EbComponentType *         svt_dec_component,
                         EbSvtAv1DecConfiguration *config_struct) {
    if (svt_dec_component == NULL || config_struct == NULL) return EB_ErrorBadParameter;

    EbDecHandle *dec_handle_ptr = (EbDecHandle *)svt_dec_component->p_component_private;

    dec_handle_ptr->dec_config = *config_struct;
    dec_handle_ptr->is_16bit_pipeline = config_struct->is_16bit_pipeline;

    return EB_ErrorNone;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_init(EbComponentType *svt_dec_component) {
    EbErrorType return_error = EB_ErrorNone;
    if (svt_dec_component == NULL) return EB_ErrorBadParameter;

    EbDecHandle *dec_handle_ptr = (EbDecHandle *)svt_dec_component->p_component_private;
#ifdef ARCH_X86
    CPU_FLAGS    cpu_flags = get_cpu_flags_to_use();
#else
    CPU_FLAGS    cpu_flags = 0;
#endif
    dec_handle_ptr->dec_cnt       = -1;
    dec_handle_ptr->num_frms_prll = 1;
    if (dec_handle_ptr->num_frms_prll > DEC_MAX_NUM_FRM_PRLL)
        dec_handle_ptr->num_frms_prll = DEC_MAX_NUM_FRM_PRLL;
    dec_handle_ptr->seq_header_done = 0;
    dec_handle_ptr->mem_init_done   = 0;

    dec_handle_ptr->seen_frame_header   = 0;
    dec_handle_ptr->show_existing_frame = 0;
    dec_handle_ptr->show_frame          = 0;
    dec_handle_ptr->showable_frame      = 0;
    dec_handle_ptr->seq_header.sb_size = 0;

    setup_common_rtcd_internal(cpu_flags);

    asm_set_convolve_asm_table();

    init_intra_dc_predictors_c_internal();

    asm_set_convolve_hbd_asm_table();

    init_intra_predictors_internal();

    av1_init_wedge_masks();

    /************************************
    * Decoder Memory Init
    ************************************/
    return_error = dec_mem_init(dec_handle_ptr);
    if (return_error != EB_ErrorNone) return return_error;

    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_frame(EbComponentType *svt_dec_component, const uint8_t *data, const size_t data_size,
                    uint32_t is_annexb) {
    EbErrorType return_error = EB_ErrorNone;
    if (svt_dec_component == NULL) return EB_ErrorBadParameter;

    EbDecHandle *dec_handle_ptr       = (EbDecHandle *)svt_dec_component->p_component_private;
    uint8_t *    data_start           = (uint8_t *)data;
    uint8_t *    data_end             = (uint8_t *)data + data_size;
    dec_handle_ptr->seen_frame_header = 0;

    while (data_start < data_end) {
        /*TODO : Remove or move. For Test purpose only */
        dec_handle_ptr->dec_cnt++;
        //SVT_LOG("\n SVT-AV1 Dec : Decoding Pic #%d", dec_handle_ptr->dec_cnt);

        uint64_t frame_size = 0;
        frame_size          = data_end - data_start;
        return_error = decode_multiple_obu(dec_handle_ptr, &data_start, frame_size, is_annexb);

        if (return_error != EB_ErrorNone) assert(0);

        dec_pic_mgr_update_ref_pic(dec_handle_ptr,
                                   (EB_ErrorNone == return_error) ? 1 : 0,
                                   dec_handle_ptr->frame_header.refresh_frame_flags);

        // Allow extra zero bytes after the frame end
        while (data < data_end) {
            const uint8_t marker = data[0];
            if (marker) break;
            ++data;
        }

        /*SVT_LOG("\nDecoding Pic #%d  frm_w : %d    frm_h : %d
            frm_typ : %d", dec_handle_ptr->dec_cnt,
            dec_handle_ptr->frame_header.frame_size.frame_width,
            dec_handle_ptr->frame_header.frame_size.frame_height,
            dec_handle_ptr->frame_header.frame_type);*/
    }

    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_get_picture(EbComponentType *svt_dec_component, EbBufferHeaderType *p_buffer,
                       EbAV1StreamInfo *stream_info, EbAV1FrameInfo *frame_info) {
    (void)stream_info;
    (void)frame_info;

    EbErrorType return_error = EB_ErrorNone;
    if (svt_dec_component == NULL) return EB_ErrorBadParameter;

    EbDecHandle *dec_handle_ptr = (EbDecHandle *)svt_dec_component->p_component_private;
    /* Copy from recon pointer and return! TODO: Should remove the memcpy! */
    if (0 == svt_dec_out_buf(dec_handle_ptr, p_buffer)) return_error = EB_DecNoOutputPicture;
    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_deinit(EbComponentType *svt_dec_component) {
    if (svt_dec_component == NULL) return EB_ErrorBadParameter;
    EbDecHandle *dec_handle_ptr = (EbDecHandle *)svt_dec_component->p_component_private;
    EbErrorType  return_error   = EB_ErrorNone;

    if (dec_handle_ptr) {
        if (dec_handle_ptr->dec_config.threads > 1) dec_sync_all_threads(dec_handle_ptr);
        if (svt_dec_memory_map) {
            // Loop through the ptr table and free all malloc'd pointers per channel
            EbMemoryMapEntry *memory_entry = svt_dec_memory_map;
            if (memory_entry) {
                do {
                    switch (memory_entry->ptr_type) {
                    case EB_N_PTR: free(memory_entry->ptr); break;
                    case EB_A_PTR:
#ifdef _WIN32
                        _aligned_free(memory_entry->ptr);
#else
                        free(memory_entry->ptr);
#endif
                        break;
                    case EB_SEMAPHORE: eb_destroy_semaphore(memory_entry->ptr); break;
                    case EB_THREAD: eb_destroy_thread(memory_entry->ptr); break;
                    case EB_MUTEX: eb_destroy_mutex(memory_entry->ptr); break;
                    default: return_error = EB_ErrorMax; break;
                    }
                    EbMemoryMapEntry *tmp_memory_entry = memory_entry;
                    memory_entry = (EbMemoryMapEntry *)tmp_memory_entry->prev_entry;
                    if (tmp_memory_entry) free(tmp_memory_entry);
                } while (memory_entry != dec_handle_ptr->memory_map_init_address && memory_entry);
                if (dec_handle_ptr->memory_map_init_address)
                    free(dec_handle_ptr->memory_map_init_address);
            }
        }
    }
    return return_error;
}

/**********************************
* Encoder Componenet DeInit
**********************************/
EbErrorType eb_dec_component_de_init(EbComponentType *svt_dec_component) {
    EbErrorType return_error = EB_ErrorNone;

    if (svt_dec_component->p_component_private)
        free((EbDecHandle *)svt_dec_component->p_component_private);
    else
        return_error = EB_ErrorUndefined;
    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType
svt_av1_dec_deinit_handle(EbComponentType *svt_dec_component) {
    EbErrorType return_error = EB_ErrorNone;

    if (svt_dec_component) {
        return_error = eb_dec_component_de_init(svt_dec_component);

        free(svt_dec_component);
    } else
        return_error = EB_ErrorInvalidComponent;
    return return_error;
}
