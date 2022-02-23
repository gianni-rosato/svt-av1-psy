// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/
// SUMMARY
//   Contains the API component functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "EbVersion.h"
#include "EbThreads.h"
#include "EbUtility.h"
#include "EbEncHandle.h"
#include "EbEncSettings.h"
#include "EbPictureControlSet.h"
#include "EbPictureOperators.h"
#include "EbReferenceObject.h"
#include "EbResourceCoordinationProcess.h"
#include "EbPictureAnalysisProcess.h"
#include "EbPictureDecisionProcess.h"
#include "EbMotionEstimationProcess.h"
#include "EbInitialRateControlProcess.h"
#include "EbSourceBasedOperationsProcess.h"
#include "EbPictureManagerProcess.h"
#include "EbRateControlProcess.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "EbEncDecProcess.h"
#include "EbEntropyCodingProcess.h"
#include "EbPacketizationProcess.h"
#include "EbResourceCoordinationResults.h"
#include "EbPictureAnalysisResults.h"
#include "EbPictureDecisionResults.h"
#include "EbMotionEstimationResults.h"
#include "EbInitialRateControlResults.h"
#include "EbPictureDemuxResults.h"
#include "EbRateControlTasks.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbPredictionStructure.h"
#include "EbRestProcess.h"
#include "EbCdefProcess.h"
#include "EbDlfProcess.h"
#include "EbRateControlResults.h"
#include "EbDefinitions.h"

#include"EbPackUnPack_C.h"

#ifdef ARCH_X86_64
#include <immintrin.h>
#endif
#include "EbLog.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#endif

#include "aom_dsp_rtcd.h"
#include "common_dsp_rtcd.h"


 EbErrorType prediction_structure_group_ctor(
     PredictionStructureGroup* pred_struct_group_ptr,
     struct SequenceControlSet* scs_ptr);

 /**************************************
  * Defines
  **************************************/
#define EB_EncodeInstancesTotalCount                    1
#define EB_ComputeSegmentInitCount                      1

  // Config Set Initial Count
#define EB_SequenceControlSetPoolInitCount              3

// Process Instantiation Initial Counts
#define EB_ResourceCoordinationProcessInitCount         1
#define EB_PictureDecisionProcessInitCount              1
#define EB_InitialRateControlProcessInitCount           1
#define EB_PictureManagerProcessInitCount               1
#define EB_RateControlProcessInitCount                  1
#define EB_PacketizationProcessInitCount                1

// Output Buffer Transfer Parameters
#define EB_OUTPUTSTREAMBUFFERSIZE                                       0x2DC6C0   //0x7D00        // match MTU Size
#define EB_OUTPUTSTATISTICSBUFFERSIZE                                   0x30            // 6X8 (8 Bytes for Y, U, V, number of bits, picture number, QP)
#define EOS_NAL_BUFFER_SIZE                                             0x0010 // Bitstream used to code EOS NAL
#define TPL_INPUT_PORT_SOP                                0
#define TPL_INPUT_PORT_TPL                                1
#define TPL_INPUT_PORT_INVALID                           -1
#define ENCDEC_INPUT_PORT_MDC                                0
#define ENCDEC_INPUT_PORT_ENCDEC                             1
#define ENCDEC_INPUT_PORT_INVALID                           -1
/**************************************
 * Globals
 **************************************/
static uint8_t                   num_groups = 0;
#ifdef _WIN32
static GROUP_AFFINITY            group_affinity;
static EbBool                    alternate_groups = 0;
#elif defined(__linux__)
static cpu_set_t                 group_affinity;
typedef struct logicalProcessorGroup {
    uint32_t num;
    uint32_t group[1024];
} processorGroup;
#define INITIAL_PROCESSOR_GROUP 16
static processorGroup           *lp_group = NULL;
#endif
uint8_t  get_tpl_synthesizer_block_size(int8_t tpl_level, uint32_t picture_width, uint32_t picture_height);

uint8_t get_disallow_nsq(EbEncMode enc_mode);
uint8_t get_disallow_4x4(EbEncMode enc_mode, EB_SLICE slice_type);
extern uint32_t tot_past_refs[];
uint32_t  get_num_refs_in_one_mg(PredictionStructure *pred_struct_ptr);

static const char *get_asm_level_name_str(CPU_FLAGS cpu_flags) {

    const struct {
        const char *name;
        CPU_FLAGS flags;
    } param_maps[] = {
        {"c",       0},
        {"mmx",     CPU_FLAGS_MMX},
        {"sse",     CPU_FLAGS_SSE},
        {"sse2",    CPU_FLAGS_SSE2},
        {"sse3",    CPU_FLAGS_SSE3},
        {"ssse3",   CPU_FLAGS_SSSE3},
        {"sse4_1",  CPU_FLAGS_SSE4_1},
        {"sse4_2",  CPU_FLAGS_SSE4_2},
        {"avx",     CPU_FLAGS_AVX},
        {"avx2",    CPU_FLAGS_AVX2},
        {"avx512",  CPU_FLAGS_AVX512F}
    };
    const uint32_t para_map_size = sizeof(param_maps) / sizeof(param_maps[0]);
    int32_t i;

    for (i = para_map_size -1; i >= 0; --i) {
        if (param_maps[i].flags & cpu_flags) {
            return param_maps[i].name;
        }
    }
    return "c";
}

//Get Number of logical processors
uint32_t get_num_processors() {
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return num_groups == 1 ? sysinfo.dwNumberOfProcessors : sysinfo.dwNumberOfProcessors << 1;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

EbErrorType init_thread_management_params() {
#ifdef _WIN32
    // Initialize group_affinity structure with Current thread info
    GetThreadGroupAffinity(GetCurrentThread(), &group_affinity);
    num_groups = (uint8_t)GetActiveProcessorGroupCount();
#elif defined(__linux__)
    memset(lp_group, 0, INITIAL_PROCESSOR_GROUP * sizeof(processorGroup));

    FILE *fin = fopen("/proc/cpuinfo", "r");
    if (fin) {
        int processor_id = 0;
        int maxSize = INITIAL_PROCESSOR_GROUP;
        char line[1024];
        while (fgets(line, sizeof(line), fin)) {
            if(strncmp(line, "processor", 9) == 0) {
                char* p = line + 9;
                while(*p < '0' || *p > '9') p++;
                processor_id = strtol(p, NULL, 0);
            }
            if(strncmp(line, "physical id", 11) == 0) {
                char* p = line + 11;
                while(*p < '0' || *p > '9') p++;
                long socket_id = strtol(p, NULL, 0);
                if (socket_id < 0) {
                    fclose(fin);
                    return EB_ErrorInsufficientResources;
                }
                if (socket_id + 1 > num_groups)
                    num_groups = socket_id + 1;
                if (socket_id >= maxSize) {
                    maxSize = maxSize * 2;
                    processorGroup *temp = realloc(lp_group, maxSize * sizeof(*temp));
                    if (temp)
                        lp_group = temp;
                    else {
                        free(lp_group);
                        fclose(fin);
                        return EB_ErrorInsufficientResources;
                    }
                }
                lp_group[socket_id].group[lp_group[socket_id].num++] = processor_id;
            }
        }
        fclose(fin);
    }
#endif
    return EB_ErrorNone;
}

#ifdef _WIN32
uint64_t get_affinity_mask(uint32_t lpnum) {
    uint64_t mask = 0x1;
    for (uint32_t i = lpnum - 1; i > 0; i--)
        mask += (uint64_t)1 << i;
    return mask;
}
#endif

void svt_set_thread_management_parameters(EbSvtAv1EncConfiguration *config_ptr)
{
#ifdef _WIN32
    const uint32_t num_logical_processors = get_num_processors();
    // For system with a single processor group(no more than 64 logic processors all together)
    // Affinity of the thread can be set to one or more logical processors
    if (num_groups == 1) {
            const uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
                config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
            group_affinity.Mask = get_affinity_mask(lps);
    }
    else if (num_groups > 1) { // For system with multiple processor group
        if (config_ptr->logical_processors == 0) {
            if (config_ptr->target_socket != -1)
                group_affinity.Group = config_ptr->target_socket;
        }
        else {
            const uint32_t num_lp_per_group = num_logical_processors / num_groups;
            if (config_ptr->target_socket == -1) {
                if (config_ptr->logical_processors > num_lp_per_group) {
                    alternate_groups = EB_TRUE;
                    SVT_WARN("-lp(logical processors) setting is ignored. Run on both sockets. \n");
                }
                else
                    group_affinity.Mask = get_affinity_mask(config_ptr->logical_processors);
            }
            else {
                const uint32_t lps =
                    config_ptr->logical_processors < num_lp_per_group ? config_ptr->logical_processors : num_lp_per_group;
                group_affinity.Mask = get_affinity_mask(lps);
                group_affinity.Group = config_ptr->target_socket;
            }
        }
    }
#elif defined(__linux__)
    uint32_t num_logical_processors = get_num_processors();
    CPU_ZERO(&group_affinity);

    if (num_groups == 1) {
        const uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
            config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
        for (uint32_t i = 0; i < lps; i++)
            CPU_SET(lp_group[0].group[i], &group_affinity);
    } else if (num_groups > 1) {
        const uint32_t num_lp_per_group = num_logical_processors / num_groups;
        if (config_ptr->logical_processors == 0) {
            if (config_ptr->target_socket != -1)
                for (uint32_t i = 0; i < lp_group[config_ptr->target_socket].num; i++)
                    CPU_SET(lp_group[config_ptr->target_socket].group[i], &group_affinity);
        } else {
            if (config_ptr->target_socket == -1) {
                const uint32_t lps =
                    config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
                if (lps > num_lp_per_group) {
                    for (uint32_t i = 0; i < lp_group[0].num; i++)
                        CPU_SET(lp_group[0].group[i], &group_affinity);
                    for (uint32_t i = 0; i < (lps - lp_group[0].num); i++)
                        CPU_SET(lp_group[1].group[i], &group_affinity);
                } else
                    for (uint32_t i = 0; i < lps; i++)
                        CPU_SET(lp_group[0].group[i], &group_affinity);
            } else {
                const uint32_t lps =
                    config_ptr->logical_processors < num_lp_per_group ? config_ptr->logical_processors : num_lp_per_group;
                for (uint32_t i = 0; i < lps; i++)
                    CPU_SET(lp_group[config_ptr->target_socket].group[i], &group_affinity);
            }
        }
    }
#else
    UNUSED(config_ptr);
    UNUSED(num_groups);
#endif
}

void asm_set_convolve_asm_table(void);
void asm_set_convolve_hbd_asm_table(void);
void init_intra_dc_predictors_c_internal(void);
void init_intra_predictors_internal(void);
void svt_av1_init_me_luts(void);

static void enc_switch_to_real_time(){
#if !defined(_WIN32)
    if (!geteuid())
        (void)pthread_setschedparam(
            pthread_self(), SCHED_FIFO, &(struct sched_param){.sched_priority = 99});
#endif
}
#define SINGLE_CORE_COUNT       1
#define CONS_CORE_COUNT         16
#define LOW_SERVER_CORE_COUNT   48
#define MED_SERVER_CORE_COUNT   128
#define HIGH_SERVER_CORE_COUNT  224
#define PARALLEL_LEVEL_2_RANGE  2
#define PARALLEL_LEVEL_4_RANGE  6
#define PARALLEL_LEVEL_8_RANGE  12
#define PARALLEL_LEVEL_16_RANGE 24
#define PARALLEL_LEVEL_32_RANGE 48

int32_t set_parent_pcs(EbSvtAv1EncConfiguration*   config, uint32_t core_count, EbInputResolution res_class) {
    if (config){
        uint32_t fps            = (uint32_t)((config->frame_rate > 1000) ?
                        config->frame_rate >> 16 :
                        config->frame_rate);
        uint32_t min_ppcs_count = (2 << config->hierarchical_levels) + 1; // min picture count to start encoding
        fps        = fps > 120 ? 120   : fps;
        fps        = fps < 24  ? 24    : fps;

        uint32_t ppcs_count = MAX(min_ppcs_count, fps);
        if (core_count <= SINGLE_CORE_COUNT)
            ppcs_count = min_ppcs_count;
        else{
            if (res_class <= INPUT_SIZE_480p_RANGE) {
                if (core_count < CONS_CORE_COUNT)
                    ppcs_count = ppcs_count * 1;                // 1 sec
                else if (core_count < LOW_SERVER_CORE_COUNT)
                    ppcs_count = (ppcs_count * 3) >> 1;         // 1.5 sec
                else if (core_count < MED_SERVER_CORE_COUNT)
                    ppcs_count = ppcs_count << 1;               // 2 sec
                else
                    ppcs_count = ppcs_count * 3;                // 3 sec
            } else if (res_class <= INPUT_SIZE_1080p_RANGE) {
                if (core_count < CONS_CORE_COUNT)
                    ppcs_count = min_ppcs_count;
                else if (core_count < LOW_SERVER_CORE_COUNT)
                    ppcs_count = (ppcs_count * 3) >> 1;         // 1.5 sec
                else if (core_count < MED_SERVER_CORE_COUNT)
                    ppcs_count = ppcs_count << 1;               // 2 sec
                else
                    ppcs_count = ppcs_count * 3;                // 3 sec
            }
            else { // 4k res and higher
                if (core_count < CONS_CORE_COUNT)
                    ppcs_count = min_ppcs_count;
                else if (core_count < LOW_SERVER_CORE_COUNT)
                    ppcs_count = ppcs_count * 1;                // 1 sec
                else if (core_count < MED_SERVER_CORE_COUNT)
                    ppcs_count = ppcs_count * 1;                // 1 sec
                else
                    ppcs_count = ppcs_count * 3;                // 3 sec
            }
        }
        return (int32_t) ppcs_count;
    }
    else{
        SVT_ERROR("Configuration struct is corrupted\n");
        return -1;
    }
}

//return max wavefronts in a given picture
uint32_t get_max_wavefronts(uint32_t width, uint32_t height, uint32_t blk_size) {

    // possible code, needs to be tested
    // return ((height + blk_size / 2) / blk_size) < ((width  + blk_size / 2) / blk_size) ? ((height + blk_size / 2) / blk_size) : ((width  + blk_size / 2) / blk_size);
    UNUSED(width);

    return (height + blk_size / 2) / blk_size;
}
/*
* When the picture width is a single SB, must use a single segment (EncDec segments
* assume a width of at least 2 SBs)
*
* Return true if the pic width is a single SB width
*/
EbBool is_pic_width_single_sb(uint32_t sb_size, uint16_t pic_width) {
    return ((pic_width + (sb_size >> 1)) / sb_size) == 1;
}
EbErrorType load_default_buffer_configuration_settings(
    EbEncHandle        *enc_handle,
    SequenceControlSet       *scs_ptr){
    EbErrorType           return_error = EB_ErrorNone;
    unsigned int lp_count   = get_num_processors();
    unsigned int core_count = lp_count;
    uint32_t me_seg_h, me_seg_w;
#if defined(_WIN32) || defined(__linux__)
    if (scs_ptr->static_config.target_socket != -1)
        core_count /= num_groups;
#endif
    if (scs_ptr->static_config.logical_processors != 0)
        core_count = scs_ptr->static_config.logical_processors < core_count ?
            scs_ptr->static_config.logical_processors: core_count;

#ifdef _WIN32
    //Handle special case on Windows
    //by default, on Windows an application is constrained to a single group
    if (scs_ptr->static_config.target_socket == -1 &&
        scs_ptr->static_config.logical_processors == 0)
        core_count /= num_groups;

    //Affininty can only be set by group on Windows.
    //Run on both sockets if -lp is larger than logical processor per group.
    if (scs_ptr->static_config.target_socket == -1 &&
        scs_ptr->static_config.logical_processors > lp_count / num_groups)
        core_count = lp_count;
#endif
    int32_t return_ppcs = set_parent_pcs(&scs_ptr->static_config,
        core_count, scs_ptr->input_resolution);
    if (return_ppcs == -1)
        return EB_ErrorInsufficientResources;

    uint32_t enc_dec_seg_h = (core_count == SINGLE_CORE_COUNT || is_pic_width_single_sb(scs_ptr->super_block_size, scs_ptr->max_input_luma_width)) ? 1 :
        (scs_ptr->super_block_size == 128) ?
        ((scs_ptr->max_input_luma_height + 64) / 128) :
        ((scs_ptr->max_input_luma_height + 32) / 64);
    uint32_t enc_dec_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
        (scs_ptr->super_block_size == 128) ?
        ((scs_ptr->max_input_luma_width + 64) / 128) :
        ((scs_ptr->max_input_luma_width + 32) / 64);

    if (scs_ptr->static_config.rate_control_mode != 0)
    {
        me_seg_h = (core_count == SINGLE_CORE_COUNT) ? 1 :
            (((scs_ptr->max_input_luma_height + 32) / BLOCK_SIZE_64) < 6) ? 1 : 8;
        me_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
            (((scs_ptr->max_input_luma_width + 32) / BLOCK_SIZE_64) < 10) ? 1 : 6;

    }
    else {
        me_seg_h = (core_count == SINGLE_CORE_COUNT) ? 1 :
            (((scs_ptr->max_input_luma_height + 32) / BLOCK_SIZE_64) < 6) ? 1 : 2;
        me_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
            (((scs_ptr->max_input_luma_width + 32) / BLOCK_SIZE_64) < 10) ? 1 : 3;
    }
    // ME segments
    scs_ptr->me_segment_row_count_array[0] = me_seg_h;
    scs_ptr->me_segment_row_count_array[1] = me_seg_h;
    scs_ptr->me_segment_row_count_array[2] = me_seg_h;
    scs_ptr->me_segment_row_count_array[3] = me_seg_h;
    scs_ptr->me_segment_row_count_array[4] = me_seg_h;
    scs_ptr->me_segment_row_count_array[5] = me_seg_h;

    scs_ptr->me_segment_column_count_array[0] = me_seg_w;
    scs_ptr->me_segment_column_count_array[1] = me_seg_w;
    scs_ptr->me_segment_column_count_array[2] = me_seg_w;
    scs_ptr->me_segment_column_count_array[3] = me_seg_w;
    scs_ptr->me_segment_column_count_array[4] = me_seg_w;
    scs_ptr->me_segment_column_count_array[5] = me_seg_w;

    // Jing:
    // A tile group can be consisted by 1 tile or NxM tiles.
    // Segments will be parallelized within a tile group
    // We can use tile group to control the threads/parallelism in ED stage
    // NOTE:1 col will have better perf for segments for large resolutions
    //by default, do not use tile prallel. to enable, one can set one tile-group per tile.
    uint8_t tile_group_col_count = 1;
    uint8_t tile_group_row_count = 1;
    scs_ptr->tile_group_col_count_array[0] = tile_group_col_count;
    scs_ptr->tile_group_col_count_array[1] = tile_group_col_count;
    scs_ptr->tile_group_col_count_array[2] = tile_group_col_count;
    scs_ptr->tile_group_col_count_array[3] = tile_group_col_count;
    scs_ptr->tile_group_col_count_array[4] = tile_group_col_count;
    scs_ptr->tile_group_col_count_array[5] = tile_group_col_count;

    scs_ptr->tile_group_row_count_array[0] = tile_group_row_count;
    scs_ptr->tile_group_row_count_array[1] = tile_group_row_count;
    scs_ptr->tile_group_row_count_array[2] = tile_group_row_count;
    scs_ptr->tile_group_row_count_array[3] = tile_group_row_count;
    scs_ptr->tile_group_row_count_array[4] = tile_group_row_count;
    scs_ptr->tile_group_row_count_array[5] = tile_group_row_count;
    // EncDec segments
    scs_ptr->enc_dec_segment_row_count_array[0] = enc_dec_seg_h;
    scs_ptr->enc_dec_segment_row_count_array[1] = enc_dec_seg_h;
    scs_ptr->enc_dec_segment_row_count_array[2] = enc_dec_seg_h;
    scs_ptr->enc_dec_segment_row_count_array[3] = enc_dec_seg_h;
    scs_ptr->enc_dec_segment_row_count_array[4] = enc_dec_seg_h;
    scs_ptr->enc_dec_segment_row_count_array[5] = enc_dec_seg_h;

    scs_ptr->enc_dec_segment_col_count_array[0] = enc_dec_seg_w;
    scs_ptr->enc_dec_segment_col_count_array[1] = enc_dec_seg_w;
    scs_ptr->enc_dec_segment_col_count_array[2] = enc_dec_seg_w;
    scs_ptr->enc_dec_segment_col_count_array[3] = enc_dec_seg_w;
    scs_ptr->enc_dec_segment_col_count_array[4] = enc_dec_seg_w;
    scs_ptr->enc_dec_segment_col_count_array[5] = enc_dec_seg_w;

    // TPL processed in 64x64 blocks, so check width against 64x64 block size (even if SB is 128x128)
    uint32_t tpl_seg_h = (core_count == SINGLE_CORE_COUNT || is_pic_width_single_sb(64, scs_ptr->max_input_luma_width)) ? 1 :
        ((scs_ptr->max_input_luma_height + 32) / 64);

    uint32_t tpl_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
        ((scs_ptr->max_input_luma_width + 32) / 64);

    scs_ptr->tpl_segment_row_count_array = tpl_seg_h;
    scs_ptr->tpl_segment_col_count_array = tpl_seg_w;

    scs_ptr->cdef_segment_column_count = me_seg_w;
    scs_ptr->cdef_segment_row_count    = me_seg_h;

    //since restoration unit size is same for Luma and Chroma, Luma segments and chroma segments do not correspond to the same area!
    //to keep proper processing, segments have to be configured based on chroma resolution.
    uint32_t unit_size                                  = 256;
    uint32_t rest_seg_w                                 = MAX((scs_ptr->max_input_luma_width /2 + (unit_size >> 1)) / unit_size, 1);
    uint32_t rest_seg_h                                 = MAX((scs_ptr->max_input_luma_height/2 + (unit_size >> 1)) / unit_size, 1);
    scs_ptr->rest_segment_column_count =  MIN(rest_seg_w, 6);
    scs_ptr->rest_segment_row_count =  MIN(rest_seg_h, 4);

    if (scs_ptr->static_config.rate_control_mode != 0)
    {
        scs_ptr->tf_segment_column_count = 15;
        scs_ptr->tf_segment_row_count = 15;

    }
    else
    {
        scs_ptr->tf_segment_column_count = me_seg_w;//1;//
        scs_ptr->tf_segment_row_count = me_seg_h;//1;//
    }

    // adjust buffer count for superres
    uint32_t superres_count = (scs_ptr->static_config.superres_mode == SUPERRES_AUTO &&
        (scs_ptr->static_config.superres_auto_search_type == SUPERRES_AUTO_DUAL ||
         scs_ptr->static_config.superres_auto_search_type == SUPERRES_AUTO_ALL)) ? 1 : 0;

    //#====================== Data Structures and Picture Buffers ======================
    // bistream buffer will be allocated at run time. app will free the buffer once written to file.
    scs_ptr->output_stream_buffer_fifo_init_count = PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH;

    uint32_t min_input, min_parent, min_child, min_paref, min_ref, min_overlay, min_recon;
    uint32_t min_me;
    uint32_t max_input, max_parent, max_child, max_paref, max_me, max_recon;
    {
        /*Look-Ahead. Picture-Decision outputs pictures by group of mini-gops so
          the needed pictures for a certain look-ahead distance (LAD) should be rounded up to the next multiple of MiniGopSize.*/
        uint32_t mg_size = 1 << scs_ptr->static_config.hierarchical_levels;
        // uint32_t needed_lad_pictures = ((mg_size - 1) / mg_size) * mg_size; // remove because always 0
        uint32_t overlay = scs_ptr->static_config.enable_overlays ? 1 : 0;

        /*To accomodate FFMPEG EOS, 1 frame delay is needed in Resource coordination.
           note that we have the option to not add 1 frame delay of Resource Coordination. In this case we have wait for first I frame
           to be released back to be able to start first base(16). Anyway poc16 needs to wait for poc0 to finish.*/
        uint32_t eos_delay = 1;

        //Minimum input pictures needed in the pipeline
        uint16_t lad_mg_pictures = (1 + mg_size + overlay) * scs_ptr->lad_mg; //Unit= 1(provision for a potential delayI) + prediction struct + potential overlay        return_ppcs = (1 + mg_size) * (scs_ptr->lad_mg + 1)  + scs_ptr->scd_delay + eos_delay;
        return_ppcs = (1 + mg_size) * (scs_ptr->lad_mg + 1) + scs_ptr->scd_delay + eos_delay;
        //scs_ptr->input_buffer_fifo_init_count = return_ppcs;
        min_input = return_ppcs;

        if (scs_ptr->static_config.enable_overlays)
            //scs_ptr->picture_control_set_pool_init_count =
            min_parent = (mg_size + eos_delay + scs_ptr->scd_delay) * 2 + //min to get a mini-gop in PD - release all and keep one.
                         mg_size + 1 + // minigop in PM (+ 1 for key frame in the first minigop)
                         overlay;      // overlay frame of minigop in PM
        else
            min_parent = return_ppcs;

        //Pic-Manager will inject one child at a time.
        min_child = 1;

        //References. Min to sustain dec order flow (RA-5L-MRP-ON) 7 pictures from previous MGs + 11 needed for curr mini-GoP
        PredictionStructure*pred_struct_ptr = get_prediction_structure(
            enc_handle->scs_instance_array[0]->encode_context_ptr->prediction_structure_group_ptr,
            enc_handle->scs_instance_array[0]->scs_ptr->static_config.pred_structure,
            4,
            scs_ptr->static_config.hierarchical_levels);

        uint16_t num_ref_from_past_mgs = tot_past_refs[scs_ptr->static_config.hierarchical_levels];
        uint16_t num_ref_from_cur_mg = get_num_refs_in_one_mg(pred_struct_ptr) + 1 ;//+1: to accomodate one for a delayed-I

        //printf("CUR_MG_REFs:%i \n", num_ref_from_cur_mg);

        uint16_t num_ref_lad_mgs = num_ref_from_cur_mg * scs_ptr->lad_mg;
        min_ref = num_ref_from_past_mgs + num_ref_from_cur_mg + num_ref_lad_mgs;
        if (scs_ptr->static_config.pass == ENC_FIRST_PASS)
            min_me = min_parent;
        else if (scs_ptr->tpl_level) {
            // PictureDecisionContext.mg_size = mg_size + overlay; see EbPictureDecisionProcess.c line 5680
            min_me = 1 +                  // potential delay I
                     lad_mg_pictures +    // 16 + 1 ME data used in store_tpl_pictures() at line 5717
                     (mg_size + overlay); // 16 + 1 ME data used in store_tpl_pictures() at line 5729
        }
        else
            min_me = 1;

        //PA REF
        uint16_t num_pa_ref_from_past_mgs = tot_past_refs[scs_ptr->static_config.hierarchical_levels];
        //printf("TOT_PAST_REFs:%i \n", num_pa_ref_from_past_mgs);
        uint16_t num_pa_ref_from_cur_mg = mg_size; //ref+nref; nRef PA buffers are processed in PicAnalysis and used in TF
        uint16_t num_pa_ref_for_cur_mg = num_pa_ref_from_past_mgs + num_pa_ref_from_cur_mg;
        min_paref = num_pa_ref_for_cur_mg + 1 + lad_mg_pictures + scs_ptr->scd_delay + eos_delay ;
        if (scs_ptr->static_config.enable_overlays) {
            // the additional paref count for overlay is mg_size + scs_ptr->scd_delay.
            // in resource_coordination, allocate 1 additional paref for each potential overlay picture in minigop. (line 1259)
            // in picture_decision, for each minigop, keep 1 paref for the real overlay picture and release others. (line 5109)
            min_paref += mg_size + scs_ptr->scd_delay; // min_paref *= 2;
        }
        //Overlays
        min_overlay = scs_ptr->static_config.enable_overlays ?
              mg_size + eos_delay + scs_ptr->scd_delay : 1;

        min_recon = min_ref;

        if (scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
            scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B)
        {
            min_input = min_parent = 1 + scs_ptr->scd_delay + eos_delay;
            min_child = 1;
            min_ref = num_ref_from_past_mgs + num_ref_from_cur_mg;
            min_me = 1;
            min_paref = num_pa_ref_for_cur_mg + scs_ptr->scd_delay + eos_delay;
            uint32_t low_delay_tf_frames = scs_ptr->tf_params_per_type[1].max_num_past_pics;
            min_input  += low_delay_tf_frames;
            min_parent += low_delay_tf_frames;
            min_ref    += low_delay_tf_frames;
            min_me     += low_delay_tf_frames;
            min_paref  += low_delay_tf_frames;

        }
        //Configure max needed buffers to process 1+n_extra_mg Mini-Gops in the pipeline. n extra MGs to feed to picMgr on top of current one.
        uint32_t n_extra_mg;
        if (scs_ptr->static_config.rate_control_mode != 0)
        {
            if ((core_count < PARALLEL_LEVEL_4_RANGE) || (scs_ptr->input_resolution > INPUT_SIZE_8K_RANGE)) {
                n_extra_mg = 0;
            }
            else if ((core_count >= PARALLEL_LEVEL_4_RANGE) && (core_count < PARALLEL_LEVEL_8_RANGE)) {
                n_extra_mg = 1;
            }
            else if ((core_count >= PARALLEL_LEVEL_8_RANGE) && (core_count < PARALLEL_LEVEL_16_RANGE)) {
                n_extra_mg = 2;
            }
            else {
                n_extra_mg = 3;
            }
        }
        else {
            if ((core_count < PARALLEL_LEVEL_4_RANGE) || (scs_ptr->input_resolution > INPUT_SIZE_8K_RANGE)) {
                n_extra_mg = 0;
            }
            else if ((core_count >= PARALLEL_LEVEL_4_RANGE) && (core_count < PARALLEL_LEVEL_8_RANGE)) {
                n_extra_mg = 1;
            }
            else if ((core_count >= PARALLEL_LEVEL_8_RANGE) && (core_count < PARALLEL_LEVEL_16_RANGE)) {
                n_extra_mg = 2;
            }
            else if ((core_count >= PARALLEL_LEVEL_16_RANGE) && (core_count < PARALLEL_LEVEL_32_RANGE)) {
                n_extra_mg = 3;
            }
            else {
                n_extra_mg = scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 6 : scs_ptr->input_resolution <= INPUT_SIZE_8K_RANGE ? 5 : 0;
            }
        }

        max_input  = min_input + (1 + mg_size) * n_extra_mg;
        max_parent = max_input;
        max_child = (mg_size / 2) * (n_extra_mg + 1);
        max_child = MAX(max_child, 1);//have at least one child for mg_size<2

        // max_ref defines here to avoid cppcheck warning
        uint32_t max_ref = min_ref   + num_ref_from_cur_mg * n_extra_mg;
        max_paref = min_paref + (1 + mg_size)       * n_extra_mg;
        max_me    = min_me    + (1 + mg_size)       * n_extra_mg;
        max_recon = max_ref;
        // if tpl_la is disabled when super-res fix/random, input speed is much faster than recon output speed,
        // recon_output_fifo might be full and freeze at recon_output()
        if (!scs_ptr->tpl_level && scs_ptr->static_config.recon_enabled)
            max_recon = min_recon = MAX(max_ref, 30);
    }


    if (core_count == SINGLE_CORE_COUNT || MIN_PIC_PARALLELIZATION) {
        scs_ptr->input_buffer_fifo_init_count                  = min_input;
        scs_ptr->picture_control_set_pool_init_count           = min_parent;
        scs_ptr->pa_reference_picture_buffer_init_count        = min_paref;
        scs_ptr->reference_picture_buffer_init_count           = min_ref;
        scs_ptr->picture_control_set_pool_init_count_child     = min_child;
        scs_ptr->enc_dec_pool_init_count                    = min_child;
        scs_ptr->overlay_input_picture_buffer_init_count       = min_overlay;

        scs_ptr->output_recon_buffer_fifo_init_count = MAX(scs_ptr->reference_picture_buffer_init_count, min_recon);
        scs_ptr->me_pool_init_count = min_me;
    }
    else if (core_count <= PARALLEL_LEVEL_2_RANGE) {
        scs_ptr->input_buffer_fifo_init_count = clamp(max_input, min_input, max_input);
        scs_ptr->picture_control_set_pool_init_count = clamp(max_parent, min_parent, max_parent);
        scs_ptr->pa_reference_picture_buffer_init_count = clamp(max_paref, min_paref, max_paref);
        scs_ptr->output_recon_buffer_fifo_init_count = scs_ptr->reference_picture_buffer_init_count = clamp(max_recon, min_recon, max_recon);
        scs_ptr->picture_control_set_pool_init_count_child = scs_ptr->enc_dec_pool_init_count = clamp(2, min_child, max_child) + superres_count;
        scs_ptr->me_pool_init_count = clamp(max_me, min_me, max_me);
        scs_ptr->overlay_input_picture_buffer_init_count = min_overlay;
    }
    else if (core_count < PARALLEL_LEVEL_4_RANGE) {
        scs_ptr->input_buffer_fifo_init_count = clamp(max_input, min_input, max_input);
        scs_ptr->picture_control_set_pool_init_count = clamp(max_parent, min_parent, max_parent);
        scs_ptr->pa_reference_picture_buffer_init_count = clamp(max_paref, min_paref, max_paref);
        scs_ptr->output_recon_buffer_fifo_init_count = scs_ptr->reference_picture_buffer_init_count = clamp(max_recon, min_recon, max_recon);
        scs_ptr->picture_control_set_pool_init_count_child = scs_ptr->enc_dec_pool_init_count = clamp(6, min_child, max_child) + superres_count;
        scs_ptr->me_pool_init_count = clamp(max_me, min_me, max_me);
        scs_ptr->overlay_input_picture_buffer_init_count = min_overlay;
    }
    else if (core_count < PARALLEL_LEVEL_32_RANGE) {
        scs_ptr->input_buffer_fifo_init_count = clamp(max_input, min_input, max_input);
        scs_ptr->picture_control_set_pool_init_count = clamp(max_parent, min_parent, max_parent);
        scs_ptr->pa_reference_picture_buffer_init_count = clamp(max_paref, min_paref, max_paref);
        scs_ptr->output_recon_buffer_fifo_init_count = scs_ptr->reference_picture_buffer_init_count = clamp(max_recon, min_recon, max_recon);
        scs_ptr->picture_control_set_pool_init_count_child = scs_ptr->enc_dec_pool_init_count = clamp(16, min_child, max_child) + superres_count;
        scs_ptr->me_pool_init_count = clamp(max_me, min_me, max_me);
        scs_ptr->overlay_input_picture_buffer_init_count = min_overlay;

    }
    else {
        scs_ptr->input_buffer_fifo_init_count = clamp(max_input, min_input, max_input);
        scs_ptr->picture_control_set_pool_init_count = clamp(max_parent, min_parent, max_parent);
        scs_ptr->pa_reference_picture_buffer_init_count = clamp(max_paref, min_paref, max_paref);
        scs_ptr->output_recon_buffer_fifo_init_count = scs_ptr->reference_picture_buffer_init_count = clamp(max_recon, min_recon, max_recon);
        scs_ptr->picture_control_set_pool_init_count_child = scs_ptr->enc_dec_pool_init_count = clamp(max_child, min_child, max_child) + superres_count;
        scs_ptr->me_pool_init_count = clamp(max_me, min_me, max_me);
        scs_ptr->overlay_input_picture_buffer_init_count = min_overlay;
    }

    //#====================== Inter process Fifos ======================
    scs_ptr->resource_coordination_fifo_init_count       = 300;
    scs_ptr->picture_analysis_fifo_init_count            = 300;
    scs_ptr->picture_decision_fifo_init_count            = 300;
    scs_ptr->initial_rate_control_fifo_init_count        = 300;
    scs_ptr->tpl_disp_fifo_init_count                    = 300;
    scs_ptr->picture_demux_fifo_init_count               = 300;
    scs_ptr->rate_control_tasks_fifo_init_count          = 300;
    scs_ptr->rate_control_fifo_init_count                = 301;
    //Jing: Too many tiles may drain the fifo
    scs_ptr->mode_decision_configuration_fifo_init_count = 300 * (MIN(9, 1<<scs_ptr->static_config.tile_rows));
    scs_ptr->motion_estimation_fifo_init_count           = 300;
    scs_ptr->entropy_coding_fifo_init_count              = 300;
    scs_ptr->enc_dec_fifo_init_count                     = 300;
    scs_ptr->dlf_fifo_init_count                         = 300;
    scs_ptr->cdef_fifo_init_count                        = 300;
    scs_ptr->rest_fifo_init_count                        = 300;
    //#====================== Processes number ======================
    scs_ptr->total_process_init_count                    = 0;

    uint32_t max_pa_proc, max_me_proc, max_tpl_proc, max_mdc_proc, max_md_proc, max_ec_proc, max_dlf_proc, max_cdef_proc, max_rest_proc;

    max_pa_proc = max_input;
    max_me_proc = max_me * me_seg_w * me_seg_h;
    max_tpl_proc = get_max_wavefronts(scs_ptr->max_input_luma_width, scs_ptr->max_input_luma_height, 64);
    max_mdc_proc = scs_ptr->picture_control_set_pool_init_count_child;
    max_md_proc = scs_ptr->picture_control_set_pool_init_count_child * get_max_wavefronts(scs_ptr->max_input_luma_width, scs_ptr->max_input_luma_height, scs_ptr->super_block_size);
    max_ec_proc = scs_ptr->picture_control_set_pool_init_count_child;
    max_dlf_proc = scs_ptr->picture_control_set_pool_init_count_child;
    max_cdef_proc = scs_ptr->picture_control_set_pool_init_count_child * scs_ptr->cdef_segment_column_count * scs_ptr->cdef_segment_row_count;
    max_rest_proc = scs_ptr->picture_control_set_pool_init_count_child * scs_ptr->rest_segment_column_count * scs_ptr->rest_segment_row_count;

    if (core_count == SINGLE_CORE_COUNT) {
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = 1);
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count           = 1);
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count                    = 1);
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = 1);
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = 1);
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = 1);
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = 1);
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = 1);
    }
    else if (core_count <= PARALLEL_LEVEL_2_RANGE) {
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = clamp(max_pa_proc, 1, max_pa_proc));
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count           = clamp(10, 1, max_me_proc));
        scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count                    = clamp(max_tpl_proc, 1, max_tpl_proc));
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = clamp(2, scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = clamp(max_ec_proc, 1, max_ec_proc));
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = clamp(max_dlf_proc, 1, max_dlf_proc));
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = clamp(max_cdef_proc, 1, max_cdef_proc));
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = clamp(max_rest_proc, 1, max_rest_proc));
    }
    else if (core_count < PARALLEL_LEVEL_4_RANGE) {
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = clamp(max_pa_proc, 1, max_pa_proc));
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count           = clamp(15, 1, max_me_proc));
        scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count                    = clamp(max_tpl_proc, 1, max_tpl_proc));
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = clamp(10, scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = clamp(max_ec_proc, 1, max_ec_proc));
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = clamp(max_dlf_proc, 1, max_dlf_proc));
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = clamp(max_cdef_proc, 1, max_cdef_proc));
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = clamp(max_rest_proc, 1, max_rest_proc));
    }
    else if (core_count < PARALLEL_LEVEL_8_RANGE) {
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = clamp(max_pa_proc, 1, max_pa_proc));
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count           = clamp(50, 1, max_me_proc));
        scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count                    = clamp(max_tpl_proc, 1, max_tpl_proc));
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = clamp(50, scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = clamp(max_ec_proc, 1, max_ec_proc));
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = clamp(max_dlf_proc, 1, max_dlf_proc));
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = clamp(max_cdef_proc, 1, max_cdef_proc));
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = clamp(max_rest_proc, 1, max_rest_proc));
    }
    else {
        if (scs_ptr->static_config.rate_control_mode != 0)
        {
            scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count = 1);
            scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count = clamp(max_pa_proc, 1, max_pa_proc));
            scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count = clamp(50, 1, max_me_proc));
            scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count = clamp(max_tpl_proc, 1, max_tpl_proc));
            scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
            scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count = clamp(150, scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
            scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count = clamp(max_ec_proc, 1, max_ec_proc));
            scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count = clamp(max_dlf_proc, 1, max_dlf_proc));
            scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count = clamp(max_cdef_proc, 1, max_cdef_proc));
            scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count = clamp(max_rest_proc, 1, max_rest_proc));
        }
        else
        {
            if (core_count < PARALLEL_LEVEL_32_RANGE) {
                scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count = 1);
                scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count = clamp(max_pa_proc, 1, max_pa_proc));
                scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count = clamp(50, 1, max_me_proc));
                scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count = clamp(max_tpl_proc, 1, max_tpl_proc));
                scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
                scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count = clamp(200, scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
                scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count = clamp(max_ec_proc, 1, max_ec_proc));
                scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count = clamp(max_dlf_proc, 1, max_dlf_proc));
                scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count = clamp(max_cdef_proc, 1, max_cdef_proc));
                scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count = clamp(max_rest_proc, 1, max_rest_proc));
            }
            else {
                scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count = 1);
                scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count = clamp(max_pa_proc, 1, max_pa_proc));
                scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count = clamp(50, 1, max_me_proc));
                scs_ptr->total_process_init_count += (scs_ptr->tpl_disp_process_init_count = clamp(max_tpl_proc, 1, max_tpl_proc));
                scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = clamp(max_mdc_proc, 1, max_mdc_proc));
                scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count = clamp(scs_ptr->picture_control_set_pool_init_count_child,
                    scs_ptr->picture_control_set_pool_init_count_child, max_md_proc));
                scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count = clamp(max_ec_proc, 1, max_ec_proc));
                scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count = clamp(max_dlf_proc, 1, max_dlf_proc));
                scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count = clamp(max_cdef_proc, 1, max_cdef_proc));
                scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count = clamp(max_rest_proc, 1, max_rest_proc));
            }
        }
    }

    scs_ptr->total_process_init_count += 6; // single processes count
    if (scs_ptr->static_config.pass == 0 || scs_ptr->static_config.pass == 3){
        SVT_INFO("Number of logical cores available: %u\n", core_count);
        SVT_INFO("Number of PPCS %u\n", scs_ptr->picture_control_set_pool_init_count);

        /******************************************************************
        * Platform detection, limit cpu flags to hardware available CPU
        ******************************************************************/
#ifdef ARCH_X86_64
        const CPU_FLAGS cpu_flags = get_cpu_flags();
        const CPU_FLAGS cpu_flags_to_use = get_cpu_flags_to_use();
        scs_ptr->static_config.use_cpu_flags &= cpu_flags_to_use;
        SVT_INFO("[asm level on system : up to %s]\n", get_asm_level_name_str(cpu_flags));
        SVT_INFO("[asm level selected : up to %s]\n", get_asm_level_name_str(scs_ptr->static_config.use_cpu_flags));
#else
        scs_ptr->static_config.use_cpu_flags &= 0;
        SVT_INFO("[asm level on system : up to %s]\n", get_asm_level_name_str(0));
        SVT_INFO("[asm level selected : up to %s]\n", get_asm_level_name_str(scs_ptr->static_config.use_cpu_flags));
#endif
    }
    return return_error;
}
 // Rate Control
static RateControlPorts rate_control_ports[] = {
    {RATE_CONTROL_INPUT_PORT_INLME,                 0},
    {RATE_CONTROL_INPUT_PORT_PACKETIZATION,         0},
    {RATE_CONTROL_INPUT_PORT_ENTROPY_CODING,        0},
    {RATE_CONTROL_INPUT_PORT_INVALID,               0}
};
// Rate Control
static uint32_t rate_control_port_lookup(
    RateControlInputPortTypes           type,
    uint32_t                                port_type_index){
    uint32_t port_index = 0;
    uint32_t port_count = 0;

    while ((type != rate_control_ports[port_index].type) && (type != RATE_CONTROL_INPUT_PORT_INVALID))
        port_count += rate_control_ports[port_index++].count;
    return (port_count + port_type_index);
}
// Rate Control
static uint32_t rate_control_port_total_count(void){
    uint32_t port_index = 0;
    uint32_t total_count = 0;

    while (rate_control_ports[port_index].type != RATE_CONTROL_INPUT_PORT_INVALID)
        total_count += rate_control_ports[port_index++].count;
    return total_count;
}

static PicMgrPorts pic_mgr_ports[] = {
    {PIC_MGR_INPUT_PORT_SOP,            0},
    {PIC_MGR_INPUT_PORT_PACKETIZATION,  0},
    {PIC_MGR_INPUT_PORT_REST,           0},
    {PIC_MGR_INPUT_PORT_INVALID,        0}
};
static uint32_t pic_mgr_port_lookup(
    PicMgrInputPortTypes           type,
    uint32_t       port_type_index) {
    uint32_t port_index = 0;
    uint32_t port_count = 0;

    while ((type != pic_mgr_ports[port_index].type) && (type != PIC_MGR_INPUT_PORT_INVALID))
        port_count += pic_mgr_ports[port_index++].count;
    return (port_count + port_type_index);
}
static uint32_t pic_mgr_port_total_count(void) {
    uint32_t port_index = 0;
    uint32_t total_count = 0;

    while (pic_mgr_ports[port_index].type != PIC_MGR_INPUT_PORT_INVALID)
        total_count += pic_mgr_ports[port_index++].count;
    return total_count;
}

// EncDec
typedef struct {
    int32_t  type;
    uint32_t  count;
} EncDecPorts_t;
static EncDecPorts_t enc_dec_ports[] = {
    {ENCDEC_INPUT_PORT_MDC,        0},
    {ENCDEC_INPUT_PORT_ENCDEC,     0},
    {ENCDEC_INPUT_PORT_INVALID,    0}
};
static EncDecPorts_t tpl_ports[] = {
    {TPL_INPUT_PORT_SOP,     0},
    {TPL_INPUT_PORT_TPL,     0},
    {TPL_INPUT_PORT_INVALID,    0}
};
// TPL
static uint32_t tpl_port_lookup(
    int32_t  type,
    uint32_t  port_type_index)
{
    uint32_t port_index = 0;
    uint32_t port_count = 0;

    while ((type != tpl_ports[port_index].type) && (type != TPL_INPUT_PORT_INVALID))
        port_count += tpl_ports[port_index++].count;
    return (port_count + port_type_index);
}

static uint32_t tpl_port_total_count(void){
    uint32_t port_index = 0;
    uint32_t total_count = 0;

    while (tpl_ports[port_index].type != TPL_INPUT_PORT_INVALID)
        total_count += tpl_ports[port_index++].count;
    return total_count;
}

/*****************************************
 * Input Port Lookup
 *****************************************/
// EncDec
static uint32_t enc_dec_port_lookup(
    int32_t  type,
    uint32_t  port_type_index)
{
    uint32_t port_index = 0;
    uint32_t port_count = 0;

    while ((type != enc_dec_ports[port_index].type) && (type != ENCDEC_INPUT_PORT_INVALID))
        port_count += enc_dec_ports[port_index++].count;
    return (port_count + port_type_index);
}
// EncDec
static uint32_t enc_dec_port_total_count(void){
    uint32_t port_index = 0;
    uint32_t total_count = 0;

    while (enc_dec_ports[port_index].type != ENCDEC_INPUT_PORT_INVALID)
        total_count += enc_dec_ports[port_index++].count;
    return total_count;
}
/*****************************************
 * Input Port Total Count
 *****************************************/

void lib_svt_encoder_send_error_exit(
    EbPtr                    hComponent,
    uint32_t                 error_code);

static void svt_enc_handle_stop_threads(EbEncHandle *enc_handle_ptr)
{
    SequenceControlSet*  control_set_ptr = enc_handle_ptr->scs_instance_array[0]->scs_ptr;
    // Resource Coordination
    EB_DESTROY_THREAD(enc_handle_ptr->resource_coordination_thread_handle);
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->picture_analysis_thread_handle_array,control_set_ptr->picture_analysis_process_init_count);

    // Picture Decision
    EB_DESTROY_THREAD(enc_handle_ptr->picture_decision_thread_handle);

    // Motion Estimation
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->motion_estimation_thread_handle_array, control_set_ptr->motion_estimation_process_init_count);

    // Initial Rate Control
    EB_DESTROY_THREAD(enc_handle_ptr->initial_rate_control_thread_handle);

    // Source Based Oprations
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->source_based_operations_thread_handle_array, control_set_ptr->source_based_operations_process_init_count);

    // TPL dispenser ME
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->tpl_disp_thread_handle_array, control_set_ptr->tpl_disp_process_init_count);

    // Picture Manager
    EB_DESTROY_THREAD(enc_handle_ptr->picture_manager_thread_handle);

    // Rate Control
    EB_DESTROY_THREAD(enc_handle_ptr->rate_control_thread_handle);

    // Mode Decision Configuration Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->mode_decision_configuration_thread_handle_array, control_set_ptr->mode_decision_configuration_process_init_count);

    // EncDec Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->enc_dec_thread_handle_array, control_set_ptr->enc_dec_process_init_count);

    // Dlf Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->dlf_thread_handle_array, control_set_ptr->dlf_process_init_count);

    // Cdef Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->cdef_thread_handle_array, control_set_ptr->cdef_process_init_count);

    // Rest Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->rest_thread_handle_array, control_set_ptr->rest_process_init_count);

    // Entropy Coding Process
    EB_DESTROY_THREAD_ARRAY(enc_handle_ptr->entropy_coding_thread_handle_array, control_set_ptr->entropy_coding_process_init_count);

    // Packetization
    EB_DESTROY_THREAD(enc_handle_ptr->packetization_thread_handle);
}
/**********************************
* Encoder Library Handle Deonstructor
**********************************/
static void svt_enc_handle_dctor(EbPtr p)
{
    EbEncHandle *enc_handle_ptr = (EbEncHandle *)p;
    svt_enc_handle_stop_threads(enc_handle_ptr);
    EB_FREE_PTR_ARRAY(enc_handle_ptr->app_callback_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->scs_pool_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->me_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_DELETE_PTR_ARRAY(enc_handle_ptr->enc_dec_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);


    EB_DELETE_PTR_ARRAY(enc_handle_ptr->pa_reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->overlay_input_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->input_cmd_resource_ptr);
    EB_DELETE(enc_handle_ptr->input_y8b_buffer_resource_ptr);

    //all buffer_y have been redirected to y8b location that just got released.
    //to prevent releasing twice, we need to reset the buffer back to NULL
    if (enc_handle_ptr->input_buffer_resource_ptr) {
        for (uint32_t w_i = 0; w_i < enc_handle_ptr->input_buffer_resource_ptr->object_total_count; ++w_i) {
            EbObjectWrapper *wrp = enc_handle_ptr->input_buffer_resource_ptr->wrapper_ptr_pool[w_i];
            EbBufferHeaderType*obj = (EbBufferHeaderType*)wrp->object_ptr;
            EbPictureBufferDesc *desc = (EbPictureBufferDesc *)obj->p_buffer;
            desc->buffer_y = 0;
        }
    }
    EB_DELETE(enc_handle_ptr->input_buffer_resource_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->output_stream_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->output_recon_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->resource_coordination_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_analysis_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_decision_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->motion_estimation_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->initial_rate_control_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_demux_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->tpl_disp_res_srm);
    EB_DELETE(enc_handle_ptr->rate_control_tasks_resource_ptr);
    EB_DELETE(enc_handle_ptr->rate_control_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->enc_dec_tasks_resource_ptr);
    EB_DELETE(enc_handle_ptr->enc_dec_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->dlf_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->cdef_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->rest_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->entropy_coding_results_resource_ptr);

    EB_DELETE(enc_handle_ptr->resource_coordination_context_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_analysis_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->motion_estimation_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->tpl_disp_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->source_based_operations_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->mode_decision_configuration_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->enc_dec_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->dlf_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->cdef_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->rest_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->entropy_coding_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->scs_instance_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->picture_decision_context_ptr);
    EB_DELETE(enc_handle_ptr->initial_rate_control_context_ptr);
    EB_DELETE(enc_handle_ptr->picture_manager_context_ptr);
    EB_DELETE(enc_handle_ptr->rate_control_context_ptr);
    EB_DELETE(enc_handle_ptr->packetization_context_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

}

/**********************************
* Encoder Library Handle Constructor
**********************************/
static EbErrorType svt_enc_handle_ctor(
    EbEncHandle *enc_handle_ptr,
    EbComponentType * ebHandlePtr)
{
    enc_handle_ptr->dctor = svt_enc_handle_dctor;

    init_thread_management_params();

    enc_handle_ptr->encode_instance_total_count                           = EB_EncodeInstancesTotalCount;
    enc_handle_ptr->compute_segments_total_count_array                    = EB_ComputeSegmentInitCount;
    // Config Set Count
    enc_handle_ptr->scs_pool_total_count                 = EB_SequenceControlSetPoolInitCount;

    // Initialize Callbacks
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->app_callback_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC(enc_handle_ptr->app_callback_ptr_array[0], sizeof(EbCallback));
    enc_handle_ptr->app_callback_ptr_array[0]->error_handler = lib_svt_encoder_send_error_exit;
    enc_handle_ptr->app_callback_ptr_array[0]->handle = ebHandlePtr;

    // Initialize Sequence Control Set Instance Array
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->scs_instance_array, enc_handle_ptr->encode_instance_total_count);
    EB_NEW(enc_handle_ptr->scs_instance_array[0], svt_sequence_control_set_instance_ctor);
    return EB_ErrorNone;
}

EbErrorType svt_input_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

EbErrorType svt_output_recon_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

EbErrorType svt_overlay_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

EbErrorType svt_output_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);

void svt_input_buffer_header_destroyer(    EbPtr p);
void svt_output_recon_buffer_header_destroyer(    EbPtr p);
void svt_output_buffer_header_destroyer(    EbPtr p);

EbErrorType svt_input_y8b_creator(EbPtr *object_dbl_ptr, EbPtr  object_init_data_ptr);
void svt_input_y8b_destroyer(EbPtr p);

EbErrorType in_cmd_ctor(
    InputCommand *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
/*
* Input Command Constructor
*/
EbErrorType svt_input_cmd_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{

    InputCommand* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, in_cmd_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType dlf_results_ctor(
    DlfResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType dlf_results_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    DlfResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, dlf_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

/*
   TPL results ctor
*/
EbErrorType tpl_disp_results_ctor(
    TplDispResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

/*
   TPL results creator
*/
EbErrorType tpl_disp_results_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    TplDispResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, tpl_disp_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType cdef_results_ctor(
    CdefResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType cdef_results_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    CdefResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, cdef_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType rest_results_ctor(
    RestResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType rest_results_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    RestResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, rest_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

static int create_pa_ref_buf_descs(EbEncHandle *enc_handle_ptr, uint32_t instance_index)
{
        SequenceControlSet* scs_ptr = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr;
        EbPaReferenceObjectDescInitData   eb_pa_ref_obj_ect_desc_init_data_structure;
        EbPictureBufferDescInitData       ref_pic_buf_desc_init_data;
        EbPictureBufferDescInitData       quart_pic_buf_desc_init_data;
        EbPictureBufferDescInitData       sixteenth_pic_buf_desc_init_data;
        // PA Reference Picture Buffers
        // Currently, only Luma samples are needed in the PA
        ref_pic_buf_desc_init_data.max_width = scs_ptr->max_input_luma_width;
        ref_pic_buf_desc_init_data.max_height = scs_ptr->max_input_luma_height;
        ref_pic_buf_desc_init_data.bit_depth = EB_8BIT;
        ref_pic_buf_desc_init_data.color_format = EB_YUV420; //use 420 for picture analysis
        //No full-resolution pixel data is allocated for PA REF,
        // it points directly to the Luma input samples of the app data
        ref_pic_buf_desc_init_data.buffer_enable_mask = 0;


        ref_pic_buf_desc_init_data.left_padding = scs_ptr->left_padding;
        ref_pic_buf_desc_init_data.right_padding = scs_ptr->right_padding;
        ref_pic_buf_desc_init_data.top_padding = scs_ptr->top_padding;
        ref_pic_buf_desc_init_data.bot_padding = scs_ptr->bot_padding;
        ref_pic_buf_desc_init_data.split_mode = EB_FALSE;
        ref_pic_buf_desc_init_data.rest_units_per_tile = scs_ptr->rest_units_per_tile;
        ref_pic_buf_desc_init_data.mfmv                = 0;
        ref_pic_buf_desc_init_data.is_16bit_pipeline   = EB_FALSE;
        ref_pic_buf_desc_init_data.enc_mode            = scs_ptr->static_config.enc_mode;

        quart_pic_buf_desc_init_data.max_width = scs_ptr->max_input_luma_width >> 1;
        quart_pic_buf_desc_init_data.max_height = scs_ptr->max_input_luma_height >> 1;
        quart_pic_buf_desc_init_data.bit_depth = EB_8BIT;
        quart_pic_buf_desc_init_data.color_format = EB_YUV420;
        quart_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        quart_pic_buf_desc_init_data.left_padding = scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.right_padding = scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.top_padding = scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.bot_padding = scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.split_mode = EB_FALSE;
        quart_pic_buf_desc_init_data.down_sampled_filtered = (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;
        quart_pic_buf_desc_init_data.rest_units_per_tile = scs_ptr->rest_units_per_tile;
        quart_pic_buf_desc_init_data.mfmv                = 0;
        quart_pic_buf_desc_init_data.is_16bit_pipeline   = EB_FALSE;
        quart_pic_buf_desc_init_data.enc_mode            = scs_ptr->static_config.enc_mode;

        sixteenth_pic_buf_desc_init_data.max_width = scs_ptr->max_input_luma_width >> 2;
        sixteenth_pic_buf_desc_init_data.max_height = scs_ptr->max_input_luma_height >> 2;
        sixteenth_pic_buf_desc_init_data.bit_depth = EB_8BIT;
        sixteenth_pic_buf_desc_init_data.color_format = EB_YUV420;
        sixteenth_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        sixteenth_pic_buf_desc_init_data.left_padding = scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.right_padding = scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.top_padding = scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.bot_padding = scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.split_mode = EB_FALSE;
        sixteenth_pic_buf_desc_init_data.down_sampled_filtered = (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;
        sixteenth_pic_buf_desc_init_data.rest_units_per_tile = scs_ptr->rest_units_per_tile;
        sixteenth_pic_buf_desc_init_data.mfmv                = 0;
        sixteenth_pic_buf_desc_init_data.is_16bit_pipeline   = EB_FALSE;
        sixteenth_pic_buf_desc_init_data.enc_mode            = scs_ptr->static_config.enc_mode;

        eb_pa_ref_obj_ect_desc_init_data_structure.reference_picture_desc_init_data = ref_pic_buf_desc_init_data;
        eb_pa_ref_obj_ect_desc_init_data_structure.quarter_picture_desc_init_data = quart_pic_buf_desc_init_data;
        eb_pa_ref_obj_ect_desc_init_data_structure.sixteenth_picture_desc_init_data = sixteenth_pic_buf_desc_init_data;
        // Reference Picture Buffers
        EB_NEW(enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index],
            svt_system_resource_ctor,
            scs_ptr->pa_reference_picture_buffer_init_count,
            EB_PictureDecisionProcessInitCount,
            0,
            svt_pa_reference_object_creator,
            &(eb_pa_ref_obj_ect_desc_init_data_structure),
            NULL);
        // Set the SequenceControlSet Picture Pool Fifo Ptrs
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->pa_reference_picture_pool_fifo_ptr =
            svt_system_resource_get_producer_fifo(enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index], 0);
#if SRM_REPORT
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->pa_reference_picture_pool_fifo_ptr->queue_ptr->log = 0;
#endif
        return 0;
}

static int create_ref_buf_descs(EbEncHandle *enc_handle_ptr, uint32_t instance_index)
{
    EbReferenceObjectDescInitData     eb_ref_obj_ect_desc_init_data_structure;
    EbPictureBufferDescInitData       ref_pic_buf_desc_init_data;
    SequenceControlSet* scs_ptr = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr;
    EbBool is_16bit = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    // Initialize the various Picture types
    ref_pic_buf_desc_init_data.max_width = scs_ptr->max_input_luma_width;
    ref_pic_buf_desc_init_data.max_height = scs_ptr->max_input_luma_height;
    ref_pic_buf_desc_init_data.bit_depth = scs_ptr->encoder_bit_depth;
    ref_pic_buf_desc_init_data.color_format = scs_ptr->static_config.encoder_color_format;
    ref_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    ref_pic_buf_desc_init_data.rest_units_per_tile = scs_ptr->rest_units_per_tile;
    ref_pic_buf_desc_init_data.sb_total_count = scs_ptr->sb_total_count;
    uint16_t padding = scs_ptr->super_block_size + 32;

    ref_pic_buf_desc_init_data.left_padding = padding;
    ref_pic_buf_desc_init_data.right_padding = padding;
    ref_pic_buf_desc_init_data.top_padding = padding;
    ref_pic_buf_desc_init_data.bot_padding = padding;
    ref_pic_buf_desc_init_data.mfmv = scs_ptr->mfmv_enabled;
    ref_pic_buf_desc_init_data.is_16bit_pipeline = scs_ptr->is_16bit_pipeline;
    // Hsan: split_mode is set @ eb_reference_object_ctor() as both unpacked reference and packed reference are needed for a 10BIT input; unpacked reference @ MD, and packed reference @ EP

    ref_pic_buf_desc_init_data.split_mode = EB_FALSE;
    ref_pic_buf_desc_init_data.down_sampled_filtered = EB_FALSE;
    ref_pic_buf_desc_init_data.enc_mode = scs_ptr->static_config.enc_mode;
    if (is_16bit)
        ref_pic_buf_desc_init_data.bit_depth = EB_10BIT;

    eb_ref_obj_ect_desc_init_data_structure.reference_picture_desc_init_data = ref_pic_buf_desc_init_data;
    eb_ref_obj_ect_desc_init_data_structure.hbd_mode_decision =
        scs_ptr->enable_hbd_mode_decision;

    // Reference Picture Buffers
    EB_NEW(
            enc_handle_ptr->reference_picture_pool_ptr_array[instance_index],
            svt_system_resource_ctor,
            scs_ptr->reference_picture_buffer_init_count,//enc_handle_ptr->ref_pic_pool_total_count,
            EB_PictureManagerProcessInitCount,
            0,
            svt_reference_object_creator,
            &(eb_ref_obj_ect_desc_init_data_structure),
            NULL);

    enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->reference_picture_pool_fifo_ptr =
        svt_system_resource_get_producer_fifo(enc_handle_ptr->reference_picture_pool_ptr_array[instance_index], 0);

#if SRM_REPORT
    enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->reference_picture_pool_fifo_ptr->queue_ptr->log = 0;
#endif

    return 0;
}

void init_fn_ptr(void);
void svt_av1_init_wedge_masks(void);
/**********************************
* Initialize Encoder Library
**********************************/
EB_API EbErrorType svt_av1_enc_init(EbComponentType *svt_enc_component)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;
    EbEncHandle *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbErrorType return_error = EB_ErrorNone;
    uint32_t instance_index;
    uint32_t process_index;
    EbColorFormat color_format = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.encoder_color_format;
    SequenceControlSet* control_set_ptr;

    setup_common_rtcd_internal(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.use_cpu_flags);
    setup_rtcd_internal(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.use_cpu_flags);

    asm_set_convolve_asm_table();

    init_intra_dc_predictors_c_internal();

    asm_set_convolve_hbd_asm_table();

    init_intra_predictors_internal();
    EbSequenceControlSetInitData scs_init;
    scs_init.sb_size = enc_handle_ptr->scs_instance_array[0]->scs_ptr->super_block_size;

    build_blk_geom(enc_handle_ptr->scs_instance_array[0]->scs_ptr->geom_idx);

    svt_av1_init_me_luts();
    init_fn_ptr();
    svt_av1_init_wedge_masks();

    /************************************
    * Sequence Control Set
    ************************************/
    EB_NEW(enc_handle_ptr->scs_pool_ptr,
        svt_system_resource_ctor,
        enc_handle_ptr->scs_pool_total_count,
        1,
        0,
        svt_sequence_control_set_creator,
        &scs_init,
        NULL);

    /************************************
    * Picture Control Set: Parent
    ************************************/
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->me_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        // The segment Width & Height Arrays are in units of SBs, not samples
        PictureControlSetInitData input_data;

        input_data.picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        input_data.picture_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
        input_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->left_padding;
        input_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->right_padding;
        input_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->top_padding;
        input_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->bot_padding;
        input_data.color_format = color_format;
        input_data.sb_sz = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz;
        input_data.max_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_sb_depth;
        input_data.ten_bit_format = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ten_bit_format;
        input_data.compressed_ten_bit_format = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.compressed_ten_bit_format;
        input_data.enc_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enc_mode;
        input_data.speed_control = (uint8_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->speed_control_flag;
        input_data.hbd_mode_decision = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enable_hbd_mode_decision;
        input_data.film_grain_noise_level = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.film_grain_denoise_strength;
        input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.encoder_bit_depth;
        input_data.ext_block_flag = (uint8_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ext_block_flag;
        input_data.log2_tile_rows = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.tile_rows;
        input_data.log2_tile_cols = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.tile_columns;
        input_data.log2_sb_sz = (scs_init.sb_size == 128) ? 5 : 4;
        input_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->is_16bit_pipeline;
        input_data.non_m8_pad_w = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_pad_right;
        input_data.non_m8_pad_h = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_pad_bottom;

        input_data.enable_tpl_la = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->tpl_level;
        input_data.in_loop_ois = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->in_loop_ois;
        input_data.enc_dec_segment_col = (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->tpl_segment_col_count_array;
        input_data.enc_dec_segment_row = (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->tpl_segment_row_count_array;
        input_data.pass = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.pass;
        input_data.skip_frame_first_pass = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.skip_frame_first_pass;
        input_data.bypass_blk_step = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.bypass_blk_step;
        input_data.ipp_ds = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.ds;
        input_data.dist_ds = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.dist_ds;
        input_data.ipp_was_ds = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_was_ds;
        input_data.final_pass_preset = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->final_pass_preset;
        input_data.bypass_zz_check = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.bypass_zz_check;
        input_data.use8blk = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.use8blk;
        input_data.reduce_me_search = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->ipp_pass_ctrls.reduce_me_search;
        input_data.rate_control_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.rate_control_mode;
        MrpCtrls* mrp_ctrl = &(enc_handle_ptr->scs_instance_array[0]->scs_ptr->mrp_ctrls);
        input_data.ref_count_used_list0 =
            MAX(mrp_ctrl->sc_base_ref_list0_count,
                MAX(mrp_ctrl->base_ref_list0_count,
                    MAX(mrp_ctrl->sc_non_base_ref_list0_count, mrp_ctrl->non_base_ref_list0_count)));

        input_data.ref_count_used_list1 =
            MAX(mrp_ctrl->sc_base_ref_list1_count,
                MAX(mrp_ctrl->base_ref_list1_count,
                    MAX(mrp_ctrl->sc_non_base_ref_list1_count, mrp_ctrl->non_base_ref_list1_count)));

        input_data.tpl_synth_size = get_tpl_synthesizer_block_size( enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->tpl_level,
            input_data.picture_width, input_data.picture_height);
        input_data.enable_adaptive_quantization = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enable_adaptive_quantization;
        input_data.calculate_variance = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->calculate_variance;
        input_data.scene_change_detection = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.scene_change_detection ||
                                            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->vq_ctrls.sharpness_ctrls.scene_transition;
        input_data.tpl_lad_mg = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->tpl_lad_mg;
        input_data.input_resolution = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->input_resolution;

        EB_NEW(
            enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index],
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->picture_control_set_pool_init_count,//enc_handle_ptr->pcs_pool_total_count,
            1,
            0,
            picture_parent_control_set_creator,
            &input_data,
            NULL);
#if SRM_REPORT
        enc_handle_ptr->picture_parent_control_set_pool_ptr_array[0]->empty_queue->log = 0;
#endif
        EB_NEW(
            enc_handle_ptr->me_pool_ptr_array[instance_index],
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->me_pool_init_count,
            1,
            0,
            me_creator,
            &input_data,
            NULL);
#if SRM_REPORT
        enc_handle_ptr->me_pool_ptr_array[instance_index]->empty_queue->log = 0;
        dump_srm_content(enc_handle_ptr->me_pool_ptr_array[instance_index], EB_FALSE);
#endif
    }



        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->enc_dec_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            // The segment Width & Height Arrays are in units of SBs, not samples
            PictureControlSetInitData input_data;
            unsigned i;
            input_data.enc_dec_segment_col = 0;
            input_data.enc_dec_segment_row = 0;
            for (i = 0; i <= enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.hierarchical_levels; ++i) {
                input_data.enc_dec_segment_col = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_col_count_array[i] > input_data.enc_dec_segment_col ?
                    (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_col_count_array[i] :
                    input_data.enc_dec_segment_col;
                input_data.enc_dec_segment_row = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_row_count_array[i] > input_data.enc_dec_segment_row ?
                    (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_row_count_array[i] :
                    input_data.enc_dec_segment_row;
            }

            input_data.picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
            input_data.picture_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
            input_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->left_padding;
            input_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->right_padding;
            input_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->top_padding;
            input_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->bot_padding;
            input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
            input_data.film_grain_noise_level = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.film_grain_denoise_strength;
            input_data.color_format = color_format;
            input_data.sb_sz = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz;
            input_data.sb_size_pix = scs_init.sb_size;
            input_data.max_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_sb_depth;
            input_data.hbd_mode_decision = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enable_hbd_mode_decision;
            input_data.cdf_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->cdf_mode;
            input_data.mfmv = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->mfmv_enabled;
            input_data.cfg_palette = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.screen_content_mode;
            //Jing: Get tile info from parent_pcs
            PictureParentControlSet *parent_pcs = (PictureParentControlSet *)enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index]->wrapper_ptr_pool[0]->object_ptr;
            input_data.tile_row_count = parent_pcs->av1_cm->tiles_info.tile_rows;
            input_data.tile_column_count = parent_pcs->av1_cm->tiles_info.tile_cols;
            input_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->is_16bit_pipeline;
            input_data.av1_cm = parent_pcs->av1_cm;
            input_data.enc_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enc_mode;

            input_data.input_resolution = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->input_resolution;

            EB_NEW(
                enc_handle_ptr->enc_dec_pool_ptr_array[instance_index],
                svt_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_pool_init_count, //EB_PictureControlSetPoolInitCountChild,
                1,
                0,
                recon_coef_creator,
                &input_data,
                NULL);
        }




        /************************************
        * Picture Control Set: Child
        ************************************/
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            // The segment Width & Height Arrays are in units of SBs, not samples
            PictureControlSetInitData input_data;
            unsigned i;

            input_data.enc_dec_segment_col = 0;
            input_data.enc_dec_segment_row = 0;
            for (i = 0; i <= enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.hierarchical_levels; ++i) {
                input_data.enc_dec_segment_col = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_col_count_array[i] > input_data.enc_dec_segment_col ?
                    (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_col_count_array[i] :
                    input_data.enc_dec_segment_col;
                input_data.enc_dec_segment_row = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_row_count_array[i] > input_data.enc_dec_segment_row ?
                    (uint16_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_segment_row_count_array[i] :
                    input_data.enc_dec_segment_row;
            }

            input_data.init_max_block_cnt = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_block_cnt;
            input_data.picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
            input_data.picture_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
            input_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->left_padding;
            input_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->right_padding;
            input_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->top_padding;
            input_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->bot_padding;
            input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
            input_data.film_grain_noise_level = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.film_grain_denoise_strength;
            input_data.color_format = color_format;
            input_data.sb_sz = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz;
            input_data.sb_size_pix = scs_init.sb_size;
            input_data.max_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_sb_depth;
            input_data.hbd_mode_decision = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enable_hbd_mode_decision;
            input_data.cdf_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->cdf_mode;
            input_data.mfmv = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->mfmv_enabled;
            input_data.cfg_palette = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.screen_content_mode;
            //Jing: Get tile info from parent_pcs
            PictureParentControlSet *parent_pcs = (PictureParentControlSet *)enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index]->wrapper_ptr_pool[0]->object_ptr;
            input_data.tile_row_count = parent_pcs->av1_cm->tiles_info.tile_rows;
            input_data.tile_column_count = parent_pcs->av1_cm->tiles_info.tile_cols;
            input_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->is_16bit_pipeline;
            input_data.av1_cm = parent_pcs->av1_cm;
            input_data.enc_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enc_mode;
            input_data.static_config = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config;

            input_data.input_resolution = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->input_resolution;

            EB_NEW(
                enc_handle_ptr->picture_control_set_pool_ptr_array[instance_index],
                svt_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->picture_control_set_pool_init_count_child, //EB_PictureControlSetPoolInitCountChild,
                1,
                0,
                picture_control_set_creator,
                &input_data,
                NULL);
        }

    /************************************
    * Picture Buffers
    ************************************/

    // Allocate Resource Arrays
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->pa_reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->overlay_input_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    pic_mgr_ports[PIC_MGR_INPUT_PORT_SOP].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count;
    pic_mgr_ports[PIC_MGR_INPUT_PORT_PACKETIZATION].count = EB_PacketizationProcessInitCount;
    pic_mgr_ports[PIC_MGR_INPUT_PORT_REST].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count;
    // Rate Control
    rate_control_ports[RATE_CONTROL_INPUT_PORT_INLME].count = EB_PictureManagerProcessInitCount;
    rate_control_ports[RATE_CONTROL_INPUT_PORT_PACKETIZATION].count = EB_PacketizationProcessInitCount;
    rate_control_ports[RATE_CONTROL_INPUT_PORT_ENTROPY_CODING].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count;

    enc_dec_ports[ENCDEC_INPUT_PORT_MDC].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count;
    enc_dec_ports[ENCDEC_INPUT_PORT_ENCDEC].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count;
    tpl_ports[TPL_INPUT_PORT_SOP].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count;
    tpl_ports[TPL_INPUT_PORT_TPL].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_process_init_count;
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {

            // Must always allocate mem b/c don't know if restoration is on or off at this point
            // The restoration assumes only 1 tile is used, so only allocate for 1 tile... see svt_av1_alloc_restoration_struct()
            PictureControlSet *pcs = (PictureControlSet *)enc_handle_ptr->picture_control_set_pool_ptr_array[instance_index]->wrapper_ptr_pool[0]->object_ptr;
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->rest_units_per_tile = pcs->rst_info[0/*Y-plane*/].units_per_tile;
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_total_count = pcs->sb_total_count;
            create_ref_buf_descs(enc_handle_ptr, instance_index);

        create_pa_ref_buf_descs(enc_handle_ptr, instance_index);

        if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.enable_overlays) {
            // Overlay Input Picture Buffers
            EB_NEW(
                enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index],
                svt_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->overlay_input_picture_buffer_init_count,
                1,
                0,
                svt_overlay_buffer_header_creator,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr,
                svt_input_buffer_header_destroyer);
            // Set the SequenceControlSet Overlay input Picture Pool Fifo Ptrs
            enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->overlay_input_picture_pool_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index], 0);
        }
    }

    /************************************
    * System Resource Managers & Fifos
    ************************************/

    //SRM to link App to Ress-Coordination via Input commands. an Input Command holds 2 picture buffers: y8bit and rest(uv8b + yuv2b)
    EB_NEW(
        enc_handle_ptr->input_cmd_resource_ptr,
        svt_system_resource_ctor,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->resource_coordination_fifo_init_count,
        1,
        EB_ResourceCoordinationProcessInitCount,
        svt_input_cmd_creator,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr,
        NULL);
    enc_handle_ptr->input_cmd_producer_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->input_cmd_resource_ptr, 0);

    //Picture Buffer SRM to hold (uv8b + yuv2b)
    EB_NEW(
        enc_handle_ptr->input_buffer_resource_ptr,
        svt_system_resource_ctor,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->input_buffer_fifo_init_count,
        1,
        0, //1/2 SRM; no consumer FIFO
        svt_input_buffer_header_creator,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr,
        svt_input_buffer_header_destroyer);
    enc_handle_ptr->input_buffer_producer_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->input_buffer_resource_ptr, 0);

    //Picture Buffer SRM to hold y8b to be shared by Pcs->enhanced and Pa_ref
    EB_NEW(
        enc_handle_ptr->input_y8b_buffer_resource_ptr,
        svt_system_resource_ctor,
        MAX(enc_handle_ptr->scs_instance_array[0]->scs_ptr->input_buffer_fifo_init_count, enc_handle_ptr->scs_instance_array[0]->scs_ptr->pa_reference_picture_buffer_init_count),
        1,
        0, //1/2 SRM; no consumer FIFO
        svt_input_y8b_creator,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr,
        svt_input_y8b_destroyer);

#if SRM_REPORT
    enc_handle_ptr->input_y8b_buffer_resource_ptr->empty_queue->log = 1;
#endif
    enc_handle_ptr->input_y8b_buffer_producer_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->input_y8b_buffer_resource_ptr, 0);

    // EbBufferHeaderType Output Stream
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_stream_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EB_NEW(
            enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index],
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->output_stream_buffer_fifo_init_count,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->total_process_init_count,//EB_PacketizationProcessInitCount,
            1,
            svt_output_buffer_header_creator,
            &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config,
            svt_output_buffer_header_destroyer);
    }
    enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr = svt_system_resource_get_consumer_fifo(enc_handle_ptr->output_stream_buffer_resource_ptr_array[0], 0);
    if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.recon_enabled) {
        // EbBufferHeaderType Output Recon
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_recon_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            EB_NEW(
                enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index],
                svt_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->output_recon_buffer_fifo_init_count,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_process_init_count,
                1,
                svt_output_recon_buffer_header_creator,
                enc_handle_ptr->scs_instance_array[0]->scs_ptr,
                svt_output_recon_buffer_header_destroyer);
        }
        enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr = svt_system_resource_get_consumer_fifo(enc_handle_ptr->output_recon_buffer_resource_ptr_array[0], 0);
    }

    // Resource Coordination Results
    {
        ResourceCoordinationResultInitData resource_coordination_result_init_data;

        EB_NEW(
            enc_handle_ptr->resource_coordination_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->resource_coordination_fifo_init_count,
            EB_ResourceCoordinationProcessInitCount,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_process_init_count,
            resource_coordination_result_creator,
            &resource_coordination_result_init_data,
            NULL);
    }

    // Picture Analysis Results
    {
        PictureAnalysisResultInitData picture_analysis_result_init_data;

        EB_NEW(
            enc_handle_ptr->picture_analysis_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_process_init_count,
            EB_PictureDecisionProcessInitCount,
            picture_analysis_result_creator,
            &picture_analysis_result_init_data,
            NULL);
    }

    // Picture Decision Results
    {
        PictureDecisionResultInitData picture_decision_result_init_data;

        EB_NEW(
            enc_handle_ptr->picture_decision_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_decision_fifo_init_count,
            EB_PictureDecisionProcessInitCount + 2,  // 1 for rate control, another 1 for packetization when superres recoding is on
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_process_init_count,
            picture_decision_result_creator,
            &picture_decision_result_init_data,
            NULL);
    }

    // Motion Estimation Results
    {
        MotionEstimationResultsInitData motion_estimation_result_init_data;

        EB_NEW(
            enc_handle_ptr->motion_estimation_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_process_init_count,
            EB_InitialRateControlProcessInitCount,
            motion_estimation_results_creator,
            &motion_estimation_result_init_data,
            NULL);
    }


    // Initial Rate Control Results
    {
        InitialRateControlResultInitData initial_rate_control_result_init_data;

        EB_NEW(
            enc_handle_ptr->initial_rate_control_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->initial_rate_control_fifo_init_count,
            EB_InitialRateControlProcessInitCount,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count,
            initial_rate_control_results_creator,
            &initial_rate_control_result_init_data,
            NULL);
    }

    // Picture Demux Results
    {
        PictureResultInitData picture_result_init_data;
        EB_NEW(
            enc_handle_ptr->picture_demux_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_demux_fifo_init_count,
            pic_mgr_port_total_count(),
            EB_PictureManagerProcessInitCount,
            picture_results_creator,
            &picture_result_init_data,
            NULL);

    }

    // TPL dispenser Results
    {
        EntropyCodingResultsInitData tpl_disp_result_init_data;
        //TPL Dispenser tasks
        EB_NEW(
            enc_handle_ptr->tpl_disp_res_srm,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_fifo_init_count,
            tpl_port_total_count(),
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_process_init_count,
            tpl_disp_results_creator,
            &tpl_disp_result_init_data,
            NULL);
    }

    // Rate Control Tasks
    {
        RateControlTasksInitData rate_control_tasks_init_data;

        EB_NEW(
            enc_handle_ptr->rate_control_tasks_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->rate_control_tasks_fifo_init_count,
            rate_control_port_total_count(),
            EB_RateControlProcessInitCount,
            rate_control_tasks_creator,
            &rate_control_tasks_init_data,
            NULL);
    }

    // Rate Control Results
    {
        RateControlResultsInitData rate_control_result_init_data;

        EB_NEW(
            enc_handle_ptr->rate_control_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->rate_control_fifo_init_count,
            EB_RateControlProcessInitCount,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count,
            rate_control_results_creator,
            &rate_control_result_init_data,
            NULL);
    }
    // EncDec Tasks
    {
        EncDecTasksInitData mode_decision_result_init_data;
        unsigned i;

        mode_decision_result_init_data.enc_dec_segment_row_count = 0;

        for (i = 0; i <= enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.hierarchical_levels; ++i) {
            mode_decision_result_init_data.enc_dec_segment_row_count = MAX(
                mode_decision_result_init_data.enc_dec_segment_row_count,
                enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_segment_row_count_array[i]);
        }

        EB_NEW(
            enc_handle_ptr->enc_dec_tasks_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_fifo_init_count,
            enc_dec_port_total_count(),
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count,
            enc_dec_tasks_creator,
            &mode_decision_result_init_data,
            NULL);
    }

    // EncDec Results
    {
        EncDecResultsInitData enc_dec_result_init_data;

        EB_NEW(
            enc_handle_ptr->enc_dec_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_process_init_count,
            enc_dec_results_creator,
            &enc_dec_result_init_data,
            NULL);
    }

    //DLF results
    {
        EntropyCodingResultsInitData delf_result_init_data;

        EB_NEW(
            enc_handle_ptr->dlf_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_process_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_process_init_count,
            dlf_results_creator,
            &delf_result_init_data,
            NULL);
    }
    //CDEF results
    {
        EntropyCodingResultsInitData cdef_result_init_data;

        EB_NEW(
            enc_handle_ptr->cdef_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_process_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count,
            cdef_results_creator,
            &cdef_result_init_data,
            NULL);
    }
    //REST results
    {
        EntropyCodingResultsInitData rest_result_init_data;

        EB_NEW(
            enc_handle_ptr->rest_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count,
            rest_results_creator,
            &rest_result_init_data,
            NULL);
    }

    // Entropy Coding Results
    {
        EntropyCodingResultsInitData entropy_coding_results_init_data;

        EB_NEW(
            enc_handle_ptr->entropy_coding_results_resource_ptr,
            svt_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count,
            EB_PacketizationProcessInitCount,
            entropy_coding_results_creator,
            &entropy_coding_results_init_data,
            NULL);
    }


    /************************************
    * App Callbacks
    ************************************/
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index)
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->app_callback_ptr = enc_handle_ptr->app_callback_ptr_array[instance_index];
    // svt Output Buffer Fifo Ptrs
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->stream_output_fifo_ptr     = svt_system_resource_get_producer_fifo(enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index], 0);
        if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.recon_enabled)
            enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->recon_output_fifo_ptr  = svt_system_resource_get_producer_fifo(enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index], 0);
    }

    /************************************
    * Contexts
    ************************************/

    // Resource Coordination Context
    EB_NEW(
        enc_handle_ptr->resource_coordination_context_ptr,
        resource_coordination_context_ctor,
        enc_handle_ptr);

    // Picture Analysis Context
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_analysis_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_process_init_count);

    for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_analysis_process_init_count; ++process_index) {

        EB_NEW(
            enc_handle_ptr->picture_analysis_context_ptr_array[process_index],
            picture_analysis_context_ctor,
            enc_handle_ptr,
            process_index);
   }

    // Picture Decision Context
    {
        // Initialize the various Picture types
        instance_index = 0;

        EB_NEW(
            enc_handle_ptr->picture_decision_context_ptr,
            picture_decision_context_ctor,
            enc_handle_ptr);
    }

    // Motion Analysis Context
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->motion_estimation_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_process_init_count);

    for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->motion_estimation_process_init_count; ++process_index) {
        EB_NEW(
            enc_handle_ptr->motion_estimation_context_ptr_array[process_index],
            motion_estimation_context_ctor,
            enc_handle_ptr,
            process_index);
    }


        // Initial Rate Control Context
        EB_NEW(
            enc_handle_ptr->initial_rate_control_context_ptr,
            initial_rate_control_context_ctor,
            enc_handle_ptr);
        // Source Based Operations Context
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->source_based_operations_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count);

        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->source_based_operations_context_ptr_array[process_index],
                source_based_operations_context_ctor,
                enc_handle_ptr,
                tpl_port_lookup(TPL_INPUT_PORT_SOP, process_index),
                pic_mgr_port_lookup(PIC_MGR_INPUT_PORT_SOP, process_index));
        }
        // TPL dispenser
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->tpl_disp_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_process_init_count);

        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->tpl_disp_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->tpl_disp_context_ptr_array[process_index],
                tpl_disp_context_ctor,
                enc_handle_ptr,
                process_index,
                tpl_port_lookup(TPL_INPUT_PORT_TPL, process_index)
            );
        }
        // Picture Manager Context
        EB_NEW(
            enc_handle_ptr->picture_manager_context_ptr,
            picture_manager_context_ctor,
            enc_handle_ptr,
            rate_control_port_lookup(RATE_CONTROL_INPUT_PORT_INLME, 0)); //Pic-Mgr uses the first Port
        // Rate Control Context
        EB_NEW(
            enc_handle_ptr->rate_control_context_ptr,
            rate_control_context_ctor,
            enc_handle_ptr,
            EB_PictureDecisionProcessInitCount);  // me_port_index

        // Mode Decision Configuration Contexts
        {
            // Mode Decision Configuration Contexts
            EB_ALLOC_PTR_ARRAY(enc_handle_ptr->mode_decision_configuration_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count);

            for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count; ++process_index) {
                EB_NEW(
                    enc_handle_ptr->mode_decision_configuration_context_ptr_array[process_index],
                    mode_decision_configuration_context_ctor,
                    enc_handle_ptr,
                    process_index,
                    enc_dec_port_lookup(ENCDEC_INPUT_PORT_MDC, process_index));
            }
        }

        uint32_t max_picture_width = 0;
        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            if (max_picture_width < enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width)
                max_picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        }

        // EncDec Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->enc_dec_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count);
        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->enc_dec_context_ptr_array[process_index],
                enc_dec_context_ctor,
                enc_handle_ptr,
                process_index,
                enc_dec_port_lookup(ENCDEC_INPUT_PORT_ENCDEC, process_index));
        }

        // Dlf Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->dlf_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_process_init_count);

        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->dlf_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->dlf_context_ptr_array[process_index],
                dlf_context_ctor,
                enc_handle_ptr,
                process_index);
        }

        //CDEF Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->cdef_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_process_init_count);

        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->cdef_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->cdef_context_ptr_array[process_index],
                cdef_context_ctor,
                enc_handle_ptr,
                process_index);
        }
        //Rest Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->rest_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count);

        EbPictureBufferDescInitData input_data;
        input_data.enc_mode = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.enc_mode;
        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->rest_context_ptr_array[process_index],
                rest_context_ctor,
                enc_handle_ptr,
                &input_data,
                process_index,
                pic_mgr_port_lookup(PIC_MGR_INPUT_PORT_REST, process_index));
        }

        // Entropy Coding Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->entropy_coding_context_ptr_array, enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count);

        for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count; ++process_index) {
            EB_NEW(
                enc_handle_ptr->entropy_coding_context_ptr_array[process_index],
                entropy_coding_context_ctor,
                enc_handle_ptr,
                process_index,
                rate_control_port_lookup(RATE_CONTROL_INPUT_PORT_ENTROPY_CODING, process_index));
        }

    // Packetization Context
    EB_NEW(
        enc_handle_ptr->packetization_context_ptr,
        packetization_context_ctor,
        enc_handle_ptr,
        rate_control_port_lookup(RATE_CONTROL_INPUT_PORT_PACKETIZATION, 0),
        pic_mgr_port_lookup(PIC_MGR_INPUT_PORT_PACKETIZATION, 0),
        EB_PictureDecisionProcessInitCount + EB_RateControlProcessInitCount);  // me_port_index
    /************************************
    * Thread Handles
    ************************************/
    EbSvtAv1EncConfiguration   *config_ptr = &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config;
    if (config_ptr->pin_threads == 1)
        svt_set_thread_management_parameters(config_ptr);

    control_set_ptr = enc_handle_ptr->scs_instance_array[0]->scs_ptr;

    // Resource Coordination
    EB_CREATE_THREAD(enc_handle_ptr->resource_coordination_thread_handle, resource_coordination_kernel, enc_handle_ptr->resource_coordination_context_ptr);
    EB_CREATE_THREAD_ARRAY(enc_handle_ptr->picture_analysis_thread_handle_array,control_set_ptr->picture_analysis_process_init_count,
        picture_analysis_kernel,
        enc_handle_ptr->picture_analysis_context_ptr_array);

    // Picture Decision
    EB_CREATE_THREAD(enc_handle_ptr->picture_decision_thread_handle, picture_decision_kernel, enc_handle_ptr->picture_decision_context_ptr);

    // Motion Estimation
    EB_CREATE_THREAD_ARRAY(enc_handle_ptr->motion_estimation_thread_handle_array, control_set_ptr->motion_estimation_process_init_count,
        motion_estimation_kernel,
        enc_handle_ptr->motion_estimation_context_ptr_array);

        // Initial Rate Control
        EB_CREATE_THREAD(enc_handle_ptr->initial_rate_control_thread_handle, initial_rate_control_kernel, enc_handle_ptr->initial_rate_control_context_ptr);

        // Source Based Oprations
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->source_based_operations_thread_handle_array, control_set_ptr->source_based_operations_process_init_count,
            source_based_operations_kernel,
            enc_handle_ptr->source_based_operations_context_ptr_array);

        // TPL dispenser
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->tpl_disp_thread_handle_array, control_set_ptr->tpl_disp_process_init_count,
            tpl_disp_kernel,//TODOOMK
            enc_handle_ptr->tpl_disp_context_ptr_array);
        // Picture Manager
        EB_CREATE_THREAD(enc_handle_ptr->picture_manager_thread_handle, picture_manager_kernel, enc_handle_ptr->picture_manager_context_ptr);
        // Rate Control
        EB_CREATE_THREAD(enc_handle_ptr->rate_control_thread_handle, rate_control_kernel, enc_handle_ptr->rate_control_context_ptr);

        // Mode Decision Configuration Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->mode_decision_configuration_thread_handle_array, control_set_ptr->mode_decision_configuration_process_init_count,
            mode_decision_configuration_kernel,
            enc_handle_ptr->mode_decision_configuration_context_ptr_array);


        // EncDec Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->enc_dec_thread_handle_array, control_set_ptr->enc_dec_process_init_count,
            mode_decision_kernel,
            enc_handle_ptr->enc_dec_context_ptr_array);

        // Dlf Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->dlf_thread_handle_array, control_set_ptr->dlf_process_init_count,
            dlf_kernel,
            enc_handle_ptr->dlf_context_ptr_array);

        // Cdef Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->cdef_thread_handle_array, control_set_ptr->cdef_process_init_count,
            cdef_kernel,
            enc_handle_ptr->cdef_context_ptr_array);

        // Rest Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->rest_thread_handle_array, control_set_ptr->rest_process_init_count,
            rest_kernel,
            enc_handle_ptr->rest_context_ptr_array);

        // Entropy Coding Process
        EB_CREATE_THREAD_ARRAY(enc_handle_ptr->entropy_coding_thread_handle_array, control_set_ptr->entropy_coding_process_init_count,
            entropy_coding_kernel,
            enc_handle_ptr->entropy_coding_context_ptr_array);

    // Packetization
    EB_CREATE_THREAD(enc_handle_ptr->packetization_thread_handle, packetization_kernel, enc_handle_ptr->packetization_context_ptr);

    svt_print_memory_usage();

    return return_error;
}

/**********************************
* DeInitialize Encoder Library
**********************************/
EB_API EbErrorType svt_av1_enc_deinit(EbComponentType *svt_enc_component){
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;

    EbEncHandle *handle = (EbEncHandle*)svt_enc_component->p_component_private;
    if (handle) {
        svt_shutdown_process(handle->input_buffer_resource_ptr);
        svt_shutdown_process(handle->input_cmd_resource_ptr);
        svt_shutdown_process(handle->resource_coordination_results_resource_ptr);
        svt_shutdown_process(handle->picture_analysis_results_resource_ptr);
        svt_shutdown_process(handle->picture_decision_results_resource_ptr);
        svt_shutdown_process(handle->motion_estimation_results_resource_ptr);
        svt_shutdown_process(handle->initial_rate_control_results_resource_ptr);
        svt_shutdown_process(handle->picture_demux_results_resource_ptr);
        svt_shutdown_process(handle->tpl_disp_res_srm);
        svt_shutdown_process(handle->rate_control_tasks_resource_ptr);
        svt_shutdown_process(handle->rate_control_results_resource_ptr);
        svt_shutdown_process(handle->enc_dec_tasks_resource_ptr);
        svt_shutdown_process(handle->enc_dec_results_resource_ptr);
        svt_shutdown_process(handle->entropy_coding_results_resource_ptr);
        svt_shutdown_process(handle->dlf_results_resource_ptr);
        svt_shutdown_process(handle->cdef_results_resource_ptr);
        svt_shutdown_process(handle->rest_results_resource_ptr);
    }

    return EB_ErrorNone;
}

EbErrorType init_svt_av1_encoder_handle(
    EbComponentType * hComponent);
/**********************************
* GetHandle
**********************************/
EB_API EbErrorType svt_av1_enc_init_handle(
    EbComponentType** p_handle,               // Function to be called in the future for manipulating the component
    void*              p_app_data,
    EbSvtAv1EncConfiguration  *config_ptr)              // pointer passed back to the client during callbacks

{
    if(p_handle == NULL)
         return EB_ErrorBadParameter;
    svt_log_init();

    #if defined(__linux__)
        if(lp_group == NULL) {
            EB_MALLOC(lp_group, INITIAL_PROCESSOR_GROUP * sizeof(processorGroup));
        }
    #endif

    *p_handle = (EbComponentType*)malloc(sizeof(EbComponentType));
    if (*p_handle == (EbComponentType*)NULL) {
        SVT_ERROR("Component Struct Malloc Failed\n");
        return EB_ErrorInsufficientResources;
    }
    // Init Component OS objects (threads, semaphores, etc.)
    // also links the various Component control functions
    EbErrorType return_error = init_svt_av1_encoder_handle(*p_handle);

    if (return_error == EB_ErrorNone) {
        ((EbComponentType*)(*p_handle))->p_application_private = p_app_data;
        return_error = svt_av1_set_default_params(config_ptr);
    }
    if (return_error != EB_ErrorNone) {
        svt_av1_enc_deinit(*p_handle);
        free(*p_handle);
        *p_handle = NULL;
        return return_error;
    }
    svt_increase_component_count();
    return return_error;
}

/**********************************
* Encoder Componenet DeInit
**********************************/
EbErrorType svt_av1_enc_component_de_init(EbComponentType  *svt_enc_component)
{
    EbErrorType       return_error = EB_ErrorNone;

    if (svt_enc_component->p_component_private) {
        EbEncHandle* handle = (EbEncHandle*)svt_enc_component->p_component_private;
        EB_DELETE(handle);
        svt_enc_component->p_component_private = NULL;
    }
    else
        return_error = EB_ErrorUndefined;
    return return_error;
}

/**********************************
* svt_av1_enc_deinit_handle
**********************************/
EB_API EbErrorType svt_av1_enc_deinit_handle(
    EbComponentType  *svt_enc_component)
{
    if (svt_enc_component) {
        EbErrorType return_error = svt_av1_enc_component_de_init(svt_enc_component);

        free(svt_enc_component);
#if  defined(__linux__)
        EB_FREE(lp_group);
#endif
        svt_decrease_component_count();
        return return_error;
    }
    return EB_ErrorInvalidComponent;
}

// Sets the default intra period the closest possible to 1 second without breaking the minigop
static int32_t compute_default_intra_period(
    SequenceControlSet       *scs_ptr){
    int32_t intra_period               = 0;
    EbSvtAv1EncConfiguration   *config = &scs_ptr->static_config;
    int32_t fps                        = config->frame_rate < 1000 ?
                                            config->frame_rate :
                                            config->frame_rate >> 16;
    int32_t mini_gop_size              = (1 << (config->hierarchical_levels));

    intra_period                       = ((int)((fps + mini_gop_size) / mini_gop_size)*(mini_gop_size));
    intra_period                       = intra_period << 1; // use a 2-sec gop by default.
    if (config->intra_refresh_type == 1)
        intra_period -= 1;

    return intra_period;
}

/*
Calculates the default LAD value
*/
static uint32_t compute_default_look_ahead(
    EbSvtAv1EncConfiguration*   config) {
    int32_t lad;
    uint32_t mg_size = 1 << config->hierarchical_levels;

    /*To accomodate FFMPEG EOS, 1 frame delay is needed in Resource coordination.
       note that we have the option to not add 1 frame delay of Resource Coordination. In this case we have wait for first I frame
       to be released back to be able to start first base(16). Anyway poc16 needs to wait for poc0 to finish.*/
    uint32_t eos_delay = 1;
    uint32_t max_tf_delay = 6;

    if (config->rate_control_mode == 0)
        lad = (1 + mg_size) * (1 + MIN_LAD_MG) + max_tf_delay + eos_delay;
    else
        lad = (1 + mg_size) * (1 + RC_DEFAULT_LAD_MG) + max_tf_delay + eos_delay;

    lad = lad > MAX_LAD ? MAX_LAD : lad; // clip to max allowed lad
    return lad;
}
/*
Updates the LAD value
*/
static void update_look_ahead(SequenceControlSet *scs_ptr) {

    /*To accomodate FFMPEG EOS, 1 frame delay is needed in Resource coordination.
           note that we have the option to not add 1 frame delay of Resource Coordination. In this case we have wait for first I frame
           to be released back to be able to start first base(16). Anyway poc16 needs to wait for poc0 to finish.*/
    uint32_t eos_delay = 1;

    uint32_t mg_size = 1 << scs_ptr->static_config.hierarchical_levels;
    if ((int32_t) (scs_ptr->static_config.look_ahead_distance - (eos_delay + scs_ptr->scd_delay)) < (int32_t) (mg_size + 1)) {
        // Not enough pictures to form the minigop. update mg_size
        scs_ptr->static_config.look_ahead_distance = mg_size + 1 + (eos_delay + scs_ptr->scd_delay);
        SVT_WARN("Minimum lookahead distance to run %dL with TF %d is %d. Force the look_ahead_distance to be %d\n",
            scs_ptr->static_config.hierarchical_levels + 1,
            scs_ptr->static_config.enable_tf,
            scs_ptr->static_config.look_ahead_distance,
            scs_ptr->static_config.look_ahead_distance);
    }

    int32_t picture_in_future = scs_ptr->static_config.look_ahead_distance;
    // Subtract pictures used for scd_delay and eos_delay
    picture_in_future = MAX(0, (int32_t)(picture_in_future - eos_delay - scs_ptr->scd_delay));
    // Subtract pictures used for minigop formation. Unit= 1(provision for a potential delayI)
    picture_in_future = MAX(0, (int32_t)(picture_in_future - (1 + mg_size)));
    // Specify the number of mini-gops to be used in the sliding window. 0: 1 mini-gop, 1: 2 mini-gops and 3: 3 mini-gops
    scs_ptr->lad_mg = (picture_in_future + (mg_size + 1) / 2) / (mg_size + 1);
    // Since TPL is tuned for 0, 1 and 2 mini-gops, we make sure lad_mg is not smaller than tpl_lad_mg
    if (scs_ptr->lad_mg < scs_ptr->tpl_lad_mg) {
        scs_ptr->lad_mg = scs_ptr->tpl_lad_mg;
        scs_ptr->static_config.look_ahead_distance = (1 + mg_size) * (scs_ptr->lad_mg + 1) + scs_ptr->scd_delay + eos_delay;
        SVT_WARN("Lookahead distance is not long enough to get best bdrate trade off. Force the look_ahead_distance to be %d\n",
            scs_ptr->static_config.look_ahead_distance);
    }
    else if (scs_ptr->lad_mg > scs_ptr->tpl_lad_mg && (scs_ptr->static_config.rate_control_mode == 0 || scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS)) {
        scs_ptr->lad_mg = scs_ptr->tpl_lad_mg;
        scs_ptr->static_config.look_ahead_distance = (1 + mg_size) * (scs_ptr->lad_mg + 1) + scs_ptr->scd_delay + eos_delay;
        SVT_WARN("For CRF or 2PASS RC mode, the maximum needed Lookahead distance is %d. Force the look_ahead_distance to be %d\n",
            scs_ptr->static_config.look_ahead_distance,
            scs_ptr->static_config.look_ahead_distance);

    }
}
/*
 * Set tf-64x64 path params
 */
void set_tf_64x64_params(
    uint8_t use_fast_filter,
    uint8_t *use_pred_64x64_only_th, uint32_t *me_exit_th,
    uint8_t use_pred_64x64_only_th_value, uint32_t me_exit_th_value) {

    *use_pred_64x64_only_th = use_fast_filter ? use_pred_64x64_only_th_value : 0;
    *me_exit_th = use_fast_filter ? me_exit_th_value : 0;
}

/*
 * Control TF
 */
void tf_controls(SequenceControlSet* scs_ptr, uint8_t tf_level) {

    switch (tf_level)
    {
    case 0:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 0;

        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 0;

        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 0;
        break;

    case 1:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 1;
        scs_ptr->tf_params_per_type[0].num_future_pics = 16;
        scs_ptr->tf_params_per_type[0].noise_adjust_future_pics = 1;
        scs_ptr->tf_params_per_type[0].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 16);
        scs_ptr->tf_params_per_type[0].hme_me_level = 0;
        scs_ptr->tf_params_per_type[0].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].do_chroma = 1;
        scs_ptr->tf_params_per_type[0].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[0].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[0].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[0].use_medium_filter = 0;
        scs_ptr->tf_params_per_type[0].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[0].use_2tap = 0;
        scs_ptr->tf_params_per_type[0].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[0].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[0].use_fast_filter,
            &scs_ptr->tf_params_per_type[0].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[0].me_exit_th, 0, 0);
        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 1;
        scs_ptr->tf_params_per_type[1].num_past_pics = 3;
        scs_ptr->tf_params_per_type[1].num_future_pics = 6;
        scs_ptr->tf_params_per_type[1].noise_adjust_past_pics = 1;
        scs_ptr->tf_params_per_type[1].noise_adjust_future_pics = 1;
        scs_ptr->tf_params_per_type[1].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 3);
        scs_ptr->tf_params_per_type[1].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 6);
        scs_ptr->tf_params_per_type[1].hme_me_level = 0;
        scs_ptr->tf_params_per_type[1].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].do_chroma = 1;
        scs_ptr->tf_params_per_type[1].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[1].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[1].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[1].use_medium_filter = 0;
        scs_ptr->tf_params_per_type[1].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[1].use_2tap = 0;
        scs_ptr->tf_params_per_type[1].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[1].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[1].use_fast_filter,
            &scs_ptr->tf_params_per_type[1].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[1].me_exit_th, 0, 0);
        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 1;
        scs_ptr->tf_params_per_type[2].num_past_pics = 1;
        scs_ptr->tf_params_per_type[2].num_future_pics = 1;
        scs_ptr->tf_params_per_type[2].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[2].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[2].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].hme_me_level = 0;
        scs_ptr->tf_params_per_type[2].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].do_chroma = 1;
        scs_ptr->tf_params_per_type[2].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[2].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[2].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[2].use_medium_filter = 0;
        scs_ptr->tf_params_per_type[2].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[2].use_2tap = 0;
        scs_ptr->tf_params_per_type[2].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[2].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[2].use_fast_filter,
            &scs_ptr->tf_params_per_type[2].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[2].me_exit_th, 0, 0);
        break;

    case 2:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 1;
        scs_ptr->tf_params_per_type[0].num_future_pics = 16;
        scs_ptr->tf_params_per_type[0].noise_adjust_future_pics = 1;
        scs_ptr->tf_params_per_type[0].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 16);
        scs_ptr->tf_params_per_type[0].hme_me_level = 1;
        scs_ptr->tf_params_per_type[0].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].do_chroma = 1;
        scs_ptr->tf_params_per_type[0].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[0].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[0].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[0].use_medium_filter = (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE) ? 1 : 0; //case 3, M4-5
        scs_ptr->tf_params_per_type[0].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[0].use_2tap = 0;
        scs_ptr->tf_params_per_type[0].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[0].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[0].use_fast_filter,
            &scs_ptr->tf_params_per_type[0].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[0].me_exit_th, 0, 0);

        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 1;
        scs_ptr->tf_params_per_type[1].num_past_pics = 3;
        scs_ptr->tf_params_per_type[1].num_future_pics = 3;
        scs_ptr->tf_params_per_type[1].noise_adjust_past_pics = 1;
        scs_ptr->tf_params_per_type[1].noise_adjust_future_pics = 1;
        scs_ptr->tf_params_per_type[1].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 3);
        scs_ptr->tf_params_per_type[1].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 6);
        scs_ptr->tf_params_per_type[1].hme_me_level = 1;
        scs_ptr->tf_params_per_type[1].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].do_chroma = 1;
        scs_ptr->tf_params_per_type[1].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[1].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[1].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[1].use_medium_filter = (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE) ? 1 : 0; //case 3, M4-5
        scs_ptr->tf_params_per_type[1].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[1].use_2tap = 0;
        scs_ptr->tf_params_per_type[1].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[1].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[1].use_fast_filter,
            &scs_ptr->tf_params_per_type[1].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[1].me_exit_th, 0, 0);

        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 1;
        scs_ptr->tf_params_per_type[2].num_past_pics = 1;
        scs_ptr->tf_params_per_type[2].num_future_pics = 1;
        scs_ptr->tf_params_per_type[2].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[2].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[2].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].hme_me_level = 1;
        scs_ptr->tf_params_per_type[2].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].eight_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].do_chroma = 1;
        scs_ptr->tf_params_per_type[2].pred_error_32x32_th = 0;
        scs_ptr->tf_params_per_type[2].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[2].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[2].use_medium_filter = (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE) ? 1 : 0; //case 3, M4-5
        scs_ptr->tf_params_per_type[2].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[2].use_2tap = 0;
        scs_ptr->tf_params_per_type[2].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[2].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[2].use_fast_filter,
            &scs_ptr->tf_params_per_type[2].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[2].me_exit_th, 0, 0);
        break;
    case 3:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 1;
        scs_ptr->tf_params_per_type[0].num_future_pics = 8;
        scs_ptr->tf_params_per_type[0].noise_adjust_future_pics = 1;
        scs_ptr->tf_params_per_type[0].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 16);
        scs_ptr->tf_params_per_type[0].hme_me_level = 2;
        scs_ptr->tf_params_per_type[0].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[0].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[0].do_chroma = 1;
        scs_ptr->tf_params_per_type[0].pred_error_32x32_th = 20 * 32 * 32;
        scs_ptr->tf_params_per_type[0].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[0].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[0].use_medium_filter = 1;
        scs_ptr->tf_params_per_type[0].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[0].use_2tap = 0;
        scs_ptr->tf_params_per_type[0].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[0].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[0].use_fast_filter,
            &scs_ptr->tf_params_per_type[0].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[0].me_exit_th, 0, 0);
        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 1;
        scs_ptr->tf_params_per_type[1].num_past_pics = 2;
        scs_ptr->tf_params_per_type[1].num_future_pics = 2;
        scs_ptr->tf_params_per_type[1].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[1].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[1].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 3);
        scs_ptr->tf_params_per_type[1].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 6);
        scs_ptr->tf_params_per_type[1].hme_me_level = 2;
        scs_ptr->tf_params_per_type[1].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[1].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[1].do_chroma = 1;
        scs_ptr->tf_params_per_type[1].pred_error_32x32_th = 20 * 32 * 32;
        scs_ptr->tf_params_per_type[1].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[1].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[1].use_medium_filter = 1;
        scs_ptr->tf_params_per_type[1].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[1].use_2tap = 0;
        scs_ptr->tf_params_per_type[1].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[1].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[1].use_fast_filter,
            &scs_ptr->tf_params_per_type[1].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[1].me_exit_th, 0, 0);
        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 1;
        scs_ptr->tf_params_per_type[2].num_past_pics = 1;
        scs_ptr->tf_params_per_type[2].num_future_pics = 1;
        scs_ptr->tf_params_per_type[2].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[2].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[2].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels) / 2, 1);
        scs_ptr->tf_params_per_type[2].hme_me_level = 2;
        scs_ptr->tf_params_per_type[2].half_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].quarter_pel_mode = 1;
        scs_ptr->tf_params_per_type[2].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[2].do_chroma = 1;
        scs_ptr->tf_params_per_type[2].pred_error_32x32_th = 20 * 32 * 32;
        scs_ptr->tf_params_per_type[2].sub_sampling_shift = 0;
        scs_ptr->tf_params_per_type[2].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[2].use_medium_filter = 1;
        scs_ptr->tf_params_per_type[2].avoid_2d_qpel = 0;
        scs_ptr->tf_params_per_type[2].use_2tap = 0;
        scs_ptr->tf_params_per_type[2].use_intra_for_noise_est = 0;
        scs_ptr->tf_params_per_type[2].use_8bit_subpel = 0;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[2].use_fast_filter,
            &scs_ptr->tf_params_per_type[2].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[2].me_exit_th, 0, 0);
        break;

    case 4:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 1;
        scs_ptr->tf_params_per_type[0].num_future_pics = 8;
        scs_ptr->tf_params_per_type[0].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[0].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 16);
        scs_ptr->tf_params_per_type[0].hme_me_level = 2;
        scs_ptr->tf_params_per_type[0].half_pel_mode = 2;
        scs_ptr->tf_params_per_type[0].quarter_pel_mode = 3;
        scs_ptr->tf_params_per_type[0].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[0].do_chroma = 0;
        scs_ptr->tf_params_per_type[0].pred_error_32x32_th = (uint64_t)~0;
        scs_ptr->tf_params_per_type[0].sub_sampling_shift = 1;
        scs_ptr->tf_params_per_type[0].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[0].use_medium_filter = 1;
        scs_ptr->tf_params_per_type[0].avoid_2d_qpel = 1;
        scs_ptr->tf_params_per_type[0].use_2tap = 1;
        scs_ptr->tf_params_per_type[0].use_intra_for_noise_est = 1;
        scs_ptr->tf_params_per_type[0].use_8bit_subpel = 1;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[0].use_fast_filter,
            &scs_ptr->tf_params_per_type[0].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[0].me_exit_th, 35, 16 * 16);
        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 1;
        scs_ptr->tf_params_per_type[1].num_past_pics = 1;
        scs_ptr->tf_params_per_type[1].num_future_pics = 1;
        scs_ptr->tf_params_per_type[1].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[1].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[1].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 3);
        scs_ptr->tf_params_per_type[1].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 6);
        scs_ptr->tf_params_per_type[1].hme_me_level = 2;
        scs_ptr->tf_params_per_type[1].half_pel_mode = 2;
        scs_ptr->tf_params_per_type[1].quarter_pel_mode = 3;
        scs_ptr->tf_params_per_type[1].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[1].do_chroma = 0;
        scs_ptr->tf_params_per_type[1].pred_error_32x32_th = (uint64_t)~0;
        scs_ptr->tf_params_per_type[1].sub_sampling_shift = 1;
        scs_ptr->tf_params_per_type[1].use_fast_filter = 0;
        scs_ptr->tf_params_per_type[1].use_medium_filter = 1;
        scs_ptr->tf_params_per_type[1].avoid_2d_qpel = 1;
        scs_ptr->tf_params_per_type[1].use_2tap = 1;
        scs_ptr->tf_params_per_type[1].use_intra_for_noise_est = 1;
        scs_ptr->tf_params_per_type[1].use_8bit_subpel = 1;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[1].use_fast_filter,
            &scs_ptr->tf_params_per_type[1].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[1].me_exit_th, 35, 16 * 16);
        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 0;

        break;

    case 5:
        // I_SLICE TF Params
        scs_ptr->tf_params_per_type[0].enabled = 1;
        scs_ptr->tf_params_per_type[0].num_future_pics = 8;
        scs_ptr->tf_params_per_type[0].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[0].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 16);
        scs_ptr->tf_params_per_type[0].hme_me_level = 2;
        scs_ptr->tf_params_per_type[0].half_pel_mode = 2;
        scs_ptr->tf_params_per_type[0].quarter_pel_mode = 3;
        scs_ptr->tf_params_per_type[0].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[0].do_chroma = 0;
        scs_ptr->tf_params_per_type[0].pred_error_32x32_th = (uint64_t)~0;
        scs_ptr->tf_params_per_type[0].sub_sampling_shift = 1;
        scs_ptr->tf_params_per_type[0].use_fast_filter = 1;
        scs_ptr->tf_params_per_type[0].use_medium_filter = 0;
        scs_ptr->tf_params_per_type[0].avoid_2d_qpel = 1;
        scs_ptr->tf_params_per_type[0].use_2tap = 1;
        scs_ptr->tf_params_per_type[0].use_intra_for_noise_est = 1;
        scs_ptr->tf_params_per_type[0].use_8bit_subpel = 1;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[0].use_fast_filter,
            &scs_ptr->tf_params_per_type[0].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[0].me_exit_th, 35, 16 * 16);
        // BASE TF Params
        scs_ptr->tf_params_per_type[1].enabled = 1;
        scs_ptr->tf_params_per_type[1].num_past_pics = 1;
        scs_ptr->tf_params_per_type[1].num_future_pics = 1;
        scs_ptr->tf_params_per_type[1].noise_adjust_past_pics = 0;
        scs_ptr->tf_params_per_type[1].noise_adjust_future_pics = 0;
        scs_ptr->tf_params_per_type[1].max_num_past_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 3);
        scs_ptr->tf_params_per_type[1].max_num_future_pics = MIN((1 << scs_ptr->static_config.hierarchical_levels), 6);
        scs_ptr->tf_params_per_type[1].hme_me_level = 2;
        scs_ptr->tf_params_per_type[1].half_pel_mode = 2;
        scs_ptr->tf_params_per_type[1].quarter_pel_mode = 3;
        scs_ptr->tf_params_per_type[1].eight_pel_mode = 0;
        scs_ptr->tf_params_per_type[1].do_chroma = 0;
        scs_ptr->tf_params_per_type[1].pred_error_32x32_th = (uint64_t)~0;
        scs_ptr->tf_params_per_type[1].sub_sampling_shift = 1;
        scs_ptr->tf_params_per_type[1].use_fast_filter = 1;
        scs_ptr->tf_params_per_type[1].use_medium_filter = 0;
        scs_ptr->tf_params_per_type[1].avoid_2d_qpel = 1;
        scs_ptr->tf_params_per_type[1].use_2tap = 1;
        scs_ptr->tf_params_per_type[1].use_intra_for_noise_est = 1;
        scs_ptr->tf_params_per_type[1].use_8bit_subpel = 1;
        set_tf_64x64_params(
            scs_ptr->tf_params_per_type[1].use_fast_filter,
            &scs_ptr->tf_params_per_type[1].use_pred_64x64_only_th,
            &scs_ptr->tf_params_per_type[1].me_exit_th, 35, 16 * 16);
        // L1 TF Params
        scs_ptr->tf_params_per_type[2].enabled = 0;
        break;
    default:
        assert(0);
        break;
    }
    // Limit the future frames used in TF for lowdelay prediction structure
    if (scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P)
        scs_ptr->tf_params_per_type[1].max_num_future_pics = 0;

    scs_ptr->tf_params_per_type[0].use_fixed_point = 1;
    scs_ptr->tf_params_per_type[1].use_fixed_point = 1;
    scs_ptr->tf_params_per_type[2].use_fixed_point = 1;

}
/*
 * Derive tune Params; 0: use objective mode params, 1: use subjective mode params
 */
void derive_vq_params(SequenceControlSet* scs_ptr) {
    VqCtrls* vq_ctrl = &scs_ptr->vq_ctrls;

    if (scs_ptr->static_config.tune == 0) {

        // Sharpness
        vq_ctrl->sharpness_ctrls.scene_transition = 1;
        vq_ctrl->sharpness_ctrls.tf               = 1;
        vq_ctrl->sharpness_ctrls.unipred_bias     = 1;
        vq_ctrl->sharpness_ctrls.ifs              = 1;
        vq_ctrl->sharpness_ctrls.cdef             = 1;
        vq_ctrl->sharpness_ctrls.restoration      = 1;
        // Sharpness
        vq_ctrl->stability_ctrls.depth_refinement = 1;
    }
    else {

        // Sharpness
        vq_ctrl->sharpness_ctrls.scene_transition = 0;
        vq_ctrl->sharpness_ctrls.unipred_bias     = 0;
        vq_ctrl->sharpness_ctrls.ifs              = 0;
        vq_ctrl->sharpness_ctrls.cdef             = 0;
        vq_ctrl->sharpness_ctrls.restoration      = 0;
        // Sharpness
        vq_ctrl->stability_ctrls.depth_refinement = 0;
    }
}
/*
 * Derive TF Params
 */
void derive_tf_params(SequenceControlSet *scs_ptr) {

    // Do not perform TF if LD or 1 Layer or 1st pass
    EbBool do_tf = scs_ptr->static_config.enable_tf && scs_ptr->static_config.hierarchical_levels >= 1 && scs_ptr->static_config.pass != ENC_FIRST_PASS;
    uint8_t tf_level = 0;
    if (do_tf == 0) {
        tf_level = 0;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M1) {
        tf_level = 1;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M5) {
        tf_level = 2;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M8) {
        tf_level = 3;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M9) {
        tf_level = 4;
    }
    else
        tf_level = 5;
    tf_controls(scs_ptr, tf_level);
}
/*
 * Set the MRP control
 */
void set_mrp_ctrl(SequenceControlSet* scs_ptr, uint8_t mrp_level) {
    MrpCtrls* mrp_ctrl = &scs_ptr->mrp_ctrls;
    switch (mrp_level)
    {
    case 0:
        mrp_ctrl->referencing_scheme          = 0;
        mrp_ctrl->sc_base_ref_list0_count     = 1;
        mrp_ctrl->sc_base_ref_list1_count     = 1;
        mrp_ctrl->sc_non_base_ref_list0_count = 1;
        mrp_ctrl->sc_non_base_ref_list1_count = 1;
        mrp_ctrl->base_ref_list0_count        = 1;
        mrp_ctrl->base_ref_list1_count        = 1;
        mrp_ctrl->non_base_ref_list0_count    = 1;
        mrp_ctrl->non_base_ref_list1_count    = 1;
        break;

    case 1:
        mrp_ctrl->referencing_scheme          = 1;
        mrp_ctrl->sc_base_ref_list0_count     = 4;
        mrp_ctrl->sc_base_ref_list1_count     = 3;
        mrp_ctrl->sc_non_base_ref_list0_count = 4;
        mrp_ctrl->sc_non_base_ref_list1_count = 3;
        mrp_ctrl->base_ref_list0_count        = 4;
        mrp_ctrl->base_ref_list1_count        = 3;
        mrp_ctrl->non_base_ref_list0_count    = 4;
        mrp_ctrl->non_base_ref_list1_count    = 3;
        break;

    case 2:
        mrp_ctrl->referencing_scheme          = 1;
        mrp_ctrl->sc_base_ref_list0_count     = 2;
        mrp_ctrl->sc_base_ref_list1_count     = 2;
        mrp_ctrl->sc_non_base_ref_list0_count = 2;
        mrp_ctrl->sc_non_base_ref_list1_count = 2;
        mrp_ctrl->base_ref_list0_count        = 4;
        mrp_ctrl->base_ref_list1_count        = 3;
        mrp_ctrl->non_base_ref_list0_count    = 4;
        mrp_ctrl->non_base_ref_list1_count    = 3;
        break;

    case 3:
        mrp_ctrl->referencing_scheme          = 1;
        mrp_ctrl->sc_base_ref_list0_count     = 2;
        mrp_ctrl->sc_base_ref_list1_count     = 2;
        mrp_ctrl->sc_non_base_ref_list0_count = 1;
        mrp_ctrl->sc_non_base_ref_list1_count = 1;
        mrp_ctrl->base_ref_list0_count        = 4;
        mrp_ctrl->base_ref_list1_count        = 3;
        mrp_ctrl->non_base_ref_list0_count    = 4;
        mrp_ctrl->non_base_ref_list1_count    = 3;
        break;

    case 4:
        mrp_ctrl->referencing_scheme          = 0;
        mrp_ctrl->sc_base_ref_list0_count     = 2;
        mrp_ctrl->sc_base_ref_list1_count     = 2;
        mrp_ctrl->sc_non_base_ref_list0_count = 1;
        mrp_ctrl->sc_non_base_ref_list1_count = 1;
        mrp_ctrl->base_ref_list0_count        = 2;
        mrp_ctrl->base_ref_list1_count        = 2;
        mrp_ctrl->non_base_ref_list0_count    = 2;
        mrp_ctrl->non_base_ref_list1_count    = 2;
        break;

    case 5:
        mrp_ctrl->referencing_scheme          = 0;
        mrp_ctrl->sc_base_ref_list0_count     = 2;
        mrp_ctrl->sc_base_ref_list1_count     = 2;
        mrp_ctrl->sc_non_base_ref_list0_count = 1;
        mrp_ctrl->sc_non_base_ref_list1_count = 1;
        mrp_ctrl->base_ref_list0_count        = 2;
        mrp_ctrl->base_ref_list1_count        = 2;
        mrp_ctrl->non_base_ref_list0_count    = 1;
        mrp_ctrl->non_base_ref_list1_count    = 1;
        break;

    default:
        assert(0);
        break;
    }
}

void set_ipp_pass_ctrls(
    SequenceControlSet* scs_ptr,
    uint8_t ipp_pass_level) {

    IppPassControls* ipp_pass_ctrls = &scs_ptr->ipp_pass_ctrls;
    switch (ipp_pass_level) {
    case 0:
        ipp_pass_ctrls->skip_frame_first_pass = 0;
        ipp_pass_ctrls->ds = 0;
        ipp_pass_ctrls->bypass_blk_step = 0;
        ipp_pass_ctrls->dist_ds = 0;
        ipp_pass_ctrls->bypass_zz_check = 0;
        ipp_pass_ctrls->use8blk = 0;
        ipp_pass_ctrls->reduce_me_search = 0;
        break;

    case 1:
        ipp_pass_ctrls->skip_frame_first_pass = 1;
        ipp_pass_ctrls->ds = 0;
        ipp_pass_ctrls->bypass_blk_step = 0;
        ipp_pass_ctrls->dist_ds = 1;
        ipp_pass_ctrls->bypass_zz_check = 1;
        ipp_pass_ctrls->use8blk = 1;
        ipp_pass_ctrls->reduce_me_search = 1;
        break;

    case 2:
        ipp_pass_ctrls->skip_frame_first_pass = 1;
        ipp_pass_ctrls->ds = 1;
        ipp_pass_ctrls->bypass_blk_step = 1;
        ipp_pass_ctrls->dist_ds = 1;
        ipp_pass_ctrls->bypass_zz_check = 1;
        ipp_pass_ctrls->use8blk = 1;
        ipp_pass_ctrls->reduce_me_search = 1;
        break;

    default:
        assert(0);
        break;
    }
}

void set_mid_pass_ctrls(
    SequenceControlSet* scs_ptr,
    uint8_t mid_pass_level) {

    MidPassControls* mid_pass_ctrls = &scs_ptr->mid_pass_ctrls;
    switch (mid_pass_level) {

    case 0:
        mid_pass_ctrls->ds = 0;
        break;

    case 1:
        mid_pass_ctrls->ds = 1;
        break;

    default:
        assert(0);
        break;
    }
}

uint8_t get_tpl_level(int8_t enc_mode, int32_t pass, int32_t lap_enabled, uint8_t pred_structure, uint8_t superres_mode) {
    uint8_t tpl_level;

    if (pred_structure == EB_PRED_LOW_DELAY_B) {
        SVT_WARN("TPL is disabled in low delay applications.\n");
        tpl_level = 0;
    }
    else if (pass == ENC_FIRST_PASS && lap_enabled == 0) {
        tpl_level = 0;
    }
    // allow TPL with auto-dual and auto-all
    else if (superres_mode > SUPERRES_NONE && superres_mode != SUPERRES_AUTO && superres_mode != SUPERRES_QTHRESH) {
        SVT_WARN("TPL will be disabled when super resolution is enabled!\n");
        tpl_level = 0;
    }
    else if (enc_mode <= ENC_M5)
        tpl_level = 1;
    else if (enc_mode <= ENC_M6)
        tpl_level = 3;
    else if (enc_mode <= ENC_M9)
        tpl_level = 4;
    else if (enc_mode <= ENC_M10)
        tpl_level = 5;
    else
        tpl_level = 7;

    return tpl_level;
}

/*
* Set multi Pass Params
*/
void set_multi_pass_params(SequenceControlSet *scs_ptr)
{
    EbSvtAv1EncConfiguration *config = &scs_ptr->static_config;

    // Update passes
    if (scs_ptr->static_config.pass != ENC_SINGLE_PASS)
        scs_ptr->passes = (scs_ptr->static_config.rate_control_mode == 1) ? 3 : 2;
    else
        scs_ptr->passes = 1;

    if (config->rate_control_mode == 0 && config->pass == 2)
        scs_ptr->static_config.pass = 3; //use last pass for 2nd pass CRF

    switch (config->pass) {

        case ENC_SINGLE_PASS: {
            set_ipp_pass_ctrls(scs_ptr, 0);
            set_mid_pass_ctrls(scs_ptr, 0);
            scs_ptr->ipp_was_ds = 0;
            scs_ptr->final_pass_preset = config->enc_mode;
            break;
        }

        case ENC_FIRST_PASS: {
            if (config->enc_mode <= ENC_M4)
                set_ipp_pass_ctrls(scs_ptr, 0);
            else if (config->enc_mode <= ENC_M10)
                set_ipp_pass_ctrls(scs_ptr, 1);
            else
                if (config->rate_control_mode == 0)
                    set_ipp_pass_ctrls(scs_ptr, 2);
                else
                    set_ipp_pass_ctrls(scs_ptr, 1);
            set_mid_pass_ctrls(scs_ptr, 0);
            scs_ptr->ipp_was_ds = 0;
            scs_ptr->final_pass_preset = config->enc_mode;
            scs_ptr->static_config.enc_mode = MAX_ENC_PRESET;
            scs_ptr->static_config.look_ahead_distance = 0;
            scs_ptr->static_config.rate_control_mode = 0;
            scs_ptr->static_config.intra_refresh_type = SVT_AV1_KF_REFRESH;
            scs_ptr->static_config.max_bit_rate = 0;
            scs_ptr->static_config.hierarchical_levels = 0;
            break;
        }

        case ENC_MIDDLE_PASS: {
            set_ipp_pass_ctrls(scs_ptr, 0);
            if (config->enc_mode <= ENC_M10)
                set_mid_pass_ctrls(scs_ptr, 0);
            else
                set_mid_pass_ctrls(scs_ptr, 1);
            scs_ptr->ipp_was_ds = 0;
            scs_ptr->final_pass_preset = config->enc_mode;
            if (scs_ptr->final_pass_preset <= ENC_M7)
                scs_ptr->static_config.enc_mode = ENC_M11;
            else if (scs_ptr->final_pass_preset <= ENC_M9)
                scs_ptr->static_config.enc_mode = ENC_M12;
            else
                scs_ptr->static_config.enc_mode = MAX_ENC_PRESET;
            scs_ptr->static_config.rate_control_mode = 0;
            scs_ptr->static_config.qp = 43;
            scs_ptr->static_config.intra_refresh_type = SVT_AV1_KF_REFRESH;
            scs_ptr->static_config.max_bit_rate = 0;
            break;
        }

        case ENC_LAST_PASS: {
            set_ipp_pass_ctrls(scs_ptr, 0);
            // Please make sure that ipp_was_ds is ON only when ipp_ctrls->ipp_ds is ON
            if (config->enc_mode <= ENC_M10)
                scs_ptr->ipp_was_ds = 0;
            else
                scs_ptr->ipp_was_ds = config->rate_control_mode == 0 ? 1 : 0;
            scs_ptr->final_pass_preset = config->enc_mode;
            scs_ptr->static_config.intra_refresh_type = SVT_AV1_KF_REFRESH;
            break;
        }
        default: {
            assert(0);
            break;
        }
    }

    int do_downsample =
        (scs_ptr->mid_pass_ctrls.ds || scs_ptr->ipp_pass_ctrls.ds) && scs_ptr->max_input_luma_width >= 128 && scs_ptr->max_input_luma_height >= 128
        ? 1
        : 0;

    if (do_downsample) {
        // To make sure the down scaled video has width and height of multiple of 2
        scs_ptr->max_input_luma_width = (scs_ptr->max_input_luma_width >> 2) << 1;
        scs_ptr->max_input_luma_height = (scs_ptr->max_input_luma_height >> 2) << 1;
    }

    if (scs_ptr->lap_enabled) {
        scs_ptr->static_config.intra_refresh_type = SVT_AV1_KF_REFRESH;
    }

    if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS && scs_ptr->final_pass_preset > ENC_M8)
        scs_ptr->rc_stat_gen_pass_mode = 1;
    else
        scs_ptr->rc_stat_gen_pass_mode = 0;

    if (scs_ptr->static_config.recode_loop > 0 &&
        (!scs_ptr->static_config.rate_control_mode && scs_ptr->static_config.max_bit_rate == 0) &&
        (!scs_ptr->static_config.rate_control_mode || (!scs_ptr->lap_enabled && !(scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS)))) {
        // Only allow re-encoding for 2pass VBR/CBR or 1 PASS LAP, otherwise force recode_loop to DISALLOW_RECODE or 0
        scs_ptr->static_config.recode_loop = DISALLOW_RECODE;
    }
    else if (scs_ptr->static_config.recode_loop == ALLOW_RECODE_DEFAULT) {
        // capped CRF has reencde enabled for base layer frames for all presets
        if (!scs_ptr->static_config.rate_control_mode && scs_ptr->static_config.max_bit_rate)
            scs_ptr->static_config.recode_loop = ALLOW_RECODE_KFARFGF;
        else
            scs_ptr->static_config.recode_loop = scs_ptr->static_config.enc_mode <= ENC_M2 ? ALLOW_RECODE_KFARFGF : ALLOW_RECODE_KFMAXBW;
    }

    scs_ptr->encode_context_ptr->recode_loop = scs_ptr->static_config.recode_loop;
}
void set_param_based_on_input(SequenceControlSet *scs_ptr)
{
    set_multi_pass_params(
        scs_ptr);

    scs_ptr->tpl_level = get_tpl_level(scs_ptr->static_config.enc_mode, scs_ptr->static_config.pass, scs_ptr->lap_enabled, scs_ptr->static_config.pred_structure, scs_ptr->static_config.superres_mode);

    uint16_t subsampling_x = scs_ptr->subsampling_x;
    uint16_t subsampling_y = scs_ptr->subsampling_y;
    // Update picture width, and picture height
    if (scs_ptr->max_input_luma_width % MIN_BLOCK_SIZE) {
        scs_ptr->max_input_pad_right = MIN_BLOCK_SIZE - (scs_ptr->max_input_luma_width % MIN_BLOCK_SIZE);
        scs_ptr->max_input_luma_width = scs_ptr->max_input_luma_width + scs_ptr->max_input_pad_right;
    } else {
        scs_ptr->max_input_pad_right = 0;
    }

    if (scs_ptr->max_input_luma_height % MIN_BLOCK_SIZE) {
        scs_ptr->max_input_pad_bottom = MIN_BLOCK_SIZE - (scs_ptr->max_input_luma_height % MIN_BLOCK_SIZE);
        scs_ptr->max_input_luma_height = scs_ptr->max_input_luma_height + scs_ptr->max_input_pad_bottom;
    } else {
        scs_ptr->max_input_pad_bottom = 0;
    }

    scs_ptr->max_input_chroma_width = scs_ptr->max_input_luma_width >> subsampling_x;
    scs_ptr->max_input_chroma_height = scs_ptr->max_input_luma_height >> subsampling_y;

    scs_ptr->chroma_width = scs_ptr->max_input_luma_width >> subsampling_x;
    scs_ptr->chroma_height = scs_ptr->max_input_luma_height >> subsampling_y;
    scs_ptr->seq_header.max_frame_width = scs_ptr->max_input_luma_width;
    scs_ptr->seq_header.max_frame_height = scs_ptr->max_input_luma_height;
    scs_ptr->static_config.source_width = scs_ptr->max_input_luma_width;
    scs_ptr->static_config.source_height = scs_ptr->max_input_luma_height;
    if (scs_ptr->static_config.superres_mode == SUPERRES_FIXED &&
        scs_ptr->static_config.superres_denom == SCALE_NUMERATOR &&
        scs_ptr->static_config.superres_kf_denom == SCALE_NUMERATOR) {
        scs_ptr->static_config.superres_mode = SUPERRES_NONE;
    }
    if (scs_ptr->static_config.superres_mode == SUPERRES_QTHRESH &&
        scs_ptr->static_config.superres_qthres == MAX_QP_VALUE &&
        scs_ptr->static_config.superres_kf_qthres == MAX_QP_VALUE) {
        scs_ptr->static_config.superres_mode = SUPERRES_NONE;
    }
    if (scs_ptr->static_config.superres_mode > SUPERRES_NONE) {
        if (scs_ptr->static_config.tile_rows || scs_ptr->static_config.tile_columns) {
            // disable tiles if super-res is on
            SVT_WARN("Tiles will be disabled when super resolution is enabled!\n");
            scs_ptr->static_config.tile_rows = 0;
            scs_ptr->static_config.tile_columns = 0;
        }
        if (scs_ptr->static_config.superres_mode == SUPERRES_RANDOM) {
            SVT_WARN("Super resolution random mode is designed for test and debugging purpose,\n"
                "it creates array of picture buffers for all scaling denominators (up to 8) of each reference frame.\n"
                "This mode retains a significant amount of memory, much more than other modes!\n");
        }
    }
#if FIX_1PVBR
    // Set initial qp for single pass vbr
    if ((scs_ptr->static_config.rate_control_mode) && (scs_ptr->static_config.pass == ENC_SINGLE_PASS)){
        if (scs_ptr->static_config.qp != DEFAULT_QP) {
            SVT_WARN("The input q vlaue is ignored in vbr mode %d", scs_ptr->static_config.qp);
        }
        const uint8_t tbr_bands[6] = { 1,2,4,6,8,10 };
        const uint64_t src_samples = (uint64_t)(scs_ptr->seq_header.max_frame_width*scs_ptr->seq_header.max_frame_height);
        const uint64_t target_bit_rate = scs_ptr->static_config.target_bit_rate / src_samples;
        if (target_bit_rate < tbr_bands[0])
            scs_ptr->static_config.qp = 50;
        else if (target_bit_rate < tbr_bands[1])
            scs_ptr->static_config.qp = 45;
        else if (target_bit_rate < tbr_bands[2])
            scs_ptr->static_config.qp = 40;
        else if (target_bit_rate < tbr_bands[3])
            scs_ptr->static_config.qp = 35;
        else if (target_bit_rate < tbr_bands[4])
            scs_ptr->static_config.qp = 30;
        else if (target_bit_rate < tbr_bands[5])
            scs_ptr->static_config.qp = 25;
        else
            scs_ptr->static_config.qp = 20;
    }
#endif
    derive_input_resolution(
        &scs_ptr->input_resolution,
        scs_ptr->seq_header.max_frame_width*scs_ptr->seq_header.max_frame_height);
    // Set tune params
    derive_vq_params(scs_ptr);

    // Set TF level
    derive_tf_params(scs_ptr);

    //Future frames window in Scene Change Detection (SCD) / TemporalFiltering
    scs_ptr->scd_delay = 0;

    // Update the scd_delay based on the the number of future frames @ ISLICE
    // This case is needed for non-delayed Intra (intra_period_length == 0)
    uint32_t scd_delay_islice = 0;
    if (scs_ptr->static_config.intra_period_length == 0)
        if (scs_ptr->tf_params_per_type[0].enabled)
            scd_delay_islice =
            MIN(scs_ptr->tf_params_per_type[0].num_future_pics + (scs_ptr->tf_params_per_type[0].noise_adjust_future_pics ? 3 : 0), // number of future picture(s) used for ISLICE + max picture(s) after noise-based adjustement (=3)
                scs_ptr->tf_params_per_type[0].max_num_future_pics);


    // Update the scd_delay based on the the number of future frames @ BASE
    uint32_t scd_delay_base = 0;
    if (scs_ptr->tf_params_per_type[1].enabled)
        scd_delay_base =
        MIN(scs_ptr->tf_params_per_type[1].num_future_pics + (scs_ptr->tf_params_per_type[1].noise_adjust_future_pics ? 3 : 0), // number of future picture(s) used for BASE + max picture(s) after noise-based adjustement (=3)
            scs_ptr->tf_params_per_type[1].max_num_future_pics);
    scs_ptr->scd_delay = MAX(scd_delay_islice, scd_delay_base);

    // Update the scd_delay based on SCD, 1first pass
    // Delay needed for SCD , 1first pass of (2pass and 1pass VBR)
    if (scs_ptr->static_config.scene_change_detection || scs_ptr->vq_ctrls.sharpness_ctrls.scene_transition || scs_ptr->static_config.pass == ENC_FIRST_PASS || scs_ptr->lap_enabled)
        scs_ptr->scd_delay = MAX(scs_ptr->scd_delay, 2);

    // no future minigop is used for lowdelay prediction structure
    if (scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P || scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) {
        scs_ptr->lad_mg = scs_ptr->tpl_lad_mg = 0;
    }
    else
     {
        uint8_t tpl_lad_mg = 1; // Specify the number of mini-gops to be used as LAD. 0: 1 mini-gop, 1: 2 mini-gops and 3: 3 mini-gops
        uint32_t mg_size = 1 << scs_ptr->static_config.hierarchical_levels;
        if (scs_ptr->static_config.look_ahead_distance < mg_size)
            tpl_lad_mg = 0;
        else
            if (scs_ptr->static_config.enc_mode <= ENC_M11 && scs_ptr->tpl_level != 0)
                if (scs_ptr->static_config.fast_decode <= 1)
                    tpl_lad_mg = 1;
                else if (scs_ptr->static_config.fast_decode <= 2)
                    tpl_lad_mg = scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? 1 : 0;
                else
                    tpl_lad_mg = scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 0;
            else
                tpl_lad_mg = 0;

        // special conditions for higher resolutions in order to decrease memory usage for tpl_lad_mg
#if !REMOVE_LP1_LPN_DIFF
        if (scs_ptr->static_config.logical_processors == 1 && scs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE && scs_ptr->static_config.hierarchical_levels >= 4) {
            tpl_lad_mg = 0;
        }
        else
#endif
       if (scs_ptr->input_resolution >= INPUT_SIZE_8K_RANGE && scs_ptr->static_config.hierarchical_levels >= 4) {
            tpl_lad_mg = 0;
        }
        scs_ptr->tpl_lad_mg = MIN(2, tpl_lad_mg);// lad_mg is capped to 2 because tpl was optimised only for 1,2 and 3 mini-gops
        if (scs_ptr->static_config.rate_control_mode == 0)
            scs_ptr->lad_mg = scs_ptr->tpl_lad_mg;
        else
            // update the look ahead size
            update_look_ahead(scs_ptr);
    }
    // In two pass encoding, the first pass uses sb size=64. Also when tpl is used
    // in 240P resolution, sb size is set to 64
    if (scs_ptr->static_config.pass == ENC_FIRST_PASS ||
        (scs_ptr->tpl_level && scs_ptr->input_resolution == INPUT_SIZE_240p_RANGE))
        scs_ptr->super_block_size = 64;
    else
        if (scs_ptr->static_config.enc_mode <= ENC_M2)
            scs_ptr->super_block_size = 128;
        else
            scs_ptr->super_block_size = 64;
    if (scs_ptr->static_config.rate_control_mode && !(scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS) && !scs_ptr->lap_enabled)
        scs_ptr->super_block_size = 64;

    // scs_ptr->static_config.hierarchical_levels = (scs_ptr->static_config.rate_control_mode > 1) ? 3 : scs_ptr->static_config.hierarchical_levels;
    if (scs_ptr->static_config.restricted_motion_vector && scs_ptr->super_block_size == 128) {
        scs_ptr->static_config.restricted_motion_vector = EB_FALSE;
        SVT_WARN("Restricted_motion_vector and SB 128x128 not supoorted, setting rmv to false\n");
    }
    if (scs_ptr->static_config.intra_refresh_type == SVT_AV1_FWDKF_REFRESH && scs_ptr->static_config.hierarchical_levels != 4){
        scs_ptr->static_config.hierarchical_levels = 4;
        SVT_WARN("Fwd key frame is only supported for hierarchical levels 4 at this point. Hierarchical levels are set to 4\n");
    }
    uint8_t disallow_nsq = get_disallow_nsq(scs_ptr->static_config.enc_mode);
    uint8_t disallow_4x4 = 1;
    for (EB_SLICE slice_type = 0; slice_type < IDR_SLICE + 1; slice_type++)
        disallow_4x4 = MIN(disallow_4x4, get_disallow_4x4(scs_ptr->static_config.enc_mode, slice_type));

    if (scs_ptr->super_block_size == 128) {
        scs_ptr->geom_idx = GEOM_2;
        scs_ptr->max_block_cnt = 4421;
    }
    else {
        //SB 64x64
        if (disallow_nsq && disallow_4x4) {
            scs_ptr->geom_idx = GEOM_0;
            scs_ptr->max_block_cnt = 85;
        }
        else {
            scs_ptr->geom_idx = GEOM_1;
            scs_ptr->max_block_cnt = 1101;
        }
    }
    //printf("\n\nGEOM:%i \n", scs_ptr->geom_idx);
    // Configure the padding
    scs_ptr->left_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->top_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->right_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->bot_padding = scs_ptr->super_block_size + 4;

    //for 10bit,  increase the pad of source from 68 to 72 (mutliple of 8) to accomodate 2bit-compression flow
    //we actually need to change the horizontal dimension only, but for simplicity/uniformity we do all directions
   // if (scs_ptr->static_config.encoder_bit_depth != EB_8BIT)
    {
        scs_ptr->left_padding  += 4;
        scs_ptr->top_padding   += 4;
        scs_ptr->right_padding += 4;
        scs_ptr->bot_padding   += 4;
    }


    scs_ptr->static_config.enable_overlays = !scs_ptr->static_config.enable_tf ||
        (scs_ptr->static_config.rate_control_mode > 0) ?
        0 : scs_ptr->static_config.enable_overlays;
    //0: ON
    //1: OFF
    // Memory Footprint reduction tool ONLY if no CDF (should be controlled using an API signal and not f(enc_mode))
    scs_ptr->cdf_mode = 0;
    // Set down-sampling method     Settings
    // 0                            0: filtering
    // 1                            1: decimation

    scs_ptr->down_sampling_method_me_search = ME_FILTERED_DOWNSAMPLED;

    // Enforce starting frame in decode order (at PicMgr)
    // Does not wait for feedback from PKT
    if (scs_ptr->static_config.logical_processors == 1 && // LP1
        scs_ptr->tpl_level)
        scs_ptr->enable_pic_mgr_dec_order = 1;
    else
        scs_ptr->enable_pic_mgr_dec_order = 0;
    // Enforce encoding frame in decode order
    // Wait for feedback from PKT
#if RC_NO_R2R
    scs_ptr->enable_dec_order = 1;
#else
    if (scs_ptr->static_config.logical_processors == 1 && // LP1
        ((scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS) || scs_ptr->lap_enabled))
        scs_ptr->enable_dec_order = 1;
    else
        scs_ptr->enable_dec_order = 0;
#endif
   // Open loop intra done with TPL, data is not stored
    scs_ptr->in_loop_ois = 1;

    // 1: Use boundary pixels in restoration filter search.
    // 0: Do not use boundary pixels in the restoration filter search.
    scs_ptr->use_boundaries_in_rest_search = 0;

    // Set over_boundary_block_mode     Settings
    // 0                            0: not allowed
    // 1                            1: allowed
    if (scs_ptr->over_bndry_blk == DEFAULT)
        scs_ptr->over_boundary_block_mode = 1;
    else
        scs_ptr->over_boundary_block_mode = scs_ptr->over_bndry_blk;
    if (scs_ptr->static_config.pass == ENC_FIRST_PASS)
        scs_ptr->over_boundary_block_mode = 0;

    if (scs_ptr->static_config.enable_mfmv == DEFAULT)
        if (scs_ptr->static_config.enc_mode <= ENC_M5)
            scs_ptr->mfmv_enabled = 1;
        else if(scs_ptr->static_config.enc_mode <= ENC_M10)
            if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                scs_ptr->mfmv_enabled = 1;
            else
                scs_ptr->mfmv_enabled = 0;
        else
            scs_ptr->mfmv_enabled = 0;
    else
        scs_ptr->mfmv_enabled = scs_ptr->static_config.enable_mfmv;
    // Set hbd_mode_decision OFF for high encode modes or bitdepth < 10
    if (scs_ptr->static_config.encoder_bit_depth < 10)
        scs_ptr->enable_hbd_mode_decision = 0;

    // Throws a warning when scene change is on, as the feature is not optimal and may produce false detections
    if (scs_ptr->static_config.scene_change_detection == 1)
        SVT_WARN("Scene Change is not optimal and may produce suboptimal keyframe placements\n");

    // MRP level
    uint8_t mrp_level;

    if (scs_ptr->static_config.enc_mode <= ENC_MRS) {

        mrp_level = 1;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M1) {

        mrp_level = 2;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M3) {

        mrp_level = 3;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M6) {

        mrp_level = 4;
    }
    else if (scs_ptr->static_config.enc_mode <= ENC_M11) {

        mrp_level = 5;
    }
    else {
        mrp_level = 0;
    }

    set_mrp_ctrl(scs_ptr, mrp_level);
    scs_ptr->is_short_clip = 0; // set to 1 if multipass and less than 200 frames in resourcecordination

    // Varaince is required for scene change detection and segmentation-based quantization
    if (scs_ptr->static_config.enable_adaptive_quantization == 1 ||
        scs_ptr->static_config.scene_change_detection == 1)
        scs_ptr->calculate_variance = 1;
    else if (scs_ptr->static_config.enc_mode <= ENC_M10)
        scs_ptr->calculate_variance = 1;
    else
        scs_ptr->calculate_variance = 0;
}
/******************************************************
 * Read Stat from File
 ******************************************************/
extern void read_stat(SequenceControlSet *scs_ptr);

extern void setup_two_pass(SequenceControlSet *scs_ptr);
void set_mini_gop_size_controls(MiniGopSizeCtrls *mgs_ctls, uint8_t mg_level,int input_resolution) {
    switch (mg_level) {
    case 0:
        mgs_ctls->adptive_enable = 0;
        break;
    case 1:
        mgs_ctls->adptive_enable = 1;
        mgs_ctls->animation_type_th = 0.40;
        mgs_ctls->hm_th = 0.95;
        mgs_ctls->hsa_th = 0.5;
        mgs_ctls->lfr_th = 50;
        mgs_ctls->lm_th = 0.0001;
        mgs_ctls->short_shot_th = 3;
        mgs_ctls->hmv_di_th = 0.75;
        mgs_ctls->lmv_di_th = input_resolution < INPUT_SIZE_360p_RANGE ? 0.5 : 0.35;
        break;
    default:
        mgs_ctls->adptive_enable = 0;
    }
}
void set_max_mini_gop_size(SequenceControlSet *scs_ptr, MiniGopSizeCtrls *mgs_ctls) {
    if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS) {
        read_stat(scs_ptr);
        setup_two_pass(scs_ptr);
        const double resolution_offset[2][INPUT_SIZE_COUNT] = { { 0.3,0.1,0.0,0.0,0.0,0.0,0.0 },{ 0.37,0.12,0.05,0.05,0.5,0.5,0.0 } };
        FIRSTPASS_STATS * stat = scs_ptr->twopass.stats_buf_ctx->total_stats;
        double low_motion_clip = (stat->pcnt_inter - stat->pcnt_motion) / (stat->count - 1);
        double content_type = (stat->intra_skip_pct / (stat->count - 1)) >= mgs_ctls->animation_type_th ? FC_GRAPHICS_ANIMATION : FC_NORMAL;
        // Avoid long gop for animation contents with low static area.
        double avoid_long_gop = content_type == FC_GRAPHICS_ANIMATION ? 1 : 0;
        // Avoid long_gop for short clips
        avoid_long_gop = stat->count < (mgs_ctls->short_shot_th * 32) ? 1 : avoid_long_gop;
        EbInputResolution input_resolution;
        if (scs_ptr->mid_pass_ctrls.ds)
            derive_input_resolution(
                &input_resolution,
                (scs_ptr->seq_header.max_frame_width << 1)*(scs_ptr->seq_header.max_frame_height << 1));
        else
            derive_input_resolution(
                &input_resolution,
                scs_ptr->seq_header.max_frame_width*scs_ptr->seq_header.max_frame_height);
        double lm_th = (0.6 + resolution_offset[scs_ptr->ipp_was_ds][input_resolution]);
        uint32_t fps = (uint32_t)((scs_ptr->static_config.frame_rate > 1000) ?
            scs_ptr->static_config.frame_rate >> 16 :
            scs_ptr->static_config.frame_rate);
        double short_shot = (stat->count < (mgs_ctls->short_shot_th * 32)) ? 1 : 0;
        double unid_motion = ((stat->mv_in_out_count / (stat->count - 1)) > mgs_ctls->lmv_di_th) && ((stat->mv_in_out_count / (stat->count - 1)) < mgs_ctls->hmv_di_th) ? 1 : 0;
        double low_frame_rate = (fps < mgs_ctls->lfr_th) ? 1 : 0;
        double lm_th_offset = (short_shot && unid_motion && low_frame_rate ? (scs_ptr->ipp_was_ds ? 0.2 : 0.065) : 0.0);
        double hm_th = (mgs_ctls->lm_th + lm_th_offset);
        uint32_t  min_gop_size = ((low_motion_clip > lm_th) && !avoid_long_gop) ? 32 : low_motion_clip > hm_th ? 16 : 8;
        double high_motion_clip = (stat->pcnt_motion) / (stat->count - 1);
        double mv_in_out_count = ABS(stat->mv_in_out_count / (stat->count - 1));
        min_gop_size = (high_motion_clip > 0.78) && (mv_in_out_count > (scs_ptr->ipp_pass_ctrls.ds ? 0.5 : 0.6)) ? min_gop_size / 2 : min_gop_size;
        min_gop_size = MIN(32, MAX(8, min_gop_size));

        switch (min_gop_size) {
        case 1:
            scs_ptr->static_config.hierarchical_levels = 0;
            break;
        case 2:
            scs_ptr->static_config.hierarchical_levels = 1;
            break;
        case 4:
            scs_ptr->static_config.hierarchical_levels = 2;
            break;
        case 8:
            scs_ptr->static_config.hierarchical_levels = 3;
            break;
        case 16:
            scs_ptr->static_config.hierarchical_levels = 4;
            break;
        case 32:
            scs_ptr->static_config.hierarchical_levels = 5;
            break;
        default:
            scs_ptr->static_config.hierarchical_levels = 4;
            break;
        }
        scs_ptr->max_temporal_layers = scs_ptr->static_config.hierarchical_levels;


    }
    scs_ptr->enable_adaptive_mini_gop = mgs_ctls->adptive_enable;
    scs_ptr->max_heirachical_level = scs_ptr->static_config.hierarchical_levels;
}
void copy_api_from_app(
    SequenceControlSet       *scs_ptr,
    EbSvtAv1EncConfiguration   *config_struct){

    scs_ptr->max_input_luma_width = config_struct->source_width;
    scs_ptr->max_input_luma_height = config_struct->source_height;
    scs_ptr->frame_rate = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate;
    // SB Definitions
    scs_ptr->static_config.pred_structure = ((EbSvtAv1EncConfiguration*)config_struct)->pred_structure;
    // Tpl is disabled in low delay applications
    if (scs_ptr->static_config.pred_structure == 0) {
        ((EbSvtAv1EncConfiguration*)config_struct)->enable_tpl_la = 0;
        SVT_WARN("TPL is disabled in low delay applications.\n");
    }
    scs_ptr->enable_qp_scaling_flag = 1;
    scs_ptr->max_blk_size = (uint8_t)64;
    scs_ptr->min_blk_size = (uint8_t)8;
    scs_ptr->max_intra_size = (uint8_t)32;
    scs_ptr->min_intra_size = (uint8_t)8;
    scs_ptr->max_ref_count = 1;
    scs_ptr->palette_level = DEFAULT;
    scs_ptr->intra_angle_delta = DEFAULT;
    scs_ptr->intrabc_mode = DEFAULT;
    scs_ptr->ten_bit_format = 0;
    scs_ptr->speed_control_flag = 0;

    // Padding Offsets
    scs_ptr->sb_sz = 64;
    scs_ptr->max_sb_depth = (uint8_t)EB_MAX_SB_DEPTH;
    scs_ptr->static_config.intra_period_length = ((EbSvtAv1EncConfiguration*)config_struct)->intra_period_length;
    scs_ptr->static_config.intra_refresh_type = ((EbSvtAv1EncConfiguration*)config_struct)->intra_refresh_type;
    scs_ptr->static_config.hierarchical_levels = ((EbSvtAv1EncConfiguration*)config_struct)->hierarchical_levels;
    scs_ptr->static_config.enc_mode = ((EbSvtAv1EncConfiguration*)config_struct)->enc_mode;
    scs_ptr->static_config.use_qp_file = ((EbSvtAv1EncConfiguration*)config_struct)->use_qp_file;
    scs_ptr->static_config.use_fixed_qindex_offsets = ((EbSvtAv1EncConfiguration*)config_struct)->use_fixed_qindex_offsets;
    scs_ptr->static_config.key_frame_chroma_qindex_offset = ((EbSvtAv1EncConfiguration*)config_struct)->key_frame_chroma_qindex_offset;
    scs_ptr->static_config.key_frame_qindex_offset = ((EbSvtAv1EncConfiguration*)config_struct)->key_frame_qindex_offset;
    if (scs_ptr->static_config.use_fixed_qindex_offsets == 1) {
        scs_ptr->enable_qp_scaling_flag = 0;
        scs_ptr->static_config.use_qp_file = 0;
        memcpy(scs_ptr->static_config.qindex_offsets, ((EbSvtAv1EncConfiguration*)config_struct)->qindex_offsets,
            MAX_TEMPORAL_LAYERS * sizeof(int32_t));
        memcpy(scs_ptr->static_config.chroma_qindex_offsets, ((EbSvtAv1EncConfiguration*)config_struct)->chroma_qindex_offsets,
            MAX_TEMPORAL_LAYERS * sizeof(int32_t));
    }
    scs_ptr->static_config.rc_stats_buffer = ((EbSvtAv1EncConfiguration*)config_struct)->rc_stats_buffer;
    scs_ptr->static_config.pass = ((EbSvtAv1EncConfiguration*)config_struct)->pass;
    // Deblock Filter
    scs_ptr->static_config.enable_dlf_flag = ((EbSvtAv1EncConfiguration*)config_struct)->enable_dlf_flag;

    // CDEF
    scs_ptr->static_config.cdef_level = ((EbSvtAv1EncConfiguration*)config_struct)->cdef_level;

    // Restoration filtering
    scs_ptr->static_config.enable_restoration_filtering = ((EbSvtAv1EncConfiguration*)config_struct)->enable_restoration_filtering;
    // motion field motion vector
    scs_ptr->static_config.enable_mfmv                  = ((EbSvtAv1EncConfiguration*)config_struct)->enable_mfmv;

    // Decoder Optimization Flag
    scs_ptr->static_config.fast_decode = ((EbSvtAv1EncConfiguration*)config_struct)->fast_decode;
    //If the set fast_decode value is in the allowable range, check that the value is supported for the current preset.
    // If the value is valid, but not supported in the current preset, change the value to one that is supported.
    if (scs_ptr->static_config.fast_decode &&
        (scs_ptr->static_config.fast_decode <=3)) {
        if (scs_ptr->static_config.enc_mode <= ENC_M4 ||
            (scs_ptr->static_config.enc_mode >= ENC_M8 && scs_ptr->static_config.fast_decode > 2) ||
            (scs_ptr->static_config.enc_mode >= ENC_M10 && scs_ptr->static_config.fast_decode > 1) ||
            scs_ptr->static_config.enc_mode >= ENC_M11) {
            SVT_WARN("Decoder speedup level %d is not supported in M%d.\n", scs_ptr->static_config.fast_decode, scs_ptr->static_config.enc_mode);
            SVT_WARN("Decoder speedup levels are supported as follows:\
                \n\t<= M4: not supported\n\tM5-M7: levels 1-3 supported\n\tM8-M9:levels 1-2 supported\n\tM10: level 1 supported\n\t>= M11: not supported\n");
            if (scs_ptr->static_config.enc_mode <= ENC_M4)
                scs_ptr->static_config.fast_decode = 0;
            else if (scs_ptr->static_config.enc_mode <= ENC_M9)
                scs_ptr->static_config.fast_decode = 2;
            else if (scs_ptr->static_config.enc_mode <= ENC_M10)
                scs_ptr->static_config.fast_decode = 1;
            else
                scs_ptr->static_config.fast_decode = 0;
            SVT_WARN("Switching to decoder speedup level %d.\n", scs_ptr->static_config.fast_decode);
        }
    }

    //Film Grain
    scs_ptr->static_config.film_grain_denoise_strength = ((EbSvtAv1EncConfiguration*)config_struct)->film_grain_denoise_strength;

    // MD Parameters
    scs_ptr->enable_hbd_mode_decision = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth > 8 ? DEFAULT : 0;
    // Adaptive Loop Filter
    scs_ptr->static_config.tile_rows = scs_ptr->static_config.pass == ENC_FIRST_PASS ? 0 : ((EbSvtAv1EncConfiguration*)config_struct)->tile_rows;
    scs_ptr->static_config.tile_columns = scs_ptr->static_config.pass == ENC_FIRST_PASS ? 0 : ((EbSvtAv1EncConfiguration*)config_struct)->tile_columns;
    scs_ptr->static_config.restricted_motion_vector = ((EbSvtAv1EncConfiguration*)config_struct)->restricted_motion_vector;

    // Rate Control
    scs_ptr->static_config.scene_change_detection = ((EbSvtAv1EncConfiguration*)config_struct)->scene_change_detection;
    scs_ptr->static_config.rate_control_mode = ((EbSvtAv1EncConfiguration*)config_struct)->rate_control_mode;

    if (scs_ptr->static_config.rate_control_mode == 2) {
        SVT_WARN("CBR Rate control is currently not supported, switching to VBR\n");
        scs_ptr->static_config.rate_control_mode = 1;
    }

    if (scs_ptr->static_config.rate_control_mode == 2 && scs_ptr->static_config.pass != ENC_FIRST_PASS && !(scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS) &&
        scs_ptr->static_config.pred_structure != 0) {
        scs_ptr->static_config.pred_structure = 0;
        SVT_WARN("Forced 1pass CBR to be always low delay mode.\n");
        if(((EbSvtAv1EncConfiguration*)config_struct)->enable_tpl_la) {
            ((EbSvtAv1EncConfiguration*)config_struct)->enable_tpl_la = 0;
            SVT_WARN("TPL is disabled in low delay applications.\n");
        }
    }
    // for 1pass CBR not real time mode
    //if (scs_ptr->static_config.rate_control_mode == 2 && !use_output_stat(scs_ptr) && !use_input_stat(scs_ptr))
    //    scs_ptr->static_config.hierarchical_levels = 0;

    scs_ptr->max_temporal_layers = scs_ptr->static_config.hierarchical_levels;
    scs_ptr->static_config.look_ahead_distance = ((EbSvtAv1EncConfiguration*)config_struct)->look_ahead_distance;
    scs_ptr->static_config.frame_rate = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate;
    scs_ptr->static_config.frame_rate_denominator = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate_denominator;
    scs_ptr->static_config.frame_rate_numerator = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate_numerator;

    scs_ptr->static_config.target_bit_rate = ((EbSvtAv1EncConfiguration*)config_struct)->target_bit_rate;
    scs_ptr->static_config.max_bit_rate = ((EbSvtAv1EncConfiguration*)config_struct)->max_bit_rate;
    if (((EbSvtAv1EncConfiguration*)config_struct)->enable_tpl_la == 0 && scs_ptr->static_config.max_bit_rate) {
        scs_ptr->static_config.max_bit_rate = 0;
        SVT_WARN("Maximum bit rate only supported with tpl on. max bit rate 0 is used instead.\n");
    }
    scs_ptr->static_config.vbv_bufsize = ((EbSvtAv1EncConfiguration*)config_struct)->vbv_bufsize;

    scs_ptr->static_config.max_qp_allowed = (scs_ptr->static_config.rate_control_mode) ?
        ((EbSvtAv1EncConfiguration*)config_struct)->max_qp_allowed :
        63;

    scs_ptr->static_config.min_qp_allowed = (scs_ptr->static_config.rate_control_mode) ?
        (((EbSvtAv1EncConfiguration*)config_struct)->min_qp_allowed > 0) ?
        ((EbSvtAv1EncConfiguration*)config_struct)->min_qp_allowed : 1 :
        1; // lossless coding not supported

    scs_ptr->static_config.vbr_bias_pct        = ((EbSvtAv1EncConfiguration*)config_struct)->vbr_bias_pct;
    scs_ptr->static_config.vbr_min_section_pct = ((EbSvtAv1EncConfiguration*)config_struct)->vbr_min_section_pct;
    scs_ptr->static_config.vbr_max_section_pct = ((EbSvtAv1EncConfiguration*)config_struct)->vbr_max_section_pct;
    scs_ptr->static_config.under_shoot_pct     = ((EbSvtAv1EncConfiguration*)config_struct)->under_shoot_pct;
    scs_ptr->static_config.over_shoot_pct      = ((EbSvtAv1EncConfiguration*)config_struct)->over_shoot_pct;
    scs_ptr->static_config.mbr_over_shoot_pct  = ((EbSvtAv1EncConfiguration*)config_struct)->mbr_over_shoot_pct;
    scs_ptr->static_config.maximum_buffer_size_ms   = ((EbSvtAv1EncConfiguration*)config_struct)->maximum_buffer_size_ms;
    scs_ptr->static_config.starting_buffer_level_ms = ((EbSvtAv1EncConfiguration*)config_struct)->starting_buffer_level_ms;
    scs_ptr->static_config.optimal_buffer_level_ms  = ((EbSvtAv1EncConfiguration*)config_struct)->optimal_buffer_level_ms;
    scs_ptr->static_config.recode_loop         = ((EbSvtAv1EncConfiguration*)config_struct)->recode_loop;
    if (scs_ptr->static_config.rate_control_mode == 1 && scs_ptr->static_config.pass != ENC_MIDDLE_PASS && scs_ptr->static_config.pass != ENC_FIRST_PASS && !(scs_ptr->static_config.pass == ENC_MIDDLE_PASS || scs_ptr->static_config.pass == ENC_LAST_PASS))
        scs_ptr->lap_enabled = 1;
    else
        scs_ptr->lap_enabled = 0;
    //Segmentation
    //TODO: check RC mode and set only when RC is enabled in the final version.
    scs_ptr->static_config.enable_adaptive_quantization = config_struct->enable_adaptive_quantization;

    // Misc
    scs_ptr->static_config.encoder_bit_depth = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth;
    scs_ptr->static_config.encoder_color_format = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_color_format;

    scs_ptr->chroma_format_idc = (uint32_t)(scs_ptr->static_config.encoder_color_format);
    scs_ptr->encoder_bit_depth = (uint32_t)(scs_ptr->static_config.encoder_bit_depth);
    //16bit pipeline
    scs_ptr->is_16bit_pipeline = ((((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth) > EB_8BIT) ? EB_TRUE: EB_FALSE;
    scs_ptr->subsampling_x = (scs_ptr->chroma_format_idc == EB_YUV444 ? 1 : 2) - 1;
    scs_ptr->subsampling_y = (scs_ptr->chroma_format_idc >= EB_YUV422 ? 1 : 2) - 1;
    scs_ptr->static_config.compressed_ten_bit_format = ((EbSvtAv1EncConfiguration*)config_struct)->compressed_ten_bit_format;

    // Thresholds
    scs_ptr->static_config.high_dynamic_range_input = ((EbSvtAv1EncConfiguration*)config_struct)->high_dynamic_range_input;
    scs_ptr->static_config.screen_content_mode = ((EbSvtAv1EncConfiguration*)config_struct)->screen_content_mode;
    // SC detection is OFF for first pass in M8
    uint8_t disable_sc_detection = (scs_ptr->static_config.pass == ENC_FIRST_PASS) ? 1 : 0;
    if (disable_sc_detection)
        scs_ptr->static_config.screen_content_mode = 0;

    // Annex A parameters
    scs_ptr->static_config.profile = ((EbSvtAv1EncConfiguration*)config_struct)->profile;
    scs_ptr->static_config.tier = ((EbSvtAv1EncConfiguration*)config_struct)->tier;
    scs_ptr->static_config.level = ((EbSvtAv1EncConfiguration*)config_struct)->level;
    scs_ptr->static_config.stat_report = ((EbSvtAv1EncConfiguration*)config_struct)->stat_report;

    // Buffers - Hardcoded(Cleanup)
    scs_ptr->static_config.use_cpu_flags = ((EbSvtAv1EncConfiguration*)config_struct)->use_cpu_flags;

    scs_ptr->static_config.channel_id = ((EbSvtAv1EncConfiguration*)config_struct)->channel_id;
    scs_ptr->static_config.active_channel_count = ((EbSvtAv1EncConfiguration*)config_struct)->active_channel_count;
    scs_ptr->static_config.logical_processors = ((EbSvtAv1EncConfiguration*)config_struct)->logical_processors;
    scs_ptr->static_config.pin_threads = ((EbSvtAv1EncConfiguration*)config_struct)->pin_threads;
    scs_ptr->static_config.target_socket = ((EbSvtAv1EncConfiguration*)config_struct)->target_socket;
    if ((scs_ptr->static_config.pin_threads == 0) && (scs_ptr->static_config.target_socket != -1)){
        SVT_WARN("threads pinning 0 and ss %d is not a valid combination: unpin will be set to 0\n", scs_ptr->static_config.target_socket);
        scs_ptr->static_config.pin_threads = 1;
    }
    scs_ptr->static_config.qp = ((EbSvtAv1EncConfiguration*)config_struct)->qp;
    scs_ptr->static_config.recon_enabled = ((EbSvtAv1EncConfiguration*)config_struct)->recon_enabled;
    scs_ptr->static_config.enable_tpl_la = ((EbSvtAv1EncConfiguration*)config_struct)->enable_tpl_la;
    if (scs_ptr->static_config.enable_tpl_la != 1){
        SVT_WARN("TPL off mode is not supported in this release, enable_tpl_la is set to 1\n");
        scs_ptr->static_config.enable_tpl_la = 1;
    }
    // Extract frame rate from Numerator and Denominator if not 0
    if (scs_ptr->static_config.frame_rate_numerator != 0 && scs_ptr->static_config.frame_rate_denominator != 0)
        scs_ptr->frame_rate = scs_ptr->static_config.frame_rate = (((scs_ptr->static_config.frame_rate_numerator << 8) / (scs_ptr->static_config.frame_rate_denominator)) << 8);
    // Get Default Intra Period if not specified
    if (scs_ptr->static_config.intra_period_length == -2)
        scs_ptr->static_config.intra_period_length = compute_default_intra_period(scs_ptr);
    if (scs_ptr->static_config.look_ahead_distance == (uint32_t)~0)
        scs_ptr->static_config.look_ahead_distance = compute_default_look_ahead(&scs_ptr->static_config);
    scs_ptr->static_config.enable_tf = config_struct->enable_tf;
    scs_ptr->static_config.enable_overlays = config_struct->enable_overlays;
    scs_ptr->static_config.tune = config_struct->tune;
    scs_ptr->static_config.superres_mode = config_struct->superres_mode;
    scs_ptr->static_config.superres_denom = config_struct->superres_denom;
    scs_ptr->static_config.superres_kf_denom = config_struct->superres_kf_denom;
    scs_ptr->static_config.superres_qthres = config_struct->superres_qthres;
    scs_ptr->static_config.superres_kf_qthres = config_struct->superres_kf_qthres;
    if (scs_ptr->static_config.superres_mode == SUPERRES_AUTO)
    {
        // TODO: set search mode based on preset
        //scs_ptr->static_config.superres_auto_search_type = SUPERRES_AUTO_SOLO;
        scs_ptr->static_config.superres_auto_search_type = SUPERRES_AUTO_DUAL;
        //scs_ptr->static_config.superres_auto_search_type = SUPERRES_AUTO_ALL;
    }

    // Prediction Structure
    scs_ptr->static_config.enable_manual_pred_struct    = config_struct->enable_manual_pred_struct;
    if(scs_ptr->static_config.enable_manual_pred_struct){
        scs_ptr->static_config.manual_pred_struct_entry_num = config_struct->manual_pred_struct_entry_num;
        EB_MEMCPY(&scs_ptr->static_config.pred_struct[0], &config_struct->pred_struct[0],config_struct->manual_pred_struct_entry_num*sizeof(PredictionStructureConfigEntry));
        switch (scs_ptr->static_config.manual_pred_struct_entry_num) {
            case 1:
                scs_ptr->static_config.hierarchical_levels =  0;
                break;
            case 2:
                scs_ptr->static_config.hierarchical_levels =  1;
                break;
            case 4:
                scs_ptr->static_config.hierarchical_levels =  2;
                break;
            case 8:
                scs_ptr->static_config.hierarchical_levels =  3;
                break;
            case 16:
                scs_ptr->static_config.hierarchical_levels =  4;
                break;
            case 32:
                scs_ptr->static_config.hierarchical_levels =  5;
                break;
            default:
                scs_ptr->static_config.hierarchical_levels =  0;
                break;
        }
    }

    // Color description
    scs_ptr->static_config.color_description_present_flag = config_struct->color_description_present_flag;
    scs_ptr->static_config.color_primaries = config_struct->color_primaries;
    scs_ptr->static_config.transfer_characteristics = config_struct->transfer_characteristics;
    scs_ptr->static_config.matrix_coefficients = config_struct->matrix_coefficients;
    scs_ptr->static_config.color_range = config_struct->color_range;
    scs_ptr->static_config.mastering_display = config_struct->mastering_display;
    scs_ptr->static_config.content_light_level = config_struct->content_light_level;

    return;
}

/**********************************

* Set Parameter
**********************************/
EB_API EbErrorType svt_av1_enc_set_parameter(
    EbComponentType              *svt_enc_component,
    EbSvtAv1EncConfiguration     *config_struct)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;

    EbEncHandle        *enc_handle  = (EbEncHandle*)svt_enc_component->p_component_private;
    uint32_t              instance_index = 0;

    // Acquire Config Mutex
    svt_block_on_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);

    copy_api_from_app(
        enc_handle->scs_instance_array[instance_index]->scs_ptr,
        (EbSvtAv1EncConfiguration*)config_struct);

    EbErrorType return_error = svt_av1_verify_settings(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    if (return_error == EB_ErrorBadParameter)
        return EB_ErrorBadParameter;

    set_param_based_on_input(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);
    MiniGopSizeCtrls *mgs_ctls = &enc_handle->scs_instance_array[instance_index]->scs_ptr->mgs_ctls;
    uint8_t mg_level = (enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.pass == ENC_MIDDLE_PASS || enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.pass == ENC_LAST_PASS) ? 1 : 0;
    enc_handle->scs_instance_array[instance_index]->scs_ptr->max_heirachical_level = enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.hierarchical_levels;
    enc_handle->scs_instance_array[instance_index]->scs_ptr->enable_adaptive_mini_gop = 0;
    set_mini_gop_size_controls(mgs_ctls, mg_level, enc_handle->scs_instance_array[instance_index]->scs_ptr->input_resolution);
    if (mgs_ctls->adptive_enable)
        set_max_mini_gop_size(
            enc_handle->scs_instance_array[instance_index]->scs_ptr, mgs_ctls);
    // Initialize the Prediction Structure Group
    EB_NO_THROW_NEW(
        enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        prediction_structure_group_ctor,
        enc_handle->scs_instance_array[instance_index]->scs_ptr);
    if (!enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr) {
        svt_release_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);
        return EB_ErrorInsufficientResources;
    }
    // Set the Prediction Structure
    enc_handle->scs_instance_array[instance_index]->scs_ptr->pred_struct_ptr = get_prediction_structure(
        enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.pred_structure,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->max_ref_count,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->max_temporal_layers);

    return_error = load_default_buffer_configuration_settings(
        enc_handle,
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    svt_av1_print_lib_params(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    // Release Config Mutex
    svt_release_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);

    return return_error;
}
EB_API EbErrorType svt_av1_enc_stream_header(
    EbComponentType           *svt_enc_component,
    EbBufferHeaderType        **output_stream_ptr)
{
    EbErrorType              return_error = EB_ErrorNone;

    if(!svt_enc_component)
        return EB_ErrorBadParameter;

    EbEncHandle             *enc_handle  = (EbEncHandle*)svt_enc_component->p_component_private;
    SequenceControlSet      *scs_ptr = enc_handle->scs_instance_array[0]->scs_ptr;
    Bitstream                bitstream;
    OutputBitstreamUnit      output_bitstream;
    EbBufferHeaderType      *output_stream_buffer;
    uint32_t output_buffer_size = get_out_buffer_size(scs_ptr->max_input_luma_width, scs_ptr->max_input_luma_height);
    memset(&bitstream, 0, sizeof(Bitstream));
    memset(&output_bitstream, 0, sizeof(OutputBitstreamUnit));
    bitstream.output_bitstream_ptr = &output_bitstream;
    output_stream_buffer = (EbBufferHeaderType *)malloc(sizeof(EbBufferHeaderType));
    if (!output_stream_buffer) {
        return EB_ErrorInsufficientResources;
    }
    output_stream_buffer->p_buffer = (uint8_t *)malloc(sizeof(uint8_t) * output_buffer_size);
    if (!output_stream_buffer->p_buffer) {
        free(output_stream_buffer);
        return EB_ErrorInsufficientResources;
    }

    output_stream_buffer->size = sizeof(EbBufferHeaderType);
    output_stream_buffer->n_alloc_len = output_buffer_size;
    output_stream_buffer->p_app_private = NULL;
    output_stream_buffer->pic_type = EB_AV1_INVALID_PICTURE;
    output_stream_buffer->n_filled_len = 0;

    ((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_begin_av1 = output_stream_buffer->p_buffer;

    output_bitstream_reset(bitstream.output_bitstream_ptr);

    // Code the SPS
    encode_sps_av1(&bitstream, scs_ptr);

    output_stream_buffer->n_filled_len = (uint32_t)(((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_av1 - ((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_begin_av1);

    *output_stream_ptr = output_stream_buffer;

    return return_error;
}
EB_API EbErrorType svt_av1_enc_stream_header_release(
    EbBufferHeaderType        *stream_header_ptr)
{
    EbErrorType           return_error = EB_ErrorNone;

    if (!stream_header_ptr || !(stream_header_ptr->p_buffer)) {
        return EB_ErrorBadParameter;
    }

    free(stream_header_ptr->p_buffer);
    free(stream_header_ptr);

    return return_error;
}

/***********************************************
**** Copy the input buffer from the
**** sample application to the library buffers
************************************************/
/*
 Down sample and Copy the input buffer
from the sample application to the library buffers
*/
/********************************************
 * downsample_2d_c_16_zero2bit
 *      downsample the input by getting the average and zero out the two LSB bit
 ********************************************/
void downsample_2d_c_16_zero2bit(uint16_t *input_samples, // input parameter, input samples Ptr
    uint32_t input_stride, // input parameter, input stride
    uint32_t input_area_width, // input parameter, input area width
    uint32_t input_area_height, // input parameter, input area height
    uint8_t *decim_8b_samples, // output parameter, decimated samples Ptr
    uint32_t decim_stride, // input parameter, output stride
    uint32_t decim_step) // input parameter, decimation amount in pixels
{
    uint32_t       horizontal_index;
    uint32_t       vertical_index;
    uint32_t       input_stripe_stride = input_stride * decim_step;
    uint32_t       decim_horizontal_index;
    const uint32_t half_decim_step = decim_step >> 1;
    for (input_samples += half_decim_step * input_stride, vertical_index = half_decim_step;
        vertical_index < input_area_height;
        vertical_index += decim_step) {
        uint16_t *prev_input_line = input_samples - input_stride;
        for (horizontal_index = half_decim_step, decim_horizontal_index = 0;
            horizontal_index < input_area_width;
            horizontal_index += decim_step, decim_horizontal_index++) {
                uint32_t sum = (uint32_t)prev_input_line[horizontal_index - 1] +
                    (uint32_t)prev_input_line[horizontal_index] +
                    (uint32_t)input_samples[horizontal_index - 1] +
                    (uint32_t)input_samples[horizontal_index];
                decim_8b_samples[decim_horizontal_index] = (uint8_t)(((sum + 2) >> 2) >> 2);

        }
        input_samples += input_stripe_stride;
        decim_8b_samples += decim_stride;
    }

    return;
}
/********************************************
 * downsample_2d_c_16_zero2bit_skipall
 *      downsample the input by skipping three pixels and zero out the two LSB bit
 ********************************************/
void downsample_2d_c_16_zero2bit_skipall(uint16_t *input_samples, // input parameter, input samples Ptr
    uint32_t input_stride, // input parameter, input stride
    uint32_t input_area_width, // input parameter, input area width
    uint32_t input_area_height, // input parameter, input area height
    uint8_t *decim_8b_samples, // output parameter, decimated samples Ptr
    uint32_t decim_stride, // input parameter, output stride
    uint32_t decim_step) // input parameter, decimation amount in pixels
{
    uint32_t       horizontal_index;
    uint32_t       vertical_index;
    uint32_t       input_stripe_stride = input_stride * decim_step;
    uint32_t       decim_horizontal_index;
    const uint32_t half_decim_step = decim_step >> 1;

    for (input_samples += half_decim_step * input_stride, vertical_index = half_decim_step;
        vertical_index < input_area_height;
        vertical_index += decim_step) {
        uint16_t *prev_input_line = input_samples - input_stride;
        for (horizontal_index = half_decim_step, decim_horizontal_index = 0;
            horizontal_index < input_area_width;
            horizontal_index += decim_step, decim_horizontal_index++) {
            decim_8b_samples[decim_horizontal_index] = (uint8_t)((prev_input_line[horizontal_index - 1]) >> 2);
        }
        input_samples += input_stripe_stride;
        decim_8b_samples += decim_stride;
    }

    return;
}
/********************************************
 * downsample_2d_c_skipall
 *      downsample the input by skipping three pixels
 ********************************************/
void downsample_2d_c_skipall(uint8_t *input_samples, // input parameter, input samples Ptr
    uint32_t input_stride, // input parameter, input stride
    uint32_t input_area_width, // input parameter, input area width
    uint32_t input_area_height, // input parameter, input area height
    uint8_t *decim_samples, // output parameter, decimated samples Ptr
    uint32_t decim_stride, // input parameter, output stride
    uint32_t decim_step) // input parameter, decimation amount in pixels
{
    uint32_t       horizontal_index;
    uint32_t       vertical_index;
    uint32_t       input_stripe_stride = input_stride * decim_step;
    uint32_t       decim_horizontal_index;
    const uint32_t half_decim_step = decim_step >> 1;

    for (input_samples += half_decim_step * input_stride, vertical_index = half_decim_step;
        vertical_index < input_area_height;
        vertical_index += decim_step) {
        uint8_t *prev_input_line = input_samples - input_stride;
        for (horizontal_index = half_decim_step, decim_horizontal_index = 0;
            horizontal_index < input_area_width;
            horizontal_index += decim_step, decim_horizontal_index++) {
            decim_samples[decim_horizontal_index] = (uint32_t)prev_input_line[horizontal_index - 1];
        }
        input_samples += input_stripe_stride;
        decim_samples += decim_stride;
    }

    return;
}
/***********************************************
 Down sample and Copy the input buffer
from the sample application to the library buffers
************************************************/
static EbErrorType downsample_copy_frame_buffer(
    SequenceControlSet            *scs_ptr,
    uint8_t                       *destination,
    uint8_t                       *destination_y8b,
    uint8_t                       *source,
    int                            pass)
{
    EbSvtAv1EncConfiguration          *config = &scs_ptr->static_config;
    EbErrorType                      return_error = EB_ErrorNone;

    EbPictureBufferDesc           *input_picture_ptr = (EbPictureBufferDesc*)destination;
    EbPictureBufferDesc           *y8b_input_picture_ptr = (EbPictureBufferDesc*)destination_y8b;
    EbSvtIOFormat                   *input_ptr = (EbSvtIOFormat*)source;
    EbBool                           is_16bit_input = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    uint8_t                         *src, *dst;

    // Need to include for Interlacing on the fly with pictureScanType = 1
    if (!is_16bit_input) {
        uint32_t     luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding) << is_16bit_input;
        uint32_t     chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1)) << is_16bit_input;
        uint16_t     luma_stride = input_picture_ptr->stride_y << is_16bit_input;
        uint16_t     chroma_stride = input_picture_ptr->stride_cb << is_16bit_input;
        uint16_t     luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);
        uint16_t     luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right);

        uint16_t     source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t     source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t     source_cb_stride = (uint16_t)(input_ptr->cb_stride);
        uint16_t source_chroma_height =
            (luma_height >> (input_picture_ptr->color_format == EB_YUV420));
        uint16_t source_chroma_width =
            (luma_width >> (input_picture_ptr->color_format == EB_YUV420));

        src = input_ptr->luma;
        dst = y8b_input_picture_ptr->buffer_y + luma_buffer_offset;
        downsample_2d_c_skipall(
            src,
            source_luma_stride,
            luma_width << 1,
            luma_height<<1,
            dst,
            luma_stride,
            2);

#define ENCODE_FIRST_PASS 1
        if (pass != ENCODE_FIRST_PASS) {
            src = input_ptr->cb;
            dst = input_picture_ptr->buffer_cb + chroma_buffer_offset;
            downsample_2d_c_skipall(
                src,
                source_cb_stride,
                source_chroma_width << 1,
                source_chroma_height << 1,
                dst,
                chroma_stride,
                2);

            src = input_ptr->cr;
            dst = input_picture_ptr->buffer_cr + chroma_buffer_offset;
            downsample_2d_c_skipall(
                src,
                source_cr_stride,
                source_chroma_width << 1,
                source_chroma_height << 1,
                dst,
                chroma_stride,
                2);
        }
    }

    else if (config->compressed_ten_bit_format == 1)
    {
        {
            SVT_WARN("Compressed_ten_bit_format not supported in downsample_copy_frame_buffer");//anaghdin
            uint32_t  luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
            uint32_t  chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
            uint16_t  luma_stride = input_picture_ptr->stride_y;
            uint16_t  chroma_stride = input_picture_ptr->stride_cb;
            uint16_t  luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

            uint16_t  source_luma_stride = (uint16_t)(input_ptr->y_stride);
            uint16_t  source_cr_stride = (uint16_t)(input_ptr->cr_stride);
            uint16_t  source_cb_stride = (uint16_t)(input_ptr->cb_stride);
            uint16_t source_chroma_height =
                (luma_height >> (input_picture_ptr->color_format == EB_YUV420));

            src = input_ptr->luma;
            dst = y8b_input_picture_ptr->buffer_y + luma_buffer_offset;
            for (unsigned i = 0; i < luma_height; i++) {
                svt_memcpy(dst, src, source_luma_stride);
                src += source_luma_stride;
                dst += luma_stride;
            }
            if (pass != ENCODE_FIRST_PASS) {
                src = input_ptr->cb;
                dst = input_picture_ptr->buffer_cb + chroma_buffer_offset;
                for (unsigned i = 0; i < source_chroma_height; i++) {
                    svt_memcpy(dst, src, source_cb_stride);
                    src += source_cb_stride;
                    dst += chroma_stride;
                }

                src = input_ptr->cr;
                dst = input_picture_ptr->buffer_cr + chroma_buffer_offset;
                for (unsigned i = 0; i < source_chroma_height; i++) {
                    svt_memcpy(dst, src, source_cr_stride);
                    src += source_cr_stride;
                    dst += chroma_stride;
                }
            }
            //efficient copy - final
            //compressed 2Bit in 1D format
            {
                uint32_t comp_stride_y = luma_stride / 4;
                uint32_t comp_luma_buffer_offset = comp_stride_y * input_picture_ptr->origin_y + input_picture_ptr->origin_x / 4;

                uint32_t comp_stride_uv = chroma_stride / 4;
                uint32_t comp_chroma_buffer_offset = comp_stride_uv * (input_picture_ptr->origin_y / 2) + input_picture_ptr->origin_x / 2 / 4;

                src = input_ptr->luma_ext;
                dst = input_picture_ptr->buffer_bit_inc_y + comp_luma_buffer_offset;
                for (unsigned i = 0; i < luma_height; i++) {
                    svt_memcpy(dst, src, source_luma_stride >> 2);
                    src += source_luma_stride >> 2;
                    dst += comp_stride_y;
                }
                if (pass != ENCODE_FIRST_PASS) {
                    src = input_ptr->cb_ext;
                    dst = input_picture_ptr->buffer_bit_inc_cb + comp_chroma_buffer_offset;
                    for (unsigned i = 0; i < source_chroma_height; i++) {
                        svt_memcpy(dst, src, source_cb_stride >> 2);
                        src += source_cb_stride >> 2;
                        dst += comp_stride_uv;
                    }

                    src = input_ptr->cr_ext;
                    dst = input_picture_ptr->buffer_bit_inc_cr + comp_chroma_buffer_offset;
                    for (unsigned i = 0; i < source_chroma_height; i++) {
                        svt_memcpy(dst, src, source_cr_stride >> 2);
                        src += source_cr_stride >> 2;
                        dst += comp_stride_uv;
                    }
                }
            }
        }
    }
    else { // 10bit packed

    uint32_t luma_offset = 0;
        uint32_t luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
        uint32_t chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
        uint16_t luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right);
        uint16_t luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        uint16_t source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t source_cb_stride = (uint16_t)(input_ptr->cb_stride);

        downsample_2d_c_16_zero2bit_skipall(
            (uint16_t*)(uint16_t*)(input_ptr->luma + luma_offset),
            source_luma_stride,
            luma_width << 1,
            luma_height << 1,
            (y8b_input_picture_ptr->buffer_y + luma_buffer_offset),
            y8b_input_picture_ptr->stride_y,
            2);

        memset(input_picture_ptr->buffer_bit_inc_y, 0, input_picture_ptr->luma_size/4);

        if (pass != ENCODE_FIRST_PASS) {
            uint32_t chroma_offset = 0;

            downsample_2d_c_16_zero2bit_skipall(
                (uint16_t*)(input_ptr->cb + chroma_offset),
                source_cb_stride,
                luma_width,
                luma_height,
                input_picture_ptr->buffer_cb + chroma_buffer_offset,
                y8b_input_picture_ptr->stride_cb,
                2);

            memset(input_picture_ptr->buffer_bit_inc_cb, 0, input_picture_ptr->chroma_size/4);

            downsample_2d_c_16_zero2bit_skipall(
                (uint16_t*)(input_ptr->cr + chroma_offset),
                source_cr_stride,
                luma_width,
                luma_height,
                input_picture_ptr->buffer_cr + chroma_buffer_offset,
                y8b_input_picture_ptr->stride_cr,
                2);

            memset(input_picture_ptr->buffer_bit_inc_cr, 0, input_picture_ptr->chroma_size/4);
        }
    }
    return return_error;
}
/*
 Copy the input buffer
from the sample application to the library buffers
*/

static EbErrorType copy_frame_buffer(
    SequenceControlSet            *scs_ptr,
    uint8_t                       *destination,
    uint8_t                       *destination_y8b,
    uint8_t                       *source,
    int                            pass)
{
    EbSvtAv1EncConfiguration          *config = &scs_ptr->static_config;
    EbErrorType                      return_error = EB_ErrorNone;

    EbPictureBufferDesc           *input_picture_ptr = (EbPictureBufferDesc*)destination;
    EbPictureBufferDesc           *y8b_input_picture_ptr = (EbPictureBufferDesc*)destination_y8b;
    EbSvtIOFormat                   *input_ptr = (EbSvtIOFormat*)source;
    EbBool                           is_16bit_input = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    uint8_t                         *src, *dst;

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is_16bit_input) {
        uint32_t     luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding) << is_16bit_input;
        uint32_t     chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1)) << is_16bit_input;
        uint16_t     luma_stride = input_picture_ptr->stride_y << is_16bit_input;
        uint16_t     chroma_stride = input_picture_ptr->stride_cb << is_16bit_input;
        uint16_t     luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        uint16_t     source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t     source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t     source_cb_stride = (uint16_t)(input_ptr->cb_stride);
        uint16_t source_chroma_height =
            (luma_height >> (input_picture_ptr->color_format == EB_YUV420));

        src = input_ptr->luma;
        dst = y8b_input_picture_ptr->buffer_y + luma_buffer_offset;
        for (unsigned i = 0; i < luma_height; i++) {
            svt_memcpy(dst, src, source_luma_stride);
            src += source_luma_stride;
            dst += luma_stride;
        }

#define ENCODE_FIRST_PASS 1
        if (pass != ENCODE_FIRST_PASS) {
            src = input_ptr->cb;
            dst = input_picture_ptr->buffer_cb + chroma_buffer_offset;
            for (unsigned i = 0; i < source_chroma_height; i++) {
                svt_memcpy(dst, src, source_cb_stride);
                src += source_cb_stride;
                dst += chroma_stride;
            }

            src = input_ptr->cr;
            dst = input_picture_ptr->buffer_cr + chroma_buffer_offset;
            for (unsigned i = 0; i < source_chroma_height; i++) {
                svt_memcpy(dst, src, source_cr_stride);
                src += source_cr_stride;
                dst += chroma_stride;
            }
        }
    }

    else if (config->compressed_ten_bit_format == 1)
    {
        {
            uint32_t  luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
            uint32_t  chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
            uint16_t  luma_stride = input_picture_ptr->stride_y;
            uint16_t  chroma_stride = input_picture_ptr->stride_cb;
            uint16_t  luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

            uint16_t  source_luma_stride = (uint16_t)(input_ptr->y_stride);
            uint16_t  source_cr_stride = (uint16_t)(input_ptr->cr_stride);
            uint16_t  source_cb_stride = (uint16_t)(input_ptr->cb_stride);
            uint16_t source_chroma_height =
                (luma_height >> (input_picture_ptr->color_format == EB_YUV420));

            src = input_ptr->luma;
            dst = y8b_input_picture_ptr->buffer_y + luma_buffer_offset;
            for (unsigned i = 0; i < luma_height; i++) {
                svt_memcpy(dst, src, source_luma_stride);
                src += source_luma_stride;
                dst += luma_stride;
            }
            if (pass != ENCODE_FIRST_PASS) {
                src = input_ptr->cb;
                dst = input_picture_ptr->buffer_cb + chroma_buffer_offset;
                for (unsigned i = 0; i < source_chroma_height; i++) {
                    svt_memcpy(dst, src, source_cb_stride);
                    src += source_cb_stride;
                    dst += chroma_stride;
                }

                src = input_ptr->cr;
                dst = input_picture_ptr->buffer_cr + chroma_buffer_offset;
                for (unsigned i = 0; i < source_chroma_height; i++) {
                    svt_memcpy(dst, src, source_cr_stride);
                    src += source_cr_stride;
                    dst += chroma_stride;
                }
            }
            //efficient copy - final
            //compressed 2Bit in 1D format
            {
                uint32_t comp_stride_y = luma_stride / 4;
                uint32_t comp_luma_buffer_offset = comp_stride_y * input_picture_ptr->origin_y + input_picture_ptr->origin_x/4;

                uint32_t comp_stride_uv = chroma_stride / 4;
                uint32_t comp_chroma_buffer_offset = comp_stride_uv * (input_picture_ptr->origin_y/2) + input_picture_ptr->origin_x /2 / 4;

                src = input_ptr->luma_ext;
                dst = input_picture_ptr->buffer_bit_inc_y + comp_luma_buffer_offset;
                for (unsigned i = 0; i < luma_height; i++) {
                    svt_memcpy(dst, src, source_luma_stride >> 2);
                    src += source_luma_stride >> 2;
                    dst += comp_stride_y;
                }
                if (pass != ENCODE_FIRST_PASS) {
                    src = input_ptr->cb_ext;
                    dst = input_picture_ptr->buffer_bit_inc_cb + comp_chroma_buffer_offset;
                    for (unsigned i = 0; i < source_chroma_height; i++) {
                        svt_memcpy(dst, src, source_cb_stride >> 2);
                        src += source_cb_stride >> 2;
                        dst += comp_stride_uv;
                    }

                    src = input_ptr->cr_ext;
                    dst = input_picture_ptr->buffer_bit_inc_cr + comp_chroma_buffer_offset;
                    for (unsigned i = 0; i < source_chroma_height; i++) {
                        svt_memcpy(dst, src, source_cr_stride >> 2);
                        src += source_cr_stride >> 2;
                        dst += comp_stride_uv;
                    }
                }
            }
        }
    }
    else { // 10bit packed
        uint32_t luma_offset = 0;
        uint32_t luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
        uint32_t chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
        uint16_t luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right);
        uint16_t luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        uint16_t source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t source_cb_stride = (uint16_t)(input_ptr->cb_stride);

        uint32_t comp_stride_y = input_picture_ptr->stride_y / 4;
        uint32_t comp_luma_buffer_offset = comp_stride_y * input_picture_ptr->origin_y + input_picture_ptr->origin_x/4;

        uint32_t comp_stride_uv = input_picture_ptr->stride_cb / 4;
        uint32_t comp_chroma_buffer_offset = comp_stride_uv * (input_picture_ptr->origin_y/2) + input_picture_ptr->origin_x /2 / 4;

        svt_unpack_and_2bcompress(
            (uint16_t*)(input_ptr->luma + luma_offset),
            source_luma_stride,
            y8b_input_picture_ptr->buffer_y + luma_buffer_offset,
            y8b_input_picture_ptr->stride_y,
            input_picture_ptr->buffer_bit_inc_y + comp_luma_buffer_offset,
            comp_stride_y,
            luma_width,
            luma_height);
        if (pass != ENCODE_FIRST_PASS) {
            uint32_t chroma_offset = 0;
        svt_unpack_and_2bcompress(
            (uint16_t*)(input_ptr->cb + chroma_offset),
            source_cb_stride,
            input_picture_ptr->buffer_cb + chroma_buffer_offset,
            input_picture_ptr->stride_cb,
            input_picture_ptr->buffer_bit_inc_cb + comp_chroma_buffer_offset,
            comp_stride_uv,
            luma_width / 2,
            luma_height / 2);

        svt_unpack_and_2bcompress(
            (uint16_t*)(input_ptr->cr + chroma_offset),
            source_cr_stride,
            input_picture_ptr->buffer_cr + chroma_buffer_offset,
            input_picture_ptr->stride_cr,
            input_picture_ptr->buffer_bit_inc_cr + comp_chroma_buffer_offset,
            comp_stride_uv,
            luma_width / 2,
            luma_height / 2);
        }
    }
    return return_error;
}
/*
 Copy the input buffer header content
from the sample application to the library buffers
*/
static void copy_input_buffer(
    SequenceControlSet*    sequenceControlSet,
    EbBufferHeaderType*     dst,
    EbBufferHeaderType*     dst_y8b,
    EbBufferHeaderType*     src,
    int                     pass
)
{
    // Copy the higher level structure
    dst->n_alloc_len = src->n_alloc_len;
    dst->n_filled_len = src->n_filled_len;
    dst->flags = src->flags;
    dst->pts = src->pts;
    dst->n_tick_count = src->n_tick_count;
    dst->size = src->size;
    dst->qp = src->qp;
    dst->pic_type = src->pic_type;

    int copy_frame = 1;
    if (sequenceControlSet->ipp_pass_ctrls.skip_frame_first_pass == 1)
        copy_frame = (((src->pts % 8) == 0) || ((src->pts % 8) == 6) || ((src->pts % 8) == 7));
    else if (sequenceControlSet->ipp_pass_ctrls.skip_frame_first_pass == 2)
        copy_frame = ((src->pts < 7) || ((src->pts % 8) == 0) || ((src->pts % 8) == 6) || ((src->pts % 8) == 7));
    if (sequenceControlSet->mid_pass_ctrls.ds || sequenceControlSet->ipp_pass_ctrls.ds) {
        // Copy the picture buffer
        if (src->p_buffer != NULL)
            downsample_copy_frame_buffer(sequenceControlSet, dst->p_buffer, dst_y8b->p_buffer, src->p_buffer, pass);
    }
    else

    // Bypass copy for the unecessary picture in IPPP pass
    if ((pass != ENCODE_FIRST_PASS) || ((pass == ENCODE_FIRST_PASS) && copy_frame)) {
        // Copy the picture buffer
        if (src->p_buffer != NULL)
            copy_frame_buffer(sequenceControlSet, dst->p_buffer, dst_y8b->p_buffer, src->p_buffer, pass);
    }

}
/***********************************************
**** Deep copy of the input metadata buffer
************************************************/
static EbErrorType copy_metadata_buffer(EbBufferHeaderType *dst, EbBufferHeaderType *src)
{
    EbErrorType return_error = EB_ErrorNone;
    for (size_t i = 0; i < src->metadata->sz; ++i) {
        SvtMetadataT *current_metadata = src->metadata->metadata_array[i];
        const uint32_t type = current_metadata->type;
        const uint8_t *payload = current_metadata->payload;
        const size_t sz = current_metadata->sz;

        if (svt_add_metadata(dst, type, payload, sz) != 0)
            SVT_WARN("Metadata of type %d could not be added to the buffer.\n", type);
    }
    return return_error;
}
/**********************************
* Empty This Buffer
**********************************/
EB_API EbErrorType svt_av1_enc_send_picture(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbEncHandle          *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr;
    EbBufferHeaderType   *app_hdr = p_buffer;

    // Get new Luma-8b buffer & a new (Chroma-8b + Luma-Chroma-2bit) buffers; Lib will release once done.
    EbObjectWrapper  *eb_y8b_wrapper_ptr;
    svt_get_empty_object(
        enc_handle_ptr->input_y8b_buffer_producer_fifo_ptr,
        &eb_y8b_wrapper_ptr);
    //set live count to 1 to be decremented at the end of the encode in RC
    svt_object_inc_live_count(eb_y8b_wrapper_ptr, 1);

   // svt_object_inc_live_count(eb_y8b_wrapper_ptr, 1);

    svt_get_empty_object(
        enc_handle_ptr->input_buffer_producer_fifo_ptr,
        &eb_wrapper_ptr);

     //set live count to 1 to be decremented at the end of the encode in RC, and released
     //this would also allow low delay TF to retain pictures
     svt_object_inc_live_count(eb_wrapper_ptr, 1);

    if (p_buffer != NULL) {


        //copy the Luma 8bit part into y8b buffer and the rest of samples into the regular buffer
        EbBufferHeaderType *lib_y8b_hdr = (EbBufferHeaderType*)eb_y8b_wrapper_ptr->object_ptr;
        EbBufferHeaderType *lib_reg_hdr = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;
        copy_input_buffer(
            enc_handle_ptr->scs_instance_array[0]->scs_ptr,
            lib_reg_hdr,
            lib_y8b_hdr,
            app_hdr,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.pass == ENC_FIRST_PASS);
    }

    //Take a new App-RessCoord command
    EbObjectWrapper *input_cmd_wrp;
    svt_get_empty_object(
        enc_handle_ptr->input_cmd_producer_fifo_ptr,
        &input_cmd_wrp);
    InputCommand *input_cmd_obj = (InputCommand*)input_cmd_wrp->object_ptr;
    //Fill the command with two picture buffers
    input_cmd_obj->eb_input_wrapper_ptr = eb_wrapper_ptr;
    input_cmd_obj->eb_y8b_wrapper_ptr = eb_y8b_wrapper_ptr;
    //Send to Lib
    svt_post_full_object(input_cmd_wrp);

    return EB_ErrorNone;
}
static void copy_output_recon_buffer(
    EbBufferHeaderType   *dst,
    EbBufferHeaderType   *src
)
{
    // copy output Bitstream fileds
    dst->size = src->size;
    dst->n_alloc_len = src->n_alloc_len;
    dst->n_filled_len = src->n_filled_len;
    dst->p_app_private = src->p_app_private;
    dst->n_tick_count = src->n_tick_count;
    dst->pts = src->pts;
    dst->dts = src->dts;
    dst->flags = src->flags;
    dst->pic_type = src->pic_type;

    // Copy the metadata array
    if (src->metadata)
        copy_metadata_buffer(dst, src);
    else
        dst->metadata = NULL;

    // Copy the picture buffer
    if (src->p_buffer)
        svt_memcpy(dst->p_buffer, src->p_buffer, src->n_filled_len);

    return;
}

/**********************************
* svt_av1_enc_get_packet sends out packet
**********************************/
EB_API EbErrorType svt_av1_enc_get_packet(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType  **p_buffer,
    unsigned char          pic_send_done)
{
    EbErrorType             return_error = EB_ErrorNone;
    EbEncHandle          *enc_handle = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr = NULL;
    EbBufferHeaderType    *packet;
    if (pic_send_done)
        svt_get_full_object(
            enc_handle->output_stream_buffer_consumer_fifo_ptr,
            &eb_wrapper_ptr);
    else
        svt_get_full_object_non_blocking(
            enc_handle->output_stream_buffer_consumer_fifo_ptr,
            &eb_wrapper_ptr);

    if (eb_wrapper_ptr) {
        packet = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;
        if ( packet->flags & 0xfffffff0 )
            return_error = EB_ErrorMax;
        // return the output stream buffer
        *p_buffer = packet;

        // save the wrapper pointer for the release
        (*p_buffer)->wrapper_ptr = (void*)eb_wrapper_ptr;
    }
    else
        return_error = EB_NoErrorEmptyQueue;
    return return_error;
}

EB_API void svt_av1_enc_release_out_buffer(
    EbBufferHeaderType  **p_buffer)
{
    if (p_buffer && (*p_buffer)->wrapper_ptr)
    {
        if((*p_buffer)->p_buffer)
           EB_FREE((*p_buffer)->p_buffer);
        // Release out put buffer back into the pool
        svt_release_object((EbObjectWrapper  *)(*p_buffer)->wrapper_ptr);
     }
    return;
}

/**********************************
* Fill This Buffer
**********************************/
EB_API EbErrorType svt_av1_get_recon(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbErrorType           return_error = EB_ErrorNone;
    EbEncHandle          *enc_handle = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr = NULL;

    if (enc_handle->scs_instance_array[0]->scs_ptr->static_config.recon_enabled) {
        svt_get_full_object_non_blocking(
            enc_handle->output_recon_buffer_consumer_fifo_ptr,
            &eb_wrapper_ptr);

        if (eb_wrapper_ptr) {
            EbBufferHeaderType* obj_ptr = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;
            copy_output_recon_buffer(
                p_buffer,
                obj_ptr);

            if (p_buffer->flags != EB_BUFFERFLAG_EOS && p_buffer->flags != 0)
                return_error = EB_ErrorMax;
            svt_release_object((EbObjectWrapper  *)eb_wrapper_ptr);
        }
        else
            return_error = EB_NoErrorEmptyQueue;
    }
    else {
        // recon is not enabled
        return_error = EB_ErrorMax;
    }

    return return_error;
}

/**********************************
* Encoder Error Handling
**********************************/
void lib_svt_encoder_send_error_exit(
    EbPtr                    hComponent,
    uint32_t                 error_code)
{
    EbComponentType      *svt_enc_component = (EbComponentType*)hComponent;
    EbEncHandle          *enc_handle = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr = NULL;
    EbBufferHeaderType    *output_packet;

    svt_get_empty_object(
        enc_handle->output_stream_buffer_consumer_fifo_ptr,
        &eb_wrapper_ptr);

    output_packet            = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;

    output_packet->size     = 0;
    output_packet->flags    = error_code;
    output_packet->p_buffer   = NULL;

    svt_post_full_object(eb_wrapper_ptr);
}

EB_API const char *svt_av1_get_version(void) {
    return SVT_AV1_CVS_VERSION;
}

EB_API void svt_av1_print_version(void) {
    SVT_INFO("-------------------------------------------\n");
    SVT_INFO("SVT [version]:\tSVT-AV1 Encoder Lib %s\n", SVT_AV1_CVS_VERSION);
    const char *compiler =
#if defined( _MSC_VER ) && (_MSC_VER >= 1930)
    "Visual Studio 2022"
#elif defined( _MSC_VER ) && (_MSC_VER >= 1920)
    "Visual Studio 2019"
#elif defined( _MSC_VER ) && (_MSC_VER >= 1910)
    "Visual Studio 2017"
#elif defined( _MSC_VER ) && (_MSC_VER >= 1900)
    "Visual Studio 2015"
#elif defined( _MSC_VER )
    "Visual Studio (old)"
#elif defined(__clang__)
    __VERSION__ "\t"
#elif defined(__GNUC__)
    "GCC " __VERSION__ "\t"
#else
    "unknown compiler"
#endif
    ;
    SVT_INFO("SVT [build]  :\t%s %zu bit\n", compiler,  sizeof(void*) * 8);
#if !REPRODUCIBLE_BUILDS
    SVT_INFO("LIB Build date: %s %s\n", __DATE__, __TIME__);
#endif
    SVT_INFO("-------------------------------------------\n");
}

/**********************************
* Encoder Handle Initialization
**********************************/
EbErrorType init_svt_av1_encoder_handle(
    EbComponentType * hComponent)
{
    EbErrorType       return_error = EB_ErrorNone;
    EbComponentType  *svt_enc_component = (EbComponentType*)hComponent;
    EbEncHandle      *handle;
    svt_av1_print_version();

    enc_switch_to_real_time();

    // Set Component Size & Version
    svt_enc_component->size = sizeof(EbComponentType);

    EB_NEW(handle, svt_enc_handle_ctor, svt_enc_component);
    svt_enc_component->p_component_private = handle;

    return return_error;
}

static EbErrorType allocate_frame_buffer(
    SequenceControlSet       *scs_ptr,
    EbBufferHeaderType       *input_buffer,
    EbBool                   noy8b)
{
    EbErrorType   return_error = EB_ErrorNone;
    EbPictureBufferDescInitData input_pic_buf_desc_init_data;
    EbSvtAv1EncConfiguration   * config = &scs_ptr->static_config;
    uint8_t is_16bit = config->encoder_bit_depth > 8 ? 1 : 0;

    input_pic_buf_desc_init_data.max_width =
        !(scs_ptr->max_input_luma_width % 8) ?
        scs_ptr->max_input_luma_width :
        scs_ptr->max_input_luma_width + (scs_ptr->max_input_luma_width % 8);

    input_pic_buf_desc_init_data.max_height =
        !(scs_ptr->max_input_luma_height % 8) ?
        scs_ptr->max_input_luma_height :
        scs_ptr->max_input_luma_height + (scs_ptr->max_input_luma_height % 8);

    input_pic_buf_desc_init_data.bit_depth = (EbBitDepthEnum)config->encoder_bit_depth;
    input_pic_buf_desc_init_data.color_format = (EbColorFormat)config->encoder_color_format;

    input_pic_buf_desc_init_data.left_padding = scs_ptr->left_padding;
    input_pic_buf_desc_init_data.right_padding = scs_ptr->right_padding;
    input_pic_buf_desc_init_data.top_padding = scs_ptr->top_padding;
    input_pic_buf_desc_init_data.bot_padding = scs_ptr->bot_padding;

    input_pic_buf_desc_init_data.split_mode = is_16bit ? EB_TRUE : EB_FALSE;

    input_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    input_pic_buf_desc_init_data.is_16bit_pipeline = 0;

    // Enhanced Picture Buffer
    {
        EbPictureBufferDesc* buf;
        if (!noy8b) {
            EB_NEW(
                buf,
                svt_picture_buffer_desc_ctor,
                (EbPtr)&input_pic_buf_desc_init_data);
        }
        else {
            EB_NEW(
                buf,
                svt_picture_buffer_desc_ctor_noy8b,
                (EbPtr)&input_pic_buf_desc_init_data);
        }
        input_buffer->p_buffer = (uint8_t*)buf;

    }

    return return_error;
}

/*
  allocate an input sample Luma-8bit buffer
*/
static EbErrorType allocate_y8b_frame_buffer(
    SequenceControlSet       *scs_ptr,
    EbBufferHeaderType        *input_buffer)
{
    EbErrorType   return_error = EB_ErrorNone;
    EbPictureBufferDescInitData input_pic_buf_desc_init_data;
    EbSvtAv1EncConfiguration   * config = &scs_ptr->static_config;
    uint8_t is_16bit = 0;

    input_pic_buf_desc_init_data.max_width =
        !(scs_ptr->max_input_luma_width % 8) ?
        scs_ptr->max_input_luma_width :
        scs_ptr->max_input_luma_width + (scs_ptr->max_input_luma_width % 8);

    input_pic_buf_desc_init_data.max_height =
        !(scs_ptr->max_input_luma_height % 8) ?
        scs_ptr->max_input_luma_height :
        scs_ptr->max_input_luma_height + (scs_ptr->max_input_luma_height % 8);
    input_pic_buf_desc_init_data.bit_depth = EB_8BIT;
    input_pic_buf_desc_init_data.color_format = (EbColorFormat)config->encoder_color_format;

    input_pic_buf_desc_init_data.left_padding = scs_ptr->left_padding;
    input_pic_buf_desc_init_data.right_padding = scs_ptr->right_padding;
    input_pic_buf_desc_init_data.top_padding = scs_ptr->top_padding;
    input_pic_buf_desc_init_data.bot_padding = scs_ptr->bot_padding;

    input_pic_buf_desc_init_data.split_mode = is_16bit ? EB_TRUE : EB_FALSE;

    input_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK; //allocate for 8bit Luma only
    input_pic_buf_desc_init_data.is_16bit_pipeline = 0;


    // Enhanced Picture Buffer
    {
        EbPictureBufferDesc* buf;
        EB_NEW(
            buf,
            svt_picture_buffer_desc_ctor,
            (EbPtr)&input_pic_buf_desc_init_data);
        input_buffer->p_buffer = (uint8_t*)buf;


    }

    return return_error;
}
/*
  create a luma 8bit buffer descriptor
*/
EbErrorType svt_input_y8b_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{
    EbBufferHeaderType* input_buffer;
    SequenceControlSet        *scs_ptr = (SequenceControlSet*)object_init_data_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(input_buffer, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)input_buffer;
    // Initialize Header
    input_buffer->size = sizeof(EbBufferHeaderType);

    EbErrorType return_error = allocate_y8b_frame_buffer(
        scs_ptr,
        input_buffer);
    if (return_error != EB_ErrorNone)
        return return_error;

    input_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}
/*
  free a luma 8bit buffer descriptor
*/
void svt_input_y8b_destroyer(EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EbPictureBufferDesc* buf = (EbPictureBufferDesc*)obj->p_buffer;
    if (buf) {
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_y);
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cb);
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cr);
    }

    EB_DELETE(buf);
    EB_FREE(obj);
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType svt_input_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{
    EbBufferHeaderType* input_buffer;
    SequenceControlSet        *scs_ptr = (SequenceControlSet*)object_init_data_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(input_buffer, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)input_buffer;
    // Initialize Header
    input_buffer->size = sizeof(EbBufferHeaderType);

    EbErrorType return_error = allocate_frame_buffer(
        scs_ptr,
        input_buffer,
        EB_TRUE);
    if (return_error != EB_ErrorNone)
        return return_error;

    input_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

void svt_input_buffer_header_destroyer(    EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EbPictureBufferDesc* buf = (EbPictureBufferDesc*)obj->p_buffer;
    if (buf) {
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_y);
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cb);
        EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cr);
    }

    EB_DELETE(buf);
    EB_FREE(obj);
}

EbErrorType svt_overlay_buffer_header_creator(
    EbPtr* object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{
    EbBufferHeaderType* input_buffer;
    SequenceControlSet* scs_ptr = (SequenceControlSet*)object_init_data_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(input_buffer, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)input_buffer;
    // Initialize Header
    input_buffer->size = sizeof(EbBufferHeaderType);

    EbErrorType return_error = allocate_frame_buffer(
        scs_ptr,
        input_buffer,
        EB_FALSE);
    if (return_error != EB_ErrorNone)
        return return_error;

    input_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType svt_output_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    (void)object_init_data_ptr;
    EbBufferHeaderType* out_buf_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(out_buf_ptr, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)out_buf_ptr;

    // Initialize Header
    out_buf_ptr->size = sizeof(EbBufferHeaderType);
    // p_buffer and n_alloc_len are dynamically set in EbPacketizationProcess
    // out_buf_ptr->n_alloc_len;
    out_buf_ptr->p_app_private = NULL;

    return EB_ErrorNone;
}

void svt_output_buffer_header_destroyer(    EbPtr p)
{
    EbBufferHeaderType* obj = (EbBufferHeaderType*)p;
    EB_FREE(obj);
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType svt_output_recon_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{
    EbBufferHeaderType         *recon_buffer;
    SequenceControlSet        *scs_ptr = (SequenceControlSet*)object_init_data_ptr;
    const uint32_t luma_size =
        scs_ptr->seq_header.max_frame_width    *
        scs_ptr->seq_header.max_frame_height;
    // both u and v
    const uint32_t chroma_size = luma_size >> 1;
    const uint32_t ten_bit = (scs_ptr->static_config.encoder_bit_depth > 8);
    const uint32_t frame_size = (luma_size + chroma_size) << ten_bit;

    *object_dbl_ptr = NULL;
    EB_CALLOC(recon_buffer, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)recon_buffer;

    // Initialize Header
    recon_buffer->size = sizeof(EbBufferHeaderType);

    // Assign the variables
    EB_MALLOC(recon_buffer->p_buffer, frame_size);

    recon_buffer->n_alloc_len = frame_size;
    recon_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

void svt_output_recon_buffer_header_destroyer(    EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EB_FREE(obj->p_buffer);
    EB_FREE(obj);
}

/**********************************
* svt_av1_enc_get_stream_info get stream information from encoder
**********************************/
EB_API EbErrorType svt_av1_enc_get_stream_info(EbComponentType *    svt_enc_component,
                                    uint32_t stream_info_id, void* info)
{
    if (stream_info_id >= SVT_AV1_STREAM_INFO_END || stream_info_id < SVT_AV1_STREAM_INFO_START) {
        return EB_ErrorBadParameter;
    }
    EbEncHandle         *enc_handle = (EbEncHandle*)svt_enc_component->p_component_private;
    if (stream_info_id == SVT_AV1_STREAM_INFO_FIRST_PASS_STATS_OUT) {
        EncodeContext*      context = enc_handle->scs_instance_array[0]->encode_context_ptr;
        SvtAv1FixedBuf*     first_pass_stats = (SvtAv1FixedBuf*)info;
        first_pass_stats->buf = context->stats_out.stat;
        first_pass_stats->sz = context->stats_out.size * sizeof(FIRSTPASS_STATS);
        return EB_ErrorNone;
    }
    return EB_ErrorBadParameter;
}
// clang-format on
