// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
// SUMMARY
//   Contains the API component functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "EbThreads.h"
#include "EbUtility.h"
#include "EbString.h"
#include "EbEncHandle.h"
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
#ifdef ARCH_X86
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
#define EB_OUTPUTRECONBUFFERSIZE                                        (MAX_PICTURE_WIDTH_SIZE*MAX_PICTURE_HEIGHT_SIZE*2)   // Recon Slice Size
#define EB_OUTPUTSTATISTICSBUFFERSIZE                                   0x30            // 6X8 (8 Bytes for Y, U, V, number of bits, picture number, QP)
#define EOS_NAL_BUFFER_SIZE                                             0x0010 // Bitstream used to code EOS NAL
#define EB_OUTPUTSTREAMBUFFERSIZE_MACRO(ResolutionSize)                ((ResolutionSize) < (INPUT_SIZE_1080i_TH) ? 0x1E8480 : (ResolutionSize) < (INPUT_SIZE_1080p_TH) ? 0x2DC6C0 : (ResolutionSize) < (INPUT_SIZE_4K_TH) ? 0x2DC6C0 : 0x2DC6C0  )

#define ENCDEC_INPUT_PORT_MDC                                0
#define ENCDEC_INPUT_PORT_ENCDEC                             1
#define ENCDEC_INPUT_PORT_INVALID                           -1

#define SCD_LAD                                              6

/**************************************
 * Globals
 **************************************/

uint8_t                          num_groups = 0;
#ifdef _WIN32
GROUP_AFFINITY                   group_affinity;
EbBool                           alternate_groups = 0;
#elif defined(__linux__)
cpu_set_t                        group_affinity;
typedef struct logicalProcessorGroup {
    uint32_t num;
    uint32_t group[1024];
}processorGroup;
#define INITIAL_PROCESSOR_GROUP 16
processorGroup                  *lp_group = NULL;
#endif

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
    const char* PROCESSORID = "processor";
    const char* PHYSICALID = "physical id";
    int processor_id_len = EB_STRLEN(PROCESSORID, 128);
    int physical_id_len = EB_STRLEN(PHYSICALID, 128);
    int maxSize = INITIAL_PROCESSOR_GROUP;
    if (processor_id_len < 0 || processor_id_len >= 128)
        return EB_ErrorInsufficientResources;
    if (physical_id_len < 0 || physical_id_len >= 128)
        return EB_ErrorInsufficientResources;
    memset(lp_group, 0, INITIAL_PROCESSOR_GROUP * sizeof(processorGroup));

    FILE *fin = fopen("/proc/cpuinfo", "r");
    if (fin) {
        int processor_id = 0, socket_id = 0;
        char line[1024];
        while (fgets(line, sizeof(line), fin)) {
            if(strncmp(line, PROCESSORID, processor_id_len) == 0) {
                char* p = line + processor_id_len;
                while(*p < '0' || *p > '9') p++;
                processor_id = strtol(p, NULL, 0);
            }
            if(strncmp(line, PHYSICALID, physical_id_len) == 0) {
                char* p = line + physical_id_len;
                while(*p < '0' || *p > '9') p++;
                socket_id = strtol(p, NULL, 0);
                if (socket_id < 0) {
                    fclose(fin);
                    return EB_ErrorInsufficientResources;
                }
                if (socket_id + 1 > num_groups)
                    num_groups = socket_id + 1;
                if (socket_id >= maxSize) {
                    maxSize = maxSize * 2;
                    lp_group = realloc(lp_group, maxSize * sizeof(processorGroup));
                    if (lp_group == NULL)
                        return EB_ErrorInsufficientResources;
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

void eb_set_thread_management_parameters(EbSvtAv1EncConfiguration *config_ptr)
{
#ifdef _WIN32
    uint32_t num_logical_processors = get_num_processors();
    // For system with a single processor group(no more than 64 logic processors all together)
    // Affinity of the thread can be set to one or more logical processors
    if (config_ptr->logical_processors == 1 && config_ptr->unpin_lp1 == 1) {
        group_affinity.Mask = get_affinity_mask(num_logical_processors);
    }
    else {
        if (num_groups == 1) {
            uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
                config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
            group_affinity.Mask = get_affinity_mask(lps);
        }
        else if (num_groups > 1) { // For system with multiple processor group
            if (config_ptr->logical_processors == 0) {
                if (config_ptr->target_socket != -1)
                    group_affinity.Group = config_ptr->target_socket;
            }
            else {
                uint32_t num_lp_per_group = num_logical_processors / num_groups;
                if (config_ptr->target_socket == -1) {
                    if (config_ptr->logical_processors > num_lp_per_group) {
                        alternate_groups = EB_TRUE;
                        SVT_LOG("SVT [WARNING]: -lp(logical processors) setting is ignored. Run on both sockets. \n");
                    }
                    else
                        group_affinity.Mask = get_affinity_mask(config_ptr->logical_processors);
                }
                else {
                    uint32_t lps = config_ptr->logical_processors == 0 ? num_lp_per_group :
                        config_ptr->logical_processors < num_lp_per_group ? config_ptr->logical_processors : num_lp_per_group;
                    group_affinity.Mask = get_affinity_mask(lps);
                    group_affinity.Group = config_ptr->target_socket;
                }
            }
        }
    }
#elif defined(__linux__)
    uint32_t num_logical_processors = get_num_processors();
    if (config_ptr->logical_processors == 1 && config_ptr->unpin_lp1 == 1) {
        pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &group_affinity);
    }
    else {
        CPU_ZERO(&group_affinity);

        if (num_groups == 1) {
            uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
                config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
            for (uint32_t i = 0; i < lps; i++)
                CPU_SET(lp_group[0].group[i], &group_affinity);
        }
        else if (num_groups > 1) {
            uint32_t num_lp_per_group = num_logical_processors / num_groups;
            if (config_ptr->logical_processors == 0) {
                if (config_ptr->target_socket != -1) {
                    for (uint32_t i = 0; i < lp_group[config_ptr->target_socket].num; i++)
                        CPU_SET(lp_group[config_ptr->target_socket].group[i], &group_affinity);
                }
            }
            else {
                if (config_ptr->target_socket == -1) {
                    uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
                        config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
                    if (lps > num_lp_per_group) {
                        for (uint32_t i = 0; i < lp_group[0].num; i++)
                            CPU_SET(lp_group[0].group[i], &group_affinity);
                        for (uint32_t i = 0; i < (lps - lp_group[0].num); i++)
                            CPU_SET(lp_group[1].group[i], &group_affinity);
                    }
                    else {
                        for (uint32_t i = 0; i < lps; i++)
                            CPU_SET(lp_group[0].group[i], &group_affinity);
                    }
                }
                else {
                    uint32_t lps = config_ptr->logical_processors == 0 ? num_lp_per_group :
                        config_ptr->logical_processors < num_lp_per_group ? config_ptr->logical_processors : num_lp_per_group;
                    for (uint32_t i = 0; i < lps; i++)
                        CPU_SET(lp_group[config_ptr->target_socket].group[i], &group_affinity);
                }
            }
        }
    }
#else
    UNUSED(config_ptr);
#endif
}

void asm_set_convolve_asm_table(void);
void asm_set_convolve_hbd_asm_table(void);
void init_intra_dc_predictors_c_internal(void);
void init_intra_predictors_internal(void);
void eb_av1_init_me_luts(void);

void switch_to_real_time(){
#ifndef _WIN32

    struct sched_param schedParam = {
        .sched_priority = 99
    };

    int32_t retValue = pthread_setschedparam(pthread_self(), SCHED_FIFO, &schedParam);
    UNUSED(retValue);
#endif
}
#define SINGLE_CORE_COUNT       1
#define CONS_CORE_COUNT         16
#define LOW_SERVER_CORE_COUNT   48
#define MED_SERVER_CORE_COUNT   128
#define HIGH_SERVER_CORE_COUNT  224

int32_t set_parent_pcs(EbSvtAv1EncConfiguration*   config, uint32_t core_count, EbInputResolution res_class) {
    if (config){
        uint32_t fps            = (uint32_t)((config->frame_rate > 1000) ?
                        config->frame_rate >> 16 :
                        config->frame_rate);
        uint32_t ppcs_count     = fps;
        uint32_t min_ppcs_count = (2 << config->hierarchical_levels) + 1; // min picture count to start encoding
        fps        = fps > 120 ? 120   : fps;
        fps        = fps < 24  ? 24    : fps;

        ppcs_count = MAX(min_ppcs_count, fps);
        if (core_count <= SINGLE_CORE_COUNT)
            ppcs_count = min_ppcs_count;
        else{
            if (res_class < INPUT_SIZE_1080i_RANGE){
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
        SVT_LOG("SVT[error]: Configuration struct is corrupted\n");
        return -1;
    }
}
EbErrorType load_default_buffer_configuration_settings(
    SequenceControlSet       *scs_ptr){
    EbErrorType           return_error = EB_ErrorNone;
#if !MD_RATE_EST_ENH
    uint32_t enc_dec_seg_h = (scs_ptr->static_config.super_block_size == 128) ?
        ((scs_ptr->max_input_luma_height + 64) / 128) :
        ((scs_ptr->max_input_luma_height + 32) / 64);
    uint32_t enc_dec_seg_w = (scs_ptr->static_config.super_block_size == 128) ?
        ((scs_ptr->max_input_luma_width + 64) / 128) :
        ((scs_ptr->max_input_luma_width + 32) / 64);
    uint32_t me_seg_h     = (((scs_ptr->max_input_luma_height + 32) / BLOCK_SIZE_64) < 6) ? 1 : 6;
    uint32_t me_seg_w     = (((scs_ptr->max_input_luma_width + 32) / BLOCK_SIZE_64) < 10) ? 1 : 10;
#endif
    unsigned int lp_count   = get_num_processors();
    unsigned int core_count = lp_count;
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
    uint32_t input_pic = (uint32_t)return_ppcs;
    scs_ptr->input_buffer_fifo_init_count = input_pic + SCD_LAD + scs_ptr->static_config.look_ahead_distance;
    scs_ptr->output_stream_buffer_fifo_init_count =
        scs_ptr->input_buffer_fifo_init_count + 4;
#if MD_RATE_EST_ENH
    uint32_t enc_dec_seg_h = (core_count == SINGLE_CORE_COUNT) ? 1 :
        (scs_ptr->static_config.super_block_size == 128) ?
        ((scs_ptr->max_input_luma_height + 64) / 128) :
        ((scs_ptr->max_input_luma_height + 32) / 64);
    uint32_t enc_dec_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
        (scs_ptr->static_config.super_block_size == 128) ?
        ((scs_ptr->max_input_luma_width + 64) / 128) :
        ((scs_ptr->max_input_luma_width + 32) / 64);
    uint32_t me_seg_h = (core_count == SINGLE_CORE_COUNT) ? 1 :
        (((scs_ptr->max_input_luma_height + 32) / BLOCK_SIZE_64) < 6) ? 1 : 6;
    uint32_t me_seg_w = (core_count == SINGLE_CORE_COUNT) ? 1 :
        (((scs_ptr->max_input_luma_width + 32) / BLOCK_SIZE_64) < 10) ? 1 : 10;
#endif
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

#if TILES_PARALLEL
    // Jing:
    // A tile group can be consisted by 1 tile or NxM tiles.
    // Segments will be parallelized within a tile group
    // We can use tile group to control the threads/parallelism in ED stage
    // NOTE:1 col will have better perf for segments for large resolutions
    uint8_t tile_group_col_count = 1;//(1 << scs_ptr->static_config.tile_columns)
    uint8_t tile_group_row_count = (1 << scs_ptr->static_config.tile_rows);

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
#endif
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

    scs_ptr->cdef_segment_column_count = me_seg_w;
    scs_ptr->cdef_segment_row_count    = me_seg_h;

    //since restoration unit size is same for Luma and Chroma, Luma segments and chroma segments do not correspond to the same area!
    //to keep proper processing, segments have to be configured based on chroma resolution.
    uint32_t unit_size                                  = 256;
    uint32_t rest_seg_w                                 = MAX((scs_ptr->max_input_luma_width /2 + (unit_size >> 1)) / unit_size, 1);
    uint32_t rest_seg_h                                 = MAX((scs_ptr->max_input_luma_height/2 + (unit_size >> 1)) / unit_size, 1);
    scs_ptr->rest_segment_column_count = MIN(rest_seg_w,6);
    scs_ptr->rest_segment_row_count    = MIN(rest_seg_h,4);

    scs_ptr->tf_segment_column_count = me_seg_w;//1;//
    scs_ptr->tf_segment_row_count =  me_seg_h;//1;//
    //#====================== Data Structures and Picture Buffers ======================
    scs_ptr->picture_control_set_pool_init_count       = input_pic + SCD_LAD + scs_ptr->static_config.look_ahead_distance;
    if (scs_ptr->static_config.enable_overlays)
        scs_ptr->picture_control_set_pool_init_count = MAX(scs_ptr->picture_control_set_pool_init_count,
            scs_ptr->static_config.look_ahead_distance + // frames in the LAD
            scs_ptr->static_config.look_ahead_distance / (1 << scs_ptr->static_config.hierarchical_levels) + 1 +  // number of overlayes in the LAD
            ((1 << scs_ptr->static_config.hierarchical_levels) + SCD_LAD) * 2 +// minigop formation in PD + SCD_LAD *(normal pictures + potential pictures )
            (1 << scs_ptr->static_config.hierarchical_levels)); // minigop in PM
    scs_ptr->picture_control_set_pool_init_count_child = MAX(MAX(MIN(3, core_count/2), core_count / 6), 1);
    scs_ptr->reference_picture_buffer_init_count       = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << scs_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          scs_ptr->static_config.look_ahead_distance + SCD_LAD;
    scs_ptr->pa_reference_picture_buffer_init_count    = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << scs_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          scs_ptr->static_config.look_ahead_distance + SCD_LAD;
    scs_ptr->output_recon_buffer_fifo_init_count       = scs_ptr->reference_picture_buffer_init_count;
    scs_ptr->overlay_input_picture_buffer_init_count   = scs_ptr->static_config.enable_overlays ?
                                                                          (2 << scs_ptr->static_config.hierarchical_levels) + SCD_LAD : 1;
    //Future frames window in Scene Change Detection (SCD) / TemporalFiltering
    scs_ptr->scd_delay =
        scs_ptr->static_config.enable_altrefs || scs_ptr->static_config.scene_change_detection ? SCD_LAD : 0;

    // bistream buffer will be allocated at run time. app will free the buffer once written to file.
    scs_ptr->output_stream_buffer_fifo_init_count = PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH;

    uint32_t min_input, min_parent, min_child, min_paref, min_ref, min_overlay;
    {
        /*Look-Ahead. Picture-Decision outputs pictures by group of mini-gops so
          the needed pictures for a certain look-ahead distance (LAD) should be rounded up to the next multiple of MiniGopSize.*/
        uint32_t mg_size = 1 << scs_ptr->static_config.hierarchical_levels;
        uint32_t needed_lad_pictures = ((scs_ptr->static_config.look_ahead_distance + mg_size - 1) / mg_size) * mg_size;

        /*To accomodate FFMPEG EOS, 1 frame delay is needed in Resource coordination.
           note that we have the option to not add 1 frame delay of Resource Coordination. In this case we have wait for first I frame
           to be released back to be able to start first base(16). Anyway poc16 needs to wait for poc0 to finish.*/
        uint32_t eos_delay = 1;

        //Minimum input pictures needed in the pipeline
        return_ppcs = (mg_size + 1) + eos_delay + scs_ptr->scd_delay + needed_lad_pictures;

        //scs_ptr->input_buffer_fifo_init_count = return_ppcs;
        min_input = return_ppcs;

        if (scs_ptr->static_config.enable_overlays)
            //scs_ptr->picture_control_set_pool_init_count =
            min_parent =  ((mg_size + 1) + eos_delay + scs_ptr->scd_delay) * 2 + //min to get a mini-gop in PD - release all and keep one.
                needed_lad_pictures +
                needed_lad_pictures / mg_size + 1;// Number of overlays in the look-ahead window
        else
           min_parent = return_ppcs;

        //Pic-Manager will inject one child at a time.
        min_child = 1;

        //References. Min to sustain flow (RA-5L-MRP-ON) 7 pictures from previous MGs + 10 needed for curr mini-GoP
        min_ref = 17;

        //Pa-References.Min to sustain flow (RA-5L-MRP-ON) -->TODO: derive numbers for other GOP Structures.
        min_paref = 25 +  scs_ptr->scd_delay + eos_delay;
        if (scs_ptr->static_config.enable_overlays)
            min_paref *= 2;

        //Overlays
        min_overlay = scs_ptr->static_config.enable_overlays ?
              mg_size + eos_delay + scs_ptr->scd_delay : 1;
    }

    if (core_count == SINGLE_CORE_COUNT) {
        scs_ptr->input_buffer_fifo_init_count                  = min_input;
        scs_ptr->picture_control_set_pool_init_count           = min_parent;
        scs_ptr->pa_reference_picture_buffer_init_count        = min_paref;
        scs_ptr->reference_picture_buffer_init_count           = min_ref;
        scs_ptr->picture_control_set_pool_init_count_child     = min_child;
        scs_ptr->overlay_input_picture_buffer_init_count       = min_overlay;

        scs_ptr->output_recon_buffer_fifo_init_count = scs_ptr->reference_picture_buffer_init_count;
    }
    else {
        scs_ptr->input_buffer_fifo_init_count              = MAX(min_input, scs_ptr->input_buffer_fifo_init_count);
        scs_ptr->picture_control_set_pool_init_count       = MAX(min_parent, scs_ptr->picture_control_set_pool_init_count);
        scs_ptr->pa_reference_picture_buffer_init_count    = MAX(min_paref, scs_ptr->pa_reference_picture_buffer_init_count);
        scs_ptr->reference_picture_buffer_init_count       = MAX(min_ref, scs_ptr->reference_picture_buffer_init_count);
        scs_ptr->picture_control_set_pool_init_count_child = MAX(min_child, scs_ptr->picture_control_set_pool_init_count_child);
        scs_ptr->overlay_input_picture_buffer_init_count   = MAX(min_overlay, scs_ptr->overlay_input_picture_buffer_init_count);
    }

    //#====================== Inter process Fifos ======================
    scs_ptr->resource_coordination_fifo_init_count       = 300;
    scs_ptr->picture_analysis_fifo_init_count            = 300;
    scs_ptr->picture_decision_fifo_init_count            = 300;
    scs_ptr->initial_rate_control_fifo_init_count        = 300;
    scs_ptr->picture_demux_fifo_init_count               = 300;
    scs_ptr->rate_control_tasks_fifo_init_count          = 300;
    scs_ptr->rate_control_fifo_init_count                = 301;
#if TILES_PARALLEL
    //Jing: Too many tiles may drain the fifo
    scs_ptr->mode_decision_configuration_fifo_init_count = 300 * (MIN(9, 1<<scs_ptr->static_config.tile_rows));
#else
    scs_ptr->mode_decision_configuration_fifo_init_count = 300;
#endif
    scs_ptr->motion_estimation_fifo_init_count           = 300;
    scs_ptr->entropy_coding_fifo_init_count              = 300;
    scs_ptr->enc_dec_fifo_init_count                     = 300;
    scs_ptr->dlf_fifo_init_count                         = 300;
    scs_ptr->cdef_fifo_init_count                        = 300;
    scs_ptr->rest_fifo_init_count                        = 300;
    //#====================== Processes number ======================
    scs_ptr->total_process_init_count                    = 0;
    if (core_count > 1){
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = MAX(MIN(15, core_count >> 1), core_count / 6));
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count =  MAX(MIN(20, core_count >> 1), core_count / 3));//1);//
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = MAX(MIN(3, core_count >> 1), core_count / 12));
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = MAX(MIN(3, core_count >> 1), core_count / 12));
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = MAX(MIN(40, core_count >> 1), core_count));
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = MAX(MIN(3, core_count >> 1), core_count / 12));
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = MAX(MIN(40, core_count >> 1), core_count));
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = MAX(MIN(40, core_count >> 1), core_count));
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = MAX(MIN(40, core_count >> 1), core_count));
    }else{
        scs_ptr->total_process_init_count += (scs_ptr->picture_analysis_process_init_count            = 1);
        scs_ptr->total_process_init_count += (scs_ptr->motion_estimation_process_init_count           = 1);
        scs_ptr->total_process_init_count += (scs_ptr->source_based_operations_process_init_count     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->mode_decision_configuration_process_init_count = 1);
        scs_ptr->total_process_init_count += (scs_ptr->enc_dec_process_init_count                     = 1);
        scs_ptr->total_process_init_count += (scs_ptr->entropy_coding_process_init_count              = 1);
        scs_ptr->total_process_init_count += (scs_ptr->dlf_process_init_count                         = 1);
        scs_ptr->total_process_init_count += (scs_ptr->cdef_process_init_count                        = 1);
        scs_ptr->total_process_init_count += (scs_ptr->rest_process_init_count                        = 1);
    }

    scs_ptr->total_process_init_count += 6; // single processes count
    SVT_LOG("Number of logical cores available: %u\nNumber of PPCS %u\n", core_count, scs_ptr->picture_control_set_pool_init_count);

    /******************************************************************
    * Platform detection, limit cpu flags to hardware available CPU
    ******************************************************************/
#ifdef ARCH_X86
    const CPU_FLAGS cpu_flags = get_cpu_flags();
    const CPU_FLAGS cpu_flags_to_use = get_cpu_flags_to_use();
    scs_ptr->static_config.use_cpu_flags &= cpu_flags_to_use;
    SVT_LOG("[asm level on system : up to %s]\n", get_asm_level_name_str(cpu_flags));
    SVT_LOG("[asm level selected : up to %s]\n", get_asm_level_name_str(scs_ptr->static_config.use_cpu_flags));
#else
    scs_ptr->static_config.use_cpu_flags &= 0;
    SVT_LOG("[asm level on system : up to %s]\n", get_asm_level_name_str(0));
    SVT_LOG("[asm level selected : up to %s]\n", get_asm_level_name_str(scs_ptr->static_config.use_cpu_flags));
#endif
    return return_error;
}
 // Rate Control
static RateControlPorts rate_control_ports[] = {
    {RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER,       0},
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

static void eb_enc_handle_stop_threads(EbEncHandle *enc_handle_ptr)
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
static void eb_enc_handle_dctor(EbPtr p)
{
    EbEncHandle *enc_handle_ptr = (EbEncHandle *)p;

    eb_enc_handle_stop_threads(enc_handle_ptr);
    EB_FREE_PTR_ARRAY(enc_handle_ptr->app_callback_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->scs_pool_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->pa_reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->overlay_input_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->input_buffer_resource_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->output_stream_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->output_recon_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->resource_coordination_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_analysis_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_decision_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->motion_estimation_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->initial_rate_control_results_resource_ptr);
    EB_DELETE(enc_handle_ptr->picture_demux_results_resource_ptr);
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
static EbErrorType eb_enc_handle_ctor(
    EbEncHandle *enc_handle_ptr,
    EbComponentType * ebHandlePtr)
{
    EbErrorType return_error = EB_ErrorNone;

    enc_handle_ptr->dctor = eb_enc_handle_dctor;

    return_error = init_thread_management_params();

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
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
    EB_NEW(enc_handle_ptr->scs_instance_array[0], eb_sequence_control_set_instance_ctor);
    return EB_ErrorNone;
}

EbErrorType eb_input_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

EbErrorType eb_output_recon_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

EbErrorType eb_output_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);

void eb_input_buffer_header_destroyer(    EbPtr p);
void eb_output_recon_buffer_header_destroyer(    EbPtr p);
void eb_output_buffer_header_destroyer(    EbPtr p);


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

void init_fn_ptr(void);
void av1_init_wedge_masks(void);
/**********************************
* Initialize Encoder Library
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_init(EbComponentType *svt_enc_component)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;
    EbEncHandle *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbErrorType return_error = EB_ErrorNone;
    uint32_t instance_index;
    uint32_t process_index;
    uint32_t max_picture_width;
    EbBool is_16bit = (EbBool)(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.encoder_color_format;
    SequenceControlSet* control_set_ptr;

    setup_common_rtcd_internal(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.use_cpu_flags);
    setup_rtcd_internal(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.use_cpu_flags);
    asm_set_convolve_asm_table();

    init_intra_dc_predictors_c_internal();

    asm_set_convolve_hbd_asm_table();

    init_intra_predictors_internal();
    EbSequenceControlSetInitData scs_init;
    scs_init.sb_size = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.super_block_size;

    build_blk_geom(scs_init.sb_size == 128);

    eb_av1_init_me_luts();
    init_fn_ptr();
    av1_init_wedge_masks();
    /************************************
    * Sequence Control Set
    ************************************/
    EB_NEW(enc_handle_ptr->scs_pool_ptr,
        eb_system_resource_ctor,
        enc_handle_ptr->scs_pool_total_count,
        1,
        0,
        eb_sequence_control_set_creator,
        &scs_init,
        NULL);

    /************************************
    * Picture Control Set: Parent
    ************************************/
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        // The segment Width & Height Arrays are in units of SBs, not samples
        PictureControlSetInitData input_data;

        input_data.picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        input_data.picture_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
        input_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->left_padding;
        input_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->right_padding;
        input_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->top_padding;
        input_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->bot_padding;
        input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        input_data.color_format = color_format;
        input_data.sb_sz = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz;
        input_data.max_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_sb_depth;
        input_data.ten_bit_format = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.ten_bit_format;
        input_data.compressed_ten_bit_format = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.compressed_ten_bit_format;
        input_data.enc_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enc_mode;
        input_data.speed_control = (uint8_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.speed_control_flag;
        input_data.hbd_mode_decision = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enable_hbd_mode_decision;
        input_data.film_grain_noise_level = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.film_grain_denoise_strength;
        input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.encoder_bit_depth;
        input_data.ext_block_flag = (uint8_t)enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.ext_block_flag;
        input_data.mrp_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->mrp_mode;
        input_data.nsq_present = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->nsq_present;
#if TILES_PARALLEL
        input_data.log2_tile_rows = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.tile_rows;
        input_data.log2_tile_cols = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.tile_columns;
        input_data.log2_sb_sz = (scs_init.sb_size == 128) ? 5 : 4;
#endif
#if OIS_MEM
        input_data.allocate_ois_struct = 0;
#endif
        input_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.is_16bit_pipeline;

        EB_NEW(
            enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->picture_control_set_pool_init_count,//enc_handle_ptr->pcs_pool_total_count,
            1,
            0,
            picture_parent_control_set_creator,
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

        input_data.picture_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        input_data.picture_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
        input_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->left_padding;
        input_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->right_padding;
        input_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->top_padding;
        input_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->bot_padding;
        input_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        input_data.film_grain_noise_level = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->film_grain_denoise_strength;
        input_data.color_format = color_format;
        input_data.sb_sz = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz;
        input_data.sb_size_pix = scs_init.sb_size;
        input_data.max_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_sb_depth;
        input_data.hbd_mode_decision = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.enable_hbd_mode_decision;
        input_data.cdf_mode = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->cdf_mode;
        input_data.mfmv = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->mfmv_enabled;
        input_data.cfg_palette = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.screen_content_mode;
#if TILES_PARALLEL
        //Jing: Get tile info from parent_pcs
        PictureParentControlSet *parent_pcs = (PictureParentControlSet *)enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index]->wrapper_ptr_pool[0]->object_ptr;
        input_data.tile_row_count = parent_pcs->av1_cm->tiles_info.tile_rows;
        input_data.tile_column_count = parent_pcs->av1_cm->tiles_info.tile_cols;
#endif
        input_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.is_16bit_pipeline;
        input_data.serial_rate_est = enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.pic_based_rate_est &&
            input_data.enc_dec_segment_col == 1 && input_data.enc_dec_segment_row == 1 ? 1 : 0;
        EB_NEW(
            enc_handle_ptr->picture_control_set_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
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

    // Rate Control
    rate_control_ports[0].count = EB_PictureManagerProcessInitCount;
    rate_control_ports[1].count = EB_PacketizationProcessInitCount;
    rate_control_ports[2].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->entropy_coding_process_init_count;
    rate_control_ports[3].count = 0;

    enc_dec_ports[ENCDEC_INPUT_PORT_MDC].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->mode_decision_configuration_process_init_count;
    enc_dec_ports[ENCDEC_INPUT_PORT_ENCDEC].count = enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count;

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EbReferenceObjectDescInitData     eb_ref_obj_ect_desc_init_data_structure;
        EbPaReferenceObjectDescInitData   eb_pa_ref_obj_ect_desc_init_data_structure;
        EbPictureBufferDescInitData       ref_pic_buf_desc_init_data;
        EbPictureBufferDescInitData       quart_pic_buf_desc_init_data;
        EbPictureBufferDescInitData       sixteenth_pic_buf_desc_init_data;
        // Initialize the various Picture types
        ref_pic_buf_desc_init_data.max_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        ref_pic_buf_desc_init_data.max_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
        ref_pic_buf_desc_init_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        ref_pic_buf_desc_init_data.color_format = color_format;
        ref_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

        ref_pic_buf_desc_init_data.left_padding = PAD_VALUE;
        ref_pic_buf_desc_init_data.right_padding = PAD_VALUE;
        ref_pic_buf_desc_init_data.top_padding = PAD_VALUE;
        ref_pic_buf_desc_init_data.bot_padding = PAD_VALUE;
        ref_pic_buf_desc_init_data.mfmv = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->mfmv_enabled;
        ref_pic_buf_desc_init_data.is_16bit_pipeline = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->static_config.is_16bit_pipeline;
        // Hsan: split_mode is set @ eb_reference_object_ctor() as both unpacked reference and packed reference are needed for a 10BIT input; unpacked reference @ MD, and packed reference @ EP

        if (is_16bit)
            ref_pic_buf_desc_init_data.bit_depth = EB_10BIT;

        eb_ref_obj_ect_desc_init_data_structure.reference_picture_desc_init_data = ref_pic_buf_desc_init_data;

        // Reference Picture Buffers
        EB_NEW(
            enc_handle_ptr->reference_picture_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->reference_picture_buffer_init_count,//enc_handle_ptr->ref_pic_pool_total_count,
            EB_PictureManagerProcessInitCount,
            0,
            eb_reference_object_creator,
            &(eb_ref_obj_ect_desc_init_data_structure),
            NULL);

        // PA Reference Picture Buffers
        // Currently, only Luma samples are needed in the PA
        ref_pic_buf_desc_init_data.max_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
        ref_pic_buf_desc_init_data.max_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
        ref_pic_buf_desc_init_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        ref_pic_buf_desc_init_data.color_format = EB_YUV420; //use 420 for picture analysis
        ref_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        ref_pic_buf_desc_init_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz + ME_FILTER_TAP;
        ref_pic_buf_desc_init_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz + ME_FILTER_TAP;
        ref_pic_buf_desc_init_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz + ME_FILTER_TAP;
        ref_pic_buf_desc_init_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz + ME_FILTER_TAP;
        ref_pic_buf_desc_init_data.split_mode = EB_FALSE;
        quart_pic_buf_desc_init_data.max_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width >> 1;
        quart_pic_buf_desc_init_data.max_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height >> 1;
        quart_pic_buf_desc_init_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        quart_pic_buf_desc_init_data.color_format = EB_YUV420;
        quart_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        quart_pic_buf_desc_init_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 1;
        quart_pic_buf_desc_init_data.split_mode = EB_FALSE;
        quart_pic_buf_desc_init_data.down_sampled_filtered = (enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;

        sixteenth_pic_buf_desc_init_data.max_width = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width >> 2;
        sixteenth_pic_buf_desc_init_data.max_height = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height >> 2;
        sixteenth_pic_buf_desc_init_data.bit_depth = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->encoder_bit_depth;
        sixteenth_pic_buf_desc_init_data.color_format = EB_YUV420;
        sixteenth_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        sixteenth_pic_buf_desc_init_data.left_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.right_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.top_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.bot_padding = enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->sb_sz >> 2;
        sixteenth_pic_buf_desc_init_data.split_mode = EB_FALSE;
        sixteenth_pic_buf_desc_init_data.down_sampled_filtered = (enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;

        eb_pa_ref_obj_ect_desc_init_data_structure.reference_picture_desc_init_data = ref_pic_buf_desc_init_data;
        eb_pa_ref_obj_ect_desc_init_data_structure.quarter_picture_desc_init_data = quart_pic_buf_desc_init_data;
        eb_pa_ref_obj_ect_desc_init_data_structure.sixteenth_picture_desc_init_data = sixteenth_pic_buf_desc_init_data;
        // Reference Picture Buffers
        EB_NEW(enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->pa_reference_picture_buffer_init_count,
            EB_PictureDecisionProcessInitCount,
            0,
            eb_pa_reference_object_creator,
            &(eb_pa_ref_obj_ect_desc_init_data_structure),
            NULL);
        // Set the SequenceControlSet Picture Pool Fifo Ptrs
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->reference_picture_pool_fifo_ptr = eb_system_resource_get_producer_fifo(enc_handle_ptr->reference_picture_pool_ptr_array[instance_index], 0);
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->pa_reference_picture_pool_fifo_ptr = eb_system_resource_get_producer_fifo(enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index], 0);

        if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.enable_overlays) {
            // Overlay Input Picture Buffers
            EB_NEW(
                enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index],
                eb_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->overlay_input_picture_buffer_init_count,
                1,
                0,
                eb_input_buffer_header_creator,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr,
                eb_input_buffer_header_destroyer);
           // Set the SequenceControlSet Overlay input Picture Pool Fifo Ptrs
            enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->overlay_input_picture_pool_fifo_ptr = eb_system_resource_get_producer_fifo(enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index], 0);
        }
    }

    /************************************
    * System Resource Managers & Fifos
    ************************************/

    // EbBufferHeaderType Input
    EB_NEW(
        enc_handle_ptr->input_buffer_resource_ptr,
        eb_system_resource_ctor,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->input_buffer_fifo_init_count,
        1,
        EB_ResourceCoordinationProcessInitCount,
        eb_input_buffer_header_creator,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr,
        eb_input_buffer_header_destroyer);

    enc_handle_ptr->input_buffer_producer_fifo_ptr = eb_system_resource_get_producer_fifo(enc_handle_ptr->input_buffer_resource_ptr, 0);


    // EbBufferHeaderType Output Stream
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_stream_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EB_NEW(
            enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->output_stream_buffer_fifo_init_count,
            enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->total_process_init_count,//EB_PacketizationProcessInitCount,
            1,
            eb_output_buffer_header_creator,
            &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config,
            eb_output_buffer_header_destroyer);
    }
    enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr = eb_system_resource_get_consumer_fifo(enc_handle_ptr->output_stream_buffer_resource_ptr_array[0], 0);
    if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.recon_enabled) {
        // EbBufferHeaderType Output Recon
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_recon_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            EB_NEW(
                enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index],
                eb_system_resource_ctor,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->output_recon_buffer_fifo_init_count,
                enc_handle_ptr->scs_instance_array[instance_index]->scs_ptr->enc_dec_process_init_count,
                1,
                eb_output_recon_buffer_header_creator,
                enc_handle_ptr->scs_instance_array[0]->scs_ptr,
                eb_output_recon_buffer_header_destroyer);
        }
        enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr = eb_system_resource_get_consumer_fifo(enc_handle_ptr->output_recon_buffer_resource_ptr_array[0], 0);
    }

    // Resource Coordination Results
    {
        ResourceCoordinationResultInitData resource_coordination_result_init_data;

        EB_NEW(
            enc_handle_ptr->resource_coordination_results_resource_ptr,
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_decision_fifo_init_count,
            EB_PictureDecisionProcessInitCount,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->picture_demux_fifo_init_count,
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count + enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count + 1, // 1 for packetization
            EB_PictureManagerProcessInitCount,
            picture_results_creator,
            &picture_result_init_data,
            NULL);

    }

    // Rate Control Tasks
    {
        RateControlTasksInitData rate_control_tasks_init_data;

        EB_NEW(
            enc_handle_ptr->rate_control_tasks_resource_ptr,
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
            eb_system_resource_ctor,
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
        enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->stream_output_fifo_ptr     = eb_system_resource_get_producer_fifo(enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index], 0);
        if (enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.recon_enabled)
            enc_handle_ptr->scs_instance_array[instance_index]->encode_context_ptr->recon_output_fifo_ptr  = eb_system_resource_get_producer_fifo(enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index], 0);
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
            process_index);
    }

    // Picture Manager Context
    EB_NEW(
        enc_handle_ptr->picture_manager_context_ptr,
        picture_manager_context_ctor,
        enc_handle_ptr,
        rate_control_port_lookup(RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER, 0));

    // Rate Control Context
    EB_NEW(
        enc_handle_ptr->rate_control_context_ptr,
        rate_control_context_ctor,
        enc_handle_ptr);

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

    max_picture_width = 0;
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
            enc_dec_port_lookup(ENCDEC_INPUT_PORT_ENCDEC, process_index),
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count + process_index);
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

    for (process_index = 0; process_index < enc_handle_ptr->scs_instance_array[0]->scs_ptr->rest_process_init_count; ++process_index) {
        EB_NEW(
            enc_handle_ptr->rest_context_ptr_array[process_index],
            rest_context_ctor,
            enc_handle_ptr,
            process_index,
            1 + process_index);
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
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->source_based_operations_process_init_count +
            enc_handle_ptr->scs_instance_array[0]->scs_ptr->enc_dec_process_init_count);

    /************************************
    * Thread Handles
    ************************************/
    EbSvtAv1EncConfiguration   *config_ptr = &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config;

    eb_set_thread_management_parameters(config_ptr);

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
        enc_dec_kernel,
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

#if DISPLAY_MEMORY
    EB_MEMORY();
#endif
    eb_print_memory_usage();

    return return_error;
}

/**********************************
* DeInitialize Encoder Library
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_deinit(EbComponentType *svt_enc_component){
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;

    EbEncHandle *handle = (EbEncHandle*)svt_enc_component->p_component_private;
    if (handle) {
        eb_shutdown_process(handle->input_buffer_resource_ptr);
        eb_shutdown_process(handle->resource_coordination_results_resource_ptr);
        eb_shutdown_process(handle->picture_analysis_results_resource_ptr);
        eb_shutdown_process(handle->picture_decision_results_resource_ptr);
        eb_shutdown_process(handle->motion_estimation_results_resource_ptr);
        eb_shutdown_process(handle->initial_rate_control_results_resource_ptr);
        eb_shutdown_process(handle->picture_demux_results_resource_ptr);
        eb_shutdown_process(handle->rate_control_tasks_resource_ptr);
        eb_shutdown_process(handle->rate_control_results_resource_ptr);
        eb_shutdown_process(handle->enc_dec_tasks_resource_ptr);
        eb_shutdown_process(handle->enc_dec_results_resource_ptr);
        eb_shutdown_process(handle->entropy_coding_results_resource_ptr);
        eb_shutdown_process(handle->dlf_results_resource_ptr);
        eb_shutdown_process(handle->cdef_results_resource_ptr);
        eb_shutdown_process(handle->rest_results_resource_ptr);
    }

    return EB_ErrorNone;
}

EbErrorType eb_svt_enc_init_parameter(
    EbSvtAv1EncConfiguration * config_ptr);

EbErrorType init_svt_av1_encoder_handle(
    EbComponentType * hComponent);
/**********************************
* GetHandle
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_init_handle(
    EbComponentType** p_handle,               // Function to be called in the future for manipulating the component
    void*              p_app_data,
    EbSvtAv1EncConfiguration  *config_ptr)              // pointer passed back to the client during callbacks

{
    EbErrorType           return_error = EB_ErrorNone;
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
        SVT_LOG("Error: Component Struct Malloc Failed\n");
        return_error = EB_ErrorInsufficientResources;
    }
    // Init Component OS objects (threads, semaphores, etc.)
    // also links the various Component control functions
    return_error = init_svt_av1_encoder_handle(*p_handle);

    if (return_error == EB_ErrorNone) {
        ((EbComponentType*)(*p_handle))->p_application_private = p_app_data;
        return_error = eb_svt_enc_init_parameter(config_ptr);
    }
    if (return_error != EB_ErrorNone) {
        svt_av1_enc_deinit(*p_handle);
        free(*p_handle);
        *p_handle = NULL;
        return return_error;
    }
    eb_increase_component_count();
    return return_error;
}

/**********************************
* Encoder Componenet DeInit
**********************************/
EbErrorType eb_av1_enc_component_de_init(EbComponentType  *svt_enc_component)
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
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_deinit_handle(
    EbComponentType  *svt_enc_component)
{
    EbErrorType return_error = EB_ErrorNone;

    if (svt_enc_component) {
        return_error = eb_av1_enc_component_de_init(svt_enc_component);

        free(svt_enc_component);
#if  defined(__linux__)
        EB_FREE(lp_group);
#endif
        eb_decrease_component_count();
    }
    else
        return_error = EB_ErrorInvalidComponent;
    return return_error;
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
    int32_t min_ip                     = ((int)((fps) / mini_gop_size)*(mini_gop_size));
    int32_t max_ip                     = ((int)((fps + mini_gop_size) / mini_gop_size)*(mini_gop_size));

    intra_period = (ABS((fps - max_ip)) > ABS((fps - min_ip))) ? min_ip : max_ip;

    if (config->intra_refresh_type == 1)
        intra_period -= 1;

    return intra_period;
}

// Set configurations for the hardcoded parameters
void set_default_configuration_parameters(
    SequenceControlSet       *scs_ptr)
{
    // SB Definitions
    scs_ptr->sb_sz = MAX_SB_SIZE;
    scs_ptr->max_sb_depth = (uint8_t)EB_MAX_SB_DEPTH;
    scs_ptr->static_config.enable_adaptive_quantization = 2;

    return;
}

static uint32_t compute_default_look_ahead(
    EbSvtAv1EncConfiguration*   config){
    int32_t lad = 0;
    if (config->rate_control_mode == 0 || config->intra_period_length < 0)
        lad = (2 << config->hierarchical_levels)+1;
    else
        lad = config->intra_period_length;

    return lad;
}

// Only use the maximum look ahead needed if
static uint32_t cap_look_ahead_distance(
    EbSvtAv1EncConfiguration*   config){
    uint32_t lad = 0;

    if(config){
        uint32_t fps = config->frame_rate < 1000 ?
                      config->frame_rate :
                      config->frame_rate >> 16;
        uint32_t max_cqp_lad = (2 << config->hierarchical_levels) + 1;
        uint32_t max_rc_lad  = fps << 1;
        lad = config->look_ahead_distance;
        if (config->rate_control_mode == 0 && lad > max_cqp_lad)
            lad = max_cqp_lad;
        else if (config->rate_control_mode != 0 && lad > max_rc_lad)
            lad = max_rc_lad;
    }

    lad = lad > MAX_LAD ? MAX_LAD: lad; // clip to max allowed lad

    return lad;
}

void set_param_based_on_input(SequenceControlSet *scs_ptr)
{
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

    derive_input_resolution(
        &scs_ptr->input_resolution,
        scs_ptr->seq_header.max_frame_width*scs_ptr->seq_header.max_frame_height);
    // In two pass encoding, the first pass uses sb size=64
    if (scs_ptr->use_output_stat_file)
        scs_ptr->static_config.super_block_size = 64;
    else
        scs_ptr->static_config.super_block_size = (scs_ptr->static_config.enc_mode <= ENC_M3) ? 128 : 64;
    scs_ptr->static_config.super_block_size = (scs_ptr->static_config.rate_control_mode > 1) ? 64 : scs_ptr->static_config.super_block_size;
   // scs_ptr->static_config.hierarchical_levels = (scs_ptr->static_config.rate_control_mode > 1) ? 3 : scs_ptr->static_config.hierarchical_levels;
    // Configure the padding
    scs_ptr->left_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->top_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->right_padding = BLOCK_SIZE_64 + 4;
    scs_ptr->bot_padding = scs_ptr->static_config.super_block_size + 4;
    scs_ptr->static_config.enable_overlays = scs_ptr->static_config.enable_altrefs == EB_FALSE ||
        (scs_ptr->static_config.altref_nframes <= 1) ||
        (scs_ptr->static_config.rate_control_mode > 0) ||
        scs_ptr->static_config.encoder_bit_depth != EB_8BIT ?
        0 : scs_ptr->static_config.enable_overlays;

    //0: MRP Mode 0 (4,3)
    //1: MRP Mode 1 (2,2)
    scs_ptr->mrp_mode = (uint8_t)(scs_ptr->static_config.enc_mode <= ENC_M7) ? 0 : 1;
    //0: ON
    //1: OFF
    scs_ptr->cdf_mode = (uint8_t)(scs_ptr->static_config.enc_mode <= ENC_M6) ? 0 : 1;

    //0: NSQ absent
    //1: NSQ present
    scs_ptr->nsq_present = (uint8_t)(scs_ptr->static_config.enc_mode <= ENC_M5) ? 1 : 0;

    // Set down-sampling method     Settings
    // 0                            0: filtering
    // 1                            1: decimation
    if (scs_ptr->static_config.screen_content_mode == 1)
        if (scs_ptr->static_config.enc_mode <= ENC_M4)
            scs_ptr->down_sampling_method_me_search = ME_FILTERED_DOWNSAMPLED;
        else
            scs_ptr->down_sampling_method_me_search = ME_DECIMATED_DOWNSAMPLED;
    else
        if (scs_ptr->static_config.enc_mode <= ENC_M4)
            scs_ptr->down_sampling_method_me_search = ME_FILTERED_DOWNSAMPLED;
        else
            scs_ptr->down_sampling_method_me_search = ME_DECIMATED_DOWNSAMPLED;


    // Set over_boundary_block_mode     Settings
    // 0                            0: not allowed
    // 1                            1: allowed
    if (scs_ptr->static_config.over_bndry_blk == DEFAULT)
        if (scs_ptr->static_config.enc_mode <= ENC_M5)
            scs_ptr->over_boundary_block_mode = 1;
        else
            scs_ptr->over_boundary_block_mode = 0;
    else
        scs_ptr->over_boundary_block_mode = scs_ptr->static_config.over_bndry_blk;
    if (scs_ptr->static_config.enable_mfmv == DEFAULT)
        if (scs_ptr->static_config.screen_content_mode == 1)
            scs_ptr->mfmv_enabled = 0;
        else
            scs_ptr->mfmv_enabled = (uint8_t)(scs_ptr->static_config.enc_mode <= ENC_M1) ? 1 : 0;
    else
        scs_ptr->mfmv_enabled = scs_ptr->static_config.enable_mfmv;

    // Set hbd_mode_decision OFF for high encode modes or bitdepth < 10
    if (scs_ptr->static_config.encoder_bit_depth < 10)
        scs_ptr->static_config.enable_hbd_mode_decision = 0;
}

void copy_api_from_app(
    SequenceControlSet       *scs_ptr,
    EbSvtAv1EncConfiguration   *config_struct){
    uint32_t                  hme_region_index = 0;

    scs_ptr->max_input_luma_width = config_struct->source_width;
    scs_ptr->max_input_luma_height = config_struct->source_height;
    scs_ptr->frame_rate = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate;
    // SB Definitions
    scs_ptr->static_config.pred_structure = 2; // Hardcoded(Cleanup)
    scs_ptr->static_config.enable_qp_scaling_flag = 1;
    scs_ptr->max_blk_size = (uint8_t)64;
    scs_ptr->min_blk_size = (uint8_t)8;
    scs_ptr->max_intra_size = (uint8_t)32;
    scs_ptr->min_intra_size = (uint8_t)8;
    scs_ptr->max_ref_count = 1;

    // Padding Offsets
    scs_ptr->sb_sz = (uint8_t)((EbSvtAv1EncConfiguration*)config_struct)->sb_sz;
    scs_ptr->max_sb_depth = (uint8_t)((EbSvtAv1EncConfiguration*)config_struct)->partition_depth;
    scs_ptr->static_config.intra_period_length = ((EbSvtAv1EncConfiguration*)config_struct)->intra_period_length;
    scs_ptr->static_config.intra_refresh_type = ((EbSvtAv1EncConfiguration*)config_struct)->intra_refresh_type;
    scs_ptr->static_config.hierarchical_levels = ((EbSvtAv1EncConfiguration*)config_struct)->hierarchical_levels;
    scs_ptr->static_config.enc_mode = ((EbSvtAv1EncConfiguration*)config_struct)->enc_mode;
    scs_ptr->static_config.snd_pass_enc_mode = ((EbSvtAv1EncConfiguration*)config_struct)->snd_pass_enc_mode;
    scs_ptr->intra_period_length = scs_ptr->static_config.intra_period_length;
    scs_ptr->intra_refresh_type = scs_ptr->static_config.intra_refresh_type;
    scs_ptr->max_temporal_layers = scs_ptr->static_config.hierarchical_levels;
    scs_ptr->static_config.use_qp_file = ((EbSvtAv1EncConfiguration*)config_struct)->use_qp_file;
    scs_ptr->static_config.input_stat_file = ((EbSvtAv1EncConfiguration*)config_struct)->input_stat_file;
    scs_ptr->static_config.output_stat_file = ((EbSvtAv1EncConfiguration*)config_struct)->output_stat_file;
    scs_ptr->use_input_stat_file = scs_ptr->static_config.input_stat_file ? 1 : 0;
    scs_ptr->use_output_stat_file = scs_ptr->static_config.output_stat_file ? 1 : 0;
    // Deblock Filter
    scs_ptr->static_config.disable_dlf_flag = ((EbSvtAv1EncConfiguration*)config_struct)->disable_dlf_flag;

    // Local Warped Motion
    scs_ptr->static_config.enable_warped_motion = ((EbSvtAv1EncConfiguration*)config_struct)->enable_warped_motion;

    // Global motion
    scs_ptr->static_config.enable_global_motion = ((EbSvtAv1EncConfiguration*)config_struct)->enable_global_motion;

    // CDEF
    scs_ptr->static_config.cdef_mode = ((EbSvtAv1EncConfiguration*)config_struct)->cdef_mode;

    // Restoration filtering
    scs_ptr->static_config.enable_restoration_filtering = ((EbSvtAv1EncConfiguration*)config_struct)->enable_restoration_filtering;
    scs_ptr->static_config.sg_filter_mode = ((EbSvtAv1EncConfiguration*)config_struct)->sg_filter_mode;
    scs_ptr->static_config.wn_filter_mode = ((EbSvtAv1EncConfiguration*)config_struct)->wn_filter_mode;

    //combine class 12
    scs_ptr->static_config.combine_class_12             = ((EbSvtAv1EncConfiguration*)config_struct)->combine_class_12;
    // motion field motion vector
    scs_ptr->static_config.enable_mfmv                  = ((EbSvtAv1EncConfiguration*)config_struct)->enable_mfmv;
    // redundant block
    scs_ptr->static_config.enable_redundant_blk         = ((EbSvtAv1EncConfiguration*)config_struct)->enable_redundant_blk;
    // spatial sse in full loop
    scs_ptr->static_config.spatial_sse_fl               = ((EbSvtAv1EncConfiguration*)config_struct)->spatial_sse_fl;
    // subpel
    scs_ptr->static_config.enable_subpel                = ((EbSvtAv1EncConfiguration*)config_struct)->enable_subpel;
    // over boundry block mode
    scs_ptr->static_config.over_bndry_blk               = ((EbSvtAv1EncConfiguration*)config_struct)->over_bndry_blk;
    // new nearest comb injection
    scs_ptr->static_config.new_nearest_comb_inject      = ((EbSvtAv1EncConfiguration*)config_struct)->new_nearest_comb_inject;
    // prune unipred at me
    scs_ptr->static_config.prune_unipred_me             = ((EbSvtAv1EncConfiguration*)config_struct)->prune_unipred_me;
    // edge skip angle intra
    scs_ptr->static_config.edge_skp_angle_intra         = ((EbSvtAv1EncConfiguration*)config_struct)->edge_skp_angle_intra;
    // intra angle delta
    scs_ptr->static_config.intra_angle_delta            = ((EbSvtAv1EncConfiguration*)config_struct)->intra_angle_delta;
    // inter intra compoound
    scs_ptr->static_config.inter_intra_compound         = ((EbSvtAv1EncConfiguration*)config_struct)->inter_intra_compound;
    //prune ref frame for rec partitions
    scs_ptr->static_config.prune_ref_rec_part           = ((EbSvtAv1EncConfiguration*)config_struct)->prune_ref_rec_part;
    // NSQ table
    scs_ptr->static_config.nsq_table                    = ((EbSvtAv1EncConfiguration*)config_struct)->nsq_table;
    // frame end cdf update mode
    scs_ptr->static_config.frame_end_cdf_update         = ((EbSvtAv1EncConfiguration*)config_struct)->frame_end_cdf_update;

    // Chroma mode
    scs_ptr->static_config.set_chroma_mode = ((EbSvtAv1EncConfiguration*)config_struct)->set_chroma_mode;

    // Chroma mode
    scs_ptr->static_config.disable_cfl_flag = ((EbSvtAv1EncConfiguration*)config_struct)->disable_cfl_flag;

    // OBMC
    scs_ptr->static_config.enable_obmc = ((EbSvtAv1EncConfiguration*)config_struct)->enable_obmc;

    // RDOQ
    scs_ptr->static_config.enable_rdoq = ((EbSvtAv1EncConfiguration*)config_struct)->enable_rdoq;

    // Predictive ME
    scs_ptr->static_config.pred_me  = ((EbSvtAv1EncConfiguration*)config_struct)->pred_me;
    // BiPred 3x3 injection
    scs_ptr->static_config.bipred_3x3_inject = ((EbSvtAv1EncConfiguration*)config_struct)->bipred_3x3_inject;
    // Compound mode
    scs_ptr->static_config.compound_level = ((EbSvtAv1EncConfiguration*)config_struct)->compound_level;

    scs_ptr->static_config.enable_paeth = ((EbSvtAv1EncConfiguration*)config_struct)->enable_paeth;
    scs_ptr->static_config.enable_smooth = ((EbSvtAv1EncConfiguration*)config_struct)->enable_smooth;

    // Filter intra prediction
    scs_ptr->static_config.enable_filter_intra = ((EbSvtAv1EncConfiguration*)config_struct)->enable_filter_intra;

    // Intra Edge Filter
    scs_ptr->static_config.enable_intra_edge_filter = ((EbSvtAv1EncConfiguration*)config_struct)->enable_intra_edge_filter;

    // Picture based rate estimation
    scs_ptr->static_config.pic_based_rate_est = ((EbSvtAv1EncConfiguration*)config_struct)->pic_based_rate_est;

    // ME Tools
    scs_ptr->static_config.use_default_me_hme = ((EbSvtAv1EncConfiguration*)config_struct)->use_default_me_hme;
    scs_ptr->static_config.enable_hme_flag = ((EbSvtAv1EncConfiguration*)config_struct)->enable_hme_flag;
    scs_ptr->static_config.enable_hme_level0_flag = ((EbSvtAv1EncConfiguration*)config_struct)->enable_hme_level0_flag;
    scs_ptr->static_config.enable_hme_level1_flag = ((EbSvtAv1EncConfiguration*)config_struct)->enable_hme_level1_flag;
    scs_ptr->static_config.enable_hme_level2_flag = ((EbSvtAv1EncConfiguration*)config_struct)->enable_hme_level2_flag;
    scs_ptr->static_config.search_area_width = ((EbSvtAv1EncConfiguration*)config_struct)->search_area_width;
    scs_ptr->static_config.search_area_height = ((EbSvtAv1EncConfiguration*)config_struct)->search_area_height;
    scs_ptr->static_config.number_hme_search_region_in_width = ((EbSvtAv1EncConfiguration*)config_struct)->number_hme_search_region_in_width;
    scs_ptr->static_config.number_hme_search_region_in_height = ((EbSvtAv1EncConfiguration*)config_struct)->number_hme_search_region_in_height;
    scs_ptr->static_config.hme_level0_total_search_area_width = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level0_total_search_area_width;
    scs_ptr->static_config.hme_level0_total_search_area_height = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level0_total_search_area_height;
    scs_ptr->static_config.ext_block_flag = ((EbSvtAv1EncConfiguration*)config_struct)->ext_block_flag;
    for (hme_region_index = 0; hme_region_index < scs_ptr->static_config.number_hme_search_region_in_width; ++hme_region_index) {
        scs_ptr->static_config.hme_level0_search_area_in_width_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level0_search_area_in_width_array[hme_region_index];
        scs_ptr->static_config.hme_level1_search_area_in_width_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level1_search_area_in_width_array[hme_region_index];
        scs_ptr->static_config.hme_level2_search_area_in_width_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level2_search_area_in_width_array[hme_region_index];
    }

    for (hme_region_index = 0; hme_region_index < scs_ptr->static_config.number_hme_search_region_in_height; ++hme_region_index) {
        scs_ptr->static_config.hme_level0_search_area_in_height_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level0_search_area_in_height_array[hme_region_index];
        scs_ptr->static_config.hme_level1_search_area_in_height_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level1_search_area_in_height_array[hme_region_index];
        scs_ptr->static_config.hme_level2_search_area_in_height_array[hme_region_index] = ((EbSvtAv1EncConfiguration*)config_struct)->hme_level2_search_area_in_height_array[hme_region_index];
    }

    //Denoise - Hardcoded(CleanUp)
    scs_ptr->static_config.enable_denoise_flag = 0;

    //Film Grain
    scs_ptr->static_config.film_grain_denoise_strength = ((EbSvtAv1EncConfiguration*)config_struct)->film_grain_denoise_strength;
    scs_ptr->film_grain_denoise_strength = scs_ptr->static_config.film_grain_denoise_strength;

    // MD Parameters
    scs_ptr->static_config.enable_hbd_mode_decision = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth > 8 ? ((EbSvtAv1EncConfiguration*)config_struct)->enable_hbd_mode_decision : 0;
    scs_ptr->static_config.enable_palette = ((EbSvtAv1EncConfiguration*)config_struct)->enable_palette;
    // Adaptive Loop Filter
    scs_ptr->static_config.tile_rows = ((EbSvtAv1EncConfiguration*)config_struct)->tile_rows;
    scs_ptr->static_config.tile_columns = ((EbSvtAv1EncConfiguration*)config_struct)->tile_columns;
    scs_ptr->static_config.unrestricted_motion_vector = ((EbSvtAv1EncConfiguration*)config_struct)->unrestricted_motion_vector;

    // Rate Control
    scs_ptr->static_config.scene_change_detection = ((EbSvtAv1EncConfiguration*)config_struct)->scene_change_detection;
    scs_ptr->static_config.rate_control_mode = ((EbSvtAv1EncConfiguration*)config_struct)->rate_control_mode;
    scs_ptr->static_config.look_ahead_distance = ((EbSvtAv1EncConfiguration*)config_struct)->look_ahead_distance;
    scs_ptr->static_config.frame_rate = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate;
    scs_ptr->static_config.frame_rate_denominator = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate_denominator;
    scs_ptr->static_config.frame_rate_numerator = ((EbSvtAv1EncConfiguration*)config_struct)->frame_rate_numerator;

    scs_ptr->static_config.target_bit_rate = ((EbSvtAv1EncConfiguration*)config_struct)->target_bit_rate;

    scs_ptr->static_config.vbv_bufsize = ((EbSvtAv1EncConfiguration*)config_struct)->vbv_bufsize;

    scs_ptr->static_config.max_qp_allowed = (scs_ptr->static_config.rate_control_mode) ?
        ((EbSvtAv1EncConfiguration*)config_struct)->max_qp_allowed :
        63;

    scs_ptr->static_config.min_qp_allowed = (scs_ptr->static_config.rate_control_mode) ?
        ((EbSvtAv1EncConfiguration*)config_struct)->min_qp_allowed :
        1; // lossless coding not supported

    //Segmentation
    //TODO: check RC mode and set only when RC is enabled in the final version.
    scs_ptr->static_config.enable_adaptive_quantization = config_struct->enable_adaptive_quantization;

    // Misc
    scs_ptr->static_config.encoder_bit_depth = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth;
    scs_ptr->static_config.encoder_color_format = ((EbSvtAv1EncConfiguration*)config_struct)->encoder_color_format;
    if (scs_ptr->static_config.encoder_color_format == EB_YUV400) {
        SVT_LOG("SVT [Warning]: Color format EB_YUV400 not supported, set to EB_YUV420\n");
        scs_ptr->static_config.encoder_color_format = EB_YUV420;
    }
    scs_ptr->chroma_format_idc = (uint32_t)(scs_ptr->static_config.encoder_color_format);
    scs_ptr->encoder_bit_depth = (uint32_t)(scs_ptr->static_config.encoder_bit_depth);
    //16bit pipeline
    scs_ptr->static_config.is_16bit_pipeline = ((((EbSvtAv1EncConfiguration*)config_struct)->encoder_bit_depth) > EB_8BIT) ? EB_TRUE: ((EbSvtAv1EncConfiguration*)config_struct)->is_16bit_pipeline;
    scs_ptr->subsampling_x = (scs_ptr->chroma_format_idc == EB_YUV444 ? 1 : 2) - 1;
    scs_ptr->subsampling_y = (scs_ptr->chroma_format_idc >= EB_YUV422 ? 1 : 2) - 1;
    scs_ptr->static_config.ten_bit_format = ((EbSvtAv1EncConfiguration*)config_struct)->ten_bit_format;
    scs_ptr->static_config.compressed_ten_bit_format = ((EbSvtAv1EncConfiguration*)config_struct)->compressed_ten_bit_format;

    // Thresholds
    scs_ptr->static_config.high_dynamic_range_input = ((EbSvtAv1EncConfiguration*)config_struct)->high_dynamic_range_input;
    scs_ptr->static_config.screen_content_mode = ((EbSvtAv1EncConfiguration*)config_struct)->screen_content_mode;

    scs_ptr->static_config.intrabc_mode = ((EbSvtAv1EncConfiguration*)config_struct)->intrabc_mode;

    // Annex A parameters
    scs_ptr->static_config.profile = ((EbSvtAv1EncConfiguration*)config_struct)->profile;
    scs_ptr->static_config.tier = ((EbSvtAv1EncConfiguration*)config_struct)->tier;
    scs_ptr->static_config.level = ((EbSvtAv1EncConfiguration*)config_struct)->level;
    scs_ptr->static_config.stat_report = ((EbSvtAv1EncConfiguration*)config_struct)->stat_report;

    scs_ptr->static_config.injector_frame_rate = ((EbSvtAv1EncConfiguration*)config_struct)->injector_frame_rate;
    scs_ptr->static_config.speed_control_flag = ((EbSvtAv1EncConfiguration*)config_struct)->speed_control_flag;

    // Buffers - Hardcoded(Cleanup)
    scs_ptr->static_config.use_cpu_flags = ((EbSvtAv1EncConfiguration*)config_struct)->use_cpu_flags;

    scs_ptr->static_config.channel_id = ((EbSvtAv1EncConfiguration*)config_struct)->channel_id;
    scs_ptr->static_config.active_channel_count = ((EbSvtAv1EncConfiguration*)config_struct)->active_channel_count;
    scs_ptr->static_config.logical_processors = ((EbSvtAv1EncConfiguration*)config_struct)->logical_processors;
    scs_ptr->static_config.unpin_lp1 = ((EbSvtAv1EncConfiguration*)config_struct)->unpin_lp1;
    scs_ptr->static_config.target_socket = ((EbSvtAv1EncConfiguration*)config_struct)->target_socket;
    scs_ptr->static_config.qp = ((EbSvtAv1EncConfiguration*)config_struct)->qp;
    scs_ptr->static_config.recon_enabled = ((EbSvtAv1EncConfiguration*)config_struct)->recon_enabled;

    // Extract frame rate from Numerator and Denominator if not 0
    if (scs_ptr->static_config.frame_rate_numerator != 0 && scs_ptr->static_config.frame_rate_denominator != 0)
        scs_ptr->frame_rate = scs_ptr->static_config.frame_rate = (((scs_ptr->static_config.frame_rate_numerator << 8) / (scs_ptr->static_config.frame_rate_denominator)) << 8);
    // Get Default Intra Period if not specified
    if (scs_ptr->static_config.intra_period_length == -2)
        scs_ptr->intra_period_length = scs_ptr->static_config.intra_period_length = compute_default_intra_period(scs_ptr);
    if (scs_ptr->static_config.look_ahead_distance == (uint32_t)~0)
        scs_ptr->static_config.look_ahead_distance = compute_default_look_ahead(&scs_ptr->static_config);
    else
        scs_ptr->static_config.look_ahead_distance = cap_look_ahead_distance(&scs_ptr->static_config);

    scs_ptr->static_config.enable_altrefs = config_struct->enable_altrefs;
    scs_ptr->static_config.altref_strength = config_struct->altref_strength;
    scs_ptr->static_config.altref_nframes = config_struct->altref_nframes;
    scs_ptr->static_config.enable_overlays = config_struct->enable_overlays;

    scs_ptr->static_config.superres_mode = config_struct->superres_mode;
    scs_ptr->static_config.superres_denom = config_struct->superres_denom;
    scs_ptr->static_config.superres_kf_denom = config_struct->superres_kf_denom;
    scs_ptr->static_config.superres_qthres = config_struct->superres_qthres;

    scs_ptr->static_config.sq_weight = config_struct->sq_weight;

    scs_ptr->static_config.md_stage_1_cand_prune_th = config_struct->md_stage_1_cand_prune_th;
    scs_ptr->static_config.md_stage_1_class_prune_th = config_struct->md_stage_1_class_prune_th;
    scs_ptr->static_config.md_stage_2_3_cand_prune_th = config_struct->md_stage_2_3_cand_prune_th;
    scs_ptr->static_config.md_stage_2_3_class_prune_th = config_struct->md_stage_2_3_class_prune_th;

    // Prediction Structure
    scs_ptr->static_config.enable_manual_pred_struct    = config_struct->enable_manual_pred_struct;
    if(scs_ptr->static_config.enable_manual_pred_struct){
        scs_ptr->static_config.manual_pred_struct_entry_num = config_struct->manual_pred_struct_entry_num;
        memcpy(&scs_ptr->static_config.pred_struct[0], &config_struct->pred_struct[0],config_struct->manual_pred_struct_entry_num*sizeof(PredictionStructureConfigEntry));
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

    return;
}

/******************************************
* Verify Settings
******************************************/
#define PowerOfTwoCheck(x) (((x) != 0) && (((x) & (~(x) + 1)) == (x)))

static int verify_hme_dimension(unsigned int index, unsigned int HmeLevel0SearchAreaInWidth, uint32_t number_hme_search_region_in_width_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT], unsigned int number_hme_search_region_in_width)
{
    int           return_error = 0;
    uint32_t        i;
    uint32_t        total_search_width = 0;

    for (i = 0; i < number_hme_search_region_in_width; i++)
        total_search_width += number_hme_search_region_in_width_array[i];
    if ((total_search_width) != (HmeLevel0SearchAreaInWidth)) {
        SVT_LOG("Error Instance %u: Invalid  HME Total Search Area. \n", index);
        return_error = -1;
        return return_error;
    }

    return return_error;
}
static int verify_hme_dimension_l1_l2(unsigned int index, uint32_t number_hme_search_region_in_width_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT], unsigned int number_hme_search_region_in_width)
{
    int             return_error = 0;
    uint32_t        i;
    uint32_t        total_search_width = 0;

    for (i = 0; i < number_hme_search_region_in_width; i++)
        total_search_width += number_hme_search_region_in_width_array[i];
    if ((total_search_width > 256) || (total_search_width == 0)) {
        SVT_LOG("Error Instance %u: Invalid  HME Total Search Area. Must be [1 - 256].\n", index);
        return_error = -1;
        return return_error;
    }

    return return_error;
}

static EbErrorType verify_settings(
    SequenceControlSet       *scs_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbSvtAv1EncConfiguration *config = &scs_ptr->static_config;
    unsigned int channel_number = config->channel_id;
    if (config->enc_mode > MAX_ENC_PRESET) {
        SVT_LOG("Error instance %u: EncoderMode must be in the range of [0-%d]\n", channel_number + 1, MAX_ENC_PRESET);
        return_error = EB_ErrorBadParameter;
    }
    if (config->snd_pass_enc_mode > MAX_ENC_PRESET + 1) {
        SVT_LOG("Error instance %u: Second pass encoder mode must be in the range of [0-%d]\n", channel_number + 1, MAX_ENC_PRESET + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->ext_block_flag > 1) {
        SVT_LOG("Error instance %u: ExtBlockFlag must be [0-1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (scs_ptr->max_input_luma_width < 64) {
        SVT_LOG("Error instance %u: Source Width must be at least 64\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (scs_ptr->max_input_luma_height < 64) {
        SVT_LOG("Error instance %u: Source Width must be at least 64\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->pred_structure != 2) {
        SVT_LOG("Error instance %u: Pred Structure must be [2]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (scs_ptr->max_input_luma_width % 8 && scs_ptr->static_config.compressed_ten_bit_format == 1) {
        SVT_LOG("Error Instance %u: Only multiple of 8 width is supported for compressed 10-bit inputs \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (scs_ptr->max_input_luma_width % 2) {
        SVT_LOG("Error Instance %u: Source Width must be even for YUV_420 colorspace\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (scs_ptr->max_input_luma_height % 2) {
        SVT_LOG("Error Instance %u: Source Height must be even for YUV_420 colorspace\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (scs_ptr->max_input_luma_width > 4096) {
        SVT_LOG("Error instance %u: Source Width must be less than 4096\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (scs_ptr->max_input_luma_height > 2160) {
        SVT_LOG("Error instance %u: Source Height must be less than 2160\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->qp > MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: QP must be [0 - %d]\n", channel_number + 1, MAX_QP_VALUE);
        return_error = EB_ErrorBadParameter;
    }
    if (config->hierarchical_levels != 3 && config->hierarchical_levels != 4 && config->hierarchical_levels != 5) {
        SVT_LOG("Error instance %u: Hierarchical Levels supported [3-5]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->intra_period_length < -2 || config->intra_period_length > 255) {
        SVT_LOG("Error Instance %u: The intra period must be [-2 - 255] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->intra_refresh_type > 2 || config->intra_refresh_type < 1) {
        SVT_LOG("Error Instance %u: Invalid intra Refresh Type [1-2]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->disable_dlf_flag > 1) {
        SVT_LOG("Error Instance %u: Invalid LoopFilterDisable. LoopFilterDisable must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_default_me_hme > 1) {
        SVT_LOG("Error Instance %u: invalid use_default_me_hme. use_default_me_hme must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->enable_hme_flag > 1) {
        SVT_LOG("Error Instance %u: invalid HME. HME must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level0_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel0. HMELevel0 must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level1_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel1. HMELevel1 must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level2_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel2. HMELevel2 must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->search_area_width > 256) || (config->search_area_width == 0)) {
        SVT_LOG("Error Instance %u: Invalid search_area_width. search_area_width must be [1 - 256]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->search_area_height > 256) || (config->search_area_height == 0)) {
        SVT_LOG("Error Instance %u: Invalid search_area_height. search_area_height must be [1 - 256]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_flag) {
        if ((config->number_hme_search_region_in_width > (uint32_t)EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT) || (config->number_hme_search_region_in_width == 0)) {
            SVT_LOG("Error Instance %u: Invalid number_hme_search_region_in_width. number_hme_search_region_in_width must be [1 - %d]\n", channel_number + 1, EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT);
            return_error = EB_ErrorBadParameter;
        }

        if ((config->number_hme_search_region_in_height > (uint32_t)EB_HME_SEARCH_AREA_ROW_MAX_COUNT) || (config->number_hme_search_region_in_height == 0)) {
            SVT_LOG("Error Instance %u: Invalid number_hme_search_region_in_height. number_hme_search_region_in_height must be [1 - %d]\n", channel_number + 1, EB_HME_SEARCH_AREA_ROW_MAX_COUNT);
            return_error = EB_ErrorBadParameter;
        }

        if ((config->hme_level0_total_search_area_height > 256) || (config->hme_level0_total_search_area_height == 0)) {
            SVT_LOG("Error Instance %u: Invalid hme_level0_total_search_area_height. hme_level0_total_search_area_height must be [1 - 256]\n", channel_number + 1);
            return_error = EB_ErrorBadParameter;
        }
        if ((config->hme_level0_total_search_area_width > 256) || (config->hme_level0_total_search_area_width == 0)) {
            SVT_LOG("Error Instance %u: Invalid hme_level0_total_search_area_width. hme_level0_total_search_area_width must be [1 - 256]\n", channel_number + 1);
            return_error = EB_ErrorBadParameter;
        }
        if (verify_hme_dimension(channel_number + 1, config->hme_level0_total_search_area_height, config->hme_level0_search_area_in_height_array, config->number_hme_search_region_in_height))
            return_error = EB_ErrorBadParameter;
        if (verify_hme_dimension(channel_number + 1, config->hme_level0_total_search_area_width, config->hme_level0_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (verify_hme_dimension_l1_l2(channel_number + 1, config->hme_level1_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (verify_hme_dimension_l1_l2(channel_number + 1, config->hme_level1_search_area_in_height_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (verify_hme_dimension_l1_l2(channel_number + 1, config->hme_level2_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (verify_hme_dimension_l1_l2(channel_number + 1, config->hme_level2_search_area_in_height_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
    }

    if (config->profile > 2) {
        SVT_LOG("Error Instance %u: The maximum allowed profile value is 2 \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    // Check if the current input video is conformant with the Level constraint
    if (config->frame_rate > (240 << 16)) {
        SVT_LOG("Error Instance %u: The maximum allowed frame rate is 240 fps\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    // Check that the frame_rate is non-zero
    if (config->frame_rate <= 0) {
        SVT_LOG("Error Instance %u: The frame rate should be greater than 0 fps \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->intra_period_length < -2 || config->intra_period_length > 255) {
        SVT_LOG("Error Instance %u: The intra period must be [-2 - 255] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->rate_control_mode > 2) {

        SVT_LOG("Error Instance %u: The rate control mode must be [0 - 2] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if ((config->rate_control_mode == 3|| config->rate_control_mode == 2) && config->look_ahead_distance != (uint32_t)config->intra_period_length && config->intra_period_length >= 0) {
        SVT_LOG("Error Instance %u: The rate control mode 2/3 LAD must be equal to intra_period \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->look_ahead_distance > MAX_LAD && config->look_ahead_distance != (uint32_t)~0) {
        SVT_LOG("Error Instance %u: The lookahead distance must be [0 - %d] \n", channel_number + 1, MAX_LAD);

        return_error = EB_ErrorBadParameter;
    }
    if (config->tile_rows < 0 || config->tile_columns < 0 || config->tile_rows > 6 || config->tile_columns > 6) {
        SVT_LOG("Error Instance %u: Log2Tile rows/cols must be [0 - 6] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
#if TILES_PARALLEL
    if ((1 << config->tile_rows) * (1 << config->tile_columns) > 128 || config->tile_columns > 4) {
        SVT_LOG("Error Instance %u: MaxTiles is 128 and MaxTileCols is 16 (Annex A.3) \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
#endif
    if (config->unrestricted_motion_vector > 1) {
        SVT_LOG("Error Instance %u : Invalid Unrestricted Motion Vector flag [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->scene_change_detection > 1) {
        SVT_LOG("Error Instance %u: The scene change detection must be [0 - 1] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->max_qp_allowed > MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: MaxQpAllowed must be [0 - %d]\n", channel_number + 1, MAX_QP_VALUE);
        return_error = EB_ErrorBadParameter;
    }
    else if (config->min_qp_allowed >= MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: MinQpAllowed must be [0 - %d]\n", channel_number + 1, MAX_QP_VALUE-1);
        return_error = EB_ErrorBadParameter;
    }
    else if ((config->min_qp_allowed) > (config->max_qp_allowed)) {
        SVT_LOG("Error Instance %u:  MinQpAllowed must be smaller than MaxQpAllowed\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->stat_report > 1) {
        SVT_LOG("Error instance %u : Invalid StatReport. StatReport must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->high_dynamic_range_input > 1) {
        SVT_LOG("Error instance %u : Invalid HighDynamicRangeInput. HighDynamicRangeInput must be [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->screen_content_mode > 2) {
        SVT_LOG("Error instance %u : Invalid screen_content_mode. screen_content_mode must be [0 - 2]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    // IntraBC
    if (config->intrabc_mode > 3 || config->intrabc_mode < -1) {
        SVT_LOG( "Error instance %u: Invalid intraBC mode [0-3, -1 for default], your input: %i\n", channel_number + 1, config->intrabc_mode);
        return_error = EB_ErrorBadParameter;
    }

    if (scs_ptr->static_config.enable_adaptive_quantization > 2) {
        SVT_LOG("Error instance %u : Invalid enable_adaptive_quantization. enable_adaptive_quantization must be [0-2]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->encoder_bit_depth != 8) &&
        (config->encoder_bit_depth != 10)
        ) {
        SVT_LOG("Error instance %u: Encoder Bit Depth shall be only 8 or 10 \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    // Check if the EncoderBitDepth is conformant with the Profile constraint
    if ((config->profile == 0 || config->profile == 1) && config->encoder_bit_depth > 10) {
        SVT_LOG("Error instance %u: The encoder bit depth shall be equal to 8 or 10 for Main/High Profile\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->encoder_color_format != EB_YUV420) {
        SVT_LOG("Error instance %u: Only support 420 now \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 0 && config->encoder_color_format > EB_YUV420) {
        SVT_LOG("Error instance %u: Non 420 color format requires profile 1 or 2\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 1 && config->encoder_color_format != EB_YUV444) {
        SVT_LOG("Error instance %u: Profile 1 requires 4:4:4 color format\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 2 && config->encoder_bit_depth <= 10 && config->encoder_color_format != EB_YUV422) {
        SVT_LOG("Error instance %u: Profile 2 bit-depth < 10 requires 4:2:2 color format\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->compressed_ten_bit_format !=0)
    {
        SVT_LOG("Error instance %u: Compressed ten bit format is not supported in this version \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->speed_control_flag > 1) {
        SVT_LOG("Error Instance %u: Invalid Speed Control flag [0 - 1]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_cpu_flags & CPU_FLAGS_INVALID) {
        SVT_LOG("Error Instance %u: param '-asm' have invalid value.\n"
            "Value should be [0 - 11, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512, max]\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->target_socket != -1 && config->target_socket != 0 && config->target_socket != 1) {
        SVT_LOG("Error instance %u: Invalid target_socket. target_socket must be [-1 - 1] \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    // alt-ref frames related
    if (config->altref_strength > ALTREF_MAX_STRENGTH ) {
        SVT_LOG("Error instance %u: invalid altref-strength, should be in the range [0 - %d] \n", channel_number + 1, ALTREF_MAX_STRENGTH);
        return_error = EB_ErrorBadParameter;
    }

    if (config->altref_nframes > ALTREF_MAX_NFRAMES ) {
        SVT_LOG("Error instance %u: invalid altref-nframes, should be in the range [0 - %d] \n", channel_number + 1, ALTREF_MAX_NFRAMES);
        return_error = EB_ErrorBadParameter;
    }

    // Local Warped Motion
    if (config->enable_warped_motion != 0 && config->enable_warped_motion != 1 && config->enable_warped_motion != -1) {
      SVT_LOG("Error instance %u: Invalid warped motion flag [0/1, -1], your input: %d\n", channel_number + 1, config->enable_warped_motion);
      return_error = EB_ErrorBadParameter;
    }

    // Global Motion
    if (config->enable_global_motion != 0 && config->enable_global_motion != 1) {
      SVT_LOG("Error instance %u: Invalid global motion flag [0 - 1], your input: %d\n", channel_number + 1, config->enable_global_motion);
      return_error = EB_ErrorBadParameter;
    }

    // OBMC
    if (config->enable_obmc != 0 && config->enable_obmc != 1) {
      SVT_LOG("Error instance %u: Invalid OBMC flag [0 - 1], your input: %d\n", channel_number + 1, config->enable_obmc);
      return_error = EB_ErrorBadParameter;
    }

    // Filter Intra prediction
    if (config->enable_filter_intra != 0 && config->enable_filter_intra != 1) {
      SVT_LOG("Error instance %u: Invalid Filter Intra flag [0 - 1], your input: %d\n", channel_number + 1, config->enable_filter_intra);
      return_error = EB_ErrorBadParameter;
    }

    // Intra Edge Filter
    if (config->enable_intra_edge_filter != 0 && config->enable_intra_edge_filter != 1 && config->enable_intra_edge_filter != -1) {
        SVT_LOG("Error instance %u: Invalid Filter Intra flag [0/1, -1], your input: %d\n", channel_number + 1, config->enable_intra_edge_filter);
        return_error = EB_ErrorBadParameter;
    }

    // Picture based rate estimation
    if (config->pic_based_rate_est != 0 && config->pic_based_rate_est != 1 && config->pic_based_rate_est != -1) {
        SVT_LOG("Error instance %u: Invalid pic_based_rate_est [0/1, -1], your input: %d\n", channel_number + 1, config->pic_based_rate_est);
        return_error = EB_ErrorBadParameter;
    }

    // HBD mode decision
    if (config->enable_hbd_mode_decision != 0 && config->enable_hbd_mode_decision != 1 && config->enable_hbd_mode_decision != 2) {
    SVT_LOG("Error instance %u: Invalid HBD mode decision flag [0 - 2], your input: %d\n", channel_number + 1, config->enable_hbd_mode_decision);
    return_error = EB_ErrorBadParameter;
    }

    // palette
    if (config->enable_palette < (int32_t)(-1) || config->enable_palette >6) {
        SVT_LOG( "Error instance %u: Invalid Palette Mode [0 .. 6], your input: %i\n", channel_number + 1, config->enable_palette);
        return_error = EB_ErrorBadParameter;
    }

    // RDOQ
    if (config->enable_rdoq != 0 && config->enable_rdoq != 1 && config->enable_rdoq != -1) {
        SVT_LOG( "Error instance %u: Invalid RDOQ parameter [-1, 0, 1], your input: %i\n", channel_number + 1, config->enable_rdoq);
        return_error = EB_ErrorBadParameter;
    }

    // Chroma Level
    if (config->set_chroma_mode > 3 || config->set_chroma_mode < -1) {
      SVT_LOG("Error instance %u: Invalid Chroma Mode [0 - 3, -1 for auto], your input: %d\n", channel_number + 1, config->set_chroma_mode);
      return_error = EB_ErrorBadParameter;
    }

    // Disable chroma from luma (CFL)
    if (config->disable_cfl_flag != 0 && config->disable_cfl_flag != 1 && config->disable_cfl_flag != -1) {
        SVT_LOG( "Error instance %u: Invalid CFL flag [0/1, -1], your input: %i\n", channel_number + 1, config->disable_cfl_flag);
        return_error = EB_ErrorBadParameter;
    }

    // CDEF
    if (config->cdef_mode > 5 || config->cdef_mode < -1) {
        SVT_LOG("Error instance %u: Invalid CDEF mode [0 - 5, -1 for auto], your input: %d\n", channel_number + 1, config->cdef_mode);
        return_error = EB_ErrorBadParameter;
    }

    // Restoration Filtering
    if (config->enable_restoration_filtering != 0 && config->enable_restoration_filtering != 1 && config->enable_restoration_filtering != -1) {
      SVT_LOG("Error instance %u: Invalid restoration flag [0 - 1, -1 for auto], your input: %d\n", channel_number + 1, config->enable_restoration_filtering);
      return_error = EB_ErrorBadParameter;
    }

    if (config->sg_filter_mode > 4 || config->sg_filter_mode < -1) {
        SVT_LOG("Error instance %u: Invalid self-guided filter mode [0 - 4, -1 for auto], your input: %d\n", channel_number + 1, config->sg_filter_mode);
        return_error = EB_ErrorBadParameter;
    }

    if (config->wn_filter_mode > 3 || config->wn_filter_mode < -1) {
        SVT_LOG("Error instance %u: Invalid Wiener filter mode [0 - 3, -1 for auto], your input: %d\n", channel_number + 1, config->wn_filter_mode);
        return_error = EB_ErrorBadParameter;
    }

    if (config->pred_me > 5 || config->pred_me < -1) {
      SVT_LOG("Error instance %u: Invalid predictive me level [0-5, -1 for auto], your input: %d\n", channel_number + 1, config->pred_me);
      return_error = EB_ErrorBadParameter;
    }

    if (config->bipred_3x3_inject > 2 || config->bipred_3x3_inject < -1) {
      SVT_LOG("Error instance %u: Invalid bipred_3x3_inject mode [0-2, -1 for auto], your input: %d\n", channel_number + 1, config->bipred_3x3_inject);
      return_error = EB_ErrorBadParameter;
    }

    if (config->compound_level > 2 || config->compound_level < -1) {
      SVT_LOG("Error instance %u: Invalid compound level [0-2, -1 for auto], your input: %d\n", channel_number + 1, config->compound_level);
      return_error = EB_ErrorBadParameter;
    }

    if (config->combine_class_12 != 0 && config->combine_class_12 != 1 && config->combine_class_12 != -1) {
      SVT_LOG("Error instance %u: Invalid combine MD Class1&2 flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->combine_class_12);
      return_error = EB_ErrorBadParameter;
    }

    if (config->edge_skp_angle_intra != 0 && config->edge_skp_angle_intra != 1 && config->edge_skp_angle_intra != -1) {
      SVT_LOG("Error instance %u: Invalid Enable skip angle intra based on edge flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->edge_skp_angle_intra);
      return_error = EB_ErrorBadParameter;
    }

    if (config->intra_angle_delta != 0 && config->intra_angle_delta != 1 && config->intra_angle_delta != -1) {
        SVT_LOG("Error instance %u: Invalid Enable intra angle delta flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->intra_angle_delta);
        return_error = EB_ErrorBadParameter;
    }

    if (config->inter_intra_compound != 0 && config->inter_intra_compound != 1 && config->inter_intra_compound != -1) {
      SVT_LOG("Error instance %u: Invalid Inter Intra Compound flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->inter_intra_compound);
      return_error = EB_ErrorBadParameter;
    }

    if (config->enable_paeth != 0 && config->enable_paeth != 1 && config->enable_paeth != -1) {
        SVT_LOG("Error instance %u: Invalid Paeth flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->enable_paeth);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_smooth != 0 && config->enable_smooth != 1 && config->enable_smooth != -1) {
        SVT_LOG("Error instance %u: Invalid Smooth flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->enable_smooth);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_mfmv != 0 && config->enable_mfmv != 1 && config->enable_mfmv != -1) {
      SVT_LOG("Error instance %u: Invalid motion field motion vector flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->enable_mfmv);
      return_error = EB_ErrorBadParameter;
    }

    if (config->enable_redundant_blk != 0 && config->enable_redundant_blk != 1 && config->enable_redundant_blk != -1) {
      SVT_LOG("Error instance %u: Invalid enable_redundant_blk  flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->enable_redundant_blk);
      return_error = EB_ErrorBadParameter;
    }

    if (config->spatial_sse_fl != 0 && config->spatial_sse_fl != 1 && config->spatial_sse_fl != -1) {
      SVT_LOG("Error instance %u: Invalid spatial_sse_fl flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->spatial_sse_fl);
      return_error = EB_ErrorBadParameter;
    }

    if (config->enable_subpel != 0 && config->enable_subpel != 1 && config->enable_subpel != -1) {
      SVT_LOG("Error instance %u: Invalid enable_subpel flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->enable_subpel);
      return_error = EB_ErrorBadParameter;
    }

    if (config->over_bndry_blk != 0 && config->over_bndry_blk != 1 && config->over_bndry_blk != -1) {
      SVT_LOG("Error instance %u: Invalid over_bndry_blk flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->over_bndry_blk);
      return_error = EB_ErrorBadParameter;
    }

    if (config->new_nearest_comb_inject != 0 && config->new_nearest_comb_inject != 1 && config->new_nearest_comb_inject != -1) {
      SVT_LOG("Error instance %u: Invalid new_nearest_comb_inject flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->new_nearest_comb_inject);
      return_error = EB_ErrorBadParameter;
    }

    if (config->prune_unipred_me != 0 && config->prune_unipred_me != 1 && config->prune_unipred_me != -1) {
      SVT_LOG("Error instance %u: Invalid prune_unipred_me flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->prune_unipred_me);
      return_error = EB_ErrorBadParameter;
    }

    if (config->prune_ref_rec_part != 0 && config->prune_ref_rec_part != 1 && config->prune_ref_rec_part != -1) {
      SVT_LOG("Error instance %u: Invalid prune_ref_rec_part flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->prune_ref_rec_part);
      return_error = EB_ErrorBadParameter;
    }

    if (config->nsq_table != 0 && config->nsq_table != 1 && config->nsq_table != -1) {
      SVT_LOG("Error instance %u: Invalid nsq_table flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->nsq_table);
      return_error = EB_ErrorBadParameter;
    }

    if (config->frame_end_cdf_update != 0 && config->frame_end_cdf_update != 1 && config->frame_end_cdf_update != -1) {
      SVT_LOG("Error instance %u: Invalid frame_end_cdf_update flag [0/1 or -1 for auto], your input: %d\n", channel_number + 1, config->frame_end_cdf_update);
      return_error = EB_ErrorBadParameter;
    }

    // prediction structure
    if(config->enable_manual_pred_struct) {
        if(config->manual_pred_struct_entry_num > (1<<(MAX_HIERARCHICAL_LEVEL-1))){
            SVT_LOG("Error instance %u: Invalid manual prediction structure entry number [1 - 32], your input: %d\n", channel_number + 1, config->manual_pred_struct_entry_num);
            return_error = EB_ErrorBadParameter;
        }
        else {
            for(int32_t i = 0; i < config->manual_pred_struct_entry_num; i++) {
                config->pred_struct[i].ref_list1[REF_LIST_MAX_DEPTH-1] = 0;
                if(config->pred_struct[i].decode_order >= (1<<(MAX_HIERARCHICAL_LEVEL-1))){
                    SVT_LOG("Error instance %u: Invalid decode order for manual prediction structure [0 - 31], your input: %d\n", channel_number + 1, config->pred_struct[i].decode_order);
                    return_error = EB_ErrorBadParameter;
                }
                if(config->pred_struct[i].temporal_layer_index >= (1<<(MAX_HIERARCHICAL_LEVEL-1))){
                    SVT_LOG("Error instance %u: Invalid temporal layer index for manual prediction structure [0 - 31], your input: %d\n", channel_number + 1, config->pred_struct[i].temporal_layer_index);
                    return_error = EB_ErrorBadParameter;
                }
                EbBool have_ref_frame_within_minigop_in_list0 = EB_FALSE;
                int32_t entry_idx = i + 1;
                for(int32_t j = 0; j < REF_LIST_MAX_DEPTH; j++) {
                    if((entry_idx - config->pred_struct[i].ref_list1[j] > config->manual_pred_struct_entry_num)) {
                        SVT_LOG("Error instance %u: Invalid ref frame %d in list1 entry%d for manual prediction structure, all ref frames in list1 should not exceed minigop end\n",
                        channel_number + 1, config->pred_struct[i].ref_list1[j], i);
                        return_error = EB_ErrorBadParameter;
                    }
                    if(config->pred_struct[i].ref_list0[j] < 0) {
                        SVT_LOG("Error instance %u: Invalid ref frame %d in list0 entry%d for manual prediction structure, only forward frames can be in list0\n",
                        channel_number + 1, config->pred_struct[i].ref_list0[j], i);
                        return_error = EB_ErrorBadParameter;
                    }
                    if(!have_ref_frame_within_minigop_in_list0 && config->pred_struct[i].ref_list0[j] && entry_idx - config->pred_struct[i].ref_list0[j] >= 0 ) {
                        have_ref_frame_within_minigop_in_list0 = EB_TRUE;
                    }
                }
                if(!have_ref_frame_within_minigop_in_list0) {
                    SVT_LOG("Error instance %u: Invalid ref frame in list0 entry%d for manual prediction structure,there should be at least one frame within minigop \n",
                    channel_number + 1, i);
                    return_error = EB_ErrorBadParameter;
                }
            }
        }
    }

    if (config->superres_mode > 2) {
        SVT_LOG("Error instance %u: invalid superres-mode %d, should be in the range [%d - %d], "
                "only SUPERRES_NONE (0), SUPERRES_FIXED (1) and SUPERRES_RANDOM (2) are currently implemented \n", channel_number + 1, config->superres_mode, 0, 2);
        return_error = EB_ErrorBadParameter;
    }

    if (config->superres_mode > 0 && (config->input_stat_file || config->output_stat_file)){
        SVT_LOG("Error instance %u: superres cannot be enabled in 2-pass mode yet \n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->superres_qthres > MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: invalid superres-qthres %d, should be in the range [%d - %d] \n", channel_number + 1, config->superres_qthres, MIN_QP_VALUE, MAX_QP_VALUE);
        return_error = EB_ErrorBadParameter;
    }

    if (config->superres_kf_denom < MIN_SUPERRES_DENOM || config->superres_kf_denom > MAX_SUPERRES_DENOM) {
        SVT_LOG("Error instance %u: invalid superres-kf-denom %d, should be in the range [%d - %d] \n", channel_number + 1, config->superres_kf_denom, MIN_SUPERRES_DENOM, MAX_SUPERRES_DENOM);
        return_error = EB_ErrorBadParameter;
    }

    if (config->superres_denom < MIN_SUPERRES_DENOM || config->superres_denom > MAX_SUPERRES_DENOM) {
        SVT_LOG("Error instance %u: invalid superres-denom %d, should be in the range [%d - %d] \n", channel_number + 1, config->superres_denom, MIN_SUPERRES_DENOM, MAX_SUPERRES_DENOM);
        return_error = EB_ErrorBadParameter;
    }
    return return_error;
}

/**********************************
Set Default Library Params
**********************************/
EbErrorType eb_svt_enc_init_parameter(
    EbSvtAv1EncConfiguration * config_ptr)
{
    EbErrorType                  return_error = EB_ErrorNone;

    if (!config_ptr) {
        SVT_LOG("The EbSvtAv1EncConfiguration structure is empty! \n");
        return EB_ErrorBadParameter;
    }

    config_ptr->frame_rate = 30 << 16;
    config_ptr->frame_rate_numerator = 0;
    config_ptr->frame_rate_denominator = 0;
    config_ptr->encoder_bit_depth = 8;
    config_ptr->ten_bit_format = 0;
    config_ptr->compressed_ten_bit_format = 0;
    config_ptr->source_width = 0;
    config_ptr->source_height = 0;
    config_ptr->stat_report = 0;
    config_ptr->tile_rows = 0;
    config_ptr->tile_columns = 0;

    config_ptr->qp = 50;
    config_ptr->use_qp_file = EB_FALSE;
    config_ptr->scene_change_detection = 0;
    config_ptr->rate_control_mode = 0;
    config_ptr->look_ahead_distance = (uint32_t)~0;
    config_ptr->target_bit_rate = 7000000;
    config_ptr->max_qp_allowed = 63;
    config_ptr->min_qp_allowed = 10;
    config_ptr->enc_mode = MAX_ENC_PRESET;
    config_ptr->snd_pass_enc_mode = MAX_ENC_PRESET + 1;
    config_ptr->intra_period_length = -2;
    config_ptr->intra_refresh_type = 1;
    config_ptr->hierarchical_levels = 4;
    config_ptr->pred_structure = EB_PRED_RANDOM_ACCESS;
    config_ptr->disable_dlf_flag = EB_FALSE;
    config_ptr->enable_warped_motion = DEFAULT;
    config_ptr->enable_global_motion = EB_TRUE;
    config_ptr->cdef_mode = DEFAULT;
    config_ptr->enable_restoration_filtering = DEFAULT;
    config_ptr->sg_filter_mode = DEFAULT;
    config_ptr->wn_filter_mode = DEFAULT;
    config_ptr->edge_skp_angle_intra = DEFAULT;
    config_ptr->intra_angle_delta = DEFAULT;
    config_ptr->combine_class_12 = DEFAULT;
    config_ptr->inter_intra_compound = DEFAULT;
    config_ptr->enable_paeth = DEFAULT;
    config_ptr->enable_smooth = DEFAULT;
    config_ptr->enable_mfmv = DEFAULT;
    config_ptr->enable_redundant_blk = DEFAULT;
    config_ptr->spatial_sse_fl = DEFAULT;
    config_ptr->enable_subpel = DEFAULT;
    config_ptr->over_bndry_blk = DEFAULT;
    config_ptr->new_nearest_comb_inject = DEFAULT;
    config_ptr->prune_unipred_me = DEFAULT;
    config_ptr->prune_ref_rec_part = DEFAULT;
    config_ptr->nsq_table = DEFAULT;
    config_ptr->frame_end_cdf_update = DEFAULT;
    config_ptr->set_chroma_mode = DEFAULT;
    config_ptr->disable_cfl_flag = DEFAULT;
    config_ptr->enable_obmc = EB_TRUE;
    config_ptr->enable_rdoq = DEFAULT;
    config_ptr->pred_me = DEFAULT;
    config_ptr->bipred_3x3_inject = DEFAULT;
    config_ptr->compound_level = DEFAULT;
    config_ptr->enable_filter_intra = EB_TRUE;
    config_ptr->enable_intra_edge_filter = DEFAULT;
    config_ptr->pic_based_rate_est = DEFAULT;
    config_ptr->ext_block_flag = EB_FALSE;
    config_ptr->use_default_me_hme = EB_TRUE;
    config_ptr->enable_hme_flag = EB_TRUE;
    config_ptr->enable_hme_level0_flag = EB_TRUE;
    config_ptr->enable_hme_level1_flag = EB_FALSE;
    config_ptr->enable_hme_level2_flag = EB_FALSE;
    config_ptr->search_area_width = 16;
    config_ptr->search_area_height = 7;
    config_ptr->number_hme_search_region_in_width = 2;
    config_ptr->number_hme_search_region_in_height = 2;
    config_ptr->hme_level0_total_search_area_width = 64;
    config_ptr->hme_level0_total_search_area_height = 25;
    config_ptr->hme_level0_search_area_in_width_array[0] = 32;
    config_ptr->hme_level0_search_area_in_width_array[1] = 32;
    config_ptr->hme_level0_search_area_in_height_array[0] = 12;
    config_ptr->hme_level0_search_area_in_height_array[1] = 13;
    config_ptr->hme_level1_search_area_in_width_array[0] = 1;
    config_ptr->hme_level1_search_area_in_width_array[1] = 1;
    config_ptr->hme_level1_search_area_in_height_array[0] = 1;
    config_ptr->hme_level1_search_area_in_height_array[1] = 1;
    config_ptr->hme_level2_search_area_in_width_array[0] = 1;
    config_ptr->hme_level2_search_area_in_width_array[1] = 1;
    config_ptr->hme_level2_search_area_in_height_array[0] = 1;
    config_ptr->hme_level2_search_area_in_height_array[1] = 1;
    config_ptr->enable_hbd_mode_decision = 1;
    config_ptr->enable_palette = -1;
    // Bitstream options
    //config_ptr->codeVpsSpsPps = 0;
    //config_ptr->codeEosNal = 0;
    config_ptr->unrestricted_motion_vector = EB_TRUE;

    config_ptr->high_dynamic_range_input = 0;
    config_ptr->screen_content_mode = 0;

    config_ptr->intrabc_mode = DEFAULT;

    // Annex A parameters
    config_ptr->profile = 0;
    config_ptr->tier = 0;
    config_ptr->level = 0;

    // Latency
    config_ptr->injector_frame_rate = 60 << 16;
    config_ptr->speed_control_flag = 0;
    config_ptr->super_block_size = 128;

    config_ptr->sb_sz = 64;
    config_ptr->partition_depth = (uint8_t)EB_MAX_SB_DEPTH;
    //config_ptr->latency_mode = 0;
    config_ptr->speed_control_flag = 0;
    config_ptr->film_grain_denoise_strength = 0;

    // CPU Flags
    config_ptr->use_cpu_flags = CPU_FLAGS_ALL;

    // Channel info
    config_ptr->logical_processors = 0;
    config_ptr->unpin_lp1 = 1;
    config_ptr->target_socket = -1;
    config_ptr->channel_id = 0;
    config_ptr->active_channel_count = 1;

    // Debug info
    config_ptr->recon_enabled = 0;

    // Alt-Ref default values
    config_ptr->enable_altrefs = EB_TRUE;
    config_ptr->altref_nframes = 7;
    config_ptr->altref_strength = 5;
    config_ptr->enable_overlays = EB_FALSE;

    // Super-resolution default values
    config_ptr->superres_mode = SUPERRES_NONE;
    config_ptr->superres_denom = 8;
    config_ptr->superres_kf_denom = 8;
    config_ptr->superres_qthres = 43; // random threshold, change

    config_ptr->sq_weight = 100;

    config_ptr->md_stage_1_cand_prune_th = 75;
    config_ptr->md_stage_1_class_prune_th = 100;
    config_ptr->md_stage_2_3_cand_prune_th = 15;
    config_ptr->md_stage_2_3_class_prune_th = 25;

    return return_error;
}
//#define DEBUG_BUFFERS
static void print_lib_params(
    SequenceControlSet* scs) {
    EbSvtAv1EncConfiguration*   config = &scs->static_config;

    SVT_LOG("------------------------------------------- ");
    if (config->profile == MAIN_PROFILE)
        SVT_LOG("\nSVT [config]: Main Profile\t");
    else if (config->profile == HIGH_PROFILE)
        SVT_LOG("\nSVT [config]: High Profile\t");
    else if (config->profile == PROFESSIONAL_PROFILE)
        SVT_LOG("\nSVT [config]: Professional Profile\t");
    else
        SVT_LOG("\nSVT [config]: Unknown Profile\t");

    if (config->tier != 0 && config->level != 0)
        SVT_LOG("Tier %d\tLevel %.1f\t", config->tier, (float)(config->level / 10));
    else {
        if (config->tier == 0)
            SVT_LOG("Tier (auto)\t");
        else
            SVT_LOG("Tier %d\t", config->tier);

        if (config->level == 0)
            SVT_LOG("Level (auto)\t");
        else
            SVT_LOG("Level %.1f\t", (float)(config->level / 10));
    }
    SVT_LOG("\nSVT [config]: EncoderMode \t\t\t\t\t\t\t: %d ", config->enc_mode);
    SVT_LOG("\nSVT [config]: EncoderBitDepth / EncoderColorFormat / CompressedTenBitFormat\t: %d / %d / %d", config->encoder_bit_depth, config->encoder_color_format, config->compressed_ten_bit_format);
    SVT_LOG("\nSVT [config]: SourceWidth / SourceHeight\t\t\t\t\t: %d / %d ", config->source_width, config->source_height);
    if (config->frame_rate_denominator != 0 && config->frame_rate_numerator != 0)
        SVT_LOG("\nSVT [config]: Fps_Numerator / Fps_Denominator / Gop Size / IntraRefreshType \t: %d / %d / %d / %d", config->frame_rate_numerator > (1 << 16) ? config->frame_rate_numerator >> 16 : config->frame_rate_numerator,
            config->frame_rate_denominator > (1 << 16) ? config->frame_rate_denominator >> 16 : config->frame_rate_denominator,
            config->intra_period_length + 1,
            config->intra_refresh_type);
    else
        SVT_LOG("\nSVT [config]: FrameRate / Gop Size\t\t\t\t\t\t: %d / %d ", config->frame_rate > 1000 ? config->frame_rate >> 16 : config->frame_rate, config->intra_period_length + 1);
    SVT_LOG("\nSVT [config]: HierarchicalLevels  / PredStructure\t\t: %d / %d", config->hierarchical_levels, config->pred_structure);
    if (config->rate_control_mode == 1)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate (kbps)/ LookaheadDistance / SceneChange\t\t: VBR / %d / %d / %d ", (int)config->target_bit_rate/1000, config->look_ahead_distance, config->scene_change_detection);
    else if (config->rate_control_mode == 2)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate (kbps)/ LookaheadDistance / SceneChange\t\t: Constraint VBR / %d / %d / %d ", (int)config->target_bit_rate/1000, config->look_ahead_distance, config->scene_change_detection);
    else
        SVT_LOG("\nSVT [config]: BRC Mode / QP  / LookaheadDistance / SceneChange\t\t\t: CQP / %d / %d / %d ", scs->static_config.qp, config->look_ahead_distance, config->scene_change_detection);
#ifdef DEBUG_BUFFERS
    SVT_LOG("\nSVT [config]: INPUT / OUTPUT \t\t\t\t\t\t\t: %d / %d", scs->input_buffer_fifo_init_count, scs->output_stream_buffer_fifo_init_count);
    SVT_LOG("\nSVT [config]: CPCS / PAREF / REF \t\t\t\t\t\t: %d / %d / %d", scs->picture_control_set_pool_init_count_child, scs->pa_reference_picture_buffer_init_count, scs->reference_picture_buffer_init_count);
    SVT_LOG("\nSVT [config]: ME_SEG_W0 / ME_SEG_W1 / ME_SEG_W2 / ME_SEG_W3 \t\t\t: %d / %d / %d / %d ",
        scs->me_segment_column_count_array[0],
        scs->me_segment_column_count_array[1],
        scs->me_segment_column_count_array[2],
        scs->me_segment_column_count_array[3]);
    SVT_LOG("\nSVT [config]: ME_SEG_H0 / ME_SEG_H1 / ME_SEG_H2 / ME_SEG_H3 \t\t\t: %d / %d / %d / %d ",
        scs->me_segment_row_count_array[0],
        scs->me_segment_row_count_array[1],
        scs->me_segment_row_count_array[2],
        scs->me_segment_row_count_array[3]);
    SVT_LOG("\nSVT [config]: ME_SEG_W0 / ME_SEG_W1 / ME_SEG_W2 / ME_SEG_W3 \t\t\t: %d / %d / %d / %d ",
        scs->enc_dec_segment_col_count_array[0],
        scs->enc_dec_segment_col_count_array[1],
        scs->enc_dec_segment_col_count_array[2],
        scs->enc_dec_segment_col_count_array[3]);
    SVT_LOG("\nSVT [config]: ME_SEG_H0 / ME_SEG_H1 / ME_SEG_H2 / ME_SEG_H3 \t\t\t: %d / %d / %d / %d ",
        scs->enc_dec_segment_row_count_array[0],
        scs->enc_dec_segment_row_count_array[1],
        scs->enc_dec_segment_row_count_array[2],
        scs->enc_dec_segment_row_count_array[3]);
    SVT_LOG("\nSVT [config]: PA_P / ME_P / SBO_P / MDC_P / ED_P / EC_P \t\t\t: %d / %d / %d / %d / %d / %d ",
        scs->picture_analysis_process_init_count,
        scs->motion_estimation_process_init_count,
        scs->source_based_operations_process_init_count,
        scs->mode_decision_configuration_process_init_count,
        scs->enc_dec_process_init_count,
        scs->entropy_coding_process_init_count);
    SVT_LOG("\nSVT [config]: DLF_P / CDEF_P / REST_P \t\t\t\t\t\t: %d / %d / %d",
        scs->dlf_process_init_count,
        scs->cdef_process_init_count,
        scs->rest_process_init_count);
#endif
    SVT_LOG("\n------------------------------------------- ");
    SVT_LOG("\n");

    fflush(stdout);
}
/**********************************

* Set Parameter
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_set_parameter(
    EbComponentType              *svt_enc_component,
    EbSvtAv1EncConfiguration     *config_struct)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;

    EbErrorType           return_error  = EB_ErrorNone;
    EbEncHandle        *enc_handle  = (EbEncHandle*)svt_enc_component->p_component_private;
    uint32_t              instance_index = 0;

    // Acquire Config Mutex
    eb_block_on_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);

    set_default_configuration_parameters(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    copy_api_from_app(
        enc_handle->scs_instance_array[instance_index]->scs_ptr,
        (EbSvtAv1EncConfiguration*)config_struct);

    return_error = (EbErrorType)verify_settings(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    if (return_error == EB_ErrorBadParameter)
        return EB_ErrorBadParameter;
    set_param_based_on_input(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    // Initialize the Prediction Structure Group
    EB_NO_THROW_NEW(
        enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        prediction_structure_group_ctor,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.enc_mode,
        &(enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config));
    if (!enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr) {
        eb_release_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);
        return EB_ErrorInsufficientResources;
    }
    // Set the Prediction Structure
    enc_handle->scs_instance_array[instance_index]->scs_ptr->pred_struct_ptr = get_prediction_structure(
        enc_handle->scs_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->static_config.pred_structure,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->max_ref_count,
        enc_handle->scs_instance_array[instance_index]->scs_ptr->max_temporal_layers);

    return_error = load_default_buffer_configuration_settings(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    print_lib_params(
        enc_handle->scs_instance_array[instance_index]->scs_ptr);

    // Release Config Mutex
    eb_release_mutex(enc_handle->scs_instance_array[instance_index]->config_mutex);

    return return_error;
}
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
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

    memset(&bitstream, 0, sizeof(Bitstream));
    memset(&output_bitstream, 0, sizeof(OutputBitstreamUnit));
    bitstream.output_bitstream_ptr = &output_bitstream;

    output_stream_buffer = (EbBufferHeaderType *)malloc(sizeof(EbBufferHeaderType));
    if (!output_stream_buffer) {
        return EB_ErrorInsufficientResources;
    }

    output_stream_buffer->p_buffer = (uint8_t *)malloc(sizeof(uint8_t) * PACKETIZATION_PROCESS_BUFFER_SIZE);
    if (!output_stream_buffer->p_buffer) {
        free(output_stream_buffer);
        return EB_ErrorInsufficientResources;
    }

    output_stream_buffer->size = sizeof(EbBufferHeaderType);
    output_stream_buffer->n_alloc_len = PACKETIZATION_PROCESS_BUFFER_SIZE;
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
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
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
//
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_eos_nal(
    EbComponentType           *svt_enc_component,
    EbBufferHeaderType       **output_stream_ptr
)
{
    EbErrorType           return_error = EB_ErrorNone;
    UNUSED(svt_enc_component);
    UNUSED(output_stream_ptr);
    return return_error;
}

/***********************************************
**** Copy the input buffer from the
**** sample application to the library buffers
************************************************/
static EbErrorType copy_frame_buffer(
    SequenceControlSet            *scs_ptr,
    uint8_t                          *dst,
    uint8_t                          *src)
{
    EbSvtAv1EncConfiguration          *config = &scs_ptr->static_config;
    EbErrorType                      return_error = EB_ErrorNone;

    EbPictureBufferDesc           *input_picture_ptr = (EbPictureBufferDesc*)dst;
    EbSvtIOFormat                   *input_ptr = (EbSvtIOFormat*)src;
    uint16_t                         input_row_index;
    EbBool                           is_16bit_input = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is_16bit_input) {
        uint32_t     luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding) << is_16bit_input;
        uint32_t     chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1)) << is_16bit_input;
        uint16_t     luma_stride = input_picture_ptr->stride_y << is_16bit_input;
        uint16_t     chroma_stride = input_picture_ptr->stride_cb << is_16bit_input;
        uint16_t     luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right) << is_16bit_input;
        uint16_t     chroma_width = (luma_width >> 1) << is_16bit_input;
        uint16_t     luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        uint16_t     source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t     source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t     source_cb_stride = (uint16_t)(input_ptr->cb_stride);

        //uint16_t     luma_height  = input_picture_ptr->max_height;
        // Y
        for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
            eb_memcpy((input_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                (input_ptr->luma + source_luma_stride * input_row_index),
                luma_width);
        }

        // U
        for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
            eb_memcpy((input_picture_ptr->buffer_cb + chroma_buffer_offset + chroma_stride * input_row_index),
                (input_ptr->cb + (source_cb_stride*input_row_index)),
                chroma_width);
        }

        // V
        for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
            eb_memcpy((input_picture_ptr->buffer_cr + chroma_buffer_offset + chroma_stride * input_row_index),
                (input_ptr->cr + (source_cr_stride*input_row_index)),
                chroma_width);
        }
    }
    else if (is_16bit_input && config->compressed_ten_bit_format == 1)
    {
        {
            uint32_t  luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
            uint32_t  chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
            uint16_t  luma_stride = input_picture_ptr->stride_y;
            uint16_t  chroma_stride = input_picture_ptr->stride_cb;
            uint16_t  luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right);
            uint16_t  chroma_width = (luma_width >> 1);
            uint16_t  luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

            uint16_t  source_luma_stride = (uint16_t)(input_ptr->y_stride);
            uint16_t  source_cr_stride = (uint16_t)(input_ptr->cr_stride);
            uint16_t  source_cb_stride = (uint16_t)(input_ptr->cb_stride);

            // Y 8bit
            for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
                eb_memcpy((input_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                    (input_ptr->luma + source_luma_stride * input_row_index),
                    luma_width);
            }

            // U 8bit
            for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
                eb_memcpy((input_picture_ptr->buffer_cb + chroma_buffer_offset + chroma_stride * input_row_index),
                    (input_ptr->cb + (source_cb_stride*input_row_index)),
                    chroma_width);
            }

            // V 8bit
            for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
                eb_memcpy((input_picture_ptr->buffer_cr + chroma_buffer_offset + chroma_stride * input_row_index),
                    (input_ptr->cr + (source_cr_stride*input_row_index)),
                    chroma_width);
            }

            //efficient copy - final
            //compressed 2Bit in 1D format
            {
                uint16_t luma_2bit_width = scs_ptr->max_input_luma_width / 4;
                uint16_t luma_height = scs_ptr->max_input_luma_height;

                uint16_t source_luma_2bit_stride = source_luma_stride / 4;
                uint16_t source_chroma_2bit_stride = source_luma_2bit_stride >> 1;

                for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
                    eb_memcpy(input_picture_ptr->buffer_bit_inc_y + luma_2bit_width * input_row_index, input_ptr->luma_ext + source_luma_2bit_stride * input_row_index, luma_2bit_width);
                }
                for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
                    eb_memcpy(input_picture_ptr->buffer_bit_inc_cb + (luma_2bit_width >> 1)*input_row_index, input_ptr->cb_ext + source_chroma_2bit_stride * input_row_index, luma_2bit_width >> 1);
                }
                for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
                    eb_memcpy(input_picture_ptr->buffer_bit_inc_cr + (luma_2bit_width >> 1)*input_row_index, input_ptr->cr_ext + source_chroma_2bit_stride * input_row_index, luma_2bit_width >> 1);
                }
            }
        }
    }
    else { // 10bit packed

        uint32_t luma_offset = 0, chroma_offset = 0;
        uint32_t luma_buffer_offset = (input_picture_ptr->stride_y*scs_ptr->top_padding + scs_ptr->left_padding);
        uint32_t chroma_buffer_offset = (input_picture_ptr->stride_cr*(scs_ptr->top_padding >> 1) + (scs_ptr->left_padding >> 1));
        uint16_t luma_width = (uint16_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right);
        uint16_t chroma_width = (luma_width >> 1);
        uint16_t luma_height = (uint16_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        uint16_t source_luma_stride = (uint16_t)(input_ptr->y_stride);
        uint16_t source_cr_stride = (uint16_t)(input_ptr->cr_stride);
        uint16_t source_cb_stride = (uint16_t)(input_ptr->cb_stride);

        un_pack2d(
            (uint16_t*)(input_ptr->luma + luma_offset),
            source_luma_stride,
            input_picture_ptr->buffer_y + luma_buffer_offset,
            input_picture_ptr->stride_y,
            input_picture_ptr->buffer_bit_inc_y + luma_buffer_offset,
            input_picture_ptr->stride_bit_inc_y,
            luma_width,
            luma_height);

        un_pack2d(
            (uint16_t*)(input_ptr->cb + chroma_offset),
            source_cb_stride,
            input_picture_ptr->buffer_cb + chroma_buffer_offset,
            input_picture_ptr->stride_cb,
            input_picture_ptr->buffer_bit_inc_cb + chroma_buffer_offset,
            input_picture_ptr->stride_bit_inc_cb,
            chroma_width,
            (luma_height >> 1));

        un_pack2d(
            (uint16_t*)(input_ptr->cr + chroma_offset),
            source_cr_stride,
            input_picture_ptr->buffer_cr + chroma_buffer_offset,
            input_picture_ptr->stride_cr,
            input_picture_ptr->buffer_bit_inc_cr + chroma_buffer_offset,
            input_picture_ptr->stride_bit_inc_cr,
            chroma_width,
            (luma_height >> 1));
    }
    return return_error;
}
static void copy_input_buffer(
    SequenceControlSet*    sequenceControlSet,
    EbBufferHeaderType*     dst,
    EbBufferHeaderType*     src
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

    // Copy the picture buffer
    if (src->p_buffer != NULL)
        copy_frame_buffer(sequenceControlSet, dst->p_buffer, src->p_buffer);
}

/**********************************
* Empty This Buffer
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_enc_send_picture(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbEncHandle          *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr;

    // Take the buffer and put it into our internal queue structure
    eb_get_empty_object(
        enc_handle_ptr->input_buffer_producer_fifo_ptr,
        &eb_wrapper_ptr);

    if (p_buffer != NULL) {
        copy_input_buffer(
            enc_handle_ptr->scs_instance_array[0]->scs_ptr,
            (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr,
            p_buffer);
    }

    eb_post_full_object(eb_wrapper_ptr);

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
    if (src->p_buffer)
        eb_memcpy(dst->p_buffer, src->p_buffer, src->n_filled_len);

    return;
}

/**********************************
* svt_av1_enc_get_packet sends out packet
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
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
        eb_get_full_object(
            enc_handle->output_stream_buffer_consumer_fifo_ptr,
            &eb_wrapper_ptr);
    else
        eb_get_full_object_non_blocking(
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

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API void svt_av1_enc_release_out_buffer(
    EbBufferHeaderType  **p_buffer)
{
    if (p_buffer && (*p_buffer)->wrapper_ptr)
    {
        if((*p_buffer)->p_buffer)
           free((*p_buffer)->p_buffer);
        // Release out put buffer back into the pool
        eb_release_object((EbObjectWrapper  *)(*p_buffer)->wrapper_ptr);
     }
    return;
}

/**********************************
* Fill This Buffer
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType svt_av1_get_recon(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbErrorType           return_error = EB_ErrorNone;
    EbEncHandle          *enc_handle = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *eb_wrapper_ptr = NULL;

    if (enc_handle->scs_instance_array[0]->scs_ptr->static_config.recon_enabled) {
        eb_get_full_object_non_blocking(
            enc_handle->output_recon_buffer_consumer_fifo_ptr,
            &eb_wrapper_ptr);

        if (eb_wrapper_ptr) {
            EbBufferHeaderType* obj_ptr = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;
            copy_output_recon_buffer(
                p_buffer,
                obj_ptr);

            if (p_buffer->flags != EB_BUFFERFLAG_EOS && p_buffer->flags != 0)
                return_error = EB_ErrorMax;
            eb_release_object((EbObjectWrapper  *)eb_wrapper_ptr);
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

    eb_get_empty_object(
        enc_handle->output_stream_buffer_consumer_fifo_ptr,
        &eb_wrapper_ptr);

    output_packet            = (EbBufferHeaderType*)eb_wrapper_ptr->object_ptr;

    output_packet->size     = 0;
    output_packet->flags    = error_code;
    output_packet->p_buffer   = NULL;

    eb_post_full_object(eb_wrapper_ptr);
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

    SVT_LOG("SVT [version]:\tSVT-AV1 Encoder Lib v%d.%d.%d\n", SVT_VERSION_MAJOR, SVT_VERSION_MINOR, SVT_VERSION_PATCHLEVEL);
#if ( defined( _MSC_VER ) && (_MSC_VER >= 1920) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2019");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1910) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2017");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1900) )
    SVT_LOG("SVT [build]  :\tVisual Studio 2015");
#elif ( defined( _MSC_VER ) && (_MSC_VER < 1900) )
    SVT_LOG("SVT [build]  :\tVisual Studio (old)");
#elif defined(__GNUC__)
    SVT_LOG("SVT [build]  :\tGCC %s\t", __VERSION__);
#else
    SVT_LOG("SVT [build]  :\tunknown compiler");
#endif
    SVT_LOG(" %u bit\n", (unsigned) sizeof(void*) * 8);
    SVT_LOG("LIB Build date: %s %s\n", __DATE__, __TIME__);
    SVT_LOG("-------------------------------------------\n");

    switch_to_real_time();

    // Set Component Size & Version
    svt_enc_component->size = sizeof(EbComponentType);

    EB_NEW(handle, eb_enc_handle_ctor, svt_enc_component);
    svt_enc_component->p_component_private = handle;

    return return_error;
}
static EbErrorType allocate_frame_buffer(
    SequenceControlSet       *scs_ptr,
    EbBufferHeaderType        *input_buffer)
{
    EbErrorType   return_error = EB_ErrorNone;
    EbPictureBufferDescInitData input_pic_buf_desc_init_data;
    EbSvtAv1EncConfiguration   * config = &scs_ptr->static_config;
    uint8_t is_16bit = config->encoder_bit_depth > 8 ? 1 : 0;
    // Init Picture Init data
    input_pic_buf_desc_init_data.max_width = (uint16_t)scs_ptr->max_input_luma_width;
    input_pic_buf_desc_init_data.max_height = (uint16_t)scs_ptr->max_input_luma_height;
    input_pic_buf_desc_init_data.bit_depth = (EbBitDepthEnum)config->encoder_bit_depth;
    input_pic_buf_desc_init_data.color_format = (EbColorFormat)config->encoder_color_format;

    if (config->compressed_ten_bit_format == 1)
        input_pic_buf_desc_init_data.buffer_enable_mask = 0;
    else
        input_pic_buf_desc_init_data.buffer_enable_mask = is_16bit ?
             PICTURE_BUFFER_DESC_FULL_MASK : 0;
    input_pic_buf_desc_init_data.left_padding = scs_ptr->left_padding;
    input_pic_buf_desc_init_data.right_padding = scs_ptr->right_padding;
    input_pic_buf_desc_init_data.top_padding = scs_ptr->top_padding;
    input_pic_buf_desc_init_data.bot_padding = scs_ptr->bot_padding;

    input_pic_buf_desc_init_data.split_mode = is_16bit ? EB_TRUE : EB_FALSE;

    input_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    if (is_16bit && config->compressed_ten_bit_format == 1)
        //do special allocation for 2bit data down below.
        input_pic_buf_desc_init_data.split_mode = EB_FALSE;

    // Enhanced Picture Buffer
    {
        EbPictureBufferDesc* buf;
        EB_NEW(
            buf,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&input_pic_buf_desc_init_data);
        input_buffer->p_buffer = (uint8_t*)buf;
        if (is_16bit && config->compressed_ten_bit_format == 1) {
            //pack 4 2bit pixels into 1Byte
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_y,
                 (input_pic_buf_desc_init_data.max_width / 4)*
                 (input_pic_buf_desc_init_data.max_height));
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_cb,
                 (input_pic_buf_desc_init_data.max_width / 8)*
                 (input_pic_buf_desc_init_data.max_height / 2));
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_cr,
                 (input_pic_buf_desc_init_data.max_width / 8)*
                 (input_pic_buf_desc_init_data.max_height / 2));
        }
    }

    return return_error;
}
/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType eb_input_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbBufferHeaderType* input_buffer;
    SequenceControlSet        *scs_ptr = (SequenceControlSet*)object_init_data_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(input_buffer, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)input_buffer;
    // Initialize Header
    input_buffer->size = sizeof(EbBufferHeaderType);

    return_error = allocate_frame_buffer(
        scs_ptr,
        input_buffer);
    if (return_error != EB_ErrorNone)
        return return_error;

    input_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

void eb_input_buffer_header_destroyer(    EbPtr p)
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
EbErrorType eb_output_buffer_header_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    EbSvtAv1EncConfiguration   * config = (EbSvtAv1EncConfiguration*)object_init_data_ptr;
    uint32_t n_stride =
        (uint32_t)(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(config->source_width * config->source_height));
    EbBufferHeaderType* out_buf_ptr;

    *object_dbl_ptr = NULL;
    EB_CALLOC(out_buf_ptr, 1, sizeof(EbBufferHeaderType));
    *object_dbl_ptr = (EbPtr)out_buf_ptr;

    // Initialize Header
    out_buf_ptr->size = sizeof(EbBufferHeaderType);
    out_buf_ptr->n_alloc_len = n_stride;
    out_buf_ptr->p_app_private = NULL;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

void eb_output_buffer_header_destroyer(    EbPtr p)
{
    EbBufferHeaderType* obj = (EbBufferHeaderType*)p;
    EB_FREE(obj);
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType eb_output_recon_buffer_header_creator(
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

void eb_output_recon_buffer_header_destroyer(    EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EB_FREE(obj->p_buffer);
    EB_FREE(obj);
}
// clang-format on
