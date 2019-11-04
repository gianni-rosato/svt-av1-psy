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
#include <immintrin.h>

#include "EbDefinitions.h"
#include "EbSvtAv1Enc.h"
#include "EbThreads.h"
#include "EbUtility.h"
#include "EbString.h"
#include "EbEncHandle.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbPictureOperators.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"
#include "EbReferenceObject.h"
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
#include "EbRateControlResults.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbPredictionStructure.h"
#include "EbDlfProcess.h"
#include "EbCdefProcess.h"
#include "EbRestProcess.h"
#include "EbObject.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#endif

#include "aom_dsp_rtcd.h"

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
        char *name;
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
uint32_t GetNumProcessors() {
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return num_groups == 1 ? sysinfo.dwNumberOfProcessors : sysinfo.dwNumberOfProcessors << 1;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

EbErrorType InitThreadManagmentParams() {
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
uint64_t GetAffinityMask(uint32_t lpnum) {
    uint64_t mask = 0x1;
    for (uint32_t i = lpnum - 1; i > 0; i--)
        mask += (uint64_t)1 << i;
    return mask;
}
#endif

void EbSetThreadManagementParameters(EbSvtAv1EncConfiguration *config_ptr)
{
#ifdef _WIN32
    uint32_t num_logical_processors = GetNumProcessors();
    // For system with a single processor group(no more than 64 logic processors all together)
    // Affinity of the thread can be set to one or more logical processors
    if (num_groups == 1) {
        uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
            config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
        group_affinity.Mask = GetAffinityMask(lps);
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
                    group_affinity.Mask = GetAffinityMask(config_ptr->logical_processors);
            }
            else {
                uint32_t lps = config_ptr->logical_processors == 0 ? num_lp_per_group :
                    config_ptr->logical_processors < num_lp_per_group ? config_ptr->logical_processors : num_lp_per_group;
                group_affinity.Mask = GetAffinityMask(lps);
                group_affinity.Group = config_ptr->target_socket;
            }
        }
    }
#elif defined(__linux__)
    uint32_t num_logical_processors = GetNumProcessors();
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
#else
    UNUSED(config_ptr);
#endif
}

void asmSetConvolveAsmTable(void);
void asmSetConvolveHbdAsmTable(void);
void init_intra_dc_predictors_c_internal(void);
void init_intra_predictors_internal(void);
void eb_av1_init_me_luts(void);

void SwitchToRealTime(){
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
    SequenceControlSet       *sequence_control_set_ptr){
    EbErrorType           return_error = EB_ErrorNone;
    uint32_t encDecSegH = (sequence_control_set_ptr->static_config.super_block_size == 128) ?
        ((sequence_control_set_ptr->max_input_luma_height + 64) / 128) :
        ((sequence_control_set_ptr->max_input_luma_height + 32) / 64);
    uint32_t encDecSegW = (sequence_control_set_ptr->static_config.super_block_size == 128) ?
        ((sequence_control_set_ptr->max_input_luma_width + 64) / 128) :
        ((sequence_control_set_ptr->max_input_luma_width + 32) / 64);

#if CABAC_SERIAL
    encDecSegH = 1;
    encDecSegW = 1;
#endif

    uint32_t meSegH     = (((sequence_control_set_ptr->max_input_luma_height + 32) / BLOCK_SIZE_64) < 6) ? 1 : 6;
    uint32_t meSegW     = (((sequence_control_set_ptr->max_input_luma_width + 32) / BLOCK_SIZE_64) < 10) ? 1 : 10;

    unsigned int lp_count   = GetNumProcessors();
    unsigned int core_count = lp_count;
#if defined(_WIN32) || defined(__linux__)
    if (sequence_control_set_ptr->static_config.target_socket != -1)
        core_count /= num_groups;
#endif
    if (sequence_control_set_ptr->static_config.logical_processors != 0)
        core_count = sequence_control_set_ptr->static_config.logical_processors < core_count ?
            sequence_control_set_ptr->static_config.logical_processors: core_count;

#ifdef _WIN32
    //Handle special case on Windows
    //By default, on Windows an application is constrained to a single group
    if (sequence_control_set_ptr->static_config.target_socket == -1 &&
        sequence_control_set_ptr->static_config.logical_processors == 0)
        core_count /= num_groups;

    //Affininty can only be set by group on Windows.
    //Run on both sockets if -lp is larger than logical processor per group.
    if (sequence_control_set_ptr->static_config.target_socket == -1 &&
        sequence_control_set_ptr->static_config.logical_processors > lp_count / num_groups)
        core_count = lp_count;
#endif
    int32_t return_ppcs = set_parent_pcs(&sequence_control_set_ptr->static_config,
                    core_count, sequence_control_set_ptr->input_resolution);
    if (return_ppcs == -1)
        return EB_ErrorInsufficientResources;
    uint32_t input_pic = (uint32_t)return_ppcs;
    sequence_control_set_ptr->input_buffer_fifo_init_count = input_pic + SCD_LAD + sequence_control_set_ptr->static_config.look_ahead_distance;
    sequence_control_set_ptr->output_stream_buffer_fifo_init_count =
        sequence_control_set_ptr->input_buffer_fifo_init_count + 4;

    // ME segments
    sequence_control_set_ptr->me_segment_row_count_array[0] = meSegH;
    sequence_control_set_ptr->me_segment_row_count_array[1] = meSegH;
    sequence_control_set_ptr->me_segment_row_count_array[2] = meSegH;
    sequence_control_set_ptr->me_segment_row_count_array[3] = meSegH;
    sequence_control_set_ptr->me_segment_row_count_array[4] = meSegH;
    sequence_control_set_ptr->me_segment_row_count_array[5] = meSegH;

    sequence_control_set_ptr->me_segment_column_count_array[0] = meSegW;
    sequence_control_set_ptr->me_segment_column_count_array[1] = meSegW;
    sequence_control_set_ptr->me_segment_column_count_array[2] = meSegW;
    sequence_control_set_ptr->me_segment_column_count_array[3] = meSegW;
    sequence_control_set_ptr->me_segment_column_count_array[4] = meSegW;
    sequence_control_set_ptr->me_segment_column_count_array[5] = meSegW;

    // EncDec segments
    sequence_control_set_ptr->enc_dec_segment_row_count_array[0] = encDecSegH;
    sequence_control_set_ptr->enc_dec_segment_row_count_array[1] = encDecSegH;
    sequence_control_set_ptr->enc_dec_segment_row_count_array[2] = encDecSegH;
    sequence_control_set_ptr->enc_dec_segment_row_count_array[3] = encDecSegH;
    sequence_control_set_ptr->enc_dec_segment_row_count_array[4] = encDecSegH;
    sequence_control_set_ptr->enc_dec_segment_row_count_array[5] = encDecSegH;

    sequence_control_set_ptr->enc_dec_segment_col_count_array[0] = encDecSegW;
    sequence_control_set_ptr->enc_dec_segment_col_count_array[1] = encDecSegW;
    sequence_control_set_ptr->enc_dec_segment_col_count_array[2] = encDecSegW;
    sequence_control_set_ptr->enc_dec_segment_col_count_array[3] = encDecSegW;
    sequence_control_set_ptr->enc_dec_segment_col_count_array[4] = encDecSegW;
    sequence_control_set_ptr->enc_dec_segment_col_count_array[5] = encDecSegW;

    sequence_control_set_ptr->cdef_segment_column_count = meSegW;
    sequence_control_set_ptr->cdef_segment_row_count    = meSegH;

    //since restoration unit size is same for Luma and Chroma, Luma segments and chroma segments do not correspond to the same area!
    //to keep proper processing, segments have to be configured based on chroma resolution.
    uint32_t unit_size                                  = 256;
    uint32_t rest_seg_w                                 = MAX((sequence_control_set_ptr->max_input_luma_width /2 + (unit_size >> 1)) / unit_size, 1);
    uint32_t rest_seg_h                                 = MAX((sequence_control_set_ptr->max_input_luma_height/2 + (unit_size >> 1)) / unit_size, 1);
    sequence_control_set_ptr->rest_segment_column_count = MIN(rest_seg_w,6);
    sequence_control_set_ptr->rest_segment_row_count    = MIN(rest_seg_h,4);

    sequence_control_set_ptr->tf_segment_column_count = meSegW;//1;//
    sequence_control_set_ptr->tf_segment_row_count =  meSegH;//1;//
    //#====================== Data Structures and Picture Buffers ======================
    sequence_control_set_ptr->picture_control_set_pool_init_count       = input_pic + SCD_LAD + sequence_control_set_ptr->static_config.look_ahead_distance;
    if (sequence_control_set_ptr->static_config.enable_overlays)
        sequence_control_set_ptr->picture_control_set_pool_init_count = MAX(sequence_control_set_ptr->picture_control_set_pool_init_count,
            sequence_control_set_ptr->static_config.look_ahead_distance + // frames in the LAD
            sequence_control_set_ptr->static_config.look_ahead_distance / (1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 1 +  // number of overlayes in the LAD
            ((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + SCD_LAD) * 2 +// minigop formation in PD + SCD_LAD *(normal pictures + potential pictures )
            (1 << sequence_control_set_ptr->static_config.hierarchical_levels)); // minigop in PM
    sequence_control_set_ptr->picture_control_set_pool_init_count_child = MAX(MAX(MIN(3, core_count/2), core_count / 6), 1);
    sequence_control_set_ptr->reference_picture_buffer_init_count       = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          sequence_control_set_ptr->static_config.look_ahead_distance + SCD_LAD;
    sequence_control_set_ptr->pa_reference_picture_buffer_init_count    = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          sequence_control_set_ptr->static_config.look_ahead_distance + SCD_LAD;
    sequence_control_set_ptr->output_recon_buffer_fifo_init_count       = sequence_control_set_ptr->reference_picture_buffer_init_count;
    sequence_control_set_ptr->overlay_input_picture_buffer_init_count   = sequence_control_set_ptr->static_config.enable_overlays ?
                                                                          (2 << sequence_control_set_ptr->static_config.hierarchical_levels) + SCD_LAD : 1;

    //#====================== Inter process Fifos ======================
    sequence_control_set_ptr->resource_coordination_fifo_init_count       = 300;
    sequence_control_set_ptr->picture_analysis_fifo_init_count            = 300;
    sequence_control_set_ptr->picture_decision_fifo_init_count            = 300;
    sequence_control_set_ptr->initial_rate_control_fifo_init_count        = 300;
    sequence_control_set_ptr->picture_demux_fifo_init_count               = 300;
    sequence_control_set_ptr->rate_control_tasks_fifo_init_count          = 300;
    sequence_control_set_ptr->rate_control_fifo_init_count                = 301;
    sequence_control_set_ptr->mode_decision_configuration_fifo_init_count = 300;
    sequence_control_set_ptr->motion_estimation_fifo_init_count           = 300;
    sequence_control_set_ptr->entropy_coding_fifo_init_count              = 300;
    sequence_control_set_ptr->enc_dec_fifo_init_count                     = 300;
    sequence_control_set_ptr->dlf_fifo_init_count                         = 300;
    sequence_control_set_ptr->cdef_fifo_init_count                        = 300;
    sequence_control_set_ptr->rest_fifo_init_count                        = 300;
    //#====================== Processes number ======================
    sequence_control_set_ptr->total_process_init_count                    = 0;
    if (core_count > 1){
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->picture_analysis_process_init_count            = MAX(MIN(15, core_count >> 1), core_count / 6));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->motion_estimation_process_init_count =  MAX(MIN(20, core_count >> 1), core_count / 3));//1);//
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->source_based_operations_process_init_count     = MAX(MIN(3, core_count >> 1), core_count / 12));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->mode_decision_configuration_process_init_count = MAX(MIN(3, core_count >> 1), core_count / 12));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->enc_dec_process_init_count                     = MAX(MIN(40, core_count >> 1), core_count));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->entropy_coding_process_init_count              = MAX(MIN(3, core_count >> 1), core_count / 12));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->dlf_process_init_count                         = MAX(MIN(40, core_count >> 1), core_count));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->cdef_process_init_count                        = MAX(MIN(40, core_count >> 1), core_count));
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->rest_process_init_count                        = MAX(MIN(40, core_count >> 1), core_count));
    }else{
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->picture_analysis_process_init_count            = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->motion_estimation_process_init_count           = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->source_based_operations_process_init_count     = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->mode_decision_configuration_process_init_count = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->enc_dec_process_init_count                     = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->entropy_coding_process_init_count              = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->dlf_process_init_count                         = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->cdef_process_init_count                        = 1);
        sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->rest_process_init_count                        = 1);
    }

    sequence_control_set_ptr->total_process_init_count += 6; // single processes count
    printf("Number of logical cores available: %u\nNumber of PPCS %u\n", core_count, sequence_control_set_ptr->picture_control_set_pool_init_count);

    /******************************************************************
    * Platform detection, limit cpu flags to hardware available CPU
    ******************************************************************/
    const CPU_FLAGS cpu_flags = get_cpu_flags();
    const CPU_FLAGS cpu_flags_to_use = get_cpu_flags_to_use();
    sequence_control_set_ptr->static_config.use_cpu_flags &= cpu_flags_to_use;
    printf("[asm level on system : up to %s]\n", get_asm_level_name_str(cpu_flags));
    printf("[asm level selected : up to %s]\n", get_asm_level_name_str(sequence_control_set_ptr->static_config.use_cpu_flags));

    return return_error;
}
 // Rate Control
static RateControlPorts rateControlPorts[] = {
    {RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER,       0},
    {RATE_CONTROL_INPUT_PORT_PACKETIZATION,         0},
    {RATE_CONTROL_INPUT_PORT_ENTROPY_CODING,        0},
    {RATE_CONTROL_INPUT_PORT_INVALID,               0}
};
// Rate Control
static uint32_t RateControlPortLookup(
    RateControlInputPortTypes           type,
    uint32_t                                portTypeIndex){
    uint32_t portIndex = 0;
    uint32_t portCount = 0;

    while ((type != rateControlPorts[portIndex].type) && (type != RATE_CONTROL_INPUT_PORT_INVALID))
        portCount += rateControlPorts[portIndex++].count;
    return (portCount + portTypeIndex);
}
// Rate Control
static uint32_t RateControlPortTotalCount(void){
    uint32_t portIndex = 0;
    uint32_t total_count = 0;

    while (rateControlPorts[portIndex].type != RATE_CONTROL_INPUT_PORT_INVALID)
        total_count += rateControlPorts[portIndex++].count;
    return total_count;
}

// EncDec
typedef struct {
    int32_t  type;
    uint32_t  count;
} EncDecPorts_t;
static EncDecPorts_t encDecPorts[] = {
    {ENCDEC_INPUT_PORT_MDC,        0},
    {ENCDEC_INPUT_PORT_ENCDEC,     0},
    {ENCDEC_INPUT_PORT_INVALID,    0}
};

/*****************************************
 * Input Port Lookup
 *****************************************/
// EncDec
static uint32_t EncDecPortLookup(
    int32_t  type,
    uint32_t  portTypeIndex)
{
    uint32_t portIndex = 0;
    uint32_t portCount = 0;

    while ((type != encDecPorts[portIndex].type) && (type != ENCDEC_INPUT_PORT_INVALID))
        portCount += encDecPorts[portIndex++].count;
    return (portCount + portTypeIndex);
}
// EncDec
static uint32_t EncDecPortTotalCount(void){
    uint32_t portIndex = 0;
    uint32_t total_count = 0;

    while (encDecPorts[portIndex].type != ENCDEC_INPUT_PORT_INVALID)
        total_count += encDecPorts[portIndex++].count;
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
    SequenceControlSet*  control_set_ptr = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr;
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
    EB_DELETE(enc_handle_ptr->sequence_control_set_pool_ptr);
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
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->picture_analysis_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->motion_estimation_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->source_based_operations_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->mode_decision_configuration_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->enc_dec_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->dlf_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->cdef_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->rest_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->entropy_coding_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->sequence_control_set_instance_array, enc_handle_ptr->encode_instance_total_count);
    EB_DELETE(enc_handle_ptr->picture_decision_context_ptr);
    EB_DELETE(enc_handle_ptr->initial_rate_control_context_ptr);
    EB_DELETE(enc_handle_ptr->picture_manager_context_ptr);
    EB_DELETE(enc_handle_ptr->rate_control_context_ptr);
    EB_DELETE(enc_handle_ptr->packetization_context_ptr);
    EB_DELETE_PTR_ARRAY(enc_handle_ptr->reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_FREE_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array);
    EB_FREE_ARRAY(enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr_dbl_array);

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

    return_error = InitThreadManagmentParams();

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
    enc_handle_ptr->encode_instance_total_count                           = EB_EncodeInstancesTotalCount;
    enc_handle_ptr->compute_segments_total_count_array                    = EB_ComputeSegmentInitCount;
    // Config Set Count
    enc_handle_ptr->sequence_control_set_pool_total_count                 = EB_SequenceControlSetPoolInitCount;

    // Initialize Callbacks
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->app_callback_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC(enc_handle_ptr->app_callback_ptr_array[0], sizeof(EbCallback));
    enc_handle_ptr->app_callback_ptr_array[0]->ErrorHandler = lib_svt_encoder_send_error_exit;
    enc_handle_ptr->app_callback_ptr_array[0]->handle = ebHandlePtr;

    // Initialize Sequence Control Set Instance Array
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->sequence_control_set_instance_array, enc_handle_ptr->encode_instance_total_count);
    EB_NEW(enc_handle_ptr->sequence_control_set_instance_array[0], eb_sequence_control_set_instance_ctor);
    return EB_ErrorNone;
}

EbErrorType EbInputBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr);

EbErrorType EbOutputReconBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr);

EbErrorType EbOutputBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr objectInitDataPtr);

void EbInputBufferHeaderDestoryer(    EbPtr p);
void EbOutputReconBufferHeaderDestoryer(    EbPtr p);
void EbOutputBufferHeaderDestoryer(    EbPtr p);


EbErrorType DlfResultsCtor(
    DlfResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType DlfResultsCreator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    DlfResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, DlfResultsCtor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType CdefResultsCtor(
    CdefResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType CdefResultsCreator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    CdefResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, CdefResultsCtor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType RestResultsCtor(
    RestResults *context_ptr,
    EbPtr object_init_data_ptr)
{
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType RestResultsCreator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    RestResults* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, RestResultsCtor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

void init_fn_ptr(void);
extern void av1_init_wedge_masks(void);
/**********************************
* Initialize Encoder Library
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_init_encoder(EbComponentType *svt_enc_component)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;
    EbEncHandle *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbErrorType return_error = EB_ErrorNone;
    uint32_t instance_index;
    uint32_t processIndex;
    uint32_t max_picture_width;
    EbBool is16bit = (EbBool)(enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.encoder_color_format;
    SequenceControlSet* control_set_ptr;

    setup_rtcd_internal(enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.use_cpu_flags);
    asmSetConvolveAsmTable();

    init_intra_dc_predictors_c_internal();

    asmSetConvolveHbdAsmTable();

    init_intra_predictors_internal();
    EbSequenceControlSetInitData scs_init;
    scs_init.sb_size = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.super_block_size;

    build_blk_geom(scs_init.sb_size == 128);

    eb_av1_init_me_luts();
    init_fn_ptr();
    av1_init_wedge_masks();
    /************************************
    * Sequence Control Set
    ************************************/
    EB_NEW(enc_handle_ptr->sequence_control_set_pool_ptr,
        eb_system_resource_ctor,
        enc_handle_ptr->sequence_control_set_pool_total_count,
        1,
        0,
        &enc_handle_ptr->sequence_control_set_pool_producer_fifo_ptr_array,
        (EbFifo ***)EB_NULL,
        EB_FALSE,
        eb_sequence_control_set_creator,
        &scs_init,
        NULL);

    /************************************
    * Picture Control Set: Parent
    ************************************/
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_MALLOC_ARRAY(enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        // The segment Width & Height Arrays are in units of LCUs, not samples
        PictureControlSetInitData inputData;

        inputData.picture_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
        inputData.picture_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height;
        inputData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->left_padding;
        inputData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->right_padding;
        inputData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->top_padding;
        inputData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->bot_padding;
        inputData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        inputData.color_format = color_format;
        inputData.sb_sz = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz;
        inputData.max_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_sb_depth;
        inputData.ten_bit_format = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.ten_bit_format;
        inputData.compressed_ten_bit_format = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.compressed_ten_bit_format;
        inputData.enc_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enc_mode;
        inputData.speed_control = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.speed_control_flag;
        inputData.hbd_mode_decision = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enable_hbd_mode_decision;
        inputData.film_grain_noise_level = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.film_grain_denoise_strength;
        inputData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.encoder_bit_depth;

        inputData.ext_block_flag = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.ext_block_flag;

        inputData.in_loop_me_flag = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.in_loop_me_flag;
        inputData.mrp_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mrp_mode;
        inputData.nsq_present = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->nsq_present;
        EB_NEW(
            enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->picture_control_set_pool_init_count,//enc_handle_ptr->picture_control_set_pool_total_count,
            1,
            0,
            &enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            picture_parent_control_set_creator,
            &inputData,
            NULL);
    }

    /************************************
    * Picture Control Set: Child
    ************************************/
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_control_set_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_MALLOC_ARRAY(enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        // The segment Width & Height Arrays are in units of LCUs, not samples
        PictureControlSetInitData inputData;
        unsigned i;

        inputData.enc_dec_segment_col = 0;
        inputData.enc_dec_segment_row = 0;
        for (i = 0; i <= enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.hierarchical_levels; ++i) {
            inputData.enc_dec_segment_col = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_segment_col_count_array[i] > inputData.enc_dec_segment_col ?
                (uint16_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_segment_col_count_array[i] :
                inputData.enc_dec_segment_col;
            inputData.enc_dec_segment_row = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_segment_row_count_array[i] > inputData.enc_dec_segment_row ?
                (uint16_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_segment_row_count_array[i] :
                inputData.enc_dec_segment_row;
        }

        inputData.picture_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
        inputData.picture_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height;
        inputData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->left_padding;
        inputData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->right_padding;
        inputData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->top_padding;
        inputData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->bot_padding;
        inputData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        inputData.film_grain_noise_level = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->film_grain_denoise_strength;
        inputData.color_format = color_format;
        inputData.sb_sz = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz;
        inputData.sb_size_pix = scs_init.sb_size;
        inputData.max_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_sb_depth;
        inputData.hbd_mode_decision = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enable_hbd_mode_decision;
        inputData.cdf_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->cdf_mode;
        inputData.mfmv = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mfmv_enabled;

#if PAL_SUP
        inputData.cfg_palette = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.screen_content_mode;
#endif
        EB_NEW(
            enc_handle_ptr->picture_control_set_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->picture_control_set_pool_init_count_child, //EB_PictureControlSetPoolInitCountChild,
            1,
            0,
            &enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            picture_control_set_creator,
            &inputData,
            NULL);
    }

    /************************************
    * Picture Buffers
    ************************************/

    // Allocate Resource Arrays
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->pa_reference_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);

    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->overlay_input_picture_pool_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC_ARRAY(enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);
    // Allocate Producer Fifo Arrays
    EB_MALLOC_ARRAY(enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC_ARRAY(enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);

    // Rate Control
    rateControlPorts[0].count = EB_PictureManagerProcessInitCount;
    rateControlPorts[1].count = EB_PacketizationProcessInitCount;
    rateControlPorts[2].count = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count;
    rateControlPorts[3].count = 0;

    encDecPorts[ENCDEC_INPUT_PORT_MDC].count = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count;
    encDecPorts[ENCDEC_INPUT_PORT_ENCDEC].count = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count;

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EbReferenceObjectDescInitData     EbReferenceObjectDescInitDataStructure;
        EbPaReferenceObjectDescInitData   EbPaReferenceObjectDescInitDataStructure;
        EbPictureBufferDescInitData       referencePictureBufferDescInitData;
        EbPictureBufferDescInitData       quarterPictureBufferDescInitData;
        EbPictureBufferDescInitData       sixteenthPictureBufferDescInitData;
        // Initialize the various Picture types
        referencePictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
        referencePictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height;
        referencePictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        referencePictureBufferDescInitData.color_format = color_format;
        referencePictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

        referencePictureBufferDescInitData.left_padding = PAD_VALUE;
        referencePictureBufferDescInitData.right_padding = PAD_VALUE;
        referencePictureBufferDescInitData.top_padding = PAD_VALUE;
        referencePictureBufferDescInitData.bot_padding = PAD_VALUE;
        referencePictureBufferDescInitData.mfmv = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mfmv_enabled;
        // Hsan: split_mode is set @ eb_reference_object_ctor() as both unpacked reference and packed reference are needed for a 10BIT input; unpacked reference @ MD, and packed reference @ EP

        if (is16bit)
            referencePictureBufferDescInitData.bit_depth = EB_10BIT;

        EbReferenceObjectDescInitDataStructure.reference_picture_desc_init_data = referencePictureBufferDescInitData;

        // Reference Picture Buffers
        EB_NEW(
            enc_handle_ptr->reference_picture_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->reference_picture_buffer_init_count,//enc_handle_ptr->reference_picture_pool_total_count,
            EB_PictureManagerProcessInitCount,
            0,
            &enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            eb_reference_object_creator,
            &(EbReferenceObjectDescInitDataStructure),
            NULL);

        // PA Reference Picture Buffers
        // Currently, only Luma samples are needed in the PA
        referencePictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
        referencePictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height;
        referencePictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        referencePictureBufferDescInitData.color_format = EB_YUV420; //use 420 for picture analysis
        referencePictureBufferDescInitData.buffer_enable_mask = 0;
        referencePictureBufferDescInitData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz + ME_FILTER_TAP;
        referencePictureBufferDescInitData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz + ME_FILTER_TAP;
        referencePictureBufferDescInitData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz + ME_FILTER_TAP;
        referencePictureBufferDescInitData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz + ME_FILTER_TAP;
        referencePictureBufferDescInitData.split_mode = EB_FALSE;
        quarterPictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width >> 1;
        quarterPictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height >> 1;
        quarterPictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        quarterPictureBufferDescInitData.color_format = EB_YUV420;
        quarterPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        quarterPictureBufferDescInitData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterPictureBufferDescInitData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterPictureBufferDescInitData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterPictureBufferDescInitData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterPictureBufferDescInitData.split_mode = EB_FALSE;
        quarterPictureBufferDescInitData.down_sampled_filtered = (enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;

        sixteenthPictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width >> 2;
        sixteenthPictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height >> 2;
        sixteenthPictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        sixteenthPictureBufferDescInitData.color_format = EB_YUV420;
        sixteenthPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        sixteenthPictureBufferDescInitData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthPictureBufferDescInitData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthPictureBufferDescInitData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthPictureBufferDescInitData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthPictureBufferDescInitData.split_mode = EB_FALSE;
        sixteenthPictureBufferDescInitData.down_sampled_filtered = (enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ? EB_TRUE : EB_FALSE;

        EbPaReferenceObjectDescInitDataStructure.reference_picture_desc_init_data = referencePictureBufferDescInitData;
        EbPaReferenceObjectDescInitDataStructure.quarter_picture_desc_init_data = quarterPictureBufferDescInitData;
        EbPaReferenceObjectDescInitDataStructure.sixteenth_picture_desc_init_data = sixteenthPictureBufferDescInitData;
        // Reference Picture Buffers
        EB_NEW(enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pa_reference_picture_buffer_init_count,
            EB_PictureDecisionProcessInitCount,
            0,
            &enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            eb_pa_reference_object_creator,
            &(EbPaReferenceObjectDescInitDataStructure),
            NULL);
        // Set the SequenceControlSet Picture Pool Fifo Ptrs
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->reference_picture_pool_fifo_ptr = (enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->pa_reference_picture_pool_fifo_ptr = (enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];

        if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.enable_overlays) {
            // Overlay Input Picture Buffers
            EB_NEW(
                enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index],
                eb_system_resource_ctor,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->overlay_input_picture_buffer_init_count,
                1,
                0,
                &enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
                (EbFifo ***)EB_NULL,
                EB_FALSE,
                EbInputBufferHeaderCreator,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr,
                EbInputBufferHeaderDestoryer);
           // Set the SequenceControlSet Overlay input Picture Pool Fifo Ptrs
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->overlay_input_picture_pool_fifo_ptr = (enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];
        }
    }

    /************************************
    * System Resource Managers & Fifos
    ************************************/

    // EbBufferHeaderType Input
    EB_NEW(
        enc_handle_ptr->input_buffer_resource_ptr,
        eb_system_resource_ctor,
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->input_buffer_fifo_init_count,
        1,
        EB_ResourceCoordinationProcessInitCount,
        &enc_handle_ptr->input_buffer_producer_fifo_ptr_array,
        &enc_handle_ptr->input_buffer_consumer_fifo_ptr_array,
        EB_TRUE,
        EbInputBufferHeaderCreator,
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr,
        EbInputBufferHeaderDestoryer);

    // EbBufferHeaderType Output Stream
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_stream_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC_ARRAY(enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);
    EB_MALLOC_ARRAY(enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EB_NEW(
            enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index],
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->output_stream_buffer_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->total_process_init_count,//EB_PacketizationProcessInitCount,
            1,
            &enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array[instance_index],
            &enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array[instance_index],
            EB_TRUE,
            EbOutputBufferHeaderCreator,
            &enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config,
            EbOutputBufferHeaderDestoryer);
    }
    if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.recon_enabled) {
        // EbBufferHeaderType Output Recon
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->output_recon_buffer_resource_ptr_array, enc_handle_ptr->encode_instance_total_count);
        EB_MALLOC_ARRAY(enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);
        EB_MALLOC_ARRAY(enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr_dbl_array, enc_handle_ptr->encode_instance_total_count);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            EB_NEW(
                enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index],
                eb_system_resource_ctor,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->output_recon_buffer_fifo_init_count,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_process_init_count,
                1,
                &enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array[instance_index],
                &enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr_dbl_array[instance_index],
                EB_TRUE,
                EbOutputReconBufferHeaderCreator,
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr,
                EbOutputReconBufferHeaderDestoryer);
        }
    }

    // Resource Coordination Results
    {
        ResourceCoordinationResultInitData resourceCoordinationResultInitData;

        EB_NEW(
            enc_handle_ptr->resource_coordination_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->resource_coordination_fifo_init_count,
            EB_ResourceCoordinationProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count,
            &enc_handle_ptr->resource_coordination_results_producer_fifo_ptr_array,
            &enc_handle_ptr->resource_coordination_results_consumer_fifo_ptr_array,
            EB_TRUE,
            resource_coordination_result_creator,
            &resourceCoordinationResultInitData,
            NULL);
    }

    // Picture Analysis Results
    {
        PictureAnalysisResultInitData pictureAnalysisResultInitData;

        EB_NEW(
            enc_handle_ptr->picture_analysis_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count,
            EB_PictureDecisionProcessInitCount,
            &enc_handle_ptr->picture_analysis_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_analysis_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_analysis_result_creator,
            &pictureAnalysisResultInitData,
            NULL);
    }

    // Picture Decision Results
    {
        PictureDecisionResultInitData pictureDecisionResultInitData;

        EB_NEW(
            enc_handle_ptr->picture_decision_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_decision_fifo_init_count,
            EB_PictureDecisionProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count,
            &enc_handle_ptr->picture_decision_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_decision_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_decision_result_creator,
            &pictureDecisionResultInitData,
            NULL);
    }

    // Motion Estimation Results
    {
        MotionEstimationResultsInitData motionEstimationResultInitData;

        EB_NEW(
            enc_handle_ptr->motion_estimation_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count,
            EB_InitialRateControlProcessInitCount,
            &enc_handle_ptr->motion_estimation_results_producer_fifo_ptr_array,
            &enc_handle_ptr->motion_estimation_results_consumer_fifo_ptr_array,
            EB_TRUE,
            motion_estimation_results_creator,
            &motionEstimationResultInitData,
            NULL);
    }

    // Initial Rate Control Results
    {
        InitialRateControlResultInitData initialRateControlResultInitData;

        EB_NEW(
            enc_handle_ptr->initial_rate_control_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->initial_rate_control_fifo_init_count,
            EB_InitialRateControlProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count,
            &enc_handle_ptr->initial_rate_control_results_producer_fifo_ptr_array,
            &enc_handle_ptr->initial_rate_control_results_consumer_fifo_ptr_array,
            EB_TRUE,
            initial_rate_control_results_creator,
            &initialRateControlResultInitData,
            NULL);
    }

    // Picture Demux Results
    {
        PictureResultInitData pictureResultInitData;
        EB_NEW(
            enc_handle_ptr->picture_demux_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_demux_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count + enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count + 1, // 1 for packetization
            EB_PictureManagerProcessInitCount,
            &enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_demux_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_results_creator,
            &pictureResultInitData,
            NULL);

    }

    // Rate Control Tasks
    {
        RateControlTasksInitData rateControlTasksInitData;

        EB_NEW(
            enc_handle_ptr->rate_control_tasks_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rate_control_tasks_fifo_init_count,
            RateControlPortTotalCount(),
            EB_RateControlProcessInitCount,
            &enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array,
            &enc_handle_ptr->rate_control_tasks_consumer_fifo_ptr_array,
            EB_TRUE,
            rate_control_tasks_creator,
            &rateControlTasksInitData,
            NULL);
    }

    // Rate Control Results
    {
        RateControlResultsInitData rateControlResultInitData;

        EB_NEW(
            enc_handle_ptr->rate_control_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rate_control_fifo_init_count,
            EB_RateControlProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count,
            &enc_handle_ptr->rate_control_results_producer_fifo_ptr_array,
            &enc_handle_ptr->rate_control_results_consumer_fifo_ptr_array,
            EB_TRUE,
            rate_control_results_creator,
            &rateControlResultInitData,
            NULL);
    }
    // EncDec Tasks
    {
        EncDecTasksInitData ModeDecisionResultInitData;
        unsigned i;

        ModeDecisionResultInitData.enc_dec_segment_row_count = 0;

        for (i = 0; i <= enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.hierarchical_levels; ++i) {
            ModeDecisionResultInitData.enc_dec_segment_row_count = MAX(
                ModeDecisionResultInitData.enc_dec_segment_row_count,
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_segment_row_count_array[i]);
        }

        EB_NEW(
            enc_handle_ptr->enc_dec_tasks_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_fifo_init_count,
            EncDecPortTotalCount(),
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count,
            &enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array,
            &enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array,
            EB_TRUE,
            enc_dec_tasks_creator,
            &ModeDecisionResultInitData,
            NULL);
    }

    // EncDec Results
    {
        EncDecResultsInitData encDecResultInitData;

        EB_NEW(
            enc_handle_ptr->enc_dec_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count,
            &enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array,
            &enc_handle_ptr->enc_dec_results_consumer_fifo_ptr_array,
            EB_TRUE,
            enc_dec_results_creator,
            &encDecResultInitData,
            NULL);
   }

    //DLF results
    {
        EntropyCodingResultsInitData dlfResultInitData;

        EB_NEW(
            enc_handle_ptr->dlf_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count,
            &enc_handle_ptr->dlf_results_producer_fifo_ptr_array,
            &enc_handle_ptr->dlf_results_consumer_fifo_ptr_array,
            EB_TRUE,
            DlfResultsCreator,
            &dlfResultInitData,
            NULL);
    }
    //CDEF results
    {
        EntropyCodingResultsInitData cdefResultInitData;

        EB_NEW(
            enc_handle_ptr->cdef_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count,
            &enc_handle_ptr->cdef_results_producer_fifo_ptr_array,
            &enc_handle_ptr->cdef_results_consumer_fifo_ptr_array,
            EB_TRUE,
            CdefResultsCreator,
            &cdefResultInitData,
            NULL);
    }
    //REST results
    {
        EntropyCodingResultsInitData restResultInitData;

        EB_NEW(
            enc_handle_ptr->rest_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count,
            &enc_handle_ptr->rest_results_producer_fifo_ptr_array,
            &enc_handle_ptr->rest_results_consumer_fifo_ptr_array,
            EB_TRUE,
            RestResultsCreator,
            &restResultInitData,
            NULL);
    }

    // Entropy Coding Results
    {
        EntropyCodingResultsInitData entropyCodingResultInitData;

        EB_NEW(
            enc_handle_ptr->entropy_coding_results_resource_ptr,
            eb_system_resource_ctor,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count,
            EB_PacketizationProcessInitCount,
            &enc_handle_ptr->entropy_coding_results_producer_fifo_ptr_array,
            &enc_handle_ptr->entropy_coding_results_consumer_fifo_ptr_array,
            EB_TRUE,
            entropy_coding_results_creator,
            &entropyCodingResultInitData,
            NULL);
    }

    /************************************
    * App Callbacks
    ************************************/
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index)
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->app_callback_ptr = enc_handle_ptr->app_callback_ptr_array[instance_index];
    // svt Output Buffer Fifo Ptrs
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->stream_output_fifo_ptr     = (enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array[instance_index])[0];
        if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.recon_enabled)
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->recon_output_fifo_ptr      = (enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array[instance_index])[0];
    }

    /************************************
    * Contexts
    ************************************/

    // Resource Coordination Context
    EB_NEW(
        enc_handle_ptr->resource_coordination_context_ptr,
        resource_coordination_context_ctor,
        enc_handle_ptr->input_buffer_consumer_fifo_ptr_array[0],
        enc_handle_ptr->resource_coordination_results_producer_fifo_ptr_array[0],
        enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array[0],//ResourceCoordination works with ParentPCS
        enc_handle_ptr->sequence_control_set_instance_array,
        enc_handle_ptr->sequence_control_set_pool_producer_fifo_ptr_array[0],
        enc_handle_ptr->app_callback_ptr_array,
        enc_handle_ptr->compute_segments_total_count_array,
        enc_handle_ptr->encode_instance_total_count);

    // Picture Analysis Context
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->picture_analysis_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count; ++processIndex) {
        EbPictureBufferDescInitData  pictureBufferDescConf;
        pictureBufferDescConf.color_format = color_format;
        pictureBufferDescConf.max_width = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width;
        pictureBufferDescConf.max_height = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height;
        pictureBufferDescConf.bit_depth = EB_8BIT;
        pictureBufferDescConf.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
        pictureBufferDescConf.left_padding = 0;
        pictureBufferDescConf.right_padding = 0;
        pictureBufferDescConf.top_padding = 0;
        pictureBufferDescConf.bot_padding = 0;
        pictureBufferDescConf.split_mode = EB_FALSE;

        EB_NEW(
            enc_handle_ptr->picture_analysis_context_ptr_array[processIndex],
            picture_analysis_context_ctor,
            &pictureBufferDescConf,
            EB_TRUE,
            enc_handle_ptr->resource_coordination_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_analysis_results_producer_fifo_ptr_array[processIndex]);
   }

    // Picture Decision Context
    {
        // Initialize the various Picture types
        instance_index = 0;

        EB_NEW(
            enc_handle_ptr->picture_decision_context_ptr,
            picture_decision_context_ctor,
            enc_handle_ptr->picture_analysis_results_consumer_fifo_ptr_array[0],
            enc_handle_ptr->picture_decision_results_producer_fifo_ptr_array[0]);
    }

    // Motion Analysis Context
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->motion_estimation_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->motion_estimation_context_ptr_array[processIndex],
            motion_estimation_context_ctor,
            enc_handle_ptr->picture_decision_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->motion_estimation_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->nsq_present,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mrp_mode);
    }

    // Initial Rate Control Context
    EB_NEW(
        enc_handle_ptr->initial_rate_control_context_ptr,
        initial_rate_control_context_ctor,
        enc_handle_ptr->motion_estimation_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->initial_rate_control_results_producer_fifo_ptr_array[0]);
    // Source Based Operations Context
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->source_based_operations_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->source_based_operations_context_ptr_array[processIndex],
            source_based_operations_context_ctor,
            enc_handle_ptr->initial_rate_control_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr);
    }

    // Picture Manager Context
    EB_NEW(
        enc_handle_ptr->picture_manager_context_ptr,
        picture_manager_context_ctor,
        enc_handle_ptr->picture_demux_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER, 0)],
        enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array[0]);//The Child PCS Pool here
    // Rate Control Context
    EB_NEW(
        enc_handle_ptr->rate_control_context_ptr,
        rate_control_context_ctor,
        enc_handle_ptr->rate_control_tasks_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_results_producer_fifo_ptr_array[0],
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->intra_period_length);
    // Mode Decision Configuration Contexts
    {
        // Mode Decision Configuration Contexts
        EB_ALLOC_PTR_ARRAY(enc_handle_ptr->mode_decision_configuration_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count);

        for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count; ++processIndex) {
            EB_NEW(
                enc_handle_ptr->mode_decision_configuration_context_ptr_array[processIndex],
                mode_decision_configuration_context_ctor,
                enc_handle_ptr->rate_control_results_consumer_fifo_ptr_array[processIndex],
                enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array[EncDecPortLookup(ENCDEC_INPUT_PORT_MDC, processIndex)],
                ((enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64) *
                ((enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64));
        }
    }

    max_picture_width = 0;
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        if (max_picture_width < enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width)
            max_picture_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
    }

    // EncDec Contexts
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->enc_dec_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count; ++processIndex) {
#if PAL_SUP
        EB_NEW(
            enc_handle_ptr->enc_dec_context_ptr_array[processIndex],
            enc_dec_context_ctor,
            enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array[EncDecPortLookup(ENCDEC_INPUT_PORT_ENCDEC, processIndex)],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count +
                    //1 +
                    processIndex], // Add port lookup logic here JMJ
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.screen_content_mode,
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.enable_hbd_mode_decision,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
            );
#else
        EB_NEW(
            enc_handle_ptr->enc_dec_context_ptr_array[processIndex],
            enc_dec_context_ctor,
            enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array[EncDecPortLookup(ENCDEC_INPUT_PORT_ENCDEC, processIndex)],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count+
                //1 +
                    processIndex], // Add port lookup logic here JMJ
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.enable_hbd_mode_decision,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );
#endif
    }

    // Dlf Contexts
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->dlf_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->dlf_context_ptr_array[processIndex],
            dlf_context_ctor,
            enc_handle_ptr->enc_dec_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->dlf_results_producer_fifo_ptr_array[processIndex],             //output to EC
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );
    }

    //CDEF Contexts
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->cdef_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->cdef_context_ptr_array[processIndex],
            cdef_context_ctor,
            enc_handle_ptr->dlf_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->cdef_results_producer_fifo_ptr_array[processIndex],
            is16bit,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );
    }
    //Rest Contexts
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->rest_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->rest_context_ptr_array[processIndex],
            rest_context_ctor,
            enc_handle_ptr->cdef_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->rest_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[
                /*enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count*/ 1+ processIndex],
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );
    }

    // Entropy Coding Contexts
    EB_ALLOC_PTR_ARRAY(enc_handle_ptr->entropy_coding_context_ptr_array, enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count; ++processIndex) {
        EB_NEW(
            enc_handle_ptr->entropy_coding_context_ptr_array[processIndex],
            entropy_coding_context_ctor,
            enc_handle_ptr->rest_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->entropy_coding_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_ENTROPY_CODING, processIndex)],
            is16bit);
    }

    // Packetization Context
    EB_NEW(
        enc_handle_ptr->packetization_context_ptr,
        packetization_context_ctor,
        enc_handle_ptr->entropy_coding_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_PACKETIZATION, 0)]
        , enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count +
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count]
    );

    /************************************
    * Thread Handles
    ************************************/
    EbSvtAv1EncConfiguration   *config_ptr = &enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config;

    EbSetThreadManagementParameters(config_ptr);

    control_set_ptr = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr;

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
EB_API EbErrorType eb_deinit_encoder(EbComponentType *svt_enc_component){
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;
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
EB_API EbErrorType eb_init_handle(
    EbComponentType** p_handle,               // Function to be called in the future for manipulating the component
    void*              p_app_data,
    EbSvtAv1EncConfiguration  *config_ptr)              // pointer passed back to the client during callbacks

{
    EbErrorType           return_error = EB_ErrorNone;
    if(p_handle == NULL)
         return EB_ErrorBadParameter;

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
        eb_deinit_encoder(*p_handle);
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
* eb_deinit_handle
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_deinit_handle(
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
    SequenceControlSet       *sequence_control_set_ptr){
    int32_t intra_period               = 0;
    EbSvtAv1EncConfiguration   *config = &sequence_control_set_ptr->static_config;
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
void SetDefaultConfigurationParameters(
    SequenceControlSet       *sequence_control_set_ptr)
{
    // LCU Definitions
    sequence_control_set_ptr->sb_sz = MAX_SB_SIZE;
    sequence_control_set_ptr->max_sb_depth = (uint8_t)EB_MAX_LCU_DEPTH;

    // No Cropping Window
    sequence_control_set_ptr->conformance_window_flag = 0;
    sequence_control_set_ptr->cropping_left_offset = 0;
    sequence_control_set_ptr->cropping_right_offset = 0;
    sequence_control_set_ptr->cropping_top_offset = 0;
    sequence_control_set_ptr->cropping_bottom_offset = 0;

    sequence_control_set_ptr->static_config.enable_adaptive_quantization = 2;

    return;
}

static uint32_t compute_default_look_ahead(
    EbSvtAv1EncConfiguration*   config){
    int32_t lad = 0;
    if (config->rate_control_mode == 0)
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

void SetParamBasedOnInput(SequenceControlSet *sequence_control_set_ptr)
{
    uint16_t subsampling_x = sequence_control_set_ptr->subsampling_x;
    uint16_t subsampling_y = sequence_control_set_ptr->subsampling_y;
    sequence_control_set_ptr->general_frame_only_constraint_flag = 0;
    sequence_control_set_ptr->general_progressive_source_flag = 1;
    sequence_control_set_ptr->general_interlaced_source_flag = 0;

    // Update picture width, and picture height
    if (sequence_control_set_ptr->max_input_luma_width % MIN_BLOCK_SIZE) {
        sequence_control_set_ptr->max_input_pad_right = MIN_BLOCK_SIZE - (sequence_control_set_ptr->max_input_luma_width % MIN_BLOCK_SIZE);
        sequence_control_set_ptr->max_input_luma_width = sequence_control_set_ptr->max_input_luma_width + sequence_control_set_ptr->max_input_pad_right;
    } else {
        sequence_control_set_ptr->max_input_pad_right = 0;
    }

    if (sequence_control_set_ptr->max_input_luma_height % MIN_BLOCK_SIZE) {
        sequence_control_set_ptr->max_input_pad_bottom = MIN_BLOCK_SIZE - (sequence_control_set_ptr->max_input_luma_height % MIN_BLOCK_SIZE);
        sequence_control_set_ptr->max_input_luma_height = sequence_control_set_ptr->max_input_luma_height + sequence_control_set_ptr->max_input_pad_bottom;
    } else {
        sequence_control_set_ptr->max_input_pad_bottom = 0;
    }

    sequence_control_set_ptr->max_input_chroma_width = sequence_control_set_ptr->max_input_luma_width >> subsampling_x;
    sequence_control_set_ptr->max_input_chroma_height = sequence_control_set_ptr->max_input_luma_height >> subsampling_y;

    sequence_control_set_ptr->chroma_width = sequence_control_set_ptr->max_input_luma_width >> subsampling_x;
    sequence_control_set_ptr->chroma_height = sequence_control_set_ptr->max_input_luma_height >> subsampling_y;
    sequence_control_set_ptr->seq_header.max_frame_width = sequence_control_set_ptr->max_input_luma_width;
    sequence_control_set_ptr->seq_header.max_frame_height = sequence_control_set_ptr->max_input_luma_height;
    sequence_control_set_ptr->static_config.source_width = sequence_control_set_ptr->max_input_luma_width;
    sequence_control_set_ptr->static_config.source_height = sequence_control_set_ptr->max_input_luma_height;

    derive_input_resolution(
        sequence_control_set_ptr,
        sequence_control_set_ptr->seq_header.max_frame_width*sequence_control_set_ptr->seq_header.max_frame_height);
#if TWO_PASS
    // In two pass encoding, the first pass uses sb size=64
    if (sequence_control_set_ptr->static_config.screen_content_mode == 1 || sequence_control_set_ptr->use_output_stat_file)
#else
    if (sequence_control_set_ptr->static_config.screen_content_mode == 1)
#endif
        sequence_control_set_ptr->static_config.super_block_size = 64;
    else
        sequence_control_set_ptr->static_config.super_block_size = (sequence_control_set_ptr->static_config.enc_mode <= ENC_M3 && sequence_control_set_ptr->input_resolution >= INPUT_SIZE_1080i_RANGE) ? 128 : 64;

    sequence_control_set_ptr->static_config.super_block_size = (sequence_control_set_ptr->static_config.rate_control_mode > 1) ? 64 : sequence_control_set_ptr->static_config.super_block_size;
   // sequence_control_set_ptr->static_config.hierarchical_levels = (sequence_control_set_ptr->static_config.rate_control_mode > 1) ? 3 : sequence_control_set_ptr->static_config.hierarchical_levels;
    // Configure the padding
    sequence_control_set_ptr->left_padding = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->top_padding = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->right_padding = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->bot_padding = sequence_control_set_ptr->static_config.super_block_size + 4;
    sequence_control_set_ptr->static_config.enable_overlays = sequence_control_set_ptr->static_config.enable_altrefs == EB_FALSE ||
        (sequence_control_set_ptr->static_config.altref_nframes <= 1) ||
        (sequence_control_set_ptr->static_config.rate_control_mode > 0) ||
        sequence_control_set_ptr->static_config.encoder_bit_depth != EB_8BIT ?
        0 : sequence_control_set_ptr->static_config.enable_overlays;

    //0: MRP Mode 0 (4,3)
    //1: MRP Mode 1 (2,2)
    sequence_control_set_ptr->mrp_mode = (uint8_t) (sequence_control_set_ptr->static_config.enc_mode == ENC_M0) ? 0 : 1;

    //0: ON
    //1: OFF
    sequence_control_set_ptr->cdf_mode = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode <= ENC_M6) ? 0 : 1;

    //0: NSQ absent
    //1: NSQ present
    sequence_control_set_ptr->nsq_present = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode <= ENC_M5) ? 1 : 0;

    // Set down-sampling method     Settings
    // 0                            0: filtering
    // 1                            1: decimation
    if (sequence_control_set_ptr->static_config.enc_mode == ENC_M0)
        sequence_control_set_ptr->down_sampling_method_me_search = ME_FILTERED_DOWNSAMPLED;
    else
        sequence_control_set_ptr->down_sampling_method_me_search = ME_DECIMATED_DOWNSAMPLED;

    // Set over_boundary_block_mode     Settings
    // 0                            0: not allowed
    // 1                            1: allowed
    if (sequence_control_set_ptr->static_config.enc_mode == ENC_M0)
        sequence_control_set_ptr->over_boundary_block_mode = 1;
    else
        sequence_control_set_ptr->over_boundary_block_mode = 0;
#if M0_OPT
    sequence_control_set_ptr->mfmv_enabled = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode == ENC_M0 && sequence_control_set_ptr->static_config.screen_content_mode != 1) ? 1 : 0;
#else
    sequence_control_set_ptr->mfmv_enabled = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode == ENC_M0) ? 1 : 0;
#endif

    // Set hbd_mode_decision OFF for high encode modes or bitdepth < 10
    if (sequence_control_set_ptr->static_config.enc_mode > ENC_M0 || sequence_control_set_ptr->static_config.encoder_bit_depth < 10)
        sequence_control_set_ptr->static_config.enable_hbd_mode_decision = 0;
}

void CopyApiFromApp(
    SequenceControlSet       *sequence_control_set_ptr,
    EbSvtAv1EncConfiguration   *pComponentParameterStructure){
    uint32_t                  hmeRegionIndex = 0;

    sequence_control_set_ptr->max_input_luma_width = pComponentParameterStructure->source_width;
    sequence_control_set_ptr->max_input_luma_height = pComponentParameterStructure->source_height;

    sequence_control_set_ptr->frame_rate = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->frame_rate;

    sequence_control_set_ptr->general_frame_only_constraint_flag = 0;
    sequence_control_set_ptr->general_progressive_source_flag = 1;
    sequence_control_set_ptr->general_interlaced_source_flag = 0;

    // SB Definitions
    sequence_control_set_ptr->static_config.pred_structure = 2; // Hardcoded(Cleanup)
    sequence_control_set_ptr->static_config.enable_qp_scaling_flag = 1;

    sequence_control_set_ptr->max_cu_size = (uint8_t)64;
    sequence_control_set_ptr->min_cu_size = (uint8_t)8;
    sequence_control_set_ptr->max_intra_size = (uint8_t)32;
    sequence_control_set_ptr->min_intra_size = (uint8_t)8;
    sequence_control_set_ptr->intra4x4_flag = 1;
    sequence_control_set_ptr->max_ref_count = 1;

    // Cropping Definitions - Hardcoded(CleanUp)
    sequence_control_set_ptr->cropping_left_offset = -1;
    sequence_control_set_ptr->cropping_right_offset = -1;
    sequence_control_set_ptr->cropping_top_offset = -1;
    sequence_control_set_ptr->cropping_bottom_offset = -1;

    // Padding Offsets
    sequence_control_set_ptr->sb_sz = (uint8_t)((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->sb_sz;
    sequence_control_set_ptr->max_sb_depth = (uint8_t)((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->partition_depth;

    if (sequence_control_set_ptr->cropping_left_offset == -1 &&
        sequence_control_set_ptr->cropping_right_offset == -1 &&
        sequence_control_set_ptr->cropping_top_offset == -1 &&
        sequence_control_set_ptr->cropping_bottom_offset == -1) {
        sequence_control_set_ptr->conformance_window_flag = 0;
    }
    else
        sequence_control_set_ptr->conformance_window_flag = 1;
    if (sequence_control_set_ptr->cropping_left_offset == -1)
        sequence_control_set_ptr->cropping_left_offset = 0;
    if (sequence_control_set_ptr->cropping_right_offset == -1)
        sequence_control_set_ptr->cropping_right_offset = 0;
    if (sequence_control_set_ptr->cropping_top_offset == -1)
        sequence_control_set_ptr->cropping_top_offset = 0;
    if (sequence_control_set_ptr->cropping_bottom_offset == -1)
        sequence_control_set_ptr->cropping_bottom_offset = 0;

    sequence_control_set_ptr->static_config.intra_period_length = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->intra_period_length;
    sequence_control_set_ptr->static_config.intra_refresh_type = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->intra_refresh_type;
    sequence_control_set_ptr->static_config.base_layer_switch_mode = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->base_layer_switch_mode;
    sequence_control_set_ptr->static_config.hierarchical_levels = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hierarchical_levels;
    sequence_control_set_ptr->static_config.enc_mode = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enc_mode;
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    sequence_control_set_ptr->static_config.snd_pass_enc_mode = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->snd_pass_enc_mode;
#endif
    sequence_control_set_ptr->intra_period_length = sequence_control_set_ptr->static_config.intra_period_length;
    sequence_control_set_ptr->intra_refresh_type = sequence_control_set_ptr->static_config.intra_refresh_type;
    sequence_control_set_ptr->max_temporal_layers = sequence_control_set_ptr->static_config.hierarchical_levels;
    sequence_control_set_ptr->static_config.use_qp_file = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->use_qp_file;
#if TWO_PASS
    sequence_control_set_ptr->static_config.input_stat_file = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->input_stat_file;
    sequence_control_set_ptr->static_config.output_stat_file = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->output_stat_file;
    sequence_control_set_ptr->use_input_stat_file = sequence_control_set_ptr->static_config.input_stat_file ? 1 : 0;
    sequence_control_set_ptr->use_output_stat_file = sequence_control_set_ptr->static_config.output_stat_file ? 1 : 0;
#endif
    // Deblock Filter
    sequence_control_set_ptr->static_config.disable_dlf_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->disable_dlf_flag;

    // Local Warped Motion
    sequence_control_set_ptr->static_config.enable_warped_motion = EB_TRUE;

    // Global motion
    sequence_control_set_ptr->static_config.enable_global_motion = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_global_motion;

    // OBMC
    sequence_control_set_ptr->static_config.enable_obmc = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_obmc;

    // RDOQ
    sequence_control_set_ptr->static_config.enable_rdoq = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_rdoq;

    // Filter intra prediction
    sequence_control_set_ptr->static_config.enable_filter_intra = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_filter_intra;

    // ME Tools
    sequence_control_set_ptr->static_config.use_default_me_hme = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->use_default_me_hme;
    sequence_control_set_ptr->static_config.enable_hme_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_hme_flag;
    sequence_control_set_ptr->static_config.enable_hme_level0_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_hme_level0_flag;
    sequence_control_set_ptr->static_config.enable_hme_level1_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_hme_level1_flag;
    sequence_control_set_ptr->static_config.enable_hme_level2_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_hme_level2_flag;
    sequence_control_set_ptr->static_config.search_area_width = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->search_area_width;
    sequence_control_set_ptr->static_config.search_area_height = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->search_area_height;
    sequence_control_set_ptr->static_config.number_hme_search_region_in_width = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->number_hme_search_region_in_width;
    sequence_control_set_ptr->static_config.number_hme_search_region_in_height = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->number_hme_search_region_in_height;
    sequence_control_set_ptr->static_config.hme_level0_total_search_area_width = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level0_total_search_area_width;
    sequence_control_set_ptr->static_config.hme_level0_total_search_area_height = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level0_total_search_area_height;
    sequence_control_set_ptr->static_config.ext_block_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->ext_block_flag;
    sequence_control_set_ptr->static_config.in_loop_me_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->in_loop_me_flag;

    for (hmeRegionIndex = 0; hmeRegionIndex < sequence_control_set_ptr->static_config.number_hme_search_region_in_width; ++hmeRegionIndex) {
        sequence_control_set_ptr->static_config.hme_level0_search_area_in_width_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level0_search_area_in_width_array[hmeRegionIndex];
        sequence_control_set_ptr->static_config.hme_level1_search_area_in_width_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level1_search_area_in_width_array[hmeRegionIndex];
        sequence_control_set_ptr->static_config.hme_level2_search_area_in_width_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level2_search_area_in_width_array[hmeRegionIndex];
    }

    for (hmeRegionIndex = 0; hmeRegionIndex < sequence_control_set_ptr->static_config.number_hme_search_region_in_height; ++hmeRegionIndex) {
        sequence_control_set_ptr->static_config.hme_level0_search_area_in_height_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level0_search_area_in_height_array[hmeRegionIndex];
        sequence_control_set_ptr->static_config.hme_level1_search_area_in_height_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level1_search_area_in_height_array[hmeRegionIndex];
        sequence_control_set_ptr->static_config.hme_level2_search_area_in_height_array[hmeRegionIndex] = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->hme_level2_search_area_in_height_array[hmeRegionIndex];
    }

    //Denoise - Hardcoded(CleanUp)
    sequence_control_set_ptr->static_config.enable_denoise_flag = 0;

    //Film Grain
    sequence_control_set_ptr->static_config.film_grain_denoise_strength = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->film_grain_denoise_strength;
    sequence_control_set_ptr->film_grain_denoise_strength = sequence_control_set_ptr->static_config.film_grain_denoise_strength;

    // MD Parameters
    sequence_control_set_ptr->static_config.enable_hbd_mode_decision = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->encoder_bit_depth > 8 ? ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_hbd_mode_decision : 0;
    sequence_control_set_ptr->static_config.constrained_intra = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->constrained_intra;
    sequence_control_set_ptr->static_config.enable_palette = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_palette;
    sequence_control_set_ptr->static_config.olpd_refinement = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->olpd_refinement;
    // Adaptive Loop Filter
    sequence_control_set_ptr->static_config.tile_rows = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->tile_rows;
    sequence_control_set_ptr->static_config.tile_columns = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->tile_columns;
    sequence_control_set_ptr->static_config.unrestricted_motion_vector = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->unrestricted_motion_vector;

    // Rate Control
    sequence_control_set_ptr->static_config.scene_change_detection = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->scene_change_detection;
    sequence_control_set_ptr->static_config.rate_control_mode = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->rate_control_mode;
    sequence_control_set_ptr->static_config.look_ahead_distance = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->look_ahead_distance;
    sequence_control_set_ptr->static_config.frames_to_be_encoded = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->frames_to_be_encoded;
    sequence_control_set_ptr->static_config.frame_rate = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->frame_rate;
    sequence_control_set_ptr->static_config.frame_rate_denominator = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->frame_rate_denominator;
    sequence_control_set_ptr->static_config.frame_rate_numerator = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->frame_rate_numerator;

    sequence_control_set_ptr->static_config.target_bit_rate = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->target_bit_rate;

    sequence_control_set_ptr->static_config.max_qp_allowed = (sequence_control_set_ptr->static_config.rate_control_mode) ?
        ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->max_qp_allowed :
        63;

    sequence_control_set_ptr->static_config.min_qp_allowed = (sequence_control_set_ptr->static_config.rate_control_mode) ?
        ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->min_qp_allowed :
        1; // lossless coding not supported

    //Segmentation
    //TODO: check RC mode and set only when RC is enabled in the final version.
    sequence_control_set_ptr->static_config.enable_adaptive_quantization = pComponentParameterStructure->enable_adaptive_quantization;

    // Misc
    sequence_control_set_ptr->static_config.encoder_bit_depth = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->encoder_bit_depth;
    sequence_control_set_ptr->static_config.encoder_color_format = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->encoder_color_format;
    if (sequence_control_set_ptr->static_config.encoder_color_format == EB_YUV400) {
        SVT_LOG("SVT [Warning]: Color format EB_YUV400 not supported, set to EB_YUV420\n");
        sequence_control_set_ptr->static_config.encoder_color_format = EB_YUV420;
    }
    sequence_control_set_ptr->chroma_format_idc = (uint32_t)(sequence_control_set_ptr->static_config.encoder_color_format);
    sequence_control_set_ptr->encoder_bit_depth = (uint32_t)(sequence_control_set_ptr->static_config.encoder_bit_depth);

    sequence_control_set_ptr->subsampling_x = (sequence_control_set_ptr->chroma_format_idc == EB_YUV444 ? 1 : 2) - 1;
    sequence_control_set_ptr->subsampling_y = (sequence_control_set_ptr->chroma_format_idc >= EB_YUV422 ? 1 : 2) - 1;
    sequence_control_set_ptr->static_config.ten_bit_format = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->ten_bit_format;
    sequence_control_set_ptr->static_config.compressed_ten_bit_format = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->compressed_ten_bit_format;

    // Thresholds
    sequence_control_set_ptr->static_config.high_dynamic_range_input = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->high_dynamic_range_input;
    sequence_control_set_ptr->static_config.screen_content_mode = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->screen_content_mode;

    // Annex A parameters
    sequence_control_set_ptr->static_config.profile = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->profile;
    sequence_control_set_ptr->static_config.tier = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->tier;
    sequence_control_set_ptr->static_config.level = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->level;
    sequence_control_set_ptr->static_config.stat_report = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->stat_report;

    sequence_control_set_ptr->static_config.injector_frame_rate = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->injector_frame_rate;
    sequence_control_set_ptr->static_config.speed_control_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->speed_control_flag;

    // Buffers - Hardcoded(Cleanup)
    sequence_control_set_ptr->static_config.use_cpu_flags = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->use_cpu_flags;

    sequence_control_set_ptr->static_config.channel_id = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->channel_id;
    sequence_control_set_ptr->static_config.active_channel_count = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->active_channel_count;
    sequence_control_set_ptr->static_config.logical_processors = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->logical_processors;
    sequence_control_set_ptr->static_config.target_socket = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->target_socket;
    sequence_control_set_ptr->qp = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->qp;
    sequence_control_set_ptr->static_config.recon_enabled = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->recon_enabled;

    // Extract frame rate from Numerator and Denominator if not 0
    if (sequence_control_set_ptr->static_config.frame_rate_numerator != 0 && sequence_control_set_ptr->static_config.frame_rate_denominator != 0)
        sequence_control_set_ptr->frame_rate = sequence_control_set_ptr->static_config.frame_rate = (((sequence_control_set_ptr->static_config.frame_rate_numerator << 8) / (sequence_control_set_ptr->static_config.frame_rate_denominator)) << 8);
    // Get Default Intra Period if not specified
    if (sequence_control_set_ptr->static_config.intra_period_length == -2)
        sequence_control_set_ptr->intra_period_length = sequence_control_set_ptr->static_config.intra_period_length = compute_default_intra_period(sequence_control_set_ptr);
    if (sequence_control_set_ptr->static_config.look_ahead_distance == (uint32_t)~0)
        sequence_control_set_ptr->static_config.look_ahead_distance = compute_default_look_ahead(&sequence_control_set_ptr->static_config);
    else
        sequence_control_set_ptr->static_config.look_ahead_distance = cap_look_ahead_distance(&sequence_control_set_ptr->static_config);

    sequence_control_set_ptr->static_config.enable_altrefs = pComponentParameterStructure->enable_altrefs;
    sequence_control_set_ptr->static_config.altref_strength = pComponentParameterStructure->altref_strength;
    sequence_control_set_ptr->static_config.altref_nframes = pComponentParameterStructure->altref_nframes;
    sequence_control_set_ptr->static_config.enable_overlays = pComponentParameterStructure->enable_overlays;

    sequence_control_set_ptr->static_config.sq_weight = pComponentParameterStructure->sq_weight;
    sequence_control_set_ptr->static_config.enable_auto_max_partition = pComponentParameterStructure->enable_auto_max_partition;

    sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th = pComponentParameterStructure->md_stage_1_cand_prune_th;
    sequence_control_set_ptr->static_config.md_stage_1_class_prune_th = pComponentParameterStructure->md_stage_1_class_prune_th;
    sequence_control_set_ptr->static_config.md_stage_2_cand_prune_th = pComponentParameterStructure->md_stage_2_cand_prune_th;
    sequence_control_set_ptr->static_config.md_stage_2_class_prune_th = pComponentParameterStructure->md_stage_2_class_prune_th;

    return;
}

/******************************************
* Verify Settings
******************************************/
#define PowerOfTwoCheck(x) (((x) != 0) && (((x) & (~(x) + 1)) == (x)))

static int VerifyHmeDimention(unsigned int index, unsigned int HmeLevel0SearchAreaInWidth, uint32_t number_hme_search_region_in_width_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT], unsigned int number_hme_search_region_in_width)
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
static int VerifyHmeDimentionL1L2(unsigned int index, uint32_t number_hme_search_region_in_width_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT], unsigned int number_hme_search_region_in_width)
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

static EbErrorType VerifySettings(
    SequenceControlSet       *sequence_control_set_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbSvtAv1EncConfiguration *config = &sequence_control_set_ptr->static_config;
    unsigned int channelNumber = config->channel_id;
    if (config->enc_mode > MAX_ENC_PRESET) {
        SVT_LOG("Error instance %u: EncoderMode must be in the range of [0-%d]\n", channelNumber + 1, MAX_ENC_PRESET);
        return_error = EB_ErrorBadParameter;
    }
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    if (config->snd_pass_enc_mode > MAX_ENC_PRESET + 1) {
        SVT_LOG("Error instance %u: Second pass encoder mode must be in the range of [0-%d]\n", channelNumber + 1, MAX_ENC_PRESET + 1);
        return_error = EB_ErrorBadParameter;
    }
#endif
    if (config->ext_block_flag > 1) {
        SVT_LOG("Error instance %u: ExtBlockFlag must be [0-1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->in_loop_me_flag > 1) {
        SVT_LOG("Error instance %u: InLoopMeFlag must be [0-1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (sequence_control_set_ptr->max_input_luma_width < 64) {
        SVT_LOG("Error instance %u: Source Width must be at least 64\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (sequence_control_set_ptr->max_input_luma_height < 64) {
        SVT_LOG("Error instance %u: Source Width must be at least 64\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->pred_structure != 2) {
        SVT_LOG("Error instance %u: Pred Structure must be [2]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->base_layer_switch_mode == 1 && config->pred_structure != 2) {
        SVT_LOG("Error Instance %u: Base Layer Switch Mode 1 only when Prediction Structure is Random Access\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (sequence_control_set_ptr->max_input_luma_width % 8 && sequence_control_set_ptr->static_config.compressed_ten_bit_format == 1) {
        SVT_LOG("Error Instance %u: Only multiple of 8 width is supported for compressed 10-bit inputs \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (sequence_control_set_ptr->max_input_luma_width % 2) {
        SVT_LOG("Error Instance %u: Source Width must be even for YUV_420 colorspace\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (sequence_control_set_ptr->max_input_luma_height % 2) {
        SVT_LOG("Error Instance %u: Source Height must be even for YUV_420 colorspace\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (sequence_control_set_ptr->max_input_luma_width > 4096) {
        SVT_LOG("Error instance %u: Source Width must be less than 4096\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (sequence_control_set_ptr->max_input_luma_height > 2160) {
        SVT_LOG("Error instance %u: Source Height must be less than 2160\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->qp > MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: QP must be [0 - %d]\n", channelNumber + 1, MAX_QP_VALUE);
        return_error = EB_ErrorBadParameter;
    }
    if (config->hierarchical_levels != 3 && config->hierarchical_levels != 4) {
        SVT_LOG("Error instance %u: Hierarchical Levels supported [3-4]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->intra_period_length < -2 || config->intra_period_length > 255) {
        SVT_LOG("Error Instance %u: The intra period must be [-2 - 255] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->intra_refresh_type > 2 || config->intra_refresh_type < 1) {
        SVT_LOG("Error Instance %u: Invalid intra Refresh Type [1-2]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->base_layer_switch_mode > 1) {
        SVT_LOG("Error Instance %u: Invalid Base Layer Switch Mode [0-1] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->disable_dlf_flag > 1) {
        SVT_LOG("Error Instance %u: Invalid LoopFilterDisable. LoopFilterDisable must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_default_me_hme > 1) {
        SVT_LOG("Error Instance %u: invalid use_default_me_hme. use_default_me_hme must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->enable_hme_flag > 1) {
        SVT_LOG("Error Instance %u: invalid HME. HME must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level0_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel0. HMELevel0 must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level1_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel1. HMELevel1 must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_level2_flag > 1) {
        SVT_LOG("Error Instance %u: invalid enable HMELevel2. HMELevel2 must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->search_area_width > 256) || (config->search_area_width == 0)) {
        SVT_LOG("Error Instance %u: Invalid search_area_width. search_area_width must be [1 - 256]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->search_area_height > 256) || (config->search_area_height == 0)) {
        SVT_LOG("Error Instance %u: Invalid search_area_height. search_area_height must be [1 - 256]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->enable_hme_flag) {
        if ((config->number_hme_search_region_in_width > (uint32_t)EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT) || (config->number_hme_search_region_in_width == 0)) {
            SVT_LOG("Error Instance %u: Invalid number_hme_search_region_in_width. number_hme_search_region_in_width must be [1 - %d]\n", channelNumber + 1, EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT);
            return_error = EB_ErrorBadParameter;
        }

        if ((config->number_hme_search_region_in_height > (uint32_t)EB_HME_SEARCH_AREA_ROW_MAX_COUNT) || (config->number_hme_search_region_in_height == 0)) {
            SVT_LOG("Error Instance %u: Invalid number_hme_search_region_in_height. number_hme_search_region_in_height must be [1 - %d]\n", channelNumber + 1, EB_HME_SEARCH_AREA_ROW_MAX_COUNT);
            return_error = EB_ErrorBadParameter;
        }

        if ((config->hme_level0_total_search_area_height > 256) || (config->hme_level0_total_search_area_height == 0)) {
            SVT_LOG("Error Instance %u: Invalid hme_level0_total_search_area_height. hme_level0_total_search_area_height must be [1 - 256]\n", channelNumber + 1);
            return_error = EB_ErrorBadParameter;
        }
        if ((config->hme_level0_total_search_area_width > 256) || (config->hme_level0_total_search_area_width == 0)) {
            SVT_LOG("Error Instance %u: Invalid hme_level0_total_search_area_width. hme_level0_total_search_area_width must be [1 - 256]\n", channelNumber + 1);
            return_error = EB_ErrorBadParameter;
        }
        if (VerifyHmeDimention(channelNumber + 1, config->hme_level0_total_search_area_height, config->hme_level0_search_area_in_height_array, config->number_hme_search_region_in_height))
            return_error = EB_ErrorBadParameter;
        if (VerifyHmeDimention(channelNumber + 1, config->hme_level0_total_search_area_width, config->hme_level0_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level1_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level1_search_area_in_height_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level2_search_area_in_width_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level2_search_area_in_height_array, config->number_hme_search_region_in_width))
            return_error = EB_ErrorBadParameter;
    }

    if (config->profile > 2) {
        SVT_LOG("Error Instance %u: The maximum allowed profile value is 2 \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    // Check if the current input video is conformant with the Level constraint
    if (config->frame_rate > (240 << 16)) {
        SVT_LOG("Error Instance %u: The maximum allowed frame rate is 240 fps\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    // Check that the frame_rate is non-zero
    if (config->frame_rate <= 0) {
        SVT_LOG("Error Instance %u: The frame rate should be greater than 0 fps \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->intra_period_length < -2 || config->intra_period_length > 255) {
        SVT_LOG("Error Instance %u: The intra period must be [-2 - 255] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->constrained_intra > 1) {
        SVT_LOG("Error Instance %u: The constrained intra must be [0 - 1] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->rate_control_mode > 2) {

        SVT_LOG("Error Instance %u: The rate control mode must be [0 - 2] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if ((config->rate_control_mode == 1 || config->rate_control_mode == 2) && config->look_ahead_distance != (uint32_t)config->intra_period_length) {
        SVT_LOG("Error Instance %u: The rate control mode 1/2 LAD must be equal to intra_period \n", channelNumber + 1);

        return_error = EB_ErrorBadParameter;
    }
    if (config->look_ahead_distance > MAX_LAD && config->look_ahead_distance != (uint32_t)~0) {
        SVT_LOG("Error Instance %u: The lookahead distance must be [0 - %d] \n", channelNumber + 1, MAX_LAD);

        return_error = EB_ErrorBadParameter;
    }
    if (config->tile_rows < 0 || config->tile_columns < 0 || config->tile_rows > 6 || config->tile_columns > 6) {
        SVT_LOG("Error Instance %u: Log2Tile rows/cols must be [0 - 6] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->unrestricted_motion_vector > 1) {
        SVT_LOG("Error Instance %u : Invalid Unrestricted Motion Vector flag [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->scene_change_detection > 1) {
        SVT_LOG("Error Instance %u: The scene change detection must be [0 - 1] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->max_qp_allowed > MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: MaxQpAllowed must be [0 - %d]\n", channelNumber + 1, MAX_QP_VALUE);
        return_error = EB_ErrorBadParameter;
    }
    else if (config->min_qp_allowed >= MAX_QP_VALUE) {
        SVT_LOG("Error instance %u: MinQpAllowed must be [0 - %d]\n", channelNumber + 1, MAX_QP_VALUE-1);
        return_error = EB_ErrorBadParameter;
    }
    else if ((config->min_qp_allowed) > (config->max_qp_allowed)) {
        SVT_LOG("Error Instance %u:  MinQpAllowed must be smaller than MaxQpAllowed\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->stat_report > 1) {
        SVT_LOG("Error instance %u : Invalid StatReport. StatReport must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->high_dynamic_range_input > 1) {
        SVT_LOG("Error instance %u : Invalid HighDynamicRangeInput. HighDynamicRangeInput must be [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->screen_content_mode > 2) {
        SVT_LOG("Error instance %u : Invalid screen_content_mode. screen_content_mode must be [0 - 2]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (sequence_control_set_ptr->static_config.enable_adaptive_quantization > 2) {
        SVT_LOG("Error instance %u : Invalid enable_adaptive_quantization. enable_adaptive_quantization must be [0-2]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if ((config->encoder_bit_depth != 8) &&
        (config->encoder_bit_depth != 10)
        ) {
        SVT_LOG("Error instance %u: Encoder Bit Depth shall be only 8 or 10 \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    // Check if the EncoderBitDepth is conformant with the Profile constraint
    if ((config->profile == 0 || config->profile == 1) && config->encoder_bit_depth > 10) {
        SVT_LOG("Error instance %u: The encoder bit depth shall be equal to 8 or 10 for Main/High Profile\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 0 && config->encoder_color_format > EB_YUV420) {
        SVT_LOG("Error instance %u: Non 420 color format requires profile 1 or 2\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 1 && config->encoder_color_format != EB_YUV444) {
        SVT_LOG("Error instance %u: Profile 1 requires 4:4:4 color format\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->profile == 2 && config->encoder_bit_depth <= 10 && config->encoder_color_format != EB_YUV422) {
        SVT_LOG("Error instance %u: Profile 2 bit-depth < 10 requires 4:2:2 color format\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->compressed_ten_bit_format !=0)
    {
        SVT_LOG("Error instance %u: Compressed ten bit format is not supported in this version \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->speed_control_flag > 1) {
        SVT_LOG("Error Instance %u: Invalid Speed Control flag [0 - 1]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_cpu_flags & CPU_FLAGS_INVALID) {
        SVT_LOG("Error Instance %u: param '-asm' have invalid value.\n"
            "Value should be [0 - 11, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512, max]\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->target_socket != -1 && config->target_socket != 0 && config->target_socket != 1) {
        SVT_LOG("Error instance %u: Invalid target_socket. target_socket must be [-1 - 1] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    // alt-ref frames related
    if (config->altref_strength > ALTREF_MAX_STRENGTH ) {
        SVT_LOG("Error instance %u: invalid altref-strength, should be in the range [0 - %d] \n", channelNumber + 1, ALTREF_MAX_STRENGTH);
        return_error = EB_ErrorBadParameter;
    }

    if (config->altref_nframes > ALTREF_MAX_NFRAMES ) {
        SVT_LOG("Error instance %u: invalid altref-nframes, should be in the range [0 - %d] \n", channelNumber + 1, ALTREF_MAX_NFRAMES);
        return_error = EB_ErrorBadParameter;
    }

    // palette
    if (config->enable_palette < (int32_t)(-1) || config->enable_palette >6) {
        SVT_LOG( "Error instance %u: Invalid Palette Mode [0 .. 6], your input: %i\n", channelNumber + 1, config->enable_palette);
        return_error = EB_ErrorBadParameter;
    }

    // RDOQ
    if (config->enable_rdoq != (int8_t)0 && config->enable_rdoq != (int8_t)1 && config->enable_rdoq != (int8_t)-1) {
        SVT_LOG( "Error instance %u: Invalid RDOQ parameter [-1, 0, 1], your input: %i\n", channelNumber + 1, config->enable_rdoq);
        return_error = EB_ErrorBadParameter;
    }

    // mdc refinement
    if (config->olpd_refinement < (int32_t)(-1) || config->olpd_refinement > 1) {
        SVT_LOG("Error instance %u: Invalid OLPD Refinement Mode [0 .. 1], your input: %i\n", channelNumber + 1, config->olpd_refinement);
        return_error = EB_ErrorBadParameter;
    }
    else if (config->olpd_refinement == 1 && config->enc_mode >= ENC_M1) {
        SVT_LOG("Error instance %u: Invalid OLPD Refinement mode for M%d [0], your input: %i\n", channelNumber + 1, config->enc_mode, config->olpd_refinement);
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
    config_ptr->frames_to_be_encoded = 0;
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
    config_ptr->base_layer_switch_mode = 0;
    config_ptr->enc_mode = MAX_ENC_PRESET;
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    config_ptr->snd_pass_enc_mode = MAX_ENC_PRESET + 1;
#endif
    config_ptr->intra_period_length = -2;
    config_ptr->intra_refresh_type = 1;
    config_ptr->hierarchical_levels = 4;
    config_ptr->pred_structure = EB_PRED_RANDOM_ACCESS;
    config_ptr->disable_dlf_flag = EB_FALSE;
    config_ptr->enable_warped_motion = EB_TRUE;
    config_ptr->enable_global_motion = EB_TRUE;
    config_ptr->enable_obmc = EB_TRUE;
    config_ptr->enable_rdoq = AUTO_MODE;
    config_ptr->enable_filter_intra = EB_TRUE;
    config_ptr->in_loop_me_flag = EB_TRUE;
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
    config_ptr->constrained_intra = EB_FALSE;
    config_ptr->enable_palette = -1;
    config_ptr->olpd_refinement = -1;
    // Bitstream options
    //config_ptr->codeVpsSpsPps = 0;
    //config_ptr->codeEosNal = 0;
    config_ptr->unrestricted_motion_vector = EB_TRUE;

    config_ptr->high_dynamic_range_input = 0;
    config_ptr->screen_content_mode = 2;

    // Annex A parameters
    config_ptr->profile = 0;
    config_ptr->tier = 0;
    config_ptr->level = 0;

    // Latency
    config_ptr->injector_frame_rate = 60 << 16;
    config_ptr->speed_control_flag = 0;
    config_ptr->super_block_size = 128;

    config_ptr->sb_sz = 64;
    config_ptr->partition_depth = (uint8_t)EB_MAX_LCU_DEPTH;
    //config_ptr->latency_mode = 0;
    config_ptr->speed_control_flag = 0;
    config_ptr->film_grain_denoise_strength = 0;

    // CPU Flags
    config_ptr->use_cpu_flags = CPU_FLAGS_ALL;

    // Channel info
    config_ptr->logical_processors = 0;
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

    config_ptr->sq_weight = 100;

    config_ptr->md_stage_1_cand_prune_th = 75;
    config_ptr->md_stage_1_class_prune_th = 100;
    config_ptr->md_stage_2_cand_prune_th = 15;
    config_ptr->md_stage_2_class_prune_th = 25;

    config_ptr->enable_auto_max_partition = 1;    return return_error;
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
    SVT_LOG("\nSVT [config]: HierarchicalLevels / BaseLayerSwitchMode / PredStructure\t\t: %d / %d / %d ", config->hierarchical_levels, config->base_layer_switch_mode, config->pred_structure);
    if (config->rate_control_mode == 1)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate (kbps)/ LookaheadDistance / SceneChange\t\t: VBR / %d / %d / %d ", (int)config->target_bit_rate/1000, config->look_ahead_distance, config->scene_change_detection);
    else if (config->rate_control_mode == 2)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate (kbps)/ LookaheadDistance / SceneChange\t\t: Constraint VBR / %d / %d / %d ", (int)config->target_bit_rate/1000, config->look_ahead_distance, config->scene_change_detection);
    else
        SVT_LOG("\nSVT [config]: BRC Mode / QP  / LookaheadDistance / SceneChange\t\t\t: CQP / %d / %d / %d ", scs->qp, config->look_ahead_distance, config->scene_change_detection);
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
EB_API EbErrorType eb_svt_enc_set_parameter(
    EbComponentType              *svt_enc_component,
    EbSvtAv1EncConfiguration     *pComponentParameterStructure)
{
    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;

    EbErrorType           return_error  = EB_ErrorNone;
    EbEncHandle        *pEncCompData  = (EbEncHandle*)svt_enc_component->p_component_private;
    uint32_t              instance_index = 0;

    // Acquire Config Mutex
    eb_block_on_mutex(pEncCompData->sequence_control_set_instance_array[instance_index]->config_mutex);

    SetDefaultConfigurationParameters(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    CopyApiFromApp(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr,
        (EbSvtAv1EncConfiguration*)pComponentParameterStructure);

    return_error = (EbErrorType)VerifySettings(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    if (return_error == EB_ErrorBadParameter)
        return EB_ErrorBadParameter;
    SetParamBasedOnInput(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    // Initialize the Prediction Structure Group
    EB_NO_THROW_NEW(
        pEncCompData->sequence_control_set_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        prediction_structure_group_ctor,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enc_mode,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.base_layer_switch_mode);
    if (!pEncCompData->sequence_control_set_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr) {
        eb_release_mutex(pEncCompData->sequence_control_set_instance_array[instance_index]->config_mutex);
        return EB_ErrorInsufficientResources;
    }
    // Set the Prediction Structure
    pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pred_struct_ptr = get_prediction_structure(
        pEncCompData->sequence_control_set_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.pred_structure,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_ref_count,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_temporal_layers);

    return_error = load_default_buffer_configuration_settings(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    print_lib_params(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    // Release Config Mutex
    eb_release_mutex(pEncCompData->sequence_control_set_instance_array[instance_index]->config_mutex);

    return return_error;
}
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_enc_stream_header(
    EbComponentType           *svt_enc_component,
    EbBufferHeaderType        **output_stream_ptr)
{
    EbErrorType              return_error = EB_ErrorNone;
    EbEncHandle             *pEncCompData  = (EbEncHandle*)svt_enc_component->p_component_private;
    SequenceControlSet      *sequenceControlSetPtr = pEncCompData->sequence_control_set_instance_array[0]->sequence_control_set_ptr;
    Bitstream                bitstream;
    OutputBitstreamUnit      output_bitstream;
    EbBufferHeaderType      *outputStreamBuffer;

    memset(&bitstream, 0, sizeof(Bitstream));
    memset(&output_bitstream, 0, sizeof(OutputBitstreamUnit));
    bitstream.output_bitstream_ptr = &output_bitstream;

    outputStreamBuffer = (EbBufferHeaderType *)malloc(sizeof(EbBufferHeaderType));
    if (!outputStreamBuffer) {
        return EB_ErrorInsufficientResources;
    }

    outputStreamBuffer->p_buffer = (uint8_t *)malloc(sizeof(uint8_t) * PACKETIZATION_PROCESS_BUFFER_SIZE);
    if (!outputStreamBuffer->p_buffer) {
        free(outputStreamBuffer);
        return EB_ErrorInsufficientResources;
    }

    outputStreamBuffer->size = sizeof(EbBufferHeaderType);
    outputStreamBuffer->n_alloc_len = PACKETIZATION_PROCESS_BUFFER_SIZE;
    outputStreamBuffer->p_app_private = NULL;
    outputStreamBuffer->pic_type = EB_AV1_INVALID_PICTURE;
    outputStreamBuffer->n_filled_len = 0;

    ((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_begin_av1 = outputStreamBuffer->p_buffer;

    output_bitstream_reset(bitstream.output_bitstream_ptr);

    // Code the SPS
    encode_sps_av1(&bitstream, sequenceControlSetPtr);

    outputStreamBuffer->n_filled_len = (uint32_t)(((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_av1 - ((OutputBitstreamUnit *)bitstream.output_bitstream_ptr)->buffer_begin_av1);

    *output_stream_ptr = outputStreamBuffer;

    return return_error;
}
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_release_enc_stream_header(
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
EB_API EbErrorType eb_svt_enc_eos_nal(
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
static EbErrorType CopyFrameBuffer(
    SequenceControlSet            *sequence_control_set_ptr,
    uint8_t                          *dst,
    uint8_t                          *src)
{
    EbSvtAv1EncConfiguration          *config = &sequence_control_set_ptr->static_config;
    EbErrorType                      return_error = EB_ErrorNone;

    EbPictureBufferDesc           *input_picture_ptr = (EbPictureBufferDesc*)dst;
    EbSvtIOFormat                   *inputPtr = (EbSvtIOFormat*)src;
    uint16_t                         inputRowIndex;
    EbBool                           is16BitInput = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is16BitInput) {
        uint32_t     lumaBufferOffset = (input_picture_ptr->stride_y*sequence_control_set_ptr->top_padding + sequence_control_set_ptr->left_padding) << is16BitInput;
        uint32_t     chromaBufferOffset = (input_picture_ptr->stride_cr*(sequence_control_set_ptr->top_padding >> 1) + (sequence_control_set_ptr->left_padding >> 1)) << is16BitInput;
        uint16_t     lumaStride = input_picture_ptr->stride_y << is16BitInput;
        uint16_t     chromaStride = input_picture_ptr->stride_cb << is16BitInput;
        uint16_t     lumaWidth = (uint16_t)(input_picture_ptr->width - sequence_control_set_ptr->max_input_pad_right) << is16BitInput;
        uint16_t     chromaWidth = (lumaWidth >> 1) << is16BitInput;
        uint16_t     lumaHeight = (uint16_t)(input_picture_ptr->height - sequence_control_set_ptr->max_input_pad_bottom);

        uint16_t     sourceLumaStride = (uint16_t)(inputPtr->y_stride);
        uint16_t     sourceCrStride = (uint16_t)(inputPtr->cr_stride);
        uint16_t     sourceCbStride = (uint16_t)(inputPtr->cb_stride);

        //uint16_t     lumaHeight  = input_picture_ptr->max_height;
        // Y
        for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
            EB_MEMCPY((input_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                (inputPtr->luma + sourceLumaStride * inputRowIndex),
                lumaWidth);
        }

        // U
        for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            EB_MEMCPY((input_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                (inputPtr->cb + (sourceCbStride*inputRowIndex)),
                chromaWidth);
        }

        // V
        for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            EB_MEMCPY((input_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                (inputPtr->cr + (sourceCrStride*inputRowIndex)),
                chromaWidth);
        }
    }
    else if (is16BitInput && config->compressed_ten_bit_format == 1)
    {
        {
            uint32_t  lumaBufferOffset = (input_picture_ptr->stride_y*sequence_control_set_ptr->top_padding + sequence_control_set_ptr->left_padding);
            uint32_t  chromaBufferOffset = (input_picture_ptr->stride_cr*(sequence_control_set_ptr->top_padding >> 1) + (sequence_control_set_ptr->left_padding >> 1));
            uint16_t  lumaStride = input_picture_ptr->stride_y;
            uint16_t  chromaStride = input_picture_ptr->stride_cb;
            uint16_t  lumaWidth = (uint16_t)(input_picture_ptr->width - sequence_control_set_ptr->max_input_pad_right);
            uint16_t  chromaWidth = (lumaWidth >> 1);
            uint16_t  lumaHeight = (uint16_t)(input_picture_ptr->height - sequence_control_set_ptr->max_input_pad_bottom);

            uint16_t  sourceLumaStride = (uint16_t)(inputPtr->y_stride);
            uint16_t  sourceCrStride = (uint16_t)(inputPtr->cr_stride);
            uint16_t  sourceCbStride = (uint16_t)(inputPtr->cb_stride);

            // Y 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
                EB_MEMCPY((input_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                    (inputPtr->luma + sourceLumaStride * inputRowIndex),
                    lumaWidth);
            }

            // U 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                EB_MEMCPY((input_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                    (inputPtr->cb + (sourceCbStride*inputRowIndex)),
                    chromaWidth);
            }

            // V 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                EB_MEMCPY((input_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                    (inputPtr->cr + (sourceCrStride*inputRowIndex)),
                    chromaWidth);
            }

            //efficient copy - final
            //compressed 2Bit in 1D format
            {
                uint16_t luma2BitWidth = sequence_control_set_ptr->max_input_luma_width / 4;
                uint16_t lumaHeight = sequence_control_set_ptr->max_input_luma_height;

                uint16_t sourceLuma2BitStride = sourceLumaStride / 4;
                uint16_t sourceChroma2BitStride = sourceLuma2BitStride >> 1;

                for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
                    EB_MEMCPY(input_picture_ptr->buffer_bit_inc_y + luma2BitWidth * inputRowIndex, inputPtr->luma_ext + sourceLuma2BitStride * inputRowIndex, luma2BitWidth);
                }
                for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                    EB_MEMCPY(input_picture_ptr->buffer_bit_inc_cb + (luma2BitWidth >> 1)*inputRowIndex, inputPtr->cb_ext + sourceChroma2BitStride * inputRowIndex, luma2BitWidth >> 1);
                }
                for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                    EB_MEMCPY(input_picture_ptr->buffer_bit_inc_cr + (luma2BitWidth >> 1)*inputRowIndex, inputPtr->cr_ext + sourceChroma2BitStride * inputRowIndex, luma2BitWidth >> 1);
                }
            }
        }
    }
    else { // 10bit packed

        uint32_t lumaOffset = 0, chromaOffset = 0;
        uint32_t lumaBufferOffset = (input_picture_ptr->stride_y*sequence_control_set_ptr->top_padding + sequence_control_set_ptr->left_padding);
        uint32_t chromaBufferOffset = (input_picture_ptr->stride_cr*(sequence_control_set_ptr->top_padding >> 1) + (sequence_control_set_ptr->left_padding >> 1));
        uint16_t lumaWidth = (uint16_t)(input_picture_ptr->width - sequence_control_set_ptr->max_input_pad_right);
        uint16_t chromaWidth = (lumaWidth >> 1);
        uint16_t lumaHeight = (uint16_t)(input_picture_ptr->height - sequence_control_set_ptr->max_input_pad_bottom);

        uint16_t sourceLumaStride = (uint16_t)(inputPtr->y_stride);
        uint16_t sourceCrStride = (uint16_t)(inputPtr->cr_stride);
        uint16_t sourceCbStride = (uint16_t)(inputPtr->cb_stride);

        un_pack2d(
            (uint16_t*)(inputPtr->luma + lumaOffset),
            sourceLumaStride,
            input_picture_ptr->buffer_y + lumaBufferOffset,
            input_picture_ptr->stride_y,
            input_picture_ptr->buffer_bit_inc_y + lumaBufferOffset,
            input_picture_ptr->stride_bit_inc_y,
            lumaWidth,
            lumaHeight);

        un_pack2d(
            (uint16_t*)(inputPtr->cb + chromaOffset),
            sourceCbStride,
            input_picture_ptr->buffer_cb + chromaBufferOffset,
            input_picture_ptr->stride_cb,
            input_picture_ptr->buffer_bit_inc_cb + chromaBufferOffset,
            input_picture_ptr->stride_bit_inc_cb,
            chromaWidth,
            (lumaHeight >> 1));

        un_pack2d(
            (uint16_t*)(inputPtr->cr + chromaOffset),
            sourceCrStride,
            input_picture_ptr->buffer_cr + chromaBufferOffset,
            input_picture_ptr->stride_cr,
            input_picture_ptr->buffer_bit_inc_cr + chromaBufferOffset,
            input_picture_ptr->stride_bit_inc_cr,
            chromaWidth,
            (lumaHeight >> 1));
    }
    return return_error;
}
static void CopyInputBuffer(
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
        CopyFrameBuffer(sequenceControlSet, dst->p_buffer, src->p_buffer);
}

/**********************************
* Empty This Buffer
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_enc_send_picture(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbEncHandle          *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *ebWrapperPtr;

    // Take the buffer and put it into our internal queue structure
    eb_get_empty_object(
        enc_handle_ptr->input_buffer_producer_fifo_ptr_array[0],
        &ebWrapperPtr);

    if (p_buffer != NULL) {
        CopyInputBuffer(
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr,
            (EbBufferHeaderType*)ebWrapperPtr->object_ptr,
            p_buffer);
    }

    eb_post_full_object(ebWrapperPtr);

    return EB_ErrorNone;
}
static void CopyOutputReconBuffer(
    EbBufferHeaderType   *dst,
    EbBufferHeaderType   *src
)
{
    // copy output bitstream fileds
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
        EB_MEMCPY(dst->p_buffer, src->p_buffer, src->n_filled_len);

    return;
}

/**********************************
* eb_svt_get_packet sends out packet
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_get_packet(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType  **p_buffer,
    unsigned char          pic_send_done)
{
    EbErrorType             return_error = EB_ErrorNone;
    EbEncHandle          *pEncCompData = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *ebWrapperPtr = NULL;
    EbBufferHeaderType    *packet;
    if (pic_send_done)
        eb_get_full_object(
        (pEncCompData->output_stream_buffer_consumer_fifo_ptr_dbl_array[0])[0],
            &ebWrapperPtr);
    else
        eb_get_full_object_non_blocking(
        (pEncCompData->output_stream_buffer_consumer_fifo_ptr_dbl_array[0])[0],
            &ebWrapperPtr);

    if (ebWrapperPtr) {
        packet = (EbBufferHeaderType*)ebWrapperPtr->object_ptr;
        if ( packet->flags & 0xfffffff0 )
            return_error = EB_ErrorMax;
        // return the output stream buffer
        *p_buffer = packet;

        // save the wrapper pointer for the release
        (*p_buffer)->wrapper_ptr = (void*)ebWrapperPtr;
    }
    else
        return_error = EB_NoErrorEmptyQueue;
    return return_error;
}

#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API void eb_svt_release_out_buffer(
    EbBufferHeaderType  **p_buffer)
{
    if (p_buffer && (*p_buffer)->wrapper_ptr)
        // Release out put buffer back into the pool
        eb_release_object((EbObjectWrapper  *)(*p_buffer)->wrapper_ptr);
    return;
}

/**********************************
* Fill This Buffer
**********************************/
#ifdef __GNUC__
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_get_recon(
    EbComponentType      *svt_enc_component,
    EbBufferHeaderType   *p_buffer)
{
    EbErrorType           return_error = EB_ErrorNone;
    EbEncHandle          *pEncCompData = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *ebWrapperPtr = NULL;

    if (pEncCompData->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.recon_enabled) {
        eb_get_full_object_non_blocking(
            (pEncCompData->output_recon_buffer_consumer_fifo_ptr_dbl_array[0])[0],
            &ebWrapperPtr);

        if (ebWrapperPtr) {
            EbBufferHeaderType* objPtr = (EbBufferHeaderType*)ebWrapperPtr->object_ptr;
            CopyOutputReconBuffer(
                p_buffer,
                objPtr);

            if (p_buffer->flags != EB_BUFFERFLAG_EOS && p_buffer->flags != 0)
                return_error = EB_ErrorMax;
            eb_release_object((EbObjectWrapper  *)ebWrapperPtr);
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
    EbEncHandle          *pEncCompData = (EbEncHandle*)svt_enc_component->p_component_private;
    EbObjectWrapper      *ebWrapperPtr = NULL;
    EbBufferHeaderType    *outputPacket;

    eb_get_empty_object(
        (pEncCompData->output_stream_buffer_producer_fifo_ptr_dbl_array[0])[0],
        &ebWrapperPtr);

    outputPacket            = (EbBufferHeaderType*)ebWrapperPtr->object_ptr;

    outputPacket->size     = 0;
    outputPacket->flags    = error_code;
    outputPacket->p_buffer   = NULL;

    eb_post_full_object(ebWrapperPtr);
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

    printf("SVT [version]:\tSVT-AV1 Encoder Lib v%d.%d.%d\n", SVT_VERSION_MAJOR, SVT_VERSION_MINOR, SVT_VERSION_PATCHLEVEL);
#if ( defined( _MSC_VER ) && (_MSC_VER < 1910) )
    printf("SVT [build]  : Visual Studio 2013");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1910) )
    printf("SVT [build]  :\tVisual Studio 2017");
#elif defined(__GNUC__)
    printf("SVT [build]  :\tGCC %s\t", __VERSION__);
#else
    printf("SVT [build]  :\tunknown compiler");
#endif
    printf(" %u bit\n", (unsigned) sizeof(void*) * 8);
    printf("LIB Build date: %s %s\n", __DATE__, __TIME__);
    printf("-------------------------------------------\n");

    SwitchToRealTime();

    // Set Component Size & Version
    svt_enc_component->size = sizeof(EbComponentType);

    EB_NEW(handle, eb_enc_handle_ctor, svt_enc_component);
    svt_enc_component->p_component_private = handle;

    return return_error;
}
static EbErrorType allocate_frame_buffer(
    SequenceControlSet       *sequence_control_set_ptr,
    EbBufferHeaderType        *inputBuffer)
{
    EbErrorType   return_error = EB_ErrorNone;
    EbPictureBufferDescInitData input_picture_buffer_desc_init_data;
    EbSvtAv1EncConfiguration   * config = &sequence_control_set_ptr->static_config;
    uint8_t is16bit = config->encoder_bit_depth > 8 ? 1 : 0;
    // Init Picture Init data
    input_picture_buffer_desc_init_data.max_width = (uint16_t)sequence_control_set_ptr->max_input_luma_width;
    input_picture_buffer_desc_init_data.max_height = (uint16_t)sequence_control_set_ptr->max_input_luma_height;
    input_picture_buffer_desc_init_data.bit_depth = (EbBitDepthEnum)config->encoder_bit_depth;
    input_picture_buffer_desc_init_data.color_format = (EbColorFormat)config->encoder_color_format;

    if (config->compressed_ten_bit_format == 1)
        input_picture_buffer_desc_init_data.buffer_enable_mask = 0;
    else
        input_picture_buffer_desc_init_data.buffer_enable_mask = is16bit ? PICTURE_BUFFER_DESC_FULL_MASK : 0;
    input_picture_buffer_desc_init_data.left_padding = sequence_control_set_ptr->left_padding;
    input_picture_buffer_desc_init_data.right_padding = sequence_control_set_ptr->right_padding;
    input_picture_buffer_desc_init_data.top_padding = sequence_control_set_ptr->top_padding;
    input_picture_buffer_desc_init_data.bot_padding = sequence_control_set_ptr->bot_padding;

    input_picture_buffer_desc_init_data.split_mode = is16bit ? EB_TRUE : EB_FALSE;

    input_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    if (is16bit && config->compressed_ten_bit_format == 1)
        input_picture_buffer_desc_init_data.split_mode = EB_FALSE;  //do special allocation for 2bit data down below.
    // Enhanced Picture Buffer
    {
        EbPictureBufferDesc* buf;
        EB_NEW(
            buf,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&input_picture_buffer_desc_init_data);
        inputBuffer->p_buffer = (uint8_t*)buf;
        if (is16bit && config->compressed_ten_bit_format == 1) {
            //pack 4 2bit pixels into 1Byte
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_y, (input_picture_buffer_desc_init_data.max_width / 4)*(input_picture_buffer_desc_init_data.max_height));
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_cb, (input_picture_buffer_desc_init_data.max_width / 8)*(input_picture_buffer_desc_init_data.max_height / 2));
            EB_MALLOC_ALIGNED_ARRAY(buf->buffer_bit_inc_cr, (input_picture_buffer_desc_init_data.max_width / 8)*(input_picture_buffer_desc_init_data.max_height / 2));
        }
    }

    return return_error;
}
/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbInputBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr)
{
    EbBufferHeaderType* inputBuffer;
    SequenceControlSet        *sequence_control_set_ptr = (SequenceControlSet*)objectInitDataPtr;

    *objectDblPtr = NULL;
    EB_CALLOC(inputBuffer, 1, sizeof(EbBufferHeaderType));
    *objectDblPtr = (EbPtr)inputBuffer;
    // Initialize Header
    inputBuffer->size = sizeof(EbBufferHeaderType);

    allocate_frame_buffer(
        sequence_control_set_ptr,
        inputBuffer);

    inputBuffer->p_app_private = NULL;

    return EB_ErrorNone;
}

void EbInputBufferHeaderDestoryer(    EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EbPictureBufferDesc* buf = (EbPictureBufferDesc*)obj->p_buffer;
    EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_y);
    EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cb);
    EB_FREE_ALIGNED_ARRAY(buf->buffer_bit_inc_cr);

    EB_DELETE(buf);
    EB_FREE(obj);
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbOutputBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr objectInitDataPtr)
{
    EbSvtAv1EncConfiguration   * config = (EbSvtAv1EncConfiguration*)objectInitDataPtr;
    uint32_t n_stride = (uint32_t)(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(config->source_width * config->source_height));  //TBC
    EbBufferHeaderType* outBufPtr;

    *objectDblPtr = NULL;
    EB_CALLOC(outBufPtr, 1, sizeof(EbBufferHeaderType));
    *objectDblPtr = (EbPtr)outBufPtr;

    // Initialize Header
    outBufPtr->size = sizeof(EbBufferHeaderType);

    EB_MALLOC(outBufPtr->p_buffer, n_stride);

    outBufPtr->n_alloc_len = n_stride;
    outBufPtr->p_app_private = NULL;

    (void)objectInitDataPtr;

    return EB_ErrorNone;
}

void EbOutputBufferHeaderDestoryer(    EbPtr p)
{
    EbBufferHeaderType* obj = (EbBufferHeaderType*)p;
    EB_FREE(obj->p_buffer);
    EB_FREE(obj);
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbOutputReconBufferHeaderCreator(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr)
{
    EbBufferHeaderType         *recon_buffer;
    SequenceControlSet        *sequence_control_set_ptr = (SequenceControlSet*)objectInitDataPtr;
    const uint32_t luma_size =
        sequence_control_set_ptr->seq_header.max_frame_width    *
        sequence_control_set_ptr->seq_header.max_frame_height;
    // both u and v
    const uint32_t chroma_size = luma_size >> 1;
    const uint32_t tenBit = (sequence_control_set_ptr->static_config.encoder_bit_depth > 8);
    const uint32_t frameSize = (luma_size + chroma_size) << tenBit;

    *objectDblPtr = NULL;
    EB_CALLOC(recon_buffer, 1, sizeof(EbBufferHeaderType));
    *objectDblPtr = (EbPtr)recon_buffer;

    // Initialize Header
    recon_buffer->size = sizeof(EbBufferHeaderType);

    // Assign the variables
    EB_MALLOC(recon_buffer->p_buffer, frameSize);

    recon_buffer->n_alloc_len = frameSize;
    recon_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

void EbOutputReconBufferHeaderDestoryer(    EbPtr p)
{
    EbBufferHeaderType *obj = (EbBufferHeaderType*)p;
    EB_FREE(obj->p_buffer);
    EB_FREE(obj);
}
