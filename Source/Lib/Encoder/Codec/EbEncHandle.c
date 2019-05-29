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
#include "EbRateControlResults.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbPredictionStructure.h"
#include "EbDlfProcess.h"
#include "EbCdefProcess.h"
#include "EbRestProcess.h"


#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#if defined(__linux__) || defined(__APPLE__)
#include <pthread.h>
#include <errno.h>
#endif


#define RTCD_C
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

EbMemoryMapEntry                 *memory_map;
uint32_t                         *memory_map_index;
uint64_t                         *total_lib_memory;

uint32_t                         lib_malloc_count = 0;
uint32_t                         lib_thread_count = 0;
uint32_t                         lib_semaphore_count = 0;
uint32_t                         lib_mutex_count = 0;

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
#define MAX_PROCESSOR_GROUP 16
processorGroup                   lp_group[MAX_PROCESSOR_GROUP];
#endif

/**************************************
* Instruction Set Support
**************************************/

#if defined(_MSC_VER)
# include <intrin.h>
#endif
// Helper Functions
void RunCpuid(uint32_t eax, uint32_t ecx, int32_t* abcd)
{
#if defined(_MSC_VER)
    __cpuidex(abcd, eax, ecx);
#else
    uint32_t ebx = 0, edx = 0;
# if defined( __i386__ ) && defined ( __PIC__ )
    /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__("movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
    __asm__("cpuid" : "+b" (ebx),
# endif
        "+a" (eax), "+c" (ecx), "=d" (edx));
    abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#endif
}
int32_t CheckXcr0Ymm()
{
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}
int32_t Check4thGenIntelCoreFeatures()
{
    int32_t abcd[4];
    int32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    int32_t avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);

    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   &&
       CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
       CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    RunCpuid(1, 0, abcd);
    if ((abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask)
        return 0;

    if (!CheckXcr0Ymm())
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
    RunCpuid(7, 0, abcd);
    if ((abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask)
        return 0;
    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    RunCpuid(0x80000001, 0, abcd);
    if ((abcd[2] & (1 << 5)) == 0)
        return 0;
    return 1;
}
static int32_t CanUseIntelCore4thGenFeatures()
{
    static int32_t the_4th_gen_features_available = -1;
    /* test is performed once */
    if (the_4th_gen_features_available < 0)
        the_4th_gen_features_available = Check4thGenIntelCoreFeatures();
    return the_4th_gen_features_available;
}
EbAsm GetCpuAsmType()
{
    EbAsm asm_type = ASM_NON_AVX2;

    if (CanUseIntelCore4thGenFeatures() == 1)
        asm_type = ASM_AVX2;
    else
        // Need to change to support lower CPU Technologies
        asm_type = ASM_NON_AVX2;
    return asm_type;
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
    if (processor_id_len < 0 || processor_id_len >= 128)
        return EB_ErrorInsufficientResources;
    if (physical_id_len < 0 || physical_id_len >= 128)
        return EB_ErrorInsufficientResources;
    memset(lp_group, 0, sizeof(lp_group));

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
                if (socket_id < 0 || socket_id > 15) {
                    fclose(fin);
                    return EB_ErrorInsufficientResources;
                }
                if (socket_id + 1 > num_groups)
                    num_groups = socket_id + 1;
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

void EbSetThreadManagementParameters(EbSvtAv1EncConfiguration   *config_ptr) {
    uint32_t num_logical_processors = GetNumProcessors();
#ifdef _WIN32
    // For system with a single processor group(no more than 64 logic processors all together)
    // Affinity of the thread can be set to one or more logical processors
    if (num_groups == 1) {
        uint32_t lps = config_ptr->logical_processors == 0 ? num_logical_processors :
            config_ptr->logical_processors < num_logical_processors ? config_ptr->logical_processors : num_logical_processors;
        group_affinity.Mask = GetAffinityMask(lps);
    }
    else if (num_groups > 1) { // For system with multiple processor group
        if (config_ptr->logical_processors == 0) {
            if (config_ptr->target_socket != -1) {
                group_affinity.Group = config_ptr->target_socket;
            }
        }
        else {
            uint32_t num_lp_per_group = num_logical_processors / num_groups;
            if (config_ptr->target_socket == -1) {
                if (config_ptr->logical_processors > num_lp_per_group) {
                    alternate_groups = EB_TRUE;
                    SVT_LOG("SVT [WARNING]: -lp(logical processors) setting is ignored. Run on both sockets. \n");
                }
                else {
                    group_affinity.Mask = GetAffinityMask(config_ptr->logical_processors);
                }
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
#endif
}
void asmSetConvolveAsmTable(void);
void asmSetConvolveHbdAsmTable(void);
void init_intra_dc_predictors_c_internal(void);
void init_intra_predictors_internal(void);
void av1_init_me_luts(void);

void SwitchToRealTime(){
#if defined(__linux__) || defined(__APPLE__)

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
#if MINI_GOP_PCS
        uint32_t min_ppcs_count = (1 << config->hierarchical_levels) + 1; // min picture count to start encoding
#else
        uint32_t min_ppcs_count = (2 << config->hierarchical_levels) + 1; // min picture count to start encoding
#endif
        fps        = fps > 120 ? 120   : fps;
        fps        = fps < 24  ? 24    : fps;

#if MINI_GOP_PCS
        if (core_count == 4) 
            return min_ppcs_count;
        else
            ppcs_count = MAX(min_ppcs_count, fps);
#else
        ppcs_count = MAX(min_ppcs_count, fps);
#endif
#if NEW_BUFF_CFG        
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
#else
        ppcs_count = ((ppcs_count * 4) >> 1);  // 2 sec worth of internal buffering
#endif
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
    if (sequence_control_set_ptr->static_config.logical_processors != 0)
        core_count = sequence_control_set_ptr->static_config.logical_processors < core_count ?
            sequence_control_set_ptr->static_config.logical_processors: core_count;
#endif

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
#if CHECK_MEM_REDUCTION
    core_count = 4;
#endif
    int32_t return_ppcs = set_parent_pcs(&sequence_control_set_ptr->static_config, 
                    core_count, sequence_control_set_ptr->input_resolution);
    if (return_ppcs == -1)
        return EB_ErrorInsufficientResources;
    uint32_t input_pic = (uint32_t)return_ppcs;
#if BUG_FIX_LOOKAHEAD
    sequence_control_set_ptr->input_buffer_fifo_init_count = input_pic + SCD_LAD + sequence_control_set_ptr->static_config.look_ahead_distance;
#else
    sequence_control_set_ptr->input_buffer_fifo_init_count         =
        input_pic + SCD_LAD + sequence_control_set_ptr->static_config.look_ahead_distance ;
#endif
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

#if ALTREF_FILTERING_SUPPORT
	sequence_control_set_ptr->tf_segment_column_count = meSegW;//1;//
	sequence_control_set_ptr->tf_segment_row_count =  meSegH;//1;//
#endif
    //#====================== Data Structures and Picture Buffers ======================
#if BUG_FIX_LOOKAHEAD
    sequence_control_set_ptr->picture_control_set_pool_init_count       = input_pic + SCD_LAD + sequence_control_set_ptr->static_config.look_ahead_distance;
#if ALT_REF_OVERLAY
    if (sequence_control_set_ptr->static_config.enable_overlays)
        sequence_control_set_ptr->picture_control_set_pool_init_count = MAX(sequence_control_set_ptr->picture_control_set_pool_init_count,
            sequence_control_set_ptr->static_config.look_ahead_distance + // frames in the LAD
            sequence_control_set_ptr->static_config.look_ahead_distance / (1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 1 +  // number of overlayes in the LAD 
            ((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + SCD_LAD) * 2 +// minigop formation in PD + SCD_LAD *(normal pictures + potential pictures ) 
            (1 << sequence_control_set_ptr->static_config.hierarchical_levels)); // minigop in PM
#endif
#else
    sequence_control_set_ptr->picture_control_set_pool_init_count       = input_pic + sequence_control_set_ptr->static_config.look_ahead_distance + SCD_LAD;
#endif
    sequence_control_set_ptr->picture_control_set_pool_init_count_child = MAX(MAX(MIN(3, core_count/2), core_count / 6), 1);
    sequence_control_set_ptr->reference_picture_buffer_init_count       = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          sequence_control_set_ptr->static_config.look_ahead_distance + SCD_LAD;
    sequence_control_set_ptr->pa_reference_picture_buffer_init_count    = MAX((uint32_t)(input_pic >> 1),
                                                                          (uint32_t)((1 << sequence_control_set_ptr->static_config.hierarchical_levels) + 2)) +
                                                                          sequence_control_set_ptr->static_config.look_ahead_distance + SCD_LAD;
    sequence_control_set_ptr->output_recon_buffer_fifo_init_count       = sequence_control_set_ptr->reference_picture_buffer_init_count;
#if ALT_REF_OVERLAY
    sequence_control_set_ptr->overlay_input_picture_buffer_init_count   = sequence_control_set_ptr->static_config.enable_overlays ? 
                                                                          (2 << sequence_control_set_ptr->static_config.hierarchical_levels) + SCD_LAD : 1;
#endif

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
#if NEW_BUFF_CFG
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
#else
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->picture_analysis_process_init_count = MAX(MIN(15, core_count), core_count / 6));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->motion_estimation_process_init_count = MAX(MIN(20, core_count), core_count / 3));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->source_based_operations_process_init_count = MAX(MIN(3, core_count), core_count / 12));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->mode_decision_configuration_process_init_count = MAX(MIN(3, core_count), core_count / 12));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->enc_dec_process_init_count = MAX(MIN(40, core_count), core_count));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->entropy_coding_process_init_count = MAX(MIN(3, core_count), core_count / 12));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->dlf_process_init_count = MAX(MIN(40, core_count), core_count));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->cdef_process_init_count = MAX(MIN(40, core_count), core_count));
    sequence_control_set_ptr->total_process_init_count += (sequence_control_set_ptr->rest_process_init_count = MAX(MIN(40, core_count), core_count));
#endif

    sequence_control_set_ptr->total_process_init_count += 6; // single processes count
    printf("Number of logical cores available: %u\nNumber of PPCS %u\n", core_count, sequence_control_set_ptr->picture_control_set_pool_init_count);

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

    while ((type != rateControlPorts[portIndex].type) && (type != RATE_CONTROL_INPUT_PORT_INVALID)) {
        portCount += rateControlPorts[portIndex++].count;
    }

    return (portCount + portTypeIndex);
}
// Rate Control
static uint32_t RateControlPortTotalCount(void){
    uint32_t portIndex = 0;
    uint32_t total_count = 0;

    while (rateControlPorts[portIndex].type != RATE_CONTROL_INPUT_PORT_INVALID) {
        total_count += rateControlPorts[portIndex++].count;
    }

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

    while ((type != encDecPorts[portIndex].type) && (type != ENCDEC_INPUT_PORT_INVALID)) {
        portCount += encDecPorts[portIndex++].count;
    }

    return (portCount + portTypeIndex);
}
// EncDec
static uint32_t EncDecPortTotalCount(void){
    uint32_t portIndex = 0;
    uint32_t total_count = 0;

    while (encDecPorts[portIndex].type != ENCDEC_INPUT_PORT_INVALID) {
        total_count += encDecPorts[portIndex++].count;
    }

    return total_count;
}
/*****************************************
 * Input Port Total Count
 *****************************************/

void lib_svt_encoder_send_error_exit(
    EbPtr                    hComponent,
    uint32_t                 error_code);

/**********************************
* Encoder Library Handle Constructor
**********************************/
static EbErrorType eb_enc_handle_ctor(
    EbEncHandle **encHandleDblPtr,
    EbComponentType * ebHandlePtr)
{
    EbErrorType return_error = EB_ErrorNone;
#if !MEM_MAP_OPT
    uint32_t instance_index  = 0;
#endif
    // Allocate Memory
    EbEncHandle *enc_handle_ptr = (EbEncHandle*)malloc(sizeof(EbEncHandle));

    *encHandleDblPtr          = enc_handle_ptr;
    if (enc_handle_ptr == (EbEncHandle*)EB_NULL)
        return EB_ErrorInsufficientResources;
#if MEM_MAP_OPT
    enc_handle_ptr->memory_map                = (EbMemoryMapEntry*)malloc(sizeof(EbMemoryMapEntry));
    enc_handle_ptr->memory_map_index          = 0;
    enc_handle_ptr->total_lib_memory          = sizeof(EbComponentType) + sizeof(EbEncHandle) + sizeof(EbMemoryMapEntry);
    enc_handle_ptr->memory_map_init_address   = enc_handle_ptr->memory_map;
#else
    enc_handle_ptr->memory_map = (EbMemoryMapEntry*)malloc(sizeof(EbMemoryMapEntry) * MAX_NUM_PTR);
    enc_handle_ptr->memory_map_index = 0;
    enc_handle_ptr->total_lib_memory = sizeof(EbEncHandle) + sizeof(EbMemoryMapEntry) * MAX_NUM_PTR;
#endif
    // Save Memory Map Pointers
    total_lib_memory                        = &enc_handle_ptr->total_lib_memory;
    memory_map                              = enc_handle_ptr->memory_map;
    memory_map_index                        = &enc_handle_ptr->memory_map_index;
    lib_malloc_count                        = 0;
    lib_thread_count                        = 0;
    lib_mutex_count                         = 0;
    lib_semaphore_count                     = 0;

    if (memory_map == (EbMemoryMapEntry*)EB_NULL)
        return EB_ErrorInsufficientResources;

    return_error = InitThreadManagmentParams();

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
#if MEM_MAP_OPT
    enc_handle_ptr->memory_map->prev_entry                                = EB_NULL;
    enc_handle_ptr->encode_instance_total_count                           = EB_EncodeInstancesTotalCount;
    enc_handle_ptr->compute_segments_total_count_array                    = EB_ComputeSegmentInitCount;
#else
    enc_handle_ptr->encode_instance_total_count = EB_EncodeInstancesTotalCount;

    EB_MALLOC(uint32_t*, enc_handle_ptr->compute_segments_total_count_array, sizeof(uint32_t) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        enc_handle_ptr->compute_segments_total_count_array[instance_index] = EB_ComputeSegmentInitCount;
    }
#endif
    // Config Set Count
    enc_handle_ptr->sequence_control_set_pool_total_count                 = EB_SequenceControlSetPoolInitCount;

    // Sequence Control Set Buffers
    enc_handle_ptr->sequence_control_set_pool_ptr                         = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->sequence_control_set_pool_producer_fifo_ptr_array     = (EbFifo**)EB_NULL;

    // Picture Buffers
    enc_handle_ptr->reference_picture_pool_ptr_array                      = (EbSystemResource**)EB_NULL;
    enc_handle_ptr->pa_reference_picture_pool_ptr_array                   = (EbSystemResource**)EB_NULL;

#if ALT_REF_OVERLAY
    // Overlay input picture
    enc_handle_ptr->overlay_input_picture_pool_ptr_array = (EbSystemResource**)EB_NULL;
    enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array = (EbFifo***)EB_NULL;
#endif
    // Picture Buffer Producer Fifos
    enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array    = (EbFifo***)EB_NULL;
    enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array = (EbFifo***)EB_NULL;

    // Threads
    enc_handle_ptr->resource_coordination_thread_handle             = (EbHandle)EB_NULL;
    enc_handle_ptr->picture_analysis_thread_handle_array            = (EbHandle*)EB_NULL;
    enc_handle_ptr->picture_decision_thread_handle                  = (EbHandle)EB_NULL;
    enc_handle_ptr->motion_estimation_thread_handle_array           = (EbHandle*)EB_NULL;
    enc_handle_ptr->initial_rate_control_thread_handle              = (EbHandle)EB_NULL;
    enc_handle_ptr->source_based_operations_thread_handle_array     = (EbHandle*)EB_NULL;
    enc_handle_ptr->picture_manager_thread_handle                   = (EbHandle)EB_NULL;
    enc_handle_ptr->rate_control_thread_handle                      = (EbHandle)EB_NULL;
    enc_handle_ptr->mode_decision_configuration_thread_handle_array = (EbHandle*)EB_NULL;
    enc_handle_ptr->enc_dec_thread_handle_array                     = (EbHandle*)EB_NULL;
    enc_handle_ptr->entropy_coding_thread_handle_array              = (EbHandle*)EB_NULL;
    enc_handle_ptr->packetization_thread_handle                     = (EbHandle)EB_NULL;
    enc_handle_ptr->dlf_thread_handle_array                         = (EbHandle*)EB_NULL;
    enc_handle_ptr->cdef_thread_handle_array                        = (EbHandle*)EB_NULL;
    enc_handle_ptr->rest_thread_handle_array                        = (EbHandle*)EB_NULL;

    // Contexts
    enc_handle_ptr->resource_coordination_context_ptr             = (EbPtr)EB_NULL;
    enc_handle_ptr->picture_analysis_context_ptr_array            = (EbPtr*)EB_NULL;
    enc_handle_ptr->picture_decision_context_ptr                  = (EbPtr)EB_NULL;
    enc_handle_ptr->motion_estimation_context_ptr_array           = (EbPtr*)EB_NULL;
    enc_handle_ptr->initial_rate_control_context_ptr              = (EbPtr)EB_NULL;
    enc_handle_ptr->source_based_operations_context_ptr_array     = (EbPtr*)EB_NULL;
    enc_handle_ptr->picture_manager_context_ptr                   = (EbPtr)EB_NULL;
    enc_handle_ptr->rate_control_context_ptr                      = (EbPtr)EB_NULL;
    enc_handle_ptr->mode_decision_configuration_context_ptr_array = (EbPtr*)EB_NULL;
    enc_handle_ptr->enc_dec_context_ptr_array                     = (EbPtr*)EB_NULL;
    enc_handle_ptr->entropy_coding_context_ptr_array              = (EbPtr*)EB_NULL;
    enc_handle_ptr->dlf_context_ptr_array                         = (EbPtr*)EB_NULL;
    enc_handle_ptr->cdef_context_ptr_array                        = (EbPtr*)EB_NULL;
    enc_handle_ptr->rest_context_ptr_array                        = (EbPtr*)EB_NULL;
    enc_handle_ptr->packetization_context_ptr                     = (EbPtr)EB_NULL;

    // System Resource Managers
    enc_handle_ptr->input_buffer_resource_ptr                  = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->output_stream_buffer_resource_ptr_array    = (EbSystemResource**)EB_NULL;
    enc_handle_ptr->resource_coordination_results_resource_ptr = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->picture_analysis_results_resource_ptr      = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->picture_decision_results_resource_ptr      = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->motion_estimation_results_resource_ptr     = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->initial_rate_control_results_resource_ptr  = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->picture_demux_results_resource_ptr         = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->rate_control_tasks_resource_ptr            = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->rate_control_results_resource_ptr          = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->enc_dec_tasks_resource_ptr                 = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->enc_dec_results_resource_ptr               = (EbSystemResource*)EB_NULL;
    enc_handle_ptr->entropy_coding_results_resource_ptr        = (EbSystemResource*)EB_NULL;

    // Inter-Process Producer Fifos
    enc_handle_ptr->input_buffer_producer_fifo_ptr_array                  = (EbFifo**)EB_NULL;
    enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array      = (EbFifo***)EB_NULL;
    enc_handle_ptr->resource_coordination_results_producer_fifo_ptr_array = (EbFifo**)EB_NULL;
    enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array         = (EbFifo**)EB_NULL;
    enc_handle_ptr->picture_manager_results_producer_fifo_ptr_array       = (EbFifo**)EB_NULL;
    enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array            = (EbFifo**)EB_NULL;
    enc_handle_ptr->rate_control_results_producer_fifo_ptr_array          = (EbFifo**)EB_NULL;
    enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array                 = (EbFifo**)EB_NULL;
    enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array               = (EbFifo**)EB_NULL;
    enc_handle_ptr->entropy_coding_results_producer_fifo_ptr_array        = (EbFifo**)EB_NULL;
    enc_handle_ptr->dlf_results_producer_fifo_ptr_array                   = (EbFifo**)EB_NULL;
    enc_handle_ptr->cdef_results_producer_fifo_ptr_array                  = (EbFifo**)EB_NULL;
    enc_handle_ptr->rest_results_producer_fifo_ptr_array                  = (EbFifo**)EB_NULL;
    enc_handle_ptr->dlf_results_consumer_fifo_ptr_array                   = (EbFifo**)EB_NULL;
    enc_handle_ptr->cdef_results_consumer_fifo_ptr_array                  = (EbFifo**)EB_NULL;
    enc_handle_ptr->rest_results_consumer_fifo_ptr_array                  = (EbFifo**)EB_NULL;

    // Inter-Process Consumer Fifos
    enc_handle_ptr->input_buffer_consumer_fifo_ptr_array                  = (EbFifo**)EB_NULL;
    enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array      = (EbFifo***)EB_NULL;
    enc_handle_ptr->resource_coordination_results_consumer_fifo_ptr_array = (EbFifo**)EB_NULL;
    enc_handle_ptr->picture_demux_results_consumer_fifo_ptr_array         = (EbFifo**)EB_NULL;
    enc_handle_ptr->rate_control_tasks_consumer_fifo_ptr_array            = (EbFifo**)EB_NULL;
    enc_handle_ptr->rate_control_results_consumer_fifo_ptr_array          = (EbFifo**)EB_NULL;
    enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array                 = (EbFifo**)EB_NULL;
    enc_handle_ptr->enc_dec_results_consumer_fifo_ptr_array               = (EbFifo**)EB_NULL;
    enc_handle_ptr->entropy_coding_results_consumer_fifo_ptr_array        = (EbFifo**)EB_NULL;

    // Initialize Callbacks
    EB_MALLOC(EbCallback**, enc_handle_ptr->app_callback_ptr_array, sizeof(EbCallback*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
#if MEM_MAP_OPT
    EB_MALLOC(EbCallback*, enc_handle_ptr->app_callback_ptr_array[0], sizeof(EbCallback), EB_N_PTR);
    enc_handle_ptr->app_callback_ptr_array[0]->ErrorHandler = lib_svt_encoder_send_error_exit;
    enc_handle_ptr->app_callback_ptr_array[0]->handle = ebHandlePtr;
#else
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        EB_MALLOC(EbCallback*, enc_handle_ptr->app_callback_ptr_array[instance_index], sizeof(EbCallback), EB_N_PTR);
        enc_handle_ptr->app_callback_ptr_array[instance_index]->ErrorHandler = lib_svt_encoder_send_error_exit;
        enc_handle_ptr->app_callback_ptr_array[instance_index]->handle = ebHandlePtr;
    }
#endif

    // Initialize Sequence Control Set Instance Array
    EB_MALLOC(EbSequenceControlSetInstance**, enc_handle_ptr->sequence_control_set_instance_array, sizeof(EbSequenceControlSetInstance*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
#if MEM_MAP_OPT
    return_error = eb_sequence_control_set_instance_ctor(&enc_handle_ptr->sequence_control_set_instance_array[0]);
    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
#else
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        return_error = eb_sequence_control_set_instance_ctor(&enc_handle_ptr->sequence_control_set_instance_array[instance_index]);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
#endif
    return EB_ErrorNone;
}

EbErrorType EbInputBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr);

EbErrorType EbOutputReconBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr);

EbErrorType EbOutputBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr objectInitDataPtr);



EbErrorType DlfResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    DlfResults *context_ptr;
    EB_MALLOC(DlfResults*, context_ptr, sizeof(DlfResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
EbErrorType CdefResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    CdefResults *context_ptr;
    EB_MALLOC(CdefResults*, context_ptr, sizeof(CdefResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType RestResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    RestResults *context_ptr;
    EB_MALLOC(RestResults*, context_ptr, sizeof(RestResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

void init_fn_ptr(void);

/**********************************
* Initialize Encoder Library
**********************************/
#if defined(__linux__) || defined(__APPLE__)
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
#if !BUG_FIX_LOOKAHEAD
    uint32_t maxLookAheadDistance = 0;
#endif
    EbBool is16bit = (EbBool)(enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.encoder_color_format;

    /************************************
    * Plateform detection
    ************************************/
    if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.asm_type == 1) {
        enc_handle_ptr->sequence_control_set_instance_array[0]->encode_context_ptr->asm_type = GetCpuAsmType();
    }
    else {
        enc_handle_ptr->sequence_control_set_instance_array[0]->encode_context_ptr->asm_type = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.asm_type;
    }

    setup_rtcd_internal(enc_handle_ptr->sequence_control_set_instance_array[0]->encode_context_ptr->asm_type);
    asmSetConvolveAsmTable();

    init_intra_dc_predictors_c_internal();

    asmSetConvolveHbdAsmTable();

    init_intra_predictors_internal();
    EbSequenceControlSetInitData scs_init;
    scs_init.sb_size = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.super_block_size;

    build_blk_geom(scs_init.sb_size == 128);

    av1_init_me_luts();
    init_fn_ptr();

    /************************************
    * Sequence Control Set
    ************************************/
    return_error = eb_system_resource_ctor(
        &enc_handle_ptr->sequence_control_set_pool_ptr,
        enc_handle_ptr->sequence_control_set_pool_total_count,
        1,
        0,
        &enc_handle_ptr->sequence_control_set_pool_producer_fifo_ptr_array,
        (EbFifo ***)EB_NULL,
        EB_FALSE,
        eb_sequence_control_set_ctor,
        &scs_init);


    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    /************************************
    * Picture Control Set: Parent
    ************************************/
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->picture_parent_control_set_pool_ptr_array, sizeof(EbSystemResource*)  * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);


    EB_MALLOC(EbFifo***, enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array, sizeof(EbSystemResource**) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
#if !BUG_FIX_LOOKAHEAD
    // Updating the picture_control_set_pool_total_count based on the maximum look ahead distance
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        maxLookAheadDistance = MAX(maxLookAheadDistance, enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.look_ahead_distance);
    }
#endif


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
#if !BUG_FIX_LOOKAHEAD
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->picture_control_set_pool_init_count += maxLookAheadDistance;
#endif
        inputData.enc_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enc_mode;
        inputData.speed_control = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.speed_control_flag;
        inputData.film_grain_noise_level = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.film_grain_denoise_strength;
        inputData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.encoder_bit_depth;

        inputData.ext_block_flag = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.ext_block_flag;

        inputData.in_loop_me_flag = (uint8_t)enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.in_loop_me_flag;
#if MEMORY_FOOTPRINT_OPT_ME_MV
        inputData.mrp_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mrp_mode;
        inputData.nsq_present = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->nsq_present;
#endif
        return_error = eb_system_resource_ctor(
            &(enc_handle_ptr->picture_parent_control_set_pool_ptr_array[instance_index]),
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->picture_control_set_pool_init_count,//enc_handle_ptr->picture_control_set_pool_total_count,
            1,
            0,
            &enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            picture_parent_control_set_ctor,
            &inputData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    /************************************
    * Picture Control Set: Child
    ************************************/
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->picture_control_set_pool_ptr_array, sizeof(EbSystemResource*)  * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

    EB_MALLOC(EbFifo***, enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array, sizeof(EbSystemResource**) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

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
        inputData.color_format = color_format;
        inputData.sb_sz = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz;
        inputData.sb_size_pix = scs_init.sb_size;
        inputData.max_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_sb_depth;
#if MEMORY_FOOTPRINT_OPT_ME_MV
        inputData.cdf_mode = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->cdf_mode;
#endif
        return_error = eb_system_resource_ctor(
            &(enc_handle_ptr->picture_control_set_pool_ptr_array[instance_index]),
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->picture_control_set_pool_init_count_child, //EB_PictureControlSetPoolInitCountChild,
            1,
            0,
            &enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            picture_control_set_ctor,
            &inputData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    /************************************
    * Picture Buffers
    ************************************/

    // Allocate Resource Arrays
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->reference_picture_pool_ptr_array, sizeof(EbSystemResource*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->pa_reference_picture_pool_ptr_array, sizeof(EbSystemResource*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

#if ALT_REF_OVERLAY
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->overlay_input_picture_pool_ptr_array, sizeof(EbSystemResource*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
    EB_MALLOC(EbFifo***, enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array, sizeof(EbFifo**) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
#endif
    // Allocate Producer Fifo Arrays
    EB_MALLOC(EbFifo***, enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array, sizeof(EbFifo**) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
    EB_MALLOC(EbFifo***, enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array, sizeof(EbFifo**) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

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
        EbPictureBufferDescInitData       quarterDecimPictureBufferDescInitData;
        EbPictureBufferDescInitData       sixteenthDecimPictureBufferDescInitData;

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
#if UNPACK_REF_POST_EP // constructor
        // Hsan: split_mode is set @ eb_reference_object_ctor() as both unpacked reference and packed reference are needed for a 10BIT input; unpacked reference @ MD, and packed reference @ EP
#else
        referencePictureBufferDescInitData.split_mode = EB_FALSE;
#endif

        if (is16bit)
            referencePictureBufferDescInitData.bit_depth = EB_10BIT;


        EbReferenceObjectDescInitDataStructure.reference_picture_desc_init_data = referencePictureBufferDescInitData;

        // Reference Picture Buffers
        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->reference_picture_pool_ptr_array[instance_index],
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->reference_picture_buffer_init_count,//enc_handle_ptr->reference_picture_pool_total_count,
            EB_PictureManagerProcessInitCount,
            0,
            &enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            eb_reference_object_ctor,
            &(EbReferenceObjectDescInitDataStructure));

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

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

        quarterDecimPictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width >> 1;
        quarterDecimPictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height >> 1;
        quarterDecimPictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        quarterDecimPictureBufferDescInitData.color_format = EB_YUV420;
        quarterDecimPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        quarterDecimPictureBufferDescInitData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterDecimPictureBufferDescInitData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterDecimPictureBufferDescInitData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterDecimPictureBufferDescInitData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 1;
        quarterDecimPictureBufferDescInitData.split_mode = EB_FALSE;


        sixteenthDecimPictureBufferDescInitData.max_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width >> 2;
        sixteenthDecimPictureBufferDescInitData.max_height = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height >> 2;
        sixteenthDecimPictureBufferDescInitData.bit_depth = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->encoder_bit_depth;
        sixteenthDecimPictureBufferDescInitData.color_format = EB_YUV420;
        sixteenthDecimPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        sixteenthDecimPictureBufferDescInitData.left_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthDecimPictureBufferDescInitData.right_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthDecimPictureBufferDescInitData.top_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthDecimPictureBufferDescInitData.bot_padding = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->sb_sz >> 2;
        sixteenthDecimPictureBufferDescInitData.split_mode = EB_FALSE;

        EbPaReferenceObjectDescInitDataStructure.reference_picture_desc_init_data = referencePictureBufferDescInitData;
        EbPaReferenceObjectDescInitDataStructure.quarter_picture_desc_init_data = quarterDecimPictureBufferDescInitData;
        EbPaReferenceObjectDescInitDataStructure.sixteenth_picture_desc_init_data = sixteenthDecimPictureBufferDescInitData;

        // Reference Picture Buffers
        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->pa_reference_picture_pool_ptr_array[instance_index],
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pa_reference_picture_buffer_init_count,
            EB_PictureDecisionProcessInitCount,
            0,
            &enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
            (EbFifo ***)EB_NULL,
            EB_FALSE,
            eb_pa_reference_object_ctor,
            &(EbPaReferenceObjectDescInitDataStructure));
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        // Set the SequenceControlSet Picture Pool Fifo Ptrs
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->reference_picture_pool_fifo_ptr = (enc_handle_ptr->reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->pa_reference_picture_pool_fifo_ptr = (enc_handle_ptr->pa_reference_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];

#if ALT_REF_OVERLAY
        if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.enable_overlays) {
            // Overlay Input Picture Buffers
            return_error = eb_system_resource_ctor(
                &enc_handle_ptr->overlay_input_picture_pool_ptr_array[instance_index],
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->overlay_input_picture_buffer_init_count,
                1,
                0,
                &enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array[instance_index],
                (EbFifo ***)EB_NULL,
                EB_FALSE,
                EbInputBufferHeaderCtor,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

            if (return_error == EB_ErrorInsufficientResources) {
                return EB_ErrorInsufficientResources;
            }

            // Set the SequenceControlSet Overlay input Picture Pool Fifo Ptrs
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->overlay_input_picture_pool_fifo_ptr = (enc_handle_ptr->overlay_input_picture_pool_producer_fifo_ptr_dbl_array[instance_index])[0];
        }
#endif
    }

    /************************************
    * System Resource Managers & Fifos
    ************************************/

    // EbBufferHeaderType Input
    return_error = eb_system_resource_ctor(
        &enc_handle_ptr->input_buffer_resource_ptr,
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->input_buffer_fifo_init_count,
        1,
        EB_ResourceCoordinationProcessInitCount,
        &enc_handle_ptr->input_buffer_producer_fifo_ptr_array,
        &enc_handle_ptr->input_buffer_consumer_fifo_ptr_array,
        EB_TRUE,
        EbInputBufferHeaderCtor,
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // EbBufferHeaderType Output Stream
    EB_MALLOC(EbSystemResource**, enc_handle_ptr->output_stream_buffer_resource_ptr_array, sizeof(EbSystemResource*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
    EB_MALLOC(EbFifo***, enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array, sizeof(EbFifo**)          * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
    EB_MALLOC(EbFifo***, enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array, sizeof(EbFifo**)          * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->output_stream_buffer_resource_ptr_array[instance_index],
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->output_stream_buffer_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->total_process_init_count,//EB_PacketizationProcessInitCount,
            1,
            &enc_handle_ptr->output_stream_buffer_producer_fifo_ptr_dbl_array[instance_index],
            &enc_handle_ptr->output_stream_buffer_consumer_fifo_ptr_dbl_array[instance_index],
            EB_TRUE,
            EbOutputBufferHeaderCtor,
            &enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    if (enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.recon_enabled) {
        // EbBufferHeaderType Output Recon
        EB_MALLOC(EbSystemResource**, enc_handle_ptr->output_recon_buffer_resource_ptr_array, sizeof(EbSystemResource*) * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
        EB_MALLOC(EbFifo***, enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array, sizeof(EbFifo**)          * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);
        EB_MALLOC(EbFifo***, enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr_dbl_array, sizeof(EbFifo**)          * enc_handle_ptr->encode_instance_total_count, EB_N_PTR);

        for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
            return_error = eb_system_resource_ctor(
                &enc_handle_ptr->output_recon_buffer_resource_ptr_array[instance_index],
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->output_recon_buffer_fifo_init_count,
                enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->enc_dec_process_init_count,
                1,
                &enc_handle_ptr->output_recon_buffer_producer_fifo_ptr_dbl_array[instance_index],
                &enc_handle_ptr->output_recon_buffer_consumer_fifo_ptr_dbl_array[instance_index],
                EB_TRUE,
                EbOutputReconBufferHeaderCtor,
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr);
            if (return_error == EB_ErrorInsufficientResources) {
                return EB_ErrorInsufficientResources;
            }
        }
    }

    // Resource Coordination Results
    {
        ResourceCoordinationResultInitData resourceCoordinationResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->resource_coordination_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->resource_coordination_fifo_init_count,
            EB_ResourceCoordinationProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count,
            &enc_handle_ptr->resource_coordination_results_producer_fifo_ptr_array,
            &enc_handle_ptr->resource_coordination_results_consumer_fifo_ptr_array,
            EB_TRUE,
            resource_coordination_result_ctor,
            &resourceCoordinationResultInitData);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }


    // Picture Analysis Results
    {
        PictureAnalysisResultInitData pictureAnalysisResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->picture_analysis_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count,
            EB_PictureDecisionProcessInitCount,
            &enc_handle_ptr->picture_analysis_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_analysis_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_analysis_result_ctor,
            &pictureAnalysisResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Decision Results
    {
        PictureDecisionResultInitData pictureDecisionResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->picture_decision_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_decision_fifo_init_count,
            EB_PictureDecisionProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count,
            &enc_handle_ptr->picture_decision_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_decision_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_decision_result_ctor,
            &pictureDecisionResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Motion Estimation Results
    {
        MotionEstimationResultsInitData motionEstimationResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->motion_estimation_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count,
            EB_InitialRateControlProcessInitCount,
            &enc_handle_ptr->motion_estimation_results_producer_fifo_ptr_array,
            &enc_handle_ptr->motion_estimation_results_consumer_fifo_ptr_array,
            EB_TRUE,
            motion_estimation_results_ctor,
            &motionEstimationResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Initial Rate Control Results
    {
        InitialRateControlResultInitData initialRateControlResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->initial_rate_control_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->initial_rate_control_fifo_init_count,
            EB_InitialRateControlProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count,
            &enc_handle_ptr->initial_rate_control_results_producer_fifo_ptr_array,
            &enc_handle_ptr->initial_rate_control_results_consumer_fifo_ptr_array,
            EB_TRUE,
            initial_rate_control_results_ctor,
            &initialRateControlResultInitData);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Demux Results
    {
        PictureResultInitData pictureResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->picture_demux_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_demux_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count + enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count,
            EB_PictureManagerProcessInitCount,
            &enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array,
            &enc_handle_ptr->picture_demux_results_consumer_fifo_ptr_array,
            EB_TRUE,
            picture_results_ctor,
            &pictureResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Rate Control Tasks
    {
        RateControlTasksInitData rateControlTasksInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->rate_control_tasks_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rate_control_tasks_fifo_init_count,
            RateControlPortTotalCount(),
            EB_RateControlProcessInitCount,
            &enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array,
            &enc_handle_ptr->rate_control_tasks_consumer_fifo_ptr_array,
            EB_TRUE,
            rate_control_tasks_ctor,
            &rateControlTasksInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Rate Control Results
    {
        RateControlResultsInitData rateControlResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->rate_control_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rate_control_fifo_init_count,
            EB_RateControlProcessInitCount,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count,
            &enc_handle_ptr->rate_control_results_producer_fifo_ptr_array,
            &enc_handle_ptr->rate_control_results_consumer_fifo_ptr_array,
            EB_TRUE,
            rate_control_results_ctor,
            &rateControlResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
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

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->enc_dec_tasks_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_fifo_init_count,
            EncDecPortTotalCount(),
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count,
            &enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array,
            &enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array,
            EB_TRUE,
            enc_dec_tasks_ctor,
            &ModeDecisionResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // EncDec Results
    {
        EncDecResultsInitData encDecResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->enc_dec_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count,
            &enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array,
            &enc_handle_ptr->enc_dec_results_consumer_fifo_ptr_array,
            EB_TRUE,
            enc_dec_results_ctor,
            &encDecResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    //DLF results
    {
        EntropyCodingResultsInitData dlfResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->dlf_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count,
            &enc_handle_ptr->dlf_results_producer_fifo_ptr_array,
            &enc_handle_ptr->dlf_results_consumer_fifo_ptr_array,
            EB_TRUE,
            DlfResultsCtor,
            &dlfResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    //CDEF results
    {
        EntropyCodingResultsInitData cdefResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->cdef_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count,
            &enc_handle_ptr->cdef_results_producer_fifo_ptr_array,
            &enc_handle_ptr->cdef_results_consumer_fifo_ptr_array,
            EB_TRUE,
            CdefResultsCtor,
            &cdefResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

    }
    //REST results
    {
        EntropyCodingResultsInitData restResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->rest_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count,
            &enc_handle_ptr->rest_results_producer_fifo_ptr_array,
            &enc_handle_ptr->rest_results_consumer_fifo_ptr_array,
            EB_TRUE,
            RestResultsCtor,
            &restResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }


    // Entropy Coding Results
    {
        EntropyCodingResultsInitData entropyCodingResultInitData;

        return_error = eb_system_resource_ctor(
            &enc_handle_ptr->entropy_coding_results_resource_ptr,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_fifo_init_count,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count,
            EB_PacketizationProcessInitCount,
            &enc_handle_ptr->entropy_coding_results_producer_fifo_ptr_array,
            &enc_handle_ptr->entropy_coding_results_consumer_fifo_ptr_array,
            EB_TRUE,
            entropy_coding_results_ctor,
            &entropyCodingResultInitData);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    /************************************
    * App Callbacks
    ************************************/
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        enc_handle_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->app_callback_ptr = enc_handle_ptr->app_callback_ptr_array[instance_index];
    }

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
    return_error = resource_coordination_context_ctor(
        (ResourceCoordinationContext**)&enc_handle_ptr->resource_coordination_context_ptr,
        enc_handle_ptr->input_buffer_consumer_fifo_ptr_array[0],
        enc_handle_ptr->resource_coordination_results_producer_fifo_ptr_array[0],
        enc_handle_ptr->picture_parent_control_set_pool_producer_fifo_ptr_dbl_array[0],//ResourceCoordination works with ParentPCS
        enc_handle_ptr->sequence_control_set_instance_array,
        enc_handle_ptr->sequence_control_set_pool_producer_fifo_ptr_array[0],
        enc_handle_ptr->app_callback_ptr_array,
        enc_handle_ptr->compute_segments_total_count_array,
        enc_handle_ptr->encode_instance_total_count);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Picture Analysis Context
    EB_MALLOC(EbPtr*, enc_handle_ptr->picture_analysis_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count; ++processIndex) {

        EbPictureBufferDescInitData  pictureBufferDescConf;
        pictureBufferDescConf.max_width = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width;
        pictureBufferDescConf.max_height = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height;
        pictureBufferDescConf.bit_depth = EB_8BIT;
        pictureBufferDescConf.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
        pictureBufferDescConf.left_padding = 0;
        pictureBufferDescConf.right_padding = 0;
        pictureBufferDescConf.top_padding = 0;
        pictureBufferDescConf.bot_padding = 0;
        pictureBufferDescConf.split_mode = EB_FALSE;

        return_error = picture_analysis_context_ctor(
            &pictureBufferDescConf,
            EB_TRUE,
            (PictureAnalysisContext**)&enc_handle_ptr->picture_analysis_context_ptr_array[processIndex],
            enc_handle_ptr->resource_coordination_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_analysis_results_producer_fifo_ptr_array[processIndex]);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Decision Context
    {
        // Initialize the various Picture types

        instance_index = 0;


        return_error = picture_decision_context_ctor(
            (PictureDecisionContext**)&enc_handle_ptr->picture_decision_context_ptr,
            enc_handle_ptr->picture_analysis_results_consumer_fifo_ptr_array[0],
            enc_handle_ptr->picture_decision_results_producer_fifo_ptr_array[0]);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Motion Analysis Context
    EB_MALLOC(EbPtr*, enc_handle_ptr->motion_estimation_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count; ++processIndex) {
#if MEMORY_FOOTPRINT_OPT_ME_MV
        return_error = motion_estimation_context_ctor(
            (MotionEstimationContext_t**)&enc_handle_ptr->motion_estimation_context_ptr_array[processIndex],
            enc_handle_ptr->picture_decision_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->motion_estimation_results_producer_fifo_ptr_array[processIndex],
#if REDUCE_ME_SEARCH_AREA
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height,
#endif
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->nsq_present,
            enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->mrp_mode);
#else
        return_error = motion_estimation_context_ctor(
            (MotionEstimationContext_t**)&enc_handle_ptr->motion_estimation_context_ptr_array[processIndex],
            enc_handle_ptr->picture_decision_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->motion_estimation_results_producer_fifo_ptr_array[processIndex]);
#endif
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Initial Rate Control Context
    return_error = initial_rate_control_context_ctor(
        (InitialRateControlContext**)&enc_handle_ptr->initial_rate_control_context_ptr,
        enc_handle_ptr->motion_estimation_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->initial_rate_control_results_producer_fifo_ptr_array[0]);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Source Based Operations Context
    EB_MALLOC(EbPtr*, enc_handle_ptr->source_based_operations_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count; ++processIndex) {

        return_error = source_based_operations_context_ctor(
            (SourceBasedOperationsContext**)&enc_handle_ptr->source_based_operations_context_ptr_array[processIndex],
            enc_handle_ptr->initial_rate_control_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Manager Context
    return_error = picture_manager_context_ctor(
        (PictureManagerContext**)&enc_handle_ptr->picture_manager_context_ptr,
        enc_handle_ptr->picture_demux_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER, 0)],
        enc_handle_ptr->picture_control_set_pool_producer_fifo_ptr_dbl_array[0]);//The Child PCS Pool here
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Rate Control Context
    return_error = rate_control_context_ctor(
        (RateControlContext**)&enc_handle_ptr->rate_control_context_ptr,
        enc_handle_ptr->rate_control_tasks_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_results_producer_fifo_ptr_array[0],
        enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->intra_period_length);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }


    // Mode Decision Configuration Contexts
    {
        // Mode Decision Configuration Contexts
        EB_MALLOC(EbPtr*, enc_handle_ptr->mode_decision_configuration_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count, EB_N_PTR);

        for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count; ++processIndex) {
            return_error = mode_decision_configuration_context_ctor(
                (ModeDecisionConfigurationContext**)&enc_handle_ptr->mode_decision_configuration_context_ptr_array[processIndex],
                enc_handle_ptr->rate_control_results_consumer_fifo_ptr_array[processIndex],

                enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array[EncDecPortLookup(ENCDEC_INPUT_PORT_MDC, processIndex)],
                ((enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64) *
                ((enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64));


            if (return_error == EB_ErrorInsufficientResources) {
                return EB_ErrorInsufficientResources;
            }
        }
    }

    max_picture_width = 0;
    for (instance_index = 0; instance_index < enc_handle_ptr->encode_instance_total_count; ++instance_index) {
        if (max_picture_width < enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width) {
            max_picture_width = enc_handle_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
        }
    }

    // EncDec Contexts
    EB_MALLOC(EbPtr*, enc_handle_ptr->enc_dec_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count; ++processIndex) {
        return_error = enc_dec_context_ctor(
            (EncDecContext**)&enc_handle_ptr->enc_dec_context_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_tasks_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_tasks_producer_fifo_ptr_array[EncDecPortLookup(ENCDEC_INPUT_PORT_ENCDEC, processIndex)],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[
                enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count+
                //1 +
                    processIndex], // Add port lookup logic here JMJ
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Dlf Contexts
    EB_MALLOC(EbPtr*, enc_handle_ptr->dlf_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count; ++processIndex) {
        return_error = dlf_context_ctor(
            (DlfContext**)&enc_handle_ptr->dlf_context_ptr_array[processIndex],
            enc_handle_ptr->enc_dec_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->dlf_results_producer_fifo_ptr_array[processIndex],             //output to EC
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    //CDEF Contexts
    EB_MALLOC(EbPtr*, enc_handle_ptr->cdef_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count; ++processIndex) {
        return_error = cdef_context_ctor(
            (CdefContext_t**)&enc_handle_ptr->cdef_context_ptr_array[processIndex],
            enc_handle_ptr->dlf_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->cdef_results_producer_fifo_ptr_array[processIndex],
            is16bit,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    //Rest Contexts
    EB_MALLOC(EbPtr*, enc_handle_ptr->rest_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count; ++processIndex) {
        return_error = rest_context_ctor(
            (RestContext**)&enc_handle_ptr->rest_context_ptr_array[processIndex],
            enc_handle_ptr->cdef_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->rest_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->picture_demux_results_producer_fifo_ptr_array[
                /*enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count*/ 1+ processIndex],
            is16bit,
            color_format,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_width,
            enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->max_input_luma_height
        );

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }


    // Entropy Coding Contexts
    EB_MALLOC(EbPtr*, enc_handle_ptr->entropy_coding_context_ptr_array, sizeof(EbPtr) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count; ++processIndex) {
        return_error = entropy_coding_context_ctor(
            (EntropyCodingContext**)&enc_handle_ptr->entropy_coding_context_ptr_array[processIndex],
            enc_handle_ptr->rest_results_consumer_fifo_ptr_array[processIndex],
            enc_handle_ptr->entropy_coding_results_producer_fifo_ptr_array[processIndex],
            enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_ENTROPY_CODING, processIndex)],
            is16bit);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Packetization Context
    return_error = packetization_context_ctor(
        (PacketizationContext**)&enc_handle_ptr->packetization_context_ptr,
        enc_handle_ptr->entropy_coding_results_consumer_fifo_ptr_array[0],
        enc_handle_ptr->rate_control_tasks_producer_fifo_ptr_array[RateControlPortLookup(RATE_CONTROL_INPUT_PORT_PACKETIZATION, 0)]);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    /************************************
    * Thread Handles
    ************************************/
    EbSvtAv1EncConfiguration   *config_ptr = &enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config;

    EbSetThreadManagementParameters(config_ptr);

    // Resource Coordination
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->resource_coordination_thread_handle, sizeof(EbHandle), EB_THREAD, resource_coordination_kernel, enc_handle_ptr->resource_coordination_context_ptr);

    // Picture Analysis
    EB_MALLOC(EbHandle*, enc_handle_ptr->picture_analysis_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->picture_analysis_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->picture_analysis_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, picture_analysis_kernel, enc_handle_ptr->picture_analysis_context_ptr_array[processIndex]);
    }

    // Picture Decision
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->picture_decision_thread_handle, sizeof(EbHandle), EB_THREAD, picture_decision_kernel, enc_handle_ptr->picture_decision_context_ptr);

    // Motion Estimation
    EB_MALLOC(EbHandle*, enc_handle_ptr->motion_estimation_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->motion_estimation_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->motion_estimation_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, motion_estimation_kernel, enc_handle_ptr->motion_estimation_context_ptr_array[processIndex]);
    }

    // Initial Rate Control
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->initial_rate_control_thread_handle, sizeof(EbHandle), EB_THREAD, initial_rate_control_kernel, enc_handle_ptr->initial_rate_control_context_ptr);

    // Source Based Oprations
    EB_MALLOC(EbHandle*, enc_handle_ptr->source_based_operations_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->source_based_operations_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->source_based_operations_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, source_based_operations_kernel, enc_handle_ptr->source_based_operations_context_ptr_array[processIndex]);
    }

    // Picture Manager
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->picture_manager_thread_handle, sizeof(EbHandle), EB_THREAD, picture_manager_kernel, enc_handle_ptr->picture_manager_context_ptr);

    // Rate Control
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->rate_control_thread_handle, sizeof(EbHandle), EB_THREAD, rate_control_kernel, enc_handle_ptr->rate_control_context_ptr);

    // Mode Decision Configuration Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->mode_decision_configuration_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->mode_decision_configuration_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->mode_decision_configuration_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, mode_decision_configuration_kernel, enc_handle_ptr->mode_decision_configuration_context_ptr_array[processIndex]);
    }

    // EncDec Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->enc_dec_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->enc_dec_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->enc_dec_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, enc_dec_kernel, enc_handle_ptr->enc_dec_context_ptr_array[processIndex]);
    }

    // Dlf Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->dlf_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->dlf_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->dlf_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, dlf_kernel, enc_handle_ptr->dlf_context_ptr_array[processIndex]);
    }


    // Cdef Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->cdef_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->cdef_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->cdef_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, cdef_kernel, enc_handle_ptr->cdef_context_ptr_array[processIndex]);
    }

    // Rest Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->rest_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->rest_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->rest_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, rest_kernel, enc_handle_ptr->rest_context_ptr_array[processIndex]);
    }

    // Entropy Coding Process
    EB_MALLOC(EbHandle*, enc_handle_ptr->entropy_coding_thread_handle_array, sizeof(EbHandle) * enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count, EB_N_PTR);

    for (processIndex = 0; processIndex < enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->entropy_coding_process_init_count; ++processIndex) {
        EB_CREATETHREAD(EbHandle, enc_handle_ptr->entropy_coding_thread_handle_array[processIndex], sizeof(EbHandle), EB_THREAD, entropy_coding_kernel, enc_handle_ptr->entropy_coding_context_ptr_array[processIndex]);
    }

    // Packetization
    EB_CREATETHREAD(EbHandle, enc_handle_ptr->packetization_thread_handle, sizeof(EbHandle), EB_THREAD, packetization_kernel, enc_handle_ptr->packetization_context_ptr);


#if DISPLAY_MEMORY
    EB_MEMORY();
#endif
    return return_error;
}

/**********************************
* DeInitialize Encoder Library
**********************************/
#if defined(__linux__) || defined(__APPLE__)
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_deinit_encoder(EbComponentType *svt_enc_component){

    if(svt_enc_component == NULL)
        return EB_ErrorBadParameter;
#if MEM_MAP_OPT
    EbEncHandle         *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbErrorType          return_error = EB_ErrorNone;

    if (enc_handle_ptr) {
        if (memory_map) {
            // Loop through the ptr table and free all malloc'd pointers per channel
            EbMemoryMapEntry*    memory_entry = memory_map;
            if (memory_entry){
                do {
                    switch (memory_entry->ptr_type) {
                        case EB_N_PTR:
                            free(memory_entry->ptr);
                            break;
                        case EB_A_PTR:
#ifdef _WIN32
                            _aligned_free(memory_entry->ptr);
#else
                            free(memory_entry->ptr);
#endif
                            break;
                        case EB_SEMAPHORE:
                            eb_destroy_semaphore(memory_entry->ptr);
                            break;
                        case EB_THREAD:
                            eb_destroy_thread(memory_entry->ptr);
                            break;
                        case EB_MUTEX:
                            eb_destroy_mutex(memory_entry->ptr);
                            break;
                        default:
                            return_error = EB_ErrorMax;
                            break;
                    }
                    EbMemoryMapEntry*    tmp_memory_entry = memory_entry;
                    memory_entry = (EbMemoryMapEntry*)tmp_memory_entry->prev_entry;
                    if (tmp_memory_entry) free(tmp_memory_entry);
                } while(memory_entry != enc_handle_ptr->memory_map_init_address && memory_entry);
                if (enc_handle_ptr->memory_map_init_address) free(enc_handle_ptr->memory_map_init_address);
            }
        }
    }
#else
    EbEncHandle *enc_handle_ptr = (EbEncHandle*)svt_enc_component->p_component_private;
    EbErrorType return_error = EB_ErrorNone;
    int32_t              ptrIndex = 0;
    EbMemoryMapEntry*   memoryEntry = (EbMemoryMapEntry*)EB_NULL;

    if (enc_handle_ptr) {
        if (enc_handle_ptr->memory_map_index) {
            // Loop through the ptr table and free all malloc'd pointers per channel
            for (ptrIndex = (enc_handle_ptr->memory_map_index) - 1; ptrIndex >= 0; --ptrIndex) {
                memoryEntry = &enc_handle_ptr->memory_map[ptrIndex];
                switch (memoryEntry->ptr_type) {
                case EB_N_PTR:
                    free(memoryEntry->ptr);
                    break;
                case EB_A_PTR:
#ifdef _WIN32
                    _aligned_free(memoryEntry->ptr);
#else
                    free(memoryEntry->ptr);
#endif
                    break;
                case EB_SEMAPHORE:
                    eb_destroy_semaphore(memoryEntry->ptr);
                    break;
                case EB_THREAD:
                    eb_destroy_thread(memoryEntry->ptr);
                    break;
                case EB_MUTEX:
                    eb_destroy_mutex(memoryEntry->ptr);
                    break;
                default:
                    return_error = EB_ErrorMax;
                    break;
                }
            }
            if (enc_handle_ptr->memory_map != (EbMemoryMapEntry*)NULL) {
                free(enc_handle_ptr->memory_map);
            }

        }
    }
#endif
    return return_error;
}

EbErrorType eb_svt_enc_init_parameter(
    EbSvtAv1EncConfiguration * config_ptr);

EbErrorType init_svt_av1_encoder_handle(
    EbComponentType * hComponent);
/**********************************
* GetHandle
**********************************/
#if defined(__linux__) || defined(__APPLE__)
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

    *p_handle = (EbComponentType*)malloc(sizeof(EbComponentType));
    if (*p_handle != (EbComponentType*)NULL) {

        // Init Component OS objects (threads, semaphores, etc.)
        // also links the various Component control functions
        return_error = init_svt_av1_encoder_handle(*p_handle);

        if (return_error == EB_ErrorNone) {
            ((EbComponentType*)(*p_handle))->p_application_private = p_app_data;

        }
        else if (return_error == EB_ErrorInsufficientResources) {
            eb_deinit_encoder((EbComponentType*)NULL);
            *p_handle = (EbComponentType*)NULL;
        }
        else
            return_error = EB_ErrorInvalidComponent;
    }
    else {
        SVT_LOG("Error: Component Struct Malloc Failed\n");
        return_error = EB_ErrorInsufficientResources;
    }
    return_error = eb_svt_enc_init_parameter(config_ptr);

    return return_error;
}

/**********************************
* Encoder Componenet DeInit
**********************************/
EbErrorType eb_h265_enc_component_de_init(EbComponentType  *svt_enc_component)
{
    EbErrorType       return_error = EB_ErrorNone;

    if (svt_enc_component->p_component_private) {
        free((EbEncHandle *)svt_enc_component->p_component_private);
    }
    else {
        return_error = EB_ErrorUndefined;
    }

    return return_error;
}

/**********************************
* eb_deinit_handle
**********************************/
#if defined(__linux__) || defined(__APPLE__)
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_deinit_handle(
    EbComponentType  *svt_enc_component)
{
    EbErrorType return_error = EB_ErrorNone;

    if (svt_enc_component) {
        return_error = eb_h265_enc_component_de_init(svt_enc_component);

        free(svt_enc_component);
    }
    else {
        return_error = EB_ErrorInvalidComponent;
    }

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
    }
    else {

        sequence_control_set_ptr->max_input_pad_right = 0;
    }
    if (sequence_control_set_ptr->max_input_luma_height % MIN_BLOCK_SIZE) {

        sequence_control_set_ptr->max_input_pad_bottom = MIN_BLOCK_SIZE - (sequence_control_set_ptr->max_input_luma_height % MIN_BLOCK_SIZE);
        sequence_control_set_ptr->max_input_luma_height = sequence_control_set_ptr->max_input_luma_height + sequence_control_set_ptr->max_input_pad_bottom;
    }
    else {
        sequence_control_set_ptr->max_input_pad_bottom = 0;
    }

    sequence_control_set_ptr->max_input_chroma_width = sequence_control_set_ptr->max_input_luma_width >> subsampling_x;
    sequence_control_set_ptr->max_input_chroma_height = sequence_control_set_ptr->max_input_luma_height >> subsampling_y;


    // Configure the padding
    sequence_control_set_ptr->left_padding  = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->top_padding = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->right_padding = BLOCK_SIZE_64 + 4;
    sequence_control_set_ptr->bot_padding = BLOCK_SIZE_64 + 4;

    sequence_control_set_ptr->chroma_width = sequence_control_set_ptr->max_input_luma_width >> subsampling_x;
    sequence_control_set_ptr->chroma_height = sequence_control_set_ptr->max_input_luma_height >> subsampling_y;
    sequence_control_set_ptr->luma_width = sequence_control_set_ptr->max_input_luma_width;
    sequence_control_set_ptr->luma_height = sequence_control_set_ptr->max_input_luma_height;
    sequence_control_set_ptr->static_config.source_width = sequence_control_set_ptr->max_input_luma_width;
    sequence_control_set_ptr->static_config.source_height = sequence_control_set_ptr->max_input_luma_height;

    derive_input_resolution(
        sequence_control_set_ptr,
        sequence_control_set_ptr->luma_width*sequence_control_set_ptr->luma_height);
#if NEW_PRESETS
    sequence_control_set_ptr->static_config.super_block_size       = (sequence_control_set_ptr->static_config.enc_mode == ENC_M0 && sequence_control_set_ptr->input_resolution >= INPUT_SIZE_1080i_RANGE) ? 128 : 64;
#else
    sequence_control_set_ptr->static_config.super_block_size       = (sequence_control_set_ptr->static_config.enc_mode <= ENC_M1 && sequence_control_set_ptr->input_resolution >= INPUT_SIZE_1080i_RANGE) ? 128 : 64;
#endif
#if RC
    sequence_control_set_ptr->static_config.super_block_size = (sequence_control_set_ptr->static_config.rate_control_mode > 1) ? 64 : sequence_control_set_ptr->static_config.super_block_size;
   // sequence_control_set_ptr->static_config.hierarchical_levels = (sequence_control_set_ptr->static_config.rate_control_mode > 1) ? 3 : sequence_control_set_ptr->static_config.hierarchical_levels;
#endif
#if ALT_REF_OVERLAY
    sequence_control_set_ptr->static_config.enable_overlays = sequence_control_set_ptr->static_config.enable_altrefs == EB_FALSE ||
        (sequence_control_set_ptr->static_config.rate_control_mode > 0) ||
        (sequence_control_set_ptr->static_config.enc_mode > ENC_M0) ||
        sequence_control_set_ptr->static_config.encoder_bit_depth != EB_8BIT ?
        0 : sequence_control_set_ptr->static_config.enable_overlays;
#endif

#if MEMORY_FOOTPRINT_OPT_ME_MV
    //0: MRP Mode 0 (4,3)
    //1: MRP Mode 1 (2,2)                            
    sequence_control_set_ptr->mrp_mode = (uint8_t) (sequence_control_set_ptr->static_config.enc_mode == ENC_M0) ? 0 : 1;

    //0: ON
    //1: OFF                            
    sequence_control_set_ptr->cdf_mode = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode <= ENC_M6) ? 0 : 1;


    //0: NSQ absent
    //1: NSQ present    
#if REDUCE_BLOCK_COUNT_ME
    sequence_control_set_ptr->nsq_present = (uint8_t)(sequence_control_set_ptr->static_config.enc_mode <= ENC_M5) ? 1 : 0;
#else
    sequence_control_set_ptr->nsq_present = 1;
#endif
#endif
}

void CopyApiFromApp(
    SequenceControlSet       *sequence_control_set_ptr,
    EbSvtAv1EncConfiguration   *pComponentParameterStructure){

    uint32_t                  hmeRegionIndex = 0;

    sequence_control_set_ptr->max_input_luma_width = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->source_width;
    sequence_control_set_ptr->max_input_luma_height = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->source_height;
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
    else {
        sequence_control_set_ptr->conformance_window_flag = 1;
    }
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
    sequence_control_set_ptr->intra_period_length = sequence_control_set_ptr->static_config.intra_period_length;
    sequence_control_set_ptr->intra_refresh_type = sequence_control_set_ptr->static_config.intra_refresh_type;
    sequence_control_set_ptr->max_temporal_layers = sequence_control_set_ptr->static_config.hierarchical_levels;
    sequence_control_set_ptr->static_config.use_qp_file = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->use_qp_file;

#if SHUT_FILTERING
    sequence_control_set_ptr->static_config.disable_dlf_flag = 1;//
#else
    // Deblock Filter
    sequence_control_set_ptr->static_config.disable_dlf_flag = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->disable_dlf_flag;
#endif

    // Local Warped Motion
#if ENABLE_WARPED_MV
    sequence_control_set_ptr->static_config.enable_warped_motion = EB_TRUE;
#else
    sequence_control_set_ptr->static_config.enable_warped_motion = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->enable_warped_motion;
#endif

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
    sequence_control_set_ptr->static_config.constrained_intra = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->constrained_intra;

    // Adaptive Loop Filter
    sequence_control_set_ptr->static_config.tile_rows = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->tile_rows;
    sequence_control_set_ptr->static_config.tile_columns = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->tile_columns;


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
    sequence_control_set_ptr->static_config.improve_sharpness = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->improve_sharpness;
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
    sequence_control_set_ptr->static_config.asm_type = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->asm_type;

    sequence_control_set_ptr->static_config.channel_id = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->channel_id;
    sequence_control_set_ptr->static_config.active_channel_count = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->active_channel_count;
    sequence_control_set_ptr->static_config.logical_processors = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->logical_processors;
    sequence_control_set_ptr->static_config.target_socket = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->target_socket;
    sequence_control_set_ptr->qp = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->qp;
    sequence_control_set_ptr->static_config.recon_enabled = ((EbSvtAv1EncConfiguration*)pComponentParameterStructure)->recon_enabled;

    // Extract frame rate from Numerator and Denominator if not 0
    if (sequence_control_set_ptr->static_config.frame_rate_numerator != 0 && sequence_control_set_ptr->static_config.frame_rate_denominator != 0) {
        sequence_control_set_ptr->frame_rate = sequence_control_set_ptr->static_config.frame_rate = (((sequence_control_set_ptr->static_config.frame_rate_numerator << 8) / (sequence_control_set_ptr->static_config.frame_rate_denominator)) << 8);
    }

    // Get Default Intra Period if not specified
    if (sequence_control_set_ptr->static_config.intra_period_length == -2) {
        sequence_control_set_ptr->intra_period_length = sequence_control_set_ptr->static_config.intra_period_length = compute_default_intra_period(sequence_control_set_ptr);
    }

    if (sequence_control_set_ptr->static_config.look_ahead_distance == (uint32_t)~0)
#if SHUT_LOOKAHEAD
        sequence_control_set_ptr->static_config.look_ahead_distance = 0;
#else
        sequence_control_set_ptr->static_config.look_ahead_distance = compute_default_look_ahead(&sequence_control_set_ptr->static_config);
#endif
    else
        sequence_control_set_ptr->static_config.look_ahead_distance = cap_look_ahead_distance(&sequence_control_set_ptr->static_config);

#if ALTREF_FILTERING_SUPPORT
    sequence_control_set_ptr->static_config.enable_altrefs = pComponentParameterStructure->enable_altrefs;
    sequence_control_set_ptr->static_config.altref_strength = pComponentParameterStructure->altref_strength;
    sequence_control_set_ptr->static_config.altref_nframes = pComponentParameterStructure->altref_nframes;
#endif
#if ALT_REF_OVERLAY
    sequence_control_set_ptr->static_config.enable_overlays = pComponentParameterStructure->enable_overlays;
#endif

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

    for (i = 0; i < number_hme_search_region_in_width; i++) {
        total_search_width += number_hme_search_region_in_width_array[i];
    }
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

    for (i = 0; i < number_hme_search_region_in_width; i++) {
        total_search_width += number_hme_search_region_in_width_array[i];
    }
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

    if ((sequence_control_set_ptr->max_input_luma_height % 8) || (sequence_control_set_ptr->max_input_luma_width % 8)) {
        SVT_LOG("Error Instance %u: Only multiple of 8 resolutions are supported \n", channelNumber + 1);
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
        if (VerifyHmeDimention(channelNumber + 1, config->hme_level0_total_search_area_height, config->hme_level0_search_area_in_height_array, config->number_hme_search_region_in_height)) {
            return_error = EB_ErrorBadParameter;
        }

        if (VerifyHmeDimention(channelNumber + 1, config->hme_level0_total_search_area_width, config->hme_level0_search_area_in_width_array, config->number_hme_search_region_in_width)) {
            return_error = EB_ErrorBadParameter;
        }
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level1_search_area_in_width_array, config->number_hme_search_region_in_width)) {
            return_error = EB_ErrorBadParameter;
        }
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level1_search_area_in_height_array, config->number_hme_search_region_in_width)) {
            return_error = EB_ErrorBadParameter;
        }
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level2_search_area_in_width_array, config->number_hme_search_region_in_width)) {
            return_error = EB_ErrorBadParameter;
        }
        if (VerifyHmeDimentionL1L2(channelNumber + 1, config->hme_level2_search_area_in_height_array, config->number_hme_search_region_in_width)) {
            return_error = EB_ErrorBadParameter;
        }
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
#if RC
    if (config->rate_control_mode > 3) {
        SVT_LOG("Error Instance %u: The rate control mode must be [0 - 3] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
#else
    if (config->rate_control_mode > 1) {
        SVT_LOG("Error Instance %u: The rate control mode must be [0 - 1] \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
#endif
    if (config->rate_control_mode == 1) {
        SVT_LOG("Error Instance %u: The rate control mode 1 is currently not supported \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
#if RC
    if ((config->rate_control_mode == 3|| config->rate_control_mode == 2) && config->look_ahead_distance != (uint32_t)config->intra_period_length) {
        SVT_LOG("Error Instance %u: The rate control mode 2/3 LAD must be equal to intra_period \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
#endif
    if (config->look_ahead_distance > MAX_LAD && config->look_ahead_distance != (uint32_t)~0) {
        SVT_LOG("Error Instance %u: The lookahead distance must be [0 - %d] \n", channelNumber + 1, MAX_LAD);

        return_error = EB_ErrorBadParameter;
    }
    if (config->tile_rows < 0 || config->tile_columns < 0 || config->tile_rows > 6 || config->tile_columns > 6) {
        SVT_LOG("Error Instance %u: Log2Tile rows/cols must be [0 - 6] \n", channelNumber + 1);
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

    if (config->improve_sharpness > 1) {
        SVT_LOG("Error instance %u : Invalid ImproveSharpness. ImproveSharpness must be [0 - 1]\n", channelNumber + 1);
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

    if (((int32_t)(config->asm_type) < -1) || ((int32_t)(config->asm_type) != 1)) {
       // SVT_LOG("Error Instance %u: Invalid asm type value [0: C Only, 1: Auto] .\n", channelNumber + 1);
        SVT_LOG("Error Instance %u: Asm 0 is not supported in this build .\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->target_socket != -1 && config->target_socket != 0 && config->target_socket != 1) {
        SVT_LOG("Error instance %u: Invalid target_socket. target_socket must be [-1 - 1] \n", channelNumber + 1);
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
#if RC
    config_ptr->min_qp_allowed = 10;
#else
    config_ptr->min_qp_allowed = 0;
#endif
    config_ptr->base_layer_switch_mode = 0;
    config_ptr->enc_mode = MAX_ENC_PRESET;
    config_ptr->intra_period_length = -2;
    config_ptr->intra_refresh_type = 1;
    config_ptr->hierarchical_levels = 4;
    config_ptr->pred_structure = EB_PRED_RANDOM_ACCESS;
    config_ptr->disable_dlf_flag = EB_FALSE;
#if ENABLE_WARPED_MV
    config_ptr->enable_warped_motion = EB_TRUE;
#else
    config_ptr->enable_warped_motion = EB_FALSE;
#endif
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
    config_ptr->constrained_intra = EB_FALSE;
    config_ptr->improve_sharpness = EB_FALSE;

    // Bitstream options
    //config_ptr->codeVpsSpsPps = 0;
    //config_ptr->codeEosNal = 0;

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

    // ASM Type
    config_ptr->asm_type = 1;

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
    config_ptr->enable_overlays = EB_TRUE;

    return return_error;
}
//#define DEBUG_BUFFERS
static void print_lib_params(
    SequenceControlSet* scs) {

    EbSvtAv1EncConfiguration*   config = &scs->static_config;

    SVT_LOG("------------------------------------------- ");
    if (config->profile == 0)
        SVT_LOG("\nSVT [config]: Main Profile\t");
    else
        SVT_LOG("\nSVT [config]: Main10 Profile\t");

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
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate / LookaheadDistance / SceneChange\t\t: ABR / %d / %d / %d ", config->target_bit_rate, config->look_ahead_distance, config->scene_change_detection);
    else if (config->rate_control_mode == 2)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate / LookaheadDistance / SceneChange\t\t: VBR / %d / %d / %d ", config->target_bit_rate, config->look_ahead_distance, config->scene_change_detection);
    else if (config->rate_control_mode == 3)
        SVT_LOG("\nSVT [config]: RCMode / TargetBitrate / LookaheadDistance / SceneChange\t\t: Constraint VBR / %d / %d / %d ", config->target_bit_rate, config->look_ahead_distance, config->scene_change_detection);
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
#if defined(__linux__) || defined(__APPLE__)
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

    if (return_error == EB_ErrorBadParameter) {
        return EB_ErrorBadParameter;
    }

    SetParamBasedOnInput(
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

    // Initialize the Prediction Structure Group
    return_error = (EbErrorType)prediction_structure_group_ctor(
#if MRP_M1
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.enc_mode,
#endif
        &pEncCompData->sequence_control_set_instance_array[instance_index]->encode_context_ptr->prediction_structure_group_ptr,
        pEncCompData->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->static_config.base_layer_switch_mode);

    if (return_error == EB_ErrorInsufficientResources) {
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
#if defined(__linux__) || defined(__APPLE__)
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_enc_stream_header(
    EbComponentType           *svt_enc_component,
    EbBufferHeaderType        **output_stream_ptr){

    EbErrorType             return_error = EB_ErrorNone;
    UNUSED(svt_enc_component);
    UNUSED(output_stream_ptr);
    return return_error;
}
//
#if defined(__linux__) || defined(__APPLE__)
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
            lumaHeight,
            config->asm_type);

        un_pack2d(
            (uint16_t*)(inputPtr->cb + chromaOffset),
            sourceCbStride,
            input_picture_ptr->buffer_cb + chromaBufferOffset,
            input_picture_ptr->stride_cb,
            input_picture_ptr->buffer_bit_inc_cb + chromaBufferOffset,
            input_picture_ptr->stride_bit_inc_cb,
            chromaWidth,
            (lumaHeight >> 1),
            config->asm_type);

        un_pack2d(
            (uint16_t*)(inputPtr->cr + chromaOffset),
            sourceCrStride,
            input_picture_ptr->buffer_cr + chromaBufferOffset,
            input_picture_ptr->stride_cr,
            input_picture_ptr->buffer_bit_inc_cr + chromaBufferOffset,
            input_picture_ptr->stride_bit_inc_cr,
            chromaWidth,
            (lumaHeight >> 1),
            config->asm_type);
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
#if defined(__linux__) || defined(__APPLE__)
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
#if defined(__linux__) || defined(__APPLE__)
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
#if ALT_REF_OVERLAY
        if ( packet->flags & 0xfffffff0 )
            return_error = EB_ErrorMax;
#else
        if (packet->flags != EB_BUFFERFLAG_EOS &&
            packet->flags != EB_BUFFERFLAG_SHOW_EXT &&
            packet->flags != EB_BUFFERFLAG_HAS_TD &&
            packet->flags != (EB_BUFFERFLAG_SHOW_EXT | EB_BUFFERFLAG_EOS) &&
            packet->flags != (EB_BUFFERFLAG_SHOW_EXT | EB_BUFFERFLAG_HAS_TD) &&
            packet->flags != (EB_BUFFERFLAG_SHOW_EXT | EB_BUFFERFLAG_HAS_TD | EB_BUFFERFLAG_EOS) &&
            packet->flags != (EB_BUFFERFLAG_HAS_TD | EB_BUFFERFLAG_EOS) &&
            packet->flags != 0) {
            return_error = EB_ErrorMax;
        }
#endif
        // return the output stream buffer
        *p_buffer = packet;

        // save the wrapper pointer for the release
        (*p_buffer)->wrapper_ptr = (void*)ebWrapperPtr;
    }
    else {
        return_error = EB_NoErrorEmptyQueue;
    }

    return return_error;
}

#if defined(__linux__) || defined(__APPLE__)
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
#if defined(__linux__) || defined(__APPLE__)
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

            if (p_buffer->flags != EB_BUFFERFLAG_EOS && p_buffer->flags != 0) {
                return_error = EB_ErrorMax;
            }
            eb_release_object((EbObjectWrapper  *)ebWrapperPtr);
        }
        else {
            return_error = EB_NoErrorEmptyQueue;
        }
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

    // Encoder Private Handle Ctor
    return_error = (EbErrorType)eb_enc_handle_ctor(
        (EbEncHandle**) &(svt_enc_component->p_component_private),
        svt_enc_component);

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

    if (config->compressed_ten_bit_format == 1) {
        input_picture_buffer_desc_init_data.buffer_enable_mask = 0;
    }
    else {
        input_picture_buffer_desc_init_data.buffer_enable_mask = is16bit ? PICTURE_BUFFER_DESC_FULL_MASK : 0;
    }

    input_picture_buffer_desc_init_data.left_padding = sequence_control_set_ptr->left_padding;
    input_picture_buffer_desc_init_data.right_padding = sequence_control_set_ptr->right_padding;
    input_picture_buffer_desc_init_data.top_padding = sequence_control_set_ptr->top_padding;
    input_picture_buffer_desc_init_data.bot_padding = sequence_control_set_ptr->bot_padding;

    input_picture_buffer_desc_init_data.split_mode = is16bit ? EB_TRUE : EB_FALSE;

    input_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    if (is16bit && config->compressed_ten_bit_format == 1) {
        input_picture_buffer_desc_init_data.split_mode = EB_FALSE;  //do special allocation for 2bit data down below.
    }

    // Enhanced Picture Buffer
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(inputBuffer->p_buffer),
        (EbPtr)&input_picture_buffer_desc_init_data);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    if (is16bit && config->compressed_ten_bit_format == 1) {
        //pack 4 2bit pixels into 1Byte
        EB_ALLIGN_MALLOC(uint8_t*, ((EbPictureBufferDesc*)(inputBuffer->p_buffer))->buffer_bit_inc_y, sizeof(uint8_t) * (input_picture_buffer_desc_init_data.max_width / 4)*(input_picture_buffer_desc_init_data.max_height), EB_A_PTR);
        EB_ALLIGN_MALLOC(uint8_t*, ((EbPictureBufferDesc*)(inputBuffer->p_buffer))->buffer_bit_inc_cb, sizeof(uint8_t) * (input_picture_buffer_desc_init_data.max_width / 8)*(input_picture_buffer_desc_init_data.max_height / 2), EB_A_PTR);
        EB_ALLIGN_MALLOC(uint8_t*, ((EbPictureBufferDesc*)(inputBuffer->p_buffer))->buffer_bit_inc_cr, sizeof(uint8_t) * (input_picture_buffer_desc_init_data.max_width / 8)*(input_picture_buffer_desc_init_data.max_height / 2), EB_A_PTR);
    }

    return return_error;
}
/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbInputBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr)
{
    EbBufferHeaderType* inputBuffer;
    SequenceControlSet        *sequence_control_set_ptr = (SequenceControlSet*)objectInitDataPtr;
    EB_MALLOC(EbBufferHeaderType*, inputBuffer, sizeof(EbBufferHeaderType), EB_N_PTR);
    *objectDblPtr = (EbPtr)inputBuffer;
    // Initialize Header
    inputBuffer->size = sizeof(EbBufferHeaderType);

    allocate_frame_buffer(
        sequence_control_set_ptr,
        inputBuffer);

    inputBuffer->p_app_private = NULL;

    return EB_ErrorNone;
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbOutputBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr objectInitDataPtr)
{
    EbSvtAv1EncConfiguration   * config = (EbSvtAv1EncConfiguration*)objectInitDataPtr;
    uint32_t n_stride = (uint32_t)(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(config->source_width * config->source_height));  //TBC
    EbBufferHeaderType* outBufPtr;

    EB_MALLOC(EbBufferHeaderType*, outBufPtr, sizeof(EbBufferHeaderType), EB_N_PTR);
    *objectDblPtr = (EbPtr)outBufPtr;

    // Initialize Header
    outBufPtr->size = sizeof(EbBufferHeaderType);

    EB_MALLOC(uint8_t*, outBufPtr->p_buffer, n_stride, EB_N_PTR);

    outBufPtr->n_alloc_len = n_stride;
    outBufPtr->p_app_private = NULL;

    (void)objectInitDataPtr;

    return EB_ErrorNone;
}

/**************************************
* EbBufferHeaderType Constructor
**************************************/
EbErrorType EbOutputReconBufferHeaderCtor(
    EbPtr *objectDblPtr,
    EbPtr  objectInitDataPtr)
{
    EbBufferHeaderType         *recon_buffer;
    SequenceControlSet        *sequence_control_set_ptr = (SequenceControlSet*)objectInitDataPtr;
    const uint32_t luma_size =
        sequence_control_set_ptr->luma_width    *
        sequence_control_set_ptr->luma_height;
    // both u and v
    const uint32_t chroma_size = luma_size >> 1;
    const uint32_t tenBit = (sequence_control_set_ptr->static_config.encoder_bit_depth > 8);
    const uint32_t frameSize = (luma_size + chroma_size) << tenBit;

    EB_MALLOC(EbBufferHeaderType*, recon_buffer, sizeof(EbBufferHeaderType), EB_N_PTR);
    *objectDblPtr = (EbPtr)recon_buffer;

    // Initialize Header
    recon_buffer->size = sizeof(EbBufferHeaderType);

    // Assign the variables
    EB_MALLOC(uint8_t*, recon_buffer->p_buffer, frameSize, EB_N_PTR);

    recon_buffer->n_alloc_len = frameSize;
    recon_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}
