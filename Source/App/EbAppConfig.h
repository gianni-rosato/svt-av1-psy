/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAppConfig_h
#define EbAppConfig_h

#include <stdio.h>

#include "EbApi.h"

#ifdef __GNUC__
#define fseeko64 fseek
#define ftello64 ftell
#endif
// Define Cross-Platform 64-bit fseek() and ftell()
#ifdef _MSC_VER
typedef __int64 off64_t;
#define fseeko64 _fseeki64
#define ftello64 _ftelli64

#elif _WIN32 // MinGW

#endif

#ifndef _RSIZE_T_DEFINED
typedef size_t rsize_t;
#define _RSIZE_T_DEFINED
#endif  /* _RSIZE_T_DEFINED */

#ifndef _ERRNO_T_DEFINED
#define _ERRNO_T_DEFINED
typedef int32_t errno_t;
#endif  /* _ERRNO_T_DEFINED */

/** The APPEXITCONDITIONTYPE type is used to define the App main loop exit
conditions.
*/
typedef enum APPEXITCONDITIONTYPE {
    APP_ExitConditionNone = 0,
    APP_ExitConditionFinished,
    APP_ExitConditionError
} APPEXITCONDITIONTYPE;

/** The APPPORTACTIVETYPE type is used to define the state of output ports in
the App.
*/
typedef enum APPPORTACTIVETYPE {
    APP_PortActive = 0,
    APP_PortInactive
} APPPORTACTIVETYPE;

//typedef enum EbPtrType {
//    EB_N_PTR = 0,                                   // malloc'd pointer
//    EB_A_PTR = 1,                                   // malloc'd pointer aligned
//    EB_MUTEX = 2,                                   // mutex
//    EB_SEMAPHORE = 3,                                   // semaphore
//    EB_THREAD = 4                                    // thread handle
//}EbPtrType;

/** The EB_PTR type is intended to be used to pass pointers to and from the svt
API.  This is a 32 bit pointer and is aligned on a 32 bit word boundary.
*/
typedef void * EB_PTR;

/** The EB_NULL type is used to define the C style NULL pointer.
*/
#define EB_NULL ((void*) 0)

// memory map to be removed and replaced by malloc / free
typedef enum EbPtrType {
    EB_N_PTR = 0,                                   // malloc'd pointer
    EB_A_PTR = 1,                                   // malloc'd pointer aligned
    EB_MUTEX = 2,                                   // mutex
    EB_SEMAPHORE = 3,                                   // semaphore
    EB_THREAD = 4                                    // thread handle
}EbPtrType;
typedef void * EbPtr;
typedef struct EbMemoryMapEntry
{
    EbPtr                     ptr;                       // points to a memory pointer
    EbPtrType                 ptrType;                   // pointer type
} EbMemoryMapEntry;

extern    EbMemoryMapEntry        *appMemoryMap;            // App Memory table
extern    uint32_t                  *appMemoryMapIndex;       // App Memory index
extern    uint64_t                  *totalAppMemory;          // App Memory malloc'd
extern    uint32_t                   appMallocCount;

#define MAX_APP_NUM_PTR                             (0x186A0 << 2)             // Maximum number of pointers to be allocated for the app

#define EB_APP_MALLOC(type, pointer, nElements, pointerClass, returnType) \
    pointer = (type)malloc(nElements); \
    if (pointer == (type)EB_NULL){ \
        return returnType; \
            } \
                else { \
        appMemoryMap[*(appMemoryMapIndex)].ptrType = pointerClass; \
        appMemoryMap[(*(appMemoryMapIndex))++].ptr = pointer; \
        if (nElements % 8 == 0) { \
            *totalAppMemory += (nElements); \
                        } \
                                else { \
            *totalAppMemory += ((nElements) + (8 - ((nElements) % 8))); \
            } \
        } \
    if (*(appMemoryMapIndex) >= MAX_APP_NUM_PTR) { \
        return returnType; \
                } \
    appMallocCount++;

#define EB_APP_MALLOC_NR(type, pointer, nElements, pointerClass,returnType) \
    (void)returnType; \
    pointer = (type)malloc(nElements); \
    if (pointer == (type)EB_NULL){ \
        returnType = EB_ErrorInsufficientResources; \
        printf("Malloc has failed due to insuffucient resources"); \
        return; \
            } \
                else { \
        appMemoryMap[*(appMemoryMapIndex)].ptrType = pointerClass; \
        appMemoryMap[(*(appMemoryMapIndex))++].ptr = pointer; \
        if (nElements % 8 == 0) { \
            *totalAppMemory += (nElements); \
                        } \
                                else { \
            *totalAppMemory += ((nElements) + (8 - ((nElements) % 8))); \
            } \
        } \
    if (*(appMemoryMapIndex) >= MAX_APP_NUM_PTR) { \
        returnType = EB_ErrorInsufficientResources; \
        printf("Malloc has failed due to insuffucient resources"); \
        return; \
                } \
    appMallocCount++;

/* string copy */
extern errno_t strcpy_ss(char *dest, rsize_t dmax, const char *src);

/* fitted string copy */
extern errno_t strncpy_ss(char *dest, rsize_t dmax, const char *src, rsize_t slen);

/* string length */
extern rsize_t strnlen_ss(const char *s, rsize_t smax);

#define EB_STRNCPY(dst, src, count) \
    strncpy_ss(dst, sizeof(dst), src, count)

#define EB_STRCPY(dst, size, src) \
    strcpy_ss(dst, size, src)

#define EB_STRCMP(target,token) \
    strcmp(target,token)

#define EB_STRLEN(target, max_size) \
    strnlen_ss(target, max_size)

#define EB_APP_MEMORY() \
    printf("Total Number of Mallocs in App: %d\n", appMallocCount); \
    printf("Total App Memory: %.2lf KB\n\n",*totalAppMemory/(double)1024);

#define MAX_CHANNEL_NUMBER      6
#define MAX_NUM_TOKENS          200

#ifdef _MSC_VER
#define FOPEN(f,s,m) fopen_s(&f,s,m)
#else
#define FOPEN(f,s,m) f=fopen(s,m)
#endif

/****************************************
* Padding
****************************************/
#define LEFT_INPUT_PADDING 0
#define RIGHT_INPUT_PADDING 0
#define TOP_INPUT_PADDING 0
#define BOTTOM_INPUT_PADDING 0


typedef struct EbPerformanceContext_s {

    /****************************************
     * Computational Performance Data
     ****************************************/
    uint64_t                  lib_start_time[2];       // [sec, micro_sec] including init time
    uint64_t                  encode_start_time[2];    // [sec, micro_sec] first frame sent

    double                    total_execution_time;    // includes init
    double                    total_encode_time;       // not including init

    uint64_t                  totalLatency;
    uint32_t                  maxLatency;

    uint64_t                  startsTime;
    uint64_t                  startuTime;
    uint64_t                  frameCount;

    double                    averageSpeed;
    double                    averageLatency;

    uint64_t                  byteCount;

}EbPerformanceContext_t;

typedef struct EbConfig_s
{
    /****************************************
     * File I/O
     ****************************************/
    FILE                    *configFile;
    FILE                    *inputFile;
    FILE                    *bitstreamFile;
    FILE                    *reconFile;
    FILE                    *errorLogFile;
    FILE                    *bufferFile;

    FILE                    *qpFile;

    EbBool                  y4mInput;
    unsigned char           y4mBuf[9];

    EbBool                  use_qp_file;

    uint32_t                 frameRate;
    uint32_t                 frameRateNumerator;
    uint32_t                 frameRateDenominator;
    uint32_t                 injector_frame_rate;
    uint32_t                 injector;
    uint32_t                 speed_control_flag;
    uint32_t                 encoderBitDepth;
    uint32_t                 compressedTenBitFormat;
    uint32_t                 sourceWidth;
    uint32_t                 sourceHeight;

    uint32_t                 inputPaddedWidth;
    uint32_t                 inputPaddedHeight;

    int64_t                  framesToBeEncoded;
    int32_t                  framesEncoded;
    int32_t                  bufferedInput;
    uint8_t                **sequenceBuffer;

    uint8_t                  latencyMode;

    /****************************************
     * // Interlaced Video
     ****************************************/
    EbBool                  interlacedVideo;
    EbBool                  separateFields;

    /*****************************************
     * Coding Structure
     *****************************************/
    uint32_t                 base_layer_switch_mode;
    uint8_t                  encMode;
    int32_t                  intraPeriod;
    uint32_t                 intraRefreshType;
    uint32_t                 hierarchicalLevels;
    uint32_t                 predStructure;


    /****************************************
     * Quantization
     ****************************************/
    uint32_t                 qp;

    /****************************************
     * Film Grain
     ****************************************/
    EbBool                  film_grain_denoise_strength;
     /****************************************
     * DLF
     ****************************************/
    EbBool                  disable_dlf_flag;

    /****************************************
     * Local Warped Motion
     ****************************************/
    EbBool                  enable_warped_motion;

    /****************************************
     * ME Tools
     ****************************************/
    EbBool                  use_default_me_hme;
    EbBool                  enableHmeFlag;
    EbBool                  enableHmeLevel0Flag;
    EbBool                  enableHmeLevel1Flag;
    EbBool                  enableHmeLevel2Flag;
    EbBool                  ext_block_flag;
    EbBool                  in_loop_me_flag;

    /****************************************
     * ME Parameters
     ****************************************/
    uint32_t                 searchAreaWidth;
    uint32_t                 searchAreaHeight;

    /****************************************
     * HME Parameters
     ****************************************/
    uint32_t                 numberHmeSearchRegionInWidth ;
    uint32_t                 numberHmeSearchRegionInHeight;
    uint32_t                 hmeLevel0TotalSearchAreaWidth;
    uint32_t                 hmeLevel0TotalSearchAreaHeight;
    uint32_t                 hmeLevel0ColumnIndex;
    uint32_t                 hmeLevel0RowIndex;
    uint32_t                 hmeLevel1ColumnIndex;
    uint32_t                 hmeLevel1RowIndex;
    uint32_t                 hmeLevel2ColumnIndex;
    uint32_t                 hmeLevel2RowIndex;
    uint32_t                 hmeLevel0SearchAreaInWidthArray[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t                 hmeLevel0SearchAreaInHeightArray[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t                 hmeLevel1SearchAreaInWidthArray[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t                 hmeLevel1SearchAreaInHeightArray[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t                 hmeLevel2SearchAreaInWidthArray[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t                 hmeLevel2SearchAreaInHeightArray[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];

    /****************************************
     * MD Parameters
     ****************************************/
    EbBool                  constrained_intra;

#if TILES
    int32_t                  tile_columns;
    int32_t                  tile_rows;
#endif

    /****************************************
     * Rate Control
     ****************************************/
    uint32_t                 scene_change_detection;
    uint32_t                 rateControlMode;
    uint32_t                 look_ahead_distance;
    uint32_t                 targetBitRate;
    uint32_t                 max_qp_allowed;
    uint32_t                 min_qp_allowed;

    /****************************************
     * Optional Features
     ****************************************/

    EbBool                   improve_sharpness;
    uint32_t                 high_dynamic_range_input;
    uint32_t                 access_unit_delimiter;
    uint32_t                 buffering_period_sei;
    uint32_t                 picture_timing_sei;
    EbBool                   registered_user_data_sei_flag;
    EbBool                   unregistered_user_data_sei_flag;
    EbBool                   recovery_point_sei_flag;
    uint32_t                 enable_temporal_id;

    /****************************************
     * Annex A Parameters
     ****************************************/
    uint32_t                 profile;
    uint32_t                 tier;
    uint32_t                 level;

    /****************************************
     * On-the-fly Testing
     ****************************************/
    uint32_t                 testUserData;
    EbBool                   eosFlag;

    /****************************************
    * Optimization Type
    ****************************************/
    uint32_t                  asmType;

    /****************************************
     * Computational Performance Data
     ****************************************/
    EbPerformanceContext_t  performanceContext;

    /****************************************
    * Instance Info
    ****************************************/
    uint32_t                channel_id;
    uint32_t                active_channel_count;
    uint32_t                logicalProcessors;
    int32_t                 targetSocket;
    EbBool                 stopEncoder;         // to signal CTRL+C Event, need to stop encoding.

    uint64_t                processedFrameCount;
    uint64_t                processedByteCount;

    uint64_t                byte_count_since_ivf;
    uint64_t                ivf_count;

} EbConfig_t;

extern void EbConfigCtor(EbConfig_t *config_ptr);
extern void EbConfigDtor(EbConfig_t *config_ptr);

extern EbErrorType    ReadCommandLine(int32_t argc, char *const argv[], EbConfig_t **config, uint32_t  numChannels,    EbErrorType *return_errors);
extern uint32_t     GetHelp(int32_t argc, char *const argv[]);
extern uint32_t        GetNumberOfChannels(int32_t argc, char *const argv[]);

#endif //EbAppConfig_h
