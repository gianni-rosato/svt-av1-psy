/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbAppConfig_h
#define EbAppConfig_h

#include <stdio.h>

#include "EbSvtAv1Enc.h"

typedef struct EbAppContext_ EbAppContext;

#ifdef _WIN32
#define fseeko _fseeki64
#define ftello _ftelli64
#endif
// Define Cross-Platform 64-bit fseek() and ftell()
/** The AppExitConditionType type is used to define the App main loop exit
conditions.
*/
typedef enum AppExitConditionType {
    APP_ExitConditionNone = 0,
    APP_ExitConditionFinished,
    APP_ExitConditionError
} AppExitConditionType;

/** The AppPortActiveType type is used to define the state of output ports in
the App.
*/
typedef enum AppPortActiveType { APP_PortActive = 0, APP_PortInactive } AppPortActiveType;

typedef enum EncodePass {
    ENCODE_SINGLE_PASS, //single pass mode
    ENCODE_FIRST_PASS,  // first pass of multi pass mode
    ENCODE_LAST_PASS,   // last pass of multi pass mode
    MAX_ENCODE_PASS = 2,
} EncodePass;

/** The EbPtr type is intended to be used to pass pointers to and from the svt
API.  This is a 32 bit pointer and is aligned on a 32 bit word boundary.
*/
typedef void *EbPtr;

#define WARNING_LENGTH 100

// memory map to be removed and replaced by malloc / free
typedef enum EbPtrType {
    EB_N_PTR     = 0, // malloc'd pointer
    EB_A_PTR     = 1, // malloc'd pointer aligned
    EB_MUTEX     = 2, // mutex
    EB_SEMAPHORE = 3, // semaphore
    EB_THREAD    = 4 // thread handle
} EbPtrType;
typedef struct EbMemoryMapEntry {
    EbPtr     ptr; // points to a memory pointer
    EbPtrType ptr_type; // pointer type
} EbMemoryMapEntry;

extern EbMemoryMapEntry *app_memory_map; // App Memory table
extern uint32_t *        app_memory_map_index; // App Memory index
extern uint64_t *        total_app_memory; // App Memory malloc'd
extern uint32_t          app_malloc_count;

#define MAX_APP_NUM_PTR (0x186A0 << 2) // Maximum number of pointers to be allocated for the app

#define EB_APP_MALLOC(type, pointer, n_elements, pointer_class, return_type) \
    pointer = (type)malloc(n_elements);                                      \
    if (pointer == (type)NULL) {                                             \
        return return_type;                                                  \
    } else {                                                                 \
        app_memory_map[*(app_memory_map_index)].ptr_type = pointer_class;    \
        app_memory_map[(*(app_memory_map_index))++].ptr  = pointer;          \
        if (n_elements % 8 == 0) {                                           \
            *total_app_memory += (n_elements);                               \
        } else {                                                             \
            *total_app_memory += ((n_elements) + (8 - ((n_elements) % 8)));  \
        }                                                                    \
    }                                                                        \
    if (*(app_memory_map_index) >= MAX_APP_NUM_PTR) { return return_type; }  \
    app_malloc_count++;

#define EB_APP_MALLOC_NR(type, pointer, n_elements, pointer_class, return_type) \
    (void)return_type;                                                          \
    pointer = (type)malloc(n_elements);                                         \
    if (pointer == (type)NULL) {                                                \
        return_type = EB_ErrorInsufficientResources;                            \
        fprintf(stderr, "Malloc has failed due to insuffucient resources");     \
        return;                                                                 \
    } else {                                                                    \
        app_memory_map[*(app_memory_map_index)].ptr_type = pointer_class;       \
        app_memory_map[(*(app_memory_map_index))++].ptr  = pointer;             \
        if (n_elements % 8 == 0) {                                              \
            *total_app_memory += (n_elements);                                  \
        } else {                                                                \
            *total_app_memory += ((n_elements) + (8 - ((n_elements) % 8)));     \
        }                                                                       \
    }                                                                           \
    if (*(app_memory_map_index) >= MAX_APP_NUM_PTR) {                           \
        return_type = EB_ErrorInsufficientResources;                            \
        fprintf(stderr, "Malloc has failed due to insuffucient resources");     \
        return;                                                                 \
    }                                                                           \
    app_malloc_count++;

#define EB_APP_MEMORY()                                                        \
    fprintf(stderr, "Total Number of Mallocs in App: %u\n", app_malloc_count); \
    fprintf(stderr, "Total App Memory: %.2lf KB\n\n", *total_app_memory / (double)1024);

#define MAX_CHANNEL_NUMBER 6U
#define MAX_NUM_TOKENS 210

#ifdef _WIN32
#define FOPEN(f, s, m) fopen_s(&f, s, m)
#else
#define FOPEN(f, s, m) f = fopen(s, m)
#endif

typedef struct EbPerformanceContext {
    /****************************************
     * Computational Performance Data
     ****************************************/
    uint64_t lib_start_time[2]; // [sec, micro_sec] including init time
    uint64_t encode_start_time[2]; // [sec, micro_sec] first frame sent

    double total_execution_time; // includes init
    double total_encode_time; // not including init

    uint64_t total_latency;
    uint32_t max_latency;

    uint64_t starts_time;
    uint64_t startu_time;
    uint64_t frame_count;

    double average_speed;
    double average_latency;

    uint64_t byte_count;
    double   sum_luma_psnr;
    double   sum_cr_psnr;
    double   sum_cb_psnr;

    double sum_luma_sse;
    double sum_cr_sse;
    double sum_cb_sse;

    double sum_luma_ssim;
    double sum_cr_ssim;
    double sum_cb_ssim;

    uint64_t sum_qp;

} EbPerformanceContext;

typedef struct EbConfig {
    /****************************************
     * File I/O
     ****************************************/
    FILE *        config_file;
    FILE *        input_file;
    EbBool        input_file_is_fifo;
    FILE *        bitstream_file;
    FILE *        recon_file;
    FILE *        error_log_file;
    FILE *        stat_file;
    FILE *        buffer_file;
    FILE *        qp_file;
    /* two pass */
    int           pass;
    const char*   stats;
    FILE *        input_stat_file;
    FILE *        output_stat_file;


    FILE *        input_pred_struct_file;
    char *        input_pred_struct_filename;
    EbBool        y4m_input;
    unsigned char y4m_buf[9];

    uint8_t       progress; // 0 = no progress output, 1 = normal, 2 = aomenc style verbose progress
    /****************************************
     * Computational Performance Data
     ****************************************/
    EbPerformanceContext performance_context;

    uint32_t input_padded_width;
    uint32_t input_padded_height;

    // -1 indicates unknown (auto-detect at earliest opportunity)
    // auto-detect is performed on load for files and at end of stream for pipes
    int64_t   frames_to_be_encoded;
    int32_t   frames_encoded;
    int32_t   buffered_input;
    uint8_t **sequence_buffer;

    uint32_t      injector_frame_rate;
    uint32_t      injector;
    uint32_t      speed_control_flag;


    uint32_t hme_level0_column_index;
    uint32_t hme_level0_row_index;
    uint32_t hme_level1_column_index;
    uint32_t hme_level1_row_index;
    uint32_t hme_level2_column_index;
    uint32_t hme_level2_row_index;
    EbBool   stop_encoder; // to signal CTRL+C Event, need to stop encoding.

    uint64_t processed_frame_count;
    uint64_t processed_byte_count;

    uint64_t byte_count_since_ivf;
    uint64_t ivf_count;
    /****************************************
     * On-the-fly Testing
     ****************************************/
    EbBool eos_flag;

    EbSvtAv1EncConfiguration config;
} EbConfig;

typedef struct EncChannel {
    EbConfig        *config; // Encoder Configuration
    EbAppContext    *app_callback; // Instances App callback date
    EbErrorType     return_error; // Error Handling
    AppExitConditionType exit_cond_output; // Processing loop exit condition
    AppExitConditionType exit_cond_recon; // Processing loop exit condition
    AppExitConditionType exit_cond_input; // Processing loop exit condition
    AppExitConditionType exit_cond; // Processing loop exit condition
    EbBool active;
} EncChannel;

typedef struct EncApp {
    SvtAv1FixedBuf rc_twopasses_stats;
} EncApp;

EbConfig * svt_config_ctor(EncodePass pass);
void svt_config_dtor(EbConfig *config_ptr);

EbErrorType enc_channel_ctor(EncChannel* c, EncodePass pass);
void enc_channel_dctor(EncChannel* c, uint32_t inst_cnt);

extern EbErrorType read_command_line(int32_t argc, char *const argv[], EncChannel *channels,
                                     uint32_t num_channels,
                                     char *warning_str[WARNING_LENGTH]);
extern uint32_t    get_help(int32_t argc, char *const argv[]);
extern uint32_t    get_number_of_channels(int32_t argc, char *const argv[]);
uint32_t get_passes(int32_t argc, char *const argv[], EncodePass pass[]);
EbErrorType set_two_passes_stats(EbConfig *config, EncodePass pass,
    const SvtAv1FixedBuf* rc_twopass_stats_in, uint32_t channel_number);
#endif //EbAppConfig_h
