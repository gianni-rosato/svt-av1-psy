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

// main.cpp
//  -Constructs the following resources needed during the encoding process
//      -memory
//      -threads
//      -semaphores
//  -Configures the encoder
//  -Calls the encoder via the API
//  -Destructs the resources

/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>
#include "EbAppConfig.h"
#include "EbAppContext.h"
#include "EbTime.h"
#include <fcntl.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h> /* _setmode() */
#else
#include <pthread.h>
#include <semaphore.h>
#include <errno.h>
#endif

#if !defined(_WIN32) || !defined(HAVE_STRNLEN_S)
#include "third_party/safestringlib/safe_str_lib.h"
#endif

#ifdef __GLIBC__
#include <malloc.h>
#endif

#ifndef _WIN32
#include <unistd.h>
#endif

#if LOG_ENC_DONE
int tot_frames_done = 0;
#endif
/***************************************
 * External Functions
 ***************************************/
void process_input_buffer(EncChannel* c);

void process_output_recon_buffer(EncChannel* c);

void process_output_stream_buffer(EncChannel* c, EncApp* enc_app, int32_t* frame_count);

void init_reader(EbConfig* app_cfg);

volatile int32_t keep_running = 1;

void event_handler(int32_t dummy) {
    (void)dummy;
    keep_running = 0;

    // restore default signal handler
    signal(SIGINT, SIG_DFL);
}

void assign_app_thread_group(uint8_t target_socket) {
#ifdef _WIN32
    if (GetActiveProcessorGroupCount() == 2) {
        GROUP_AFFINITY svt_aom_group_affinity;
        GetThreadGroupAffinity(GetCurrentThread(), &svt_aom_group_affinity);
        svt_aom_group_affinity.Group = target_socket;
        SetThreadGroupAffinity(GetCurrentThread(), &svt_aom_group_affinity, NULL);
    }
#else
    (void)target_socket;
    return;
#endif
}

typedef struct EncContext {
    uint32_t   num_channels;
    EncChannel channels[MAX_CHANNEL_NUMBER];
    char*      warning[MAX_NUM_TOKENS];
    EncPass    enc_pass;
    int32_t    passes;
    int32_t    total_frames;
} EncContext;

//initilize memory mapped file handler
static void init_memory_file_map(EbConfig* app_cfg) {
    app_cfg->mmap.enable = app_cfg->buffered_input == -1 && !app_cfg->input_file_is_fifo;

    if (!app_cfg->mmap.enable)
        return;
    // app_cfg->mmap.y4m_seq_hdr   = 0; // already initialized to 0 or some value in read_y4m_header()
    app_cfg->mmap.y4m_frm_hdr   = 0;
    app_cfg->mmap.file_frame_it = 0;
    // The cur_offset shows the current offset to be read. Initially set to 0 or y4m hdr
    // After each frame, the offset is updated. The changes were made to support changing the resolution on the fly
    app_cfg->mmap.cur_offset = app_cfg->mmap.y4m_seq_hdr;
    const int64_t curr_loc   = ftello(app_cfg->input_file); // get current fp location
    fseeko(app_cfg->input_file, 0L, SEEK_END); // seek to end of file
    app_cfg->mmap.file_size = ftello(app_cfg->input_file); // get file size
    fseeko(app_cfg->input_file, curr_loc, SEEK_SET); // seek back to that location
#ifdef _WIN32
    SYSTEM_INFO sys_info;
    GetSystemInfo(&sys_info);
    app_cfg->mmap.fd         = _fileno(app_cfg->input_file);
    app_cfg->mmap.align_mask = sys_info.dwAllocationGranularity > sys_info.dwPageSize
        ? sys_info.dwAllocationGranularity - 1
        : sys_info.dwPageSize - 1;
    HANDLE fhandle           = (HANDLE)_get_osfhandle(app_cfg->mmap.fd);
    app_cfg->mmap.map_handle = CreateFileMapping(fhandle, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!app_cfg->mmap.map_handle) {
        app_cfg->mmap.enable = false;
        return;
    }
#else
    app_cfg->mmap.fd = fileno(app_cfg->input_file);
#if defined(__linux__)
    posix_fadvise(app_cfg->mmap.fd, 0, 0, POSIX_FADV_NOREUSE);
#endif
    app_cfg->mmap.align_mask = sysconf(_SC_PAGESIZE) - 1;
#endif
}

static void deinit_memory_file_map(EbConfig* app_cfg) {
    if (!app_cfg->mmap.enable)
        return;
#ifdef _WIN32
    CloseHandle(app_cfg->mmap.map_handle);
#endif
}

static int compar_uint64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static EbErrorType enc_context_ctor(EncApp* enc_app, EncContext* enc_context, int32_t argc, char* argv[],
                                    EncPass enc_pass, int32_t passes) {
#if LOG_ENC_DONE
    tot_frames_done = 0;
#endif

    memset(enc_context, 0, sizeof(*enc_context));
    uint32_t num_channels = get_number_of_channels(argc, argv);
    if (num_channels == 0)
        return EB_ErrorBadParameter;
    enc_context->enc_pass    = enc_pass;
    enc_context->passes      = passes;
    EbErrorType return_error = EB_ErrorNone;

    enc_context->num_channels = num_channels;
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        return_error = enc_channel_ctor(enc_context->channels + inst_cnt);
        if (return_error != EB_ErrorNone)
            return return_error;
    }

    char** warning = enc_context->warning;

    for (int token_id = 0; token_id < MAX_NUM_TOKENS; token_id++) {
        warning[token_id] = (char*)malloc(WARNING_LENGTH);
        if (!warning[token_id]) {
            return EB_ErrorInsufficientResources;
        }
        strcpy_s(warning[token_id], WARNING_LENGTH, "");
    }
    // Process any command line options, including the configuration file
    // Read all configuration files.
    return_error = read_command_line(argc, argv, enc_context->channels, num_channels);
    if (return_error != EB_ErrorNone) {
        fprintf(stderr, "Error in configuration, could not begin encoding! ... \n");
        fprintf(stderr, "Run %s --help for a list of options\n", argv[0]);
        return return_error;
    }
    // Set main thread affinity
    if (enc_context->channels[0].app_cfg->config.target_socket != -1)
        assign_app_thread_group(enc_context->channels[0].app_cfg->config.target_socket);

    // Init the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        EncChannel* c = enc_context->channels + inst_cnt;
        if (c->return_error == EB_ErrorNone) {
            EbConfig* app_cfg                    = c->app_cfg;
            app_cfg->config.active_channel_count = num_channels;
            app_cfg->config.channel_id           = inst_cnt;
            app_cfg->config.recon_enabled        = app_cfg->recon_file ? TRUE : FALSE;

            // set force_key_frames frames
            if (app_cfg->config.force_key_frames) {
                const double fps = (double)app_cfg->config.frame_rate_numerator /
                    app_cfg->config.frame_rate_denominator;
                struct forced_key_frames* forced_keyframes = &app_cfg->forced_keyframes;

                for (size_t i = 0; i < forced_keyframes->count; ++i) {
                    char*  p;
                    double val = strtod(forced_keyframes->specifiers[i], &p);
                    switch (*p) {
                    case 'f':
                    case 'F': break;
                    case 's':
                    case 'S':
                    default: val *= fps; break;
                    }
                    forced_keyframes->frames[i] = (uint64_t)val;
                }
                qsort(forced_keyframes->frames,
                      forced_keyframes->count,
                      sizeof(forced_keyframes->frames[0]),
                      compar_uint64);
            }
            init_memory_file_map(app_cfg);
            init_reader(app_cfg);

            app_svt_av1_get_time(&app_cfg->performance_context.lib_start_time[0],
                                 &app_cfg->performance_context.lib_start_time[1]);
            // Update pass
            app_cfg->config.pass = passes == 1 ? app_cfg->config.pass // Single-Pass
                                               : (int)enc_pass; // Multi-Pass

            c->return_error = handle_stats_file(app_cfg, enc_pass, &enc_app->rc_twopasses_stats, num_channels);
            if (c->return_error == EB_ErrorNone) {
                c->return_error = init_encoder(app_cfg, inst_cnt);
            }
            return_error = (EbErrorType)(return_error | c->return_error);
        } else
            c->active = FALSE;
    }
    return return_error;
}

static void enc_context_dctor(EncContext* enc_context) {
    // DeInit Encoder
    for (int32_t inst_cnt = enc_context->num_channels - 1; inst_cnt >= 0; --inst_cnt) {
        EncChannel* c = enc_context->channels + inst_cnt;
        deinit_memory_file_map(c->app_cfg);
        enc_channel_dctor(c, inst_cnt);
    }

    for (uint32_t warning_id = 0; warning_id < MAX_NUM_TOKENS; warning_id++) free(enc_context->warning[warning_id]);
}

double      get_psnr(double sse, double max);
static void print_summary(const EncContext* const enc_context) {
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        const EncChannel* const c       = &enc_context->channels[inst_cnt];
        const EbConfig*         app_cfg = c->app_cfg;
        if (c->exit_cond == APP_ExitConditionFinished && c->return_error == EB_ErrorNone &&
            (app_cfg->config.pass == 0 || app_cfg->config.pass == 2)) {
#if LOG_ENC_DONE
            tot_frames_done = (int)app_cfg->performance_context.frame_count;
#endif

            uint64_t frame_count    = (uint32_t)app_cfg->performance_context.frame_count;
            uint32_t max_luma_value = (app_cfg->config.encoder_bit_depth == 8) ? 255 : 1023;
            double   max_luma_sse   = (double)max_luma_value * max_luma_value *
                (app_cfg->config.source_width * app_cfg->config.source_height);
            double max_chroma_sse = (double)max_luma_value * max_luma_value *
                (app_cfg->config.source_width / 2 * app_cfg->config.source_height / 2);

            const double frame_rate = (double)app_cfg->config.frame_rate_numerator /
                (double)app_cfg->config.frame_rate_denominator;

            if (app_cfg->config.stat_report && app_cfg->stat_file) {
                fprintf(app_cfg->stat_file,
                        "\nSUMMARY "
                        "------------------------------------------------------"
                        "---------------\n");
                fprintf(app_cfg->stat_file,
                        "\n\t\t\t\tAverage PSNR (using per-frame "
                        "PSNR)\t\t|\tOverall PSNR (using per-frame MSE)\t\t|"
                        "\tAverage SSIM\n");
                fprintf(app_cfg->stat_file,
                        "Total Frames\tAverage QP  \tY-PSNR   \tU-PSNR   "
                        "\tV-PSNR\t\t| \tY-PSNR   \tU-PSNR   \tV-PSNR   \t|"
                        "\tY-SSIM   \tU-SSIM   \tV-SSIM   "
                        "\t|\tBitrate\n");
                fprintf(app_cfg->stat_file,
                        "%10ld  \t   %2.2f    \t%3.2f dB\t%3.2f dB\t%3.2f dB  "
                        "\t|\t%3.2f dB\t%3.2f dB\t%3.2f dB \t|\t%1.5f \t%1.5f "
                        "\t%1.5f\t\t|\t%.2f kbps\n",
                        (long int)frame_count,
                        (float)app_cfg->performance_context.sum_qp / frame_count,
                        (float)app_cfg->performance_context.sum_luma_psnr / frame_count,
                        (float)app_cfg->performance_context.sum_cb_psnr / frame_count,
                        (float)app_cfg->performance_context.sum_cr_psnr / frame_count,
                        (float)(get_psnr((app_cfg->performance_context.sum_luma_sse / frame_count), max_luma_sse)),
                        (float)(get_psnr((app_cfg->performance_context.sum_cb_sse / frame_count), max_chroma_sse)),
                        (float)(get_psnr((app_cfg->performance_context.sum_cr_sse / frame_count), max_chroma_sse)),
                        (float)app_cfg->performance_context.sum_luma_ssim / frame_count,
                        (float)app_cfg->performance_context.sum_cb_ssim / frame_count,
                        (float)app_cfg->performance_context.sum_cr_ssim / frame_count,
                        ((double)(app_cfg->performance_context.byte_count << 3) * frame_rate /
                         (app_cfg->frames_encoded * 1000)));
            }

            fprintf(stderr,
                    "\nSUMMARY --------------------------------- Channel %u  "
                    "--------------------------------\n",
                    inst_cnt + 1);
            fprintf(stderr, "Total Frames\t\tFrame Rate\t\tByte Count\t\tBitrate\n");
            fprintf(stderr,
                    "%12d\t\t%4.2f fps\t\t%10.0f\t\t%5.2f kbps\n",
                    (int32_t)frame_count,
                    frame_rate,
                    (double)app_cfg->performance_context.byte_count,
                    ((double)(app_cfg->performance_context.byte_count << 3) * frame_rate /
                     (app_cfg->frames_encoded * 1000)));

            if (app_cfg->config.stat_report) {
                fprintf(stderr,
                        "\n\t\tAverage PSNR (using per-frame "
                        "PSNR)\t\t|\tOverall PSNR (using per-frame MSE)\t\t|\t"
                        "Average SSIM\n");
                fprintf(stderr,
                        "Average "
                        "QP\tY-PSNR\t\tU-PSNR\t\tV-PSNR\t\t|\tY-PSNR\t\tU-"
                        "PSNR\t\tV-PSNR\t\t|\tY-SSIM\tU-SSIM\tV-SSIM\n");
                fprintf(stderr,
                        "%11.2f\t%4.2f dB\t%4.2f dB\t%4.2f dB\t|\t%4.2f "
                        "dB\t%4.2f dB\t%4.2f dB\t|\t%1.5f\t%1.5f\t%1.5f\n",
                        (float)app_cfg->performance_context.sum_qp / frame_count,
                        (float)app_cfg->performance_context.sum_luma_psnr / frame_count,
                        (float)app_cfg->performance_context.sum_cb_psnr / frame_count,
                        (float)app_cfg->performance_context.sum_cr_psnr / frame_count,
                        (float)(get_psnr((app_cfg->performance_context.sum_luma_sse / frame_count), max_luma_sse)),
                        (float)(get_psnr((app_cfg->performance_context.sum_cb_sse / frame_count), max_chroma_sse)),
                        (float)(get_psnr((app_cfg->performance_context.sum_cr_sse / frame_count), max_chroma_sse)),
                        (float)app_cfg->performance_context.sum_luma_ssim / frame_count,
                        (float)app_cfg->performance_context.sum_cb_ssim / frame_count,
                        (float)app_cfg->performance_context.sum_cr_ssim / frame_count);
            }

            fflush(stdout);
        }
    }
    fprintf(stderr, "\n");
    fflush(stdout);
}

static void print_performance(const EncContext* const enc_context) {
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        const EncChannel* c = enc_context->channels + inst_cnt;
        if (c->exit_cond == APP_ExitConditionFinished && c->return_error == EB_ErrorNone) {
            EbConfig* app_cfg = c->app_cfg;
            if (app_cfg->stop_encoder == FALSE) {
                if ((app_cfg->config.pass == 0 ||
                     (app_cfg->config.pass == 2 && app_cfg->config.rate_control_mode == SVT_AV1_RC_MODE_CQP_OR_CRF) ||
                     app_cfg->config.pass == 3))
                    fprintf(stderr,
                            "\nChannel %u\nAverage Speed:\t\t%.3f fps\nTotal Encoding Time:\t%.0f "
                            "ms\nTotal Execution Time:\t%.0f ms\nAverage Latency:\t%.0f ms\nMax "
                            "Latency:\t\t%u ms\n",
                            (uint32_t)(inst_cnt + 1),
                            app_cfg->performance_context.average_speed,
                            app_cfg->performance_context.total_encode_time * 1000,
                            app_cfg->performance_context.total_execution_time * 1000,
                            app_cfg->performance_context.average_latency,
                            (uint32_t)(app_cfg->performance_context.max_latency));
            } else
                fprintf(stderr, "\nChannel %u Encoding Interrupted\n", (uint32_t)(inst_cnt + 1));
        } else if (c->return_error == EB_ErrorInsufficientResources)
            fprintf(stderr, "Could not allocate enough memory for channel %u\n", inst_cnt + 1);
        else
            fprintf(stderr,
                    "Error encoding at channel %u! Check error log file for more details "
                    "... \n",
                    inst_cnt + 1);
    }
}

static void print_warnnings(const EncContext* const enc_context) {
    char* const* warning = enc_context->warning;
    for (uint32_t warning_id = 0;; warning_id++) {
        if (*warning[warning_id] == '-')
            fprintf(stderr, "warning: %s\n", warning[warning_id]);
        else if (*warning[warning_id + 1] != '-')
            break;
    }
}

static Bool is_active(const EncChannel* c) { return c->active; }

static Bool has_active_channel(const EncContext* const enc_context) {
    // check if all channels are inactive
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        if (is_active(enc_context->channels + inst_cnt))
            return TRUE;
    }
    return FALSE;
}
bool process_skip(EbConfig* app_cfg, EbBufferHeaderType* header_ptr);

static void enc_channel_step(EncChannel* c, EncApp* enc_app, EncContext* enc_context) {
    EbConfig* app_cfg = c->app_cfg;

    if (app_cfg->need_to_skip) {
        bool skip   = !process_skip(app_cfg, app_cfg->input_buffer_pool);
        int  next_c = fgetc(app_cfg->input_file);
        if (!skip && next_c == EOF) {
            fputs("\n[SVT-Error]: Skipped all available frames!\n", stderr);
            c->exit_cond_input = APP_ExitConditionFinished;
            c->active          = FALSE;
            return;
        }
        ungetc(next_c, app_cfg->input_file);
    }

    process_input_buffer(c);
    process_output_recon_buffer(c);
    process_output_stream_buffer(c, enc_app, &enc_context->total_frames);

    if (((c->exit_cond_recon == APP_ExitConditionFinished || !app_cfg->recon_file) &&
         c->exit_cond_output == APP_ExitConditionFinished && c->exit_cond_input == APP_ExitConditionFinished) ||
        ((c->exit_cond_recon == APP_ExitConditionError && app_cfg->recon_file) ||
         c->exit_cond_output == APP_ExitConditionError || c->exit_cond_input == APP_ExitConditionError)) {
        c->active = FALSE;
        if (app_cfg->recon_file)
            c->exit_cond = (AppExitConditionType)(c->exit_cond_recon | c->exit_cond_output | c->exit_cond_input);
        else
            c->exit_cond = (AppExitConditionType)(c->exit_cond_output | c->exit_cond_input);
    }
}
static const char* get_pass_name(EncPass enc_pass) {
    switch (enc_pass) {
    case ENC_FIRST_PASS: return "Pass 1/2 ";
    case ENC_SECOND_PASS: return "Pass 2/2 ";
    default: return "";
    }
}
static void enc_channel_start(EncChannel* c) {
    if (c->return_error == EB_ErrorNone) {
        EbConfig* app_cfg   = c->app_cfg;
        c->exit_cond        = APP_ExitConditionNone;
        c->exit_cond_output = APP_ExitConditionNone;
        c->exit_cond_recon  = app_cfg->recon_file ? APP_ExitConditionNone : APP_ExitConditionError;
        c->exit_cond_input  = APP_ExitConditionNone;
        c->active           = TRUE;
        app_svt_av1_get_time(&app_cfg->performance_context.encode_start_time[0],
                             &app_cfg->performance_context.encode_start_time[1]);
    }
}

static EbErrorType encode(EncApp* enc_app, EncContext* enc_context) {
    EbErrorType return_error = EB_ErrorNone;

    // Get num_channels
    uint32_t num_channels = enc_context->num_channels;
    EncPass  enc_pass     = enc_context->enc_pass;
    // Start the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt)
        enc_channel_start(enc_context->channels + inst_cnt);
    print_warnnings(enc_context);
    fprintf(stderr, "%sEncoding          ", get_pass_name(enc_pass));

    while (has_active_channel(enc_context)) {
        for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
            EncChannel* c = enc_context->channels + inst_cnt;
            if (is_active(c)) {
                enc_channel_step(c, enc_app, enc_context);
            }
        }
    }
    print_summary(enc_context);
    print_performance(enc_context);
    return return_error;
}

EbErrorType enc_app_ctor(EncApp* enc_app) {
    memset(enc_app, 0, sizeof(*enc_app));
    return EB_ErrorNone;
}

void enc_app_dctor(EncApp* enc_app) { free(enc_app->rc_twopasses_stats.buf); }

/***************************************
 * Encoder App Main
 ***************************************/
int32_t main(int32_t argc, char* argv[]) {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // GLOBAL VARIABLES
    EbErrorType return_error = EB_ErrorNone; // Error Handling
    uint32_t    passes;
    EncPass     enc_pass[MAX_ENC_PASS];
    EncApp      enc_app;
    EncContext  enc_context;

    signal(SIGINT, event_handler);
    if (get_version(argc, argv))
        return 0;

    if (get_help(argc, argv))
        return 0;

    if (get_color_help(argc, argv))
        return 0;

    enc_app_ctor(&enc_app);
    passes = get_passes(argc, argv, enc_pass);
    for (uint8_t pass_idx = 0; pass_idx < passes; pass_idx++) {
        return_error = enc_context_ctor(&enc_app, &enc_context, argc, argv, enc_pass[pass_idx], passes);

        if (return_error == EB_ErrorNone)
            return_error = encode(&enc_app, &enc_context);

        enc_context_dctor(&enc_context);
        if (return_error != EB_ErrorNone)
            break;

#ifdef __GLIBC__
        malloc_trim(0);
#endif
    }
    enc_app_dctor(&enc_app);

#if LOG_ENC_DONE
    fprintf(stderr, "all_done_encoding  %i frames \n", tot_frames_done);
#endif
    return return_error != EB_ErrorNone;
}
