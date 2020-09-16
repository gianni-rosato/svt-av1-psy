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
#ifdef _WIN32
#include <windows.h>
#include <io.h> /* _setmode() */
#include <fcntl.h> /* _O_BINARY */
#else
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#include <errno.h>
#endif

#if !defined(_WIN32) || !defined(HAVE_STRNLEN_S)
#include "safe_str_lib.h"
#endif

/***************************************
 * External Functions
 ***************************************/
extern AppExitConditionType process_input_buffer(EbConfig *config, EbAppContext *app_call_back);

extern AppExitConditionType process_output_recon_buffer(EbConfig *    config,
                                                        EbAppContext *app_call_back);

extern AppExitConditionType process_output_stream_buffer(EncApp* enc_app, EbConfig *    config,
                                                         EbAppContext *app_call_back,
                                                         uint8_t       pic_send_done,
                                                         int32_t      *frame_count);

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
        GROUP_AFFINITY group_affinity;
        GetThreadGroupAffinity(GetCurrentThread(), &group_affinity);
        group_affinity.Group = target_socket;
        SetThreadGroupAffinity(GetCurrentThread(), &group_affinity, NULL);
    }
#else
    (void)target_socket;
    return;
#endif
}

typedef struct EncContext {
    uint32_t        num_channels;
    EbConfig        *configs[MAX_CHANNEL_NUMBER]; // Encoder Configuration
    EbAppContext    *app_callbacks[MAX_CHANNEL_NUMBER]; // Instances App callback date
    char            *warning[MAX_NUM_TOKENS];
    EbErrorType     return_errors[MAX_CHANNEL_NUMBER]; // Error Handling
    AppExitConditionType exit_cond_output[MAX_CHANNEL_NUMBER]; // Processing loop exit condition
    AppExitConditionType exit_cond_recon[MAX_CHANNEL_NUMBER]; // Processing loop exit condition
    AppExitConditionType exit_cond_input[MAX_CHANNEL_NUMBER]; // Processing loop exit condition

    AppExitConditionType exit_cond[MAX_CHANNEL_NUMBER]; // Processing loop exit condition
    EbBool          channel_active[MAX_CHANNEL_NUMBER];
    EncodePass      pass;
    int32_t         total_frames;
} EncContext;

static EbErrorType enc_context_ctor(EncApp* enc_app, EncContext* enc_context, int32_t argc, char *argv[], EncodePass pass)
{
    memset(enc_context, 0, sizeof(*enc_context));
    uint32_t num_channels = get_number_of_channels(argc, argv);
    if (num_channels == 0)
        return EB_ErrorBadParameter;

    enc_context->pass = pass;

    EbConfig        **configs = enc_context->configs;
    EbAppContext    **app_callbacks = enc_context->app_callbacks;
    EbErrorType     return_error   = EB_ErrorNone;

    enc_context->num_channels = num_channels;
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        configs[inst_cnt] = eb_config_ctor(pass);
        if (!configs[inst_cnt]) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Initialize appCallback
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        app_callbacks[inst_cnt] = (EbAppContext *)malloc(sizeof(EbAppContext));
        if (!app_callbacks[inst_cnt]) {
            return EB_ErrorInsufficientResources;
        }
    }

    char **warning = enc_context->warning;

    for (int token_id = 0; token_id < MAX_NUM_TOKENS; token_id++) {
        warning[token_id] = (char *)malloc(WARNING_LENGTH);
        strcpy_s(warning[token_id], WARNING_LENGTH, "");
    }
    // Process any command line options, including the configuration file
    // Read all configuration files.
    return_error = read_command_line(argc, argv, configs, num_channels, enc_context->return_errors, warning);
    if (return_error != EB_ErrorNone) {
        fprintf(stderr, "Error in configuration, could not begin encoding! ... \n");
        fprintf(stderr, "Run %s --help for a list of options\n", argv[0]);
        return return_error;
    }
    // Set main thread affinity
    if (configs[0]->target_socket != -1)
        assign_app_thread_group(configs[0]->target_socket);

    EbErrorType*     return_errors = enc_context->return_errors;
    // Init the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        if (return_errors[inst_cnt] == EB_ErrorNone) {
            configs[inst_cnt]->active_channel_count = num_channels;
            configs[inst_cnt]->channel_id           = inst_cnt;

            start_time((uint64_t *)&configs[inst_cnt]->performance_context.lib_start_time[0],
                        (uint64_t *)&configs[inst_cnt]->performance_context.lib_start_time[1]);

            if (pass == ENCODE_FIRST_PASS) {
                //TODO: remove this if we can use a quick first pass in svt av1 library.
                configs[inst_cnt]->enc_mode = MAX_ENC_PRESET;
                configs[inst_cnt]->look_ahead_distance = 1;
                configs[inst_cnt]->enable_tpl_la = 0;
                configs[inst_cnt]->rate_control_mode = 0;
            }
            return_errors[inst_cnt] = set_two_passes_stats(configs[inst_cnt], pass,
                &enc_app->rc_twopasses_stats, num_channels);
            if (return_errors[inst_cnt] == EB_ErrorNone) {
                return_errors[inst_cnt] = init_encoder(
                    configs[inst_cnt], app_callbacks[inst_cnt], inst_cnt);
            }
            return_error = (EbErrorType)(return_error | return_errors[inst_cnt]);
        } else
            enc_context->channel_active[inst_cnt] = EB_FALSE;
    }
    return EB_ErrorNone;
}

static void enc_context_dctor(EncContext* enc_context){
    // DeInit Encoder
    for (int32_t inst_cnt = enc_context->num_channels - 1; inst_cnt >= 0; --inst_cnt) {
        de_init_encoder(enc_context->app_callbacks[inst_cnt],
                        inst_cnt);
        eb_config_dtor(enc_context->configs[inst_cnt]);
        free(enc_context->app_callbacks[inst_cnt]);
    }

    for (uint32_t warning_id = 0; warning_id < MAX_NUM_TOKENS; warning_id++)
        free(enc_context->warning[warning_id]);
}

double get_psnr(double sse, double max);
static void print_summary(const EncContext* const enc_context)
{
    EbConfig*   const* configs = enc_context->configs;
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        if (enc_context->exit_cond[inst_cnt] == APP_ExitConditionFinished &&
            enc_context->return_errors[inst_cnt] == EB_ErrorNone) {
            uint64_t frame_count =
                (uint32_t)configs[inst_cnt]->performance_context.frame_count;
            uint32_t max_luma_value = (configs[inst_cnt]->encoder_bit_depth == 8) ? 255
                                                                                    : 1023;
            double max_luma_sse = (double)max_luma_value * max_luma_value *
                (configs[inst_cnt]->source_width * configs[inst_cnt]->source_height);
            double max_chroma_sse = (double)max_luma_value * max_luma_value *
                (configs[inst_cnt]->source_width / 2 * configs[inst_cnt]->source_height /
                    2);

            if ((configs[inst_cnt]->frame_rate_numerator != 0 &&
                    configs[inst_cnt]->frame_rate_denominator != 0) ||
                configs[inst_cnt]->frame_rate != 0) {
                double frame_rate = configs[inst_cnt]->frame_rate_numerator &&
                        configs[inst_cnt]->frame_rate_denominator
                    ? (double)configs[inst_cnt]->frame_rate_numerator /
                        (double)configs[inst_cnt]->frame_rate_denominator
                    : configs[inst_cnt]->frame_rate > 1000
                        // Correct for 16-bit fixed-point fractional precision
                        ? (double)configs[inst_cnt]->frame_rate / (1 << 16)
                        : (double)configs[inst_cnt]->frame_rate;

                if (configs[inst_cnt]->stat_report) {
                    if (configs[inst_cnt]->stat_file) {
                        fprintf(configs[inst_cnt]->stat_file,
                                "\nSUMMARY "
                                "------------------------------------------------------"
                                "---------------\n");
                        fprintf(configs[inst_cnt]->stat_file,
                                "\n\t\t\t\tAverage PSNR (using per-frame "
                                "PSNR)\t\t|\tOverall PSNR (using per-frame MSE)\t\t|"
                                "\tAverage SSIM\n");
                        fprintf(configs[inst_cnt]->stat_file,
                                "Total Frames\tAverage QP  \tY-PSNR   \tU-PSNR   "
                                "\tV-PSNR\t\t| \tY-PSNR   \tU-PSNR   \tV-PSNR   \t|"
                                "\tY-SSIM   \tU-SSIM   \tV-SSIM   "
                                "\t|\tBitrate\n");
                        fprintf(
                            configs[inst_cnt]->stat_file,
                            "%10ld  \t   %2.2f    \t%3.2f dB\t%3.2f dB\t%3.2f dB  "
                            "\t|\t%3.2f dB\t%3.2f dB\t%3.2f dB \t|\t%1.5f \t%1.5f "
                            "\t%1.5f\t\t|\t%.2f kbps\n",
                            (long int)frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_qp /
                                frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_luma_psnr /
                                frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_cb_psnr /
                                frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_cr_psnr /
                                frame_count,
                            (float)(get_psnr(
                                (configs[inst_cnt]->performance_context.sum_luma_sse /
                                    frame_count),
                                max_luma_sse)),
                            (float)(get_psnr(
                                (configs[inst_cnt]->performance_context.sum_cb_sse /
                                    frame_count),
                                max_chroma_sse)),
                            (float)(get_psnr(
                                (configs[inst_cnt]->performance_context.sum_cr_sse /
                                    frame_count),
                                max_chroma_sse)),
                            (float)configs[inst_cnt]->performance_context.sum_luma_ssim /
                                frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_cb_ssim /
                                frame_count,
                            (float)configs[inst_cnt]->performance_context.sum_cr_ssim /
                                frame_count,
                            ((double)(configs[inst_cnt]->performance_context.byte_count
                                        << 3) *
                                frame_rate / (configs[inst_cnt]->frames_encoded * 1000)));
                    }
                }

                fprintf(stderr,
                        "\nSUMMARY --------------------------------- Channel %u  "
                        "--------------------------------\n",
                        inst_cnt + 1);
                {
                    fprintf(stderr,
                            "Total Frames\t\tFrame Rate\t\tByte Count\t\tBitrate\n");
                    fprintf(
                        stderr,
                        "%12d\t\t%4.2f fps\t\t%10.0f\t\t%5.2f kbps\n",
                        (int32_t)frame_count,
                        (double)frame_rate,
                        (double)configs[inst_cnt]->performance_context.byte_count,
                        ((double)(configs[inst_cnt]->performance_context.byte_count << 3) *
                            frame_rate / (configs[inst_cnt]->frames_encoded * 1000)));
                }

                if (configs[inst_cnt]->stat_report) {
                    fprintf(stderr,
                            "\n\t\tAverage PSNR (using per-frame "
                            "PSNR)\t\t|\tOverall PSNR (using per-frame MSE)\t\t|\t"
                            "Average SSIM\n");
                    fprintf(stderr,
                            "Average "
                            "QP\tY-PSNR\t\tU-PSNR\t\tV-PSNR\t\t|\tY-PSNR\t\tU-"
                            "PSNR\t\tV-PSNR\t\t|\tY-SSIM\tU-SSIM\tV-SSIM\n");
                    fprintf(
                        stderr,
                        "%11.2f\t%4.2f dB\t%4.2f dB\t%4.2f dB\t|\t%4.2f "
                        "dB\t%4.2f dB\t%4.2f dB\t|\t%1.5f\t%1.5f\t%1.5f\n",
                        (float)configs[inst_cnt]->performance_context.sum_qp / frame_count,
                        (float)configs[inst_cnt]->performance_context.sum_luma_psnr /
                            frame_count,
                        (float)configs[inst_cnt]->performance_context.sum_cb_psnr /
                            frame_count,
                        (float)configs[inst_cnt]->performance_context.sum_cr_psnr /
                            frame_count,
                        (float)(get_psnr(
                            (configs[inst_cnt]->performance_context.sum_luma_sse /
                                frame_count),
                            max_luma_sse)),
                        (float)(get_psnr(
                            (configs[inst_cnt]->performance_context.sum_cb_sse /
                                frame_count),
                            max_chroma_sse)),
                        (float)(get_psnr(
                            (configs[inst_cnt]->performance_context.sum_cr_sse /
                                frame_count),
                            max_chroma_sse)),
                        (float)configs[inst_cnt]->performance_context.sum_luma_ssim /
                            frame_count,
                        (float)configs[inst_cnt]->performance_context.sum_cb_ssim /
                            frame_count,
                        (float)configs[inst_cnt]->performance_context.sum_cr_ssim /
                            frame_count);
                }

                fflush(stdout);
            }
        }
    }
    fprintf(stderr, "\n");
    fflush(stdout);
}

static void print_performance(const EncContext* const enc_context)
{
    EbConfig*   const* configs = enc_context->configs;
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        if (enc_context->exit_cond[inst_cnt] == APP_ExitConditionFinished &&
            enc_context->return_errors[inst_cnt] == EB_ErrorNone) {
            if (configs[inst_cnt]->stop_encoder == EB_FALSE) {
                fprintf(stderr,
                        "\nChannel %u\nAverage Speed:\t\t%.3f fps\nTotal Encoding Time:\t%.0f "
                        "ms\nTotal Execution Time:\t%.0f ms\nAverage Latency:\t%.0f ms\nMax "
                        "Latency:\t\t%u ms\n",
                        (uint32_t)(inst_cnt + 1),
                        configs[inst_cnt]->performance_context.average_speed,
                        configs[inst_cnt]->performance_context.total_encode_time * 1000,
                        configs[inst_cnt]->performance_context.total_execution_time * 1000,
                        configs[inst_cnt]->performance_context.average_latency,
                        (uint32_t)(configs[inst_cnt]->performance_context.max_latency));
            } else
                fprintf(
                    stderr, "\nChannel %u Encoding Interrupted\n", (uint32_t)(inst_cnt + 1));
        } else if (enc_context->return_errors[inst_cnt] == EB_ErrorInsufficientResources)
            fprintf(stderr, "Could not allocate enough memory for channel %u\n", inst_cnt + 1);
        else
            fprintf(stderr,
                    "Error encoding at channel %u! Check error log file for more details "
                    "... \n",
                    inst_cnt + 1);
    }
}

static void print_warnnings(const EncContext* const enc_context)
{
    char* const* warning = enc_context->warning;
    for (uint32_t warning_id = 0;; warning_id++) {
        if (*warning[warning_id] == '-')
            fprintf(stderr, "warning: %s\n", warning[warning_id]);
        else if (*warning[warning_id + 1] != '-')
                break;
    }

}

EbBool has_active_channel(const EncContext* const enc_context)
{
    // check if all channels are inactive
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        if (enc_context->channel_active[inst_cnt])
            return EB_TRUE;
    }
    return EB_FALSE;
}

void enc_channel_step(EncApp* enc_app, EncContext* enc_context, uint32_t inst_cnt )
{
    EbAppContext **app_callbacks = enc_context->app_callbacks; // Instances App callback date

    AppExitConditionType *exit_cond_output = enc_context->exit_cond_output;// Processing loop exit condition
    AppExitConditionType *exit_cond_recon = enc_context->exit_cond_recon; // Processing loop exit condition
    AppExitConditionType *exit_cond_input = enc_context->exit_cond_input; // Processing loop exit condition
    AppExitConditionType* exit_cond = enc_context->exit_cond;
    EbBool* channel_active = enc_context->channel_active;
    EbConfig **configs = enc_context->configs; // Encoder Configuration
    if (exit_cond_input[inst_cnt] == APP_ExitConditionNone)
        exit_cond_input[inst_cnt] = process_input_buffer(
            configs[inst_cnt], app_callbacks[inst_cnt]);
    if (exit_cond_recon[inst_cnt] == APP_ExitConditionNone)
        exit_cond_recon[inst_cnt] = process_output_recon_buffer(
            configs[inst_cnt], app_callbacks[inst_cnt]);
    if (exit_cond_output[inst_cnt] == APP_ExitConditionNone)
        exit_cond_output[inst_cnt] = process_output_stream_buffer(
            enc_app,
            configs[inst_cnt],
            app_callbacks[inst_cnt],
            (exit_cond_input[inst_cnt] == APP_ExitConditionNone) ||
                    (exit_cond_recon[inst_cnt] == APP_ExitConditionNone)
                ? 0
                : 1, &enc_context->total_frames);

    if (((exit_cond_recon[inst_cnt] == APP_ExitConditionFinished ||
            !configs[inst_cnt]->recon_file) &&
            exit_cond_output[inst_cnt] == APP_ExitConditionFinished &&
            exit_cond_input[inst_cnt] == APP_ExitConditionFinished) ||
        ((exit_cond_recon[inst_cnt] == APP_ExitConditionError &&
            configs[inst_cnt]->recon_file) ||
            exit_cond_output[inst_cnt] == APP_ExitConditionError ||
            exit_cond_input[inst_cnt] == APP_ExitConditionError)) {
        channel_active[inst_cnt] = EB_FALSE;
        if (configs[inst_cnt]->recon_file)
            exit_cond[inst_cnt] = (AppExitConditionType)(
                exit_cond_recon[inst_cnt] | exit_cond_output[inst_cnt] |
                exit_cond_input[inst_cnt]);
        else
            exit_cond[inst_cnt] = (AppExitConditionType)(
                exit_cond_output[inst_cnt] | exit_cond_input[inst_cnt]);
    }

}

static const char *get_pass_name(EncodePass pass) {
    switch (pass) {
    case ENCODE_FIRST_PASS: return "Pass 1/2 ";
    case ENCODE_LAST_PASS: return "Pass 2/2 ";
    default: return "";
    }
}

static EbErrorType encode(EncApp* enc_app, EncContext* enc_context) {
    EbErrorType          return_error   = EB_ErrorNone;

    AppExitConditionType *exit_cond_output = enc_context->exit_cond_output;// Processing loop exit condition
    AppExitConditionType *exit_cond_recon = enc_context->exit_cond_recon; // Processing loop exit condition
    AppExitConditionType *exit_cond_input = enc_context->exit_cond_input; // Processing loop exit condition
    AppExitConditionType* exit_cond = enc_context->exit_cond;
    EbBool* channel_active = enc_context->channel_active;
    EbConfig **configs = enc_context->configs; // Encoder Configuration
    // Get num_channels
    uint32_t num_channels = enc_context->num_channels;

    for (uint32_t inst_cnt = 0; inst_cnt < MAX_CHANNEL_NUMBER; ++inst_cnt) {
        exit_cond[inst_cnt]        = APP_ExitConditionError; // Processing loop exit condition
        exit_cond_output[inst_cnt] = APP_ExitConditionError; // Processing loop exit condition
        exit_cond_recon[inst_cnt]  = APP_ExitConditionError; // Processing loop exit condition
        exit_cond_input[inst_cnt]  = APP_ExitConditionError; // Processing loop exit condition
        channel_active[inst_cnt]   = EB_FALSE;
    }

    EbErrorType* return_errors = enc_context->return_errors;
    EncodePass pass = enc_context->pass;
    // Start the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        if (return_errors[inst_cnt] == EB_ErrorNone) {
            return_error        = (EbErrorType)(return_error & return_errors[inst_cnt]);
            exit_cond[inst_cnt] = APP_ExitConditionNone;
            exit_cond_output[inst_cnt] = APP_ExitConditionNone;
            exit_cond_recon[inst_cnt]  = configs[inst_cnt]->recon_file
                ? APP_ExitConditionNone
                : APP_ExitConditionError;
            exit_cond_input[inst_cnt] = APP_ExitConditionNone;
            channel_active[inst_cnt]  = EB_TRUE;
            start_time(
                (uint64_t *)&configs[inst_cnt]->performance_context.encode_start_time[0],
                (uint64_t *)&configs[inst_cnt]->performance_context.encode_start_time[1]);
        } else {
            exit_cond[inst_cnt]        = APP_ExitConditionError;
            exit_cond_output[inst_cnt] = APP_ExitConditionError;
            exit_cond_recon[inst_cnt]  = APP_ExitConditionError;
            exit_cond_input[inst_cnt]  = APP_ExitConditionError;
        }

#if DISPLAY_MEMORY
        EB_APP_MEMORY();
#endif
    }
    print_warnnings(enc_context);

    fprintf(stderr, "%sEncoding          ", get_pass_name(pass));

    while (has_active_channel(enc_context)) {
        for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
            if (channel_active[inst_cnt] == EB_TRUE) {
                enc_channel_step(enc_app, enc_context, inst_cnt);
            }
        }
    }
    print_summary(enc_context);
    print_performance(enc_context);
    return return_error;
}

EbErrorType enc_app_ctor(EncApp* enc_app)
{
    memset(enc_app, 0, sizeof(*enc_app));
    return EB_ErrorNone;
}

void enc_app_dctor(EncApp* enc_app)
{
    free(enc_app->rc_twopasses_stats.buf);
}

/***************************************
 * Encoder App Main
 ***************************************/
int32_t main(int32_t argc, char *argv[]) {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // GLOBAL VARIABLES
    EbErrorType return_error = EB_ErrorNone; // Error Handling
    uint32_t    passes;
    EncodePass  pass[MAX_ENCODE_PASS];
    EncApp      enc_app;
    EncContext   enc_context;

    signal(SIGINT, event_handler);
    if (get_help(argc, argv))
        return 0;

    enc_app_ctor(&enc_app);
    passes = get_passes(argc, argv, pass);
    for (uint32_t i = 0; i < passes; i++) {
        return_error = enc_context_ctor(&enc_app, &enc_context, argc, argv, pass[i]);

        if (return_error == EB_ErrorNone)
            return_error = encode(&enc_app, &enc_context);

        enc_context_dctor(&enc_context);
        if (return_error != EB_ErrorNone)
            break;
    }
    enc_app_dctor(&enc_app);
    return return_error != EB_ErrorNone;
}
