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
#include <errno.h>
#endif

#if !defined(_WIN32) || !defined(HAVE_STRNLEN_S)
#include "safe_str_lib.h"
#endif

/***************************************
 * External Functions
 ***************************************/
void process_input_buffer(EncChannel* c);

void process_output_recon_buffer(EncChannel* c);

void process_output_stream_buffer(EncChannel* c, EncApp* enc_app, int32_t *frame_count);

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
    EncChannel      channels[MAX_CHANNEL_NUMBER];
    char            *warning[MAX_NUM_TOKENS];

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
    EbErrorType     return_error   = EB_ErrorNone;

    enc_context->num_channels = num_channels;
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        return_error = enc_channel_ctor(enc_context->channels + inst_cnt, pass);
        if (return_error != EB_ErrorNone)
            return return_error;
    }

    char **warning = enc_context->warning;

    for (int token_id = 0; token_id < MAX_NUM_TOKENS; token_id++) {
        warning[token_id] = (char *)malloc(WARNING_LENGTH);
        if (!warning[token_id]) {
            return EB_ErrorInsufficientResources;
        }
        strcpy_s(warning[token_id], WARNING_LENGTH, "");
    }
    // Process any command line options, including the configuration file
    // Read all configuration files.
    return_error = read_command_line(argc, argv, enc_context->channels, num_channels, warning);
    if (return_error != EB_ErrorNone) {
        fprintf(stderr, "Error in configuration, could not begin encoding! ... \n");
        fprintf(stderr, "Run %s --help for a list of options\n", argv[0]);
        return return_error;
    }
    // Set main thread affinity
    if (enc_context->channels[0].config->target_socket != -1)
        assign_app_thread_group(enc_context->channels[0].config->target_socket);

    // Init the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        EncChannel* c = enc_context->channels + inst_cnt;
        if (c->return_error == EB_ErrorNone) {
            EbConfig* config = c->config;
            config->active_channel_count = num_channels;
            config->channel_id           = inst_cnt;

            app_svt_av1_get_time(&config->performance_context.lib_start_time[0],
                                 &config->performance_context.lib_start_time[1]);

            c->return_error = set_two_passes_stats(config, pass,
                &enc_app->rc_twopasses_stats, num_channels);
            if (c->return_error == EB_ErrorNone) {
                c->return_error = init_encoder(
                    config, c->app_callback, inst_cnt);
            }
            return_error = (EbErrorType)(return_error | c->return_error);
        } else
            c->active = EB_FALSE;
    }
    return EB_ErrorNone;
}

static void enc_context_dctor(EncContext* enc_context){
    // DeInit Encoder
    for (int32_t inst_cnt = enc_context->num_channels - 1; inst_cnt >= 0; --inst_cnt) {
        EncChannel* c = enc_context->channels + inst_cnt;
        de_init_encoder(c->app_callback,
                        inst_cnt);
        enc_channel_dctor(c);
    }

    for (uint32_t warning_id = 0; warning_id < MAX_NUM_TOKENS; warning_id++)
        free(enc_context->warning[warning_id]);
}

double get_psnr(double sse, double max);
static void print_summary(const EncContext* const enc_context)
{
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        const EncChannel* const c = &enc_context->channels[inst_cnt];
        const EbConfig*  config = c->config;
        if (c->exit_cond == APP_ExitConditionFinished &&
            c->return_error == EB_ErrorNone) {
            uint64_t frame_count =
                (uint32_t)config->performance_context.frame_count;
            uint32_t max_luma_value = (config->encoder_bit_depth == 8) ? 255
                                                                                    : 1023;
            double max_luma_sse = (double)max_luma_value * max_luma_value *
                (config->source_width * config->source_height);
            double max_chroma_sse = (double)max_luma_value * max_luma_value *
                (config->source_width / 2 * config->source_height /
                    2);

            if ((config->frame_rate_numerator != 0 &&
                    config->frame_rate_denominator != 0) ||
                config->frame_rate != 0) {
                double frame_rate = config->frame_rate_numerator &&
                        config->frame_rate_denominator
                    ? (double)config->frame_rate_numerator /
                        (double)config->frame_rate_denominator
                    : config->frame_rate > 1000
                        // Correct for 16-bit fixed-point fractional precision
                        ? (double)config->frame_rate / (1 << 16)
                        : (double)config->frame_rate;

                if (config->stat_report) {
                    if (config->stat_file) {
                        fprintf(config->stat_file,
                                "\nSUMMARY "
                                "------------------------------------------------------"
                                "---------------\n");
                        fprintf(config->stat_file,
                                "\n\t\t\t\tAverage PSNR (using per-frame "
                                "PSNR)\t\t|\tOverall PSNR (using per-frame MSE)\t\t|"
                                "\tAverage SSIM\n");
                        fprintf(config->stat_file,
                                "Total Frames\tAverage QP  \tY-PSNR   \tU-PSNR   "
                                "\tV-PSNR\t\t| \tY-PSNR   \tU-PSNR   \tV-PSNR   \t|"
                                "\tY-SSIM   \tU-SSIM   \tV-SSIM   "
                                "\t|\tBitrate\n");
                        fprintf(
                            config->stat_file,
                            "%10ld  \t   %2.2f    \t%3.2f dB\t%3.2f dB\t%3.2f dB  "
                            "\t|\t%3.2f dB\t%3.2f dB\t%3.2f dB \t|\t%1.5f \t%1.5f "
                            "\t%1.5f\t\t|\t%.2f kbps\n",
                            (long int)frame_count,
                            (float)config->performance_context.sum_qp /
                                frame_count,
                            (float)config->performance_context.sum_luma_psnr /
                                frame_count,
                            (float)config->performance_context.sum_cb_psnr /
                                frame_count,
                            (float)config->performance_context.sum_cr_psnr /
                                frame_count,
                            (float)(get_psnr(
                                (config->performance_context.sum_luma_sse /
                                    frame_count),
                                max_luma_sse)),
                            (float)(get_psnr(
                                (config->performance_context.sum_cb_sse /
                                    frame_count),
                                max_chroma_sse)),
                            (float)(get_psnr(
                                (config->performance_context.sum_cr_sse /
                                    frame_count),
                                max_chroma_sse)),
                            (float)config->performance_context.sum_luma_ssim /
                                frame_count,
                            (float)config->performance_context.sum_cb_ssim /
                                frame_count,
                            (float)config->performance_context.sum_cr_ssim /
                                frame_count,
                            ((double)(config->performance_context.byte_count
                                        << 3) *
                                frame_rate / (config->frames_encoded * 1000)));
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
                        (double)config->performance_context.byte_count,
                        ((double)(config->performance_context.byte_count << 3) *
                            frame_rate / (config->frames_encoded * 1000)));
                }

                if (config->stat_report) {
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
                        (float)config->performance_context.sum_qp / frame_count,
                        (float)config->performance_context.sum_luma_psnr /
                            frame_count,
                        (float)config->performance_context.sum_cb_psnr /
                            frame_count,
                        (float)config->performance_context.sum_cr_psnr /
                            frame_count,
                        (float)(get_psnr(
                            (config->performance_context.sum_luma_sse /
                                frame_count),
                            max_luma_sse)),
                        (float)(get_psnr(
                            (config->performance_context.sum_cb_sse /
                                frame_count),
                            max_chroma_sse)),
                        (float)(get_psnr(
                            (config->performance_context.sum_cr_sse /
                                frame_count),
                            max_chroma_sse)),
                        (float)config->performance_context.sum_luma_ssim /
                            frame_count,
                        (float)config->performance_context.sum_cb_ssim /
                            frame_count,
                        (float)config->performance_context.sum_cr_ssim /
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
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        const EncChannel* c = enc_context->channels + inst_cnt;
        if (c->exit_cond == APP_ExitConditionFinished &&
            c->return_error == EB_ErrorNone) {
            EbConfig* config = c->config;
            if (config->stop_encoder == EB_FALSE) {
                fprintf(stderr,
                        "\nChannel %u\nAverage Speed:\t\t%.3f fps\nTotal Encoding Time:\t%.0f "
                        "ms\nTotal Execution Time:\t%.0f ms\nAverage Latency:\t%.0f ms\nMax "
                        "Latency:\t\t%u ms\n",
                        (uint32_t)(inst_cnt + 1),
                        config->performance_context.average_speed,
                        config->performance_context.total_encode_time * 1000,
                        config->performance_context.total_execution_time * 1000,
                        config->performance_context.average_latency,
                        (uint32_t)(config->performance_context.max_latency));
            } else
                fprintf(
                    stderr, "\nChannel %u Encoding Interrupted\n", (uint32_t)(inst_cnt + 1));
        } else if (c->return_error == EB_ErrorInsufficientResources)
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

static EbBool is_active(const EncChannel* c)
{
    return c->active;
}

static EbBool has_active_channel(const EncContext* const enc_context)
{
    // check if all channels are inactive
    for (uint32_t inst_cnt = 0; inst_cnt < enc_context->num_channels; ++inst_cnt) {
        if (is_active(enc_context->channels + inst_cnt))
            return EB_TRUE;
    }
    return EB_FALSE;
}

static void enc_channel_step(EncChannel* c, EncApp* enc_app, EncContext* enc_context)
{
    EbConfig* config = c->config;
    process_input_buffer(c);
    process_output_recon_buffer(c);
    process_output_stream_buffer(c,
            enc_app,
            &enc_context->total_frames);

    if (((c->exit_cond_recon == APP_ExitConditionFinished ||
            !config->recon_file) &&
            c->exit_cond_output == APP_ExitConditionFinished &&
            c->exit_cond_input == APP_ExitConditionFinished) ||
        ((c->exit_cond_recon == APP_ExitConditionError &&
            config->recon_file) ||
            c->exit_cond_output == APP_ExitConditionError ||
            c->exit_cond_input == APP_ExitConditionError)) {
        c->active = EB_FALSE;
        if (config->recon_file)
            c->exit_cond = (AppExitConditionType)(
                c->exit_cond_recon | c->exit_cond_output |
                c->exit_cond_input);
        else
            c->exit_cond = (AppExitConditionType)(
                c->exit_cond_output | c->exit_cond_input);
    }

}

static const char *get_pass_name(EncodePass pass) {
    switch (pass) {
    case ENCODE_FIRST_PASS: return "Pass 1/2 ";
    case ENCODE_LAST_PASS: return "Pass 2/2 ";
    default: return "";
    }
}

static void enc_channel_start(EncChannel* c)
{
    if (c->return_error == EB_ErrorNone) {
        EbConfig* config = c->config;
        c->exit_cond = APP_ExitConditionNone;
        c->exit_cond_output = APP_ExitConditionNone;
        c->exit_cond_recon  = config->recon_file
            ? APP_ExitConditionNone
            : APP_ExitConditionError;
        c->exit_cond_input = APP_ExitConditionNone;
        c->active  = EB_TRUE;
        app_svt_av1_get_time(&config->performance_context.encode_start_time[0],
                             &config->performance_context.encode_start_time[1]);
    }
}

static EbErrorType encode(EncApp* enc_app, EncContext* enc_context) {
    EbErrorType          return_error   = EB_ErrorNone;

    // Get num_channels
    uint32_t num_channels = enc_context->num_channels;

    EncodePass pass = enc_context->pass;
    // Start the Encoder
    for (uint32_t inst_cnt = 0; inst_cnt < num_channels; ++inst_cnt) {
        enc_channel_start(enc_context->channels + inst_cnt);

#if DISPLAY_MEMORY
        EB_APP_MEMORY();
#endif
    }
    print_warnnings(enc_context);

    fprintf(stderr, "%sEncoding          ", get_pass_name(pass));

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
