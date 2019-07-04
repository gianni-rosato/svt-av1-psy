/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "EbSvtAv1Dec.h"
#include "EbDecParamParser.h"
#include "EbMD5Utility.h"

#ifdef _MSC_VER
#include <io.h>     /* _setmode() */
#include <fcntl.h>  /* _O_BINARY */
#endif

int init_pic_buffer(EbSvtIOFormat *pic_buffer, CLInput *cli) {
    switch (cli->fmt) {
    case EB_YUV420:
        pic_buffer->cb_stride = cli->width / 2;
        pic_buffer->cr_stride = cli->width / 2;
        pic_buffer->y_stride = cli->width;
        break;
    default:
        printf("Unsupported colour format. \n");
        return 0;
    }
    pic_buffer->width = cli->width;
    pic_buffer->height = cli->width;
    pic_buffer->luma_ext = NULL;
    pic_buffer->cb_ext = NULL;
    pic_buffer->cr_ext = NULL;
    pic_buffer->origin_x = 0;
    pic_buffer->origin_y = 0;
    return 0;
}

int read_input_frame(CLInput *cli, uint8_t **buffer, size_t *bytes_read,
                     size_t *buffer_size, int64_t *pts)
{
    switch (cli->inFileType)
    {
    case FILE_TYPE_IVF:
        return read_ivf_frame(cli->inFile, buffer, bytes_read, buffer_size, pts);
        break;
    default:
        printf("Unsupported bitstream type. \n");
        return 0;
    }
}

void write_frame(EbBufferHeaderType *recon_buffer, CLInput *cli) {
    EbSvtIOFormat* img = (EbSvtIOFormat*)recon_buffer->p_buffer;

    // Support only for 420 images
    assert(cli->fmt == EB_YUV420);

    const int bytes_per_sample = (cli->bit_depth == EB_EIGHT_BIT) ? 1 : 2;

    // Write luma plane
    unsigned char *buf = img->luma;
    int stride = img->y_stride;
    int w = cli->width;
    int h = cli->height;
    int y = 0;
    for (y = 0; y < h; ++y) {
        fwrite(buf, bytes_per_sample, w, cli->outFile);
        buf += (stride* bytes_per_sample);
    }

    //Write chroma planes
    buf = img->cb;
    stride = img->cb_stride;
    w /= 2;
    h /= 2;
    for (y = 0; y < h; ++y) {
        fwrite(buf, bytes_per_sample, w, cli->outFile);
        buf += (stride* bytes_per_sample);
    }

    buf = img->cr;
    stride = img->cr_stride;
    for (y = 0; y < h; ++y) {
        fwrite(buf, bytes_per_sample, w, cli->outFile);
        buf += (stride* bytes_per_sample);
    }

    fflush(cli->outFile);
}

/***************************************
 * Decoder App Main
 ***************************************/
int32_t main(int32_t argc, char* argv[])
{
#ifdef _MSC_VER
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // GLOBAL VARIABLES
    EbErrorType            return_error = EB_ErrorNone; // Error Handling
    EbSvtAv1DecConfiguration  *config_ptr = (EbSvtAv1DecConfiguration*)
                                            malloc(sizeof(EbSvtAv1DecConfiguration));
    CLInput cli;
    cli.inFile = NULL;
    cli.outFile = NULL;
    cli.enable_md5 = 0;

    uint64_t stop_after = 0;
    uint32_t in_frame = 0;

    MD5Context md5_ctx;
    unsigned char md5_digest[16];

    uint8_t *buf = NULL;
    size_t bytes_in_buffer = 0, buffer_size = 0;

    // Print Decoder Info
    printf("\n**WARNING** decoder is not feature complete\n");
    printf("Current support: intra & inter(no OBMC, no Compund, no Wedge, ");
    printf("no Temporal Scan) tools (no SCC, no loop filters)\n\n");

    printf("-------------------------------------\n");
    printf("SVT-AV1 Decoder Sample Application v1.2.0\n");
    printf("Platform:   %u bit\n", (unsigned) sizeof(void*) * 8);
#if ( defined( _MSC_VER ) && (_MSC_VER < 1910) )
    printf("Compiler: VS13\n");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1910) )
    printf("Compiler: VS17\n");
#elif defined(__INTEL_COMPILER)
    printf("Compiler: Intel\n");
#elif defined(__GNUC__)
    printf("Compiler: GCC\n");
#else
    printf("Compiler: unknown\n");
#endif

    printf("APP Build date: %s %s\n", __DATE__, __TIME__);

    // Initialize config
    if (!config_ptr)
        return EB_ErrorInsufficientResources;
    EbComponentType* p_handle;
    void* p_app_data = NULL;

    return_error |= eb_dec_init_handle(&p_handle, p_app_data, config_ptr);
    if (return_error != EB_ErrorNone) goto fail;

    if (read_command_line(argc, argv, config_ptr, &cli) == 0 &&
        !eb_svt_dec_set_parameter(p_handle, config_ptr)) {
        return_error = eb_init_decoder(p_handle);
        if (return_error != EB_ErrorNone) {
            return_error |= eb_dec_deinit_handle(p_handle);
            goto fail;
        }

        assert(config_ptr->max_color_format == EB_YUV420);
        assert(config_ptr->max_bit_depth < EB_TWELVE_BIT);

        int enable_md5 = cli.enable_md5;

        EbBufferHeaderType *recon_buffer = NULL;
        recon_buffer = (EbBufferHeaderType*)malloc(sizeof(EbBufferHeaderType));
        recon_buffer->p_buffer = (uint8_t *)malloc(sizeof(EbSvtIOFormat));

        int size = (config_ptr->max_bit_depth == EB_EIGHT_BIT) ?
                                sizeof(uint8_t) : sizeof(uint16_t);
        size = size * cli.height * cli.width;

        ((EbSvtIOFormat *)recon_buffer->p_buffer)->luma = (uint8_t*)malloc(size);
        ((EbSvtIOFormat *)recon_buffer->p_buffer)->cb = (uint8_t*)malloc(size >> 2);
        ((EbSvtIOFormat *)recon_buffer->p_buffer)->cr = (uint8_t*)malloc(size >> 2);
        if (!init_pic_buffer((EbSvtIOFormat*)recon_buffer->p_buffer, &cli)) {
            printf("Decoding \n");
            EbAV1StreamInfo *stream_info = (EbAV1StreamInfo*)malloc(sizeof(EbAV1StreamInfo));
            EbAV1FrameInfo *frame_info = (EbAV1FrameInfo*)malloc(sizeof(EbAV1FrameInfo));

            if (config_ptr->skip_frames)
                fprintf(stderr, "Skipping first %" PRIu64 " frames.\n", config_ptr->skip_frames);
            uint64_t skip_frame = config_ptr->skip_frames;
            while (skip_frame) {
                if (!read_input_frame(&cli, &buf, &bytes_in_buffer, &buffer_size, NULL)) break;
                skip_frame--;
            }
            stop_after = config_ptr->frames_to_be_decoded;
            if (enable_md5)
                md5_init(&md5_ctx);
            // Input Loop Thread
            while (read_input_frame(&cli, &buf, &bytes_in_buffer, &buffer_size, NULL)) {
                if (!stop_after || in_frame < stop_after) {
                    return_error |= eb_svt_decode_frame(p_handle, buf, bytes_in_buffer);

                    in_frame++;

                    if (eb_svt_dec_get_picture(p_handle, recon_buffer, stream_info, frame_info) != EB_DecNoOutputPicture) {
                        if (enable_md5)
                            write_md5(recon_buffer, &cli, &md5_ctx);
                        else
                            write_frame(recon_buffer, &cli);
                    }
                }
                else break;
            }

            if (enable_md5) {
                md5_final(md5_digest, &md5_ctx);
                print_md5(md5_digest);
            }

            return_error |= eb_deinit_decoder(p_handle);

            free(frame_info);
            free(stream_info);
        }
        free(((EbSvtIOFormat *)recon_buffer->p_buffer)->cr);
        free(((EbSvtIOFormat *)recon_buffer->p_buffer)->cb);
        free(((EbSvtIOFormat *)recon_buffer->p_buffer)->luma);
        free(recon_buffer->p_buffer);
        free(recon_buffer);
        free(buf);
    }
    else
        printf("Error in configuration. \n");
    return_error |= eb_dec_deinit_handle(p_handle);

fail:
    if (cli.inFile) fclose(cli.inFile);
    if (cli.outFile) fclose(cli.outFile);

    return return_error;
}
