/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2EFramework.cc
 *
 * @brief Impelmentation of End to End test framework
 *
 * @author Cidana-Edmond Cidana-Ryan
 *
 ******************************************************************************/

#include "EbSvtAv1Enc.h"
#include "Y4mVideoSource.h"
#include "YuvVideoSource.h"
#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "RefDecoder.h"
#include "SvtAv1E2EFramework.h"
#include "CompareTools.h"

#if _WIN32
#define fseeko64 _fseeki64
#define ftello64 _ftelli64
#else
#define fseeko64 fseek
#define ftello64 ftell
#endif

using namespace svt_av1_e2e_test;
using namespace svt_av1_e2e_tools;

VideoSource *SvtAv1E2ETestFramework::prepare_video_src(
    const TestVideoVector &vector) {
    VideoSource *video_src = nullptr;
    switch (std::get<1>(vector)) {
    case YUV_VIDEO_FILE:
        video_src = new YuvVideoSource(std::get<0>(vector),
                                       std::get<2>(vector),
                                       std::get<3>(vector),
                                       std::get<4>(vector),
                                       std::get<5>(vector));
        break;
    case Y4M_VIDEO_FILE:
        video_src = new Y4MVideoSource(std::get<0>(vector),
                                       std::get<2>(vector),
                                       std::get<3>(vector),
                                       std::get<4>(vector),
                                       std::get<5>(vector),
                                       std::get<6>(vector));
        break;
    default: assert(0); break;
    }
    return video_src;
}

void SvtAv1E2ETestFramework::trans_src_param(const VideoSource *source,
                                             EbSvtAv1EncConfiguration &config) {
    VideoColorFormat fmt = source->get_image_format();
    switch (fmt) {
    case IMG_FMT_420:
    case IMG_FMT_420P10_PACKED: config.encoder_color_format = EB_YUV420; break;
    case IMG_FMT_422:
    case IMG_FMT_422P10_PACKED: config.encoder_color_format = EB_YUV422; break;
    case IMG_FMT_444:
    case IMG_FMT_444P10_PACKED: config.encoder_color_format = EB_YUV444; break;
    default: assert(0); break;
    }
    config.source_width = source->get_width_with_padding();
    config.source_height = source->get_height_with_padding();
    config.encoder_bit_depth = source->get_bit_depth();
    config.compressed_ten_bit_format = source->get_compressed_10bit_mode();
    config.frames_to_be_encoded = source->get_frame_count();
}

SvtAv1E2ETestFramework::SvtAv1E2ETestFramework()
    : video_src_(SvtAv1E2ETestFramework::prepare_video_src(GetParam())),
      psnr_src_(SvtAv1E2ETestFramework::prepare_video_src(GetParam())) {
    memset(&av1enc_ctx_, 0, sizeof(av1enc_ctx_));
    recon_queue_ = nullptr;
    refer_dec_ = nullptr;
    output_file_ = nullptr;
    obu_frame_header_size_ = 0;
    collect_ = nullptr;
    ref_compare_ = nullptr;
    start_pos_ = std::get<7>(GetParam());
    frames_to_test_ = std::get<8>(GetParam());
    printf("start: %d, count: %d\n", start_pos_, frames_to_test_);
}

SvtAv1E2ETestFramework::~SvtAv1E2ETestFramework() {
    if (video_src_) {
        delete video_src_;
        video_src_ = nullptr;
    }
    if (recon_queue_) {
        delete recon_queue_;
        recon_queue_ = nullptr;
    }
    if (refer_dec_) {
        delete refer_dec_;
        refer_dec_ = nullptr;
    }
    if (output_file_) {
        delete output_file_;
        output_file_ = nullptr;
    }
    if (collect_) {
        delete collect_;
        collect_ = nullptr;
    }
    if (psnr_src_) {
        psnr_src_->close_source();
        delete psnr_src_;
        psnr_src_ = nullptr;
    }
    if (ref_compare_) {
        delete ref_compare_;
        ref_compare_ = nullptr;
    }
}

void SvtAv1E2ETestFramework::SetUp() {
    EbErrorType return_error = EB_ErrorNone;

    // check for video source
    ASSERT_NE(video_src_, nullptr) << "video source create failed!";
    return_error = video_src_->open_source(start_pos_, frames_to_test_);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "open_source return error:" << return_error;
    // Check input parameters
    uint32_t width = video_src_->get_width_with_padding();
    uint32_t height = video_src_->get_height_with_padding();
    uint32_t bit_depth = video_src_->get_bit_depth();
    ASSERT_GT(width, 0) << "Video vector width error.";
    ASSERT_GT(height, 0) << "Video vector height error.";
    ASSERT_TRUE(bit_depth == 10 || bit_depth == 8)
        << "Video vector bitDepth error.";

    //
    // Init handle
    //
    return_error = eb_init_handle(
        &av1enc_ctx_.enc_handle, &av1enc_ctx_, &av1enc_ctx_.enc_params);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_init_handle return error:" << return_error;
    ASSERT_NE(av1enc_ctx_.enc_handle, nullptr)
        << "eb_init_handle return null handle.";
    trans_src_param(video_src_, av1enc_ctx_.enc_params);
    av1enc_ctx_.enc_params.recon_enabled = 0;

    //
    // Prepare input and output buffer
    //
    // Input Buffer
    av1enc_ctx_.input_picture_buffer = new EbBufferHeaderType;
    ASSERT_NE(av1enc_ctx_.input_picture_buffer, nullptr)
        << "Malloc memory for inputPictureBuffer failed.";
    av1enc_ctx_.input_picture_buffer->p_buffer = nullptr;
    av1enc_ctx_.input_picture_buffer->size = sizeof(EbBufferHeaderType);
    av1enc_ctx_.input_picture_buffer->p_app_private = nullptr;
    av1enc_ctx_.input_picture_buffer->pic_type = EB_AV1_INVALID_PICTURE;

    // Output buffer
    av1enc_ctx_.output_stream_buffer = new EbBufferHeaderType;
    ASSERT_NE(av1enc_ctx_.output_stream_buffer, nullptr)
        << "Malloc memory for outputStreamBuffer failed.";
    av1enc_ctx_.output_stream_buffer->p_buffer =
        new uint8_t[EB_OUTPUTSTREAMBUFFERSIZE_MACRO(width * height)];
    ASSERT_NE(av1enc_ctx_.output_stream_buffer->p_buffer, nullptr)
        << "Malloc memory for outputStreamBuffer->p_buffer failed.";
    av1enc_ctx_.output_stream_buffer->size = sizeof(EbBufferHeaderType);
    av1enc_ctx_.output_stream_buffer->n_alloc_len =
        EB_OUTPUTSTREAMBUFFERSIZE_MACRO(width * height);
    av1enc_ctx_.output_stream_buffer->p_app_private = nullptr;
    av1enc_ctx_.output_stream_buffer->pic_type = EB_AV1_INVALID_PICTURE;

    // initialize test
    init_test();
}

void SvtAv1E2ETestFramework::TearDown() {
    // close test before teardown
    close_test();

    EbErrorType return_error = EB_ErrorNone;
    // Destruct the component
    return_error = eb_deinit_handle(av1enc_ctx_.enc_handle);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_deinit_handle return error:" << return_error;
    av1enc_ctx_.enc_handle = nullptr;

    // Clear
    if (av1enc_ctx_.output_stream_buffer != nullptr) {
        if (av1enc_ctx_.output_stream_buffer->p_buffer != nullptr) {
            delete[] av1enc_ctx_.output_stream_buffer->p_buffer;
        }
        delete av1enc_ctx_.output_stream_buffer;
        av1enc_ctx_.output_stream_buffer = nullptr;
    }
    if (av1enc_ctx_.input_picture_buffer != nullptr) {
        delete av1enc_ctx_.input_picture_buffer;
        av1enc_ctx_.input_picture_buffer = nullptr;
    }

    ASSERT_NE(video_src_, nullptr);
    video_src_->close_source();
}

/** initialization for test */
void SvtAv1E2ETestFramework::init_test() {
    EbErrorType return_error = eb_svt_enc_set_parameter(
        av1enc_ctx_.enc_handle, &av1enc_ctx_.enc_params);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_svt_enc_set_parameter return error:" << return_error;

    return_error = eb_init_encoder(av1enc_ctx_.enc_handle);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_init_encoder return error:" << return_error;

    // Get ivf header
    return_error = eb_svt_enc_stream_header(av1enc_ctx_.enc_handle,
                                            &av1enc_ctx_.output_stream_buffer);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_svt_enc_stream_header return error:" << return_error;
    ASSERT_NE(av1enc_ctx_.output_stream_buffer, nullptr)
        << "eb_svt_enc_stream_header return null output buffer."
        << return_error;

#if TILES
    EbBool has_tiles = (EbBool)(av1enc_ctx_.enc_params.tile_columns ||
                                av1enc_ctx_.enc_params.tile_rows);
#else
    EbBool has_tiles = (EbBool)EB_FALSE;
#endif
    obu_frame_header_size_ =
        has_tiles ? OBU_FRAME_HEADER_SIZE + 1 : OBU_FRAME_HEADER_SIZE;

    ASSERT_NE(psnr_src_, nullptr) << "PSNR source create failed!";
    EbErrorType err = psnr_src_->open_source(start_pos_, frames_to_test_);
    ASSERT_EQ(err, EB_ErrorNone) << "open_source return error:" << err;
    total_enc_out_ = 0;
}

void SvtAv1E2ETestFramework::close_test() {
    EbErrorType return_error = EB_ErrorNone;
    // Deinit
    return_error = eb_deinit_encoder(av1enc_ctx_.enc_handle);
    ASSERT_EQ(return_error, EB_ErrorNone)
        << "eb_deinit_encoder return error:" << return_error;
}

void SvtAv1E2ETestFramework::run_encode_process() {
    static const char READ_SRC[] = "read_src";
    static const char ENCODING[] = "encoding";
    static const char RECON[] = "recon";
    static const char CONFORMANCE[] = "conformance";

    EbErrorType return_error = EB_ErrorNone;

    uint32_t frame_count = video_src_->get_frame_count();
    ASSERT_GT(frame_count, 0) << "video srouce file does not contain frame!!";
    if (recon_queue_) {
        recon_queue_->set_frame_count(frame_count);
    }

    if (output_file_) {
        write_output_header();
    }

    uint8_t *frame = nullptr;
    bool src_file_eos = false;
    bool enc_file_eos = false;
    bool rec_file_eos = recon_queue_ ? false : true;
    do {
        if (!src_file_eos) {
            {
                TimeAutoCount counter(READ_SRC, collect_);
                frame = (uint8_t *)video_src_->get_next_frame();
            }
            {
                TimeAutoCount counter(ENCODING, collect_);
                if (frame != nullptr && frame_count) {
                    frame_count--;
                    // Fill in Buffers Header control data
                    av1enc_ctx_.input_picture_buffer->p_buffer = frame;
                    av1enc_ctx_.input_picture_buffer->n_filled_len =
                        video_src_->get_frame_size();
                    av1enc_ctx_.input_picture_buffer->flags = 0;
                    av1enc_ctx_.input_picture_buffer->p_app_private = nullptr;
                    av1enc_ctx_.input_picture_buffer->pts =
                        video_src_->get_frame_index();
                    av1enc_ctx_.input_picture_buffer->pic_type =
                        EB_AV1_INVALID_PICTURE;
                    // Send the picture
                    EXPECT_EQ(EB_ErrorNone,
                              return_error = eb_svt_enc_send_picture(
                                  av1enc_ctx_.enc_handle,
                                  av1enc_ctx_.input_picture_buffer))
                        << "eb_svt_enc_send_picture error at: "
                        << av1enc_ctx_.input_picture_buffer->pts;
                }
                if (frame_count == 0 || frame == nullptr) {
                    src_file_eos = true;
                    EbBufferHeaderType headerPtrLast;
                    headerPtrLast.n_alloc_len = 0;
                    headerPtrLast.n_filled_len = 0;
                    headerPtrLast.n_tick_count = 0;
                    headerPtrLast.p_app_private = nullptr;
                    headerPtrLast.flags = EB_BUFFERFLAG_EOS;
                    headerPtrLast.p_buffer = nullptr;
                    headerPtrLast.pic_type = EB_AV1_INVALID_PICTURE;
                    av1enc_ctx_.input_picture_buffer->flags = EB_BUFFERFLAG_EOS;
                    EXPECT_EQ(EB_ErrorNone,
                              return_error = eb_svt_enc_send_picture(
                                  av1enc_ctx_.enc_handle, &headerPtrLast))
                        << "eb_svt_enc_send_picture EOS error";
                }
            }
        }

        // recon
        if (recon_queue_ && !rec_file_eos) {
            TimeAutoCount counter(RECON, collect_);
            if (!rec_file_eos)
                get_recon_frame(av1enc_ctx_, recon_queue_, rec_file_eos);
        }

        if (!enc_file_eos) {
            do {
                // non-blocking call
                EbBufferHeaderType *enc_out = nullptr;
                {
                    TimeAutoCount counter(ENCODING, collect_);
                    int pic_send_done = (src_file_eos && rec_file_eos) ? 1 : 0;
                    return_error = eb_svt_get_packet(
                        av1enc_ctx_.enc_handle, &enc_out, pic_send_done);
                    ASSERT_NE(return_error, EB_ErrorMax)
                        << "Error while encoding, code:" << enc_out->flags;
                }

                // process the output buffer
                if (return_error != EB_NoErrorEmptyQueue && enc_out) {
                    TimeAutoCount counter(CONFORMANCE, collect_);
                    process_compress_data(enc_out);
                    if (enc_out->flags & EB_BUFFERFLAG_EOS) {
                        enc_file_eos = true;
                        printf("Encoder EOS\n");
                        break;
                    }
                } else {
                    if (return_error != EB_NoErrorEmptyQueue) {
                        enc_file_eos = true;
                        GTEST_FAIL() << "decoder return: " << return_error;
                    }
                    break;
                }

                // Release the output buffer
                if (enc_out != nullptr) {
                    eb_svt_release_out_buffer(&enc_out);
                    // EXPECT_EQ(enc_out, nullptr)
                    //    << "enc_out buffer is not well released";
                }
            } while (src_file_eos);
        }
    } while (!rec_file_eos || !src_file_eos || !enc_file_eos);

    /** complete the reference buffers in list comparison with recon */
    if (ref_compare_) {
        TimeAutoCount counter(CONFORMANCE, collect_);
        ASSERT_TRUE(ref_compare_->flush_video());
        delete ref_compare_;
        ref_compare_ = nullptr;
    }

    /** PSNR report */
    int count = 0;
    double psnr[4];
    pnsr_statistics_.get_statistics(count, psnr[0], psnr[1], psnr[2], psnr[3]);
    if (count > 0) {
        printf(
            "PSNR: %d frames, total: %0.4f, luma: %0.4f, cb: %0.4f, cr: "
            "%0.4f\r\n",
            count,
            psnr[0],
            psnr[1],
            psnr[2],
            psnr[3]);
    }
    pnsr_statistics_.reset();

    /** performance report */
    if (collect_) {
        frame_count = video_src_->get_frame_count();
        uint64_t total_enc_time = collect_->read_count(ENCODING);
        if (total_enc_time) {
            printf("Enc Performance: %.2fsec/frame (%.4fFPS)\n",
                   (double)total_enc_time / frame_count / 1000,
                   (double)frame_count * 1000 / total_enc_time);
        }
    }
}

void SvtAv1E2ETestFramework::write_output_header() {
    char header[IVF_STREAM_HEADER_SIZE];
    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4, 0);           // version
    mem_put_le16(header + 6, 32);          // header size
    mem_put_le32(header + 8, AV1_FOURCC);  // fourcc
    mem_put_le16(header + 12, av1enc_ctx_.enc_params.source_width);   // width
    mem_put_le16(header + 14, av1enc_ctx_.enc_params.source_height);  // height
    if (av1enc_ctx_.enc_params.frame_rate_denominator != 0 &&
        av1enc_ctx_.enc_params.frame_rate_numerator != 0) {
        mem_put_le32(header + 16,
                     av1enc_ctx_.enc_params.frame_rate_numerator);  // rate
        mem_put_le32(header + 20,
                     av1enc_ctx_.enc_params.frame_rate_denominator);  // scale
    } else {
        mem_put_le32(header + 16,
                     (av1enc_ctx_.enc_params.frame_rate >> 16) * 1000);  // rate
        mem_put_le32(header + 20, 1000);  // scale
    }
    mem_put_le32(header + 24, 0);  // length
    mem_put_le32(header + 28, 0);  // unused
    if (output_file_ && output_file_->file)
        fwrite(header, 1, IVF_STREAM_HEADER_SIZE, output_file_->file);
}

static void update_prev_ivf_header(
    svt_av1_e2e_test::SvtAv1E2ETestFramework::IvfFile *ivf) {
    char header[4];  // only for the number of bytes
    if (ivf && ivf->file && ivf->byte_count_since_ivf != 0) {
        fseeko64(
            ivf->file,
            (-(int32_t)(ivf->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE)),
            SEEK_CUR);
        mem_put_le32(&header[0], (int32_t)(ivf->byte_count_since_ivf));
        fwrite(header, 1, 4, ivf->file);
        fseeko64(ivf->file,
                 (ivf->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE - 4),
                 SEEK_CUR);
        ivf->byte_count_since_ivf = 0;
    }
}

static void write_ivf_frame_header(
    svt_av1_e2e_test::SvtAv1E2ETestFramework::IvfFile *ivf,
    uint32_t byte_count) {
    char header[IVF_FRAME_HEADER_SIZE];
    int32_t write_location = 0;

    mem_put_le32(&header[write_location], (int32_t)byte_count);
    write_location = write_location + 4;
    mem_put_le32(&header[write_location],
                 (int32_t)((ivf->ivf_count) & 0xFFFFFFFF));
    write_location = write_location + 4;
    mem_put_le32(&header[write_location], (int32_t)((ivf->ivf_count) >> 32));
    write_location = write_location + 4;

    ivf->byte_count_since_ivf = (byte_count);

    ivf->ivf_count++;
    fflush(stdout);

    if (ivf->file)
        fwrite(header, 1, IVF_FRAME_HEADER_SIZE, ivf->file);
}

void SvtAv1E2ETestFramework::write_compress_data(
    const EbBufferHeaderType *output) {
    // Check for the flags EB_BUFFERFLAG_HAS_TD and
    // EB_BUFFERFLAG_SHOW_EXT
    switch (output->flags & 0x00000006) {
    case (EB_BUFFERFLAG_HAS_TD | EB_BUFFERFLAG_SHOW_EXT):
        // terminate previous ivf packet, update the combined size of
        // packets sent
        update_prev_ivf_header(output_file_);

        // Write a new IVF frame header to file as a TD is in the packet
        write_ivf_frame_header(
            output_file_,
            output->n_filled_len - (obu_frame_header_size_ + TD_SIZE));
        fwrite(output->p_buffer,
               1,
               output->n_filled_len - (obu_frame_header_size_ + TD_SIZE),
               output_file_->file);

        // An EB_BUFFERFLAG_SHOW_EXT means that another TD has been added to
        // the packet to show another frame, a new IVF is needed
        write_ivf_frame_header(output_file_,
                               (obu_frame_header_size_ + TD_SIZE));
        fwrite(output->p_buffer + output->n_filled_len -
                   (obu_frame_header_size_ + TD_SIZE),
               1,
               (obu_frame_header_size_ + TD_SIZE),
               output_file_->file);

        break;
    case (EB_BUFFERFLAG_HAS_TD):
        // terminate previous ivf packet, update the combined size of
        // packets sent
        update_prev_ivf_header(output_file_);

        // Write a new IVF frame header to file as a TD is in the packet
        write_ivf_frame_header(output_file_, output->n_filled_len);
        fwrite(output->p_buffer, 1, output->n_filled_len, output_file_->file);
        break;
    case (EB_BUFFERFLAG_SHOW_EXT):
        // this case means that there's only one TD in this packet and is
        // relater
        fwrite(output->p_buffer,
               1,
               output->n_filled_len - (obu_frame_header_size_ + TD_SIZE),
               output_file_->file);
        // this packet will be part of the previous IVF header
        output_file_->byte_count_since_ivf +=
            (output->n_filled_len - (obu_frame_header_size_ + TD_SIZE));

        // terminate previous ivf packet, update the combined size of
        // packets sent
        update_prev_ivf_header(output_file_);

        // An EB_BUFFERFLAG_SHOW_EXT means that another TD has been added to
        // the packet to show another frame, a new IVF is needed
        write_ivf_frame_header(output_file_,
                               (obu_frame_header_size_ + TD_SIZE));
        fwrite(output->p_buffer + output->n_filled_len -
                   (obu_frame_header_size_ + TD_SIZE),
               1,
               (obu_frame_header_size_ + TD_SIZE),
               output_file_->file);

        break;
    default:
        // This is a packet without a TD, write it straight to file
        fwrite(output->p_buffer, 1, output->n_filled_len, output_file_->file);

        // this packet will be part of the previous IVF header
        output_file_->byte_count_since_ivf += (output->n_filled_len);
        break;
    }
}

void SvtAv1E2ETestFramework::process_compress_data(
    const EbBufferHeaderType *data) {
    ASSERT_NE(data, nullptr);
    if (refer_dec_ == nullptr) {
        if (output_file_) {
            write_compress_data(data);
        }
        return;
    }

    if (data->flags & EB_BUFFERFLAG_SHOW_EXT) {
        uint32_t first_part_size =
            data->n_filled_len - obu_frame_header_size_ - TD_SIZE;
        decode_compress_data(data->p_buffer, first_part_size);
        decode_compress_data(data->p_buffer + first_part_size,
                             obu_frame_header_size_ + TD_SIZE);
    } else {
        decode_compress_data(data->p_buffer, data->n_filled_len);
    }
}

void SvtAv1E2ETestFramework::decode_compress_data(const uint8_t *data,
                                                  const uint32_t size) {
    ASSERT_NE(data, nullptr);
    ASSERT_GT(size, 0);

    // input the compressed data into decoder
    ASSERT_EQ(refer_dec_->process_data(data, size), RefDecoder::REF_CODEC_OK);

    VideoFrame ref_frame;
    memset(&ref_frame, 0, sizeof(ref_frame));
    while (refer_dec_->get_frame(ref_frame) == RefDecoder::REF_CODEC_OK) {
        if (recon_queue_) {
            // compare tools
            if (ref_compare_ == nullptr) {
                ref_compare_ =
                    create_ref_compare_queue(ref_frame, recon_queue_);
                ASSERT_NE(ref_compare_, nullptr);
            }
            // Compare ref decode output with recon output.
            ASSERT_TRUE(ref_compare_->compare_video(ref_frame))
                << "image compare failed on " << ref_frame.timestamp;

            // PSNR tool
            check_psnr(ref_frame);
        }
    }
}

void SvtAv1E2ETestFramework::check_psnr(const VideoFrame &frame) {
    // Calculate psnr with input frame and
    EbSvtIOFormat *src_frame =
        psnr_src_->get_frame_by_index((const uint32_t)frame.timestamp);
    if (src_frame) {
        double luma_psnr = 0.0;
        double cb_psnr = 0.0;
        double cr_psnr = 0.0;
        psnr_frame(src_frame,
                   video_src_->get_bit_depth(),
                   frame,
                   luma_psnr,
                   cb_psnr,
                   cr_psnr);
        pnsr_statistics_.add(luma_psnr, cb_psnr, cr_psnr);
        // TODO: Check PSNR value reasonable here?
    }
}

static bool transfer_frame_planes(VideoFrame *frame) {
    uint32_t luma_size = frame->stride[0] * frame->height;
    if (frame->buf_size != 4 * luma_size) {
        printf("buffer size doesn't match!\n");
        return false;
    }
    bool half_height = false;
    if (frame->format == IMG_FMT_420P10_PACKED ||
        frame->format == IMG_FMT_422P10_PACKED ||
        frame->format == IMG_FMT_444P10_PACKED) {
        luma_size *= 2;
    }
    if (frame->format == IMG_FMT_420 ||
        frame->format == IMG_FMT_420P10_PACKED ||
        frame->format == IMG_FMT_NV12) {
        half_height = true;
    }
    if (frame->stride[3]) {
        uint32_t offset = frame->stride[0] * frame->height +
                          ((frame->stride[1] + frame->stride[2]) *
                           (half_height ? frame->height / 2 : frame->height));
        memcpy(frame->planes[3], frame->buffer + offset, frame->buf_size >> 2);
    }
    if (frame->stride[2]) {
        uint32_t offset = frame->stride[0] * frame->height +
                          (frame->stride[1] *
                           (half_height ? frame->height / 2 : frame->height));
        memcpy(frame->planes[2], frame->buffer + offset, frame->buf_size >> 2);
    }
    if (frame->stride[1]) {
        uint32_t offset = frame->stride[0] * frame->height;
        memcpy(frame->planes[1], frame->buffer + offset, frame->buf_size >> 2);
    }
    return true;
}

void SvtAv1E2ETestFramework::get_recon_frame(const SvtAv1Context &ctxt,
                                             FrameQueue *recon, bool &is_eos) {
    do {
        VideoFrame *new_frame = recon->get_empty_frame();
        ASSERT_NE(new_frame, nullptr)
            << "can not get new frame for recon frame!!";
        ASSERT_NE(new_frame->buffer, nullptr)
            << "can not get new buffer of recon frame!!";

        EbBufferHeaderType recon_frame = {0};
        recon_frame.size = sizeof(EbBufferHeaderType);
        recon_frame.p_buffer = new_frame->buffer;
        recon_frame.n_alloc_len = new_frame->buf_size;
        recon_frame.p_app_private = nullptr;
        // non-blocking call until all input frames are sent
        EbErrorType recon_status =
            eb_svt_get_recon(ctxt.enc_handle, &recon_frame);
        ASSERT_NE(recon_status, EB_ErrorMax)
            << "Error while outputing recon, code:" << recon_frame.flags;
        if (recon_status == EB_NoErrorEmptyQueue) {
            recon->delete_frame(new_frame);
            break;
        } else {
            ASSERT_LE(recon_frame.n_filled_len, new_frame->buf_size)
                << "recon frame size incorrect@" << recon_frame.pts;
            // mark the recon eos flag
            if (recon_frame.flags & EB_BUFFERFLAG_EOS)
                is_eos = true;
            transfer_frame_planes(new_frame);
            new_frame->timestamp = recon_frame.pts;
            recon->add_frame(new_frame);
        }
    } while (true);
}

SvtAv1E2ETestFramework::IvfFile::IvfFile(std::string path) {
    FOPEN(file, path.c_str(), "wb");
    byte_count_since_ivf = 0;
    ivf_count = 0;
}
