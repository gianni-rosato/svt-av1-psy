/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file DummyVideoSource.h
 *
 * @brief A dummy video source for testing different resolution without real
 *        source file.
 *        Dummy source will generate a color bars and keep moving to right side,
 *        FRAME_PER_LOOP defined the frame count in on loop.
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#ifndef _SVT_TEST_DUMMY_VIDEO_SOURCE_H_
#define _SVT_TEST_DUMMY_VIDEO_SOURCE_H_
#include "VideoSource.h"

namespace svt_av1_video_source {

#define FRAME_PER_LOOP 300

static const uint8_t color_bar_luma[8] = {255, 173, 131, 115, 96, 83, 34, 0};
static const uint8_t color_bar_cb[8] = {128, 52, 161, 87, 151, 80, 214, 128};
static const uint8_t color_bar_cr[8] = {128, 137, 34, 45, 198, 205, 134, 128};

class DummyVideoSource : public VideoSource {
  public:
    DummyVideoSource(const VideoColorFormat format, const uint32_t width,
                     const uint32_t height, const uint8_t bit_depth,
                     const bool use_compressed_2bit_plan_output)
        : VideoSource(format, width, height, bit_depth,
                      use_compressed_2bit_plan_output) {
        if (width_ % 16 != 0)
            width_with_padding_ = ((width_ >> 4) + 1) << 4;
        if (height_ % 16 != 0)
            height_with_padding_ = ((height_ >> 4) + 1) << 4;
    }

    virtual ~DummyVideoSource() {
        deinit_frame_buffer();
    }

    EbErrorType open_source(const uint32_t init_pos,
                            const uint32_t frame_count) override {
        if (image_format_ != IMG_FMT_420) {
            printf("Open dummy source error, support YUV420 8bit only\r\n");
            return EB_ErrorBadParameter;
        }

        init_pos_ = init_pos;
        if (frame_count == 0)
            frame_count_ = FRAME_PER_LOOP;
        else
            frame_count_ = frame_count;

        current_frame_index_ = -1;
        init_frame_buffer();

        return EB_ErrorNone;
    }

    /*!\brief Close stream. */
    void close_source() override {
        deinit_frame_buffer();
    }

    /*!\brief Get next frame. */
    EbSvtIOFormat *get_next_frame() override {
        if ((uint32_t)(current_frame_index_ + 1) >= frame_count_)
            return nullptr;
        current_frame_index_++;
        generate_frame(current_frame_index_ + init_pos_);
        return frame_buffer_;
    }

    /*!\brief Get frame by index. */
    EbSvtIOFormat *get_frame_by_index(const uint32_t index) override {
        if (index >= frame_count_)
            return nullptr;
        current_frame_index_ = index;
        generate_frame(current_frame_index_ + init_pos_);
        return frame_buffer_;
    }

  protected:
    void generate_frame(const uint32_t index) {
        // generate a color bar pattern, moving with FRAME_PER_LOOP frame as one
        // loop.
        const uint32_t single_width = width_with_padding_ / 8;
        uint8_t *src_p = nullptr;
        uint8_t *p = nullptr;
        uint32_t offset = 0;

        offset =
            (index % FRAME_PER_LOOP) * width_with_padding_ / FRAME_PER_LOOP;
        offset -= offset % 2;
        // Support yuv420 8bit only.
        // luma
        src_p = frame_buffer_->luma;
        p = src_p + offset;

        for (int r = 0; r < 8; r++) {
            memset(p, color_bar_luma[r], single_width);
            p += single_width;
        }
        memcpy(src_p, src_p + width_with_padding_, offset);
        p = src_p;
        for (uint32_t l = 1; l < height_with_padding_; l++) {
            memcpy(p, src_p, width_with_padding_);
            p += width_with_padding_;
        }

        // cb
        src_p = frame_buffer_->cb;
        p = src_p + offset / 2;
        for (int r = 0; r < 8; r++) {
            memset(p, color_bar_cb[r], single_width / 2);
            p += single_width / 2;
        }
        memcpy(src_p, src_p + width_with_padding_ / 2, offset / 2);
        p = src_p;
        for (uint32_t l = 1; l < height_with_padding_ / 2; l++) {
            memcpy(p, src_p, width_with_padding_ / 2);
            p += width_with_padding_ / 2;
        }

        // cr
        src_p = frame_buffer_->cr;
        p = src_p + offset / 2;
        for (int r = 0; r < 8; r++) {
            memset(p, color_bar_cr[r], single_width / 2);
            p += single_width / 2;
        }
        memcpy(src_p, src_p + width_with_padding_ / 2, offset / 2);
        p = src_p;
        for (uint32_t l = 0; l < height_with_padding_ / 2; l++) {
            memcpy(p, src_p, width_with_padding_ / 2);
            p += width_with_padding_ / 2;
        }
    }
};
}  // namespace svt_av1_video_source
#endif  //_SVT_TEST_DUMMY_VIDEO_SOURCE_H_
