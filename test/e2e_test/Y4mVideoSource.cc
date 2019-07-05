/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file Y4mVideoSource.h
 *
 * @brief Impelmentation of Y4mVideoSource class for reading y4m file.
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include "Y4mVideoSource.h"
extern "C" {
#include "EbAppInputy4m.h"
}
using namespace svt_av1_video_source;

Y4MVideoSource::Y4MVideoSource(const std::string& file_name,
                               const VideoColorFormat format,
                               const uint32_t width, const uint32_t height,
                               const uint8_t bit_depth,
                               const bool use_compressed_2bit_plane_output)
    : VideoFileSource(file_name, format, width, height, bit_depth,
                      use_compressed_2bit_plane_output) {
    frame_length_ = 0;
    header_length_ = 0;
#ifdef ENABLE_DEBUG_MONITOR
    monitor = nullptr;
#endif
}

Y4MVideoSource::~Y4MVideoSource() {
    if (file_handle_ != nullptr) {
        fclose(file_handle_);
    }
#ifdef ENABLE_DEBUG_MONITOR
    if (monitor != nullptr)
        delete monitor;
#endif
}

#define SKIP_TAG                                      \
    {                                                 \
        char tmp;                                     \
        do {                                          \
            if (1 != fread(&tmp, 1, 1, file_handle_)) \
                break;                                \
        } while ((tmp != 0xA) && (tmp != ' '));       \
    }

EbErrorType Y4MVideoSource::parse_file_info() {
    if (file_handle_ == nullptr)
        return EB_ErrorBadParameter;

    // Get file length
    fseek(file_handle_, 0, SEEK_END);
    file_length_ = ftell(file_handle_);

    // Seek to begin
    fseek(file_handle_, 0, SEEK_SET);

    EbConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.input_file = file_handle_;
    if (check_if_y4m(&cfg) != EB_TRUE)
        return EB_ErrorBadParameter;
    read_y4m_header(&cfg);

    width_ = cfg.source_width;
    height_ = cfg.source_height;
    bit_depth_ = cfg.encoder_bit_depth;

    // Workaround, check_if_y4m can not output color format
    // and svt-av1 encoder support yuv420 color format only in current version,
    // so use bit_depth to get iamge format.
    if (bit_depth_ > 8) {
        image_format_ = IMG_FMT_420P10_PACKED;
    } else {
        image_format_ = IMG_FMT_420;
    }

    // Get header lenght
    header_length_ = ftell(file_handle_);
    frame_length_ = width_ * height_;

    // Calculate frame length
    switch (image_format_) {
    case IMG_FMT_420P10_PACKED:
    case IMG_FMT_420: {
        frame_length_ = frame_length_ * 3 / 2 * (bit_depth_ > 8 ? 2 : 1);
    } break;
    case IMG_FMT_422P10_PACKED:
    case IMG_FMT_422: {
        frame_length_ = frame_length_ * 2 * (bit_depth_ > 8 ? 2 : 1);
    } break;
    case IMG_FMT_444P10_PACKED:
    case IMG_FMT_444: {
        frame_length_ = frame_length_ * 3 * (bit_depth_ > 8 ? 2 : 1);
    } break;
    default: break;
    }
    frame_length_ += 6;  // FRAME header

    // Calculate frame count
    file_frames_ = (file_length_ - header_length_) / frame_length_;

    printf("File len:%d; frame w:%d, h:%d, frame len:%d; frame count:%d\r\n",
           file_length_,
           width_,
           height_,
           frame_length_,
           file_frames_);

    return EB_ErrorNone;
}

uint32_t Y4MVideoSource::read_input_frame() {
    char frame_header[6] = {0};
    if (file_handle_ == nullptr) {
        printf("Error file handle\r\n");
        return 0;
    }

    if (feof(file_handle_) != 0) {
        printf("Reach file end\r\n");
        return 0;
    }

    if (6 != fread(frame_header, 1, 6, file_handle_)) {
        printf("can not found frame header\r\n");
        return 0;
    }

    // Check frame header
    if (!((strncmp("FRAME", frame_header, 5) == 0) && frame_header[5] == 0xA)) {
        printf("Read frame error\n");
        return 0;
    }
    return VideoFileSource::read_input_frame();
}

EbErrorType Y4MVideoSource::seek_to_frame(const uint32_t index) {
    if (file_handle_ == nullptr) {
        printf("Error file handle\r\n");
        return EB_ErrorInsufficientResources;
    }
    fseek(file_handle_, index * frame_length_ + header_length_, SEEK_SET);
    return EB_ErrorNone;
}
