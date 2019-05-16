/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file Y4mVideoSource.h
 *
 * @brief Impelmentation of YuvVideoSource class for reading y4m file.
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#include "YuvVideoSource.h"
using namespace svt_av1_video_source;
#define SIZE_OF_ONE_FRAME_IN_BYTES(width, height) (((width) * (height)*3) >> 1)

YuvVideoSource::YuvVideoSource(const std::string &file_name,
                               const VideoColorFormat format,
                               const uint32_t width, const uint32_t height,
                               const uint8_t bit_depth)
    : VideoFileSource(file_name, format, width, height, bit_depth, false) {
#ifdef ENABLE_DEBUG_MONITOR
    monitor = nullptr;
#endif
}
YuvVideoSource::~YuvVideoSource() {
    if (file_handle_ != nullptr) {
        fclose(file_handle_);
        file_handle_ = nullptr;
    }
#ifdef ENABLE_DEBUG_MONITOR
    if (monitor != nullptr)
        delete monitor;
#endif
}

EbErrorType YuvVideoSource::parse_file_info() {
    if (file_handle_ == nullptr)
        return EB_ErrorBadParameter;

    // Prepare buffer
    if (EB_ErrorNone != init_frame_buffer()) {
        fclose(file_handle_);
        file_handle_ = nullptr;
        return EB_ErrorInsufficientResources;
    }

    // Get file length
    fseek(file_handle_, 0, SEEK_END);
    file_length_ = ftell(file_handle_);

    // Seek to begin
    fseek(file_handle_, 0, SEEK_SET);

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

    // Calculate frame count
    file_frames_ = file_length_ / frame_length_;
    return EB_ErrorNone;
}

EbErrorType YuvVideoSource::seek_to_frame(const uint32_t index) {
    if (file_handle_ == nullptr)
        return EB_ErrorBadParameter;
    if (fseek(file_handle_, (init_pos_ + index) * frame_length_, SEEK_SET) != 0)
        return EB_ErrorInsufficientResources;
    return EB_ErrorNone;
}
