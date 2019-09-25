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
    src_name_ = "YUV Source";
}

YuvVideoSource::~YuvVideoSource() {
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

    cal_yuv_plane_param();

    // Get file length
    fseek(file_handle_, 0, SEEK_END);
    file_length_ = ftell(file_handle_);

    // Seek to begin
    fseek(file_handle_, 0, SEEK_SET);

    uint32_t luma_size = width_ * height_ * bytes_per_sample_;
    frame_length_ =
        luma_size + ((luma_size >> (width_downsize_ + height_downsize_)) * 2);

    // Calculate frame count
    file_frames_ = file_length_ / frame_length_;
    return EB_ErrorNone;
}

EbErrorType YuvVideoSource::seek_to_frame(const uint32_t index) {
    if (file_handle_ == nullptr)
        return EB_ErrorBadParameter;
    uint32_t real_index = init_pos_ + index;
    if (real_index >= file_frames_)
        real_index = real_index % file_frames_;
    if (fseek(file_handle_, real_index * frame_length_, SEEK_SET) != 0)
        return EB_ErrorInsufficientResources;
    return EB_ErrorNone;
}
