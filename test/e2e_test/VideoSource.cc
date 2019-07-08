/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file VideoSource.cc
 *
 * @brief Impelmentation class for video source & video file source.
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#include <stdio.h>
#include "EbDefinitions.h"
#include "VideoSource.h"
#include "random.h"

using namespace svt_av1_video_source;
VideoSource::VideoSource(const VideoColorFormat format, const uint32_t width,
                         const uint32_t height, const uint8_t bit_depth,
                         const bool use_compressed_2bit_plane_output)
    : width_(width),
      height_(height),
      width_with_padding_(width),
      height_with_padding_(height),
      bit_depth_(bit_depth),
      init_pos_(0),
      frame_count_(0),
      current_frame_index_(-1),
      frame_size_(0),
      frame_buffer_(nullptr),
      image_format_(format) {
    if (bit_depth_ > 8 && use_compressed_2bit_plane_output)
        svt_compressed_2bit_plane_ = true;
    else
        svt_compressed_2bit_plane_ = false;
    frame_qp_list_.clear();
}

VideoSource::~VideoSource() {
    deinit_frame_buffer();
}

bool VideoSource::is_10bit_mode() {
    if (image_format_ == IMG_FMT_420P10_PACKED ||
        image_format_ == IMG_FMT_422P10_PACKED ||
        image_format_ == IMG_FMT_444P10_PACKED) {
        return true;
    }
    return false;
}

void VideoSource::deinit_frame_buffer() {
    if (frame_buffer_ == nullptr)
        return;

    if (frame_buffer_->luma != nullptr) {
        free(frame_buffer_->luma);
        frame_buffer_->luma = nullptr;
    }
    if (frame_buffer_->cb != nullptr) {
        free(frame_buffer_->cb);
        frame_buffer_->cb = nullptr;
    }
    if (frame_buffer_->cr != nullptr) {
        free(frame_buffer_->cr);
        frame_buffer_->cr = nullptr;
    }
    if (frame_buffer_->luma_ext != nullptr) {
        free(frame_buffer_->luma_ext);
        frame_buffer_->luma_ext = nullptr;
    }
    if (frame_buffer_->cb_ext != nullptr) {
        free(frame_buffer_->cb_ext);
        frame_buffer_->cb_ext = nullptr;
    }
    if (frame_buffer_->cr_ext != nullptr) {
        free(frame_buffer_->cr_ext);
        frame_buffer_->cr_ext = nullptr;
    }
    free(frame_buffer_);
    frame_buffer_ = nullptr;
}
EbErrorType VideoSource::init_frame_buffer() {
    // Determine size of each plane
    uint32_t luma_size = width_with_padding_ * height_with_padding_;
    uint32_t chroma_size = 0;

    switch (image_format_) {
    case IMG_FMT_420P10_PACKED:
    case IMG_FMT_420: {
        chroma_size = luma_size >> 2;
    } break;
    case IMG_FMT_422P10_PACKED:
    case IMG_FMT_422: {
        chroma_size = luma_size >> 1;
    } break;
    case IMG_FMT_444P10_PACKED:
    case IMG_FMT_444: {
        chroma_size = luma_size;
    } break;
    default: { chroma_size = luma_size >> 2; } break;
    }

    // Determine
    if (frame_buffer_ == nullptr)
        frame_buffer_ = (EbSvtIOFormat *)malloc(sizeof(EbSvtIOFormat));

    if (frame_buffer_ == nullptr) {
        deinit_frame_buffer();
        return EB_ErrorInsufficientResources;
    }

    memset(frame_buffer_, 0, sizeof(EbSvtIOFormat));
    frame_buffer_->width = width_with_padding_;
    frame_buffer_->height = height_with_padding_;
    frame_buffer_->origin_x = 0;
    frame_buffer_->origin_y = 0;

    // SVT-AV1 use pixel size as stride?
    frame_buffer_->y_stride = luma_size;
    frame_buffer_->cb_stride = chroma_size;
    frame_buffer_->cr_stride = chroma_size;

    if (is_10bit_mode() && !svt_compressed_2bit_plane_) {
        luma_size *= 2;
        chroma_size *= 2;
    }

    frame_buffer_->luma = (uint8_t *)malloc(luma_size);
    if (frame_buffer_->luma == nullptr) {
        deinit_frame_buffer();
        return EB_ErrorInsufficientResources;
    }

    frame_buffer_->cb = (uint8_t *)malloc(chroma_size);
    if (frame_buffer_->cb == nullptr) {
        deinit_frame_buffer();
        return EB_ErrorInsufficientResources;
    }

    frame_buffer_->cr = (uint8_t *)malloc(chroma_size);
    if (frame_buffer_->cr == nullptr) {
        deinit_frame_buffer();
        return EB_ErrorInsufficientResources;
    }

    if (is_10bit_mode() && svt_compressed_2bit_plane_) {
        frame_buffer_->luma_ext = (uint8_t *)malloc(luma_size / 4);
        if (frame_buffer_->luma_ext == nullptr) {
            deinit_frame_buffer();
            return EB_ErrorInsufficientResources;
        }

        frame_buffer_->cb_ext = (uint8_t *)malloc(chroma_size / 4);
        if (frame_buffer_->cb_ext == nullptr) {
            deinit_frame_buffer();
            return EB_ErrorInsufficientResources;
        }

        frame_buffer_->cr_ext = (uint8_t *)malloc(chroma_size / 4);
        if (frame_buffer_->cr_ext == nullptr) {
            deinit_frame_buffer();
            return EB_ErrorInsufficientResources;
        }
    } else {
        frame_buffer_->luma_ext = nullptr;
        frame_buffer_->cb_ext = nullptr;
        frame_buffer_->cr_ext = nullptr;
    }

    return EB_ErrorNone;
}

VideoFileSource::VideoFileSource(const std::string &file_name,
                                 const VideoColorFormat format,
                                 const uint32_t width, const uint32_t height,
                                 const uint8_t bit_depth,
                                 const bool use_compressed_2bit_plane_output)
    : VideoSource(format, width, height, bit_depth,
                  use_compressed_2bit_plane_output),
      file_name_(file_name),
      file_handle_(nullptr),
      file_length_(0),
      file_frames_(0) {
    if (width_ % 8 != 0)
        width_with_padding_ += (8 - width_ % 8);
    if (height_ % 8 != 0)
        height_with_padding_ += (8 - height_ % 8);
}

VideoFileSource::~VideoFileSource() {
}
/**
 * @brief      Use this function to get vector path defined by envrionment
 * variable SVT_AV1_TEST_VECTOR_PATH, or it will return a default path.
 *
 * @return     The vectors path.
 */
std::string VideoFileSource::get_vector_dir() {
    const char *const path = getenv("SVT_AV1_TEST_VECTOR_PATH");
    if (path == nullptr) {
#ifdef _WIN32
        return "../../../../test/vectors";
#else
        return "../../test/vectors";
#endif  // _WIN32
    }
    return path;
}

uint32_t VideoFileSource::read_input_frame() {
    uint32_t filled_len = 0;
    if (file_handle_ == nullptr) {
        printf("Error file handle\r\n");
        return 0;
    }

    if (feof(file_handle_) != 0) {
        printf("Reach file end\r\n");
        return 0;
    }
    int width_downsize = 1;
    int height_downsize = 1;
    int pixel_byte_size = 1;

    switch (image_format_) {
    case IMG_FMT_420: {
        width_downsize = 1;
        height_downsize = 1;
        pixel_byte_size = 1;
        break;
    }
    case IMG_FMT_422: {
        width_downsize = 1;
        height_downsize = 0;
        pixel_byte_size = 1;
        break;
    }
    case IMG_FMT_444: {
        width_downsize = 0;
        height_downsize = 0;
        pixel_byte_size = 1;
        break;
    }
    case IMG_FMT_420P10_PACKED: {
        width_downsize = 1;
        height_downsize = 1;
        pixel_byte_size = 2;
        break;
    }
    case IMG_FMT_422P10_PACKED: {
        width_downsize = 1;
        height_downsize = 0;
        pixel_byte_size = 2;
        break;
    }
    case IMG_FMT_444P10_PACKED: {
        width_downsize = 0;
        height_downsize = 0;
        pixel_byte_size = 2;
        break;
    }
    default: break;
    }

    // SVT-AV1 use pixel size as stride?
    frame_buffer_->y_stride = width_with_padding_;
    frame_buffer_->cb_stride = (width_with_padding_ >> width_downsize);
    frame_buffer_->cr_stride = (width_with_padding_ >> width_downsize);

    // Read raw data from file
    const uint32_t bottom_padding = height_with_padding_ - height_;
    const uint32_t righ_padding = width_with_padding_ - width_;
    size_t read_len = 0;
    uint32_t i;
    if (bit_depth_ <= 8 || (bit_depth_ > 8 && !svt_compressed_2bit_plane_)) {
        uint8_t *eb_input_ptr = nullptr;
        // Y
        eb_input_ptr = frame_buffer_->luma;
        for (i = 0; i < height_; ++i) {
            read_len =
                fread(eb_input_ptr, 1, width_ * pixel_byte_size, file_handle_);
            if (read_len != width_ * pixel_byte_size)
                return 0;  // read file error.

            memset(eb_input_ptr + width_ * pixel_byte_size,
                   0,
                   righ_padding * pixel_byte_size);
            eb_input_ptr += frame_buffer_->y_stride * pixel_byte_size;
            filled_len += frame_buffer_->y_stride * pixel_byte_size;
        }
        for (i = 0; i < bottom_padding; ++i) {
            memset(eb_input_ptr, 0, width_with_padding_ * pixel_byte_size);
            eb_input_ptr += frame_buffer_->y_stride * pixel_byte_size;
            filled_len += frame_buffer_->y_stride * pixel_byte_size;
        }
        // Cb
        eb_input_ptr = frame_buffer_->cb;
        for (i = 0; i < (height_ >> height_downsize); ++i) {
            read_len = fread(eb_input_ptr,
                             1,
                             (width_ >> width_downsize) * pixel_byte_size,
                             file_handle_);
            if (read_len != (width_ >> width_downsize) * pixel_byte_size)
                return 0;  // read file error.

            memset(eb_input_ptr + (width_ >> width_downsize) * pixel_byte_size,
                   0,
                   (righ_padding >> width_downsize) * pixel_byte_size);
            eb_input_ptr += frame_buffer_->cb_stride * pixel_byte_size;
            filled_len += frame_buffer_->cb_stride * pixel_byte_size;
        }
        for (i = 0; i < (bottom_padding >> height_downsize); ++i) {
            memset(eb_input_ptr,
                   0,
                   (width_with_padding_ >> width_downsize) * pixel_byte_size);
            eb_input_ptr += frame_buffer_->cb_stride * pixel_byte_size;
            filled_len += frame_buffer_->cb_stride * pixel_byte_size;
        }

        // Cr
        eb_input_ptr = frame_buffer_->cr;

        for (i = 0; i < (height_ >> height_downsize); ++i) {
            read_len = fread(eb_input_ptr,
                             1,
                             (width_ >> width_downsize) * pixel_byte_size,
                             file_handle_);
            if (read_len != (width_ >> width_downsize) * pixel_byte_size)
                return 0;  // read file error.

            memset(eb_input_ptr + (width_ >> width_downsize) * pixel_byte_size,
                   0,
                   (righ_padding >> width_downsize) * pixel_byte_size);
            eb_input_ptr += frame_buffer_->cr_stride * pixel_byte_size;
            filled_len += frame_buffer_->cr_stride * pixel_byte_size;
        }
        for (i = 0; i < (bottom_padding >> height_downsize); ++i) {
            memset(eb_input_ptr,
                   0,
                   (width_with_padding_ >> width_downsize) * pixel_byte_size);
            eb_input_ptr += frame_buffer_->cr_stride * pixel_byte_size;
            filled_len += frame_buffer_->cr_stride * pixel_byte_size;
        }
    } else if (bit_depth_ > 8 && svt_compressed_2bit_plane_) {
        uint8_t *eb_input_ptr = nullptr;
        uint8_t *eb_ext_input_ptr = nullptr;
        // Y
        uint16_t pix = 0;
        eb_input_ptr = frame_buffer_->luma;
        eb_ext_input_ptr = frame_buffer_->luma_ext;
        for (i = 0; i < height_; ++i) {
            uint32_t j = 0;
            for (j = 0; j < width_; ++j) {
                // Get one pixel
                if (2 != fread(&pix, 1, 2, file_handle_)) {
                    return 0;
                }
                eb_input_ptr[j] = (uint8_t)(pix >> 2);
                eb_ext_input_ptr[j / 4] &= (pix & 0x3) << (j % 4);
            }

            for (; j < width_ + righ_padding; ++j) {
                eb_input_ptr[j] = 0;
                eb_ext_input_ptr[j / 4] = 0;
            }
            eb_input_ptr += frame_buffer_->y_stride;
            eb_ext_input_ptr += frame_buffer_->y_stride / 4;
            filled_len += frame_buffer_->y_stride * 5 / 4;
        }
        for (i = 0; i < bottom_padding; ++i) {
            memset(eb_input_ptr, 0, width_with_padding_);
            memset(eb_ext_input_ptr, 0, width_with_padding_ / 4);
            eb_input_ptr += frame_buffer_->y_stride;
            eb_ext_input_ptr += frame_buffer_->y_stride / 4;
            filled_len += frame_buffer_->y_stride * 5 / 4;
        }
        // Cb
        eb_input_ptr = frame_buffer_->cb;
        eb_ext_input_ptr = frame_buffer_->cb_ext;
        for (i = 0; i < (height_ >> height_downsize); ++i) {
            uint32_t j = 0;
            for (j = 0; j < (width_ >> width_downsize); ++j) {
                // Get one pixel
                if (2 != fread(&pix, 1, 2, file_handle_)) {
                    return 0;
                }
                eb_input_ptr[j] = (uint8_t)(pix >> 2);
                eb_ext_input_ptr[j / 4] &= (pix & 0x3) << (j % 4);
            }

            for (; j < ((width_ + righ_padding) >> width_downsize); ++j) {
                eb_input_ptr[j] = 0;
                eb_ext_input_ptr[j / 4] = 0;
            }
            eb_input_ptr += frame_buffer_->cb_stride;
            eb_ext_input_ptr += frame_buffer_->cb_stride / 4;
            filled_len += frame_buffer_->cb_stride * 5 / 4;
        }
        for (i = 0; i<bottom_padding>> height_downsize; ++i) {
            memset(eb_input_ptr, 0, (width_with_padding_ >> width_downsize));
            memset(eb_ext_input_ptr,
                   0,
                   (width_with_padding_ >> width_downsize) / 4);
            eb_input_ptr += frame_buffer_->cb_stride;
            eb_ext_input_ptr += frame_buffer_->cb_stride / 4;
            filled_len += frame_buffer_->cb_stride * 5 / 4;
        }

        // Cr
        eb_input_ptr = frame_buffer_->cr;
        eb_ext_input_ptr = frame_buffer_->cr_ext;
        for (i = 0; i < (height_ >> height_downsize); ++i) {
            uint32_t j = 0;
            for (j = 0; j < (width_ >> width_downsize); ++j) {
                // Get one pixel
                if (2 != fread(&pix, 1, 2, file_handle_)) {
                    return 0;
                }
                eb_input_ptr[j] = (uint8_t)(pix >> 2);
                eb_ext_input_ptr[j / 4] &= (pix & 0x3) << (j % 4);
            }

            for (; j < ((width_ + righ_padding) >> width_downsize); ++j) {
                eb_input_ptr[j] = 0;
                eb_ext_input_ptr[j / 4] = 0;
            }
            eb_input_ptr += frame_buffer_->cr_stride;
            eb_ext_input_ptr += frame_buffer_->cr_stride / 4;
            filled_len += frame_buffer_->cr_stride * 5 / 4;
        }
        for (i = 0; i<bottom_padding>> height_downsize; ++i) {
            memset(eb_input_ptr, 0, (width_with_padding_ >> width_downsize));
            memset(eb_ext_input_ptr,
                   0,
                   (width_with_padding_ >> width_downsize) / 4);
            eb_input_ptr += frame_buffer_->cr_stride;
            eb_ext_input_ptr += frame_buffer_->cr_stride / 4;
            filled_len += frame_buffer_->cr_stride * 5 / 4;
        }
    }

    return filled_len;
}

/*!\brief Prepare stream. */
EbErrorType VideoFileSource::open_source(const uint32_t init_pos,
                                         const uint32_t frame_count) {
    EbErrorType return_error = EB_ErrorNone;
    if (file_handle_ != nullptr)
        return EB_ErrorNone;
    std::string full_path = get_vector_dir() + "/" + file_name_.c_str();

    FOPEN(file_handle_, full_path.c_str(), "rb");
    if (file_handle_ == nullptr) {
        printf(">>> Open video source %s failed!\r\n", full_path.c_str());
        printf(
            "    You can use CMake generated build target TestVectors to get "
            "correct test vectors before run this test.\r\n");
        printf(
            "    For unix like system, run 'make TestVectors', or build "
            "TestVectors project on "
            "VisualStudio.\r\n");
        return EB_ErrorBadParameter;
    }

    // Seek to begin
    fseek(file_handle_, 0, SEEK_SET);

    // Get file info before prepare buffer
    return_error = parse_file_info();
    if (return_error != EB_ErrorNone) {
        fclose(file_handle_);
        file_handle_ = nullptr;
        return return_error;
    }

    // Prepare buffer
    if (EB_ErrorNone != init_frame_buffer()) {
        fclose(file_handle_);
        file_handle_ = nullptr;
        return EB_ErrorInsufficientResources;
    }
    if (file_frames_ <= init_pos || init_pos + frame_count >= file_frames_) {
        printf(
            "setup of initial position(%u) and output frame count(%u) is out "
            "of bound!\n",
            init_pos,
            frame_count);
        fclose(file_handle_);
        file_handle_ = nullptr;
        return EB_ErrorInsufficientResources;
    }
    init_pos_ = init_pos;
    if (frame_count == 0)
        frame_count_ = file_frames_ - init_pos_;
    else
        frame_count_ = (file_frames_ - init_pos_) > frame_count
                           ? frame_count
                           : (file_frames_ - init_pos_);

    // generate frame qp from random if failed from file
    if (true /* TODO: first from qp file */) {
        svt_av1_test_tool::SVTRandom rnd(1, 63);
        for (uint32_t i = 0; i < frame_count_; i++) {
            frame_qp_list_.push_back(rnd.random());
        }
    }

    if (seek_to_frame(init_pos_) != EB_ErrorNone) {
        fclose(file_handle_);
        file_handle_ = nullptr;
        return EB_ErrorInsufficientResources;
    }

    current_frame_index_ = -1;
#ifdef ENABLE_DEBUG_MONITOR
    if (monitor == nullptr) {
        monitor = new VideoMonitor(
            width_with_padding_,
            height_with_padding_,
            (bit_depth_ > 8) ? width_with_padding_ * 2 : width_with_padding_,
            bit_depth_,
            svt_compressed_2bit_plane_,
            "Y4M Source");
    }
#endif
    return EB_ErrorNone;
}
/*!\brief Close stream. */
void VideoFileSource::close_source() {
    if (file_handle_) {
        fclose(file_handle_);
        file_handle_ = nullptr;
    }

#ifdef ENABLE_DEBUG_MONITOR
    if (monitor) {
        delete monitor;
        monitor = nullptr;
    }
#endif
    deinit_frame_buffer();
    init_pos_ = 0;
    frame_count_ = 0;
    file_frames_ = 0;
}

EbSvtIOFormat *VideoFileSource::get_frame_by_index(const uint32_t index) {
    if (index >= frame_count_) {
        return nullptr;
    }
    // Seek to frame by index
    if (seek_to_frame(index) != EB_ErrorNone)
        return nullptr;

    frame_size_ = read_input_frame();
    if (frame_size_ == 0)
        return nullptr;

    current_frame_index_ = index;
#ifdef ENABLE_DEBUG_MONITOR
    if (monitor != nullptr) {
        monitor->draw_frame(
            frame_buffer_->luma, frame_buffer_->cb, frame_buffer_->cr);
    }
#endif
    return frame_buffer_;
}

EbSvtIOFormat *VideoFileSource::get_next_frame() {
    // printf("Get Next Frame:%d\r\n", current_frame_index_ + 1);
    frame_size_ = read_input_frame();
    if (frame_size_ == 0)
        return nullptr;
    ++current_frame_index_;

#ifdef ENABLE_DEBUG_MONITOR
    if (monitor != nullptr) {
        monitor->draw_frame(
            frame_buffer_->luma, frame_buffer_->cb, frame_buffer_->cr);
    }
#endif
    return frame_buffer_;
}
