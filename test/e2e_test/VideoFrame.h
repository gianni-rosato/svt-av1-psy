/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file VideoFrame.h
 *
 * @brief Defines video frame types and structure
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/
#ifndef _SVT_TEST_VIDEO_FRAME_H_
#define _SVT_TEST_VIDEO_FRAME_H_

#include <memory.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#define INVALID_QP (0xFF)

/** VideoColorFormat defines the format of YUV video */
typedef enum VideoColorFormat {
    IMG_FMT_YV12,
    IMG_FMT_420 = IMG_FMT_YV12,
    IMG_FMT_422,
    IMG_FMT_444,
    IMG_FMT_420P10_PACKED,
    IMG_FMT_422P10_PACKED,
    IMG_FMT_444P10_PACKED,
    IMG_FMT_NV12,
    IMG_FMT_YV12_CUSTOM_COLOR_SPACE,
    IMG_FMT_NV12_CUSTOM_COLOR_SPACE,
    IMG_FMT_444A,
} VideoColorFormat;

static inline bool is_ten_bits(VideoColorFormat fmt) {
    switch (fmt) {
    case IMG_FMT_420P10_PACKED:
    case IMG_FMT_422P10_PACKED:
    case IMG_FMT_444P10_PACKED: return true;
    default: break;
    }
    return false;
}

/** VideoFrameParam defines the basic parameters of video frame */
typedef struct VideoFrameParam {
    VideoColorFormat format; /**< video format type */
    uint32_t width;          /**< width of video frame in pixel */
    uint32_t height;         /**< height of video frame in pixel */
    VideoFrameParam() {
        format = IMG_FMT_YV12;
        width = 0;
        height = 0;
    }
} VideoFrameParam;

/** VideoFrame defines the full parameters of video frame */
typedef struct VideoFrame : public VideoFrameParam {
    uint32_t disp_width;      /**< width to display*/
    uint32_t disp_height;     /**< height to display*/
    uint32_t stride[4];       /**< stride in array */
    uint8_t *planes[4];       /**< plane pointer to buffer address*/
    uint32_t bits_per_sample; /**< bit for each sample, default is 8, more for
                                 packed formats */
    void *context;            /**< private data context from creator */
    uint64_t timestamp;       /**< timestamp(index) of this frame */
    uint8_t *buffer;          /**< self own buffer */
    uint32_t buf_size;        /**< buffer size in bytes */
    uint32_t qp;              /**< qp of this frame */
    VideoFrame() {
        disp_width = 0;
        disp_height = 0;
        memset(&stride, 0, sizeof(stride));
        memset(&planes, 0, sizeof(planes));
        bits_per_sample = 8;
        context = nullptr;
        timestamp = 0;
        buffer = nullptr;
        buf_size = 0;
        qp = INVALID_QP;
    }
    VideoFrame(const VideoFrameParam &param) {
        // copy basic info from param
        *(VideoFrameParam *)this = param;
        disp_width = param.width;
        disp_height = param.height;
        memset(&stride, 0, sizeof(stride));
        memset(&planes, 0, sizeof(planes));
        bits_per_sample = is_ten_bits(param.format) ? 10 : 8;
        context = nullptr;
        timestamp = 0;
        buffer = nullptr;
        buf_size = 0;
        qp = INVALID_QP;

        // allocate memory for new frame
        uint32_t max_size = calculate_max_frame_size(param);
        buffer = new uint8_t[max_size];
        if (buffer) {
            memset(buffer, 0, max_size);
            planes[0] = buffer;
            planes[1] = buffer + (max_size >> 2);
            planes[2] = buffer + (max_size >> 1);
            planes[3] = buffer + (max_size * 3 / 4);
            calculate_strides(param, stride);
            buf_size = max_size;
        } else
            printf("video frame buffer is out of memory!!\n");
    }
    VideoFrame(const VideoFrame &origin) {
        // copy from origin
        *this = origin;
        // maintain own buffer
        const uint32_t luma_len = stride[0] * height;
        /** create video frame buffer with maximun size in 4 planes */
        uint32_t max_size = calculate_max_frame_size(origin);
        buffer = new uint8_t[max_size];
        if (buffer) {
            memset(buffer, 0, max_size);
            for (size_t i = 0; i < 4; ++i) {
                if (i != 3 || origin.planes[i]) {
                    const int buffer_len =
                        (i == 1 || i == 2) ? luma_len >> 2 : luma_len;
                    planes[i] = buffer + (i * luma_len);
                    memcpy(planes[i], origin.planes[i], buffer_len);
                }
            }
            buf_size = max_size;
        } else
            printf("video frame buffer is out of memory!!\n");
    }
    ~VideoFrame() {
        if (buf_size) {
            delete[] buffer;
            buffer = nullptr;
            memset(&planes, 0, sizeof(planes));
            buf_size = 0;
        }
    }
    /** Trim video frame buffer size for memory useage */
    void trim_buffer() {
        if (buf_size) {
            delete[] buffer;
            buffer = nullptr;
            buf_size = 0;
            memset(planes, 0, sizeof(planes));
        }
    }
    static uint32_t calculate_max_frame_size(const VideoFrameParam &param) {
        uint32_t luma_size = param.width * param.height;
        return 4 * luma_size;
    }
    static uint32_t calculate_max_frame_size(const VideoFrame &frame) {
        uint32_t luma_size = frame.stride[0] * frame.height;
        return 4 * luma_size;
    }
    static void calculate_strides(const VideoFrameParam &param,
                                  uint32_t strides[4]) {
        strides[0] = param.width;
        strides[3] = 0;
        switch (param.format) {
        case IMG_FMT_420:
        case IMG_FMT_420P10_PACKED:
        case IMG_FMT_NV12:
        case IMG_FMT_YV12_CUSTOM_COLOR_SPACE:
        case IMG_FMT_NV12_CUSTOM_COLOR_SPACE:
        case IMG_FMT_422:
        case IMG_FMT_422P10_PACKED:
            strides[1] = strides[2] = param.width >> 1;
            break;
        case IMG_FMT_444A: strides[3] = param.width;
        case IMG_FMT_444:
        case IMG_FMT_444P10_PACKED:
            strides[1] = strides[2] = param.width;
            break;
        default: assert(0); break;
        }
    }

} VideoFrame;

#endif  //_SVT_TEST_VIDEO_FRAME_H_
