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

#ifndef EbSvtAv1_h
#define EbSvtAv1_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include <stdint.h>
#include "EbSvtAv1Formats.h"
#include "EbDebugMacros.h"

struct SvtMetadataArray;

// API Version
#define SVT_AV1_VERSION_MAJOR 1
#define SVT_AV1_VERSION_MINOR 6
#define SVT_AV1_VERSION_PATCHLEVEL 0

#define SVT_AV1_CHECK_VERSION(major, minor, patch)                                                               \
    (SVT_AV1_VERSION_MAJOR > (major) || (SVT_AV1_VERSION_MAJOR == (major) && SVT_AV1_VERSION_MINOR > (minor)) || \
     (SVT_AV1_VERSION_MAJOR == (major) && SVT_AV1_VERSION_MINOR == (minor) && SVT_AV1_VERSION_PATCHLEVEL >= (patch)))

#if defined(_WIN32)
#define EB_HELPER_EXPORT __declspec(dllexport)
#define EB_HELPER_IMPORT __declspec(dllimport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#define EB_HELPER_EXPORT __attribute__((visibility("default")))
#define EB_HELPER_IMPORT
#else
#define EB_HELPER_EXPORT
#define EB_HELPER_IMPORT
#endif

#if defined(EB_DLL)
#if defined(EB_BUILDING_SHARED_LIBS)
#define EB_API EB_HELPER_EXPORT
#else
#define EB_API EB_HELPER_IMPORT
#endif // if defined(EB_BUILDING_SHARED_LIBS)
#else
#define EB_API
#endif //if defined(EB_DLL)

#define EB_MAX_NUM_OPERATING_POINTS 32

#define MAX_TEMPORAL_LAYERS 6

#define EB_MAX_TEMPORAL_LAYERS MAX_TEMPORAL_LAYERS

/********************************
* Defines
********************************/
#define EB_PICTURE uint32_t

typedef enum EbAv1PictureType {
    EB_AV1_INTER_PICTURE         = 0,
    EB_AV1_ALT_REF_PICTURE       = 1,
    EB_AV1_INTRA_ONLY_PICTURE    = 2,
    EB_AV1_KEY_PICTURE           = 3,
    EB_AV1_NON_REF_PICTURE       = 4,
    EB_AV1_SHOW_EXISTING_PICTURE = 6,
    EB_AV1_FW_KEY_PICTURE        = 5,
    EB_AV1_SWITCH_PICTURE        = 7,
    EB_AV1_INVALID_PICTURE       = 0xFF
} EbAv1PictureType;

/** The Bool type is intended to be used to represent a true or a false
value when passing parameters to and from the eBrisk API.  The
Bool is an 8 bit quantity.
*/
typedef uint8_t Bool;
#define FALSE 0
#define TRUE 1

typedef struct EbBufferHeaderType {
    // EbBufferHeaderType size
    uint32_t size;

    // picture (input or output) buffer
    uint8_t *p_buffer;
    uint32_t n_filled_len;
    uint32_t n_alloc_len;

    // pic private data
    void *p_app_private;
    void *wrapper_ptr;

    // pic timing param
    uint32_t n_tick_count;
    int64_t  dts;
    int64_t  pts;

    // pic info
    uint32_t         qp;
    EbAv1PictureType pic_type;
    uint64_t         luma_sse;
    uint64_t         cr_sse;
    uint64_t         cb_sse;
    // pic flags
    uint32_t flags;

    double luma_ssim;
    double cr_ssim;
    double cb_ssim;

    struct SvtMetadataArray *metadata;
} EbBufferHeaderType;

typedef struct EbComponentType {
    uint32_t size;
    void    *p_component_private;
    void    *p_application_private;
} EbComponentType;

typedef enum EbErrorType {
    EB_ErrorNone                   = 0,
    EB_DecUnsupportedBitstream     = (int32_t)0x40001000,
    EB_DecNoOutputPicture          = (int32_t)0x40001004,
    EB_DecDecodingError            = (int32_t)0x40001008,
    EB_Corrupt_Frame               = (int32_t)0x4000100C,
    EB_ErrorInsufficientResources  = (int32_t)0x80001000,
    EB_ErrorUndefined              = (int32_t)0x80001001,
    EB_ErrorInvalidComponent       = (int32_t)0x80001004,
    EB_ErrorBadParameter           = (int32_t)0x80001005,
    EB_ErrorDestroyThreadFailed    = (int32_t)0x80002012,
    EB_ErrorSemaphoreUnresponsive  = (int32_t)0x80002021,
    EB_ErrorDestroySemaphoreFailed = (int32_t)0x80002022,
    EB_ErrorCreateMutexFailed      = (int32_t)0x80002030,
    EB_ErrorMutexUnresponsive      = (int32_t)0x80002031,
    EB_ErrorDestroyMutexFailed     = (int32_t)0x80002032,
    EB_NoErrorEmptyQueue           = (int32_t)0x80002033,
    EB_NoErrorFifoShutdown         = (int32_t)0x80002034,
    EB_ErrorMax                    = 0x7FFFFFFF
} EbErrorType;

/* AV1 bistream profile (seq_profile syntax element) */
typedef enum EbAv1SeqProfile { MAIN_PROFILE = 0, HIGH_PROFILE = 1, PROFESSIONAL_PROFILE = 2 } EbAv1SeqProfile;

// For 8-bit and 10-bit packed inputs and outputs, the luma, cb, and cr fields should be used
//   for the three input picture planes.  However, for 10-bit unpacked planes the
//   lumaExt, cbExt, and crExt fields should be used hold the extra 2-bits of
//   precision while the luma, cb, and cr fields hold the 8-bit data.
typedef struct EbSvtIOFormat //former EbSvtEncInput
{
    // Hosts 8 bit or 16 bit input YUV420p / YUV420p10le
    uint8_t *luma;
    uint8_t *cb;
    uint8_t *cr;

    // Hosts LSB 2 bits of 10bit input/output when the compressed 10bit format is used
#if !SVT_AV1_CHECK_VERSION(1, 5, 0)
    /* DEPRECATED: to be removed in 1.5.0. */
    void *luma_ext;
    void *cb_ext;
    void *cr_ext;
#endif

    uint32_t y_stride;
    uint32_t cr_stride;
    uint32_t cb_stride;

    uint32_t width;
    uint32_t height;

    uint32_t org_x;
    uint32_t org_y;

    EbColorFormat color_fmt;
    EbBitDepth    bit_depth;
} EbSvtIOFormat;

typedef struct EbOperatingParametersInfo {
    /*!<Specifies the time interval between the arrival of the first bit in the
     * smoothing buffer and the subsequent removal of the data that belongs to
     * the first coded frame for operating point*/
    uint32_t decoder_buffer_delay;

    /*!<Specifies, in combination with decoder_buffer_delay[op] syntax element,
     * the first bit arrival time of frames to be decoded to the smoothing
     * buffer */
    uint32_t encoder_buffer_delay;

    /*!< Equal to 1 indicates that the smoothing buffer operates in low-delay
     * mode for operating point*/
    uint8_t low_delay_mode_flag;

} EbOperatingParametersInfo;

typedef struct EbAV1OperatingPoint {
    uint32_t op_idc;
    uint32_t seq_level_idx;
    uint32_t seq_tier;

    /*!< 1 -> Indicates that there is a decoder model associated with operating
             point,
     *   0 -> Indicates that there is not a decoder model associated with
             operating point*/
    uint8_t decoder_model_present_for_this_op;

    /*!< Operating Parameters Information structure*/
    EbOperatingParametersInfo operating_parameters_info;

    uint32_t initial_display_delay_present_for_this_op;
    uint32_t initial_display_delay;

} EbAv1OperatingPoint;

typedef struct EbColorConfig {
    /*!< bit depth */
    EbBitDepth bit_depth;

    /*!< 1: Indicates that the video does not contain U and V color planes.
     *   0: Indicates that the video contains Y, U, and V color planes. */
    Bool mono_chrome;

    /*!< Specify the chroma subsampling format */
    uint8_t subsampling_x;

    /*!< Specify the chroma subsampling format */
    uint8_t subsampling_y;

    /*!< 1: Specifies that color_primaries, transfer_characteristics, and
            matrix_coefficients are present. color_description_present_flag
     *   0: Specifies that color_primaries, transfer_characteristics and
            matrix_coefficients are not present */
    Bool color_description_present_flag;

    /*!< An integer that is defined by the "Color primaries" section of
     * ISO/IEC 23091-4/ITU-T H.273 */
    EbColorPrimaries color_primaries;

    /*!< An integer that is defined by the "Transfer characteristics" section
     * of ISO/IEC 23091-4/ITU-T H.273 */
    EbTransferCharacteristics transfer_characteristics;

    /*!< An integer that is defined by the "Matrix coefficients" section of
     * ISO/IEC 23091-4/ITU-T H.273 */
    EbMatrixCoefficients matrix_coefficients;

    /*!< 0: shall be referred to as the studio swing representation
     *   1: shall be referred to as the full swing representation */
    EbColorRange color_range;

    /*!< Specifies the sample position for subsampled streams */
    EbChromaSamplePosition chroma_sample_position;

    /*!< 1: Indicates that the U and V planes may have separate delta quantizer
     *   0: Indicates that the U and V planes will share the same delta
            quantizer value */
    Bool separate_uv_delta_q;

} EbColorConfig;

typedef struct EbTimingInfo {
    /*!< Timing info present flag */
    Bool timing_info_present;

    /*!< Number of time units of a clock operating at the frequency time_scale
     * Hz that corresponds to one increment of a clock tick counter*/
    uint32_t num_units_in_display_tick;

    /*!< Number of time units that pass in one second*/
    uint32_t time_scale;

    /*!< Equal to 1 indicates that pictures should be displayed according to
     * their output order with the number of ticks between two consecutive
     * pictures specified by num_ticks_per_picture.*/
    uint8_t equal_picture_interval;

    /*!< Specifies the number of clock ticks corresponding to output time
     * between two consecutive pictures in the output order.
     * Range - [0 to (1 << 32) - 2]*/
    uint32_t num_ticks_per_picture;

} EbTimingInfo;

// structure to be allocated at the sample application and passed to the library
// on a per picture basis through the p_app_private field in the EbBufferHeaderType structure
// this structure and the data inside would be casted, validated, then copied at the
// svt_av1_enc_send_picture API call
typedef enum {
    PRIVATE_DATA, // data to be passed through and written to the bitstream
    //FILM_GRAIN_PARAM,        // passing film grain parameters per picture
    REF_FRAME_SCALING_EVENT, // reference frame scaling data per picture
    ROI_MAP_EVENT, // ROI map data per picture
    PRIVATE_DATA_TYPES // end of private data types
} PrivDataType;
typedef struct EbPrivDataNode {
    PrivDataType           node_type;
    void                  *data; // pointer to data structure e.g. EbRefFrameScale or AomFilmGrain
    uint32_t               size; // size of data being sent for the library to know how much to copy
    struct EbPrivDataNode *next; // pointer to the next node, NULL if done.
} EbPrivDataNode;
typedef struct EbRefFrameScale {
    uint8_t  scale_mode; // scaling mode, support for RESIZE_NONE, RESIZE_FIXED and RESIZE_RANDOM
    uint32_t scale_denom; // scaling denominator for non-key frame, from 8~16
    uint32_t scale_kf_denom; // scaling denominator for key frame, from 8~16
} EbRefFrameScale;
typedef struct SvtAv1RoiMapEvt {
    uint64_t                start_picture_number;
    uint8_t                *b64_seg_map;
    int16_t                 seg_qp[8]; // 8: MAX_SEGMENTS
    int8_t                  max_seg_id;
    struct SvtAv1RoiMapEvt *next;
} SvtAv1RoiMapEvt;
typedef struct SvtAv1RoiMap {
    uint32_t         evt_num;
    SvtAv1RoiMapEvt *evt_list;
    SvtAv1RoiMapEvt *cur_evt;
    int16_t         *qp_map;
    char            *buf;
} SvtAv1RoiMap;
/**
CPU FLAGS
*/
typedef uint64_t EbCpuFlags;
#define EB_CPU_FLAGS_MMX (1 << 0)
#define EB_CPU_FLAGS_SSE (1 << 1)
#define EB_CPU_FLAGS_SSE2 (1 << 2)
#define EB_CPU_FLAGS_SSE3 (1 << 3)
#define EB_CPU_FLAGS_SSSE3 (1 << 4)
#define EB_CPU_FLAGS_SSE4_1 (1 << 5)
#define EB_CPU_FLAGS_SSE4_2 (1 << 6)
#define EB_CPU_FLAGS_AVX (1 << 7)
#define EB_CPU_FLAGS_AVX2 (1 << 8)
#define EB_CPU_FLAGS_AVX512F (1 << 9)
#define EB_CPU_FLAGS_AVX512CD (1 << 10)
#define EB_CPU_FLAGS_AVX512DQ (1 << 11)
#define EB_CPU_FLAGS_AVX512ER (1 << 12)
#define EB_CPU_FLAGS_AVX512PF (1 << 13)
#define EB_CPU_FLAGS_AVX512BW (1 << 14)
#define EB_CPU_FLAGS_AVX512VL (1 << 15)
#define EB_CPU_FLAGS_ALL ((EB_CPU_FLAGS_AVX512VL << 1) - 1)
#define EB_CPU_FLAGS_INVALID (1ULL << (sizeof(EbCpuFlags) * 8ULL - 1ULL))

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbSvtAv1_h
