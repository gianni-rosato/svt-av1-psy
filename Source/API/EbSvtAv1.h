/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbSvtAv1_h
#define EbSvtAv1_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "stdint.h"
#include "EbSvtAv1Formats.h"

// API Version
#define SVT_VERSION_MAJOR 0
#define SVT_VERSION_MINOR 8
#define SVT_VERSION_PATCHLEVEL 2

#ifdef _WIN32
#define EB_API __declspec(dllexport)
#else
#define EB_API
#endif

#define EB_MAX_NUM_OPERATING_POINTS 32

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

/** The EbBool type is intended to be used to represent a true or a false
value when passing parameters to and from the eBrisk API.  The
EbBool is a 32 bit quantity and is aligned on a 32 bit word boundary.
*/

#define EbBool uint8_t
#define EB_FALSE 0
#define EB_TRUE 1

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
    uint32_t qp;
    uint32_t pic_type;
    uint32_t luma_sse;
    uint32_t cr_sse;
    uint32_t cb_sse;

    // pic flags
    uint32_t flags;
} EbBufferHeaderType;

typedef struct EbComponentType {
    uint32_t size;
    void *   p_component_private;
    void *   p_application_private;
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
typedef enum EbAv1SeqProfile {
    MAIN_PROFILE         = 0,
    HIGH_PROFILE         = 1,
    PROFESSIONAL_PROFILE = 2
} EbAv1SeqProfile;

typedef enum AomBitDepth {
    AOM_BITS_8  = 8, /**<  8 bits */
    AOM_BITS_10 = 10, /**< 10 bits */
    AOM_BITS_12 = 12, /**< 12 bits */
} AomBitDepth;

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
    uint8_t *luma_ext;
    uint8_t *cb_ext;
    uint8_t *cr_ext;

    uint32_t y_stride;
    uint32_t cr_stride;
    uint32_t cb_stride;

    uint32_t width;
    uint32_t height;

    uint32_t origin_x;
    uint32_t origin_y;

    EbColorFormat color_fmt;
    EbBitDepth    bit_depth;
} EbSvtIOFormat;

typedef struct BitstreamLevel {
    uint8_t major;
    uint8_t minor;
} BitstreamLevel;

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
    uint32_t bit_depth;

    /*!< 1: Indicates that the video does not contain U and V color planes.
     *   0: Indicates that the video contains Y, U, and V color planes. */
    EbBool mono_chrome;

    /*!< Specify the chroma subsampling format */
    uint8_t subsampling_x;

    /*!< Specify the chroma subsampling format */
    uint8_t subsampling_y;

    /*!< 1: Specifies that color_primaries, transfer_characteristics, and
            matrix_coefficients are present. color_description_present_flag
     *   0: Specifies that color_primaries, transfer_characteristics and
            matrix_coefficients are not present */
    EbBool color_description_present_flag;

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
    EbBool separate_uv_delta_q;

} EbColorConfig;

typedef struct EbTimingInfo {
    /*!< Timing info present flag */
    EbBool timing_info_present;

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

/**
CPU FLAGS
*/
typedef uint64_t CPU_FLAGS;
#define CPU_FLAGS_MMX (1 << 0)
#define CPU_FLAGS_SSE (1 << 1)
#define CPU_FLAGS_SSE2 (1 << 2)
#define CPU_FLAGS_SSE3 (1 << 3)
#define CPU_FLAGS_SSSE3 (1 << 4)
#define CPU_FLAGS_SSE4_1 (1 << 5)
#define CPU_FLAGS_SSE4_2 (1 << 6)
#define CPU_FLAGS_AVX (1 << 7)
#define CPU_FLAGS_AVX2 (1 << 8)
#define CPU_FLAGS_AVX512F (1 << 9)
#define CPU_FLAGS_AVX512CD (1 << 10)
#define CPU_FLAGS_AVX512DQ (1 << 11)
#define CPU_FLAGS_AVX512ER (1 << 12)
#define CPU_FLAGS_AVX512PF (1 << 13)
#define CPU_FLAGS_AVX512BW (1 << 14)
#define CPU_FLAGS_AVX512VL (1 << 15)
#define CPU_FLAGS_ALL ((CPU_FLAGS_AVX512VL << 1) - 1)
#define CPU_FLAGS_INVALID (1ULL << (sizeof(CPU_FLAGS) * 8ULL - 1ULL))

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbSvtAv1_h
