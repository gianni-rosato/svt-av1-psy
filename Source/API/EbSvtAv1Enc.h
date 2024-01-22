﻿/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbSvtAv1Enc_h
#define EbSvtAv1Enc_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include <stdint.h>
#include "EbSvtAv1.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief SVT-AV1 encoder ABI version
 *
 * Should be increased by 1 everytime a public struct in the encoder
 * has been modified, and reset anytime the major API version has
 * been changed. Used to keep track if a field has been added or not.
 */
#define SVT_AV1_ENC_ABI_VERSION 9

//***HME***

#define MAX_HIERARCHICAL_LEVEL 6
#define REF_LIST_MAX_DEPTH 4
/*!\brief Decorator indicating that given struct/union/enum is packed */
#ifndef ATTRIBUTE_PACKED
#if defined(__GNUC__) && __GNUC__
#define ATTRIBUTE_PACKED __attribute__((packed))
#else
#define ATTRIBUTE_PACKED
#endif
#endif /* ATTRIBUTE_PACKED */
typedef enum ATTRIBUTE_PACKED {
    ENC_MRS        = -2, // Highest quality research mode (slowest)
    ENC_MR         = -1, //Research mode with higher quality than M0
    ENC_M0         = 0,
    ENC_M1         = 1,
    ENC_M2         = 2,
    ENC_M3         = 3,
    ENC_M4         = 4,
    ENC_M5         = 5,
    ENC_M6         = 6,
    ENC_M7         = 7,
    ENC_M8         = 8,
    ENC_M9         = 9,
    ENC_M10        = 10,
    ENC_M11        = 11,
    ENC_M12        = 12,
    ENC_M13        = 13,
    MAX_ENC_PRESET = ENC_M13
} EncMode;

#define DEFAULT -1

#define EB_BUFFERFLAG_EOS 0x00000001 // signals the last packet of the stream
#define EB_BUFFERFLAG_SHOW_EXT 0x00000002 // signals that the packet contains a show existing frame at the end
#define EB_BUFFERFLAG_HAS_TD 0x00000004 // signals that the packet contains a TD
#define EB_BUFFERFLAG_IS_ALT_REF 0x00000008 // signals that the packet contains an ALT_REF frame
#define EB_BUFFERFLAG_ERROR_MASK \
    0xFFFFFFF0 // mask for signalling error assuming top flags fit in 4 bits. To be changed, if more flags are added.

/*
 * Struct for storing content light level information
 * Values are stored in BE format
 * Refer to the AV1 specification 6.7.3 for more details
 */
struct EbContentLightLevel {
    uint16_t max_cll;
    uint16_t max_fall;
};

/*
 * Struct for storing x and y chroma points, values are stored in BE format
 */
struct EbSvtAv1ChromaPoints {
    uint16_t x;
    uint16_t y;
};

/*
 * Struct for storing mastering-display information
 * values are stored in BE format
 * Refer to the AV1 specification 6.7.4 for more details
 */
struct EbSvtAv1MasteringDisplayInfo {
    struct EbSvtAv1ChromaPoints r;
    struct EbSvtAv1ChromaPoints g;
    struct EbSvtAv1ChromaPoints b;
    struct EbSvtAv1ChromaPoints white_point;
    uint32_t                    max_luma;
    uint32_t                    min_luma;
};

/************************************************
 * Prediction Structure Config Entry
 *   Contains the basic reference lists and
 *   configurations for each Prediction Structure
 *   Config Entry.
 ************************************************/
typedef struct PredictionStructureConfigEntry {
    uint32_t temporal_layer_index;
    uint32_t decode_order;
} PredictionStructureConfigEntry;

// super-res modes
typedef enum {
    SUPERRES_NONE, // No frame superres allowed.
    SUPERRES_FIXED, // All frames are coded at the specified scale, and super-resolved.
    SUPERRES_RANDOM, // All frames are coded at a random scale, and super-resolved.
    SUPERRES_QTHRESH, // Superres scale for a frame is determined based on q_index.
    SUPERRES_AUTO, // Automatically select superres for appropriate frames.
    SUPERRES_MODES
} SUPERRES_MODE;

// super-res auto search type
typedef enum {
    SUPERRES_AUTO_ALL, // Tries all possible superres ratios
    SUPERRES_AUTO_DUAL, // Tries no superres and q-based superres ratios
    SUPERRES_AUTO_SOLO, // Only apply the q-based superres ratio
    SUPERRES_AUTO_SEARCH_TYPES
} SUPERRES_AUTO_SEARCH_TYPE;

// reference scaling modes
typedef enum {
    RESIZE_NONE, // No frame resize allowed.
    RESIZE_FIXED, // All frames are coded at the specified scale.
    RESIZE_RANDOM, // All frames are coded at a random scale.
    RESIZE_DYNAMIC, // Resize scale for a frame in dynamic.
    RESIZE_RANDOM_ACCESS, // Resize scale frame by event in random access
    RESIZE_MODES
} RESIZE_MODE;

/** The SvtAv1IntraRefreshType is used to describe the intra refresh type.
*/
typedef enum SvtAv1IntraRefreshType {
    SVT_AV1_FWDKF_REFRESH = 1,
    SVT_AV1_KF_REFRESH    = 2,
} SvtAv1IntraRefreshType;

typedef enum {
    SVT_AV1_STREAM_INFO_START                = 1,
    SVT_AV1_STREAM_INFO_FIRST_PASS_STATS_OUT = SVT_AV1_STREAM_INFO_START,

    SVT_AV1_STREAM_INFO_END,
} SVT_AV1_STREAM_INFO_ID;

/*!\brief Generic fixed size buffer structure
 *
 * This structure is able to hold a reference to any fixed size buffer.
 */
typedef struct SvtAv1FixedBuf {
    void    *buf; /**< Pointer to the data. Does NOT own the data! */
    uint64_t sz; /**< Length of the buffer, in chars */
} SvtAv1FixedBuf; /**< alias for struct aom_fixed_buf */

/** Indicates how an S-Frame should be inserted.
*/
typedef enum EbSFrameMode {
    SFRAME_STRICT_BASE =
        1, /**< The considered frame will be made into an S-Frame only if it is a base layer inter frame */
    SFRAME_NEAREST_BASE =
        2, /**< If the considered frame is not an altref frame, the next base layer inter frame will be made into an S-Frame */
} EbSFrameMode;

/* Indicates what prediction structure to use
 * was PredStructure in EbDefinitions.h
 * Only SVT_AV1_PRED_LOW_DELAY_B and SVT_AV1_PRED_RANDOM_ACCESS are valid
 */
typedef enum SvtAv1PredStructure {
    SVT_AV1_PRED_LOW_DELAY_P   = 0, // No longer active
    SVT_AV1_PRED_LOW_DELAY_B   = 1,
    SVT_AV1_PRED_RANDOM_ACCESS = 2,
    SVT_AV1_PRED_TOTAL_COUNT   = 3,
    SVT_AV1_PRED_INVALID       = 0xFF,
} SvtAv1PredStructure;

/* Indicates what rate control mode is used.
 * Currently, cqp is distinguised by setting enable_adaptive_quantization to 0
 */
typedef enum SvtAv1RcMode {
    SVT_AV1_RC_MODE_CQP_OR_CRF = 0, // constant quantization parameter/constant rate factor
    SVT_AV1_RC_MODE_VBR        = 1, // variable bit rate
    SVT_AV1_RC_MODE_CBR        = 2, // constant bit rate
} SvtAv1RcMode;

typedef enum SvtAv1FrameUpdateType {
    SVT_AV1_KF_UPDATE,
    SVT_AV1_LF_UPDATE,
    SVT_AV1_GF_UPDATE,
    SVT_AV1_ARF_UPDATE,
    SVT_AV1_OVERLAY_UPDATE,
    SVT_AV1_INTNL_OVERLAY_UPDATE, // Internal Overlay Frame
    SVT_AV1_INTNL_ARF_UPDATE, // Internal Altref Frame
    SVT_AV1_FRAME_UPDATE_TYPES
} SvtAv1FrameUpdateType;

typedef struct SvtAv1FrameScaleEvts {
    uint32_t  evt_num;
    uint64_t *start_frame_nums;
    uint32_t *resize_kf_denoms;
    uint32_t *resize_denoms;
} SvtAv1FrameScaleEvts;

// Will contain the EbEncApi which will live in the EncHandle class
// Only modifiable during config-time.
typedef struct EbSvtAv1EncConfiguration {
    /**
     * @brief Encoder preset used.
     * -2 and -1 are for debug purposes and should not be used.
     * 0 is the highest quality mode but is the slowest,
     * 13 is the fastest mode but is not as high quality.
     *
     * Min value is -2.
     * Max value is 13.
     * Default is 12.
     */
    int8_t enc_mode;

    // GOP Structure

    /* The intra period defines the interval of frames after which you insert an
     * Intra refresh. It is strongly recommended to set the value to multiple of
     * 2^(hierarchical_levels), subtracting one if using open GOP (intra_refresh_type == 1).
     * For instance, to get a 5-second GOP (default being >=5 seconds)
     * with hierarchical_levels = 3 and open GOP you could use 319, 279, 159
     * for 60, 50, or 30 respectively.
     *
     * -1 = no intra update.
     * -2 = auto.
     *
     * Default is -2. */
    int32_t intra_period_length;

    /* Random access.
     *
     * 1 = CRA, open GOP.
     * 2 = IDR, closed GOP.
     *
     * Default is 1. */
    SvtAv1IntraRefreshType intra_refresh_type;

    /* Number of hierarchical layers used to construct GOP.
     * Minigop size = 2^HierarchicalLevels.
     *
     * Default is 5 upt to M12 4, for M13. */
    uint32_t hierarchical_levels;

    /* Prediction structure used to construct GOP. There are two main structures
     * supported, which are: Low Delay (P or B) and Random Access.
     *
     * In Low Delay structure, pictures within a mini GOP refer to the previously
     * encoded pictures in display order. In other words, pictures with display
     * order N can only be referenced by pictures with display order greater than
     * N, and it can only refer pictures with picture order lower than N. The Low
     * Delay structure can be flat structured (e.g. IPPPPPPP...) or hierarchically
     * structured. B/b pictures can be used instead of P/p pictures. However, the
     * reference picture list 0 and the reference picture list 1 will contain the
     * same reference picture.
     *
     * In Random Access structure, the B/b pictures can refer to reference pictures
     * from both directions (past and future).
     *
     * Refer to SvtAv1PredStructure enum for valid values.
     *
     * Default is SVT_AV1_PRED_RANDOM_ACCESS. */
    uint8_t pred_structure;

    // Input Info

    /**
     * @brief Frame width in pixels.
     *
     * Min is 64.
     * Max is 16384.
     * Default is 0.
     */
    uint32_t source_width;

    /**
     * @brief Frame height in pixels
     *
     * Min is 64.
     * Max is 8704.
     * Default is 0.
     */
    uint32_t source_height;

    /* Specifies the maximum frame width/height for the frames represented by the sequence header
     * (max_frame_width_minus_1 and max_frame_height_minus_1, spec 5.5.1).
     * Actual frame height could be equal to or less than this value. E.g. Use this value to indicate
     * the maximum height between renditions when switch frame feature is on.
     */
    uint32_t forced_max_frame_width;
    uint32_t forced_max_frame_height;
    /* Frame rate numerator. When zero, the encoder will use -fps if
     * FrameRateDenominator is also zero, otherwise an error is returned.
     *
     * Default is 0. */
    uint32_t frame_rate_numerator;
    /* Frame rate denominator. When zero, the encoder will use -fps if
     * FrameRateNumerator is also zero, otherwise an error is returned.
     *
     * Default is 0. */
    uint32_t frame_rate_denominator;
    /* Specifies the bit depth of input video.
     *
     * 8 = 8 bit.
     * 10 = 10 bit.
     *
     * Default is 8. */
    uint32_t encoder_bit_depth;

    /**
     * @brief Encoder color format.
     * Only YUV420 is supported for now.
     *
     * Min is YUV400.
     * Max is YUV444.
     * Default is YUV420.
     */
    EbColorFormat encoder_color_format;

    /**
     * @brief Currently unused.
     *
     * Default is 0.
     */
    uint8_t high_dynamic_range_input;

    /**
     * @brief Bitstream profile to use.
     * 0: main, 1: high, 2: professional.
     *
     * Min is MAIN_PROFILE.
     * Max is PROFESSIONAL_PROFILE.
     * Default is MAIN_PROFILE.
     */
    EbAv1SeqProfile profile;
    /* Constraints for bitstream in terms of max bitrate and max buffer size.
     *
     * 0 = Main, for most applications.
     * 1 = High, for demanding applications.
     *
     * Default is 0. */
    uint32_t tier;

    /**
     * @brief Bitstream level.
     * 0: autodetect from bitstream, 20: level 2.0, 63: level 6.3, only levels 2.0-6.3 are properly defined.
     * The levels are defined at https://aomediacodec.github.io/av1-spec/av1-spec.pdf
     * under "A.3. Levels".
     *
     * Min is 0.
     * Max is 73.
     * Default is 0.
     */
    uint32_t level;

    /* Color description present flag
    *
    * It is not necessary to set this parameter manually.
    * It is set internally to true once one of the color_primaries, transfer_characteristics or
    * matrix coefficients is set to non-default value.
    *
    Default is false. */
    Bool color_description_present_flag;
    /* Color primaries
    * values are from EbColorPrimaries
    Default is 2 (CP_UNSPECIFIED). */
    EbColorPrimaries color_primaries;
    /* Transfer characteristics
    * values are from EbTransferCharacteristics
    Default is 2 (TC_UNSPECIFIED). */
    EbTransferCharacteristics transfer_characteristics;
    /* Matrix coefficients
    * values are from EbMatrixCoefficients
    Default is 2 (MC_UNSPECIFIED). */
    EbMatrixCoefficients matrix_coefficients;
    /* Color range
    * values are from EbColorRange
    * 0: studio swing.
    * 1: full swing.
    Default is 0. */
    EbColorRange color_range;
    /* Mastering display metadata
    * values are from set using svt_aom_parse_mastering_display()
    */
    struct EbSvtAv1MasteringDisplayInfo mastering_display;
    /* Content light level
    * values are from set using svt_aom_parse_content_light_level()
    */
    struct EbContentLightLevel content_light_level;

    /* Chroma sample position
     * Values as per 6.4.2 of the specification:
     * EB_CSP_UNKNOWN:   default
     * EB_CSP_VERTICAL:  value 0 from H.273 AKA "left"
     * EB_CSP_COLOCATED: value 2 from H.273 AKA "top left"
     */
    EbChromaSamplePosition chroma_sample_position;

    // End input info

    /* Rate control mode.
     *
     * Refer to the SvtAv1RcMode enum for valid values
     * Default is 0. */
    uint32_t rate_control_mode;

    // Rate control tuning

    // Quantization
    /* Initial quantization parameter for the Intra pictures used under constant
     * qp rate control mode.
     *
     * Default is 50. */
    uint32_t qp;

    /* force qp values for every picture that are passed in the header pointer
    *
    * Default is 0.*/
    Bool use_qp_file;

    /* Target bitrate in bits/second, only applicable when rate control mode is
     * set to 1 (VBR) or 2 (CBR).
     *
     * Default is 2000000. */
    uint32_t target_bit_rate;
    /* maximum bitrate in bits/second, only apllicable when rate control mode is
     * set to 0.
     *
     * Default is 0. */
    uint32_t max_bit_rate;

#if !SVT_AV1_CHECK_VERSION(1, 5, 0)
    /* DEPRECATED: to be removed in 1.5.0. */
    uint32_t vbv_bufsize;
#endif

    /* Maxium QP value allowed for rate control use, only applicable when rate
     * control mode is set to 1. It has to be greater or equal to minQpAllowed.
     *
     * Default is 63. */
    uint32_t max_qp_allowed;
    /* Minimum QP value allowed for rate control use, only applicable when rate
     * control mode is set to 1. It has to be smaller or equal to maxQpAllowed.
     *
     * Default is 0. */
    uint32_t min_qp_allowed;

    // DATARATE CONTROL OPTIONS
#if !SVT_AV1_CHECK_VERSION(2, 0, 0)
    /* DEPRECATED: to be removed in 2.0.0. */
    /**
     * @brief Variable Bit Rate Bias Percentage
     *
     * Indicates the bias for determining target size for the current frame.
     * A value 0 indicates the optimal CBR mode value should be used, and 100
     * indicates the optimal VBR mode value should be used.
     *
     * Min is 0.
     * Max is 100.
     * Default is 100
     */
    uint32_t vbr_bias_pct;
#endif

    /**
     * @brief Variable Bit Rate Minimum Section Percentage
     *
     * Indicates the minimum bitrate to be used for a single GOP as a percentage
     * of the target bitrate.
     *
     * Min is 0.
     * Max is 100.
     * Default is 0.
     */
    uint32_t vbr_min_section_pct;

    /**
     * @brief Variable Bit Rate Maximum Section Percentage
     *
     * Indicates the maximum bitrate to be used for a single GOP as a percentage
     * of the target bitrate.
     *
     * Min is 0.
     * Max is 10000.
     * Default is 2000.
     */
    uint32_t vbr_max_section_pct;

    /**
     * @brief UnderShoot Percentage
     *
     * Only applicable for VBR and CBR.
     *
     * Indicates the tolerance of the VBR algorithm to undershoot and is used
     * as a trigger threshold for more agressive adaptation of Quantization.
     *
     * Min is 0.
     * Max is 100.
     * Default is 25 for CBR and 50 for VBR.
     */
    uint32_t under_shoot_pct;

    /**
     * @brief OverShoot Percentage
     *
     * Only applicable for VBR and CBR
     *
     * Indicates the tolerance of the VBR algorithm to overshoot and is used as
     * a trigger threshold for more agressive adaptation of Quantization.
     *
     * Min is 0.
     * Max is 100.
     * Default is 25.
     */
    uint32_t over_shoot_pct;

    /**
     * @brief MaxBitRate OverShoot Percentage
     *
     * Only applicable for Capped CRF.
     *
     * Indicates the tolerance of the Capped CRF algorithm to overshoot
     * and is used as a trigger threshold for more agressive adaptation of
     * Quantization.
     *
     * Min is 0.
     * Max is 100.
     * Default is 50.
     */
    uint32_t mbr_over_shoot_pct;

    /**
     * @brief Starting Buffer Level in MilliSeconds
     *
     * Only applicable for CBR.
     *
     * Indicates the amount of data that will be buffered by the decoding
     * application prior to beginning playback, and is expressed in units of
     * time. Must be less than maximum_buffer_size_ms.
     *
     * Min is 20.
     * Max is 10000.
     * Default is 600.
     */
    int64_t starting_buffer_level_ms;

    /**
     * @brief Optimal Buffer Level in MilliSeconds
     *
     * Only applicable for CBR.
     *
     * indicates the amount of data that the encoder should try to maintain in the
     * decoder's buffer, and is expressed in units of time. Must be less than
     * maximum_buffer_size_ms.
     *
     * Min is 20.
     * Max is 10000.
     * Default is 600.
     */
    int64_t optimal_buffer_level_ms;

    /**
     * @brief Maximum Buffer Size in MilliSeconds
     *
     * Only applicable for CBR.
     *
     * indicates the maximum amount of data that may be buffered by the
     * decoding application, and is expressed in units of time.
     *
     * Min is 20.
     * Max is 10000.
     * Default is 1000.
     */
    int64_t maximum_buffer_size_ms;

    // input / output buffer to be used for multi-pass encoding
    SvtAv1FixedBuf rc_stats_buffer;
    int            pass;

    // End rate control tuning

    // Individual tuning flags
    /* use fixed qp offset for every picture based on temporal layer index
    * 0: off (use the auto mode QP)
    * 1: on (the offset is applied on top of the user QP)
    * 2: on (the offset is applied on top of the auto mode QP)
    *
    * Default is 0.*/
    uint8_t use_fixed_qindex_offsets;
    int32_t qindex_offsets[EB_MAX_TEMPORAL_LAYERS];
    int32_t key_frame_chroma_qindex_offset;
    int32_t key_frame_qindex_offset;
    int32_t chroma_qindex_offsets[EB_MAX_TEMPORAL_LAYERS];

    int32_t luma_y_dc_qindex_offset;
    int32_t chroma_u_dc_qindex_offset;
    int32_t chroma_u_ac_qindex_offset;
    int32_t chroma_v_dc_qindex_offset;
    int32_t chroma_v_ac_qindex_offset;

    /**
     * @brief Deblocking loop filter control
     *
     * Default is true.
     */
    Bool enable_dlf_flag;

    /* Film grain denoising the input picture
    * Flag to enable the denoising
    *
    * Default is 0. */
    uint32_t film_grain_denoise_strength;

    /**
    * @brief Determines how much denoising is used.
    * Only applicable when film grain is ON.
    *
    * 0 is no denoising
    * 1 is full denoising
    */
    uint8_t film_grain_denoise_apply;

    /* CDEF Level
    *
    * Default is -1. */
    int cdef_level;

    /* Restoration filtering
    *  enable/disable
    *  set Self-Guided (sg) mode
    *  set Wiener (wn) mode
    *
    * Default is -1. */
    int enable_restoration_filtering;

    /* motion field motion vector
    *
    *  Default is -1. */
    int enable_mfmv;

    /* Flag to enable the scene change detection algorithm.
     *
     * Default is 1. */
    uint32_t scene_change_detection;

    /**
     * @brief API signal to constrain motion vectors.
     *
     * Default is false.
     */
    Bool restricted_motion_vector;

    /* Log 2 Tile Rows and columns . 0 means no tiling,1 means that we split the dimension
        * into 2
        * Default is 0. */
    int32_t tile_columns;
    int32_t tile_rows;

    /* When RateControlMode is set to 1 it's best to set this parameter to be
     * equal to the Intra period value (such is the default set by the encoder).
     * When CQP is chosen, then a (2 * minigopsize +1) look ahead is recommended.
     *
     * Default depends on rate control mode.*/
    uint32_t look_ahead_distance;

    /* Enable TPL in look ahead
     * 0 = disable TPL in look ahead
     * 1 = enable TPL in look ahead
     * Default is 0  */
    uint8_t enable_tpl_la;

    /* recode_loop indicates the recode levels,
     * DISALLOW_RECODE = 0, No recode.
     * ALLOW_RECODE_KFMAXBW = 1, Allow recode for KF and exceeding maximum frame bandwidth.
     * ALLOW_RECODE_KFARFGF = 2, Allow recode only for KF/ARF/GF frames.
     * ALLOW_RECODE = 3, Allow recode for all frames based on bitrate constraints.
     * ALLOW_RECODE_DEFAULT = 4, Default setting, ALLOW_RECODE_KFARFGF for M0~5 and
     *                                            ALLOW_RECODE_KFMAXBW for M6~8.
     * default is 4
     */
    uint32_t recode_loop;

    /* Flag to signal the content being a screen sharing content type
    *
    * Default is 0. */
    uint32_t screen_content_mode;

    /* Enable adaptive quantization within a frame using segmentation.
     *
     * For rate control mode 0, setting this to 0 will use CQP mode, else CRF mode will be used.
     * Default is 2. */
    uint8_t enable_adaptive_quantization;

    /**
     * @brief Enable use of ALT-REF (temporally filtered) frames.
     *
     * Default is true.
     */
    Bool enable_tf;

    Bool enable_overlays;
    /**
     * @brief Tune for a particular metric; 0: VQ, 1: PSNR, 2: SSIM.
     *
     * Default is 1.
     */
    uint8_t tune;

    // super-resolution parameters
    uint8_t superres_mode;
    uint8_t superres_denom;
    uint8_t superres_kf_denom;
    uint8_t superres_qthres;
    uint8_t superres_kf_qthres;
    uint8_t superres_auto_search_type;

#if !SVT_AV1_CHECK_VERSION(1, 5, 0)
    /* DEPRECATED: to be removed in 1.5.0. */
    PredictionStructureConfigEntry pred_struct[1 << (MAX_HIERARCHICAL_LEVEL - 1)];

    /* DEPRECATED: to be removed in 1.5.0. */
    Bool enable_manual_pred_struct;

    /* DEPRECATED: to be removed in 1.5.0. */
    int32_t manual_pred_struct_entry_num;
#endif
    /* Decoder-speed-targeted encoder optimization level (produce bitstreams that can be decoded faster).
    * 0: No decoder speed optimization
    * 1: Decoder speed optimization enabled (fast decode)
    */
    Bool fast_decode;

    /* S-Frame interval (frames)
    * 0: S-Frame off
    * >0: S-Frame on and indicates the number of frames after which a frame may be coded as an S-Frame
    */
    int32_t sframe_dist;
    /* Indicates how an S-Frame should be inserted
    * values are from EbSFrameMode
    * SFRAME_STRICT_ARF: the considered frame will be made into an S-Frame only if it is an altref frame
    * SFRAME_NEAREST_ARF: if the considered frame is not an altref frame, the next altref frame will be made into an S-Frame
    */
    EbSFrameMode sframe_mode;

    // End of individual tuning flags

    // Application Specific parameters

    /**
     * @brief API signal for the library to know the channel ID (used for pinning to cores).
     *
     * Min value is 0.
     * Max value is 0xFFFFFFFF.
     * Default is 0.
     */
    uint32_t channel_id;

    /**
     * @brief API signal for the library to know the active number of channels being encoded simultaneously.
     *
     * Min value is 1.
     * Max value is 0xFFFFFFFF.
     * Default is 1.
     */
    uint32_t active_channel_count;

    // Threads management

    /* The number of logical processor which encoder threads run on. If
     * LogicalProcessors and TargetSocket are not set, threads are managed by
     * OS thread scheduler. */
    uint32_t logical_processors;

    /* Unpin the execution .This option does not
    * set the execution to be pinned to a specific number of cores when set to 1. this allows the execution
    * of multiple encodes on the CPU without having to pin them to a specific mask
    * 1: pinned threads
    * 0: unpinned
    * default 0 */
    uint32_t pin_threads;

    /* Target socket to run on. For dual socket systems, this can specify which
     * socket the encoder runs on.
     *
     * -1 = Both Sockets.
     *  0 = Socket 0.
     *  1 = Socket 1.
     *
     * Default is -1. */
    int32_t target_socket;

    /* CPU FLAGS to limit assembly instruction set used by encoder.
    * Default is EB_CPU_FLAGS_ALL. */
    EbCpuFlags use_cpu_flags;

    // Debug tools
    /* Instruct the library to calculate the recon to source for PSNR calculation
    *
    * Default is 0.*/
    uint32_t stat_report;

    /**
     * @brief API Signal to output reconstructed yuv used for debug purposes.
     * Using this will affect the speed of encoder.
     *
     * Default is false.
     */
    Bool recon_enabled;
    // 1.0.0: Any additional fields shall go after here

    /**
     * @brief Signal that force-key-frames is enabled.
     *
     */
    Bool force_key_frames;

    /**
     * @brief Signal to the library to treat intra_period_length as seconds and
     * multiply by fps_num/fps_den.
     */
    Bool multiply_keyint;

    // reference scaling parameters
    /**
     * @brief Reference scaling mode
     * the available modes are defined in RESIZE_MODE
     */
    uint8_t resize_mode;
    /**
     * @brief Resize denominator
     * this value can be from 8 to 16, means downscaling to 8/8-8/16 of original
     * resolution in both width and height
     */
    uint8_t resize_denom;
    /**
     * @brief Resize denominator of key frames
     * this value can be from 8 to 16, means downscaling to 8/8-8/16 of original
     * resolution in both width and height
     */
    uint8_t resize_kf_denom;
    /**
     * @brief Signal to the library to enable quantisation matrices
     *
     * Default is false.
     */
    Bool enable_qm;
    /**
     * @brief Min quant matrix flatness. Applicable when enable_qm is true.
     * Min value is 0.
     * Max value is 15.
     * Default is 8.
     */
    uint8_t min_qm_level;
    /**
     * @brief Max quant matrix flatness. Applicable when enable_qm is true.
     * Min value is 0.
     * Max value is 15.
     * Default is 15.
     */
    uint8_t max_qm_level;

    /**
     * @brief gop_constraint_rc
     *
     * Currently, only applicable for VBR and  when GoP size is greater than 119 frames.
     *
     * When enabled, the rate control matches the target rate for each GoP.
     *
     * 0: off
     * 1: on
     * Default is 0.
     */
    Bool gop_constraint_rc;

    /**
     * @brief scale factors for lambda value for different frame update types
     * factor >> 7 (/ 128) is the actual value in float
     */
    int32_t lambda_scale_factors[SVT_AV1_FRAME_UPDATE_TYPES];

    /* Dynamic gop
    *
    * 0 = disable Dynamic GoP
    * 1 = enable Dynamic GoP
    *  Default is 1. */
    Bool enable_dg;

    /**
     * @brief startup_mg_size
     *
     * When enabled, a MG with specified size will be inserted after the key frame.
     * The MG size is determined by 2^startup_mg_size.
     *
     * 0: off
     * 2: set hierarchical levels to 2 (MG size 4)
     * 3: set hierarchical levels to 3 (MG size 8)
     * 4: set hierarchical levels to 4 (MG size 16)
     * Default is 0.
     */
    uint8_t startup_mg_size;

    /* @brief reference scaling events for random access mode (resize-mode = 4)
     *
     * evt_num:          total count of events
     * start_frame_nums: array of scaling start frame numbers
     * resize_kf_denoms: array of scaling denominators of key-frame
     * resize_denoms:    array of scaling denominators of non-key-frame
     */
    SvtAv1FrameScaleEvts frame_scale_evts;

    /*Add 64 Byte Padding to Struct to avoid changing the size of the public configuration struct*/
    /* ROI map
    *
    * 0 = disable ROI
    * 1 = enable ROI
    *  Default is 0. */
    Bool    enable_roi_map;
    uint8_t padding[64 - sizeof(Bool)];

    AomFilmGrain *fgs_table;

    /* @brief Boost low variance regions according to a fast-growing formula
        0: no boost
        1: mild boost
        2: gentle boost
        3: medium boost
        4: aggressive boost
        5: maximum boost
       Default is 3 (medium curve) */
    uint8_t variance_boost_strength;

    /* @brief Enables the new 8x8-based variance algorithm, and picks an 8x8 variance value per superblock to determine boost
       Lower values enable detecting more blocks that need boosting, at the expense of more possible false positives (overall bitrate increase)
        0: disabled, uses classic 64x64 based variance algorithm instead
        1: enabled, 1st octile
        4: enabled, median
        8: enabled, maximum
       Default is 4
    */
    uint8_t new_variance_octile;
} EbSvtAv1EncConfiguration;

/**
 * Returns a string containing "v$tag-$commit_count-g$hash${dirty:+-dirty}"
 * @param[out] SVT_AV1_CVS_VERSION
 */
EB_API const char *svt_av1_get_version(void);

/**
 * Prints the version header and build information to the file
 * specified by the SVT_LOG_FILE environment variable or stderr
 */
EB_API void svt_av1_print_version(void);

/* STEP 1: Call the library to construct a Component Handle.
     *
     * Parameter:
     * @ **p_handle      Handle to be called in the future for manipulating the
     *                  component.
     * @ *p_app_data      Callback data.
     * @ *config_ptr     Pointer passed back to the client during callbacks, it will be
     *                  loaded with default params from the library. */
EB_API EbErrorType svt_av1_enc_init_handle(
    EbComponentType **p_handle, void *p_app_data,
    EbSvtAv1EncConfiguration *config_ptr); // config_ptr will be loaded with default params from the library

/* STEP 2: Set all configuration parameters.
     *
     * Parameter:
     * @ *svt_enc_component              Encoder handler.
     * @ *pComponentParameterStructure  Encoder and buffer configurations will be copied to the library. */
EB_API EbErrorType svt_av1_enc_set_parameter(
    EbComponentType *svt_enc_component,
    EbSvtAv1EncConfiguration
        *pComponentParameterStructure); // pComponentParameterStructure contents will be copied to the library

/* OPTIONAL: Set a single configuration parameter.
     *
     * Parameter:
     * @ *pComponentParameterStructure  Encoder parameters structure.
     * @ *name                          Null terminated string containing the parameter name
     * @ *value                         Null terminated string containing the parameter value */
EB_API EbErrorType svt_av1_enc_parse_parameter(EbSvtAv1EncConfiguration *pComponentParameterStructure, const char *name,
                                               const char *value);

/* STEP 3: Initialize encoder and allocates memory to necessary buffers.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_init(EbComponentType *svt_enc_component);

/* OPTIONAL: Get stream headers at init time.
     *
     * Parameter:
     * @ *svt_enc_component   Encoder handler.
     * @ **output_stream_ptr  Output buffer. */
EB_API EbErrorType svt_av1_enc_stream_header(EbComponentType     *svt_enc_component,
                                             EbBufferHeaderType **output_stream_ptr);

/* OPTIONAL: Release stream headers at init time.
     *
     * Parameter:
     * @ *stream_header_ptr  stream header buffer. */
EB_API EbErrorType svt_av1_enc_stream_header_release(EbBufferHeaderType *stream_header_ptr);

/* STEP 4: Send the picture.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *p_buffer           Header pointer, picture buffer. */
EB_API EbErrorType svt_av1_enc_send_picture(EbComponentType *svt_enc_component, EbBufferHeaderType *p_buffer);

/* STEP 5: Receive packet.
     * Parameter:
    * @ *svt_enc_component  Encoder handler.
     * @ **p_buffer          Header pointer to return packet with.
     * @ pic_send_done       Flag to signal that all input pictures have been sent, this call becomes locking one this signal is 1.
     * Non-locking call, returns EB_ErrorMax for an encode error, EB_NoErrorEmptyQueue when the library does not have any available packets.*/
EB_API EbErrorType svt_av1_enc_get_packet(EbComponentType *svt_enc_component, EbBufferHeaderType **p_buffer,
                                          uint8_t pic_send_done);

/* STEP 5-1: Release output buffer back into the pool.
     *
     * Parameter:
     * @ **p_buffer          Header pointer that contains the output packet to be released. */
EB_API void svt_av1_enc_release_out_buffer(EbBufferHeaderType **p_buffer);

/* OPTIONAL: Fill buffer with reconstructed picture.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *p_buffer           Output buffer. */
EB_API EbErrorType svt_av1_get_recon(EbComponentType *svt_enc_component, EbBufferHeaderType *p_buffer);

/* OPTIONAL: get stream information
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *stream_info_id SVT_AV1_STREAM_INFO_ID.
     * @ *info         output, the type depends on id */
EB_API EbErrorType svt_av1_enc_get_stream_info(EbComponentType *svt_enc_component, uint32_t stream_info_id, void *info);

/* STEP 6: Deinitialize encoder library.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_deinit(EbComponentType *svt_enc_component);

/* STEP 7: Deconstruct encoder handler.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_deinit_handle(EbComponentType *svt_enc_component);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbSvtAv1Enc_h
