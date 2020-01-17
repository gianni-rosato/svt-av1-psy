/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file params.h
 *
 * @brief Define the default values, valid values and invalid values for
 * individual params of encoder configuration.
 *
 * @author Cidana-Edmond
 *
 * @version 1.0 <br>
 * 04-2019 Cidana-Edmond
 *
 ******************************************************************************/

#ifndef _TEST_PARAMS_H_
#define _TEST_PARAMS_H_

#include <vector>
#include "EbSvtAv1Enc.h"
#include "gtest/gtest.h"
#include "EbDefinitions.h"

using std::vector;

/** @brief Definitions of Marco to help get test value from test vectors for
 * individual parameter of encoder configuration
 *
 */
#define SIZE_VALID_PARAM(p) (svt_av1_test_params::valid_##p.size())
#define SIZE_INVALID_PARAM(p) (svt_av1_test_params::invalid_##p.size())
#define SIZE_SPECIAL_PARAM(p) (svt_av1_test_params::special_##p.size())
#define GET_DEFAULT_PARAM(p) (svt_av1_test_params::default_##p[0])
#define GET_VALID_PARAM(p, n) (svt_av1_test_params::valid_##p[n])
#define GET_INVALID_PARAM(p, n) (svt_av1_test_params::invalid_##p[n])
#define GET_SPECIAL_PARAM(p, n) (svt_av1_test_params::special_##p[n])

/** @defgroup svt_av1_test_params Test vectors of individual params test
 *  Defines test values of encoder params in test vectors, includes default,
 * valid, invalid and special case in vector
 *  @{
 */

namespace svt_av1_test_params {

/** @brief Test vectors for parameter enc_mode, inludes default, valid and
 * invalid <br>
 * paramter reference: <br>
 * A preset defining the quality vs density tradeoff point that the
 * encoding is to be performed at. 0 is the highest quality mode, 3 is
 * the highest density mode. <br>
 *
 * Default is defined as MAX_ENC_PRESET.
 *
 */
static const vector<uint8_t> default_enc_mode = {
    MAX_ENC_PRESET,
};
static const vector<uint8_t> valid_enc_mode = {
    ENC_M0, /**< highest quality mode */
    ENC_M1,
    ENC_M2,
    ENC_M3,
    ENC_M4,
    ENC_M5,
    ENC_M6,
    ENC_M7,
    MAX_ENC_PRESET,
};
static const vector<uint8_t> invalid_enc_mode = {
    MAX_ENC_PRESET + 1,
};

/* The intra period defines the interval of frames after which you insert an
 * Intra refresh. It is strongly recommended to set the value to multiple of
 * 8 minus 1 the closest to 1 second (e.g. 55, 47, 31, 23 should be used for
 * 60, 50, 30, (24 or 25) respectively.
 *
 * -1 = no intra update.
 * -2 = auto.
 *
 * Deault is -2. */
static const vector<int32_t> default_intra_period_length = {
    -2,
};
static const vector<int32_t> valid_intra_period_length = {
    23,  // 24, 25 fps
    31,  // 30 fps
    47,  // 50 fps
    55,  // 60 fps
    -1,  // no intra update
    -2,  // auto
    0,   // all I frame
};
static const vector<int32_t> invalid_intra_period_length = {
    -3,   // < -2
    256,  // > 255
};

/* Random access.
 *
 * 1 = CRA, open GOP.
 * 2 = IDR, closed GOP.
 *
 * Default is 1. */
static const vector<uint32_t> default_intra_refresh_type = {
    1,
};
static const vector<uint32_t> valid_intra_refresh_type = {
    1,  // CRA, open GOP
    2,  // IDR, closed GOP
};
static const vector<uint32_t> invalid_intra_refresh_type = {
    0,
    3,
};

/* Number of hierarchical layers used to construct GOP.
 * Minigop size = 2^HierarchicalLevels.
 *
 * Default is 3. */
static const vector<uint32_t> default_hierarchical_levels = {
    4,
};
static const vector<uint32_t> valid_hierarchical_levels = {3, 4};
static const vector<uint32_t> invalid_hierarchical_levels = {
    0, 1, 2, 5,  // ...
};

/* Prediction structure used to construct GOP. There are two main structures
 * supported, which are: Low Delay (P or b) and Random Access.
 *
 * In Low Delay structure, pictures within a mini GOP refer to the previously
 * encoded pictures in display order. In other words, pictures with display
 * order N can only be referenced by pictures with display order greater than
 * N, and it can only refer pictures with picture order lower than N. The Low
 * Delay structure can be flat structured (e.g. IPPPPPPP...) or hierarchically
 * structured. b/b pictures can be used instead of P/p pictures. However, the
 * reference picture list 0 and the reference picture list 1 will contain the
 * same reference picture.
 *
 * Following values are supported and defined in EbDefinitions.h
 * #define EB_PRED_LOW_DELAY_P     0
 * #define EB_PRED_LOW_DELAY_B     1
 * #define EB_PRED_RANDOM_ACCESS   2
 * #define EB_PRED_TOTAL_COUNT     3

 * In Random Access structure, the b/b pictures can refer to reference pictures
 * from both directions (past and future).
 *
 * Default is 2. */
static const vector<uint8_t> default_pred_structure = {
    EB_PRED_RANDOM_ACCESS,
};
static const vector<uint8_t> valid_pred_structure = {
    EB_PRED_LOW_DELAY_P, EB_PRED_LOW_DELAY_B, EB_PRED_RANDOM_ACCESS};
static const vector<uint8_t> invalid_pred_structure = {
    /* _pred_structure override in code
    EB_PRED_TOTAL_COUNT, EB_PRED_TOTAL_COUNT + 1, EB_PRED_INVALID*/};

/* Decides whether to use b picture or P picture in the base layer.
 *
 * 0 = b Picture.
 * 1 = P Picture.
 *
 * Default is 0. */
static const vector<uint32_t> default_base_layer_switch_mode = {
    0,
};
static const vector<uint32_t> valid_base_layer_switch_mode = {
    0,  // b Picture.
    1,  // P Picture
};
static const vector<uint32_t> invalid_base_layer_switch_mode = {
    2,  // > 1
};

// Input Info
/* The width of input source in units of picture luma pixels.
 *
 * Default is 0. */
static const vector<uint32_t> default_source_width = {
    0,
};
static const vector<uint32_t> valid_source_width = {
    64, 320, 640, 800, 1280, 1920, 2560, 3840, 4096,  // ...
};
static const vector<uint32_t> invalid_source_width = {
    0, 1, 2, 4, 8, 16, 32, 63, 65, 4097,  // ...
};

/* The height of input source in units of picture luma pixels.
 *
 * Default is 0. */
static const vector<uint32_t> default_source_height = {
    0,
};
static const vector<uint32_t> valid_source_height = {
    64, 240, 480, 720, 1080, 1440, 1600, 2160,  // ...
};
static const vector<uint32_t> invalid_source_height = {
    0, 1, 2, 4, 8, 16, 32, 63, 65, 2161,  // ...
};

/* The frequecy of images being displayed. If the number is less than 1000,
 * the input frame rate is an integer number between 1 and 60, else the input
 * number is in Q16 format, shifted by 16 bits, where max allowed is 240 fps.
 * If FrameRateNumerator and FrameRateDenominator are both not equal to zero,
 * the encoder will ignore this parameter.
 *
 * Default is 25. */
static const vector<uint32_t> default_frame_rate = {
    30 << 16,
};
static const vector<uint32_t> valid_frame_rate = {
    1,
    24,
    25,
    30,
    50,
    60,
    1 << 16,
    24 << 16,
    25 << 16,
    30 << 16,
    50 << 16,
    60 << 16,
    120 << 16,
    240 << 16,
};
static const vector<uint32_t> invalid_frame_rate = {
    0, 241 << 16, 0xFFFFFFFF  // ...
};

// TODO: follwoing two parameters should be a combination test
/* Frame rate numerator. When zero, the encoder will use -fps if
 * FrameRateDenominator is also zero, otherwise an error is returned.
 *
 * Default is 0. */
static const vector<uint32_t> default_frame_rate_numerator = {
    0,
};
static const vector<uint32_t> valid_frame_rate_numerator = {
    1,
    2,
    3,
    4,
    5,
    6,
    8,
    10,
    20,
    30,
    50,
    60,
    90,
    120,
    240,
};
static const vector<uint32_t> invalid_frame_rate_numerator = {
    0xFFFFFFFF,
};
static const vector<uint32_t> special_frame_rate_numerator = {
    0,
};

/* Frame rate denominator. When zero, the encoder will use -fps if
 * FrameRateNumerator is also zero, otherwise an error is returned.
 *
 * Default is 0. */
static const vector<uint32_t> default_frame_rate_denominator = {
    0,
};
static const vector<uint32_t> valid_frame_rate_denominator = {
    1,
    2,
    3,
    4,
    5,
    6,
    10,
    30,
    60,
};
static const vector<uint32_t> invalid_frame_rate_denominator = {
    0xFFFFFFFF,
};
static const vector<uint32_t> special_frame_rate_denominator = {
    0,
};

/* Specifies the bit depth of input video.
 *
 * 8 = 8 bit.
 * 10 = 10 bit.
 *
 * Default is 8. */
static const vector<uint32_t> default_encoder_bit_depth = {
    8,
};
static const vector<uint32_t> valid_encoder_bit_depth = {
    8,   // 8 bit.
    10,  // 10 bit.
};
static const vector<uint32_t> invalid_encoder_bit_depth = {
    0, 1, 2, 11,  // ...
};

/* Offline packing of the 2bits: requires two bits packed input.
 *
 * Default is 0. */
static const vector<uint32_t> default_compressed_ten_bit_format = {
    0,
};
static const vector<uint32_t> valid_compressed_ten_bit_format = {
    0,
    //1, Not supported in this version
};
static const vector<uint32_t> invalid_compressed_ten_bit_format = {
    2, 10,  // ...
};

/* Number of frames of sequence to be encoded. If number of frames is greater
 * than the number of frames in file, the encoder will loop to the beginning
 * and continue the encode.
 *
 * 0 = encodes the full clip.
 *
 * Default is 0. */
static const vector<uint64_t> default_frames_to_be_encoded = {
    0,
};
static const vector<uint64_t> valid_frames_to_be_encoded = {
    0, 1, 10, 100, 10000, (uint64_t)0xFFFFFFFFFFFFFFFF,  // ...
};
static const vector<uint64_t> invalid_frames_to_be_encoded = {
    // none
};

/* Super block size for motion estimation
 *
 * Default is 64. */
static const vector<uint32_t> default_sb_sz = {
    64,
};
static const vector<uint32_t> valid_sb_sz = {
    4,
    8,
    16,
    32,
    64,
    MAX_SB_SIZE,
};
static const vector<uint32_t> invalid_sb_sz = {
    /* sb_sz override in code
    0,
    (MAX_SB_SIZE + 1),
    (MAX_SB_SIZE << 1),
    */
};

/* Super block size
 * Only support 64 and 128
 * Default is 128. */
static const vector<uint32_t> default_super_block_size = {
    128,
};
static const vector<uint32_t> valid_super_block_size = {
    64,
    MAX_SB_SIZE,
};
static const vector<uint32_t> invalid_super_block_size = {
    /* super_block_size override in code
    0,
    (MAX_SB_SIZE + 1),
    (MAX_SB_SIZE << 1),
    */
};

/* The maximum partitioning depth with 0 being the superblock depth
 *
 * Default is 4. */
static const vector<uint32_t> default_partition_depth = {
    4,
};
static const vector<uint32_t> valid_partition_depth = {
    0,
    1,
    2,
    3,
    EB_MAX_SB_DEPTH,
};
static const vector<uint32_t> invalid_partition_depth = {
    /* partition_depth override in code
    (EB_MAX_SB_DEPTH + 1),
    */
};

// Quantization
/* Initial quantization parameter for the Intra pictures used under constant
 * qp rate control mode.
 *
 * Default is 50. */
static const vector<uint32_t> default_qp = {
    50,
};
static const vector<uint32_t> valid_qp = {
    MIN_QP_VALUE,
    1,
    2,
    10,
    25,
    32,
    50,
    56,
    62,
    MAX_QP_VALUE,
};
static const vector<uint32_t> invalid_qp = {
    (MAX_QP_VALUE + 1),
};

/* force qp values for every picture that are passed in the header pointer
 *
 * Default is 0.*/
static const vector<EbBool> default_use_qp_file = {
    EB_FALSE,
};
static const vector<EbBool> valid_use_qp_file = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_use_qp_file = {
    // none
};

/* Enable picture QP scaling between hierarchical levels
 *
 * Default is null.*/
// TODO: the description of this parameter is incorrect, refer to source code,
// it should be a EbBool
static const vector<uint32_t> default_enable_qp_scaling_flag = {
    EB_FALSE,
};
static const vector<uint32_t> valid_enable_qp_scaling_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<uint32_t> invalid_enable_qp_scaling_flag = {
    // none
};

// Deblock Filter
/* Flag to disable the Deblocking Loop Filtering.
 *
 * Default is 0. */
static const vector<EbBool> default_disable_dlf_flag = {
    EB_FALSE,
};
static const vector<EbBool> valid_disable_dlf_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_disable_dlf_flag = {
    // none
};

/* Denoise the input picture when noise levels are too high
 * Flag to enable the denoising
 *
 * Default is 0. */
static const vector<EbBool> default_enable_denoise_flag = {
    EB_FALSE,
};
static const vector<EbBool> valid_enable_denoise_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_enable_denoise_flag = {
    // none
};

/* Film grain denoising the input picture
 * Flag to enable the denoising
 *
 * Default is 0. */
// TODO: the description of this parameter is incorrect, refer to source code,
// it should be a EbBool
static const vector<uint32_t> default_film_grain_denoise_strength = {
    EB_FALSE,
};
static const vector<uint32_t> valid_film_grain_denoise_strength = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<uint32_t> invalid_film_grain_denoise_strength = {
    // none
};

/* Warped motion
 *
 * Default is 0. */
static const vector<EbBool> default_enable_warped_motion = {
    EB_TRUE,
};
static const vector<EbBool> valid_enable_warped_motion = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_enable_warped_motion = {
    // none
};

/* Global motion
 *
 * Default is 1. */
static const vector<EbBool> default_enable_global_motion = {
    EB_TRUE,
};
static const vector<EbBool> valid_enable_global_motion = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_enable_global_motion = {
    // none
};

/* Flag to enable the use of default ME HME parameters.
 *
 * Default is 1. */
static const vector<EbBool> default_use_default_me_hme = {
    EB_TRUE,
};
static const vector<EbBool> valid_use_default_me_hme = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_use_default_me_hme = {
    // none
};

/* Flag to enable Hierarchical Motion Estimation.
 *
 * Default is 1. */
static const vector<EbBool> default_enable_hme_flag = {
    EB_TRUE,
};
static const vector<EbBool> valid_enable_hme_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_enable_hme_flag = {
    // none
};

/* Flag to enable the use of non-swaure partitions
 *
 * Default is 0. */
static const vector<EbBool> default_ext_block_flag = {
    EB_FALSE,
};
static const vector<EbBool> valid_ext_block_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_ext_block_flag = {
    // none
};

/* Flag to enable the use of recon pictures for motion estimation
 *
 * Default is 0. */
static const vector<EbBool> default_in_loop_me_flag = {
    EB_FALSE,
};
static const vector<EbBool> valid_in_loop_me_flag = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_in_loop_me_flag = {
    // none
};

// ME Parameters
/* Number of search positions in the horizontal direction.
 *
 * Default depends on input resolution. */
static const vector<uint32_t> default_search_area_width = {
    16,  // 0,
};
static const vector<uint32_t> valid_search_area_width = {
    1, 2, 3, 4, 8, 10, 16, 32, 64, 128, 256,  // ...
};
static const vector<uint32_t> invalid_search_area_width = {
    0, 257, 1000,  // ...
};

/* Number of search positions in the vertical direction.
 *
 * Default depends on input resolution. */
static const vector<uint32_t> default_search_area_height = {
    7,  // 0,
};
static const vector<uint32_t> valid_search_area_height = {
    1, 2, 3, 4, 8, 10, 16, 32, 64, 128, 256,  // ...
};
static const vector<uint32_t> invalid_search_area_height = {
    0, 257, 1000,  // ...
};

// MD Parameters
/* Palette Mode
 *-1:Auto Mode(ON at level6 when SC is detected)
 * 0:OFF
 * 1:Slow    NIC=7/4/4
 * 2:        NIC=7/2/2
 * 3:        NIC=7/2/2 + No K means for non ref
 * 4:        NIC=4/2/1
 * 5:        NIC=4/2/1 + No K means for Inter frame
 * 6:Fastest NIC=4/2/1 + No K means for non base + step for non base for
 * most dominant
 * Default is -1. */
static const vector<int32_t> default_enable_palette = {-1};
static const vector<int32_t> valid_enable_palette = {-1, 0, 1, 2, 3, 4, 5, 6};
static const vector<int32_t> invalid_enable_palette = {-2, 7};

/* Enable the use of Constrained Intra, which yields sending two picture
 * parameter sets in the elementary streams .
 *
 * Default is 0. */
static const vector<EbBool> default_constrained_intra = {
    EB_FALSE,
};
static const vector<EbBool> valid_constrained_intra = {
    EB_FALSE,
    EB_TRUE,
};
static const vector<EbBool> invalid_constrained_intra = {
    // none
};

// Rate Control
/* Rate control mode.
 *
 * 0 = Constant QP.
 * 1 = Average BitRate.
 *
 * Default is 0. */
static const vector<uint32_t> default_rate_control_mode = {0};
static const vector<uint32_t> valid_rate_control_mode = {0, 1, 2};
static const vector<uint32_t> invalid_rate_control_mode = {3, 4};

/* Flag to enable the scene change detection algorithm.
 *
 * Default is 0. */
static const vector<uint32_t> default_scene_change_detection = {
    0,
};
static const vector<uint32_t> valid_scene_change_detection = {
    0,
    1,
};
static const vector<uint32_t> invalid_scene_change_detection = {
    2,
};

/* When RateControlMode is set to 1 it's best to set this parameter to be
 * equal to the Intra period value (such is the default set by the encoder).
 * When CQP is chosen, then a (2 * minigopsize +1) look ahead is recommended.
 *
 * Default depends on rate control mode.*/
static const vector<uint32_t> default_look_ahead_distance = {
    (uint32_t)~0,
};
static const vector<uint32_t> valid_look_ahead_distance = {
    (uint32_t)~0,
    0,
    1,
    10,
    24,
    25,
    30,
    60,
    MAX_LAD,
};
static const vector<uint32_t> invalid_look_ahead_distance = {
    /* look_ahead_distance override in code
    MAX_LAD + 1, ((uint32_t)~0 - 1)
    */
    };

/* Target bitrate in bits/second, only apllicable when rate control mode is
 * set to 1.
 *
 * Default is 7000000. */
static const vector<uint32_t> default_target_bit_rate = {
    7000000,
};
static const vector<uint32_t> valid_target_bit_rate = {
    0,
    1,
    100,
    1000,
    10000,
    100000,
    1000000,
    7000000,
    10000000,
    0xFFFFFFFF,
};
static const vector<uint32_t> invalid_target_bit_rate = {
    // none
};

/* Maxium QP value allowed for rate control use, only applicable when rate
 * control mode is set to 1. It has to be greater or equal to minQpAllowed.
 *
 * Default is 63. */
static const vector<uint32_t> default_max_qp_allowed = {
    MAX_QP_VALUE,
};
static const vector<uint32_t> valid_max_qp_allowed = {
    MIN_QP_VALUE,
    1,
    2,
    10,
    25,
    32,
    50,
    56,
    62,
    MAX_QP_VALUE,
};
static const vector<uint32_t> invalid_max_qp_allowed = {
    (MAX_QP_VALUE + 1),
};

/* Minimum QP value allowed for rate control use, only applicable when rate
 * control mode is set to 1. It has to be smaller or equal to maxQpAllowed.
 *
 * Default is 0. */
/*
 * There is a value check for min_qp_allowed in EbEncHandle.c :
 * else if (config->min_qp_allowed >= MAX_QP_VALUE) {
 *     SVT_LOG("Error instance %u: MinQpAllowed must be [0 - %d]\n",
 *         channel_number + 1, MAX_QP_VALUE-1); return_error =
 *         EB_ErrorBadParameter;
 * }
 * The maximum valid value should be MAX_QP_VALUE - 10.
 */
static const vector<uint32_t> default_min_qp_allowed = {
    10,
};
static const vector<uint32_t> valid_min_qp_allowed = {
    MIN_QP_VALUE,
    1,
    2,
    10,
    25,
    32,
    50,
    56,
    MAX_QP_VALUE - 1,
};
static const vector<uint32_t> invalid_min_qp_allowed = {
    MAX_QP_VALUE,
};

// Tresholds
/* Flag to signal that the input yuv is HDR10 BT2020 using SMPTE ST2048,
 * requires
 *
 * Default is 0. */
static const vector<uint32_t> default_high_dynamic_range_input = {
    0,
};
static const vector<uint32_t> valid_high_dynamic_range_input = {
    0,
    1,
};
static const vector<uint32_t> invalid_high_dynamic_range_input = {
    2,
};

/* Defined set of coding tools to create Bitstream.
 *
 * 1 = Main, allows bit depth of 8.
 * 2 = Main 10, allows bit depth of 8 to 10.
 *
 * Default is 0. */
static const vector<uint32_t> default_profile = {
    PROFILE_0,
};
static const vector<uint32_t> valid_profile = {
    PROFILE_0,
    PROFILE_1,
    PROFILE_2,
};
static const vector<uint32_t> invalid_profile = {
    MAX_PROFILES,
};

/* Constraints for Bitstream in terms of max bitrate and max buffer size.
 *
 * 0 = Main, for most applications.
 * 1 = High, for demanding applications.
 *
 * Default is 0. */
static const vector<uint32_t> default_tier = {
    0,
};
static const vector<uint32_t> valid_tier = {
    0,
    1,
};
static const vector<uint32_t> invalid_tier = {
    /* tier override in code
    2,
    */
};

/* Constraints for Bitstream in terms of max bitrate and max buffer size.
 *
 * 0 = auto determination.
 *
 * Default is 0. */
static const vector<uint32_t> default_level = {
    0,
};
static const vector<uint32_t> valid_level = {0, 1, 2, 10, 64, 100, 0xFFFFFFFF};
static const vector<uint32_t> invalid_level = {
    // none
};

/* Assembly instruction set used by encoder.
 *
 * 0 = non-AVX2, C only.
 * 1 = up to AVX512, auto-select highest assembly instruction set supported.
 *
 * Default is 1. */
static const vector<CPU_FLAGS> default_use_cpu_flags = {
    CPU_FLAGS_ALL,
};
static const vector<CPU_FLAGS> valid_use_cpu_flags = {
    CPU_FLAGS_MMX,
    CPU_FLAGS_SSE,
    CPU_FLAGS_SSE2,
    CPU_FLAGS_SSE3,
    CPU_FLAGS_SSSE3,
    CPU_FLAGS_SSE4_1,
    CPU_FLAGS_SSE4_2,
    CPU_FLAGS_AVX,
    CPU_FLAGS_AVX2,
    CPU_FLAGS_AVX512F,
    CPU_FLAGS_AVX512CD,
    CPU_FLAGS_AVX512DQ,
    CPU_FLAGS_AVX512ER,
    CPU_FLAGS_AVX512PF,
    CPU_FLAGS_AVX512BW,
    CPU_FLAGS_AVX512VL,
};
static const vector<CPU_FLAGS> invalid_use_cpu_flags = {
    CPU_FLAGS_INVALID
};

// Application Specific parameters
/* ID assigned to each channel when multiple instances are running within the
 * same application. */
static const vector<uint32_t> default_channel_id = {
    0,
};
static const vector<uint32_t> valid_channel_id = {
    0,
    1,
    2,
    3,
    10,
    100,
    0xFFFFFFFF,
};
static const vector<uint32_t> invalid_channel_id = {
    // none
};

static const vector<uint32_t> default_active_channel_count = {
    1,
};
static const vector<uint32_t> valid_active_channel_count = {
    1,
    2,
    3,
    10,
    100,
    0xFFFFFFFF,
};
static const vector<uint32_t> invalid_active_channel_count = {
    /* active_channel_count override in code
    0,
    */
};

/* Flag to enable the Speed Control functionality to achieve the real-time
 * encoding speed defined by dynamically changing the encoding preset to meet
 * the average speed defined in injectorFrameRate. When this parameter is set
 * to 1 it forces -inj to be 1 -inj-frm-rt to be set to the -fps.
 *
 * Default is 0. */
static const vector<uint32_t> default_speed_control_flag = {
    0,
};
static const vector<uint32_t> valid_speed_control_flag = {
    0,
    1,
};
static const vector<uint32_t> invalid_speed_control_flag = {
    2,
};

/* Frame Rate used for the injector. Recommended to match the encoder speed.
 *
 * Default is 60. */
static const vector<int32_t> default_injector_frame_rate = {
    60<<16,
};
static const vector<int32_t> valid_injector_frame_rate = {
    24,
    25,
    30,
    50,
    60,
    120,
    240,
    24 << 16,
    25 << 16,
    30 << 16,
    50 << 16,
    60 << 16,
    120 << 16,
    240 << 16,  // ...
};
static const vector<int32_t> invalid_injector_frame_rate = {
    /* injector_frame_rate override in code
    0, 1, 2, 10, 15, 29, 241,  // ...
    */
};

// Threads management

/* The number of logical processor which encoder threads run on. If
 * LogicalProcessorNumber and TargetSocket are not set, threads are managed by
 * OS thread scheduler. */
static const vector<uint32_t> default_logical_processors = {
    0,
};
static const vector<uint32_t> valid_logical_processors = {
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    20,
    40,
    1000,
    0xFFFFFFFF,
};
static const vector<uint32_t> invalid_logical_processors = {
    // ...
};

/* Target socket to run on. For dual socket systems, this can specify which
 * socket the encoder runs on.
 *
 * -1 = Both Sockets.
 *  0 = Socket 0.
 *  1 = Socket 1.
 *
 * Default is -1. */
static const vector<int32_t> default_target_socket = {
    -1,
};
static const vector<int32_t> valid_target_socket = {
    -1,
    0,
    1,
};
static const vector<int32_t> invalid_target_socket = {
    2,
};

// Debug tools

/* Output reconstructed yuv used for debug purposes. The value is set through
 * ReconFile token (-o) and using the feature will affect the speed of encoder.
 *
 * Default is 0. */
static const vector<uint32_t> default_recon_enabled = {EB_FALSE};
static const vector<uint32_t> valid_recon_enabled = {EB_FALSE, EB_TRUE};
static const vector<uint32_t> invalid_recon_enabled = {/** none */};

#if TILES
/* Log 2 Tile Rows and colums . 0 means no tiling,1 means that we split the
 * dimension into 2 Default is 0. */
static const vector<int32_t> default_tile_columns = {
    0,
};
static const vector<int32_t> valid_tile_columns = {
    0,
    1,
    2,
    3,
    4,
    5,
    6,
};
static const vector<int32_t> invalid_tile_columns = {
    -1,
    7,
};

static const vector<int32_t> default_tile_rows = {
    0,
};
static const vector<int32_t> valid_tile_rows = {
    0,
    1,
    2,
    3,
    4,
    5,
    6,
};
static const vector<int32_t> invalid_tile_rows = {
    -1,
    7,
};
#endif

/* Flag to signal the content being a screen sharing content type
 *
 * Default is 2. */
static const vector<uint32_t> default_screen_content_mode = {2};
static const vector<uint32_t> valid_screen_content_mode = {0, 1, 2};
static const vector<uint32_t> invalid_screen_content_mode = {3};

/* Variables to control the use of ALT-REF (temporally filtered frames)
 */
static const vector<EbBool> default_enable_altrefs = {EB_TRUE};
static const vector<EbBool> valid_enable_altrefs = {EB_FALSE, EB_TRUE};
static const vector<EbBool> invalid_enable_altrefs = {/*none*/};

static const vector<uint8_t> default_altref_strength = {5};
static const vector<uint8_t> valid_altref_strength = {0, 1, 2, 3, 4, 5, 6};
static const vector<uint8_t> invalid_altref_strength = {7};

static const vector<uint8_t> default_altref_nframes = {7};
static const vector<uint8_t> valid_altref_nframes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
static const vector<uint8_t> invalid_altref_nframes = {11};

static const vector<EbBool> default_enable_overlays = {EB_FALSE};
static const vector<EbBool> valid_enable_overlays = {EB_FALSE, EB_TRUE};
static const vector<EbBool> invalid_enable_overlays = {/*none*/};

}  // namespace svt_av1_test_params

/** @} */  // end of svt_av1_test_params

#endif  // _TEST_PARAMS_H_
