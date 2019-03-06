/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef EbFormats_h
#define EbFormats_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


/*!\brief List of supported color primaries */
typedef enum aom_color_primaries {
    AOM_CICP_CP_RESERVED_0 = 0,  /**< For future use */
    AOM_CICP_CP_BT_709 = 1,      /**< BT.709 */
    AOM_CICP_CP_UNSPECIFIED = 2, /**< Unspecified */
    AOM_CICP_CP_RESERVED_3 = 3,  /**< For future use */
    AOM_CICP_CP_BT_470_M = 4,    /**< BT.470 System M (historical) */
    AOM_CICP_CP_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
    AOM_CICP_CP_BT_601 = 6,      /**< BT.601 */
    AOM_CICP_CP_SMPTE_240 = 7,   /**< SMPTE 240 */
    AOM_CICP_CP_GENERIC_FILM =
    8, /**< Generic film (color filters using illuminant C) */
    AOM_CICP_CP_BT_2020 = 9,      /**< BT.2020, BT.2100 */
    AOM_CICP_CP_XYZ = 10,         /**< SMPTE 428 (CIE 1921 XYZ) */
    AOM_CICP_CP_SMPTE_431 = 11,   /**< SMPTE RP 431-2 */
    AOM_CICP_CP_SMPTE_432 = 12,   /**< SMPTE EG 432-1  */
    AOM_CICP_CP_RESERVED_13 = 13, /**< For future use (values 13 - 21)  */
    AOM_CICP_CP_EBU_3213 = 22,    /**< EBU Tech. 3213-E  */
    AOM_CICP_CP_RESERVED_23 = 23  /**< For future use (values 23 - 255)  */
} aom_color_primaries_t;        /**< alias for enum aom_color_primaries */

/*!\brief List of supported transfer functions */
typedef enum aom_transfer_characteristics {
    AOM_CICP_TC_RESERVED_0 = 0,  /**< For future use */
    AOM_CICP_TC_BT_709 = 1,      /**< BT.709 */
    AOM_CICP_TC_UNSPECIFIED = 2, /**< Unspecified */
    AOM_CICP_TC_RESERVED_3 = 3,  /**< For future use */
    AOM_CICP_TC_BT_470_M = 4,    /**< BT.470 System M (historical)  */
    AOM_CICP_TC_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
    AOM_CICP_TC_BT_601 = 6,      /**< BT.601 */
    AOM_CICP_TC_SMPTE_240 = 7,   /**< SMPTE 240 M */
    AOM_CICP_TC_LINEAR = 8,      /**< Linear */
    AOM_CICP_TC_LOG_100 = 9,     /**< Logarithmic (100 : 1 range) */
    AOM_CICP_TC_LOG_100_SQRT10 =
    10,                     /**< Logarithmic (100 * Sqrt(10) : 1 range) */
    AOM_CICP_TC_IEC_61966 = 11, /**< IEC 61966-2-4 */
    AOM_CICP_TC_BT_1361 = 12,   /**< BT.1361 */
    AOM_CICP_TC_SRGB = 13,      /**< sRGB or sYCC*/
    AOM_CICP_TC_BT_2020_10_BIT = 14, /**< BT.2020 10-bit systems */
    AOM_CICP_TC_BT_2020_12_BIT = 15, /**< BT.2020 12-bit systems */
    AOM_CICP_TC_SMPTE_2084 = 16,     /**< SMPTE ST 2084, ITU BT.2100 PQ */
    AOM_CICP_TC_SMPTE_428 = 17,      /**< SMPTE ST 428 */
    AOM_CICP_TC_HLG = 18,            /**< BT.2100 HLG, ARIB STD-B67 */
    AOM_CICP_TC_RESERVED_19 = 19     /**< For future use (values 19-255) */
} aom_transfer_characteristics_t;  /**< alias for enum aom_transfer_function */

/*!\brief List of supported matrix coefficients */
typedef enum aom_matrix_coefficients {
    AOM_CICP_MC_IDENTITY = 0,    /**< Identity matrix */
    AOM_CICP_MC_BT_709 = 1,      /**< BT.709 */
    AOM_CICP_MC_UNSPECIFIED = 2, /**< Unspecified */
    AOM_CICP_MC_RESERVED_3 = 3,  /**< For future use */
    AOM_CICP_MC_FCC = 4,         /**< US FCC 73.628 */
    AOM_CICP_MC_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
    AOM_CICP_MC_BT_601 = 6,      /**< BT.601 */
    AOM_CICP_MC_SMPTE_240 = 7,   /**< SMPTE 240 M */
    AOM_CICP_MC_SMPTE_YCGCO = 8, /**< YCgCo */
    AOM_CICP_MC_BT_2020_NCL =
    9, /**< BT.2020 non-constant luminance, BT.2100 YCbCr  */
    AOM_CICP_MC_BT_2020_CL = 10, /**< BT.2020 constant luminance */
    AOM_CICP_MC_SMPTE_2085 = 11, /**< SMPTE ST 2085 YDzDx */
    AOM_CICP_MC_CHROMAT_NCL =
    12, /**< Chromaticity-derived non-constant luminance */
    AOM_CICP_MC_CHROMAT_CL = 13, /**< Chromaticity-derived constant luminance */
    AOM_CICP_MC_ICTCP = 14,      /**< BT.2100 ICtCp */
    AOM_CICP_MC_RESERVED_15 = 15 /**< For future use (values 15-255)  */
} aom_matrix_coefficients_t;

/*!\brief List of supported color range */
typedef enum aom_color_range {
    AOM_CR_STUDIO_RANGE = 0, /**< Y [16..235], UV [16..240] */
    AOM_CR_FULL_RANGE = 1    /**< YUV/RGB [0..255] */
} aom_color_range_t;       /**< alias for enum aom_color_range */

/* AV1 bit depth */
typedef enum EbBitDepth {
    EB_EIGHT_BIT = 8,
    EB_TEN_BIT  = 10,
    EB_TWELVE_BIT = 12    
} EbBitDepth;

/* AV1 Chroma Format */
typedef enum EbColorFormat {
    EB_YUV400,
    EB_YUV420,
    EB_YUV422,
    EB_YUV444
} EbColorFormat;

/*!\brief List of chroma sample positions */
typedef enum aom_chroma_sample_position {
    AOM_CSP_UNKNOWN = 0,          /**< Unknown */
    AOM_CSP_VERTICAL = 1,         /**< Horizontally co-located with luma(0, 0)*/
    /**< sample, between two vertical samples */
    AOM_CSP_COLOCATED = 2,        /**< Co-located with luma(0, 0) sample */
    AOM_CSP_RESERVED = 3          /**< Reserved value */
} aom_chroma_sample_position_t; /**< alias for enum aom_transfer_function */


#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbFormats_h