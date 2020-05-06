/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimationContext_h
#define EbMotionEstimationContext_h

#include "EbDefinitions.h"
#include "EbMdRateEstimation.h"
#include "EbCodingUnit.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif

// Max Search Area
#define MAX_SEARCH_AREA_WIDTH 1280
#define MAX_SEARCH_AREA_HEIGHT 1280
#define MAX_SEARCH_AREA_WIDTH_CH MAX_SEARCH_AREA_WIDTH + PAD_VALUE
#define MAX_SEARCH_AREA_HEIGHT_CH MAX_SEARCH_AREA_HEIGHT + PAD_VALUE

// 1-D interpolation shift value
#define if_shift 6
#define NUMBER_OF_SB_QUAD 4
#define VARIANCE_PRECISION 16
#define MEAN_PRECISION (VARIANCE_PRECISION >> 1)
#define HME_RECTANGULAR 0
#define HME_SPARSE 1
#define HME_DECIM_FILTER_TAP 9

// Quater pel refinement methods
typedef enum EbQuarterPelRefinementMethod {
    EB_QUARTER_IN_FULL,
    EB_QUARTER_IN_HALF_HORIZONTAL,
    EB_QUARTER_IN_HALF_VERTICAL,
    EB_QUARTER_IN_HALF_DIAGONAL
} EbQuarterPelInterpolationMethod;

typedef struct MePredictionUnit {
    uint64_t distortion;
    int16_t  x_mv;
    int16_t  y_mv;
    uint32_t sub_pel_direction;
} MePredictionUnit;

typedef enum EbMeTierZeroPu {
    // 2Nx2N [85 partitions]
    ME_TIER_ZERO_PU_64x64    = 0,
    ME_TIER_ZERO_PU_32x32_0  = 1,
    ME_TIER_ZERO_PU_32x32_1  = 2,
    ME_TIER_ZERO_PU_32x32_2  = 3,
    ME_TIER_ZERO_PU_32x32_3  = 4,
    ME_TIER_ZERO_PU_16x16_0  = 5,
    ME_TIER_ZERO_PU_16x16_1  = 6,
    ME_TIER_ZERO_PU_16x16_2  = 7,
    ME_TIER_ZERO_PU_16x16_3  = 8,
    ME_TIER_ZERO_PU_16x16_4  = 9,
    ME_TIER_ZERO_PU_16x16_5  = 10,
    ME_TIER_ZERO_PU_16x16_6  = 11,
    ME_TIER_ZERO_PU_16x16_7  = 12,
    ME_TIER_ZERO_PU_16x16_8  = 13,
    ME_TIER_ZERO_PU_16x16_9  = 14,
    ME_TIER_ZERO_PU_16x16_10 = 15,
    ME_TIER_ZERO_PU_16x16_11 = 16,
    ME_TIER_ZERO_PU_16x16_12 = 17,
    ME_TIER_ZERO_PU_16x16_13 = 18,
    ME_TIER_ZERO_PU_16x16_14 = 19,
    ME_TIER_ZERO_PU_16x16_15 = 20,
    ME_TIER_ZERO_PU_8x8_0    = 21,
    ME_TIER_ZERO_PU_8x8_1    = 22,
    ME_TIER_ZERO_PU_8x8_2    = 23,
    ME_TIER_ZERO_PU_8x8_3    = 24,
    ME_TIER_ZERO_PU_8x8_4    = 25,
    ME_TIER_ZERO_PU_8x8_5    = 26,
    ME_TIER_ZERO_PU_8x8_6    = 27,
    ME_TIER_ZERO_PU_8x8_7    = 28,
    ME_TIER_ZERO_PU_8x8_8    = 29,
    ME_TIER_ZERO_PU_8x8_9    = 30,
    ME_TIER_ZERO_PU_8x8_10   = 31,
    ME_TIER_ZERO_PU_8x8_11   = 32,
    ME_TIER_ZERO_PU_8x8_12   = 33,
    ME_TIER_ZERO_PU_8x8_13   = 34,
    ME_TIER_ZERO_PU_8x8_14   = 35,
    ME_TIER_ZERO_PU_8x8_15   = 36,
    ME_TIER_ZERO_PU_8x8_16   = 37,
    ME_TIER_ZERO_PU_8x8_17   = 38,
    ME_TIER_ZERO_PU_8x8_18   = 39,
    ME_TIER_ZERO_PU_8x8_19   = 40,
    ME_TIER_ZERO_PU_8x8_20   = 41,
    ME_TIER_ZERO_PU_8x8_21   = 42,
    ME_TIER_ZERO_PU_8x8_22   = 43,
    ME_TIER_ZERO_PU_8x8_23   = 44,
    ME_TIER_ZERO_PU_8x8_24   = 45,
    ME_TIER_ZERO_PU_8x8_25   = 46,
    ME_TIER_ZERO_PU_8x8_26   = 47,
    ME_TIER_ZERO_PU_8x8_27   = 48,
    ME_TIER_ZERO_PU_8x8_28   = 49,
    ME_TIER_ZERO_PU_8x8_29   = 50,
    ME_TIER_ZERO_PU_8x8_30   = 51,
    ME_TIER_ZERO_PU_8x8_31   = 52,
    ME_TIER_ZERO_PU_8x8_32   = 53,
    ME_TIER_ZERO_PU_8x8_33   = 54,
    ME_TIER_ZERO_PU_8x8_34   = 55,
    ME_TIER_ZERO_PU_8x8_35   = 56,
    ME_TIER_ZERO_PU_8x8_36   = 57,
    ME_TIER_ZERO_PU_8x8_37   = 58,
    ME_TIER_ZERO_PU_8x8_38   = 59,
    ME_TIER_ZERO_PU_8x8_39   = 60,
    ME_TIER_ZERO_PU_8x8_40   = 61,
    ME_TIER_ZERO_PU_8x8_41   = 62,
    ME_TIER_ZERO_PU_8x8_42   = 63,
    ME_TIER_ZERO_PU_8x8_43   = 64,
    ME_TIER_ZERO_PU_8x8_44   = 65,
    ME_TIER_ZERO_PU_8x8_45   = 66,
    ME_TIER_ZERO_PU_8x8_46   = 67,
    ME_TIER_ZERO_PU_8x8_47   = 68,
    ME_TIER_ZERO_PU_8x8_48   = 69,
    ME_TIER_ZERO_PU_8x8_49   = 70,
    ME_TIER_ZERO_PU_8x8_50   = 71,
    ME_TIER_ZERO_PU_8x8_51   = 72,
    ME_TIER_ZERO_PU_8x8_52   = 73,
    ME_TIER_ZERO_PU_8x8_53   = 74,
    ME_TIER_ZERO_PU_8x8_54   = 75,
    ME_TIER_ZERO_PU_8x8_55   = 76,
    ME_TIER_ZERO_PU_8x8_56   = 77,
    ME_TIER_ZERO_PU_8x8_57   = 78,
    ME_TIER_ZERO_PU_8x8_58   = 79,
    ME_TIER_ZERO_PU_8x8_59   = 80,
    ME_TIER_ZERO_PU_8x8_60   = 81,
    ME_TIER_ZERO_PU_8x8_61   = 82,
    ME_TIER_ZERO_PU_8x8_62   = 83,
    ME_TIER_ZERO_PU_8x8_63   = 84,
    // H  [42 partitions]
    ME_TIER_ZERO_PU_64x32_0 = 85,
    ME_TIER_ZERO_PU_64x32_1 = 86,
    ME_TIER_ZERO_PU_32x16_0 = 87,
    ME_TIER_ZERO_PU_32x16_1 = 88,
    ME_TIER_ZERO_PU_32x16_2 = 89,
    ME_TIER_ZERO_PU_32x16_3 = 90,
    ME_TIER_ZERO_PU_32x16_4 = 91,
    ME_TIER_ZERO_PU_32x16_5 = 92,
    ME_TIER_ZERO_PU_32x16_6 = 93,
    ME_TIER_ZERO_PU_32x16_7 = 94,
    ME_TIER_ZERO_PU_16x8_0  = 95,
    ME_TIER_ZERO_PU_16x8_1  = 96,
    ME_TIER_ZERO_PU_16x8_2  = 97,
    ME_TIER_ZERO_PU_16x8_3  = 98,
    ME_TIER_ZERO_PU_16x8_4  = 99,
    ME_TIER_ZERO_PU_16x8_5  = 100,
    ME_TIER_ZERO_PU_16x8_6  = 101,
    ME_TIER_ZERO_PU_16x8_7  = 102,
    ME_TIER_ZERO_PU_16x8_8  = 103,
    ME_TIER_ZERO_PU_16x8_9  = 104,
    ME_TIER_ZERO_PU_16x8_10 = 105,
    ME_TIER_ZERO_PU_16x8_11 = 106,
    ME_TIER_ZERO_PU_16x8_12 = 107,
    ME_TIER_ZERO_PU_16x8_13 = 108,
    ME_TIER_ZERO_PU_16x8_14 = 109,
    ME_TIER_ZERO_PU_16x8_15 = 110,
    ME_TIER_ZERO_PU_16x8_16 = 111,
    ME_TIER_ZERO_PU_16x8_17 = 112,
    ME_TIER_ZERO_PU_16x8_18 = 113,
    ME_TIER_ZERO_PU_16x8_19 = 114,
    ME_TIER_ZERO_PU_16x8_20 = 115,
    ME_TIER_ZERO_PU_16x8_21 = 116,
    ME_TIER_ZERO_PU_16x8_22 = 117,
    ME_TIER_ZERO_PU_16x8_23 = 118,
    ME_TIER_ZERO_PU_16x8_24 = 119,
    ME_TIER_ZERO_PU_16x8_25 = 120,
    ME_TIER_ZERO_PU_16x8_26 = 121,
    ME_TIER_ZERO_PU_16x8_27 = 122,
    ME_TIER_ZERO_PU_16x8_28 = 123,
    ME_TIER_ZERO_PU_16x8_29 = 124,
    ME_TIER_ZERO_PU_16x8_30 = 125,
    ME_TIER_ZERO_PU_16x8_31 = 126,
    // V  [42 partitions]
    ME_TIER_ZERO_PU_32x64_0 = 127,
    ME_TIER_ZERO_PU_32x64_1 = 128,
    ME_TIER_ZERO_PU_16x32_0 = 129,
    ME_TIER_ZERO_PU_16x32_1 = 130,
    ME_TIER_ZERO_PU_16x32_2 = 131,
    ME_TIER_ZERO_PU_16x32_3 = 132,
    ME_TIER_ZERO_PU_16x32_4 = 133,
    ME_TIER_ZERO_PU_16x32_5 = 134,
    ME_TIER_ZERO_PU_16x32_6 = 135,
    ME_TIER_ZERO_PU_16x32_7 = 136,
    ME_TIER_ZERO_PU_8x16_0  = 137,
    ME_TIER_ZERO_PU_8x16_1  = 138,
    ME_TIER_ZERO_PU_8x16_2  = 139,
    ME_TIER_ZERO_PU_8x16_3  = 140,
    ME_TIER_ZERO_PU_8x16_4  = 141,
    ME_TIER_ZERO_PU_8x16_5  = 142,
    ME_TIER_ZERO_PU_8x16_6  = 143,
    ME_TIER_ZERO_PU_8x16_7  = 144,
    ME_TIER_ZERO_PU_8x16_8  = 145,
    ME_TIER_ZERO_PU_8x16_9  = 146,
    ME_TIER_ZERO_PU_8x16_10 = 147,
    ME_TIER_ZERO_PU_8x16_11 = 148,
    ME_TIER_ZERO_PU_8x16_12 = 149,
    ME_TIER_ZERO_PU_8x16_13 = 150,
    ME_TIER_ZERO_PU_8x16_14 = 151,
    ME_TIER_ZERO_PU_8x16_15 = 152,
    ME_TIER_ZERO_PU_8x16_16 = 153,
    ME_TIER_ZERO_PU_8x16_17 = 154,
    ME_TIER_ZERO_PU_8x16_18 = 155,
    ME_TIER_ZERO_PU_8x16_19 = 156,
    ME_TIER_ZERO_PU_8x16_20 = 157,
    ME_TIER_ZERO_PU_8x16_21 = 158,
    ME_TIER_ZERO_PU_8x16_22 = 159,
    ME_TIER_ZERO_PU_8x16_23 = 160,
    ME_TIER_ZERO_PU_8x16_24 = 161,
    ME_TIER_ZERO_PU_8x16_25 = 162,
    ME_TIER_ZERO_PU_8x16_26 = 163,
    ME_TIER_ZERO_PU_8x16_27 = 164,
    ME_TIER_ZERO_PU_8x16_28 = 165,
    ME_TIER_ZERO_PU_8x16_29 = 166,
    ME_TIER_ZERO_PU_8x16_30 = 167,
    ME_TIER_ZERO_PU_8x16_31 = 168,
    // H4 [16 partitions]
    ME_TIER_ZERO_PU_32x8_0  = 169,
    ME_TIER_ZERO_PU_32x8_1  = 170,
    ME_TIER_ZERO_PU_32x8_2  = 171,
    ME_TIER_ZERO_PU_32x8_3  = 172,
    ME_TIER_ZERO_PU_32x8_4  = 173,
    ME_TIER_ZERO_PU_32x8_5  = 174,
    ME_TIER_ZERO_PU_32x8_6  = 175,
    ME_TIER_ZERO_PU_32x8_7  = 176,
    ME_TIER_ZERO_PU_32x8_8  = 177,
    ME_TIER_ZERO_PU_32x8_9  = 178,
    ME_TIER_ZERO_PU_32x8_10 = 179,
    ME_TIER_ZERO_PU_32x8_11 = 180,
    ME_TIER_ZERO_PU_32x8_12 = 181,
    ME_TIER_ZERO_PU_32x8_13 = 182,
    ME_TIER_ZERO_PU_32x8_14 = 183,
    ME_TIER_ZERO_PU_32x8_15 = 184,
    // V4 [16 partitions]
    ME_TIER_ZERO_PU_8x32_0  = 185,
    ME_TIER_ZERO_PU_8x32_1  = 186,
    ME_TIER_ZERO_PU_8x32_2  = 187,
    ME_TIER_ZERO_PU_8x32_3  = 188,
    ME_TIER_ZERO_PU_8x32_4  = 189,
    ME_TIER_ZERO_PU_8x32_5  = 190,
    ME_TIER_ZERO_PU_8x32_6  = 191,
    ME_TIER_ZERO_PU_8x32_7  = 192,
    ME_TIER_ZERO_PU_8x32_8  = 193,
    ME_TIER_ZERO_PU_8x32_9  = 194,
    ME_TIER_ZERO_PU_8x32_10 = 195,
    ME_TIER_ZERO_PU_8x32_11 = 196,
    ME_TIER_ZERO_PU_8x32_12 = 197,
    ME_TIER_ZERO_PU_8x32_13 = 198,
    ME_TIER_ZERO_PU_8x32_14 = 199,
    ME_TIER_ZERO_PU_8x32_15 = 200,
    ME_TIER_ZERO_PU_64x16_0 = 201,
    ME_TIER_ZERO_PU_64x16_1 = 202,
    ME_TIER_ZERO_PU_64x16_2 = 203,
    ME_TIER_ZERO_PU_64x16_3 = 204,
    ME_TIER_ZERO_PU_16x64_0 = 205,
    ME_TIER_ZERO_PU_16x64_1 = 206,
    ME_TIER_ZERO_PU_16x64_2 = 207,
    ME_TIER_ZERO_PU_16x64_3 = 208
} EbMeTierZeroPu;

typedef struct MeTierZero {
    MePredictionUnit pu[MAX_ME_PU_COUNT];
} MeTierZero;

typedef struct IntraReferenceSamplesOpenLoop {
    EbDctor  dctor;
    uint8_t *y_intra_reference_array_reverse;

    // Scratch buffers used in the interpolaiton process
    uint8_t reference_above_line_y[MAX_INTRA_REFERENCE_SAMPLES];
    uint8_t reference_left_line_y[MAX_INTRA_REFERENCE_SAMPLES];
    EbBool  above_ready_flag_y;
    EbBool  left_ready_flag_y;
} IntraReferenceSamplesOpenLoop;

typedef struct MePredUnit {
    uint8_t         ref_index[MAX_NUM_OF_REF_PIC_LIST];
    uint8_t         ref0_list;
    uint8_t         ref1_list;
    uint32_t        distortion;
    EbPredDirection prediction_direction;
} MePredUnit;

typedef struct MotionEstimationTierZero {
    MePredUnit pu[MAX_ME_PU_COUNT];
} MotionEstimationTierZero;
typedef struct HmeResults {
    uint8_t  list_i;   // list index of this ref
    uint8_t  ref_i;    // ref list index of this ref
    int16_t  hme_sc_x; // hme search centre x
    int16_t  hme_sc_y; // hme search centre y
    uint64_t hme_sad;  // hme sad
    uint8_t  do_ref;   // to process this ref in ME or not
} HmeResults;
typedef struct MeContext {
    EbDctor dctor;
    // Search region stride
    uint32_t                  interpolated_stride;
    uint32_t                  interpolated_full_stride[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    MotionEstimationTierZero *me_candidate;
    // Intermediate SB-sized buffer to retain the input samples
    uint8_t * sb_buffer;
    uint8_t * sb_buffer_ptr;
    uint32_t  sb_buffer_stride;
    uint8_t * sb_src_ptr;
    uint32_t  sb_src_stride;
    uint8_t * quarter_sb_buffer;
    uint32_t  quarter_sb_buffer_stride;
    uint8_t * sixteenth_sb_buffer;
    uint32_t  sixteenth_sb_buffer_stride;
    uint8_t * integer_buffer_ptr[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    uint8_t * pos_b_buffer[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    uint8_t * pos_h_buffer[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    uint8_t * pos_j_buffer[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    uint8_t * one_d_intermediate_results_buf0;
    uint8_t * one_d_intermediate_results_buf1;
    int16_t   x_search_area_origin[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    int16_t   y_search_area_origin[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    // ME Parameters
    /* Number of search positions in the horizontal direction.
    *
    * Default depends on input resolution. */
    uint32_t sa_width[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    /* Number of search positions in the vertical direction.
    *
    * Default depends on input resolution. */
    uint32_t sa_height[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    uint8_t * avctemp_buffer;
    uint32_t *p_best_sad_8x8;
    uint32_t *p_best_sad_16x16;
    uint32_t *p_best_sad_32x32;
    uint32_t *p_best_sad_64x64;
    uint32_t *p_best_mv8x8;
    uint32_t *p_best_mv16x16;
    uint32_t *p_best_mv32x32;
    uint32_t *p_best_mv64x64;
    uint32_t *p_best_sad_64x32;
    uint32_t *p_best_sad_32x16;
    uint32_t *p_best_sad_16x8;
    uint32_t *p_best_sad_32x64;
    uint32_t *p_best_sad_16x32;
    uint32_t *p_best_sad_8x16;
    uint32_t *p_best_sad_32x8;
    uint32_t *p_best_sad_8x32;
    uint32_t *p_best_sad_64x16;
    uint32_t *p_best_sad_16x64;
    uint32_t *p_best_mv64x32;
    uint32_t *p_best_mv32x16;
    uint32_t *p_best_mv16x8;
    uint32_t *p_best_mv32x64;
    uint32_t *p_best_mv16x32;
    uint32_t *p_best_mv8x16;
    uint32_t *p_best_mv32x8;
    uint32_t *p_best_mv8x32;
    uint32_t *p_best_mv64x16;
    uint32_t *p_best_mv16x64;
    EB_ALIGN(16) uint32_t p_sad32x32[4];
    EB_ALIGN(64) uint32_t p_sad16x16[16];
    EB_ALIGN(64) uint32_t p_sad8x8[64];

    uint8_t psub_pel_direction64x64;
    uint8_t psub_pel_direction32x32[4];
    uint8_t psub_pel_direction16x16[16];
    uint8_t psub_pel_direction8x8[64];
    uint8_t psub_pel_direction64x32[2];
    uint8_t psub_pel_direction32x16[8];
    uint8_t psub_pel_direction16x8[32];
    uint8_t psub_pel_direction32x64[2];
    uint8_t psub_pel_direction16x32[8];
    uint8_t psub_pel_direction8x16[32];
    uint8_t psub_pel_direction32x8[16];
    uint8_t psub_pel_direction8x32[16];
    uint8_t psub_pel_direction64x16[4];
    uint8_t psub_pel_direction16x64[4];

    uint32_t  p_sb_best_sad[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][MAX_ME_PU_COUNT];
    uint32_t  p_sb_best_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][MAX_ME_PU_COUNT];
    uint32_t  p_sb_bipred_sad[MAX_ME_PU_COUNT]; //needs to be upgraded to 209 pus
    uint32_t  p_sb_best_full_pel_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][MAX_ME_PU_COUNT];
    uint32_t *p_best_full_pel_mv8x8;
    uint32_t *p_best_full_pel_mv16x16;
    uint32_t *p_best_full_pel_mv32x32;
    uint32_t *p_best_full_pel_mv64x64;
    uint32_t *p_best_full_pel_mv64x32;
    uint32_t *p_best_full_pel_mv32x16;
    uint32_t *p_best_full_pel_mv16x8;
    uint32_t *p_best_full_pel_mv32x64;
    uint32_t *p_best_full_pel_mv16x32;
    uint32_t *p_best_full_pel_mv8x16;
    uint32_t *p_best_full_pel_mv32x8;
    uint32_t *p_best_full_pel_mv8x32;
    uint32_t *p_best_full_pel_mv64x16;
    uint32_t *p_best_full_pel_mv16x64;
    uint8_t   full_quarter_pel_refinement;
    uint32_t  p_sb_best_ssd[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][MAX_ME_PU_COUNT];
    uint32_t *p_best_ssd8x8;
    uint32_t *p_best_ssd16x16;
    uint32_t *p_best_ssd32x32;
    uint32_t *p_best_ssd64x64;
    uint32_t *p_best_ssd64x32;
    uint32_t *p_best_ssd32x16;
    uint32_t *p_best_ssd16x8;
    uint32_t *p_best_ssd32x64;
    uint32_t *p_best_ssd16x32;
    uint32_t *p_best_ssd8x16;
    uint32_t *p_best_ssd32x8;
    uint32_t *p_best_ssd8x32;
    uint32_t *p_best_ssd64x16;
    uint32_t *p_best_ssd16x64;

    uint8_t   p_sb_best_nsq[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][MAX_ME_PU_COUNT];
    uint8_t * p_best_nsq8x8;
    uint8_t * p_best_nsq16x16;
    uint8_t * p_best_nsq32x32;
    uint8_t * p_best_nsq64x64;
    uint16_t *p_eight_pos_sad16x16;
    EB_ALIGN(64) uint32_t p_eight_sad32x32[4][8];
    EB_ALIGN(64) uint32_t p_eight_sad16x16[16][8];
    EB_ALIGN(64) uint32_t p_eight_sad8x8[64][8];
    EbBitFraction *mvd_bits_array;
    uint64_t       lambda;
    uint8_t        hme_search_type;

    uint8_t fractional_search_method;

    uint8_t fractional_search_model;
    uint8_t hme_search_method;
    uint8_t me_search_method;

    EbBool enable_hme_flag;
    EbBool enable_hme_level0_flag;
    EbBool enable_hme_level1_flag;
    EbBool enable_hme_level2_flag;

    EbBool use_subpel_flag;
    EbBool half_pel_mode;

    EbBool compute_global_motion;

    // ME
    uint16_t search_area_width;
    uint16_t search_area_height;
    uint16_t max_me_search_width;
    uint16_t max_me_search_height;
    uint8_t inherit_rec_mv_from_sq_block;
    uint8_t best_list_idx;
    uint8_t best_ref_idx;
    // HME
    uint16_t number_hme_search_region_in_width;
    uint16_t number_hme_search_region_in_height;
    uint16_t hme_level0_total_search_area_width;
    uint16_t hme_level0_total_search_area_height;
    uint16_t hme_level0_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint16_t hme_level0_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint16_t hme_level1_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint16_t hme_level1_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint16_t hme_level2_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint16_t hme_level2_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint8_t  update_hme_search_center_flag;
    HmeResults hme_results[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint32_t reduce_me_sr_flag[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    EbBool local_hp_mode[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // ------- Context for Alt-Ref ME ------
    uint16_t adj_search_area_width;
    uint16_t adj_search_area_height;
    EbBool   me_alt_ref;
    void *   alt_ref_reference_ptr;
    // tf
    int tf_frame_index;
    int tf_index_center;
    uint8_t h_pel_search_wind;
    signed short tf_16x16_mv_x[16];
    signed short tf_16x16_mv_y[16];
    uint64_t tf_16x16_block_error[16];

    signed short tf_32x32_mv_x[4];
    signed short tf_32x32_mv_y[4];
    uint64_t tf_32x32_block_error[4];

    int tf_32x32_block_split_flag[4];
    // -------
} MeContext;

typedef uint64_t (*EB_ME_DISTORTION_FUNC)(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                          uint32_t ref_stride, uint32_t width, uint32_t height);

extern EbErrorType me_context_ctor(MeContext *object_ptr, uint16_t max_input_luma_width,
                                   uint16_t max_input_luma_height, uint8_t nsq_present,
                                   uint8_t mrp_mode);

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationContext_h
