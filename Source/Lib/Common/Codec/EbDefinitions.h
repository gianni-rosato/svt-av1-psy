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

#ifndef EbDefinitions_h
#define EbDefinitions_h
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "EbSvtAv1.h"
#ifdef _WIN32
#define inline __inline
#elif __GNUC__
#define inline __inline__
#else
#define inline
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MULTI_PASS_PD                1 // Multi-Pass Partitioning Depth (Multi-Pass PD) performs multiple PD stages for the same SB towards 1 final Partitioning Structure. As we go from PDn to PDn + 1, the prediction accuracy of the MD feature(s) increases while the number of block(s) decreases

#define HBD_CLEAN_UP                 1

#define IFS_8BIT_MD                  1

#define COMP_HBD                     1
#define INTERINTRA_HBD               1
#define ATB_HBD                      1
#define ATB_FIX                      1

#define M0_OPT                       1
/* Note: shutting the macro PAL_SUP will not give SS as pcs->palette_mode = 0
   rate estimation is changed for I frame + enabled sc for P (rate estimation
   is a result changed for P frames)
*/
#define PAL_SUP                      1 //Palette prediction support.
#if PAL_SUP
#define PAL_CLASS   1
#endif

#define AUTO_MAX_PARTITION           1 // Shortcut to skip search depths depending on motion estimation info and statistics

#define LESS_RECTANGULAR_CHECK_LEVEL 1 // Shortcut to skip a/b shapes depending on SQ/H/V shape costs

#define INTER_INTRA_CLASS_PRUNING    1

#define FIX_ALTREF                   1 // Address ALTREF mismatch between rtime-m0-test and master: fixed actual_future_pics derivation, shut padding of the central frame, fixed end past frame index prior to window shrinking
#define FIX_NEAREST_NEW              1 // Address NEAREST_NEW mismatch between rtime-m0-test and master: fixed injection and fixed settings
#define FIX_ESTIMATE_INTRA           1 // Address ESTIMATE_INTRA mismatch between rtime-m0-test and master: fixed settings
#define FIX_SKIP_REDUNDANT_BLOCK     1 // Address SKIP_REDUNDANT_BLOCK mismatch between rtime-m0-test and master: fixed the action to bypass MD
#define FIX_COEF_BASED_ATB_SKIP      1 // Address COEF_BASED_ATB_SKIP mismatch between rtime-m0-test and master: reset coeff_based_skip_atb @ each SB
#define FIX_WM_SETTINGS              1 // Address WM_SETTINGS mismatch between rtime-m0-test and master: fixed settings
#define FIX_ENABLE_CDF_UPDATE        1 // Address ENABLE_CDF_UPDATE mismatch between rtime-m0-test and master: removed useless update
#define FIX_SORTING_METHOD           1 // Address SORTING mismatch between rtime-m0-test and master: used same method
#define FIX_SETTINGS_RESET           1 // Address SEGMENT_RESET mismatch between rtime-m0-test and master: only @ 1st segment
#define FIX_COMPOUND                 1 // Address COMPOUND mismatch between rtime-m0-test and master: used block size @ the derivation of compound count

#define OBMC_FLAG            1 // OBMC motion mode flag

#define INJECT_NEW_NEAR_NEAR_NEW   1   // Inject NEW_NEAR / NEAR_NEW inter prediction
#define FILTER_INTRA_FLAG    1 // Filter intra prediction


#define II_COMP_FLAG                 1 // InterIntra compound
#define PAETH_HBD                    1 // Enbale Intra PAETH for 10bit
#define INTER_INTER_HBD              1 // Upgrade InterInter compound 10bit
#define INTER_INTRA_HBD              1 // Upgrade InterIntra compound 10bit
#define ATB_10_BIT                   1 // Upgrade ATB  10bit

#define PRED_CHANGE                  1 // Change the MRP in 4L Pictures 3, 5 , 7 and 9 use 1 as the reference
#define PRED_CHANGE_5L               1 // Change the MRP in 5L Pictures 3, 5 , 7 and 9 use 1 as the reference, 11, 13, 15 and 17 use 9 as the reference
#define PRED_CHANGE_MOD              1 // Reorder the references for MRP
#define SPEED_OPT                    1 // Speed optimization(s)
#define GLOBAL_WARPED_MOTION         1 // Global warped motion detection and insertion
#define GM_OPT                       1 // Perform global motion estimation on a down-sampled version of the input picture

#ifndef NON_AVX512_SUPPORT
#define NON_AVX512_SUPPORT
#endif

#define MR_MODE                           0

#define WARP_UPDATE                       1 // Modified Warp settings: ON for MR mode. ON for ref frames in M0
#define UPDATE_CDEF                       1 // Update bit cost estimation for CDEF filter
#define EIGTH_PEL_MV                      1
#define EIGHT_PEL_PREDICTIVE_ME           1
#define EIGHT_PEL_FIX                     1 // Improve the 8th pel performance by shutting it based on the qindex and bug fixes
#define HIGH_PRECISION_MV_QTHRESH         150
#define COMP_INTERINTRA                   1 // InterIntra mode support

#define ENHANCE_ATB                       1

#define RDOQ_CHROMA                       1


#define TWO_PASS                          1 // Two pass encoding. For now, the encoder is called two times and data transfered using file.
                                            // Actions in the second pass: Frame and SB QP assignment and temporal filtering strenght change
#define TWO_PASS_USE_2NDP_ME_IN_1STP      1 // Add a config parameter to the first pass to use the ME settings of the second pass

#define REMOVE_MD_STAGE_1                 1 // Simplified MD Staging; removed md_stage_1
#define NON_KF_INTRA_TF_FIX               0 // Fix temporal filtering for non-key Intra frames

#define TWO_PASS_IMPROVEMENT              1 // Tune 2 pass for better Luma by adjusting the reference area and the actions
//FOR DEBUGGING - Do not remove
#define NO_ENCDEC                         0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch

#define ADP_STATS_PER_LAYER                             0
#define AOM_INTERP_EXTEND                               4
#define OPTIMISED_EX_SUBPEL                             1

#define AOM_LEFT_TOP_MARGIN_PX(subsampling) \
  ((AOM_BORDER_IN_PIXELS >> subsampling) - AOM_INTERP_EXTEND)
#define AOM_LEFT_TOP_MARGIN_SCALED(subsampling) \
  (AOM_LEFT_TOP_MARGIN_PX(subsampling) << SCALE_SUBPEL_BITS)

#if OPTIMISED_EX_SUBPEL
#define H_PEL_SEARCH_WIND 3  // 1/2-pel serach window
#else
#define H_PEL_SEARCH_WIND 4  // 1/2-pel serach window
#endif
#define Q_PEL_SEARCH_WIND 2  // 1/4-pel serach window
#define HP_REF_OPT        1  // Remove redundant positions.
typedef enum ME_HP_MODE {
    EX_HP_MODE = 0,       // Exhaustive  1/2-pel serach mode.
    REFINMENT_HP_MODE = 1 // Refinement 1/2-pel serach mode.
} ME_HP_MODE;
typedef enum ME_QP_MODE {
    EX_QP_MODE = 0,       // Exhaustive  1/4-pel serach mode.
    REFINMENT_QP_MODE = 1 // Refinement 1/4-pel serach mode.
} ME_QP_MODE;
#if GM_OPT
typedef enum GM_LEVEL {
    GM_FULL         = 0,       // Exhaustive search mode.
    GM_DOWN         = 1,       // Downsampled search mode, with a downsampling factor of 2 in each dimension
    GM_TRAN_ONLY    = 2        // Translation only using ME MV.
} GM_LEVEL;
#endif
struct Buf2D
{
    uint8_t *buf;
    uint8_t *buf0;
    int width;
    int height;
    int stride;
};
typedef struct MvLimits
{
    int col_min;
    int col_max;
    int row_min;
    int row_max;
} MvLimits;
typedef struct {
    uint8_t by;
    uint8_t bx;
    uint8_t skip;
} cdef_list;

/*!\brief force enum to be unsigned 1 byte*/
#define UENUM1BYTE(enumvar) \
  ;                         \
  typedef uint8_t enumvar

enum {
    DIAMOND = 0,
    NSTEP = 1,
    HEX = 2,
    BIGDIA = 3,
    SQUARE = 4,
    FAST_HEX = 5,
    FAST_DIAMOND = 6
} UENUM1BYTE(SEARCH_METHODS);

/********************************************************/
/****************** Pre-defined Values ******************/
/********************************************************/

#define ALTREF_MAX_NFRAMES 10 // maximum number of frames allowed for the Alt-ref picture computation
                              // this number can be increased by increasing the constant
                              // FUTURE_WINDOW_WIDTH defined in EbPictureDecisionProcess.c
#define ALTREF_MAX_STRENGTH 6

#define PAD_VALUE                                (128+32)

/* Use open-loop data to predict the NSQ partitions. */
#if MULTI_PASS_PD
#define PREDICT_NSQ_SHAPE                               0 // to set to 0
#else
#define PREDICT_NSQ_SHAPE                               1
#endif
#if PREDICT_NSQ_SHAPE
#define NUMBER_OF_DEPTH                                 6
#define NUMBER_OF_SHAPES                                10
#define ADD_SAD_FOR_128X128                             1
#define ADJUST_NSQ_RANK_BASED_ON_NEIGH                  1
#define COMBINE_MDC_NSQ_TABLE                           1
#define ADD_SUPPORT_TO_SKIP_PART_N                      1
#define ADD_MDC_REFINEMENT_LOOP                         1
#define ADD_MDC_FULL_COST                               1
#define NSQ_TAB_SIZE                                    8
#define MAX_MDC_LEVEL                                   8
#define MDC_ADAPTIVE_LEVEL                              1
#else
#if MULTI_PASS_PD
#define NSQ_TAB_SIZE                                    8
#define NUMBER_OF_DEPTH                                 6
#define NUMBER_OF_SHAPES                                10
#else
#define NSQ_TAB_SIZE                                    6
#endif
#endif

#define AUTO_MODE                                 -1

//  Delta QP support
#define ADD_DELTA_QP_SUPPORT                      1  // Add delta QP support
#define BLOCK_MAX_COUNT_SB_128                    4421  // TODO: reduce alloction for 64x64
#define BLOCK_MAX_COUNT_SB_64                     1101  // TODO: reduce alloction for 64x64
#define MAX_TXB_COUNT                             4 // Maximum number of transform blocks.
#if II_COMP_FLAG
#if OBMC_FLAG
#if FILTER_INTRA_FLAG
#define MAX_NFL                                 125 // Maximum number of candidates MD can support
#else
#define MAX_NFL                                 120 // Maximum number of candidates MD can support
#endif
#else
#define MAX_NFL                                  95
#endif
#else
#define MAX_NFL                                   80
#endif
#define MAX_NFL_BUFF                              (MAX_NFL + CAND_CLASS_TOTAL)  //need one extra temp buffer for each fast loop call
#define MAX_LAD                                   120 // max lookahead-distance 2x60fps
#define ROUND_UV(x) (((x)>>3)<<3)
#define AV1_PROB_COST_SHIFT 9
#define AOMINNERBORDERINPIXELS 160
#define SWITCHABLE_FILTER_CONTEXTS ((SWITCHABLE_FILTERS + 1) * 4)
#define MAX_MB_PLANE   3
#define CFL_MAX_BlockSize (BLOCK_32X32)
#define CFL_BUF_LINE (32)
#define CFL_BUF_LINE_I128 (CFL_BUF_LINE >> 3)
#define CFL_BUF_LINE_I256 (CFL_BUF_LINE >> 4)
#define CFL_BUF_SQUARE (CFL_BUF_LINE * CFL_BUF_LINE)
/***********************************    AV1_OBU     ********************************/
#define INVALID_NEIGHBOR_DATA 0xFFu
#define CONFIG_BITSTREAM_DEBUG 0
#define CONFIG_BUFFER_MODEL 1
#define CONFIG_COEFFICIENT_RANGE_CHECKING 0
#define CONFIG_ENTROPY_STATS 0
#define CONFIG_FP_MB_STATS 0
#define CONFIG_INTERNAL_STATS 0
#define CONFIG_RD_DEBUG 0

// Max superblock size
#define MAX_SB_SIZE_LOG2 7
#define MAX_SB_SIZE (1 << MAX_SB_SIZE_LOG2)
#define MAX_SB_SQUARE (MAX_SB_SIZE * MAX_SB_SIZE)
#define SB_STRIDE_Y    MAX_SB_SIZE
#define SB_STRIDE_UV  (MAX_SB_SIZE>>1)

// Min superblock size
#define MIN_SB_SIZE_LOG2 6

// Pixels per Mode Info (MI) unit
#define MI_SIZE_LOG2 2
#define MI_SIZE (1 << MI_SIZE_LOG2)

// MI-units per max superblock (MI Block - MIB)
#define MAX_MIB_SIZE_LOG2 (MAX_SB_SIZE_LOG2 - MI_SIZE_LOG2)
#define MAX_MIB_SIZE (1 << MAX_MIB_SIZE_LOG2)

// MI-units per min superblock
#define SB64_MIB_SIZE 16

// MI-units per min superblock
#define MIN_MIB_SIZE_LOG2 (MIN_SB_SIZE_LOG2 - MI_SIZE_LOG2)

// Mask to extract MI offset within max MIB
#define MAX_MIB_MASK (MAX_MIB_SIZE - 1)

// Maximum size of a loop restoration tile
#define RESTORATION_TILESIZE_MAX 256
// Maximum number of tile rows and tile columns
#define MAX_TILE_ROWS 64
#define MAX_TILE_COLS 64
#if ENHANCE_ATB
#define MAX_VARTX_DEPTH 2
#else
#define MAX_VARTX_DEPTH 1
#endif
#define MI_SIZE_64X64 (64 >> MI_SIZE_LOG2)
#define MI_SIZE_128X128 (128 >> MI_SIZE_LOG2)
#define MAX_PALETTE_SQUARE (64 * 64)
// Maximum number of colors in a palette.
#define PALETTE_MAX_SIZE 8
// Minimum number of colors in a palette.
#define PALETTE_MIN_SIZE 2
#define FRAME_OFFSET_BITS 5
#define MAX_FRAME_DISTANCE ((1 << FRAME_OFFSET_BITS) - 1)

// 4 frame filter levels: y plane vertical, y plane horizontal,
// u plane, and v plane
#define FRAME_LF_COUNT 4
#define DEFAULT_DELTA_LF_MULTI 0
#define MAX_MODE_LF_DELTAS 2
#define LEVEL_MAJOR_BITS 3
#define LEVEL_MINOR_BITS 2
#define LEVEL_BITS (LEVEL_MAJOR_BITS + LEVEL_MINOR_BITS)

#define LEVEL_MAJOR_MIN 2
#define LEVEL_MAJOR_MAX ((1 << LEVEL_MAJOR_BITS) - 1 + LEVEL_MAJOR_MIN)
#define LEVEL_MINOR_MIN 0
#define LEVEL_MINOR_MAX ((1 << LEVEL_MINOR_BITS) - 1)

#define OP_POINTS_CNT_MINUS_1_BITS 5
#define OP_POINTS_IDC_BITS 12
#define TX_SIZE_LUMA_MIN (TX_4X4)
/* We don't need to code a transform size unless the allowed size is at least
one more than the minimum. */
#define TX_SIZE_CTX_MIN (TX_SIZE_LUMA_MIN + 1)

// Maximum tx_size categories
#define MAX_TX_CATS (TX_SIZES - TX_SIZE_CTX_MIN)
#define MAX_TX_DEPTH 2

#define MAX_TX_SIZE_LOG2 (6)
#define MAX_TX_SIZE (1 << MAX_TX_SIZE_LOG2)
#define MIN_TX_SIZE_LOG2 2
#define MIN_TX_SIZE (1 << MIN_TX_SIZE_LOG2)
#define MAX_TX_SQUARE (MAX_TX_SIZE * MAX_TX_SIZE)

// Pad 4 extra columns to remove horizontal availability check.
#define TX_PAD_HOR_LOG2 2
#define TX_PAD_HOR 4
// Pad 6 extra rows (2 on top and 4 on bottom) to remove vertical availability
// check.
#define TX_PAD_TOP 2
#define TX_PAD_BOTTOM 4
#define TX_PAD_VER (TX_PAD_TOP + TX_PAD_BOTTOM)
// Pad 16 extra bytes to avoid reading overflow in SIMD optimization.
#define TX_PAD_END 16
#define TX_PAD_2D \
((MAX_TX_SIZE + TX_PAD_HOR) * (MAX_TX_SIZE + TX_PAD_VER) + TX_PAD_END)
#define COMPOUND_WEIGHT_MODE DIST
#define DIST_PRECISION_BITS 4
#define DIST_PRECISION (1 << DIST_PRECISION_BITS)  // 16

// TODO(chengchen): Temporal flag serve as experimental flag for WIP
// bitmask construction.
// Shall be removed when bitmask code is completely checkedin
#define LOOP_FILTER_BITMASK 0
#define PROFILE_BITS 3

// AV1 Loop Filter
#define AV1_LF                                    1  // AV1 Loop Filter
#if AV1_LF
#define LF_SHARPNESS 0
#endif

#define FILTER_BITS 7
#define SUBPEL_BITS 4
#define SUBPEL_MASK ((1 << SUBPEL_BITS) - 1)
#define SUBPEL_SHIFTS (1 << SUBPEL_BITS)
#define SUBPEL_TAPS 8

#define SCALE_SUBPEL_BITS 10
#define SCALE_SUBPEL_SHIFTS (1 << SCALE_SUBPEL_BITS)
#define SCALE_SUBPEL_MASK (SCALE_SUBPEL_SHIFTS - 1)
#define SCALE_EXTRA_BITS (SCALE_SUBPEL_BITS - SUBPEL_BITS)
#define SCALE_EXTRA_OFF ((1 << SCALE_EXTRA_BITS) / 2)

typedef int16_t InterpKernel[SUBPEL_TAPS];

/***************************************************/
/****************** Helper Macros ******************/
/***************************************************/
void aom_reset_mmx_state(void);
extern void RunEmms();
#define aom_clear_system_state() RunEmms() //aom_reset_mmx_state()

/* Shift down with rounding for use when n >= 0, value >= 0 */
#define ROUND_POWER_OF_TWO(value, n) (((value)+(((1 << (n)) >> 1))) >> (n))

/* Shift down with rounding for signed integers, for use when n >= 0 */
#define ROUND_POWER_OF_TWO_SIGNED(value, n)           \
    (((value) < 0) ? -ROUND_POWER_OF_TWO(-(value), (n)) \
                 : ROUND_POWER_OF_TWO((value), (n)))

/* Shift down with rounding for use when n >= 0, value >= 0 for (64 bit) */
#define ROUND_POWER_OF_TWO_64(value, n) \
    (((value) + ((((int64_t)1 << (n)) >> 1))) >> (n))

/* Shift down with rounding for signed integers, for use when n >= 0 (64 bit) */
#define ROUND_POWER_OF_TWO_SIGNED_64(value, n)           \
    (((value) < 0) ? -ROUND_POWER_OF_TWO_64(-(value), (n)) \
                 : ROUND_POWER_OF_TWO_64((value), (n)))

#define IS_POWER_OF_TWO(x) (((x) & ((x)-1)) == 0)

#ifdef __cplusplus
#define EB_EXTERN extern "C"
#else
#define EB_EXTERN
#endif // __cplusplus

#define INLINE __inline
#define RESTRICT
#ifdef _WIN32
#define FOPEN(f,s,m) fopen_s(&f,s,m)
#else
#define FOPEN(f,s,m) f=fopen(s,m)
#endif

#define IMPLIES(a, b) (!(a) || (b))  //  Logical 'a implies b' (or 'a -> b')
#if (defined(__GNUC__) && __GNUC__) || defined(__SUNPRO_C)
#define DECLARE_ALIGNED(n, typ, val) typ val __attribute__((aligned(n)))
#elif defined(_WIN32)
#define DECLARE_ALIGNED(n, typ, val) __declspec(align(n)) typ val
#else
#warning No alignment directives known for this compiler.
#define DECLARE_ALIGNED(n, typ, val) typ val
#endif

#ifdef _MSC_VER
#define EB_ALIGN(n) __declspec(align(n))
#elif defined(__GNUC__)
#define EB_ALIGN(n) __attribute__((__aligned__(n)))
#else
#define EB_ALIGN(n)
#endif

#ifdef _MSC_VER
#define AOM_FORCE_INLINE __forceinline
#define AOM_INLINE __inline
#else
#define AOM_FORCE_INLINE __inline__ __attribute__((always_inline))
    // TODO(jbb): Allow a way to force inline off for older compilers.
#define AOM_INLINE inline
#endif

#define SIMD_INLINE static AOM_FORCE_INLINE

    //*********************************************************************************************************************//
    // mem.h
    /* shift right or left depending on sign of n */
#define RIGHT_SIGNED_SHIFT(value, n) \
((n) < 0 ? ((value) << (-(n))) : ((value) >> (n)))
    //*********************************************************************************************************************//
    // cpmmom.h
    // Only need this for fixed-size arrays, for structs just assign.
#define av1_copy(dest, src)              \
{                                      \
    assert(sizeof(dest) == sizeof(src)); \
    memcpy(dest, src, sizeof(src));      \
}

    // mem_ops.h
#ifndef MAU_T
    /* Minimum Access Unit for this target */
#define MAU_T uint8_t
#endif

#ifndef MEM_VALUE_T
#define MEM_VALUE_T int32_t
#endif

#undef MEM_VALUE_T_SZ_BITS
#define MEM_VALUE_T_SZ_BITS (sizeof(MEM_VALUE_T) << 3)

static __inline void mem_put_le16(void *vmem, MEM_VALUE_T val) {
    MAU_T *mem = (MAU_T *)vmem;

    mem[0] = (MAU_T)((val >> 0) & 0xff);
    mem[1] = (MAU_T)((val >> 8) & 0xff);
}

static __inline void mem_put_le32(void *vmem, MEM_VALUE_T val) {
    MAU_T *mem = (MAU_T *)vmem;

    mem[0] = (MAU_T)((val >> 0) & 0xff);
    mem[1] = (MAU_T)((val >> 8) & 0xff);
    mem[2] = (MAU_T)((val >> 16) & 0xff);
    mem[3] = (MAU_T)((val >> 24) & 0xff);
}
/* clang-format on */
//#endif  // AOM_PORTS_MEM_OPS_H_

typedef uint16_t ConvBufType;

typedef struct ConvolveParams
{
    int32_t ref;
    int32_t do_average;
    ConvBufType *dst;
    int32_t dst_stride;
    int32_t round_0;
    int32_t round_1;
    int32_t plane;
    int32_t is_compound;
    int32_t use_jnt_comp_avg;
    int32_t fwd_offset;
    int32_t bck_offset;
    int32_t use_dist_wtd_comp_avg;
} ConvolveParams;

// texture component type
typedef enum ATTRIBUTE_PACKED
{
    COMPONENT_LUMA = 0,            // luma
    COMPONENT_CHROMA = 1,            // chroma (Cb+Cr)
    COMPONENT_CHROMA_CB = 2,            // chroma Cb
    COMPONENT_CHROMA_CR = 3,            // chroma Cr
    COMPONENT_ALL = 4,            // Y+Cb+Cr
    COMPONENT_NONE = 15
}COMPONENT_TYPE;

static INLINE int32_t clamp(int32_t value, int32_t low, int32_t high) {
    return value < low ? low : (value > high ? high : value);
}

static INLINE int64_t clamp64(int64_t value, int64_t low, int64_t high) {
    return value < low ? low : (value > high ? high : value);
}

static INLINE uint8_t clip_pixel(int32_t val) {
    return (uint8_t)((val > 255) ? 255 : (val < 0) ? 0 : val);
}

static INLINE uint16_t clip_pixel_highbd(int32_t val, int32_t bd) {
    switch (bd) {
    case 8:
    default: return (uint16_t)clamp(val, 0, 255);
    case 10: return (uint16_t)clamp(val, 0, 1023);
    case 12: return (uint16_t)clamp(val, 0, 4095);
    }
}

static INLINE unsigned int negative_to_zero(int value) {
    return value & ~(value >> (sizeof(value) * 8 - 1));
}

static INLINE int av1_num_planes(EbColorConfig   *color_info) {
    return color_info->mono_chrome ? 1 : MAX_MB_PLANE;
}

//*********************************************************************************************************************//
// enums.h
/*!\brief Decorator indicating that given struct/union/enum is packed */
#ifndef ATTRIBUTE_PACKED
#if defined(__GNUC__) && __GNUC__
#define ATTRIBUTE_PACKED __attribute__((packed))
#else
#define ATTRIBUTE_PACKED
#endif
#endif /* ATTRIBUTE_PACKED */

#if MULTI_PASS_PD
typedef enum PD_PASS {
    PD_PASS_0,
    PD_PASS_1,
    PD_PASS_2,
    PD_PASS_TOTAL,
} PD_PASS;
#endif

typedef enum CAND_CLASS {
    CAND_CLASS_0,
    CAND_CLASS_1,
    CAND_CLASS_2,
    CAND_CLASS_3,
#if II_COMP_FLAG
    CAND_CLASS_4,
#endif
#if OBMC_FLAG
    CAND_CLASS_5,
#endif
#if FILTER_INTRA_FLAG
    CAND_CLASS_6,
#endif
#if PAL_CLASS
    CAND_CLASS_7,
#endif
    CAND_CLASS_8,
    CAND_CLASS_TOTAL
} CAND_CLASS;

typedef enum MD_STAGE {
    MD_STAGE_0,
    MD_STAGE_1,
    MD_STAGE_2,
    MD_STAGE_3,
    MD_STAGE_TOTAL
} MD_STAGE;
#if REMOVE_MD_STAGE_1
#define MD_STAGING_MODE_0    0
#define MD_STAGING_MODE_1    1
#else
#define MD_STAGING_MODE_0    0
#define MD_STAGING_MODE_1    1
#define MD_STAGING_MODE_2    2
#define MD_STAGING_MODE_3    3
#endif
#define INTRA_NFL           16
#define INTER_NEW_NFL       16
#define INTER_PRED_NFL      16


#define BEST_CANDIDATE_COUNT 4
#define MAX_REF_TYPE_CAND   30
#define PRUNE_REC_TH         5
#define PRUNE_REF_ME_TH      2
#if !SPEED_OPT
#define MD_EXIT_THSL         0 // MD_EXIT_THSL -->0 is lossless 100 is maximum. Increase with a step of 10-20.
#endif
typedef enum
{
    EIGHTTAP_REGULAR,
    EIGHTTAP_SMOOTH,
    MULTITAP_SHARP,
    BILINEAR,
    INTERP_FILTERS_ALL,
    SWITCHABLE_FILTERS = BILINEAR,
    SWITCHABLE = SWITCHABLE_FILTERS + 1, /* the last switchable one */
    EXTRA_FILTERS = INTERP_FILTERS_ALL - SWITCHABLE_FILTERS,
}InterpFilter;
#if OBMC_FLAG

#define AV1_COMMON Av1Common
enum {
  USE_2_TAPS_ORIG = 0,  // This is used in temporal filtering.
  USE_2_TAPS,
  USE_4_TAPS,
  USE_8_TAPS,
} UENUM1BYTE(SUBPEL_SEARCH_TYPE);
#endif
typedef struct InterpFilterParams
{
    const int16_t *filter_ptr;
    uint16_t taps;
    uint16_t subpel_shifts;
    InterpFilter interp_filter;
} InterpFilterParams;

typedef enum TxSearchLevel
{
    TX_SEARCH_OFF,
    TX_SEARCH_ENC_DEC,
    TX_SEARCH_INTER_DEPTH,
    TX_SEARCH_FULL_LOOP
} TxSearchLevel;

typedef enum InterpolationSearchLevel
{
#if MULTI_PASS_PD
    IT_SEARCH_OFF,
    IT_SEARCH_FAST_LOOP_UV_BLIND,
    IT_SEARCH_FAST_LOOP,
#else
    IT_SEARCH_OFF,
    IT_SEARCH_INTER_DEPTH,
    IT_SEARCH_FULL_LOOP,
    IT_SEARCH_FAST_LOOP_UV_BLIND,
    IT_SEARCH_FAST_LOOP,
#endif
} InterpolationSearchLevel;

typedef enum NsqSearchLevel
{
    NSQ_SEARCH_OFF,
    NSQ_SEARCH_LEVEL1,
    NSQ_SEARCH_LEVEL2,
    NSQ_SEARCH_LEVEL3,
    NSQ_SEARCH_LEVEL4,
    NSQ_SEARCH_LEVEL5,
    NSQ_SEARCH_LEVEL6,
#if PREDICT_NSQ_SHAPE
    NSQ_SEARCH_LEVEL7,
#endif
    NSQ_SEARCH_FULL
} NsqSearchLevel;

#define MAX_PARENT_SQ     6
typedef enum CompoundDistWeightMode {
    DIST,
} CompoundDistWeightMode;

// Profile 0.  8-bit and 10-bit 4:2:0 and 4:0:0 only.
// Profile 1.  8-bit and 10-bit 4:4:4
// Profile 2.  8-bit and 10-bit 4:2:2
//            12 bit  4:0:0, 4:2:2 and 4:4:4
typedef enum BitstreamProfile {
    PROFILE_0,
    PROFILE_1,
    PROFILE_2,
    MAX_PROFILES
} BitstreamProfile;
// Note: Some enums use the attribute 'packed' to use smallest possible integer
// type, so that we can save memory when they are used in structs/arrays.

typedef enum ATTRIBUTE_PACKED {
    BLOCK_4X4,
    BLOCK_4X8,
    BLOCK_8X4,
    BLOCK_8X8,
    BLOCK_8X16,
    BLOCK_16X8,
    BLOCK_16X16,
    BLOCK_16X32,
    BLOCK_32X16,
    BLOCK_32X32,
    BLOCK_32X64,
    BLOCK_64X32,
    BLOCK_64X64,
    BLOCK_64X128,
    BLOCK_128X64,
    BLOCK_128X128,
    BLOCK_4X16,
    BLOCK_16X4,
    BLOCK_8X32,
    BLOCK_32X8,
    BLOCK_16X64,
    BLOCK_64X16,
    BlockSizeS_ALL,
    BlockSizeS = BLOCK_4X16,
    BLOCK_INVALID = 255,
    BLOCK_LARGEST = (BlockSizeS - 1)
} BlockSize;

typedef enum ATTRIBUTE_PACKED {
    PARTITION_NONE,
    PARTITION_HORZ,
    PARTITION_VERT,
    PARTITION_SPLIT,
    PARTITION_HORZ_A,  // HORZ split and the top partition is split again
    PARTITION_HORZ_B,  // HORZ split and the bottom partition is split again
    PARTITION_VERT_A,  // VERT split and the left partition is split again
    PARTITION_VERT_B,  // VERT split and the right partition is split again
    PARTITION_HORZ_4,  // 4:1 horizontal partition
    PARTITION_VERT_4,  // 4:1 vertical partition
    EXT_PARTITION_TYPES,
    PARTITION_TYPES = PARTITION_SPLIT + 1,
    PARTITION_INVALID = 255
} PartitionType;

#define MAX_NUM_BLOCKS_ALLOC  7493  //max number of blocks assuming 128x128-4x4 all partitions.

typedef enum ATTRIBUTE_PACKED {
    PART_N,
    PART_H,
    PART_V,
    PART_HA,
    PART_HB,
    PART_VA,
    PART_VB,
    PART_H4,
    PART_V4,
    PART_S
} PART;

static const uint8_t mi_size_wide[BlockSizeS_ALL] = {
    1, 1, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 1, 4, 2, 8, 4, 16
};
static const uint8_t mi_size_high[BlockSizeS_ALL] = {
    1, 2, 1, 2, 4, 2, 4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 4, 1, 8, 2, 16, 4
};

typedef char PartitionContextType;
#define PARTITION_PLOFFSET 4  // number of probability models per block size
#define PARTITION_BlockSizeS 5
#define PARTITION_CONTEXTS (PARTITION_BlockSizeS * PARTITION_PLOFFSET)

// block transform size
#ifdef _MSC_VER
typedef uint8_t TxSize;
enum ATTRIBUTE_PACKED {
#else
typedef enum ATTRIBUTE_PACKED {
#endif
    TX_4X4,             // 4x4 transform
    TX_8X8,             // 8x8 transform
    TX_16X16,           // 16x16 transform
    TX_32X32,           // 32x32 transform
    TX_64X64,           // 64x64 transform
    TX_4X8,             // 4x8 transform
    TX_8X4,             // 8x4 transform
    TX_8X16,            // 8x16 transform
    TX_16X8,            // 16x8 transform
    TX_16X32,           // 16x32 transform
    TX_32X16,           // 32x16 transform
    TX_32X64,           // 32x64 transform
    TX_64X32,           // 64x32 transform
    TX_4X16,            // 4x16 transform
    TX_16X4,            // 16x4 transform
    TX_8X32,            // 8x32 transform
    TX_32X8,            // 32x8 transform
    TX_16X64,           // 16x64 transform
    TX_64X16,           // 64x16 transform
    TX_SIZES_ALL,       // Includes rectangular transforms
    TX_SIZES = TX_4X8,  // Does NOT include rectangular transforms
    TX_SIZES_LARGEST = TX_64X64,
    TX_INVALID = 255  // Invalid transform size

#ifdef _MSC_VER
};
#else
} TxSize;
#endif
static const TxSize tx_depth_to_tx_size[3][BlockSizeS_ALL] = {
    // tx_depth 0
    {
        TX_4X4,
        TX_4X8,
        TX_8X4,
        TX_8X8,
        TX_8X16,
        TX_16X8,
        TX_16X16,
        TX_16X32,
        TX_32X16,
        TX_32X32,
        TX_32X64,
        TX_64X32,
        TX_64X64,
        TX_64X64,//TX_64X128,
        TX_64X64,//TX_128X64,
        TX_64X64,//TX_128X128,
        TX_4X16,
        TX_16X4,
        TX_8X32,
        TX_32X8,
        TX_16X64,
        TX_64X16
    },
    // tx_depth 1:
    {
        TX_4X4,
        TX_4X8,
        TX_8X4,
        TX_4X4,
        TX_8X8,
        TX_8X8,
        TX_8X8,
        TX_16X16,
        TX_16X16,
        TX_16X16,
        TX_32X32,
        TX_32X32,
        TX_32X32,
        TX_64X64,//TX_64X128,
        TX_64X64,//TX_128X64,
        TX_64X64,//TX_128X128,
        TX_4X4,
        TX_4X4,
        TX_8X8,
        TX_8X8,
        TX_16X16,
        TX_16X16
    },
    // tx_depth 2
    {
        TX_4X4,
        TX_4X8,
        TX_8X4,
        TX_8X8,
        TX_4X4,
        TX_4X4,
        TX_4X4,
        TX_8X8,
        TX_8X8,
        TX_8X8,
        TX_16X16,
        TX_16X16,
        TX_16X16,
        TX_64X64,//TX_64X128,
        TX_64X64,//TX_128X64,
        TX_64X64,//TX_128X128,
        TX_4X16, // No depth 2
        TX_16X4, // No depth 2
        TX_4X4,
        TX_4X4,
        TX_8X8,
        TX_8X8
    }
};
static const int32_t tx_size_wide[TX_SIZES_ALL] = {
    4, 8, 16, 32, 64, 4, 8, 8, 16, 16, 32, 32, 64, 4, 16, 8, 32, 16, 64,
};
// Transform block height in pixels
static const int32_t tx_size_high[TX_SIZES_ALL] = {
    4, 8, 16, 32, 64, 8, 4, 16, 8, 32, 16, 64, 32, 16, 4, 32, 8, 64, 16,
};

 // TranLow  is the datatype used for final transform coefficients.
typedef int32_t TranLow;
typedef uint8_t QmVal;

typedef enum TxClass
{
    TX_CLASS_2D = 0,
    TX_CLASS_HORIZ = 1,
    TX_CLASS_VERT = 2,
    TX_CLASSES = 3,
} TxClass;

static INLINE TxSize av1_get_adjusted_tx_size(TxSize tx_size)
{
    switch (tx_size) {
    case TX_64X64:
    case TX_64X32:
    case TX_32X64: return TX_32X32;
    case TX_64X16: return TX_32X16;
    case TX_16X64: return TX_16X32;
    default: return tx_size;
    }
}

// Transform block width in log2
static const int32_t tx_size_wide_log2[TX_SIZES_ALL] = {
    2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 5, 6, 2, 4, 3, 5, 4, 6,
};

// Transform block height in log2
static const int32_t tx_size_high_log2[TX_SIZES_ALL] = {
    2, 3, 4, 5, 6, 3, 2, 4, 3, 5, 4, 6, 5, 4, 2, 5, 3, 6, 4,
};
#define ALIGN_POWER_OF_TWO(value, n) \
(((value) + ((1 << (n)) - 1)) & ~((1 << (n)) - 1))
#define AOM_PLANE_Y 0       /**< Y (Luminance) plane */
#define AOM_PLANE_U 1       /**< U (Chroma) plane */
#define AOM_PLANE_V 2       /**< V (Chroma) plane */

#define CONVERT_TO_SHORTPTR(x) ((uint16_t *)(((uintptr_t)(x)) << 1))
#define CONVERT_TO_BYTEPTR(x) ((uint8_t *)(((uintptr_t)(x)) >> 1))

#define AOMMIN(x, y) (((x) < (y)) ? (x) : (y))
#define AOMMAX(x, y) (((x) > (y)) ? (x) : (y))

// frame transform mode
typedef enum ATTRIBUTE_PACKED
{
    ONLY_4X4,         // use only 4x4 transform
    TX_MODE_LARGEST,  // transform size is the largest possible for pu size
    TX_MODE_SELECT,   // transform specified for each block
    TX_MODES,
} TxMode;

// 1D tx types
typedef enum ATTRIBUTE_PACKED
{
    DCT_1D,
    ADST_1D,
    FLIPADST_1D,
    IDTX_1D,
    // TODO(sarahparker) need to eventually put something here for the
    // mrc experiment to make this work with the ext-tx pruning functions
    TX_TYPES_1D,
} TxType1D;

typedef enum ATTRIBUTE_PACKED
{
    DCT_DCT,    // DCT  in both horizontal and vertical
    ADST_DCT,   // ADST in vertical, DCT in horizontal
    DCT_ADST,   // DCT  in vertical, ADST in horizontal
    ADST_ADST,  // ADST in both directions
    FLIPADST_DCT,
    DCT_FLIPADST,
    FLIPADST_FLIPADST,
    ADST_FLIPADST,
    FLIPADST_ADST,
    IDTX,
    V_DCT,
    H_DCT,
    V_ADST,
    H_ADST,
    V_FLIPADST,
    H_FLIPADST,
    TX_TYPES,
} TxType;

typedef enum ATTRIBUTE_PACKED
{
    // DCT only
    EXT_TX_SET_DCTONLY,
    // DCT + Identity only
    EXT_TX_SET_DCT_IDTX,
    // Discrete Trig transforms w/o flip (4) + Identity (1)
    EXT_TX_SET_DTT4_IDTX,
    // Discrete Trig transforms w/o flip (4) + Identity (1) + 1D Hor/vert DCT (2)
    EXT_TX_SET_DTT4_IDTX_1DDCT,
    // Discrete Trig transforms w/ flip (9) + Identity (1) + 1D Hor/Ver DCT (2)
    EXT_TX_SET_DTT9_IDTX_1DDCT,
    // Discrete Trig transforms w/ flip (9) + Identity (1) + 1D Hor/Ver (6)
    EXT_TX_SET_ALL16,
    EXT_TX_SET_TYPES
} TxSetType;

typedef struct TxfmParam
{
    // for both forward and inverse transforms
    TxType tx_type;
    TxSize tx_size;
    int32_t lossless;
    int32_t bd;
    // are the pixel buffers octets or shorts?  This should collapse to
    // bd==8 implies !is_hbd, but that's not certain right now.
    int32_t is_hbd;
    TxSetType tx_set_type;
    // for inverse transforms only
    int32_t eob;
} TxfmParam;

#define IS_2D_TRANSFORM(tx_type) (tx_type < IDTX)
#define EXT_TX_SIZES 4       // number of sizes that use extended transforms
#define EXT_TX_SETS_INTER 4  // Sets of transform selections for INTER
#define EXT_TX_SETS_INTRA 3  // Sets of transform selections for INTRA

typedef enum ATTRIBUTE_PACKED
{
    UNIDIR_COMP_REFERENCE,
    BIDIR_COMP_REFERENCE,
    COMP_REFERENCE_TYPES,
} CompReferenceType;

typedef enum ATTRIBUTE_PACKED
{
    PLANE_TYPE_Y,
    PLANE_TYPE_UV,
    PLANE_TYPES
} PlaneType;

#define CFL_ALPHABET_SIZE_LOG2 4
#define CFL_ALPHABET_SIZE (1 << CFL_ALPHABET_SIZE_LOG2)
#define CFL_MAGS_SIZE ((2 << CFL_ALPHABET_SIZE_LOG2) + 1)
#define CFL_IDX_U(idx) (idx >> CFL_ALPHABET_SIZE_LOG2)
#define CFL_IDX_V(idx) (idx & (CFL_ALPHABET_SIZE - 1))

typedef enum ATTRIBUTE_PACKED
{
    CFL_PRED_U,
    CFL_PRED_V,
    CFL_PRED_PLANES
} CflPredType;

typedef enum ATTRIBUTE_PACKED
{
    CFL_SIGN_ZERO,
    CFL_SIGN_NEG,
    CFL_SIGN_POS,
    CFL_SIGNS
} CflSignType;

typedef enum ATTRIBUTE_PACKED
{
    CFL_DISALLOWED,
    CFL_ALLOWED,
    CFL_ALLOWED_TYPES
} CflAllowedType;

// CFL_SIGN_ZERO,CFL_SIGN_ZERO is invalid
#define CFL_JOINT_SIGNS (CFL_SIGNS * CFL_SIGNS - 1)
// CFL_SIGN_U is equivalent to (js + 1) / 3 for js in 0 to 8
#define CFL_SIGN_U(js) (((js + 1) * 11) >> 5)
// CFL_SIGN_V is equivalent to (js + 1) % 3 for js in 0 to 8
#define CFL_SIGN_V(js) ((js + 1) - CFL_SIGNS * CFL_SIGN_U(js))

// There is no context when the alpha for a given plane is zero.
// So there are 2 fewer contexts than joint signs.
#define CFL_ALPHA_CONTEXTS (CFL_JOINT_SIGNS + 1 - CFL_SIGNS)
#define CFL_CONTEXT_U(js) (js + 1 - CFL_SIGNS)
// Also, the contexts are symmetric under swapping the planes.
#define CFL_CONTEXT_V(js) \
(CFL_SIGN_V(js) * CFL_SIGNS + CFL_SIGN_U(js) - CFL_SIGNS)

typedef enum ATTRIBUTE_PACKED {
    PALETTE_MAP,
    COLOR_MAP_TYPES,
} COLOR_MAP_TYPE;

typedef enum ATTRIBUTE_PACKED {
    TWO_COLORS,
    THREE_COLORS,
    FOUR_COLORS,
    FIVE_COLORS,
    SIX_COLORS,
    SEVEN_COLORS,
    EIGHT_COLORS,
    PALETTE_SIZES
} PaletteSize;

typedef enum ATTRIBUTE_PACKED
{
    PALETTE_COLOR_ONE,
    PALETTE_COLOR_TWO,
    PALETTE_COLOR_THREE,
    PALETTE_COLOR_FOUR,
    PALETTE_COLOR_FIVE,
    PALETTE_COLOR_SIX,
    PALETTE_COLOR_SEVEN,
    PALETTE_COLOR_EIGHT,
    PALETTE_COLORS
} PaletteColor;

// Note: All directional predictors must be between V_PRED and D67_PRED (both
// inclusive).
typedef enum ATTRIBUTE_PACKED
{
    DC_PRED,        // Average of above and left pixels
    V_PRED,         // Vertical
    H_PRED,         // Horizontal
    D45_PRED,       // Directional 45  degree
    D135_PRED,      // Directional 135 degree
    D113_PRED,      // Directional 113 degree
    D157_PRED,      // Directional 157 degree
    D203_PRED,      // Directional 203 degree
    D67_PRED,       // Directional 67  degree
    SMOOTH_PRED,    // Combination of horizontal and vertical interpolation
    SMOOTH_V_PRED,  // Vertical interpolation
    SMOOTH_H_PRED,  // Horizontal interpolation
    PAETH_PRED,     // Predict from the direction of smallest gradient
    NEARESTMV,
    NEARMV,
    GLOBALMV,
    NEWMV,
    // Compound ref compound modes
    NEAREST_NEARESTMV,
    NEAR_NEARMV,
    NEAREST_NEWMV,
    NEW_NEARESTMV,
    NEAR_NEWMV,
    NEW_NEARMV,
    GLOBAL_GLOBALMV,
    NEW_NEWMV,
    MB_MODE_COUNT,
    INTRA_MODE_START = DC_PRED,
    INTRA_MODE_END = NEARESTMV,
    INTRA_MODE_NUM = INTRA_MODE_END - INTRA_MODE_START,
    SINGLE_INTER_MODE_START = NEARESTMV,
    SINGLE_INTER_MODE_END = NEAREST_NEARESTMV,
    SINGLE_INTER_MODE_NUM = SINGLE_INTER_MODE_END - SINGLE_INTER_MODE_START,
    COMP_INTER_MODE_START = NEAREST_NEARESTMV,
    COMP_INTER_MODE_END = MB_MODE_COUNT,
    COMP_INTER_MODE_NUM = COMP_INTER_MODE_END - COMP_INTER_MODE_START,
    INTRA_MODES = PAETH_PRED + 1,  // PAETH_PRED has to be the last intra mode.
    INTRA_INVALID = MB_MODE_COUNT,  // For uv_mode in inter blocks
    INTRA_MODE_4x4
} PredictionMode;

// TODO(ltrudeau) Do we really want to pack this?
// TODO(ltrudeau) Do we match with PredictionMode?
typedef enum ATTRIBUTE_PACKED
{
    UV_DC_PRED,        // Average of above and left pixels
    UV_V_PRED,         // Vertical
    UV_H_PRED,         // Horizontal
    UV_D45_PRED,       // Directional 45  degree
    UV_D135_PRED,      // Directional 135 degree
    UV_D113_PRED,      // Directional 113 degree
    UV_D157_PRED,      // Directional 157 degree
    UV_D203_PRED,      // Directional 203 degree
    UV_D67_PRED,       // Directional 67  degree
    UV_SMOOTH_PRED,    // Combination of horizontal and vertical interpolation
    UV_SMOOTH_V_PRED,  // Vertical interpolation
    UV_SMOOTH_H_PRED,  // Horizontal interpolation
    UV_PAETH_PRED,     // Predict from the direction of smallest gradient
    UV_CFL_PRED,       // Chroma-from-Luma
    UV_INTRA_MODES,
    UV_MODE_INVALID,  // For uv_mode in inter blocks
} UvPredictionMode;

typedef enum ATTRIBUTE_PACKED
{
    SIMPLE_TRANSLATION,
    OBMC_CAUSAL,    // 2-sided OBMC
    WARPED_CAUSAL,  // 2-sided WARPED
    MOTION_MODES
} MotionMode;

typedef enum ATTRIBUTE_PACKED
{
    II_DC_PRED,
    II_V_PRED,
    II_H_PRED,
    II_SMOOTH_PRED,
    INTERINTRA_MODES
} InterIntraMode;

typedef enum ATTRIBUTE_PACKED
{
    COMPOUND_AVERAGE,
    COMPOUND_DISTWTD,
    COMPOUND_WEDGE,
    COMPOUND_DIFFWTD,
    COMPOUND_TYPES,
    MASKED_COMPOUND_TYPES = 2,
} CompoundType;

#define   COMPOUND_INTRA  4//just for the decoder
#define AOM_BLEND_A64_ROUND_BITS 6
#define AOM_BLEND_A64_MAX_ALPHA (1 << AOM_BLEND_A64_ROUND_BITS)  // 64

#define AOM_BLEND_A64(a, v0, v1)                                          \
  ROUND_POWER_OF_TWO((a) * (v0) + (AOM_BLEND_A64_MAX_ALPHA - (a)) * (v1), \
                     AOM_BLEND_A64_ROUND_BITS)
#define DIFF_FACTOR_LOG2 4
#define DIFF_FACTOR (1 << DIFF_FACTOR_LOG2)
#define AOM_BLEND_AVG(v0, v1) ROUND_POWER_OF_TWO((v0) + (v1), 1)
typedef uint16_t CONV_BUF_TYPE;
#define MAX_WEDGE_TYPES (1 << 4)
#define MAX_WEDGE_SIZE_LOG2 5  // 32x32
#define MAX_WEDGE_SIZE (1 << MAX_WEDGE_SIZE_LOG2)
#define MAX_WEDGE_SQUARE (MAX_WEDGE_SIZE * MAX_WEDGE_SIZE)
#define WEDGE_WEIGHT_BITS 6
#define WEDGE_NONE -1
#define MASK_MASTER_SIZE ((MAX_WEDGE_SIZE) << 1)
#define MASK_MASTER_STRIDE (MASK_MASTER_SIZE)
typedef struct {
    int enable_order_hint;           // 0 - disable order hint, and related tools
    int order_hint_bits_minus_1;     // dist_wtd_comp, ref_frame_mvs,
    int enable_dist_wtd_comp;        // 0 - disable dist-wtd compound modes
    int enable_ref_frame_mvs;        // 0 - disable ref frame mvs
} OrderHintInfoEnc;
enum {
    MD_COMP_AVG,
    MD_COMP_DIST,
    MD_COMP_DIFF0,
    MD_COMP_WEDGE,
    MD_COMP_TYPES,
} UENUM1BYTE(MD_COMP_TYPE);
#define COMPOUND_TYPE  CompoundType
#define MAX_DIFFWTD_MASK_BITS 1
enum {
    DIFFWTD_38 = 0,
    DIFFWTD_38_INV,
    DIFFWTD_MASK_TYPES,
} UENUM1BYTE(DIFFWTD_MASK_TYPE);
typedef struct {

    /*!< Specifies how the two predictions should be blended together. */
    CompoundType type;

    /*!< Used to derive the direction and offset of the wedge mask used during blending. */
    uint8_t wedge_index;

    /*!< Specifies the sign of the wedge blend. */
    uint8_t wedge_sign;

    /*!< Specifies the type of mask to be used during blending. */
    DIFFWTD_MASK_TYPE mask_type;
} InterInterCompoundData;

#if II_COMP_FLAG
#define INTERINTRA_MODE  InterIntraMode
#endif
typedef enum ATTRIBUTE_PACKED
{
    FILTER_DC_PRED,
    FILTER_V_PRED,
    FILTER_H_PRED,
    FILTER_D157_PRED,
    FILTER_PAETH_PRED,
    FILTER_INTRA_MODES,
} FilterIntraMode;

#if FILTER_INTRA_FLAG
static const PredictionMode fimode_to_intramode[FILTER_INTRA_MODES] = {
  DC_PRED, V_PRED, H_PRED, D157_PRED, PAETH_PRED
};
#endif
#define DIRECTIONAL_MODES 8
#define MAX_ANGLE_DELTA 3
#define ANGLE_STEP 3

#define INTER_MODES (1 + NEWMV - NEARESTMV)

#define INTER_COMPOUND_MODES (1 + NEW_NEWMV - NEAREST_NEARESTMV)

#define SKIP_CONTEXTS 3
#define SKIP_MODE_CONTEXTS 3

#define COMP_INDEX_CONTEXTS 6
#define COMP_GROUP_IDX_CONTEXTS 6

#define NMV_CONTEXTS 3

#define NEWMV_MODE_CONTEXTS 6
#define GLOBALMV_MODE_CONTEXTS 2
#define REFMV_MODE_CONTEXTS 6
#define DRL_MODE_CONTEXTS 3

#define GLOBALMV_OFFSET 3
#define REFMV_OFFSET 4

#define NEWMV_CTX_MASK ((1 << GLOBALMV_OFFSET) - 1)
#define GLOBALMV_CTX_MASK ((1 << (REFMV_OFFSET - GLOBALMV_OFFSET)) - 1)
#define REFMV_CTX_MASK ((1 << (8 - REFMV_OFFSET)) - 1)

#define COMP_NEWMV_CTXS 5
#define INTER_MODE_CONTEXTS 8

#define DELTA_Q_SMALL 3
#define DELTA_Q_PROBS (DELTA_Q_SMALL)
#define DEFAULT_DELTA_Q_RES 1
#define DELTA_LF_SMALL 3
#define DELTA_LF_PROBS (DELTA_LF_SMALL)
#define DEFAULT_DELTA_LF_RES 2

/* Segment Feature Masks */
#define MAX_MV_REF_CANDIDATES 2

#define MAX_REF_MV_STACK_SIZE 8
#define REF_CAT_LEVEL 640

#define INTRA_INTER_CONTEXTS 4
#define COMP_INTER_CONTEXTS 5
#define REF_CONTEXTS 3

#define COMP_REF_TYPE_CONTEXTS 5
#define UNI_COMP_REF_CONTEXTS 3

#define TXFM_PARTITION_CONTEXTS ((TX_SIZES - TX_8X8) * 6 - 3)
typedef uint8_t TXFM_CONTEXT;

// frame types
#define NONE_FRAME -1
#define INTRA_FRAME 0
#define LAST_FRAME 1
#define LAST2_FRAME 2
#define LAST3_FRAME 3
#define GOLDEN_FRAME 4
#define BWDREF_FRAME 5
#define ALTREF2_FRAME 6
#define ALTREF_FRAME 7
#define LAST_REF_FRAMES (LAST3_FRAME - LAST_FRAME + 1)

#define REFS_PER_FRAME 7
#define INTER_REFS_PER_FRAME (ALTREF_FRAME - LAST_FRAME + 1)
#define TOTAL_REFS_PER_FRAME (ALTREF_FRAME - INTRA_FRAME + 1)

#define FWD_REFS (GOLDEN_FRAME - LAST_FRAME + 1)
#define FWD_RF_OFFSET(ref) (ref - LAST_FRAME)
#define BWD_REFS (ALTREF_FRAME - BWDREF_FRAME + 1)
#define BWD_RF_OFFSET(ref) (ref - BWDREF_FRAME)

#define SINGLE_REFS (FWD_REFS + BWD_REFS)

typedef enum ATTRIBUTE_PACKED
{
    LAST_LAST2_FRAMES,      // { LAST_FRAME, LAST2_FRAME }
    LAST_LAST3_FRAMES,      // { LAST_FRAME, LAST3_FRAME }
    LAST_GOLDEN_FRAMES,     // { LAST_FRAME, GOLDEN_FRAME }
    BWDREF_ALTREF_FRAMES,   // { BWDREF_FRAME, ALTREF_FRAME }
    LAST2_LAST3_FRAMES,     // { LAST2_FRAME, LAST3_FRAME }
    LAST2_GOLDEN_FRAMES,    // { LAST2_FRAME, GOLDEN_FRAME }
    LAST3_GOLDEN_FRAMES,    // { LAST3_FRAME, GOLDEN_FRAME }
    BWDREF_ALTREF2_FRAMES,  // { BWDREF_FRAME, ALTREF2_FRAME }
    ALTREF2_ALTREF_FRAMES,  // { ALTREF2_FRAME, ALTREF_FRAME }
    TOTAL_UNIDIR_COMP_REFS,
    // NOTE: UNIDIR_COMP_REFS is the number of uni-directional reference pairs
    //       that are explicitly signaled.
    UNIDIR_COMP_REFS = BWDREF_ALTREF_FRAMES + 1,
} UniDirCompRef;

#define TOTAL_COMP_REFS (FWD_REFS * BWD_REFS + TOTAL_UNIDIR_COMP_REFS)

#define COMP_REFS (FWD_REFS * BWD_REFS + UNIDIR_COMP_REFS)

// NOTE: A limited number of unidirectional reference pairs can be signalled for
//       compound prediction. The use of skip mode, on the other hand, makes it
//       possible to have a reference pair not listed for explicit signaling.
#define MODE_CTX_REF_FRAMES (TOTAL_REFS_PER_FRAME + TOTAL_COMP_REFS)

typedef enum ATTRIBUTE_PACKED
{
    RESTORE_NONE,
    RESTORE_WIENER,
    RESTORE_SGRPROJ,
    RESTORE_SWITCHABLE,
    RESTORE_SWITCHABLE_TYPES = RESTORE_SWITCHABLE,
    RESTORE_TYPES = 4,
} RestorationType;

#define SCALE_NUMERATOR 8
#define SUPERRES_SCALE_BITS 3
#define SUPERRES_SCALE_DENOMINATOR_MIN (SCALE_NUMERATOR + 1)

//*********************************************************************************************************************//
// assert.h
#undef assert

#ifdef NDEBUG

#define assert(expression) ((void)0)

#else
#define assert(expression) ((void)0)

#endif
//**********************************************************************************************************************//
// onyxc_int.h
#define CDEF_MAX_STRENGTHS 16

#define REF_FRAMES_LOG2 3
#define REF_FRAMES (1 << REF_FRAMES_LOG2)

// 4 scratch frames for the new frames to support a maximum of 4 cores decoding
// in parallel, 3 for scaled references on the encoder.
// TODO(hkuang): Add ondemand frame buffers instead of hardcoding the number
// of framebuffers.
// TODO(jkoleszar): These 3 extra references could probably come from the
// normal reference pool.
#define FRAME_BUFFERS (REF_FRAMES + 7)

/* Constant values while waiting for the sequence header */
#define FRAME_ID_LENGTH 15
#define DELTA_FRAME_ID_LENGTH 14

#define FRAME_CONTEXTS (FRAME_BUFFERS + 1)
// Extra frame context which is always kept at default values
#define FRAME_CONTEXT_DEFAULTS (FRAME_CONTEXTS - 1)
#define PRIMARY_REF_BITS 3
#define PRIMARY_REF_NONE 7

#define NUM_PING_PONG_BUFFERS 2

#define MAX_NUM_TEMPORAL_LAYERS 8
#define MAX_NUM_SPATIAL_LAYERS 4
/* clang-format off */
// clang-format seems to think this is a pointer dereference and not a
// multiplication.
#define MAX_NUM_OPERATING_POINTS \
MAX_NUM_TEMPORAL_LAYERS * MAX_NUM_SPATIAL_LAYERS

static INLINE int32_t is_valid_seq_level_idx(uint8_t seq_level_idx) {
    return seq_level_idx < 24 || seq_level_idx == 31;
}
// TODO(jingning): Turning this on to set up transform coefficient
// processing timer.
#define TXCOEFF_TIMER 0
#define TXCOEFF_COST_TIMER 0

typedef enum
{
    SINGLE_REFERENCE = 0,
    COMPOUND_REFERENCE = 1,
    REFERENCE_MODE_SELECT = 2,
    REFERENCE_MODES = 3,
} ReferenceMode;

typedef enum RefreshFrameContextMode
{
    /**
    * Frame context updates are disabled
    */
    REFRESH_FRAME_CONTEXT_DISABLED,
    /**
    * Update frame context to values resulting from backward probability
    * updates based on entropy/counts in the decoded frame
    */
    REFRESH_FRAME_CONTEXT_BACKWARD,
} RefreshFrameContextMode;

//**********************************************************************************************************************//
// aom_codec.h
/*!\brief Algorithm return codes */
typedef enum AomCodecErr
{
    /*!\brief Operation completed without error */
    AOM_CODEC_OK,
    /*!\brief Unspecified error */
    AOM_CODEC_ERROR,
    /*!\brief Memory operation failed */
    AOM_CODEC_MEM_ERROR,
    /*!\brief ABI version mismatch */
    AOM_CODEC_ABI_MISMATCH,
    /*!\brief Algorithm does not have required capability */
    AOM_CODEC_INCAPABLE,
    /*!\brief The given bitstream is not supported.
    *
    * The bitstream was unable to be parsed at the highest level. The decoder
    * is unable to proceed. This error \ref SHOULD be treated as fatal to the
    * stream. */
    AOM_CODEC_UNSUP_BITSTREAM,
    /*!\brief Encoded bitstream uses an unsupported feature
    *
    * The decoder does not implement a feature required by the encoder. This
    * return code should only be used for features that prevent future
    * pictures from being properly decoded. This error \ref MAY be treated as
    * fatal to the stream or \ref MAY be treated as fatal to the current GOP.
    */
    AOM_CODEC_UNSUP_FEATURE,
    /*!\brief The coded data for this stream is corrupt or incomplete
    *
    * There was a problem decoding the current frame.  This return code
    * should only be used for failures that prevent future pictures from
    * being properly decoded. This error \ref MAY be treated as fatal to the
    * stream or \ref MAY be treated as fatal to the current GOP. If decoding
    * is continued for the current GOP, artifacts may be present.
    */
    AOM_CODEC_CORRUPT_FRAME,
    /*!\brief An application-supplied parameter is not valid.
    *
    */
    AOM_CODEC_INVALID_PARAM,
    /*!\brief An iterator reached the end of list.
    *
    */
    AOM_CODEC_LIST_END
} AomCodecErr;

//**********************************************************************************************************************//
// Common_data.h
static const int32_t intra_mode_context[INTRA_MODES] = {
    0, 1, 2, 3, 4, 4, 4, 4, 3, 0, 1, 2, 0,
};

static const TxSize txsize_sqr_map[TX_SIZES_ALL] = {
    TX_4X4,    // TX_4X4
    TX_8X8,    // TX_8X8
    TX_16X16,  // TX_16X16
    TX_32X32,  // TX_32X32
    TX_64X64,  // TX_64X64
    TX_4X4,    // TX_4X8
    TX_4X4,    // TX_8X4
    TX_8X8,    // TX_8X16
    TX_8X8,    // TX_16X8
    TX_16X16,  // TX_16X32
    TX_16X16,  // TX_32X16
    TX_32X32,  // TX_32X64
    TX_32X32,  // TX_64X32
    TX_4X4,    // TX_4X16
    TX_4X4,    // TX_16X4
    TX_8X8,    // TX_8X32
    TX_8X8,    // TX_32X8
    TX_16X16,  // TX_16X64
    TX_16X16,  // TX_64X16
};
static const TxSize txsize_sqr_up_map[TX_SIZES_ALL] = {
    TX_4X4,    // TX_4X4
    TX_8X8,    // TX_8X8
    TX_16X16,  // TX_16X16
    TX_32X32,  // TX_32X32
    TX_64X64,  // TX_64X64
    TX_8X8,    // TX_4X8
    TX_8X8,    // TX_8X4
    TX_16X16,  // TX_8X16
    TX_16X16,  // TX_16X8
    TX_32X32,  // TX_16X32
    TX_32X32,  // TX_32X16
    TX_64X64,  // TX_32X64
    TX_64X64,  // TX_64X32
    TX_16X16,  // TX_4X16
    TX_16X16,  // TX_16X4
    TX_32X32,  // TX_8X32
    TX_32X32,  // TX_32X8
    TX_64X64,  // TX_16X64
    TX_64X64,  // TX_64X16
};

// above and left partition
typedef struct PartitionContext
{
    PartitionContextType above;
    PartitionContextType left;
}PartitionContext;
// Generates 5 bit field in which each bit set to 1 represents
// a BlockSize partition  11111 means we split 128x128, 64x64, 32x32, 16x16
// and 8x8.  10000 means we just split the 128x128 to 64x64
/* clang-format off */
static const struct
{
    PartitionContextType above;
    PartitionContextType left;
} partition_context_lookup[BlockSizeS_ALL] = {
{ 31, 31 },  // 4X4   - {0b11111, 0b11111}
{ 31, 30 },  // 4X8   - {0b11111, 0b11110}
{ 30, 31 },  // 8X4   - {0b11110, 0b11111}
{ 30, 30 },  // 8X8   - {0b11110, 0b11110}
{ 30, 28 },  // 8X16  - {0b11110, 0b11100}
{ 28, 30 },  // 16X8  - {0b11100, 0b11110}
{ 28, 28 },  // 16X16 - {0b11100, 0b11100}
{ 28, 24 },  // 16X32 - {0b11100, 0b11000}
{ 24, 28 },  // 32X16 - {0b11000, 0b11100}
{ 24, 24 },  // 32X32 - {0b11000, 0b11000}
{ 24, 16 },  // 32X64 - {0b11000, 0b10000}
{ 16, 24 },  // 64X32 - {0b10000, 0b11000}
{ 16, 16 },  // 64X64 - {0b10000, 0b10000}
{ 16, 0 },   // 64X128- {0b10000, 0b00000}
{ 0, 16 },   // 128X64- {0b00000, 0b10000}
{ 0, 0 },    // 128X128-{0b00000, 0b00000}
{ 31, 28 },  // 4X16  - {0b11111, 0b11100}
{ 28, 31 },  // 16X4  - {0b11100, 0b11111}
{ 30, 24 },  // 8X32  - {0b11110, 0b11000}
{ 24, 30 },  // 32X8  - {0b11000, 0b11110}
{ 28, 16 },  // 16X64 - {0b11100, 0b10000}
{ 16, 28 },  // 64X16 - {0b10000, 0b11100}
};
/* clang-format on */

// Width/height lookup tables in units of various block sizes
static const uint8_t block_size_wide[BlockSizeS_ALL] = {
    4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32,
    64, 64, 64, 128, 128, 4, 16, 8, 32, 16, 64
};

static const uint8_t block_size_high[BlockSizeS_ALL] = {
    4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 64,
    32, 64, 128, 64, 128, 16, 4, 32, 8, 64, 16
};

// AOMMIN(3, AOMMIN(b_width_log2(bsize), b_height_log2(bsize)))
static const uint8_t size_group_lookup[BlockSizeS_ALL] = {
    0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 2, 2
};

static const uint8_t num_pels_log2_lookup[BlockSizeS_ALL] = {
    4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12, 13, 13, 14, 6, 6, 8, 8, 10, 10
};
static const TxSize max_txsize_lookup[BlockSizeS_ALL] = {
    //                   4X4
    TX_4X4,
    // 4X8,    8X4,      8X8
    TX_4X4, TX_4X4, TX_8X8,
    // 8X16,   16X8,     16X16
    TX_8X8, TX_8X8, TX_16X16,
    // 16X32,  32X16,    32X32
    TX_16X16, TX_16X16, TX_32X32,
    // 32X64,  64X32,
    TX_32X32, TX_32X32,
    // 64X64
    TX_64X64,
    // 64x128, 128x64,   128x128
    TX_64X64, TX_64X64, TX_64X64,
    // 4x16,   16x4,     8x32
    TX_4X4, TX_4X4, TX_8X8,
    // 32x8,   16x64     64x16
    TX_8X8, TX_16X16, TX_16X16
};

static const TxSize max_txsize_rect_lookup[BlockSizeS_ALL] = {
    // 4X4
    TX_4X4,
    // 4X8,    8X4,      8X8
    TX_4X8, TX_8X4, TX_8X8,
    // 8X16,   16X8,     16X16
    TX_8X16, TX_16X8, TX_16X16,
    // 16X32,  32X16,    32X32
    TX_16X32, TX_32X16, TX_32X32,
    // 32X64,  64X32,
    TX_32X64, TX_64X32,
    // 64X64
    TX_64X64,
    // 64x128, 128x64,   128x128
    TX_64X64, TX_64X64, TX_64X64,
    // 4x16,   16x4,
    TX_4X16, TX_16X4,
    // 8x32,   32x8
    TX_8X32, TX_32X8,
    // 16x64,  64x16
    TX_16X64, TX_64X16
};

// Transform block width in unit
static const int32_t tx_size_wide_unit[TX_SIZES_ALL] = {
    1, 2, 4, 8, 16, 1, 2, 2, 4, 4, 8, 8, 16, 1, 4, 2, 8, 4, 16,
};
// Transform block height in unit
static const int32_t tx_size_high_unit[TX_SIZES_ALL] = {
    1, 2, 4, 8, 16, 2, 1, 4, 2, 8, 4, 16, 8, 4, 1, 8, 2, 16, 4,
};

static const TxSize sub_tx_size_map[TX_SIZES_ALL] = {
    TX_4X4,    // TX_4X4
    TX_4X4,    // TX_8X8
    TX_8X8,    // TX_16X16
    TX_16X16,  // TX_32X32
    TX_32X32,  // TX_64X64
    TX_4X4,    // TX_4X8
    TX_4X4,    // TX_8X4
    TX_8X8,    // TX_8X16
    TX_8X8,    // TX_16X8
    TX_16X16,  // TX_16X32
    TX_16X16,  // TX_32X16
    TX_32X32,  // TX_32X64
    TX_32X32,  // TX_64X32
    TX_4X8,    // TX_4X16
    TX_8X4,    // TX_16X4
    TX_8X16,   // TX_8X32
    TX_16X8,   // TX_32X8
    TX_16X32,  // TX_16X64
    TX_32X16,  // TX_64X16
};
static const TxSize txsize_horz_map[TX_SIZES_ALL] = {
    TX_4X4,    // TX_4X4
    TX_8X8,    // TX_8X8
    TX_16X16,  // TX_16X16
    TX_32X32,  // TX_32X32
    TX_64X64,  // TX_64X64
    TX_4X4,    // TX_4X8
    TX_8X8,    // TX_8X4
    TX_8X8,    // TX_8X16
    TX_16X16,  // TX_16X8
    TX_16X16,  // TX_16X32
    TX_32X32,  // TX_32X16
    TX_32X32,  // TX_32X64
    TX_64X64,  // TX_64X32
    TX_4X4,    // TX_4X16
    TX_16X16,  // TX_16X4
    TX_8X8,    // TX_8X32
    TX_32X32,  // TX_32X8
    TX_16X16,  // TX_16X64
    TX_64X64,  // TX_64X16
};

static const TxSize txsize_vert_map[TX_SIZES_ALL] = {
    TX_4X4,    // TX_4X4
    TX_8X8,    // TX_8X8
    TX_16X16,  // TX_16X16
    TX_32X32,  // TX_32X32
    TX_64X64,  // TX_64X64
    TX_8X8,    // TX_4X8
    TX_4X4,    // TX_8X4
    TX_16X16,  // TX_8X16
    TX_8X8,    // TX_16X8
    TX_32X32,  // TX_16X32
    TX_16X16,  // TX_32X16
    TX_64X64,  // TX_32X64
    TX_32X32,  // TX_64X32
    TX_16X16,  // TX_4X16
    TX_4X4,    // TX_16X4
    TX_32X32,  // TX_8X32
    TX_8X8,    // TX_32X8
    TX_64X64,  // TX_16X64
    TX_16X16,  // TX_64X16
};
static const uint8_t mi_size_wide_log2[BlockSizeS_ALL] = {
    0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 0, 2, 1, 3, 2, 4
};
static const uint8_t mi_size_high_log2[BlockSizeS_ALL] = {
    0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 2, 0, 3, 1, 4, 2
};

typedef struct SgrParamsType
{
    int32_t r[2];  // radii
    int32_t s[2];  // sgr parameters for r[0] and r[1], based on GenSgrprojVtable()
} SgrParamsType;

//**********************************************************************************************************************//
// blockd.h
typedef enum FrameType
{
    KEY_FRAME = 0,
    INTER_FRAME = 1,
    INTRA_ONLY_FRAME = 2,  // replaces intra-only
    S_FRAME = 3,
    FRAME_TYPES,
} FrameType;

typedef int8_t MvReferenceFrame;

// Number of transform types in each set type

static const int32_t av1_num_ext_tx_set[EXT_TX_SET_TYPES] = {
    1, 2, 5, 7, 12, 16,
};

static const int32_t av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
{ 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
{ 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
};

static INLINE TxSetType get_ext_tx_set_type(TxSize tx_size, int32_t is_inter,
    int32_t use_reduced_set) {
    const TxSize tx_size_sqr_up = txsize_sqr_up_map[tx_size];

    if (tx_size_sqr_up > TX_32X32) return EXT_TX_SET_DCTONLY;
    if (tx_size_sqr_up == TX_32X32)
        return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DCTONLY;
    if (use_reduced_set)
        return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DTT4_IDTX;
    const TxSize tx_size_sqr = txsize_sqr_map[tx_size];
    if (is_inter) {
        return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT9_IDTX_1DDCT
            : EXT_TX_SET_ALL16);
    }
    else {
        return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT4_IDTX
            : EXT_TX_SET_DTT4_IDTX_1DDCT);
    }
}
static INLINE int32_t get_ext_tx_types(TxSize tx_size, int32_t is_inter,
    int32_t use_reduced_set) {
    const int32_t set_type = get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
    return av1_num_ext_tx_set[set_type];
}
// Maps tx set types to the indices.
static const int32_t ext_tx_set_index[2][EXT_TX_SET_TYPES] = {
    { // Intra
        0, -1, 2, 1, -1, -1 },
        { // Inter
            0, 3, -1, -1, 2, 1 },
};

static INLINE int32_t get_ext_tx_set(TxSize tx_size, int32_t is_inter,
    int32_t use_reduced_set) {
    const TxSetType set_type =
        get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
    return ext_tx_set_index[is_inter][set_type];
}

static INLINE int32_t is_inter_compound_mode(PredictionMode mode) {
    return mode >= NEAREST_NEARESTMV && mode <= NEW_NEWMV;
}
static INLINE int is_inter_singleref_mode(PredictionMode mode) {
    return mode >= SINGLE_INTER_MODE_START && mode < SINGLE_INTER_MODE_END;
}

//**********************************************************************************************************************//
// encoder.h
typedef enum FrameContextIndex
{
    // regular inter frame
    REGULAR_FRAME = 0,
    // alternate reference frame
    ARF_FRAME = 1,
    // overlay frame
    OVERLAY_FRAME = 2,
    // golden frame
    GLD_FRAME = 3,
    // backward reference frame
    BRF_FRAME = 4,
    // extra alternate reference frame
    EXT_ARF_FRAME = 5,
    FRAME_CONTEXT_INDEXES
} FrameContextIndex;

//**********************************************************************************************************************//
// common.h
#define av1_zero(dest) memset(&(dest), 0, sizeof(dest))
//**********************************************************************************************************************//
// alloccommon.h
#define INVALID_IDX -1  // Invalid buffer index.

//**********************************************************************************************************************//
// quant_common.h
#define MINQ 0
#define MAXQ 255
#define QINDEX_RANGE (MAXQ - MINQ + 1)
#define QINDEX_BITS 8
// Total number of QM sets stored
#define QM_LEVEL_BITS 4
#define NUM_QM_LEVELS (1 << QM_LEVEL_BITS)
/* Range of QMS is between first and last value, with offset applied to inter
* blocks*/
#define DEFAULT_QM_Y 10
#define DEFAULT_QM_U 11
#define DEFAULT_QM_V 12
#define DEFAULT_QM_FIRST 5
#define DEFAULT_QM_LAST 9

//**********************************************************************************************************************//
// blockd.h
#define NO_FILTER_FOR_IBC 1  // Disable in-loop filters for frame with intrabc
//**********************************************************************************************************************//
// av1_loopfilter.h
#define MAX_LOOP_FILTER 63
#define MAX_SHARPNESS 7
#define SIMD_WIDTH 16

struct LoopFilter {
    int32_t filter_level[2];
    int32_t filter_level_u;
    int32_t filter_level_v;

    int32_t sharpness_level;

    uint8_t mode_ref_delta_enabled;
    uint8_t mode_ref_delta_update;

    // 0 = Intra, Last, Last2+Last3,
    // GF, BRF, ARF2, ARF
    int8_t ref_deltas[REF_FRAMES];

    // 0 = ZERO_MV, MV
    int8_t mode_deltas[MAX_MODE_LF_DELTAS];
    int32_t combine_vert_horz_lf;
};

#define MAX_SEGMENTS 8
#define MAX_MB_PLANE 3

#define MAX_LOOP_FILTER 63
#define MAX_SHARPNESS 7

#define SIMD_WIDTH 16
// Need to align this structure so when it is declared and
// passed it can be loaded into vector registers.
typedef struct LoopFilterThresh
{
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, mblim[SIMD_WIDTH]);
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, lim[SIMD_WIDTH]);
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, hev_thr[SIMD_WIDTH]);
} LoopFilterThresh;

typedef struct LoopFilterInfoN
{
    LoopFilterThresh lfthr[MAX_LOOP_FILTER + 1];
    uint8_t lvl[MAX_MB_PLANE][MAX_SEGMENTS][2][REF_FRAMES][MAX_MODE_LF_DELTAS];
} LoopFilterInfoN;

//**********************************************************************************************************************//
// cdef.h
#define CDEF_STRENGTH_BITS 6

#define CDEF_PRI_STRENGTHS 16
#define CDEF_SEC_STRENGTHS 4

// Bits of precision used for the model
#define WARPEDMODEL_PREC_BITS 16
// The following constants describe the various precisions
// of different parameters in the global motion experiment.
//
// Given the general homography:
//      [x'     (a  b  c   [x
//  z .  y'  =   d  e  f *  y
//       1]      g  h  i)    1]
//
// Constants using the name ALPHA here are related to parameters
// a, b, d, e. Constants using the name TRANS are related
// to parameters c and f.
//
// Anything ending in PREC_BITS is the number of bits of precision
// to maintain when converting from double to integer.
//
// The ABS parameters are used to create an upper and lower bound
// for each parameter. In other words, after a parameter is integerized
// it is clamped between -(1 << ABS_XXX_BITS) and (1 << ABS_XXX_BITS).
//
// XXX_PREC_DIFF and XXX_DECODE_FACTOR
// are computed once here to prevent repetitive
// computation on the decoder side. These are
// to allow the global motion parameters to be encoded in a lower
// precision than the warped model precision. This means that they
// need to be changed to warped precision when they are decoded.
//
// XX_MIN, XX_MAX are also computed to avoid repeated computation

#define SUBEXPFIN_K 3
#define GM_TRANS_PREC_BITS 6
#define GM_ABS_TRANS_BITS 12
#define GM_ABS_TRANS_ONLY_BITS (GM_ABS_TRANS_BITS - GM_TRANS_PREC_BITS + 3)
#define GM_TRANS_PREC_DIFF (WARPEDMODEL_PREC_BITS - GM_TRANS_PREC_BITS)
#define GM_TRANS_ONLY_PREC_DIFF (WARPEDMODEL_PREC_BITS - 3)
#define GM_TRANS_DECODE_FACTOR (1 << GM_TRANS_PREC_DIFF)
#define GM_TRANS_ONLY_DECODE_FACTOR (1 << GM_TRANS_ONLY_PREC_DIFF)
#define GM_TRANS_ONLY_PREC_BITS 3

#define GM_ALPHA_PREC_BITS 15
#define GM_ABS_ALPHA_BITS 12
#define GM_ALPHA_PREC_DIFF (WARPEDMODEL_PREC_BITS - GM_ALPHA_PREC_BITS)
#define GM_ALPHA_DECODE_FACTOR (1 << GM_ALPHA_PREC_DIFF)

#define GM_ROW3HOMO_PREC_BITS 16
#define GM_ABS_ROW3HOMO_BITS 11
#define GM_ROW3HOMO_PREC_DIFF \
(WARPEDMODEL_ROW3HOMO_PREC_BITS - GM_ROW3HOMO_PREC_BITS)
#define GM_ROW3HOMO_DECODE_FACTOR (1 << GM_ROW3HOMO_PREC_DIFF)

#define GM_TRANS_MAX (1 << GM_ABS_TRANS_BITS)
#define GM_ALPHA_MAX (1 << GM_ABS_ALPHA_BITS)
#define GM_ROW3HOMO_MAX (1 << GM_ABS_ROW3HOMO_BITS)

#define GM_TRANS_MIN -GM_TRANS_MAX
#define GM_ALPHA_MIN -GM_ALPHA_MAX
#define GM_ROW3HOMO_MIN -GM_ROW3HOMO_MAX

#define USE_CUR_GM_REFMV 1

/* clang-format off */
typedef enum TransformationType
{
    IDENTITY = 0,      // identity transformation, 0-parameter
    TRANSLATION = 1,   // translational motion 2-parameter
    ROTZOOM = 2,       // simplified affine with rotation + zoom only, 4-parameter
    AFFINE = 3,        // affine, 6-parameter
    TRANS_TYPES,
} TransformationType;
// The order of values in the wmmat matrix below is best described
// by the homography:
//      [x'     (m2 m3 m0   [x
//  z .  y'  =   m4 m5 m1 *  y
//       1]      m6 m7 1)    1]
typedef struct EbWarpedMotionParams
{
    TransformationType wmtype;
    int32_t wmmat[8];
    int16_t alpha, beta, gamma, delta;
    int8_t invalid;
} EbWarpedMotionParams;

/*! Scale factors and scaling function pointers  when reference and current frame dimensions are not equal */
typedef struct ScaleFactors {
    int32_t x_scale_fp;  // horizontal fixed point scale factor
    int32_t y_scale_fp;  // vertical fixed point scale factor
    int32_t x_step_q4;
    int32_t y_step_q4;

    int32_t(*scale_value_x)(int32_t val, const struct ScaleFactors *sf);
    int32_t(*scale_value_y)(int32_t val, const struct ScaleFactors *sf);
} ScaleFactors;

/* clang-format off */
static const EbWarpedMotionParams default_warp_params = {
    IDENTITY,
{ 0, 0, (1 << WARPEDMODEL_PREC_BITS), 0, 0, (1 << WARPEDMODEL_PREC_BITS), 0,
0 },
0, 0, 0, 0,
0,
};

/***********************************    AV1_OBU     ********************************/

//**********************************************************************************************************************//
//**********************************************************************************************************************//

#define YBITS_THSHLD                        50
#define YDC_THSHLD                          5
#define M6_YBITS_THSHLD                     80
#define M6_YDC_THSHLD                       10

#ifdef _WIN32
#define NOINLINE                __declspec ( noinline )
#define FORCE_INLINE            __forceinline
#else
#define NOINLINE                __attribute__(( noinline ))
#define FORCE_INLINE            __attribute__((always_inline))
#endif

#define EB_STRINGIZE( L )       #L
#define EB_MAKESTRING( M, L )   M( L )
#define $Line                   EB_MAKESTRING( EB_STRINGIZE, __LINE__ )
#define EB_SRC_LINE             __FILE__ "(" $Line ") : message "

// ***************************** Definitions *****************************
#define PM_DC_TRSHLD1                       10 // The threshold for DC to disable masking for DC

#define MAX_BITS_PER_FRAME            8000000
#define VAR_BASED_STAT_AREA_THRSLHD         (32*32)

#define ANTI_TRAILING_VAR_THRSLD         1000
#define MAX_VAR_BIAS               100
#define MEAN_DIFF_THRSHOLD         10
#define VAR_DIFF_THRSHOLD          10

#define HME_BIAS_X_THRSHLD1       64
#define HME_BIAS_Y_THRSHLD1       64
#define HME_BIAS_X_THRSHLD2       32
#define HME_BIAS_Y_THRSHLD2       32

#define ASPECT_RATIO_4_3    13           // Limit Ration to detect VGA resolutiosn
#define ASPECT_RATIO_16_9   17           // Limit Ration to detect UHD,1080p,720p ... or similar resolutions

#define ASPECT_RATIO_CLASS_0  0           // 4:3 aspect ratios
#define ASPECT_RATIO_CLASS_1  1           // 16:9 aspect ratios
#define ASPECT_RATIO_CLASS_2  2           // Other aspect ratios

#define SC_FRAMES_TO_IGNORE     1000 // The speed control algorith starts after SC_FRAMES_TO_IGNORE number frames.
#define SC_FRAMES_INTERVAL_SPEED      60 // The speed control Interval To Check the speed
#define SC_FRAMES_INTERVAL_T1         60 // The speed control Interval Threshold1
#define SC_FRAMES_INTERVAL_T2        180 // The speed control Interval Threshold2
#define SC_FRAMES_INTERVAL_T3        120 // The speed control Interval Threshold3

#define SC_SPEED_T2             1250 // speed level thershold. If speed is higher than target speed x SC_SPEED_T2, a slower mode is selected (+25% x 1000 (for precision))
#define SC_SPEED_T1              750 // speed level thershold. If speed is less than target speed x SC_SPEED_T1, a fast mode is selected (-25% x 1000 (for precision))
#define EB_CMPLX_CLASS           uint8_t
#define CMPLX_LOW                0
#define CMPLX_MEDIUM             1
#define CMPLX_HIGH               2
#define CMPLX_VHIGH              3
#define CMPLX_NOISE              4
#define EB_NORMAL_LATENCY        0
#define EB_LOW_LATENCY           1

typedef enum EbBitFieldMasks
{
    BITMASK_0 = 1,
    BITMASK_1 = 2,
    BITMASK_2 = 4,
    BITMASK_3 = 8
} EbBitFieldMasks;

// CLEAN_BASIS_FUNCTIONS
#define CLEAN_BASIS_FUNCTIONS_VAR_TRSHLD 10
#define CLEAN_BASIS_FUNCTIONS_NZCOEF_TRSHLD0 10
#define CLEAN_BASIS_FUNCTIONS_NZCOEF_TRSHLD1 15
#define CLEAN_BASIS_FUNCTIONS_NZCOEF_TRSHLD2 20
// Anti-contouring
#define C3_TRSHLF_N                                    45
#define C3_TRSHLF_D                                    10
#define C4_TRSHLF_N                                    35
#define C4_TRSHLF_D                                    10

#define C1_TRSHLF_4K_N                                45
#define C1_TRSHLF_4K_D                                10
#define C2_TRSHLF_4K_N                                35
#define C2_TRSHLF_4K_D                                10

#define AC_ENERGY_BASED_4K_ANTI_CONTOURING_QP_DELTA     3
#define AC_ENERGY_BASED_4K_ANTI_CONTOURING_MIN_QP       22

#define C1_TRSHLF_N       1
#define C1_TRSHLF_D       1
#define C2_TRSHLF_N       16
#define C2_TRSHLF_D       10

#define CHANGE_LAMBDA_FOR_AURA   0x01
#define RESTRICT_CUS_AND_MODIFY_COST  0x02

#define ANTI_CONTOURING_TH_0     16 * 16
#define ANTI_CONTOURING_TH_1     32 * 32
#define ANTI_CONTOURING_TH_2 2 * 32 * 32

#define ANTI_CONTOURING_DELTA_QP_0  -3
#define ANTI_CONTOURING_DELTA_QP_1  -9
#define ANTI_CONTOURING_DELTA_QP_2  -11

#define AC_ENERGY_BASED_ANTI_CONTOURING_QP_DELTA 11
#define AC_ENERGY_BASED_ANTI_CONTOURING_MIN_QP 20
#define ANTI_CONTOURING_LUMA_T1                40
#define ANTI_CONTOURING_LUMA_T2                180

#define VAR_BASED_DETAIL_PRESERVATION_SELECTOR_THRSLHD         (64*64)

#define LAST_BWD_FRAME     8
#define LAST_ALT_FRAME    16

#define MAX_NUM_TOKENS          200

#define LAD_DISABLE                       0
#define INIT_RC_OPT_G1                    1
#define INIT_RC_OPT_G2                    1
#define HIST_OPT                          2 // 1 is intrinsic, 2 is C
#define ENABLE_8x8                        0

#define    Log2f                              Log2f_SSE2

#define INPUT_SIZE_576p_TH                  0x90000        // 0.58 Million
#define INPUT_SIZE_1080i_TH                 0xB71B0        // 0.75 Million
#define INPUT_SIZE_1080p_TH                 0x1AB3F0    // 1.75 Million
#define INPUT_SIZE_4K_TH                    0x29F630    // 2.75 Million
#define INPUT_SIZE_8K_TH                    0xA7D8C0    // 11 Million

/** Redefine ASSERT() to avoid warnings
*/
#if defined _DEBUG || _DEBUG_
#include <assert.h>
#define ASSERT assert
#elif defined _DEBUG
#define ASSERT assert
#else
#define ASSERT(exp) ((void)sizeof(exp))
#endif

#define    INTERPOLATION_NEED  4
#define    BUFF_PITCH          (INTERPOLATION_NEED*2+64)
#define    ME_FILTER_TAP       4
#define    SUB_SAD_SEARCH      0
#define    FULL_SAD_SEARCH     1
#define    SSD_SEARCH          2
/************************ INPUT CLASS **************************/

#define EbInputResolution             uint8_t
#define INPUT_SIZE_576p_RANGE_OR_LOWER     0
#define INPUT_SIZE_1080i_RANGE             1
#define INPUT_SIZE_1080p_RANGE             2
#define INPUT_SIZE_4K_RANGE                3
#define INPUT_SIZE_COUNT                   INPUT_SIZE_4K_RANGE + 1

/** The EbPtr type is intended to be used to pass pointers to and from the eBrisk
API.  This is a 32 bit pointer and is aligned on a 32 bit word boundary.
*/
typedef void *EbPtr;

/** The EbString type is intended to be used to pass "C" type strings to and
from the eBrisk API.  The EbString type is a 32 bit pointer to a zero terminated
string.  The pointer is word aligned and the string is byte aligned.
*/
typedef char * EbString;

/** The EbByte type is intended to be used to pass arrays of bytes such as
buffers to and from the eBrisk API.  The EbByte type is a 32 bit pointer.
The pointer is word aligned and the buffer is byte aligned.
*/
typedef uint8_t * EbByte;

/** The EB_SAMPLE type is intended to be used to pass arrays of bytes such as
buffers to and from the eBrisk API.  The EbByte type is a 32 bit pointer.
The pointer is word aligned and the buffer is byte aligned.
*/

/** The EbBitDepthEnum type is used to describe the bitdepth of video data.
*/
typedef enum EbBitDepthEnum
{
    EB_8BIT = 8,
    EB_10BIT = 10,
    EB_12BIT = 12,
    EB_14BIT = 14,
    EB_16BIT = 16,
    EB_32BIT = 32
} EbBitDepthEnum;
#if HBD_CLEAN_UP
/** The MD_BIT_DEPTH_MODE type is used to describe the bitdepth of MD path.
*/

typedef enum MD_BIT_DEPTH_MODE
{
    EB_8_BIT_MD     = 0,    // 8bit mode decision
    EB_10_BIT_MD    = 1,    // 10bit mode decision
    EB_DUAL_BIT_MD  = 2     // Auto: 8bit & 10bit mode decision
} MD_BIT_DEPTH_MODE;
#endif
/** The EB_GOP type is used to describe the hierarchical coding structure of
Groups of Pictures (GOP) units.
*/
#define EbPred                 uint8_t
#define EB_PRED_LOW_DELAY_P     0
#define EB_PRED_LOW_DELAY_B     1
#define EB_PRED_RANDOM_ACCESS   2
#define EB_PRED_TOTAL_COUNT     3
#define EB_PRED_INVALID         0xFF

/** The EB_SLICE type is used to describe the slice prediction type.
*/

#define EB_SLICE        uint8_t
#define B_SLICE         0
#define P_SLICE         1
#define I_SLICE         2
#define IDR_SLICE       3
#define INVALID_SLICE   0xFF

/** The EbPictStruct type is used to describe the picture structure.
*/
#define EbPictStruct           uint8_t
#define PROGRESSIVE_PICT_STRUCT  0
#define TOP_FIELD_PICT_STRUCT    1
#define BOTTOM_FIELD_PICT_STRUCT 2

/** The EbModeType type is used to describe the PU type.
*/
typedef uint8_t EbModeType;
#define INTER_MODE 1
#define INTRA_MODE 2

#define INVALID_MODE 0xFFu

/** INTRA_4x4 offsets
*/
static const uint8_t INTRA_4x4_OFFSET_X[4] = { 0, 4, 0, 4 };
static const uint8_t INTRA_4x4_OFFSET_Y[4] = { 0, 0, 4, 4 };

/** The EbPartMode type is used to describe the CU partition size.
*/
typedef uint8_t EbPartMode;
#define SIZE_2Nx2N 0
#define SIZE_2NxN  1
#define SIZE_Nx2N  2
#define SIZE_NxN   3
#define SIZE_2NxnU 4
#define SIZE_2NxnD 5
#define SIZE_nLx2N 6
#define SIZE_nRx2N 7
#define SIZE_PART_MODE 8

/** The EbIntraRefreshType is used to describe the intra refresh type.
*/
typedef enum EbIntraRefreshType
{
    NO_REFRESH = 0,
    CRA_REFRESH = 1,
    IDR_REFRESH = 2
}EbIntraRefreshType;

#define SIZE_2Nx2N_PARTITION_MASK   (1 << SIZE_2Nx2N)
#define SIZE_2NxN_PARTITION_MASK    (1 << SIZE_2NxN)
#define SIZE_Nx2N_PARTITION_MASK    (1 << SIZE_Nx2N)
#define SIZE_NxN_PARTITION_MASK     (1 << SIZE_NxN)
#define SIZE_2NxnU_PARTITION_MASK   (1 << SIZE_2NxnU)
#define SIZE_2NxnD_PARTITION_MASK   (1 << SIZE_2NxnD)
#define SIZE_nLx2N_PARTITION_MASK   (1 << SIZE_nLx2N)
#define SIZE_nRx2N_PARTITION_MASK   (1 << SIZE_nRx2N)

/** The EbEncMode type is used to describe the encoder mode .
*/

#define EbEncMode     uint8_t
#define ENC_M0          0
#define ENC_M1          1
#define ENC_M2          2
#define ENC_M3          3
#define ENC_M4          4
#define ENC_M5          5
#define ENC_M6          6
#define ENC_M7          7
#define ENC_M8          8
#define ENC_M9          9
#define ENC_M10         10
#define ENC_M11         11
#define ENC_M12         12

#define MAX_SUPPORTED_MODES 13

#define SPEED_CONTROL_INIT_MOD ENC_M4;
/** The EB_TUID type is used to identify a TU within a CU.
*/
typedef enum EbTuSize
{
    TU_2Nx2N       = 0,
    TU_NxN_0       = 1,
    TU_NxN_1       = 2,
    TU_NxN_2       = 3,
    TU_NxN_3       = 4,
    TU_N2xN2_0     = 5,
    TU_N2xN2_1     = 6,
    TU_N2xN2_2     = 7,
    TU_N2xN2_3     = 8,
    INVALID_TUSIZE = ~0
}EbTuSize;

#define TU_2Nx2N_PARTITION_MASK     (1 << TU_2Nx2N)
#define TU_NxN_0_PARTITION_MASK     (1 << TU_NxN_0)
#define TU_NxN_1_PARTITION_MASK     (1 << TU_NxN_1)
#define TU_NxN_2_PARTITION_MASK     (1 << TU_NxN_2)
#define TU_NxN_3_PARTITION_MASK     (1 << TU_NxN_3)
#define TU_N2xN2_0_PARTITION_MASK   (1 << TU_N2xN2_0)
#define TU_N2xN2_1_PARTITION_MASK   (1 << TU_N2xN2_1)
#define TU_N2xN2_2_PARTITION_MASK   (1 << TU_N2xN2_2)
#define TU_N2xN2_3_PARTITION_MASK   (1 << TU_N2xN2_3)

#define EbReflist            uint8_t
#define REF_LIST_0             0
#define REF_LIST_1             1
#define TOTAL_NUM_OF_REF_LISTS 2
#define INVALID_LIST           0xFF

#define EbPredDirection         uint8_t
#define UNI_PRED_LIST_0          0
#define UNI_PRED_LIST_1          1
#define BI_PRED                  2
#define EB_PREDDIRECTION_TOTAL   3
#define INVALID_PRED_DIRECTION   0xFF

#define UNI_PRED_LIST_0_MASK    (1 << UNI_PRED_LIST_0)
#define UNI_PRED_LIST_1_MASK    (1 << UNI_PRED_LIST_1)
#define BI_PRED_MASK            (1 << BI_PRED)

// The EB_QP_OFFSET_MODE type is used to describe the QP offset
#define EB_FRAME_CARACTERICTICS uint8_t
#define EB_FRAME_CARAC_0           0
#define EB_FRAME_CARAC_1           1
#define EB_FRAME_CARAC_2           2
#define EB_FRAME_CARAC_3           3
#define EB_FRAME_CARAC_4           4

static const uint8_t QP_OFFSET_WEIGHT[3][4] = { // [Slice Type][QP Offset Weight Level]
    { 9, 8, 7, 6 },
    { 9, 8, 7, 6 },
    { 10, 9, 8, 7 }
};

#if PAL_SUP
#define  MAX_PAL_CAND   14
typedef struct {
    // Value of base colors for Y, U, and V
    uint16_t palette_colors[3 * PALETTE_MAX_SIZE];
    // Number of base colors for Y (0) and UV (1)
    uint8_t palette_size[2];

} PaletteModeInfo;

typedef struct {
    PaletteModeInfo pmi;
    uint8_t  *color_idx_map;
} PaletteInfo;
#endif
/** The EB_NULL type is used to define the C style NULL pointer.
*/
#define EB_NULL ((void*) 0)

/** The EbHandle type is used to define OS object handles for threads,
semaphores, mutexs, etc.
*/
typedef void * EbHandle;

/**
object_ptr is a EbPtr to the object being constructed.
object_init_data_ptr is a EbPtr to a data structure used to initialize the object.
*/
typedef EbErrorType(*EbCreator)(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);

#define INVALID_MV            0x80008000 //0xFFFFFFFF    //ICOPY They changed this to 0x80008000
#define BLKSIZE 64

/***************************************
* Generic linked list data structure for passing data into/out from the library
***************************************/
// Reserved types for lib's internal use. Must be less than EB_EXT_TYPE_BASE
#define       EB_TYPE_PIC_TIMING_SEI         0
#define       EB_TYPE_BUFFERING_PERIOD_SEI   1
#define       EB_TYPE_RECOVERY_POINT_SEI     2
#define       EB_TYPE_UNREG_USER_DATA_SEI    3
#define       EB_TYPE_REG_USER_DATA_SEI      4
#define       EB_TYPE_PIC_STRUCT             5             // It is a requirement (for the application) that if pictureStruct is present for 1 picture it shall be present for every picture
#define       EB_TYPE_INPUT_PICTURE_DEF      6

#define       EB_TYPE_HIERARCHICAL_LEVELS  100
#define       EB_TYPE_PRED_STRUCTURE       101

typedef int32_t EbLinkedListType;

typedef struct EbLinkedListNode
{
    void*                     app;                       // points to an application object this node is associated
                                                            // with. this is an opaque pointer to the encoder lib, but
                                                            // release_cb_fnc_ptr may need to access it.
    EbLinkedListType       type;                      // type of data pointed by "data" member variable
    uint32_t                    size;                      // size of (data)
    EbBool                   passthrough;               // whether this is passthrough data from application
    void(*release_cb_fnc_ptr)(struct EbLinkedListNode*); // callback to be executed by encoder when picture reaches end of pipeline, or
                                                        // when aborting. However, at end of pipeline encoder shall
                                                        // NOT invoke this callback if passthrough is TRUE (but
                                                        // still needs to do so when aborting)
    void                     *data;                      // pointer to application's data
    struct EbLinkedListNode  *next;                      // pointer to next node (null when last)
} EbLinkedListNode;

typedef enum DistCalcType
{
    DIST_CALC_RESIDUAL = 0,    // SSE(Coefficients - ReconCoefficients)
    DIST_CALC_PREDICTION = 1,    // SSE(Coefficients) *Note - useful in modes that don't send residual coeff bits
    DIST_CALC_TOTAL = 2
} DistCalcType;

typedef enum EbPtrType
{
    EB_N_PTR        = 0,     // malloc'd pointer
    EB_C_PTR        = 1,     // calloc'd pointer
    EB_A_PTR        = 2,     // malloc'd pointer aligned
    EB_MUTEX        = 3,     // mutex
    EB_SEMAPHORE    = 4,     // semaphore
    EB_THREAD       = 5,      // thread handle
    EB_PTR_TYPE_TOTAL,
} EbPtrType;

typedef struct EbMemoryMapEntry
{
    EbPtr                    ptr;            // points to a memory pointer
    EbPtrType                ptr_type;       // pointer type
    EbPtr                    prev_entry;     // pointer to the prev entry
} EbMemoryMapEntry;

// Rate Control
#define THRESHOLD1QPINCREASE     1
#define THRESHOLD2QPINCREASE     2
#define EB_IOS_POINT            uint8_t
#define OIS_VERY_FAST_MODE       0
#define OIS_FAST_MODE            1
#define OIS_MEDUIM_MODE          2
#define OIS_COMPLEX_MODE         3
#define OIS_VERY_COMPLEX_MODE    4
// Display Total Memory at the end of the memory allocations
#define DISPLAY_MEMORY                              0

extern    EbMemoryMapEntry          *app_memory_map;            // App Memory table
extern    uint32_t                  *app_memory_map_index;       // App Memory index
extern    uint64_t                  *total_app_memory;          // App Memory malloc'd

extern    EbMemoryMapEntry          *memory_map;               // library Memory table
extern    uint32_t                  *memory_map_index;          // library memory index
extern    uint64_t                  *total_lib_memory;          // library Memory malloc'd

extern    uint32_t                   lib_malloc_count;
extern    uint32_t                   lib_thread_count;
extern    uint32_t                   lib_semaphore_count;
extern    uint32_t                   lib_mutex_count;

extern    uint32_t                   app_malloc_count;

#define ALVALUE 64

#define EB_ADD_APP_MEM(pointer, size, pointer_class, count, release, return_type) \
    do { \
        if (!pointer) return return_type; \
        if (*(app_memory_map_index) >= MAX_APP_NUM_PTR) { \
            printf("Malloc has failed due to insuffucient resources"); \
            release(pointer); \
            return return_type; \
        } \
        app_memory_map[*(app_memory_map_index)].ptr_type = pointer_class; \
        app_memory_map[(*(app_memory_map_index))++].ptr = pointer; \
        *total_app_memory += (size + 7) / 8; \
        count++; \
    } while (0)

#define EB_APP_MALLOC(type, pointer, n_elements, pointer_class, return_type) \
    pointer = (type)malloc(n_elements); \
    EB_ADD_APP_MEM(pointer, n_elements, pointer_class, app_malloc_count, return_type);


#define EB_APP_MALLOC_NR(type, pointer, n_elements, pointer_class,return_type) \
    pointer = (type)malloc(n_elements); \
    EB_ADD_APP_MEM(pointer, n_elements, pointer_class, app_malloc_count, return_type);

#define ALVALUE 64

#define EB_CREATE_SEMAPHORE(pointer, initial_count, max_count) \
    do { \
        pointer = eb_create_semaphore(initial_count, max_count); \
        EB_ADD_MEM(pointer, 1, EB_SEMAPHORE); \
    }while (0)

#define EB_DESTROY_SEMAPHORE(pointer) \
    do { \
        if (pointer) { \
            eb_destroy_semaphore(pointer); \
            EB_REMOVE_MEM_ENTRY(pointer, EB_SEMAPHORE); \
            pointer = NULL; \
        } \
    }while (0)

#define EB_CREATE_MUTEX(pointer) \
    do { \
        pointer = eb_create_mutex(); \
        EB_ADD_MEM(pointer, 1, EB_MUTEX); \
    } while (0)

#define EB_DESTROY_MUTEX(pointer) \
    do { \
        if (pointer) { \
            eb_destroy_mutex(pointer); \
            EB_REMOVE_MEM_ENTRY(pointer, EB_MUTEX); \
            pointer = NULL; \
        } \
    } while (0)


#define EB_MEMORY() \
printf("Total Number of Mallocs in Library: %d\n", lib_malloc_count); \
printf("Total Number of Threads in Library: %d\n", lib_thread_count); \
printf("Total Number of Semaphore in Library: %d\n", lib_semaphore_count); \
printf("Total Number of Mutex in Library: %d\n", lib_mutex_count); \
printf("Total Library Memory: %.2lf KB\n\n",*total_lib_memory/(double)1024);

#define EB_APP_MEMORY() \
printf("Total Number of Mallocs in App: %d\n", app_malloc_count); \
printf("Total App Memory: %.2lf KB\n\n",*total_app_memory/(double)1024);

#define RSIZE_MAX_MEM      ( 256UL << 20 )     /* 256MB */

#define EXPORT_SYMBOL(sym)

#ifndef _ERRNO_T_DEFINED
#define _ERRNO_T_DEFINED
typedef int32_t errno_t;
#endif  /* _ERRNO_T_DEFINED */

extern void
    eb_memcpy(void  *dst_ptr, void  *src_ptr, size_t size);

#define EB_MEMCPY(dst, src, size) \
    eb_memcpy(dst, src, size)

#define EB_MEMSET(dst, val, count) \
memset(dst, val, count)

//#ifdef __cplusplus
//}
//#endif // __cplusplus

/**************************************
* Callback Functions
**************************************/
typedef struct EbCallback
{
EbPtr appPrivateData;
EbPtr handle;
void(*ErrorHandler)(
    EbPtr handle,
    uint32_t errorCode);
} EbCallback;

// DEBUG MACROS
#define LIB_PRINTF_ENABLE                1

#if LIB_PRINTF_ENABLE
#define SVT_LOG printf
#else
#ifdef _MSC_VER
#define SVT_LOG(s, ...) printf("")
#else
#define SVT_LOG(s, ...) printf("",##__VA_ARGS__)
#endif
#endif

// Common Macros
#define UNUSED(x) (void)(x)

//***Profile, tier, level***
#define TOTAL_LEVEL_COUNT                           13

//***Encoding Parameters***
#define MAX_PICTURE_WIDTH_SIZE                      4672u
#define MAX_PICTURE_HEIGHT_SIZE                     2560u
#define MAX_PICTURE_WIDTH_SIZE_CH                   2336u
#define MAX_PICTURE_HEIGHT_SIZE_CH                  1280u
#define INTERNAL_BIT_DEPTH                          8 // to be modified
#define MAX_SAMPLE_VALUE                            ((1 << INTERNAL_BIT_DEPTH) - 1)
#define MAX_SAMPLE_VALUE_10BIT                      0x3FF
#define BLOCK_SIZE_64                                64u
#define LOG2F_MAX_LCU_SIZE                          6u
#define LOG2_64_SIZE                                6 // log2(BLOCK_SIZE_64)
#define MAX_LEVEL_COUNT                             5 // log2(BLOCK_SIZE_64) - log2(MIN_BLOCK_SIZE)
#if !ENHANCE_ATB
#define MAX_TU_DEPTH                                2
#endif
#define LOG_MIN_BLOCK_SIZE                          3
#define MIN_BLOCK_SIZE                              (1 << LOG_MIN_BLOCK_SIZE)
#define LOG_MIN_PU_SIZE                             2
#define MIN_PU_SIZE                                 (1 << LOG_MIN_PU_SIZE)
#define MAX_NUM_OF_PU_PER_CU                        1
#define MAX_NUM_OF_REF_PIC_LIST                     2
#define MAX_NUM_OF_PART_SIZE                        8
#define EB_MAX_LCU_DEPTH                            (((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 1) ? 1 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 2) ? 2 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 4) ? 3 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 8) ? 4 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 16) ? 5 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 32) ? 6 : 7)
#define MIN_CU_BLK_COUNT                            ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) * (BLOCK_SIZE_64 / MIN_BLOCK_SIZE))
#define MAX_NUM_OF_TU_PER_CU                        21
#define MIN_NUM_OF_TU_PER_CU                        5
#define MAX_LCU_ROWS                                ((MAX_PICTURE_HEIGHT_SIZE) / (BLOCK_SIZE_64))

#define MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE       ((MAX_PICTURE_WIDTH_SIZE + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64) * \
                                                ((MAX_PICTURE_HEIGHT_SIZE + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64)

//***Prediction Structure***
#define REF_LIST_MAX_DEPTH                          4 // NM - To be specified
#define MAX_TEMPORAL_LAYERS                         6
#define MAX_HIERARCHICAL_LEVEL                      6
#define MAX_REF_IDX                                 4
#define INVALID_POC                                 (((uint32_t) (~0)) - (((uint32_t) (~0)) >> 1))
#define MAX_ELAPSED_IDR_COUNT                       1024

typedef enum DownSamplingMethod
{
    ME_FILTERED_DOWNSAMPLED  = 0,
    ME_DECIMATED_DOWNSAMPLED = 1
} DownSamplingMethod;

//***Segments***
#define EB_SEGMENT_MIN_COUNT                        1
#define EB_SEGMENT_MAX_COUNT                        64
#define CU_MAX_COUNT                                85

#define EB_EVENT_MAX_COUNT                          20

#define MAX_INTRA_REFERENCE_SAMPLES                 (BLOCK_SIZE_64 << 2) + 1

#define MAX_INTRA_MODES                             35

#define _MVXT(mv) ( (int16_t)((mv) &  0xFFFF) )
#define _MVYT(mv) ( (int16_t)((mv) >> 16    ) )

//***MCP***
#define MaxChromaFilterTag          4
#define MaxVerticalLumaFliterTag    8
#define MaxHorizontalLumaFliterTag  8

#define MCPXPaddingOffset           16                                    // to be modified
#define MCPYPaddingOffset           16                                    // to be modified

#define InternalBitDepth            8                                     // to be modified
#define MAX_Sample_Value            ((1 << InternalBitDepth) - 1)
#define IF_Shift                    6                                     // to be modified
#define IF_Prec                     14                                    // to be modified
#define IF_Negative_Offset          (IF_Prec - 1)                         // to be modified
#define InternalBitDepthIncrement   (InternalBitDepth - 8)

#define MIN_QP_VALUE                     0
#define MAX_QP_VALUE                    63
#define MAX_CHROMA_MAP_QP_VALUE         63

//***Transforms***
#define TRANSFORMS_LUMA_FLAG        0
#define TRANSFORMS_CHROMA_FLAG      1
#define TRANSFORMS_COLOR_LEN        2
#define TRANSFORMS_LUMA_MASK        (1 << TRANSFORMS_LUMA_FLAG)
#define TRANSFORMS_CHROMA_MASK      (1 << TRANSFORMS_CHROMA_FLAG)
#define TRANSFORMS_FULL_MASK        ((1 << TRANSFORMS_LUMA_FLAG) | (1 << TRANSFORMS_CHROMA_FLAG))

#define TRANSFORMS_SIZE_32_FLAG     0
#define TRANSFORMS_SIZE_16_FLAG     1
#define TRANSFORMS_SIZE_8_FLAG      2
#define TRANSFORMS_SIZE_4_FLAG      3
#define TRANSFORMS_SIZE_LEN         4
#define TRANSFORM_MAX_SIZE          64
#define TRANSFORM_MIN_SIZE          4

#define BIT_INCREMENT_10BIT    2
#define BIT_INCREMENT_8BIT     0

#define TRANS_BIT_INCREMENT    0
#define QUANT_IQUANT_SHIFT     20 // Q(QP%6) * IQ(QP%6) = 2^20
#define QUANT_SHIFT            14 // Q(4) = 2^14
#define SCALE_BITS             15 // Inherited from TMuC, pressumably for fractional bit estimates in RDOQ
#define MAX_TR_DYNAMIC_RANGE   15 // Maximum transform dynamic range (excluding sign bit)
#define MAX_POS_16BIT_NUM      32767
#define MIN_NEG_16BIT_NUM      -32768
#define QUANT_OFFSET_I         171
#define QUANT_OFFSET_P         85
#define LOW_LCU_VARIANCE        10
#define MEDIUM_LCU_VARIANCE        50

/*********************************************************
* used for the first time, but not the last time interpolation filter
*********************************************************/
#define Shift1       InternalBitDepthIncrement
#define MinusOffset1 (1 << (IF_Negative_Offset + InternalBitDepthIncrement))
#if (InternalBitDepthIncrement == 0)
#define ChromaMinusOffset1 0
#else
#define ChromaMinusOffset1 MinusOffset1
#endif

/*********************************************************
* used for neither the first time nor the last time interpolation filter
*********************************************************/
#define Shift2       IF_Shift

/*********************************************************
* used for the first time, and also the last time interpolation filter
*********************************************************/
#define Shift3       IF_Shift
#define Offset3      (1<<(Shift3-1))

/*********************************************************
* used for not the first time, but the last time interpolation filter
*********************************************************/
#define Shift4       (IF_Shift + IF_Shift - InternalBitDepthIncrement)
#define Offset4      ((1 << (IF_Shift + IF_Negative_Offset)) + (1 << (Shift4 - 1)))
#if (InternalBitDepthIncrement == 0)
#define ChromaOffset4 (1 << (Shift4 - 1))
#else
#define ChromaOffset4 Offset4
#endif

/*********************************************************
* used for weighted sample prediction
*********************************************************/
#define Shift5       (IF_Shift - InternalBitDepthIncrement + 1)
#define Offset5      ((1 << (Shift5 - 1)) + (1 << (IF_Negative_Offset + 1)))
#if (InternalBitDepthIncrement == 0)
#define ChromaOffset5 (1 << (Shift5 - 1))
#else
#define ChromaOffset5 Offset5
#endif

/*********************************************************
* used for biPredCopy()
*********************************************************/
#define Shift6       (IF_Shift - InternalBitDepthIncrement)
#define MinusOffset6 (1 << IF_Negative_Offset)
#if (InternalBitDepthIncrement == 0)
#define ChromaMinusOffset6 0
#else
#define ChromaMinusOffset6 MinusOffset6
#endif

/*********************************************************
* 10bit case
*********************************************************/

#define  SHIFT1D_10BIT      6
#define  OFFSET1D_10BIT     32

#define  SHIFT2D1_10BIT     2
#define  OFFSET2D1_10BIT    (-32768)

#define  SHIFT2D2_10BIT     10
#define  OFFSET2D2_10BIT    524800

//BIPRED
#define  BI_SHIFT_10BIT         4
#define  BI_OFFSET_10BIT        8192//2^(14-1)

#define  BI_AVG_SHIFT_10BIT     5
#define  BI_AVG_OFFSET_10BIT    16400

#define  BI_SHIFT2D2_10BIT      6
#define  BI_OFFSET2D2_10BIT     0

// Noise detection
#define  NOISE_VARIANCE_TH                390

#define  EbPicnoiseClass    uint8_t
#define  PIC_NOISE_CLASS_INV  0 //not computed
#define  PIC_NOISE_CLASS_1    1 //No Noise
#define  PIC_NOISE_CLASS_2    2
#define  PIC_NOISE_CLASS_3    3
#define  PIC_NOISE_CLASS_3_1  4
#define  PIC_NOISE_CLASS_4    5
#define  PIC_NOISE_CLASS_5    6
#define  PIC_NOISE_CLASS_6    7
#define  PIC_NOISE_CLASS_7    8
#define  PIC_NOISE_CLASS_8    9
#define  PIC_NOISE_CLASS_9    10
#define  PIC_NOISE_CLASS_10   11 //Extreme Noise

// Intrinisc
#define INTRINSIC_SSE2                                1

// Enhance background macros for decimated 64x64
#define BEA_CLASS_0_0_DEC_TH 16 * 16    // 16x16 block size * 1
#define BEA_CLASS_0_DEC_TH     16 * 16 * 2    // 16x16 block size * 2
#define BEA_CLASS_1_DEC_TH     16 * 16 * 4    // 16x16 block size * 4
#define BEA_CLASS_2_DEC_TH     16 * 16 * 8    // 16x16 block size * 8

// Enhance background macros
#define BEA_CLASS_0_0_TH 8 * 8        // 8x8 block size * 1

#define BEA_CLASS_0_TH    8 * 8 * 2    // 8x8 block size * 2
#define BEA_CLASS_1_TH    8 * 8 * 4    // 8x8 block size * 4
#define BEA_CLASS_2_TH    8 * 8 * 8    // 8x8 block size * 8

#define UNCOVERED_AREA_ZZ_TH 4 * 4 * 14

#define BEA_CLASS_0_ZZ_COST     0
#define BEA_CLASS_0_1_ZZ_COST     3

#define BEA_CLASS_1_ZZ_COST    10
#define BEA_CLASS_2_ZZ_COST    20
#define BEA_CLASS_3_ZZ_COST    30
#define INVALID_ZZ_COST    (uint8_t) ~0

#define PM_NON_MOVING_INDEX_TH 23

#define QP_OFFSET_LCU_SCORE_0    0
#define QP_OFFSET_LCU_SCORE_1    50
#define QP_OFFSET_LCU_SCORE_2    100
#define UNCOVERED_AREA_ZZ_COST_TH 8
#define BEA_MIN_DELTA_QP_T00 1
#define BEA_MIN_DELTA_QP_T0  3
#define BEA_MIN_DELTA_QP_T1  5
#define BEA_MIN_DELTA_QP_T2  5
#define BEA_DISTANSE_RATIO_T0 900
#define BEA_DISTANSE_RATIO_T1 600
#define ACTIVE_PICTURE_ZZ_COST_TH 29

#define BEA_MAX_DELTA_QP 1

#define FAILING_MOTION_DELTA_QP            -5
#define FAILING_MOTION_VAR_THRSLHD        50
static const uint8_t intra_area_th_class_1[MAX_HIERARCHICAL_LEVEL][MAX_TEMPORAL_LAYERS] = { // [Highest Temporal Layer] [Temporal Layer Index]
    { 20 },
    { 30, 20 },
    { 40, 30, 20 },
    { 50, 40, 30, 20 },
    { 50, 40, 30, 20, 10 },
    { 50, 40, 30, 20, 10, 10 }
};

#define NON_MOVING_SCORE_0     0
#define NON_MOVING_SCORE_1    10
#define NON_MOVING_SCORE_2    20
#define NON_MOVING_SCORE_3    30
#define INVALID_NON_MOVING_SCORE (uint8_t) ~0

// Picture split into regions for analysis (SCD, Dynamic GOP)
#define CLASS_SUB_0_REGION_SPLIT_PER_WIDTH    1
#define CLASS_SUB_0_REGION_SPLIT_PER_HEIGHT    1

#define CLASS_1_REGION_SPLIT_PER_WIDTH        2
#define CLASS_1_REGION_SPLIT_PER_HEIGHT        2

#define HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_WIDTH        4
#define HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_HEIGHT        4

// Dynamic GOP activity TH - to tune

#define DYNAMIC_GOP_SUB_1080P_L6_VS_L5_COST_TH        11
#define DYNAMIC_GOP_SUB_1080P_L5_VS_L4_COST_TH        19
#define DYNAMIC_GOP_SUB_1080P_L4_VS_L3_COST_TH        30    // No L4_VS_L3 - 25 is the TH after 1st round of tuning

#define DYNAMIC_GOP_ABOVE_1080P_L6_VS_L5_COST_TH    15//25//5//
#define DYNAMIC_GOP_ABOVE_1080P_L5_VS_L4_COST_TH    25//28//9//
#define DYNAMIC_GOP_ABOVE_1080P_L4_VS_L3_COST_TH    30    // No L4_VS_L3 - 28 is the TH after 1st round of tuning
#define DYNAMIC_GOP_SUB_480P_L6_VS_L5_COST_TH        9

#define SB_COMPLEXITY_NON_MOVING_INDEX_TH_0 30
#define SB_COMPLEXITY_NON_MOVING_INDEX_TH_1 29
#define SB_COMPLEXITY_NON_MOVING_INDEX_TH_2 23

#define GRADUAL_LUMINOSITY_CHANGE_TH                        3
#define FADED_LCU_PERCENTAGE_TH                             10
#define FADED_PICTURES_TH                                   15
#define CLASS_SUB_0_PICTURE_ACTIVITY_REGIONS_TH             1
#define CLASS_1_SIZE_PICTURE_ACTIVITY_REGIONS_TH            2
#define HIGHER_THAN_CLASS_1_PICTURE_ACTIVITY_REGIONS_TH     8

#define IS_COMPLEX_LCU_VARIANCE_TH                          100
#define IS_COMPLEX_LCU_FLAT_VARIANCE_TH                     10
#define IS_COMPLEX_LCU_VARIANCE_DEVIATION_TH                13
#define IS_COMPLEX_LCU_ZZ_SAD_FACTOR_TH                     25

#define MAX_SUPPORTED_SEGMENTS                            7
#define NUM_QPS                                           52

// The EbAuraStatus type is used to describe the aura status
#define EbAuraStatus       uint8_t
#define AURA_STATUS_0        0
#define AURA_STATUS_1        1
#define AURA_STATUS_2        2
#define AURA_STATUS_3        3
#define INVALID_AURA_STATUS  128

// Aura detection definitions
#define    AURA_4K_DISTORTION_TH    25
#define    AURA_4K_DISTORTION_TH_6L 20

// The EB_4L_PRED_ERROR_CLASS type is used to inform about the prediction error compared to 4L
#define EB_4L_PRED_ERROR_CLASS    uint8_t
#define PRED_ERROR_CLASS_0          0
#define PRED_ERROR_CLASS_1          1
#define INVALID_PRED_ERROR_CLASS    128

#define EbScdMode uint8_t
#define SCD_MODE_0  0     // SCD OFF
#define SCD_MODE_1   1     // Light SCD (histograms generation on the 1/16 decimated input)
#define SCD_MODE_2   2     // Full SCD

#define EbBlockMeanPrec uint8_t
#define BLOCK_MEAN_PREC_FULL 0
#define BLOCK_MEAN_PREC_SUB  1

#define EbPmMode uint8_t
#define PM_MODE_0  0     // 1-stage PM
#define PM_MODE_1  1     // 2-stage PM 4K
#define PM_MODE_2  2     // 2-stage PM Sub 4K

#define EB_ZZ_SAD_MODE uint8_t
#define ZZ_SAD_MODE_0  0        // ZZ SAD on Decimated resolution
#define ZZ_SAD_MODE_1  1        // ZZ SAD on Full resolution

#define EbPfMode uint8_t
#define PF_OFF  0
#define PF_N2   1
#define PF_N4   2
#define STAGE uint8_t
#define ED_STAGE  1      // ENCDEC stage

#define EB_TRANS_COEFF_SHAPE uint8_t
#define DEFAULT_SHAPE 0
#define N2_SHAPE      1
#define N4_SHAPE      2
#define ONLY_DC_SHAPE 3

#define EB_CHROMA_LEVEL uint8_t
#define CHROMA_MODE_0  0 // Full chroma search @ MD
#define CHROMA_MODE_1  1 // Fast chroma search @ MD
#define CHROMA_MODE_2  2 // Chroma blind @ MD + CFL @ EP
#define CHROMA_MODE_3  3 // Chroma blind @ MD + no CFL @ EP

typedef enum EbSbComplexityStatus
{
    SB_COMPLEXITY_STATUS_0 = 0,
    SB_COMPLEXITY_STATUS_1 = 1,
    SB_COMPLEXITY_STATUS_2 = 2,
    SB_COMPLEXITY_STATUS_INVALID = (uint8_t)~0
} EbSbComplexityStatus;

typedef enum EbCleanUpMode
{
    CLEAN_UP_MODE_0 = 0,
    CLEAN_UP_MODE_1 = 1
} EbCleanUpMode;

typedef enum EbSaoMode
{
    SAO_MODE_0 = 0,
    SAO_MODE_1 = 1
} EbSaoMode;

typedef enum EbCu8x8Mode
{
    CU_8x8_MODE_0 = 0,  // Perform OIS, Full_Search, Fractional_Search & Bipred for CU_8x8
    CU_8x8_MODE_1 = 1   // Perform OIS and only Full_Search for CU_8x8
} EbCu8x8Mode;


#if MULTI_PASS_PD
// Multi-Pass Partitioning Depth(Multi - Pass PD) performs multiple PD stages for the same SB towards 1 final Partitioning Structure
// As we go from PDn to PDn + 1, the prediction accuracy of the MD feature(s) increases while the number of block(s) decreases
typedef enum EbPictureDepthMode
{
    PIC_MULTI_PASS_PD_MODE_0    = 0, // Multi-Pass PD Mode 0: PD0 | PD0_REFINEMENT
    PIC_MULTI_PASS_PD_MODE_1    = 1, // Multi-Pass PD Mode 1: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ only
    PIC_MULTI_PASS_PD_MODE_2    = 2, // Multi-Pass PD Mode 2: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ and SQ coeff info
    PIC_MULTI_PASS_PD_MODE_3    = 3, // Multi-Pass PD Mode 3: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ and both SQ and NSQ coeff info
    PIC_ALL_DEPTH_MODE          = 4, // ALL sq and nsq:  SB size -> 4x4
    PIC_ALL_C_DEPTH_MODE        = 5, // ALL sq and nsq with control :  SB size -> 4x4
    PIC_SQ_DEPTH_MODE           = 6, // ALL sq:  SB size -> 4x4
    PIC_SQ_NON4_DEPTH_MODE      = 7, // SQ:  SB size -> 8x8
    PIC_OPEN_LOOP_DEPTH_MODE    = 8, // Early Inter Depth Decision:  SB size -> 8x8
    PIC_SB_SWITCH_DEPTH_MODE    = 9  // Adaptive Depth Partitioning

} EbPictureDepthMode;
#else
typedef enum EbPictureDepthMode
{
    PIC_ALL_DEPTH_MODE          = 0, // ALL sq and nsq:  SB size -> 4x4
    PIC_ALL_C_DEPTH_MODE        = 1, // ALL sq and nsq with control :  SB size -> 4x4
    PIC_SQ_DEPTH_MODE           = 2, // ALL sq:  SB size -> 4x4
    PIC_SQ_NON4_DEPTH_MODE      = 3, // SQ:  SB size -> 8x8
    PIC_OPEN_LOOP_DEPTH_MODE    = 4, // Early Inter Depth Decision:  SB size -> 8x8
    PIC_SB_SWITCH_DEPTH_MODE    = 5  // Adaptive Depth Partitioning
} EbPictureDepthMode;
#endif
#define EB_SB_DEPTH_MODE              uint8_t
#define SB_SQ_BLOCKS_DEPTH_MODE             1
#define SB_SQ_NON4_BLOCKS_DEPTH_MODE        2
#define SB_OPEN_LOOP_DEPTH_MODE             3
#define SB_FAST_OPEN_LOOP_DEPTH_MODE        4
#define SB_PRED_OPEN_LOOP_DEPTH_MODE        5

typedef enum EbIntrA4x4SearchMethod
{
    INTRA4x4_OFF = 0,
    INTRA4x4_INLINE_SEARCH = 1,
    INTRA4x4_REFINEMENT_SEARCH = 2,
} EbIntrA4x4SearchMethod;

static const int32_t global_motion_threshold[MAX_HIERARCHICAL_LEVEL][MAX_TEMPORAL_LAYERS] = { // [Highest Temporal Layer] [Temporal Layer Index]
    { 2 },
    { 4, 2 },
    { 8, 4, 2 },
    { 16, 8, 4, 2 },
    { 32, 16, 8, 4, 2 },    // Derived by analogy from 4-layer settings
    { 64, 32, 16, 8, 4, 2 }
};

static const int32_t hme_level_0_search_area_multiplier_x[MAX_HIERARCHICAL_LEVEL][MAX_TEMPORAL_LAYERS] = { // [Highest Temporal Layer] [Temporal Layer Index]
    { 100 },
    { 100, 100 },
    { 100, 100, 100 },
    { 200, 140, 100,  70 },
    { 350, 200, 100, 100, 100 },
    { 525, 350, 200, 100, 100, 100 }
};

static const int32_t hme_level_0_search_area_multiplier_y[MAX_HIERARCHICAL_LEVEL][MAX_TEMPORAL_LAYERS] = { // [Highest Temporal Layer] [Temporal Layer Index]
    { 100 },
    { 100, 100 },
    { 100, 100, 100 },
    { 200, 140, 100, 70 },
    { 350, 200, 100, 100, 100 },
    { 525, 350, 200, 100, 100, 100 }
};

typedef enum RasterScanCuIndex
{
    // 2Nx2N [85 partitions]
    RASTER_SCAN_CU_INDEX_64x64 = 0,
    RASTER_SCAN_CU_INDEX_32x32_0 = 1,
    RASTER_SCAN_CU_INDEX_32x32_1 = 2,
    RASTER_SCAN_CU_INDEX_32x32_2 = 3,
    RASTER_SCAN_CU_INDEX_32x32_3 = 4,
    RASTER_SCAN_CU_INDEX_16x16_0 = 5,
    RASTER_SCAN_CU_INDEX_16x16_1 = 6,
    RASTER_SCAN_CU_INDEX_16x16_2 = 7,
    RASTER_SCAN_CU_INDEX_16x16_3 = 8,
    RASTER_SCAN_CU_INDEX_16x16_4 = 9,
    RASTER_SCAN_CU_INDEX_16x16_5 = 10,
    RASTER_SCAN_CU_INDEX_16x16_6 = 11,
    RASTER_SCAN_CU_INDEX_16x16_7 = 12,
    RASTER_SCAN_CU_INDEX_16x16_8 = 13,
    RASTER_SCAN_CU_INDEX_16x16_9 = 14,
    RASTER_SCAN_CU_INDEX_16x16_10 = 15,
    RASTER_SCAN_CU_INDEX_16x16_11 = 16,
    RASTER_SCAN_CU_INDEX_16x16_12 = 17,
    RASTER_SCAN_CU_INDEX_16x16_13 = 18,
    RASTER_SCAN_CU_INDEX_16x16_14 = 19,
    RASTER_SCAN_CU_INDEX_16x16_15 = 20,
    RASTER_SCAN_CU_INDEX_8x8_0 = 21,
    RASTER_SCAN_CU_INDEX_8x8_1 = 22,
    RASTER_SCAN_CU_INDEX_8x8_2 = 23,
    RASTER_SCAN_CU_INDEX_8x8_3 = 24,
    RASTER_SCAN_CU_INDEX_8x8_4 = 25,
    RASTER_SCAN_CU_INDEX_8x8_5 = 26,
    RASTER_SCAN_CU_INDEX_8x8_6 = 27,
    RASTER_SCAN_CU_INDEX_8x8_7 = 28,
    RASTER_SCAN_CU_INDEX_8x8_8 = 29,
    RASTER_SCAN_CU_INDEX_8x8_9 = 30,
    RASTER_SCAN_CU_INDEX_8x8_10 = 31,
    RASTER_SCAN_CU_INDEX_8x8_11 = 32,
    RASTER_SCAN_CU_INDEX_8x8_12 = 33,
    RASTER_SCAN_CU_INDEX_8x8_13 = 34,
    RASTER_SCAN_CU_INDEX_8x8_14 = 35,
    RASTER_SCAN_CU_INDEX_8x8_15 = 36,
    RASTER_SCAN_CU_INDEX_8x8_16 = 37,
    RASTER_SCAN_CU_INDEX_8x8_17 = 38,
    RASTER_SCAN_CU_INDEX_8x8_18 = 39,
    RASTER_SCAN_CU_INDEX_8x8_19 = 40,
    RASTER_SCAN_CU_INDEX_8x8_20 = 41,
    RASTER_SCAN_CU_INDEX_8x8_21 = 42,
    RASTER_SCAN_CU_INDEX_8x8_22 = 43,
    RASTER_SCAN_CU_INDEX_8x8_23 = 44,
    RASTER_SCAN_CU_INDEX_8x8_24 = 45,
    RASTER_SCAN_CU_INDEX_8x8_25 = 46,
    RASTER_SCAN_CU_INDEX_8x8_26 = 47,
    RASTER_SCAN_CU_INDEX_8x8_27 = 48,
    RASTER_SCAN_CU_INDEX_8x8_28 = 49,
    RASTER_SCAN_CU_INDEX_8x8_29 = 50,
    RASTER_SCAN_CU_INDEX_8x8_30 = 51,
    RASTER_SCAN_CU_INDEX_8x8_31 = 52,
    RASTER_SCAN_CU_INDEX_8x8_32 = 53,
    RASTER_SCAN_CU_INDEX_8x8_33 = 54,
    RASTER_SCAN_CU_INDEX_8x8_34 = 55,
    RASTER_SCAN_CU_INDEX_8x8_35 = 56,
    RASTER_SCAN_CU_INDEX_8x8_36 = 57,
    RASTER_SCAN_CU_INDEX_8x8_37 = 58,
    RASTER_SCAN_CU_INDEX_8x8_38 = 59,
    RASTER_SCAN_CU_INDEX_8x8_39 = 60,
    RASTER_SCAN_CU_INDEX_8x8_40 = 61,
    RASTER_SCAN_CU_INDEX_8x8_41 = 62,
    RASTER_SCAN_CU_INDEX_8x8_42 = 63,
    RASTER_SCAN_CU_INDEX_8x8_43 = 64,
    RASTER_SCAN_CU_INDEX_8x8_44 = 65,
    RASTER_SCAN_CU_INDEX_8x8_45 = 66,
    RASTER_SCAN_CU_INDEX_8x8_46 = 67,
    RASTER_SCAN_CU_INDEX_8x8_47 = 68,
    RASTER_SCAN_CU_INDEX_8x8_48 = 69,
    RASTER_SCAN_CU_INDEX_8x8_49 = 70,
    RASTER_SCAN_CU_INDEX_8x8_50 = 71,
    RASTER_SCAN_CU_INDEX_8x8_51 = 72,
    RASTER_SCAN_CU_INDEX_8x8_52 = 73,
    RASTER_SCAN_CU_INDEX_8x8_53 = 74,
    RASTER_SCAN_CU_INDEX_8x8_54 = 75,
    RASTER_SCAN_CU_INDEX_8x8_55 = 76,
    RASTER_SCAN_CU_INDEX_8x8_56 = 77,
    RASTER_SCAN_CU_INDEX_8x8_57 = 78,
    RASTER_SCAN_CU_INDEX_8x8_58 = 79,
    RASTER_SCAN_CU_INDEX_8x8_59 = 80,
    RASTER_SCAN_CU_INDEX_8x8_60 = 81,
    RASTER_SCAN_CU_INDEX_8x8_61 = 82,
    RASTER_SCAN_CU_INDEX_8x8_62 = 83,
    RASTER_SCAN_CU_INDEX_8x8_63 = 84
} RasterScanCuIndex;

static const uint32_t raster_scan_cu_x[CU_MAX_COUNT] =
{
    0,
    0, 32,
    0, 32,
    0, 16, 32, 48,
    0, 16, 32, 48,
    0, 16, 32, 48,
    0, 16, 32, 48,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56
};

static const uint32_t raster_scan_cu_y[CU_MAX_COUNT] =
{
    0,
    0, 0,
    32, 32,
    0, 0, 0, 0,
    16, 16, 16, 16,
    32, 32, 32, 32,
    48, 48, 48, 48,
    0, 0, 0, 0, 0, 0, 0, 0,
    8, 8, 8, 8, 8, 8, 8, 8,
    16, 16, 16, 16, 16, 16, 16, 16,
    24, 24, 24, 24, 24, 24, 24, 24,
    32, 32, 32, 32, 32, 32, 32, 32,
    40, 40, 40, 40, 40, 40, 40, 40,
    48, 48, 48, 48, 48, 48, 48, 48,
    56, 56, 56, 56, 56, 56, 56, 56
};

static const uint32_t raster_scan_cu_size[CU_MAX_COUNT] =
{   64,
    32, 32,
    32, 32,
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8
};

static const uint32_t RASTER_SCAN_CU_DEPTH[CU_MAX_COUNT] =
{   0,
    1, 1,
    1, 1,
    2, 2, 2, 2,
    2, 2, 2, 2,
    2, 2, 2, 2,
    2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3
};

static const uint32_t raster_scan_to_md_scan[CU_MAX_COUNT] =
{
    0,
    1, 22,
    43, 64,
    2, 7, 23, 28,
    12, 17, 33, 38,
    44, 49, 65, 70,
    54, 59, 75, 80,
    3, 4, 8, 9, 24, 25, 29, 30,
    5, 6, 10, 11, 26, 27, 31, 32,
    13, 14, 18, 19, 34, 35, 39, 40,
    15, 16, 20, 21, 36, 37, 41, 42,
    45, 46, 50, 51, 66, 67, 71, 72,
    47, 48, 52, 53, 68, 69, 73, 74,
    55, 56, 60, 61, 76, 77, 81, 82,
    57, 58, 62, 63, 78, 79, 83, 84
};

static const uint32_t ParentBlockIndex[85] = { 0, 0, 0, 2, 2, 2, 2, 0, 7, 7, 7, 7, 0, 12, 12, 12, 12, 0, 17, 17, 17, 17, 0, 0,
    23, 23, 23, 23, 0, 28, 28, 28, 28, 0, 33, 33, 33, 33, 0, 38, 38, 38, 38, 0, 0,
    44, 44, 44, 44, 0, 49, 49, 49, 49, 0, 54, 54, 54, 54, 0, 59, 59, 59, 59, 0, 0,
    65, 65, 65, 65, 0, 70, 70, 70, 70, 0, 75, 75, 75, 75, 0, 80, 80, 80, 80 };

static const uint32_t md_scan_to_raster_scan[CU_MAX_COUNT] =
{
    0,
    1,
    5, 21, 22, 29, 30,
    6, 23, 24, 31, 32,
    9, 37, 38, 45, 46,
    10, 39, 40, 47, 48,
    2,
    7, 25, 26, 33, 34,
    8, 27, 28, 35, 36,
    11, 41, 42, 49, 50,
    12, 43, 44, 51, 52,
    3,
    13, 53, 54, 61, 62,
    14, 55, 56, 63, 64,
    17, 69, 70, 77, 78,
    18, 71, 72, 79, 80,
    4,
    15, 57, 58, 65, 66,
    16, 59, 60, 67, 68,
    19, 73, 74, 81, 82,
    20, 75, 76, 83, 84
};

static const uint32_t raster_scan_cu_parent_index[CU_MAX_COUNT] =
{   0,
    0, 0,
    0, 0,
    1, 1, 2, 2,
    1, 1, 2, 2,
    3, 3, 4, 4,
    3, 3, 4, 4,
    5, 5, 6, 6, 7, 7, 8, 8,
    5, 5, 6, 6, 7, 7, 8, 8,
    9, 9, 10, 10, 11, 11, 12, 12,
    9, 9, 10, 10, 11, 11, 12, 12,
    13, 13, 14, 14, 15, 15, 16, 16,
    13, 13, 14, 14, 15, 15, 16, 16,
    17, 17, 18, 18, 19, 19, 20, 20,
    17, 17, 18, 18, 19, 19, 20, 20
};

#define UNCOMPRESS_SAD(x) ( ((x) & 0x1FFF)<<(((x)>>13) & 7) )

static const uint32_t MD_SCAN_TO_OIS_32x32_SCAN[CU_MAX_COUNT] =
{
    /*0  */0,
    /*1  */0,
    /*2  */0,
    /*3  */0,
    /*4  */0,
    /*5  */0,
    /*6  */0,
    /*7  */0,
    /*8  */0,
    /*9  */0,
    /*10 */0,
    /*11 */0,
    /*12 */0,
    /*13 */0,
    /*14 */0,
    /*15 */0,
    /*16 */0,
    /*17 */0,
    /*18 */0,
    /*19 */0,
    /*20 */0,
    /*21 */0,
    /*22 */1,
    /*23 */1,
    /*24 */1,
    /*25 */1,
    /*26 */1,
    /*27 */1,
    /*28 */1,
    /*29 */1,
    /*30 */1,
    /*31 */1,
    /*32 */1,
    /*33 */1,
    /*34 */1,
    /*35 */1,
    /*36 */1,
    /*37 */1,
    /*38 */1,
    /*39 */1,
    /*40 */1,
    /*41 */1,
    /*42 */1,
    /*43 */2,
    /*44 */2,
    /*45 */2,
    /*46 */2,
    /*47 */2,
    /*48 */2,
    /*49 */2,
    /*50 */2,
    /*51 */2,
    /*52 */2,
    /*53 */2,
    /*54 */2,
    /*55 */2,
    /*56 */2,
    /*57 */2,
    /*58 */2,
    /*59 */2,
    /*60 */2,
    /*61 */2,
    /*62 */2,
    /*63 */2,
    /*64 */3,
    /*65 */3,
    /*66 */3,
    /*67 */3,
    /*68 */3,
    /*69 */3,
    /*70 */3,
    /*71 */3,
    /*72 */3,
    /*73 */3,
    /*74 */3,
    /*75 */3,
    /*76 */3,
    /*77 */3,
    /*78 */3,
    /*79 */3,
    /*80 */3,
    /*81 */3,
    /*82 */3,
    /*83 */3,
    /*84 */3,
};

#if TWO_PASS
typedef struct stat_struct_t
{
    uint32_t                        referenced_area[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
} stat_struct_t;
#if TWO_PASS_IMPROVEMENT
#define TWO_PASS_IR_THRSHLD 40  // Intra refresh threshold used to reduce the reference area.
                                // If the periodic Intra refresh is less than the threshold,
                                // the referenced area is normalized
#endif
#endif
#define SC_MAX_LEVEL 2 // 2 sets of HME/ME settings are used depending on the scene content mode

/******************************************************************************
                            ME/HME settings
*******************************************************************************/
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    },{
        {   0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    }
};
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_level0_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    }
};

static const uint16_t hme_level0_total_search_area_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        {  96,  96,    96,   96,  112,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 },
     } , {
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 }
    }
};

static const uint16_t hme_level0_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  48,   48,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
    }
};
static const uint16_t hme_level0_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
    }
};
static const uint16_t hme_level0_total_search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    } , {
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 }
    }
};
static const uint16_t hme_level0_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    } , {
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
    }
};
static const uint16_t hme_level0_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }, {
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
    }
};

// HME LEVEL 1
   //      M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_level1_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }, {
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }
};
static const uint16_t hme_level1_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
    }
};
static const uint16_t hme_level1_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
    }
};
static const uint16_t hme_level1_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
    }
};
static const uint16_t hme_level1_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
    }
};
// HME LEVEL 2
    //     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_level2_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }
};
static const uint16_t hme_level2_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t hme_level2_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t hme_level2_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t hme_level2_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};

static const uint16_t search_area_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        { 128,  128,  128,  128,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 160,  160,  160,  160,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 192,  192,  192,  192,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 192,  192,  192,  192,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
    } , {
#if M0_OPT
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {640 ,  640,  640,  288,  288,  168,  128,  128,   64,   80,   80,    80,   80 },
        {640 ,  640,  640,  288,  288,  168,  128,  128,   64,   80,   80,    80,   80 }
#else
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {960 ,  640,  640,  288,  288,  168,  128,  128,   64,   80,   80,    80,   80 },
        {960 ,  640,  640,  288,  288,  168,  128,  128,   64,   80,   80,    80,   80 }
#endif
    }
};
static const uint16_t search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        { 128,  128,  128,  128,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 160,  160,  160,  160,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 192,  192,  192,  192,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 192,  192,  192,  192,   32,   32,   32,   32,   16,    9,    9,     9,    9 }
    } , {
#if M0_OPT
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {640 ,  640,  640,  288,  288,  168,   80,   80,   48,   80,   80,    80,   80 },
        {640 ,  640,  640,  288,  288,  168,   80,   80,   48,   80,   80,    80,   80 }
#else
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {480 ,  480,  480,  144,  144,   88,   48,   48,   48,   48,   48,    48,   48 },
        {960 ,  640,  640,  288,  288,  168,   80,   80,   48,   80,   80,    80,   80 },
        {960 ,  640,  640,  288,  288,  168,   80,   80,   48,   80,   80,    80,   80 }
#endif
    }

    //     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
};

/******************************************************************************
                            ME/HME settings for Altref Temporal Filtering
*******************************************************************************/
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    }
};
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_level0_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_4K_RANGE
    }
};

static const uint16_t tf_hme_level0_total_search_area_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  112,  112,  112,  112,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 },
     } , {
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  112,  112,  112,  112,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 },
    }
};

static const uint16_t tf_hme_level0_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    }
};
static const uint16_t tf_hme_level0_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    }
};
static const uint16_t tf_hme_level0_total_search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    } , {
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    }
};
static const uint16_t tf_hme_level0_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    } , {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }
};
static const uint16_t tf_hme_level0_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }, {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }
};

// HME LEVEL 1
   //      M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_level1_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }, {
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }
};
static const uint16_t tf_hme_level1_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    }
};
static const uint16_t tf_hme_level1_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    }
};
static const uint16_t tf_hme_level1_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    }
};
static const uint16_t tf_hme_level1_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    }
};
// HME LEVEL 2
    //     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_level2_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }
};
static const uint16_t tf_hme_level2_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t tf_hme_level2_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t tf_hme_level2_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};
static const uint16_t tf_hme_level2_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    }
};

static const uint16_t tf_search_area_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  64,   64,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 112,  112,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 128,  128,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 128,  128,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 }
    } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 112,  112,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 128,  128,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 },
        { 128,  128,   64,   64,   64,   64,   64,   64,   48,   16,   16,    16,   16 }
    }
};
static const uint16_t tf_search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {  64,   64,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 112,  112,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 128,  128,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 128,  128,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 }
    } , {
        {  64,   64,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 112,  112,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 128,  128,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 },
        { 128,  128,   64,   64,   32,   32,   32,   32,   16,    9,    9,     9,    9 }
    }

    //     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
};

static const uint16_t ep_to_pa_block_index[BLOCK_MAX_COUNT_SB_64] = {
    0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    2 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    3 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    5 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    6 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    7 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    8 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    9 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    10,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    11,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    12,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    13,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    14,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    15,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    16,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    17,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    18,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    19,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    20,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    21,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    22,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    23,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    24,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    25,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    26,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    27,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    28,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    29,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    30,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    31,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    32,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    33,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    34,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    35,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    36,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    37,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    38,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    39,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    40,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    41,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    42,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    43,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    44,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    45,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    46,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    47,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    48,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    49,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    50,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    51,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    52,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    53,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    54,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    55,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    56,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    57,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    58,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    59,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    60,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    61,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    62,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    63,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    64,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    65,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    66,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    67,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    68,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    69,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    70,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    71,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    72,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    73,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    74,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    75,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 , 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    76,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    77,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    78,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    79,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    80,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    81,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    82,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    83,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
    84,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
};
#ifdef __cplusplus
}
#endif
#endif // EbDefinitions_h
/* File EOF */
