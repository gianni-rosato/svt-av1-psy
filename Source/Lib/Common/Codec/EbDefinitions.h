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
//
#ifndef EbDefinitions_h
#define EbDefinitions_h
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>
#include "EbSvtAv1.h"
#include "EbSvtAv1Enc.h"
#ifdef _WIN32
#define inline __inline
#elif __GNUC__
#define inline __inline__
#else
#define inline
#endif

#if defined(_MSC_VER)
#if defined(_M_IX86) || defined(_M_X64)
#define ARCH_X86
#endif
#endif

#if __GNUC__
#if defined(__i386__) || defined(__x86_64__)
#define ARCH_X86
#endif
#endif
#ifdef __cplusplus
extern "C" {
#endif

#ifndef NON_AVX512_SUPPORT
#define NON_AVX512_SUPPORT
#endif

#define MASTER_SYNCH   1  // based on master d9b5494b9d68df19c5d8ceb62d203cdf9b1e1b5c


///////// START MASTER_SYNCH

#if MASTER_SYNCH
#define FIX_WARNINGS                1
#define FIXED_SQ_WEIGHT_PER_QP      1
#define SQ_WEIGHT_PATCH_0           1
#define MAR2_M8_ADOPTIONS           1
#define MAR2_M7_ADOPTIONS           1
#define MAR3_M2_ADOPTIONS           1
#define MAR3_M6_ADOPTIONS           1
#define IBC_SYNCH                   1
#define MRP_CTRL                    1  //add control to inject smaller number of references.
#define MAR4_M8_ADOPTIONS           1
#define MAR4_M3_ADOPTIONS           1
#define PALETTE_SYNCH               1
#define MAR4_M6_ADOPTIONS           1
#define SHUT_ME_DISTORTION          1 //Removed the ME distortions (209 elements), and the HEVC-legacy early inter-depth decision.
#define MULTI_STAGE_HME             1
#if MULTI_STAGE_HME
#define DISABLE_HME_PRE_CHECK             1
#define ENABLE_HME_AT_INC_SB              1
#define NEW_HME_DISTANCE_ALGORITHM        1
#define DISABLE_HME_OF_SAME_POC           1
#define DISABLE_HME_L0_FOR_240P           1
#endif
#define HME_PRUNE_BUG_FIX                 1
#define MAR10_ADOPTIONS                   1 // Adoptions for all presets
#define NICS_SYNCH                        1
#define MAR11_ADOPTIONS                   1 // Adoptions for M2, M3, M4, M5
#define WEDGE_SYNCH                       1
#define DEPTH_PART_CLEAN_UP               1 // sb_128x128 if NSC, sb_64x64 if SC, and multi-pass PD till M8

#define REMOVE_COMBINE_CLASS12            1 // remove code associated with combine_class12 feature
#define REMOVE_OLD_NICS                   1 // Remove code for old NICS levels
#define ADD_ME_SIGNAL_FOR_PRUNING_TH      1 // Add signals for mode-dependent ME thresholds
#define ADD_HME_MIN_MAX_MULTIPLIER_SIGNAL 1 // Add ME signal for the max HME search area multiplier
#define MAR12_M8_ADOPTIONS                1
#define MAR12_ADOPTIONS                   1 // Adoptions for all modes
#define REMOVED_MEM_OPT_CDF               1
#define FIX_MR_PD1                        1 // Disable PD1 refinement changes for MR.
#define PME_SORT_REF                      1 //add reference sorting of pme results
#define OBMC_FAST                         1 //faster obmc mode (3). cleaner obmc signalling.
#define REMOVE_MD_EXIT                    1 // remove md_exit
#define MAR16_M8_ADOPTIONS                1 // M8 adoption for TH value
#define ADDED_CFL_OFF                     1
#define ADOPT_CHROMA_MODE1_CFL_OFF        1
#define PIC_BASED_RE_OFF                  0 //-------------------------------------------------------
#define MR_MODE_FOR_PIC_MULTI_PASS_PD_MODE_1 1 // shut SQ vs. NSQ if MR (for multi_pass_pd_level = PIC_MULTI_PASS_PD_MODE_1 or PIC_MULTI_PASS_PD_MODE_2 or PIC_MULTI_PASS_PD_MODE_3)
#define ADD_SAD_AT_PME_SIGNAL      1 // Add signal for using SAD at PME
#define MAR17_ADOPTIONS            1 // Push features with bad slope to M8 & beyond.
    // it3
#define M5_CHROMA_NICS             1
#define INTER_COMP_REDESIGN        1 // new fast mode, cleaner signaling and code
#define MAR18_MR_TESTS_ADOPTIONS   1 // adoptions for MR, M0, and M2
#define SKIP_DEPTH_SYNCH           1
#define MAR18_ADOPTIONS            1 // adoptions in M5/M8
#define REU_UPDATE                 1 // use top right instead of top SB for CDF calculation
#define ADD_NEW_MPPD_LEVEL         1 // add a new MPPD level with PD0 | PD1 | PD2 w/o sq/nsq decision
#define LOG_MV_VALIDITY            1 //report error message if MV is beyond av1 limits
#define MD_CFL                     1 // Modified cfl search in MD
#if MD_CFL
#define CFL_REDUCED_ALPHA    1 // use faster libaom_short_cuts_ths
#endif
#define UV_SEARCH_MODE_INJCECTION  1 // use the luma mode ijection method in chroma independent mode search
#define MAR19_ADOPTIONS            1 // Adoptions for all modes
#define MAR20_M4_ADOPTIONS         1 // Adoptions in M4
#define ADOPT_SQ_ME_SEARCH_AREA    1 // Adopt a square search area for ME (all modes)
#define MAR20_ADOPTIONS            1 // Adoptions affecting all modes
#define MD_CONFIG_SB               0
#endif

///////// END MASTER_SYNCH

#define MAX_TILE_CNTS 128 // Annex A.3
        // Add support for other prediction structure
#define MR_MODE 0

#define ALT_REF_QP_THRESH 20
#define HIGH_PRECISION_MV_QTHRESH 150

// Actions in the second pass: Frame and SB QP assignment and temporal filtering strenght change
//FOR DEBUGGING - Do not remove
#define NO_ENCDEC         0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define AOM_INTERP_EXTEND 4
#define AOM_LEFT_TOP_MARGIN_PX(subsampling) \
    ((AOM_BORDER_IN_PIXELS >> subsampling) - AOM_INTERP_EXTEND)
#define AOM_LEFT_TOP_MARGIN_SCALED(subsampling) \
    (AOM_LEFT_TOP_MARGIN_PX(subsampling) << SCALE_SUBPEL_BITS)

#define H_PEL_SEARCH_WIND_3 3 // 1/2-pel serach window 3
#define H_PEL_SEARCH_WIND_2 2 // 1/2-pel serach window 2
#define Q_PEL_SEARCH_WIND 2 // 1/4-pel search window
#define HP_REF_OPT 1 // Remove redundant positions.

#define ENABLE_PME_SAD 0
#define SWITCH_XY_LOOPS_PME_SAD_SSD 0
#if SWITCH_XY_LOOPS_PME_SAD_SSD
#define RESTRUCTURE_SAD 1
#endif

#define FIX_HBD_R2R 1 // Fix 10bit error in over-boundaries CUs (incomplete SB)

typedef enum MeHpMode {
    EX_HP_MODE         = 0, // Exhaustive  1/2-pel serach mode.
    REFINEMENT_HP_MODE = 1,// Refinement 1/2-pel serach mode.
    SWITCHABLE_HP_MODE = 2 // Switch between EX_HP_MODE and REFINEMENT_HP_MODE mode.
} MeHpMode;
typedef enum GM_LEVEL {
    GM_FULL      = 0, // Exhaustive search mode.
    GM_DOWN      = 1, // Downsampled search mode, with a downsampling factor of 2 in each dimension
    GM_TRAN_ONLY = 2 // Translation only using ME MV.
} GM_LEVEL;
typedef enum SqWeightOffsets {
    CONSERVATIVE_OFFSET_0 = 5,
    CONSERVATIVE_OFFSET_1 = 10,
    AGGRESSIVE_OFFSET_0   = -5,
    AGGRESSIVE_OFFSET_1   = -10
} SqWeightOffsets;

#if !FIXED_SQ_WEIGHT_PER_QP
typedef enum Qp { QP_20 = 20, QP_32 = 32, QP_43 = 43, QP_55 = 55, QP_63 = 63 } Qp;
#endif
struct Buf2D {
    uint8_t *buf;
    uint8_t *buf0;
    int      width;
    int      height;
    int      stride;
};
typedef struct MvLimits {
    int col_min;
    int col_max;
    int row_min;
    int row_max;
} MvLimits;
typedef struct {
    uint8_t by;
    uint8_t bx;
    uint8_t skip;
} CdefList;

/*!\brief force enum to be unsigned 1 byte*/
#define UENUM1BYTE(enumvar) \
    ;                       \
    typedef uint8_t enumvar

enum {
    DIAMOND      = 0,
    NSTEP        = 1,
    HEX          = 2,
    BIGDIA       = 3,
    SQUARE       = 4,
    FAST_HEX     = 5,
    FAST_DIAMOND = 6
} UENUM1BYTE(SEARCH_METHODS);

/********************************************************/
/****************** Pre-defined Values ******************/
/********************************************************/

/* maximum number of frames allowed for the Alt-ref picture computation
 * this number can be increased by increasing the constant
 * FUTURE_WINDOW_WIDTH defined in EbPictureDecisionProcess.c
 */
#define ALTREF_MAX_NFRAMES 10
#define ALTREF_MAX_STRENGTH 6
#define PAD_VALUE (128 + 32)
#define PAD_VALUE_SCALED (128 + 128 + 32)
#define NSQ_TAB_SIZE 8
#define NUMBER_OF_DEPTH 6
#define NUMBER_OF_SHAPES 10
//  Delta QP support
#define ADD_DELTA_QP_SUPPORT 1 // Add delta QP support
#define BLOCK_MAX_COUNT_SB_128 4421
#define BLOCK_MAX_COUNT_SB_64 1101
#define MAX_TXB_COUNT 16 // Maximum number of transform blocks per depth
#if REMOVE_MR_MACRO
#define MAX_NFL 250 // Maximum number of candidates MD can support
#else
#if MAR12_ADOPTIONS && MR_MODE
#define MAX_NFL 250 // Maximum number of candidates MD can support
#else
#if M0_NIC
#define MAX_NFL 150 // Maximum number of candidates MD can support
#else
#define MAX_NFL 125 // Maximum number of candidates MD can support
#endif
#endif
#endif
#define MAX_NFL_BUFF_Y \
    (MAX_NFL + CAND_CLASS_TOTAL) //need one extra temp buffer for each fast loop call
#define MAX_NFL_BUFF (MAX_NFL_BUFF_Y + 84) //need one extra temp buffer for each fast loop call
#define MAX_LAD 120 // max lookahead-distance 2x60fps
#define ROUND_UV(x) (((x) >> 3) << 3)
#define AV1_PROB_COST_SHIFT 9
#define AOMINNERBORDERINPIXELS 160
#define SWITCHABLE_FILTER_CONTEXTS ((SWITCHABLE_FILTERS + 1) * 4)
#define MAX_MB_PLANE 3
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
#define SB_STRIDE_Y MAX_SB_SIZE
#define SB_STRIDE_UV (MAX_SB_SIZE >> 1)

// Min superblock size
#define MIN_SB_SIZE 64
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

#define MAX_VARTX_DEPTH 2
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
#define TX_PAD_2D ((MAX_TX_SIZE + TX_PAD_HOR) * (MAX_TX_SIZE + TX_PAD_VER) + TX_PAD_END)
#define COMPOUND_WEIGHT_MODE DIST
#define DIST_PRECISION_BITS 4
#define DIST_PRECISION (1 << DIST_PRECISION_BITS) // 16

#define LOOP_FILTER_BITMASK 0
#define PROFILE_BITS 3

// AV1 Loop Filter
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
#ifdef ARCH_X86
extern void RunEmms();
#define aom_clear_system_state() RunEmms()
#endif

/* Shift down with rounding for use when n >= 0, value >= 0 */
#define ROUND_POWER_OF_TWO(value, n) (((value) + (((1 << (n)) >> 1))) >> (n))

/* Shift down with rounding for signed integers, for use when n >= 0 */
#define ROUND_POWER_OF_TWO_SIGNED(value, n) \
    (((value) < 0) ? -ROUND_POWER_OF_TWO(-(value), (n)) : ROUND_POWER_OF_TWO((value), (n)))

/* Shift down with rounding for use when n >= 0, value >= 0 for (64 bit) */
#define ROUND_POWER_OF_TWO_64(value, n) (((value) + ((((int64_t)1 << (n)) >> 1))) >> (n))

/* Shift down with rounding for signed integers, for use when n >= 0 (64 bit) */
#define ROUND_POWER_OF_TWO_SIGNED_64(value, n) \
    (((value) < 0) ? -ROUND_POWER_OF_TWO_64(-(value), (n)) : ROUND_POWER_OF_TWO_64((value), (n)))

#define IS_POWER_OF_TWO(x) (((x) & ((x)-1)) == 0)

#ifdef __cplusplus
#define EB_EXTERN extern "C"
#else
#define EB_EXTERN
#endif // __cplusplus

#define INLINE __inline
#define RESTRICT
#ifdef _WIN32
#define FOPEN(f, s, m) fopen_s(&f, s, m)
#else
#define FOPEN(f, s, m) f = fopen(s, m)
#endif

#define IMPLIES(a, b) (!(a) || (b)) //  Logical 'a implies b' (or 'a -> b')
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
#define AOM_INLINE inline
#endif

#define SIMD_INLINE static AOM_FORCE_INLINE

//*********************************************************************************************************************//
// mem.h
/* shift right or left depending on sign of n */
#define RIGHT_SIGNED_SHIFT(value, n) ((n) < 0 ? ((value) << (-(n))) : ((value) >> (n)))
//*********************************************************************************************************************//
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

static __inline void mem_put_le24(void *vmem, MEM_VALUE_T val) {
    MAU_T *mem = (MAU_T *)vmem;

    mem[0] = (MAU_T)((val >> 0) & 0xff);
    mem[1] = (MAU_T)((val >> 8) & 0xff);
    mem[2] = (MAU_T)((val >> 16) & 0xff);
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

typedef struct ConvolveParams {
    int32_t      ref;
    int32_t      do_average;
    ConvBufType *dst;
    int32_t      dst_stride;
    int32_t      round_0;
    int32_t      round_1;
    int32_t      plane;
    int32_t      is_compound;
    int32_t      use_jnt_comp_avg;
    int32_t      fwd_offset;
    int32_t      bck_offset;
    int32_t      use_dist_wtd_comp_avg;
} ConvolveParams;

// texture component type
typedef enum ATTRIBUTE_PACKED {
    COMPONENT_LUMA      = 0, // luma
    COMPONENT_CHROMA    = 1, // chroma (Cb+Cr)
    COMPONENT_CHROMA_CB = 2, // chroma Cb
    COMPONENT_CHROMA_CR = 3, // chroma Cr
    COMPONENT_ALL       = 4, // Y+Cb+Cr
    COMPONENT_NONE      = 15
} COMPONENT_TYPE;

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

static INLINE unsigned int negative_to_zero(int value) { return (value < 0) ? 0 : value; }

static INLINE int av1_num_planes(EbColorConfig *color_info) {
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
typedef enum PdPass {
    PD_PASS_0,
    PD_PASS_1,
    PD_PASS_2,
    PD_PASS_TOTAL,
} PdPass;
typedef enum CandClass {
    CAND_CLASS_0,
    CAND_CLASS_1,
    CAND_CLASS_2,
    CAND_CLASS_3,
    CAND_CLASS_4,
    CAND_CLASS_5,
    CAND_CLASS_6,
    CAND_CLASS_7,
    CAND_CLASS_8,
    CAND_CLASS_TOTAL
} CandClass;

typedef enum MdStage { MD_STAGE_0, MD_STAGE_1, MD_STAGE_2, MD_STAGE_3, MD_STAGE_TOTAL } MdStage;

typedef enum MdStagingMode {
    MD_STAGING_MODE_0,
    MD_STAGING_MODE_1,
    MD_STAGING_MODE_2,
    MD_STAGING_MODE_TOTAL
} MdStagingMode;

// NICS
#define MAX_FRAME_TYPE 3 // Max number of frame type allowed for nics
#define ALL_S0 -1 // Allow all candidates from stage0
#define BEST_CANDIDATE_COUNT 4
#define MAX_REF_TYPE_CAND 30
#define PRUNE_REC_TH 5
#define PRUNE_REF_ME_TH 2
typedef enum {
    EIGHTTAP_REGULAR,
    EIGHTTAP_SMOOTH,
    MULTITAP_SHARP,
    BILINEAR,
    INTERP_FILTERS_ALL,
    SWITCHABLE_FILTERS = BILINEAR,
    SWITCHABLE         = SWITCHABLE_FILTERS + 1, /* the last switchable one */
    EXTRA_FILTERS      = INTERP_FILTERS_ALL - SWITCHABLE_FILTERS,
} InterpFilter;

#define AV1_COMMON Av1Common
enum {
    USE_2_TAPS_ORIG = 0, // This is used in temporal filtering.
    USE_2_TAPS,
    USE_4_TAPS,
    USE_8_TAPS,
} UENUM1BYTE(SUBPEL_SEARCH_TYPE);

typedef struct InterpFilterParams {
    const int16_t *filter_ptr;
    uint16_t       taps;
    uint16_t       subpel_shifts;
    InterpFilter   interp_filter;
} InterpFilterParams;

typedef enum TxSearchLevel {
    TX_SEARCH_OFF,
    TX_SEARCH_ENC_DEC,
    TX_SEARCH_INTER_DEPTH,
    TX_SEARCH_FULL_LOOP
} TxSearchLevel;

typedef enum InterpolationSearchLevel {
    IT_SEARCH_OFF,
    IT_SEARCH_FAST_LOOP_UV_BLIND,
    IT_SEARCH_FAST_LOOP,
} InterpolationSearchLevel;

typedef enum NsqSearchLevel {
    NSQ_SEARCH_OFF,
    NSQ_SEARCH_LEVEL1,
    NSQ_SEARCH_LEVEL2,
    NSQ_SEARCH_LEVEL3,
    NSQ_SEARCH_LEVEL4,
    NSQ_SEARCH_LEVEL5,
    NSQ_SEARCH_LEVEL6,
    NSQ_SEARCH_FULL
} NsqSearchLevel;

#define MAX_PARENT_SQ 6
typedef enum CompoundDistWeightMode {
    DIST,
} CompoundDistWeightMode;

// Profile 0.  8-bit and 10-bit 4:2:0 and 4:0:0 only.
// Profile 1.  8-bit and 10-bit 4:4:4
// Profile 2.  8-bit and 10-bit 4:2:2
//            12 bit  4:0:0, 4:2:2 and 4:4:4
typedef enum BitstreamProfile { PROFILE_0, PROFILE_1, PROFILE_2, MAX_PROFILES } BitstreamProfile;
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
    BlockSizeS    = BLOCK_4X16,
    BLOCK_INVALID = 255,
    BLOCK_LARGEST = (BlockSizeS - 1)
} BlockSize;

typedef enum ATTRIBUTE_PACKED {
    PARTITION_NONE,
    PARTITION_HORZ,
    PARTITION_VERT,
    PARTITION_SPLIT,
    PARTITION_HORZ_A, // HORZ split and the top partition is split again
    PARTITION_HORZ_B, // HORZ split and the bottom partition is split again
    PARTITION_VERT_A, // VERT split and the left partition is split again
    PARTITION_VERT_B, // VERT split and the right partition is split again
    PARTITION_HORZ_4, // 4:1 horizontal partition
    PARTITION_VERT_4, // 4:1 vertical partition
    EXT_PARTITION_TYPES,
    PARTITION_TYPES   = PARTITION_SPLIT + 1,
    PARTITION_INVALID = 255
} PartitionType;

#define MAX_NUM_BLOCKS_ALLOC 7493 //max number of blocks assuming 128x128-4x4 all partitions.

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
} Part;

static const uint8_t mi_size_wide[BlockSizeS_ALL] = {1,  1,  2,  2,  2,  4, 4, 4, 8, 8, 8,
                                                     16, 16, 16, 32, 32, 1, 4, 2, 8, 4, 16};
static const uint8_t mi_size_high[BlockSizeS_ALL] = {1, 2,  1,  2,  4,  2, 4, 8, 4, 8,  16,
                                                     8, 16, 32, 16, 32, 4, 1, 8, 2, 16, 4};

typedef char PartitionContextType;
#define PARTITION_PLOFFSET 4 // number of probability models per block size
#define PARTITION_BlockSizeS 5
#define PARTITION_CONTEXTS (PARTITION_BlockSizeS * PARTITION_PLOFFSET)

// block transform size
#ifdef _MSC_VER
typedef uint8_t TxSize;
enum ATTRIBUTE_PACKED {
#else
typedef enum ATTRIBUTE_PACKED {
#endif
    TX_4X4, // 4x4 transform
    TX_8X8, // 8x8 transform
    TX_16X16, // 16x16 transform
    TX_32X32, // 32x32 transform
    TX_64X64, // 64x64 transform
    TX_4X8, // 4x8 transform
    TX_8X4, // 8x4 transform
    TX_8X16, // 8x16 transform
    TX_16X8, // 16x8 transform
    TX_16X32, // 16x32 transform
    TX_32X16, // 32x16 transform
    TX_32X64, // 32x64 transform
    TX_64X32, // 64x32 transform
    TX_4X16, // 4x16 transform
    TX_16X4, // 16x4 transform
    TX_8X32, // 8x32 transform
    TX_32X8, // 32x8 transform
    TX_16X64, // 16x64 transform
    TX_64X16, // 64x16 transform
    TX_SIZES_ALL, // Includes rectangular transforms
    TX_SIZES         = TX_4X8, // Does NOT include rectangular transforms
    TX_SIZES_LARGEST = TX_64X64,
    TX_INVALID       = 255 // Invalid transform size

#ifdef _MSC_VER
};
#else
} TxSize;
#endif
static const TxSize tx_depth_to_tx_size[3][BlockSizeS_ALL] = {
    // tx_depth 0
    {TX_4X4,   TX_4X8,   TX_8X4,   TX_8X8,   TX_8X16,  TX_16X8,  TX_16X16,
     TX_16X32, TX_32X16, TX_32X32, TX_32X64, TX_64X32, TX_64X64,
     TX_64X64, //TX_64X128,
     TX_64X64, //TX_128X64,
     TX_64X64, //TX_128X128,
     TX_4X16,  TX_16X4,  TX_8X32,  TX_32X8,  TX_16X64, TX_64X16},
    // tx_depth 1:
    {TX_4X4,   TX_4X8,   TX_8X4,   TX_4X4,   TX_8X8,   TX_8X8,   TX_8X8,
     TX_16X16, TX_16X16, TX_16X16, TX_32X32, TX_32X32, TX_32X32,
     TX_64X64, //TX_64X128,
     TX_64X64, //TX_128X64,
     TX_64X64, //TX_128X128,
     TX_4X4,   TX_4X4,   TX_8X8,   TX_8X8,   TX_16X16, TX_16X16},
    // tx_depth 2
    {TX_4X4,   TX_4X8, TX_8X4, TX_8X8,   TX_4X4,   TX_4X4,   TX_4X4,
     TX_8X8,   TX_8X8, TX_8X8, TX_16X16, TX_16X16, TX_16X16,
     TX_64X64, //TX_64X128,
     TX_64X64, //TX_128X64,
     TX_64X64, //TX_128X128,
     TX_4X16, // No depth 2
     TX_16X4, // No depth 2
     TX_4X4,   TX_4X4, TX_8X8, TX_8X8}};
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

typedef enum TxClass {
    TX_CLASS_2D    = 0,
    TX_CLASS_HORIZ = 1,
    TX_CLASS_VERT  = 2,
    TX_CLASSES     = 3,
} TxClass;

static INLINE TxSize av1_get_adjusted_tx_size(TxSize tx_size) {
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
#define ALIGN_POWER_OF_TWO(value, n) (((value) + ((1 << (n)) - 1)) & ~((1 << (n)) - 1))
#define AOM_PLANE_Y 0 /**< Y (Luminance) plane */
#define AOM_PLANE_U 1 /**< U (Chroma) plane */
#define AOM_PLANE_V 2 /**< V (Chroma) plane */

#define CONVERT_TO_SHORTPTR(x) ((uint16_t *)(((uintptr_t)(x)) << 1))
#define CONVERT_TO_BYTEPTR(x) ((uint8_t *)(((uintptr_t)(x)) >> 1))

#define AOMMIN(x, y) (((x) < (y)) ? (x) : (y))
#define AOMMAX(x, y) (((x) > (y)) ? (x) : (y))

// frame transform mode
typedef enum ATTRIBUTE_PACKED {
    ONLY_4X4, // use only 4x4 transform
    TX_MODE_LARGEST, // transform size is the largest possible for pu size
    TX_MODE_SELECT, // transform specified for each block
    TX_MODES,
} TxMode;

// 1D tx types
typedef enum ATTRIBUTE_PACKED {
    DCT_1D,
    ADST_1D,
    FLIPADST_1D,
    IDTX_1D,
    TX_TYPES_1D,
} TxType1D;

typedef enum ATTRIBUTE_PACKED {
    DCT_DCT, // DCT  in both horizontal and vertical
    ADST_DCT, // ADST in vertical, DCT in horizontal
    DCT_ADST, // DCT  in vertical, ADST in horizontal
    ADST_ADST, // ADST in both directions
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

typedef enum ATTRIBUTE_PACKED {
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

typedef struct TxfmParam {
    // for both forward and inverse transforms
    TxType  tx_type;
    TxSize  tx_size;
    int32_t lossless;
    int32_t bd;
    // are the pixel buffers octets or shorts?  This should collapse to
    // bd==8 implies !is_hbd, but that's not certain right now.
    int32_t   is_hbd;
    TxSetType tx_set_type;
    // for inverse transforms only
    int32_t eob;
} TxfmParam;

#define IS_2D_TRANSFORM(tx_type) (tx_type < IDTX)
#define EXT_TX_SIZES 4 // number of sizes that use extended transforms
#define EXT_TX_SETS_INTER 4 // Sets of transform selections for INTER
#define EXT_TX_SETS_INTRA 3 // Sets of transform selections for INTRA

typedef enum ATTRIBUTE_PACKED {
    UNIDIR_COMP_REFERENCE,
    BIDIR_COMP_REFERENCE,
    COMP_REFERENCE_TYPES,
} CompReferenceType;

typedef enum ATTRIBUTE_PACKED { PLANE_TYPE_Y, PLANE_TYPE_UV, PLANE_TYPES } PlaneType;

#define CFL_ALPHABET_SIZE_LOG2 4
#define CFL_ALPHABET_SIZE (1 << CFL_ALPHABET_SIZE_LOG2)
#define CFL_MAGS_SIZE ((2 << CFL_ALPHABET_SIZE_LOG2) + 1)
#define CFL_IDX_U(idx) (idx >> CFL_ALPHABET_SIZE_LOG2)
#define CFL_IDX_V(idx) (idx & (CFL_ALPHABET_SIZE - 1))

typedef enum ATTRIBUTE_PACKED { CFL_PRED_U, CFL_PRED_V, CFL_PRED_PLANES } CflPredType;

typedef enum ATTRIBUTE_PACKED { CFL_SIGN_ZERO, CFL_SIGN_NEG, CFL_SIGN_POS, CFL_SIGNS } CflSignType;

typedef enum ATTRIBUTE_PACKED { CFL_DISALLOWED, CFL_ALLOWED, CFL_ALLOWED_TYPES } CflAllowedType;

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
#define CFL_CONTEXT_V(js) (CFL_SIGN_V(js) * CFL_SIGNS + CFL_SIGN_U(js) - CFL_SIGNS)

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

typedef enum ATTRIBUTE_PACKED {
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
typedef enum ATTRIBUTE_PACKED {
    DC_PRED, // Average of above and left pixels
    V_PRED, // Vertical
    H_PRED, // Horizontal
    D45_PRED, // Directional 45  degree
    D135_PRED, // Directional 135 degree
    D113_PRED, // Directional 113 degree
    D157_PRED, // Directional 157 degree
    D203_PRED, // Directional 203 degree
    D67_PRED, // Directional 67  degree
    SMOOTH_PRED, // Combination of horizontal and vertical interpolation
    SMOOTH_V_PRED, // Vertical interpolation
    SMOOTH_H_PRED, // Horizontal interpolation
    PAETH_PRED, // Predict from the direction of smallest gradient
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
    INTRA_MODE_START        = DC_PRED,
    INTRA_MODE_END          = NEARESTMV,
    INTRA_MODE_NUM          = INTRA_MODE_END - INTRA_MODE_START,
    SINGLE_INTER_MODE_START = NEARESTMV,
    SINGLE_INTER_MODE_END   = NEAREST_NEARESTMV,
    SINGLE_INTER_MODE_NUM   = SINGLE_INTER_MODE_END - SINGLE_INTER_MODE_START,
    COMP_INTER_MODE_START   = NEAREST_NEARESTMV,
    COMP_INTER_MODE_END     = MB_MODE_COUNT,
    COMP_INTER_MODE_NUM     = COMP_INTER_MODE_END - COMP_INTER_MODE_START,
    INTRA_MODES             = PAETH_PRED + 1, // PAETH_PRED has to be the last intra mode.
    INTRA_INVALID           = MB_MODE_COUNT, // For uv_mode in inter blocks
    INTRA_MODE_4x4
} PredictionMode;

#define MAX_UPSAMPLE_SZ 16

typedef enum ATTRIBUTE_PACKED {
    UV_DC_PRED, // Average of above and left pixels
    UV_V_PRED, // Vertical
    UV_H_PRED, // Horizontal
    UV_D45_PRED, // Directional 45  degree
    UV_D135_PRED, // Directional 135 degree
    UV_D113_PRED, // Directional 113 degree
    UV_D157_PRED, // Directional 157 degree
    UV_D203_PRED, // Directional 203 degree
    UV_D67_PRED, // Directional 67  degree
    UV_SMOOTH_PRED, // Combination of horizontal and vertical interpolation
    UV_SMOOTH_V_PRED, // Vertical interpolation
    UV_SMOOTH_H_PRED, // Horizontal interpolation
    UV_PAETH_PRED, // Predict from the direction of smallest gradient
    UV_CFL_PRED, // Chroma-from-Luma
    UV_INTRA_MODES,
    UV_MODE_INVALID, // For uv_mode in inter blocks
} UvPredictionMode;

typedef enum ATTRIBUTE_PACKED {
    SIMPLE_TRANSLATION,
    OBMC_CAUSAL, // 2-sided OBMC
    WARPED_CAUSAL, // 2-sided WARPED
    MOTION_MODES
} MotionMode;

typedef enum ATTRIBUTE_PACKED {
    II_DC_PRED,
    II_V_PRED,
    II_H_PRED,
    II_SMOOTH_PRED,
    INTERINTRA_MODES
} InterIntraMode;

typedef enum ATTRIBUTE_PACKED {
    COMPOUND_AVERAGE,
    COMPOUND_DISTWTD,
    COMPOUND_WEDGE,
    COMPOUND_DIFFWTD,
    COMPOUND_TYPES,
    MASKED_COMPOUND_TYPES = 2,
} CompoundType;

#define COMPOUND_INTRA 4 //just for the decoder
#define AOM_BLEND_A64_ROUND_BITS 6
#define AOM_BLEND_A64_MAX_ALPHA (1 << AOM_BLEND_A64_ROUND_BITS) // 64

#define AOM_BLEND_A64(a, v0, v1)                                            \
    ROUND_POWER_OF_TWO((a) * (v0) + (AOM_BLEND_A64_MAX_ALPHA - (a)) * (v1), \
                       AOM_BLEND_A64_ROUND_BITS)
#define DIFF_FACTOR_LOG2 4
#define DIFF_FACTOR (1 << DIFF_FACTOR_LOG2)
#define AOM_BLEND_AVG(v0, v1) ROUND_POWER_OF_TWO((v0) + (v1), 1)
typedef uint16_t CONV_BUF_TYPE;
#define MAX_WEDGE_TYPES (1 << 4)
#define MAX_WEDGE_SIZE_LOG2 5 // 32x32
#define MAX_WEDGE_SIZE (1 << MAX_WEDGE_SIZE_LOG2)
#define MAX_WEDGE_SQUARE (MAX_WEDGE_SIZE * MAX_WEDGE_SIZE)
#define WEDGE_WEIGHT_BITS 6
#define WEDGE_NONE -1
#define MASK_MASTER_SIZE ((MAX_WEDGE_SIZE) << 1)
#define MASK_MASTER_STRIDE (MASK_MASTER_SIZE)
typedef struct {
    int enable_order_hint; // 0 - disable order hint, and related tools
    int order_hint_bits_minus_1; // dist_wtd_comp, ref_frame_mvs,
    int enable_dist_wtd_comp; // 0 - disable dist-wtd compound modes
    int enable_ref_frame_mvs; // 0 - disable ref frame mvs
} OrderHintInfoEnc;
enum {
    MD_COMP_AVG,
    MD_COMP_DIST,
    MD_COMP_DIFF0,
    MD_COMP_WEDGE,
    MD_COMP_TYPES,
} UENUM1BYTE(MD_COMP_TYPE);
#define COMPOUND_TYPE CompoundType
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

#define InterIntraMode InterIntraMode
typedef enum ATTRIBUTE_PACKED {
    FILTER_DC_PRED,
    FILTER_V_PRED,
    FILTER_H_PRED,
    FILTER_D157_PRED,
    FILTER_PAETH_PRED,
    FILTER_INTRA_MODES,
} FilterIntraMode;

static const PredictionMode fimode_to_intramode[FILTER_INTRA_MODES] = {
    DC_PRED, V_PRED, H_PRED, D157_PRED, PAETH_PRED};
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

typedef enum ATTRIBUTE_PACKED {
    LAST_LAST2_FRAMES, // { LAST_FRAME, LAST2_FRAME }
    LAST_LAST3_FRAMES, // { LAST_FRAME, LAST3_FRAME }
    LAST_GOLDEN_FRAMES, // { LAST_FRAME, GOLDEN_FRAME }
    BWDREF_ALTREF_FRAMES, // { BWDREF_FRAME, ALTREF_FRAME }
    LAST2_LAST3_FRAMES, // { LAST2_FRAME, LAST3_FRAME }
    LAST2_GOLDEN_FRAMES, // { LAST2_FRAME, GOLDEN_FRAME }
    LAST3_GOLDEN_FRAMES, // { LAST3_FRAME, GOLDEN_FRAME }
    BWDREF_ALTREF2_FRAMES, // { BWDREF_FRAME, ALTREF2_FRAME }
    ALTREF2_ALTREF_FRAMES, // { ALTREF2_FRAME, ALTREF_FRAME }
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

typedef enum ATTRIBUTE_PACKED {
    RESTORE_NONE,
    RESTORE_WIENER,
    RESTORE_SGRPROJ,
    RESTORE_SWITCHABLE,
    RESTORE_SWITCHABLE_TYPES = RESTORE_SWITCHABLE,
    RESTORE_TYPES            = 4,
} RestorationType;

#define SCALE_NUMERATOR 8
#define SUPERRES_SCALE_BITS 3
#define SUPERRES_SCALE_DENOMINATOR_MIN (SCALE_NUMERATOR + 1)
#define NUM_SCALES 8

//**********************************************************************************************************************//
// onyxc_int.h
#define CDEF_MAX_STRENGTHS 16

#define REF_FRAMES_LOG2 3
#define REF_FRAMES (1 << REF_FRAMES_LOG2)

// 4 scratch frames for the new frames to support a maximum of 4 cores decoding
// in parallel, 3 for scaled references on the encoder.
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
    /*!\brief The given Bitstream is not supported.
    *
    * The Bitstream was unable to be parsed at the highest level. The decoder
    * is unable to proceed. This error \ref SHOULD be treated as fatal to the
    * stream. */
    AOM_CODEC_UNSUP_BITSTREAM,
    /*!\brief Encoded Bitstream uses an unsupported feature
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
    4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32, 64, 64, 64, 128, 128, 4, 16, 8, 32, 16, 64};

static const uint8_t block_size_high[BlockSizeS_ALL] = {
    4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 64, 32, 64, 128, 64, 128, 16, 4, 32, 8, 64, 16};

// AOMMIN(3, AOMMIN(b_width_log2(bsize), b_height_log2(bsize)))
static const uint8_t size_group_lookup[BlockSizeS_ALL] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3,
                                                          3, 3, 3, 3, 3, 0, 0, 1, 1, 2, 2};

static const uint8_t num_pels_log2_lookup[BlockSizeS_ALL] = {
    4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12, 13, 13, 14, 6, 6, 8, 8, 10, 10};
static const TxSize max_txsize_lookup[BlockSizeS_ALL] = {
    //                   4X4
    TX_4X4,
    // 4X8,    8X4,      8X8
    TX_4X4,
    TX_4X4,
    TX_8X8,
    // 8X16,   16X8,     16X16
    TX_8X8,
    TX_8X8,
    TX_16X16,
    // 16X32,  32X16,    32X32
    TX_16X16,
    TX_16X16,
    TX_32X32,
    // 32X64,  64X32,
    TX_32X32,
    TX_32X32,
    // 64X64
    TX_64X64,
    // 64x128, 128x64,   128x128
    TX_64X64,
    TX_64X64,
    TX_64X64,
    // 4x16,   16x4,     8x32
    TX_4X4,
    TX_4X4,
    TX_8X8,
    // 32x8,   16x64     64x16
    TX_8X8,
    TX_16X16,
    TX_16X16};

static const TxSize max_txsize_rect_lookup[BlockSizeS_ALL] = {
    // 4X4
    TX_4X4,
    // 4X8,    8X4,      8X8
    TX_4X8,
    TX_8X4,
    TX_8X8,
    // 8X16,   16X8,     16X16
    TX_8X16,
    TX_16X8,
    TX_16X16,
    // 16X32,  32X16,    32X32
    TX_16X32,
    TX_32X16,
    TX_32X32,
    // 32X64,  64X32,
    TX_32X64,
    TX_64X32,
    // 64X64
    TX_64X64,
    // 64x128, 128x64,   128x128
    TX_64X64,
    TX_64X64,
    TX_64X64,
    // 4x16,   16x4,
    TX_4X16,
    TX_16X4,
    // 8x32,   32x8
    TX_8X32,
    TX_32X8,
    // 16x64,  64x16
    TX_16X64,
    TX_64X16};

// Transform block width in unit
static const int32_t tx_size_wide_unit[TX_SIZES_ALL] = {
    1, 2, 4, 8, 16, 1, 2, 2, 4, 4, 8, 8, 16, 1, 4, 2, 8, 4, 16,
};
// Transform block height in unit
static const int32_t tx_size_high_unit[TX_SIZES_ALL] = {
    1, 2, 4, 8, 16, 2, 1, 4, 2, 8, 4, 16, 8, 4, 1, 8, 2, 16, 4,
};

static const TxSize sub_tx_size_map[TX_SIZES_ALL] = {
    TX_4X4, // TX_4X4
    TX_4X4, // TX_8X8
    TX_8X8, // TX_16X16
    TX_16X16, // TX_32X32
    TX_32X32, // TX_64X64
    TX_4X4, // TX_4X8
    TX_4X4, // TX_8X4
    TX_8X8, // TX_8X16
    TX_8X8, // TX_16X8
    TX_16X16, // TX_16X32
    TX_16X16, // TX_32X16
    TX_32X32, // TX_32X64
    TX_32X32, // TX_64X32
    TX_4X8, // TX_4X16
    TX_8X4, // TX_16X4
    TX_8X16, // TX_8X32
    TX_16X8, // TX_32X8
    TX_16X32, // TX_16X64
    TX_32X16, // TX_64X16
};
static const TxSize txsize_horz_map[TX_SIZES_ALL] = {
    TX_4X4, // TX_4X4
    TX_8X8, // TX_8X8
    TX_16X16, // TX_16X16
    TX_32X32, // TX_32X32
    TX_64X64, // TX_64X64
    TX_4X4, // TX_4X8
    TX_8X8, // TX_8X4
    TX_8X8, // TX_8X16
    TX_16X16, // TX_16X8
    TX_16X16, // TX_16X32
    TX_32X32, // TX_32X16
    TX_32X32, // TX_32X64
    TX_64X64, // TX_64X32
    TX_4X4, // TX_4X16
    TX_16X16, // TX_16X4
    TX_8X8, // TX_8X32
    TX_32X32, // TX_32X8
    TX_16X16, // TX_16X64
    TX_64X64, // TX_64X16
};

static const TxSize txsize_vert_map[TX_SIZES_ALL] = {
    TX_4X4, // TX_4X4
    TX_8X8, // TX_8X8
    TX_16X16, // TX_16X16
    TX_32X32, // TX_32X32
    TX_64X64, // TX_64X64
    TX_8X8, // TX_4X8
    TX_4X4, // TX_8X4
    TX_16X16, // TX_8X16
    TX_8X8, // TX_16X8
    TX_32X32, // TX_16X32
    TX_16X16, // TX_32X16
    TX_64X64, // TX_32X64
    TX_32X32, // TX_64X32
    TX_16X16, // TX_4X16
    TX_4X4, // TX_16X4
    TX_32X32, // TX_8X32
    TX_8X8, // TX_32X8
    TX_64X64, // TX_16X64
    TX_16X16, // TX_64X16
};
static const uint8_t mi_size_wide_log2[BlockSizeS_ALL] = {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
                                                          4, 4, 4, 5, 5, 0, 2, 1, 3, 2, 4};
static const uint8_t mi_size_high_log2[BlockSizeS_ALL] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4,
                                                          3, 4, 5, 4, 5, 2, 0, 3, 1, 4, 2};

typedef struct SgrParamsType {
    int32_t r[2]; // radii
    int32_t s[2]; // sgr parameters for r[0] and r[1], based on GenSgrprojVtable()
} SgrParamsType;

//**********************************************************************************************************************//
// blockd.h
typedef enum FrameType {
    KEY_FRAME        = 0,
    INTER_FRAME      = 1,
    INTRA_ONLY_FRAME = 2, // replaces intra-only
    S_FRAME          = 3,
    FRAME_TYPES,
} FrameType;

typedef int8_t MvReferenceFrame;

// Number of transform types in each set type

static const int32_t av1_num_ext_tx_set[EXT_TX_SET_TYPES] = {
    1,
    2,
    5,
    7,
    12,
    16,
};

static const int32_t av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
};

static INLINE TxSetType get_ext_tx_set_type(TxSize tx_size, int32_t is_inter,
                                            int32_t use_reduced_set) {
    const TxSize tx_size_sqr_up = txsize_sqr_up_map[tx_size];

    if (tx_size_sqr_up > TX_32X32)
        return EXT_TX_SET_DCTONLY;
    if (tx_size_sqr_up == TX_32X32)
        return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DCTONLY;
    if (use_reduced_set)
        return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DTT4_IDTX;
    const TxSize tx_size_sqr = txsize_sqr_map[tx_size];
    if (is_inter) {
        return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT9_IDTX_1DDCT : EXT_TX_SET_ALL16);
    } else {
        return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT4_IDTX : EXT_TX_SET_DTT4_IDTX_1DDCT);
    }
}
static INLINE int32_t get_ext_tx_types(TxSize tx_size, int32_t is_inter, int32_t use_reduced_set) {
    const int32_t set_type = get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
    return av1_num_ext_tx_set[set_type];
}
// Maps tx set types to the indices.
static const int32_t ext_tx_set_index[2][EXT_TX_SET_TYPES] = {
    {// Intra
     0,
     -1,
     2,
     1,
     -1,
     -1},
    {// Inter
     0,
     3,
     -1,
     -1,
     2,
     1},
};

static INLINE int32_t get_ext_tx_set(TxSize tx_size, int32_t is_inter, int32_t use_reduced_set) {
    const TxSetType set_type = get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
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
typedef enum FrameContextIndex {
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
#define INVALID_IDX -1 // Invalid buffer index.

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
#define NO_FILTER_FOR_IBC 1 // Disable in-loop filters for frame with intrabc
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
    int8_t  mode_deltas[MAX_MODE_LF_DELTAS];
    int32_t combine_vert_horz_lf;
};

#define MAX_SEGMENTS 8
#define MAX_MB_PLANE 3

#define MAX_LOOP_FILTER 63
#define MAX_SHARPNESS 7

#define SIMD_WIDTH 16
// Need to align this structure so when it is declared and
// passed it can be loaded into vector registers.
typedef struct LoopFilterThresh {
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, mblim[SIMD_WIDTH]);
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, lim[SIMD_WIDTH]);
    DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, hev_thr[SIMD_WIDTH]);
} LoopFilterThresh;

typedef struct LoopFilterInfoN {
    LoopFilterThresh lfthr[MAX_LOOP_FILTER + 1];
    uint8_t          lvl[MAX_MB_PLANE][MAX_SEGMENTS][2][REF_FRAMES][MAX_MODE_LF_DELTAS];
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
#define GM_ROW3HOMO_PREC_DIFF (WARPEDMODEL_ROW3HOMO_PREC_BITS - GM_ROW3HOMO_PREC_BITS)
#define GM_ROW3HOMO_DECODE_FACTOR (1 << GM_ROW3HOMO_PREC_DIFF)

#define GM_TRANS_MAX (1 << GM_ABS_TRANS_BITS)
#define GM_ALPHA_MAX (1 << GM_ABS_ALPHA_BITS)
#define GM_ROW3HOMO_MAX (1 << GM_ABS_ROW3HOMO_BITS)

#define GM_TRANS_MIN -GM_TRANS_MAX
#define GM_ALPHA_MIN -GM_ALPHA_MAX
#define GM_ROW3HOMO_MIN -GM_ROW3HOMO_MAX
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

#ifdef _WIN32
#define NOINLINE                __declspec ( noinline )
#define FORCE_INLINE            __forceinline
#else
#define NOINLINE                __attribute__(( noinline ))
#define FORCE_INLINE            __attribute__((always_inline))
#endif

// ***************************** Definitions *****************************

#define MAX_BITS_PER_FRAME            8000000

#define SC_FRAMES_TO_IGNORE     1000 // The speed control algorith starts after SC_FRAMES_TO_IGNORE number frames.
#define SC_FRAMES_INTERVAL_SPEED      60 // The speed control Interval To Check the speed
#define SC_FRAMES_INTERVAL_T1         60 // The speed control Interval Threshold1
#define SC_FRAMES_INTERVAL_T2        180 // The speed control Interval Threshold2
#define SC_FRAMES_INTERVAL_T3        120 // The speed control Interval Threshold3

#define LAST_BWD_FRAME     8
#define LAST_ALT_FRAME    16

#define MAX_NUM_TOKENS          200


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
/** The MD_BIT_DEPTH_MODE type is used to describe the bitdepth of MD path.
*/

typedef enum MD_BIT_DEPTH_MODE
{
    EB_8_BIT_MD     = 0,    // 8bit mode decision
    EB_10_BIT_MD    = 1,    // 10bit mode decision
    EB_DUAL_BIT_MD  = 2     // Auto: 8bit & 10bit mode decision
} MD_BIT_DEPTH_MODE;

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
            SVT_LOG("Malloc has failed due to insuffucient resources"); \
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
SVT_LOG("Total Number of Mallocs in Library: %d\n", lib_malloc_count); \
SVT_LOG("Total Number of Threads in Library: %d\n", lib_thread_count); \
SVT_LOG("Total Number of Semaphore in Library: %d\n", lib_semaphore_count); \
SVT_LOG("Total Number of Mutex in Library: %d\n", lib_mutex_count); \
SVT_LOG("Total Library Memory: %.2lf KB\n\n",*total_lib_memory/(double)1024);

#define EB_APP_MEMORY() \
SVT_LOG("Total Number of Mallocs in App: %d\n", app_malloc_count); \
SVT_LOG("Total App Memory: %.2lf KB\n\n",*total_app_memory/(double)1024);

#ifndef _ERRNO_T_DEFINED
#define _ERRNO_T_DEFINED
typedef int32_t errno_t;
#endif  /* _ERRNO_T_DEFINED */

extern void
    eb_memcpy_app(void  *dst_ptr, void  *src_ptr, size_t size);
#ifdef ARCH_X86
#define EB_MEMCPY(dst, src, size) \
    eb_memcpy_app(dst, src, size)
#else
#define EB_MEMCPY(dst, src, size) \
    memcpy(dst, src, size)
#endif

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
EbPtr app_private_data;
EbPtr handle;
void(*error_handler)(
    EbPtr handle,
    uint32_t errorCode);
} EbCallback;

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
#define LOG2F_MAX_SB_SIZE                          6u
#define LOG2_64_SIZE                                6 // log2(BLOCK_SIZE_64)
#define MAX_LEVEL_COUNT                             5 // log2(BLOCK_SIZE_64) - log2(MIN_BLOCK_SIZE)
#define LOG_MIN_BLOCK_SIZE                          3
#define MIN_BLOCK_SIZE                              (1 << LOG_MIN_BLOCK_SIZE)
#define LOG_MIN_PU_SIZE                             2
#define MIN_PU_SIZE                                 (1 << LOG_MIN_PU_SIZE)
#define MAX_NUM_OF_PU_PER_CU                        1
#define MAX_NUM_OF_REF_PIC_LIST                     2
#define MAX_NUM_OF_PART_SIZE                        8
#define EB_MAX_SB_DEPTH                            (((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 1) ? 1 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 2) ? 2 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 4) ? 3 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 8) ? 4 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 16) ? 5 : \
                                                    ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) == 32) ? 6 : 7)
#define MIN_CU_BLK_COUNT                            ((BLOCK_SIZE_64 / MIN_BLOCK_SIZE) * (BLOCK_SIZE_64 / MIN_BLOCK_SIZE))
#define MAX_NUM_OF_TU_PER_CU                        21
#define MIN_NUM_OF_TU_PER_CU                        5
#define MAX_SB_ROWS                                ((MAX_PICTURE_HEIGHT_SIZE) / (BLOCK_SIZE_64))

#define MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE       ((MAX_PICTURE_WIDTH_SIZE + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64) * \
                                                ((MAX_PICTURE_HEIGHT_SIZE + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64)
// super-resolution definitions
#define MIN_SUPERRES_DENOM                          8
#define MAX_SUPERRES_DENOM                          16

//***Prediction Structure***
#define MAX_TEMPORAL_LAYERS                         6
#define MAX_REF_IDX                                 4
#define INVALID_POC                                 (((uint32_t) (~0)) - (((uint32_t) (~0)) >> 1))
#define MAX_ELAPSED_IDR_COUNT                       1024

typedef enum DownSamplingMethod
{
    ME_FILTERED_DOWNSAMPLED  = 0,
    ME_DECIMATED_DOWNSAMPLED = 1
} DownSamplingMethod;

//***Segments***
#define CU_MAX_COUNT                                85

#define MAX_INTRA_REFERENCE_SAMPLES                 (BLOCK_SIZE_64 << 2) + 1

#define _MVXT(mv) ( (int16_t)((mv) &  0xFFFF) )
#define _MVYT(mv) ( (int16_t)((mv) >> 16    ) )

//***MCP***

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

#define TRANS_BIT_INCREMENT    0
#define QUANT_IQUANT_SHIFT     20 // Q(QP%6) * IQ(QP%6) = 2^20
#define QUANT_SHIFT            14 // Q(4) = 2^14
#define SCALE_BITS             15 // Inherited from TMuC, pressumably for fractional bit estimates in RDOQ
#define MAX_TR_DYNAMIC_RANGE   15 // Maximum transform dynamic range (excluding sign bit)
#define MAX_POS_16BIT_NUM      32767
#define MIN_NEG_16BIT_NUM      -32768
#define QUANT_OFFSET_I         171
#define QUANT_OFFSET_P         85
#define LOW_SB_VARIANCE        10
#define MEDIUM_SB_VARIANCE        50

/*********************************************************
* used for the first time, but not the last time interpolation filter
*********************************************************/
#define Shift1       InternalBitDepthIncrement
#define MinusOffset1 (1 << (IF_Negative_Offset + InternalBitDepthIncrement))

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

#define QP_OFFSET_SB_SCORE_0    0
#define QP_OFFSET_SB_SCORE_1    50
#define QP_OFFSET_SB_SCORE_2    100
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
#define GRADUAL_LUMINOSITY_CHANGE_TH                        3
#define FADED_SB_PERCENTAGE_TH                             10
#define FADED_PICTURES_TH                                   15
#define CLASS_SUB_0_PICTURE_ACTIVITY_REGIONS_TH             1
#define CLASS_1_SIZE_PICTURE_ACTIVITY_REGIONS_TH            2
#define HIGHER_THAN_CLASS_1_PICTURE_ACTIVITY_REGIONS_TH     8

#define IS_COMPLEX_SB_VARIANCE_TH                          100
#define IS_COMPLEX_SB_FLAT_VARIANCE_TH                     10
#define IS_COMPLEX_SB_VARIANCE_DEVIATION_TH                13
#define IS_COMPLEX_SB_ZZ_SAD_FACTOR_TH                     25

#define MAX_SUPPORTED_SEGMENTS                            7
#define NUM_QPS                                           52


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

// Multi-Pass Partitioning Depth(Multi - Pass PD) performs multiple PD stages for the same SB towards 1 final Partitioning Structure
// As we go from PDn to PDn + 1, the prediction accuracy of the MD feature(s) increases while the number of block(s) decreases
#if DEPTH_PART_CLEAN_UP
#if ADD_NEW_MPPD_LEVEL
typedef enum MultiPassPdLevel
{
    MULTI_PASS_PD_OFF     = 0, // Multi-Pass PD OFF = 1-single PD Pass (e.g. I_SLICE, SC)
    MULTI_PASS_PD_LEVEL_0 = 1, // Multi-Pass PD Mode 0: PD0 | PD0_REFINEMENT
    MULTI_PASS_PD_LEVEL_1 = 2, // Multi-Pass PD Mode 1: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT
    MULTI_PASS_PD_LEVEL_2 = 3, // Multi-Pass PD Mode 1: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ only
    MULTI_PASS_PD_LEVEL_3 = 4, // Multi-Pass PD Mode 2: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ and SQ coeff info
    MULTI_PASS_PD_LEVEL_4 = 5, // reserved = MULTI_PASS_PD_LEVEL_3
    MULTI_PASS_PD_INVALID = 6, // Invalid Multi-Pass PD Mode
} MultiPassPdLevel;
#else
typedef enum MultiPassPdLevel
{
    MULTI_PASS_PD_OFF = 0, // Multi-Pass PD OFF = 1-single PD Pass (e.g. I_SLICE, SC)
    MULTI_PASS_PD_LEVEL_0 = 1, // Multi-Pass PD Mode 0: PD0 | PD0_REFINEMENT
    MULTI_PASS_PD_LEVEL_1 = 2, // Multi-Pass PD Mode 1: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ only
    MULTI_PASS_PD_LEVEL_2 = 3, // Multi-Pass PD Mode 2: PD0 | PD0_REFINEMENT | PD1 | PD1_REFINEMENT using SQ vs. NSQ and SQ coeff info
    MULTI_PASS_PD_LEVEL_3 = 4, // reserved = MULTI_PASS_PD_LEVEL_2
    MULTI_PASS_PD_INVALID = 5, // Invalid Multi-Pass PD Mode
} MultiPassPdLevel;
#endif

typedef enum AdpLevel
{
    ADP_OFF = 0, // All SBs use the same Multi-Pass PD level
    ADP_LEVEL_1 = 1, // read @ ADP budget derivation (e.g. high budget_boost)
    ADP_LEVEL_2 = 2, // read @ ADP budget derivation (e.g. moderate budget_boost)
    ADP_LEVEL_3 = 3, // read @ ADP budget derivation (e.g. low budget_boost)
} AdpLevel;
#else
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
#if SHUT_ME_DISTORTION
    PIC_SB_SWITCH_DEPTH_MODE    = 8  // Adaptive Depth Partitioning
#else
    PIC_OPEN_LOOP_DEPTH_MODE    = 8, // Early Inter Depth Decision:  SB size -> 8x8
    PIC_SB_SWITCH_DEPTH_MODE    = 9  // Adaptive Depth Partitioning
#endif

} EbPictureDepthMode;
#endif
#define EB_SB_DEPTH_MODE              uint8_t
#define SB_SQ_BLOCKS_DEPTH_MODE             1
#define SB_SQ_NON4_BLOCKS_DEPTH_MODE        2
#if !SHUT_ME_DISTORTION
#define SB_OPEN_LOOP_DEPTH_MODE             3
#define SB_FAST_OPEN_LOOP_DEPTH_MODE        4
#define SB_PRED_OPEN_LOOP_DEPTH_MODE        5
#endif
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

static const uint32_t raster_scan_blk_x[CU_MAX_COUNT] =
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

static const uint32_t raster_scan_blk_y[CU_MAX_COUNT] =
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

static const uint32_t raster_scan_blk_size[CU_MAX_COUNT] =
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

static const uint32_t raster_scan_blk_parent_index[CU_MAX_COUNT] =
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

static const uint32_t md_scan_to_ois_32x32_scan[CU_MAX_COUNT] =
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

typedef struct StatStruct
{
    uint32_t                        referenced_area[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
} StatStruct;
#define TWO_PASS_IR_THRSHLD 40  // Intra refresh threshold used to reduce the reference area.
                                // If the periodic Intra refresh is less than the threshold,
                                // the referenced area is normalized
#define SC_MAX_LEVEL 2 // 2 sets of HME/ME settings are used depending on the scene content mode

#if !REFACTOR_ME_HME
/******************************************************************************
                            ME/HME settings
*******************************************************************************/
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 }       // INPUT_SIZE_4K_RANGE
    }
};
//     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t enable_hme_level0_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 }       // INPUT_SIZE_4K_RANGE
    }
};

static const uint16_t hme_level0_total_search_area_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 }
#else
#if MAR25_ADOPTIONS
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 }
#else
#if MAR23_ADOPTIONS
        { 200,  200,   48,   48,   48,    48,    48,   48,   48,   48,   48,   48,   48 },
        { 256,  256,   96,   96,   96,    96,    96,   96,   96,   48,   48,   48,   48 },
        { 320,  320,  128,  128,  128,   128,   128,  128,  128,   48,   48,   48,   48 },
        { 400,  400,  128,  128,  128,   128,   128,  128,  128,   96,   96,   96,   96 }
     } , {
        { 400,  400,  400,  400,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        { 512,  512,  512,  512,   96,   96,   96,   96,   96,   96,   96,   96,   96 },
        { 640,  640,  640,  640,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 640,  640,  640,  640,  128,  128,  128,  128,  128,  128,  128,  128,  128 }
#else
#if MAR17_ADOPTIONS
        {  64,   64,   48,   48,   48,    48,    48,   48,   48,   48,   48,   48,   48 },
        {  96,   96,   96,   96,   96,    96,    96,   96,   96,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   128,   128,  128,  128,   48,   48,   48,   48 },
        { 192,  192,  128,  128,  128,   128,   128,  128,  128,   96,   96,   96,   96 }
     } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  96,   96,   96,   96,   96,   96,   96,   96,   96,   96,   96,   96,   96 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 192,  192,  192,  192,  128,  128,  128,  128,  128,  128,  128,  128,  128 }
#else
#if MAR12_ADOPTIONS
        {  64,   64,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        {  96,   96,   96,   96,   96,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 192,  192,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
     } , {
        {  64,   64,   64,   64,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        {  96,   96,   96,   96,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 192,  192,  192,  192,  128,  128,  128,  128,  128,  128,  128,  128,  128 }
#else
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        {  96,   96,   96,   96,   96,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
     } , {
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 },
        { 128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128,  128 }
#endif
#endif
#endif
#endif
#endif
    }
};

static const uint16_t hme_level0_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        { 100,   100,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 128,   128,   48,   48,   48,   48,   48,   48,   48,   24,   24,   24,   24 },
        { 160,   160,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 200,   200,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        { 200,  200,   200,  200,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        { 256,  256,   256,  256,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR17_ADOPTIONS
        { 32,   32,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 48,   48,   56,   56,   56,   56,   56,   56,   56,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 96,   96,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  96,   96,   96,   96,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR12_ADOPTIONS
        { 32,   32,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 48,   48,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 96,   96,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }

    } , {
        {  32,   32,   32,   32,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  48,   48,   48,   48,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  96,   96,   96,   96,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
        { 24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }

    } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t hme_level0_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        { 100,   100,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 128,   128,   48,   48,   48,   48,   48,   48,   48,   24,   24,   24,   24 },
        { 160,   160,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 200,   200,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        { 200,  200,   200,  200,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        { 256,  256,   256,  256,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR17_ADOPTIONS
        { 32,   32,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 48,   48,   56,   56,   56,   56,   56,   56,   56,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 96,   96,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  96,   96,   96,   96,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR12_ADOPTIONS
        { 32,   32,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 48,   48,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 96,   96,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }

    } , {
        {  32,   32,   32,   32,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  48,   48,   48,   48,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  96,   96,   96,   96,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
        { 24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }

    } , {
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {  64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t hme_level0_total_search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 }
#else
#if MAR25_ADOPTIONS
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 }
#else
#if MAR23_ADOPTIONS
        {  200,   200,   48,   48,   48,   48,   48,   48,    48,   48,  48,  48,   48 },
        {  256,   256,   96,   96,   96,   96,   96,   96,    96,   48,  48,  48,   48 },
        {  320,   320,  128,  128,  128,  128,  128,  128,   128,   48,  48,  48,   48 },
        {  400,   400,  128,  128,  128,  128,  128,  128,   128,   96,  96,  96,   96 }
    } , {
        {  400,  400,  400,   400,   64,   64,   64,   64,   64,    64,   64,   64,   64 },
        {  512,  512,  512,   512,   96,   96,   96,   96,   96,    96,   96,   96,   96 },
        {  640,  640,  640,   640,  128,  128,  128,  128,  128,   128,  128,  128,  128 },
        {  640,  640,  640,   640,  128,  128,  128,  128,  128,   128,  128,  128,  128 }
#else
#if MAR17_ADOPTIONS
        {   32,    32,   40,   40,   40,   40,   40,   40,   40,   24,   24,   24,   24 },
        {   64,    64,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        {   80,    80,   80,   80,   80,   80,   80,   80,   80,   24,   24,   24,   24 },
        {  128,   128,   80,   80,   80,   80,   80,   80,   80,   24,   24,   24,   24 }
    } , {
        {   32,    32,   32,    32,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        {   64,    64,   64,    64,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        {   80,    80,   80,    80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  128,   128,  128,   128,   80,   80,   80,   80,   80,   80,   80,   80,   80 }
#else
#if MAR12_ADOPTIONS
        {   32,    32,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {   64,    64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {   80,    80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  128,   128,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    } , {
        {   32,    32,   32,    32,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {   64,    64,   64,    64,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {   80,    80,   80,    80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  128,   128,  128,   128,   80,   80,   80,   80,   80,   80,   80,   80,   80 }
#else
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    } , {
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 },
        {  80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80,   80 }
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t hme_level0_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        { 100,   100,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 128,   128,   48,   48,   48,   48,   48,   48,   48,   24,   24,   24,   24 },
        { 160,   160,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 200,   200,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        { 200,  200,   200,  200,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        { 256,  256,   256,  256,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR17_ADOPTIONS
        {  16,   16,   20,   20,   20,   20,   20,   20,   20,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   12,   12,   12,   12 },
        {  64,   64,   40,   40,   40,   40,   40,   40,   40,   12,   12,   12,   12 }
    } , {
        {  16,   16,  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  32,   32,  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        {  40,   40,  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  64,   64,  64,   64,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#else
#if MAR12_ADOPTIONS
        {  16,   16,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  64,   64,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    } , {
        {  16,   16,  16,   16,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  32,   32,  32,   32,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  64,   64,  64,   64,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#else
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    } , {
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t hme_level0_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        { 100,   100,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        { 128,   128,   48,   48,   48,   48,   48,   48,   48,   24,   24,   24,   24 },
        { 160,   160,   64,   64,   64,   64,   64,   64,   64,   24,   24,   24,   24 },
        { 200,   200,   64,   64,   64,   64,   64,   64,   64,   48,   48,   48,   48 }

    } , {
        { 200,  200,   200,  200,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        { 256,  256,   256,  256,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 },
        { 320,  320,   320,  320,   64,   64,   64,   64,   64,   64,   64,   64,   64 }
#else
#if MAR17_ADOPTIONS
        {  16,   16,   20,   20,   20,   20,   20,   20,   20,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   12,   12,   12,   12 },
        {  64,   64,   40,   40,   40,   40,   40,   40,   40,   12,   12,   12,   12 }
    }, {
        {  16,   16,  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  32,   32,  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32 },
        {  40,   40,  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  64,   64,  64,   64,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#else
#if MAR12_ADOPTIONS
        {  16,   16,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  64,   64,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }, {
        {  16,   16,  16,   16,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  32,   32,  32,   32,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  64,   64,  64,   64,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#else
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }, {
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40 }
#endif
#endif
#endif
#endif
#endif
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
#if MAR12_ADOPTIONS
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#else
        {   1,    1,    1,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#endif
    }
};
static const uint16_t hme_level1_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR17_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }

#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
#if MAR12_ADOPTIONS
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
#endif
#endif
    }
};
static const uint16_t hme_level1_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR17_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }

#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
#if MAR12_ADOPTIONS
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
#endif
#endif
    }
};
static const uint16_t hme_level1_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR17_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }

#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
#if MAR12_ADOPTIONS
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
#endif
#endif
    }
};
static const uint16_t hme_level1_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR17_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }

#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
#if MAR12_ADOPTIONS
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   8,    8,    8,    8,     8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,   32,    8,    8,    8,     8 }
#endif
#endif
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
#if MAR12_ADOPTIONS
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#else
        {   1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#endif
    }
};
static const uint16_t hme_level2_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
#if MAR17_ADOPTIONS
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
#endif
    }
};
static const uint16_t hme_level2_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
#if MAR17_ADOPTIONS
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
#endif
    }
};
static const uint16_t hme_level2_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
#if MAR17_ADOPTIONS
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
#endif
    }
};
static const uint16_t hme_level2_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
#if MAR17_ADOPTIONS
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    8,    8,    8,    8,    4,    4,     4,    4 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
#endif
    }
};
/******************************************************************************
                          MAX & MIN  ME search region
*******************************************************************************/
    //    M0      M1      M2      M3      M4    M5       M6      M7      M8      M9      M10      M11    M12
static const uint16_t max_me_search_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 }
#else
#if MAR25_ADOPTIONS
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 200,    200,    114,    114,    114,    114,    114,    114,    114,    32,    32,    32,    32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 256,    256,    170,    170,    136,    136,    136,    136,    136,    48,    48,    48,    48 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 320,    320,    212,    212,    186,    186,    186,    186,    186,    64,    64,    64,    64 },  // INPUT_SIZE_1080p_RANGE
        { 400,    400,    256,    256,    186,    186,    186,    186,    186,    96,    96,    96,    96 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 400,   400,   400,    400,   368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 512,   512,   512,    512,   368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 640,   640,   640,    640,   640,    640,    640,    640,    504,    360,    360,   360,  360 },
        { 640,   640,   640,    640,   640,    640,    640,    640,    540,    540,    540,   540,  540 }
#else
#if ADOPT_SQ_ME_SEARCH_AREA
        // NSC
        { 226,    226,    114,    114,    114,    114,    114,    114,    114,    32,    32 ,   32,    32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 340,    340,    170,    170,    136,    136,    136,    136,    136,    48,    48,    48,    48 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 424,    424,    212,    212,    186,    186,    186,    186,    186,    64,    64,    64,    64 },  // INPUT_SIZE_1080p_RANGE
        { 510,    510,    256,    256,    186,    186,    186,    186,    186,    96,    96,    96,    96 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 600,    600,    600,     600,    368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 900,    900,    900,     900,    368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 1200,   1200,   1200,    1200,   784,    784,    784,    784,    504,    360,    360,   360,  360 },
        { 1800,   1800,   1800,    1800,   784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR19_ADOPTIONS
        // NSC
        { 320,    320,    256,    256,    160,    160,    160,    160,    160,    32,    32 ,   32,    32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 480,    480,    384,    384,    192,    192,    192,    192,    192,    48,    48,    48,    48 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 600,    600,    480,    480,    240,    240,    240,    240,    240,    64,    64,    64,    64 },  // INPUT_SIZE_1080p_RANGE
        { 720,    720,    720,    720,    240,    240,    240,    240,    240,    96,    96,    96,    96 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 600,    600,    600,     600,    368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 900,    900,    900,     900,    368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 1200,   1200,   1200,    1200,   784,    784,    784,    784,    504,    360,    360,   360,  360 },
        { 1800,   1800,   1800,    1800,   784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 256,    256,    256,    256,     96,     96,     96,     96,    32,    32,    32 ,   32,    32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 384,    384,    384,    384,    200,    200,    200,    200,    48,    48,    48,    48,    48 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 480,    480,    480,    480,    296,    296,    296,    296,    64,    64,    64,    64,    64 },  // INPUT_SIZE_1080p_RANGE
        { 720,    720,    720,    720,    296,    296,    296,    296,    96,    96,    96,    96,    96 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 512,    512,    512,    512,    368,    368,    368,    368,    184,    184,    184,   184,  184 },
        { 800,    800,    800,    800,    368,    368,    368,    368,    280,    280,    280,   280,  280 },
        { 1024,   1024,   1024,   1024,   784,    784,    784,    784,    400,    360,    360,   360,  360 },
        { 1536,   1536,   1536,   1536,   784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 256,    256,    256,    256,    96 ,    96,     72 ,    72,     32,     72,     72 ,    72,     72 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 384,    384,    384,    384,    200,    200,    152,    152,    48,    152,    152,    152,    152 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 480,    480,    480,    480,    296,    296,    224,    224,    64,    224,    224,    224,    224 },  // INPUT_SIZE_1080p_RANGE
        { 720,    720,    720,    720,    296,    296,    224,    224,    96,    224,    224,    224,    224 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 512,    512,    512,    512,    368,    368,    280,    280,    128,    280,    280,    280,    280 },
        { 800,    800,    800,    800,    368,    368,    280,    280,    200,    280,    280,    280,    280 },
        { 1024,   1024,   1024,   1024,   784,    784,    600,    600,    256,    600,    600,    600,    600 },
        { 1536,   1536,   1536,   1536,   784,    784,    600,    600,    384,    600,    600,    600,    600 }
#else
#if MAR10_ADOPTIONS
        // NSC
        { 256,    256,    256,    128,    96 ,    96,     72 ,    72,     32,     72,     72 ,    72,     72 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 384,    384,    384,    256,    200,    200,    152,    152,    48,    152,    152,    152,    152 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 480,    480,    480,    384,    296,    296,    224,    224,    64,    224,    224,    224,    224 },  // INPUT_SIZE_1080p_RANGE
        { 720,    720,    720,    384,    296,    296,    224,    224,    96,    224,    224,    224,    224 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 512,    512,    512,    480,    368,    368,    280,    280,    128,    280,    280,    280,    280 },
        { 800,    800,    800,    480,    368,    368,    280,    280,    200,    280,    280,    280,    280 },
        { 1024,   1024,   1024,   1024,   784,    784,    600,    600,    256,    600,    600,    600,    600 },
        { 1536,   1536,   1536,   1024,   784,    784,    600,    600,    384,    600,    600,    600,    600 }
#else
        // NSC
        { 256,   256,    128,    128,    96 ,    96,     72 ,    72,     72,     72,     72 ,    72,     72  }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 384,   384,    256,    256,    200,    200,    152,    152,    152,    152,    152,    152,    152 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 480,   480,    384,    384,    296,    296,    224,    224,    224,    224,    224,    224,    224 }, // INPUT_SIZE_1080p_RANGE
        { 720,   720,    384,    384,    296,    296,    224,    224,    224,    224,    224,    224,    224 }  // INPUT_SIZE_4K_RANGE
    },{
    // SC
        {  512,    512,    480,    480,    368,    368,    280,    280,    280,    280,    280,    280,    280 },
        {  800,    800,    480,    480,    368,    368,    280,    280,    280,    280,    280,    280,    280 },
        {  1024,   800 ,   1024,   1024,   784,    784,    600,    600,    600,    600,    600,    600,    600 },
        {  1536,   1536,   1024,   1024,   784,    784,    600,    600,    600,    600,    600,    600,    600 }
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t max_me_search_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  150,   150,  150,  150,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  450,   450,  450,  450,  450,  450,  450,  450,  450 }
#else
#if MAR25_ADOPTIONS
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 },
        { 300,  300,  300,  300,  300,   300,   300,  300,  150,  150,  150,  150,  150 }
     } , {
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 },
        { 780,  780,  780,  780,  780,   780,   780,  780,  450,  450,  450,  450,  450 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 200,    200,    114,    114,    114,    114,    114,    114,    114,    32,    32,    32,    32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 256,    256,    170,    170,    136,    136,    136,    136,    136,    48,    48,    48,    48 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 320,    320,    212,    212,    186,    186,    186,    186,    186,    64,    64,    64,    64 },  // INPUT_SIZE_1080p_RANGE
        { 400,    400,    256,    256,    186,    186,    186,    186,    186,    96,    96,    96,    96 }   // INPUT_SIZE_4K_RANGE
    },{
        // SC
        { 400,   400,   400,    400,   368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 512,   512,   512,    512,   368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 640,   640,   640,    640,   640,    640,    640,    640,    504,    360,    360,   360,  360 },
        { 640,   640,   640,    640,   640,    640,    640,    640,    540,    540,    540,   540,  540 }
#else
#if ADOPT_SQ_ME_SEARCH_AREA
        // NSC
        { 226,    226,    114,    114,   114,    114,    114,    114,    114,    16,   16,   16,   16 },
        { 340,    340,    170,    170,   136,    136,    136,    136,    136,    24,   24,   24,   24 },
        { 424,    424,    212,    212,   186,    186,    186,    186,    170,    32,   32,   32,   32 },
        { 510,    510,    256,    256,   186,    186,    186,    186,    170,    48,   48,   48,   48 }
    },{
        // SC
        { 600,   600,     600,     600,     368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 900,   900,     900,     900,     368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 1200,  1200,    1200,    1200,    784,    784,    784,    784,    504,    360,    360,   360,  360 },
        { 1800,  1800,    1800,    1800,    784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR19_ADOPTIONS
        // NSC
        { 160,    160,    128,    128,     80,     80,     80,     80,    80,    16,   16,   16,   16 },
        { 240,    240,    192,    192,     96,     96,     96,     96,    96,    24,   24,   24,   24 },
        { 300,    300,    240,    240,    144,    144,    144,    144,    120,    32,   32,   32,   32 },
        { 360,    360,    360,    360,    144,    144,    144,    144,    120,    48,   48,   48,   48 }
    },{
        // SC
        { 600,   600,     600,     600,     368,    368,    368,    368,    336,    184,    184,   184,  184 },
        { 900,   900,     900,     900,     368,    368,    368,    368,    408,    280,    280,   280,  280 },
        { 1200,  1200,    1200,    1200,    784,    784,    784,    784,    504,    360,    360,   360,  360 },
        { 1800,  1800,    1800,    1800,    784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 128,    128,    128,    128,     48,     48,     48,     48,    16,    16,   16,   16,   16 },
        { 192,    192,    192,    192,     96,     96,     96,     96,    24,    24,   24,   24,   24 },
        { 240,    240,    240,    240,    144,    144,    144,    144,    32,    32,   32,   32,   32 },
        { 360,    360,    360,    360,    144,    144,    144,    144,    48,    48,   48,   48,   48 }
    },{
        // SC
        { 512,    512,     512,     512,     368,    368,    368,    368,    184,    184,    184,   184,  184 },
        { 800,    800,     800,     800,     368,    368,    368,    368,    280,    280,    280,   280,  280 },
        { 1024,   1024,    1024,    1024,    784,    784,    784,    784,    400,    360,    360,   360,  360 },
        { 1536,   1536,    1536,    1536,    784,    784,    784,    784,    540,    540,    540,   540,  540 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 128,    128,    128,    128,    48,    48,      40,     40,     16,     40,     40,     40,     40 },
        { 192,    192,    192,    192,    96,    96,      72,     72,     24,     72,     72,     72,     72 },
        { 240,    240,    240,    240,    144,    144,    112,    112,    32,    112,    112,    112,    112 },
        { 360,    360,    360,    360,    144,    144,    112,    112,    48,    112,    112,    112,    112 }
    },{
        // SC
        { 512,    512,     512,     512,     368,    368,    280,    280,    128,    280,    280,    280,    280 },
        { 800,    800,     800,     800,     368,    368,    280,    280,    200,    280,    280,    280,    280 },
        { 1024,   1024,    1024,    1024,    784,    784,    600,    600,    256,    600,    600,    600,    600 },
        { 1536,   1536,    1536,    1536,    784,    784,    600,    600,    384,    600,    600,    600,    600 }
#else
#if MAR10_ADOPTIONS
        // NSC
        { 128,    128,    128,    64,     48,    48,      40,     40,     16,     40,     40,     40,     40 },
        { 192,    192,    192,    128,    96,    96,      72,     72,     24,     72,     72,     72,     72 },
        { 240,    240,    240,    192,    144,    144,    112,    112,    32,    112,    112,    112,    112 },
        { 360,    360,    360,    192,    144,    144,    112,    112,    48,    112,    112,    112,    112 }
    },{
        // SC
        { 512,    512,     512,     480,     368,    368,    280,    280,    128,    280,    280,    280,    280 },
        { 800,    800,     800,     480,     368,    368,    280,    280,    200,    280,    280,    280,    280 },
        { 1024,   1024,    1024,    1024,    784,    784,    600,    600,    256,    600,    600,    600,    600 },
        { 1536,   1536,    1536,    1024,    784,    784,    600,    600,    384,    600,    600,    600,    600 }
#else
        // NSC
        { 128,    128,     64,     64,    48,    48,      40,     40,     40,     40,     40,     40,     40  },
        { 192,    192,    128,    128,    96,    96,      72,     72,     72,     72,     72,     72,     72  },
        { 240,    240,    192,    192,    144,    144,    112,    112,    112,    112,    112,    112,    112 },
        { 360,    360,    192,    192,    144,    144,    112,    112,    112,    112,    112,    112,    112 }
    },{
        // SC
        {  512,    512,     480,     480,     368,    368,    280,    280,    280,    280,    280,    280,    280 },
        {  800,    800,     480,     480,     368,    368,    280,    280,    280,    280,    280,    280,    280 },
        {  1024,   800 ,    1024,    1024,    784,    784,    600,    600,    600,    600,    600,    600,    600 },
        {  1536,   1536,    1024,    1024,    784,    784,    600,    600,    600,    600,    600,    600,    600 }
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    }
};
    //    M0      M1      M2      M3   M4     M5     M6     M7     M8     M9     M10    M11    M12
static const uint16_t min_me_search_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 100,  100,    58,     58,     56,    56,    56,    56,    56,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 128,  128,    86,     86,     68,    68,    68,    68,    68,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 160,  160,   106,    106,     86,    86,    86,    86,    86,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 200,  200,   128,    128,     86,    86,    86,    86,    86,    48,    48,    48,    48 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 200,   200,   200,   200,   112,   112,   112,  112,     112,    56,    56,    56,    56 },
        { 256,   256,   256,   256,   136,   136,   136,  136,     136,    56,    56,    56,    56 },
        { 320,   320,   320,   320,   168,   168,   168,  168,     168,    56,    56,    56,    56 },
        { 320,   320,   320,   320,   168,   168,   168,  168,     168,    56,    56,    56,    56 }
#else
#if ADOPT_SQ_ME_SEARCH_AREA
        // NSC
        { 114,  114,    58,     58,     56,    56,    56,    56,    56,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 170,  170,    86,     86,     68,    68,    68,    68,    68,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 212,  212,   106,    106,     86,    86,    86,    86,    86,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 256,  256,   128,    128,     86,    86,    86,    86,    86,    48,    48,    48,    48 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 200,   200,   200,    200,   112,   112,   112,  112,     112,    56,    56,    56,    56 },
        { 300,   300,   300,    300,   136,   136,   136,  136,     136,    56,    56,    56,    56 },
        { 400,   400,   400,    400,   168,   168,   168,  168,     168,    56,    56,    56,    56 },
        { 600,   600,   600,    600,   168,   168,   168,  168,     168,    56,    56,    56,    56 }
#else
#if MAR19_ADOPTIONS
        // NSC
        { 160,  160,   128,     128,    80,     80,     80,     80,     80,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 240,  240,   192,     192,    96,     96,     96,     96,     96,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 300,  300,   240,     240,    120,    120,    120,    120,    120,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 360,  360,   360,     360,    120,    120,    120,    120,    120,    48,    48,    48,    48 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 200,   200,   200,    200,   112,   112,   112,  112,     112,    56,    56,    56,    56 },
        { 300,   300,   300,    300,   136,   136,   136,  136,     136,    56,    56,    56,    56 },
        { 400,   400,   400,    400,   168,   168,   168,  168,     168,    56,    56,    56,    56 },
        { 600,   600,   600,    600,   168,   168,   168,  168,     168,    56,    56,    56,    56 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 128,  128,   128,     128,    24,    24,    24,    24,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 192,  192,   192,     192,    48,    48,    48,    48,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 240,  240,   240,     240,    72,    72,    72,    72,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 360,  360,   360,     360,    72,    72,    72,    72,    48,    48,    48,    48,    48 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 168,   168,   168,    168,    72,    72,    72,   72,     56,    56,    56,    56,    56 },
        { 256,   256,   256,    256,    88,    88,    88,   88,     88,    56,    56,    56,    56 },
        { 320,   320,   320,    320,   112,   112,   112,  112,    112,    56,    56,    56,    56 },
        { 480,   480,   480,    480,   168,   168,   168,  168,    168,    56,    56,    56,    56 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 128,  128,   128,     128,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 192,  192,   192,     192,    48,    48,    40,    40,    24,    40,    40,    40,    40 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 240,  240,   240,     240,    72,    72,    56,    56,    32,    56,    56,    56,    56 }, // INPUT_SIZE_1080p_RANGE
        { 360,  360,   360,     360,    72,    72,    56,    56,    48,    56,    56,    56,    56 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 168,   168,   168,    168,    72,    72,    56,    56,    48,    56,    56,    56,    56 },
        { 256,   256,   256,    256,    72,    72,    56,    56,    64,    56,    56,    56,    56 },
        { 320,   320,   320,    320,    72,    72,    56,    56,    80,    56,    56,    56,    56 },
        { 480,   480,   480,    480,    72,    72,    56,    56,    120,   56,    56,    56,    56 }
#else
#if MAR10_ADOPTIONS
        // NSC
        { 128,  128,   128,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 192,  192,   192,    64,    48,    48,    40,    40,    24,    40,    40,    40,    40 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 240,  240,   240,    96,    72,    72,    56,    56,    32,    56,    56,    56,    56 }, // INPUT_SIZE_1080p_RANGE
        { 360,  360,   360,    96,    72,    72,    56,    56,    48,    56,    56,    56,    56 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 168,   168,   168,    96,    72,    72,    56,    56,    48,    56,    56,    56,    56 },
        { 256,   256,   256,    96,    72,    72,    56,    56,    64,    56,    56,    56,    56 },
        { 320,   320,   320,    96,    72,    72,    56,    56,    80,    56,    56,    56,    56 },
        { 480,   480,   480,    96,    72,    72,    56,    56,    120,   56,    56,    56,    56 }
#else
        // NSC
        { 128,    128,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 192,    192,    64,    64,    48,    48,    40,    40,    40,    40,    40,    40,    40 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 240,    240,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 }, // INPUT_SIZE_1080p_RANGE
        { 360,    360,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 168,    168,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 },
        { 256,    256,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 },
        { 320,    256,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 },
        { 480,    480,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 }
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    }
};
static const uint16_t min_me_search_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  75,   75,   75,   75,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  225,  225,  225,  225,  225,  225,  225,  225,  225 }
#else
#if MAR25_ADOPTIONS
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 },
        { 150,  150,  150,  150,  150,   150,   150,  150,   75,   75,   75,   75,   75 }
     } , {
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 },
        { 390,  390,  390,  390,  390,   390,   390,  390,  225,  225,  225,  225,  225 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 100,  100,    58,     58,     56,    56,    56,    56,    56,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 128,  128,    86,     86,     68,    68,    68,    68,    68,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 160,  160,   106,    106,     86,    86,    86,    86,    86,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 200,  200,   128,    128,     86,    86,    86,    86,    86,    48,    48,    48,    48 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 200,   200,   200,   200,   112,   112,   112,  112,     112,    56,    56,    56,    56 },
        { 256,   256,   256,   256,   136,   136,   136,  136,     136,    56,    56,    56,    56 },
        { 320,   320,   320,   320,   168,   168,   168,  168,     168,    56,    56,    56,    56 },
        { 320,   320,   320,   320,   168,   168,   168,  168,     168,    56,    56,    56,    56 }
#else
#if ADOPT_SQ_ME_SEARCH_AREA
        // NSC
        { 114,  114,   58,    58,      56,    56,   56,    56,    56,     8,     8,    8,     8 },
        { 170,  170,   86,    86,      68,    68,   68,    68,    68,    12,    12,   12,    12 },
        { 212,  212,   106,   106,     86,    86,   86,    86,    86,    16,    16,   16,    16 },
        { 256,  256,   128,   128,     86,    86,   86,    86,    86,    24,    24,   24,    24 }
    } , {
        // SC
        { 200,  200,   200,    200,   112,   112,   112,  112,   112,   16,    16,    16,    16 },
        { 300,  300,   300,    300,   136,   136,   136,  136,   136,   40,    40,    40,    40 },
        { 400,  400,   400,    400,   168,   168,   168,  168,   168,   56,    56,    56,    56 },
        { 600,  600,   600,    600,   168,   168,   168,  168,   168,   56,    56,    56,    56 }
#else
#if MAR19_ADOPTIONS
        // NSC
        { 80,   80,     64,     64,     40,    40,   40,    40,    40,     8,     8,    8,     8 },
        { 120,  120,    96,     96,     48,    48,   48,    48,    48,    12,    12,   12,    12 },
        { 150,  150,   120,    120,     60,    60,   60,    60,    60,    16,    16,   16,    16 },
        { 180,  180,   180,    180,     60,    60,   60,    60,    60,    24,    24,   24,    24 }
    } , {
        // SC
        { 200,  200,   200,    200,   112,   112,   112,  112,   112,   16,    16,    16,    16 },
        { 300,  300,   300,    300,   136,   136,   136,  136,   136,   40,    40,    40,    40 },
        { 400,  400,   400,    400,   168,   168,   168,  168,   168,   56,    56,    56,    56 },
        { 600,  600,   600,    600,   168,   168,   168,  168,   168,   56,    56,    56,    56 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 64,   64,    64,     64,     16,    16,   16,    16,     8,     8,     8,    8,     8 },
        { 96,   96,    96,     96,     24,    24,   24,    24,    12,    12,    12,   12,    12 },
        { 120,  120,   120,    120,    40,    40,   40,    40,    16,    16,    16,   16,    16 },
        { 180,  180,   180,    180,    40,    40,   40,    40,    24,    24,    24,   24,    24 }
    } , {
        // SC
        { 168,  168,   168,    168,    72,    72,    72,   72,    56,   16,    16,    16,    16 },
        { 256,  256,   256,    256,    88,    88,    88,   88,    88,   40,    40,    40,    40 },
        { 320,  320,   320,    320,   112,   112,   112,  112,   112,   56,    56,    56,    56 },
        { 480,  480,   480,    480,   168,   168,   168,  168,   168,   56,    56,    56,    56 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 64,   64,    64,     64,     16,    16,    16,    16,    8,     16,    16,    16,    16 },
        { 96,   96,    96,     96,     24,    24,    16,    16,    12,    16,    16,    16,    16 },
        { 120,  120,   120,    120,    40,    40,    32,    32,    16,    32,    32,    32,    32 },
        { 180,  180,   180,    180,    40,    40,    32,    32,    24,    32,    32,    32,    32 }
    } , {
        // SC
        { 168,  168,   168,    168,    24,    24,    16,    16,    48,    16,    16,    16,    16 },
        { 256,  256,   256,    256,    48,    48,    40,    40,    64,    40,    40,    40,    40 },
        { 320,  320,   320,    320,    72,    72,    56,    56,    80,    56,    56,    56,    56 },
        { 480,  480,   480,    480,    72,    72,    56,    56,    120,   56,    56,    56,    56 }
#else
#if MAR10_ADOPTIONS
        // NSC
        { 64,   64,    64,     16,    16,    16,    16,    16,    8,     16,    16,    16,    16 },
        { 96,   96,    96,     32,    24,    24,    16,    16,    12,    16,    16,    16,    16 },
        { 120,  120,   120,    48,    40,    40,    32,    32,    16,    32,    32,    32,    32 },
        { 180,  180,   180,    48,    40,    40,    32,    32,    24,    32,    32,    32,    32 }
    } , {
        // SC
        { 168,  168,   168,    96,    24,    24,    16,    16,    48,    16,    16,    16,    16 },
        { 256,  256,   256,    96,    48,    48,    40,    40,    64,    40,    40,    40,    40 },
        { 320,  320,   320,    96,    72,    72,    56,    56,    80,    56,    56,    56,    56 },
        { 480,  480,   480,    96,    72,    72,    56,    56,    120,   56,    56,    56,    56 }
#else
        // NSC
        {  64,    64,     16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        {  96,    96,     32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 },
        {  120,   120,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 },
        {  180,   180,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 }
    } , {
        // SC
        { 168,   168,    96,    96,    24,    24,    16,    16,    16,    16,    16,    16,    16 },
        { 256,   256,    96,    96,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 320,   256,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 },
        { 480,   480,    96,    96,    72,    72,    56,    56,    56,    56,    56,    56,    56 }
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

    }
};

/******************************************************************************
                          MAX & MIN  ME search region for temporal filtering
*******************************************************************************/
    //    M0      M1      M2      M3   M4     M5     M6    M7     M8     M9     M10   M11    M12
static const uint16_t max_metf_search_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },  // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },  // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },  // INPUT_SIZE_1080p_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }   // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 64,    64,    64,    64,    24,    16,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    80,    80,    32,    24,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    96,    96,    40,    32,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    96,    96,    40,    32,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 60,    60,   60,     60,    48,    40,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,   60,     60,    48,    40,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,  128,    128,   96,    72,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,  128,    128,   96,    72,    72,    72,    72,    72,    72,    72,    72 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 64,    64,    64,    64,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    80,    80,    32,    32,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    96,    96,    40,    40,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    96,    96,    40,    40,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 60,    60,   60,     60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,   60,     60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,  128,    128,   96,    96,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,  128,    128,   96,    96,    72,    72,    72,    72,    72,    72,    72 }
#else
        // NSC
        { 64,    64,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    40,    40,    32,    32,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 60,    60,   60,     60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,   60,     60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,  128,    128,   96,    96,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,  128,    128,   96,    96,    72,    72,    72,    72,    72,    72,    72 }
#endif
#endif
#endif
    }
};
static const uint16_t max_metf_search_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_1080p_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }
#else
#if MAR17_ADOPTIONS
        //NSC
        { 64,    64,    64,   64,    24,    16,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    80,   80,    32,    24,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    96,   96,    40,    32,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    96,   96,    40,    32,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 60,    60,    60,    60,    48,    40,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,    60,    60,    48,    40,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,   128,   128,   96,    72,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,   128,   128,   96,    72,    72,    72,    72,    72,    72,    72,    72 }
#else
#if MAR12_ADOPTIONS
        //NSC
        { 64,    64,    64,   64,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    80,   80,    32,    32,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    96,   96,    40,    40,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    96,   96,    40,    40,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        { 60,    60,    60,    60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,    60,    60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,   128,   128,   96,    96,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,   128,   128,   96,    96,    72,    72,    72,    72,    72,    72,    72 }
#else
    //NSC
        { 64,    64,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        { 80,    80,    40,    40,    32,    32,    24,    24,    24,    24,    24,    24,    24 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        { 96,    96,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 }, // INPUT_SIZE_1080p_RANGE
        { 96,    96,    48,    48,    40,    40,    32,    32,    32,    32,    32,    32,    32 }  // INPUT_SIZE_4K_RANGE
    } , {
    // SC
        { 60,    60,    60,    60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 60,    60,    60,    60,    48,    48,    40,    40,    40,    40,    40,    40,    40 },
        { 128,   128,   128,   128,   96,    96,    72,    72,    72,    72,    72,    72,    72 },
        { 128,   128,   128,   128,   96,    96,    72,    72,    72,    72,    72,    72,    72 }
#endif
#endif
#endif
    }
};
    //    M0      M1      M2      M3   M4     M5     M6    M7     M8     M9     M10   M11    M12
static const uint16_t min_metf_search_width[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 12,    12,    12,    12,    12,    8,     8,     8,     8,     8,     8,     8,     8  },
        { 18,    18,    18,    18,    18,    8,     8,     8,     8,     8,     8,     8,     8  },
        { 24,    24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16 }
    } , {
        // SC
        { 12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12 },
        { 18,    18,    18,    18,    18,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 24,    24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16 }
#else
#if MAR17_ADOPTIONS
        // NSC
        { 16,    16,    16,    16,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 24,    24,    24,    24,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16,    16 }
    } , {
        // SC
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    16,    16,    16,    16,    16,    16,    16,    16 }
#else
#if MAR12_ADOPTIONS
        // NSC
        { 16,    16,    16,    16,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 24,    24,    24,    24,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16,    16 }
    } , {
        // SC
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }
#else
        // NSC
        { 16,   16,    8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 24,   24,    12,    12,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 32,   32,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,   32,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 }
    } , {
        // SC
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    24,    24,    16,    16,    16,    16,    16,    16,    16 }
#endif
#endif
#endif
#endif
    }
};
static const uint16_t min_metf_search_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        // NSC
        { 12,    12,    12,    12,    12,    8,     8,     8,     8,     8,     8,     8,     8  },
        { 18,    18,    18,    18,    18,    8,     8,     8,     8,     8,     8,     8,     8  },
        { 24,    24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16 }
    } , {
        // SC
        { 12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12 },
        { 18,    18,    18,    18,    18,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 24,    24,    24,    24,    24,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 32,    32,    32,    32,    32,    16,    16,    16,    16,    16,    16,    16,    16 }
#else
#if MAR12_ADOPTIONS
        // NSC
        {  8,    8,     8,     8,      4,     4,     4,     4,     4,     4,    4,     4,     4 },
        {  12,   12,    12,    12,     6,     6,     6,     6,     6,     6,    6,     6,     6 },
        {  16,   16,    16,    16,     8,     8,     8,     8,     8,     8,    8,     8,     8 },
        {  16,   16,    16,    16,     8,     8,     8,     8,     8,     8,    8,     8,     8 }
    } , {
        // SC
        { 8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 12,    12,    12,    12,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 }
#else
        // NSC
        {  8,     8,      4,     4,     4,     4,     4,     4,     4,     4,    4,     4,     4 },
        {  12,    12,     6,     6,     6,     6,     6,     6,     6,     6,    6,     6,     6 },
        {  16,    16,     8,     8,     8,     8,     8,     8,     8,     8,    8,     8,     8 },
        {  16,    16,     8,     8,     8,     8,     8,     8,     8,     8,    8,     8,     8 }
    } , {
        // SC
        { 8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 12,    12,    12,    12,    8,     8,     8,     8,     8,     8,     8,     8,     8  },
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 },
        { 16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16,    16 }
#endif
#endif
#endif
    }
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
#if MAR26_ADOPTIONS
        // NSC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_1080p_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }
#else
#if MAR23_ADOPTIONS
        {  48,   48,   48,   48,   48,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  96,   96,   96,   96,   96,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
     } , {
        {  48,   48,   48,   48,   48,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  96,   96,   96,   96,   96,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
#else
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  112,  112,  112,  112,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 },
     } , {
        {  48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 112,  112,  112,  112,  112,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   48,   48,   48,   48,   48,   48,   48,   48 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 },
#endif
#endif
    }
};

static const uint16_t tf_hme_level0_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#else
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#endif
#endif
    }
};
static const uint16_t tf_hme_level0_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#else
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  56,   56,   56,   56,   56,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#endif
#endif
    }
};
static const uint16_t tf_hme_level0_total_search_area_height[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }, // INPUT_SIZE_1080p_RANGE
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 },
        {  32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,  32 }
#else
#if MAR23_ADOPTIONS
        {  48,   48,   48,   48,   48,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  96,   96,   96,   96,   96,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
     } , {
        {  48,   48,   48,   48,   48,   40,   40,   40,   40,   40,   40,   40,   40 },
        {  96,   96,   96,   96,   96,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   40,   40,   40,   40,   40,   40,   40,   40 },
        { 128,  128,  128,  128,  128,   96,   96,   96,   96,   96,   96,   96,   96 }
#else
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
    } , {
        {  40,   40,   40,   40,   40,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  64,   64,   64,   64,   64,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 },
        {  80,   80,   80,   80,   80,   24,   24,   24,   24,   24,   24,   24,   24 }
#endif
#endif
    }
};
static const uint16_t tf_hme_level0_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#else
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    } , {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
#endif
#endif
    }
};
static const uint16_t tf_hme_level0_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR26_ADOPTIONS
        // NSC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_576p_RANGE_OR_LOWER
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }, // INPUT_SIZE_1080p_RANGE
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }  // INPUT_SIZE_4K_RANGE
    } , {
        // SC
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 },
        {  16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16 }
#else
#if MAR23_ADOPTIONS
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
    } , {
        {  24,   24,   24,   24,   24,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  48,   48,   48,   48,   48,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   20,   20,   20,   20,   20,   20,   20,   20 },
        {  64,   64,   64,   64,   64,   48,   48,   48,   48,   48,   48,   48,   48 }
#else
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
    }, {
        {  20,   20,   20,   20,   20,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  32,   32,   32,   32,   32,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 },
        {  40,   40,   40,   40,   40,   12,   12,   12,   12,   12,   12,   12,   12 }
#endif
#endif
    }
};

// HME LEVEL 1
   //      M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_level1_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR23_ADOPTIONS
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }, {
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#else
#if MAR17_ADOPTIONS
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }, {
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#else
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    }, {
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#endif
#endif
    }
};
static const uint16_t tf_hme_level1_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
#endif
    }
};
static const uint16_t tf_hme_level1_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
#endif
    }
};
static const uint16_t tf_hme_level1_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
#endif
    }
};
static const uint16_t tf_hme_level1_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
    } , {
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 },
        {  16,   16,   16,   16,   16,    8,    8,    8,    8,    8,    8,    8,     8 }
#endif
    }
};
// HME LEVEL 2
    //     M0    M1    M2    M3    M4    M5    M6    M7    M8    M9    M10    M11    M12
static const uint8_t tf_enable_hme_level2_flag[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR23_ADOPTIONS
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#else
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
    },{
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_576p_RANGE_OR_LOWER
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_720P_RANGE/INPUT_SIZE_1080i_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 },      // INPUT_SIZE_1080p_RANGE
        {   1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,     0,    0 }       // INPUT_SIZE_4K_RANGE
#endif
    }
};
static const uint16_t tf_hme_level2_search_area_in_width_array_left[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
    }
};
static const uint16_t tf_hme_level2_search_area_in_width_array_right[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
    }
};
static const uint16_t tf_hme_level2_search_area_in_height_array_top[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
    }
};
static const uint16_t tf_hme_level2_search_area_in_height_array_bottom[SC_MAX_LEVEL][INPUT_SIZE_COUNT][MAX_SUPPORTED_MODES] = {
    {
#if MAR25_ADOPTIONS
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
    } , {
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 },
        {  16,   16,   16,   16,   16,    16,   16,   16,   16,    8,    8,    8,  8 }
#else
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
    } , {
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 },
        {   8,    8,    8,    8,    8,    4,    4,    4,    4,    4,    4,     4,    4 }
#endif
    }
};

#endif
#if !MAR17_ADOPTIONS
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

#endif
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
typedef struct _EbEncHandle EbEncHandle;
typedef struct _EbThreadContext EbThreadContext;
#ifdef __cplusplus
}
#endif
#endif // EbDefinitions_h
/* File EOF */
