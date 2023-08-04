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

#ifndef EbUtility_h
#define EbUtility_h

#include "EbDefinitions.h"
#include "common_dsp_rtcd.h"
#ifdef __cplusplus
extern "C" {
#endif
#include <limits.h>
/****************************
     * UTILITY FUNCTIONS
     ****************************/
typedef struct BlockList {
    uint8_t  list_size;
    uint16_t blk_mds_table[3]; //stores a max of 3 redundant blocks
} BlockList_t;
typedef enum GeomIndex {
    GEOM_0, //64x64  ->8x8  NSQ:OFF
    GEOM_1, //64x64  ->8x8  NSQ:ON (only H & V shapes, but not 8x4 and 4x8 and not 16x8 and 8x16)
    GEOM_2, //64x64  ->8x8  NSQ:ON (only H & V shapes, but not 8x4 and 4x8)
    GEOM_3, //64x64  ->8x8  NSQ:ON (only H & V shapes)
    GEOM_4, //64x64  ->4x4  NSQ:ON (only H & V shapes)
    GEOM_5, //64x64  ->4x4  NSQ:ON (all shapes)
    GEOM_6, //128x128->4x4  NSQ:ON (all shapes
    GEOM_TOT
} GeomIndex;
void svt_aom_build_blk_geom(GeomIndex geom);

typedef struct BlockGeom {
    Part    shape; // P_N..P_V4 . P_S is not used.
    uint8_t org_x; // orgin x from topleft of sb
    uint8_t org_y; // orgin x from topleft of sb

    uint8_t has_uv;
    uint8_t sq_size; // size of parent square

    uint8_t   bwidth; // block width
    uint8_t   bheight; // block height
    uint8_t   bwidth_uv; // block width for Chroma 4:2:0
    uint8_t   bheight_uv; // block height for Chroma 4:2:0
    BlockSize bsize; // bloc size
    BlockSize bsize_uv; // bloc size for Chroma 4:2:0
    uint8_t   txb_count[MAX_VARTX_DEPTH + 1]; //4-2-1
    TxSize    txsize[MAX_VARTX_DEPTH + 1];
    TxSize    txsize_uv[MAX_VARTX_DEPTH + 1];
    uint8_t   tx_width[MAX_VARTX_DEPTH + 1]; //tx_size_wide
    uint8_t   tx_height[MAX_VARTX_DEPTH + 1]; //tx_size_high
    uint8_t   tx_width_uv[MAX_VARTX_DEPTH + 1]; //tx_size_wide
    uint8_t   tx_height_uv[MAX_VARTX_DEPTH + 1]; //tx_size_high

    //origin is SB - separate tables for INTRA (idx 0) and INTER (idx 1)
    uint8_t tx_org_x[2][MAX_VARTX_DEPTH + 1][MAX_TXB_COUNT];
    uint8_t tx_org_y[2][MAX_VARTX_DEPTH + 1][MAX_TXB_COUNT];

    uint16_t blkidx_mds; // block index in md scan
    // index of the block in d1 dimension 0..24  (0 is parent square, 1 top half of H , ...., 24:last quarter of V4)
    uint8_t  d1i;
    uint16_t sqi_mds; // index of the parent square in md  scan.
    uint16_t parent_depth_idx_mds; // index of the parent block of a given depth
    // max number of ns blocks within one partition 1..4 (N:1,H:2,V:2,HA:3,HB:3,VA:3,VB:3,H4:4,V4:4)
    uint8_t     totns;
    uint8_t     nsi; // non square index within a partition  0..totns-1
    uint8_t     quadi; // parent square is in which quadrant 0..3
    GeomIndex   svt_aom_geom_idx; //type of geom this block belongs
    uint8_t     depth; // depth of the block
    uint16_t    d1_depth_offset; // offset to the next d1 sq block
    uint16_t    ns_depth_offset; // offset to the next nsq block (skip remaining d2 blocks)
    uint8_t     is_last_quadrant; // only for square bloks, is this the fourth quadrant block?
    uint8_t     redund; // 1: means that this block is redundant to another
    BlockList_t redund_list; // the list where the block is redundant
} BlockGeom;

static const BlockSize ss_size_lookup[BlockSizeS_ALL][2][2] = {
    //  ss_x == 0    ss_x == 0        ss_x == 1      ss_x == 1
    //  ss_y == 0    ss_y == 1        ss_y == 0      ss_y == 1
    {{BLOCK_4X4, BLOCK_4X4}, {BLOCK_4X4, BLOCK_4X4}},
    {{BLOCK_4X8, BLOCK_4X4}, {BLOCK_INVALID, BLOCK_4X4}},
    {{BLOCK_8X4, BLOCK_INVALID}, {BLOCK_4X4, BLOCK_4X4}},
    {{BLOCK_8X8, BLOCK_8X4}, {BLOCK_4X8, BLOCK_4X4}},
    {{BLOCK_8X16, BLOCK_8X8}, {BLOCK_INVALID, BLOCK_4X8}},
    {{BLOCK_16X8, BLOCK_INVALID}, {BLOCK_8X8, BLOCK_8X4}},
    {{BLOCK_16X16, BLOCK_16X8}, {BLOCK_8X16, BLOCK_8X8}},
    {{BLOCK_16X32, BLOCK_16X16}, {BLOCK_INVALID, BLOCK_8X16}},
    {{BLOCK_32X16, BLOCK_INVALID}, {BLOCK_16X16, BLOCK_16X8}},
    {{BLOCK_32X32, BLOCK_32X16}, {BLOCK_16X32, BLOCK_16X16}},
    {{BLOCK_32X64, BLOCK_32X32}, {BLOCK_INVALID, BLOCK_16X32}},
    {{BLOCK_64X32, BLOCK_INVALID}, {BLOCK_32X32, BLOCK_32X16}},
    {{BLOCK_64X64, BLOCK_64X32}, {BLOCK_32X64, BLOCK_32X32}},
    {{BLOCK_64X128, BLOCK_64X64}, {BLOCK_INVALID, BLOCK_32X64}},
    {{BLOCK_128X64, BLOCK_INVALID}, {BLOCK_64X64, BLOCK_64X32}},
    {{BLOCK_128X128, BLOCK_128X64}, {BLOCK_64X128, BLOCK_64X64}},
    {{BLOCK_4X16, BLOCK_4X8}, {BLOCK_INVALID, BLOCK_4X8}},
    {{BLOCK_16X4, BLOCK_INVALID}, {BLOCK_8X4, BLOCK_8X4}},
    {{BLOCK_8X32, BLOCK_8X16}, {BLOCK_INVALID, BLOCK_4X16}},
    {{BLOCK_32X8, BLOCK_INVALID}, {BLOCK_16X8, BLOCK_16X4}},
    {{BLOCK_16X64, BLOCK_16X32}, {BLOCK_INVALID, BLOCK_8X32}},
    {{BLOCK_64X16, BLOCK_INVALID}, {BLOCK_32X16, BLOCK_32X8}}};
static INLINE BlockSize get_plane_block_size(BlockSize bsize, int32_t subsampling_x, int32_t subsampling_y) {
    if (bsize == BLOCK_INVALID)
        return BLOCK_INVALID;
    return ss_size_lookup[bsize][subsampling_x][subsampling_y];
}

static INLINE TxSize av1_get_max_uv_txsize(BlockSize bsize, int32_t subsampling_x, int32_t subsampling_y) {
    const BlockSize plane_bsize = get_plane_block_size(bsize, subsampling_x, subsampling_y);
    TxSize          uv_tx       = TX_INVALID;
    if (plane_bsize < BlockSizeS_ALL)
        uv_tx = max_txsize_rect_lookup[plane_bsize];
    return av1_get_adjusted_tx_size(uv_tx);
}

#define NOT_USED_VALUE 0
//gives the index of parent from the last qudrant child
static const uint32_t parent_depth_offset[GEOM_TOT][6] = {{NOT_USED_VALUE, 64, 16, 4, NOT_USED_VALUE, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 80, 20, 4, NOT_USED_VALUE, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 128, 32, 8, NOT_USED_VALUE, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 320, 80, 20, NOT_USED_VALUE, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 512, 128, 32, 8, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 832, 208, 52, 8, NOT_USED_VALUE},
                                                          {NOT_USED_VALUE, 3320, 832, 208, 52, 8}};
//gives the index of next quadrant child within a depth
static const uint32_t ns_depth_offset[GEOM_TOT][6] = {{85, 21, 5, 1, NOT_USED_VALUE, NOT_USED_VALUE},
                                                      {105, 25, 5, 1, NOT_USED_VALUE, NOT_USED_VALUE},
                                                      {169, 41, 9, 1, NOT_USED_VALUE, NOT_USED_VALUE},
                                                      {425, 105, 25, 5, NOT_USED_VALUE, NOT_USED_VALUE},
                                                      {681, 169, 41, 9, 1, NOT_USED_VALUE},
                                                      {1101, 269, 61, 9, 1, NOT_USED_VALUE},
                                                      {4421, 1101, 269, 61, 9, 1}};
//gives the next depth block(first qudrant child) from a given parent square
static const uint32_t d1_depth_offset[GEOM_TOT][6] = {{1, 1, 1, 1, 1, NOT_USED_VALUE},
                                                      {5, 5, 1, 1, 1, NOT_USED_VALUE},
                                                      {5, 5, 5, 1, 1, NOT_USED_VALUE},
                                                      {5, 5, 5, 5, 1, NOT_USED_VALUE},
                                                      {5, 5, 5, 5, 1, NOT_USED_VALUE},
                                                      {25, 25, 25, 5, 1, NOT_USED_VALUE},
                                                      {17, 25, 25, 25, 5, 1}};
extern BlockGeom      svt_aom_blk_geom_mds[MAX_NUM_BLOCKS_ALLOC];

static INLINE const BlockGeom* get_blk_geom_mds(uint32_t bidx_mds) { return &svt_aom_blk_geom_mds[bidx_mds]; }
// CU Stats Helper Functions
typedef struct CodedBlockStats {
    uint8_t depth;
    uint8_t size;
    uint8_t size_log2;
    uint8_t org_x;
    uint8_t org_y;
    uint8_t cu_num_in_depth;
    uint8_t parent32x32_index;
} CodedBlockStats;

extern void*                  svt_aom_memalign(size_t align, size_t size);
extern void*                  svt_aom_malloc(size_t size);
extern void                   svt_aom_free(void* memblk);
extern void*                  svt_aom_memset16(void* dest, int32_t val, size_t length);
extern const CodedBlockStats* svt_aom_get_coded_blk_stats(const uint32_t cu_idx);

#define PU_ORIGIN_ADJUST(cu_origin, cu_size, offset) ((((cu_size) * (offset)) >> 2) + (cu_origin))
#define PU_SIZE_ADJUST(cu_size, puSize) (((cu_size) * (puSize)) >> 2)

#define TU_ORIGIN_ADJUST(cu_origin, cu_size, offset) ((((cu_size) * (offset)) >> 2) + (cu_origin))
#define TU_SIZE_ADJUST(cu_size, tuDepth) ((cu_size) >> (tuDepth))

/****************************
     * MACROS
     ****************************/

#define MULTI_LINE_MACRO_BEGIN do {
#define MULTI_LINE_MACRO_END \
    }                        \
    while (0)

//**************************************************
// MACROS
//**************************************************
#define DIVIDE_AND_ROUND(x, y) (((x) + ((y) >> 1)) / (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MEDIAN(a, b, c)                   ((a)>(b)?(a)>?(b)>?(b)::(a):(b)>?(a)>?(a)::(b))
#define CLIP3(min_val, max_val, a) (((a) < (min_val)) ? (min_val) : (((a) > (max_val)) ? (max_val) : (a)))
#define CLIP3EQ(min_val, max_val, a) (((a) <= (min_val)) ? (min_val) : (((a) >= (max_val)) ? (max_val) : (a)))
#define BITDEPTH_MIDRANGE_VALUE(precision) (1 << ((precision)-1))
#define SWAP(a, b)                    \
    MULTI_LINE_MACRO_BEGIN(a) ^= (b); \
    (b) ^= (a);                       \
    (a) ^= (b);                       \
    MULTI_LINE_MACRO_END
#define ABS(a) (((a) < 0) ? (-(a)) : (a))
#define EB_ABS_DIFF(a, b) ((a) > (b) ? ((a) - (b)) : ((b) - (a)))
#define EB_DIFF_SQR(a, b) (((a) - (b)) * ((a) - (b)))
#define SQR(x) ((x) * (x))
#define POW2(x) (1 << (x))
#define SIGN(a, b) (((a - b) < 0) ? (-1) : ((a - b) > 0) ? (1) : 0)
#define ROUND(a) (a >= 0) ? (a + 1 / 2) : (a - 1 / 2);
#define UNSIGNED_DEC(x)                                    \
    MULTI_LINE_MACRO_BEGIN(x) = (((x) > 0) ? ((x)-1) : 0); \
    MULTI_LINE_MACRO_END
#define CIRCULAR_ADD(x, max) (((x) >= (max)) ? ((x) - (max)) : ((x) < 0) ? ((max) + (x)) : (x))
#define CIRCULAR_ADD_UNSIGNED(x, max) (((x) >= (max)) ? ((x) - (max)) : (x))
#define CEILING(x, base) ((((x) + (base)-1) / (base)) * (base))
#define POW2_CHECK(x) ((x) == ((x) & (-((int32_t)(x)))))
#define ROUND_UP_MUL_8(x) ((x) + ((8 - ((x)&0x7)) & 0x7))
#define ROUND_UP_MULT(x, mult) ((x) + (((mult) - ((x) & ((mult)-1))) & ((mult)-1)))

// rounds down to the next power of two
#define FLOOR_POW2(x)                        \
    MULTI_LINE_MACRO_BEGIN(x) |= ((x) >> 1); \
    (x) |= ((x) >> 2);                       \
    (x) |= ((x) >> 4);                       \
    (x) |= ((x) >> 8);                       \
    (x) |= ((x) >> 16);                      \
    (x) -= ((x) >> 1);                       \
    MULTI_LINE_MACRO_END

// rounds up to the next power of two
#define CEIL_POW2(x)                \
    MULTI_LINE_MACRO_BEGIN(x) -= 1; \
    (x) |= ((x) >> 1);              \
    (x) |= ((x) >> 2);              \
    (x) |= ((x) >> 4);              \
    (x) |= ((x) >> 8);              \
    (x) |= ((x) >> 16);             \
    (x) += 1;                       \
    MULTI_LINE_MACRO_END
#define LOG2F_8(x)              \
    (((x) < 0x0002u)       ? 0u \
         : ((x) < 0x0004u) ? 1u \
         : ((x) < 0x0008u) ? 2u \
         : ((x) < 0x0010u) ? 3u \
         : ((x) < 0x0020u) ? 4u \
         : ((x) < 0x0040u) ? 5u \
                           : 6u)

#define TWO_D_INDEX(x, y, stride) (((y) * (stride)) + (x))

// MAX_CU_COUNT is used to find the total number of partitions for the max partition depth and for
// each parent partition up to the root partition level (i.e. SB level).

// MAX_CU_COUNT is given by SUM from k=1 to n (4^(k-1)), reduces by using the following finite sum
// SUM from k=1 to n (q^(k-1)) = (q^n - 1)/(q-1) => (4^n - 1) / 3
#define MAX_CU_COUNT(max_depth_count) ((((1 << (max_depth_count)) * (1 << (max_depth_count))) - 1) / 3)

//**************************************************
// CONSTANTS
//**************************************************
#define MIN_UNSIGNED_VALUE 0
#define MAX_UNSIGNED_VALUE ~0u
#define MIN_SIGNED_VALUE ~0 - ((signed)(~0u >> 1))
#define MAX_SIGNED_VALUE ((signed)(~0u >> 1))
#define CONST_SQRT2 (1.4142135623730950488016887242097) /*sqrt(2)*/

// Helper functions for EbLinkedListNode.

// concatenate two linked list, and return the pointer to the new concatenated list
EbLinkedListNode* svt_aom_concat_eb_linked_list(EbLinkedListNode* a, EbLinkedListNode* b);

// split a linked list into two. return the pointer to a linked list whose nodes meets the condition
// predicate_func(node) == TRUE, the rest of the nodes will be collected into another linked list to which (*restLL) is
// set. Does not gaurantee the original order of the nodes.

EbLinkedListNode* svt_aom_split_eb_linked_list(EbLinkedListNode* input, EbLinkedListNode** restLL,
                                               Bool (*predicate_func)(EbLinkedListNode*));

#if OPT_ENABLE_2L_INCOMP
#define MINI_GOP_MAX_COUNT 31
#define MINI_GOP_WINDOW_MAX_COUNT 16 // window subdivision: 16 x 2L

#define MIN_HIERARCHICAL_LEVEL 1
static const uint8_t mini_gop_offset[MAX_HIERARCHICAL_LEVEL - MIN_HIERARCHICAL_LEVEL] = {1, 3, 7, 15, 31};

typedef struct MiniGopStats {
    uint8_t hierarchical_levels;
    uint8_t start_index;
    uint8_t end_index;
    uint8_t length;
} MiniGopStats;
extern const MiniGopStats* svt_aom_get_mini_gop_stats(const uint32_t mini_gop_index);
typedef enum MinigopIndex {
    L6_INDEX    = 0,
    L5_0_INDEX  = 1,
    L4_0_INDEX  = 2,
    L3_0_INDEX  = 3,
    L2_0_INDEX  = 4,
    L2_1_INDEX  = 5,
    L3_1_INDEX  = 6,
    L2_2_INDEX  = 7,
    L2_3_INDEX  = 8,
    L4_1_INDEX  = 9,
    L3_2_INDEX  = 10,
    L2_4_INDEX  = 11,
    L2_5_INDEX  = 12,
    L3_3_INDEX  = 13,
    L2_6_INDEX  = 14,
    L2_7_INDEX  = 15,
    L5_1_INDEX  = 16,
    L4_2_INDEX  = 17,
    L3_4_INDEX  = 18,
    L2_8_INDEX  = 19,
    L2_9_INDEX  = 20,
    L3_5_INDEX  = 21,
    L2_10_INDEX = 22,
    L2_11_INDEX = 23,
    L4_3_INDEX  = 24,
    L3_6_INDEX  = 25,
    L2_12_INDEX = 26,
    L2_13_INDEX = 27,
    L3_7_INDEX  = 28,
    L2_14_INDEX = 29,
    L2_15_INDEX = 30
} MinigopIndex;
#else
#define MINI_GOP_MAX_COUNT 15
#define MINI_GOP_WINDOW_MAX_COUNT 8 // widow subdivision: 8 x 3L

#define MIN_HIERARCHICAL_LEVEL 2
static const uint32_t mini_gop_offset[4] = {1, 3, 7, 31};

typedef struct MiniGopStats {
    uint32_t hierarchical_levels;
    uint32_t start_index;
    uint32_t end_index;
    uint32_t lenght;
} MiniGopStats;
extern const MiniGopStats* svt_aom_get_mini_gop_stats(const uint32_t mini_gop_index);
typedef enum MinigopIndex {
    L6_INDEX   = 0,
    L5_0_INDEX = 1,
    L4_0_INDEX = 2,
    L3_0_INDEX = 3,
    L3_1_INDEX = 4,
    L4_1_INDEX = 5,
    L3_2_INDEX = 6,
    L3_3_INDEX = 7,
    L5_1_INDEX = 8,
    L4_2_INDEX = 9,
    L3_4_INDEX = 10,
    L3_5_INDEX = 11,
    L4_3_INDEX = 12,
    L3_6_INDEX = 13,
    L3_7_INDEX = 14
} MinigopIndex;
#endif
// Right shift that replicates gcc's implementation

static inline int gcc_right_shift(int a, unsigned shift) {
    if (!a)
        return 0;
    if (a > 0)
        return a >> shift;
    static const unsigned sbit = 1u << (sizeof(sbit) * CHAR_BIT - 1);
    a                          = (unsigned)a >> shift;
    while (shift) a |= sbit >> shift--;
    return a ^ sbit;
}

static INLINE int convert_to_trans_prec(int allow_hp, int coor) {
    if (allow_hp)
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 3);
    else
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 2) * 2;
}

/* Convert Floating Point to Fixed Point example: int32_t val_fp8 = FLOAT2FP(val_float, 8, int32_t) */
#define FLOAT2FP(x_float, base_move, fp_type) ((fp_type)((x_float) * (((fp_type)(1)) << (base_move))))

/* Convert Fixed Point to Floating Point example: double val = FP2FLOAT(val_fp8, 8, int32_t, double) */
#define FP2FLOAT(x_fp, base_move, fp_type, float_type) \
    ((((float_type)((fp_type)(x_fp))) / ((float_type)(((fp_type)1) << (base_move)))))

#ifndef FIXED_POINT_ASSERT_TEST
#if NDEBUG
#define FIXED_POINT_ASSERT_TEST 0
#else
#define FIXED_POINT_ASSERT_TEST 1
#endif
#endif

#if FIXED_POINT_ASSERT_TEST
void svt_fixed_point_test_breakpoint(char* file, unsigned line);
#define FP_ASSERT(expression)                                                                             \
    if (!(expression)) {                                                                                  \
        fprintf(stderr, "ERROR: FP_ASSERT Fixed Point overload %s:%u\n", __FILE__, (unsigned)(__LINE__)); \
        svt_fixed_point_test_breakpoint(__FILE__, (unsigned)(__LINE__));                                  \
        assert(0);                                                                                        \
        abort();                                                                                          \
    }
#else /*FIXED_POINT_ASSERT_TEST*/
#define FP_ASSERT(expression)
#endif

#ifdef __cplusplus
}
#endif

#endif // EbUtility_h
