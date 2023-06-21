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

#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <stdlib.h>
#include <sys/time.h>
#endif

#include "EbUtility.h"
#include "EbLog.h"
#include <math.h>

/* assert a certain condition and report err if condition not met */
void svt_aom_assert_err(uint32_t condition, char* err_msg) {
    assert(condition);
    if (!condition)
        SVT_ERROR("\n %s \n", err_msg);
}

/********************************************************************************************
* faster memcopy for <= 64B blocks, great w/ inlining and size known at compile time (or w/ PGO)
* THIS NEEDS TO STAY IN A HEADER FOR BEST PERFORMANCE
********************************************************************************************/
#ifdef ARCH_X86_64
#include <immintrin.h>
#if defined(__GNUC__) && !defined(__clang__) && !defined(__ICC__)
__attribute__((optimize("unroll-loops")))
#endif
static void
svt_memcpy_small(void* dst_ptr, const void* src_ptr, size_t size) {
    const unsigned char* src = src_ptr;
    unsigned char*       dst = dst_ptr;
    size_t               i   = 0;

#ifdef _INTEL_COMPILER
#pragma unroll
#endif
    while ((i + 16) <= size) {
        _mm_storeu_ps((float*)(void*)(dst + i), _mm_loadu_ps((const float*)(const void*)(src + i)));
        i += 16;
    }

    if ((i + 8) <= size) {
        _mm_store_sd((double*)(void*)(dst + i), _mm_load_sd((const double*)(const void*)(src + i)));
        i += 8;
    }

    for (; i < size; ++i) dst[i] = src[i];
}
#define EB_MIN(a, b) (((a) < (b)) ? (a) : (b))
static void svt_memcpy_sse(void* dst_ptr, const void* src_ptr, size_t size) {
    const unsigned char* src       = src_ptr;
    unsigned char*       dst       = dst_ptr;
    size_t               i         = 0;
    size_t               align_cnt = EB_MIN((64 - ((size_t)dst & 63)), size);

    // align dest to a $line
    if (align_cnt != 64) {
        svt_memcpy_small(dst, src, align_cnt);
        dst += align_cnt;
        src += align_cnt;
        size -= align_cnt;
    }

    // copy a $line at a time
    // dst aligned to a $line
    size_t cline_cnt = (size & ~(size_t)63);
    for (i = 0; i < cline_cnt; i += 64) {
        __m128 c0 = _mm_loadu_ps((const float*)(const void*)(src + i));
        __m128 c1 = _mm_loadu_ps((const float*)(const void*)(src + i + sizeof(c0)));
        __m128 c2 = _mm_loadu_ps((const float*)(const void*)(src + i + sizeof(c0) * 2));
        __m128 c3 = _mm_loadu_ps((const float*)(const void*)(src + i + sizeof(c0) * 3));

        _mm_storeu_ps((float*)(void*)(dst + i), c0);
        _mm_storeu_ps((float*)(void*)(dst + i + sizeof(c0)), c1);
        _mm_storeu_ps((float*)(void*)(dst + i + sizeof(c0) * 2), c2);
        _mm_storeu_ps((float*)(void*)(dst + i + sizeof(c0) * 3), c3);
    }

    // copy the remainder
    if (i < size)
        svt_memcpy_small(dst + i, src + i, size - i);
}
void svt_memcpy_app(void* dst_ptr, const void* src_ptr, size_t size) {
    if (size > 64)
        svt_memcpy_sse(dst_ptr, src_ptr, size);
    else
        svt_memcpy_small(dst_ptr, src_ptr, size);
}
#endif
/*****************************************
 * Z-Order
 *****************************************/
static TxSize blocksize_to_txsize[BlockSizeS_ALL] = {
    TX_4X4, // BLOCK_4X4
    TX_4X8, // BLOCK_4X8
    TX_8X4, // BLOCK_8X4
    TX_8X8, // BLOCK_8X8
    TX_8X16, // BLOCK_8X16
    TX_16X8, // BLOCK_16X8
    TX_16X16, // BLOCK_16X16
    TX_16X32, // BLOCK_16X32
    TX_32X16, // BLOCK_32X16
    TX_32X32, // BLOCK_32X32
    TX_32X64, // BLOCK_32X64
    TX_64X32, // BLOCK_64X32
    TX_64X64, // BLOCK_64X64
    TX_64X64, // BLOCK_64X128
    TX_64X64, // BLOCK_128X64
    TX_64X64, // BLOCK_128X128
    TX_4X16, // BLOCK_4X16
    TX_16X4, // BLOCK_16X4
    TX_8X32, // BLOCK_8X32
    TX_32X8, // BLOCK_32X8
    TX_16X64, // BLOCK_16X64
    TX_64X16 // BLOCK_64X16
};

static CodedBlockStats coded_unit_stats_array[] = {
    //   Depth       Size      SizeLog2     OriginX    OriginY   cu_num_in_depth   Index
    {0, 64, 6, 0, 0, 0, 0}, // 0
    {1, 32, 5, 0, 0, 0, 1}, // 1
    {2, 16, 4, 0, 0, 0, 1}, // 2
    {3, 8, 3, 0, 0, 0, 1}, // 3
    {3, 8, 3, 8, 0, 1, 1}, // 4
    {3, 8, 3, 0, 8, 8, 1}, // 5
    {3, 8, 3, 8, 8, 9, 1}, // 6
    {2, 16, 4, 16, 0, 1, 1}, // 7
    {3, 8, 3, 16, 0, 2, 1}, // 8
    {3, 8, 3, 24, 0, 3, 1}, // 9
    {3, 8, 3, 16, 8, 10, 1}, // 10
    {3, 8, 3, 24, 8, 11, 1}, // 11
    {2, 16, 4, 0, 16, 4, 1}, // 12
    {3, 8, 3, 0, 16, 16, 1}, // 13
    {3, 8, 3, 8, 16, 17, 1}, // 14
    {3, 8, 3, 0, 24, 24, 1}, // 15
    {3, 8, 3, 8, 24, 25, 1}, // 16
    {2, 16, 4, 16, 16, 5, 1}, // 17
    {3, 8, 3, 16, 16, 18, 1}, // 18
    {3, 8, 3, 24, 16, 19, 1}, // 19
    {3, 8, 3, 16, 24, 26, 1}, // 20
    {3, 8, 3, 24, 24, 27, 1}, // 21
    {1, 32, 5, 32, 0, 1, 2}, // 22
    {2, 16, 4, 32, 0, 2, 2}, // 23
    {3, 8, 3, 32, 0, 4, 2}, // 24
    {3, 8, 3, 40, 0, 5, 2}, // 25
    {3, 8, 3, 32, 8, 12, 2}, // 26
    {3, 8, 3, 40, 8, 13, 2}, // 27
    {2, 16, 4, 48, 0, 3, 2}, // 28
    {3, 8, 3, 48, 0, 6, 2}, // 29
    {3, 8, 3, 56, 0, 7, 2}, // 30
    {3, 8, 3, 48, 8, 14, 2}, // 31
    {3, 8, 3, 56, 8, 15, 2}, // 32
    {2, 16, 4, 32, 16, 6, 2}, // 33
    {3, 8, 3, 32, 16, 20, 2}, // 34
    {3, 8, 3, 40, 16, 21, 2}, // 35
    {3, 8, 3, 32, 24, 28, 2}, // 36
    {3, 8, 3, 40, 24, 29, 2}, // 37
    {2, 16, 4, 48, 16, 7, 2}, // 38
    {3, 8, 3, 48, 16, 22, 2}, // 39
    {3, 8, 3, 56, 16, 23, 2}, // 40
    {3, 8, 3, 48, 24, 30, 2}, // 41
    {3, 8, 3, 56, 24, 31, 2}, // 42
    {1, 32, 5, 0, 32, 2, 3}, // 43
    {2, 16, 4, 0, 32, 8, 3}, // 44
    {3, 8, 3, 0, 32, 32, 3}, // 45
    {3, 8, 3, 8, 32, 33, 3}, // 46
    {3, 8, 3, 0, 40, 40, 3}, // 47
    {3, 8, 3, 8, 40, 41, 3}, // 48
    {2, 16, 4, 16, 32, 9, 3}, // 49
    {3, 8, 3, 16, 32, 34, 3}, // 50
    {3, 8, 3, 24, 32, 35, 3}, // 51
    {3, 8, 3, 16, 40, 42, 3}, // 52
    {3, 8, 3, 24, 40, 43, 3}, // 53
    {2, 16, 4, 0, 48, 12, 3}, // 54
    {3, 8, 3, 0, 48, 48, 3}, // 55
    {3, 8, 3, 8, 48, 49, 3}, // 56
    {3, 8, 3, 0, 56, 56, 3}, // 57
    {3, 8, 3, 8, 56, 57, 3}, // 58
    {2, 16, 4, 16, 48, 13, 3}, // 59
    {3, 8, 3, 16, 48, 50, 3}, // 60
    {3, 8, 3, 24, 48, 51, 3}, // 61
    {3, 8, 3, 16, 56, 58, 3}, // 62
    {3, 8, 3, 24, 56, 59, 3}, // 63
    {1, 32, 5, 32, 32, 3, 4}, // 64
    {2, 16, 4, 32, 32, 10, 4}, // 65
    {3, 8, 3, 32, 32, 36, 4}, // 66
    {3, 8, 3, 40, 32, 37, 4}, // 67
    {3, 8, 3, 32, 40, 44, 4}, // 68
    {3, 8, 3, 40, 40, 45, 4}, // 69
    {2, 16, 4, 48, 32, 11, 4}, // 70
    {3, 8, 3, 48, 32, 38, 4}, // 71
    {3, 8, 3, 56, 32, 39, 4}, // 72
    {3, 8, 3, 48, 40, 46, 4}, // 73
    {3, 8, 3, 56, 40, 47, 4}, // 74
    {2, 16, 4, 32, 48, 14, 4}, // 75
    {3, 8, 3, 32, 48, 52, 4}, // 76
    {3, 8, 3, 40, 48, 53, 4}, // 77
    {3, 8, 3, 32, 56, 60, 4}, // 78
    {3, 8, 3, 40, 56, 61, 4}, // 79
    {2, 16, 4, 48, 48, 15, 4}, // 80
    {3, 8, 3, 48, 48, 54, 4}, // 81
    {3, 8, 3, 56, 48, 55, 4}, // 82
    {3, 8, 3, 48, 56, 62, 4}, // 83
    {3, 8, 3, 56, 56, 63, 4} // 84
};

/**************************************************************
 * Get Coded Unit Statistics
 **************************************************************/
const CodedBlockStats* svt_aom_get_coded_blk_stats(const uint32_t cu_idx) { return &coded_unit_stats_array[cu_idx]; }

/*****************************************
  * Long Log 2
  *  This is a quick adaptation of a Number
  *  Leading Zeros (NLZ) algorithm to get
  *  the log2f of a 32-bit number
  *****************************************/
uint32_t svt_aom_log2f_32(uint32_t x) {
    uint32_t log = 0;
    int32_t  i;
    for (i = 4; i >= 0; --i) {
        const uint32_t shift = (1 << i);
        const uint32_t n     = x >> shift;
        if (n != 0) {
            x = n;
            log += shift;
        }
    }
    return log;
}
// concatenate two linked list, and return the pointer to the new concatenated list
EbLinkedListNode* svt_aom_concat_eb_linked_list(EbLinkedListNode* a, EbLinkedListNode* b) {
    if (a) {
        while (a->next) a = a->next;
        a->next = b;
        return a;
    } else
        return b;
}

// split a linked list
EbLinkedListNode* svt_aom_split_eb_linked_list(EbLinkedListNode* input, EbLinkedListNode** restLL,
                                               Bool (*predicate_func)(EbLinkedListNode*)) {
    EbLinkedListNode* ll_true_ptr = (EbLinkedListNode*)NULL; // list of nodes satifying predicate_func(node) == TRUE
    EbLinkedListNode* ll_rest_ptr = (EbLinkedListNode*)NULL; // list of nodes satifying predicate_func(node) != TRUE

    while (input) {
        EbLinkedListNode* next = input->next;
        input->next            = (EbLinkedListNode*)NULL;
        if (predicate_func(input))
            ll_true_ptr = svt_aom_concat_eb_linked_list(input, ll_true_ptr);
        else
            ll_rest_ptr = svt_aom_concat_eb_linked_list(input, ll_rest_ptr);
        input = next;
    }

    *restLL = ll_rest_ptr;
    return ll_true_ptr;
}

#if OPT_ENABLE_2L_INCOMP
static const MiniGopStats mini_gop_stats_array[] = {
    // hierarchical_levels    start_index    end_index    Length
    {5, 0, 31, 32}, // 0
    {4, 0, 15, 16}, // 1
    {3, 0, 7, 8}, // 2
    {2, 0, 3, 4}, // 3
    {1, 0, 1, 2}, // 4
    {1, 2, 3, 2}, // 5
    {2, 4, 7, 4}, // 6
    {1, 4, 5, 2}, // 7
    {1, 6, 7, 2}, // 8
    {3, 8, 15, 8}, // 9
    {2, 8, 11, 4}, // 10
    {1, 8, 9, 2}, // 11
    {1, 10, 11, 2}, // 12
    {2, 12, 15, 4}, // 13
    {1, 12, 13, 2}, // 14
    {1, 14, 15, 2}, // 15
    {4, 16, 31, 16}, // 16
    {3, 16, 23, 8}, // 17
    {2, 16, 19, 4}, // 18
    {1, 16, 17, 2}, // 19
    {1, 18, 19, 2}, // 20
    {2, 20, 23, 4}, // 21
    {1, 20, 21, 2}, // 22
    {1, 22, 23, 2}, // 23
    {3, 24, 31, 8}, // 24
    {2, 24, 27, 4}, // 25
    {1, 24, 25, 2}, // 26
    {1, 26, 27, 2}, // 27
    {2, 28, 31, 4}, // 28
    {1, 28, 29, 2}, // 29
    {1, 30, 31, 2} // 30
};
#else
static const MiniGopStats mini_gop_stats_array[] = {
    //    hierarchical_levels    start_index    end_index    Lenght    mini_gop_index
    {5, 0, 31, 32}, // 0
    {4, 0, 15, 16}, // 1
    {3, 0, 7, 8}, // 2
    {2, 0, 3, 4}, // 3
    {2, 4, 7, 4}, // 4
    {3, 8, 15, 8}, // 5
    {2, 8, 11, 4}, // 6
    {2, 12, 15, 4}, // 7
    {4, 16, 31, 16}, // 8
    {3, 16, 23, 8}, // 9
    {2, 16, 19, 4}, // 10
    {2, 20, 23, 4}, // 11
    {3, 24, 31, 8}, // 12
    {2, 24, 27, 4}, // 13
    {2, 28, 31, 4} // 14
};
#endif
/**************************************************************
* Get Mini GOP Statistics
**************************************************************/
const MiniGopStats* svt_aom_get_mini_gop_stats(const uint32_t mini_gop_index) {
    return &mini_gop_stats_array[mini_gop_index];
}

static uint32_t ns_quarter_off_mult[9 /*Up to 9 part*/][2 /*x+y*/][4 /*Up to 4 ns blocks per part*/] = {
    //9 means not used.

    //          |   x   |     |   y   |

    /*P=0*/ {{0, 9, 9, 9}, {0, 9, 9, 9}},
    /*P=1*/ {{0, 0, 9, 9}, {0, 2, 9, 9}},
    /*P=2*/ {{0, 2, 9, 9}, {0, 0, 9, 9}},
    /*P=3*/ {{0, 2, 0, 9}, {0, 0, 2, 9}},
    /*P=4*/ {{0, 0, 2, 9}, {0, 2, 2, 9}},
    /*P=5*/ {{0, 0, 2, 9}, {0, 2, 0, 9}},
    /*P=6*/ {{0, 2, 2, 9}, {0, 0, 2, 9}},
    /*P=7*/ {{0, 0, 0, 0}, {0, 1, 2, 3}},
    /*P=8*/ {{0, 1, 2, 3}, {0, 0, 0, 0}}};

static uint32_t ns_quarter_size_mult[9 /*Up to 9 part*/][2 /*h+v*/][4 /*Up to 4 ns blocks per part*/] = {
    //9 means not used.

    //          |   h   |     |   v   |

    /*P=0*/ {{4, 9, 9, 9}, {4, 9, 9, 9}},
    /*P=1*/ {{4, 4, 9, 9}, {2, 2, 9, 9}},
    /*P=2*/ {{2, 2, 9, 9}, {4, 4, 9, 9}},
    /*P=3*/ {{2, 2, 4, 9}, {2, 2, 2, 9}},
    /*P=4*/ {{4, 2, 2, 9}, {2, 2, 2, 9}},
    /*P=5*/ {{2, 2, 2, 9}, {2, 2, 4, 9}},
    /*P=6*/ {{2, 2, 2, 9}, {4, 2, 2, 9}},
    /*P=7*/ {{4, 4, 4, 4}, {1, 1, 1, 1}},
    /*P=8*/ {{1, 1, 1, 1}, {4, 4, 4, 4}}};

static BlockSize hvsize_to_bsize[/*H*/ 6][/*V*/ 6] = {
    {BLOCK_4X4, BLOCK_4X8, BLOCK_4X16, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_8X4, BLOCK_8X8, BLOCK_8X16, BLOCK_8X32, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_16X4, BLOCK_16X8, BLOCK_16X16, BLOCK_16X32, BLOCK_16X64, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_32X8, BLOCK_32X16, BLOCK_32X32, BLOCK_32X64, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X16, BLOCK_64X32, BLOCK_64X64, BLOCK_64X128},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, BLOCK_128X128}};

static uint32_t max_sb    = 64;
static uint32_t max_depth = 5;
static uint32_t max_part  = 9;
static uint32_t max_num_active_blocks;

GeomIndex svt_aom_geom_idx;
//TODO need to remove above globals for multi-channel support

BlockGeom svt_aom_blk_geom_mds
    [MAX_NUM_BLOCKS_ALLOC]; // to access geom info of a particular block; use this table if you have the block index in md scan
static INLINE TxSize av1_get_tx_size(BlockSize bsize, int32_t plane /*, const MacroBlockD *xd*/) {
    UNUSED(plane);
    //const MbModeInfo *mbmi = xd->mi[0];
    // if (xd->lossless[mbmi->segment_id]) return TX_4X4;
    if (plane == 0)
        return blocksize_to_txsize[bsize];
    // const MacroblockdPlane *pd = &xd->plane[plane];

    uint32_t subsampling_x = plane > 0 ? 1 : 0;
    uint32_t subsampling_y = plane > 0 ? 1 : 0;
    return av1_get_max_uv_txsize(/*mbmi->*/ bsize, subsampling_x, subsampling_y);
}
static void md_scan_all_blks(uint32_t* idx_mds, uint32_t sq_size, uint32_t x, uint32_t y, int32_t is_last_quadrant,
                             uint8_t quad_it, uint8_t min_nsq_bsize) {
    //the input block is the parent square block of size sq_size located at pos (x,y)

    assert(quad_it <= 3);
    uint32_t part_it, nsq_it, d1_it, sqi_mds;

    uint32_t halfsize  = sq_size / 2;
    uint32_t quartsize = sq_size / 4;

    uint32_t max_part_updated = sq_size == 128 ? MIN(max_part, 7)
        : sq_size == 8                         ? MIN(max_part, 3)
        :

        sq_size == 4 ? 1
                     : max_part;
    if (sq_size <= min_nsq_bsize)
        max_part_updated = 1;
    d1_it   = 0;
    sqi_mds = *idx_mds;

    for (part_it = 0; part_it < max_part_updated; part_it++) {
        uint32_t tot_num_ns_per_part = part_it < 1 ? 1 : part_it < 3 ? 2 : part_it < 7 ? 3 : 4;

        for (nsq_it = 0; nsq_it < tot_num_ns_per_part; nsq_it++) {
            svt_aom_blk_geom_mds[*idx_mds].depth = sq_size == max_sb / 1 ? 0
                : sq_size == max_sb / 2                                  ? 1
                : sq_size == max_sb / 4                                  ? 2
                : sq_size == max_sb / 8                                  ? 3
                : sq_size == max_sb / 16                                 ? 4
                                                                         : 5;

            svt_aom_blk_geom_mds[*idx_mds].sq_size          = sq_size;
            svt_aom_blk_geom_mds[*idx_mds].is_last_quadrant = is_last_quadrant;
            svt_aom_blk_geom_mds[*idx_mds].quadi            = quad_it;

            svt_aom_blk_geom_mds[*idx_mds].shape = (Part)part_it;
            svt_aom_blk_geom_mds[*idx_mds].org_x = x + quartsize * ns_quarter_off_mult[part_it][0][nsq_it];
            svt_aom_blk_geom_mds[*idx_mds].org_y = y + quartsize * ns_quarter_off_mult[part_it][1][nsq_it];

            svt_aom_blk_geom_mds[*idx_mds].d1i     = d1_it++;
            svt_aom_blk_geom_mds[*idx_mds].sqi_mds = sqi_mds;

            svt_aom_blk_geom_mds[*idx_mds].svt_aom_geom_idx = svt_aom_geom_idx;

            svt_aom_blk_geom_mds[*idx_mds].parent_depth_idx_mds = sqi_mds == 0
                ? 0
                : (sqi_mds + (3 - quad_it) * ns_depth_offset[svt_aom_geom_idx][svt_aom_blk_geom_mds[*idx_mds].depth]) -
                    parent_depth_offset[svt_aom_geom_idx][svt_aom_blk_geom_mds[*idx_mds].depth];
            svt_aom_blk_geom_mds[*idx_mds].d1_depth_offset =
                d1_depth_offset[svt_aom_geom_idx][svt_aom_blk_geom_mds[*idx_mds].depth];
            svt_aom_blk_geom_mds[*idx_mds].ns_depth_offset =
                ns_depth_offset[svt_aom_geom_idx][svt_aom_blk_geom_mds[*idx_mds].depth];
            svt_aom_blk_geom_mds[*idx_mds].totns   = tot_num_ns_per_part;
            svt_aom_blk_geom_mds[*idx_mds].nsi     = nsq_it;
            svt_aom_blk_geom_mds[*idx_mds].bwidth  = quartsize * ns_quarter_size_mult[part_it][0][nsq_it];
            svt_aom_blk_geom_mds[*idx_mds].bheight = quartsize * ns_quarter_size_mult[part_it][1][nsq_it];
            svt_aom_blk_geom_mds[*idx_mds].bsize =
                hvsize_to_bsize[svt_log2f(svt_aom_blk_geom_mds[*idx_mds].bwidth) - 2]
                               [svt_log2f(svt_aom_blk_geom_mds[*idx_mds].bheight) - 2];
            svt_aom_blk_geom_mds[*idx_mds].bwidth_uv  = MAX(4, svt_aom_blk_geom_mds[*idx_mds].bwidth >> 1);
            svt_aom_blk_geom_mds[*idx_mds].bheight_uv = MAX(4, svt_aom_blk_geom_mds[*idx_mds].bheight >> 1);
            svt_aom_blk_geom_mds[*idx_mds].has_uv     = 1;

            if (svt_aom_blk_geom_mds[*idx_mds].bwidth == 4 && svt_aom_blk_geom_mds[*idx_mds].bheight == 4)
                svt_aom_blk_geom_mds[*idx_mds].has_uv = is_last_quadrant ? 1 : 0;

            else if ((svt_aom_blk_geom_mds[*idx_mds].bwidth >> 1) < svt_aom_blk_geom_mds[*idx_mds].bwidth_uv ||
                     (svt_aom_blk_geom_mds[*idx_mds].bheight >> 1) < svt_aom_blk_geom_mds[*idx_mds].bheight_uv) {
                int32_t num_blk_same_uv = 1;
                if (svt_aom_blk_geom_mds[*idx_mds].bwidth >> 1 < 4)
                    num_blk_same_uv *= 2;
                if (svt_aom_blk_geom_mds[*idx_mds].bheight >> 1 < 4)
                    num_blk_same_uv *= 2;
                //if (svt_aom_blk_geom_mds[*idx_mds].nsi % 2 == 0)
                //if (svt_aom_blk_geom_mds[*idx_mds].nsi != (svt_aom_blk_geom_mds[*idx_mds].totns-1) )
                if (svt_aom_blk_geom_mds[*idx_mds].nsi != (num_blk_same_uv - 1) &&
                    svt_aom_blk_geom_mds[*idx_mds].nsi != (2 * num_blk_same_uv - 1))
                    svt_aom_blk_geom_mds[*idx_mds].has_uv = 0;
            }

            svt_aom_blk_geom_mds[*idx_mds].bsize_uv = get_plane_block_size(svt_aom_blk_geom_mds[*idx_mds].bsize, 1, 1);
            uint16_t txb_itr                        = 0;
            // tx_depth 1 geom settings
            uint8_t tx_depth                                   = 0;
            svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128
                ? 4
                : svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64 ||
                    svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128
                ? 2
                : 1;
            for (txb_itr = 0; txb_itr < svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth]; txb_itr++) {
                svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(svt_aom_blk_geom_mds[*idx_mds].bsize,
                                                                                  0);
                svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = av1_get_tx_size(
                    svt_aom_blk_geom_mds[*idx_mds].bsize, 1);
                if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128) {
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0 || txb_itr == 2)
                        ? svt_aom_blk_geom_mds[*idx_mds].org_x
                        : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0 || txb_itr == 1)
                        ? svt_aom_blk_geom_mds[*idx_mds].org_y
                        : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64) {
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0)
                        ? svt_aom_blk_geom_mds[*idx_mds].org_x
                        : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128) {
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0)
                        ? svt_aom_blk_geom_mds[*idx_mds].org_y
                        : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                } else {
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y;
                }
                /*if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X8)
                    SVT_LOG("");*/
                svt_aom_blk_geom_mds[*idx_mds].tx_width[tx_depth] =
                    tx_size_wide[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_height[tx_depth] =
                    tx_size_high[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_width_uv[tx_depth] =
                    tx_size_wide[svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_height_uv[tx_depth] =
                    tx_size_high[svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth]];
            }
            // tx_depth 1 geom settings
            tx_depth                                           = 1;
            svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128
                ? 4
                : svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64 ||
                    svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128
                ? 2
                : 1;

            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X16 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X8) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 4;
            }

            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X16 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X8 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X16) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 2;
            }
            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X16 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X8 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X4 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_4X16) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 2;
            }
            for (txb_itr = 0; txb_itr < svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth]; txb_itr++) {
                if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_32X32, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[4]                                 = {0, 32, 0, 32};
                    uint8_t offsety[4]                                 = {0, 0, 32, 32};
                    //   0  1
                    //   2  3
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_32X32, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 32};
                    uint8_t offsety[2]                                 = {0, 0};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_32X32, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 0};
                    uint8_t offsety[2]                                 = {0, 32};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[4]                                 = {0, 16, 0, 16};
                    uint8_t offsety[4]                                 = {0, 0, 16, 16};
                    //   0  1
                    //   2  3
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 16};
                    uint8_t offsety[2]                                 = {0, 0};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 0};
                    uint8_t offsety[2]                                 = {0, 16};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[4]                                 = {0, 8, 0, 8};
                    uint8_t offsety[4]                                 = {0, 0, 8, 8};
                    //   0  1
                    //   2  3
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X8) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 8};
                    uint8_t offsety[2]                                 = {0, 0};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[2]                                 = {0, 0};
                    uint8_t offsety[2]                                 = {0, 8};
                    //   0  1
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X8) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    uint8_t offsetx[4]                                 = {0, 4, 0, 4};
                    uint8_t offsety[4]                                 = {0, 0, 4, 4};
                    //   0  1
                    //   2  3
                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_32X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx[2] = {0, 32};
                    uint8_t offsety[2] = {0, 0};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X32, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx[2] = {0, 0};
                    uint8_t offsety[2] = {0, 32};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X8) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx[2] = {0, 16};
                    uint8_t offsety[2] = {0, 0};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[2] = {0, 0};
                    uint8_t offsety[2] = {0, 16};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X4) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx[2] = {0, 8};
                    uint8_t offsety[2] = {0, 0};

                    uint8_t tbx = offsetx[txb_itr];
                    uint8_t tby = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_4X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx[2] = {0, 0};
                    uint8_t offsety[2] = {0, 8};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else {
                    if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0 ||
                                                                                             txb_itr == 2)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_x
                            : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0 ||
                                                                                             txb_itr == 1)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_y
                            : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                    } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_x
                            : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_y;
                    } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_x;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_y
                            : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                    } else {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_x;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_y;
                    }
                }
                svt_aom_blk_geom_mds[*idx_mds].tx_width[tx_depth] =
                    tx_size_wide[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_height[tx_depth] =
                    tx_size_high[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_width_uv[tx_depth]  = svt_aom_blk_geom_mds[*idx_mds].tx_width_uv[0];
                svt_aom_blk_geom_mds[*idx_mds].tx_height_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].tx_height_uv[0];
            }
            // tx_depth 2 geom settings
            tx_depth = 2;

            svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128
                ? 4
                : svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64 ||
                    svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128
                ? 2
                : 1;

            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X16) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 16;
            }
            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X16 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X8 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X16) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 8;
            }
            if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X16 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X64 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X8 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X32 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X4 ||
                svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_4X16) {
                svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth] = 4;
            }

            for (txb_itr = 0; txb_itr < svt_aom_blk_geom_mds[*idx_mds].txb_count[tx_depth]; txb_itr++) {
                if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[16] = {0, 16, 32, 48, 0, 16, 32, 48, 0, 16, 32, 48, 0, 16, 32, 48};
                    uint8_t offsety_intra[16] = {0, 0, 0, 0, 16, 16, 16, 16, 32, 32, 32, 32, 48, 48, 48, 48};

                    uint8_t offsetx_inter[16] = {0, 16, 0, 16, 32, 48, 32, 48, 0, 16, 0, 16, 32, 48, 32, 48};
                    uint8_t offsety_inter[16] = {0, 0, 16, 16, 0, 0, 16, 16, 32, 32, 48, 48, 32, 32, 48, 48};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];

                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 16, 32, 48, 0, 16, 32, 48};
                    uint8_t offsety_intra[8] = {0, 0, 0, 0, 16, 16, 16, 16};

                    uint8_t offsetx_inter[8] = {0, 16, 0, 16, 32, 48, 32, 48};
                    uint8_t offsety_inter[8] = {0, 0, 16, 16, 0, 0, 16, 16};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 16, 0, 16, 0, 16, 0, 16};
                    uint8_t offsety_intra[8] = {0, 0, 16, 16, 32, 32, 48, 48};

                    uint8_t offsetx_inter[8] = {0, 16, 0, 16, 0, 16, 0, 16};
                    uint8_t offsety_inter[8] = {0, 0, 16, 16, 32, 32, 48, 48};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];

                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[16] = {0, 8, 16, 24, 0, 8, 16, 24, 0, 8, 16, 24, 0, 8, 16, 24};
                    uint8_t offsety_intra[16] = {0, 0, 0, 0, 8, 8, 8, 8, 16, 16, 16, 16, 24, 24, 24, 24};

                    uint8_t offsetx_inter[16] = {0, 8, 0, 8, 16, 24, 16, 24, 0, 8, 0, 8, 16, 24, 16, 24};
                    uint8_t offsety_inter[16] = {0, 0, 8, 8, 0, 0, 8, 8, 16, 16, 24, 24, 16, 16, 24, 24};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 8, 16, 24, 0, 8, 16, 24};
                    uint8_t offsety_intra[8] = {0, 0, 0, 0, 8, 8, 8, 8};

                    uint8_t offsetx_inter[8] = {0, 8, 0, 8, 16, 24, 16, 24};
                    uint8_t offsety_inter[8] = {0, 0, 8, 8, 0, 0, 8, 8};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 8, 0, 8, 0, 8, 0, 8};
                    uint8_t offsety_intra[8] = {0, 0, 8, 8, 16, 16, 24, 24};

                    uint8_t offsetx_inter[8] = {0, 8, 0, 8, 0, 8, 0, 8};
                    uint8_t offsety_inter[8] = {0, 0, 8, 8, 16, 16, 24, 24};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X8) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 4, 8, 12, 0, 4, 8, 12};
                    uint8_t offsety_intra[8] = {0, 0, 0, 0, 4, 4, 4, 4};

                    uint8_t offsetx_inter[8] = {0, 4, 0, 4, 8, 12, 8, 12};
                    uint8_t offsety_inter[8] = {0, 0, 4, 4, 0, 0, 4, 4};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[8] = {0, 4, 0, 4, 0, 4, 0, 4};
                    uint8_t offsety_intra[8] = {0, 0, 4, 4, 8, 8, 12, 12};

                    uint8_t offsetx_inter[8] = {0, 4, 0, 4, 0, 4, 0, 4};
                    uint8_t offsety_inter[8] = {0, 0, 4, 4, 8, 8, 12, 12};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];

                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];

                    uint8_t offsetx_intra[16] = {0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12};
                    uint8_t offsety_intra[16] = {0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12};

                    uint8_t offsetx_inter[16] = {0, 4, 0, 4, 8, 12, 8, 12, 0, 4, 0, 4, 8, 12, 8, 12};
                    uint8_t offsety_inter[16] = {0, 0, 4, 4, 0, 0, 4, 4, 8, 8, 12, 12, 8, 8, 12, 12};

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_intra[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_intra[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_x + offsetx_inter[txb_itr];
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].org_y + offsety_inter[txb_itr];
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 16, 32, 48};
                    uint8_t offsety[4] = {0, 0, 0, 0};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X64) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_16X16, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 0, 0, 0};
                    uint8_t offsety[4] = {0, 16, 32, 48};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_32X8) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 8, 16, 24};
                    uint8_t offsety[4] = {0, 0, 0, 0};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_8X32) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_8X8, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 0, 0, 0};
                    uint8_t offsety[4] = {0, 8, 16, 24};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_16X4) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 4, 8, 12};
                    uint8_t offsety[4] = {0, 0, 0, 0};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_4X16) {
                    svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]    = av1_get_tx_size(BLOCK_4X4, 0);
                    svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                    //   0  1 2 3
                    uint8_t offsetx[4] = {0, 0, 0, 0};
                    uint8_t offsety[4] = {0, 4, 8, 12};
                    uint8_t tbx        = offsetx[txb_itr];
                    uint8_t tby        = offsety[txb_itr];

                    svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_x + tbx;
                    svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].org_y + tby;
                } else {
                    if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X128) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0 ||
                                                                                             txb_itr == 2)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_x
                            : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0 ||
                                                                                             txb_itr == 1)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_y
                            : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                    } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_128X64) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] = (txb_itr == 0)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_x
                            : svt_aom_blk_geom_mds[*idx_mds].org_x + 64;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_y;
                    } else if (svt_aom_blk_geom_mds[*idx_mds].bsize == BLOCK_64X128) {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_x;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] = (txb_itr == 0)
                            ? svt_aom_blk_geom_mds[*idx_mds].org_y
                            : svt_aom_blk_geom_mds[*idx_mds].org_y + 64;
                    } else {
                        svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth] = av1_get_tx_size(
                            svt_aom_blk_geom_mds[*idx_mds].bsize, 0);
                        svt_aom_blk_geom_mds[*idx_mds].txsize_uv[tx_depth] =
                            svt_aom_blk_geom_mds[*idx_mds].txsize_uv[0];
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_x[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_x[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_x;
                        svt_aom_blk_geom_mds[*idx_mds].tx_org_y[0][tx_depth][txb_itr] =
                            svt_aom_blk_geom_mds[*idx_mds].tx_org_y[1][tx_depth][txb_itr] =
                                svt_aom_blk_geom_mds[*idx_mds].org_y;
                    }
                }
                svt_aom_blk_geom_mds[*idx_mds].tx_width[tx_depth] =
                    tx_size_wide[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_height[tx_depth] =
                    tx_size_high[svt_aom_blk_geom_mds[*idx_mds].txsize[tx_depth]];
                svt_aom_blk_geom_mds[*idx_mds].tx_width_uv[tx_depth]  = svt_aom_blk_geom_mds[*idx_mds].tx_width_uv[0];
                svt_aom_blk_geom_mds[*idx_mds].tx_height_uv[tx_depth] = svt_aom_blk_geom_mds[*idx_mds].tx_height_uv[0];
            }
            svt_aom_blk_geom_mds[*idx_mds].blkidx_mds = (*idx_mds);
            (*idx_mds)                                = (*idx_mds) + 1;
        }
    }

    uint32_t min_size = max_sb >> (max_depth - 1);
    if (halfsize >= min_size) {
        md_scan_all_blks(idx_mds, halfsize, x, y, 0, 0, min_nsq_bsize);
        md_scan_all_blks(idx_mds, halfsize, x + halfsize, y, 0, 1, min_nsq_bsize);
        md_scan_all_blks(idx_mds, halfsize, x, y + halfsize, 0, 2, min_nsq_bsize);
        md_scan_all_blks(idx_mds, halfsize, x + halfsize, y + halfsize, 1, 3, min_nsq_bsize);
    }
}
static uint32_t count_total_num_of_active_blks(uint8_t min_nsq_bsize) {
    uint32_t depth_it, sq_it_y, sq_it_x, part_it, nsq_it;

    uint32_t depth_scan_idx = 0;

    for (depth_it = 0; depth_it < max_depth; depth_it++) {
        uint32_t tot_num_sq = 1 << depth_it;
        uint32_t sq_size    = depth_it == 0 ? max_sb
               : depth_it == 1              ? max_sb / 2
               : depth_it == 2              ? max_sb / 4
               : depth_it == 3              ? max_sb / 8
               : depth_it == 4              ? max_sb / 16
                                            : max_sb / 32;

        uint32_t max_part_updated = sq_size == 128 ? MIN(max_part, 7)
            : sq_size == 8                         ? MIN(max_part, 3)
            : sq_size == 4                         ? 1
                                                   : max_part;
        if (sq_size <= min_nsq_bsize)
            max_part_updated = 1;
        for (sq_it_y = 0; sq_it_y < tot_num_sq; sq_it_y++) {
            for (sq_it_x = 0; sq_it_x < tot_num_sq; sq_it_x++) {
                for (part_it = 0; part_it < max_part_updated; part_it++) {
                    uint32_t tot_num_ns_per_part = part_it < 1 ? 1 : part_it < 3 ? 2 : part_it < 7 ? 3 : 4;

                    for (nsq_it = 0; nsq_it < tot_num_ns_per_part; nsq_it++) depth_scan_idx++;
                }
            }
        }
    }

    return depth_scan_idx;
}
static void log_redundancy_similarity(uint32_t max_block_count) {
    uint32_t blk_it, s_it;

    for (blk_it = 0; blk_it < max_block_count; blk_it++) {
        BlockGeom* cur_geom             = &svt_aom_blk_geom_mds[blk_it];
        cur_geom->redund                = 0;
        cur_geom->redund_list.list_size = 0;

        for (s_it = 0; s_it < max_block_count; s_it++) {
            BlockGeom* search_geom = &svt_aom_blk_geom_mds[s_it];

            if (cur_geom->bsize == search_geom->bsize && cur_geom->org_x == search_geom->org_x &&
                cur_geom->org_y == search_geom->org_y && s_it != blk_it) {
                if (cur_geom->nsi == 0 && search_geom->nsi == 0 && cur_geom->redund_list.list_size < 3) {
                    cur_geom->redund                                                     = 1;
                    cur_geom->redund_list.blk_mds_table[cur_geom->redund_list.list_size] = search_geom->blkidx_mds;
                    cur_geom->redund_list.list_size++;
                }
            }
        }
    }
}

/*
  Build Block Geometry
*/
void svt_aom_build_blk_geom(GeomIndex geom) {
    uint32_t max_block_count;
    svt_aom_geom_idx = geom;
    uint32_t min_nsq_bsize;
    if (geom == GEOM_0) {
        max_sb          = 64;
        max_depth       = 4;
        max_part        = 1;
        max_block_count = 85;
        min_nsq_bsize   = 16;
    } else if (geom == GEOM_1) {
        max_sb          = 64;
        max_depth       = 4;
        max_part        = 3;
        max_block_count = 105;
        min_nsq_bsize   = 16;
    } else if (geom == GEOM_2) {
        max_sb          = 64;
        max_depth       = 4;
        max_part        = 3;
        max_block_count = 169;
        min_nsq_bsize   = 8;
    } else if (geom == GEOM_3) {
        max_sb          = 64;
        max_depth       = 4;
        max_part        = 3;
        max_block_count = 425;
        min_nsq_bsize   = 0;
    } else if (geom == GEOM_4) {
        max_sb          = 64;
        max_depth       = 5;
        max_part        = 3;
        max_block_count = 681;
        min_nsq_bsize   = 0;
    } else if (geom == GEOM_5) {
        max_sb          = 64;
        max_depth       = 5;
        max_part        = 9;
        max_block_count = 1101;
        min_nsq_bsize   = 0;
    } else {
        max_sb          = 128;
        max_depth       = 6;
        max_part        = 9;
        max_block_count = 4421;
        min_nsq_bsize   = 0;
    }
    //(0)compute total number of blocks using the information provided
    max_num_active_blocks = count_total_num_of_active_blks(min_nsq_bsize);
    if (max_num_active_blocks != max_block_count)
        SVT_LOG(" \n\n Error %i blocks\n\n ", max_num_active_blocks);
    //(2) Construct md scan blk_geom_mds:  use info from dps
    uint32_t idx_mds = 0;
    md_scan_all_blks(&idx_mds, max_sb, 0, 0, 0, 0, min_nsq_bsize);
    log_redundancy_similarity(max_block_count);
}
uint32_t get_mds_idx(uint32_t orgx, uint32_t orgy, uint32_t size, uint32_t use_128x128) {
    (void)use_128x128;
    uint32_t max_block_count = max_num_active_blocks;
    uint32_t mds             = 0;

    for (uint32_t blk_it = 0; blk_it < max_block_count; blk_it++) {
        BlockGeom* cur_geom = &svt_aom_blk_geom_mds[blk_it];

        if ((uint32_t)cur_geom->sq_size == size && cur_geom->org_x == orgx && cur_geom->org_y == orgy &&
            cur_geom->shape == PART_N) {
            mds = cur_geom->blkidx_mds;
            break;
        }
    }
    return mds;
}

#if FIXED_POINT_ASSERT_TEST
void svt_fixed_point_test_breakpoint(char* file, unsigned line) {
    printf("ERROR: Fixed Point Test Assert:  %s:%u", file, line);
}
#endif
