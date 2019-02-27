/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbUtility_h
#define EbUtility_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
    /****************************
     * UTILITY FUNCTIONS
     ****************************/

    void build_blk_geom();
    typedef struct BlockGeom
    {
        uint8_t    depth;                       // depth of the block
        PART       shape;                       // P_N..P_V4 . P_S is not used.
        uint8_t    origin_x;                    // orgin x from topleft of sb
        uint8_t    origin_y;                    // orgin x from topleft of sb
                   
        uint8_t    d1i;                         // index of the block in d1 dimension 0..24  (0 is parent square, 1 top half of H , ...., 24:last quarter of V4)
        uint16_t   sqi_mds;                     // index of the parent square in md  scan.
        uint8_t    totns;                       // max number of ns blocks within one partition 1..4 (N:1,H:2,V:2,HA:3,HB:3,VA:3,VB:3,H4:4,V4:4)
        uint8_t    nsi;                         // non square index within a partition  0..totns-1
                   
                   
        uint8_t    bwidth;                      // block width
        uint8_t    bheight;                     // block height
        uint8_t    bwidth_uv;                   // block width for Chroma 4:2:0
        uint8_t    bheight_uv;                  // block height for Chroma 4:2:0
        uint8_t    bwidth_log2;                 // block width log2
        uint8_t    bheight_log2;                // block height log2
        block_size bsize;                       // bloc size
        block_size bsize_uv;                    // bloc size for Chroma 4:2:0
        uint16_t   txb_count;                   //4-2-1
        TxSize     txsize[MAX_TXB_COUNT];
        TxSize     txsize_uv[MAX_TXB_COUNT];
        uint16_t   tx_org_x[MAX_TXB_COUNT];     //orgin is SB
        uint16_t   tx_org_y[MAX_TXB_COUNT];     //origin is SB
        uint16_t   tx_boff_x[MAX_TXB_COUNT];    //block offset , origin is block
        uint16_t   tx_boff_y[MAX_TXB_COUNT];    //block offset , origin is block
        uint8_t    tx_width[MAX_TXB_COUNT];     //tx_size_wide
        uint8_t    tx_height[MAX_TXB_COUNT];    //tx_size_wide
        uint8_t    tx_width_uv[MAX_TXB_COUNT];  //tx_size_wide
        uint8_t    tx_height_uv[MAX_TXB_COUNT]; //tx_size_wide
                   
                   
        uint16_t   blkidx_mds;                  // block index in md scan
        uint16_t   blkidx_dps;                  // block index in depth scan
        int32_t    has_uv;
        int32_t    sq_size;
        int32_t    is_last_quadrant;            // only for square bloks, is this the fourth quadrant block?


    } BlockGeom;

    static const block_size ss_size_lookup[BlockSizeS_ALL][2][2] = {
        //  ss_x == 0    ss_x == 0        ss_x == 1      ss_x == 1
        //  ss_y == 0    ss_y == 1        ss_y == 0      ss_y == 1
        { { BLOCK_4X4, BLOCK_4X4 }, { BLOCK_4X4, BLOCK_4X4 } },
        { { BLOCK_4X8, BLOCK_4X4 }, { BLOCK_4X4, BLOCK_4X4 } },
        { { BLOCK_8X4, BLOCK_4X4 }, { BLOCK_4X4, BLOCK_4X4 } },
        { { BLOCK_8X8, BLOCK_8X4 }, { BLOCK_4X8, BLOCK_4X4 } },
        { { BLOCK_8X16, BLOCK_8X8 }, { BLOCK_4X16, BLOCK_4X8 } },
        { { BLOCK_16X8, BLOCK_16X4 }, { BLOCK_8X8, BLOCK_8X4 } },
        { { BLOCK_16X16, BLOCK_16X8 }, { BLOCK_8X16, BLOCK_8X8 } },
        { { BLOCK_16X32, BLOCK_16X16 }, { BLOCK_8X32, BLOCK_8X16 } },
        { { BLOCK_32X16, BLOCK_32X8 }, { BLOCK_16X16, BLOCK_16X8 } },
        { { BLOCK_32X32, BLOCK_32X16 }, { BLOCK_16X32, BLOCK_16X16 } },
        { { BLOCK_32X64, BLOCK_32X32 }, { BLOCK_16X64, BLOCK_16X32 } },
        { { BLOCK_64X32, BLOCK_64X16 }, { BLOCK_32X32, BLOCK_32X16 } },
        { { BLOCK_64X64, BLOCK_64X32 }, { BLOCK_32X64, BLOCK_32X32 } },
        { { BLOCK_64X128, BLOCK_64X64 }, { BLOCK_INVALID, BLOCK_32X64 } },
        { { BLOCK_128X64, BLOCK_INVALID }, { BLOCK_64X64, BLOCK_64X32 } },
        { { BLOCK_128X128, BLOCK_128X64 }, { BLOCK_64X128, BLOCK_64X64 } },
        { { BLOCK_4X16, BLOCK_4X8 }, { BLOCK_4X16, BLOCK_4X8 } },
        { { BLOCK_16X4, BLOCK_16X4 }, { BLOCK_8X4, BLOCK_8X4 } },
        { { BLOCK_8X32, BLOCK_8X16 }, { BLOCK_INVALID, BLOCK_4X16 } },
        { { BLOCK_32X8, BLOCK_INVALID }, { BLOCK_16X8, BLOCK_16X4 } },
        { { BLOCK_16X64, BLOCK_16X32 }, { BLOCK_INVALID, BLOCK_8X32 } },
        { { BLOCK_64X16, BLOCK_INVALID }, { BLOCK_32X16, BLOCK_32X8 } }
    };
    static INLINE block_size get_plane_block_size(block_size bsize,
        int32_t subsampling_x,
        int32_t subsampling_y) {
        if (bsize == BLOCK_INVALID) return BLOCK_INVALID;
        return ss_size_lookup[bsize][subsampling_x][subsampling_y];
    }

#define NOT_USED_VALUE  0
    static const uint32_t parent_depth_offset[2][6] =
    { /*64x64*/ { NOT_USED_VALUE, 832, 208, 52, 8 ,NOT_USED_VALUE},
        /*128x128*/{ NOT_USED_VALUE,3320 , 832 ,208,52,8 }
    };
    static const uint32_t ns_depth_offset[2][6] =
    { /*64x64*/  { 1101 , 269  ,61,9,1 ,NOT_USED_VALUE},
        /*128x128*/{ 4421 , 1101 ,269 ,61,9,1 }
    };
    static  const uint32_t d1_depth_offset[2][6] =
    { /*64x64*/ { 25,25,25, 5 ,1 ,NOT_USED_VALUE} ,
        /*128x128*/{ 17 , 25,25,25, 5 ,1 }
    };


    const BlockGeom * get_blk_geom_mds(uint32_t bidx_mds);




    // CU Stats Helper Functions
    typedef struct CodedUnitStats_s
    {
        uint8_t   depth;
        uint8_t   size;
        uint8_t   size_log2;
        uint16_t  origin_x;
        uint16_t  origin_y;
        uint8_t   cuNumInDepth;
        uint8_t   parent32x32Index;

    } CodedUnitStats_t;

    // PU Stats Helper Functions
    typedef struct PredictionUnitStats_t
    {
        uint8_t  width;
        uint8_t  height;
        uint8_t  offsetX;
        uint8_t  offsetY;

    } PredictionUnitStats_t;

    // TU Stats Helper Functions
    typedef struct TransformUnitStats_s
    {
        uint8_t  depth;
        uint8_t  offsetX;
        uint8_t  offsetY;

    } TransformUnitStats_t;

    extern uint64_t Log2fHighPrecision(uint64_t x, uint8_t precision);

    extern const CodedUnitStats_t* GetCodedUnitStats(const uint32_t cuIdx);
    extern const TransformUnitStats_t* GetTransformUnitStats(const uint32_t tuIdx);

#define PU_ORIGIN_ADJUST(cuOrigin, cu_size, offset) ((((cu_size) * (offset)) >> 2) + (cuOrigin))
#define PU_SIZE_ADJUST(cu_size, puSize) (((cu_size) * (puSize)) >> 2)

#define TU_ORIGIN_ADJUST(cuOrigin, cu_size, offset) ((((cu_size) * (offset)) >> 2) + (cuOrigin))
#define TU_SIZE_ADJUST(cu_size, tuDepth) ((cu_size) >> (tuDepth))

    extern EbErrorType ZOrderIncrement(uint32_t *xLoc, uint32_t *yLoc);
    extern void ZOrderIncrementWithLevel(
        uint32_t *xLoc,
        uint32_t *yLoc,
        uint32_t *level,
        uint32_t *index);

    extern uint32_t Log2f(uint32_t x);
    extern uint64_t Log2f64(uint64_t x);
    extern uint32_t EndianSwap(uint32_t ui);

    /****************************
     * MACROS
     ****************************/

#ifdef _MSC_VER
#define MULTI_LINE_MACRO_BEGIN do {
#define MULTI_LINE_MACRO_END \
    __pragma(warning(push)) \
    __pragma(warning(disable:4127)) \
    } while(0) \
    __pragma(warning(pop))
#else
#define MULTI_LINE_MACRO_BEGIN do {
#define MULTI_LINE_MACRO_END } while(0)
#endif

     //**************************************************
     // MACROS
     //**************************************************
#define MAX(x, y)                       ((x)>(y)?(x):(y))
#define MIN(x, y)                       ((x)<(y)?(x):(y))
#define MEDIAN(a,b,c)                   ((a)>(b)?(a)>©?(b)>©?(b):©:(a):(b)>©?(a)>©?(a):©:(b))
#define CLIP3(MinVal, MaxVal, a)        (((a)<(MinVal)) ? (MinVal) : (((a)>(MaxVal)) ? (MaxVal) :(a)))
#define CLIP3EQ(MinVal, MaxVal, a)        (((a)<=(MinVal)) ? (MinVal) : (((a)>=(MaxVal)) ? (MaxVal) :(a)))
#define BITDEPTH_MIDRANGE_VALUE(precision)  (1 << ((precision) - 1))
#define SWAP(a, b)                      MULTI_LINE_MACRO_BEGIN (a) ^= (b); (b) ^= (a); (a) ^= (b); MULTI_LINE_MACRO_END
#define ABS(a)                          (((a) < 0) ? (-(a)) : (a))
#define EB_ABS_DIFF(a,b)                ((a) > (b) ? ((a) - (b)) : ((b) - (a)))
#define EB_DIFF_SQR(a,b)                (((a) - (b)) * ((a) - (b)))
#define SQR(x)                          ((x)*(x))
#define POW2(x)                         (1 << (x))
#define SIGN(a,b)                       (((a - b) < 0) ? (-1) : ((a - b) > 0) ? (1) : 0)
#define ROUND(a)                        (a >= 0) ? ( a + 1/2) : (a - 1/2);
#define UNSIGNED_DEC(x)                 MULTI_LINE_MACRO_BEGIN (x) = (((x) > 0) ? ((x) - 1) : 0); MULTI_LINE_MACRO_END
#define CIRCULAR_ADD(x,max)             (((x) >= (max)) ? ((x) - (max)) : ((x) < 0) ? ((max) + (x)) : (x))
#define CIRCULAR_ADD_UNSIGNED(x,max)    (((x) >= (max)) ? ((x) - (max)) : (x))
#define CEILING(x,base)                 ((((x) + (base) - 1) / (base)) * (base))
#define POW2_CHECK(x)                   ((x) == ((x) & (-((int32_t)(x)))))
#define ROUND_UP_MUL_8(x)               ((x) + ((8 - ((x) & 0x7)) & 0x7))
#define ROUND_UP_MULT(x,mult)           ((x) + (((mult) - ((x) & ((mult)-1))) & ((mult)-1)))

// rounds down to the next power of two
#define FLOOR_POW2(x)       \
    MULTI_LINE_MACRO_BEGIN  \
        (x) |= ((x) >> 1);  \
        (x) |= ((x) >> 2);  \
        (x) |= ((x) >> 4);  \
        (x) |= ((x) >> 8);  \
        (x) |= ((x) >>16);  \
        (x) -= ((x) >> 1);  \
    MULTI_LINE_MACRO_END

// rounds up to the next power of two
#define CEIL_POW2(x)        \
    MULTI_LINE_MACRO_BEGIN  \
        (x) -= 1;           \
        (x) |= ((x) >> 1);  \
        (x) |= ((x) >> 2);  \
        (x) |= ((x) >> 4);  \
        (x) |= ((x) >> 8);  \
        (x) |= ((x) >>16);  \
        (x) += 1;           \
    MULTI_LINE_MACRO_END

// Calculates the Log2 floor of the integer 'x'
//   Intended to only be used for macro definitions
#define LOG2F Log2f_SSE2

#define LOG2F_8(x) (              \
    ((x) < 0x0002u) ? 0u :      \
    ((x) < 0x0004u) ? 1u :      \
    ((x) < 0x0008u) ? 2u :      \
    ((x) < 0x0010u) ? 3u :      \
    ((x) < 0x0020u) ? 4u :      \
    ((x) < 0x0040u) ? 5u :6u )

#define TWO_D_INDEX(x, y, stride)   \
    (((y) * (stride)) + (x))

// MAX_CU_COUNT is used to find the total number of partitions for the max partition depth and for
// each parent partition up to the root partition level (i.e. SB level).

// MAX_CU_COUNT is given by SUM from k=1 to n (4^(k-1)), reduces by using the following finite sum
// SUM from k=1 to n (q^(k-1)) = (q^n - 1)/(q-1) => (4^n - 1) / 3
#define MAX_CU_COUNT(max_depth_count) ((((1 << (max_depth_count)) * (1 << (max_depth_count))) - 1)/3)

//**************************************************
// CONSTANTS
//**************************************************
#define MIN_UNSIGNED_VALUE      0
#define MAX_UNSIGNED_VALUE     ~0u
#define MIN_SIGNED_VALUE       ~0 - ((signed) (~0u >> 1))
#define MAX_SIGNED_VALUE       ((signed) (~0u >> 1))


// Helper functions for EbLinkedListNode.

// concatenate two linked list, and return the pointer to the new concatenated list
    EbLinkedListNode* concatEbLinkedList(EbLinkedListNode* a, EbLinkedListNode* b);

    // split a linked list into two. return the pointer to a linked list whose nodes meets the condition
    // predicateFunc(node) == TRUE, the rest of the nodes will be collected into another linked list to which (*restLL) is
    // set. Does not gaurantee the original order of the nodes.

    EbLinkedListNode* splitEbLinkedList(EbLinkedListNode* input, EbLinkedListNode** restLL, EbBool(*predicateFunc)(EbLinkedListNode*));

#define MINI_GOP_MAX_COUNT            15
#define MINI_GOP_WINDOW_MAX_COUNT     8    // widow subdivision: 8 x 3L

#define MIN_HIERARCHICAL_LEVEL         2
    static const uint32_t MiniGopOffset[4] = { 1, 3, 7, 31 };

    typedef struct MiniGopStats_s
    {
        uint32_t  hierarchical_levels;
        uint32_t  startIndex;
        uint32_t  endIndex;
        uint32_t  lenght;

    } MiniGopStats_t;
    extern const MiniGopStats_t* GetMiniGopStats(const uint32_t miniGopIndex);
    typedef enum MINI_GOP_INDEX {
        L6_INDEX = 0,
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
    } MINI_GOP_INDEX;

#ifdef __cplusplus
}
#endif

#endif // EbUtility_h
