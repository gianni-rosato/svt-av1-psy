// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimation_h
#define EbMotionEstimation_h

#include "EbDefinitions.h"
#include "EbCodingUnit.h"

#include "EbMotionEstimationProcess.h"
#include "EbMotionEstimationContext.h"
#include "EbPictureBufferDesc.h"
#include "EbReferenceObject.h"
#include "EbPictureDecisionProcess.h"
#include "EbUtility.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern EbErrorType motion_estimate_sb(
        PictureParentControlSet   *pcs_ptr,
        uint32_t                       sb_index,
        uint32_t                       sb_origin_x,
        uint32_t                       sb_origin_y,
        MeContext                 *context_ptr,
        EbPictureBufferDesc       *input_ptr);

    extern void decimation_2d(
        uint8_t                   *input_samples,
        uint32_t                   input_stride,
        uint32_t                   input_area_width,
        uint32_t                   input_area_height,
        uint8_t                   *decim_samples,
        uint32_t                   decim_stride,
        uint32_t                   decim_step);

    extern void downsample_2d(
        uint8_t                   *input_samples,
        uint32_t                   input_stride,
        uint32_t                   input_area_width,
        uint32_t                   input_area_height,
        uint8_t                   *decim_samples,
        uint32_t                   decim_stride,
        uint32_t                   decim_step);

    extern EbErrorType open_loop_intra_search_sb(
        PictureParentControlSet   *pcs_ptr,
        uint32_t                       sb_index,
        MotionEstimationContext_t   *context_ptr,
        EbPictureBufferDesc       *input_ptr);
    

#if TPL_LA
    extern EbErrorType av1_open_loop_intra_search(
        PictureParentControlSet   *picture_control_set_ptr,
        MotionEstimationContext_t *context_ptr,
        EbPictureBufferDesc       *input_ptr);
    extern EbErrorType open_loop_intra_search_mb(
        PictureParentControlSet   *picture_control_set_ptr,
        uint32_t                   sb_index,
        EbPictureBufferDesc       *input_ptr);
#endif
#define a_b_c  0
#define a_c_b  1
#define b_a_c  2
#define b_c_a  3
#define c_a_b  4
#define c_b_a  5

#define TOP_LEFT_POSITION       0
#define TOP_POSITION            1
#define TOP_RIGHT_POSITION      2
#define RIGHT_POSITION          3
#define BOTTOM_RIGHT_POSITION   4
#define BOTTOM_POSITION         5
#define BOTTOM_LEFT_POSITION    6
#define LEFT_POSITION           7

    // The interpolation is performed using a set of three 4 tap filters
#define IFShiftAvcStyle         1
#define F0 0
#define F1 1
#define F2 2
#define MAX_SSE_VALUE 128 * 128 * 255 * 255
#define  MAX_SAD_VALUE 128*128*255

// Interpolation Filters
    static const int32_t me_if_coeff[3][4] = {
        { -4, 54, 16, -2 }, // F0
        { -4, 36, 36, -4 }, // F1
        { -2, 16, 54, -4 }, // F2
    };

    static const uint32_t tab16x16[16] = {
        0, 1, 4, 5,
        2, 3, 6, 7,
        8, 9, 12, 13,
        10, 11, 14, 15
    };

    static const uint32_t tab32x32[16] = {
        0, 1,
        2, 3
    };

    static const uint32_t tab64x32[16] = {
        0, 1,
    };

    static const uint32_t tab32x64[16] = {
        0, 1,
    };

    static const uint32_t tab8x8[64] = {
        0, 1, 4, 5, 16, 17, 20, 21,
        2, 3, 6, 7, 18, 19, 22, 23,
        8, 9, 12, 13, 24, 25, 28, 29,
        10, 11, 14, 15, 26, 27, 30, 31,
        32, 33, 36, 37, 48, 49, 52, 53,
        34, 35, 38, 39, 50, 51, 54, 55,
        40, 41, 44, 45, 56, 57, 60, 61,
        42, 43, 46, 47, 58, 59, 62, 63
    };

    static const uint32_t tab32x16[8] = {
        0, 2,
        1, 3,
        4, 6,
        5, 7,
    };

    static const uint32_t tab16x32[8] = {
        0, 1, 2, 3,
        4, 5, 6, 7,
    };

    static const uint32_t tab32x8[16] = {
        0, 4,
        1, 5,
        2, 6,
        3, 7,
        8, 12,
        9, 13,
        10, 14,
        11, 15
    };

    static const uint32_t tab8x32[16] = {
        0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15
    };

    static const uint32_t tab16x8[32] = {
        0, 2, 8, 10,
        1, 3, 9, 11,
        4, 6, 12, 14,
        5, 7, 13, 15,
        16, 18, 24, 26,
        17, 19, 25, 27,
        20, 22, 28, 30,
        21, 23, 29, 31,
    };

    static const uint32_t tab8x16[32] = {
        0, 1, 2, 3, 8, 9, 10, 11,
        4, 5, 6, 7, 12, 13, 14, 15,
        16, 17, 18, 19, 24, 25, 26, 27,
        20, 21, 22, 23, 28, 29, 30, 31,
    };
#if NSQ_ME_CONTEXT_CLEAN_UP
    static const uint32_t partition_width[SQUARE_PU_COUNT] = {
#else
    static const uint32_t partition_width[MAX_ME_PU_COUNT] = {
#endif
        64,                                                                          // (1)
        32, 32, 32, 32,                                                              // (4)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
#if !NSQ_ME_CONTEXT_CLEAN_UP

        //H  Partitions
        64, 64,                                                                      // (2)
        32, 32, 32, 32, 32, 32, 32, 32,                                              // (8)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)

        //V  Partitions
        32, 32,                                                                      // (2)
        16, 16, 16, 16, 16, 16, 16, 16,                                              // (8)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)

        // H4 Partitions
        32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,              // (16)

        // V4 Partitions
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8                               // (16)

        ,64,64,64,64,
         16,16,16,16
#endif
    };
#if NSQ_ME_CONTEXT_CLEAN_UP
    static const uint32_t partition_height[SQUARE_PU_COUNT] = {
#else
    static const uint32_t partition_height[MAX_ME_PU_COUNT] = {
#endif
        64,                                                                          // (1)
        32, 32, 32, 32,                                                              // (4)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
#if !NSQ_ME_CONTEXT_CLEAN_UP

        // H Partitions
        32, 32,                                                                      // (2)
        16, 16, 16, 16, 16, 16, 16, 16,                                              // (8)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)

        //V  Partitions
        64, 64,                                                                      // (2)
        32, 32, 32, 32, 32, 32, 32, 32,                                              // (8)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,              // (16)

        // H4 Partitions
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,                              // (16)

        // V4 Partitions
        32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32               // (16)

        ,16,16,16,16,
         64,64,64,64
#endif
    };
#if NSQ_ME_CONTEXT_CLEAN_UP
    static const uint32_t pu_search_index_map[SQUARE_PU_COUNT][2] = {
#else
    static const uint32_t pu_search_index_map[MAX_ME_PU_COUNT][2] = {
#endif
        { 0, 0 },
        { 0, 0 }, { 32, 0 }, { 0, 32 }, { 32, 32 },
        { 0, 0 }, { 16, 0 }, { 32, 0 }, { 48, 0 },
        { 0, 16 }, { 16, 16 }, { 32, 16 }, { 48, 16 },
        { 0, 32 }, { 16, 32 }, { 32, 32 }, { 48, 32 },
        { 0, 48 }, { 16, 48 }, { 32, 48 }, { 48, 48 },
        { 0, 0 }, { 8, 0 }, { 16, 0 }, { 24, 0 }, { 32, 0 }, { 40, 0 }, { 48, 0 }, { 56, 0 },
        { 0, 8 }, { 8, 8 }, { 16, 8 }, { 24, 8 }, { 32, 8 }, { 40, 8 }, { 48, 8 }, { 56, 8 },
        { 0, 16 }, { 8, 16 }, { 16, 16 }, { 24, 16 }, { 32, 16 }, { 40, 16 }, { 48, 16 }, { 56, 16 },
        { 0, 24 }, { 8, 24 }, { 16, 24 }, { 24, 24 }, { 32, 24 }, { 40, 24 }, { 48, 24 }, { 56, 24 },
        { 0, 32 }, { 8, 32 }, { 16, 32 }, { 24, 32 }, { 32, 32 }, { 40, 32 }, { 48, 32 }, { 56, 32 },
        { 0, 40 }, { 8, 40 }, { 16, 40 }, { 24, 40 }, { 32, 40 }, { 40, 40 }, { 48, 40 }, { 56, 40 },
        { 0, 48 }, { 8, 48 }, { 16, 48 }, { 24, 48 }, { 32, 48 }, { 40, 48 }, { 48, 48 }, { 56, 48 },
        { 0, 56 }, { 8, 56 }, { 16, 56 }, { 24, 56 }, { 32, 56 }, { 40, 56 }, { 48, 56 }, { 56, 56 },
#if !NSQ_ME_CONTEXT_CLEAN_UP
        //H  Partitions
        { 0, 0 },
        { 0, 32 },
        { 0, 0 }, { 32, 0 },
        { 0, 16 }, { 32, 16 },
        { 0, 32 }, { 32, 32 },
        { 0, 48 }, { 32, 48 },
        { 0, 0 }, { 16, 0 }, { 32, 0 }, { 48, 0 },
        { 0, 8 }, { 16, 8 }, { 32, 8 }, { 48, 8 },
        { 0, 16 }, { 16, 16 }, { 32, 16 }, { 48, 16 },
        { 0, 24 }, { 16, 24 }, { 32, 24 }, { 48, 24 },
        { 0, 32 }, { 16, 32 }, { 32, 32 }, { 48, 32 },
        { 0, 40 }, { 16, 40 }, { 32, 40 }, { 48, 40 },
        { 0, 48 }, { 16, 48 }, { 32, 48 }, { 48, 48 },
        { 0, 56 }, { 16, 56 }, { 32, 56 }, { 48, 56 },
        //V  Partitions
        { 0, 0 }, { 32, 0 },
        { 0, 0 }, { 16, 0 }, { 32, 0 }, { 48, 0 },
        { 0, 32 }, { 16, 32 }, { 32, 32 }, { 48, 32 },
        { 0, 0 }, { 8, 0 }, { 16, 0 }, { 24, 0 }, { 32, 0 }, { 40, 0 }, { 48, 0 }, { 56, 0 },
        { 0, 16 }, { 8, 16 }, { 16, 16 }, { 24, 16 }, { 32, 16 }, { 40, 16 }, { 48, 16 }, { 56, 16 },
        { 0, 32 }, { 8, 32 }, { 16, 32 }, { 24, 32 }, { 32, 32 }, { 40, 32 }, { 48, 32 }, { 56, 32 },
        { 0, 48 }, { 8, 48 }, { 16, 48 }, { 24, 48 }, { 32, 48 }, { 40, 48 }, { 48, 48 }, { 56, 48 },
        // H4 Partitions
        { 0, 0 },  { 32, 0 },
        { 0, 8 },  { 32, 8 },
        { 0, 16 }, { 32, 16 },
        { 0, 24 }, { 32, 24 },
        { 0, 32 }, { 32, 32 },
        { 0, 40 }, { 32, 40 },
        { 0, 48 }, { 32, 48 },
        { 0, 56 }, { 32, 56 },
        // V4 Partitions
        { 0, 0  }, { 8, 0  }, { 16, 0  }, { 24, 0  }, { 32, 0  }, { 40, 0  }, { 48, 0  }, { 56, 0 },
        { 0, 32 }, { 8, 32 }, { 16, 32 }, { 24, 32 }, { 32, 32 }, { 40, 32 }, { 48, 32 }, { 56, 32 },
        { 0, 0},
        { 0, 16},
        { 0, 32},
        { 0, 48},
        { 0, 0}, { 16, 0},{ 32, 0}, { 48, 0}
#endif
    };

    static const uint8_t sub_position_type[16] = { 0, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2 };

    extern uint32_t compute8x4_sad_kernel_c(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride);
#if !REMOVE_UNUSED_CODE
    /*******************************************
    * GetEightHorizontalSearchPointResults_8x8_16x16_PU
    *******************************************/
    extern void get_eight_horizontal_search_point_results_8x8_16x16_pu_c(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad_8x8,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_sad_16x16,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint16_t  *p_sad16x16,
        EbBool     sub_sad);

    /*******************************************
    Calculate SAD for 32x32,64x64 from 16x16
    and check if there is improvement, if yes keep
    the best SAD+MV
    *******************************************/
    extern void get_eight_horizontal_search_point_results_32x32_64x64_pu_c(
        uint16_t  *p_sad16x16,
        uint32_t  *p_best_sad_32x32,
        uint32_t  *p_best_sad_64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    /*******************************************
    Calculate SAD for 16x16 and its 8x8 sublcoks
    and check if there is improvment, if yes keep
    the best SAD+MV
    *******************************************/
    extern void sad_calculation_8x8_16x16_c(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad_8x8,
        uint32_t  *p_best_sad_16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint32_t  *p_sad16x16,
        EbBool     sub_sad);

    /*******************************************
    Calculate SAD for 32x32,64x64 from 16x16
    and check if there is improvment, if yes keep
    the best SAD+MV
    *******************************************/
    extern void sad_calculation_32x32_64x64_c(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad_32x32,
        uint32_t  *p_best_sad_64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

#endif
    extern void ext_all_sad_calculation_8x8_16x16_c(
        uint8_t *src,
        uint32_t src_stride,
        uint8_t *ref,
        uint32_t ref_stride,
        uint32_t mv,
        uint32_t *p_best_sad_8x8,
        uint32_t *p_best_sad_16x16,
        uint32_t *p_best_mv8x8,
        uint32_t *p_best_mv16x16,
        uint32_t p_eight_sad16x16[16][8],
        uint32_t p_eight_sad8x8[64][8]);

#if !SHUT_ME_NSQ_SEARCH
    /****************************************************
    Calculate SAD for Rect H, V and H4, V4 partitions
    and update its Motion info if the result SAD is better
    ****************************************************/
    extern void ext_eigth_sad_calculation_nsq_c(
        uint32_t p_sad8x8[64][8],
        uint32_t p_sad16x16[16][8],
        uint32_t p_sad32x32[4][8],
        uint32_t *p_best_sad_64x32,
        uint32_t *p_best_mv64x32,
        uint32_t *p_best_sad_32x16,
        uint32_t *p_best_mv32x16,
        uint32_t *p_best_sad_16x8,
        uint32_t *p_best_mv16x8,
        uint32_t *p_best_sad_32x64,
        uint32_t *p_best_mv32x64,
        uint32_t *p_best_sad_16x32,
        uint32_t *p_best_mv16x32,
        uint32_t *p_best_sad_8x16,
        uint32_t *p_best_mv8x16,
        uint32_t *p_best_sad_32x8,
        uint32_t *p_best_mv32x8,
        uint32_t *p_best_sad_8x32,
        uint32_t *p_best_mv8x32,
        uint32_t *p_best_sad_64x16,
        uint32_t *p_best_mv64x16,
        uint32_t *p_best_sad_16x64,
        uint32_t *p_best_mv16x64,
        uint32_t mv);
#endif

    /*******************************************
    Calculate SAD for 32x32,64x64 from 16x16
    and check if there is improvment, if yes keep
    the best SAD+MV
    *******************************************/
    extern void ext_eight_sad_calculation_32x32_64x64_c(
        uint32_t p_sad16x16[16][8],
        uint32_t *p_best_sad_32x32,
        uint32_t *p_best_sad_64x64,
        uint32_t *p_best_mv32x32,
        uint32_t *p_best_mv64x64,
        uint32_t mv,
        uint32_t p_sad32x32[4][8]);

    // Nader - to be replaced by loock-up table
    /*******************************************
    * get_me_info_index
    *   search the correct index of the motion
    *   info that corresponds to the input
    *   md candidate
    *******************************************/
    extern uint32_t get_me_info_index(
        uint32_t         max_me_block,
        const BlockGeom *blk_geom,
        uint32_t         geom_offset_x,
        uint32_t         geom_offset_y);
#if !REMOVE_ME_SUBPEL_CODE
    void half_pel_refinement_sb(
#if !SHUT_ME_NSQ_SEARCH
        PictureParentControlSet *pcs_ptr,
#endif
        MeContext *context_ptr,  // input/output parameter, ME context Ptr, used
                                 // to get/update ME results
        uint8_t *refBuffer, uint32_t ref_stride,
        uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated
                                // search area Ptr
        uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated
                                // search area Ptr
        uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated
                                // search area Ptr
        int16_t x_search_area_origin,  // input parameter, search area origin in
                                       // the horizontal direction, used to
                                       // point to reference samples
        int16_t y_search_area_origin,  // input parameter, search area origin in
                                       // the vertical direction, used to point
                                       // to reference samples
        uint32_t search_area_height,  // input parameter, search area height
        uint32_t search_area_width,  // input parameter, search area width
#if !SHUT_ME_NSQ_SEARCH
        uint8_t list_index, // reference picture list
        uint8_t ref_pic_index, // reference picture index
#endif
        uint32_t integer_mv);         // input parameter, integer MV
#endif

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimation_h
// clang-format on
