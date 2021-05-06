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

#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"

#include "EbSourceBasedOperationsProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbPictureDemuxResults.h"
#ifdef ARCH_X86_64
#include <emmintrin.h>
#endif
#include "EbEncHandle.h"
#include "EbUtility.h"
#include "EbPictureManagerProcess.h"
#include "EbReferenceObject.h"
#include "EbTransforms.h"
#include "aom_dsp_rtcd.h"
#include "EbLog.h"
#include "EbIntraPrediction.h"
#include "EbMotionEstimation.h"
#include "EbEncDecResults.h"
#include "EbRateDistortionCost.h"

/**************************************
 * Context
 **************************************/

typedef struct SourceBasedOperationsContext {
    EbDctor dctor;
    EbFifo *initial_rate_control_results_input_fifo_ptr;
    EbFifo *picture_demux_results_output_fifo_ptr;
    EbFifo *sbo_output_fifo_ptr;
    // local zz cost array
    uint32_t complete_sb_count;
    uint8_t *y_mean_ptr;
    uint8_t *cr_mean_ptr;
    uint8_t *cb_mean_ptr;
} SourceBasedOperationsContext;
typedef struct TplDispenserContext {
    EbDctor dctor;
    EbFifo *tpl_disp_input_fifo_ptr;
    EbFifo *tpl_disp_fb_fifo_ptr;
    uint32_t sb_index;
    uint32_t coded_sb_count;
} TplDispenserContext;

static void source_based_operations_context_dctor(EbPtr p) {
    EbThreadContext *             thread_context_ptr = (EbThreadContext *)p;
    SourceBasedOperationsContext *obj = (SourceBasedOperationsContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/************************************************
* Source Based Operation Context Constructor
************************************************/
EbErrorType source_based_operations_context_ctor(EbThreadContext *  thread_context_ptr,
                                                 const EbEncHandle *enc_handle_ptr, int index) {
    SourceBasedOperationsContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = source_based_operations_context_dctor;

    context_ptr->initial_rate_control_results_input_fifo_ptr =
        svt_system_resource_get_consumer_fifo(
            enc_handle_ptr->initial_rate_control_results_resource_ptr, index);

    context_ptr->sbo_output_fifo_ptr= svt_system_resource_get_producer_fifo(
        enc_handle_ptr->tpl_disp_res_srm, index);
    context_ptr->picture_demux_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, index);
    return EB_ErrorNone;
}

/***************************************************
* Derives BEA statistics and set activity flags
***************************************************/
void derive_picture_activity_statistics(PictureParentControlSet *pcs_ptr)

{
    uint64_t non_moving_index_min = ~0u;
    uint64_t non_moving_index_max = 0;
    uint64_t non_moving_index_sum = 0;
    uint32_t complete_sb_count    = 0;
    uint32_t non_moving_sb_count  = 0;
    uint32_t sb_total_count       = pcs_ptr->sb_total_count;

    for (uint32_t sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        SbParams *sb_params = &pcs_ptr->sb_params_array[sb_index];
        if (sb_params->is_complete_sb) {
            non_moving_index_min = pcs_ptr->non_moving_index_array[sb_index] < non_moving_index_min
                ? pcs_ptr->non_moving_index_array[sb_index]
                : non_moving_index_min;

            non_moving_index_max = pcs_ptr->non_moving_index_array[sb_index] > non_moving_index_max
                ? pcs_ptr->non_moving_index_array[sb_index]
                : non_moving_index_max;
            if (pcs_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1)
                non_moving_sb_count++;
            complete_sb_count++;

            non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_index];
        }
    }

    if (complete_sb_count > 0) {
        pcs_ptr->non_moving_index_average = (uint16_t)(non_moving_index_sum / complete_sb_count);
        pcs_ptr->kf_zeromotion_pct        = (non_moving_sb_count * 100) / complete_sb_count;
    }
    pcs_ptr->non_moving_index_min_distance = (uint16_t)(
        ABS((int32_t)(pcs_ptr->non_moving_index_average) - (int32_t)non_moving_index_min));
    pcs_ptr->non_moving_index_max_distance = (uint16_t)(
        ABS((int32_t)(pcs_ptr->non_moving_index_average) - (int32_t)non_moving_index_max));
    return;
}

/*
     TPL dispenser context dctor
*/
static void tpl_disp_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    TplDispenserContext *obj = (TplDispenserContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}
/*
     TPL dispenser context cctor
*/
EbErrorType tpl_disp_context_ctor(EbThreadContext *  thread_context_ptr,
        const EbEncHandle *enc_handle_ptr, int index, int tasks_index) {
    TplDispenserContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);

    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = tpl_disp_context_dctor;

    context_ptr->tpl_disp_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->tpl_disp_res_srm, index);

    context_ptr->tpl_disp_fb_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->tpl_disp_res_srm, tasks_index);

    return EB_ErrorNone;
}


void tpl_prep_info(PictureParentControlSet    *pcs) ;


// Generate lambda factor to tune lambda based on TPL stats
static void generate_lambda_scaling_factor(PictureParentControlSet *pcs_ptr,
                                           int64_t                  mc_dep_cost_base) {
    Av1Common *cm         = pcs_ptr->av1_cm;
    const int  step       = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int  mi_cols_sr = ((pcs_ptr->aligned_width + 15) / 16) << 2;

    const int    block_size = BLOCK_16X16;
    const int    num_mi_w   = mi_size_wide[block_size];
    const int    num_mi_h   = mi_size_high[block_size];
    const int    num_cols   = (mi_cols_sr + num_mi_w - 1) / num_mi_w;
    const int    num_rows   = (cm->mi_rows + num_mi_h - 1) / num_mi_h;
    const int    stride     = mi_cols_sr >> (1 + pcs_ptr->is_720p_or_larger);
    const double c          = 1.2;

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            double    intra_cost  = 0.0;
            double    mc_dep_cost = 0.0;
            const int index       = row * num_cols + col;
            for (int mi_row = row * num_mi_h; mi_row < (row + 1) * num_mi_h; mi_row += step) {
                for (int mi_col = col * num_mi_w; mi_col < (col + 1) * num_mi_w; mi_col += step) {
                    if (mi_row >= cm->mi_rows || mi_col >= mi_cols_sr)
                        continue;

                    const int index1 = (mi_row >> (1 + pcs_ptr->is_720p_or_larger)) * stride +
                        (mi_col >> (1 + pcs_ptr->is_720p_or_larger));
                    TplStats *tpl_stats_ptr = pcs_ptr->tpl_stats[index1];
                    int64_t   mc_dep_delta  = RDCOST(pcs_ptr->base_rdmult,
                                                  tpl_stats_ptr->mc_dep_rate,
                                                  tpl_stats_ptr->mc_dep_dist);
                    intra_cost += (double)(tpl_stats_ptr->recrf_dist << RDDIV_BITS);
                    mc_dep_cost += (double)(tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
                }
            }
            double rk = 0;
            if (mc_dep_cost > 0 && intra_cost > 0) {
                rk = intra_cost / mc_dep_cost;
            }

            pcs_ptr->tpl_rdmult_scaling_factors[index] = (mc_dep_cost_base) ? rk / pcs_ptr->r0 + c
                                                                            : c;
        }
    }

    return;
}

static AOM_INLINE void get_quantize_error(MacroblockPlane *p, const TranLow *coeff, TranLow *qcoeff,
                                          TranLow *dqcoeff, TxSize tx_size, uint16_t *eob,
                                          int64_t *recon_error, int64_t *sse) {
    const ScanOrder *const scan_order =
        &av1_scan_orders[tx_size][DCT_DCT]; //&av1_default_scan_orders[tx_size]
    int       pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
    const int shift   = tx_size == TX_32X32 ? 0 : 2;

    svt_av1_quantize_fp(coeff,
                        pix_num,
                        p->zbin_qtx,
                        p->round_fp_qtx,
                        p->quant_fp_qtx,
                        p->quant_shift_qtx,
                        qcoeff,
                        dqcoeff,
                        p->dequant_qtx,
                        eob,
                        scan_order->scan,
                        scan_order->iscan);

    *recon_error = svt_av1_block_error(coeff, dqcoeff, pix_num, sse) >> shift;
    *recon_error = AOMMAX(*recon_error, 1);

    *sse = (*sse) >> shift;
    *sse = AOMMAX(*sse, 1);
}

static int rate_estimator(TranLow *qcoeff, int eob, TxSize tx_size) {
    const ScanOrder *const scan_order =
        &av1_scan_orders[tx_size][DCT_DCT]; //&av1_default_scan_orders[tx_size]

    assert((1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]]) >= eob);

    int rate_cost = 1;

    for (int idx = 0; idx < eob; ++idx) {
        int abs_level = abs(qcoeff[scan_order->scan[idx]]);
        rate_cost += (int)(log1p(abs_level) / log(2.0)) + 1;
    }

    return (rate_cost << AV1_PROB_COST_SHIFT);
}




static void result_model_store(PictureParentControlSet *pcs_ptr, TplStats  *tpl_stats_ptr,
                               uint32_t mb_origin_x, uint32_t mb_origin_y) {
    const int mi_height       = mi_size_high[BLOCK_16X16];
    const int mi_width        = mi_size_wide[BLOCK_16X16];
    const int step            = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int shift           = 3 + pcs_ptr->is_720p_or_larger;
    const int aligned16_width = ((pcs_ptr->aligned_width + 15) / 16) << 4;

    int64_t srcrf_dist = tpl_stats_ptr->srcrf_dist / (mi_height * mi_width);
    int64_t recrf_dist = tpl_stats_ptr->recrf_dist / (mi_height * mi_width);
    int64_t srcrf_rate = tpl_stats_ptr->srcrf_rate / (mi_height * mi_width);
    int64_t recrf_rate = tpl_stats_ptr->recrf_rate / (mi_height * mi_width);

    srcrf_dist = AOMMAX(1, srcrf_dist);
    recrf_dist = AOMMAX(1, recrf_dist);
    srcrf_rate = AOMMAX(1, srcrf_rate);
    recrf_rate = AOMMAX(1, recrf_rate);

    for (int idy = 0; idy < mi_height; idy += step) {
        TplStats *dst_ptr =
            pcs_ptr->tpl_stats[((mb_origin_y >> shift) + (idy >> 1)) * (aligned16_width >> shift) +
                               (mb_origin_x >> shift)];
        for (int idx = 0; idx < mi_width; idx += step) {
            dst_ptr->srcrf_dist    = srcrf_dist;
            dst_ptr->recrf_dist    = recrf_dist;
            dst_ptr->srcrf_rate    = srcrf_rate;
            dst_ptr->recrf_rate    = recrf_rate;
            dst_ptr->mv            = tpl_stats_ptr->mv;
            dst_ptr->ref_frame_poc = tpl_stats_ptr->ref_frame_poc;
            ++dst_ptr;
        }
    }
}


static const int16_t dc_qlookup_QTX[QINDEX_RANGE] = {
    4,   8,   8,   9,   10,  11,  12,  12,  13,   14,   15,   16,   17,   18,   19,   19,
    20,  21,  22,  23,  24,  25,  26,  26,  27,   28,   29,   30,   31,   32,   32,   33,
    34,  35,  36,  37,  38,  38,  39,  40,  41,   42,   43,   43,   44,   45,   46,   47,
    48,  48,  49,  50,  51,  52,  53,  53,  54,   55,   56,   57,   57,   58,   59,   60,
    61,  62,  62,  63,  64,  65,  66,  66,  67,   68,   69,   70,   70,   71,   72,   73,
    74,  74,  75,  76,  77,  78,  78,  79,  80,   81,   81,   82,   83,   84,   85,   85,
    87,  88,  90,  92,  93,  95,  96,  98,  99,   101,  102,  104,  105,  107,  108,  110,
    111, 113, 114, 116, 117, 118, 120, 121, 123,  125,  127,  129,  131,  134,  136,  138,
    140, 142, 144, 146, 148, 150, 152, 154, 156,  158,  161,  164,  166,  169,  172,  174,
    177, 180, 182, 185, 187, 190, 192, 195, 199,  202,  205,  208,  211,  214,  217,  220,
    223, 226, 230, 233, 237, 240, 243, 247, 250,  253,  257,  261,  265,  269,  272,  276,
    280, 284, 288, 292, 296, 300, 304, 309, 313,  317,  322,  326,  330,  335,  340,  344,
    349, 354, 359, 364, 369, 374, 379, 384, 389,  395,  400,  406,  411,  417,  423,  429,
    435, 441, 447, 454, 461, 467, 475, 482, 489,  497,  505,  513,  522,  530,  539,  549,
    559, 569, 579, 590, 602, 614, 626, 640, 654,  668,  684,  700,  717,  736,  755,  775,
    796, 819, 843, 869, 896, 925, 955, 988, 1022, 1058, 1098, 1139, 1184, 1232, 1282, 1336,
};

static const int16_t dc_qlookup_10_QTX[QINDEX_RANGE] = {
    4,    9,    10,   13,   15,   17,   20,   22,   25,   28,   31,   34,   37,   40,   43,   47,
    50,   53,   57,   60,   64,   68,   71,   75,   78,   82,   86,   90,   93,   97,   101,  105,
    109,  113,  116,  120,  124,  128,  132,  136,  140,  143,  147,  151,  155,  159,  163,  166,
    170,  174,  178,  182,  185,  189,  193,  197,  200,  204,  208,  212,  215,  219,  223,  226,
    230,  233,  237,  241,  244,  248,  251,  255,  259,  262,  266,  269,  273,  276,  280,  283,
    287,  290,  293,  297,  300,  304,  307,  310,  314,  317,  321,  324,  327,  331,  334,  337,
    343,  350,  356,  362,  369,  375,  381,  387,  394,  400,  406,  412,  418,  424,  430,  436,
    442,  448,  454,  460,  466,  472,  478,  484,  490,  499,  507,  516,  525,  533,  542,  550,
    559,  567,  576,  584,  592,  601,  609,  617,  625,  634,  644,  655,  666,  676,  687,  698,
    708,  718,  729,  739,  749,  759,  770,  782,  795,  807,  819,  831,  844,  856,  868,  880,
    891,  906,  920,  933,  947,  961,  975,  988,  1001, 1015, 1030, 1045, 1061, 1076, 1090, 1105,
    1120, 1137, 1153, 1170, 1186, 1202, 1218, 1236, 1253, 1271, 1288, 1306, 1323, 1342, 1361, 1379,
    1398, 1416, 1436, 1456, 1476, 1496, 1516, 1537, 1559, 1580, 1601, 1624, 1647, 1670, 1692, 1717,
    1741, 1766, 1791, 1817, 1844, 1871, 1900, 1929, 1958, 1990, 2021, 2054, 2088, 2123, 2159, 2197,
    2236, 2276, 2319, 2363, 2410, 2458, 2508, 2561, 2616, 2675, 2737, 2802, 2871, 2944, 3020, 3102,
    3188, 3280, 3375, 3478, 3586, 3702, 3823, 3953, 4089, 4236, 4394, 4559, 4737, 4929, 5130, 5347,
};

static const int16_t dc_qlookup_12_QTX[QINDEX_RANGE] = {
    4,     12,    18,    25,    33,    41,    50,    60,    70,    80,    91,    103,   115,
    127,   140,   153,   166,   180,   194,   208,   222,   237,   251,   266,   281,   296,
    312,   327,   343,   358,   374,   390,   405,   421,   437,   453,   469,   484,   500,
    516,   532,   548,   564,   580,   596,   611,   627,   643,   659,   674,   690,   706,
    721,   737,   752,   768,   783,   798,   814,   829,   844,   859,   874,   889,   904,
    919,   934,   949,   964,   978,   993,   1008,  1022,  1037,  1051,  1065,  1080,  1094,
    1108,  1122,  1136,  1151,  1165,  1179,  1192,  1206,  1220,  1234,  1248,  1261,  1275,
    1288,  1302,  1315,  1329,  1342,  1368,  1393,  1419,  1444,  1469,  1494,  1519,  1544,
    1569,  1594,  1618,  1643,  1668,  1692,  1717,  1741,  1765,  1789,  1814,  1838,  1862,
    1885,  1909,  1933,  1957,  1992,  2027,  2061,  2096,  2130,  2165,  2199,  2233,  2267,
    2300,  2334,  2367,  2400,  2434,  2467,  2499,  2532,  2575,  2618,  2661,  2704,  2746,
    2788,  2830,  2872,  2913,  2954,  2995,  3036,  3076,  3127,  3177,  3226,  3275,  3324,
    3373,  3421,  3469,  3517,  3565,  3621,  3677,  3733,  3788,  3843,  3897,  3951,  4005,
    4058,  4119,  4181,  4241,  4301,  4361,  4420,  4479,  4546,  4612,  4677,  4742,  4807,
    4871,  4942,  5013,  5083,  5153,  5222,  5291,  5367,  5442,  5517,  5591,  5665,  5745,
    5825,  5905,  5984,  6063,  6149,  6234,  6319,  6404,  6495,  6587,  6678,  6769,  6867,
    6966,  7064,  7163,  7269,  7376,  7483,  7599,  7715,  7832,  7958,  8085,  8214,  8352,
    8492,  8635,  8788,  8945,  9104,  9275,  9450,  9639,  9832,  10031, 10245, 10465, 10702,
    10946, 11210, 11482, 11776, 12081, 12409, 12750, 13118, 13501, 13913, 14343, 14807, 15290,
    15812, 16356, 16943, 17575, 18237, 18949, 19718, 20521, 21387,
};

int16_t av1_dc_quant_qtx(int qindex, int delta, AomBitDepth bit_depth) {
    const int q_clamped = clamp(qindex + delta, 0, MAXQ);
    switch (bit_depth) {
    case AOM_BITS_8: return dc_qlookup_QTX[q_clamped];
    case AOM_BITS_10: return dc_qlookup_10_QTX[q_clamped];
    case AOM_BITS_12: return dc_qlookup_12_QTX[q_clamped];
    default: assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12"); return -1;
    }
}

int svt_av1_compute_rd_mult_based_on_qindex(AomBitDepth bit_depth, int qindex) {
    const int q = av1_dc_quant_qtx(qindex, 0, bit_depth);
    //const int q = svt_av1_dc_quant_Q3(qindex, 0, bit_depth);
    int rdmult = q * q;
    rdmult     = rdmult * 3 + (rdmult * 2 / 3);
    switch (bit_depth) {
    case AOM_BITS_8: break;
    case AOM_BITS_10: rdmult = ROUND_POWER_OF_TWO(rdmult, 4); break;
    case AOM_BITS_12: rdmult = ROUND_POWER_OF_TWO(rdmult, 8); break;
    default: assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12"); return -1;
    }
    return rdmult > 0 ? rdmult : 1;
}

void svt_av1_build_quantizer(AomBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q,
                             int32_t u_ac_delta_q, int32_t v_dc_delta_q, int32_t v_ac_delta_q,
                             Quants *const quants, Dequants *const deq);

double svt_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);

int32_t svt_av1_compute_qdelta(double qstart, double qtarget, AomBitDepth bit_depth);

extern void filter_intra_edge(OisMbResults *ois_mb_results_ptr, uint8_t mode,
                              uint16_t max_frame_width, uint16_t max_frame_height, int32_t p_angle,
                              int32_t cu_origin_x, int32_t cu_origin_y, uint8_t *above_row,
                              uint8_t *left_col);

//Given one reference frame identified by the pair (list_index,ref_index)
//indicate if ME data is valid
static uint8_t is_me_data_valid(const MeSbResults *me_results, uint32_t me_mb_offset,
                                uint8_t list_idx, uint8_t ref_idx) {
    uint8_t            total_me_cnt = me_results->total_me_candidate_index[me_mb_offset];
    const MeCandidate *me_block_results =
        &me_results->me_candidate_array[me_mb_offset * MAX_PA_ME_CAND];

    for (uint32_t me_cand_i = 0; me_cand_i < total_me_cnt; ++me_cand_i) {
        const MeCandidate *me_cand = &me_block_results[me_cand_i];
        assert(/*me_cand->direction >= 0 && */ me_cand->direction <= 2);
        if (me_cand->direction == 0 || me_cand->direction == 2) {
            if (list_idx == me_cand->ref0_list && ref_idx == me_cand->ref_idx_l0)
                return 1;
        }
        if (me_cand->direction == 1 || me_cand->direction == 2) {
            if (list_idx == me_cand->ref1_list && ref_idx == me_cand->ref_idx_l1)
                return 1;
        }
    }
    return 0;
}


void clip_mv_in_pad(
    EbPictureBufferDesc *ref_pic_ptr,
    uint32_t mb_origin_x,
    uint32_t mb_origin_y,
    int16_t *x_curr_mv,
    int16_t *y_curr_mv)
{
    // Search area adjustment
    int16_t blk_origin_x = mb_origin_x;
    int16_t blk_origin_y = mb_origin_y;
    int16_t bwidth = 16;
    int16_t bheight = 16;
    int16_t mvx = *x_curr_mv;
    int16_t mvy = *y_curr_mv;
    int16_t padx = TPL_PADX;
    int16_t pady = TPL_PADY;

    if ((blk_origin_x + (mvx >> 3)) < -padx)
        mvx = (-padx - blk_origin_x) << 3;

    if ((blk_origin_x + bwidth + (mvx >> 3)) > (padx + ref_pic_ptr->max_width - 1))
        mvx = ((padx + ref_pic_ptr->max_width - 1) - (blk_origin_x + bwidth)) << 3;

    if ((blk_origin_y + (mvy >> 3)) < -pady)
        mvy = (-pady - blk_origin_y) << 3;

    if ((blk_origin_y + bheight + (mvy >> 3)) > (pady + ref_pic_ptr->max_height - 1))
        mvy = ((pady + ref_pic_ptr->max_height - 1) - (blk_origin_y + bheight)) << 3;

    *x_curr_mv = mvx;
    *y_curr_mv = mvy;
}
// Reference pruning, Loop over all available references and get the best reference idx based on SAD
void get_best_reference(
    PictureParentControlSet *pcs_ptr,
    uint32_t sb_index,
    uint32_t   me_mb_offset,
    uint32_t mb_origin_x,
    uint32_t mb_origin_y,
    uint32_t *best_reference )
{
    EbPictureBufferDesc *input_ptr     = pcs_ptr->enhanced_picture_ptr;
    uint32_t             max_inter_ref = MAX_PA_ME_MV;
    EbPictureBufferDesc *ref_pic_ptr;
    int16_t              x_curr_mv          = 0;
    int16_t              y_curr_mv          = 0;
    uint32_t             best_reference_sad = UINT32_MAX;
    uint32_t             reference_sad;
    uint8_t *            src_mb = input_ptr->buffer_y + input_ptr->origin_x + mb_origin_x +
        (input_ptr->origin_y + mb_origin_y) * input_ptr->stride_y;

    for (uint32_t rf_idx = 0; rf_idx < max_inter_ref; rf_idx++) {
        uint32_t list_index    = rf_idx < 4 ? 0 : 1;
        uint32_t ref_pic_index = rf_idx >= 4 ? (rf_idx - 4) : rf_idx;
        if ((list_index == 0 && (ref_pic_index + 1) > pcs_ptr->tpl_data.tpl_ref0_count) ||
            (list_index == 1 && (ref_pic_index + 1) > pcs_ptr->tpl_data.tpl_ref1_count))
            continue;
        if (!is_me_data_valid(
                pcs_ptr->pa_me_data->me_results[sb_index], me_mb_offset, list_index, ref_pic_index))
            continue;
        ref_pic_ptr = (EbPictureBufferDesc *)pcs_ptr->tpl_data
                          .tpl_ref_ds_ptr_array[list_index][ref_pic_index]
                          .picture_ptr;

        const MeSbResults *me_results = pcs_ptr->pa_me_data->me_results[sb_index];
        x_curr_mv =
            me_results
                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index]
                .x_mv
            << 1;
        y_curr_mv =
            me_results
                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index]
                .y_mv
            << 1;
        clip_mv_in_pad(ref_pic_ptr,mb_origin_x,mb_origin_y,&x_curr_mv,&y_curr_mv);
        MV      best_mv          = {y_curr_mv, x_curr_mv};
        int32_t ref_origin_index = ref_pic_ptr->origin_x + (mb_origin_x + (best_mv.col >> 3)) +
            (mb_origin_y + (best_mv.row >> 3) + ref_pic_ptr->origin_y) * ref_pic_ptr->stride_y;
        reference_sad = svt_nxm_sad_kernel_sub_sampled(src_mb,
                                                       input_ptr->stride_y,
                                                       ref_pic_ptr->buffer_y + ref_origin_index,
                                                       ref_pic_ptr->stride_y,
                                                       16,
                                                       16);
        if (reference_sad < best_reference_sad) {
            best_reference_sad = reference_sad;
            *best_reference    = rf_idx;
        }
    }
    return;
}





/*
    TPL Dispenser SB based (sz 64x64)
*/
void tpl_mc_flow_dispenser_sb(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *scs_ptr,
    PictureParentControlSet        *pcs_ptr,
    int32_t                          frame_idx,
    uint32_t                        sb_index,
    int32_t                         qIndex)
{
{
    uint32_t             picture_width_in_mb = (pcs_ptr->enhanced_picture_ptr->width + 16 - 1) / 16;
    int16_t              x_curr_mv           = 0;
    int16_t              y_curr_mv           = 0;
    uint32_t             me_mb_offset        = 0;
    TxSize               tx_size             = TX_16X16;
    EbPictureBufferDesc *ref_pic_ptr;
    BlockGeom            blk_geom;
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc *recon_picture_ptr =
        encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx];
    TplStats tpl_stats;

    DECLARE_ALIGNED(32, uint8_t, predictor8[256 * 2]);
    DECLARE_ALIGNED(32, int16_t, src_diff[256]);
    DECLARE_ALIGNED(32, TranLow, coeff[256]);
    DECLARE_ALIGNED(32, TranLow, qcoeff[256]);
    DECLARE_ALIGNED(32, TranLow, dqcoeff[256]);
    DECLARE_ALIGNED(32, TranLow, best_coeff[256]);
    uint8_t *predictor = predictor8;

    blk_geom.bwidth  = 16;
    blk_geom.bheight = 16;

    MacroblockPlane mb_plane;
    mb_plane.quant_qtx       = scs_ptr->quants_8bit.y_quant[qIndex];
    mb_plane.quant_fp_qtx    = scs_ptr->quants_8bit.y_quant_fp[qIndex];
    mb_plane.round_fp_qtx    = scs_ptr->quants_8bit.y_round_fp[qIndex];
    mb_plane.quant_shift_qtx = scs_ptr->quants_8bit.y_quant_shift[qIndex];
    mb_plane.zbin_qtx        = scs_ptr->quants_8bit.y_zbin[qIndex];
    mb_plane.round_qtx       = scs_ptr->quants_8bit.y_round[qIndex];
    mb_plane.dequant_qtx     = scs_ptr->deq_8bit.y_dequant_qtx[qIndex];

    EbPictureBufferDesc *input_ptr = pcs_ptr->enhanced_picture_ptr;
    const uint8_t tpl_opt_flag = pcs_ptr->tpl_ctrls.tpl_opt_flag;


    SbParams *sb_params    = &scs_ptr->sb_params_array[sb_index];
    uint32_t  pa_blk_index = 0;
    while (pa_blk_index < CU_MAX_COUNT) {
        const CodedBlockStats *blk_stats_ptr;
        blk_stats_ptr              = get_coded_blk_stats(pa_blk_index);
        uint8_t bsize              = blk_stats_ptr->size;
        EbBool  small_boundary_blk = EB_FALSE;

        {
            uint32_t cu_origin_x = sb_params->origin_x + blk_stats_ptr->origin_x;
            uint32_t cu_origin_y = sb_params->origin_y + blk_stats_ptr->origin_y;
            if ((blk_stats_ptr->origin_x % 16) == 0 &&
                (blk_stats_ptr->origin_y % 16) == 0 &&
                ((pcs_ptr->enhanced_picture_ptr->width - cu_origin_x) < 16 ||
                    (pcs_ptr->enhanced_picture_ptr->height - cu_origin_y) < 16))
                small_boundary_blk = EB_TRUE;
        }
        if (bsize != 16 && !small_boundary_blk) {
            pa_blk_index++;
            continue;
        }
        if (sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[pa_blk_index]]) {
            uint32_t  mb_origin_x       = sb_params->origin_x + blk_stats_ptr->origin_x;
            uint32_t  mb_origin_y       = sb_params->origin_y + blk_stats_ptr->origin_y;
            const int dst_buffer_stride = recon_picture_ptr->stride_y;
            const int dst_mb_offset     = mb_origin_y * dst_buffer_stride + mb_origin_x;
            const int dst_basic_offset  = recon_picture_ptr->origin_y *
                    recon_picture_ptr->stride_y +
                recon_picture_ptr->origin_x;
            uint8_t *dst_buffer = recon_picture_ptr->buffer_y + dst_basic_offset +
                dst_mb_offset;

            int64_t  inter_cost;
            int64_t  recon_error = 1, sse = 1;
            uint64_t best_ref_poc    = 0;
            int32_t  best_rf_idx     = -1;
            int64_t  best_inter_cost = INT64_MAX;
            MV       final_best_mv   = {0, 0};
            uint32_t max_inter_ref   = MAX_PA_ME_MV;

            PredictionMode best_intra_mode = DC_PRED;
            int64_t        best_intra_cost = INT64_MAX;
            // Disable intra prediction
            uint8_t disable_intra_pred  = tpl_opt_flag && (pcs_ptr->tpl_ctrls.disable_intra_pred_nref ||
                pcs_ptr->tpl_ctrls.disable_intra_pred_nbase);
            if (!disable_intra_pred ||
                (pcs_ptr->tpl_ctrls.disable_intra_pred_nref && pcs_ptr->tpl_data.is_used_as_reference_flag) ||
                (pcs_ptr->tpl_ctrls.disable_intra_pred_nbase && pcs_ptr->tpl_data.tpl_temporal_layer_index == 0)){
                if (scs_ptr->in_loop_ois == 0) {
                    OisMbResults *ois_mb_results_ptr =
                        pcs_ptr->ois_mb_results[(mb_origin_y >> 4) * picture_width_in_mb +
                                                (mb_origin_x >> 4)];
                    best_intra_mode = ois_mb_results_ptr->intra_mode;
                    best_intra_cost = ois_mb_results_ptr->intra_cost;

                } else { // ois
                    // always process as block16x16 even bsize or tx_size is 8x8
                    bsize = 16;
                    DECLARE_ALIGNED(16, uint8_t, left0_data[MAX_TX_SIZE * 2 + 32]);
                    DECLARE_ALIGNED(16, uint8_t, above0_data[MAX_TX_SIZE * 2 + 32]);
                    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
                    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);

                    uint8_t *above_row;
                    uint8_t *left_col;
                    uint8_t *above0_row;
                    uint8_t *left0_col;
                    above0_row = above0_data + 16;
                    left0_col  = left0_data + 16;
                    above_row  = above_data + 16;
                    left_col   = left_data + 16;

                    uint8_t *src = input_ptr->buffer_y +
                        pcs_ptr->enhanced_picture_ptr->origin_x + mb_origin_x +
                        (pcs_ptr->enhanced_picture_ptr->origin_y + mb_origin_y) *
                            input_ptr->stride_y;

                    // Fill Neighbor Arrays
                    update_neighbor_samples_array_open_loop_mb(
                                                                1, // use_top_righ_bottom_left
                                                                1, // update_top_neighbor
                                                                above0_row - 1,
                                                                left0_col - 1,
                                                                input_ptr,
                                                                input_ptr->stride_y,
                                                                mb_origin_x,
                                                                mb_origin_y,
                                                                bsize,
                                                                bsize);

                    uint8_t ois_intra_mode;
                    uint8_t intra_mode_start = DC_PRED;
                    EbBool  enable_paeth  = pcs_ptr->scs_ptr->static_config.enable_paeth ==
                            DEFAULT
                            ? EB_TRUE
                            : (EbBool)pcs_ptr->scs_ptr->static_config.enable_paeth;
                    EbBool  enable_smooth = pcs_ptr->scs_ptr->static_config.enable_smooth ==
                            DEFAULT
                            ? EB_TRUE
                            : (EbBool)pcs_ptr->scs_ptr->static_config.enable_smooth;
                    uint8_t intra_mode_end =
                    pcs_ptr->tpl_ctrls.tpl_opt_flag

                        ? DC_PRED
                        : enable_paeth      ? PAETH_PRED
                            : enable_smooth ? SMOOTH_H_PRED
                                            : D67_PRED;

                    for (ois_intra_mode = intra_mode_start;
                            ois_intra_mode <= intra_mode_end;
                            ++ois_intra_mode) {
                        int32_t p_angle = av1_is_directional_mode(
                                                (PredictionMode)ois_intra_mode)
                            ? mode_to_angle_map[(PredictionMode)ois_intra_mode]
                            : 0;
                        // Edge filter
                        if (av1_is_directional_mode((PredictionMode)ois_intra_mode) &&
                            1 /*scs_ptr->seq_header.enable_intra_edge_filter*/) {
                            EB_MEMCPY(left_data,
                                        left0_data,
                                        sizeof(uint8_t) * (MAX_TX_SIZE * 2 + 32));
                            EB_MEMCPY(above_data,
                                        above0_data,
                                        sizeof(uint8_t) * (MAX_TX_SIZE * 2 + 32));
                            above_row = above_data + 16;
                            left_col  = left_data + 16;
                            filter_intra_edge(NULL,
                                                ois_intra_mode,
                                                scs_ptr->seq_header.max_frame_width,
                                                scs_ptr->seq_header.max_frame_height,
                                                p_angle,
                                                (int32_t)mb_origin_x,
                                                (int32_t)mb_origin_y,
                                                above_row,
                                                left_col);
                        } else {
                            above_row = above0_row;
                            left_col  = left0_col;
                        }
                        // PRED
                        intra_prediction_open_loop_mb(p_angle,
                                                        ois_intra_mode,
                                                        mb_origin_x,
                                                        mb_origin_y,
                                                        tx_size,
                                                        above_row,
                                                        left_col,
                                                        predictor,
                                                        16);

                        // Distortion
                        int64_t intra_cost;
                        if (pcs_ptr->tpl_ctrls.tpl_opt_flag && pcs_ptr->tpl_ctrls.use_pred_sad_in_intra_search) {
                            intra_cost = svt_nxm_sad_kernel_sub_sampled(
                                src,
                                input_ptr->stride_y,
                                predictor,
                                16,
                                16,
                                16);
                        }
                        else {
                            svt_aom_subtract_block(
                                16, 16, src_diff, 16, src, input_ptr->stride_y, predictor, 16);
                            EB_TRANS_COEFF_SHAPE pf_shape = pcs_ptr->tpl_ctrls.tpl_opt_flag ? pcs_ptr->tpl_ctrls.pf_shape : DEFAULT_SHAPE;
                            svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size, pf_shape, 8, 0);
                             intra_cost = svt_aom_satd(coeff, 16 * 16);
                        }

                        if (intra_cost < best_intra_cost) {
                            best_intra_cost = intra_cost;
                            best_intra_mode = ois_intra_mode;
                        }
                    }
                }
            }
            uint8_t  best_mode = DC_PRED;
            uint8_t *src_mb    = input_picture_ptr->buffer_y + input_picture_ptr->origin_x +
                mb_origin_x +
                (input_picture_ptr->origin_y + mb_origin_y) * input_picture_ptr->stride_y;
            memset(&tpl_stats, 0, sizeof(tpl_stats));
            blk_geom.origin_x = blk_stats_ptr->origin_x;
            blk_geom.origin_y = blk_stats_ptr->origin_y;
            me_mb_offset      = get_me_info_index(
                pcs_ptr->max_number_of_pus_per_sb, &blk_geom, 0, 0);

            uint32_t best_reference = 0;
            if (pcs_ptr->tpl_ctrls.tpl_opt_flag && pcs_ptr->tpl_ctrls.get_best_ref)
                // Reference pruning
                get_best_reference(pcs_ptr,
                                    sb_index,
                                    me_mb_offset,
                                    mb_origin_x,
                                    mb_origin_y,
                                    &best_reference);

            for (uint32_t rf_idx = 0; rf_idx < max_inter_ref; rf_idx++) {
                if (pcs_ptr->tpl_ctrls.get_best_ref)
                    if (rf_idx != best_reference)
                        continue;
                uint32_t list_index    = rf_idx < 4 ? 0 : 1;
                uint32_t ref_pic_index = rf_idx >= 4 ? (rf_idx - 4) : rf_idx;
                if ((list_index == 0 &&
                        (ref_pic_index + 1) > pcs_ptr->tpl_data.tpl_ref0_count) ||
                    (list_index == 1 &&
                        (ref_pic_index + 1) > pcs_ptr->tpl_data.tpl_ref1_count))
                    continue;
                if (!is_me_data_valid(pcs_ptr->pa_me_data->me_results[sb_index],
                                        me_mb_offset,
                                        list_index,
                                        ref_pic_index))
                    continue;
                ref_pic_ptr = (EbPictureBufferDesc *)pcs_ptr->tpl_data
                                    .tpl_ref_ds_ptr_array[list_index][ref_pic_index]
                                    .picture_ptr;
                const MeSbResults *me_results = pcs_ptr->pa_me_data->me_results[sb_index];
                x_curr_mv                     = me_results
                                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV +
                                                (list_index ? 4 : 0) + ref_pic_index]
                                .x_mv
                    << 1;
                y_curr_mv = me_results
                                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV +
                                                (list_index ? 4 : 0) + ref_pic_index]
                                .y_mv
                    << 1;
                clip_mv_in_pad(ref_pic_ptr,mb_origin_x,mb_origin_y,&x_curr_mv,&y_curr_mv);
                MV      best_mv          = {y_curr_mv, x_curr_mv};
                if (pcs_ptr->tpl_ctrls.tpl_opt_flag && pcs_ptr->tpl_ctrls.use_pred_sad_in_inter_search) {
                    int32_t ref_origin_index = ref_pic_ptr->origin_x +
                        (mb_origin_x + (best_mv.col >> 3)) +
                        (mb_origin_y + (best_mv.row >> 3) +
                            ref_pic_ptr->origin_y) * ref_pic_ptr->stride_y;
                    //sad_1
                    inter_cost = svt_nxm_sad_kernel_sub_sampled(
                        src_mb,
                        input_ptr->stride_y,
                        ref_pic_ptr->buffer_y + ref_origin_index,
                        ref_pic_ptr->stride_y,
                        16,
                        16);
                }
                else {
                    int32_t ref_origin_index = ref_pic_ptr->origin_x +
                        (mb_origin_x + (best_mv.col >> 3)) +
                        (mb_origin_y + (best_mv.row >> 3) + ref_pic_ptr->origin_y) *
                        ref_pic_ptr->stride_y;

                    svt_aom_subtract_block(16,
                        16,
                        src_diff,
                        16,
                        src_mb,
                        input_picture_ptr->stride_y,
                        ref_pic_ptr->buffer_y + ref_origin_index,
                        ref_pic_ptr->stride_y);
                    EB_TRANS_COEFF_SHAPE pf_shape = pcs_ptr->tpl_ctrls.tpl_opt_flag ? pcs_ptr->tpl_ctrls.pf_shape : DEFAULT_SHAPE;
                    svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size, pf_shape, 8, 0);

                    inter_cost = svt_aom_satd(coeff, 256);
                }
                if (inter_cost < best_inter_cost) {
                    if (!(pcs_ptr->tpl_ctrls.tpl_opt_flag && pcs_ptr->tpl_ctrls.use_pred_sad_in_inter_search))
                    EB_MEMCPY(best_coeff, coeff, sizeof(best_coeff));
                    best_ref_poc = pcs_ptr->tpl_data
                                        .tpl_ref_ds_ptr_array[list_index][ref_pic_index]
                                        .picture_number;

                    best_rf_idx     = rf_idx;
                    best_inter_cost = inter_cost;
                    final_best_mv   = best_mv;

                    if (best_inter_cost < best_intra_cost)
                        best_mode = NEWMV;
                }
            } // rf_idx

            if (best_mode == NEWMV) {
                uint16_t eob = 0;
                if (pcs_ptr->tpl_ctrls.tpl_opt_flag && pcs_ptr->tpl_ctrls.use_pred_sad_in_inter_search) {
                    uint32_t list_index = best_rf_idx < 4 ? 0 : 1;
                    uint32_t ref_pic_index = best_rf_idx >= 4 ? (best_rf_idx - 4) : best_rf_idx;

                    ref_pic_ptr = (EbPictureBufferDesc*)pcs_ptr->tpl_data.tpl_ref_ds_ptr_array[list_index][ref_pic_index].picture_ptr;

                    int32_t ref_origin_index = ref_pic_ptr->origin_x +
                        (mb_origin_x + (final_best_mv.col >> 3)) +
                        (mb_origin_y + (final_best_mv.row >> 3) +
                            ref_pic_ptr->origin_y) * ref_pic_ptr->stride_y;
                    svt_aom_subtract_block(16, 16, src_diff, 16, src_mb, input_picture_ptr->stride_y,
                        ref_pic_ptr->buffer_y + ref_origin_index, ref_pic_ptr->stride_y);
                    EB_TRANS_COEFF_SHAPE pf_shape = pcs_ptr->tpl_ctrls.tpl_opt_flag ? pcs_ptr->tpl_ctrls.pf_shape : DEFAULT_SHAPE;
                    svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size, pf_shape, 8, 0);
                    memcpy(best_coeff, coeff, sizeof(best_coeff));
                }
                get_quantize_error(&mb_plane,
                                    best_coeff,
                                    qcoeff,
                                    dqcoeff,
                                    tx_size,
                                    &eob,
                                    &recon_error,
                                    &sse);
                int rate_cost = pcs_ptr->tpl_ctrls.tpl_opt_flag ? 0 : rate_estimator(qcoeff, eob, tx_size);

                tpl_stats.srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                tpl_stats.srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
            }
            if (best_mode == NEWMV) {
                // inter recon with rec_picture as reference pic
                uint64_t ref_poc       = best_ref_poc;
                uint32_t list_index    = best_rf_idx < 4 ? 0 : 1;
                uint32_t ref_pic_index = best_rf_idx >= 4 ? (best_rf_idx - 4) : best_rf_idx;
                if (pcs_ptr->tpl_data.ref_in_slide_window[list_index][ref_pic_index]) {
                    uint32_t ref_frame_idx = 0;
                    while (ref_frame_idx < MAX_TPL_LA_SW &&
                            encode_context_ptr->poc_map_idx[ref_frame_idx] != ref_poc)
                        ref_frame_idx++;
                    assert(ref_frame_idx != MAX_TPL_LA_SW);
                    ref_pic_ptr =
                        encode_context_ptr->mc_flow_rec_picture_buffer[ref_frame_idx];
                } else
                    ref_pic_ptr = (EbPictureBufferDesc *)pcs_ptr->tpl_data
                                        .tpl_ref_ds_ptr_array[list_index][ref_pic_index]
                                        .picture_ptr;
                int32_t ref_origin_index = ref_pic_ptr->origin_x +
                    (mb_origin_x + (final_best_mv.col >> 3)) +
                    (mb_origin_y + (final_best_mv.row >> 3) + ref_pic_ptr->origin_y) *
                        ref_pic_ptr->stride_y;
                for (int i = 0; i < 16; ++i)
                    EB_MEMCPY(dst_buffer + i * dst_buffer_stride,
                                ref_pic_ptr->buffer_y + ref_origin_index +
                                    i * ref_pic_ptr->stride_y,
                                sizeof(uint8_t) * (16));
            } else {
                // intra recon

                uint8_t *above_row;
                uint8_t *left_col;
                DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
                DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);

                above_row             = above_data + 16;
                left_col              = left_data + 16;
                uint8_t *recon_buffer = recon_picture_ptr->buffer_y + dst_basic_offset;

                update_neighbor_samples_array_open_loop_mb_recon(
                                                                    1, // use_top_righ_bottom_left
                                                                    1, // update_top_neighbor
                                                                    above_row - 1,
                                                                    left_col - 1,
                                                                    recon_buffer,
                                                                    dst_buffer_stride,
                                                                    mb_origin_x,
                                                                    mb_origin_y,
                                                                    16,
                                                                    16,
                                                                    input_picture_ptr->width,
                                                                    input_picture_ptr->height);

                uint8_t ois_intra_mode = best_intra_mode; // ois_mb_results_ptr->intra_mode;
                int32_t p_angle = av1_is_directional_mode((PredictionMode)ois_intra_mode)
                    ? mode_to_angle_map[(PredictionMode)ois_intra_mode]
                    : 0;
                // Edge filter
                if (av1_is_directional_mode((PredictionMode)ois_intra_mode) &&
                    1 /*scs_ptr->seq_header.enable_intra_edge_filter*/) {
                    filter_intra_edge(NULL,
                                        ois_intra_mode,
                                        scs_ptr->seq_header.max_frame_width,
                                        scs_ptr->seq_header.max_frame_height,
                                        p_angle,
                                        mb_origin_x,
                                        mb_origin_y,
                                        above_row,
                                        left_col);
                }
                // PRED
                intra_prediction_open_loop_mb(p_angle,
                                                ois_intra_mode,
                                                mb_origin_x,
                                                mb_origin_y,
                                                tx_size,
                                                above_row,
                                                left_col,
                                                dst_buffer,
                                                dst_buffer_stride);
            }

            svt_aom_subtract_block(16,
                                    16,
                                    src_diff,
                                    16,
                                    src_mb,
                                    input_picture_ptr->stride_y,
                                    dst_buffer,
                                    dst_buffer_stride);
            EB_TRANS_COEFF_SHAPE pf_shape = pcs_ptr->tpl_ctrls.tpl_opt_flag ? pcs_ptr->tpl_ctrls.pf_shape : DEFAULT_SHAPE;
            svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size,pf_shape, 8, 0);

            uint16_t eob = 0;

            get_quantize_error(
                &mb_plane, coeff, qcoeff, dqcoeff, tx_size, &eob, &recon_error, &sse);
            int rate_cost = pcs_ptr->tpl_ctrls.tpl_opt_flag ? 0 : rate_estimator(qcoeff, eob, tx_size);
            // Disable intra prediction
            disable_intra_pred  = tpl_opt_flag && (pcs_ptr->tpl_ctrls.disable_intra_pred_nref ||
                pcs_ptr->tpl_ctrls.disable_intra_pred_nbase);
            if (!disable_intra_pred || (pcs_ptr->tpl_data.is_used_as_reference_flag))
                if (eob) {
                    av1_inv_transform_recon8bit((int32_t *)dqcoeff,
                                                dst_buffer,
                                                dst_buffer_stride,
                                                dst_buffer,
                                                dst_buffer_stride,
                                                TX_16X16,
                                                DCT_DCT,
                                                PLANE_TYPE_Y,
                                                eob,
                                                0);
                }

            tpl_stats.recrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
            tpl_stats.recrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
            if (best_mode != NEWMV) {
                tpl_stats.srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
                tpl_stats.srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
            }
            tpl_stats.recrf_dist = AOMMAX(tpl_stats.srcrf_dist, tpl_stats.recrf_dist);
            tpl_stats.recrf_rate = AOMMAX(tpl_stats.srcrf_rate, tpl_stats.recrf_rate);
            if (pcs_ptr->tpl_data.tpl_slice_type != I_SLICE && best_rf_idx != -1) {
                tpl_stats.mv            = final_best_mv;
                tpl_stats.ref_frame_poc = best_ref_poc;
            }
            // Motion flow dependency dispenser.
            result_model_store(pcs_ptr, &tpl_stats, mb_origin_x, mb_origin_y);
        }
        pa_blk_index++;
    }

    }

}

#define TPL_TASKS_MDC_INPUT 0
#define TPL_TASKS_ENCDEC_INPUT 1
#define TPL_TASKS_CONTINUE 2
/*
   Assign TPL dispenser segments
*/
EbBool assign_tpl_segments(EncDecSegments *segmentPtr, uint16_t *segmentInOutIndex,
                               TplDispResults * taskPtr,
    int32_t                          frame_idx, EbFifo *srmFifoPtr) {
    EbBool           continue_processing_flag = EB_FALSE;
    uint32_t row_segment_index = 0;
    uint32_t segment_index;
    uint32_t right_segment_index;
    uint32_t bottom_left_segment_index;

    int16_t feedback_row_index = -1;

    uint32_t self_assigned = EB_FALSE;

    //static FILE *trace = 0;
    //
    //if(trace == 0) {
    //    trace = fopen("seg-trace.txt","w");
    //}

    switch (taskPtr->input_type) {
    case TPL_TASKS_MDC_INPUT:

        // The entire picture is provided by the MDC process, so
        //   no logic is necessary to clear input dependencies.
        for (uint32_t row_index = 0; row_index < segmentPtr->segment_row_count; ++row_index) {
            segmentPtr->row_array[row_index].current_seg_index =
                segmentPtr->row_array[row_index].starting_seg_index;
        }


        // Start on Segment 0 immediately
        *segmentInOutIndex  = segmentPtr->row_array[0].current_seg_index;
        taskPtr->input_type = TPL_TASKS_CONTINUE;
        ++segmentPtr->row_array[0].current_seg_index;
        continue_processing_flag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case TPL_TASKS_ENCDEC_INPUT:

        // Setup row_segment_index to release the in_progress token
        //row_segment_index = taskPtr->encDecSegmentRowArray[0];

        // Start on the assigned row immediately
        *segmentInOutIndex  = segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        taskPtr->input_type = TPL_TASKS_CONTINUE;
        ++segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        continue_processing_flag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case TPL_TASKS_CONTINUE:

        // Update the Dependency List for Right and Bottom Neighbors
        segment_index     = *segmentInOutIndex;
        row_segment_index = segment_index / segmentPtr->segment_band_count;

        right_segment_index       = segment_index + 1;
        bottom_left_segment_index = segment_index + segmentPtr->segment_band_count;

        // Right Neighbor
        if (segment_index < segmentPtr->row_array[row_segment_index].ending_seg_index) {
            svt_block_on_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[right_segment_index];

            if (segmentPtr->dep_map.dependency_map[right_segment_index] == 0) {
                *segmentInOutIndex = segmentPtr->row_array[row_segment_index].current_seg_index;
                ++segmentPtr->row_array[row_segment_index].current_seg_index;
                self_assigned            = EB_TRUE;
                continue_processing_flag = EB_TRUE;

                //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
                //    *segmentInOutIndex);
            }

            svt_release_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);
        }

        // Bottom-left Neighbor
        if (row_segment_index < segmentPtr->segment_row_count - 1 &&
            bottom_left_segment_index >=
                segmentPtr->row_array[row_segment_index + 1].starting_seg_index) {
            svt_block_on_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[bottom_left_segment_index];

            if (segmentPtr->dep_map.dependency_map[bottom_left_segment_index] == 0) {
                if (self_assigned == EB_TRUE)
                    feedback_row_index = (int16_t)row_segment_index + 1;
                else {
                    *segmentInOutIndex =
                        segmentPtr->row_array[row_segment_index + 1].current_seg_index;
                    ++segmentPtr->row_array[row_segment_index + 1].current_seg_index;
                    continue_processing_flag = EB_TRUE;

                    //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                    //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
                    //    *segmentInOutIndex);
                }
            }
            svt_release_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);
        }

        if (feedback_row_index > 0) {

            EbObjectWrapper *out_results_wrapper_ptr;

            svt_get_empty_object(
                    srmFifoPtr ,
                    &out_results_wrapper_ptr);

            TplDispResults *out_results_ptr = (TplDispResults*)out_results_wrapper_ptr->object_ptr;
            out_results_ptr->input_type          = TPL_TASKS_ENCDEC_INPUT;

            out_results_ptr->enc_dec_segment_row = feedback_row_index;
            out_results_ptr->tile_group_index = taskPtr->tile_group_index;
            out_results_ptr->qIndex = taskPtr->qIndex;

            out_results_ptr->pcs_wrapper_ptr = taskPtr->pcs_wrapper_ptr;
            out_results_ptr->pcs_ptr = taskPtr->pcs_ptr;
            out_results_ptr->frame_index = frame_idx;
            svt_post_full_object(out_results_wrapper_ptr);
        }

        break;

    default: break;
    }

    return continue_processing_flag;
}





/************************************************
* Genrate TPL MC Flow Dispenser  Based on Lookahead
** LAD Window: sliding window size
************************************************/


void tpl_mc_flow_dispenser(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *scs_ptr,
    int32_t                         *base_rdmult,
    PictureParentControlSet        *pcs_ptr,
    int32_t                          frame_idx,
    SourceBasedOperationsContext    *context_ptr)
{
    EbPictureBufferDesc *recon_picture_ptr = encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx];


    int32_t         qIndex = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
    if (pcs_ptr->tpl_ctrls.enable_tpl_qps){
        const double delta_rate_new[7][6] = {
            {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, // 1L
            {0.6, 1.0, 1.0, 1.0, 1.0, 1.0}, // 2L
            {0.6, 0.8, 1.0, 1.0, 1.0, 1.0}, // 3L
            {0.6, 0.8, 0.9, 1.0, 1.0, 1.0}, // 4L
            {0.35, 0.6, 0.8, 0.9, 1.0, 1.0}, //5L
            {0.35, 0.6, 0.8, 0.9, 0.95, 1.0} //6L
        };
        double q_val;
        q_val = svt_av1_convert_qindex_to_q(qIndex, 8);
        int32_t delta_qindex;
        if (pcs_ptr->tpl_data.tpl_slice_type == I_SLICE)
            delta_qindex = svt_av1_compute_qdelta(q_val, q_val * 0.25, 8);
        else
            delta_qindex = svt_av1_compute_qdelta(
                q_val,
                q_val *
                    delta_rate_new[pcs_ptr->hierarchical_levels]
                                  [pcs_ptr->tpl_data.tpl_temporal_layer_index],
                8);
        qIndex = (qIndex + delta_qindex);
    }
    *base_rdmult = svt_av1_compute_rd_mult_based_on_qindex((AomBitDepth)8/*scs_ptr->static_config.encoder_bit_depth*/, qIndex) / 6;

    {
        {


        // reset number of TPLed sbs per pic
        pcs_ptr->tpl_disp_coded_sb_count = 0;

        EbObjectWrapper *out_results_wrapper_ptr;

        // TPL dispenser kernel
        svt_get_empty_object(
                context_ptr->sbo_output_fifo_ptr,
                &out_results_wrapper_ptr);

        TplDispResults *out_results_ptr = (TplDispResults*)out_results_wrapper_ptr->object_ptr;
       // out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
        out_results_ptr->pcs_ptr = pcs_ptr;
        out_results_ptr->input_type       = TPL_TASKS_MDC_INPUT;
        out_results_ptr->tile_group_index = /*tile_group_idx*/0;

        out_results_ptr->frame_index = frame_idx;
        out_results_ptr->qIndex = qIndex;

        svt_post_full_object(out_results_wrapper_ptr);

        svt_block_on_semaphore(pcs_ptr->tpl_disp_done_semaphore); // we can do all in // ?


        }
    }

    // padding current recon picture
    generate_padding(recon_picture_ptr->buffer_y,
                     recon_picture_ptr->stride_y,
                     recon_picture_ptr->width,
                     recon_picture_ptr->height,
                     recon_picture_ptr->origin_x,
                     recon_picture_ptr->origin_y);

    return;
}


static int get_overlap_area(int grid_pos_row, int grid_pos_col, int ref_pos_row, int ref_pos_col,
                            int block, int /*BLOCK_SIZE*/ bsize) {
    int width = 0, height = 0;
    int bw = 4 << mi_size_wide_log2[bsize];
    int bh = 4 << mi_size_high_log2[bsize];

    switch (block) {
    case 0:
        width  = grid_pos_col + bw - ref_pos_col;
        height = grid_pos_row + bh - ref_pos_row;
        break;
    case 1:
        width  = ref_pos_col + bw - grid_pos_col;
        height = grid_pos_row + bh - ref_pos_row;
        break;
    case 2:
        width  = grid_pos_col + bw - ref_pos_col;
        height = ref_pos_row + bh - grid_pos_row;
        break;
    case 3:
        width  = ref_pos_col + bw - grid_pos_col;
        height = ref_pos_row + bh - grid_pos_row;
        break;
    default: assert(0);
    }

    return width * height;
}

static int round_floor(int ref_pos, int bsize_pix) {
    int round;
    if (ref_pos < 0)
        round = -(1 + (-ref_pos - 1) / bsize_pix);
    else
        round = ref_pos / bsize_pix;

    return round;
}

static int64_t delta_rate_cost(int64_t delta_rate, int64_t recrf_dist, int64_t srcrf_dist,
                               int pix_num) {
    double  beta      = (double)srcrf_dist / recrf_dist;
    int64_t rate_cost = delta_rate;

    if (srcrf_dist <= 128)
        return rate_cost;

    double dr = (double)(delta_rate >> (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT)) / pix_num;

    double log_den = log(beta) / log(2.0) + 2.0 * dr;

    if (log_den > log(10.0) / log(2.0)) {
        rate_cost = (int64_t)((log(1.0 / beta) * pix_num) / log(2.0) / 2.0);
        rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);
        return rate_cost;
    }

    double num = pow(2.0, log_den);
    double den = num * beta + (1 - beta) * beta;

    rate_cost = (int64_t)((pix_num * log(num / den)) / log(2.0) / 2.0);

    rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);

    return rate_cost;
}
/************************************************
* Genrate TPL MC Flow Synthesizer
************************************************/


static AOM_INLINE void tpl_model_update_b(PictureParentControlSet *ref_pcs_ptr, PictureParentControlSet *pcs_ptr,
    TplStats *tpl_stats_ptr,
    int mi_row, int mi_col,
    const int/*BLOCK_SIZE*/ bsize) {
    Av1Common *ref_cm = ref_pcs_ptr->av1_cm;
    TplStats * ref_tpl_stats_ptr;

    const FULLPEL_MV full_mv     = get_fullmv_from_mv(&tpl_stats_ptr->mv);
    const int        ref_pos_row = mi_row * MI_SIZE + full_mv.row;
    const int        ref_pos_col = mi_col * MI_SIZE + full_mv.col;

    const int bw         = 4 << mi_size_wide_log2[bsize];
    const int bh         = 4 << mi_size_high_log2[bsize];
    const int mi_height  = mi_size_high[bsize];
    const int mi_width   = mi_size_wide[bsize];
    const int pix_num    = bw * bh;
    const int shift      = pcs_ptr->is_720p_or_larger ? 2 : 1;
    const int mi_cols_sr = ((ref_pcs_ptr->aligned_width + 15) / 16) << 2;

    // top-left on grid block location in pixel
    int grid_pos_row_base = round_floor(ref_pos_row, bh) * bh;
    int grid_pos_col_base = round_floor(ref_pos_col, bw) * bw;
    int block;

    int64_t cur_dep_dist = tpl_stats_ptr->recrf_dist - tpl_stats_ptr->srcrf_dist;
    int64_t mc_dep_dist  = (int64_t)(
        tpl_stats_ptr->mc_dep_dist *
        ((double)(tpl_stats_ptr->recrf_dist - tpl_stats_ptr->srcrf_dist) /
         tpl_stats_ptr->recrf_dist));
    int64_t delta_rate  = tpl_stats_ptr->recrf_rate - tpl_stats_ptr->srcrf_rate;
    int64_t mc_dep_rate = pcs_ptr->tpl_ctrls.tpl_opt_flag ? 0

        : delta_rate_cost(tpl_stats_ptr->mc_dep_rate,
                          tpl_stats_ptr->recrf_dist,
                          tpl_stats_ptr->srcrf_dist,
                          pix_num);

    for (block = 0; block < 4; ++block) {
        int grid_pos_row = grid_pos_row_base + bh * (block >> 1);
        int grid_pos_col = grid_pos_col_base + bw * (block & 0x01);

        if (grid_pos_row >= 0 && grid_pos_row < ref_cm->mi_rows * MI_SIZE && grid_pos_col >= 0 &&
            grid_pos_col < ref_cm->mi_cols * MI_SIZE) {
            int overlap_area = get_overlap_area(
                grid_pos_row, grid_pos_col, ref_pos_row, ref_pos_col, block, bsize);
            int       ref_mi_row = round_floor(grid_pos_row, bh) * mi_height;
            int       ref_mi_col = round_floor(grid_pos_col, bw) * mi_width;
            const int step       = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);

            for (int idy = 0; idy < mi_height; idy += step) {
                for (int idx = 0; idx < mi_width; idx += step) {
                    ref_tpl_stats_ptr = ref_pcs_ptr->tpl_stats[((ref_mi_row + idy) >> shift) *
                                                                   (mi_cols_sr >> shift) +
                                                               ((ref_mi_col + idx) >> shift)];
                    ref_tpl_stats_ptr->mc_dep_dist += ((cur_dep_dist + mc_dep_dist) *
                                                       overlap_area) /
                        pix_num;
                    ref_tpl_stats_ptr->mc_dep_rate += ((delta_rate + mc_dep_rate) * overlap_area) /
                        pix_num;
                    assert(overlap_area >= 0);
                }
            }
        }
    }
}

/************************************************
* Genrate TPL MC Flow Synthesizer
************************************************/


static AOM_INLINE void tpl_model_update(
    PictureParentControlSet     *pcs_array[MAX_TPL_LA_SW],
    int32_t frame_idx, int mi_row, int mi_col,
    const int/*BLOCK_SIZE*/ bsize, uint8_t frames_in_sw) {
    const int                mi_height  = mi_size_high[bsize];
    const int                mi_width   = mi_size_wide[bsize];
    PictureParentControlSet  *pcs_ptr = pcs_array[frame_idx];
    const int /*BLOCK_SIZE*/ block_size = pcs_ptr->is_720p_or_larger ? BLOCK_16X16 : BLOCK_8X8;
    const int                step       = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int                shift      = pcs_ptr->is_720p_or_larger ? 2 : 1;
    const int                mi_cols_sr = ((pcs_ptr->aligned_width + 15) / 16) << 2;
    int                      i          = 0;

    for (int idy = 0; idy < mi_height; idy += step) {
        for (int idx = 0; idx < mi_width; idx += step) {
            TplStats *tpl_stats_ptr =
                pcs_ptr->tpl_stats[(((mi_row + idy) >> shift) * (mi_cols_sr >> shift)) +
                                   ((mi_col + idx) >> shift)];

            while (i < frames_in_sw && pcs_array[i]->picture_number != tpl_stats_ptr->ref_frame_poc)
                i++;
            if (i < frames_in_sw)
                tpl_model_update_b(
                    pcs_array[i], pcs_ptr, tpl_stats_ptr, mi_row + idy, mi_col + idx, block_size);
        }
    }
}



/************************************************
* Genrate TPL MC Flow Synthesizer Based on Lookahead
** LAD Window: sliding window size
************************************************/


void tpl_mc_flow_synthesizer(
    PictureParentControlSet         *pcs_array[MAX_TPL_LA_SW],
    int32_t                          frame_idx,
    uint8_t                          frames_in_sw)
{
    Av1Common *              cm        = pcs_array[frame_idx]->av1_cm;
    const int /*BLOCK_SIZE*/ bsize     = BLOCK_16X16;
    const int                mi_height = mi_size_high[bsize];
    const int                mi_width  = mi_size_wide[bsize];

    for (int mi_row = 0; mi_row < cm->mi_rows; mi_row += mi_height) {
        for (int mi_col = 0; mi_col < cm->mi_cols; mi_col += mi_width) {
            tpl_model_update(pcs_array, frame_idx, mi_row, mi_col, bsize, frames_in_sw);
        }
    }
    return;
}


static void generate_r0beta(PictureParentControlSet *pcs_ptr) {
    Av1Common *         cm               = pcs_ptr->av1_cm;
    SequenceControlSet *scs_ptr          = pcs_ptr->scs_ptr;
    int64_t             intra_cost_base  = 0;
    int64_t             mc_dep_cost_base = 0;
    const int           step             = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int           mi_cols_sr       = ((pcs_ptr->aligned_width + 15) / 16) << 2;
    const int           shift            = pcs_ptr->is_720p_or_larger ? 2 : 1;

    for (int row = 0; row < cm->mi_rows; row += step) {
        for (int col = 0; col < mi_cols_sr; col += step) {
            TplStats *tpl_stats_ptr =
                pcs_ptr->tpl_stats[(row >> shift) * (mi_cols_sr >> shift) + (col >> shift)];
            int64_t mc_dep_delta = RDCOST(
                pcs_ptr->base_rdmult, tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->mc_dep_dist);
            intra_cost_base += (tpl_stats_ptr->recrf_dist << RDDIV_BITS);
            mc_dep_cost_base += (tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
        }
    }

    if (mc_dep_cost_base != 0) {
        pcs_ptr->r0 = (double)intra_cost_base / mc_dep_cost_base;
        pcs_ptr->tpl_is_valid = 1;
    }
    else {
        pcs_ptr->tpl_is_valid = 0;
    }

#if DEBUG_TPL
    SVT_LOG("generate_r0beta ------> poc %ld\t%.0f\t%.0f \t%.5f base_rdmult=%d\n",
            pcs_ptr->picture_number,
            (double)intra_cost_base,
            (double)mc_dep_cost_base,
            pcs_ptr->r0,
            pcs_ptr->base_rdmult);
#endif
    generate_lambda_scaling_factor(pcs_ptr, mc_dep_cost_base);

    const uint32_t sb_sz            = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 128 : 64;
    const uint32_t picture_sb_width = (uint32_t)((scs_ptr->seq_header.max_frame_width + sb_sz - 1) /
                                                 sb_sz);
    const uint32_t picture_sb_height = (uint32_t)(
        (scs_ptr->seq_header.max_frame_height + sb_sz - 1) / sb_sz);
    const uint32_t picture_width_in_mb  = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t picture_height_in_mb = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    const uint32_t blks                 = scs_ptr->seq_header.sb_size == BLOCK_128X128
                        ? (128 >> (3 + pcs_ptr->is_720p_or_larger))
                        : (64 >> (3 + pcs_ptr->is_720p_or_larger));
    for (uint32_t sb_y = 0; sb_y < picture_sb_height; ++sb_y) {
        for (uint32_t sb_x = 0; sb_x < picture_sb_width; ++sb_x) {
            int64_t intra_cost  = 0;
            int64_t mc_dep_cost = 0;
            for (uint32_t blky_offset = 0; blky_offset < blks; blky_offset++) {
                for (uint32_t blkx_offset = 0; blkx_offset < blks; blkx_offset++) {
                    uint32_t blkx = ((sb_x * sb_sz) >> (3 + pcs_ptr->is_720p_or_larger)) +
                        blkx_offset;
                    uint32_t blky = ((sb_y * sb_sz) >> (3 + pcs_ptr->is_720p_or_larger)) +
                        blky_offset;
                    if ((blkx >> (1 - pcs_ptr->is_720p_or_larger)) >= picture_width_in_mb ||
                        (blky >> (1 - pcs_ptr->is_720p_or_larger)) >= picture_height_in_mb)
                        continue;
                    TplStats *tpl_stats_ptr =
                        pcs_ptr->tpl_stats[blky * (mi_cols_sr >> shift) + blkx];
                    int64_t mc_dep_delta = RDCOST(pcs_ptr->base_rdmult,
                                                  tpl_stats_ptr->mc_dep_rate,
                                                  tpl_stats_ptr->mc_dep_dist);
                    intra_cost += (tpl_stats_ptr->recrf_dist << RDDIV_BITS);
                    mc_dep_cost += (tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
                }
            }
            double beta = 1.0;
            if (mc_dep_cost > 0 && intra_cost > 0) {
                double rk = (double)intra_cost / mc_dep_cost;
                beta      = (pcs_ptr->r0 / rk);
                assert(beta > 0.0);
            }
            pcs_ptr->tpl_beta[sb_y * picture_sb_width + sb_x] = beta;
        }
    }
    return;
}
/************************************************
* Allocate and initialize buffers needed for tpl
************************************************/
EbErrorType init_tpl_buffers(
    EncodeContext                   *encode_context_ptr,
    PictureParentControlSet         *pcs_ptr){
    int32_t frames_in_sw = MIN(MAX_TPL_LA_SW, pcs_ptr->tpl_group_size);
    int32_t frame_idx;

    for (frame_idx = 0; frame_idx < MAX_TPL_LA_SW; frame_idx++) {
        encode_context_ptr->poc_map_idx[frame_idx]                = -1;
        encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] = NULL;
    }
    EbPictureBufferDescInitData picture_buffer_desc_init_data;
    picture_buffer_desc_init_data.max_width          = pcs_ptr->enhanced_picture_ptr->max_width;
    picture_buffer_desc_init_data.max_height         = pcs_ptr->enhanced_picture_ptr->max_height;
    picture_buffer_desc_init_data.bit_depth          = pcs_ptr->enhanced_picture_ptr->bit_depth;
    picture_buffer_desc_init_data.color_format       = pcs_ptr->enhanced_picture_ptr->color_format;
    picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
    picture_buffer_desc_init_data.left_padding       = TPL_PADX;
    picture_buffer_desc_init_data.right_padding      = TPL_PADX;
    picture_buffer_desc_init_data.top_padding        = TPL_PADY;
    picture_buffer_desc_init_data.bot_padding        = TPL_PADY;
    picture_buffer_desc_init_data.split_mode         = EB_FALSE;

    EB_NEW(encode_context_ptr->mc_flow_rec_picture_buffer_noref,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&picture_buffer_desc_init_data);

    for (frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
        if (pcs_ptr->tpl_valid_pic[frame_idx]) {
            EB_NEW(encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx],
                   svt_picture_buffer_desc_ctor,
                   (EbPtr)&picture_buffer_desc_init_data);
        } else {
            encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] =
                encode_context_ptr->mc_flow_rec_picture_buffer_noref;
        }
    }
    return EB_ErrorNone;
}




/************************************************
* init tpl tpl_disp_segment_ctrl
************************************************/
void init_tpl_segments(
    SequenceControlSet              *scs_ptr,
    PictureParentControlSet         *pcs_ptr,
    PictureParentControlSet        **pcs_array,
    int32_t                         frames_in_sw) {

    for (int32_t frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
        uint32_t enc_dec_seg_col_cnt = scs_ptr->tpl_segment_col_count_array;
        uint32_t enc_dec_seg_row_cnt = scs_ptr->tpl_segment_row_count_array;

        const int tile_cols = pcs_ptr->av1_cm->tiles_info.tile_cols;
        const int tile_rows = pcs_ptr->av1_cm->tiles_info.tile_rows;
        uint8_t   tile_group_cols = MIN(
            tile_cols,
            scs_ptr->tile_group_col_count_array[pcs_ptr->temporal_layer_index]);
        uint8_t tile_group_rows = MIN(
            tile_rows,
            scs_ptr->tile_group_row_count_array[pcs_ptr->temporal_layer_index]);

        // Valid when only one tile used
        // TPL segments + tiles (not working)
        // TPL segments are 64x64 SB based
        uint16_t                        pic_width_in_sb;
        uint16_t                        pic_height_in_sb;
        pic_width_in_sb = (pcs_ptr->aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
        pic_height_in_sb   = (pcs_ptr->aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;

        if (tile_group_cols * tile_group_rows > 1) {
            enc_dec_seg_col_cnt = MIN(enc_dec_seg_col_cnt,
                (uint8_t)(pic_width_in_sb / tile_group_cols));
            enc_dec_seg_row_cnt = MIN(
                enc_dec_seg_row_cnt,
                (uint8_t)(pic_height_in_sb / tile_group_rows));
        }
        // Init segments within the tile group
        int      sb_size_log2 = scs_ptr->seq_header.sb_size_log2;

        uint8_t tile_group_col_start_tile_idx[1024];
        uint8_t tile_group_row_start_tile_idx[1024];

        // Get the tile start index for tile group
        for (uint8_t c = 0; c <= tile_group_cols; c++) {
            tile_group_col_start_tile_idx[c] = c * tile_cols / tile_group_cols;
        }
        for (uint8_t r = 0; r <= tile_group_rows; r++) {
            tile_group_row_start_tile_idx[r] = r * tile_rows / tile_group_rows;
        }

        for (uint8_t r = 0; r < tile_group_rows; r++) {
            for (uint8_t c = 0; c < tile_group_cols; c++) {
                uint16_t tile_group_idx = r * tile_group_cols + c;
                uint16_t top_left_tile_col_idx = tile_group_col_start_tile_idx[c];
                uint16_t top_left_tile_row_idx = tile_group_row_start_tile_idx[r];
                uint16_t bottom_right_tile_col_idx =
                    tile_group_col_start_tile_idx[c + 1];
                uint16_t bottom_right_tile_row_idx =
                    tile_group_row_start_tile_idx[r + 1];

                TileGroupInfo *tg_info_ptr =
                    &pcs_array[frame_idx]->tile_group_info[tile_group_idx];

                tg_info_ptr->tile_group_tile_start_x = top_left_tile_col_idx;
                tg_info_ptr->tile_group_tile_end_x = bottom_right_tile_col_idx;

                tg_info_ptr->tile_group_tile_start_y = top_left_tile_row_idx;
                tg_info_ptr->tile_group_tile_end_y = bottom_right_tile_row_idx;

                tg_info_ptr->tile_group_sb_start_x =
                    pcs_ptr->av1_cm->tiles_info.tile_col_start_mi[top_left_tile_col_idx] >>
                    sb_size_log2;
                tg_info_ptr->tile_group_sb_start_y =
                    pcs_ptr->av1_cm->tiles_info.tile_row_start_mi[top_left_tile_row_idx] >>
                    sb_size_log2;




                // Get the SB end of the bottom right tile
                tg_info_ptr->tile_group_sb_end_x = pic_width_in_sb ;
                    //(pcs_ptr->av1_cm->tiles_info.tile_col_start_mi[bottom_right_tile_col_idx] >>
                    //    sb_size_log2);
                tg_info_ptr->tile_group_sb_end_y = pic_height_in_sb;
                    //(pcs_ptr->av1_cm->tiles_info.tile_row_start_mi[bottom_right_tile_row_idx] >>
                    //    sb_size_log2);

                // Get the width/height of tile group in SB
                tg_info_ptr->tile_group_height_in_sb =
                    tg_info_ptr->tile_group_sb_end_y -
                    tg_info_ptr->tile_group_sb_start_y;
                tg_info_ptr->tile_group_width_in_sb =
                    tg_info_ptr->tile_group_sb_end_x -
                    tg_info_ptr->tile_group_sb_start_x;

                enc_dec_segments_init(
                    pcs_array[frame_idx]->tpl_disp_segment_ctrl[tile_group_idx],
                    enc_dec_seg_col_cnt,
                    enc_dec_seg_row_cnt,
                    tg_info_ptr->tile_group_width_in_sb,
                    tg_info_ptr->tile_group_height_in_sb);
            }
        }
    }
}


/************************************************
* Genrate TPL MC Flow Based on frames in the tpl group
************************************************/
EbErrorType tpl_mc_flow(EncodeContext *encode_context_ptr, SequenceControlSet *scs_ptr,
                        PictureParentControlSet *pcs_ptr,  SourceBasedOperationsContext    *context_ptr) {

    int32_t  frames_in_sw = MIN(MAX_TPL_LA_SW, pcs_ptr->tpl_group_size);
    int32_t  frame_idx;
    uint32_t shift                = pcs_ptr->is_720p_or_larger ? 0 : 1;
    uint32_t picture_width_in_mb  = (pcs_ptr->enhanced_picture_ptr->width + 16 - 1) / 16;
    uint32_t picture_height_in_mb = (pcs_ptr->enhanced_picture_ptr->height + 16 - 1) / 16;

    //wait for PA ME to be done.
    for (uint32_t i = 1; i < pcs_ptr->tpl_group_size; i++) {
        svt_wait_cond_var(&pcs_ptr->tpl_group[i]->me_ready, 0);
    }
    pcs_ptr->tpl_is_valid = 0;
    init_tpl_buffers(encode_context_ptr, pcs_ptr);

    if (pcs_ptr->tpl_group[0]->tpl_data.tpl_temporal_layer_index == 0) {


        // no Tiles path
        if (scs_ptr->static_config.tile_rows == 0 && scs_ptr->static_config.tile_columns == 0 )
            init_tpl_segments(
                scs_ptr,
                pcs_ptr,
                pcs_ptr->tpl_group,
                frames_in_sw) ;



        uint8_t tpl_on;
        encode_context_ptr->poc_map_idx[0] = pcs_ptr->tpl_group[0]->picture_number;
        for (frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
            encode_context_ptr->poc_map_idx[frame_idx] = pcs_ptr->tpl_group[frame_idx]->picture_number;
            for (uint32_t blky = 0; blky < (picture_height_in_mb << shift); blky++) {
                memset(pcs_ptr->tpl_group[frame_idx]->tpl_stats[blky * (picture_width_in_mb << shift)],
                        0,
                        (picture_width_in_mb << shift) * sizeof(TplStats));
            }
            if(scs_ptr->lad_mg)
                tpl_on = pcs_ptr->tpl_valid_pic[frame_idx];
            else {
                tpl_on = !(pcs_ptr->tpl_group[0]->tpl_ctrls.disable_tpl_nref);
                tpl_on = (pcs_ptr->tpl_group[0]->slice_type == I_SLICE) ? 1 : tpl_on;
                if (tpl_on == 0) {
                    tpl_on = pcs_ptr->tpl_group[frame_idx]->tpl_data.is_used_as_reference_flag ? 1 :
                        (ABS((int64_t)pcs_ptr->tpl_group[0]->picture_number -
                        (int64_t)pcs_ptr->tpl_group[frame_idx]->picture_number)
                        <= pcs_ptr->tpl_group[0]->tpl_ctrls.disable_tpl_pic_dist) ? 1 : tpl_on;
                }
            }
            if (tpl_on)
                tpl_mc_flow_dispenser(encode_context_ptr, scs_ptr, &pcs_ptr->base_rdmult, pcs_ptr->tpl_group[frame_idx], frame_idx,context_ptr);
        }

        // synthesizer
        for (frame_idx = frames_in_sw - 1; frame_idx >= 0; frame_idx--) {
            if(scs_ptr->lad_mg)
                tpl_on = pcs_ptr->tpl_valid_pic[frame_idx];
            else {
                tpl_on = !(pcs_ptr->tpl_group[0]->tpl_ctrls.disable_tpl_nref);
                tpl_on = (pcs_ptr->tpl_group[0]->slice_type == I_SLICE) ? 1 : tpl_on;
                if (tpl_on == 0) {
                    tpl_on = pcs_ptr->tpl_group[frame_idx]->tpl_data.is_used_as_reference_flag ? 1 :
                        (ABS((int64_t)pcs_ptr->tpl_group[0]->picture_number -
                        (int64_t)pcs_ptr->tpl_group[frame_idx]->picture_number)
                            <= pcs_ptr->tpl_group[0]->tpl_ctrls.disable_tpl_pic_dist) ? 1 : tpl_on;
                }
            }
            if (tpl_on)
                tpl_mc_flow_synthesizer(pcs_ptr->tpl_group, frame_idx, frames_in_sw);
        }

        // generate tpl stats
        generate_r0beta(pcs_ptr);
#if DEBUG_TPL
        SVT_LOG("LOG displayorder:%ld\n",
            pcs_array[0]->picture_number);
        for (frame_idx = 0; frame_idx < frames_in_sw; frame_idx++)
        {
            PictureParentControlSet *pcs_ptr_tmp = pcs_array[frame_idx];
            Av1Common *cm = pcs_ptr->av1_cm;
            SequenceControlSet *scs_ptr = pcs_ptr_tmp->scs_ptr;
            int64_t intra_cost_base = 0;
            int64_t mc_dep_cost_base = 0;
            const int step = 1 << (pcs_ptr_tmp->is_720p_or_larger ? 2 : 1);
            const int mi_cols_sr = ((pcs_ptr_tmp->aligned_width + 15) / 16) << 2;
            const int shift = pcs_ptr_tmp->is_720p_or_larger ? 2 : 1;

            for (int row = 0; row < cm->mi_rows; row += step) {
                for (int col = 0; col < mi_cols_sr; col += step) {
                    TplStats *tpl_stats_ptr = pcs_ptr_tmp->tpl_stats[(row >> shift) * (mi_cols_sr >> shift) + (col >> shift)];
                    int64_t mc_dep_delta =
                        RDCOST(pcs_ptr->base_rdmult, tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->mc_dep_dist);
                    intra_cost_base += (tpl_stats_ptr->recrf_dist << RDDIV_BITS);
                    mc_dep_cost_base += (tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
                }
            }

            SVT_LOG("After mc_flow_synthesizer:\tframe_indx:%d\tdisplayorder:%ld\tIntra:%lld\tmc_dep:%lld rdmult:%i\n",
                frame_idx, pcs_ptr_tmp->picture_number, intra_cost_base, mc_dep_cost_base, pcs_ptr->base_rdmult);
        }
#endif



    }

    for (frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
        if (encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] &&
            encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] !=
                encode_context_ptr->mc_flow_rec_picture_buffer_noref)
            EB_DELETE(encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx]);
    }
    EB_DELETE(encode_context_ptr->mc_flow_rec_picture_buffer_noref);

    for (uint32_t i = 0; i < pcs_ptr->tpl_group_size; i++) {
        if (pcs_ptr->tpl_group[i]->slice_type == P_SLICE) {
            if (pcs_ptr->tpl_group[i]->ext_mg_id == pcs_ptr->ext_mg_id + 1)
                release_pa_reference_objects(scs_ptr, pcs_ptr->tpl_group[i]);
        }
        else {
            if (pcs_ptr->tpl_group[i]->ext_mg_id == pcs_ptr->ext_mg_id)
                release_pa_reference_objects(scs_ptr, pcs_ptr->tpl_group[i]);
        }
        if (pcs_ptr->tpl_group[i]->non_tf_input)
            EB_DELETE(pcs_ptr->tpl_group[i]->non_tf_input);
    }

    return EB_ErrorNone;
}


/*
   TPL dispenser kernel
   process one picture of TPL group
*/


void *tpl_disp_kernel(void *input_ptr) {
    EbThreadContext *             thread_context_ptr = (EbThreadContext *)input_ptr;
    TplDispenserContext *context_ptr =
        (TplDispenserContext *)thread_context_ptr->priv;
    EbObjectWrapper *          in_results_wrapper_ptr;
    TplDispResults *in_results_ptr;
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->tpl_disp_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (TplDispResults *)in_results_wrapper_ptr->object_ptr;

        PictureParentControlSet* pcs_ptr = in_results_ptr->pcs_ptr;

        SequenceControlSet* scs_ptr = (SequenceControlSet *)pcs_ptr->scs_ptr;

        int32_t frame_idx =in_results_ptr->frame_index;
        context_ptr->coded_sb_count   = 0;

        uint16_t tile_group_width_in_sb = pcs_ptr->tile_group_info[0/*context_ptr->tile_group_index*/] //  1 tile
                                              .tile_group_width_in_sb;
        EncDecSegments *segments_ptr;

        segments_ptr = pcs_ptr->tpl_disp_segment_ctrl[0/*context_ptr->tile_group_index*/]; //  1 tile
    // Segments
    uint16_t        segment_index;

    uint8_t sb_sz      = (uint8_t)scs_ptr->sb_sz ;
    uint8_t sb_size_log2 = (uint8_t)svt_log2f(sb_sz);
    uint32_t pic_width_in_sb = (pcs_ptr->aligned_width + sb_sz - 1) >> sb_size_log2;

    segment_index = 0;
    // no Tiles path
    if (scs_ptr->static_config.tile_rows == 0 && scs_ptr->static_config.tile_columns == 0 ){
        // segments loop
        while (
            assign_tpl_segments(
                segments_ptr,
                &segment_index,
                in_results_ptr,
                frame_idx,
                context_ptr->tpl_disp_fb_fifo_ptr)
            == EB_TRUE) {

            uint32_t        x_sb_start_index;
            uint32_t        y_sb_start_index;
            uint32_t        sb_start_index;
            uint32_t        sb_segment_count;
            uint32_t        sb_segment_index;
            uint32_t        segment_row_index;
            uint32_t        segment_band_index;
            uint32_t        segment_band_size;
            // SB Loop variables
            uint32_t        x_sb_index;
            uint32_t        y_sb_index;

            x_sb_start_index = segments_ptr->x_start_array[segment_index];
            y_sb_start_index = segments_ptr->y_start_array[segment_index];
            sb_start_index = y_sb_start_index * tile_group_width_in_sb + x_sb_start_index;
            sb_segment_count = segments_ptr->valid_sb_count_array[segment_index];

            segment_row_index = segment_index / segments_ptr->segment_band_count;
            segment_band_index =
                segment_index - segment_row_index * segments_ptr->segment_band_count;
            segment_band_size = (segments_ptr->sb_band_count * (segment_band_index + 1) +
                segments_ptr->segment_band_count - 1) /
                segments_ptr->segment_band_count;


            for (y_sb_index = y_sb_start_index, sb_segment_index = sb_start_index;
                sb_segment_index < sb_start_index + sb_segment_count;
                ++y_sb_index) {
                for (x_sb_index = x_sb_start_index;
                    x_sb_index < tile_group_width_in_sb &&
                    (x_sb_index + y_sb_index < segment_band_size) &&
                    sb_segment_index < sb_start_index + sb_segment_count;
                    ++x_sb_index, ++sb_segment_index) {
                    uint16_t tile_group_y_sb_start =
                        pcs_ptr->tile_group_info[0/*context_ptr->tile_group_index*/] //  1 tile
                        .tile_group_sb_start_y;
                    uint16_t tile_group_x_sb_start =
                        pcs_ptr->tile_group_info[0/*context_ptr->tile_group_index*/] //  1 tile
                        .tile_group_sb_start_x;

                    context_ptr->sb_index = (uint16_t)((y_sb_index + tile_group_y_sb_start) * pic_width_in_sb +
                        x_sb_index + tile_group_x_sb_start);

                    // TPL dispenser per SB (64)
                    tpl_mc_flow_dispenser_sb(
                        pcs_ptr->scs_ptr->encode_context_ptr,
                        scs_ptr,
                        pcs_ptr,
                        frame_idx,
                        context_ptr->sb_index,
                        in_results_ptr->qIndex);

                    context_ptr->coded_sb_count++;

                }

                x_sb_start_index = (x_sb_start_index > 0) ? x_sb_start_index - 1 : 0;
            }
        }

        svt_block_on_mutex(pcs_ptr->tpl_disp_mutex);
        pcs_ptr->tpl_disp_coded_sb_count += (uint32_t)context_ptr->coded_sb_count;
        EbBool last_sb_flag = (pcs_ptr->sb_total_count == pcs_ptr->tpl_disp_coded_sb_count);

        svt_release_mutex(pcs_ptr->tpl_disp_mutex);
        if (last_sb_flag)
            svt_post_semaphore(pcs_ptr->tpl_disp_done_semaphore);
    }
    else {
        // Tiles path does not suupport segments
        for (uint32_t sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {

            tpl_mc_flow_dispenser_sb(
                pcs_ptr->scs_ptr->encode_context_ptr,
                scs_ptr,
                pcs_ptr,
                frame_idx,
                sb_index,
                in_results_ptr->qIndex);
        }
        svt_post_semaphore(pcs_ptr->tpl_disp_done_semaphore);

    }
        svt_release_object(in_results_wrapper_ptr);

    }
    return NULL;
}




/************************************************
 * Source Based Operations Kernel
 * Source-based operations process involves a number of analysis algorithms
 * to identify spatiotemporal characteristics of the input pictures.
 ************************************************/
void *source_based_operations_kernel(void *input_ptr) {
    EbThreadContext *             thread_context_ptr = (EbThreadContext *)input_ptr;
    SourceBasedOperationsContext *context_ptr        = (SourceBasedOperationsContext *)
                                                    thread_context_ptr->priv;
    PictureParentControlSet *  pcs_ptr;
    EbObjectWrapper *          in_results_wrapper_ptr;
    InitialRateControlResults *in_results_ptr;
    EbObjectWrapper *          out_results_wrapper_ptr;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->initial_rate_control_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (InitialRateControlResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        context_ptr->complete_sb_count = 0;
        uint32_t sb_total_count        = pcs_ptr->sb_total_count;
        uint32_t sb_index;

        SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        // Get TPL ME

        if (scs_ptr->static_config.enable_tpl_la) {

            if (scs_ptr->static_config.enable_tpl_la &&
                pcs_ptr->temporal_layer_index == 0) {

                tpl_prep_info(pcs_ptr);
                tpl_mc_flow(scs_ptr->encode_context_ptr, scs_ptr, pcs_ptr,context_ptr);
            }
        }

        /***********************************************SB-based operations************************************************************/
        for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
            SbParams *sb_params      = &pcs_ptr->sb_params_array[sb_index];
            EbBool    is_complete_sb = sb_params->is_complete_sb;
            if (is_complete_sb) {
                context_ptr->complete_sb_count++;
            }
        }
        /*********************************************Picture-based operations**********************************************************/

        // Activity statistics derivation
        derive_picture_activity_statistics(pcs_ptr);

        // Get Empty Results Object
        svt_get_empty_object(context_ptr->picture_demux_results_output_fifo_ptr,
                             &out_results_wrapper_ptr);

        PictureDemuxResults *out_results_ptr = (PictureDemuxResults *)
                                                   out_results_wrapper_ptr->object_ptr;
        out_results_ptr->pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
        out_results_ptr->picture_type    = EB_PIC_INPUT;

        // Release the Input Results
        svt_release_object(in_results_wrapper_ptr);

        // Post the Full Results Object
        svt_post_full_object(out_results_wrapper_ptr);
    }
    return NULL;
}
