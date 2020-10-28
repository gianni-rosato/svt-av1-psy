/*
* Copyright(c) 2019 Netflix, Inc.
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "EbDefinitions.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecInverseQuantize.h"
#include "EbObuParse.h"
#include "EbDecProcessFrame.h"
#include "EbCoefficients.h"
#include "EbQMatrices.h"
#include "EbInvTransforms.h"

// Same wrapper(av1_ac/dc_quant_qtx) available in .c file of encoder
static INLINE int16_t get_dc_quant(int32_t qindex, int32_t delta, AomBitDepth bit_depth) {
    return svt_av1_dc_quant_q3(qindex, delta, bit_depth);
}

static INLINE int16_t get_ac_quant(int32_t qindex, int32_t delta, AomBitDepth bit_depth) {
    return svt_av1_ac_quant_q3(qindex, delta, bit_depth);
}

// Called in read_frame_header_obu() -> av1_decode_frame_headers_and_setup() -> read_uncompressed_header()
void setup_segmentation_dequant(DecModCtxt *dec_mod_ctxt) {
    SeqHeader *  seq_header = dec_mod_ctxt->seq_header;
    FrameHeader *frame_info = dec_mod_ctxt->frame_header;
    int          bit_depth  = seq_header->color_config.bit_depth;
    /*int max_segments = frame_info->segmentation_params.segmentation_enabled ?
        MAX_SEGMENTS : 1;*/
    int32_t dc_delta_q, ac_delta_q;
    for (int i = 0; i < MAX_SEGMENTS; i++) {
        const int qindex = get_qindex(
            &frame_info->segmentation_params, i, frame_info->quantization_params.base_q_idx);

        for (int plane = 0; plane < MAX_MB_PLANE; plane++) {
            dc_delta_q = frame_info->quantization_params.delta_q_dc[plane];
            ac_delta_q = frame_info->quantization_params.delta_q_ac[plane];

            dec_mod_ctxt->dequants.dequant_qtx[i][plane][0] =
                get_dc_quant(qindex, dc_delta_q, bit_depth);
            dec_mod_ctxt->dequants.dequant_qtx[i][plane][1] =
                get_ac_quant(qindex, ac_delta_q, bit_depth);
        }
    }
}

void av1_inverse_qm_init(DecModCtxt *dec_mod_ctxt, SeqHeader *seq_header) {
    const int num_planes = av1_num_planes(&seq_header->color_config);
    int       q, c;
    uint8_t   t;
    int       current;
    for (q = 0; q < NUM_QM_LEVELS; ++q) {
        for (c = 0; c < num_planes; ++c) {
            current = 0;
            for (t = 0; t < TX_SIZES_ALL; ++t) {
                const int     size       = tx_size_2d[t];
                const uint8_t qm_tx_size = av1_get_adjusted_tx_size(t);
                if (q == NUM_QM_LEVELS - 1)
                    dec_mod_ctxt->giqmatrix[q][c][t] = NULL;
                else if (t != qm_tx_size) { // Reuse matrices for 'qm_tx_size'
                    dec_mod_ctxt->giqmatrix[q][c][t] = dec_mod_ctxt->giqmatrix[q][c][qm_tx_size];
                } else {
                    assert(current + size <= QM_TOTAL_SIZE);
                    dec_mod_ctxt->giqmatrix[q][c][t] = &iwt_matrix_ref[q][c >= 1][current];
                    current += size;
                }
            }
        }
    }
}

// Called in decode_tiles()
// Default initilaization of dequant and iquant
/*
void av1_init_sb(FrameHeader *frame)
{
    frame->dequants_delta_q = &frame->dequants;
}
*/

// Called in parse_decode_block()
// Update de-quantization parameter based on delta qp param
void update_dequant(DecModCtxt *dec_mod_ctxt, SBInfo *sb_info) {
    SeqHeader *  seq_header = dec_mod_ctxt->seq_header;
    FrameHeader *frame      = dec_mod_ctxt->frame_header;

    int bit_depth                  = seq_header->color_config.bit_depth;
    dec_mod_ctxt->dequants_delta_q = &dec_mod_ctxt->dequants;
    if (frame->delta_q_params.delta_q_present) {
        for (int i = 0; i < MAX_SEGMENTS; i++) {
            const int current_qindex = get_qindex(
                &frame->segmentation_params, i, sb_info->sb_delta_q[0]);

            for (int plane = 0; plane < MAX_MB_PLANE; plane++) {
                const int dc_delta_q = frame->quantization_params.delta_q_dc[plane];
                const int ac_delta_q = frame->quantization_params.delta_q_ac[plane];

                dec_mod_ctxt->dequants_delta_q->dequant_qtx[i][plane][0] =
                    get_dc_quant(current_qindex, dc_delta_q, bit_depth);
                dec_mod_ctxt->dequants_delta_q->dequant_qtx[i][plane][1] =
                    get_ac_quant(current_qindex, ac_delta_q, bit_depth);
            }
        }
    }
}

static INLINE int get_dqv(const int16_t dequant, int coeff_idx, const QmVal *iqmatrix) {
    int dqv = dequant;
    if (iqmatrix != NULL)
        dqv = ((iqmatrix[coeff_idx] * dqv) + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
    return dqv;
}

int32_t inverse_quantize(DecModCtxt *dec_mod_ctxt, PartitionInfo *part, BlockModeInfo *mode,
                         int32_t *level, int32_t *qcoeffs, TxType tx_type, TxSize tx_size,
                         int plane) {
    (void)part;
    SeqHeader *            seq   = dec_mod_ctxt->seq_header;
    FrameHeader *          frame = dec_mod_ctxt->frame_header;
    const ScanOrder *const scan_order =
        &av1_scan_orders[tx_size][tx_type]; //get_scan(tx_size, tx_type);
    const int16_t *scan      = scan_order->scan;
    const int32_t  max_value = (1 << (7 + seq->color_config.bit_depth)) - 1;
    const int32_t  min_value = -(1 << (7 + seq->color_config.bit_depth));
    int            n_coeffs, i, pos, qmlevel;
    int16_t *      dequant;
    const QmVal *  iqmatrix;
    const TxSize   qm_tx_size = av1_get_adjusted_tx_size(tx_size);

    int using_qm = frame->quantization_params.using_qmatrix;
    int lossless = frame->lossless_array[mode->segment_id];
    dequant      = dec_mod_ctxt->dequants_delta_q->dequant_qtx[mode->segment_id][plane];
    qmlevel =
        (lossless || using_qm == 0) ? NUM_QM_LEVELS - 1 : frame->quantization_params.qm[plane];
    iqmatrix = IS_2D_TRANSFORM(tx_type) ? dec_mod_ctxt->giqmatrix[qmlevel][plane][qm_tx_size]
                                        : dec_mod_ctxt->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
    const int shift = av1_get_tx_scale(tx_size);

    // Level is 1D array with eob length as first value then continued by
    // coeffs value to the length of eob.
#if SVT_DEC_COEFF_DEBUG
    int16_t *cur_coeff = (int16_t *)level;
    n_coeffs           = cur_coeff[1];
#else
    n_coeffs = level[0]; // coeffs length
#endif
    level++;

    TranLow q_coeff;
    int32_t lev;
    lev = level[0];
    if (lev) {
        pos     = scan[0];
        q_coeff = (TranLow)((int64_t)abs(lev) * get_dqv(dequant[0], pos, iqmatrix) & 0xffffff);
        q_coeff = q_coeff >> shift;

        if (lev < 0) q_coeff = -q_coeff;
        qcoeffs[0] = clamp(q_coeff, min_value, max_value);
    }

    for (i = 1; i < n_coeffs; i++) {
        lev = level[i];
        if (lev != 0) {
            pos     = scan[i];
            q_coeff = (TranLow)((int64_t)abs(lev) * get_dqv(dequant[1], pos, iqmatrix) & 0xffffff);
            q_coeff = q_coeff >> shift;

            if (lev < 0) q_coeff = -q_coeff;

            qcoeffs[pos] = clamp(q_coeff, min_value, max_value);
        }
    }
    return n_coeffs;
}
