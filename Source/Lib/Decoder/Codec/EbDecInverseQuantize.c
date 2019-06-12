/*
* Copyright(c) 2019 Netflix, Inc.
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

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecInverseQuantize.h"
#include "EbObuParse.h"
#include "EbDecProcessFrame.h"

#include "EbTransforms.h"

// Same wrapper(av1_ac/dc_quant_QTX) available in .c file of encoder
int16_t get_dc_quant(int32_t qindex, int32_t delta, AomBitDepth bit_depth)
{
    return av1_dc_quant_Q3(qindex, delta, bit_depth);
}

int16_t get_ac_quant(int32_t qindex, int32_t delta, AomBitDepth bit_depth)
{
    return av1_ac_quant_Q3(qindex, delta, bit_depth);
}

// Called in read_frame_header_obu() -> av1_decode_frame_headers_and_setup() -> read_uncompressed_header()
void setup_segmentation_dequant(EbDecHandle *dec_handle_ptr, EbColorConfig *color_config)
{
    (void)color_config;
    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    FrameHeader *frame_info = &dec_handle_ptr->frame_header;
    DecModCtxt *dec_mod_ctxt = (DecModCtxt*)dec_handle_ptr->pv_dec_mod_ctxt;
    int bit_depth = seq_header->color_config.bit_depth;
    /*int max_segments = frame_info->segmentation_params.segmentation_enabled ?
        MAX_SEGMENTS : 1;*/
    int32_t qindex;
    for (int i = 0; i < MAX_SEGMENTS; i++) {
        qindex = get_qindex(&frame_info->segmentation_params, i,
            frame_info->quantization_params.base_q_idx);

        // Y plane: DC and AC
        dec_mod_ctxt->dequants.y_dequant_QTX[i][0] = get_dc_quant(qindex,
            frame_info->quantization_params.delta_q_y_dc, bit_depth);
        dec_mod_ctxt->dequants.y_dequant_QTX[i][1] = get_ac_quant(qindex,
            0, bit_depth);

        // U plane: DC and AC
        dec_mod_ctxt->dequants.u_dequant_QTX[i][0] = get_dc_quant(qindex,
            frame_info->quantization_params.delta_q_u_dc, bit_depth);
        dec_mod_ctxt->dequants.u_dequant_QTX[i][1] = get_ac_quant(qindex,
            frame_info->quantization_params.delta_q_u_ac, bit_depth);

        // V plane: DC and AC
        dec_mod_ctxt->dequants.v_dequant_QTX[i][0] = get_dc_quant(qindex,
            frame_info->quantization_params.delta_q_v_dc, bit_depth);
        dec_mod_ctxt->dequants.v_dequant_QTX[i][1] = get_ac_quant(qindex,
            frame_info->quantization_params.delta_q_v_ac, bit_depth);
    }
}

void av1_inverse_qm_init(EbDecHandle *dec_handle_ptr)
{
    DecModCtxt *dec_mod_ctxt = (DecModCtxt*)dec_handle_ptr->pv_dec_mod_ctxt;
    const int num_planes = av1_num_planes(&dec_handle_ptr->seq_header.color_config);
    int q, c;
    uint8_t t;
    int current;
    for (q = 0; q < NUM_QM_LEVELS; ++q) {
        for (c = 0; c < num_planes; ++c) {
            current = 0;
            for (t = 0; t < TX_SIZES_ALL; ++t) {
                const int size = tx_size_2d[t];
                const uint8_t qm_tx_size = av1_get_adjusted_tx_size(t);
                if (q == NUM_QM_LEVELS - 1)
                    dec_mod_ctxt->giqmatrix[q][c][t] = NULL;
                else if (t != qm_tx_size) {  // Reuse matrices for 'qm_tx_size'
                    dec_mod_ctxt->giqmatrix[q][c][t] =
                        dec_mod_ctxt->giqmatrix[q][c][qm_tx_size];
                }
                else {
                    assert(current + size <= QM_TOTAL_SIZE);
                    dec_mod_ctxt->giqmatrix[q][c][t] =
                        &iwt_matrix_ref[q][c >= 1][current];
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
void update_dequant(EbDecHandle *dec_handle, SBInfo *sb_info)
{
    int32_t current_qindex;
    int dc_delta_q, ac_delta_q;
    SeqHeader *seq_header = &dec_handle->seq_header;
    FrameHeader *frame = &dec_handle->frame_header;
    DecModCtxt *dec_mod_ctxt = (DecModCtxt*)dec_handle->pv_dec_mod_ctxt;

    if (!frame->delta_q_params.delta_q_present)
        dec_mod_ctxt->dequants_delta_q = &dec_mod_ctxt->dequants;
    else {
        for (int i = 0; i < MAX_SEGMENTS; i++) {
            current_qindex =
                get_qindex(&frame->segmentation_params, i,
                    clamp(frame->quantization_params.base_q_idx +
                    sb_info->sb_delta_q[0], 1, MAXQ));

            // Y Plane: AC and DC
            dc_delta_q = frame->quantization_params.delta_q_y_dc;
            ac_delta_q = 0;
            dec_mod_ctxt->dequants_delta_q->y_dequant_QTX[i][0] = get_dc_quant(
                current_qindex, dc_delta_q, seq_header->color_config.bit_depth);
            dec_mod_ctxt->dequants_delta_q->y_dequant_QTX[i][1] = get_ac_quant(
                current_qindex, ac_delta_q, seq_header->color_config.bit_depth);

            // U Plane: AC and DC
            dc_delta_q = frame->quantization_params.delta_q_u_dc;
            ac_delta_q = frame->quantization_params.delta_q_u_ac;
            dec_mod_ctxt->dequants_delta_q->u_dequant_QTX[i][0] = get_dc_quant(
                current_qindex, dc_delta_q, seq_header->color_config.bit_depth);
            dec_mod_ctxt->dequants_delta_q->u_dequant_QTX[i][1] = get_ac_quant(
                current_qindex, ac_delta_q, seq_header->color_config.bit_depth);

            // V Plane: AC and DC
            dc_delta_q = frame->quantization_params.delta_q_v_dc;
            ac_delta_q = frame->quantization_params.delta_q_v_ac;
            dec_mod_ctxt->dequants_delta_q->v_dequant_QTX[i][0] = get_dc_quant(
                current_qindex, dc_delta_q, seq_header->color_config.bit_depth);
            dec_mod_ctxt->dequants_delta_q->v_dequant_QTX[i][1] = get_ac_quant(
                current_qindex, ac_delta_q, seq_header->color_config.bit_depth);
        }
    }
}

int get_dqv(const int16_t *dequant, int coeff_idx, const QmVal *iqmatrix) {
    int dqv = dequant[!!coeff_idx];
    if (iqmatrix != NULL)
        dqv =
        ((iqmatrix[coeff_idx] * dqv) + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
    return dqv;
}

int32_t inverse_quantize(EbDecHandle * dec_handle, PartitionInfo_t *part, ModeInfo_t *mode,
    int32_t *level, int32_t *qcoeffs, TxType tx_type, TxSize tx_size, int plane)
{
    (void)part;
    SeqHeader *seq = &dec_handle->seq_header;
    FrameHeader *frame = &dec_handle->frame_header;
    DecModCtxt *dec_mod_ctxt = (DecModCtxt*)dec_handle->pv_dec_mod_ctxt;
    const ScanOrder *const scan_order = &av1_scan_orders[tx_size][tx_type]; //get_scan(tx_size, tx_type);
    const int16_t *scan = scan_order->scan;
    const int32_t max_value = (1 << (7 + seq->color_config.bit_depth)) - 1;
    const int32_t min_value = -(1 << (7 + seq->color_config.bit_depth));
    int n_coeffs, i, pos, qmlevel;
    int16_t *dequant;
    const QmVal *iqmatrix;
    const TxSize qm_tx_size = av1_get_adjusted_tx_size(tx_size);

    int using_qm = frame->quantization_params.using_qmatrix;
    int lossless = frame->lossless_array[mode->segment_id];
    if (plane == 0) {
        qmlevel = (lossless || using_qm == 0) ? NUM_QM_LEVELS - 1 :
            frame->quantization_params.qm_y;
        iqmatrix = IS_2D_TRANSFORM(tx_type)
            ? dec_mod_ctxt->giqmatrix[qmlevel][AOM_PLANE_Y][qm_tx_size]
            : dec_mod_ctxt->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
        dequant = dec_mod_ctxt->dequants_delta_q->y_dequant_QTX[mode->segment_id];
    }
    else if (plane == 1) {
        qmlevel = (lossless || using_qm == 0) ? NUM_QM_LEVELS - 1 :
            frame->quantization_params.qm_u;
        iqmatrix = IS_2D_TRANSFORM(tx_type)
            ? dec_mod_ctxt->giqmatrix[qmlevel][AOM_PLANE_U][qm_tx_size]
            : dec_mod_ctxt->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
        dequant = dec_mod_ctxt->dequants_delta_q->u_dequant_QTX[mode->segment_id];
    }
    else {
        qmlevel = (lossless || using_qm == 0) ? NUM_QM_LEVELS - 1 :
            frame->quantization_params.qm_v;
        iqmatrix = IS_2D_TRANSFORM(tx_type)
            ? dec_mod_ctxt->giqmatrix[qmlevel][AOM_PLANE_V][qm_tx_size]
            : dec_mod_ctxt->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
        dequant = dec_mod_ctxt->dequants_delta_q->v_dequant_QTX[mode->segment_id];
    }

    const int shift = av1_get_tx_scale(tx_size);

    // Level is 1D array with eob length as first value then continued by
    // coeffs value to the length of eob.
#if SVT_DEC_COEFF_DEBUG
    int16_t *cur_coeff = (int16_t *)level;
    n_coeffs = cur_coeff[1];
#else
    n_coeffs = level[0]; // coeffs length
#endif
    level++;

    memset(qcoeffs, 0, (1 << seq->sb_size_log2) * (1 << seq->sb_size_log2) * sizeof(*qcoeffs));

    for (i = 0; i < n_coeffs; i++) {
        pos = scan[i];
        qcoeffs[pos] = (TranLow)((int64_t)abs(level[i]) *
            get_dqv(dequant, pos, iqmatrix) & 0xffffff);
        qcoeffs[pos] = qcoeffs[pos] >> shift;

        if (level[i] < 0)
            qcoeffs[pos] = -qcoeffs[pos];

        qcoeffs[pos] = clamp(qcoeffs[pos], min_value, max_value);
    }
    return n_coeffs;
}
