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

#include "EbDefinitions.h"
#include "EbModeDecisionProcess.h"
#include "EbTransforms.h"
#include "EbFullLoop.h"
#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"
#if QT_10BIT_SUPPORT
#ifdef __GNUC__
#define LIKELY(v) __builtin_expect(v, 1)
#define UNLIKELY(v) __builtin_expect(v, 0)
#else
#define LIKELY(v) (v)
#define UNLIKELY(v) (v)
#endif
#endif
static PartitionType from_shape_to_part[] =
{
PARTITION_NONE,
PARTITION_HORZ,
PARTITION_VERT,
PARTITION_HORZ_A,
PARTITION_HORZ_B,
PARTITION_VERT_A,
PARTITION_VERT_B,
PARTITION_HORZ_4,
PARTITION_VERT_4,
PARTITION_SPLIT

};
void quantize_b_helper_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
    uint16_t *eob_ptr, const int16_t *scan,
    const int16_t *iscan, const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr, const int32_t log_scale) {
    const int32_t zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
                           ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale) };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t i, non_zero_count = (int32_t)n_coeffs, eob = -1;
    (void)iscan;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        // Pre-scan pass
        for (i = (int32_t)n_coeffs - 1; i >= 0; i--) {
            const int32_t rc = scan[i];
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;

            if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
                coeff >(nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                non_zero_count--;
            else
                break;
        }

        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < non_zero_count; i++) {
            const int32_t rc = scan[i];
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            int32_t tmp32;

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
                int64_t tmp =
                    clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale),
                        INT16_MIN, INT16_MAX);
                tmp *= wt;
                tmp32 = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
                    quant_shift_ptr[rc != 0]) >>
                    (16 - log_scale + AOM_QM_BITS));  // quantization
                qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
                const int32_t iwt = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
                const int32_t dequant =
                    (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                    AOM_QM_BITS;
                const tran_low_t abs_dqcoeff = (tmp32 * dequant) >> log_scale;
                dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);

                if (tmp32) eob = i;
            }
        }
    }
    *eob_ptr = (uint16_t)(eob + 1);
}
void aom_quantize_b_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
    uint16_t *eob_ptr, const int16_t *scan,
    const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 0);
}

void aom_quantize_b_32x32_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 1);
}

void aom_quantize_b_64x64_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 2);
}

void quantize_b_helper_c(
    const tran_low_t *coeff_ptr,
    int32_t stride,
#
    int32_t width,
    int32_t height,

    intptr_t n_coeffs,
    int32_t skip_block,
    const int16_t *zbin_ptr,
    const int16_t *round_ptr,
    const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr,
    uint16_t *eob_ptr,
    const int16_t *scan,
    const int16_t *iscan,
    const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr,
    const int32_t log_scale)
{

    const int32_t zbins[2] = {
        ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
        ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale)
    };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t i, non_zero_count = (int32_t)n_coeffs, eob = -1;
    (void)iscan;

    // Nader quantisation
    for (int32_t x = 0; x < height; x++) {
        memset(qcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*qcoeff_ptr));
        memset(dqcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*dqcoeff_ptr));
    }


    if (!skip_block) {
        // Pre-scan pass
        for (i = (int32_t)n_coeffs - 1; i >= 0; i--) {
            const int32_t mapRc = scan[i];

            const int32_t rc = ((mapRc / MIN(32, height))  * stride) + (mapRc % MIN(32, width));

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;


            ////if (mapRc != NewTab[rc])
            //printf("%d\n", coeff);

            if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
                coeff >(nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                non_zero_count--;
            else
                break;
        }
        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < non_zero_count; i++) {
            const int32_t mapRc = scan[i];

            const int32_t rc = ((mapRc / MIN(32, height))  * stride) + (mapRc % MIN(32, width));
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            int32_t tmp32;

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[mapRc] : (1 << AOM_QM_BITS);

            if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {

                int64_t tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale), INT16_MIN, INT16_MAX);

                tmp *= wt;

                tmp32 = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *    quant_shift_ptr[rc != 0]) >> (16 - log_scale + AOM_QM_BITS));  // quantization

                qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;

                const int32_t iwt = iqm_ptr != NULL ? iqm_ptr[mapRc] : (1 << AOM_QM_BITS);

                const int32_t dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

                dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant / (1 << log_scale);

                if (tmp32) eob = i;
            }
        }
    }


    *eob_ptr = (uint16_t)(eob + 1);
}
#if QT_10BIT_SUPPORT
void highbd_quantize_b_helper_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan, const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr, const int32_t log_scale) {
    int32_t i, eob = -1;
    const int32_t zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
        ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale) };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t dequant;
    int32_t idx_arr[4096];
    (void)iscan;
    int32_t idx = 0;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        // Pre-scan pass
        for (i = 0; i < n_coeffs; i++) {
            const int32_t rc = scan[i];
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;

            // If the coefficient is out of the base ZBIN range, keep it for
            // quantization.
            if (coeff >= (zbins[rc != 0] * (1 << AOM_QM_BITS)) ||
                coeff <= (nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                idx_arr[idx++] = i;
        }

        // Quantization pass: only process the coefficients selected in
        // pre-scan pass. Note: idx can be zero.
        for (i = 0; i < idx; i++) {
            const int32_t rc = scan[idx_arr[i]];
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const qm_val_t iwt = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            const int64_t tmp1 =
                abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale);
            const int64_t tmpw = tmp1 * wt;
            const int64_t tmp2 = ((tmpw * quant_ptr[rc != 0]) >> 16) + tmpw;
            const int32_t abs_qcoeff = (int32_t)((tmp2 * quant_shift_ptr[rc != 0]) >>
                (16 - log_scale + AOM_QM_BITS));
            qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
            dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;
            const tran_low_t abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
            dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
            if (abs_qcoeff) eob = idx_arr[i];
        }
    }
    *eob_ptr = (uint16_t)(eob + 1);
}

void aom_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 0);
}

void aom_highbd_quantize_b_32x32_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 1);
}

void aom_highbd_quantize_b_64x64_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 2);
}

void av1_highbd_quantize_b_facade(const tran_low_t *coeff_ptr,
    intptr_t n_coeffs, const MacroblockPlane *p,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
    const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam) {
    // obsolete skip_block
    const int32_t skip_block = 0;
    const qm_val_t *qm_ptr = qparam->qmatrix;
    const qm_val_t *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
            p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
            qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
            sc->scan, sc->iscan, qm_ptr, iqm_ptr,
            qparam->log_scale);
    }
    else {
        switch (qparam->log_scale) {
        case 0:
            if (LIKELY(n_coeffs >= 8)) {

                aom_highbd_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                    p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                    qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX,
                    eob_ptr, sc->scan, sc->iscan);

            }
            else {

                aom_highbd_quantize_b_c(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                    p->round_QTX, p->quant_QTX,
                    p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                    p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            }
            break;
        case 1:
            aom_highbd_quantize_b_32x32(
                coeff_ptr, n_coeffs, skip_block, p->zbin_QTX, p->round_QTX,
                p->quant_QTX, p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            break;
        case 2:
            aom_highbd_quantize_b_64x64(
                coeff_ptr, n_coeffs, skip_block, p->zbin_QTX, p->round_QTX,
                p->quant_QTX, p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            break;
        default: assert(0);
        }
    }
}

/*
static INLINE void highbd_quantize_dc(
    const tran_low_t *coeff_ptr, int32_t n_coeffs, int32_t skip_block,
    const int16_t *round_ptr, const int16_t quant, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t dequant_ptr, uint16_t *eob_ptr,
    const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr, const int32_t log_scale) {
    int32_t eob = -1;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        const qm_val_t wt = qm_ptr != NULL ? qm_ptr[0] : (1 << AOM_QM_BITS);
        const qm_val_t iwt = iqm_ptr != NULL ? iqm_ptr[0] : (1 << AOM_QM_BITS);
        const int32_t coeff = coeff_ptr[0];
        const int32_t coeff_sign = (coeff >> 31);
        const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
        const int64_t tmp = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[0], log_scale);
        const int64_t tmpw = tmp * wt;
        const int32_t abs_qcoeff =
            (int32_t)((tmpw * quant) >> (16 - log_scale + AOM_QM_BITS));
        qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        const int32_t dequant =
            (dequant_ptr * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

        const tran_low_t abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
        dqcoeff_ptr[0] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
        if (abs_qcoeff) eob = 0;
    }
    *eob_ptr = (uint16_t)(eob + 1);
}
*/

#endif
#if !QT_10BIT_SUPPORT
/* These functions should only be called when quantisation matrices
are not used. */
void aom_quantize_b_c(const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t txb_size,
    intptr_t n_coeffs,
    int32_t skip_block,
    const int16_t *zbin_ptr,
    const int16_t *round_ptr,
    const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr,
    uint16_t *eob_ptr,
    const int16_t *scan,
    const int16_t *iscan)
{
    quantize_b_helper_c(
        coeff_ptr,
        stride, txb_size,
        n_coeffs,
        skip_block,
        zbin_ptr,
        round_ptr,
        quant_ptr,
        quant_shift_ptr,
        qcoeff_ptr,
        dqcoeff_ptr,
        dequant_ptr,
        eob_ptr, scan,
        iscan,
        NULL,
        NULL,
        0);
}

void aom_quantize_b_32x32_c(const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t txb_size,
    intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c(coeff_ptr,
        stride,
        txb_size,
        n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 1);
}

void aom_quantize_b_64x64_c(const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t txb_size,
    intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c(coeff_ptr,
        stride,
        txb_size,
        n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 2);
}
void av1_quantize_b_facade(
    const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t txb_size,
    intptr_t n_coeffs,

    const MacroblockPlane *p,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    uint16_t *eob_ptr,
    const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam)
{
    // obsolete skip_block
    const int32_t skip_block = 0;
    const qm_val_t *qm_ptr = qparam->qmatrix;
    const qm_val_t *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        quantize_b_helper_c(
            coeff_ptr,
            stride,
            txb_size,
            n_coeffs,
            skip_block,
            p->zbin_QTX,
            p->round_QTX,
            p->quant_QTX,
            p->quant_shift_QTX,
            qcoeff_ptr,
            dqcoeff_ptr,
            p->dequant_QTX,
            eob_ptr,
            sc->scan,
            sc->iscan,
            qm_ptr,
            iqm_ptr,
            qparam->log_scale);
    }
    else {
        switch (qparam->log_scale) {
        case 0:
            aom_quantize_b_c(
                coeff_ptr,
                stride,
                txb_size,
                n_coeffs,
                skip_block,
                p->zbin_QTX,
                p->round_QTX,
                p->quant_QTX,
                p->quant_shift_QTX,
                qcoeff_ptr,
                dqcoeff_ptr,
                p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);
            break;
        case 1:
            aom_quantize_b_32x32_c(
                coeff_ptr,
                stride,
                txb_size,
                n_coeffs,
                skip_block,
                p->zbin_QTX,
                p->round_QTX,
                p->quant_QTX,
                p->quant_shift_QTX,
                qcoeff_ptr,
                dqcoeff_ptr,
                p->dequant_QTX,
                eob_ptr,
                sc->scan,
                sc->iscan);
            break;
        case 2:
            aom_quantize_b_64x64_c(
                coeff_ptr,
                stride,
                txb_size,
                n_coeffs,
                skip_block,
                p->zbin_QTX,
                p->round_QTX,
                p->quant_QTX,
                p->quant_shift_QTX,
                qcoeff_ptr,
                dqcoeff_ptr,
                p->dequant_QTX,
                eob_ptr,
                sc->scan,
                sc->iscan);
            break;
        default: assert(0);
        }
    }
}

#endif
void av1_quantize_b_facade_II(
    const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t                width,
    int32_t                height,
    intptr_t n_coeffs,

    const MacroblockPlane *p,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    uint16_t *eob_ptr,
    const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam)
{
    // obsolete skip_block
    const int32_t skip_block = 0;
    const qm_val_t *qm_ptr = qparam->qmatrix;
    const qm_val_t *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        quantize_b_helper_c(
            coeff_ptr,
            stride,
            width,
            height,
            n_coeffs,
            skip_block,
            p->zbin_QTX,
            p->round_QTX,
            p->quant_QTX,
            p->quant_shift_QTX,
            qcoeff_ptr,
            dqcoeff_ptr,
            p->dequant_QTX,
            eob_ptr,
            sc->scan,
            sc->iscan,
            qm_ptr,
            iqm_ptr,
            qparam->log_scale);
    }
    else {
        switch (qparam->log_scale) {
        case 0:
            aom_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        case 1:

            aom_quantize_b_32x32(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        case 2:

            aom_quantize_b_64x64(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        default: assert(0);
        }
    }
}

/*********************************************************************
* UnifiedQuantizeInvQuantize
*
*  Unified Quant +iQuant
*********************************************************************/
void Av1QuantizeInvQuantize_II(
    PictureControlSet_t  *picture_control_set_ptr,
    int32_t               *coeff,
    const uint32_t          coeffStride,
    int32_t               *quantCoeff,
    int32_t               *reconCoeff,
    uint32_t                qp,
    uint32_t                width,
    uint32_t                height,
    TxSize              transform_size,
    uint16_t                *eob,
    //  MacroblockPlane      candidate_plane,
    EbAsm                asm_type,
    uint32_t                *y_count_non_zero_coeffs,
    EbPfMode              pfMode,
    uint8_t                 enableContouringQCUpdateFlag,
    uint32_t                componentType,
#if QT_10BIT_SUPPORT
    uint32_t                bitIncrement,
#endif
    TxType               tx_type,
    EbBool               cleanSparseCoeffFlag)
{
    (void)cleanSparseCoeffFlag;
    (void)enableContouringQCUpdateFlag;
    (void)pfMode;
    (void)asm_type;
#if !ADD_DELTA_QP_SUPPORT
    (void) qp;
#endif
#if MACRO_BLOCK_CLEANUP
    MacroblockPlane      candidate_plane ;
#else
    MacroblockPlane      candidate_plane = { 0 };
#endif

    //    EB_SLICE          slice_type = picture_control_set_ptr->slice_type;
    //    uint32_t            temporal_layer_index = picture_control_set_ptr->temporal_layer_index;

    // Use a flat matrix (i.e. no weighting) for 1D and Identity transforms
    //const qm_val_t *qmatrix =
    //    IS_2D_TRANSFORM(tx_type) ? pd->seg_qmatrix[seg_id][qm_tx_size]
    //    : cm->gqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
    //const qm_val_t *iqmatrix =
    //    IS_2D_TRANSFORM(tx_type)
    //    ? pd->seg_iqmatrix[seg_id][qm_tx_size]
    //    : cm->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];

    const qm_val_t *qMatrix = picture_control_set_ptr->parent_pcs_ptr->gqmatrix[NUM_QM_LEVELS - 1][0][transform_size];
    const qm_val_t *iqMatrix = picture_control_set_ptr->parent_pcs_ptr->giqmatrix[NUM_QM_LEVELS - 1][0][transform_size];
#if ADD_DELTA_QP_SUPPORT
    uint32_t qIndex = qp;
#else
    uint32_t qIndex = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
#endif
#if MD_10BIT_FIX
    if (bitIncrement == 0) {
        if (componentType == COMPONENT_LUMA) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.y_dequant_QTX[qIndex];
        }

        if (componentType == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.u_dequant_QTX[qIndex];

        }

        if (componentType == COMPONENT_CHROMA_CR) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.v_dequant_QTX[qIndex];

        }

    }
    else {
        if (componentType == COMPONENT_LUMA) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[qIndex];
        }

        if (componentType == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.u_dequant_QTX[qIndex];

        }

        if (componentType == COMPONENT_CHROMA_CR) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.v_dequant_QTX[qIndex];
        }
    }
#else
    if (componentType == COMPONENT_LUMA) {
        candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant[qIndex];
        candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_fp[qIndex];
        candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round_fp[qIndex];
        candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_shift[qIndex];
        candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_zbin[qIndex];
        candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round[qIndex];
        candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[qIndex];
    }

    if (componentType == COMPONENT_CHROMA_CB) {
        candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant[qIndex];
        candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_fp[qIndex];
        candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round_fp[qIndex];
        candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_shift[qIndex];
        candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_zbin[qIndex];
        candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round[qIndex];
        candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.u_dequant_QTX[qIndex];

    }

    if (componentType == COMPONENT_CHROMA_CR) {
        candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant[qIndex];
        candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_fp[qIndex];
        candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round_fp[qIndex];
        candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_shift[qIndex];
        candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_zbin[qIndex];
        candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round[qIndex];
        candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.v_dequant_QTX[qIndex];

        }
#endif
    const SCAN_ORDER *const scan_order = &av1_scan_orders[transform_size][tx_type];  //get_scan(tx_size, tx_type);

    const int32_t n_coeffs = av1_get_max_eob(transform_size);

    QUANT_PARAM qparam;

    qparam.log_scale = av1_get_tx_scale(transform_size);
    qparam.tx_size = transform_size;
    qparam.qmatrix = qMatrix;
    qparam.iqmatrix = iqMatrix;

#if QT_10BIT_SUPPORT
    if (bitIncrement)
        av1_highbd_quantize_b_facade(
        (tran_low_t*)coeff,
            n_coeffs,
            &candidate_plane,
            quantCoeff,
            (tran_low_t*)reconCoeff,
            eob,
            scan_order,
            &qparam);
    else
        av1_quantize_b_facade_II(
        (tran_low_t*)coeff,
            coeffStride,
            width,
            height,
            n_coeffs,
            &candidate_plane,
            quantCoeff,
            (tran_low_t*)reconCoeff,
            eob,
            scan_order,
            &qparam);
#else
    av1_quantize_b_facade_II(
        (tran_low_t*)coeff,
        coeffStride,
        areaSize,
        n_coeffs,
        &candidate_plane,
        quantCoeff,
        (tran_low_t*)reconCoeff,
        eob,
        scan_order,
        &qparam);
#endif
    *y_count_non_zero_coeffs = *eob;

    }
void Av1QuantizeInvQuantize(
    PictureControlSet_t  *picture_control_set_ptr,
    int32_t               *coeff,
    const uint32_t          coeffStride,
    int32_t               *quantCoeff,
    int32_t               *reconCoeff,
    uint32_t                qp,
    uint32_t                width,
    uint32_t                height,
    TxSize               txsize,
    uint16_t                *eob,
    MacroblockPlane      candidate_plane,
    EbAsm                asm_type,
    uint32_t                *y_count_non_zero_coeffs,
    EbPfMode              pfMode,
    uint8_t                 enableContouringQCUpdateFlag,
    uint32_t                componentType,
#if QT_10BIT_SUPPORT
    uint32_t                bitIncrement,
#endif
    TxType               tx_type,
    EbBool               cleanSparseCoeffFlag)
{
    (void)coeffStride;
    (void)candidate_plane;
    (void)enableContouringQCUpdateFlag;
    (void)pfMode;
    //Note: Transformed, Quantized, iQuantized coeff are stored in 1D fashion. 64x64 is hence in the first 32x32 corner.

    uint32_t i;

    for (i = 0; i < height; i++)
    {
        memset(quantCoeff + i * width, 0, width * sizeof(int32_t));
        memset(reconCoeff + i * width, 0, width * sizeof(int32_t));
    }




    Av1QuantizeInvQuantize_II(
        picture_control_set_ptr,
        coeff,
        0,
        quantCoeff,
        reconCoeff,
        qp,
        width,
        height,
        txsize,
        &eob[0],
        asm_type,
        y_count_non_zero_coeffs,
        0,
        0,
        componentType,
#if QT_10BIT_SUPPORT
        bitIncrement,
#endif
        tx_type,
        cleanSparseCoeffFlag);



}

/****************************************
 ************  Full loop ****************
****************************************/
void ProductFullLoop(
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint32_t                          qp,
    uint32_t                          *y_count_non_zero_coeffs,
    uint64_t                         *y_coeff_bits,
    uint64_t                         *y_full_distortion)
{
    uint32_t                       tuOriginIndex;
    uint64_t                      yFullCost;
    SequenceControlSet_t        *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    EbAsm                         asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    EbBool                      cleanSparseCoeffFlag = EB_FALSE;

    //    uint32_t   currentTuIndex,tuIt;
    uint64_t   yTuCoeffBits;
    uint64_t   tuFullDistortion[3][DIST_CALC_TOTAL];
    context_ptr->three_quad_energy = 0;
    uint32_t  txb_1d_offset = 0;
    uint32_t txb_itr = 0;
    for (txb_itr = 0; txb_itr < context_ptr->blk_geom->txb_count; txb_itr++)
    {
        uint16_t tx_org_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        uint16_t tx_org_y = context_ptr->blk_geom->tx_org_y[txb_itr];
        tuOriginIndex = tx_org_x + (tx_org_y * candidateBuffer->residual_ptr->strideY);
        yTuCoeffBits = 0;

        // Y: T Q iQ
        Av1EstimateTransform(
            &(((int16_t*)candidateBuffer->residual_ptr->bufferY)[tuOriginIndex]),
            candidateBuffer->residual_ptr->strideY,
            &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferY)[txb_1d_offset]),
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[txb_itr],
            &context_ptr->three_quad_energy,
            context_ptr->transform_inner_array_ptr,
            0,
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
            asm_type,
            PLANE_TYPE_Y,
            context_ptr->pf_md_mode);

#if SIMULATE_PF_N2     
        int32_t* transCoeffBuffer = &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferY)[txb_1d_offset]);
        int32_t tu_size;
        switch (context_ptr->blk_geom->txsize[txb_itr]) {
        case TX_64X64:
            tu_size = 64;
            break;
        case TX_32X32:
            tu_size = 32;
            break;
        case TX_16X16:
            tu_size = 16;
            break;
        case TX_8X8:
            tu_size = 8;
            break;
        case TX_4X4:
            tu_size = 4;
            break;
        default: assert(0); break;}

        if ((tu_size == 32 || tu_size == 16 || tu_size == 8 || tu_size == 4) && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0) {
            for (int i = 0; i < (tu_size* tu_size); i++) {
                if (i % tu_size >= (tu_size >> 1) || i / tu_size >= (tu_size >> 1)) {
                    transCoeffBuffer[i] = 0;
                }
            }
        }
#endif
        Av1QuantizeInvQuantize(
            picture_control_set_ptr,
            &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferY)[txb_1d_offset]),
            NOT_USED_VALUE,
            &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferY)[txb_1d_offset]),
            &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferY)[txb_1d_offset]),
            qp,
            context_ptr->blk_geom->tx_width[txb_itr],
            context_ptr->blk_geom->tx_height[txb_itr],
            context_ptr->blk_geom->txsize[txb_itr],
            &candidateBuffer->candidate_ptr->eob[0][txb_itr],
            candidateBuffer->candidate_ptr->candidate_plane[0],
            asm_type,
            &(y_count_non_zero_coeffs[txb_itr]),
            context_ptr->pf_md_mode,
            0,
            COMPONENT_LUMA,
#if QT_10BIT_SUPPORT
            BIT_INCREMENT_8BIT,
#endif
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
            cleanSparseCoeffFlag);

        candidateBuffer->candidate_ptr->quantized_dc[0] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferY)[txb_1d_offset]);



        // LUMA DISTORTION
        PictureFullDistortion32Bits(
            context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr,
            txb_1d_offset,
            0,
            candidateBuffer->reconCoeffPtr,
            txb_1d_offset,
            0,
            context_ptr->blk_geom->tx_width[txb_itr],
            context_ptr->blk_geom->tx_height[txb_itr],
            NOT_USED_VALUE,
            NOT_USED_VALUE,
            tuFullDistortion[0],
            NOT_USED_VALUE,
            NOT_USED_VALUE,
            y_count_non_zero_coeffs[txb_itr],
            0,
            0,
            COMPONENT_LUMA,
            asm_type);


        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;
        //assert(context_ptr->three_quad_energy == 0 && context_ptr->cu_stats->size < 64);
        TxSize    txSize = context_ptr->blk_geom->txsize[0];
        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);


        //LUMA-ONLY
        Av1TuEstimateCoeffBits(
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            txb_1d_offset,
            0,
            context_ptr->coeff_est_entropy_coder_ptr,
            candidateBuffer->residualQuantCoeffPtr,
            y_count_non_zero_coeffs[txb_itr],
            0,
            0,
            &yTuCoeffBits,
            &yTuCoeffBits,
            &yTuCoeffBits,
            context_ptr->blk_geom->txsize[0],
            context_ptr->blk_geom->txsize_uv[0],
            COMPONENT_LUMA,
            asm_type);



        //TODO: fix cbf decision
        Av1TuCalcCostLuma(
            context_ptr->cu_ptr->luma_txb_skip_context,//this should be updated here.
            candidateBuffer->candidate_ptr,
            txb_itr,
            context_ptr->blk_geom->txsize[0],
            y_count_non_zero_coeffs[txb_itr],
            tuFullDistortion[0],      //gets updated inside based on cbf decision
            &yTuCoeffBits,            //gets updated inside based on cbf decision
            &yFullCost,
            context_ptr->full_lambda);


        (*y_coeff_bits) += yTuCoeffBits;

        y_full_distortion[DIST_CALC_RESIDUAL] += tuFullDistortion[0][DIST_CALC_RESIDUAL];
        y_full_distortion[DIST_CALC_PREDICTION] += tuFullDistortion[0][DIST_CALC_PREDICTION];
        txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];



    }
}
#if FAST_TX_SEARCH
// T1
uint8_t allowed_tx_set_a[TX_SIZES_ALL][TX_TYPES] = {
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    1,    0,    1},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    1,    0,    1},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    1,    0,    1,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0} };

uint8_t allowed_tx_set_b[TX_SIZES_ALL][TX_TYPES] = {
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{0,    0,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    0,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0}
};
#endif
void ProductFullLoopTxSearch(
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    PictureControlSet_t            *picture_control_set_ptr)
{
    uint32_t                       tuOriginIndex;
    SequenceControlSet_t          *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    EbAsm                          asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    EbBool                         cleanSparseCoeffFlag = EB_FALSE;
    uint64_t                       yTuCoeffBits;
    uint64_t                       tuFullDistortion[3][DIST_CALC_TOTAL];
    int32_t                        plane = 0;
    const int32_t                  is_inter = candidateBuffer->candidate_ptr->type == INTER_MODE ? EB_TRUE : EB_FALSE;//is_inter_block(mbmi);
    uint64_t                       bestFullCost = UINT64_MAX;
    uint64_t                       yFullCost = MAX_CU_COST;
    uint32_t                       yCountNonZeroCoeffsTemp;
    TxType                         txk_start = DCT_DCT;
    TxType                         txk_end = TX_TYPES;
    TxType                         tx_type;
    int32_t                        txb_itr = 0;
    TxSize                         txSize = context_ptr->blk_geom->txsize[txb_itr];
    const TxSetType                tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);


    int32_t allowed_tx_mask[TX_TYPES] = { 0 };  // 1: allow; 0: skip.
    int32_t allowed_tx_num = 0;
    TxType uv_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        allowed_tx_mask[tx_type] = 1;
        if (plane == 0) {
            if (allowed_tx_mask[tx_type]) {
                const TxType ref_tx_type = ((!av1_ext_tx_used[tx_set_type][tx_type]) || txsize_sqr_up_map[txSize] > TX_32X32) ? DCT_DCT : tx_type;
                if (tx_type != ref_tx_type) {

                    allowed_tx_mask[tx_type] = 0;
                }
            }
        }

        allowed_tx_num += allowed_tx_mask[tx_type];
    }
    // Need to have at least one transform type allowed.
    if (allowed_tx_num == 0) {
        allowed_tx_mask[plane ? uv_tx_type : DCT_DCT] = 1;
    }
#if BUG_FIX
    TxType best_tx_type = DCT_DCT;
#endif
    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        if (!allowed_tx_mask[tx_type]) continue;
#if TX_SEARCH_LEVELS
        if (picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;
#else
#if FAST_TX_SEARCH
#if ENCODER_MODE_CLEANUP
        if (picture_control_set_ptr->enc_mode == ENC_M1)
#endif
         if (!allowed_tx_set_a[txSize][tx_type]) continue;
#endif
#endif
        context_ptr->three_quad_energy = 0;
        uint32_t txb_itr = 0;
        for (txb_itr = 0; txb_itr < context_ptr->blk_geom->txb_count; txb_itr++)


        {
            tuOriginIndex = context_ptr->blk_geom->origin_x + (context_ptr->blk_geom->origin_y * candidateBuffer->residual_ptr->strideY);
            yTuCoeffBits = 0;
            // Y: T Q iQ
            Av1EstimateTransform(
                &(((int16_t*)candidateBuffer->residual_ptr->bufferY)[tuOriginIndex]),
                candidateBuffer->residual_ptr->strideY,

                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferY)[tuOriginIndex]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                tx_type,
                asm_type,
                PLANE_TYPE_Y,
                context_ptr->pf_md_mode);

            Av1QuantizeInvQuantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferY)[tuOriginIndex]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferY)[tuOriginIndex]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferY)[tuOriginIndex]),
                context_ptr->cu_ptr->qp,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight,
                context_ptr->blk_geom->txsize[txb_itr],
                &candidateBuffer->candidate_ptr->eob[0][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[0],
                asm_type,
                &yCountNonZeroCoeffsTemp,
                context_ptr->pf_md_mode,
                0,
                COMPONENT_LUMA,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_8BIT,
#endif
                tx_type,
                cleanSparseCoeffFlag);

            candidateBuffer->candidate_ptr->quantized_dc[0] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferY)[tuOriginIndex]);


#if TX_TYPE_FIX
            //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
            if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {

                continue;
            }
#endif


            // LUMA DISTORTION
            PictureFullDistortion32Bits(
                context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr,
                tuOriginIndex,
                0,
                candidateBuffer->reconCoeffPtr,
                tuOriginIndex,
                0,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv,
                tuFullDistortion[0],
                tuFullDistortion[0],
                tuFullDistortion[0],
                yCountNonZeroCoeffsTemp,
                0,
                0,
                COMPONENT_LUMA,
                asm_type);

            tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
            tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

            int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
            tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
            tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
#if BUG_FIX
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = tx_type;
#endif
            //LUMA-ONLY
            Av1TuEstimateCoeffBits(
                picture_control_set_ptr,
                candidateBuffer,
                context_ptr->cu_ptr,
                tuOriginIndex,
                0,
                context_ptr->coeff_est_entropy_coder_ptr,
                candidateBuffer->residualQuantCoeffPtr,
                yCountNonZeroCoeffsTemp,
                0,
                0,
                &yTuCoeffBits,
                &yTuCoeffBits,
                &yTuCoeffBits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                COMPONENT_LUMA,
                asm_type);

            Av1TuCalcCostLuma(
                context_ptr->cu_ptr->luma_txb_skip_context,
                candidateBuffer->candidate_ptr,
                txb_itr,
                context_ptr->blk_geom->txsize[txb_itr],
                yCountNonZeroCoeffsTemp,
                tuFullDistortion[0],
                &yTuCoeffBits,
                &yFullCost,
                context_ptr->full_lambda);

        }

        if (yFullCost < bestFullCost) {
            bestFullCost = yFullCost;
#if BUG_FIX
            best_tx_type = tx_type;
#else
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = tx_type;
#endif
        }

        //if (cpi->sf.adaptive_txb_search_level) {
        //    if ((best_rd - (best_rd >> cpi->sf.adaptive_txb_search_level)) >
        //        ref_best_rd) {
        //        break;
        //    }
        //}
        //// Skip transform type search when we found the block has been quantized to
        //// all zero and at the same time, it has better rdcost than doing transform.
        //if (cpi->sf.tx_type_search.skip_tx_search && !best_eob) break;


    }
#if BUG_FIX
    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;
#endif
    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
#if TX_TYPE_FIX
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y];
#else
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV];
#endif
}

#if ENCDEC_TX_SEARCH
void encode_pass_tx_search(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,          
    EbPictureBufferDesc_t          *residual16bit,           
    EbPictureBufferDesc_t          *transform16bit,          
    EbPictureBufferDesc_t          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane){

    (void)dZoffset;
    (void)use_delta_qp;
    (void)cbQp;
    UNUSED(count_non_zero_coeffs);
    UNUSED(component_mask);

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t               qp = cu_ptr->qp;
    const uint32_t         scratchLumaOffset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t         coeff1dOffset = context_ptr->coded_area_sb;

    EbBool                 cleanSparseCoeffFlag = EB_FALSE;
    uint64_t               yTuCoeffBits;
    uint64_t               tuFullDistortion[3][DIST_CALC_TOTAL];
    const int32_t          is_inter = context_ptr->is_inter;
    uint64_t               bestFullCost = UINT64_MAX;
    uint64_t               yFullCost;
    uint32_t               yCountNonZeroCoeffsTemp;
    TxType                 txk_start = DCT_DCT;
    TxType                 txk_end = TX_TYPES;
    TxType                 tx_type;
    TxSize                 txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];
    const TxSetType        tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);

    TxType best_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;

#if TX_SEARCH_LEVELS
        if(picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;
#else
#if FAST_TX_SEARCH
#if ENCODER_MODE_CLEANUP
        if (picture_control_set_ptr->enc_mode == ENC_M1)
#endif
            if (!allowed_tx_set_a[txSize][tx_type]) continue;
#endif
#endif
        const int32_t eset = get_ext_tx_set(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
        // eset == 0 should correspond to a set with only DCT_DCT and there
        // is no need to send the tx_type
        if (eset <= 0) continue;
        if (av1_ext_tx_used[tx_set_type][tx_type] == 0) continue;

        context_ptr->three_quad_energy = 0;

        yTuCoeffBits = 0;

        Av1EstimateTransform(
            ((int16_t*)residual16bit->bufferY) + scratchLumaOffset,
            residual16bit->strideY,
            ((tran_low_t*)transform16bit->bufferY) + coeff1dOffset,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            tx_type,
            asm_type,
            PLANE_TYPE_Y,
            context_ptr->trans_coeff_shape_luma);

        Av1QuantizeInvQuantize(
            sb_ptr->picture_control_set_ptr,
            ((tran_low_t*)transform16bit->bufferY) + coeff1dOffset,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->bufferY) + coeff1dOffset,
            ((int32_t*)inverse_quant_buffer->bufferY) + coeff1dOffset,
            qp,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &eob[0],
            candidate_plane[0],
            asm_type,
            &yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
            BIT_INCREMENT_8BIT,
            tx_type,
            cleanSparseCoeffFlag);

        //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
        if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {
            continue;
        }


        // LUMA DISTORTION
        PictureFullDistortion32Bits(
            transform16bit,
            coeff1dOffset,
            0,
            inverse_quant_buffer,
            coeff1dOffset,
            0,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv,
            tuFullDistortion[0],
            tuFullDistortion[0],
            tuFullDistortion[0],
            yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
            asm_type);

        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
        txb_ptr->transform_type[PLANE_TYPE_Y] = tx_type;

        //LUMA-ONLY

        ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[context_ptr->md_context->buffer_depth_index_start[0]]);
        ModeDecisionCandidateBuffer_t          *candidateBuffer;

        // Set the Candidate Buffer
        candidateBuffer = candidate_buffer_ptr_array[0];
        // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
        EntropyCoder_t  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;
        candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
        candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

        const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

        Av1TuEstimateCoeffBits(
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            coeff1dOffset,
            0,
            coeff_est_entropy_coder_ptr,
            coeffSamplesTB,
            yCountNonZeroCoeffsTemp,
            0,
            0,
            &yTuCoeffBits,
            &yTuCoeffBits,
            &yTuCoeffBits,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            COMPONENT_LUMA,
            asm_type);

        Av1TuCalcCostLuma(
            context_ptr->cu_ptr->luma_txb_skip_context,
            candidateBuffer->candidate_ptr,
            context_ptr->txb_itr,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            yCountNonZeroCoeffsTemp,
            tuFullDistortion[0],
            &yTuCoeffBits,
            &yFullCost,
            context_ptr->full_lambda);

        if (yFullCost < bestFullCost) {
            bestFullCost = yFullCost;
            best_tx_type = tx_type;
        }


    }

    txb_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;

    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        txb_ptr->transform_type[PLANE_TYPE_UV] = txb_ptr->transform_type[PLANE_TYPE_Y];

}

void encode_pass_tx_search_hbd(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,
    EbPictureBufferDesc_t          *residual16bit,
    EbPictureBufferDesc_t          *transform16bit,
    EbPictureBufferDesc_t          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane){

    (void)dZoffset;
    (void)use_delta_qp;
    (void)cbQp;
    UNUSED(component_mask);
    UNUSED(count_non_zero_coeffs);

    CodingUnit_t    *cu_ptr               = context_ptr->cu_ptr;
    TransformUnit_t *txb_ptr              = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t         qp                   = cu_ptr->qp;
    const uint32_t   scratchLumaOffset    = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t   coeff1dOffset        = context_ptr->coded_area_sb;
    EbBool           cleanSparseCoeffFlag = EB_FALSE;

    //Update QP for Quant
    qp += QP_BD_OFFSET;
    uint64_t                    yTuCoeffBits;
    uint64_t                    tuFullDistortion[3][DIST_CALC_TOTAL];
    const int32_t               is_inter = context_ptr->is_inter;
    uint64_t                    bestFullCost = UINT64_MAX;
    uint64_t                    yFullCost;
    uint32_t                    yCountNonZeroCoeffsTemp;
    TxType                      txk_start = DCT_DCT;
    TxType                      txk_end = TX_TYPES;
    TxType                      tx_type;
    TxSize                      txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];
    const TxSetType             tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);

    TxType best_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        ////if (!allowed_tx_mask[tx_type]) continue;
#if TX_SEARCH_LEVELS
        if (picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;
#else
#if FAST_TX_SEARCH
#if ENCODER_MODE_CLEANUP
        if (picture_control_set_ptr->enc_mode == ENC_M1 )
#endif
            if (!allowed_tx_set_a[txSize][tx_type]) continue;
#endif
#endif
        const int32_t eset = get_ext_tx_set(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
        // eset == 0 should correspond to a set with only DCT_DCT and there
        // is no need to send the tx_type
        if (eset <= 0) continue;
        if (av1_ext_tx_used[tx_set_type][tx_type] == 0) continue;

        context_ptr->three_quad_energy = 0;


        yTuCoeffBits = 0;

        Av1EstimateTransform(
            ((int16_t*)residual16bit->bufferY) + scratchLumaOffset,
            residual16bit->strideY,
            ((tran_low_t*)transform16bit->bufferY) + coeff1dOffset,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_10BIT,
            tx_type,
            asm_type,
            PLANE_TYPE_Y,
            context_ptr->trans_coeff_shape_luma);

        Av1QuantizeInvQuantize(
            sb_ptr->picture_control_set_ptr,
            ((int32_t*)transform16bit->bufferY) + coeff1dOffset,
            NOT_USED_VALUE,
#if QT_10BIT_SUPPORT
            ((int32_t*)coeffSamplesTB->bufferY) + coeff1dOffset,
#else
            ((int32_t*)coeffSamplesTB->bufferY) + scratchLumaOffset,
#endif
            ((int32_t*)inverse_quant_buffer->bufferY) + coeff1dOffset,
            qp,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &eob[0],
            candidate_plane[0],
            asm_type,
            &yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
#if QT_10BIT_SUPPORT
            BIT_INCREMENT_10BIT,
#endif
            tx_type,
            cleanSparseCoeffFlag);

        //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
        if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {
            continue;
        }


        // LUMA DISTORTION
        PictureFullDistortion32Bits(
            transform16bit,
            coeff1dOffset,
            0,
            inverse_quant_buffer,
            coeff1dOffset,
            0,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv,
            tuFullDistortion[0],
            tuFullDistortion[0],
            tuFullDistortion[0],
            yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
            asm_type);

        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
        txb_ptr->transform_type[PLANE_TYPE_Y] = tx_type;

        //LUMA-ONLY

        ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[context_ptr->md_context->buffer_depth_index_start[0]]);
        ModeDecisionCandidateBuffer_t          *candidateBuffer;

        // Set the Candidate Buffer
        candidateBuffer = candidate_buffer_ptr_array[0];
        // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
        EntropyCoder_t  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;
        candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
        candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

        const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

        Av1TuEstimateCoeffBits(
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            coeff1dOffset,
            0,
            coeff_est_entropy_coder_ptr,
            coeffSamplesTB,
            yCountNonZeroCoeffsTemp,
            0,
            0,
            &yTuCoeffBits,
            &yTuCoeffBits,
            &yTuCoeffBits,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            COMPONENT_LUMA,
            asm_type);

        Av1TuCalcCostLuma(
            context_ptr->cu_ptr->luma_txb_skip_context,
            candidateBuffer->candidate_ptr,
            context_ptr->txb_itr,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            yCountNonZeroCoeffsTemp,
            tuFullDistortion[0],
            &yTuCoeffBits,
            &yFullCost,
            context_ptr->full_lambda);

        if (yFullCost < bestFullCost) {
            bestFullCost = yFullCost;
            best_tx_type = tx_type;
        }


    }


    txb_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;

    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        txb_ptr->transform_type[PLANE_TYPE_UV] = txb_ptr->transform_type[PLANE_TYPE_Y];

}
#endif
/****************************************
 ************  Full loop ****************
****************************************/
void FullLoop_R(
    LargestCodingUnit_t            *sb_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    EbPictureBufferDesc_t          *inputPicturePtr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint32_t                          component_mask,
    uint32_t                          cbQp,
    uint32_t                          crQp,
    uint32_t                          *cbCountNonZeroCoeffs,
    uint32_t                          *crCountNonZeroCoeffs)
{
    (void)sb_ptr;
    (void)crQp;
    (void)inputPicturePtr;
    int16_t                *chromaResidualPtr;
    uint32_t                 tuOriginIndex;
    UNUSED(tuOriginIndex);
    uint32_t                 tuCbOriginIndex;
    uint32_t                 tuCrOriginIndex;
    uint32_t                 tuCount;
    uint32_t                 txb_itr;
    uint32_t                 txb_origin_x;
    uint32_t                 txb_origin_y;

    // EbPictureBufferDesc_t         * tuTransCoeffTmpPtr;
     //EbPictureBufferDesc_t         * tuQuantCoeffTmpPtr;

    SequenceControlSet_t    *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    EbAsm     asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    EbBool cleanSparseCoeffFlag = EB_FALSE;
    context_ptr->three_quad_energy = 0;

    tuCount = context_ptr->blk_geom->txb_count;
    uint32_t  txb_1d_offset = 0;

    txb_itr = 0;
    do {

        txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];


        // NADER - TU
        tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->residualQuantCoeffPtr->strideY;
        tuCbOriginIndex = (((txb_origin_x >> 3) << 3) + (((txb_origin_y >> 3) << 3) * candidateBuffer->residualQuantCoeffPtr->strideCb)) >> 1;
        tuCrOriginIndex = (((txb_origin_x >> 3) << 3) + (((txb_origin_y >> 3) << 3) * candidateBuffer->residualQuantCoeffPtr->strideCr)) >> 1;

        //    This function replaces the previous Intra Chroma mode if the LM fast
            //    cost is better.
            //    *Note - this might require that we have inv transform in the loop
        EbPfMode    correctedPFMode = PF_OFF;

        if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
            // Configure the Chroma Residual Ptr

            chromaResidualPtr = //(candidateBuffer->candidate_ptr->type  == INTRA_MODE )?
                  //&(((int16_t*) candidateBuffer->intraChromaResidualPtr->bufferCb)[tuChromaOriginIndex]):
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCb)[tuCbOriginIndex]);


            // Cb Transform
            Av1EstimateTransform(
                chromaResidualPtr,
                candidateBuffer->residual_ptr->strideCb,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCb)[txb_1d_offset]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
                correctedPFMode);


#if SIMULATE_PF_N2     
            int32_t* transCoeffBuffer = &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCb)[txb_1d_offset]);
            int32_t tu_size;
            switch (context_ptr->blk_geom->txsize_uv[txb_itr]) {
            case TX_64X64:
                tu_size = 64;
                break;
            case TX_32X32:
                tu_size = 32;
                break;
            case TX_16X16:
                tu_size = 16;
                break;
            case TX_8X8:
                tu_size = 8;
                break;
            case TX_4X4:
                tu_size = 4;
                break;
            default: assert(0); break;
            }

            if ((tu_size == 32 || tu_size == 16 || tu_size == 8 || tu_size == 4) && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0) {
                for (int i = 0; i < (tu_size* tu_size); i++) {
                    if (i % tu_size >= (tu_size >> 1) || i / tu_size >= (tu_size >> 1)) {
                        transCoeffBuffer[i] = 0;
                    }
                }
            }
#endif
            Av1QuantizeInvQuantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCb)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCb)[txb_1d_offset]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCb)[txb_1d_offset]),
                cbQp,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &candidateBuffer->candidate_ptr->eob[1][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[1],
                asm_type,
                &(cbCountNonZeroCoeffs[txb_itr]),
                context_ptr->pf_md_mode,
                0,
                COMPONENT_CHROMA_CB,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_8BIT,
#endif
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                cleanSparseCoeffFlag);
            candidateBuffer->candidate_ptr->quantized_dc[1] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCb)[txb_1d_offset]);
        }


        if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
            // Configure the Chroma Residual Ptr

            chromaResidualPtr = //(candidateBuffer->candidate_ptr->type  == INTRA_MODE )?
                //&(((int16_t*) candidateBuffer->intraChromaResidualPtr->bufferCr)[tuChromaOriginIndex]):
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCr)[tuCrOriginIndex]);

            // Cr Transform
            Av1EstimateTransform(
                chromaResidualPtr,
                candidateBuffer->residual_ptr->strideCr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCr)[txb_1d_offset]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
                correctedPFMode);

#if SIMULATE_PF_N2     
            int32_t* transCoeffBuffer = &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCr)[txb_1d_offset]);
            int32_t tu_size;
            switch (context_ptr->blk_geom->txsize_uv[txb_itr]) {
            case TX_64X64:
                tu_size = 64;
                break;
            case TX_32X32:
                tu_size = 32;
                break;
            case TX_16X16:
                tu_size = 16;
                break;
            case TX_8X8:
                tu_size = 8;
                break;
            case TX_4X4:
                tu_size = 4;
                break;
            default: assert(0); break;
            }

            if ((tu_size == 32 || tu_size == 16 || tu_size == 8 || tu_size == 4) && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0) {
                    
                    for (int i = 0; i < (tu_size* tu_size); i++) {
                    if (i % tu_size >= (tu_size >> 1) || i / tu_size >= (tu_size >> 1)) {
                        transCoeffBuffer[i] = 0;
                    }
                }
            }
#endif
            Av1QuantizeInvQuantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCr)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCr)[txb_1d_offset]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCr)[txb_1d_offset]),
                cbQp,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &candidateBuffer->candidate_ptr->eob[2][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[2],
                asm_type,
                &(crCountNonZeroCoeffs[txb_itr]),
                context_ptr->pf_md_mode,
                0,
                COMPONENT_CHROMA_CR,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_8BIT,
#endif
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                cleanSparseCoeffFlag);
            candidateBuffer->candidate_ptr->quantized_dc[2] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCr)[txb_1d_offset]);
        }

        txb_1d_offset += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];


        ++txb_itr;

    } while (txb_itr < tuCount);

}

//****************************************
// ************ CuFullDistortionFastTuMode ****************
//****************************************/
void CuFullDistortionFastTuMode_R(
    LargestCodingUnit_t            *sb_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t            *context_ptr,
    ModeDecisionCandidate_t           *candidate_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint64_t                          cbFullDistortion[DIST_CALC_TOTAL],
    uint64_t                          crFullDistortion[DIST_CALC_TOTAL],
    uint32_t                          count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU],
    COMPONENT_TYPE                  componentType,
    uint64_t                         *cb_coeff_bits,
    uint64_t                         *cr_coeff_bits,
    EbAsm                            asm_type)
{
    (void)sb_ptr;

    uint64_t                          yTuCoeffBits;
    uint64_t                          cbTuCoeffBits;
    uint64_t                          crTuCoeffBits;
    uint32_t                          tuOriginIndex;
    uint32_t                          txb_origin_x;
    uint32_t                          txb_origin_y;
    uint32_t                          currentTuIndex;
    int32_t                          chromaShift;
    uint32_t                          tuChromaOriginIndex;
    uint64_t                          tuFullDistortion[3][DIST_CALC_TOTAL];
    EbPictureBufferDesc_t          *transform_buffer;
    uint32_t                          tuTotalCount;
    uint32_t                          txb_itr = 0;
    //    SequenceControlSet_t           *sequence_control_set_ptr=((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr);

    tuTotalCount = context_ptr->blk_geom->txb_count;
    currentTuIndex = 0;
    transform_buffer = context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr;


    uint32_t  txb_1d_offset = 0;
    candidate_ptr->u_has_coeff = 0;
    candidate_ptr->v_has_coeff = 0;

    do {

        txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];



        tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->residualQuantCoeffPtr->strideY;
        tuChromaOriginIndex = txb_1d_offset;
        // Reset the Bit Costs
        yTuCoeffBits = 0;
        cbTuCoeffBits = 0;
        crTuCoeffBits = 0;

        if (componentType == COMPONENT_CHROMA_CB || componentType == COMPONENT_CHROMA_CR || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {
            uint32_t countNonZeroCoeffsAll[3];
            countNonZeroCoeffsAll[0] = count_non_zero_coeffs[0][currentTuIndex];
            countNonZeroCoeffsAll[1] = count_non_zero_coeffs[1][currentTuIndex];
            countNonZeroCoeffsAll[2] = count_non_zero_coeffs[2][currentTuIndex];


            // *Full Distortion (SSE)
            // *Note - there are known issues with how this distortion metric is currently
            //    calculated.  The amount of scaling between the two arrays is not
            //    equivalent.


            PictureFullDistortion32Bits(
                transform_buffer,
                NOT_USED_VALUE,
                tuChromaOriginIndex,
                candidateBuffer->reconCoeffPtr,
                NOT_USED_VALUE,
                tuChromaOriginIndex,
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                tuFullDistortion[0],
                tuFullDistortion[1],
                tuFullDistortion[2],
                countNonZeroCoeffsAll[0],
                countNonZeroCoeffsAll[1],
                countNonZeroCoeffsAll[2],
                componentType,
                asm_type);

            TxSize    txSize = context_ptr->blk_geom->txsize_uv[txb_itr];
            chromaShift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
            tuFullDistortion[1][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_RESIDUAL], chromaShift);
            tuFullDistortion[1][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_PREDICTION], chromaShift);
            tuFullDistortion[2][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_RESIDUAL], chromaShift);
            tuFullDistortion[2][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_PREDICTION], chromaShift);





            //CHROMA-ONLY
            Av1TuEstimateCoeffBits(
                picture_control_set_ptr,
                candidateBuffer,
                context_ptr->cu_ptr,
                tuOriginIndex,
                tuChromaOriginIndex,
                context_ptr->coeff_est_entropy_coder_ptr,
                candidateBuffer->residualQuantCoeffPtr,
                count_non_zero_coeffs[0][currentTuIndex],
                count_non_zero_coeffs[1][currentTuIndex],
                count_non_zero_coeffs[2][currentTuIndex],
                &yTuCoeffBits,
                &cbTuCoeffBits,
                &crTuCoeffBits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],

                componentType,
                asm_type);


            // OMK Useless ? We don't calculate Chroma CBF here
            Av1TuCalcCost(
                candidate_ptr,
                context_ptr->cu_ptr->luma_txb_skip_context,
                currentTuIndex,
                count_non_zero_coeffs[0][currentTuIndex],
                count_non_zero_coeffs[1][currentTuIndex],
                count_non_zero_coeffs[2][currentTuIndex],
                tuFullDistortion[0],
                tuFullDistortion[1],
                tuFullDistortion[2],
                componentType,
                &yTuCoeffBits,
                &cbTuCoeffBits,
                &crTuCoeffBits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->full_lambda);



            *cb_coeff_bits += cbTuCoeffBits;
            *cr_coeff_bits += crTuCoeffBits;
            cbFullDistortion[DIST_CALC_RESIDUAL] += tuFullDistortion[1][DIST_CALC_RESIDUAL];
            crFullDistortion[DIST_CALC_RESIDUAL] += tuFullDistortion[2][DIST_CALC_RESIDUAL];
            cbFullDistortion[DIST_CALC_PREDICTION] += tuFullDistortion[1][DIST_CALC_PREDICTION];
            crFullDistortion[DIST_CALC_PREDICTION] += tuFullDistortion[2][DIST_CALC_PREDICTION];

        }


        txb_1d_offset += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];
        currentTuIndex++;


        ++txb_itr;


    } while (txb_itr < tuTotalCount);
}


void  d1_non_square_block_decision(
    ModeDecisionContext_t               *context_ptr
)
{
    //compute total cost for the whole block partition
    uint64_t tot_cost = 0;
    uint32_t first_blk_idx = context_ptr->cu_ptr->mds_idx - (context_ptr->blk_geom->totns - 1);//index of first block in this partition
    uint32_t blk_it;
    for (blk_it = 0; blk_it < context_ptr->blk_geom->totns; blk_it++)
    {
        tot_cost += context_ptr->md_local_cu_unit[first_blk_idx + blk_it].cost;
    }

    if (context_ptr->blk_geom->shape == PART_N || tot_cost < context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].cost)
    {
        //store best partition cost in parent square
        context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].cost = tot_cost;
        context_ptr->md_cu_arr_nsq[context_ptr->blk_geom->sqi_mds].part = from_shape_to_part[context_ptr->blk_geom->shape];
        context_ptr->md_cu_arr_nsq[context_ptr->blk_geom->sqi_mds].best_d1_blk = first_blk_idx;
    }


}


#if FIX_INTER_DEPTH
/// compute the cost of curr depth, and the depth above
void   compute_depth_costs(
    ModeDecisionContext_t    *context_ptr,
    SequenceControlSet_t     *sequence_control_set_ptr,
    uint32_t                  curr_depth_mds,
    uint32_t                  above_depth_mds,
    uint32_t                  step,
#if FIX_47
    uint64_t                 *above_depth_cost,
    uint64_t                 *curr_depth_cost)
#else
    uint64_t                 *depthNCost,
    uint64_t                 *depthNPlusOneCost)
#endif
{
#if FIX_47
    uint64_t       above_non_split_rate = 0;
    uint64_t       above_split_rate = 0;
#else
    uint64_t       depthNRate = 0;
    uint64_t       depthNPlusOneRate = 0;
#endif


    /*
    ___________
    |     |     |
    |blk0 |blk1 |
    |-----|-----|
    |blk2 |blk3 |
    |_____|_____|
    */
    // current depth blocks
    uint32_t       curr_depth_blk0_mds = curr_depth_mds - 3 * step;
    uint32_t       curr_depth_blk1_mds = curr_depth_mds - 2 * step;
    uint32_t       curr_depth_blk2_mds = curr_depth_mds - 1 * step;
    uint32_t       curr_depth_blk3_mds = curr_depth_mds;

    // Rate of not spliting the current depth (Depth != 4) in case the children were omitted by MDC
    uint64_t       curr_non_split_rate_blk0 = 0;
    uint64_t       curr_non_split_rate_blk1 = 0;
    uint64_t       curr_non_split_rate_blk2 = 0;
    uint64_t       curr_non_split_rate_blk3 = 0;


    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_mode = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_depth = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_mode = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].top_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_depth = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].top_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_partition = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_partition;
    context_ptr->md_local_cu_unit[above_depth_mds].above_neighbor_partition = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].above_neighbor_partition;

    // Compute above depth  cost
    if (context_ptr->md_local_cu_unit[above_depth_mds].tested_cu_flag == EB_TRUE)
    {
        Av1SplitFlagRate(
            sequence_control_set_ptr,
            context_ptr,
            &context_ptr->md_cu_arr_nsq[above_depth_mds],
            0,
            PARTITION_NONE,//shouldn't this be final partition for above depth?
            &above_non_split_rate,
            context_ptr->full_lambda,
            context_ptr->md_rate_estimation_ptr,
            sequence_control_set_ptr->max_sb_depth);

        *above_depth_cost = context_ptr->md_local_cu_unit[above_depth_mds].cost + above_non_split_rate;
    }
    else {
        *above_depth_cost = MAX_MODE_COST;
    }

    // Compute curr depth  cost
    Av1SplitFlagRate(
        sequence_control_set_ptr,
        context_ptr,
        &context_ptr->md_cu_arr_nsq[above_depth_mds],
        0,
        PARTITION_SPLIT,
        &above_split_rate,
        context_ptr->full_lambda,
        context_ptr->md_rate_estimation_ptr,
        sequence_control_set_ptr->max_sb_depth);

    if (context_ptr->blk_geom->bsize > BLOCK_4X4) {

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk0_mds].mdc_split_flag == 0)
            Av1SplitFlagRate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk0_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk0,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk1_mds].mdc_split_flag == 0)
            Av1SplitFlagRate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk1_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk1,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk2_mds].mdc_split_flag == 0)
            Av1SplitFlagRate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk2_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk2,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk3_mds].mdc_split_flag == 0)
            Av1SplitFlagRate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk3_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk3,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);
    }
    //curr_non_split_rate_344 = splitflag_mdc_344 || 4x4 ? 0 : compute; 


    *curr_depth_cost =
        context_ptr->md_local_cu_unit[curr_depth_mds].cost + curr_non_split_rate_blk3 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 1 * step].cost + curr_non_split_rate_blk2 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 2 * step].cost + curr_non_split_rate_blk1 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 3 * step].cost + curr_non_split_rate_blk0 +
        above_split_rate;

}
#else
/// compute the cost of curr depth, and the depth above
void   compute_depth_costs(
    ModeDecisionContext_t    *context_ptr,
    SequenceControlSet_t     *sequence_control_set_ptr,
    uint32_t                  curr_depth_mds,
    uint32_t                  above_depth_mds,
    uint32_t                  step,
#if FIX_47
    uint64_t                 *above_depth_cost,
    uint64_t                 *curr_depth_cost)
#else
    uint64_t                 *depthNCost,
    uint64_t                 *depthNPlusOneCost)
#endif
{
#if FIX_47
    uint64_t       above_rate = 0;
    uint64_t       curr_rate = 0;
#else
    uint64_t       depthNRate = 0;
    uint64_t       depthNPlusOneRate = 0;
#endif
    uint32_t     top_left_idx_mds = curr_depth_mds - 3 * step;

    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_mode = context_ptr->md_local_cu_unit[top_left_idx_mds].left_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_depth = context_ptr->md_local_cu_unit[top_left_idx_mds].left_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_mode = context_ptr->md_local_cu_unit[top_left_idx_mds].top_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_depth = context_ptr->md_local_cu_unit[top_left_idx_mds].top_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_partition = context_ptr->md_local_cu_unit[top_left_idx_mds].left_neighbor_partition;
    context_ptr->md_local_cu_unit[above_depth_mds].above_neighbor_partition = context_ptr->md_local_cu_unit[top_left_idx_mds].above_neighbor_partition;
#if FIX_47
    // Compute above depth  cost
    if (context_ptr->md_local_cu_unit[above_depth_mds].tested_cu_flag == EB_TRUE)
    {
        Av1SplitFlagRate(
            sequence_control_set_ptr,
            context_ptr,
            &context_ptr->md_cu_arr_nsq[above_depth_mds],
            0,
            PARTITION_NONE,//shouldn't this be final partition for above depth?
            &above_rate,
            context_ptr->full_lambda,
            context_ptr->md_rate_estimation_ptr,
            sequence_control_set_ptr->max_sb_depth);

        *above_depth_cost = context_ptr->md_local_cu_unit[above_depth_mds].cost + above_rate;
    }
    else {
        *above_depth_cost = MAX_MODE_COST;
    }

    // Compute curr depth  cost
    Av1SplitFlagRate(
        sequence_control_set_ptr,
        context_ptr,
        &context_ptr->md_cu_arr_nsq[above_depth_mds],
        0,
        PARTITION_SPLIT,
        &curr_rate,
        context_ptr->full_lambda,
        context_ptr->md_rate_estimation_ptr,
        sequence_control_set_ptr->max_sb_depth);

    *curr_depth_cost =
        context_ptr->md_local_cu_unit[curr_depth_mds].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 1 * step].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 2 * step].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 3 * step].cost +
        curr_rate;
#else  
    // Compute depth N cost
    Av1SplitFlagRate(
        sequence_control_set_ptr,
        context_ptr,
        &context_ptr->md_cu_arr_nsq[above_depth_mds],
        0,
        PARTITION_NONE,
        &depthNRate,
        context_ptr->full_lambda,
        context_ptr->md_rate_estimation_ptr,
        sequence_control_set_ptr->max_sb_depth);

    if (context_ptr->md_local_cu_unit[above_depth_mds].tested_cu_flag == EB_FALSE)
        context_ptr->md_local_cu_unit[above_depth_mds].cost = MAX_MODE_COST;

    *depthNCost = context_ptr->md_local_cu_unit[above_depth_mds].cost + depthNRate;

    // Compute depth N+1 cost
    Av1SplitFlagRate(
        sequence_control_set_ptr,
        context_ptr,
        &context_ptr->md_cu_arr_nsq[above_depth_mds],
        0,
        PARTITION_SPLIT,
        &depthNPlusOneRate,
        context_ptr->full_lambda,
        context_ptr->md_rate_estimation_ptr,
        sequence_control_set_ptr->max_sb_depth);


    *depthNPlusOneCost =
        context_ptr->md_local_cu_unit[curr_depth_mds].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 1 * step].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 2 * step].cost +
        context_ptr->md_local_cu_unit[curr_depth_mds - 3 * step].cost +
        depthNPlusOneRate;
#endif

}
#endif
uint32_t d2_inter_depth_block_decision(
    ModeDecisionContext_t          *context_ptr,
    uint32_t                        blk_mds,
    LargestCodingUnit_t            *tbPtr,
    uint32_t                          lcuAddr,
    uint32_t                          tbOriginX,
    uint32_t                          tbOriginY,
    uint64_t                          full_lambda,
    MdRateEstimationContext_t      *md_rate_estimation_ptr,
    PictureControlSet_t            *picture_control_set_ptr)
{
    UNUSED(tbPtr);
    UNUSED(lcuAddr);
    UNUSED(tbOriginX);
    UNUSED(tbOriginY);
    UNUSED(full_lambda);
    UNUSED(md_rate_estimation_ptr);

    uint32_t                  lastCuIndex, d0_idx_mds, d1_idx_mds, d2_idx_mds, top_left_idx_mds;
    UNUSED(top_left_idx_mds);
    UNUSED(d2_idx_mds);
    UNUSED(d1_idx_mds);
    UNUSED(d0_idx_mds);
    uint64_t                    parent_depth_cost = 0, current_depth_cost = 0;
    SequenceControlSet_t     *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    EbBool                    lastDepthFlag;
    const BlockGeom          * blk_geom;

    lastDepthFlag = context_ptr->md_cu_arr_nsq[blk_mds].split_flag == EB_FALSE ? EB_TRUE : EB_FALSE;
    d1_idx_mds = blk_mds;
    d2_idx_mds = blk_mds;
    lastCuIndex = blk_mds;
    blk_geom = Get_blk_geom_mds(blk_mds);
    uint32_t    parent_depth_idx_mds = blk_mds;
    uint32_t    current_depth_idx_mds = blk_mds;

    if (lastDepthFlag) {
        while (blk_geom->is_last_quadrant) {

            //get parent idx
            parent_depth_idx_mds = current_depth_idx_mds - parent_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];
#if FIX_DEBUG_CRASH
            if (picture_control_set_ptr->slice_type == I_SLICE && parent_depth_idx_mds == 0) {
                parent_depth_cost = MAX_MODE_COST;
            }
            else {
#endif
                compute_depth_costs(context_ptr, sequence_control_set_ptr, current_depth_idx_mds, parent_depth_idx_mds, ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth], &parent_depth_cost, &current_depth_cost);
#if FIX_DEBUG_CRASH
            }
#endif
            if (parent_depth_cost <= current_depth_cost) {
                context_ptr->md_cu_arr_nsq[parent_depth_idx_mds].split_flag = EB_FALSE;
                context_ptr->md_local_cu_unit[parent_depth_idx_mds].cost = parent_depth_cost;
                lastCuIndex = parent_depth_idx_mds;
            }
            else {
                context_ptr->md_local_cu_unit[parent_depth_idx_mds].cost = current_depth_cost;
                context_ptr->md_cu_arr_nsq[parent_depth_idx_mds].part = PARTITION_SPLIT;
            }

            //setup next parent inter depth
            blk_geom = Get_blk_geom_mds(parent_depth_idx_mds);
            current_depth_idx_mds = parent_depth_idx_mds;
        }
    }


    return lastCuIndex;
}



