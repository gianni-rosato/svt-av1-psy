/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "EbDefinitions.h"
#include "EbFullLoop.h"
#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"

#ifdef __GNUC__
#define LIKELY(v) __builtin_expect(v, 1)
#define UNLIKELY(v) __builtin_expect(v, 0)
#else
#define LIKELY(v) (v)
#define UNLIKELY(v) (v)
#endif
void residual_kernel(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *pred,
                     uint32_t pred_offset, uint32_t pred_stride, int16_t *residual,
                     uint32_t residual_offset, uint32_t residual_stride, Bool hbd,
                     uint32_t area_width, uint32_t area_height);

void svt_aom_quantize_b_c_ii(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,
                             const int16_t *round_ptr, const int16_t *quant_ptr,
                             const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                             TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
                             const int16_t *scan, const int16_t *iscan, const QmVal *qm_ptr,
                             const QmVal *iqm_ptr, const int32_t log_scale) {
    const int32_t zbins[2]       = {ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
                                    ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale)};
    const int32_t nzbins[2]      = {zbins[0] * -1, zbins[1] * -1};
    intptr_t      non_zero_count = n_coeffs, eob = -1;
    (void)iscan;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    // Pre-scan pass
    for (intptr_t i = n_coeffs - 1; i >= 0; i--) {
        const int32_t rc    = scan[i];
        const QmVal   wt    = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
        const int32_t coeff = coeff_ptr[rc] * wt;

        if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
            coeff > (nzbins[rc != 0] * (1 << AOM_QM_BITS)))
            non_zero_count--;
        else
            break;
    }

    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (intptr_t i = 0; i < non_zero_count; i++) {
        const int32_t rc         = scan[i];
        const int32_t coeff      = coeff_ptr[rc];
        const int     coeff_sign = coeff < 0 ? -1 : 0;
        const int32_t abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;

        const QmVal wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
        if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
            int64_t tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale),
                                INT16_MIN,
                                INT16_MAX);
            tmp *= wt;
            int32_t tmp32         = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
                                       quant_shift_ptr[rc != 0]) >>
                                      (16 - log_scale + AOM_QM_BITS)); // quantization
            qcoeff_ptr[rc]        = (tmp32 ^ coeff_sign) - coeff_sign;
            const int32_t iwt     = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;
            const TranLow abs_dqcoeff = (tmp32 * dequant) >> log_scale;
            dqcoeff_ptr[rc]           = (TranLow)((abs_dqcoeff ^ coeff_sign) - coeff_sign);

            if (tmp32)
                eob = i;
        }
    }
    *eob_ptr = (uint16_t)(eob + 1);
}

void svt_aom_quantize_b_c(const TranLow *coeff_ptr, int32_t stride, int32_t width, int32_t height,
                          intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr,
                          const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                          TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr,
                          uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan,
                          const QmVal *qm_ptr, const QmVal *iqm_ptr, const int32_t log_scale) {
    const int32_t zbins[2]       = {ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
                                    ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale)};
    const int32_t nzbins[2]      = {zbins[0] * -1, zbins[1] * -1};
    intptr_t      non_zero_count = n_coeffs, eob = -1;
    (void)iscan;

    // Nader quantisation
    for (int32_t x = 0; x < height; x++) {
        memset(qcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*qcoeff_ptr));
        memset(dqcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*dqcoeff_ptr));
    }

    // Pre-scan pass
    for (intptr_t i = n_coeffs - 1; i >= 0; i--) {
        const int32_t map_rc = scan[i];

        const int32_t rc = ((map_rc / MIN(32, height)) * stride) + (map_rc % MIN(32, width));

        const QmVal   wt    = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
        const int32_t coeff = coeff_ptr[rc] * wt;

        ////if (map_rc != NewTab[rc])
        //SVT_LOG("%d\n", coeff);

        if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
            coeff > (nzbins[rc != 0] * (1 << AOM_QM_BITS)))
            non_zero_count--;
        else
            break;
    }
    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (intptr_t i = 0; i < non_zero_count; i++) {
        const int32_t map_rc = scan[i];

        const int32_t rc    = ((map_rc / MIN(32, height)) * stride) + (map_rc % MIN(32, width));
        const int32_t coeff = coeff_ptr[rc];
        const int     coeff_sign = coeff < 0 ? -1 : 0;
        const int32_t abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;

        const QmVal wt = qm_ptr != NULL ? qm_ptr[map_rc] : (1 << AOM_QM_BITS);

        if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
            int64_t tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale),
                                INT16_MIN,
                                INT16_MAX);

            tmp *= wt;

            int32_t tmp32 = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
                                       quant_shift_ptr[rc != 0]) >>
                                      (16 - log_scale + AOM_QM_BITS)); // quantization

            qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;

            const int32_t iwt = iqm_ptr != NULL ? iqm_ptr[map_rc] : (1 << AOM_QM_BITS);

            const int32_t dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;

            dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant / (1 << log_scale);

            if (tmp32)
                eob = i;
        }
    }

    *eob_ptr = (uint16_t)(eob + 1);
}
void svt_aom_highbd_quantize_b_c(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                 const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                 TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                 const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan,
                                 const int16_t *iscan, const QmVal *qm_ptr, const QmVal *iqm_ptr,
                                 const int32_t log_scale) {
    intptr_t eob = -1;
    (void)iscan;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    const int32_t zbins[2]  = {ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
                               ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale)};
    const int32_t nzbins[2] = {zbins[0] * -1, zbins[1] * -1};
    intptr_t      idx_arr[4096];
    int           idx = 0;
    // Pre-scan pass
    for (intptr_t i = 0; i < n_coeffs; i++) {
        const int32_t rc    = scan[i];
        const QmVal   wt    = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
        const int32_t coeff = coeff_ptr[rc] * wt;

        // If the coefficient is out of the base ZBIN range, keep it for
        // quantization.
        if (coeff >= (zbins[rc != 0] * (1 << AOM_QM_BITS)) ||
            coeff <= (nzbins[rc != 0] * (1 << AOM_QM_BITS)))
            idx_arr[idx++] = i;
    }

    // Quantization pass: only process the coefficients selected in
    // pre-scan pass. Note: idx can be zero.
    for (int i = 0; i < idx; i++) {
        const int32_t rc         = scan[idx_arr[i]];
        const int32_t coeff      = coeff_ptr[rc];
        const int     coeff_sign = coeff < 0 ? -1 : 0;
        const QmVal   wt         = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
        const QmVal   iwt        = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
        const int32_t abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;
        const int64_t tmp1       = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale);
        const int64_t tmpw       = tmp1 * wt;
        const int64_t tmp2       = ((tmpw * quant_ptr[rc != 0]) >> 16) + tmpw;
        const int32_t abs_qcoeff = (int32_t)((tmp2 * quant_shift_ptr[rc != 0]) >>
                                             (16 - log_scale + AOM_QM_BITS));
        qcoeff_ptr[rc]           = (TranLow)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        int32_t dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
        const TranLow abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
        dqcoeff_ptr[rc]           = (TranLow)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
        if (abs_qcoeff)
            eob = idx_arr[i];
    }

    *eob_ptr = (uint16_t)(eob + 1);
}

void svt_av1_highbd_quantize_b_facade(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                      const MacroblockPlane *p, TranLow *qcoeff_ptr,
                                      TranLow *dqcoeff_ptr, uint16_t *eob_ptr, const ScanOrder *sc,
                                      const QuantParam *qparam) {
    const QmVal *qm_ptr  = qparam->qmatrix;
    const QmVal *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr || iqm_ptr) {
        svt_av1_highbd_quantize_b_qm(coeff_ptr,
                                     n_coeffs,
                                     p->zbin_qtx,
                                     p->round_qtx,
                                     p->quant_qtx,
                                     p->quant_shift_qtx,
                                     qcoeff_ptr,
                                     dqcoeff_ptr,
                                     p->dequant_qtx,
                                     eob_ptr,
                                     sc->scan,
                                     sc->iscan,
                                     qm_ptr,
                                     iqm_ptr,
                                     qparam->log_scale);
    } else {
        svt_aom_highbd_quantize_b(coeff_ptr,
                                  n_coeffs,
                                  p->zbin_qtx,
                                  p->round_qtx,
                                  p->quant_qtx,
                                  p->quant_shift_qtx,
                                  qcoeff_ptr,
                                  dqcoeff_ptr,
                                  p->dequant_qtx,
                                  eob_ptr,
                                  sc->scan,
                                  sc->iscan,
                                  NULL,
                                  NULL,
                                  qparam->log_scale);
    }
    assert(qparam->log_scale <= 2);
}

void av1_quantize_b_facade_ii(const TranLow *coeff_ptr, intptr_t n_coeffs, const MacroblockPlane *p,
                              TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, uint16_t *eob_ptr,
                              const ScanOrder *sc, const QuantParam *qparam) {
    const QmVal *qm_ptr  = qparam->qmatrix;
    const QmVal *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr || iqm_ptr) {
        svt_av1_quantize_b_qm(coeff_ptr,
                              n_coeffs,
                              p->zbin_qtx,
                              p->round_qtx,
                              p->quant_qtx,
                              p->quant_shift_qtx,
                              qcoeff_ptr,
                              dqcoeff_ptr,
                              p->dequant_qtx,
                              eob_ptr,
                              sc->scan,
                              sc->iscan,
                              qm_ptr,
                              iqm_ptr,
                              qparam->log_scale);
    } else {
        svt_aom_quantize_b(coeff_ptr,
                           n_coeffs,
                           p->zbin_qtx,
                           p->round_qtx,
                           p->quant_qtx,
                           p->quant_shift_qtx,
                           qcoeff_ptr,
                           dqcoeff_ptr,
                           p->dequant_qtx,
                           eob_ptr,
                           sc->scan,
                           sc->iscan,
                           NULL,
                           NULL,
                           qparam->log_scale);
    }
    assert(qparam->log_scale <= 2);
}

static void quantize_fp_helper_c(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                 const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                 TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                 const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan,
                                 const int16_t *iscan, const QmVal *qm_ptr, const QmVal *iqm_ptr,
                                 int log_scale) {
    int       i, eob = -1;
    const int rounding[2] = {ROUND_POWER_OF_TWO(round_ptr[0], log_scale),
                             ROUND_POWER_OF_TWO(round_ptr[1], log_scale)};
    (void)zbin_ptr;
    (void)quant_shift_ptr;
    (void)iscan;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (qm_ptr == NULL && iqm_ptr == NULL) {
        for (i = 0; i < n_coeffs; i++) {
            const int     rc         = scan[i];
            const int32_t thresh     = (int32_t)(dequant_ptr[rc != 0]);
            const int     coeff      = coeff_ptr[rc];
            const int     coeff_sign = coeff < 0 ? -1 : 0;
            int64_t       abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;
            int           tmp32      = 0;
            if ((abs_coeff << (1 + log_scale)) >= thresh) {
                abs_coeff = clamp64(abs_coeff + rounding[rc != 0], INT16_MIN, INT16_MAX);
                tmp32     = (int)((abs_coeff * quant_ptr[rc != 0]) >> (16 - log_scale));
                if (tmp32) {
                    qcoeff_ptr[rc]            = (tmp32 ^ coeff_sign) - coeff_sign;
                    const TranLow abs_dqcoeff = (tmp32 * dequant_ptr[rc != 0]) >> log_scale;
                    dqcoeff_ptr[rc]           = (abs_dqcoeff ^ coeff_sign) - coeff_sign;
                }
            }
            if (tmp32)
                eob = i;
        }
    } else {
        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < n_coeffs; i++) {
            const int   rc      = scan[i];
            const int   coeff   = coeff_ptr[rc];
            const QmVal wt      = qm_ptr ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const QmVal iwt     = iqm_ptr ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
            const int   dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;
            const int coeff_sign = coeff < 0 ? -1 : 0;
            int64_t   abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;
            int       tmp32      = 0;
            if (abs_coeff * wt >= (dequant_ptr[rc != 0] << (AOM_QM_BITS - (1 + log_scale)))) {
                abs_coeff += rounding[rc != 0];
                abs_coeff                 = clamp64(abs_coeff, INT16_MIN, INT16_MAX);
                tmp32                     = (int)((abs_coeff * wt * quant_ptr[rc != 0]) >>
                              (16 - log_scale + AOM_QM_BITS));
                qcoeff_ptr[rc]            = (tmp32 ^ coeff_sign) - coeff_sign;
                const TranLow abs_dqcoeff = (tmp32 * dequant) >> log_scale;
                dqcoeff_ptr[rc]           = (abs_dqcoeff ^ coeff_sign) - coeff_sign;
            }

            if (tmp32)
                eob = i;
        }
    }
    *eob_ptr = eob + 1;
}

void svt_av1_quantize_fp_c(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,
                           const int16_t *round_ptr, const int16_t *quant_ptr,
                           const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                           TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
                           const int16_t *scan, const int16_t *iscan) {
    quantize_fp_helper_c(coeff_ptr,
                         n_coeffs,
                         zbin_ptr,
                         round_ptr,
                         quant_ptr,
                         quant_shift_ptr,
                         qcoeff_ptr,
                         dqcoeff_ptr,
                         dequant_ptr,
                         eob_ptr,
                         scan,
                         iscan,
                         NULL,
                         NULL,
                         0);
}

void svt_av1_quantize_fp_qm_c(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,
                              const int16_t *round_ptr, const int16_t *quant_ptr,
                              const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                              TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
                              const int16_t *scan, const int16_t *iscan, const QmVal *qm_ptr,
                              const QmVal *iqm_ptr, int16_t log_scale) {
    quantize_fp_helper_c(coeff_ptr,
                         n_coeffs,
                         zbin_ptr,
                         round_ptr,
                         quant_ptr,
                         quant_shift_ptr,
                         qcoeff_ptr,
                         dqcoeff_ptr,
                         dequant_ptr,
                         eob_ptr,
                         scan,
                         iscan,
                         qm_ptr,
                         iqm_ptr,
                         log_scale);
}

static void highbd_quantize_fp_helper_c(
    const TranLow *coeff_ptr, intptr_t count, const int16_t *zbin_ptr, const int16_t *round_ptr,
    const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
    TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan,
    const int16_t *iscan, const QmVal *qm_ptr, const QmVal *iqm_ptr, int16_t log_scale) {
    int       i;
    int       eob   = -1;
    const int shift = 16 - log_scale;
    (void)zbin_ptr;
    (void)quant_shift_ptr;
    (void)iscan;

    if (qm_ptr || iqm_ptr) {
        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < count; i++) {
            const int   rc      = scan[i];
            const int   coeff   = coeff_ptr[rc];
            const QmVal wt      = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const QmVal iwt     = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
            const int   dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;
            const int     coeff_sign = coeff < 0 ? -1 : 0;
            const int64_t abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;
            if (abs_coeff * wt >= (dequant_ptr[rc != 0] << (AOM_QM_BITS - (1 + log_scale)))) {
                const int64_t tmp = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale);
                const int     abs_qcoeff  = (int)((tmp * quant_ptr[rc != 0] * wt) >>
                                             (shift + AOM_QM_BITS));
                qcoeff_ptr[rc]            = (TranLow)((abs_qcoeff ^ coeff_sign) - coeff_sign);
                const TranLow abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
                dqcoeff_ptr[rc]           = (TranLow)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
                if (abs_qcoeff)
                    eob = i;
            } else {
                qcoeff_ptr[rc]  = 0;
                dqcoeff_ptr[rc] = 0;
            }
        }
    } else {
        const int log_scaled_round_arr[2] = {
            ROUND_POWER_OF_TWO(round_ptr[0], log_scale),
            ROUND_POWER_OF_TWO(round_ptr[1], log_scale),
        };
        for (i = 0; i < count; i++) {
            const int rc               = scan[i];
            const int coeff            = coeff_ptr[rc];
            const int rc01             = (rc != 0);
            const int coeff_sign       = coeff < 0 ? -1 : 0;
            const int abs_coeff        = (coeff ^ coeff_sign) - coeff_sign;
            const int log_scaled_round = log_scaled_round_arr[rc01];
            if ((abs_coeff << (1 + log_scale)) >= dequant_ptr[rc01]) {
                const int     quant       = quant_ptr[rc01];
                const int     dequant     = dequant_ptr[rc01];
                const int64_t tmp         = (int64_t)abs_coeff + log_scaled_round;
                const int     abs_qcoeff  = (int)((tmp * quant) >> shift);
                qcoeff_ptr[rc]            = (TranLow)((abs_qcoeff ^ coeff_sign) - coeff_sign);
                const TranLow abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
                if (abs_qcoeff)
                    eob = i;
                dqcoeff_ptr[rc] = (TranLow)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
            } else {
                qcoeff_ptr[rc]  = 0;
                dqcoeff_ptr[rc] = 0;
            }
        }
    }
    *eob_ptr = eob + 1;
}

void svt_av1_highbd_quantize_fp_c(const TranLow *coeff_ptr, intptr_t count, const int16_t *zbin_ptr,
                                  const int16_t *round_ptr, const int16_t *quant_ptr,
                                  const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                                  TranLow *dqcoeff_ptr, const int16_t *dequant_ptr,
                                  uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan,
                                  int16_t log_scale) {
    highbd_quantize_fp_helper_c(coeff_ptr,
                                count,
                                zbin_ptr,
                                round_ptr,
                                quant_ptr,
                                quant_shift_ptr,
                                qcoeff_ptr,
                                dqcoeff_ptr,
                                dequant_ptr,
                                eob_ptr,
                                scan,
                                iscan,
                                NULL,
                                NULL,
                                log_scale);
}

void svt_av1_quantize_fp_32x32_c(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                 const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                 TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                 const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan,
                                 const int16_t *iscan) {
    quantize_fp_helper_c(coeff_ptr,
                         n_coeffs,
                         zbin_ptr,
                         round_ptr,
                         quant_ptr,
                         quant_shift_ptr,
                         qcoeff_ptr,
                         dqcoeff_ptr,
                         dequant_ptr,
                         eob_ptr,
                         scan,
                         iscan,
                         NULL,
                         NULL,
                         1);
}

void svt_av1_quantize_fp_64x64_c(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                 const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                 TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                 const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan,
                                 const int16_t *iscan) {
    quantize_fp_helper_c(coeff_ptr,
                         n_coeffs,
                         zbin_ptr,
                         round_ptr,
                         quant_ptr,
                         quant_shift_ptr,
                         qcoeff_ptr,
                         dqcoeff_ptr,
                         dequant_ptr,
                         eob_ptr,
                         scan,
                         iscan,
                         NULL,
                         NULL,
                         2);
}

void svt_av1_quantize_fp_facade(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                const MacroblockPlane *p, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                uint16_t *eob_ptr, const ScanOrder *sc, const QuantParam *qparam) {
    const QmVal *qm_ptr  = qparam->qmatrix;
    const QmVal *iqm_ptr = qparam->iqmatrix;

    if (qm_ptr || iqm_ptr)
        svt_av1_quantize_fp_qm(coeff_ptr,
                               n_coeffs,
                               p->zbin_qtx,
                               p->round_fp_qtx,
                               p->quant_fp_qtx,
                               p->quant_shift_qtx,
                               qcoeff_ptr,
                               dqcoeff_ptr,
                               p->dequant_qtx,
                               eob_ptr,
                               sc->scan,
                               sc->iscan,
                               qm_ptr,
                               iqm_ptr,
                               qparam->log_scale);
    else {
        switch (qparam->log_scale) {
        case 0:
            svt_av1_quantize_fp(coeff_ptr,
                                n_coeffs,
                                p->zbin_qtx,
                                p->round_fp_qtx,
                                p->quant_fp_qtx,
                                p->quant_shift_qtx,
                                qcoeff_ptr,
                                dqcoeff_ptr,
                                p->dequant_qtx,
                                eob_ptr,
                                sc->scan,
                                sc->iscan);
            break;
        case 1:
            svt_av1_quantize_fp_32x32(coeff_ptr,
                                      n_coeffs,
                                      p->zbin_qtx,
                                      p->round_fp_qtx,
                                      p->quant_fp_qtx,
                                      p->quant_shift_qtx,
                                      qcoeff_ptr,
                                      dqcoeff_ptr,
                                      p->dequant_qtx,
                                      eob_ptr,
                                      sc->scan,
                                      sc->iscan);
            break;
        case 2:
            svt_av1_quantize_fp_64x64(coeff_ptr,
                                      n_coeffs,
                                      p->zbin_qtx,
                                      p->round_fp_qtx,
                                      p->quant_fp_qtx,
                                      p->quant_shift_qtx,
                                      qcoeff_ptr,
                                      dqcoeff_ptr,
                                      p->dequant_qtx,
                                      eob_ptr,
                                      sc->scan,
                                      sc->iscan);
            break;
        default: assert(0);
        }
    }
}

void svt_av1_highbd_quantize_fp_facade(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                       const MacroblockPlane *p, TranLow *qcoeff_ptr,
                                       TranLow *dqcoeff_ptr, uint16_t *eob_ptr, const ScanOrder *sc,
                                       const QuantParam *qparam) {
    const QmVal *qm_ptr  = qparam->qmatrix;
    const QmVal *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        svt_av1_highbd_quantize_fp_qm(coeff_ptr,
                                      n_coeffs,
                                      p->zbin_qtx,
                                      p->round_fp_qtx,
                                      p->quant_fp_qtx,
                                      p->quant_shift_qtx,
                                      qcoeff_ptr,
                                      dqcoeff_ptr,
                                      p->dequant_qtx,
                                      eob_ptr,
                                      sc->scan,
                                      sc->iscan,
                                      qm_ptr,
                                      iqm_ptr,
                                      qparam->log_scale);
    } else {
        svt_av1_highbd_quantize_fp(coeff_ptr,
                                   n_coeffs,
                                   p->zbin_qtx,
                                   p->round_fp_qtx,
                                   p->quant_fp_qtx,
                                   p->quant_shift_qtx,
                                   qcoeff_ptr,
                                   dqcoeff_ptr,
                                   p->dequant_qtx,
                                   eob_ptr,
                                   sc->scan,
                                   sc->iscan,
                                   qparam->log_scale);
    }
}

void svt_av1_highbd_quantize_fp_qm_c(const TranLow *coeff_ptr, intptr_t count,
                                     const int16_t *zbin_ptr, const int16_t *round_ptr,
                                     const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                     TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                     const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                     const int16_t *scan, const int16_t *iscan, const QmVal *qm_ptr,
                                     const QmVal *iqm_ptr, int16_t log_scale) {
    highbd_quantize_fp_helper_c(coeff_ptr,
                                count,
                                zbin_ptr,
                                round_ptr,
                                quant_ptr,
                                quant_shift_ptr,
                                qcoeff_ptr,
                                dqcoeff_ptr,
                                dequant_ptr,
                                eob_ptr,
                                scan,
                                iscan,
                                qm_ptr,
                                iqm_ptr,
                                log_scale);
}

// Hsan: code clean up; from static to extern as now used @ more than 1 file

static const int8_t eob_to_pos_small[33] = {
    0, 1, 2, // 0-2
    3, 3, // 3-4
    4, 4, 4, 4, // 5-8
    5, 5, 5, 5, 5, 5, 5, 5, // 9-16
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 // 17-32
};

static const int8_t eob_to_pos_large[17] = {
    6, // place holder
    7, // 33-64
    8,
    8, // 65-128
    9,
    9,
    9,
    9, // 129-256
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    10, // 257-512
    11 // 513-
};

static INLINE int32_t get_eob_pos_token(const int32_t eob, int32_t *const extra) {
    int32_t t;

    if (eob < 33)
        t = eob_to_pos_small[eob];
    else {
        const int32_t e = AOMMIN((eob - 1) >> 5, 16);
        t               = eob_to_pos_large[e];
    }

    *extra = eob - eb_k_eob_group_start[t];

    return t;
}
// Transform end of block bit estimation
static int get_eob_cost(int eob, const LvMapEobCost *txb_eob_costs, const LvMapCoeffCost *txb_costs,
                        TxClass tx_class) {
    int       eob_extra;
    const int eob_pt        = get_eob_pos_token(eob, &eob_extra);
    int       eob_cost      = 0;
    const int eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;
    eob_cost                = txb_eob_costs->eob_cost[eob_multi_ctx][eob_pt - 1];

    if (eb_k_eob_offset_bits[eob_pt] > 0) {
        const int eob_ctx   = eob_pt - 3;
        const int eob_shift = eb_k_eob_offset_bits[eob_pt] - 1;
        const int bit       = (eob_extra & (1 << eob_shift)) ? 1 : 0;
        eob_cost += txb_costs->eob_extra_cost[eob_ctx][bit];
        const int offset_bits = eb_k_eob_offset_bits[eob_pt];
        if (offset_bits > 1)
            eob_cost += av1_cost_literal(offset_bits - 1);
    }
    return eob_cost;
}
static INLINE int get_lower_levels_ctx_general(int is_last, int scan_idx, int bwl, int height,
                                               const uint8_t *levels, int coeff_idx, TxSize tx_size,
                                               TxClass tx_class) {
    if (is_last) {
        if (scan_idx == 0)
            return 0;
        if (scan_idx <= (height << bwl) >> 3)
            return 1;
        if (scan_idx <= (height << bwl) >> 2)
            return 2;
        return 3;
    }
    return get_lower_levels_ctx(levels, coeff_idx, bwl, tx_size, tx_class);
}

static INLINE int32_t get_golomb_cost(int32_t abs_qc) {
    if (abs_qc >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
        const int32_t r      = abs_qc - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
        const int32_t length = get_msb(r) + 1;
        return av1_cost_literal(2 * length - 1);
    }
    return 0;
}
static INLINE int get_br_cost(TranLow level, const int *coeff_lps) {
    const int base_range = AOMMIN(level - 1 - NUM_BASE_LEVELS, COEFF_BASE_RANGE);
    return coeff_lps[base_range] + get_golomb_cost(level);
}
static INLINE int get_coeff_cost_general(int is_last, int ci, TranLow abs_qc, int sign,
                                         int coeff_ctx, int dc_sign_ctx,
                                         const LvMapCoeffCost *txb_costs, int bwl, TxClass tx_class,
                                         const uint8_t *levels) {
    int cost = 0;
    if (is_last)
        cost += txb_costs->base_eob_cost[coeff_ctx][AOMMIN(abs_qc, 3) - 1];
    else
        cost += txb_costs->base_cost[coeff_ctx][AOMMIN(abs_qc, 3)];
    if (abs_qc != 0) {
        if (ci == 0)
            cost += txb_costs->dc_sign_cost[dc_sign_ctx][sign];
        else
            cost += av1_cost_literal(1);
        if (abs_qc > NUM_BASE_LEVELS) {
            int br_ctx;
            if (is_last)
                br_ctx = get_br_ctx_eob(ci, bwl, tx_class);
            else
                br_ctx = get_br_ctx(levels, ci, bwl, tx_class);
            cost += get_br_cost(abs_qc, txb_costs->lps_cost[br_ctx]);
        }
    }
    return cost;
}
static INLINE int64_t get_coeff_dist(TranLow tcoeff, TranLow dqcoeff, int shift) {
    const int64_t diff = ((int64_t)tcoeff - dqcoeff) * (int64_t)(1lu << shift);
    return diff * diff;
}
static INLINE void get_qc_dqc_low(TranLow abs_qc, int sign, int dqv, int shift, TranLow *qc_low,
                                  TranLow *dqc_low) {
    TranLow abs_qc_low = abs_qc - 1;
    *qc_low            = (-sign ^ abs_qc_low) + sign;
    assert((sign ? -abs_qc_low : abs_qc_low) == *qc_low);
    TranLow abs_dqc_low = (abs_qc_low * dqv) >> shift;
    *dqc_low            = (-sign ^ abs_dqc_low) + sign;
    assert((sign ? -abs_dqc_low : abs_dqc_low) == *dqc_low);
}
static const int golomb_bits_cost[32] = {
    0,       512,     512 * 3, 512 * 3, 512 * 5, 512 * 5, 512 * 5, 512 * 5,
    512 * 7, 512 * 7, 512 * 7, 512 * 7, 512 * 7, 512 * 7, 512 * 7, 512 * 7,
    512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9,
    512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9, 512 * 9};
static const int  golomb_cost_diff[32] = {0, 512, 512 * 2, 0, 512 * 2, 0,       0, 0, 512 * 2, 0, 0,
                                          0, 0,   0,       0, 0,       512 * 2, 0, 0, 0,       0, 0,
                                          0, 0,   0,       0, 0,       0,       0, 0, 0,       0};
static INLINE int get_br_cost_with_diff(TranLow level, const int *coeff_lps, int *diff) {
    const int base_range  = AOMMIN(level - 1 - NUM_BASE_LEVELS, COEFF_BASE_RANGE);
    int       golomb_bits = 0;
    if (level <= COEFF_BASE_RANGE + 1 + NUM_BASE_LEVELS)
        *diff += coeff_lps[base_range + COEFF_BASE_RANGE + 1];

    if (level >= COEFF_BASE_RANGE + 1 + NUM_BASE_LEVELS) {
        int r = level - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
        if (r < 32) {
            golomb_bits = golomb_bits_cost[r];
            *diff += golomb_cost_diff[r];
        } else {
            golomb_bits = get_golomb_cost(level);
            *diff += (r & (r - 1)) == 0 ? 1024 : 0;
        }
    }

    return coeff_lps[base_range] + golomb_bits;
}
static AOM_FORCE_INLINE int get_two_coeff_cost_simple(int ci, TranLow abs_qc, int coeff_ctx,
                                                      const LvMapCoeffCost *txb_costs, int bwl,
                                                      TxClass tx_class, const uint8_t *levels,
                                                      int *cost_low) {
    // this simple version assumes the coeff's scan_idx is not DC (scan_idx != 0)
    // and not the last (scan_idx != eob - 1)
    assert(ci > 0);
    //assert(abs_qc + 4 < 4);
    int cost = txb_costs->base_cost[coeff_ctx][AOMMIN(abs_qc, 3)];
    int diff = 0;
    if (abs_qc <= 3)
        diff = txb_costs->base_cost[coeff_ctx][abs_qc + 4];
    if (abs_qc) {
        cost += av1_cost_literal(1);
        if (abs_qc > NUM_BASE_LEVELS) {
            const int br_ctx      = get_br_ctx(levels, ci, bwl, tx_class);
            int       brcost_diff = 0;
            cost += get_br_cost_with_diff(abs_qc, txb_costs->lps_cost[br_ctx], &brcost_diff);
            diff += brcost_diff;
        }
    }
    *cost_low = cost - diff;

    return cost;
}
static INLINE int get_coeff_cost_eob(int ci, TranLow abs_qc, int sign, int coeff_ctx,
                                     int dc_sign_ctx, const LvMapCoeffCost *txb_costs, int bwl,
                                     TxClass tx_class) {
    int cost = 0;
    cost += txb_costs->base_eob_cost[coeff_ctx][AOMMIN(abs_qc, 3) - 1];
    if (abs_qc != 0) {
        if (ci == 0)
            cost += txb_costs->dc_sign_cost[dc_sign_ctx][sign];
        else
            cost += av1_cost_literal(1);
        if (abs_qc > NUM_BASE_LEVELS) {
            int br_ctx;
            br_ctx = get_br_ctx_eob(ci, bwl, tx_class);
            cost += get_br_cost(abs_qc, txb_costs->lps_cost[br_ctx]);
        }
    }
    return cost;
}

static INLINE int get_dqv(const int16_t *dequant, int coeff_idx, const QmVal *iqm_ptr) {
    int dqv = dequant[!!coeff_idx];
    if (iqm_ptr != NULL)
        dqv = ((iqm_ptr[coeff_idx] * dqv) + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
    return dqv;
}

static AOM_FORCE_INLINE void update_coeff_eob(
    int *accu_rate, int64_t *accu_dist, uint16_t *eob, int *nz_num, int *nz_ci, int si,
    TxSize tx_size, TxClass tx_class, int bwl, int height, int dc_sign_ctx, int64_t rdmult,
    int shift, const int16_t *dequant, const int16_t *scan, const LvMapEobCost *txb_eob_costs,
    const LvMapCoeffCost *txb_costs, const TranLow *tcoeff, TranLow *qcoeff, TranLow *dqcoeff,
    uint8_t *levels, int sharpness, const QmVal *iqm_ptr) {
    assert(si != *eob - 1);
    const int     ci        = scan[si];
    const int     dqv       = get_dqv(dequant, ci, iqm_ptr);
    const TranLow qc        = qcoeff[ci];
    const int     coeff_ctx = get_lower_levels_ctx(levels, ci, bwl, tx_size, tx_class);
    if (qc == 0)
        *accu_rate += txb_costs->base_cost[coeff_ctx][0];
    else {
        int           lower_level = 0;
        const TranLow abs_qc      = abs(qc);
        const TranLow tqc         = tcoeff[ci];
        const TranLow dqc         = dqcoeff[ci];
        const int     sign        = (qc < 0) ? 1 : 0;
        const int64_t dist0       = get_coeff_dist(tqc, 0, shift);
        int64_t       dist        = get_coeff_dist(tqc, dqc, shift) - dist0;
        int           rate        = get_coeff_cost_general(
            0, ci, abs_qc, sign, coeff_ctx, dc_sign_ctx, txb_costs, bwl, tx_class, levels);
        int64_t rd = RDCOST(rdmult, *accu_rate + rate, *accu_dist + dist);

        TranLow qc_low, dqc_low;
        TranLow abs_qc_low;
        int64_t dist_low, rd_low;
        int     rate_low;
        if (abs_qc == 1) {
            abs_qc_low = 0;
            dqc_low = qc_low = 0;
            dist_low         = 0;
            rate_low         = txb_costs->base_cost[coeff_ctx][0];
            rd_low           = RDCOST(rdmult, *accu_rate + rate_low, *accu_dist);
        } else {
            get_qc_dqc_low(abs_qc, sign, dqv, shift, &qc_low, &dqc_low);
            abs_qc_low = abs_qc - 1;
            dist_low   = get_coeff_dist(tqc, dqc_low, shift) - dist0;
            rate_low   = get_coeff_cost_general(
                0, ci, abs_qc_low, sign, coeff_ctx, dc_sign_ctx, txb_costs, bwl, tx_class, levels);
            rd_low = RDCOST(rdmult, *accu_rate + rate_low, *accu_dist + dist_low);
        }

        int       lower_level_new_eob = 0;
        const int new_eob             = si + 1;
        const int coeff_ctx_new_eob   = get_lower_levels_ctx_eob(bwl, height, si);
        const int new_eob_cost        = get_eob_cost(new_eob, txb_eob_costs, txb_costs, tx_class);
        int       rate_coeff_eob =
            new_eob_cost +
            get_coeff_cost_eob(
                ci, abs_qc, sign, coeff_ctx_new_eob, dc_sign_ctx, txb_costs, bwl, tx_class);
        int64_t dist_new_eob = dist;
        int64_t rd_new_eob   = RDCOST(rdmult, rate_coeff_eob, dist_new_eob);

        if (abs_qc_low > 0) {
            const int rate_coeff_eob_low =
                new_eob_cost +
                get_coeff_cost_eob(
                    ci, abs_qc_low, sign, coeff_ctx_new_eob, dc_sign_ctx, txb_costs, bwl, tx_class);
            const int64_t dist_new_eob_low = dist_low;
            const int64_t rd_new_eob_low   = RDCOST(rdmult, rate_coeff_eob_low, dist_new_eob_low);
            if (rd_new_eob_low < rd_new_eob) {
                lower_level_new_eob = 1;
                rd_new_eob          = rd_new_eob_low;
                rate_coeff_eob      = rate_coeff_eob_low;
                dist_new_eob        = dist_new_eob_low;
            }
        }

        if (rd_low < rd) {
            lower_level = 1;
            rd          = rd_low;
            rate        = rate_low;
            dist        = dist_low;
        }

        if (sharpness == 0 && rd_new_eob < rd) {
            for (int ni = 0; ni < *nz_num; ++ni) {
                int last_ci                          = nz_ci[ni];
                levels[get_padded_idx(last_ci, bwl)] = 0;
                qcoeff[last_ci]                      = 0;
                dqcoeff[last_ci]                     = 0;
            }
            *eob        = new_eob;
            *nz_num     = 0;
            *accu_rate  = rate_coeff_eob;
            *accu_dist  = dist_new_eob;
            lower_level = lower_level_new_eob;
        } else {
            *accu_rate += rate;
            *accu_dist += dist;
        }

        if (lower_level) {
            qcoeff[ci]                      = qc_low;
            dqcoeff[ci]                     = dqc_low;
            levels[get_padded_idx(ci, bwl)] = AOMMIN(abs_qc_low, INT8_MAX);
        }
        if (qcoeff[ci]) {
            nz_ci[*nz_num] = ci;
            ++*nz_num;
        }
    }
}
static INLINE void update_coeff_general(int *accu_rate, int64_t *accu_dist, int si, int eob,
                                        TxSize tx_size, TxClass tx_class, int bwl, int height,
                                        int64_t rdmult, int shift, int dc_sign_ctx,
                                        const int16_t *dequant, const int16_t *scan,
                                        const LvMapCoeffCost *txb_costs, const TranLow *tcoeff,
                                        TranLow *qcoeff, TranLow *dqcoeff, uint8_t *levels,
                                        const QmVal *iqm_ptr) {
    const int     ci        = scan[si];
    const int     dqv       = get_dqv(dequant, ci, iqm_ptr);
    const TranLow qc        = qcoeff[ci];
    const int     is_last   = si == (eob - 1);
    const int     coeff_ctx = get_lower_levels_ctx_general(
        is_last, si, bwl, height, levels, ci, tx_size, tx_class);
    if (qc == 0)
        *accu_rate += txb_costs->base_cost[coeff_ctx][0];
    else {
        const int     sign   = (qc < 0) ? 1 : 0;
        const TranLow abs_qc = abs(qc);
        const TranLow tqc    = tcoeff[ci];
        const TranLow dqc    = dqcoeff[ci];
        const int64_t dist   = get_coeff_dist(tqc, dqc, shift);
        const int64_t dist0  = get_coeff_dist(tqc, 0, shift);
        const int     rate   = get_coeff_cost_general(
            is_last, ci, abs_qc, sign, coeff_ctx, dc_sign_ctx, txb_costs, bwl, tx_class, levels);
        const int64_t rd = RDCOST(rdmult, rate, dist);

        TranLow qc_low, dqc_low;
        TranLow abs_qc_low;
        int64_t dist_low, rd_low;
        int     rate_low;
        if (abs_qc == 1) {
            abs_qc_low = qc_low = dqc_low = 0;
            dist_low                      = dist0;
            rate_low                      = txb_costs->base_cost[coeff_ctx][0];
        } else {
            get_qc_dqc_low(abs_qc, sign, dqv, shift, &qc_low, &dqc_low);
            abs_qc_low = abs_qc - 1;
            dist_low   = get_coeff_dist(tqc, dqc_low, shift);
            rate_low   = get_coeff_cost_general(is_last,
                                              ci,
                                              abs_qc_low,
                                              sign,
                                              coeff_ctx,
                                              dc_sign_ctx,
                                              txb_costs,
                                              bwl,
                                              tx_class,
                                              levels);
        }

        rd_low = RDCOST(rdmult, rate_low, dist_low);
        if (rd_low < rd) {
            qcoeff[ci]                      = qc_low;
            dqcoeff[ci]                     = dqc_low;
            levels[get_padded_idx(ci, bwl)] = AOMMIN(abs_qc_low, INT8_MAX);
            *accu_rate += rate_low;
            *accu_dist += dist_low - dist0;
        } else {
            *accu_rate += rate;
            *accu_dist += dist - dist0;
        }
    }
}

static AOM_FORCE_INLINE void update_coeff_simple(
    int *accu_rate, int si, int eob, TxSize tx_size, TxClass tx_class, int bwl, int64_t rdmult,
    int shift, const int16_t *dequant, const int16_t *scan, const LvMapCoeffCost *txb_costs,
    const TranLow *tcoeff, TranLow *qcoeff, TranLow *dqcoeff, uint8_t *levels,
    const QmVal *iqm_ptr) {
    const int dqv = get_dqv(dequant, scan[si], iqm_ptr);
    (void)eob;
    // this simple version assumes the coeff's scan_idx is not DC (scan_idx != 0)
    // and not the last (scan_idx != eob - 1)
    assert(si != eob - 1);
    assert(si > 0);
    const int     ci        = scan[si];
    const TranLow qc        = qcoeff[ci];
    const int     coeff_ctx = get_lower_levels_ctx(levels, ci, bwl, tx_size, tx_class);
    if (qc == 0)
        *accu_rate += txb_costs->base_cost[coeff_ctx][0];
    else {
        const TranLow abs_qc   = abs(qc);
        const TranLow abs_tqc  = abs(tcoeff[ci]);
        const TranLow abs_dqc  = abs(dqcoeff[ci]);
        int           rate_low = 0;
        const int     rate     = get_two_coeff_cost_simple(
            ci, abs_qc, coeff_ctx, txb_costs, bwl, tx_class, levels, &rate_low);
        if (abs_dqc < abs_tqc) {
            *accu_rate += rate;
            return;
        }

        const int64_t dist = get_coeff_dist(abs_tqc, abs_dqc, shift);
        const int64_t rd   = RDCOST(rdmult, rate, dist);

        const TranLow abs_qc_low  = abs_qc - 1;
        const TranLow abs_dqc_low = (abs_qc_low * dqv) >> shift;
        const int64_t dist_low    = get_coeff_dist(abs_tqc, abs_dqc_low, shift);
        const int64_t rd_low      = RDCOST(rdmult, rate_low, dist_low);

        if (rd_low < rd) {
            const int sign                  = (qc < 0) ? 1 : 0;
            qcoeff[ci]                      = (-sign ^ abs_qc_low) + sign;
            dqcoeff[ci]                     = (-sign ^ abs_dqc_low) + sign;
            levels[get_padded_idx(ci, bwl)] = AOMMIN(abs_qc_low, INT8_MAX);
            *accu_rate += rate_low;
        } else
            *accu_rate += rate;
    }
}
static INLINE void update_skip(int *accu_rate, int64_t accu_dist, uint16_t *eob, int nz_num,
                               int *nz_ci, int64_t rdmult, int skip_cost, int non_skip_cost,
                               TranLow *qcoeff, TranLow *dqcoeff, int sharpness) {
    const int64_t rd         = RDCOST(rdmult, *accu_rate + non_skip_cost, accu_dist);
    const int64_t rd_new_eob = RDCOST(rdmult, skip_cost, 0);
    if (sharpness == 0 && rd_new_eob < rd) {
        for (int i = 0; i < nz_num; ++i) {
            const int ci = nz_ci[i];
            qcoeff[ci]   = 0;
            dqcoeff[ci]  = 0;
            // no need to set up levels because this is the last step
            // levels[get_padded_idx(ci, bwl)] = 0;
        }
        *accu_rate = 0;
        *eob       = 0;
    }
}
enum {
    NO_AQ             = 0,
    VARIANCE_AQ       = 1,
    COMPLEXITY_AQ     = 2,
    CYCLIC_REFRESH_AQ = 3,
    AQ_MODE_COUNT // This should always be the last member of the enum
} UENUM1BYTE(AQ_MODE);
enum {
    NO_DELTA_Q   = 0,
    DELTA_Q_ONLY = 1,
    DELTA_Q_LF   = 2,
    DELTAQ_MODE_COUNT // This should always be the last member of the enum
} UENUM1BYTE(DELTAQ_MODE);

// These numbers are empirically obtained.
static const int plane_rd_mult[REF_TYPES][PLANE_TYPES] = {
#if TUNE_CHROMA_SSIM
    {17, 13},
    {16, 10},
#else
    {17, 20},
    {16, 20},
#endif
};

/*
 * Reduce the number of non-zero quantized coefficients before getting to the main/complex RDOQ stage
 * (it performs an early check of whether to zero out each of the non-zero quantized coefficients,
 * and updates the quantized coeffs if it is determined it can be zeroed out).
 */
static INLINE void update_coeff_eob_fast(uint16_t *eob, int shift, const int16_t *dequant_ptr,
                                         const int16_t *scan, const TranLow *coeff_ptr,
                                         TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr) {
    int eob_out = *eob;
    int zbin[2] = {dequant_ptr[0] + ROUND_POWER_OF_TWO(dequant_ptr[0] * 70, 7),
                   dequant_ptr[1] + ROUND_POWER_OF_TWO(dequant_ptr[1] * 70, 7)};
    for (int i = *eob - 1; i >= 0; i--) {
        const int rc         = scan[i];
        const int qcoeff     = qcoeff_ptr[rc];
        const int coeff      = coeff_ptr[rc];
        const int coeff_sign = -(coeff < 0);
        int64_t   abs_coeff  = (coeff ^ coeff_sign) - coeff_sign;
        if (((abs_coeff << (1 + shift)) < zbin[rc != 0]) || (qcoeff == 0)) {
            eob_out--;
            qcoeff_ptr[rc]  = 0;
            dqcoeff_ptr[rc] = 0;
        } else
            break;
    }
    *eob = eob_out;
}

// look-up table for sqrt of number of pixels in a transform block
// rounded up to the nearest integer.
static const int sqrt_tx_pixels_2d[TX_SIZES_ALL] = {
    4, 8, 16, 32, 32, 6, 6, 12, 12, 23, 23, 32, 32, 8, 8, 16, 16, 23, 23};

void svt_fast_optimize_b(const TranLow *coeff_ptr, const MacroblockPlane *p, TranLow *qcoeff_ptr,
                         TranLow *dqcoeff_ptr, uint16_t *eob, TxSize tx_size, TxType tx_type)

{
    const ScanOrder *const scan_order = &av1_scan_orders[tx_size][tx_type];
    const int16_t         *scan       = scan_order->scan;
    const int              shift      = av1_get_tx_scale_tab[tx_size];
    update_coeff_eob_fast(eob, shift, p->dequant_qtx, scan, coeff_ptr, qcoeff_ptr, dqcoeff_ptr);
}
void svt_av1_optimize_b(ModeDecisionContext *md_context, int16_t txb_skip_context,
                        int16_t dc_sign_context, const TranLow *coeff_ptr, int32_t stride,
                        intptr_t n_coeffs, const MacroblockPlane *p, TranLow *qcoeff_ptr,
                        TranLow *dqcoeff_ptr, uint16_t *eob, const ScanOrder *sc,
                        const QuantParam *qparam, TxSize tx_size, TxType tx_type, Bool is_inter,
                        uint8_t use_sharpness, uint8_t delta_q_present, uint8_t picture_qp,
                        uint32_t lambda, int plane)

{
    (void)stride;
    (void)n_coeffs;
    (void)sc;
    int sharpness = 0; // No Sharpness
    int fast_mode = (md_context->rdoq_ctrls.eob_fast_y_inter && is_inter && !plane) ||
            (md_context->rdoq_ctrls.eob_fast_y_intra && !is_inter && !plane) ||
            (md_context->rdoq_ctrls.eob_fast_uv_inter && is_inter && plane) ||
            (md_context->rdoq_ctrls.eob_fast_uv_intra && !is_inter && plane)
        ? 1
        : 0;
    const ScanOrder *const scan_order = &av1_scan_orders[tx_size][tx_type];
    const int16_t         *scan       = scan_order->scan;
    const int              shift      = av1_get_tx_scale_tab[tx_size];
    const PlaneType        plane_type = plane;
    const TxSize           txs_ctx    = get_txsize_entropy_ctx_tab[tx_size];
    const TxClass          tx_class   = tx_type_to_class[tx_type];
    const int              bwl        = get_txb_bwl_tab[tx_size];
    const int              width      = get_txb_wide_tab[tx_size];
    const int              height     = get_txb_high_tab[tx_size];
    assert(width == (1 << bwl));
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *txb_costs =
        &md_context->md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][plane_type];
    const int           eob_multi_size = txsize_log2_minus4[tx_size];
    const LvMapEobCost *txb_eob_costs =
        &md_context->md_rate_estimation_ptr->eob_frac_bits[eob_multi_size][plane_type];
    const int non_skip_cost = txb_costs->txb_skip_cost[txb_skip_context][0];
    const int skip_cost     = txb_costs->txb_skip_cost[txb_skip_context][1];
    const int eob_cost      = get_eob_cost(*eob, txb_eob_costs, txb_costs, tx_class);

    int sq_size_idx = 7 - (int)svt_log2f(md_context->blk_geom->sq_size);
    if (eob_cost < (int)(width * height * sq_size_idx * md_context->rdoq_ctrls.early_exit_th)) {
        if (skip_cost < non_skip_cost) {
            return;
        }
    }

    if (fast_mode) {
        update_coeff_eob_fast(eob, shift, p->dequant_qtx, scan, coeff_ptr, qcoeff_ptr, dqcoeff_ptr);
        if (*eob == 0)
            return;
    }
    int       rweight = 100;
    const int rshift  = 2;
    if (use_sharpness && delta_q_present && plane == 0) {
        int diff = md_context->sb_ptr->qindex - quantizer_to_qindex[picture_qp];
        if (diff < 0) {
            sharpness = 1;
            rweight   = 0;
        }
    }
    const int64_t rdmult =
        (((((int64_t)lambda * plane_rd_mult[is_inter][plane_type]) * rweight) / 100) + 2) >> rshift;
    uint8_t        levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);

    if (*eob > 1)
        svt_av1_txb_init_levels(qcoeff_ptr, width, height, levels);
    int accu_rate = eob_cost;

    int64_t       accu_dist  = 0;
    int           si         = *eob - 1;
    const int     ci         = scan[si];
    const TranLow qc         = qcoeff_ptr[ci];
    const TranLow abs_qc     = abs(qc);
    const int     sign       = qc < 0;
    const int     max_nz_num = 2;
    int           nz_num     = 1;
    int           nz_ci[3]   = {ci, 0, 0};

    if (abs_qc >= 2) {
        update_coeff_general(&accu_rate,
                             &accu_dist,
                             si,
                             *eob,
                             tx_size,
                             tx_class,
                             bwl,
                             height,
                             rdmult,
                             shift,
                             dc_sign_context,
                             p->dequant_qtx,
                             scan,
                             txb_costs,
                             coeff_ptr,
                             qcoeff_ptr,
                             dqcoeff_ptr,
                             levels,
                             qparam->iqmatrix);
        --si;
    } else {
        assert(abs_qc == 1);
        const int coeff_ctx = get_lower_levels_ctx_eob(bwl, height, si);
        accu_rate += get_coeff_cost_eob(
            ci, abs_qc, sign, coeff_ctx, dc_sign_context, txb_costs, bwl, tx_class);

        const TranLow tqc   = coeff_ptr[ci];
        const TranLow dqc   = dqcoeff_ptr[ci];
        const int64_t dist  = get_coeff_dist(tqc, dqc, shift);
        const int64_t dist0 = get_coeff_dist(tqc, 0, shift);
        accu_dist += dist - dist0;
        --si;
    }

#define UPDATE_COEFF_EOB_CASE(tx_class_literal)                       \
    case tx_class_literal:                                            \
        for (; si >= 0 && nz_num <= max_nz_num && !fast_mode; --si) { \
            update_coeff_eob(&accu_rate,                              \
                             &accu_dist,                              \
                             eob,                                     \
                             &nz_num,                                 \
                             nz_ci,                                   \
                             si,                                      \
                             tx_size,                                 \
                             tx_class_literal,                        \
                             bwl,                                     \
                             height,                                  \
                             dc_sign_context,                         \
                             rdmult,                                  \
                             shift,                                   \
                             p->dequant_qtx,                          \
                             scan,                                    \
                             txb_eob_costs,                           \
                             txb_costs,                               \
                             coeff_ptr,                               \
                             qcoeff_ptr,                              \
                             dqcoeff_ptr,                             \
                             levels,                                  \
                             sharpness,                               \
                             qparam->iqmatrix);                       \
        }                                                             \
        break;
    switch (tx_class) {
        UPDATE_COEFF_EOB_CASE(TX_CLASS_2D);
        UPDATE_COEFF_EOB_CASE(TX_CLASS_HORIZ);
        UPDATE_COEFF_EOB_CASE(TX_CLASS_VERT);
#undef UPDATE_COEFF_EOB_CASE
    default: assert(false);
    }

    if (si == -1 && nz_num <= max_nz_num) {
        update_skip(&accu_rate,
                    accu_dist,
                    eob,
                    nz_num,
                    nz_ci,
                    rdmult,
                    skip_cost,
                    non_skip_cost,
                    qcoeff_ptr,
                    dqcoeff_ptr,
                    sharpness);
    }

#define UPDATE_COEFF_SIMPLE_CASE(tx_class_literal) \
    case tx_class_literal:                         \
        for (; si >= 1; --si) {                    \
            update_coeff_simple(&accu_rate,        \
                                si,                \
                                *eob,              \
                                tx_size,           \
                                tx_class_literal,  \
                                bwl,               \
                                rdmult,            \
                                shift,             \
                                p->dequant_qtx,    \
                                scan,              \
                                txb_costs,         \
                                coeff_ptr,         \
                                qcoeff_ptr,        \
                                dqcoeff_ptr,       \
                                levels,            \
                                qparam->iqmatrix); \
        }                                          \
        break;
    switch (tx_class) {
        UPDATE_COEFF_SIMPLE_CASE(TX_CLASS_2D);
        UPDATE_COEFF_SIMPLE_CASE(TX_CLASS_HORIZ);
        UPDATE_COEFF_SIMPLE_CASE(TX_CLASS_VERT);
#undef UPDATE_COEFF_SIMPLE_CASE
    default: assert(false);
    }

    // DC position
    if (si == 0) {
        // no need to update accu_dist because it's not used after this point
        int64_t dummy_dist = 0;
        update_coeff_general(&accu_rate,
                             &dummy_dist,
                             si,
                             *eob,
                             tx_size,
                             tx_class,
                             bwl,
                             height,
                             rdmult,
                             shift,
                             dc_sign_context,
                             p->dequant_qtx,
                             scan,
                             txb_costs,
                             coeff_ptr,
                             qcoeff_ptr,
                             dqcoeff_ptr,
                             levels,
                             qparam->iqmatrix);
    }
}

static INLINE void set_dc_sign(int32_t *cul_level, int32_t dc_val) {
    if (dc_val < 0)
        *cul_level |= 1 << COEFF_CONTEXT_BITS;
    else if (dc_val > 0)
        *cul_level += 2 << COEFF_CONTEXT_BITS;
}
static INLINE TxSize aom_av1_get_adjusted_tx_size(TxSize tx_size) {
    switch (tx_size) {
    case TX_64X64:
    case TX_64X32:
    case TX_32X64: return TX_32X32;
    case TX_64X16: return TX_32X16;
    case TX_16X64: return TX_16X32;
    default: return tx_size;
    }
}
int32_t av1_quantize_inv_quantize_light(PictureControlSet *pcs_ptr, int32_t *coeff,
                                        int32_t *quant_coeff, int32_t *recon_coeff, uint32_t qindex,
                                        TxSize txsize, uint16_t *eob,
                                        uint32_t *count_non_zero_coeffs, uint32_t bit_depth,
                                        TxType tx_type) {
    SequenceControlSet    *scs_ptr          = pcs_ptr->scs_ptr;
    uint32_t               q_index          = qindex;
    const ScanOrder *const scan_order       = &av1_scan_orders[txsize][tx_type];
    const int32_t          n_coeffs         = av1_get_max_eob(txsize);
    int32_t                qmatrix_level    = (IS_2D_TRANSFORM(tx_type) &&
                             pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.using_qmatrix)
                          ? pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.qm[AOM_PLANE_Y]
                          : NUM_QM_LEVELS - 1;
    TxSize                 adjusted_tx_size = aom_av1_get_adjusted_tx_size(txsize);
    const QmVal           *q_matrix =
        pcs_ptr->parent_pcs_ptr->gqmatrix[qmatrix_level][AOM_PLANE_Y][adjusted_tx_size];
    const QmVal *iq_matrix =
        pcs_ptr->parent_pcs_ptr->giqmatrix[qmatrix_level][AOM_PLANE_Y][adjusted_tx_size];

    if (q_matrix == NULL && iq_matrix == NULL) {
        if (bit_depth > EB_EIGHT_BIT) {
            svt_aom_highbd_quantize_b((TranLow *)coeff,
                                      n_coeffs,
                                      scs_ptr->quants_bd.y_zbin[q_index],
                                      scs_ptr->quants_bd.y_round[q_index],
                                      scs_ptr->quants_bd.y_quant[q_index],
                                      scs_ptr->quants_bd.y_quant_shift[q_index],
                                      quant_coeff,
                                      (TranLow *)recon_coeff,
                                      scs_ptr->deq_bd.v_dequant_qtx[q_index],
                                      eob,
                                      scan_order->scan,
                                      scan_order->iscan,
                                      q_matrix,
                                      iq_matrix,
                                      av1_get_tx_scale_tab[txsize]);
        } else {
            //  svt_aom_quantize_b_c_ii((TranLow *)coeff,
            svt_aom_quantize_b((TranLow *)coeff,
                               n_coeffs,
                               scs_ptr->quants_8bit.v_zbin[q_index],
                               scs_ptr->quants_8bit.v_round[q_index],
                               scs_ptr->quants_8bit.v_quant[q_index],
                               scs_ptr->quants_8bit.v_quant_shift[q_index],
                               quant_coeff,
                               (TranLow *)recon_coeff,
                               scs_ptr->deq_8bit.y_dequant_qtx[q_index],
                               eob,
                               scan_order->scan,
                               scan_order->iscan,
                               q_matrix,
                               iq_matrix,
                               av1_get_tx_scale_tab[txsize]);
        }
    } else {
        if (bit_depth > EB_EIGHT_BIT) {
            svt_av1_highbd_quantize_b_qm((TranLow *)coeff,
                                         n_coeffs,
                                         scs_ptr->quants_bd.y_zbin[q_index],
                                         scs_ptr->quants_bd.y_round[q_index],
                                         scs_ptr->quants_bd.y_quant[q_index],
                                         scs_ptr->quants_bd.y_quant_shift[q_index],
                                         quant_coeff,
                                         (TranLow *)recon_coeff,
                                         scs_ptr->deq_bd.v_dequant_qtx[q_index],
                                         eob,
                                         scan_order->scan,
                                         scan_order->iscan,
                                         q_matrix,
                                         iq_matrix,
                                         av1_get_tx_scale_tab[txsize]);
        } else {
            svt_av1_quantize_b_qm((TranLow *)coeff,
                                  n_coeffs,
                                  scs_ptr->quants_8bit.v_zbin[q_index],
                                  scs_ptr->quants_8bit.v_round[q_index],
                                  scs_ptr->quants_8bit.v_quant[q_index],
                                  scs_ptr->quants_8bit.v_quant_shift[q_index],
                                  quant_coeff,
                                  (TranLow *)recon_coeff,
                                  scs_ptr->deq_8bit.y_dequant_qtx[q_index],
                                  eob,
                                  scan_order->scan,
                                  scan_order->iscan,
                                  q_matrix,
                                  iq_matrix,
                                  av1_get_tx_scale_tab[txsize]);
        }
    }

    *count_non_zero_coeffs = *eob;
    return 0;
}

int32_t svt_av1_compute_cul_level_c(const int16_t *const scan, const int32_t *const quant_coeff,
                                    uint16_t *eob) {
    int32_t cul_level = 0;
    for (int32_t c = 0; c < *eob; ++c) {
        const int16_t pos   = scan[c];
        const int32_t v     = quant_coeff[pos];
        int32_t       level = ABS(v);
        cul_level += level;
        // Early exit the loop if cul_level reaches COEFF_CONTEXT_MASK
        if (cul_level >= COEFF_CONTEXT_MASK)
            break;
    }

    cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);
    // DC value
    set_dc_sign(&cul_level, quant_coeff[0]);
    return cul_level;
}

int32_t av1_quantize_inv_quantize(PictureControlSet *pcs_ptr, ModeDecisionContext *md_context,
                                  int32_t *coeff, const uint32_t coeff_stride, int32_t *quant_coeff,
                                  int32_t *recon_coeff, uint32_t qindex,
                                  int32_t segmentation_qp_offset, uint32_t width, uint32_t height,
                                  TxSize txsize, uint16_t *eob, uint32_t *count_non_zero_coeffs,

                                  uint32_t component_type, uint32_t bit_depth, TxType tx_type,
                                  ModeDecisionCandidateBuffer *candidate_buffer,
                                  int16_t txb_skip_context, int16_t dc_sign_context,
                                  PredictionMode pred_mode, Bool is_intra_bc, uint32_t lambda,
                                  Bool is_encode_pass) {
    (void)candidate_buffer;
    (void)is_encode_pass;
    (void)coeff_stride;
    (void)is_intra_bc;
    SequenceControlSet *scs_ptr          = pcs_ptr->scs_ptr;
    int32_t             plane            = component_type == COMPONENT_LUMA ? AOM_PLANE_Y
                               : COMPONENT_CHROMA_CB                        ? AOM_PLANE_U
                                                                            : AOM_PLANE_V;
    int32_t             qmatrix_level    = (IS_2D_TRANSFORM(tx_type) &&
                             pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.using_qmatrix)
                       ? pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.qm[plane]
                       : NUM_QM_LEVELS - 1;
    TxSize              adjusted_tx_size = aom_av1_get_adjusted_tx_size(txsize);
    MacroblockPlane     candidate_plane;
    const QmVal        *q_matrix =
        pcs_ptr->parent_pcs_ptr->gqmatrix[qmatrix_level][plane][adjusted_tx_size];
    const QmVal *iq_matrix =
        pcs_ptr->parent_pcs_ptr->giqmatrix[qmatrix_level][plane][adjusted_tx_size];
    uint32_t q_index = pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
        ? qindex
        : (uint32_t)CLIP3(0,
                          255,
                          (int32_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx +
                              segmentation_qp_offset);
    if (component_type != COMPONENT_LUMA) {
        const int8_t offset = (component_type == COMPONENT_CHROMA_CB)
            ? pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params
                  .delta_q_dc[1] // we are assuming delta_q_ac == delta_q_dc
            : pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.delta_q_dc[2];
        q_index += offset;
        q_index = (uint32_t)CLIP3(0, 255, (int32_t)q_index);
    }
    if (bit_depth == EB_EIGHT_BIT) {
        if (component_type == COMPONENT_LUMA) {
            candidate_plane.quant_qtx       = scs_ptr->quants_8bit.y_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_8bit.y_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_8bit.y_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_8bit.y_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_8bit.y_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_8bit.y_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_8bit.y_dequant_qtx[q_index];
        }

        else if (component_type == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_qtx       = scs_ptr->quants_8bit.u_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_8bit.u_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_8bit.u_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_8bit.u_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_8bit.u_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_8bit.u_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_8bit.u_dequant_qtx[q_index];
        }

        else {
            candidate_plane.quant_qtx       = scs_ptr->quants_8bit.v_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_8bit.v_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_8bit.v_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_8bit.v_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_8bit.v_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_8bit.v_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_8bit.v_dequant_qtx[q_index];
        }
    } else {
        if (component_type == COMPONENT_LUMA) {
            candidate_plane.quant_qtx       = scs_ptr->quants_bd.y_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_bd.y_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_bd.y_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_bd.y_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_bd.y_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_bd.y_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_bd.y_dequant_qtx[q_index];
        }

        else if (component_type == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_qtx       = scs_ptr->quants_bd.u_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_bd.u_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_bd.u_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_bd.u_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_bd.u_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_bd.u_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_bd.u_dequant_qtx[q_index];
        }

        else {
            candidate_plane.quant_qtx       = scs_ptr->quants_bd.v_quant[q_index];
            candidate_plane.quant_fp_qtx    = scs_ptr->quants_bd.v_quant_fp[q_index];
            candidate_plane.round_fp_qtx    = scs_ptr->quants_bd.v_round_fp[q_index];
            candidate_plane.quant_shift_qtx = scs_ptr->quants_bd.v_quant_shift[q_index];
            candidate_plane.zbin_qtx        = scs_ptr->quants_bd.v_zbin[q_index];
            candidate_plane.round_qtx       = scs_ptr->quants_bd.v_round[q_index];
            candidate_plane.dequant_qtx     = scs_ptr->deq_bd.v_dequant_qtx[q_index];
        }
    }

    const ScanOrder *const scan_order =
        &av1_scan_orders[txsize][tx_type]; //get_scan(tx_size, tx_type);

    const int32_t n_coeffs = av1_get_max_eob(txsize);

    QuantParam qparam;

    qparam.log_scale = av1_get_tx_scale_tab[txsize];
    qparam.tx_size   = txsize;
    qparam.qmatrix   = q_matrix;
    qparam.iqmatrix  = iq_matrix;

    Bool is_inter = (pred_mode >= NEARESTMV);
    Bool perform_rdoq;

    // If rdoq_level is specified in the command line instruction, set perform_rdoq accordingly.
    perform_rdoq            = ((md_context->md_staging_skip_rdoq == FALSE || is_encode_pass) &&
                    md_context->rdoq_level);
    const int dequant_shift = md_context->hbd_mode_decision
        ? pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth - 5
        : 3;
    const int qstep         = candidate_plane.dequant_qtx[1] /*[AC]*/ >> dequant_shift;
    if (!is_encode_pass) {
        if ((md_context->rdoq_ctrls.dct_dct_only && tx_type != DCT_DCT) ||
            (md_context->rdoq_ctrls.skip_uv && component_type != COMPONENT_LUMA))
            perform_rdoq = 0;
    }
    if (perform_rdoq && md_context->rdoq_ctrls.satd_factor != ((uint8_t)~0)) {
        int       satd  = svt_aom_satd(coeff, n_coeffs);
        const int shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[txsize]);

        satd = RIGHT_SIGNED_SHIFT(satd, shift);
        satd >>= (pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth - 8);
        const int skip_block_trellis = ((uint64_t)satd >
                                        (uint64_t)md_context->rdoq_ctrls.satd_factor * qstep *
                                            sqrt_tx_pixels_2d[txsize]);
        if (skip_block_trellis)
            perform_rdoq = 0;
    }

    if (perform_rdoq &&
        ((!component_type && md_context->rdoq_ctrls.fp_q_y) ||
         (component_type && md_context->rdoq_ctrls.fp_q_uv))) {
        if ((bit_depth > EB_EIGHT_BIT) || (is_encode_pass && scs_ptr->is_16bit_pipeline)) {
            svt_av1_highbd_quantize_fp_facade((TranLow *)coeff,
                                              n_coeffs,
                                              &candidate_plane,
                                              quant_coeff,
                                              (TranLow *)recon_coeff,
                                              eob,
                                              scan_order,
                                              &qparam);
        } else {
            svt_av1_quantize_fp_facade((TranLow *)coeff,
                                       n_coeffs,
                                       &candidate_plane,
                                       quant_coeff,
                                       (TranLow *)recon_coeff,
                                       eob,
                                       scan_order,
                                       &qparam);
        }
    } else {
        if ((bit_depth > EB_EIGHT_BIT) || (is_encode_pass && scs_ptr->is_16bit_pipeline)) {
            svt_av1_highbd_quantize_b_facade((TranLow *)coeff,
                                             n_coeffs,
                                             &candidate_plane,
                                             quant_coeff,
                                             (TranLow *)recon_coeff,
                                             eob,
                                             scan_order,
                                             &qparam);
        } else {
            av1_quantize_b_facade_ii((TranLow *)coeff,
                                     n_coeffs,
                                     &candidate_plane,
                                     quant_coeff,
                                     (TranLow *)recon_coeff,
                                     eob,
                                     scan_order,
                                     &qparam);
        }
    }
    if (perform_rdoq && *eob != 0) {
        int eob_perc = (*eob) * 100 / (width * height);
        if (eob_perc >= md_context->rdoq_ctrls.eob_th) {
            perform_rdoq = 0;
        }
        if (perform_rdoq && (eob_perc >= md_context->rdoq_ctrls.eob_fast_th)) {
            svt_fast_optimize_b((TranLow *)coeff,
                                &candidate_plane,
                                quant_coeff,
                                (TranLow *)recon_coeff,
                                eob,
                                txsize,
                                tx_type);
        }
        if (perform_rdoq == 0) {
            if ((bit_depth > EB_EIGHT_BIT) || (is_encode_pass && scs_ptr->is_16bit_pipeline)) {
                svt_av1_highbd_quantize_b_facade((TranLow *)coeff,
                                                 n_coeffs,
                                                 &candidate_plane,
                                                 quant_coeff,
                                                 (TranLow *)recon_coeff,
                                                 eob,
                                                 scan_order,
                                                 &qparam);
            } else {
                av1_quantize_b_facade_ii((TranLow *)coeff,
                                         n_coeffs,
                                         &candidate_plane,
                                         quant_coeff,
                                         (TranLow *)recon_coeff,
                                         eob,
                                         scan_order,
                                         &qparam);
            }
        }
    }
    if (perform_rdoq && *eob != 0) {
        // Perform rdoq
        svt_av1_optimize_b(md_context,
                           txb_skip_context,
                           dc_sign_context,
                           (TranLow *)coeff,
                           coeff_stride,
                           n_coeffs,
                           &candidate_plane,
                           quant_coeff,
                           (TranLow *)recon_coeff,
                           eob,
                           scan_order,
                           &qparam,
                           txsize,
                           tx_type,
                           is_inter,
                           scs_ptr->vq_ctrls.sharpness_ctrls.rdoq,
                           pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present,
                           pcs_ptr->picture_qp,
                           lambda,
                           (component_type == COMPONENT_LUMA) ? 0 : 1);
    }

    *count_non_zero_coeffs = *eob;
    if (!md_context->rate_est_ctrls.update_skip_ctx_dc_sign_ctx)
        return 0;

    // Derive cul_level
    return svt_av1_compute_cul_level(scan_order->scan, quant_coeff, eob);
}

void inv_transform_recon_wrapper(uint8_t *pred_buffer, uint32_t pred_offset, uint32_t pred_stride,
                                 uint8_t *rec_buffer, uint32_t rec_offset, uint32_t rec_stride,
                                 int32_t *rec_coeff_buffer, uint32_t coeff_offset, Bool hbd,
                                 TxSize txsize, TxType transform_type, PlaneType component_type,
                                 uint32_t eob) {
    if (hbd) {
        av1_inv_transform_recon(rec_coeff_buffer + coeff_offset,
                                CONVERT_TO_BYTEPTR(((uint16_t *)pred_buffer) + pred_offset),
                                pred_stride,
                                CONVERT_TO_BYTEPTR(((uint16_t *)rec_buffer) + rec_offset),
                                rec_stride,
                                txsize,
                                EB_TEN_BIT,
                                transform_type,
                                component_type,
                                eob,
                                0 /*lossless*/);
    } else {
        av1_inv_transform_recon8bit(rec_coeff_buffer + coeff_offset,
                                    pred_buffer + pred_offset,
                                    pred_stride,
                                    rec_buffer + rec_offset,
                                    rec_stride,
                                    txsize,
                                    transform_type,
                                    component_type,
                                    eob,
                                    0 /*lossless*/);
    }
}
/*
  tx path for light PD1 chroma
*/
void full_loop_chroma_light_pd1(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                ModeDecisionCandidateBuffer *candidate_buffer,
                                EbPictureBufferDesc         *input_picture_ptr,
                                uint32_t input_cb_origin_in_index, uint32_t blk_chroma_origin_index,
                                COMPONENT_TYPE component_type, uint32_t chroma_qindex,
                                uint64_t  cb_full_distortion[DIST_CALC_TOTAL],
                                uint64_t  cr_full_distortion[DIST_CALC_TOTAL],
                                uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits) {
    uint32_t     full_lambda = context_ptr->hbd_mode_decision
            ? context_ptr->full_lambda_md[EB_10_BIT_MD]
            : context_ptr->full_lambda_md[EB_8_BIT_MD];
    uint32_t     nz_count_dummy;
    const TxSize tx_size_uv = context_ptr->blk_geom->txsize_uv[0][0];

    EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;
    // If Cb component not detected as complex, can use TX shortcuts
    if (context_ptr->use_tx_shortcuts_mds3 &&
        (context_ptr->chroma_complexity == COMPONENT_LUMA ||
         context_ptr->chroma_complexity == COMPONENT_CHROMA_CR)) {
        pf_shape = N4_SHAPE;
    } else {
        uint8_t use_pfn4_cond = 0;
        if (context_ptr->lpd1_tx_ctrls.use_uv_shortcuts_on_y_coeffs &&
            (context_ptr->chroma_complexity == COMPONENT_LUMA ||
             context_ptr->chroma_complexity == COMPONENT_CHROMA_CR)) {
            const uint16_t th = ((context_ptr->blk_geom->tx_width_uv[0][0] >> 4) *
                                 (context_ptr->blk_geom->tx_height_uv[0][0] >> 4));
            use_pfn4_cond     = (candidate_buffer->count_non_zero_coeffs < th) ||
                    !candidate_buffer->block_has_coeff
                    ? 1
                    : 0;
        }
        if (use_pfn4_cond)
            pf_shape = N4_SHAPE;
    }
    assert(tx_size_uv < TX_SIZES_ALL);
    const int32_t chroma_shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size_uv]) * 2;
    uint32_t      bwidth       = context_ptr->blk_geom->tx_width_uv[0][0];
    uint32_t      bheight      = context_ptr->blk_geom->tx_height_uv[0][0];
    if (pf_shape) {
        bwidth  = MAX((bwidth >> pf_shape), 4);
        bheight = (bheight >> pf_shape);
    }
    if (component_type == COMPONENT_CHROMA || component_type == COMPONENT_CHROMA_CB) {
        residual_kernel(input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cb,
                        candidate_buffer->prediction_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cb,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cb,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);

        // Cb Transform
        av1_estimate_transform(
            &(((int16_t *)candidate_buffer->residual_ptr->buffer_cb)[blk_chroma_origin_index]),
            candidate_buffer->residual_ptr->stride_cb,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[0]),
            NOT_USED_VALUE,
            tx_size_uv,
            &context_ptr->three_quad_energy,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            candidate_buffer->candidate_ptr->transform_type_uv,
            PLANE_TYPE_UV,
            pf_shape);
        candidate_buffer->quantized_dc[1][0] = av1_quantize_inv_quantize(
            pcs_ptr,
            context_ptr,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[0]),
            NOT_USED_VALUE,
            &(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_cb)[0]),
            &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb)[0]),
            chroma_qindex,
            0,
            context_ptr->blk_geom->tx_width_uv[0][0],
            context_ptr->blk_geom->tx_height_uv[0][0],
            tx_size_uv,
            &candidate_buffer->eob[1][0],
            &nz_count_dummy,
            COMPONENT_CHROMA_CB,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            candidate_buffer->candidate_ptr->transform_type_uv,
            candidate_buffer,
            0,
            0,
            candidate_buffer->candidate_ptr->pred_mode,
            0,
            full_lambda,
            FALSE);

        picture_full_distortion32_bits_single(
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[0]),
            &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb)[0]),
            context_ptr->blk_geom->tx_width_uv[0][0],
            bwidth,
            bheight,
            cb_full_distortion,
            candidate_buffer->eob[1][0]);
        cb_full_distortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
            cb_full_distortion[DIST_CALC_RESIDUAL], chroma_shift);
        cb_full_distortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
            cb_full_distortion[DIST_CALC_PREDICTION], chroma_shift);
        candidate_buffer->u_has_coeff = (candidate_buffer->eob[1][0] > 0);
    }

    pf_shape = context_ptr->pf_ctrls.pf_shape;
    // If Cr component not detected as complex, can use TX shortcuts
    if (context_ptr->use_tx_shortcuts_mds3 &&
        (context_ptr->chroma_complexity == COMPONENT_LUMA ||
         context_ptr->chroma_complexity == COMPONENT_CHROMA_CB)) {
        pf_shape = N4_SHAPE;
    } else {
        uint8_t use_pfn4_cond = 0;
        if (context_ptr->lpd1_tx_ctrls.use_uv_shortcuts_on_y_coeffs &&
            (context_ptr->chroma_complexity == COMPONENT_LUMA ||
             context_ptr->chroma_complexity == COMPONENT_CHROMA_CB)) {
            const uint16_t th = ((context_ptr->blk_geom->tx_width_uv[0][0] >> 4) *
                                 (context_ptr->blk_geom->tx_height_uv[0][0] >> 4));
            use_pfn4_cond     = (candidate_buffer->count_non_zero_coeffs < th) ||
                    !candidate_buffer->block_has_coeff
                    ? 1
                    : 0;
        }
        if (use_pfn4_cond)
            pf_shape = N4_SHAPE;
    }
    bwidth  = context_ptr->blk_geom->tx_width_uv[0][0];
    bheight = context_ptr->blk_geom->tx_height_uv[0][0];
    if (pf_shape) {
        bwidth  = MAX((bwidth >> pf_shape), 4);
        bheight = (bheight >> pf_shape);
    }

    if (component_type == COMPONENT_CHROMA || component_type == COMPONENT_CHROMA_CR) {
        //Cr Residual
        residual_kernel(input_picture_ptr->buffer_cr,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cr,
                        candidate_buffer->prediction_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cr,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cr,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);
        // Cr Transform
        av1_estimate_transform(
            &(((int16_t *)candidate_buffer->residual_ptr->buffer_cr)[blk_chroma_origin_index]),
            candidate_buffer->residual_ptr->stride_cr,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[0]),
            NOT_USED_VALUE,
            tx_size_uv,
            &context_ptr->three_quad_energy,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            candidate_buffer->candidate_ptr->transform_type_uv,
            PLANE_TYPE_UV,
            pf_shape);
        candidate_buffer->quantized_dc[2][0] = av1_quantize_inv_quantize(
            pcs_ptr,
            context_ptr,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[0]),
            NOT_USED_VALUE,
            &(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_cr)[0]),
            &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr)[0]),
            chroma_qindex,
            0,
            context_ptr->blk_geom->tx_width_uv[0][0],
            context_ptr->blk_geom->tx_height_uv[0][0],
            tx_size_uv,
            &candidate_buffer->eob[2][0],
            &nz_count_dummy,
            COMPONENT_CHROMA_CR,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            candidate_buffer->candidate_ptr->transform_type_uv,
            candidate_buffer,
            0,
            0,
            candidate_buffer->candidate_ptr->pred_mode,
            0,
            full_lambda,
            FALSE);

        picture_full_distortion32_bits_single(
            &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[0]),
            &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr)[0]),
            context_ptr->blk_geom->tx_width_uv[0][0],
            bwidth,
            bheight,
            cr_full_distortion,
            candidate_buffer->eob[2][0]);

        cr_full_distortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
            cr_full_distortion[DIST_CALC_RESIDUAL], chroma_shift);
        cr_full_distortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
            cr_full_distortion[DIST_CALC_PREDICTION], chroma_shift);
        candidate_buffer->v_has_coeff = (candidate_buffer->eob[2][0] > 0);
    }

    //CHROMA-ONLY
    av1_txb_estimate_coeff_bits(context_ptr,
                                0,
                                NULL,
                                pcs_ptr,
                                candidate_buffer,
                                NOT_USED_VALUE,
                                0,
                                candidate_buffer->quant_coeff_ptr,
                                NOT_USED_VALUE,
                                candidate_buffer->eob[1][0],
                                candidate_buffer->eob[2][0],
                                NOT_USED_VALUE,
                                cb_coeff_bits,
                                cr_coeff_bits,
                                NOT_USED_VALUE,
                                tx_size_uv,
                                NOT_USED_VALUE,
                                candidate_buffer->candidate_ptr->transform_type_uv,
                                component_type);
}
/****************************************
 ************  Full loop ****************
****************************************/
void svt_aom_full_loop_uv(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          ModeDecisionCandidateBuffer *candidate_buffer,
                          EbPictureBufferDesc *input_picture_ptr, COMPONENT_TYPE component_type,
                          uint32_t chroma_qindex,
                          uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU],
                          uint64_t cb_full_distortion[DIST_CALC_TOTAL],
                          uint64_t cr_full_distortion[DIST_CALC_TOTAL], uint64_t *cb_coeff_bits,
                          uint64_t *cr_coeff_bits, Bool is_full_loop) {
    EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
        ? svt_full_distortion_kernel16_bits
        : svt_spatial_full_distortion_kernel;
    uint32_t             *cb_count_non_zero_coeffs   = count_non_zero_coeffs[1];
    uint32_t             *cr_count_non_zero_coeffs   = count_non_zero_coeffs[2];
    EB_ALIGN(16) uint64_t txb_full_distortion[3][DIST_CALC_TOTAL];
    candidate_buffer->u_has_coeff = 0;
    candidate_buffer->v_has_coeff = 0;
    int16_t *chroma_residual_ptr;
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];

    context_ptr->three_quad_energy = 0;

    const uint8_t tx_depth      = candidate_buffer->candidate_ptr->tx_depth;
    const Bool    is_inter      = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                           candidate_buffer->candidate_ptr->use_intrabc)
                ? TRUE
                : FALSE;
    const int     tu_count      = tx_depth
                 ? 1
                 : context_ptr->blk_geom
              ->txb_count[candidate_buffer->candidate_ptr->tx_depth]; //NM: 128x128 exeption
    uint32_t      txb_1d_offset = 0;

    int txb_itr = 0;
    do {
        const uint32_t txb_origin_x = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
        const uint32_t txb_origin_y = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
        int32_t        cropped_tx_width_uv = MIN(
            context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
            pcs_ptr->parent_pcs_ptr->aligned_width / 2 -
                ((context_ptr->sb_origin_x + ((txb_origin_x >> 3) << 3)) >> 1));
        int32_t cropped_tx_height_uv = MIN(
            context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr],
            pcs_ptr->parent_pcs_ptr->aligned_height / 2 -
                ((context_ptr->sb_origin_y + ((txb_origin_y >> 3) << 3)) >> 1));
        // NADER - TU
        uint32_t tu_cb_origin_index = (((txb_origin_x >> 3) << 3) +
                                       (((txb_origin_y >> 3) << 3) *
                                        candidate_buffer->residual_ptr->stride_cb)) >>
            1;
        uint32_t tu_cr_origin_index = (((txb_origin_x >> 3) << 3) +
                                       (((txb_origin_y >> 3) << 3) *
                                        candidate_buffer->residual_ptr->stride_cr)) >>
            1;
        EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;
        if (context_ptr->md_stage == MD_STAGE_3 && context_ptr->use_tx_shortcuts_mds3 &&
            context_ptr->chroma_complexity == COMPONENT_LUMA) {
            pf_shape = N4_SHAPE;
        }
        // for chroma path, use luma coeff info to make shortcut decisions (available even if MDS1 is skipped)
        else if (context_ptr->tx_shortcut_ctrls.apply_pf_on_coeffs &&
                 context_ptr->md_stage == MD_STAGE_3 &&
                 context_ptr->chroma_complexity == COMPONENT_LUMA) {
            uint8_t use_pfn4_cond = 0;

            const uint16_t th = ((context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr] >> 4) *
                                 (context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr] >> 4));
            use_pfn4_cond     = (candidate_buffer->count_non_zero_coeffs < th) ||
                    !candidate_buffer->block_has_coeff
                    ? 1
                    : 0;

            if (use_pfn4_cond)
                pf_shape = N4_SHAPE;
        }
        //    This function replaces the previous Intra Chroma mode if the LM fast
        //    cost is better.
        //    *Note - this might require that we have inv transform in the loop
        if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA ||
            component_type == COMPONENT_ALL) {
            context_ptr->cb_txb_skip_context = 0;
            context_ptr->cb_dc_sign_context  = 0;
            if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx)
                get_txb_ctx(pcs_ptr,
                            COMPONENT_CHROMA,
                            context_ptr->cb_dc_sign_level_coeff_neighbor_array,
                            ROUND_UV(context_ptr->sb_origin_x + txb_origin_x) >> 1,
                            ROUND_UV(context_ptr->sb_origin_y + txb_origin_y) >> 1,
                            context_ptr->blk_geom->bsize_uv,
                            context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                            &context_ptr->cb_txb_skip_context,
                            &context_ptr->cb_dc_sign_context);
            // Configure the Chroma Residual Ptr

            chroma_residual_ptr = //(candidate_buffer->candidate_ptr->type  == INTRA_MODE )?
                //&(((int16_t*) candidate_buffer->intraChromaResidualPtr->buffer_cb)[txb_chroma_origin_index]):
                &(((int16_t *)candidate_buffer->residual_ptr->buffer_cb)[tu_cb_origin_index]);

            // Cb Transform
            av1_estimate_transform(chroma_residual_ptr,
                                   candidate_buffer->residual_ptr->stride_cb,
                                   &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[txb_1d_offset]),
                                   NOT_USED_VALUE,
                                   context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                                   &context_ptr->three_quad_energy,
                                   context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                                   candidate_buffer->candidate_ptr->transform_type_uv,
                                   PLANE_TYPE_UV,
                                   pf_shape);

            int32_t seg_qp =
                pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled
                ? pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params
                      .feature_data[context_ptr->blk_ptr->segment_id][SEG_LVL_ALT_Q]
                : 0;
            candidate_buffer->quantized_dc[1][0] = av1_quantize_inv_quantize(
                pcs_ptr,
                context_ptr,
                &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_cb)[txb_1d_offset]),
                &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb)[txb_1d_offset]),
                chroma_qindex,
                seg_qp,
                context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr],
                context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                &candidate_buffer->eob[1][txb_itr],
                &(cb_count_non_zero_coeffs[txb_itr]),
                COMPONENT_CHROMA_CB,
                context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                candidate_buffer->candidate_ptr->transform_type_uv,
                candidate_buffer,
                context_ptr->cb_txb_skip_context,
                context_ptr->cb_dc_sign_context,
                candidate_buffer->candidate_ptr->pred_mode,
                candidate_buffer->candidate_ptr->use_intrabc,
                full_lambda,
                FALSE);

            if (is_full_loop && context_ptr->md_staging_spatial_sse_full_loop_level) {
                uint32_t cb_has_coeff = cb_count_non_zero_coeffs[txb_itr] > 0;

                if (cb_has_coeff)
                    inv_transform_recon_wrapper(
                        candidate_buffer->prediction_ptr->buffer_cb,
                        tu_cb_origin_index,
                        candidate_buffer->prediction_ptr->stride_cb,
                        candidate_buffer->recon_ptr->buffer_cb,
                        tu_cb_origin_index,
                        candidate_buffer->recon_ptr->stride_cb,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb,
                        txb_1d_offset,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->eob[1][txb_itr]);
                else
                    svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                         0,
                                         tu_cb_origin_index,
                                         candidate_buffer->recon_ptr,
                                         0,
                                         tu_cb_origin_index,
                                         0,
                                         0,
                                         context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                                         context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr],
                                         PICTURE_BUFFER_DESC_Cb_FLAG,
                                         context_ptr->hbd_mode_decision);

                uint32_t input_chroma_txb_origin_index =
                    (((context_ptr->sb_origin_y + ((txb_origin_y >> 3) << 3)) >> 1) +
                     (input_picture_ptr->origin_y >> 1)) *
                        input_picture_ptr->stride_cb +
                    (((context_ptr->sb_origin_x + ((txb_origin_x >> 3) << 3)) >> 1) +
                     (input_picture_ptr->origin_x >> 1));

                int32_t txb_uv_origin_index = (((txb_origin_x >> 3) << 3) +
                                               (((txb_origin_y >> 3) << 3) *
                                                candidate_buffer->quant_coeff_ptr->stride_cb)) >>
                    1;

                txb_full_distortion[1][DIST_CALC_PREDICTION] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_cb,
                    input_chroma_txb_origin_index,
                    input_picture_ptr->stride_cb,
                    candidate_buffer->prediction_ptr->buffer_cb,
                    txb_uv_origin_index,
                    candidate_buffer->prediction_ptr->stride_cb,
                    cropped_tx_width_uv,
                    cropped_tx_height_uv);

                txb_full_distortion[1][DIST_CALC_RESIDUAL] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_cb,
                    input_chroma_txb_origin_index,
                    input_picture_ptr->stride_cb,
                    candidate_buffer->recon_ptr->buffer_cb,
                    txb_uv_origin_index,
                    candidate_buffer->recon_ptr->stride_cb,
                    cropped_tx_width_uv,
                    cropped_tx_height_uv);

                txb_full_distortion[1][DIST_CALC_PREDICTION] <<= 4;
                txb_full_distortion[1][DIST_CALC_RESIDUAL] <<= 4;
            } else {
                // *Full Distortion (SSE)
                // *Note - there are known issues with how this distortion metric is currently
                //    calculated.  The amount of scaling between the two arrays is not
                //    equivalent.

                uint32_t bwidth  = context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr];
                uint32_t bheight = context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr];
                if (pf_shape) {
                    bwidth  = MAX((bwidth >> pf_shape), 4);
                    bheight = (bheight >> pf_shape);
                }
                picture_full_distortion32_bits_single(
                    &(((int32_t *)context_ptr->tx_coeffs->buffer_cb)[txb_1d_offset]),
                    &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb)[txb_1d_offset]),
                    context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                    bwidth,
                    bheight,
                    txb_full_distortion[1],
                    cb_count_non_zero_coeffs[txb_itr]);

                TxSize        tx_size      = context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr];
                const int32_t chroma_shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
                txb_full_distortion[1][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion[1][DIST_CALC_RESIDUAL], chroma_shift);
                txb_full_distortion[1][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion[1][DIST_CALC_PREDICTION], chroma_shift);
            }
            candidate_buffer->u_has_coeff |= ((cb_count_non_zero_coeffs[txb_itr] != 0) << txb_itr);
            cb_full_distortion[DIST_CALC_RESIDUAL] += txb_full_distortion[1][DIST_CALC_RESIDUAL];
            cb_full_distortion[DIST_CALC_PREDICTION] +=
                txb_full_distortion[1][DIST_CALC_PREDICTION];
        }

        if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA ||
            component_type == COMPONENT_ALL) {
            context_ptr->cr_txb_skip_context = 0;
            context_ptr->cr_dc_sign_context  = 0;
            if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx)
                get_txb_ctx(pcs_ptr,
                            COMPONENT_CHROMA,
                            context_ptr->cr_dc_sign_level_coeff_neighbor_array,
                            ROUND_UV(context_ptr->sb_origin_x + txb_origin_x) >> 1,
                            ROUND_UV(context_ptr->sb_origin_y + txb_origin_y) >> 1,
                            context_ptr->blk_geom->bsize_uv,
                            context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                            &context_ptr->cr_txb_skip_context,
                            &context_ptr->cr_dc_sign_context);
            // Configure the Chroma Residual Ptr

            chroma_residual_ptr = //(candidate_buffer->candidate_ptr->type  == INTRA_MODE )?
                //&(((int16_t*) candidate_buffer->intraChromaResidualPtr->buffer_cr)[txb_chroma_origin_index]):
                &(((int16_t *)candidate_buffer->residual_ptr->buffer_cr)[tu_cr_origin_index]);

            // Cr Transform
            av1_estimate_transform(chroma_residual_ptr,
                                   candidate_buffer->residual_ptr->stride_cr,
                                   &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[txb_1d_offset]),
                                   NOT_USED_VALUE,
                                   context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                                   &context_ptr->three_quad_energy,
                                   context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                                   candidate_buffer->candidate_ptr->transform_type_uv,
                                   PLANE_TYPE_UV,
                                   pf_shape);
            int32_t seg_qp =
                pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled
                ? pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params
                      .feature_data[context_ptr->blk_ptr->segment_id][SEG_LVL_ALT_Q]
                : 0;
            candidate_buffer->quantized_dc[2][0] = av1_quantize_inv_quantize(
                pcs_ptr,
                context_ptr,
                &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_cr)[txb_1d_offset]),
                &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr)[txb_1d_offset]),
                chroma_qindex,
                seg_qp,
                context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr],
                context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                &candidate_buffer->eob[2][txb_itr],
                &(cr_count_non_zero_coeffs[txb_itr]),
                COMPONENT_CHROMA_CR,
                context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                candidate_buffer->candidate_ptr->transform_type_uv,
                candidate_buffer,
                context_ptr->cr_txb_skip_context,
                context_ptr->cr_dc_sign_context,
                candidate_buffer->candidate_ptr->pred_mode,
                candidate_buffer->candidate_ptr->use_intrabc,
                full_lambda,
                FALSE);
            if (is_full_loop && context_ptr->md_staging_spatial_sse_full_loop_level) {
                uint32_t cr_has_coeff = cr_count_non_zero_coeffs[txb_itr] > 0;

                if (cr_has_coeff)
                    inv_transform_recon_wrapper(
                        candidate_buffer->prediction_ptr->buffer_cr,
                        tu_cr_origin_index,
                        candidate_buffer->prediction_ptr->stride_cr,
                        candidate_buffer->recon_ptr->buffer_cr,
                        tu_cr_origin_index,
                        candidate_buffer->recon_ptr->stride_cr,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr,
                        txb_1d_offset,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->eob[2][txb_itr]);
                else
                    svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                         0,
                                         tu_cb_origin_index,
                                         candidate_buffer->recon_ptr,
                                         0,
                                         tu_cb_origin_index,
                                         0,
                                         0,
                                         context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                                         context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr],
                                         PICTURE_BUFFER_DESC_Cr_FLAG,
                                         context_ptr->hbd_mode_decision);
                uint32_t input_chroma_txb_origin_index =
                    (((context_ptr->sb_origin_y + ((txb_origin_y >> 3) << 3)) >> 1) +
                     (input_picture_ptr->origin_y >> 1)) *
                        input_picture_ptr->stride_cb +
                    (((context_ptr->sb_origin_x + ((txb_origin_x >> 3) << 3)) >> 1) +
                     (input_picture_ptr->origin_x >> 1));

                int32_t txb_uv_origin_index = (((txb_origin_x >> 3) << 3) +
                                               (((txb_origin_y >> 3) << 3) *
                                                candidate_buffer->quant_coeff_ptr->stride_cb)) >>
                    1;

                txb_full_distortion[2][DIST_CALC_PREDICTION] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_cr,
                    input_chroma_txb_origin_index,
                    input_picture_ptr->stride_cr,
                    candidate_buffer->prediction_ptr->buffer_cr,
                    txb_uv_origin_index,
                    candidate_buffer->prediction_ptr->stride_cr,
                    cropped_tx_width_uv,
                    cropped_tx_height_uv);

                txb_full_distortion[2][DIST_CALC_RESIDUAL] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_cr,
                    input_chroma_txb_origin_index,
                    input_picture_ptr->stride_cr,
                    candidate_buffer->recon_ptr->buffer_cr,
                    txb_uv_origin_index,
                    candidate_buffer->recon_ptr->stride_cr,
                    cropped_tx_width_uv,
                    cropped_tx_height_uv);

                txb_full_distortion[2][DIST_CALC_PREDICTION] <<= 4;
                txb_full_distortion[2][DIST_CALC_RESIDUAL] <<= 4;
            } else {
                // *Full Distortion (SSE)
                // *Note - there are known issues with how this distortion metric is currently
                //    calculated.  The amount of scaling between the two arrays is not
                //    equivalent.

                uint32_t bwidth  = context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr];
                uint32_t bheight = context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr];
                if (pf_shape) {
                    bwidth  = MAX((bwidth >> pf_shape), 4);
                    bheight = (bheight >> pf_shape);
                }
                picture_full_distortion32_bits_single(
                    &(((int32_t *)context_ptr->tx_coeffs->buffer_cr)[txb_1d_offset]),
                    &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr)[txb_1d_offset]),
                    context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr],
                    bwidth,
                    bheight,
                    txb_full_distortion[2],
                    cr_count_non_zero_coeffs[txb_itr]);

                TxSize        tx_size      = context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr];
                const int32_t chroma_shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
                txb_full_distortion[2][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion[2][DIST_CALC_RESIDUAL], chroma_shift);
                txb_full_distortion[2][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion[2][DIST_CALC_PREDICTION], chroma_shift);
            }
            candidate_buffer->v_has_coeff |= ((cr_count_non_zero_coeffs[txb_itr] != 0) << txb_itr);
            cr_full_distortion[DIST_CALC_RESIDUAL] += txb_full_distortion[2][DIST_CALC_RESIDUAL];
            cr_full_distortion[DIST_CALC_PREDICTION] +=
                txb_full_distortion[2][DIST_CALC_PREDICTION];
        }

        const uint32_t txb_origin_index = txb_origin_x +
            txb_origin_y * candidate_buffer->quant_coeff_ptr->stride_y;

        // Reset the Bit Costs
        uint64_t y_txb_coeff_bits  = 0;
        uint64_t cb_txb_coeff_bits = 0;
        uint64_t cr_txb_coeff_bits = 0;

        //CHROMA-ONLY
        av1_txb_estimate_coeff_bits(context_ptr,
                                    0, //allow_update_cdf,
                                    NULL, //FRAME_CONTEXT *ec_ctx,
                                    pcs_ptr,
                                    candidate_buffer,
                                    txb_origin_index,
                                    txb_1d_offset,
                                    candidate_buffer->quant_coeff_ptr,
                                    count_non_zero_coeffs[0][txb_itr],
                                    count_non_zero_coeffs[1][txb_itr],
                                    count_non_zero_coeffs[2][txb_itr],
                                    &y_txb_coeff_bits,
                                    &cb_txb_coeff_bits,
                                    &cr_txb_coeff_bits,
                                    context_ptr->blk_geom->txsize[tx_depth][txb_itr],
                                    context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                                    candidate_buffer->candidate_ptr->transform_type[txb_itr],
                                    candidate_buffer->candidate_ptr->transform_type_uv,
                                    component_type);

        *cb_coeff_bits += cb_txb_coeff_bits;
        *cr_coeff_bits += cr_txb_coeff_bits;
        txb_1d_offset += context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr] *
            context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr];

        ++txb_itr;
    } while (txb_itr < tu_count);
}
uint64_t d1_non_square_block_decision(ModeDecisionContext *context_ptr, uint32_t d1_block_itr) {
    //compute total cost for the whole block partition
    uint64_t tot_cost      = 0;
    uint32_t first_blk_idx = context_ptr->blk_ptr->mds_idx -
        (context_ptr->blk_geom->totns - 1); //index of first block in this partition
    uint32_t blk_it;
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_sb_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_sb_lambda_md[EB_8_BIT_MD];
    for (blk_it = 0; blk_it < context_ptr->blk_geom->totns; blk_it++) {
        tot_cost += context_ptr->md_local_blk_unit[first_blk_idx + blk_it].cost;
    }
    uint64_t split_cost = svt_aom_partition_rate_cost(
        context_ptr->sb_ptr->pcs_ptr->parent_pcs_ptr,
        context_ptr,
        &context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds],
        from_shape_to_part[context_ptr->blk_geom->shape],
        full_lambda,
        context_ptr->md_rate_estimation_ptr);

    tot_cost += split_cost;
    if ((d1_block_itr == 0) ||
        (tot_cost < context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].cost)) {
        //store best partition cost in parent square
        context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].cost = tot_cost;
        context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds].part =
            from_shape_to_part[context_ptr->blk_geom->shape];
        context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds].best_d1_blk = first_blk_idx;
    }
    return tot_cost;
}

/// compute the cost of curr depth, and the depth above
void compute_depth_costs(ModeDecisionContext *context_ptr, PictureParentControlSet *pcs_ptr,
                         uint32_t curr_depth_mds, uint32_t above_depth_mds, uint32_t step,
                         uint64_t *above_depth_cost, uint64_t *curr_depth_cost) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_sb_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_sb_lambda_md[EB_8_BIT_MD];

    uint64_t above_split_rate = 0;

    /*
    ___________
    |     |     |
    |blk0 |blk1 |
    |-----|-----|
    |blk2 |blk3 |
    |_____|_____|
    */
    // current depth blocks
    uint32_t curr_depth_blk0_mds = curr_depth_mds - 3 * step;
    uint32_t curr_depth_blk1_mds = curr_depth_mds - 2 * step;
    uint32_t curr_depth_blk2_mds = curr_depth_mds - 1 * step;
    uint32_t curr_depth_blk3_mds = curr_depth_mds;
    context_ptr->md_local_blk_unit[above_depth_mds].left_neighbor_partition =
        context_ptr->md_local_blk_unit[curr_depth_blk0_mds].left_neighbor_partition;
    context_ptr->md_local_blk_unit[above_depth_mds].above_neighbor_partition =
        context_ptr->md_local_blk_unit[curr_depth_blk0_mds].above_neighbor_partition;

    // Get split rate for current depth
    above_split_rate = svt_aom_partition_rate_cost(pcs_ptr,
                                                   context_ptr,
                                                   &context_ptr->md_blk_arr_nsq[above_depth_mds],
                                                   PARTITION_SPLIT,
                                                   full_lambda,
                                                   context_ptr->md_rate_estimation_ptr);

    // Compute current depth cost
    *curr_depth_cost = context_ptr->md_local_blk_unit[curr_depth_blk0_mds].cost +
        context_ptr->md_local_blk_unit[curr_depth_blk1_mds].cost +
        context_ptr->md_local_blk_unit[curr_depth_blk2_mds].cost +
        context_ptr->md_local_blk_unit[curr_depth_blk3_mds].cost + above_split_rate;

    // Compute above depth cost
    *above_depth_cost = context_ptr->tested_blk_flag[above_depth_mds]
        ? context_ptr->md_local_blk_unit[above_depth_mds].cost
        : MAX_MODE_COST;
}

/*
 * Compare costs between depths, then update cost/splitting info in the parent blocks
 * to reflect chosen partition.  Cost comparison only performed when the all quadrants
 * of a given depth have been evaluted.
 */
uint32_t d2_inter_depth_block_decision(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                       ModeDecisionContext *context_ptr, uint32_t blk_mds,
                                       uint32_t sb_addr) {
    uint64_t         parent_depth_cost = 0, current_depth_cost = 0;
    Bool             last_depth_flag = (context_ptr->md_blk_arr_nsq[blk_mds].split_flag == FALSE);
    uint32_t         last_blk_index = blk_mds, current_depth_idx_mds = blk_mds;
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_mds);
    if (last_depth_flag) {
        while (blk_geom->is_last_quadrant) {
            //get parent idx
            uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
            if (pcs_ptr->slice_type == I_SLICE && parent_depth_idx_mds == 0 &&
                scs_ptr->seq_header.sb_size == BLOCK_128X128)
                parent_depth_cost = MAX_MODE_COST;
            else
                compute_depth_costs(context_ptr,
                                    pcs_ptr->parent_pcs_ptr,
                                    current_depth_idx_mds,
                                    parent_depth_idx_mds,
                                    blk_geom->ns_depth_offset,
                                    &parent_depth_cost,
                                    &current_depth_cost);
            if (!pcs_ptr->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[parent_depth_idx_mds])
                parent_depth_cost = MAX_MODE_COST;
            if (context_ptr->inter_depth_bias) {
                current_depth_cost = (current_depth_cost * context_ptr->inter_depth_bias) / 1000;
            }
            if (parent_depth_cost <= current_depth_cost) {
                context_ptr->md_blk_arr_nsq[parent_depth_idx_mds].split_flag = FALSE;
                context_ptr->md_local_blk_unit[parent_depth_idx_mds].cost    = parent_depth_cost;
                last_blk_index                                               = parent_depth_idx_mds;
            } else {
                context_ptr->md_local_blk_unit[parent_depth_idx_mds].cost    = current_depth_cost;
                context_ptr->md_blk_arr_nsq[parent_depth_idx_mds].part       = PARTITION_SPLIT;
                context_ptr->md_blk_arr_nsq[parent_depth_idx_mds].split_flag = TRUE;
            }

            //setup next parent inter depth
            blk_geom              = get_blk_geom_mds(parent_depth_idx_mds);
            current_depth_idx_mds = parent_depth_idx_mds;
        }
    }

    return last_blk_index;
}
void compute_depth_costs_md_skip_light_pd0(ModeDecisionContext *context_ptr,
                                           uint32_t above_depth_mds, uint32_t step,
                                           uint64_t *above_depth_cost, uint64_t *curr_depth_cost) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_sb_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_sb_lambda_md[EB_8_BIT_MD];

    *curr_depth_cost = 0;
    // sum the previous ones
    for (int i = 1; i < context_ptr->blk_geom->quadi + 1; i++) {
        uint32_t curr_depth_cur_blk_mds = context_ptr->blk_geom->sqi_mds - i * step;
        *curr_depth_cost += context_ptr->md_local_blk_unit[curr_depth_cur_blk_mds].cost;
    }
    // Add split rate to the cost of the current depth
    // Use context index 0 for the split rate as an approximation to skip call to av1_partition_rate_cost
    *curr_depth_cost += RDCOST(
        full_lambda,
        (uint64_t)context_ptr->md_rate_estimation_ptr->partition_fac_bits[0][PARTITION_SPLIT],
        0);

    *above_depth_cost = context_ptr->tested_blk_flag[above_depth_mds]
        ? context_ptr->md_local_blk_unit[above_depth_mds].cost
        : MAX_MODE_COST;
}
void compute_depth_costs_md_skip(ModeDecisionContext *context_ptr, PictureParentControlSet *pcs_ptr,
                                 uint32_t above_depth_mds, uint32_t step,
                                 uint64_t *above_depth_cost, uint64_t *curr_depth_cost) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_sb_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_sb_lambda_md[EB_8_BIT_MD];

    uint64_t above_split_rate = 0;
    *curr_depth_cost          = 0;
    // sum the previous ones
    for (int i = 1; i < context_ptr->blk_geom->quadi + 1; i++) {
        uint32_t curr_depth_cur_blk_mds = context_ptr->blk_geom->sqi_mds - i * step;

        *curr_depth_cost += context_ptr->md_local_blk_unit[curr_depth_cur_blk_mds].cost;
    }
    /*
    ___________
    |     |     |
    |blk0 |blk1 |
    |-----|-----|
    |blk2 |blk3 |
    |_____|_____|
    */
    // current depth blocks
    uint32_t curr_depth_blk0_mds = context_ptr->blk_geom->sqi_mds -
        context_ptr->blk_geom->quadi * step;
    context_ptr->md_local_blk_unit[above_depth_mds].left_neighbor_partition =
        context_ptr->md_local_blk_unit[curr_depth_blk0_mds].left_neighbor_partition;
    context_ptr->md_local_blk_unit[above_depth_mds].above_neighbor_partition =
        context_ptr->md_local_blk_unit[curr_depth_blk0_mds].above_neighbor_partition;

    above_split_rate = svt_aom_partition_rate_cost(pcs_ptr,
                                                   context_ptr,
                                                   &context_ptr->md_blk_arr_nsq[above_depth_mds],
                                                   PARTITION_SPLIT,
                                                   full_lambda,
                                                   context_ptr->md_rate_estimation_ptr);

    *curr_depth_cost += above_split_rate;

    *above_depth_cost = context_ptr->tested_blk_flag[above_depth_mds]
        ? context_ptr->md_local_blk_unit[above_depth_mds].cost
        : MAX_MODE_COST;
}
