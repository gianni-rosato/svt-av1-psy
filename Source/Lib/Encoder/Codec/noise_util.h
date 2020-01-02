/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_NOISE_UTIL_H_
#define AOM_AOM_DSP_NOISE_UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// aom_noise_tx_t is an abstraction of a transform that is used for denoising.
// It is meant to be lightweight and does hold the transformed data (as
// the user should not be manipulating the transformed data directly).
struct aom_noise_tx_t;

// Allocates and returns a aom_noise_tx_t useful for denoising the given
// block_size. The resulting aom_noise_tx_t should be free'd with
// eb_aom_noise_tx_free.
struct aom_noise_tx_t *eb_aom_noise_tx_malloc(int32_t block_size);
void                   eb_aom_noise_tx_free(struct aom_noise_tx_t *aom_noise_tx);

// Transforms the internal data and holds it in the aom_noise_tx's internal
// buffer. For compatibility with existing SIMD implementations, "data" must
// be 32-byte aligned.
void eb_aom_noise_tx_forward(struct aom_noise_tx_t *aom_noise_tx, const float *data);

// Filters aom_noise_tx's internal data using the provided noise power spectral
// density. The PSD must be at least BlockSize * BlockSize and should be
// populated with a constant or via estimates taken from
// aom_noise_tx_add_energy.
void eb_aom_noise_tx_filter(struct aom_noise_tx_t *aom_noise_tx, const float *psd);

// Performs an inverse transform using the internal transform data.
// For compatibility with existing SIMD implementations, "data" must be 32-byte
// aligned.
void eb_aom_noise_tx_inverse(struct aom_noise_tx_t *aom_noise_tx, float *data);

// Returns a default value suitable for denosing a transform of the given
// BlockSize. The noise "factor" determines the strength of the noise to
// be removed. A value of about 2.5 can be used for moderate denoising,
// where a value of 5.0 can be used for a high level of denoising.
float eb_aom_noise_psd_get_default_value(int32_t block_size, float factor);
#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus

#endif // AOM_AOM_DSP_NOISE_UTIL_H_
