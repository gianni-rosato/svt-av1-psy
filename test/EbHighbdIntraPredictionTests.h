/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#pragma once
#ifndef EbHighbdIntraPredictionTests_h
#define EbHighbdIntraPredictionTests_h

#include "aom_dsp_rtcd.h"

#ifndef NON_AVX512_SUPPORT
typedef void(*aom_highbd_dc_top_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_dc_top_predictor_func aom_highbd_dc_top_funcptr_array_opt[7] = {
                            aom_highbd_dc_top_predictor_32x8_avx512,
                            aom_highbd_dc_top_predictor_32x16_avx512,
                            aom_highbd_dc_top_predictor_32x32_avx512,
                            aom_highbd_dc_top_predictor_32x64_avx512,
                            aom_highbd_dc_top_predictor_64x16_avx512,
                            aom_highbd_dc_top_predictor_64x32_avx512,
                            aom_highbd_dc_top_predictor_64x64_avx512 };

static aom_highbd_dc_top_predictor_func aom_highbd_dc_top_funcptr_array_base[7] = {
                            eb_aom_highbd_dc_top_predictor_32x8_avx2,
                            eb_aom_highbd_dc_top_predictor_32x16_avx2,
                            eb_aom_highbd_dc_top_predictor_32x32_avx2,
                            eb_aom_highbd_dc_top_predictor_32x64_avx2,
                            eb_aom_highbd_dc_top_predictor_64x16_avx2,
                            eb_aom_highbd_dc_top_predictor_64x32_avx2,
                            eb_aom_highbd_dc_top_predictor_64x64_avx2 };

static aom_highbd_dc_top_predictor_func aom_highbd_dc_top_funcptr_array_naive[7] = {
                            eb_aom_highbd_dc_top_predictor_32x8_c,
                            eb_aom_highbd_dc_top_predictor_32x16_c,
                            eb_aom_highbd_dc_top_predictor_32x32_c,
                            eb_aom_highbd_dc_top_predictor_32x64_c,
                            eb_aom_highbd_dc_top_predictor_64x16_c,
                            eb_aom_highbd_dc_top_predictor_64x32_c,
                            eb_aom_highbd_dc_top_predictor_64x64_c};

typedef void(*aom_highbd_dc_left_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_dc_left_predictor_func aom_highbd_dc_left_funcptr_array_opt[7] = {
                            aom_highbd_dc_left_predictor_32x8_avx512,
                            aom_highbd_dc_left_predictor_32x16_avx512,
                            aom_highbd_dc_left_predictor_32x32_avx512,
                            aom_highbd_dc_left_predictor_32x64_avx512,
                            aom_highbd_dc_left_predictor_64x16_avx512,
                            aom_highbd_dc_left_predictor_64x32_avx512,
                            aom_highbd_dc_left_predictor_64x64_avx512 };

static aom_highbd_dc_left_predictor_func aom_highbd_dc_left_funcptr_array_base[7] = {
                            eb_aom_highbd_dc_left_predictor_32x8_avx2,
                            eb_aom_highbd_dc_left_predictor_32x16_avx2,
                            eb_aom_highbd_dc_left_predictor_32x32_avx2,
                            eb_aom_highbd_dc_left_predictor_32x64_avx2,
                            eb_aom_highbd_dc_left_predictor_64x16_avx2,
                            eb_aom_highbd_dc_left_predictor_64x32_avx2,
                            eb_aom_highbd_dc_left_predictor_64x64_avx2 };

static aom_highbd_dc_left_predictor_func aom_highbd_dc_left_funcptr_array_naive[7] = {
                            eb_aom_highbd_dc_left_predictor_32x8_c,
                            eb_aom_highbd_dc_left_predictor_32x16_c,
                            eb_aom_highbd_dc_left_predictor_32x32_c,
                            eb_aom_highbd_dc_left_predictor_32x64_c,
                            eb_aom_highbd_dc_left_predictor_64x16_c,
                            eb_aom_highbd_dc_left_predictor_64x32_c,
                            eb_aom_highbd_dc_left_predictor_64x64_c };


typedef void(*aom_highbd_dc_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_dc_predictor_func aom_highbd_dc_pred_funcptr_array_opt[7] = {
                            aom_highbd_dc_predictor_32x8_avx512,
                            aom_highbd_dc_predictor_32x16_avx512,
                            aom_highbd_dc_predictor_32x32_avx512,
                            aom_highbd_dc_predictor_32x64_avx512,
                            aom_highbd_dc_predictor_64x16_avx512,
                            aom_highbd_dc_predictor_64x32_avx512,
                            aom_highbd_dc_predictor_64x64_avx512 };

static aom_highbd_dc_predictor_func aom_highbd_dc_pred_funcptr_array_base[7] = {
                            eb_aom_highbd_dc_predictor_32x8_avx2,
                            eb_aom_highbd_dc_predictor_32x16_avx2,
                            eb_aom_highbd_dc_predictor_32x32_avx2,
                            eb_aom_highbd_dc_predictor_32x64_avx2,
                            eb_aom_highbd_dc_predictor_64x16_avx2,
                            eb_aom_highbd_dc_predictor_64x32_avx2,
                            eb_aom_highbd_dc_predictor_64x64_avx2 };

static aom_highbd_dc_predictor_func aom_highbd_dc_pred_funcptr_array_naive[7] = {
                            eb_aom_highbd_dc_predictor_32x8_c,
                            eb_aom_highbd_dc_predictor_32x16_c,
                            eb_aom_highbd_dc_predictor_32x32_c,
                            eb_aom_highbd_dc_predictor_32x64_c,
                            eb_aom_highbd_dc_predictor_64x16_c,
                            eb_aom_highbd_dc_predictor_64x32_c,
                            eb_aom_highbd_dc_predictor_64x64_c };


typedef void(*aom_highbd_h_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_h_predictor_func aom_highbd_h_pred_funcptr_array_opt[7] = {
                            aom_highbd_h_predictor_32x8_avx512,
                            aom_highbd_h_predictor_32x16_avx512,
                            aom_highbd_h_predictor_32x32_avx512,
                            aom_highbd_h_predictor_32x64_avx512,
                            aom_highbd_h_predictor_64x16_avx512,
                            aom_highbd_h_predictor_64x32_avx512,
                            aom_highbd_h_predictor_64x64_avx512 };

static aom_highbd_h_predictor_func aom_highbd_h_pred_funcptr_array_base[7] = {
                            eb_aom_highbd_h_predictor_32x8_avx2,
                            eb_aom_highbd_h_predictor_32x16_sse2,
                            eb_aom_highbd_h_predictor_32x32_sse2,
                            eb_aom_highbd_h_predictor_32x64_avx2,
                            eb_aom_highbd_h_predictor_64x16_avx2,
                            eb_aom_highbd_h_predictor_64x32_avx2,
                            eb_aom_highbd_h_predictor_64x64_avx2 };

static aom_highbd_h_predictor_func aom_highbd_h_pred_funcptr_array_naive[7] = {
                            eb_aom_highbd_h_predictor_32x8_c,
                            eb_aom_highbd_h_predictor_32x16_c,
                            eb_aom_highbd_h_predictor_32x32_c,
                            eb_aom_highbd_h_predictor_32x64_c,
                            eb_aom_highbd_h_predictor_64x16_c,
                            eb_aom_highbd_h_predictor_64x32_c,
                            eb_aom_highbd_h_predictor_64x64_c };

typedef void(*aom_highbd_v_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_v_predictor_func aom_highbd_v_pred_funcptr_array_opt[7] = {
                            aom_highbd_v_predictor_32x8_avx512,
                            aom_highbd_v_predictor_32x16_avx512,
                            aom_highbd_v_predictor_32x32_avx512,
                            aom_highbd_v_predictor_32x64_avx512,
                            aom_highbd_v_predictor_64x16_avx512,
                            aom_highbd_v_predictor_64x32_avx512,
                            aom_highbd_v_predictor_64x64_avx512 };

static aom_highbd_v_predictor_func aom_highbd_v_pred_funcptr_array_base[7] = {
                            eb_aom_highbd_v_predictor_32x8_avx2,
                            eb_aom_highbd_v_predictor_32x16_avx2,
                            eb_aom_highbd_v_predictor_32x32_avx2,
                            eb_aom_highbd_v_predictor_32x64_avx2,
                            eb_aom_highbd_v_predictor_64x16_avx2,
                            eb_aom_highbd_v_predictor_64x32_avx2,
                            eb_aom_highbd_v_predictor_64x64_avx2 };

static aom_highbd_v_predictor_func aom_highbd_v_pred_funcptr_array_naive[7] = {
                            eb_aom_highbd_v_predictor_32x8_c,
                            eb_aom_highbd_v_predictor_32x16_c,
                            eb_aom_highbd_v_predictor_32x32_c,
                            eb_aom_highbd_v_predictor_32x64_c,
                            eb_aom_highbd_v_predictor_64x16_c,
                            eb_aom_highbd_v_predictor_64x32_c,
                            eb_aom_highbd_v_predictor_64x64_c };


typedef void(*aom_highbd_smooth_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_smooth_predictor_func aom_smooth_predictor_funcptr_array_opt[7] = {
                            aom_highbd_smooth_predictor_32x8_avx512,
                            aom_highbd_smooth_predictor_32x16_avx512,
                            aom_highbd_smooth_predictor_32x32_avx512,
                            aom_highbd_smooth_predictor_32x64_avx512,
                            aom_highbd_smooth_predictor_64x16_avx512,
                            aom_highbd_smooth_predictor_64x32_avx512,
                            aom_highbd_smooth_predictor_64x64_avx512 };

static aom_highbd_smooth_predictor_func aom_smooth_predictor_funcptr_array_base[7] = {
                            eb_aom_highbd_smooth_predictor_32x8_avx2,
                            eb_aom_highbd_smooth_predictor_32x16_avx2,
                            eb_aom_highbd_smooth_predictor_32x32_avx2,
                            eb_aom_highbd_smooth_predictor_32x64_avx2,
                            eb_aom_highbd_smooth_predictor_64x16_avx2,
                            eb_aom_highbd_smooth_predictor_64x32_avx2,
                            eb_aom_highbd_smooth_predictor_64x64_avx2 };

typedef void(*aom_highbd_smooth_v_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_smooth_v_predictor_func aom_highbd_smooth_v_predictor_funcptr_array_opt[7] = {
                            aom_highbd_smooth_v_predictor_32x8_avx512,
                            aom_highbd_smooth_v_predictor_32x16_avx512,
                            aom_highbd_smooth_v_predictor_32x32_avx512,
                            aom_highbd_smooth_v_predictor_32x64_avx512,
                            aom_highbd_smooth_v_predictor_64x16_avx512,
                            aom_highbd_smooth_v_predictor_64x32_avx512,
                            aom_highbd_smooth_v_predictor_64x64_avx512 };

static aom_highbd_smooth_v_predictor_func aom_highbd_smooth_v_predictor_funcptr_array_base[7] = {
                            eb_aom_highbd_smooth_v_predictor_32x8_avx2,
                            eb_aom_highbd_smooth_v_predictor_32x16_avx2,
                            eb_aom_highbd_smooth_v_predictor_32x32_avx2,
                            eb_aom_highbd_smooth_v_predictor_32x64_avx2,
                            eb_aom_highbd_smooth_v_predictor_64x16_avx2,
                            eb_aom_highbd_smooth_v_predictor_64x32_avx2,
                            eb_aom_highbd_smooth_v_predictor_64x64_avx2 };

typedef void(*aom_highbd_smooth_h_predictor_func)(uint16_t *dst, ptrdiff_t stride, const uint16_t *above, const uint16_t *left, int32_t bd);

static aom_highbd_smooth_h_predictor_func aom_highbd_smooth_h_predictor_funcptr_array_opt[7] = {
                            aom_highbd_smooth_h_predictor_32x8_avx512,
                            aom_highbd_smooth_h_predictor_32x16_avx512,
                            aom_highbd_smooth_h_predictor_32x32_avx512,
                            aom_highbd_smooth_h_predictor_32x64_avx512,
                            aom_highbd_smooth_h_predictor_64x16_avx512,
                            aom_highbd_smooth_h_predictor_64x32_avx512,
                            aom_highbd_smooth_h_predictor_64x64_avx512 };

static aom_highbd_smooth_h_predictor_func aom_highbd_smooth_h_predictor_funcptr_array_base[7] = {
                            eb_aom_highbd_smooth_h_predictor_32x8_avx2,
                            eb_aom_highbd_smooth_h_predictor_32x16_avx2,
                            eb_aom_highbd_smooth_h_predictor_32x32_avx2,
                            eb_aom_highbd_smooth_h_predictor_32x64_avx2,
                            eb_aom_highbd_smooth_h_predictor_64x16_avx2,
                            eb_aom_highbd_smooth_h_predictor_64x32_avx2,
                            eb_aom_highbd_smooth_h_predictor_64x64_avx2 };

#endif
#endif
