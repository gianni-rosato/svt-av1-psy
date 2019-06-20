/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbRestoration.h"
#include "EbUnitTestUtility.h"

#include <immintrin.h>  // AVX2
#include "synonyms.h"
#include "synonyms_avx2.h"
#include "EbPictureOperators_AVX2.h"
#include "EbRestorationPick.h"
#include "EbUtility.h"

static const int32_t sizes[2] = {RESTORATION_UNITSIZE_MAX,
                                 RESTORATION_UNITSIZE_MAX * 3 / 2};

static const int32_t wins[3] = {WIENER_WIN_3TAP, WIENER_WIN_CHROMA, WIENER_WIN};

static void init_data(uint8_t **dgd, int32_t *dgd_stride, uint8_t **src,
                      int32_t *src_stride, const int idx) {
    *dgd_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *src_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *dgd = (uint8_t *)malloc(sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX *
                             *dgd_stride);
    *src = (uint8_t *)malloc(sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX *
                             *src_stride);

    if (!idx) {
        memset(*dgd,
               0,
               sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (1 == idx) {
        eb_buf_random_u8_to_255(*dgd,
                                2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        eb_buf_random_u8_to_255(*src,
                                2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (2 == idx) {
        eb_buf_random_u8_to_255(*dgd,
                                2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (3 == idx) {
        memset(*dgd,
               0,
               sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        eb_buf_random_u8_to_255(*src,
                                2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (4 == idx) {
        eb_buf_random_u8_to_0_or_255(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        eb_buf_random_u8_to_0_or_255(
            *src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else {
        eb_buf_random_u8(*dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        eb_buf_random_u8(*src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    }
}

static void init_data_highbd(uint16_t **dgd, int32_t *dgd_stride,
                             uint16_t **src, int32_t *src_stride,
                             const AomBitDepth bd, const int idx) {
    *dgd_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *src_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *dgd = (uint16_t *)malloc(sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX *
                              *dgd_stride);
    *src = (uint16_t *)malloc(sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX *
                              *src_stride);

    if (!idx) {
        memset(*dgd,
               0,
               sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (1 == idx) {
        eb_buf_random_u16_to_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        eb_buf_random_u16_to_bd(
            *src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, bd);
    } else if (2 == idx) {
        eb_buf_random_u16_to_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (3 == idx) {
        memset(*dgd,
               0,
               sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride);
        eb_buf_random_u16_to_bd(
            *src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, bd);
    } else if (4 == idx) {
        eb_buf_random_u16_to_0_or_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        eb_buf_random_u16_to_0_or_bd(
            *src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, bd);
    } else if (5 == idx) {
        // Trigger the 32-bit overflow in Step 3 and 4 for AOM_BITS_12.
        eb_buf_random_u16_to_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        for (int i = 0; i < 2 * RESTORATION_UNITSIZE_MAX; i++) {
            memset(*dgd + i * *dgd_stride, 0, sizeof(**dgd) * 20);
        }
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else if (6 == idx) {
        // Trigger the 32-bit overflow in Step 5 and 6 for AOM_BITS_12.
        eb_buf_random_u16_to_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        memset(*dgd, 0, sizeof(**dgd) * 2 * RESTORATION_UNITSIZE_MAX * 20);
        memset(*src,
               0,
               sizeof(**src) * 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    } else {
        eb_buf_random_u16_with_bd(
            *dgd, 2 * RESTORATION_UNITSIZE_MAX * *dgd_stride, bd);
        eb_buf_random_u16_with_bd(
            *src, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, bd);
    }
}

static void uninit_data(uint8_t *dgd, uint8_t *src) {
    free(dgd);
    free(src);
}

static void uninit_data_highbd(uint16_t *dgd, uint16_t *src) {
    free(dgd);
    free(src);
}

static void init_output(int64_t M_org[WIENER_WIN2], int64_t M_opt[WIENER_WIN2],
                        int64_t H_org[WIENER_WIN2 * WIENER_WIN2],
                        int64_t H_opt[WIENER_WIN2 * WIENER_WIN2]) {
    eb_buf_random_s64(M_org, WIENER_WIN2);
    eb_buf_random_s64(M_opt, WIENER_WIN2);
    eb_buf_random_s64(H_org, WIENER_WIN2 * WIENER_WIN2);
    eb_buf_random_s64(H_opt, WIENER_WIN2 * WIENER_WIN2);
    memcpy(M_opt, M_org, sizeof(*M_org) * WIENER_WIN2);
    memcpy(H_opt, H_org, sizeof(*H_org) * WIENER_WIN2 * WIENER_WIN2);
}

TEST(EbRestorationPick, compute_stats) {
    uint8_t *dgd, *src;
    int32_t dgd_stride, src_stride;
    int64_t M_org[WIENER_WIN2], M_opt[WIENER_WIN2];
    int64_t H_org[WIENER_WIN2 * WIENER_WIN2], H_opt[WIENER_WIN2 * WIENER_WIN2];

    init_output(M_org, M_opt, H_org, H_opt);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
            init_data(&dgd, &dgd_stride, &src, &src_stride, j);
            uint8_t *const d = dgd + WIENER_WIN * dgd_stride;
            uint8_t *const s = src + WIENER_WIN * src_stride;

            for (int wn_filter_mode = 0; wn_filter_mode <= 2;
                 wn_filter_mode++) {
                for (int plane = 0; plane <= 1; plane++) {
                    const int32_t wn_luma = wn_filter_mode == 1
                                                ? WIENER_WIN_3TAP
                                                : wn_filter_mode == 2
                                                      ? WIENER_WIN_CHROMA
                                                      : WIENER_WIN;
                    const int32_t wiener_win = wn_filter_mode == 1
                                                   ? WIENER_WIN_3TAP
                                                   : (plane == AOM_PLANE_Y)
                                                         ? wn_luma
                                                         : WIENER_WIN_CHROMA;

                    for (int start = 0; start <= 2; start += 2) {
                        const int32_t h_start = start;
                        const int32_t v_start = start;

                        for (int end = 0; end < 32; end += 2) {
                            const int32_t h_end =
                                RESTORATION_UNITSIZE_MAX * 3 / 2 + end;
                            const int32_t v_end =
                                RESTORATION_UNITSIZE_MAX * 3 / 2 + end;

                            av1_compute_stats_c(wiener_win,
                                                d,
                                                s,
                                                h_start,
                                                h_end,
                                                v_start,
                                                v_end,
                                                dgd_stride,
                                                src_stride,
                                                M_org,
                                                H_org);
                            av1_compute_stats_avx2(wiener_win,
                                                   d,
                                                   s,
                                                   h_start,
                                                   h_end,
                                                   v_start,
                                                   v_end,
                                                   dgd_stride,
                                                   src_stride,
                                                   M_opt,
                                                   H_opt);

                            for (int idx = 0; idx < WIENER_WIN2; idx++) {
                                if (M_org[idx] != M_opt[idx]) {
                                    printf(
                                        "\n M: width = %d, height = %d, "
                                        "wiener_win=%d, idx=%d, [%d, %d]",
                                        h_end - h_start,
                                        v_end - v_start,
                                        wiener_win,
                                        idx,
                                        idx / (wiener_win * wiener_win),
                                        idx % (wiener_win * wiener_win));
                                }
                                EXPECT_EQ(M_org[idx], M_opt[idx]);
                            }
                            for (int idx = 0; idx < WIENER_WIN2 * WIENER_WIN2;
                                 idx++) {
                                if (H_org[idx] != H_opt[idx]) {
                                    printf(
                                        "\n H: width = %d, height = %d, "
                                        "wiener_win=%d, idx=%d, [%d, %d]",
                                        h_end - h_start,
                                        v_end - v_start,
                                        wiener_win,
                                        idx,
                                        idx / (wiener_win * wiener_win),
                                        idx % (wiener_win * wiener_win));
                                }
                                EXPECT_EQ(H_org[idx], H_opt[idx]);
                            }
                        }
                    }
                }
            }

            uninit_data(dgd, src);
        }
    }
}

TEST(EbRestorationPick, compute_stats_highbd) {
    uint16_t *dgd, *src;
    int32_t dgd_stride, src_stride;
    int64_t M_org[WIENER_WIN2], M_opt[WIENER_WIN2];
    int64_t H_org[WIENER_WIN2 * WIENER_WIN2], H_opt[WIENER_WIN2 * WIENER_WIN2];

    init_output(M_org, M_opt, H_org, H_opt);

    for (int bd = AOM_BITS_8; bd <= AOM_BITS_12; bd += 2) {
        const AomBitDepth bit_depth = (AomBitDepth)bd;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 8; j++) {
                init_data_highbd(
                    &dgd, &dgd_stride, &src, &src_stride, bit_depth, j);
                const uint16_t *const d =
                    dgd + WIENER_WIN * dgd_stride + WIENER_WIN;
                const uint16_t *const s =
                    src + WIENER_WIN * src_stride + WIENER_WIN;
                const uint8_t *const dgd8 = CONVERT_TO_BYTEPTR(d);
                const uint8_t *const src8 = CONVERT_TO_BYTEPTR(s);

                for (int wn_filter_mode = 0; wn_filter_mode <= 2;
                     wn_filter_mode++) {
                    for (int plane = 0; plane <= 1; plane++) {
                        const int32_t wn_luma = wn_filter_mode == 1
                                                    ? WIENER_WIN_3TAP
                                                    : wn_filter_mode == 2
                                                          ? WIENER_WIN_CHROMA
                                                          : WIENER_WIN;
                        const int32_t wiener_win =
                            wn_filter_mode == 1
                                ? WIENER_WIN_3TAP
                                : (plane == AOM_PLANE_Y) ? wn_luma
                                                         : WIENER_WIN_CHROMA;

                        for (int size_mode = 0; size_mode < 2; size_mode++) {
                            for (int start = 0; start <= 2; start += 2) {
                                const int32_t h_start = start;
                                const int32_t v_start = start;

                                for (int end = 0; end <= 2; end += 2) {
                                    const int32_t h_end =
                                        sizes[size_mode] + end;
                                    const int32_t v_end =
                                        sizes[size_mode] + end;

                                    av1_compute_stats_highbd_c(wiener_win,
                                                               dgd8,
                                                               src8,
                                                               h_start,
                                                               h_end,
                                                               v_start,
                                                               v_end,
                                                               dgd_stride,
                                                               src_stride,
                                                               M_org,
                                                               H_org,
                                                               bit_depth);
                                    av1_compute_stats_highbd_avx2(wiener_win,
                                                                  dgd8,
                                                                  src8,
                                                                  h_start,
                                                                  h_end,
                                                                  v_start,
                                                                  v_end,
                                                                  dgd_stride,
                                                                  src_stride,
                                                                  M_opt,
                                                                  H_opt,
                                                                  bit_depth);

                                    for (int idx = 0; idx < WIENER_WIN2;
                                         idx++) {
                                        if (M_org[idx] != M_opt[idx]) {
                                            printf(
                                                "\n M: bd=%d, width = %d, "
                                                "height = %d, wiener_win=%d, "
                                                "idx=%d, [%d, %d]",
                                                bd,
                                                h_end - h_start,
                                                v_end - v_start,
                                                wiener_win,
                                                idx,
                                                idx / (wiener_win * wiener_win),
                                                idx %
                                                    (wiener_win * wiener_win));
                                        }
                                        EXPECT_EQ(M_org[idx], M_opt[idx]);
                                    }
                                    for (int idx = 0;
                                         idx < WIENER_WIN2 * WIENER_WIN2;
                                         idx++) {
                                        if (H_org[idx] != H_opt[idx]) {
                                            printf(
                                                "\n H: bd=%d, width = %d, "
                                                "height = %d, wiener_win=%d, "
                                                "idx=%d, [%d, %d]",
                                                bd,
                                                h_end - h_start,
                                                v_end - v_start,
                                                wiener_win,
                                                idx,
                                                idx / (wiener_win * wiener_win),
                                                idx %
                                                    (wiener_win * wiener_win));
                                        }
                                        EXPECT_EQ(H_org[idx], H_opt[idx]);
                                    }
                                }
                            }
                        }
                    }
                }

                uninit_data_highbd(dgd, src);
            }
        }
    }
}

TEST(EbRestorationPick, DISABLED_compute_stats_speed) {
    uint8_t *dgd, *src;
    int32_t dgd_stride, src_stride;
    int64_t M_org[WIENER_WIN2], M_opt[WIENER_WIN2];
    int64_t H_org[WIENER_WIN2 * WIENER_WIN2], H_opt[WIENER_WIN2 * WIENER_WIN2];
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    init_output(M_org, M_opt, H_org, H_opt);
    init_data(&dgd, &dgd_stride, &src, &src_stride, 5);
    uint8_t *const d = dgd + WIENER_WIN * dgd_stride;
    uint8_t *const s = src + WIENER_WIN * src_stride;
    const int32_t h_start = 0;
    const int32_t v_start = 0;
    const int32_t h_end = RESTORATION_UNITSIZE_MAX;
    const int32_t v_end = RESTORATION_UNITSIZE_MAX;

    for (int j = 0; j < 3; j++) {
        const int32_t wiener_win = wins[j];
        const uint64_t num_loop = 100000 / (wiener_win * wiener_win);

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            av1_compute_stats_c(wiener_win,
                                d,
                                s,
                                h_start,
                                h_end,
                                v_start,
                                v_end,
                                dgd_stride,
                                src_stride,
                                M_org,
                                H_org);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            av1_compute_stats_avx2(wiener_win,
                                   d,
                                   s,
                                   h_start,
                                   h_end,
                                   v_start,
                                   v_end,
                                   dgd_stride,
                                   src_stride,
                                   M_opt,
                                   H_opt);
        }

        EbStartTime(&finish_time_seconds, &finish_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);
        EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &time_o);

        for (int idx = 0; idx < WIENER_WIN2; idx++) {
            EXPECT_EQ(M_org[idx], M_opt[idx]);
        }
        for (int idx = 0; idx < WIENER_WIN2 * WIENER_WIN2; idx++) {
            EXPECT_EQ(H_org[idx], H_opt[idx]);
        }

        printf("Average Nanoseconds per Function Call\n");
        printf("    av1_compute_stats_avx2(%d)   : %6.2f\n",
               wiener_win,
               1000000 * time_c / num_loop);
        printf(
            "    av1_compute_stats_avx512(%d) : %6.2f   (Comparison: "
            "%5.2fx)\n",
            wiener_win,
            1000000 * time_o / num_loop,
            time_c / time_o);
    }
}

TEST(EbRestorationPick, DISABLED_compute_stats_highbd_speed) {
    uint16_t *dgd, *src;
    int32_t dgd_stride, src_stride;
    int64_t M_org[WIENER_WIN2], M_opt[WIENER_WIN2];
    int64_t H_org[WIENER_WIN2 * WIENER_WIN2], H_opt[WIENER_WIN2 * WIENER_WIN2];
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    init_output(M_org, M_opt, H_org, H_opt);
    const int32_t h_start = 0;
    const int32_t v_start = 0;
    const int32_t h_end = RESTORATION_UNITSIZE_MAX;
    const int32_t v_end = RESTORATION_UNITSIZE_MAX;

    for (int bd = AOM_BITS_8; bd <= AOM_BITS_12; bd += 2) {
        const AomBitDepth bit_depth = (AomBitDepth)bd;
        init_data_highbd(&dgd, &dgd_stride, &src, &src_stride, bit_depth, 7);
        uint16_t *const d = dgd + WIENER_WIN * dgd_stride;
        uint16_t *const s = src + WIENER_WIN * src_stride;
        const uint8_t *const dgd8 = CONVERT_TO_BYTEPTR(d);
        const uint8_t *const src8 = CONVERT_TO_BYTEPTR(s);

        for (int j = 0; j < 3; j++) {
            const int32_t wiener_win = wins[j];
            const uint64_t num_loop = 1000000 / (wiener_win * wiener_win);

            EbStartTime(&start_time_seconds, &start_time_useconds);

            for (uint64_t i = 0; i < num_loop; i++) {
                av1_compute_stats_highbd_c(wiener_win,
                                           dgd8,
                                           src8,
                                           h_start,
                                           h_end,
                                           v_start,
                                           v_end,
                                           dgd_stride,
                                           src_stride,
                                           M_org,
                                           H_org,
                                           bit_depth);
            }

            EbStartTime(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t i = 0; i < num_loop; i++) {
                av1_compute_stats_highbd_avx2(wiener_win,
                                              dgd8,
                                              src8,
                                              h_start,
                                              h_end,
                                              v_start,
                                              v_end,
                                              dgd_stride,
                                              src_stride,
                                              M_opt,
                                              H_opt,
                                              bit_depth);
            }

            EbStartTime(&finish_time_seconds, &finish_time_useconds);
            EbComputeOverallElapsedTimeMs(start_time_seconds,
                                          start_time_useconds,
                                          middle_time_seconds,
                                          middle_time_useconds,
                                          &time_c);
            EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                          middle_time_useconds,
                                          finish_time_seconds,
                                          finish_time_useconds,
                                          &time_o);

            for (int idx = 0; idx < WIENER_WIN2; idx++) {
                EXPECT_EQ(M_org[idx], M_opt[idx]);
            }
            for (int idx = 0; idx < WIENER_WIN2 * WIENER_WIN2; idx++) {
                EXPECT_EQ(H_org[idx], H_opt[idx]);
            }

            printf("Average Nanoseconds per Function Call\n");
            printf("    av1_compute_stats_highbd_avx2(%d)   : %6.2f\n",
                   wiener_win,
                   1000000 * time_c / num_loop);
            printf(
                "    av1_compute_stats_highbd_avx512(%d) : %6.2f   "
                "(Comparison: "
                "%5.2fx)\n",
                wiener_win,
                1000000 * time_o / num_loop,
                time_c / time_o);
        }
    }
}

#if 0
// To test integral_images() and integral_images_highbd(), make them not static,
// and add declarations to header file.
static void init_data_integral_images(uint8_t **src8, uint16_t **src16,
                                      int32_t *src_stride) {
    *src_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *src8 = (uint8_t *)malloc(sizeof(**src8) * 2 * RESTORATION_UNITSIZE_MAX *
                              *src_stride);
    *src16 = (uint16_t *)malloc(sizeof(**src16) * 2 * RESTORATION_UNITSIZE_MAX *
                                *src_stride);
    eb_buf_random_u8(*src8, 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    eb_buf_random_u16_with_bd(
        *src16, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, 12);
}

static void uninit_data_integral_images(uint8_t *src8, uint16_t *src16) {
    free(src8);
    free(src16);
}

static void integral_images_c(const uint8_t *src, int32_t src_stride,
                              int32_t width, int32_t height, int32_t *A,
                              int32_t *B, int32_t buf_stride) {
    for (int x = 0; x < width; x++) {
        A[x + 1] = 0;
        B[x + 1] = 0;
    }

    for (int y = 0; y < height; y++) {
        A[(y + 1) * buf_stride] = 0;
        B[(y + 1) * buf_stride] = 0;
        for (int x = 0; x < width; x++) {
            A[(y + 1) * buf_stride + x + 1] =
                A[(y + 1) * buf_stride + x] +
                src[y * src_stride + x] * src[y * src_stride + x];
            B[(y + 1) * buf_stride + x + 1] =
                B[(y + 1) * buf_stride + x] + src[y * src_stride + x];
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            A[(y + 1) * buf_stride + x + 1] += A[y * buf_stride + x + 1];
            B[(y + 1) * buf_stride + x + 1] += B[y * buf_stride + x + 1];
        }
    }
}

static void integral_images_highbd_c(const uint16_t *src, int32_t src_stride,
                                     int32_t width, int32_t height, int32_t *A,
                                     int32_t *B, int32_t buf_stride) {
    for (int x = 0; x < width; x++) {
        A[x + 1] = 0;
        B[x + 1] = 0;
    }

    for (int y = 0; y < height; y++) {
        A[(y + 1) * buf_stride] = 0;
        B[(y + 1) * buf_stride] = 0;
        for (int x = 0; x < width; x++) {
            A[(y + 1) * buf_stride + x + 1] =
                A[(y + 1) * buf_stride + x] +
                src[y * src_stride + x] * src[y * src_stride + x];
            B[(y + 1) * buf_stride + x + 1] =
                B[(y + 1) * buf_stride + x] + src[y * src_stride + x];
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            A[(y + 1) * buf_stride + x + 1] += A[y * buf_stride + x + 1];
            B[(y + 1) * buf_stride + x + 1] += B[y * buf_stride + x + 1];
        }
    }
}

TEST(EbRestorationPick, integral_images) {
    uint8_t *src8;
    uint16_t *src16;
    int32_t src_stride;

    // The ALIGN_POWER_OF_TWO macro here ensures that column 1 of Atl, Btl,
    // Ctl and Dtl is 32-byte aligned.
    const int32_t buf_elts = ALIGN_POWER_OF_TWO(RESTORATION_PROC_UNIT_PELS, 3);
    int32_t *buf_c, *buf_o;

    buf_c = (int32_t *)aom_memalign(32, sizeof(*buf_c) * 4 * buf_elts);
    buf_o = (int32_t *)aom_memalign(32, sizeof(*buf_o) * 4 * buf_elts);

    for (int i = 0; i < 10; i++) {
        init_data_integral_images(&src8, &src16, &src_stride);

        for (int height = 2; height < 32; height += 2) {
            for (int width = 2; width < 32; width += 2) {
                // From av1_selfguided_restoration_avx2()
                const int32_t width_ext = width + 2 * SGRPROJ_BORDER_HORZ;
                const int32_t height_ext = height + 2 * SGRPROJ_BORDER_VERT;

                // Adjusting the stride of A and B here appears to avoid bad
                // cache effects, leading to a significant speed
                // improvement. We also align the stride to a multiple of 32
                // bytes for efficiency.
                int32_t buf_stride = ALIGN_POWER_OF_TWO(width_ext + 16, 3);

                // The "tl" pointers point at the top-left of the
                // initialised data for the array.
                int32_t *Atl_c = buf_c + 0 * buf_elts + 7;
                int32_t *Btl_c = buf_c + 1 * buf_elts + 7;
                int32_t *Ctl_c = buf_c + 2 * buf_elts + 7;
                int32_t *Dtl_c = buf_c + 3 * buf_elts + 7;
                int32_t *Atl_o = buf_o + 0 * buf_elts + 7;
                int32_t *Btl_o = buf_o + 1 * buf_elts + 7;
                int32_t *Ctl_o = buf_o + 2 * buf_elts + 7;
                int32_t *Dtl_o = buf_o + 3 * buf_elts + 7;

                // The "0" pointers are (- SGRPROJ_BORDER_VERT,
                // -SGRPROJ_BORDER_HORZ). Note there's a zero row and column
                // in A, B (integral images), so we move down and right one
                // for them.
                const int32_t buf_diag_border =
                    SGRPROJ_BORDER_HORZ + buf_stride * SGRPROJ_BORDER_VERT;

                int32_t *A0_c = Atl_c + 1 + buf_stride;
                int32_t *B0_c = Btl_c + 1 + buf_stride;
                int32_t *C0_c = Ctl_c + 1 + buf_stride;
                int32_t *D0_c = Dtl_c + 1 + buf_stride;
                int32_t *A0_o = Atl_o + 1 + buf_stride;
                int32_t *B0_o = Btl_o + 1 + buf_stride;
                int32_t *C0_o = Ctl_o + 1 + buf_stride;
                int32_t *D0_o = Dtl_o + 1 + buf_stride;

                // Finally, A, B, C, D point at position (0, 0).
                int32_t *A_c = A0_c + buf_diag_border;
                int32_t *B_c = B0_c + buf_diag_border;
                int32_t *C_c = C0_c + buf_diag_border;
                int32_t *D_c = D0_c + buf_diag_border;
                int32_t *A_o = A0_o + buf_diag_border;
                int32_t *B_o = B0_o + buf_diag_border;
                int32_t *C_o = C0_o + buf_diag_border;
                int32_t *D_o = D0_o + buf_diag_border;

                (void)A_c;
                (void)A_o;
                (void)B_c;
                (void)B_o;
                (void)C_c;
                (void)C_o;
                (void)D_c;
                (void)D_o;

                const int32_t dgd_diag_border =
                    SGRPROJ_BORDER_HORZ + src_stride * SGRPROJ_BORDER_VERT;
                // Done from av1_selfguided_restoration_avx2()

                integral_images_c(src8,
                                  src_stride,
                                  width_ext,
                                  height_ext,
                                  Ctl_c,
                                  Dtl_c,
                                  buf_stride);

                integral_images(src8,
                                src_stride,
                                width_ext,
                                height_ext,
                                Ctl_o,
                                Dtl_o,
                                buf_stride);

                for (int y = 0; y < height_ext; y++) {
                    for (int x = 0; x < width_ext; x++) {
                        const int idx = y * buf_stride + x;
                        EXPECT_EQ(C0_c[idx], C0_o[idx]);
                        EXPECT_EQ(D0_c[idx], D0_o[idx]);
                    }
                }

                integral_images_highbd_c(src16,
                                         src_stride,
                                         width_ext,
                                         height_ext,
                                         Ctl_c,
                                         Dtl_c,
                                         buf_stride);

                integral_images_highbd(src16,
                                       src_stride,
                                       width_ext,
                                       height_ext,
                                       Ctl_o,
                                       Dtl_o,
                                       buf_stride);

                for (int y = 0; y < height_ext; y++) {
                    for (int x = 0; x < width_ext; x++) {
                        const int idx = y * buf_stride + x;
                        EXPECT_EQ(C0_c[idx], C0_o[idx]);
                        EXPECT_EQ(D0_c[idx], D0_o[idx]);
                    }
                }
            }
        }

        uninit_data_integral_images(src8, src16);
    }

    aom_free(buf_c);
    aom_free(buf_o);
}

TEST(EbRestorationPick, DISABLED_integral_images_speed) {
    uint8_t *src8;
    uint16_t *src16;
    int32_t src_stride;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    // The ALIGN_POWER_OF_TWO macro here ensures that column 1 of Atl, Btl,
    // Ctl and Dtl is 32-byte aligned.
    const int32_t buf_elts = ALIGN_POWER_OF_TWO(RESTORATION_PROC_UNIT_PELS, 3);
    int32_t *buf_c, *buf_o;

    buf_c = (int32_t *)aom_memalign(32, sizeof(*buf_c) * 4 * buf_elts);
    buf_o = (int32_t *)aom_memalign(32, sizeof(*buf_o) * 4 * buf_elts);

    init_data_integral_images(&src8, &src16, &src_stride);
    const int32_t width = 32;
    const int32_t height = 32;

    // From av1_selfguided_restoration_avx2()
    const int32_t width_ext = width + 2 * SGRPROJ_BORDER_HORZ;
    const int32_t height_ext = height + 2 * SGRPROJ_BORDER_VERT;

    // Adjusting the stride of A and B here appears to avoid bad
    // cache effects, leading to a significant speed
    // improvement. We also align the stride to a multiple of 32
    // bytes for efficiency.
    int32_t buf_stride = ALIGN_POWER_OF_TWO(width_ext + 16, 3);

    // The "tl" pointers point at the top-left of the
    // initialised data for the array.
    int32_t *Atl_c = buf_c + 0 * buf_elts + 7;
    int32_t *Btl_c = buf_c + 1 * buf_elts + 7;
    int32_t *Ctl_c = buf_c + 2 * buf_elts + 7;
    int32_t *Dtl_c = buf_c + 3 * buf_elts + 7;
    int32_t *Atl_o = buf_o + 0 * buf_elts + 7;
    int32_t *Btl_o = buf_o + 1 * buf_elts + 7;
    int32_t *Ctl_o = buf_o + 2 * buf_elts + 7;
    int32_t *Dtl_o = buf_o + 3 * buf_elts + 7;

    // The "0" pointers are (- SGRPROJ_BORDER_VERT,
    // -SGRPROJ_BORDER_HORZ). Note there's a zero row and column
    // in A, B (integral images), so we move down and right one
    // for them.
    const int32_t buf_diag_border =
        SGRPROJ_BORDER_HORZ + buf_stride * SGRPROJ_BORDER_VERT;

    int32_t *A0_c = Atl_c + 1 + buf_stride;
    int32_t *B0_c = Btl_c + 1 + buf_stride;
    int32_t *C0_c = Ctl_c + 1 + buf_stride;
    int32_t *D0_c = Dtl_c + 1 + buf_stride;
    int32_t *A0_o = Atl_o + 1 + buf_stride;
    int32_t *B0_o = Btl_o + 1 + buf_stride;
    int32_t *C0_o = Ctl_o + 1 + buf_stride;
    int32_t *D0_o = Dtl_o + 1 + buf_stride;

    // Finally, A, B, C, D point at position (0, 0).
    int32_t *A_c = A0_c + buf_diag_border;
    int32_t *B_c = B0_c + buf_diag_border;
    int32_t *C_c = C0_c + buf_diag_border;
    int32_t *D_c = D0_c + buf_diag_border;
    int32_t *A_o = A0_o + buf_diag_border;
    int32_t *B_o = B0_o + buf_diag_border;
    int32_t *C_o = C0_o + buf_diag_border;
    int32_t *D_o = D0_o + buf_diag_border;

    (void)A_c;
    (void)A_o;
    (void)B_c;
    (void)B_o;
    (void)C_c;
    (void)C_o;
    (void)D_c;
    (void)D_o;

    const int32_t dgd_diag_border =
        SGRPROJ_BORDER_HORZ + src_stride * SGRPROJ_BORDER_VERT;
    // Done from av1_selfguided_restoration_avx2()

    const uint64_t num_loop = 1000000;

    EbStartTime(&start_time_seconds, &start_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_c(
            src8, src_stride, width_ext, height_ext, Ctl_c, Dtl_c, buf_stride);
    }

    EbStartTime(&middle_time_seconds, &middle_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images(
            src8, src_stride, width_ext, height_ext, Ctl_o, Dtl_o, buf_stride);
    }

    EbStartTime(&finish_time_seconds, &finish_time_useconds);
    EbComputeOverallElapsedTimeMs(start_time_seconds,
                                  start_time_useconds,
                                  middle_time_seconds,
                                  middle_time_useconds,
                                  &time_c);
    EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                  middle_time_useconds,
                                  finish_time_seconds,
                                  finish_time_useconds,
                                  &time_o);

    for (int y = 0; y < height_ext; y++) {
        for (int x = 0; x < width_ext; x++) {
            const int idx = y * buf_stride + x;
            EXPECT_EQ(C0_c[idx], C0_o[idx]);
            EXPECT_EQ(D0_c[idx], D0_o[idx]);
        }
    }

    printf("Average Nanoseconds per Function Call\n");
    printf("    integral_images_c    : %6.2f\n", 1000000 * time_c / num_loop);
    printf("    integral_images_avx2 : %6.2f   (Comparison: %5.2fx)\n",
           1000000 * time_o / num_loop,
           time_c / time_o);

    EbStartTime(&start_time_seconds, &start_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_highbd_c(
            src16, src_stride, width_ext, height_ext, Ctl_c, Dtl_c, buf_stride);
    }

    EbStartTime(&middle_time_seconds, &middle_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_highbd(
            src16, src_stride, width_ext, height_ext, Ctl_o, Dtl_o, buf_stride);
    }

    EbStartTime(&finish_time_seconds, &finish_time_useconds);
    EbComputeOverallElapsedTimeMs(start_time_seconds,
                                  start_time_useconds,
                                  middle_time_seconds,
                                  middle_time_useconds,
                                  &time_c);
    EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                  middle_time_useconds,
                                  finish_time_seconds,
                                  finish_time_useconds,
                                  &time_o);

    for (int y = 0; y < height_ext; y++) {
        for (int x = 0; x < width_ext; x++) {
            const int idx = y * buf_stride + x;
            EXPECT_EQ(C0_c[idx], C0_o[idx]);
            EXPECT_EQ(D0_c[idx], D0_o[idx]);
        }
    }

    printf("Average Nanoseconds per Function Call\n");
    printf("    integral_images_highbd_c    : %6.2f\n",
           1000000 * time_c / num_loop);
    printf("    integral_images_highbd_avx2 : %6.2f   (Comparison: %5.2fx)\n",
           1000000 * time_o / num_loop,
           time_c / time_o);

    uninit_data_integral_images(src8, src16);
    aom_free(buf_c);
    aom_free(buf_o);
}
#endif
