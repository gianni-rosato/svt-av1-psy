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

typedef void (*av1_compute_stats_func)(int32_t wiener_win, const uint8_t *dgd8,
                                       const uint8_t *src8, int32_t h_start,
                                       int32_t h_end, int32_t v_start,
                                       int32_t v_end, int32_t dgd_stride,
                                       int32_t src_stride, int64_t *M,
                                       int64_t *H);

typedef void (*av1_compute_stats_highbd_func)(
    int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8,
    int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end,
    int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H,
    AomBitDepth bit_depth);

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

void match_test(av1_compute_stats_func func) {
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

                            eb_av1_compute_stats_c(wiener_win,
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

                            func(wiener_win,
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

void highbd_match_test(av1_compute_stats_highbd_func func) {
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

                                    eb_av1_compute_stats_highbd_c(wiener_win,
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

                                    func(wiener_win,
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

void speed_test(av1_compute_stats_func func) {
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
            eb_av1_compute_stats_c(wiener_win,
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
            func(wiener_win,
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
        printf("    eb_av1_compute_stats_c(%d)   : %6.2f\n",
               wiener_win,
               1000000 * time_c / num_loop);
        printf(
            "    eb_av1_compute_stats_opt(%d) : %6.2f   (Comparison: %5.2fx)\n",
            wiener_win,
            1000000 * time_o / num_loop,
            time_c / time_o);
    }
}

void highbd_speed_test(av1_compute_stats_highbd_func func) {
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
                eb_av1_compute_stats_highbd_c(wiener_win,
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
                func(wiener_win,
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
            printf("    eb_av1_compute_stats_highbd_c(%d)   : %6.2f\n",
                   wiener_win,
                   1000000 * time_c / num_loop);
            printf(
                "    av1_compute_stats_highbd_opt(%d) : %6.2f   "
                "(Comparison: "
                "%5.2fx)\n",
                wiener_win,
                1000000 * time_o / num_loop,
                time_c / time_o);
        }
    }
}

TEST(EbRestorationPick_avx2, match) {
    match_test(eb_av1_compute_stats_avx2);
}

TEST(EbRestorationPick_avx2, match_highbd) {
    highbd_match_test(eb_av1_compute_stats_highbd_avx2);
}

TEST(EbRestorationPick_avx2, DISABLED_speed) {
    speed_test(eb_av1_compute_stats_avx2);
}

TEST(EbRestorationPick_avx2, DISABLED_speed_highbd) {
    highbd_speed_test(eb_av1_compute_stats_highbd_avx2);
}

#ifndef NON_AVX512_SUPPORT

TEST(EbRestorationPick_avx512, match) {
    if (get_cpu_flags_to_use() & CPU_FLAGS_AVX512F) {
        match_test(eb_av1_compute_stats_avx512);
    }
}

TEST(EbRestorationPick_avx512, match_highbd) {
    if (get_cpu_flags_to_use() & CPU_FLAGS_AVX512F) {
        highbd_match_test(eb_av1_compute_stats_highbd_avx512);
    }
}

TEST(EbRestorationPick_avx512, DISABLED_speed) {
    if (get_cpu_flags_to_use() & CPU_FLAGS_AVX512F) {
        speed_test(eb_av1_compute_stats_avx512);
    }
}

TEST(EbRestorationPick_avx512, DISABLED_speed_highbd) {
    if (get_cpu_flags_to_use() & CPU_FLAGS_AVX512F) {
        highbd_speed_test(eb_av1_compute_stats_highbd_avx512);
    }
}

#endif
