/*
 * Copyright(c) 2019 Netflix, Inc.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the
 * Alliance for Open Media Patent License 1.0 was not distributed with this
 * source code in the PATENTS file, you can obtain it at
 * https://www.aomedia.org/license/patent-license.
 */

/******************************************************************************
 * @file AdaptiveScanTest.cc
 *
 * @brief Unit test for adaptive scan:
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "gtest/gtest.h"
#include <stdlib.h>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "EbTransforms.h"
#include "TxfmCommon.h"
#include "EbCabacContextModel.h"  // use tx_type_to_class
#include "util.h"
#include "random.h"

static bool valid_scan(const int16_t *scan, const int16_t *iscan, int si,
                       int expected_pos) {
    if (scan[si] != expected_pos || iscan[expected_pos] != si)
        return false;
    else
        return true;
}

/**
 * @brief Unit test for scan tables:
 *
 * Test strategy:
 * Iterate the positions according to the scan type and verify the pos from
 * scan table and si from iscan table.
 *
 * Expect result:
 * Position from the scan table should be identical with the pos calculated
 * from x and y.
 *
 * Test coverage:
 * Test cases:
 * All valid scan, including zig-zag scan, vertical scan, horizontal scan
 * and rectangle scan
 *
 */
TEST(AdaptiveScanTest, scan_tables_test) {
    for (int tx_size = TX_4X4; tx_size < TX_SIZES_ALL; ++tx_size) {
        const int txb_height = get_txb_high((TxSize)tx_size);
        const int txb_width = get_txb_wide((TxSize)tx_size);
        const int bwl = get_txb_bwl((TxSize)tx_size);

        for (int tx_type = DCT_DCT; tx_type < TX_TYPES; ++tx_type) {
            if (!is_txfm_allowed((TxType)tx_type, (TxSize)tx_size))
                continue;

            TxClass tx_class = tx_type_to_class[tx_type];
            const ScanOrder *scan_order = &av1_scan_orders[tx_size][tx_type];
            const int16_t *scan = scan_order->scan;
            const int16_t *iscan = scan_order->iscan;
            const int sum_xy = txb_width + txb_height - 1;
            int si = 0;
            if (tx_class == TX_CLASS_VERT) {
                // vertical transform will scan the rows horizontally,
                // so si is the same with pos
                // process horizontal scan first
                for (int y = 0; y < txb_height; ++y) {
                    for (int x = 0; x < txb_width; ++x) {
                        EXPECT_TRUE(valid_scan(scan, iscan, si, x + (y << bwl)))
                            << "tx_size " << tx_size << " tx_type " << tx_type
                            << " scan table test fail";
                        ++si;
                    }
                }
            } else if (tx_class == TX_CLASS_HORIZ) {
                // horizontal transform will scan vertically.
                for (int x = 0; x < txb_width; ++x) {
                    for (int y = 0; y < txb_height; ++y) {
                        EXPECT_TRUE(valid_scan(scan, iscan, si, x + (y << bwl)))
                            << "tx_size " << tx_size << " tx_type " << tx_type
                            << " scan table test fail";
                        ++si;
                    }
                }
            } else if (txb_width < txb_height) {
                // cols are smaller than rows
                for (int i = 0; i < sum_xy; ++i) {
                    for (int y = 0; y < txb_height; ++y) {
                        int x = i - y;
                        if (x >= 0 && x < txb_width) {
                            EXPECT_TRUE(
                                valid_scan(scan, iscan, si, x + (y << bwl)))
                                << "tx_size " << tx_size << " tx_type "
                                << tx_type << " scan table test fail";
                            ++si;
                        }
                    }
                }
            } else if (txb_width > txb_height) {
                // col diag
                for (int i = 0; i < sum_xy; ++i) {
                    for (int x = 0; x < txb_width; ++x) {
                        int y = i - x;
                        if (y >= 0 && y < txb_height) {
                            EXPECT_TRUE(
                                valid_scan(scan, iscan, si, x + (y << bwl)))
                                << "tx_size " << tx_size << " tx_type "
                                << tx_type << " scan table test fail";
                            ++si;
                        }
                    }
                }
            } else {
                // zig-zag scan
                int swap = 1;
                for (int i = 0; i < sum_xy; ++i) {
                    swap ^= 1;
                    for (int j = 0; j < txb_width; ++j) {
                        const int x = swap ? i - j : j;
                        const int y = swap ? j : i - j;
                        if (i - j >= 0 && i - j < txb_width) {
                            EXPECT_TRUE(
                                valid_scan(scan, iscan, si, x + (y << bwl)))
                                << "tx_size " << tx_size << " tx_type "
                                << tx_type << " scan table test fail";
                            ++si;
                        }
                    }
                }
            }
        }
    }
}

using svt_av1_test_tool::SVTRandom;
TEST(CopyMiMapGrid, avx2) {
    SVTRandom rnd(0, (1 << 10) - 1);
    const int max_size = 100;
    ModeInfo *mi_grid_ref[max_size * max_size];
    ModeInfo *mi_grid_tst[max_size * max_size];

    for (int tx_size = TX_4X4; tx_size < TX_SIZES_ALL; ++tx_size) {
        const int txb_height = get_txb_high((TxSize)tx_size);
        const int txb_width = get_txb_wide((TxSize)tx_size);
        const uint32_t stride = txb_width + rnd.random() % 10;

        memset(mi_grid_ref, 0xcd, sizeof(mi_grid_ref));
        memset(mi_grid_tst, 0xcd, sizeof(mi_grid_ref));

        mi_grid_ref[0] = mi_grid_tst[0] = (ModeInfo *)((uint64_t)rnd.random());

        svt_copy_mi_map_grid_c(mi_grid_ref, stride, txb_height, txb_width);
        svt_copy_mi_map_grid_avx2(mi_grid_tst, stride, txb_height, txb_width);

        EXPECT_TRUE(memcmp(mi_grid_ref, mi_grid_tst, sizeof(mi_grid_ref)) == 0);
    }

    // special case non power of 2 cols
    for (int cols = 0; cols < 15; ++cols) {
        const uint32_t stride = cols + rnd.random() % 10;

        memset(mi_grid_ref, 0xcd, sizeof(mi_grid_ref));
        memset(mi_grid_tst, 0xcd, sizeof(mi_grid_ref));

        mi_grid_ref[0] = mi_grid_tst[0] = (ModeInfo *)((uint64_t)rnd.random());

        svt_copy_mi_map_grid_c(mi_grid_ref, stride, cols, cols);
        svt_copy_mi_map_grid_avx2(mi_grid_tst, stride, cols, cols);

        EXPECT_TRUE(memcmp(mi_grid_ref, mi_grid_tst, sizeof(mi_grid_ref)) == 0);
    }
}
