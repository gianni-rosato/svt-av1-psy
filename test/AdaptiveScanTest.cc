/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
