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

/******************************************************************************
 * @file ParseUtil.cc
 *
 * @brief Impelmentation of sequence header parser.
 *
 ******************************************************************************/

#include "EbDefinitions.h"
#include "EbDecParseObuUtil.h"
#include "ParseUtil.h"
#include "gtest/gtest.h"

namespace svt_av1_e2e_tools {

void SequenceHeaderParser::input_obu_data(const uint8_t* obu_data,
                                          const uint32_t size,
                                          RefDecoder::StreamInfo* stream_info) {
    SeqHeader seg_header;
    if (eb_get_sequence_info(obu_data, size, &seg_header) == EB_ErrorNone) {
        profile_ = seg_header.seq_profile;
        switch (seg_header.sb_size) {
        case BLOCK_64X64: sb_size_ = 64; break;
        case BLOCK_128X128: sb_size_ = 128; break;
        default:
            ASSERT_TRUE(false)
                << "super block size is invalid in sequence header!";
            break;
        }
        printf("SPS header: profile(%u), sb_size(%u)\n", profile_, sb_size_);

        // update stream info
        stream_info->profile = seg_header.seq_profile;
        stream_info->still_pic = seg_header.still_picture == 1;
        stream_info->sb_size = sb_size_;
        stream_info->force_integer_mv = seg_header.seq_force_integer_mv;
        stream_info->enable_filter_intra = seg_header.enable_filter_intra;
        stream_info->enable_intra_edge_filter =
            seg_header.enable_intra_edge_filter;
        stream_info->enable_masked_compound = seg_header.enable_masked_compound;
        stream_info->enable_dual_filter = seg_header.enable_dual_filter;
        stream_info->enable_jnt_comp =
            seg_header.order_hint_info.enable_jnt_comp;
        stream_info->enable_ref_frame_mvs =
            seg_header.order_hint_info.enable_ref_frame_mvs;
        stream_info->enable_warped_motion = seg_header.enable_warped_motion;
        stream_info->enable_cdef = seg_header.enable_cdef;
        stream_info->enable_restoration = seg_header.enable_restoration;
        stream_info->enable_superres = seg_header.enable_superres;
    }
}
}  // namespace svt_av1_e2e_tools
