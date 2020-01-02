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

#ifndef SVT_AV1_EBSEGMENTATION_H
#define SVT_AV1_EBSEGMENTATION_H

#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbRateControlProcess.h"
#include "EbUtility.h"

void apply_segmentation_based_quantization(const BlockGeom *blk_geom, PictureControlSet *pcs_ptr,
                                           SuperBlock *sb_ptr, BlkStruct *blk_ptr);

void setup_segmentation(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                        RateControlLayerContext *rateControlLayerPtr);

void find_segment_qps(SegmentationParams *segmentation_params, PictureControlSet *pcs_ptr);

void temporally_update_qps(int32_t *segment_qp_ptr, int32_t *prev_segment_qp_ptr,
                           EbBool temporal_update);

void calculate_segmentation_data(SegmentationParams *segmentation_params);

#endif //SVT_AV1_EBSEGMENTATIONS_H
