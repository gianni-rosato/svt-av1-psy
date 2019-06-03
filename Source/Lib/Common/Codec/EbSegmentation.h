

#ifndef SVT_AV1_EBSEGMENTATION_H
#define SVT_AV1_EBSEGMENTATION_H


#include <stdint.h>
#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbRateControlProcess.h"



void ApplySegmentationBasedQuantization(
        const BlockGeom * blk_geom,
        PictureControlSet *picture_control_set_ptr,
        LargestCodingUnit *sb_ptr,
        CodingUnit *cu_ptr
);

void SetupSegmentation(
        PictureControlSet *picture_control_set_ptr,
        SequenceControlSet *sequence_control_set_ptr,
        RateControlContext *context_ptr,
        RateControlLayerContext *rateControlLayerPtr,
        RateControlIntervalParamContext *rateControlParamPtr
);

void FindSegmentQps(
        SegmentationParams *segmentation_params,
        PictureControlSet *picture_control_set_ptr
);

void TemporallyUpdateQPs(
        int32_t *segment_qp_ptr,
        int32_t *prev_segment_qp_ptr,
        EbBool temporal_update
);

void CalculateSegmentationData(
        SegmentationParams *segmentation_params
);

#endif //SVT_AV1_EBSEGMENTATIONS_H


