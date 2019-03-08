/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbUtility.h"

#include "EbPictureControlSet.h"
#include "RateControlModel.h"
#include "RateControlGopInfo.h"

/*
 * @private
 * @function get_inter_qp_for_size. Helper to extract the suggested QP from the
 * prediction model for a desired size
 * @param {EbRateControlModel*} model_ptr.
 * @param {uint32_t} desired_size.
 * @return {uint32_t}.
 */
static uint32_t    get_inter_qp_for_size(EbRateControlModel *model_ptr, uint32_t desired_size);

/*
 * @private
 * @function record_new_gop. Take into account a new group of picture in the
 * model
 * @param {EbRateControlModel*} model_ptr.
 * @param {PictureParentControlSet_t*} picture_ptr. Picture holding the intra frame.
 * @return {void}.
 */
static void record_new_gop(EbRateControlModel *model_ptr, PictureParentControlSet_t *picture_ptr);

/*
 * Average size in bits for and intra frame per QP for a 1920x1080 reference video clip
 */
static const size_t DEFAULT_REF_INTRA_PICTURE_COMPRESSION_RATIO[64] = {
    14834056, 14589237, 14344418, 14099599, 13854780, 13609961, 13365142,
    13120323, 12875503, 12492037, 12108571, 11725105, 11341639, 10958173,
    10574707, 10191241, 9942206, 9693170, 9444135, 9195099, 8946064, 8697028,
    8447993, 8198957, 7949922, 7700886, 7440985, 7181083, 6921182, 6661280,
    6401379, 6141478, 5881576, 5621675, 5361773, 5101872, 4841971, 4582069,
    4322168, 4062266, 3802365, 3542464, 3282562, 3022661, 2762759, 2502858,
    2398196, 2293535, 2188873, 2084211, 1979550, 1874888, 1770226, 1665565,
    1560903, 1456241, 1351580, 1246918, 1142256, 1037595, 932933, 828271,
    723610, 618948
};

/*
 * Average size in bits for and inter frame per QP for a 1920x1080 reference video clip
 */
static const size_t DEFAULT_REF_INTER_PICTURE_COMPRESSION_RATIO[64] = {
     12459127, 10817116, 9175105, 7533094, 5891083, 4249072, 3819718, 3390366,
     2961014, 2531662, 2102309, 1931469, 1760631, 1589793, 1418955, 1248117,
     1165879, 1083646, 1001413, 919180, 836947, 783885, 730823, 677761,
     624699, 571634, 531742, 491850, 451958, 412066, 372174, 344015,
     315858, 287701, 259544, 231387, 213604, 195823, 178042, 160261, 142480,
     131278, 120077, 108876, 97675, 86474, 79681, 72892, 66103, 59314,
     52525, 48203, 43882, 39561, 35240, 30919, 28284, 25651, 23018,
     20385, 17752, 15465, 13178, 10891
};

EbErrorType    rate_control_model_ctor(EbRateControlModel **object_doubble_ptr) {
    EbRateControlModel  *model_ptr;

    EB_MALLOC(EbRateControlModel*, model_ptr, sizeof(EbRateControlModel), EB_N_PTR);
    *object_doubble_ptr = (void*)model_ptr;

    EB_MEMSET(model_ptr, 0, sizeof(EbRateControlModel));
    EB_MEMCPY(model_ptr->intra_size_predictions, (void*)DEFAULT_REF_INTRA_PICTURE_COMPRESSION_RATIO, sizeof(DEFAULT_REF_INTRA_PICTURE_COMPRESSION_RATIO));
    EB_MEMCPY(model_ptr->inter_size_predictions, (void*)DEFAULT_REF_INTER_PICTURE_COMPRESSION_RATIO, sizeof(DEFAULT_REF_INTER_PICTURE_COMPRESSION_RATIO));

    return EB_ErrorNone;
}

EbErrorType rate_control_model_init(EbRateControlModel *model_ptr, SequenceControlSet_t *sequenceControlSetPtr) {
    uint32_t                number_of_frame = sequenceControlSetPtr->static_config.frames_to_be_encoded;
    EbRateControlGopInfo    *gop_infos;

    EB_MALLOC(EbRateControlGopInfo*, gop_infos, sizeof(EbRateControlModel) * number_of_frame, EB_N_PTR);
    memset(gop_infos, 0, sizeof(EbRateControlModel) * number_of_frame);

    model_ptr->desired_bitrate = sequenceControlSetPtr->static_config.target_bit_rate;
    model_ptr->frame_rate = sequenceControlSetPtr->static_config.frame_rate >> 16;
    model_ptr->width = sequenceControlSetPtr->luma_width;
    model_ptr->height = sequenceControlSetPtr->luma_height;
    model_ptr->pixels = model_ptr->width * model_ptr->height;
    model_ptr->gop_infos = gop_infos;
    model_ptr->intra_period = sequenceControlSetPtr->static_config.intra_period_length;

    return EB_ErrorNone;
}

EbErrorType    rate_control_update_model(EbRateControlModel *model_ptr, PictureParentControlSet_t *picture_ptr) {
    uint64_t                size = picture_ptr->total_num_bits;
    EbRateControlGopInfo    *gop = get_gop_infos(model_ptr->gop_infos, picture_ptr->picture_number);

    model_ptr->total_bytes += size;
    model_ptr->reported_frames++;
    gop->actual_size += size;
    gop->reported_frames++;

    if (gop->reported_frames == gop->length) {
      float      variation = 1;

      model_ptr->model_variation_reported++;
      if (gop->actual_size > gop->desired_size) {
        variation = -((int64_t)gop->actual_size / (int64_t)gop->desired_size);
      } else if (gop->desired_size > gop->actual_size) {
        variation = ((int64_t)(gop->desired_size / (int64_t)gop->actual_size));
      }
      variation = CLIP3(-100, 100, variation);
      variation -= gop->model_variation;

      model_ptr->model_variation = model_ptr->model_variation * (model_ptr->model_variation_reported - 1) / model_ptr->model_variation_reported + variation / model_ptr->model_variation_reported;
    }
    
    return EB_ErrorNone;
}

uint8_t    rate_control_get_quantizer(EbRateControlModel *model_ptr, PictureParentControlSet_t *picture_ptr) {
    FRAME_TYPE  type = picture_ptr->av1FrameType;

    if (type == INTRA_ONLY_FRAME || type == KEY_FRAME) {
        record_new_gop(model_ptr, picture_ptr);
    }

    EbRateControlGopInfo *gop = get_gop_infos(model_ptr->gop_infos, picture_ptr->picture_number);

    return gop->qp;
}

uint32_t get_inter_qp_for_size(EbRateControlModel *model_ptr, uint32_t desired_size) {
    uint8_t     qp;
    
    for (qp = 0; qp < MAX_QP_VALUE; qp++) {
        float    size = model_ptr->inter_size_predictions[qp];

        size = (size / (1920 * 1080)) * model_ptr->pixels;
        // Scale size for current resolution
        if (desired_size > size) {
            break;
        }
    }

    return qp;
}

static void record_new_gop(EbRateControlModel *model_ptr, PictureParentControlSet_t *picture_ptr) {
    uint64_t                pictureNumber = picture_ptr->picture_number;
    EbRateControlGopInfo    *gop = &model_ptr->gop_infos[pictureNumber];

    gop->index = pictureNumber;
    gop->exists = EB_TRUE;
    gop->desired_size = get_gop_size_in_bytes(model_ptr);
    gop->model_variation = model_ptr->model_variation;

    uint32_t                size = gop->desired_size / model_ptr->intra_period;

    gop->qp = get_inter_qp_for_size(model_ptr, size);

    // Update length in gopinfos
    if (pictureNumber != 0) {
        EbRateControlGopInfo *previousGop = get_gop_infos(model_ptr->gop_infos, pictureNumber - 1);

        previousGop->length = gop->index - previousGop->index;
    } 
}

uint32_t get_gop_size_in_bytes(EbRateControlModel *model_ptr) {
    uint32_t    gop_per_second = (model_ptr->frame_rate << 8)  / model_ptr->intra_period;
    uint32_t    gop_size = ((model_ptr->desired_bitrate << 8) / gop_per_second);
    uint32_t    desired_total_bytes = (model_ptr->desired_bitrate / model_ptr->frame_rate) * model_ptr->reported_frames;
    int64_t     delta_bytes = desired_total_bytes - model_ptr->total_bytes;
    float       extra = 1;

    if (delta_bytes < 0) {
      extra = (uint64_t)((float)delta_bytes / (float)gop_size) + 1;
    } else if (delta_bytes > 0) {
      extra = (uint64_t)((float)delta_bytes / (float)gop_size) + 1;
    }

    if (model_ptr->model_variation < 0) {
      gop_size = gop_size / -(model_ptr->model_variation);
    } else if (model_ptr->model_variation > 0) {
      gop_size = gop_size * model_ptr->model_variation;
    }

    if (model_ptr->model_variation_reported > 5) {
      if (extra < 0) {
        gop_size = gop_size / -extra;
      } else if (extra > 0) {
        gop_size = gop_size * extra;
      }
    }

    return gop_size + 1;
}
