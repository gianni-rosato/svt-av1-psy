/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef RateControlModel_h
#define RateControlModel_h

#include "EbSequenceControlSet.h"
#include "RateControlGopInfo.h"
#include "EbObject.h"

/*
 * @struct Holds a prediction model for a sequence
 */
typedef struct    RateControlModel {
    EbDctor     dctor;
    /*
     * @variable uint32_t[64]. Predicted intra frame size for a given QP.
     * Indexed by QP. intraSizePrediction[15] references the estimated size of an
     * intra frame at QP = 15
     */
    size_t      intra_size_predictions[64];

    /*
     * @variable uint32_t[64]. Predicted inter frame size for a given QP.
     * Indexed by QP. interSizePrediction[15] references the estimated size of an
     * inter frame at QP = 15
     */
    size_t      inter_size_predictions[64];

    /*
     * @variable uint32_t. Desired bitrate set in the configuration
     */
    uint32_t    desired_bitrate;

    /*
     * @variable EB_U32. Desired frame rate set in the configuration
     */
    uint32_t    frame_rate;

    /*
     * @variable EB_S32. Intra period length set in the configuration
     */
    int32_t     intra_period;

    /*
     * @variable float. Variation from model. -2 means frames tend to be twice as small. 2 means twice bigger than expeted.
     */
    float       model_variation;

    /*
     * @variable int32_t. Number of variation from model calculated so far
     */
    int32_t     model_variation_reported;

    /*
     * @variable uint32_t. Sum of the bytes of all the encoded frames
     */
    uint64_t    total_bytes;

    /*
     * @variable uint64_t. Number of frames encoded so far
     */
    uint64_t    reported_frames;

    /*
     * @variable uint32_t. Video width in pixels
     */
    uint32_t    width;

    /*
     * @variable uint32_t. Video height in pixels
     */
    uint32_t    height;

    /*
     * @variable uint32_t. Number on pixels per frame. width * height
     */
    uint32_t    pixels;

    /*
     * @variable EbRateControlGopInfo[]. Information about group of picture in
     * the current sequence. Dynamically allocated in RateControlInit.
     * Indexed by pictureNumber.
     */
    EbRateControlGopInfo    *gop_infos;
} EbRateControlModel;

/*
 * @function rate_control_model_ctor. Allocate and initialize a new EbRateControlModel
 * with default values.
 * @param {EbRateControlModel*} model_ptr. Address of the pointer to EbRateControlModel
 * @return {EbErrorType}.
 */
EbErrorType rate_control_model_ctor(EbRateControlModel *model_ptr);

/*
 * @function rate_control_model_init. Initialize a model with data specific to a sequence.
 * Must be called before RateControlUpdateModel and RateControlGetQuantizer
 * @param {EbRateControlModel*} model_ptr.
 * @param {SequenceControlSet*} sequence_control_set_ptr. First frame used to initialize the model
 * @return {EbErrorType}.
 */
EbErrorType    rate_control_model_init(EbRateControlModel *model_ptr,
                                    SequenceControlSet *sequence_control_set_ptr);

/*
 * @function rate_control_update_model. Update a model with information from an encoded frame.
 * @param {EbRateControlModel*} model_ptr.
 * @param {PictureParentControlSet*} picture_ptr. Encoded frame.
 * @return {EbErrorType}.
 */
EbErrorType    rate_control_update_model(EbRateControlModel *model_ptr,
                                      PictureParentControlSet *picture_ptr);

/*
 * @function rate_control_get_quantizer. Return a QP for the given frame to be encoded.
 * Uses data from the model and information from the given frame to make a decision
 * @param {EbRateControlModel*} model_ptr.
 * @param {PictureParentControlSet*} picture_ptr. Frame to be encoded.
 * @return {uint8_t}. Suggested QP for the given frame
 */
uint8_t    rate_control_get_quantizer(EbRateControlModel *model_ptr,
                                   PictureParentControlSet *picture_ptr);

/*
 * @function get_gop_size_in_bytes. Return the size in bytes a new gop should take
 * to fit closer to the rate control constraints.
 * @param {EbRateControlModel*} model_ptr.
 * @return {uint32_t}. Ideal size the gop should take.
 */
uint32_t  get_gop_size_in_bytes(EbRateControlModel *model_ptr);

#endif // RateControlModel_h
