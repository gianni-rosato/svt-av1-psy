/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#ifndef RateControlGopInfo_h
#define RateControlGopInfo_h

struct  RateControlGopInfo_s;

/*
 * @struct Holds rate control information for a group of picture
 */
typedef struct  RateControlGopInfo_s {
    /*
     * @variable uint8_t. Represents an intra if not zero.
     */
    uint8_t     exists;

    /*
     * @variable uint64_t. Frame index.
     */
    uint64_t    index;

    /*
     * @variable uint32_t. Estimated size for the frame in bytes.
     */
    uint32_t    desired_size;

    /*
     * @variable uint32_t. Actual size of the frame once encoded.
     */
    uint32_t    actual_size;

    /*
     * @variable uint8_t. Number of encoded frames in this GOP.
     */
    uint8_t     reported_frames;

    /*
     * @variable uint8_t. Number of frames in the GOP.
     */
    uint8_t     length;

    /*
     * @variable uint8_t. Assigned QP.
     */
    uint8_t     qp;

    /*
     * @variable int32_t. Variation from the model taken into account when the intra for this GOP started encoding.
     */
    int32_t     model_variation;
} RateControlGopInfo_t;

/*
 * @function get_gop_infos. Retreive the gop structure of the frame at the given
 * position from the list of GOP already recorded. Because intra frames are not
 * guaranteed to be processed in order, the previous GOP returned at the
 * moment of the call *may* not be the actual previous one in the final
 * bitsteam
 * @param {RateControlGopInfo_t*} gopInfo. Typically RateControlModel->gopInfos
 * @param {uint64_t} position. Position of the GOP to start the search from.
 * @return {RateControlGopInfo_t*}. Pointer to the GOP structure.
 * or EB_NULL if not found (unlikely).
 */
RateControlGopInfo_t *get_gop_infos(RateControlGopInfo_t *gop_info,
                                    uint64_t position);

#endif /* RateControlGopInfo_h */
