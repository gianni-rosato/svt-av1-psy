/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPacketizationReorderQueue_h
#define EbPacketizationReorderQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPredictionStructure.h"

#ifdef __cplusplus
extern "C" {
#endif
    /************************************************
     * Packetization Reorder Queue Entry
     ************************************************/
    typedef struct PacketizationReorderEntry
    {
        uint64_t                          picture_number;
        EbObjectWrapper              *output_stream_wrapper_ptr;
        EbObjectWrapper              *outputStatisticsWrapperPtr;

        EbLinkedListNode               *out_meta_data;

        uint64_t                          start_time_seconds;
        uint64_t                          start_time_u_seconds;

        uint8_t                                 slice_type;
        uint64_t                                ref_poc_list0;
        uint64_t                                ref_poc_list1;
#if REF_ORDER
        uint64_t                                ref_poc_array[7];
#endif
        uint64_t                                 poc;
#if RC
        uint64_t                                total_num_bits;
#endif
        FrameType                            av1_frame_type;
        Av1RpsNode                          av1_ref_signal;
        EbBool                               show_frame;
        EbBool                               has_show_existing;
        uint8_t                              show_existing_loc;
#if ALT_REF_OVERLAY
        uint8_t                              is_alt_ref;
#endif
    } PacketizationReorderEntry;

    extern EbErrorType packetization_reorder_entry_ctor(
        PacketizationReorderEntry **entry_dbl_ptr,
        uint32_t                      picture_number);

#ifdef __cplusplus
}
#endif
#endif //EbPacketizationReorderQueue_h
