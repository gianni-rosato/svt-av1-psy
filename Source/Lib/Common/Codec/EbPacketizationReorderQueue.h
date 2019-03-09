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
    typedef struct PacketizationReorderEntry_s {
        uint64_t                          picture_number;
        EbObjectWrapper_t              *output_stream_wrapper_ptr;
        EbObjectWrapper_t              *outputStatisticsWrapperPtr;

        EbLinkedListNode               *outMetaData;

        uint64_t                          start_time_seconds;
        uint64_t                          start_time_u_seconds;

        uint8_t                                 slice_type;
        uint64_t                                refPOCList0;
        uint64_t                                refPOCList1;
        uint64_t                                 poc;
        FRAME_TYPE                            av1FrameType;
        Av1RpsNode_t                          av1RefSignal;
        EbBool                               showFrame;
        EbBool                               hasShowExisting;
        uint8_t                                 showExistingLoc;


    } PacketizationReorderEntry_t;

    extern EbErrorType packetization_reorder_entry_ctor(
        PacketizationReorderEntry_t **entry_dbl_ptr,
        uint32_t                      picture_number);


#ifdef __cplusplus
}
#endif
#endif //EbPacketizationReorderQueue_h
