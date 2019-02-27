/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPacketization_h
#define EbPacketization_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#ifdef __cplusplus
extern "C" {
#endif

    /**************************************
     * Type Declarations
     **************************************/
    typedef struct EbPPSConfig_s
    {
        uint8_t ppsId;
        uint8_t constrainedFlag;

    } EbPPSConfig_t;

    /**************************************
     * Context
     **************************************/
    typedef struct PacketizationContext_s
    {
        EbFifo_t      *entropy_coding_input_fifo_ptr;
        EbFifo_t      *rate_control_tasks_output_fifo_ptr;
        EbPPSConfig_t *ppsConfig;

        uint64_t   dpbDispOrder[8], dpbDecOrder[8];
        uint64_t   totShownFrames;
        uint64_t   dispOrderContinuityCount;

    } PacketizationContext_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType packetization_context_ctor(
        PacketizationContext_t **context_dbl_ptr,
        EbFifo_t                *entropy_coding_input_fifo_ptr,
        EbFifo_t                *rate_control_tasks_output_fifo_ptr);



    extern void* PacketizationKernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbPacketization_h
