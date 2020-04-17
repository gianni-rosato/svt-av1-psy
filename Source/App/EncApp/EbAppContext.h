/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAppContext_h
#define EbAppContext_h

#include "EbSvtAv1Enc.h"
#include "EbAppConfig.h"

/***************************************

 * App Callback data struct
 ***************************************/
typedef struct EbAppContext {
    EbSvtAv1EncConfiguration eb_enc_parameters;

    // Output Ports Active Flags
    AppPortActiveType output_stream_port_active;

    // Component Handle
    EbComponentType *svt_encoder_handle;

    // Buffer Pools
    EbBufferHeaderType *input_buffer_pool;
    EbBufferHeaderType *recon_buffer;

    // Instance Index
    uint8_t instance_idx;
} EbAppContext;

/********************************
 * External Function
 ********************************/
extern EbErrorType init_encoder(EbConfig *config, EbAppContext *callback_data,
                                uint32_t instance_idx);
extern EbErrorType de_init_encoder(EbAppContext *callback_data_ptr, uint32_t instance_index);

#endif // EbAppContext_h
