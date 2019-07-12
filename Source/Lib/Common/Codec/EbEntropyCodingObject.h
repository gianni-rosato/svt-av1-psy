/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEntropyCodingObject_h
#define EbEntropyCodingObject_h

#include "EbDefinitions.h"
#include "EbCabacContextModel.h"
#include "EbBitstreamUnit.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
    typedef struct Bitstream {
        EbDctor dctor;
        EbPtr output_bitstream_ptr;
    } Bitstream;

    typedef struct EntropyCoder
    {
        EbDctor dctor;
        EbPtr cabac_encode_context_ptr;
        FRAME_CONTEXT   *fc;              /* this frame entropy */
        AomWriter       ec_writer;
        EbPtr           ec_output_bitstream_ptr;
        uint64_t   ec_frame_size;
    } EntropyCoder;

    extern EbErrorType bitstream_ctor(
        Bitstream *bitstream_ptr,
        uint32_t buffer_size);

    extern EbErrorType entropy_coder_ctor(
        EntropyCoder *entropy_coder_ptr,
        uint32_t buffer_size);

    extern EbPtr entropy_coder_get_bitstream_ptr(
        EntropyCoder *entropy_coder_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEntropyCodingObject_h
