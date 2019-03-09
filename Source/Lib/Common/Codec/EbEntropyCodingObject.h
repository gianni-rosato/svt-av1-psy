/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEntropyCodingObject_h
#define EbEntropyCodingObject_h

#include "EbDefinitions.h"
#include "EbCabacContextModel.h"
#include "EbBitstreamUnit.h"
#ifdef __cplusplus
extern "C" {
#endif
    typedef struct Bitstream_s {
        EbPtr outputBitstreamPtr;
    } Bitstream_t;

    typedef struct EntropyCoder_s {
        EbPtr cabacEncodeContextPtr;
        FRAME_CONTEXT   *fc;              /* this frame entropy */
        aom_writer       ecWriter;
        EbPtr           ecOutputBitstreamPtr;
#if TILES
        uint64_t   ec_frame_size;
#endif
    } EntropyCoder_t;

    extern EbErrorType BitstreamCtor(
        Bitstream_t **bitstreamDblPtr,
        uint32_t bufferSize);


    extern EbErrorType EntropyCoderCtor(
        EntropyCoder_t **entropyCoderDblPtr,
        uint32_t bufferSize);



    extern EbPtr EntropyCoderGetBitstreamPtr(
        EntropyCoder_t *entropy_coder_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbEntropyCodingObject_h