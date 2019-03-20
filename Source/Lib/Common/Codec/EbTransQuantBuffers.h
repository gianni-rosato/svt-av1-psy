/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTransQuantBuffers_h
#define EbTransQuantBuffers_h

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct EbTransQuantBuffers_s
    {
        EbPictureBufferDesc_t         *tuTransCoeff2Nx2NPtr;
        EbPictureBufferDesc_t         *tuTransCoeffNxNPtr;
        EbPictureBufferDesc_t         *tuTransCoeffN2xN2Ptr;
        EbPictureBufferDesc_t         *tuQuantCoeffNxNPtr;
        EbPictureBufferDesc_t         *tuQuantCoeffN2xN2Ptr;

    } EbTransQuantBuffers_t;


    extern EbErrorType EbTransQuantBuffersCtor(
        EbTransQuantBuffers_t            *trans_quant_buffers_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbTransQuantBuffers_h