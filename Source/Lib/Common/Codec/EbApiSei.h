/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAPISEI_h
#define EBAPISEI_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// EbPicturePlane defines the data formatting of a singple plane of picture data.
typedef struct EbPicturePlane
{
    // "start" is the starting position of the first
    //   valid pixel in the picture plane.
    uint8_t* start;

    // "stride" is the number of bytes required to increment
    //   an address by one line without changing the horizontal
    //   positioning.
    uint32_t stride;
} EbPicturePlane;

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // EBAPISEI_h
