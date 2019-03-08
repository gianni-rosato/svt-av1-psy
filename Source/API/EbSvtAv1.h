/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCodec_h
#define EbCodec_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "stdint.h"
#include "EbSvtAv1Formats.h"

    // API Version
#define SVT_VERSION_MAJOR       0
#define SVT_VERSION_MINOR       4
#define SVT_VERSION_PATCHLEVEL  0

#ifdef _WIN32
#define EB_API __declspec(dllexport)
#else
#define EB_API
#endif

#define EB_MAX_NUM_OPERATING_POINTS            32

#define EB_MAX_TEMPORAL_LAYERS              MAX_TEMPORAL_LAYERS


/********************************
* Defines
********************************/
#define EB_PICTURE           uint32_t

typedef enum EbAv1PictureType
{
    EB_AV1_INTER_PICTURE         = 0,
    EB_AV1_ALT_REF_PICTURE       = 1,
    EB_AV1_INTRA_ONLY_PICTURE    = 2,
    EB_AV1_KEY_PICTURE           = 3,
    EB_AV1_NON_REF_PICTURE       = 4,
    EB_AV1_SHOW_EXISTING_PICTURE = 6,
    EB_AV1_FW_KEY_PICTURE        = 5,
    EB_AV1_SWITCH_PICTURE        = 7,
    EB_AV1_INVALID_PICTURE       = 0xFF
} EbAv1PictureType;

/** The EbBool type is intended to be used to represent a true or a false
value when passing parameters to and from the eBrisk API.  The
EbBool is a 32 bit quantity and is aligned on a 32 bit word boundary.
*/

#define EbBool   uint8_t
#define EB_FALSE  0
#define EB_TRUE   1

typedef struct EbBufferHeaderType
{
    // EbBufferHeaderType size
    uint32_t size;

    // picture (input or output) buffer
    uint8_t* p_buffer;
    uint32_t n_filled_len;
    uint32_t n_alloc_len;

    // pic private data
    void*    p_app_private;
    void*    wrapper_ptr;

    // pic timing param
    uint32_t n_tick_count;
    int64_t  dts;
    int64_t  pts;

    // pic info
    uint32_t qp;
    uint32_t pic_type;

    // pic flags
    uint32_t flags;
} EbBufferHeaderType;

typedef struct EbComponentType
{
    uint32_t size;
    void* pComponentPrivate;
    void* pApplicationPrivate;
} EbComponentType;

typedef enum EbErrorType
{
    EB_ErrorNone = 0,

    EB_DecUnsupportedBitstream = (int32_t)0x40001000,
    EB_DecNoOutputPicture = (int32_t)0x40001004,
    EB_DecDecodingError = (int32_t)0x40001008,

    EB_ErrorInsufficientResources = (int32_t)0x80001000,
    EB_ErrorUndefined = (int32_t)0x80001001,
    EB_ErrorInvalidComponent = (int32_t)0x80001004,
    EB_ErrorBadParameter = (int32_t)0x80001005,
    EB_ErrorDestroyThreadFailed = (int32_t)0x80002012,
    EB_ErrorSemaphoreUnresponsive = (int32_t)0x80002021,
    EB_ErrorDestroySemaphoreFailed = (int32_t)0x80002022,
    EB_ErrorCreateMutexFailed = (int32_t)0x80002030,
    EB_ErrorMutexUnresponsive = (int32_t)0x80002031,
    EB_ErrorDestroyMutexFailed = (int32_t)0x80002032,
    EB_NoErrorEmptyQueue = (int32_t)0x80002033,

    EB_ErrorMax = 0x7FFFFFFF
} EbErrorType;

/* AV1 bistream profile (seq_profile syntax element) */
typedef enum EbAv1SeqProfile 
{
    MAIN_PROFILE = 0,
    HIGH_PROFILE  = 1,
    PROFESSIONAL_PROFILE = 2
} EbAv1SeqProfile;


// For 8-bit and 10-bit packed inputs and outputs, the luma, cb, and cr fields should be used
//   for the three input picture planes.  However, for 10-bit unpacked planes the
//   lumaExt, cbExt, and crExt fields should be used hold the extra 2-bits of
//   precision while the luma, cb, and cr fields hold the 8-bit data.
typedef struct EbSvtIOFormat            //former EbSvtEncInput
{
    // Hosts 8 bit or 16 bit input YUV420p / YUV420p10le
    uint8_t *luma;
    uint8_t *cb;
    uint8_t *cr;

    // Hosts LSB 2 bits of 10bit input/output when the compressed 10bit format is used
    uint8_t *lumaExt;
    uint8_t *cbExt;
    uint8_t *crExt;

    uint32_t yStride;
    uint32_t crStride;
    uint32_t cbStride;

    uint32_t width;
    uint32_t height;

    uint32_t origin_x;
    uint32_t origin_y;

} EbSvtIOFormat;

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbCodec_h
