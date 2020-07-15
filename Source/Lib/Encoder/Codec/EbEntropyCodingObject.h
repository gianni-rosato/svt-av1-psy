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
    OutputBitstreamUnit* output_bitstream_ptr;
} Bitstream;

typedef struct EntropyCoder {
    EbDctor        dctor;
#if !EC_MEM_OPT
    EbPtr          cabac_encode_context_ptr;
#endif
    FRAME_CONTEXT *fc; /* this frame entropy */
    AomWriter      ec_writer;
    EbPtr          ec_output_bitstream_ptr;
    uint64_t       ec_frame_size;
} EntropyCoder;

typedef struct EntropyTileInfo
{
    EbDctor           dctor;
    EntropyCoder     *entropy_coder_ptr;
    int8_t            entropy_coding_current_available_row;
    EbBool            entropy_coding_row_array[MAX_SB_ROWS];
    int8_t            entropy_coding_current_row;
    int8_t            entropy_coding_row_count;
    EbHandle          entropy_coding_mutex;
    EbBool            entropy_coding_in_progress;
    EbBool            entropy_coding_tile_done;
} EntropyTileInfo;

extern EbErrorType entropy_tile_info_ctor(
        EntropyTileInfo *entropy_tile_info_ptr,
        uint32_t buf_size);

extern EbErrorType bitstream_ctor(Bitstream *bitstream_ptr, uint32_t buffer_size);

void bitstream_reset(Bitstream* bitstream_ptr);

int bitstream_get_bytes_count(const Bitstream* bitstream_ptr);

//copy size bytes from bistream_ptr to dst
void bitstream_copy(const Bitstream* bitstream_ptr, void* dest, int size);

extern EbErrorType entropy_coder_ctor(EntropyCoder *entropy_coder_ptr, uint32_t buffer_size);

#if !EC_MEM_OPT
extern OutputBitstreamUnit* entropy_coder_get_bitstream_ptr(EntropyCoder *entropy_coder_ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbEntropyCodingObject_h
