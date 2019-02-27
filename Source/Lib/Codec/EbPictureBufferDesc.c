/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdlib.h>

#include "EbPictureBufferDesc.h"


/*****************************************
 * eb_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_picture_buffer_desc_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{
    EbPictureBufferDesc_t          *pictureBufferDescPtr;
    EbPictureBufferDescInitData_t  *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData_t*)object_init_data_ptr;

    uint32_t bytesPerPixel = (pictureBufferDescInitDataPtr->bit_depth == EB_8BIT) ? 1 : (pictureBufferDescInitDataPtr->bit_depth <= EB_16BIT) ? 2 : 4;

    if (pictureBufferDescInitDataPtr->bit_depth > EB_8BIT && pictureBufferDescInitDataPtr->bit_depth <= EB_16BIT && pictureBufferDescInitDataPtr->splitMode == EB_TRUE)
        bytesPerPixel = 1;

    EB_MALLOC(EbPictureBufferDesc_t*, pictureBufferDescPtr, sizeof(EbPictureBufferDesc_t), EB_N_PTR);

    // Allocate the PictureBufferDesc Object
    *object_dbl_ptr = (EbPtr)pictureBufferDescPtr;

    // Set the Picture Buffer Static variables
    pictureBufferDescPtr->maxWidth = pictureBufferDescInitDataPtr->maxWidth;
    pictureBufferDescPtr->maxHeight = pictureBufferDescInitDataPtr->maxHeight;
    pictureBufferDescPtr->width = pictureBufferDescInitDataPtr->maxWidth;
    pictureBufferDescPtr->height = pictureBufferDescInitDataPtr->maxHeight;
    pictureBufferDescPtr->bit_depth = pictureBufferDescInitDataPtr->bit_depth;
    pictureBufferDescPtr->stride_y = pictureBufferDescInitDataPtr->maxWidth + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding;
    pictureBufferDescPtr->strideCb = pictureBufferDescPtr->strideCr = pictureBufferDescPtr->stride_y >> 1;
    pictureBufferDescPtr->origin_x = pictureBufferDescInitDataPtr->left_padding;
    pictureBufferDescPtr->origin_y = pictureBufferDescInitDataPtr->top_padding;

    pictureBufferDescPtr->lumaSize = (pictureBufferDescInitDataPtr->maxWidth + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding) *
        (pictureBufferDescInitDataPtr->maxHeight + pictureBufferDescInitDataPtr->top_padding + pictureBufferDescInitDataPtr->bot_padding);
    pictureBufferDescPtr->chromaSize = pictureBufferDescPtr->lumaSize >> 2;
    pictureBufferDescPtr->packedFlag = EB_FALSE;

    if (pictureBufferDescInitDataPtr->splitMode == EB_TRUE) {
        pictureBufferDescPtr->strideBitIncY = pictureBufferDescPtr->stride_y;
        pictureBufferDescPtr->strideBitIncCb = pictureBufferDescPtr->strideCb;
        pictureBufferDescPtr->strideBitIncCr = pictureBufferDescPtr->strideCr;
    }
    else {
        pictureBufferDescPtr->strideBitIncY = 0;
        pictureBufferDescPtr->strideBitIncCb = 0;
        pictureBufferDescPtr->strideBitIncCr = 0;
    }

    // Allocate the Picture Buffers (luma & chroma)
    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->buffer_y, pictureBufferDescPtr->lumaSize      * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
        //pictureBufferDescPtr->buffer_y = (EbByte) EB_aligned_malloc( pictureBufferDescPtr->lumaSize      * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        pictureBufferDescPtr->bufferBitIncY = 0;
        if (pictureBufferDescInitDataPtr->splitMode == EB_TRUE) {
            EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferBitIncY, pictureBufferDescPtr->lumaSize      * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
            //pictureBufferDescPtr->bufferBitIncY = (EbByte) EB_aligned_malloc( pictureBufferDescPtr->lumaSize      * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        }

        memset(pictureBufferDescPtr->buffer_y, 0, pictureBufferDescPtr->lumaSize * bytesPerPixel * sizeof(uint8_t));

    }
    else {
        pictureBufferDescPtr->buffer_y = 0;
        pictureBufferDescPtr->bufferBitIncY = 0;
    }

    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferCb, pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
        //pictureBufferDescPtr->bufferCb = (EbByte) EB_aligned_malloc(pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        pictureBufferDescPtr->bufferBitIncCb = 0;
        if (pictureBufferDescInitDataPtr->splitMode == EB_TRUE) {
            EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferBitIncCb, pictureBufferDescPtr->chromaSize      * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
            //pictureBufferDescPtr->bufferBitIncCb = (EbByte) EB_aligned_malloc(pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        }

        memset(pictureBufferDescPtr->bufferCb, 0, pictureBufferDescPtr->chromaSize * bytesPerPixel * sizeof(uint8_t));

    }
    else {
        pictureBufferDescPtr->bufferCb = 0;
        pictureBufferDescPtr->bufferBitIncCb = 0;
    }

    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferCr, pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
        //pictureBufferDescPtr->bufferCr = (EbByte) EB_aligned_malloc(pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        pictureBufferDescPtr->bufferBitIncCr = 0;
        if (pictureBufferDescInitDataPtr->splitMode == EB_TRUE) {
            EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferBitIncCr, pictureBufferDescPtr->chromaSize      * bytesPerPixel * sizeof(uint8_t), EB_A_PTR);
            //pictureBufferDescPtr->bufferBitIncCr = (EbByte) EB_aligned_malloc(pictureBufferDescPtr->chromaSize    * bytesPerPixel * sizeof(uint8_t),ALVALUE);
        }

        memset(pictureBufferDescPtr->bufferCr, 0, pictureBufferDescPtr->chromaSize * bytesPerPixel * sizeof(uint8_t));

    }
    else {
        pictureBufferDescPtr->bufferCr = 0;
        pictureBufferDescPtr->bufferBitIncCr = 0;
    }

    return EB_ErrorNone;
}


/*****************************************
 * eb_recon_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_recon_picture_buffer_desc_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{
    EbPictureBufferDesc_t          *pictureBufferDescPtr;
    EbPictureBufferDescInitData_t  *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData_t*)object_init_data_ptr;

    uint32_t bytesPerPixel = (pictureBufferDescInitDataPtr->bit_depth == EB_8BIT) ? 1 : 2;

    EB_MALLOC(EbPictureBufferDesc_t*, pictureBufferDescPtr, sizeof(EbPictureBufferDesc_t), EB_N_PTR);

    // Allocate the PictureBufferDesc Object
    *object_dbl_ptr = (EbPtr)pictureBufferDescPtr;

    // Set the Picture Buffer Static variables
    pictureBufferDescPtr->maxWidth = pictureBufferDescInitDataPtr->maxWidth;
    pictureBufferDescPtr->maxHeight = pictureBufferDescInitDataPtr->maxHeight;
    pictureBufferDescPtr->width = pictureBufferDescInitDataPtr->maxWidth;
    pictureBufferDescPtr->height = pictureBufferDescInitDataPtr->maxHeight;
    pictureBufferDescPtr->bit_depth = pictureBufferDescInitDataPtr->bit_depth;
    pictureBufferDescPtr->stride_y = pictureBufferDescInitDataPtr->maxWidth + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding;
    pictureBufferDescPtr->strideCb = pictureBufferDescPtr->strideCr = pictureBufferDescPtr->stride_y >> 1;
    pictureBufferDescPtr->origin_x = pictureBufferDescInitDataPtr->left_padding;
    pictureBufferDescPtr->origin_y = pictureBufferDescInitDataPtr->top_padding;

    pictureBufferDescPtr->lumaSize = (pictureBufferDescInitDataPtr->maxWidth + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding) *
        (pictureBufferDescInitDataPtr->maxHeight + pictureBufferDescInitDataPtr->top_padding + pictureBufferDescInitDataPtr->bot_padding);
    pictureBufferDescPtr->chromaSize = pictureBufferDescPtr->lumaSize >> 2;
    pictureBufferDescPtr->packedFlag = EB_FALSE;

    pictureBufferDescPtr->strideBitIncY = 0;
    pictureBufferDescPtr->strideBitIncCb = 0;
    pictureBufferDescPtr->strideBitIncCr = 0;

    // Allocate the Picture Buffers (luma & chroma)
    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Y_FLAG) {

        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->buffer_y, pictureBufferDescPtr->lumaSize * bytesPerPixel, EB_A_PTR);
        memset(pictureBufferDescPtr->buffer_y, 0, pictureBufferDescPtr->lumaSize      * bytesPerPixel);

    }
    else {
        pictureBufferDescPtr->buffer_y = 0;
    }

    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Cb_FLAG) {

        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferCb, pictureBufferDescPtr->chromaSize * bytesPerPixel, EB_A_PTR);
        memset(pictureBufferDescPtr->bufferCb, 0, pictureBufferDescPtr->chromaSize      * bytesPerPixel);

    }
    else {
        pictureBufferDescPtr->bufferCb = 0;
    }

    if (pictureBufferDescInitDataPtr->bufferEnableMask & PICTURE_BUFFER_DESC_Cr_FLAG) {

        EB_ALLIGN_MALLOC(EbByte, pictureBufferDescPtr->bufferCr, pictureBufferDescPtr->chromaSize * bytesPerPixel, EB_A_PTR);
        memset(pictureBufferDescPtr->bufferCr, 0, pictureBufferDescPtr->chromaSize      * bytesPerPixel);

    }
    else {
        pictureBufferDescPtr->bufferCr = 0;
    }

    return EB_ErrorNone;
}


void LinkEbToAomBufferDesc(
    EbPictureBufferDesc_t          *picBuffDsc,
    Yv12BufferConfig             *aomBuffDsc
)
{

    //NOTe:  Not all fileds are connected. add more connections as needed.
    if (picBuffDsc->bit_depth == EB_8BIT) {

        aomBuffDsc->y_buffer = picBuffDsc->buffer_y + picBuffDsc->origin_x + (picBuffDsc->origin_y     * picBuffDsc->stride_y);
        aomBuffDsc->u_buffer = picBuffDsc->bufferCb + picBuffDsc->origin_x / 2 + (picBuffDsc->origin_y / 2 * picBuffDsc->strideCb);
        aomBuffDsc->v_buffer = picBuffDsc->bufferCr + picBuffDsc->origin_x / 2 + (picBuffDsc->origin_y / 2 * picBuffDsc->strideCb);


        aomBuffDsc->y_width = picBuffDsc->width;
        aomBuffDsc->uv_width = picBuffDsc->width / 2;

        aomBuffDsc->y_height = picBuffDsc->height;
        aomBuffDsc->uv_height = picBuffDsc->height / 2;


        aomBuffDsc->y_stride = picBuffDsc->stride_y;
        aomBuffDsc->uv_stride = picBuffDsc->strideCb;

        aomBuffDsc->border = picBuffDsc->origin_x;

        aomBuffDsc->subsampling_x = 1;
        aomBuffDsc->subsampling_y = 1;

        aomBuffDsc->y_crop_width = aomBuffDsc->y_width;
        aomBuffDsc->uv_crop_width = aomBuffDsc->uv_width;
        aomBuffDsc->y_crop_height = aomBuffDsc->y_height;
        aomBuffDsc->uv_crop_height = aomBuffDsc->uv_height;


        aomBuffDsc->flags = 0;
    }
    else {

        /*
        Moving within a 16bit memory area: 2 possible mecanisms:

        1. to move from one location to another by an offset x, using 16bit pointers
           int32_t x;
           U16* Base16b;
           U16* NewAdd16b = Base16b + x
           int32_t data = NewAdd16b[0];

         2. to move from one location to another by an offset x, using 8bit pointers

            int32_t x;
            U16* Base16b;

            U16* baseAd8b = Base16b/2; //convert the base address into 8bit
            U16* newAd8b  = baseAd8b + x;

            then before reading the data, we need to convert the pointer back to 16b
            U16* NewAdd16b = newAd8b*2 ;
            int32_t data = NewAdd16b[0];

            NewAdd16b = Base16b + off
                      = Base16b_asInt + 2*off
                      =(Base16b_asInt/2 +off)*2
        */

        aomBuffDsc->y_buffer = CONVERT_TO_BYTEPTR(picBuffDsc->buffer_y);
        aomBuffDsc->u_buffer = CONVERT_TO_BYTEPTR(picBuffDsc->bufferCb);
        aomBuffDsc->v_buffer = CONVERT_TO_BYTEPTR(picBuffDsc->bufferCr);

        aomBuffDsc->y_buffer += picBuffDsc->origin_x + (picBuffDsc->origin_y     * picBuffDsc->stride_y);
        aomBuffDsc->u_buffer += picBuffDsc->origin_x / 2 + (picBuffDsc->origin_y / 2 * picBuffDsc->strideCb);
        aomBuffDsc->v_buffer += picBuffDsc->origin_x / 2 + (picBuffDsc->origin_y / 2 * picBuffDsc->strideCb);


        aomBuffDsc->y_width = picBuffDsc->width;
        aomBuffDsc->uv_width = picBuffDsc->width / 2;

        aomBuffDsc->y_height = picBuffDsc->height;
        aomBuffDsc->uv_height = picBuffDsc->height / 2;


        aomBuffDsc->y_stride = picBuffDsc->stride_y;
        aomBuffDsc->uv_stride = picBuffDsc->strideCb;

        aomBuffDsc->border = picBuffDsc->origin_x;

        aomBuffDsc->subsampling_x = 1;
        aomBuffDsc->subsampling_y = 1;

        aomBuffDsc->y_crop_width = aomBuffDsc->y_width;
        aomBuffDsc->uv_crop_width = aomBuffDsc->uv_width;
        aomBuffDsc->y_crop_height = aomBuffDsc->y_height;
        aomBuffDsc->uv_crop_height = aomBuffDsc->uv_height;
        aomBuffDsc->flags = YV12_FLAG_HIGHBITDEPTH;

    }

}



void *aom_memalign(size_t align, size_t size);
void aom_free(void *memblk);

#define yv12_align_addr(addr, align) \
  (void *)(((size_t)(addr) + ((align)-1)) & (size_t) - (align))



int32_t aom_realloc_frame_buffer(Yv12BufferConfig *ybf, int32_t width, int32_t height,
    int32_t ss_x, int32_t ss_y, int32_t use_highbitdepth,
    int32_t border, int32_t byte_alignment,
    aom_codec_frame_buffer_t *fb,
    aom_get_frame_buffer_cb_fn_t cb, void *cb_priv) {
    if (ybf) {
        const int32_t aom_byte_align = (byte_alignment == 0) ? 1 : byte_alignment;
        const int32_t aligned_width = (width + 7) & ~7;
        const int32_t aligned_height = (height + 7) & ~7;
        const int32_t y_stride = ((aligned_width + 2 * border) + 31) & ~31;
        const uint64_t yplane_size =
            (aligned_height + 2 * border) * (uint64_t)y_stride + byte_alignment;
        const int32_t uv_width = aligned_width >> ss_x;
        const int32_t uv_height = aligned_height >> ss_y;
        const int32_t uv_stride = y_stride >> ss_x;
        const int32_t uv_border_w = border >> ss_x;
        const int32_t uv_border_h = border >> ss_y;
        const uint64_t uvplane_size =
            (uv_height + 2 * uv_border_h) * (uint64_t)uv_stride + byte_alignment;

        const uint64_t frame_size =
            (1 + use_highbitdepth) * (yplane_size + 2 * uvplane_size);

        uint8_t *buf = NULL;

        if (cb != NULL) {
            const int32_t align_addr_extra_size = 31;
            const uint64_t external_frame_size = frame_size + align_addr_extra_size;

            assert(fb != NULL);

            if (external_frame_size != (size_t)external_frame_size) return -1;

            // Allocation to hold larger frame, or first allocation.
            if (cb(cb_priv, (size_t)external_frame_size, fb) < 0) return -1;

            if (fb->data == NULL || fb->size < external_frame_size) return -1;

            ybf->buffer_alloc = (uint8_t *)yv12_align_addr(fb->data, 32);

#if defined(__has_feature)
#if __has_feature(memory_sanitizer)
            // This memset is needed for fixing the issue of using uninitialized
            // value in msan test. It will cause a perf loss, so only do this for
            // msan test.
            memset(ybf->buffer_alloc, 0, (int32_t)frame_size);
#endif
#endif
        }
        else if (frame_size > (size_t)ybf->buffer_alloc_sz) {
            // Allocation to hold larger frame, or first allocation.
            aom_free(ybf->buffer_alloc);
            ybf->buffer_alloc = NULL;

            if (frame_size != (size_t)frame_size) return -1;


            EB_MALLOC(uint8_t *, ybf->buffer_alloc, frame_size, EB_N_PTR);

            if (!ybf->buffer_alloc) return -1;

            ybf->buffer_alloc_sz = (size_t)frame_size;

            // This memset is needed for fixing valgrind error from C loop filter
            // due to access uninitialized memory in frame border. It could be
            // removed if border is totally removed.
            memset(ybf->buffer_alloc, 0, ybf->buffer_alloc_sz);
        }

        /* Only support allocating buffers that have a border that's a multiple
        * of 32. The border restriction is required to get 16-byte alignment of
        * the start of the chroma rows without introducing an arbitrary gap
        * between planes, which would break the semantics of things like
        * aom_img_set_rect(). */
        if (border & 0x1f) return -3;

        ybf->y_crop_width = width;
        ybf->y_crop_height = height;
        ybf->y_width = aligned_width;
        ybf->y_height = aligned_height;
        ybf->y_stride = y_stride;

        ybf->uv_crop_width = (width + ss_x) >> ss_x;
        ybf->uv_crop_height = (height + ss_y) >> ss_y;
        ybf->uv_width = uv_width;
        ybf->uv_height = uv_height;
        ybf->uv_stride = uv_stride;

        ybf->border = border;
        ybf->frame_size = (size_t)frame_size;
        ybf->subsampling_x = ss_x;
        ybf->subsampling_y = ss_y;

        buf = ybf->buffer_alloc;
        if (use_highbitdepth) {
            // Store uint16 addresses when using 16bit framebuffers
            buf = CONVERT_TO_BYTEPTR(ybf->buffer_alloc);
            ybf->flags = YV12_FLAG_HIGHBITDEPTH;
        }
        else {
            ybf->flags = 0;
        }

        ybf->y_buffer = (uint8_t *)yv12_align_addr(
            buf + (border * y_stride) + border, aom_byte_align);
        ybf->u_buffer = (uint8_t *)yv12_align_addr(
            buf + yplane_size + (uv_border_h * uv_stride) + uv_border_w,
            aom_byte_align);
        ybf->v_buffer =
            (uint8_t *)yv12_align_addr(buf + yplane_size + uvplane_size +
            (uv_border_h * uv_stride) + uv_border_w,
                aom_byte_align);

        ybf->use_external_refernce_buffers = 0;

        //if (use_highbitdepth) {
        //    if (ybf->y_buffer_8bit) aom_free(ybf->y_buffer_8bit);
        //    ybf->y_buffer_8bit = (uint8_t *)aom_memalign(32, (size_t)yplane_size);
        //    if (!ybf->y_buffer_8bit) return -1;
        //}
        //else {
        //    assert(!ybf->y_buffer_8bit);
        //}

        ybf->corrupted = 0; /* assume not corrupted by errors */
        return 0;
    }
    return -2;
}
