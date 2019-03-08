/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbExtFrameBuf_h
#define EbExtFrameBuf_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "stdint.h"

/*!\brief External frame buffer
 *
 * This structure holds allocated frame buffers used by the decoder.
 */
typedef struct EbExtFrameBuf {
    
   /* Pointer to the memory allocates externally for the codec
    * picture buffer */
  uint8_t       *buffer;
  
  /* Size of the memory allocates externally for the codec
   * picture buffer */  
  uint32_t       buffer_size;   

  /* Pointer to private data associated with allocating the memory */  
  void          *private_data;    

} EbExtFrameBuf;


/* Callback function to allocate frame buffer. 
 * 
 * This function is called by the decoder to allocate the
 * data for the frame buffer.
 * Parameters:
 * @  *frame_buf pointer to the frame buffer structure to be allocated
 * @  min_size  requested data size in bytes.
 * @  private_data is the private data that can be used by the allocator*/
typedef int (*eb_allocate_frame_buffer)(
            EbExtFrameBuf   *frame_buf,
            uint32_t        min_size,
            void            *private_data);

/* Callback function to release frame buffer. 
 * 
 * This function is called by the decoder to release the
 * data for the frame buffer. THe buffer should not be used
 * by the decoder. 
 *
 * Parameters:
 * @  *frame_buf pointer to the frame buffer structure to be allocated
 * @  min_size  requested data size in bytes.
 * @  private_data is the private data that can be used by the allocator*/
typedef int (*eb_release_frame_buffer)(
            EbExtFrameBuf   *fb,
            void            *private_data);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDecApi_h