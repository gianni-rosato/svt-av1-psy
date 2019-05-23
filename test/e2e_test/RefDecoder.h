/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file RefDecoder.h
 *
 * @brief Defines a reference decoder wrapped AOM decoder
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/
#ifndef _REF_DECODER_H_
#define _REF_DECODER_H_

#include <memory.h>
#include <stdio.h>
#include "VideoFrame.h"

/** RefDecoder is a class designed for a refenece tool of conformance
 * test. It provides decoding AV1 compressed data with OBU frames, its output is
 * the YUV frame in display order. User should call get_frame right after
 * process_data to avoid missing any video frame
 */
class RefDecoder {
  public:
    /** RefDecoderErr is enumerate type of errors from decoder, refered to
     * errors in AOM */
    typedef enum {
        /*!\brief Operation completed without error */
        REF_CODEC_OK,

        /*!\brief Unspecified error */
        REF_CODEC_ERROR = 0 - AOM_CODEC_ERROR,

        /*!\brief Memory operation failed */
        REF_CODEC_MEM_ERROR = 0 - AOM_CODEC_MEM_ERROR,

        /*!\brief ABI version mismatch */
        REF_CODEC_ABI_MISMATCH = 0 - AOM_CODEC_ABI_MISMATCH,

        /*!\brief Algorithm does not have required capability */
        REF_CODEC_INCAPABLE = 0 - AOM_CODEC_INCAPABLE,

        /*!\brief The given bitstream is not supported.
         *
         * The bitstream was unable to be parsed at the highest level. The
         * decoder is unable to proceed. This error \ref SHOULD be treated as
         * fatal to the stream. */
        REF_CODEC_UNSUP_BITSTREAM = 0 - AOM_CODEC_UNSUP_BITSTREAM,

        /*!\brief Encoded bitstream uses an unsupported feature
         *
         * The decoder does not implement a feature required by the encoder.
         * This return code should only be used for features that prevent future
         * pictures from being properly decoded. This error \ref MAY be treated
         * as fatal to the stream or \ref MAY be treated as fatal to the current
         * GOP.
         */
        REF_CODEC_UNSUP_FEATURE = 0 - AOM_CODEC_UNSUP_FEATURE,

        /*!\brief The coded data for this stream is corrupt or incomplete
         *
         * There was a problem decoding the current frame.  This return code
         * should only be used for failures that prevent future pictures from
         * being properly decoded. This error \ref MAY be treated as fatal to
         * the stream or \ref MAY be treated as fatal to the current GOP. If
         * decoding is continued for the current GOP, artifacts may be present.
         */
        REF_CODEC_CORRUPT_FRAME = 0 - AOM_CODEC_CORRUPT_FRAME,

        /*!\brief An application-supplied parameter is not valid.
         *
         */
        REF_CODEC_INVALID_PARAM = 0 - AOM_CODEC_INVALID_PARAM,

        /*!\brief An iterator reached the end of list.
         *
         */
        REF_CODEC_LIST_END = 0 - AOM_CODEC_LIST_END,

        /*!\brief Decoder need more input data to generate frame
         *
         */
        REF_CODEC_NEED_MORE_INPUT = -100,

    } RefDecoderErr;

  public:
    /** Constructor of RefDecoder
     * @param ret the error code found in construction
     */
    RefDecoder(RefDecoderErr &ret);
    /** Destructor of RefDecoder	  */
    virtual ~RefDecoder();

  public:
    /** Setup decoder
     * @param init_ts  initial timestamp of input stream
     * @param interval  the time interval in two frames
     * @return
     * REF_CODEC_OK -- no error found in processing
     * others -- errors found in setup
     */
    RefDecoderErr setup(const uint64_t init_ts, const uint32_t interval);

    /** Process compressed data
     * @param data  the memory buffer of a frame of compressed data
     * @param size  the size of data
     * @return
     * REF_CODEC_OK -- no error found in processing
     * others -- errors found in process, refer to RefDecoderErr
     */
    RefDecoderErr process_data(const uint8_t *data, const uint32_t size);
    /** Get a video frame after data proceed
     * @param frame  the video frame with its attributes
     * @return
     * REF_CODEC_OK -- no error found in processing
     * others -- errors found in process, refer to RefDecoderErr
     */
    RefDecoderErr get_frame(VideoFrame &frame);

  private:
    /** Tool of translation from AOM image info to a video frame
     * @param image  the video image from AOM decoder
     * @param frame  the video frame to output
     */
    void trans_video_frame(const void *image, VideoFrame &frame);

  protected:
    void *codec_handle_;      /**< AOM codec context */
    uint32_t dec_frame_cnt_;  /**< count of decoded frame in processing */
    uint64_t init_timestamp_; /**< initial timestamp of stream */
    uint32_t frame_interval_; /**< time interval of two frame in miliseconds */
};

/** Interface of reference decoder creation
 * @return
 * RefDecoder -- decoder handle created
 * nullptr -- creation failed
 */
RefDecoder *create_reference_decoder();

#endif  // !_REF_DECODER_H_
