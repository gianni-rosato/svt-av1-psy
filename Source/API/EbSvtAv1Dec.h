/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbSvtAv1Dec_h
#define EbSvtAv1Dec_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "stdint.h"
#include "EbSvtAv1.h"
#include "EbSvtAv1ExtFrameBuf.h"

typedef struct EbAV1StreamInfo {
    /*seq_profile*/
    EbAv1SeqProfile seq_profile;

    /* Max Picture Dimensions */
    uint32_t max_picture_width;
    uint32_t max_picture_height;

    /* Operating points present in the bitstream */
    uint32_t            num_operating_points;
    EbAv1OperatingPoint op_points[EB_MAX_NUM_OPERATING_POINTS];

    /* Display Timing Info*/
    EbTimingInfo timing_info;

    /* Color description */
    EbColorConfig color_config;

    /* Film Grain Synthesis Present */
    EbBool film_grain_params_present;

    /* The stream is in annex_b format */
    EbBool is_annex_b;
} EbAV1StreamInfo;

typedef struct EbAV1FrameInfo {
    /* Layer to which the current frame belong */
    uint32_t layer;

    /* Frame presentation time */
    uint64_t frame_presentation_time;
} EbAV1FrameInfo;

typedef struct EbSvtAv1DecConfiguration {
    /* Bitstream operating point to decode.
     *
     * Default is -1, the highest operating point present in the bitstream
     * A value higher than the maximum number of operating points present
     * returns the highest available operating point. */

    int32_t operating_point; // Operating point to decode

    /* When set to 1, returns output pictures from all scalable layers present in the bitstream.
     *
     * Default is 0, only one output layer is returned, defined by operating_point parameter */
    uint32_t output_all_layers;

    /* Skip film grain synthesis if it is present in the bitstream. Can be used for debugging purpose.
     *
     * Default is 0 */
    EbBool skip_film_grain;

    /* Skip N output frames in the display order.
     *
     * 0 = decodes from the start of the bitstream.
     *
     * Default is 0. */
    uint64_t skip_frames;

    /* Maximum number of frames in the sequence to be decoded.
     *
     * 0 = decodes the full bitstream.
     *
     * Default is 0. */
    uint64_t frames_to_be_decoded;

    /* Offline packing of the 2bits: requires two bits packed input.
     *
     * Default is 0. */
    uint32_t compressed_ten_bit_format; //remove?

    /* Outputs 8-bit pictures even if the bitstream has higher bit depth.
     * Ignored if the bitstream is 8-bit
     *
     * Default is 0. */

    EbBool eight_bit_output;

    /* Picture parameters */
    uint32_t max_picture_width;
    uint32_t max_picture_height;

    EbBitDepth max_bit_depth;

    EbColorFormat max_color_format;

    // Application Specific parameters

    /* Number of threads used by decoder.
     * Default is 1. */
    uint32_t threads;

    /* Number of frames that can be processed
       in parallel. Default is 1 */
    uint32_t num_p_frames;

    // Application Specific parameters

    /* ID assigned to each channel when multiple instances are running within the
     * same application. */
    uint32_t channel_id;
    uint32_t active_channel_count;

    uint32_t stat_report;
    /* Decoder internal bit-depth is set to 16-bit even if the bitstream is 8-bit
 *
 * Default is 0. */
    EbBool is_16bit_pipeline;
} EbSvtAv1DecConfiguration;

/* STEP 1: Call the library to construct a Component Handle.
     *
     * Parameter:
     * @ **p_handle      Handle to be called in the future for manipulating the
     *                   component.
     * @ *p_app_data     Callback data.
     * @ *config_ptr     Pointer passed back to the client during callbacks, it will be
     *                   loaded with default parameters from the library. */
EB_API EbErrorType svt_av1_dec_init_handle(EbComponentType **p_handle, void *p_app_data,
                                      EbSvtAv1DecConfiguration *config_ptr);

/* STEP 2: Set configuration parameters.
     *
     * Parameter:
     * @ *svt_dec_component             Decoder handle.
     * @ *pComponentParameterStructure  Decoder and buffer configurations will be copied to the library. */
EB_API EbErrorType svt_av1_dec_set_parameter(
    EbComponentType *svt_dec_component,
    EbSvtAv1DecConfiguration *
        pComponentParameterStructure); // pComponentParameterStructure contents will be copied to the library

/* STEP 3: Initialize decoder and allocate memory to necessary buffers.
     *
     * Parameter:
     * @ *svt_dec_component  Decoder handle. */
EB_API EbErrorType svt_av1_dec_init(EbComponentType *svt_dec_component);

/*!\brief STEP 4: Decodes a frame with associated data. The data in *data
     * should belong to one frame, possibly with sequence header and metadata.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle
     * @ *data                  Buffer with data
     * @ data_size              Data size in bytes
     *
     *  Returns EB_ErrorNone if the coded data has been processed successfully. */
EB_API EbErrorType svt_av1_dec_frame(EbComponentType *svt_dec_component, const uint8_t *data,
                                       const size_t data_size, uint32_t is_annexb);

/* STEP 5: Get the next decoded picture. When several output pictures
     * have been generated, calling this function multiple times will
     * iterate over the decoded pictures. The previous output picture becomes
     * unavailable after the svt_av1_dec_get_picture() or one of the decoding
     * functions is called. The pictures are returned in their display order.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle.
     * @ *p_buffer              Header pointer, picture buffer.
     * @ *stream_info           Sequence header info
     * @ *frame_info            Last decoded frame info
     *
     *  Returns EB_ErrorNone if the picture has been returned successfully.
     *  Returns EB_DecNoOutputPicture if the next output picture has not
     *  been generated yet. Calling a decoding function is needed to generate more pictures. */
EB_API EbErrorType svt_av1_dec_get_picture(EbComponentType *   svt_dec_component,
                                          EbBufferHeaderType *p_buffer,
                                          EbAV1StreamInfo *stream_info, EbAV1FrameInfo *frame_info);

/* STEP 6: Deinitialize decoder library.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle */
EB_API EbErrorType svt_av1_dec_deinit(EbComponentType *svt_dec_component);

/* STEP 7: Deconstruct decoder handler.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle */
EB_API EbErrorType svt_av1_dec_deinit_handle(EbComponentType *svt_dec_component);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // EbSvtAv1Dec_h
