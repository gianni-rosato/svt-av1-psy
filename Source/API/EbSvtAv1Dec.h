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


typedef struct EbOperatingParametersInfo {

    /*!<Specifies the time interval between the arrival of the first bit in the
     * smoothing buffer and the subsequent removal of the data that belongs to
     * the first coded frame for operating point*/
    uint32_t    decoder_buffer_delay;

    /*!<Specifies, in combination with decoder_buffer_delay[op] syntax element,
     * the first bit arrival time of frames to be decoded to the smoothing 
     * buffer */
    uint32_t    encoder_buffer_delay;

    /*!< Equal to 1 indicates that the smoothing buffer operates in low-delay 
     * mode for operating point*/
    uint8_t     low_delay_mode_flag;

} EbOperatingParametersInfo;

typedef struct EbAV1OperatingPoint {

    uint32_t    op_idc;
    uint32_t    seq_level_idx;
    uint32_t    seq_tier;

    /*!< 1 -> Indicates that there is a decoder model associated with operating
             point, 
     *   0 -> Indicates that there is not a decoder model associated with 
             operating point*/
    uint8_t     decoder_model_present_for_this_op;

    /*!< Operating Parameters Information structure*/
    EbOperatingParametersInfo operating_parameters_info;

    uint32_t    initial_display_delay_present_for_this_op;
    uint32_t    initial_display_delay;
} EbAv1OperatingPoint;

typedef struct EbColorConfig {

    /*!< bit depth */
    AomBitDepth                     bit_depth;

    /*!< 1: Indicates that the video does not contain U and V color planes.
     *   0: Indicates that the video contains Y, U, and V color planes. */
    EbBool                          mono_chrome;

    /*!< Specify the chroma subsampling format */
    uint8_t                         subsampling_x;

    /*!< Specify the chroma subsampling format */
    uint8_t                         subsampling_y;
    
    /*!< 1: Specifies that color_primaries, transfer_characteristics, and 
            matrix_coefficients are present. color_description_present_flag 
     *   0: Specifies that color_primaries, transfer_characteristics and 
            matrix_coefficients are not present */
    EbBool                         color_description_present_flag;

    /*!< An integer that is defined by the "Color primaries" section of 
     * ISO/IEC 23091-4/ITU-T H.273 */
    EbColorPrimaries                   color_primaries;

    /*!< An integer that is defined by the "Transfer characteristics" section 
     * of ISO/IEC 23091-4/ITU-T H.273 */
    EbTransferCharacteristics  transfer_characteristics;

    /*!< An integer that is defined by the "Matrix coefficients" section of 
     * ISO/IEC 23091-4/ITU-T H.273 */
    EbMatrixCoefficients       matrix_coefficients;

    /*!< 0: shall be referred to as the studio swing representation
     *   1: shall be referred to as the full swing representation */
    EbColorRange                         color_range;

    /*!< Specifies the sample position for subsampled streams */
    EbChromaSamplePosition    chroma_sample_position;

    /*!< 1: Indicates that the U and V planes may have separate delta quantizer
     *   0: Indicates that the U and V planes will share the same delta 
            quantizer value */
    EbBool                         separate_uv_delta_q;

} EbColorConfig;

typedef struct EbTimingInfo {
    /*!< Timing info present flag */
    EbBool      timing_info_present;
    
    /*!< Number of time units of a clock operating at the frequency time_scale 
     * Hz that corresponds to one increment of a clock tick counter*/
    uint32_t    num_units_in_display_tick; 

    /*!< Number of time units that pass in one second*/
    uint32_t    time_scale;

    /*!< Equal to 1 indicates that pictures should be displayed according to
     * their output order with the number of ticks between two consecutive 
     * pictures specified by num_ticks_per_picture.*/
    uint8_t     equal_picture_interval;

    /*!< Specifies the number of clock ticks corresponding to output time 
     * between two consecutive pictures in the output order. 
     * Range - [0 to (1 << 32) - 2]*/
    uint32_t    num_ticks_per_picture;

} EbTimingInfo;

typedef struct EbAV1StreamInfo
{
    /*seq_profile*/
    EbAv1SeqProfile seq_profile;

    /* Max Picture Dimensions */
    uint32_t    max_picture_width;
    uint32_t    max_picture_height;

    /* Operating points present in the bitstream */    
    uint32_t    num_operating_points;
    EbAv1OperatingPoint op_points[EB_MAX_NUM_OPERATING_POINTS];

    /* Display Timing Info*/
    EbTimingInfo    timing_info;

    /* Color description */
    EbColorConfig    color_config;

    /* Film Grain Synthesis Present */
    EbBool    film_grain_params_present;

    /* The stream is in annex_b format */
    EbBool    is_annex_b;

} EbAV1StreamInfo;

typedef struct EbAV1FrameInfo
{

    /* Layer to which the current frame belong */
    uint32_t    layer;

    /* Frame presentation time */
    uint64_t    frame_presentation_time;

} EbAV1FrameInfo;

typedef struct EbSvtAv1DecConfiguration 
{

    /* Bitstream operating point to decode. 
     * 
     * Default is -1, the highest operating point present in the bitstream
     * A value higher than the maximum number of operating points present 
     * returns the highest available operating point. */

    int32_t                operating_point;   // Operating point to decode

    /* When set to 1, returns output pictures from all scalable layers present in the bitstream. 
     * 
     * Default is 0, only one output layer is returned, defined by operating_point parameter */
    uint32_t                output_all_layers; 

    /* Skip film grain synthesis if it is present in the bitstream. Can be used for debugging purpose.  
     * 
     * Default is 0 */
    EbBool                  skip_film_grain; 

    /* Skip N output frames in the display order.
     *
     * 0 = decodes from the start of the bitstream.
     *
     * Default is 0. */
    uint64_t                 skip_frames;

    /* Maximum number of frames in the sequence to be decoded.
     *
     * 0 = decodes the full bitstream.
     *
     * Default is 0. */
    uint64_t                 frames_to_be_decoded;


    /* Offline packing of the 2bits: requires two bits packed input.
     *
     * Default is 0. */
    uint32_t                 compressed_ten_bit_format;  //remove?

    /* Outputs 8-bit pictures even if the bitstream has higher bit depth.
     * Ignored if the bitstream is 8-bit
     *
     * Default is 0. */

    EbBool                  eight_bit_output;

    /* Picture parameters */
    uint32_t                max_picture_width;
    uint32_t                max_picture_height;

    EbBitDepth              max_bit_depth;

    EbColorFormat           max_color_format;

    /* Assembly instruction set used by encoder.
    *
    * 0 = non-AVX2, C only.
    * 1 = up to AVX512, auto-select highest assembly instruction set supported.
    *
    * Default is 1. */
    uint32_t                 asm_type;
    // Application Specific parameters

    /* Number of threads used by decoder.
    *
    * 0 = System default.
    * 1 = Single thread decoding.
    *
    * Default is 0. */
    uint32_t                 threads;
    // Application Specific parameters

    /* ID assigned to each channel when multiple instances are running within the
     * same application. */
    uint32_t                 channel_id;
    uint32_t                 active_channel_count;

    uint32_t                 stat_report;

} EbSvtAv1DecConfiguration;


    /* STEP 1: Call the library to construct a Component Handle.
     *
     * Parameter:
     * @ **p_handle      Handle to be called in the future for manipulating the
     *                   component.
     * @ *p_app_data     Callback data.
     * @ *config_ptr     Pointer passed back to the client during callbacks, it will be
     *                   loaded with default parameters from the library. */
    EB_API EbErrorType eb_dec_init_handle(
        EbComponentType** p_handle,
        void* p_app_data,
        EbSvtAv1DecConfiguration  *config_ptr);

    /* STEP 2: Peek into sequence header.
     *
     * The function is used to find the sequence header
     * before the decoder is initialized. The decoder needs
     * to read the sequence header again when starting decoding.
     * The function can be called multiple times before it returns
     * the EB_ErrorNone, which means the the sequence header is found.
     * When the OBU is not a valid sequence header, EB_DecUnsupportedBitstream
     * is returned.
     *
     * Parameter:
     * @ *header            Sequence header info.
     * @ *data              Input buffer pointer.
     * @ data_size          Input data size in bytes */
    EB_API EbErrorType eb_peek_sequence_header(
        EbAV1StreamInfo *header,
        const uint8_t   *data,
        const uint32_t  data_size);

    /* STEP 3: Set configuration parameters.
     *
     * Parameter:
     * @ *svt_dec_component             Decoder handle.
     * @ *pComponentParameterStructure  Decoder and buffer configurations will be copied to the library. */
    EB_API EbErrorType eb_svt_dec_set_parameter(
        EbComponentType             *svt_dec_component,
        EbSvtAv1DecConfiguration    *pComponentParameterStructure); // pComponentParameterStructure contents will be copied to the library

    /* STEP 4: Initialize decoder and allocate memory to necessary buffers.
     *
     * Parameter:
     * @ *svt_dec_component  Decoder handle. */
    EB_API EbErrorType eb_init_decoder(
        EbComponentType         *svt_dec_component);

    /*!\brief STEP 5: Decodes one OBU.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle
     * @ *data                  Buffer with data
     * @ data_size              Data size in bytes
     * @ *user_priv             pointer to the private user data
     *
     *  Returns EB_ErrorNone if the coded data has been processed successfully. */
    EB_API EbErrorType eb_svt_decode_obu(
        EbComponentType     *svt_dec_component,
        const uint8_t       *data,
        const uint32_t       data_size);

    /*!\brief STEP 5-alt-1: Decodes a frame with associated data. The data in *data
     * should belong to one frame, possibly with sequence header and metadata.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle
     * @ *data                  Buffer with data
     * @ data_size              Data size in bytes
     *
     *  Returns EB_ErrorNone if the coded data has been processed successfully. */
    EB_API EbErrorType eb_svt_decode_frame(
        EbComponentType     *svt_dec_component,
        const uint8_t       *data,
        const uint32_t       data_size);

    /*!\brief STEP 5-alt-2: Decodes a temporal unit (TU). Decoding a TU
     * may result in several output pictures generated if output_all_layers
     * was set to 1 in the  EbSvtAv1DecConfiguration.
     * In this case, calling eb_svt_dec_get_picture() multiple times
     * would output pictures that belong to the corresponding quality layers
     * in the increasing order.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle
     * @ *data                  Buffer with data
     * @ data_size              Data size in bytes
     *
     *  Returns EB_ErrorNone if the coded data has been processed successfully. */
    EB_API EbErrorType eb_svt_decode_tu(
        EbComponentType     *svt_dec_component,
        const uint8_t       *data,
        const uint32_t       data_size);

    /* STEP 6: Get the next decoded picture. When several output pictures
     * have been generated, calling this function multiple times will
     * iterate over the decoded pictures. The previous output picture becomes
     * unavailable after the eb_svt_dec_get_picture() or one of the decoding
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
    EB_API EbErrorType eb_svt_dec_get_picture(
        EbComponentType      *svt_dec_component,
        EbBufferHeaderType   *p_buffer,
        EbAV1StreamInfo      *stream_info,
        EbAV1FrameInfo       *frame_info);

    /* STEP 7: Deinitialize decoder library.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle */
    EB_API EbErrorType eb_deinit_decoder(
        EbComponentType     *svt_dec_component);

    /* STEP 8: Deconstruct decoder handler.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle */
    EB_API EbErrorType eb_dec_deinit_handle(
        EbComponentType     *svt_dec_component);

    /*  Flush a decoder
     *
     *  Clears the decoder frame buffers. 
     *  The decoder is ready to parse a new sequence header.
     *
     *  Parameter:
     *  @ *svt_dec_component     Decoder handle
     *
     *  Returns EB_ErrorNone if the decode buffer has been flushed successfully.
     */    
    EB_API EbErrorType eb_dec_flush(
        EbComponentType     *svt_dec_component);
    
    /* Initialize callback functions.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle
     * @ allocate_buffer        callback function to allocate frame buffer
     * @ release_buffer         callback function to release frame buffer
     * @ priv_data              private data used by the allocator */
    
    EB_API EbErrorType eb_dec_set_frame_buffer_callbacks(
        EbComponentType             *svt_dec_component,
        eb_allocate_frame_buffer    allocate_buffer,
        eb_release_frame_buffer     release_buffer,
        void                        *priv_data);

    /* Returns information about the bitstream and
     * the last decoded frame.
     *
     * Parameter:
     * @ *svt_dec_component     Decoder handle.
     * @ *stream_info           Sequence header info
     * @ *frame_info            Last decoded frame info */
    EB_API EbErrorType eb_get_stream_info(
        EbComponentType             *svt_dec_component,
        EbAV1StreamInfo             *stream_info,
        EbAV1FrameInfo              *frame_info);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // EbSvtAv1Dec_h