/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAV1Structs_h
#define EbAV1Structs_h

#ifdef __cplusplus
extern "C" {
#endif

/*!\brief OBU types. */
typedef enum ATTRIBUTE_PACKED {
    OBU_SEQUENCE_HEADER = 1,
    OBU_TEMPORAL_DELIMITER = 2,
    OBU_FRAME_HEADER = 3,
    OBU_TILE_GROUP = 4,
    OBU_METADATA = 5,
    OBU_FRAME = 6,
    OBU_REDUNDANT_FRAME_HEADER = 7,
    OBU_PADDING = 15,
} obuType;

typedef struct ObuHeader {

   /*!<Size (1 or 2 B) of OBU header (including optional OBU extension header) */
    size_t size;

    /*!<Must be set to 0*/
    uint8_t obu_forbidden_bit;

    /*!<Specifies the type of data structure contained in the OBU payload*/
    obuType obu_type;

    /*!<Indicates if the optional obu_extension_header is present*/
    uint8_t obu_extension_flag;

    /*!< 1: indicates that the obu_size syntax element will be present
     *   0: Indicates that the obu_size syntax element will not be present*/
    uint8_t obu_has_size_field;

    /*!<Specifies the temporal level of the data contained in the OBU*/
    uint8_t temporal_id;

    /*!<Specifies the spatial level of the data contained in the OBU*/
    uint8_t spatial_id;

    size_t payload_size;

} ObuHeader;


typedef struct DecoderModelInfo {

    /*!< Specifies the length of the decoder_buffer_delay and the
     * encoder_buffer_delay syntax elements, in bits*/
    uint8_t     buffer_delay_length_minus_1;

    /*!< Number of time units of a decoding clock operating at the frequency
     *time_scale Hz that corresponds to one increment of a clock tick counter*/
    uint32_t    num_units_in_decoding_tick;

    /*!<Specifies the length of the buffer_removal_time syntax element,in bits*/
    uint8_t     buffer_removal_time_length_minus_1;

    /*!< Specifies the length of the frame_presentation_time syntax element,
     * in bits*/
    uint8_t     frame_presentation_time_length_minus_1;


} DecoderModelInfo;


typedef struct OrderHintInfo {

    /*!<1: Indicates that tools based on the values of order hints may be used.
     *  0: Indicates that tools based on order hints are disabled */
    uint8_t             enable_order_hint;

    /*!< 1: Indicates that the distance weights process may be used for
     * inter prediction */
    uint8_t             enable_jnt_comp;

    /*!< 1: Indicates that the use_ref_frame_mvs syntax element may be present.
     *   0: Indicates that the use_ref_frame_mvs syntax element will not be
            present */
    uint8_t             enable_ref_frame_mvs;

    /*!< Used to compute OrderHintBits*/
    uint8_t             order_hint_bits;

} OrderHintInfo;


typedef struct SeqHeader {

    /*!< Specifies the features that can be used in the coded video sequence */
    EbAv1SeqProfile   seq_profile;

    /*!<1: Specifies that the coded video sequence contains only one coded frame
     *  0: Specifies that the coded video sequence contains one or more coded
           frames */
    uint8_t             still_picture;

    /*!< Specifies that the syntax elements not needed by a still picture are
     * omitted */
    uint8_t             reduced_still_picture_header;

    /*!< Timing Information structure*/
    EbTimingInfo        timing_info;

    /*!< Specifies whether decoder model information is present in the coded
     * video sequence */
    uint8_t             decoder_model_info_present_flag;

    /*!< Decoder Mode Information structure*/
    DecoderModelInfo  decoder_model_info;

    /*!< Specifies whether initial display delay information is present in the
     * coded video sequence.*/
    uint8_t             initial_display_delay_present_flag;

    /*!< Indicates the number of operating points minus 1 present in the coded video sequence*/
    uint8_t             operating_points_cnt_minus_1;

    /*!< Operating Point Param structure*/
    EbAv1OperatingPoint operating_point[MAX_NUM_OPERATING_POINTS];

    /*!< Specifies the number of bits minus 1 used for transmitting the frame
     * width syntax elements */
    uint8_t             frame_width_bits;

    /*!< Specifies the number of bits minus 1 used for transmitting the frame
     * height syntax elements*/
    uint8_t             frame_height_bits;

    /*!< Specifies the maximum frame width minus 1 for the frames represented
     * by this sequence header */
    uint16_t            max_frame_width;

    /*!< Specifies the maximum frame height minus 1 for the frames represented
     * by this sequence header */
    uint16_t            max_frame_height;

    /*!< Specifies whether frame id numbers are present in the coded video
     * sequence*/
    uint8_t             frame_id_numbers_present_flag;

    /*!< Specifies the number of bits used to encode delta_frame_id
     * syntax elements*/
    uint8_t             delta_frame_id_length;

    /*!<Used to calculate the number of bits used to encode the frame_id syntax
     * element.*/
    uint8_t             frame_id_length;

    /*!< 1: Indicates that superblocks contain 128x128 luma samples
     *   0: Indicates that superblocks contain 64x64 luma samples.*/
    uint8_t             use_128x128_superblock;

    /*Size of the superblock used for this frame*/
    BlockSize         sb_size;

    /*!< superblock size in 4x4 MI unit */
    uint8_t             sb_mi_size;

    /*!< superblock size inlog2 unit */
    uint8_t             sb_size_log2;

    /*!< 1: Specifies that the use_filter_intra syntax element may be present.
     *   0: Specifies that the use_filter_intra syntax element will not be
     *       present*/
    uint8_t             enable_filter_intra;

    /*!< Specifies whether the intra edge filtering process should be enabled */
    uint8_t             enable_intra_edge_filter;

    /*!<1: Specifies that the mode info for inter blocks may contain the syntax
     *     element interintra.
     *  0: Specifies that the syntax element interintra will not be present */
    uint8_t             enable_interintra_compound;

    /*!<1: Specifies that the mode info for inter blocks may contain the syntax
     *     element compound_type
     *  0: Specifies that the syntax element compound_type will not be present*/
    uint8_t             enable_masked_compound;

    /*!<1: Indicates that the allow_warped_motion syntax element may be present
     *  0: Indicates that the allow_warped_motion syntax element will not be
     *     present*/
    uint8_t             enable_warped_motion;

    /*!< 1: Indicates that the inter prediction filter type may be specified
     *      independently in the horizontal and vertical directions.
     *   0: Indicates only one filter type may be specified, which is then used
     *      in both directions.*/
    uint8_t             enable_dual_filter;

    /*!< Order Hint Information structure */
    OrderHintInfo     order_hint_info;

    /*!<Equal to SELECT_SCREEN_CONTENT_TOOLS, indicates that the
     * allow_screen_content_tools syntax element will be present in the frame
     * header. Otherwise, seq_force_screen_content_tools contains the value for
     * allow_screen_content_tools*/
    uint8_t             seq_force_screen_content_tools;

    /*!< Equal to SELECT_INTEGER_MV indicates that the force_integer_mv syntax
     * element will be present in the frame header (providing
     * allow_screen_content_tools is equal to 1). Otherwise, seq_force_integer_mv
     * contains the value for force_integer_mv */
    uint8_t             seq_force_integer_mv;

    /*!< 1: Specifies that the use_superres syntax element will be present in
     *      the uncompressed header.
     *   0: Specifies that the use_superres syntax element will not be present*/
    uint8_t             enable_superres;

    /*!< 1: Specifies that cdef filtering may be enabled.
         0: specifies that cdef filtering is disabled */
    uint8_t             enable_cdef;

    /*!< 1: Specifies that loop restoration filtering may be enabled.
         0: Specifies that loop restoration filtering is disabled*/
    uint8_t             enable_restoration;

    /*!< Colour Configuration structure*/
    EbColorConfig       color_config;

    /*!< Specifies whether film grain parameters are present in the coded video sequence*/
    uint8_t             film_grain_params_present;

} SeqHeader;


#ifdef __cplusplus
}
#endif
#endif // EbAV1Structs_h
