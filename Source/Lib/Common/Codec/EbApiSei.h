/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAPISEI_h
#define EBAPISEI_h

#include "EbDefinitions.h"

#define  MAX_DECODING_UNIT_COUNT                  64   // picture timing SEI

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    // These are temporary definitions, should be changed in the future according to user's requirements
#define MAX_CPB_COUNT    4

// User defined structures for passing data from application to the library should be added here

    typedef struct RegistedUserData
    {
        uint8_t   *user_data;    // First byte is itu_t_t35_country_code.
                              // If itu_t_t35_country_code  ==  0xFF, second byte is itu_t_t35_country_code_extension_byte.
                              // the rest are the payloadByte
        uint32_t   user_data_size;
    } RegistedUserData;

    // SEI structures
    typedef struct AppHrdParameters
    {
        EbBool                            nal_hrd_parameters_present_flag;
        EbBool                            vcl_hrd_parameters_present_flag;
        EbBool                            sub_pic_cpb_params_present_flag;

        uint32_t                          tick_divisor_minus2;
        uint32_t                          du_cpb_removal_delay_length_minus1;

        EbBool                            sub_pic_cpb_params_pic_timing_sei_flag;

        uint32_t                          dpb_output_delay_du_length_minus1;

        uint32_t                          bit_rate_scale;
        uint32_t                          cpb_size_scale;
        uint32_t                          du_cpb_size_scale;

        uint32_t                          initial_cpb_removal_delay_length_minus1;
        uint32_t                          au_cpb_removal_delay_length_minus1;
        uint32_t                          dpb_output_delay_length_minus1;

        EbBool                            fixed_pic_rate_general_flag[MAX_TEMPORAL_LAYERS];
        EbBool                            fixed_pic_rate_within_cvs_flag[MAX_TEMPORAL_LAYERS];

        uint32_t                          elemental_duration_tc_minus1[MAX_TEMPORAL_LAYERS];

        EbBool                            low_delay_hrd_flag[MAX_TEMPORAL_LAYERS];

        uint32_t                          cpb_count_minus1[MAX_TEMPORAL_LAYERS];

        uint32_t                          bit_rate_value_minus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          cpb_size_value_minus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          bit_rate_du_value_minus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          cpb_size_du_value_minus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];

        EbBool                            cbr_flag[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];

        EbBool                            cpb_dpb_delays_present_flag;
    } AppHrdParameters;

    typedef struct AppVideoUsabilityInfo
    {
        EbBool            aspect_ratio_info_present_flag;
        uint32_t          aspect_ratio_idc;
        uint32_t          sar_width;
        uint32_t          sar_height;

        EbBool            overscan_info_present_flag;
        EbBool            overscan_approriate_flag;
        EbBool            video_signal_type_present_flag;

        uint32_t          video_format;
        EbBool            video_full_range_flag;
        EbBool            color_description_present_flag;

        uint32_t          color_primaries;
        uint32_t          transfer_characteristics;
        uint32_t          matrix_coeffs;

        EbBool            chroma_loc_info_present_flag;
        uint32_t          chroma_sample_loc_type_top_field;
        uint32_t          chroma_sample_loc_type_bottom_field;

        EbBool            neutral_chroma_indication_flag;
        EbBool            field_seq_flag;
        EbBool            frame_field_info_present_flag;

        EbBool            default_display_window_flag;
        uint32_t          default_display_win_left_offset;
        uint32_t          default_display_win_right_offset;
        uint32_t          default_display_win_top_offset;
        uint32_t          default_display_win_bottom_offset;

        EbBool            vui_timing_info_present_flag;
        uint32_t          vui_num_units_in_tick;
        uint32_t          vui_time_scale;

        EbBool            vui_poc_propotional_timing_flag;
        uint32_t          vui_num_ticks_poc_diff_one_minus1;

        EbBool            vui_hrd_parameters_present_flag;

        EbBool            bitstream_restriction_flag;

        EbBool            motion_vectors_over_pic_boundaries_flag;
        EbBool            restricted_ref_pic_lists_flag;

        uint32_t          min_spatial_segmentation_idc;
        uint32_t          max_bytes_per_pic_denom;
        uint32_t          max_bits_per_min_cu_denom;
        uint32_t          log2_max_mv_length_horizontal;
        uint32_t          log2_max_mv_length_vertical;

        AppHrdParameters *hrd_parameters_ptr;
    } AppVideoUsabilityInfo;

    typedef struct AppPictureTimingSei
    {
        uint32_t   pic_struct;
        uint32_t   source_scan_type;
        EbBool     duplicate_flag;
        uint32_t   au_cpb_removal_delay_minus1;
        uint32_t   pic_dpb_output_delay;
        uint32_t   pic_dpb_output_du_delay;
        uint32_t   num_decoding_units_minus1;
        EbBool     du_common_cpb_removal_delay_flag;
        uint32_t   du_common_cpb_removal_delay_minus1;
        uint32_t   num_nalus_in_du_minus1;
        uint32_t   du_cpb_removal_delay_minus1[MAX_DECODING_UNIT_COUNT];
    } AppPictureTimingSei;

    typedef struct AppBufferingPeriodSei
    {
        uint32_t  bp_seq_parameter_set_id;
        EbBool    rap_cpb_params_present_flag;
        EbBool    concatenation_flag;
        uint32_t  au_cpb_removal_delay_delta_minus1;
        uint32_t  cpb_delay_offset;
        uint32_t  dpb_delay_offset;
        uint32_t  initial_cpb_removal_delay[2][MAX_CPB_COUNT];
        uint32_t  initial_cpb_removal_delay_offset[2][MAX_CPB_COUNT];
        uint32_t  initial_alt_cpb_removal_delay[2][MAX_CPB_COUNT];
        uint32_t  initial_alt_cpb_removal_delay_offset[2][MAX_CPB_COUNT];
    } AppBufferingPeriodSei;

    typedef struct AppRecoveryPoint
    {
        uint32_t recovery_poc_cnt;
        EbBool   exact_matching_flag;
        EbBool   broken_link_flag;
    } AppRecoveryPoint;

    // Below is an example of PanScanRectangle SEI data structure
    // Other SEI messages can have data structure in this format
    typedef struct AppPanScanRectangleSei
    {
        uint32_t pan_scan_rect_id;
        EbBool   pan_scan_rect_cancel_flag;

        uint32_t pan_scan_count_minus1;
        uint32_t pan_scan_rect_left_offset[3];
        uint32_t pan_scan_rect_right_offset[3];
        uint32_t pan_scan_rect_top_offset[3];
        uint32_t pan_scan_rect_bottom_offset[3];

        EbBool   pan_scan_rect_persist_flag;
    }AppPanScanRectangleSei;

    typedef struct EbFrameRateCfg
    {
        uint32_t num;
        uint32_t den;
    } EbFrameRateCfg;

    typedef struct EbLatencyCalc
    {
        uint64_t start_times_seconds;
        uint64_t finish_times_seconds;
        uint64_t start_timesu_seconds;
        uint64_t finish_timesu_seconds;
        uint64_t poc;
    } EbLatencyCalc;

    typedef struct EbCuQpCfg
    {
        uint32_t      qp_array_count;
        int8_t       *qp_array;

        uint32_t      cu_stride;
        uint32_t      sb_count;
        uint32_t      sb_stride;

        uint32_t      num_cus_per_last_row_lcu;
        uint32_t      num_cus_per_typical_lcu;

        // Hold the DifCUDeltaQpDepth in this struct so that it is simpler to access it in the Lib,
        // rather than having to reference something deep
    } EbCuQpCfg;

    // Signals that the default prediction structure and controls are to be
    //   overwritten and manually controlled. Manual control should be active
    //   for an entire encode, from beginning to termination.  Mixing of default
    //   prediction structure control and override prediction structure control
    //   is not supported.
    //
    // The ref_list0_count and ref_list1_count variables control the picture/slice type.
    //   I_SLICE: ref_list0_count == 0 && ref_list1_count == 0
    //   P_SLICE: ref_list0_count  > 0 && ref_list1_count == 0
    //   B_SLICE: ref_list0_count  > 0 && ref_list1_count  > 0

    typedef struct EbPredStructureCfg
    {
        uint64_t  picture_number;       // Corresponds to the display order of the picture
        uint64_t  decode_order_number;  // Corresponds to the decode order of the picture
        uint32_t  temporal_layer_index; // Corresponds to the temporal layer index of the picture
        uint32_t  nal_unit_type;        // Pictures NAL Unit Type
        uint32_t  ref_list0_count;      // A count of zero indicates the list is inactive
        uint32_t  ref_list1_count;      // A count of zero indicates the list is inactive
        EbBool    is_referenced;        // Indicates whether or not the picture is used as
                                        //   future reference.
        uint32_t  future_reference_count;
        int32_t  *future_reference_list;// Contains a list of delta POCs whose references shall
                                        //   be saved for future reference.  This signalling must
                                        //   be done with respect to decode picture order. Must
                                        //   be conformant with the DPB rules.
    } EbPredStructureCfg;

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

    // EbInputPictureDef defines the data formatting of an input picture. Note, each
    //   element can change independently of the previously coded pictures.  This allows
    //   for sub-picture coding and de-interlacing.
    typedef struct EbInputPictureDef
    {
        uint32_t width;        // Y plane width  (in units of samples)
        uint32_t height;       // Y plane height (in units of samples)

        // Padding (in lines/pixels of Y plane)
        uint32_t top_padding;
        uint32_t bot_padding;
        uint32_t left_padding;
        uint32_t right_padding;

        // Y, Cb, Cr planes. Note, for bitdepths greater than 8 using the "unpacked" format,
        //  the following elements contain the 8 MSB bits of the sample.
        EbPicturePlane y_plane;
        EbPicturePlane cb_plane;
        EbPicturePlane cr_plane;

        // Auxiliary planes used when bitdepths are greater than 8 and using the "unpacked"
        //   format.for The LSB bits of the sample should be located at the MSB bit positions
        //   of each byte. Must be set to NULL for 8-bit input or "packed" 10 or 12-bit input.
        EbPicturePlane y_aux_plane;
        EbPicturePlane cb_aux_plane;
        EbPicturePlane cr_aux_plane;
    } EbInputPictureDef;

    typedef struct EbEosDataDef
    {
        EbBool             data_valid;          // Indicates if the data attached with the last frame in the bitstream is valid or not.
                                                //   If the last frame is valid, the data will be included in the bistream
                                                //   If the last frame is NOT valid, the frame will be coded as IDR in the encoder, but
                                                //   not included in the bitstream.
    } EbEosDataDef;

    extern EbErrorType eb_app_video_usability_info_ctor(
        AppVideoUsabilityInfo *vui_ptr);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // EBAPISEI_h
