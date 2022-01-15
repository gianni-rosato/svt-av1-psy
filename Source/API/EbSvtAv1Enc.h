/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbSvtAv1Enc_h
#define EbSvtAv1Enc_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include <stdint.h>
#include "EbSvtAv1.h"
#include <stdlib.h>
#include <stdio.h>

//***HME***

#define MAX_HIERARCHICAL_LEVEL 6
#define REF_LIST_MAX_DEPTH 4
#define MAX_ENC_PRESET 13
#define NUM_MV_COMPONENTS 2
#define NUM_MV_HIST 2
#define MAX_MV_HIST_SIZE 2 * REF_LIST_MAX_DEPTH *NUM_MV_COMPONENTS *NUM_MV_HIST
#define DEFAULT -1

#define EB_BUFFERFLAG_EOS 0x00000001 // signals the last packet of the stream
#define EB_BUFFERFLAG_SHOW_EXT \
    0x00000002 // signals that the packet contains a show existing frame at the end
#define EB_BUFFERFLAG_HAS_TD 0x00000004 // signals that the packet contains a TD
#define EB_BUFFERFLAG_IS_ALT_REF 0x00000008 // signals that the packet contains an ALT_REF frame
#define EB_BUFFERFLAG_ERROR_MASK \
    0xFFFFFFF0 // mask for signalling error assuming top flags fit in 4 bits. To be changed, if more flags are added.

/*
 * Struct for storing content light level information
 * Values are stored in BE format
 * Refer to the AV1 specification 6.7.3 for more details
 */
struct EbContentLightLevel {
    uint16_t max_cll;
    uint16_t max_fall;
};

/*
 * Struct for storing x and y chroma points, values are stored in BE format
 */
struct EbSvtAv1ChromaPoints {
    uint16_t x;
    uint16_t y;
};

/*
 * Struct for storing mastering-display information
 * values are stored in BE format
 * Refer to the AV1 specification 6.7.4 for more details
 */
struct EbSvtAv1MasteringDisplayInfo {
    struct EbSvtAv1ChromaPoints r;
    struct EbSvtAv1ChromaPoints g;
    struct EbSvtAv1ChromaPoints b;
    struct EbSvtAv1ChromaPoints white_point;
    uint32_t                    max_luma;
    uint32_t                    min_luma;
};

/************************************************
 * Prediction Structure Config Entry
 *   Contains the basic reference lists and
 *   configurations for each Prediction Structure
 *   Config Entry.
 ************************************************/
typedef struct PredictionStructureConfigEntry {
    uint32_t temporal_layer_index;
    uint32_t decode_order;
    int32_t  ref_list0[REF_LIST_MAX_DEPTH];
    int32_t  ref_list1[REF_LIST_MAX_DEPTH];
} PredictionStructureConfigEntry;

// super-res modes
typedef enum {
    SUPERRES_NONE, // No frame superres allowed.
    SUPERRES_FIXED, // All frames are coded at the specified scale, and super-resolved.
    SUPERRES_RANDOM, // All frames are coded at a random scale, and super-resolved.
    SUPERRES_QTHRESH, // Superres scale for a frame is determined based on q_index.
    SUPERRES_AUTO, // Automatically select superres for appropriate frames.
    SUPERRES_MODES
} SUPERRES_MODE;

// super-res auto search type
typedef enum {
    SUPERRES_AUTO_ALL, // Tries all possible superres ratios
    SUPERRES_AUTO_DUAL, // Tries no superres and q-based superres ratios
    SUPERRES_AUTO_SOLO, // Only apply the q-based superres ratio
    SUPERRES_AUTO_SEARCH_TYPES
} SUPERRES_AUTO_SEARCH_TYPE;

/** The SvtAv1IntraRefreshType is used to describe the intra refresh type.
*/
typedef enum SvtAv1IntraRefreshType {
    SVT_AV1_FWDKF_REFRESH = 1,
    SVT_AV1_KF_REFRESH = 2,
} SvtAv1IntraRefreshType;

typedef enum {
    SVT_AV1_STREAM_INFO_START                = 1,
    SVT_AV1_STREAM_INFO_FIRST_PASS_STATS_OUT = SVT_AV1_STREAM_INFO_START,

    SVT_AV1_STREAM_INFO_END,
} SVT_AV1_STREAM_INFO_ID;

/*!\brief Generic fixed size buffer structure
 *
 * This structure is able to hold a reference to any fixed size buffer.
 */
typedef struct SvtAv1FixedBuf {
    void *   buf; /**< Pointer to the data. Does NOT own the data! */
    uint64_t sz; /**< Length of the buffer, in chars */
} SvtAv1FixedBuf; /**< alias for struct aom_fixed_buf */


// Will contain the EbEncApi which will live in the EncHandle class
// Only modifiable during config-time.
typedef struct EbSvtAv1EncConfiguration {
    // Encoding preset

    /**
     * @brief Encoder preset used.
     * -2 and -1 are for debug purposes and should not be used.
     * 0 is the highest quality mode but is the slowest,
     * 13 is the fastest mode but is not as high quality.
     *
     * Min value is -2.
     * Max value is 13.
     * Default is 12.
     */
    int8_t enc_mode;

    // GOP Structure

    /* The intra period defines the interval of frames after which you insert an
     * Intra refresh. It is strongly recommended to set the value to multiple of
     * 8 minus 1 the closest to 1 second (e.g. 55, 47, 31, 23 should be used for
     * 60, 50, 30, (24 or 25) respectively.
     *
     * -1 = no intra update.
     * -2 = auto.
     *
     * Default is -2. */
    int32_t intra_period_length;
    /* Random access.
     *
     * 1 = CRA, open GOP.
     * 2 = IDR, closed GOP.
     *
     * Default is 1. */
    SvtAv1IntraRefreshType intra_refresh_type;
    /* Number of hierarchical layers used to construct GOP.
     * Minigop size = 2^HierarchicalLevels.
     *
     * Default is 3. */
    uint32_t hierarchical_levels;

    /* Prediction structure used to construct GOP. There are two main structures
     * supported, which are: Low Delay (P or B) and Random Access.
     *
     * In Low Delay structure, pictures within a mini GOP refer to the previously
     * encoded pictures in display order. In other words, pictures with display
     * order N can only be referenced by pictures with display order greater than
     * N, and it can only refer pictures with picture order lower than N. The Low
     * Delay structure can be flat structured (e.g. IPPPPPPP...) or hierarchically
     * structured. B/b pictures can be used instead of P/p pictures. However, the
     * reference picture list 0 and the reference picture list 1 will contain the
     * same reference picture.
     *
     * In Random Access structure, the B/b pictures can refer to reference pictures
     * from both directions (past and future).
     *
     * Default is 2. */
    uint8_t pred_structure;

    // Input Info

    /**
     * @brief Frame width in pixels.
     *
     * Min is 64.
     * Max is 16384.
     * Default is 0.
     */
    uint32_t source_width;

    /**
     * @brief Frame height in pixels
     *
     * Min is 64.
     * Max is 8704.
     * Default is 0.
     */
    uint32_t source_height;

    /* The frequecy of images being displayed. If the number is less than 1000,
     * the input frame rate is an integer number between 1 and 60, else the input
     * number is in Q16 format, shifted by 16 bits, where max allowed is 240 fps.
     * If FrameRateNumerator and FrameRateDenominator are both not equal to zero,
     * the encoder will ignore this parameter.
     *
     * Default is 25. */
    uint32_t frame_rate;

    /* Frame rate numerator. When zero, the encoder will use -fps if
     * FrameRateDenominator is also zero, otherwise an error is returned.
     *
     * Default is 0. */
    uint32_t frame_rate_numerator;
    /* Frame rate denominator. When zero, the encoder will use -fps if
     * FrameRateNumerator is also zero, otherwise an error is returned.
     *
     * Default is 0. */
    uint32_t frame_rate_denominator;
    /* Specifies the bit depth of input video.
     *
     * 8 = 8 bit.
     * 10 = 10 bit.
     *
     * Default is 8. */
    uint32_t encoder_bit_depth;

    /**
     * @brief Encoder color format.
     * Only YUV420 is supported for now.
     *
     * Min is YUV400.
     * Max is YUV444.
     * Default is YUV420.
     */
    EbColorFormat encoder_color_format;
    /* Offline packing of the 2bits: requires two bits packed input.
     *
     * Default is 0. */
    uint32_t compressed_ten_bit_format;

    /* Instruct the library to calculate the recon to source for PSNR calculation
    *
    * Default is 0.*/
    uint32_t stat_report;

    // Quantization
    /* Initial quantization parameter for the Intra pictures used under constant
     * qp rate control mode.
     *
     * Default is 50. */
    uint32_t qp;

    /* force qp values for every picture that are passed in the header pointer
    *
    * Default is 0.*/
    EbBool use_qp_file;

    /* use fixed qp offset for every picture based on temporal layer index
    *
    * Default is 0.*/
    EbBool  use_fixed_qindex_offsets;
    int32_t qindex_offsets[EB_MAX_TEMPORAL_LAYERS];
    int32_t key_frame_chroma_qindex_offset;
    int32_t key_frame_qindex_offset;
    int32_t chroma_qindex_offsets[EB_MAX_TEMPORAL_LAYERS];

    // input / output buffer to be used for multi-pass encoding
    SvtAv1FixedBuf rc_stats_buffer;
    int            pass;

    // Deblock Filter

    /**
     * @brief Deblocking loop filter control
     *
     * Default is true.
     */
    EbBool enable_dlf_flag;

    /* Film grain denoising the input picture
    * Flag to enable the denoising
    *
    * Default is 0. */
    uint32_t film_grain_denoise_strength;

    /* CDEF Level
    *
    * Default is -1. */
    int cdef_level;

    /* Restoration filtering
    *  enable/disable
    *  set Self-Guided (sg) mode
    *  set Wiener (wn) mode
    *
    * Default is -1. */
    int enable_restoration_filtering;
    /* motion field motion vector
    *
    *  Default is -1. */
    int enable_mfmv;

    // Rate Control

    /* Rate control mode.
     *
     * 0 = Constant QP.
     * 1 = Variable Bit Rate, achieve the target bitrate at entire stream.
     * 2 = Constrained Variable Bit Rate, achieve the target bitrate at each gop
     * Default is 0. */
    uint32_t rate_control_mode;
    /* Flag to enable the scene change detection algorithm.
     *
     * Default is 1. */
    uint32_t scene_change_detection;

    /* When RateControlMode is set to 1 it's best to set this parameter to be
     * equal to the Intra period value (such is the default set by the encoder).
     * When CQP is chosen, then a (2 * minigopsize +1) look ahead is recommended.
     *
     * Default depends on rate control mode.*/
    uint32_t look_ahead_distance;

    /* Enable TPL in look ahead
     * 0 = disable TPL in look ahead
     * 1 = enable TPL in look ahead
     * Default is 0  */
    uint8_t enable_tpl_la;

    /* Target bitrate in bits/second, only apllicable when rate control mode is
     * set to 2 or 3.
     *
     * Default is 7000000. */
    uint32_t target_bit_rate;
    /* maximum bitrate in bits/second, only apllicable when rate control mode is
     * set to 0.
     *
     * Default is 0. */
    uint32_t max_bit_rate;
    /* VBV Buffer size */
    uint32_t vbv_bufsize;

    /* Maxium QP value allowed for rate control use, only applicable when rate
     * control mode is set to 1. It has to be greater or equal to minQpAllowed.
     *
     * Default is 63. */
    uint32_t max_qp_allowed;
    /* Minimum QP value allowed for rate control use, only applicable when rate
     * control mode is set to 1. It has to be smaller or equal to maxQpAllowed.
     *
     * Default is 0. */
    uint32_t min_qp_allowed;

    /* TWO PASS DATARATE CONTROL OPTIONS.
     * Indicates the bias (expressed on a scale of 0 to 100) for determining
     * target size for the current frame. The value 0 indicates the optimal CBR
     * mode value should be used, and 100 indicates the optimal VBR mode value
     * should be used. */
    uint32_t vbr_bias_pct;
    /* Indicates the minimum bitrate to be used for a single GOP as a percentage
     * of the target bitrate. */
    uint32_t vbr_min_section_pct;
    /* Indicates the maximum bitrate to be used for a single GOP as a percentage
     * of the target bitrate. */
    uint32_t vbr_max_section_pct;
    /* under_shoot_pct indicates the tolerance of the VBR algorithm to undershoot
     * and is used as a trigger threshold for more agressive adaptation of Q. Its
     * value can range from 0-100. */
    uint32_t under_shoot_pct;
    /* over_shoot_pct indicates the tolerance of the VBR algorithm to overshoot
     * and is used as a trigger threshold for more agressive adaptation of Q. Its
     * value can range from 0-1000. */
    uint32_t over_shoot_pct;
    /* over_shoot_pct indicates the tolerance of the Capped CRF algorithm to overshoot
     * and is used as a trigger threshold for more agressive adaptation of Q. Its
     * value can range from 0-1000. */
    uint32_t mbr_over_shoot_pct;
    /* Indicates the amount of data that will be buffered by the decoding
     * application prior to beginning playback, and is expressed in units of
     * time(milliseconds). */
    int64_t starting_buffer_level_ms;
    /* Indicates the amount of data that the encoder should try to maintain in the
     * decoder's buffer, and is expressed in units of time(milliseconds). */
    int64_t optimal_buffer_level_ms;
    /* Indicates the maximum amount of data that may be buffered by the decoding
     * application, and is expressed in units of time(milliseconds).*/
    int64_t maximum_buffer_size_ms;

    /* recode_loop indicates the recode levels,
     * DISALLOW_RECODE = 0, No recode.
     * ALLOW_RECODE_KFMAXBW = 1, Allow recode for KF and exceeding maximum frame bandwidth.
     * ALLOW_RECODE_KFARFGF = 2, Allow recode only for KF/ARF/GF frames.
     * ALLOW_RECODE = 3, Allow recode for all frames based on bitrate constraints.
     * ALLOW_RECODE_DEFAULT = 4, Default setting, ALLOW_RECODE_KFARFGF for M0~5 and
     *                                            ALLOW_RECODE_KFMAXBW for M6~8.
     * default is 4
     */
    uint32_t recode_loop;

    /* Flag to signal the content being a screen sharing content type
    *
    * Default is 0. */
    uint32_t screen_content_mode;

    /* Enable adaptive quantization within a frame using segmentation.
     *
     * Default is 2. */
    uint8_t enable_adaptive_quantization;

    // Tresholds

    /**
     * @brief Enable writing of HDR metadata in the bitstream
     *
     * Default is false.
     */
    EbBool high_dynamic_range_input;

    /**
     * @brief Bitstream profile to use.
     * 0: main, 1: high, 2: professional.
     *
     * Min is MAIN_PROFILE.
     * Max is PROFESSIONAL_PROFILE.
     * Default is MAIN_PROFILE.
     */
    EbAv1SeqProfile profile;
    /* Constraints for bitstream in terms of max bitrate and max buffer size.
     *
     * 0 = Main, for most applications.
     * 1 = High, for demanding applications.
     *
     * Default is 0. */
    uint32_t tier;

    /**
     * @brief Bitstream level.
     * 0: autodetect from bitstream, 20: level 2.0, 63: level 6.3, only levels 2.0-6.3 are properly defined.
     * The levels are defined at https://aomediacodec.github.io/av1-spec/av1-spec.pdf
     * under "A.3. Levels".
     *
     * Min is 0.
     * Max is 73.
     * Default is 0.
     */
    uint32_t level;

    /* CPU FLAGS to limit assembly instruction set used by encoder.
    * Default is CPU_FLAGS_ALL. */
    CPU_FLAGS use_cpu_flags;

    // Application Specific parameters

    /**
     * @brief API signal for the library to know the channel ID (used for pinning to cores).
     *
     * Min value is 0.
     * Max value is 0xFFFFFFFF.
     * Default is 0.
     */
    uint32_t channel_id;

    /**
     * @brief API signal for the library to know the active number of channels being encoded simultaneously.
     *
     * Min value is 1.
     * Max value is 0xFFFFFFFF.
     * Default is 1.
     */
    uint32_t active_channel_count;

    /**
     * @brief API signal to constrain motion vectors.
     *
     * Default is false.
     */
    EbBool restricted_motion_vector;

    // Threads management

    /* The number of logical processor which encoder threads run on. If
     * LogicalProcessors and TargetSocket are not set, threads are managed by
     * OS thread scheduler. */
    uint32_t logical_processors;

    /* Unpin the execution .This option does not
    * set the execution to be pinned to a specific number of cores when set to 1. this allows the execution
    * of multiple encodes on the CPU without having to pin them to a specific mask
    * 1: pinned threads
    * 0: unpinned
    * default 0 */
    uint32_t pin_threads;

    /* Target socket to run on. For dual socket systems, this can specify which
     * socket the encoder runs on.
     *
     * -1 = Both Sockets.
     *  0 = Socket 0.
     *  1 = Socket 1.
     *
     * Default is -1. */
    int32_t target_socket;

    // Debug tools

    /**
     * @brief API Signal to output reconstructed yuv used for debug purposes.
     * Using this will affect the speed of encoder.
     *
     * Default is false.
     */
    EbBool recon_enabled;

    /* Log 2 Tile Rows and colums . 0 means no tiling,1 means that we split the dimension
        * into 2
        * Default is 0. */
    int32_t tile_columns;
    int32_t tile_rows;

    /**
     * @brief Enable use of ALT-REF (temporally filtered) frames.
     *
     * Default is true.
     */
    EbBool enable_tf;

    EbBool enable_overlays;
    // super-resolution parameters
    uint8_t superres_mode;
    uint8_t superres_denom;
    uint8_t superres_kf_denom;
    uint8_t superres_qthres;
    uint8_t superres_kf_qthres;
    uint8_t superres_auto_search_type;

    /**
     * @brief API signal containing the manual prediction structure parameters.
     * Only used when enable_manual_pred_struct is enabled. This list is copied
     * into internal buffers after svt_av1_enc_set_parameter().
     */
    PredictionStructureConfigEntry pred_struct[1 << (MAX_HIERARCHICAL_LEVEL - 1)];

    /**
     * @brief API signal to overwrite the encoder's default prediction structure.
     *
     * Default is false.
     */
    EbBool enable_manual_pred_struct;

    /**
     * @brief API signal specifying the size (number of entries) of the manual prediction structure buffer.
     * Only checked and used when enable_manual_pred_struct is enabled.
     *
     * Min is 1.
     * Max is 32.
     * Default is 0.
     */
    int32_t manual_pred_struct_entry_num;

    // Color description
    /* Color description present flag
    *
    * It is not necessary to set this parameter manually.
    * It is set internally to true once one of the color_primaries, transfer_characteristics or
    * matrix coefficients is set to non-default value.
    *
    Default is false. */
    EbBool color_description_present_flag;
    /* Color primaries
    * values are from EbColorPrimaries
    Default is 2 (CP_UNSPECIFIED). */
    uint8_t color_primaries;
    /* Transfer characteristics
    * values are from EbTransferCharacteristics
    Default is 2 (TC_UNSPECIFIED). */
    uint8_t transfer_characteristics;
    /* Matrix coefficients
    * values are from EbMatrixCoefficients
    Default is 2 (MC_UNSPECIFIED). */
    uint8_t matrix_coefficients;
    /* Color range
    * values are from EbColorRange
    * 0: studio swing.
    * 1: full swing.
    Default is 0. */
    uint8_t color_range;
    /* Mastering display metadata
    * values are from set using svt_aom_parse_mastering_display()
    */
    struct EbSvtAv1MasteringDisplayInfo mastering_display;
    /* Content light level
    * values are from set using svt_aom_parse_content_light_level()
    */
    struct EbContentLightLevel content_light_level;
} EbSvtAv1EncConfiguration;

/**
 * Returns a string containing "v$tag-$commit_count-g$hash${dirty:+-dirty}"
 * @param[out] SVT_AV1_CVS_VERSION
 */
EB_API const char *svt_av1_get_version(void);

/**
 * Prints the version header and build information to the file
 * specified by the SVT_LOG_FILE environment variable or stderr
 */
EB_API void svt_av1_print_version(void);

/* STEP 1: Call the library to construct a Component Handle.
     *
     * Parameter:
     * @ **p_handle      Handle to be called in the future for manipulating the
     *                  component.
     * @ *p_app_data      Callback data.
     * @ *config_ptr     Pointer passed back to the client during callbacks, it will be
     *                  loaded with default params from the library. */
EB_API EbErrorType svt_av1_enc_init_handle(
    EbComponentType **p_handle, void *p_app_data,
    EbSvtAv1EncConfiguration
        *config_ptr); // config_ptr will be loaded with default params from the library

/* STEP 2: Set all configuration parameters.
     *
     * Parameter:
     * @ *svt_enc_component              Encoder handler.
     * @ *pComponentParameterStructure  Encoder and buffer configurations will be copied to the library. */
EB_API EbErrorType svt_av1_enc_set_parameter(
    EbComponentType *svt_enc_component,
    EbSvtAv1EncConfiguration *
        pComponentParameterStructure); // pComponentParameterStructure contents will be copied to the library

/* STEP 3: Initialize encoder and allocates memory to necessary buffers.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_init(EbComponentType *svt_enc_component);

/* OPTIONAL: Get stream headers at init time.
     *
     * Parameter:
     * @ *svt_enc_component   Encoder handler.
     * @ **output_stream_ptr  Output buffer. */
EB_API EbErrorType svt_av1_enc_stream_header(EbComponentType *    svt_enc_component,
                                             EbBufferHeaderType **output_stream_ptr);

/* OPTIONAL: Release stream headers at init time.
     *
     * Parameter:
     * @ *stream_header_ptr  stream header buffer. */
EB_API EbErrorType svt_av1_enc_stream_header_release(EbBufferHeaderType *stream_header_ptr);

/* STEP 4: Send the picture.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *p_buffer           Header pointer, picture buffer. */
EB_API EbErrorType svt_av1_enc_send_picture(EbComponentType *   svt_enc_component,
                                            EbBufferHeaderType *p_buffer);

/* STEP 5: Receive packet.
     * Parameter:
    * @ *svt_enc_component  Encoder handler.
     * @ **p_buffer          Header pointer to return packet with.
     * @ pic_send_done       Flag to signal that all input pictures have been sent, this call becomes locking one this signal is 1.
     * Non-locking call, returns EB_ErrorMax for an encode error, EB_NoErrorEmptyQueue when the library does not have any available packets.*/
EB_API EbErrorType svt_av1_enc_get_packet(EbComponentType *    svt_enc_component,
                                          EbBufferHeaderType **p_buffer, uint8_t pic_send_done);

/* STEP 5-1: Release output buffer back into the pool.
     *
     * Parameter:
     * @ **p_buffer          Header pointer that contains the output packet to be released. */
EB_API void svt_av1_enc_release_out_buffer(EbBufferHeaderType **p_buffer);

/* OPTIONAL: Fill buffer with reconstructed picture.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *p_buffer           Output buffer. */
EB_API EbErrorType svt_av1_get_recon(EbComponentType *   svt_enc_component,
                                     EbBufferHeaderType *p_buffer);

/* OPTIONAL: get stream information
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ *stream_info_id SVT_AV1_STREAM_INFO_ID.
     * @ *info         output, the type depends on id */
EB_API EbErrorType svt_av1_enc_get_stream_info(EbComponentType *svt_enc_component,
                                               uint32_t stream_info_id, void *info);

/* STEP 6: Deinitialize encoder library.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_deinit(EbComponentType *svt_enc_component);

/* STEP 7: Deconstruct encoder handler.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler. */
EB_API EbErrorType svt_av1_enc_deinit_handle(EbComponentType *svt_enc_component);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbSvtAv1Enc_h
