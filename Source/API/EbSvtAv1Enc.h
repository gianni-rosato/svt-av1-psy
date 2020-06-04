/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbSvtAv1Enc_h
#define EbSvtAv1Enc_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "stdint.h"
#include "EbSvtAv1.h"
#include <stdlib.h>
#include <stdio.h>
//***HME***
#define EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT 2
#define EB_HME_SEARCH_AREA_ROW_MAX_COUNT 2
#define MAX_HIERARCHICAL_LEVEL 6
#define REF_LIST_MAX_DEPTH 4

#define MAX_ENC_PRESET 8

#define DEFAULT -1

#define EB_BUFFERFLAG_EOS 0x00000001 // signals the last packet of the stream
#define EB_BUFFERFLAG_SHOW_EXT \
    0x00000002 // signals that the packet contains a show existing frame at the end
#define EB_BUFFERFLAG_HAS_TD 0x00000004 // signals that the packet contains a TD
#define EB_BUFFERFLAG_IS_ALT_REF 0x00000008 // signals that the packet contains an ALT_REF frame
#define EB_BUFFERFLAG_ERROR_MASK \
    0xFFFFFFF0 // mask for signalling error assuming top flags fit in 4 bits. To be changed, if more flags are added.

/************************************************
 * Prediction Structure Config Entry
 *   Contains the basic reference lists and
 *   configurations for each Prediction Structure
 *   Config Entry.
 ************************************************/
typedef struct PredictionStructureConfigEntry {
  uint32_t temporal_layer_index;
  uint32_t decode_order;
  int32_t ref_list0[REF_LIST_MAX_DEPTH];
  int32_t ref_list1[REF_LIST_MAX_DEPTH];
} PredictionStructureConfigEntry;

// super-res modes
typedef enum {
    SUPERRES_NONE,     // No frame superres allowed.
    SUPERRES_FIXED,    // All frames are coded at the specified scale, and super-resolved.
    SUPERRES_RANDOM,   // All frames are coded at a random scale, and super-resolved.
    SUPERRES_QTHRESH,  // Superres scale for a frame is determined based on q_index.
    SUPERRES_AUTO,     // Automatically select superres for appropriate frames.
    SUPERRES_MODES
} SUPERRES_MODE;

// Will contain the EbEncApi which will live in the EncHandle class
// Only modifiable during config-time.
typedef struct EbSvtAv1EncConfiguration {
    // Encoding preset

    /* A preset defining the quality vs density tradeoff point that the encoding
     * is to be performed at. 0 is the highest quality mode, 3 is the highest
     * density mode.
     *
     * Default is defined as MAX_ENC_PRESET. */
    uint8_t enc_mode;
    /* For two pass encoding, the enc_mod of the second pass is passed in the first pass.
    * First pass has the option to run with second pass ME settings.
    *
    * Default is defined as MAX_ENC_PRESET. */
    uint8_t snd_pass_enc_mode;
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
    uint32_t intra_refresh_type;
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
    /* The width of input source in units of picture luma pixels.
     *
     * Default is 0. */
    uint32_t source_width;
    /* The height of input source in units of picture luma pixels.
     *
     * Default is 0. */
    uint32_t source_height;

    uint32_t render_width, render_height;

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
    /* Specifies whether to use 16bit pipeline.
     *
     * 0: 8 bit pipeline.
     * 1: 16 bit pipeline.
     * Now 16bit pipeline is only enabled in filter
     * Default is 0. */
    EbBool is_16bit_pipeline;
    /* Specifies the chroma subsampleing format of input video.
     *
     * 0 = mono.
     * 1 = 420.
     * 2 = 422.
     * 3 = 444.
     *
     * Default is 1. */
    EbColorFormat encoder_color_format;
    /* Offline packing of the 2bits: requires two bits packed input.
     *
     * Default is 0. */
    uint32_t compressed_ten_bit_format;

    /* Super block size for motion estimation
    *
    * Default is 64. */
    uint32_t sb_sz;

    /* Super block size (mm-signal)
    *
    * Default is 128. */
    uint32_t super_block_size;
    /* The maximum partitioning depth with 0 being the superblock depth
    *
    * Default is 4. */
    uint32_t partition_depth;

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
    /* Input stats file */
    FILE *input_stat_file;
    /* output stats file */
    FILE *output_stat_file;
    /* Enable picture QP scaling between hierarchical levels
    *
    * Default is null.*/
    uint32_t enable_qp_scaling_flag;

    // Deblock Filter
    /* Flag to disable the Deblocking Loop Filtering.
     *
     * Default is 0. */
    EbBool disable_dlf_flag;

    /* Denoise the input picture when noise levels are too high
    * Flag to enable the denoising
    *
    * Default is 0. */
    EbBool enable_denoise_flag;

    /* Film grain denoising the input picture
    * Flag to enable the denoising
    *
    * Default is 0. */
    uint32_t film_grain_denoise_strength;

    /* Warped motion
    *
    * Default is -1. */
    int enable_warped_motion;

    /* Global motion
    *
    * Default is 1. */
    EbBool enable_global_motion;

    /* CDEF mode
    *
    * Default is -1. */
    int cdef_mode;

    /* Restoration filtering
    *  enable/disable
    *  set Self-Guided (sg) mode
    *  set Wiener (wn) mode
    *
    * Default is -1. */
    int enable_restoration_filtering;
    int sg_filter_mode;
    int wn_filter_mode;

    /* edge based skip angle intra
    *
    * Default is -1. */
    int edge_skp_angle_intra;

    /* enable angle intra
    *
    * Default is -1. */
    int intra_angle_delta;

    /* inter intra compound
    *
    * Default is -1. */
    int inter_intra_compound;

    /* enable paeth
    *
    * Default is -1. */
    int enable_paeth;

    /* enable smooth
    *
    * Default is -1. */
    int enable_smooth;

    /* combine class 12
    *
    * Default is -1. */
    int combine_class_12;

    /* motion field motion vector
    *
    *  Default is -1. */
    int enable_mfmv;
    /* redundant block
    *
    * Default is -1. */
    int enable_redundant_blk;
    /* spatial sse in full loop
    *
    * Default is -1. */
    int spatial_sse_fl;
    /* subpel
    *
    * Default is -1. */
    int enable_subpel;
    /* over boundry block
    *
    * Default is -1. */
    int over_bndry_blk;
    /* new nearest comb injection
    *
    * Default is -1. */
    int new_nearest_comb_inject;
    /* prune unipred at me
    *
    * Default is -1. */
    int prune_unipred_me;
    /* prune ref frame for rec partitions
    *
    * Default is -1. */
    int prune_ref_rec_part;
    /* nsq table
    *
    * Default is -1. */
    int nsq_table;
    /* frame end cdf update
    *
    * Default is -1. */
    int frame_end_cdf_update;

    /* Predictive Me
    *
    * Default is -1. */
    int pred_me;

    /* Bipred 3x3 Injection
    *
    * Default is -1. */
    int bipred_3x3_inject;

    /* Compound Mode
    *
    * Default is -1. */
    int compound_level;

    /* Chroma mode
    *
    * Level                Settings
    * CHROMA_MODE_0  0     Full chroma search @ MD
    * CHROMA_MODE_1  1     Fast chroma search @ MD
    * CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    * CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
    *
    * Default is -1 (AUTO) */
    int set_chroma_mode;

    /* Disable chroma from luma (CFL)
     *
     * Default is -1 (auto) */
    int disable_cfl_flag;

    /* OBMC
    *
    * Default is 1. */
    EbBool enable_obmc;

    /* RDOQ
    *
    * Default is -1. */
    int enable_rdoq;

    /* Filter intra prediction
    *
    * Default is 1. */
    EbBool enable_filter_intra;

    /* Intra Edge Filter
    *
    * Default is -1. */
    int enable_intra_edge_filter;

    /* Picture based rate estimation
    *
    * Default is - 1. */
    int pic_based_rate_est;

    /* Flag to enable the use of default ME HME parameters.
    *
    * Default is 1. */
    EbBool use_default_me_hme;

    /* Flag to enable Hierarchical Motion Estimation.
    *
    * Default is 1. */
    EbBool enable_hme_flag;

    /* Flag to enable the use of non-swaure partitions
    *
    * Default is 1. */
    EbBool ext_block_flag;

    /* Flag to enable the use of recon pictures for motion estimation
    *
    * Default is 1. */
    EbBool in_loop_me_flag;

    // ME Parameters
    /* Number of search positions in the horizontal direction.
     *
     * Default depends on input resolution. */
    uint32_t search_area_width;
    /* Number of search positions in the vertical direction.
     *
     * Default depends on input resolution. */
    uint32_t search_area_height;

    // MD Parameters
    /* Enable the use of HBD (10-bit) for 10 bit content at the mode decision step
     *
     * 0 = 8bit mode decision
     * 1 = 10bit mode decision
     * 2 = Auto: 8bit & 10bit mode decision
     *
    * Default is 1. */
    uint8_t enable_hbd_mode_decision;

    /* Palette Mode
    *
    * Default is -1. */
    int32_t enable_palette;

    // Rate Control

    /* Rate control mode.
     *
     * 0 = Constant QP.
     * 1 = Average BitRate.
     *
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

    /* Target bitrate in bits/second, only apllicable when rate control mode is
     * set to 2 or 3.
     *
     * Default is 7000000. */
    uint32_t target_bit_rate;

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

    /* Flag to signal the content being a screen sharing content type
    *
    * Default is 0. */
    uint32_t screen_content_mode;

    /* Flag to control intraBC mode
    *  0      OFF
    *  1      slow
    *  2      faster
    *  3      fastest
    *
    * Default is -1 (DEFAULT behavior). */
    int intrabc_mode;

    /* Enable adaptive quantization within a frame using segmentation.
     *
     * Default is FALSE. */
    EbBool enable_adaptive_quantization;

    // Tresholds
    /* Flag to signal that the input yuv is HDR10 BT2020 using SMPTE ST2048, requires
     *
     * Default is 0. */
    uint32_t high_dynamic_range_input;

    /* Defined set of coding tools to create bitstream.
     *
     * 1 = Main, allows bit depth of 8.
     * 2 = Main 10, allows bit depth of 8 to 10.
     *
     * Default is 2. */
    uint32_t profile;
    /* Constraints for bitstream in terms of max bitrate and max buffer size.
     *
     * 0 = Main, for most applications.
     * 1 = High, for demanding applications.
     *
     * Default is 0. */
    uint32_t tier;
    /* Constraints for bitstream in terms of max bitrate and max buffer size.
     *
     * 0 = auto determination.
     *
     * Default is 0. */
    uint32_t level;

    /* CPU FLAGS to limit assembly instruction set used by encoder.
    * Default is CPU_FLAGS_ALL. */
    CPU_FLAGS use_cpu_flags;

    // Application Specific parameters

    /* ID assigned to each channel when multiple instances are running within the
     * same application. */
    uint32_t channel_id;
    uint32_t active_channel_count;

    /* Flag to enable the Speed Control functionality to achieve the real-time
    * encoding speed defined by dynamically changing the encoding preset to meet
    * the average speed defined in injectorFrameRate. When this parameter is set
    * to 1 it forces -inj to be 1 -inj-frm-rt to be set to the -fps.
    *
    * Default is 0. */
    uint32_t speed_control_flag;

    /* Frame Rate used for the injector. Recommended to match the encoder speed.
    *
    * Default is 60. */
    int32_t injector_frame_rate;

    /* Flag to constrain motion vectors.
     *
     * 1: Motion vectors are allowed to point outside frame boundary.
     * 0: Motion vectors are NOT allowed to point outside frame boundary.
     *
     * Default is 1. */
    uint8_t unrestricted_motion_vector;

    // Threads management

    /* The number of logical processor which encoder threads run on. If
     * LogicalProcessorNumber and TargetSocket are not set, threads are managed by
     * OS thread scheduler. */
    uint32_t logical_processors;

    /* Unpin the execution .This option does not
    * set the execution to be pinned to a specific number of cores when set to 1. this allows the execution
    * of multiple encodes on the CPU wihtout having to pin them to a specific mask
    * 1: unpinned
    * 0: pinned
    *
    * If logical_processors is set to 1 default is 1 otherwise it is 0. */
    uint32_t unpin;

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

    /* Output reconstructed yuv used for debug purposes. The value is set through
     * ReconFile token (-o) and using the feature will affect the speed of encoder.
     *
     * Default is 0. */
    uint32_t recon_enabled;
    /* Log 2 Tile Rows and colums . 0 means no tiling,1 means that we split the dimension
        * into 2
        * Default is 0. */
    int32_t tile_columns;
    int32_t tile_rows;

    /* To be deprecated.
 * Encoder configuration parameters below this line are to be deprecated. */

    /* Flag to enable Hierarchical Motion Estimation 1/16th of the picture
    *
    * Default is 1. */
    EbBool enable_hme_level0_flag;

    /* Flag to enable Hierarchical Motion Estimation 1/4th of the picture
    *
    * Default is 1. */
    EbBool enable_hme_level1_flag;

    /* Flag to enable Hierarchical Motion Estimation full sample of the picture
    *
    * Default is 1. */
    EbBool enable_hme_level2_flag;

    // HME Parameters
    /* Number of search positions in width and height for the HME
    *
    * Default depends on input resolution. */
    uint32_t number_hme_search_region_in_width;
    uint32_t number_hme_search_region_in_height;
    uint32_t hme_level0_total_search_area_width;
    uint32_t hme_level0_total_search_area_height;
    uint32_t hme_level0_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level0_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t hme_level1_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level1_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t hme_level2_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level2_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];

    uint32_t ten_bit_format;

    /* Variables to control the use of ALT-REF (temporally filtered frames)
    */
    EbBool  enable_altrefs;
    uint8_t altref_strength;
    uint8_t altref_nframes;
    EbBool  enable_overlays;

    // super-resolution parameters
    uint8_t superres_mode;
    uint8_t superres_denom;
    uint8_t superres_kf_denom;
    uint8_t superres_qthres;

    uint32_t sq_weight;

    uint64_t md_stage_1_cand_prune_th;
    uint64_t md_stage_1_class_prune_th;
    uint64_t md_stage_2_3_cand_prune_th;
    uint64_t md_stage_2_3_class_prune_th;

  /* Prediction Structure user defined
   */
  PredictionStructureConfigEntry pred_struct[1 << (MAX_HIERARCHICAL_LEVEL - 1)];
  /* Flag to enable use prediction structure user defined
   *
   * Default is false. */
  EbBool enable_manual_pred_struct;
  /* The minigop size of prediction structure user defined
   *
   * Default is 0. */
  int32_t manual_pred_struct_entry_num;
} EbSvtAv1EncConfiguration;

/* STEP 1: Call the library to construct a Component Handle.
     *
     * Parameter:
     * @ **p_handle      Handle to be called in the future for manipulating the
     *                  component.
     * @ *p_app_data      Callback data.
     * @ *config_ptr     Pointer passed back to the client during callbacks, it will be
     *                  loaded with default params from the library. */
EB_API EbErrorType
svt_av1_enc_init_handle(EbComponentType **p_handle, void *p_app_data,
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


/* OPTIONAL: Get the end of sequence Network Abstraction Layer.
     *
     * Parameter:
     * @ *svt_enc_component  Encoder handler.
     * @ **output_stream_ptr  Output stream. */
EB_API EbErrorType svt_av1_enc_eos_nal(EbComponentType *    svt_enc_component,
                                      EbBufferHeaderType **output_stream_ptr);

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
