/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureBuffer_h
#define EbPictureBuffer_h

#include <stdio.h>

#include "EbDefinitions.h"
#include "grainSynthesis.h"


#ifdef __cplusplus
extern "C" {
#endif
#define PICTURE_BUFFER_DESC_Y_FLAG              (1 << 0)
#define PICTURE_BUFFER_DESC_Cb_FLAG             (1 << 1)
#define PICTURE_BUFFER_DESC_Cr_FLAG             (1 << 2)
#define PICTURE_BUFFER_DESC_LUMA_MASK           PICTURE_BUFFER_DESC_Y_FLAG
#define PICTURE_BUFFER_DESC_CHROMA_MASK         (PICTURE_BUFFER_DESC_Cb_FLAG | PICTURE_BUFFER_DESC_Cr_FLAG)
#define PICTURE_BUFFER_DESC_FULL_MASK           (PICTURE_BUFFER_DESC_Y_FLAG | PICTURE_BUFFER_DESC_Cb_FLAG | PICTURE_BUFFER_DESC_Cr_FLAG)

    /************************************
     * EbPictureBufferDesc
     ************************************/
    typedef struct EbPictureBufferDesc_s
    {
        // Buffer Ptrs
        EbByte         buffer_y;             // pointer to the Y luma buffer
        EbByte         bufferCb;            // pointer to the U chroma buffer
        EbByte         bufferCr;            // pointer to the V chroma buffer
        //Bit increment
        EbByte         bufferBitIncY;       // pointer to the Y luma buffer Bit increment
        EbByte         bufferBitIncCb;      // pointer to the U chroma buffer Bit increment
        EbByte         bufferBitIncCr;      // pointer to the V chroma buffer Bit increment

        uint16_t          stride_y;          // pointer to the Y luma buffer
        uint16_t          strideCb;         // pointer to the U chroma buffer
        uint16_t          strideCr;         // pointer to the V chroma buffer

        uint16_t          strideBitIncY;    // pointer to the Y luma buffer Bit increment
        uint16_t          strideBitIncCb;   // pointer to the U chroma buffer Bit increment
        uint16_t          strideBitIncCr;   // pointer to the V chroma buffer Bit increment

        // Picture Parameters
        uint16_t          origin_x;         // Horizontal padding distance
        uint16_t          origin_y;         // Vertical padding distance
        uint16_t          width;            // Luma picture width which excludes the padding
        uint16_t          height;           // Luma picture height which excludes the padding
        uint16_t          maxWidth;         // input Luma picture width
        uint16_t          maxHeight;        // input Luma picture height
        EB_BITDEPTH       bit_depth;        // Pixel Bit Depth

        // Buffer Parameters
        uint32_t          lumaSize;         // Size of the luma buffer
        uint32_t          chromaSize;       // Size of the chroma buffers
        EbBool            packedFlag;       // Indicates if sample buffers are packed or not

        EbBool            film_grain_flag;  // Indicates if film grain parameters are present for the frame

    } EbPictureBufferDesc_t;



#define YV12_FLAG_HIGHBITDEPTH 8

    /*!\brief List of supported color primaries */
    typedef enum aom_color_primaries {
        AOM_CICP_CP_RESERVED_0 = 0,  /**< For future use */
        AOM_CICP_CP_BT_709 = 1,      /**< BT.709 */
        AOM_CICP_CP_UNSPECIFIED = 2, /**< Unspecified */
        AOM_CICP_CP_RESERVED_3 = 3,  /**< For future use */
        AOM_CICP_CP_BT_470_M = 4,    /**< BT.470 System M (historical) */
        AOM_CICP_CP_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
        AOM_CICP_CP_BT_601 = 6,      /**< BT.601 */
        AOM_CICP_CP_SMPTE_240 = 7,   /**< SMPTE 240 */
        AOM_CICP_CP_GENERIC_FILM =
        8, /**< Generic film (color filters using illuminant C) */
        AOM_CICP_CP_BT_2020 = 9,      /**< BT.2020, BT.2100 */
        AOM_CICP_CP_XYZ = 10,         /**< SMPTE 428 (CIE 1921 XYZ) */
        AOM_CICP_CP_SMPTE_431 = 11,   /**< SMPTE RP 431-2 */
        AOM_CICP_CP_SMPTE_432 = 12,   /**< SMPTE EG 432-1  */
        AOM_CICP_CP_RESERVED_13 = 13, /**< For future use (values 13 - 21)  */
        AOM_CICP_CP_EBU_3213 = 22,    /**< EBU Tech. 3213-E  */
        AOM_CICP_CP_RESERVED_23 = 23  /**< For future use (values 23 - 255)  */
    } aom_color_primaries_t;        /**< alias for enum aom_color_primaries */

    /*!\brief List of supported transfer functions */
    typedef enum aom_transfer_characteristics {
        AOM_CICP_TC_RESERVED_0 = 0,  /**< For future use */
        AOM_CICP_TC_BT_709 = 1,      /**< BT.709 */
        AOM_CICP_TC_UNSPECIFIED = 2, /**< Unspecified */
        AOM_CICP_TC_RESERVED_3 = 3,  /**< For future use */
        AOM_CICP_TC_BT_470_M = 4,    /**< BT.470 System M (historical)  */
        AOM_CICP_TC_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
        AOM_CICP_TC_BT_601 = 6,      /**< BT.601 */
        AOM_CICP_TC_SMPTE_240 = 7,   /**< SMPTE 240 M */
        AOM_CICP_TC_LINEAR = 8,      /**< Linear */
        AOM_CICP_TC_LOG_100 = 9,     /**< Logarithmic (100 : 1 range) */
        AOM_CICP_TC_LOG_100_SQRT10 =
        10,                     /**< Logarithmic (100 * Sqrt(10) : 1 range) */
        AOM_CICP_TC_IEC_61966 = 11, /**< IEC 61966-2-4 */
        AOM_CICP_TC_BT_1361 = 12,   /**< BT.1361 */
        AOM_CICP_TC_SRGB = 13,      /**< sRGB or sYCC*/
        AOM_CICP_TC_BT_2020_10_BIT = 14, /**< BT.2020 10-bit systems */
        AOM_CICP_TC_BT_2020_12_BIT = 15, /**< BT.2020 12-bit systems */
        AOM_CICP_TC_SMPTE_2084 = 16,     /**< SMPTE ST 2084, ITU BT.2100 PQ */
        AOM_CICP_TC_SMPTE_428 = 17,      /**< SMPTE ST 428 */
        AOM_CICP_TC_HLG = 18,            /**< BT.2100 HLG, ARIB STD-B67 */
        AOM_CICP_TC_RESERVED_19 = 19     /**< For future use (values 19-255) */
    } aom_transfer_characteristics_t;  /**< alias for enum aom_transfer_function */

    /*!\brief List of supported matrix coefficients */
    typedef enum aom_matrix_coefficients {
        AOM_CICP_MC_IDENTITY = 0,    /**< Identity matrix */
        AOM_CICP_MC_BT_709 = 1,      /**< BT.709 */
        AOM_CICP_MC_UNSPECIFIED = 2, /**< Unspecified */
        AOM_CICP_MC_RESERVED_3 = 3,  /**< For future use */
        AOM_CICP_MC_FCC = 4,         /**< US FCC 73.628 */
        AOM_CICP_MC_BT_470_B_G = 5,  /**< BT.470 System B, G (historical) */
        AOM_CICP_MC_BT_601 = 6,      /**< BT.601 */
        AOM_CICP_MC_SMPTE_240 = 7,   /**< SMPTE 240 M */
        AOM_CICP_MC_SMPTE_YCGCO = 8, /**< YCgCo */
        AOM_CICP_MC_BT_2020_NCL =
        9, /**< BT.2020 non-constant luminance, BT.2100 YCbCr  */
        AOM_CICP_MC_BT_2020_CL = 10, /**< BT.2020 constant luminance */
        AOM_CICP_MC_SMPTE_2085 = 11, /**< SMPTE ST 2085 YDzDx */
        AOM_CICP_MC_CHROMAT_NCL =
        12, /**< Chromaticity-derived non-constant luminance */
        AOM_CICP_MC_CHROMAT_CL = 13, /**< Chromaticity-derived constant luminance */
        AOM_CICP_MC_ICTCP = 14,      /**< BT.2100 ICtCp */
        AOM_CICP_MC_RESERVED_15 = 15 /**< For future use (values 15-255)  */
    } aom_matrix_coefficients_t;

    /*!\brief List of supported color range */
    typedef enum aom_color_range {
        AOM_CR_STUDIO_RANGE = 0, /**< Y [16..235], UV [16..240] */
        AOM_CR_FULL_RANGE = 1    /**< YUV/RGB [0..255] */
    } aom_color_range_t;       /**< alias for enum aom_color_range */

    /*!\brief List of chroma sample positions */
    typedef enum aom_chroma_sample_position {
        AOM_CSP_UNKNOWN = 0,          /**< Unknown */
        AOM_CSP_VERTICAL = 1,         /**< Horizontally co-located with luma(0, 0)*/
        /**< sample, between two vertical samples */
        AOM_CSP_COLOCATED = 2,        /**< Co-located with luma(0, 0) sample */
        AOM_CSP_RESERVED = 3          /**< Reserved value */
    } aom_chroma_sample_position_t; /**< alias for enum aom_transfer_function */
    typedef struct yv12_buffer_config {
        union {
            struct {
                int32_t y_width;
                int32_t uv_width;
                int32_t alpha_width;
            };
            int32_t widths[3];
        };
        union {
            struct {
                int32_t y_height;
                int32_t uv_height;
                int32_t alpha_height;
            };
            int32_t heights[3];
        };
        union {
            struct {
                int32_t y_crop_width;
                int32_t uv_crop_width;
            };
            int32_t crop_widths[2];
        };
        union {
            struct {
                int32_t y_crop_height;
                int32_t uv_crop_height;
            };
            int32_t crop_heights[2];
        };
        union {
            struct {
                int32_t y_stride;
                int32_t uv_stride;
                int32_t alpha_stride;
            };
            int32_t strides[3];
        };
        union {
            struct {
                uint8_t *y_buffer;
                uint8_t *u_buffer;
                uint8_t *v_buffer;
                uint8_t *alpha_buffer;
            };
            uint8_t *buffers[4];
        };

        // Indicate whether y_buffer, u_buffer, and v_buffer points to the internally
        // allocated memory or external buffers.
        int32_t use_external_refernce_buffers;
        // This is needed to store y_buffer, u_buffer, and v_buffer when set reference
        // uses an external refernece, and restore those buffer pointers after the
        // external reference frame is no longer used.
        uint8_t *store_buf_adr[3];

        // If the frame is stored in a 16-bit buffer, this stores an 8-bit version
        // for use in global motion detection. It is allocated on-demand.
        uint8_t *y_buffer_8bit;
        int32_t buf_8bit_valid;

        uint8_t *buffer_alloc;
        size_t buffer_alloc_sz;
        int32_t border;
        size_t frame_size;
        int32_t subsampling_x;
        int32_t subsampling_y;
        uint32_t bit_depth;
        aom_color_primaries_t color_primaries;
        aom_transfer_characteristics_t transfer_characteristics;
        aom_matrix_coefficients_t matrix_coefficients;
        int32_t monochrome;
        aom_chroma_sample_position_t chroma_sample_position;
        aom_color_range_t color_range;
        int32_t render_width;
        int32_t render_height;

        int32_t corrupted;
        int32_t flags;
    } Yv12BufferConfig;

    void LinkEbToAomBufferDesc(
        EbPictureBufferDesc_t          *picBuffDsc,
        Yv12BufferConfig             *aomBuffDsc
    );

    typedef struct aom_codec_frame_buffer {
        uint8_t *data; /**< pointer to the data buffer */
        size_t size;   /**< Size of data in bytes */
        void *priv;    /**< Frame's private data */
    } aom_codec_frame_buffer_t;

    /*!\brief get frame buffer callback prototype
    *
    * This callback is invoked by the decoder to retrieve data for the frame
    * buffer in order for the decode call to complete. The callback must
    * allocate at least min_size in bytes and assign it to fb->data. The callback
    * must zero out all the data allocated. Then the callback must set fb->size
    * to the allocated size. The application does not need to align the allocated
    * data. The callback is triggered when the decoder needs a frame buffer to
    * decode a compressed image into. This function may be called more than once
    * for every call to aom_codec_decode. The application may set fb->priv to
    * some data which will be passed back in the ximage and the release function
    * call. |fb| is guaranteed to not be NULL. On success the callback must
    * return 0. Any failure the callback must return a value less than 0.
    *
    * \param[in] priv         Callback's private data
    * \param[in] new_size     Size in bytes needed by the buffer
    * \param[in,out] fb       pointer to aom_codec_frame_buffer_t
    */
    typedef int32_t(*aom_get_frame_buffer_cb_fn_t)(void *priv, size_t min_size,
        aom_codec_frame_buffer_t *fb);

#define ADDRESS_STORAGE_SIZE sizeof(size_t)

    /*returns an addr aligned to the byte boundary specified by align*/
#define align_addr(addr, align) \
  (void *)(((size_t)(addr) + ((align)-1)) & ~(size_t)((align)-1))


#define AOM_BORDER_IN_PIXELS 288


/************************************
 * EbPictureBufferDesc Init Data
 ************************************/
    typedef struct EbPictureBufferDescInitData_s
    {
        uint16_t          maxWidth;
        uint16_t          maxHeight;
        EB_BITDEPTH       bit_depth;
        uint32_t          bufferEnableMask;
        uint16_t          left_padding;
        uint16_t          right_padding;
        uint16_t          top_padding;
        uint16_t          bot_padding;
        EbBool            splitMode;         //ON: allocate 8bit data seperately from nbit data

    } EbPictureBufferDescInitData_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType eb_picture_buffer_desc_ctor(
        EbPtr *object_dbl_ptr,
        EbPtr  object_init_data_ptr);

    extern EbErrorType eb_recon_picture_buffer_desc_ctor(
        EbPtr *object_dbl_ptr,
        EbPtr  object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureBuffer_h