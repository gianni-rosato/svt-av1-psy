/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPackUnPack_C.h"
#include "EbSvtAv1.h"

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
************************************************/
void eb_enc_msb_pack2_d(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                        uint16_t *out16_bit_buffer, uint32_t inn_stride, uint32_t out_stride,
                        uint32_t width, uint32_t height) {
    uint64_t j, k;
    uint16_t out_pixel;
    uint8_t  n_bit_pixel;

    //SIMD hint: use _mm_unpacklo_epi8 +_mm_unpackhi_epi8 to do the concatenation

    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            out_pixel   = (in8_bit_buffer[k + j * in8_stride]) << 2;
            n_bit_pixel = (inn_bit_buffer[k + j * inn_stride] >> 6) & 3;

            out16_bit_buffer[k + j * out_stride] = out_pixel | n_bit_pixel;
        }
    }
}

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void compressed_packmsb_c(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                          uint16_t *out16_bit_buffer, uint32_t inn_stride, uint32_t out_stride,
                          uint32_t width, uint32_t height) {
    uint64_t row, k_idx;
    uint16_t out_pixel;
    uint8_t  n_bit_pixel;
    uint8_t  four_2bit_pels;

    for (row = 0; row < height; row++) {
        for (k_idx = 0; k_idx < width / 4; k_idx++) {
            four_2bit_pels = inn_bit_buffer[k_idx + row * inn_stride];

            n_bit_pixel = (four_2bit_pels >> 6) & 3;

            out_pixel = in8_bit_buffer[k_idx * 4 + 0 + row * in8_stride] << 2;
            out16_bit_buffer[k_idx * 4 + 0 + row * out_stride] = out_pixel | n_bit_pixel;

            n_bit_pixel = (four_2bit_pels >> 4) & 3;
            out_pixel   = in8_bit_buffer[k_idx * 4 + 1 + row * in8_stride] << 2;
            out16_bit_buffer[k_idx * 4 + 1 + row * out_stride] = out_pixel | n_bit_pixel;

            n_bit_pixel = (four_2bit_pels >> 2) & 3;
            out_pixel   = in8_bit_buffer[k_idx * 4 + 2 + row * in8_stride] << 2;
            out16_bit_buffer[k_idx * 4 + 2 + row * out_stride] = out_pixel | n_bit_pixel;

            n_bit_pixel = (four_2bit_pels >> 0) & 3;
            out_pixel   = in8_bit_buffer[k_idx * 4 + 3 + row * in8_stride] << 2;
            out16_bit_buffer[k_idx * 4 + 3 + row * out_stride] = out_pixel | n_bit_pixel;
        }
    }
}

/************************************************
* convert unpacked nbit (n=2) data to compressedPAcked
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void c_pack_c(const uint8_t *inn_bit_buffer, uint32_t inn_stride, uint8_t *in_compn_bit_buffer,
              uint32_t out_stride, uint8_t *local_cache, uint32_t width, uint32_t height) {
    uint32_t row_index, col_index;
    (void)local_cache;

    for (row_index = 0; row_index < height; row_index++) {
        for (col_index = 0; col_index < width; col_index += 4) {
            uint32_t i = col_index + row_index * inn_stride;

            uint8_t compressed_unpacked_pixel = 0;
            compressed_unpacked_pixel =
                compressed_unpacked_pixel | ((inn_bit_buffer[i + 0] >> 0) & 0xC0); //1100.0000
            compressed_unpacked_pixel =
                compressed_unpacked_pixel | ((inn_bit_buffer[i + 1] >> 2) & 0x30); //0011.0000
            compressed_unpacked_pixel =
                compressed_unpacked_pixel | ((inn_bit_buffer[i + 2] >> 4) & 0x0C); //0000.1100
            compressed_unpacked_pixel =
                compressed_unpacked_pixel | ((inn_bit_buffer[i + 3] >> 6) & 0x03); //0000.0011

            uint32_t j             = col_index / 4 + row_index * out_stride;
            in_compn_bit_buffer[j] = compressed_unpacked_pixel;
        }
    }
}

/************************************************
* unpack 10 bit data into  8 and 2 bit 2D data
************************************************/
void eb_enc_msb_un_pack2_d(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer,
                           uint8_t *outn_bit_buffer, uint32_t out8_stride, uint32_t outn_stride,
                           uint32_t width, uint32_t height) {
    uint64_t j, k;
    uint16_t in_pixel;
    uint8_t  tmp_pixel;
    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            in_pixel                             = in16_bit_buffer[k + j * in_stride];
            out8_bit_buffer[k + j * out8_stride] = (uint8_t)(in_pixel >> 2);
            tmp_pixel                            = (uint8_t)(in_pixel << 6);
            outn_bit_buffer[k + j * outn_stride] = tmp_pixel;
        }
    }
}
void un_pack8_bit_data_c(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer,
                         uint32_t out8_stride, uint32_t width, uint32_t height) {
    uint64_t j, k;
    uint16_t in_pixel;
    //uint8_t    tmp_pixel;
    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            in_pixel                             = in16_bit_buffer[k + j * in_stride];
            out8_bit_buffer[k + j * out8_stride] = (uint8_t)(in_pixel >> 2);
            //tmp_pixel = (uint8_t)(in_pixel << 6);
            //outn_bit_buffer[k + j*outn_stride] = tmp_pixel;
        }
    }
}
void unpack_avg_c(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1,
                  uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride, uint32_t width,
                  uint32_t height) {
    uint64_t j, k;
    uint8_t  in_pixel_l0, in_pixel_l1;

    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            in_pixel_l0                 = (uint8_t)(ref16_l0[k + j * ref_l0_stride] >> 2);
            in_pixel_l1                 = (uint8_t)(ref16_l1[k + j * ref_l1_stride] >> 2);
            dst_ptr[k + j * dst_stride] = (in_pixel_l0 + in_pixel_l1 + 1) >> 1;
        }
    }
}

void unpack_avg_safe_sub_c(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1,
                           uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride,
                           EbBool sub_pred, uint32_t width, uint32_t height) {
    uint64_t j, k;
    uint8_t  in_pixel_l0, in_pixel_l1;

    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            in_pixel_l0                 = (uint8_t)(ref16_l0[k + j * ref_l0_stride] >> 2);
            in_pixel_l1                 = (uint8_t)(ref16_l1[k + j * ref_l1_stride] >> 2);
            dst_ptr[k + j * dst_stride] = (in_pixel_l0 + in_pixel_l1 + 1) >> 1;
        }
    }

    if (sub_pred) {
        //Last row
        j = height * 2 - 1;
        for (k = 0; k < width; k++) {
            in_pixel_l0                     = (uint8_t)(ref16_l0[k + j * ref_l0_stride / 2] >> 2);
            in_pixel_l1                     = (uint8_t)(ref16_l1[k + j * ref_l1_stride / 2] >> 2);
            dst_ptr[k + j * dst_stride / 2] = (in_pixel_l0 + in_pixel_l1 + 1) >> 1;
        }
    }
}
void convert_8bit_to_16bit_c(uint8_t* src, uint32_t src_stride, uint16_t* dst, uint32_t dst_stride,
    uint32_t width, uint32_t height) {
    for (uint32_t j = 0; j < height; j++) {
        for (uint32_t k = 0; k < width; k++) {
            dst[k + j * dst_stride] = src[k + j * src_stride];
        }
    }
}

void convert_16bit_to_8bit_c(uint16_t *src, uint32_t src_stride, uint8_t *dst, uint32_t dst_stride,
    uint32_t width, uint32_t height) {
    for (uint32_t j = 0; j < height; j++) {
        for (uint32_t k = 0; k < width; k++) {
            dst[k + j * dst_stride] = (uint8_t)(src[k + j * src_stride]);
        }
    }
}
