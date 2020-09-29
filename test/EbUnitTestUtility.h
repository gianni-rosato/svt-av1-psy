/*
 * Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbUnitTestUtility_h
#define EbUnitTestUtility_h

#include <stddef.h>
#include <stdint.h>
#include "EbDefinitions.h"
#include "EbTime.h"

#ifdef __cplusplus
extern "C" {
#endif

#define Eb_UNIT_TEST_BUF_SIZE 0x04000000
#define Eb_UNIT_TEST_BUF_ALIGN_BYTE 256
#define Eb_UNIT_TEST_BUF_UNALIGN_BYTE 3
#define NELEMENTS(x) (int)(sizeof(x) / sizeof(x[0]))

extern void eb_buf_random_void(void *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_s8(int8_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_s8_with_max(int8_t *const buf, const uint32_t sizeBuf,
                                      const uint8_t max);
extern void eb_buf_random_u8(uint8_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_u8_to_0_or_255(uint8_t *const buf,
                                         const uint32_t sizeBuf);
extern void eb_buf_random_u8_to_255(uint8_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_u8_to_large(uint8_t *const buf,
                                      const uint32_t sizeBuf);
extern void eb_buf_random_u8_to_near_value(uint8_t *const buf,
                                           const uint32_t sizeBuf,
                                           const uint8_t val,
                                           const uint32_t range);
extern void eb_buf_random_u8_to_small(uint8_t *const buf,
                                      const uint32_t sizeBuf);
extern void eb_buf_random_u8_to_small_or_large(uint8_t *const buf,
                                               const uint32_t sizeBuf);
extern void eb_buf_random_u8_with_max(uint8_t *const buf,
                                      const uint32_t sizeBuf,
                                      const uint8_t max);
extern void eb_buf_random_s16(int16_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_s16_to_bd(int16_t *const buf, const uint32_t sizeBuf,
                                    const uint32_t bd);
extern void eb_buf_random_s16_with_bd(int16_t *const buf,
                                      const uint32_t sizeBuf,
                                      const uint32_t bd);
extern void eb_buf_random_s16_with_max(int16_t *const buf,
                                       const uint32_t sizeBuf,
                                       const uint16_t max);
extern void eb_buf_random_u16(uint16_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_u16_to_0_or_1023(uint16_t *const buf,
                                           const uint32_t sizeBuf);
extern void eb_buf_random_u16_to_0_or_bd(uint16_t *const buf,
                                         const uint32_t sizeBuf,
                                         const uint32_t bd);
extern void eb_buf_random_u16_to_bd(uint16_t *const buf, const uint32_t sizeBuf,
                                    const uint32_t bd);
extern void eb_buf_random_u16_with_bd(uint16_t *const buf,
                                      const uint32_t sizeBuf,
                                      const uint32_t bd);
extern void eb_buf_random_u16_with_max(uint16_t *const buf,
                                       const uint32_t sizeBuf,
                                       const uint32_t max);
extern void eb_buf_random_s32(int32_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_s32_with_max(int32_t *const buf,
                                       const uint32_t sizeBuf,
                                       const int32_t max);
extern void eb_buf_random_u32(uint32_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_u32_with_max(uint32_t *const buf,
                                       const uint32_t sizeBuf,
                                       const uint32_t max);
extern void eb_buf_random_s64(int64_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_u64(uint64_t *const buf, const uint32_t sizeBuf);
extern void eb_buf_random_ptrdiff_with_max(ptrdiff_t *const buf,
                                           const uint32_t sizeBuf,
                                           const ptrdiff_t max);

extern uint32_t eb_create_random_aligned_stride(const uint32_t width,
                                                const uint32_t align);

extern EbBool eb_buf_compare_u8(const uint8_t *const buf1,
                                const uint8_t *const buf2,
                                const size_t bufSize);
extern EbBool eb_buf_compare_s8(const int8_t *const buf1,
                                const int8_t *const buf2, const size_t bufSize);
extern EbBool eb_buf_compare_u16(const uint16_t *const buf1,
                                 const uint16_t *const buf2,
                                 const size_t bufSize);
extern EbBool eb_buf_compare_s16(const int16_t *const buf1,
                                 const int16_t *const buf2,
                                 const size_t bufSize);
extern EbBool eb_buf_compare_u32(const uint32_t *const buf1,
                                 const uint32_t *const buf2,
                                 const size_t bufSize);
extern EbBool eb_buf_compare_s32(const int32_t *const buf1,
                                 const int32_t *const buf2,
                                 const size_t bufSize);

extern EbBool eb_image_compare_u8(const uint8_t *buf1, const uint8_t *buf2,
                                  const size_t width, const size_t height,
                                  const size_t stride);
extern EbBool eb_image_compare_s8(const int8_t *buf1, const int8_t *buf2,
                                  const size_t width, const size_t height,
                                  const size_t stride);
extern EbBool eb_image_compare_u16(const uint16_t *buf1, const uint16_t *buf2,
                                   const size_t width, const size_t height,
                                   const size_t stride);
extern EbBool eb_image_compare_s16(const int16_t *buf1, const int16_t *buf2,
                                   const size_t width, const size_t height,
                                   const size_t stride);
extern EbBool eb_image_compare_u32(const uint32_t *buf1, const uint32_t *buf2,
                                   const size_t width, const size_t height,
                                   const size_t stride);
extern EbBool eb_image_compare_s32(const int32_t *buf1, const int32_t *buf2,
                                   const size_t width, const size_t height,
                                   const size_t stride);
extern EbBool eb_image_compare_u64(const uint64_t *buf1, const uint64_t *buf2,
                                   const size_t width, const size_t height,
                                   const size_t stride);

extern void eb_unit_test_log_s8(const char *const nameBuf,
                                const int8_t *const buf,
                                const uint32_t sizeBuf);
extern void eb_unit_test_log_u8(const char *const nameBuf,
                                const uint8_t *const buf,
                                const uint32_t sizeBuf);
extern void eb_unit_test_log_s16(const char *const nameBuf,
                                 const int16_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_u16(const char *const nameBuf,
                                 const uint16_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_s32(const char *const nameBuf,
                                 const int32_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_u32(const char *const nameBuf,
                                 const uint32_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_s64(const char *const nameBuf,
                                 const int64_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_u64(const char *const nameBuf,
                                 const uint64_t *const buf,
                                 const uint32_t sizeBuf);
extern void eb_unit_test_log_ptrdiff(const char *const nameBuf,
                                     const ptrdiff_t *const buf,
                                     const uint32_t sizeBuf);
extern void eb_unit_test_log_bool(const char *const nameBuf,
                                  const EbBool *const buf,
                                  const uint32_t sizeBuf);
extern void eb_unit_test_log_double(const char *const nameBuf,
                                    const double *const buf,
                                    const uint32_t sizeBuf);

extern void eb_unit_test_log_image_u8(const char *const nameBuf,
                                      const uint8_t *const buf,
                                      const uint32_t width,
                                      const uint32_t height,
                                      const ptrdiff_t stride);
extern void eb_unit_test_log_image_s16(const char *const nameBuf,
                                       const int16_t *const buf,
                                       const uint32_t width,
                                       const uint32_t height,
                                       const ptrdiff_t stride);
extern void eb_unit_test_log_image_u16(const char *const nameBuf,
                                       const uint16_t *const buf,
                                       const uint32_t width,
                                       const uint32_t height,
                                       const ptrdiff_t stride);
extern void eb_unit_test_log_image_s32(const char *const nameBuf,
                                       const int32_t *const buf,
                                       const uint32_t width,
                                       const uint32_t height,
                                       const ptrdiff_t stride);
extern void eb_unit_test_log_image_u32(const char *const nameBuf,
                                       const uint32_t *const buf,
                                       const uint32_t width,
                                       const uint32_t height,
                                       const ptrdiff_t stride);

#ifdef __cplusplus
}
#endif

#endif  // EbUnitTestUtility_h
