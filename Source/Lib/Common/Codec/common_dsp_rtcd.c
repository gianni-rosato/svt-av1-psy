/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#define RTCD_C
#include "common_dsp_rtcd.h"
#include "EbPictureOperators.h"
#include "EbPackUnPack_C.h"
#include "EbAvcStyleMcp.h"

/*
 * DSP deprecated flags
 */
#define HAS_MMX CPU_FLAGS_MMX
#define HAS_SSE CPU_FLAGS_SSE
#define HAS_SSE2 CPU_FLAGS_SSE2
#define HAS_SSE3 CPU_FLAGS_SSE3
#define HAS_SSSE3 CPU_FLAGS_SSSE3
#define HAS_SSE4_1 CPU_FLAGS_SSE4_1
#define HAS_SSE4_2 CPU_FLAGS_SSE4_2
#define HAS_AVX CPU_FLAGS_AVX
#define HAS_AVX2 CPU_FLAGS_AVX2
#define HAS_AVX512F CPU_FLAGS_AVX512F
#define HAS_AVX512CD CPU_FLAGS_AVX512CD
#define HAS_AVX512DQ CPU_FLAGS_AVX512DQ
#define HAS_AVX512ER CPU_FLAGS_AVX512ER
#define HAS_AVX512PF CPU_FLAGS_AVX512PF
#define HAS_AVX512BW CPU_FLAGS_AVX512BW
#define HAS_AVX512VL CPU_FLAGS_AVX512VL

// coeff: 16 bits, dynamic range [-32640, 32640].
// length: value range {16, 64, 256, 1024}.
int svt_aom_satd_c(const TranLow *coeff, int length) {
  int i;
  int satd = 0;
  for (i = 0; i < length; ++i) satd += abs(coeff[i]);

  // satd: 26 bits, dynamic range [-32640 * 1024, 32640 * 1024]
  return satd;
}

int64_t svt_av1_block_error_c(const TranLow *coeff, const TranLow *dqcoeff,
                          intptr_t block_size, int64_t *ssz) {
  int i;
  int64_t error = 0, sqcoeff = 0;

  for (i = 0; i < block_size; i++) {
    const int diff = coeff[i] - dqcoeff[i];
    error += diff * diff;
    sqcoeff += coeff[i] * coeff[i];
  }

  *ssz = sqcoeff;
  return error;
}

/**************************************
 * Instruction Set Support
 **************************************/
#ifdef ARCH_X86_64
// Helper Functions
static INLINE void run_cpuid(int eax, int ecx, int* abcd) {
#ifdef _WIN32
    __cpuidex(abcd, eax, ecx);
#else
    int ebx = 0, edx = 0;
#if defined(__i386__) && defined(__PIC__)
    /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__("movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi"
            : "=D"(ebx),
#else
    __asm__("cpuid"
    : "+b"(ebx),
#endif
    "+a"(eax),
    "+c"(ecx),
    "=d"(edx));
    abcd[0] = eax;
    abcd[1] = ebx;
    abcd[2] = ecx;
    abcd[3] = edx;
#endif
}

static INLINE int32_t check_xcr0_ymm() {
    uint32_t xcr0;
#ifdef _WIN32
    xcr0 = (uint32_t)_xgetbv(0); /* min VS2010 SP1 compiler is required */
#else
    __asm__("xgetbv" : "=a"(xcr0) : "c"(0) : "%edx");
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

static int32_t check_4thgen_intel_core_features() {
    int       abcd[4];
    const int fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    const int avx2_bmi12_mask        = (1 << 5) | (1 << 3) | (1 << 8);

    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   &&
    CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
    CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid(1, 0, abcd);
    if ((abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask) return 0;

    if (!check_xcr0_ymm()) return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
    CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
    CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
    run_cpuid(7, 0, abcd);
    if ((abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask) return 0;
    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    run_cpuid(0x80000001, 0, abcd);
    if ((abcd[2] & (1 << 5)) == 0) return 0;
    return 1;
}

static INLINE int check_xcr0_zmm() {
    uint32_t xcr0;
    uint32_t zmm_ymm_xmm = (7 << 5) | (1 << 2) | (1 << 1);
#ifdef _WIN32
    xcr0 = (uint32_t)_xgetbv(0); /* min VS2010 SP1 compiler is required */
#else
    __asm__("xgetbv" : "=a"(xcr0) : "c"(0) : "%edx");
#endif
    return ((xcr0 & zmm_ymm_xmm) ==
            zmm_ymm_xmm); /* check if xmm, ymm and zmm state are enabled in XCR0 */
}

static int32_t can_use_intel_avx512() {
    int abcd[4];

    /*  CPUID.(EAX=07H, ECX=0):EBX[bit 16]==1 AVX512F
    CPUID.(EAX=07H, ECX=0):EBX[bit 17] AVX512DQ
    CPUID.(EAX=07H, ECX=0):EBX[bit 28] AVX512CD
    CPUID.(EAX=07H, ECX=0):EBX[bit 30] AVX512BW
    CPUID.(EAX=07H, ECX=0):EBX[bit 31] AVX512VL */

    const int avx512_ebx_mask = (1u << 16) // AVX-512F
        | (1u << 17) // AVX-512DQ
        | (1u << 28) // AVX-512CD
        | (1u << 30) // AVX-512BW
        | (1u << 31); // AVX-512VL

    if (!check_4thgen_intel_core_features()) return 0;

    // ensure OS supports ZMM registers (and YMM, and XMM)
    if (!check_xcr0_zmm()) return 0;

    run_cpuid(7, 0, abcd);
    if ((abcd[1] & avx512_ebx_mask) != avx512_ebx_mask) return 0;

    return 1;
}

CPU_FLAGS get_cpu_flags() {
    CPU_FLAGS flags = 0;

    /* To detail tests CPU features, requires more accurate implementation.
        Documentation help:
        https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?redirectedfrom=MSDN&view=vs-2019
    */

    if (check_4thgen_intel_core_features()) {
        flags |= CPU_FLAGS_MMX | CPU_FLAGS_SSE | CPU_FLAGS_SSE2 | CPU_FLAGS_SSE3 | CPU_FLAGS_SSSE3 |
                 CPU_FLAGS_SSE4_1 | CPU_FLAGS_SSE4_2 | CPU_FLAGS_AVX | CPU_FLAGS_AVX2;
    }

    if (can_use_intel_avx512()) {
        flags |= CPU_FLAGS_AVX512F | CPU_FLAGS_AVX512DQ | CPU_FLAGS_AVX512CD | CPU_FLAGS_AVX512BW |
                 CPU_FLAGS_AVX512VL;
    }

    return flags;
}

CPU_FLAGS get_cpu_flags_to_use() {
    CPU_FLAGS flags = get_cpu_flags();
#ifdef NON_AVX512_SUPPORT
    /* Remove AVX512 flags. */
    flags &= (CPU_FLAGS_AVX512F - 1);
#endif
    return flags;
}
#endif

#ifdef ARCH_X86_64
#ifndef NON_AVX512_SUPPORT
#define SET_FUNCTIONS_AVX512(ptr, avx512)                                                         \
    if (((uintptr_t)NULL != (uintptr_t)avx512) && (flags & HAS_AVX512F)) ptr = avx512;
#else /* NON_AVX512_SUPPORT */
#define SET_FUNCTIONS_AVX512(ptr, avx512)
#endif /* NON_AVX512_SUPPORT */

#define SET_FUNCTIONS_X86(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    if (((uintptr_t)NULL != (uintptr_t)mmx)    && (flags & HAS_MMX))    ptr = mmx;                \
    if (((uintptr_t)NULL != (uintptr_t)sse)    && (flags & HAS_SSE))    ptr = sse;                \
    if (((uintptr_t)NULL != (uintptr_t)sse2)   && (flags & HAS_SSE2))   ptr = sse2;               \
    if (((uintptr_t)NULL != (uintptr_t)sse3)   && (flags & HAS_SSE3))   ptr = sse3;               \
    if (((uintptr_t)NULL != (uintptr_t)ssse3)  && (flags & HAS_SSSE3))  ptr = ssse3;              \
    if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1;             \
    if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2;             \
    if (((uintptr_t)NULL != (uintptr_t)avx)    && (flags & HAS_AVX))    ptr = avx;                \
    if (((uintptr_t)NULL != (uintptr_t)avx2)   && (flags & HAS_AVX2))   ptr = avx2;               \
    SET_FUNCTIONS_AVX512(ptr, avx512)
#else /* ARCH_X86_64 */
#define SET_FUNCTIONS_X86(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512)
#endif /* ARCH_X86_64 */

#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512)     \
    do {                                                                                          \
        if (ptr != 0) {                                                                           \
            printf("Error: %s:%i: Pointer \"%s\" is set before!\n", __FILE__, __LINE__, #ptr);    \
            assert(0);                                                                            \
        }                                                                                         \
        if ((uintptr_t)NULL == (uintptr_t)c) {                                                    \
            printf("Error: %s:%i: Pointer \"%s\" on C is NULL!\n", __FILE__, __LINE__, #ptr);     \
            assert(0);                                                                            \
        }                                                                                         \
        ptr = c;                                                                                  \
        SET_FUNCTIONS_X86(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    } while (0)

/* Macros SET_* use local variable CPU_FLAGS flags */
#define SET_ONLY_C(ptr, c)                                  SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#define SET_SSE2(ptr, c, sse2)                              SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, 0, 0)
#define SET_SSE2_AVX2(ptr, c, sse2, avx2)                   SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, avx2, 0)
#define SET_SSE2_AVX512(ptr, c, sse2, avx512)               SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, 0, avx512)
#define SET_SSSE3(ptr, c, ssse3)                            SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, ssse3, 0, 0, 0, 0, 0)
#define SET_SSSE3_AVX2(ptr, c, ssse3, avx2)                 SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, ssse3, 0, 0, 0, avx2, 0)
#define SET_SSE41(ptr, c, sse4_1)                           SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41_AVX2(ptr, c, sse4_1, avx2)                SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, 0)
#define SET_SSE41_AVX2_AVX512(ptr, c, sse4_1, avx2, avx512) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, avx512)
#define SET_AVX2(ptr, c, avx2)                              SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, 0)
#define SET_AVX2_AVX512(ptr, c, avx2, avx512)               SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, avx512)


void setup_common_rtcd_internal(CPU_FLAGS flags) {
#ifdef ARCH_X86_64
    /** Should be done during library initialization,
        but for safe limiting cpu flags again. */
    flags &= get_cpu_flags_to_use();
    //to use C: flags=0
#else
    (void)flags;
#endif

    SET_SSE41_AVX2(svt_aom_blend_a64_mask, svt_aom_blend_a64_mask_c, svt_aom_blend_a64_mask_sse4_1, svt_aom_blend_a64_mask_avx2);
    SET_SSE41(svt_aom_blend_a64_hmask, svt_aom_blend_a64_hmask_c, svt_aom_blend_a64_hmask_sse4_1);
    SET_SSE41(svt_aom_blend_a64_vmask, svt_aom_blend_a64_vmask_c, svt_aom_blend_a64_vmask_sse4_1);
    SET_AVX2(svt_aom_lowbd_blend_a64_d16_mask, svt_aom_lowbd_blend_a64_d16_mask_c, svt_aom_lowbd_blend_a64_d16_mask_avx2);
    SET_SSE41(svt_aom_highbd_blend_a64_mask, svt_aom_highbd_blend_a64_mask_c, svt_aom_highbd_blend_a64_mask_8bit_sse4_1);
    SET_SSE41(svt_aom_highbd_blend_a64_hmask_8bit, svt_aom_highbd_blend_a64_hmask_8bit_c, svt_aom_highbd_blend_a64_hmask_8bit_sse4_1);
    SET_SSE41(svt_aom_highbd_blend_a64_vmask_8bit, svt_aom_highbd_blend_a64_vmask_8bit_c, svt_aom_highbd_blend_a64_vmask_8bit_sse4_1);
    SET_SSE41(svt_aom_highbd_blend_a64_vmask_16bit, svt_aom_highbd_blend_a64_vmask_16bit_c, svt_aom_highbd_blend_a64_vmask_16bit_sse4_1);
    SET_SSE41(svt_aom_highbd_blend_a64_hmask_16bit, svt_aom_highbd_blend_a64_hmask_16bit_c, svt_aom_highbd_blend_a64_hmask_16bit_sse4_1);
    SET_AVX2(svt_aom_highbd_blend_a64_d16_mask, svt_aom_highbd_blend_a64_d16_mask_c, svt_aom_highbd_blend_a64_d16_mask_avx2);
    SET_AVX2(svt_cfl_predict_lbd, svt_cfl_predict_lbd_c, svt_cfl_predict_lbd_avx2);
    SET_AVX2(svt_cfl_predict_hbd, svt_cfl_predict_hbd_c, svt_cfl_predict_hbd_avx2);
    SET_SSE41(svt_av1_filter_intra_predictor, svt_av1_filter_intra_predictor_c, svt_av1_filter_intra_predictor_sse4_1);
    SET_SSE41(svt_av1_filter_intra_edge_high, svt_av1_filter_intra_edge_high_c, svt_av1_filter_intra_edge_high_sse4_1);
    SET_SSE41(svt_av1_filter_intra_edge, svt_av1_filter_intra_edge_c, svt_av1_filter_intra_edge_sse4_1);
    SET_SSE41(svt_av1_upsample_intra_edge, svt_av1_upsample_intra_edge_c, svt_av1_upsample_intra_edge_sse4_1);
    SET_AVX2(svt_av1_build_compound_diffwtd_mask_d16, svt_av1_build_compound_diffwtd_mask_d16_c, svt_av1_build_compound_diffwtd_mask_d16_avx2);
    SET_AVX2(svt_av1_highbd_wiener_convolve_add_src, svt_av1_highbd_wiener_convolve_add_src_c, svt_av1_highbd_wiener_convolve_add_src_avx2);
    SET_AVX2(svt_apply_selfguided_restoration, svt_apply_selfguided_restoration_c, svt_apply_selfguided_restoration_avx2);
    SET_AVX2(svt_av1_selfguided_restoration, svt_av1_selfguided_restoration_c, svt_av1_selfguided_restoration_avx2);
    SET_AVX2(svt_av1_inv_txfm2d_add_4x4, svt_av1_inv_txfm2d_add_4x4_c, svt_av1_inv_txfm2d_add_4x4_avx2);
    SET_SSE41(svt_av1_inv_txfm2d_add_4x8, svt_av1_inv_txfm2d_add_4x8_c, svt_av1_inv_txfm2d_add_4x8_sse4_1);
    SET_SSE41(svt_av1_inv_txfm2d_add_4x16, svt_av1_inv_txfm2d_add_4x16_c, svt_av1_inv_txfm2d_add_4x16_sse4_1);
    SET_SSE41(svt_av1_inv_txfm2d_add_8x4, svt_av1_inv_txfm2d_add_8x4_c, svt_av1_inv_txfm2d_add_8x4_sse4_1);
    SET_AVX2(svt_av1_inv_txfm2d_add_8x8, svt_av1_inv_txfm2d_add_8x8_c, svt_av1_inv_txfm2d_add_8x8_avx2);
    SET_AVX2(svt_av1_inv_txfm2d_add_8x16, svt_av1_inv_txfm2d_add_8x16_c, svt_av1_highbd_inv_txfm_add_avx2);
    SET_AVX2(svt_av1_inv_txfm2d_add_8x32, svt_av1_inv_txfm2d_add_8x32_c, svt_av1_highbd_inv_txfm_add_avx2);
    SET_SSE41(svt_av1_inv_txfm2d_add_16x4, svt_av1_inv_txfm2d_add_16x4_c, svt_av1_inv_txfm2d_add_16x4_sse4_1);
    SET_AVX2(svt_av1_inv_txfm2d_add_16x8, svt_av1_inv_txfm2d_add_16x8_c, svt_av1_highbd_inv_txfm_add_avx2);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_16x16, svt_av1_inv_txfm2d_add_16x16_c, svt_av1_inv_txfm2d_add_16x16_avx2, svt_av1_inv_txfm2d_add_16x16_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_16x32, svt_av1_inv_txfm2d_add_16x32_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_16x32_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_16x64, svt_av1_inv_txfm2d_add_16x64_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_16x64_avx512);
    SET_AVX2(svt_av1_inv_txfm2d_add_32x8, svt_av1_inv_txfm2d_add_32x8_c, svt_av1_highbd_inv_txfm_add_avx2);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_32x16, svt_av1_inv_txfm2d_add_32x16_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_32x16_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_32x32, svt_av1_inv_txfm2d_add_32x32_c, svt_av1_inv_txfm2d_add_32x32_avx2, svt_av1_inv_txfm2d_add_32x32_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_32x64, svt_av1_inv_txfm2d_add_32x64_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_32x64_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_64x16, svt_av1_inv_txfm2d_add_64x16_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_64x16_avx512);
    SET_AVX2_AVX512(svt_av1_inv_txfm2d_add_64x32, svt_av1_inv_txfm2d_add_64x32_c, svt_av1_highbd_inv_txfm_add_avx2, svt_av1_inv_txfm2d_add_64x32_avx512);
    SET_SSE41_AVX2_AVX512(svt_av1_inv_txfm2d_add_64x64, svt_av1_inv_txfm2d_add_64x64_c, svt_av1_inv_txfm2d_add_64x64_sse4_1, svt_av1_inv_txfm2d_add_64x64_avx2, svt_av1_inv_txfm2d_add_64x64_avx512);
    SET_SSSE3_AVX2(svt_av1_inv_txfm_add, svt_av1_inv_txfm_add_c, svt_av1_inv_txfm_add_ssse3, svt_av1_inv_txfm_add_avx2);
    SET_AVX2(svt_compressed_packmsb, svt_compressed_packmsb_c, svt_compressed_packmsb_avx2_intrin);
    SET_AVX2(svt_c_pack, svt_c_pack_c, svt_c_pack_avx2_intrin);
    SET_SSE2_AVX2(svt_unpack_avg, svt_unpack_avg_c, svt_unpack_avg_sse2_intrin, svt_unpack_avg_avx2_intrin);
    SET_AVX2(svt_unpack_avg_safe_sub, svt_unpack_avg_safe_sub_c, svt_unpack_avg_safe_sub_avx2_intrin);
    SET_AVX2(svt_un_pack8_bit_data, svt_un_pack8_bit_data_c, svt_enc_un_pack8_bit_data_avx2_intrin);
    SET_AVX2(svt_cfl_luma_subsampling_420_lbd, svt_cfl_luma_subsampling_420_lbd_c, svt_cfl_luma_subsampling_420_lbd_avx2);
    SET_AVX2(svt_cfl_luma_subsampling_420_hbd, svt_cfl_luma_subsampling_420_hbd_c, svt_cfl_luma_subsampling_420_hbd_avx2);
    SET_AVX2(svt_convert_8bit_to_16bit, svt_convert_8bit_to_16bit_c, svt_convert_8bit_to_16bit_avx2);
    SET_AVX2(svt_convert_16bit_to_8bit, svt_convert_16bit_to_8bit_c, svt_convert_16bit_to_8bit_avx2);
    SET_SSE2_AVX2(svt_pack2d_16_bit_src_mul4, svt_enc_msb_pack2_d, svt_enc_msb_pack2d_sse2_intrin, svt_enc_msb_pack2d_avx2_intrin_al);
    SET_SSE2(svt_un_pack2d_16_bit_src_mul4, svt_enc_msb_un_pack2_d, svt_enc_msb_un_pack2d_sse2_intrin);
    SET_AVX2(svt_full_distortion_kernel_cbf_zero32_bits, svt_full_distortion_kernel_cbf_zero32_bits_c, svt_full_distortion_kernel_cbf_zero32_bits_avx2);
    SET_AVX2(svt_full_distortion_kernel32_bits, svt_full_distortion_kernel32_bits_c, svt_full_distortion_kernel32_bits_avx2);

    SET_AVX2_AVX512(svt_spatial_full_distortion_kernel, svt_spatial_full_distortion_kernel_c, svt_spatial_full_distortion_kernel_avx2, svt_spatial_full_distortion_kernel_avx512);
    SET_AVX2(svt_full_distortion_kernel16_bits, svt_full_distortion_kernel16_bits_c, svt_full_distortion_kernel16_bits_avx2);
    SET_AVX2_AVX512(svt_residual_kernel8bit, svt_residual_kernel8bit_c, svt_residual_kernel8bit_avx2, svt_residual_kernel8bit_avx512);
    SET_SSE2_AVX2(svt_residual_kernel16bit, svt_residual_kernel16bit_c, svt_residual_kernel16bit_sse2_intrin, svt_residual_kernel16bit_avx2);
    SET_SSE2(svt_picture_average_kernel, svt_picture_average_kernel_c, svt_picture_average_kernel_sse2_intrin);
    SET_SSE2(svt_picture_average_kernel1_line, svt_picture_average_kernel1_line_c, svt_picture_average_kernel1_line_sse2_intrin);
    SET_AVX2_AVX512(svt_av1_wiener_convolve_add_src, svt_av1_wiener_convolve_add_src_c, svt_av1_wiener_convolve_add_src_avx2, svt_av1_wiener_convolve_add_src_avx512);
    SET_ONLY_C(svt_av1_convolve_2d_scale, svt_av1_convolve_2d_scale_c);
    SET_AVX2(svt_av1_highbd_convolve_y_sr, svt_av1_highbd_convolve_y_sr_c, svt_av1_highbd_convolve_y_sr_avx2);
    SET_AVX2(svt_av1_highbd_convolve_2d_sr, svt_av1_highbd_convolve_2d_sr_c, svt_av1_highbd_convolve_2d_sr_avx2);
    SET_ONLY_C(svt_av1_highbd_convolve_2d_scale, svt_av1_highbd_convolve_2d_scale_c);
    SET_AVX2(svt_av1_highbd_convolve_2d_copy_sr, svt_av1_highbd_convolve_2d_copy_sr_c, svt_av1_highbd_convolve_2d_copy_sr_avx2);
    SET_AVX2(svt_av1_highbd_jnt_convolve_2d, svt_av1_highbd_jnt_convolve_2d_c, svt_av1_highbd_jnt_convolve_2d_avx2);
    SET_AVX2(svt_av1_highbd_jnt_convolve_2d_copy, svt_av1_highbd_jnt_convolve_2d_copy_c, svt_av1_highbd_jnt_convolve_2d_copy_avx2);
    SET_AVX2(svt_av1_highbd_jnt_convolve_x, svt_av1_highbd_jnt_convolve_x_c, svt_av1_highbd_jnt_convolve_x_avx2);
    SET_AVX2(svt_av1_highbd_jnt_convolve_y, svt_av1_highbd_jnt_convolve_y_c, svt_av1_highbd_jnt_convolve_y_avx2);
    SET_AVX2(svt_av1_highbd_convolve_x_sr, svt_av1_highbd_convolve_x_sr_c, svt_av1_highbd_convolve_x_sr_avx2);
    SET_AVX2_AVX512(svt_av1_convolve_2d_sr, svt_av1_convolve_2d_sr_c, svt_av1_convolve_2d_sr_avx2, svt_av1_convolve_2d_sr_avx512);
    SET_AVX2_AVX512(svt_av1_convolve_2d_copy_sr, svt_av1_convolve_2d_copy_sr_c, svt_av1_convolve_2d_copy_sr_avx2, svt_av1_convolve_2d_copy_sr_avx512);
    SET_AVX2_AVX512(svt_av1_convolve_x_sr, svt_av1_convolve_x_sr_c, svt_av1_convolve_x_sr_avx2, svt_av1_convolve_x_sr_avx512);
    SET_AVX2_AVX512(svt_av1_convolve_y_sr, svt_av1_convolve_y_sr_c, svt_av1_convolve_y_sr_avx2, svt_av1_convolve_y_sr_avx512);
    SET_AVX2_AVX512(svt_av1_jnt_convolve_2d, svt_av1_jnt_convolve_2d_c, svt_av1_jnt_convolve_2d_avx2, svt_av1_jnt_convolve_2d_avx512);
    SET_AVX2_AVX512(svt_av1_jnt_convolve_2d_copy, svt_av1_jnt_convolve_2d_copy_c, svt_av1_jnt_convolve_2d_copy_avx2, svt_av1_jnt_convolve_2d_copy_avx512);
    SET_AVX2_AVX512(svt_av1_jnt_convolve_x, svt_av1_jnt_convolve_x_c, svt_av1_jnt_convolve_x_avx2, svt_av1_jnt_convolve_x_avx512);
    SET_AVX2_AVX512(svt_av1_jnt_convolve_y, svt_av1_jnt_convolve_y_c, svt_av1_jnt_convolve_y_avx2, svt_av1_jnt_convolve_y_avx512);
    SET_AVX2(svt_aom_convolve8_horiz, svt_aom_convolve8_horiz_c, svt_aom_convolve8_horiz_avx2);
    SET_AVX2(svt_aom_convolve8_vert, svt_aom_convolve8_vert_c, svt_aom_convolve8_vert_avx2);
    SET_AVX2(svt_av1_build_compound_diffwtd_mask, svt_av1_build_compound_diffwtd_mask_c, svt_av1_build_compound_diffwtd_mask_avx2);
    SET_AVX2(svt_av1_build_compound_diffwtd_mask_highbd, svt_av1_build_compound_diffwtd_mask_highbd_c, svt_av1_build_compound_diffwtd_mask_highbd_avx2);
    SET_AVX2(svt_av1_wedge_sse_from_residuals, svt_av1_wedge_sse_from_residuals_c, svt_av1_wedge_sse_from_residuals_avx2);
    SET_AVX2(svt_aom_subtract_block, svt_aom_subtract_block_c, svt_aom_subtract_block_avx2);
    SET_SSE2(svt_aom_highbd_subtract_block, svt_aom_highbd_subtract_block_c, svt_aom_highbd_subtract_block_sse2);
    SET_ONLY_C(svt_aom_highbd_smooth_v_predictor_2x2, svt_aom_highbd_smooth_v_predictor_2x2_c);
    SET_SSSE3(svt_aom_highbd_smooth_v_predictor_4x4, svt_aom_highbd_smooth_v_predictor_4x4_c, svt_aom_highbd_smooth_v_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_v_predictor_4x8, svt_aom_highbd_smooth_v_predictor_4x8_c, svt_aom_highbd_smooth_v_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_v_predictor_4x16, svt_aom_highbd_smooth_v_predictor_4x16_c, svt_aom_highbd_smooth_v_predictor_4x16_ssse3);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_8x4, svt_aom_highbd_smooth_v_predictor_8x4_c, svt_aom_highbd_smooth_v_predictor_8x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_8x8, svt_aom_highbd_smooth_v_predictor_8x8_c, svt_aom_highbd_smooth_v_predictor_8x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_8x16, svt_aom_highbd_smooth_v_predictor_8x16_c, svt_aom_highbd_smooth_v_predictor_8x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_8x32, svt_aom_highbd_smooth_v_predictor_8x32_c, svt_aom_highbd_smooth_v_predictor_8x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_16x4, svt_aom_highbd_smooth_v_predictor_16x4_c, svt_aom_highbd_smooth_v_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_16x8, svt_aom_highbd_smooth_v_predictor_16x8_c, svt_aom_highbd_smooth_v_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_16x16, svt_aom_highbd_smooth_v_predictor_16x16_c, svt_aom_highbd_smooth_v_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_16x32, svt_aom_highbd_smooth_v_predictor_16x32_c, svt_aom_highbd_smooth_v_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_v_predictor_16x64, svt_aom_highbd_smooth_v_predictor_16x64_c, svt_aom_highbd_smooth_v_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_32x8, svt_aom_highbd_smooth_v_predictor_32x8_c, svt_aom_highbd_smooth_v_predictor_32x8_avx2, aom_highbd_smooth_v_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_32x16, svt_aom_highbd_smooth_v_predictor_32x16_c, svt_aom_highbd_smooth_v_predictor_32x16_avx2, aom_highbd_smooth_v_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_32x32, svt_aom_highbd_smooth_v_predictor_32x32_c, svt_aom_highbd_smooth_v_predictor_32x32_avx2, aom_highbd_smooth_v_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_32x64, svt_aom_highbd_smooth_v_predictor_32x64_c, svt_aom_highbd_smooth_v_predictor_32x64_avx2, aom_highbd_smooth_v_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_64x16, svt_aom_highbd_smooth_v_predictor_64x16_c, svt_aom_highbd_smooth_v_predictor_64x16_avx2, aom_highbd_smooth_v_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_64x32, svt_aom_highbd_smooth_v_predictor_64x32_c, svt_aom_highbd_smooth_v_predictor_64x32_avx2, aom_highbd_smooth_v_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_v_predictor_64x64, svt_aom_highbd_smooth_v_predictor_64x64_c, svt_aom_highbd_smooth_v_predictor_64x64_avx2, aom_highbd_smooth_v_predictor_64x64_avx512);
    SET_AVX2(svt_av1_dr_prediction_z1, svt_av1_dr_prediction_z1_c, svt_av1_dr_prediction_z1_avx2);
    SET_AVX2(svt_av1_dr_prediction_z2, svt_av1_dr_prediction_z2_c, svt_av1_dr_prediction_z2_avx2);
    SET_AVX2(svt_av1_dr_prediction_z3, svt_av1_dr_prediction_z3_c, svt_av1_dr_prediction_z3_avx2);
    SET_AVX2(svt_av1_highbd_dr_prediction_z1, svt_av1_highbd_dr_prediction_z1_c, svt_av1_highbd_dr_prediction_z1_avx2);
    SET_AVX2(svt_av1_highbd_dr_prediction_z2, svt_av1_highbd_dr_prediction_z2_c, svt_av1_highbd_dr_prediction_z2_avx2);
    SET_AVX2(svt_av1_highbd_dr_prediction_z3, svt_av1_highbd_dr_prediction_z3_c, svt_av1_highbd_dr_prediction_z3_avx2);
    SET_SSSE3(svt_aom_paeth_predictor_4x4, svt_aom_paeth_predictor_4x4_c, svt_aom_paeth_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_4x8, svt_aom_paeth_predictor_4x8_c, svt_aom_paeth_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_4x16, svt_aom_paeth_predictor_4x16_c, svt_aom_paeth_predictor_4x16_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_8x4, svt_aom_paeth_predictor_8x4_c, svt_aom_paeth_predictor_8x4_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_8x8, svt_aom_paeth_predictor_8x8_c, svt_aom_paeth_predictor_8x8_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_8x16, svt_aom_paeth_predictor_8x16_c, svt_aom_paeth_predictor_8x16_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_8x32, svt_aom_paeth_predictor_8x32_c, svt_aom_paeth_predictor_8x32_ssse3);
    SET_SSSE3(svt_aom_paeth_predictor_16x4, svt_aom_paeth_predictor_16x4_c, svt_aom_paeth_predictor_16x4_ssse3);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_16x8, svt_aom_paeth_predictor_16x8_c, svt_aom_paeth_predictor_16x8_ssse3, svt_aom_paeth_predictor_16x8_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_16x16, svt_aom_paeth_predictor_16x16_c, svt_aom_paeth_predictor_16x16_ssse3, svt_aom_paeth_predictor_16x16_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_16x32, svt_aom_paeth_predictor_16x32_c, svt_aom_paeth_predictor_16x32_ssse3, svt_aom_paeth_predictor_16x32_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_16x64, svt_aom_paeth_predictor_16x64_c, svt_aom_paeth_predictor_16x64_ssse3, svt_aom_paeth_predictor_16x64_avx2);
    SET_SSSE3(svt_aom_paeth_predictor_32x8, svt_aom_paeth_predictor_32x8_c, svt_aom_paeth_predictor_32x8_ssse3);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_32x16, svt_aom_paeth_predictor_32x16_c, svt_aom_paeth_predictor_32x16_ssse3, svt_aom_paeth_predictor_32x16_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_32x32, svt_aom_paeth_predictor_32x32_c, svt_aom_paeth_predictor_32x32_ssse3, svt_aom_paeth_predictor_32x32_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_32x64, svt_aom_paeth_predictor_32x64_c, svt_aom_paeth_predictor_32x64_ssse3, svt_aom_paeth_predictor_32x64_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_64x16, svt_aom_paeth_predictor_64x16_c, svt_aom_paeth_predictor_64x16_ssse3, svt_aom_paeth_predictor_64x16_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_64x32, svt_aom_paeth_predictor_64x32_c, svt_aom_paeth_predictor_64x32_ssse3, svt_aom_paeth_predictor_64x32_avx2);
    SET_SSSE3_AVX2(svt_aom_paeth_predictor_64x64, svt_aom_paeth_predictor_64x64_c, svt_aom_paeth_predictor_64x64_ssse3, svt_aom_paeth_predictor_64x64_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_2x2, svt_aom_highbd_paeth_predictor_2x2_c, svt_aom_highbd_paeth_predictor_2x2_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_4x4, svt_aom_highbd_paeth_predictor_4x4_c, svt_aom_highbd_paeth_predictor_4x4_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_4x8, svt_aom_highbd_paeth_predictor_4x8_c, svt_aom_highbd_paeth_predictor_4x8_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_4x16, svt_aom_highbd_paeth_predictor_4x16_c, svt_aom_highbd_paeth_predictor_4x16_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_8x4, svt_aom_highbd_paeth_predictor_8x4_c, svt_aom_highbd_paeth_predictor_8x4_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_8x8, svt_aom_highbd_paeth_predictor_8x8_c, svt_aom_highbd_paeth_predictor_8x8_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_8x16, svt_aom_highbd_paeth_predictor_8x16_c, svt_aom_highbd_paeth_predictor_8x16_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_8x32, svt_aom_highbd_paeth_predictor_8x32_c, svt_aom_highbd_paeth_predictor_8x32_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_16x4, svt_aom_highbd_paeth_predictor_16x4_c, svt_aom_highbd_paeth_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_16x8, svt_aom_highbd_paeth_predictor_16x8_c, svt_aom_highbd_paeth_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_16x16, svt_aom_highbd_paeth_predictor_16x16_c, svt_aom_highbd_paeth_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_16x32, svt_aom_highbd_paeth_predictor_16x32_c, svt_aom_highbd_paeth_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_16x64, svt_aom_highbd_paeth_predictor_16x64_c, svt_aom_highbd_paeth_predictor_16x64_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_32x8, svt_aom_highbd_paeth_predictor_32x8_c, svt_aom_highbd_paeth_predictor_32x8_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_32x16, svt_aom_highbd_paeth_predictor_32x16_c, svt_aom_highbd_paeth_predictor_32x16_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_32x32, svt_aom_highbd_paeth_predictor_32x32_c, svt_aom_highbd_paeth_predictor_32x32_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_32x64, svt_aom_highbd_paeth_predictor_32x64_c, svt_aom_highbd_paeth_predictor_32x64_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_64x16, svt_aom_highbd_paeth_predictor_64x16_c, svt_aom_highbd_paeth_predictor_64x16_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_64x32, svt_aom_highbd_paeth_predictor_64x32_c, svt_aom_highbd_paeth_predictor_64x32_avx2);
    SET_AVX2(svt_aom_highbd_paeth_predictor_64x64, svt_aom_highbd_paeth_predictor_64x64_c, svt_aom_highbd_paeth_predictor_64x64_avx2);
    SET_SSE2(aom_sum_squares_i16, svt_aom_sum_squares_i16_c, svt_aom_sum_squares_i16_sse2);
    SET_SSE2(svt_aom_dc_predictor_4x4, svt_aom_dc_predictor_4x4_c, svt_aom_dc_predictor_4x4_sse2);
    SET_SSE2(svt_aom_dc_predictor_4x8, svt_aom_dc_predictor_4x8_c, svt_aom_dc_predictor_4x8_sse2);
    SET_SSE2(svt_aom_dc_predictor_4x16, svt_aom_dc_predictor_4x16_c, svt_aom_dc_predictor_4x16_sse2);
    SET_SSE2(svt_aom_dc_predictor_8x4, svt_aom_dc_predictor_8x4_c, svt_aom_dc_predictor_8x4_sse2);
    SET_SSE2(svt_aom_dc_predictor_8x8, svt_aom_dc_predictor_8x8_c, svt_aom_dc_predictor_8x8_sse2);
    SET_SSE2(svt_aom_dc_predictor_8x16, svt_aom_dc_predictor_8x16_c, svt_aom_dc_predictor_8x16_sse2);
    SET_SSE2(svt_aom_dc_predictor_8x32, svt_aom_dc_predictor_8x32_c, svt_aom_dc_predictor_8x32_sse2);
    SET_SSE2(svt_aom_dc_predictor_16x4, svt_aom_dc_predictor_16x4_c, svt_aom_dc_predictor_16x4_sse2);
    SET_SSE2(svt_aom_dc_predictor_16x8, svt_aom_dc_predictor_16x8_c, svt_aom_dc_predictor_16x8_sse2);
    SET_SSE2(svt_aom_dc_predictor_16x16, svt_aom_dc_predictor_16x16_c, svt_aom_dc_predictor_16x16_sse2);
    SET_SSE2(svt_aom_dc_predictor_16x32, svt_aom_dc_predictor_16x32_c, svt_aom_dc_predictor_16x32_sse2);
    SET_SSE2(svt_aom_dc_predictor_16x64, svt_aom_dc_predictor_16x64_c, svt_aom_dc_predictor_16x64_sse2);
    SET_SSE2(svt_aom_dc_predictor_32x8, svt_aom_dc_predictor_32x8_c, svt_aom_dc_predictor_32x8_sse2);
    SET_AVX2(svt_aom_dc_predictor_32x16, svt_aom_dc_predictor_32x16_c, svt_aom_dc_predictor_32x16_avx2);
    SET_AVX2(svt_aom_dc_predictor_32x32, svt_aom_dc_predictor_32x32_c, svt_aom_dc_predictor_32x32_avx2);
    SET_AVX2(svt_aom_dc_predictor_32x64, svt_aom_dc_predictor_32x64_c, svt_aom_dc_predictor_32x64_avx2);
    SET_AVX2(svt_aom_dc_predictor_64x16, svt_aom_dc_predictor_64x16_c, svt_aom_dc_predictor_64x16_avx2);
    SET_AVX2(svt_aom_dc_predictor_64x32, svt_aom_dc_predictor_64x32_c, svt_aom_dc_predictor_64x32_avx2);
    SET_AVX2(svt_aom_dc_predictor_64x64, svt_aom_dc_predictor_64x64_c, svt_aom_dc_predictor_64x64_avx2);

    SET_SSE2(svt_aom_dc_top_predictor_4x4, svt_aom_dc_top_predictor_4x4_c, svt_aom_dc_top_predictor_4x4_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_4x8, svt_aom_dc_top_predictor_4x8_c, svt_aom_dc_top_predictor_4x8_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_4x16, svt_aom_dc_top_predictor_4x16_c, svt_aom_dc_top_predictor_4x16_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_8x4, svt_aom_dc_top_predictor_8x4_c, svt_aom_dc_top_predictor_8x4_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_8x8, svt_aom_dc_top_predictor_8x8_c, svt_aom_dc_top_predictor_8x8_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_8x16, svt_aom_dc_top_predictor_8x16_c, svt_aom_dc_top_predictor_8x16_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_8x32, svt_aom_dc_top_predictor_8x32_c, svt_aom_dc_top_predictor_8x32_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_16x4, svt_aom_dc_top_predictor_16x4_c, svt_aom_dc_top_predictor_16x4_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_16x8, svt_aom_dc_top_predictor_16x8_c, svt_aom_dc_top_predictor_16x8_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_16x16, svt_aom_dc_top_predictor_16x16_c, svt_aom_dc_top_predictor_16x16_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_16x32, svt_aom_dc_top_predictor_16x32_c, svt_aom_dc_top_predictor_16x32_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_16x64, svt_aom_dc_top_predictor_16x64_c, svt_aom_dc_top_predictor_16x64_sse2);
    SET_SSE2(svt_aom_dc_top_predictor_32x8, svt_aom_dc_top_predictor_32x8_c, svt_aom_dc_top_predictor_32x8_sse2);
    SET_AVX2(svt_aom_dc_top_predictor_32x16, svt_aom_dc_top_predictor_32x16_c, svt_aom_dc_top_predictor_32x16_avx2);
    SET_AVX2(svt_aom_dc_top_predictor_32x32, svt_aom_dc_top_predictor_32x32_c, svt_aom_dc_top_predictor_32x32_avx2);
    SET_AVX2(svt_aom_dc_top_predictor_32x64, svt_aom_dc_top_predictor_32x64_c, svt_aom_dc_top_predictor_32x64_avx2);
    SET_AVX2(svt_aom_dc_top_predictor_64x16, svt_aom_dc_top_predictor_64x16_c, svt_aom_dc_top_predictor_64x16_avx2);
    SET_AVX2(svt_aom_dc_top_predictor_64x32, svt_aom_dc_top_predictor_64x32_c, svt_aom_dc_top_predictor_64x32_avx2);
    SET_AVX2(svt_aom_dc_top_predictor_64x64, svt_aom_dc_top_predictor_64x64_c, svt_aom_dc_top_predictor_64x64_avx2);

    SET_SSE2(svt_aom_dc_left_predictor_4x4, svt_aom_dc_left_predictor_4x4_c, svt_aom_dc_left_predictor_4x4_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_4x8, svt_aom_dc_left_predictor_4x8_c, svt_aom_dc_left_predictor_4x8_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_4x16, svt_aom_dc_left_predictor_4x16_c, svt_aom_dc_left_predictor_4x16_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_8x4, svt_aom_dc_left_predictor_8x4_c, svt_aom_dc_left_predictor_8x4_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_8x8, svt_aom_dc_left_predictor_8x8_c, svt_aom_dc_left_predictor_8x8_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_8x16, svt_aom_dc_left_predictor_8x16_c, svt_aom_dc_left_predictor_8x16_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_8x32, svt_aom_dc_left_predictor_8x32_c, svt_aom_dc_left_predictor_8x32_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_16x4, svt_aom_dc_left_predictor_16x4_c, svt_aom_dc_left_predictor_16x4_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_16x8, svt_aom_dc_left_predictor_16x8_c, svt_aom_dc_left_predictor_16x8_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_16x16, svt_aom_dc_left_predictor_16x16_c, svt_aom_dc_left_predictor_16x16_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_16x32, svt_aom_dc_left_predictor_16x32_c, svt_aom_dc_left_predictor_16x32_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_16x64, svt_aom_dc_left_predictor_16x64_c, svt_aom_dc_left_predictor_16x64_sse2);
    SET_SSE2(svt_aom_dc_left_predictor_32x8, svt_aom_dc_left_predictor_32x8_c, svt_aom_dc_left_predictor_32x8_sse2);
    SET_AVX2(svt_aom_dc_left_predictor_32x16, svt_aom_dc_left_predictor_32x16_c, svt_aom_dc_left_predictor_32x16_avx2);
    SET_AVX2(svt_aom_dc_left_predictor_32x32, svt_aom_dc_left_predictor_32x32_c, svt_aom_dc_left_predictor_32x32_avx2);
    SET_AVX2(svt_aom_dc_left_predictor_32x64, svt_aom_dc_left_predictor_32x64_c, svt_aom_dc_left_predictor_32x64_avx2);
    SET_AVX2(svt_aom_dc_left_predictor_64x16, svt_aom_dc_left_predictor_64x16_c, svt_aom_dc_left_predictor_64x16_avx2);
    SET_AVX2(svt_aom_dc_left_predictor_64x32, svt_aom_dc_left_predictor_64x32_c, svt_aom_dc_left_predictor_64x32_avx2);
    SET_AVX2(svt_aom_dc_left_predictor_64x64, svt_aom_dc_left_predictor_64x64_c, svt_aom_dc_left_predictor_64x64_avx2);

    SET_SSE2(svt_aom_dc_128_predictor_4x4, svt_aom_dc_128_predictor_4x4_c, svt_aom_dc_128_predictor_4x4_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_4x8, svt_aom_dc_128_predictor_4x8_c, svt_aom_dc_128_predictor_4x8_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_4x16, svt_aom_dc_128_predictor_4x16_c, svt_aom_dc_128_predictor_4x16_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_8x4, svt_aom_dc_128_predictor_8x4_c, svt_aom_dc_128_predictor_8x4_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_8x8, svt_aom_dc_128_predictor_8x8_c, svt_aom_dc_128_predictor_8x8_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_8x16, svt_aom_dc_128_predictor_8x16_c, svt_aom_dc_128_predictor_8x16_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_8x32, svt_aom_dc_128_predictor_8x32_c, svt_aom_dc_128_predictor_8x32_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_16x4, svt_aom_dc_128_predictor_16x4_c, svt_aom_dc_128_predictor_16x4_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_16x8, svt_aom_dc_128_predictor_16x8_c, svt_aom_dc_128_predictor_16x8_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_16x16, svt_aom_dc_128_predictor_16x16_c, svt_aom_dc_128_predictor_16x16_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_16x32, svt_aom_dc_128_predictor_16x32_c, svt_aom_dc_128_predictor_16x32_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_16x64, svt_aom_dc_128_predictor_16x64_c, svt_aom_dc_128_predictor_16x64_sse2);
    SET_SSE2(svt_aom_dc_128_predictor_32x8, svt_aom_dc_128_predictor_32x8_c, svt_aom_dc_128_predictor_32x8_sse2);
    SET_AVX2(svt_aom_dc_128_predictor_32x16, svt_aom_dc_128_predictor_32x16_c, svt_aom_dc_128_predictor_32x16_avx2);
    SET_AVX2(svt_aom_dc_128_predictor_32x32, svt_aom_dc_128_predictor_32x32_c, svt_aom_dc_128_predictor_32x32_avx2);
    SET_AVX2(svt_aom_dc_128_predictor_32x64, svt_aom_dc_128_predictor_32x64_c, svt_aom_dc_128_predictor_32x64_avx2);
    SET_AVX2(svt_aom_dc_128_predictor_64x16, svt_aom_dc_128_predictor_64x16_c, svt_aom_dc_128_predictor_64x16_avx2);
    SET_AVX2(svt_aom_dc_128_predictor_64x32, svt_aom_dc_128_predictor_64x32_c, svt_aom_dc_128_predictor_64x32_avx2);
    SET_AVX2(svt_aom_dc_128_predictor_64x64, svt_aom_dc_128_predictor_64x64_c, svt_aom_dc_128_predictor_64x64_avx2);

    SET_SSSE3(svt_aom_smooth_h_predictor_4x4, svt_aom_smooth_h_predictor_4x4_c, svt_aom_smooth_h_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_4x8, svt_aom_smooth_h_predictor_4x8_c, svt_aom_smooth_h_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_4x16, svt_aom_smooth_h_predictor_4x16_c, svt_aom_smooth_h_predictor_4x16_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_8x4, svt_aom_smooth_h_predictor_8x4_c, svt_aom_smooth_h_predictor_8x4_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_8x8, svt_aom_smooth_h_predictor_8x8_c, svt_aom_smooth_h_predictor_8x8_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_8x16, svt_aom_smooth_h_predictor_8x16_c, svt_aom_smooth_h_predictor_8x16_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_8x32, svt_aom_smooth_h_predictor_8x32_c, svt_aom_smooth_h_predictor_8x32_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_16x4, svt_aom_smooth_h_predictor_16x4_c, svt_aom_smooth_h_predictor_16x4_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_16x8, svt_aom_smooth_h_predictor_16x8_c, svt_aom_smooth_h_predictor_16x8_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_16x16, svt_aom_smooth_h_predictor_16x16_c, svt_aom_smooth_h_predictor_16x16_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_16x32, svt_aom_smooth_h_predictor_16x32_c, svt_aom_smooth_h_predictor_16x32_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_16x64, svt_aom_smooth_h_predictor_16x64_c, svt_aom_smooth_h_predictor_16x64_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_32x8, svt_aom_smooth_h_predictor_32x8_c, svt_aom_smooth_h_predictor_32x8_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_32x16, svt_aom_smooth_h_predictor_32x16_c, svt_aom_smooth_h_predictor_32x16_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_32x32, svt_aom_smooth_h_predictor_32x32_c, svt_aom_smooth_h_predictor_32x32_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_32x64, svt_aom_smooth_h_predictor_32x64_c, svt_aom_smooth_h_predictor_32x64_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_64x16, svt_aom_smooth_h_predictor_64x16_c, svt_aom_smooth_h_predictor_64x16_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_64x32, svt_aom_smooth_h_predictor_64x32_c, svt_aom_smooth_h_predictor_64x32_ssse3);
    SET_SSSE3(svt_aom_smooth_h_predictor_64x64, svt_aom_smooth_h_predictor_64x64_c, svt_aom_smooth_h_predictor_64x64_ssse3);

    SET_SSSE3(svt_aom_smooth_v_predictor_4x4, svt_aom_smooth_v_predictor_4x4_c, svt_aom_smooth_v_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_4x8, svt_aom_smooth_v_predictor_4x8_c, svt_aom_smooth_v_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_4x16, svt_aom_smooth_v_predictor_4x16_c, svt_aom_smooth_v_predictor_4x16_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_8x4, svt_aom_smooth_v_predictor_8x4_c, svt_aom_smooth_v_predictor_8x4_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_8x8, svt_aom_smooth_v_predictor_8x8_c, svt_aom_smooth_v_predictor_8x8_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_8x16, svt_aom_smooth_v_predictor_8x16_c, svt_aom_smooth_v_predictor_8x16_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_8x32, svt_aom_smooth_v_predictor_8x32_c, svt_aom_smooth_v_predictor_8x32_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_16x4, svt_aom_smooth_v_predictor_16x4_c, svt_aom_smooth_v_predictor_16x4_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_16x8, svt_aom_smooth_v_predictor_16x8_c, svt_aom_smooth_v_predictor_16x8_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_16x16, svt_aom_smooth_v_predictor_16x16_c, svt_aom_smooth_v_predictor_16x16_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_16x32, svt_aom_smooth_v_predictor_16x32_c, svt_aom_smooth_v_predictor_16x32_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_16x64, svt_aom_smooth_v_predictor_16x64_c, svt_aom_smooth_v_predictor_16x64_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_32x8, svt_aom_smooth_v_predictor_32x8_c, svt_aom_smooth_v_predictor_32x8_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_32x16, svt_aom_smooth_v_predictor_32x16_c, svt_aom_smooth_v_predictor_32x16_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_32x32, svt_aom_smooth_v_predictor_32x32_c, svt_aom_smooth_v_predictor_32x32_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_32x64, svt_aom_smooth_v_predictor_32x64_c, svt_aom_smooth_v_predictor_32x64_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_64x16, svt_aom_smooth_v_predictor_64x16_c, svt_aom_smooth_v_predictor_64x16_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_64x32, svt_aom_smooth_v_predictor_64x32_c, svt_aom_smooth_v_predictor_64x32_ssse3);
    SET_SSSE3(svt_aom_smooth_v_predictor_64x64, svt_aom_smooth_v_predictor_64x64_c, svt_aom_smooth_v_predictor_64x64_ssse3);

    SET_SSSE3(svt_aom_smooth_predictor_4x4, svt_aom_smooth_predictor_4x4_c, svt_aom_smooth_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_4x8, svt_aom_smooth_predictor_4x8_c, svt_aom_smooth_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_4x16, svt_aom_smooth_predictor_4x16_c, svt_aom_smooth_predictor_4x16_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_8x4, svt_aom_smooth_predictor_8x4_c, svt_aom_smooth_predictor_8x4_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_8x8, svt_aom_smooth_predictor_8x8_c, svt_aom_smooth_predictor_8x8_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_8x16, svt_aom_smooth_predictor_8x16_c, svt_aom_smooth_predictor_8x16_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_8x32, svt_aom_smooth_predictor_8x32_c, svt_aom_smooth_predictor_8x32_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_16x4, svt_aom_smooth_predictor_16x4_c, svt_aom_smooth_predictor_16x4_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_16x8, svt_aom_smooth_predictor_16x8_c, svt_aom_smooth_predictor_16x8_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_16x16, svt_aom_smooth_predictor_16x16_c, svt_aom_smooth_predictor_16x16_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_16x32, svt_aom_smooth_predictor_16x32_c, svt_aom_smooth_predictor_16x32_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_16x64, svt_aom_smooth_predictor_16x64_c, svt_aom_smooth_predictor_16x64_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_32x8, svt_aom_smooth_predictor_32x8_c, svt_aom_smooth_predictor_32x8_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_32x16, svt_aom_smooth_predictor_32x16_c, svt_aom_smooth_predictor_32x16_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_32x32, svt_aom_smooth_predictor_32x32_c, svt_aom_smooth_predictor_32x32_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_32x64, svt_aom_smooth_predictor_32x64_c, svt_aom_smooth_predictor_32x64_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_64x16, svt_aom_smooth_predictor_64x16_c, svt_aom_smooth_predictor_64x16_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_64x32, svt_aom_smooth_predictor_64x32_c, svt_aom_smooth_predictor_64x32_ssse3);
    SET_SSSE3(svt_aom_smooth_predictor_64x64, svt_aom_smooth_predictor_64x64_c, svt_aom_smooth_predictor_64x64_ssse3);

    SET_SSE2(svt_aom_v_predictor_4x4, svt_aom_v_predictor_4x4_c, svt_aom_v_predictor_4x4_sse2);
    SET_SSE2(svt_aom_v_predictor_4x8, svt_aom_v_predictor_4x8_c, svt_aom_v_predictor_4x8_sse2);
    SET_SSE2(svt_aom_v_predictor_4x16, svt_aom_v_predictor_4x16_c, svt_aom_v_predictor_4x16_sse2);
    SET_SSE2(svt_aom_v_predictor_8x4, svt_aom_v_predictor_8x4_c, svt_aom_v_predictor_8x4_sse2);
    SET_SSE2(svt_aom_v_predictor_8x8, svt_aom_v_predictor_8x8_c, svt_aom_v_predictor_8x8_sse2);
    SET_SSE2(svt_aom_v_predictor_8x16, svt_aom_v_predictor_8x16_c, svt_aom_v_predictor_8x16_sse2);
    SET_SSE2(svt_aom_v_predictor_8x32, svt_aom_v_predictor_8x32_c, svt_aom_v_predictor_8x32_sse2);
    SET_SSE2(svt_aom_v_predictor_16x4, svt_aom_v_predictor_16x4_c, svt_aom_v_predictor_16x4_sse2);
    SET_SSE2(svt_aom_v_predictor_16x8, svt_aom_v_predictor_16x8_c, svt_aom_v_predictor_16x8_sse2);
    SET_SSE2(svt_aom_v_predictor_16x16, svt_aom_v_predictor_16x16_c, svt_aom_v_predictor_16x16_sse2);
    SET_SSE2(svt_aom_v_predictor_16x32, svt_aom_v_predictor_16x32_c, svt_aom_v_predictor_16x32_sse2);
    SET_SSE2(svt_aom_v_predictor_16x64, svt_aom_v_predictor_16x64_c, svt_aom_v_predictor_16x64_sse2);
    SET_SSE2(svt_aom_v_predictor_32x8, svt_aom_v_predictor_32x8_c, svt_aom_v_predictor_32x8_sse2);
    SET_AVX2(svt_aom_v_predictor_32x16, svt_aom_v_predictor_32x16_c, svt_aom_v_predictor_32x16_avx2);
    SET_AVX2(svt_aom_v_predictor_32x32, svt_aom_v_predictor_32x32_c, svt_aom_v_predictor_32x32_avx2);
    SET_AVX2(svt_aom_v_predictor_32x64, svt_aom_v_predictor_32x64_c, svt_aom_v_predictor_32x64_avx2);
    SET_AVX2(svt_aom_v_predictor_64x16, svt_aom_v_predictor_64x16_c, svt_aom_v_predictor_64x16_avx2);
    SET_AVX2(svt_aom_v_predictor_64x32, svt_aom_v_predictor_64x32_c, svt_aom_v_predictor_64x32_avx2);
    SET_AVX2(svt_aom_v_predictor_64x64, svt_aom_v_predictor_64x64_c, svt_aom_v_predictor_64x64_avx2);

    SET_SSE2(svt_aom_h_predictor_4x4, svt_aom_h_predictor_4x4_c, svt_aom_h_predictor_4x4_sse2);
    SET_SSE2(svt_aom_h_predictor_4x8, svt_aom_h_predictor_4x8_c, svt_aom_h_predictor_4x8_sse2);
    SET_SSE2(svt_aom_h_predictor_4x16, svt_aom_h_predictor_4x16_c, svt_aom_h_predictor_4x16_sse2);
    SET_SSE2(svt_aom_h_predictor_8x4, svt_aom_h_predictor_8x4_c, svt_aom_h_predictor_8x4_sse2);
    SET_SSE2(svt_aom_h_predictor_8x8, svt_aom_h_predictor_8x8_c, svt_aom_h_predictor_8x8_sse2);
    SET_SSE2(svt_aom_h_predictor_8x16, svt_aom_h_predictor_8x16_c, svt_aom_h_predictor_8x16_sse2);
    SET_SSE2(svt_aom_h_predictor_8x32, svt_aom_h_predictor_8x32_c, svt_aom_h_predictor_8x32_sse2);
    SET_SSE2(svt_aom_h_predictor_16x4, svt_aom_h_predictor_16x4_c, svt_aom_h_predictor_16x4_sse2);
    SET_SSE2(svt_aom_h_predictor_16x8, svt_aom_h_predictor_16x8_c, svt_aom_h_predictor_16x8_sse2);
    SET_SSE2(svt_aom_h_predictor_16x16, svt_aom_h_predictor_16x16_c, svt_aom_h_predictor_16x16_sse2);
    SET_SSE2(svt_aom_h_predictor_16x32, svt_aom_h_predictor_16x32_c, svt_aom_h_predictor_16x32_sse2);
    SET_SSE2(svt_aom_h_predictor_16x64, svt_aom_h_predictor_16x64_c, svt_aom_h_predictor_16x64_sse2);
    SET_SSE2(svt_aom_h_predictor_32x8, svt_aom_h_predictor_32x8_c, svt_aom_h_predictor_32x8_sse2);
    SET_SSE2(svt_aom_h_predictor_32x16, svt_aom_h_predictor_32x16_c, svt_aom_h_predictor_32x16_sse2);
    SET_AVX2(svt_aom_h_predictor_32x32, svt_aom_h_predictor_32x32_c, svt_aom_h_predictor_32x32_avx2);
    SET_SSE2(svt_aom_h_predictor_32x64, svt_aom_h_predictor_32x64_c, svt_aom_h_predictor_32x64_sse2);
    SET_SSE2(svt_aom_h_predictor_64x16, svt_aom_h_predictor_64x16_c, svt_aom_h_predictor_64x16_sse2);
    SET_SSE2(svt_aom_h_predictor_64x32, svt_aom_h_predictor_64x32_c, svt_aom_h_predictor_64x32_sse2);
    SET_SSE2(svt_aom_h_predictor_64x64, svt_aom_h_predictor_64x64_c, svt_aom_h_predictor_64x64_sse2);

    SET_AVX2(svt_cdef_find_dir, svt_cdef_find_dir_c, svt_cdef_find_dir_avx2);
    SET_AVX2(svt_cdef_filter_block, svt_cdef_filter_block_c, svt_cdef_filter_block_avx2);
    /* No C version, use only internal in kerneal: svt_cdef_filter_block_avx2() */
#ifdef ARCH_X86_64
    if (flags & HAS_AVX2)    svt_cdef_filter_block_8x8_16 = svt_cdef_filter_block_8x8_16_avx2;
#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) svt_cdef_filter_block_8x8_16 = svt_cdef_filter_block_8x8_16_avx512;
#endif
#endif

    SET_AVX2(svt_copy_rect8_8bit_to_16bit, svt_copy_rect8_8bit_to_16bit_c, svt_copy_rect8_8bit_to_16bit_avx2);
    SET_AVX2(svt_av1_highbd_warp_affine, svt_av1_highbd_warp_affine_c, svt_av1_highbd_warp_affine_avx2);
    SET_AVX2(svt_av1_warp_affine, svt_av1_warp_affine_c, svt_av1_warp_affine_avx2);

    SET_SSE2(svt_aom_highbd_lpf_horizontal_4, svt_aom_highbd_lpf_horizontal_4_c, svt_aom_highbd_lpf_horizontal_4_sse2);
    SET_SSE2(svt_aom_highbd_lpf_horizontal_6, svt_aom_highbd_lpf_horizontal_6_c, svt_aom_highbd_lpf_horizontal_6_sse2);
    SET_SSE2(svt_aom_highbd_lpf_horizontal_8, svt_aom_highbd_lpf_horizontal_8_c, svt_aom_highbd_lpf_horizontal_8_sse2);
    SET_SSE2(svt_aom_highbd_lpf_horizontal_14, svt_aom_highbd_lpf_horizontal_14_c, svt_aom_highbd_lpf_horizontal_14_sse2);
    SET_SSE2(svt_aom_highbd_lpf_vertical_4, svt_aom_highbd_lpf_vertical_4_c, svt_aom_highbd_lpf_vertical_4_sse2);
    SET_SSE2(svt_aom_highbd_lpf_vertical_6, svt_aom_highbd_lpf_vertical_6_c, svt_aom_highbd_lpf_vertical_6_sse2);
    SET_SSE2(svt_aom_highbd_lpf_vertical_8, svt_aom_highbd_lpf_vertical_8_c, svt_aom_highbd_lpf_vertical_8_sse2);
    SET_SSE2(svt_aom_highbd_lpf_vertical_14, svt_aom_highbd_lpf_vertical_14_c, svt_aom_highbd_lpf_vertical_14_sse2);
    SET_SSE2(svt_aom_lpf_horizontal_4, svt_aom_lpf_horizontal_4_c, svt_aom_lpf_horizontal_4_sse2);
    SET_SSE2(svt_aom_lpf_horizontal_6, svt_aom_lpf_horizontal_6_c, svt_aom_lpf_horizontal_6_sse2);
    SET_SSE2(svt_aom_lpf_horizontal_8, svt_aom_lpf_horizontal_8_c, svt_aom_lpf_horizontal_8_sse2);
    SET_SSE2(svt_aom_lpf_horizontal_14, svt_aom_lpf_horizontal_14_c, svt_aom_lpf_horizontal_14_sse2);
    SET_SSE2(svt_aom_lpf_vertical_4, svt_aom_lpf_vertical_4_c, svt_aom_lpf_vertical_4_sse2);
    SET_SSE2(svt_aom_lpf_vertical_6, svt_aom_lpf_vertical_6_c, svt_aom_lpf_vertical_6_sse2);
    SET_SSE2(svt_aom_lpf_vertical_8, svt_aom_lpf_vertical_8_c, svt_aom_lpf_vertical_8_sse2);
    SET_SSE2(svt_aom_lpf_vertical_14, svt_aom_lpf_vertical_14_c, svt_aom_lpf_vertical_14_sse2);

    // svt_aom_highbd_v_predictor
    SET_SSE2(svt_aom_highbd_v_predictor_4x4, svt_aom_highbd_v_predictor_4x4_c, svt_aom_highbd_v_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_4x8, svt_aom_highbd_v_predictor_4x8_c, svt_aom_highbd_v_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_4x16, svt_aom_highbd_v_predictor_4x16_c, svt_aom_highbd_v_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_8x4, svt_aom_highbd_v_predictor_8x4_c, svt_aom_highbd_v_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_8x8, svt_aom_highbd_v_predictor_8x8_c, svt_aom_highbd_v_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_8x16, svt_aom_highbd_v_predictor_8x16_c, svt_aom_highbd_v_predictor_8x16_sse2);
    SET_SSE2(svt_aom_highbd_v_predictor_8x32, svt_aom_highbd_v_predictor_8x32_c, svt_aom_highbd_v_predictor_8x32_sse2);
    SET_AVX2(svt_aom_highbd_v_predictor_16x4, svt_aom_highbd_v_predictor_16x4_c, svt_aom_highbd_v_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_v_predictor_16x8, svt_aom_highbd_v_predictor_16x8_c, svt_aom_highbd_v_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_v_predictor_16x16, svt_aom_highbd_v_predictor_16x16_c, svt_aom_highbd_v_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_v_predictor_16x32, svt_aom_highbd_v_predictor_16x32_c, svt_aom_highbd_v_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_v_predictor_16x64, svt_aom_highbd_v_predictor_16x64_c, svt_aom_highbd_v_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_32x8, svt_aom_highbd_v_predictor_32x8_c, svt_aom_highbd_v_predictor_32x8_avx2, aom_highbd_v_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_32x16, svt_aom_highbd_v_predictor_32x16_c, svt_aom_highbd_v_predictor_32x16_avx2, aom_highbd_v_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_32x32, svt_aom_highbd_v_predictor_32x32_c, svt_aom_highbd_v_predictor_32x32_avx2, aom_highbd_v_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_32x64, svt_aom_highbd_v_predictor_32x64_c, svt_aom_highbd_v_predictor_32x64_avx2, aom_highbd_v_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_64x16, svt_aom_highbd_v_predictor_64x16_c, svt_aom_highbd_v_predictor_64x16_avx2, aom_highbd_v_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_64x32, svt_aom_highbd_v_predictor_64x32_c, svt_aom_highbd_v_predictor_64x32_avx2, aom_highbd_v_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_v_predictor_64x64, svt_aom_highbd_v_predictor_64x64_c, svt_aom_highbd_v_predictor_64x64_avx2, aom_highbd_v_predictor_64x64_avx512);

    //aom_highbd_smooth_predictor
    SET_ONLY_C(svt_aom_highbd_smooth_predictor_2x2, svt_aom_highbd_smooth_predictor_2x2_c);
    SET_SSSE3(svt_aom_highbd_smooth_predictor_4x4, svt_aom_highbd_smooth_predictor_4x4_c, svt_aom_highbd_smooth_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_predictor_4x8, svt_aom_highbd_smooth_predictor_4x8_c, svt_aom_highbd_smooth_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_predictor_4x16, svt_aom_highbd_smooth_predictor_4x16_c, svt_aom_highbd_smooth_predictor_4x16_ssse3);
    SET_AVX2(svt_aom_highbd_smooth_predictor_8x4, svt_aom_highbd_smooth_predictor_8x4_c, svt_aom_highbd_smooth_predictor_8x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_8x8, svt_aom_highbd_smooth_predictor_8x8_c, svt_aom_highbd_smooth_predictor_8x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_8x16, svt_aom_highbd_smooth_predictor_8x16_c, svt_aom_highbd_smooth_predictor_8x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_8x32, svt_aom_highbd_smooth_predictor_8x32_c, svt_aom_highbd_smooth_predictor_8x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_16x4, svt_aom_highbd_smooth_predictor_16x4_c, svt_aom_highbd_smooth_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_16x8, svt_aom_highbd_smooth_predictor_16x8_c, svt_aom_highbd_smooth_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_16x16, svt_aom_highbd_smooth_predictor_16x16_c, svt_aom_highbd_smooth_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_16x32, svt_aom_highbd_smooth_predictor_16x32_c, svt_aom_highbd_smooth_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_predictor_16x64, svt_aom_highbd_smooth_predictor_16x64_c, svt_aom_highbd_smooth_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_32x8, svt_aom_highbd_smooth_predictor_32x8_c, svt_aom_highbd_smooth_predictor_32x8_avx2, aom_highbd_smooth_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_32x16, svt_aom_highbd_smooth_predictor_32x16_c, svt_aom_highbd_smooth_predictor_32x16_avx2, aom_highbd_smooth_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_32x32, svt_aom_highbd_smooth_predictor_32x32_c, svt_aom_highbd_smooth_predictor_32x32_avx2, aom_highbd_smooth_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_32x64, svt_aom_highbd_smooth_predictor_32x64_c, svt_aom_highbd_smooth_predictor_32x64_avx2, aom_highbd_smooth_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_64x16, svt_aom_highbd_smooth_predictor_64x16_c, svt_aom_highbd_smooth_predictor_64x16_avx2, aom_highbd_smooth_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_64x32, svt_aom_highbd_smooth_predictor_64x32_c, svt_aom_highbd_smooth_predictor_64x32_avx2, aom_highbd_smooth_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_predictor_64x64, svt_aom_highbd_smooth_predictor_64x64_c, svt_aom_highbd_smooth_predictor_64x64_avx2, aom_highbd_smooth_predictor_64x64_avx512);

    //aom_highbd_smooth_h_predictor
    SET_SSSE3(svt_aom_highbd_smooth_h_predictor_4x4, svt_aom_highbd_smooth_h_predictor_4x4_c, svt_aom_highbd_smooth_h_predictor_4x4_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_h_predictor_4x8, svt_aom_highbd_smooth_h_predictor_4x8_c, svt_aom_highbd_smooth_h_predictor_4x8_ssse3);
    SET_SSSE3(svt_aom_highbd_smooth_h_predictor_4x16, svt_aom_highbd_smooth_h_predictor_4x16_c, svt_aom_highbd_smooth_h_predictor_4x16_ssse3);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_8x4, svt_aom_highbd_smooth_h_predictor_8x4_c, svt_aom_highbd_smooth_h_predictor_8x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_8x8, svt_aom_highbd_smooth_h_predictor_8x8_c, svt_aom_highbd_smooth_h_predictor_8x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_8x16, svt_aom_highbd_smooth_h_predictor_8x16_c, svt_aom_highbd_smooth_h_predictor_8x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_8x32, svt_aom_highbd_smooth_h_predictor_8x32_c, svt_aom_highbd_smooth_h_predictor_8x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_16x4, svt_aom_highbd_smooth_h_predictor_16x4_c, svt_aom_highbd_smooth_h_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_16x8, svt_aom_highbd_smooth_h_predictor_16x8_c, svt_aom_highbd_smooth_h_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_16x16, svt_aom_highbd_smooth_h_predictor_16x16_c, svt_aom_highbd_smooth_h_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_16x32, svt_aom_highbd_smooth_h_predictor_16x32_c, svt_aom_highbd_smooth_h_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_smooth_h_predictor_16x64, svt_aom_highbd_smooth_h_predictor_16x64_c, svt_aom_highbd_smooth_h_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_32x8, svt_aom_highbd_smooth_h_predictor_32x8_c, svt_aom_highbd_smooth_h_predictor_32x8_avx2, aom_highbd_smooth_h_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_32x16, svt_aom_highbd_smooth_h_predictor_32x16_c, svt_aom_highbd_smooth_h_predictor_32x16_avx2, aom_highbd_smooth_h_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_32x32, svt_aom_highbd_smooth_h_predictor_32x32_c, svt_aom_highbd_smooth_h_predictor_32x32_avx2, aom_highbd_smooth_h_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_32x64, svt_aom_highbd_smooth_h_predictor_32x64_c, svt_aom_highbd_smooth_h_predictor_32x64_avx2, aom_highbd_smooth_h_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_64x16, svt_aom_highbd_smooth_h_predictor_64x16_c, svt_aom_highbd_smooth_h_predictor_64x16_avx2, aom_highbd_smooth_h_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_64x32, svt_aom_highbd_smooth_h_predictor_64x32_c, svt_aom_highbd_smooth_h_predictor_64x32_avx2, aom_highbd_smooth_h_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_smooth_h_predictor_64x64, svt_aom_highbd_smooth_h_predictor_64x64_c, svt_aom_highbd_smooth_h_predictor_64x64_avx2, aom_highbd_smooth_h_predictor_64x64_avx512);

    //aom_highbd_dc_128_predictor
    SET_SSE2(svt_aom_highbd_dc_128_predictor_4x4, svt_aom_highbd_dc_128_predictor_4x4_c, svt_aom_highbd_dc_128_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_4x8, svt_aom_highbd_dc_128_predictor_4x8_c, svt_aom_highbd_dc_128_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_4x16, svt_aom_highbd_dc_128_predictor_4x16_c, svt_aom_highbd_dc_128_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_8x4, svt_aom_highbd_dc_128_predictor_8x4_c, svt_aom_highbd_dc_128_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_8x8, svt_aom_highbd_dc_128_predictor_8x8_c, svt_aom_highbd_dc_128_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_8x16, svt_aom_highbd_dc_128_predictor_8x16_c, svt_aom_highbd_dc_128_predictor_8x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_128_predictor_8x32, svt_aom_highbd_dc_128_predictor_8x32_c, svt_aom_highbd_dc_128_predictor_8x32_sse2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_16x4, svt_aom_highbd_dc_128_predictor_16x4_c, svt_aom_highbd_dc_128_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_16x8, svt_aom_highbd_dc_128_predictor_16x8_c, svt_aom_highbd_dc_128_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_16x16, svt_aom_highbd_dc_128_predictor_16x16_c, svt_aom_highbd_dc_128_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_16x32, svt_aom_highbd_dc_128_predictor_16x32_c, svt_aom_highbd_dc_128_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_16x64, svt_aom_highbd_dc_128_predictor_16x64_c, svt_aom_highbd_dc_128_predictor_16x64_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_32x8, svt_aom_highbd_dc_128_predictor_32x8_c, svt_aom_highbd_dc_128_predictor_32x8_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_32x16, svt_aom_highbd_dc_128_predictor_32x16_c, svt_aom_highbd_dc_128_predictor_32x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_32x32, svt_aom_highbd_dc_128_predictor_32x32_c, svt_aom_highbd_dc_128_predictor_32x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_32x64, svt_aom_highbd_dc_128_predictor_32x64_c, svt_aom_highbd_dc_128_predictor_32x64_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_64x16, svt_aom_highbd_dc_128_predictor_64x16_c, svt_aom_highbd_dc_128_predictor_64x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_64x32, svt_aom_highbd_dc_128_predictor_64x32_c, svt_aom_highbd_dc_128_predictor_64x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_128_predictor_64x64, svt_aom_highbd_dc_128_predictor_64x64_c, svt_aom_highbd_dc_128_predictor_64x64_avx2);

    //aom_highbd_dc_left_predictor
    SET_ONLY_C(svt_aom_highbd_dc_left_predictor_2x2, svt_aom_highbd_dc_left_predictor_2x2_c);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_4x4, svt_aom_highbd_dc_left_predictor_4x4_c, svt_aom_highbd_dc_left_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_4x8, svt_aom_highbd_dc_left_predictor_4x8_c, svt_aom_highbd_dc_left_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_4x16, svt_aom_highbd_dc_left_predictor_4x16_c, svt_aom_highbd_dc_left_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_8x4, svt_aom_highbd_dc_left_predictor_8x4_c, svt_aom_highbd_dc_left_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_8x8, svt_aom_highbd_dc_left_predictor_8x8_c, svt_aom_highbd_dc_left_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_8x16, svt_aom_highbd_dc_left_predictor_8x16_c, svt_aom_highbd_dc_left_predictor_8x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_left_predictor_8x32, svt_aom_highbd_dc_left_predictor_8x32_c, svt_aom_highbd_dc_left_predictor_8x32_sse2);
    SET_AVX2(svt_aom_highbd_dc_left_predictor_16x4, svt_aom_highbd_dc_left_predictor_16x4_c, svt_aom_highbd_dc_left_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_dc_left_predictor_16x8, svt_aom_highbd_dc_left_predictor_16x8_c, svt_aom_highbd_dc_left_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_dc_left_predictor_16x16, svt_aom_highbd_dc_left_predictor_16x16_c, svt_aom_highbd_dc_left_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_left_predictor_16x32, svt_aom_highbd_dc_left_predictor_16x32_c, svt_aom_highbd_dc_left_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_left_predictor_16x64, svt_aom_highbd_dc_left_predictor_16x64_c, svt_aom_highbd_dc_left_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_32x8, svt_aom_highbd_dc_left_predictor_32x8_c, svt_aom_highbd_dc_left_predictor_32x8_avx2, aom_highbd_dc_left_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_32x16, svt_aom_highbd_dc_left_predictor_32x16_c, svt_aom_highbd_dc_left_predictor_32x16_avx2, aom_highbd_dc_left_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_32x32, svt_aom_highbd_dc_left_predictor_32x32_c, svt_aom_highbd_dc_left_predictor_32x32_avx2, aom_highbd_dc_left_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_32x64, svt_aom_highbd_dc_left_predictor_32x64_c, svt_aom_highbd_dc_left_predictor_32x64_avx2, aom_highbd_dc_left_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_64x16, svt_aom_highbd_dc_left_predictor_64x16_c, svt_aom_highbd_dc_left_predictor_64x16_avx2, aom_highbd_dc_left_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_64x32, svt_aom_highbd_dc_left_predictor_64x32_c, svt_aom_highbd_dc_left_predictor_64x32_avx2, aom_highbd_dc_left_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_left_predictor_64x64, svt_aom_highbd_dc_left_predictor_64x64_c, svt_aom_highbd_dc_left_predictor_64x64_avx2, aom_highbd_dc_left_predictor_64x64_avx512);

    SET_ONLY_C(svt_aom_highbd_dc_predictor_2x2, svt_aom_highbd_dc_predictor_2x2_c);
    SET_SSE2(svt_aom_highbd_dc_predictor_4x4, svt_aom_highbd_dc_predictor_4x4_c, svt_aom_highbd_dc_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_4x8, svt_aom_highbd_dc_predictor_4x8_c, svt_aom_highbd_dc_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_4x16, svt_aom_highbd_dc_predictor_4x16_c, svt_aom_highbd_dc_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_8x4, svt_aom_highbd_dc_predictor_8x4_c, svt_aom_highbd_dc_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_8x8, svt_aom_highbd_dc_predictor_8x8_c, svt_aom_highbd_dc_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_8x16, svt_aom_highbd_dc_predictor_8x16_c, svt_aom_highbd_dc_predictor_8x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_predictor_8x32, svt_aom_highbd_dc_predictor_8x32_c, svt_aom_highbd_dc_predictor_8x32_sse2);
    SET_AVX2(svt_aom_highbd_dc_predictor_16x4, svt_aom_highbd_dc_predictor_16x4_c, svt_aom_highbd_dc_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_dc_predictor_16x8, svt_aom_highbd_dc_predictor_16x8_c, svt_aom_highbd_dc_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_dc_predictor_16x16, svt_aom_highbd_dc_predictor_16x16_c, svt_aom_highbd_dc_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_predictor_16x32, svt_aom_highbd_dc_predictor_16x32_c, svt_aom_highbd_dc_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_predictor_16x64, svt_aom_highbd_dc_predictor_16x64_c, svt_aom_highbd_dc_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_32x8, svt_aom_highbd_dc_predictor_32x8_c, svt_aom_highbd_dc_predictor_32x8_avx2, aom_highbd_dc_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_32x16, svt_aom_highbd_dc_predictor_32x16_c, svt_aom_highbd_dc_predictor_32x16_avx2, aom_highbd_dc_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_32x32, svt_aom_highbd_dc_predictor_32x32_c, svt_aom_highbd_dc_predictor_32x32_avx2, aom_highbd_dc_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_32x64, svt_aom_highbd_dc_predictor_32x64_c, svt_aom_highbd_dc_predictor_32x64_avx2, aom_highbd_dc_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_64x16, svt_aom_highbd_dc_predictor_64x16_c, svt_aom_highbd_dc_predictor_64x16_avx2, aom_highbd_dc_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_64x32, svt_aom_highbd_dc_predictor_64x32_c, svt_aom_highbd_dc_predictor_64x32_avx2, aom_highbd_dc_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_predictor_64x64, svt_aom_highbd_dc_predictor_64x64_c, svt_aom_highbd_dc_predictor_64x64_avx2, aom_highbd_dc_predictor_64x64_avx512);

    //aom_highbd_dc_top_predictor
    SET_SSE2(svt_aom_highbd_dc_top_predictor_4x4, svt_aom_highbd_dc_top_predictor_4x4_c, svt_aom_highbd_dc_top_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_top_predictor_4x8, svt_aom_highbd_dc_top_predictor_4x8_c, svt_aom_highbd_dc_top_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_top_predictor_4x16, svt_aom_highbd_dc_top_predictor_4x16_c, svt_aom_highbd_dc_top_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_dc_top_predictor_8x4, svt_aom_highbd_dc_top_predictor_8x4_c, svt_aom_highbd_dc_top_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_dc_top_predictor_8x8, svt_aom_highbd_dc_top_predictor_8x8_c, svt_aom_highbd_dc_top_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_dc_top_predictor_8x16, svt_aom_highbd_dc_top_predictor_8x16_c, svt_aom_highbd_dc_top_predictor_8x16_sse2);
    SET_ONLY_C(svt_aom_highbd_dc_top_predictor_8x32, svt_aom_highbd_dc_top_predictor_8x32_c);
    SET_AVX2(svt_aom_highbd_dc_top_predictor_16x4, svt_aom_highbd_dc_top_predictor_16x4_c, svt_aom_highbd_dc_top_predictor_16x4_avx2);
    SET_AVX2(svt_aom_highbd_dc_top_predictor_16x8, svt_aom_highbd_dc_top_predictor_16x8_c, svt_aom_highbd_dc_top_predictor_16x8_avx2);
    SET_AVX2(svt_aom_highbd_dc_top_predictor_16x16, svt_aom_highbd_dc_top_predictor_16x16_c, svt_aom_highbd_dc_top_predictor_16x16_avx2);
    SET_AVX2(svt_aom_highbd_dc_top_predictor_16x32, svt_aom_highbd_dc_top_predictor_16x32_c, svt_aom_highbd_dc_top_predictor_16x32_avx2);
    SET_AVX2(svt_aom_highbd_dc_top_predictor_16x64, svt_aom_highbd_dc_top_predictor_16x64_c, svt_aom_highbd_dc_top_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_32x8, svt_aom_highbd_dc_top_predictor_32x8_c, svt_aom_highbd_dc_top_predictor_32x8_avx2, aom_highbd_dc_top_predictor_32x8_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_32x16, svt_aom_highbd_dc_top_predictor_32x16_c, svt_aom_highbd_dc_top_predictor_32x16_avx2, aom_highbd_dc_top_predictor_32x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_32x32, svt_aom_highbd_dc_top_predictor_32x32_c, svt_aom_highbd_dc_top_predictor_32x32_avx2, aom_highbd_dc_top_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_32x64, svt_aom_highbd_dc_top_predictor_32x64_c, svt_aom_highbd_dc_top_predictor_32x64_avx2, aom_highbd_dc_top_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_64x16, svt_aom_highbd_dc_top_predictor_64x16_c, svt_aom_highbd_dc_top_predictor_64x16_avx2, aom_highbd_dc_top_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_64x32, svt_aom_highbd_dc_top_predictor_64x32_c, svt_aom_highbd_dc_top_predictor_64x32_avx2, aom_highbd_dc_top_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_dc_top_predictor_64x64, svt_aom_highbd_dc_top_predictor_64x64_c, svt_aom_highbd_dc_top_predictor_64x64_avx2, aom_highbd_dc_top_predictor_64x64_avx512);

    // svt_aom_highbd_h_predictor
    SET_SSE2(svt_aom_highbd_h_predictor_4x4, svt_aom_highbd_h_predictor_4x4_c, svt_aom_highbd_h_predictor_4x4_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_4x8, svt_aom_highbd_h_predictor_4x8_c, svt_aom_highbd_h_predictor_4x8_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_4x16, svt_aom_highbd_h_predictor_4x16_c, svt_aom_highbd_h_predictor_4x16_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_8x4, svt_aom_highbd_h_predictor_8x4_c, svt_aom_highbd_h_predictor_8x4_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_8x8, svt_aom_highbd_h_predictor_8x8_c, svt_aom_highbd_h_predictor_8x8_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_8x16, svt_aom_highbd_h_predictor_8x16_c, svt_aom_highbd_h_predictor_8x16_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_8x32, svt_aom_highbd_h_predictor_8x32_c, svt_aom_highbd_h_predictor_8x32_sse2);
    SET_AVX2(svt_aom_highbd_h_predictor_16x4, svt_aom_highbd_h_predictor_16x4_c, svt_aom_highbd_h_predictor_16x4_avx2);
    SET_SSE2(svt_aom_highbd_h_predictor_16x8, svt_aom_highbd_h_predictor_16x8_c, svt_aom_highbd_h_predictor_16x8_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_16x16, svt_aom_highbd_h_predictor_16x16_c, svt_aom_highbd_h_predictor_16x16_sse2);
    SET_SSE2(svt_aom_highbd_h_predictor_16x32, svt_aom_highbd_h_predictor_16x32_c, svt_aom_highbd_h_predictor_16x32_sse2);
    SET_AVX2(svt_aom_highbd_h_predictor_16x64, svt_aom_highbd_h_predictor_16x64_c, svt_aom_highbd_h_predictor_16x64_avx2);
    SET_AVX2_AVX512(svt_aom_highbd_h_predictor_32x8, svt_aom_highbd_h_predictor_32x8_c, svt_aom_highbd_h_predictor_32x8_avx2, aom_highbd_h_predictor_32x8_avx512);
    SET_SSE2_AVX512(svt_aom_highbd_h_predictor_32x16, svt_aom_highbd_h_predictor_32x16_c, svt_aom_highbd_h_predictor_32x16_sse2, aom_highbd_h_predictor_32x16_avx512);
    SET_SSE2_AVX512(svt_aom_highbd_h_predictor_32x32, svt_aom_highbd_h_predictor_32x32_c, svt_aom_highbd_h_predictor_32x32_sse2, aom_highbd_h_predictor_32x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_h_predictor_32x64, svt_aom_highbd_h_predictor_32x64_c, svt_aom_highbd_h_predictor_32x64_avx2, aom_highbd_h_predictor_32x64_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_h_predictor_64x16, svt_aom_highbd_h_predictor_64x16_c, svt_aom_highbd_h_predictor_64x16_avx2, aom_highbd_h_predictor_64x16_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_h_predictor_64x32, svt_aom_highbd_h_predictor_64x32_c, svt_aom_highbd_h_predictor_64x32_avx2, aom_highbd_h_predictor_64x32_avx512);
    SET_AVX2_AVX512(svt_aom_highbd_h_predictor_64x64, svt_aom_highbd_h_predictor_64x64_c, svt_aom_highbd_h_predictor_64x64_avx2, aom_highbd_h_predictor_64x64_avx512);
    SET_SSE2(svt_log2f, log2f_32, Log2f_ASM);
    SET_SSE2(svt_memcpy, svt_memcpy_c, svt_memcpy_intrin_sse);

}
