/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#define RTCD_C
#include "common_dsp_rtcd.h"
#include "EbPictureOperators.h"
#include "EbPackUnPack_C.h"
#include "EbPackUnPack_SSE2.h"
#include "EbPackUnPack_AVX2.h"
#include "EbMcp_SSE2.h"
#include "EbAvcStyleMcp_SSSE3.h"
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

/**************************************
 * Instruction Set Support
 **************************************/

// Helper Functions
static INLINE void run_cpuid(uint32_t eax, uint32_t ecx, int32_t* abcd) {
#ifdef _WIN32
    __cpuidex(abcd, eax, ecx);
#else
    uint32_t ebx = 0, edx = 0;
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
    int32_t abcd[4];
    int32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    int32_t avx2_bmi12_mask        = (1 << 5) | (1 << 3) | (1 << 8);

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

    int avx512_ebx_mask = (1 << 16) // AVX-512F
                          | (1 << 17) // AVX-512DQ
                          | (1 << 28) // AVX-512CD
                          | (1 << 30) // AVX-512BW
                          | (1 << 31); // AVX-512VL

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

#ifndef NON_AVX512_SUPPORT
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    do {                                                                                      \
        ptr = c;                                                                              \
        if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx;              \
        if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse;              \
        if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2;           \
        if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3;           \
        if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3;        \
        if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1;     \
        if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2;     \
        if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx;              \
        if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2;           \
        if (((uintptr_t)NULL != (uintptr_t)avx512) && (flags & HAS_AVX512F)) ptr = avx512;    \
    } while (0)
#else
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    do {                                                                                      \
        ptr = c;                                                                              \
        if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx;              \
        if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse;              \
        if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2;           \
        if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3;           \
        if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3;        \
        if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1;     \
        if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2;     \
        if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx;              \
        if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2;           \
    } while (0)
#endif

#define SET_SSE2(ptr, c, sse2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, 0, 0)
#define SET_SSE2_AVX2(ptr, c, sse2, avx2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, avx2, 0)
#define SET_SSSE3(ptr, c, ssse3) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, ssse3, 0, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41_AVX2(ptr, c, sse4_1, avx2) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, 0)
#define SET_SSE41_AVX2_AVX512(ptr, c, sse4_1, avx2, avx512) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, avx512)
#define SET_AVX2(ptr, c, avx2) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, 0)
#define SET_AVX2_AVX512(ptr, c, avx2, avx512) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, avx512)


void setup_common_rtcd_internal(CPU_FLAGS flags) {
    /** Should be done during library initialization,
        but for safe limiting cpu flags again. */
    flags &= get_cpu_flags_to_use();

    //to use C: flags=0

    aom_blend_a64_mask = aom_blend_a64_mask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_mask = aom_blend_a64_mask_sse4_1;
    if (flags & HAS_AVX2) aom_blend_a64_mask = aom_blend_a64_mask_avx2;
    aom_blend_a64_hmask = aom_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_hmask = aom_blend_a64_hmask_sse4_1;
    aom_blend_a64_vmask = aom_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_vmask = aom_blend_a64_vmask_sse4_1;

    aom_highbd_blend_a64_mask = aom_highbd_blend_a64_mask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_mask = aom_highbd_blend_a64_mask_sse4_1;
    aom_highbd_blend_a64_hmask = aom_highbd_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_hmask = aom_highbd_blend_a64_hmask_sse4_1;
    aom_highbd_blend_a64_vmask = aom_highbd_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_vmask = aom_highbd_blend_a64_vmask_sse4_1;

    eb_aom_highbd_blend_a64_vmask = eb_aom_highbd_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) eb_aom_highbd_blend_a64_vmask = eb_aom_highbd_blend_a64_vmask_sse4_1;
    eb_aom_highbd_blend_a64_hmask = eb_aom_highbd_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) eb_aom_highbd_blend_a64_hmask = eb_aom_highbd_blend_a64_hmask_sse4_1;

    eb_cfl_predict_lbd = eb_cfl_predict_lbd_c;
    if (flags & HAS_AVX2) eb_cfl_predict_lbd = eb_cfl_predict_lbd_avx2;
    eb_cfl_predict_hbd = eb_cfl_predict_hbd_c;
    if (flags & HAS_AVX2) eb_cfl_predict_hbd = eb_cfl_predict_hbd_avx2;

    eb_av1_filter_intra_predictor = eb_av1_filter_intra_predictor_c;
    if (flags & HAS_SSE4_1) eb_av1_filter_intra_predictor = eb_av1_filter_intra_predictor_sse4_1;

    eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_c;
    if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_sse4_1;

    eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_high_c_old;
    if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_sse4_1;

    eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_c;
    if (flags & HAS_SSE4_1) eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_sse4_1;

    av1_build_compound_diffwtd_mask_d16 = av1_build_compound_diffwtd_mask_d16_c;
    if (flags & HAS_AVX2)
        av1_build_compound_diffwtd_mask_d16 = av1_build_compound_diffwtd_mask_d16_avx2;


    eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_c;
    eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_c;
    eb_av1_inv_txfm2d_add_4x4   = eb_av1_inv_txfm2d_add_4x4_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_4x4 = eb_av1_inv_txfm2d_add_4x4_avx2;
    eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_c;
    eb_av1_inv_txfm2d_add_8x8   = eb_av1_inv_txfm2d_add_8x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x8 = eb_av1_inv_txfm2d_add_8x8_avx2;

    eb_av1_inv_txfm2d_add_8x16 = eb_av1_inv_txfm2d_add_8x16_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x16 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_16x8 = eb_av1_inv_txfm2d_add_16x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x8 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_16x32 = eb_av1_inv_txfm2d_add_16x32_c;
    eb_av1_inv_txfm2d_add_32x16 = eb_av1_inv_txfm2d_add_32x16_c;
    eb_av1_inv_txfm2d_add_32x8  = eb_av1_inv_txfm2d_add_32x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x8 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_8x32 = eb_av1_inv_txfm2d_add_8x32_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x32 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_32x64 = eb_av1_inv_txfm2d_add_32x64_c;
    eb_av1_inv_txfm2d_add_64x32 = eb_av1_inv_txfm2d_add_64x32_c;
    eb_av1_inv_txfm2d_add_16x64 = eb_av1_inv_txfm2d_add_16x64_c;
    eb_av1_inv_txfm2d_add_64x16 = eb_av1_inv_txfm2d_add_64x16_c;
    eb_av1_inv_txfm2d_add_4x8   = eb_av1_inv_txfm2d_add_4x8_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x8 = eb_av1_inv_txfm2d_add_4x8_sse4_1;
    eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_sse4_1;
    eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_sse4_1;
    eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_sse4_1;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_avx2;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_sse4_1;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x64 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x16 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x64 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x32 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x32 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x16 = eb_av1_highbd_inv_txfm_add_avx2;


#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_av1_inv_txfm2d_add_16x16   = eb_av1_inv_txfm2d_add_16x16_avx512;
        eb_av1_inv_txfm2d_add_32x32   = eb_av1_inv_txfm2d_add_32x32_avx512;
        eb_av1_inv_txfm2d_add_64x64   = eb_av1_inv_txfm2d_add_64x64_avx512;
        eb_av1_inv_txfm2d_add_16x64   = eb_av1_inv_txfm2d_add_16x64_avx512;
        eb_av1_inv_txfm2d_add_64x16   = eb_av1_inv_txfm2d_add_64x16_avx512;
        eb_av1_inv_txfm2d_add_32x64   = eb_av1_inv_txfm2d_add_32x64_avx512;
        eb_av1_inv_txfm2d_add_64x32   = eb_av1_inv_txfm2d_add_64x32_avx512;
        eb_av1_inv_txfm2d_add_16x32   = eb_av1_inv_txfm2d_add_16x32_avx512;
        eb_av1_inv_txfm2d_add_32x16   = eb_av1_inv_txfm2d_add_32x16_avx512;
    }
#endif

    eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_c;
    if (flags & HAS_SSSE3) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_ssse3;
    if (flags & HAS_AVX2) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_avx2;

    SET_AVX2(compressed_packmsb, compressed_packmsb_c, compressed_packmsb_avx2_intrin);
    SET_AVX2(c_pack, c_pack_c, c_pack_avx2_intrin);
    SET_SSE2_AVX2(unpack_avg, unpack_avg_c, unpack_avg_sse2_intrin, unpack_avg_avx2_intrin);
    SET_AVX2(unpack_avg_safe_sub, unpack_avg_safe_sub_c, unpack_avg_safe_sub_avx2_intrin);
    SET_AVX2(un_pack8_bit_data, un_pack8_bit_data_c, eb_enc_un_pack8_bit_data_avx2_intrin);

    SET_SSE2_AVX2(pack2d_16_bit_src_mul4,
                  eb_enc_msb_pack2_d,
                  eb_enc_msb_pack2d_sse2_intrin,
                  eb_enc_msb_pack2d_avx2_intrin_al);
    SET_SSE2(un_pack2d_16_bit_src_mul4, eb_enc_msb_un_pack2_d, eb_enc_msb_un_pack2d_sse2_intrin);

    SET_AVX2(full_distortion_kernel_cbf_zero32_bits,
             full_distortion_kernel_cbf_zero32_bits_c,
             full_distortion_kernel_cbf_zero32_bits_avx2);
    SET_AVX2(full_distortion_kernel32_bits,
             full_distortion_kernel32_bits_c,
             full_distortion_kernel32_bits_avx2);

    SET_AVX2_AVX512(spatial_full_distortion_kernel,
                    spatial_full_distortion_kernel_c,
                    spatial_full_distortion_kernel_avx2,
                    spatial_full_distortion_kernel_avx512);
    SET_AVX2(full_distortion_kernel16_bits,
             full_distortion_kernel16_bits_c,
             full_distortion_kernel16_bits_avx2);
    SET_AVX2_AVX512(residual_kernel8bit,
                    residual_kernel8bit_c,
                    residual_kernel8bit_avx2,
                    residual_kernel8bit_avx512);

    SET_SSE2(residual_kernel16bit, residual_kernel16bit_c, residual_kernel16bit_sse2_intrin);

    SET_SSE2(picture_average_kernel, picture_average_kernel_c, picture_average_kernel_sse2_intrin);
    SET_SSE2(picture_average_kernel1_line,
             picture_average_kernel1_line_c,
             picture_average_kernel1_line_sse2_intrin);

    SET_SSSE3(avc_style_luma_interpolation_filter,
              avc_style_luma_interpolation_filter_helper_c,
              avc_style_luma_interpolation_filter_helper_ssse3);
}
