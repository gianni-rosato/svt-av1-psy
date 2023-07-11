#
# Copyright(c) 2019 Intel Corporation
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
#

# Encoder/ASM_AVX2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(../../../API
        ../../Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/)
link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/)

check_both_flags_add(-mavx2)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
        # Intel Windows (*Note - The Warning level /W0 should be made to /W4 at some point)
        check_both_flags_add(/Qdiag-disable:10010,10148,10157 /W0)
    else()
        check_both_flags_add(-static-intel -w)
    endif()
endif()

if(MSVC)
    check_both_flags_add(/arch:AVX2)
endif()

set(all_files
    EbCdef_AVX2.c
    EbCombinedAveragingSAD_Intrinsic_AVX2.c
    EbComputeSAD_AVX2.h
    EbComputeSAD_Intrinsic_AVX2.c
    EbNoiseModel_avx2.c
    EbRestorationPick_AVX2.c
    EbTemporalFiltering_AVX2.c
    EbTransforms_Intrinsic_AVX2.c
    av1_quantize_avx2.c
    corner_match_avx2.c
    dwt_avx2.c
    enc_warp_plane_avx2.c
    encodetxb_avx2.c
    fft_avx2.c
    highbd_fwd_txfm_avx2.c
    highbd_quantize_intrin_avx2.c
    highbd_variance_avx2.c
    obmc_sad_avx2.c
    obmc_variance_avx2.c
    palette_avx2.c
    pickrst_avx2.c
    pickrst_avx2.h
    resize_avx2.c
    sse_avx2.c
    variance_avx2.c
    variance_impl_avx2.c
    ssim_avx2.c
    )

add_library(ENCODER_ASM_AVX2 OBJECT ${all_files})
