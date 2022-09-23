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

# Common/ASM_AVX2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/)
link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)

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
    EbBlend_a64_mask_avx2.c
    EbHighbdIntraPrediction_AVX2.c
    EbInterPredictionCom_avx2.c
    EbIntraPrediction_AVX2.h
    EbIntraPrediction_Intrinsic_AVX2.c
    EbMemory_AVX2.h
    EbPackUnPack_Intrinsic_AVX2.c
    EbPictureOperators_Inline_AVX2.h
    EbPictureOperators_Intrinsic_AVX2.c
    aom_subpixel_8t_intrin_avx2.c
    av1_inv_txfm_avx2.c
    av1_inv_txfm_avx2.h
    cdef_block_avx2.c
    cfl_avx2.c
    convolve_2d_avx2.c
    convolve_avx2.c
    convolve_avx2.h
    highbd_convolve_2d_avx2.c
    highbd_convolve_avx2.c
    highbd_inv_txfm_avx2.c
    highbd_jnt_convolve_avx2.c
    highbd_warp_affine_avx2.c
    jnt_convolve_2d_avx2.c
    jnt_convolve_avx2.c
    selfguided_avx2.c
    synonyms_avx2.h
    transpose_avx2.h
    txfm_common_avx2.h
    warp_plane_avx2.c
    wiener_convolve_avx2.c
    wiener_convolve_avx2.h
    itx_lbd.h
    itx_hbd.h
    mc.h
    )

set(asm_files
    itx_avx2.asm
    itx16_avx2.asm
    mc_avx2.asm
    mc16_avx2.asm)

add_library(COMMON_ASM_AVX2 OBJECT ${all_files})

asm_compile_to_target(COMMON_ASM_AVX2 ${asm_files})
