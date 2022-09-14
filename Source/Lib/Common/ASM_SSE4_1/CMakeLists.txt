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

# ASM_SSE4.1 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/)

check_both_flags_add(-msse4.1)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel" AND NOT WIN32)
    check_both_flags_add(-static-intel -w)
endif()

set(all_files
    EbBlend_a64_mask_sse4.c
    EbBlend_sse4.h
    EbIntraPrediction16bit_Intrinsic_SSE4_1.c
    EbMemory_SSE4_1.h
    EbPictureOperators_Intrinsic_SSE4_1.c
    av1_convolve_scale_sse4.c
    cdef_block_sse4_1.c
    filterintra_sse4.c
    highbd_convolve_2d_sse4.c
    highbd_inv_txfm_sse4.c
    highbd_jnt_convolve_sse4.c
    highbd_txfm_utility_sse4.h
    reconinter_sse4.c
    selfguided_sse4.c
    warp_plane_sse4.c
    )

add_library(COMMON_ASM_SSE4_1 OBJECT ${all_files})
