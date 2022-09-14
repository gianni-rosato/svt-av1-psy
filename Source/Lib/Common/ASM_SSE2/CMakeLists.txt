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

# ASM_SSE2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)

check_both_flags_add(-msse2)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel" AND NOT WIN32)
    check_both_flags_add(-static-intel -w)
endif()

set(all_files
    EbAvcStyleMcp_Intrinsic_SSE2.c
    EbDeblockingFilter_Intrinsic_SSE2.c
    EbDeblockingFilter_SSE2.h
    EbHighbdIntraPrediction_SSE2.h
    EbIntraPrediction_AV1_Intrinsic_SSE2.c
    EbMcp_SSE2.h
    EbPackUnPack_Intrinsic_SSE2.c
    EbPictureOperators_Intrinsic_SSE2.c
    EbPictureOperators_SSE2.h
    av1_txfm_sse2.h
    convolve_2d_sse2.c
    convolve_sse2.c
    highbd_intrapred_sse2.c
    highbd_subtract_sse2.c
    jnt_convolve_2d_sse2.c
    jnt_convolve_sse2.c
    lpf_common_sse2.h
    synonyms.h
    transpose_sse2.h
    wiener_convolve_sse2.c
    )

set(asm_files
    intrapred_sse2.asm
    EbPictureOperators_SSE2.asm
    highbd_intrapred_sse2_.asm
    aom_subpixel_8t_sse2.asm
    subtract_sse2.asm)

add_library(COMMON_ASM_SSE2 OBJECT ${all_files})

asm_compile_to_target(COMMON_ASM_SSE2 ${asm_files})
