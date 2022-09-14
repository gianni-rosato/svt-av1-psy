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
include_directories(../../../API
        ../../Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
     )

check_both_flags_add(-msse2)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel" AND NOT WIN32)
    check_both_flags_add(-static-intel -w)
endif()

set(all_files
    EbComputeMean_Intrinsic_SSE2.c
    EbComputeMean_SSE2.h
    EbMeSadCalculation_Intrinsic_SSE2.c
    EbVariance_SSE2.h
    encodetxb_sse2.c
    fft_sse2.c
    highbd_variance_sse2.c
    variance_sse2.c
    )

set(asm_files
    highbd_variance_impl_sse2.asm
    subpel_variance_sse2.asm
    x64RegisterUtil.asm)

add_library(ENCODER_ASM_SSE2 OBJECT ${all_files})

asm_compile_to_target(ENCODER_ASM_SSE2 ${asm_files})
