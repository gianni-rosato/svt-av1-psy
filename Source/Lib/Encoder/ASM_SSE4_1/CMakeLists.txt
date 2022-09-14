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
include_directories(../../../API
        ../../Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/)

check_both_flags_add(-msse4.1)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel" AND NOT WIN32)
    check_both_flags_add(-static-intel -w)
endif()

set(all_files
    EbCdef_sse4.c
    EbComputeSAD_Intrinsic_SSE4_1.c
    EbTemporalFiltering_SSE4_1.c
    EbTemporalFiltering_constants.h
    av1_quantize_sse4_1.c
    av1_txfm1d_sse4.h
    av1_txfm_sse4.h
    corner_match_sse4.c
    encodetxb_sse4.c
    highbd_fwd_txfm_sse4.c
    obmc_variance_sse4.c
    pickrst_sse4.c
    )

add_library(ENCODER_ASM_SSE4_1 OBJECT ${all_files})
