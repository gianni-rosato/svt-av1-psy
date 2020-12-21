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

# Common/ASM_AVX512 Directory CMakeLists.txt

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

check_both_flags_add(
    -mavx512f
    -mavx512bw
    -mavx512dq
    -mavx512vl)

if(MSVC)
    check_both_flags_add(/arch:AVX2)
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
        # Intel Windows (*Note - The Warning level /W0 should be made to /W4 at some point)
        check_both_flags_add(/Qdiag-disable:10010,10148,10157 /W0)
    else()
        check_both_flags_add(-static-intel -w)
    endif()
endif()

file(GLOB all_files
    "*.h"
    "*.asm"
    "*.c")

add_library(ENCODER_ASM_AVX512 OBJECT ${all_files})
