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

# Globals Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(../../../API
        ${PROJECT_BINARY_DIR}/Source/Lib/Common/Codec/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
        ${PROJECT_SOURCE_DIR}/third_party/fastfeat/
        )

if(X86)
# Include Encoder Subdirectories
include_directories(
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
        )
endif()

file(GLOB all_files
    "*.h"
    "*.c")

add_library(ENCODER_GLOBALS OBJECT ${all_files})
add_dependencies(ENCODER_GLOBALS EbVersionHeaderGen)
