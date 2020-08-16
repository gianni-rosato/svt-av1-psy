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

# C_DEFAULT Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(../../../API
        ../../Encoder/Codec
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
        ${PROJECT_SOURCE_DIR}/third_party/fastfeat/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Globals/
        )

if(X86)
# Include Encoder Subdirectories
include_directories(
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/
        )
endif ()
file(GLOB all_files
    "*.h"
    "*.c")

add_library(ENCODER_CODEC OBJECT ${all_files})



