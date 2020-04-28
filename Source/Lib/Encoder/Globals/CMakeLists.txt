#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
