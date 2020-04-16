#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
# Common/Codec Directory CMakeLists.txt
if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set (X86 TRUE)
else ()
    set (X86 FALSE)
endif ()

if(X86)
# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)
else ()
# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)
endif ()
file(GLOB all_files
    "*.h"
    "*.c")

add_library(COMMON_CODEC OBJECT ${all_files})
