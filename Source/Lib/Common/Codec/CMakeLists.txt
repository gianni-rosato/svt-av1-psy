#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
# Common/Codec Directory CMakeLists.txt

add_custom_target(EbVersionHeaderGen
    ${CMAKE_COMMAND}
        -D INPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/EbVersion.h.in
        -D OUTPUT_FILE=EbVersion.h
        -D GIT_ROOT_DIR=${CMAKE_SOURCE_DIR}
        -D PACKAGE_VERSION_STRING=${PROJECT_VERSION}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/ConfigureGitVersion.cmake
    BYPRODUCTS EbVersion.h
    COMMENT "Generating version header"
    VERBATIM)

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
