#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# Common/ASM_AVX512 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/)
link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)

set(flags_to_test
    -mavx512f
    -mavx512bw
    -mavx512dq
    -mavx512vl
    -static-intel
    /Qwd10148
    /Qwd10010
    /Qwd10157)

foreach(cflag ${flags_to_test})
    string(REGEX REPLACE "[^A-Za-z0-9]" "_" cflag_var "${cflag}")
    set(test_c_flag "C_FLAG${cflag_var}")
    check_c_compiler_flag(${cflag} "${test_c_flag}")
    if(${test_c_flag})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${cflag}")
    endif()
endforeach()

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
        # Intel Windows (*Note - The Warning level /W0 should be made to /W4 at some point)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    endif()
endif()

file(GLOB all_files
    "*.h"
    "*.asm"
    "*.c")

add_library(COMMON_ASM_AVX512 OBJECT ${all_files})
