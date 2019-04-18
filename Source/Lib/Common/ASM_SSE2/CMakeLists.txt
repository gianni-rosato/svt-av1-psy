# 
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 

# ASM_SSE2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/Source/API/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_C_FLAGS "/MP")    
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -static-intel -w")
endif()

file(GLOB all_files
    "*.h"
    "*.asm"
    "*.c")

add_library(COMMON_ASM_SSE2 OBJECT ${all_files})




