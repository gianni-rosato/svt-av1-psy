#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
# Encoder Library directory CMakeLists.txt

# Shared Encoder Version
set(ENC_VERSION_MAJOR 0)
set(ENC_VERSION_MINOR 8)
set(ENC_VERSION_PATCH 2)
set(ENC_VERSION ${ENC_VERSION_MAJOR}.${ENC_VERSION_MINOR}.${ENC_VERSION_PATCH})

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
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

link_directories(
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Globals/
)

add_subdirectory(C_DEFAULT)
add_subdirectory(ASM_SSE2)
add_subdirectory(ASM_SSSE3)
add_subdirectory(ASM_SSE4_1)
add_subdirectory(ASM_AVX2)
add_subdirectory(ASM_AVX512)
add_subdirectory(Codec)
add_subdirectory(Globals)

#file(GLOB all_files
#    "Codec/*.h"
#    "Codec/*.c")

file(GLOB all_files
    "*.h"
    "*.c")

if(UNIX)
    if(NOT APPLE)
        find_library(M_LIB name m)
        if(M_LIB)
            list(APPEND PLATFORM_LIBS m)
        endif()
    endif()
    set(LIBS_PRIVATE "-lpthread -lm")
endif()

if(COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        set(LIBS_PRIVATE "-lgcov ${LIBS_PRIVATE}")
    else(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(LIBS_PRIVATE "--coverage ${LIBS_PRIVATE}")
    endif()
endif()

# Encoder Lib Source Files
add_library(SvtAv1Enc
    ${all_files}
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:FASTFEAT>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    $<TARGET_OBJECTS:COMMON_ASM_AVX512>
    $<TARGET_OBJECTS:ENCODER_GLOBALS>
    $<TARGET_OBJECTS:ENCODER_CODEC>
    $<TARGET_OBJECTS:ENCODER_C_DEFAULT>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE2>
    $<TARGET_OBJECTS:ENCODER_ASM_SSSE3>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE4_1>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX2>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX512>)
else ()
# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

link_directories(
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Globals/
)

add_subdirectory(C_DEFAULT)
add_subdirectory(Codec)
add_subdirectory(Globals)

#file(GLOB all_files
#    "Codec/*.h"
#    "Codec/*.c")

file(GLOB all_files
    "*.h"
    "*.c")

if(UNIX)
    if(NOT APPLE)
        find_library(M_LIB name m)
        if(M_LIB)
            list(APPEND PLATFORM_LIBS m)
        endif()
    endif()
    set(LIBS_PRIVATE "-lpthread -lm")
endif()

if(COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        set(LIBS_PRIVATE "-lgcov ${LIBS_PRIVATE}")
    else(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(LIBS_PRIVATE "--coverage ${LIBS_PRIVATE}")
    endif()
endif()

# Encoder Lib Source Files
add_library(SvtAv1Enc
    ${all_files}
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:FASTFEAT>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:ENCODER_GLOBALS>
    $<TARGET_OBJECTS:ENCODER_CODEC>
    $<TARGET_OBJECTS:ENCODER_C_DEFAULT>)
endif ()

set_target_properties(SvtAv1Enc PROPERTIES VERSION ${ENC_VERSION})
set_target_properties(SvtAv1Enc PROPERTIES SOVERSION ${ENC_VERSION_MAJOR})
target_link_libraries(SvtAv1Enc ${PLATFORM_LIBS})
install(TARGETS SvtAv1Enc DESTINATION "${CMAKE_INSTALL_LIBDIR}")

configure_file(pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
