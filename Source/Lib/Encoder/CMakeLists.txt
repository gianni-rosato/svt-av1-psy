#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
# Encoder Library directory CMakeLists.txt

# Shared Encoder Version
set(ENC_VERSION_MAJOR 0)
set(ENC_VERSION_MINOR 7)
set(ENC_VERSION_PATCH 5)
set(ENC_VERSION ${ENC_VERSION_MAJOR}.${ENC_VERSION_MINOR}.${ENC_VERSION_PATCH})

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/)

link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

file(GLOB all_files
    "Codec/*.h"
    "Codec/*.c")

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
    $<TARGET_OBJECTS:COMMON_ASM_AVX512>)
set_target_properties(SvtAv1Enc PROPERTIES VERSION ${ENC_VERSION})
set_target_properties(SvtAv1Enc PROPERTIES SOVERSION ${ENC_VERSION_MAJOR})
target_link_libraries(SvtAv1Enc ${PLATFORM_LIBS})
install(TARGETS SvtAv1Enc DESTINATION "${CMAKE_INSTALL_LIBDIR}")

configure_file(pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
