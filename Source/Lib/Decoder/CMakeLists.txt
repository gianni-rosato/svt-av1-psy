#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
# Decoder Library directory CMakeLists.txt

# Shared Decoder Version
set(DEC_VERSION_MAJOR 0)
set(DEC_VERSION_MINOR 8)
set(DEC_VERSION_PATCH 3)
set(DEC_VERSION ${DEC_VERSION_MAJOR}.${DEC_VERSION_MINOR}.${DEC_VERSION_PATCH})
if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set(X86 TRUE)
else()
    set(X86 FALSE)
endif()

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

include_directories(${PROJECT_BINARY_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)
link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

if(X86)
    # Include Decoder Subdirectories
    include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/)

    link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/)
endif()

# Required for cmake to be able to tell Xcode how to link all of the object files
if(CMAKE_GENERATOR STREQUAL "Xcode")
    file(WRITE ${PROJECT_BINARY_DIR}/dummy.c "")
endif()

file(GLOB all_files
    "Codec/*.h"
    "Codec/*.c")

get_target_property(common_lib_source common_lib INTERFACE_SOURCES)
if(common_lib_source)
    list(APPEND all_files ${common_lib_source} ${PROJECT_BINARY_DIR}/dummy.c)
endif()

# Decoder Lib Source Files
if(X86)
    add_library(SvtAv1Dec
        ${all_files}
        $<TARGET_OBJECTS:COMMON_CODEC>
        $<TARGET_OBJECTS:FASTFEAT>
        $<TARGET_OBJECTS:COMMON_C_DEFAULT>
        $<TARGET_OBJECTS:COMMON_ASM_SSE2>
        $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
        $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
        $<TARGET_OBJECTS:COMMON_ASM_AVX2>
        $<TARGET_OBJECTS:COMMON_ASM_AVX512>)
else()
    add_library(SvtAv1Dec
        ${all_files}
        $<TARGET_OBJECTS:COMMON_CODEC>
        $<TARGET_OBJECTS:FASTFEAT>
        $<TARGET_OBJECTS:COMMON_C_DEFAULT>)
endif()

set_target_properties(SvtAv1Dec PROPERTIES VERSION ${DEC_VERSION})
set_target_properties(SvtAv1Dec PROPERTIES SOVERSION ${DEC_VERSION_MAJOR})
add_dependencies(SvtAv1Dec EbVersionHeaderGen)
target_link_libraries(SvtAv1Dec ${PLATFORM_LIBS})
install(TARGETS SvtAv1Dec DESTINATION "${CMAKE_INSTALL_LIBDIR}")

configure_file(pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Dec.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Dec.pc DESTINATION "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig")
