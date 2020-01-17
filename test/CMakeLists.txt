#
# Copyright(c) 2019 Netflix, Inc.
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# APP Directory CMakeLists.txt

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/test/
    ${PROJECT_SOURCE_DIR}/test/ref
    ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/
    ${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src
    ${PROJECT_SOURCE_DIR}/Source/API
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec
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
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Globals
    ${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/App/EncApp
    ${PROJECT_SOURCE_DIR}/Source/API
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

# Define helper functions and macros used by Google Test.
include(../third_party/googletest/cmake/internal_utils.cmake)

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_ALLOW_KEYWORD_MACROS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALLOW_KEYWORD_MACROS")
endif()

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(gtest_force_shared_crt
    "Use shared (DLL) run-time lib even when Google Test is built as static lib."
    OFF)

add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1)

enable_testing()

if(UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()

file(GLOB all_files
    "*.h"
    "*.c"
    "*.cc"
    "ref/*.h"
    "ref/*.cc"
    "../Source/Lib/Encoder/Codec/*.c"
    "../Source/Lib/Decoder/Codec/EbDecBitReader.c"
    "../Source/Lib/Decoder/Codec/EbDecBitstreamUnit.c")

set(lib_list
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:FASTFEAT>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    $<TARGET_OBJECTS:COMMON_ASM_AVX512>
    $<TARGET_OBJECTS:ENCODER_C_DEFAULT>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE2>
    $<TARGET_OBJECTS:ENCODER_ASM_SSSE3>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE4_1>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX2>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX512>
    $<TARGET_OBJECTS:ENCODER_GLOBALS>
    gtest_all)
if(UNIX)
  # App Source Files
    add_executable(SvtAv1UnitTests
      ${all_files})
    # Link the Encoder App
    target_link_libraries(SvtAv1UnitTests
        ${lib_list}
        pthread
        m)
else()
    cxx_executable_with_flags(SvtAv1UnitTests
        "${cxx_default}"
        "${lib_list}"
        ${all_files})

    set_target_properties(SvtAv1UnitTests
        PROPERTIES
        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif()

install(TARGETS SvtAv1UnitTests RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(SvtAv1UnitTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1UnitTests)

add_subdirectory(api_test)
add_subdirectory(e2e_test)
