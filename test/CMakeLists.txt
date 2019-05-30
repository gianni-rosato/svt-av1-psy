# 
# Copyright(c) 2019 Netflix, Inc.
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 

# APP Directory CMakeLists.txt
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

# Include Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/test/)
include_directories (${PROJECT_SOURCE_DIR}/test/ref)
include_directories (${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src)
include_directories(${PROJECT_SOURCE_DIR}/Source/API )
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec )
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec)
include_directories (${PROJECT_SOURCE_DIR}/Source/App/EncApp)
include_directories (${PROJECT_SOURCE_DIR}/Source/API)

# Define helper functions and macros used by Google Test.
include(../third_party/googletest/cmake/internal_utils.cmake)

# BUILD_SHARED_LIBS is a standard CMake variable, but we declare it here to
# make it prominent in the GUI.
if (UNIX)
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)
endif()

if (MSVC OR MSYS OR MINGW OR WIN32)
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." ON)
endif()

if(MSVC)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 -D_ALLOW_KEYWORD_MACROS")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 -D_ALLOW_KEYWORD_MACROS")
endif (MSVC)

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(
  gtest_force_shared_crt
  "Use shared (DLL) run-time lib even when Google Test is built as static lib."
  OFF)

add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1)

enable_testing()

file(GLOB all_files
    "*.h"
    "*.cc"
    "ref/*.h"
    "ref/*.cc"
    "../Source/Lib/Encoder/Codec/*.c"
    "../Source/Lib/Decoder/Codec/EbDecBitReader.c"
    "../Source/Lib/Decoder/Codec/EbDecBitstreamUnit.c"
    )

if (UNIX)
  # App Source Files
    add_executable (SvtAv1UnitTests
      ${all_files})

    # Link the Encoder App
     target_link_libraries (SvtAv1UnitTests
        $<TARGET_OBJECTS:COMMON_CODEC>
        $<TARGET_OBJECTS:COMMON_C_DEFAULT>
        $<TARGET_OBJECTS:COMMON_ASM_SSE2>
        $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
        $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
        $<TARGET_OBJECTS:COMMON_ASM_AVX2>
        gtest_all
        pthread
        m)
endif(UNIX)

if (MSVC OR MSYS OR MINGW OR WIN32)
    set (lib_list
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    gtest_all)
    cxx_executable_with_flags(SvtAv1UnitTests "${cxx_default}"
      "${lib_list}" ${all_files})

    set_target_properties(SvtAv1UnitTests
                        PROPERTIES
                        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")

endif(MSVC OR MSYS OR MINGW OR WIN32)

install(TARGETS SvtAv1UnitTests RUNTIME DESTINATION bin)

add_test(SvtAv1UnitTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1UnitTests)

add_subdirectory (api_test)
add_subdirectory (e2e_test)

