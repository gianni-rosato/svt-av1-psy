# 
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 

# APP Directory CMakeLists.txt
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

# Include Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/test/)
include_directories (${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src)
include_directories(${PROJECT_SOURCE_DIR}/Source/API )
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Codec )

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
    "../Source/Lib/Codec/grainSynthesis.h"
    "../Source/Lib/Codec/grainSynthesis.c")

if (UNIX)
  # App Source Files
    add_executable (SvtAv1UnitTests 
      ${all_files})

    # Link the Encoder App
     target_link_libraries (SvtAv1UnitTests
        SvtAv1Enc
        gtest_all 
        pthread
        m)
endif(UNIX)

if (MSVC OR MSYS OR MINGW OR WIN32)
  
    set (lib_list SvtAv1Enc gtest_all)
    cxx_executable_with_flags(SvtAv1UnitTests "${cxx_default}"
      "${lib_list}" ${all_files})

    set_target_properties(SvtAv1UnitTests
                        PROPERTIES
                        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")

endif(MSVC OR MSYS OR MINGW OR WIN32)

install(TARGETS SvtAv1UnitTests RUNTIME DESTINATION bin)

add_test(SvtAv1UnitTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1UnitTests)
