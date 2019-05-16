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
include_directories (${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src)
include_directories(${PROJECT_SOURCE_DIR}/Source/API)

# Define helper functions and macros used by Google Test.
include(../../third_party/googletest/cmake/internal_utils.cmake)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 -D_ALLOW_KEYWORD_MACROS")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 -D_ALLOW_KEYWORD_MACROS")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

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
    "*.cc")

if (UNIX)
        
  # App Source Files
    add_executable (SvtAv1ApiTests 
      ${all_files})

    # Link the Encoder App
     target_link_libraries (SvtAv1ApiTests
        SvtAv1Enc
        gtest_all 
        pthread
        m)
        
endif(UNIX)

if (MSVC OR MSYS OR MINGW OR WIN32)

    set (lib_list SvtAv1Enc gtest_all)
                        
    cxx_executable_with_flags(SvtAv1ApiTests "${cxx_default}"
      "${lib_list}" ${all_files})

    set_target_properties(SvtAv1ApiTests
                        PROPERTIES
                        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif(MSVC OR MSYS OR MINGW OR WIN32)

install(TARGETS SvtAv1ApiTests RUNTIME DESTINATION bin)

add_test(SvtAv1ApiTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1ApiTests)

