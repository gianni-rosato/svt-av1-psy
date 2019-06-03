# 
# Copyright(c) 2019 Netflix, Inc
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 

include(ProcessorCount)
set(SVT_AV1_E2E_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
include("${SVT_AV1_E2E_ROOT}/test_data_util.cmake")
include("${SVT_AV1_E2E_ROOT}/libaom_dependency.cmake")

# APP Directory CMakeLists.txt
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

# Include Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/test/)
include_directories (${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/SDL2-2.0.9/include)
include_directories(${PROJECT_SOURCE_DIR}/third_party/aom/inc)
include_directories(${PROJECT_SOURCE_DIR}/Source/API)

# Link Subdirectories
if (MSVC OR MSYS OR MINGW OR WIN32)
link_directories(${PROJECT_SOURCE_DIR}/third_party/SDL2-2.0.9/lib/x64)
link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/msvc)
endif(MSVC OR MSYS OR MINGW OR WIN32)
if (UNIX AND NOT APPLE)
link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/linux)
endif(UNIX AND NOT APPLE)

if (APPLE)
link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/mac)
endif(APPLE)


# Define helper functions and macros used by Google Test.
include(../../third_party/googletest/cmake/internal_utils.cmake)

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
    "../../Source/App/EncApp/EbAppInputy4m.c"
    "../../Source/App/EncApp/EbAppString.c")

if (UNIX)
        
  # App Source Files
    add_executable (SvtAv1E2ETests 
      ${all_files})

    if (APPLE)
      add_custom_command(TARGET SvtAv1E2ETests POST_BUILD
                      COMMAND ${CMAKE_COMMAND} -E copy
                           ${PROJECT_SOURCE_DIR}/third_party/aom/lib/mac/* ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
    endif(APPLE)
    if (UNIX AND NOT APPLE)
      add_custom_command(TARGET SvtAv1E2ETests POST_BUILD
                      COMMAND ${CMAKE_COMMAND} -E copy
                           ${PROJECT_SOURCE_DIR}/third_party/aom/lib/linux/* ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
    endif(UNIX AND NOT APPLE)

    # Link the Encoder App
     target_link_libraries (SvtAv1E2ETests
        aom
        SvtAv1Enc
        gtest_all 
        pthread
        m)
        
endif(UNIX)

if (MSVC OR MSYS OR MINGW OR WIN32)

    set (lib_list SvtAv1Enc gtest_all)
                        
    cxx_executable_with_flags(SvtAv1E2ETests "${cxx_default}"
      "${lib_list}" ${all_files})
    target_link_libraries(SvtAv1E2ETests 
    #SDL2.lib 
    aom.lib)
    set_target_properties(SvtAv1E2ETests
                        PROPERTIES
                        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif(MSVC OR MSYS OR MINGW OR WIN32)

install(TARGETS SvtAv1E2ETests RUNTIME DESTINATION bin)

add_test(SvtAv1E2ETests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1E2ETests)
add_dependencies(SvtAv1E2ETests DepLibAom)

# testdata download targets.
set(TEST_DATA_LIST_FILE "${SVT_AV1_E2E_ROOT}/test_vector_list.txt")
make_test_data_lists("${TEST_DATA_LIST_FILE}" test_files
                         test_file_checksums)
list(LENGTH test_files num_test_files)
list(LENGTH test_file_checksums num_test_file_checksums)

math(EXPR max_file_index "${num_test_files} - 1")
foreach(test_index RANGE ${max_file_index})
    list(GET test_files ${test_index} test_file)
    list(GET test_file_checksums ${test_index} test_file_checksum)
    add_custom_target(testdata_${test_file}
                        COMMAND
                          ${CMAKE_COMMAND}
                          -DSVT_AV1_E2E_ROOT="${SVT_AV1_E2E_ROOT}"
                          -DSVT_AV1_TEST_FILE="${test_file}"
                          -DAVT_AV1_TEST_CHECKSUM=${test_file_checksum} -P
                          "${SVT_AV1_E2E_ROOT}/test_data_download_worker.cmake")
    list(APPEND testdata_targets testdata_${test_file})
endforeach()

# Create a custom build target for running each test data download target.
add_custom_target(TestVectors)
add_dependencies(TestVectors ${testdata_targets})
