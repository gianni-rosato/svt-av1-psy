#
# Copyright(c) 2019 Netflix, Inc
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# APP Directory CMakeLists.txt

include(ProcessorCount)
set(SVT_AV1_E2E_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
include("${SVT_AV1_E2E_ROOT}/test_data_util.cmake")
include("${SVT_AV1_E2E_ROOT}/libaom_dependency.cmake")

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/test/
    ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/
    ${PROJECT_SOURCE_DIR}/third_party/googletest/include
    ${PROJECT_SOURCE_DIR}/third_party/googletest/src
    ${PROJECT_SOURCE_DIR}/third_party/aom/inc
    ${PROJECT_SOURCE_DIR}/Source/API
    ${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec)

# Link Subdirectories
if(APPLE)
    link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/mac)
elseif(UNIX)
    link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/linux)
else()
    link_directories(${PROJECT_SOURCE_DIR}/third_party/aom/lib/msvc)
endif()

file(GLOB all_files
    "*.h"
    "*.cc"
    "../../Source/App/EncApp/EbAppInputy4m.c"
    "../../Source/App/EncApp/EbAppString.c")

set(lib_list
    SvtAv1Enc
    SvtAv1Dec
    gtest_all
    #SDL2.lib
    aom)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_WARNINGS")
endif()

if(UNIX)
  # App Source Files
    add_executable(SvtAv1E2ETests
        ${all_files})

    if(APPLE)
        add_custom_command(TARGET SvtAv1E2ETests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/third_party/aom/lib/mac/* ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
    else()
        add_custom_command(TARGET SvtAv1E2ETests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/third_party/aom/lib/linux/* ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
    endif()

    # Link the Encoder App
     target_link_libraries(SvtAv1E2ETests
        ${lib_list}
        pthread
        m)
else()
    cxx_executable_with_flags(SvtAv1E2ETests
        "${cxx_default}"
        "${lib_list}"
        ${all_files})
    set_target_properties(SvtAv1E2ETests
        PROPERTIES
        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif()

install(TARGETS SvtAv1E2ETests RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(SvtAv1E2ETests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1E2ETests)
add_dependencies(SvtAv1E2ETests DepLibAom)

# testdata download targets.
set(TEST_DATA_LIST_FILE "${SVT_AV1_E2E_ROOT}/test_vector_list.txt")
make_test_data_lists("${TEST_DATA_LIST_FILE}"
    test_files
    test_file_checksums)
list(LENGTH test_files num_test_files)
list(LENGTH test_file_checksums num_test_file_checksums)

math(EXPR max_file_index "${num_test_files} - 1")
foreach(test_index RANGE ${max_file_index})
    list(GET test_files ${test_index} test_file)
    list(GET test_file_checksums ${test_index} test_file_checksum)
    add_custom_target(testdata_${test_file}
        COMMAND ${CMAKE_COMMAND}
            -DSVT_AV1_E2E_ROOT="${SVT_AV1_E2E_ROOT}"
            -DSVT_AV1_TEST_FILE="${test_file}"
            -DAVT_AV1_TEST_CHECKSUM=${test_file_checksum} -P
            "${SVT_AV1_E2E_ROOT}/test_data_download_worker.cmake")
    list(APPEND testdata_targets testdata_${test_file})
endforeach()

# Create a custom build target for running each test data download target.
add_custom_target(TestVectors)
add_dependencies(TestVectors ${testdata_targets})
