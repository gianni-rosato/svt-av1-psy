#
# Copyright(c) 2019 Netflix, Inc.
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# APP Directory CMakeLists.txt

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/test/
    ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/
    ${PROJECT_SOURCE_DIR}/third_party/googletest/include
    ${PROJECT_SOURCE_DIR}/third_party/googletest/src
    ${PROJECT_SOURCE_DIR}/Source/API)

file(GLOB all_files
    "*.h"
    "*.cc")

set(lib_list
    SvtAv1Enc
    gtest_all)

if(UNIX)
  # App Source Files
    add_executable(SvtAv1ApiTests
      ${all_files})

    # Link the Encoder App
     target_link_libraries(SvtAv1ApiTests
        ${lib_list}
        pthread
        m)

else()
    cxx_executable_with_flags(SvtAv1ApiTests
        "${cxx_default}"
        "${lib_list}"
        ${all_files})

    set_target_properties(SvtAv1ApiTests
        PROPERTIES
        COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif()

install(TARGETS SvtAv1ApiTests RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(SvtAv1ApiTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1ApiTests)
