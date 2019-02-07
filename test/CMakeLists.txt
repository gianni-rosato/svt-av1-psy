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

enable_testing()

file(GLOB all_files
    "*.h"
    "*.cc")

# App Source Files
add_executable (SvtAv1UnitTests 
    ${all_files}
)

#********** SET COMPILE FLAGS************

if (UNIX)
    # Link the Encoder App
     target_link_libraries (SvtAv1UnitTests
        SvtAv1Enc
        gtest 
        gtest_main
        pthread
        m)
endif()

if (MSVC OR MSYS OR MINGW OR WIN32)
    
    # Link the Encoder App
    target_link_libraries (SvtAv1UnitTests
                           SvtAv1Enc
                           gtest 
                           gtest_main)
endif()

add_test(SvtAv1UnitTests ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/SvtAv1UnitTests)

install(TARGETS SvtAv1UnitTests RUNTIME DESTINATION bin)

