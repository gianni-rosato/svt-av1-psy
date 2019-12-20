#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# APP Directory CMakeLists.txt
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

file(GLOB all_files
    "*.h"
    "../API/*.h"
    "*.c")

# App Source Files
add_executable(SvtAv1DecApp
    ${all_files})

#********** SET COMPILE FLAGS************
# Link the Decoder App
target_link_libraries(SvtAv1DecApp
    SvtAv1Dec SvtAv1Enc)

if(UNIX)
    target_link_libraries(SvtAv1DecApp
        pthread
        m)
endif()

install(TARGETS SvtAv1DecApp RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
