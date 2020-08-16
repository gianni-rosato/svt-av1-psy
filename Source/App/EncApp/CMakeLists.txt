#
# Copyright(c) 2019 Intel Corporation
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
#

# APP Directory CMakeLists.txt

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

file(GLOB all_files
    "*.h"
    "../API/*.h"
    "*.c")

# App Source Files
add_executable(SvtAv1EncApp
    ${all_files})

#********** SET COMPILE FLAGS************
# Link the Encoder App
target_link_libraries(SvtAv1EncApp
        SvtAv1Enc)
if(UNIX)
    target_link_libraries(SvtAv1EncApp
        pthread
        m)
endif()

install(TARGETS SvtAv1EncApp RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
