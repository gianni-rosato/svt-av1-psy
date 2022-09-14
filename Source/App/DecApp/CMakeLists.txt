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
include_directories(${PROJECT_SOURCE_DIR}/Source/API/)

set(all_files
    ../../API/EbDebugMacros.h
    ../../API/EbSvtAv1.h
    ../../API/EbSvtAv1Dec.h
    ../../API/EbSvtAv1Enc.h
    ../../API/EbSvtAv1ErrorCodes.h
    ../../API/EbSvtAv1ExtFrameBuf.h
    ../../API/EbSvtAv1Formats.h
    ../../API/EbSvtAv1Metadata.h
    EbDecAppMain.c
    EbDecParamParser.c
    EbDecParamParser.h
    EbDecTime.c
    EbDecTime.h
    EbFileUtils.c
    EbFileUtils.h
    EbMD5Utility.c
    EbMD5Utility.h
    )

# App Source Files
add_executable(SvtAv1DecApp
    ${all_files})

#********** SET COMPILE FLAGS************
# Link the Decoder App
target_link_libraries(SvtAv1DecApp
    SvtAv1Dec)

if(UNIX)
    target_link_libraries(SvtAv1DecApp
        pthread
        m)
endif()

install(TARGETS SvtAv1DecApp RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
