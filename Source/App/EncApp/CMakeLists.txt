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

if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_LIBDOVI libdovi QUIET)
endif()

find_library(LIBDOVI_LIBRARY NAMES dovi libdovi
                             PATHS ${PC_LIBDOVI_LIBDIR}
)
find_path(LIBDOVI_INCLUDE_DIR NAMES libdovi/rpu_parser.h
                              PATHS ${PC_LIBDOVI_INCLUDEDIR})

set(all_files
    ../../API/EbDebugMacros.h
    ../../API/EbSvtAv1.h
    ../../API/EbSvtAv1Dec.h
    ../../API/EbSvtAv1Enc.h
    ../../API/EbSvtAv1ErrorCodes.h
    ../../API/EbSvtAv1ExtFrameBuf.h
    ../../API/EbSvtAv1Formats.h
    ../../API/EbSvtAv1Metadata.h
    EbAppConfig.c
    EbAppConfig.h
    EbAppContext.c
    EbAppContext.h
    EbAppInputy4m.c
    EbAppInputy4m.h
    EbAppMain.c
    EbAppOutputivf.c
    EbAppOutputivf.h
    EbAppProcessCmd.c
    EbTime.c
    EbTime.h
    )

if(LIBDOVI_FOUND)
  set(LIBDOVI_INCLUDE_DIRS ${LIBDOVI_INCLUDE_DIR})
  set(LIBDOVI_LIBRARIES ${LIBDOVI_LIBRARY})
endif()

# App Source Files
add_executable(SvtAv1EncApp
    ${all_files})

if(TARGET safestringlib)
    target_sources(SvtAv1EncApp PRIVATE $<TARGET_OBJECTS:safestringlib>)
endif()

#********** SET COMPILE FLAGS************
# Link the Encoder App
target_link_libraries(SvtAv1EncApp
        SvtAv1Enc
        ${LIBDOVI_LIBRARY})

if(UNIX)
    target_link_libraries(SvtAv1EncApp
        pthread
        m)
endif()

install(TARGETS SvtAv1EncApp RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
