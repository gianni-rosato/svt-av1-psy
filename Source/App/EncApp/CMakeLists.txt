#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# APP Directory CMakeLists.txt
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

# Include Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/Source/API/)
include_directories (${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

file(GLOB all_files
    "*.h"
    "../API/*.h"
    "*.c")

# App Source Files
add_executable (SvtAv1EncApp
    ${all_files}
)

#********** SET COMPILE FLAGS************

if(BUILD_SHARED_LIBS)
    set (SvtAv1EncLib SvtAv1EncShared)
else()
    set (SvtAv1EncLib SvtAv1EncStatic)
endif()

if (UNIX)
    # Link the Encoder App
     target_link_libraries (SvtAv1EncApp
        ${SvtAv1EncLib}
        pthread
        m)
endif()

if (MSVC OR MSYS OR MINGW OR WIN32)

    # Link the Encoder App
    target_link_libraries (SvtAv1EncApp
                           ${SvtAv1EncLib})
endif()

install(TARGETS SvtAv1EncApp RUNTIME DESTINATION bin)
