#!/usr/bin/cmake -P
# cmake -P Build/obj1fastdownloader.cmake $PWD/Build

if(CMAKE_ARGC LESS 4)
    message(FATAL_ERROR "Usage: cmake -P ${CMAKE_ARGV2} build_dir")
endif()

set(BUILD_DIRECTORY ${CMAKE_ARGV3})

# We only need to encode one short video during profiling, in order to execute
# the code paths we want to optimise. Encoding 29 more files provides no
# benefits and slows down the process to an unacceptable level.
set(PROFILING_VIDEO "stefan_sif.y4m")

set(EXPECTED_HASH
    28c0dc1afdd9081d8004cc35f839f8d34cabc0fc8b677f9f70286be6bdbd9f57)
set(FILEHASH "")
if(EXISTS ${BUILD_DIRECTORY}/objective-1-fast/${PROFILING_VIDEO})
    file(SHA256 "${BUILD_DIRECTORY}/objective-1-fast/${PROFILING_VIDEO}" FILEHASH)
endif()

if(NOT FILEHASH STREQUAL EXPECTED_HASH)
    file(
        DOWNLOAD https://media.xiph.org/video/derf/y4m/${PROFILING_VIDEO}
        ${BUILD_DIRECTORY}/objective-1-fast/${PROFILING_VIDEO}
        SHOW_PROGRESS
        EXPECTED_HASH SHA256=${EXPECTED_HASH})
endif()
