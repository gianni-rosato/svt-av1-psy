#!/usr/bin/cmake -P
# cmake -P Build/obj1fastdownloader.cmake $PWD/Build

if(CMAKE_ARGC LESS 4)
    message(FATAL_ERROR "Usage: cmake -P ${CMAKE_ARGV2} build_dir")
endif()

set(BUILD_DIRECTORY ${CMAKE_ARGV3})

set(EXPECTED_HASH
    a4aab832a5dd67eb42687e2935c1044b87eabe0a96fc73550339150a8f05c5b1)
set(TARHASH "")
if(EXISTS ${BUILD_DIRECTORY}/objective-1-fast.tar.gz)
    file(SHA256 "${BUILD_DIRECTORY}/objective-1-fast.tar.gz" TARHASH)
endif()

if(NOT TARHASH STREQUAL EXPECTED_HASH)
    file(
        DOWNLOAD https://media.xiph.org/video/derf/objective-1-fast.tar.gz
        ${BUILD_DIRECTORY}/objective-1-fast.tar.gz
        SHOW_PROGRESS
        EXPECTED_HASH SHA256=${EXPECTED_HASH})
endif()
if(NOT EXISTS ${BUILD_DIRECTORY}/TAR_DONE
   OR NOT EXISTS ${BUILD_DIRECTORY}/objective-1-fast/DOTA2_60f_420.y4m)
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf objective-1-fast.tar.gz
                    WORKING_DIRECTORY ${BUILD_DIRECTORY})
    execute_process(COMMAND ${CMAKE_COMMAND} -E touch TAR_DONE
                    WORKING_DIRECTORY ${BUILD_DIRECTORY})
endif()
