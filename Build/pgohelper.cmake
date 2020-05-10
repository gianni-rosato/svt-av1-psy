#!/usr/bin/cmake -P
# cmake -P Build/pgohelper.cmake $PWD/Build $PWD/objective-2-fast::$PWD/objective-1-fast $PWD/Bin/Release/SvtAv1EncApp $PWD/Bin/Release/SvtAv1DecApp

if(CMAKE_ARGC LESS 7)
    message(
        FATAL_ERROR
            "Usage: cmake -P ${CMAKE_ARGV2} build_dir /path/to/videofiles::additional_paths /path/to/SvtAv1EncApp /path/to/SvtAv1DecApp"
    )
endif()

set(BUILD_DIRECTORY "${CMAKE_ARGV3}")
set(VIDEO_DIRECTORY "${CMAKE_ARGV4}")
set(SvtAv1EncApp "${CMAKE_ARGV5}")
set(SvtAv1DecApp "${CMAKE_ARGV6}")

if(NOT EXISTS ${SvtAv1EncApp} OR NOT EXISTS ${SvtAv1DecApp})
    message(
        FATAL_ERROR
            "Can't run pgo if the binaries don't exist. Looked at ${SvtAv1EncApp} and ${SvtAv1DecApp}"
    )
endif()

unset(videofiles)
string(REPLACE "::" ";" VIDEO_DIRECTORY "${VIDEO_DIRECTORY}")

foreach(dir IN LISTS VIDEO_DIRECTORY)
    file(GLOB videofile ${dir}/*.y4m)
    list(APPEND videofiles ${videofile})
endforeach()

foreach(video IN LISTS videofiles)
    get_filename_component(videoname "${video}" NAME_WE)

    message(
        STATUS
            "Running ${SvtAv1EncApp} -i ${video} -b ${BUILD_DIRECTORY}/${videoname}.ivf"
    )
    execute_process(COMMAND ${SvtAv1EncApp} -i ${video} -b
                            ${BUILD_DIRECTORY}/${videoname}.ivf)
    message(
        STATUS "Running ${SvtAv1DecApp} -i ${BUILD_DIRECTORY}/${videoname}.ivf")
    execute_process(COMMAND ${SvtAv1DecApp} -i
                            ${BUILD_DIRECTORY}/${videoname}.ivf)
endforeach()
