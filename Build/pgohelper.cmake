#!/usr/bin/cmake -P
# cmake -P Build/pgohelper.cmake $PWD/Build $PWD/objective-2-fast::$PWD/objective-1-fast $PWD/Bin/Release/SvtAv1EncApp

if(CMAKE_ARGC LESS 6)
    message(
        FATAL_ERROR
            "Usage: cmake -P ${CMAKE_ARGV2} build_dir /path/to/videofiles::additional_paths /path/to/SvtAv1EncApp"
    )
endif()

set(BUILD_DIRECTORY "${CMAKE_ARGV3}")
set(VIDEO_DIRECTORY "${CMAKE_ARGV4}")
set(SvtAv1EncApp "${CMAKE_ARGV5}")

if(NOT EXISTS ${SvtAv1EncApp})
    message(
        FATAL_ERROR
            "Can't run pgo if the binaries don't exist. Looked at ${SvtAv1EncApp}"
    )
endif()

# Delete any existing Clang profiling data.

file(GLOB OLD_FILES
    "${BUILD_DIRECTORY}/*.profraw"
    "${BUILD_DIRECTORY}/*.profdata"
)
if(OLD_FILES)
    file(REMOVE ${OLD_FILES})
endif()

unset(videofiles)
string(REPLACE "::" ";" VIDEO_DIRECTORY "${VIDEO_DIRECTORY}")

foreach(dir IN LISTS VIDEO_DIRECTORY)
    file(GLOB videofile ${dir}/*.y4m)
    list(APPEND videofiles ${videofile})
endforeach()

foreach(video IN LISTS videofiles)
    get_filename_component(videoname "${video}" NAME_WE)

    SET(ENCODING_COMMAND ${SvtAv1EncApp} -i ${video} -b "${BUILD_DIRECTORY}/${videoname}.ivf" --preset 2 --film-grain 8 --tune 0)
    list(JOIN ENCODING_COMMAND " " ENCODING_COMMAND_STR)
    message(
        STATUS
            "Running ${ENCODING_COMMAND_STR}"
    )
    execute_process(COMMAND ${ENCODING_COMMAND})
endforeach()
