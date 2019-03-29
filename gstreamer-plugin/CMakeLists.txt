cmake_minimum_required(VERSION 2.8)
project(gst-svt-av1 C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

FIND_PACKAGE(PkgConfig REQUIRED)
pkg_check_modules(GSTREAMER REQUIRED gstreamer-1.0>=1.8)
include_directories(${GSTREAMER_INCLUDE_DIRS})

pkg_check_modules(GSTREAMER_BASE REQUIRED gstreamer-base-1.0>=1.8)
include_directories(${GSTREAMER_BASE_INCLUDE_DIRS})

pkg_check_modules(GSTREAMER_VIDEO REQUIRED gstreamer-video-1.0>=1.8)
include_directories(${GSTREAMER_VIDEO_INCLUDE_DIRS})

pkg_check_modules(SVT_AV1 REQUIRED SvtAv1Enc)
include_directories (${SVT_AV1_INCLUDE_DIRS})

set(CMAKE_C_FLAGS "-fstack-protector-strong -fPIE -fPIC -D_FORTIFY_SOURCE=2 -O2 -flto -Wall -Wextra -Wformat -Wformat-security")
set(CMAKE_SHARED_LINKER_FLAGS "-z noexecstack -z relro -z now -pie ")

add_library(gstsvtav1enc SHARED gstsvtav1enc.c)
target_link_libraries(gstsvtav1enc
    ${GSTREAMER_LIBRARIES}
    ${GSTREAMER_BASE_LIBRARIES}
    ${GSTREAMER_VIDEO_LIBRARIES}
    ${SVT_AV1_LIBRARIES})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${GSTREAMER_LIBDIR}")
    message(STATUS "   Install to: " ${GSTREAMER_LIBDIR})
endif()

install(TARGETS gstsvtav1enc
    LIBRARY DESTINATION gstreamer-1.0)

message(STATUS "   Build: " ${CMAKE_BUILD_TYPE})

