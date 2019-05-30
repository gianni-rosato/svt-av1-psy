#
# Copyright(c) 2018 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

cmake_minimum_required(VERSION 3.5)
set(project_name "svt-av1")

if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(WARNING "Building in-source is highly not recommended\n"
                    "Please use the Build folder or create your own.")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_program(YASM_EXE yasm)
if(YASM_EXE AND NOT CMAKE_ASM_NASM_COMPILER MATCHES "yasm")
    set(CMAKE_ASM_NASM_COMPILER ${YASM_EXE})
    message(STATUS "Found YASM: ${YASM_EXE}")
endif()

project(${project_name} C CXX ASM_NASM)

include(GNUInstallDirs)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR "lib")
endif()

if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
    set(CMAKE_INSTALL_INCLUDEDIR "include")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CAN_USE_ASSEMBLER TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR ON)


    endif()
endif()

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
#if(USE_CODE_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
#    include(CodeCoverage)
#    setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
#endif()

if (BUILD_TESTING)
    option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." ON)
endif()

# Add Subdirectories
if (BUILD_TESTING)
    option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." ON)
endif()
add_subdirectory (Source/Lib/Common)
add_subdirectory (Source/Lib/Encoder)
add_subdirectory (Source/App/EncApp)
add_subdirectory (Source/Lib/Decoder)
add_subdirectory (Source/App/DecApp)
if (BUILD_TESTING)
    add_subdirectory (test)
    add_subdirectory (third_party/googletest)
endif()
