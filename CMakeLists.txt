# 
# Copyright(c) 2018 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 

cmake_minimum_required(VERSION 3.5)
set(project_name "svt-av1")

project(${project_name} C CXX ASM_NASM)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_STANDARD 99)
set(CAN_USE_ASSEMBLER TRUE)
set(CMAKE_NASM_LINK_EXECUTABLE "ld -o <TARGET> <OBJECTS>")
set(CMAKE_ASM_NASM_SOURCE_FILE_EXTENSIONS asm)


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_ASM_NASM_OBJECT_FORMAT macho64)
    set(CMAKE_ASM_NASM_FLAGS "-DUNIX64")
    set(CMAKE_MACOSX_RPATH 1)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_C_FLAGS "-Wall -Wextra -Wformat -Wformat-security -mavx -fstack-protector-strong -fPIE -fPIC -D_FORTIFY_SOURCE=2 ")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
	SET(CMAKE_EXE_LINKER_FLAGS "-z noexecstack -z relro -z now -pie ")
	SET(CMAKE_SHARED_LINKER_FLAGS "-z noexecstack -z relro -z now -pie ")
    set(CMAKE_ASM_NASM_FLAGS "-DUNIX64")
    set(CMAKE_ASM_NASM_FLAGS_DEBUG "-g null")
    set(CMAKE_ASM_NASM_OBJECT_FORMAT elf64)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_EXE_LINKER_FLAGS " ")
    set(CMAKE_SHARED_LINKER_FLAGS " ")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_ASM_NASM_FLAGS "-DWIN64")
    set(CMAKE_ASM_NASM_FLAGS_DEBUG "")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

include(CTest)

set(USE_CODE_COVERAGE 0)

# Prepare for Coveralls
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
if(USE_CODE_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
endif()

# Add Subdirectories
add_subdirectory (Source/Lib)
add_subdirectory (Source/App)
add_subdirectory (Source/SimpleApp)
if (BUILD_TESTING)
    add_subdirectory (test)
    add_subdirectory (third_party/googletest)
endif()
