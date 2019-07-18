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


if(WIN32)
    set(CMAKE_ASM_NASM_FLAGS "${CMAKE_ASM_NASM_FLAGS} -DWIN64")
else()
    set(CMAKE_ASM_NASM_FLAGS "${CMAKE_ASM_NASM_FLAGS} -DUNIX64")
endif()

if(UNIX)
    if(APPLE)
        set(CMAKE_MACOSX_RPATH 1)
        set(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
        set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
        set(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
        set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pie -z noexecstack -z relro -z now")
    endif()
endif()

set(flags_to_test
    -Wextra
    -Wformat
    -Wformat-security
    -fstack-protector-strong
    -fPIE
    -fPIC
    -D_FORTIFY_SOURCE=2)
if(NOT MSVC)
    list(INSERT flags_to_test 0 -Wall)
        option(NATIVE "Build for native performance (march=native)" ON)
        if(NATIVE)
        list(APPEND flags_to_test -march=native)
    endif()
else()
    list(APPEND flags_to_test /MP)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /Od")
endif()
set(release_flags_to_test
    -O2
    -O3)
set(debug_flags_to_test
    -g
    -O0)

foreach(flag ${flags_to_test};${release_flags_to_test};${debug_flags_to_test})
    string(REGEX REPLACE "[^A-Za-z0-9]" "_" flag_var "${flag}")
    set(test_c_flag "C_FLAG${flag_var}")
    set(test_cxx_flag "CXX_FLAG${flag_var}")
    check_c_compiler_flag(${flag} "${test_c_flag}")
    check_cxx_compiler_flag(${flag} "${test_cxx_flag}")
    if(${flag} IN_LIST release_flags_to_test)
        if(${test_c_flag})
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${flag}")
        endif()
        if(${test_cxx_flag})
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag}")
        endif()
    elseif(${flag} IN_LIST debug_flags_to_test)
        if(${test_c_flag})
            set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${flag}")
        endif()
        if(${test_cxx_flag})
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag}")
        endif()
    else()
        if(${test_c_flag})
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
        endif()
        if(${test_cxx_flag})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        endif()
    endif()
endforeach()

if(CMAKE_ASM_NASM_OBJECT_FORMAT MATCHES "win")
    set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -gcv8")
elseif(CMAKE_ASM_NASM_OBJECT_FORMAT MATCHES "macho")
    if(CMAKE_ASM_NASM_COMPILER MATCHES "nasm")
        set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -gdwarf")
    elseif(CMAKE_ASM_NASM_COMPILER MATCHES "yasm")
        set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -gnull")
    endif()
else()
    if(CMAKE_ASM_NASM_COMPILER MATCHES "nasm")
        set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -gdwarf")
    elseif(CMAKE_ASM_NASM_COMPILER MATCHES "yasm")
        set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -gdwarf2")
    endif()
endif()

# BUILD_SHARED_LIBS is a standard CMake variable, but we declare it here to
# make it prominent in the GUI.
if(WIN32 OR MSYS OR BUILD_TESTING)
    option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." ON)
else()
    option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)
endif()

message(STATUS "BUILD_SHARED_LIBS : " ${BUILD_SHARED_LIBS})

option(BUILD_TESTING "Build SvtAv1UnitTests, SvtAv1ApiTests, and SvtAv1E2ETests unit tests")
option(COVERAGE "Generate coverage report")
option(BUILD_APPS "Build Enc and Dec Apps" ON)

# Prepare for Coveralls
if(COVERAGE AND NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

include(CheckSymbolExists)
check_symbol_exists("_mm512_extracti64x4_epi64" "immintrin.h" HAS_AVX512)
if(NOT HAS_AVX512)
    add_definitions(-DNON_AVX512_SUPPORT)
endif()

# Add Subdirectories
add_subdirectory(Source/Lib/Common)
add_subdirectory(Source/Lib/Encoder)
add_subdirectory(Source/Lib/Decoder)
if(BUILD_APPS)
    add_subdirectory(Source/App/EncApp)
    add_subdirectory(Source/App/DecApp)
endif()
if(BUILD_TESTING)
    include(CTest)
    message(STATUS "Building UnitTests")
    add_subdirectory(test)
    add_subdirectory(third_party/googletest)
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/Source/API/ DESTINATION "${CMAKE_INSTALL_FULL_INCLUDEDIR}/svt-av1" FILES_MATCHING PATTERN "*.h")
