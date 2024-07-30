#
# Copyright(c) 2019 Intel Corporation
#
# This source code is subject to the terms of the BSD 3 Clause Clear License and
# the Alliance for Open Media Patent License 1.0. If the BSD 3 Clause Clear License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at https://www.aomedia.org/license/license. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
#

# Encoder Library directory CMakeLists.txt

# Shared Encoder Version
set(ENC_VERSION_MAJOR 2)
set(ENC_VERSION_MINOR 2)
set(ENC_VERSION_PATCH 0)
set(ENC_VERSION ${ENC_VERSION_MAJOR}.${ENC_VERSION_MINOR}.${ENC_VERSION_PATCH})

add_library(common_lib INTERFACE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Android")
    find_library(M_LIB NAMES m)
    if(M_LIB)
        list(APPEND PLATFORM_LIBS m)
    endif()
    list(APPEND PLATFORM_LIBS Threads::Threads)
elseif(UNIX)
    if(NOT APPLE)
        find_library(M_LIB NAMES m)
        if(M_LIB)
            list(APPEND PLATFORM_LIBS m)
        endif()
        list(APPEND PLATFORM_LIBS Threads::Threads)
    endif()
    set(LIBS_PRIVATE "-lpthread -lm")
endif()

if(COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        set(LIBS_PRIVATE "-lgcov ${LIBS_PRIVATE}")
    else(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(LIBS_PRIVATE "--coverage ${LIBS_PRIVATE}")
    endif()
endif()

include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

add_library(SvtAv1Enc)
# Required for cmake to be able to tell Xcode how to link all of the object files
if(CMAKE_GENERATOR STREQUAL "Xcode")
    file(WRITE ${PROJECT_BINARY_DIR}/dummy.c "")
    target_sources(SvtAv1Enc PRIVATE ${PROJECT_BINARY_DIR}/dummy.c)
endif()

add_subdirectory(C_DEFAULT)
add_subdirectory(Codec)
add_subdirectory(Globals)
if(NOT COMPILE_C_ONLY AND HAVE_X86_PLATFORM)
    # Include Encoder Subdirectories
    target_include_directories(SvtAv1Enc PRIVATE
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_SSE2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_SSSE3/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_SSE4_1/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_AVX2/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_AVX512/)
    add_subdirectory(ASM_SSE2)
    add_subdirectory(ASM_SSSE3)
    add_subdirectory(ASM_SSE4_1)
    add_subdirectory(ASM_AVX2)
    if(ENABLE_AVX512)
        add_subdirectory(ASM_AVX512)
    endif()
elseif(NOT COMPILE_C_ONLY AND HAVE_ARM_PLATFORM)
    target_include_directories(SvtAv1Enc PRIVATE
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_NEON/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_NEON_DOTPROD/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_NEON_I8MM/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_SVE/
        ${PROJECT_SOURCE_DIR}/Source/Lib/ASM_SVE2/)
    add_subdirectory(ASM_NEON)
    if(ENABLE_NEON_DOTPROD)
        add_subdirectory(ASM_NEON_DOTPROD)
    endif()
    if(ENABLE_NEON_I8MM)
        add_subdirectory(ASM_NEON_I8MM)
    endif()
    if(ENABLE_SVE)
      add_subdirectory(ASM_SVE)
    endif()
    if(ENABLE_SVE2)
      add_subdirectory(ASM_SVE2)
    endif()
endif()

# Encoder Lib Source Files
target_sources(SvtAv1Enc PRIVATE
    $<TARGET_OBJECTS:FASTFEAT>
    $<TARGET_OBJECTS:GLOBALS>
    $<TARGET_OBJECTS:CODEC>
    $<TARGET_OBJECTS:C_DEFAULT>)
if(NOT COMPILE_C_ONLY AND HAVE_X86_PLATFORM)
    target_sources(SvtAv1Enc PRIVATE
        $<TARGET_OBJECTS:ASM_SSE2>
        $<TARGET_OBJECTS:ASM_SSSE3>
        $<TARGET_OBJECTS:ASM_SSE4_1>
        $<TARGET_OBJECTS:ASM_AVX2>)
    if(ENABLE_AVX512)
        target_sources(SvtAv1Enc PRIVATE
            $<TARGET_OBJECTS:ASM_AVX512>)
    endif()
elseif(NOT COMPILE_C_ONLY AND HAVE_ARM_PLATFORM)
    target_sources(SvtAv1Enc PRIVATE $<TARGET_OBJECTS:ASM_NEON>)
    if(ENABLE_NEON_DOTPROD)
        target_sources(SvtAv1Enc PRIVATE $<TARGET_OBJECTS:ASM_NEON_DOTPROD>)
    endif()
    if(ENABLE_NEON_I8MM)
        target_sources(SvtAv1Enc PRIVATE $<TARGET_OBJECTS:ASM_NEON_I8MM>)
    endif()
    if(ENABLE_SVE)
        target_sources(SvtAv1Enc PRIVATE $<TARGET_OBJECTS:ASM_SVE>)
    endif()
    if(ENABLE_SVE2)
        target_sources(SvtAv1Enc PRIVATE $<TARGET_OBJECTS:ASM_SVE2>)
    endif()
endif()

get_target_property(common_lib_source common_lib INTERFACE_SOURCES)
if(common_lib_source)
    target_sources(SvtAv1Enc PRIVATE ${common_lib_source})
endif()

set_target_properties(SvtAv1Enc PROPERTIES VERSION ${ENC_VERSION})
set_target_properties(SvtAv1Enc PROPERTIES SOVERSION ${ENC_VERSION_MAJOR})
set_target_properties(SvtAv1Enc PROPERTIES C_VISIBILITY_PRESET hidden)
target_link_libraries(SvtAv1Enc PUBLIC ${PLATFORM_LIBS})
if(USE_EXTERNAL_CPUINFO)
    target_link_libraries(SvtAv1Enc PRIVATE cpuinfo::cpuinfo)
elseif(NOT COMPILE_C_ONLY AND HAVE_X86_PLATFORM)
    target_link_libraries(SvtAv1Enc PRIVATE cpuinfo_public)
endif()

set(ENC_PKG_CONFIG_EXTRA_CFLAGS "")
if (BUILD_SHARED_LIBS)
    target_compile_definitions(SvtAv1Enc
        PUBLIC EB_DLL
        PRIVATE EB_BUILDING_SHARED_LIBS)
    set(ENC_PKG_CONFIG_EXTRA_CFLAGS " -DEB_DLL")
endif()
install(TARGETS SvtAv1Enc
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

configure_file(pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
