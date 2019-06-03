# 
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 
# EncoderLib Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/Source/API/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/)


link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/Bin/${CMAKE_BUILD_TYPE}/)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   SET(CMAKE_C_FLAGS "/MP")    
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR "lib")
endif()

file(GLOB all_files
    "*.h"
    "*.c")
 
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    list(APPEND PLATFORM_LIBS m)
endif()

if(BUILD_SHARED_LIBS)
    # Encoder Lib Source Files
    add_library(SvtAv1EncShared SHARED
        ${all_files}
        $<TARGET_OBJECTS:COMMON_CODEC>
        $<TARGET_OBJECTS:COMMON_C_DEFAULT>
        $<TARGET_OBJECTS:COMMON_ASM_SSE2>
        $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
        $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
        $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    )
    set_target_properties(SvtAv1EncShared PROPERTIES OUTPUT_NAME SvtAv1Enc)
    target_link_libraries(SvtAv1EncShared ${PLATFORM_LIBS})
    install(TARGETS SvtAv1EncShared DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endif(BUILD_SHARED_LIBS)

add_library(SvtAv1EncStatic STATIC
    ${all_files}
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
)
set_target_properties(SvtAv1EncStatic PROPERTIES OUTPUT_NAME SvtAv1Enc)
target_link_libraries(SvtAv1EncStatic ${PLATFORM_LIBS})
install(TARGETS SvtAv1EncStatic DESTINATION "${CMAKE_INSTALL_LIBDIR}")

if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
    set(CMAKE_INSTALL_INCLUDEDIR "include/svt-av1")
endif()
    
configure_file(../pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Enc.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
install(DIRECTORY ../../../API/ DESTINATION "${CMAKE_INSTALL_PREFIX}/include/svt-av1" FILES_MATCHING PATTERN "*.h")
