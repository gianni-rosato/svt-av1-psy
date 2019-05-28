# 
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
# 
# DecoderLib Directory CMakeLists.txt

# Include Decoder Subdirectories
include_directories (${PROJECT_SOURCE_DIR}/Source/API/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
include_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec/)


link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
link_directories (${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec/)

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
    # Decoder Lib Source Files
    add_library(SvtAv1DecShared SHARED
        ${all_files}
        $<TARGET_OBJECTS:COMMON_CODEC>
        $<TARGET_OBJECTS:COMMON_C_DEFAULT>
        $<TARGET_OBJECTS:COMMON_ASM_SSE2>
        $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
        $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
        $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    )
    set_target_properties(SvtAv1DecShared PROPERTIES OUTPUT_NAME SvtAv1Dec)
    target_link_libraries(SvtAv1DecShared ${PLATFORM_LIBS})
    install(TARGETS SvtAv1DecShared DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endif(BUILD_SHARED_LIBS)

add_library(SvtAv1DecStatic STATIC
    ${all_files}
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
)
set_target_properties(SvtAv1DecStatic PROPERTIES OUTPUT_NAME SvtAv1Dec)
target_link_libraries(SvtAv1DecStatic ${PLATFORM_LIBS})
install(TARGETS SvtAv1DecStatic DESTINATION "${CMAKE_INSTALL_LIBDIR}")

if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
    set(CMAKE_INSTALL_INCLUDEDIR "include/svt-av1")
endif()
    
configure_file(../pkg-config.pc.in ${CMAKE_BINARY_DIR}/SvtAv1Dec.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/SvtAv1Dec.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
install(DIRECTORY ../../../API/ DESTINATION "${CMAKE_INSTALL_PREFIX}/include/svt-av1" FILES_MATCHING PATTERN "*.h")
