#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# Libraries Directory CMakeLists.txt

macro(test_apply_compiler_flags)
    foreach(cflag ${ARGV})
        string(REGEX REPLACE "[^A-Za-z0-9]" "_" cflag_var "${cflag}")
        set(test_c_flag "C_FLAG${cflag_var}")
        check_c_compiler_flag(${cflag} "${test_c_flag}")
        if(${test_c_flag})
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${cflag}")
        endif()
    endforeach()
    if(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
        message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}: ${CMAKE_C_FLAGS}")
    endif()
endmacro()
if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set (X86 TRUE)
else ()
    set (X86 FALSE)
endif ()

if(X86)
add_subdirectory(Codec)
add_subdirectory(C_DEFAULT)
add_subdirectory(ASM_SSE2)
add_subdirectory(ASM_SSSE3)
add_subdirectory(ASM_SSE4_1)
add_subdirectory(ASM_AVX2)
add_subdirectory(ASM_AVX512)
else ()
add_subdirectory(Codec)
add_subdirectory(C_DEFAULT)
endif ()
