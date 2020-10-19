file(WRITE ${CMAKE_BINARY_DIR}/conftest.c [=[
#include <string.h>
#include <stdint.h>
int main() {
    return (uintptr_t)strnlen_s;
}
]=])
try_compile(HAVE_STRNLEN_S ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/conftest.c)

file(WRITE ${CMAKE_BINARY_DIR}/conftest.c [=[
#include <string.h>
#include <stdint.h>
int main() {
    return (uintptr_t)strcpy_s;
}
]=])
try_compile(HAVE_STRNCPY_S ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/conftest.c)

file(WRITE ${CMAKE_BINARY_DIR}/conftest.c [=[
#include <string.h>
#include <stdint.h>
int main() {
    return (uintptr_t)strcpy_s;
}
]=])
try_compile(HAVE_STRCPY_S ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/conftest.c)

file(WRITE ${CMAKE_BINARY_DIR}/dummy.c "")
add_library(safestringlib OBJECT ${CMAKE_BINARY_DIR}/dummy.c)
if(NOT HAVE_STRCPY_S OR NOT HAVE_STRNCPY_S OR NOT HAVE_STRNLEN_S)
    target_sources(safestringlib PRIVATE
        safeclib_private.h
        safe_lib.h
        safe_types.h
        safe_lib_errno.h
        safe_str_lib.h
        safe_str_constraint.h
        safe_str_constraint.c
        ignore_handler_s.c
        $<$<NOT:$<BOOL:${HAVE_STRCPY_S}>>:strcpy_s.c>
        $<$<NOT:$<BOOL:${HAVE_STRNCPY_S}>>:strncpy_s.c>
        $<$<NOT:$<BOOL:${HAVE_STRNLEN_S}>>:strnlen_s.c>)
endif()

target_compile_definitions(safestringlib PUBLIC
    SAFECLIB_STR_NULL_SLACK=1
    $<$<BOOL:${HAVE_STRCPY_S}>:HAVE_STRCPY_S=1>
    $<$<BOOL:${HAVE_STRNCPY_S}>:HAVE_STRNCPY_S=1>
    $<$<BOOL:${HAVE_STRNLEN_S}>:HAVE_STRNLEN_S=1>
    $<$<BOOL:${WIN32}>:_WIN32_WINNT=0x0601>)

set(HAVE_STRCPY_S PARENT_SCOPE)
set(HAVE_STRNCPY_S PARENT_SCOPE)
set(HAVE_STRNLEN_S PARENT_SCOPE)
