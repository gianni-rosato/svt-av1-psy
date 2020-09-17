if(WIN32)
    file(WRITE dummy.c "")
    add_library(safestringlib OBJECT dummy.c)
else()
    add_library(safestringlib OBJECT
        safe_lib_errno.h
        safe_str_lib.h
        safe_lib.h
        safe_types.h
        safeclib_private.h
        strcpy_s.c
        safe_str_constraint.h
        safe_str_constraint.c
        ignore_handler_s.c
        strncpy_s.c
        strnlen_s.c)
    target_compile_definitions(safestringlib PUBLIC SAFECLIB_STR_NULL_SLACK=1)
endif()
