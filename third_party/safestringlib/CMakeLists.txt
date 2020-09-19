if(WIN32)
    if(HAVE_STRNLEN_S)
        file(WRITE dummy.c "")
        add_library(safestringlib OBJECT dummy.c)
    else()
        add_library(safestringlib OBJECT
            safeclib_private.h
            safe_lib.h
            safe_types.h
            safe_lib_errno.h
            safe_str_lib.h
            safe_str_constraint.h
            strnlen_s.c
            safe_str_constraint.c
            ignore_handler_s.c)
    endif()
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
endif()

target_compile_definitions(safestringlib PUBLIC SAFECLIB_STR_NULL_SLACK=1)
