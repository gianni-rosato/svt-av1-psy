include(ExternalProject)
#/property:Configuration=${CMAKE_BUILD_TYPE} ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/
if(APPLE)
    set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/mac)
    set(TARGET "libaom*dylib")
    install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.dylib DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.0.dylib DESTINATION ${CMAKE_INSTALL_LIBDIR})
elseif(UNIX)
    set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/linux)
    set(TARGET "libaom.so*")
    install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.so DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.so.0 DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
    set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/msvc)
    set(CUSTOM_POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/Release/
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/Release/
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/MinSizeRel/
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/MinSizeRel/
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/RelWithDebInfo/
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/RelWithDebInfo/
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/Debug/
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/Debug/aom.lib ${TARGET_OUTPUT_PATH}/Debug/)
endif()

if(UNIX)
    set(CUSTOM_CONFIG "-DBUILD_SHARED_LIBS=1")
    set(CUSTOM_POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/${TARGET} ${TARGET_OUTPUT_PATH})
else()
    set(AOM_BUILD_COMMAND --config MinSizeRel COMMAND ${CMAKE_COMMAND} --build . --target aom --config Debug)
endif()

ExternalProject_Add(DepLibAom
    PREFIX "${CMAKE_BINARY_DIR}/libaom"
    GIT_REPOSITORY "https://github.com/Cidana-Developers/aom.git"
    GIT_TAG av1-normative
    GIT_SHALLOW 1
    CMAKE_ARGS
        ${CUSTOM_CONFIG}
        -DCONFIG_INSPECTION=1
        -DENABLE_TESTS=0
        -DCONFIG_AV1_ENCODER=0
        -DENABLE_DOCS=0
        -DENABLE_EXAMPLES=0
        -DENABLE_TESTDATA=0
        -DENABLE_TOOLS=0
    BUILD_COMMAND ${CMAKE_COMMAND} --build . --target aom ${AOM_BUILD_COMMAND}
    BUILD_ALWAYS 0
    INSTALL_COMMAND "")

add_custom_command(TARGET DepLibAom POST_BUILD
                  ${CUSTOM_POST_BUILD})
