include(ExternalProject)
#/property:Configuration=${CMAKE_BUILD_TYPE} ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/
if (MSVC OR MSYS OR MINGW OR WIN32)
set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/msvc)
set(CUSTOM_CONFIG -G "Visual Studio 15 2017 Win64")
set(CUSTOM_BUILD_CMD msbuild /p:Configuration=MinSizeRel AOM.sln
                     COMMAND msbuild /p:Configuration=Debug AOM.sln)

set(CUSTOM_POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/Release/
                      COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/Release/
                      COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/MinSizeRel/
                      COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/MinSizeRel/
                      COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/RelWithDebInfo/
                      COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/MinSizeRel/aom.lib ${TARGET_OUTPUT_PATH}/RelWithDebInfo/
                      COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}/Debug/
                      COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/Debug/aom.lib ${TARGET_OUTPUT_PATH}/Debug/)
endif(MSVC OR MSYS OR MINGW OR WIN32)


if (UNIX AND NOT APPLE)
set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/linux)
set(TARGET "libaom.so*")
install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.so DESTINATION lib)
install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.so.0 DESTINATION lib)
endif(UNIX AND NOT APPLE)

if (APPLE)
set(TARGET_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/aom/lib/mac)
set(TARGET "libaom*dylib")
install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.dylib DESTINATION lib)
install(FILES ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/libaom.0.dylib DESTINATION lib)
endif(APPLE)

if (UNIX)
set(CUSTOM_CONFIG "-DBUILD_SHARED_LIBS=1")
set(CUSTOM_BUILD_CMD make aom)
set(CUSTOM_POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET_OUTPUT_PATH}
                      COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_BINARY_DIR}/libaom/src/DepLibAom-build/${TARGET} ${TARGET_OUTPUT_PATH})


endif(UNIX)


ExternalProject_Add(DepLibAom
  PREFIX "${CMAKE_BINARY_DIR}/libaom"
  GIT_REPOSITORY "https://github.com/Cidana-Developers/aom.git"
  GIT_TAG av1-normative
  CMAKE_ARGS
    ${CUSTOM_CONFIG}
    -DCONFIG_INSPECTION=1
    -DENABLE_TESTS=0
    -DCONFIG_AV1_ENCODER=0
    -DENABLE_DOCS=0
    -DENABLE_EXAMPLES=0
    -DENABLE_TESTDATA=0
    -DENABLE_TOOLS=0
  BUILD_COMMAND ${CUSTOM_BUILD_CMD}
  BUILD_ALWAYS 0
  INSTALL_COMMAND ""
  )

add_custom_command(TARGET DepLibAom POST_BUILD
                  ${CUSTOM_POST_BUILD})
