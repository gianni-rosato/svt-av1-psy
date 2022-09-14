#
# Copyright(c) 2019 Netflix, Inc.
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
#

# APP Directory CMakeLists.txt

# Include Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/test/
    ${PROJECT_SOURCE_DIR}/test/ref
    ${PROJECT_SOURCE_DIR}/third_party/googletest/include third_party/googletest/src
    ${PROJECT_SOURCE_DIR}/Source/API
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Globals
    ${PROJECT_SOURCE_DIR}/Source/Lib/Decoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/App/EncApp
    ${PROJECT_SOURCE_DIR}/Source/API
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/)

# Define helper functions and macros used by Google Test.
include(../third_party/googletest/cmake/internal_utils.cmake)

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_ALLOW_KEYWORD_MACROS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALLOW_KEYWORD_MACROS")
endif()

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(gtest_force_shared_crt
    "Use shared (DLL) run-time lib even when Google Test is built as static lib."
    OFF)

add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1)

enable_testing()

if(UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()

set(all_files
    AdaptiveScanTest.cc
    BitstreamWriterTest.cc
    CdefTest.cc
    CompoundUtilTest.cc
    DeblockTest.cc
    EbBlend_a64_mask_1d_test.cc
    EbBlend_a64_mask_test.cc
    EbHighbdIntraPredictionTests.cc
    EbHighbdIntraPredictionTests.h
    EbNoiseModel_test.cc
    EbUnitTest.h
    EbUnitTestUtility.c
    EbUnitTestUtility.h
    EncodeTxbAsmTest.cc
    FFTTest.cc
    FilmGrainExpectedResult.h
    FilmGrainTest.cc
    FilterIntraPredTest.cc
    ForwardtransformTests.cc
    FwdTxfm1dTest.cc
    FwdTxfm2dAsmTest.cc
    FwdTxfm2dTest.cc
    GlobalMotionUtilTest.cc
    HbdVarianceTest.cc
    IntraBcUtilTest.cc
    InvTxfm1dTest.cc
    InvTxfm2dAsmTest.cc
    MotionEstimationTest.cc
    OBMCSadTest.cc
    OBMCVarianceTest.cc
    PackUnPackTest.cc
    PaletteModeUtilTest.cc
    PictureOperatorTest.cc
    PsnrTest.cc
    QuantAsmTest.cc
    ResidualTest.cc
    ResizeTest.cc
    RestorationPickTest.cc
    SadTest.cc
    SelfGuidedUtilTest.cc
    SpatialFullDistortionTest.cc
    TemporalFilterTestPlanewise.cc
    TestEnv.c
    TxfmCommon.h
    VarianceTest.cc
    WedgeUtilTest.cc
    acm_random.h
    av1_convolve_scale_test.cc
    compute_mean_test.cc
    convolve_2d_test.cc
    corner_match_test.cc
    dwt_test.cc
    frame_error_test.cc
    hadamard_test.cc
    intrapred_cfl_test.cc
    intrapred_dr_test.cc
    intrapred_edge_filter_test.cc
    intrapred_test.cc
    quantize_func_test.cc
    random.h
    selfguided_filter_test.cc
    subtract_avg_cfl_test.cc
    util.h
    warp_filter_test.cc
    warp_filter_test_util.cc
    warp_filter_test_util.h
    wiener_convolve_test.cc
    ref/TxfmRef.cc
    ref/TxfmRef.h

    # Decoder sources
    ../Source/Lib/Decoder/Codec/EbDecBitstreamUnit.c
    )

set(lib_list
    $<TARGET_OBJECTS:COMMON_CODEC>
    $<TARGET_OBJECTS:FASTFEAT>
    $<TARGET_OBJECTS:COMMON_C_DEFAULT>
    $<TARGET_OBJECTS:COMMON_ASM_SSE2>
    $<TARGET_OBJECTS:COMMON_ASM_SSSE3>
    $<TARGET_OBJECTS:COMMON_ASM_SSE4_1>
    $<TARGET_OBJECTS:COMMON_ASM_AVX2>
    $<TARGET_OBJECTS:COMMON_ASM_AVX512>
    $<TARGET_OBJECTS:ENCODER_C_DEFAULT>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE2>
    $<TARGET_OBJECTS:ENCODER_ASM_SSSE3>
    $<TARGET_OBJECTS:ENCODER_ASM_SSE4_1>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX2>
    $<TARGET_OBJECTS:ENCODER_ASM_AVX512>
    $<TARGET_OBJECTS:ENCODER_GLOBALS>
    $<TARGET_OBJECTS:ENCODER_CODEC>
    cpuinfo_public
    gtest_all)
if(UNIX)
  # App Source Files
    add_executable(SvtAv1UnitTests
      ${all_files})
    # Link the Encoder App
    target_link_libraries(SvtAv1UnitTests
        ${lib_list}
        pthread
        m)
else()
    cxx_executable_with_flags(SvtAv1UnitTests
        "${cxx_default}"
        "${lib_list}"
        ${all_files})
endif()

add_dependencies(SvtAv1UnitTests SvtAv1Enc SvtAv1Dec)

install(TARGETS SvtAv1UnitTests RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(SvtAv1UnitTests ${CMAKE_OUTPUT_DIRECTORY}/SvtAv1UnitTests)

add_subdirectory(api_test)
add_subdirectory(e2e_test)
