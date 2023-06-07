#
# Copyright(c) 2019 Intel Corporation
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
#
# Common/Codec Directory CMakeLists.txt

if(EXCLUDE_HASH)
    set(PACKAGE_VERSION_STRING ${PROJECT_VERSION})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/EbVersion.h.in EbVersion.h @ONLY)
    add_custom_target(EbVersionHeaderGen)
else()
    add_custom_target(EbVersionHeaderGen
        ${CMAKE_COMMAND}
            -D INPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/EbVersion.h.in
            -D OUTPUT_FILE=EbVersion.h
            -D GIT_ROOT_DIR=${PROJECT_SOURCE_DIR}
            -D PACKAGE_VERSION_STRING=${PROJECT_VERSION}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/ConfigureGitVersion.cmake
        BYPRODUCTS EbVersion.h
        COMMENT "Generating version header"
        VERBATIM)
endif()

if(NOT COMPILE_C_ONLY AND HAVE_X86_PLATFORM)
# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX512/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/
    ${PROJECT_SOURCE_DIR}/third_party/cpuinfo/include/)
else ()
# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/third_party/fastfeat/
    ${PROJECT_SOURCE_DIR}/third_party/cpuinfo/include/)
endif ()

set(all_files
    Av1Common.h
    EbAv1Structs.h
    EbAvcStyleMcp.h
    EbBitstreamUnit.c
    EbBitstreamUnit.h
    EbBlend_a64_mask.c
    EbBlockStructures.c
    EbBlockStructures.h
    EbCabacContextModel.c
    EbCabacContextModel.h
    EbCdef.c
    EbCdef.h
    EbCoefficients.h
    EbCommonUtils.h
    EbDeblockingCommon.c
    EbDeblockingCommon.h
    EbDefinitions.h
    EbInterPrediction.c
    EbInterPrediction.h
    EbIntraCommon.h
    EbIntraPrediction.c
    EbIntraPrediction.h
    EbInvTransforms.c
    EbInvTransforms.h
    EbLog.c
    EbLog.h
    EbMalloc.c
    EbMalloc.h
    EbMcp.c
    EbMcp.h
    EbObject.h
    EbPackUnPack.h
    EbPictureBufferDesc.c
    EbPictureBufferDesc.h
    EbPictureOperators.c
    EbPictureOperators.h
    EbQMatrices.h
    EbRestoration.c
    EbRestoration.h
    EbSegmentationParams.c
    EbSegmentationParams.h
    EbSuperRes.c
    EbSuperRes.h
    EbSystemResourceManager.c
    EbSystemResourceManager.h
    EbThreads.c
    EbThreads.h
    EbTime.c
    EbTime.h
    EbUtility.c
    EbUtility.h
    EbWarpedMotion.c
    EbWarpedMotion.h
    common_dsp_rtcd.c
    common_dsp_rtcd.h
    convolve.c
    convolve.h
    filter.h
    grainSynthesis.c
    grainSynthesis.h
    )

add_library(COMMON_CODEC OBJECT ${all_files})
