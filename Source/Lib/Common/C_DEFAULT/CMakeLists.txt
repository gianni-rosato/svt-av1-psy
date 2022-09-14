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
# C_DEFAULT Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)

set(all_files
    EbBlend_a64_mask_c.c
    EbInterPrediction_c.c
    EbIntraPrediction_c.c
    EbPackUnPack_C.c
    EbPackUnPack_C.h
    EbPictureOperators_C.c
    EbPictureOperators_C.h
    cfl_c.c
    filterintra_c.c
    )

add_library(COMMON_C_DEFAULT OBJECT ${all_files})
