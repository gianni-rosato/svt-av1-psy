#
# Copyright (c) 2023, Alliance for Open Media. All rights reserved
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at www.aomedia.org/license/software. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at www.aomedia.org/license/patent.
#

# Encoder/ASM_NEON Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(../../../API
        ../../Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_NEON/)

set(neon_files
    pickrst_neon.c
    sad_neon.c
    sum_neon.h
    EbComputeSAD_neon.c
    variance_neon.c
    obmc_sad_neon.c
    obmc_variance_neon.c
    )

add_library(ENCODER_ASM_NEON OBJECT ${neon_files})
