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

# Common/ASM_NEON Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
    ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_NEON/)

set(neon_files
    aom_convolve8_neon.c
    compound_convolve_neon.c
    compound_convolve_neon.h
    convolve_neon.c
    convolve_neon.h
    mem_neon.h
    sse_neon.c
    transpose_neon.h
    warp_plane_neon.c
    warp_plane_neon.h
    wiener_convolve_neon.c
    convolve_2d_neon.c
    )

add_library(COMMON_ASM_NEON OBJECT ${neon_files})
