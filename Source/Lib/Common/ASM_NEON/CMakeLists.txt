#
# Copyright (c) 2023, Alliance for Open Media. All rights reserved
#
# This source code is subject to the terms of the BSD 2 Clause License and the
# Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License was
# not distributed with this source code in the LICENSE file, you can obtain it
# at www.aomedia.org/license/software. If the Alliance for Open Media Patent
# License 1.0 was not distributed with this source code in the PATENTS file, you
# can obtain it at www.aomedia.org/license/patent.
#

# Common/ASM_NEON Directory CMakeLists.txt

# Include Common Subdirectories
enable_language(ASM) # needed for .S files
add_library(COMMON_ASM_NEON OBJECT)
target_sources(
  COMMON_ASM_NEON
  PUBLIC aom_convolve8_neon.c
  PUBLIC av1_inv_txfm_neon.c
  PUBLIC av1_inv_txfm_neon.c
  PUBLIC cdef_block_neon.c
  PUBLIC cdef_filter_block_neon.c
  PUBLIC cfl_neon.c
  PUBLIC compound_convolve_neon.c
  PUBLIC convolve_2d_neon.c
  PUBLIC convolve_neon.c
  PUBLIC highbd_inv_txfm_neon.c
  PUBLIC EbBlend_a64_mask_neon.c
  PUBLIC EbDeblockingFilter_Intrinsic_neon.c
  PUBLIC EbIntraPrediction_neon.c
  PUBLIC EbPictureOperators_Intrinsic_neon.c
  PUBLIC highbd_convolve_2d_neon.c
  PUBLIC selfguided_neon.c
  PUBLIC sse_neon.c
  PUBLIC warp_plane_neon.c
  PUBLIC wiener_convolve_neon.c
  PUBLIC dav1d_asm.S
  PUBLIC dav1d_util.S
  PUBLIC itx.S)
target_include_directories(
  COMMON_ASM_NEON
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/API/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_NEON/)

# Causes the files to be explicitly compiled with `-x c`, which is wrong.
# https://cmake.org/cmake/help/git-master/policy/CMP0119.html
# set_source_files_properties(dav1d_asm.S dav1d_util.S itx.S PROPERTIES LANGUAGE
#                                                                       C)
