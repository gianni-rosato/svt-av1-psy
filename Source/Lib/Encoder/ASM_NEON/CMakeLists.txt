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

# Encoder/ASM_NEON Directory CMakeLists.txt

# Include Encoder Subdirectories
add_library(ENCODER_ASM_NEON OBJECT)
target_sources(
  ENCODER_ASM_NEON
  PUBLIC pickrst_neon.c
  PUBLIC sad_neon.c
  PUBLIC EbComputeSAD_neon.c
  PUBLIC variance_neon.c
  PUBLIC obmc_sad_neon.c
  PUBLIC obmc_variance_neon.c
  PUBLIC highbd_fwd_txfm_neon.c
  PUBLIC av1_quantize_neon.c)
target_include_directories(
  ENCODER_ASM_NEON
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/API/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/Codec/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/
  PRIVATE ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/ASM_NEON/)
