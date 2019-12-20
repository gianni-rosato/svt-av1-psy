#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# C_DEFAULT Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(../../../API
        ../../Encoder/Codec
    ${PROJECT_SOURCE_DIR}/Source/Lib/Encoder/C_DEFAULT/)

file(GLOB all_files
    "*.h"
    "*.c")

add_library(ENCODER_C_DEFAULT OBJECT ${all_files})
