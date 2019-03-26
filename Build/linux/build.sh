#!/bin/bash
#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

function clean {
    rm -R -f debug
    rm -R -f release
    rm -R -f ../../Bin/Debug
    rm -R -f ../../Bin/Release
}

function build {
    build_type=$1
    lowercase_build_type=$(echo ${build_type} | tr '[:upper:]' '[:lower:]')
    mkdir -p $lowercase_build_type
    mkdir -p ../../Bin/$build_type
    cd $lowercase_build_type
    PATH=$PATH:/usr/local/bin/
    cmake ../../..                                  \
        -DCMAKE_BUILD_TYPE=$build_type              \
        -DCMAKE_C_COMPILER=$CMAKE_COMPILER          \
        -DCMAKE_ASM_NASM_COMPILER=$CMAKE_ASSEMBLER  \

    # Compile the Library
    make -j $(if [ "$(uname -s)" = "Darwin" ]; then sysctl -n hw.ncpu; else nproc; fi) SvtAv1EncApp SvtAv1DecApp SvtAv1UnitTests
    cd ..
}

# Defines
CMAKE_ASSEMBLER=yasm
GCC_COMPILER=gcc
ICC_COMPILER=/opt/intel/composerxe/bin/icc

if [ ! -e $ICC_COMPILER ]; then
    CMAKE_COMPILER=$GCC_COMPILER
elif [ "$1" == "gcc" ]; then
    CMAKE_COMPILER=$GCC_COMPILER
elif [ "$2" == "gcc" ]; then
    CMAKE_COMPILER=$GCC_COMPILER
else
    CMAKE_COMPILER=$ICC_COMPILER
fi

cd $(dirname $(if [ "$(uname -s)" = "Darwin" ]; then perl -e 'use Cwd "abs_path";print abs_path(shift)' $0; else realpath $0; fi))

if [ $# -eq 0 ]; then
    build Debug
    build Release
elif [ "$1" = "clean" ]; then
    clean
elif [ "$1" = "debug" ]; then
    build Debug
elif [ "$1" = "release" ]; then
    build Release
elif [ "$1" = "cpp" ]; then
    build Debug
    build Release
elif [ "$1" = "all" ]; then
    build Debug
    build Release
elif [ "$1" = "gcc" ]; then
    build Debug
    build Release
else
    echo "build.sh <clean|all|debug|release|help>"
fi

exit
