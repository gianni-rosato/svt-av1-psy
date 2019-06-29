#!/bin/bash
#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

function build {
    build_type=$1
    lowercase_build_type=$(echo ${build_type} | tr '[:upper:]' '[:lower:]')
    mkdir -p $lowercase_build_type
    mkdir -p ../../Bin/$build_type
    cd $lowercase_build_type

    cmake_env=" -DCMAKE_BUILD_TYPE=${build_type}"
    cmake_env+=" -DCMAKE_C_COMPILER=${compiler}"
    cmake_env+=" -DCMAKE_CXX_COMPILER=${compilerxx}"
    cmake_env+=" -DCMAKE_INSTALL_PREFIX=${install_prefix}"
    cmake_env+=" -DCMAKE_ASM_NASM_COMPILER=${assembler}"
    if [ "$2" = "test" ]; then
        cmake_env+=" -DBUILD_SHARED_LIBS=ON"
        cmake_env+=" -DBUILD_TESTING=ON"
    else
        cmake_env+=" -DBUILD_SHARED_LIBS=${enable_shared}"
    fi

    if [ "x${target_system}" != "x" ] ; then
        cmake_env+=" -DCMAKE_SYSTEM_NAME=${target_system}"
        if [ "${target_system}" == "Windows" ] ; then
            cmake_env+=" -DCMAKE_RC_COMPILER=${windres}"
        fi
    fi
    PATH=$PATH:/usr/local/bin/
    cmake ${verbose:+"-v"} ../../.. ${cmake_env}

    # Compile the Library
    make VERBOSE=${verbose} -j ${make_threads} SvtAv1EncApp SvtAv1DecApp
    if [ "$2" = "test" ]; then
        make VERBOSE=${verbose} -j ${make_threads} SvtAv1UnitTests SvtAv1ApiTests SvtAv1E2ETests TestVectors
    fi
    cd ..
}

function clean {
    rm -R -f debug
    rm -R -f release
    rm -R -f ../../Bin/Debug
    rm -R -f ../../Bin/Release
}

function ncpu {
    if [ "$(uname -s)" = "Darwin" ]; then
        sysctl -n hw.ncpu
    else
        nproc
    fi
}

function absdir {
    if [ "$(uname -s)" = "Darwin" ]; then
        dirname $(perl -e 'use Cwd "abs_path";print abs_path(shift)' $1)
    else
        dirname $(realpath $1)
    fi
}

function probe_tools {
    compiler="${GCC_COMPILER}"
    if [ "x${CC}" != "x" ]; then
        compiler=${CC}
    elif [ -n "${cross_prefix}" ]; then
        if hash "${cross_prefix}${GCC_COMPILER}" 2> /dev/null; then
            compiler="${cross_prefix}${GCC_COMPILER}"
        elif hash "${cross_prefix}${CC_COMPILER}" 2> /dev/null; then
            compiler="${cross_prefix}${CC_COMPILER}"
        else
            echo "No suitable CC for ${cross_prefix}"
            usage
            exit 1
        fi
    elif [ -e "${ICC_COMPILER_PATH}/${ICC_COMPILER}" ]; then
        compiler="${ICC_COMPILER_PATH}/${ICC_COMPILER}"
    elif hash "${ICC_COMPILER}" 2> /dev/null; then
        compiler="${ICC_COMPILER}"
    elif hash "${GCC_COMPILER}" 2> /dev/null; then
        compiler="${GCC_COMPILER}"
    elif hash "${CLANG_COMPILER}" 2> /dev/null; then
        compiler="${CLANG_COMPILER}"
    fi

    # Choose CXX
    compilerxx="${GXX_COMPILER}"
    if [ "x${CXX}" != "x" ]; then
        compilerxx=${CXX}
    elif [ -n "${cross_prefix}" ]; then
        if hash "${cross_prefix}${GXX_COMPILER}" 2> /dev/null; then
            compilerxx="${cross_prefix}${GXX_COMPILER}"
        elif hash "${cross_prefix}${CXX_COMPILER}" 2> /dev/null; then
            compilerxx="${cross_prefix}${CXX_COMPILER}"
        else
            echo "No suitable CXX for ${cross_prefix}"
            usage
            exit 1
        fi
    elif [ -e "${ICC_COMPILER_PATH}/${ICPC_COMPILER}" ]; then
        compilerxx="${ICC_COMPILER_PATH}/${ICPC_COMPILER}"
    elif hash "${ICPC_COMPILER}" 2> /dev/null; then
        compilerxx="${ICPC_COMPILER}"
    elif hash "${GXX_COMPILER}" 2> /dev/null; then
        compilerxx="${GXX_COMPILER}"
    elif hash "${CLANGXX_COMPILER}" 2> /dev/null; then
        compilerxx="${CLANGXX_COMPILER}"
    fi

    if hash "${YASM_ASSEMBLER}" 2> /dev/null; then
        assembler="${YASM_ASSEMBLER}"
    else
            echo "No suitable assembler found"
            usage
            exit 1
    fi

    if [ "${target_system}" == "Windows" ] ; then
        windres="${cross_prefix}windres"
    fi
}

function run {
    # Run command
    if [ "$1" = "clean" ]; then
        clean
    elif [ "$1" = "debug" ]; then
        build Debug $2
    elif [ "$1" = "release" ]; then
        build Release $2
    elif [ "$1" = "all" ]; then
        build Debug $2
        build Release $2
    elif [ "$1" = "test" ]; then
        build Debug "test"
        build Release "test"
    else
        usage
        exit 1
    fi
}

function usage {
    cat <<EOF
Usage: $0 [-p install_prefix] [-c cross_prefix] [-s system] [-j threads] [-x] [-v] [-h] command
Options:
    -p install_prefix - Installation directory prefix (default ${install_prefix})
    -c cross_prefix   - cross-tools prefix (default "${cross_prefix}")
    -s system         - Target system (default "${target_system}")
    -j threads        - Parallel make threads
    -x                - Build static library. (default shared. Test target need shared library,
                        this parameter has no effect with test target on. )
    -v                - Verbose output
    -h                - Show this message
    command           - clean, all, all [test], release [test],
                        debug [test], test (default ${build_command})

Environment:
    CC                - Override C compiler
    CXX               - Override C++ compiler
EOF
}


####################
# main() starts here
####################

# CLI option defaults
cross_prefix=""
target_system=""
install_prefix="/usr/local"
enable_shared=ON
build_command=all
unit_test=
make_threads=$(ncpu)
verbose=

# Compilers, auto-detected
# May be overridden using CC and CXX environment variables
CC_COMPILER=cc
CLANG_COMPILER=clang
GCC_COMPILER=gcc
ICC_COMPILER=icc
ICC_COMPILER_PATH=/opt/intel/composerxe/bin

CXX_COMPILER=c++
CLANGXX_COMPILER=clang++
GXX_COMPILER=g++
ICPC_COMPILER=icpc

YASM_ASSEMBLER=yasm

windres=""

# Parse CLI options
while getopts "p:c:s:j:xvh" opt; do
    case ${opt} in
        p)
            install_prefix=${OPTARG}
            ;;
        c)
            cross_prefix=${OPTARG}
            ;;
        s)
            target_system=${OPTARG}
            ;;
        j)
            make_threads=${OPTARG}
            ;;
        x)
            enable_shared=OFF
            ;;
        v)
            set -x
            verbose=1
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

# Parse non-option arguments
shift $((OPTIND-1))
if [ $# -eq 1 ]; then
    build_command=$1
elif [ $# -eq 2 ]; then
    build_command=$1
    unit_test=$2
elif [ $# != 0 ]; then
    usage
    exit 1
fi

# Select suitable build tools
probe_tools

# Build directories subdirectories of build script's directory
cd $(absdir $0)

run ${build_command} ${unit_test}

exit
