#!/bin/sh
#
# Copyright(c) 2019 Intel Corporation
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at www.aomedia.org/license/software. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at www.aomedia.org/license/patent.
#

# Sets context of if running in the script.
# Helpful when copying and pasting functions and debuging.
if printf '%s' "$0" | grep -q '\.sh'; then
    IN_SCRIPT=true
fi

# Fails the script if any of the commands error (Other than if and some others)
set -e

# Function to print to stderr or stdout while stripping preceding tabs(spaces)
# If tabs are desired, use cat
print_message() {
    if [ "$1" = "stdout" ] || [ "$1" = "-" ]; then
        shift
        printf '%b\n' "${@:-Unknown Message}" | sed 's/^[ \t]*//'
    else
        printf '%b\n' "${@:-Unknown Message}" | sed 's/^[ \t]*//' >&2
    fi
}

# Convenient function to exit with message
die() {
    print_message "${@:-Unknown error}"
    if ${IN_SCRIPT:-false}; then
        print_message "The script will now exit."
        exit 1
    fi
}

# Function for changing directory.
cd_safe() {
    if (cd "$1"); then
        cd "$1"
    else
        _dir="$1"
        shift
        die "${@:-Failed cd to $_dir.}"
    fi
}

# used for resolving certain paths relative to original caller like -t
prev_pwd=
${IN_SCRIPT:-false} && {
    prev_pwd=$PWD
    cd_safe "$(cd "$(dirname "$0")" > /dev/null 2>&1 && pwd -P)"
}

# Help message
echo_help() {
    cat << EOF
Usage: $0 [OPTION] ... -- [OPTIONS FOR CMAKE]

-a, --all, all          Builds release and debug
    --asm, asm=*        Set assembly compiler [$ASM]
-b, --bindir, bindir=*  Directory to install binaries
    --cc, cc=*          Set C compiler [$CC]
    --cxx, cxx=*        Set CXX compiler [$CXX]
    --c-only, c-only    Compile only C code
    --clean, clean      Remove build and Bin folders
    --debug, debug      Build debug
    --shared, shared    Build shared libs
-x, --static, static    Build static libs
-g, --gen, gen=*        Set CMake generator
-i, --install, install  Install build [Default release]
-j, --jobs, jobs=*      Set number of jobs for make/CMake [$jobs]
    --no-enc, no-enc    Don't build the encoder app and libs
    --no-dec, no-dec    Don't build the decoder app and libs
-p, --prefix, prefix=*  Set installation prefix
    --release, release  Build release
-s, --target_system,    Set CMake target system
    target_system=*
    --test, test        Build Unit Tests
-t, --toolchain,        Set CMake toolchain file
    toolchain=*
-v, --verbose, verbose  Print out commands

Example usage:
    build.sh -xi debug test
    build.sh jobs=8 all cc=clang cxx=clang++
    build.sh -j 4 all -t "https://gist.githubusercontent.com/peterspackman/8cf73f7f12ba270aa8192d6911972fe8/raw/mingw-w64-x86_64.cmake"
    build.sh generator=Xcode cc=clang
EOF
}

# Usage: build <release|debug> [test]
build() (
    build_type=Release
    while [ -n "$*" ]; do
        case $(printf %s "$1" | tr '[:upper:]' '[:lower:]') in
        release) build_type=Release && shift ;;
        debug) build_type=Debug && shift ;;
        *) break ;;
        esac
    done

    rm -rf $build_type
    mkdir -p $build_type > /dev/null 2>&1
    cd_safe $build_type

    cmake ../../.. -DCMAKE_BUILD_TYPE=$build_type $CMAKE_EXTRA_FLAGS "$@"

    if [ -f Makefile ]; then
        make -j "$jobs"
        return
    fi

    set --
    if cmake --build 2>&1 | grep -q parallel; then
        set -- --parallel "$jobs"
    fi

    # Compile the Library
    cmake --build . --config $build_type "$@"
)

check_executable() (
    print_exec=false
    while true; do
        case $1 in
        -p) print_exec=true && shift ;;
        *) break ;;
        esac
    done
    [ -n "$1" ] && command_to_check="$1" || return 1
    shift
    if [ -x "$command_to_check" ]; then
        $print_exec && printf '%s\n' "$command_to_check"
        return 0
    fi
    for d in "$@" $(printf '%s ' "$PATH" | tr ':' ' '); do
        if [ -x "$d/$command_to_check" ]; then
            $print_exec && printf '%s\n' "$d/$command_to_check"
            return 0
        fi
    done
    return 127
)

install_build() (
    build_type=Release
    sudo=$(check_executable -p sudo)
    while [ -n "$*" ]; do
        case $(printf %s "$1" | tr '[:upper:]' '[:lower:]') in
        release) build_type="Release" && shift ;;
        debug) build_type="Debug" && shift ;;
        *) break ;;
        esac
    done

    { [ -d $build_type ] && cd_safe $build_type; } ||
        die "Unable to find the build folder. Did the build command run?"
    cmake --build . --target install --config $build_type || {
        test -n "$sudo" &&
            eval ${sudo:+echo cmake failed to install, trying with sudo && $sudo cmake --build . --target install --config $build_type}
    } || die "Unable to run install"
)

if [ -z "$CC" ] && [ "$(uname -a | cut -c1-5)" != "MINGW" ]; then
    ! CC=$(check_executable -p icc /opt/intel/bin) &&
        ! CC=$(check_executable -p gcc) &&
        ! CC=$(check_executable -p clang) &&
        ! CC=$(check_executable -p cc) &&
        die "No suitable c compiler found in path" \
            "Please either install one or set it via cc=*"
    export CC
fi

if [ -z "$CXX" ] && [ "$(uname -a | cut -c1-5)" != "MINGW" ]; then
    ! CXX=$(check_executable -p icpc "/opt/intel/bin") &&
        ! CXX=$(check_executable -p g++) &&
        ! CXX=$(check_executable -p clang++) &&
        ! CXX=$(check_executable -p c++) &&
        die "No suitable cpp compiler found in path" \
            "Please either install one or set it via cxx=*"
    export CXX
fi

build_release=false
build_debug=false
build_install=false

parse_options() {
    while true; do
        [ -z "$1" ] && break
        case $(printf %s "$1" | tr '[:upper:]' '[:lower:]') in
        help) echo_help && ${IN_SCRIPT:-false} && exit ;;
        all) build_debug=true build_release=true && shift ;;
        asm=*)
            check_executable "${1#*=}" &&
                CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_ASM_NASM_COMPILER=$(check_executable -p "${1#*=}")" &&
                case $1 in
                *nasm*) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DENABLE_NASM=ON" ;;
                esac
            shift
            ;;
        bindir=*)
            CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_INSTALL_BINDIR=${1#*=}"
            shift
            ;;
        cc=*)
            if check_executable "${1#*=}"; then
                CC=$(check_executable -p "${1#*=}")
                export CC
            fi
            shift
            ;;
        cxx=*)
            if check_executable "${1#*=}"; then
                CXX=$(check_executable -p "${1#*=}")
                export CXX
            fi
            shift
            ;;
        c-only) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCOMPILE_C_ONLY=ON" && shift ;;
        clean)
            for d in *; do
                [ -d "$d" ] && rm -rf "$d"
            done
            for d in ../../Bin/*; do
                [ -d "$d" ] && rm -rf "$d"
            done
            shift && ${IN_SCRIPT:-false} && exit
            ;;
        debug) build_debug=true && shift ;;
        shared) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_SHARED_LIBS=ON" && shift ;;
        static) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_SHARED_LIBS=OFF" && shift ;;
        gen=*) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -G${1#*=}" && shift ;;
        install) build_install=true && shift ;;
        jobs=*) jobs="${1#*=}" && shift ;;
        no-enc) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_ENC=OFF" && shift ;;
        no-dec) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_DEC=OFF" && shift ;;
        prefix=*) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_INSTALL_PREFIX=${1#*=}" && shift ;;
        release) build_release=true && shift ;;
        target_system=*)
            CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_SYSTEM_NAME=${1#*=}"
            shift
            ;;
        tests) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_TESTING=ON" && shift ;;
        toolchain=*)
            toolchain=''
            url=${1#*=}
            case $(echo "$url" | cut -c1-4) in
            http*)
                toolchain=${url%%\?*}
                toolchain=${toolchain##*/}
                toolchain=$PWD/$toolchain
                curl --connect-timeout 15 --retry 3 --retry-delay 5 -sfLk -o "$toolchain" "$url"
                ;;
            *) toolchain=$(
                case ${url%/*} in
                */*)
                    [ -n "$prev_pwd" ] && cd "$prev_pwd"
                    cd "${url%/*}"
                    ;;
                esac
                pwd -P 2> /dev/null || pwd
            )/${url##*/} ;;
            esac
            CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_TOOLCHAIN_FILE=$toolchain" && shift
            ;;
        verbose) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_VERBOSE_MAKEFILE=1" && shift ;;
        *) print_message "Unknown option: $1" && shift ;;
        esac
    done
}

parse_equal_option() {
    case $1 in
    *=*) parse_options "$(printf %s "$1" | cut -c3- | cut -d= -f1 | tr '[:upper:]' '[:lower:]')=${1#*=}" ;;
    *) parse_options "$(printf %s "$1" | cut -c3- | cut -d= -f1 | tr '[:upper:]' '[:lower:]')=$2" ;;
    esac
}

if [ -z "$*" ]; then
    build_release=true
else
    while [ -n "$*" ]; do
        case $1 in
        --*) # Handle --* based args
            case $(printf %s "$1" | cut -c3- | cut -d= -f1 | tr '[:upper:]' '[:lower:]') in
            # Stop on "--", pass the rest to cmake
            "") shift && break ;;
            help) parse_options help && shift ;;
            all) parse_options debug release && shift ;;
            c-only) parse_options c-only && shift ;;
            clean) parse_options clean && shift ;;
            debug) parse_options debug && shift ;;
            install) parse_options install && shift ;;
            no-enc) parse_options no-enc && shift ;;
            no-dec) parse_options no-dec && shift ;;
            release) parse_options release && shift ;;
            shared) parse_options shared && shift ;;
            static) parse_options static && shift ;;
            toolchain) parse_options toolchain="$2" && shift ;;
            test) parse_options tests && shift ;;
            verbose) parse_options verbose && shift ;;
            asm | bindir | cc | cxx | gen | jobs | prefix | target_system)
                parse_equal_option "$1" "$2"
                case $1 in
                *=*) shift ;;
                *) shift 2 ;;
                esac
                ;;
            *) die "Error, unknown option: $1" ;;
            esac
            ;;
        -*) # Handle -* based args. Currently doesn't differentiate upper and lower since there's not need at the momment.
            i=2
            match=$1
            shift
            while [ $i -ne $((${#match} + 1)) ]; do
                case $(echo "$match" | cut -c$i | tr '[:upper:]' '[:lower:]') in
                h) parse_options help ;;
                a)
                    parse_options all
                    i=$((i + 1))
                    ;;
                b)
                    parse_options bindir="$1"
                    i=$((i + 1))
                    shift
                    ;;
                g)
                    case $(echo "$match" | cut -c$((i + 1))-) in
                    "")
                        # if it's -g Ninja
                        parse_options gen="$1"
                        i=$((i + 1))
                        shift
                        ;;
                    *)
                        # if it's -GNinja
                        # Just put everything past -g as the generator
                        parse_options gen="$(echo "$match" | cut -c$((i + 1))-)"
                        # go ahead and skip this block
                        i=$((${#match} + 1))
                        ;;
                    esac
                    ;;
                i)
                    parse_options install
                    i=$((i + 1))
                    ;;
                j)
                    case $(echo "$match" | cut -c$((i + 1))-) in
                    *[!0-9]*)
                        parse_options jobs="$1"
                        i=$((i + 1))
                        shift
                        ;;
                    *)
                        # Found number right after
                        parse_options jobs="$(echo "$match" | cut -c$((i + 1))-)"
                        # go ahead and skip this block
                        i=$((${#match} + 1))
                        ;;
                    esac
                    ;;
                p)
                    parse_options prefix="$1"
                    i=$((i + 1))
                    shift
                    ;;
                s)
                    parse_options target_system="$1"
                    i=$((i + 1))
                    shift
                    ;;
                t)
                    parse_options toolchain="$1"
                    i=$((i + 1))
                    shift
                    ;;
                x)
                    parse_options static
                    i=$((i + 1))
                    ;;
                v)
                    parse_options verbose
                    i=$((i + 1))
                    ;;
                *) die "Error, unknown option: -$(echo "$match" | cut -c$i | tr '[:upper:]' '[:lower:]')" ;;
                esac
            done
            ;;
        *) # Handle single word args
            case $(printf %s "$1" | tr '[:upper:]' '[:lower:]') in
            all) parse_options release debug && shift ;;
            asm=*) parse_options asm="${1#*=}" && shift ;;
            bindir=*) parse_options bindir="${1#*=}" && shift ;;
            cc=*) parse_options cc="${1#*=}" && shift ;;
            cxx=*) parse_options cxx="${1#*=}" && shift ;;
            c-only) parse_options c-only && shift ;;
            clean) parse_options clean && shift ;;
            debug) parse_options debug && shift ;;
            gen=*) parse_options gen="${1#*=}" && shift ;;
            help) parse_options help && shift ;;
            install) parse_options install && shift ;;
            jobs=*) parse_options jobs="${1#*=}" && shift ;;
            prefix=*) parse_options prefix="${1#*=}" && shift ;;
            no-enc) parse_options no-enc && shift ;;
            no-dec) parse_options no-dec && shift ;;
            target_system=*) parse_options target_system="${1#*=}" && shift ;;
            shared) parse_options shared && shift ;;
            static) parse_options static && shift ;;
            release) parse_options release && shift ;;
            test) parse_options tests && shift ;;
            toolchain=*) parse_options toolchain="${1#*=}" && shift ;;
            verbose) parse_options verbose && shift ;;
            end) ${IN_SCRIPT:-false} && exit ;;
            *) die "Error, unknown option: $1" ;;
            esac
            ;;
        esac
    done
fi

case $jobs in
*[!0-9]* | "") jobs=$(getconf _NPROCESSORS_ONLN 2> /dev/null || nproc 2> /dev/null || sysctl -n hw.ncpu 2> /dev/null) ;;
esac

[ "${PATH#*\/usr\/local\/bin}" = "$PATH" ] && PATH=$PATH:/usr/local/bin

if $build_debug && $build_release; then
    build release "$@"
    build debug "$@"
elif $build_debug; then
    build debug "$@"
else
    build_release=true
    build release "$@"
fi

if $build_install; then
    if $build_release; then
        install_build release
    else
        install_build debug
    fi
fi
