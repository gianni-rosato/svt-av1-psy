#!/bin/sh
#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
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

${IN_SCRIPT:-false} && cd_safe "$(cd "$(dirname "$0")" > /dev/null 2>&1 && pwd -P)"

# Help message
echo_help() {
    cat << EOF
Usage: $0 [OPTION] ... -- [OPTIONS FOR CMAKE]

-a, --all, all          Builds release and debug
    --asm, asm=*        Set assembly compiler [$ASM]
-b, --bindir, bindir=*  Directory to install binaries
    --cc, cc=*          Set C compiler [$CC]
    --cxx, cxx=*        Set CXX compiler [$CXX]
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
    build.sh jobs=8 all cc=clang cpp=clang++
    build.sh -j 4 all -t "https://gist.githubusercontent.com/peterspackman/8cf73f7f12ba270aa8192d6911972fe8/raw/mingw-w64-x86_64.cmake"
    build.sh generator=Xcode cc=clang
EOF
}

# Usage: build <release|debug> [test]
build() (
    build_type=Release
    while [ -n "$*" ]; do
        case "$(printf %s "$1" | tr '[:upper:]' '[:lower:]')" in
        release) build_type="Release" && shift ;;
        debug) build_type="Debug" && shift ;;
        *) break ;;
        esac
    done

    mkdir -p "$build_type" > /dev/null 2>&1
    cd_safe "$build_type"

    for file in *; do
        rm -rf "$file"
    done

    cmake ../../.. -DCMAKE_BUILD_TYPE="$build_type" $CMAKE_EXTRA_FLAGS "$@"

    # Compile the Library
    if [ -f Makefile ]; then
        make -j "${jobs:-4}"
    elif cmake --build 2>&1 | grep -q -- '-j'; then
        cmake --build . --config "$build_type" -j "${jobs:-4}"
    else
        cmake --build . --config "$build_type"
    fi
)

check_executable() (
    print_exec=false
    while true; do
        case "$1" in
        -p) print_exec=true && shift ;;
        *) break ;;
        esac
    done
    [ -n "$1" ] && command_to_check="$1" || return 1
    shift
    if [ -e "$command_to_check" ]; then
        $print_exec && printf '%s\n' "$command_to_check"
        return 0
    fi
    for d in "$@" $(printf '%s ' "$PATH" | tr ':' ' '); do
        if [ -e "$d/$command_to_check" ]; then
            $print_exec && printf '%s\n' "$d/$command_to_check"
            return 0
        fi
    done
    return 127
)

install_build() (
    build_type="Release"
    sudo=
    while [ -n "$*" ]; do
        case $(printf %s "$1" | tr '[:upper:]' '[:lower:]') in
        release) build_type="Release" && shift ;;
        debug) build_type="Debug" && shift ;;
        *) break ;;
        esac
    done
    check_executable sudo && sudo -v > /dev/null 2>&1 && sudo=sudo
    { [ -d "$build_type" ] && cd_safe "$build_type"; } ||
        die "Unable to find the build folder. Did the build command run?"
    $sudo cmake --build . --target install --config "$build_type" ||
        die "Unable to run install"
)

if [ -z "$CC" ] && [ "$(uname -a | cut -c1-5)" != "MINGW" ]; then
    if check_executable icc /opt/intel/bin; then
        CC=$(check_executable -p icc /opt/intel/bin)
    elif check_executable gcc; then
        CC=$(check_executable -p gcc)
    elif check_executable clang; then
        CC=$(check_executable -p clang)
    elif check_executable cc; then
        CC=$(check_executable -p cc)
    else
        die "No suitable c compiler found in path" \
            "Please either install one or set it via cc=*"
    fi
    export CC
fi

if [ -z "$CXX" ] && [ "$(uname -a | cut -c1-5)" != "MINGW" ]; then
    if check_executable icpc "/opt/intel/bin"; then
        CXX=$(check_executable -p icpc "/opt/intel/bin")
    elif check_executable g++; then
        CXX=$(check_executable -p g++)
    elif check_executable clang++; then
        CXX=$(check_executable -p clang++)
    elif check_executable c++; then
        CC=$(check_executable -p c++)
    else
        die "No suitable cpp compiler found in path" \
            "Please either install one or set it via cxx=*"
    fi
    export CXX
fi

if [ -z "$jobs" ]; then
    jobs=$(getconf _NPROCESSORS_ONLN) ||
        jobs=$(nproc) ||
        jobs=$(sysctl -n hw.ncpu) ||
        jobs=2
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
                CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_ASM_NASM_COMPILER=$(check_executable -p "${1#*=}")"
            shift
            ;;
        bindir=*)
            CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -CMAKE_INSTALL_BINDIR=${1#*=}"
            shift
            ;;
        cc=*)
            if check_executable "${1#*=}"; then
                CC="$(check_executable -p "${1#*=}")"
                export CC
            fi
            shift
            ;;
        cxx=*)
            if check_executable "${1#*=}"; then
                CXX="$(check_executable -p "${1#*=}")"
                export CXX
            fi
            shift
            ;;
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
        build_shared) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_SHARED_LIBS=ON" && shift ;;
        build_static) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DBUILD_SHARED_LIBS=OFF" && shift ;;
        generator=*) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -G${1#*=}" && shift ;;
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
            toolchain=""
            url="${1#*=}"
            if [ "$(echo "$url" | cut -c1-4)" = "http" ]; then
                toolchain="${url%%\?*}"
                toolchain="${toolchain##*/}"
                toolchain="$PWD/$toolchain"
                curl --connect-timeout 15 --retry 3 --retry-delay 5 -sfLk -o "$toolchain" "$url"
            else
                toolchain="$url"
            fi
            CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_TOOLCHAIN_FILE=../$toolchain" && shift
            ;;
        verbose) CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS -DCMAKE_VERBOSE_MAKEFILE=1" && shift ;;
        *) print_message "Unknown option: $1" && shift ;;
        esac
    done
}

if [ -z "$*" ]; then
    build_release=true
else
    while [ -n "$*" ]; do
        # Handle --* based args
        if [ "$(printf %s "$1" | cut -c1-2)" = "--" ]; then
            # Stop on "--", pass the rest to cmake
            [ -z "$(printf %s "$1" | cut -c3-)" ] && shift && break
            case "$(printf %s "$1" | cut -c3- | tr '[:upper:]' '[:lower:]')" in
            help) parse_options help && shift ;;
            all) parse_options debug release && shift ;;
            asm) parse_options asm="$2" && shift 2 ;;
            bindir) parse_options bindir="$2" && shift 2 ;;
            cc) parse_options cc="$2" && shift 2 ;;
            cxx) parse_options cxx="$2" && shift 2 ;;
            clean) parse_options clean && shift ;;
            debug) parse_options debug && shift ;;
            gen) parse_options generator="$2" && shift 2 ;;
            install) parse_options install && shift ;;
            jobs) parse_options jobs="$2" && shift 2 ;;
            no-enc) parse_options no-enc && shift ;;
            no-dec) parse_options no-dec && shift ;;
            prefix) parse_options prefix="$2" && shift 2 ;;
            release) parse_options release && shift ;;
            shared) parse_options build_shared && shift ;;
            static) parse_options build_static && shift ;;
            target_system) parse_options target_system="$2" && shift 2 ;;
            toolchain) parse_options toolchain="$2" && shift ;;
            test) parse_options tests && shift ;;
            verbose) parse_options verbose && shift ;;
            *) die "Error, unknown option: $1" ;;
            esac
        # Handle -* based args. Currently doesn't differentiate upper and lower since there's not need at the momment.
        elif [ "$(printf %s "$1" | cut -c1)" = "-" ]; then
            i=2
            match="$1"
            shift
            while [ $i -ne $((${#match} + 1)) ]; do
                case "$(echo "$match" | cut -c$i | tr '[:upper:]' '[:lower:]')" in
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
                    parse_options generator="$1"
                    i=$((i + 1))
                    shift
                    ;;
                i)
                    parse_options install
                    i=$((i + 1))
                    ;;
                j)
                    parse_options jobs="$1"
                    i=$((i + 1))
                    shift
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
                    parse_options build_static
                    i=$((i + 1))
                    ;;
                v)
                    parse_options verbose
                    i=$((i + 1))
                    ;;
                *) die "Error, unknown option: -$(echo "$match" | cut -c$i | tr '[:upper:]' '[:lower:]')" ;;
                esac
            done
        # Handle single word args
        else
            case "$(printf %s "$1" | tr '[:upper:]' '[:lower:]')" in
            all) parse_options release debug && shift ;;
            asm=*) parse_options asm="${1#*=}" && shift ;;
            bindir=*) parse_options bindir="${1#*=}" && shift ;;
            cc=*) parse_options cc="${1#*=}" && shift ;;
            cxx=*) parse_options cxx="${1#*=}" && shift ;;
            clean) parse_options clean && shift ;;
            debug) parse_options debug && shift ;;
            gen=*) parse_options generator="${1#*=}" && shift ;;
            help) parse_options help && shift ;;
            install) parse_options install && shift ;;
            jobs=*) parse_options jobs="${1#*=}" && shift ;;
            prefix=*) parse_options prefix="${1#*=}" && shift ;;
            no-enc) parse_options no-enc && shift ;;
            no-dec) parse_options no-dec && shift ;;
            target_system=*) parse_options target_system="${1#*=}" && shift ;;
            shared) parse_options build_shared && shift ;;
            static) parse_options build_static && shift ;;
            release) parse_options release && shift ;;
            test) parse_options tests && shift ;;
            toolchain=*) parse_options toolchain="${1#*=}" && shift ;;
            verbose) parse_options verbose && shift ;;
            end) ${IN_SCRIPT:-false} && exit ;;
            *) die "Error, unknown option: $1" ;;
            esac
        fi
    done
fi

if [ "${PATH#*\/usr\/local\/bin}" = "$PATH" ]; then
    PATH="$PATH:/usr/local/bin"
fi

if $build_debug && $build_release; then
    build release "$@"
    build debug "$@"
elif $build_debug; then
    build debug "$@"
else
    build release "$@"
fi

if $build_install; then
    if $build_release; then
        install_build release
    elif $build_debug; then
        install_build debug
    else
        build release "$@"
        install_build release
    fi
fi
