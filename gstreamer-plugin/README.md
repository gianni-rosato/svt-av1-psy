# GStreamer-SVT-AV1

## Overview

This plugin provides svtav1enc element to GStreamer in order to use the Scalable Video Technology for AV1 Encoder ([SVT-AV1](https://github.com/OpenVisualCloud/SVT-AV1)).

## Requirements

* GStreamer 1.13.1 or later (https://github.com/GStreamer/gst-plugins-good/commit/de8f080a5c04db505dd754820aea1fc479302ee8)
* Scalable Video Technology for AV1 Encoder
  * `SvtAv1Enc.dll` or `libSvtAv1Enc.so` has to be in the `PATH` or next to the plugin's DLL/.so.
* Windows or Linux operating system
* A 64-bit CPU with AVX2 support

## Usage

Make sure that the SvtAv1Enc library is in a path the OS looks for when loading dynamic libraries. If using default install locations, this means for example:

For *nix:

``` bash
export LD_LIBRARY_PATH+=":/usr/local/lib"
```

for Windows:

``` batch
set "PATH=C:\svt-encoders\lib;%PATH%"
```

Then a sample GStreamer pipeline is:

``` bash
gst-launch-1.0 -e videotestsrc ! video/x-raw ! svtav1enc ! matroskamux ! filesink location=out.mkv
```

If you're not familiar with GStreamer, gst-launch-1.0 is part of GStreamer tools, and mpegtsmux is part of GStreamer Bad plugins, `-e` option allows CTRL+C to translate to an EOS (end of stream) signal on the pipeline.

## Compiling and Installing

### Build Dependencies

* GStreamer and GStreamer-plugins-base dev 1.13.1 or later
  * on Debian/Ubuntu: `sudo apt install libgstreamer1.0-dev libgstreamer-plugins-base1.0-dev`
  * on Windows: [gstreamer-1.0-devel-x86_64-*.msi](https://gstreamer.freedesktop.org/data/pkg/windows/)
* Scalable Video Technology for AV1 Encoder ([SVT-AV1](https://github.com/OpenVisualCloud/SVT-AV1))
* meson 0.29 or later
  * install python3 and run `pip3 install meson`
* pkg-config
  * on Debian/Ubuntu: `sudo apt install pkg-config`
  * on Windows, we recommend [pkgconfiglite](https://sourceforge.net/projects/pkgconfiglite/)
* *(optional on Windows)* ninja
  * install python3 and run `pip3 install ninja`
  * or on Ubuntu: `sudo apt install ninja-build`

This plugin uses `meson` build tools and the dependency on SVT-AV1 library is set-up using `pkg-config`.

### Linux specific instructions

Make sure first that SVT-AV1 library is installed and can be found using pkg-config. You can do that using CMake:

``` bash
cmake -P SVT-AV1/Build/linux/Release/Source/Lib/Encoder/cmake_install.cmake
```

Then you can compile and install the plugin the following way:

``` bash
meson -Dprefix=/usr build
ninja -C build
sudo ninja -C build install
```

### Windows specific instructions

Make sure first that SVT-AV1 library is installed and can be found using pkg-config. You can do that using CMake:

``` batch
cmake -P SVT-AV1\Build\Windows\Source\Lib\Encoder\cmake_install.cmake
```

The following commands should be run from a Visual Studio command prompt or another build environment like MinGW, not Windows built-in command prompt.

Specify the path to pkgconfig configuration files for GStreamer and SVT-AV1, and installation prefix to `%GSTREAMER_1_0_ROOT_X86_64%`:

``` batch
set "PKG_CONFIG_PATH=%GSTREAMER_1_0_ROOT_X86_64%lib\pkgconfig;C:\svt-encoders\lib\pkgconfig"
```

Then the plugin can be compiled and installed using Ninja:

``` batch
meson -Dprefix=%GSTREAMER_1_0_ROOT_X86_64% build
ninja -C build
ninja -C build install
```

Or made available as a Visual Studio project:

``` batch
meson -Dprefix=%GSTREAMER_1_0_ROOT_X86_64% build --backend=vs2017
```
