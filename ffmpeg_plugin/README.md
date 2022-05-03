# FFmpeg patches for SVT-AV1

This directory contains patches and backported commits that might be of interest
to users of libsvtav1 as a single patch file.

For the original commits the patch was generated from, see <https://gitlab.com/1480c1/FFmpeg/>

## Changes

### n5.0 - [svt-av1/n5.0](https://gitlab.com/1480c1/FFmpeg/-/tree/svt-av1/n5.0)

- [c33b404885](https://gitlab.com/1480c1/FFmpeg/-/commit/c33b404885): Backport `-svtav1-params:v`
- [1dddb930aa](https://gitlab.com/1480c1/FFmpeg/-/commit/1dddb930aa): Backport `-crf:v`, remove `-rc:v` and instead use `-b:v`, `-maxrate:v`, `-crf:v`, and `-qp:v` to set rc mode
- [50bc872635](https://gitlab.com/1480c1/FFmpeg/-/commit/50bc872635): Backport patch for using aq-mode to determine crf or qp
- [51c0b9e829](https://gitlab.com/1480c1/FFmpeg/-/commit/51c0b9e829): Backport patch for passing color description info
- [e3c4442b24](https://gitlab.com/1480c1/FFmpeg/-/commit/f579c1aca1): Backport patch for parsing svtav1-params last

Additionally, the following patch is included, but is not in upstream:

- b7da10c27b: avcodec/libsvtav1: pass pict_type to library

## How to build

Assuming `$PWD` == the root of your SVT-AV1 clone and you have already built
and installed SVT-AV1 and your `PKG_CONFIG_PATH` environment variable is setup
so that `pkg-config --libs SvtAv1Enc` works properly, this may require exporting
`PKG_CONFIG_PATH` to `/usr/local/lib/pkgconfig` or where your prefix is setup

For n5.0:

```bash
git clone --branch n5.0 https://github.com/FFmpeg/FFmpeg.git
git -C FFmpeg am "$PWD/ffmpeg_plugin/n5.0"/*.patch
```

```bash
mkdir -p ffmpeg-build
(
    cd ffmpeg-build
    ../FFmpeg/configure --enable-libsvtav1 # Append other options as needed
)
make -C ffmpeg-build -j$(($(nproc) + 2))
```

Adapt as needed for depending on your setup

## Sample command lines

Basic ffmpeg line

```bash
ffmpeg -y -i input.mkv -c:v libsvtav1 -crf 30 output.webm
```

FFmpeg line with crf+maxrate for capped CRF

```bash
ffmpeg -y -i input.mkv -c:v libsvtav1 -crf 30 -maxrate 6M output.webm
```

FFmpeg line with svtav1-params setting lp and asm

```bash
ffmpeg -y -i input.mkv -c:v libsvtav1 -crf 30 -svtav1-params lp=4:asm=sse4_1 output.webm
```
