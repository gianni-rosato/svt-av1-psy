# Scalable Video Technology for AV1 Encoder (SVT-AV1 Encoder) User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Sample Application Guide](#sample-application-guide)
    - [Input Video Format](#input-video-format)
    - [Compressed 10-bit format](#compressed-10-bit-format)
    - [Running the encoder](#running-the-encoder)
    - [Sample command lines](#sample-command-lines)
    - [List of all configuration parameters](#list-of-all-configuration-parameters)

## Introduction

This document describes how to use the Scalable Video Technology for AV1 Encoder (SVT-AV1). In particular, this user guide describes how to run the sample application with the respective dynamically linked library.

## Sample Application Guide

This section describes how to run the sample encoder application that uses the SVT-AV1 Encoder library. It describes the input video format, the command line input parameters and the resulting outputs.

### Input Video Format

The SVT-AV1 Encoder supports the following input formats:

_8-bit yuv420p_\
![8-bit yuv420p](img/8bit_yuv420p.png "8-bit yuv420p")

_10-bit yuv420p10le_\
![10-bit yuv420p10le](img/10bit_yuv420p.png "10-bit yuv420p10le")

### Compressed 10-bit format

In order to reduce the size of the input original YUV file, the SVT-AV1 Encoder uses a compressed 10-bit format allowing the software to achieve a higher speed and channel density levels. The conversion between the 10-bit yuv420p10le and the compressed 10-bit format is a lossless operation and is performed using the following steps.

#### Unpack the 10-bit picture

This step consists of separating the 10 bit video samples into 8 bit and 2 bit planes so that each 10-bit picture will be represented as two separate pictures as shown in the figure below. As a result of the operation, the 2 least significant bits of the 10 bits will be written into a full byte.

_10-bit yuv420p10le unpacked_\
![10-bit yuv420p10le unpacked](img/10bit_unpacked.png "10-bit yuv420p10le unpacked")

#### Compress the 2 bit Plane

The unpacking steps separates the 10bits into a group of 8 bits and a group of 2 bits, where the 2 bits are stored in a byte. In this step, every group of consecutive 4 bytes, each containing 2bits from the unpacking step, are compressed into one byte. As a result, each 10bit picture will be represented as two separate pictures as shown in the figure below.

_10-bit yuv420p10le compressed_\
![10-bit yuv420p10le compressed](img/10bit_packed.png "10-bit yuv420p10le compressed")

#### Unroll the 64x64

Now for a faster read of the samples, every 64x64 block of the 2 bit picture should be written into a one dimensional array. Therefore, the top left 64x64 sample block which is now written into a 16 bytes x 64 bytes after the compression of the 2bit samples, will be written into a 1024 bytes x 1 byte array as shown in the picture below.

_64x64 block after 2 bit compression_\
![64x64 block after 2 bit compression](img/64x64_after_2bit_compression.png "64x64 block after 2 bit compression")

_64x64 block after unrolling_\
![64x64 block after unrolling](img/64x64_after_unrolling.png "64x64 block after unrolling")

### Running the encoder

This section describes how to run the sample encoder application `SvtAv1EncApp.exe` (on Windows\*) or `SvtAv1EncApp` (on Linux\*) from the command line, including descriptions of the most commonly used input parameters and outputs.

The sample application typically takes the following command line parameters:

`-c filename` **[Optional]**

A text file that contains encoder parameters such as input file name, quantization parameter etc. Refer to the comments in the Config/Sample.cfg for specific details. The list of encoder parameters are also listed below. Note that command line parameters take precedence over the parameters included in the configuration file when there is a conflict.

`-i filename` **[Required]**

A YUV file (e.g. 8 bit 4:2:0 planar) containing the video sequence that will be encoded. The dimensions of each image are specified by `-w` and `-h` as indicated below.

`-b filename` **[Optional]**

The resulting encoded bit stream file in binary format. If none specified, no output bit stream will be produced by the encoder.

`-w integer` **[Required]**

The width of each input image in units of picture luma pixels, e.g. 1920

`-h integer` **[Required]**]

The height of each input image in units of picture luma pixels, e.g. 1080

`-n integer` **[Optional]**

The number of frames of the sequence to encode. e.g. 100. If the input frame count is larger than the number of frames in the input video, the encoder will loopback to the first frame when it is done.

`--keyint integer` **[Optional]**

The intra period defines the interval of frames after which you insert an Intra refresh. It is strongly recommended to use (multiple of 8) -1 the closest to 1 second (e.g. 55, 47, 31, 23 should be used for 60, 50, 30, (24 or 25) respectively). When using closed gop (-irefresh-type 2) add 1 to the value above (e.g. 56 instead of 55).

`--rc integer` **[Optional]**

This token sets the bitrate control encoding mode [1: Variable Bitrate, 0: Constant QP OR Constant Rate Factor]. When `--rc` is set to 1.

With `--rc` set to 0, if `--crf` is used then enable-tpl-la is forced to 1, however, if `-q`/`--qp` is used then the encoder will work in CRF mode if `--enable-tpl-la` is set to 1 and in CQP mode (fixed qp offsets regardless of the content) when `--enable-tpl-la` is set to 0.

If a qp/crf value is not specified, a default value is assigned (50).

For example, the following command encodes 100 frames of the YUV video sequence into the bin bit stream file. The picture is 1920 luma pixels wide and 1080 pixels high using the `Sample.cfg` configuration. The QP equals 30 and the md5 checksum is not included in the bit stream.

`SvtAv1EncApp.exe -c Sample.cfg -i CrowdRun_1920x1080.yuv -w 1920 -h 1080 -n 100 -q 30 --keyint 31 -b CrowdRun_1920x1080_qp30.bin`

It should be noted that not all the encoder parameters present in the `Sample.cfg` can be changed using the command line.

### Sample command lines

Here are some sample encode command lines

#### 1 pass CRF at maximum speed from 24fps yuv 1920x1080 input
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --crf 30 --preset 8 -b output.ivf`

#### 1 pass VBR 1000 Kbps at medium speed from 24fps yuv 1920x1080 input
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --rc 1 --tbr 1000 --preset 5 -b output.ivf`

#### 2 pass VBR 1000 Kbps at maximum quality from 24fps yuv 1920x1080 input
1 command line :

`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --rc 1 --tbr 1000 --preset 0 --irefresh-type 2 --passes 2 --stats stat_file.stat -b output.ivf`

or

2 command lines :

`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --rc 1 --tbr 1000 --preset 8 --irefresh-type 2 --pass 1 --stats stat_file.stat`
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --rc 1 --tbr 1000 --preset 0 --irefresh-type 2 --pass 2 --stats stat_file.stat -b output.ivf`

#### 1 pass CRF at maximum speed from 24fps yuv 1920x1080 input with full range video signal
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --crf 30 --preset 8 --color-range 1 -b output.ivf`

#### 1 pass CRF at maximum speed from 24fps yuv 1920x1080 input with colorimetry set to BT.709
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 --fps 24 --crf 30 --preset 8 --color-primaries 1 --transfer-characteristics 1 --matrix-coefficients 1 -b output.ivf`

### List of all configuration parameters

The encoder parameters present in the `Sample.cfg` file are listed in this table below along with their status of support, command line parameter and the range of values that the parameters can take.

#### Options

| **Configuration file parameter** | **Command line**   | **Range**  | **Default** | **Description**                                                                                                 |
|----------------------------------|--------------------|------------|-------------|-----------------------------------------------------------------------------------------------------------------|
|                                  | --help             |            |             | Shows the command line options currently available                                                              |
|                                  | --version          |            |             | Shows the version of the library that's linked to the library                                                   |
| **InputFile**                    | -i                 | any string | None        | Input raw video (y4m and yuv) file path, use `stdin` to read from pipe                                          |
| **StreamFile**                   | -b                 | any string | None        | Output compressed (ivf) file path, use `stdout` to write to pipe                                                |
|                                  | -c                 | any string | None        | Configuration file path                                                                                         |
| **ErrorFile**                    | --errlog           | any string | `stderr`    | Error file path                                                                                                 |
| **ReconFile**                    | -o                 | any string | None        | Reconstructed yuv file path                                                                                     |
| **StatFile**                     | --stat-file        | any string | None        | PSNR / SSIM per picture stat output file path, requires `--enable-stat-report 1`                                |
| **PredStructFile**               | --pred-struct-file | any string | None        | Manual prediction structure file path                                                                           |
| **Progress**                     | --progress         | [0-2]      | 1           | Verbosity of the output [0: no progress is printed, 2: aomenc style output]                                     |
| **NoProgress**                   | --no-progress      | [0-1]      | 0           | Do not print out progress [1: `--progress 0`, 0: `--progress 1`]                                                |
| **EncoderMode**                  | --preset           | [-2-13]    | 12          | Encoder preset, presets < 0 are for debugging. Higher presets means faster encodes, but with a quality tradeoff |
| **SvtAv1Params**                 | --svtav1-params    | any string | None        | Colon-separated list of `key=value` pairs of parameters with keys based on command line options without `--`    |
|                                  | --nch              | [1-6]      | 1           | Number of channels (library instance) that will be instantiated                                                 |

##### Usage of **SvtAv1Params**

To use the `--svtav1-params` option, the syntax is `--svtav1-params option1=value1:option2=value2...`.

An example is:

```bash
SvtAv1EncApp \
  -i input.y4m \
  -b output.ivf \
  --svtav1-params \
  "preset=10:crf=30:irefresh-type=kf:matrix-coefficients=bt709:mastering-display=G(0.2649,0.6900)B(0.1500,0.0600)R(0.6800,0.3200)WP(0.3127,0.3290)L(1000.0,1)"
```

This will set `--preset` to 10 and `--crf` to 30 inside the API along with some other parameters.

Do note however, that there is no error checking for duplicate keys and only for invalid keys or values.

For more information on valid values for specific keys, refer to the [EbEncSettings](../Source/Lib/Encoder/Globals/EbEncSettings.c) file.

#### Encoder Global Options

| **Configuration file parameter** | **Command line**            | **Range**                      | **Default** | **Description**                                                                                               |
|----------------------------------|-----------------------------|--------------------------------|-------------|---------------------------------------------------------------------------------------------------------------|
| **SourceWidth**                  | -w                          | [64-16384]                     | None        | Frame width in pixels, inferred if y4m.                                                                       |
| **SourceHeight**                 | -h                          | [64-8704]                      | None        | Frame height in pixels, inferred if y4m.                                                                      |
| **FrameToBeEncoded**             | -n                          | [0-`(2^63)-1`]                 | 0           | Number of frames to encode. If `n` is larger than the input, the encoder will loop back and continue encoding |
| **BufferedInput**                | --nb                        | [-1, 1-`(2^31)-1`]             | -1          | Buffer `n` input frames into memory and use them to encode                                                    |
| **EncoderColorFormat**           | --color-format              | [0-3]                          | 1           | Color format, only yuv420 is supported at this time [0: yuv400, 1: yuv420, 2: yuv422, 3: yuv444]              |
| **Profile**                      | --profile                   | [0-2]                          | 0           | Bitstream profile [0: main, 1: high, 2: professional]                                                         |
| **Level**                        | --level                     | [0,2.0-7.3]                    | 0           | Bitstream level, defined in A.3 of the av1 spec [0: auto]                                                     |
| **HighDynamicRangeInput**        | --enable-hdr                | [0-1]                          | 0           | Enable writing of HDR metadata in the bitstream                                                               |
| **FrameRate**                    | --fps                       | [1-240]                        | 25          | Input video frame rate, integer values only, inferred if y4m                                                  |
| **FrameRateNumerator**           | --fps-num                   | [0-2^32-1]                     | 25000       | Input video frame rate numerator                                                                              |
| **FrameRateDenominator**         | --fps-denom                 | [0-2^32-1]                     | 1000        | Input video frame rate denominator                                                                            |
| **EncoderBitDepth**              | --input-depth               | [8, 10]                        | 8           | Input video file and output bitstream bit-depth                                                               |
| **CompressedTenBitFormat**       | --compressed-ten-bit-format | [0-1]                          | 0           | Pack 10bit video, handled between the app and library                                                         |
| **Injector**                     | --inj                       | [0-1]                          | 0           | Inject pictures to the library at defined frame rate                                                          |
| **InjectorFrameRate**            | --inj-frm-rt                | [0-240]                        | 60          | Set injector frame rate, only applicable with `--inj 1`                                                       |
| **StatReport**                   | --enable-stat-report        | [0-1]                          | 0           | Calculates and outputs PSNR SSIM metrics at the end of encoding                                               |
| **Asm**                          | --asm                       | [0-11, c-max]                  | max         | Limit assembly instruction set [c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512, max]       |
| **LogicalProcessors**            | --lp                        | [0, core count of the machine] | 0           | Target (best effort) number of logical cores to be used. 0 means all. Refer to Appendix A.1                   |
| **PinnedExecution**              | --pin                       | [0-1]                          | 0           | Pin the execution to the first --lp cores. Overwritten to 0 when `--ss` is set. Refer to Appendix A.1         |
| **TargetSocket**                 | --ss                        | [-1,1]                         | -1          | Specifies which socket to run on, assumes a max of two sockets. Refer to Appendix A.1                         |
| **FastDecode**                   | --fast-decode               | [0,3]                          | 0           | Tune settings to output bitstreams that can be decoded faster, higher values for faster decoding              |
| **Tune**                         | --tune                      | [0,1]                          | 1           | Specifies whether to use PSNR or VQ as the tuning metric [0 = VQ, 1 = PSNR]                                   |

#### Rate Control Options

| **Configuration file parameter** | **Command line**                 | **Range**      | **Default**     | **Description**                                                                                                      |
|----------------------------------|----------------------------------|----------------|-----------------|----------------------------------------------------------------------------------------------------------------------|
| **RateControlMode**              | --rc                             | [0-2]          | 0               | Rate control mode [0: CRF or CQP (if `--enable-tpl-la` is 0) [Default], 1: VBR, 2: CBR]                              |
| **QP**                           | --qp                             | [1-63]         | 50              | Initial QP level value                                                                                               |
| **CRF**                          | --crf                            | [1-63]         | 50              | Constant Rate Factor value, setting this value is equal to `--rc 0 --enable-tpl-la 1 --qp x`                         |
| **TargetBitRate**                | --tbr                            | [1-4294967]    | 2000            | Target Bitrate (kbps), only applicable for VBR and CBR encoding                                                      |
| **MaxBitRate**                   | --mbr                            | [1-4294967]    | 0               | Maximum Bitrate (kbps) only applicable for CRF and VBR encoding                                                      |
| **UseQpFile**                    | --use-q-file                     | [0-1]          | 0               | Overwrite the encoder default picture based QP assignments and use QP values from `--qp-file`                        |
| **QpFile**                       | --qpfile                         | any string     | Null            | Path to a file containing per picture QP value                                                                       |
| **MaxQpAllowed**                 | --max-qp                         | [1-63]         | 63              | Maximum (highest) quantizer, only applicable for VBR and CBR                                                         |
| **MinQpAllowed**                 | --min-qp                         | [1-63]         | 1               | Minimum (lowest) quantizer, only applicable for VBR and CBR                                                          |
| **AdaptiveQuantization**         | --aq-mode                        | [0-2]          | 2               | Set adaptive QP level [0: off, 1: variance base using AV1 segments, 2: deltaq pred efficiency]                       |
| **VBVBufSize**                   | --vbv-bufsize                    | [1-4294967]    | `TargetBitRate` | VBV buffer size.                                                                                                     |
| **UseFixedQIndexOffsets**        | --use-fixed-qindex-offsets       | [0-1]          | 0               | Overwrite the encoder default hierarchical layer based QP assignment and use fixed Q index offsets                   |
| **KeyFrameQIndexOffset**         | --key-frame-qindex-offset        | [-256-255]     | 0               | Overwrite the encoder default keyframe Q index assignment                                                            |
| **KeyFrameChromaQIndexOffset**   | --key-frame-chroma-qindex-offset | [-256-255]     | 0               | Overwrite the encoder default chroma keyframe Q index assignment                                                     |
| **QIndexOffsets**                | --qindex-offsets                 | any string     | `0,0,..,0`      | list of luma Q index offsets per hierarchical layer, separated by `,` with each offset in the range of [-256-255]    |
| **ChromaQIndexOffsets**          | --chroma-qindex-offsets          | any string     | `0,0,..,0`      | list of chroma Q index offsets per hierarchical layer, separated by `,` with each offset in the range of [-256-255]  |
| **UnderShootPct**                | --undershoot-pct                 | [0-100]        | 25              | Allowable datarate undershoot (min) target (%), default depends on the rate control mode                             |
| **OverShootPct**                 | --overshoot-pct                  | [0-100]        | 25              | Allowable datarate overshoot (max) target (%), default depends on the rate control mode                              |
| **BufSz**                        | --buf-sz                         | [0-`(2^63)-1`] | 6000            | Client buffer size (ms), only applicable for CBR                                                                     |
| **BufInitialSz**                 | --buf-initial-sz                 | [0-`(2^63)-1`] | 4000            | Client initial buffer size (ms), only applicable for CBR                                                             |
| **BufOptimalSz**                 | --buf-optimal-sz                 | [0-`(2^63)-1`] | 5000            | Client optimal buffer size (ms), only applicable for CBR                                                             |
| **RecodeLoop**                   | --recode-loop                    | [0-4]          | 4               | Recode loop level, look at the "Recode loop level table" in the user's guide for more info [0: off, 4: preset based] |
| **VBRBiasPct**                   | --bias-pct                       | [0-100]        | 50              | CBR/VBR bias [0: CBR-like, 100: VBR-like]                                                                            |
| **MinSectionPct**                | --minsection-pct                 | [0-`(2^32)-1`] | 0               | GOP min bitrate (expressed as a percentage of the target rate)                                                       |
| **MaxSectionPct**                | --maxsection-pct                 | [0-`(2^32)-1`] | 2000            | GOP max bitrate (expressed as a percentage of the target rate)                                                       |

##### **UseFixedQIndexOffsets** and more information

`UseFixedQIndexOffsets` and its associated arguments (`HierarchicalLevels`, `QIndexOffsets`, `ChromaQIndexOffsets`, `KeyFrameQIndexOffset`, `KeyFrameChromaQIndexOffset`)
are used together to specify the qindex offsets based on frame type and temporal layer when rc is set to 0.

QP value specified by the `--qp` argument is assigned to the pictures at the highest temporal layer.
It is first converted to a qindex, then the corresponding qindex offsets are added on top of it based on the frame types (Key/Inter) and temporal layer id.

Qindex offset can be negative.
The final qindex value will be clamped within the valid min/max qindex range.

For chroma plane, after deciding the qindex for the luma plane, the corresponding chroma qindex offsets are added on top of the luma plane qindex based on frame types and temporal layer id.

`--qindex-offsets` and `--chroma-qindex-offsets` have to be used after the `--hierarchical-levels` parameter.
The number of qindex offsets should be `HierarchicalLevels` plus 1, and they can be enclosed in `[]` to separate the list.

An example command line is:

```bash
SvtAv1EncApp -i in.y4m -b out.ivf --rc 0 -q 42 --hierarchical-levels 3 --use-fixed-qindex-offsets 1 --qindex-offsets [-12,-8,-4,0] --key-frame-qindex-offset -20 --key-frame-chroma-qindex-offset -6 --chroma-qindex-offsets [-6,0,12,24]
```

For this command line, corresponding qindex values are:

| **Frame Type**   | **Luma qindex** | **Chroma qindex** |
|------------------|-----------------|-------------------|
| **Key Frame**    | 148 (42x4 - 20) | 142 (148 - 6)     |
| **Layer0 Frame** | 156 (42x4 - 12) | 150 (156 - 6)     |
| **Layer1 Frame** | 160 (42x4 - 8)  | 160 (160 + 0)     |
| **Layer2 Frame** | 164 (42x4 - 4)  | 176 (164 + 12)    |
| **Layer3 Frame** | 168 (42x4 + 0)  | 192 (168 + 24)    |

##### Recode loop level table

| level | description                                                                     |
|-------|---------------------------------------------------------------------------------|
| 0     | Off                                                                             |
| 1     | Allow recode for KF and exceeding maximum frame bandwidth                       |
| 2     | Allow recode only for key frames, alternate reference frames, and Golden frames |
| 3     | Allow recode for all frame types based on bitrate constraints                   |
| 4     | Preset based decision                                                           |


#### Multi-pass Options

| **Configuration file parameter** | **Command line** | **Range**      | **Default**        | **Description**                                                                                   |
|----------------------------------|------------------|----------------|--------------------|---------------------------------------------------------------------------------------------------|
| **Pass**                         | --pass           | [0-3]          | 0                  | Multi-pass selection [0: single pass encode, 1: first pass, 2: second pass, 3: third pass]        |
| **Stats**                        | --stats          | any string     | "svtav1_2pass.log" | Filename for multi-pass encoding                                                                  |
| **Passes**                       | --passes         | [1-2]          | 1                  | Number of encoding passes, default is preset dependent [1: one pass encode, 2: multi-pass encode] |

##### **Pass** information

| **Pass** | **Stats** io            |
|----------|-------------------------|
| 0        | ""                      |
| 1        | "w"                     |
| 2        | "rw" if 3-pass else "r" |
| 3        | "r"                     |

`--pass 3` is only available for non-crf modes and all passes except single-pass requires the `--stats` parameter to point to a valid path

#### GOP size and type Options

| **Configuration file parameter** | **Command line**      | **Range**       | **Default** | **Description**                                                                                                 |
|----------------------------------|-----------------------|-----------------|-------------|-----------------------------------------------------------------------------------------------------------------|
| **Keyint**                       | --keyint              | [-2-`(2^31)-1`] | -2          | GOP size (frames) [-2: ~2 seconds, -1: "infinite" and only applicable for CRF, 0: same as -1]                   |
| **IntraRefreshType**             | --irefresh-type       | [1-2]           | 2           | Intra refresh type [1: FWD Frame (Open GOP), 2: KEY Frame (Closed GOP)]                                         |
| **SceneChangeDetection**         | --scd                 | [0-1]           | 0           | Scene change detection control                                                                                  |
| **Lookahead**                    | --lookahead           | [-1,0-120]      | -1          | Number of frames in the future to look ahead, beyond minigop, temporal filtering, and rate control [-1: auto]   |
| **HierarchicalLevels**           | --hierarchical-levels | [3-5]           | 4           | Set hierarchical levels beyond the base layer [3: 4 temporal layers, 5: 6 temporal layers]                      |
| **PredStructure**                | --pred-struct         | [0-2]           | 2           | Set prediction structure [0: low delay P-frames, 1: low delay B-frames, 2: random access]                       |

#### AV1 Specific Options

| **Configuration file parameter** | **Command line**     | **Range** | **Default** | **Description**                                                                                                           |
|----------------------------------|----------------------|-----------|-------------|---------------------------------------------------------------------------------------------------------------------------|
| **TileRow**                      | --tile-rows          | [0-6]     | 0           | Number of tile rows to use, `TileRow == log2(x)`, default changes per resolution                                          |
| **TileCol**                      | --tile-columns       | [0-4]     | 0           | Number of tile columns to use, `TileCol == log2(x)`, default changes per resolution                                       |
| **LoopFilterEnable**             | --enable-dlf         | [0-1]     | 1           | Deblocking loop filter control                                                                                            |
| **CDEFLevel**                    | --enable-cdef        | [0-1]     | 1           | Enable Constrained Directional Enhancement Filter                                                                         |
| **EnableRestoration**            | --enable-restoration | [0-1]     | 1           | Enable loop restoration filter                                                                                            |
| **EnableTPLModel**               | --enable-tpl-la      | [0-1]     | 1           | Temporal Dependency model control, currently forced on library side, only applicable for CRF/CQP                          |
| **Mfmv**                         | --enable-mfmv        | [-1-1]    | -1          | Motion Field Motion Vector control [-1: auto]                                                                             |
| **EnableTF**                     | --enable-tf          | [0-1]     | 1           | Enable ALT-REF (temporally filtered) frames                                                                               |
| **EnableOverlays**               | --enable-overlays    | [0-1]     | 0           | Enable the insertion of overlayer pictures which will be used as an additional reference frame for the base layer picture |
| **ScreenContentMode**            | --scm                | [0-2]     | 2           | Set screen content detection level [0: off, 1: on, 2: content adaptive]                                                   |
| **RestrictedMotionVector**       | --rmv                | [0-1]     | 0           | Restrict motion vectors from reaching outside the picture boundary                                                        |
| **FilmGrain**                    | --film-grain         | [0-50]    | 0           | Enable film grain [0: off, 1-50: level of denoising for film grain]                                                       |
| **SuperresMode**                 | --superres-mode      | [0-4]     | 0           | Enable super-resolution mode, refer to the super-resolution section below for more info                                   |
| **SuperresDenom**                | --superres-denom     | [8-16]    | 8           | Super-resolution denominator, only applicable for mode == 1 [8: no scaling, 16: half-scaling]                             |
| **SuperresKfDenom**              | --superres-kf-denom  | [8-16]    | 8           | Super-resolution denominator for key frames, only applicable for mode == 1 [8: no scaling, 16: half-scaling]              |
| **SuperresQthres**               | --superres-qthres    | [0-63]    | 43          | Super-resolution q-threshold, only applicable for mode == 3                                                               |
| **SuperresKfQthres**             | --superres-kf-qthres | [0-63]    | 43          | Super-resolution q-threshold for key frames, only applicable for mode == 3                                                |

##### **Super-Resolution**

Super resolution is better described in [the Super-Resolution documentation](./Appendix-Super-Resolution.md), but
this basically allows the input to be encoded at a lower resolution, horizontally, but then later upscaled back to
the original resolution by the decoder.

| **SuperresMode** | **Value**                                                                                                                   |
|------------------|-----------------------------------------------------------------------------------------------------------------------------|
| 0                | None, no frame super-resolution allowed                                                                                     |
| 1                | All frames are encoded at the specified scale of 8/`denom`, thus a `denom` of 8 means no scaling, and 16 means half-scaling |
| 2                | All frames are coded at a random scale                                                                                      |
| 3                | Super-resolution scale for a frame is determined based on the q_index, a qthreshold of 63 means no scaling                  |
| 4                | Automatically select the super-resolution mode for appropriate frames                                                       |

The performance of the encoder will be affected for all modes other than mode 0. And for mode 4, it should be noted that
the encoder will run at least twice, one for down scaling, and another with no scaling, and then it will choose the best
one for each of the appropriate frames.

For more information on the decision-making process,
please look at [section 2.2 of the super-resolution doc](./Appendix-Super-Resolution.md#22-determination-of-the-downscaling-factor)

#### Color Description Options

| **Configuration file parameter** | **Command line**           | **Range**  | **Default** | **Description**                                                                                                                          |
|----------------------------------|----------------------------|------------|-------------|------------------------------------------------------------------------------------------------------------------------------------------|
| **ColorPrimaries**               | --color-primaries          | [0-12, 22] | 2           | Color primaries, refer to the user guide Appendix A.2 for full details                                                                   |
| **TransferCharacteristics**      | --transfer-characteristics | [0-22]     | 2           | Transfer characteristics, refer to the user guide Appendix A.2 for full details                                                          |
| **MatrixCoefficients**           | --matrix-coefficients      | [0-14]     | 2           | Matrix coefficients, refer to the user guide Appendix A.2 for full details                                                               |
| **ColorRange**                   | --color-range              | [0-1]      | 0           | Color range [0: Studio, 1: Full]                                                                                                         |
| **MasteringDisplay**             | --mastering-display        | any string | none        | Mastering display metadata in the format of "G(x,y)B(x,y)R(x,y)WP(x,y)L(max,min)", refer to the user guide Appendix A.2 for full details |
| **ContentLightLevel**            | --content-light            | any string | none        | Set content light level in the format of "max_cll,max_fall", refer to the user guide Appendix A.2 for full details                       |

## Appendix A Encoder Parameters

### 1. Thread management parameters

`LogicalProcessors` (`--lp`) and `TargetSocket` (`--ss`) parameters are used to management thread affinity on Windows and Ubuntu OS. These are some examples how you use them together.

If `LogicalProcessors` and `TargetSocket` are not set, threads are managed by OS thread scheduler.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 --lp 40`

If only `LogicalProcessors` is set, threads run on 40 logical processors. Threads may run on dual sockets if 40 is larger than logical processor number of a socket.

NOTE: On Windows, thread affinity can be set only by group on system with more than 64 logical processors. So, if 40 is larger than logical processor number of a single socket, threads run on all logical processors of both sockets.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 --ss 1`

If only `TargetSocket` is set, threads run on all the logical processors of socket 1.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 --lp 20 --ss 0`

If both `LogicalProcessors` and `TargetSocket` are set, threads run on 20 logical processors of socket 0. Threads guaranteed to run only on socket 0 if 20 is larger than logical processor number of socket 0.

The (`--pin`) option allows the user to pin/unpin the execution to/from a specific number of cores.

The combinational use of (`--pin`)  with (`--lp`) results in memory reduction while allowing the execution to work on any of the cores and not restrict it to specific cores.

This is an example on how to use them together.

so -lp 4 with --pin 1 would restrict the encoder to work on cpu0-3 and reduce the resource allocation to only what's needed to using 4 cores. --lp 4 with --pin 1, would reduce the allocation to what's needed for 4 cores but not restrict the encoder to run on cpu 0-3, in this case the encoder might end up using more than 4 cores due to the multi-threading nature of the encoder, but would at least allow for more multiple -lp4 encodes to run on the same machine without them being all restricted to run on cpu 0-3 or overflow the memory usage.

Example: 72 core machine:

72 jobs x --lp 1 --pin 0 (In order to maximize the CPU utilization 72 jobs are run simultaneously with each job utilitizing 1 core without being pined to a specific core)

36 jobs x --lp 2 --pin 1

18 jobs x --lp 4 --pin 1

(`--ss`) and (`--pin 0`) is not a valid combination.(`--pin`) is overwritten to 1 when (`-ss`) is used.

### 2. AV1 metadata

Please see the subsection 6.4.2, 6.7.3, and 6.7.4 of the [AV1 Bitstream & Decoding Process Specification](https://aomediacodec.github.io/av1-spec/av1-spec.pdf) for more details on some expected values.

`MasteringDisplay` (`--mastering-display`) and `ContentLightLevel` (`--content-light`) parameters are used to set the mastering display and content light level in the AV1 bitstream.

`MasteringDisplay` takes the format of `G(x,y)B(x,y)R(x,y)WP(x,y)L(max,min)` where

- `G(x,y)` is the green channel of the mastering display
- `B(x,y)` is the blue channel of the mastering display
- `R(x,y)` is the red channel of the mastering display
- `WP(x,y)` is the white point of the mastering display
- `L(max,min)` is the light level of the mastering display

The `x` and `y` values can be coordinates from 0.0 to 1.0, as specified in CIE 1931 while the min,max values can be floating point values representing candelas per square meter, or nits.
For the `max,min` values, they are generally specified in the range of 0.0 to 1.0, but there are no constraints on the provided values.
Invalid values will be clipped accordingly.

`ContentLightLevel` takes the format of `max_cll,max_fall` where both values are integers clipped into a range of 0 to 65535.

Examples:

```bash
SvtAv1EncApp -i int.y4m -b out.ivf \
    --mastering-display "G(0.2649,0.6900)B(0.1500,0.0600)R(0.6800,0.3200)WP(0.3127,0.3290)L(1000.0,1)" \
    --content-light 100,50 \
    --color-primaries 9 \
    --transfer-characteristics 16 \
    --matrix-coefficients 9
    # Color primary 9 is BT.2020, BT.2100
    # Transfer characteristic 16 is SMPTE ST 2084, ITU BT.2100 PQ
    # matrix coefficients 9 is BT.2020 non-constant luminance, BT.2100 YCbCr

# or

ffmpeg -i in.mp4 -strict -1 -f yuv4mpegpipe - |
  SvtAv1EncApp -i stdin -b stdout \
    --mastering-display "G(0.2649,0.6900)B(0.1500,0.0600)R(0.6800,0.3200)WP(0.3127,0.3290)L(1000.0,1)" \
    --content-light 100,50 \
    --color-primaries 9 \
    --transfer-characteristics 16 \
    --matrix-coefficients 9 |
  ffmpeg -y -i - -i audio.ogg -c copy out.mp4
```
