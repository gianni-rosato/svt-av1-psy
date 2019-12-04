# Scalable Video Technology for AV1 Encoder (SVT-AV1 Encoder) User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Sample Application Guide](#sample-application-guide)
    - [Input Video Format](#input-video-format)
    - [Compressed 10-bit format](#compressed-10-bit-format)
    - [Running the encoder](#running-the-encoder)
    - [Sample command lines](#sample-command-lines)
    - [List of all configuration parameters](#list-of-all-configuration-parameters)
3. [Legal Disclaimer](#legal-disclaimer)

## Introduction

This document describes how to use the Scalable Video Technology for AV1 Encoder (SVT-AV1). In particular, this user guide describes how to run the sample application with the respective dynamically linked library.

## Sample Application Guide

This section describes how to run the sample encoder application that uses the SVT-AV1 Encoder library. It describes the input video format, the command line input parameters and the resulting outputs.

### Input Video Format

The SVT-AV1 Encoder supports the following input formats:

_8-bit yuv420p_\
![8-bit yuv420p](8bit_yuv420p.png "8-bit yuv420p")

_10-bit yuv420p10le_\
![10-bit yuv420p10le](10bit_yuv420p.png "10-bit yuv420p10le")

### Compressed 10-bit format

In order to reduce the size of the input original YUV file, the SVT-AV1 Encoder uses a compressed 10-bit format allowing the software to achieve a higher speed and channel density levels. The conversion between the 10-bit yuv420p10le and the compressed 10-bit format is a lossless operation and is performed using the following steps.

#### Unpack the 10-bit picture

This step consists of separating the 10 bit video samples into 8 bit and 2 bit planes so that each 10-bit picture will be represented as two separate pictures as shown in the figure below. As a result of the operation, the 2 least significant bits of the 10 bits will be written into a full byte.

_10-bit yuv420p10le unpacked_\
![10-bit yuv420p10le unpacked](10bit_unpacked.png "10-bit yuv420p10le unpacked")

#### Compress the 2 bit Plane

The unpacking steps separates the 10bits into a group of 8 bits and a group of 2 bits, where the 2 bits are stored in a byte. In this step, every group of consecutive 4 bytes, each containing 2bits from the unpacking step, are compressed into one byte. As a result, each 10bit picture will be represented as two separate pictures as shown in the figure below.

_10-bit yuv420p10le compressed_\
![10-bit yuv420p10le compressed](10bit_packed.png "10-bit yuv420p10le compressed")

#### Unroll the 64x64

Now for a faster read of the samples, every 64x64 block of the 2 bit picture should be written into a one dimensional array. Therefore, the top left 64x64 sample block which is now written into a 16 bytes x 64 bytes after the compression of the 2bit samples, will be written into a 1024 bytes x 1 byte array as shown in the picture below.

_64x64 block after 2 bit compression_\
![64x64 block after 2 bit compression](64x64_after_2bit_compression.png "64x64 block after 2 bit compression")

_64x64 block after unrolling_\
![64x64 block after unrolling](64x64_after_unrolling.png "64x64 block after unrolling")

### Running the encoder

This section describes how to run the sample encoder application `SvtAv1EncApp.exe` (on Windows\*) or `SvtAv1EncApp` (on Linux\*) from the command line, including descriptions of the most commonly used input parameters and outputs.

The sample application typically takes the following command line parameters:

`-c filename` **[Optional]**

A text file that contains encoder parameters such as input file name, quantization parameter etc. Refer to the comments in the Config/Sample.cfg for specific details. The list of encoder parameters are also listed below. Note that command line parameters take precedence over the parameters included in the configuration file when there is a conflict.

`-i filename` **[Required]**

A YUV file (e.g. 8 bit 4:2:0 planar) containing the video sequence that will be encoded. The dimensions of each image are specified by `–w` and `–h` as indicated below.

`-b filename` **[Optional]**

The resulting encoded bit stream file in binary format. If none specified, no output bit stream will be produced by the encoder.

`-w integer` **[Required]**

The width of each input image in units of picture luma pixels, e.g. 1920

`-h integer` **[Required]**]

The height of each input image in units of picture luma pixels, e.g. 1080

`-n integer` **[Optional]**

The number of frames of the sequence to encode. e.g. 100. If the input frame count is larger than the number of frames in the input video, the encoder will loopback to the first frame when it is done.

`-intra-period integer` **[Optional]**

The intra period defines the interval of frames after which you insert an Intra refresh. It is strongly recommended to use (multiple of 8) -1 the closest to 1 second (e.g. 55, 47, 31, 23 should be used for 60, 50, 30, (24 or 25) respectively). When using closed gop (-irefresh-type 2) add 1 to the value above (e.g. 56 instead of 55).

`-rc integer` **[Optional]**

This token sets the bitrate control encoding mode [1: Variable Bitrate, 0: Constant QP]. When `-rc` is set to 1, it is best to match the `–lad` (lookahead distance described in the next section) parameter to the `-intra-period`. When `–rc` is set to 0, a qp value is expected with the use of the `–q` command line option otherwise a default value is assigned (25).

For example, the following command encodes 100 frames of the YUV video sequence into the bin bit stream file. The picture is 1920 luma pixels wide and 1080 pixels high using the `Sample.cfg` configuration. The QP equals 30 and the md5 checksum is not included in the bit stream.

`SvtAv1EncApp.exe -c Sample.cfg -i CrowdRun\_1920x1080.yuv -w 1920 -h 1080 -n 100 -q 30 -intra-period 31 -b CrowdRun\_1920x1080\_qp30.bin`

It should be noted that not all the encoder parameters present in the `Sample.cfg` can be changed using the command line.

### Sample command lines

Here are some sample encode command lines

#### 1 pass fixed QP at maximum speed from 24fps yuv 1920x1080 input
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 -fps 24 -rc 0 -q 30 -enc-mode 8 -b output.ivf`

#### 1 pass VBR 10000 Kbps at medium speed from 24fps yuv 1920x1080 input
`SvtAv1EncApp -i input.yuv -w 1920 -h 1080 -fps 24 -rc 2 -tbr 10000 -enc-mode 5 -b output.ivf`

### List of all configuration parameters

The encoder parameters present in the `Sample.cfg` file are listed in this table below along with their status of support, command line parameter and the range of values that the parameters can take.

| **Configuration file parameter** | **Command line** | **Range** | **Default** | **Description** |
| --- | --- | --- | --- | --- |
| **ChannelNumber** | -nch | [1 - 6] | 1 | Number of encode instances |
| **ConfigFile** | -c | any string | null | Configuration file path |
| **InputFile** | -i | any string | None | Input file path |
| **StreamFile** | -b | any string | null | output bitstream file path |
| **ErrorFile** | -errlog | any string | stderr | error log displaying configuration or encode errors |
| **UseQpFile** | -use-q-file | [0 - 1] | 0 | When set to 1, overwrite the picture qp assignment using qp values in QpFile |
| **QpFile** | -qp-file | any string | Null | Path to qp file |
| **StatReport** | -stat-report | [0 - 1] | 0 | When set to 1, calculate and display PSNR values |
| **StatFile** | -stat-file | any string | Null | Path to statistics file if specified and StatReport is set to 1, per picture statistics are outputted in the file|
| **EncoderMode2p** | -enc-mode-2p | [0 - 8] | 8 | Encoder Preset [0,1,2,3,4,5,6,7,8] 0 = highest quality, 8 = highest speed. Passed to encoder's first pass to use the ME settings of the second pass to achieve better bdRate|
| **InputStatFile** | -input-stat-file | any string | Null | Input stat file for second pass|
| **OutputStatFile** | -output-stat-file | any string | Null | Output stat file for first pass|
| **EncoderMode** | -enc-mode | [0 - 8] | 8 | Encoder Preset [0,1,2,3,4,5,6,7,8] 0 = highest quality, 8 = highest speed |
| **EncoderBitDepth** | -bit-depth | [8 , 10] | 8 | specifies the bit depth of the input video |
| **CompressedTenBitFormat** | -compressed-ten-bit-format | [0 - 1] | 0 | Offline packing of the 2bits: requires two bits packed input (0: OFF, 1: ON) |
| **SourceWidth** | -w | [64 - 4096] | None | Input source width |
| **SourceHeight** | -h | [0 - 2304] | None | Input source height |
| **FrameToBeEncoded** | -n | [0 - 2^64 -1] | 0 | Number of frames to be encoded, if number of frames is > number of frames in file, the encoder will loop to the beginning and continue the encode. Use -1 to not buffer. |
| **BufferedInput** | -nb | [-1, 1 to 2^31 -1] | -1 | number of frames to preload to the RAM before the start of the encode If -nb = 100 and –n 1000 -- > the encoder will encode the first 100 frames of the video 10 times |
| **FrameRate** | -fps | [0 - 2^64 -1] | 25 | If the number is less than 1000, the input frame rate is an integer number between 1 and 60, else the input number is in Q16 format (shifted by 16 bits) [Max allowed is 240 fps] |
| **FrameRateNumerator** | -fps-num | [0 - 2^64 -1] | 0 | Frame rate numerator e.g. 6000 |
| **FrameRateDenominator** | -fps-denom | [0 - 2^64 -1] | 0 | Frame rate denominator e.g. 100 |
| **HierarchicalLevels** | -hierarchical-levels | [3 – 4] | 4 | 0 : Flat4: 5-Level HierarchyMinigop Size = (2^HierarchicalLevels) (e.g. 3 == > 7B pyramid, 4 == > 15B Pyramid) |
| **IntraPeriod** | -intra-period | [-2 - 255] | -2 | Distance Between Intra Frame inserted. -1 denotes no intra update. -2 denotes auto. |
| **IntraRefreshType** | -irefresh-type | [1 – 2] | 1 | 1: CRA (Open GOP)2: IDR (Closed GOP) |
| **TargetBitRate** | -tbr | [1 - 4294967] | 7000 | Target bitrate in kilobits per second when RateControlMode is set to 2, or 3 |
| **QP** | -q | [0 - 63] | 50 | Quantization parameter used when RateControl is set to 0 |
| **RateControlMode** | -rc | [0 - 2] | 0 | 0 = CQP , 1 = VBR , 2 = CVBR |
| **AdaptiveQuantization** | -adaptive-quantization | [0 - 2] | 0 | 0 = OFF , 1 = variance base using segments , 2 = Deltaq pred efficiency (default) |
| **UseDefaultMeHme** | -use-default-me-hme | [0 - 1] | 1 | 0 : Overwrite Default ME HME parameters1 : Use default ME HME parameters, dependent on width and height |
| **HME** | -hme | [0 - 1] | 1 | Enable HME, 0 = OFF, 1 = ON |
| **HMELevel0** | -hme-l0 | [0 - 1] | 1 | Enable HME Level 0 , 0 = OFF, 1 = ON |
| **HMELevel1** | -hme-l1 | [0 - 1] | Depends on input resolution | Enable HME Level 1 , 0 = OFF, 1 = ON |
| **HMELevel2** | -hme-l2 | [0 - 1] | Depends on input resolution | Enable HME Level 2 , 0 = OFF, 1 = ON |
| **InLoopMeFlag** | -in-loop-me | [0 - 1] | Depends on –enc-mode | 0=ME on source samples, 1= ME on recon samples |
| **LocalWarpedMotion** | -local-warp | [0 - 1] | 0 | Enable warped motion use , 0 = OFF, 1 = ON |
| **RDOQ** | -rdoq | [0/1, -1 for auto] | AUTO | Enable RDOQ, 0 = OFF, 1 = ON, -1 = AUTO |
| **ExtBlockFlag** | -ext-block | [0 - 1] | Depends on –enc-mode | Enable the non-square block 0=OFF, 1= ON |
| **ScreenContentMode** | -scm | [0 - 2] | 2 | Enable Screen Content Optimization mode (0: OFF, 1: ON, 2: Content Based Detection) |
| **SearchAreaWidth** | -search-w | [1 - 256] | Depends on input resolution | Search Area in Width |
| **SearchAreaHeight** | -search-h | [1 - 256] | Depends on input resolution | Search Area in Height |
| **NumberHmeSearchRegionInWidth** | -num-hme-w | [1 - 2] | Depends on input resolution | Search Regions in Width |
| **NumberHmeSearchRegionInHeight** | -num-hme-h | [1 - 2] | Depends on input resolution | Search Regions in Height |
| **HmeLevel0TotalSearchAreaWidth** | -hme-tot-l0-w | [1 - 256] | Depends on input resolution | Total HME Level 0 Search Area in Width |
| **HmeLevel0TotalSearchAreaHeight** | -hme-tot-l0-h | [1 - 256] | Depends on input resolution | Total HME Level 1 Search Area in Width |
| **HmeLevel0SearchAreaInWidth** | -hme-l0-w | [1 - 256] | Depends on input resolution | HME Level 0 Search Area in Width for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInWidth, and the sum must equal toHmeLevel0TotalSearchAreaWidth |
| **HmeLevel0SearchAreaInHeight** | -hme-l0-h | [1 - 256] | Depends on input resolution | HME Level 0 Search Area in Height for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInHeight, and the sum must equal toHmeLevel0TotalSearchAreaHeight |
| **HmeLevel1SearchAreaInWidth** | -hme-l1-w | [1 - 256] | Depends on input resolution | HME Level 1 Search Area in Width for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInWidth |
| **HmeLevel1SearchAreaInHeight** | -hme-l1-h | [1 - 256] | Depends on input resolution | HME Level 1 Search Area in Height for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInHeight |
| **HmeLevel2SearchAreaInWidth** | -hme-l2-w | [1 - 256] | Depends on input resolution | HME Level 2 Search Area in Width for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInWidth |
| **HmeLevel2SearchAreaInHeight** | -hme-l2-h | [1 - 256] | Depends on input resolution | HME Level 2 Search Area in Height for each region, separated in spaces, the number of input search areas must equal to NumberHmeSearchRegionInHeight |
| **LookAheadDistance** | -lad | [0 - 120] | 33 | When Rate Control is set to 1 it&#39;s best to set this parameter to be equal to the Intra period value (such is the default set by the encoder) [this value is capped by the encoder to its maximum need e.g. 33 for CQP, 2*fps for rate control] |
| **SceneChangeDetection** | -scd | [0 - 1] | 1 | Enables or disables the scene change detection algorithm |
| **AsmType** | -asm | [0 - 1] | 1 | Assembly instruction set (0: Automatically select lowest assembly instruction set supported, 1: Automatically select highest assembly instruction set supported,) |
| **LogicalProcessorNumber** | -lp | [0, total number of logical processor] | 0 | The number of logical processor which encoder threads run on.Refer to Appendix A.1 |
| **TargetSocket** | -ss | [-1,1] | -1 | For dual socket systems, this can specify which socket the encoder runs on.Refer to Appendix A.1 |
| **ReconFile** | -o | any string | null | Recon file path. Optional output of recon. |
| **TileRow** | -tile-rows | [0-6] | 0 | log2 of tile rows |
| **TileCol** | -tile-columns | [0-6] | 0 | log2 of tile columns |
| **UnrestrictedMotionVector** | -umv | [0-1] | 1 | Enables or disables unrestriced motion vectors, 0 = OFF(motion vectors are constrained within tile boundary), 1 = ON. For MCTS support, set -umv 0 |
| **PaletteMode** | -palette | [0 - 6] | -1 | Enable Palette mode (-1: Auto Mode(ON at level6 when SC is detected), 0: OFF 1: ON Level 1, ...6: ON Level6 ) |
| **OlpdRefinement** | -olpd-refinement | [0 - 1] | -1 | Enable open loop partitioning decision refinement (-1: Auto Mode(ON for M0, no SC, OFF otherwise), 0: OFF 1: ON for M0, error otherwise ) |
| **SquareWeight** | -sqw | 0 for off and any whole number percentage | 100 | Weighting applied to square/h/v shape costs when deciding if a and b shapes could be skipped. Set to 100 for neutral weighting, lesser than 100 for faster encode and BD-Rate loss, and greater than 100 for slower encode and BD-Rate gain|
| **MDStage1PruneClassThreshold** | -mds1p-class-th | 0 for off and any whole number percentage | 100 | Deviation threshold (expressed as a percentage) of an inter-class class pruning mechanism before MD Stage 1 |
| **MDStage1PruneCandThreshold** | -mds1p-cand-th | 0 for off and any whole number percentage | 75 | Deviation threshold (expressed as a percentage) of an intra-class candidate pruning mechanism before MD Stage 1 |
| **MDStage2PruneClassThreshold** | -mds2p-class-th | 0 for off and any whole number percentage | 25 | Deviation threshold (expressed as a percentage) of an inter-class class pruning mechanism before MD Stage 2 |
| **MDStage2PruneCandThreshold** | -mds2p-cand-th | 0 for off and any whole number percentage | 15 | Deviation threshold (expressed as a percentage) of an intra-class candidate pruning mechanism before MD Stage 2 |
## Appendix A Encoder Parameters

### 1. Thread management parameters

LogicalProcessorNumber (`-lp`) and TargetSocket (`-ss`) parameters are used to management thread affinity on Windows and Ubuntu OS. These are some examples how you use them together.

If LogicalProcessorNumber and TargetSocket are not set, threads are managed by OS thread scheduler.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 –lp 40`

If only LogicalProcessorNumber is set, threads run on 40 logical processors. Threads may run on dual sockets if 40 is larger than logical processor number of a socket.

NOTE: On Windows, thread affinity can be set only by group on system with more than 64 logical processors. So, if 40 is larger than logical processor number of a single socket, threads run on all logical processors of both sockets.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 –ss 1`

If only TargetSocket is set, threads run on all the logical processors of socket 1.

`SvtAv1EncApp.exe -i in.yuv -w 3840 -h 2160 –lp 20 –ss 0`

If both LogicalProcessorNumber and TargetSocket are set, threads run on 20 logical processors of socket 0. Threads guaranteed to run only on socket 0 if 20 is larger than logical processor number of socket 0.

## Legal Disclaimer

### Optimization Notice

Intel compilers may or may not optimize to the same degree for non-Intel microprocessors for optimizations that are not unique to Intel microprocessors. These optimizations include SSE2, SSE3, and SSSE3 instruction sets and other optimizations. Intel does not guarantee the availability, functionality, or effectiveness of any optimization on microprocessors not manufactured by Intel. Microprocessor-dependent optimizations in this product are intended for use with Intel microprocessors. Certain optimizations not specific to Intel microarchitecture are reserved for Intel microprocessors. Please refer to the applicable product User and Reference Guides for more information regarding the specific instruction sets covered by this notice.

### Notice Revision #20110804

Intel technologies features and benefits depend on system configuration and may require enabled hardware, software or service activation. Performance varies depending on system configuration. No computer system can be absolutely secure. Check with your system manufacturer or retailer.

No license (express or implied, by estoppel or otherwise) to any intellectual property rights is granted by this document.

Intel disclaims all express and implied warranties, including without limitation, the implied warranties of merchantability, fitness for a particular purpose, and non-infringement, as well as any warranty arising from course of performance, course of dealing, or usage in trade.

The products and services described may contain defects or errors known as errata which may cause deviations from published specifications. Current characterized errata are available on request.  ** ** No product or component can be absolutely secure.

This document contains information on products, services and/or processes in development. All information provided here is subject to change without notice. Contact your Intel representative to obtain the latest forecast, schedule, specifications and roadmaps.

Intel, Intel Xeon, Intel Core, the Intel logo and others are trademarks of Intel Corporation and its subsidiaries in the U.S. and/or other countries.

\*Other names and brands may be claimed as the property of others.

Copyright 2019 Intel Corporation.
