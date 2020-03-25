# Scalable Video Technology for AV1 Decoder (SVT-AV1 Decoder) User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Sample Application Guide](#sample-application-guide)
    - [Running the decoder](#running-the-decoder)
3. [Legal Disclaimer](#legal-disclaimer)

## Introduction

This document describes how to use the Scalable Video Technology for AV1 Decoder (SVT-AV1). In particular, this user guide describes how to run the sample application with the respective dynamically linked library.

## Sample Application Guide

This section describes how to run the sample decoder application that uses the SVT-AV1 Decoder library. It describes the command line input parameters and the resulting outputs.

### Running the decoder

This section describes how to run the sample decoder application `SvtAv1DecApp.exe` (on Windows\*) or `SvtAv1DecApp` (on Linux\*) from the command line, including descriptions of the most commonly used input parameters and outputs.

The sample application typically takes the following command line parameters:

``` none
 -help                     Show usage options and exit
 -i <arg>                  Input file name
 -o <arg>                  Output file name
 -skip <arg>               Skip the first n input frames
 -limit <arg>              Stop decoding after n frames
 -bit-depth <arg>          Input bitdepth. [8, 10]
 -w <arg>                  Input picture width
 -h <arg>                  Input picture height
 -colour-space <arg>       Input picture colour space. [400, 420, 422, 444]
 -threads <arg>            Number of threads to be launched
 -parallel-frames <arg>    Number of frames to be processed in parallel
 -md5                      MD5 support flag
 -fps-frm                  Show fps after each frame decoded
 -fps-summary              Show fps summary -skip-film-grain
```

Sample usage: `SvtAv1DecApp.exe -i test.ivf -o out.yuv`

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
