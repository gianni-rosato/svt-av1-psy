# Scalable Video Technology for AV1 Decoder (SVT-AV1 Decoder) User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Sample Application Guide](#sample-application-guide)
    - [Running the decoder](#running-the-decoder)

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
