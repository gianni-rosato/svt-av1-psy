# Scalable Video Technology for AV1 Encoder (SVT-AV1 Encoder)

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/OpenVisualCloud/SVT-AV1?branch=master&svg=true)](https://ci.appveyor.com/project/OpenVisualCloud/SVT-AV1)
[![Travis Build Status](https://travis-ci.org/OpenVisualCloud/SVT-AV1.svg?branch=master)](https://travis-ci.org/OpenVisualCloud/SVT-AV1)
[![Coverage Status](https://coveralls.io/repos/github/OpenVisualCloud/SVT-AV1/badge.svg?branch=master)](https://coveralls.io/github/OpenVisualCloud/SVT-AV1?branch=master)

The Scalable Video Technology for AV1 Encoder (SVT-AV1 Encoder) is an AV1-compliant encoder library core. The SVT-AV1 development is a work-in-progress targeting performance levels applicable to both VOD and Live encoding / transcoding video applications.

# License

SVT-AV1 Encoder is licensed under the OSI-approved BSD+Patent license. See [LICENSE](LICENSE.md) for details.

# Documentation

More details about the SVT-AV1 Encoder usage can be found under:
-   [svt-av1-encoder-user-guide](Docs/svt-av1_encoder_user_guide.md)

# System Requirements

## Operating System

SVT-AV1 Encoder may run on any Windows* or Linux* 64 bit operating systems. The list below represents the operating systems that the encoder application and library were tested and validated on:

* __Windows* Operating Systems (64-bit):__

    -  Windows* Server 2016

* __Linux* Operating Systems (64-bit):__

    -  Ubuntu* 16.04 Server LTS

    -  Ubuntu* 18.04 Server LTS

* __Unix* Operating Systems (64-bit):__

    -  MacOS

## Hardware

The SVT-AV1 Encoder library supports the x86 architecture

* __CPU Requirements__

In order to achieve the performance targeted by the SVT-AV1 Encoder, the specific CPU model listed above would need to be used when running the encoder. Otherwise, the encoder runs on any 5th Generation Intel® Core™ processor, (Intel® Xeon® CPUs, E5-v4 or newer).

* __RAM Requirements__
The SVT-AV1 Encoder adapts to the system that is being ran on. The memory requirements depend on the number of cores the system contains, the input frame rate of the input sequence (-fps) and the look ahead distance passed to the encoder (-lad). The SVT-AV1 Encoder application will display an error if the system does not have enough RAM to support the encode prior to the start of the encode. The following table shows the minimum amount of RAM required for some standard resolutions of 10bit video per stream:


|       Resolution      | 8-vCPU Commit Size (GB)| 40-vCPU Commit Size (GB)|
|-----------------------|------------------------|-------------------------|
|       4k              |           14           |           24            |
|       1080p           |            6           |           10            |
|       720p            |            4           |            7            |
|       480p            |            3           |            5            |


# Build and Install

## Windows* Operating Systems (64-bit):

* __Build Requirements__
    -    Visual Studio* 2017 (can be downloaded [here](https://www.visualstudio.com/vs/older-downloads/))
    -    CMake 3.5 or later (can be downloaded [here](https://github.com/Kitware/CMake/releases/download/v3.13.0/cmake-3.13.0-win64-x64.msi))
    -   YASM Assembler version 1.2.0 or later
    -    Download the yasm exe from the following [link](http://www.tortall.net/projects/yasm/releases/yasm-1.3.0-win64.exe)
    -    Rename yasm-1.3.0-win64.exe to yasm.exe
    -   Copy yasm.exe into a location that is in the PATH environment variable

* __Build Instructions__
    -    Generate the Visual Studio* 2017 project files by following the steps below cd Build\windows
        -    run generate_vs17.bat [such would generate the visual studio project files]
    -    Open the "svt-av1.sln" using Visual Studio* 2017 and click on Build -- > Build Solution

* __Binaries Location__
    -   Binaries can be found under <repo dir>\Bin/Release or <repo dir>\Bin/Debug, depending on whether Debug or Release were selected in the build mode

* __Installation__
-    For the binaries to operate properly on your system, the following conditions have to be met:
    -    On any of the Windows* Operating Systems listed in the OS requirements section, install Visual Studio* 2017
    -    Once the installation is complete, copy the binaries to a location making sure that both the sample application "SvtAv1EncApp.exe” and library "SvtAv1Enc.dll” are in the same folder.
    -    Open the command prompt window at the chosen location and run the sample application to encode. SvtAV1EncApp.exe -i [in.yuv] -w [width] -h [height] -b [out.ivf].
    -    Sample application supports reading from pipe. E.g. ffmpeg -i [input.mp4] -nostdin -f rawvideo -pix_fmt yuv420p - | SvtAv1EncApp.exe -i stdin -n [number_of_frames_to_encode] -w [width] -h [height].

## Linux* Operating Systems (64-bit):

* __Build Requirements__
     -    GCC 5.4.0 or later
     -    CMake 3.5.1 or later
     -    YASM Assembler version 1.2.0 or later

* __Build Instructions__
	 -	./Build/linux/build.sh <release | debug> (if none specified, both release and debug will be built)


* __Sample Binaries location__
     -    Binaries can be found under Bin/Release and / or Bin/Debug

* __Installation__
For the binaries to operate properly on your system, the following conditions have to be met:
    -    On any of the Linux* Operating Systems listed above, copy the binaries under a location of your choice.
    -    Change the permissions on the sample application “SvtAV1EncApp” executable by running the command:                 chmod +x SvtAv1EncApp
    -    cd into your chosen location
    -    Run the sample application to encode. ./SvtAv1EncApp -i [in.yuv] -w [width] -h [height] -b [out.ivf].
    -    Sample application supports reading from pipe. E.g. ffmpeg -i [input.mp4] -nostdin -f rawvideo -pix_fmt yuv420p - | ./SvtAv1EncApp -i stdin -n [number_of_frames_to_encode] -w [width] -h [height].

# Demo features and limitations

-  **Multi-instance support:** The multi-instance functionality is a demo feature implemented in the SVT-AV1 Encoder sample application as an example of one sample application using multiple encoding libraries. Encoding using the multi-instance support is limited to only 6 simultaneous streams. For example two channels encoding on Windows: SvtAV1EncApp.exe -nch 2 -c firstchannel.cfg secondchannel.cfg
-  **Features enabled:** The library will display an error message any feature combination that is not currently supported. 

# How to Contribute

We welcome community contributions to the SVT-AV1 Encoder. Thank you for your time! By contributing to the project, you agree to the license and copyright terms in the OSI-approved BSD+Patent license and to the release of your contribution under these terms. See [LICENSE](LICENSE.md) for details.

## Contribution process

-  Follow the [coding guidelines](STYLE.md)

-  Validate that your changes do not break a build

-  Perform smoke tests and ensure they pass

-  Submit a pull request for review to the maintainer

# How to Report Bugs and Provide Feedback

Use the "Issues" tab on Github. To avoid duplicate issues, please make sure you go through the existing issues before logging a new one.

# Notices and Disclaimers

The notices and disclaimers can be found [here](NOTICES.md)
