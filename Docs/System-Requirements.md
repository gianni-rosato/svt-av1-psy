[Top level](../README.md)

# System Requirements

## Operating System

SVT-AV1 Encoder may run on any Windows* or Linux* 64 bit operating systems. The
list below represents the operating systems that the encoder application and
library were tested and validated on:

- __Windows* Operating Systems (64-bit):__
  - Windows* Server 2016
- __Linux* Operating Systems (64-bit):__
  - Ubuntu* 16.04 Server LTS
  - Ubuntu* 18.04 Server LTS
- __Unix* Operating Systems (64-bit):__
  - MacOS

## Hardware

The SVT-AV1 Encoder library supports the x86 architecture

- __CPU Requirements__

    In order to achieve the performance targeted by the SVT-AV1 Encoder, the
    specific CPU model listed above would need to be used when running the
    encoder. Otherwise, the encoder runs on any 5th Generation Intel® Core™
    processor, (Intel® Xeon® CPUs, E5-v4 or newer).

- __RAM Requirements__

    The SVT-AV1 Encoder adapts to the system on which it is being run. The
    memory requirements depend on the number of cores the system contains, the
    input frame rate of the input sequence (`--fps`) and the look ahead
    distance passed to the encoder (`--lad`). The SVT-AV1 Encoder application
    will display an error if the system does not have enough RAM to support the
    encode prior to the start of the encode.



