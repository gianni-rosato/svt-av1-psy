# Changelog

## [upcoming]

Decoder
- Decoder SuperRes support
- Decoder Reference Frame Scaling support

## [0.7.0] - 2019-09-26

Encoder
- Enhanced MRP Reference Frames
- Intra Inter Compound
- QP Modulation support
- MFMV Support
- MD Staging design [Up to 4 MD stages and 3 prediction classes: Intra / Inter / Compound]
- Compound Motion prediction
- 10-bit Mode Decision support for Intra
- Thread safe resource allocation
- Added AVX512 Optimizations
- Added AVX2 Optimizations

Decoder
- Screen Content Tools
- Temporal MV scan support
- Inter support
- Screen Content Tools support
- Post Processing Filters support
- Compound Mode (InterInter & InterIntra) Tool support
- Decoder Film Grain support

Build and Testing
- Improve CI
- Improve build scripts
- Improve cmake lists
- Improve Unit Test Coverage
- API update
- Bug fixes

## [0.6.0] - 2019-06-28

- Inital decoder implementation
- Static library support
- Alt-ref pictures - temporal filtering
- Adaptive Transform Block for INTRA
- Adaptive QP scaling
- Decoder - Support for Tiles and 10bit
- API - add option to calculate / report PSNR values
- Support for segmentation
- SIMD Optimizations
- Downsampling 2x2 filtering
- Handle incomplete SBs
- MACROS / trailing / tabs-spaces cleanup

## [0.5.0] - 2019-05-17

- 8 bit / 10 bit 4:2:0 up to 4k60p resolutions
- Presets 0-8
- New API, FFmpeg, GStreamer plugins
- Rate control support  (VBR, CVBR)
- Block sizes from 4x4 to 128x128
- Non-square blocks
- Tiles
- Deblocking / CDEF / Restoration filters
- Film Grain
- Warped motion estimation
- Intra block copy
- Trellis quantized coefficient optimization
- Support for 4 and 5 layers prediction structures
- Chroma search in MD
- Multi-reference picture support
