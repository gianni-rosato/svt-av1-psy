# Changelog

## [1.3.0] - 2022-10-18

Encoder
- Port SIMD optimizations from libDav1D making the conformant path (Inv. Transform) faster
- Enabling smaller mini-GOP size configurations and tuning it for the low delay mode
- Tuning the low-latency mode in random access targeting latencies from 250ms to 1s
- Adding GOP-constrained Rate Control targeting low-latency streaming applications
- Optimize mode decision features levels for depth partitioning, RDOQ, MD stage0 pruning in-loop filtering temporal filtering and TPL adding more granularity and gaining further quality
- Preset tuning M0-M13 to smooth the spacing and utilize the quality improvements towards better tradeoffs

Build, Cleanup and Documentation
- Update preset and API documentation
- Various functional bug fixes
- Remove the use of GLOB in cmake and use file names

## [1.2.1] - 2022-08-15

- Fix a crash at the end of the encode that may occur when an invalid metadata packet is sent with the EOS packet

Build, Cleanup
- y4m header pasring code cleanup
- API cleanup and enhancements adding string options for RC mode
- Added option to build without app / dec / enc using the build.sh / build.bat scripts

## [1.2.0] - 2022-08-02

Encoder
- Improve CRF preset tradeoffs for both the default and fast-decode modes
- Improve the SSIM-based tradeoffs of the presets without impacting those of PSNR / VMAF
- Improve CBR mode by enhancing the bit-distribution within the gop
- Added support for reference frame scaling
- Added support for quantization matrices
- Added svtparams patches applicable to ffmpeg 4.4
- AVX2 optimizations for low-delay mode
- TPL-based VBR mode improvements
- Improved Chroma RDOQ
- Improve TPL QP Scaling
- Add length info to ivf header
- Fix support for metadata pass-through
- Add ability to specify Chroma and Luma qindex offsets independently on top of CRF qp assignments

Build, Cleanup and Documentation
- Fix multiple API documentation mismatches
- Updated features documentation
- Various functional bug fixes

## [1.1.0] - 2022-05-17

Encoder
- TPL tradeoff optimizations for 4L pred structure
- Quality-vs-cycles tradeoff improvements across all presets
- Add ability to force key_frame positions through ffmpeg for CRF mode
- Minimize the quality impact of fast-decode while maintaining the decoder speedup
- AVX2 optimizations for low delay mode
- Fix VQ issues #1896 #1857 and #1819

Build, Cleanup and Documentation
- API / ABI cleanup and implement independent versioning
- Add UEB_DLL for static linking with pkgconf
- Update system requirements docs
- Rate control code refactoring
- Fix AVX512 vs AVX2 mismatch

## [1.0.0] - 2022-04-22

Encoder
- Added S-frames support
- CBR Rate control mode for low delay
- Added support for chroma position signalling
- Added support for skipping denoising pictures after film grain synthesis
- Extend fast-decode support to cover presets M0-M10
- Simplified --fast-decode to have only one level
- Optimized --fast-decode level 1 for better tradeoffs
- Visual quality improvements addressing issues #1819 / #1297
- Visual quality fixes and improvements for both tune 0 and 1
- Quality vs density tradeoffs tuning across all presets in CRF mode with TPL improvements
- Update default settings to use a longer gop / higher quality preset and lower CRF value
- Various code cleanups and memory optimizations
- Additional AVX2 optimizations
- Fixed all known functional bugs
- More robust rate control parameter verification

Build and Documentation
- Major documentation update and re-structure
- Added more user guides, preset guides and common questions section
- Improve CI coverage
- Reduced unnecessary warnings
- Improved the documentation of the configuration parameters
- Improve Unit Test Coverage
- Address C vs asm mismatches

## [0.9.1] - 2022-02-23

Encoder
- New `key=val` API for setting encoder options along with `--svtav1-params` appside
- New `--fast-decode` option for producing bitstreams tuned for faster decoding M5- M10
- New `--tune` support for a subjectively optimized encoding mode
- Quality vs density tradeoffs improvements for 4k resolutions and reduction or resolution checks
- New SSE kernels for better encoding speed on older hardware
- Bugfix: Removed `DISABLE_REALTIME` and instead implemented better checks for realtime encoding
- All open library bugs resolved and closed

Build and Testing
- Windows: Moved .dll files to the binary directory next to the exe to fix dll loading issues
- Windows: General VS 2017 compilation speedup
- Windows: Added VS 2022
- BSD: Fixed compilation issues and errors surounding GLIBC specific interfaces and conflicting names
- Dockerfile added

## [0.9] - 2022-01-19

Encoder
- New faster presets M9-M12, M12 reaching similar complexity level to that of x264 veryfast
- New multi-pass and single pass VBR implementation minimizing the quality difference vs CRF while reducing the cycle overhead
- Quality vs density tradeoffs improvements across all presets in CRF mode
- Added support for CRF with capped bitrate
- Added support for overlay frames and super resolution
- Fixed film grain synthesis bugs
- Added experimental support for > 4k resolutions
- Added experimental support for the low delay prediction structure
- Significant memory reduction especially for faster presets in a multi-threaded environment
- API configuration structure cleanup removing all invalid or out of date parameters
- Speedup legacy CPUs for faster development by adding SSE code for corresponding C-kernels
- Updated the code license from BSD 2-clause to BSD 3-clause clear
- Cleaned up the code for various kernels
- Updated the user guide and feature documentation

Build and Testing
- Bug fixes
- Improve CI coverage
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Fix static analysis warnings / errors


## [0.8.7] - 2021-05-08

Encoder
- Feature optimizations: creating new mode decision / encode-decode feature levels allowing better speed / quality trade-off granularity
- Preset repositioning after adopting new levels
- Preset 8 achieving similar speed levels as x265 medium tuned for VOD operating mode
- New 1-pass and 2-pass VBR implementation ported from libaom and adapted to the SVT architecture - still a WIP
- Cleaned up old VBR and CVBR RC code along with the lookahead mechanism associated with them
- 2-pass encoding is only available for VBR mode and 1-pass with lookahead is used for CRF
- Improvements for TPL algorithm to handle long clips and easy content
- Memory optimizations, cleaning up data structures to reduce memory usage up to 2x memory reduction in multi-threaded VBR environment
- Additional AVX2 and AVX512 optimizations
- Cleaned up unused command line parameters and left the config params that are linked to ffmpeg / gst
- Update documentation accordingly
- Added HDR support and color primaries signalling (off by default until integrated with ffmpeg)

Build and Testing
- Bug fixes
- Improve CI coverage
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Fix static analysis warnings / errors

## [0.8.6] - 2020-11-28

Encoder
- Further quality-speed tradeoffs tuning for VOD use cases
- Improved TPL support within 1-pass and 2-pass CRF moode
- Continued non-optimized support for 2pass VBR and CRF
- Align kernel nomenclature to prefix svt_aom for kernels brough from libaom to avoid symbol conflicts

Build and Testing
- Bug fixes
- Improve CI
- Added CI support for gitlab
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Fix static analysis warnings / errors
- Add address sanitizer
- Fix symbol conflicts with libaom and libvpx when staticly linked to ffmpeg

## [0.8.5] - 2020-09-04
Relicensing notice
- Change the outbound license from BSD+Patent to the AOM license / patent

Encoder
- Added tpl support to adaptively change lambda and quantization parameters within the frame
- Added multi staged hme support
- Quality speed trade-offs tuned to VOD use cases
- Added first level non-optimized support for 2pass VBR and CRF
- Added combined cli two pass support with options for stats being written to a memory buffer and a specified file
- Added non square partitioning optimizations
- Improved lambda generation

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Fix static analysis warnings / errors
- Add address sanitizer
- Fix symbol conflicts with libaom and libvpx when staticly linked to ffmpeg


## [0.8.4] - 2020-06-26

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Fix static analysis warnings / errors
- Add address sanitizer
- Various ffmpeg patch fixes

## [0.8.3] - 2020-04-24

Encoder
- Presets optimization

Build and Testing
- Bug fixes
- Xcode build support


## [0.8.2] - 2020-04-18

Encoder
- Initial Super resolution support
- Mode decision rate estimation enhancements
- Altref improvements
- Manual prediction structure support
- Enhanced RDOQ support
- Improved warp motion
- Improved documentation and help menu
- New command line parameters
- Fix non-multiple of 8 resolution video corruption for 8bit
- Improved multi-stage mode decision support
- Added multi-stage motion estimation support
- Chroma mode decision optimizations
- Added 16bit pipeline support
- Added Mode decision non-square partition weights
- Enhanced qp-scaling assignment
- Memory optimizations
- Added AVX512 Optimizations
- Added AVX2 Optimizations

Decoder
- Improved multi-threading support
- Memory optimizations
- Updated documentation
- Loop filter flow optimization

Encoder and Decoder
- Encoder-Decoder-Common folder structure improvements

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches
- Support C only builds (for platforms other than x86)

## [0.8.1] - 2020-01-28

Encoder
- Palette support for 10-bit
- Added the first batch of detailed documentation to help developers

Encoder and Decoder
- Code cleanup and refactoring

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches


## [0.8.0] - 2019-12-20

Encoder
- Preset Optimizations
- Single-core execution memory optimization [-lp 1 -lad 0]
- Rate estimation update enhancements
- Added on / off flags for feature run-time switching
- Added auto-max partitioning algorithm
- Multi-pass partitioning depth support
- Remove deprecated RC mode 1 and shifter RC mode 2 and mode 3 to mode 1 and mode 2 respectively
- Update Cost Calculation for CDEF Filtering
- Intra-Inter Compound for 10-bit
- Eigth-pel optimization
- Added AVX512 Optimizations
- Added AVX2 Optimizations

Decoder
- Initial multi-threading support
- Decoder optimizations / cleanup

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches

## [0.7.5] - 2019-11-24

Encoder
- RDOQ for 10-bit
- Inter Intra Class pruning at MD-Staging
- Global Motion Vector support for 8-bit and 10-bit
- Interpolation Filter Search support for 10-bit
- Palette Prediction support
- 2-pass encoding support
- ATB 10-bit support at the encode pass
- Simplified MD Staging [only 3 stages]
- Inter-Inter and Inter-Intra Compound for 10-bit
- Intra Paeth for 10-bit
- Filter Intra Prediction
- New-Near and Near-New support
- OBMC Support for 8-bit and 10-bit
- RDOQ Chroma
- ATB Support for Inter Blocks
- Temporal Filtering for 10-bit
- Eight-pel support in predictive ME
- MCTS Tiles support
- Added AVX512 Optimizations
- Added AVX2 Optimizations

Decoder
- SuperRes support
- Reference Frame Scaling support
- 12-bit support
- Annex B support

Build and Testing
- Bug fixes
- Improve CI
- Improve Unit Test Coverage
- Address C vs asm mismatches

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
