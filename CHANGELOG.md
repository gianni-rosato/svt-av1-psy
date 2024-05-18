# Changelog

## [2.1.0] - 2024-05-17

API updates
- One config parameter added within the padding size. Config param structure size remains unchanged
- Presets 6 and 12 are now pointing to presets 7 and 13 respectively due to the lack of spacing between the presets
- Further preset shuffling is being discussed in #2152

Encoder
- Added variance boost support to improve visual quality for the tune vq mode
- Improve the tradeoffs for the random access mode across presets:
-   Speedup of 12-40% presets M0, M3, M5 and M6 while maintaining similar quality levels
-   Improved the compression efficiency of presets M11-M13 by 1-2% (!2213)
- Added ARM optimizations for functions with c_only equivalent

Cleanup Build and bug fixes and documentation
- Use nasm as a default assembler and yasm as a fallback
- Fix performance regression for systems with multiple processor groups
- Enable building SvtAv1ApiTests and SvtAv1E2ETests for arm
- Added variance boost documentation
- Added a mailmap file to map duplicate git generated emails to the appropriate author

## [2.0.0] - 2024-03-13

Major API updates
- Changed the API signaling the End Of Stream (EOS) with the last frame vs with an empty frame
- OPT_LD_LATENCY2 making the change above is kept in the code to help devs with integration
- Removed the 3-pass VBR mode which changed the calling mechanism of multi-pass VBR
- Moved to a new versioning scheme where the project major version will be updated everytime API/ABI is changed

Encoder
- Improve the tradeoffs for the random access mode across presets:
-   Speedup presets MR by ~100% and improved quality along with tradeoff improvements across the higher quality presets (!2179)
-   Improved the compression efficiency of presets M9-M13 by 1-4% (!2179)
-   Simplified VBR multi-pass to use 2 passes to allow integration with ffmpeg
- Continued adding ARM optimizations for functions with c_only equivalent
- Replaced the 3-pass VBR with a 2-pass VBR to ease the multi-pass integration with ffmpeg
- Memory savings of 20-35% for LP 8 mode in preset M6 and below and 1-5% in other modes / presets

Cleanup and bug fixes and documentation
- Various cleanups and functional bug fixes
- Update the documentation to reflect the rate control changes

## [1.8.0] - 2023-12-11

Encoder
- Improve the tradeoffs for the random access mode across presets:
- Speedup CRF presets M6 to M0 by 17-53% while maintaining similar quality levels
- Re-adjust CRF presets M7 to M13 for better quality with BD-rate gains ranging from 1-4%
- Improve the quality and speed of the 1-pass VBR mode
- More details on the per preset improvements can be found in MR !2143
- Add API allowing to update bitrate / CRF and Key_frame placement during the encoding session for CBR lowdelay mode and CRF Random Access mode
- ARM Neon SIMD optimizations for most critical kernels allowing for a 4.5-8x fps speedup vs the c implementation

Cleanup and bug fixes and documentation
- Various cleanups and functional bug fixes
- Update the documentation for preset options and individual features

## [1.7.0] - 2023-08-23

Encoder
- Improve the tradeoffs for the random access mode across presets MR-M13:
 - Quality improvements across all presets and metrics ranging from 0.3% to 4.5% in BD-rate (!2129)
 - Spacing between presets [M1-M6] has been adjusted to account for the tradeoff improvements achieved
  - As a user guidance when comparing v1.7 vs v1.6 in a convexhull encoding setup:
   - v1.7.0 M2 is now at similar quality levels as v1.6.0 M1 while being ~50% faster
   - v1.7.0 M3 is now at similar quality levels as v1.6.0 M2 while being ~50% faster
   - v1.7.0 M4 is now at similar quality levels as v1.6.0 M3 while being ~40% faster
   - v1.7.0 M5 is now at similar quality levels as v1.6.0 M4 while being ~30% faster
   - v1.7.0 M6 is now at similar quality levels as v1.6.0 M5 while being ~25% faster
- Added an experimental tune SSIM mode yielding ~3-4% additional SSIM BD-rate gains (!2109)

Build, cleanup and bug fixes
- Various cleanups and functional bug fixes
- Fix build conflict with libaom

## [1.6.0] - 2023-06-18

Encoder
- Improve the tradeoffs for the random access mode across presets M1-M13:
 - Speeding up the higher quality presets by 30-40%
 - Improving the BD-rate by 1-2% for the faster presets
- Improve the tradeoffs for the low delay mode for both scrren content and non-screen content encoding modes
- Add a toggle to remove the legacy one-frame buffer at the input of the pipeline alowing the low delay mode to operate at sub-frame processing latencies
- Add a new API allowing the user to specify quantization offsets for a region of interest per frame

Build, cleanup and bug fixes
- Various cleanups and functional bug fixes
- Fix the startup minigop size BD-rate loss
- Add ability to run the ci-testing offline

## [1.5.0] - 2023-04-25

Encoder
- Optimize the tradeoffs for M0-M13 speeding up M5-M1 by 15-30% and improving the BDR of M6-M13 by 1-3%
- Create a new preset MR (--preset -1) to be the quality reference
- Optimize the tradeoffs for M8-M13 in the low delay encoding mode (!2052 !2096 and !2102) for SC and non-SC modes
- Add dynamic minigop support for the random access configuration enabled by default in M9 and below
- Add support to allow users to specify lambda scaling factors through the commandline
- Rewrite the gstreamer plugin and updating it to be uptodate with the latest API changes
- Add skip frames feature allowing the user to start encoding after n frames in the file
- Add ability to specify a smaller startup minigop size for every gop to enable faster prefetching
- Fix segmentation support and re-enable it with --aq-mode 1 to allow work on the region of interest API
- Add padding bytes to the EbSvtAv1EncConfiguration configuration structure keep its size unchanged until v2.0

Build, Cleanup and Documentation
- Major cleanups for unused variables, static functions, and comments formatting
- Reduce the size of variables
- Refine app level parsing and reference scaling API calls in the application
- Add dynamic minigop documentation along with updating the documentation accordingly

## [1.4.1] - 2022-12-12

Bugfixes:
- Fix CRF with maxrate bug causing bitrate to be significantly limited for CRF encodings
- Fix command line parsing forcing 1-pass in a 2-pass encoding mode when the --keyint=`x`s format is used
- Fix decoder segfault due to assuming aligned buffers in the inverse transform assembly

## [1.4.0] - 2022-11-30

Encoder
- Adopted the 6L / 32-picture mini-GOP configuraion as default and adjusted MD feature levels accordingly yielding higher compression efficiency gains
- Update the TPL model to account for the larger mini-gop size
- Re-tune presets M0-M12 using the gained coding efficiency for improved quality vs cycles tradeoffs
- Allow duplicate commandline parameters to be parsed and take into consideration the latter to allow AWCY integration

Build, Cleanup and Documentation
- Make include and lib paths friendly to abs and rel paths
- Update preset and API documentation
- Various functional bug fixes
- Remove manual prediction structure support

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
