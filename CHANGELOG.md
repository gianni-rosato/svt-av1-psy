# Changelog

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
