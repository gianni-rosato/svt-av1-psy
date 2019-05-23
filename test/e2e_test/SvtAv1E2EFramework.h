/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2EFramework.h
 *
 * @brief Defines a test framework for End to End test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#ifndef _SVT_AV1_E2E_FRAMEWORK_H_
#define _SVT_AV1_E2E_FRAMEWORK_H_

#include "E2eTestVectors.h"
#include "FrameQueue.h"
#include "PerformanceCollect.h"
#include "CompareTools.h"

class RefDecoder;
extern RefDecoder *create_reference_decoder();

#define INPUT_SIZE_576p_TH 0x90000    // 0.58 Million
#define INPUT_SIZE_1080i_TH 0xB71B0   // 0.75 Million
#define INPUT_SIZE_1080p_TH 0x1AB3F0  // 1.75 Million
#define INPUT_SIZE_4K_TH 0x29F630     // 2.75 Million
#define EB_OUTPUTSTREAMBUFFERSIZE_MACRO(resolution_size) \
    ((resolution_size) < (INPUT_SIZE_1080i_TH)           \
         ? 0x1E8480                                      \
         : (resolution_size) < (INPUT_SIZE_1080p_TH)     \
               ? 0x2DC6C0                                \
               : (resolution_size) < (INPUT_SIZE_4K_TH) ? 0x2DC6C0 : 0x2DC6C0)

// Copied from EbAppProcessCmd.c
#define LONG_ENCODE_FRAME_ENCODE 4000
#define SPEED_MEASUREMENT_INTERVAL 2000
#define START_STEADY_STATE 1000
#define AV1_FOURCC 0x31305641  // used for ivf header
#define IVF_STREAM_HEADER_SIZE 32
#define IVF_FRAME_HEADER_SIZE 12
#define OBU_FRAME_HEADER_SIZE 3
#define TD_SIZE 2

/** @defgroup svt_av1_e2e_test Test framework for E2E test
 *  Defines the framework body of E2E test for the mainly test progress
 *  @{
 */
namespace svt_av1_e2e_test {

using namespace svt_av1_e2e_test_vector;
using namespace svt_av1_e2e_tools;
using namespace svt_av1_video_source;

/** SvtAv1Context is a set of test contexts in whole test progress */
typedef struct {
    EbComponentType
        *enc_handle; /**< encoder handle, created from encoder library */
    EbSvtAv1EncConfiguration enc_params; /**< encoder parameter set */
    EbBufferHeaderType
        *output_stream_buffer; /**< output buffer of encoder in test */
    EbBufferHeaderType
        *input_picture_buffer; /**< input buffer of encoder in test */
} SvtAv1Context;

/** SvtAv1E2ETestFramework is a class with impelmention of video source control,
 * encoding progress, decoding progress, data collection and data comparision */
class SvtAv1E2ETestFramework
    : public ::testing::TestWithParam<TestVideoVector> {
  public:
    typedef struct IvfFile {
        FILE *file;
        uint64_t byte_count_since_ivf;
        uint64_t ivf_count;
        IvfFile(std::string path);
        ~IvfFile() {
            if (file) {
                fclose(file);
                file = nullptr;
            }
            byte_count_since_ivf = 0;
            ivf_count = 0;
        }
    } IvfFile;

  protected:
    SvtAv1E2ETestFramework();
    virtual ~SvtAv1E2ETestFramework();

  protected:
    void SetUp() override;
    void TearDown() override;
    /** initialization for test */
    virtual void init_test();
    /** close for test */
    virtual void close_test();
    /** test processing body */
    virtual void run_encode_process();

  public:
    static VideoSource *prepare_video_src(const TestVideoVector &vector);
    static void trans_src_param(const VideoSource *source,
                                EbSvtAv1EncConfiguration &config);
    /** get reconstructed frame from encoder, it should call after send data
     * @param ctxt  context of encoder
     * @param recon  video frame queue of reconstructed
     * @param is_eos  flag of recon frames is eos
     * into decoder */
    static void get_recon_frame(const SvtAv1Context &ctxt, FrameQueue *recon,
                                bool &is_eos);

  private:
    /** write ivf header to output file */
    void write_output_header();
    /** write compressed data into file
     * @param output  compressed data from encoder
     */
    void write_compress_data(const EbBufferHeaderType *output);
    /** process compressed data by write to file for send to decoder
     * @param data  compressed data from encoder
     */
    void process_compress_data(const EbBufferHeaderType *data);
    /** send compressed data to decoder
     * @param data  compressed data from encoder, single OBU
     * @param size  size of compressed data
     */
    void decode_compress_data(const uint8_t *data, const uint32_t size);
    /** check video frame psnr with source
     * @param frame  video frame from reference decoder
     */
    void check_psnr(const VideoFrame &frame);

  protected:
    VideoSource *video_src_;   /**< video source context */
    SvtAv1Context av1enc_ctx_; /**< AV1 encoder context */
    uint32_t start_pos_;       /**< start position of video frame */
    uint32_t frames_to_test_;  /**< frame count for this test */
    FrameQueue *recon_queue_;  /**< reconstructed frame collection */
    RefDecoder *refer_dec_;    /**< reference decoder context */
    IvfFile *output_file_;     /**< file handle for save encoder output data */
    uint8_t obu_frame_header_size_; /**< size of obu frame header */
    PerformanceCollect *collect_;   /**< performance and time collection*/
    VideoSource *psnr_src_;         /**< video source context for psnr */
    ICompareQueue *ref_compare_; /**< sink of reference to compare with recon*/
    PsnrStatistics pnsr_statistics_; /**< psnr statistics recorder.*/
    uint64_t total_enc_out_;
};

}  // namespace svt_av1_e2e_test
/** @} */  // end of svt_av1_e2e_test_vector

#endif  //_SVT_AV1_E2E_FRAMEWORK_H_
