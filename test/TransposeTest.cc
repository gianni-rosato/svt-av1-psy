/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbUnitTestUtility.h"
#include "EbTransforms.h"
#include "random.h"
#include "util.h"
#include "transpose_avx2.h"
#include "transpose_sse2.h"

#if !REMOVE_UNUSED_CODE
namespace {

static INLINE void transpose_8bit_16x16_c(const uint8_t *const in,
                                          uint8_t *const out) {
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            out[16 * j + i] = in[16 * i + j];
        }
    }
}

class TransposeTest : public ::testing::Test {
  public:
    ~TransposeTest();

    void SetUp() {
        rnd_ = new svt_av1_test_tool::SVTRandom(0, (1 << 16) - 1);
    }
    void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunCheckOutput();
    void RunSpeedTest();

    uint8_t in_[256];
    uint8_t out_c_[256];
    uint8_t out_o_[256];
    svt_av1_test_tool::SVTRandom *rnd_;
};

TransposeTest::~TransposeTest() {
}

void TransposeTest::RunCheckOutput() {
    __m128i in[16], out[16];

    for (int i = 0; i < 256; ++i) {
        in_[i] = rnd_->Rand8();
    }

    transpose_8bit_16x16_c(in_, out_c_);

    for (int i = 0; i < 16; i++)
        in[i] = _mm_loadu_si128((__m128i *)(in_ + 16 * i));
    transpose_8bit_16x16_reg128bit_instance_avx2(in, out);
    for (int i = 0; i < 16; i++)
        _mm_storeu_si128((__m128i *)(out_o_ + 16 * i), out[i]);
    for (int i = 0; i < 256; ++i)
        ASSERT_EQ(out_c_[i], out_o_[i]) << "[" << i << "]";
}

void TransposeTest::RunSpeedTest() {
    const int num_loops = 100000000;
    __m128i in[16], out_c[16], out_o[16];
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    for (int i = 0; i < 256; ++i) {
        in_[i] = rnd_->Rand8();
    }

    for (int i = 0; i < 16; i++)
        in[i] = _mm_loadu_si128((__m128i *)(in_ + 16 * i));

    eb_start_time(&start_time_seconds, &start_time_useconds);

    for (int i = 0; i < num_loops; ++i) {
        transpose_8bit_16x16_sse2(in, out_c);
    }

    eb_start_time(&middle_time_seconds, &middle_time_useconds);

    for (int i = 0; i < num_loops; ++i) {
        transpose_8bit_16x16_reg128bit_instance_avx2(in, out_o);
    }
    eb_start_time(&finish_time_seconds, &finish_time_useconds);

    for (int i = 0; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(out_c_ + 16 * i), out_c[i]);
        _mm_storeu_si128((__m128i *)(out_o_ + 16 * i), out_o[i]);
    }

    for (int i = 0; i < 256; ++i)
        ASSERT_EQ(out_c_[i], out_o_[i]) << "[" << i << "]";

    eb_compute_overall_elapsed_time_ms(start_time_seconds,
                                  start_time_useconds,
                                  middle_time_seconds,
                                  middle_time_useconds,
                                  &time_c);
    eb_compute_overall_elapsed_time_ms(middle_time_seconds,
                                  middle_time_useconds,
                                  finish_time_seconds,
                                  finish_time_useconds,
                                  &time_o);
    printf("Average Nanoseconds per Function Call\n");
    printf("    transpose_8bit_16x16 (SSE2) : %6.2f\n",
           1000000 * time_c / num_loops);
    printf(
        "    transpose_8bit_16x16 (AVX2) : %6.2f   (Comparison: "
        "%5.2fx)\n",
        1000000 * time_o / num_loops,
        time_c / time_o);
}

TEST_F(TransposeTest, CheckOutput) {
    RunCheckOutput();
}

TEST_F(TransposeTest, DISABLED_Speed) {
    RunSpeedTest();
}

/**
 * @brief Unit test for transpose:
 * - transpose_8bit_4x4_reg128bit_instance_sse2
 * - transpose_8bit_8x8_reg128bit_instance_sse2
 * - transpose_8bit_16x8_reg128bit_instance_sse2
 * - transpose_8bit_16x16_reg128bit_instance_sse2
 * - transpose_8bit_16x16_reg128bit_instance_avx2
 * - transpose_16bit_4x4_reg128bit_instance_sse2
 * - transpose_16bit_4x8_reg128bit_instance_sse2
 * - transpose_16bit_8x4_reg128bit_instance_sse2
 * - transpose_16bit_8x8_reg128bit_instance_sse2
 * - transpose_16bit_16x16_reg128bit_instance_sse2
 * - transpose_32bit_4x4_reg128bit_instance_sse2
 * - transpose_32bit_4x4x2_reg128bit_instance_sse2
 * - transpose_32bit_8x4_reg128bit_instance_sse2
 * - transpose_32bit_8x8_reg256bit_instance_avx2
 * - transpose_64bit_4x4_reg256bit_instance_avx2
 * - transpose_64bit_4x6_reg256bit_instance_avx2
 * - transpose_64bit_4x8_reg256bit_instance_avx2
 * - partial_transpose_8bit_8x8_reg128bit_instance_sse2
 *
 * Test strategy:
 * Since the input parameters for assembly functions are vector type,
 * the input data should be packed beforehand and output data should
 * be unpacked after the assembly function is invoked. Then unpacked
 * output should match output from c implementation.
 *
 * Expected result:
 * The result from reference and test should be the same in value
 *
 */
using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

/** reference for transpose function in c */
template <typename Sample>
static void transpose_xbit_wxh_c(uint32_t width, uint32_t height,
                                 const Sample *const in, Sample *const out) {
    for (uint32_t h = 0; h < height; h++) {
        for (uint32_t w = 0; w < width; w++) {
            out[height * w + h] = in[width * h + w];
        }
    }
}

/* intrinsic function wrapper since purely intrinsic
 * functions have no address, can not be passed as
 * function pointer.
 */

static __m256i m256_pack(const __m256i *p) {
    return _mm256_loadu_si256(p);
}

static void m256_unpack(__m256i *p, __m256i b) {
    _mm256_storeu_si256(p, b);
}

/** two function types of transpose using __m128i and __m256i */
using TransposeFunc = void (*)(const __m128i *const in, __m128i *const out);
using TransposeM256Func = void (*)(const __m256i *const in, __m256i *const out);

/** parameter types of using different function types */
using TransposeParam = std::tuple<uint32_t, uint32_t, TransposeFunc>;
using TransposeM256Param = std::tuple<uint32_t, uint32_t, TransposeM256Func>;

// maximum transpose size supported is 16x16
static const size_t mem_size = 16 * 16 * sizeof(__m256i);

/** basic template class with following types:
 * Sample: uint8_t, uint16_t, uint32_t, uint64_t
 * FuncType: TransposeFunc, TransposeM256Func
 * PackType: __m128i, __m256i
 * ParamType: TransposeParam, TransposeM256Param
 */
template <typename Sample, typename FuncType, typename PackType,
          typename ParamType>
class TransposeTestBase : public ::testing::TestWithParam<ParamType> {
  protected:
    TransposeTestBase() : rnd_(32, false), bd_(8) {
        width_ = 0;
        height_ = 0;
        func_tst_ = nullptr;
        input_tst_ = reinterpret_cast<Sample *>(eb_aom_memalign(32, mem_size));
        memset(input_tst_, 0, mem_size);
        input_ref_ = reinterpret_cast<Sample *>(eb_aom_memalign(32, mem_size));
        memset(input_ref_, 0, mem_size);
        output_tst_ = reinterpret_cast<Sample *>(eb_aom_memalign(32, mem_size));
        memset(output_tst_, 0, mem_size);
        output_tmp_ = reinterpret_cast<Sample *>(eb_aom_memalign(32, mem_size));
        memset(output_tmp_, 0, mem_size);
        output_ref_ = reinterpret_cast<Sample *>(eb_aom_memalign(32, mem_size));
        memset(output_ref_, 0, mem_size);
        pack_tool = nullptr;
        unpack_tool = nullptr;
    }

    virtual ~TransposeTestBase() {
        if (input_tst_) {
            eb_aom_free(input_tst_);
            input_tst_ = nullptr;
        }
        if (input_ref_) {
            eb_aom_free(input_ref_);
            input_ref_ = nullptr;
        }
        if (output_tst_) {
            eb_aom_free(output_tst_);
            output_tst_ = nullptr;
        }
        if (output_tmp_) {
            eb_aom_free(output_tmp_);
            output_tmp_ = nullptr;
        }
        if (output_ref_) {
            eb_aom_free(output_ref_);
            output_ref_ = nullptr;
        }
        aom_clear_system_state();
    }

    void run_check_output() {
        // prepare input data for c implementation
        size_t size = prepare_input_data();

        // prepare input data for assembly implementation
        pack_input_data();

        transpose_xbit_wxh_c(width_, height_, input_ref_, output_ref_);
        if (func_tst_) {
            func_tst_((PackType *)input_tst_, (PackType *)output_tmp_);
            unpack_output_data();
        }

        for (size_t i = 0; i < size; i++) {
            ASSERT_EQ(output_tst_[i], output_ref_[i])
                << "output mismatches at [" << i << "]";
        }
    }

    virtual size_t prepare_input_data() {
        size_t size = width_ * height_;
        for (size_t i = 0; i < size; ++i)
            input_ref_[i] = (Sample)(rnd_.random() % (1 << bd_));
        return size;
    }

    virtual void pack_input_data() {
        for (size_t i = 0; i < height_; i++)
            ((PackType *)input_tst_)[i] =
                pack_tool((PackType *)(&input_ref_[width_ * i]));
    }

    virtual void unpack_output_data() {
        for (size_t i = 0; i < width_; i++) {
            unpack_tool((PackType *)(&output_tst_[height_ * i]),
                        ((PackType *)output_tmp_)[i]);
        }
    }

  protected:
    SVTRandom rnd_;      /**< random tool */
    uint32_t bd_;        /**< bit-depth of value */
    uint32_t width_;     /**< width of test matrix */
    uint32_t height_;    /**< height of test matrix */
    FuncType func_tst_;  /**< function to test */
    Sample *input_tst_;  /**< input buffer for test */
    Sample *input_ref_;  /**< input buffer for reference */
    Sample *output_tmp_; /**< output buffer for temp */
    Sample *output_tst_; /**< output buffer for test */
    Sample *output_ref_; /**< output buffer for reference */

    typedef PackType (*PackTool)(PackType const *);
    typedef void (*UnpackTool)(PackType *, PackType);
    PackTool pack_tool;
    UnpackTool unpack_tool;
};

/**
 * @brief Unit test for 32-bit transpose with __m256i packed:
 * - transpose_32bit_8x8_reg256bit_instance_avx2
 */
class Transpose32Bit8x8Test
    : public TransposeTestBase<uint32_t, TransposeM256Func, __m256i,
                               TransposeM256Param> {
  public:
    Transpose32Bit8x8Test() {
        bd_ = 31;
        width_ = TEST_GET_PARAM(0);
        height_ = TEST_GET_PARAM(1);
        func_tst_ = TEST_GET_PARAM(2);
        pack_tool = m256_pack;
        unpack_tool = m256_unpack;
    }

    static ::testing::internal::ParamGenerator<TransposeM256Param>
    BuildParams() {
        const TransposeM256Param params[] = {
            make_tuple(8, 8, transpose_32bit_8x8_reg256bit_instance_avx2)};
        return ::testing::ValuesIn(params);
    }
};

TEST_P(Transpose32Bit8x8Test, RunCheckOutput) {
    run_check_output();
}

INSTANTIATE_TEST_CASE_P(TRANSPOSE, Transpose32Bit8x8Test,
                        Transpose32Bit8x8Test::BuildParams());

/**
 * @brief Unit test for 64-bit transpose:
 * - transpose_64bit_4x4_reg256bit_instance_avx2
 * - transpose_64bit_4x6_reg256bit_instance_avx2
 * - transpose_64bit_4x8_reg256bit_instance_avx2
 */
class Transpose64BitTest
    : public TransposeTestBase<uint64_t, TransposeM256Func, __m256i,
                               TransposeM256Param> {
  public:
    Transpose64BitTest() {
        bd_ = 31;
        width_ = TEST_GET_PARAM(0);
        height_ = TEST_GET_PARAM(1);
        func_tst_ = TEST_GET_PARAM(2);
        pack_tool = m256_pack;
        unpack_tool = m256_unpack;
    }

    size_t prepare_input_data() override {
        size_t size = width_ * height_;
        for (size_t i = 0; i < size; ++i)
            input_ref_[i] = (((uint64_t)rnd_.random()) << 32) + rnd_.random();
        return size;
    }

    void unpack_output_data() override {
        if (width_ == 4 && (height_ == 6 || height_ == 8)) {
            for (size_t i = 0; i < width_; i++) {
                unpack_tool((__m256i *)(&output_tst_[height_ * i]),
                            ((__m256i *)output_tmp_)[2 * i]);
                unpack_tool((__m256i *)(&output_tst_[height_ * i + 4]),
                            ((__m256i *)output_tmp_)[2 * i + 1]);
            }
        } else
            TransposeTestBase::unpack_output_data();
    }

    static ::testing::internal::ParamGenerator<TransposeM256Param>
    BuildParams() {
        const TransposeM256Param params[] = {
            make_tuple(4, 4, transpose_64bit_4x4_reg256bit_instance_avx2),
            make_tuple(4, 6, transpose_64bit_4x6_reg256bit_instance_avx2),
            make_tuple(4, 8, transpose_64bit_4x8_reg256bit_instance_avx2)};
        return ::testing::ValuesIn(params);
    }
};

TEST_P(Transpose64BitTest, RunCheckOutput) {
    run_check_output();
}

INSTANTIATE_TEST_CASE_P(TRANSPOSE, Transpose64BitTest,
                        Transpose64BitTest::BuildParams());
};  // namespace
#endif
