/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file PaletteModeUtilTest.cc
 *
 * @brief Unit test for util functions in palette mode:
 * - eb_av1_count_colors
 * - av1_count_colors_highbd
 * - av1_k_means_dim1
 * - av1_k_means_dim2
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/
#include <vector>
#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "random.h"
#include "util.h"

using std::tuple;
using std::vector;
using svt_av1_test_tool::SVTRandom;
namespace {

extern "C" int eb_av1_count_colors(const uint8_t *src, int stride, int rows,
                                   int cols, int *val_count);
extern "C" int av1_count_colors_highbd(uint16_t *src, int stride, int rows,
                                       int cols, int bit_depth, int *val_count);

/**
 * @brief Unit test for counting colors:
 * - eb_av1_count_colors
 * - av1_count_colors_highbd
 *
 * Test strategy:
 * Feeds the random value both into test function and the vector without
 * duplicated, then compares the count of result and the individual item count
 * in vector.
 *
 * Expected result:
 * The count numbers from test function and vector are the same.
 *
 * Test coverage:
 * The input can be 8-bit and 8-bit/10-bit/12-bit for HBD cases
 */
template <typename Sample>
class ColorCountTest : public ::testing::Test {
  protected:
    ColorCountTest() : rnd_(16, false) {
        input_ =
            (Sample *)eb_aom_memalign(32, MAX_PALETTE_SQUARE * sizeof(Sample));
        bd_ = 8;
        ref_.clear();
        val_count_ = nullptr;
    }

    ~ColorCountTest() {
        if (input_) {
            eb_aom_free(input_);
            input_ = nullptr;
        }
        aom_clear_system_state();
    }

    void prepare_data() {
        memset(input_, 0, MAX_PALETTE_SQUARE * sizeof(Sample));
        ref_.clear();
        const int32_t mask = (1 << bd_) - 1;
        for (size_t i = 0; i < MAX_PALETTE_SQUARE; i++) {
            input_[i] = rnd_.random() & mask;
            /** put same value into a vector for reference */
            ref_.push_back(input_[i]);
        }
        /** remove all duplicated items */
        std::sort(ref_.begin(), ref_.end());
        vector<int>::iterator it = std::unique(ref_.begin(), ref_.end());
        ref_.erase(it, ref_.end());
    }

    void run_test(size_t times) {
        const int max_colors = (1 << bd_);
        val_count_ = (int *)eb_aom_memalign(32, max_colors * sizeof(int));
        for (size_t i = 0; i < times; i++) {
            prepare_data();
            ASSERT_EQ(count_color(), ref_.size())
                << "color count failed at: " << i;
        }
        if (val_count_) {
            eb_aom_free(val_count_);
            val_count_ = nullptr;
        }
    }

    virtual unsigned int count_color() = 0;

  protected:
    SVTRandom rnd_;
    Sample *input_;
    uint8_t bd_;
    vector<int> ref_;
    int *val_count_;
};

class ColorCountLbdTest : public ColorCountTest<uint8_t> {
  protected:
    unsigned int count_color() override {
        const int max_colors = (1 << bd_);
        memset(val_count_, 0, max_colors * sizeof(int));
        unsigned int colors =
            (unsigned int)eb_av1_count_colors(input_, 64, 64, 64, val_count_);
        return colors;
    }
};

TEST_F(ColorCountLbdTest, MatchTest) {
    run_test(1000);
}

class ColorCountHbdTest : public ColorCountTest<uint16_t> {
  protected:
    unsigned int count_color() override {
        const int max_colors = (1 << bd_);
        memset(val_count_, 0, max_colors * sizeof(int));
        unsigned int colors = (unsigned int)av1_count_colors_highbd(
            input_, 64, 64, 64, bd_, val_count_);
        return colors;
    }
};

TEST_F(ColorCountHbdTest, MatchTest8Bit) {
    bd_ = 8;
    run_test(1000);
}

TEST_F(ColorCountHbdTest, MatchTest10Bit) {
    bd_ = 10;
    run_test(1000);
}

TEST_F(ColorCountHbdTest, MatchTest12Bit) {
    bd_ = 12;
    run_test(1000);
}

extern "C" void av1_k_means_dim1(const int *data, int *centroids,
                                 uint8_t *indices, int n, int k, int max_itr);
extern "C" void av1_k_means_dim2(const int *data, int *centroids,
                                 uint8_t *indices, int n, int k, int max_itr);
static const int MaxItr = 50;
/**
 * @brief Unit test for kmeans functions:
 * - av1_k_means_dim1
 * - av1_k_means_dim2
 *
 * Test strategy:
 * Feeds the plane buffer with random colors into kmeans function and get the
 * centroids and indices, verifies each color being the closest to the centroid
 * in all candidates.
 *
 * Expected result:
 * Every pixels are closest to their centroid in all candidates
 *
 * Test coverage:
 * Tests for K from PALETTE_MIN_SIZE to PALETTE_MAX_SIZE
 */
class KMeansTest : public ::testing::TestWithParam<int> {
  protected:
    KMeansTest() : rnd_(8, false), palette_rnd_(2, 64) {
        k_ = GetParam();
        data_ = new int[2 * MAX_PALETTE_SQUARE];
    }

    ~KMeansTest() {
        if (data_) {
            delete[] data_;
            data_ = nullptr;
        }
    }

    /** functions for 1d test */
    int prepare_data(const int max_colors) {
        memset(data_, 0, MAX_PALETTE_SQUARE * sizeof(int));
        uint8_t *palette = new uint8_t[max_colors];
        for (int i = 0; i < max_colors; i++)
            palette[i] = rnd_.random();
        uint8_t tmp[MAX_PALETTE_SQUARE] = {0};
        for (size_t i = 0; i < MAX_PALETTE_SQUARE; i++)
            data_[i] = tmp[i] = palette[rnd_.random() % max_colors];
        delete[] palette;
        int val_count[MAX_PALETTE_SQUARE] = {0};
        return eb_av1_count_colors(tmp, 64, 64, 64, val_count);
    }

    void run_test(size_t times) {
        uint8_t indices[MAX_PALETTE_SQUARE] = {0};
        for (size_t i = 0; i < times; i++) {
            const int max_colors = palette_rnd_.random();
            const int colors = prepare_data(max_colors);
            int centroids[PALETTE_MAX_SIZE] = {0};
            int k = AOMMIN(colors, k_);
            av1_k_means_dim1(
                data_, centroids, indices, MAX_PALETTE_SQUARE, k, MaxItr);
            check_output(centroids, k, data_, indices, MAX_PALETTE_SQUARE);
        }
    }

    void check_output(const int *centroids, const int k, const int *data,
                      const uint8_t *indices, const int n) {
        for (int i = 0; i < n; i++) {
            const int min_delta = abs(data[i] - centroids[indices[i]]);
            for (int j = 0; j < k; j++) {
                const int delta = abs(data[i] - centroids[j]);
                ASSERT_GE(delta, min_delta)
                    << "index error at " << i << ", value is " << data[i]
                    << ", distance to centroid( " << centroids[indices[i]]
                    << ") is greater than to " << centroids[j];
            }
        }
    }

    /** functions for 2d test */
    int prepare_data_2d(const int max_colors) {
        memset(data_, 0, 2 * MAX_PALETTE_SQUARE * sizeof(int));
        uint16_t *palette = new uint16_t[max_colors];
        for (int i = 0; i < max_colors; i++)
            palette[i] = (rnd_.random() << 8) + rnd_.random();
        vector<uint16_t> val_vec;
        for (size_t i = 0; i < MAX_PALETTE_SQUARE; i++) {
            uint16_t tmp = palette[rnd_.random() % max_colors];
            data_[2 * i] = tmp >> 8;
            data_[2 * i + 1] = tmp & 0xFF;
            val_vec.push_back(tmp);
        }
        delete[] palette;
        std::sort(val_vec.begin(), val_vec.end());
        vector<uint16_t>::iterator it =
            std::unique(val_vec.begin(), val_vec.end());
        val_vec.erase(it, val_vec.end());
        return (const int)val_vec.size();
    }

    void run_test_2d(size_t times) {
        uint8_t indices[2 * MAX_PALETTE_SQUARE] = {0};
        for (size_t i = 0; i < times; i++) {
            const int max_colors = palette_rnd_.random();
            const int colors = prepare_data_2d(max_colors);
            int centroids[2 * PALETTE_MAX_SIZE] = {0};
            int k = AOMMIN(colors, k_);
            av1_k_means_dim2(
                data_, centroids, indices, MAX_PALETTE_SQUARE, k, MaxItr);
            check_output_2d(centroids, k, data_, indices, MAX_PALETTE_SQUARE);
        }
    }

    static double distance_2d(int x1, int y1, int x2, int y2) {
        int x_d = x1 - x2;
        int y_d = y1 - y2;
        return sqrt(x_d * x_d + y_d * y_d);
    }

    void check_output_2d(const int *centroids, const int k, const int *data,
                         const uint8_t *indices, const int n) {
        for (int i = 0; i < n; i++) {
            const double min_delta = distance_2d(data[2 * i],
                                                 data[2 * i + 1],
                                                 centroids[2 * indices[i]],
                                                 centroids[2 * indices[i] + 1]);
            for (int j = 0; j < k; j++) {
                const double delta = distance_2d(data[2 * i],
                                                 data[2 * i + 1],
                                                 centroids[2 * j],
                                                 centroids[2 * j + 1]);
                ASSERT_GE(delta, min_delta)
                    << "index error at " << i << ", value is " << data[i]
                    << ", distance to centroid( " << centroids[indices[i]]
                    << ") is greater than to " << centroids[j];
            }
        }
    }

  protected:
    int *data_;
    int k_;
    SVTRandom rnd_;
    SVTRandom palette_rnd_;
};

TEST_P(KMeansTest, CheckOutput) {
    run_test(1000);
};

TEST_P(KMeansTest, CheckOutput2D) {
    run_test_2d(1000);
};

INSTANTIATE_TEST_CASE_P(PalleteMode, KMeansTest,
                        ::testing::Range(PALETTE_MIN_SIZE, PALETTE_MAX_SIZE));

}  // namespace
