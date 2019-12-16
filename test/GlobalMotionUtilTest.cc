/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file GlobalMotionUtilTest.cc
 *
 * @brief Unit test for utility functions in global motion:
 * - ransac_affine
 * - ransac_rotzoom
 * - ransac_translation
 * - ransac_affine_double_prec
 * - ransac_rotzoom_double_prec
 * - ransac_translation_double_prec
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
extern "C" {
#include "ransac.h"
}
#include "random.h"
#include "util.h"

using std::tuple;
using std::vector;
using svt_av1_test_tool::SVTRandom;
namespace {

using Point = struct {
    double x;
    double y;
    bool is_inlier;
};
using AffineMat = struct {
    /* clang-format off */
    /* affine matrix is defined as
     * | m0  m1  m2 |
     * | m3  m4  m5 |
     * |  0   0   1 |
     */
    /* clang-format on */
    double m0;
    double m1;
    double m2;
    double m3;
    double m4;
    double m5;
};

static const int CoordinateMax = (1 << 16) - 1;
static const int PointCountMin = 15; /**< 3*MINPTS_MULTIPLIER */

/**
 * @brief Unit test for RANSAC functions:
 * - ransac_affine
 * - ransac_rotzoom
 * - ransac_translation
 * - ransac_affine_double_prec
 * - ransac_rotzoom_double_prec
 * - ransac_translation_double_prec
 *
 * Test strategy:
 * Create a pair of 2D point sets by the matrix of affine transform
 * (translation, zoom, rotate and affine), then add some noise in test data;
 * check the motion parameters in the result of RANSAC function
 *
 * Expected result:
 * The difference between affine transform matrix and the parameters in motion
 * should be less than threshold
 *
 */
template <typename Sample, typename FuncType>
class RansacTest : public ::testing::TestWithParam<TransformationType> {
  protected:
    RansacTest() : rnd_(0, CoordinateMax) {
        data_.clear();
        ref_.clear();
        memset(&mat_, 0, sizeof(mat_));
    }

    void generate_data(int count, int &inliers) {
        /** ransac function requires more than 15 samples */
        ASSERT_GE(count, PointCountMin);
        inliers = rnd_.random() % count;
        inliers = inliers < 10 ? 10 : inliers;
        for (int i = 0; i < inliers; i++)
            data_.push_back(
                {(double)rnd_.random(), (double)rnd_.random(), true});

        TransDataFunc trans_func = get_trans_data_func(GetParam());
        ASSERT_NE(trans_func, nullptr);
        const int data_count = trans_func(rnd_, data_, ref_, mat_);
        ASSERT_EQ(data_count, inliers);

        /** add noise, less than 25% */
        const int noise_count = AOMMAX(AOMMIN(count - inliers, inliers >> 2),
                                       PointCountMin - inliers);
        for (int i = 0; i < noise_count; i++) {
            size_t insert_pos = rnd_.random() % data_.size();
            data_.insert(data_.begin() + insert_pos,
                         {(double)rnd_.random(), (double)rnd_.random(), false});
            ref_.insert(ref_.begin() + insert_pos,
                        {(double)rnd_.random(), (double)rnd_.random(), false});
        }
    }

    void run_test(size_t times) {
        for (size_t i = 0; i < times; i++) {
            data_.clear();
            ref_.clear();

            const int max_data_count = MAX_CORNERS;
            int inliers = 0;
            generate_data(max_data_count, inliers);
            ASSERT_NE(inliers, 0) << "generate data failed!";
            do_ransac_check();
        }
    }

    virtual void prepare_input(Sample *input, int npoints) = 0;

    virtual FuncType get_ransac_func(TransformationType type) = 0;

    void do_ransac_check() {
        const int npoints = (int)data_.size();
        Sample *points = new Sample[npoints * 4];
        ASSERT_NE(points, nullptr);
        prepare_input(points, npoints);

        const int num_motions = RANSAC_NUM_MOTIONS;
        int *num_inliers_by_motion = new int[num_motions];
        ASSERT_NE(num_inliers_by_motion, nullptr);
        memset(num_inliers_by_motion, 0, sizeof(*num_inliers_by_motion));

        MotionModel *motions = new MotionModel[num_motions];
        ASSERT_NE(motions, nullptr);
        for (int i = 0; i < num_motions; i++) {
            memset(&motions[i], 0, sizeof(MotionModel));
            motions[i].inliers = new int[2 * MAX_CORNERS];
            ASSERT_NE(motions[i].inliers, nullptr);
            memset(motions[i].inliers, 0, sizeof(int) * 2 * npoints);
        }

        FuncType ransac_func = (FuncType)get_ransac_func(GetParam());
        ASSERT_NE(ransac_func, nullptr);
        int ret = ransac_func(
            points, npoints, num_inliers_by_motion, motions, num_motions);
        ASSERT_EQ(ret, 0);

        /** check for the number of inlier */
        ASSERT_NE(num_inliers_by_motion[0], 0);

        /** check for the transform matrix of motion */
        check_transform_matrix(mat_, motions[0].params);

        if (points)
            delete[] points;
        if (num_inliers_by_motion)
            delete[] num_inliers_by_motion;
        if (motions) {
            for (int i = 0; i < num_motions; i++) {
                if (motions[i].inliers)
                    delete[] motions[i].inliers;
            }
            delete[] motions;
        }
    }

    /* clang-format off */
    /* the transform matrix in motion parameters is defined:
     * | m2  m3  m0 |
     * | m4  m5  m1 |
     * | m6  m7   1 |
     */
    /* clang-format on */
    void check_transform_matrix(const AffineMat &mat, const double *params) {
        ASSERT_NEAR(mat.m0, params[2], 0.0001f);
        ASSERT_NEAR(mat.m1, params[3], 0.0001f);
        ASSERT_NEAR(mat.m2, params[0], 1.0f);
        ASSERT_NEAR(mat.m3, params[4], 0.0001f);
        ASSERT_NEAR(mat.m4, params[5], 0.0001f);
        ASSERT_NEAR(mat.m5, params[1], 1.0f);
    }

    using TransDataFunc = int (*)(SVTRandom &, vector<Point> &, vector<Point> &,
                                  AffineMat &);
    TransDataFunc get_trans_data_func(TransformationType type) {
        switch (type) {
        case TRANSLATION: return transform_data_translation;
        case ROTZOOM: return transform_data_zoom_rotate;
        case AFFINE: return transform_data_affine;
        default: return nullptr;
        }
    }

    static int transform_data_translation(SVTRandom &rnd, vector<Point> &data,
                                          vector<Point> &ref, AffineMat &mat) {
        ref.clear();
        const int offset_x =
            (rnd.random() >> 2) * ((rnd.random() % 2) ? -1 : 1);
        const int offset_y =
            (rnd.random() >> 2) * ((rnd.random() % 2) ? -1 : 1);
        for (vector<Point>::iterator it = data.begin(); it != data.end();
             ++it) {
            Point point = translate(*it, offset_x, offset_y, mat);
            ref.push_back(point);
        }
        return (int)ref.size();
    }

    static int transform_data_zoom_rotate(SVTRandom &rnd, vector<Point> &data,
                                          vector<Point> &ref, AffineMat &mat) {
        ref.clear();
        /** limit zoom rate from 50% to 150% */
        const double zoom_scale =
            (double)(rnd.random() % 1000) / 1000.0f + 0.5f;
        const double theta = PI * (rnd.random() % 360) / 360;
        AffineMat mat_zoom, mat_rotate;
        memset(&mat_zoom, 0, sizeof(mat_zoom));
        memset(&mat_rotate, 0, sizeof(mat_rotate));
        for (vector<Point>::iterator it = data.begin(); it != data.end();
             ++it) {
            Point point = zoom(*it, zoom_scale, mat_zoom);
            point = rotate(point, theta, mat_rotate);
            ref.push_back(point);
        }
        mat.m0 = mat_rotate.m0 * mat_zoom.m0;
        mat.m1 = mat_rotate.m1 * mat_zoom.m0;
        mat.m2 = 0;
        mat.m3 = mat_rotate.m3 * mat_zoom.m4;
        mat.m4 = mat_rotate.m4 * mat_zoom.m4;
        mat.m5 = 0;
        return (int)ref.size();
    }

    static int transform_data_affine(SVTRandom &rnd, vector<Point> &data,
                                     vector<Point> &ref, AffineMat &mat) {
        ref.clear();
        mat = {(double)(rnd.random() % 500) / 1000.0f *
                   (((rnd.random() % 2) ? -1 : 1)),
               (double)(rnd.random() % 500) / 1000.0f *
                   (((rnd.random() % 2) ? -1 : 1)),
               (double)(rnd.random() >> 2) * (((rnd.random() % 2) ? -1 : 1)),
               (double)(rnd.random() % 500) / 1000.0f *
                   (((rnd.random() % 2) ? -1 : 1)),
               (double)(rnd.random() % 500) / 1000.0f *
                   (((rnd.random() % 2) ? -1 : 1)),
               (double)(rnd.random() >> 2) * (((rnd.random() % 2) ? -1 : 1))};
        for (vector<Point>::iterator it = data.begin(); it != data.end();
             ++it) {
            Point point = affine(*it, mat);
            ref.push_back(point);
        }
        return (int)ref.size();
    }

    /* clang-format off */
    /* translation transform
     * | x |   | 1  0  offset_x |   | x+offset_x |
     * | y | X | 0  1  offset_y | = | y+offset_y |
     * | 1 |   | 0  0     1     |   |     1      |
     */
    /* clang-format on */
    static Point translate(const Point &pt, double offset_x, double offset_y,
                           AffineMat &mat) {
        mat = {1.0f, 0.0f, offset_x, 0.0f, 1.0f, offset_y};
        return affine(pt, mat);
    }

    /* clang-format off */
    /* zoom transform
     * | x |   | scale  0    0 |   | x*scale |
     * | y | X |   0  scale  0 | = | y*scale |
     * | 1 |   |   0    0    1 |   |    1    |
     */
    /* clang-format on */
    static Point zoom(const Point &pt, double scale, AffineMat &mat) {
        mat = {scale, 0.0f, 0.0f, 0.0f, scale, 0.0f};
        return affine(pt, mat);
    }

    /* clang-format off */
    /* rotate transform
     * | x |   | cos(t)  -sin(t)  0 |   | x*cos(t) - y*sin(t) |
     * | y | X | sin(t)  cos(t)   0 | = | x*sin(t) + y*cos(t) |
     * | 1 |   |   0       0      1 |   |          1          |
     */
    /* clang-format on */
    static Point rotate(const Point &pt, double theta, AffineMat &mat) {
        double sin_theta = sin(theta);
        double cos_theta = cos(theta);
        mat = {cos_theta, -1 * sin_theta, 0.0f, sin_theta, cos_theta, 0.0f};
        return affine(pt, mat);
    }

    /* clang-format off */
    /* affine transform
     * | x |   | m0  m1  m2 |   | x*m0 + y*m1 + m2 |
     * | y | X | m3  m4  m5 | = | x*m3 + y*m4 + m5 |
     * | 1 |   |  0   0   1 |   |         1        |
     */
    /* clang-format on */
    static Point affine(const Point &pt, AffineMat mat) {
        Point dst;
        dst.x = (pt.x * mat.m0) + (pt.y * mat.m1) + mat.m2;
        dst.y = (pt.x * mat.m3) + (pt.y * mat.m4) + mat.m5;
        return dst;
    }

  protected:
    SVTRandom rnd_;
    vector<Point> data_;
    vector<Point> ref_;
    AffineMat mat_;
};

class RansacIntTest : public RansacTest<int, RansacFunc> {
  protected:
    void prepare_input(int *input, int npoints) {
        for (int i = 0; i < npoints; i++) {
            input[4 * i] = (int)round(data_.at(i).x);
            input[4 * i + 1] = (int)round(data_.at(i).y);
            input[4 * i + 2] = (int)round(ref_.at(i).x);
            input[4 * i + 3] = (int)round(ref_.at(i).y);
        }
    }

    RansacFunc get_ransac_func(TransformationType type) {
        return av1_get_ransac_type(type);
    }
};

class RansacDoubleTest : public RansacTest<double, RansacFuncDouble> {
  protected:
    void prepare_input(double *input, int npoints) {
        for (int i = 0; i < npoints; i++) {
            input[4 * i] = data_.at(i).x;
            input[4 * i + 1] = data_.at(i).y;
            input[4 * i + 2] = ref_.at(i).x;
            input[4 * i + 3] = ref_.at(i).y;
        }
    }

    RansacFuncDouble get_ransac_func(TransformationType type) {
        return av1_get_ransac_double_prec_type(type);
    }
};

static const TransformationType transform_table[] = {
    TRANSLATION, ROTZOOM, AFFINE};

TEST_P(RansacIntTest, CheckOutput) {
    run_test(1000);
};

INSTANTIATE_TEST_CASE_P(GlobalMotion, RansacIntTest,
                        ::testing::ValuesIn(transform_table));

TEST_P(RansacDoubleTest, CheckOutput) {
    run_test(1000);
};

INSTANTIATE_TEST_CASE_P(GlobalMotion, RansacDoubleTest,
                        ::testing::ValuesIn(transform_table));

}  // namespace
