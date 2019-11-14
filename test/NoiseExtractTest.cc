/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file NoiseExtractTest.cc
 *
 * @brief Unit test for noise extract functions
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "gtest/gtest.h"
#include <stdlib.h>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "random.h"
extern "C" {
#include "EbPictureControlSet.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureAnalysisProcess.h"
}

void eb_picture_buffer_desc_dctor(EbPtr p) {
    EbPictureBufferDesc *obj = (EbPictureBufferDesc *)p;
    if (obj->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_FREE_ALIGNED_ARRAY(obj->buffer_y);
        EB_FREE_ALIGNED_ARRAY(obj->buffer_bit_inc_y);
    }
    if (obj->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_FREE_ALIGNED_ARRAY(obj->buffer_cb);
        EB_FREE_ALIGNED_ARRAY(obj->buffer_bit_inc_cb);
    }
    if (obj->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_FREE_ALIGNED_ARRAY(obj->buffer_cr);
        EB_FREE_ALIGNED_ARRAY(obj->buffer_bit_inc_cr);
    }
}

using svt_av1_test_tool::SVTRandom;
class ExtractFilterTest : public ::testing::Test {
  public:
    ExtractFilterTest() {
        // set the height and width not multiple of 64 to test the incomplete sb
        width_ = 336;
        height_ = 240;

        EbPictureBufferDescInitData pbd_init_data;
        pbd_init_data.max_width = width_;
        pbd_init_data.max_height = height_;
        pbd_init_data.bit_depth = EB_8BIT;
        // allocate all the components
        pbd_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        pbd_init_data.left_padding = 0;
        pbd_init_data.right_padding = 0;
        pbd_init_data.top_padding = 0;
        pbd_init_data.bot_padding = 0;
        pbd_init_data.color_format = EB_YUV420;
        pbd_init_data.split_mode = EB_FALSE;

        subsampling_x_ = (pbd_init_data.color_format == EB_YUV444 ? 1 : 2) - 1;
        subsampling_y_ = (pbd_init_data.color_format >= EB_YUV422 ? 1 : 2) - 1;

        eb_picture_buffer_desc_ctor(&in_pic_, &pbd_init_data);
        eb_picture_buffer_desc_ctor(&denoised_pic_tst_, &pbd_init_data);
        eb_picture_buffer_desc_ctor(&denoised_pic_ref_, &pbd_init_data);
        eb_picture_buffer_desc_ctor(&noise_pic_tst_, &pbd_init_data);
        eb_picture_buffer_desc_ctor(&noise_pic_ref_, &pbd_init_data);
    }

    ~ExtractFilterTest() {
        eb_picture_buffer_desc_dctor(&in_pic_);
        eb_picture_buffer_desc_dctor(&denoised_pic_tst_);
        eb_picture_buffer_desc_dctor(&denoised_pic_ref_);
        eb_picture_buffer_desc_dctor(&noise_pic_tst_);
        eb_picture_buffer_desc_dctor(&noise_pic_ref_);
    }

    virtual void init_pic(EbPictureBufferDesc *pic) {
        SVTRandom rnd(0, 255);
        // Y data
        uint8_t *buf = pic->buffer_y;
        int32_t stride = pic->stride_y;
        for (uint32_t i = 0; i < height_; ++i)
            for (uint32_t j = 0; j < width_; ++j)
                buf[i * stride + j] = rnd.random();

        buf = pic->buffer_cb;
        stride >>= subsampling_x_;
        for (uint32_t i = 0; i < (height_ >> subsampling_y_); ++i)
            for (uint32_t j = 0; j < (width_ >> subsampling_x_); ++j)
                buf[i * stride + j] = rnd.random();

        buf = pic->buffer_cr;
        for (uint32_t i = 0; i < (height_ >> subsampling_y_); ++i)
            for (uint32_t j = 0; j < (width_ >> subsampling_x_); ++j)
                buf[i * stride + j] = rnd.random();
    }

    EbBool check_pic_content(EbPictureBufferDesc *pic1,
                             EbPictureBufferDesc *pic2) {
        int mismatch_samples_y = 0;
        int mismatch_samples_cb = 0, mismatch_samples_cr = 0;

        // check y data
        uint8_t *buf1 = pic1->buffer_y;
        uint8_t *buf2 = pic2->buffer_y;
        int32_t stride = pic1->stride_y;
        for (uint32_t i = 0; i < height_; ++i) {
            for (uint32_t j = 0; j < width_; ++j) {
                if (buf1[i * stride + j] != buf2[i * stride + j])
                    ++mismatch_samples_y;
            }
        }

        // check cb data
        buf1 = pic1->buffer_cb;
        buf2 = pic2->buffer_cb;
        stride >>= subsampling_x_;
        for (uint32_t i = 0; i < (height_ >> subsampling_y_); ++i) {
            for (uint32_t j = 0; j < (width_ >> subsampling_x_); ++j) {
                if (buf1[i * stride + j] != buf2[i * stride + j])
                    ++mismatch_samples_cb;
            }
        }

        // check cr data
        buf1 = pic1->buffer_cr;
        buf2 = pic2->buffer_cr;
        for (uint32_t i = 0; i < (height_ >> subsampling_y_); ++i) {
            for (uint32_t j = 0; j < (width_ >> subsampling_x_); ++j) {
                if (buf1[i * stride + j] != buf2[i * stride + j])
                    ++mismatch_samples_cr;
            }
        }

        EXPECT_EQ(0, mismatch_samples_y)
            << mismatch_samples_y << " y samples are different.";
        EXPECT_EQ(0, mismatch_samples_cb)
            << mismatch_samples_cb << " cb samples are different.";
        EXPECT_EQ(0, mismatch_samples_cr)
            << mismatch_samples_cr << " cr samples are different.";

        return HasFailure() ? EB_FALSE : EB_TRUE;
    }

    void run_luma_test() {
        init_pic(&in_pic_);

        for (uint32_t sb_y = 0; sb_y < height_; sb_y += 64) {
            for (uint32_t sb_x = 0; sb_x < width_; sb_x += 64) {
                noise_extract_luma_weak_c(
                    &in_pic_, &denoised_pic_ref_, &noise_pic_ref_, sb_y, sb_x);
                noise_extract_luma_weak_avx2_intrin(
                    &in_pic_, &denoised_pic_tst_, &noise_pic_tst_, sb_y, sb_x);

                EbBool ret =
                    check_pic_content(&denoised_pic_ref_, &denoised_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_weak test: Denoised pic different "
                       "with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";

                ret = check_pic_content(&noise_pic_ref_, &noise_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_weak test: noise pic different with "
                       "dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";
            }
        }
    }

    void run_luma_lcu_test() {
        init_pic(&in_pic_);

        for (uint32_t sb_y = 0; sb_y < height_; sb_y += 64) {
            for (uint32_t sb_x = 0; sb_x < width_; sb_x += 64) {
                noise_extract_luma_weak_lcu_c(
                    &in_pic_, &denoised_pic_ref_, &noise_pic_ref_, 0, 0);
                noise_extract_luma_weak_lcu_avx2_intrin(
                    &in_pic_, &denoised_pic_tst_, &noise_pic_tst_, 0, 0);

                EbBool ret =
                    check_pic_content(&denoised_pic_ref_, &denoised_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_weak_lcu test: Denoised pic "
                       "different with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";

                ret = check_pic_content(&noise_pic_ref_, &noise_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_weak_lcu test: noise pic different "
                       "with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";
            }
        }
    }

    void run_luma_strong_test() {
        init_pic(&in_pic_);

        for (uint32_t sb_y = 0; sb_y < height_; sb_y += 64) {
            for (uint32_t sb_x = 0; sb_x < width_; sb_x += 64) {
                noise_extract_luma_strong_c(
                    &in_pic_, &denoised_pic_ref_, sb_y, sb_x);
                noise_extract_luma_strong_avx2_intrin(
                    &in_pic_, &denoised_pic_tst_, sb_y, sb_x);

                EbBool ret =
                    check_pic_content(&denoised_pic_ref_, &denoised_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_strong test: Denoised pic different "
                       "with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";

                ret = check_pic_content(&noise_pic_ref_, &noise_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_luma_strong test: noise pic different "
                       "with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";
            }
        }
    }

    void run_chroma_test() {
        init_pic(&in_pic_);

        for (uint32_t sb_y = 0; sb_y < height_; sb_y += 64) {
            for (uint32_t sb_x = 0; sb_x < width_; sb_x += 64) {
                noise_extract_chroma_weak_c(&in_pic_, &denoised_pic_ref_,
                    sb_y >> subsampling_y_, sb_x >> subsampling_x_);
                noise_extract_chroma_weak_avx2_intrin(&in_pic_, &denoised_pic_tst_,
                    sb_y >> subsampling_y_, sb_x >> subsampling_x_);

                EbBool ret =
                    check_pic_content(&denoised_pic_ref_, &denoised_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "Denoised pic different with dim [" << width_ << " x "
                    << height_ << "] sb_origin [" << sb_x << " x " << sb_y
                    << "]";

                ret = check_pic_content(&noise_pic_ref_, &noise_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_chroma_weak noise pic different with dim "
                       "["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";
            }
        }
    }

    void run_chroma_strong_test() {
        init_pic(&in_pic_);

        for (uint32_t sb_y = 0; sb_y < height_; sb_y += 64) {
            for (uint32_t sb_x = 0; sb_x < width_; sb_x += 64) {
                noise_extract_chroma_strong_c(&in_pic_, &denoised_pic_ref_,
                    sb_y >> subsampling_y_, sb_x >> subsampling_x_);
                noise_extract_chroma_strong_avx2_intrin(&in_pic_, &denoised_pic_tst_,
                    sb_y >> subsampling_y_, sb_x >> subsampling_x_);

                EbBool ret =
                    check_pic_content(&denoised_pic_ref_, &denoised_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "Denoised pic different with dim [" << width_ << " x "
                    << height_ << "] sb_origin [" << sb_x << " x " << sb_y
                    << "]";

                ret = check_pic_content(&noise_pic_ref_, &noise_pic_tst_);
                EXPECT_EQ(ret, EB_TRUE)
                    << "noise_extract_chroma_strong test, noise pic different "
                       "with dim ["
                    << width_ << " x " << height_ << "] sb_origin [" << sb_x
                    << " x " << sb_y << "]";
            }
        }
    }

  protected:
    uint32_t width_;
    uint32_t height_;

    int subsampling_x_;
    int subsampling_y_;
    EbPictureBufferDesc in_pic_;
    EbPictureBufferDesc denoised_pic_tst_;
    EbPictureBufferDesc denoised_pic_ref_;
    EbPictureBufferDesc noise_pic_tst_;
    EbPictureBufferDesc noise_pic_ref_;
};

TEST_F(ExtractFilterTest, LumaWeakNoiseTest) {
    run_luma_test();
}

TEST_F(ExtractFilterTest, LumaWeakNoiseLcuTest) {
    run_luma_lcu_test();
}

TEST_F(ExtractFilterTest, LumaStrongNoiseTest) {
    run_luma_strong_test();
}

TEST_F(ExtractFilterTest, ChromaWeakNoiseTest) {
    run_chroma_test();
}

TEST_F(ExtractFilterTest, ChromaStrongNoiseTest) {
    run_chroma_strong_test();
}
