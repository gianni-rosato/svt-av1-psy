/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>
#include <stdlib.h>
#include "EbDefinitions.h"
#include "EbModeDecisionProcess.h"
#include "aom_dsp_rtcd.h"

#define DIVIDE_AND_ROUND(x, y) (((x) + ((y) >> 1)) / (y))

// Generate a random number in the range [0, 32768).
static INLINE unsigned int lcg_rand16(unsigned int *state) {
    *state = (unsigned int)(*state * 1103515245ULL + 12345);
    return *state / 65536 % 32768;
}

#define AV1_K_MEANS_RENAME(func, dim) func##_dim##dim##_c

void AV1_K_MEANS_RENAME(av1_calc_indices, 1)(const int *data, const int *centroids,
                                             uint8_t *indices, int n, int k);
void AV1_K_MEANS_RENAME(av1_calc_indices, 2)(const int *data, const int *centroids,
                                             uint8_t *indices, int n, int k);
void AV1_K_MEANS_RENAME(av1_k_means, 1)(const int *data, int *centroids, uint8_t *indices, int n,
                                        int k, int max_itr);
void AV1_K_MEANS_RENAME(av1_k_means, 2)(const int *data, int *centroids, uint8_t *indices, int n,
                                        int k, int max_itr);

// Given 'n' 'data' points and 'k' 'centroids' each of dimension 'dim',
// calculate the centroid 'indices' for the data points.
static inline void av1_calc_indices(const int *data, const int *centroids, uint8_t *indices, int n,
                                    int k, int dim) {
    if (dim == 1) {
        av1_calc_indices_dim1(data, centroids, indices, n, k);
    } else if (dim == 2) {
        av1_calc_indices_dim2(data, centroids, indices, n, k);
    } else {
        assert(0 && "Untemplated k means dimension");
    }
}

// Given 'n' 'data' points and an initial guess of 'k' 'centroids' each of
// dimension 'dim', runs up to 'max_itr' iterations of k-means algorithm to get
// updated 'centroids' and the centroid 'indices' for elements in 'data'.
// Note: the output centroids are rounded off to nearest integers.
static inline void av1_k_means(const int *data, int *centroids, uint8_t *indices, int n, int k,
                               int dim, int max_itr) {
    if (dim == 1) {
        av1_k_means_dim1(data, centroids, indices, n, k, max_itr);
    } else if (dim == 2) {
        av1_k_means_dim2(data, centroids, indices, n, k, max_itr);
    } else {
        assert(0 && "Untemplated k means dimension");
    }
}

#define AV1_K_MEANS_DIM 1
#include "k_means_template.h"
#undef AV1_K_MEANS_DIM
#define AV1_K_MEANS_DIM 2
#include "k_means_template.h"
#undef AV1_K_MEANS_DIM

static int int_comparer(const void *a, const void *b) { return (*(int *)a - *(int *)b); }

int av1_remove_duplicates(int *centroids, int num_centroids) {
    int num_unique; // number of unique centroids
    int i;
    qsort(centroids, num_centroids, sizeof(*centroids), int_comparer);
    // Remove duplicates.
    num_unique = 1;
    for (i = 1; i < num_centroids; ++i) {
        if (centroids[i] != centroids[i - 1]) { // found a new unique centroid
            centroids[num_unique++] = centroids[i];
        }
    }
    return num_unique;
}

static int delta_encode_cost(const int *colors, int num, int bit_depth, int min_val) {
    if (num <= 0) return 0;
    int bits_cost = bit_depth;
    if (num == 1) return bits_cost;
    bits_cost += 2;
    int       max_delta = 0;
    int       deltas[PALETTE_MAX_SIZE];
    const int min_bits = bit_depth - 3;
    for (int i = 1; i < num; ++i) {
        const int delta = colors[i] - colors[i - 1];
        deltas[i - 1]   = delta;
        assert(delta >= min_val);
        if (delta > max_delta) max_delta = delta;
    }
    int bits_per_delta = AOMMAX(av1_ceil_log2(max_delta + 1 - min_val), min_bits);
    assert(bits_per_delta <= bit_depth);
    int range = (1 << bit_depth) - colors[0] - min_val;
    for (int i = 0; i < num - 1; ++i) {
        bits_cost += bits_per_delta;
        range -= deltas[i];
        bits_per_delta = AOMMIN(bits_per_delta, av1_ceil_log2(range));
    }
    return bits_cost;
}

int av1_index_color_cache(const uint16_t *color_cache, int n_cache, const uint16_t *colors,
                          int n_colors, uint8_t *cache_color_found, int *out_cache_colors) {
    if (n_cache <= 0) {
        for (int i = 0; i < n_colors; ++i) out_cache_colors[i] = colors[i];
        return n_colors;
    }
    memset(cache_color_found, 0, n_cache * sizeof(*cache_color_found));
    int n_in_cache = 0;
    int in_cache_flags[PALETTE_MAX_SIZE];
    memset(in_cache_flags, 0, sizeof(in_cache_flags));
    for (int i = 0; i < n_cache && n_in_cache < n_colors; ++i) {
        for (int j = 0; j < n_colors; ++j) {
            if (colors[j] == color_cache[i]) {
                in_cache_flags[j]    = 1;
                cache_color_found[i] = 1;
                ++n_in_cache;
                break;
            }
        }
    }
    int j = 0;
    for (int i = 0; i < n_colors; ++i)
        if (!in_cache_flags[i]) out_cache_colors[j++] = colors[i];
    assert(j == n_colors - n_in_cache);
    return j;
}

int av1_palette_color_cost_y(const PaletteModeInfo *const pmi, uint16_t *color_cache, int n_cache,
                             int bit_depth) {
    const int n = pmi->palette_size[0];
    int       out_cache_colors[PALETTE_MAX_SIZE];
    uint8_t   cache_color_found[2 * PALETTE_MAX_SIZE];
    const int n_out_cache = av1_index_color_cache(
        color_cache, n_cache, pmi->palette_colors, n, cache_color_found, out_cache_colors);
    const int total_bits = n_cache + delta_encode_cost(out_cache_colors, n_out_cache, bit_depth, 1);
    return av1_cost_literal(total_bits);
}

static void palette_add_to_cache(uint16_t *cache, int *n, uint16_t val) {
    // Do not add an already existing value
    if (*n > 0 && val == cache[*n - 1]) return;

    cache[(*n)++] = val;
}

int eb_get_palette_cache(const MacroBlockD *const xd, int plane, uint16_t *cache) {
    const int row = -xd->mb_to_top_edge >> 3;
    // Do not refer to above SB row when on SB boundary.
    const MbModeInfo *const above_mi = (row % (1 << MIN_SB_SIZE_LOG2)) ? xd->above_mbmi : NULL;
    const MbModeInfo *const left_mi  = xd->left_mbmi;
    int                     above_n = 0, left_n = 0;
    if (above_mi) above_n = above_mi->palette_mode_info.palette_size[plane != 0];
    if (left_mi) left_n = left_mi->palette_mode_info.palette_size[plane != 0];
    if (above_n == 0 && left_n == 0) return 0;
    int             above_idx    = plane * PALETTE_MAX_SIZE;
    int             left_idx     = plane * PALETTE_MAX_SIZE;
    int             n            = 0;
    const uint16_t *above_colors = above_mi ? above_mi->palette_mode_info.palette_colors : NULL;
    const uint16_t *left_colors  = left_mi ? left_mi->palette_mode_info.palette_colors : NULL;
    // Merge the sorted lists of base colors from above and left to get
    // combined sorted color cache.
    while (above_n > 0 && left_n > 0) {
        uint16_t v_above = above_colors[above_idx];
        uint16_t v_left  = left_colors[left_idx];
        if (v_left < v_above) {
            palette_add_to_cache(cache, &n, v_left);
            ++left_idx, --left_n;
        } else {
            palette_add_to_cache(cache, &n, v_above);
            ++above_idx, --above_n;
            if (v_left == v_above) ++left_idx, --left_n;
        }
    }
    while (above_n-- > 0) {
        uint16_t val = above_colors[above_idx++];
        palette_add_to_cache(cache, &n, val);
    }
    while (left_n-- > 0) {
        uint16_t val = left_colors[left_idx++];
        palette_add_to_cache(cache, &n, val);
    }
    assert(n <= 2 * PALETTE_MAX_SIZE);
    return n;
}
// Returns sub-sampled dimensions of the given block.
// The output values for 'rows_within_bounds' and 'cols_within_bounds' will
// differ from 'height' and 'width' when part of the block is outside the
// right
// and/or bottom image boundary.
void av1_get_block_dimensions(BlockSize bsize, int plane, const MacroBlockD *xd, int *width,
                              int *height, int *rows_within_bounds, int *cols_within_bounds) {
    const int block_height = block_size_high[bsize];
    const int block_width  = block_size_wide[bsize];
    const int block_rows =
        (xd->mb_to_bottom_edge >= 0) ? block_height : (xd->mb_to_bottom_edge >> 3) + block_height;
    const int block_cols =
        (xd->mb_to_right_edge >= 0) ? block_width : (xd->mb_to_right_edge >> 3) + block_width;

    uint8_t subsampling_x = plane == 0 ? 0 : 1;
    uint8_t subsampling_y = plane == 0 ? 0 : 1;

    assert(block_width >= block_cols);
    assert(block_height >= block_rows);
    const int plane_block_width  = block_width >> subsampling_x;
    const int plane_block_height = block_height >> subsampling_y;
    // Special handling for chroma sub8x8.
    const int is_chroma_sub8_x = plane > 0 && plane_block_width < 4;
    const int is_chroma_sub8_y = plane > 0 && plane_block_height < 4;
    if (width) *width = plane_block_width + 2 * is_chroma_sub8_x;
    if (height) *height = plane_block_height + 2 * is_chroma_sub8_y;
    if (rows_within_bounds) {
        *rows_within_bounds = (block_rows >> subsampling_y) + 2 * is_chroma_sub8_y;
    }
    if (cols_within_bounds) {
        *cols_within_bounds = (block_cols >> subsampling_x) + 2 * is_chroma_sub8_x;
    }
}

int av1_remove_duplicates(int *centroids, int num_centroids);
// Bias toward using colors in the cache.
// TODO: Try other schemes to improve compression.
static AOM_INLINE void optimize_palette_colors(uint16_t *color_cache, int n_cache, int n_colors,
                                               int stride, int *centroids) {
    if (n_cache <= 0) return;
    for (int i = 0; i < n_colors * stride; i += stride) {
        int min_diff = abs(centroids[i] - (int)color_cache[0]);
        int idx      = 0;
        for (int j = 1; j < n_cache; ++j) {
            const int this_diff = abs(centroids[i] - color_cache[j]);
            if (this_diff < min_diff) {
                min_diff = this_diff;
                idx      = j;
            }
        }
        if (min_diff <= 1) centroids[i] = color_cache[idx];
    }
}
// Extends 'color_map' array from 'orig_width x orig_height' to 'new_width x
// new_height'. Extra rows and columns are filled in by copying last valid
// row/column.
static AOM_INLINE void extend_palette_color_map(uint8_t *const color_map, int orig_width,
                                                int orig_height, int new_width, int new_height) {
    int j;
    assert(new_width >= orig_width);
    assert(new_height >= orig_height);
    if (new_width == orig_width && new_height == orig_height) return;

    for (j = orig_height - 1; j >= 0; --j) {
        memmove(color_map + j * new_width, color_map + j * orig_width, orig_width);
        // Copy last column to extra columns.
        memset(color_map + j * new_width + orig_width,
               color_map[j * new_width + orig_width - 1],
               new_width - orig_width);
    }
    // Copy last row to extra rows.
    for (j = orig_height; j < new_height; ++j) {
        memcpy(color_map + j * new_width, color_map + (orig_height - 1) * new_width, new_width);
    }
}
void palette_rd_y(PaletteInfo *palette_info, ModeDecisionContext *context_ptr, BlockSize bsize,
                  const int *data, int *centroids, int n, uint16_t *color_cache, int n_cache,
                  int bit_depth) {
    optimize_palette_colors(color_cache, n_cache, n, 1, centroids);
    int k = av1_remove_duplicates(centroids, n);
    if (k < PALETTE_MIN_SIZE) {
        // Too few unique colors to create a palette. And DC_PRED will work
        // well for that case anyway. So skip.
        palette_info->pmi.palette_size[0] = 0;
        return;
    }

    if (bit_depth > EB_8BIT) {
       for (int i = 0; i < k; ++i)
           palette_info->pmi.palette_colors[i] =
                   clip_pixel_highbd((int)centroids[i], bit_depth);
    } else {
        for (int i = 0; i < k; ++i)
            palette_info->pmi.palette_colors[i] = clip_pixel(centroids[i]);
    }
    palette_info->pmi.palette_size[0] = k;

    uint8_t *const color_map = palette_info->color_idx_map;
    int            block_width, block_height, rows, cols;
    av1_get_block_dimensions(
        bsize, 0, context_ptr->blk_ptr->av1xd, &block_width, &block_height, &rows, &cols);
    av1_calc_indices(data, centroids, color_map, rows * cols, k, 1);
    extend_palette_color_map(color_map, cols, rows, block_width, block_height);
}

int eb_av1_count_colors(const uint8_t *src, int stride, int rows, int cols, int *val_count);
int av1_count_colors_highbd(uint16_t *src, int stride, int rows, int cols, int bit_depth,
                            int *val_count);
/****************************************
   determine all palette luma candidates
 ****************************************/
void search_palette_luma(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                         PaletteInfo *palette_cand, uint32_t *tot_palette_cands) {
    int colors, n;
    EbBool is16bit = context_ptr->hbd_mode_decision > 0;

    EbPictureBufferDesc *src_pic = is16bit ?
        pcs_ptr->input_frame16bit :
        pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const int src_stride = src_pic->stride_y;

    const uint8_t *const src = src_pic->buffer_y +
          (((context_ptr->blk_origin_x + src_pic->origin_x) +
            (context_ptr->blk_origin_y + src_pic->origin_y) *
            src_pic->stride_y) << is16bit);
    int     block_width, block_height, rows, cols;

    Av1Common *  cm       = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD *xd       = context_ptr->blk_ptr->av1xd;
    TileInfo *   tile     = &context_ptr->sb_ptr->tile_info;
    BlockSize    bsize    = context_ptr->blk_geom->bsize;
    int32_t      mirow    = context_ptr->blk_origin_y >> 2;
    int32_t      micol    = context_ptr->blk_origin_x >> 2;
    xd->up_available      = (mirow > tile->mi_row_start);
    xd->left_available    = (micol > tile->mi_col_start);
    const int32_t bw      = mi_size_wide[bsize];
    const int32_t bh      = mi_size_high[bsize];
    xd->mb_to_top_edge    = -((mirow * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mirow) * MI_SIZE) * 8;
    xd->mb_to_left_edge   = -((micol * MI_SIZE) * 8);
    xd->mb_to_right_edge  = ((cm->mi_cols - bw - micol) * MI_SIZE) * 8;
    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end   = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end   = tile->mi_row_end;
    xd->n8_h              = bh;
    xd->n8_w              = bw;
    xd->is_sec_rect       = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((micol + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mirow & (xd->n8_w - 1)) xd->is_sec_rect = 1;

    int           mi_stride = cm->mi_stride;
    const int32_t offset    = mirow * mi_stride + micol;
    xd->mi                  = pcs_ptr->mi_grid_base + offset;
    ModeInfo *mi_ptr        = *xd->mi;
    if (xd->up_available) {
        xd->above_mbmi = &mi_ptr[-mi_stride].mbmi;
    } else
        xd->above_mbmi = NULL;
    if (xd->left_available) {
        xd->left_mbmi = &mi_ptr[-1].mbmi;
    } else
        xd->left_mbmi = NULL;

    av1_get_block_dimensions(context_ptr->blk_geom->bsize,
                             0,
                             context_ptr->blk_ptr->av1xd,
                             &block_width,
                             &block_height,
                             &rows,
                             &cols);

    int count_buf[1 << 12]; // Maximum (1 << 12) color levels.

    unsigned bit_depth = pcs_ptr->parent_pcs_ptr->scs_ptr->encoder_bit_depth;
    if (is16bit)
        colors = av1_count_colors_highbd((uint16_t *)src, src_stride, rows, cols,
            bit_depth, count_buf);
    else
        colors = eb_av1_count_colors(src, src_stride, rows, cols, count_buf);

    if (colors > 1 && colors <= 64) {
        int        r, c, i;
        const int  max_itr = 50;
        int *const data    = context_ptr->palette_buffer.kmeans_data_buf;
        int        centroids[PALETTE_MAX_SIZE];
        int lb, ub;

        #define GENERATE_KMEANS_DATA(src_data_type)                           \
            do {                                                              \
                lb = ub = ((src_data_type)src)[0];                            \
                for (r = 0; r < rows; ++r) {                                  \
                    for (c = 0; c < cols; ++c) {                              \
                        int val = ((src_data_type)src)[r * src_stride + c];   \
                        data[r * cols + c] = val;                             \
                        if (val < lb)                                         \
                            lb = val;                                         \
                        else if (val > ub)                                    \
                            ub = val;                                         \
                    }                                                         \
                }                                                             \
            } while (0)

        if (is16bit)
            GENERATE_KMEANS_DATA(uint16_t *);
        else
            GENERATE_KMEANS_DATA(uint8_t *);

        uint16_t  color_cache[2 * PALETTE_MAX_SIZE];
        const int n_cache = eb_get_palette_cache(xd, 0, color_cache);

        // Find the dominant colors, stored in top_colors[].
        int top_colors[PALETTE_MAX_SIZE] = {0};
        for (i = 0; i < AOMMIN(colors, PALETTE_MAX_SIZE); ++i) {
            int max_count = 0;
            for (int j = 0; j < (1 << bit_depth); ++j) {
                if (count_buf[j] > max_count) {
                    max_count     = count_buf[j];
                    top_colors[i] = j;
                }
            }
            assert(max_count > 0);
            count_buf[top_colors[i]] = 0;
        }

        // Try the dominant colors directly.
        // TODO: Try to avoid duplicate computation in cases
        // where the dominant colors and the k-means results are similar.

#if CS2_ADOPTIONS_1
        int step = (pcs_ptr->parent_pcs_ptr->palette_mode == 6)
                       ? 2
                       : 1;
#else
        int step = (pcs_ptr->parent_pcs_ptr->palette_mode == 6 && pcs_ptr->temporal_layer_index > 0)
                       ? 2
                       : 1;
#endif
        for (n = AOMMIN(colors, PALETTE_MAX_SIZE); n >= 2; n -= step) {
            for (i = 0; i < n; ++i) centroids[i] = top_colors[i];

            palette_rd_y(&palette_cand[*tot_palette_cands], context_ptr, bsize, data,
                         centroids, n, color_cache, n_cache, bit_depth);

            //consider this candidate if it has some non zero palette
            if (palette_cand[*tot_palette_cands].pmi.palette_size[0] > 2) (*tot_palette_cands)++;
            assert((*tot_palette_cands) <= 14);
        }

#if !CS2_ADOPTIONS_1
        if (pcs_ptr->parent_pcs_ptr->palette_mode == 3)
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == 0) return;

        if (pcs_ptr->parent_pcs_ptr->palette_mode == 5 ||
            pcs_ptr->parent_pcs_ptr->palette_mode == 6)
            if (pcs_ptr->temporal_layer_index > 0) return;
#endif
        // K-means clustering.
        for (n = AOMMIN(colors, PALETTE_MAX_SIZE); n >= 2; --n) {
            if (colors == PALETTE_MIN_SIZE) {
                // Special case: These colors automatically become the centroids.
                assert(colors == n);
                assert(colors == 2);
                centroids[0] = lb;
                centroids[1] = ub;
            } else {
                for (i = 0; i < n; ++i) { centroids[i] = lb + (2 * i + 1) * (ub - lb) / n / 2; }
                uint8_t *const color_map = palette_cand[*tot_palette_cands].color_idx_map;
                av1_k_means(data, centroids, color_map, rows * cols, n, 1, max_itr);
            }

            palette_rd_y(&palette_cand[*tot_palette_cands], context_ptr, bsize, data,
                         centroids, n, color_cache, n_cache, bit_depth);

            //consider this candidate if it has some non zero palette
            if (palette_cand[*tot_palette_cands].pmi.palette_size[0] > 2) (*tot_palette_cands)++;

            assert((*tot_palette_cands) <= 14);
        }
    }
}

/* clang-format off */
 typedef AomCdfProb(*MapCdf)[PALETTE_COLOR_INDEX_CONTEXTS]
     [CDF_SIZE(PALETTE_COLORS)];
 typedef const int(*ColorCost)[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS]
     [PALETTE_COLORS];
/* clang-format on */

typedef struct {
    int       rows;
    int       cols;
    int       n_colors;
    int       plane_width;
    uint8_t * color_map;
    MapCdf    map_cdf;
    ColorCost color_cost;
} Av1ColorMapParam;

static void get_palette_params(FRAME_CONTEXT *frame_context, BlkStruct *blk_ptr, int plane,
                               BlockSize bsize, Av1ColorMapParam *params) {
    const MacroBlockD *const     xd   = blk_ptr->av1xd;
    MbModeInfo *                 mbmi = &(xd->mi[0]->mbmi);
    const PaletteModeInfo *const pmi  = &mbmi->palette_mode_info;
    params->color_map                 = blk_ptr->palette_info.color_idx_map;
    params->map_cdf                   = plane ? frame_context->palette_uv_color_index_cdf
                            : frame_context->palette_y_color_index_cdf;
    params->color_cost = NULL;
    params->n_colors   = pmi->palette_size[plane];
    av1_get_block_dimensions(
        bsize, plane, xd, &params->plane_width, NULL, &params->rows, &params->cols);
}

static void get_color_map_params(FRAME_CONTEXT *frame_context, BlkStruct *blk_ptr, int plane,
                                 BlockSize bsize, TxSize tx_size, COLOR_MAP_TYPE type,
                                 Av1ColorMapParam *params) {
    (void)tx_size;
    memset(params, 0, sizeof(*params));
    switch (type) {
    case PALETTE_MAP: get_palette_params(frame_context, blk_ptr, plane, bsize, params); break;
    default: assert(0 && "Invalid color map type"); return;
    }
}
static void get_palette_params_rate(PaletteInfo *palette_info, MdRateEstimationContext *rate_table,
                                    BlkStruct *blk_ptr, int plane, BlockSize bsize,
                                    Av1ColorMapParam *params) {
    const MacroBlockD *const     xd  = blk_ptr->av1xd;
    const PaletteModeInfo *const pmi = &palette_info->pmi;

    params->color_map  = palette_info->color_idx_map;
    params->map_cdf    = NULL;
    params->color_cost = plane ? NULL : &rate_table->palette_ycolor_fac_bitss;
    params->n_colors   = pmi->palette_size[plane];

    av1_get_block_dimensions(
        bsize, plane, xd, &params->plane_width, NULL, &params->rows, &params->cols);
}

static void get_color_map_params_rate(PaletteInfo *                             palette_info,
                                      MdRateEstimationContext *                 rate_table,
                                      /*const MACROBLOCK *const x*/ BlkStruct *blk_ptr, int plane,
                                      BlockSize bsize, COLOR_MAP_TYPE type,
                                      Av1ColorMapParam *params) {
    memset(params, 0, sizeof(*params));
    switch (type) {
    case PALETTE_MAP:
        get_palette_params_rate(palette_info, rate_table, blk_ptr, plane, bsize, params);
        break;
    default: assert(0 && "Invalid color map type"); return;
    }
}

static int cost_and_tokenize_map(Av1ColorMapParam *param, TOKENEXTRA **t, int plane, int calc_rate,
                                 int allow_update_cdf, MapCdf map_pb_cdf) {
    const uint8_t *const color_map         = param->color_map;
    MapCdf               map_cdf           = param->map_cdf;
    ColorCost            color_cost        = param->color_cost;
    const int            plane_block_width = param->plane_width;
    const int            rows              = param->rows;
    const int            cols              = param->cols;
    const int            n                 = param->n_colors;
    const int            palette_size_idx  = n - PALETTE_MIN_SIZE;
    int                  this_rate         = 0;
#if !PALETTE_SPEEDUP
    uint8_t              color_order[PALETTE_MAX_SIZE];
#endif

    (void)plane;

    for (int k = 1; k < rows + cols - 1; ++k) {
        for (int j = AOMMIN(k, cols - 1); j >= AOMMAX(0, k - rows + 1); --j) {
            int       i = k - j;
            int       color_new_idx;
#if PALETTE_SPEEDUP
            const int color_ctx = av1_get_palette_color_index_context_optimized(
                color_map, plane_block_width, i, j, n, &color_new_idx);
#else
            const int color_ctx = av1_get_palette_color_index_context(
                color_map, plane_block_width, i, j, n, color_order, &color_new_idx);
#endif
            assert(color_new_idx >= 0 && color_new_idx < n);
            if (calc_rate) {
                this_rate += (*color_cost)[palette_size_idx][color_ctx][color_new_idx];
            } else {
                (*t)->token         = color_new_idx;
                (*t)->color_map_cdf = map_pb_cdf[palette_size_idx][color_ctx];
                ++(*t);
                if (allow_update_cdf)
                    update_cdf(map_cdf[palette_size_idx][color_ctx], color_new_idx, n);
#if CONFIG_ENTROPY_STATS
                if (plane) {
                    ++counts->palette_uv_color_index[palette_size_idx][color_ctx][color_new_idx];
                } else {
                    ++counts->palette_y_color_index[palette_size_idx][color_ctx][color_new_idx];
                }
#endif
            }
        }
    }
#if PALETTE_SPEEDUP
    return this_rate;
#else
    if (calc_rate) return this_rate;
    return 0;
#endif
}

void av1_tokenize_color_map(FRAME_CONTEXT *frame_context, BlkStruct *blk_ptr, int plane,
                            TOKENEXTRA **t, BlockSize bsize, TxSize tx_size, COLOR_MAP_TYPE type,
                            int allow_update_cdf) {
    assert(plane == 0 || plane == 1);
    Av1ColorMapParam color_map_params;
    get_color_map_params(frame_context, blk_ptr, plane, bsize, tx_size, type, &color_map_params);
    // The first color index does not use context or entropy.
    (*t)->token         = color_map_params.color_map[0];
    (*t)->color_map_cdf = NULL;
    ++(*t);
    MapCdf map_pb_cdf = plane ? frame_context->palette_uv_color_index_cdf
                              : frame_context->palette_y_color_index_cdf;
    cost_and_tokenize_map(&color_map_params, t, plane, 0, allow_update_cdf, map_pb_cdf);
}
int av1_cost_color_map(PaletteInfo *palette_info, MdRateEstimationContext *rate_table,
                       BlkStruct *blk_ptr, int plane, BlockSize bsize, COLOR_MAP_TYPE type) {
    assert(plane == 0 || plane == 1);
    Av1ColorMapParam color_map_params;
    get_color_map_params_rate(
        palette_info, rate_table, blk_ptr, plane, bsize, type, &color_map_params);
    MapCdf map_pb_cdf = NULL;
    return cost_and_tokenize_map(&color_map_params, NULL, plane, 1, 0, map_pb_cdf);
}
