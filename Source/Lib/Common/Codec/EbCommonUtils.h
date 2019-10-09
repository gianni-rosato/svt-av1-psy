#include "EbCodingUnit.h"
#include "EbDefinitions.h"


#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units
#define MAX_OFFSET_WIDTH 64
#define MAX_OFFSET_HEIGHT 0

#define INTRABC_DELAY_PIXELS 256  //  Delay of 256 pixels
#define INTRABC_DELAY_SB64  (INTRABC_DELAY_PIXELS / 64)

#define NELEMENTS(x) (int)(sizeof(x) / sizeof(x[0]))

#define CHECK_BACKWARD_REFS(ref_frame) \
  (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)


#define NZ_MAP_CTX_0 SIG_COEF_CONTEXTS_2D
#define NZ_MAP_CTX_5 (NZ_MAP_CTX_0 + 5)
#define NZ_MAP_CTX_10 (NZ_MAP_CTX_0 + 10)

static const int nz_map_ctx_offset_1d[32] = {
  NZ_MAP_CTX_0,  NZ_MAP_CTX_5,  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10,
};

static const int16_t eb_k_eob_group_start[12] = { 0, 1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513 };
static const int16_t eb_k_eob_offset_bits[12] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

// Although we assign 32 bit integers, all the values are strictly under 14
// bits.
static int div_mult[32] = { 0,    16384, 8192, 5461, 4096, 3276, 2730, 2340,
                            2048, 1820,  1638, 1489, 1365, 1260, 1170, 1092,
                            1024, 963,   910,  862,  819,  780,  744,  712,
                            682,  655,   630,  606,  585,  564,  546,  528 };

static const PredictionMode fimode_to_intradir[FILTER_INTRA_MODES] = {
    DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};

static INLINE uint8_t *set_levels(uint8_t *const levels_buf, const int32_t width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

static INLINE int get_padded_idx(const int idx, const int bwl) {
    return idx + ((idx >> bwl) << TX_PAD_HOR_LOG2);
}

static INLINE int get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide_log2[tx_size];
}

static INLINE int get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}

static INLINE int get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}

static INLINE int bsize_to_max_depth(BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    int depth = 0;
    while (depth < MAX_TX_DEPTH && tx_size != TX_4X4) {
        depth++;
        tx_size = sub_tx_size_map[tx_size];
    }
    return depth;
}

static INLINE int bsize_to_tx_size_cat(BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    assert(tx_size != TX_4X4);
    int depth = 0;
    while (tx_size != TX_4X4) {
        depth++;
        tx_size = sub_tx_size_map[tx_size];
        assert(depth < 10);
    }
    assert(depth <= MAX_TX_CATS);
    return depth - 1;
}

static INLINE void integer_mv_precision(MV *mv) {
    int mod = (mv->row % 8);
    if (mod != 0) {
        mv->row -= mod;
        if (abs(mod) > 4) {
            if (mod > 0)
                mv->row += 8;
            else
                mv->row -= 8;
        }
    }

    mod = (mv->col % 8);
    if (mod != 0) {
        mv->col -= mod;
        if (abs(mod) > 4) {
            if (mod > 0)
                mv->col += 8;
            else
                mv->col -= 8;
        }
    }
}

static INLINE void lower_mv_precision(MV *mv, int allow_hp, int is_integer) {
    if (is_integer)
        integer_mv_precision(mv);
    else {
        if (!allow_hp) {
            if (mv->row & 1) mv->row += (mv->row > 0 ? -1 : 1);
            if (mv->col & 1) mv->col += (mv->col > 0 ? -1 : 1);
        }
    }
}

static INLINE int get_palette_bsize_ctx(BlockSize bsize) {
    return num_pels_log2_lookup[bsize] - num_pels_log2_lookup[BLOCK_8X8];
}

static INLINE PredictionMode get_uv_mode(UvPredictionMode mode) {
    assert(mode < UV_INTRA_MODES);
    static const PredictionMode uv2y[] = {
      DC_PRED,        // UV_DC_PRED
      V_PRED,         // UV_V_PRED
      H_PRED,         // UV_H_PRED
      D45_PRED,       // UV_D45_PRED
      D135_PRED,      // UV_D135_PRED
      D113_PRED,      // UV_D113_PRED
      D157_PRED,      // UV_D157_PRED
      D203_PRED,      // UV_D203_PRED
      D67_PRED,       // UV_D67_PRED
      SMOOTH_PRED,    // UV_SMOOTH_PRED
      SMOOTH_V_PRED,  // UV_SMOOTH_V_PRED
      SMOOTH_H_PRED,  // UV_SMOOTH_H_PRED
      PAETH_PRED,     // UV_PAETH_PRED
      DC_PRED,        // UV_CFL_PRED
      INTRA_INVALID,  // UV_INTRA_MODES
      INTRA_INVALID,  // UV_MODE_INVALID
    };
    return uv2y[mode];
}

static AOM_FORCE_INLINE int get_br_ctx_eob(const int c,  // raster order
    const int bwl, const TxClass tx_class)
{
    const int row = c >> bwl;
    const int col = c - (row << bwl);
    if (c == 0) return 0;
    if ((tx_class == TX_CLASS_2D && row < 2 && col < 2) ||
        (tx_class == TX_CLASS_HORIZ && col == 0) ||
        (tx_class == TX_CLASS_VERT && row == 0))
        return 7;
    return 14;
}

static INLINE void get_mv_projection(MV *output, MV ref, int num, int den) {
    den = AOMMIN(den, MAX_FRAME_DISTANCE);
    num = num > 0 ? AOMMIN(num, MAX_FRAME_DISTANCE)
        : AOMMAX(num, -MAX_FRAME_DISTANCE);
    const int mv_row =
        ROUND_POWER_OF_TWO_SIGNED(ref.row * num * div_mult[den], 14);
    const int mv_col =
        ROUND_POWER_OF_TWO_SIGNED(ref.col * num * div_mult[den], 14);
    const int clamp_max = MV_UPP - 1;
    const int clamp_min = MV_LOW + 1;
    output->row = (int16_t)clamp(mv_row, clamp_min, clamp_max);
    output->col = (int16_t)clamp(mv_col, clamp_min, clamp_max);
}

static INLINE int32_t get_br_ctx(const uint8_t *const levels,
    const int32_t c,  // raster order
    const int32_t bwl, const TxType tx_type) {
    const int32_t row = c >> bwl;
    const int32_t col = c - (row << bwl);
    const int32_t stride = (1 << bwl) + TX_PAD_HOR;
    const TxClass tx_class = tx_type_to_class[tx_type];
    const int32_t pos = row * stride + col;
    int32_t mag = levels[pos + 1];
    mag += levels[pos + stride];
    switch (tx_class) {
    case TX_CLASS_2D:
        mag += levels[pos + stride + 1];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if ((row < 2) && (col < 2)) return mag + 7;
        break;
    case TX_CLASS_HORIZ:
        mag += levels[pos + 2];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (col == 0) return mag + 7;
        break;
    case TX_CLASS_VERT:
        mag += levels[pos + (stride << 1)];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (row == 0) return mag + 7;
        break;
    default: break;
    }

    return mag + 14;
}

static INLINE int get_lower_levels_ctx_eob(int bwl, int height, int scan_idx) {
    if (scan_idx == 0) return 0;
    if (scan_idx <= (height << bwl) / 8) return 1;
    if (scan_idx <= (height << bwl) / 4) return 2;
    return 3;
}

static AomCdfProb cdf_element_prob(const AomCdfProb *const cdf,
    size_t element) {
    assert(cdf != NULL);
    return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
}


static INLINE void partition_gather_horz_alike(AomCdfProb *out,
    const AomCdfProb *const in,
    BlockSize bsize) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_HORZ);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
    out[2] = 0;
}

static INLINE void partition_gather_vert_alike(AomCdfProb *out,
    const AomCdfProb *const in,
    BlockSize bsize) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_VERT);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
    out[2] = 0;
}

static INLINE int check_sb_border(const int mi_row, const int mi_col,
    const int row_offset, const int col_offset) {
    const int sb_mi_size = mi_size_wide[BLOCK_64X64];
    const int row = mi_row & (sb_mi_size - 1);
    const int col = mi_col & (sb_mi_size - 1);

    if (row + row_offset < 0 || row + row_offset >= sb_mi_size ||
        col + col_offset < 0 || col + col_offset >= sb_mi_size)
        return 0;

    return 1;
}

static INLINE PredictionMode compound_ref0_mode(PredictionMode mode) {
    static PredictionMode lut[] = {
      MB_MODE_COUNT,  // DC_PRED
      MB_MODE_COUNT,  // V_PRED
      MB_MODE_COUNT,  // H_PRED
      MB_MODE_COUNT,  // D45_PRED
      MB_MODE_COUNT,  // D135_PRED
      MB_MODE_COUNT,  // D113_PRED
      MB_MODE_COUNT,  // D157_PRED
      MB_MODE_COUNT,  // D203_PRED
      MB_MODE_COUNT,  // D67_PRED
      MB_MODE_COUNT,  // SMOOTH_PRED
      MB_MODE_COUNT,  // SMOOTH_V_PRED
      MB_MODE_COUNT,  // SMOOTH_H_PRED
      MB_MODE_COUNT,  // PAETH_PRED
      MB_MODE_COUNT,  // NEARESTMV
      MB_MODE_COUNT,  // NEARMV
      MB_MODE_COUNT,  // GLOBALMV
      MB_MODE_COUNT,  // NEWMV
      NEARESTMV,      // NEAREST_NEARESTMV
      NEARMV,         // NEAR_NEARMV
      NEARESTMV,      // NEAREST_NEWMV
      NEWMV,          // NEW_NEARESTMV
      NEARMV,         // NEAR_NEWMV
      NEWMV,          // NEW_NEARMV
      GLOBALMV,       // GLOBAL_GLOBALMV
      NEWMV,          // NEW_NEWMV
    };
    assert(NELEMENTS(lut) == MB_MODE_COUNT);
    assert(is_inter_compound_mode(mode));
    return lut[mode];
}

static INLINE PredictionMode compound_ref1_mode(PredictionMode mode) {
    static PredictionMode lut[] = {
      MB_MODE_COUNT,  // DC_PRED
      MB_MODE_COUNT,  // V_PRED
      MB_MODE_COUNT,  // H_PRED
      MB_MODE_COUNT,  // D45_PRED
      MB_MODE_COUNT,  // D135_PRED
      MB_MODE_COUNT,  // D113_PRED
      MB_MODE_COUNT,  // D157_PRED
      MB_MODE_COUNT,  // D203_PRED
      MB_MODE_COUNT,  // D67_PRED
      MB_MODE_COUNT,  // SMOOTH_PRED
      MB_MODE_COUNT,  // SMOOTH_V_PRED
      MB_MODE_COUNT,  // SMOOTH_H_PRED
      MB_MODE_COUNT,  // PAETH_PRED
      MB_MODE_COUNT,  // NEARESTMV
      MB_MODE_COUNT,  // NEARMV
      MB_MODE_COUNT,  // GLOBALMV
      MB_MODE_COUNT,  // NEWMV
      NEARESTMV,      // NEAREST_NEARESTMV
      NEARMV,         // NEAR_NEARMV
      NEWMV,          // NEAREST_NEWMV
      NEARESTMV,      // NEW_NEARESTMV
      NEWMV,          // NEAR_NEWMV
      NEARMV,         // NEW_NEARMV
      GLOBALMV,       // GLOBAL_GLOBALMV
      NEWMV,          // NEW_NEWMV
    };
    assert(NELEMENTS(lut) == MB_MODE_COUNT);
    assert(is_inter_compound_mode(mode));
    return lut[mode];
}

static INLINE int32_t is_chroma_reference(int32_t mi_row, int32_t mi_col, BlockSize bsize,
    int32_t subsampling_x, int32_t subsampling_y) {
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];
    int32_t ref_pos = ((mi_row & 0x01) || !(bh & 0x01) || !subsampling_y) &&
        ((mi_col & 0x01) || !(bw & 0x01) || !subsampling_x);
    return ref_pos;
}