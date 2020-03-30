/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInterPrediction_h
#define EbInterPrediction_h

//#include "EbModeDecision.h"
#include "EbBlockStructures.h"
#include "EbAv1Structs.h"
#include "EbObject.h"
#include "EbMcp.h"
#include "filter.h"
#include "convolve.h"
#include "EbCabacContextModel.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REF_SCALE_SHIFT 14
#define REF_NO_SCALE (1 << REF_SCALE_SHIFT)
#define REF_INVALID_SCALE -1

#define RS_SUBPEL_BITS 6
#define RS_SUBPEL_MASK ((1 << RS_SUBPEL_BITS) - 1)
#define RS_SCALE_SUBPEL_BITS 14
#define RS_SCALE_SUBPEL_MASK ((1 << RS_SCALE_SUBPEL_BITS) - 1)
#define RS_SCALE_EXTRA_BITS (RS_SCALE_SUBPEL_BITS - RS_SUBPEL_BITS)
#define RS_SCALE_EXTRA_OFF (1 << (RS_SCALE_EXTRA_BITS - 1))
#define MV_BORDER (16 << 3) // Allow 16 pels in 1/8th pel units


#define NELEMENTS(x) (int)(sizeof(x) / sizeof(x[0]))

#define INTRABC_DELAY_PIXELS 256 //  Delay of 256 pixels
#define INTRABC_DELAY_SB64 (INTRABC_DELAY_PIXELS / 64)

// HW does not support < 4x4 prediction. To limit the bandwidth requirement, if
// block-size of current plane is smaller than 8x8, always only blend with the
// left neighbor(s) (skip blending with the above side).
#define DISABLE_CHROMA_U8X8_OBMC 0 // 0: one-sided obmc; 1: disable

    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_4[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8sharp[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8smooth[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, bilinear_filters[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_4smooth[SUBPEL_SHIFTS]);

    typedef struct SubpelParams {
        int32_t xs;
        int32_t ys;
        int32_t subpel_x;
        int32_t subpel_y;
    } SubpelParams;

    typedef struct PadBlock {
        int x0;
        int x1;
        int y0;
        int y1;
    } PadBlock;

#define INTERINTRA_WEDGE_SIGN 0

typedef uint8_t *WedgeMasksType[MAX_WEDGE_TYPES];

// Angles are with respect to horizontal anti-clockwise
typedef enum WedgeDirectionType
{
    WEDGE_HORIZONTAL = 0,
    WEDGE_VERTICAL = 1,
    WEDGE_OBLIQUE27 = 2,
    WEDGE_OBLIQUE63 = 3,
    WEDGE_OBLIQUE117 = 4,
    WEDGE_OBLIQUE153 = 5,
    WEDGE_DIRECTIONS
} WedgeDirectionType;

static INLINE void clamp_mv(MV *mv, int32_t min_col, int32_t max_col, int32_t min_row,
                            int32_t max_row) {
    mv->col = (int16_t)clamp(mv->col, min_col, max_col);
    mv->row = (int16_t)clamp(mv->row, min_row, max_row);
}

// 3-tuple: {direction, x_offset, y_offset}
typedef struct WedgeCodeType
{
    WedgeDirectionType direction;
    int32_t x_offset;
    int32_t y_offset;
} WedgeCodeType;

typedef struct WedgeParamsType
{
    int32_t bits;
    const WedgeCodeType *codebook;
    uint8_t *signflip;
    WedgeMasksType *masks;
} WedgeParamsType;


    struct ModeDecisionContext;

    typedef struct InterPredictionContext {
        EbDctor                               dctor;
        MotionCompensationPredictionContext  *mcp_context;
    } InterPredictionContext;

    void svt_inter_predictor(const uint8_t *src, int32_t src_stride,
        uint8_t *dst, int32_t dst_stride, const SubpelParams *subpel_params,
        const ScaleFactors *sf, int32_t w, int32_t h, ConvolveParams *conv_params,
        InterpFilters interp_filters, int32_t is_intrabc);

    void svt_highbd_inter_predictor(const uint16_t *src, int32_t src_stride,
        uint16_t *dst, int32_t dst_stride, const SubpelParams *subpel_params,
        const ScaleFactors *sf, int32_t w, int32_t h, ConvolveParams *conv_params,
        InterpFilters interp_filters, int32_t is_intrabc, int32_t bd);


    void av1_dist_wtd_comp_weight_assign(
        SeqHeader *seq_header,
        int cur_frame_index,
        int bck_frame_index,
        int fwd_frame_index,
        int compound_idx,
        int order_idx,
        int *fwd_offset, int *bck_offset,
        int *use_dist_wtd_comp_avg,
        int is_compound);

    void build_masked_compound_no_round(uint8_t *dst, int dst_stride,
        const CONV_BUF_TYPE *src0, int src0_stride,
        const CONV_BUF_TYPE *src1, int src1_stride,
        const InterInterCompoundData *const comp_data, uint8_t *seg_mask,
        BlockSize sb_type, int h, int w, ConvolveParams *conv_params,
        uint8_t bd, EbBool is_16bit);

    void av1_get_convolve_filter_params(uint32_t interp_filters,
        InterpFilterParams *params_x, InterpFilterParams *params_y,
        int32_t w, int32_t h);

    /* Mapping of interintra to intra mode for use in the intra component */
    static const PredictionMode interintra_to_intra_mode[INTERINTRA_MODES] = {
      DC_PRED, V_PRED, H_PRED, SMOOTH_PRED
    };

    void combine_interintra(InterIntraMode mode,
        int8_t use_wedge_interintra, int wedge_index,
        int wedge_sign, BlockSize bsize,
        BlockSize plane_bsize, uint8_t *comppred,
        int compstride, const uint8_t *interpred,
        int interstride, const uint8_t *intrapred,
        int intrastride);

    void combine_interintra_highbd(
        InterIntraMode mode, uint8_t use_wedge_interintra, uint8_t wedge_index,
        uint8_t wedge_sign, BlockSize bsize, BlockSize plane_bsize,
        uint8_t *comppred8, int compstride, const uint8_t *interpred8,
        int interstride, const uint8_t *intrapred8, int intrastride, int bd);


    void av1_setup_scale_factors_for_frame(ScaleFactors *sf, int other_w,
        int other_h, int this_w, int this_h);

    static INLINE int av1_is_valid_scale(const struct ScaleFactors *sf) {
        return sf->x_scale_fp != REF_INVALID_SCALE &&
            sf->y_scale_fp != REF_INVALID_SCALE;
    }
    static INLINE int av1_is_scaled(const struct ScaleFactors *sf) {
        return av1_is_valid_scale(sf) &&
            (sf->x_scale_fp != REF_NO_SCALE || sf->y_scale_fp != REF_NO_SCALE);
    }
    static INLINE int valid_ref_frame_size(int ref_width, int ref_height,
        int this_width, int this_height)
    {
        return 2 * this_width >= ref_width && 2 * this_height >= ref_height &&
            this_width <= 16 * ref_width && this_height <= 16 * ref_height;
    }
    MV32 av1_scale_mv(const MV *mvq4, int x, int y,
        const ScaleFactors *sf);


void build_smooth_interintra_mask(uint8_t *mask, int stride, BlockSize plane_bsize,
                                  InterIntraMode mode);

void highbd_convolve_2d_for_intrabc(const uint16_t *src, int src_stride, uint16_t *dst,
                                    int dst_stride, int w, int h, int subpel_x_q4,
                                    int subpel_y_q4, ConvolveParams *conv_params, int bd);

void convolve_2d_for_intrabc(const uint8_t *src, int src_stride, uint8_t *dst,
                             int dst_stride, int w, int h, int subpel_x_q4, int subpel_y_q4,
                             ConvolveParams *conv_params);

extern aom_highbd_convolve_fn_t convolve_hbd[/*sub_x*/2][/*sub_y*/2][/*bi*/2];

extern AomConvolveFn convolve[/*sub_x*/2][/*sub_y*/2][/*bi*/2];

int is_interintra_wedge_used(BlockSize sb_type);

int32_t get_wedge_bits_lookup(BlockSize sb_type);

const uint8_t *av1_get_contiguous_soft_mask(int wedge_index, int wedge_sign, BlockSize sb_type);

int get_wedge_params_bits(BlockSize sb_type);

int is_masked_compound_type(COMPOUND_TYPE type);

// Although we assign 32 bit integers, all the values are strictly under 14
// bits.
static int div_mult[32] = {0,    16384, 8192, 5461, 4096, 3276, 2730, 2340, 2048, 1820, 1638,
                           1489, 1365,  1260, 1170, 1092, 1024, 963,  910,  862,  819,  780,
                           744,  712,   682,  655,  630,  606,  585,  564,  546,  528};

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

static INLINE void get_mv_projection(MV *output, MV ref, int num, int den) {
    den              = AOMMIN(den, MAX_FRAME_DISTANCE);
    num              = num > 0 ? AOMMIN(num, MAX_FRAME_DISTANCE) : AOMMAX(num, -MAX_FRAME_DISTANCE);
    const int mv_row = ROUND_POWER_OF_TWO_SIGNED(ref.row * num * div_mult[den], 14);
    const int mv_col = ROUND_POWER_OF_TWO_SIGNED(ref.col * num * div_mult[den], 14);
    const int clamp_max = MV_UPP - 1;
    const int clamp_min = MV_LOW + 1;
    output->row         = (int16_t)clamp(mv_row, clamp_min, clamp_max);
    output->col         = (int16_t)clamp(mv_col, clamp_min, clamp_max);
}

static INLINE int check_sb_border(const int mi_row, const int mi_col, const int row_offset,
                                  const int col_offset) {
    const int sb_mi_size = mi_size_wide[BLOCK_64X64];
    const int row        = mi_row & (sb_mi_size - 1);
    const int col        = mi_col & (sb_mi_size - 1);

    if (row + row_offset < 0 || row + row_offset >= sb_mi_size || col + col_offset < 0 ||
        col + col_offset >= sb_mi_size)
        return 0;

    return 1;
}

static INLINE int is_neighbor_overlappable(const MbModeInfo *mbmi) {
    return /*is_intrabc_block(mbmi) ||*/ mbmi->block_mi.ref_frame[0] >
                                         INTRA_FRAME; // TODO: modify when add intra_bc
}

static INLINE int32_t is_mv_valid(const MV *mv) {
    return mv->row > MV_LOW && mv->row < MV_UPP && mv->col > MV_LOW && mv->col < MV_UPP;
}

#define CHECK_BACKWARD_REFS(ref_frame) \
    (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)

void av1_find_ref_dv(IntMv *ref_dv, const TileInfo *const tile, int mib_size, int mi_row,
                     int mi_col);

static INLINE int32_t is_comp_ref_allowed(BlockSize bsize) {
    return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}

static INLINE int is_interinter_compound_used(CompoundType type, BlockSize sb_type) {
    const int comp_allowed = is_comp_ref_allowed(sb_type);
    switch (type) {
        case COMPOUND_AVERAGE:
        case COMPOUND_DISTWTD:
        case COMPOUND_DIFFWTD: return comp_allowed;
        case COMPOUND_WEDGE: return comp_allowed && get_wedge_params_bits(sb_type) > 0;
        default: assert(0); return 0;
    }
}

static INLINE int is_any_masked_compound_used(BlockSize sb_type) {
    CompoundType comp_type;
    int          i;
    if (!is_comp_ref_allowed(sb_type)) return 0;
    for (i = 0; i < COMPOUND_TYPES; i++) {
        comp_type = (CompoundType)i;
        if (is_masked_compound_type(comp_type) && is_interinter_compound_used(comp_type, sb_type))
            return 1;
    }
    return 0;
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

static INLINE int bsize_to_max_depth(BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    int    depth   = 0;
    while (depth < MAX_TX_DEPTH && tx_size != TX_4X4) {
        depth++;
        tx_size = sub_tx_size_map[tx_size];
    }
    return depth;
}

static INLINE PredictionMode compound_ref0_mode(PredictionMode mode) {
    static PredictionMode lut[] = {
            MB_MODE_COUNT, // DC_PRED
            MB_MODE_COUNT, // V_PRED
            MB_MODE_COUNT, // H_PRED
            MB_MODE_COUNT, // D45_PRED
            MB_MODE_COUNT, // D135_PRED
            MB_MODE_COUNT, // D113_PRED
            MB_MODE_COUNT, // D157_PRED
            MB_MODE_COUNT, // D203_PRED
            MB_MODE_COUNT, // D67_PRED
            MB_MODE_COUNT, // SMOOTH_PRED
            MB_MODE_COUNT, // SMOOTH_V_PRED
            MB_MODE_COUNT, // SMOOTH_H_PRED
            MB_MODE_COUNT, // PAETH_PRED
            MB_MODE_COUNT, // NEARESTMV
            MB_MODE_COUNT, // NEARMV
            MB_MODE_COUNT, // GLOBALMV
            MB_MODE_COUNT, // NEWMV
            NEARESTMV, // NEAREST_NEARESTMV
            NEARMV, // NEAR_NEARMV
            NEARESTMV, // NEAREST_NEWMV
            NEWMV, // NEW_NEARESTMV
            NEARMV, // NEAR_NEWMV
            NEWMV, // NEW_NEARMV
            GLOBALMV, // GLOBAL_GLOBALMV
            NEWMV, // NEW_NEWMV
    };
    assert(NELEMENTS(lut) == MB_MODE_COUNT);
    assert(is_inter_compound_mode(mode));
    return lut[mode];
}

static INLINE PredictionMode compound_ref1_mode(PredictionMode mode) {
    static PredictionMode lut[] = {
            MB_MODE_COUNT, // DC_PRED
            MB_MODE_COUNT, // V_PRED
            MB_MODE_COUNT, // H_PRED
            MB_MODE_COUNT, // D45_PRED
            MB_MODE_COUNT, // D135_PRED
            MB_MODE_COUNT, // D113_PRED
            MB_MODE_COUNT, // D157_PRED
            MB_MODE_COUNT, // D203_PRED
            MB_MODE_COUNT, // D67_PRED
            MB_MODE_COUNT, // SMOOTH_PRED
            MB_MODE_COUNT, // SMOOTH_V_PRED
            MB_MODE_COUNT, // SMOOTH_H_PRED
            MB_MODE_COUNT, // PAETH_PRED
            MB_MODE_COUNT, // NEARESTMV
            MB_MODE_COUNT, // NEARMV
            MB_MODE_COUNT, // GLOBALMV
            MB_MODE_COUNT, // NEWMV
            NEARESTMV, // NEAREST_NEARESTMV
            NEARMV, // NEAR_NEARMV
            NEWMV, // NEAREST_NEWMV
            NEARESTMV, // NEW_NEARESTMV
            NEWMV, // NEAR_NEWMV
            NEARMV, // NEW_NEARMV
            GLOBALMV, // GLOBAL_GLOBALMV
            NEWMV, // NEW_NEWMV
    };
    assert(NELEMENTS(lut) == MB_MODE_COUNT);
    assert(is_inter_compound_mode(mode));
    return lut[mode];
}

static INLINE EbBool is_motion_variation_allowed_bsize(const BlockSize bsize) {
    return (block_size_wide[bsize] >= 8 && block_size_high[bsize] >= 8);
}

static INLINE int is_global_mv_block(const PredictionMode mode, const BlockSize bsize,
                                     TransformationType type) {
    return (mode == GLOBALMV || mode == GLOBAL_GLOBALMV) && type > TRANSLATION &&
           is_motion_variation_allowed_bsize(bsize);
}

static INLINE uint32_t have_nearmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE int is_intrabc_block(const BlockModeInfo *block_mi) { return block_mi->use_intrabc; }

static INLINE int is_inter_block(const BlockModeInfo *bloc_mi) {
    return is_intrabc_block(bloc_mi) || bloc_mi->ref_frame[0] > INTRA_FRAME;
}

void av1_set_ref_frame(MvReferenceFrame *rf, int8_t ref_frame_type);

int av1_skip_u4x4_pred_in_obmc(BlockSize bsize, int dir, int subsampling_x, int subsampling_y);

#ifdef __cplusplus
}
#endif
#endif //EbInterPrediction_h
