#include <stdio.h>
#include <stdlib.h>

#include "EbPictureControlSet.h"
#include "EbResize.h"
#include "EbEncDecProcess.h"
#include "EbPictureDecisionProcess.h"
#include "EbPictureBufferDesc.h"

uint16_t svt_aom_get_max_can_count(EncMode enc_mode);
void     svt_aom_md_pme_search_controls(ModeDecisionContext *ctx, uint8_t md_pme_level);
void     svt_aom_set_inter_intra_ctrls(ModeDecisionContext *ctx, uint8_t inter_intra_level);

void    svt_aom_set_txt_controls(ModeDecisionContext *ctx, uint8_t txt_level);
void    svt_aom_set_obmc_controls(ModeDecisionContext *ctx, uint8_t obmc_mode);
void    svt_aom_set_wm_controls(ModeDecisionContext *ctx, uint8_t wm_level);
uint8_t svt_aom_set_nic_controls(ModeDecisionContext *ctx, uint8_t nic_level);
uint8_t svt_aom_set_chroma_controls(ModeDecisionContext *ctx, uint8_t uv_level);
uint8_t svt_aom_get_update_cdf_level(EncMode enc_mode, SliceType is_islice, uint8_t is_base);
uint8_t svt_aom_get_chroma_level(EncMode enc_mode);
uint8_t svt_aom_get_bypass_encdec(EncMode enc_mode, uint8_t hbd_md, uint8_t encoder_bit_depth);
uint8_t svt_aom_get_nic_level(EncMode enc_mode, uint8_t is_base, uint8_t hierarchical_levels, bool rtc_tune);

void    svt_aom_set_depth_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t depth_level);
uint8_t svt_aom_get_enable_me_16x16(EncMode enc_mode, bool rtc_tune);

Bool    svt_aom_is_ref_same_size(PictureControlSet *pcs, uint8_t list_idx, uint8_t ref_idx);
uint8_t svt_aom_get_enable_me_8x8(EncMode enc_mode, bool rtc_tune);
void    svt_aom_set_tpl_extended_controls(PictureParentControlSet *pcs, uint8_t tpl_level);
void    svt_aom_sig_deriv_mode_decision_config(SequenceControlSet *scs, PictureControlSet *pcs);
void    svt_aom_first_pass_sig_deriv_mode_decision_config(PictureControlSet *pcs);
void    svt_aom_sig_deriv_block(PictureControlSet *pcs, ModeDecisionContext *ctx);
void    svt_aom_sig_deriv_pre_analysis_pcs(PictureParentControlSet *pcs);
void    svt_aom_first_pass_sig_deriv_pre_analysis_pcs(PictureParentControlSet *pcs);
void    svt_aom_sig_deriv_pre_analysis_scs(SequenceControlSet *scs);
void    svt_aom_first_pass_sig_deriv_pre_analysis_scs(SequenceControlSet *scs);
void    svt_aom_sig_deriv_multi_processes(SequenceControlSet *scs, PictureParentControlSet *pcs,
                                          PictureDecisionContext *context_ptr);
void    svt_aom_first_pass_sig_deriv_multi_processes(SequenceControlSet *scs, PictureParentControlSet *pcs);
void    svt_aom_sig_deriv_me_tf(PictureParentControlSet *pcs, MeContext *me_ctx);
void    svt_aom_first_pass_sig_deriv_me(SequenceControlSet *scs, PictureParentControlSet *pcs, MeContext *me_ctx);

void svt_aom_sig_deriv_enc_dec_light_pd1(PictureControlSet *pcs, ModeDecisionContext *ctx);
void svt_aom_sig_deriv_enc_dec_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx);
void svt_aom_sig_deriv_enc_dec_common(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx);

void svt_aom_sig_deriv_me(SequenceControlSet *scs, PictureParentControlSet *pcs, MeContext *me_ctx);

void svt_aom_sig_deriv_enc_dec(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx);
#if TUNE_M6
bool svt_aom_need_gm_ref_info(EncMode enc_mode, bool super_res_off);
#else
bool    svt_aom_need_gm_ref_info(EncMode enc_mode, uint8_t is_base, bool super_res_off);
#endif
uint8_t svt_aom_derive_gm_level(PictureParentControlSet *pcs, bool super_res_off);

void svt_aom_set_gm_controls(PictureParentControlSet *pcs, uint8_t gm_level);

uint8_t svt_aom_get_enable_sg(EncMode enc_mode, uint8_t input_resolution, Bool fast_decode);

uint8_t svt_aom_get_enable_restoration(EncMode enc_mode, int8_t config_enable_restoration, uint8_t input_resolution,
                                       Bool fast_decode);

void svt_aom_set_dist_based_ref_pruning_controls(ModeDecisionContext *ctx, uint8_t dist_based_ref_pruning_level);

bool svt_aom_get_disallow_4x4(EncMode enc_mode, SliceType slice_type);

uint8_t svt_aom_get_nsq_level(EncMode enc_mode, uint8_t is_islice, uint8_t is_base, InputCoeffLvl coeff_lvl);
uint8_t get_inter_compound_level(EncMode enc_mode);
uint8_t get_filter_intra_level(EncMode enc_mode);
uint8_t svt_aom_get_inter_intra_level(EncMode enc_mode, uint8_t is_base, uint8_t transition_present);
#if TUNE_M6
uint8_t svt_aom_get_obmc_level(EncMode enc_mode, uint8_t fast_decode, EbInputResolution input_resolution);
#else
uint8_t svt_aom_get_obmc_level(EncMode enc_mode, uint8_t is_ref, uint8_t fast_decode,
                               EbInputResolution input_resolution);
#endif

void    svt_aom_set_nsq_ctrls(ModeDecisionContext *ctx, uint8_t nsq_level, uint8_t *allow_HVA_HVB, uint8_t *allow_HV4,
                              uint8_t *min_nsq_bsize);
uint8_t svt_aom_get_tpl_synthesizer_block_size(int8_t tpl_level, uint32_t picture_width, uint32_t picture_height);
