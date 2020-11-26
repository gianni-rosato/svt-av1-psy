/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "EbAppContext.h"
#include "EbAppConfig.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbAppInputy4m.h"
#include "EbTime.h"
/***************************************
 * Macros
 ***************************************/
#define CLIP3(min_val, max_val, a) \
    (((a) < (min_val)) ? (min_val) : (((a) > (max_val)) ? (max_val) : (a)))
#define FUTURE_WINDOW_WIDTH 4
#define SIZE_OF_ONE_FRAME_IN_BYTES(width, height, csp, is_16bit) \
    ((((width) * (height)) + 2 * (((width) * (height)) >> (3 - csp))) << is_16bit)
#define YUV4MPEG2_IND_SIZE 9
extern volatile int32_t keep_running;

/***************************************
* Process Error Log
***************************************/
void log_error_output(FILE *error_log_file, uint32_t error_code) {
    switch (error_code) {
        // EB_ENC_AMVP_ERRORS:
    case EB_ENC_AMVP_ERROR1:
        fprintf(error_log_file,
                "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR2:
        fprintf(error_log_file,
                "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

    case EB_ENC_AMVP_ERROR3:
        fprintf(error_log_file,
                "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR4:
        fprintf(error_log_file,
                "Error: The availability parameter in GetSpatialMVPPosAx() function can not be > 3 "
                "!\n");
        break;

    case EB_ENC_AMVP_ERROR5:
        fprintf(error_log_file,
                "Error: The availability parameter in GetSpatialMVPPosBx() function can not be > 7 "
                "!\n");
        break;

    case EB_ENC_AMVP_NULL_REF_ERROR:
        fprintf(error_log_file, "Error: The reference_object can not be NULL!\n");
        break;

    case EB_ENC_AMVP_SPATIAL_NA_ERROR:
        fprintf(error_log_file,
                "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

        // EB_ENC_CL_ERRORS:
    case EB_ENC_CL_ERROR1:
        fprintf(error_log_file, "Error: Unknown Inter Prediction Direction!\n");
        break;

    case EB_ENC_CL_ERROR2: fprintf(error_log_file, "Error: Unknown coding mode!\n"); break;

    case EB_ENC_CL_ERROR3:
        fprintf(error_log_file, "Error: Mode Decision Candidate Buffer Overflow!\n");
        break;

    case EB_ENC_CL_ERROR4:
        fprintf(error_log_file, "Error: Too many Mode Decision Fast Candidates!\n");
        break;

    case EB_ENC_CL_ERROR5:
        fprintf(error_log_file,
                "Error: Too many buffers chosen for this level by PreModeDecision!\n");
        break;

    case EB_ENC_CL_ERROR6:
        fprintf(error_log_file,
                "Error: Ping-Pong structure needs at least two buffers to work properly!\n");
        break;

    case EB_ENC_CL_ERROR7: fprintf(error_log_file, "Error: Invalid Intra Partition\n"); break;

    case EB_ENC_CL_ERROR8: fprintf(error_log_file, "Error: Invalid TU Configuration\n"); break;
    case EB_ENC_CL_ERROR9:
        fprintf(error_log_file, "Error: Invalid Prediction Mode\n");
        break;

        // EB_ENC_DLF_ERRORS:
    case EB_ENC_DLF_ERROR1:
        fprintf(error_log_file, "Error: While calculating bS for DLF!\n");
        break;

    case EB_ENC_DLF_ERROR2:
        fprintf(error_log_file, "Error: Unknown Inter Prediction Direction Combination!\n");
        break;

    case EB_ENC_DLF_ERROR3:
        fprintf(error_log_file, "Error: If any PU is in I_MODE, the bS will be 2!\n");
        break;

    case EB_ENC_DLF_ERROR4:
        fprintf(error_log_file,
                "Error: The x/y location of the CU must be the multiple of minmum CU size!");
        break;

    case EB_ENC_DLF_ERROR5: fprintf(error_log_file, "Error: Unknown Slice Type!"); break;

    case EB_ENC_DLF_ERROR6:
        fprintf(error_log_file,
                "Error: While calculating the bS for the PU bounday, the 4x4 block must be "
                "guaranteed at the PU boundary!");
        break;

    case EB_ENC_DLF_ERROR7: fprintf(error_log_file, "Error: SB size must be power of 2!"); break;

    case EB_ENC_DLF_ERROR8:
        fprintf(error_log_file,
                "Error: Deblocking filter can not support the picture whose width or height is not "
                "the multiple of 8!");
        break;

    case EB_ENC_DLF_ERROR9: fprintf(error_log_file, "Error: Neighbor PU must be available!"); break;

    case EB_ENC_DLF_ERROR10:
        fprintf(error_log_file,
                "Error: Deblocking filter can not support the picture whose width or height is not "
                "the multiple of 8!");
        break;

        // EB_ENC_EC_ERRORS:
    case EB_ENC_EC_ERROR1:
        fprintf(error_log_file, "Error: EncodeCodedBlockFlags: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR10:
        fprintf(error_log_file, "Error: EncodeTuSplitCoeff: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR11:
        fprintf(error_log_file,
                "Error: CodeSPS: Long term reference pictures are not currently handled!\n");
        break;

    case EB_ENC_EC_ERROR12:
        fprintf(error_log_file,
                "Error: CodeProfileTierLevel: The maximum sublayers must be equal to 1!\n");
        break;

    case EB_ENC_EC_ERROR13:
        fprintf(error_log_file, "Error: EncodeRootCodedBlockFlag: rootCbf too large!\n");
        break;

    case EB_ENC_EC_ERROR14:
        fprintf(error_log_file,
                "Error: cpbCountMinus1 in HRD parameter exceeds the upper limit 4!\n");
        break;

    case EB_ENC_EC_ERROR15:
        fprintf(
            error_log_file,
            "Error: numDecodingUnitsMinus1 in picture timeing SEI exceeds the upper limit 64!\n");
        break;

    case EB_ENC_EC_ERROR16:
        fprintf(error_log_file,
                "Error: The size of the unregistered user data SEI payload is not allowed!\n");
        break;

    case EB_ENC_EC_ERROR2:
        fprintf(error_log_file,
                "Error: copy_payload: output buffer too small!\n");
        break;

    case EB_ENC_EC_ERROR3: fprintf(error_log_file, "Error: encode_sb: Unknown mode type!\n"); break;

    case EB_ENC_EC_ERROR4:
        fprintf(error_log_file, "Error: 8x4 & 4x8 PU should not have Bi-pred mode!\n");
        break;

    case EB_ENC_EC_ERROR5:
        fprintf(error_log_file, "Error: EncodeMergeIndex: value too large!\n");
        break;

    case EB_ENC_EC_ERROR6:
        fprintf(error_log_file, "Error: EncodeSkipFlag: context too large!\n");
        break;

    case EB_ENC_EC_ERROR7:
        fprintf(error_log_file, "Error: EncodeBypassBins: binsLength must be less than 32!\n");
        break;

    case EB_ENC_EC_ERROR8:
        fprintf(error_log_file, "Error: EncodeQuantizedCoefficients: Invalid block size!\n");
        break;

    case EB_ENC_EC_ERROR9:
        fprintf(error_log_file, "Error: EncodeSplitFlag: context too large!\n");
        break;

    case EB_ENC_EC_ERROR26: fprintf(error_log_file, "Error: Level not recognized!\n"); break;

    case EB_ENC_EC_ERROR27:
        fprintf(error_log_file, "Error: EncodeOneBin:  BinaryValue must be less than 2\n");
        break;

    case EB_ENC_EC_ERROR28: fprintf(error_log_file, "Error: No more than 6 SAO types\n"); break;

    case EB_ENC_EC_ERROR29:
        fprintf(error_log_file, "Error: No more than 6 SAO types\n");
        break;

        // EB_ENC_FL_ERRORS:
    case EB_ENC_FL_ERROR1: fprintf(error_log_file, "Error: Uncovered area inside Cu!\n"); break;

    case EB_ENC_FL_ERROR2:
        fprintf(error_log_file, "Error: Depth 2 is not allowed for 8x8 CU!\n");
        break;

    case EB_ENC_FL_ERROR3:
        fprintf(error_log_file, "Error: Depth 0 is not allowed for 64x64 CU!\n");
        break;

    case EB_ENC_FL_ERROR4:
        fprintf(error_log_file, "Error: Max CU Depth Exceeded!\n");
        break;

        // EB_ENC_HANDLE_ERRORS:
    case EB_ENC_HANDLE_ERROR1:
        fprintf(error_log_file, "Error: Only one Resource Coordination Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR10:
        fprintf(error_log_file, "Error: Need at least one Entropy Coding Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR11:
        fprintf(error_log_file, "Error: Only one Packetization Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR12:
        fprintf(error_log_file,
                "Error: RC Results Fifo Size should be greater than RC Tasks Fifo Size in order to "
                "avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR13:
        fprintf(error_log_file,
                "Error: RC Tasks Fifo Size should be greater than EC results Fifo Size in order to "
                "avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR14:
        fprintf(error_log_file,
                "Error: RC Tasks Fifo Size should be greater than Picture Manager results Fifo "
                "Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR18:
        fprintf(error_log_file, "Error: Intra period setting breaks mini-gop!\n");
        break;

    case EB_ENC_HANDLE_ERROR2:
        fprintf(error_log_file, "Error: Only one Picture Enhancement Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR3:
        fprintf(error_log_file, "Error: Only one Picture Manager Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR4:
        fprintf(error_log_file, "Error: Need at least one ME Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR5:
        fprintf(error_log_file, "Error: Only one Rate-Control Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR6:
        fprintf(error_log_file, "Error: Need at least one Mode Decision Configuration Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR7:
        fprintf(error_log_file, "Error: Need at least one Coding Loop Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR8:
        fprintf(error_log_file, "Error: Only one Second Pass Deblocking Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR9:
        fprintf(error_log_file, "Error: Only one ALF Process allowed!\n");
        break;

        // EB_ENC_INTER_ERRORS:
    case EB_ENC_INTER_INVLD_MCP_ERROR:
        fprintf(error_log_file,
                "Error: Motion compensation prediction is out of the picture boundary!\n");
        break;

    case EB_ENC_INTER_PRED_ERROR0:
        fprintf(error_log_file, "Error: Unkown Inter Prediction Direction!\n");
        break;

    case EB_ENC_INTER_PRED_ERROR1:
        fprintf(error_log_file, "Error: Inter prediction can not support more than 2 MVPs!\n");
        break;

        // EB_ENC_INTRA_ERRORS:
    case EB_ENC_INTRA_PRED_ERROR1:
        fprintf(error_log_file, "Error: IntraPrediction does not support 2Nx2N partition size!\n");
        break;

    case EB_ENC_INTRA_PRED_ERROR2:
        fprintf(error_log_file,
                "Error: IntraPrediction: intra prediction only supports square PU!\n");
        break;

    case EB_ENC_INTRA_PRED_ERROR3:
        fprintf(error_log_file, "Error: IntraPredictionChroma: Only Planar!\n");
        break;

    case EB_ENC_INVLD_PART_SIZE_ERROR:
        fprintf(error_log_file,
                "Error: IntraPrediction: only PU sizes of 8 or largers are currently supported!\n");
        break;

        // EB_ENC_MD_ERRORS:
    case EB_ENC_MD_ERROR1:
        fprintf(error_log_file, "Error: Unknown AMVP Mode Decision Candidate Type!\n");
        break;

    case EB_ENC_MD_ERROR2:
        fprintf(error_log_file, "Error: PreModeDecision: need at least one buffer!\n");
        break;

    case EB_ENC_MD_ERROR3:
        fprintf(error_log_file, "Error: Unknow Inter Prediction Direction!\n");
        break;

    case EB_ENC_MD_ERROR4: fprintf(error_log_file, "Error: Unknown ME SAD Level!\n"); break;

    case EB_ENC_MD_ERROR5:
        fprintf(error_log_file,
                "Error: Invalid encoder mode. The encoder mode should be 0, 1 or 2!\n");
        break;

    case EB_ENC_MD_ERROR6: fprintf(error_log_file, "Error: Invalid TU size!\n"); break;

    case EB_ENC_MD_ERROR7: fprintf(error_log_file, "Error: Unknown depth!\n"); break;

    case EB_ENC_MD_ERROR8: fprintf(error_log_file, "Error: Depth not supported!\n"); break;

    case EB_ENC_MD_ERROR9:
        fprintf(error_log_file,
                "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

    case EB_ENC_MD_ERROR10:
        fprintf(error_log_file,
                "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

        // EB_ENC_ME_ERRORS:
    case EB_ENC_ME_ERROR1:
        fprintf(error_log_file,
                "Error: Motion Estimation: non valid value of the subPelDirection !\n");
        break;

    case EB_ENC_ME_ERROR2:
        fprintf(error_log_file,
                "Error: FillMvMergeCandidate() method only supports P or b slices!\n");
        break;

        // EB_ENC_ERRORS:
    case EB_ENC_ROB_OF_ERROR:
        fprintf(error_log_file, "Error: Recon Output Buffer Overflow!\n");
        break;

        // EB_ENC_PACKETIZATION_ERRORS:
    case EB_ENC_PACKETIZATION_ERROR1:
        fprintf(error_log_file,
                "Error: PacketizationProcess: Picture Number does not match entry. "
                "PacketizationReorderQueue overflow!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR2:
        fprintf(error_log_file,
                "Error: Entropy Coding Result can not be outputed by processes other than entropy "
                "coder and ALF!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR3:
        fprintf(error_log_file,
                "Error: The encoder can not support the SliceMode other than 0 and 1!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR4:
        fprintf(error_log_file, "Error: Statistics Output Buffer Overflow!\n");
        break;
    case EB_ENC_PACKETIZATION_ERROR5:
        fprintf(error_log_file,
                "Error: Stream Fifo is starving..deadlock, increase "
                "EB_outputStreamBufferFifoInitCount APP_ENCODERSTREAMBUFFERCOUNT \n");
        break;

        // EB_ENC_PM_ERRORS:
    case EB_ENC_PM_ERROR0:
        fprintf(error_log_file, "Error: PictureManagerProcess: Unknown Slice Type!\n");
        break;

    case EB_ENC_PM_ERROR1:
        fprintf(error_log_file, "Error: EbPictureManager: dependent_count underflow!\n");
        break;

    case EB_ENC_PM_ERROR10:
        fprintf(error_log_file,
                "Error: picture_manager_kernel: reference_entry_ptr should never be null!\n");
        break;

    case EB_ENC_PM_ERROR2:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: The Reference Structure period must be less than "
                "the MAX_ELAPSED_IDR_COUNT or false-IDR boundary logic will be activated!\n");
        break;

    case EB_ENC_PM_ERROR3:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: The dependent_count underflow detected!\n");
        break;

    case EB_ENC_PM_ERROR4:
        fprintf(error_log_file, "Error: PictureManagerProcess: Empty input queue!\n");
        break;

    case EB_ENC_PM_ERROR5:
        fprintf(error_log_file, "Error: PictureManagerProcess: Empty reference queue!\n");
        break;

    case EB_ENC_PM_ERROR6:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: The capped elaspedNonIdrCount must be larger than "
                "the maximum supported delta ref poc!\n");
        break;

    case EB_ENC_PM_ERROR7:
        fprintf(error_log_file, "Error: PictureManagerProcess: Reference Picture Queue Full!\n");
        break;

    case EB_ENC_PM_ERROR8:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: No reference match found - this will lead to a "
                "memory leak!\n");
        break;

    case EB_ENC_PM_ERROR9:
        fprintf(error_log_file, "Error: PictureManagerProcess: Unknown picture type!\n");
        break;

    case EB_ENC_PM_ERROR12:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: prediction structure configuration API has too many "
                "reference pictures\n");
        break;

    case EB_ENC_PM_ERROR13:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: The maximum allowed frame rate is 60 fps\n");
        break;

    case EB_ENC_PM_ERROR14:
        fprintf(error_log_file,
                "Error: PictureManagerProcess: The minimum allowed frame rate is 1 fps\n");
        break;

        // EB_ENC_PRED_STRC_ERRORS:
    case EB_ENC_PRED_STRC_ERROR1:
        fprintf(error_log_file, "Error: prediction_structure_ctor: DecodeOrder lut too small!\n");
        break;

    case EB_ENC_PRED_STRC_ERROR2:
        fprintf(error_log_file,
                "Error: prediction_structure_ctor: prediction structure improperly configured!\n");
        break;

        // EB_ENC_PU_ERRORS:
    case EB_ENC_PU_ERROR1: fprintf(error_log_file, "Error: Unknown partition size!\n"); break;

    case EB_ENC_PU_ERROR2:
        fprintf(error_log_file, "Error: The target area is not inside the CU!\n");
        break;

        // EB_ENC_RC_ERRORS:
    case EB_ENC_RC_ERROR1:
        fprintf(error_log_file, "Error: RateControlProcess: Unknown input task_type!\n");
        break;

    case EB_ENC_RC_ERROR2:
        fprintf(error_log_file, "Error: RateControlProcess: No RC interval found!\n");
        break;

    case EB_ENC_RC_ERROR3:
        fprintf(error_log_file, "Error: RateControlProcess: RC input Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR4:
        fprintf(error_log_file, "Error: RateControlProcess: RC feedback Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR5:
        fprintf(error_log_file, "Error: RateControlProcess: RC feedback Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR6:
        fprintf(error_log_file,
                "Error: RateControlProcess: No feedback frame match found - this will lead to a "
                "memory leak!\n");
        break;

    case EB_ENC_RC_ERROR7:
        fprintf(error_log_file, "Error: remainingBytes has to be multiple of 2 for 16 bit input\n");
        break;

    case EB_ENC_RC_ERROR8:
        fprintf(error_log_file, "Error: hlRateControlHistorgramQueue Overflow\n");
        break;

        // EB_ENC_RD_COST_ERRORS:
    case EB_ENC_RD_COST_ERROR1:
        fprintf(error_log_file, "Error: Skip mode only exists in 2Nx2N partition type!\n");
        break;

    case EB_ENC_RD_COST_ERROR2:
        fprintf(error_log_file, "Error: IntraChromaCost: Unknown slice type!\n");
        break;

    case EB_ENC_RD_COST_ERROR3:
        fprintf(error_log_file,
                "Error: intra2_nx2_n_fast_cost_islice can only support 2Nx2N partition type!\n");
        break;

        // EB_ENC_SAO_ERRORS:
    case EB_ENC_SAO_ERROR1: fprintf(error_log_file, "Error: No more than 6 SAO types!\n"); break;

    case EB_ENC_SAO_ERROR2:
        fprintf(error_log_file, "Error: No more than 5 EO SAO categories!\n");
        break;
        // EB_ENC_SCS_ERRORS:
    case EB_ENC_SCS_ERROR1:
        fprintf(error_log_file,
                "Error: SequenceControlSetCopy: Not all SequenceControlSet members are being "
                "copied!\n");
        break;

        // EB_ENC_BITSTREAM_ERRORS:
    case EB_ENC_BITSTREAM_ERROR1:
        fprintf(error_log_file,
                "Error: OutputBitstreamRBSPToPayload: Bitstream payload buffer empty!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR2:
        fprintf(error_log_file, "Error: OutputBitstreamWrite: Empty Bitstream!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR3:
        fprintf(error_log_file,
                "Error: OutputBitstreamRBSPToPayload: Buffer index more than buffer size!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR4:
        fprintf(error_log_file,
                "Error: OutputBitstreamRBSPToPayload: Start Location in not inside the buffer!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR5:
        fprintf(error_log_file,
                "Error: OutputBitstreamWrite: Trying to write more than one word!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR6:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Expecting Start code!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR7:
        fprintf(
            error_log_file,
            "Error: OutputBitstreamRBSPToPayload: Bitstream not flushed (i.e. byte-aligned)!\n");
        break;

    case EB_ENC_RESS_COOR_ERRORS1:
        fprintf(error_log_file,
                "Error: ResourceCoordinationProcess: The received input data should be equal to "
                "the buffer size - only complete frame transfer is supported\n");
        break;

    case EB_ENC_RES_COORD_InvalidQP:
        fprintf(error_log_file,
                "Error: ResourceCoordinationProcess: The QP value in the QP file is invalid\n");
        break;

    case EB_ENC_RES_COORD_InvalidSliceType:
        fprintf(error_log_file, "Error: ResourceCoordinationProcess: Slice Type Invalid\n");
        break;

        // picture decision Errors
    case EB_ENC_PD_ERROR8:
        fprintf(error_log_file,
                "Error: PictureDecisionProcess: Picture Decision Reorder Queue overflow\n");
        break;

    default: fprintf(error_log_file, "Error: Others!\n"); break;
    }

    return;
}

void read_input_frames(EbConfig *config, uint8_t is_16bit, EbBufferHeaderType *header_ptr) {
    const uint32_t input_padded_width  = config->input_padded_width;
    const uint32_t input_padded_height = config->input_padded_height;
    FILE *         input_file          = config->input_file;
    EbSvtIOFormat *input_ptr = (EbSvtIOFormat *)header_ptr->p_buffer;

    const uint8_t color_format  = config->config.encoder_color_format;
    const uint8_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;

    input_ptr->y_stride  = input_padded_width;
    input_ptr->cr_stride = input_padded_width >> subsampling_x;
    input_ptr->cb_stride = input_padded_width >> subsampling_x;

    if (config->buffered_input == -1) {
        uint64_t read_size;
        if (is_16bit == 0 || (is_16bit == 1 && config->config.compressed_ten_bit_format == 0)) {
            read_size = (uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(
                input_padded_width, input_padded_height, color_format, is_16bit);

            header_ptr->n_filled_len = 0;
            /* if input is a y4m file, read next line which contains "FRAME" */
            if (config->y4m_input == EB_TRUE) read_y4m_frame_delimiter(config);
            uint64_t luma_read_size = (uint64_t)input_padded_width * input_padded_height
                                      << is_16bit;
            uint8_t *eb_input_ptr = input_ptr->luma;
            if (!config->y4m_input && config->processed_frame_count == 0 &&
                (config->input_file == stdin || config->input_file_is_fifo)) {
                /* 9 bytes were already buffered during the the YUV4MPEG2 header probe */
                memcpy(eb_input_ptr, config->y4m_buf, YUV4MPEG2_IND_SIZE);
                header_ptr->n_filled_len += YUV4MPEG2_IND_SIZE;
                eb_input_ptr += YUV4MPEG2_IND_SIZE;
                header_ptr->n_filled_len += (uint32_t)fread(
                    eb_input_ptr, 1, luma_read_size - YUV4MPEG2_IND_SIZE, input_file);
            } else {
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->luma, 1, luma_read_size, input_file);
            }
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cb, 1, luma_read_size >> (3 - color_format), input_file);
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cr, 1, luma_read_size >> (3 - color_format), input_file);

            if (read_size != header_ptr->n_filled_len) {
                fseek(input_file, 0, SEEK_SET);
                if (config->y4m_input == EB_TRUE) {
                    read_and_skip_y4m_header(config);
                    read_y4m_frame_delimiter(config);
                }
                header_ptr->n_filled_len =
                    (uint32_t)fread(input_ptr->luma, 1, luma_read_size, input_file);
                header_ptr->n_filled_len += (uint32_t)fread(
                    input_ptr->cb, 1, luma_read_size >> (3 - color_format), input_file);
                header_ptr->n_filled_len += (uint32_t)fread(
                    input_ptr->cr, 1, luma_read_size >> (3 - color_format), input_file);
            }
        } else {
            assert(is_16bit == 1 && config->config.compressed_ten_bit_format == 1);
            // 10-bit Compressed Unpacked Mode
            const uint32_t luma_read_size        = input_padded_width * input_padded_height;
            const uint32_t chroma_read_size      = luma_read_size >> (3 - color_format);
            const uint32_t nbit_luma_read_size   = (input_padded_width / 4) * input_padded_height;
            const uint32_t nbit_chroma_read_size = nbit_luma_read_size >> (3 - color_format);

            // Fill the buffer with a complete frame
            header_ptr->n_filled_len = 0;

            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->luma, 1, luma_read_size, input_file);
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cb, 1, chroma_read_size, input_file);
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cr, 1, chroma_read_size, input_file);

            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->luma_ext, 1, nbit_luma_read_size, input_file);
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cb_ext, 1, nbit_chroma_read_size, input_file);
            header_ptr->n_filled_len +=
                (uint32_t)fread(input_ptr->cr_ext, 1, nbit_chroma_read_size, input_file);

            read_size = luma_read_size + nbit_luma_read_size +
                2 * (chroma_read_size + nbit_chroma_read_size);

            if (read_size != header_ptr->n_filled_len) {
                fseek(input_file, 0, SEEK_SET);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->luma, 1, luma_read_size, input_file);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->cb, 1, chroma_read_size, input_file);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->cr, 1, chroma_read_size, input_file);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->luma_ext, 1, nbit_luma_read_size, input_file);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->cb_ext, 1, nbit_chroma_read_size, input_file);
                header_ptr->n_filled_len +=
                    (uint32_t)fread(input_ptr->cr_ext, 1, nbit_chroma_read_size, input_file);
            }
        }

        if (feof(input_file) != 0) {
            if ((input_file == stdin) || (config->input_file_is_fifo)) {
                //for a fifo, we only know this when we reach eof
                config->frames_to_be_encoded = config->frames_encoded;
                if (header_ptr->n_filled_len != read_size) {
                    // not a completed frame
                    header_ptr->n_filled_len = 0;
                }
            } else {
                // If we reached the end of file, loop over again
                fseek(input_file, 0, SEEK_SET);
            }
        }

    } else {
        if (is_16bit && config->config.compressed_ten_bit_format == 1) {
            // Determine size of each plane
            const size_t luma_8bit_size   = input_padded_width * input_padded_height;
            const size_t chroma_8bit_size = luma_8bit_size >> (3 - color_format);
            const size_t luma_2bit_size   = luma_8bit_size / 4; //4-2bit pixels into 1 byte
            const size_t chroma_2bit_size = luma_2bit_size >> (3 - color_format);

            input_ptr = (EbSvtIOFormat *)header_ptr->p_buffer;
            input_ptr->y_stride      = input_padded_width;
            input_ptr->cr_stride     = input_padded_width >> subsampling_x;
            input_ptr->cb_stride     = input_padded_width >> subsampling_x;

            input_ptr->luma =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input];
            input_ptr->cb =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_8bit_size;
            input_ptr->cr =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_8bit_size + chroma_8bit_size;

            input_ptr->luma_ext =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_8bit_size + 2 * chroma_8bit_size;
            input_ptr->cb_ext =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_8bit_size + 2 * chroma_8bit_size + luma_2bit_size;
            input_ptr->cr_ext =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_8bit_size + 2 * chroma_8bit_size + luma_2bit_size + chroma_2bit_size;

            header_ptr->n_filled_len = (uint32_t)(luma_8bit_size + luma_2bit_size +
                                                  2 * (chroma_8bit_size + chroma_2bit_size));
        } else {
            //Normal unpacked mode:yuv420p10le yuv422p10le yuv444p10le
            const size_t luma_size   = (input_padded_width * input_padded_height) << is_16bit;
            const size_t chroma_size = luma_size >> (3 - color_format);

            input_ptr = (EbSvtIOFormat *)header_ptr->p_buffer;

            input_ptr->y_stride  = input_padded_width;
            input_ptr->cr_stride = input_padded_width >> subsampling_x;
            input_ptr->cb_stride = input_padded_width >> subsampling_x;

            input_ptr->luma =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input];
            input_ptr->cb =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_size;
            input_ptr->cr =
                config->sequence_buffer[config->processed_frame_count % config->buffered_input] +
                luma_size + chroma_size;

            header_ptr->n_filled_len = (uint32_t)(luma_size + 2 * chroma_size);
        }
    }

    return;
}

/**
 * Reads and extracts one qp from the qp_file
 * @param qp_file file to read a value from
 * @param qp_read_from_file boolean value to check if a value was read
 * @return long value of qp. -1 is returned if eof or eol is reached without reading anything.
 * A 0 may also be returned if a line starting with '#', '/', or '-' is found
 */
static long get_next_qp_from_qp_file(FILE *const qp_file, int *const qp_read_from_file) {
    long qp = 0;
    char line[512], *pos = line;
    // Read single line until \n
    if (!fgets(line, 512, qp_file))
        // eof
        return -1;
    // Clear out beginning spaces
    while (isspace(*pos)) ++pos;
    if (!*pos)
        // eol
        return -1;
    switch (*pos) {
    case '#':
    case '/':
    case '-': return 0;
    }
    if (isdigit(*pos))
        qp = strtol(pos, NULL, 0);
    if (qp > 0)
        *qp_read_from_file = 1;
    return qp;
}

static unsigned char send_qp_on_the_fly(FILE *const qp_file, uint8_t *use_qp_file) {
    long tmp_qp            = 0;
    int  qp_read_from_file = 0;

    while (tmp_qp == 0 || (tmp_qp == -1 && qp_read_from_file))
        // get next qp
        tmp_qp = get_next_qp_from_qp_file(qp_file, &qp_read_from_file);

    if (tmp_qp == -1) {
        *use_qp_file = EB_FALSE;
        fprintf(stderr, "\nWarning: QP File did not contain any valid QPs");
    }
    return (unsigned)CLIP3(0, 63, tmp_qp);
}

static void injector(uint64_t processed_frame_count, uint32_t injector_frame_rate) {
    static uint64_t start_times_seconds;
    static uint64_t start_timesu_seconds;
    static int      first_time = 0;

    if (first_time == 0) {
        first_time = 1;
        app_svt_av1_get_time(&start_times_seconds, &start_timesu_seconds);
    } else {
        uint64_t current_times_seconds, current_timesu_seconds;
        app_svt_av1_get_time(&current_times_seconds, &current_timesu_seconds);
        const double elapsed_time = app_svt_av1_compute_overall_elapsed_time(
            start_times_seconds,
            start_timesu_seconds,
            current_times_seconds,
            current_timesu_seconds);
        const int    buffer_frames     = 1; // How far ahead of time should we let it get
        const double injector_interval = (double)(1 << 16) /
            injector_frame_rate; // 1.0 / injector frame rate (in this
        // case, 1.0/encodRate)
        const double predicted_time  = (processed_frame_count - buffer_frames) * injector_interval;
        const int    milli_sec_ahead = (int)(1000 * (predicted_time - elapsed_time));
        if (milli_sec_ahead > 0)
            app_svt_av1_sleep(milli_sec_ahead);
    }
}

//************************************/
// process_input_buffer
// Reads yuv frames from file and copy
// them into the input buffer
/************************************/
void process_input_buffer(EncChannel* channel) {
    EbConfig *config = channel->config;
    EbAppContext *app_call_back = channel->app_callback;
    uint8_t             is_16bit         = (uint8_t)(config->config.encoder_bit_depth > 8);
    EbBufferHeaderType *header_ptr       = app_call_back->input_buffer_pool;
    EbComponentType *   component_handle = (EbComponentType *)app_call_back->svt_encoder_handle;

    AppExitConditionType return_value = APP_ExitConditionNone;

    const uint8_t color_format         = config->config.encoder_color_format;
    const int64_t input_padded_width   = config->input_padded_width;
    const int64_t input_padded_height  = config->input_padded_height;
    const int64_t frames_to_be_encoded = config->frames_to_be_encoded;
    int64_t       total_bytes_to_process_count;
    int64_t       remaining_byte_count;
    uint32_t      compressed10bit_frame_size =
        (uint32_t)((input_padded_width * input_padded_height) +
                   2 * ((input_padded_width * input_padded_width) >> (3 - color_format)));
    compressed10bit_frame_size += compressed10bit_frame_size / 4;

    if (channel->exit_cond_input != APP_ExitConditionNone)
        return;
    if (config->injector && config->processed_frame_count)
        injector(config->processed_frame_count, config->injector_frame_rate);
    total_bytes_to_process_count =
        (frames_to_be_encoded < 0)
            ? -1
            : (config->config.encoder_bit_depth == 10 && config->config.compressed_ten_bit_format == 1)
                  ? frames_to_be_encoded * (int64_t)compressed10bit_frame_size
                  : frames_to_be_encoded *
                        SIZE_OF_ONE_FRAME_IN_BYTES(
                            input_padded_width, input_padded_height, color_format, is_16bit);

    remaining_byte_count =
        (total_bytes_to_process_count < 0)
            ? -1
            : total_bytes_to_process_count - (int64_t)config->processed_byte_count;

    // If there are bytes left to encode, configure the header
    if (remaining_byte_count != 0 && config->stop_encoder == EB_FALSE) {
        read_input_frames(config, is_16bit, header_ptr);
        if (header_ptr->n_filled_len) {
            // Update the context parameters
            config->processed_byte_count += header_ptr->n_filled_len;
            header_ptr->p_app_private = (EbPtr)NULL;
            config->frames_encoded    = (int32_t)(++config->processed_frame_count);

            // Configuration parameters changed on the fly
            if (config->config.use_qp_file && config->qp_file)
                header_ptr->qp = send_qp_on_the_fly(config->qp_file, &config->config.use_qp_file);

            if (keep_running == 0 && !config->stop_encoder) config->stop_encoder = EB_TRUE;
            // Fill in Buffers Header control data
            header_ptr->pts      = config->processed_frame_count - 1;
            header_ptr->pic_type = EB_AV1_INVALID_PICTURE;
            header_ptr->flags    = 0;

            // Send the picture
            svt_av1_enc_send_picture(component_handle, header_ptr);
        }

        if ((config->processed_frame_count == (uint64_t)config->frames_to_be_encoded) ||
            config->stop_encoder) {
            header_ptr->n_alloc_len   = 0;
            header_ptr->n_filled_len  = 0;
            header_ptr->n_tick_count  = 0;
            header_ptr->p_app_private = NULL;
            header_ptr->flags         = EB_BUFFERFLAG_EOS;
            header_ptr->p_buffer      = NULL;
            header_ptr->pic_type      = EB_AV1_INVALID_PICTURE;

            svt_av1_enc_send_picture(component_handle, header_ptr);
        }

        return_value =
            (header_ptr->flags == EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : return_value;
    }

    channel->exit_cond_input = return_value;
}

#define LONG_ENCODE_FRAME_ENCODE 4000
#define SPEED_MEASUREMENT_INTERVAL 2000
#define START_STEADY_STATE 1000
#define AV1_FOURCC 0x31305641 // used for ivf header
#define IVF_STREAM_HEADER_SIZE 32
#define IVF_FRAME_HEADER_SIZE 12
#define OBU_FRAME_HEADER_SIZE 3
#define TD_SIZE 2
static __inline void mem_put_le32(void *vmem, int32_t val) {
    uint8_t *mem = (uint8_t *)vmem;

    mem[0] = (uint8_t)((val >> 0) & 0xff);
    mem[1] = (uint8_t)((val >> 8) & 0xff);
    mem[2] = (uint8_t)((val >> 16) & 0xff);
    mem[3] = (uint8_t)((val >> 24) & 0xff);
}
#define MEM_VALUE_T_SZ_BITS (sizeof(MEM_VALUE_T) << 3)

static __inline void mem_put_le16(void *vmem, int32_t val) {
    uint8_t *mem = (uint8_t *)vmem;

    mem[0] = (uint8_t)((val >> 0) & 0xff);
    mem[1] = (uint8_t)((val >> 8) & 0xff);
}

static void write_ivf_stream_header(EbConfig *config) {
    char header[IVF_STREAM_HEADER_SIZE];
    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4, 0); // version
    mem_put_le16(header + 6, 32); // header size
    mem_put_le32(header + 8, AV1_FOURCC); // fourcc
    mem_put_le16(header + 12, config->input_padded_width); // width
    mem_put_le16(header + 14, config->input_padded_height); // height
    if (config->config.frame_rate_denominator != 0 && config->config.frame_rate_numerator != 0) {
        mem_put_le32(header + 16, config->config.frame_rate_numerator); // rate
        mem_put_le32(header + 20, config->config.frame_rate_denominator); // scale
            //mem_put_le32(header + 16, config->frame_rate_denominator);  // rate
            //mem_put_le32(header + 20, config->frame_rate_numerator);  // scale
    } else {
        mem_put_le32(header + 16, (config->config.frame_rate >> 16) * 1000); // rate
        mem_put_le32(header + 20, 1000); // scale
            //mem_put_le32(header + 16, config->frame_rate_denominator);  // rate
            //mem_put_le32(header + 20, config->frame_rate_numerator);  // scale
    }
    mem_put_le32(header + 24, 0); // length
    mem_put_le32(header + 28, 0); // unused
    //config->performance_context.byte_count += 32;
    if (config->bitstream_file) fwrite(header, 1, IVF_STREAM_HEADER_SIZE, config->bitstream_file);

    return;
}

static void write_ivf_frame_header(EbConfig *config, uint32_t byte_count) {
    char    header[IVF_FRAME_HEADER_SIZE];
    int32_t write_location = 0;

    mem_put_le32(&header[write_location], (int32_t)byte_count);
    write_location = write_location + 4;
    mem_put_le32(&header[write_location], (int32_t)((config->ivf_count) & 0xFFFFFFFF));
    write_location = write_location + 4;
    mem_put_le32(&header[write_location], (int32_t)((config->ivf_count) >> 32));

    config->byte_count_since_ivf = (byte_count);

    config->ivf_count++;
    fflush(stdout);

    if (config->bitstream_file) fwrite(header, 1, IVF_FRAME_HEADER_SIZE, config->bitstream_file);
}
double get_psnr(double sse, double max) {
    double psnr;
    if (sse == 0)
        psnr = 10 * log10(max / (double)0.1);
    else
        psnr = 10 * log10(max / sse);

    return psnr;
}

/***************************************
* Process Output STATISTICS Buffer
***************************************/
void process_output_statistics_buffer(EbBufferHeaderType *header_ptr, EbConfig *config) {
    uint32_t max_luma_value = (config->config.encoder_bit_depth == 8) ? 255 : 1023;
    uint64_t picture_stream_size, luma_sse, cr_sse, cb_sse, picture_number, picture_qp;
    double   luma_ssim, cr_ssim, cb_ssim;
    double   temp_var, luma_psnr, cb_psnr, cr_psnr;
    uint32_t source_width = config->config.source_width;
    uint32_t source_height = config->config.source_height;

    picture_stream_size = header_ptr->n_filled_len;
    luma_sse            = header_ptr->luma_sse;
    cr_sse              = header_ptr->cr_sse;
    cb_sse              = header_ptr->cb_sse;
    picture_number      = header_ptr->pts;
    picture_qp          = header_ptr->qp;
    luma_ssim           = header_ptr->luma_ssim;
    cr_ssim             = header_ptr->cr_ssim;
    cb_ssim             = header_ptr->cb_ssim;

    temp_var =
        (double)max_luma_value * max_luma_value * (source_width * source_height);

    luma_psnr = get_psnr((double)luma_sse, temp_var);

    temp_var = (double)max_luma_value * max_luma_value *
               (source_width / 2 * source_height / 2);

    cb_psnr = get_psnr((double)cb_sse, temp_var);

    cr_psnr = get_psnr((double)cr_sse, temp_var);

    config->performance_context.sum_luma_psnr += luma_psnr;
    config->performance_context.sum_cr_psnr += cr_psnr;
    config->performance_context.sum_cb_psnr += cb_psnr;

    config->performance_context.sum_luma_sse += luma_sse;
    config->performance_context.sum_cr_sse += cr_sse;
    config->performance_context.sum_cb_sse += cb_sse;

    config->performance_context.sum_qp += picture_qp;
    config->performance_context.sum_luma_ssim    += luma_ssim;
    config->performance_context.sum_cr_ssim      += cr_ssim;
    config->performance_context.sum_cb_ssim      += cb_ssim;

    // Write statistic Data to file
    if (config->stat_file) {

        fprintf(config->stat_file,
                "Picture Number: %4d\t QP: %4d  [ "
                "PSNR-Y: %.2f dB,\tPSNR-U: %.2f dB,\tPSNR-V: %.2f "
                "dB,\tMSE-Y: %.2f,\tMSE-U: %.2f,\tMSE-V: %.2f,\t"
                "SSIM-Y: %.5f,\tSSIM-U: %.5f,\tSSIM-V: %.5f"
                " ]\t %6d bytes\n",
                (int)picture_number,
                (int)picture_qp,
                luma_psnr,
                cb_psnr,
                cr_psnr,
                (double)luma_sse / (source_width * source_height),
                (double)cb_sse / (source_width / 2 * source_height / 2),
                (double)cr_sse / (source_width / 2 * source_height / 2),
                luma_ssim,
                cr_ssim,
                cb_ssim,
                (int)picture_stream_size);
    }

    return;
}

void process_output_stream_buffer(EncChannel* channel, EncApp* enc_app,
                                                int32_t *frame_count) {
    EbConfig *config = channel->config;
    EbAppContext *app_call_back = channel->app_callback;
    AppPortActiveType *  port_state = &app_call_back->output_stream_port_active;
    EbBufferHeaderType * header_ptr;
    EbComponentType *    component_handle = (EbComponentType *)app_call_back->svt_encoder_handle;
    AppExitConditionType return_value     = APP_ExitConditionNone;
    // Per channel variables
    FILE *stream_file = config->bitstream_file;

    uint64_t *total_latency = &config->performance_context.total_latency;
    uint32_t *max_latency   = &config->performance_context.max_latency;

    // Local variables
    uint64_t finish_s_time = 0;
    uint64_t finish_u_time = 0;
    uint8_t  is_alt_ref    = 1;
    if (channel->exit_cond_output != APP_ExitConditionNone)
        return;
    uint8_t pic_send_done = (channel->exit_cond_input == APP_ExitConditionNone) ||
            (channel->exit_cond_recon == APP_ExitConditionNone)
            ? 0
            : 1;
    while (is_alt_ref) {
        is_alt_ref = 0;
        // non-blocking call until all input frames are sent
        EbErrorType stream_status = svt_av1_enc_get_packet(
            component_handle, &header_ptr, pic_send_done);

        if (stream_status == EB_ErrorMax) {
            fprintf(stderr, "\n");
            log_error_output(config->error_log_file, header_ptr->flags);
            channel->exit_cond_output = APP_ExitConditionError;
            return;
        } else if (stream_status != EB_NoErrorEmptyQueue) {
            uint32_t flags = header_ptr->flags;
            is_alt_ref        = (flags & EB_BUFFERFLAG_IS_ALT_REF);
            if (!(flags & EB_BUFFERFLAG_IS_ALT_REF))
                ++(config->performance_context.frame_count);
            *total_latency += (uint64_t)header_ptr->n_tick_count;
            *max_latency =
                (header_ptr->n_tick_count > *max_latency) ? header_ptr->n_tick_count : *max_latency;

            app_svt_av1_get_time(&finish_s_time, &finish_u_time);

            // total execution time, inc init time
            config->performance_context.total_execution_time =
                app_svt_av1_compute_overall_elapsed_time(
                    config->performance_context.lib_start_time[0],
                    config->performance_context.lib_start_time[1],
                    finish_s_time,
                    finish_u_time);

            // total encode time
            config->performance_context.total_encode_time =
                app_svt_av1_compute_overall_elapsed_time(
                    config->performance_context.encode_start_time[0],
                    config->performance_context.encode_start_time[1],
                    finish_s_time,
                    finish_u_time);

            // Write Stream Data to file
            if (stream_file) {
                if (config->performance_context.frame_count == 1 &&
                    !(flags & EB_BUFFERFLAG_IS_ALT_REF)) {
                    write_ivf_stream_header(config);
                }
                write_ivf_frame_header(config, header_ptr->n_filled_len);
                fwrite(header_ptr->p_buffer, 1, header_ptr->n_filled_len, stream_file);
            }

            config->performance_context.byte_count += header_ptr->n_filled_len;

            if (config->config.stat_report && !(flags & EB_BUFFERFLAG_IS_ALT_REF))
                process_output_statistics_buffer(header_ptr, config);

            // Update Output Port Activity State
            *port_state  = (flags & EB_BUFFERFLAG_EOS) ? APP_PortInactive : *port_state;
            return_value = (flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished
                                                                   : APP_ExitConditionNone;
            // Release the output buffer
            svt_av1_enc_release_out_buffer(&header_ptr);

            if (flags & EB_BUFFERFLAG_EOS) {
                if (config->config.rc_firstpass_stats_out) {
                    SvtAv1FixedBuf first_pass_stat;
                    EbErrorType ret = svt_av1_enc_get_stream_info(component_handle,
                        SVT_AV1_STREAM_INFO_FIRST_PASS_STATS_OUT, &first_pass_stat);
                    if (ret == EB_ErrorNone) {
                        if (config->output_stat_file) {
                            fwrite(first_pass_stat.buf,
                                1, first_pass_stat.sz, config->output_stat_file);
                        }
                        enc_app->rc_twopasses_stats.buf = realloc(enc_app->rc_twopasses_stats.buf, first_pass_stat.sz);
                        if (enc_app->rc_twopasses_stats.buf) {
                            memcpy(enc_app->rc_twopasses_stats.buf, first_pass_stat.buf, first_pass_stat.sz);
                            enc_app->rc_twopasses_stats.sz = first_pass_stat.sz;
                        }
                    }
                }

            }

            ++*frame_count;
            const double fps = (double)*frame_count / config->performance_context.total_encode_time;
            const double frame_rate = config->config.frame_rate_numerator && config->config.frame_rate_denominator
                ? (double)config->config.frame_rate_numerator / (double)config->config.frame_rate_denominator
                : config->config.frame_rate > 1000
                    // Correct for 16-bit fixed-point fractional precision
                    ? (double)config->config.frame_rate / (1 << 16)
                    : (double)config->config.frame_rate;
            switch (config->progress) {
            case 0: break;
            case 1:
                if (!(flags & EB_BUFFERFLAG_IS_ALT_REF))
                    fprintf(stderr, "\b\b\b\b\b\b\b\b\b%9d", *frame_count);
                break;
            case 2:
                fprintf(stderr,
                        "\rEncoding frame %4d %.2f kbps %.2f fp%c  ",
                        *frame_count,
                        ((double)(config->performance_context.byte_count << 3) * frame_rate /
                         (config->frames_encoded * 1000)),
                        fps >= 1.0 ? fps : fps * 60,
                        fps >= 1.0 ? 's' : 'm');
            default: break;
            }
            fflush(stderr);

            config->performance_context.average_speed =
                (double)config->performance_context.frame_count /
                config->performance_context.total_encode_time;
            config->performance_context.average_latency =
                (double)config->performance_context.total_latency /
                config->performance_context.frame_count;

            if (config->progress == 1 && !(*frame_count % SPEED_MEASUREMENT_INTERVAL))
                fprintf(stderr,
                        "\nAverage System Encoding Speed:        %.2f\n",
                        (double)*frame_count / config->performance_context.total_encode_time);
        }
    }
    channel->exit_cond_output = return_value;
}
void process_output_recon_buffer(EncChannel* channel) {
    EbConfig *config = channel->config;
    EbAppContext *app_call_back = channel->app_callback;
    EbBufferHeaderType *header_ptr =
        app_call_back->recon_buffer; // needs to change for buffered input
    EbComponentType *    component_handle = (EbComponentType *)app_call_back->svt_encoder_handle;
    AppExitConditionType return_value     = APP_ExitConditionNone;
    int32_t              fseek_return_val;
    if (channel->exit_cond_recon != APP_ExitConditionNone) {
        return;
    }
    // non-blocking call until all input frames are sent
    EbErrorType recon_status = svt_av1_get_recon(component_handle, header_ptr);

    if (recon_status == EB_ErrorMax) {
        fprintf(stderr, "\n");
        log_error_output(config->error_log_file, header_ptr->flags);
        channel->exit_cond_recon = APP_ExitConditionError;
        return;
    } else if (recon_status != EB_NoErrorEmptyQueue) {
        //Sets the File position to the beginning of the file.
        rewind(config->recon_file);
        uint64_t frame_num = header_ptr->pts;
        while (frame_num > 0) {
            fseek_return_val = fseeko(config->recon_file, header_ptr->n_filled_len, SEEK_CUR);

            if (fseek_return_val != 0) {
                fprintf(stderr, "Error in fseeko  returnVal %i\n", fseek_return_val);
                channel->exit_cond_recon = APP_ExitConditionError;
                return;
            }
            frame_num = frame_num - 1;
        }

        fwrite(header_ptr->p_buffer, 1, header_ptr->n_filled_len, config->recon_file);

        // Update Output Port Activity State
        return_value = (header_ptr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished
                                                               : APP_ExitConditionNone;
    }
    channel->exit_cond_recon =  return_value;
}
