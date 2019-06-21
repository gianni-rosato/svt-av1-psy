/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "EbAppContext.h"
#include "EbAppConfig.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbAppInputy4m.h"

#include "EbTime.h"

#define IVF_FRAME_HEADER_IN_LIB                     0

/***************************************
 * Macros
 ***************************************/
#define CLIP3(min_val, max_val, a)        (((a)<(min_val)) ? (min_val) : (((a)>(max_val)) ? (max_val) :(a)))
#define FUTURE_WINDOW_WIDTH                 4
#define SIZE_OF_ONE_FRAME_IN_BYTES(width, height, csp, is16bit) \
    ( (((width)*(height)) + 2*(((width)*(height))>>(3-csp)))<<is16bit)
#define YUV4MPEG2_IND_SIZE 9
extern volatile int32_t keepRunning;

/***************************************
* Process Error Log
***************************************/
void LogErrorOutput(
    FILE                     *error_log_file,
    uint32_t                  error_code)
{
    switch (error_code) {
        // EB_ENC_AMVP_ERRORS:
    case EB_ENC_AMVP_ERROR1:
        fprintf(error_log_file, "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR2:
        fprintf(error_log_file, "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

    case EB_ENC_AMVP_ERROR3:
        fprintf(error_log_file, "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR4:
        fprintf(error_log_file, "Error: The availability parameter in GetSpatialMVPPosAx() function can not be > 3 !\n");
        break;

    case EB_ENC_AMVP_ERROR5:
        fprintf(error_log_file, "Error: The availability parameter in GetSpatialMVPPosBx() function can not be > 7 !\n");
        break;

    case EB_ENC_AMVP_ERROR6:
        fprintf(error_log_file, "Error: GetTemporalMVP: tmvpMapLcuIndex must be either 0 or 1!\n");
        break;

    case EB_ENC_AMVP_ERROR7:
        fprintf(error_log_file, "Error: the input PU to GetNonScalingSpatialAMVP() must be available!");
        break;

    case EB_ENC_AMVP_ERROR8:
        fprintf(error_log_file, "Error: GetTemporalMVP: tmvpMapLcuIndex must be either 0 or 1");
        break;

    case EB_ENC_AMVP_NULL_REF_ERROR:
        fprintf(error_log_file, "Error: The referenceObject can not be NULL!\n");
        break;

    case EB_ENC_AMVP_SPATIAL_NA_ERROR:
        fprintf(error_log_file, "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

        // EB_ENC_CL_ERRORS:
    case EB_ENC_CL_ERROR1:
        fprintf(error_log_file, "Error: Unknown Inter Prediction Direction!\n");
        break;

    case EB_ENC_CL_ERROR2:
        fprintf(error_log_file, "Error: Unknown coding mode!\n");
        break;

    case EB_ENC_CL_ERROR3:
        fprintf(error_log_file, "Error: Mode Decision Candidate Buffer Overflow!\n");
        break;

    case EB_ENC_CL_ERROR4:
        fprintf(error_log_file, "Error: Too many Mode Decision Fast Candidates!\n");
        break;

    case EB_ENC_CL_ERROR5:
        fprintf(error_log_file, "Error: Too many buffers chosen for this level by PreModeDecision!\n");
        break;

    case EB_ENC_CL_ERROR6:
        fprintf(error_log_file, "Error: Ping-Pong structure needs at least two buffers to work properly!\n");
        break;

    case EB_ENC_CL_ERROR7:
        fprintf(error_log_file, "Error: Invalid Intra Partition\n");
        break;

    case EB_ENC_CL_ERROR8:
        fprintf(error_log_file, "Error: Invalid TU Configuration\n");
        break;
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
        fprintf(error_log_file, "Error: The x/y location of the CU must be the multiple of minmum CU size!");
        break;

    case EB_ENC_DLF_ERROR5:
        fprintf(error_log_file, "Error: Unknown Slice Type!");
        break;

    case EB_ENC_DLF_ERROR6:
        fprintf(error_log_file, "Error: While calculating the bS for the PU bounday, the 4x4 block must be guaranteed at the PU boundary!");
        break;

    case EB_ENC_DLF_ERROR7:
        fprintf(error_log_file, "Error: LCU size must be power of 2!");
        break;

    case EB_ENC_DLF_ERROR8:
        fprintf(error_log_file, "Error: Deblocking filter can not support the picture whose width or height is not the multiple of 8!");
        break;

    case EB_ENC_DLF_ERROR9:
        fprintf(error_log_file, "Error: Neighbor PU must be available!");
        break;

    case EB_ENC_DLF_ERROR10:
        fprintf(error_log_file, "Error: Deblocking filter can not support the picture whose width or height is not the multiple of 8!");
        break;

        // EB_ENC_EC_ERRORS:
    case EB_ENC_EC_ERROR1:
        fprintf(error_log_file, "Error: EncodeCodedBlockFlags: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR10:
        fprintf(error_log_file, "Error: EncodeTuSplitCoeff: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR11:
        fprintf(error_log_file, "Error: CodeSPS: Long term reference pictures are not currently handled!\n");
        break;

    case EB_ENC_EC_ERROR12:
        fprintf(error_log_file, "Error: CodeProfileTierLevel: The maximum sublayers must be equal to 1!\n");
        break;

    case EB_ENC_EC_ERROR13:
        fprintf(error_log_file, "Error: EncodeRootCodedBlockFlag: rootCbf too large!\n");
        break;

    case EB_ENC_EC_ERROR14:
        fprintf(error_log_file, "Error: cpbCountMinus1 in HRD parameter exceeds the upper limit 4!\n");
        break;

    case EB_ENC_EC_ERROR15:
        fprintf(error_log_file, "Error: numDecodingUnitsMinus1 in picture timeing SEI exceeds the upper limit 64!\n");
        break;

    case EB_ENC_EC_ERROR16:
        fprintf(error_log_file, "Error: The size of the unregistered user data SEI payload is not allowed!\n");
        break;

    case EB_ENC_EC_ERROR2:
        fprintf(error_log_file, "Error: copy_rbsp_bitstream_to_payload: output buffer too small!\n");
        break;

    case EB_ENC_EC_ERROR3:
        fprintf(error_log_file, "Error: EncodeLcu: Unknown mode type!\n");
        break;

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

    case EB_ENC_EC_ERROR26:
        fprintf(error_log_file, "Error: Level not recognized!\n");
        break;

    case EB_ENC_EC_ERROR27:
        fprintf(error_log_file, "Error: EncodeOneBin:  BinaryValue must be less than 2\n");
        break;

    case EB_ENC_EC_ERROR28:
        fprintf(error_log_file, "Error: No more than 6 SAO types\n");
        break;

    case EB_ENC_EC_ERROR29:
        fprintf(error_log_file, "Error: No more than 6 SAO types\n");
        break;

        // EB_ENC_FL_ERRORS:
    case EB_ENC_FL_ERROR1:
        fprintf(error_log_file, "Error: Uncovered area inside Cu!\n");
        break;

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
        fprintf(error_log_file, "Error: RC Results Fifo Size should be greater than RC Tasks Fifo Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR13:
        fprintf(error_log_file, "Error: RC Tasks Fifo Size should be greater than EC results Fifo Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR14:
        fprintf(error_log_file, "Error: RC Tasks Fifo Size should be greater than Picture Manager results Fifo Size in order to avoid deadlock!\n");
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
        fprintf(error_log_file, "Error: Motion compensation prediction is out of the picture boundary!\n");
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
        fprintf(error_log_file, "Error: IntraPrediction: intra prediction only supports square PU!\n");
        break;

    case EB_ENC_INTRA_PRED_ERROR3:
        fprintf(error_log_file, "Error: IntraPredictionChroma: Only Planar!\n");
        break;

    case EB_ENC_INVLD_PART_SIZE_ERROR:
        fprintf(error_log_file, "Error: IntraPrediction: only PU sizes of 8 or largers are currently supported!\n");
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

    case EB_ENC_MD_ERROR4:
        fprintf(error_log_file, "Error: Unknown ME SAD Level!\n");
        break;

    case EB_ENC_MD_ERROR5:
        fprintf(error_log_file, "Error: Invalid encoder mode. The encoder mode should be 0, 1 or 2!\n");
        break;

    case EB_ENC_MD_ERROR6:
        fprintf(error_log_file, "Error: Invalid TU size!\n");
        break;

    case EB_ENC_MD_ERROR7:
        fprintf(error_log_file, "Error: Unknown depth!\n");
        break;

    case EB_ENC_MD_ERROR8:
        fprintf(error_log_file, "Error: Depth not supported!\n");
        break;

    case EB_ENC_MD_ERROR9:
        fprintf(error_log_file, "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

    case EB_ENC_MD_ERROR10:
        fprintf(error_log_file, "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

        // EB_ENC_ME_ERRORS:
    case EB_ENC_ME_ERROR1:
        fprintf(error_log_file, "Error: Motion Estimation: non valid value of the subPelDirection !\n");
        break;

    case EB_ENC_ME_ERROR2:
        fprintf(error_log_file, "Error: FillMvMergeCandidate() method only supports P or B slices!\n");
        break;

        // EB_ENC_ERRORS:
    case EB_ENC_ROB_OF_ERROR:
        fprintf(error_log_file, "Error: Recon Output Buffer Overflow!\n");
        break;

        // EB_ENC_PACKETIZATION_ERRORS:
    case EB_ENC_PACKETIZATION_ERROR1:
        fprintf(error_log_file, "Error: PacketizationProcess: Picture Number does not match entry. PacketizationReorderQueue overflow!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR2:
        fprintf(error_log_file, "Error: Entropy Coding Result can not be outputed by processes other than entropy coder and ALF!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR3:
        fprintf(error_log_file, "Error: The encoder can not support the SliceMode other than 0 and 1!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR4:
        fprintf(error_log_file, "Error: Statistics Output Buffer Overflow!\n");
        break;
    case EB_ENC_PACKETIZATION_ERROR5:
        fprintf(error_log_file, "Error: Stream Fifo is starving..deadlock, increase EB_outputStreamBufferFifoInitCount APP_ENCODERSTREAMBUFFERCOUNT \n");
        break;

        // EB_ENC_PM_ERRORS:
    case EB_ENC_PM_ERROR0:
        fprintf(error_log_file, "Error: PictureManagerProcess: Unknown Slice Type!\n");
        break;

    case EB_ENC_PM_ERROR1:
        fprintf(error_log_file, "Error: EbPictureManager: dependent_count underflow!\n");
        break;

    case EB_ENC_PM_ERROR10:
        fprintf(error_log_file, "Error: picture_manager_kernel: referenceEntryPtr should never be null!\n");
        break;

    case EB_ENC_PM_ERROR2:
        fprintf(error_log_file, "Error: PictureManagerProcess: The Reference Structure period must be less than the MAX_ELAPSED_IDR_COUNT or false-IDR boundary logic will be activated!\n");
        break;

    case EB_ENC_PM_ERROR3:
        fprintf(error_log_file, "Error: PictureManagerProcess: The dependent_count underflow detected!\n");
        break;

    case EB_ENC_PM_ERROR4:
        fprintf(error_log_file, "Error: PictureManagerProcess: Empty input queue!\n");
        break;

    case EB_ENC_PM_ERROR5:
        fprintf(error_log_file, "Error: PictureManagerProcess: Empty reference queue!\n");
        break;

    case EB_ENC_PM_ERROR6:
        fprintf(error_log_file, "Error: PictureManagerProcess: The capped elaspedNonIdrCount must be larger than the maximum supported delta ref poc!\n");
        break;

    case EB_ENC_PM_ERROR7:
        fprintf(error_log_file, "Error: PictureManagerProcess: Reference Picture Queue Full!\n");
        break;

    case EB_ENC_PM_ERROR8:
        fprintf(error_log_file, "Error: PictureManagerProcess: No reference match found - this will lead to a memory leak!\n");
        break;

    case EB_ENC_PM_ERROR9:
        fprintf(error_log_file, "Error: PictureManagerProcess: Unknown picture type!\n");
        break;

    case EB_ENC_PM_ERROR12:
        fprintf(error_log_file, "Error: PictureManagerProcess: prediction structure configuration API has too many reference pictures\n");
        break;

    case EB_ENC_PM_ERROR13:
        fprintf(error_log_file, "Error: PictureManagerProcess: The maximum allowed frame rate is 60 fps\n");
        break;

    case EB_ENC_PM_ERROR14:
        fprintf(error_log_file, "Error: PictureManagerProcess: The minimum allowed frame rate is 1 fps\n");
        break;

        // EB_ENC_PRED_STRC_ERRORS:
    case EB_ENC_PRED_STRC_ERROR1:
        fprintf(error_log_file, "Error: PredictionStructureCtor: DecodeOrder LUT too small!\n");
        break;

    case EB_ENC_PRED_STRC_ERROR2:
        fprintf(error_log_file, "Error: PredictionStructureCtor: prediction structure improperly configured!\n");
        break;

        // EB_ENC_PU_ERRORS:
    case EB_ENC_PU_ERROR1:
        fprintf(error_log_file, "Error: Unknown partition size!\n");
        break;

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
        fprintf(error_log_file, "Error: RateControlProcess: No feedback frame match found - this will lead to a memory leak!\n");
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
        fprintf(error_log_file, "Error: intra2_nx2_n_fast_cost_islice can only support 2Nx2N partition type!\n");
        break;

        // EB_ENC_SAO_ERRORS:
    case EB_ENC_SAO_ERROR1:
        fprintf(error_log_file, "Error: No more than 6 SAO types!\n");
        break;

    case EB_ENC_SAO_ERROR2:
        fprintf(error_log_file, "Error: No more than 5 EO SAO categories!\n");
        break;
        // EB_ENC_SCS_ERRORS:
    case EB_ENC_SCS_ERROR1:
        fprintf(error_log_file, "Error: SequenceControlSetCopy: Not all SequenceControlSet members are being copied!\n");
        break;

        // EB_ENC_BITSTREAM_ERRORS:
    case EB_ENC_BITSTREAM_ERROR1:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Bitstream payload buffer empty!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR2:
        fprintf(error_log_file, "Error: OutputBitstreamWrite: Empty bitstream!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR3:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Buffer index more than buffer size!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR4:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Start Location in not inside the buffer!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR5:
        fprintf(error_log_file, "Error: OutputBitstreamWrite: Trying to write more than one word!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR6:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Expecting Start code!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR7:
        fprintf(error_log_file, "Error: OutputBitstreamRBSPToPayload: Bitstream not flushed (i.e. byte-aligned)!\n");
        break;

    case EB_ENC_RESS_COOR_ERRORS1:
        fprintf(error_log_file, "Error: ResourceCoordinationProcess: The received input data should be equal to the buffer size - only complete frame transfer is supported\n");
        break;

    case EB_ENC_RES_COORD_InvalidQP:
        fprintf(error_log_file, "Error: ResourceCoordinationProcess: The QP value in the QP file is invalid\n");
        break;

    case EB_ENC_RES_COORD_InvalidSliceType:
        fprintf(error_log_file, "Error: ResourceCoordinationProcess: Slice Type Invalid\n");
        break;

        // picture decision Errors
    case EB_ENC_PD_ERROR8:
        fprintf(error_log_file, "Error: PictureDecisionProcess: Picture Decision Reorder Queue overflow\n");
        break;

    default:
        fprintf(error_log_file, "Error: Others!\n");
        break;
    }

    return;
}

/******************************************************
* Copy fields from the stream to the input buffer
    Input   : stream
    Output  : valid input buffer
******************************************************/
void ProcessInputFieldStandardMode(
    EbConfig               *config,
    EbBufferHeaderType      *headerPtr,
    FILE                     *input_file,
    uint8_t                    *lumaInputPtr,
    uint8_t                    *cbInputPtr,
    uint8_t                    *crInputPtr,
    uint8_t                   is16bit) {
    const int64_t input_padded_width  = config->input_padded_width;
    const int64_t input_padded_height = config->input_padded_height;
    const uint8_t color_format = config->encoder_color_format;
    const uint8_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint8_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;
    const uint64_t source_luma_row_size = (uint64_t)(input_padded_width << is16bit);
    const uint64_t source_chroma_row_size = source_luma_row_size >> subsampling_x;
    uint8_t  *ebInputPtr;
    uint32_t  inputRowIndex;

    // Y
    ebInputPtr = lumaInputPtr;
    // Skip 1 luma row if bottom field (point to the bottom field)
    if (config->processed_frame_count % 2 != 0)
        fseeko64(input_file, (long)source_luma_row_size, SEEK_CUR);

    for (inputRowIndex = 0; inputRowIndex < input_padded_height; inputRowIndex++) {
        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, source_luma_row_size, input_file);
        // Skip 1 luma row (only fields)
        fseeko64(input_file, (long)source_luma_row_size, SEEK_CUR);
        ebInputPtr += source_luma_row_size;
    }

    // U
    ebInputPtr = cbInputPtr;
    // Step back 1 luma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    if (config->processed_frame_count % 2 != 0) {
        fseeko64(input_file, -(long)source_luma_row_size, SEEK_CUR);
        fseeko64(input_file, (long)source_chroma_row_size, SEEK_CUR);
    }

    for (inputRowIndex = 0; inputRowIndex < input_padded_height >> subsampling_y; inputRowIndex++) {
        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, source_chroma_row_size, input_file);
        // Skip 1 chroma row (only fields)
        fseeko64(input_file, (long)source_chroma_row_size, SEEK_CUR);
        ebInputPtr += source_chroma_row_size;
    }

    // V
    ebInputPtr = crInputPtr;
    // Step back 1 chroma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    // => no action

    for (inputRowIndex = 0; inputRowIndex < input_padded_height >> subsampling_y; inputRowIndex++) {
        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, source_chroma_row_size, input_file);
        // Skip 1 chroma row (only fields)
        fseeko64(input_file, (long)source_chroma_row_size, SEEK_CUR);
        ebInputPtr += source_chroma_row_size;
    }

    // Step back 1 chroma row if bottom field (undo the previous jump)
    if (config->processed_frame_count % 2 != 0)
        fseeko64(input_file, -(long)source_chroma_row_size, SEEK_CUR);
}

//************************************/
// GetNextQpFromQpFile
// Reads and extracts one qp from the qp_file
// Input  : QP file
// Output : QP value
/************************************/
static int32_t qpReadFromFile = 0;

int32_t GetNextQpFromQpFile(
    EbConfig  *config
)
{
    uint8_t *line;
    int32_t qp = 0;
    uint32_t readsize = 0, eof = 0;
    EB_APP_MALLOC(uint8_t*, line, 8, EB_N_PTR, EB_ErrorInsufficientResources);
    memset(line,0,8);
    readsize = (uint32_t)fread(line, 1, 2, config->qp_file);

    if (readsize == 0) {
        // end of file
        return -1;
    }
    else if (readsize == 1) {
        qp = strtol((const char*)line, NULL, 0);
        if (qp == 0) // eof
            qp = -1;
    }
    else if (readsize == 2 && (line[0] == '\n')) {
        // new line
        fseek(config->qp_file, -1, SEEK_CUR);
        qp = 0;
    }
    else if (readsize == 2 && (line[1] == '\n')) {
        // new line
        qp = strtol((const char*)line, NULL, 0);
    }
    else if (readsize == 2 && (line[0] == '#' || line[0] == '/' || line[0] == '-' || line[0] == ' ')) {
        // Backup one step to not miss the new line char
        fseek(config->qp_file, -1, SEEK_CUR);
        do {
            readsize = (uint32_t)fread(line, 1, 1, config->qp_file);
            if (readsize != 1)
                break;
        } while (line[0] != '\n');

        if (eof != 0)
            // end of file
            qp = -1;
        else
            // skip line
            qp = 0;
    }
    else if (readsize == 2) {
        qp = strtol((const char*)line, NULL, 0);
        do {
            readsize = (uint32_t)fread(line, 1, 1, config->qp_file);
            if (readsize != 1)
                break;
        } while (line[0] != '\n');
    }

    if (qp > 0)
        qpReadFromFile |= 1;

    return qp;
}

void ReadInputFrames(
    EbConfig                  *config,
    uint8_t                      is16bit,
    EbBufferHeaderType         *headerPtr)
{
    uint64_t  readSize;
    const uint32_t  input_padded_width = config->input_padded_width;
    const uint32_t  input_padded_height = config->input_padded_height;
    FILE   *input_file = config->input_file;
    uint8_t  *ebInputPtr;
    EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;

    const uint8_t color_format = config->encoder_color_format;
    const uint8_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;

    inputPtr->y_stride  = input_padded_width;
    inputPtr->cr_stride = input_padded_width >> subsampling_x;
    inputPtr->cb_stride = input_padded_width >> subsampling_x;

    if (config->buffered_input == -1) {
        if (is16bit == 0 || (is16bit == 1 && config->compressed_ten_bit_format == 0)) {
            readSize = (uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(input_padded_width, input_padded_height, color_format, is16bit);

            headerPtr->n_filled_len = 0;

            // Interlaced Video
            if (config->separate_fields) {
                ProcessInputFieldStandardMode(
                    config,
                    headerPtr,
                    input_file,
                    inputPtr->luma,
                    inputPtr->cb,
                    inputPtr->cr,
                    is16bit);

                if (readSize != headerPtr->n_filled_len) {
                    fseek(input_file, 0, SEEK_SET);
                    headerPtr->n_filled_len = 0;

                    ProcessInputFieldStandardMode(
                        config,
                        headerPtr,
                        input_file,
                        inputPtr->luma,
                        inputPtr->cb,
                        inputPtr->cr,
                        is16bit);
                }

                // Reset the pointer position after a top field
                if (config->processed_frame_count % 2 == 0)
                    fseek(input_file, -(long)(readSize << 1), SEEK_CUR);
            }
            else {
                /* if input is a y4m file, read next line which contains "FRAME" */
                if(config->y4m_input==EB_TRUE)
                    read_y4m_frame_delimiter(config);
                uint64_t lumaReadSize = (uint64_t)input_padded_width*input_padded_height << is16bit;
                ebInputPtr = inputPtr->luma;
                if(config->y4m_input==EB_FALSE && config->processed_frame_count == 0 && config->input_file == stdin) {
                    /* if not a y4m file and input is read from stdin, 9 bytes were already read when checking
                        or the YUV4MPEG2 string in the stream, so copy those bytes over */
                    memcpy(ebInputPtr,config->y4m_buf,YUV4MPEG2_IND_SIZE);
                    headerPtr->n_filled_len += YUV4MPEG2_IND_SIZE;
                    ebInputPtr += YUV4MPEG2_IND_SIZE;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize-YUV4MPEG2_IND_SIZE, input_file);
                }else {
                    headerPtr->n_filled_len += (uint32_t)fread(inputPtr->luma, 1, lumaReadSize, input_file);
                }
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb, 1, lumaReadSize >> (3 - color_format), input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr, 1, lumaReadSize >> (3 - color_format), input_file);

                if (readSize != headerPtr->n_filled_len) {
                    fseek(input_file, 0, SEEK_SET);
                    headerPtr->n_filled_len = (uint32_t)fread(inputPtr->luma, 1, lumaReadSize, input_file);
                    headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb, 1, lumaReadSize >> (3 - color_format), input_file);
                    headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr, 1, lumaReadSize >> (3 - color_format), input_file);
                }
            }
        } else if (is16bit == 1 && config->compressed_ten_bit_format == 1) {
            // 10-bit Compressed Unpacked Mode
            const uint32_t lumaReadSize = input_padded_width * input_padded_height;
            const uint32_t chromaReadSize = lumaReadSize >> (3 - color_format);
            const uint32_t nbitLumaReadSize = (input_padded_width / 4) * input_padded_height;
            const uint32_t nbitChromaReadSize = nbitLumaReadSize >> (3 - color_format);

            // Fill the buffer with a complete frame
            headerPtr->n_filled_len = 0;

            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->luma, 1, lumaReadSize, input_file);
            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb, 1, chromaReadSize, input_file);
            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr, 1, chromaReadSize, input_file);

            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->luma_ext, 1, nbitLumaReadSize, input_file);
            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb_ext, 1, nbitChromaReadSize, input_file);
            headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr_ext, 1, nbitChromaReadSize, input_file);

            readSize = lumaReadSize + nbitLumaReadSize + 2 * (chromaReadSize + nbitChromaReadSize);

            if (readSize != headerPtr->n_filled_len) {
                fseek(input_file, 0, SEEK_SET);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->luma, 1, lumaReadSize, input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb, 1, chromaReadSize, input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr, 1, chromaReadSize, input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->luma_ext, 1, nbitLumaReadSize, input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cb_ext, 1, nbitChromaReadSize, input_file);
                headerPtr->n_filled_len += (uint32_t)fread(inputPtr->cr_ext, 1, nbitChromaReadSize, input_file);
            }
        }
    } else {
        if (is16bit && config->compressed_ten_bit_format == 1) {
            // Determine size of each plane
            const size_t luma8bitSize = input_padded_width * input_padded_height;
            const size_t chroma8bitSize = luma8bitSize >> (3 - color_format);
            const size_t luma2bitSize = luma8bitSize / 4; //4-2bit pixels into 1 byte
            const size_t chroma2bitSize = luma2bitSize >> (3 - color_format);

            EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;
            inputPtr->y_stride = input_padded_width;
            inputPtr->cr_stride = input_padded_width >> subsampling_x;
            inputPtr->cb_stride = input_padded_width >> subsampling_x;

            inputPtr->luma = config->sequence_buffer[config->processed_frame_count % config->buffered_input];
            inputPtr->cb = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + luma8bitSize;
            inputPtr->cr = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + luma8bitSize + chroma8bitSize;

            inputPtr->luma_ext = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + luma8bitSize + 2 * chroma8bitSize;
            inputPtr->cb_ext = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + luma8bitSize + 2 * chroma8bitSize + luma2bitSize;
            inputPtr->cr_ext = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + luma8bitSize + 2 * chroma8bitSize + luma2bitSize + chroma2bitSize;

            headerPtr->n_filled_len = (uint32_t)(luma8bitSize + luma2bitSize + 2 * (chroma8bitSize + chroma2bitSize));
        } else {
            //Normal unpacked mode:yuv420p10le yuv422p10le yuv444p10le
            const size_t lumaSize = (input_padded_width * input_padded_height) << is16bit;
            const size_t chromaSize = lumaSize >> (3 - color_format);

            EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;

            inputPtr->y_stride = input_padded_width;
            inputPtr->cr_stride = input_padded_width >> subsampling_x;
            inputPtr->cb_stride = input_padded_width >> subsampling_x;

            inputPtr->luma = config->sequence_buffer[config->processed_frame_count % config->buffered_input];
            inputPtr->cb = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + lumaSize;
            inputPtr->cr = config->sequence_buffer[config->processed_frame_count % config->buffered_input] + lumaSize+ chromaSize;

            headerPtr->n_filled_len = (uint32_t)(lumaSize + 2 * chromaSize);
        }
    }

    // If we reached the end of file, loop over again
    if (feof(input_file) != 0)
        fseek(input_file, 0, SEEK_SET);
    return;
}

void SendQpOnTheFly(
    EbConfig                  *config,
    EbBufferHeaderType        *headerPtr){
    {
        uint32_t           qpPtr;
        int32_t           tmpQp = 0;

        do {
            // get next qp
            tmpQp = GetNextQpFromQpFile(config);

            if (tmpQp == (int32_t)EB_ErrorInsufficientResources) {
                printf("Malloc has failed due to insuffucient resources");
                return;
            }

            // check if eof
            if ((tmpQp == -1) && (qpReadFromFile != 0))
                fseek(config->qp_file, 0, SEEK_SET);

            // check if the qp read is valid
            else if (tmpQp > 0)
                break;
        } while (tmpQp == 0 || ((tmpQp == -1) && (qpReadFromFile != 0)));

        if (tmpQp == -1) {
            config->use_qp_file = EB_FALSE;
            printf("\nWarning: QP File did not contain any valid QPs");
        }

        qpPtr = CLIP3(0, 51, tmpQp);

        headerPtr->qp = qpPtr;
    }
    return;
}

//************************************/
// ProcessInputBuffer
// Reads yuv frames from file and copy
// them into the input buffer
/************************************/
AppExitConditionType ProcessInputBuffer(
    EbConfig             *config,
    EbAppContext         *appCallBack)
{
    uint8_t                 is16bit = (uint8_t)(config->encoder_bit_depth > 8);
    EbBufferHeaderType     *headerPtr = appCallBack->input_buffer_pool;
    EbComponentType        *componentHandle = (EbComponentType*)appCallBack->svt_encoder_handle;

    AppExitConditionType    return_value = APP_ExitConditionNone;

    const uint8_t color_format = config->encoder_color_format;
    const int64_t input_padded_width = config->input_padded_width;
    const int64_t input_padded_height = config->input_padded_height;
    const int64_t frames_to_be_encoded = config->frames_to_be_encoded;
    int64_t totalBytesToProcessCount;
    int64_t remainingByteCount;
    uint32_t compressed10bitFrameSize = (uint32_t)((input_padded_width*input_padded_height) + 2 * ((input_padded_width*input_padded_width) >> (3 - color_format)));
    compressed10bitFrameSize += compressed10bitFrameSize / 4;

    if (config->injector && config->processed_frame_count)
        Injector(config->processed_frame_count, config->injector_frame_rate);
    totalBytesToProcessCount = (frames_to_be_encoded < 0) ? -1 : (config->encoder_bit_depth == 10 && config->compressed_ten_bit_format == 1) ?
        frames_to_be_encoded * (int64_t)compressed10bitFrameSize:
        frames_to_be_encoded * SIZE_OF_ONE_FRAME_IN_BYTES(input_padded_width, input_padded_height, color_format, is16bit);

    remainingByteCount       = (totalBytesToProcessCount < 0) ?   -1 :  totalBytesToProcessCount - (int64_t)config->processed_byte_count;

    // If there are bytes left to encode, configure the header
    if (remainingByteCount != 0 && config->stop_encoder == EB_FALSE) {
        ReadInputFrames(
            config,
            is16bit,
            headerPtr);

        // Update the context parameters
        config->processed_byte_count += headerPtr->n_filled_len;
        headerPtr->p_app_private          = (EbPtr)EB_NULL;
        config->frames_encoded           = (int32_t)(++config->processed_frame_count);

        // Configuration parameters changed on the fly
        if (config->use_qp_file && config->qp_file)
            SendQpOnTheFly(
                config,
                headerPtr);

        if (keepRunning == 0 && !config->stop_encoder)
            config->stop_encoder = EB_TRUE;
        // Fill in Buffers Header control data
        headerPtr->pts          = config->processed_frame_count-1;
        headerPtr->pic_type    = EB_AV1_INVALID_PICTURE;

        headerPtr->flags = 0;

        // Send the picture
        eb_svt_enc_send_picture(componentHandle, headerPtr);

        if ((config->processed_frame_count == (uint64_t)config->frames_to_be_encoded) || config->stop_encoder) {
            headerPtr->n_alloc_len    = 0;
            headerPtr->n_filled_len   = 0;
            headerPtr->n_tick_count   = 0;
            headerPtr->p_app_private  = NULL;
            headerPtr->flags       = EB_BUFFERFLAG_EOS;
            headerPtr->p_buffer      = NULL;
            headerPtr->pic_type    = EB_AV1_INVALID_PICTURE;

            eb_svt_enc_send_picture(componentHandle, headerPtr);
        }

        return_value = (headerPtr->flags == EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : return_value;
    }

    return return_value;
}

#define LONG_ENCODE_FRAME_ENCODE    4000
#define SPEED_MEASUREMENT_INTERVAL  2000
#define START_STEADY_STATE          1000
#define AV1_FOURCC                  0x31305641 // used for ivf header
#define IVF_STREAM_HEADER_SIZE      32
#define IVF_FRAME_HEADER_SIZE       12
#define OBU_FRAME_HEADER_SIZE       3
#define TD_SIZE                     2
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

static void write_ivf_stream_header(EbConfig *config)
{
    char header[IVF_STREAM_HEADER_SIZE];
    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4, 0);                     // version
    mem_put_le16(header + 6, 32);                    // header size
    mem_put_le32(header + 8, AV1_FOURCC);                // fourcc
    mem_put_le16(header + 12, config->input_padded_width);  // width
    mem_put_le16(header + 14, config->input_padded_height); // height
    if (config->frame_rate_denominator != 0 && config->frame_rate_numerator != 0){
        mem_put_le32(header + 16, config->frame_rate_numerator);  // rate
        mem_put_le32(header + 20, config->frame_rate_denominator);            // scale
                                                    //mem_put_le32(header + 16, config->frame_rate_denominator);  // rate
                                                    //mem_put_le32(header + 20, config->frame_rate_numerator);  // scale
    }
    else {
        mem_put_le32(header + 16, (config->frame_rate >> 16) * 1000);  // rate
        mem_put_le32(header + 20, 1000);            // scale
                                                    //mem_put_le32(header + 16, config->frame_rate_denominator);  // rate
                                                    //mem_put_le32(header + 20, config->frame_rate_numerator);  // scale
    }
    mem_put_le32(header + 24, 0);               // length
    mem_put_le32(header + 28, 0);               // unused
    //config->performance_context.byte_count += 32;
    if (config->bitstream_file)
        fwrite(header, 1, IVF_STREAM_HEADER_SIZE, config->bitstream_file);

    return;
}
static void update_prev_ivf_header(EbConfig *config){
    char header[4]; // only for the number of bytes
    if (config && config->bitstream_file && config->byte_count_since_ivf != 0){
        fseeko64(config->bitstream_file, (-(int32_t)(config->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE)),SEEK_CUR);
        mem_put_le32(&header[0], (int32_t)(config->byte_count_since_ivf));
        fwrite(header, 1, 4, config->bitstream_file);
        fseeko64(config->bitstream_file, (config->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE - 4), SEEK_CUR);
        config->byte_count_since_ivf = 0;
    }
}

static void write_ivf_frame_header(EbConfig *config, uint32_t byte_count){
    char header[IVF_FRAME_HEADER_SIZE];
    int32_t write_location = 0;

    mem_put_le32(&header[write_location], (int32_t)byte_count);
    write_location = write_location + 4;
    mem_put_le32(&header[write_location], (int32_t)((config->ivf_count) & 0xFFFFFFFF));
    write_location = write_location + 4;
    mem_put_le32(&header[write_location], (int32_t)((config->ivf_count) >> 32));
    write_location = write_location + 4;

    config->byte_count_since_ivf = (byte_count);

    config->ivf_count++;
    fflush(stdout);

    if (config->bitstream_file)
        fwrite(header, 1, IVF_FRAME_HEADER_SIZE, config->bitstream_file);
}

/***************************************
* Process Output STATISTICS Buffer
***************************************/
void process_output_statistics_buffer(
    EbBufferHeaderType      *header_ptr,
    EbConfig                *config){

    uint32_t    max_luma_value = (config->encoder_bit_depth == 8) ? 255 : 1023;
    uint64_t    picture_stream_size, luma_sse, cr_sse, cb_sse, picture_number, picture_qp;
    double      temp_var, luma_psnr, cb_psnr, cr_psnr;

    picture_stream_size  = header_ptr->n_filled_len;
    luma_sse             = header_ptr->luma_sse;
    cr_sse               = header_ptr->cr_sse;
    cb_sse               = header_ptr->cb_sse;
    picture_number       = header_ptr->pts;
    picture_qp           = header_ptr->qp;

    temp_var = (double)max_luma_value*max_luma_value *
             (config->source_width*config->source_height);

    if (luma_sse == 0)
        luma_psnr = 100;
    else
        luma_psnr = 10 * log10((double)temp_var / (double)luma_sse);

    temp_var = (double) max_luma_value * max_luma_value *
            (config->source_width / 2 * config->source_height / 2);

    if (cb_sse == 0)
        cb_psnr = 100;
    else
        cb_psnr = 10 * log10((double)temp_var / (double)cb_sse);

    if (cr_sse == 0)
        cr_psnr = 100;
    else
        cr_psnr = 10 * log10((double)temp_var / (double)cr_sse);

    config->performance_context.sum_luma_psnr += luma_psnr;
    config->performance_context.sum_cr_psnr   += cr_psnr;
    config->performance_context.sum_cb_psnr   += cb_psnr;
    config->performance_context.sum_qp        += picture_qp;

    // Write statistic Data to file
    if (config->stat_file)
        fprintf(config->stat_file, "Picture Number: %4d\t QP: %4d [Y: %.2f dB, U: %.2f dB, V: %.2f dB]\t %6d bits\n",
        (int)picture_number,
        (int)picture_qp,
        luma_psnr,
        cr_psnr,
        cb_psnr,
        (int)picture_stream_size);

    return;
}


AppExitConditionType ProcessOutputStreamBuffer(
    EbConfig             *config,
    EbAppContext         *appCallBack,
    uint8_t                 pic_send_done)
{
    AppPortActiveType      *portState       = &appCallBack->output_stream_port_active;
    EbBufferHeaderType     *headerPtr;
    EbComponentType        *componentHandle = (EbComponentType*)appCallBack->svt_encoder_handle;
    AppExitConditionType    return_value    = APP_ExitConditionNone;
    EbErrorType             stream_status   = EB_ErrorNone;
    // Per channel variables
    FILE                   *streamFile       = config->bitstream_file;

    uint64_t               *total_latency     = &config->performance_context.total_latency;
    uint32_t               *max_latency       = &config->performance_context.max_latency;

    // System performance variables
    static int32_t         frame_count                = 0;

    // Local variables
    uint64_t                finishsTime     = 0;
    uint64_t                finishuTime     = 0;
    uint8_t is_alt_ref = 1;
    while (is_alt_ref) {
        is_alt_ref = 0;
        // non-blocking call until all input frames are sent
        stream_status = eb_svt_get_packet(componentHandle, &headerPtr, pic_send_done);

        if (stream_status == EB_ErrorMax) {
            printf("\n");
            LogErrorOutput(
                config->error_log_file,
                headerPtr->flags);
            return APP_ExitConditionError;
        }
        else if (stream_status != EB_NoErrorEmptyQueue) {
            is_alt_ref = (headerPtr->flags & EB_BUFFERFLAG_IS_ALT_REF);
            EbBool   has_tiles = (EbBool)(appCallBack->eb_enc_parameters.tile_columns || appCallBack->eb_enc_parameters.tile_rows);
            uint8_t  obu_frame_header_size = has_tiles ? OBU_FRAME_HEADER_SIZE + 1 : OBU_FRAME_HEADER_SIZE;
            if (!(headerPtr->flags & EB_BUFFERFLAG_IS_ALT_REF))
                ++(config->performance_context.frame_count);
            *total_latency += (uint64_t)headerPtr->n_tick_count;
            *max_latency = (headerPtr->n_tick_count > *max_latency) ? headerPtr->n_tick_count : *max_latency;

            FinishTime((uint64_t*)&finishsTime, (uint64_t*)&finishuTime);

            // total execution time, inc init time
            ComputeOverallElapsedTime(
                config->performance_context.lib_start_time[0],
                config->performance_context.lib_start_time[1],
                finishsTime,
                finishuTime,
                &config->performance_context.total_execution_time);

            // total encode time
            ComputeOverallElapsedTime(
                config->performance_context.encode_start_time[0],
                config->performance_context.encode_start_time[1],
                finishsTime,
                finishuTime,
                &config->performance_context.total_encode_time);

            // Write Stream Data to file
            if (streamFile) {
                if (config->performance_context.frame_count ==  1 && !(headerPtr->flags & EB_BUFFERFLAG_IS_ALT_REF)){
                    write_ivf_stream_header(config);
                }

                switch (headerPtr->flags & 0x00000006) { // Check for the flags EB_BUFFERFLAG_HAS_TD and EB_BUFFERFLAG_SHOW_EXT
                case (EB_BUFFERFLAG_HAS_TD | EB_BUFFERFLAG_SHOW_EXT):

                    // terminate previous ivf packet, update the combined size of packets sent
                    update_prev_ivf_header(config);

                    // Write a new IVF frame header to file as a TD is in the packet
                    write_ivf_frame_header(config, headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE));
                    fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE), streamFile);

                    // An EB_BUFFERFLAG_SHOW_EXT means that another TD has been added to the packet to show another frame, a new IVF is needed
                    write_ivf_frame_header(config, (obu_frame_header_size + TD_SIZE));
                    fwrite(headerPtr->p_buffer + headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE), 1, (obu_frame_header_size + TD_SIZE), streamFile);

                    break;

                case (EB_BUFFERFLAG_HAS_TD):

                    // terminate previous ivf packet, update the combined size of packets sent
                    update_prev_ivf_header(config);

                    // Write a new IVF frame header to file as a TD is in the packet
                    write_ivf_frame_header(config, headerPtr->n_filled_len);
                    fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, streamFile);

                    break;

                case (EB_BUFFERFLAG_SHOW_EXT):

                    // this case means that there's only one TD in this packet and is relater
                    fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE), streamFile);
                    // this packet will be part of the previous IVF header
                    config->byte_count_since_ivf += (headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE));

                    // terminate previous ivf packet, update the combined size of packets sent
                    update_prev_ivf_header(config);

                    // An EB_BUFFERFLAG_SHOW_EXT means that another TD has been added to the packet to show another frame, a new IVF is needed
                    write_ivf_frame_header(config, (obu_frame_header_size + TD_SIZE));
                    fwrite(headerPtr->p_buffer + headerPtr->n_filled_len - (obu_frame_header_size + TD_SIZE), 1, (obu_frame_header_size + TD_SIZE), streamFile);

                    break;

                default:

                    // This is a packet without a TD, write it straight to file
                    fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, streamFile);

                    // this packet will be part of the previous IVF header
                    config->byte_count_since_ivf += (headerPtr->n_filled_len);
                    break;
                }
            }
            config->performance_context.byte_count += headerPtr->n_filled_len;

            if (config->stat_report && !(headerPtr->flags & EB_BUFFERFLAG_IS_ALT_REF))
                process_output_statistics_buffer(headerPtr, config);

            // Update Output Port Activity State
            *portState = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_PortInactive : *portState;
            return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;

            // Release the output buffer
            eb_svt_release_out_buffer(&headerPtr);

#if DEADLOCK_DEBUG
            ++frame_count;
#else
            //++frame_count;
            if (!(headerPtr->flags & EB_BUFFERFLAG_IS_ALT_REF))
                printf("\b\b\b\b\b\b\b\b\b%9d", ++frame_count);
#endif

            //++frame_count;
            fflush(stdout);

            {
                config->performance_context.average_speed = (config->performance_context.frame_count) / config->performance_context.total_encode_time;
                config->performance_context.average_latency = config->performance_context.total_latency / (double)(config->performance_context.frame_count);
            }

            if (!(frame_count % SPEED_MEASUREMENT_INTERVAL)) {
                {
                    printf("\n");
                    printf("Average System Encoding Speed:        %.2f\n", (double)(frame_count) / config->performance_context.total_encode_time);
                }
            }
        }
    }
    return return_value;
}
AppExitConditionType ProcessOutputReconBuffer(
    EbConfig             *config,
    EbAppContext         *appCallBack)
{
    EbBufferHeaderType    *headerPtr = appCallBack->recon_buffer; // needs to change for buffered input
    EbComponentType       *componentHandle = (EbComponentType*)appCallBack->svt_encoder_handle;
    AppExitConditionType    return_value = APP_ExitConditionNone;
    EbErrorType            recon_status = EB_ErrorNone;
    int32_t fseekReturnVal;
    // non-blocking call until all input frames are sent
    recon_status = eb_svt_get_recon(componentHandle, headerPtr);

    if (recon_status == EB_ErrorMax) {
        printf("\n");
        LogErrorOutput(
            config->error_log_file,
            headerPtr->flags);
        return APP_ExitConditionError;
    }
    else if (recon_status != EB_NoErrorEmptyQueue) {
        //Sets the File position to the beginning of the file.
        rewind(config->recon_file);
        uint64_t frameNum = headerPtr->pts;
        while (frameNum>0) {
            fseekReturnVal = fseeko64(config->recon_file, headerPtr->n_filled_len, SEEK_CUR);

            if (fseekReturnVal != 0) {
                printf("Error in fseeko64  returnVal %i\n", fseekReturnVal);
                return APP_ExitConditionError;
            }
            frameNum = frameNum - 1;
        }

        fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, config->recon_file);

        // Update Output Port Activity State
        return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;
    }
    return return_value;
}
