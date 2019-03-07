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

#include "EbAppContext.h"
#include "EbAppConfig.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbAppInputy4m.h"

#include "EbSvtAv1Time.h"


#define IVF_FRAME_HEADER_IN_LIB                     0

/***************************************
 * Macros
 ***************************************/
#define CLIP3(MinVal, MaxVal, a)        (((a)<(MinVal)) ? (MinVal) : (((a)>(MaxVal)) ? (MaxVal) :(a)))
#define FUTURE_WINDOW_WIDTH                 4
#define SIZE_OF_ONE_FRAME_IN_BYTES(width, height,is16bit) ( ( ((width)*(height)*3)>>1 )<<is16bit)
#define YUV4MPEG2_IND_SIZE 9
extern volatile int32_t keepRunning;

/***************************************
* Process Error Log
***************************************/
void LogErrorOutput(
    FILE                     *errorLogFile,
    uint32_t                  errorCode)
{

    switch (errorCode) {

        // EB_ENC_AMVP_ERRORS:
    case EB_ENC_AMVP_ERROR1:
        fprintf(errorLogFile, "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR2:
        fprintf(errorLogFile, "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

    case EB_ENC_AMVP_ERROR3:
        fprintf(errorLogFile, "Error: The input PU to GetNonScalingSpatialAMVP() can not be I_MODE!\n");
        break;

    case EB_ENC_AMVP_ERROR4:
        fprintf(errorLogFile, "Error: The availability parameter in GetSpatialMVPPosAx() function can not be > 3 !\n");
        break;

    case EB_ENC_AMVP_ERROR5:
        fprintf(errorLogFile, "Error: The availability parameter in GetSpatialMVPPosBx() function can not be > 7 !\n");
        break;

    case EB_ENC_AMVP_ERROR6:
        fprintf(errorLogFile, "Error: GetTemporalMVP: tmvpMapLcuIndex must be either 0 or 1!\n");
        break;

    case EB_ENC_AMVP_ERROR7:
        fprintf(errorLogFile, "Error: the input PU to GetNonScalingSpatialAMVP() must be available!");
        break;

    case EB_ENC_AMVP_ERROR8:
        fprintf(errorLogFile, "Error: GetTemporalMVP: tmvpMapLcuIndex must be either 0 or 1");
        break;

    case EB_ENC_AMVP_NULL_REF_ERROR:
        fprintf(errorLogFile, "Error: The referenceObject can not be NULL!\n");
        break;

    case EB_ENC_AMVP_SPATIAL_NA_ERROR:
        fprintf(errorLogFile, "Error: The input PU to GetNonScalingSpatialAMVP() must be available!\n");
        break;

        // EB_ENC_CL_ERRORS:
    case EB_ENC_CL_ERROR1:
        fprintf(errorLogFile, "Error: Unknown Inter Prediction Direction!\n");
        break;

    case EB_ENC_CL_ERROR2:
        fprintf(errorLogFile, "Error: Unknown coding mode!\n");
        break;

    case EB_ENC_CL_ERROR3:
        fprintf(errorLogFile, "Error: Mode Decision Candidate Buffer Overflow!\n");
        break;

    case EB_ENC_CL_ERROR4:
        fprintf(errorLogFile, "Error: Too many Mode Decision Fast Candidates!\n");
        break;

    case EB_ENC_CL_ERROR5:
        fprintf(errorLogFile, "Error: Too many buffers chosen for this level by PreModeDecision!\n");
        break;

    case EB_ENC_CL_ERROR6:
        fprintf(errorLogFile, "Error: Ping-Pong structure needs at least two buffers to work properly!\n");
        break;

    case EB_ENC_CL_ERROR7:
        fprintf(errorLogFile, "Error: Invalid Intra Partition\n");
        break;

    case EB_ENC_CL_ERROR8:
        fprintf(errorLogFile, "Error: Invalid TU Configuration\n");
        break;
    case EB_ENC_CL_ERROR9:
        fprintf(errorLogFile, "Error: Invalid Prediction Mode\n");
        break;


        // EB_ENC_DLF_ERRORS:
    case EB_ENC_DLF_ERROR1:
        fprintf(errorLogFile, "Error: While calculating bS for DLF!\n");
        break;

    case EB_ENC_DLF_ERROR2:
        fprintf(errorLogFile, "Error: Unknown Inter Prediction Direction Combination!\n");
        break;

    case EB_ENC_DLF_ERROR3:
        fprintf(errorLogFile, "Error: If any PU is in I_MODE, the bS will be 2!\n");
        break;

    case EB_ENC_DLF_ERROR4:
        fprintf(errorLogFile, "Error: The x/y location of the CU must be the multiple of minmum CU size!");
        break;

    case EB_ENC_DLF_ERROR5:
        fprintf(errorLogFile, "Error: Unknown Slice Type!");
        break;

    case EB_ENC_DLF_ERROR6:
        fprintf(errorLogFile, "Error: While calculating the bS for the PU bounday, the 4x4 block must be guaranteed at the PU boundary!");
        break;

    case EB_ENC_DLF_ERROR7:
        fprintf(errorLogFile, "Error: LCU size must be power of 2!");
        break;

    case EB_ENC_DLF_ERROR8:
        fprintf(errorLogFile, "Error: Deblocking filter can not support the picture whose width or height is not the multiple of 8!");
        break;

    case EB_ENC_DLF_ERROR9:
        fprintf(errorLogFile, "Error: Neighbor PU must be available!");
        break;

    case EB_ENC_DLF_ERROR10:
        fprintf(errorLogFile, "Error: Deblocking filter can not support the picture whose width or height is not the multiple of 8!");
        break;



        // EB_ENC_EC_ERRORS:
    case EB_ENC_EC_ERROR1:
        fprintf(errorLogFile, "Error: EncodeCodedBlockFlags: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR10:
        fprintf(errorLogFile, "Error: EncodeTuSplitCoeff: context value too large!\n");
        break;

    case EB_ENC_EC_ERROR11:
        fprintf(errorLogFile, "Error: CodeSPS: Long term reference pictures are not currently handled!\n");
        break;

    case EB_ENC_EC_ERROR12:
        fprintf(errorLogFile, "Error: CodeProfileTierLevel: The maximum sublayers must be equal to 1!\n");
        break;

    case EB_ENC_EC_ERROR13:
        fprintf(errorLogFile, "Error: EncodeRootCodedBlockFlag: rootCbf too large!\n");
        break;

    case EB_ENC_EC_ERROR14:
        fprintf(errorLogFile, "Error: cpbCountMinus1 in HRD parameter exceeds the upper limit 4!\n");
        break;

    case EB_ENC_EC_ERROR15:
        fprintf(errorLogFile, "Error: numDecodingUnitsMinus1 in picture timeing SEI exceeds the upper limit 64!\n");
        break;

    case EB_ENC_EC_ERROR16:
        fprintf(errorLogFile, "Error: The size of the unregistered user data SEI payload is not allowed!\n");
        break;

    case EB_ENC_EC_ERROR2:
        fprintf(errorLogFile, "Error: CopyRbspBitstreamToPayload: output buffer too small!\n");
        break;

    case EB_ENC_EC_ERROR3:
        fprintf(errorLogFile, "Error: EncodeLcu: Unknown mode type!\n");
        break;

    case EB_ENC_EC_ERROR4:
        fprintf(errorLogFile, "Error: 8x4 & 4x8 PU should not have Bi-pred mode!\n");
        break;

    case EB_ENC_EC_ERROR5:
        fprintf(errorLogFile, "Error: EncodeMergeIndex: value too large!\n");
        break;

    case EB_ENC_EC_ERROR6:
        fprintf(errorLogFile, "Error: EncodeSkipFlag: context too large!\n");
        break;

    case EB_ENC_EC_ERROR7:
        fprintf(errorLogFile, "Error: EncodeBypassBins: binsLength must be less than 32!\n");
        break;

    case EB_ENC_EC_ERROR8:
        fprintf(errorLogFile, "Error: EncodeQuantizedCoefficients: Invalid block size!\n");
        break;

    case EB_ENC_EC_ERROR9:
        fprintf(errorLogFile, "Error: EncodeSplitFlag: context too large!\n");
        break;

    case EB_ENC_EC_ERROR26:
        fprintf(errorLogFile, "Error: Level not recognized!\n");
        break;

    case EB_ENC_EC_ERROR27:
        fprintf(errorLogFile, "Error: EncodeOneBin:  BinaryValue must be less than 2\n");
        break;

    case EB_ENC_EC_ERROR28:
        fprintf(errorLogFile, "Error: No more than 6 SAO types\n");
        break;

    case EB_ENC_EC_ERROR29:
        fprintf(errorLogFile, "Error: No more than 6 SAO types\n");
        break;

        // EB_ENC_FL_ERRORS:
    case EB_ENC_FL_ERROR1:
        fprintf(errorLogFile, "Error: Uncovered area inside Cu!\n");
        break;

    case EB_ENC_FL_ERROR2:
        fprintf(errorLogFile, "Error: Depth 2 is not allowed for 8x8 CU!\n");
        break;

    case EB_ENC_FL_ERROR3:
        fprintf(errorLogFile, "Error: Depth 0 is not allowed for 64x64 CU!\n");
        break;

    case EB_ENC_FL_ERROR4:
        fprintf(errorLogFile, "Error: Max CU Depth Exceeded!\n");
        break;

        // EB_ENC_HANDLE_ERRORS:
    case EB_ENC_HANDLE_ERROR1:
        fprintf(errorLogFile, "Error: Only one Resource Coordination Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR10:
        fprintf(errorLogFile, "Error: Need at least one Entropy Coding Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR11:
        fprintf(errorLogFile, "Error: Only one Packetization Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR12:
        fprintf(errorLogFile, "Error: RC Results Fifo Size should be greater than RC Tasks Fifo Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR13:
        fprintf(errorLogFile, "Error: RC Tasks Fifo Size should be greater than EC results Fifo Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR14:
        fprintf(errorLogFile, "Error: RC Tasks Fifo Size should be greater than Picture Manager results Fifo Size in order to avoid deadlock!\n");
        break;

    case EB_ENC_HANDLE_ERROR18:
        fprintf(errorLogFile, "Error: Intra period setting breaks mini-gop!\n");
        break;

    case EB_ENC_HANDLE_ERROR2:
        fprintf(errorLogFile, "Error: Only one Picture Enhancement Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR3:
        fprintf(errorLogFile, "Error: Only one Picture Manager Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR4:
        fprintf(errorLogFile, "Error: Need at least one ME Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR5:
        fprintf(errorLogFile, "Error: Only one Rate-Control Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR6:
        fprintf(errorLogFile, "Error: Need at least one Mode Decision Configuration Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR7:
        fprintf(errorLogFile, "Error: Need at least one Coding Loop Process!\n");
        break;

    case EB_ENC_HANDLE_ERROR8:
        fprintf(errorLogFile, "Error: Only one Second Pass Deblocking Process allowed!\n");
        break;

    case EB_ENC_HANDLE_ERROR9:
        fprintf(errorLogFile, "Error: Only one ALF Process allowed!\n");
        break;

        // EB_ENC_INTER_ERRORS:
    case EB_ENC_INTER_INVLD_MCP_ERROR:
        fprintf(errorLogFile, "Error: Motion compensation prediction is out of the picture boundary!\n");
        break;

    case EB_ENC_INTER_PRED_ERROR0:
        fprintf(errorLogFile, "Error: Unkown Inter Prediction Direction!\n");
        break;

    case EB_ENC_INTER_PRED_ERROR1:
        fprintf(errorLogFile, "Error: Inter prediction can not support more than 2 MVPs!\n");
        break;

        // EB_ENC_INTRA_ERRORS:
    case EB_ENC_INTRA_PRED_ERROR1:
        fprintf(errorLogFile, "Error: IntraPrediction does not support 2Nx2N partition size!\n");
        break;

    case EB_ENC_INTRA_PRED_ERROR2:
        fprintf(errorLogFile, "Error: IntraPrediction: intra prediction only supports square PU!\n");
        break;

    case EB_ENC_INTRA_PRED_ERROR3:
        fprintf(errorLogFile, "Error: IntraPredictionChroma: Only Planar!\n");
        break;

    case EB_ENC_INVLD_PART_SIZE_ERROR:
        fprintf(errorLogFile, "Error: IntraPrediction: only PU sizes of 8 or largers are currently supported!\n");
        break;


        // EB_ENC_MD_ERRORS:
    case EB_ENC_MD_ERROR1:
        fprintf(errorLogFile, "Error: Unknown AMVP Mode Decision Candidate Type!\n");
        break;

    case EB_ENC_MD_ERROR2:
        fprintf(errorLogFile, "Error: PreModeDecision: need at least one buffer!\n");
        break;

    case EB_ENC_MD_ERROR3:
        fprintf(errorLogFile, "Error: Unknow Inter Prediction Direction!\n");
        break;

    case EB_ENC_MD_ERROR4:
        fprintf(errorLogFile, "Error: Unknown ME SAD Level!\n");
        break;

    case EB_ENC_MD_ERROR5:
        fprintf(errorLogFile, "Error: Invalid encoder mode. The encoder mode should be 0, 1 or 2!\n");
        break;

    case EB_ENC_MD_ERROR6:
        fprintf(errorLogFile, "Error: Invalid TU size!\n");
        break;

    case EB_ENC_MD_ERROR7:
        fprintf(errorLogFile, "Error: Unknown depth!\n");
        break;

    case EB_ENC_MD_ERROR8:
        fprintf(errorLogFile, "Error: Depth not supported!\n");
        break;

    case EB_ENC_MD_ERROR9:
        fprintf(errorLogFile, "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

    case EB_ENC_MD_ERROR10:
        fprintf(errorLogFile, "Error: Ping-Pong structure needs at least two buffers to work properly\n");
        break;

        // EB_ENC_ME_ERRORS:
    case EB_ENC_ME_ERROR1:
        fprintf(errorLogFile, "Error: Motion Estimation: non valid value of the subPelDirection !\n");
        break;

    case EB_ENC_ME_ERROR2:
        fprintf(errorLogFile, "Error: FillMvMergeCandidate() method only supports P or B slices!\n");
        break;

        // EB_ENC_ERRORS:
    case EB_ENC_ROB_OF_ERROR:
        fprintf(errorLogFile, "Error: Recon Output Buffer Overflow!\n");
        break;

        // EB_ENC_PACKETIZATION_ERRORS:
    case EB_ENC_PACKETIZATION_ERROR1:
        fprintf(errorLogFile, "Error: PacketizationProcess: Picture Number does not match entry. PacketizationReorderQueue overflow!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR2:
        fprintf(errorLogFile, "Error: Entropy Coding Result can not be outputed by processes other than entropy coder and ALF!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR3:
        fprintf(errorLogFile, "Error: The encoder can not support the SliceMode other than 0 and 1!\n");
        break;

    case EB_ENC_PACKETIZATION_ERROR4:
        fprintf(errorLogFile, "Error: Statistics Output Buffer Overflow!\n");
        break;
    case EB_ENC_PACKETIZATION_ERROR5:
        fprintf(errorLogFile, "Error: Stream Fifo is starving..deadlock, increase EB_outputStreamBufferFifoInitCount APP_ENCODERSTREAMBUFFERCOUNT \n");
        break;

        // EB_ENC_PM_ERRORS:
    case EB_ENC_PM_ERROR0:
        fprintf(errorLogFile, "Error: PictureManagerProcess: Unknown Slice Type!\n");
        break;

    case EB_ENC_PM_ERROR1:
        fprintf(errorLogFile, "Error: EbPictureManager: dependentCount underflow!\n");
        break;

    case EB_ENC_PM_ERROR10:
        fprintf(errorLogFile, "Error: picture_manager_kernel: referenceEntryPtr should never be null!\n");
        break;

    case EB_ENC_PM_ERROR2:
        fprintf(errorLogFile, "Error: PictureManagerProcess: The Reference Structure period must be less than the MAX_ELAPSED_IDR_COUNT or false-IDR boundary logic will be activated!\n");
        break;

    case EB_ENC_PM_ERROR3:
        fprintf(errorLogFile, "Error: PictureManagerProcess: The dependentCount underflow detected!\n");
        break;

    case EB_ENC_PM_ERROR4:
        fprintf(errorLogFile, "Error: PictureManagerProcess: Empty input queue!\n");
        break;

    case EB_ENC_PM_ERROR5:
        fprintf(errorLogFile, "Error: PictureManagerProcess: Empty reference queue!\n");
        break;

    case EB_ENC_PM_ERROR6:
        fprintf(errorLogFile, "Error: PictureManagerProcess: The capped elaspedNonIdrCount must be larger than the maximum supported delta ref poc!\n");
        break;

    case EB_ENC_PM_ERROR7:
        fprintf(errorLogFile, "Error: PictureManagerProcess: Reference Picture Queue Full!\n");
        break;

    case EB_ENC_PM_ERROR8:
        fprintf(errorLogFile, "Error: PictureManagerProcess: No reference match found - this will lead to a memory leak!\n");
        break;

    case EB_ENC_PM_ERROR9:
        fprintf(errorLogFile, "Error: PictureManagerProcess: Unknown picture type!\n");
        break;

    case EB_ENC_PM_ERROR12:
        fprintf(errorLogFile, "Error: PictureManagerProcess: prediction structure configuration API has too many reference pictures\n");
        break;

    case EB_ENC_PM_ERROR13:
        fprintf(errorLogFile, "Error: PictureManagerProcess: The maximum allowed frame rate is 60 fps\n");
        break;

    case EB_ENC_PM_ERROR14:
        fprintf(errorLogFile, "Error: PictureManagerProcess: The minimum allowed frame rate is 1 fps\n");
        break;

        // EB_ENC_PRED_STRC_ERRORS:
    case EB_ENC_PRED_STRC_ERROR1:
        fprintf(errorLogFile, "Error: PredictionStructureCtor: DecodeOrder LUT too small!\n");
        break;

    case EB_ENC_PRED_STRC_ERROR2:
        fprintf(errorLogFile, "Error: PredictionStructureCtor: prediction structure improperly configured!\n");
        break;

        // EB_ENC_PU_ERRORS:
    case EB_ENC_PU_ERROR1:
        fprintf(errorLogFile, "Error: Unknown partition size!\n");
        break;

    case EB_ENC_PU_ERROR2:
        fprintf(errorLogFile, "Error: The target area is not inside the CU!\n");
        break;

        // EB_ENC_RC_ERRORS:
    case EB_ENC_RC_ERROR1:
        fprintf(errorLogFile, "Error: RateControlProcess: Unknown input tasktype!\n");
        break;

    case EB_ENC_RC_ERROR2:
        fprintf(errorLogFile, "Error: RateControlProcess: No RC interval found!\n");
        break;

    case EB_ENC_RC_ERROR3:
        fprintf(errorLogFile, "Error: RateControlProcess: RC input Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR4:
        fprintf(errorLogFile, "Error: RateControlProcess: RC feedback Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR5:
        fprintf(errorLogFile, "Error: RateControlProcess: RC feedback Picture Queue Full!\n");
        break;

    case EB_ENC_RC_ERROR6:
        fprintf(errorLogFile, "Error: RateControlProcess: No feedback frame match found - this will lead to a memory leak!\n");
        break;

    case EB_ENC_RC_ERROR7:
        fprintf(errorLogFile, "Error: remainingBytes has to be multiple of 2 for 16 bit input\n");
        break;

    case EB_ENC_RC_ERROR8:
        fprintf(errorLogFile, "Error: hlRateControlHistorgramQueue Overflow\n");
        break;

        // EB_ENC_RD_COST_ERRORS:
    case EB_ENC_RD_COST_ERROR1:
        fprintf(errorLogFile, "Error: Skip mode only exists in 2Nx2N partition type!\n");
        break;

    case EB_ENC_RD_COST_ERROR2:
        fprintf(errorLogFile, "Error: IntraChromaCost: Unknown slice type!\n");
        break;

    case EB_ENC_RD_COST_ERROR3:
        fprintf(errorLogFile, "Error: intra2_nx2_n_fast_cost_islice can only support 2Nx2N partition type!\n");
        break;

        // EB_ENC_SAO_ERRORS:
    case EB_ENC_SAO_ERROR1:
        fprintf(errorLogFile, "Error: No more than 6 SAO types!\n");
        break;

    case EB_ENC_SAO_ERROR2:
        fprintf(errorLogFile, "Error: No more than 5 EO SAO categories!\n");
        break;
        // EB_ENC_SCS_ERRORS:
    case EB_ENC_SCS_ERROR1:
        fprintf(errorLogFile, "Error: SequenceControlSetCopy: Not all SequenceControlSet_t members are being copied!\n");
        break;

        // EB_ENC_BITSTREAM_ERRORS:
    case EB_ENC_BITSTREAM_ERROR1:
        fprintf(errorLogFile, "Error: OutputBitstreamRBSPToPayload: Bitstream payload buffer empty!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR2:
        fprintf(errorLogFile, "Error: OutputBitstreamWrite: Empty bitstream!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR3:
        fprintf(errorLogFile, "Error: OutputBitstreamRBSPToPayload: Buffer index more than buffer size!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR4:
        fprintf(errorLogFile, "Error: OutputBitstreamRBSPToPayload: Start Location in not inside the buffer!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR5:
        fprintf(errorLogFile, "Error: OutputBitstreamWrite: Trying to write more than one word!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR6:
        fprintf(errorLogFile, "Error: OutputBitstreamRBSPToPayload: Expecting Start code!\n");
        break;

    case EB_ENC_BITSTREAM_ERROR7:
        fprintf(errorLogFile, "Error: OutputBitstreamRBSPToPayload: Bitstream not flushed (i.e. byte-aligned)!\n");
        break;

    case EB_ENC_RESS_COOR_ERRORS1:
        fprintf(errorLogFile, "Error: ResourceCoordinationProcess: The received input data should be equal to the buffer size - only complete frame transfer is supported\n");
        break;

    case EB_ENC_RES_COORD_InvalidQP:
        fprintf(errorLogFile, "Error: ResourceCoordinationProcess: The QP value in the QP file is invalid\n");
        break;

    case EB_ENC_RES_COORD_InvalidSliceType:
        fprintf(errorLogFile, "Error: ResourceCoordinationProcess: Slice Type Invalid\n");
        break;

        // picture decision Errors
    case EB_ENC_PD_ERROR8:
        fprintf(errorLogFile, "Error: PictureDecisionProcess: Picture Decision Reorder Queue overflow\n");
        break;

    default:
        fprintf(errorLogFile, "Error: Others!\n");
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

    EbConfig_t               *config,
    EbBufferHeaderType      *headerPtr,
    FILE                     *inputFile,
    uint8_t                    *lumaInputPtr,
    uint8_t                    *cbInputPtr,
    uint8_t                    *crInputPtr,
    uint8_t                   is16bit) {


    int64_t  inputPaddedWidth  = config->inputPaddedWidth;
    int64_t  inputPaddedHeight = config->inputPaddedHeight;
    uint64_t  sourceLumaRowSize = (uint64_t)(inputPaddedWidth << is16bit);
    uint64_t  sourceChromaRowSize = sourceLumaRowSize >> 1;
    uint8_t  *ebInputPtr;
    uint32_t  inputRowIndex;

    // Y
    ebInputPtr = lumaInputPtr;
    // Skip 1 luma row if bottom field (point to the bottom field)
    if (config->processedFrameCount % 2 != 0)
        fseeko64(inputFile, (long)sourceLumaRowSize, SEEK_CUR);

    for (inputRowIndex = 0; inputRowIndex < inputPaddedHeight; inputRowIndex++) {

        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, sourceLumaRowSize, inputFile);
        // Skip 1 luma row (only fields)
        fseeko64(inputFile, (long)sourceLumaRowSize, SEEK_CUR);
        ebInputPtr += sourceLumaRowSize;
    }

    // U
    ebInputPtr = cbInputPtr;
    // Step back 1 luma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    if (config->processedFrameCount % 2 != 0) {
        fseeko64(inputFile, -(long)sourceLumaRowSize, SEEK_CUR);
        fseeko64(inputFile, (long)sourceChromaRowSize, SEEK_CUR);
    }

    for (inputRowIndex = 0; inputRowIndex < inputPaddedHeight >> 1; inputRowIndex++) {

        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, sourceChromaRowSize, inputFile);
        // Skip 1 chroma row (only fields)
        fseeko64(inputFile, (long)sourceChromaRowSize, SEEK_CUR);
        ebInputPtr += sourceChromaRowSize;
    }

    // V
    ebInputPtr = crInputPtr;
    // Step back 1 chroma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    // => no action


    for (inputRowIndex = 0; inputRowIndex < inputPaddedHeight >> 1; inputRowIndex++) {

        headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, sourceChromaRowSize, inputFile);
        // Skip 1 chroma row (only fields)
        fseeko64(inputFile, (long)sourceChromaRowSize, SEEK_CUR);
        ebInputPtr += sourceChromaRowSize;
    }

    // Step back 1 chroma row if bottom field (undo the previous jump)
    if (config->processedFrameCount % 2 != 0) {
        fseeko64(inputFile, -(long)sourceChromaRowSize, SEEK_CUR);
    }
}


//************************************/
// GetNextQpFromQpFile
// Reads and extracts one qp from the qpfile
// Input  : QP file
// Output : QP value
/************************************/
static int32_t qpReadFromFile = 0;

int32_t GetNextQpFromQpFile(
    EbConfig_t  *config
)
{
    uint8_t *line;
    int32_t qp = 0;
    uint32_t readsize = 0, eof = 0;
    EB_APP_MALLOC(uint8_t*, line, 8, EB_N_PTR, EB_ErrorInsufficientResources);
    memset(line,0,8);
    readsize = (uint32_t)fread(line, 1, 2, config->qpFile);

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
        fseek(config->qpFile, -1, SEEK_CUR);
        qp = 0;
    }
    else if (readsize == 2 && (line[1] == '\n')) {
        // new line
        qp = strtol((const char*)line, NULL, 0);
    }
    else if (readsize == 2 && (line[0] == '#' || line[0] == '/' || line[0] == '-' || line[0] == ' ')) {
        // Backup one step to not miss the new line char
        fseek(config->qpFile, -1, SEEK_CUR);
        do {
            readsize = (uint32_t)fread(line, 1, 1, config->qpFile);
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
            readsize = (uint32_t)fread(line, 1, 1, config->qpFile);
            if (readsize != 1)
                break;
        } while (line[0] != '\n');
    }

    if (qp > 0)
        qpReadFromFile |= 1;

    return qp;
}

void ReadInputFrames(
    EbConfig_t                  *config,
    uint8_t                      is16bit,
    EbBufferHeaderType         *headerPtr){

    uint64_t  readSize;
    uint32_t  inputPaddedWidth = config->inputPaddedWidth;
    uint32_t  inputPaddedHeight = config->inputPaddedHeight;
    FILE   *inputFile = config->inputFile;
    uint8_t  *ebInputPtr;
    EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;

    uint64_t frameSize = (uint64_t)((inputPaddedWidth*inputPaddedHeight * 3) / 2 + (inputPaddedWidth / 4 * inputPaddedHeight * 3) / 2);
    inputPtr->yStride  = inputPaddedWidth;
    inputPtr->crStride = inputPaddedWidth >> 1;
    inputPtr->cbStride = inputPaddedWidth >> 1;

    if (config->bufferedInput == -1) {

        if (is16bit == 0 || (is16bit == 1 && config->compressedTenBitFormat == 0)) {

            readSize = (uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(inputPaddedWidth, inputPaddedHeight, is16bit);

            headerPtr->n_filled_len = 0;

            // Interlaced Video
            if (config->separateFields) {

                ProcessInputFieldStandardMode(
                    config,
                    headerPtr,
                    inputFile,
                    inputPtr->luma,
                    inputPtr->cb,
                    inputPtr->cr,
                    is16bit);

                if (readSize != headerPtr->n_filled_len) {

                    fseek(inputFile, 0, SEEK_SET);
                    headerPtr->n_filled_len = 0;

                    ProcessInputFieldStandardMode(
                        config,
                        headerPtr,
                        inputFile,
                        inputPtr->luma,
                        inputPtr->cb,
                        inputPtr->cr,
                        is16bit);
                }

                // Reset the pointer position after a top field
                if (config->processedFrameCount % 2 == 0) {
                    fseek(inputFile, -(long)(readSize << 1), SEEK_CUR);
                }
            }
            else {

                /* if input is a y4m file, read next line which contains "FRAME" */
                if(config->y4mInput==EB_TRUE) {
                    readY4mFrameDelimiter(config);
                }

                uint64_t lumaReadSize = (uint64_t)inputPaddedWidth*inputPaddedHeight << is16bit;
                ebInputPtr = inputPtr->luma;
                if(config->y4mInput==EB_FALSE && config->processedFrameCount == 0 && config->inputFile == stdin) {
                    /* if not a y4m file and input is read from stdin, 9 bytes were already read when checking
                        or the YUV4MPEG2 string in the stream, so copy those bytes over */
                    memcpy(ebInputPtr,config->y4mBuf,YUV4MPEG2_IND_SIZE);
                    headerPtr->n_filled_len += YUV4MPEG2_IND_SIZE;
                    ebInputPtr += YUV4MPEG2_IND_SIZE;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize-YUV4MPEG2_IND_SIZE, inputFile);
                }else {
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                }
                ebInputPtr = inputPtr->cb;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->cr;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING) << is16bit);
                inputPtr->cb   = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << is16bit);
                inputPtr->cr   = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << is16bit);


                if (readSize != headerPtr->n_filled_len) {

                    fseek(inputFile, 0, SEEK_SET);
                    ebInputPtr = inputPtr->luma;
                    headerPtr->n_filled_len = (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                    ebInputPtr = inputPtr->cb;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                    ebInputPtr = inputPtr->cr;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                    inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
                    inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
                    inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

                }
            }
        }
        // 10-bit Compressed Unpacked Mode
        else if (is16bit == 1 && config->compressedTenBitFormat == 1) {

            // Fill the buffer with a complete frame
            headerPtr->n_filled_len = 0;


            uint64_t lumaReadSize = (uint64_t)inputPaddedWidth*inputPaddedHeight;
            uint64_t nbitlumaReadSize = (uint64_t)(inputPaddedWidth / 4)*inputPaddedHeight;

            ebInputPtr = inputPtr->luma;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
            ebInputPtr = inputPtr->cb;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
            ebInputPtr = inputPtr->cr;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

            inputPtr->luma = inputPtr->luma + config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING;
            inputPtr->cb = inputPtr->cb + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);
            inputPtr->cr = inputPtr->cr + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);


            ebInputPtr = inputPtr->lumaExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize, inputFile);
            ebInputPtr = inputPtr->cbExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);
            ebInputPtr = inputPtr->crExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);

            inputPtr->lumaExt = inputPtr->lumaExt + ((config->inputPaddedWidth >> 2)*TOP_INPUT_PADDING + (LEFT_INPUT_PADDING >> 2));
            inputPtr->cbExt = inputPtr->cbExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));
            inputPtr->crExt = inputPtr->crExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));

            readSize = ((lumaReadSize * 3) >> 1) + ((nbitlumaReadSize * 3) >> 1);

            if (readSize != headerPtr->n_filled_len) {

                fseek(inputFile, 0, SEEK_SET);
                ebInputPtr = inputPtr->luma;
                headerPtr->n_filled_len = (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                ebInputPtr = inputPtr->cb;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->cr;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                inputPtr->luma = inputPtr->luma + config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING;
                inputPtr->cb = inputPtr->cb + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);
                inputPtr->cr = inputPtr->cr + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);

                ebInputPtr = inputPtr->lumaExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize, inputFile);
                ebInputPtr = inputPtr->cbExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->crExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);

                inputPtr->lumaExt = inputPtr->lumaExt + ((config->inputPaddedWidth >> 2)*TOP_INPUT_PADDING + (LEFT_INPUT_PADDING >> 2));
                inputPtr->cbExt = inputPtr->cbExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));
                inputPtr->crExt = inputPtr->crExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));

            }

        }

        // 10-bit Unpacked Mode
        else {

            readSize = (uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(inputPaddedWidth, inputPaddedHeight, 1);

            headerPtr->n_filled_len = 0;

            // Interlaced Video
            if (config->separateFields) {

                ProcessInputFieldStandardMode(
                    config,
                    headerPtr,
                    inputFile,
                    inputPtr->luma,
                    inputPtr->cb,
                    inputPtr->cr,
                    0);

                ProcessInputFieldStandardMode(
                    config,
                    headerPtr,
                    inputFile,
                    inputPtr->lumaExt,
                    inputPtr->cbExt,
                    inputPtr->crExt,
                    0);

                if (readSize != headerPtr->n_filled_len) {

                    fseek(inputFile, 0, SEEK_SET);
                    headerPtr->n_filled_len = 0;

                    ProcessInputFieldStandardMode(
                        config,
                        headerPtr,
                        inputFile,
                        inputPtr->luma,
                        inputPtr->cb,
                        inputPtr->cr,
                        0);

                    ProcessInputFieldStandardMode(
                        config,
                        headerPtr,
                        inputFile,
                        inputPtr->lumaExt,
                        inputPtr->cbExt,
                        inputPtr->crExt,
                        0);
                }

                // Reset the pointer position after a top field
                if (config->processedFrameCount % 2 == 0) {
                    fseek(inputFile, -(long)(readSize << 1), SEEK_CUR);
                }

            }
            else {


                uint64_t lumaReadSize = (uint64_t)inputPaddedWidth*inputPaddedHeight;

                ebInputPtr = inputPtr->luma;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                ebInputPtr = inputPtr->cb;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->cr;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
                inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
                inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

                ebInputPtr = inputPtr->lumaExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                ebInputPtr = inputPtr->cbExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->crExt;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                inputPtr->lumaExt = inputPtr->lumaExt + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
                inputPtr->cbExt = inputPtr->cbExt + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
                inputPtr->crExt = inputPtr->crExt + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

                if (readSize != headerPtr->n_filled_len) {

                    fseek(inputFile, 0, SEEK_SET);
                    ebInputPtr = inputPtr->luma;
                    headerPtr->n_filled_len = (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                    ebInputPtr = inputPtr->cb;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                    ebInputPtr = inputPtr->cr;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                    inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
                    inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
                    inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

                    ebInputPtr = inputPtr->lumaExt;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                    ebInputPtr = inputPtr->cbExt;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                    ebInputPtr = inputPtr->crExt;
                    headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

                    inputPtr->lumaExt = inputPtr->lumaExt + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
                    inputPtr->cbExt = inputPtr->cbExt + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
                    inputPtr->crExt = inputPtr->crExt + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

                }

            }




        }

    }
    else {
        if (config->encoderBitDepth == 10 && config->compressedTenBitFormat == 1)
        {
            // Determine size of each plane

            const size_t luma8bitSize = config->inputPaddedWidth * config->inputPaddedHeight;
            const size_t chroma8bitSize = luma8bitSize >> 2;

            const size_t luma2bitSize = luma8bitSize / 4; //4-2bit pixels into 1 byte
            const size_t chroma2bitSize = luma2bitSize >> 2;

            EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;
            inputPtr->yStride = config->inputPaddedWidth;
            inputPtr->crStride = config->inputPaddedWidth >> 1;
            inputPtr->cbStride = config->inputPaddedWidth >> 1;

            inputPtr->luma = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput];
            inputPtr->cb = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize;
            inputPtr->cr = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + chroma8bitSize;

            inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING));
            inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));
            inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)));

            if (is16bit) {
                inputPtr->lumaExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize;
                inputPtr->cbExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize + luma2bitSize;
                inputPtr->crExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize + luma2bitSize + chroma2bitSize;

                inputPtr->lumaExt = inputPtr->lumaExt + config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING;
                inputPtr->cbExt = inputPtr->cbExt + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);
                inputPtr->crExt = inputPtr->crExt + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);

            }

            headerPtr->n_filled_len = (uint32_t)frameSize;
        }
        else
        {
            const int32_t tenBitPackedMode = (config->encoderBitDepth > 8) && (config->compressedTenBitFormat == 0) ? 1 : 0;

            // Determine size of each plane
            const size_t luma8bitSize =
                (config->inputPaddedWidth) *
                (config->inputPaddedHeight) *
                (1 << tenBitPackedMode);

            const size_t chroma8bitSize = luma8bitSize >> 2;

            const size_t luma10bitSize = (config->encoderBitDepth > 8 && tenBitPackedMode == 0) ? luma8bitSize : 0;
            const size_t chroma10bitSize = (config->encoderBitDepth > 8 && tenBitPackedMode == 0) ? chroma8bitSize : 0;

            EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)headerPtr->p_buffer;

            inputPtr->yStride = config->inputPaddedWidth;
            inputPtr->crStride = config->inputPaddedWidth >> 1;
            inputPtr->cbStride = config->inputPaddedWidth >> 1;

            inputPtr->luma = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput];
            inputPtr->cb = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize;
            inputPtr->cr = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + chroma8bitSize;
            inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING) << tenBitPackedMode);
            inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << tenBitPackedMode);
            inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << tenBitPackedMode);


            if (is16bit) {
                inputPtr->lumaExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize;
                inputPtr->cbExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize + luma10bitSize;
                inputPtr->crExt = config->sequenceBuffer[config->processedFrameCount % config->bufferedInput] + luma8bitSize + 2 * chroma8bitSize + luma10bitSize + chroma10bitSize;
                inputPtr->lumaExt = inputPtr->lumaExt + config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING;
                inputPtr->cbExt = inputPtr->cbExt + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);
                inputPtr->crExt = inputPtr->crExt + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);

            }

            headerPtr->n_filled_len = (uint32_t)(uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(inputPaddedWidth, inputPaddedHeight, is16bit);

        }


    }

    // If we reached the end of file, loop over again
    if (feof(inputFile) != 0) {
        fseek(inputFile, 0, SEEK_SET);
    }

    return;
}

void SendQpOnTheFly(
    EbConfig_t                  *config,
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
                fseek(config->qpFile, 0, SEEK_SET);

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
APPEXITCONDITIONTYPE ProcessInputBuffer(
    EbConfig_t             *config,
    EbAppContext_t         *appCallBack)
{
    uint8_t                 is16bit = (uint8_t)(config->encoderBitDepth > 8);
    EbBufferHeaderType     *headerPtr = appCallBack->inputBufferPool;
    EbComponentType        *componentHandle = (EbComponentType*)appCallBack->svtEncoderHandle;

    APPEXITCONDITIONTYPE    return_value = APP_ExitConditionNone;

    int64_t                  inputPaddedWidth           = config->inputPaddedWidth;
    int64_t                  inputPaddedHeight          = config->inputPaddedHeight;
    int64_t                  frames_to_be_encoded          = config->frames_to_be_encoded;
    uint64_t                 frameSize                  = (uint64_t)((inputPaddedWidth*inputPaddedHeight * 3) / 2 + (inputPaddedWidth / 4 * inputPaddedHeight * 3) / 2);
    int64_t                  totalBytesToProcessCount;
    int64_t                  remainingByteCount;

    if (config->injector && config->processedFrameCount)
    {
        EbInjector(config->processedFrameCount, config->injector_frame_rate);
    }

    totalBytesToProcessCount = (frames_to_be_encoded < 0) ? -1 : (config->encoderBitDepth == 10 && config->compressedTenBitFormat == 1) ?
            frames_to_be_encoded * (int64_t)frameSize :
            frames_to_be_encoded * SIZE_OF_ONE_FRAME_IN_BYTES(inputPaddedWidth, inputPaddedHeight, is16bit);


    remainingByteCount       = (totalBytesToProcessCount < 0) ?   -1 :  totalBytesToProcessCount - (int64_t)config->processedByteCount;

    // If there are bytes left to encode, configure the header
    if (remainingByteCount != 0 && config->stopEncoder == EB_FALSE) {
        ReadInputFrames(
            config,
            is16bit,
            headerPtr);

        // Update the context parameters
        config->processedByteCount += headerPtr->n_filled_len;
        headerPtr->p_app_private          = (EB_PTR)EB_NULL;
        config->framesEncoded           = (int32_t)(++config->processedFrameCount);

        // Configuration parameters changed on the fly
        if (config->use_qp_file && config->qpFile)
            SendQpOnTheFly(
                config,
                headerPtr);

        if (keepRunning == 0 && !config->stopEncoder) {
            config->stopEncoder = EB_TRUE;
        }

        // Fill in Buffers Header control data
        headerPtr->pts          = config->processedFrameCount-1;
        headerPtr->pic_type    = EB_AV1_INVALID_PICTURE;

        headerPtr->flags = 0;

        // Send the picture
        eb_svt_enc_send_picture(componentHandle, headerPtr);

        if ((config->processedFrameCount == (uint64_t)config->frames_to_be_encoded) || config->stopEncoder) {

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

static void write_ivf_stream_header(EbConfig_t *config)
{
    char header[IVF_STREAM_HEADER_SIZE];
    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4, 0);                     // version
    mem_put_le16(header + 6, 32);                    // header size
    mem_put_le32(header + 8, AV1_FOURCC);                // fourcc
    mem_put_le16(header + 12, config->inputPaddedWidth);  // width
    mem_put_le16(header + 14, config->inputPaddedHeight); // height
    if (config->frameRateDenominator != 0 && config->frameRateNumerator != 0){
        mem_put_le32(header + 16, config->frameRateNumerator);  // rate
        mem_put_le32(header + 20, config->frameRateDenominator);            // scale
                                                    //mem_put_le32(header + 16, config->frameRateDenominator);  // rate
                                                    //mem_put_le32(header + 20, config->frameRateNumerator);  // scale
    }
    else {
        mem_put_le32(header + 16, (config->frameRate >> 16) * 1000);  // rate
        mem_put_le32(header + 20, 1000);            // scale
                                                    //mem_put_le32(header + 16, config->frameRateDenominator);  // rate
                                                    //mem_put_le32(header + 20, config->frameRateNumerator);  // scale
    }
    mem_put_le32(header + 24, 0);               // length
    mem_put_le32(header + 28, 0);               // unused
    //config->performanceContext.byteCount += 32;
    if (config->bitstreamFile)
        fwrite(header, 1, IVF_STREAM_HEADER_SIZE, config->bitstreamFile);

    return;
}
static void update_prev_ivf_header(EbConfig_t *config){

    char header[4]; // only for the number of bytes
    if (config && config->bitstreamFile && config->byte_count_since_ivf != 0){
        fseeko64(config->bitstreamFile, (-(int32_t)(config->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE)),SEEK_CUR);
        mem_put_le32(&header[0], (int32_t)(config->byte_count_since_ivf));
        fwrite(header, 1, 4, config->bitstreamFile);
        fseeko64(config->bitstreamFile, (config->byte_count_since_ivf + IVF_FRAME_HEADER_SIZE - 4), SEEK_CUR);
        config->byte_count_since_ivf = 0;
    }
}

static void write_ivf_frame_header(EbConfig_t *config, uint32_t byte_count){
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

    if (config->bitstreamFile)
        fwrite(header, 1, IVF_FRAME_HEADER_SIZE, config->bitstreamFile);
}

APPEXITCONDITIONTYPE ProcessOutputStreamBuffer(
    EbConfig_t             *config,
    EbAppContext_t         *appCallBack,
    uint8_t                 pic_send_done)
{
    APPPORTACTIVETYPE      *portState       = &appCallBack->outputStreamPortActive;
    EbBufferHeaderType     *headerPtr;
    EbComponentType        *componentHandle = (EbComponentType*)appCallBack->svtEncoderHandle;
    APPEXITCONDITIONTYPE    return_value    = APP_ExitConditionNone;
    EbErrorType             stream_status   = EB_ErrorNone;
    // Per channel variables
    FILE                   *streamFile       = config->bitstreamFile;

    uint64_t               *totalLatency     = &config->performanceContext.totalLatency;
    uint32_t               *maxLatency       = &config->performanceContext.maxLatency;

    // System performance variables
    static int32_t         frameCount                = 0;

    // Local variables
    uint64_t                finishsTime     = 0;
    uint64_t                finishuTime     = 0;

    // non-blocking call until all input frames are sent
    stream_status = eb_svt_get_packet(componentHandle, &headerPtr, pic_send_done);

    if (stream_status == EB_ErrorMax) {
        printf("\n");
        LogErrorOutput(
            config->errorLogFile,
            headerPtr->flags);
        return APP_ExitConditionError;
    }
    else if (stream_status != EB_NoErrorEmptyQueue) {
#if TILES
        EbBool   has_tiles                = (EbBool)(appCallBack->ebEncParameters.tile_columns || appCallBack->ebEncParameters.tile_rows);
#else
        EbBool   has_tiles                = (EbBool)EB_FALSE;
#endif
        uint8_t  obu_frame_header_size    = has_tiles ? OBU_FRAME_HEADER_SIZE + 1 : OBU_FRAME_HEADER_SIZE;
        ++(config->performanceContext.frameCount);
        *totalLatency += (uint64_t)headerPtr->n_tick_count;
        *maxLatency = (headerPtr->n_tick_count > *maxLatency) ? headerPtr->n_tick_count : *maxLatency;

        EbFinishTime((uint64_t*)&finishsTime, (uint64_t*)&finishuTime);

        // total execution time, inc init time
        EbComputeOverallElapsedTime(
            config->performanceContext.lib_start_time[0],
            config->performanceContext.lib_start_time[1],
            finishsTime,
            finishuTime,
            &config->performanceContext.total_execution_time);

        // total encode time
        EbComputeOverallElapsedTime(
            config->performanceContext.encode_start_time[0],
            config->performanceContext.encode_start_time[1],
            finishsTime,
            finishuTime,
            &config->performanceContext.total_encode_time);

        // Write Stream Data to file
        if (streamFile) {
            if (config->performanceContext.frameCount == 1){
                write_ivf_stream_header(config);
            }

            switch(headerPtr->flags & 0x00000006){ // Check for the flags EB_BUFFERFLAG_HAS_TD and EB_BUFFERFLAG_SHOW_EXT

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
        config->performanceContext.byteCount += headerPtr->n_filled_len;

        // Update Output Port Activity State
        *portState = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_PortInactive : *portState;
        return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;

        // Release the output buffer
        eb_svt_release_out_buffer(&headerPtr);

#if DEADLOCK_DEBUG
        ++frameCount;
#else
        //++frameCount;
        printf("\b\b\b\b\b\b\b\b\b%9d", ++frameCount);
#endif

        //++frameCount;
        fflush(stdout);

        {
            config->performanceContext.averageSpeed = (config->performanceContext.frameCount) / config->performanceContext.total_encode_time;
            config->performanceContext.averageLatency = config->performanceContext.totalLatency / (double)(config->performanceContext.frameCount);
        }

        if (!(frameCount % SPEED_MEASUREMENT_INTERVAL)) {
            {
                printf("\n");
                printf("Average System Encoding Speed:        %.2f\n", (double)(frameCount) / config->performanceContext.total_encode_time);
            }
        }
    }
    return return_value;
}
APPEXITCONDITIONTYPE ProcessOutputReconBuffer(
    EbConfig_t             *config,
    EbAppContext_t         *appCallBack)
{
    EbBufferHeaderType    *headerPtr = appCallBack->recon_buffer; // needs to change for buffered input
    EbComponentType       *componentHandle = (EbComponentType*)appCallBack->svtEncoderHandle;
    APPEXITCONDITIONTYPE    return_value = APP_ExitConditionNone;
    EbErrorType            recon_status = EB_ErrorNone;
    int32_t fseekReturnVal;
    // non-blocking call until all input frames are sent
    recon_status = eb_svt_get_recon(componentHandle, headerPtr);

    if (recon_status == EB_ErrorMax) {
        printf("\n");
        LogErrorOutput(
            config->errorLogFile,
            headerPtr->flags);
        return APP_ExitConditionError;
    }
    else if (recon_status != EB_NoErrorEmptyQueue) {
        //Sets the File position to the beginning of the file.
        rewind(config->reconFile);
        uint64_t frameNum = headerPtr->pts;
        while (frameNum>0) {
            fseekReturnVal = fseeko64(config->reconFile, headerPtr->n_filled_len, SEEK_CUR);

            if (fseekReturnVal != 0) {
                printf("Error in fseeko64  returnVal %i\n", fseekReturnVal);
                return APP_ExitConditionError;
            }
            frameNum = frameNum - 1;
        }

        fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, config->reconFile);

        // Update Output Port Activity State
        return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;
    }
    return return_value;
}

