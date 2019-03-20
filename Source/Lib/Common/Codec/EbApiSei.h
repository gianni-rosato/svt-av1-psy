/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAPISEI_h
#define EBAPISEI_h

#include "EbDefinitions.h"

#define  MAX_DECODING_UNIT_COUNT                  64   // picture timing SEI

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    // These are temporary definitions, should be changed in the future according to user's requirements
#define MAX_CPB_COUNT    4

// User defined structures for passing data from application to the library should be added here

    typedef struct RegistedUserData_s {

        uint8_t   *userData;    // First byte is itu_t_t35_country_code.
                              // If itu_t_t35_country_code  ==  0xFF, second byte is itu_t_t35_country_code_extension_byte.
                              // the rest are the payloadByte
        uint32_t   userDataSize;

    } RegistedUserData_t;

    // SEI structures
    typedef struct AppHrdParameters_s {

        EbBool                            nalHrdParametersPresentFlag;
        EbBool                            vclHrdParametersPresentFlag;
        EbBool                            subPicCpbParamsPresentFlag;
                                          
        uint32_t                          tickDivisorMinus2;
        uint32_t                          duCpbRemovalDelayLengthMinus1;
                                          
        EbBool                            subPicCpbParamsPicTimingSeiFlag;

        uint32_t                          dpbOutputDelayDuLengthMinus1;

        uint32_t                          bitRateScale;
        uint32_t                          cpbSizeScale;
        uint32_t                          duCpbSizeScale;

        uint32_t                          initialCpbRemovalDelayLengthMinus1;
        uint32_t                          auCpbRemovalDelayLengthMinus1;
        uint32_t                          dpbOutputDelayLengthMinus1;

        EbBool                            fixedPicRateGeneralFlag[MAX_TEMPORAL_LAYERS];
        EbBool                            fixedPicRateWithinCvsFlag[MAX_TEMPORAL_LAYERS];

        uint32_t                          elementalDurationTcMinus1[MAX_TEMPORAL_LAYERS];

        EbBool                            lowDelayHrdFlag[MAX_TEMPORAL_LAYERS];

        uint32_t                          cpbCountMinus1[MAX_TEMPORAL_LAYERS];

        uint32_t                          bitRateValueMinus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          cpbSizeValueMinus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          bitRateDuValueMinus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];
        uint32_t                          cpbSizeDuValueMinus1[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];

        EbBool                            cbrFlag[MAX_TEMPORAL_LAYERS][2][MAX_CPB_COUNT];

        EbBool                            cpbDpbDelaysPresentFlag;

    } AppHrdParameters_t;

    typedef struct AppVideoUsabilityInfo_s {

        EbBool                    aspectRatioInfoPresentFlag;
        uint32_t                  aspectRatioIdc;
        uint32_t                  sarWidth;
        uint32_t                  sarHeight;

        EbBool                    overscanInfoPresentFlag;
        EbBool                    overscanApproriateFlag;
        EbBool                    videoSignalTypePresentFlag;

        uint32_t                  videoFormat;
        EbBool                    videoFullRangeFlag;
        EbBool                    colorDescriptionPresentFlag;

        uint32_t                  colorPrimaries;
        uint32_t                  transferCharacteristics;
        uint32_t                  matrixCoeffs;

        EbBool                    chromaLocInfoPresentFlag;
        uint32_t                  chromaSampleLocTypeTopField;
        uint32_t                  chromaSampleLocTypeBottomField;


        EbBool                    neutralChromaIndicationFlag;
        EbBool                    fieldSeqFlag;
        EbBool                    frameFieldInfoPresentFlag;

        EbBool                    defaultDisplayWindowFlag;
        uint32_t                  defaultDisplayWinLeftOffset;
        uint32_t                  defaultDisplayWinRightOffset;
        uint32_t                  defaultDisplayWinTopOffset;
        uint32_t                  defaultDisplayWinBottomOffset;

        EbBool                    vuiTimingInfoPresentFlag;
        uint32_t                  vuiNumUnitsInTick;
        uint32_t                  vuiTimeScale;

        EbBool                    vuiPocPropotionalTimingFlag;
        uint32_t                  vuiNumTicksPocDiffOneMinus1;

        EbBool                    vuiHrdParametersPresentFlag;
                                  
        EbBool                    bitstreamRestrictionFlag;
                                  
        EbBool                    motionVectorsOverPicBoundariesFlag;
        EbBool                    restrictedRefPicListsFlag;

        uint32_t                  minSpatialSegmentationIdc;
        uint32_t                  maxBytesPerPicDenom;
        uint32_t                  maxBitsPerMinCuDenom;
        uint32_t                  log2MaxMvLengthHorizontal;
        uint32_t                  log2MaxMvLengthVertical;

        AppHrdParameters_t     *hrdParametersPtr;

    } AppVideoUsabilityInfo_t;


    typedef struct AppPictureTimingSei_s {

        uint32_t   picStruct;
        uint32_t   sourceScanType;
        EbBool  duplicateFlag;
        uint32_t   auCpbRemovalDelayMinus1;
        uint32_t   picDpbOutputDelay;
        uint32_t   picDpbOutputDuDelay;
        uint32_t   numDecodingUnitsMinus1;
        EbBool  duCommonCpbRemovalDelayFlag;
        uint32_t   duCommonCpbRemovalDelayMinus1;
        uint32_t   numNalusInDuMinus1;
        uint32_t   duCpbRemovalDelayMinus1[MAX_DECODING_UNIT_COUNT];

    } AppPictureTimingSei_t;

    typedef struct AppBufferingPeriodSei_s {

        uint32_t  bpSeqParameterSetId;
        EbBool rapCpbParamsPresentFlag;
        EbBool concatenationFlag;
        uint32_t  auCpbRemovalDelayDeltaMinus1;
        uint32_t  cpbDelayOffset;
        uint32_t  dpbDelayOffset;
        uint32_t  initialCpbRemovalDelay[2][MAX_CPB_COUNT];
        uint32_t  initialCpbRemovalDelayOffset[2][MAX_CPB_COUNT];
        uint32_t  initialAltCpbRemovalDelay[2][MAX_CPB_COUNT];
        uint32_t  initialAltCpbRemovalDelayOffset[2][MAX_CPB_COUNT];

    } AppBufferingPeriodSei_t;


    typedef struct AppRecoveryPoint_s {

        uint32_t  recoveryPocCnt;
        EbBool exactMatchingFlag;
        EbBool brokenLinkFlag;

    } AppRecoveryPoint_t;

    // Below is an example of PanScanRectangle SEI data structure
    // Other SEI messages can have data structure in this format
    typedef struct AppPanScanRectangleSei_s {

        uint32_t      panScanRectId;
        EbBool     panScanRectCancelFlag;

        uint32_t      panScanCountMinus1;
        uint32_t      panScanRectLeftOffset[3];
        uint32_t      panScanRectRightOffset[3];
        uint32_t      panScanRectTopOffset[3];
        uint32_t      panScanRectBottomOffset[3];

        EbBool     panScanRectPersistFlag;

    }AppPanScanRectangleSei_t;

    //


    typedef struct EB_FRAME_RATE_CFG
    {
        uint32_t num;
        uint32_t den;
    } EB_FRAME_RATE_CFG;

    typedef struct EB_LATENCY_CALC
    {
        uint64_t    startTimesSeconds;
        uint64_t    finishTimesSeconds;
        uint64_t    startTimesuSeconds;
        uint64_t    finishTimesuSeconds;
        uint64_t  poc;
    } EB_LATENCY_CALC;


    typedef struct EB_CU_QP_CFG_s
    {
        uint32_t      qpArrayCount;
        int8_t       *qp_array;

        uint32_t      cuStride;
        uint32_t      sb_count;
        uint32_t      sb_stride;

        uint32_t      numCUsPerLastRowLCU;
        uint32_t      numCUsPerTypicalLCU;

        // Hold the DifCUDeltaQpDepth in this struct so that it is simpler to access it in the Lib,
        // rather than having to reference something deep

    } EB_CU_QP_CFG_t;



    // Signals that the default prediction structure and controls are to be
    //   overwritten and manually controlled. Manual control should be active
    //   for an entire encode, from beginning to termination.  Mixing of default
    //   prediction structure control and override prediction structure control
    //   is not supported.
    //
    // The ref_list0_count and ref_list1_count variables control the picture/slice type.
    //   I_SLICE: ref_list0_count == 0 && ref_list1_count == 0
    //   P_SLICE: ref_list0_count  > 0 && ref_list1_count == 0
    //   B_SLICE: ref_list0_count  > 0 && ref_list1_count  > 0

    typedef struct EB_PRED_STRUCTURE_CFG
    {
        uint64_t              picture_number;      // Corresponds to the display order of the picture
        uint64_t              decodeOrderNumber;  // Corresponds to the decode order of the picture
        uint32_t              temporal_layer_index; // Corresponds to the temporal layer index of the picture
        uint32_t              nalUnitType;        // Pictures NAL Unit Type
        uint32_t              ref_list0_count;      // A count of zero indicates the list is inactive
        uint32_t              ref_list1_count;      // A count of zero indicates the list is inactive
        EbBool             isReferenced;       // Indicates whether or not the picture is used as
                                                //   future reference.
        uint32_t              futureReferenceCount;
        int32_t             *futureReferenceList;// Contains a list of delta POCs whose references shall
                                                //   be saved for future reference.  This signalling must
                                                //   be done with respect to decode picture order. Must
                                                //   be conformant with the DPB rules.
    } EB_PRED_STRUCTURE_CFG;

    // EB_PICTURE_PLANE defines the data formatting of a singple plane of picture data.
    typedef struct EB_PICTURE_PLANE
    {
        // "start" is the starting position of the first
        //   valid pixel in the picture plane.
        uint8_t* start;

        // "stride" is the number of bytes required to increment
        //   an address by one line without changing the horizontal
        //   positioning.
        uint32_t stride;

    } EB_PICTURE_PLANE;



    // EB_INPUT_PICTURE_DEF defines the data formatting of an input picture. Note, each
    //   element can change independently of the previously coded pictures.  This allows
    //   for sub-picture coding and de-interlacing.
    typedef struct EB_INPUT_PICTURE_DEF
    {
        uint32_t width;        // Y plane width  (in units of samples)
        uint32_t height;       // Y plane height (in units of samples)

        // Padding (in lines/pixels of Y plane)
        uint32_t top_padding;
        uint32_t bot_padding;
        uint32_t left_padding;
        uint32_t right_padding;

        // Y, Cb, Cr planes. Note, for bitdepths greater than 8 using the "unpacked" format,
        //  the following elements contain the 8 MSB bits of the sample.
        EB_PICTURE_PLANE yPlane;
        EB_PICTURE_PLANE cbPlane;
        EB_PICTURE_PLANE crPlane;

        // Auxiliary planes used when bitdepths are greater than 8 and using the "unpacked"
        //   format.for The LSB bits of the sample should be located at the MSB bit positions
        //   of each byte. Must be set to NULL for 8-bit input or "packed" 10 or 12-bit input.
        EB_PICTURE_PLANE yAuxPlane;
        EB_PICTURE_PLANE cbAuxPlane;
        EB_PICTURE_PLANE crAuxPlane;

    } EB_INPUT_PICTURE_DEF;

    typedef struct EB_EOS_DATA_DEF
    {
        EbBool             dataValid;          // Indicates if the data attached with the last frame in the bitstream is valid or not.
                                                //   If the last frame is valid, the data will be included in the bistream
                                                //   If the last frame is NOT valid, the frame will be coded as IDR in the encoder, but
                                                //   not included in the bitstream.

    } EB_EOS_DATA_DEF;




    extern EbErrorType EbAppVideoUsabilityInfoCtor(
        AppVideoUsabilityInfo_t *vui_ptr);





#ifdef __cplusplus
}
#endif // __cplusplus
#endif // EBAPISEI_h
