/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"

#include "EbPictureManagerProcess.h"
#include "EbReferenceObject.h"

#include "EbPictureDemuxResults.h"
#include "EbPictureManagerQueue.h"
#include "EbPredictionStructure.h"
#include "EbRateControlTasks.h"
#include "EbSvtAv1ErrorCodes.h"

#if TILES
void av1_tile_set_col(TileInfo *tile, PictureParentControlSet_t * pcsPtr, int col);
void av1_tile_set_row(TileInfo *tile, PictureParentControlSet_t * pcsPtr, int row);
void set_tile_info(PictureParentControlSet_t * pcsPtr);
#endif

/************************************************
 * Defines
 ************************************************/
#define POC_CIRCULAR_ADD(base, offset)             (((base) + (offset)))

 /************************************************
  * Configure Picture edges
  ************************************************/
static void ConfigurePictureEdges(
    SequenceControlSet_t *scsPtr,
    PictureControlSet_t  *ppsPtr)
{
    // Tiles Initialisation
    const uint16_t picture_width_in_sb = (scsPtr->luma_width + scsPtr->sb_size_pix - 1) / scsPtr->sb_size_pix;
    const uint16_t picture_height_in_sb = (scsPtr->luma_height + scsPtr->sb_size_pix - 1) / scsPtr->sb_size_pix;
    unsigned xLcuIndex, yLcuIndex, sb_index;

    // LCU-loops
    for (yLcuIndex = 0; yLcuIndex < picture_height_in_sb; ++yLcuIndex) {
        for (xLcuIndex = 0; xLcuIndex < picture_width_in_sb; ++xLcuIndex) {
            sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);
            ppsPtr->sb_ptr_array[sb_index]->picture_left_edge_flag = (xLcuIndex == 0) ? EB_TRUE : EB_FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_top_edge_flag = (yLcuIndex == 0) ? EB_TRUE : EB_FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_right_edge_flag = (xLcuIndex == (unsigned)(picture_width_in_sb - 1)) ? EB_TRUE : EB_FALSE;
        }
    }

    return;
}

/************************************************
 * Picture Manager Context Constructor
 ************************************************/
EbErrorType picture_manager_context_ctor(
    PictureManagerContext_t **context_dbl_ptr,
    EbFifo_t                 *picture_input_fifo_ptr,
    EbFifo_t                 *pictureManagerOutputFifoPtr,
    EbFifo_t                **picture_control_set_fifo_ptr_array)
{
    PictureManagerContext_t *context_ptr;
    EB_MALLOC(PictureManagerContext_t*, context_ptr, sizeof(PictureManagerContext_t), EB_N_PTR);

    *context_dbl_ptr = context_ptr;

    context_ptr->picture_input_fifo_ptr = picture_input_fifo_ptr;
    context_ptr->pictureManagerOutputFifoPtr = pictureManagerOutputFifoPtr;
    context_ptr->picture_control_set_fifo_ptr_array = picture_control_set_fifo_ptr_array;

    return EB_ErrorNone;
}




/***************************************************************************************************
 * Picture Manager Kernel
 *
 * Notes on the Picture Manager:
 *
 * The Picture Manager Process performs the function of managing both the Input Picture buffers and
 * the Reference Picture buffers and subdividing the Input Picture into Tiles. Both the Input Picture
 * and Reference Picture buffers particular management depends on the GoP structure already implemented in
 * the Picture decision. Also note that the Picture Manager sets up the RPS for Entropy Coding as well.
 *
 * Inputs:
 * Input Picture
 *   -Input Picture Data
 *
 *  Reference Picture
 *   -Reference Picture Data
 *
 *  Outputs:
 *   -Picture Control Set with fully available Reference List
 *
 ***************************************************************************************************/
void* picture_manager_kernel(void *input_ptr)
{
    PictureManagerContext_t         *context_ptr = (PictureManagerContext_t*)input_ptr;

    EbObjectWrapper_t               *ChildPictureControlSetWrapperPtr;
    PictureControlSet_t             *ChildPictureControlSetPtr;
    PictureParentControlSet_t       *picture_control_set_ptr;
    SequenceControlSet_t            *sequence_control_set_ptr;
    EncodeContext_t                 *encode_context_ptr;


    EbObjectWrapper_t               *inputPictureDemuxWrapperPtr;
    PictureDemuxResults_t           *inputPictureDemuxPtr;

    EbObjectWrapper_t               *outputWrapperPtr;
    RateControlTasks_t              *rateControlTasksPtr;

    EbBool                           availabilityFlag;

    PredictionStructureEntry_t      *predPositionPtr;
    InputQueueEntry_t               *inputEntryPtr;
    uint32_t                         inputQueueIndex;
    uint64_t                         current_input_poc;
    ReferenceQueueEntry_t           *referenceEntryPtr;
    uint32_t                         referenceQueueIndex;
    uint64_t                         refPoc;
    uint32_t                         depIdx;
    uint64_t                         depPoc;
    uint32_t                         depListCount;
    PictureParentControlSet_t       *entryPictureControlSetPtr;
    SequenceControlSet_t            *entrySequenceControlSetPtr;

    // Initialization
    uint8_t                          picture_width_in_sb;
    uint8_t                          picture_height_in_sb;
    PictureManagerReorderEntry_t    *queueEntryPtr;
    int32_t                          queueEntryIndex;

    // Debug
    uint32_t loopCount = 0;

    for (;;) {

        // Get Input Full Object
        eb_get_full_object(
            context_ptr->picture_input_fifo_ptr,
            &inputPictureDemuxWrapperPtr);

        inputPictureDemuxPtr = (PictureDemuxResults_t*)inputPictureDemuxWrapperPtr->object_ptr;

        // *Note - This should be overhauled and/or replaced when we
        //   need hierarchical support.
        loopCount++;

        switch (inputPictureDemuxPtr->pictureType) {

        case EB_PIC_INPUT:

            picture_control_set_ptr = (PictureParentControlSet_t*)inputPictureDemuxPtr->pictureControlSetWrapperPtr->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

            //printf("\nPicture Manager Process @ %d \n ", picture_control_set_ptr->picture_number);

            queueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index]->picture_number);
            queueEntryIndex += encode_context_ptr->picture_manager_reorder_queue_head_index;
            queueEntryIndex = (queueEntryIndex > PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH : queueEntryIndex;
            queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[queueEntryIndex];
            if (queueEntryPtr->parentPcsWrapperPtr != NULL) {
                CHECK_REPORT_ERROR_NC(
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PD_ERROR8);
            }
            else {
                queueEntryPtr->parentPcsWrapperPtr = inputPictureDemuxPtr->pictureControlSetWrapperPtr;
                queueEntryPtr->picture_number = picture_control_set_ptr->picture_number;
            }
            // Process the head of the Picture Manager Reorder Queue
            queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index];

            while (queueEntryPtr->parentPcsWrapperPtr != EB_NULL) {

                picture_control_set_ptr = (PictureParentControlSet_t*)queueEntryPtr->parentPcsWrapperPtr->object_ptr;

                predPositionPtr = picture_control_set_ptr->pred_struct_ptr->predStructEntryPtrArray[picture_control_set_ptr->pred_struct_index];
#if NEW_PRED_STRUCT
                // If there was a change in the number of temporal layers, then cleanup the Reference Queue's Dependent Counts
                if (picture_control_set_ptr->hierarchical_layers_diff != 0) {

                    // Dynamic GOP
                    PredictionStructure_t          *next_pred_struct_ptr;
                    PredictionStructureEntry_t     *next_base_layer_pred_position_ptr;
                    
                    uint32_t                        dependant_list_positive_entries;
                    uint32_t                        dependant_list_removed_entries;

                    referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;

                    while (referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) {

                        referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                        if (referenceEntryPtr->picture_number == (picture_control_set_ptr->picture_number - 1)) { // Picture where the change happened 

                            // Get the prediction struct entry of the next GOP structure
                            next_pred_struct_ptr = GetPredictionStructure(
                                encode_context_ptr->prediction_structure_group_ptr,
                                picture_control_set_ptr->pred_structure,
                                1,
                                picture_control_set_ptr->hierarchical_levels);

                            // Get the prediction struct of a picture in temporal layer 0 (from the new GOP structure)
                            next_base_layer_pred_position_ptr = next_pred_struct_ptr->predStructEntryPtrArray[next_pred_struct_ptr->predStructEntryCount - 1];


                            // Remove all positive entries from the dependant lists
                            dependant_list_positive_entries = 0;
                            for (depIdx = 0; depIdx < referenceEntryPtr->list0.listCount; ++depIdx) {
                                if (referenceEntryPtr->list0.list[depIdx] >= 0) {
                                    dependant_list_positive_entries++;
                                }
                            }
                            referenceEntryPtr->list0.listCount = referenceEntryPtr->list0.listCount - dependant_list_positive_entries;

                            dependant_list_positive_entries = 0;
                            for (depIdx = 0; depIdx < referenceEntryPtr->list1.listCount; ++depIdx) {
                                if (referenceEntryPtr->list1.list[depIdx] >= 0) {
                                    dependant_list_positive_entries++;
                                }
                            }
                            referenceEntryPtr->list1.listCount = referenceEntryPtr->list1.listCount - dependant_list_positive_entries;

                            for (depIdx = 0; depIdx < next_base_layer_pred_position_ptr->depList0.listCount; ++depIdx) {
                                if (next_base_layer_pred_position_ptr->depList0.list[depIdx] >= 0) {
                                    referenceEntryPtr->list0.list[referenceEntryPtr->list0.listCount++] = next_base_layer_pred_position_ptr->depList0.list[depIdx];
                                }
                            }


                            for (depIdx = 0; depIdx < next_base_layer_pred_position_ptr->depList1.listCount; ++depIdx) {
                                if (next_base_layer_pred_position_ptr->depList1.list[depIdx] >= 0) {
                                    referenceEntryPtr->list1.list[referenceEntryPtr->list1.listCount++] = next_base_layer_pred_position_ptr->depList1.list[depIdx];
                                }
                            }

                            // Update the dependant count update
                            dependant_list_removed_entries = referenceEntryPtr->depList0Count + referenceEntryPtr->depList1Count - referenceEntryPtr->dependentCount;

                            referenceEntryPtr->depList0Count = referenceEntryPtr->list0.listCount;
                            referenceEntryPtr->depList1Count = referenceEntryPtr->list1.listCount;
                            referenceEntryPtr->dependentCount = referenceEntryPtr->depList0Count + referenceEntryPtr->depList1Count - dependant_list_removed_entries;

                        }
                        else {

                            // Modify Dependent List0
                            depListCount = referenceEntryPtr->list0.listCount;
                            for (depIdx = 0; depIdx < depListCount; ++depIdx) {


                                // Adjust the latest currentInputPoc in case we're in a POC rollover scenario 
                                // currentInputPoc += (currentInputPoc < referenceEntryPtr->pocNumber) ? (1 << sequence_control_set_ptr->bitsForPictureOrderCount) : 0;

                                depPoc = POC_CIRCULAR_ADD(
                                    referenceEntryPtr->picture_number, // can't use a value that gets reset
                                    referenceEntryPtr->list0.list[depIdx]/*,
                                    sequence_control_set_ptr->bitsForPictureOrderCount*/);

                                    // If Dependent POC is greater or equal to the IDR POC
                                if (depPoc >= picture_control_set_ptr->picture_number && referenceEntryPtr->list0.list[depIdx]) {

                                    referenceEntryPtr->list0.list[depIdx] = 0;

                                    // Decrement the Reference's referenceCount
                                    --referenceEntryPtr->dependentCount;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependentCount != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR3);
                                }
                            }

                            // Modify Dependent List1
                            depListCount = referenceEntryPtr->list1.listCount;
                            for (depIdx = 0; depIdx < depListCount; ++depIdx) {

                                // Adjust the latest currentInputPoc in case we're in a POC rollover scenario 
                                // currentInputPoc += (currentInputPoc < referenceEntryPtr->pocNumber) ? (1 << sequence_control_set_ptr->bitsForPictureOrderCount) : 0;

                                depPoc = POC_CIRCULAR_ADD(
                                    referenceEntryPtr->picture_number,
                                    referenceEntryPtr->list1.list[depIdx]/*,
                                    sequence_control_set_ptr->bitsForPictureOrderCount*/);

                                    // If Dependent POC is greater or equal to the IDR POC
                                if ((depPoc >= picture_control_set_ptr->picture_number) && referenceEntryPtr->list1.list[depIdx]) {
                                    referenceEntryPtr->list1.list[depIdx] = 0;

                                    // Decrement the Reference's referenceCount
                                    --referenceEntryPtr->dependentCount;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependentCount != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR3);
                                }
                            }
                        }

                        // Increment the referenceQueueIndex Iterator
                        referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;
                    }
                }
#endif
                // If there was an I-frame or Scene Change, then cleanup the Reference Queue's Dependent Counts
                if (picture_control_set_ptr->slice_type == I_SLICE)
                {

                    referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;
                    while (referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) {

                        referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                        // Modify Dependent List0
                        depListCount = referenceEntryPtr->list0.listCount;
                        for (depIdx = 0; depIdx < depListCount; ++depIdx) {

                            current_input_poc = picture_control_set_ptr->picture_number;

                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                            // current_input_poc += (current_input_poc < referenceEntryPtr->picture_number) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                            depPoc = POC_CIRCULAR_ADD(
                                referenceEntryPtr->picture_number, // can't use a value that gets reset
                                referenceEntryPtr->list0.list[depIdx]/*,
                                sequence_control_set_ptr->bits_for_picture_order_count*/);

                                // If Dependent POC is greater or equal to the IDR POC
                            if (depPoc >= current_input_poc && referenceEntryPtr->list0.list[depIdx]) {

                                referenceEntryPtr->list0.list[depIdx] = 0;

                                // Decrement the Reference's referenceCount
                                --referenceEntryPtr->dependentCount;
                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependentCount != ~0u),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR3);
                            }
                        }

                        // Modify Dependent List1
                        depListCount = referenceEntryPtr->list1.listCount;
                        for (depIdx = 0; depIdx < depListCount; ++depIdx) {

                            current_input_poc = picture_control_set_ptr->picture_number;

                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                            // current_input_poc += (current_input_poc < referenceEntryPtr->picture_number) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                            depPoc = POC_CIRCULAR_ADD(
                                referenceEntryPtr->picture_number,
                                referenceEntryPtr->list1.list[depIdx]/*,
                                sequence_control_set_ptr->bits_for_picture_order_count*/);

                                // If Dependent POC is greater or equal to the IDR POC or if we inserted trailing Ps
                            if (((depPoc >= current_input_poc) || (((picture_control_set_ptr->pre_assignment_buffer_count != picture_control_set_ptr->pred_struct_ptr->predStructPeriod) || (picture_control_set_ptr->idr_flag == EB_TRUE)) && (depPoc > (current_input_poc - picture_control_set_ptr->pre_assignment_buffer_count)))) && referenceEntryPtr->list1.list[depIdx]) {

                                referenceEntryPtr->list1.list[depIdx] = 0;

                                // Decrement the Reference's referenceCount
                                --referenceEntryPtr->dependentCount;
                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependentCount != ~0u),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR3);
                            }

                        }

                        // Increment the referenceQueueIndex Iterator
                        referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;
                    }

                }
                else if (picture_control_set_ptr->idr_flag == EB_TRUE) {

                    // Set Reference Entry pointer
                    referenceEntryPtr = (ReferenceQueueEntry_t*)EB_NULL;
                }

                // Check if the EnhancedPictureQueue is full.
                // *Note - Having the number of Enhanced Pictures less than the InputQueueSize should ensure this never gets hit

                CHECK_REPORT_ERROR(
                    (((encode_context_ptr->input_picture_queue_head_index != encode_context_ptr->input_picture_queue_tail_index) || (encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]->inputObjectPtr == EB_NULL))),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR4);

                // Place Picture in input queue
                inputEntryPtr = encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_tail_index];
                inputEntryPtr->inputObjectPtr = queueEntryPtr->parentPcsWrapperPtr;
                inputEntryPtr->referenceEntryIndex = encode_context_ptr->reference_picture_queue_tail_index;
                encode_context_ptr->input_picture_queue_tail_index =
                    (encode_context_ptr->input_picture_queue_tail_index == INPUT_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->input_picture_queue_tail_index + 1;

                // Copy the reference lists into the inputEntry and
                // set the Reference Counts Based on Temporal Layer and how many frames are active
                picture_control_set_ptr->ref_list0_count = (picture_control_set_ptr->slice_type == I_SLICE) ? 0 : (uint8_t)predPositionPtr->refList0.referenceListCount;
                picture_control_set_ptr->ref_list1_count = (picture_control_set_ptr->slice_type == I_SLICE) ? 0 : (uint8_t)predPositionPtr->refList1.referenceListCount;
                inputEntryPtr->list0Ptr = &predPositionPtr->refList0;
                inputEntryPtr->list1Ptr = &predPositionPtr->refList1;

                // Check if the ReferencePictureQueue is full.
                CHECK_REPORT_ERROR(
                    (((encode_context_ptr->reference_picture_queue_head_index != encode_context_ptr->reference_picture_queue_tail_index) || (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->referenceObjectPtr == EB_NULL))),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR5);

                // Create Reference Queue Entry even if picture will not be referenced
                referenceEntryPtr = encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_tail_index];
                referenceEntryPtr->picture_number = picture_control_set_ptr->picture_number;
                referenceEntryPtr->referenceObjectPtr = (EbObjectWrapper_t*)EB_NULL;
                referenceEntryPtr->releaseEnable = EB_TRUE;
                referenceEntryPtr->referenceAvailable = EB_FALSE;
                referenceEntryPtr->is_used_as_reference_flag = picture_control_set_ptr->is_used_as_reference_flag;
                encode_context_ptr->reference_picture_queue_tail_index =
                    (encode_context_ptr->reference_picture_queue_tail_index == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->reference_picture_queue_tail_index + 1;

                // Copy the Dependent Lists
                // *Note - we are removing any leading picture dependencies for now
                referenceEntryPtr->list0.listCount = 0;
                for (depIdx = 0; depIdx < predPositionPtr->depList0.listCount; ++depIdx) {
                    if (predPositionPtr->depList0.list[depIdx] >= 0) {
                        referenceEntryPtr->list0.list[referenceEntryPtr->list0.listCount++] = predPositionPtr->depList0.list[depIdx];
                    }
                }

                referenceEntryPtr->list1.listCount = predPositionPtr->depList1.listCount;
                for (depIdx = 0; depIdx < predPositionPtr->depList1.listCount; ++depIdx) {
                    referenceEntryPtr->list1.list[depIdx] = predPositionPtr->depList1.list[depIdx];
                }

                referenceEntryPtr->depList0Count = referenceEntryPtr->list0.listCount;
                referenceEntryPtr->depList1Count = referenceEntryPtr->list1.listCount;
                referenceEntryPtr->dependentCount = referenceEntryPtr->depList0Count + referenceEntryPtr->depList1Count;

                CHECK_REPORT_ERROR(
                    (picture_control_set_ptr->pred_struct_ptr->predStructPeriod < MAX_ELAPSED_IDR_COUNT),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR6);

                // Release the Reference Buffer once we know it is not a reference
                if (picture_control_set_ptr->is_used_as_reference_flag == EB_FALSE) {
                    // Release the nominal liveCount value
                    eb_release_object(picture_control_set_ptr->reference_picture_wrapper_ptr);
                    picture_control_set_ptr->reference_picture_wrapper_ptr = (EbObjectWrapper_t*)EB_NULL;
                }

                // Release the Picture Manager Reorder Queue
                queueEntryPtr->parentPcsWrapperPtr = (EbObjectWrapper_t*)EB_NULL;
                queueEntryPtr->picture_number += PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH;

                // Increment the Picture Manager Reorder Queue
                encode_context_ptr->picture_manager_reorder_queue_head_index = (encode_context_ptr->picture_manager_reorder_queue_head_index == PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_manager_reorder_queue_head_index + 1;

                // Get the next entry from the Picture Manager Reorder Queue (Entry N+1)
                queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index];

            }
            break;

        case EB_PIC_REFERENCE:

            sequence_control_set_ptr = (SequenceControlSet_t*)inputPictureDemuxPtr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

            // Check if Reference Queue is full
            CHECK_REPORT_ERROR(
                (encode_context_ptr->reference_picture_queue_head_index != encode_context_ptr->reference_picture_queue_tail_index),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR7);

            referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;

            // Find the Reference in the Reference Queue
            do {

                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                if (referenceEntryPtr->picture_number == inputPictureDemuxPtr->picture_number) {

                    // Assign the reference object if there is a match
                    referenceEntryPtr->referenceObjectPtr = inputPictureDemuxPtr->reference_picture_wrapper_ptr;

                    // Set the reference availability
                    referenceEntryPtr->referenceAvailable = EB_TRUE;
                }

                // Increment the referenceQueueIndex Iterator
                referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;

            } while ((referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) && (referenceEntryPtr->picture_number != inputPictureDemuxPtr->picture_number));


            CHECK_REPORT_ERROR(
                (referenceEntryPtr->picture_number == inputPictureDemuxPtr->picture_number),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR8);

            //keep the relase of SCS here because we still need the encodeContext strucutre here
            // Release the Reference's SequenceControlSet
            eb_release_object(inputPictureDemuxPtr->sequence_control_set_wrapper_ptr);

            break;

        default:

            sequence_control_set_ptr = (SequenceControlSet_t*)inputPictureDemuxPtr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

            CHECK_REPORT_ERROR_NC(
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR9);

            picture_control_set_ptr = (PictureParentControlSet_t*)EB_NULL;
            encode_context_ptr = (EncodeContext_t*)EB_NULL;

            break;
        }

        // ***********************************
        //  Common Code
        // *************************************

        // Walk the input queue and start all ready pictures.  Mark entry as null after started.  Increment the head as you go.
        if (encode_context_ptr != (EncodeContext_t*)EB_NULL) {
            inputQueueIndex = encode_context_ptr->input_picture_queue_head_index;
            while (inputQueueIndex != encode_context_ptr->input_picture_queue_tail_index) {

                inputEntryPtr = encode_context_ptr->input_picture_queue[inputQueueIndex];

                if (inputEntryPtr->inputObjectPtr != EB_NULL) {

                    entryPictureControlSetPtr = (PictureParentControlSet_t*)inputEntryPtr->inputObjectPtr->object_ptr;
                    entrySequenceControlSetPtr = (SequenceControlSet_t*)entryPictureControlSetPtr->sequence_control_set_wrapper_ptr->object_ptr;

                    availabilityFlag = EB_TRUE;

                    // Check RefList0 Availability
                    if (entryPictureControlSetPtr->ref_list0_count) {
                        referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                            ((int32_t)inputEntryPtr->referenceEntryIndex) -     // Base
                            inputEntryPtr->list0Ptr->referenceList,     // Offset
                            REFERENCE_QUEUE_MAX_DEPTH);                         // Max

                        referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                        CHECK_REPORT_ERROR(
                            (referenceEntryPtr),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PM_ERROR10);

                        refPoc = POC_CIRCULAR_ADD(
                            entryPictureControlSetPtr->picture_number,
                            -inputEntryPtr->list0Ptr->referenceList/*,
                            entrySequenceControlSetPtr->bits_for_picture_order_count*/);

                            // Increment the current_input_poc is the case of POC rollover
                        current_input_poc = encode_context_ptr->current_input_poc;
                        //current_input_poc += ((current_input_poc < refPoc) && (inputEntryPtr->list0Ptr->referenceList[refIdx] > 0)) ?
                        //    (1 << entrySequenceControlSetPtr->bits_for_picture_order_count) :
                        //    0;

                        availabilityFlag =
                            (availabilityFlag == EB_FALSE) ? EB_FALSE :   // Don't update if already False
                            (refPoc > current_input_poc) ? EB_FALSE :   // The Reference has not been received as an Input Picture yet, then its availability is false
                            (referenceEntryPtr->referenceAvailable) ? EB_TRUE :   // The Reference has been completed
                            EB_FALSE;     // The Reference has not been completed
                    }

                    // Check RefList1 Availability
                    if (entryPictureControlSetPtr->slice_type == B_SLICE) {
                        if (entryPictureControlSetPtr->ref_list1_count) {
                            // If Reference is valid (non-zero), update the availability
                            if (inputEntryPtr->list1Ptr->referenceList != (int32_t)INVALID_POC) {

                                referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                    ((int32_t)inputEntryPtr->referenceEntryIndex) -     // Base
                                    inputEntryPtr->list1Ptr->referenceList,     // Offset
                                    REFERENCE_QUEUE_MAX_DEPTH);                         // Max

                                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR10);

                                refPoc = POC_CIRCULAR_ADD(
                                    entryPictureControlSetPtr->picture_number,
                                    -inputEntryPtr->list1Ptr->referenceList/*,
                                    entrySequenceControlSetPtr->bits_for_picture_order_count*/);

                                    // Increment the current_input_poc is the case of POC rollover
                                current_input_poc = encode_context_ptr->current_input_poc;
                                //current_input_poc += ((current_input_poc < refPoc && inputEntryPtr->list1Ptr->referenceList[refIdx] > 0)) ?
                                //    (1 << entrySequenceControlSetPtr->bits_for_picture_order_count) :
                                //    0;

                                availabilityFlag =
                                    (availabilityFlag == EB_FALSE) ? EB_FALSE :   // Don't update if already False
                                    (refPoc > current_input_poc) ? EB_FALSE :   // The Reference has not been received as an Input Picture yet, then its availability is false
                                    (referenceEntryPtr->referenceAvailable) ? EB_TRUE :   // The Reference has been completed
                                    EB_FALSE;     // The Reference has not been completed
                            }
                        }
                    }

                    if (availabilityFlag == EB_TRUE) {

                        // Get New  Empty Child PCS from PCS Pool
                        eb_get_empty_object(
                            context_ptr->picture_control_set_fifo_ptr_array[0],
                            &ChildPictureControlSetWrapperPtr);

                        // Child PCS is released by Packetization
                        eb_object_inc_live_count(
                            ChildPictureControlSetWrapperPtr,
                            1);

                        ChildPictureControlSetPtr = (PictureControlSet_t*)ChildPictureControlSetWrapperPtr->object_ptr;

                        //1.Link The Child PCS to its Parent
                        ChildPictureControlSetPtr->picture_parent_control_set_wrapper_ptr = inputEntryPtr->inputObjectPtr;
                        ChildPictureControlSetPtr->parent_pcs_ptr = entryPictureControlSetPtr;


                        ChildPictureControlSetPtr->parent_pcs_ptr->childPcs = ChildPictureControlSetPtr;


                        //2. Have some common information between  ChildPCS and ParentPCS.
                        ChildPictureControlSetPtr->sequence_control_set_wrapper_ptr = entryPictureControlSetPtr->sequence_control_set_wrapper_ptr;
                        ChildPictureControlSetPtr->picture_qp = entryPictureControlSetPtr->picture_qp;
                        ChildPictureControlSetPtr->picture_number = entryPictureControlSetPtr->picture_number;
                        ChildPictureControlSetPtr->slice_type = entryPictureControlSetPtr->slice_type;
                        ChildPictureControlSetPtr->temporal_layer_index = entryPictureControlSetPtr->temporal_layer_index;

                        ChildPictureControlSetPtr->parent_pcs_ptr->total_num_bits = 0;
                        ChildPictureControlSetPtr->parent_pcs_ptr->picture_qp = entryPictureControlSetPtr->picture_qp;
                        ChildPictureControlSetPtr->parent_pcs_ptr->sad_me = 0;
                        ChildPictureControlSetPtr->parent_pcs_ptr->quantized_coeff_num_bits = 0;
                        ChildPictureControlSetPtr->enc_mode = entryPictureControlSetPtr->enc_mode;

                        //3.make all  init for ChildPCS
                        picture_width_in_sb = (uint8_t)((entrySequenceControlSetPtr->luma_width + entrySequenceControlSetPtr->sb_size_pix - 1) / entrySequenceControlSetPtr->sb_size_pix);
                        picture_height_in_sb = (uint8_t)((entrySequenceControlSetPtr->luma_height + entrySequenceControlSetPtr->sb_size_pix - 1) / entrySequenceControlSetPtr->sb_size_pix);

                        // EncDec Segments
                        EncDecSegmentsInit(
                            ChildPictureControlSetPtr->enc_dec_segment_ctrl,
                            entrySequenceControlSetPtr->enc_dec_segment_col_count_array[entryPictureControlSetPtr->temporal_layer_index],
                            entrySequenceControlSetPtr->enc_dec_segment_row_count_array[entryPictureControlSetPtr->temporal_layer_index],
                            picture_width_in_sb,
                            picture_height_in_sb);

                        // Entropy Coding Rows
                        {
                            unsigned row_index;

                            ChildPictureControlSetPtr->entropy_coding_current_row = 0;
                            ChildPictureControlSetPtr->entropy_coding_current_available_row = 0;
                            ChildPictureControlSetPtr->entropy_coding_row_count = picture_height_in_sb;
                            ChildPictureControlSetPtr->entropy_coding_in_progress = EB_FALSE;

                            for (row_index = 0; row_index < MAX_LCU_ROWS; ++row_index) {
                                ChildPictureControlSetPtr->entropy_coding_row_array[row_index] = EB_FALSE;
                            }
                        }

#if TILES             
                        set_tile_info(ChildPictureControlSetPtr->parent_pcs_ptr);

                        struct PictureParentControlSet_s     *ppcs_ptr = ChildPictureControlSetPtr->parent_pcs_ptr;
                        Av1Common *const cm = ppcs_ptr->av1_cm;
                        int tile_row, tile_col;
                        uint32_t  x_lcu_index,  y_lcu_index;
                        const int tile_cols = ppcs_ptr->av1_cm->tile_cols;
                        const int tile_rows = ppcs_ptr->av1_cm->tile_rows;                        
                        TileInfo tile_info;
                        //Tile Loop
                        for (tile_row = 0; tile_row < tile_rows; tile_row++)
                        {                           
                            av1_tile_set_row(&tile_info, ppcs_ptr, tile_row);

                            for (tile_col = 0; tile_col < tile_cols; tile_col++)
                            {
                                av1_tile_set_col(&tile_info, ppcs_ptr, tile_col);

                                for (y_lcu_index = cm->tile_row_start_sb[tile_row]; y_lcu_index < (uint32_t)cm->tile_row_start_sb[tile_row + 1]; ++y_lcu_index)
                                {
                                    for (x_lcu_index = cm->tile_col_start_sb[tile_col]; x_lcu_index < (uint32_t)cm->tile_col_start_sb[tile_col + 1]; ++x_lcu_index)
                                    {
                                        int sb_index = (uint16_t)(x_lcu_index + y_lcu_index * picture_width_in_sb);
                                        ChildPictureControlSetPtr->sb_ptr_array[sb_index]->tile_info = tile_info;                                       
                                    }
                                }                               
                            }
                        }                       

#endif

                        // Picture edges
                        ConfigurePictureEdges(entrySequenceControlSetPtr, ChildPictureControlSetPtr);

                        // Reset the qp array for DLF
                        EB_MEMSET(ChildPictureControlSetPtr->qp_array, 0, sizeof(uint8_t)*ChildPictureControlSetPtr->qp_array_size);


                        // Error resilience related
                        ChildPictureControlSetPtr->colocated_pu_ref_list = REF_LIST_0;     // to be modified

                        ChildPictureControlSetPtr->is_low_delay = (EbBool)(
                            ChildPictureControlSetPtr->parent_pcs_ptr->pred_struct_ptr->predStructEntryPtrArray[ChildPictureControlSetPtr->parent_pcs_ptr->pred_struct_index]->positiveRefPicsTotalCount == 0);

                        // Rate Control
                        ChildPictureControlSetPtr->use_delta_qp = (uint8_t)entrySequenceControlSetPtr->static_config.improve_sharpness;
                        ChildPictureControlSetPtr->dif_cu_delta_qp_depth = (uint8_t)entrySequenceControlSetPtr->input_resolution == INPUT_SIZE_4K_RANGE ? 3 : 2;

                        // Reset the Reference Lists
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_ptr_array, 0, 2 * sizeof(EbObjectWrapper_t*));

                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_qp_array, 0, 2 * sizeof(uint8_t));

                        EB_MEMSET(ChildPictureControlSetPtr->ref_slice_type_array, 0, 2 * sizeof(EB_SLICE));

                        // Configure List0
                        if ((entryPictureControlSetPtr->slice_type == P_SLICE) || (entryPictureControlSetPtr->slice_type == B_SLICE)) {

                            if (entryPictureControlSetPtr->ref_list0_count) {
                                referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                    ((int32_t)inputEntryPtr->referenceEntryIndex) - inputEntryPtr->list0Ptr->referenceList,
                                    REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                // Set the Reference Object
                                ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0] = referenceEntryPtr->referenceObjectPtr;

#if ADD_DELTA_QP_SUPPORT
                                ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_0] = (uint8_t)((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->qp;
                                ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_0] = (uint8_t)((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->slice_type;
#else
                                ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_0] = ((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->qp;
                                ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_0] = ((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->slice_type;
#endif
                                // Increment the Reference's liveCount by the number of tiles in the input picture
                                eb_object_inc_live_count(
                                    referenceEntryPtr->referenceObjectPtr,
                                    1);

                                // Decrement the Reference's dependentCount Count
                                --referenceEntryPtr->dependentCount;

                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependentCount != ~0u),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR1);

                            }
                        }

                        // Configure List1
                        if (entryPictureControlSetPtr->slice_type == B_SLICE) {

                            if (entryPictureControlSetPtr->ref_list1_count) {
                                referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                    ((int32_t)inputEntryPtr->referenceEntryIndex) - inputEntryPtr->list1Ptr->referenceList,
                                    REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                // Set the Reference Object
                                ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1] = referenceEntryPtr->referenceObjectPtr;

                                ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_1] = (uint8_t)((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->qp;
                                ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_1] = ((EbReferenceObject_t*)referenceEntryPtr->referenceObjectPtr->object_ptr)->slice_type;



                                // Increment the Reference's liveCount by the number of tiles in the input picture
                                eb_object_inc_live_count(
                                    referenceEntryPtr->referenceObjectPtr,
                                    1);

                                // Decrement the Reference's dependentCount Count
                                --referenceEntryPtr->dependentCount;

                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependentCount != ~0u),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR1);
                            }
                        }

                        // Adjust the Slice-type if the Lists are Empty, but don't reset the Prediction Structure
                        entryPictureControlSetPtr->slice_type =
                            (entryPictureControlSetPtr->ref_list1_count > 0) ? B_SLICE :
                            (entryPictureControlSetPtr->ref_list0_count > 0) ? P_SLICE :
                            I_SLICE;


                        // Increment the sequenceControlSet Wrapper's live count by 1 for only the pictures which are used as reference
                        if (ChildPictureControlSetPtr->parent_pcs_ptr->is_used_as_reference_flag) {
                            eb_object_inc_live_count(
                                ChildPictureControlSetPtr->parent_pcs_ptr->sequence_control_set_wrapper_ptr,
                                1);
                        }


                        // Get Empty Results Object
                        eb_get_empty_object(
                            context_ptr->pictureManagerOutputFifoPtr,
                            &outputWrapperPtr);

                        rateControlTasksPtr = (RateControlTasks_t*)outputWrapperPtr->object_ptr;
                        rateControlTasksPtr->pictureControlSetWrapperPtr = ChildPictureControlSetWrapperPtr;
                        rateControlTasksPtr->taskType = RC_PICTURE_MANAGER_RESULT;

                        // Post the Full Results Object
                        eb_post_full_object(outputWrapperPtr);

                        // Remove the Input Entry from the Input Queue
                        inputEntryPtr->inputObjectPtr = (EbObjectWrapper_t*)EB_NULL;

                    }
                }

                // Increment the HeadIndex if the head is null
                encode_context_ptr->input_picture_queue_head_index =
                    (encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]->inputObjectPtr) ? encode_context_ptr->input_picture_queue_head_index :
                    (encode_context_ptr->input_picture_queue_head_index == INPUT_QUEUE_MAX_DEPTH - 1) ? 0
                    : encode_context_ptr->input_picture_queue_head_index + 1;

                // Increment the inputQueueIndex Iterator
                inputQueueIndex = (inputQueueIndex == INPUT_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;

            }

            // Walk the reference queue and remove entries that have been completely referenced.
            referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;
            while (referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) {

                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                // Remove the entry & release the reference if there are no remaining references
                if ((referenceEntryPtr->dependentCount == 0) &&
                    (referenceEntryPtr->referenceAvailable) &&
                    (referenceEntryPtr->releaseEnable) &&
                    (referenceEntryPtr->referenceObjectPtr))
                {
                    // Release the nominal liveCount value
                    eb_release_object(referenceEntryPtr->referenceObjectPtr);

                    referenceEntryPtr->referenceObjectPtr = (EbObjectWrapper_t*)EB_NULL;
                    referenceEntryPtr->referenceAvailable = EB_FALSE;
                    referenceEntryPtr->is_used_as_reference_flag = EB_FALSE;
                }

                // Increment the HeadIndex if the head is empty
                encode_context_ptr->reference_picture_queue_head_index =
                    (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->releaseEnable == EB_FALSE) ? encode_context_ptr->reference_picture_queue_head_index :
                    (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->referenceAvailable == EB_FALSE &&
                        encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->is_used_as_reference_flag == EB_TRUE) ? encode_context_ptr->reference_picture_queue_head_index :
                        (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->dependentCount > 0) ? encode_context_ptr->reference_picture_queue_head_index :
                    (encode_context_ptr->reference_picture_queue_head_index == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0
                    : encode_context_ptr->reference_picture_queue_head_index + 1;
                // Increment the referenceQueueIndex Iterator
                referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;
            }
        }

        // Release the Input Picture Demux Results
        eb_release_object(inputPictureDemuxWrapperPtr);

    }
    return EB_NULL;
}
