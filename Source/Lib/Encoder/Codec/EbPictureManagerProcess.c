// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureManagerProcess.h"
#include "EbReferenceObject.h"
#include "EbPictureDemuxResults.h"
#include "EbPictureManagerQueue.h"
#include "EbPredictionStructure.h"
#include "EbRateControlTasks.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbEntropyCoding.h"

void set_tile_info(PictureParentControlSet * pcs_ptr);
extern MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx);
/************************************************
 * Defines
 ************************************************/
#define POC_CIRCULAR_ADD(base, offset)             (((base) + (offset)))

 /************************************************
  * Configure Picture edges
  ************************************************/
static void ConfigurePictureEdges(
    SequenceControlSet *scs_ptr,
    PictureControlSet  *ppsPtr)
{
    // Tiles Initialisation
    const uint16_t picture_width_in_sb = (scs_ptr->seq_header.max_frame_width + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    const uint16_t picture_height_in_sb = (scs_ptr->seq_header.max_frame_height + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    unsigned x_lcu_index, y_lcu_index, sb_index;

    // LCU-loops
    for (y_lcu_index = 0; y_lcu_index < picture_height_in_sb; ++y_lcu_index) {
        for (x_lcu_index = 0; x_lcu_index < picture_width_in_sb; ++x_lcu_index) {
            sb_index = (uint16_t)(x_lcu_index + y_lcu_index * picture_width_in_sb);
            ppsPtr->sb_ptr_array[sb_index]->picture_left_edge_flag = (x_lcu_index == 0) ? EB_TRUE : EB_FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_top_edge_flag = (y_lcu_index == 0) ? EB_TRUE : EB_FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_right_edge_flag = (x_lcu_index == (unsigned)(picture_width_in_sb - 1)) ? EB_TRUE : EB_FALSE;
        }
    }

    return;
}
#if TWO_PASS
void write_stat_to_file(
    SequenceControlSet    *sequence_control_set_ptr,
    stat_struct_t          stat_struct,
    uint64_t               ref_poc);
#endif
/************************************************
 * Picture Manager Context Constructor
 ************************************************/
EbErrorType picture_manager_context_ctor(
    PictureManagerContext  *context_ptr,
    EbFifo                 *picture_input_fifo_ptr,
    EbFifo                 *picture_manager_output_fifo_ptr,
    EbFifo                **picture_control_set_fifo_ptr_array)
{
    context_ptr->picture_input_fifo_ptr = picture_input_fifo_ptr;
    context_ptr->picture_manager_output_fifo_ptr = picture_manager_output_fifo_ptr;
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
    PictureManagerContext         *context_ptr = (PictureManagerContext*)input_ptr;

    EbObjectWrapper               *ChildPictureControlSetWrapperPtr;
    PictureControlSet             *ChildPictureControlSetPtr;
    PictureParentControlSet       *picture_control_set_ptr;
    SequenceControlSet            *sequence_control_set_ptr;
    EncodeContext                 *encode_context_ptr;

    EbObjectWrapper               *inputPictureDemuxWrapperPtr;
    PictureDemuxResults           *inputPictureDemuxPtr;

    EbObjectWrapper               *outputWrapperPtr;
    RateControlTasks              *rateControlTasksPtr;

    EbBool                           availabilityFlag;

    PredictionStructureEntry      *predPositionPtr;
    InputQueueEntry               *inputEntryPtr;
    uint32_t                         inputQueueIndex;
    uint64_t                         current_input_poc;
    ReferenceQueueEntry           *referenceEntryPtr;
    uint32_t                         referenceQueueIndex;
    uint64_t                         ref_poc;
    uint32_t                         depIdx;
    uint64_t                         depPoc;
    uint32_t                         depListCount;
    PictureParentControlSet       *entryPictureControlSetPtr;
    SequenceControlSet            *entrySequenceControlSetPtr;

    // Initialization
    uint8_t                          picture_width_in_sb;
    uint8_t                          picture_height_in_sb;
    PictureManagerReorderEntry    *queueEntryPtr;
    int32_t                          queueEntryIndex;

    // Debug
    uint32_t loopCount = 0;

    for (;;) {
        // Get Input Full Object
        eb_get_full_object(
            context_ptr->picture_input_fifo_ptr,
            &inputPictureDemuxWrapperPtr);

        inputPictureDemuxPtr = (PictureDemuxResults*)inputPictureDemuxWrapperPtr->object_ptr;

        // *Note - This should be overhauled and/or replaced when we
        //   need hierarchical support.
        loopCount++;

        switch (inputPictureDemuxPtr->picture_type) {
        case EB_PIC_INPUT:

            picture_control_set_ptr = (PictureParentControlSet*)inputPictureDemuxPtr->picture_control_set_wrapper_ptr->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

            //printf("\nPicture Manager Process @ %d \n ", picture_control_set_ptr->picture_number);

            queueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number_alt - encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index]->picture_number);
            queueEntryIndex += encode_context_ptr->picture_manager_reorder_queue_head_index;
            queueEntryIndex = (queueEntryIndex > PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH : queueEntryIndex;
            queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[queueEntryIndex];
            if (queueEntryPtr->parent_pcs_wrapper_ptr != NULL) {
                CHECK_REPORT_ERROR_NC(
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PD_ERROR8);
            }
            else {
                queueEntryPtr->parent_pcs_wrapper_ptr = inputPictureDemuxPtr->picture_control_set_wrapper_ptr;
                queueEntryPtr->picture_number = picture_control_set_ptr->picture_number_alt;
            }
            // Process the head of the Picture Manager Reorder Queue
            queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index];

            while (queueEntryPtr->parent_pcs_wrapper_ptr != EB_NULL) {
                picture_control_set_ptr = (PictureParentControlSet*)queueEntryPtr->parent_pcs_wrapper_ptr->object_ptr;

                predPositionPtr = picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];
                // If there was a change in the number of temporal layers, then cleanup the Reference Queue's Dependent Counts
                if (picture_control_set_ptr->hierarchical_layers_diff != 0) {
                    // Dynamic GOP
                    PredictionStructure          *next_pred_struct_ptr;
                    PredictionStructureEntry     *next_base_layer_pred_position_ptr;

                    uint32_t                        dependant_list_positive_entries;
                    uint32_t                        dependant_list_removed_entries;

                    referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;

                    while (referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) {
                        referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                        if (referenceEntryPtr->picture_number == (picture_control_set_ptr->picture_number - 1)) { // Picture where the change happened

                            // Get the prediction struct entry of the next GOP structure
                            next_pred_struct_ptr = get_prediction_structure(
                                encode_context_ptr->prediction_structure_group_ptr,
                                picture_control_set_ptr->pred_structure,
                                sequence_control_set_ptr->reference_count,
                                picture_control_set_ptr->hierarchical_levels);

                            // Get the prediction struct of a picture in temporal layer 0 (from the new GOP structure)
                            next_base_layer_pred_position_ptr = next_pred_struct_ptr->pred_struct_entry_ptr_array[next_pred_struct_ptr->pred_struct_entry_count - 1];

                            // Remove all positive entries from the dependant lists
                            dependant_list_positive_entries = 0;
                            for (depIdx = 0; depIdx < referenceEntryPtr->list0.list_count; ++depIdx) {
                                if (referenceEntryPtr->list0.list[depIdx] >= 0)
                                    dependant_list_positive_entries++;
                            }
                            referenceEntryPtr->list0.list_count = referenceEntryPtr->list0.list_count - dependant_list_positive_entries;

                            dependant_list_positive_entries = 0;
                            for (depIdx = 0; depIdx < referenceEntryPtr->list1.list_count; ++depIdx) {
                                if (referenceEntryPtr->list1.list[depIdx] >= 0)
                                    dependant_list_positive_entries++;
                            }
                            referenceEntryPtr->list1.list_count = referenceEntryPtr->list1.list_count - dependant_list_positive_entries;

                            for (depIdx = 0; depIdx < next_base_layer_pred_position_ptr->dep_list0.list_count; ++depIdx) {
                                if (next_base_layer_pred_position_ptr->dep_list0.list[depIdx] >= 0)
                                    referenceEntryPtr->list0.list[referenceEntryPtr->list0.list_count++] = next_base_layer_pred_position_ptr->dep_list0.list[depIdx];
                            }

                            for (depIdx = 0; depIdx < next_base_layer_pred_position_ptr->dep_list1.list_count; ++depIdx) {
                                if (next_base_layer_pred_position_ptr->dep_list1.list[depIdx] >= 0)
                                    referenceEntryPtr->list1.list[referenceEntryPtr->list1.list_count++] = next_base_layer_pred_position_ptr->dep_list1.list[depIdx];
                            }

                            // Update the dependant count update
                            dependant_list_removed_entries = referenceEntryPtr->dep_list0_count + referenceEntryPtr->dep_list1_count - referenceEntryPtr->dependent_count;
                            referenceEntryPtr->dep_list0_count = (referenceEntryPtr->is_alt_ref) ? referenceEntryPtr->list0.list_count + 1 : referenceEntryPtr->list0.list_count;
                            referenceEntryPtr->dep_list1_count = referenceEntryPtr->list1.list_count;
                            referenceEntryPtr->dependent_count = referenceEntryPtr->dep_list0_count + referenceEntryPtr->dep_list1_count - dependant_list_removed_entries;
                        }
                        else {
                            // Modify Dependent List0
                            depListCount = referenceEntryPtr->list0.list_count;
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
                                    --referenceEntryPtr->dependent_count;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR3);
                                }
                            }

                            // Modify Dependent List1
                            depListCount = referenceEntryPtr->list1.list_count;
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
                                    --referenceEntryPtr->dependent_count;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR3);
                                }
                            }
                        }

                        // Increment the referenceQueueIndex Iterator
                        referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;
                    }
                }
                // If there was an I-frame or Scene Change, then cleanup the Reference Queue's Dependent Counts
                if (picture_control_set_ptr->slice_type == I_SLICE)
                {
                    referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;
                    while (referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) {
                        referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                        // Modify Dependent List0
                        depListCount = referenceEntryPtr->list0.list_count;
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
                                --referenceEntryPtr->dependent_count;
                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependent_count != ~0u),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR3);
                            }
                        }

                        // Modify Dependent List1
                        depListCount = referenceEntryPtr->list1.list_count;
                        for (depIdx = 0; depIdx < depListCount; ++depIdx) {
                            current_input_poc = picture_control_set_ptr->picture_number;

                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                            // current_input_poc += (current_input_poc < referenceEntryPtr->picture_number) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                            depPoc = POC_CIRCULAR_ADD(
                                referenceEntryPtr->picture_number,
                                referenceEntryPtr->list1.list[depIdx]/*,
                                sequence_control_set_ptr->bits_for_picture_order_count*/);

                                // If Dependent POC is greater or equal to the IDR POC or if we inserted trailing Ps
                            if (((depPoc >= current_input_poc) || (((picture_control_set_ptr->pre_assignment_buffer_count != picture_control_set_ptr->pred_struct_ptr->pred_struct_period) || (picture_control_set_ptr->idr_flag == EB_TRUE)) && (depPoc > (current_input_poc - picture_control_set_ptr->pre_assignment_buffer_count)))) && referenceEntryPtr->list1.list[depIdx]) {
                                referenceEntryPtr->list1.list[depIdx] = 0;

                                // Decrement the Reference's referenceCount
                                --referenceEntryPtr->dependent_count;
                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr->dependent_count != ~0u),
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
                    referenceEntryPtr = (ReferenceQueueEntry*)EB_NULL;
                }

                // Check if the EnhancedPictureQueue is full.
                // *Note - Having the number of Enhanced Pictures less than the InputQueueSize should ensure this never gets hit

                CHECK_REPORT_ERROR(
                    (((encode_context_ptr->input_picture_queue_head_index != encode_context_ptr->input_picture_queue_tail_index) || (encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]->input_object_ptr == EB_NULL))),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR4);

                // Place Picture in input queue
                inputEntryPtr = encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_tail_index];
                inputEntryPtr->input_object_ptr = queueEntryPtr->parent_pcs_wrapper_ptr;
                // Since the overlay picture is not added to the reference queue, reference_entry_index points to the previous picture which is the alt ref
                inputEntryPtr->reference_entry_index = (!picture_control_set_ptr->is_overlay) ? encode_context_ptr->reference_picture_queue_tail_index :
                    (encode_context_ptr->reference_picture_queue_tail_index == 0) ? REFERENCE_QUEUE_MAX_DEPTH - 1 : encode_context_ptr->reference_picture_queue_tail_index - 1;
                encode_context_ptr->input_picture_queue_tail_index =
                    (encode_context_ptr->input_picture_queue_tail_index == INPUT_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->input_picture_queue_tail_index + 1;

                // Copy the reference lists into the inputEntry and
                // set the Reference Counts Based on Temporal Layer and how many frames are active
                inputEntryPtr->list0_ptr = &predPositionPtr->ref_list0;
                inputEntryPtr->list1_ptr = &predPositionPtr->ref_list1;
                if (!picture_control_set_ptr->is_overlay) {
                    // Check if the ReferencePictureQueue is full.
                    CHECK_REPORT_ERROR(
                        (((encode_context_ptr->reference_picture_queue_head_index != encode_context_ptr->reference_picture_queue_tail_index) || (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->reference_object_ptr == EB_NULL))),
                        encode_context_ptr->app_callback_ptr,
                        EB_ENC_PM_ERROR5);

                    // Create Reference Queue Entry even if picture will not be referenced
                    referenceEntryPtr = encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_tail_index];
                    referenceEntryPtr->picture_number = picture_control_set_ptr->picture_number;
                    referenceEntryPtr->reference_object_ptr = (EbObjectWrapper*)EB_NULL;
                    referenceEntryPtr->release_enable = EB_TRUE;
                    referenceEntryPtr->reference_available = EB_FALSE;
                    referenceEntryPtr->slice_type = picture_control_set_ptr->slice_type;
                    referenceEntryPtr->temporal_layer_index = picture_control_set_ptr->temporal_layer_index;
                    referenceEntryPtr->frame_context_updated = EB_FALSE;
                    referenceEntryPtr->is_alt_ref = picture_control_set_ptr->is_alt_ref;
                    referenceEntryPtr->feedback_arrived = EB_FALSE;
                    referenceEntryPtr->is_used_as_reference_flag = picture_control_set_ptr->is_used_as_reference_flag;
                    encode_context_ptr->reference_picture_queue_tail_index =
                        (encode_context_ptr->reference_picture_queue_tail_index == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->reference_picture_queue_tail_index + 1;

                    // Copy the Dependent Lists
                    // *Note - we are removing any leading picture dependencies for now
                    referenceEntryPtr->list0.list_count = 0;
                    for (depIdx = 0; depIdx < predPositionPtr->dep_list0.list_count; ++depIdx) {
                        if (predPositionPtr->dep_list0.list[depIdx] >= 0)
                            referenceEntryPtr->list0.list[referenceEntryPtr->list0.list_count++] = predPositionPtr->dep_list0.list[depIdx];
                    }

                    referenceEntryPtr->list1.list_count = predPositionPtr->dep_list1.list_count;
                    for (depIdx = 0; depIdx < predPositionPtr->dep_list1.list_count; ++depIdx)
                        referenceEntryPtr->list1.list[depIdx] = predPositionPtr->dep_list1.list[depIdx];
                    referenceEntryPtr->dep_list0_count = (picture_control_set_ptr->is_alt_ref) ? referenceEntryPtr->list0.list_count + 1 : referenceEntryPtr->list0.list_count;

                    referenceEntryPtr->dep_list1_count = referenceEntryPtr->list1.list_count;
                    referenceEntryPtr->dependent_count = referenceEntryPtr->dep_list0_count + referenceEntryPtr->dep_list1_count;

                    CHECK_REPORT_ERROR(
                        (picture_control_set_ptr->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH < MAX_ELAPSED_IDR_COUNT),
                        encode_context_ptr->app_callback_ptr,
                        EB_ENC_PM_ERROR6);
                }
                // Release the Reference Buffer once we know it is not a reference
                if (picture_control_set_ptr->is_used_as_reference_flag == EB_FALSE) {
                    // Release the nominal live_count value
                    eb_release_object(picture_control_set_ptr->reference_picture_wrapper_ptr);
                    picture_control_set_ptr->reference_picture_wrapper_ptr = (EbObjectWrapper*)EB_NULL;
                }

                // Release the Picture Manager Reorder Queue
                queueEntryPtr->parent_pcs_wrapper_ptr = (EbObjectWrapper*)EB_NULL;
                queueEntryPtr->picture_number += PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH;

                // Increment the Picture Manager Reorder Queue
                encode_context_ptr->picture_manager_reorder_queue_head_index = (encode_context_ptr->picture_manager_reorder_queue_head_index == PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_manager_reorder_queue_head_index + 1;

                // Get the next entry from the Picture Manager Reorder Queue (Entry N+1)
                queueEntryPtr = encode_context_ptr->picture_manager_reorder_queue[encode_context_ptr->picture_manager_reorder_queue_head_index];
            }

            break;

        case EB_PIC_REFERENCE:

            sequence_control_set_ptr = (SequenceControlSet*)inputPictureDemuxPtr->sequence_control_set_wrapper_ptr->object_ptr;
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
                    referenceEntryPtr->reference_object_ptr = inputPictureDemuxPtr->reference_picture_wrapper_ptr;

                    // Set the reference availability
                    referenceEntryPtr->reference_available = EB_TRUE;
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
        case EB_PIC_FEEDBACK:
            sequence_control_set_ptr = (SequenceControlSet*)inputPictureDemuxPtr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;
            referenceQueueIndex = encode_context_ptr->reference_picture_queue_head_index;
            // Find the Reference in the Reference Queue
            do {
                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];
                if (referenceEntryPtr->picture_number == inputPictureDemuxPtr->picture_number) {
                    // Set the feedback arrived
                    referenceEntryPtr->feedback_arrived = EB_TRUE;
                    referenceEntryPtr->frame_context_updated = EB_TRUE;
                }
                // Increment the referenceQueueIndex Iterator
                referenceQueueIndex = (referenceQueueIndex == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : referenceQueueIndex + 1;

            } while ((referenceQueueIndex != encode_context_ptr->reference_picture_queue_tail_index) && (referenceEntryPtr->picture_number != inputPictureDemuxPtr->picture_number));

            //keep the release of SCS here because we still need the encodeContext structure here
            // Release the Reference's SequenceControlSet
            eb_release_object(inputPictureDemuxPtr->sequence_control_set_wrapper_ptr);

            break;
        default:

            sequence_control_set_ptr = (SequenceControlSet*)inputPictureDemuxPtr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

            CHECK_REPORT_ERROR_NC(
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR9);

            picture_control_set_ptr = (PictureParentControlSet*)EB_NULL;
            encode_context_ptr = (EncodeContext*)EB_NULL;

            break;
        }

        // ***********************************
        //  Common Code
        // *************************************

        // Walk the input queue and start all ready pictures.  Mark entry as null after started.  Increment the head as you go.
        if (encode_context_ptr != (EncodeContext*)EB_NULL) {
            inputQueueIndex = encode_context_ptr->input_picture_queue_head_index;
            while (inputQueueIndex != encode_context_ptr->input_picture_queue_tail_index) {
                inputEntryPtr = encode_context_ptr->input_picture_queue[inputQueueIndex];

                if (inputEntryPtr->input_object_ptr != EB_NULL) {
                    entryPictureControlSetPtr = (PictureParentControlSet*)inputEntryPtr->input_object_ptr->object_ptr;
                    entrySequenceControlSetPtr = (SequenceControlSet*)entryPictureControlSetPtr->sequence_control_set_wrapper_ptr->object_ptr;

                    availabilityFlag = EB_TRUE;

                    // Check RefList0 Availability
                    uint8_t refIdx;
                    for (refIdx = 0; refIdx < entryPictureControlSetPtr->ref_list0_count; ++refIdx) {
                        //if (entryPictureControlSetPtr->ref_list0_count)  // NM: to double check.
                        {
                            // hardcode the reference for the overlay frame
                            if (entryPictureControlSetPtr->is_overlay) {
                                referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                    ((int32_t)inputEntryPtr->reference_entry_index),
                                    REFERENCE_QUEUE_MAX_DEPTH);

                                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR10);

                                ref_poc = POC_CIRCULAR_ADD(
                                    entryPictureControlSetPtr->picture_number,
                                    0);
                            }
                            else {
                                referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                    ((int32_t)inputEntryPtr->reference_entry_index) -     // Base
                                    inputEntryPtr->list0_ptr->reference_list[refIdx],     // Offset
                                    REFERENCE_QUEUE_MAX_DEPTH);                         // Max

                                referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                CHECK_REPORT_ERROR(
                                    (referenceEntryPtr),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PM_ERROR10);

                                ref_poc = POC_CIRCULAR_ADD(
                                    entryPictureControlSetPtr->picture_number,
                                    -inputEntryPtr->list0_ptr->reference_list[refIdx]/*,
                                    entrySequenceControlSetPtr->bits_for_picture_order_count*/);
                            }
                                // Increment the current_input_poc is the case of POC rollover
                            current_input_poc = encode_context_ptr->current_input_poc;
                            //current_input_poc += ((current_input_poc < ref_poc) && (inputEntryPtr->list0_ptr->reference_list[refIdx] > 0)) ?
                            //    (1 << entrySequenceControlSetPtr->bits_for_picture_order_count) :
                            //    0;

                            availabilityFlag =
                                (availabilityFlag == EB_FALSE) ? EB_FALSE :   // Don't update if already False
                                (ref_poc > current_input_poc) ? EB_FALSE :   // The Reference has not been received as an Input Picture yet, then its availability is false
                                (!encode_context_ptr->terminating_sequence_flag_received &&
                                (sequence_control_set_ptr->static_config.rate_control_mode && entryPictureControlSetPtr->slice_type != I_SLICE
                                    && entryPictureControlSetPtr->temporal_layer_index == 0 && !referenceEntryPtr->feedback_arrived)) ? EB_FALSE :
                                    (entryPictureControlSetPtr->frame_end_cdf_update_mode && !referenceEntryPtr->frame_context_updated) ? EB_FALSE :
                                (referenceEntryPtr->reference_available) ? EB_TRUE :   // The Reference has been completed
                                EB_FALSE;     // The Reference has not been completed
                        }
                    }
                    // Check RefList1 Availability
                    if (entryPictureControlSetPtr->slice_type == B_SLICE) {
                        uint8_t refIdx;
                        for (refIdx = 0; refIdx < entryPictureControlSetPtr->ref_list1_count; ++refIdx) {
                            // if (entryPictureControlSetPtr->ref_list1_count) // NM: To double check
                            {
                                // If Reference is valid (non-zero), update the availability
                                if (inputEntryPtr->list1_ptr->reference_list[refIdx] != (int32_t)INVALID_POC) {
                                    referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                        ((int32_t)inputEntryPtr->reference_entry_index) -     // Base
                                        inputEntryPtr->list1_ptr->reference_list[refIdx],     // Offset
                                        REFERENCE_QUEUE_MAX_DEPTH);                         // Max

                                    referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PM_ERROR10);

                                    ref_poc = POC_CIRCULAR_ADD(
                                        entryPictureControlSetPtr->picture_number,
                                        -inputEntryPtr->list1_ptr->reference_list[refIdx]/*,
                                        entrySequenceControlSetPtr->bits_for_picture_order_count*/);

                                        // Increment the current_input_poc is the case of POC rollover
                                    current_input_poc = encode_context_ptr->current_input_poc;
                                    //current_input_poc += ((current_input_poc < ref_poc && inputEntryPtr->list1_ptr->reference_list[refIdx] > 0)) ?
                                    //    (1 << entrySequenceControlSetPtr->bits_for_picture_order_count) :
                                    //    0;

                                    availabilityFlag =
                                        (availabilityFlag == EB_FALSE) ? EB_FALSE :   // Don't update if already False
                                        (ref_poc > current_input_poc) ? EB_FALSE :   // The Reference has not been received as an Input Picture yet, then its availability is false
                                        (!encode_context_ptr->terminating_sequence_flag_received &&
                                        (sequence_control_set_ptr->static_config.rate_control_mode && entryPictureControlSetPtr->slice_type != I_SLICE
                                            && entryPictureControlSetPtr->temporal_layer_index == 0 && !referenceEntryPtr->feedback_arrived)) ? EB_FALSE :
                                            (entryPictureControlSetPtr->frame_end_cdf_update_mode && !referenceEntryPtr->frame_context_updated) ? EB_FALSE :
                                        (referenceEntryPtr->reference_available) ? EB_TRUE :   // The Reference has been completed
                                        EB_FALSE;     // The Reference has not been completed
                                }
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

                        ChildPictureControlSetPtr = (PictureControlSet*)ChildPictureControlSetWrapperPtr->object_ptr;

                        //1.Link The Child PCS to its Parent
                        ChildPictureControlSetPtr->picture_parent_control_set_wrapper_ptr = inputEntryPtr->input_object_ptr;
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
                        picture_width_in_sb = (uint8_t)((entrySequenceControlSetPtr->seq_header.max_frame_width + entrySequenceControlSetPtr->sb_size_pix - 1) / entrySequenceControlSetPtr->sb_size_pix);
                        picture_height_in_sb = (uint8_t)((entrySequenceControlSetPtr->seq_header.max_frame_height + entrySequenceControlSetPtr->sb_size_pix - 1) / entrySequenceControlSetPtr->sb_size_pix);

                        // EncDec Segments
                        enc_dec_segments_init(
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

                            for (row_index = 0; row_index < MAX_LCU_ROWS; ++row_index)
                                ChildPictureControlSetPtr->entropy_coding_row_array[row_index] = EB_FALSE;
                        }

                        ChildPictureControlSetPtr->parent_pcs_ptr->av1_cm->pcs_ptr = ChildPictureControlSetPtr;

                        struct PictureParentControlSet     *ppcs_ptr = ChildPictureControlSetPtr->parent_pcs_ptr;
                        Av1Common *const cm = ppcs_ptr->av1_cm;
                        int tile_row, tile_col;
                        uint32_t  x_lcu_index,  y_lcu_index;
                        const int tile_cols = ppcs_ptr->av1_cm->tiles_info.tile_cols;
                        const int tile_rows = ppcs_ptr->av1_cm->tiles_info.tile_rows;
                        TileInfo tile_info;
                        int sb_size_log2 = sequence_control_set_ptr->seq_header.sb_size_log2;
                        //Tile Loop
                        for (tile_row = 0; tile_row < tile_rows; tile_row++)
                        {
                            eb_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                            for (tile_col = 0; tile_col < tile_cols; tile_col++)
                            {
                                eb_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);

                                for ((y_lcu_index = cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                                     (y_lcu_index < (uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >> sb_size_log2);
                                     y_lcu_index++)
                                {
                                    for ((x_lcu_index = cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                                         (x_lcu_index < (uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >> sb_size_log2);
                                         x_lcu_index++)
                                    {
                                        int sb_index = (uint16_t)(x_lcu_index + y_lcu_index * picture_width_in_sb);
                                        ChildPictureControlSetPtr->sb_ptr_array[sb_index]->tile_info = tile_info;
                                    }
                                }
                            }
                        }
                        cm->mi_stride = ChildPictureControlSetPtr->mi_stride;
                        // Picture edges
                        ConfigurePictureEdges(entrySequenceControlSetPtr, ChildPictureControlSetPtr);

                        // Reset the qp array for DLF
                        EB_MEMSET(ChildPictureControlSetPtr->qp_array, 0, sizeof(uint8_t)*ChildPictureControlSetPtr->qp_array_size);

                        // Error resilience related
                        ChildPictureControlSetPtr->colocated_pu_ref_list = REF_LIST_0;     // to be modified

                        ChildPictureControlSetPtr->is_low_delay = (EbBool)(
                            ChildPictureControlSetPtr->parent_pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[ChildPictureControlSetPtr->parent_pcs_ptr->pred_struct_index]->positive_ref_pics_total_count == 0);

                        // Rate Control
                        ChildPictureControlSetPtr->dif_cu_delta_qp_depth = (uint8_t)entrySequenceControlSetPtr->input_resolution == INPUT_SIZE_4K_RANGE ? 3 : 2;

                        // Reset the Reference Lists
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));

                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint8_t));
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint8_t));

                        EB_MEMSET(ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EB_SLICE));
                        EB_MEMSET(ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EB_SLICE));
#if TWO_PASS
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_referenced_area_avg_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                        EB_MEMSET(ChildPictureControlSetPtr->ref_pic_referenced_area_avg_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
#endif
                        int8_t max_temporal_index = -1, ref_index = 0;
                        // Configure List0
                        if ((entryPictureControlSetPtr->slice_type == P_SLICE) || (entryPictureControlSetPtr->slice_type == B_SLICE)) {
                            uint8_t refIdx;
                            for (refIdx = 0; refIdx < entryPictureControlSetPtr->ref_list0_count; ++refIdx) {
                                if (entryPictureControlSetPtr->ref_list0_count) {
                                    // hardcode the reference for the overlay frame
                                    if (entryPictureControlSetPtr->is_overlay)
                                        referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                            ((int32_t)inputEntryPtr->reference_entry_index),
                                            REFERENCE_QUEUE_MAX_DEPTH);
                                    else
                                    referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                        ((int32_t)inputEntryPtr->reference_entry_index) - inputEntryPtr->list0_ptr->reference_list[refIdx],
                                        REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                    referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];
                                    if (entryPictureControlSetPtr->frame_end_cdf_update_mode) {
                                        ChildPictureControlSetPtr->ref_frame_context[svt_get_ref_frame_type(REF_LIST_0, refIdx) - LAST_FRAME] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->frame_context;
                                        if (max_temporal_index < (int8_t)referenceEntryPtr->temporal_layer_index) {
                                            max_temporal_index = (int8_t)referenceEntryPtr->temporal_layer_index;
                                            ref_index = svt_get_ref_frame_type(REF_LIST_0, refIdx) - LAST_FRAME;
                                            for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame)
                                                ChildPictureControlSetPtr->ref_global_motion[frame] =
                                                ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->global_motion[frame];
                                        }
                                    }
                                    // Set the Reference Object
                                    ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0][refIdx] = referenceEntryPtr->reference_object_ptr;
                                    ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_0][refIdx] = (uint8_t)((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->qp;
                                    ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_0][refIdx] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->slice_type;
#if TWO_PASS
                                    ChildPictureControlSetPtr->ref_pic_referenced_area_avg_array[REF_LIST_0][refIdx] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->referenced_area_avg;
#endif
                                    // Increment the Reference's liveCount by the number of tiles in the input picture
                                    eb_object_inc_live_count(
                                        referenceEntryPtr->reference_object_ptr,
                                        1);

                                    // Decrement the Reference's dependent_count Count
                                    --referenceEntryPtr->dependent_count;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PM_ERROR1);
                                }
                            }
                            //fill the non used spots to be used in MFMV.
                            for (refIdx = entryPictureControlSetPtr->ref_list0_count; refIdx < 4; ++refIdx)
                                ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0][refIdx] = ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0][0];

                            if (entryPictureControlSetPtr->ref_list1_count == 0) {
                                for (refIdx = entryPictureControlSetPtr->ref_list1_count; refIdx < 3; ++refIdx)
                                    ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1][refIdx] = ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_0][0];
                            }
                        }

                        // Configure List1
                        if (entryPictureControlSetPtr->slice_type == B_SLICE) {
                            uint8_t refIdx;
                            for (refIdx = 0; refIdx < entryPictureControlSetPtr->ref_list1_count; ++refIdx) {
                                if (entryPictureControlSetPtr->ref_list1_count) {
                                    referenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                        ((int32_t)inputEntryPtr->reference_entry_index) - inputEntryPtr->list1_ptr->reference_list[refIdx],
                                        REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                    referenceEntryPtr = encode_context_ptr->reference_picture_queue[referenceQueueIndex];
                                    if (entryPictureControlSetPtr->frame_end_cdf_update_mode) {
                                        ChildPictureControlSetPtr->ref_frame_context[svt_get_ref_frame_type(REF_LIST_1, refIdx) - LAST_FRAME] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->frame_context;
                                        if (max_temporal_index < (int8_t)referenceEntryPtr->temporal_layer_index && referenceEntryPtr->slice_type != I_SLICE/* && ChildPictureControlSetPtr->temporal_layer_index != 0*/) {
                                            max_temporal_index = (int8_t)referenceEntryPtr->temporal_layer_index;
                                            ref_index = svt_get_ref_frame_type(REF_LIST_1, refIdx) - LAST_FRAME;
                                            for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame)
                                                ChildPictureControlSetPtr->ref_global_motion[frame] =
                                                ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->global_motion[frame];
                                        }
                                    }
                                    // Set the Reference Object
                                    ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1][refIdx] = referenceEntryPtr->reference_object_ptr;

                                    ChildPictureControlSetPtr->ref_pic_qp_array[REF_LIST_1][refIdx] = (uint8_t)((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->qp;
                                    ChildPictureControlSetPtr->ref_slice_type_array[REF_LIST_1][refIdx] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->slice_type;
#if TWO_PASS
                                    ChildPictureControlSetPtr->ref_pic_referenced_area_avg_array[REF_LIST_1][refIdx] = ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->referenced_area_avg;
#endif
                                    // Increment the Reference's liveCount by the number of tiles in the input picture
                                    eb_object_inc_live_count(
                                        referenceEntryPtr->reference_object_ptr,
                                        1);

                                    // Decrement the Reference's dependent_count Count
                                    --referenceEntryPtr->dependent_count;

                                    CHECK_REPORT_ERROR(
                                        (referenceEntryPtr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PM_ERROR1);
                                }
                            }
                            //fill the non used spots to be used in MFMV.
                            if (entryPictureControlSetPtr->ref_list1_count) {
                                for (refIdx = entryPictureControlSetPtr->ref_list1_count; refIdx < 3; ++refIdx)
                                    ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1][refIdx] = ChildPictureControlSetPtr->ref_pic_ptr_array[REF_LIST_1][0];
                            }
                        }

                        // Adjust the Slice-type if the Lists are Empty, but don't reset the Prediction Structure
                        entryPictureControlSetPtr->slice_type =
                            (entryPictureControlSetPtr->ref_list1_count > 0) ? B_SLICE :
                            (entryPictureControlSetPtr->ref_list0_count > 0) ? P_SLICE :
                            I_SLICE;
                        if (entryPictureControlSetPtr->frame_end_cdf_update_mode) {
                            if (entryPictureControlSetPtr->slice_type != I_SLICE)
                                ChildPictureControlSetPtr->parent_pcs_ptr->frm_hdr.primary_ref_frame = ref_index;
                            else
                                ChildPictureControlSetPtr->parent_pcs_ptr->frm_hdr.primary_ref_frame = PRIMARY_REF_NONE;
                            ChildPictureControlSetPtr->parent_pcs_ptr->refresh_frame_context = REFRESH_FRAME_CONTEXT_BACKWARD;

                        }
                        else {
                            ChildPictureControlSetPtr->parent_pcs_ptr->frm_hdr.primary_ref_frame = PRIMARY_REF_NONE;
                            ChildPictureControlSetPtr->parent_pcs_ptr->refresh_frame_context = REFRESH_FRAME_CONTEXT_DISABLED;
                        }
                        // Increment the sequenceControlSet Wrapper's live count by 1 for only the pictures which are used as reference
                        if (ChildPictureControlSetPtr->parent_pcs_ptr->is_used_as_reference_flag) {
                            eb_object_inc_live_count(
                                ChildPictureControlSetPtr->parent_pcs_ptr->sequence_control_set_wrapper_ptr,
                                1);
                        }

                        // Get Empty Results Object
                        eb_get_empty_object(
                            context_ptr->picture_manager_output_fifo_ptr,
                            &outputWrapperPtr);

                        rateControlTasksPtr = (RateControlTasks*)outputWrapperPtr->object_ptr;
                        rateControlTasksPtr->picture_control_set_wrapper_ptr = ChildPictureControlSetWrapperPtr;
                        rateControlTasksPtr->task_type = RC_PICTURE_MANAGER_RESULT;

                        // Post the Full Results Object
                        eb_post_full_object(outputWrapperPtr);

                        // Remove the Input Entry from the Input Queue
                        inputEntryPtr->input_object_ptr = (EbObjectWrapper*)EB_NULL;
                    }
                }

                // Increment the head_index if the head is null
                encode_context_ptr->input_picture_queue_head_index =
                    (encode_context_ptr->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]->input_object_ptr) ? encode_context_ptr->input_picture_queue_head_index :
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
                if ((referenceEntryPtr->dependent_count == 0) &&
                    (referenceEntryPtr->reference_available) &&
                    (referenceEntryPtr->release_enable) &&
                    (referenceEntryPtr->reference_object_ptr))
                {
                    // Release the nominal live_count value
#if TWO_PASS
                    if (sequence_control_set_ptr->use_output_stat_file &&
                        referenceEntryPtr->reference_object_ptr->live_count == 1)
                        write_stat_to_file(
                            sequence_control_set_ptr,
                            ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->stat_struct,
                            ((EbReferenceObject*)referenceEntryPtr->reference_object_ptr->object_ptr)->ref_poc);
#endif
                    eb_release_object(referenceEntryPtr->reference_object_ptr);
                    referenceEntryPtr->reference_object_ptr = (EbObjectWrapper*)EB_NULL;
                    referenceEntryPtr->reference_available = EB_FALSE;
                    referenceEntryPtr->is_used_as_reference_flag = EB_FALSE;
                }

                // Increment the head_index if the head is empty
                encode_context_ptr->reference_picture_queue_head_index =
                    (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->release_enable == EB_FALSE) ? encode_context_ptr->reference_picture_queue_head_index :
                    (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->reference_available == EB_FALSE &&
                        encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->is_used_as_reference_flag == EB_TRUE) ? encode_context_ptr->reference_picture_queue_head_index :
                        (encode_context_ptr->reference_picture_queue[encode_context_ptr->reference_picture_queue_head_index]->dependent_count > 0) ? encode_context_ptr->reference_picture_queue_head_index :
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
// clang-format on
