/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbEncodeContext.h"
#include "EbPictureManagerQueue.h"
#include "EbCabacContextModel.h"
#include "EbSvtAv1ErrorCodes.h"

EbErrorType encode_context_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    uint32_t pictureIndex;
    EbErrorType return_error = EB_ErrorNone;

    EncodeContext_t *encode_context_ptr;
    EB_MALLOC(EncodeContext_t*, encode_context_ptr, sizeof(EncodeContext_t), EB_N_PTR);
    *object_dbl_ptr = (EbPtr)encode_context_ptr;

    object_init_data_ptr = 0;
    CHECK_REPORT_ERROR(
        (object_init_data_ptr == 0),
        encode_context_ptr->app_callback_ptr,
        EB_ENC_EC_ERROR29);

    // Callback Functions
    encode_context_ptr->app_callback_ptr = (EbCallback_t*)EB_NULL;

    EB_CREATEMUTEX(EbHandle, encode_context_ptr->total_number_of_recon_frame_mutex, sizeof(EbHandle), EB_MUTEX);
    encode_context_ptr->total_number_of_recon_frames = 0;
    encode_context_ptr->statistics_port_active = EB_FALSE;
    
    // Output Buffer Fifos
    encode_context_ptr->stream_output_fifo_ptr = (EbFifo_t*)EB_NULL;
    encode_context_ptr->recon_output_fifo_ptr = (EbFifo_t*)EB_NULL;

    // Picture Buffer Fifos
    encode_context_ptr->reference_picture_pool_fifo_ptr = (EbFifo_t*)EB_NULL;
    encode_context_ptr->pa_reference_picture_pool_fifo_ptr = (EbFifo_t*)EB_NULL;

    // Picture Decision Reordering Queue
    encode_context_ptr->picture_decision_reorder_queue_head_index = 0;
    EB_MALLOC(PictureDecisionReorderEntry_t**, encode_context_ptr->picture_decision_reorder_queue, sizeof(PictureDecisionReorderEntry_t*) * PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = picture_decision_reorder_entry_ctor(
            &(encode_context_ptr->picture_decision_reorder_queue[pictureIndex]),
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Manager Reordering Queue
    encode_context_ptr->picture_manager_reorder_queue_head_index = 0;
    EB_MALLOC(PictureManagerReorderEntry_t**, encode_context_ptr->picture_manager_reorder_queue, sizeof(PictureManagerReorderEntry_t*) * PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = picture_manager_reorder_entry_ctor(
            &(encode_context_ptr->picture_manager_reorder_queue[pictureIndex]),
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }


    // Picture Manager Pre-Assignment Buffer
    encode_context_ptr->pre_assignment_buffer_intra_count = 0;
    encode_context_ptr->pre_assignment_buffer_idr_count = 0;
    encode_context_ptr->pre_assignment_buffer_scene_change_count = 0;
    encode_context_ptr->pre_assignment_buffer_scene_change_index = 0;
    encode_context_ptr->pre_assignment_buffer_eos_flag = EB_FALSE;
    encode_context_ptr->decode_base_number = 0;

    encode_context_ptr->pre_assignment_buffer_count = 0;

    EB_MALLOC(EbObjectWrapper_t**, encode_context_ptr->pre_assignment_buffer, sizeof(EbObjectWrapper_t*) * PRE_ASSIGNMENT_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < PRE_ASSIGNMENT_MAX_DEPTH; ++pictureIndex) {
        encode_context_ptr->pre_assignment_buffer[pictureIndex] = (EbObjectWrapper_t*)EB_NULL;
    }

    // Picture Manager Input Queue
    encode_context_ptr->input_picture_queue_head_index = 0;
    encode_context_ptr->input_picture_queue_tail_index = 0;
    EB_MALLOC(InputQueueEntry_t**, encode_context_ptr->input_picture_queue, sizeof(InputQueueEntry_t*) * INPUT_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < INPUT_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = input_queue_entry_ctor(
            &(encode_context_ptr->input_picture_queue[pictureIndex]));
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Manager Reference Queue
    encode_context_ptr->reference_picture_queue_head_index = 0;
    encode_context_ptr->reference_picture_queue_tail_index = 0;
    EB_MALLOC(ReferenceQueueEntry_t**, encode_context_ptr->reference_picture_queue, sizeof(ReferenceQueueEntry_t*) * REFERENCE_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < REFERENCE_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = reference_queue_entry_ctor(
            &(encode_context_ptr->reference_picture_queue[pictureIndex]));
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Picture Decision PA Reference Queue
    encode_context_ptr->picture_decision_pa_reference_queue_head_index = 0;
    encode_context_ptr->picture_decision_pa_reference_queue_tail_index = 0;
    EB_MALLOC(PaReferenceQueueEntry_t**, encode_context_ptr->picture_decision_pa_reference_queue, sizeof(PaReferenceQueueEntry_t*) * PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = pa_reference_queue_entry_ctor(
            &(encode_context_ptr->picture_decision_pa_reference_queue[pictureIndex]));
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // Initial Rate Control Reordering Queue
    encode_context_ptr->initial_rate_control_reorder_queue_head_index = 0;
    EB_MALLOC(InitialRateControlReorderEntry_t**, encode_context_ptr->initial_rate_control_reorder_queue, sizeof(InitialRateControlReorderEntry_t*) * INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = InitialRateControlReorderEntryCtor(
            &(encode_context_ptr->initial_rate_control_reorder_queue[pictureIndex]),
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    // High level Rate Control histogram Queue
    encode_context_ptr->hl_rate_control_historgram_queue_head_index = 0;

    EB_MALLOC(HlRateControlHistogramEntry_t**, encode_context_ptr->hl_rate_control_historgram_queue, sizeof(HlRateControlHistogramEntry_t*) * HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = HlRateControlHistogramEntryCtor(
            &(encode_context_ptr->hl_rate_control_historgram_queue[pictureIndex]),
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    // HLRateControl Historgram Queue Mutex
    EB_CREATEMUTEX(EbHandle, encode_context_ptr->hl_rate_control_historgram_queue_mutex, sizeof(EbHandle), EB_MUTEX);

    // Packetization Reordering Queue
    encode_context_ptr->packetization_reorder_queue_head_index = 0;
    EB_MALLOC(PacketizationReorderEntry_t**, encode_context_ptr->packetization_reorder_queue, sizeof(PacketizationReorderEntry_t*) * PACKETIZATION_REORDER_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < PACKETIZATION_REORDER_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = packetization_reorder_entry_ctor(
            &(encode_context_ptr->packetization_reorder_queue[pictureIndex]),
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    encode_context_ptr->intra_period_position = 0;
    encode_context_ptr->pred_struct_position = 0;
    encode_context_ptr->current_input_poc = -1;
    encode_context_ptr->elapsed_non_idr_count = 0;
    encode_context_ptr->elapsed_non_cra_count = 0;
    encode_context_ptr->initial_picture = EB_TRUE;

    encode_context_ptr->last_idr_picture = 0;

    // Sequence Termination Flags
    encode_context_ptr->terminating_picture_number = ~0u;
    encode_context_ptr->terminating_sequence_flag_received = EB_FALSE;

    // Signalling the need for a td structure to be written in the bitstream - on when the sequence starts
    encode_context_ptr->td_needed = EB_TRUE;

    // Prediction Structure Group
    encode_context_ptr->prediction_structure_group_ptr = (PredictionStructureGroup_t*)EB_NULL;

    // MD Rate Estimation Array
    EB_MALLOC(MdRateEstimationContext_t*, encode_context_ptr->md_rate_estimation_array, sizeof(MdRateEstimationContext_t) * TOTAL_NUMBER_OF_MD_RATE_ESTIMATION_CASE_BUFFERS, EB_N_PTR);

    memset(encode_context_ptr->md_rate_estimation_array, 0, sizeof(MdRateEstimationContext_t) * TOTAL_NUMBER_OF_MD_RATE_ESTIMATION_CASE_BUFFERS);

    return_error = MdRateEstimationContextCtor(encode_context_ptr->md_rate_estimation_array);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Temporal Filter

    // Rate Control Bit Tables
    EB_MALLOC(RateControlTables_t*, encode_context_ptr->rate_control_tables_array, sizeof(RateControlTables_t) * TOTAL_NUMBER_OF_INITIAL_RC_TABLES_ENTRY, EB_N_PTR);

    return_error = rate_control_tables_ctor(encode_context_ptr->rate_control_tables_array);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // RC Rate Table Update Mutex
    EB_CREATEMUTEX(EbHandle, encode_context_ptr->rate_table_update_mutex, sizeof(EbHandle), EB_MUTEX);

    encode_context_ptr->rate_control_tables_array_updated = EB_FALSE;

    EB_CREATEMUTEX(EbHandle, encode_context_ptr->sc_buffer_mutex, sizeof(EbHandle), EB_MUTEX);
    encode_context_ptr->sc_buffer = 0;
    encode_context_ptr->sc_frame_in = 0;
    encode_context_ptr->sc_frame_out = 0;

    encode_context_ptr->enc_mode = SPEED_CONTROL_INIT_MOD;

    encode_context_ptr->previous_selected_ref_qp = 32;
    encode_context_ptr->max_coded_poc = 0;
    encode_context_ptr->max_coded_poc_selected_ref_qp = 32;

    encode_context_ptr->shared_reference_mutex = eb_create_mutex();
    if (encode_context_ptr->shared_reference_mutex == (EbHandle)EB_NULL) {
        return EB_ErrorInsufficientResources;
    }
    else {
        memory_map[*(memory_map_index)].ptrType = EB_MUTEX;
        memory_map[(*(memory_map_index))++].ptr = encode_context_ptr->shared_reference_mutex;
        *total_lib_memory += (sizeof(EbHandle));
    }


    return EB_ErrorNone;
}

