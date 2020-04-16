/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecMemInit_h
#define EbDecMemInit_h

#ifdef __cplusplus
extern "C" {
#endif

extern EbMemoryMapEntry *svt_dec_memory_map;
extern uint32_t         *svt_dec_memory_map_index;
extern uint64_t         *svt_dec_total_lib_memory;
extern uint32_t          svt_dec_lib_malloc_count;

extern EbMemoryMapEntry                 *memory_map_start_address;
extern EbMemoryMapEntry                 *memory_map_end_address;

#ifdef _WIN32
#define EB_ALLIGN_MALLOC_DEC(type, pointer, n_elements, pointer_class)                  \
    pointer = (type)_aligned_malloc(n_elements, ALVALUE);                               \
    if (pointer == (type)NULL)                                                       \
        return EB_ErrorInsufficientResources;                                           \
    else {                                                                              \
        EbMemoryMapEntry *node = malloc(sizeof(EbMemoryMapEntry));                      \
        if (node == (EbMemoryMapEntry *)NULL) return EB_ErrorInsufficientResources;  \
        node->ptr_type     = pointer_class;                                             \
        node->ptr          = (EbPtr)pointer;                                            \
        node->prev_entry   = (EbPtr)svt_dec_memory_map;                                 \
        svt_dec_memory_map = node;                                                      \
        (*svt_dec_memory_map_index)++;                                                  \
        if (n_elements % 8 == 0)                                                        \
            *svt_dec_total_lib_memory += ((n_elements) + sizeof(EbMemoryMapEntry));     \
        else                                                                            \
            *svt_dec_total_lib_memory +=                                                \
                (((n_elements) + (8 - ((n_elements) % 8))) + sizeof(EbMemoryMapEntry)); \
        svt_dec_lib_malloc_count++;                                                     \
    }
#else
#define EB_ALLIGN_MALLOC_DEC(type, pointer, n_elements, pointer_class)                  \
    if (posix_memalign((void **)(&(pointer)), ALVALUE, n_elements) != 0)                \
        return EB_ErrorInsufficientResources;                                           \
    else {                                                                              \
        pointer                = (type)pointer;                                         \
        EbMemoryMapEntry *node = malloc(sizeof(EbMemoryMapEntry));                      \
        if (node == (EbMemoryMapEntry *)NULL) return EB_ErrorInsufficientResources;  \
        node->ptr_type     = pointer_class;                                             \
        node->ptr          = (EbPtr)pointer;                                            \
        node->prev_entry   = (EbPtr)svt_dec_memory_map;                                 \
        svt_dec_memory_map = node;                                                      \
        (*svt_dec_memory_map_index)++;                                                  \
        if (n_elements % 8 == 0)                                                        \
            *svt_dec_total_lib_memory += ((n_elements) + sizeof(EbMemoryMapEntry));     \
        else                                                                            \
            *svt_dec_total_lib_memory +=                                                \
                (((n_elements) + (8 - ((n_elements) % 8))) + sizeof(EbMemoryMapEntry)); \
        svt_dec_lib_malloc_count++;                                                     \
    }
#endif
#define EB_MALLOC_DEC(type, pointer, n_elements, pointer_class)                         \
    pointer = (type)malloc(n_elements);                                                 \
    if (pointer == (type)NULL)                                                       \
        return EB_ErrorInsufficientResources;                                           \
    else {                                                                              \
        EbMemoryMapEntry *node = malloc(sizeof(EbMemoryMapEntry));                      \
        if (node == (EbMemoryMapEntry *)NULL) return EB_ErrorInsufficientResources;  \
        node->ptr_type     = pointer_class;                                             \
        node->ptr          = (EbPtr)pointer;                                            \
        node->prev_entry   = (EbPtr)svt_dec_memory_map;                                 \
        svt_dec_memory_map = node;                                                      \
        (*svt_dec_memory_map_index)++;                                                  \
        if (n_elements % 8 == 0)                                                        \
            *svt_dec_total_lib_memory += ((n_elements) + sizeof(EbMemoryMapEntry));     \
        else                                                                            \
            *svt_dec_total_lib_memory +=                                                \
                (((n_elements) + (8 - ((n_elements) % 8))) + sizeof(EbMemoryMapEntry)); \
        svt_dec_lib_malloc_count++;                                                     \
    }

EbErrorType dec_eb_recon_picture_buffer_desc_ctor(EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr,
    EbBool is_16bit_pipeline);

EbErrorType dec_mem_init(EbDecHandle *dec_handle_ptr);

EbErrorType init_dec_mod_ctxt(EbDecHandle *dec_handle_ptr, void **dec_mod_ctxt);

#ifdef __cplusplus
}
#endif
#endif // EbDecMemInit_h
