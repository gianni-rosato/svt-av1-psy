/*
* Copyright (c) 2021, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>

#include "EbSvtAv1Metadata.h"

EB_API SvtMetadataT *svt_metadata_alloc(
    const uint32_t type, const uint8_t *data, const size_t sz) {
    if (!data || sz == 0) return NULL;
    SvtMetadataT *metadata = (SvtMetadataT *)malloc(sizeof(SvtMetadataT));
    if (!metadata) return NULL;
    metadata->type = type;
    metadata->payload = (uint8_t *)malloc(sz);
    if (!metadata->payload) {
        free(metadata);
        return NULL;
    }
    memcpy(metadata->payload, data, sz);
    metadata->sz = sz;
    return metadata;
}

EB_API void svt_metadata_free(void *ptr) {
    SvtMetadataT **metadata = (SvtMetadataT **)ptr;
    if (*metadata) {
        if ((*metadata)->payload) {
            free((*metadata)->payload);
            (*metadata)->payload = NULL;
        }
        free(*metadata);
        *metadata = NULL;
    }
}

EB_API SvtMetadataArrayT *svt_metadata_array_alloc(const size_t sz) {
    SvtMetadataArrayT *arr = (SvtMetadataArrayT *)calloc(1, sizeof(SvtMetadataArrayT));
    if (!arr) return NULL;
    if (sz > 0) {
        arr->metadata_array = (SvtMetadataT **)calloc(sz, sizeof(SvtMetadataT *));
        if (!arr->metadata_array) {
            svt_metadata_array_free(&arr);
            return NULL;
        }
        arr->sz = sz;
    }
    return arr;
}

EB_API void svt_metadata_array_free(void *arr) {
    SvtMetadataArrayT **metadata = (SvtMetadataArrayT **)arr;
    if (*metadata) {
        if ((*metadata)->metadata_array) {
            for (size_t i = 0; i < (*metadata)->sz; i++) {
                svt_metadata_free(&((*metadata)->metadata_array[i]));
            }
            free((*metadata)->metadata_array);
        }
        free(*metadata);
    }
    *metadata = NULL;
}

EB_API int svt_add_metadata(EbBufferHeaderType *buffer, const uint32_t type, const uint8_t *data, const size_t sz) {
    if (!buffer) return -1;
    if (!buffer->metadata) {
        buffer->metadata = svt_metadata_array_alloc(0);
        if (!buffer->metadata) return -1;
    }
    SvtMetadataT *metadata = svt_metadata_alloc(type, data, sz);
    if (!metadata) return -1;
    SvtMetadataT **metadata_array =
        (SvtMetadataT **)realloc(buffer->metadata->metadata_array, (buffer->metadata->sz + 1) * sizeof(metadata));
    if (!metadata_array) {
        svt_metadata_free(&metadata);
        return -1;
    }
    buffer->metadata->metadata_array = metadata_array;
    buffer->metadata->metadata_array[buffer->metadata->sz] = metadata;
    buffer->metadata->sz++;
    return 0;
}

EB_API size_t svt_metadata_size(SvtMetadataArrayT *metadata, const EbAv1MetadataType type) {
    size_t sz = 0;
    if (!metadata || !metadata->metadata_array || metadata->sz == 0) {
        return 0;
    } else {
        for (size_t i = 0; i < metadata->sz; i++) {
            SvtMetadataT *current_metadata = metadata->metadata_array[i];
            if (current_metadata && current_metadata->payload && current_metadata->type == type) {
                sz = current_metadata->sz + 1 //obu type
                    + 1 //trailing byte
                    + 1 //header size
                    + 1; //length field size
            }
        }
    }
    return sz;
}
