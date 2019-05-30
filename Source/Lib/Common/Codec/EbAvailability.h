/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAvailability_h
#define EbAvailability_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    extern EbBool is_bottom_left_available(
        uint32_t                depth,
        uint32_t                part_index);

    extern EbBool is_upper_right_available(
        uint32_t                depth,
        uint32_t                part_index);
#ifdef __cplusplus
}
#endif
#endif // EbAvailability_h
