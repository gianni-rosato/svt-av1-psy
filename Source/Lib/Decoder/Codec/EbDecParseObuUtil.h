/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#ifndef EbDecParseObuUtil_h
#define EbDecParseObuUtil_h

#include "EbAv1Structs.h"

/* Returns information about the sequence header
 *
 * Parameter:
 * @ *obu_data             OBU data.
 * @ size                  size of the OBU data
 * @ *sequence_info        information about sequence header */
#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
EB_API EbErrorType eb_get_sequence_info(const uint8_t *obu_data, size_t size,
                                        SeqHeader *sequence_info);
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDecParseObuUtil_h
