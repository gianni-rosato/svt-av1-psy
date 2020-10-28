/*
* Copyright(c) 2019 Netflix, Inc.
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
EB_API EbErrorType svt_get_sequence_info(const uint8_t *obu_data, size_t size,
                                         SeqHeader *sequence_info);
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDecParseObuUtil_h
