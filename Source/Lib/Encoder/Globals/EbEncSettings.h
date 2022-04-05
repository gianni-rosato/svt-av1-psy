/*
* Copyright(c) 2022 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbEncSettings_h
#define EbEncSettings_h

#include "EbSvtAv1Enc.h"
#include "EbDebugMacros.h"
#include "EbSequenceControlSet.h"

/**************************************
 * Defines
 **************************************/
#define DEFAULT_QP 35

void svt_av1_print_lib_params(SequenceControlSet* scs);

EbErrorType svt_av1_verify_settings(struct SequenceControlSet* scs_ptr);

EbErrorType svt_av1_set_default_params(EbSvtAv1EncConfiguration* config_ptr);

#endif // EbEncSettings_h
