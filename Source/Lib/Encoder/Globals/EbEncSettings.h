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

struct SvtAv1PrivOptions {
    /* This structure is meant to allow for the extension and shortening of features.
    * This allows library and end users to use these capabilities through svt_av1_enc_parse_parameter without requiring a wrapper extension.
    * New feature applications would no longer need recompilation and relinking. */

    /* Dynamic gop
    *
    * 0 = disable Dynamic GoP
    * 1 = enable Dynamic GoP
    *  Default is 1. */
    Bool enable_dg;

#if FTR_STARTUP_MG_SIZE
    /**
     * @brief startup_mg_size
     *
     * When enabled, a MG with specified size will be inserted after the key frame.
     * The MG size is determined by 2^startup_mg_size.
     *
     * 0: off
     * 2: set hierarchical levels to 2 (MG size 4)
     * 3: set hierarchical levels to 3 (MG size 8)
     * 4: set hierarchical levels to 4 (MG size 16)
     * Default is 0.
     */
    uint8_t startup_mg_size;
#endif
};

void svt_av1_print_lib_params(SequenceControlSet* scs);

EbErrorType svt_av1_verify_settings(struct SequenceControlSet* scs);

EbErrorType svt_av1_set_default_params(EbSvtAv1EncConfiguration* config_ptr);

#endif // EbEncSettings_h
