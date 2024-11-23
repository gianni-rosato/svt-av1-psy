/*
* Copyright(c) 2024 Gianni Rosato
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "definitions.h"

 /**************************************
 * Instruction Set Support
 **************************************/

#ifdef _WIN32
# include <intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* */
int64_t calculate_psy_cost_c(const TranLow *orig_coeff, const TranLow *recon_coeff,
                             int width, int height);


#ifdef __cplusplus
}  // extern "C"
#endif

// clang-format on
