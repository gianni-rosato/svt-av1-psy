/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file TestEnv.c
 *
 * @brief Environment setup for unit test, since C++ compiler is not friendly
 * with macro "RTCD_C" opened, we have to separate it in a C file
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#define RTCD_C
#include "aom_dsp_rtcd.h"

/** setup_test_env is a util for unit test setup environment without create a
 * encoder */
void setup_test_env(EbAsm asm_type) {
    setup_rtcd_internal(asm_type);
}
