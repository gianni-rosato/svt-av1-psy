/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file TestEnv.c
 *
 * @brief Environment setup for unit test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "aom_dsp_rtcd.h"

/** setup_test_env is a util for unit test setup environment without create a
 * encoder */
void setup_test_env(EbAsm asm_type) {
    setup_rtcd_internal(asm_type);
}
