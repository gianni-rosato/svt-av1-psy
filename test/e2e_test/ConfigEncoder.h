/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#ifndef _SVT_AV1_CONFIG_ENCODER_H_
#define _SVT_AV1_CONFIG_ENCODER_H_
#include "EbSvtAv1Enc.h"

void *create_enc_config();
void release_enc_config(void *config_ptr);
void set_enc_config(void *config, const char *name, const char *value);
int copy_enc_param(EbSvtAv1EncConfiguration *dst_enc_config, void *config_ptr);
std::string get_enc_token(const char* name);

#endif
