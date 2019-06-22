/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

extern "C" {
#include "EbAppConfig.c"
#include "EbAppContext.c"
}

void set_enc_config(void *config_ptr, const char *name, const char *value) {
    EbConfig *config = (EbConfig *)config_ptr;
    SetConfigValue(config, name, value);
}

void *create_enc_config() {
    EbConfig *config = (EbConfig *)malloc(sizeof(EbConfig));
    eb_config_ctor(config);
    return config;
}

void release_enc_config(void *config_ptr) {
    if (config_ptr) {
        EbConfig *config = (EbConfig *)config_ptr;
        eb_config_dtor(config);
        free(config_ptr);
    }
}

int copy_enc_param(EbSvtAv1EncConfiguration *dst_enc_config, void *config_ptr) {
    EbAppContext app_ctx;
    EbConfig *config = (EbConfig *)config_ptr;
    // Initial app_ctx.eb_enc_parameters by copying from dst_enc_config,
    // which should be initialed in eb_init_handle
    memcpy(&app_ctx.eb_enc_parameters,
           dst_enc_config,
           sizeof(EbSvtAv1EncConfiguration));

    CopyConfigurationParameters(
        config, &app_ctx, 0);  // instance_idx is not used;

    // copy back to dst_enc_config
    memcpy(dst_enc_config,
           &app_ctx.eb_enc_parameters,
           sizeof(EbSvtAv1EncConfiguration));
    return 0;
}
