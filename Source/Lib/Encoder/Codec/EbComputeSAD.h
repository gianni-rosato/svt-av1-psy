/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeSAD_h
#define EbComputeSAD_h

#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "EbComputeSAD_C.h"
#include "EbUtility.h"
#ifdef __cplusplus
extern "C" {
#endif

uint32_t sad_16b_kernel_c(uint16_t *src, // input parameter, source samples Ptr
                          uint32_t  src_stride, // input parameter, source stride
                          uint16_t *ref, // input parameter, reference samples Ptr
                          uint32_t  ref_stride, // input parameter, reference stride
                          uint32_t  height, // input parameter, block height (M)
                          uint32_t  width); // input parameter, block width (N)

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_h
