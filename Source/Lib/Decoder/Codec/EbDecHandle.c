/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbSvtAv1Dec.h"

/**********************************

* Set Parameter
**********************************/
#if defined(__linux__) || defined(__APPLE__)
__attribute__((visibility("default")))
#endif
EB_API EbErrorType eb_svt_dec_set_parameter(
    EbComponentType              *svt_dec_component,
    EbSvtAv1DecConfiguration     *pComponentParameterStructure)
{
    svt_dec_component->size = 0;
    pComponentParameterStructure->compressed_ten_bit_format=0;
    
    return EB_ErrorNone;
}