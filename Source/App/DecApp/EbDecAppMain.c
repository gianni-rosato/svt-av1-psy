/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

//  Dummy decoder application file
//    -- not functioning yet


/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include "EbSvtAv1Dec.h"


#ifdef _MSC_VER
#include <io.h>     /* _setmode() */
#include <fcntl.h>  /* _O_BINARY */
#endif

/***************************************
 * Decoder App Main
 ***************************************/
int32_t main()
{
#ifdef _MSC_VER
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // GLOBAL VARIABLES
    EbErrorType            return_error = EB_ErrorNone;            // Error Handling
    
    EbComponentType svt_dec;
    EbSvtAv1DecConfiguration eb_dec_parameters;
    
    return_error = eb_svt_dec_set_parameter(
                       &svt_dec,
                       &eb_dec_parameters);
    
    printf ("The decoder has not been implemented yet.");

    return (return_error == 0) ? 0 : 1;    
}
