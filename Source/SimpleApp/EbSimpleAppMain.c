/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// main.cpp
//  -Contructs the following resources needed during the encoding process
//      -memory
//      -threads
//  -Configures the encoder
//  -Calls the encoder via the API
//  -Destructs the resources

/***************************************
 * Includes
 ***************************************/
#include "EbSimpleAppContext.h"
#include "EbApi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#if _WIN32
#define fseeko64 _fseeki64
#define ftello64 _ftelli64
#define FOPEN(f,s,m) fopen_s(&f,s,m)
#else
#define fseeko64 fseek
#define ftello64 ftell
#define FOPEN(f,s,m) f=fopen(s,m)
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#include <errno.h>
#endif

/** The APPEXITCONDITIONTYPE type is used to define the App main loop exit
conditions.
*/
typedef enum APPEXITCONDITIONTYPE {
    APP_ExitConditionNone = 0,
    APP_ExitConditionFinished,
    APP_ExitConditionError
} APPEXITCONDITIONTYPE;

/****************************************
* Padding
****************************************/
#define LEFT_INPUT_PADDING 0
#define RIGHT_INPUT_PADDING 0
#define TOP_INPUT_PADDING 0
#define BOTTOM_INPUT_PADDING 0

 /**********************************
 * Constructor
 **********************************/
static void EbConfigCtor(EbConfig_t *config_ptr)
{
    config_ptr->inputFile = NULL;
    config_ptr->bitstreamFile = NULL;
    config_ptr->reconFile = NULL;
    config_ptr->encoderBitDepth = 8;
    config_ptr->compressedTenBitFormat = 0;
    config_ptr->sourceWidth = 0;
    config_ptr->sourceHeight = 0;
    config_ptr->framesToBeEncoded = 0;
    config_ptr->channel_id = 0;
    config_ptr->stopEncoder = 0;

    return;
}

/**********************************
* Destructor
**********************************/
static void EbConfigDtor(EbConfig_t *config_ptr)
{

    if (config_ptr->inputFile) {
        fclose(config_ptr->inputFile);
        config_ptr->inputFile = (FILE *)NULL;
    }

    if (config_ptr->reconFile) {
        fclose(config_ptr->reconFile);
        config_ptr->reconFile = (FILE *)NULL;
    }

    if (config_ptr->bitstreamFile) {
        fclose(config_ptr->bitstreamFile);
        config_ptr->bitstreamFile = (FILE *)NULL;
    }

    return;
}

APPEXITCONDITIONTYPE ProcessOutputReconBuffer(
    EbConfig_t             *config,
    EbAppContext_t         *appCallBack)
{
    EbBufferHeaderType    *headerPtr = appCallBack->recon_buffer; // needs to change for buffered input
    EbComponentType       *componentHandle = (EbComponentType*)appCallBack->svtEncoderHandle;
    APPEXITCONDITIONTYPE    return_value = APP_ExitConditionNone;
    EbErrorType            recon_status = EB_ErrorNone;
    int32_t fseekReturnVal;
    // non-blocking call until all input frames are sent
    recon_status = eb_svt_get_recon(componentHandle, headerPtr);

    if (recon_status == EB_ErrorMax) {
        printf("\nError while outputing recon, code 0x%x\n", headerPtr->flags);
        return APP_ExitConditionError;
    }
    else if (recon_status != EB_NoErrorEmptyQueue) {
        //Sets the File position to the beginning of the file.
        rewind(config->reconFile);
        uint64_t frameNum = headerPtr->pts;
        while (frameNum>0) {
            fseekReturnVal = fseeko64(config->reconFile, headerPtr->n_filled_len, SEEK_CUR);

            if (fseekReturnVal != 0) {
                printf("Error in fseeko64  returnVal %i\n", fseekReturnVal);
                return APP_ExitConditionError;
            }
            frameNum = frameNum - 1;
        }

        fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, config->reconFile);

        // Update Output Port Activity State
        return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;
    }
    return return_value;
}
APPEXITCONDITIONTYPE ProcessOutputStreamBuffer(
    EbConfig_t             *config,
    EbAppContext_t         *appCallback,
    uint8_t           pic_send_done
)
{
    EbBufferHeaderType    *headerPtr;
    EbComponentType       *componentHandle = (EbComponentType*)appCallback->svtEncoderHandle;
    APPEXITCONDITIONTYPE    return_value = APP_ExitConditionNone;
    EbErrorType            stream_status = EB_ErrorNone;
    // System performance variables
    static int64_t          frameCount = 0;

    // non-blocking call
    stream_status = eb_svt_get_packet(componentHandle, &headerPtr, pic_send_done);

    if (stream_status == EB_ErrorMax) {
        printf("\nError while encoding, code 0x%x\n", headerPtr->flags);
        return APP_ExitConditionError;
    }else if (stream_status != EB_NoErrorEmptyQueue) {
        fwrite(headerPtr->p_buffer, 1, headerPtr->n_filled_len, config->bitstreamFile);

        // Update Output Port Activity State
        return_value = (headerPtr->flags & EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : APP_ExitConditionNone;
        //printf("\b\b\b\b\b\b\b\b\b%9d", ++frameCount);
        printf("\nDecode Order:\t%ld\tdts:\t%ld\tpts:\t%ld\tSliceType:\t%d", (long int)frameCount++, (long int)headerPtr->dts , (long int)headerPtr->pts, (int)headerPtr->pic_type);

        fflush(stdout);

        // Release the output buffer
        eb_svt_release_out_buffer(&headerPtr);
    }
    return return_value;
}

#define SIZE_OF_ONE_FRAME_IN_BYTES(width, height,is16bit) ( ( ((width)*(height)*3)>>1 )<<is16bit)
void ReadInputFrames(
    EbConfig_t                  *config,
    uint8_t                      is16bit,
    EbBufferHeaderType         *headerPtr)
{

    uint64_t  readSize;
    uint32_t  inputPaddedWidth = config->inputPaddedWidth;
    uint32_t  inputPaddedHeight = config->inputPaddedHeight;
    FILE   *inputFile = config->inputFile;
    uint8_t  *ebInputPtr;
    EbSvtEncInput* inputPtr = (EbSvtEncInput*)headerPtr->p_buffer;
    inputPtr->yStride  = inputPaddedWidth;
    inputPtr->cbStride = inputPaddedWidth >> 1;
    inputPtr->crStride = inputPaddedWidth >> 1;
    {
        if (is16bit == 0 || (is16bit == 1 && config->compressedTenBitFormat == 0)) {

            readSize = (uint64_t)SIZE_OF_ONE_FRAME_IN_BYTES(inputPaddedWidth, inputPaddedHeight, is16bit);

            headerPtr->n_filled_len = 0;

            {
                uint64_t lumaReadSize = (uint64_t)inputPaddedWidth*inputPaddedHeight << is16bit;
                ebInputPtr = inputPtr->luma;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
                ebInputPtr = inputPtr->cb;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                ebInputPtr = inputPtr->cr;
                headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
                inputPtr->luma = inputPtr->luma + ((config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING) << is16bit);
                inputPtr->cb = inputPtr->cb + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << is16bit);
                inputPtr->cr = inputPtr->cr + (((config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1)) << is16bit);

                 if (readSize != headerPtr->n_filled_len) {
                    config->stopEncoder = 1;
                 }
            }
        }
        // 10-bit Compressed Unpacked Mode
        else if (is16bit == 1 && config->compressedTenBitFormat == 1) {

            // Fill the buffer with a complete frame
            headerPtr->n_filled_len = 0;

            uint64_t lumaReadSize = (uint64_t)inputPaddedWidth*inputPaddedHeight;
            uint64_t nbitlumaReadSize = (uint64_t)(inputPaddedWidth / 4)*inputPaddedHeight;

            ebInputPtr = inputPtr->luma;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize, inputFile);
            ebInputPtr = inputPtr->cb;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);
            ebInputPtr = inputPtr->cr;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, lumaReadSize >> 2, inputFile);

            inputPtr->luma = inputPtr->luma + config->inputPaddedWidth*TOP_INPUT_PADDING + LEFT_INPUT_PADDING;
            inputPtr->cb = inputPtr->cb + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);
            inputPtr->cr = inputPtr->cr + (config->inputPaddedWidth >> 1)*(TOP_INPUT_PADDING >> 1) + (LEFT_INPUT_PADDING >> 1);


            ebInputPtr = inputPtr->lumaExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize, inputFile);
            ebInputPtr = inputPtr->cbExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);
            ebInputPtr = inputPtr->crExt;
            headerPtr->n_filled_len += (uint32_t)fread(ebInputPtr, 1, nbitlumaReadSize >> 2, inputFile);

            inputPtr->lumaExt = inputPtr->lumaExt + ((config->inputPaddedWidth >> 2)*TOP_INPUT_PADDING + (LEFT_INPUT_PADDING >> 2));
            inputPtr->cbExt = inputPtr->cbExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));
            inputPtr->crExt = inputPtr->crExt + (((config->inputPaddedWidth >> 1) >> 2)*(TOP_INPUT_PADDING >> 1) + ((LEFT_INPUT_PADDING >> 1) >> 2));

            readSize = ((lumaReadSize * 3) >> 1) + ((nbitlumaReadSize * 3) >> 1);

            if (readSize != headerPtr->n_filled_len) {
                config->stopEncoder = 1;
            }

        }
    }
    // If we reached the end of file, loop over again
    if (feof(inputFile) != 0) {
        //fseek(inputFile, 0, SEEK_SET);
        config->stopEncoder = 1;
    }

    return;
}
#define  TEST_IDR 0
APPEXITCONDITIONTYPE ProcessInputBuffer(
    EbConfig_t                  *config,
    EbAppContext_t              *appCallBack)
{
    uint8_t            is16bit = (uint8_t)(config->encoderBitDepth > 8);
    EbBufferHeaderType     *headerPtr = appCallBack->inputPictureBuffer; // needs to change for buffered input
    EbComponentType        *componentHandle = (EbComponentType*)appCallBack->svtEncoderHandle;
    APPEXITCONDITIONTYPE     return_value = APP_ExitConditionNone;
    static int32_t               frameCount = 0;

    if (config->stopEncoder == 0) {
        ReadInputFrames(
            config,
            is16bit,
            headerPtr);

        if (config->stopEncoder == 0) {
            //printf ("DISP: %d", frameCount);
            // Fill in Buffers Header control data
            headerPtr->flags = 0;
            headerPtr->p_app_private = NULL;
            headerPtr->pts         = frameCount++;
            headerPtr->pic_type   = EB_INVALID_PICTURE;
#if TEST_IDR
            if (frameCount == 200)
                headerPtr->pic_type = IDR_SLICE;
            if (frameCount == 150)
                headerPtr->pic_type = I_SLICE;
#endif
            // Send the picture
            eb_svt_enc_send_picture(componentHandle, headerPtr);
        }
        else {
            EbBufferHeaderType headerPtrLast;
            headerPtrLast.n_alloc_len = 0;
            headerPtrLast.n_filled_len = 0;
            headerPtrLast.n_tick_count = 0;
            headerPtrLast.p_app_private = NULL;
            headerPtrLast.flags = EB_BUFFERFLAG_EOS;
            headerPtrLast.p_buffer = NULL;
            headerPtr->flags = EB_BUFFERFLAG_EOS;

            eb_svt_enc_send_picture(componentHandle, &headerPtrLast);
        }
        return_value = (headerPtr->flags == EB_BUFFERFLAG_EOS) ? APP_ExitConditionFinished : return_value;
    }
    return return_value;
}

/***************************************
 * Encoder App Main
 ***************************************/
int32_t main(int32_t argc, char* argv[])
{
    EbErrorType            return_error = EB_ErrorNone;            // Error Handling
    APPEXITCONDITIONTYPE    exitConditionOutput = APP_ExitConditionNone , exitConditionInput = APP_ExitConditionNone , exitConditionRecon = APP_ExitConditionNone;    // Processing loop exit condition
    EbConfig_t             *config;        // Encoder Configuration
    EbAppContext_t         *appCallback;   // Instances App callback data

    // Print Encoder Info
    printf("-------------------------------------\n");
    printf("SVT-AV1 Encoder Simple Sample Application v1.2.0\n");
    printf("Platform:   %u bit\n", (unsigned) sizeof(void*)*8);
#if ( defined( _MSC_VER ) && (_MSC_VER < 1910) )
    printf("Compiler: VS13\n");
#elif ( defined( _MSC_VER ) && (_MSC_VER >= 1910) )
    printf("Compiler: VS17\n");
#elif defined(__INTEL_COMPILER)
    printf("Compiler: Intel\n");
#elif defined(__GNUC__)
    printf("Compiler: GCC\n");
#else
    printf("Compiler: unknown\n");
#endif

    printf("APP Build date: %s %s\n",__DATE__,__TIME__);
    fflush(stdout);
    {
        // Initialize config
        config = (EbConfig_t*)malloc(sizeof(EbConfig_t));
        if (config){
            EbConfigCtor(config);
            if (argc != 6 && argc != 7) {
                printf("Usage: ./SvtHevcEncSimpleApp in.yuv out.ivf width height bitdepth recon.yuv(optional)\n");
                return_error = EB_ErrorBadParameter;
            }
            else {
                // Get info for config
                FILE * fin;
                FOPEN(fin,argv[1], "rb");
                if (!fin) {
                    printf("Invalid input file \n");
                    return_error = EB_ErrorBadParameter;
                }
                else
                    config->inputFile = fin;

                FILE * fout;
                FOPEN(fout,argv[2], "wb");
                if (!fout) {
                    printf("Invalid input file \n");
                    return_error = EB_ErrorBadParameter;
                }
                else
                    config->bitstreamFile = fout;

                uint32_t width = 0, height = 0;

                width = strtoul(argv[3], NULL, 0);
                height = strtoul(argv[4], NULL, 0);
                if ((width&&height) == 0){
                    printf("Invalid video dimensions\n");
                    return_error = EB_ErrorBadParameter;
                }

                config->inputPaddedWidth  = config->sourceWidth = width;
                config->inputPaddedHeight = config->sourceHeight = height;

                uint32_t bdepth = width = strtoul(argv[5], NULL, 0);
                if ((bdepth != 8) && (bdepth != 10)){
                    printf("Invalid bit depth\n");
                    return_error = EB_ErrorBadParameter;
                }
                config->encoderBitDepth = bdepth;

                if (argc == 7) {
                    FILE * frec;
                    FOPEN(frec, argv[6], "wb");
                    if (!frec) {
                        printf("Invalid recon file \n");
                        return_error = EB_ErrorBadParameter;
                    }
                    else
                        config->reconFile = frec;
                }

            }
        }
        if (return_error == EB_ErrorNone && (config != NULL)) {

            // Initialize appCallback
            appCallback = (EbAppContext_t*)malloc(sizeof(EbAppContext_t));
            if (appCallback){
                EbAppContextCtor(appCallback,config);

                return_error = InitEncoder(config, appCallback, 0);

                printf("Encoding          ");
                fflush(stdout);

                // Input Loop Thread
                exitConditionOutput = APP_ExitConditionNone;
                exitConditionRecon = APP_ExitConditionNone;
                while (exitConditionOutput == APP_ExitConditionNone) {
                    exitConditionInput = ProcessInputBuffer(config, appCallback);
                    if (config->reconFile) {
                        exitConditionRecon = ProcessOutputReconBuffer(config, appCallback);
                    }
                    exitConditionOutput = ProcessOutputStreamBuffer(config, appCallback, (exitConditionInput == APP_ExitConditionNone || (exitConditionRecon == APP_ExitConditionNone && config->reconFile) ? 0 : 1));
                }
                return_error = DeInitEncoder(appCallback, 0);
                // Destruct the App memory variables
                EbAppContextDtor(appCallback);
                free(appCallback);

            }else
                printf("Error allocating EbAppContext structure");

            printf("\n");
            fflush(stdout);

        }
        else {
            printf("Error in configuration, could not begin encoding! ... \n");
        }

        if (config){
            EbConfigDtor(config);
            free(config);
        }
    }
    printf("Encoder finished\n");

    return (return_error == 0) ? 0 : 1;
}
