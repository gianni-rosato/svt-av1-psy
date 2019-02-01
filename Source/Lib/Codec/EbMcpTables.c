/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbMcp.h"

/**************************************************
* Function Pointer Tables
**************************************************/

// Luma
const InterpolationFilterNew uniPredLumaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][16] = {     //[ASM type][Interpolation position]
    // NON_AVX2
    {
        LumaInterpolationCopy_SSSE3,                        //A
        LumaInterpolationFilterPosa_SSSE3,                    //a
        LumaInterpolationFilterPosb_SSSE3,                   //b
        LumaInterpolationFilterPosc_SSSE3,                   //c
        LumaInterpolationFilterPosd_SSSE3,                   //d
        LumaInterpolationFilterPose_SSSE3,                   //e
        LumaInterpolationFilterPosf_SSSE3,                   //f
        LumaInterpolationFilterPosg_SSSE3,                   //g
        LumaInterpolationFilterPosh_SSSE3,                   //h
        LumaInterpolationFilterPosi_SSSE3,                   //i
        LumaInterpolationFilterPosj_SSSE3,                   //j
        LumaInterpolationFilterPosk_SSSE3,                   //k
        LumaInterpolationFilterPosn_SSSE3,                   //n
        LumaInterpolationFilterPosp_SSSE3,                   //p
        LumaInterpolationFilterPosq_SSSE3,                   //q
        LumaInterpolationFilterPosr_SSSE3,                   //r
    },
    // AVX2
    {
        LumaInterpolationCopy_SSSE3,                        //A
        LumaInterpolationFilterPosa_SSSE3,                    //a
        LumaInterpolationFilterPosb_SSSE3,                   //b
        LumaInterpolationFilterPosc_SSSE3,                   //c
        LumaInterpolationFilterPosd_SSSE3,                   //d
        LumaInterpolationFilterPose_SSSE3,                   //e
        LumaInterpolationFilterPosf_SSSE3,                   //f
        LumaInterpolationFilterPosg_SSSE3,                   //g
        LumaInterpolationFilterPosh_SSSE3,                   //h
        LumaInterpolationFilterPosi_SSSE3,                   //i
        LumaInterpolationFilterPosj_SSSE3,                   //j
        LumaInterpolationFilterPosk_SSSE3,                   //k
        LumaInterpolationFilterPosn_SSSE3,                   //n
        LumaInterpolationFilterPosp_SSSE3,                   //p
        LumaInterpolationFilterPosq_SSSE3,                   //q
        LumaInterpolationFilterPosr_SSSE3,                   //r
    },
};

const InterpolationFilterOutRaw biPredLumaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][16] = {     //[ASM type][Interpolation position]
        // NON_AVX2
        {
            LumaInterpolationCopyOutRaw_SSSE3,                   //A
            LumaInterpolationFilterPosaOutRaw_SSSE3,             //a
            LumaInterpolationFilterPosbOutRaw_SSSE3,             //b
            LumaInterpolationFilterPoscOutRaw_SSSE3,             //c
            LumaInterpolationFilterPosdOutRaw_SSSE3,             //d
            LumaInterpolationFilterPoseOutRaw_SSSE3,             //e
            LumaInterpolationFilterPosfOutRaw_SSSE3,             //f
            LumaInterpolationFilterPosgOutRaw_SSSE3,             //g
            LumaInterpolationFilterPoshOutRaw_SSSE3,             //h
            LumaInterpolationFilterPosiOutRaw_SSSE3,             //i
            LumaInterpolationFilterPosjOutRaw_SSSE3,             //j
            LumaInterpolationFilterPoskOutRaw_SSSE3,             //k
            LumaInterpolationFilterPosnOutRaw_SSSE3,             //n
            LumaInterpolationFilterPospOutRaw_SSSE3,             //p
            LumaInterpolationFilterPosqOutRaw_SSSE3,             //q
            LumaInterpolationFilterPosrOutRaw_SSSE3,             //r
        },
        // AVX2
        {
            LumaInterpolationCopyOutRaw_SSSE3,                   //A
            LumaInterpolationFilterPosaOutRaw_SSSE3,             //a
            LumaInterpolationFilterPosbOutRaw_SSSE3,             //b
            LumaInterpolationFilterPoscOutRaw_SSSE3,             //c
            LumaInterpolationFilterPosdOutRaw_SSSE3,             //d
            LumaInterpolationFilterPoseOutRaw_SSSE3,             //e
            LumaInterpolationFilterPosfOutRaw_SSSE3,             //f
            LumaInterpolationFilterPosgOutRaw_SSSE3,             //g
            LumaInterpolationFilterPoshOutRaw_SSSE3,             //h
            LumaInterpolationFilterPosiOutRaw_SSSE3,             //i
            LumaInterpolationFilterPosjOutRaw_SSSE3,             //j
            LumaInterpolationFilterPoskOutRaw_SSSE3,             //k
            LumaInterpolationFilterPosnOutRaw_SSSE3,             //n
            LumaInterpolationFilterPospOutRaw_SSSE3,             //p
            LumaInterpolationFilterPosqOutRaw_SSSE3,             //q
            LumaInterpolationFilterPosrOutRaw_SSSE3,             //r
        },
};

// Chroma
const ChromaFilterNew uniPredChromaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][64] = {
    // NON_AVX2
    {
        ChromaInterpolationCopy_SSSE3,                       //B
        ChromaInterpolationFilterOneDHorizontal_SSSE3,         //ab
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ac
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ad
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ae
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //af
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ag
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ah
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ba
        ChromaInterpolationFilterTwoD_SSSE3,                 //bb
        ChromaInterpolationFilterTwoD_SSSE3,                 //bc
        ChromaInterpolationFilterTwoD_SSSE3,                 //bd
        ChromaInterpolationFilterTwoD_SSSE3,                 //be
        ChromaInterpolationFilterTwoD_SSSE3,                 //bf
        ChromaInterpolationFilterTwoD_SSSE3,                 //bg
        ChromaInterpolationFilterTwoD_SSSE3,                 //bh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ca
        ChromaInterpolationFilterTwoD_SSSE3,                 //cb
        ChromaInterpolationFilterTwoD_SSSE3,                 //cc
        ChromaInterpolationFilterTwoD_SSSE3,                 //cd
        ChromaInterpolationFilterTwoD_SSSE3,                 //ce
        ChromaInterpolationFilterTwoD_SSSE3,                 //cf
        ChromaInterpolationFilterTwoD_SSSE3,                 //cg
        ChromaInterpolationFilterTwoD_SSSE3,                 //ch
        ChromaInterpolationFilterOneDVertical_SSSE3,         //da
        ChromaInterpolationFilterTwoD_SSSE3,                 //db
        ChromaInterpolationFilterTwoD_SSSE3,                 //dc
        ChromaInterpolationFilterTwoD_SSSE3,                 //dd
        ChromaInterpolationFilterTwoD_SSSE3,                 //de
        ChromaInterpolationFilterTwoD_SSSE3,                 //df
        ChromaInterpolationFilterTwoD_SSSE3,                 //dg
        ChromaInterpolationFilterTwoD_SSSE3,                 //dh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ea
        ChromaInterpolationFilterTwoD_SSSE3,                 //eb
        ChromaInterpolationFilterTwoD_SSSE3,                 //ec
        ChromaInterpolationFilterTwoD_SSSE3,                 //ed
        ChromaInterpolationFilterTwoD_SSSE3,                 //ee
        ChromaInterpolationFilterTwoD_SSSE3,                 //ef
        ChromaInterpolationFilterTwoD_SSSE3,                 //eg
        ChromaInterpolationFilterTwoD_SSSE3,                 //eh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //fa
        ChromaInterpolationFilterTwoD_SSSE3,                 //fb
        ChromaInterpolationFilterTwoD_SSSE3,                 //fc
        ChromaInterpolationFilterTwoD_SSSE3,                 //fd
        ChromaInterpolationFilterTwoD_SSSE3,                 //fe
        ChromaInterpolationFilterTwoD_SSSE3,                 //ff
        ChromaInterpolationFilterTwoD_SSSE3,                 //fg
        ChromaInterpolationFilterTwoD_SSSE3,                 //fh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ga
        ChromaInterpolationFilterTwoD_SSSE3,                 //gb
        ChromaInterpolationFilterTwoD_SSSE3,                 //gc
        ChromaInterpolationFilterTwoD_SSSE3,                 //gd
        ChromaInterpolationFilterTwoD_SSSE3,                 //ge
        ChromaInterpolationFilterTwoD_SSSE3,                 //gf
        ChromaInterpolationFilterTwoD_SSSE3,                 //gg
        ChromaInterpolationFilterTwoD_SSSE3,                 //gh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ha
        ChromaInterpolationFilterTwoD_SSSE3,                 //hb
        ChromaInterpolationFilterTwoD_SSSE3,                 //hc
        ChromaInterpolationFilterTwoD_SSSE3,                 //hd
        ChromaInterpolationFilterTwoD_SSSE3,                 //he
        ChromaInterpolationFilterTwoD_SSSE3,                 //hf
        ChromaInterpolationFilterTwoD_SSSE3,                 //hg
        ChromaInterpolationFilterTwoD_SSSE3,                 //hh
    },
    // AVX2
    {

        ChromaInterpolationCopy_SSSE3,                       //B
        ChromaInterpolationFilterOneDHorizontal_SSSE3,         //ab
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ac
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ad
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ae
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //af
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ag
        ChromaInterpolationFilterOneDHorizontal_SSSE3,       //ah
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ba
        ChromaInterpolationFilterTwoD_SSSE3,                 //bb
        ChromaInterpolationFilterTwoD_SSSE3,                 //bc
        ChromaInterpolationFilterTwoD_SSSE3,                 //bd
        ChromaInterpolationFilterTwoD_SSSE3,                 //be
        ChromaInterpolationFilterTwoD_SSSE3,                 //bf
        ChromaInterpolationFilterTwoD_SSSE3,                 //bg
        ChromaInterpolationFilterTwoD_SSSE3,                 //bh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ca
        ChromaInterpolationFilterTwoD_SSSE3,                 //cb
        ChromaInterpolationFilterTwoD_SSSE3,                 //cc
        ChromaInterpolationFilterTwoD_SSSE3,                 //cd
        ChromaInterpolationFilterTwoD_SSSE3,                 //ce
        ChromaInterpolationFilterTwoD_SSSE3,                 //cf
        ChromaInterpolationFilterTwoD_SSSE3,                 //cg
        ChromaInterpolationFilterTwoD_SSSE3,                 //ch
        ChromaInterpolationFilterOneDVertical_SSSE3,         //da
        ChromaInterpolationFilterTwoD_SSSE3,                 //db
        ChromaInterpolationFilterTwoD_SSSE3,                 //dc
        ChromaInterpolationFilterTwoD_SSSE3,                 //dd
        ChromaInterpolationFilterTwoD_SSSE3,                 //de
        ChromaInterpolationFilterTwoD_SSSE3,                 //df
        ChromaInterpolationFilterTwoD_SSSE3,                 //dg
        ChromaInterpolationFilterTwoD_SSSE3,                 //dh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ea
        ChromaInterpolationFilterTwoD_SSSE3,                 //eb
        ChromaInterpolationFilterTwoD_SSSE3,                 //ec
        ChromaInterpolationFilterTwoD_SSSE3,                 //ed
        ChromaInterpolationFilterTwoD_SSSE3,                 //ee
        ChromaInterpolationFilterTwoD_SSSE3,                 //ef
        ChromaInterpolationFilterTwoD_SSSE3,                 //eg
        ChromaInterpolationFilterTwoD_SSSE3,                 //eh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //fa
        ChromaInterpolationFilterTwoD_SSSE3,                 //fb
        ChromaInterpolationFilterTwoD_SSSE3,                 //fc
        ChromaInterpolationFilterTwoD_SSSE3,                 //fd
        ChromaInterpolationFilterTwoD_SSSE3,                 //fe
        ChromaInterpolationFilterTwoD_SSSE3,                 //ff
        ChromaInterpolationFilterTwoD_SSSE3,                 //fg
        ChromaInterpolationFilterTwoD_SSSE3,                 //fh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ga
        ChromaInterpolationFilterTwoD_SSSE3,                 //gb
        ChromaInterpolationFilterTwoD_SSSE3,                 //gc
        ChromaInterpolationFilterTwoD_SSSE3,                 //gd
        ChromaInterpolationFilterTwoD_SSSE3,                 //ge
        ChromaInterpolationFilterTwoD_SSSE3,                 //gf
        ChromaInterpolationFilterTwoD_SSSE3,                 //gg
        ChromaInterpolationFilterTwoD_SSSE3,                 //gh
        ChromaInterpolationFilterOneDVertical_SSSE3,         //ha
        ChromaInterpolationFilterTwoD_SSSE3,                 //hb
        ChromaInterpolationFilterTwoD_SSSE3,                 //hc
        ChromaInterpolationFilterTwoD_SSSE3,                 //hd
        ChromaInterpolationFilterTwoD_SSSE3,                 //he
        ChromaInterpolationFilterTwoD_SSSE3,                 //hf
        ChromaInterpolationFilterTwoD_SSSE3,                 //hg
        ChromaInterpolationFilterTwoD_SSSE3,                 //hh
    },
};

const ChromaFilterOutRaw biPredChromaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][64] = {
    // NON_AVX2
    {
        ChromaInterpolationCopyOutRaw_SSSE3,                 //B
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ab
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ac
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ad
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ae
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //af
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ag
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ah
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,   //ba
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //bb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //be
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,      //ca
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //ce
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //ch
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //da
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //db
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //de
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //df
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ea
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ec
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ed
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ee
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ef
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //fa
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fe
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ff
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ga
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ge
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ha
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //he
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hh
    },
    // AVX2
    {
        ChromaInterpolationCopyOutRaw_SSSE3,                 //B
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ab
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ac
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ad
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ae
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //af
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ag
        ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3, //ah
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,   //ba
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //bb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //be
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //bh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,      //ca
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //ce
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //cg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,              //ch
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //da
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //db
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //de
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //df
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //dh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ea
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ec
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ed
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ee
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ef
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //eh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //fa
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fe
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ff
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //fh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ga
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //ge
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //gh
        ChromaInterpolationFilterOneDOutRawVertical_SSSE3,    //ha
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hb
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hc
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hd
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //he
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hf
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hg
        ChromaInterpolationFilterTwoDOutRaw_SSSE3,            //hh
    },
};

