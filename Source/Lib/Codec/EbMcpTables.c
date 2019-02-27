/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbMcp.h"

/**************************************************
* Function pointer Tables
**************************************************/

// Luma
const InterpolationFilterNew uniPredLumaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][16] = {     //[ASM type][Interpolation position]
    // NON_AVX2
    {
        luma_interpolation_copy_ssse3,                        //A
        luma_interpolation_filter_posa_ssse3,                    //a
        luma_interpolation_filter_posb_ssse3,                   //b
        luma_interpolation_filter_posc_ssse3,                   //c
        luma_interpolation_filter_posd_ssse3,                   //d
        luma_interpolation_filter_pose_ssse3,                   //e
        luma_interpolation_filter_posf_ssse3,                   //f
        luma_interpolation_filter_posg_ssse3,                   //g
        luma_interpolation_filter_posh_ssse3,                   //h
        luma_interpolation_filter_posi_ssse3,                   //i
        luma_interpolation_filter_posj_ssse3,                   //j
        luma_interpolation_filter_posk_ssse3,                   //k
        luma_interpolation_filter_posn_ssse3,                   //n
        luma_interpolation_filter_posp_ssse3,                   //p
        luma_interpolation_filter_posq_ssse3,                   //q
        luma_interpolation_filter_posr_ssse3,                   //r
    },
    // AVX2
    {
        luma_interpolation_copy_ssse3,                        //A
        luma_interpolation_filter_posa_ssse3,                    //a
        luma_interpolation_filter_posb_ssse3,                   //b
        luma_interpolation_filter_posc_ssse3,                   //c
        luma_interpolation_filter_posd_ssse3,                   //d
        luma_interpolation_filter_pose_ssse3,                   //e
        luma_interpolation_filter_posf_ssse3,                   //f
        luma_interpolation_filter_posg_ssse3,                   //g
        luma_interpolation_filter_posh_ssse3,                   //h
        luma_interpolation_filter_posi_ssse3,                   //i
        luma_interpolation_filter_posj_ssse3,                   //j
        luma_interpolation_filter_posk_ssse3,                   //k
        luma_interpolation_filter_posn_ssse3,                   //n
        luma_interpolation_filter_posp_ssse3,                   //p
        luma_interpolation_filter_posq_ssse3,                   //q
        luma_interpolation_filter_posr_ssse3,                   //r
    },
};

const InterpolationFilterOutRaw biPredLumaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][16] = {     //[ASM type][Interpolation position]
        // NON_AVX2
        {
            luma_interpolation_copy_out_raw_ssse3,                   //A
            luma_interpolation_filter_posa_out_raw_ssse3,             //a
            luma_interpolation_filter_posb_out_raw_ssse3,             //b
            luma_interpolation_filter_posc_out_raw_ssse3,             //c
            luma_interpolation_filter_posd_out_raw_ssse3,             //d
            luma_interpolation_filter_pose_out_raw_ssse3,             //e
            luma_interpolation_filter_posf_out_raw_ssse3,             //f
            luma_interpolation_filter_posg_out_raw_ssse3,             //g
            luma_interpolation_filter_posh_out_raw_ssse3,             //h
            luma_interpolation_filter_posi_out_raw_ssse3,             //i
            luma_interpolation_filter_posj_out_raw_ssse3,             //j
            luma_interpolation_filter_posk_out_raw_ssse3,             //k
            luma_interpolation_filter_posn_out_raw_ssse3,             //n
            luma_interpolation_filter_posp_out_raw_ssse3,             //p
            luma_interpolation_filter_posq_out_raw_ssse3,             //q
            luma_interpolation_filter_posr_out_raw_ssse3,             //r
        },
        // AVX2
        {
            luma_interpolation_copy_out_raw_ssse3,                   //A
            luma_interpolation_filter_posa_out_raw_ssse3,             //a
            luma_interpolation_filter_posb_out_raw_ssse3,             //b
            luma_interpolation_filter_posc_out_raw_ssse3,             //c
            luma_interpolation_filter_posd_out_raw_ssse3,             //d
            luma_interpolation_filter_pose_out_raw_ssse3,             //e
            luma_interpolation_filter_posf_out_raw_ssse3,             //f
            luma_interpolation_filter_posg_out_raw_ssse3,             //g
            luma_interpolation_filter_posh_out_raw_ssse3,             //h
            luma_interpolation_filter_posi_out_raw_ssse3,             //i
            luma_interpolation_filter_posj_out_raw_ssse3,             //j
            luma_interpolation_filter_posk_out_raw_ssse3,             //k
            luma_interpolation_filter_posn_out_raw_ssse3,             //n
            luma_interpolation_filter_posp_out_raw_ssse3,             //p
            luma_interpolation_filter_posq_out_raw_ssse3,             //q
            luma_interpolation_filter_posr_out_raw_ssse3,             //r
        },
};

// Chroma
const ChromaFilterNew uniPredChromaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][64] = {
    // NON_AVX2
    {
        chroma_interpolation_copy_ssse3,                       //B
        chroma_interpolation_filter_one_d_horizontal_ssse3,         //ab
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ac
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ad
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ae
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //af
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ag
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ah
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ba
        chroma_interpolation_filter_two_d_ssse3,                 //bb
        chroma_interpolation_filter_two_d_ssse3,                 //bc
        chroma_interpolation_filter_two_d_ssse3,                 //bd
        chroma_interpolation_filter_two_d_ssse3,                 //be
        chroma_interpolation_filter_two_d_ssse3,                 //bf
        chroma_interpolation_filter_two_d_ssse3,                 //bg
        chroma_interpolation_filter_two_d_ssse3,                 //bh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ca
        chroma_interpolation_filter_two_d_ssse3,                 //cb
        chroma_interpolation_filter_two_d_ssse3,                 //cc
        chroma_interpolation_filter_two_d_ssse3,                 //cd
        chroma_interpolation_filter_two_d_ssse3,                 //ce
        chroma_interpolation_filter_two_d_ssse3,                 //cf
        chroma_interpolation_filter_two_d_ssse3,                 //cg
        chroma_interpolation_filter_two_d_ssse3,                 //ch
        chroma_interpolation_filter_one_d_vertical_ssse3,         //da
        chroma_interpolation_filter_two_d_ssse3,                 //db
        chroma_interpolation_filter_two_d_ssse3,                 //dc
        chroma_interpolation_filter_two_d_ssse3,                 //dd
        chroma_interpolation_filter_two_d_ssse3,                 //de
        chroma_interpolation_filter_two_d_ssse3,                 //df
        chroma_interpolation_filter_two_d_ssse3,                 //dg
        chroma_interpolation_filter_two_d_ssse3,                 //dh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ea
        chroma_interpolation_filter_two_d_ssse3,                 //eb
        chroma_interpolation_filter_two_d_ssse3,                 //ec
        chroma_interpolation_filter_two_d_ssse3,                 //ed
        chroma_interpolation_filter_two_d_ssse3,                 //ee
        chroma_interpolation_filter_two_d_ssse3,                 //ef
        chroma_interpolation_filter_two_d_ssse3,                 //eg
        chroma_interpolation_filter_two_d_ssse3,                 //eh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //fa
        chroma_interpolation_filter_two_d_ssse3,                 //fb
        chroma_interpolation_filter_two_d_ssse3,                 //fc
        chroma_interpolation_filter_two_d_ssse3,                 //fd
        chroma_interpolation_filter_two_d_ssse3,                 //fe
        chroma_interpolation_filter_two_d_ssse3,                 //ff
        chroma_interpolation_filter_two_d_ssse3,                 //fg
        chroma_interpolation_filter_two_d_ssse3,                 //fh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ga
        chroma_interpolation_filter_two_d_ssse3,                 //gb
        chroma_interpolation_filter_two_d_ssse3,                 //gc
        chroma_interpolation_filter_two_d_ssse3,                 //gd
        chroma_interpolation_filter_two_d_ssse3,                 //ge
        chroma_interpolation_filter_two_d_ssse3,                 //gf
        chroma_interpolation_filter_two_d_ssse3,                 //gg
        chroma_interpolation_filter_two_d_ssse3,                 //gh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ha
        chroma_interpolation_filter_two_d_ssse3,                 //hb
        chroma_interpolation_filter_two_d_ssse3,                 //hc
        chroma_interpolation_filter_two_d_ssse3,                 //hd
        chroma_interpolation_filter_two_d_ssse3,                 //he
        chroma_interpolation_filter_two_d_ssse3,                 //hf
        chroma_interpolation_filter_two_d_ssse3,                 //hg
        chroma_interpolation_filter_two_d_ssse3,                 //hh
    },
    // AVX2
    {

        chroma_interpolation_copy_ssse3,                       //B
        chroma_interpolation_filter_one_d_horizontal_ssse3,         //ab
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ac
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ad
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ae
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //af
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ag
        chroma_interpolation_filter_one_d_horizontal_ssse3,       //ah
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ba
        chroma_interpolation_filter_two_d_ssse3,                 //bb
        chroma_interpolation_filter_two_d_ssse3,                 //bc
        chroma_interpolation_filter_two_d_ssse3,                 //bd
        chroma_interpolation_filter_two_d_ssse3,                 //be
        chroma_interpolation_filter_two_d_ssse3,                 //bf
        chroma_interpolation_filter_two_d_ssse3,                 //bg
        chroma_interpolation_filter_two_d_ssse3,                 //bh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ca
        chroma_interpolation_filter_two_d_ssse3,                 //cb
        chroma_interpolation_filter_two_d_ssse3,                 //cc
        chroma_interpolation_filter_two_d_ssse3,                 //cd
        chroma_interpolation_filter_two_d_ssse3,                 //ce
        chroma_interpolation_filter_two_d_ssse3,                 //cf
        chroma_interpolation_filter_two_d_ssse3,                 //cg
        chroma_interpolation_filter_two_d_ssse3,                 //ch
        chroma_interpolation_filter_one_d_vertical_ssse3,         //da
        chroma_interpolation_filter_two_d_ssse3,                 //db
        chroma_interpolation_filter_two_d_ssse3,                 //dc
        chroma_interpolation_filter_two_d_ssse3,                 //dd
        chroma_interpolation_filter_two_d_ssse3,                 //de
        chroma_interpolation_filter_two_d_ssse3,                 //df
        chroma_interpolation_filter_two_d_ssse3,                 //dg
        chroma_interpolation_filter_two_d_ssse3,                 //dh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ea
        chroma_interpolation_filter_two_d_ssse3,                 //eb
        chroma_interpolation_filter_two_d_ssse3,                 //ec
        chroma_interpolation_filter_two_d_ssse3,                 //ed
        chroma_interpolation_filter_two_d_ssse3,                 //ee
        chroma_interpolation_filter_two_d_ssse3,                 //ef
        chroma_interpolation_filter_two_d_ssse3,                 //eg
        chroma_interpolation_filter_two_d_ssse3,                 //eh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //fa
        chroma_interpolation_filter_two_d_ssse3,                 //fb
        chroma_interpolation_filter_two_d_ssse3,                 //fc
        chroma_interpolation_filter_two_d_ssse3,                 //fd
        chroma_interpolation_filter_two_d_ssse3,                 //fe
        chroma_interpolation_filter_two_d_ssse3,                 //ff
        chroma_interpolation_filter_two_d_ssse3,                 //fg
        chroma_interpolation_filter_two_d_ssse3,                 //fh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ga
        chroma_interpolation_filter_two_d_ssse3,                 //gb
        chroma_interpolation_filter_two_d_ssse3,                 //gc
        chroma_interpolation_filter_two_d_ssse3,                 //gd
        chroma_interpolation_filter_two_d_ssse3,                 //ge
        chroma_interpolation_filter_two_d_ssse3,                 //gf
        chroma_interpolation_filter_two_d_ssse3,                 //gg
        chroma_interpolation_filter_two_d_ssse3,                 //gh
        chroma_interpolation_filter_one_d_vertical_ssse3,         //ha
        chroma_interpolation_filter_two_d_ssse3,                 //hb
        chroma_interpolation_filter_two_d_ssse3,                 //hc
        chroma_interpolation_filter_two_d_ssse3,                 //hd
        chroma_interpolation_filter_two_d_ssse3,                 //he
        chroma_interpolation_filter_two_d_ssse3,                 //hf
        chroma_interpolation_filter_two_d_ssse3,                 //hg
        chroma_interpolation_filter_two_d_ssse3,                 //hh
    },
};

const ChromaFilterOutRaw biPredChromaIFFunctionPtrArrayNew[ASM_TYPE_TOTAL][64] = {
    // NON_AVX2
    {
        chroma_interpolation_copy_out_raw_ssse3,                 //B
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ab
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ac
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ad
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ae
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //af
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ag
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ah
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,   //ba
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //bb
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bc
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bd
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //be
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bf
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bg
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,      //ca
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cb
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cc
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cd
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //ce
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cf
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cg
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //ch
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //da
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //db
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //de
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //df
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ea
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ec
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ed
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ee
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ef
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //fa
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fe
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ff
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ga
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ge
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gf
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ha
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //he
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hf
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hh
    },
    // AVX2
    {
        chroma_interpolation_copy_out_raw_ssse3,                 //B
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ab
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ac
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ad
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ae
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //af
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ag
        chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3, //ah
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,   //ba
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //bb
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bc
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bd
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //be
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bf
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bg
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //bh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,      //ca
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cb
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cc
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cd
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //ce
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cf
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //cg
        chroma_interpolation_filter_two_d_out_raw_ssse3,              //ch
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //da
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //db
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //de
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //df
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //dh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ea
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ec
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ed
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ee
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ef
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //eh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //fa
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fe
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ff
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //fh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ga
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //ge
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gf
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //gh
        chroma_interpolation_filter_one_d_out_raw_vertical_ssse3,    //ha
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hb
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hc
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hd
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //he
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hf
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hg
        chroma_interpolation_filter_two_d_out_raw_ssse3,            //hh
    },
};

