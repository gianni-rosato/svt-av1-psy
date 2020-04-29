;
; Copyright(c) 2019 Intel Corporation
; SPDX - License - Identifier: BSD - 2 - Clause - Patent
;

%include "x64inc.asm"
%include "x64Macro.asm"
section .text
; ----------------------------------------------------------------------------------------

cglobal picture_copy_kernel_sse2

; Requirement: areaWidthInBytes = 4, 8, 12, 16, 24, 32, 48, 64 or 128
; Requirement: area_height % 2 = 0

%define src              r0
%define srcStrideInBytes r1
%define dst              r2
%define dstStrideInBytes r3
%define areaWidthInBytes r4
%define area_height       r5

    GET_PARAM_5UXD
    GET_PARAM_6UXD
    XMM_SAVE

    cmp             areaWidthInBytes,        16
    jg              Label_PictureCopyKernel_SSE2_WIDTH_Big
    je              Label_PictureCopyKernel_SSE2_WIDTH16

    cmp             areaWidthInBytes,        4
    je              Label_PictureCopyKernel_SSE2_WIDTH4

    cmp             areaWidthInBytes,        8
    je              Label_PictureCopyKernel_SSE2_WIDTH8

Label_PictureCopyKernel_SSE2_WIDTH12:
    movq            xmm0,           [src]
    movd            xmm1,           [src+8]
    movq            xmm2,           [src+srcStrideInBytes]
    movd            xmm3,           [src+srcStrideInBytes+8]
    lea             src,            [src+2*srcStrideInBytes]
    movq            [dst],          xmm0
    movd            [dst+8],        xmm1
    movq            [dst+dstStrideInBytes], xmm2
    movd            [dst+dstStrideInBytes+8], xmm3
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH12

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH8:
    movq            xmm0,           [src]
    movq            xmm1,           [src+srcStrideInBytes]
    lea             src,            [src+2*srcStrideInBytes]
    movq            [dst],          xmm0
    movq            [dst+dstStrideInBytes], xmm1
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH8

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH4:
    movd            xmm0,           [src]
    movd            xmm1,           [src+srcStrideInBytes]
    lea             src,            [src+2*srcStrideInBytes]
    movd            [dst],          xmm0
    movd            [dst+dstStrideInBytes], xmm1
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH4

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH_Big:
    cmp             areaWidthInBytes,        24
    je              Label_PictureCopyKernel_SSE2_WIDTH24

    cmp             areaWidthInBytes,        32
    je              Label_PictureCopyKernel_SSE2_WIDTH32

    cmp             areaWidthInBytes,        48
    je              Label_PictureCopyKernel_SSE2_WIDTH48

    cmp             areaWidthInBytes,        64
    je              Label_PictureCopyKernel_SSE2_WIDTH64

Label_PictureCopyKernel_SSE2_WIDTH128:
    movdqu          xmm0,           [src]
    movdqu          xmm1,           [src+16]
    movdqu          xmm2,           [src+32]
    movdqu          xmm3,           [src+48]
    movdqu          xmm4,           [src+64]
    movdqu          xmm5,           [src+80]
    movdqu          xmm6,           [src+96]
    movdqu          xmm7,           [src+112]
    lea             src,            [src+srcStrideInBytes]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+32],       xmm2
    movdqu          [dst+48],       xmm3
    movdqu          [dst+64],       xmm4
    movdqu          [dst+80],       xmm5
    movdqu          [dst+96],       xmm6
    movdqu          [dst+112],      xmm7
    lea             dst,            [dst+dstStrideInBytes]
    sub             area_height,     1
    jne             Label_PictureCopyKernel_SSE2_WIDTH128

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH64:
    movdqu          xmm0,           [src]
    movdqu          xmm1,           [src+16]
    movdqu          xmm2,           [src+32]
    movdqu          xmm3,           [src+48]
    movdqu          xmm4,           [src+srcStrideInBytes]
    movdqu          xmm5,           [src+srcStrideInBytes+16]
    movdqu          xmm6,           [src+srcStrideInBytes+32]
    movdqu          xmm7,           [src+srcStrideInBytes+48]
    lea             src,            [src+2*srcStrideInBytes]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+32],       xmm2
    movdqu          [dst+48],       xmm3
    movdqu          [dst+dstStrideInBytes],    xmm4
    movdqu          [dst+dstStrideInBytes+16], xmm5
    movdqu          [dst+dstStrideInBytes+32], xmm6
    movdqu          [dst+dstStrideInBytes+48], xmm7
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH64

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH48:
    movdqu          xmm0,           [src]
    movdqu          xmm1,           [src+16]
    movdqu          xmm2,           [src+32]
    movdqu          xmm3,           [src+srcStrideInBytes]
    movdqu          xmm4,           [src+srcStrideInBytes+16]
    movdqu          xmm5,           [src+srcStrideInBytes+32]
    lea             src,            [src+2*srcStrideInBytes]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+32],       xmm2
    movdqu          [dst+dstStrideInBytes], xmm3
    movdqu          [dst+dstStrideInBytes+16], xmm4
    movdqu          [dst+dstStrideInBytes+32], xmm5
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH48

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH32:
    movdqu          xmm0,           [src]
    movdqu          xmm1,           [src+16]
    movdqu          xmm2,           [src+srcStrideInBytes]
    movdqu          xmm3,           [src+srcStrideInBytes+16]
    lea             src,            [src+2*srcStrideInBytes]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+dstStrideInBytes], xmm2
    movdqu          [dst+dstStrideInBytes+16], xmm3
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH32

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH24:
    movdqu          xmm0,           [src]
    movq            xmm1,           [src+16]
    movdqu          xmm2,           [src+srcStrideInBytes]
    movq            xmm3,           [src+srcStrideInBytes+16]
    lea             src,            [src+2*srcStrideInBytes]
    movdqu          [dst],          xmm0
    movq            [dst+16],       xmm1
    movdqu          [dst+dstStrideInBytes], xmm2
    movq            [dst+dstStrideInBytes+16], xmm3
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH24

    XMM_RESTORE
    ret

Label_PictureCopyKernel_SSE2_WIDTH16:
    movdqu          xmm0,           [src]
    movdqu          xmm1,           [src+srcStrideInBytes]
    lea             src,            [src+2*srcStrideInBytes]
    movdqu          [dst],          xmm0
    movdqu          [dst+dstStrideInBytes], xmm1
    lea             dst,            [dst+2*dstStrideInBytes]
    sub             area_height,     2
    jne             Label_PictureCopyKernel_SSE2_WIDTH16

    XMM_RESTORE
    ret

; ----------------------------------------------------------------------------------------

cglobal picture_average_kernel_sse2

; Requirement: pu_width         = 4, 8, 12, 16, 24, 32, 48 or 64
; Requirement: pu_height   %  2 = 0
; Requirement: src0       % 16 = 0 when pu_width >= 16
; Requirement: src1       % 16 = 0 when pu_width >= 16
; Requirement: dst        % 16 = 0 when pu_width >= 16
; Requirement: src0_stride % 16 = 0 when pu_width >= 16
; Requirement: src1_stride % 16 = 0 when pu_width >= 16
; Requirement: dst_stride  % 16 = 0 when pu_width >= 16

%define src0       r0
%define src0_stride r1
%define src1       r2
%define src1_stride r3
%define dst        r4
%define dst_stride  r5
%define area_width  r6
%define area_height r7

    GET_PARAM_5Q
    GET_PARAM_6UXD
    GET_PARAM_7UXD
    PUSH_REG 7
    XMM_SAVE
    GET_PARAM_8UXD  r7d                         ; area_height

    cmp             area_width,      16
    jg              Label_PictureAverageKernel_SSE2_WIDTH_Big
    je              Label_PictureAverageKernel_SSE2_WIDTH16

    cmp             area_width,      4
    je              Label_PictureAverageKernel_SSE2_WIDTH4

    cmp             area_width,      8
    je              Label_PictureAverageKernel_SSE2_WIDTH8

Label_PictureAverageKernel_SSE2_WIDTH12:
    movq            mm0,            [src0]
    movd            mm1,            [src0+8]
    movq            mm2,            [src0+src0_stride]
    movd            mm3,            [src0+src0_stride+8]
    pavgb           mm0,            [src1]
    pavgb           mm1,            [src1+8]
    pavgb           mm2,            [src1+src1_stride]
    pavgb           mm3,            [src1+src1_stride+8]
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movq            [dst],          mm0
    movd            [dst+8],        mm1
    movq            [dst+dst_stride], mm2
    movd            [dst+dst_stride+8], mm3
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH12

    XMM_RESTORE
    POP_REG 7

%if NEED_EMMS
    emms
%endif
    ret

Label_PictureAverageKernel_SSE2_WIDTH8:
    movq            mm0,            [src0]
    movq            mm1,            [src0+src0_stride]
    pavgb           mm0,            [src1]
    pavgb           mm1,            [src1+src1_stride]
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movq            [dst],          mm0
    movq            [dst+dst_stride], mm1
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH8

    XMM_RESTORE
    POP_REG 7

%if NEED_EMMS
    emms
%endif
    ret

Label_PictureAverageKernel_SSE2_WIDTH4:
    movd            mm0,            [src0]
    movd            mm1,            [src0+src0_stride]
    pavgb           mm0,            [src1]
    pavgb           mm1,            [src1+src1_stride]
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movd            [dst],          mm0
    movd            [dst+dst_stride], mm1
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH4

    XMM_RESTORE
    POP_REG 7

%if NEED_EMMS
    emms
%endif
    ret

Label_PictureAverageKernel_SSE2_WIDTH_Big:
    cmp             area_width,      24
    je              Label_PictureAverageKernel_SSE2_WIDTH24

    cmp             area_width,      32
    je              Label_PictureAverageKernel_SSE2_WIDTH32

    cmp             area_width,      48
    je              Label_PictureAverageKernel_SSE2_WIDTH48

Label_PictureAverageKernel_SSE2_WIDTH64:
    movdqu          xmm0,           [src0]
    movdqu          xmm1,           [src0+16]
    movdqu          xmm2,           [src0+32]
    movdqu          xmm3,           [src0+48]
    movdqu          xmm4,           [src0+src0_stride]
    movdqu          xmm5,           [src0+src0_stride+16]
    movdqu          xmm6,           [src0+src0_stride+32]
    movdqu          xmm7,           [src0+src0_stride+48]
    movdqu            xmm8,            [src1]
    pavgb           xmm0,            xmm8
    movdqu            xmm8,            [src1+16]
    pavgb           xmm1,            xmm8
    movdqu            xmm8,            [src1+32]
    pavgb           xmm2,            xmm8
    movdqu            xmm8,            [src1+48]
    pavgb           xmm3,            xmm8
    movdqu            xmm8,            [src1+src1_stride]
    pavgb           xmm4,            xmm8
    movdqu            xmm8,            [src1+src1_stride+16]
    pavgb           xmm5,            xmm8
    movdqu            xmm8,            [src1+src1_stride+32]
    pavgb           xmm6,            xmm8
    movdqu            xmm8,            [src1+src1_stride+48]
    pavgb           xmm7,            xmm8
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+32],       xmm2
    movdqu          [dst+48],       xmm3
    movdqu          [dst+dst_stride], xmm4
    movdqu          [dst+dst_stride+16], xmm5
    movdqu          [dst+dst_stride+32], xmm6
    movdqu          [dst+dst_stride+48], xmm7
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH64

    XMM_RESTORE
    POP_REG 7
    ret

Label_PictureAverageKernel_SSE2_WIDTH48:
    movdqu          xmm0,           [src0]
    movdqu          xmm1,           [src0+16]
    movdqu          xmm2,           [src0+32]
    movdqu          xmm3,           [src0+src0_stride]
    movdqu          xmm4,           [src0+src0_stride+16]
    movdqu          xmm5,           [src0+src0_stride+32]
    movdqu          xmm6,            [src1]
    pavgb           xmm0,            xmm6
    movdqu          xmm6,            [src1+16]
    pavgb           xmm1,            xmm6
    movdqu          xmm6,            [src1+32]
    pavgb           xmm2,            xmm6
    movdqu          xmm6,            [src1+src1_stride]
    pavgb           xmm3,            xmm6
    movdqu          xmm6,            [src1+src1_stride+16]
    pavgb           xmm4,            xmm6
    movdqu          xmm6,            [src1+src1_stride+32]
    pavgb           xmm5,            xmm6
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+32],       xmm2
    movdqu          [dst+dst_stride], xmm3
    movdqu          [dst+dst_stride+16], xmm4
    movdqu          [dst+dst_stride+32], xmm5
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH48

    XMM_RESTORE
    POP_REG 7
    ret

Label_PictureAverageKernel_SSE2_WIDTH32:
    movdqu          xmm0,           [src0]
    movdqu          xmm1,           [src0+16]
    movdqu          xmm2,           [src0+src0_stride]
    movdqu          xmm3,           [src0+src0_stride+16]
    movdqu            xmm4,            [src1]
    pavgb           xmm0,           xmm4
    movdqu            xmm4,            [src1+16]
    pavgb           xmm1,           xmm4
    movdqu            xmm4,            [src1+src1_stride]
    pavgb           xmm2,           xmm4
    movdqu            xmm4,            [src1+src1_stride+16]
    pavgb           xmm3,           xmm4
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movdqu          [dst],          xmm0
    movdqu          [dst+16],       xmm1
    movdqu          [dst+dst_stride], xmm2
    movdqu          [dst+dst_stride+16], xmm3
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH32

    XMM_RESTORE
    POP_REG 7
    ret

Label_PictureAverageKernel_SSE2_WIDTH24:
    movdqu          xmm0,           [src0]
    movq            mm0,            [src0+16]
    movdqu          xmm1,           [src0+src0_stride]
    movq            mm1,            [src0+src0_stride+16]
    movdqu          xmm2,           [src1]
    pavgb           xmm0,           xmm2
    pavgb           mm0,            [src1+16]
    movdqu          xmm3,           [src1+src1_stride]
    pavgb           xmm1,           xmm3
    pavgb           mm1,            [src1+src1_stride+16]
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movdqu          [dst],          xmm0
    movq            [dst+16],       mm0
    movdqu          [dst+dst_stride], xmm1
    movq            [dst+dst_stride+16], mm1
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH24

    XMM_RESTORE
    POP_REG 7

%if NEED_EMMS
    emms
%endif
    ret

Label_PictureAverageKernel_SSE2_WIDTH16:
    movdqu          xmm0,           [src0]
    movdqu          xmm1,           [src0+src0_stride]
    movdqu          xmm4,           [src1]
    pavgb           xmm0,           xmm4
    movdqu          xmm4,           [src1+src1_stride]
    pavgb           xmm1,           xmm4
    lea             src0,           [src0+2*src0_stride]
    lea             src1,           [src1+2*src1_stride]
    movdqu          [dst],          xmm0
    movdqu          [dst+dst_stride], xmm1
    lea             dst,            [dst+2*dst_stride]
    sub             area_height,     2
    jne             Label_PictureAverageKernel_SSE2_WIDTH16

    XMM_RESTORE
    POP_REG 7
    ret

; ----------------------------------------------------------------------------------------
    cglobal Log2f_ASM
;   If (r0 == 0) then bsr return undefined behavior. For 0 return 0
    or r0, 1
    bsr rax, r0
    ret
