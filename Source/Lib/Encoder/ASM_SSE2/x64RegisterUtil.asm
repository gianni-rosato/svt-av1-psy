;
; Copyright(c) 2019 Intel Corporation
;
; This source code is subject to the terms of the BSD 2 Clause License and
; the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
; was not distributed with this source code in the LICENSE file, you can
; obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
; Media Patent License 1.0 was not distributed with this source code in the
; PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
;

%include "../../Common/ASM_SSE2/x64inc.asm"

section .text

; (Copied from Intel 64 and IA-32 Architectures Software Developer's Manual)
; ("9.6.3 Using the EMMS Instruction")
; The EMMS instruction should be used in each of the following cases:
;  1. When an application using the x87 FPU instructions calls an MMX technology
;     library/DLL (use the EMMS instruction at the end of the MMX code).
;  2. When an application using MMX instructions calls a x87 FPU floating-point
;     library/DLL (use the EMMS instruction before calling the x87 FPU code).
;  3. When a switch is made between MMX code in a task or thread and other tasks or
;     threads in cooperative operating systems, unless it is certain that more MMX
;     instructions will be executed before any x87 FPU code.

; So if RunEmms() is called according to the above cases,
; then the "emms" instruction in all other assembly functions can be removed.

cglobal RunEmms
    emms
    ret

; ----------------------------------------------------------------------------------------


