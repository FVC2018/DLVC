/*!
*******************************************************************************
*                               JEM_ITT Encoder
*                      ITTIAM SYSTEMS PVT LTD, BANGALORE
*                               COPYRIGHT(C) 2017
*
*  This program is proprietary to Ittiam Systems Pvt. Ltd. and is protected
*  under Indian Copyright Act as an unpublished work. Its use and disclosure is
*  limited by the terms and conditions of a license agreement. It may not be
*  copied or otherwise reproduced or disclosed to persons outside the licensee's
*  organization except in accordance with the terms and conditions of such an
*  agreement. All copies and reproductions shall be the property of Ittiam
*******************************************************************************
*/

/** \file     TComFxnVariantSelector.h
    \brief    Set of macro definitions used for selecting between the C and SIMD variants of a function
*/

#pragma once

#define C_FXN 0
#define SIMD_FXN 1
#define COARSE_SIMD_FXN 2

#if ITH_OPT_C3_SIMD

#define ITH_C3_TWODIMBUFCOPIER                                  SIMD_FXN    ///< C and SIMD

#define ITH_C3_RESIDUECOMPUTER                                  SIMD_FXN    ///< C and SIMD

#define ITH_C3_RECONCOMPUTER                                    SIMD_FXN    ///< C and SIMD

#define ITH_C3_SCALEADDCLIP                                     SIMD_FXN    ///< C and SIMD

#define ITH_C3_SCALESUB                                         SIMD_FXN    ///< C and SIMD

#define ITH_C3_INTERP_FILT                                      SIMD_FXN    ///< C, SIMD and COARSE_SIMD

#define ITH_C3_GRADFILTERS                                      SIMD_FXN    ///< C and SIMD

#define ITH_C3_DEQUANT                                          SIMD_FXN    ///< C and SIMD

#define ITH_C3_HSADCOMPUTER                                     SIMD_FXN    ///< C, SIMD and COARSE_SIMD

#define ITH_C3_BIOPREDCOMPUTER                                  (!!(ITH_C4_BIO_BUG_FIX) ? SIMD_FXN : C_FXN)    ///< C and SIMD   BIO相关的汇编，JEM7.0已对BIO进行修改，该汇编可能不适用，请移植时确认一下。

#define ITH_C3_SADCOMPUTER                                      SIMD_FXN    ///< C, SIMD and COARSE_SIMD

#define ITH_C3_LEFTSHIFTSUB                                     SIMD_FXN    ///< C and SIMD

#define ITH_C3_RIGHTSHIFTCLIP                                   SIMD_FXN    ///< C and SIMD

#define ITH_C3_SMT_BL_FILT                                      (!!(ITH_C2_DIV2MUL_LUT_MOD) ? SIMD_FXN : C_FXN)    ///< C and SIMD

#define ITH_C3_EQUAL_COEFF_COMPUTER                             SIMD_FXN

#define ITH_C3_HORIZ_SOBEL_FILT                                 SIMD_FXN

#define ITH_C3_VERT_SOBEL_FILT                                  SIMD_FXN

#else

#define ITH_C3_TWODIMBUFCOPIER                                  C_FXN

#define ITH_C3_RESIDUECOMPUTER                                  C_FXN

#define ITH_C3_RECONCOMPUTER                                    C_FXN

#define ITH_C3_SCALEADDCLIP                                     C_FXN

#define ITH_C3_SCALESUB                                         C_FXN

#define ITH_C3_GRADFILTERS                                      C_FXN

#define ITH_C3_DEQUANT                                          C_FXN

#define ITH_C3_BIOPREDCOMPUTER                                  C_FXN

#define ITH_C3_LEFTSHIFTSUB                                     C_FXN

#define ITH_C3_RIGHTSHIFTCLIP                                   C_FXN

#define ITH_C3_SMT_BL_FILT                                      C_FXN

#if COM16_C806_SIMD_OPT

#define ITH_C3_INTERP_FILT                                      COARSE_SIMD_FXN

#define ITH_C3_HSADCOMPUTER                                     COARSE_SIMD_FXN

#define ITH_C3_SADCOMPUTER                                      COARSE_SIMD_FXN

#else

#define ITH_C3_INTERP_FILT                                      C_FXN

#define ITH_C3_HSADCOMPUTER                                     C_FXN

#define ITH_C3_SADCOMPUTER                                      C_FXN

#define ITH_C3_EQUAL_COEFF_COMPUTER                             C_FXN

#define ITH_C3_HORIZ_SOBEL_FILT                                 C_FXN

#define ITH_C3_VERT_SOBEL_FILT                                  C_FXN

#endif

#endif