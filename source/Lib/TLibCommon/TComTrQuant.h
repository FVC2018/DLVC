/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2015, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TComTrQuant.h
    \brief    transform and quantization class (header)
*/

#ifndef __TCOMTRQUANT__
#define __TCOMTRQUANT__

#include "CommonDef.h"
#include "TComYuv.h"
#include "TComDataCU.h"
#include "TComChromaFormat.h"
#include "ContextTables.h"
#if JVET_F0096_BILATERAL_FILTER
#include "TLibCommon/TComBilateralFilter.h"
#endif
#if PVC_AQS
#include "CommonImageProcessFunctions.h"
#endif
//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constants
// ====================================================================================================================

#define QP_BITS                 15
#if VCEG_AZ08_KLT_COMMON
#define MAX_1DTRANS_LEN         (1 << (((USE_MORE_BLOCKSIZE_DEPTH_MAX) + 1) << 1)) ///< 4x4 = 16, 8x8 = 64, 16x16=256, 32x32 = 1024
extern UInt g_uiDepth2Width[5];
#if VCEG_AZ08_INTER_KLT
extern UInt g_uiDepth2InterTempSize[5];
#endif
#if VCEG_AZ08_INTRA_KLT
extern UInt g_uiDepth2IntraTempSize[5];
#endif
extern UInt g_uiDepth2MaxCandiNum[5];
extern UInt g_uiDepth2MinCandiNum[5];
#endif
// ====================================================================================================================
// Type definition
// ====================================================================================================================
#if VCEG_AZ05_ROT_TR
 __inline static Short  xMult ( Int i, UInt uiShift ) { return ((i)>>uiShift); }
#endif
#if COM16_C806_EMT
#if JVET_D0077_TRANSFORM_OPT
#if JVET_C0024_QTBT
void fastForwardDCT2_B2    (TCoeff *block, TCoeff *coeff, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B2    (TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B128  (TCoeff *block, TCoeff *coeff, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B128  (TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif
void fastForwardDCT2_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT2_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);

void fastForwardDST7_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDST7_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST7_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDCT5_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDCT5_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT5_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDCT8_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDCT8_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDCT8_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDST1_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDST1_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use);
void fastInverseDST1_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

typedef void FwdTrans (TCoeff *, TCoeff *, Int, Int, Int, Int, Int);
typedef void InvTrans (TCoeff *, TCoeff *, Int, Int, Int, Int, Int, const TCoeff, const TCoeff);

#else//JVET_D0077_TRANSFORM_OPT

#if JVET_C0024_QTBT
void fastForwardDCT2_B2 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT2_B2 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT2_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif
void fastForwardDCT2_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT2_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT2_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT2_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT2_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT2_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT2_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);

void fastForwardDST7_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST7_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);

void fastForwardDST7_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST7_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST7_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST7_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDST7_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDST7_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST7_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDST7_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDCT5_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT5_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT5_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT5_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT5_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDCT5_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT5_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT5_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT5_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDCT8_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT8_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT8_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT8_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDCT8_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDCT8_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT8_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDCT8_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDCT8_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

void fastForwardDST1_B4 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST1_B4 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B8 (TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST1_B8 (TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B16(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST1_B16(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B32(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int zo, Int use);
void fastInverseDST1_B32(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum);
#if JVET_C0024_QTBT
void fastForwardDST1_B64(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDST1_B64(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
void fastForwardDST1_B128(TCoeff *block,TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2);
void fastInverseDST1_B128(TCoeff *coeff,TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum);
#endif

typedef void FwdTrans (TCoeff *, TCoeff *, Int, Int, Int, Int);
typedef void InvTrans (TCoeff *, TCoeff *, Int, Int, Int, Int, const TCoeff, const TCoeff);
#endif//JVET_D0077_TRANSFORM_OPT
#endif

typedef struct
{
  Int significantCoeffGroupBits[NUM_SIG_CG_FLAG_CTX][2 /*Flag = [0|1]*/];
  Int significantBits[NUM_SIG_FLAG_CTX][2 /*Flag = [0|1]*/];
  Int lastXBits[MAX_NUM_CHANNEL_TYPE][LAST_SIGNIFICANT_GROUPS];
  Int lastYBits[MAX_NUM_CHANNEL_TYPE][LAST_SIGNIFICANT_GROUPS];
  Int m_greaterOneBits[NUM_ONE_FLAG_CTX][2 /*Flag = [0|1]*/];
#if !VCEG_AZ07_CTX_RESIDUALCODING
  Int m_levelAbsBits[NUM_ABS_FLAG_CTX][2 /*Flag = [0|1]*/];
#endif
  Int blockCbpBits[NUM_QT_CBF_CTX_SETS * NUM_QT_CBF_CTX_PER_SET][2 /*Flag = [0|1]*/];
  Int blockRootCbpBits[4][2 /*Flag = [0|1]*/];

  Int golombRiceAdaptationStatistics[RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS];
} estBitsSbacStruct;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// QP struct
struct QpParam
{
  Int Qp;
  Int per;
  Int rem;

  QpParam(const Int           qpy,
          const ChannelType   chType,
          const Int           qpBdOffset,
          const Int           chromaQPOffset,
          const ChromaFormat  chFmt );

  QpParam(const TComDataCU   &cu, const ComponentID compID);

}; // END STRUCT DEFINITION QpParam

#if VCEG_AZ08_KLT_COMMON
class TempLibFast
{
public:
    Int *m_pX;    //offset X
    Int *m_pY;    //offset Y
    Int *m_pXInteger;    //offset X for integer pixel search
    Int *m_pYInteger;    //offset Y for integer pixel search
    DistType *m_pDiffInteger;
    Int* getXInteger() { return m_pXInteger; }
    Int* getYInteger() { return m_pYInteger; }
    DistType* getDiffInteger() { return m_pDiffInteger; }
    Short *m_pIdInteger; //frame id
    Short* getIdInteger() { return m_pIdInteger; }
    DistType *m_pDiff; //mse
    Short *m_pId; //frame id
    Int m_iSize;

    TempLibFast();
    ~TempLibFast();
    Void init(UInt iSize);
    Int* getX() { return m_pX; }
    Int* getY() { return m_pY; }
    DistType* getDiff() { return m_pDiff; }
    Short* getId() { return m_pId; }
    Void initDiff(UInt uiPatchSize, Int bitDepth);
    Void initDiff(UInt uiPatchSize, Int bitDepth, Int iCandiNumber);
#if VCEG_AZ08_INTRA_KLT
    Void initTemplateDiff(UInt uiPatchSize, UInt uiBlkSize, Int bitDepth, Int iCandiNumber);
#endif
    Int m_diffMax;
    Int getDiffMax() { return m_diffMax; }
};

typedef Short TrainDataType; 
typedef Double TrainDataTypeD;
#endif

/// transform and quantization class
class TComTrQuant
{
public:
  TComTrQuant();
  ~TComTrQuant();

#if PVC_AQS
  Int    m_iEstQSFactor;
  Int    getEstQSFactor() { return m_iEstQSFactor; }
  Void   setEstQSFactor(Pel* pPred, UInt uiStride, UInt uiWidth, UInt uiHeight, SliceType eSliceType, TComDataCU* pcCU, UInt uiAbsPartIdx, ComponentID comId = COMPONENT_Y);
#endif

  // initialize class
  Void init                 ( UInt  uiMaxTrSize,
#if VCEG_AZ08_USE_KLT
                              UInt  uiUseKLT              = 0,
#endif
                              Bool useRDOQ                = false,
                              Bool useRDOQTS              = false,
#if T0196_SELECTIVE_RDOQ
                              Bool useSelectiveRDOQ       = false,
#endif
                              Bool bEnc                   = false,
                              Bool useTransformSkipFast   = false
#if ADAPTIVE_QP_SELECTION
                            , Bool bUseAdaptQpSelect      = false
#endif
                              );

#if COM16_C806_EMT
  UChar getEmtTrIdx ( TComTU &rTu, const ComponentID compID );
  UChar getEmtMode  ( TComTU &rTu, const ComponentID compID );
#endif
#if VCEG_AZ08_USE_KLT
  UInt getUseKLT() { return m_useKLT;}
  UInt getUseIntraKLT() { return m_useKLT & 1; }
  UInt getUseInterKLT() { return (m_useKLT >> 1) & 1; }
#endif
#if VCEG_AZ05_ROT_TR
Void InvRotTransform4I(  Int* matrix, UChar index );
Void RotTransform4I( Int* matrix, UChar index );
#elif COM16_C1044_NSST
Void FwdNsst4x4( Int* src, UInt uiMode, UChar index );
Void InvNsst4x4( Int* src, UInt uiMode, UChar index );
#if JVET_D0120_NSST_IMPROV
Void FwdNsst8x8( Int* src, UInt uiMode, UChar index );
Void InvNsst8x8( Int* src, UInt uiMode, UChar index );
#endif
#endif 
  // transform & inverse transform functions
  Void transformNxN(       TComTU         & rTu,
                     const ComponentID      compID,
                           Pel           *  pcResidual,
                     const UInt             uiStride,
                           TCoeff        *  rpcCoeff,
#if ADAPTIVE_QP_SELECTION
                           TCoeff        * rpcArlCoeff,
#endif
                           TCoeff         & uiAbsSum,
                     const QpParam        & cQP 
#if VCEG_AZ08_KLT_COMMON
                     , Bool useKLT = false
#endif
                     );


  Void invTransformNxN(      TComTU       & rTu,
                       const ComponentID    compID,
                             Pel         *pcResidual,
                       const UInt           uiStride,
                             TCoeff      *  pcCoeff,
                       const QpParam      & cQP
#if VCEG_AZ08_KLT_COMMON
                       , Bool useKLT = false
#endif
                             DEBUG_STRING_FN_DECLAREP(psDebug) 
                             );

  Void invRecurTransformNxN ( const ComponentID compID, TComYuv *pResidual, TComTU &rTu 
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
  , TComYuv *pPred
#endif
  );

  Void rdpcmNxN   ( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride, const QpParam& cQP, TCoeff* pcCoeff, TCoeff &uiAbsSum, RDPCMMode& rdpcmMode );
  Void invRdpcmNxN( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride );

  Void applyForwardRDPCM( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride, const QpParam& cQP, TCoeff* pcCoeff, TCoeff &uiAbsSum, const RDPCMMode mode );

  // Misc functions

#if RDOQ_CHROMA_LAMBDA
  Void setLambdas(const Double lambdas[MAX_NUM_COMPONENT]) { for (UInt component = 0; component < MAX_NUM_COMPONENT; component++) m_lambdas[component] = lambdas[component]; }
  Void selectLambda(const ComponentID compIdx) { m_dLambda = m_lambdas[compIdx]; }
#if COM16_C806_CR_FROM_CB_LAMBDA_ADJUSTMENT
  Void setLambda( Double dLambda) { m_dLambda = dLambda; }
#endif
#else
  Void setLambda(Double dLambda) { m_dLambda = dLambda;}
#endif

#if COM16_C806_CR_FROM_CB_LAMBDA_ADJUSTMENT
  Double getlambda () { return m_dLambda; }
#endif

  Void setRDOQOffset( UInt uiRDOQOffset ) { m_uiRDOQOffset = uiRDOQOffset; }

  estBitsSbacStruct* m_pcEstBitsSbac;
#if VCEG_AZ07_CTX_RESIDUALCODING
  static Int getGrtZeroCtxInc   (  TCoeff*                         pcCoeff,
                                   const Int                       scanPosition,
                                   Int                             width,
                                   Int                             height,
                                   const ChannelType               chanType
                                );
  static Int getGrtOneCtxInc    ( TCoeff*                         pcCoeff,
                                  const Int                       scanPosition,
                                  Int                             width,
                                  Int                             height,
                                  const ChannelType               chanType
                                );
  static Int getGrtTwoCtxInc    ( TCoeff*                         pcCoeff,
                                  const Int                       scanPosition,
                                  Int                             width,
                                  Int                             height,
                                  const ChannelType               chanType
                                 );
  static Int getRemainCoeffCtxInc( TCoeff*                         pcCoeff,
                                   const Int                       scanPosition,
                                   Int                             width,
                                   Int                             height
                                 );
  static Int      getSigCtxInc   ( TCoeff*                         pcCoeff,
                                   const Int                       scanPosition,
                                   Int                             width,
                                   Int                             height,
                                   const ChannelType               chanType,
                                   UInt&                           sumOne,
                                   UInt&                           sumTwo,
                                   UInt&                           sumAbs
                                  );
#else  
  static Int      calcPatternSigCtx( const UInt* sigCoeffGroupFlag, UInt uiCGPosX, UInt uiCGPosY, UInt widthInGroups, UInt heightInGroups );

  static Int      getSigCtxInc     ( Int                              patternSigCtx,
                                     const TUEntropyCodingParameters &codingParameters,
                                     const Int                        scanPosition,
                                     const Int                        log2BlockWidth,
                                     const Int                        log2BlockHeight,
                                     const ChannelType                chanType
                                    );
#endif
  static UInt getSigCoeffGroupCtxInc  (const UInt*  uiSigCoeffGroupFlag,
                                       const UInt   uiCGPosX,
                                       const UInt   uiCGPosY,
#if VCEG_AZ07_CTX_RESIDUALCODING
                                             UInt   widthInGroups,
                                             UInt   heightInGroups,
                                       COEFF_SCAN_TYPE   scanIdx
#else
                                       const UInt   widthInGroups,
                                       const UInt   heightInGroups
#endif
                                       );

  Void initScalingList                      ();
  Void destroyScalingList                   ();
#if JVET_C0024_QTBT
  Void setErrScaleCoeff    ( UInt list, UInt w, UInt h, Int qp, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths );
  Double* getErrScaleCoeff              ( UInt list, UInt w, UInt h, Int qp ) { return m_errScale             [w][h][list][qp]; };  //!< get Error Scale Coefficent
  Double& getErrScaleCoeffNoScalingList ( UInt list, UInt w, UInt h, Int qp ) { return m_errScaleNoScalingList[w][h][list][qp]; };  //!< get Error Scale Coefficent
  Int* getQuantCoeff                    ( UInt list, Int qp, UInt w, UInt h ) { return m_quantCoef            [w][h][list][qp]; };  //!< get Quant Coefficent
  Int* getDequantCoeff                  ( UInt list, Int qp, UInt w, UInt h ) { return m_dequantCoef          [w][h][list][qp]; };  //!< get DeQuant Coefficent
#else
  Void setErrScaleCoeff    ( UInt list, UInt size, Int qp, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths );
  Double* getErrScaleCoeff              ( UInt list, UInt size, Int qp ) { return m_errScale             [size][list][qp]; };  //!< get Error Scale Coefficent
  Double& getErrScaleCoeffNoScalingList ( UInt list, UInt size, Int qp ) { return m_errScaleNoScalingList[size][list][qp]; };  //!< get Error Scale Coefficent
  Int* getQuantCoeff                    ( UInt list, Int qp, UInt size ) { return m_quantCoef            [size][list][qp]; };  //!< get Quant Coefficent
  Int* getDequantCoeff                  ( UInt list, Int qp, UInt size ) { return m_dequantCoef          [size][list][qp]; };  //!< get DeQuant Coefficent
#endif
  Void setUseScalingList   ( Bool bUseScalingList){ m_scalingListEnabledFlag = bUseScalingList; };
  Bool getUseScalingList   (const UInt width, const UInt height, const Bool isTransformSkip){ return m_scalingListEnabledFlag && (!isTransformSkip || ((width == 4) && (height == 4))); };
  Void setFlatScalingList  (const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths);
#if JVET_C0024_QTBT
  Void xsetFlatScalingList ( UInt list, UInt w, UInt h, Int qp);
#else
  Void xsetFlatScalingList ( UInt list, UInt size, Int qp);
#endif
  Void xSetScalingListEnc  ( TComScalingList *scalingList, UInt list, UInt size, Int qp);
  Void xSetScalingListDec  ( const TComScalingList &scalingList, UInt list, UInt size, Int qp);
  Void setScalingList      ( TComScalingList *scalingList, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths);
  Void setScalingListDec   ( const TComScalingList &scalingList);
  Void processScalingListEnc( Int *coeff, Int *quantcoeff, Int quantScales, UInt height, UInt width, UInt ratio, Int sizuNum, UInt dc);
  Void processScalingListDec( const Int *coeff, Int *dequantcoeff, Int invQuantScales, UInt height, UInt width, UInt ratio, Int sizuNum, UInt dc);
#if ADAPTIVE_QP_SELECTION
  Void    initSliceQpDelta() ;
  Void    storeSliceQpNext(TComSlice* pcSlice);
  Void    clearSliceARLCnt();
  Int     getQpDelta(Int qp) { return m_qpDelta[qp]; }
  Int*    getSliceNSamples(){ return m_sliceNsamples ;}
  Double* getSliceSumC()    { return m_sliceSumC; }
#endif
  Void transformSkipQuantOneSample(TComTU &rTu, const ComponentID compID, const TCoeff resiDiff, TCoeff* pcCoeff, const UInt uiPos, const QpParam &cQP, const Bool bUseHalfRoundingPoint);
  Void invTrSkipDeQuantOneSample(TComTU &rTu, ComponentID compID, TCoeff pcCoeff, Pel &reconSample, const QpParam &cQP, UInt uiPos );

#if VCEG_AZ08_KLT_COMMON
  Void calcCovMatrix(TrainDataType **pData, UInt uiSampleNum, covMatrixType *pCovMatrix, UInt uiDim, DistType *pDiff);
  Void calcCovMatrixXXt(TrainDataType **pData, UInt uiSampleNum, covMatrixType *pCovMatrix, UInt uiDim);
  DistType calcTemplateDiff(Pel *ref, UInt uiStride, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, DistType iMax);
  Void calcCovMatrix(TrainDataType **pData, UInt uiSampleNum, covMatrixType *pCovMatrix, UInt uiDim);
  Bool deriveKLT(UInt uiBlkSize, UInt uiUseCandiNumber);
  Bool derive1DimKLT_Fast(UInt uiBlkSize, UInt uiUseCandiNumber);
  Bool derive1DimKLT(UInt uiBlkSize, UInt uiUseCandiNumber);
  Bool derive2DimKLT(UInt uiBlkSize, DistType *pDiff);
  Pel  **getTargetPatch(UInt uiDepth) { return m_pppTarPatch[uiDepth]; }
  Pel* getRefPicUsed(UInt uiId) { return m_refPicUsed[uiId]; }
  Void setRefPicUsed(UInt uiId, Pel *ref) { m_refPicUsed[uiId] = ref; }
  UInt getStride() { return m_uiPicStride; }
  Void setStride(UInt uiPicStride) { m_uiPicStride = uiPicStride; }
#endif
#if VCEG_AZ08_INTRA_KLT
  Void searchCandidateFromOnePicIntra(TComDataCU *pcCU, UInt uiPartAddr, TComPic* refPicSrc, TComPicYuv *refPic, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt setId);
  Void candidateSearchIntra(TComDataCU *pcCU, UInt uiPartAddr, UInt uiBlkSize, UInt uiTempSize);
  Bool generateTMPrediction(Pel *piPred, UInt uiStride, UInt uiBlkSize, UInt uiTempSize, Int &foundCandiNum);
  Void getTargetTemplate(TComDataCU *pcCU, UInt uiAbsPartIdx, UInt uiBlkSize, UInt uiTempSize);
  Bool calcKLTIntra(Pel *piPred, UInt uiStride, UInt uiBlkSize);
  Bool prepareKLTSamplesIntra(Pel *piPred, UInt uiStride, UInt uiBlkSize);
#endif
#if VCEG_AZ08_INTER_KLT
  Void getTargetPatch(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absTUPartIdx, TComYuv* pcPred, UInt uiBlkSize, UInt uiTempSize);
  Void candidateSearch(TComDataCU *pcCU, UInt uiPartAddr, UInt uiBlkSize, UInt uiTempSize);
  Void searchCandidateFromOnePicInteger(TComDataCU *pcCU, UInt uiPartAddr, TComPicYuv *refPic, TComMv  cMv, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt setId);
  Void searchCandidateFraBasedOnInteger(TComDataCU *pcCU, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt uiPartAddr, Short setIdFraStart);
  Void RecordPosition(UInt uiTargetCandiNum);
  Bool candidateTrain(UInt uiBlkSize, UInt uiTempSize);
  Bool prepareKLTSamplesInter(UInt uiBlkSize, UInt uiTempSize);
  Void setRefPicBuf(UInt uiId, TComPic *refPic) { m_refPicBuf[uiId] = refPic; }
  TComPic* getRefPicBuf(UInt uiId) { return m_refPicBuf[uiId]; }
  DistType calcPatchDiff(Pel *ref, UInt uiStride, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, DistType iMax);
  Void xSetSearchRange(TComDataCU* pcCU, TComMv& cMvPred, Int iSrchRng, TComMv& rcMvSrchRngLT, TComMv& rcMvSrchRngRB);
#endif

protected:
#if ADAPTIVE_QP_SELECTION
  Int     m_qpDelta[MAX_QP+1];
  Int     m_sliceNsamples[LEVEL_RANGE+1];
  Double  m_sliceSumC[LEVEL_RANGE+1] ;
#endif
  TCoeff* m_plTempCoeff;

//  QpParam  m_cQP; - removed - placed on the stack.
#if RDOQ_CHROMA_LAMBDA
  Double   m_lambdas[MAX_NUM_COMPONENT];
#endif
  Double   m_dLambda;
  UInt     m_uiRDOQOffset;
  UInt     m_uiMaxTrSize;
  Bool     m_bEnc;
  Bool     m_useRDOQ;
  Bool     m_useRDOQTS;
#if VCEG_AZ08_USE_KLT
  UInt     m_useKLT;
#endif
#if T0196_SELECTIVE_RDOQ
  Bool     m_useSelectiveRDOQ;
#endif
#if ADAPTIVE_QP_SELECTION
  Bool     m_bUseAdaptQpSelect;
#endif
  Bool     m_useTransformSkipFast;

  Bool     m_scalingListEnabledFlag;

#if JVET_C0024_QTBT
  Int      *m_quantCoef            [SCALING_LIST_SIZE_NUM][SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
  Int      *m_dequantCoef          [SCALING_LIST_SIZE_NUM][SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of dequantization matrix coefficient 4x4
  Double   *m_errScale             [SCALING_LIST_SIZE_NUM][SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
  Double    m_errScaleNoScalingList[SCALING_LIST_SIZE_NUM][SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
#else
  Int      *m_quantCoef            [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
  Int      *m_dequantCoef          [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of dequantization matrix coefficient 4x4
  Double   *m_errScale             [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
  Double    m_errScaleNoScalingList[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM]; ///< array of quantization matrix coefficient 4x4
#endif

#if JVET_C0024_QTBT && ITH_C1_RDOQ_MEM_OPT
																									/* scratch buffers for RDOQ */
  Int       *m_piRateIncUp;
  Int       *m_piRateIncDown;
  Int       *m_piSigRateDelta;
  Double    *m_pdCostCoeff;
  Double    *m_pdCostSig;
  Double    *m_pdCostCoeff0;
  TCoeff    *m_pcDeltaU;
  /*Note: Although not directly related to QTBT, this buffers were removed from stack to avoid stackoverflow, possibly due to extreme stack usage during QTBT recursion */
#endif

#if VCEG_AZ08_KLT_COMMON
  Int m_uiPartLibSize;
  TempLibFast m_tempLibFast;
  Pel *m_refPicUsed[MAX_NUM_REF_IDS];
  TComPic *m_refPicBuf[MAX_NUM_REF_IDS];
  UInt m_uiPicStride;
  TrainDataType *m_pData[MAX_CANDI_NUM];
#if VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY
  TrainDataType *m_pDataT[MAX_1DTRANS_LEN];
#endif
  UInt m_uiVaildCandiNum;
  Double m_pEigenValues[MAX_1DTRANS_LEN];
  Int m_pIDTmp[MAX_1DTRANS_LEN];
  EigenType ***m_pppdEigenVector;
  Short ***m_pppsEigenVector;
  covMatrixType **m_pCovMatrix;
  Pel ***m_pppTarPatch;
#if VCEG_AZ08_FAST_DERIVE_KLT
  EigenType **m_pppdTmpEigenVector;
#endif
#endif

private:
  // forward Transform
  Void xT   ( const Int channelBitDepth, Bool useDST, Pel* piBlkResi, UInt uiStride, TCoeff* psCoeff, Int iWidth, Int iHeight, const Int maxLog2TrDynamicRange 
#if COM16_C806_EMT
    , UChar ucMode
    , UChar ucTrIdx
#endif
#if VCEG_AZ08_KLT_COMMON
    , Bool useKLT = false
#endif
    );

  // skipping Transform
  Void xTransformSkip ( Pel* piBlkResi, UInt uiStride, TCoeff* psCoeff, TComTU &rTu, const ComponentID component );

  Void signBitHidingHDQ( TCoeff* pQCoef, TCoeff* pCoef, TCoeff* deltaU, const TUEntropyCodingParameters &codingParameters, const Int maxLog2TrDynamicRange );

  // quantization
  Void xQuant(       TComTU       &rTu,
                     TCoeff      * pSrc,
                     TCoeff      * pDes,
#if ADAPTIVE_QP_SELECTION
                     TCoeff      *pArlDes,
#endif
                     TCoeff       &uiAbsSum,
               const ComponentID   compID,
               const QpParam      &cQP 
               );

#if T0196_SELECTIVE_RDOQ
  Bool xNeedRDOQ(    TComTU       &rTu,
                     TCoeff      * pSrc,
               const ComponentID   compID,
               const QpParam      &cQP );
#endif

  // RDOQ functions

  Void           xRateDistOptQuant (       TComTU       &rTu,
                                           TCoeff      * plSrcCoeff,
                                           TCoeff      * piDstCoeff,
#if ADAPTIVE_QP_SELECTION
                                           TCoeff      *piArlDstCoeff,
#endif
                                           TCoeff       &uiAbsSum,
                                     const ComponentID   compID,
                                     const QpParam      &cQP 
                                     );

__inline UInt              xGetCodedLevel  ( Double&          rd64CodedCost,
                                             Double&          rd64CodedCost0,
                                             Double&          rd64CodedCostSig,
                                             Intermediate_Int lLevelDouble,
                                             UInt             uiMaxAbsLevel,
                                             UShort           ui16CtxNumSig,
                                             UShort           ui16CtxNumOne,
                                             UShort           ui16CtxNumAbs,
                                             UShort           ui16AbsGoRice,
                                             UInt             c1Idx,
                                             UInt             c2Idx,
                                             Int              iQBits,
                                             Double           errorScale,
                                             Bool             bLast,
                                             Bool             useLimitedPrefixLength,
                                             const Int        maxLog2TrDynamicRange
                                             ) const;


  __inline Int xGetICRate  ( const UInt   uiAbsLevel,
                             const UShort ui16CtxNumOne,
                             const UShort ui16CtxNumAbs,
                             const UShort ui16AbsGoRice,
                             const UInt   c1Idx,
                             const UInt   c2Idx,
                             const Bool   useLimitedPrefixLength,
                             const Int maxLog2TrDynamicRange
                           ) const;

  __inline Double xGetRateLast         ( const UInt uiPosX, const UInt uiPosY, const ComponentID component ) const;
  __inline Double xGetRateSigCoeffGroup( UShort uiSignificanceCoeffGroup, UShort ui16CtxNumSig             ) const;
  __inline Double xGetRateSigCoef      ( UShort uiSignificance,           UShort ui16CtxNumSig             ) const;
  __inline Double xGetICost            ( Double dRate                                                      ) const;
  __inline Double xGetIEPRate          (                                                                   ) const;


  // dequantization
  Void xDeQuant(       TComTU       &rTu,
                 const TCoeff      * pSrc,
                       TCoeff      * pDes,
                 const ComponentID   compID,
                 const QpParam      &cQP );

  // inverse transform
  Void xIT    ( const Int channelBitDepth, Bool useDST, TCoeff* plCoef, Pel* pResidual, UInt uiStride, Int iWidth, Int iHeight, const Int maxLog2TrDynamicRange 
#if COM16_C806_EMT
    , UChar ucMode
    , UChar ucTrIdx
#endif
#if JVET_C0024_ITSKIP
    , UInt uiSkipWidth
    , UInt uiSkipHeight
#endif
#if VCEG_AZ08_KLT_COMMON
    , Bool useKLT
#endif
    );

  // inverse skipping transform
  Void xITransformSkip ( TCoeff* plCoef, Pel* pResidual, UInt uiStride, TComTU &rTu, const ComponentID component );

public:
  static Void crossComponentPrediction(      TComTU      &rTu,
                                       const ComponentID  compID,
                                       const Pel         *piResiL,
                                       const Pel         *piResiC,
                                             Pel         *piResiT,
                                       const Int          width,
                                       const Int          height,
                                       const Int          strideL,
                                       const Int          strideC,
                                       const Int          strideT,
                                       const Bool         reverse);

};// END CLASS DEFINITION TComTrQuant

//! \}

#endif // __TCOMTRQUANT__
