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

/** \file     TComRdCost.cpp
    \brief    RD cost computation class
*/

#include <math.h>
#include <assert.h>
#include <limits>
#include "TComRom.h"
#include "TComRdCost.h"

#if COM16_C806_SIMD_OPT
#include <emmintrin.h>  
#include <xmmintrin.h>
#endif

//! \ingroup TLibCommon
//! \{
#if WCG_LUMA_DQP_CM_SCALE
Double  TComRdCost::m_lumaLevelToWeightPLUT[LUMA_LEVEL_TO_DQP_LUT_MAXSIZE];
#endif

TComRdCost::TComRdCost()
{
  init();
}

TComRdCost::~TComRdCost()
{
}

// Calculate RD functions
Double TComRdCost::calcRdCost( UInt uiBits, Distortion uiDistortion, Bool bFlag, DFunc eDFunc )
{
  Double dRdCost = 0.0;
  Double dLambda = 0.0;

  switch ( eDFunc )
  {
    case DF_SSE:
      assert(0);
      break;
    case DF_SAD:
#if RExt__HIGH_BIT_DEPTH_SUPPORT
      dLambda = m_dLambdaMotionSAD[0]; // 0 is valid, because for lossless blocks, the cost equation is modified to compensate.
#else
      dLambda = (Double)m_uiLambdaMotionSAD[0]; // 0 is valid, because for lossless blocks, the cost equation is modified to compensate.
#endif
      break;
#if !WCG_LUMA_DQP_CM_SCALE
    case DF_DEFAULT:
      dLambda =         m_dLambda;
      break;
#else
    case DF_DEFAULT_ORI:
      dLambda =         m_dLambda;
      break;
    case DF_DEFAULT:
      dLambda =         m_dLambda_unadjusted;  // use lamda without the luma QP adjustment; when luma adaptive QP is not enabled, m_dLambda_unadjusted is same as m_dLambda
      break;
#endif
    case DF_SSE_FRAME:
      dLambda =         m_dFrameLambda;
      break;
    default:
      assert (0);
      break;
  }

  if (bFlag) //NOTE: this "bFlag" is never true
  {
    // Intra8x8, Intra4x4 Block only...
    if (m_costMode != COST_STANDARD_LOSSY)
    {
      dRdCost = (Double(uiDistortion) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
    }
    else
    {
      dRdCost = (((Double)uiDistortion) + ((Double)uiBits * dLambda));
    }
  }
  else
  {
    if (eDFunc == DF_SAD)
    {
      if (m_costMode != COST_STANDARD_LOSSY)
      {
        dRdCost = ((Double(uiDistortion) * 65536) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
      }
      else
      {
        dRdCost = floor(Double(uiDistortion) + (floor((Double(uiBits) * dLambda) + 0.5) / 65536.0));
      }
    }
    else
    {
      if (m_costMode != COST_STANDARD_LOSSY)
      {
        dRdCost = (Double(uiDistortion) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
      }
      else
      {
        dRdCost = floor(Double(uiDistortion) + (Double(uiBits) * dLambda) + 0.5);
      }
    }
  }

  return dRdCost;
}

Double TComRdCost::calcRdCost64( UInt64 uiBits, UInt64 uiDistortion, Bool bFlag, DFunc eDFunc )
{
  Double dRdCost = 0.0;
  Double dLambda = 0.0;

  switch ( eDFunc )
  {
    case DF_SSE:
      assert(0);
      break;
    case DF_SAD:
#if RExt__HIGH_BIT_DEPTH_SUPPORT
      dLambda = m_dLambdaMotionSAD[0]; // 0 is valid, because for lossless blocks, the cost equation is modified to compensate.
#else
      dLambda = (Double)m_uiLambdaMotionSAD[0]; // 0 is valid, because for lossless blocks, the cost equation is modified to compensate.
#endif
      break;
    case DF_DEFAULT:
      dLambda =         m_dLambda;
      break;
    case DF_SSE_FRAME:
      dLambda =         m_dFrameLambda;
      break;
    default:
      assert (0);
      break;
  }

  if (bFlag) //NOTE: this "bFlag" is never true
  {
    // Intra8x8, Intra4x4 Block only...
    if (m_costMode != COST_STANDARD_LOSSY)
    {
      dRdCost = (Double(uiDistortion) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
    }
    else
    {
      dRdCost = (((Double)(Int64)uiDistortion) + ((Double)(Int64)uiBits * dLambda));
    }
  }
  else
  {
    if (eDFunc == DF_SAD)
    {
      if (m_costMode != COST_STANDARD_LOSSY)
      {
        dRdCost = ((Double(uiDistortion) * 65536) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
      }
      else
      {
        dRdCost = floor(Double(uiDistortion) + (floor((Double(uiBits) * dLambda) + 0.5) / 65536.0));
      }
    }
    else
    {
      if (m_costMode != COST_STANDARD_LOSSY)
      {
        dRdCost = (Double(uiDistortion) / dLambda) + Double(uiBits); // all lossless costs would have uiDistortion=0, and therefore this cost function can be used.
      }
      else
      {
        dRdCost = floor(Double(uiDistortion) + (Double(uiBits) * dLambda) + 0.5);
      }
    }
  }

  return dRdCost;
}

Void TComRdCost::setLambda( Double dLambda, const BitDepths &bitDepths )
{
  m_dLambda           = dLambda;
  m_sqrtLambda        = sqrt(m_dLambda);
#if RExt__HIGH_BIT_DEPTH_SUPPORT
  m_dLambdaMotionSAD[0] = 65536.0 * m_sqrtLambda;
  m_dLambdaMotionSSE[0] = 65536.0 * m_dLambda;
#if FULL_NBIT
  dLambda = 0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12) / 3.0));
#else
  dLambda = 0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12 - 6 * (bitDepths.recon[CHANNEL_TYPE_LUMA] - 8)) / 3.0));
#endif
  m_dLambdaMotionSAD[1] = 65536.0 * sqrt(dLambda);
  m_dLambdaMotionSSE[1] = 65536.0 * dLambda;
#else
  m_uiLambdaMotionSAD[0] = (UInt)floor(65536.0 * m_sqrtLambda);
  m_uiLambdaMotionSSE[0] = (UInt)floor(65536.0 * m_dLambda   );
#if FULL_NBIT
  dLambda = 0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12) / 3.0));
#else
  dLambda = 0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12 - 6 * (bitDepths.recon[CHANNEL_TYPE_LUMA] - 8)) / 3.0));
#endif
  m_uiLambdaMotionSAD[1] = (UInt)floor(65536.0 * sqrt(dLambda));
  m_uiLambdaMotionSSE[1] = (UInt)floor(65536.0 * dLambda   );
#endif
}


// Initalize Function Pointer by [eDFunc]
Void TComRdCost::init()
{
  m_afpDistortFunc[DF_DEFAULT] = NULL;                  // for DF_DEFAULT

  m_afpDistortFunc[DF_SSE    ] = TComRdCost::xGetSSE;
  m_afpDistortFunc[DF_SSE4   ] = TComRdCost::xGetSSE4;
  m_afpDistortFunc[DF_SSE8   ] = TComRdCost::xGetSSE8;
  m_afpDistortFunc[DF_SSE16  ] = TComRdCost::xGetSSE16;
  m_afpDistortFunc[DF_SSE32  ] = TComRdCost::xGetSSE32;
  m_afpDistortFunc[DF_SSE64  ] = TComRdCost::xGetSSE64;
  m_afpDistortFunc[DF_SSE16N ] = TComRdCost::xGetSSE16N;

  m_afpDistortFunc[DF_SAD    ] = TComRdCost::xGetSAD;
  m_afpDistortFunc[DF_SAD4   ] = TComRdCost::xGetSAD4;
  m_afpDistortFunc[DF_SAD8   ] = TComRdCost::xGetSAD8;
  m_afpDistortFunc[DF_SAD16  ] = TComRdCost::xGetSAD16;
  m_afpDistortFunc[DF_SAD32  ] = TComRdCost::xGetSAD32;
  m_afpDistortFunc[DF_SAD64  ] = TComRdCost::xGetSAD64;
  m_afpDistortFunc[DF_SAD16N ] = TComRdCost::xGetSAD16N;

  m_afpDistortFunc[DF_SADS   ] = TComRdCost::xGetSAD;
  m_afpDistortFunc[DF_SADS4  ] = TComRdCost::xGetSAD4;
  m_afpDistortFunc[DF_SADS8  ] = TComRdCost::xGetSAD8;
  m_afpDistortFunc[DF_SADS16 ] = TComRdCost::xGetSAD16;
  m_afpDistortFunc[DF_SADS32 ] = TComRdCost::xGetSAD32;
  m_afpDistortFunc[DF_SADS64 ] = TComRdCost::xGetSAD64;
  m_afpDistortFunc[DF_SADS16N] = TComRdCost::xGetSAD16N;

  m_afpDistortFunc[DF_SAD12  ] = TComRdCost::xGetSAD12;
  m_afpDistortFunc[DF_SAD24  ] = TComRdCost::xGetSAD24;
  m_afpDistortFunc[DF_SAD48  ] = TComRdCost::xGetSAD48;

  m_afpDistortFunc[DF_SADS12 ] = TComRdCost::xGetSAD12;
  m_afpDistortFunc[DF_SADS24 ] = TComRdCost::xGetSAD24;
  m_afpDistortFunc[DF_SADS48 ] = TComRdCost::xGetSAD48;

  m_afpDistortFunc[DF_HADS   ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS4  ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS8  ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS16 ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS32 ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS64 ] = TComRdCost::xGetHADs;
  m_afpDistortFunc[DF_HADS16N] = TComRdCost::xGetHADs;

#if WCG_LUMA_DQP_CM_SCALE
  m_afpDistortFunc[DF_SSE_WTD    ] = TComRdCost::xGetSSE_WTD;
  m_afpDistortFunc[DF_SSE4_WTD   ] = TComRdCost::xGetSSE4_WTD;
  m_afpDistortFunc[DF_SSE8_WTD   ] = TComRdCost::xGetSSE8_WTD;
  m_afpDistortFunc[DF_SSE16_WTD  ] = TComRdCost::xGetSSE16_WTD;
  m_afpDistortFunc[DF_SSE32_WTD  ] = TComRdCost::xGetSSE32_WTD;
  m_afpDistortFunc[DF_SSE64_WTD  ] = TComRdCost::xGetSSE64_WTD;
  m_afpDistortFunc[DF_SSE16N_WTD ] = TComRdCost::xGetSSE16N_WTD;
#endif
  m_costMode                   = COST_STANDARD_LOSSY;

#if RExt__HIGH_BIT_DEPTH_SUPPORT
  m_dCost                      = 0;
#else
  m_uiCost                     = 0;
#endif
  m_iCostScale                 = 0;
}

// Static member function
UInt TComRdCost::xGetExpGolombNumberOfBits( Int iVal 
#if VCEG_AZ07_IMV
  , Int iMvFlag 
#endif
  )
{
#if VCEG_AZ07_IMV
  if( iMvFlag )
  {
    iVal >>= 2;
#if JVET_E0076_MULTI_PEL_MVD
    if (iMvFlag == 2)
    {
      iVal >>= MULTI_PEL_MVD_BITS;
    }
#endif
  }
#endif

  assert(iVal != std::numeric_limits<Int>::min());
  UInt uiLength = 1;
  UInt uiTemp   = ( iVal <= 0) ? (UInt(-iVal)<<1)+1: UInt(iVal<<1);

  while ( 1 != uiTemp )
  {
    uiTemp >>= 1;
    uiLength += 2;
  }

  return uiLength;
}

Void TComRdCost::setDistParam( UInt uiBlkWidth, UInt uiBlkHeight, DFunc eDFunc, DistParam& rcDistParam )
{
  // set Block Width / Height
  rcDistParam.iCols    = uiBlkWidth;
  rcDistParam.iRows    = uiBlkHeight;
#if COM16_C806_LARGE_CTU
  if( rcDistParam.iCols > 64 )
    rcDistParam.DistFunc = m_afpDistortFunc[eDFunc];
  else
#endif
#if JVET_C0024_QTBT
  rcDistParam.DistFunc = m_afpDistortFunc[eDFunc + g_aucConvertToBit[ rcDistParam.iCols ] + MIN_CU_LOG2 - 1 ];
#else
  rcDistParam.DistFunc = m_afpDistortFunc[eDFunc + g_aucConvertToBit[ rcDistParam.iCols ] + 1 ];
#endif

  // initialize
  rcDistParam.iSubShift  = 0;
}

// Setting the Distortion Parameter for Inter (ME)
Void TComRdCost::setDistParam( TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, DistParam& rcDistParam )
{
  // set Original & Curr Pointer / Stride
  rcDistParam.pOrg = pcPatternKey->getROIY();
  rcDistParam.pCur = piRefY;

  rcDistParam.iStrideOrg = pcPatternKey->getPatternLStride();
  rcDistParam.iStrideCur = iRefStride;

  // set Block Width / Height
  rcDistParam.iCols    = pcPatternKey->getROIYWidth();
  rcDistParam.iRows    = pcPatternKey->getROIYHeight();
#if COM16_C806_LARGE_CTU
  if( rcDistParam.iCols > 64 )
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD];
  else
#endif
#if JVET_C0024_QTBT
  rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD + g_aucConvertToBit[ rcDistParam.iCols ] + MIN_CU_LOG2 - 1 ];
#else
  rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD + g_aucConvertToBit[ rcDistParam.iCols ] + 1 ];
#endif

  if (rcDistParam.iCols == 12)
  {
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD12];
  }
  else if (rcDistParam.iCols == 24)
  {
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD24];
  }
  else if (rcDistParam.iCols == 48)
  {
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SAD48];
  }

  // initialize
  rcDistParam.iSubShift  = 0;
}

// Setting the Distortion Parameter for Inter (subpel ME with step)
Void TComRdCost::setDistParam( TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, Int iStep, DistParam& rcDistParam, Bool bHADME )
{
  // set Original & Curr Pointer / Stride
  rcDistParam.pOrg = pcPatternKey->getROIY();
  rcDistParam.pCur = piRefY;

  rcDistParam.iStrideOrg = pcPatternKey->getPatternLStride();
  rcDistParam.iStrideCur = iRefStride * iStep;

  // set Step for interpolated buffer
  rcDistParam.iStep = iStep;

  // set Block Width / Height
  rcDistParam.iCols    = pcPatternKey->getROIYWidth();
  rcDistParam.iRows    = pcPatternKey->getROIYHeight();

  // set distortion function
  if ( !bHADME )
  {
#if COM16_C806_LARGE_CTU
    if( rcDistParam.iCols > 64 )
      rcDistParam.DistFunc =  m_afpDistortFunc[DF_SADS];
    else
#endif
#if JVET_C0024_QTBT
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SADS + g_aucConvertToBit[ rcDistParam.iCols ] + MIN_CU_LOG2 - 1 ];
#else
    rcDistParam.DistFunc = m_afpDistortFunc[DF_SADS + g_aucConvertToBit[ rcDistParam.iCols ] + 1 ];
#endif
    if (rcDistParam.iCols == 12)
    {
      rcDistParam.DistFunc = m_afpDistortFunc[DF_SADS12];
    }
    else if (rcDistParam.iCols == 24)
    {
      rcDistParam.DistFunc = m_afpDistortFunc[DF_SADS24];
    }
    else if (rcDistParam.iCols == 48)
    {
      rcDistParam.DistFunc = m_afpDistortFunc[DF_SADS48];
    }
  }
  else
  {
#if COM16_C806_LARGE_CTU
    if( rcDistParam.iCols > 64 )
      rcDistParam.DistFunc =  m_afpDistortFunc[DF_HADS];
    else
#endif
    rcDistParam.DistFunc = m_afpDistortFunc[DF_HADS + g_aucConvertToBit[ rcDistParam.iCols ] + 1 ];
  }

  // initialize
  rcDistParam.iSubShift  = 0;
}

Void TComRdCost::setDistParam( DistParam& rcDP, Int bitDepth, Pel* p1, Int iStride1, Pel* p2, Int iStride2, Int iWidth, Int iHeight, Bool bHadamard )
{
  rcDP.pOrg       = p1;
  rcDP.pCur       = p2;
  rcDP.iStrideOrg = iStride1;
  rcDP.iStrideCur = iStride2;
  rcDP.iCols      = iWidth;
  rcDP.iRows      = iHeight;
  rcDP.iStep      = 1;
  rcDP.iSubShift  = 0;
  rcDP.bitDepth   = bitDepth;
#if COM16_C806_LARGE_CTU
  if( rcDP.iCols > 64 )
    rcDP.DistFunc   = m_afpDistortFunc[ ( bHadamard ? DF_HADS : DF_SADS ) ];
  else
#endif
#if JVET_C0024_QTBT
  rcDP.DistFunc   = m_afpDistortFunc[ ( bHadamard ? DF_HADS : DF_SADS ) + g_aucConvertToBit[ iWidth ] + MIN_CU_LOG2 - 1 ];
#else
  rcDP.DistFunc   = m_afpDistortFunc[ ( bHadamard ? DF_HADS : DF_SADS ) + g_aucConvertToBit[ iWidth ] + 1 ];
#endif
}

Distortion TComRdCost::calcHAD( Int bitDepth, Pel* pi0, Int iStride0, Pel* pi1, Int iStride1, Int iWidth, Int iHeight )
{
  Distortion uiSum = 0;
  Int x, y;

  if ( ( (iWidth % 8) == 0 ) && ( (iHeight % 8) == 0 ) )
  {
    for ( y=0; y<iHeight; y+= 8 )
    {
      for ( x=0; x<iWidth; x+= 8 )
      {
        uiSum += xCalcHADs8x8( &pi0[x], &pi1[x], iStride0, iStride1, 1 
#if COM16_C806_SIMD_OPT
          , bitDepth
#endif
          );
      }
      pi0 += iStride0*8;
      pi1 += iStride1*8;
    }
  }
  else
  {
    assert ( ( (iWidth % 4) == 0 ) && ( (iHeight % 4) == 0 ) );

    for ( y=0; y<iHeight; y+= 4 )
    {
      for ( x=0; x<iWidth; x+= 4 )
      {
        uiSum += xCalcHADs4x4( &pi0[x], &pi1[x], iStride0, iStride1, 1, bitDepth);
      }
      pi0 += iStride0*4;
      pi1 += iStride1*4;
    }
  }

  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(bitDepth-8) );
}

Distortion TComRdCost::getDistPart( Int bitDepth, Pel* piCur, Int iCurStride,  Pel* piOrg, Int iOrgStride, UInt uiBlkWidth, UInt uiBlkHeight, const ComponentID compID, DFunc eDFunc 
#if WCG_LUMA_DQP_CM_SCALE
  , Pel* piOrgLuma, 
  Int iOrgStrideLuma
#endif
#if VCEG_AZ06_IC
  , Bool bMRFlag
#endif
  )
{
  DistParam cDtParam;
  setDistParam( uiBlkWidth, uiBlkHeight, eDFunc, cDtParam );
  cDtParam.pOrg       = piOrg;
  cDtParam.pCur       = piCur;
  cDtParam.iStrideOrg = iOrgStride;
  cDtParam.iStrideCur = iCurStride;
  cDtParam.iStep      = 1;

  cDtParam.bApplyWeight = false;
#if WCG_LUMA_DQP_CM_SCALE
  cDtParam.compIdx      = compID;
#else
  cDtParam.compIdx      = MAX_NUM_COMPONENT; // just for assert: to be sure it was set before use
#endif
  cDtParam.bitDepth     = bitDepth;
#if VCEG_AZ06_IC
  cDtParam.bMRFlag      = bMRFlag;
#endif
  if (isChroma(compID))
  {
#if WCG_LUMA_DQP_CM_SCALE
    cDtParam.pOrgLuma = piOrgLuma;
    cDtParam.iStrideOrgLuma = iOrgStrideLuma;
#endif
    return ((Distortion) (m_distortionWeight[compID] * cDtParam.DistFunc( &cDtParam )));
  }
  else
  {
#if WCG_LUMA_DQP_CM_SCALE
    cDtParam.pOrgLuma = piOrg;
    cDtParam.iStrideOrgLuma = iOrgStride;
#endif
    return cDtParam.DistFunc( &cDtParam );
  }
}

// ====================================================================================================================
// Distortion functions
// ====================================================================================================================

// --------------------------------------------------------------------------------------------------------------------
// SAD
// --------------------------------------------------------------------------------------------------------------------

#if COM16_C806_SIMD_OPT
inline Int simdSADLine4n16b( const Pel * piOrg , const Pel * piCur , Int nWidth )
{
  // internal bit-depth must be 12-bit or lower
  assert( !( nWidth & 0x03 ) );
  __m128i org , cur , abs , sum;
  sum = _mm_setzero_si128();
  for( Int n = 0 ; n < nWidth ; n += 4 )
  {
    org = _mm_loadl_epi64( ( __m128i* )( piOrg + n ) );
    cur = _mm_loadl_epi64( ( __m128i* )( piCur + n ) );
    abs = _mm_subs_epi16( _mm_max_epi16( org , cur )  , _mm_min_epi16( org , cur ) );
    sum = _mm_adds_epu16( abs , sum );
  }
  __m128i zero =  _mm_setzero_si128();
  sum = _mm_unpacklo_epi16( sum , zero );
  sum = _mm_add_epi32( sum , _mm_shuffle_epi32( sum , _MM_SHUFFLE( 2 , 3 , 0 , 1 ) ) );
  sum = _mm_add_epi32( sum , _mm_shuffle_epi32( sum , _MM_SHUFFLE( 1 , 0 , 3 , 2 ) ) );
  return( _mm_cvtsi128_si32( sum ) );
}

inline Int simdSADLine8n16b( const Pel * piOrg , const Pel * piCur , Int nWidth )
{
  // internal bit-depth must be 12-bit or lower
  assert( !( nWidth & 0x07 ) );
  __m128i org , cur , abs , sum;
  sum = _mm_setzero_si128();
  for( Int n = 0 ; n < nWidth ; n += 8 )
  {
    org = _mm_loadu_si128( ( __m128i* )( piOrg + n ) );
    cur = _mm_loadu_si128( ( __m128i* )( piCur + n ) );
    abs = _mm_subs_epi16( _mm_max_epi16( org , cur )  , _mm_min_epi16( org , cur ) );
    sum = _mm_adds_epu16( abs , sum );
  }
  __m128i zero =  _mm_setzero_si128();
  __m128i hi = _mm_unpackhi_epi16( sum , zero );
  __m128i lo = _mm_unpacklo_epi16( sum , zero );
  sum = _mm_add_epi32( lo , hi );
  sum = _mm_add_epi32( sum , _mm_shuffle_epi32( sum , _MM_SHUFFLE( 2 , 3 , 0 , 1 ) ) );
  sum = _mm_add_epi32( sum , _mm_shuffle_epi32( sum , _MM_SHUFFLE( 1 , 0 , 3 , 2 ) ) );
  return( _mm_cvtsi128_si32( sum ) );
}

inline Void simd8x8Transpose32b( __m128i * pBuffer )
{
  __m128 tmp[16];
  for( Int n = 0 ; n < 16 ; n++ )
  {
    tmp[n] = _mm_castsi128_ps( pBuffer[n] );
  }
  _MM_TRANSPOSE4_PS( tmp[0] , tmp[2] , tmp[4] , tmp[6] );
  _MM_TRANSPOSE4_PS( tmp[1] , tmp[3] , tmp[5] , tmp[7] );
  _MM_TRANSPOSE4_PS( tmp[8] , tmp[10] , tmp[12] , tmp[14] );
  _MM_TRANSPOSE4_PS( tmp[9] , tmp[11] , tmp[13] , tmp[15] );
  for( Int n = 0 ; n < 8 ; n += 2 )
  {
    pBuffer[n] = _mm_castps_si128( tmp[n] );
    pBuffer[n+1]  = _mm_castps_si128( tmp[n+8] );
    pBuffer[n+8] = _mm_castps_si128( tmp[n+1] );
    pBuffer[n+9]  = _mm_castps_si128( tmp[n+9] );
  }
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
Void simd8x8HAD1D32b( __m128i * pInput , __m128i * pOutput )
{
  __m128i m1[8][2] , m2[8][2];

  m2[0][0] = _mm_add_epi32( pInput[0] ,pInput[8 ] );  m2[0][1] = _mm_add_epi32( pInput[1] ,pInput[9 ] );
  m2[1][0] = _mm_add_epi32( pInput[2] ,pInput[10] );  m2[1][1] = _mm_add_epi32( pInput[3] ,pInput[11] );
  m2[2][0] = _mm_add_epi32( pInput[4] ,pInput[12] );  m2[2][1] = _mm_add_epi32( pInput[5] ,pInput[13] );
  m2[3][0] = _mm_add_epi32( pInput[6] ,pInput[14] );  m2[3][1] = _mm_add_epi32( pInput[7] ,pInput[15] );
  m2[4][0] = _mm_sub_epi32( pInput[0] ,pInput[8 ] );  m2[4][1] = _mm_sub_epi32( pInput[1] ,pInput[9 ] );
  m2[5][0] = _mm_sub_epi32( pInput[2] ,pInput[10] );  m2[5][1] = _mm_sub_epi32( pInput[3] ,pInput[11] );
  m2[6][0] = _mm_sub_epi32( pInput[4] ,pInput[12] );  m2[6][1] = _mm_sub_epi32( pInput[5] ,pInput[13] );
  m2[7][0] = _mm_sub_epi32( pInput[6] ,pInput[14] );  m2[7][1] = _mm_sub_epi32( pInput[7] ,pInput[15] );

  m1[0][0] = _mm_add_epi32( m2[0][0] , m2[2][0] );  m1[0][1] = _mm_add_epi32( m2[0][1] , m2[2][1] );
  m1[1][0] = _mm_add_epi32( m2[1][0] , m2[3][0] );  m1[1][1] = _mm_add_epi32( m2[1][1] , m2[3][1] );
  m1[2][0] = _mm_sub_epi32( m2[0][0] , m2[2][0] );  m1[2][1] = _mm_sub_epi32( m2[0][1] , m2[2][1] );
  m1[3][0] = _mm_sub_epi32( m2[1][0] , m2[3][0] );  m1[3][1] = _mm_sub_epi32( m2[1][1] , m2[3][1] );
  m1[4][0] = _mm_add_epi32( m2[4][0] , m2[6][0] );  m1[4][1] = _mm_add_epi32( m2[4][1] , m2[6][1] );
  m1[5][0] = _mm_add_epi32( m2[5][0] , m2[7][0] );  m1[5][1] = _mm_add_epi32( m2[5][1] , m2[7][1] );
  m1[6][0] = _mm_sub_epi32( m2[4][0] , m2[6][0] );  m1[6][1] = _mm_sub_epi32( m2[4][1] , m2[6][1] );
  m1[7][0] = _mm_sub_epi32( m2[5][0] , m2[7][0] );  m1[7][1] = _mm_sub_epi32( m2[5][1] , m2[7][1] );

  pInput[0 ] = _mm_add_epi32( m1[0][0] , m1[1][0] );  pInput[1 ] = _mm_add_epi32( m1[0][1] , m1[1][1] );
  pInput[2 ] = _mm_sub_epi32( m1[0][0] , m1[1][0] );  pInput[3 ] = _mm_sub_epi32( m1[0][1] , m1[1][1] );
  pInput[4 ] = _mm_add_epi32( m1[2][0] , m1[3][0] );  pInput[5 ] = _mm_add_epi32( m1[2][1] , m1[3][1] );
  pInput[6 ] = _mm_sub_epi32( m1[2][0] , m1[3][0] );  pInput[7 ] = _mm_sub_epi32( m1[2][1] , m1[3][1] );
  pInput[8 ] = _mm_add_epi32( m1[4][0] , m1[5][0] );  pInput[9 ] = _mm_add_epi32( m1[4][1] , m1[5][1] );
  pInput[10] = _mm_sub_epi32( m1[4][0] , m1[5][0] );  pInput[11] = _mm_sub_epi32( m1[4][1] , m1[5][1] );
  pInput[12] = _mm_add_epi32( m1[6][0] , m1[7][0] );  pInput[13] = _mm_add_epi32( m1[6][1] , m1[7][1] );
  pInput[14] = _mm_sub_epi32( m1[6][0] , m1[7][0] );  pInput[15] = _mm_sub_epi32( m1[6][1] , m1[7][1] );
}

inline __m128i simdAbs32b( __m128i m )
{
  const __m128i zero = _mm_setzero_si128();
  __m128i tmp = _mm_sub_epi32( zero , m );
  __m128i mask = _mm_cmpgt_epi32( m , tmp );
  return( _mm_or_si128( _mm_and_si128( mask , m ) , _mm_andnot_si128( mask , tmp ) ) );
}

UInt simdHADs8x8( Pel * piOrg, Pel * piCur, Int iStrideOrg, Int iStrideCur )
{
  __m128i mmDiff[8][2];
  __m128i mmZero = _mm_setzero_si128();
  for( Int n = 0 ; n < 8 ; n++ , piOrg += iStrideOrg , piCur += iStrideCur )
  {
    __m128i diff = _mm_sub_epi16( _mm_loadu_si128( ( __m128i* )piOrg ) , _mm_loadu_si128( ( __m128i* )piCur ) );
    // sign extension
    __m128i mask = _mm_cmplt_epi16( diff , mmZero );
    mmDiff[n][0] = _mm_unpacklo_epi16( diff , mask );
    mmDiff[n][1] = _mm_unpackhi_epi16( diff , mask );
  }

  // transpose
  simd8x8Transpose32b( &mmDiff[0][0] );

  // horizontal
  simd8x8HAD1D32b( &mmDiff[0][0] , &mmDiff[0][0] );

  // transpose
  simd8x8Transpose32b( &mmDiff[0][0] );

  // vertical
  simd8x8HAD1D32b( &mmDiff[0][0] , &mmDiff[0][0] );

  __m128i mmSum = _mm_setzero_si128();
  for( Int n = 0 ; n < 8 ; n++ )
  {
    mmSum = _mm_add_epi32( mmSum , simdAbs32b( mmDiff[n][0] ) );
    mmSum = _mm_add_epi32( mmSum , simdAbs32b( mmDiff[n][1] ) );
  }
  mmSum = _mm_add_epi32( mmSum , _mm_shuffle_epi32( mmSum , _MM_SHUFFLE( 2 , 3 , 0 , 1 ) ) );
  mmSum = _mm_add_epi32( mmSum , _mm_shuffle_epi32( mmSum , _MM_SHUFFLE( 1 , 0 , 3 , 2 ) ) );

  UInt sad = _mm_cvtsi128_si32( mmSum );
  sad = ( sad + 2 ) >> 2;

  return( sad );
}
#endif

Distortion TComRdCost::xGetSAD( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideCur = pcDtParam->iStrideCur;
  Int  iStrideOrg = pcDtParam->iStrideOrg;

  Distortion uiSum = 0;

#if ITH_C3_SADCOMPUTER!=SIMD_FXN
#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    if( ( iCols & 0x07 ) == 0 )
    {
      for( ; iRows != 0; iRows-- )
      {
        uiSum += simdSADLine8n16b( piOrg , piCur , iCols );
        piOrg += iStrideOrg;
        piCur += iStrideCur;
      }
    }
    else
    {
      for( ; iRows != 0; iRows-- )
      {
        uiSum += simdSADLine4n16b( piOrg , piCur , iCols );
        piOrg += iStrideOrg;
        piCur += iStrideCur;
      }
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n++ )
    {
      uiSum += abs( piOrg[n] - piCur[n] );
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
#elif ITH_C3_SADCOMPUTER==SIMD_FXN
  int iSubShift = pcDtParam->iSubShift;
  Int iBitDepth = pcDtParam->bitDepth;
  assert(iBitDepth <= 14);

  /*widths and height has to be more than 4*/
  assert(iCols >= 4);
  assert(iRows >= 4);

  /*widths & heights has to be even*/
  assert(!(iCols & 1));
  assert(!(iRows & 1));

  /*iSubShift has to be either  or 1 in this function*/
  assert((iSubShift == 0) || (iSubShift == 1));

  __m128i src_8x16b;
  __m128i src_8x16b_1;
  __m128i src_8x16b_2;
  __m128i src_8x16b_3;

  __m128i pred_8x16b;
  __m128i pred_8x16b_1;
  __m128i pred_8x16b_2;
  __m128i pred_8x16b_3;

  __m128i temp;
  __m128i temp_1;
  __m128i temp_2;
  __m128i temp_3;
  __m128i temp_4;
  __m128i temp_5;

  __m128i temp_dummy;
  __m128i result, result_1;

  Int ht = iRows;

  Int  iSubStep = (1 << iSubShift);

  Short *pu2_inp;
  Short *pu2_ref;
  pu2_inp = pcDtParam->pOrg;
  pu2_ref = pcDtParam->pCur;

  Int  i, j, iRemWd;
  Int sad = 0;
  Int  iSrcStride = iStrideOrg;
  Int  iRefStride = iStrideCur;


  iRemWd = iCols;
  temp_dummy = _mm_setzero_si128();
  result = _mm_setzero_si128();
  result_1 = _mm_setzero_si128();

  if (1 == iSubStep)
  {
	  if (iRemWd >= 8)
	  {
		  for (i = 0; i < ht / 4; i++)
		  {
			  Int iCount = 0;

			  for (j = iRemWd; j > 7; j -= 8)
			  {
				  src_8x16b = _mm_loadu_si128((__m128i *) (pu2_inp + iCount));
				  src_8x16b_1 = _mm_loadu_si128((__m128i *) (pu2_inp + iCount + iSrcStride));
				  src_8x16b_2 = _mm_loadu_si128((__m128i *) (pu2_inp + iCount + (iSrcStride * 2)));
				  src_8x16b_3 = _mm_loadu_si128((__m128i *) (pu2_inp + iCount + (iSrcStride * 3)));

				  pred_8x16b = _mm_loadu_si128((__m128i *) (pu2_ref + iCount));
				  pred_8x16b_1 = _mm_loadu_si128((__m128i *) (pu2_ref + iCount + iRefStride));
				  pred_8x16b_2 = _mm_loadu_si128((__m128i *) (pu2_ref + iCount + (iRefStride * 2)));
				  pred_8x16b_3 = _mm_loadu_si128((__m128i *) (pu2_ref + iCount + (iRefStride * 3)));

				  temp = _mm_sub_epi16(src_8x16b, pred_8x16b);
				  temp_1 = _mm_sub_epi16(src_8x16b_1, pred_8x16b_1);
				  temp_2 = _mm_sub_epi16(src_8x16b_2, pred_8x16b_2);
				  temp_3 = _mm_sub_epi16(src_8x16b_3, pred_8x16b_3);

				  temp = _mm_abs_epi16(temp);
				  temp_1 = _mm_abs_epi16(temp_1);
				  temp_2 = _mm_abs_epi16(temp_2);
				  temp_3 = _mm_abs_epi16(temp_3);

				  temp = _mm_add_epi16(temp, temp_1);
				  temp_2 = _mm_add_epi16(temp_2, temp_3);

				  temp_1 = _mm_unpackhi_epi16(temp, temp_dummy);
				  temp_3 = _mm_unpacklo_epi16(temp, temp_dummy);
				  temp_4 = _mm_unpackhi_epi16(temp_2, temp_dummy);
				  temp_5 = _mm_unpacklo_epi16(temp_2, temp_dummy);

				  temp = _mm_add_epi32(temp_1, temp_3);
				  temp_2 = _mm_add_epi32(temp_4, temp_5);

				  result = _mm_add_epi32(result, temp);
				  result_1 = _mm_add_epi32(result_1, temp_2);

				  iCount += 8;
			  }

			  pu2_inp += (4 * iSrcStride);
			  pu2_ref += (4 * iRefStride);
		  }
		  result = _mm_add_epi32(result, result_1);
	  }

	  iRemWd &= 0x7;
	  pu2_inp = pcDtParam->pOrg + ((iCols / 8) * 8);
	  pu2_ref = pcDtParam->pCur + ((iCols / 8) * 8);

	  if (iRemWd > 3)
	  {
		  for (i = 0; i < (ht / 4); i++)
		  {
			  src_8x16b = _mm_loadl_epi64((__m128i *) (pu2_inp));
			  src_8x16b_1 = _mm_loadl_epi64((__m128i *) (pu2_inp + iSrcStride));
			  src_8x16b_2 = _mm_loadl_epi64((__m128i *) (pu2_inp + (iSrcStride * 2)));
			  src_8x16b_3 = _mm_loadl_epi64((__m128i *) (pu2_inp + (iSrcStride * 3)));

			  pred_8x16b = _mm_loadl_epi64((__m128i *) (pu2_ref));
			  pred_8x16b_1 = _mm_loadl_epi64((__m128i *) (pu2_ref + iRefStride));
			  pred_8x16b_2 = _mm_loadl_epi64((__m128i *) (pu2_ref + (iRefStride * 2)));
			  pred_8x16b_3 = _mm_loadl_epi64((__m128i *) (pu2_ref + (iRefStride * 3)));

			  src_8x16b = _mm_unpacklo_epi16(src_8x16b, src_8x16b_1);
			  pred_8x16b = _mm_unpacklo_epi16(pred_8x16b, pred_8x16b_1);
			  src_8x16b_2 = _mm_unpacklo_epi16(src_8x16b_2, src_8x16b_3);
			  pred_8x16b_2 = _mm_unpacklo_epi16(pred_8x16b_2, pred_8x16b_3);

			  temp = _mm_sub_epi16(src_8x16b, pred_8x16b);
			  temp_1 = _mm_sub_epi16(src_8x16b_2, pred_8x16b_2);

			  temp = _mm_abs_epi16(temp);
			  temp_1 = _mm_abs_epi16(temp_1);

			  temp = _mm_adds_epu16(temp, temp_1);

			  temp_1 = _mm_unpackhi_epi16(temp, temp_dummy);
			  temp_2 = _mm_unpacklo_epi16(temp, temp_dummy);

			  temp_1 = _mm_add_epi32(temp_2, temp_1);

			  result = _mm_add_epi32(result, temp_1);

			  pu2_inp += 4 * (iSrcStride);
			  pu2_ref += 4 * (iRefStride);

		  }
	  }

	  {
		  Int *val = (Int*)&result;
		  sad = val[0] + val[1] + val[2] + val[3];
	  }
	  sad <<= iSubShift;
	  return (sad >> DISTORTION_PRECISION_ADJUSTMENT(iBitDepth - 8));
  }
  else	/* iSubStep == 2 case */	/*we can do only 2 rows at a time*/
  {
	  if (iRemWd >= 8)
	  {
		  for (i = 0; i < ht / 4; i++)
		  {
			  Int iCount = 0;
			  for (j = iRemWd; j > 7; j -= 8)
			  {
				  src_8x16b = _mm_loadu_si128((__m128i *) (pu2_inp + iCount));
				  src_8x16b_1 = _mm_loadu_si128((__m128i *) (pu2_inp + iCount + iSrcStride));

				  pred_8x16b = _mm_loadu_si128((__m128i *) (pu2_ref + iCount));
				  pred_8x16b_1 = _mm_loadu_si128((__m128i *) (pu2_ref + iCount + iRefStride));

				  temp = _mm_sub_epi16(src_8x16b, pred_8x16b);
				  temp_1 = _mm_sub_epi16(src_8x16b_1, pred_8x16b_1);

				  temp = _mm_abs_epi16(temp);
				  temp_1 = _mm_abs_epi16(temp_1);

				  temp = _mm_adds_epu16(temp, temp_1);

				  temp_1 = _mm_unpackhi_epi16(temp, temp_dummy);
				  temp_3 = _mm_unpacklo_epi16(temp, temp_dummy);

				  temp = _mm_add_epi32(temp_1, temp_3);

				  result = _mm_add_epi32(result, temp);

				  iCount += 8;
			  }
			  pu2_inp += (2 * iSrcStride);
			  pu2_ref += (2 * iRefStride);
		  }
	  }

	  iRemWd &= 0x7;
	  pu2_inp = pcDtParam->pOrg + ((iCols / 8) * 8);
	  pu2_ref = pcDtParam->pCur + ((iCols / 8) * 8);

	  if (iRemWd > 3)
	  {
		  for (i = 0; i < (ht / 4); i++)
		  {
			  src_8x16b = _mm_loadl_epi64((__m128i *) (pu2_inp));
			  src_8x16b_1 = _mm_loadl_epi64((__m128i *) (pu2_inp + iSrcStride));

			  pred_8x16b = _mm_loadl_epi64((__m128i *) (pu2_ref));
			  pred_8x16b_1 = _mm_loadl_epi64((__m128i *) (pu2_ref + iRefStride));

			  src_8x16b = _mm_unpacklo_epi16(src_8x16b, src_8x16b_1);
			  pred_8x16b = _mm_unpacklo_epi16(pred_8x16b, pred_8x16b_1);

			  temp = _mm_sub_epi16(src_8x16b, pred_8x16b);

			  temp = _mm_abs_epi16(temp);

			  temp_1 = _mm_unpackhi_epi16(temp, temp_dummy);
			  temp_2 = _mm_unpacklo_epi16(temp, temp_dummy);

			  temp_1 = _mm_add_epi32(temp_1, temp_2);
			  result = _mm_add_epi32(result, temp_1);

			  pu2_inp += 2 * (iSrcStride);
			  pu2_ref += 2 * (iRefStride);

		  }
	  }

	  {
		  Int *val = (Int*)&result;
		  sad = val[0] + val[1] + val[2] + val[3];
	  }

	  sad <<= iSubShift;
	  return (sad >> DISTORTION_PRECISION_ADJUSTMENT(iBitDepth - 8));
  }
#endif
}

Distortion TComRdCost::xGetSAD4( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD4( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine4n16b( piOrg , piCur , 4 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD8( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD8( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg      = pcDtParam->pOrg;
  const Pel* piCur      = pcDtParam->pCur;
  Int  iRows      = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 8 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD16( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD16( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 16 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );
    uiSum += abs( piOrg[12] - piCur[12] );
    uiSum += abs( piOrg[13] - piCur[13] );
    uiSum += abs( piOrg[14] - piCur[14] );
    uiSum += abs( piOrg[15] - piCur[15] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD12( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD12( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD16N( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD16N( pcDtParam );
  }
#endif
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , iCols );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    for (Int n = 0; n < iCols; n+=16 )
    {
      uiSum += abs( piOrg[n+ 0] - piCur[n+ 0] );
      uiSum += abs( piOrg[n+ 1] - piCur[n+ 1] );
      uiSum += abs( piOrg[n+ 2] - piCur[n+ 2] );
      uiSum += abs( piOrg[n+ 3] - piCur[n+ 3] );
      uiSum += abs( piOrg[n+ 4] - piCur[n+ 4] );
      uiSum += abs( piOrg[n+ 5] - piCur[n+ 5] );
      uiSum += abs( piOrg[n+ 6] - piCur[n+ 6] );
      uiSum += abs( piOrg[n+ 7] - piCur[n+ 7] );
      uiSum += abs( piOrg[n+ 8] - piCur[n+ 8] );
      uiSum += abs( piOrg[n+ 9] - piCur[n+ 9] );
      uiSum += abs( piOrg[n+10] - piCur[n+10] );
      uiSum += abs( piOrg[n+11] - piCur[n+11] );
      uiSum += abs( piOrg[n+12] - piCur[n+12] );
      uiSum += abs( piOrg[n+13] - piCur[n+13] );
      uiSum += abs( piOrg[n+14] - piCur[n+14] );
      uiSum += abs( piOrg[n+15] - piCur[n+15] );
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD32( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD32( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 32 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );
    uiSum += abs( piOrg[12] - piCur[12] );
    uiSum += abs( piOrg[13] - piCur[13] );
    uiSum += abs( piOrg[14] - piCur[14] );
    uiSum += abs( piOrg[15] - piCur[15] );
    uiSum += abs( piOrg[16] - piCur[16] );
    uiSum += abs( piOrg[17] - piCur[17] );
    uiSum += abs( piOrg[18] - piCur[18] );
    uiSum += abs( piOrg[19] - piCur[19] );
    uiSum += abs( piOrg[20] - piCur[20] );
    uiSum += abs( piOrg[21] - piCur[21] );
    uiSum += abs( piOrg[22] - piCur[22] );
    uiSum += abs( piOrg[23] - piCur[23] );
    uiSum += abs( piOrg[24] - piCur[24] );
    uiSum += abs( piOrg[25] - piCur[25] );
    uiSum += abs( piOrg[26] - piCur[26] );
    uiSum += abs( piOrg[27] - piCur[27] );
    uiSum += abs( piOrg[28] - piCur[28] );
    uiSum += abs( piOrg[29] - piCur[29] );
    uiSum += abs( piOrg[30] - piCur[30] );
    uiSum += abs( piOrg[31] - piCur[31] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD24( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD24( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 24 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );
    uiSum += abs( piOrg[12] - piCur[12] );
    uiSum += abs( piOrg[13] - piCur[13] );
    uiSum += abs( piOrg[14] - piCur[14] );
    uiSum += abs( piOrg[15] - piCur[15] );
    uiSum += abs( piOrg[16] - piCur[16] );
    uiSum += abs( piOrg[17] - piCur[17] );
    uiSum += abs( piOrg[18] - piCur[18] );
    uiSum += abs( piOrg[19] - piCur[19] );
    uiSum += abs( piOrg[20] - piCur[20] );
    uiSum += abs( piOrg[21] - piCur[21] );
    uiSum += abs( piOrg[22] - piCur[22] );
    uiSum += abs( piOrg[23] - piCur[23] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD64( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD64( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 64 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );
    uiSum += abs( piOrg[12] - piCur[12] );
    uiSum += abs( piOrg[13] - piCur[13] );
    uiSum += abs( piOrg[14] - piCur[14] );
    uiSum += abs( piOrg[15] - piCur[15] );
    uiSum += abs( piOrg[16] - piCur[16] );
    uiSum += abs( piOrg[17] - piCur[17] );
    uiSum += abs( piOrg[18] - piCur[18] );
    uiSum += abs( piOrg[19] - piCur[19] );
    uiSum += abs( piOrg[20] - piCur[20] );
    uiSum += abs( piOrg[21] - piCur[21] );
    uiSum += abs( piOrg[22] - piCur[22] );
    uiSum += abs( piOrg[23] - piCur[23] );
    uiSum += abs( piOrg[24] - piCur[24] );
    uiSum += abs( piOrg[25] - piCur[25] );
    uiSum += abs( piOrg[26] - piCur[26] );
    uiSum += abs( piOrg[27] - piCur[27] );
    uiSum += abs( piOrg[28] - piCur[28] );
    uiSum += abs( piOrg[29] - piCur[29] );
    uiSum += abs( piOrg[30] - piCur[30] );
    uiSum += abs( piOrg[31] - piCur[31] );
    uiSum += abs( piOrg[32] - piCur[32] );
    uiSum += abs( piOrg[33] - piCur[33] );
    uiSum += abs( piOrg[34] - piCur[34] );
    uiSum += abs( piOrg[35] - piCur[35] );
    uiSum += abs( piOrg[36] - piCur[36] );
    uiSum += abs( piOrg[37] - piCur[37] );
    uiSum += abs( piOrg[38] - piCur[38] );
    uiSum += abs( piOrg[39] - piCur[39] );
    uiSum += abs( piOrg[40] - piCur[40] );
    uiSum += abs( piOrg[41] - piCur[41] );
    uiSum += abs( piOrg[42] - piCur[42] );
    uiSum += abs( piOrg[43] - piCur[43] );
    uiSum += abs( piOrg[44] - piCur[44] );
    uiSum += abs( piOrg[45] - piCur[45] );
    uiSum += abs( piOrg[46] - piCur[46] );
    uiSum += abs( piOrg[47] - piCur[47] );
    uiSum += abs( piOrg[48] - piCur[48] );
    uiSum += abs( piOrg[49] - piCur[49] );
    uiSum += abs( piOrg[50] - piCur[50] );
    uiSum += abs( piOrg[51] - piCur[51] );
    uiSum += abs( piOrg[52] - piCur[52] );
    uiSum += abs( piOrg[53] - piCur[53] );
    uiSum += abs( piOrg[54] - piCur[54] );
    uiSum += abs( piOrg[55] - piCur[55] );
    uiSum += abs( piOrg[56] - piCur[56] );
    uiSum += abs( piOrg[57] - piCur[57] );
    uiSum += abs( piOrg[58] - piCur[58] );
    uiSum += abs( piOrg[59] - piCur[59] );
    uiSum += abs( piOrg[60] - piCur[60] );
    uiSum += abs( piOrg[61] - piCur[61] );
    uiSum += abs( piOrg[62] - piCur[62] );
    uiSum += abs( piOrg[63] - piCur[63] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

Distortion TComRdCost::xGetSAD48( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRSAD48( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSADw( pcDtParam );
  }

  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  Distortion uiSum = 0;

#if COM16_C806_SIMD_OPT
  if( pcDtParam->bitDepth <= 10 )
  {
    for( ; iRows != 0; iRows-=iSubStep )
    {
      uiSum += simdSADLine8n16b( piOrg , piCur , 48 );
      piOrg += iStrideOrg;
      piCur += iStrideCur;
    }
  }
  else
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] );
    uiSum += abs( piOrg[1] - piCur[1] );
    uiSum += abs( piOrg[2] - piCur[2] );
    uiSum += abs( piOrg[3] - piCur[3] );
    uiSum += abs( piOrg[4] - piCur[4] );
    uiSum += abs( piOrg[5] - piCur[5] );
    uiSum += abs( piOrg[6] - piCur[6] );
    uiSum += abs( piOrg[7] - piCur[7] );
    uiSum += abs( piOrg[8] - piCur[8] );
    uiSum += abs( piOrg[9] - piCur[9] );
    uiSum += abs( piOrg[10] - piCur[10] );
    uiSum += abs( piOrg[11] - piCur[11] );
    uiSum += abs( piOrg[12] - piCur[12] );
    uiSum += abs( piOrg[13] - piCur[13] );
    uiSum += abs( piOrg[14] - piCur[14] );
    uiSum += abs( piOrg[15] - piCur[15] );
    uiSum += abs( piOrg[16] - piCur[16] );
    uiSum += abs( piOrg[17] - piCur[17] );
    uiSum += abs( piOrg[18] - piCur[18] );
    uiSum += abs( piOrg[19] - piCur[19] );
    uiSum += abs( piOrg[20] - piCur[20] );
    uiSum += abs( piOrg[21] - piCur[21] );
    uiSum += abs( piOrg[22] - piCur[22] );
    uiSum += abs( piOrg[23] - piCur[23] );
    uiSum += abs( piOrg[24] - piCur[24] );
    uiSum += abs( piOrg[25] - piCur[25] );
    uiSum += abs( piOrg[26] - piCur[26] );
    uiSum += abs( piOrg[27] - piCur[27] );
    uiSum += abs( piOrg[28] - piCur[28] );
    uiSum += abs( piOrg[29] - piCur[29] );
    uiSum += abs( piOrg[30] - piCur[30] );
    uiSum += abs( piOrg[31] - piCur[31] );
    uiSum += abs( piOrg[32] - piCur[32] );
    uiSum += abs( piOrg[33] - piCur[33] );
    uiSum += abs( piOrg[34] - piCur[34] );
    uiSum += abs( piOrg[35] - piCur[35] );
    uiSum += abs( piOrg[36] - piCur[36] );
    uiSum += abs( piOrg[37] - piCur[37] );
    uiSum += abs( piOrg[38] - piCur[38] );
    uiSum += abs( piOrg[39] - piCur[39] );
    uiSum += abs( piOrg[40] - piCur[40] );
    uiSum += abs( piOrg[41] - piCur[41] );
    uiSum += abs( piOrg[42] - piCur[42] );
    uiSum += abs( piOrg[43] - piCur[43] );
    uiSum += abs( piOrg[44] - piCur[44] );
    uiSum += abs( piOrg[45] - piCur[45] );
    uiSum += abs( piOrg[46] - piCur[46] );
    uiSum += abs( piOrg[47] - piCur[47] );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if COM16_C806_SIMD_OPT
  }
#endif

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

#if VCEG_AZ06_IC
UInt TComRdCost::xGetMRSAD( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideCur = pcDtParam->iStrideCur;
  Int  iStrideOrg = pcDtParam->iStrideOrg;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows--)
    {
      for (Int n = 0; n < iCols; n++)
      {
        iCurAvg += piCur[n];
      }
      piCur += iStrideCur;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n++ )
    {
      iOrigAvg += piOrg[n];
      iCurAvg  += piCur[n];
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = (iOrigAvg - iCurAvg)/iCols/iRows;

  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n++ )
    {
      uiSum += abs( piOrg[n] - piCur[n] - iDeltaC );
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD4( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/4) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD8( DistParam* pcDtParam )
{
  Pel* piOrg      = pcDtParam->pOrg;
  Pel* piCur      = pcDtParam->pCur;
  Int  iRows      = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];
    iOrigAvg += piOrg[4];
    iOrigAvg += piOrg[5];
    iOrigAvg += piOrg[6];
    iOrigAvg += piOrg[7];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];
    iCurAvg  += piCur[4];
    iCurAvg  += piCur[5];
    iCurAvg  += piCur[6];
    iCurAvg  += piCur[7];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/8) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD16( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];
      iCurAvg += piCur[12];
      iCurAvg += piCur[13];
      iCurAvg += piCur[14];
      iCurAvg += piCur[15];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];
    iOrigAvg += piOrg[4];
    iOrigAvg += piOrg[5];
    iOrigAvg += piOrg[6];
    iOrigAvg += piOrg[7];
    iOrigAvg += piOrg[8];
    iOrigAvg += piOrg[9];
    iOrigAvg += piOrg[10];
    iOrigAvg += piOrg[11];
    iOrigAvg += piOrg[12];
    iOrigAvg += piOrg[13];
    iOrigAvg += piOrg[14];
    iOrigAvg += piOrg[15];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];
    iCurAvg  += piCur[4];
    iCurAvg  += piCur[5];
    iCurAvg  += piCur[6];
    iCurAvg  += piCur[7];
    iCurAvg  += piCur[8];
    iCurAvg  += piCur[9];
    iCurAvg  += piCur[10];
    iCurAvg  += piCur[11];
    iCurAvg  += piCur[12];
    iCurAvg  += piCur[13];
    iCurAvg  += piCur[14];
    iCurAvg  += piCur[15];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/16) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );
    uiSum += abs( piOrg[12] - piCur[12] - iDeltaC );
    uiSum += abs( piOrg[13] - piCur[13] - iDeltaC );
    uiSum += abs( piOrg[14] - piCur[14] - iDeltaC );
    uiSum += abs( piOrg[15] - piCur[15] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD12( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];
    iOrigAvg += piOrg[4];
    iOrigAvg += piOrg[5];
    iOrigAvg += piOrg[6];
    iOrigAvg += piOrg[7];
    iOrigAvg += piOrg[8];
    iOrigAvg += piOrg[9];
    iOrigAvg += piOrg[10];
    iOrigAvg += piOrg[11];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];
    iCurAvg  += piCur[4];
    iCurAvg  += piCur[5];
    iCurAvg  += piCur[6];
    iCurAvg  += piCur[7];
    iCurAvg  += piCur[8];
    iCurAvg  += piCur[9];
    iCurAvg  += piCur[10];
    iCurAvg  += piCur[11];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/12) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD16N( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0, uiColCnt = (iCols-1)/16 + 1;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      for (Int n = 0; n < iCols; n += 16)
      {
        iCurAvg += piCur[n + 0];
        iCurAvg += piCur[n + 1];
        iCurAvg += piCur[n + 2];
        iCurAvg += piCur[n + 3];
        iCurAvg += piCur[n + 4];
        iCurAvg += piCur[n + 5];
        iCurAvg += piCur[n + 6];
        iCurAvg += piCur[n + 7];
        iCurAvg += piCur[n + 8];
        iCurAvg += piCur[n + 9];
        iCurAvg += piCur[n + 10];
        iCurAvg += piCur[n + 11];
        iCurAvg += piCur[n + 12];
        iCurAvg += piCur[n + 13];
        iCurAvg += piCur[n + 14];
        iCurAvg += piCur[n + 15];
      }
      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    for (Int n = 0; n < iCols; n+=16 )
    {
      iOrigAvg += piOrg[n + 0];
      iOrigAvg += piOrg[n + 1];
      iOrigAvg += piOrg[n + 2];
      iOrigAvg += piOrg[n + 3];
      iOrigAvg += piOrg[n + 4];
      iOrigAvg += piOrg[n + 5];
      iOrigAvg += piOrg[n + 6];
      iOrigAvg += piOrg[n + 7];
      iOrigAvg += piOrg[n + 8];
      iOrigAvg += piOrg[n + 9];
      iOrigAvg += piOrg[n + 10];
      iOrigAvg += piOrg[n + 11];
      iOrigAvg += piOrg[n + 12];
      iOrigAvg += piOrg[n + 13];
      iOrigAvg += piOrg[n + 14];
      iOrigAvg += piOrg[n + 15];

      iCurAvg  += piCur[n + 0];
      iCurAvg  += piCur[n + 1];
      iCurAvg  += piCur[n + 2];
      iCurAvg  += piCur[n + 3];
      iCurAvg  += piCur[n + 4];
      iCurAvg  += piCur[n + 5];
      iCurAvg  += piCur[n + 6];
      iCurAvg  += piCur[n + 7];
      iCurAvg  += piCur[n + 8];
      iCurAvg  += piCur[n + 9];
      iCurAvg  += piCur[n + 10];
      iCurAvg  += piCur[n + 11];
      iCurAvg  += piCur[n + 12];
      iCurAvg  += piCur[n + 13];
      iCurAvg  += piCur[n + 14];
      iCurAvg  += piCur[n + 15];
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif
  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = (uiRowCnt && uiColCnt) ? ((iOrigAvg - iCurAvg)/uiRowCnt/uiColCnt/16) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    for (Int n = 0; n < iCols; n+=16 )
    {
      uiSum += abs( piOrg[n+ 0] - piCur[n+ 0] - iDeltaC );
      uiSum += abs( piOrg[n+ 1] - piCur[n+ 1] - iDeltaC );
      uiSum += abs( piOrg[n+ 2] - piCur[n+ 2] - iDeltaC );
      uiSum += abs( piOrg[n+ 3] - piCur[n+ 3] - iDeltaC );
      uiSum += abs( piOrg[n+ 4] - piCur[n+ 4] - iDeltaC );
      uiSum += abs( piOrg[n+ 5] - piCur[n+ 5] - iDeltaC );
      uiSum += abs( piOrg[n+ 6] - piCur[n+ 6] - iDeltaC );
      uiSum += abs( piOrg[n+ 7] - piCur[n+ 7] - iDeltaC );
      uiSum += abs( piOrg[n+ 8] - piCur[n+ 8] - iDeltaC );
      uiSum += abs( piOrg[n+ 9] - piCur[n+ 9] - iDeltaC );
      uiSum += abs( piOrg[n+10] - piCur[n+10] - iDeltaC );
      uiSum += abs( piOrg[n+11] - piCur[n+11] - iDeltaC );
      uiSum += abs( piOrg[n+12] - piCur[n+12] - iDeltaC );
      uiSum += abs( piOrg[n+13] - piCur[n+13] - iDeltaC );
      uiSum += abs( piOrg[n+14] - piCur[n+14] - iDeltaC );
      uiSum += abs( piOrg[n+15] - piCur[n+15] - iDeltaC );
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD32( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];
      iCurAvg += piCur[12];
      iCurAvg += piCur[13];
      iCurAvg += piCur[14];
      iCurAvg += piCur[15];
      iCurAvg += piCur[16];
      iCurAvg += piCur[17];
      iCurAvg += piCur[18];
      iCurAvg += piCur[19];
      iCurAvg += piCur[20];
      iCurAvg += piCur[21];
      iCurAvg += piCur[22];
      iCurAvg += piCur[23];
      iCurAvg += piCur[24];
      iCurAvg += piCur[25];
      iCurAvg += piCur[26];
      iCurAvg += piCur[27];
      iCurAvg += piCur[28];
      iCurAvg += piCur[29];
      iCurAvg += piCur[30];
      iCurAvg += piCur[31];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];
    iOrigAvg += piOrg[4];
    iOrigAvg += piOrg[5];
    iOrigAvg += piOrg[6];
    iOrigAvg += piOrg[7];
    iOrigAvg += piOrg[8];
    iOrigAvg += piOrg[9];
    iOrigAvg += piOrg[10];
    iOrigAvg += piOrg[11];
    iOrigAvg += piOrg[12];
    iOrigAvg += piOrg[13];
    iOrigAvg += piOrg[14];
    iOrigAvg += piOrg[15];
    iOrigAvg += piOrg[16];
    iOrigAvg += piOrg[17];
    iOrigAvg += piOrg[18];
    iOrigAvg += piOrg[19];
    iOrigAvg += piOrg[20];
    iOrigAvg += piOrg[21];
    iOrigAvg += piOrg[22];
    iOrigAvg += piOrg[23];
    iOrigAvg += piOrg[24];
    iOrigAvg += piOrg[25];
    iOrigAvg += piOrg[26];
    iOrigAvg += piOrg[27];
    iOrigAvg += piOrg[28];
    iOrigAvg += piOrg[29];
    iOrigAvg += piOrg[30];
    iOrigAvg += piOrg[31];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];
    iCurAvg  += piCur[4];
    iCurAvg  += piCur[5];
    iCurAvg  += piCur[6];
    iCurAvg  += piCur[7];
    iCurAvg  += piCur[8];
    iCurAvg  += piCur[9];
    iCurAvg  += piCur[10];
    iCurAvg  += piCur[11];
    iCurAvg  += piCur[12];
    iCurAvg  += piCur[13];
    iCurAvg  += piCur[14];
    iCurAvg  += piCur[15];
    iCurAvg  += piCur[16];
    iCurAvg  += piCur[17];
    iCurAvg  += piCur[18];
    iCurAvg  += piCur[19];
    iCurAvg  += piCur[20];
    iCurAvg  += piCur[21];
    iCurAvg  += piCur[22];
    iCurAvg  += piCur[23];
    iCurAvg  += piCur[24];
    iCurAvg  += piCur[25];
    iCurAvg  += piCur[26];
    iCurAvg  += piCur[27];
    iCurAvg  += piCur[28];
    iCurAvg  += piCur[29];
    iCurAvg  += piCur[30];
    iCurAvg  += piCur[31];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/32) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );
    uiSum += abs( piOrg[12] - piCur[12] - iDeltaC );
    uiSum += abs( piOrg[13] - piCur[13] - iDeltaC );
    uiSum += abs( piOrg[14] - piCur[14] - iDeltaC );
    uiSum += abs( piOrg[15] - piCur[15] - iDeltaC );
    uiSum += abs( piOrg[16] - piCur[16] - iDeltaC );
    uiSum += abs( piOrg[17] - piCur[17] - iDeltaC );
    uiSum += abs( piOrg[18] - piCur[18] - iDeltaC );
    uiSum += abs( piOrg[19] - piCur[19] - iDeltaC );
    uiSum += abs( piOrg[20] - piCur[20] - iDeltaC );
    uiSum += abs( piOrg[21] - piCur[21] - iDeltaC );
    uiSum += abs( piOrg[22] - piCur[22] - iDeltaC );
    uiSum += abs( piOrg[23] - piCur[23] - iDeltaC );
    uiSum += abs( piOrg[24] - piCur[24] - iDeltaC );
    uiSum += abs( piOrg[25] - piCur[25] - iDeltaC );
    uiSum += abs( piOrg[26] - piCur[26] - iDeltaC );
    uiSum += abs( piOrg[27] - piCur[27] - iDeltaC );
    uiSum += abs( piOrg[28] - piCur[28] - iDeltaC );
    uiSum += abs( piOrg[29] - piCur[29] - iDeltaC );
    uiSum += abs( piOrg[30] - piCur[30] - iDeltaC );
    uiSum += abs( piOrg[31] - piCur[31] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD24( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];
      iCurAvg += piCur[12];
      iCurAvg += piCur[13];
      iCurAvg += piCur[14];
      iCurAvg += piCur[15];
      iCurAvg += piCur[16];
      iCurAvg += piCur[17];
      iCurAvg += piCur[18];
      iCurAvg += piCur[19];
      iCurAvg += piCur[20];
      iCurAvg += piCur[21];
      iCurAvg += piCur[22];
      iCurAvg += piCur[23];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0];
    iOrigAvg += piOrg[1];
    iOrigAvg += piOrg[2];
    iOrigAvg += piOrg[3];
    iOrigAvg += piOrg[4];
    iOrigAvg += piOrg[5];
    iOrigAvg += piOrg[6];
    iOrigAvg += piOrg[7];
    iOrigAvg += piOrg[8];
    iOrigAvg += piOrg[9];
    iOrigAvg += piOrg[10];
    iOrigAvg += piOrg[11];
    iOrigAvg += piOrg[12];
    iOrigAvg += piOrg[13];
    iOrigAvg += piOrg[14];
    iOrigAvg += piOrg[15];
    iOrigAvg += piOrg[16];
    iOrigAvg += piOrg[17];
    iOrigAvg += piOrg[18];
    iOrigAvg += piOrg[19];
    iOrigAvg += piOrg[20];
    iOrigAvg += piOrg[21];
    iOrigAvg += piOrg[22];
    iOrigAvg += piOrg[23];

    iCurAvg  += piCur[0];
    iCurAvg  += piCur[1];
    iCurAvg  += piCur[2];
    iCurAvg  += piCur[3];
    iCurAvg  += piCur[4];
    iCurAvg  += piCur[5];
    iCurAvg  += piCur[6];
    iCurAvg  += piCur[7];
    iCurAvg  += piCur[8];
    iCurAvg  += piCur[9];
    iCurAvg  += piCur[10];
    iCurAvg  += piCur[11];
    iCurAvg  += piCur[12];
    iCurAvg  += piCur[13];
    iCurAvg  += piCur[14];
    iCurAvg  += piCur[15];
    iCurAvg  += piCur[16];
    iCurAvg  += piCur[17];
    iCurAvg  += piCur[18];
    iCurAvg  += piCur[19];
    iCurAvg  += piCur[20];
    iCurAvg  += piCur[21];
    iCurAvg  += piCur[22];
    iCurAvg  += piCur[23];

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/24) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );
    uiSum += abs( piOrg[12] - piCur[12] - iDeltaC );
    uiSum += abs( piOrg[13] - piCur[13] - iDeltaC );
    uiSum += abs( piOrg[14] - piCur[14] - iDeltaC );
    uiSum += abs( piOrg[15] - piCur[15] - iDeltaC );
    uiSum += abs( piOrg[16] - piCur[16] - iDeltaC );
    uiSum += abs( piOrg[17] - piCur[17] - iDeltaC );
    uiSum += abs( piOrg[18] - piCur[18] - iDeltaC );
    uiSum += abs( piOrg[19] - piCur[19] - iDeltaC );
    uiSum += abs( piOrg[20] - piCur[20] - iDeltaC );
    uiSum += abs( piOrg[21] - piCur[21] - iDeltaC );
    uiSum += abs( piOrg[22] - piCur[22] - iDeltaC );
    uiSum += abs( piOrg[23] - piCur[23] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD64( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];
      iCurAvg += piCur[12];
      iCurAvg += piCur[13];
      iCurAvg += piCur[14];
      iCurAvg += piCur[15];
      iCurAvg += piCur[16];
      iCurAvg += piCur[17];
      iCurAvg += piCur[18];
      iCurAvg += piCur[19];
      iCurAvg += piCur[20];
      iCurAvg += piCur[21];
      iCurAvg += piCur[22];
      iCurAvg += piCur[23];
      iCurAvg += piCur[24];
      iCurAvg += piCur[25];
      iCurAvg += piCur[26];
      iCurAvg += piCur[27];
      iCurAvg += piCur[28];
      iCurAvg += piCur[29];
      iCurAvg += piCur[30];
      iCurAvg += piCur[31];
      iCurAvg += piCur[32];
      iCurAvg += piCur[33];
      iCurAvg += piCur[34];
      iCurAvg += piCur[35];
      iCurAvg += piCur[36];
      iCurAvg += piCur[37];
      iCurAvg += piCur[38];
      iCurAvg += piCur[39];
      iCurAvg += piCur[40];
      iCurAvg += piCur[41];
      iCurAvg += piCur[42];
      iCurAvg += piCur[43];
      iCurAvg += piCur[44];
      iCurAvg += piCur[45];
      iCurAvg += piCur[46];
      iCurAvg += piCur[47];
      iCurAvg += piCur[48];
      iCurAvg += piCur[49];
      iCurAvg += piCur[50];
      iCurAvg += piCur[51];
      iCurAvg += piCur[52];
      iCurAvg += piCur[53];
      iCurAvg += piCur[54];
      iCurAvg += piCur[55];
      iCurAvg += piCur[56];
      iCurAvg += piCur[57];
      iCurAvg += piCur[58];
      iCurAvg += piCur[59];
      iCurAvg += piCur[60];
      iCurAvg += piCur[61];
      iCurAvg += piCur[62];
      iCurAvg += piCur[63];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0] ;
    iOrigAvg += piOrg[1] ;
    iOrigAvg += piOrg[2] ;
    iOrigAvg += piOrg[3] ;
    iOrigAvg += piOrg[4] ;
    iOrigAvg += piOrg[5] ;
    iOrigAvg += piOrg[6] ;
    iOrigAvg += piOrg[7] ;
    iOrigAvg += piOrg[8] ;
    iOrigAvg += piOrg[9] ;
    iOrigAvg += piOrg[10] ;
    iOrigAvg += piOrg[11] ;
    iOrigAvg += piOrg[12] ;
    iOrigAvg += piOrg[13] ;
    iOrigAvg += piOrg[14] ;
    iOrigAvg += piOrg[15] ;
    iOrigAvg += piOrg[16] ;
    iOrigAvg += piOrg[17] ;
    iOrigAvg += piOrg[18] ;
    iOrigAvg += piOrg[19] ;
    iOrigAvg += piOrg[20] ;
    iOrigAvg += piOrg[21] ;
    iOrigAvg += piOrg[22] ;
    iOrigAvg += piOrg[23] ;
    iOrigAvg += piOrg[24] ;
    iOrigAvg += piOrg[25] ;
    iOrigAvg += piOrg[26] ;
    iOrigAvg += piOrg[27] ;
    iOrigAvg += piOrg[28] ;
    iOrigAvg += piOrg[29] ;
    iOrigAvg += piOrg[30] ;
    iOrigAvg += piOrg[31] ;
    iOrigAvg += piOrg[32] ;
    iOrigAvg += piOrg[33] ;
    iOrigAvg += piOrg[34] ;
    iOrigAvg += piOrg[35] ;
    iOrigAvg += piOrg[36] ;
    iOrigAvg += piOrg[37] ;
    iOrigAvg += piOrg[38] ;
    iOrigAvg += piOrg[39] ;
    iOrigAvg += piOrg[40] ;
    iOrigAvg += piOrg[41] ;
    iOrigAvg += piOrg[42] ;
    iOrigAvg += piOrg[43] ;
    iOrigAvg += piOrg[44] ;
    iOrigAvg += piOrg[45] ;
    iOrigAvg += piOrg[46] ;
    iOrigAvg += piOrg[47] ;
    iOrigAvg += piOrg[48] ;
    iOrigAvg += piOrg[49] ;
    iOrigAvg += piOrg[50] ;
    iOrigAvg += piOrg[51] ;
    iOrigAvg += piOrg[52] ;
    iOrigAvg += piOrg[53] ;
    iOrigAvg += piOrg[54] ;
    iOrigAvg += piOrg[55] ;
    iOrigAvg += piOrg[56] ;
    iOrigAvg += piOrg[57] ;
    iOrigAvg += piOrg[58] ;
    iOrigAvg += piOrg[59] ;
    iOrigAvg += piOrg[60] ;
    iOrigAvg += piOrg[61] ;
    iOrigAvg += piOrg[62] ;
    iOrigAvg += piOrg[63] ;

    iCurAvg += piCur[0] ;
    iCurAvg += piCur[1] ;
    iCurAvg += piCur[2] ;
    iCurAvg += piCur[3] ;
    iCurAvg += piCur[4] ;
    iCurAvg += piCur[5] ;
    iCurAvg += piCur[6] ;
    iCurAvg += piCur[7] ;
    iCurAvg += piCur[8] ;
    iCurAvg += piCur[9] ;
    iCurAvg += piCur[10] ;
    iCurAvg += piCur[11] ;
    iCurAvg += piCur[12] ;
    iCurAvg += piCur[13] ;
    iCurAvg += piCur[14] ;
    iCurAvg += piCur[15] ;
    iCurAvg += piCur[16] ;
    iCurAvg += piCur[17] ;
    iCurAvg += piCur[18] ;
    iCurAvg += piCur[19] ;
    iCurAvg += piCur[20] ;
    iCurAvg += piCur[21] ;
    iCurAvg += piCur[22] ;
    iCurAvg += piCur[23] ;
    iCurAvg += piCur[24] ;
    iCurAvg += piCur[25] ;
    iCurAvg += piCur[26] ;
    iCurAvg += piCur[27] ;
    iCurAvg += piCur[28] ;
    iCurAvg += piCur[29] ;
    iCurAvg += piCur[30] ;
    iCurAvg += piCur[31] ;
    iCurAvg += piCur[32] ;
    iCurAvg += piCur[33] ;
    iCurAvg += piCur[34] ;
    iCurAvg += piCur[35] ;
    iCurAvg += piCur[36] ;
    iCurAvg += piCur[37] ;
    iCurAvg += piCur[38] ;
    iCurAvg += piCur[39] ;
    iCurAvg += piCur[40] ;
    iCurAvg += piCur[41] ;
    iCurAvg += piCur[42] ;
    iCurAvg += piCur[43] ;
    iCurAvg += piCur[44] ;
    iCurAvg += piCur[45] ;
    iCurAvg += piCur[46] ;
    iCurAvg += piCur[47] ;
    iCurAvg += piCur[48] ;
    iCurAvg += piCur[49] ;
    iCurAvg += piCur[50] ;
    iCurAvg += piCur[51] ;
    iCurAvg += piCur[52] ;
    iCurAvg += piCur[53] ;
    iCurAvg += piCur[54] ;
    iCurAvg += piCur[55] ;
    iCurAvg += piCur[56] ;
    iCurAvg += piCur[57] ;
    iCurAvg += piCur[58] ;
    iCurAvg += piCur[59] ;
    iCurAvg += piCur[60] ;
    iCurAvg += piCur[61] ;
    iCurAvg += piCur[62] ;
    iCurAvg += piCur[63] ;

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
  }
#endif

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/64) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );
    uiSum += abs( piOrg[12] - piCur[12] - iDeltaC );
    uiSum += abs( piOrg[13] - piCur[13] - iDeltaC );
    uiSum += abs( piOrg[14] - piCur[14] - iDeltaC );
    uiSum += abs( piOrg[15] - piCur[15] - iDeltaC );
    uiSum += abs( piOrg[16] - piCur[16] - iDeltaC );
    uiSum += abs( piOrg[17] - piCur[17] - iDeltaC );
    uiSum += abs( piOrg[18] - piCur[18] - iDeltaC );
    uiSum += abs( piOrg[19] - piCur[19] - iDeltaC );
    uiSum += abs( piOrg[20] - piCur[20] - iDeltaC );
    uiSum += abs( piOrg[21] - piCur[21] - iDeltaC );
    uiSum += abs( piOrg[22] - piCur[22] - iDeltaC );
    uiSum += abs( piOrg[23] - piCur[23] - iDeltaC );
    uiSum += abs( piOrg[24] - piCur[24] - iDeltaC );
    uiSum += abs( piOrg[25] - piCur[25] - iDeltaC );
    uiSum += abs( piOrg[26] - piCur[26] - iDeltaC );
    uiSum += abs( piOrg[27] - piCur[27] - iDeltaC );
    uiSum += abs( piOrg[28] - piCur[28] - iDeltaC );
    uiSum += abs( piOrg[29] - piCur[29] - iDeltaC );
    uiSum += abs( piOrg[30] - piCur[30] - iDeltaC );
    uiSum += abs( piOrg[31] - piCur[31] - iDeltaC );
    uiSum += abs( piOrg[32] - piCur[32] - iDeltaC );
    uiSum += abs( piOrg[33] - piCur[33] - iDeltaC );
    uiSum += abs( piOrg[34] - piCur[34] - iDeltaC );
    uiSum += abs( piOrg[35] - piCur[35] - iDeltaC );
    uiSum += abs( piOrg[36] - piCur[36] - iDeltaC );
    uiSum += abs( piOrg[37] - piCur[37] - iDeltaC );
    uiSum += abs( piOrg[38] - piCur[38] - iDeltaC );
    uiSum += abs( piOrg[39] - piCur[39] - iDeltaC );
    uiSum += abs( piOrg[40] - piCur[40] - iDeltaC );
    uiSum += abs( piOrg[41] - piCur[41] - iDeltaC );
    uiSum += abs( piOrg[42] - piCur[42] - iDeltaC );
    uiSum += abs( piOrg[43] - piCur[43] - iDeltaC );
    uiSum += abs( piOrg[44] - piCur[44] - iDeltaC );
    uiSum += abs( piOrg[45] - piCur[45] - iDeltaC );
    uiSum += abs( piOrg[46] - piCur[46] - iDeltaC );
    uiSum += abs( piOrg[47] - piCur[47] - iDeltaC );
    uiSum += abs( piOrg[48] - piCur[48] - iDeltaC );
    uiSum += abs( piOrg[49] - piCur[49] - iDeltaC );
    uiSum += abs( piOrg[50] - piCur[50] - iDeltaC );
    uiSum += abs( piOrg[51] - piCur[51] - iDeltaC );
    uiSum += abs( piOrg[52] - piCur[52] - iDeltaC );
    uiSum += abs( piOrg[53] - piCur[53] - iDeltaC );
    uiSum += abs( piOrg[54] - piCur[54] - iDeltaC );
    uiSum += abs( piOrg[55] - piCur[55] - iDeltaC );
    uiSum += abs( piOrg[56] - piCur[56] - iDeltaC );
    uiSum += abs( piOrg[57] - piCur[57] - iDeltaC );
    uiSum += abs( piOrg[58] - piCur[58] - iDeltaC );
    uiSum += abs( piOrg[59] - piCur[59] - iDeltaC );
    uiSum += abs( piOrg[60] - piCur[60] - iDeltaC );
    uiSum += abs( piOrg[61] - piCur[61] - iDeltaC );
    uiSum += abs( piOrg[62] - piCur[62] - iDeltaC );
    uiSum += abs( piOrg[63] - piCur[63] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}

UInt TComRdCost::xGetMRSAD48( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iSubShift  = pcDtParam->iSubShift;
  Int  iSubStep   = ( 1 << iSubShift );
  Int  iStrideCur = pcDtParam->iStrideCur*iSubStep;
  Int  iStrideOrg = pcDtParam->iStrideOrg*iSubStep;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0, uiRowCnt = 0;
  Int  iDeltaC;

#if DMVR_MRSAD_SPEEDUP
  if (pcDtParam->CalcStrategy == 1)
  {
    iOrigAvg = pcDtParam->OrgMean;
    for (; iRows != 0; iRows -= iSubStep)
    {
      iCurAvg += piCur[0];
      iCurAvg += piCur[1];
      iCurAvg += piCur[2];
      iCurAvg += piCur[3];
      iCurAvg += piCur[4];
      iCurAvg += piCur[5];
      iCurAvg += piCur[6];
      iCurAvg += piCur[7];
      iCurAvg += piCur[8];
      iCurAvg += piCur[9];
      iCurAvg += piCur[10];
      iCurAvg += piCur[11];
      iCurAvg += piCur[12];
      iCurAvg += piCur[13];
      iCurAvg += piCur[14];
      iCurAvg += piCur[15];
      iCurAvg += piCur[16];
      iCurAvg += piCur[17];
      iCurAvg += piCur[18];
      iCurAvg += piCur[19];
      iCurAvg += piCur[20];
      iCurAvg += piCur[21];
      iCurAvg += piCur[22];
      iCurAvg += piCur[23];
      iCurAvg += piCur[24];
      iCurAvg += piCur[25];
      iCurAvg += piCur[26];
      iCurAvg += piCur[27];
      iCurAvg += piCur[28];
      iCurAvg += piCur[29];
      iCurAvg += piCur[30];
      iCurAvg += piCur[31];
      iCurAvg += piCur[32];
      iCurAvg += piCur[33];
      iCurAvg += piCur[34];
      iCurAvg += piCur[35];
      iCurAvg += piCur[36];
      iCurAvg += piCur[37];
      iCurAvg += piCur[38];
      iCurAvg += piCur[39];
      iCurAvg += piCur[40];
      iCurAvg += piCur[41];
      iCurAvg += piCur[42];
      iCurAvg += piCur[43];
      iCurAvg += piCur[44];
      iCurAvg += piCur[45];
      iCurAvg += piCur[46];
      iCurAvg += piCur[47];

      piCur += iStrideCur;
      uiRowCnt++;
    }
  }
  if (pcDtParam->CalcStrategy == 0)
  {
#endif
  for( ; iRows != 0; iRows-=iSubStep )
  {
    iOrigAvg += piOrg[0] ;
    iOrigAvg += piOrg[1] ;
    iOrigAvg += piOrg[2] ;
    iOrigAvg += piOrg[3] ;
    iOrigAvg += piOrg[4] ;
    iOrigAvg += piOrg[5] ;
    iOrigAvg += piOrg[6] ;
    iOrigAvg += piOrg[7] ;
    iOrigAvg += piOrg[8] ;
    iOrigAvg += piOrg[9] ;
    iOrigAvg += piOrg[10] ;
    iOrigAvg += piOrg[11] ;
    iOrigAvg += piOrg[12] ;
    iOrigAvg += piOrg[13] ;
    iOrigAvg += piOrg[14] ;
    iOrigAvg += piOrg[15] ;
    iOrigAvg += piOrg[16] ;
    iOrigAvg += piOrg[17] ;
    iOrigAvg += piOrg[18] ;
    iOrigAvg += piOrg[19] ;
    iOrigAvg += piOrg[20] ;
    iOrigAvg += piOrg[21] ;
    iOrigAvg += piOrg[22] ;
    iOrigAvg += piOrg[23] ;
    iOrigAvg += piOrg[24] ;
    iOrigAvg += piOrg[25] ;
    iOrigAvg += piOrg[26] ;
    iOrigAvg += piOrg[27] ;
    iOrigAvg += piOrg[28] ;
    iOrigAvg += piOrg[29] ;
    iOrigAvg += piOrg[30] ;
    iOrigAvg += piOrg[31] ;
    iOrigAvg += piOrg[32] ;
    iOrigAvg += piOrg[33] ;
    iOrigAvg += piOrg[34] ;
    iOrigAvg += piOrg[35] ;
    iOrigAvg += piOrg[36] ;
    iOrigAvg += piOrg[37] ;
    iOrigAvg += piOrg[38] ;
    iOrigAvg += piOrg[39] ;
    iOrigAvg += piOrg[40] ;
    iOrigAvg += piOrg[41] ;
    iOrigAvg += piOrg[42] ;
    iOrigAvg += piOrg[43] ;
    iOrigAvg += piOrg[44] ;
    iOrigAvg += piOrg[45] ;
    iOrigAvg += piOrg[46] ;
    iOrigAvg += piOrg[47] ;

    iCurAvg += piCur[0] ;
    iCurAvg += piCur[1] ;
    iCurAvg += piCur[2] ;
    iCurAvg += piCur[3] ;
    iCurAvg += piCur[4] ;
    iCurAvg += piCur[5] ;
    iCurAvg += piCur[6] ;
    iCurAvg += piCur[7] ;
    iCurAvg += piCur[8] ;
    iCurAvg += piCur[9] ;
    iCurAvg += piCur[10] ;
    iCurAvg += piCur[11] ;
    iCurAvg += piCur[12] ;
    iCurAvg += piCur[13] ;
    iCurAvg += piCur[14] ;
    iCurAvg += piCur[15] ;
    iCurAvg += piCur[16] ;
    iCurAvg += piCur[17] ;
    iCurAvg += piCur[18] ;
    iCurAvg += piCur[19] ;
    iCurAvg += piCur[20] ;
    iCurAvg += piCur[21] ;
    iCurAvg += piCur[22] ;
    iCurAvg += piCur[23] ;
    iCurAvg += piCur[24] ;
    iCurAvg += piCur[25] ;
    iCurAvg += piCur[26] ;
    iCurAvg += piCur[27] ;
    iCurAvg += piCur[28] ;
    iCurAvg += piCur[29] ;
    iCurAvg += piCur[30] ;
    iCurAvg += piCur[31] ;
    iCurAvg += piCur[32] ;
    iCurAvg += piCur[33] ;
    iCurAvg += piCur[34] ;
    iCurAvg += piCur[35] ;
    iCurAvg += piCur[36] ;
    iCurAvg += piCur[37] ;
    iCurAvg += piCur[38] ;
    iCurAvg += piCur[39] ;
    iCurAvg += piCur[40] ;
    iCurAvg += piCur[41] ;
    iCurAvg += piCur[42] ;
    iCurAvg += piCur[43] ;
    iCurAvg += piCur[44] ;
    iCurAvg += piCur[45] ;
    iCurAvg += piCur[46] ;
    iCurAvg += piCur[47] ;

    piOrg += iStrideOrg;
    piCur += iStrideCur;
    uiRowCnt++;
  }
#if DMVR_MRSAD_SPEEDUP
  pcDtParam->OrgMean = (iOrigAvg);
  pcDtParam->CurrMean = (iCurAvg);
    }
#endif
  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;
  iRows   = pcDtParam->iRows;

  iDeltaC = uiRowCnt ? ((iOrigAvg - iCurAvg)/uiRowCnt/48) : 0;

  for( ; iRows != 0; iRows-=iSubStep )
  {
    uiSum += abs( piOrg[0] - piCur[0] - iDeltaC );
    uiSum += abs( piOrg[1] - piCur[1] - iDeltaC );
    uiSum += abs( piOrg[2] - piCur[2] - iDeltaC );
    uiSum += abs( piOrg[3] - piCur[3] - iDeltaC );
    uiSum += abs( piOrg[4] - piCur[4] - iDeltaC );
    uiSum += abs( piOrg[5] - piCur[5] - iDeltaC );
    uiSum += abs( piOrg[6] - piCur[6] - iDeltaC );
    uiSum += abs( piOrg[7] - piCur[7] - iDeltaC );
    uiSum += abs( piOrg[8] - piCur[8] - iDeltaC );
    uiSum += abs( piOrg[9] - piCur[9] - iDeltaC );
    uiSum += abs( piOrg[10] - piCur[10] - iDeltaC );
    uiSum += abs( piOrg[11] - piCur[11] - iDeltaC );
    uiSum += abs( piOrg[12] - piCur[12] - iDeltaC );
    uiSum += abs( piOrg[13] - piCur[13] - iDeltaC );
    uiSum += abs( piOrg[14] - piCur[14] - iDeltaC );
    uiSum += abs( piOrg[15] - piCur[15] - iDeltaC );
    uiSum += abs( piOrg[16] - piCur[16] - iDeltaC );
    uiSum += abs( piOrg[17] - piCur[17] - iDeltaC );
    uiSum += abs( piOrg[18] - piCur[18] - iDeltaC );
    uiSum += abs( piOrg[19] - piCur[19] - iDeltaC );
    uiSum += abs( piOrg[20] - piCur[20] - iDeltaC );
    uiSum += abs( piOrg[21] - piCur[21] - iDeltaC );
    uiSum += abs( piOrg[22] - piCur[22] - iDeltaC );
    uiSum += abs( piOrg[23] - piCur[23] - iDeltaC );
    uiSum += abs( piOrg[24] - piCur[24] - iDeltaC );
    uiSum += abs( piOrg[25] - piCur[25] - iDeltaC );
    uiSum += abs( piOrg[26] - piCur[26] - iDeltaC );
    uiSum += abs( piOrg[27] - piCur[27] - iDeltaC );
    uiSum += abs( piOrg[28] - piCur[28] - iDeltaC );
    uiSum += abs( piOrg[29] - piCur[29] - iDeltaC );
    uiSum += abs( piOrg[30] - piCur[30] - iDeltaC );
    uiSum += abs( piOrg[31] - piCur[31] - iDeltaC );
    uiSum += abs( piOrg[32] - piCur[32] - iDeltaC );
    uiSum += abs( piOrg[33] - piCur[33] - iDeltaC );
    uiSum += abs( piOrg[34] - piCur[34] - iDeltaC );
    uiSum += abs( piOrg[35] - piCur[35] - iDeltaC );
    uiSum += abs( piOrg[36] - piCur[36] - iDeltaC );
    uiSum += abs( piOrg[37] - piCur[37] - iDeltaC );
    uiSum += abs( piOrg[38] - piCur[38] - iDeltaC );
    uiSum += abs( piOrg[39] - piCur[39] - iDeltaC );
    uiSum += abs( piOrg[40] - piCur[40] - iDeltaC );
    uiSum += abs( piOrg[41] - piCur[41] - iDeltaC );
    uiSum += abs( piOrg[42] - piCur[42] - iDeltaC );
    uiSum += abs( piOrg[43] - piCur[43] - iDeltaC );
    uiSum += abs( piOrg[44] - piCur[44] - iDeltaC );
    uiSum += abs( piOrg[45] - piCur[45] - iDeltaC );
    uiSum += abs( piOrg[46] - piCur[46] - iDeltaC );
    uiSum += abs( piOrg[47] - piCur[47] - iDeltaC );

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  uiSum <<= iSubShift;
  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}
#endif
// --------------------------------------------------------------------------------------------------------------------
// SSE
// --------------------------------------------------------------------------------------------------------------------

Distortion TComRdCost::xGetSSE( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int iTemp;

  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n++ )
    {
      iTemp = piOrg[n  ] - piCur[n  ];
      uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE4( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 4 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {

    iTemp = piOrg[0] - piCur[0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[1] - piCur[1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[2] - piCur[2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[3] - piCur[3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE8( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 8 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {
    iTemp = piOrg[0] - piCur[0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[1] - piCur[1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[2] - piCur[2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[3] - piCur[3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[4] - piCur[4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[5] - piCur[5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[6] - piCur[6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[7] - piCur[7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE16( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 16 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {

    iTemp = piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE16N( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n+=16 )
    {

      iTemp = piOrg[n+ 0] - piCur[n+ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 1] - piCur[n+ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 2] - piCur[n+ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 3] - piCur[n+ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 4] - piCur[n+ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 5] - piCur[n+ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 6] - piCur[n+ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 7] - piCur[n+ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 8] - piCur[n+ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+ 9] - piCur[n+ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+10] - piCur[n+10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+11] - piCur[n+11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+12] - piCur[n+12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+13] - piCur[n+13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+14] - piCur[n+14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      iTemp = piOrg[n+15] - piCur[n+15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE32( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 32 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {

    iTemp = piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[16] - piCur[16]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[17] - piCur[17]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[18] - piCur[18]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[19] - piCur[19]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[20] - piCur[20]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[21] - piCur[21]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[22] - piCur[22]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[23] - piCur[23]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[24] - piCur[24]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[25] - piCur[25]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[26] - piCur[26]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[27] - piCur[27]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[28] - piCur[28]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[29] - piCur[29]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[30] - piCur[30]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[31] - piCur[31]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

Distortion TComRdCost::xGetSSE64( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 64 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);

  Intermediate_Int  iTemp;

  for( ; iRows != 0; iRows-- )
  {
    iTemp = piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[16] - piCur[16]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[17] - piCur[17]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[18] - piCur[18]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[19] - piCur[19]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[20] - piCur[20]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[21] - piCur[21]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[22] - piCur[22]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[23] - piCur[23]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[24] - piCur[24]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[25] - piCur[25]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[26] - piCur[26]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[27] - piCur[27]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[28] - piCur[28]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[29] - piCur[29]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[30] - piCur[30]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[31] - piCur[31]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[32] - piCur[32]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[33] - piCur[33]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[34] - piCur[34]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[35] - piCur[35]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[36] - piCur[36]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[37] - piCur[37]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[38] - piCur[38]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[39] - piCur[39]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[40] - piCur[40]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[41] - piCur[41]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[42] - piCur[42]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[43] - piCur[43]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[44] - piCur[44]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[45] - piCur[45]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[46] - piCur[46]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[47] - piCur[47]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[48] - piCur[48]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[49] - piCur[49]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[50] - piCur[50]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[51] - piCur[51]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[52] - piCur[52]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[53] - piCur[53]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[54] - piCur[54]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[55] - piCur[55]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[56] - piCur[56]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[57] - piCur[57]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[58] - piCur[58]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[59] - piCur[59]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[60] - piCur[60]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[61] - piCur[61]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[62] - piCur[62]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    iTemp = piOrg[63] - piCur[63]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);

    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  return ( uiSum );
}

#if WCG_LUMA_DQP_CM_SCALE
Void TComRdCost::initLumaLevelToWeightTable(Bool isSDR) {
  for (Int i = 0; i < LUMA_LEVEL_TO_DQP_LUT_MAXSIZE; i++) {
    Double x = i;
    Double y;

    if (isSDR)  // set SDR weight table
    {
      y = 0.03*x - 3.0;        // this is the Equation used to derive the luma qp LUT for SDR in ST-2084
      y = y<0 ? 0 : (y>12 ? 12 : y);

    }
    else
    { // set SDR weight table
      y = 0.015*x - 1.5 - 6;   // this is the Equation used to derive the luma qp LUT for HDR in MPEG HDR anchor3.2 (JCTCX-X1020) 
      y = y<-3 ? -3 : (y>6 ? 6 : y);
    }
    
    m_lumaLevelToWeightPLUT[i] = pow(2.0, y / 3.0);      // or power(10, dQp/10)      they are almost equal       

  }
}

Distortion TComRdCost::getWeightedMSE(Int compIdx, Pel org, Pel cur, UInt       uiShift, Pel orgLuma) {
  Distortion distortionVal = 0;
  Intermediate_Int iTemp = org - cur;
  assert(org>=0);

  if (compIdx == COMPONENT_Y)
  {
     assert(org==orgLuma);
  }
  // use luma to get weight
  Double weight = m_lumaLevelToWeightPLUT[orgLuma];
  Intermediate_Int mse = Intermediate_Int(weight*(Double)iTemp*(Double)iTemp+0.5);
  distortionVal = Distortion( mse >> uiShift);
  return distortionVal;
}
Distortion TComRdCost::xGetSSE_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );  // ignore it for now
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420, could use getComponentScaleX, getComponentScaleY

  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n++ )
    {
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n  ], piCur[n  ], uiShift, piOrgLuma[n<<cShift]); 
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE4_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 4 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam ); // ignore it for now
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[0  ], piCur[0  ], uiShift, piOrgLuma[0<<cShift]);   // piOrg[0] - piCur[0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[1  ], piCur[1  ], uiShift, piOrgLuma[1<<cShift] );   // piOrg[1] - piCur[1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[2  ], piCur[2  ], uiShift, piOrgLuma[2<<cShift] );   // piOrg[2] - piCur[2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[3  ], piCur[3  ], uiShift, piOrgLuma[3<<cShift] );   // piOrg[3] - piCur[3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE8_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 8 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
 
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[0  ], piCur[0  ], uiShift, piOrgLuma[0  ]);   // piOrg[0] - piCur[0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[1  ], piCur[1  ], uiShift, piOrgLuma[1<<cShift  ]);  // piOrg[1] - piCur[1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[2  ], piCur[2  ], uiShift, piOrgLuma[2<<cShift  ]);  //piOrg[2] - piCur[2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[3  ], piCur[3  ], uiShift, piOrgLuma[3<<cShift  ]);  // piOrg[3] - piCur[3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[4  ], piCur[4  ], uiShift, piOrgLuma[4<<cShift  ]);  // piOrg[4] - piCur[4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[5  ], piCur[5  ], uiShift, piOrgLuma[5<<cShift  ]);  // piOrg[5] - piCur[5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[6  ], piCur[6  ], uiShift, piOrgLuma[6<<cShift  ]);  // piOrg[6] - piCur[6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[7  ], piCur[7  ], uiShift, piOrgLuma[7<<cShift  ]);  // piOrg[7] - piCur[7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE16_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 16 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
  
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[0  ], piCur[0  ], uiShift, piOrgLuma[0  ]);  // piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[1  ], piCur[1  ], uiShift, piOrgLuma[1<<cShift  ]);  //piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[2  ], piCur[2  ], uiShift, piOrgLuma[2<<cShift  ]);  //piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[3  ], piCur[3  ], uiShift, piOrgLuma[3<<cShift  ]);  //piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[4  ], piCur[4  ], uiShift, piOrgLuma[4<<cShift  ]);  //piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[5  ], piCur[5  ], uiShift, piOrgLuma[5<<cShift  ]);  //piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[6  ], piCur[6  ], uiShift, piOrgLuma[6<<cShift  ]);  //piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[7  ], piCur[7  ], uiShift, piOrgLuma[7<<cShift  ]);  //piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[8  ], piCur[8  ], uiShift, piOrgLuma[8<<cShift  ]);  //piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[9  ], piCur[9  ], uiShift, piOrgLuma[9<<cShift  ]);  //piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[10 ], piCur[10 ], uiShift, piOrgLuma[10<<cShift  ]);  //piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[11 ], piCur[11 ], uiShift, piOrgLuma[11<<cShift  ]);  //piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[12 ], piCur[12 ], uiShift, piOrgLuma[12<<cShift  ]);  //piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[13 ], piCur[13 ], uiShift, piOrgLuma[13<<cShift  ]);  //piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[14 ], piCur[14 ], uiShift, piOrgLuma[14<<cShift  ]);  //piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[15 ], piCur[15 ], uiShift, piOrgLuma[15<<cShift  ]);  //piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE16N_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    for (Int n = 0; n < iCols; n+=16 )
    {
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+0 ], piCur[n+0 ], uiShift, piOrgLuma[(n+0)<<cShift ]);  // iTemp = piOrg[n+ 0] - piCur[n+ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+1 ], piCur[n+1 ], uiShift, piOrgLuma[(n+1)<<cShift ]);  // iTemp = piOrg[n+ 1] - piCur[n+ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+2 ], piCur[n+2 ], uiShift, piOrgLuma[(n+2)<<cShift ]);  // iTemp = piOrg[n+ 2] - piCur[n+ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+3 ], piCur[n+3 ], uiShift, piOrgLuma[(n+3)<<cShift ]);  // iTemp = piOrg[n+ 3] - piCur[n+ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+4 ], piCur[n+4 ], uiShift, piOrgLuma[(n+4)<<cShift ]);  // iTemp = piOrg[n+ 4] - piCur[n+ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+5 ], piCur[n+5 ], uiShift, piOrgLuma[(n+5)<<cShift ]);  // iTemp = piOrg[n+ 5] - piCur[n+ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+6 ], piCur[n+6 ], uiShift, piOrgLuma[(n+6)<<cShift ]);  // iTemp = piOrg[n+ 6] - piCur[n+ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+7 ], piCur[n+7 ], uiShift, piOrgLuma[(n+7)<<cShift ]);  // iTemp = piOrg[n+ 7] - piCur[n+ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+8 ], piCur[n+8 ], uiShift, piOrgLuma[(n+8)<<cShift ]);  // iTemp = piOrg[n+ 8] - piCur[n+ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+9 ], piCur[n+9 ], uiShift, piOrgLuma[(n+9)<<cShift ]);  // iTemp = piOrg[n+ 9] - piCur[n+ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+10], piCur[n+10], uiShift, piOrgLuma[(n+10)<<cShift ]);  // iTemp = piOrg[n+10] - piCur[n+10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+11], piCur[n+11], uiShift, piOrgLuma[(n+11)<<cShift ]);  // iTemp = piOrg[n+11] - piCur[n+11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+12], piCur[n+12], uiShift, piOrgLuma[(n+12)<<cShift]);  // iTemp = piOrg[n+12] - piCur[n+12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+13], piCur[n+13], uiShift, piOrgLuma[(n+13)<<cShift ]);  // iTemp = piOrg[n+13] - piCur[n+13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+14], piCur[n+14], uiShift, piOrgLuma[(n+14)<<cShift ]);  // iTemp = piOrg[n+14] - piCur[n+14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
      uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[n+15], piCur[n+15], uiShift, piOrgLuma[(n+15)<<cShift ]);  // iTemp = piOrg[n+15] - piCur[n+15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE32_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 32 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
  
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[0 ], piCur[0 ], uiShift, piOrgLuma[0 ]);  // iTemp = piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[1 ], piCur[1 ], uiShift, piOrgLuma[1<<cShift ]);  // iTemp = piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[2 ], piCur[2 ], uiShift, piOrgLuma[2<<cShift ]);  // iTemp = piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[3 ], piCur[3 ], uiShift, piOrgLuma[3<<cShift ]);  // iTemp = piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[4 ], piCur[4 ], uiShift, piOrgLuma[4<<cShift ]);  // iTemp = piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[5 ], piCur[5 ], uiShift, piOrgLuma[5<<cShift ]);  // iTemp = piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[6 ], piCur[6 ], uiShift, piOrgLuma[6<<cShift ]);  // iTemp = piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[7 ], piCur[7 ], uiShift, piOrgLuma[7<<cShift ]);  // iTemp = piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[8 ], piCur[8 ], uiShift, piOrgLuma[8<<cShift ]);  // iTemp = piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[9 ], piCur[9 ], uiShift, piOrgLuma[9<<cShift ]);  // iTemp = piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[10], piCur[10], uiShift, piOrgLuma[10<<cShift ]);  // iTemp = piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[11], piCur[11], uiShift, piOrgLuma[11<<cShift ]);  // iTemp = piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[12], piCur[12], uiShift, piOrgLuma[12<<cShift ]);  // iTemp = piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[13], piCur[13], uiShift, piOrgLuma[13<<cShift ]);  // iTemp = piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[14], piCur[14], uiShift, piOrgLuma[14<<cShift ]);  // iTemp = piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[15], piCur[15], uiShift, piOrgLuma[15<<cShift ]);  // iTemp = piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[16], piCur[16], uiShift, piOrgLuma[16<<cShift ]);  //  iTemp = piOrg[16] - piCur[16]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[17], piCur[17], uiShift, piOrgLuma[17<<cShift ]);  //  iTemp = piOrg[17] - piCur[17]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[18], piCur[18], uiShift, piOrgLuma[18<<cShift ]);  //  iTemp = piOrg[18] - piCur[18]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[19], piCur[19], uiShift, piOrgLuma[19<<cShift ]);  //  iTemp = piOrg[19] - piCur[19]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[20], piCur[20], uiShift, piOrgLuma[20<<cShift ]);  //  iTemp = piOrg[20] - piCur[20]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[21], piCur[21], uiShift, piOrgLuma[21<<cShift ]);  //  iTemp = piOrg[21] - piCur[21]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[22], piCur[22], uiShift, piOrgLuma[22<<cShift ]);  //  iTemp = piOrg[22] - piCur[22]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[23], piCur[23], uiShift, piOrgLuma[23<<cShift ]);  //  iTemp = piOrg[23] - piCur[23]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[24], piCur[24], uiShift, piOrgLuma[24<<cShift ]);  //  iTemp = piOrg[24] - piCur[24]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[25], piCur[25], uiShift, piOrgLuma[25<<cShift ]);  //  iTemp = piOrg[25] - piCur[25]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[26], piCur[26], uiShift, piOrgLuma[26<<cShift ]);  //  iTemp = piOrg[26] - piCur[26]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[27], piCur[27], uiShift, piOrgLuma[27<<cShift ]);  //  iTemp = piOrg[27] - piCur[27]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[28], piCur[28], uiShift, piOrgLuma[28<<cShift ]);  //  iTemp = piOrg[28] - piCur[28]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[29], piCur[29], uiShift, piOrgLuma[29<<cShift ]);  //  iTemp = piOrg[29] - piCur[29]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[30], piCur[30], uiShift, piOrgLuma[30<<cShift ]);  //  iTemp = piOrg[30] - piCur[30]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[31], piCur[31], uiShift, piOrgLuma[31<<cShift ]);  //  iTemp = piOrg[31] - piCur[31]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
Distortion TComRdCost::xGetSSE64_WTD( DistParam* pcDtParam )
{
  if ( pcDtParam->bApplyWeight )
  {
    assert( pcDtParam->iCols == 64 );
    return TComRdCostWeightPrediction::xGetSSEw( pcDtParam );
  }
  const Pel* piOrg   = pcDtParam->pOrg;
  const Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStrideCur = pcDtParam->iStrideCur;
  const Pel* piOrgLuma   = pcDtParam->pOrgLuma;
  Int  iStrideOrgLuma = pcDtParam->iStrideOrgLuma;
  Int cShift = (pcDtParam->compIdx==COMPONENT_Y) ? 0 : 1; // assume 420
 
  Distortion uiSum   = 0;
  UInt       uiShift = DISTORTION_PRECISION_ADJUSTMENT((pcDtParam->bitDepth-8) << 1);
  for( ; iRows != 0; iRows-- )
  {
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[0 ], piCur[0 ], uiShift, piOrgLuma[0 ]);  // iTemp = piOrg[ 0] - piCur[ 0]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[1 ], piCur[1 ], uiShift, piOrgLuma[1<<cShift ]);  // iTemp = piOrg[ 1] - piCur[ 1]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[2 ], piCur[2 ], uiShift, piOrgLuma[2<<cShift ]);  // iTemp = piOrg[ 2] - piCur[ 2]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[3 ], piCur[3 ], uiShift, piOrgLuma[3<<cShift ]);  // iTemp = piOrg[ 3] - piCur[ 3]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[4 ], piCur[4 ], uiShift, piOrgLuma[4<<cShift ]);  // iTemp = piOrg[ 4] - piCur[ 4]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[5 ], piCur[5 ], uiShift, piOrgLuma[5<<cShift ]);  // iTemp = piOrg[ 5] - piCur[ 5]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[6 ], piCur[6 ], uiShift, piOrgLuma[6<<cShift ]);  // iTemp = piOrg[ 6] - piCur[ 6]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[7 ], piCur[7 ], uiShift, piOrgLuma[7<<cShift ]);  // iTemp = piOrg[ 7] - piCur[ 7]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[8 ], piCur[8 ], uiShift, piOrgLuma[8<<cShift ]);  // iTemp = piOrg[ 8] - piCur[ 8]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[9 ], piCur[9 ], uiShift, piOrgLuma[9<<cShift ]);  // iTemp = piOrg[ 9] - piCur[ 9]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[10], piCur[10], uiShift, piOrgLuma[10<<cShift]);  // iTemp = piOrg[10] - piCur[10]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[11], piCur[11], uiShift, piOrgLuma[11<<cShift]);  // iTemp = piOrg[11] - piCur[11]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[12], piCur[12], uiShift, piOrgLuma[12<<cShift]);  // iTemp = piOrg[12] - piCur[12]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[13], piCur[13], uiShift, piOrgLuma[13<<cShift]);  // iTemp = piOrg[13] - piCur[13]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[14], piCur[14], uiShift, piOrgLuma[14<<cShift]);  // iTemp = piOrg[14] - piCur[14]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[15], piCur[15], uiShift, piOrgLuma[15<<cShift]);  // iTemp = piOrg[15] - piCur[15]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[16], piCur[16], uiShift, piOrgLuma[16<<cShift]);  //  iTemp = piOrg[16] - piCur[16]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[17], piCur[17], uiShift, piOrgLuma[17<<cShift]);  //  iTemp = piOrg[17] - piCur[17]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[18], piCur[18], uiShift, piOrgLuma[18<<cShift]);  //  iTemp = piOrg[18] - piCur[18]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[19], piCur[19], uiShift, piOrgLuma[19<<cShift]);  //  iTemp = piOrg[19] - piCur[19]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[20], piCur[20], uiShift, piOrgLuma[20<<cShift]);  //  iTemp = piOrg[20] - piCur[20]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[21], piCur[21], uiShift, piOrgLuma[21<<cShift]);  //  iTemp = piOrg[21] - piCur[21]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[22], piCur[22], uiShift, piOrgLuma[22<<cShift]);  //  iTemp = piOrg[22] - piCur[22]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[23], piCur[23], uiShift, piOrgLuma[23<<cShift]);  //  iTemp = piOrg[23] - piCur[23]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[24], piCur[24], uiShift, piOrgLuma[24<<cShift]);  //  iTemp = piOrg[24] - piCur[24]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[25], piCur[25], uiShift, piOrgLuma[25<<cShift]);  //  iTemp = piOrg[25] - piCur[25]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[26], piCur[26], uiShift, piOrgLuma[26<<cShift]);  //  iTemp = piOrg[26] - piCur[26]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[27], piCur[27], uiShift, piOrgLuma[27<<cShift]);  //  iTemp = piOrg[27] - piCur[27]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[28], piCur[28], uiShift, piOrgLuma[28<<cShift]);  //  iTemp = piOrg[28] - piCur[28]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[29], piCur[29], uiShift, piOrgLuma[29<<cShift]);  //  iTemp = piOrg[29] - piCur[29]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[30], piCur[30], uiShift, piOrgLuma[30<<cShift]);  //  iTemp = piOrg[30] - piCur[30]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[31], piCur[31], uiShift, piOrgLuma[31<<cShift]);  //  iTemp = piOrg[31] - piCur[31]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[32], piCur[32], uiShift, piOrgLuma[32<<cShift]);  // iTemp = piOrg[32] - piCur[32]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[33], piCur[33], uiShift, piOrgLuma[33<<cShift]);  // iTemp = piOrg[33] - piCur[33]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[34], piCur[34], uiShift, piOrgLuma[34<<cShift]);  // iTemp = piOrg[34] - piCur[34]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[35], piCur[35], uiShift, piOrgLuma[35<<cShift]);  // iTemp = piOrg[35] - piCur[35]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[36], piCur[36], uiShift, piOrgLuma[36<<cShift]);  // iTemp = piOrg[36] - piCur[36]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[37], piCur[37], uiShift, piOrgLuma[37<<cShift]);  // iTemp = piOrg[37] - piCur[37]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[38], piCur[38], uiShift, piOrgLuma[38<<cShift]);  // iTemp = piOrg[38] - piCur[38]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[39], piCur[39], uiShift, piOrgLuma[39<<cShift]);  // iTemp = piOrg[39] - piCur[39]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[40], piCur[40], uiShift, piOrgLuma[40<<cShift]);  // iTemp = piOrg[40] - piCur[40]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[41], piCur[41], uiShift, piOrgLuma[41<<cShift]);  // iTemp = piOrg[41] - piCur[41]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[42], piCur[42], uiShift, piOrgLuma[42<<cShift]);  // iTemp = piOrg[42] - piCur[42]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[43], piCur[43], uiShift, piOrgLuma[43<<cShift]);  // iTemp = piOrg[43] - piCur[43]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[44], piCur[44], uiShift, piOrgLuma[44<<cShift]);  // iTemp = piOrg[44] - piCur[44]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[45], piCur[45], uiShift, piOrgLuma[45<<cShift]);  // iTemp = piOrg[45] - piCur[45]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[46], piCur[46], uiShift, piOrgLuma[46<<cShift]);  // iTemp = piOrg[46] - piCur[46]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[47], piCur[47], uiShift, piOrgLuma[47<<cShift]);  // iTemp = piOrg[47] - piCur[47]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[48], piCur[48], uiShift, piOrgLuma[48<<cShift]);  // iTemp = piOrg[48] - piCur[48]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[49], piCur[49], uiShift, piOrgLuma[49<<cShift]);  // iTemp = piOrg[49] - piCur[49]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[50], piCur[50], uiShift, piOrgLuma[50<<cShift]);  // iTemp = piOrg[50] - piCur[50]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[51], piCur[51], uiShift, piOrgLuma[51<<cShift]);  // iTemp = piOrg[51] - piCur[51]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[52], piCur[52], uiShift, piOrgLuma[52<<cShift]);  // iTemp = piOrg[52] - piCur[52]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[53], piCur[53], uiShift, piOrgLuma[53<<cShift]);  // iTemp = piOrg[53] - piCur[53]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[54], piCur[54], uiShift, piOrgLuma[54<<cShift]);  // iTemp = piOrg[54] - piCur[54]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[55], piCur[55], uiShift, piOrgLuma[55<<cShift]);  // iTemp = piOrg[55] - piCur[55]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[56], piCur[56], uiShift, piOrgLuma[56<<cShift]);  // iTemp = piOrg[56] - piCur[56]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[57], piCur[57], uiShift, piOrgLuma[57<<cShift]);  // iTemp = piOrg[57] - piCur[57]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[58], piCur[58], uiShift, piOrgLuma[58<<cShift]);  // iTemp = piOrg[58] - piCur[58]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[59], piCur[59], uiShift, piOrgLuma[59<<cShift]);  // iTemp = piOrg[59] - piCur[59]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[60], piCur[60], uiShift, piOrgLuma[60<<cShift]);  // iTemp = piOrg[60] - piCur[60]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[61], piCur[61], uiShift, piOrgLuma[61<<cShift]);  // iTemp = piOrg[61] - piCur[61]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[62], piCur[62], uiShift, piOrgLuma[62<<cShift]);  // iTemp = piOrg[62] - piCur[62]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    uiSum += getWeightedMSE(pcDtParam->compIdx, piOrg[63], piCur[63], uiShift, piOrgLuma[63<<cShift]);  // iTemp = piOrg[63] - piCur[63]; uiSum += Distortion(( iTemp * iTemp ) >> uiShift);
    piOrg += iStrideOrg;
    piCur += iStrideCur;
    piOrgLuma += iStrideOrgLuma<<cShift;
  }
  return ( uiSum );
}
#endif
// --------------------------------------------------------------------------------------------------------------------
// HADAMARD with step (used in fractional search)
// --------------------------------------------------------------------------------------------------------------------

Distortion TComRdCost::xCalcHADs2x2( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep )
{
  Distortion satd = 0;
  TCoeff diff[4], m[4];
  assert( iStep == 1 );
  diff[0] = piOrg[0             ] - piCur[0];
  diff[1] = piOrg[1             ] - piCur[1];
  diff[2] = piOrg[iStrideOrg    ] - piCur[0 + iStrideCur];
  diff[3] = piOrg[iStrideOrg + 1] - piCur[1 + iStrideCur];
  m[0] = diff[0] + diff[2];
  m[1] = diff[1] + diff[3];
  m[2] = diff[0] - diff[2];
  m[3] = diff[1] - diff[3];

  satd += abs(m[0] + m[1]);
  satd += abs(m[0] - m[1]);
  satd += abs(m[2] + m[3]);
  satd += abs(m[2] - m[3]);

  return satd;
}
#if ITH_C3_HSADCOMPUTER==C_FXN
Distortion TComRdCost::xCalcHADs4x4( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth)
{
  Int k;
  Distortion satd = 0;
  TCoeff diff[16], m[16], d[16];

  assert( iStep == 1 );
  for( k = 0; k < 16; k+=4 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  /*===== hadamard transform =====*/
  m[ 0] = diff[ 0] + diff[12];
  m[ 1] = diff[ 1] + diff[13];
  m[ 2] = diff[ 2] + diff[14];
  m[ 3] = diff[ 3] + diff[15];
  m[ 4] = diff[ 4] + diff[ 8];
  m[ 5] = diff[ 5] + diff[ 9];
  m[ 6] = diff[ 6] + diff[10];
  m[ 7] = diff[ 7] + diff[11];
  m[ 8] = diff[ 4] - diff[ 8];
  m[ 9] = diff[ 5] - diff[ 9];
  m[10] = diff[ 6] - diff[10];
  m[11] = diff[ 7] - diff[11];
  m[12] = diff[ 0] - diff[12];
  m[13] = diff[ 1] - diff[13];
  m[14] = diff[ 2] - diff[14];
  m[15] = diff[ 3] - diff[15];

  d[ 0] = m[ 0] + m[ 4];
  d[ 1] = m[ 1] + m[ 5];
  d[ 2] = m[ 2] + m[ 6];
  d[ 3] = m[ 3] + m[ 7];
  d[ 4] = m[ 8] + m[12];
  d[ 5] = m[ 9] + m[13];
  d[ 6] = m[10] + m[14];
  d[ 7] = m[11] + m[15];
  d[ 8] = m[ 0] - m[ 4];
  d[ 9] = m[ 1] - m[ 5];
  d[10] = m[ 2] - m[ 6];
  d[11] = m[ 3] - m[ 7];
  d[12] = m[12] - m[ 8];
  d[13] = m[13] - m[ 9];
  d[14] = m[14] - m[10];
  d[15] = m[15] - m[11];

  m[ 0] = d[ 0] + d[ 3];
  m[ 1] = d[ 1] + d[ 2];
  m[ 2] = d[ 1] - d[ 2];
  m[ 3] = d[ 0] - d[ 3];
  m[ 4] = d[ 4] + d[ 7];
  m[ 5] = d[ 5] + d[ 6];
  m[ 6] = d[ 5] - d[ 6];
  m[ 7] = d[ 4] - d[ 7];
  m[ 8] = d[ 8] + d[11];
  m[ 9] = d[ 9] + d[10];
  m[10] = d[ 9] - d[10];
  m[11] = d[ 8] - d[11];
  m[12] = d[12] + d[15];
  m[13] = d[13] + d[14];
  m[14] = d[13] - d[14];
  m[15] = d[12] - d[15];

  d[ 0] = m[ 0] + m[ 1];
  d[ 1] = m[ 0] - m[ 1];
  d[ 2] = m[ 2] + m[ 3];
  d[ 3] = m[ 3] - m[ 2];
  d[ 4] = m[ 4] + m[ 5];
  d[ 5] = m[ 4] - m[ 5];
  d[ 6] = m[ 6] + m[ 7];
  d[ 7] = m[ 7] - m[ 6];
  d[ 8] = m[ 8] + m[ 9];
  d[ 9] = m[ 8] - m[ 9];
  d[10] = m[10] + m[11];
  d[11] = m[11] - m[10];
  d[12] = m[12] + m[13];
  d[13] = m[12] - m[13];
  d[14] = m[14] + m[15];
  d[15] = m[15] - m[14];

  for (k=0; k<16; ++k)
  {
    satd += abs(d[k]);
  }
  satd = ((satd+1)>>1);

  return satd;
}

Distortion TComRdCost::xCalcHADs8x8( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep 
#if COM16_C806_SIMD_OPT
  , Int bitDepth
#endif
  )
{
#if COM16_C806_SIMD_OPT
  if( bitDepth <= 10 )
  {
    return( simdHADs8x8( piOrg , piCur , iStrideOrg , iStrideCur ) );
  }
#endif
  Int k, i, j, jj;
  Distortion sad = 0;
  TCoeff diff[64], m1[8][8], m2[8][8], m3[8][8];
  assert( iStep == 1 );
  for( k = 0; k < 64; k += 8 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];
    diff[k+4] = piOrg[4] - piCur[4];
    diff[k+5] = piOrg[5] - piCur[5];
    diff[k+6] = piOrg[6] - piCur[6];
    diff[k+7] = piOrg[7] - piCur[7];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 3;
    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i=0; i < 8; i++)
  {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      sad += abs(m2[i][j]);
    }
  }

  sad=((sad+2)>>2);

  return sad;
}

#if JVET_C0024_QTBT
Distortion TComRdCost::xCalcHADs16x8( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{   //need to add SIMD implementation ,JCA
  Int k, i, j, jj, sad=0;
  Int diff[128], m1[8][16], m2[8][16];
  for( k = 0; k < 128; k += 16 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];
    diff[k+4] = piOrg[4] - piCur[4];
    diff[k+5] = piOrg[5] - piCur[5];
    diff[k+6] = piOrg[6] - piCur[6];
    diff[k+7] = piOrg[7] - piCur[7];

    diff[k+8]  = piOrg[8]  - piCur[8] ;
    diff[k+9]  = piOrg[9]  - piCur[9] ;
    diff[k+10] = piOrg[10] - piCur[10];
    diff[k+11] = piOrg[11] - piCur[11];
    diff[k+12] = piOrg[12] - piCur[12];
    diff[k+13] = piOrg[13] - piCur[13];
    diff[k+14] = piOrg[14] - piCur[14];
    diff[k+15] = piOrg[15] - piCur[15];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 4;

    m2[j][0]  = diff[jj  ] + diff[jj+8];
    m2[j][1]  = diff[jj+1] + diff[jj+9];
    m2[j][2]  = diff[jj+2] + diff[jj+10];
    m2[j][3]  = diff[jj+3] + diff[jj+11];
    m2[j][4]  = diff[jj+4] + diff[jj+12];
    m2[j][5]  = diff[jj+5] + diff[jj+13];
    m2[j][6]  = diff[jj+6] + diff[jj+14];
    m2[j][7]  = diff[jj+7] + diff[jj+15];
    m2[j][8]  = diff[jj  ] - diff[jj+8];
    m2[j][9]  = diff[jj+1] - diff[jj+9];
    m2[j][10] = diff[jj+2] - diff[jj+10];
    m2[j][11] = diff[jj+3] - diff[jj+11];
    m2[j][12] = diff[jj+4] - diff[jj+12];
    m2[j][13] = diff[jj+5] - diff[jj+13];
    m2[j][14] = diff[jj+6] - diff[jj+14];
    m2[j][15] = diff[jj+7] - diff[jj+15];

    m1[j][0]  = m2[j][0]  + m2[j][4];
    m1[j][1]  = m2[j][1]  + m2[j][5];
    m1[j][2]  = m2[j][2]  + m2[j][6];
    m1[j][3]  = m2[j][3]  + m2[j][7];
    m1[j][4]  = m2[j][0]  - m2[j][4];
    m1[j][5]  = m2[j][1]  - m2[j][5];
    m1[j][6]  = m2[j][2]  - m2[j][6];
    m1[j][7]  = m2[j][3]  - m2[j][7];
    m1[j][8]  = m2[j][8]  + m2[j][12];
    m1[j][9]  = m2[j][9]  + m2[j][13];
    m1[j][10] = m2[j][10] + m2[j][14];
    m1[j][11] = m2[j][11] + m2[j][15];
    m1[j][12] = m2[j][8]  - m2[j][12];
    m1[j][13] = m2[j][9]  - m2[j][13];
    m1[j][14] = m2[j][10] - m2[j][14];
    m1[j][15] = m2[j][11] - m2[j][15];

    m2[j][0]  = m1[j][0]  + m1[j][2];
    m2[j][1]  = m1[j][1]  + m1[j][3];
    m2[j][2]  = m1[j][0]  - m1[j][2];
    m2[j][3]  = m1[j][1]  - m1[j][3];
    m2[j][4]  = m1[j][4]  + m1[j][6];
    m2[j][5]  = m1[j][5]  + m1[j][7];
    m2[j][6]  = m1[j][4]  - m1[j][6];
    m2[j][7]  = m1[j][5]  - m1[j][7];
    m2[j][8]  = m1[j][8]  + m1[j][10];
    m2[j][9]  = m1[j][9]  + m1[j][11];
    m2[j][10] = m1[j][8]  - m1[j][10];
    m2[j][11] = m1[j][9]  - m1[j][11];
    m2[j][12] = m1[j][12] + m1[j][14];
    m2[j][13] = m1[j][13] + m1[j][15];
    m2[j][14] = m1[j][12] - m1[j][14];
    m2[j][15] = m1[j][13] - m1[j][15];

    m1[j][0]  = m2[j][0]  + m2[j][1];
    m1[j][1]  = m2[j][0]  - m2[j][1];
    m1[j][2]  = m2[j][2]  + m2[j][3];
    m1[j][3]  = m2[j][2]  - m2[j][3];
    m1[j][4]  = m2[j][4]  + m2[j][5];
    m1[j][5]  = m2[j][4]  - m2[j][5];
    m1[j][6]  = m2[j][6]  + m2[j][7];
    m1[j][7]  = m2[j][6]  - m2[j][7];
    m1[j][8]  = m2[j][8]  + m2[j][9];
    m1[j][9]  = m2[j][8]  - m2[j][9];
    m1[j][10] = m2[j][10] + m2[j][11];
    m1[j][11] = m2[j][10] - m2[j][11];
    m1[j][12] = m2[j][12] + m2[j][13];
    m1[j][13] = m2[j][12] - m2[j][13];
    m1[j][14] = m2[j][14] + m2[j][15];
    m1[j][15] = m2[j][14] - m2[j][15];
  }

  //vertical
  for (i=0; i < 16; i++)
  {    
    m2[0][i] = m1[0][i] + m1[4][i];
    m2[1][i] = m1[1][i] + m1[5][i];
    m2[2][i] = m1[2][i] + m1[6][i];
    m2[3][i] = m1[3][i] + m1[7][i];
    m2[4][i] = m1[0][i] - m1[4][i];
    m2[5][i] = m1[1][i] - m1[5][i];
    m2[6][i] = m1[2][i] - m1[6][i];
    m2[7][i] = m1[3][i] - m1[7][i];

    m1[0][i] = m2[0][i] + m2[2][i];
    m1[1][i] = m2[1][i] + m2[3][i];
    m1[2][i] = m2[0][i] - m2[2][i];
    m1[3][i] = m2[1][i] - m2[3][i];
    m1[4][i] = m2[4][i] + m2[6][i];
    m1[5][i] = m2[5][i] + m2[7][i];
    m1[6][i] = m2[4][i] - m2[6][i];
    m1[7][i] = m2[5][i] - m2[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 16; j++)
    {
      sad += abs(m2[i][j]);
    }
  }

  sad=(Int)(sad/sqrt(16.0*8)*2);

  return sad;
}

Distortion TComRdCost::xCalcHADs8x16( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
  Int k, i, j, jj, sad=0;
  Int diff[128], m1[16][8], m2[16][8];
  for( k = 0; k < 128; k += 8 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];
    diff[k+4] = piOrg[4] - piCur[4];
    diff[k+5] = piOrg[5] - piCur[5];
    diff[k+6] = piOrg[6] - piCur[6];
    diff[k+7] = piOrg[7] - piCur[7];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 16; j++)
  {
    jj = j << 3;

    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i=0; i < 8; i++)
  {    
    m1[0][i]  = m2[0][i] + m2[8][i];
    m1[1][i]  = m2[1][i] + m2[9][i];
    m1[2][i]  = m2[2][i] + m2[10][i];
    m1[3][i]  = m2[3][i] + m2[11][i];
    m1[4][i]  = m2[4][i] + m2[12][i];
    m1[5][i]  = m2[5][i] + m2[13][i];
    m1[6][i]  = m2[6][i] + m2[14][i];
    m1[7][i]  = m2[7][i] + m2[15][i];
    m1[8][i]  = m2[0][i] - m2[8][i];
    m1[9][i]  = m2[1][i] - m2[9][i];
    m1[10][i]  = m2[2][i] - m2[10][i];
    m1[11][i]  = m2[3][i] - m2[11][i];
    m1[12][i]  = m2[4][i] - m2[12][i];
    m1[13][i]  = m2[5][i] - m2[13][i];
    m1[14][i]  = m2[6][i] - m2[14][i];
    m1[15][i]  = m2[7][i] - m2[15][i];

    m2[0][i]  = m1[0][i]  + m1[4][i];
    m2[1][i]  = m1[1][i]  + m1[5][i];
    m2[2][i]  = m1[2][i]  + m1[6][i];
    m2[3][i]  = m1[3][i]  + m1[7][i];
    m2[4][i]  = m1[0][i]  - m1[4][i];
    m2[5][i]  = m1[1][i]  - m1[5][i];
    m2[6][i]  = m1[2][i]  - m1[6][i];
    m2[7][i]  = m1[3][i]  - m1[7][i];
    m2[8][i]  = m1[8][i]  + m1[12][i];
    m2[9][i]  = m1[9][i]  + m1[13][i];
    m2[10][i] = m1[10][i] + m1[14][i];
    m2[11][i] = m1[11][i] + m1[15][i];
    m2[12][i] = m1[8][i]  - m1[12][i];
    m2[13][i] = m1[9][i]  - m1[13][i];
    m2[14][i] = m1[10][i] - m1[14][i];
    m2[15][i] = m1[11][i] - m1[15][i];

    m1[0][i]  = m2[0][i]  + m2[2][i];
    m1[1][i]  = m2[1][i]  + m2[3][i];
    m1[2][i]  = m2[0][i]  - m2[2][i];
    m1[3][i]  = m2[1][i]  - m2[3][i];
    m1[4][i]  = m2[4][i]  + m2[6][i];
    m1[5][i]  = m2[5][i]  + m2[7][i];
    m1[6][i]  = m2[4][i]  - m2[6][i];
    m1[7][i]  = m2[5][i]  - m2[7][i];
    m1[8][i]  = m2[8][i]  + m2[10][i];
    m1[9][i]  = m2[9][i]  + m2[11][i];
    m1[10][i] = m2[8][i]  - m2[10][i];
    m1[11][i] = m2[9][i]  - m2[11][i];
    m1[12][i] = m2[12][i] + m2[14][i];
    m1[13][i] = m2[13][i] + m2[15][i];
    m1[14][i] = m2[12][i] - m2[14][i];
    m1[15][i] = m2[13][i] - m2[15][i];

    m2[0][i]  = m1[0][i]  + m1[1][i];
    m2[1][i]  = m1[0][i]  - m1[1][i];
    m2[2][i]  = m1[2][i]  + m1[3][i];
    m2[3][i]  = m1[2][i]  - m1[3][i];
    m2[4][i]  = m1[4][i]  + m1[5][i];
    m2[5][i]  = m1[4][i]  - m1[5][i];
    m2[6][i]  = m1[6][i]  + m1[7][i];
    m2[7][i]  = m1[6][i]  - m1[7][i];
    m2[8][i]  = m1[8][i]  + m1[9][i];
    m2[9][i]  = m1[8][i]  - m1[9][i];
    m2[10][i] = m1[10][i] + m1[11][i];
    m2[11][i] = m1[10][i] - m1[11][i];
    m2[12][i] = m1[12][i] + m1[13][i];
    m2[13][i] = m1[12][i] - m1[13][i];
    m2[14][i] = m1[14][i] + m1[15][i];
    m2[15][i] = m1[14][i] - m1[15][i];
  }

  for (i = 0; i < 16; i++)
  {
    for (j = 0; j < 8; j++)
    {
      sad += abs(m2[i][j]);
    }
  }

  sad=(Int)(sad/sqrt(16.0*8)*2);

  return sad;
}
Distortion TComRdCost::xCalcHADs4x8( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
  Int k, i, j, jj, sad=0;
  Int diff[32], m1[8][4], m2[8][4];
  for( k = 0; k < 32; k += 4 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 2;
    m2[j][0] = diff[jj  ] + diff[jj+2];
    m2[j][1] = diff[jj+1] + diff[jj+3];
    m2[j][2] = diff[jj  ] - diff[jj+2];
    m2[j][3] = diff[jj+1] - diff[jj+3];

    m1[j][0] = m2[j][0] + m2[j][1];
    m1[j][1] = m2[j][0] - m2[j][1];
    m1[j][2] = m2[j][2] + m2[j][3];
    m1[j][3] = m2[j][2] - m2[j][3];
  }

  //vertical
  for (i=0; i<4; i++)
  {
    m2[0][i] = m1[0][i] + m1[4][i];
    m2[1][i] = m1[1][i] + m1[5][i];
    m2[2][i] = m1[2][i] + m1[6][i];
    m2[3][i] = m1[3][i] + m1[7][i];
    m2[4][i] = m1[0][i] - m1[4][i];
    m2[5][i] = m1[1][i] - m1[5][i];
    m2[6][i] = m1[2][i] - m1[6][i];
    m2[7][i] = m1[3][i] - m1[7][i];

    m1[0][i] = m2[0][i] + m2[2][i];
    m1[1][i] = m2[1][i] + m2[3][i];
    m1[2][i] = m2[0][i] - m2[2][i];
    m1[3][i] = m2[1][i] - m2[3][i];
    m1[4][i] = m2[4][i] + m2[6][i];
    m1[5][i] = m2[5][i] + m2[7][i];
    m1[6][i] = m2[4][i] - m2[6][i];
    m1[7][i] = m2[5][i] - m2[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 4; j++)
    {
      sad += abs(m2[i][j]);
    }
  }

  sad=(Int)(sad/sqrt(4.0*8)*2);

  return sad;
}

Distortion TComRdCost::xCalcHADs8x4( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
  Int k, i, j, jj, sad=0;
  Int diff[32], m1[4][8], m2[4][8];
  for( k = 0; k < 32; k += 8 )
  {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];
    diff[k+4] = piOrg[4] - piCur[4];
    diff[k+5] = piOrg[5] - piCur[5];
    diff[k+6] = piOrg[6] - piCur[6];
    diff[k+7] = piOrg[7] - piCur[7];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 4; j++)
  {
    jj = j << 3;

    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i=0; i < 8; i++)
  {    
    m1[0][i] = m2[0][i] + m2[2][i];
    m1[1][i] = m2[1][i] + m2[3][i];
    m1[2][i] = m2[0][i] - m2[2][i];
    m1[3][i] = m2[1][i] - m2[3][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
  }

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 8; j++)
    {
      sad += abs(m2[i][j]);
    }
  }

  sad=(Int)(sad/sqrt(4.0*8)*2);

  return sad;
}

#endif
#elif ITH_C3_HSADCOMPUTER==SIMD_FXN
Distortion TComRdCost::xCalcHADs4x4(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src1_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src3_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;

		pred0_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred1_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred3_8x16b = _mm_loadl_epi64((__m128i *) piCur);

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		/**********************Residue Calculation********************************/

		/************************* 4x4 Vertical Transform*************************/
		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);

		/* r0 + r1 + r2 + r3 */
		out0_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 + r1 - r2 - r3 */
		out2_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);

		/* r0 - r1 + r2 - r3 */
		out1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 - r1 - r2 + r3 */
		out3_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/************************* 4x4 Vertical Transform*************************/

		/**************** 4x4 horizontal transform *******************************/
		/***********************    4x4 16 bit Transpose  ************************/
		src3_8x16b = _mm_unpacklo_epi16(out0_8x16b, out1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(out2_8x16b, out3_8x16b);

		out0_8x16b = _mm_unpacklo_epi32(src3_8x16b, pred0_8x16b);
		out1_8x16b = _mm_srli_si128(out0_8x16b, 8);
		out2_8x16b = _mm_unpackhi_epi32(src3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_srli_si128(out2_8x16b, 8);
		/**********************   4x4 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);

		/* r0 + r1 + r2 + r3 */
		src0_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 + r1 - r2 - r3 */
		src2_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);

		/* r0 - r1 + r2 - r3 */
		src1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 - r1 - r2 + r3 */
		src3_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/**************** 4x4 horizontal transform *******************************/

		/****************************SATD calculation ****************************/
		src0_8x16b = _mm_abs_epi16(src0_8x16b);
		src1_8x16b = _mm_abs_epi16(src1_8x16b);
		src2_8x16b = _mm_abs_epi16(src2_8x16b);
		src3_8x16b = _mm_abs_epi16(src3_8x16b);

		src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
		src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
		src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
		src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);

		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

		sad += _mm_cvtsi128_si32(src0_8x16b);

		sad = ((sad + 1) >> 1);

		return sad;
	}
	else
	{
		Int k;
		Distortion satd = 0;
		TCoeff diff[16], m[16], d[16];

		assert(iStep == 1);
		for (k = 0; k < 16; k += 4)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		/*===== hadamard transform =====*/
		m[0] = diff[0] + diff[12];
		m[1] = diff[1] + diff[13];
		m[2] = diff[2] + diff[14];
		m[3] = diff[3] + diff[15];
		m[4] = diff[4] + diff[8];
		m[5] = diff[5] + diff[9];
		m[6] = diff[6] + diff[10];
		m[7] = diff[7] + diff[11];
		m[8] = diff[4] - diff[8];
		m[9] = diff[5] - diff[9];
		m[10] = diff[6] - diff[10];
		m[11] = diff[7] - diff[11];
		m[12] = diff[0] - diff[12];
		m[13] = diff[1] - diff[13];
		m[14] = diff[2] - diff[14];
		m[15] = diff[3] - diff[15];

		d[0] = m[0] + m[4];
		d[1] = m[1] + m[5];
		d[2] = m[2] + m[6];
		d[3] = m[3] + m[7];
		d[4] = m[8] + m[12];
		d[5] = m[9] + m[13];
		d[6] = m[10] + m[14];
		d[7] = m[11] + m[15];
		d[8] = m[0] - m[4];
		d[9] = m[1] - m[5];
		d[10] = m[2] - m[6];
		d[11] = m[3] - m[7];
		d[12] = m[12] - m[8];
		d[13] = m[13] - m[9];
		d[14] = m[14] - m[10];
		d[15] = m[15] - m[11];

		m[0] = d[0] + d[3];
		m[1] = d[1] + d[2];
		m[2] = d[1] - d[2];
		m[3] = d[0] - d[3];
		m[4] = d[4] + d[7];
		m[5] = d[5] + d[6];
		m[6] = d[5] - d[6];
		m[7] = d[4] - d[7];
		m[8] = d[8] + d[11];
		m[9] = d[9] + d[10];
		m[10] = d[9] - d[10];
		m[11] = d[8] - d[11];
		m[12] = d[12] + d[15];
		m[13] = d[13] + d[14];
		m[14] = d[13] - d[14];
		m[15] = d[12] - d[15];

		d[0] = m[0] + m[1];
		d[1] = m[0] - m[1];
		d[2] = m[2] + m[3];
		d[3] = m[3] - m[2];
		d[4] = m[4] + m[5];
		d[5] = m[4] - m[5];
		d[6] = m[6] + m[7];
		d[7] = m[7] - m[6];
		d[8] = m[8] + m[9];
		d[9] = m[8] - m[9];
		d[10] = m[10] + m[11];
		d[11] = m[11] - m[10];
		d[12] = m[12] + m[13];
		d[13] = m[12] - m[13];
		d[14] = m[14] + m[15];
		d[15] = m[15] - m[14];

		for (k = 0; k < 16; ++k)
		{
			satd += abs(d[k]);
		}
		satd = ((satd + 1) >> 1);

		return satd;
	}
}

Distortion TComRdCost::xCalcHADs8x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i src4_8x16b, src5_8x16b, src6_8x16b, src7_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i pred4_8x16b, pred5_8x16b, pred6_8x16b, pred7_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;
		__m128i out4_8x16b, out5_8x16b, out6_8x16b, out7_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src1_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src3_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src4_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src5_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src6_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src7_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;

		pred0_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred1_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred3_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred4_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred5_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred6_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred7_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_sub_epi16(src4_8x16b, pred4_8x16b);
		src5_8x16b = _mm_sub_epi16(src5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_sub_epi16(src6_8x16b, pred6_8x16b);
		src7_8x16b = _mm_sub_epi16(src7_8x16b, pred7_8x16b);
		/**********************Residue Calculation********************************/

		/**************** 8x8 horizontal transform *******************************/
		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		out0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		out1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		out2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		out3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		out4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		out5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		out6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		out7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		src5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		src7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/
		/**************** 8x8 horizontal transform *******************************/

		if (bitDepth == 8)
		{
			/************************* 8x8 Vertical Transform*************************/
			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi16(src4_8x16b, src5_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi16(src6_8x16b, src7_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi16(src4_8x16b, src5_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi16(src6_8x16b, src7_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi16(out0_8x16b);
			src1_8x16b = _mm_abs_epi16(out1_8x16b);
			src2_8x16b = _mm_abs_epi16(out2_8x16b);
			src3_8x16b = _mm_abs_epi16(out3_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out4_8x16b);
			src1_8x16b = _mm_abs_epi16(out5_8x16b);
			src2_8x16b = _mm_abs_epi16(out6_8x16b);
			src3_8x16b = _mm_abs_epi16(out7_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (sad + 2) >> 2;

			return sad;
		}
		else
		{
			__m128i out0a_8x16b, out1a_8x16b, out2a_8x16b, out3a_8x16b;
			__m128i out4a_8x16b, out5a_8x16b, out6a_8x16b, out7a_8x16b;
			__m128i tmp0_8x16b, tmp1_8x16b, tmp2_8x16b, tmp3_8x16b;
			__m128i tmp4_8x16b, tmp5_8x16b, tmp6_8x16b, tmp7_8x16b;

			/************************* 8x8 Vertical Transform*************************/
			tmp0_8x16b = _mm_srli_si128(src0_8x16b, 8);
			tmp1_8x16b = _mm_srli_si128(src1_8x16b, 8);
			tmp2_8x16b = _mm_srli_si128(src2_8x16b, 8);
			tmp3_8x16b = _mm_srli_si128(src3_8x16b, 8);
			tmp4_8x16b = _mm_srli_si128(src4_8x16b, 8);
			tmp5_8x16b = _mm_srli_si128(src5_8x16b, 8);
			tmp6_8x16b = _mm_srli_si128(src6_8x16b, 8);
			tmp7_8x16b = _mm_srli_si128(src7_8x16b, 8);

			/*************************First 4 pixels ********************************/
			src0_8x16b = _mm_cvtepi16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepi16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(src7_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out4_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out2_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out6_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src1_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src2_8x16b, src3_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src4_8x16b, src5_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src6_8x16b, src7_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out1_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out5_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out3_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out7_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/*************************First 4 pixels ********************************/

			/**************************Next 4 pixels *******************************/
			src0_8x16b = _mm_cvtepi16_epi32(tmp0_8x16b);
			src1_8x16b = _mm_cvtepi16_epi32(tmp1_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(tmp2_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(tmp3_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(tmp4_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(tmp5_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(tmp6_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(tmp7_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out4a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out2a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out6a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src1_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src2_8x16b, src3_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src4_8x16b, src5_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src6_8x16b, src7_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out1a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out5a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out3a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out7a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/**************************Next 4 pixels *******************************/
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi32(out0_8x16b);
			src1_8x16b = _mm_abs_epi32(out1_8x16b);
			src2_8x16b = _mm_abs_epi32(out2_8x16b);
			src3_8x16b = _mm_abs_epi32(out3_8x16b);
			src4_8x16b = _mm_abs_epi32(out4_8x16b);
			src5_8x16b = _mm_abs_epi32(out5_8x16b);
			src6_8x16b = _mm_abs_epi32(out6_8x16b);
			src7_8x16b = _mm_abs_epi32(out7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out0a_8x16b);
			src1_8x16b = _mm_abs_epi32(out1a_8x16b);
			src2_8x16b = _mm_abs_epi32(out2a_8x16b);
			src3_8x16b = _mm_abs_epi32(out3a_8x16b);
			src4_8x16b = _mm_abs_epi32(out4a_8x16b);
			src5_8x16b = _mm_abs_epi32(out5a_8x16b);
			src6_8x16b = _mm_abs_epi32(out6a_8x16b);
			src7_8x16b = _mm_abs_epi32(out7a_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (sad + 2) >> 2;

			return sad;
		}
	}
	else
	{
		Int k, i, j, jj;
		Distortion sad = 0;
		TCoeff diff[64], m1[8][8], m2[8][8], m3[8][8];
		assert(iStep == 1);
		for (k = 0; k < 64; k += 8)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];
			diff[k + 4] = piOrg[4] - piCur[4];
			diff[k + 5] = piOrg[5] - piCur[5];
			diff[k + 6] = piOrg[6] - piCur[6];
			diff[k + 7] = piOrg[7] - piCur[7];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		//horizontal
		for (j = 0; j < 8; j++)
		{
			jj = j << 3;
			m2[j][0] = diff[jj] + diff[jj + 4];
			m2[j][1] = diff[jj + 1] + diff[jj + 5];
			m2[j][2] = diff[jj + 2] + diff[jj + 6];
			m2[j][3] = diff[jj + 3] + diff[jj + 7];
			m2[j][4] = diff[jj] - diff[jj + 4];
			m2[j][5] = diff[jj + 1] - diff[jj + 5];
			m2[j][6] = diff[jj + 2] - diff[jj + 6];
			m2[j][7] = diff[jj + 3] - diff[jj + 7];

			m1[j][0] = m2[j][0] + m2[j][2];
			m1[j][1] = m2[j][1] + m2[j][3];
			m1[j][2] = m2[j][0] - m2[j][2];
			m1[j][3] = m2[j][1] - m2[j][3];
			m1[j][4] = m2[j][4] + m2[j][6];
			m1[j][5] = m2[j][5] + m2[j][7];
			m1[j][6] = m2[j][4] - m2[j][6];
			m1[j][7] = m2[j][5] - m2[j][7];

			m2[j][0] = m1[j][0] + m1[j][1];
			m2[j][1] = m1[j][0] - m1[j][1];
			m2[j][2] = m1[j][2] + m1[j][3];
			m2[j][3] = m1[j][2] - m1[j][3];
			m2[j][4] = m1[j][4] + m1[j][5];
			m2[j][5] = m1[j][4] - m1[j][5];
			m2[j][6] = m1[j][6] + m1[j][7];
			m2[j][7] = m1[j][6] - m1[j][7];
		}

		//vertical
		for (i = 0; i < 8; i++)
		{
			m3[0][i] = m2[0][i] + m2[4][i];
			m3[1][i] = m2[1][i] + m2[5][i];
			m3[2][i] = m2[2][i] + m2[6][i];
			m3[3][i] = m2[3][i] + m2[7][i];
			m3[4][i] = m2[0][i] - m2[4][i];
			m3[5][i] = m2[1][i] - m2[5][i];
			m3[6][i] = m2[2][i] - m2[6][i];
			m3[7][i] = m2[3][i] - m2[7][i];

			m1[0][i] = m3[0][i] + m3[2][i];
			m1[1][i] = m3[1][i] + m3[3][i];
			m1[2][i] = m3[0][i] - m3[2][i];
			m1[3][i] = m3[1][i] - m3[3][i];
			m1[4][i] = m3[4][i] + m3[6][i];
			m1[5][i] = m3[5][i] + m3[7][i];
			m1[6][i] = m3[4][i] - m3[6][i];
			m1[7][i] = m3[5][i] - m3[7][i];

			m2[0][i] = m1[0][i] + m1[1][i];
			m2[1][i] = m1[0][i] - m1[1][i];
			m2[2][i] = m1[2][i] + m1[3][i];
			m2[3][i] = m1[2][i] - m1[3][i];
			m2[4][i] = m1[4][i] + m1[5][i];
			m2[5][i] = m1[4][i] - m1[5][i];
			m2[6][i] = m1[6][i] + m1[7][i];
			m2[7][i] = m1[6][i] - m1[7][i];
		}

		for (i = 0; i < 8; i++)
		{
			for (j = 0; j < 8; j++)
			{
				sad += abs(m2[i][j]);
			}
		}

		sad = ((sad + 2) >> 2);

		return sad;
	}
}

#if JVET_C0024_QTBT
#define SAD_DIV_16x8 sqrt(16.0 * 8) * 2
Distortion TComRdCost::xCalcHADs16x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i src4_8x16b, src5_8x16b, src6_8x16b, src7_8x16b;
		__m128i src8_8x16b, src9_8x16b, src10_8x16b, src11_8x16b;
		__m128i src12_8x16b, src13_8x16b, src14_8x16b, src15_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i pred4_8x16b, pred5_8x16b, pred6_8x16b, pred7_8x16b;
		__m128i pred8_8x16b, pred9_8x16b, pred10_8x16b, pred11_8x16b;
		__m128i pred12_8x16b, pred13_8x16b, pred14_8x16b, pred15_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;
		__m128i out4_8x16b, out5_8x16b, out6_8x16b, out7_8x16b;
		__m128i out8_8x16b, out9_8x16b, out10_8x16b, out11_8x16b;
		__m128i out12_8x16b, out13_8x16b, out14_8x16b, out15_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src1_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src3_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src4_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src5_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src6_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src7_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;

		pred0_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred1_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred3_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred4_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred5_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred6_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred7_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_sub_epi16(src4_8x16b, pred4_8x16b);
		src5_8x16b = _mm_sub_epi16(src5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_sub_epi16(src6_8x16b, pred6_8x16b);
		src7_8x16b = _mm_sub_epi16(src7_8x16b, pred7_8x16b);

		src8_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src9_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src10_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src11_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src12_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src13_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;
		src14_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		src15_8x16b = _mm_loadu_si128((__m128i *)(piOrg + 8));
		piOrg = piOrg + iStrideOrg;

		pred8_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred9_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred10_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred11_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred12_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred13_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;
		pred14_8x16b = _mm_loadu_si128((__m128i *) piCur);
		pred15_8x16b = _mm_loadu_si128((__m128i *)(piCur + 8));
		piCur = piCur + iStrideCur;

		src8_8x16b = _mm_sub_epi16(src8_8x16b, pred8_8x16b);
		src9_8x16b = _mm_sub_epi16(src9_8x16b, pred9_8x16b);
		src10_8x16b = _mm_sub_epi16(src10_8x16b, pred10_8x16b);
		src11_8x16b = _mm_sub_epi16(src11_8x16b, pred11_8x16b);
		src12_8x16b = _mm_sub_epi16(src12_8x16b, pred12_8x16b);
		src13_8x16b = _mm_sub_epi16(src13_8x16b, pred13_8x16b);
		src14_8x16b = _mm_sub_epi16(src14_8x16b, pred14_8x16b);
		src15_8x16b = _mm_sub_epi16(src15_8x16b, pred15_8x16b);
		/**********************Residue Calculation********************************/

		/**************** 8x8 horizontal transform *******************************/
		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		out0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		out1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		out2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		out3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		out4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		out5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		out6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		out7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		src5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		src7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/
		/**************** 8x8 horizontal transform *******************************/

		/**************** 8x8 horizontal transform *******************************/
		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src8_8x16b, src9_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src10_8x16b, src11_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src12_8x16b, src13_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src14_8x16b, src15_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src8_8x16b, src9_8x16b);
		src10_8x16b = _mm_unpackhi_epi16(src10_8x16b, src11_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src12_8x16b, src13_8x16b);
		src14_8x16b = _mm_unpackhi_epi16(src14_8x16b, src15_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src10_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src10_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src14_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src14_8x16b);

		out0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		out1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		out2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		out3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		out4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		out5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		out6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		out7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src8_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src12_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src10_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src14_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src9_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src13_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src11_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src15_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src8_8x16b, src9_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src10_8x16b, src11_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src12_8x16b, src13_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src14_8x16b, src15_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src8_8x16b, src9_8x16b);
		src10_8x16b = _mm_unpackhi_epi16(src10_8x16b, src11_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src12_8x16b, src13_8x16b);
		src14_8x16b = _mm_unpackhi_epi16(src14_8x16b, src15_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src10_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src10_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src14_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src14_8x16b);

		src8_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src9_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src10_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src11_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		src12_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		src13_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		src14_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		src15_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/
		/**************** 8x8 horizontal transform *******************************/

		/****************Horizontal Transform Addition****************************/
		out0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
		out1_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);

		out2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);
		out3_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);

		out4_8x16b = _mm_add_epi16(src4_8x16b, src5_8x16b);
		out5_8x16b = _mm_sub_epi16(src4_8x16b, src5_8x16b);

		out6_8x16b = _mm_add_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_sub_epi16(src6_8x16b, src7_8x16b);

		out8_8x16b = _mm_add_epi16(src8_8x16b, src9_8x16b);
		out9_8x16b = _mm_sub_epi16(src8_8x16b, src9_8x16b);

		out10_8x16b = _mm_add_epi16(src10_8x16b, src11_8x16b);
		out11_8x16b = _mm_sub_epi16(src10_8x16b, src11_8x16b);

		out12_8x16b = _mm_add_epi16(src12_8x16b, src13_8x16b);
		out13_8x16b = _mm_sub_epi16(src12_8x16b, src13_8x16b);

		out14_8x16b = _mm_add_epi16(src14_8x16b, src15_8x16b);
		out15_8x16b = _mm_sub_epi16(src14_8x16b, src15_8x16b);
		/****************Horizontal Transform Addition****************************/

		src0_8x16b = out0_8x16b;
		src1_8x16b = out1_8x16b;
		src2_8x16b = out2_8x16b;
		src3_8x16b = out3_8x16b;
		src4_8x16b = out4_8x16b;
		src5_8x16b = out5_8x16b;
		src6_8x16b = out6_8x16b;
		src7_8x16b = out7_8x16b;
		src8_8x16b = out8_8x16b;
		src9_8x16b = out9_8x16b;
		src10_8x16b = out10_8x16b;
		src11_8x16b = out11_8x16b;
		src12_8x16b = out12_8x16b;
		src13_8x16b = out13_8x16b;
		src14_8x16b = out14_8x16b;
		src15_8x16b = out15_8x16b;

		if (bitDepth == 8)
		{
			/************************* 8x8 Vertical Transform*************************/
			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi16(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi16(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi16(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi16(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi16(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi16(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi16(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi16(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/************************* 8x8 Vertical Transform*************************/
			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi16(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi16(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi16(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi16(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi16(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi16(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi16(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi16(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi16(out0_8x16b);
			src1_8x16b = _mm_abs_epi16(out1_8x16b);
			src2_8x16b = _mm_abs_epi16(out2_8x16b);
			src3_8x16b = _mm_abs_epi16(out3_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out4_8x16b);
			src1_8x16b = _mm_abs_epi16(out5_8x16b);
			src2_8x16b = _mm_abs_epi16(out6_8x16b);
			src3_8x16b = _mm_abs_epi16(out7_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out8_8x16b);
			src1_8x16b = _mm_abs_epi16(out9_8x16b);
			src2_8x16b = _mm_abs_epi16(out10_8x16b);
			src3_8x16b = _mm_abs_epi16(out11_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out12_8x16b);
			src1_8x16b = _mm_abs_epi16(out13_8x16b);
			src2_8x16b = _mm_abs_epi16(out14_8x16b);
			src3_8x16b = _mm_abs_epi16(out15_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (Int)(sad / SAD_DIV_16x8);

			return sad;
		}
		else
		{
			__m128i out0a_8x16b, out1a_8x16b, out2a_8x16b, out3a_8x16b;
			__m128i out4a_8x16b, out5a_8x16b, out6a_8x16b, out7a_8x16b;
			__m128i out8a_8x16b, out9a_8x16b, out10a_8x16b, out11a_8x16b;
			__m128i out12a_8x16b, out13a_8x16b, out14a_8x16b, out15a_8x16b;
			__m128i tmp0_8x16b, tmp1_8x16b, tmp2_8x16b, tmp3_8x16b;
			__m128i tmp4_8x16b, tmp5_8x16b, tmp6_8x16b, tmp7_8x16b;
			__m128i tmp8_8x16b, tmp9_8x16b, tmp10_8x16b, tmp11_8x16b;
			__m128i tmp12_8x16b, tmp13_8x16b, tmp14_8x16b, tmp15_8x16b;

			/************************* 8x8 Vertical Transform*************************/
			tmp0_8x16b = _mm_srli_si128(src0_8x16b, 8);
			tmp2_8x16b = _mm_srli_si128(src2_8x16b, 8);
			tmp4_8x16b = _mm_srli_si128(src4_8x16b, 8);
			tmp6_8x16b = _mm_srli_si128(src6_8x16b, 8);
			tmp8_8x16b = _mm_srli_si128(src8_8x16b, 8);
			tmp10_8x16b = _mm_srli_si128(src10_8x16b, 8);
			tmp12_8x16b = _mm_srli_si128(src12_8x16b, 8);
			tmp14_8x16b = _mm_srli_si128(src14_8x16b, 8);

			/*************************First 4 pixels ********************************/
			src0_8x16b = _mm_cvtepi16_epi32(src0_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(src2_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(src4_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(src6_8x16b);
			src8_8x16b = _mm_cvtepi16_epi32(src8_8x16b);
			src10_8x16b = _mm_cvtepi16_epi32(src10_8x16b);
			src12_8x16b = _mm_cvtepi16_epi32(src12_8x16b);
			src14_8x16b = _mm_cvtepi16_epi32(src14_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/*************************First 4 pixels ********************************/

			/**************************Next 4 pixels *******************************/
			src0_8x16b = _mm_cvtepi16_epi32(tmp0_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(tmp2_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(tmp4_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(tmp6_8x16b);
			src8_8x16b = _mm_cvtepi16_epi32(tmp8_8x16b);
			src10_8x16b = _mm_cvtepi16_epi32(tmp10_8x16b);
			src12_8x16b = _mm_cvtepi16_epi32(tmp12_8x16b);
			src14_8x16b = _mm_cvtepi16_epi32(tmp14_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/**************************Next 4 pixels *******************************/
			/************************* 8x8 Vertical Transform*************************/

			/************************* 8x8 Vertical Transform*************************/
			tmp1_8x16b = _mm_srli_si128(src1_8x16b, 8);
			tmp3_8x16b = _mm_srli_si128(src3_8x16b, 8);
			tmp5_8x16b = _mm_srli_si128(src5_8x16b, 8);
			tmp7_8x16b = _mm_srli_si128(src7_8x16b, 8);
			tmp9_8x16b = _mm_srli_si128(src9_8x16b, 8);
			tmp11_8x16b = _mm_srli_si128(src11_8x16b, 8);
			tmp13_8x16b = _mm_srli_si128(src13_8x16b, 8);
			tmp15_8x16b = _mm_srli_si128(src15_8x16b, 8);

			/*************************First 4 pixels ********************************/
			src1_8x16b = _mm_cvtepi16_epi32(src1_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(src3_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(src5_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(src7_8x16b);
			src9_8x16b = _mm_cvtepi16_epi32(src9_8x16b);
			src11_8x16b = _mm_cvtepi16_epi32(src11_8x16b);
			src13_8x16b = _mm_cvtepi16_epi32(src13_8x16b);
			src15_8x16b = _mm_cvtepi16_epi32(src15_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/*************************First 4 pixels ********************************/

			/*************************Next 4 pixels ********************************/
			src1_8x16b = _mm_cvtepi16_epi32(tmp1_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(tmp3_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(tmp5_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(tmp7_8x16b);
			src9_8x16b = _mm_cvtepi16_epi32(tmp9_8x16b);
			src11_8x16b = _mm_cvtepi16_epi32(tmp11_8x16b);
			src13_8x16b = _mm_cvtepi16_epi32(tmp13_8x16b);
			src15_8x16b = _mm_cvtepi16_epi32(tmp15_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi32(out0_8x16b);
			src1_8x16b = _mm_abs_epi32(out1_8x16b);
			src2_8x16b = _mm_abs_epi32(out2_8x16b);
			src3_8x16b = _mm_abs_epi32(out3_8x16b);
			src4_8x16b = _mm_abs_epi32(out4_8x16b);
			src5_8x16b = _mm_abs_epi32(out5_8x16b);
			src6_8x16b = _mm_abs_epi32(out6_8x16b);
			src7_8x16b = _mm_abs_epi32(out7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out8_8x16b);
			src1_8x16b = _mm_abs_epi32(out9_8x16b);
			src2_8x16b = _mm_abs_epi32(out10_8x16b);
			src3_8x16b = _mm_abs_epi32(out11_8x16b);
			src4_8x16b = _mm_abs_epi32(out12_8x16b);
			src5_8x16b = _mm_abs_epi32(out13_8x16b);
			src6_8x16b = _mm_abs_epi32(out14_8x16b);
			src7_8x16b = _mm_abs_epi32(out15_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out0a_8x16b);
			src1_8x16b = _mm_abs_epi32(out1a_8x16b);
			src2_8x16b = _mm_abs_epi32(out2a_8x16b);
			src3_8x16b = _mm_abs_epi32(out3a_8x16b);
			src4_8x16b = _mm_abs_epi32(out4a_8x16b);
			src5_8x16b = _mm_abs_epi32(out5a_8x16b);
			src6_8x16b = _mm_abs_epi32(out6a_8x16b);
			src7_8x16b = _mm_abs_epi32(out7a_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out8a_8x16b);
			src1_8x16b = _mm_abs_epi32(out9a_8x16b);
			src2_8x16b = _mm_abs_epi32(out10a_8x16b);
			src3_8x16b = _mm_abs_epi32(out11a_8x16b);
			src4_8x16b = _mm_abs_epi32(out12a_8x16b);
			src5_8x16b = _mm_abs_epi32(out13a_8x16b);
			src6_8x16b = _mm_abs_epi32(out14a_8x16b);
			src7_8x16b = _mm_abs_epi32(out15a_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (Int)(sad / SAD_DIV_16x8);

			return sad;
		}
	}
	else
	{
		Int k, i, j, jj, sad = 0;
		Int diff[128], m1[8][16], m2[8][16];
		for (k = 0; k < 128; k += 16)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];
			diff[k + 4] = piOrg[4] - piCur[4];
			diff[k + 5] = piOrg[5] - piCur[5];
			diff[k + 6] = piOrg[6] - piCur[6];
			diff[k + 7] = piOrg[7] - piCur[7];

			diff[k + 8] = piOrg[8] - piCur[8];
			diff[k + 9] = piOrg[9] - piCur[9];
			diff[k + 10] = piOrg[10] - piCur[10];
			diff[k + 11] = piOrg[11] - piCur[11];
			diff[k + 12] = piOrg[12] - piCur[12];
			diff[k + 13] = piOrg[13] - piCur[13];
			diff[k + 14] = piOrg[14] - piCur[14];
			diff[k + 15] = piOrg[15] - piCur[15];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		//horizontal
		for (j = 0; j < 8; j++)
		{
			jj = j << 4;

			m2[j][0] = diff[jj] + diff[jj + 8];
			m2[j][1] = diff[jj + 1] + diff[jj + 9];
			m2[j][2] = diff[jj + 2] + diff[jj + 10];
			m2[j][3] = diff[jj + 3] + diff[jj + 11];
			m2[j][4] = diff[jj + 4] + diff[jj + 12];
			m2[j][5] = diff[jj + 5] + diff[jj + 13];
			m2[j][6] = diff[jj + 6] + diff[jj + 14];
			m2[j][7] = diff[jj + 7] + diff[jj + 15];
			m2[j][8] = diff[jj] - diff[jj + 8];
			m2[j][9] = diff[jj + 1] - diff[jj + 9];
			m2[j][10] = diff[jj + 2] - diff[jj + 10];
			m2[j][11] = diff[jj + 3] - diff[jj + 11];
			m2[j][12] = diff[jj + 4] - diff[jj + 12];
			m2[j][13] = diff[jj + 5] - diff[jj + 13];
			m2[j][14] = diff[jj + 6] - diff[jj + 14];
			m2[j][15] = diff[jj + 7] - diff[jj + 15];

			m1[j][0] = m2[j][0] + m2[j][4];
			m1[j][1] = m2[j][1] + m2[j][5];
			m1[j][2] = m2[j][2] + m2[j][6];
			m1[j][3] = m2[j][3] + m2[j][7];
			m1[j][4] = m2[j][0] - m2[j][4];
			m1[j][5] = m2[j][1] - m2[j][5];
			m1[j][6] = m2[j][2] - m2[j][6];
			m1[j][7] = m2[j][3] - m2[j][7];
			m1[j][8] = m2[j][8] + m2[j][12];
			m1[j][9] = m2[j][9] + m2[j][13];
			m1[j][10] = m2[j][10] + m2[j][14];
			m1[j][11] = m2[j][11] + m2[j][15];
			m1[j][12] = m2[j][8] - m2[j][12];
			m1[j][13] = m2[j][9] - m2[j][13];
			m1[j][14] = m2[j][10] - m2[j][14];
			m1[j][15] = m2[j][11] - m2[j][15];

			m2[j][0] = m1[j][0] + m1[j][2];
			m2[j][1] = m1[j][1] + m1[j][3];
			m2[j][2] = m1[j][0] - m1[j][2];
			m2[j][3] = m1[j][1] - m1[j][3];
			m2[j][4] = m1[j][4] + m1[j][6];
			m2[j][5] = m1[j][5] + m1[j][7];
			m2[j][6] = m1[j][4] - m1[j][6];
			m2[j][7] = m1[j][5] - m1[j][7];
			m2[j][8] = m1[j][8] + m1[j][10];
			m2[j][9] = m1[j][9] + m1[j][11];
			m2[j][10] = m1[j][8] - m1[j][10];
			m2[j][11] = m1[j][9] - m1[j][11];
			m2[j][12] = m1[j][12] + m1[j][14];
			m2[j][13] = m1[j][13] + m1[j][15];
			m2[j][14] = m1[j][12] - m1[j][14];
			m2[j][15] = m1[j][13] - m1[j][15];

			m1[j][0] = m2[j][0] + m2[j][1];
			m1[j][1] = m2[j][0] - m2[j][1];
			m1[j][2] = m2[j][2] + m2[j][3];
			m1[j][3] = m2[j][2] - m2[j][3];
			m1[j][4] = m2[j][4] + m2[j][5];
			m1[j][5] = m2[j][4] - m2[j][5];
			m1[j][6] = m2[j][6] + m2[j][7];
			m1[j][7] = m2[j][6] - m2[j][7];
			m1[j][8] = m2[j][8] + m2[j][9];
			m1[j][9] = m2[j][8] - m2[j][9];
			m1[j][10] = m2[j][10] + m2[j][11];
			m1[j][11] = m2[j][10] - m2[j][11];
			m1[j][12] = m2[j][12] + m2[j][13];
			m1[j][13] = m2[j][12] - m2[j][13];
			m1[j][14] = m2[j][14] + m2[j][15];
			m1[j][15] = m2[j][14] - m2[j][15];
		}

		//vertical
		for (i = 0; i < 16; i++)
		{
			m2[0][i] = m1[0][i] + m1[4][i];
			m2[1][i] = m1[1][i] + m1[5][i];
			m2[2][i] = m1[2][i] + m1[6][i];
			m2[3][i] = m1[3][i] + m1[7][i];
			m2[4][i] = m1[0][i] - m1[4][i];
			m2[5][i] = m1[1][i] - m1[5][i];
			m2[6][i] = m1[2][i] - m1[6][i];
			m2[7][i] = m1[3][i] - m1[7][i];

			m1[0][i] = m2[0][i] + m2[2][i];
			m1[1][i] = m2[1][i] + m2[3][i];
			m1[2][i] = m2[0][i] - m2[2][i];
			m1[3][i] = m2[1][i] - m2[3][i];
			m1[4][i] = m2[4][i] + m2[6][i];
			m1[5][i] = m2[5][i] + m2[7][i];
			m1[6][i] = m2[4][i] - m2[6][i];
			m1[7][i] = m2[5][i] - m2[7][i];

			m2[0][i] = m1[0][i] + m1[1][i];
			m2[1][i] = m1[0][i] - m1[1][i];
			m2[2][i] = m1[2][i] + m1[3][i];
			m2[3][i] = m1[2][i] - m1[3][i];
			m2[4][i] = m1[4][i] + m1[5][i];
			m2[5][i] = m1[4][i] - m1[5][i];
			m2[6][i] = m1[6][i] + m1[7][i];
			m2[7][i] = m1[6][i] - m1[7][i];
		}

		for (i = 0; i < 8; i++)
		{
			for (j = 0; j < 16; j++)
			{
				sad += abs(m2[i][j]);
			}
		}

		sad = (Int)(sad / SAD_DIV_16x8);

		return sad;
	}
}

#define SAD_DIV_8x16 sqrt(16.0 * 8) * 2
Distortion TComRdCost::xCalcHADs8x16(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i src4_8x16b, src5_8x16b, src6_8x16b, src7_8x16b;
		__m128i src8_8x16b, src9_8x16b, src10_8x16b, src11_8x16b;
		__m128i src12_8x16b, src13_8x16b, src14_8x16b, src15_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i pred4_8x16b, pred5_8x16b, pred6_8x16b, pred7_8x16b;
		__m128i pred8_8x16b, pred9_8x16b, pred10_8x16b, pred11_8x16b;
		__m128i pred12_8x16b, pred13_8x16b, pred14_8x16b, pred15_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;
		__m128i out4_8x16b, out5_8x16b, out6_8x16b, out7_8x16b;
		__m128i out8_8x16b, out9_8x16b, out10_8x16b, out11_8x16b;
		__m128i out12_8x16b, out13_8x16b, out14_8x16b, out15_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src1_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src3_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src4_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src5_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src6_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src7_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;

		pred0_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred1_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred3_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred4_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred5_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred6_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred7_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_sub_epi16(src4_8x16b, pred4_8x16b);
		src5_8x16b = _mm_sub_epi16(src5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_sub_epi16(src6_8x16b, pred6_8x16b);
		src7_8x16b = _mm_sub_epi16(src7_8x16b, pred7_8x16b);

		src8_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src9_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src10_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src11_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src12_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src13_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src14_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src15_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;

		pred8_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred9_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred10_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred11_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred12_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred13_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred14_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred15_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;

		src8_8x16b = _mm_sub_epi16(src8_8x16b, pred8_8x16b);
		src9_8x16b = _mm_sub_epi16(src9_8x16b, pred9_8x16b);
		src10_8x16b = _mm_sub_epi16(src10_8x16b, pred10_8x16b);
		src11_8x16b = _mm_sub_epi16(src11_8x16b, pred11_8x16b);
		src12_8x16b = _mm_sub_epi16(src12_8x16b, pred12_8x16b);
		src13_8x16b = _mm_sub_epi16(src13_8x16b, pred13_8x16b);
		src14_8x16b = _mm_sub_epi16(src14_8x16b, pred14_8x16b);
		src15_8x16b = _mm_sub_epi16(src15_8x16b, pred15_8x16b);
		/**********************Residue Calculation********************************/

		/**************** 8x8 horizontal transform *******************************/
		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		out0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		out1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		out2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		out3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		out4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		out5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		out6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		out7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_unpackhi_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src6_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src6_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		src5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		src7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/
		/**************** 8x8 horizontal transform *******************************/

		/**************** 8x8 horizontal transform *******************************/
		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src8_8x16b, src9_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src10_8x16b, src11_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src12_8x16b, src13_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src14_8x16b, src15_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src8_8x16b, src9_8x16b);
		src10_8x16b = _mm_unpackhi_epi16(src10_8x16b, src11_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src12_8x16b, src13_8x16b);
		src14_8x16b = _mm_unpackhi_epi16(src14_8x16b, src15_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src10_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src10_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src14_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src14_8x16b);

		out0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		out1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		out2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		out3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		out4_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		out5_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		out6_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		out7_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src8_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src12_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src10_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src14_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src9_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src13_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src11_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src15_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x8 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src8_8x16b, src9_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src10_8x16b, src11_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src12_8x16b, src13_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src14_8x16b, src15_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src8_8x16b, src9_8x16b);
		src10_8x16b = _mm_unpackhi_epi16(src10_8x16b, src11_8x16b);
		pred7_8x16b = _mm_unpackhi_epi16(src12_8x16b, src13_8x16b);
		src14_8x16b = _mm_unpackhi_epi16(src14_8x16b, src15_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);
		out5_8x16b = _mm_unpacklo_epi32(out7_8x16b, src10_8x16b);
		out7_8x16b = _mm_unpackhi_epi32(out7_8x16b, src10_8x16b);
		pred5_8x16b = _mm_unpacklo_epi32(pred7_8x16b, src14_8x16b);
		pred7_8x16b = _mm_unpackhi_epi32(pred7_8x16b, src14_8x16b);

		src8_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src9_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src10_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src11_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		src12_8x16b = _mm_unpacklo_epi64(out5_8x16b, pred5_8x16b);
		src13_8x16b = _mm_unpackhi_epi64(out5_8x16b, pred5_8x16b);
		src14_8x16b = _mm_unpacklo_epi64(out7_8x16b, pred7_8x16b);
		src15_8x16b = _mm_unpackhi_epi64(out7_8x16b, pred7_8x16b);
		/**********************   8x8 16 bit Transpose End   *********************/
		/**************** 8x8 horizontal transform *******************************/

		/****************Horizontal Transform Addition****************************/
		out0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
		out1_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);

		out2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);
		out3_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);

		out4_8x16b = _mm_add_epi16(src4_8x16b, src5_8x16b);
		out5_8x16b = _mm_sub_epi16(src4_8x16b, src5_8x16b);

		out6_8x16b = _mm_add_epi16(src6_8x16b, src7_8x16b);
		out7_8x16b = _mm_sub_epi16(src6_8x16b, src7_8x16b);

		out8_8x16b = _mm_add_epi16(src8_8x16b, src9_8x16b);
		out9_8x16b = _mm_sub_epi16(src8_8x16b, src9_8x16b);

		out10_8x16b = _mm_add_epi16(src10_8x16b, src11_8x16b);
		out11_8x16b = _mm_sub_epi16(src10_8x16b, src11_8x16b);

		out12_8x16b = _mm_add_epi16(src12_8x16b, src13_8x16b);
		out13_8x16b = _mm_sub_epi16(src12_8x16b, src13_8x16b);

		out14_8x16b = _mm_add_epi16(src14_8x16b, src15_8x16b);
		out15_8x16b = _mm_sub_epi16(src14_8x16b, src15_8x16b);
		/****************Horizontal Transform Addition****************************/

		src0_8x16b = out0_8x16b;
		src1_8x16b = out1_8x16b;
		src2_8x16b = out2_8x16b;
		src3_8x16b = out3_8x16b;
		src4_8x16b = out4_8x16b;
		src5_8x16b = out5_8x16b;
		src6_8x16b = out6_8x16b;
		src7_8x16b = out7_8x16b;
		src8_8x16b = out8_8x16b;
		src9_8x16b = out9_8x16b;
		src10_8x16b = out10_8x16b;
		src11_8x16b = out11_8x16b;
		src12_8x16b = out12_8x16b;
		src13_8x16b = out13_8x16b;
		src14_8x16b = out14_8x16b;
		src15_8x16b = out15_8x16b;

		if (bitDepth == 8)
		{
			/************************* 8x8 Vertical Transform*************************/
			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi16(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi16(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi16(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi16(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi16(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi16(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi16(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi16(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/************************* 8x8 Vertical Transform*************************/
			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi16(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi16(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi16(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi16(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi16(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi16(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi16(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi16(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi16(out0_8x16b);
			src1_8x16b = _mm_abs_epi16(out1_8x16b);
			src2_8x16b = _mm_abs_epi16(out2_8x16b);
			src3_8x16b = _mm_abs_epi16(out3_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out4_8x16b);
			src1_8x16b = _mm_abs_epi16(out5_8x16b);
			src2_8x16b = _mm_abs_epi16(out6_8x16b);
			src3_8x16b = _mm_abs_epi16(out7_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out8_8x16b);
			src1_8x16b = _mm_abs_epi16(out9_8x16b);
			src2_8x16b = _mm_abs_epi16(out10_8x16b);
			src3_8x16b = _mm_abs_epi16(out11_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi16(out12_8x16b);
			src1_8x16b = _mm_abs_epi16(out13_8x16b);
			src2_8x16b = _mm_abs_epi16(out14_8x16b);
			src3_8x16b = _mm_abs_epi16(out15_8x16b);

			src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
			src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
			src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
			src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

			src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
			src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
			src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
			src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
			src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
			src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
			src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
			src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (Int)(sad / SAD_DIV_8x16);

			return sad;
		}
		else
		{
			__m128i out0a_8x16b, out1a_8x16b, out2a_8x16b, out3a_8x16b;
			__m128i out4a_8x16b, out5a_8x16b, out6a_8x16b, out7a_8x16b;
			__m128i out8a_8x16b, out9a_8x16b, out10a_8x16b, out11a_8x16b;
			__m128i out12a_8x16b, out13a_8x16b, out14a_8x16b, out15a_8x16b;
			__m128i tmp0_8x16b, tmp1_8x16b, tmp2_8x16b, tmp3_8x16b;
			__m128i tmp4_8x16b, tmp5_8x16b, tmp6_8x16b, tmp7_8x16b;
			__m128i tmp8_8x16b, tmp9_8x16b, tmp10_8x16b, tmp11_8x16b;
			__m128i tmp12_8x16b, tmp13_8x16b, tmp14_8x16b, tmp15_8x16b;

			/************************* 8x8 Vertical Transform*************************/
			tmp0_8x16b = _mm_srli_si128(src0_8x16b, 8);
			tmp2_8x16b = _mm_srli_si128(src2_8x16b, 8);
			tmp4_8x16b = _mm_srli_si128(src4_8x16b, 8);
			tmp6_8x16b = _mm_srli_si128(src6_8x16b, 8);
			tmp8_8x16b = _mm_srli_si128(src8_8x16b, 8);
			tmp10_8x16b = _mm_srli_si128(src10_8x16b, 8);
			tmp12_8x16b = _mm_srli_si128(src12_8x16b, 8);
			tmp14_8x16b = _mm_srli_si128(src14_8x16b, 8);

			/*************************First 4 pixels ********************************/
			src0_8x16b = _mm_cvtepi16_epi32(src0_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(src2_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(src4_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(src6_8x16b);
			src8_8x16b = _mm_cvtepi16_epi32(src8_8x16b);
			src10_8x16b = _mm_cvtepi16_epi32(src10_8x16b);
			src12_8x16b = _mm_cvtepi16_epi32(src12_8x16b);
			src14_8x16b = _mm_cvtepi16_epi32(src14_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/*************************First 4 pixels ********************************/

			/**************************Next 4 pixels *******************************/
			src0_8x16b = _mm_cvtepi16_epi32(tmp0_8x16b);
			src2_8x16b = _mm_cvtepi16_epi32(tmp2_8x16b);
			src4_8x16b = _mm_cvtepi16_epi32(tmp4_8x16b);
			src6_8x16b = _mm_cvtepi16_epi32(tmp6_8x16b);
			src8_8x16b = _mm_cvtepi16_epi32(tmp8_8x16b);
			src10_8x16b = _mm_cvtepi16_epi32(tmp10_8x16b);
			src12_8x16b = _mm_cvtepi16_epi32(tmp12_8x16b);
			src14_8x16b = _mm_cvtepi16_epi32(tmp14_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src8_8x16b, src10_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src12_8x16b, src14_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out0a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out8a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out4a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out12a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src0_8x16b, src2_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src4_8x16b, src6_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src8_8x16b, src10_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src12_8x16b, src14_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out2a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out10a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out6a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out14a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/**************************Next 4 pixels *******************************/
			/************************* 8x8 Vertical Transform*************************/

			/************************* 8x8 Vertical Transform*************************/
			tmp1_8x16b = _mm_srli_si128(src1_8x16b, 8);
			tmp3_8x16b = _mm_srli_si128(src3_8x16b, 8);
			tmp5_8x16b = _mm_srli_si128(src5_8x16b, 8);
			tmp7_8x16b = _mm_srli_si128(src7_8x16b, 8);
			tmp9_8x16b = _mm_srli_si128(src9_8x16b, 8);
			tmp11_8x16b = _mm_srli_si128(src11_8x16b, 8);
			tmp13_8x16b = _mm_srli_si128(src13_8x16b, 8);
			tmp15_8x16b = _mm_srli_si128(src15_8x16b, 8);

			/*************************First 4 pixels ********************************/
			src1_8x16b = _mm_cvtepi16_epi32(src1_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(src3_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(src5_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(src7_8x16b);
			src9_8x16b = _mm_cvtepi16_epi32(src9_8x16b);
			src11_8x16b = _mm_cvtepi16_epi32(src11_8x16b);
			src13_8x16b = _mm_cvtepi16_epi32(src13_8x16b);
			src15_8x16b = _mm_cvtepi16_epi32(src15_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/*************************First 4 pixels ********************************/

			/*************************Next 4 pixels ********************************/
			src1_8x16b = _mm_cvtepi16_epi32(tmp1_8x16b);
			src3_8x16b = _mm_cvtepi16_epi32(tmp3_8x16b);
			src5_8x16b = _mm_cvtepi16_epi32(tmp5_8x16b);
			src7_8x16b = _mm_cvtepi16_epi32(tmp7_8x16b);
			src9_8x16b = _mm_cvtepi16_epi32(tmp9_8x16b);
			src11_8x16b = _mm_cvtepi16_epi32(tmp11_8x16b);
			src13_8x16b = _mm_cvtepi16_epi32(tmp13_8x16b);
			src15_8x16b = _mm_cvtepi16_epi32(tmp15_8x16b);

			/* r0 + r1 */
			pred0_8x16b = _mm_add_epi32(src1_8x16b, src3_8x16b);
			/* r2 + r3 */
			pred2_8x16b = _mm_add_epi32(src5_8x16b, src7_8x16b);
			/* r4 + r5 */
			pred4_8x16b = _mm_add_epi32(src9_8x16b, src11_8x16b);
			/* r6 + r7 */
			pred6_8x16b = _mm_add_epi32(src13_8x16b, src15_8x16b);

			/* r0 + r1 + r2 + r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 + r6 + r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
			out1a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
			out9a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 + r1 - r2 - r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 + r5 - r6 - r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
			out5a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
			out13a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 */
			pred0_8x16b = _mm_sub_epi32(src1_8x16b, src3_8x16b);
			/* r2 - r3 */
			pred2_8x16b = _mm_sub_epi32(src5_8x16b, src7_8x16b);
			/* r4 - r5 */
			pred4_8x16b = _mm_sub_epi32(src9_8x16b, src11_8x16b);
			/* r6 - r7 */
			pred6_8x16b = _mm_sub_epi32(src13_8x16b, src15_8x16b);

			/* r0 - r1 + r2 - r3 */
			pred1_8x16b = _mm_add_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 + r6 - r7 */
			pred5_8x16b = _mm_add_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
			out3a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
			out11a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);

			/* r0 - r1 - r2 + r3 */
			pred1_8x16b = _mm_sub_epi32(pred0_8x16b, pred2_8x16b);
			/* r4 - r5 - r6 + r7 */
			pred5_8x16b = _mm_sub_epi32(pred4_8x16b, pred6_8x16b);
			/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
			out7a_8x16b = _mm_add_epi32(pred1_8x16b, pred5_8x16b);
			/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
			out15a_8x16b = _mm_sub_epi32(pred1_8x16b, pred5_8x16b);
			/************************* 8x8 Vertical Transform*************************/

			/****************************SATD calculation ****************************/
			src0_8x16b = _mm_abs_epi32(out0_8x16b);
			src1_8x16b = _mm_abs_epi32(out1_8x16b);
			src2_8x16b = _mm_abs_epi32(out2_8x16b);
			src3_8x16b = _mm_abs_epi32(out3_8x16b);
			src4_8x16b = _mm_abs_epi32(out4_8x16b);
			src5_8x16b = _mm_abs_epi32(out5_8x16b);
			src6_8x16b = _mm_abs_epi32(out6_8x16b);
			src7_8x16b = _mm_abs_epi32(out7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out8_8x16b);
			src1_8x16b = _mm_abs_epi32(out9_8x16b);
			src2_8x16b = _mm_abs_epi32(out10_8x16b);
			src3_8x16b = _mm_abs_epi32(out11_8x16b);
			src4_8x16b = _mm_abs_epi32(out12_8x16b);
			src5_8x16b = _mm_abs_epi32(out13_8x16b);
			src6_8x16b = _mm_abs_epi32(out14_8x16b);
			src7_8x16b = _mm_abs_epi32(out15_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out0a_8x16b);
			src1_8x16b = _mm_abs_epi32(out1a_8x16b);
			src2_8x16b = _mm_abs_epi32(out2a_8x16b);
			src3_8x16b = _mm_abs_epi32(out3a_8x16b);
			src4_8x16b = _mm_abs_epi32(out4a_8x16b);
			src5_8x16b = _mm_abs_epi32(out5a_8x16b);
			src6_8x16b = _mm_abs_epi32(out6a_8x16b);
			src7_8x16b = _mm_abs_epi32(out7a_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			src0_8x16b = _mm_abs_epi32(out8a_8x16b);
			src1_8x16b = _mm_abs_epi32(out9a_8x16b);
			src2_8x16b = _mm_abs_epi32(out10a_8x16b);
			src3_8x16b = _mm_abs_epi32(out11a_8x16b);
			src4_8x16b = _mm_abs_epi32(out12a_8x16b);
			src5_8x16b = _mm_abs_epi32(out13a_8x16b);
			src6_8x16b = _mm_abs_epi32(out14a_8x16b);
			src7_8x16b = _mm_abs_epi32(out15a_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
			src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
			src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
			src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

			src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
			src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

			sad += _mm_cvtsi128_si32(src0_8x16b);

			sad = (Int)(sad / SAD_DIV_8x16);

			return sad;
		}
	}
	else
	{
		Int k, i, j, jj, sad = 0;
		Int diff[128], m1[16][8], m2[16][8];
		for (k = 0; k < 128; k += 8)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];
			diff[k + 4] = piOrg[4] - piCur[4];
			diff[k + 5] = piOrg[5] - piCur[5];
			diff[k + 6] = piOrg[6] - piCur[6];
			diff[k + 7] = piOrg[7] - piCur[7];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		//horizontal
		for (j = 0; j < 16; j++)
		{
			jj = j << 3;

			m2[j][0] = diff[jj] + diff[jj + 4];
			m2[j][1] = diff[jj + 1] + diff[jj + 5];
			m2[j][2] = diff[jj + 2] + diff[jj + 6];
			m2[j][3] = diff[jj + 3] + diff[jj + 7];
			m2[j][4] = diff[jj] - diff[jj + 4];
			m2[j][5] = diff[jj + 1] - diff[jj + 5];
			m2[j][6] = diff[jj + 2] - diff[jj + 6];
			m2[j][7] = diff[jj + 3] - diff[jj + 7];

			m1[j][0] = m2[j][0] + m2[j][2];
			m1[j][1] = m2[j][1] + m2[j][3];
			m1[j][2] = m2[j][0] - m2[j][2];
			m1[j][3] = m2[j][1] - m2[j][3];
			m1[j][4] = m2[j][4] + m2[j][6];
			m1[j][5] = m2[j][5] + m2[j][7];
			m1[j][6] = m2[j][4] - m2[j][6];
			m1[j][7] = m2[j][5] - m2[j][7];

			m2[j][0] = m1[j][0] + m1[j][1];
			m2[j][1] = m1[j][0] - m1[j][1];
			m2[j][2] = m1[j][2] + m1[j][3];
			m2[j][3] = m1[j][2] - m1[j][3];
			m2[j][4] = m1[j][4] + m1[j][5];
			m2[j][5] = m1[j][4] - m1[j][5];
			m2[j][6] = m1[j][6] + m1[j][7];
			m2[j][7] = m1[j][6] - m1[j][7];
		}

		//vertical
		for (i = 0; i < 8; i++)
		{
			m1[0][i] = m2[0][i] + m2[8][i];
			m1[1][i] = m2[1][i] + m2[9][i];
			m1[2][i] = m2[2][i] + m2[10][i];
			m1[3][i] = m2[3][i] + m2[11][i];
			m1[4][i] = m2[4][i] + m2[12][i];
			m1[5][i] = m2[5][i] + m2[13][i];
			m1[6][i] = m2[6][i] + m2[14][i];
			m1[7][i] = m2[7][i] + m2[15][i];
			m1[8][i] = m2[0][i] - m2[8][i];
			m1[9][i] = m2[1][i] - m2[9][i];
			m1[10][i] = m2[2][i] - m2[10][i];
			m1[11][i] = m2[3][i] - m2[11][i];
			m1[12][i] = m2[4][i] - m2[12][i];
			m1[13][i] = m2[5][i] - m2[13][i];
			m1[14][i] = m2[6][i] - m2[14][i];
			m1[15][i] = m2[7][i] - m2[15][i];

			m2[0][i] = m1[0][i] + m1[4][i];
			m2[1][i] = m1[1][i] + m1[5][i];
			m2[2][i] = m1[2][i] + m1[6][i];
			m2[3][i] = m1[3][i] + m1[7][i];
			m2[4][i] = m1[0][i] - m1[4][i];
			m2[5][i] = m1[1][i] - m1[5][i];
			m2[6][i] = m1[2][i] - m1[6][i];
			m2[7][i] = m1[3][i] - m1[7][i];
			m2[8][i] = m1[8][i] + m1[12][i];
			m2[9][i] = m1[9][i] + m1[13][i];
			m2[10][i] = m1[10][i] + m1[14][i];
			m2[11][i] = m1[11][i] + m1[15][i];
			m2[12][i] = m1[8][i] - m1[12][i];
			m2[13][i] = m1[9][i] - m1[13][i];
			m2[14][i] = m1[10][i] - m1[14][i];
			m2[15][i] = m1[11][i] - m1[15][i];

			m1[0][i] = m2[0][i] + m2[2][i];
			m1[1][i] = m2[1][i] + m2[3][i];
			m1[2][i] = m2[0][i] - m2[2][i];
			m1[3][i] = m2[1][i] - m2[3][i];
			m1[4][i] = m2[4][i] + m2[6][i];
			m1[5][i] = m2[5][i] + m2[7][i];
			m1[6][i] = m2[4][i] - m2[6][i];
			m1[7][i] = m2[5][i] - m2[7][i];
			m1[8][i] = m2[8][i] + m2[10][i];
			m1[9][i] = m2[9][i] + m2[11][i];
			m1[10][i] = m2[8][i] - m2[10][i];
			m1[11][i] = m2[9][i] - m2[11][i];
			m1[12][i] = m2[12][i] + m2[14][i];
			m1[13][i] = m2[13][i] + m2[15][i];
			m1[14][i] = m2[12][i] - m2[14][i];
			m1[15][i] = m2[13][i] - m2[15][i];

			m2[0][i] = m1[0][i] + m1[1][i];
			m2[1][i] = m1[0][i] - m1[1][i];
			m2[2][i] = m1[2][i] + m1[3][i];
			m2[3][i] = m1[2][i] - m1[3][i];
			m2[4][i] = m1[4][i] + m1[5][i];
			m2[5][i] = m1[4][i] - m1[5][i];
			m2[6][i] = m1[6][i] + m1[7][i];
			m2[7][i] = m1[6][i] - m1[7][i];
			m2[8][i] = m1[8][i] + m1[9][i];
			m2[9][i] = m1[8][i] - m1[9][i];
			m2[10][i] = m1[10][i] + m1[11][i];
			m2[11][i] = m1[10][i] - m1[11][i];
			m2[12][i] = m1[12][i] + m1[13][i];
			m2[13][i] = m1[12][i] - m1[13][i];
			m2[14][i] = m1[14][i] + m1[15][i];
			m2[15][i] = m1[14][i] - m1[15][i];
		}

		for (i = 0; i < 16; i++)
		{
			for (j = 0; j < 8; j++)
			{
				sad += abs(m2[i][j]);
			}
		}

		sad = (Int)(sad / SAD_DIV_8x16);

		return sad;
	}
}

#define SAD_DIV_4x8 sqrt(4.0*8)*2
Distortion TComRdCost::xCalcHADs4x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i src4_8x16b, src5_8x16b, src6_8x16b, src7_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i pred4_8x16b, pred5_8x16b, pred6_8x16b, pred7_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;
		__m128i out4_8x16b, out5_8x16b, out6_8x16b, out7_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src1_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src3_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src4_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src5_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src6_8x16b = _mm_loadl_epi64((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src7_8x16b = _mm_loadl_epi64((__m128i *) piOrg);

		pred0_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred1_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred3_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred4_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred5_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred6_8x16b = _mm_loadl_epi64((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred7_8x16b = _mm_loadl_epi64((__m128i *) piCur);

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		src4_8x16b = _mm_sub_epi16(src4_8x16b, pred4_8x16b);
		src5_8x16b = _mm_sub_epi16(src5_8x16b, pred5_8x16b);
		src6_8x16b = _mm_sub_epi16(src6_8x16b, pred6_8x16b);
		src7_8x16b = _mm_sub_epi16(src7_8x16b, pred7_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(src0_8x16b, src1_8x16b);
		src1_8x16b = _mm_unpacklo_epi64(src2_8x16b, src3_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(src4_8x16b, src5_8x16b);
		src3_8x16b = _mm_unpacklo_epi64(src6_8x16b, src7_8x16b);
		/**********************Residue Calculation********************************/

		/**************** 8x4 horizontal transform *******************************/
		/***********************    8x4 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);

		out0_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out1_8x16b = _mm_srli_si128(out0_8x16b, 8);
		out2_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_srli_si128(out2_8x16b, 8);
		out4_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out5_8x16b = _mm_srli_si128(out4_8x16b, 8);
		out6_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_srli_si128(out6_8x16b, 8);
		/**********************   8x4 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x4 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		/**********************   8x4 16 bit Transpose End   *********************/
		/**************** 8x4 horizontal transform *******************************/

		/************************* 4x8 Vertical Transform*************************/
		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);

		/* r0 + r1 + r2 + r3 */
		out0_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 + r1 - r2 - r3 */
		out2_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);

		/* r0 - r1 + r2 - r3 */
		out1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 - r1 - r2 + r3 */
		out3_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/************************* 4x8 Vertical Transform*************************/

		/****************************SATD calculation ****************************/
		src0_8x16b = _mm_abs_epi16(out0_8x16b);
		src1_8x16b = _mm_abs_epi16(out1_8x16b);
		src2_8x16b = _mm_abs_epi16(out2_8x16b);
		src3_8x16b = _mm_abs_epi16(out3_8x16b);

		src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
		src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
		src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
		src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

		src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
		src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
		src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
		src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
		src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
		src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
		src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
		src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
		src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
		src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

		sad += _mm_cvtsi128_si32(src0_8x16b);

		sad = (Int)(sad / SAD_DIV_4x8);

		return sad;
	}
	else
	{
		Int k, i, j, jj, sad = 0;
		Int diff[32], m1[4][8], m2[4][8];
		for (k = 0; k < 32; k += 8)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];
			diff[k + 4] = piOrg[4] - piCur[4];
			diff[k + 5] = piOrg[5] - piCur[5];
			diff[k + 6] = piOrg[6] - piCur[6];
			diff[k + 7] = piOrg[7] - piCur[7];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		//horizontal
		for (j = 0; j < 4; j++)
		{
			jj = j << 3;

			m2[j][0] = diff[jj] + diff[jj + 4];
			m2[j][1] = diff[jj + 1] + diff[jj + 5];
			m2[j][2] = diff[jj + 2] + diff[jj + 6];
			m2[j][3] = diff[jj + 3] + diff[jj + 7];
			m2[j][4] = diff[jj] - diff[jj + 4];
			m2[j][5] = diff[jj + 1] - diff[jj + 5];
			m2[j][6] = diff[jj + 2] - diff[jj + 6];
			m2[j][7] = diff[jj + 3] - diff[jj + 7];

			m1[j][0] = m2[j][0] + m2[j][2];
			m1[j][1] = m2[j][1] + m2[j][3];
			m1[j][2] = m2[j][0] - m2[j][2];
			m1[j][3] = m2[j][1] - m2[j][3];
			m1[j][4] = m2[j][4] + m2[j][6];
			m1[j][5] = m2[j][5] + m2[j][7];
			m1[j][6] = m2[j][4] - m2[j][6];
			m1[j][7] = m2[j][5] - m2[j][7];

			m2[j][0] = m1[j][0] + m1[j][1];
			m2[j][1] = m1[j][0] - m1[j][1];
			m2[j][2] = m1[j][2] + m1[j][3];
			m2[j][3] = m1[j][2] - m1[j][3];
			m2[j][4] = m1[j][4] + m1[j][5];
			m2[j][5] = m1[j][4] - m1[j][5];
			m2[j][6] = m1[j][6] + m1[j][7];
			m2[j][7] = m1[j][6] - m1[j][7];
		}

		//vertical
		for (i = 0; i < 8; i++)
		{
			m1[0][i] = m2[0][i] + m2[2][i];
			m1[1][i] = m2[1][i] + m2[3][i];
			m1[2][i] = m2[0][i] - m2[2][i];
			m1[3][i] = m2[1][i] - m2[3][i];

			m2[0][i] = m1[0][i] + m1[1][i];
			m2[1][i] = m1[0][i] - m1[1][i];
			m2[2][i] = m1[2][i] + m1[3][i];
			m2[3][i] = m1[2][i] - m1[3][i];
		}

		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++)
			{
				sad += abs(m2[i][j]);
			}
		}

		sad = (Int)(sad / SAD_DIV_4x8);

		return sad;
	}
}

#define SAD_DIV_8x4 sqrt(4.0*8)*2
Distortion TComRdCost::xCalcHADs8x4(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	if (bitDepth == 8 || bitDepth == 10)
	{
		Int sad = 0;
		/* all 128 bit registers are named with a suffix mxnb, where m is the */
		/* number of n bits packed in the register                            */
		__m128i src0_8x16b, src1_8x16b, src2_8x16b, src3_8x16b;
		__m128i src4_8x16b, src5_8x16b, src6_8x16b, src7_8x16b;
		__m128i pred0_8x16b, pred1_8x16b, pred2_8x16b, pred3_8x16b;
		__m128i pred4_8x16b, pred5_8x16b, pred6_8x16b;
		__m128i out0_8x16b, out1_8x16b, out2_8x16b, out3_8x16b;
		__m128i out4_8x16b, out5_8x16b, out6_8x16b, out7_8x16b;

		/**********************Residue Calculation********************************/
		src0_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src1_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src2_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;
		src3_8x16b = _mm_loadu_si128((__m128i *) piOrg);
		piOrg = piOrg + iStrideOrg;

		pred0_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred1_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred2_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;
		pred3_8x16b = _mm_loadu_si128((__m128i *) piCur);
		piCur = piCur + iStrideCur;

		src0_8x16b = _mm_sub_epi16(src0_8x16b, pred0_8x16b);
		src1_8x16b = _mm_sub_epi16(src1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_sub_epi16(src2_8x16b, pred2_8x16b);
		src3_8x16b = _mm_sub_epi16(src3_8x16b, pred3_8x16b);
		/**********************Residue Calculation********************************/

		/**************** 8x4 horizontal transform *******************************/
		/***********************    8x4 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out7_8x16b = _mm_unpackhi_epi16(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_unpackhi_epi16(src2_8x16b, src3_8x16b);

		out0_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out1_8x16b = _mm_srli_si128(out0_8x16b, 8);
		out2_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_srli_si128(out2_8x16b, 8);
		out4_8x16b = _mm_unpacklo_epi32(out7_8x16b, src2_8x16b);
		out5_8x16b = _mm_srli_si128(out4_8x16b, 8);
		out6_8x16b = _mm_unpackhi_epi32(out7_8x16b, src2_8x16b);
		out7_8x16b = _mm_srli_si128(out6_8x16b, 8);
		/**********************   8x4 16 bit Transpose End   *********************/

		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(out0_8x16b, out1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(out2_8x16b, out3_8x16b);
		/* r4 + r5 */
		pred4_8x16b = _mm_add_epi16(out4_8x16b, out5_8x16b);
		/* r6 + r7 */
		pred6_8x16b = _mm_add_epi16(out6_8x16b, out7_8x16b);

		/* r0 + r1 + r2 + r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 + r6 + r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 */
		src0_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 + r2 + r3 - r4 - r5 - r6 - r7 */
		src4_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 + r1 - r2 - r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 + r5 - r6 - r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 + r1 - r2 - r3 + r4 + r5 - r6 - r7 */
		src2_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 + r1 - r2 - r3 - r4 - r5 + r6 + r7 */
		src6_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(out0_8x16b, out1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(out2_8x16b, out3_8x16b);
		/* r4 - r5 */
		pred4_8x16b = _mm_sub_epi16(out4_8x16b, out5_8x16b);
		/* r6 - r7 */
		pred6_8x16b = _mm_sub_epi16(out6_8x16b, out7_8x16b);

		/* r0 - r1 + r2 - r3 */
		pred1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 + r6 - r7 */
		pred5_8x16b = _mm_add_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 + r2 - r3 + r4 - r5 + r6 - r7 */
		src1_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 + r2 - r3 - r4 + r5 - r6 + r7 */
		src5_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/* r0 - r1 - r2 + r3 */
		pred1_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/* r4 - r5 - r6 + r7 */
		pred5_8x16b = _mm_sub_epi16(pred4_8x16b, pred6_8x16b);
		/* r0 - r1 - r2 + r3 + r4 - r5 - r6 + r7 */
		src3_8x16b = _mm_add_epi16(pred1_8x16b, pred5_8x16b);
		/* r0 - r1 - r2 + r3 - r4 + r5 + r6 - r7 */
		src7_8x16b = _mm_sub_epi16(pred1_8x16b, pred5_8x16b);

		/***********************    8x4 16 bit Transpose  ************************/
		out3_8x16b = _mm_unpacklo_epi16(src0_8x16b, src1_8x16b);
		pred0_8x16b = _mm_unpacklo_epi16(src2_8x16b, src3_8x16b);
		out2_8x16b = _mm_unpacklo_epi16(src4_8x16b, src5_8x16b);
		pred3_8x16b = _mm_unpacklo_epi16(src6_8x16b, src7_8x16b);

		out1_8x16b = _mm_unpacklo_epi32(out3_8x16b, pred0_8x16b);
		out3_8x16b = _mm_unpackhi_epi32(out3_8x16b, pred0_8x16b);
		pred1_8x16b = _mm_unpacklo_epi32(out2_8x16b, pred3_8x16b);
		pred3_8x16b = _mm_unpackhi_epi32(out2_8x16b, pred3_8x16b);

		src0_8x16b = _mm_unpacklo_epi64(out1_8x16b, pred1_8x16b);
		src1_8x16b = _mm_unpackhi_epi64(out1_8x16b, pred1_8x16b);
		src2_8x16b = _mm_unpacklo_epi64(out3_8x16b, pred3_8x16b);
		src3_8x16b = _mm_unpackhi_epi64(out3_8x16b, pred3_8x16b);
		/**********************   8x4 16 bit Transpose End   *********************/
		/**************** 8x4 horizontal transform *******************************/

		/************************* 8x4 Vertical Transform*************************/
		/* r0 + r1 */
		pred0_8x16b = _mm_add_epi16(src0_8x16b, src1_8x16b);
		/* r2 + r3 */
		pred2_8x16b = _mm_add_epi16(src2_8x16b, src3_8x16b);

		/* r0 + r1 + r2 + r3 */
		out0_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 + r1 - r2 - r3 */
		out2_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);

		/* r0 - r1 */
		pred0_8x16b = _mm_sub_epi16(src0_8x16b, src1_8x16b);
		/* r2 - r3 */
		pred2_8x16b = _mm_sub_epi16(src2_8x16b, src3_8x16b);

		/* r0 - r1 + r2 - r3 */
		out1_8x16b = _mm_add_epi16(pred0_8x16b, pred2_8x16b);
		/* r0 - r1 - r2 + r3 */
		out3_8x16b = _mm_sub_epi16(pred0_8x16b, pred2_8x16b);
		/************************* 8x4 Vertical Transform*************************/

		/****************************SATD calculation ****************************/
		src0_8x16b = _mm_abs_epi16(out0_8x16b);
		src1_8x16b = _mm_abs_epi16(out1_8x16b);
		src2_8x16b = _mm_abs_epi16(out2_8x16b);
		src3_8x16b = _mm_abs_epi16(out3_8x16b);

		src4_8x16b = _mm_srli_si128(src0_8x16b, 8);
		src5_8x16b = _mm_srli_si128(src1_8x16b, 8);
		src6_8x16b = _mm_srli_si128(src2_8x16b, 8);
		src7_8x16b = _mm_srli_si128(src3_8x16b, 8);

		src0_8x16b = _mm_cvtepu16_epi32(src0_8x16b);
		src1_8x16b = _mm_cvtepu16_epi32(src1_8x16b);
		src2_8x16b = _mm_cvtepu16_epi32(src2_8x16b);
		src3_8x16b = _mm_cvtepu16_epi32(src3_8x16b);
		src4_8x16b = _mm_cvtepu16_epi32(src4_8x16b);
		src5_8x16b = _mm_cvtepu16_epi32(src5_8x16b);
		src6_8x16b = _mm_cvtepu16_epi32(src6_8x16b);
		src7_8x16b = _mm_cvtepu16_epi32(src7_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src1_8x16b);
		src2_8x16b = _mm_add_epi32(src2_8x16b, src3_8x16b);
		src4_8x16b = _mm_add_epi32(src4_8x16b, src5_8x16b);
		src6_8x16b = _mm_add_epi32(src6_8x16b, src7_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src2_8x16b);
		src4_8x16b = _mm_add_epi32(src4_8x16b, src6_8x16b);

		src0_8x16b = _mm_add_epi32(src0_8x16b, src4_8x16b);

		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);
		src0_8x16b = _mm_hadd_epi32(src0_8x16b, src0_8x16b);

		sad += _mm_cvtsi128_si32(src0_8x16b);

		sad = (Int)(sad / SAD_DIV_8x4);

		return sad;
	}
	else
	{
		Int k, i, j, jj, sad = 0;
		Int diff[32], m1[4][8], m2[4][8];
		for (k = 0; k < 32; k += 8)
		{
			diff[k + 0] = piOrg[0] - piCur[0];
			diff[k + 1] = piOrg[1] - piCur[1];
			diff[k + 2] = piOrg[2] - piCur[2];
			diff[k + 3] = piOrg[3] - piCur[3];
			diff[k + 4] = piOrg[4] - piCur[4];
			diff[k + 5] = piOrg[5] - piCur[5];
			diff[k + 6] = piOrg[6] - piCur[6];
			diff[k + 7] = piOrg[7] - piCur[7];

			piCur += iStrideCur;
			piOrg += iStrideOrg;
		}

		//horizontal
		for (j = 0; j < 4; j++)
		{
			jj = j << 3;

			m2[j][0] = diff[jj] + diff[jj + 4];
			m2[j][1] = diff[jj + 1] + diff[jj + 5];
			m2[j][2] = diff[jj + 2] + diff[jj + 6];
			m2[j][3] = diff[jj + 3] + diff[jj + 7];
			m2[j][4] = diff[jj] - diff[jj + 4];
			m2[j][5] = diff[jj + 1] - diff[jj + 5];
			m2[j][6] = diff[jj + 2] - diff[jj + 6];
			m2[j][7] = diff[jj + 3] - diff[jj + 7];

			m1[j][0] = m2[j][0] + m2[j][2];
			m1[j][1] = m2[j][1] + m2[j][3];
			m1[j][2] = m2[j][0] - m2[j][2];
			m1[j][3] = m2[j][1] - m2[j][3];
			m1[j][4] = m2[j][4] + m2[j][6];
			m1[j][5] = m2[j][5] + m2[j][7];
			m1[j][6] = m2[j][4] - m2[j][6];
			m1[j][7] = m2[j][5] - m2[j][7];

			m2[j][0] = m1[j][0] + m1[j][1];
			m2[j][1] = m1[j][0] - m1[j][1];
			m2[j][2] = m1[j][2] + m1[j][3];
			m2[j][3] = m1[j][2] - m1[j][3];
			m2[j][4] = m1[j][4] + m1[j][5];
			m2[j][5] = m1[j][4] - m1[j][5];
			m2[j][6] = m1[j][6] + m1[j][7];
			m2[j][7] = m1[j][6] - m1[j][7];
		}

		//vertical
		for (i = 0; i < 8; i++)
		{
			m1[0][i] = m2[0][i] + m2[2][i];
			m1[1][i] = m2[1][i] + m2[3][i];
			m1[2][i] = m2[0][i] - m2[2][i];
			m1[3][i] = m2[1][i] - m2[3][i];

			m2[0][i] = m1[0][i] + m1[1][i];
			m2[1][i] = m1[0][i] - m1[1][i];
			m2[2][i] = m1[2][i] + m1[3][i];
			m2[3][i] = m1[2][i] - m1[3][i];
		}

		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++)
			{
				sad += abs(m2[i][j]);
			}
		}

		sad = (Int)(sad / SAD_DIV_8x4);

		return sad;
	}
}
#endif

#elif ITH_C3_HSADCOMPUTER==COARSE_SIMD_FXN

Distortion TComRdCost::xCalcHADs4x4(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth)
{
	Int k;
	Distortion satd = 0;
	TCoeff diff[16], m[16], d[16];

	assert(iStep == 1);
	for (k = 0; k < 16; k += 4)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	/*===== hadamard transform =====*/
	m[0] = diff[0] + diff[12];
	m[1] = diff[1] + diff[13];
	m[2] = diff[2] + diff[14];
	m[3] = diff[3] + diff[15];
	m[4] = diff[4] + diff[8];
	m[5] = diff[5] + diff[9];
	m[6] = diff[6] + diff[10];
	m[7] = diff[7] + diff[11];
	m[8] = diff[4] - diff[8];
	m[9] = diff[5] - diff[9];
	m[10] = diff[6] - diff[10];
	m[11] = diff[7] - diff[11];
	m[12] = diff[0] - diff[12];
	m[13] = diff[1] - diff[13];
	m[14] = diff[2] - diff[14];
	m[15] = diff[3] - diff[15];

	d[0] = m[0] + m[4];
	d[1] = m[1] + m[5];
	d[2] = m[2] + m[6];
	d[3] = m[3] + m[7];
	d[4] = m[8] + m[12];
	d[5] = m[9] + m[13];
	d[6] = m[10] + m[14];
	d[7] = m[11] + m[15];
	d[8] = m[0] - m[4];
	d[9] = m[1] - m[5];
	d[10] = m[2] - m[6];
	d[11] = m[3] - m[7];
	d[12] = m[12] - m[8];
	d[13] = m[13] - m[9];
	d[14] = m[14] - m[10];
	d[15] = m[15] - m[11];

	m[0] = d[0] + d[3];
	m[1] = d[1] + d[2];
	m[2] = d[1] - d[2];
	m[3] = d[0] - d[3];
	m[4] = d[4] + d[7];
	m[5] = d[5] + d[6];
	m[6] = d[5] - d[6];
	m[7] = d[4] - d[7];
	m[8] = d[8] + d[11];
	m[9] = d[9] + d[10];
	m[10] = d[9] - d[10];
	m[11] = d[8] - d[11];
	m[12] = d[12] + d[15];
	m[13] = d[13] + d[14];
	m[14] = d[13] - d[14];
	m[15] = d[12] - d[15];

	d[0] = m[0] + m[1];
	d[1] = m[0] - m[1];
	d[2] = m[2] + m[3];
	d[3] = m[3] - m[2];
	d[4] = m[4] + m[5];
	d[5] = m[4] - m[5];
	d[6] = m[6] + m[7];
	d[7] = m[7] - m[6];
	d[8] = m[8] + m[9];
	d[9] = m[8] - m[9];
	d[10] = m[10] + m[11];
	d[11] = m[11] - m[10];
	d[12] = m[12] + m[13];
	d[13] = m[12] - m[13];
	d[14] = m[14] + m[15];
	d[15] = m[15] - m[14];

	for (k = 0; k < 16; ++k)
	{
		satd += abs(d[k]);
	}
	satd = ((satd + 1) >> 1);

	return satd;
}

Distortion TComRdCost::xCalcHADs8x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth)
{
	if (bitDepth <= 10)
	{
		return(simdHADs8x8(piOrg, piCur, iStrideOrg, iStrideCur));
	}
	Int k, i, j, jj;
	Distortion sad = 0;
	TCoeff diff[64], m1[8][8], m2[8][8], m3[8][8];
	assert(iStep == 1);
	for (k = 0; k < 64; k += 8)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];
		diff[k + 4] = piOrg[4] - piCur[4];
		diff[k + 5] = piOrg[5] - piCur[5];
		diff[k + 6] = piOrg[6] - piCur[6];
		diff[k + 7] = piOrg[7] - piCur[7];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 8; j++)
	{
		jj = j << 3;
		m2[j][0] = diff[jj] + diff[jj + 4];
		m2[j][1] = diff[jj + 1] + diff[jj + 5];
		m2[j][2] = diff[jj + 2] + diff[jj + 6];
		m2[j][3] = diff[jj + 3] + diff[jj + 7];
		m2[j][4] = diff[jj] - diff[jj + 4];
		m2[j][5] = diff[jj + 1] - diff[jj + 5];
		m2[j][6] = diff[jj + 2] - diff[jj + 6];
		m2[j][7] = diff[jj + 3] - diff[jj + 7];

		m1[j][0] = m2[j][0] + m2[j][2];
		m1[j][1] = m2[j][1] + m2[j][3];
		m1[j][2] = m2[j][0] - m2[j][2];
		m1[j][3] = m2[j][1] - m2[j][3];
		m1[j][4] = m2[j][4] + m2[j][6];
		m1[j][5] = m2[j][5] + m2[j][7];
		m1[j][6] = m2[j][4] - m2[j][6];
		m1[j][7] = m2[j][5] - m2[j][7];

		m2[j][0] = m1[j][0] + m1[j][1];
		m2[j][1] = m1[j][0] - m1[j][1];
		m2[j][2] = m1[j][2] + m1[j][3];
		m2[j][3] = m1[j][2] - m1[j][3];
		m2[j][4] = m1[j][4] + m1[j][5];
		m2[j][5] = m1[j][4] - m1[j][5];
		m2[j][6] = m1[j][6] + m1[j][7];
		m2[j][7] = m1[j][6] - m1[j][7];
	}

	//vertical
	for (i = 0; i < 8; i++)
	{
		m3[0][i] = m2[0][i] + m2[4][i];
		m3[1][i] = m2[1][i] + m2[5][i];
		m3[2][i] = m2[2][i] + m2[6][i];
		m3[3][i] = m2[3][i] + m2[7][i];
		m3[4][i] = m2[0][i] - m2[4][i];
		m3[5][i] = m2[1][i] - m2[5][i];
		m3[6][i] = m2[2][i] - m2[6][i];
		m3[7][i] = m2[3][i] - m2[7][i];

		m1[0][i] = m3[0][i] + m3[2][i];
		m1[1][i] = m3[1][i] + m3[3][i];
		m1[2][i] = m3[0][i] - m3[2][i];
		m1[3][i] = m3[1][i] - m3[3][i];
		m1[4][i] = m3[4][i] + m3[6][i];
		m1[5][i] = m3[5][i] + m3[7][i];
		m1[6][i] = m3[4][i] - m3[6][i];
		m1[7][i] = m3[5][i] - m3[7][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
		m2[4][i] = m1[4][i] + m1[5][i];
		m2[5][i] = m1[4][i] - m1[5][i];
		m2[6][i] = m1[6][i] + m1[7][i];
		m2[7][i] = m1[6][i] - m1[7][i];
	}

	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			sad += abs(m2[i][j]);
		}
	}

	sad = ((sad + 2) >> 2);

	return sad;
}

#if JVET_C0024_QTBT
Distortion TComRdCost::xCalcHADs16x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{   //need to add SIMD implementation ,JCA
	Int k, i, j, jj, sad = 0;
	Int diff[128], m1[8][16], m2[8][16];
	for (k = 0; k < 128; k += 16)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];
		diff[k + 4] = piOrg[4] - piCur[4];
		diff[k + 5] = piOrg[5] - piCur[5];
		diff[k + 6] = piOrg[6] - piCur[6];
		diff[k + 7] = piOrg[7] - piCur[7];

		diff[k + 8] = piOrg[8] - piCur[8];
		diff[k + 9] = piOrg[9] - piCur[9];
		diff[k + 10] = piOrg[10] - piCur[10];
		diff[k + 11] = piOrg[11] - piCur[11];
		diff[k + 12] = piOrg[12] - piCur[12];
		diff[k + 13] = piOrg[13] - piCur[13];
		diff[k + 14] = piOrg[14] - piCur[14];
		diff[k + 15] = piOrg[15] - piCur[15];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 8; j++)
	{
		jj = j << 4;

		m2[j][0] = diff[jj] + diff[jj + 8];
		m2[j][1] = diff[jj + 1] + diff[jj + 9];
		m2[j][2] = diff[jj + 2] + diff[jj + 10];
		m2[j][3] = diff[jj + 3] + diff[jj + 11];
		m2[j][4] = diff[jj + 4] + diff[jj + 12];
		m2[j][5] = diff[jj + 5] + diff[jj + 13];
		m2[j][6] = diff[jj + 6] + diff[jj + 14];
		m2[j][7] = diff[jj + 7] + diff[jj + 15];
		m2[j][8] = diff[jj] - diff[jj + 8];
		m2[j][9] = diff[jj + 1] - diff[jj + 9];
		m2[j][10] = diff[jj + 2] - diff[jj + 10];
		m2[j][11] = diff[jj + 3] - diff[jj + 11];
		m2[j][12] = diff[jj + 4] - diff[jj + 12];
		m2[j][13] = diff[jj + 5] - diff[jj + 13];
		m2[j][14] = diff[jj + 6] - diff[jj + 14];
		m2[j][15] = diff[jj + 7] - diff[jj + 15];

		m1[j][0] = m2[j][0] + m2[j][4];
		m1[j][1] = m2[j][1] + m2[j][5];
		m1[j][2] = m2[j][2] + m2[j][6];
		m1[j][3] = m2[j][3] + m2[j][7];
		m1[j][4] = m2[j][0] - m2[j][4];
		m1[j][5] = m2[j][1] - m2[j][5];
		m1[j][6] = m2[j][2] - m2[j][6];
		m1[j][7] = m2[j][3] - m2[j][7];
		m1[j][8] = m2[j][8] + m2[j][12];
		m1[j][9] = m2[j][9] + m2[j][13];
		m1[j][10] = m2[j][10] + m2[j][14];
		m1[j][11] = m2[j][11] + m2[j][15];
		m1[j][12] = m2[j][8] - m2[j][12];
		m1[j][13] = m2[j][9] - m2[j][13];
		m1[j][14] = m2[j][10] - m2[j][14];
		m1[j][15] = m2[j][11] - m2[j][15];

		m2[j][0] = m1[j][0] + m1[j][2];
		m2[j][1] = m1[j][1] + m1[j][3];
		m2[j][2] = m1[j][0] - m1[j][2];
		m2[j][3] = m1[j][1] - m1[j][3];
		m2[j][4] = m1[j][4] + m1[j][6];
		m2[j][5] = m1[j][5] + m1[j][7];
		m2[j][6] = m1[j][4] - m1[j][6];
		m2[j][7] = m1[j][5] - m1[j][7];
		m2[j][8] = m1[j][8] + m1[j][10];
		m2[j][9] = m1[j][9] + m1[j][11];
		m2[j][10] = m1[j][8] - m1[j][10];
		m2[j][11] = m1[j][9] - m1[j][11];
		m2[j][12] = m1[j][12] + m1[j][14];
		m2[j][13] = m1[j][13] + m1[j][15];
		m2[j][14] = m1[j][12] - m1[j][14];
		m2[j][15] = m1[j][13] - m1[j][15];

		m1[j][0] = m2[j][0] + m2[j][1];
		m1[j][1] = m2[j][0] - m2[j][1];
		m1[j][2] = m2[j][2] + m2[j][3];
		m1[j][3] = m2[j][2] - m2[j][3];
		m1[j][4] = m2[j][4] + m2[j][5];
		m1[j][5] = m2[j][4] - m2[j][5];
		m1[j][6] = m2[j][6] + m2[j][7];
		m1[j][7] = m2[j][6] - m2[j][7];
		m1[j][8] = m2[j][8] + m2[j][9];
		m1[j][9] = m2[j][8] - m2[j][9];
		m1[j][10] = m2[j][10] + m2[j][11];
		m1[j][11] = m2[j][10] - m2[j][11];
		m1[j][12] = m2[j][12] + m2[j][13];
		m1[j][13] = m2[j][12] - m2[j][13];
		m1[j][14] = m2[j][14] + m2[j][15];
		m1[j][15] = m2[j][14] - m2[j][15];
	}

	//vertical
	for (i = 0; i < 16; i++)
	{
		m2[0][i] = m1[0][i] + m1[4][i];
		m2[1][i] = m1[1][i] + m1[5][i];
		m2[2][i] = m1[2][i] + m1[6][i];
		m2[3][i] = m1[3][i] + m1[7][i];
		m2[4][i] = m1[0][i] - m1[4][i];
		m2[5][i] = m1[1][i] - m1[5][i];
		m2[6][i] = m1[2][i] - m1[6][i];
		m2[7][i] = m1[3][i] - m1[7][i];

		m1[0][i] = m2[0][i] + m2[2][i];
		m1[1][i] = m2[1][i] + m2[3][i];
		m1[2][i] = m2[0][i] - m2[2][i];
		m1[3][i] = m2[1][i] - m2[3][i];
		m1[4][i] = m2[4][i] + m2[6][i];
		m1[5][i] = m2[5][i] + m2[7][i];
		m1[6][i] = m2[4][i] - m2[6][i];
		m1[7][i] = m2[5][i] - m2[7][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
		m2[4][i] = m1[4][i] + m1[5][i];
		m2[5][i] = m1[4][i] - m1[5][i];
		m2[6][i] = m1[6][i] + m1[7][i];
		m2[7][i] = m1[6][i] - m1[7][i];
	}

	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 16; j++)
		{
			sad += abs(m2[i][j]);
		}
	}

	sad = (Int)(sad / sqrt(16.0 * 8) * 2);

	return sad;
}

Distortion TComRdCost::xCalcHADs8x16(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	Int k, i, j, jj, sad = 0;
	Int diff[128], m1[16][8], m2[16][8];
	for (k = 0; k < 128; k += 8)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];
		diff[k + 4] = piOrg[4] - piCur[4];
		diff[k + 5] = piOrg[5] - piCur[5];
		diff[k + 6] = piOrg[6] - piCur[6];
		diff[k + 7] = piOrg[7] - piCur[7];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 16; j++)
	{
		jj = j << 3;

		m2[j][0] = diff[jj] + diff[jj + 4];
		m2[j][1] = diff[jj + 1] + diff[jj + 5];
		m2[j][2] = diff[jj + 2] + diff[jj + 6];
		m2[j][3] = diff[jj + 3] + diff[jj + 7];
		m2[j][4] = diff[jj] - diff[jj + 4];
		m2[j][5] = diff[jj + 1] - diff[jj + 5];
		m2[j][6] = diff[jj + 2] - diff[jj + 6];
		m2[j][7] = diff[jj + 3] - diff[jj + 7];

		m1[j][0] = m2[j][0] + m2[j][2];
		m1[j][1] = m2[j][1] + m2[j][3];
		m1[j][2] = m2[j][0] - m2[j][2];
		m1[j][3] = m2[j][1] - m2[j][3];
		m1[j][4] = m2[j][4] + m2[j][6];
		m1[j][5] = m2[j][5] + m2[j][7];
		m1[j][6] = m2[j][4] - m2[j][6];
		m1[j][7] = m2[j][5] - m2[j][7];

		m2[j][0] = m1[j][0] + m1[j][1];
		m2[j][1] = m1[j][0] - m1[j][1];
		m2[j][2] = m1[j][2] + m1[j][3];
		m2[j][3] = m1[j][2] - m1[j][3];
		m2[j][4] = m1[j][4] + m1[j][5];
		m2[j][5] = m1[j][4] - m1[j][5];
		m2[j][6] = m1[j][6] + m1[j][7];
		m2[j][7] = m1[j][6] - m1[j][7];
	}

	//vertical
	for (i = 0; i < 8; i++)
	{
		m1[0][i] = m2[0][i] + m2[8][i];
		m1[1][i] = m2[1][i] + m2[9][i];
		m1[2][i] = m2[2][i] + m2[10][i];
		m1[3][i] = m2[3][i] + m2[11][i];
		m1[4][i] = m2[4][i] + m2[12][i];
		m1[5][i] = m2[5][i] + m2[13][i];
		m1[6][i] = m2[6][i] + m2[14][i];
		m1[7][i] = m2[7][i] + m2[15][i];
		m1[8][i] = m2[0][i] - m2[8][i];
		m1[9][i] = m2[1][i] - m2[9][i];
		m1[10][i] = m2[2][i] - m2[10][i];
		m1[11][i] = m2[3][i] - m2[11][i];
		m1[12][i] = m2[4][i] - m2[12][i];
		m1[13][i] = m2[5][i] - m2[13][i];
		m1[14][i] = m2[6][i] - m2[14][i];
		m1[15][i] = m2[7][i] - m2[15][i];

		m2[0][i] = m1[0][i] + m1[4][i];
		m2[1][i] = m1[1][i] + m1[5][i];
		m2[2][i] = m1[2][i] + m1[6][i];
		m2[3][i] = m1[3][i] + m1[7][i];
		m2[4][i] = m1[0][i] - m1[4][i];
		m2[5][i] = m1[1][i] - m1[5][i];
		m2[6][i] = m1[2][i] - m1[6][i];
		m2[7][i] = m1[3][i] - m1[7][i];
		m2[8][i] = m1[8][i] + m1[12][i];
		m2[9][i] = m1[9][i] + m1[13][i];
		m2[10][i] = m1[10][i] + m1[14][i];
		m2[11][i] = m1[11][i] + m1[15][i];
		m2[12][i] = m1[8][i] - m1[12][i];
		m2[13][i] = m1[9][i] - m1[13][i];
		m2[14][i] = m1[10][i] - m1[14][i];
		m2[15][i] = m1[11][i] - m1[15][i];

		m1[0][i] = m2[0][i] + m2[2][i];
		m1[1][i] = m2[1][i] + m2[3][i];
		m1[2][i] = m2[0][i] - m2[2][i];
		m1[3][i] = m2[1][i] - m2[3][i];
		m1[4][i] = m2[4][i] + m2[6][i];
		m1[5][i] = m2[5][i] + m2[7][i];
		m1[6][i] = m2[4][i] - m2[6][i];
		m1[7][i] = m2[5][i] - m2[7][i];
		m1[8][i] = m2[8][i] + m2[10][i];
		m1[9][i] = m2[9][i] + m2[11][i];
		m1[10][i] = m2[8][i] - m2[10][i];
		m1[11][i] = m2[9][i] - m2[11][i];
		m1[12][i] = m2[12][i] + m2[14][i];
		m1[13][i] = m2[13][i] + m2[15][i];
		m1[14][i] = m2[12][i] - m2[14][i];
		m1[15][i] = m2[13][i] - m2[15][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
		m2[4][i] = m1[4][i] + m1[5][i];
		m2[5][i] = m1[4][i] - m1[5][i];
		m2[6][i] = m1[6][i] + m1[7][i];
		m2[7][i] = m1[6][i] - m1[7][i];
		m2[8][i] = m1[8][i] + m1[9][i];
		m2[9][i] = m1[8][i] - m1[9][i];
		m2[10][i] = m1[10][i] + m1[11][i];
		m2[11][i] = m1[10][i] - m1[11][i];
		m2[12][i] = m1[12][i] + m1[13][i];
		m2[13][i] = m1[12][i] - m1[13][i];
		m2[14][i] = m1[14][i] + m1[15][i];
		m2[15][i] = m1[14][i] - m1[15][i];
	}

	for (i = 0; i < 16; i++)
	{
		for (j = 0; j < 8; j++)
		{
			sad += abs(m2[i][j]);
		}
	}

	sad = (Int)(sad / sqrt(16.0 * 8) * 2);

	return sad;
}

Distortion TComRdCost::xCalcHADs4x8(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	Int k, i, j, jj, sad = 0;
	Int diff[32], m1[8][4], m2[8][4];
	for (k = 0; k < 32; k += 4)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 8; j++)
	{
		jj = j << 2;
		m2[j][0] = diff[jj] + diff[jj + 2];
		m2[j][1] = diff[jj + 1] + diff[jj + 3];
		m2[j][2] = diff[jj] - diff[jj + 2];
		m2[j][3] = diff[jj + 1] - diff[jj + 3];

		m1[j][0] = m2[j][0] + m2[j][1];
		m1[j][1] = m2[j][0] - m2[j][1];
		m1[j][2] = m2[j][2] + m2[j][3];
		m1[j][3] = m2[j][2] - m2[j][3];
	}

	//vertical
	for (i = 0; i < 4; i++)
	{
		m2[0][i] = m1[0][i] + m1[4][i];
		m2[1][i] = m1[1][i] + m1[5][i];
		m2[2][i] = m1[2][i] + m1[6][i];
		m2[3][i] = m1[3][i] + m1[7][i];
		m2[4][i] = m1[0][i] - m1[4][i];
		m2[5][i] = m1[1][i] - m1[5][i];
		m2[6][i] = m1[2][i] - m1[6][i];
		m2[7][i] = m1[3][i] - m1[7][i];

		m1[0][i] = m2[0][i] + m2[2][i];
		m1[1][i] = m2[1][i] + m2[3][i];
		m1[2][i] = m2[0][i] - m2[2][i];
		m1[3][i] = m2[1][i] - m2[3][i];
		m1[4][i] = m2[4][i] + m2[6][i];
		m1[5][i] = m2[5][i] + m2[7][i];
		m1[6][i] = m2[4][i] - m2[6][i];
		m1[7][i] = m2[5][i] - m2[7][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
		m2[4][i] = m1[4][i] + m1[5][i];
		m2[5][i] = m1[4][i] - m1[5][i];
		m2[6][i] = m1[6][i] + m1[7][i];
		m2[7][i] = m1[6][i] - m1[7][i];
	}

	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 4; j++)
		{
			sad += abs(m2[i][j]);
		}
	}

	sad = (Int)(sad / sqrt(4.0 * 8) * 2);

	return sad;
}

Distortion TComRdCost::xCalcHADs8x4(Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth)
{
	Int k, i, j, jj, sad = 0;
	Int diff[32], m1[4][8], m2[4][8];
	for (k = 0; k < 32; k += 8)
	{
		diff[k + 0] = piOrg[0] - piCur[0];
		diff[k + 1] = piOrg[1] - piCur[1];
		diff[k + 2] = piOrg[2] - piCur[2];
		diff[k + 3] = piOrg[3] - piCur[3];
		diff[k + 4] = piOrg[4] - piCur[4];
		diff[k + 5] = piOrg[5] - piCur[5];
		diff[k + 6] = piOrg[6] - piCur[6];
		diff[k + 7] = piOrg[7] - piCur[7];

		piCur += iStrideCur;
		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 4; j++)
	{
		jj = j << 3;

		m2[j][0] = diff[jj] + diff[jj + 4];
		m2[j][1] = diff[jj + 1] + diff[jj + 5];
		m2[j][2] = diff[jj + 2] + diff[jj + 6];
		m2[j][3] = diff[jj + 3] + diff[jj + 7];
		m2[j][4] = diff[jj] - diff[jj + 4];
		m2[j][5] = diff[jj + 1] - diff[jj + 5];
		m2[j][6] = diff[jj + 2] - diff[jj + 6];
		m2[j][7] = diff[jj + 3] - diff[jj + 7];

		m1[j][0] = m2[j][0] + m2[j][2];
		m1[j][1] = m2[j][1] + m2[j][3];
		m1[j][2] = m2[j][0] - m2[j][2];
		m1[j][3] = m2[j][1] - m2[j][3];
		m1[j][4] = m2[j][4] + m2[j][6];
		m1[j][5] = m2[j][5] + m2[j][7];
		m1[j][6] = m2[j][4] - m2[j][6];
		m1[j][7] = m2[j][5] - m2[j][7];

		m2[j][0] = m1[j][0] + m1[j][1];
		m2[j][1] = m1[j][0] - m1[j][1];
		m2[j][2] = m1[j][2] + m1[j][3];
		m2[j][3] = m1[j][2] - m1[j][3];
		m2[j][4] = m1[j][4] + m1[j][5];
		m2[j][5] = m1[j][4] - m1[j][5];
		m2[j][6] = m1[j][6] + m1[j][7];
		m2[j][7] = m1[j][6] - m1[j][7];
	}

	//vertical
	for (i = 0; i < 8; i++)
	{
		m1[0][i] = m2[0][i] + m2[2][i];
		m1[1][i] = m2[1][i] + m2[3][i];
		m1[2][i] = m2[0][i] - m2[2][i];
		m1[3][i] = m2[1][i] - m2[3][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
	}

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 8; j++)
		{
			sad += abs(m2[i][j]);
		}
	}

	sad = (Int)(sad / sqrt(4.0 * 8) * 2);

	return sad;
}
#endif
#endif
Distortion TComRdCost::xGetHADs( DistParam* pcDtParam )
{
#if VCEG_AZ06_IC
  if( pcDtParam->bMRFlag )
  {
    return xGetMRHADs( pcDtParam );
  }
#endif
  if ( pcDtParam->bApplyWeight )
  {
    return TComRdCostWeightPrediction::xGetHADsw( pcDtParam );
  }

  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideCur = pcDtParam->iStrideCur;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStep  = pcDtParam->iStep;

  Int  x, y;

  Distortion uiSum = 0;

#if JVET_C0024_QTBT
  if ( iCols > iRows && iRows>=8) 
  {
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 16 )
      {
        uiSum += xCalcHADs16x8( &piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*8;
      piCur += iStrideCur*8;
    }
  }
  else if (  iCols < iRows && iCols>=8) 
  {
    for ( y=0; y<iRows; y+= 16 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x16(&piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*16;
      piCur += iStrideCur*16;
    }
  }
  else if ( iCols > iRows && iRows==4) 
  {
    for ( y=0; y<iRows; y+= 4 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x4( &piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*4;
      piCur += iStrideCur*4;
    }
  }
  else if (  iCols < iRows && iCols==4) 
  {
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 4 )
      {
        uiSum += xCalcHADs4x8(&piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*8;
      piCur += iStrideCur*8;
    }
  }
  else
#endif
  if( ( iRows % 8 == 0) && (iCols % 8 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<3;
    Int  iOffsetCur = iStrideCur<<3;
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x8( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep 
#if COM16_C806_SIMD_OPT
          , pcDtParam->bitDepth
#endif
          );
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else if( ( iRows % 4 == 0) && (iCols % 4 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<2;
    Int  iOffsetCur = iStrideCur<<2;

    for ( y=0; y<iRows; y+= 4 )
    {
      for ( x=0; x<iCols; x+= 4 )
      {
        uiSum += xCalcHADs4x4( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep, pcDtParam->bitDepth);
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else if( ( iRows % 2 == 0) && (iCols % 2 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<1;
    Int  iOffsetCur = iStrideCur<<1;
    for ( y=0; y<iRows; y+=2 )
    {
      for ( x=0; x<iCols; x+=2 )
      {
        uiSum += xCalcHADs2x2( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep );
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else
  {
    assert(false);
  }

  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT(pcDtParam->bitDepth-8) );
}

#if VCEG_AZ06_IC
UInt TComRdCost::xGetMRHADs( DistParam* pcDtParam )
{
  Pel* piOrg   = pcDtParam->pOrg;
  Pel* piCur   = pcDtParam->pCur;
  Int  iRows   = pcDtParam->iRows;
  Int  iCols   = pcDtParam->iCols;
  Int  iStrideCur = pcDtParam->iStrideCur;
  Int  iStrideOrg = pcDtParam->iStrideOrg;
  Int  iStep  = pcDtParam->iStep;

  Int  x, y;

  UInt uiSum = 0;
  Int  iOrigAvg = 0, iCurAvg = 0;
  Int  iDeltaC;

  for ( y=0; y<iRows; y++ )
  {
    for ( x=0; x<iCols; x++ )
    {        
      iOrigAvg += piOrg[x];
      iCurAvg  += piCur[x];
    }
    piOrg += iStrideOrg;
    piCur += iStrideCur;
  }

  piOrg   = pcDtParam->pOrg;
  piCur   = pcDtParam->pCur;

  iDeltaC = (iOrigAvg - iCurAvg)/iRows/iCols;

  for ( y=0; y<iRows; y++ )
  {
    for ( x=0; x<iCols; x++ )
    {        
      piOrg[x] -= iDeltaC;
    }
    piOrg += iStrideOrg;
  }

  piOrg   = pcDtParam->pOrg;

#if JVET_C0024_QTBT 
  if ( iCols > iRows && iRows>=8) 
  {
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 16 )
      {
        uiSum += xCalcHADs16x8( &piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*8;
      piCur += iStrideCur*8;
    }
  }
  else if (  iCols < iRows && iCols>=8) 
  {
    for ( y=0; y<iRows; y+= 16 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x16(&piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*16;
      piCur += iStrideCur*16;
    }
  }
  else if ( iCols > iRows && iRows==4) 
  {
    for ( y=0; y<iRows; y+= 4 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x4( &piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*4;
      piCur += iStrideCur*4;
    }
  }
  else if (  iCols < iRows && iCols==4) 
  {
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 4 )
      {
        uiSum += xCalcHADs4x8(&piOrg[x], &piCur[x], iStrideOrg, iStrideCur, pcDtParam->bitDepth);
      }
      piOrg += iStrideOrg*8;
      piCur += iStrideCur*8;
    }
  }
  else
#endif
  if( ( iRows % 8 == 0) && (iCols % 8 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<3;
    Int  iOffsetCur = iStrideCur<<3;
    for ( y=0; y<iRows; y+= 8 )
    {
      for ( x=0; x<iCols; x+= 8 )
      {
        uiSum += xCalcHADs8x8( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep 
#if COM16_C806_SIMD_OPT
          , pcDtParam->bitDepth
#endif
          );
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else if( ( iRows % 4 == 0) && (iCols % 4 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<2;
    Int  iOffsetCur = iStrideCur<<2;

    for ( y=0; y<iRows; y+= 4 )
    {
      for ( x=0; x<iCols; x+= 4 )
      {
        uiSum += xCalcHADs4x4( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep, pcDtParam->bitDepth);
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else if( ( iRows % 2 == 0) && (iCols % 2 == 0) )
  {
    Int  iOffsetOrg = iStrideOrg<<1;
    Int  iOffsetCur = iStrideCur<<1;
    for ( y=0; y<iRows; y+=2 )
    {
      for ( x=0; x<iCols; x+=2 )
      {
        uiSum += xCalcHADs2x2( &piOrg[x], &piCur[x*iStep], iStrideOrg, iStrideCur, iStep );
      }
      piOrg += iOffsetOrg;
      piCur += iOffsetCur;
    }
  }
  else
  {
    assert(false);
  }

  piOrg   = pcDtParam->pOrg;

  for ( y=0; y<iRows; y++ )
  {
    for ( x=0; x<iCols; x++ )
    {        
      piOrg[x] += iDeltaC;
    }
    piOrg += iStrideOrg;
  }

  return ( uiSum >> DISTORTION_PRECISION_ADJUSTMENT( pcDtParam->bitDepth - 8 ) );
}
#endif
//! \}
