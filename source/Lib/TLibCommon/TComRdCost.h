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

/** \file     TComRdCost.h
    \brief    RD cost computation classes (header)
*/

#ifndef __TCOMRDCOST__
#define __TCOMRDCOST__


#include "CommonDef.h"
#include "TComPattern.h"
#include "TComMv.h"

#include "TComSlice.h"
#include "TComRdCostWeightPrediction.h"

//! \ingroup TLibCommon
//! \{

class DistParam;
class TComPattern;

// ====================================================================================================================
// Type definition
// ====================================================================================================================

// for function pointer
typedef Distortion (*FpDistFunc) (DistParam*); // TODO: can this pointer be replaced with a reference? - there are no NULL checks on pointer.

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// distortion parameter class
class DistParam
{
public:
#if DMVR_MRSAD_SPEEDUP
  Int CalcStrategy;
  Int OrgMean;
  Int CurrMean;
#endif
  Pel*  pOrg;
  Pel*  pCur;
  Int   iStrideOrg;
  Int   iStrideCur;
#if WCG_LUMA_DQP_CM_SCALE
  Pel*  pOrgLuma;          // use luma to get weighting
  Int   iStrideOrgLuma; 
#endif

  Int   iRows;
  Int   iCols;
  Int   iStep;
#if VCEG_AZ06_IC
  Bool  bMRFlag;
#endif
  FpDistFunc DistFunc;
  Int   bitDepth;

  Bool            bApplyWeight;     // whether weighted prediction is used or not
  WPScalingParam  *wpCur;           // weighted prediction scaling parameters for current ref
  ComponentID     compIdx;

  // (vertical) subsampling shift (for reducing complexity)
  // - 0 = no subsampling, 1 = even rows, 2 = every 4th, etc.
  Int   iSubShift;

  DistParam()
  {
#if DMVR_MRSAD_SPEEDUP
    CalcStrategy = 0;
    OrgMean = 0;
    CurrMean = 0;
#endif
    pOrg = NULL;
    pCur = NULL;
    iStrideOrg = 0;
    iStrideCur = 0;
#if WCG_LUMA_DQP_CM_SCALE
    pOrgLuma = NULL;          // use luma to get weighting
    iStrideOrgLuma = 0; 
    compIdx = COMPONENT_Y;
#endif
    iRows = 0;
    iCols = 0;
    iStep = 1;
    DistFunc = NULL;
    iSubShift = 0;
    bitDepth = 0;
#if VCEG_AZ06_IC
    bMRFlag = false;
#endif
  }
};

/// RD cost computation class
class TComRdCost
{
private:
  // for distortion

  FpDistFunc              m_afpDistortFunc[DF_TOTAL_FUNCTIONS]; // [eDFunc]
  CostMode                m_costMode;
  Double                  m_distortionWeight[MAX_NUM_COMPONENT]; // only chroma values are used.
  Double                  m_dLambda;
#if WCG_LUMA_DQP_CM_SCALE
  Double                  m_dLambda_unadjusted;
  static Double           m_lumaLevelToWeightPLUT[LUMA_LEVEL_TO_DQP_LUT_MAXSIZE];
#endif
  Double                  m_sqrtLambda;
#if RExt__HIGH_BIT_DEPTH_SUPPORT
  Double                  m_dLambdaMotionSAD[2 /* 0=standard, 1=for transquant bypass when mixed-lossless cost evaluation enabled*/];
  Double                  m_dLambdaMotionSSE[2 /* 0=standard, 1=for transquant bypass when mixed-lossless cost evaluation enabled*/];
#else
  UInt                    m_uiLambdaMotionSAD[2 /* 0=standard, 1=for transquant bypass when mixed-lossless cost evaluation enabled*/];
  UInt                    m_uiLambdaMotionSSE[2 /* 0=standard, 1=for transquant bypass when mixed-lossless cost evaluation enabled*/];
#endif
  Double                  m_dFrameLambda;

  // for motion cost
  TComMv                  m_mvPredictor;
#if RExt__HIGH_BIT_DEPTH_SUPPORT
  Double                  m_dCost;
#else
  UInt                    m_uiCost;
#endif
  Int                     m_iCostScale;

public:
  TComRdCost();
  virtual ~TComRdCost();

  Double  calcRdCost  ( UInt   uiBits, Distortion uiDistortion, Bool bFlag = false, DFunc eDFunc = DF_DEFAULT );
  Double  calcRdCost64( UInt64 uiBits, UInt64 uiDistortion, Bool bFlag = false, DFunc eDFunc = DF_DEFAULT );

  Void    setDistortionWeight  ( const ComponentID compID, const Double distortionWeight ) { m_distortionWeight[compID] = distortionWeight; }
  Void    setLambda      ( Double dLambda, const BitDepths &bitDepths );
  Void    setFrameLambda ( Double dLambda ) { m_dFrameLambda = dLambda; }
#if WCG_LUMA_DQP_CM_SCALE
  Void    saveUnadjustedLambda() {m_dLambda_unadjusted = m_dLambda;}
  Void    initLumaLevelToWeightTable(Bool isSDR);
  static  Distortion getWeightedMSE(Int compIdx, Pel org, Pel cur, UInt uiShift, Pel orgLuma);
  inline  Double     getWPSNRLumaLevelWeight(Int val) { return m_lumaLevelToWeightPLUT[val]; }
#endif
  Double  getSqrtLambda ()   { return m_sqrtLambda; }

  Double  getLambda() { return m_dLambda; }
  Double  getChromaWeight () { return ((m_distortionWeight[COMPONENT_Cb] + m_distortionWeight[COMPONENT_Cr]) / 2.0); }

  Void    setCostMode(CostMode   m )    { m_costMode = m; }

  // Distortion Functions
  Void    init();

  Void    setDistParam( UInt uiBlkWidth, UInt uiBlkHeight, DFunc eDFunc, DistParam& rcDistParam );
  Void    setDistParam( TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride,            DistParam& rcDistParam );
  Void    setDistParam( TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, Int iStep, DistParam& rcDistParam, Bool bHADME=false );
  Void    setDistParam( DistParam& rcDP, Int bitDepth, Pel* p1, Int iStride1, Pel* p2, Int iStride2, Int iWidth, Int iHeight, Bool bHadamard = false );

  Distortion calcHAD(Int bitDepth, Pel* pi0, Int iStride0, Pel* pi1, Int iStride1, Int iWidth, Int iHeight );

  // for motion cost
  static UInt    xGetExpGolombNumberOfBits( Int iVal 
#if VCEG_AZ07_IMV
    , Int iMvFlag = false
#endif
    );
#if RExt__HIGH_BIT_DEPTH_SUPPORT
  Void    getMotionCost( Bool bSad, Int iAdd, Bool bIsTransquantBypass ) { m_dCost = (bSad ? m_dLambdaMotionSAD[(bIsTransquantBypass && m_costMode==COST_MIXED_LOSSLESS_LOSSY_CODING) ?1:0] + iAdd : m_dLambdaMotionSSE[(bIsTransquantBypass && m_costMode==COST_MIXED_LOSSLESS_LOSSY_CODING)?1:0] + iAdd); }
#else
  Void    getMotionCost( Bool bSad, Int iAdd, Bool bIsTransquantBypass ) { m_uiCost = (bSad ? m_uiLambdaMotionSAD[(bIsTransquantBypass && m_costMode==COST_MIXED_LOSSLESS_LOSSY_CODING) ?1:0] + iAdd : m_uiLambdaMotionSSE[(bIsTransquantBypass && m_costMode==COST_MIXED_LOSSLESS_LOSSY_CODING)?1:0] + iAdd); }
#endif
  Void    setPredictor( TComMv& rcMv )
  {
    m_mvPredictor = rcMv;
  }
  Void    setCostScale( Int iCostScale )    { m_iCostScale = iCostScale; }
  __inline Distortion getCost( Int x, Int y 
#if VCEG_AZ07_IMV
    , Int iMvFlag = false
#endif
    )
  {
#if RExt__HIGH_BIT_DEPTH_SUPPORT
    return Distortion((m_dCost * getBits(x, y
#if VCEG_AZ07_IMV
      , iMvFlag
#endif
      )) / 65536.0);
#else
    return m_uiCost * getBits(x, y
#if VCEG_AZ07_IMV
      , iMvFlag
#endif
      ) >> 16;
#endif
  }
#if JVET_D0123_ME_CTX_LUT_BITS
  Distortion getCostSearch (UInt uiBits)
  {
    return (Distortion) (m_sqrtLambda * uiBits / ((Double)EPBIT));
  }
#endif
#if RExt__HIGH_BIT_DEPTH_SUPPORT
  Distortion getCost( UInt b )                 { return Distortion(( m_dCost * b ) / 65536.0); }
#else
  Distortion getCost( UInt b )                 { return ( m_uiCost * b ) >> 16; }
#endif
  UInt    getBits( Int x, Int y 
#if VCEG_AZ07_IMV
    , Int iMvFlag
#endif
    )
  {
    return xGetExpGolombNumberOfBits((x << m_iCostScale) 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      - ( m_mvPredictor.getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE )
#else
      - m_mvPredictor.getHor()
#endif
#if VCEG_AZ07_IMV
      , iMvFlag
#endif
      )
    +      xGetExpGolombNumberOfBits((y << m_iCostScale) 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      - ( m_mvPredictor.getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE )
#else
    - m_mvPredictor.getVer()
#endif
#if VCEG_AZ07_IMV
      , iMvFlag
#endif
      );
  }

private:

  static Distortion xGetSSE           ( DistParam* pcDtParam );
  static Distortion xGetSSE4          ( DistParam* pcDtParam );
  static Distortion xGetSSE8          ( DistParam* pcDtParam );
  static Distortion xGetSSE16         ( DistParam* pcDtParam );
  static Distortion xGetSSE32         ( DistParam* pcDtParam );
  static Distortion xGetSSE64         ( DistParam* pcDtParam );
  static Distortion xGetSSE16N        ( DistParam* pcDtParam );
#if WCG_LUMA_DQP_CM_SCALE
  static Distortion xGetSSE_WTD           ( DistParam* pcDtParam );
  static Distortion xGetSSE4_WTD          ( DistParam* pcDtParam );
  static Distortion xGetSSE8_WTD          ( DistParam* pcDtParam );
  static Distortion xGetSSE16_WTD         ( DistParam* pcDtParam );
  static Distortion xGetSSE32_WTD         ( DistParam* pcDtParam );
  static Distortion xGetSSE64_WTD         ( DistParam* pcDtParam );
  static Distortion xGetSSE16N_WTD        ( DistParam* pcDtParam );
#endif
  static Distortion xGetSAD           ( DistParam* pcDtParam );
  static Distortion xGetSAD4          ( DistParam* pcDtParam );
  static Distortion xGetSAD8          ( DistParam* pcDtParam );
  static Distortion xGetSAD16         ( DistParam* pcDtParam );
  static Distortion xGetSAD32         ( DistParam* pcDtParam );
  static Distortion xGetSAD64         ( DistParam* pcDtParam );
  static Distortion xGetSAD16N        ( DistParam* pcDtParam );

  static Distortion xGetSAD12         ( DistParam* pcDtParam );
  static Distortion xGetSAD24         ( DistParam* pcDtParam );
  static Distortion xGetSAD48         ( DistParam* pcDtParam );

  static Distortion xGetHADs          ( DistParam* pcDtParam );
  static Distortion xCalcHADs2x2      ( Pel *piOrg, Pel *piCurr, Int iStrideOrg, Int iStrideCur, Int iStep );
  static Distortion xCalcHADs4x4      ( Pel *piOrg, Pel *piCurr, Int iStrideOrg, Int iStrideCur, Int iStep, Int bitDepth);
  static Distortion xCalcHADs8x8      ( Pel *piOrg, Pel *piCurr, Int iStrideOrg, Int iStrideCur, Int iStep 
#if COM16_C806_SIMD_OPT
    , Int bitDepth
#endif
    );
#if JVET_C0024_QTBT
  static Distortion xCalcHADs16x8     ( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth);
  static Distortion xCalcHADs8x16     ( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth);
  static Distortion xCalcHADs4x8     ( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth);
  static Distortion xCalcHADs8x4     ( Pel *piOrg, Pel *piCur, Int iStrideOrg, Int iStrideCur, Int bitDepth);
#endif

#if VCEG_AZ06_IC
  static Distortion xGetMRSAD         ( DistParam* pcDtParam );
  static Distortion xGetMRSAD4        ( DistParam* pcDtParam );
  static Distortion xGetMRSAD8        ( DistParam* pcDtParam );
  static Distortion xGetMRSAD16       ( DistParam* pcDtParam );
  static Distortion xGetMRSAD32       ( DistParam* pcDtParam );
  static Distortion xGetMRSAD64       ( DistParam* pcDtParam );
  static Distortion xGetMRSAD16N      ( DistParam* pcDtParam );
  static Distortion xGetMRSAD12       ( DistParam* pcDtParam );
  static Distortion xGetMRSAD24       ( DistParam* pcDtParam );
  static Distortion xGetMRSAD48       ( DistParam* pcDtParam );
  static Distortion xGetMRHADs        ( DistParam* pcDtParam );
#endif

public:

  Distortion   getDistPart(Int bitDepth, Pel* piCur, Int iCurStride,  Pel* piOrg, Int iOrgStride, UInt uiBlkWidth, UInt uiBlkHeight, const ComponentID compID, DFunc eDFunc = DF_SSE 
#if WCG_LUMA_DQP_CM_SCALE
  , Pel* piOrgLuma = NULL, 
  Int iOrgStrideLuma = 0
#endif
#if VCEG_AZ06_IC
    , Bool bMRFlag = false
#endif
    );

};// END CLASS DEFINITION TComRdCost

//! \}

#endif // __TCOMRDCOST__
