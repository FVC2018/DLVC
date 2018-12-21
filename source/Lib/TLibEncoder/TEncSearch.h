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

/** \file     TEncSearch.h
    \brief    encoder search class (header)
*/

#ifndef __TENCSEARCH__
#define __TENCSEARCH__

// Include files
#include "TLibCommon/TComYuv.h"
#include "TLibCommon/TComMotionInfo.h"
#include "TLibCommon/TComPattern.h"
#include "TLibCommon/TComPrediction.h"
#include "TLibCommon/TComTrQuant.h"
#include "TLibCommon/TComPic.h"
#include "TLibCommon/TComRectangle.h"
#include "TEncEntropy.h"
#include "TEncSbac.h"
#include "TEncCfg.h"
#if JVET_F0096_BILATERAL_FILTER
#include "TLibCommon/TComBilateralFilter.h"
#endif

//! \ingroup TLibEncoder
//! \{

class TEncCu;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

static const UInt MAX_NUM_REF_LIST_ADAPT_SR=2;
static const UInt MAX_IDX_ADAPT_SR=33;
static const UInt NUM_MV_PREDICTORS=3;

/// encoder search class
class TEncSearch : public TComPrediction
{
private:
#if JVET_C0024_QTBT
  TCoeff***       m_pppcQTTempCoeff[MAX_NUM_COMPONENT /* 0->Y, 1->Cb, 2->Cr*/];
#else
  TCoeff**        m_ppcQTTempCoeff[MAX_NUM_COMPONENT /* 0->Y, 1->Cb, 2->Cr*/];
#endif
  TCoeff*         m_pcQTTempCoeff[MAX_NUM_COMPONENT];
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
  TCoeff***       m_pppcQTTempArlCoeff[MAX_NUM_COMPONENT];
#else
  TCoeff**        m_ppcQTTempArlCoeff[MAX_NUM_COMPONENT];
#endif
  TCoeff*         m_pcQTTempArlCoeff[MAX_NUM_COMPONENT];
#endif
#if !JVET_C0024_QTBT
  UChar*          m_puhQTTempTrIdx;
#endif
  UChar*          m_puhQTTempCbf[MAX_NUM_COMPONENT];

#if JVET_C0024_QTBT
  TComYuv**       m_ppcQTTempTComYuv;
#else
  TComYuv*        m_pcQTTempTComYuv;
#endif
#if VCEG_AZ08_INTER_KLT
#if JVET_C0024_QTBT
  TComYuv**       m_ppcQTTempTComYuvRec;
#else
  TComYuv*        m_pcQTTempTComYuvRec;
#endif
#endif  
  TComYuv         m_tmpYuvPred; // To be used in xGetInterPredictionError() to avoid constant memory allocation/deallocation

  Char*           m_phQTTempCrossComponentPredictionAlpha[MAX_NUM_COMPONENT];
  Pel*            m_pSharedPredTransformSkip[MAX_NUM_COMPONENT];
  TCoeff*         m_pcQTTempTUCoeff[MAX_NUM_COMPONENT];
  UChar*          m_puhQTTempTransformSkipFlag[MAX_NUM_COMPONENT];
#if VCEG_AZ08_KLT_COMMON
  UChar*          m_puhQTTempKLTFlag[MAX_NUM_COMPONENT];
#endif
#if COM16_C806_EMT
  UChar*          m_puhQTTempExplicitRdpcmMode[MAX_NUM_COMPONENT];
#endif
  TComYuv         m_pcQTTempTransformSkipTComYuv;

#if COM16_C806_LMCHROMA
  TComYuv         m_pcQTTempResiTComYuv;
#endif

#if ADAPTIVE_QP_SELECTION
  TCoeff*         m_ppcQTTempTUArlCoeff[MAX_NUM_COMPONENT];
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComMvField   * m_pMvFieldSP[NUM_MGR_TYPE];
  UChar         * m_phInterDirSP[NUM_MGR_TYPE];
#else
  TComMvField   * m_pMvFieldSP[2];
  UChar         * m_phInterDirSP[2];
#endif
#endif
#if VCEG_AZ07_FRUC_MERGE
  TComMvField   * m_pMvFieldFRUC;
  UChar         * m_phInterDirFRUC;
  UChar         * m_phFRUCRefineDist[2];
  UChar         * m_phFRUCSBlkRefineDist[2];
#endif


#if   JVET_E0077_LM_MF
  Pel             **m_pLMMFPredSaved;
#endif

#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
  Pel*            m_resiSingleBuffer[MAX_CU_DEPTH][MAX_CU_DEPTH][NUMBER_OF_STORED_RESIDUAL_TYPES];
#else
  Pel*            m_resiSplitBuffer[MAX_CU_DEPTH][NUMBER_OF_STORED_RESIDUAL_TYPES];
  Pel*            m_resiSingleBuffer[MAX_CU_DEPTH][NUMBER_OF_STORED_RESIDUAL_TYPES];
#endif
  Pel*            m_resiPUBuffer[NUMBER_OF_STORED_RESIDUAL_TYPES];
#endif

#if COM16_C1016_AFFINE
#if !ITH_C2_FXD_PT_AFF_PARAMS
  Int*            m_tmpError;
  Double*         m_tmpDerivate[2];
#else
  Int*            m_tmpDerivate[2];
  Pel*            m_tmpError;
#endif
#endif

#if DIMD_INTRA_PRED
  UChar           m_ucDIMDIntraMode;
  UChar           m_ucParentDIMDIntraMode;
#endif


#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  UInt*           m_SaveLoadNumRouteSaved[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UInt*           m_SaveLoadRouteCodingOrder[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UInt*           m_SaveLoadRouteLeftNeigbor[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UInt*           m_SaveLoadRouteUpNeigbor[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];

  UChar*          m_SaveLoadTag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];     ///< 0: not used, 1: save, 2: load
#if COM16_C806_EMT
  UChar*          m_SaveLoadEmtFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
  UChar*          m_SaveLoadEmtIdx[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
  Char*          m_SaveLoadRotIdx[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if COM16_C1046_PDPC_INTRA
  Char*          m_SaveLoadPdpcIdx[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if VCEG_AZ07_FRUC_MERGE
  UChar*          m_SaveLoadFrucMode[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
  UChar*          m_SaveLoadIMVFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#else
  Bool*           m_SaveLoadIMVFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#endif
#if VCEG_AZ06_IC
  Bool*           m_SaveLoadICFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#if COM16_C1016_AFFINE
  Bool*           m_SaveLoadAffineFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
  Bool*           m_SaveLoadMergeFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UChar*          m_SaveLoadInterDir[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UChar*          m_SaveLoadSplit[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#else
  UInt            m_SaveLoadPartIdx[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; ///< partition index of the block for save/load encoder decision 
  UChar           m_SaveLoadTag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];     ///< 0: not used, 1: save, 2: load
#if COM16_C806_EMT
  UChar           m_SaveLoadEmtFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
  UChar           m_SaveLoadEmtIdx [MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
   Char           m_SaveLoadRotIdx [MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if RRC_INTRA_SHARPENING && INTRA_SHARP_D0077
   Char           m_SaveLoadSharpIdx [MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if COM16_C1046_PDPC_INTRA
   Char           m_SaveLoadPdpcIdx[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; 
#endif
#if VCEG_AZ07_FRUC_MERGE
  UChar           m_SaveLoadFrucMode[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
  UChar           m_SaveLoadIMVFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#else
  Bool            m_SaveLoadIMVFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#endif
#if VCEG_AZ06_IC
  Bool            m_SaveLoadICFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#if COM16_C1016_AFFINE
  Bool            m_SaveLoadAffineFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
  Bool            m_SaveLoadMergeFlag[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UChar           m_SaveLoadInterDir[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
  UChar           m_SaveLoadSplit[MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1];
#endif
#endif


#if JVET_D0123_ME_CTX_LUT_BITS
  Int             iCostScale;
#endif

#if ITH_C2_FXD_PT_AFF_PARAMS
#if JVECT_C0062_AFFINE_SIX_PARAM
  static const Int iEqualCoeffBufDimension = 7;
#else
  static const Int iEqualCoeffBufDimension = 5;
#endif
#endif

protected:
  // interface to option
  TEncCfg*        m_pcEncCfg;

  // interface to classes
  TComTrQuant*    m_pcTrQuant;
  TComRdCost*     m_pcRdCost;
  TEncEntropy*    m_pcEntropyCoder;

  // ME parameters
  Int             m_iSearchRange;
  Int             m_bipredSearchRange; // Search range for bi-prediction
  Int             m_iFastSearch;
  Int             m_aaiAdaptSR[MAX_NUM_REF_LIST_ADAPT_SR][MAX_IDX_ADAPT_SR];
  TComMv          m_cSrchRngLT;
  TComMv          m_cSrchRngRB;
  TComMv          m_acMvPredictors[NUM_MV_PREDICTORS]; // Left, Above, AboveRight. enum MVP_DIR first NUM_MV_PREDICTORS entries are suitable for accessing.

  // RD computation
#if JVET_C0024_QTBT
  TEncSbac****    m_ppppcRDSbacCoder;
#else
  TEncSbac***     m_pppcRDSbacCoder;
#endif
  TEncSbac*       m_pcRDGoOnSbacCoder;
  DistParam       m_cDistParam;

#if VCEG_AZ07_IMV && JVET_D0123_ME_CTX_LUT_BITS
  UInt        m_uiBitsIMVFlag[2];
#endif
  // Misc.
  Pel*            m_pTempPel;
  const UInt*     m_puiDFilter;

  // AMVP cost computation
  // UInt            m_auiMVPIdxCost[AMVP_MAX_NUM_CANDS+1][AMVP_MAX_NUM_CANDS];
  UInt            m_auiMVPIdxCost[AMVP_MAX_NUM_CANDS+1][AMVP_MAX_NUM_CANDS+1]; //th array bounds

  TComMv          m_integerMv2Nx2N[NUM_REF_PIC_LIST_01][MAX_NUM_REF];

  Bool            m_isInitialized;

#if COM16_C806_EMT
  UChar*          m_puhQTTempEmtTuIdx;
  UChar*          m_puhQTTempEmtCuFlag;
#endif

public:
#if JVET_D0123_ME_CTX_LUT_BITS
  estPuMeBitsSbacStruct* m_pcPuMeEstBitsSbac;
#endif
  TEncSearch();
  virtual ~TEncSearch();

  Void init(TEncCfg*      pcEncCfg,
            TComTrQuant*  pcTrQuant,
            Int           iSearchRange,
            Int           bipredSearchRange,
            Int           iFastSearch,
            const UInt    maxCUWidth,
            const UInt    maxCUHeight,
            const UInt    maxTotalCUDepth,
            TEncEntropy*  pcEntropyCoder,
            TComRdCost*   pcRdCost,
#if JVET_C0024_QTBT
            TEncSbac****  ppppcRDSbacCoder,
#else
            TEncSbac***   pppcRDSbacCoder,
#endif
            TEncSbac*     pcRDGoOnSbacCoder );

  Void destroy();

#if DIMD_INTRA_PRED
  UChar getDIMDIntraMode() { return m_ucDIMDIntraMode; }
  Void  setDIMDIntraMode(UChar uc) { m_ucDIMDIntraMode = uc; }
  UChar getParentDIMDIntraMode() { return m_ucParentDIMDIntraMode; }
  Void  setParentDIMDIntraMode(UChar uc) { m_ucParentDIMDIntraMode = uc; }
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT // search RouteInfo function
  Void  cleanSaveLoadInfo(UInt uiPartIdx, Int uiWIdx, Int uiHIdx);
  UChar searchSaveLoadIdx( UInt uiPartIdx, UInt uiWIdx, UInt uiHIdx, UInt uiBTDepth, UInt* puiCURouteInfo);
  
  UChar getSaveLoadTag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadTag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }
  Void  setSaveLoadTag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadTag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; }
#if COM16_C806_EMT
  UChar getSaveLoadEmtFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadEmtFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }
  Void  setSaveLoadEmtFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadEmtFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; }
  UChar getSaveLoadEmtIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadEmtIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadEmtIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadEmtIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
  Char  getSaveLoadRotIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadRotIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadRotIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Char c ) { m_SaveLoadRotIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
#endif
#if COM16_C1046_PDPC_INTRA
  Char  getSaveLoadPdpcIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadPdpcIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadPdpcIdx( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Char c ) { m_SaveLoadPdpcIdx[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
#endif
#if VCEG_AZ07_FRUC_MERGE
  UChar getSaveLoadFrucMode( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadFrucMode[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadFrucMode( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadFrucMode[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
  UChar getSaveLoadIMVFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadIMVFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadIMVFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar b ) { m_SaveLoadIMVFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = b; };
#else
  Bool  getSaveLoadIMVFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadIMVFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadIMVFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadIMVFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = b; };
#endif
#endif
#if VCEG_AZ06_IC
  Bool  getSaveLoadICFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadICFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadICFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadICFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = b; };
#endif
#if COM16_C1016_AFFINE
  Bool  getSaveLoadAffineFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadAffineFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadAffineFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadAffineFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = b; };
#endif
  Bool  getSaveLoadMergeFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadMergeFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadMergeFlag( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadMergeFlag[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = b; };
  UChar getSaveLoadInterDir( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadInterDir[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadInterDir( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadInterDir[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
  UChar getSaveLoadSplit( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadSplit[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx]; }; 
  Void  setSaveLoadSplit( UInt uiPartIdx, UInt uiRouteIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadSplit[uiWIdx][uiHIdx][uiPartIdx+MAX_PART_IN_CTU*uiRouteIdx] = c; };
#else
  UChar getSaveLoadTag( UInt uiPartIdx, UInt uiWIdx, UInt uiHIdx ) {  return uiPartIdx == m_SaveLoadPartIdx[uiWIdx][uiHIdx] ? m_SaveLoadTag[uiWIdx][uiHIdx] : SAVE_LOAD_INIT; };
  Void  setSaveLoadTag( UInt uiPartIdx, UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadPartIdx[uiWIdx][uiHIdx] = uiPartIdx; m_SaveLoadTag[uiWIdx][uiHIdx] = c; };
#if COM16_C806_EMT
  UChar getSaveLoadEmtFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadEmtFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadEmtFlag( UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadEmtFlag[uiWIdx][uiHIdx] = c; };
  UChar getSaveLoadEmtIdx( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadEmtIdx[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadEmtIdx( UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadEmtIdx[uiWIdx][uiHIdx] = c; };
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
  Char  getSaveLoadRotIdx( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadRotIdx[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadRotIdx( UInt uiWIdx, UInt uiHIdx, Char c ) { m_SaveLoadRotIdx[uiWIdx][uiHIdx] = c; };
#endif
#if RRC_INTRA_SHARPENING && INTRA_SHARP_D0077
  Char  getSaveLoadSharpIdx(UInt uiWIdx, UInt uiHIdx) { return m_SaveLoadSharpIdx[uiWIdx][uiHIdx]; };
  Void  setSaveLoadSharpIdx(UInt uiWIdx, UInt uiHIdx, Char c) { m_SaveLoadSharpIdx[uiWIdx][uiHIdx] = c; };
#endif
#if COM16_C1046_PDPC_INTRA
  Char  getSaveLoadPdpcIdx( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadPdpcIdx[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadPdpcIdx( UInt uiWIdx, UInt uiHIdx, Char c ) { m_SaveLoadPdpcIdx[uiWIdx][uiHIdx] = c; };
#endif
#if VCEG_AZ07_FRUC_MERGE
  UChar getSaveLoadFrucMode( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadFrucMode[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadFrucMode( UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadFrucMode[uiWIdx][uiHIdx] = c; };
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
  UChar  getSaveLoadIMVFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadIMVFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadIMVFlag( UInt uiWIdx, UInt uiHIdx, UChar b ) { m_SaveLoadIMVFlag[uiWIdx][uiHIdx] = b; };
#else
  Bool  getSaveLoadIMVFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadIMVFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadIMVFlag( UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadIMVFlag[uiWIdx][uiHIdx] = b; };
#endif
#endif
#if VCEG_AZ06_IC
  Bool  getSaveLoadICFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadICFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadICFlag( UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadICFlag[uiWIdx][uiHIdx] = b; };
#endif
#if COM16_C1016_AFFINE
  Bool  getSaveLoadAffineFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadAffineFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadAffineFlag( UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadAffineFlag[uiWIdx][uiHIdx] = b; };
#endif
  Bool  getSaveLoadMergeFlag( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadMergeFlag[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadMergeFlag( UInt uiWIdx, UInt uiHIdx, Bool b ) { m_SaveLoadMergeFlag[uiWIdx][uiHIdx] = b; };
  UChar getSaveLoadInterDir( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadInterDir[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadInterDir( UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadInterDir[uiWIdx][uiHIdx] = c; };
  UChar getSaveLoadSplit( UInt uiWIdx, UInt uiHIdx ) {  return m_SaveLoadSplit[uiWIdx][uiHIdx]; }; 
  Void  setSaveLoadSplit( UInt uiWIdx, UInt uiHIdx, UChar c ) { m_SaveLoadSplit[uiWIdx][uiHIdx] = c; };
#endif
#endif

#if MRG_CONSIDER_COMPLEX
  Void xMergeComplexEstimation(TComDataCU*     pcCU,
    TComYuv*        pcYuvOrg,
    Int             iPartIdx,
    UInt&           uiComplexInterDir,
    NodeInfo&       NodeMotionInfoP,
    UInt&           uiMergeComplexIndex,
    UInt&           ruiCost
  );
  Int getComplexIndexCodeBits(UInt numValidMergeComplexCand, UInt uiMergeComplexIndex);
#endif

protected:

  /// sub-function for motion vector refinement used in fractional-pel accuracy
  Distortion  xPatternRefinement( TComPattern* pcPatternKey,
                                  TComMv baseRefMv,
                                  Int iFrac, TComMv& rcMvFrac, Bool bAllowUseOfHadamard
#if RRC_SELECTABLE_IF
                                  , TComDataCU* pcCU
                                  , TComPic *pcRefPic
                                  , UInt uiPartAddr
#endif
                                 );

  typedef struct
  {
    Pel*        piRefY;
    Int         iYStride;
    Int         iBestX;
    Int         iBestY;
    UInt        uiBestRound;
    UInt        uiBestDistance;
    Distortion  uiBestSad;
#if JVET_D0123_ME_CTX_LUT_BITS
    Distortion  uiBestSadOnly;
#endif
    UChar       ucPointNr;
  } IntTZSearchStruct;

  // sub-functions for ME
#if JVET_D0123_ME_CTX_LUT_BITS
  __inline UInt xGetBitsEpExGolomb( UInt uiSymbol, UInt uiCount );
  __inline UInt xGetMvdBitsLut (TComDataCU* pcCU, UInt uiAbsPartIdx, Int iCurVal);
  __inline UInt xMvdBits(TComPattern* pcPatternKey, Int iX, Int iY);
#endif
  __inline Void xTZSearchHelp         ( TComPattern* pcPatternKey, IntTZSearchStruct& rcStruct, const Int iSearchX, const Int iSearchY, const UChar ucPointNr, const UInt uiDistance );
  __inline Void xTZ2PointSearch       ( TComPattern* pcPatternKey, IntTZSearchStruct& rcStrukt, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB );
  __inline Void xTZ8PointSquareSearch ( TComPattern* pcPatternKey, IntTZSearchStruct& rcStrukt, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist );
  __inline Void xTZ8PointDiamondSearch( TComPattern* pcPatternKey, IntTZSearchStruct& rcStrukt, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist );

  Void xGetInterPredictionError( TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, Distortion& ruiSAD, Bool Hadamard );

public:
  Void  estIntraPredLumaQT      ( TComDataCU* pcCU,
                                  TComYuv*    pcOrgYuv,
                                  TComYuv*    pcPredYuv,
                                  TComYuv*    pcResiYuv,
                                  TComYuv*    pcRecoYuv,
// #if DIMD_MLIP_HARM
//       Bool bUseDIMD,
// #endif
#if COM16_C806_LARGE_CTU
                                  Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
                                  Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
                                  DEBUG_STRING_FN_DECLARE(sDebug));

#if DIMD_INTRA_PRED
  Void  estDIMDIntraPredLumaQT(TComDataCU* pcCU,
	  TComYuv*    pcOrgYuv,
	  TComYuv*    pcPredYuv,
	  TComYuv*    pcResiYuv,
	  TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
	  Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
	  Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
	  );

  Void  estDIMDIntraPredChromaQT(TComDataCU* pcCU,
	  TComYuv*    pcOrgYuv,
	  TComYuv*    pcPredYuv,
	  TComYuv*    pcResiYuv,
	  TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
	  Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
	  Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
	  );
#endif

  Void  estIntraPredChromaQT    ( TComDataCU* pcCU,
                                  TComYuv*    pcOrgYuv,
                                  TComYuv*    pcPredYuv,
                                  TComYuv*    pcResiYuv,
                                  TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
                                  Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
                                  Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
                                  DEBUG_STRING_FN_DECLARE(sDebug));

  /// encoder estimation - inter prediction (non-skip)
  Void predInterSearch          ( TComDataCU* pcCU,
#if AFF_EARLY_TERMINATION
		                              TComDataCU* pcBestCU,
#endif
                                  TComYuv*    pcOrgYuv,
                                  TComYuv*    pcPredYuv,
                                  TComYuv*    pcResiYuv,
                                  TComYuv*    pcRecoYuv
                                  DEBUG_STRING_FN_DECLARE(sDebug),
#if COM16_C806_OBMC
                                  TComYuv*    pcPredYuvWoOBMC, 
                                  TComYuv*    pcTmpYuv1,
                                  TComYuv*    pcTmpYuv2,
#endif
                                  Bool        bUseRes = false
#if AMP_MRG
                                 ,Bool        bUseMRG = false
#endif
                                );

  /// encode residual and compute rd-cost for inter mode
  Void encodeResAndCalcRdInterCU( TComDataCU* pcCU,
                                  TComYuv*    pcYuvOrg,
                                  TComYuv*    pcYuvPred,
                                  TComYuv*    pcYuvResi,
                                  TComYuv*    pcYuvResiBest,
                                  TComYuv*    pcYuvRec,
                                  Bool        bSkipResidual
#if COM16_C806_EMT
                                  , Double dBestCost
#endif
                                  DEBUG_STRING_FN_DECLARE(sDebug) );

  /// set ME search range
  Void setAdaptiveSearchRange   ( Int iDir, Int iRefIdx, Int iSearchRange) { assert(iDir < MAX_NUM_REF_LIST_ADAPT_SR && iRefIdx<Int(MAX_IDX_ADAPT_SR)); m_aaiAdaptSR[iDir][iRefIdx] = iSearchRange; }

  Void xEncPCM    (TComDataCU* pcCU, UInt uiAbsPartIdx, Pel* piOrg, Pel* piPCM, Pel* piPred, Pel* piResi, Pel* piReco, UInt uiStride, UInt uiWidth, UInt uiHeight, const ComponentID compID );
  Void IPCMSearch (TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* rpcPredYuv, TComYuv* rpcResiYuv, TComYuv* rpcRecoYuv );
#if JVET_C0024_FAST_MRG
  static UInt updateCandList( UInt uiMode, Double uiCost, UInt uiFastCandNum, UInt * CandModeList, Double * CandCostList
#if MLIP_EXTEND_REF_LINE
      , UInt * ExternRefList = NULL, UInt ExternRef = 0
#endif
  );
#if JVET_D0123_ME_CTX_LUT_BITS && JVET_C0024_FAST_MRG 
  Void getMrgCandBits( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiMrgIdxBits[MRG_MAX_NUM_CANDS]);
#endif
#endif
protected:

  // -------------------------------------------------------------------------------------------------------------------
  // Intra search
  // -------------------------------------------------------------------------------------------------------------------

  Void  xEncSubdivCbfQT           ( TComTU      &rTu,
                                    Bool         bLuma,
                                    Bool         bChroma );

  Void  xEncCoeffQT               ( TComTU &rTu,
                                    ComponentID  component,
                                    Bool         bRealCoeff );
  Void  xEncIntraHeader           ( TComDataCU*  pcCU,
                                    UInt         uiTrDepth,
                                    UInt         uiAbsPartIdx,
                                    Bool         bLuma,
                                    Bool         bChroma );
  UInt  xGetIntraBitsQT           ( TComTU &rTu,
                                    Bool         bLuma,
                                    Bool         bChroma,
                                    Bool         bRealCoeff );

  UInt  xGetIntraBitsQTChroma    ( TComTU &rTu,
                                   ComponentID compID,
                                   Bool          bRealCoeff );

  Void  xIntraCodingTUBlock       (       TComYuv*      pcOrgYuv,
                                          TComYuv*      pcPredYuv,
                                          TComYuv*      pcResiYuv,
#if COM16_C806_LARGE_CTU
                                          Pel*          resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
                                          Pel           resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
                                    const Bool          checkCrossCPrediction,
                                          Distortion&   ruiDist,
                                    const ComponentID   compID,
                                          TComTU        &rTu
                                    DEBUG_STRING_FN_DECLARE(sTest)
                                         ,Int           default0Save1Load2 = 0
#if COM16_C806_EMT
                                         ,UInt*         puiSigNum = NULL
#endif
                                   );
#if VCEG_AZ08_INTRA_KLT
  Bool xIntraCodingTUBlockTM(TComYuv*    pcOrgYuv,
                                      TComYuv*    pcPredYuv,
                                      TComYuv*    pcResiYuv,
                                      Distortion& ruiDist,
                                      const ComponentID compID,
                                      TComTU&     rTu
                                      DEBUG_STRING_FN_DECLARE(sDebug)
#if COM16_C806_EMT
                                      , UInt*      puiSigNum
#endif
                                      , Int tmpred0_tmpredklt1_ori2
                                      );
#endif

  Void  xRecurIntraCodingLumaQT   ( TComYuv*    pcOrgYuv,
                                    TComYuv*    pcPredYuv,
                                    TComYuv*    pcResiYuv,
#if COM16_C806_LARGE_CTU
                                    Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
                                    Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
                                    Distortion& ruiDistY,
#if HHI_RQT_INTRA_SPEEDUP
                                    Bool         bCheckFirst,
#endif
                                    Double&      dRDCost,
                                    TComTU      &rTu
                                    DEBUG_STRING_FN_DECLARE(sDebug));

  Void  xSetIntraResultLumaQT     ( TComYuv*     pcRecoYuv,
                                    TComTU &rTu);

  Void xStoreCrossComponentPredictionResult  (       Pel    *pResiLuma,
                                               const Pel    *pBestLuma,
                                                     TComTU &rTu,
                                               const Int     xOffset,
                                               const Int     yOffset,
                                               const Int     strideResi,
                                               const Int     strideBest );

  Char xCalcCrossComponentPredictionAlpha    (       TComTU &rTu,
                                               const ComponentID compID,
                                               const Pel*        piResiL,
                                               const Pel*        piResiC,
                                               const Int         width,
                                               const Int         height,
                                               const Int         strideL,
                                               const Int         strideC );

  Void  xRecurIntraChromaCodingQT ( TComYuv*    pcOrgYuv,
                                    TComYuv*    pcPredYuv,
                                    TComYuv*    pcResiYuv,
#if COM16_C806_LARGE_CTU
                                    Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
                                    Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
                                    Distortion& ruiDist,
                                    TComTU      &rTu
                                    DEBUG_STRING_FN_DECLARE(sDebug));

  Void  xSetIntraResultChromaQT   ( TComYuv*    pcRecoYuv, TComTU &rTu);

#if COM16_C806_LMCHROMA
  Void  xStoreIntraResultQT       ( const ComponentID compID, TComTU &rTu, TComYuv* pcResiYuv = NULL);
  Void  xLoadIntraResultQT        ( const ComponentID compID, TComTU &rTu, TComYuv* pcResiYuv = NULL);
#else
  Void  xStoreIntraResultQT       ( const ComponentID compID, TComTU &rTu);
  Void  xLoadIntraResultQT        ( const ComponentID compID, TComTU &rTu);
#endif


  // -------------------------------------------------------------------------------------------------------------------
  // Inter search (AMP)
  // -------------------------------------------------------------------------------------------------------------------

  Void xEstimateMvPredAMVP        ( TComDataCU* pcCU,
                                    TComYuv*    pcOrgYuv,
                                    UInt        uiPartIdx,
                                    RefPicList  eRefPicList,
                                    Int         iRefIdx,
                                    TComMv&     rcMvPred,
                                    Bool        bFilled = false
                                  , Distortion* puiDistBiP = NULL
                                     );

  Void xCheckBestMVP              ( TComDataCU* pcCU,
                                    RefPicList  eRefPicList,
                                    TComMv      cMv,
                                    TComMv&     rcMvPred,
                                    Int&        riMVPIdx,
                                    UInt&       ruiBits,
                                    Distortion& ruiCost 
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
                                    , UInt uiPartAddr
#endif
                                    );

#if JVET_E0076_MULTI_PEL_MVD
  Bool intMvRefineNeeded (TComDataCU* pcCU, Int iPartIdx) { return pcCU->getiMVFlag(iPartIdx) != 0 ? true : false; }
  Void intMvRefine ( TComDataCU* pcCU, RefPicList eRefPicList, TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, TComMv& rcMv, TComMv& rcMvPred, Int& riMVPIdx, UInt& ruiBits, Distortion& ruiCost, Double fWeight  );
#endif

  Distortion xGetTemplateCost    ( TComDataCU*  pcCU,
                                    UInt        uiPartAddr,
                                    TComYuv*    pcOrgYuv,
                                    TComYuv*    pcTemplateCand,
                                    TComMv      cMvCand,
                                    Int         iMVPIdx,
                                    Int         iMVPNum,
                                    RefPicList  eRefPicList,
                                    Int         iRefIdx,
                                    Int         iSizeX,
                                    Int         iSizeY
                                   );


  Void xCopyAMVPInfo              ( AMVPInfo*   pSrc, AMVPInfo* pDst );
  UInt xGetMvpIdxBits             ( Int iIdx, Int iNum );
  Void xGetBlkBits                ( PartSize  eCUMode, Bool bPSlice, Int iPartIdx,  UInt uiLastMode, UInt uiBlkBit[3]
#if JVET_D0123_ME_CTX_LUT_BITS
  ,TComDataCU* pcCU, UInt uiAbsPartIdx
#endif
    );

#if JVET_D0123_ME_CTX_LUT_BITS
#if VCEG_AZ07_IMV
  Void xGetIMVFlagBits        (TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
  Void loadCtxMe(TComDataCU* pcCU);
  Void loadCtxMeSbacCoder(TComDataCU* pcCU);
  Void storeCtxMeSbacCoder(TComDataCU* pcCU);
  UInt getInterDirPUBits(UInt uiInterDir, UInt uiCtx, TComDataCU* pcCU, UInt uiAbsPartIdx);
  Void xGetMrgFlagBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt mrgFlagBits[2]);
  UInt xGetMrgIdxBits( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiMergeCand);
  Void xUpdateMvpIdxBits(TComDataCU* pcCU);
  UInt xRefFrameIdxBits(TComDataCU* pcCU, Int iRefFrame, UInt uiRefNum);
#if COM16_C1016_AFFINE
  Void xGetAffineFlagBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt iPUIdx, UInt uiAffineFlagBits[2]);
#endif
#if VCEG_AZ07_FRUC_MERGE
  UInt xGetFRUCMgrModeBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt iPUIdx, UChar uhFRUCMode);
#endif
#endif

  Void xMergeEstimation           ( TComDataCU*  pcCU,
                                    TComYuv*     pcYuvOrg,
                                    Int          iPartIdx,
                                    UInt&        uiInterDir,
                                    TComMvField* pacMvField,
                                    UInt&        uiMergeIndex,
                                    Distortion&  ruiCost,
                                    TComMvField* cMvFieldNeighbours,
                                    UChar*       uhInterDirNeighbours,
                                    Int&         numValidMergeCand
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
                                  , UChar*      pMergeTypeNeighbor 
#if JVET_C0035_ATMVP_SIMPLIFICATION
                                  , TComMvField*    pcMvFieldSP[NUM_MGR_TYPE]
                                  , UChar*          puhInterDirSP[NUM_MGR_TYPE]
#else
                                  , TComMvField*    pcMvFieldSP[2]
                                  , UChar*          puhInterDirSP[2]
#endif
#endif
                                   );

#if VCEG_AZ07_FRUC_MERGE
  Void xFRUCMgrEstimation         ( TComDataCU* pcCU, 
                                    TComYuv* pcYuvOrg, 
                                    Int iPUIdx, 
                                    TComMvField* pacMvField, 
                                    UChar * phInterDir , 
                                    UChar ** phFRUCRefineDist , 
                                    UChar ** phFRUCSBlkRefineDist ,
                                    UInt& ruiMinCost , 
                                    UChar & ruhFRUCMode );
#endif

  Void xRestrictBipredMergeCand   ( TComDataCU*     pcCU,
                                    UInt            puIdx,
                                    TComMvField*    mvFieldNeighbours,
                                    UChar*          interDirNeighbours,
                                    Int             numValidMergeCand );


  // -------------------------------------------------------------------------------------------------------------------
  // motion estimation
  // -------------------------------------------------------------------------------------------------------------------

  Void xMotionEstimation          ( TComDataCU*  pcCU,
                                    TComYuv*     pcYuvOrg,
                                    Int          iPartIdx,
                                    RefPicList   eRefPicList,
                                    TComMv*      pcMvPred,
                                    Int          iRefIdxPred,
                                    TComMv&      rcMv,
#if JVET_E0076_MULTI_PEL_MVD
                                    Int& riMVPIdx, 
#endif
                                    UInt&        ruiBits,
                                    Distortion&  ruiCost,
                                    Bool         bBi = false  );

  Void xTZSearch                  ( TComDataCU*  pcCU,
                                    TComPattern* pcPatternKey,
                                    Pel*         piRefY,
                                    Int          iRefStride,
                                    TComMv*      pcMvSrchRngLT,
                                    TComMv*      pcMvSrchRngRB,
                                    TComMv&      rcMv,
                                    Distortion&  ruiSAD,
                                    const TComMv *pIntegerMv2Nx2NPred
                                    );

  Void xTZSearchSelective         ( TComDataCU*  pcCU,
                                    TComPattern* pcPatternKey,
                                    Pel*         piRefY,
                                    Int          iRefStride,
                                    TComMv*      pcMvSrchRngLT,
                                    TComMv*      pcMvSrchRngRB,
                                    TComMv&      rcMv,
                                    Distortion&  ruiSAD,
                                    const TComMv *pIntegerMv2Nx2NPred
                                    );

  Void xSetSearchRange            ( TComDataCU*  pcCU,
                                    TComMv&      cMvPred,
                                    Int          iSrchRng,
                                    TComMv&      rcMvSrchRngLT,
                                    TComMv&      rcMvSrchRngRB );

  Void xPatternSearchFast         ( TComDataCU*  pcCU,
                                    TComPattern* pcPatternKey,
                                    Pel*         piRefY,
                                    Int          iRefStride,
                                    TComMv*      pcMvSrchRngLT,
                                    TComMv*      pcMvSrchRngRB,
                                    TComMv&      rcMv,
                                    Distortion&  ruiSAD,
                                    const TComMv* pIntegerMv2Nx2NPred
                                  );

  Void xPatternSearch             ( TComPattern* pcPatternKey,
                                    Pel*         piRefY,
                                    Int          iRefStride,
                                    TComMv*      pcMvSrchRngLT,
                                    TComMv*      pcMvSrchRngRB,
                                    TComMv&      rcMv,
                                    Distortion&  ruiSAD );

  Void xPatternSearchFracDIF      (
                                    Bool         bIsLosslessCoded,
                                    TComPattern* pcPatternKey,
                                    Pel*         piRefY,
                                    Int          iRefStride,
                                    TComMv*      pcMvInt,
                                    TComMv&      rcMvHalf,
                                    TComMv&      rcMvQter,
                                    Distortion&  ruiCost
#if RRC_SELECTABLE_IF
                                    , TComDataCU* pcCU
                                    , TComPic *pcRefPic
                                    , UInt uiPartAddr
#endif
                                   );

  #if COM16_C1016_AFFINE
  Void predAffineInterSearch    ( TComDataCU* pcCU,
                                  TComYuv*    pcOrgYuv,
                                  Int         iPartIdx,
                                  UInt&       ruiLastMode,
                                  Distortion& ruiAffineCost,
                                  TComMv      cHevcMv[2][33]
#if JVECT_C0062_AFFINE_SIX_PARAM
                                  , UInt affineParam
                                  , UInt&uiBestAffineParam
                                  , TComMv cMVBest[2][3], Int iBestRefIdx[2], TComMv bestMVD[2][3]
                                  , Int&bestDir, Int bestMVPIdx[2], Int bestMVPNum[2]
#endif
#if AFFINE_INHERIT_FOUR_PARAM
		, TComMv cMvAffineParam[2][33][3]
#endif
                                );

  Void xAffineMotionEstimation  ( TComDataCU*   pcCU,
                                  TComYuv*      pcYuvOrg,
                                  Int           iPartIdx,
                                  RefPicList    eRefPicList,
                                  TComMv        acMvPred[3],
                                  Int           iRefIdxPred,
                                  TComMv        acMv[3],
                                  UInt&         ruiBits,
                                  Distortion&   ruiCost,
#if JVECT_C0062_AFFINE_SIX_PARAM
                                  Bool          b6Param,
#endif
                                  Bool          bBi = false  );
  
  Bool xEstimateAffineAMVP      ( TComDataCU*   pcCU, 
                                  TComYuv*      pcOrgYuv, 
                                  UInt          uiPartIdx, 
                                  RefPicList    eRefPicList, 
                                  Int           iRefIdx,
                                  TComMv        acMvPred[3], 
                                  Distortion*   puiDistBiP
#if JVECT_C0062_AFFINE_SIX_PARAM
                                  , Bool b6Param
#endif
                                );

  Distortion xGetAffineTemplateCost ( TComDataCU* pcCU,
                                      UInt        uiPartAddr,
                                      TComYuv*    pcOrgYuv,
                                      TComYuv*    pcTemplateCand,
                                      TComMv      acMvCand[3],
                                      Int         iMVPIdx,
                                      Int         iMVPNum,
                                      RefPicList  eRefPicList,
                                      Int         iRefIdx,
                                      Int         iSizeX,
                                      Int         iSizeY
#if JVECT_C0062_AFFINE_SIX_PARAM
                                      , Bool b6Param
#endif
                                    );

  Void xCopyAffineAMVPInfo     ( AffineAMVPInfo* pSrc, AffineAMVPInfo* pDst );
#if ITH_C2_FXD_PT_AFF_PARAMS
  Void solveEqual(Double(*dEqualCoeff)[TEncSearch::iEqualCoeffBufDimension], Int iOrder, Double* dAffinePara);
#else
  Void solveEqual( Double** dEqualCoeff, Int iOrder, Double* dAffinePara );
#endif
  Void xCheckBestAffineMVP      ( TComDataCU* pcCU,
                                  RefPicList  eRefPicList,
                                  TComMv acMv[3], TComMv acMvPred[3],
                                  Int&        riMVPIdx,
                                  UInt&       ruiBits,
                                  Distortion& ruiCost
#if JVECT_C0062_AFFINE_SIX_PARAM
                                  , Bool b6Param
#endif
#if JVET_D0123_ME_CTX_LUT_BITS
                                  , UInt uiPartAddr
#endif
                                  );
#endif

  Void xExtDIFUpSamplingH( TComPattern* pcPattern );
  Void xExtDIFUpSamplingQ( TComPattern* pcPatternKey, TComMv halfPelRef 
#if SIF_GAUSS
    , TComDataCU* pcCU
    , TComPic *pcRefPic
    , UInt uiPartAddr
#endif
    );

  // -------------------------------------------------------------------------------------------------------------------
  // T & Q & Q-1 & T-1
  // -------------------------------------------------------------------------------------------------------------------

#if !JVET_C0024_QTBT
  Void xEncodeInterResidualQT( const ComponentID compID, TComTU &rTu );
#endif  
  Void xEstimateInterResidualQT( TComYuv* pcResi, Double &rdCost, UInt &ruiBits, Distortion &ruiDist, Distortion *puiZeroDist, TComTU &rTu 
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER
      , TComYuv* pcPred
#endif
      DEBUG_STRING_FN_DECLARE(sDebug) );
  Void xSetInterResidualQTData( TComYuv* pcResi, Bool bSpatial, TComTU &rTu  );

  UInt  xModeBitsIntra ( TComDataCU* pcCU, UInt uiMode, UInt uiPartOffset, UInt uiDepth, const ChannelType compID 
#if VCEG_AZ07_INTRA_65ANG_MODES
    , Int* piModes = NULL
    , Int  iAboveLeftCase = -1
#endif
    );
#if !JVET_C0024_FAST_MRG
  UInt  xUpdateCandList( UInt uiMode, Double uiCost, UInt uiFastCandNum, UInt * CandModeList, Double * CandCostList );
#endif

  // -------------------------------------------------------------------------------------------------------------------
  // compute symbol bits
  // -------------------------------------------------------------------------------------------------------------------

  Void xAddSymbolBitsInter       ( TComDataCU*   pcCU,
                                   UInt&         ruiBits);

  Void  setWpScalingDistParam( TComDataCU* pcCU, Int iRefIdx, RefPicList eRefPicListCur );
  inline  Void  setDistParamComp( ComponentID compIdx )  { m_cDistParam.compIdx = compIdx; }

	private:
#if ITH_C2_FXD_PT_AFF_PARAMS
		Void equalCoeffComputer
			(
				Pel *pResidue,
				Int iResidueStride,
				Int **ppiDerivate,
				Int iDerivateBufStride,
				Int64(*pi64EqualCoeff)[iEqualCoeffBufDimension],
				Int iWidth,
				Int iHeight,
				Bool b6Param
				);
#endif
};// END CLASS DEFINITION TEncSearch

//! \}

#endif // __TENCSEARCH__