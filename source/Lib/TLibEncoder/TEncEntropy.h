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

/** \file     TEncEntropy.h
    \brief    entropy encoder class (header)
*/

#ifndef __TENCENTROPY__
#define __TENCENTROPY__

#include "TLibCommon/TComSlice.h"
#include "TLibCommon/TComDataCU.h"
#include "TLibCommon/TComBitStream.h"
#include "TLibCommon/ContextModel.h"
#include "TLibCommon/TComPic.h"
#include "TLibCommon/TComTrQuant.h"
#include "TLibCommon/TComSampleAdaptiveOffset.h"
#include "TLibCommon/TComChromaFormat.h"
#if ALF_HM3_REFACTOR
#include "TLibCommon//TComAdaptiveLoopFilter.h"
#endif

class TEncSbac;
class TEncCavlc;
class SEI;
#if JVET_D0123_ME_CTX_LUT_BITS
typedef struct
{
  Int mvdBits[NUM_MV_RES_CTX][2 /*Flag = [0|1]*/];
  Int mvpIdxBits[NUM_MVP_IDX_CTX][2 /*Flag = [0|1]*/];
  Int mrgFlagBits[NUM_MERGE_FLAG_EXT_CTX][2 /*Flag = [0|1]*/];
  Int mrgIdxBits[NUM_MERGE_IDX_EXT_CTX][2 /*Flag = [0|1]*/];
  Int refIdxBits[NUM_REF_NO_CTX][2 /*Flag = [0|1]*/];
#if VCEG_AZ07_FRUC_MERGE
  Int frucMrgBits[NUM_FRUCMGRMODE_CTX][2 /*Flag = [0|1]*/];
  Int frucMeBits [NUM_FRUCME_CTX][2 /*Flag = [0|1]*/];
#if THIRD_FRUC_MODE
  Int frucMeBitsUni[NUM_FRUCME_CTX][2 /*Flag = [0|1]*/];
#endif
#endif
  Int interDirBits[NUM_INTER_DIR_CTX][2 /*Flag = [0|1]*/];
#if COM16_C1016_AFFINE
  Int affineFlagBits[NUM_AFFINE_FLAG_CTX][2 /*Flag = [0|1]*/];
#endif
#if MRG_CONSIDER_COMPLEX
  Int mrgComplexFlagBits[NUM_MERGE_FLAG_EXT_CTX][2 /*Flag = [0|1]*/];
  Int mrgComplexIdxBits[NUM_MERGE_COMPLEX_IDX_EXT_CTX][2 /*Flag = [0|1]*/];
#endif
#if VCEG_AZ07_IMV
  Int iMVFlagBits[NUM_IMV_FLAG_CTX][2/*Flag = [0|1]*/];
#endif
} estPuMeBitsSbacStruct;
#endif
// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// entropy encoder pure class
class TEncEntropyIf
{
public:
  virtual Void  resetEntropy          (const TComSlice *pSlice)                = 0;
  virtual SliceType determineCabacInitIdx (const TComSlice *pSlice)                = 0;
  virtual Void  setBitstream          ( TComBitIf* p )  = 0;
  virtual Void  resetBits             ()                = 0;
  virtual UInt  getNumberOfWrittenBits()                = 0;

  virtual Void  codeVPS                 ( const TComVPS* pcVPS )                                      = 0;
  virtual Void  codeSPS                 ( const TComSPS* pcSPS )                                      = 0;
  virtual Void  codePPS                 ( const TComPPS* pcPPS )                                      = 0;
  virtual Void  codeSliceHeader         ( TComSlice* pcSlice )                                  = 0;

  virtual Void  codeTilesWPPEntryPoint  ( TComSlice* pSlice )     = 0;
  virtual Void  codeTerminatingBit      ( UInt uilsLast )                                       = 0;
  virtual Void  codeSliceFinish         ()                                                      = 0;
  virtual Void  codeMVPIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList ) = 0;

public:
  virtual Void codeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeSkipFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#if BLOCK_BASED_BACKGROUND_REFERENCE
	virtual Void codeLCUBgFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
#if BARC
  virtual Void codeBARCFlag      ( TComDataCU* pcCU, ChannelType eCtype ) = 0;
#endif
#if BARC_UP_SAMPLE_RDO
  virtual Void codeUpCNNFlag     (TComDataCU* pcCU, UInt uiAbsPartIdx) = 0;
#endif
#if VCEG_AZ05_INTRA_MPI
  virtual Void codeMPIIdx        (TComDataCU* pcCU, UInt uiAbsPartIdx)   = 0;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  virtual Void codePDPCIdx       (TComDataCU* pcCU, UInt uiAbsPartIdx) = 0;
#endif
#if DIMD_INTRA_PRED
  virtual Void codeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth) = 0;
  virtual Void codeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth) = 0;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
   virtual Void codeROTIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth ) = 0;
#if JVET_C0024_QTBT
   virtual Void codeROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth ) = 0;
#endif
#endif
#if VCEG_AZ07_IMV
  virtual Void codeiMVFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
#if COM16_C806_OBMC
  virtual Void codeOBMCFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
#if VCEG_AZ06_IC
  virtual Void codeICFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
  virtual Void codeMergeFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeMergeIndex    ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#if CNN_BASED_INLOOP_FILTER
  virtual Void codeCNNLFBlkParam (ComponentID compIdx, CNNLFBlkParam* &ctuParam, TComDataCU *pcCU) = 0;
#endif
#if VCEG_AZ07_FRUC_MERGE
  virtual Void codeFRUCMgrMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx ) = 0;
#endif
  virtual Void codeSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#if JVET_C0024_QTBT
  virtual Void codeBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight) = 0;
#else
  virtual Void codePartSize      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
  virtual Void codePredMode      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;

  virtual Void codeIPCMInfo      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;

  virtual Void codeTransformSubdivFlag( UInt uiSymbol, UInt uiCtx ) = 0;
  virtual Void codeQtCbf         ( TComTU &rTu, const ComponentID compID, const Bool lowestLevel ) = 0;
  virtual Void codeQtRootCbf     ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeQtCbfZero     ( TComTU &rTu, const ChannelType chType ) = 0;
  virtual Void codeQtRootCbfZero ( ) = 0;
  virtual Void codeIntraDirLumaAng( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool isMultiplePU 
#if VCEG_AZ07_INTRA_65ANG_MODES
    , Int* piModes = NULL, Int iAboveLeftCase = -1
#endif
    ) = 0;

  virtual Void codeIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeInterDir      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeRefFrmIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )      = 0;
  virtual Void codeMvd           ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )      = 0;
  virtual Void codeCrossComponentPrediction( TComTU &rTu, ComponentID compID ) = 0;

  virtual Void codeDeltaQP       ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeChromaQpAdjustment( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeCoeffNxN      ( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID 
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , int& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    ) = 0;
  virtual Void codeTransformSkipFlags ( TComTU &rTu, ComponentID component ) = 0;
#if VCEG_AZ08_KLT_COMMON
  virtual Void codeKLTFlags      (TComTU &rTu, ComponentID component) = 0;
#endif
  virtual Void codeSAOBlkParam   (SAOBlkParam& saoBlkParam, const BitDepths &bitDepths, Bool* sliceEnabled, Bool leftMergeAvail, Bool aboveMergeAvail, Bool onlyEstMergeInfo = false)    =0;
  virtual Void estBit               (estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , UInt uiScanIdx
#endif
    ) = 0;
#if JVET_D0123_ME_CTX_LUT_BITS
  virtual Void estPuMeBit           (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac) = 0;
#endif
  virtual Void codeExplicitRdpcmMode ( TComTU &rTu, const ComponentID compID ) = 0;

  virtual ~TEncEntropyIf() {}

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  virtual TComStats* getStatesHandle () = 0;
  virtual Void setStatesHandle   ( TComStats* pcStats ) = 0;
#if VCEG_AZ07_BAC_ADAPT_WDOW 
  virtual Void codeCtxUpdateInfo ( TComSlice* pcSlice,  TComStats* apcStats ) = 0;
#endif
#endif

#if ALF_HM3_REFACTOR
#if JVET_C0038_GALF  
  virtual Void xWriteTruncBinCode   (UInt iSymbol, UInt iMaxSymbol) = 0;
  virtual Void codeALFPrevFiltType( UInt uiCode) = 0;
  virtual Void codeALFPrevFiltFlag( Int uiCode)  = 0;
#endif
  virtual Bool getAlfCtrl()                = 0;
  virtual UInt getMaxAlfCtrlDepth()                = 0;
  virtual Void setAlfCtrl(Bool bAlfCtrl)                = 0;
  virtual Void setMaxAlfCtrlDepth(UInt uiMaxAlfCtrlDepth)                = 0;
  virtual Void codeAlfCtrlDepth     ( UInt uiMaxTotalCUDepth ) = 0;
#if !JVET_C0024_QTBT
  virtual Void codeAlfCtrlFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
  virtual Void codeAlfFlag          ( UInt uiCode ) = 0;
  virtual Void codeAlfUvlc          ( UInt uiCode ) = 0;
  virtual Void codeAlfSvlc          ( Int   iCode ) = 0;
  virtual Void codeAlfFlagNum       ( UInt uiCode, UInt minValue ) = 0;
  virtual Void codeAlfCtrlFlag      ( UInt uiSymbol ) = 0;
#endif

#if COM16_C806_EMT
  virtual Void codeEmtTuIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
  virtual Void codeEmtCuFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf ) = 0;
#endif

#if COM16_C1016_AFFINE
  virtual Void codeAffineFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
#if MRG_CONSIDER_COMPLEX
  virtual Void codeMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
  virtual Void codeMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  virtual Void codeAffineParamFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx ) = 0;
#endif

#if RRC_NS
	virtual Void codeNsBlockSize(UInt   uiCode) = 0;
	virtual Void codeNsFlag(UInt   uiCode) = 0;
	virtual Void codeAdaptiveNsFlag(UInt   uiCode) = 0;
	virtual Void codeNsFWholeFrameFlag(UInt   uiCode) = 0;
#endif

#if RRC_INTRA_SHARPENING
  virtual Void codeIntraSharpFlag     ( UInt uiCode) = 0;
#endif

#if RRC_SELECTABLE_IF
  virtual Void codeAdaptiveIfIndex(Int idx) = 0;
#endif

};

/// entropy encoder class
class TEncEntropy
{
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
private:
  TComStats*                      m_pcStats;
#endif

public:
  Void    setEntropyCoder           ( TEncEntropyIf* e );
  Void    setBitstream              ( TComBitIf* p )          { m_pcEntropyCoderIf->setBitstream(p);  }
  Void    resetBits                 ()                        { m_pcEntropyCoderIf->resetBits();      }
  UInt    getNumberOfWrittenBits    ()                        { return m_pcEntropyCoderIf->getNumberOfWrittenBits(); }
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Void    resetEntropy              (const TComSlice *pSlice)                        
  {
    m_pcEntropyCoderIf->setStatesHandle (m_pcStats);
    m_pcEntropyCoderIf->resetEntropy(pSlice);  
  }
#else
  Void    resetEntropy              (const TComSlice *pSlice) { m_pcEntropyCoderIf->resetEntropy(pSlice);  }
#endif

  SliceType determineCabacInitIdx   (const TComSlice *pSlice) { return m_pcEntropyCoderIf->determineCabacInitIdx(pSlice); }

  Void    encodeSliceHeader         ( TComSlice* pcSlice );
  Void    encodeTilesWPPEntryPoint( TComSlice* pSlice );
  Void    encodeTerminatingBit      ( UInt uiIsLast );
  Void    encodeSliceFinish         ();
  TEncEntropyIf*      m_pcEntropyCoderIf;
#if VCEG_AZ07_BAC_ADAPT_WDOW
  TEncSbac* getCABACCoder  ( ) { return (TEncSbac*)m_pcEntropyCoderIf; }
  Void encodeCtxUpdateInfo ( TComSlice* pcSlice,  TComStats* apcStats );
#endif

public:
  Void encodeVPS               ( const TComVPS* pcVPS);
  // SPS
  Void encodeSPS               ( const TComSPS* pcSPS );
  Void encodePPS               ( const TComPPS* pcPPS );
  Void encodeSplitFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRD = false );
#if JVET_C0024_QTBT
  Void encodeBTSplitMode       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, Bool bRD = false );
#endif
  Void encodeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
  Void encodeSkipFlag          ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void encodeLCUBgFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#endif
#if BARC_UP_SAMPLE_RDO
	Void encodeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false);
#endif
#if BARC
  Void encodeBARCFlag          ( TComDataCU* pcCU, ChannelType eCtype );
#endif
#if COM16_C806_OBMC
  Void encodeOBMCFlag          ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#endif
#if VCEG_AZ07_IMV
  Void encodeiMVFlag           ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if VCEG_AZ06_IC
  Void encodeICFlag            ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if JVET_D0123_ME_CTX_LUT_BITS
  Void encodePuMotionInfo     (TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPartIdx);
#endif
  Void encodePUWise       ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void encodeInterDirPU   ( TComDataCU* pcSubCU, UInt uiAbsPartIdx  );
  Void encodeRefFrmIdxPU  ( TComDataCU* pcSubCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void encodeMvdPU        ( TComDataCU* pcSubCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void encodeMVPIdxPU     ( TComDataCU* pcSubCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void encodeMergeFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void encodeMergeIndex   ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#if VCEG_AZ05_INTRA_MPI
  Void encodeMPIIdx       ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false);
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void encodePDPCIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false);
#endif
#if DIMD_INTRA_PRED
  Void encodeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, Bool bRD = false);
  Void encodeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, Bool bRD = false);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    Void encodeROTIdx    ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth, Bool bRD = false );
#if JVET_C0024_QTBT
    Void encodeROTIdxChroma( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth, Bool bRD = false );
#endif
#endif
#if VCEG_AZ07_FRUC_MERGE
  Void encodeFRUCMgrMode  ( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx );
#endif
  Void encodePredMode          ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#if !JVET_C0024_QTBT
  Void encodePartSize          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRD = false );
#endif
  Void encodeIPCMInfo          ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
  Void encodePredInfo          ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void encodeIntraDirModeLuma  ( TComDataCU* pcCU, UInt absPartIdx, Bool isMultiplePU = false 
#if VCEG_AZ07_INTRA_65ANG_MODES
    , Int* piModes = NULL, Int iAboveLeftCase = -1
#endif
    );

  Void encodeIntraDirModeChroma( TComDataCU* pcCU, UInt uiAbsPartIdx );

  Void encodeTransformSubdivFlag( UInt uiSymbol, UInt uiCtx );
  Void encodeQtCbf             ( TComTU &rTu, const ComponentID compID, const Bool lowestLevel );

  Void encodeQtCbfZero         ( TComTU &rTu, const ChannelType chType );
  Void encodeQtRootCbfZero     ( );
  Void encodeQtRootCbf         ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void encodeQP                ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
  Void encodeChromaQpAdjustment ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );

  Void encodeCrossComponentPrediction( TComTU &rTu, ComponentID compID );

private:
#if JVET_C0024_QTBT
  Void xEncodeTransform( Bool& bCodeDQP, Bool& codeChromaQpAdj, TComTU &rTu, ComponentID compID);
#else
  Void xEncodeTransform        ( Bool& bCodeDQP, Bool& codeChromaQpAdj, TComTU &rTu
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , Int& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );
#endif

public:
  Void encodeCoeff             ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool& bCodeDQP, Bool& codeChromaQpAdj
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , Int& bNonZeroCoeff
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );

  Void encodeCoeffNxN         ( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID );

  Void estimateBit             ( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , UInt uiScanIdx
#endif
    );
#if JVET_D0123_ME_CTX_LUT_BITS
  Void estimatePuMeBit        ( estPuMeBitsSbacStruct* pcMePuEstBitsSbac );
#endif
#if CNN_BASED_INLOOP_FILTER
	Void encodeCNNLFBlkParam(ComponentID compIdx, CNNLFBlkParam* &ctuParam, TComDataCU *pcCU) { m_pcEntropyCoderIf->codeCNNLFBlkParam(compIdx, ctuParam, pcCU); }
#endif
  Void encodeSAOBlkParam(SAOBlkParam& saoBlkParam, const BitDepths &bitDepths, Bool* sliceEnabled, Bool leftMergeAvail, Bool aboveMergeAvail){m_pcEntropyCoderIf->codeSAOBlkParam(saoBlkParam, bitDepths, sliceEnabled, leftMergeAvail, aboveMergeAvail, false);}

  static Int countNonZeroCoeffs( TCoeff* pcCoef, UInt uiSize );

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Void setStatsHandle       ( TComStats*  pcStats)  { m_pcStats=pcStats; }
  TComStats* getStatsHandle ()                      { return m_pcStats;  }
#endif 

#if ALF_HM3_REFACTOR
  Void encodeAlfParam(ALFParam* pAlfParam, UInt uiMaxTotalCUDepth
#if FIX_TICKET12
        ,const TComSlice * pSlice
#endif
    );
#if !JVET_C0024_QTBT
  Void encodeAlfCtrlFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#endif
  Void encodeAlfCtrlParam      ( ALFParam *pAlfParam );
  Bool getAlfCtrl() {return m_pcEntropyCoderIf->getAlfCtrl();}
  UInt getMaxAlfCtrlDepth() {return m_pcEntropyCoderIf->getMaxAlfCtrlDepth();}
  Void setAlfCtrl(Bool bAlfCtrl) {m_pcEntropyCoderIf->setAlfCtrl(bAlfCtrl);}
  Void setMaxAlfCtrlDepth(UInt uiMaxAlfCtrlDepth) {m_pcEntropyCoderIf->setMaxAlfCtrlDepth(uiMaxAlfCtrlDepth);}
  Void codeAuxCountBit(ALFParam* pAlfParam, Int64* ruiRate, const TComSlice * pSlice);
  Void codeFiltCountBit(ALFParam* pAlfParam, Int64* ruiRate, const TComSlice * pSlice);
  Void codeAux (ALFParam* pAlfParam);
  Void codeFilt (ALFParam* pAlfParam);
  Int codeFilterCoeff(ALFParam* ALFp
#if JVET_C0038_GALF
    , Bool bChroma = false
#endif
    );
#if JVET_C0038_GALF
  Void codeFilterCoeffForce0(ALFParam* ALFp);
  Int writeFilterCodingParams(Int minKStart, Int maxScanVal, Int kMinTab[], Bool forceCoeff0, Int filters_per_group, Bool codedVarBins[]);
  Int writeFilterCoeffs(Int sqrFiltLength, Int filters_per_group, const Int pDepthInt[], Int **FilterCoeff, Int kMinTab[], Bool codedVarBins[]);
#else
  Int writeFilterCodingParams(int minKStart, int maxScanVal, int kMinTab[]);
  Int writeFilterCoeffs(int sqrFiltLength, int filters_per_group, const int pDepthInt[], 
    int **FilterCoeff, int kMinTab[]);
#endif
#if JVET_C0038_GALF
  Int writeFilterCoeffsForChroma(Int sqrFiltLength, const Int pDepthInt[], Int *FilterCoeff, Int kMinTab[]);
#endif
  Int golombEncode(int coeff, int k);
  Int lengthGolomb(int coeffVal, int k);
#endif
#if COM16_C806_EMT
  Void encodeEmtCuFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bCodeCuFlag );
#endif

#if COM16_C1016_AFFINE
  Void encodeAffineFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPuIdx );
#endif
#if MRG_CONSIDER_COMPLEX
  Void encodeMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void encodeMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD = false );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  Void encodeAffineParamFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPuIdx );
#endif

#if RRC_INTRA_SHARPENING
  void encodeIntraSharpFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif

#if RRC_SELECTABLE_IF
  Void encodeAdaptiveIfIndex(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif

#if RRC_NS
	Void encodeNsBlockSize(Int size) { { m_pcEntropyCoderIf->codeNsBlockSize(size); } }
	Void encodeNsFlag(Int flag) { m_pcEntropyCoderIf->codeNsFlag(flag); }
	Void encodeAdaptiveNsFlag(Int flag) { m_pcEntropyCoderIf->codeAdaptiveNsFlag(flag); }
	Void encodeNsFWholeFrameFlag(Int flag) { m_pcEntropyCoderIf->codeNsFWholeFrameFlag(flag); }
#endif

};// END CLASS DEFINITION TEncEntropy

//! \}

#endif // __TENCENTROPY__

