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

/** \file     TEncCavlc.h
    \brief    CAVLC encoder class (header)
*/

#ifndef __TENCCAVLC__
#define __TENCCAVLC__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComBitStream.h"
#include "TLibCommon/TComRom.h"
#include "TEncEntropy.h"
#include "SyntaxElementWriter.h"

//! \ingroup TLibEncoder
//! \{

class TEncTop;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CAVLC encoder class
class TEncCavlc : public SyntaxElementWriter, public TEncEntropyIf
{
#if ALF_HM3_REFACTOR
protected:
  Bool          m_bAlfCtrl;
  UInt          m_uiMaxAlfCtrlDepth;
#endif

public:
  TEncCavlc();
  virtual ~TEncCavlc();

protected:
  Void codeShortTermRefPicSet              ( const TComReferencePictureSet* pcRPS, Bool calledFromSliceHeader, Int idx );
  Bool findMatchingLTRP ( TComSlice* pcSlice, UInt *ltrpsIndex, Int ltrpPOC, Bool usedFlag );

public:

  Void  resetEntropy          (const TComSlice *pSlice);
  SliceType determineCabacInitIdx  (const TComSlice* /*pSlice*/) { assert(0); return I_SLICE; };

  Void  setBitstream          ( TComBitIf* p )  { m_pcBitIf = p;  }
  Void  resetBits             ()                { m_pcBitIf->resetBits(); }
  UInt  getNumberOfWrittenBits()                { return  m_pcBitIf->getNumberOfWrittenBits();  }
  Void  codeVPS                 ( const TComVPS* pcVPS );
  Void  codeVUI                 ( const TComVUI *pcVUI, const TComSPS* pcSPS );
  Void  codeSPS                 ( const TComSPS* pcSPS );
  Void  codePPS                 ( const TComPPS* pcPPS );
  Void  codeSliceHeader         ( TComSlice* pcSlice );
  Void  codePTL                 ( const TComPTL* pcPTL, Bool profilePresentFlag, Int maxNumSubLayersMinus1);
  Void  codeProfileTier         ( const ProfileTierLevel* ptl, const Bool bIsSubLayer );
  Void  codeHrdParameters       ( const TComHRD *hrd, Bool commonInfPresentFlag, UInt maxNumSubLayersMinus1 );
  Void  codeTilesWPPEntryPoint( TComSlice* pSlice );
  Void  codeTerminatingBit      ( UInt uilsLast );
  Void  codeSliceFinish         ();
#if VCEG_AZ07_BAC_ADAPT_WDOW
  Void codeCtxUpdateInfo        ( TComSlice* pcSlice,  TComStats* apcStats );
  Void xRunCoding               ( Bool * uiCtxMAP, UInt uiNumCtx );
  Void xCtxCodewordCoding       ( Bool * uiCtxMAP, UChar * uiCtxCodeIdx, UInt uiNumCtx );
#endif

  Void codeMVPIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );
#if CNN_BASED_INLOOP_FILTER
	Void codeCNNLFBlkParam(ComponentID compIdx, CNNLFBlkParam* &ctuParam, TComDataCU *pcCU) { printf("only supported in CABAC"); assert(0); exit(-1); }
#endif
  Void codeSAOBlkParam(SAOBlkParam& /*saoBlkParam*/, const BitDepths& /*bitDepths*/, Bool* /*sliceEnabled*/, Bool /*leftMergeAvail*/, Bool /*aboveMergeAvail*/, Bool /*onlyEstMergeInfo*/ = false){printf("only supported in CABAC"); assert(0); exit(-1);}
  Void codeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeSkipFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void codeLCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if BARC_UP_SAMPLE_RDO
	Void codeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if BARC
  Void codeBARCFlag      ( TComDataCU* pcCU, ChannelType eCtype );
#endif
#if VCEG_AZ05_INTRA_MPI
  Void codeMPIIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void codePDPCIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if DIMD_INTRA_PRED
  Void codeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
  Void codeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Void codeROTIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth );
#if JVET_C0024_QTBT
  Void codeROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth );
#endif
#endif
#if VCEG_AZ07_IMV
  Void codeiMVFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if COM16_C806_OBMC
  Void codeOBMCFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if VCEG_AZ06_IC
  Void codeICFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
  Void codeMergeFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeMergeIndex    ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#if VCEG_AZ07_FRUC_MERGE
  Void codeFRUCMgrMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx );
#endif
#if !JVET_C0038_GALF
  Void codeAlfCtrlFlag   ( ComponentID /*component*/, UInt /*code*/ ) {printf("Not supported\n"); assert(0);}
#endif
  Void codeInterModeFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiEncMode );
  Void codeSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void codeBTSplitMode   (TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight);
#else
  Void codePartSize      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void codePredMode      ( TComDataCU* pcCU, UInt uiAbsPartIdx );

  Void codeIPCMInfo      ( TComDataCU* pcCU, UInt uiAbsPartIdx );

  Void codeTransformSubdivFlag( UInt uiSymbol, UInt uiCtx );
  Void codeQtCbf         ( TComTU &rTu, const ComponentID compID, const Bool lowestLevel );
  Void codeQtRootCbf     ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeQtCbfZero     ( TComTU &rTu, const ChannelType chType );
  Void codeQtRootCbfZero ( );
  Void codeIntraDirLumaAng( TComDataCU* pcCU, UInt absPartIdx, Bool isMultiple
#if VCEG_AZ07_INTRA_65ANG_MODES
    , Int* piModes = NULL, Int  iAboveLeftCase = -1
#endif
    );
  Void codeIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeInterDir      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeRefFrmIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void codeMvd           ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void codeCrossComponentPrediction( TComTU &rTu, ComponentID compID );

  Void codeDeltaQP       ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeChromaQpAdjustment( TComDataCU* pcCU, UInt uiAbsPartIdx );

  Void codeCoeffNxN      ( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID 
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , Int& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );
  Void codeTransformSkipFlags ( TComTU &rTu, ComponentID component );

#if VCEG_AZ08_KLT_COMMON
  Void codeKLTFlags      (TComTU &rTu, ComponentID component);
#endif

#if JVET_D0123_ME_CTX_LUT_BITS
  Void estPuMeBit        ( estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
#endif

  Void estBit            ( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , UInt uiScanIdx
#endif
    );

  Void xCodePredWeightTable          ( TComSlice* pcSlice );

  Void codeScalingList  ( const TComScalingList &scalingList );
  Void xCodeScalingList ( const TComScalingList* scalingList, UInt sizeId, UInt listId);

  Void codeExplicitRdpcmMode( TComTU &rTu, const ComponentID compID );
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  TComStats* m_pcStats;
  TComStats* getStatesHandle ( )                    { return m_pcStats;    }
  Void setStatesHandle       ( TComStats* pcStats ) { m_pcStats = pcStats; }
#endif

#if ALF_HM3_REFACTOR
  Void xWriteUnaryMaxSymbol( UInt uiSymbol, UInt uiMaxSymbol );
  Bool getAlfCtrl() {return m_bAlfCtrl;}
  UInt getMaxAlfCtrlDepth() {return m_uiMaxAlfCtrlDepth;}
  Void setAlfCtrl(Bool bAlfCtrl) {m_bAlfCtrl = bAlfCtrl;}
  Void setMaxAlfCtrlDepth(UInt uiMaxAlfCtrlDepth) {m_uiMaxAlfCtrlDepth = uiMaxAlfCtrlDepth;}
  Void codeAlfFlag       ( UInt uiCode );
#if JVET_C0038_GALF  
  Void xWriteTruncBinCode(UInt uiSymbol, UInt uiMaxSymbol);
  Void codeALFPrevFiltType( UInt uiCode);
  Void codeALFPrevFiltFlag( Int uiCode);
#endif
  Void codeAlfUvlc       ( UInt uiCode );
  Void codeAlfSvlc       ( Int   iCode );
  Void codeAlfCtrlDepth  ( UInt uiMaxTotalCUDepth );
#if !JVET_C0024_QTBT
  Void codeAlfCtrlFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
  Void codeAlfFlagNum    ( UInt uiCode, UInt minValue );
  Void codeAlfCtrlFlag   ( UInt uiSymbol );
#endif

#if COM16_C806_EMT
  Void codeEmtTuIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
  Void codeEmtCuFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf );
#endif

#if COM16_C1016_AFFINE
  Void codeAffineFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif

#if MRG_CONSIDER_COMPLEX
  Void codeMergeComplexFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
  Void codeMergeComplexIndex(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  Void codeAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif

#if RRC_NS
	Void codeNsBlockSize(UInt uiCode);
	Void codeNsFlag(UInt uiCode);
	Void codeAdaptiveNsFlag(UInt uiCode);
	Void codeNsFWholeFrameFlag(UInt uiCode);
#endif

#if RRC_INTRA_SHARPENING
  Void codeIntraSharpFlag    ( UInt uiCode );
#endif

#if RRC_SELECTABLE_IF
  Void codeAdaptiveIfIndex(int idx);
#endif
};

//! \}

#endif // !defined(AFX_TENCCAVLC_H__EE8A0B30_945B_4169_B290_24D3AD52296F__INCLUDED_)

