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

/** \file     TEncSbac.h
    \brief    Context-adaptive entropy encoder class (header)
*/

#ifndef __TENCSBAC__
#define __TENCSBAC__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TLibCommon/TComBitStream.h"
#include "TLibCommon/ContextTables.h"
#include "TLibCommon/ContextModel.h"
#include "TLibCommon/ContextModel3DBuffer.h"
#include "TEncEntropy.h"
#include "TEncBinCoder.h"
#include "TEncBinCoderCABAC.h"
#if FAST_BIT_EST
#include "TEncBinCoderCABACCounter.h"
#endif

class TEncTop;

//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// SBAC encoder class
class TEncSbac : public TEncEntropyIf
{
public:
  TEncSbac();
  virtual ~TEncSbac();

  Void  init                   ( TEncBinIf* p )  { m_pcBinIf = p; }
  Void  uninit                 ()                { m_pcBinIf = 0; }

  //  Virtual list
  Void  resetEntropy           (const TComSlice *pSlice);
  SliceType determineCabacInitIdx  (const TComSlice *pSlice);
  Void  setBitstream           ( TComBitIf* p )  { m_pcBitIf = p; m_pcBinIf->init( p ); }

  Void  load                   ( const TEncSbac* pSrc  );
  Void  loadIntraDirMode       ( const TEncSbac* pScr, const ChannelType chType  );
  Void  store                  ( TEncSbac* pDest ) const;
  Void  loadContexts           ( const TEncSbac* pSrc  );
  Void  resetBits              ()                { m_pcBinIf->resetBits(); m_pcBitIf->resetBits(); }
  UInt  getNumberOfWrittenBits ()                { return m_pcBinIf->getNumWrittenBits(); }
  //--SBAC RD

  Void  codeVPS                ( const TComVPS* pcVPS );
  Void  codeSPS                ( const TComSPS* pcSPS     );
  Void  codePPS                ( const TComPPS* pcPPS     );
  Void  codeSliceHeader        ( TComSlice* pcSlice );
  Void  codeTilesWPPEntryPoint ( TComSlice* pSlice );
  Void  codeTerminatingBit     ( UInt uilsLast      );
  Void  codeSliceFinish        ();
#if CNN_BASED_INLOOP_FILTER
  Void  codeCNNLFBlkParam      ( ComponentID compIdx, CNNLFBlkParam* &ctuParam, TComDataCU *pcCU);
#endif
  Void  codeSaoMaxUvlc       ( UInt code, UInt maxSymbol );
  Void  codeSaoMerge         ( UInt  uiCode );
  Void  codeSaoTypeIdx       ( UInt  uiCode);
  Void  codeSaoUflc          ( UInt uiLength, UInt  uiCode );
  Void  codeSAOSign          ( UInt  uiCode);  //<! code SAO offset sign

  Void codeSAOOffsetParam(ComponentID compIdx, SAOOffset& ctbParam, Bool sliceEnabled, const Int channelBitDepth);
  Void codeSAOBlkParam(SAOBlkParam& saoBlkParam, const BitDepths &bitDepths
                    , Bool* sliceEnabled
                    , Bool leftMergeAvail
                    , Bool aboveMergeAvail
                    , Bool onlyEstMergeInfo = false
                    );
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Int  getCtxNumber()     { return m_numContextModels; }
  Void codeCtxUpdateInfo  ( TComSlice* pcSlice,  TComStats* apcStats );
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW
  Void xUpdateWindowSize  ( SliceType eSliceType, Int iQPIdx, TComStats* apcStats);
#endif
#if VCEG_AZ07_INIT_PREVFRAME
  Void loadContextsFromPrev ( TComStats* apcStats, SliceType eSliceType, Int iQPIdx, Bool bFromGloble, Int iQPIdxRst =-1, Bool bAfterLastISlice= false );
#endif

private:
#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
  Void  xWriteTruncBinCode      ( UInt uiSymbol, UInt uiMaxSymbol );
#endif
  Void  xWriteUnarySymbol    ( UInt uiSymbol, ContextModel* pcSCModel, Int iOffset );
  Void  xWriteUnaryMaxSymbol ( UInt uiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol );
  Void  xWriteEpExGolomb     ( UInt uiSymbol, UInt uiCount );
  Void  xWriteCoefRemainExGolomb ( UInt symbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange );

  Void  xCopyFrom            ( const TEncSbac* pSrc );
  Void  xCopyContextsFrom    ( const TEncSbac* pSrc );

protected:
  TComBitIf*    m_pcBitIf;
  TEncBinIf*    m_pcBinIf;

  //--Adaptive loop filter

public:
  Void codeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeSkipFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#if BLOCK_BASED_BACKGROUND_REFERENCE
  Void codeLCUBgFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if BARC_UP_SAMPLE_RDO
  Void codeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if BARC
  Void codeBARCFlag      ( TComDataCU* pcCU, ChannelType eCtype);
#endif
#if VCEG_AZ05_INTRA_MPI
  Void codeMPIIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void codePDPCIdx       ( TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif

#if DIMD_INTRA_PRED
  Void codeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
  Void codeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Void codeROTIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth  );
#if JVET_C0024_QTBT
  Void codeROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth  );
#endif
#endif
#if COM16_C806_OBMC
  Void codeOBMCFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if VCEG_AZ07_IMV
  Void codeiMVFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#if VCEG_AZ06_IC
  Void codeICFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif 
  Void codeMergeFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeMergeIndex    ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#if VCEG_AZ07_FRUC_MERGE
  Void codeFRUCMgrMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx );
#endif
  Void codeSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void codeBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight);
#endif
  Void codeMVPIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );

#if !JVET_C0024_QTBT
  Void codePartSize      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void codePredMode      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeIPCMInfo      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeTransformSubdivFlag ( UInt uiSymbol, UInt uiCtx );
  Void codeQtCbf               ( TComTU & rTu, const ComponentID compID, const Bool lowestLevel );
  Void codeQtRootCbf           ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeQtCbfZero           ( TComTU &rTu, const ChannelType chType );
  Void codeQtRootCbfZero       ( );
  Void codeIntraDirLumaAng     ( TComDataCU* pcCU, UInt absPartIdx, Bool isMultiple
#if VCEG_AZ07_INTRA_65ANG_MODES
    , Int* piModes = NULL, Int  iAboveLeftCase = -1
#endif
    );

  Void codeIntraDirChroma      ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeInterDir            ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeRefFrmIdx           ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );
  Void codeMvd                 ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList );

  Void codeCrossComponentPrediction( TComTU &rTu, ComponentID compID );

  Void codeDeltaQP             ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeChromaQpAdjustment  ( TComDataCU* cu, UInt absPartIdx );

  Void codeLastSignificantXY ( UInt uiPosX, UInt uiPosY, Int width, Int height, ComponentID component, UInt uiScanIdx );
  Void codeCoeffNxN            ( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID
#if VCEG_AZ05_ROT_TR   || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , Int& bCbfCU
#endif  
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );
  Void codeTransformSkipFlags ( TComTU &rTu, ComponentID component );
#if VCEG_AZ08_KLT_COMMON
  Void codeKLTFlags           ( TComTU &rTu, ComponentID component );
#endif

#if RRC_NS
	Void codeNsBlockSize(UInt  uiCode);
	Void codeNsFlag(UInt  uiCode);
	Void codeAdaptiveNsFlag(UInt uiCode);
	Void codeNsFWholeFrameFlag(UInt uiCode);
#endif

#if JVET_E0077_ENHANCED_LM
  Void codeLMModes(TComDataCU* pcCU, UInt uiAbsPartIdx, Int iMode);
#endif

  // -------------------------------------------------------------------------------------------------------------------
  // for RD-optimizatioon
  // -------------------------------------------------------------------------------------------------------------------
#if JVET_D0123_ME_CTX_LUT_BITS
  Void estPuMeBit           (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estMvdBit            (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estMvpIdxBit         (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estRefIdxBit         (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estInterDirBit       (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estMrgFlagBit        (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estMrgIdxBit         (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
#if VCEG_AZ07_FRUC_MERGE
  Void estFrucModeBit       (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);  
#endif
#if COM16_C1016_AFFINE
  Void estAffineFlagBit     (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
#endif
#if MRG_CONSIDER_COMPLEX
  Void estMrgComplexFlagBit (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
  Void estMrgComplexIdxBit  (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
#endif
#if VCEG_AZ07_IMV
  Void estIMVFlagBit        (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac);
#endif
#endif
  Void estBit               (estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , UInt uiScanIdx
#endif
    );
  Void estCBFBit                     ( estBitsSbacStruct* pcEstBitsSbac );
  Void estSignificantCoeffGroupMapBit( estBitsSbacStruct* pcEstBitsSbac, ChannelType chType );
  Void estSignificantMapBit          ( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType );
  Void estLastSignificantPositionBit ( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , UInt uiScanIdx
#endif
    );
  Void estSignificantCoefficientsBit ( estBitsSbacStruct* pcEstBitsSbac, ChannelType chType );

  Void codeExplicitRdpcmMode            ( TComTU &rTu, const ComponentID compID );


  TEncBinIf* getEncBinIf()  { return m_pcBinIf; }

#if ALF_HM3_REFACTOR
#if JVET_C0038_GALF  
  Void codeALFPrevFiltType( UInt uiCode);
  Void codeALFPrevFiltFlag( Int uiCode);
#endif
  Bool  getAlfCtrl             ()                         { return m_bAlfCtrl;          }
  UInt  getMaxAlfCtrlDepth     ()                         { return m_uiMaxAlfCtrlDepth; }
  Void  setAlfCtrl             ( Bool bAlfCtrl          ) { m_bAlfCtrl          = bAlfCtrl;          }
  Void  setMaxAlfCtrlDepth     ( UInt uiMaxAlfCtrlDepth ) { m_uiMaxAlfCtrlDepth = uiMaxAlfCtrlDepth; }
  Void  codeAlfFlag       ( UInt uiCode );
  Void  codeAlfUvlc       ( UInt uiCode );
  Void  codeAlfSvlc       ( Int  uiCode );
  Void  codeAlfCtrlDepth  ( UInt uiMaxTotalCUDepth );
  Void  codeAlfFlagNum    ( UInt uiCode, UInt minValue );
  Void  codeAlfCtrlFlag   ( UInt uiSymbol );
#if !JVET_C0024_QTBT
  Void  codeAlfCtrlFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif
#endif

#if COM16_C806_EMT
  Void codeEmtTuIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void codeEmtCuFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf );
#endif

#if COM16_C1016_AFFINE
  Void codeAffineFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  Void codeAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if MRG_CONSIDER_COMPLEX
  Void codeMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void codeMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#endif

#if RRC_INTRA_SHARPENING
  Void codeIntraSharpFlag     ( UInt uiCode );
#endif

#if RRC_SELECTABLE_IF
  Void codeAdaptiveIfIndex(Int idx);
#endif

private:
  ContextModel         m_contextModels[MAX_NUM_CTX_MOD];
  Int                  m_numContextModels;
  ContextModel3DBuffer m_cCUSplitFlagSCModel;
#if JVET_C0024_QTBT
  ContextModel3DBuffer m_cBTSplitFlagSCModel;
#endif
  ContextModel3DBuffer m_cCUSkipFlagSCModel;
#if BARC_UP_SAMPLE_RDO
  ContextModel3DBuffer m_cCUUpCNNFlagSCModel;
#endif
#if BARC
  ContextModel3DBuffer m_cCUBARCFlagSCModel;
#endif
#if VCEG_AZ05_INTRA_MPI
  ContextModel3DBuffer m_cMPIIdxSCModel;
#endif
#if COM16_C1046_PDPC_INTRA
  ContextModel3DBuffer m_cPDPCIdxSCModel;
#endif

#if DIMD_INTRA_PRED
  ContextModel3DBuffer m_cDIMDEnabledSCModel;
  ContextModel3DBuffer m_cDIMDNoBTFlagSCModel;
#endif

#if RRC_DWDIP
  ContextModel3DBuffer m_cDWDIPIdxSCModel;
#endif
#if RRC_DWDIP_MLR
  ContextModel3DBuffer m_cDWDIP_GPIdxSCModel;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  ContextModel3DBuffer m_cROTidxSCModel;
#endif
#if FAST_EMT_v1
  ContextModel3DBuffer m_cEMTIntraidxSCModel;
#endif
  ContextModel3DBuffer m_cCUMergeFlagExtSCModel;
  ContextModel3DBuffer m_cCUMergeIdxExtSCModel;
#if VCEG_AZ07_FRUC_MERGE
  ContextModel3DBuffer m_cCUFRUCMgrModeSCModel;
  ContextModel3DBuffer m_cCUFRUCMESCModel;
#if THIRD_FRUC_MODE
  ContextModel3DBuffer m_cCUFRUCMESCModel_uni;
#endif
#endif
  ContextModel3DBuffer m_cCUPartSizeSCModel;
  ContextModel3DBuffer m_cCUPredModeSCModel;
  ContextModel3DBuffer m_cCUIntraPredSCModel;
#if MLIP_EXTEND_REF_LINE
  ContextModel3DBuffer m_cCUExternRefSCModel;
#endif
  ContextModel3DBuffer m_cCUChromaPredSCModel;
  ContextModel3DBuffer m_cCUDeltaQpSCModel;
  ContextModel3DBuffer m_cCUInterDirSCModel;
  ContextModel3DBuffer m_cCURefPicSCModel;
  ContextModel3DBuffer m_cCUMvdSCModel;
  ContextModel3DBuffer m_cCUQtCbfSCModel;
  ContextModel3DBuffer m_cCUTransSubdivFlagSCModel;
  ContextModel3DBuffer m_cCUQtRootCbfSCModel;

  ContextModel3DBuffer m_cCUSigCoeffGroupSCModel;
  ContextModel3DBuffer m_cCUSigSCModel;
  ContextModel3DBuffer m_cCuCtxLastX;
  ContextModel3DBuffer m_cCuCtxLastY;
  ContextModel3DBuffer m_cCUOneSCModel;
#if !VCEG_AZ07_CTX_RESIDUALCODING
  ContextModel3DBuffer m_cCUAbsSCModel;
#endif

  ContextModel3DBuffer m_cMVPIdxSCModel;
#if CNN_BASED_INLOOP_FILTER
  ContextModel3DBuffer m_cCUCNNLoopFilterFlagSCModel;
#endif
  ContextModel3DBuffer m_cSaoMergeSCModel;
  ContextModel3DBuffer m_cSaoTypeIdxSCModel;
  ContextModel3DBuffer m_cTransformSkipSCModel;
#if VCEG_AZ08_KLT_COMMON
  ContextModel3DBuffer m_cKLTFlagSCModel;
#endif
  ContextModel3DBuffer m_CUTransquantBypassFlagSCModel;
  ContextModel3DBuffer m_explicitRdpcmFlagSCModel;
  ContextModel3DBuffer m_explicitRdpcmDirSCModel;
  ContextModel3DBuffer m_cCrossComponentPredictionSCModel;

  ContextModel3DBuffer m_ChromaQpAdjFlagSCModel;
  ContextModel3DBuffer m_ChromaQpAdjIdcSCModel;
#if COM16_C806_OBMC
  ContextModel3DBuffer m_cCUOBMCFlagSCModel;
#endif
#if VCEG_AZ07_IMV
  ContextModel3DBuffer m_cCUiMVFlagSCModel;
#endif
#if VCEG_AZ06_IC
  ContextModel3DBuffer m_cCUICFlagSCModel;
#endif
#if ALF_HM3_REFACTOR
  Bool          m_bAlfCtrl;
  UInt          m_uiMaxAlfCtrlDepth;
  ContextModel3DBuffer m_cCUAlfCtrlFlagSCModel;
#if !JVET_C0038_GALF
  ContextModel3DBuffer m_cALFFlagSCModel;
#endif
  ContextModel3DBuffer m_cALFUvlcSCModel;
#if !JVET_C0038_GALF
  ContextModel3DBuffer m_cALFSvlcSCModel;
#endif
#endif

#if COM16_C806_EMT
  ContextModel3DBuffer m_cEmtTuIdxSCModel;
  ContextModel3DBuffer m_cEmtCuFlagSCModel;
#endif

#if COM16_C1016_AFFINE
  ContextModel3DBuffer m_cCUAffineFlagSCModel;
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
	ContextModel3DBuffer m_CUBgFlagSCModel;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  ContextModel3DBuffer m_cCUAffineParamFlagSCModel;
#endif
#if MRG_CONSIDER_COMPLEX
  ContextModel3DBuffer m_cCUMergeComplexFlagExtSCModel;
  ContextModel3DBuffer m_cCUMergeComplexIdxExtSCModel;
#endif
#if RRC_NS
	ContextModel3DBuffer m_cNsFlagSCModel;
#endif
#if RRC_INTRA_SHARPENING
  ContextModel3DBuffer m_cSharpFlagSCModel;
#endif

#if RRC_SELECTABLE_IF
  ContextModel3DBuffer m_cCUAdaptiveIfIndex;
#endif


  UInt m_golombRiceAdaptationStatistics[RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS];

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
public:
  TComStats* m_pcStats;
  TComStats* getStatesHandle ()    { return m_pcStats;             }
  Void setStatesHandle ( TComStats* pcStats ) { m_pcStats = pcStats;}
  Int getContextModelNum()         { return m_numContextModels;    }
#if VCEG_AZ07_BAC_ADAPT_WDOW
  ContextModel*  getContextModel() { return m_contextModels;       }
  TEncBinIf* getBinIf  ()          { return m_pcBinIf;             }
#endif
#endif
};

//! \}

#endif // !defined(AFX_TENCSBAC_H__DDA7CDC4_EDE3_4015_9D32_2156249C82AA__INCLUDED_)
