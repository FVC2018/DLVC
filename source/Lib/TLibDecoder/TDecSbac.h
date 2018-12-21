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

/** \file     TDecSbac.h
    \brief    SBAC decoder class (header)
*/

#ifndef __TDECSBAC__
#define __TDECSBAC__


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "TDecEntropy.h"
#include "TDecBinCoder.h"
#include "TLibCommon/ContextTables.h"
#include "TLibCommon/ContextModel.h"
#include "TLibCommon/ContextModel3DBuffer.h"

//! \ingroup TLibDecoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

class SEImessages;

/// SBAC decoder class
class TDecSbac : public TDecEntropyIf
{
public:
  TDecSbac();
  virtual ~TDecSbac();
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  TComStats* m_pcStats;
  TComStats* getStatesHandle      ()                    { return m_pcStats;          }
  Void setStatesHandle            ( TComStats* pcStats) { m_pcStats = pcStats;       }
  Int  getCtxNumber               ()                    { return m_numContextModels; }
#endif

  Void  init                      ( TDecBinIf* p )    { m_pcTDecBinIf = p; }
  Void  uninit                    (              )    { m_pcTDecBinIf = 0; }

  Void load                       ( const TDecSbac* pSrc );
  Void loadContexts               ( const TDecSbac* pSrc );
  Void xCopyFrom                  ( const TDecSbac* pSrc );
  Void xCopyContextsFrom          ( const TDecSbac* pSrc );

  Void  resetEntropy (TComSlice* pSlice );
  Void  setBitstream              ( TComInputBitstream* p  ) { m_pcBitstream = p; m_pcTDecBinIf->init( p ); }
  Void  parseVPS                  ( TComVPS* /*pcVPS*/ ) {}
  Void  parseSPS                  ( TComSPS* /*pcSPS*/ ) {}
  Void  parsePPS                  ( TComPPS* /*pcPPS*/ ) {}
#if VCEG_AZ07_BAC_ADAPT_WDOW
  Void parseCtxUpdateInfo         ( TComSlice*& rpcSlice,  TComStats* apcStats )   {}
  Void xUpdateWindowSize          ( SliceType eSliceType, Int uiQPIdx, TComStats* apcStats );  
#endif

  Void  parseSliceHeader          ( TComSlice* /*pcSlice*/, ParameterSetManager* /*parameterSetManager*/, const Int /*prevTid0POC*/) {}
  Void  parseTerminatingBit       ( UInt& ruiBit );
  Void  parseRemainingBytes       ( Bool noTrailingBytesExpected);
  Void  parseMVPIdx               ( Int& riMVPIdx          );
  Void  parseSaoMaxUvlc           ( UInt& val, UInt maxSymbol );
  Void  parseSaoMerge             ( UInt&  ruiVal   );
  Void  parseSaoTypeIdx           ( UInt&  ruiVal  );
  Void  parseSaoUflc              ( UInt uiLength, UInt& ruiVal     );
  Void parseSAOBlkParam           (SAOBlkParam& saoBlkParam, Bool* sliceEnabled, Bool leftMergeAvail, Bool aboveMergeAvail, const BitDepths &bitDepths);
  Void parseSaoSign               (UInt& val);
#if CNN_BASED_INLOOP_FILTER
	Void parseCNNLoopFilterTypeIdx  (TComDataCU *pcCU, UInt&  ruiVal, CNNLFBlkParam* &ctuParam, ComponentID compIdx);
	Void parseConvNetBlkParam       (CNNLFBlkParam*& convNetBlkParam, TComDataCU *pcCU);
#endif
private:
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
  Void  xReadTruncBinCode   (UInt& ruiSymbol, UInt uiMaxSymbol, const class TComCodingStatisticsClassType &whichStat);
#endif
  Void  xReadUnarySymbol    ( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, const class TComCodingStatisticsClassType &whichStat );
  Void  xReadUnaryMaxSymbol ( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol, const class TComCodingStatisticsClassType &whichStat );
  Void  xReadEpExGolomb     ( UInt& ruiSymbol, UInt uiCount, const class TComCodingStatisticsClassType &whichStat );
  Void  xReadCoefRemainExGolomb ( UInt &rSymbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange, const class TComCodingStatisticsClassType &whichStat );
#else
#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
 Void  xReadTruncBinCode       (UInt& ruiSymbol, UInt uiMaxSymbol);
#endif
  Void  xReadUnarySymbol    ( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset );
  Void  xReadUnaryMaxSymbol ( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol );
  Void  xReadEpExGolomb     ( UInt& ruiSymbol, UInt uiCount );
  Void  xReadCoefRemainExGolomb ( UInt &rSymbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange );
#endif
private:
  TComInputBitstream* m_pcBitstream;
  TDecBinIf*        m_pcTDecBinIf;

public:

  Void parseSkipFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if BARC_UP_SAMPLE_RDO
  Void parseUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if BARC
  Void parseBARCFlag      ( TComDataCU* pcCU, ChannelType eCtype );
#endif
#if VCEG_AZ07_IMV
  Void parseiMVFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if COM16_C806_OBMC
  Void parseOBMCFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if VCEG_AZ06_IC
  Void parseICFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void parseCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void parseCUBgFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void parseSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void parseBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
    , UInt uiBTDepth
#endif
    );
#endif
#if VCEG_AZ05_INTRA_MPI
  Void parseMPIIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void parsePDPCIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if DIMD_INTRA_PRED
  Void parseDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
  Void parseDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Void parseROTIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void parseROTIdxChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#endif
  Void parseMergeFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void parseMergeIndex    ( TComDataCU* pcCU, UInt& ruiMergeIndex );
#if VCEG_AZ07_FRUC_MERGE
  Void parseFRUCMgrMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
#endif
#if !JVET_C0024_QTBT
  Void parsePartSize      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void parsePredMode      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void parseIntraDirLumaAng( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void parseInterDir      ( TComDataCU* pcCU, UInt& ruiInterDir, UInt uiAbsPartIdx );
  Void parseRefFrmIdx     ( TComDataCU* pcCU, Int& riRefFrmIdx, RefPicList eRefList );
  Void parseMvd           ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList );
  Void parseCrossComponentPrediction ( class TComTU &rTu, ComponentID compID );

  Void parseTransformSubdivFlag( UInt& ruiSubdivFlag, UInt uiLog2TransformBlockSize );
  Void parseQtCbf         ( TComTU &rTu, const ComponentID compID, const Bool lowestLevel );
  Void parseQtRootCbf     ( UInt uiAbsPartIdx, UInt& uiQtRootCbf );

  Void parseDeltaQP       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseChromaQpAdjustment( TComDataCU* cu, UInt absPartIdx, UInt depth );

  Void parseIPCMInfo      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);

  Void parseLastSignificantXY( UInt& uiPosLastX, UInt& uiPosLastY, Int width, Int height, ComponentID component, UInt uiScanIdx );
  Void parseCoeffNxN      ( class TComTU &rTu, ComponentID compID  
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );
  Void parseTransformSkipFlags ( class TComTU &rTu, ComponentID component );

#if VCEG_AZ08_KLT_COMMON
  Void parseKLTFlags      ( TComTU &rTu, ComponentID component);
#endif

  Void  parseScalingList ( TComScalingList* /*scalingList*/ ) {}

#if ALF_HM3_REFACTOR
  Void  parseAlfFlag          ( UInt& ruiVal           );
#if JVET_C0038_GALF
  Void parseALFTruncBinVal  ( UInt& ruiSymbol, UInt uiMaxSymbol );
  Void parseALFPrevFiltType ( UInt& uiCode );
  Void parseALFPrevFiltFlag ( UInt& uiCode );
#endif
  Void  parseAlfUvlc          ( UInt& ruiVal           );
  Void  parseAlfSvlc          ( Int&  riVal            );
  Void  parseAlfCtrlDepth     ( UInt& ruiAlfCtrlDepth , UInt uiMaxTotalCUDepth ); 
#if !JVET_C0024_QTBT
  Void  parseAlfCtrlFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth , UInt uiMaxAlfCtrlDepth );
#endif
  Void  parseAlfFlagNum       ( UInt& ruiVal, UInt minValue, UInt depth );
  Void  parseAlfCtrlFlag      ( UInt &ruiAlfCtrlFlag );
#endif
#if RRC_SELECTABLE_IF
  Int parseAdaptiveIfIndex();
#endif

#if RRC_NS
	UInt parseNsFlag();
	UInt parseAdaptiveNsFlag();
	UInt parseNsFWholeFrameFlag();
	Int  parseNsBlockSize();
#endif

#if RRC_INTRA_SHARPENING
  UInt parseIntraSharpFlag();
#endif

#if COM16_C806_EMT
  Void parseEmtTuIdx          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseEmtCuFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf );
#endif

#if VCEG_AZ07_INIT_PREVFRAME
  Void  loadContextsFromPrev  ( TComStats* apcStats, SliceType eSliceType, Int iQPIdx, Bool bFromGloble, Int iQPIdxRst =-1, Bool bAfterLastISlice= false );
#endif

#if COM16_C1016_AFFINE
  Void  parseAffineFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
  Void  parseAffineMvd      ( TComDataCU* pcCU, UInt uiAbsPartAddr, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  Void  parseAffineParamFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
#endif

#if MRG_CONSIDER_COMPLEX
  Void parseMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void parseMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt& ruiMergeComplexIndex );
#endif

  Void  parseExplicitRdpcmMode( TComTU &rTu, ComponentID compID );


#if JVET_E0077_ENHANCED_LM
  Int  parseLMMode(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
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
};

//! \}

#endif // !defined(AFX_TDECSBAC_H__CFCAAA19_8110_47F4_9A16_810C4B5499D5__INCLUDED_)
