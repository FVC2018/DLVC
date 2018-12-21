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

/** \file     TDecCAVLC.h
    \brief    CAVLC decoder class (header)
*/

#ifndef __TDECCAVLC__
#define __TDECCAVLC__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TDecEntropy.h"
#include "SyntaxElementParser.h"

//! \ingroup TLibDecoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CAVLC decoder class
class TDecCavlc : public SyntaxElementParser, public TDecEntropyIf
{
public:
  TDecCavlc();
  virtual ~TDecCavlc();

protected:

  Void  parseShortTermRefPicSet            (TComSPS* pcSPS, TComReferencePictureSet* pcRPS, Int idx);

public:

  /// rest entropy coder by intial QP and IDC in CABAC
  Void  resetEntropy        ( TComSlice* /*pcSlice*/  )     { assert(0); };
  Void  setBitstream        ( TComInputBitstream* p )   { m_pcBitstream = p; }
  Void  parseTransformSubdivFlag( UInt& ruiSubdivFlag, UInt uiLog2TransformBlockSize );
  Void  parseQtCbf          ( class TComTU &rTu, const ComponentID compID, const Bool lowestLevel );
  Void  parseQtRootCbf      ( UInt uiAbsPartIdx, UInt& uiQtRootCbf );
  Void  parseVPS            ( TComVPS* pcVPS );
  Void  parseSPS            ( TComSPS* pcSPS );
  Void  parsePPS            ( TComPPS* pcPPS );
  Void  parseVUI            ( TComVUI* pcVUI, TComSPS* pcSPS );
  Void  parseSEI            ( SEIMessages& );
  Void  parsePTL            ( TComPTL *rpcPTL, Bool profilePresentFlag, Int maxNumSubLayersMinus1 );
  Void  parseProfileTier    (ProfileTierLevel *ptl, const Bool bIsSubLayer);
  Void  parseHrdParameters  (TComHRD *hrd, Bool cprms_present_flag, UInt tempLevelHigh);
  Void  parseSliceHeader    ( TComSlice* pcSlice, ParameterSetManager *parameterSetManager, const Int prevTid0POC);
  Void  parseTerminatingBit ( UInt& ruiBit );
  Void  parseRemainingBytes ( Bool noTrailingBytesExpected );

  Void parseMVPIdx          ( Int& riMVPIdx );

  Void parseSkipFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if BARC_UP_SAMPLE_RDO
  Void parseUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx);
#endif
#if BARC
  Void parseBARCFlag        ( TComDataCU* pcCU, ChannelType eCtype );
#endif
#if COM16_C806_OBMC
  Void  parseOBMCFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if VCEG_AZ07_IMV
  Void  parseiMVFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if VCEG_AZ06_IC
  Void  parseICFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void parseCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void  parseCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if VCEG_AZ05_INTRA_MPI
  Void parseMPIIdx          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void parsePDPCIdx         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if MLIP_EXTEND_REF_LINE
  Void parseMLIPIdx           (TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if DIMD_INTRA_PRED
  Void parseDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
  Void parseDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Void parseROTIdx       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void parseROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#endif
  Void parseMergeFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void parseMergeIndex      ( TComDataCU* pcCU, UInt& ruiMergeIndex );
#if VCEG_AZ07_FRUC_MERGE
  Void parseFRUCMgrMode     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
#endif
  Void parseSplitFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void parseBTSplitMode     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
    , UInt uiBTDepth
#endif
    );
#else
  Void parsePartSize        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void parsePredMode        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void parseIntraDirLumaAng ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseIntraDirChroma  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void parseInterDir        ( TComDataCU* pcCU, UInt& ruiInterDir, UInt uiAbsPartIdx );
  Void parseRefFrmIdx       ( TComDataCU* pcCU, Int& riRefFrmIdx, RefPicList eRefList );
  Void parseMvd             ( TComDataCU* pcCU, UInt uiAbsPartAddr,UInt uiPartIdx,    UInt uiDepth, RefPicList eRefList );
  Void parseCrossComponentPrediction( class TComTU &rTu, ComponentID compID );

  Void parseDeltaQP         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseChromaQpAdjustment( TComDataCU* cu, UInt absPartIdx, UInt depth);

  Void parseCoeffNxN        ( class TComTU &rTu, ComponentID compID 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );

  Void parseTransformSkipFlags ( class TComTU &rTu, ComponentID component );

#if VCEG_AZ08_KLT_COMMON
  Void parseKLTFlags        ( TComTU &rTu, ComponentID component);
#endif

  Void parseIPCMInfo        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);

  Void xParsePredWeightTable ( TComSlice* pcSlice, const TComSPS *sps );
  Void  parseScalingList     ( TComScalingList* scalingList );
  Void xDecodeScalingList    ( TComScalingList *scalingList, UInt sizeId, UInt listId);

  Void  parseExplicitRdpcmMode( TComTU &rTu, ComponentID compID );
#if VCEG_AZ07_BAC_ADAPT_WDOW
  Void parseCtxUpdateInfo   ( TComSlice*& rpcSlice, TComStats* apcStats );
  Void xRunDecoding         ( Bool * uiCtxMAP, UInt uiNumCtx );
  Void xLevelDecoding       ( Bool * uiCtxMAP, UChar *uiCtxCodeIdx, UInt uiNumCtx );
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Void setStatesHandle      ( TComStats* pcStats) { m_pcStats = pcStats; }
  TComStats* getStatesHandle()                    { return m_pcStats; }
#endif

#if ALF_HM3_REFACTOR
  Void  xReadUnaryMaxSymbol ( UInt& ruiSymbol, UInt uiMaxSymbol );
  Void  parseAlfFlag        ( UInt& ruiVal );
#if JVET_C0038_GALF  
  Void parseALFTruncBinVal  (UInt& ruiSymbol, UInt uiMaxSymbol);
  Void parseALFPrevFiltType ( UInt& uiCode );
  Void parseALFPrevFiltFlag ( UInt& uiCode );
#endif
  Void  parseAlfUvlc        ( UInt& ruiVal );
  Void  parseAlfSvlc        ( Int&  riVal  );
  Void  parseAlfCtrlDepth   ( UInt& ruiAlfCtrlDepth , UInt uiMaxTotalCUDepth ); 
#if !JVET_C0024_QTBT
  Void  parseAlfCtrlFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth , UInt uiMaxAlfCtrlDepth );
#endif
  Void  parseAlfFlagNum     ( UInt& ruiVal, UInt minValue, UInt depth );
  Void  parseAlfCtrlFlag    ( UInt &ruiAlfCtrlFlag );
#endif
#if COM16_C806_EMT
  Void parseEmtTuIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void parseEmtCuFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf );
#endif
#if RRC_NS
	UInt parseNsFlag();
	UInt parseAdaptiveNsFlag();
	UInt parseNsFWholeFrameFlag();
	Int  parseNsBlockSize();
#endif

#if COM16_C1016_AFFINE
  Void  parseAffineFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
  Void  parseAffineMvd      ( TComDataCU* pcCU, UInt uiAbsPartAddr, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  Void  parseAffineParamFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
#endif
#if MRG_CONSIDER_COMPLEX
  Void parseMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void parseMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt& ruiMergeComplexIndex );
#endif

#if RRC_SELECTABLE_IF
  Int parseAdaptiveIfIndex();
#endif

#if RRC_INTRA_SHARPENING
  UInt parseIntraSharpFlag();
#endif

protected:
  Bool  xMoreRbspData();
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  TComStats* m_pcStats;
#endif

};

//! \}

#endif // !defined(AFX_TDECCAVLC_H__9732DD64_59B0_4A41_B29E_1A5B18821EAD__INCLUDED_)
