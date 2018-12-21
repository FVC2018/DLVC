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

/** \file     TDecEntropy.h
    \brief    entropy decoder class (header)
*/

#ifndef __TDECENTROPY__
#define __TDECENTROPY__

#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComBitStream.h"
#include "TLibCommon/TComSlice.h"
#include "TLibCommon/TComPic.h"
#include "TLibCommon/TComSampleAdaptiveOffset.h"
#include "TLibCommon/TComRectangle.h"
#if ALF_HM3_REFACTOR
#include "TLibCommon//TComAdaptiveLoopFilter.h"
#endif

class TDecSbac;
class TDecCavlc;
class ParameterSetManagerDecoder;
class TComPrediction;

//! \ingroup TLibDecoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// entropy decoder pure class
class TDecEntropyIf
{
public:
  //  Virtual list for SBAC/CAVLC
  virtual Void  resetEntropy          ( TComSlice* pcSlice )     = 0;
  virtual Void  setBitstream          ( TComInputBitstream* p )  = 0;

  virtual Void  parseVPS                  ( TComVPS* pcVPS )     = 0;
  virtual Void  parseSPS                  ( TComSPS* pcSPS )     = 0;
  virtual Void  parsePPS                  ( TComPPS* pcPPS )     = 0;

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  virtual TComStats* getStatesHandle      ( )                    = 0;
  virtual Void setStatesHandle            ( TComStats* pcStats)  = 0;
#if VCEG_AZ07_BAC_ADAPT_WDOW 
  virtual Void parseCtxUpdateInfo         ( TComSlice*& rpcSlice, TComStats* apcStats )   = 0;
#endif
#endif

  virtual Void parseSliceHeader          ( TComSlice* pcSlice, ParameterSetManager *parameterSetManager, const Int prevTid0POC)       = 0;

  virtual Void parseTerminatingBit       ( UInt& ruilsLast )                                     = 0;
  virtual Void parseRemainingBytes( Bool noTrailingBytesExpected ) = 0;

  virtual Void parseMVPIdx        ( Int& riMVPIdx ) = 0;

public:
  virtual Void parseSkipFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#if BLOCK_BASED_BACKGROUND_REFERENCE
	virtual Void parseCUBgFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
#if BARC_UP_SAMPLE_RDO
	virtual Void parseUpCNNFlag    (TComDataCU* pcCU, UInt uiAbsPartIdx) = 0;
#endif
#if BARC
  virtual Void parseBARCFlag      ( TComDataCU* pcCU, ChannelType eCtype ) = 0;
#endif
#if COM16_C806_OBMC
  virtual Void parseOBMCFlag      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
#if VCEG_AZ07_IMV
  virtual Void parseiMVFlag       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
#if VCEG_AZ06_IC
  virtual Void parseICFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
  virtual Void parseCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
  virtual Void parseSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#if JVET_C0024_QTBT
  virtual Void parseBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
    , UInt uiBTDepth
#endif
    ) = 0;
#endif
  virtual Void parseMergeFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx ) = 0;
#if VCEG_AZ05_INTRA_MPI
  virtual Void parseMPIIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  virtual Void parsePDPCIdx       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth) = 0;
#endif

#if DIMD_INTRA_PRED
  virtual Void parseDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight) = 0;
  virtual Void parseDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight) = 0;
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  virtual Void parseROTIdx     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#if JVET_C0024_QTBT
  virtual Void parseROTIdxChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
#endif
  virtual Void parseMergeIndex    ( TComDataCU* pcCU, UInt& ruiMergeIndex ) = 0;
#if VCEG_AZ07_FRUC_MERGE
  virtual Void parseFRUCMgrMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx ) = 0;
#endif
#if !JVET_C0024_QTBT
  virtual Void parsePartSize      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
#endif
  virtual Void parsePredMode      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;

  virtual Void parseIntraDirLumaAng( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
  virtual Void parseIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;

  virtual Void parseInterDir      ( TComDataCU* pcCU, UInt& ruiInterDir, UInt uiAbsPartIdx ) = 0;
  virtual Void parseRefFrmIdx     ( TComDataCU* pcCU, Int& riRefFrmIdx, RefPicList eRefList ) = 0;
  virtual Void parseMvd           ( TComDataCU* pcCU, UInt uiAbsPartAddr, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList ) = 0;
  virtual Void parseCrossComponentPrediction ( class TComTU &rTu, ComponentID compID ) = 0;

  virtual Void parseTransformSubdivFlag( UInt& ruiSubdivFlag, UInt uiLog2TransformBlockSize ) = 0;
  virtual Void parseQtCbf         ( TComTU &rTu, const ComponentID compID, const Bool lowestLevel ) = 0;
  virtual Void parseQtRootCbf     ( UInt uiAbsPartIdx, UInt& uiQtRootCbf ) = 0;

  virtual Void parseDeltaQP       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
  virtual Void parseChromaQpAdjustment( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;

  virtual Void parseIPCMInfo     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth) = 0;

  virtual Void parseCoeffNxN( class TComTU &rTu, ComponentID compID  
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    ) = 0;

  virtual Void parseTransformSkipFlags ( class TComTU &rTu, ComponentID component ) = 0;

#if VCEG_AZ08_KLT_COMMON
  virtual Void parseKLTFlags(TComTU &rTu, ComponentID component) = 0;
#endif

  virtual Void parseExplicitRdpcmMode ( TComTU &rTu, ComponentID compID ) = 0;

  virtual ~TDecEntropyIf() {}

#if ALF_HM3_REFACTOR
  virtual Void parseAlfFlag       ( UInt& ruiVal           ) = 0;
#if JVET_C0038_GALF  
  virtual Void parseALFTruncBinVal  ( UInt& ruiSymbol, UInt uiMaxSymbol ) = 0;
  virtual Void parseALFPrevFiltType ( UInt& uiCode ) = 0;
  virtual Void parseALFPrevFiltFlag ( UInt& uiCode ) = 0;
#endif
  virtual Void parseAlfUvlc       ( UInt& ruiVal           ) = 0;
  virtual Void parseAlfSvlc       ( Int&  riVal            ) = 0;
  virtual Void parseAlfCtrlDepth   ( UInt& ruiAlfCtrlDepth , UInt uiMaxTotalCUDepth ) = 0; 
#if !JVET_C0024_QTBT
  virtual Void parseAlfCtrlFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth , UInt uiMaxAlfCtrlDepth ) = 0;
#endif
  virtual Void parseAlfFlagNum    ( UInt& ruiVal, UInt minValue, UInt depth ) = 0;
  virtual Void parseAlfCtrlFlag   ( UInt &ruiAlfCtrlFlag ) = 0;
#endif

#if RRC_NS
	virtual UInt parseNsFlag() = 0;
	virtual UInt parseAdaptiveNsFlag() = 0;
	virtual UInt parseNsFWholeFrameFlag() = 0;
	virtual Int  parseNsBlockSize() = 0;
#endif

#if COM16_C806_EMT
  virtual Void parseEmtTuIdx      ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth ) = 0;
  virtual Void parseEmtCuFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRootCbf ) = 0;
#endif

#if COM16_C1016_AFFINE
  virtual Void parseAffineFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx ) = 0;
  virtual Void parseAffineMvd     ( TComDataCU* pcCU, UInt uiAbsPartAddr, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList ) = 0;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  virtual Void parseAffineParamFlag ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx ) = 0;
#endif

#if MRG_CONSIDER_COMPLEX
  virtual Void parseMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx ) = 0;
  virtual Void parseMergeComplexIndex ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt& ruiMergeComplexIndex ) = 0;
#endif

#if RRC_INTRA_SHARPENING
  virtual UInt parseIntraSharpFlag() = 0;
#endif

#if RRC_SELECTABLE_IF
  virtual Int parseAdaptiveIfIndex() = 0;
#endif
};

/// entropy decoder class
class TDecEntropy
{
private:
  TDecEntropyIf*  m_pcEntropyDecoderIf;
  TComPrediction* m_pcPrediction;
  //UInt    m_uiBakAbsPartIdx;
  //UInt    m_uiBakChromaOffset;
  //UInt    m_bakAbsPartIdxCU;

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComMvField*        m_pMvFieldSP[NUM_MGR_TYPE];
  UChar*              m_phInterDirSP[NUM_MGR_TYPE];
#else
  TComMvField*        m_pMvFieldSP[2];
  UChar*              m_phInterDirSP[2];
#endif
#endif

public:

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  TDecEntropy();
  ~TDecEntropy();
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  TComStats* m_pcStats;
  Void       updateStates   ( SliceType uiSliceType, UInt uiSliceQP, TComStats*  apcStats);
  Void       setStatsHandle ( TComStats*  pcStats)  { m_pcStats = pcStats; }
  TComStats* getStatsHandle ( )                     { return m_pcStats; }
#endif 

  Void init (TComPrediction* p) {m_pcPrediction = p;}
  Void decodePUWise       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, TComDataCU* pcSubCU );
  Void decodeInterDirPU   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx );
  Void decodeRefFrmIdxPU  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList );
  Void decodeMvdPU        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList );
  Void decodeMVPIdxPU     ( TComDataCU* pcSubCU, UInt uiPartAddr, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList );

  Void    setEntropyDecoder           ( TDecEntropyIf* p );
  Void    setBitstream                ( TComInputBitstream* p ) { m_pcEntropyDecoderIf->setBitstream(p);                    }
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Void    resetEntropy                ( TComSlice* p )       
  {
    m_pcEntropyDecoderIf->setStatesHandle (m_pcStats);
    m_pcEntropyDecoderIf->resetEntropy(p);  
  }
#if VCEG_AZ07_BAC_ADAPT_WDOW
  Void    decodeCtxUpdateInfo         ( TComSlice*& rpcSlice, TComStats* apcStats )  { m_pcEntropyDecoderIf->parseCtxUpdateInfo( rpcSlice, apcStats); }
#endif
#else
  Void    resetEntropy                ( TComSlice* p)           { m_pcEntropyDecoderIf->resetEntropy(p);                    }
#endif

  Void    decodeVPS                   ( TComVPS* pcVPS ) { m_pcEntropyDecoderIf->parseVPS(pcVPS); }
  Void    decodeSPS                   ( TComSPS* pcSPS ) { m_pcEntropyDecoderIf->parseSPS(pcSPS); }
  Void    decodePPS                   ( TComPPS* pcPPS ) { m_pcEntropyDecoderIf->parsePPS(pcPPS); }
  Void    decodeSliceHeader           ( TComSlice* pcSlice, ParameterSetManager *parameterSetManager, const Int prevTid0POC)  { m_pcEntropyDecoderIf->parseSliceHeader(pcSlice, parameterSetManager, prevTid0POC);         }

  Void    decodeTerminatingBit        ( UInt& ruiIsLast )       { m_pcEntropyDecoderIf->parseTerminatingBit(ruiIsLast);     }
  Void    decodeRemainingBytes( Bool noTrailingBytesExpected ) { m_pcEntropyDecoderIf->parseRemainingBytes(noTrailingBytesExpected); }

  TDecEntropyIf* getEntropyDecoder() { return m_pcEntropyDecoderIf; }

public:
#if JVET_C0024_QTBT
  Void decodeBTSplitMode       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
    , UInt uiBTDepth
#endif
    );
#endif
  Void decodeSplitFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void decodeSkipFlag          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void decodeCUBgFlag          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if BARC_UP_SAMPLE_RDO
	Void decodeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif
#if BARC
  Void decodeBARCFlag          ( TComDataCU* pcCU, ChannelType eCtype );
#endif
#if VCEG_AZ05_INTRA_MPI
  Void decodeMPIIdx            ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  Void decodePDPCIdx           ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif

#if DIMD_INTRA_PRED
  Void decodeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
  Void decodeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Void decodeROTIdx        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if JVET_C0024_QTBT
  Void decodeROTIdxChroma  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#endif
#if COM16_C806_OBMC
  Void decodeOBMCFlag          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if VCEG_AZ07_IMV
  Void decodeiMVFlag           ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
#if VCEG_AZ06_IC
  Void decodeICFlag            ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif
  Void decodeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void decodeMergeFlag         ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void decodeMergeIndex        ( TComDataCU* pcSubCU, UInt uiPartIdx, UInt uiPartAddr, UInt uiDepth );
#if VCEG_AZ07_FRUC_MERGE
  Void decodeFRUCMgrMode       ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
#endif
  Void decodePredMode          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#if !JVET_C0024_QTBT
  Void decodePartSize          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif

  Void decodeIPCMInfo          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void decodePredInfo          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, TComDataCU* pcSubCU );

  Void decodeIntraDirModeLuma  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
  Void decodeIntraDirModeChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );

  Void decodeQP                ( TComDataCU* pcCU, UInt uiAbsPartIdx );
  Void decodeChromaQpAdjustment( TComDataCU* pcCU, UInt uiAbsPartIdx );

private:
#if JVET_C0024_QTBT
  Void xDecodeTransform        ( Bool& bCodeDQP, Bool& isChromaQpAdjCoded, TComTU &rTu, ComponentID compID);
#else
  Void xDecodeTransform        ( Bool& bCodeDQP, Bool& isChromaQpAdjCoded, TComTU &rTu, const Int quadtreeTULog2MinSizeInCU 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    );
#endif

public:

  Void decodeCoeff             ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool& bCodeDQP, Bool& isChromaQpAdjCoded );

#if QTBT_NSST
  Int countNonZeroCoeffs       ( TCoeff* pcCoef, UInt uiSize );
#endif

#if ALF_HM3_REFACTOR
  // Adaptive Loop filter
  Void decodeAlfParam(ALFParam* pAlfParam, UInt uiMaxTotalCUDepth
#if FIX_TICKET12
        , TComSlice * pSlice
#endif
    );
  Void decodeAux(ALFParam* pAlfParam);
  Void decodeFilt(ALFParam* pAlfParam);
#if JVET_C0038_GALF
  Void readFilterCodingParams(ALFParam* pAlfParam, Bool bChroma = false);
  Void readFilterCoeffs(ALFParam* pAlfParam, Bool bChroma = false);
  Void decodeFilterCoeff (ALFParam* pAlfParam, Bool bChroma = false);
#else
  Void readFilterCodingParams(ALFParam* pAlfParam);
  Void readFilterCoeffs(ALFParam* pAlfParam);
  Void decodeFilterCoeff (ALFParam* pAlfParam);
#endif
  Int  golombDecode(Int k);
  Void decodeAlfCtrlParam      ( ALFParam *pAlfParam );
#endif
#if RRC_NS
	UInt decodeNsFlag();
	UInt decodeAdaptiveNsFlag();
	UInt decodeNsFWholeFrameFlag();
	Int  decodeNsBlockSize();
#endif
#if COM16_C1016_AFFINE
  Void decodeAffineFlag        ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  Void decodeAffineParamFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx );
#endif

#if MRG_CONSIDER_COMPLEX
  Void decodeMergeComplexFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx );
  Void decodeMergeComplexIndex ( TComDataCU* pcSubCU, UInt uiPartIdx, UInt uiPartAddr, UInt uiDepth );
#endif

#if RRC_INTRA_SHARPENING
  Void decodeIntraSharpFlag    ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth );
#endif

#if RRC_SELECTABLE_IF
  Void decodeAdaptiveIfIndex(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth);
#endif
};// END CLASS DEFINITION TDecEntropy

//! \}

#endif // __TDECENTROPY__
