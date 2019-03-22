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

/** \file     TEncCu.h
    \brief    Coding Unit (CU) encoder class (header)
*/

#ifndef __TENCCU__
#define __TENCCU__

// Include files
#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComYuv.h"
#include "TLibCommon/TComPrediction.h"
#include "TLibCommon/TComTrQuant.h"
#include "TLibCommon/TComBitCounter.h"
#include "TLibCommon/TComDataCU.h"

#include "TEncEntropy.h"
#include "TEncSearch.h"
#include "TEncRateCtrl.h"
#if VCEG_AZ07_IMV
#include <set>
#endif
//! \ingroup TLibEncoder
//! \{

class TEncTop;
class TEncSbac;
class TEncCavlc;
class TEncSlice;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CU encoder class
class TEncCu
{
private:

#if JVET_C0024_QTBT
  TComDataCU***            m_pppcBestCU;      ///< Best CUs in each depth
  TComDataCU***            m_pppcTempCU;      ///< Temporary CUs in each depth
#if COM16_C806_OBMC
  TComDataCU***            m_pppcTempCUWoOBMC; ///< Temporary CUs in each depth
#endif
#if VCEG_AZ07_FRUC_MERGE
  TComDataCU***            m_pppcFRUCBufferCU;      
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
  TComDataCU***            m_pppcTempCUIMVCache[NUMBER_OF_PART_SIZES]; 
#endif
#else
  TComDataCU**            m_ppcBestCU;      ///< Best CUs in each depth
  TComDataCU**            m_ppcTempCU;      ///< Temporary CUs in each depth
#if COM16_C806_OBMC
  TComDataCU**            m_ppcTempCUWoOBMC; ///< Temporary CUs in each depth
#endif
#if VCEG_AZ07_FRUC_MERGE
  TComDataCU**            m_ppcFRUCBufferCU;      
#endif
#if VCEG_AZ07_IMV
  TComDataCU**            m_ppcTempCUIMVCache[NUMBER_OF_PART_SIZES]; 
#endif
#endif
  UChar                   m_uhTotalDepth;

#if JVET_C0024_QTBT
#if JVET_C0024_FAST_MRG
  TComYuv*                 m_pcMrgPredTempYuv[MRG_MAX_NUM_CANDS]; ///< Best Prediction Yuv for each depth
#endif
  TComYuv***               m_pppcPredYuvBest; ///< Best Prediction Yuv for each depth
  TComYuv***               m_pppcResiYuvBest; ///< Best Residual Yuv for each depth
  TComYuv***               m_pppcRecoYuvBest; ///< Best Reconstruction Yuv for each depth
  TComYuv***               m_pppcPredYuvTemp; ///< Temporary Prediction Yuv for each depth
  TComYuv***               m_pppcResiYuvTemp; ///< Temporary Residual Yuv for each depth
  TComYuv***               m_pppcRecoYuvTemp; ///< Temporary Reconstruction Yuv for each depth
  TComYuv***               m_pppcOrigYuv;     ///< Original Yuv for each depth
#if COM16_C806_OBMC
  TComYuv***               m_pppcTmpYuv1;     ///< Temporary Yuv used for OBMC
  TComYuv***               m_pppcTmpYuv2;     ///< Temporary Yuv used for OBMC, seems no use? JCA
  TComYuv***               m_pppcPredYuvWoOBMC; ///< Temporary Prediction Yuv for each depth
#endif
#else
  TComYuv**               m_ppcPredYuvBest; ///< Best Prediction Yuv for each depth
  TComYuv**               m_ppcResiYuvBest; ///< Best Residual Yuv for each depth
  TComYuv**               m_ppcRecoYuvBest; ///< Best Reconstruction Yuv for each depth
  TComYuv**               m_ppcPredYuvTemp; ///< Temporary Prediction Yuv for each depth
  TComYuv**               m_ppcResiYuvTemp; ///< Temporary Residual Yuv for each depth
  TComYuv**               m_ppcRecoYuvTemp; ///< Temporary Reconstruction Yuv for each depth
  TComYuv**               m_ppcOrigYuv;     ///< Original Yuv for each depth
#if COM16_C806_OBMC
  TComYuv**               m_ppcTmpYuv1;     ///< Temporary Yuv used for OBMC
  TComYuv**               m_ppcTmpYuv2;     ///< Temporary Yuv used for OBMC
  TComYuv**               m_ppcPredYuvWoOBMC; ///< Temporary Prediction Yuv for each depth
#endif
#endif
#if COM16_C806_LARGE_CTU
  Pel*                    m_resiBuffer[NUMBER_OF_STORED_RESIDUAL_TYPES];
#endif
#if BARC_FUNC
#if BARC_ONLY_FOR_LUMA
  Pel*                    m_reconLumaBufferHM;
#else
  Pel*                    m_reconBuffer[MAX_NUM_COMPONENT];
  Pel*                    m_reconLumaBuffer;
#endif
#if COMPACT_RESOLUTION_CNN
  Pel*                    m_reconBufferBicubic;
  Pel*                    m_DScaleBufferBicubic;
#endif
#endif
#if JVET_E0076_MULTI_PEL_MVD
  Double                 m_dBestMvDPelCost[3]; //0: 1/4, 1: Int, 2: 4-pel
#endif
  //  Data : encoder control
  Bool                    m_bEncodeDQP;
  Bool                    m_bFastDeltaQP;
  Bool                    m_stillToCodeChromaQpOffsetFlag; //indicates whether chroma QP offset flag needs to coded at this particular CU granularity.
  Int                     m_cuChromaQpOffsetIdxPlus1; // if 0, then cu_chroma_qp_offset_flag will be 0, otherwise cu_chroma_qp_offset_flag will be 1.
#if WCG_LUMA_DQP_CM_SCALE  
  Int                     m_lumaLevelToDeltaQPLUT[LUMA_LEVEL_TO_DQP_LUT_MAXSIZE];
  Int                     m_LumaQPOffset;
#endif
  //  Access channel
  TEncCfg*                m_pcEncCfg;
  TEncSearch*             m_pcPredSearch;
  TComTrQuant*            m_pcTrQuant;
  TComRdCost*             m_pcRdCost;

  TEncEntropy*            m_pcEntropyCoder;
  TEncBinCABAC*           m_pcBinCABAC;

  // SBAC RD
#if JVET_C0024_QTBT
  TEncSbac****            m_ppppcRDSbacCoder;
#else
  TEncSbac***             m_pppcRDSbacCoder;
#endif
  TEncSbac*               m_pcRDGoOnSbacCoder;
  TEncRateCtrl*           m_pcRateCtrl;

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  //ATMVP 
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComMvField  * m_pMvFieldSP[NUM_MGR_TYPE];
  UChar        * m_phInterDirSP[NUM_MGR_TYPE];
#else
  TComMvField  * m_pMvFieldSP[2];
  UChar        * m_phInterDirSP[2];
#endif
#endif

public:
  /// copy parameters from encoder class
  Void  init                ( TEncTop* pcEncTop );

#if WCG_LUMA_DQP_CM_SCALE
  TEncSlice*             m_pcSliceEncoder;
  Void                   setSliceEncoder( TEncSlice* pSliceEncoder) {m_pcSliceEncoder=pSliceEncoder;}
  TEncSlice*             getSliceEncoder() { return m_pcSliceEncoder;}
  Void                   initLumaDeltaQpLUT(); 
  Int                    calculateLumaDQP(TComDataCU *pcCU, const UInt uiAbsPartIdx, const TComYuv * pOrgYuv);
#endif
  /// create internal buffers
  Void  create              ( UChar uhTotalDepth, UInt iMaxWidth, UInt iMaxHeight, ChromaFormat chromaFormat );

  /// destroy internal buffers
  Void  destroy             ();

  /// CTU analysis function
  Void  compressCtu         ( TComDataCU*  pCtu );

  /// CTU encoding function
  Void  encodeCtu           ( TComDataCU*  pCtu );

  Int   updateCtuDataISlice ( TComDataCU* pCtu, Int width, Int height );

  Void setFastDeltaQp       ( Bool b)                 { m_bFastDeltaQP = b;         }

protected:
  Void  finishCU            ( TComDataCU*  pcCU, UInt uiAbsPartIdx );

#if AMP_ENC_SPEEDUP
#if JVET_C0024_QTBT
  Void  xCompressCU         ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug), UInt uiSplitConstrain=0 );
  Void  CompressBT          ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug), Int* piParam, Double* pdParam, UInt uiSplitConstrain=0 );
#if TT_ENC_OPT
  Void  CompressTT          ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug), Int* piParam, Double* pdParam, UInt uiSplitConstrain=0 );
#endif
#else
  Void  xCompressCU         ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth DEBUG_STRING_FN_DECLARE(sDebug), PartSize eParentPartSize = NUMBER_OF_PART_SIZES );
#endif
#else
  Void  xCompressCU         ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth        );
#endif
#if JVET_C0024_QTBT
  Void  xEncodeCU           ( TComDataCU*  pcCU, UInt uiAbsPartIdx,           UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiSplitConstrain=0);
#else
  Void  xEncodeCU           ( TComDataCU*  pcCU, UInt uiAbsPartIdx,           UInt uiDepth        );
#endif

  Int   xComputeQP          ( TComDataCU* pcCU, UInt uiDepth );
#if JVET_C0024_QTBT
  Void  xCheckBestMode      ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth, UInt uiWidth=0, UInt uiHeight=0 DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo=true));
#else
  Void  xCheckBestMode      ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo=true));
#endif

  Void  xCheckRDCostMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU DEBUG_STRING_FN_DECLARE(sDebug), Bool *earlyDetectionSkipMode );
#if VCEG_AZ07_FRUC_MERGE
  Void  xCheckRDCostMerge2Nx2NFRUC( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU , Bool *earlyDetectionSkipMode );
#endif
#if COM16_C1016_AFFINE
  Void  xCheckRDCostAffineMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU );
#endif

#if AMP_MRG
  Void  xCheckRDCostInter   ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseMRG = false 
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
    , UChar bIMV = 0 , TComDataCU * pcCUInfo2Reuse = NULL 
#else
    , Bool bIMV = false , TComDataCU * pcCUInfo2Reuse = NULL 
#endif
#endif
    );
#else
  Void  xCheckRDCostInter   ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize 
#if VCEG_AZ07_IMV
    , Bool bIMV = false , TComDataCU * pcCUInfo2Reuse = NULL 
#endif
    );
#endif

#if VCEG_AZ08_INTER_KLT
  Void  xCheckRDCostInterKLT(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize);
#endif

  Void  xCheckRDCostIntra   ( TComDataCU *&rpcBestCU,
                              TComDataCU *&rpcTempCU,
                              Double      &cost,
                              PartSize     ePartSize
                              DEBUG_STRING_FN_DECLARE(sDebug)
#if VCEG_AZ05_ROT_TR   || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                              , Int& bNonZeroCoeff
#endif
#if DIMD_INTRA_PRED
    , UChar ucDIMDIntraRD
#endif
// #if DIMD_MLIP_HARM
//       , Bool bUseDIMD
// #endif
                            );

  Void  xCheckDQP           ( TComDataCU*  pcCU );

  Void  xCheckIntraPCM      ( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU                      );
  Void  xCopyAMVPInfo       ( AMVPInfo* pSrc, AMVPInfo* pDst );
#if JVET_C0024_QTBT //uiSplitMethod: 0: quadtree; 1: hor; 2: ver
  Void  xCopyYuv2Pic        (TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth, UInt uiWidth, UInt uiHeight 
#if BARC
    ,Bool isBARC = false
#endif
  );
#if TT_TREE_DS1
  Void  xCopyYuv2Tmp        ( UInt uiPartUnitIdxSrc, UInt uiPartUnitIdxDst, UInt uiWidth, UInt uiHeight, UInt uiSplitMethod=0 );
#else
  Void  xCopyYuv2Tmp        ( UInt uiPartUnitIdx, UInt uiWidth, UInt uiHeight, UInt uiSplitMethod=0 );
#endif
#else
  Void  xCopyYuv2Pic        (TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth );
  Void  xCopyYuv2Tmp        ( UInt uhPartUnitIdx, UInt uiDepth );
#endif

  Bool getdQPFlag           ()                        { return m_bEncodeDQP;        }
  Void setdQPFlag           ( Bool b )                { m_bEncodeDQP = b;           }

  Bool getFastDeltaQp       () const                  { return m_bFastDeltaQP;      }

  Bool getCodeChromaQpAdjFlag() { return m_stillToCodeChromaQpOffsetFlag; }
  Void setCodeChromaQpAdjFlag( Bool b ) { m_stillToCodeChromaQpOffsetFlag = b; }

#if ADAPTIVE_QP_SELECTION
  // Adaptive reconstruction level (ARL) statistics collection functions
  Void xCtuCollectARLStats(TComDataCU* pCtu);
  Int  xTuCollectARLStats(TCoeff* rpcCoeff, TCoeff* rpcArlCoeff, Int NumCoeffInCU, Double* cSum, UInt* numSamples );
#endif

#if !JVET_C0024_QTBT
#if AMP_ENC_SPEEDUP
#if AMP_MRG
  Void deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver, Bool &bTestMergeAMP_Hor, Bool &bTestMergeAMP_Ver
    );
#else
  Void deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver
    );
#endif
#endif
#endif

  Void  xFillPCMBuffer     ( TComDataCU* pCU, TComYuv* pOrgYuv );
#if BARC_FUNC
  Void reComputeBARCDistortion(TComDataCU* pCtu, ChannelType chType);
  Void upSampleRec( TComDataCU*& pCtu, ChannelType chType);
  Void downSampleRec( TComDataCU*& pCtu, ChannelType chType );
  Void initEncSliceBARC( ChannelType chType, TComPic* pcPic, Int iGOPid, TComSlice*& rpcSlice, Bool isBARC );
  Void copyLumaDscaletoRec(TComDataCU* pCtu, Bool isCopyDscaletoRec);
  Void copyRecToTmp(TComDataCU* pCtu, Bool isCopyRecToTmp);
  UInt getDist(Pel* piCur, Int iCurStride, Pel* piOrg, Int iOrgStride, UInt uiBlkWidth, UInt uiBlkHeight, UInt bitDepth = 10);
#if COMPACT_RESOLUTION_CNN
  Void copyBicubicRecToTmp(TComDataCU* pCtu, Bool isCopyRecToTmp);
  Void upSampleRecDeep(TComDataCU*& pCtu, ChannelType chType);
#endif
#endif

#if VCEG_AZ07_IMV
private:
  typedef struct
  {
    Double    dRDCost;
    Bool      bUseMrg;
    PartSize  eInterPartSize;
    TComDataCU * pcCUMode;
  }SModeCand;

  class cmpModeCand
  {
  public:
    bool operator()(const SModeCand & r1 , const SModeCand & r2 ) const
    {
      return r1.dRDCost < r2.dRDCost;
    }
  };
  std::set <SModeCand, cmpModeCand> m_setInterCand;
#endif
};

//! \}

#endif // __TENCMB__
