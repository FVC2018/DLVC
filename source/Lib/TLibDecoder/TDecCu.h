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

/** \file     TDecCu.h
    \brief    CU decoder class (header)
*/

#ifndef __TDECCU__
#define __TDECCU__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TLibCommon/TComTrQuant.h"
#include "TDecEntropy.h"
#if VCEG_AZ08_INTER_KLT
#include "TLibCommon/TComPrediction.h"
#endif
#if JVET_F0096_BILATERAL_FILTER
#include "TLibCommon/TComBilateralFilter.h"
#endif
//! \ingroup TLibDecoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CU decoder class
class TDecCu
{
private:
  UInt                m_uiMaxDepth;       ///< max. number of depth
#if JVET_C0024_QTBT
  TComYuv***          m_pppcYuvResi;
  TComYuv***          m_pppcYuvReco;
  TComDataCU***       m_pppcCU;
#if COM16_C806_OBMC
  TComYuv***           m_pppcTmpYuv1;       ///< array of OBMC prediction buffer
  TComYuv***           m_pppcTmpYuv2;
#endif
#else
  TComYuv**           m_ppcYuvResi;       ///< array of residual buffer
  TComYuv**           m_ppcYuvReco;       ///< array of prediction & reconstruction buffer
  TComDataCU**        m_ppcCU;            ///< CU data array
#if COM16_C806_OBMC
  TComYuv**           m_ppcTmpYuv1;       ///< array of OBMC prediction buffer
  TComYuv**           m_ppcTmpYuv2;
#endif
#endif
  // access channel
  TComTrQuant*        m_pcTrQuant;
  TComPrediction*     m_pcPrediction;
  TDecEntropy*        m_pcEntropyDecoder;

  Bool                m_bDecodeDQP;
  Bool                m_IsChromaQpAdjCoded;

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION 
  TComMvField*        m_pMvFieldSP[NUM_MGR_TYPE];
  UChar*              m_phInterDirSP[NUM_MGR_TYPE];
#else
  TComMvField*        m_pMvFieldSP[2];
  UChar*              m_phInterDirSP[2];
#endif
#endif
#if BARC_FUNC
  Pel*               m_reconLumaBuffer;
#endif
public:
  TDecCu();
  virtual ~TDecCu();

  /// initialize access channels
  Void  init                    ( TDecEntropy* pcEntropyDecoder, TComTrQuant* pcTrQuant, TComPrediction* pcPrediction );

  /// create internal buffers
  Void  create                  ( UInt uiMaxDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormatIDC );

  /// destroy internal buffers
  Void  destroy                 ();

  /// decode Ctu information
  Void  decodeCtu               ( TComDataCU* pCtu, Bool &isLastCtuOfSliceSegment );

  /// reconstruct Ctu information
  Void  decompressCtu           ( TComDataCU* pCtu );

#if VCEG_AZ08_INTER_KLT
  TComPrediction* getPointerPrediction() { return m_pcPrediction; }
#endif

protected:

#if JVET_C0024_QTBT
  Void xDecodeCU                ( TComDataCU* const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, const UInt uiWidth, const UInt uiHeight
#if TT_TREE_DS1
    , UInt uiBTDepth
#endif
    , Bool &isLastCtuOfSliceSegment, UInt uiSplitConstrain=0);
#else
  Void xDecodeCU                ( TComDataCU* const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, Bool &isLastCtuOfSliceSegment);
#endif
  Void xFinishDecodeCU          ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool &isLastCtuOfSliceSegment);
  Bool xDecodeSliceEnd          ( TComDataCU* pcCU, UInt uiAbsPartIdx );
#if JVET_C0024_QTBT
  Void xDecompressCU            ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight );
#else
  Void xDecompressCU            ( TComDataCU* pCtu, UInt uiAbsPartIdx, UInt uiDepth );
#endif

  Void xReconInter              ( TComDataCU* pcCU, UInt uiDepth );

  Void xReconIntraQT            ( TComDataCU* pcCU, UInt uiDepth );
  Void xIntraRecBlk             ( TComYuv* pcRecoYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, const ComponentID component, TComTU &rTu );
#if VCEG_AZ08_INTRA_KLT
  Void xIntraRecBlkTM           ( TComYuv* pcRecoYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, const ComponentID compID, TComTU &rTu, Int tmpred0_tmpredklt1_ori2 );
#endif
  Void xIntraRecQT              ( TComYuv* pcRecoYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, const ChannelType chType, TComTU &rTu );

  Void xReconPCM                ( TComDataCU* pcCU, UInt uiDepth );

  Void xDecodeInterTexture      ( TComDataCU* pcCU, UInt uiDepth );
  Void xDecodePCMTexture        ( TComDataCU* pcCU, const UInt uiPartIdx, const Pel *piPCM, Pel* piReco, const UInt uiStride, const UInt uiWidth, const UInt uiHeight, const ComponentID compID);

#if !JVET_C0024_QTBT
  Void xCopyToPic               ( TComDataCU* pcCU, TComPic* pcPic, UInt uiZorderIdx, UInt uiDepth );
#endif

  Bool getdQPFlag               ()                        { return m_bDecodeDQP;        }
  Void setdQPFlag               ( Bool b )                { m_bDecodeDQP = b;           }
  Bool getIsChromaQpAdjCoded    ()                        { return m_IsChromaQpAdjCoded;        }
  Void setIsChromaQpAdjCoded    ( Bool b )                { m_IsChromaQpAdjCoded = b;           }

  Void xFillPCMBuffer           (TComDataCU* pCU, UInt depth);

#if VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
  Void xDeriveCUMV              ( TComDataCU * pcCU , UInt uiAbsPartIdx , UInt uiDepth );
#endif
#if BARC
  Void upSampleRec( TComDataCU*& pCtu, ChannelType chType);
  Void downSampleRec( TComDataCU*& pCtu, ChannelType chType);
  Void copyLumaDscaletoRec(TComDataCU* pCtu, Bool isCopyDscaletoRec);
#endif
};

//! \}

#endif

