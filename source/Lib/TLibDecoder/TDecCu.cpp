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

/** \file     TDecCu.cpp
    \brief    CU decoder class
*/

#include "TDecCu.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/TComPrediction.h"


//! \ingroup TLibDecoder
//! \{
#if VCEG_AZ08_INTER_KLT
extern UInt g_uiDepth2InterTempSize[5];
#endif
#if VCEG_AZ08_INTRA_KLT
extern UInt g_uiDepth2IntraTempSize[5];
#endif
// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TDecCu::TDecCu()
{
#if JVET_C0024_QTBT
  m_pppcYuvResi = NULL;
  m_pppcYuvReco = NULL;
  m_pppcCU      = NULL;
#else
  m_ppcYuvResi = NULL;
  m_ppcYuvReco = NULL;
  m_ppcCU      = NULL;
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
  for (Int i=0; i< NUM_MGR_TYPE; i++)
  {
    m_pMvFieldSP[i] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
    m_phInterDirSP[i] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
    assert( m_pMvFieldSP[i] != NULL && m_phInterDirSP[i] != NULL );
  }
#else
  m_pMvFieldSP[0] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
  m_pMvFieldSP[1] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
  m_phInterDirSP[0] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
  m_phInterDirSP[1] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
  assert( m_pMvFieldSP[0] != NULL && m_phInterDirSP[0] != NULL );
  assert( m_pMvFieldSP[1] != NULL && m_phInterDirSP[1] != NULL );
#endif
#endif
#if COM16_C806_OBMC
#if JVET_C0024_QTBT
  m_pppcTmpYuv1 = NULL;
  m_pppcTmpYuv2 = NULL;
#else
  m_ppcTmpYuv1 = NULL;
  m_ppcTmpYuv2 = NULL;
#endif
#endif
#if BARC_FUNC
  m_reconLumaBuffer = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#endif
}

TDecCu::~TDecCu()
{
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION 
  for (UInt ui=0;ui<NUM_MGR_TYPE;ui++) 
#else 
  for (UInt ui=0;ui<2;ui++)
#endif 
  {
    if( m_pMvFieldSP[ui] != NULL )
    {
      delete [] m_pMvFieldSP[ui];
      m_pMvFieldSP[ui] = NULL;
    }
    if( m_phInterDirSP[ui] != NULL )
    {
      delete [] m_phInterDirSP[ui];
      m_phInterDirSP[ui] = NULL;
    }
  }
#endif
#if BARC_FUNC
  delete[] m_reconLumaBuffer;
  m_reconLumaBuffer = NULL;
#endif
}

Void TDecCu::init( TDecEntropy* pcEntropyDecoder, TComTrQuant* pcTrQuant, TComPrediction* pcPrediction)
{
  m_pcEntropyDecoder  = pcEntropyDecoder;
  m_pcTrQuant         = pcTrQuant;
  m_pcPrediction      = pcPrediction;
}

/**
 \param    uiMaxDepth      total number of allowable depth
 \param    uiMaxWidth      largest CU width
 \param    uiMaxHeight     largest CU height
 \param    chromaFormatIDC chroma format
 */
Void TDecCu::create( UInt uiMaxDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormatIDC )
{
  m_uiMaxDepth = uiMaxDepth+1;

#if JVET_C0024_QTBT
  assert((1<<MIN_CU_LOG2<<uiMaxDepth)==uiMaxWidth);

  UInt uiNumWidthIdx = g_aucConvertToBit[uiMaxWidth] + 1;
  UInt uiNumHeightIdx = g_aucConvertToBit[uiMaxHeight] + 1;

  m_pppcCU      = new TComDataCU**[uiNumWidthIdx];

  m_pppcYuvResi = new TComYuv**[uiNumWidthIdx];
  m_pppcYuvReco = new TComYuv**[uiNumWidthIdx];
#if COM16_C806_OBMC
  m_pppcTmpYuv1 = new TComYuv**[uiNumWidthIdx];
  m_pppcTmpYuv2 = new TComYuv**[uiNumWidthIdx];
#endif

  for (UInt wIdx=0; wIdx<uiNumWidthIdx; wIdx++)
  {
    m_pppcCU[wIdx] = new TComDataCU* [uiNumHeightIdx];
    m_pppcYuvResi[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcYuvReco[wIdx] = new TComYuv*[uiNumHeightIdx];
#if COM16_C806_OBMC
    m_pppcTmpYuv1[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcTmpYuv2[wIdx] = new TComYuv*[uiNumHeightIdx];
#endif
    for (UInt hIdx=0; hIdx<uiNumHeightIdx; hIdx++)
    {
      UInt uiWidth = 1<<(wIdx+MIN_CU_LOG2);
      UInt uiHeight = 1<<(hIdx+MIN_CU_LOG2);

      m_pppcCU[wIdx][hIdx]      = new TComDataCU; m_pppcCU[wIdx][hIdx]->create( chromaFormatIDC, 1<<( ( m_uiMaxDepth - 1 )<<1 ), uiMaxWidth, uiMaxHeight, true, uiMaxWidth >> (m_uiMaxDepth - 1), uiWidth, uiHeight );
      m_pppcYuvResi[wIdx][hIdx] = new TComYuv;    m_pppcYuvResi[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormatIDC );
      m_pppcYuvReco[wIdx][hIdx] = new TComYuv;    m_pppcYuvReco[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormatIDC );
#if COM16_C806_OBMC
      m_pppcTmpYuv1[wIdx][hIdx] = new TComYuv;    m_pppcTmpYuv1[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormatIDC );
      m_pppcTmpYuv2[wIdx][hIdx] = new TComYuv;    m_pppcTmpYuv2[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormatIDC );
#endif
    }
  }
#else
  m_ppcYuvResi = new TComYuv*[m_uiMaxDepth-1];
  m_ppcYuvReco = new TComYuv*[m_uiMaxDepth-1];
  m_ppcCU      = new TComDataCU*[m_uiMaxDepth-1];
#if COM16_C806_OBMC
  m_ppcTmpYuv1 = new TComYuv*[m_uiMaxDepth-1];
  m_ppcTmpYuv2 = new TComYuv*[m_uiMaxDepth-1];
#endif

  for ( UInt ui = 0; ui < m_uiMaxDepth-1; ui++ )
  {
    UInt uiNumPartitions = 1<<( ( m_uiMaxDepth - ui - 1 )<<1 );
    UInt uiWidth  = uiMaxWidth  >> ui;
    UInt uiHeight = uiMaxHeight >> ui;

    // The following arrays (m_ppcYuvResi, m_ppcYuvReco and m_ppcCU) are only required for CU depths
    // although data is allocated for all possible depths of the CU/TU tree except the last.
    // Since the TU tree will always include at least one additional depth greater than the CU tree,
    // there will be enough entries for these arrays.
    // (Section 7.4.3.2: "The CVS shall not contain data that result in (Log2MinTrafoSize) MinTbLog2SizeY
    //                    greater than or equal to MinCbLog2SizeY")
    // TODO: tidy the array allocation given the above comment.

    m_ppcYuvResi[ui] = new TComYuv;    m_ppcYuvResi[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
    m_ppcYuvReco[ui] = new TComYuv;    m_ppcYuvReco[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
    m_ppcCU     [ui] = new TComDataCU; m_ppcCU     [ui]->create( chromaFormatIDC, uiNumPartitions, uiWidth, uiHeight, true, uiMaxWidth >> (m_uiMaxDepth - 1) );
#if COM16_C806_OBMC
    m_ppcTmpYuv1[ui] = new TComYuv;    m_ppcTmpYuv1[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
    m_ppcTmpYuv2[ui] = new TComYuv;    m_ppcTmpYuv2[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
#endif
  }
#endif

  m_bDecodeDQP = false;
  m_IsChromaQpAdjCoded = false;

  // initialize partition order.
  UInt* piTmp = &g_auiZscanToRaster[0];
  initZscanToRaster(m_uiMaxDepth, 1, 0, piTmp);
  initRasterToZscan( uiMaxWidth, uiMaxHeight, m_uiMaxDepth );

  // initialize conversion matrix from partition index to pel
  initRasterToPelXY( uiMaxWidth, uiMaxHeight, m_uiMaxDepth );

#if PVC_AQS
  initFactorValue();
#endif
}

Void TDecCu::destroy()
{
#if JVET_C0024_QTBT
  UInt uiMaxWidth = 1<<MIN_CU_LOG2<<(m_uiMaxDepth-1);
  UInt uiMaxHeight = uiMaxWidth;
  UInt uiNumWidthIdx = g_aucConvertToBit[uiMaxWidth] + 1;
  UInt uiNumHeightIdx = g_aucConvertToBit[uiMaxHeight] + 1;

  for (UInt wIdx=0; wIdx<uiNumWidthIdx; wIdx++)
  {
    for (UInt hIdx=0; hIdx<uiNumHeightIdx; hIdx++)
    {
      m_pppcCU[wIdx][hIdx]->destroy(); delete m_pppcCU[wIdx][hIdx] ; m_pppcCU[wIdx][hIdx] = NULL;
      m_pppcYuvResi[wIdx][hIdx]->destroy(); delete m_pppcYuvResi[wIdx][hIdx]; m_pppcYuvResi[wIdx][hIdx] = NULL;
      m_pppcYuvReco[wIdx][hIdx]->destroy(); delete m_pppcYuvReco[wIdx][hIdx]; m_pppcYuvReco[wIdx][hIdx] = NULL;
#if COM16_C806_OBMC
      m_pppcTmpYuv1[wIdx][hIdx]->destroy(); delete m_pppcTmpYuv1[wIdx][hIdx]; m_pppcTmpYuv1[wIdx][hIdx] = NULL;
      m_pppcTmpYuv2[wIdx][hIdx]->destroy(); delete m_pppcTmpYuv2[wIdx][hIdx]; m_pppcTmpYuv2[wIdx][hIdx] = NULL;
#endif
    }
    delete [] m_pppcCU[wIdx]; m_pppcCU[wIdx] = NULL;
    delete [] m_pppcYuvResi[wIdx]; m_pppcYuvResi[wIdx] = NULL;
    delete [] m_pppcYuvReco[wIdx]; m_pppcYuvReco[wIdx] = NULL;
#if COM16_C806_OBMC
    delete [] m_pppcTmpYuv1[wIdx]; m_pppcTmpYuv1[wIdx] = NULL;
    delete [] m_pppcTmpYuv2[wIdx]; m_pppcTmpYuv2[wIdx] = NULL;
#endif
  }
  delete [] m_pppcCU; m_pppcCU = NULL;
  delete [] m_pppcYuvResi; m_pppcYuvResi = NULL;
  delete [] m_pppcYuvReco; m_pppcYuvReco = NULL;
#if COM16_C806_OBMC
  delete [] m_pppcTmpYuv1; m_pppcTmpYuv1 = NULL;
  delete [] m_pppcTmpYuv2; m_pppcTmpYuv2 = NULL;
#endif
#else
  for ( UInt ui = 0; ui < m_uiMaxDepth-1; ui++ )
  {
    m_ppcYuvResi[ui]->destroy(); delete m_ppcYuvResi[ui]; m_ppcYuvResi[ui] = NULL;
    m_ppcYuvReco[ui]->destroy(); delete m_ppcYuvReco[ui]; m_ppcYuvReco[ui] = NULL;
    m_ppcCU     [ui]->destroy(); delete m_ppcCU     [ui]; m_ppcCU     [ui] = NULL;
#if COM16_C806_OBMC
    m_ppcTmpYuv1[ui]->destroy(); delete m_ppcTmpYuv1[ui]; m_ppcTmpYuv1[ui] = NULL;
    m_ppcTmpYuv2[ui]->destroy(); delete m_ppcTmpYuv2[ui]; m_ppcTmpYuv2[ui] = NULL;
#endif
  }

  delete [] m_ppcYuvResi; m_ppcYuvResi = NULL;
  delete [] m_ppcYuvReco; m_ppcYuvReco = NULL;
  delete [] m_ppcCU     ; m_ppcCU      = NULL;
#if COM16_C806_OBMC
  delete [] m_ppcTmpYuv1; m_ppcTmpYuv1 = NULL;
  delete [] m_ppcTmpYuv2; m_ppcTmpYuv2 = NULL;
#endif
#endif
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/** 
 Parse a CTU.
 \param    pCtu                      [in/out] pointer to CTU data structure
 \param    isLastCtuOfSliceSegment   [out]    true, if last CTU of the slice segment
 */
Void TDecCu::decodeCtu( TComDataCU* pCtu, Bool& isLastCtuOfSliceSegment )
{
#if JVET_C0024_QTBT
  pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
#endif
  if ( pCtu->getSlice()->getPPS()->getUseDQP() )
  {
    setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX // init current quantization unit info.
    Char qp = pCtu->getCtuLastCodedQP();
    pCtu->setQuPartIdx(0); 
    pCtu->setQuLastCodedQP( qp );
    pCtu->setCodedQP( qp );
#endif
  }

  if ( pCtu->getSlice()->getUseChromaQpAdj() )
  {
    setIsChromaQpAdjCoded(true);
  }

  // start from the top level CU
#if BARC
  if (pCtu->getSlice()->isIntra())
  {
	  m_pcEntropyDecoder->decodeBARCFlag(pCtu, CHANNEL_TYPE_LUMA);
#if BARC_UP_SAMPLE_RDO
#if BARC_INTRA_LUMA_SPS
	  if (pCtu->getSlice()->getSPS()->getUseBarcIntraLuma())
	  {
#endif
		  if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
		  {
			  m_pcEntropyDecoder->decodeUpCNNFlag(pCtu, 0, 0);
		  }
#if BARC_INTRA_LUMA_SPS
	  }
#endif
#endif
  }
#endif
#if JVET_C0024_QTBT
  UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize();
  pCtu->getPic()->setCodedAreaInCTU(0);
  pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2);  //only used for affine merge code or not

#if TT_TREE_DS1
#if BARC
  UInt uiBARCCTUSize = uiCTUSize;
  if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
  {
    uiBARCCTUSize >>= BARC_DScaleRatio;
    pCtu->initDscaleCtuDec(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
    xDecodeCU( pCtu, 0, BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize, 0, isLastCtuOfSliceSegment);
  }
  else
#endif
  xDecodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, 0, isLastCtuOfSliceSegment);
#else
#if BARC
  UInt uiBARCCTUSize = uiCTUSize;
  if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
  {
    uiBARCCTUSize >>= BARC_DScaleRatio;
    pCtu->initDscaleCtuDec(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
    xDecodeCU(pCtu, 0, BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize,isLastCtuOfSliceSegment);
   }
  else
#endif
  xDecodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, isLastCtuOfSliceSegment);
#endif
#else
  xDecodeCU( pCtu, 0, 0, isLastCtuOfSliceSegment);
#endif
#if JVET_C0024_QTBT
  if (pCtu->getSlice()->isIntra())
  {
#if TRACE_FLAGS
	  gNumCu++;
#endif
    pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
    pCtu->getPic()->setCodedAreaInCTU(0);
#if JVET_E0077_ENHANCED_LM || DIMD_INTRA_PRED
    pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize >> MIN_CU_LOG2, uiCTUSize >> MIN_CU_LOG2);
#endif
#if JVET_C0024_DELTA_QP_FIX
    if ( pCtu->getSlice()->getPPS()->getUseDQP() )
    {
      setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX // init current quantization unit info.
      Char qp = pCtu->getCtuLastCodedQP();
      pCtu->setQuPartIdx(0); 
      pCtu->setQuLastCodedQP( qp );
      pCtu->setCodedQP( qp );
#endif
    }
#endif
#if JVET_E0062_MULTI_DMS
    pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2); 
#endif

    // start from the top level CU
#if TT_TREE_DS1
    xDecodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, 0, isLastCtuOfSliceSegment);
#else
    xDecodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, isLastCtuOfSliceSegment);
#endif
  }
  xFinishDecodeCU( pCtu, 0, 0, isLastCtuOfSliceSegment );
#endif
}

/** 
 Decoding process for a CTU.
 \param    pCtu                      [in/out] pointer to CTU data structure
 */
Void TDecCu::decompressCtu( TComDataCU* pCtu )
{
#if JVET_C0024_QTBT
  pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
  UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize();
  pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2);
  pCtu->getPic()->setCodedAreaInCTU(0);
#if BARC
  Bool isBARCFlag = pCtu->getBARCFlag(CHANNEL_TYPE_LUMA);
  if (isBARCFlag)
  {
    UInt uiBARCCTUSize = uiCTUSize>> BARC_DScaleRatio;
    xDecompressCU( pCtu, 0, BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize );
    upSampleRec(pCtu, CHANNEL_TYPE_LUMA);
  }
  else
  {
    xDecompressCU(pCtu, 0, 0, uiCTUSize, uiCTUSize);
    if (pCtu->getSlice()->isIntra())
    {
      downSampleRec(pCtu, CHANNEL_TYPE_LUMA);
    }
  }
#else
  xDecompressCU( pCtu, 0,  0, uiCTUSize, uiCTUSize ); 
#endif

  if (pCtu->getSlice()->isIntra())
  {
#if TRACE_FLAGS
	  gNumCu++;
#endif
    pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
    pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2);
    pCtu->getPic()->setCodedAreaInCTU(0);
    xDecompressCU( pCtu, 0,  0, uiCTUSize, uiCTUSize );
  }
#else
  xDecompressCU( pCtu, 0,  0 );
#endif
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

//! decode end-of-slice flag
Bool TDecCu::xDecodeSliceEnd( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiIsLastCtuOfSliceSegment;

#if JVET_C0024_QTBT
  m_pcEntropyDecoder->decodeTerminatingBit( uiIsLastCtuOfSliceSegment );
#else
  if (pcCU->isLastSubCUOfCtu(uiAbsPartIdx))
  {
    m_pcEntropyDecoder->decodeTerminatingBit( uiIsLastCtuOfSliceSegment );
  }
  else
  {
    uiIsLastCtuOfSliceSegment=0;
  }
#endif

  return uiIsLastCtuOfSliceSegment>0;
}

//! decode CU block recursively
#if JVET_C0024_QTBT
#if TT_TREE_DS1
Void TDecCu::xDecodeCU(TComDataCU*const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, const UInt uiWidth, const UInt uiHeight, UInt uiBTDepth, Bool &isLastCtuOfSliceSegment, UInt uiSplitConstrain)
#else
Void TDecCu::xDecodeCU(TComDataCU*const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, const UInt uiWidth, const UInt uiHeight, Bool &isLastCtuOfSliceSegment, UInt uiSplitConstrain)
#endif
#else
Void TDecCu::xDecodeCU(TComDataCU*const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, Bool &isLastCtuOfSliceSegment)
#endif
{
  TComPic* pcPic = pcCU->getPic();
  const TComSPS &sps = pcPic->getPicSym()->getSPS();
  const TComPPS &pps = pcPic->getPicSym()->getPPS();
#if !JVET_C0024_QTBT
  const UInt maxCuWidth = sps.getMaxCUWidth();
  const UInt maxCuHeight = sps.getMaxCUHeight();
#endif
  UInt uiCurNumParts = pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiQNumParts = uiCurNumParts >> 2;

  Bool bBoundary = false;
  UInt uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
#if JVET_C0024_QTBT
  UInt uiRPelX = uiLPelX + uiWidth - 1;
#else
  UInt uiRPelX = uiLPelX + (maxCuWidth >> uiDepth) - 1;
#endif
  UInt uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if JVET_C0024_QTBT
  UInt uiBPelY = uiTPelY + uiHeight - 1;
#else
  UInt uiBPelY = uiTPelY + (maxCuHeight >> uiDepth) - 1;
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
	if (pcCU->getSlice()->getSPS()->getBackgroundRefsPresent() && pcCU->getSlice()->getSliceType() != I_SLICE && uiDepth == 0 && uiAbsPartIdx == 0)
	{
		m_pcEntropyDecoder->decodeCUBgFlag(pcCU, 0, 0);
	}
#endif

#if BT_FORCE
#if BARC_BTFORCE
  UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
  UInt deltaY = (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#else
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#endif

#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(pcCU->getTextType()) ? pcCU->getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcCU->getSlice()->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif

#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  const UInt uiDIMDMaxBtDepht = (pcCU->getBoundaryCtuLocation() && pcCU->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (pcCU->getBoundaryCtuLocation() && pcCU->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
#endif

#if JVET_C0024_BT_RMV_REDUNDANT
  pcCU->setSplitConstrain(uiSplitConstrain);
  Bool bQTreeValid = false;
#endif

#if JVET_C0024_QTBT
  Bool bForceQT = uiWidth > MAX_TU_SIZE;
  if (bForceQT)
  {
    assert(uiWidth == uiHeight);
  }
  UInt uiCTUSize = pcCU->getSlice()->getSPS()->getCTUSize();
#if JVET_C0024_DELTA_QP_FIX
  UInt uiQTWidth = uiCTUSize >> uiDepth;
  UInt uiQTHeight = uiCTUSize >> uiDepth;
  const UInt uiQTBTDepth = (uiDepth << 1) + (g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight]);
  const UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif
  if (uiCTUSize >> uiDepth == uiWidth && uiWidth == uiHeight)
  {
#endif
#if BARC_Boundary
    UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
    UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

    UInt deltaX =  (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
    UInt deltaY =  (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());
    if( ( uiRPelX < sps.getPicWidthInLumaSamples()- deltaX ) && ( uiBPelY < sps.getPicHeightInLumaSamples() - deltaY) )
#else
    if ((uiRPelX < sps.getPicWidthInLumaSamples()) && (uiBPelY < sps.getPicHeightInLumaSamples()))
#endif
    {
#if JVET_C0024_QTBT
      if (bForceQT)
      {
        pcCU->setDepthSubParts(uiDepth + 1, uiAbsPartIdx);
      }
      else
#endif
        m_pcEntropyDecoder->decodeSplitFlag(pcCU, uiAbsPartIdx, uiDepth);
    }
    else
    {
      bBoundary = true;
#if BT_FORCE
      if (
#if !BT_FORCE_INTRA_ON
        !pcCU->getSlice()->isIntra() && 
#endif
        ((bBBoundary || bRBoundary) && (uiDepth >= iMaxQtLv)))
        pcCU->setDepthSubParts(uiDepth, uiAbsPartIdx);
#endif
    }

#if JVET_C0024_BT_RMV_REDUNDANT
    bQTreeValid = true;

    UInt uiMinQTSize = sps.getMinQTSize(pcCU->getSlice()->getSliceType(), pcCU->getTextType());

    if ((uiCTUSize >> uiDepth) <= uiMinQTSize)
    {
      bQTreeValid = false;
    }
#endif

#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL
    if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiDepth))
    {
#if JVET_C0024_QTBT
      pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth); // set size for correct flags assignment
#endif
      m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
    }
#endif

#if JVET_C0024_QTBT
    if (((uiDepth < pcCU->getDepth(uiAbsPartIdx)) || bBoundary)
#if BT_FORCE
      && (
#if !BT_FORCE_INTRA_ON
        pcCU->getSlice()->isIntra() || 
#endif
        (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
      )
#else
    if (((uiDepth < pcCU->getDepth(uiAbsPartIdx)) && (uiDepth < sps.getLog2DiffMaxMinCodingBlockSize())) || bBoundary)
#endif
    {
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL
      if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiDepth))
      {
#if JVET_C0024_QTBT
        pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth); // set size for correct flags assignment
#endif
        m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
      }
#endif
      UInt uiIdx = uiAbsPartIdx;
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
      if (uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx(uiAbsPartIdx);
        pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#else
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiDepth); // set QP to default QP
#endif
      }

      if (uiDepth == pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcCU->getSlice()->getUseChromaQpAdj())
      {
        setIsChromaQpAdjCoded(true);
      }

      for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++)
      {
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiIdx]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiIdx]];
#if BARC_Boundary
        UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX =  uiPelXInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY =  uiPelYInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        if ( !isLastCtuOfSliceSegment && ( uiLPelX + deltaX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY + deltaY < sps.getPicHeightInLumaSamples() ) )
#else
        if (!isLastCtuOfSliceSegment && (uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#if JVET_C0024_QTBT
#if TT_TREE_DS1
          xDecodeCU(pcCU, uiIdx, uiDepth + 1, uiWidth >> 1, uiHeight >> 1, uiBTDepth, isLastCtuOfSliceSegment);
#else
          xDecodeCU(pcCU, uiIdx, uiDepth + 1, uiWidth >> 1, uiHeight >> 1, isLastCtuOfSliceSegment);
#endif
#else
          xDecodeCU(pcCU, uiIdx, uiDepth + 1, isLastCtuOfSliceSegment);
#endif
        }
        else
        {
          pcCU->setOutsideCUPart(uiIdx, uiDepth + 1);
#if JVET_C0024_QTBT
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 2);
#endif
        }

        uiIdx += uiQNumParts;
      }
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
      if (uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
      {
        if (getdQPFlag())
        {
          UInt uiQPSrcPartIdx = uiAbsPartIdx;
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#else
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiDepth); // set QP to default QP
#endif
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiQPSrcPartIdx));
#endif
        }
      }
      return;
    }
#if JVET_C0024_QTBT
  }
#if JVET_C0024_BT_RMV_REDUNDANT  
  Bool bBTHorRmvEnable = false;
  Bool bBTVerRmvEnable = false;
  if (pcCU->getSlice()->getSliceType() != I_SLICE)
  {
    bBTHorRmvEnable = true;
    bBTVerRmvEnable = bQTreeValid;
  }
#endif

#if JVET_C0024_SPS_MAX_BT_DEPTH
  UInt uiMaxBTD = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? sps.getMaxBTDepthISliceL() : sps.getMaxBTDepthISliceC()) : sps.getMaxBTDepth();
#if BT_FORCE_EXT
  Int iMaxELv =
#if BT_FORCE_INTRA_ON
    pcCU->getSlice()->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;
  if (
#if !BT_FORCE_INTRA_ON
    !pcCU->getSlice()->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx)&&
#endif
    pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
#else
  UInt uiMaxBTD = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? MAX_BT_DEPTH : MAX_BT_DEPTH_C) : MAX_BT_DEPTH_INTER;
#if BT_FORCE_EXT
  Int iMaxELv = 
#if BT_FORCE_INTRA_ON
    pcCU->getSlice()->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;
  if (
#if !BT_FORCE_INTRA_ON
    !pcCU->getSlice()->isIntra() &&
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx) &&
#endif
    pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
#endif
#if JVET_C0024_SPS_MAX_BT_SIZE
  UInt uiMaxBTSize = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? sps.getMaxBTSizeISliceL() : sps.getMaxBTSizeISliceC()) : sps.getMaxBTSize();
#else
  UInt uiMaxBTSize = isLuma(pcCU->getTextType()) ? pcCU->getSlice()->getMaxBTSize() : MAX_BT_SIZE_C;
#endif
  UInt uiMinBTSize = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? MIN_BT_SIZE : MIN_BT_SIZE_C) : MIN_BT_SIZE_INTER;
#if !TT_TREE_DS1
  UInt uiBTDepth = pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);
#endif

#if TT_SYNTAX //extend conditional check
  Bool bAdditionalCond = true;
#if TT_LIMIT_SIZEIDX
  UInt uiDeltaSizeIdx = g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight];
  bAdditionalCond = bAdditionalCond && (uiDeltaSizeIdx < (uiMaxBTD + 1));
#endif
#if TT_LEAF_NODE
  UInt uiBTSplitModeLastBtDepth = (uiBTDepth == 0) ? 0 : pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth - 1);
  if (uiBTSplitModeLastBtDepth & 2) //if tt
    bAdditionalCond = false;
#endif
  if (((uiHeight > uiMinBTSize || uiWidth > uiMinBTSize)
    && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD && bAdditionalCond)
#if BT_FORCE
    || (
#if !BT_FORCE_INTRA_ON
      !pcCU->getSlice()->isIntra() &&
#endif
       ((bBBoundary || bRBoundary) && (uiDepth >= iMaxQtLv)))
#endif
    )
#else
  if (((uiHeight > uiMinBTSize || uiWidth > uiMinBTSize)
    && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD)
#if BT_FORCE
    || (
#if !BT_FORCE_INTRA_ON
      !pcCU->getSlice()->isIntra() && 
#endif
      ((bBBoundary || bRBoundary) && (uiDepth >= iMaxQtLv)))
#endif
    )
#endif
  {
    UInt uiIdxBT = uiAbsPartIdx;
#if BT_FORCE
#if !BT_FORCE_INTRA_ON
    if (pcCU->getSlice()->isIntra())
    {
      m_pcEntropyDecoder->decodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight
#if TT_TREE_DS1
        , uiBTDepth
#endif
        );
    }
    else
    {
#endif
      if (bBBoundary)
      {
#if TT_TREE_DS1
#if BT_FORCE_TT
        if (pcCU->isForcedBT(uiAbsPartIdx, uiWidth, uiHeight))
        {
#if (BT_FORCE_AR>0)
          if (uiWidth >= uiHeight*BT_FORCE_AR)
            pcCU->setBTSplitModeSubParts(ST_VERT_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
          else
#endif
          pcCU->setBTSplitModeSubParts(ST_HORI_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        }
        else
        {
#if (BT_FORCE_AR_TT>0)
          if (uiWidth >= uiHeight*BT_FORCE_AR_TT)
            pcCU->setBTSplitModeSubParts(ST_VERT_TT_L, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
          else
#endif
          pcCU->setBTSplitModeSubParts(ST_HORI_TT_L, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        }
#else
#if (BT_FORCE_AR>0)
        if (uiWidth >= uiHeight*BT_FORCE_AR)
          pcCU->setBTSplitModeSubParts(ST_VERT_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        else
#endif
          pcCU->setBTSplitModeSubParts(ST_HORI_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
#endif
#else
#if (BT_FORCE_AR>0)
        if (uiWidth >= uiHeight*BT_FORCE_AR)
          pcCU->setBTSplitModeSubParts(2, uiAbsPartIdx, uiWidth, uiHeight);
        else
#endif
        pcCU->setBTSplitModeSubParts(1, uiAbsPartIdx, uiWidth, uiHeight);
#endif
      }
      else if (bRBoundary)
      {
#if TT_TREE_DS1
#if BT_FORCE_TT
        if (pcCU->isForcedBT(uiAbsPartIdx, uiWidth, uiHeight))
        {
#if (BT_FORCE_AR>0)
          if (uiHeight >= uiWidth*BT_FORCE_AR)
            pcCU->setBTSplitModeSubParts(ST_HORI_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
          else
#endif
          pcCU->setBTSplitModeSubParts(ST_VERT_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        }
        else
        {
#if (BT_FORCE_AR_TT>0)
          if (uiHeight >= uiWidth*BT_FORCE_AR_TT)
            pcCU->setBTSplitModeSubParts(ST_HORI_TT_L, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
          else
#endif
          pcCU->setBTSplitModeSubParts(ST_VERT_TT_L, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        }
#else
#if (BT_FORCE_AR>0)
        if (uiHeight >= uiWidth*BT_FORCE_AR)
          pcCU->setBTSplitModeSubParts(ST_HORI_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
        else
#endif
          pcCU->setBTSplitModeSubParts(ST_VERT_BT, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
#endif
#else
#if (BT_FORCE_AR>0)
        if (uiHeight >= uiWidth*BT_FORCE_AR)
          pcCU->setBTSplitModeSubParts(1, uiAbsPartIdx, uiWidth, uiHeight);
        else
#endif
        pcCU->setBTSplitModeSubParts(2, uiAbsPartIdx, uiWidth, uiHeight);
#endif
      }
      else
      {
        m_pcEntropyDecoder->decodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight
#if TT_TREE_DS1
          , uiBTDepth
#endif
          );
      }
#if !BT_FORCE_INTRA_ON
    }
#endif
#else
    m_pcEntropyDecoder->decodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight
#if TT_TREE_DS1
      , uiBTDepth
#endif
      );
#endif
#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
    if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiQTBTDepth))
    {
      pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth); // set size for correct flags assignment
      m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
    }
#endif
#if JVET_C0024_BT_RMV_REDUNDANT
    uiSplitConstrain = 0;
#endif

#if TT_TREE_DS1
    if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth) == ST_HORI_BT)
#else
    if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth) == 1)
#endif
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyDecoder->decodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);

        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
        if (!pcCU->getSlice()->isIntra())
        {
          pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);

          if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
          {
            pcCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
          }
        }

      }
#endif

#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx(uiAbsPartIdx);
        pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
      }
#endif
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
      if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, /*uiDepth*/uiQTBTDepth))
      {
        pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth); // set size for correct flags assignment
        UInt uiBlkX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
        UInt uiBlkY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
        pcCU->getPic()->setCodedBlkInCTU(true, uiBlkX, uiBlkY, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
        m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
      }
#endif
      for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++)
      {
        if (uiPartUnitIdx == 1)
        {
          uiIdxBT = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT]
            + (uiHeight >> 1) / pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
        }
#if JVET_C0024_BT_RMV_REDUNDANT
#if BT_FORCE
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiIdxBT]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiIdxBT]];
#if BARC_BTFORCE
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if (!isLastCtuOfSliceSegment && (uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if (!isLastCtuOfSliceSegment && (uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
#if TT_TREE_DS1
          xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth, uiHeight >> 1, uiBTDepth + 1, isLastCtuOfSliceSegment, uiSplitConstrain);
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == ST_VERT_BT && bBTHorRmvEnable && uiPartUnitIdx == 0)
#else
          xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth, uiHeight >> 1, isLastCtuOfSliceSegment, uiSplitConstrain);
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == 2 && bBTHorRmvEnable && uiPartUnitIdx == 0)
#endif
          {
            uiSplitConstrain = 2;
          }
#if TT_RMV_PART1_BTBT && !TT_LEAF_NODE
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == ST_HORI_BT && uiPartUnitIdx == 0)
            uiSplitConstrain = 1;
#endif
#if BT_FORCE
        }
        else
        {
          pcCU->setSizeSubParts(uiWidth, uiHeight >> 1, uiIdxBT, uiDepth);
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
        }
#endif
#else
        xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth, uiHeight >> 1, isLastCtuOfSliceSegment);
#endif
      }
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        if (getdQPFlag())
        {
          UInt uiQPSrcPartIdx = uiAbsPartIdx;
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiQPSrcPartIdx));
#endif
        }
      }
#endif
      return;
    }
#if TT_TREE_DS1
    else if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth) == ST_VERT_BT)
#else
    else if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth) == 2)
#endif
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyDecoder->decodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);

        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
        if (!pcCU->getSlice()->isIntra())
        {
          pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);

          if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
          {
            pcCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
          }
        }
      }
#endif

#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx(uiAbsPartIdx);
        pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
      }
#endif
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
      if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, /*uiDepth*/uiQTBTDepth))
      {
        pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth); // set size for correct flags assignment
        UInt uiBlkX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
        UInt uiBlkY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
        pcCU->getPic()->setCodedBlkInCTU(true, uiBlkX, uiBlkY, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
        m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
      }
#endif
      for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++)
      {
        if (uiPartUnitIdx == 1)
        {
          uiIdxBT = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT]
            + (uiWidth >> 1) / pcCU->getPic()->getMinCUWidth()];
        }
#if JVET_C0024_BT_RMV_REDUNDANT
#if BT_FORCE
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiIdxBT]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiIdxBT]];
#if BARC_BTFORCE
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if (!isLastCtuOfSliceSegment && (uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if (!isLastCtuOfSliceSegment && (uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
#if TT_TREE_DS1
          xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth >> 1, uiHeight, uiBTDepth + 1, isLastCtuOfSliceSegment, uiSplitConstrain);
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == ST_HORI_BT && bBTVerRmvEnable && uiPartUnitIdx == 0)
#else
          xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth >> 1, uiHeight, isLastCtuOfSliceSegment, uiSplitConstrain);
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == 1 && bBTVerRmvEnable && uiPartUnitIdx == 0)
#endif
          {
            uiSplitConstrain = 1;
          }
#if TT_RMV_PART1_BTBT && !TT_LEAF_NODE
          if (pcCU->getBTSplitModeForBTDepth(uiIdxBT, uiBTDepth + 1) == ST_VERT_BT && uiPartUnitIdx == 0)
            uiSplitConstrain = 2;
#endif
#if BT_FORCE
        }
        else
        {
          pcCU->setSizeSubParts(uiWidth >> 1, uiHeight, uiIdxBT, uiDepth);
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
        }
#endif
#else
        xDecodeCU(pcCU, uiIdxBT, uiDepth, uiWidth >> 1, uiHeight, isLastCtuOfSliceSegment);
#endif
      }
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        if (getdQPFlag())
        {
          UInt uiQPSrcPartIdx = uiAbsPartIdx;
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiQPSrcPartIdx));
#endif
        }
      }
#endif
      return;
    }
#if TT_DEC_OPT //xDecodeCU
    UInt uiHeightSubPart, uiWidthSubPart, uiIdxBTSubPart;
    UInt uiBTSplitMode = pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth);
    if (uiBTSplitMode == ST_HORI_TT_L || uiBTSplitMode == ST_HORI_TT_S)
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyDecoder->decodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);

        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
        if (!pcCU->getSlice()->isIntra())
        {
          pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);

          if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
          {
            pcCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
          }
        }
      }
#endif

#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx(uiAbsPartIdx);
        pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
      }
#endif

      for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++)
      {
        if (uiPartUnitIdx == 1)
        {
          uiIdxBTSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT] + (uiHeight >> 2) / pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
          uiHeightSubPart = uiHeight >> 1;
        }
        else
        {
          if (uiPartUnitIdx == 2)
            uiIdxBTSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT] + ((uiHeight >> 2) + (uiHeight >> 1)) / pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
          else
            uiIdxBTSubPart = uiIdxBT;
          uiHeightSubPart = uiHeight >> 2;
        }
#if BT_FORCE_TT
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiIdxBTSubPart]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiIdxBTSubPart]];
#if BARC_Boundary
        UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        if (!isLastCtuOfSliceSegment && (uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if (!isLastCtuOfSliceSegment && (uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xDecodeCU(pcCU, uiIdxBTSubPart, uiDepth, uiWidth, uiHeightSubPart, uiBTDepth + 1, isLastCtuOfSliceSegment, 0);
#if BT_FORCE_TT
        }
        else
        {
          pcCU->setSizeSubParts(uiWidth, uiHeightSubPart, uiIdxBTSubPart, uiDepth);
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeightSubPart);
        }
#endif
      }
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        if (getdQPFlag())
        {
          UInt uiQPSrcPartIdx = uiAbsPartIdx;
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiQPSrcPartIdx));
#endif
        }
      }
#endif
      return;
    }

    if (uiBTSplitMode == ST_VERT_TT_L || uiBTSplitMode == ST_VERT_TT_S)
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyDecoder->decodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);

        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
        if (!pcCU->getSlice()->isIntra())
        {
          pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);

          if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
          {
            pcCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
          }
        }
      }
#endif

#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx(uiAbsPartIdx);
        pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
        pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
      }
#endif

      for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++)
      {
        if (uiPartUnitIdx == 1)
        {
          uiIdxBTSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT] + (uiWidth >> 2) / pcCU->getPic()->getMinCUWidth()];
          uiWidthSubPart = uiWidth >> 1;
        }
        else
        {
          if (uiPartUnitIdx == 2)
            uiIdxBTSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiIdxBT] + ((uiWidth >> 2) + (uiWidth >> 1)) / pcCU->getPic()->getMinCUWidth()];
          else
            uiIdxBTSubPart = uiIdxBT;
          uiWidthSubPart = uiWidth >> 2;
        }
#if BT_FORCE_TT
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiIdxBTSubPart]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiIdxBTSubPart]];
#if BARC_Boundary
        UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTU*pcCU->getBARCFlag(pcCU->getTextType());
        if (!isLastCtuOfSliceSegment && (uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if (!isLastCtuOfSliceSegment && (uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xDecodeCU(pcCU, uiIdxBTSubPart, uiDepth, uiWidthSubPart, uiHeight, uiBTDepth + 1, isLastCtuOfSliceSegment, 0);
#if BT_FORCE_TT
        }
        else
        {
          pcCU->setSizeSubParts(uiWidthSubPart, uiHeight, uiIdxBTSubPart, uiDepth);
          pcCU->getPic()->addCodedAreaInCTU(uiWidthSubPart*uiHeight);
        }
#endif
      }
#if JVET_C0024_DELTA_QP_FIX
      if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        if (getdQPFlag())
        {
          UInt uiQPSrcPartIdx = uiAbsPartIdx;
          pcCU->setQPSubParts(pcCU->getRefQP(uiQPSrcPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiQPSrcPartIdx));
#endif
        }
      }
#endif
      return;
    }
#endif
  }
#if  TRACE_FLAGS 
  if (pcCU->getSlice()->getTextType() == CHANNEL_TYPE_LUMA)
    gNumCu++;
#endif
  pcCU->setSizeSubParts(uiWidth, uiHeight, uiAbsPartIdx, uiDepth);
  UInt uiBlkX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
  UInt uiBlkY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]] >> MIN_CU_LOG2;
  pcCU->getPic()->setCodedBlkInCTU(true, uiBlkX, uiBlkY, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
  pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight);

  UInt uiWidthIdx = g_aucConvertToBit[uiWidth];
  UInt uiHeightIdx = g_aucConvertToBit[uiHeight];
  assert(g_aucConvertToBit[uiWidth] >= 0 && g_aucConvertToBit[uiHeight] >= 0);
#endif

#if JVET_C0024_DELTA_QP_FIX
  if (uiQTBTDepth <= uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
  if (uiDepth <= pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
  {
    setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
    pcCU->setQuPartIdx(uiAbsPartIdx);
    pcCU->setQuLastCodedQP(pcCU->getCodedQP());
#endif
#if JVET_C0024_DELTA_QP_FIX
    pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
#else
    pcCU->setQPSubParts(pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiDepth); // set QP to default QP
#endif
  }

  if (uiDepth <= pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcCU->getSlice()->getUseChromaQpAdj())
  {
    setIsChromaQpAdjCoded(true);
  }

  if (pps.getTransquantBypassEnableFlag())
  {
    m_pcEntropyDecoder->decodeCUTransquantBypassFlag(pcCU, uiAbsPartIdx, uiDepth);
  }

  // decode CU mode and the partition size
  if (!pcCU->getSlice()->isIntra())
  {
    m_pcEntropyDecoder->decodeSkipFlag(pcCU, uiAbsPartIdx, uiDepth);
  }
#if COM16_C806_OBMC
  pcCU->setOBMCFlagSubParts(true, uiAbsPartIdx, uiDepth);
#endif

  if (pcCU->isSkipped(uiAbsPartIdx))
  {
#if JVET_C0024_QTBT
    m_pppcCU[uiWidthIdx][uiHeightIdx]->copyInterPredInfoFrom(pcCU, uiAbsPartIdx, REF_PIC_LIST_0);
    m_pppcCU[uiWidthIdx][uiHeightIdx]->copyInterPredInfoFrom(pcCU, uiAbsPartIdx, REF_PIC_LIST_1);
#else
    m_ppcCU[uiDepth]->copyInterPredInfoFrom(pcCU, uiAbsPartIdx, REF_PIC_LIST_0);
    m_ppcCU[uiDepth]->copyInterPredInfoFrom(pcCU, uiAbsPartIdx, REF_PIC_LIST_1);
#endif
#if VCEG_AZ07_FRUC_MERGE
    m_pcEntropyDecoder->decodeFRUCMgrMode(pcCU, uiAbsPartIdx, uiDepth, 0);
    if (!pcCU->getFRUCMgrMode(uiAbsPartIdx))
    {
#endif
#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT
      TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
#if COM16_C1016_AFFINE
      TComMvField cAffineMvField[2][3]; // double length for mv of both lists, 3 mv for affine
#endif

      UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP 
      UChar eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
      memset(eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, sizeof(UChar)*MRG_MAX_NUM_CANDS);
#endif
#if VCEG_AZ06_IC
      Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif
      Int numValidMergeCand = 0;
#if JVET_C0024_QTBT
      for (UInt ui = 0; ui < m_pppcCU[uiWidthIdx][uiHeightIdx]->getSlice()->getMaxNumMergeCand(); ++ui)
#else
      for (UInt ui = 0; ui < m_ppcCU[uiDepth]->getSlice()->getMaxNumMergeCand(); ++ui)
#endif
      {
        uhInterDirNeighbours[ui] = 0;
      }
#endif

#if MRG_CONSIDER_COMPLEX
      if ((pcCU->getWidth(uiAbsPartIdx) >= 16) || (pcCU->getHeight(uiAbsPartIdx) >= 16))
      {
        m_pcEntropyDecoder->decodeMergeComplexFlag(pcCU, uiAbsPartIdx, uiDepth, 0);
        if (pcCU->getMergeComplexFlag(uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeMergeComplexIndex(pcCU, 0, uiAbsPartIdx, uiDepth);
        }
      }

      if (!pcCU->getMergeComplexFlag(uiAbsPartIdx))
      {
#endif

#if COM16_C1016_AFFINE
        if (pcCU->isAffineMrgFlagCoded(uiAbsPartIdx, 0))
        {
          m_pcEntropyDecoder->decodeAffineFlag(pcCU, uiAbsPartIdx, uiDepth, 0);
        }

        if (!pcCU->isAffine(uiAbsPartIdx))
        {
          m_pcEntropyDecoder->decodeMergeIndex(pcCU, 0, uiAbsPartIdx, uiDepth);
        }
#else
        m_pcEntropyDecoder->decodeMergeIndex(pcCU, 0, uiAbsPartIdx, uiDepth);
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM && VCEG_AZ07_FRUC_MERGE
        if (pcCU->isAffine(uiAbsPartIdx))
        {
          TComMvField cAffineMvField[2][3]; // double length for mv of both lists, 3 mv for affine
          UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
          Int numValidMergeCand = 0;
          UInt neighborAffineParam = 0;
#if JVET_C0024_QTBT
          for (UInt ui = 0; ui < m_pppcCU[uiWidthIdx][uiHeightIdx]->getSlice()->getMaxNumMergeCand(); ++ui)
#else
          for (UInt ui = 0; ui < m_ppcCU[uiDepth]->getSlice()->getMaxNumMergeCand(); ++ui)
#endif
          {
            uhInterDirNeighbours[ui] = 0;
          }
#if JVET_C0024_QTBT
          m_pppcCU[uiWidthIdx][uiHeightIdx]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#else
          m_ppcCU[uiDepth]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#endif
          assert(neighborAffineParam == 0 || neighborAffineParam == 1);
          pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiAbsPartIdx, 0, uiDepth);
        }
#endif

#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT
        UInt uiMergeIndex = pcCU->getMergeIndex(uiAbsPartIdx);

#if COM16_C1016_AFFINE
        if (pcCU->isAffine(uiAbsPartIdx))
        {
#if JVECT_C0062_AFFINE_SIX_PARAM
          UInt neighborAffineParam = 0;
#if JVET_C0024_QTBT
          m_pppcCU[uiWidthIdx][uiHeightIdx]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#else
          m_ppcCU[uiDepth]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#endif
          assert(neighborAffineParam == 0 || neighborAffineParam == 1);
          pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiAbsPartIdx, 0, uiDepth);
#else
#if JVET_C0024_QTBT
          m_pppcCU[uiWidthIdx][uiHeightIdx]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand);
#else
          m_ppcCU[uiDepth]->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand);
#endif
#endif
          pcCU->setInterDirSubParts(uhInterDirNeighbours[uiMergeIndex], uiAbsPartIdx, 0, uiDepth);

#if JVECT_C0062_AFFINE_SIX_PARAM
          Bool b6Param = neighborAffineParam;//ZF
#endif
          TComMv cTmpMv(0, 0);
          for (UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
          {
            if (pcCU->getSlice()->getNumRefIdx(RefPicList(uiRefListIdx)) > 0)
            {
              TComMvField* pcMvField = cAffineMvField[2 * uiMergeIndex + uiRefListIdx];

              pcCU->setMVPIdxSubParts(0, RefPicList(uiRefListIdx), uiAbsPartIdx, 0, uiDepth);
              pcCU->setMVPNumSubParts(0, RefPicList(uiRefListIdx), uiAbsPartIdx, 0, uiDepth);
              pcCU->getCUMvField(RefPicList(uiRefListIdx))->setAllMvd(cTmpMv, SIZE_2Nx2N, uiAbsPartIdx, uiDepth);
#if JVECT_C0062_AFFINE_SIX_PARAM
              pcCU->setAllAffineMvField(uiAbsPartIdx, 0, pcMvField, RefPicList(uiRefListIdx), uiDepth, b6Param);
#else
              pcCU->setAllAffineMvField(uiAbsPartIdx, 0, pcMvField, RefPicList(uiRefListIdx), uiDepth);
#endif
            }
          }
          xFinishDecodeCU(pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment);
          return;
        }
#endif
#if JVET_C0024_QTBT
        m_pppcCU[uiWidthIdx][uiHeightIdx]->getInterMergeCandidates(0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#else
        m_ppcCU[uiDepth]->getInterMergeCandidates(0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#endif
#if VCEG_AZ06_IC
          , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          , eMergeCandTypeNieghors
          , m_pMvFieldSP
          , m_phInterDirSP
          , uiAbsPartIdx
          , pcCU
#endif
          , uiMergeIndex);

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
        pcCU->setMergeTypeSubParts(eMergeCandTypeNieghors[uiMergeIndex], uiAbsPartIdx, 0, uiDepth);
        if (eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_DEFAULT_N)
#endif
          pcCU->setInterDirSubParts(uhInterDirNeighbours[uiMergeIndex], uiAbsPartIdx, 0, uiDepth);
#if VCEG_AZ06_IC
        pcCU->setICFlagSubParts(pcCU->getSlice()->getApplyIC() ? abICFlag[uiMergeIndex] : 0, uiAbsPartIdx, uiDepth);
#endif
        TComMv cTmpMv(0, 0);
        for (UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
        {
          if (pcCU->getSlice()->getNumRefIdx(RefPicList(uiRefListIdx)) > 0)
          {
            pcCU->setMVPIdxSubParts(0, RefPicList(uiRefListIdx), uiAbsPartIdx, 0, uiDepth);
            pcCU->setMVPNumSubParts(0, RefPicList(uiRefListIdx), uiAbsPartIdx, 0, uiDepth);
            pcCU->getCUMvField(RefPicList(uiRefListIdx))->setAllMvd(cTmpMv, SIZE_2Nx2N, uiAbsPartIdx, uiDepth);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
            if (eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_DEFAULT_N)
            {
#endif
              pcCU->getCUMvField(RefPicList(uiRefListIdx))->setAllMvField(cMvFieldNeighbours[2 * uiMergeIndex + uiRefListIdx], SIZE_2Nx2N, uiAbsPartIdx, uiDepth);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
            }
#endif
          }
        }
#endif

#if MRG_CONSIDER_COMPLEX
      }
#endif

#if VCEG_AZ07_FRUC_MERGE
    }
#if VCEG_AZ06_IC
    m_pcEntropyDecoder->decodeICFlag(pcCU, uiAbsPartIdx, uiDepth);
#endif
#endif
#if SIF_CU_LEVEL
    m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
#endif
#if JVET_C0024_QTBT
    if (pcCU->getSlice()->getPPS()->getUseDQP())
    {
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiWidth, uiHeight); // set QP
#else
      pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiDepth); // set QP
#endif
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setCodedQP(pcCU->getQP(uiAbsPartIdx));
#endif
    }

    if (pcCU->getSlice()->getUseChromaQpAdj() && !getIsChromaQpAdjCoded())
    {
      pcCU->setChromaQpAdjSubParts(pcCU->getCodedChromaQpAdj(), uiAbsPartIdx, uiDepth); // set QP
    }
#else
    xFinishDecodeCU(pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment);
#endif
    return;
  }
#if BARC_AdpQP
  if (pcCU->getBARCFlag(pcCU->getTextType()))
  {
    if (pcCU->getTextType() == CHANNEL_TYPE_LUMA)
      pcCU->setQPSubParts(pcCU->getQP(uiAbsPartIdx) - BARC_AdpQPSTEP, uiAbsPartIdx, uiWidth, uiHeight); // set QP to default QP
  }
#endif
  m_pcEntropyDecoder->decodePredMode(pcCU, uiAbsPartIdx, uiDepth);

#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
  if (isLuma(pcCU->getTextType()) && (uiBTDepth <= uiDIMDMaxBtDepht || !pcCU->getSlice()->isIntra()))
#else
  if (isLuma(pcCU->getTextType()) && (uiBTDepth <= DIMD_MAX_BT_DEPTH_SIG_AI || !pcCU->getSlice()->isIntra()))
#endif
  {
    m_pcEntropyDecoder->decodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
  }
#endif

#if JVET_C0024_QTBT
  if (isLuma(pcCU->getTextType()))
  {
#endif
#if VCEG_AZ05_INTRA_MPI
    m_pcEntropyDecoder->decodeMPIIdx(pcCU, uiAbsPartIdx, uiDepth);
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    m_pcEntropyDecoder->decodePDPCIdx(pcCU, uiAbsPartIdx, uiDepth);
#endif
#if JVET_C0024_QTBT
  }
#else
    m_pcEntropyDecoder->decodePartSize(pcCU, uiAbsPartIdx, uiDepth);
#endif

#if DIMD_INTRA_PRED
    if (pcCU->isIntra(uiAbsPartIdx))
    {
      if (isLuma(pcCU->getTextType()))  // luma of I-slice, and luma and chroma of P/B-slice
      {
        UChar ucDIMDFlag = pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
        if (!pcCU->getSlice()->isIntra())
        {
          pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
        }
        UChar ucSecIntraMode = pcCU->generateSecIntraDir(uiAbsPartIdx, uiWidth, uiHeight);
        pcCU->setSecIntraDirSubParts(ucSecIntraMode, uiAbsPartIdx, uiDepth);
      }
      else  // chroma of I-slice
      {
        pcCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, uiAbsPartIdx, uiDepth);
      }
    }
#endif

#if JVET_C0024_QTBT
    if (pcCU->isIntra(uiAbsPartIdx))
#else
    if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N)
#endif
    {
      m_pcEntropyDecoder->decodeIPCMInfo(pcCU, uiAbsPartIdx, uiDepth);

      if (pcCU->getIPCMFlag(uiAbsPartIdx))
      {
#if JVET_C0024_QTBT
        if (pcCU->getSlice()->getPPS()->getUseDQP())
        {
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiWidth, uiHeight); // set QP
#else
          pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiDepth); // set QP
#endif
#if JVET_C0024_DELTA_QP_FIX
          pcCU->setCodedQP(pcCU->getQP(uiAbsPartIdx));
#endif
        }

        if (pcCU->getSlice()->getUseChromaQpAdj() && !getIsChromaQpAdjCoded())
        {
          pcCU->setChromaQpAdjSubParts(pcCU->getCodedChromaQpAdj(), uiAbsPartIdx, uiDepth); // set QP
        }
#else
        xFinishDecodeCU(pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment);
#endif
        return;
      }
    }

    // prediction mode ( Intra : direction mode, Inter : Mv, reference idx )
#if JVET_C0024_QTBT
    m_pcEntropyDecoder->decodePredInfo(pcCU, uiAbsPartIdx, uiDepth, m_pppcCU[uiWidthIdx][uiHeightIdx]);
#else
    m_pcEntropyDecoder->decodePredInfo(pcCU, uiAbsPartIdx, uiDepth, m_ppcCU[uiDepth]);
#endif
#if COM16_C806_OBMC
    m_pcEntropyDecoder->decodeOBMCFlag(pcCU, uiAbsPartIdx, uiDepth);
#endif
#if SIF_CU_LEVEL
    m_pcEntropyDecoder->decodeAdaptiveIfIndex(pcCU, uiAbsPartIdx, uiDepth);
#endif
#if VCEG_AZ06_IC
    m_pcEntropyDecoder->decodeICFlag(pcCU, uiAbsPartIdx, uiDepth);
#endif
    // Coefficient decoding
    Bool bCodeDQP = getdQPFlag();
    Bool isChromaQpAdjCoded = getIsChromaQpAdjCoded();
    m_pcEntropyDecoder->decodeCoeff(pcCU, uiAbsPartIdx, uiDepth, bCodeDQP, isChromaQpAdjCoded);
    setIsChromaQpAdjCoded(isChromaQpAdjCoded);
    setdQPFlag(bCodeDQP);

#if JVET_C0024_QTBT
    if (pcCU->getSlice()->getPPS()->getUseDQP())
    {
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiWidth, uiHeight); // set QP
#else
      pcCU->setQPSubParts(getdQPFlag() ? pcCU->getRefQP(uiAbsPartIdx) : pcCU->getCodedQP(), uiAbsPartIdx, uiDepth); // set QP
#endif
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setCodedQP(pcCU->getQP(uiAbsPartIdx));
#endif
    }

    if (pcCU->getSlice()->getUseChromaQpAdj() && !getIsChromaQpAdjCoded())
    {
      pcCU->setChromaQpAdjSubParts(pcCU->getCodedChromaQpAdj(), uiAbsPartIdx, uiDepth); // set QP
    }
#else
    xFinishDecodeCU(pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment);
#endif
}

Void TDecCu::xFinishDecodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool &isLastCtuOfSliceSegment)
{
#if !JVET_C0024_QTBT
  if(  pcCU->getSlice()->getPPS()->getUseDQP())
  {
    pcCU->setQPSubParts( getdQPFlag()?pcCU->getRefQP(uiAbsPartIdx):pcCU->getCodedQP(), uiAbsPartIdx, uiDepth ); // set QP
  }

  if (pcCU->getSlice()->getUseChromaQpAdj() && !getIsChromaQpAdjCoded())
  {
    pcCU->setChromaQpAdjSubParts( pcCU->getCodedChromaQpAdj(), uiAbsPartIdx, uiDepth ); // set QP
  }
#endif

  isLastCtuOfSliceSegment = xDecodeSliceEnd( pcCU, uiAbsPartIdx );
}

#if JVET_C0024_QTBT
Void TDecCu::xDecompressCU( TComDataCU* pCtu, UInt uiAbsPartIdx,  UInt uiDepth, UInt uiWidth, UInt uiHeight )
#else
Void TDecCu::xDecompressCU( TComDataCU* pCtu, UInt uiAbsPartIdx,  UInt uiDepth )
#endif
{
#if !JVET_C0024_QTBT
  TComPic* pcPic = pCtu->getPic();
#endif
  TComSlice * pcSlice = pCtu->getSlice();
  const TComSPS &sps=*(pcSlice->getSPS());
#if !JVET_C0024_QTBT
  Bool bBoundary = false;
#endif

#if JVET_C0024_QTBT
  UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize();
#if BT_FORCE
  UInt uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
#if JVET_C0024_QTBT
  UInt uiRPelX = uiLPelX + uiWidth - 1;
#else
  UInt uiRPelX = uiLPelX + (sps.getMaxCUWidth() >> uiDepth) - 1;
#endif
  UInt uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if JVET_C0024_QTBT
  UInt uiBPelY = uiTPelY + uiHeight - 1;
#else
  UInt uiBPelY = uiTPelY + (sps.getMaxCUHeight() >> uiDepth) - 1;
#endif

#if BARC_BTFORCE
  UInt uiPelXInCTU = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTU + uiWidth)*pCtu->getBARCFlag(pCtu->getTextType());
  UInt deltaY = (uiPelYInCTU + uiHeight)*pCtu->getBARCFlag(pCtu->getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#else
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#endif

#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(pCtu->getTextType()) ? pCtu->getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  const UInt uiDIMDMaxBtDepht = (pCtu->getBoundaryCtuLocation() && pCtu->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (pCtu->getBoundaryCtuLocation() && pCtu->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
#endif

  if (uiCTUSize>>uiDepth == uiWidth && uiWidth==uiHeight)
  {
#endif
#if !BT_FORCE
  UInt uiLPelX   = pCtu->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiRPelX   = uiLPelX + (sps.getCTUSize()>>uiDepth)  - 1;
#else
  UInt uiRPelX   = uiLPelX + (sps.getMaxCUWidth()>>uiDepth)  - 1;
#endif
  UInt uiTPelY   = pCtu->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiBPelY   = uiTPelY + (sps.getCTUSize()>>uiDepth) - 1;
#else
  UInt uiBPelY   = uiTPelY + (sps.getMaxCUHeight()>>uiDepth) - 1;
#endif
#endif
#if BARC_Boundary
  UInt uiPelXInCTU = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

  UInt deltaX =  (uiPelXInCTU+(sps.getCTUSize()>>uiDepth))*pCtu->getBARCFlag(pCtu->getTextType());
  UInt deltaY =  (uiPelYInCTU +(sps.getCTUSize()>>uiDepth))*pCtu->getBARCFlag(pCtu->getTextType());
  if( ( uiRPelX + deltaX >= sps.getPicWidthInLumaSamples() ) || ( uiBPelY + deltaY >= sps.getPicHeightInLumaSamples() ) )
#else
  if( ( uiRPelX >= sps.getPicWidthInLumaSamples() ) || ( uiBPelY >= sps.getPicHeightInLumaSamples() ) )
#endif
  {
#if JVET_C0024_QTBT
#if BT_FORCE
   if (
#if !BT_FORCE_INTRA_ON
     pcSlice->isIntra() || 
#endif
     (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
      assert(uiDepth < pCtu->getDepth(uiAbsPartIdx));
#else
    bBoundary = true;
#endif
  }

#if JVET_C0024_QTBT  
  if( (uiDepth < pCtu->getDepth( uiAbsPartIdx ) )  
#if BT_FORCE
    && (
#if !BT_FORCE_INTRA_ON
      pcSlice->isIntra() || 
#endif
      (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
    )
#else
  if( ( ( uiDepth < pCtu->getDepth( uiAbsPartIdx ) ) && ( uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
#endif
  {
    UInt uiNextDepth = uiDepth + 1;
    UInt uiQNumParts = pCtu->getTotalNumPart() >> (uiNextDepth<<1);
    UInt uiIdx = uiAbsPartIdx;
    for ( UInt uiPartIdx = 0; uiPartIdx < 4; uiPartIdx++ )
    {
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiIdx] ];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiIdx] ];
#if BARC_Boundary
      UInt uiPelXInCTU = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTU = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      UInt deltaY = uiPelYInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      if( ( uiLPelX + deltaX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY + deltaY < sps.getPicHeightInLumaSamples() ) )
#else
      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
#endif
      {
#if JVET_C0024_QTBT
        xDecompressCU(pCtu, uiIdx, uiNextDepth, uiWidth>>1, uiHeight>>1 );
#else
        xDecompressCU(pCtu, uiIdx, uiNextDepth );
#endif
      }
#if JVET_C0024_QTBT
      else
      {
        pCtu->getPic()->addCodedAreaInCTU(uiWidth*uiHeight>>2);
      }
#endif

      uiIdx += uiQNumParts;
    }
    return;
  }

#if JVET_C0024_QTBT
  }
  UInt uiBTDepth = pCtu->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);

#if TT_TREE_DS1
  if (pCtu->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==ST_HORI_BT)
#else
  if (pCtu->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==1)
#endif
  {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
    {
      UInt uiDIMDIntraMode = DC_IDX;
      UInt uiOrgWidth = pCtu->getWidth(uiAbsPartIdx);
      UInt uiOrgHeight = pCtu->getHeight(uiAbsPartIdx);
      pCtu->setWidth(uiAbsPartIdx, uiWidth);
      pCtu->setHeight(uiAbsPartIdx, uiHeight);

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
      UInt uiOrgPDPCIdx = pCtu->getPDPCIdx(uiAbsPartIdx);
      pCtu->setPDPCIdx(uiAbsPartIdx, 0);
#endif
      if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);

#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
      }
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#endif
        if (!pCtu->getSlice()->isIntra())
        {
          assert(pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) == pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
          if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
          {
            assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
          }
        }
        pCtu->setWidth(uiAbsPartIdx, uiOrgWidth);
        pCtu->setHeight(uiAbsPartIdx, uiOrgHeight);

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC 
        pCtu->setPDPCIdx(uiAbsPartIdx, uiOrgPDPCIdx);
#endif
    }
#endif

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] 
        + (uiHeight>>1)/pCtu->getPic()->getMinCUHeight()*pCtu->getPic()->getNumPartInCtuWidth()];
      }
#if BT_FORCE 
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pCtu->getBARCFlag(pCtu->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pCtu->getBARCFlag(pCtu->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDecompressCU(pCtu, uiAbsPartIdx, uiDepth, uiWidth, uiHeight >> 1);
#if BT_FORCE 
      }
      else
      {
        pCtu->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
      }
#endif
    }
    return;
  }
#if TT_TREE_DS1
  else if (pCtu->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==ST_VERT_BT)
#else
  else if (pCtu->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==2)
#endif
  {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
    {
      UInt uiDIMDIntraMode = DC_IDX;
      UInt uiOrgWidth = pCtu->getWidth(uiAbsPartIdx);
      UInt uiOrgHeight = pCtu->getHeight(uiAbsPartIdx);
      pCtu->setWidth(uiAbsPartIdx, uiWidth);
      pCtu->setHeight(uiAbsPartIdx, uiHeight);

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
      UInt uiOrgPDPCIdx = pCtu->getPDPCIdx(uiAbsPartIdx);
      pCtu->setPDPCIdx(uiAbsPartIdx, 0);
#endif
      if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);

#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);

#endif
      }
      if (!pCtu->getSlice()->isIntra())
      {
        assert(pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) == pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
        {
          assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        }
      }

      pCtu->setWidth(uiAbsPartIdx, uiOrgWidth);
      pCtu->setHeight(uiAbsPartIdx, uiOrgHeight);
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC 
      pCtu->setPDPCIdx(uiAbsPartIdx, uiOrgPDPCIdx);
#endif  
    }
#endif

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] 
        + (uiWidth>>1)/pCtu->getPic()->getMinCUWidth()];
      }
#if BT_FORCE
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pCtu->getBARCFlag(pCtu->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pCtu->getBARCFlag(pCtu->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDecompressCU(pCtu, uiAbsPartIdx, uiDepth, uiWidth >> 1, uiHeight);
#if BT_FORCE
      }
      else
      {
        pCtu->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
      }
#endif
    }
    return;
  }
#if TT_DEC_OPT //xDecompressCU
  UInt uiHeightSubPart, uiWidthSubPart, uiAbsPartIdxSubPart;
  UInt uiBTSplitMode = pCtu->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth);
  if (uiBTSplitMode==ST_HORI_TT_L || uiBTSplitMode==ST_HORI_TT_S)
  {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
    {
      UInt uiDIMDIntraMode = DC_IDX;
      UInt uiOrgWidth = pCtu->getWidth(uiAbsPartIdx);
      UInt uiOrgHeight = pCtu->getHeight(uiAbsPartIdx);
      pCtu->setWidth(uiAbsPartIdx, uiWidth);
      pCtu->setHeight(uiAbsPartIdx, uiHeight);

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
      UInt uiOrgPDPCIdx = pCtu->getPDPCIdx(uiAbsPartIdx);
      pCtu->setPDPCIdx(uiAbsPartIdx, 0);
#endif
      if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);
#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);

#endif
      }
      if (!pCtu->getSlice()->isIntra())
      {
        assert(pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) == pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
        {
          assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        }
      }
      pCtu->setWidth(uiAbsPartIdx, uiOrgWidth);
      pCtu->setHeight(uiAbsPartIdx, uiOrgHeight);
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC 
      pCtu->setPDPCIdx(uiAbsPartIdx, uiOrgPDPCIdx);
#endif  
    }

#endif
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + (uiHeight>>2)/pCtu->getPic()->getMinCUHeight()*pCtu->getPic()->getNumPartInCtuWidth()];
        uiHeightSubPart = uiHeight>>1;
      }
      else
      {
        if(uiPartUnitIdx==2)
          uiAbsPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + ((uiHeight>>2)+(uiHeight>>1))/pCtu->getPic()->getMinCUHeight()*pCtu->getPic()->getNumPartInCtuWidth()];
        else
          uiAbsPartIdxSubPart = uiAbsPartIdx;
        uiHeightSubPart = uiHeight>>2;
      }
#if BT_FORCE_TT
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdxSubPart]];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdxSubPart]];
#if BARC_Boundary
      UInt uiPelXInCTU = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTU = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      UInt deltaY = uiPelYInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDecompressCU(pCtu, uiAbsPartIdxSubPart, uiDepth, uiWidth, uiHeightSubPart);
#if BT_FORCE_TT
      }
      else
      {
        pCtu->getPic()->addCodedAreaInCTU(uiWidth*uiHeightSubPart);
      }
#endif
    }
    return;
  }

  if (uiBTSplitMode==ST_VERT_TT_L || uiBTSplitMode==ST_VERT_TT_S)
  {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
    if (isLuma(pCtu->getTextType()) && pCtu->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
    {
      UInt uiDIMDIntraMode = DC_IDX;
      UInt uiOrgWidth = pCtu->getWidth(uiAbsPartIdx);
      UInt uiOrgHeight = pCtu->getHeight(uiAbsPartIdx);
      pCtu->setWidth(uiAbsPartIdx, uiWidth);
      pCtu->setHeight(uiAbsPartIdx, uiHeight);

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
      UInt uiOrgPDPCIdx = pCtu->getPDPCIdx(uiAbsPartIdx);
      pCtu->setPDPCIdx(uiAbsPartIdx, 0);
#endif
      if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);
#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);

#endif
      }
      if (!pCtu->getSlice()->isIntra())
      {
        assert(pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) == pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        if (pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx) && pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
        {
          assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
        }
      }
      pCtu->setWidth(uiAbsPartIdx, uiOrgWidth);
      pCtu->setHeight(uiAbsPartIdx, uiOrgHeight);
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC 
      pCtu->setPDPCIdx(uiAbsPartIdx, uiOrgPDPCIdx);
#endif  
    }
#endif

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + (uiWidth>>2)/pCtu->getPic()->getMinCUWidth()];
        uiWidthSubPart = uiWidth>>1;
      }
      else
      {
        if(uiPartUnitIdx==2)
          uiAbsPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + ((uiWidth>>2)+(uiWidth>>1))/pCtu->getPic()->getMinCUWidth()];
        else
          uiAbsPartIdxSubPart = uiAbsPartIdx;
        uiWidthSubPart = uiWidth>>2;
      }
#if BT_FORCE_TT
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdxSubPart]];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdxSubPart]];
#if BARC_Boundary
      UInt uiPelXInCTU = uiLPelX - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTU = uiTPelY - pCtu->getPic()->getCtu(pCtu->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      UInt deltaY = uiPelYInCTU*pCtu->getBARCFlag(pCtu->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDecompressCU(pCtu, uiAbsPartIdxSubPart, uiDepth, uiWidthSubPart, uiHeight);
#if BT_FORCE_TT
      }
      else
      {
        pCtu->getPic()->addCodedAreaInCTU(uiWidthSubPart*uiHeight);
      }
#endif
    }
    return;
  }
#endif

  UInt uiWidthIdx = g_aucConvertToBit[uiWidth];
  UInt uiHeightIdx = g_aucConvertToBit[uiHeight];

#if TT_TREE_AFFINE
  pCtu->setCUAbsZIdxInCtuSubParts(uiAbsPartIdx);
#endif

#if DIMD_INTRA_PRED
  if (isLuma(pCtu->getTextType()))  // luma of I-slice, and luma and chroma of P/B-slice
  {
#if BT_FORCE_DIMD_HARM
    if (uiBTDepth <= uiDIMDMaxBtDepht || !pCtu->getSlice()->isIntra()) // BTdepth<=1 || P,B slice
#else
    if (uiBTDepth <= DIMD_MAX_BT_DEPTH_SIG_AI || !pCtu->getSlice()->isIntra()) // BTdepth<=1 || P,B slice
#endif
    {
      UChar ucDIMDFlag = pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
      UInt uiDIMDIntraMode = DC_IDX;

      if (ucDIMDFlag)
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);
#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#endif
      }
      if (!pCtu->getSlice()->isIntra())
      {
        assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
      }
    }
    else
    {
      UInt uiDIMDIntraMode = DC_IDX;
      UChar ucDIMDFlag = pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
      if (ucDIMDFlag && !pCtu->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      {
        uiDIMDIntraMode = m_pcPrediction->deriveNeighborIntraDirs(pCtu, uiAbsPartIdx);
#if DIMD_NUM_INTRA_DIR_INC
        pCtu->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiDIMDIntraMode), uiAbsPartIdx, uiDepth);
#endif
#else
        pCtu->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiDIMDIntraMode, uiAbsPartIdx, uiDepth);
#endif
      }

      if (!pCtu->getSlice()->isIntra())
      {
        assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
      }
    }
  }
  else  // chroma of I-slice
  {
    assert(!pCtu->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
  }
#endif

  // Residual reconstruction
  m_pppcYuvResi[uiWidthIdx][uiHeightIdx]->clear();

  m_pppcCU[uiWidthIdx][uiHeightIdx]->copySubCU( pCtu, uiAbsPartIdx );
#if  TRACE_FLAGS 
  if (pCtu->getTextType() == CHANNEL_TYPE_LUMA)
	  gNumCu++;
#endif

  switch( m_pppcCU[uiWidthIdx][uiHeightIdx]->getPredictionMode(0) )
  {
    case MODE_INTER:
#if VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
      xDeriveCUMV( pCtu , uiAbsPartIdx , uiDepth );
#if TRACE_FLAGS && JVET_C0024_QTBT
	  DUMP_CU_FLAGS_INFO("decoder", pCtu, uiAbsPartIdx, uiAbsPartIdx, uiAbsPartIdx);
#endif
#endif
      xReconInter( m_pppcCU[uiWidthIdx][uiHeightIdx], uiDepth );
      break;
    case MODE_INTRA:
      xReconIntraQT( m_pppcCU[uiWidthIdx][uiHeightIdx], uiDepth );
      break;
    default:
      assert(0);
      break;
  }
#else
  // Residual reconstruction
  m_ppcYuvResi[uiDepth]->clear();

  m_ppcCU[uiDepth]->copySubCU( pCtu, uiAbsPartIdx );

  switch( m_ppcCU[uiDepth]->getPredictionMode(0) )
  {
    case MODE_INTER:
#if TRACE_FLAGS
	  DUMP_PU_FLAGS_INTER("decoder", pCtu, uiAbsPartIdx, -1);
#endif
#if VCEG_AZ07_FRUC_MERGE
      xDeriveCUMV( pCtu , uiAbsPartIdx , uiDepth );
#endif
      xReconInter( m_ppcCU[uiDepth], uiDepth );
      break;
    case MODE_INTRA:
      xReconIntraQT( m_ppcCU[uiDepth], uiDepth );
      break;
    default:
      assert(0);
      break;
  }
#endif

#if DEBUG_STRING
  const PredMode predMode=m_ppcCU[uiDepth]->getPredictionMode(0);
  if (DebugOptionList::DebugString_Structure.getInt()&DebugStringGetPredModeMask(predMode))
  {
    PartSize eSize=m_ppcCU[uiDepth]->getPartitionSize(0);
    std::ostream &ss(std::cout);

    ss <<"###: " << (predMode==MODE_INTRA?"Intra   ":"Inter   ") << partSizeToString[eSize] << " CU at " << m_ppcCU[uiDepth]->getCUPelX() << ", " << m_ppcCU[uiDepth]->getCUPelY() << " width=" << UInt(m_ppcCU[uiDepth]->getWidth(0)) << std::endl;
  }
#endif

#if JVET_C0024_QTBT
  if ( m_pppcCU[uiWidthIdx][uiHeightIdx]->isLosslessCoded(0) && (m_pppcCU[uiWidthIdx][uiHeightIdx]->getIPCMFlag(0) == false))
  {
    xFillPCMBuffer(m_pppcCU[uiWidthIdx][uiHeightIdx], uiDepth);
  }

  UInt uiBlkX = g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ] >> MIN_CU_LOG2;
  UInt uiBlkY = g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ] >> MIN_CU_LOG2;
  pCtu->getPic()->setCodedBlkInCTU(true, uiBlkX, uiBlkY, uiWidth>> MIN_CU_LOG2, uiHeight>> MIN_CU_LOG2);
  pCtu->getPic()->addCodedAreaInCTU(uiWidth*uiHeight);
#if TRACE_YUV
  if (!pCtu->getSlice()->isIntra() || pCtu->getTextType() == CHANNEL_TYPE_LUMA) // add me
  {

    TRACE_CU_PART(m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcCU[uiWidthIdx][uiHeightIdx], "t_rec_dec.yuv");
    TRACE_CU_PART(m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcCU[uiWidthIdx][uiHeightIdx], "t_rec_dec_no_mark.yuv");
    TRACE_CU_PART_DIFF(m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcYuvResi[uiWidthIdx][uiHeightIdx], m_pppcCU[uiWidthIdx][uiHeightIdx], "t_predicted_dec.yuv");
    MARK_PARTIONING(m_pppcCU[uiWidthIdx][uiHeightIdx], uiAbsPartIdx, uiWidth, uiHeight, "t_rec_dec.yuv");
    MARK_PARTIONING(m_pppcCU[uiWidthIdx][uiHeightIdx], uiAbsPartIdx, uiWidth, uiHeight, "t_predicted_dec.yuv");
  }
#endif

#else
  if ( m_ppcCU[uiDepth]->isLosslessCoded(0) && (m_ppcCU[uiDepth]->getIPCMFlag(0) == false))
  {
    xFillPCMBuffer(m_ppcCU[uiDepth], uiDepth);
  }

  xCopyToPic( m_ppcCU[uiDepth], pcPic, uiAbsPartIdx, uiDepth );
#if TRACE_YUV
  TRACE_CU_PART(m_ppcYuvReco[uiDepth], m_ppcCU[uiDepth], "t_rec_dec.yuv");
  TRACE_CU_PART(m_ppcYuvReco[uiDepth], m_ppcCU[uiDepth], "t_rec_no_mark.yuv");
  TRACE_CU_PART_DIFF(m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], m_ppcCU[uiDepth], "t_predicted_dec.yuv");
  MARK_PARTIONING(m_ppcCU[uiDepth], 0, m_ppcCU[uiDepth]->getWidth(uiAbsPartIdx), m_ppcCU[uiDepth]->getHeight(uiAbsPartIdx), "t_rec_dec.yuv");
  MARK_PARTIONING(m_ppcCU[uiDepth], 0, m_ppcCU[uiDepth]->getWidth(uiAbsPartIdx), m_ppcCU[uiDepth]->getHeight(uiAbsPartIdx), "t_predicted_dec.yuv");
#endif
#endif
}

Void TDecCu::xReconInter( TComDataCU* pcCU, UInt uiDepth )
{

#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[pcCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[pcCU->getHeight(0)];
  m_pcPrediction->motionCompensation( pcCU, m_pppcYuvReco[uiWidthIdx][uiHeightIdx] );
#if COM16_C806_OBMC
  m_pcPrediction->subBlockOBMC( pcCU, 0, m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv1[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv2[uiWidthIdx][uiHeightIdx] );
#endif
#else
  // inter prediction
  m_pcPrediction->motionCompensation( pcCU, m_ppcYuvReco[uiDepth] );
#if COM16_C806_OBMC
  m_pcPrediction->subBlockOBMC( pcCU, 0, m_ppcYuvReco[uiDepth], m_ppcTmpYuv1[uiDepth], m_ppcTmpYuv2[uiDepth] );
#endif
#endif
#if DEBUG_STRING
  const Int debugPredModeMask=DebugStringGetPredModeMask(MODE_INTER);
  if (DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-pred: ", *(m_ppcYuvReco[uiDepth]));
  }
#endif

  // inter recon
  xDecodeInterTexture( pcCU, uiDepth );

#if DEBUG_STRING
  if (DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-resi: ", *(m_ppcYuvResi[uiDepth]));
  }
#endif

  // clip for only non-zero cbp case
  if  ( pcCU->getQtRootCbf( 0) )
  {
#if JVET_C0024_QTBT
    m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->addClip( m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcYuvResi[uiWidthIdx][uiHeightIdx], 0, pcCU->getWidth( 0 ), pcCU->getHeight(0), pcCU->getSlice()->getSPS()->getBitDepths() );
#else
    m_ppcYuvReco[uiDepth]->addClip( m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], 0, pcCU->getWidth( 0 ), pcCU->getSlice()->getSPS()->getBitDepths() );
#endif
  }
#if !JVET_C0024_QTBT
  else
  {
#if JVET_D0033_ADAPTIVE_CLIPPING // decoder, inter no res
      if (g_ClipParam.isActive)
      {
          m_ppcYuvReco[uiDepth]->clipPartToPartYuv( m_ppcYuvReco[uiDepth], 0, pcCU->getWidth(0), pcCU->getHeight(0));
      } 
      else 
      {
    m_ppcYuvReco[uiDepth]->copyPartToPartYuv( m_ppcYuvReco[uiDepth],0, pcCU->getWidth( 0 ),pcCU->getHeight( 0 ));
  }
#else
    m_ppcYuvReco[uiDepth]->copyPartToPartYuv( m_ppcYuvReco[uiDepth],0, pcCU->getWidth( 0 ),pcCU->getHeight( 0 ));
#endif
  }
#endif
#if DEBUG_STRING
  if (DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-reco: ", *(m_ppcYuvReco[uiDepth]));
  }
#endif

#if JVET_C0024_QTBT
#if JVET_D0033_ADAPTIVE_CLIPPING
  m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->clipToPicYuv(pcCU->getPic()->getPicYuvRec(), pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
#else
  m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->copyToPicYuv(pcCU->getPic()->getPicYuvRec(), pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
#endif
#endif
}


Void
TDecCu::xIntraRecBlk(       TComYuv*    pcRecoYuv,
                            TComYuv*    pcPredYuv,
                            TComYuv*    pcResiYuv,
                      const ComponentID compID,
                            TComTU     &rTu)
{
  if (!rTu.ProcessComponentSection(compID))
  {
    return;
  }
  const Bool       bIsLuma = isLuma(compID);


  TComDataCU *pcCU = rTu.getCU();
#if !JVET_D0033_ADAPTIVE_CLIPPING
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());
#endif
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

  const TComRectangle &tuRect  =rTu.getRect(compID);
  const UInt uiWidth           = tuRect.width;
  const UInt uiHeight          = tuRect.height;
  const UInt uiStride          = pcRecoYuv->getStride (compID);
        Pel* piPred            = pcPredYuv->getAddr( compID, uiAbsPartIdx );
  const ChromaFormat chFmt     = rTu.GetChromaFormat();

#if !JVET_C0024_QTBT
  if (uiWidth != uiHeight)
  {
    //------------------------------------------------

    //split at current level if dividing into square sub-TUs

    TComTURecurse subTURecurse(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

    //recurse further
    do
    {
      xIntraRecBlk(pcRecoYuv, pcPredYuv, pcResiYuv, compID, subTURecurse);
    } while (subTURecurse.nextSection(rTu));

    //------------------------------------------------

    return;
  }
#endif
#if TRACE_FLAGS
#if !JVET_C0024_QTBT
  if (COMPONENT_Y == compID)
  {
	  DUMP_TU_FLAGS_INFO("decoder", pcCU, uiAbsPartIdx, pcCU->getZorderIdxInCtu() + uiAbsPartIdx, rTu);
  }
#else
  if (COMPONENT_Y == compID)
  {
	  DUMP_CU_FLAGS_INFO("decoder", pcCU, uiAbsPartIdx, pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
  }
#endif
#endif
  const UInt uiChPredMode  = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx );
#if JVET_C0024_QTBT
#if JVET_E0062_MULTI_DMS
  const UInt uiChCodedMode = uiChPredMode; 
#else
  const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) : uiChPredMode;
#endif
#else
  const UInt partsPerMinCU = 1<<(2*(sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
  const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
#endif
  const UInt uiChFinalMode = ((chFmt == CHROMA_422)       && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;

#if RRC_INTRA_SHARPENING
  Bool doSharpening = pcCU->getIntraSharpFlag(uiAbsPartIdx);
#if !INTRA_SHARP_CHROMA
  doSharpening = (compID == COMPONENT_Y) ? doSharpening : 0;
#elif JVET_C0024_QTBT
  if (isChroma(pcCU->getTextType()))
  {
    UInt absPartIdx = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
    absPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + (pcCU->getHeight(uiAbsPartIdx) / pcCU->getPic()->getMinCUHeight()) / 2 * pcCU->getPic()->getNumPartInCtuWidth() + (pcCU->getWidth(uiAbsPartIdx) / pcCU->getPic()->getMinCUWidth()) / 2];
    doSharpening = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->isIntra(absPartIdx) ? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraSharpFlag(absPartIdx) : 0;
    if (uiChFinalMode <= DC_IDX || uiChFinalMode >= LM_CHROMA_IDX)
      doSharpening = 0;
  }
#endif
#endif

  //===== init availability pattern =====
#if !RRC_INTRA_SHARPENING
  const 
#endif
        Bool bUseFilteredPredictions=TComPrediction::filteringIntraReferenceSamples(compID, 
                                                                                    uiChFinalMode, 
                                                                                    uiWidth, 
                                                                                    uiHeight, 
                                                                                    chFmt, 
                                                                                    pcCU->getSlice()->getSPS()->getSpsRangeExtension().getIntraSmoothingDisabledFlag()
#if DIMD_NUM_INTRA_DIR_INC
			, pcCU, uiAbsPartIdx
#endif
                                                                                   );
#if DEBUG_STRING
  std::ostream &ss(std::cout);
#endif

  DEBUG_STRING_NEW(sTemp)
#if BARC_CORE
  Bool isCTUBoundary = 0;
  if (pcCU->getBARCFlag(toChannelType(compID)))
  {
#if JVET_G0104_PLANAR_PDPC
		UInt leftBoundary = (pcCU->getCUPelX() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
		UInt aboveBoundary = (pcCU->getCUPelY() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
#else
    UInt leftBoundary = (pcCU->getCUPelX() % sps.getCTUSize()) == 0 ? 1 : 0;
    UInt aboveBoundary = (pcCU->getCUPelY() % sps.getCTUSize()) == 0 ? 1 : 0;
#endif
    isCTUBoundary = leftBoundary + aboveBoundary > 0 ? 1 : 0;
  }
#endif
  m_pcPrediction->initIntraPatternChType( rTu, compID, bUseFilteredPredictions  
#if RRC_INTRA_SHARPENING
                                        , doSharpening
#endif
#if BARC_CORE
      , isCTUBoundary
#endif
                                          DEBUG_STRING_PASS_INTO(sTemp) );


  //===== get prediction signal =====

#if COM16_C806_LMCHROMA
  if( uiChFinalMode == LM_CHROMA_IDX 
#if JVET_E0077_ENHANCED_LM
      || IsLMMode(uiChFinalMode)
#endif
      )
  {
    m_pcPrediction->getLumaRecPixels( rTu, uiWidth, uiHeight );
#if JVET_E0077_ENHANCED_LM
    Int iLMType = uiChFinalMode;
    m_pcPrediction->predLMIntraChroma(rTu, compID, piPred, uiStride, uiWidth, uiHeight, iLMType);
#else
    m_pcPrediction->predLMIntraChroma( rTu, compID, piPred, uiStride, uiWidth, uiHeight );
#endif
  }
  else
  {
#endif 

#if RRC_INTRA_SHARPENING
  bUseFilteredPredictions |= doSharpening;
#endif
#if RRC_DWDIP && !RRC_DWDIP_2Nx2N_CHECK_DISABLE
    if (pcCU->getDWDIPIdx(uiAbsPartIdx))
      assert(pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N);
#endif
    m_pcPrediction->predIntraAng( compID,   
                                  uiChFinalMode, 
                                 0 /* Decoder does not have an original image */, 0, piPred, uiStride, rTu, bUseFilteredPredictions
                                );

#if COM16_C806_LMCHROMA
    if( compID == COMPONENT_Cr && pcCU->getSlice()->getSPS()->getUseLMChroma() )
    { 
      m_pcPrediction->addCrossColorResi( rTu, compID, piPred, uiStride, uiWidth, uiHeight, pcResiYuv->getAddr( COMPONENT_Cb, uiAbsPartIdx ), pcResiYuv->getStride(COMPONENT_Cb) );
    }
  }
#endif

#if DEBUG_STRING
  ss << sTemp;
#endif

  //===== inverse transform =====
  Pel*      piResi            = pcResiYuv->getAddr( compID, uiAbsPartIdx );
  TCoeff*   pcCoeff           = pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID);//( uiNumCoeffInc * uiAbsPartIdx );

  const QpParam cQP(*pcCU, compID);


  DEBUG_STRING_NEW(sDebug);
#if DEBUG_STRING
  const Int debugPredModeMask=DebugStringGetPredModeMask(MODE_INTRA);
  std::string *psDebug=(DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask) ? &sDebug : 0;
#endif

  if (pcCU->getCbf(uiAbsPartIdx, compID, rTu.GetTransformDepthRel()) != 0)
  {
#if PVC_AQS
    m_pcTrQuant->setEstQSFactor(piPred, uiStride, uiWidth, uiHeight, pcCU->getSlice()->getSliceType(), pcCU, uiAbsPartIdx, COMPONENT_Y);
#endif
    m_pcTrQuant->invTransformNxN( rTu, compID, piResi, uiStride, pcCoeff, cQP DEBUG_STRING_PASS_INTO(psDebug) );
  }
  else
  {
    for (UInt y = 0; y < uiHeight; y++)
    {
      for (UInt x = 0; x < uiWidth; x++)
      {
        piResi[(y * uiStride) + x] = 0;
      }
    }
  }

#if DEBUG_STRING
  if (psDebug)
  {
    ss << (*psDebug);
  }
#endif

  //===== reconstruction =====
  const UInt uiRecIPredStride  = pcCU->getPic()->getPicYuvRec()->getStride(compID);

  const Bool useCrossComponentPrediction = isChroma(compID) && (pcCU->getCrossComponentPredictionAlpha(uiAbsPartIdx, compID) != 0);
  const Pel* pResiLuma  = pcResiYuv->getAddr( COMPONENT_Y, uiAbsPartIdx );
  const Int  strideLuma = pcResiYuv->getStride( COMPONENT_Y );

        Pel* pPred      = piPred;
        Pel* pResi      = piResi;
        Pel* pReco      = pcRecoYuv->getAddr( compID, uiAbsPartIdx );
        Pel* pRecIPred  = pcCU->getPic()->getPicYuvRec()->getAddr( compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx );
#if BARC
        Pel* pDscaleRec = pcCU->getPic()->getPicYuvDScaleRef()->getAddrDscale( compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx );
#endif

#if DEBUG_STRING
  const Bool bDebugPred=((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  const Bool bDebugResi=((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  const Bool bDebugReco=((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  if (bDebugPred || bDebugResi || bDebugReco)
  {
    ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << std::endl;
  }
#endif

#if !JVET_D0033_ADAPTIVE_CLIPPING
  const Int clipbd = sps.getBitDepth(toChannelType(compID));
#endif
#if O0043_BEST_EFFORT_DECODING
  const Int bitDepthDelta = sps.getStreamBitDepth(toChannelType(compID)) - clipbd;
#endif

  if( useCrossComponentPrediction )
  {
    TComTrQuant::crossComponentPrediction( rTu, compID, pResiLuma, piResi, piResi, uiWidth, uiHeight, strideLuma, uiStride, uiStride, true );
  }

  for( UInt uiY = 0; uiY < uiHeight; uiY++ )
  {
#if DEBUG_STRING
    if (bDebugPred || bDebugResi || bDebugReco)
    {
      ss << "###: ";
    }

    if (bDebugPred)
    {
      ss << " - pred: ";
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        ss << pPred[ uiX ] << ", ";
      }
    }
    if (bDebugResi)
    {
      ss << " - resi: ";
    }
#endif

    for( UInt uiX = 0; uiX < uiWidth; uiX++ )
    {
#if DEBUG_STRING
      if (bDebugResi)
      {
        ss << pResi[ uiX ] << ", ";
      }
#endif
#if O0043_BEST_EFFORT_DECODING
#if JVET_D0033_ADAPTIVE_CLIPPING
      pReco    [ uiX ] = ClipA( rightShiftEvenRounding<Pel>(pPred[ uiX ] + pResi[ uiX ], bitDepthDelta),  compID);
#else
      pReco    [ uiX ] = ClipBD( rightShiftEvenRounding<Pel>(pPred[ uiX ] + pResi[ uiX ], bitDepthDelta), clipbd );
#endif
#else

#if JVET_D0033_ADAPTIVE_CLIPPING // decoder intra rec
      pReco    [ uiX ] = ClipA(pPred[ uiX ] + pResi[ uiX ], compID);
#else
      pReco    [ uiX ] = ClipBD( pPred[ uiX ] + pResi[ uiX ], clipbd );
#endif
#endif
      pRecIPred[ uiX ] = pReco[ uiX ];
#if BARC
      pDscaleRec[ uiX ] = pReco[ uiX ];
#endif
    }
#if DEBUG_STRING
    if (bDebugReco)
    {
      ss << " - reco: ";
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        ss << pReco[ uiX ] << ", ";
      }
    }

    if (bDebugPred || bDebugResi || bDebugReco)
    {
      ss << "\n";
    }
#endif
    pPred     += uiStride;
    pResi     += uiStride;
    pReco     += uiStride;
    pRecIPred += uiRecIPredStride;
#if BARC
    pDscaleRec += uiRecIPredStride;
#endif
  }
  #if JVET_F0096_BILATERAL_FILTER
	if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
	{
		Pel* piReco = pcRecoYuv->getAddr(compID, uiAbsPartIdx);
		Pel* piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
#if BARC
		Pel* piDscaleRec = pcCU->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
#endif

		if (isLuma(compID))
		{
			if ((pcCU->getCbf(uiAbsPartIdx, compID, rTu.GetTransformDepthRel()) != 0) && (pcCU->getQP(COMPONENT_Y) > 17))
			{
				TComBilateralFilter::instance()->bilateralFilterIntra(pcCU, uiWidth, uiHeight, piReco, uiStride, pcCU->getQP(COMPONENT_Y));
				
				twoDimBufCopier
					(
						(UChar *)piRecIPred,
						uiRecIPredStride * sizeof(Pel),
						(UChar *)piReco,
						uiStride * sizeof(Pel),
						uiWidth * sizeof(Short),
						uiHeight
						);
#if BARC
				twoDimBufCopier
					(
						(UChar *)piDscaleRec,
						uiRecIPredStride * sizeof(Pel),
						(UChar *)piReco,
						uiStride * sizeof(Pel),
						uiWidth * sizeof(Short),
						uiHeight
						);

#endif
				
//				for (UInt uiY = 0; uiY < uiHeight; uiY++)
//				{
//					memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
//#if BARC
//					memcpy(piDscaleRec + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
//#endif
////					uiY++;
////					memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
////#if BARC
////					memcpy(piDscaleRec + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
////#endif
//				}
			}
		}
	}
#endif
}

#if VCEG_AZ08_INTRA_KLT
Void
TDecCu::xIntraRecBlkTM( TComYuv*    pcRecoYuv,
                      TComYuv*    pcPredYuv,
                      TComYuv*    pcResiYuv,
                      const ComponentID compID,
                      TComTU     &rTu,
                      Int tmpred0_tmpredklt1_ori2
                      )
{
    if (!rTu.ProcessComponentSection(compID))
    {
        return;
    }

    TComDataCU *pcCU = rTu.getCU();
#if !JVET_D0033_ADAPTIVE_CLIPPING
    const TComSPS &sps = *(pcCU->getSlice()->getSPS());
#endif
    const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();

    const TComRectangle &tuRect = rTu.getRect(compID);
    const UInt uiWidth = tuRect.width;
    const UInt uiHeight = tuRect.height;
    const UInt uiStride = pcRecoYuv->getStride(compID);
    Pel* piPred = pcPredYuv->getAddr(compID, uiAbsPartIdx);
    assert(uiWidth == uiHeight);

    Bool useKLT = false;
    UInt uiBlkSize = uiWidth;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiTempSize = g_uiDepth2IntraTempSize[uiTarDepth];
    m_pcTrQuant->getTargetTemplate(pcCU, uiAbsPartIdx, uiBlkSize, uiTempSize);
    m_pcTrQuant->candidateSearchIntra(pcCU, uiAbsPartIdx, uiBlkSize, uiTempSize);
    Int foundCandiNum;
    Bool bSuccessful = m_pcTrQuant->generateTMPrediction(piPred, uiStride, uiBlkSize, uiTempSize, foundCandiNum);
    if (1 == tmpred0_tmpredklt1_ori2 && bSuccessful)
    {
        useKLT = m_pcTrQuant->calcKLTIntra(piPred, uiStride, uiBlkSize);
    }
    assert(foundCandiNum >= 1);

    //===== inverse transform =====
    Pel*      piResi = pcResiYuv->getAddr(compID, uiAbsPartIdx);
    TCoeff*   pcCoeff = pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID);

    const QpParam cQP(*pcCU, compID);


    DEBUG_STRING_NEW(sDebug);
#if DEBUG_STRING
    const Int debugPredModeMask = DebugStringGetPredModeMask(MODE_INTRA);
    std::string *psDebug = (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask) ? &sDebug : 0;
#endif

    if (pcCU->getCbf(uiAbsPartIdx, compID, rTu.GetTransformDepthRel()) != 0)
    {
        const UInt *scan = NULL;
        if (useKLT)
        {
            TUEntropyCodingParameters codingParameters;
            getTUEntropyCodingParameters(codingParameters, rTu, compID);
            scan = codingParameters.scan;
            recoverOrderCoeff(pcCoeff, scan, uiWidth, uiHeight);
        }
        m_pcTrQuant->invTransformNxN(rTu, compID, piResi, uiStride, pcCoeff, cQP, useKLT DEBUG_STRING_PASS_INTO(psDebug));
        if (useKLT)
        {
            reOrderCoeff(pcCoeff, scan, uiWidth, uiHeight);
        }
    }
    else
    {
        for (UInt y = 0; y < uiHeight; y++)
        {
            for (UInt x = 0; x < uiWidth; x++)
            {
                piResi[(y * uiStride) + x] = 0;
            }
        }
    }

    //===== reconstruction =====
    const UInt uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);

    const Bool useCrossComponentPrediction = isChroma(compID) && (pcCU->getCrossComponentPredictionAlpha(uiAbsPartIdx, compID) != 0);
    assert(useCrossComponentPrediction == false);
    Pel* pPred = piPred;
    Pel* pResi = piResi;
    Pel* pReco = pcRecoYuv->getAddr(compID, uiAbsPartIdx);
    Pel* pRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);

#if DEBUG_STRING
    const Bool bDebugPred = ((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
    const Bool bDebugResi = ((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
    const Bool bDebugReco = ((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
    if (bDebugPred || bDebugResi || bDebugReco)
    {
        ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << std::endl;
    }
#endif

#if !JVET_D0033_ADAPTIVE_CLIPPING
    const Int clipbd = sps.getBitDepth(toChannelType(compID));
#endif
#if O0043_BEST_EFFORT_DECODING
    const Int bitDepthDelta = sps.getStreamBitDepth(toChannelType(compID)) - clipbd;
#endif

    for (UInt uiY = 0; uiY < uiHeight; uiY++)
    {
#if DEBUG_STRING
        if (bDebugPred || bDebugResi || bDebugReco)
        {
            ss << "###: ";
        }

        if (bDebugPred)
        {
            ss << " - pred: ";
            for (UInt uiX = 0; uiX < uiWidth; uiX++)
            {
                ss << pPred[uiX] << ", ";
            }
        }
        if (bDebugResi)
        {
            ss << " - resi: ";
        }
#endif

        for (UInt uiX = 0; uiX < uiWidth; uiX++)
        {
#if DEBUG_STRING
            if (bDebugResi)
            {
                ss << pResi[uiX] << ", ";
            }
#endif
#if JVET_D0033_ADAPTIVE_CLIPPING
#if O0043_BEST_EFFORT_DECODING

            pReco[uiX] = ClipA(rightShiftEvenRounding<Pel>(pPred[uiX] + pResi[uiX], bitDepthDelta), compID);
#else
            pReco[uiX] = ClipA(pPred[uiX] + pResi[uiX], compID);
#endif
#else
#if O0043_BEST_EFFORT_DECODING

            pReco[uiX] = ClipBD(rightShiftEvenRounding<Pel>(pPred[uiX] + pResi[uiX], bitDepthDelta), clipbd);
#else
            pReco[uiX] = ClipBD(pPred[uiX] + pResi[uiX], clipbd);
#endif
#endif
            pRecIPred[uiX] = pReco[uiX];
        }
#if DEBUG_STRING
        if (bDebugReco)
        {
            ss << " - reco: ";
            for (UInt uiX = 0; uiX < uiWidth; uiX++)
            {
                ss << pReco[uiX] << ", ";
            }
        }

        if (bDebugPred || bDebugResi || bDebugReco)
        {
            ss << "\n";
        }
#endif
        pPred += uiStride;
        pResi += uiStride;
        pReco += uiStride;
        pRecIPred += uiRecIPredStride;
    }
#if JVET_F0096_BILATERAL_FILTER
		if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
		{
			Pel* piReco = pcRecoYuv->getAddr(compID, uiAbsPartIdx);
			Pel* piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
			if (isLuma(compID))
			{
				if ((pcCU->getCbf(uiAbsPartIdx, compID, rTu.GetTransformDepthRel()) != 0) && (pcCU->getQP(COMPONENT_Y) > 17))
				{
					TComBilateralFilter::instance()->bilateralFilterIntra(pcCU, uiWidth, uiHeight, piReco, uiStride, pcCU->getQP(COMPONENT_Y));
					
					twoDimBufCopier
						(
							(UChar *)piRecIPred,
							uiRecIPredStride * sizeof(Pel),
							(UChar *)piReco,
							uiStride * sizeof(Pel),
							uiWidth * sizeof(Short),
							uiHeight
							);
					
					//for (UInt uiY = 0; uiY < uiHeight; uiY++)
					//{
					//	memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
					//	//uiY++;
					//	//memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
					//}
				}
			}
		}
#endif
}
#endif


Void
TDecCu::xReconIntraQT( TComDataCU* pcCU, UInt uiDepth )
{
  if (pcCU->getIPCMFlag(0))
  {
    xReconPCM( pcCU, uiDepth );
    return;
  }
#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[pcCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[pcCU->getHeight(0)];

  const UInt numChType = (pcCU->getPic()->getChromaFormat()==CHROMA_400 
    || (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra())) ? 1 : 2;
  UInt stChType = isChroma(pcCU->getTextType()) ? CHANNEL_TYPE_CHROMA: CHANNEL_TYPE_LUMA;
  for (UInt chType=stChType; chType<numChType; chType++)
#else
  const UInt numChType = pcCU->getPic()->getChromaFormat()!=CHROMA_400 ? 2 : 1;
  for (UInt chType=CHANNEL_TYPE_LUMA; chType<numChType; chType++)
#endif
  {
    const ChannelType chanType=ChannelType(chType);
#if JVET_C0024_QTBT
    const UInt uiInitTrDepth = 0;
#else
    const Bool NxNPUHas4Parts = ::isChroma(chanType) ? enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) : true;
    const UInt uiInitTrDepth = ( pcCU->getPartitionSize(0) != SIZE_2Nx2N && NxNPUHas4Parts ? 1 : 0 );
#endif

    TComTURecurse tuRecurseCU(pcCU, 0);
    TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth==0)?TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

    do
    {
#if JVET_C0024_QTBT
      xIntraRecQT( m_pppcYuvReco[uiWidthIdx][uiHeightIdx], m_pppcYuvReco[uiWidthIdx][uiHeightIdx]
      , m_pppcYuvResi[uiWidthIdx][uiHeightIdx], chanType, tuRecurseWithPU );
#else
      xIntraRecQT( m_ppcYuvReco[uiDepth], m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], chanType, tuRecurseWithPU );
#endif
    } while (tuRecurseWithPU.nextSection(tuRecurseCU));
  }
}



/** Function for deriving reconstructed PU/CU chroma samples with QTree structure
 * \param pcRecoYuv pointer to reconstructed sample arrays
 * \param pcPredYuv pointer to prediction sample arrays
 * \param pcResiYuv pointer to residue sample arrays
 * \param chType    texture channel type (luma/chroma)
 * \param rTu       reference to transform data
 *
 \ This function derives reconstructed PU/CU chroma samples with QTree recursive structure
 */

Void
TDecCu::xIntraRecQT(TComYuv*    pcRecoYuv,
                    TComYuv*    pcPredYuv,
                    TComYuv*    pcResiYuv,
                    const ChannelType chType,
                    TComTU     &rTu)
{
  UInt uiTrDepth    = rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
  assert(uiTrDepth==0);
#else
  TComDataCU *pcCU  = rTu.getCU();
  UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if( uiTrMode == uiTrDepth )
  {
#endif
    if (isLuma(chType))
    {
#if VCEG_AZ08_INTRA_KLT
#if VCEG_AZ08_USE_KLT
#if JVET_C0024_QTBT
      TComDataCU *pcCU  = rTu.getCU();
      UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#endif
        if (pcCU->getSlice()->getSPS()->getUseIntraKLT() && (Bool)pcCU->getKLTFlag(uiAbsPartIdx, COMPONENT_Y))
#else
        Bool bTMFlag = (Bool)pcCU->getKLTFlag(uiAbsPartIdx, COMPONENT_Y);
        if (bTMFlag)
#endif
        {
            xIntraRecBlkTM(pcRecoYuv, pcPredYuv, pcResiYuv, COMPONENT_Y, rTu, TMPRED0_TMPREDKLT1_ORI2);
        }
        else
        {
      xIntraRecBlk( pcRecoYuv, pcPredYuv, pcResiYuv, COMPONENT_Y,  rTu );
    }
#else
      xIntraRecBlk( pcRecoYuv, pcPredYuv, pcResiYuv, COMPONENT_Y,  rTu );
#endif
    }
    else
    {
      const UInt numValidComp=getNumberValidComponents(rTu.GetChromaFormat());
      for(UInt compID=COMPONENT_Cb; compID<numValidComp; compID++)
      {
        xIntraRecBlk( pcRecoYuv, pcPredYuv, pcResiYuv, ComponentID(compID), rTu );
      }
    }
#if !JVET_C0024_QTBT
  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xIntraRecQT( pcRecoYuv, pcPredYuv, pcResiYuv, chType, tuRecurseChild );
    } while (tuRecurseChild.nextSection(rTu));
  }
#endif
}

#if !JVET_C0024_QTBT
Void TDecCu::xCopyToPic( TComDataCU* pcCU, TComPic* pcPic, UInt uiZorderIdx, UInt uiDepth )
{
  UInt uiCtuRsAddr = pcCU->getCtuRsAddr();

  m_ppcYuvReco[uiDepth]->copyToPicYuv  ( pcPic->getPicYuvRec (), uiCtuRsAddr, uiZorderIdx );

  return;
}
#endif

Void TDecCu::xDecodeInterTexture ( TComDataCU* pcCU, UInt uiDepth )
{
#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[pcCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#endif
  TComTURecurse tuRecur(pcCU, 0, uiDepth);

#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER
#if JVET_C0024_QTBT
  TComYuv* pcPred = m_pppcYuvReco[uiWidthIdx][uiHeightIdx];
#else
  TComYuv* pcPred = m_ppcYuvReco[uiDepth];
#endif
#endif
  for(UInt ch=0; ch<pcCU->getPic()->getNumberValidComponents(); ch++)
  {
    const ComponentID compID=ComponentID(ch);
    DEBUG_STRING_OUTPUT(std::cout, debug_reorder_data_inter_token[compID])
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER
#if JVET_C0024_QTBT
    m_pcTrQuant->invRecurTransformNxN ( compID, m_pppcYuvResi[uiWidthIdx][uiHeightIdx], tuRecur, pcPred);
#else
    m_pcTrQuant->invRecurTransformNxN ( compID, m_ppcYuvResi[uiDepth], tuRecur, pcPred);
#endif
#else
#if JVET_C0024_QTBT
    m_pcTrQuant->invRecurTransformNxN ( compID, m_pppcYuvResi[uiWidthIdx][uiHeightIdx], tuRecur
#if PVC_AQS
    ,m_ppcYuvReco[uiDepth]
#else
	);
#endif

#else
    m_pcTrQuant->invRecurTransformNxN ( compID, m_ppcYuvResi[uiDepth], tuRecur 
#if PVC_AQS
    ,m_ppcYuvReco[uiDepth]
#else
	);
#endif
#endif
#endif
  }

  DEBUG_STRING_OUTPUT(std::cout, debug_reorder_data_inter_token[MAX_NUM_COMPONENT])
}

/** Function for deriving reconstructed luma/chroma samples of a PCM mode CU.
 * \param pcCU pointer to current CU
 * \param uiPartIdx part index
 * \param piPCM pointer to PCM code arrays
 * \param piReco pointer to reconstructed sample arrays
 * \param uiStride stride of reconstructed sample arrays
 * \param uiWidth CU width
 * \param uiHeight CU height
 * \param compID colour component ID
 * \returns Void
 */
Void TDecCu::xDecodePCMTexture( TComDataCU* pcCU, const UInt uiPartIdx, const Pel *piPCM, Pel* piReco, const UInt uiStride, const UInt uiWidth, const UInt uiHeight, const ComponentID compID)
{
        Pel* piPicReco         = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu()+uiPartIdx);
  const UInt uiPicStride       = pcCU->getPic()->getPicYuvRec()->getStride(compID);
  const TComSPS &sps           = *(pcCU->getSlice()->getSPS());
  const UInt uiPcmLeftShiftBit = sps.getBitDepth(toChannelType(compID)) - sps.getPCMBitDepth(toChannelType(compID));

  for(UInt uiY = 0; uiY < uiHeight; uiY++ )
  {
    for(UInt uiX = 0; uiX < uiWidth; uiX++ )
    {
      piReco[uiX] = (piPCM[uiX] << uiPcmLeftShiftBit);
      piPicReco[uiX] = piReco[uiX];
    }
    piPCM += uiWidth;
    piReco += uiStride;
    piPicReco += uiPicStride;
  }
}

/** Function for reconstructing a PCM mode CU.
 * \param pcCU pointer to current CU
 * \param uiDepth CU Depth
 * \returns Void
 */
Void TDecCu::xReconPCM( TComDataCU* pcCU, UInt uiDepth )
{
#if JVET_C0024_QTBT 
  UInt uiWidthIdx = g_aucConvertToBit[pcCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
  const UInt maxCuWidth     = pcCU->getSlice()->getSPS()->getMaxCUWidth();
  const UInt maxCuHeight    = pcCU->getSlice()->getSPS()->getMaxCUHeight();
#endif
  for (UInt ch=0; ch < pcCU->getPic()->getNumberValidComponents(); ch++)
  {
    const ComponentID compID = ComponentID(ch);
#if JVET_C0024_QTBT
    const UInt width  = pcCU->getWidth(0) >> m_pppcYuvResi[uiWidthIdx][uiHeightIdx]->getComponentScaleX(compID);
    const UInt height = pcCU->getHeight(0) >> m_pppcYuvResi[uiWidthIdx][uiHeightIdx]->getComponentScaleY(compID);
    const UInt stride = m_pppcYuvResi[uiWidthIdx][uiHeightIdx]->getStride(compID);
    Pel * pPCMChannel = pcCU->getPCMSample(compID);
    Pel * pRecChannel = m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->getAddr(compID);
#else
    const UInt width  = (maxCuWidth >>(uiDepth+m_ppcYuvResi[uiDepth]->getComponentScaleX(compID)));
    const UInt height = (maxCuHeight>>(uiDepth+m_ppcYuvResi[uiDepth]->getComponentScaleY(compID)));
    const UInt stride = m_ppcYuvResi[uiDepth]->getStride(compID);
    Pel * pPCMChannel = pcCU->getPCMSample(compID);
    Pel * pRecChannel = m_ppcYuvReco[uiDepth]->getAddr(compID);
#endif
    xDecodePCMTexture( pcCU, 0, pPCMChannel, pRecChannel, stride, width, height, compID );
  }
}

/** Function for filling the PCM buffer of a CU using its reconstructed sample array
 * \param pCU   pointer to current CU
 * \param depth CU Depth
 */
Void TDecCu::xFillPCMBuffer(TComDataCU* pCU, UInt depth)
{
  const ChromaFormat format = pCU->getPic()->getChromaFormat();
  const UInt numValidComp   = getNumberValidComponents(format);
#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[pCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[pCU->getHeight(0)];
#else
  const UInt maxCuWidth     = pCU->getSlice()->getSPS()->getMaxCUWidth();
  const UInt maxCuHeight    = pCU->getSlice()->getSPS()->getMaxCUHeight();
#endif

  for (UInt componentIndex = 0; componentIndex < numValidComp; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

#if JVET_C0024_QTBT
    const UInt width  = pCU->getWidth(0) >> getComponentScaleX(component, format);
    const UInt height = pCU->getHeight(0) >> getComponentScaleY(component, format);
    Pel *source      = m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->getAddr(component, 0, width);
    Pel *destination = pCU->getPCMSample(component);

    const UInt sourceStride = m_pppcYuvReco[uiWidthIdx][uiHeightIdx]->getStride(component);
#else
    const UInt width  = maxCuWidth  >> (depth + getComponentScaleX(component, format));
    const UInt height = maxCuHeight >> (depth + getComponentScaleY(component, format));

    Pel *source      = m_ppcYuvReco[depth]->getAddr(component, 0, width);
    Pel *destination = pCU->getPCMSample(component);

    const UInt sourceStride = m_ppcYuvReco[depth]->getStride(component);
#endif

    for (Int line = 0; line < height; line++)
    {
      for (Int column = 0; column < width; column++)
      {
        destination[column] = source[column];
      }

      source      += sourceStride;
      destination += width;
    }
  }
}

#if VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
Void TDecCu::xDeriveCUMV( TComDataCU * pcCU , UInt uiAbsPartIdx , UInt uiDepth )
{
#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[pcCU->getWidth(uiAbsPartIdx)];
  UInt uiHeightIdx = g_aucConvertToBit[pcCU->getHeight(uiAbsPartIdx)];

  PartSize ePartSize = SIZE_2Nx2N;
  UInt uiNumPU = 1;
  UInt uiPUOffset = 0;
  TComDataCU * pcSubCU = m_pppcCU[uiWidthIdx][uiHeightIdx];
#else
  PartSize ePartSize = pcCU->getPartitionSize( uiAbsPartIdx );
  UInt uiNumPU = ( ePartSize == SIZE_2Nx2N ? 1 : ( ePartSize == SIZE_NxN ? 4 : 2 ) );
  UInt uiPUOffset = ( g_auiPUOffset[UInt( ePartSize )] << ( ( pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth ) << 1 ) ) >> 4;
  TComDataCU * pcSubCU = m_ppcCU[uiDepth];
#endif

  TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  memset( uhInterDirNeighbours , 0 , sizeof( uhInterDirNeighbours ) );
  Int numValidMergeCand = 0;
  Bool isMerged = false;

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UChar    eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
  memset( eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, sizeof(UChar)*MRG_MAX_NUM_CANDS);
#endif
#if VCEG_AZ06_IC
  Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif
  for ( UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset )
  {
    if ( pcCU->getMergeFlag( uiSubPartIdx ) )
    {
#if VCEG_AZ07_FRUC_MERGE
      if( pcCU->getFRUCMgrMode( uiSubPartIdx ) )
      {
        Bool bAvailable = m_pcPrediction->deriveFRUCMV( pcSubCU , uiDepth , uiSubPartIdx - pcSubCU->getZorderIdxInCtu() , uiPartIdx );
        if (!bAvailable)
        {
          int a = 0;
        }
        assert( bAvailable );
      }
      else
#endif
      {
#if MRG_CONSIDER_COMPLEX
        if (pcCU->getMergeComplexFlag(uiSubPartIdx))
        {
          TComSlice* pSlice = pcCU->getSlice();
          TComMvField  cMvFieldNeighbourNodes[MRG_NODE_NUM << 1];
          UChar InterDirNeighbourNodes[MRG_NODE_NUM];
          Int numExistMergeNode[2] = { 0 }; // numExistMergeNode[i] represents number of Exist merge nodes  of list i (i=0,1) without removing duplication
          Int numValidMergeNode[MRG_COMPLEX_NUM_CANDS << 1] = { 0 };// numValidMergeNode[i] represents number of valid merge nodes of list i  (i=0,1) after removing duplication,
          UChar numValidMergeComplexCand = 0;
          Int RealComplexIndex[11] = { 0 };
          for (UInt ui = 0; ui < MRG_NODE_NUM; ++ui)
          {
            InterDirNeighbourNodes[ui] = 0;
          }

          UInt MergeComplexIndex = pcCU->getMergeComplexIndex(uiSubPartIdx);
          pcSubCU->getInterMergeComplexCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, numValidMergeComplexCand, cMvFieldNeighbourNodes, InterDirNeighbourNodes, numExistMergeNode, numValidMergeNode, RealComplexIndex, MergeComplexIndex);

          UChar uhInterDirComplex = 0;
          // pre process the control motion information
          pcCU->preproccessMotionFieldNodes(RealComplexIndex, MergeComplexIndex, InterDirNeighbourNodes, uhInterDirComplex, numValidMergeNode, cMvFieldNeighbourNodes);
          pcCU->setInterDirSubParts(uhInterDirComplex, uiSubPartIdx, uiPartIdx, uiDepth);

          NodeInfo NodeMotionInfo;
          for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
          {
            NodeMotionInfo.MvFieldNbNodes[i] = cMvFieldNeighbourNodes[i];
          }
          for (Int i = 0; i < MRG_NODE_NUM; i++)
          {
            NodeMotionInfo.InterDirNbNodes[i] = InterDirNeighbourNodes[i];
          }
          for (Int i = 0; i < (MRG_COMPLEX_NUM_CANDS << 1); i++)
          {
            NodeMotionInfo.numValidNodes[i] = numValidMergeNode[i];
          }
          for (Int i = 0; i < MRG_COMPLEX_NUM_CANDS; i++)
          {
            NodeMotionInfo.RealComplexIdx[i] = RealComplexIndex[i];
          }
#if JVET_C0024_QTBT
          pcCU->setNodeMotionInfoSubParts(NodeMotionInfo, uiSubPartIdx);
#else
          pcCU->setNodeMotionInfoSubParts(NodeMotionInfo, uiSubPartIdx, uiPartIdx, uiDepth);
#endif

          TComMv cTmpMv(0, 0);
          if (uhInterDirComplex == 1)
          {
            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0) > 0) //if ( ref. frame list0 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_0, uiSubPartIdx, uiDepth, uiPartIdx);
            }

            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 0) //if ( ref. frame list1 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(TComMvField(), ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
            }
          }

          if (uhInterDirComplex == 2)
          {
            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0) > 0) //if ( ref. frame list0 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(TComMvField(), ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
            }
            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 0) //if ( ref. frame list1 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_1, uiSubPartIdx, uiDepth, uiPartIdx);
            }
          }

          if (uhInterDirComplex == 3)
          {
            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0) > 0) //if ( ref. frame list0 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_0, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_0, uiSubPartIdx, uiDepth, uiPartIdx);

            }
            if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 0) //if ( ref. frame list1 has at least 1 entry )
            {
              pcCU->setMVPIdxSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts(0, REF_PIC_LIST_1, uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
              pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_1, uiSubPartIdx, uiDepth, uiPartIdx);
            }
          }
        }
        else
        {
#endif

#if COM16_C1016_AFFINE
        if ( pcCU->isAffine( uiSubPartIdx ) )
        {
          TComMvField  cAffineMvField[2][3]; // double length for mv of both lists, 3 mv for affine
          UInt uiMergeIndex = 0;

#if JVECT_C0062_AFFINE_SIX_PARAM
          UInt neighborAffineParam = 0;
#if JVET_C0024_QTBT
          m_pppcCU[uiWidthIdx][uiHeightIdx]->getAffineMergeCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#else
          m_ppcCU[uiDepth]->getAffineMergeCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
#endif
          pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiSubPartIdx, uiPartIdx, uiDepth);
#else
#if JVET_C0024_QTBT
          m_pppcCU[uiWidthIdx][uiHeightIdx]->getAffineMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx,cAffineMvField, uhInterDirNeighbours, numValidMergeCand );
#else
          m_ppcCU[uiDepth]->getAffineMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx,cAffineMvField, uhInterDirNeighbours, numValidMergeCand );
#endif
#endif
          pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth );

          TComMv cTmpMv( 0, 0 );
          for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              TComMvField* pcMvField = cAffineMvField[ 2 * uiMergeIndex + uiRefListIdx ];

              pcCU->setMVPIdxSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx, uiPartIdx, uiDepth);
              pcCU->setMVPNumSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx , uiPartIdx, uiDepth);
              pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cTmpMv, SIZE_2Nx2N, uiSubPartIdx, uiDepth );
#if JVECT_C0062_AFFINE_SIX_PARAM
              m_pppcCU[uiWidthIdx][uiHeightIdx]->setAllAffineMvField(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, pcMvField, RefPicList(uiRefListIdx), uiDepth, neighborAffineParam);
#else
              pcCU->setAllAffineMvField( uiSubPartIdx , uiPartIdx, pcMvField, RefPicList(uiRefListIdx), uiDepth );
#endif
            }
          }
        }
        else
        {
#endif
        UInt uiMergeIndex = pcCU->getMergeIndex(uiSubPartIdx);
        if ( pcCU->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() && ePartSize != SIZE_2Nx2N && pcSubCU->getWidth( 0 ) <= 8 ) 
        {
#if !JVET_C0024_QTBT
          pcSubCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth );
#endif
          if ( !isMerged )
          {
            pcSubCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand 
#if VCEG_AZ06_IC
              , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
              , eMergeCandTypeNieghors , m_pMvFieldSP , m_phInterDirSP
#endif
              );
            isMerged = true;
          }
#if !JVET_C0024_QTBT
          pcSubCU->setPartSizeSubParts( ePartSize, 0, uiDepth );
#endif
        }
        else
        {
          uiMergeIndex = pcCU->getMergeIndex(uiSubPartIdx);
          pcSubCU->getInterMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#if VCEG_AZ06_IC
            , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
            , eMergeCandTypeNieghors , m_pMvFieldSP , m_phInterDirSP  
#endif
            );
        }
        pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth );
#if VCEG_AZ06_IC
        if( ePartSize == SIZE_2Nx2N )
        {
          pcCU->setICFlagSubParts( pcCU->getSlice()->getApplyIC() ? abICFlag[uiMergeIndex] : 0, uiSubPartIdx, uiDepth );
        }
#endif
        TComMv cTmpMv( 0, 0 );
        for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
        {        
          if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
          {
            pcCU->setMVPIdxSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx, uiPartIdx, uiDepth);
            pcCU->setMVPNumSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx, uiPartIdx, uiDepth);
            pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );
            pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvField( cMvFieldNeighbours[ 2*uiMergeIndex + uiRefListIdx ], ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );
          }
        }
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
        pcCU->setMergeTypeSubParts( eMergeCandTypeNieghors[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth ); 
#if JVET_C0035_ATMVP_SIMPLIFICATION
        if( eMergeCandTypeNieghors[uiMergeIndex]==MGR_TYPE_SUBPU_ATMVP || eMergeCandTypeNieghors[uiMergeIndex]==MGR_TYPE_SUBPU_ATMVP_EXT )
#else
        if( eMergeCandTypeNieghors[uiMergeIndex]==MGR_TYPE_SUBPU_TMVP || eMergeCandTypeNieghors[uiMergeIndex]==MGR_TYPE_SUBPU_TMVP_EXT )
#endif
        {
          Int iWidth, iHeight;
          UInt uiIdx;
          pcCU->getPartIndexAndSize( uiPartIdx, uiIdx, iWidth, iHeight, uiSubPartIdx, true );

          UInt uiSPAddr;

          Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
          UInt uiSPListIndex =  eMergeCandTypeNieghors[uiMergeIndex];
#else
          UInt uiSPListIndex =  eMergeCandTypeNieghors[uiMergeIndex]==MGR_TYPE_SUBPU_TMVP ? 0:1;
#endif
          pcCU->getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

          for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
          {
            pcCU->getSPAbsPartIdx(uiSubPartIdx, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
            pcCU->setInterDirSP(m_phInterDirSP[uiSPListIndex][iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
            pcCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx], iSPWidth, iSPHeight);
            pcCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx + 1], iSPWidth, iSPHeight);
          }
        }
#endif
#if COM16_C1016_AFFINE
        }
#endif
#if MRG_CONSIDER_COMPLEX
        }
#endif
      }
    }
    else
    {
      for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
      {
        RefPicList eRefList = RefPicList( uiRefListIdx );
        if ( pcCU->getSlice()->getNumRefIdx( eRefList ) > 0 && ( pcCU->getInterDir( uiSubPartIdx ) & ( 1 << uiRefListIdx ) ) )
        {
#if COM16_C1016_AFFINE
          if ( pcCU->isAffine(uiSubPartIdx) )
          {
            TComMv acMv[3];
            memset( acMv, 0, sizeof(TComMv) * 3 );

#if JVECT_C0062_AFFINE_SIX_PARAM
            Bool b6Param = pcCU->getAffineParamFlag(uiSubPartIdx);
#endif

            UInt uiSubCUPartIdx = uiSubPartIdx - uiAbsPartIdx;
            TComCUMvField* pcSubCUMvField = pcSubCU->getCUMvField( eRefList );
            Int iRefIdx = pcSubCUMvField->getRefIdx( uiSubCUPartIdx ); 

#if JVECT_C0062_AFFINE_SIX_PARAM
            assert(iRefIdx >= 0);
            AffineAMVPInfo* pAffineAMVPInfo = pcSubCUMvField->getAffineAMVPInfo(b6Param);
            pcSubCU->fillAffineMvpCand(uiPartIdx, uiSubCUPartIdx, eRefList, iRefIdx, pAffineAMVPInfo, b6Param);
            m_pcPrediction->getMvPredAffineAMVP(pcSubCU, uiPartIdx, uiSubCUPartIdx, eRefList, acMv, b6Param);
#else
            AffineAMVPInfo* pAffineAMVPInfo = pcSubCUMvField->getAffineAMVPInfo();
            pcSubCU->fillAffineMvpCand( uiPartIdx, uiSubCUPartIdx, eRefList, iRefIdx, pAffineAMVPInfo );
            m_pcPrediction->getMvPredAffineAMVP( pcSubCU, uiPartIdx, uiSubCUPartIdx, eRefList, acMv );
#endif


#if JVECT_C0062_AFFINE_SIX_PARAM
            UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB, uiAbsIndexInLCU;
            uiAbsIndexInLCU = pcSubCU->getZorderIdxInCtu();
            pcSubCU->deriveLeftRightTopIdx(uiPartIdx, uiPartIdxLT, uiPartIdxRT);
            pcSubCU->deriveLeftBottomIdx(uiPartIdx, uiPartIdxLB);
#if TT_TREE_AFFINE // deriveZIdxInCU
            assert(pcSubCU->getUseAbsoluteZIdxToCtu()==true);
#endif
            if (b6Param)
            {
              acMv[0] += pcSubCUMvField->getMvd(uiPartIdxLT - uiAbsIndexInLCU);
              acMv[1] += pcSubCUMvField->getMvd(uiPartIdxRT - uiAbsIndexInLCU);
              acMv[2] += pcSubCUMvField->getMvd(uiPartIdxLB - uiAbsIndexInLCU);
            }
            else
            {
              Int iWidth = pcSubCU->getWidth(uiSubCUPartIdx);
              Int iHeight = pcSubCU->getHeight(uiSubCUPartIdx);

              Bool bRecShape = (iHeight / iWidth) >= 2;
              if (bRecShape)
              {
                acMv[0] += pcSubCUMvField->getMvd(uiPartIdxLT - uiAbsIndexInLCU);
                acMv[2] += pcSubCUMvField->getMvd(uiPartIdxLB - uiAbsIndexInLCU);
#if !AFFINE_CLIP_REMOVE
                pcSubCU->clipMv(acMv[0]);//ZF
                pcSubCU->clipMv(acMv[2]);//ZF
#endif
                Int vx1 = (acMv[2].getVer() - acMv[0].getVer())*iWidth / iHeight + acMv[0].getHor();
                Int vy1 = -(acMv[2].getHor() - acMv[0].getHor())*iWidth / iHeight + acMv[0].getVer();
                acMv[1].set(vx1, vy1);
              }
              else
              {
                acMv[0] += pcSubCUMvField->getMvd(uiPartIdxLT - uiAbsIndexInLCU);
                acMv[1] += pcSubCUMvField->getMvd(uiPartIdxRT - uiAbsIndexInLCU);
#if !AFFINE_CLIP_REMOVE
                pcSubCU->clipMv(acMv[0]);//ZF
                pcSubCU->clipMv(acMv[1]);//ZF
#endif
                Int vx2 = -(acMv[1].getVer() - acMv[0].getVer()) * iHeight / iWidth + acMv[0].getHor();
                Int vy2 = (acMv[1].getHor() - acMv[0].getHor()) * iHeight / iWidth + acMv[0].getVer();
                acMv[2].set(vx2, vy2);
              }
            }
#else
            UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB, uiAbsIndexInLCU;
            uiAbsIndexInLCU = pcSubCU->getZorderIdxInCtu();
            pcSubCU->deriveLeftRightTopIdx( uiPartIdx, uiPartIdxLT, uiPartIdxRT );
            pcSubCU->deriveLeftBottomIdx( uiPartIdx, uiPartIdxLB );

#if TT_TREE_AFFINE // deriveZIdxInCU
            assert(pcSubCU->getUseAbsoluteZIdxToCtu()==true);
#endif
            acMv[0] += pcSubCUMvField->getMvd( uiPartIdxLT - uiAbsIndexInLCU );
            acMv[1] += pcSubCUMvField->getMvd( uiPartIdxRT - uiAbsIndexInLCU );

            Int iWidth = pcSubCU->getWidth(uiSubCUPartIdx);
            Int iHeight = pcSubCU->getHeight(uiSubCUPartIdx);
            Int vx2 =  - ( acMv[1].getVer() - acMv[0].getVer() ) * iHeight / iWidth + acMv[0].getHor();
            Int vy2 =    ( acMv[1].getHor() - acMv[0].getHor() ) * iHeight / iWidth + acMv[0].getVer();
            acMv[2].set( vx2, vy2 );
#endif
#if !AFFINE_CLIP_REMOVE
            pcSubCU->clipMv(acMv[0]);
            pcSubCU->clipMv(acMv[1]);
            pcSubCU->clipMv(acMv[2]);
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
            pcSubCU->setAllAffineMv(uiSubCUPartIdx, uiPartIdx, acMv, eRefList, 0, b6Param);
#else
            pcSubCU->setAllAffineMv( uiSubCUPartIdx, uiPartIdx, acMv, eRefList, 0 );
#endif
          }
          else
          {
#endif

          Short iMvdHor = pcCU->getCUMvField( eRefList )->getMvd( uiSubPartIdx ).getHor();
          Short iMvdVer = pcCU->getCUMvField( eRefList )->getMvd( uiSubPartIdx ).getVer();
#if VCEG_AZ07_IMV
          if( pcCU->getiMVFlag( uiSubPartIdx ) )
          {
            iMvdHor <<= 2;
            iMvdVer <<= 2;

#if  JVET_E0076_MULTI_PEL_MVD
            if( pcCU->getiMVFlag( uiSubPartIdx ) == 2)
            {
              iMvdHor <<= MULTI_PEL_MVD_BITS;
              iMvdVer <<= MULTI_PEL_MVD_BITS;
            }
#endif
          }
#endif
          TComMv cMv( iMvdHor, iMvdVer );

          pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );

          AMVPInfo* pAMVPInfo = pcSubCU->getCUMvField( eRefList )->getAMVPInfo();
#if JVET_C0024_QTBT && !VCEG_AZ07_FRUC_MERGE
          pcSubCU->fillMvpCand(uiPartIdx, uiSubPartIdx - uiAbsPartIdx, eRefList, pcSubCU->getCUMvField( eRefList )->getRefIdx( uiSubPartIdx - uiAbsPartIdx ), pAMVPInfo );
#else
          pcSubCU->fillMvpCand(uiPartIdx, uiSubPartIdx - uiAbsPartIdx, eRefList, pcSubCU->getCUMvField( eRefList )->getRefIdx( uiSubPartIdx - uiAbsPartIdx ), pAMVPInfo, m_pcPrediction );
#endif
          m_pcPrediction->getMvPredAMVP( pcSubCU, uiPartIdx, uiSubPartIdx - uiAbsPartIdx, RefPicList( uiRefListIdx ), cMv );
          cMv += TComMv( iMvdHor, iMvdVer );
         
#if TT_TREE_DS3 //important fix when TT is ON
          pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMv( cMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
#else
          pcSubCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMv( cMv, ePartSize, uiSubPartIdx - uiAbsPartIdx, 0, uiPartIdx);
#endif
#if COM16_C1016_AFFINE
          }
#endif
        }
      }
    }

#if JVET_C0024_QTBT && !VCEG_AZ07_FRUC_MERGE
    if ( (pcCU->getInterDir(uiSubPartIdx) == 3) && pcSubCU->isBipredRestriction(uiPartIdx) ) 
#else
    if ( (pcCU->getInterDir(uiSubPartIdx) == 3) && pcSubCU->isBipredRestriction(uiPartIdx) && !pcCU->getFRUCMgrMode( uiSubPartIdx ) ) 
#endif
    {
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMv( TComMv(0,0), ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllRefIdx( -1, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->setInterDirSubParts( 1, uiSubPartIdx, uiPartIdx, uiDepth);
    }
  }
}
#endif

#if BARC
Void TDecCu::upSampleRec ( TComDataCU*& pCtu, ChannelType chType)
{
  Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false : true;
  Bool isLeftCUAvail = pCtu->getCUPelX() == 0 ? false : true;

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

  UInt upWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt upHeight = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiWidth = upWidth >> BARC_DScaleRatio;
  UInt uiHeight = upHeight >> BARC_DScaleRatio;
  UInt iSrcStride =  pCtu->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
  UInt starIdx, endIdx;

  if (m_uiCUPelX + upWidth > picWidth)
  {
    upWidth = picWidth - m_uiCUPelX;
    uiWidth = upWidth >> BARC_DScaleRatio;
  }
  if (m_uiCUPelY + upHeight > picHeight)
  {
    upHeight = picHeight - m_uiCUPelY;
    uiHeight = upHeight >> BARC_DScaleRatio;
  }

  if (chType == CHANNEL_TYPE_LUMA)
  {
    starIdx = 0;
    endIdx = 1;
  }
  else
  {
    starIdx = 1;
    endIdx = 3;
    upWidth >>= 1;
    upHeight >>= 1;
    uiWidth >>= 1;
    uiHeight >>= 1;
  }
  for (UInt comp=starIdx; comp<endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID,pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);

    for (Int i = 0; i < uiHeight; i++)
    memcpy(pRefScale+i*iSrcStride, pSrc+i*iSrcStride, uiWidth*sizeof(Pel));

    pCtu->upScaleComponent2x2(pRefScale, pSrc, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail,isLeftCUAvail, compID);
  }
#if BARC_REF_UP_SAMPLE_CNN
#if BARC_INTRA_LUMA_SPS
  if (pCtu->getSlice()->getSPS()->getUseBarcIntraLuma())
  {
#endif
	  if ((chType == CHANNEL_TYPE_LUMA) && pCtu->getUpCNNFlag())
	  {
		  const ComponentID compID = ComponentID(starIdx);

		  Pel* pSrc = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pCtu->getCtuRsAddr(), 0);
		  Pel* pSrcRec = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
		  UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
		  Double* pDstNet = new Double[MAX_CU_SIZE*MAX_CU_SIZE];
		  UInt iDstStride = MAX_CU_SIZE;
		  pCtu->getPic()->CNNUpScaleLumaComponent2x2(true, pSrc, pDstNet, iSrcStride, iDstStride, uiWidth, uiHeight);
		  Double* tpDstNet = pDstNet;
		  for (Int i = 0; i < upHeight; i++)
		  {
			  for (Int j = 0; j < upWidth; j++)
			  {

#if COMPACT_RESOLUTION_CNN
				  pSrcRec[j] = ClipA(int(tpDstNet[j] * PixelRangeInt + 0.5), compID);
#else
				  pSrcRec[j] = ClipA(int(pSrcRec[j] + tpDstNet[j] * PixelRangeInt + 0.5), compID);
#endif
			  }
			  pSrcRec += iSrcStride;
			  tpDstNet += iDstStride;
		  }
		  if (pDstNet)
		  {
			  delete[] pDstNet;
			  pDstNet = NULL;
		  }
	  }
#if BARC_INTRA_LUMA_SPS
  }
#endif
#endif
}

Void TDecCu::downSampleRec ( TComDataCU*& pCtu, ChannelType chType)
{
  UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

  Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false:true;
  Bool isLeftCUAvail = pCtu->getCUPelX() == 0? false:true;

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();
  Double distortion = 0;
  UInt startIdx, endIdx;		

  if (m_uiCUPelX + uiWidth > picWidth)
  {
    uiWidth = picWidth - m_uiCUPelX;
  }
  if (m_uiCUPelY + uiHeight > picHeight)
  {
    uiHeight = picHeight - m_uiCUPelY;
  }

  if (chType == CHANNEL_TYPE_LUMA)
  {
    startIdx = 0;
    endIdx = 1;
  }
  else
  {
    startIdx = 1;
    endIdx = 3;
    uiWidth >>= 1;
    uiHeight >>= 1;
  }

  for (UInt comp=startIdx; comp<endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID,pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
#if COMPACT_RESOLUTION_CNN
	if (chType == CHANNEL_TYPE_LUMA)
	{
#if BARC_INTRA_LUMA_SPS
		if (pCtu->getSlice()->getSPS()->getUseBarcIntraLuma())
		{
#endif
			pCtu->getPic()->CNNDownScaleLumaComponent2x2(false, pSrc, pRefScale, iSrcStride, iSrcStride, uiWidth, uiHeight);
#if BARC_INTRA_LUMA_SPS
		}
		else
		{
			pCtu->downScaleComponent2x2ZeroPhaseBoundary(pSrc, pRefScale, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, compID);
		}
#endif
	}
	else
	{
		pCtu->downScaleComponent2x2ZeroPhaseBoundary(pSrc, pRefScale, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, compID);
	}
#else
    pCtu->downScaleComponent2x2ZeroPhaseBoundary(pSrc, pRefScale, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, compID);
#endif
  }
}
Void TDecCu::copyLumaDscaletoRec(TComDataCU* pCtu, Bool isCopyDscaletoRec)
{
    UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
    UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

    UInt m_uiCUPelX = pCtu->getCUPelX();
    UInt m_uiCUPelY = pCtu->getCUPelY();
    UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
    UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

    if (m_uiCUPelX + uiWidth > picWidth)
    {
        uiWidth = picWidth - m_uiCUPelX;
    }
    if (m_uiCUPelY + uiHeight > picHeight)
    {
        uiHeight = picHeight - m_uiCUPelY;
    }

    for (UInt comp = 0; comp < 1; comp++)
    {
        const ComponentID compID = ComponentID(comp);

        Pel* pRec = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
        Pel* pRecScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pCtu->getCtuRsAddr(), 0);
        Pel* pTmp = m_reconLumaBuffer;
  
        UInt iRecStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
        UInt iTmpStride = MAX_CU_SIZE;

        if (isCopyDscaletoRec)
        {
          for (Int y = 0; y < uiHeight; y++)
            memcpy(pTmp + y*iTmpStride, pRec + y*iRecStride, sizeof(Pel)*uiWidth);
          for (Int y = 0; y < uiHeight / 2; y++)
            memcpy(pRec + y*iRecStride, pRecScale + y*iRecStride, sizeof(Pel)*uiWidth / 2);
        }
        else
        {
          for (Int y = 0; y < uiHeight; y++)
            memcpy(pRec + y*iRecStride, pTmp + y*iTmpStride, sizeof(Pel)*uiWidth);
        }
    }
}
#endif
//! \}
