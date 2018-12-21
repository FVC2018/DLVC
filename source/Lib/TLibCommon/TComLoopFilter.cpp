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

/** \file     TComLoopFilter.cpp
    \brief    deblocking filter
*/

#include "TComLoopFilter.h"
#include "TComSlice.h"
#include "TComMv.h"
#include "TComTU.h"

#if JVET_D0033_ADAPTIVE_CLIPPING
#include "CommonDef.h"
#endif
//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constants
// ====================================================================================================================

//#define   EDGE_VER    0
//#define   EDGE_HOR    1

#define DEFAULT_INTRA_TC_OFFSET 2 ///< Default intra TC offset

// ====================================================================================================================
// Tables
// ====================================================================================================================

const UChar TComLoopFilter::sm_tcTable[MAX_QP + 1 + DEFAULT_INTRA_TC_OFFSET] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,7,8,9,10,11,13,14,16,18,20,22,24
};

const UChar TComLoopFilter::sm_betaTable[MAX_QP + 1] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64
};

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TComLoopFilter::TComLoopFilter()
: m_uiNumPartitions(0)
, m_bLFCrossTileBoundary(true)
{
  for( Int edgeDir = 0; edgeDir < NUM_EDGE_DIR; edgeDir++ )
  {
    m_aapucBS       [edgeDir] = NULL;
    m_aapbEdgeFilter[edgeDir] = NULL;
  }
}

TComLoopFilter::~TComLoopFilter()
{
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================
Void TComLoopFilter::setCfg( Bool bLFCrossTileBoundary )
{
  m_bLFCrossTileBoundary = bLFCrossTileBoundary;
}

Void TComLoopFilter::create( UInt uiMaxCUDepth )
{
  destroy();
  m_uiNumPartitions = 1 << ( uiMaxCUDepth<<1 );
  for( Int edgeDir = 0; edgeDir < NUM_EDGE_DIR; edgeDir++ )
  {
    m_aapucBS       [edgeDir] = new UChar[m_uiNumPartitions];
    m_aapbEdgeFilter[edgeDir] = new Bool [m_uiNumPartitions];
  }
}

Void TComLoopFilter::destroy()
{
  for( Int edgeDir = 0; edgeDir < NUM_EDGE_DIR; edgeDir++ )
  {
    if (m_aapucBS[edgeDir] != NULL)
    {
      delete [] m_aapucBS[edgeDir];
      m_aapucBS[edgeDir] = NULL;
    }

    if (m_aapbEdgeFilter[edgeDir])
    {
      delete [] m_aapbEdgeFilter[edgeDir];
      m_aapbEdgeFilter[edgeDir] = NULL;
    }
  }
}

/**
 - call deblocking function for every CU
 .
 \param  pcPic   picture class (TComPic) pointer
 */
Void TComLoopFilter::loopFilterPic( TComPic* pcPic )
{
  // Horizontal filtering
  for ( UInt ctuRsAddr = 0; ctuRsAddr < pcPic->getNumberOfCtusInFrame(); ctuRsAddr++ )
  {
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );

    ::memset( m_aapucBS       [EDGE_VER], 0, sizeof( UChar ) * m_uiNumPartitions );
    ::memset( m_aapbEdgeFilter[EDGE_VER], 0, sizeof( Bool  ) * m_uiNumPartitions );

    // CU-based deblocking
#if JVET_C0024_QTBT
    pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
    UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize() ;
#if BARC_DBK
    if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
      xDeblockCU(pCtu, 0, BARC_DScaleRatio, uiCTUSize >> BARC_DScaleRatio, uiCTUSize >> BARC_DScaleRatio, EDGE_VER);
    else
#endif
    xDeblockCU( pCtu, 0, 0, uiCTUSize, uiCTUSize,EDGE_VER );
    if (pCtu->getSlice()->isIntra())
    {
      ::memset( m_aapucBS       [EDGE_VER], 0, sizeof( UChar ) * m_uiNumPartitions );
      ::memset( m_aapbEdgeFilter[EDGE_VER], 0, sizeof( Bool  ) * m_uiNumPartitions );

      pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
      xDeblockCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, EDGE_VER );
    }
#else
    xDeblockCU( pCtu, 0, 0, EDGE_VER );
#endif
  }

  // Vertical filtering
  for ( UInt ctuRsAddr = 0; ctuRsAddr < pcPic->getNumberOfCtusInFrame(); ctuRsAddr++ )
  {
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );

    ::memset( m_aapucBS       [EDGE_HOR], 0, sizeof( UChar ) * m_uiNumPartitions );
    ::memset( m_aapbEdgeFilter[EDGE_HOR], 0, sizeof( Bool  ) * m_uiNumPartitions );

    // CU-based deblocking
#if JVET_C0024_QTBT
    pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
    UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize() ;
#if BARC_DBK
    if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
      xDeblockCU(pCtu, 0, BARC_DScaleRatio, uiCTUSize >> BARC_DScaleRatio, uiCTUSize >> BARC_DScaleRatio, EDGE_HOR);
    else
#endif
    xDeblockCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, EDGE_HOR );
    if (pCtu->getSlice()->isIntra())
    {
      ::memset( m_aapucBS       [EDGE_HOR], 0, sizeof( UChar ) * m_uiNumPartitions );
      ::memset( m_aapbEdgeFilter[EDGE_HOR], 0, sizeof( Bool  ) * m_uiNumPartitions );

      pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
      xDeblockCU( pCtu, 0, 0, uiCTUSize, uiCTUSize, EDGE_HOR );
    }
#else
    xDeblockCU( pCtu, 0, 0, EDGE_HOR );
#endif
  }
}


// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

/**
 Deblocking filter process in CU-based (the same function as conventional's)

 \param pcCU             Pointer to CTU/CU structure
 \param uiAbsZorderIdx   Position in CU
 \param uiDepth          Depth in CU
 \param edgeDir          the direction of the edge in block boundary (horizontal/vertical), which is added newly
*/
#if JVET_C0024_QTBT
Void TComLoopFilter::xDeblockCU( TComDataCU* pcCU, UInt uiAbsZorderIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, DeblockEdgeDir edgeDir )
#else
Void TComLoopFilter::xDeblockCU( TComDataCU* pcCU, UInt uiAbsZorderIdx, UInt uiDepth, DeblockEdgeDir edgeDir )
#endif
{
#if JVET_C0024_QTBT 
  if(pcCU->getPic()==0)
#else
  if(pcCU->getPic()==0||pcCU->getPartitionSize(uiAbsZorderIdx)==NUMBER_OF_PART_SIZES)
#endif
  {
    return;
  }
  TComPic* pcPic     = pcCU->getPic();
  UInt uiCurNumParts = pcPic->getNumPartitionsInCtu() >> (uiDepth<<1);
  UInt uiQNumParts   = uiCurNumParts>>2;
  const TComSPS &sps = *(pcCU->getSlice()->getSPS());
#if BT_FORCE
  UInt uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdx]];
#if JVET_C0024_QTBT
  UInt uiRPelX = uiLPelX + uiWidth - 1;
#else
  UInt uiRPelX = uiLPelX + (maxCuWidth >> uiDepth) - 1;
#endif
  UInt uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdx]];
#if JVET_C0024_QTBT
  UInt uiBPelY = uiTPelY + uiHeight - 1;
#else
  UInt uiBPelY = uiTPelY + (maxCuHeight >> uiDepth) - 1;
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
#endif
#endif
  if ((pcCU->getDepth(uiAbsZorderIdx) > uiDepth)
#if BT_FORCE
    && (
#if !BT_FORCE_INTRA_ON
      pcCU->getSlice()->isIntra() || 
#endif
      (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
    )
  {
    for ( UInt uiPartIdx = 0; uiPartIdx < 4; uiPartIdx++, uiAbsZorderIdx+=uiQNumParts )
    {
#if BT_FORCE
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdx]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdx]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
#endif
#else
      UInt uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdx]];
      UInt uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdx]];
#if BARC_Boundary
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
#endif
#endif
      {
#if JVET_C0024_QTBT
        xDeblockCU( pcCU, uiAbsZorderIdx, uiDepth+1, uiWidth>>1, uiHeight>>1, edgeDir );
#else
        xDeblockCU( pcCU, uiAbsZorderIdx, uiDepth+1, edgeDir );
#endif
      }
    }
    return;
  }

#if JVET_C0024_QTBT
  UInt uiBTDepth = pcCU->getBTDepth(uiAbsZorderIdx, uiWidth, uiHeight);
  UInt uiMinCUW = pcCU->getPic()->getMinCUWidth();
  UInt uiMinCUH = pcCU->getPic()->getMinCUHeight();

#if TT_TREE_DS1
  if (pcCU->getBTSplitModeForBTDepth(uiAbsZorderIdx, uiBTDepth)==ST_HORI_BT)
#else
  if (pcCU->getBTSplitModeForBTDepth(uiAbsZorderIdx, uiBTDepth)==1)
#endif
  {
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsZorderIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] 
        + (uiHeight>>1)/uiMinCUH*pcCU->getPic()->getNumPartInCtuWidth()];
      }
#if BT_FORCE
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdx]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdx]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDeblockCU(pcCU, uiAbsZorderIdx, uiDepth, uiWidth, uiHeight >> 1, edgeDir);
#if BT_FORCE
      }
#endif
    }
    return;
  }
#if TT_TREE_DS1
  else if (pcCU->getBTSplitModeForBTDepth(uiAbsZorderIdx, uiBTDepth)==ST_VERT_BT)
#else
  else if (pcCU->getBTSplitModeForBTDepth(uiAbsZorderIdx, uiBTDepth)==2)
#endif
  {
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsZorderIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] 
        + (uiWidth>>1)/uiMinCUW];
      }
#if BT_FORCE
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdx]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdx]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDeblockCU(pcCU, uiAbsZorderIdx, uiDepth, uiWidth >> 1, uiHeight, edgeDir);
#if BT_FORCE
      }
#endif
    }
    return;
  }
#if TT_ENC_OPT || TT_DEC_OPT //xDeblockCU
  UInt uiHeightSubPart, uiWidthSubPart, uiAbsZorderIdxSubPart;
  UInt uiSTSplitMode = pcCU->getBTSplitModeForBTDepth(uiAbsZorderIdx, uiBTDepth);
  if (uiSTSplitMode==ST_HORI_TT_L || uiSTSplitMode==ST_HORI_TT_S )
  {
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsZorderIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + (uiHeight>>2)/uiMinCUH*pcCU->getPic()->getNumPartInCtuWidth()];
        uiHeightSubPart = uiHeight>>1;
      }
      else
      {
        if(uiPartUnitIdx==2)
          uiAbsZorderIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + ((uiHeight>>2)+(uiHeight>>1))/uiMinCUH*pcCU->getPic()->getNumPartInCtuWidth()];
        else
          uiAbsZorderIdxSubPart = uiAbsZorderIdx;
        uiHeightSubPart = uiHeight>>2;
      }
#if BT_FORCE_TT
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdxSubPart]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdxSubPart]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDeblockCU(pcCU, uiAbsZorderIdxSubPart, uiDepth, uiWidth, uiHeightSubPart, edgeDir);
#if BT_FORCE_TT
      }
#endif
    }
    return;
  }
  if(uiSTSplitMode == ST_VERT_TT_L || uiSTSplitMode == ST_VERT_TT_S )
  {
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
    {
      if (uiPartUnitIdx==1)
      {
        uiAbsZorderIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + (uiWidth>>2)/uiMinCUW];
        uiWidthSubPart = uiWidth>>1;
      }
      else
      {
        if(uiPartUnitIdx==2)
          uiAbsZorderIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsZorderIdx] + ((uiWidth>>2)+(uiWidth>>1))/uiMinCUW];
        else
          uiAbsZorderIdxSubPart = uiAbsZorderIdx;
        uiWidthSubPart = uiWidth>>2;
      }
#if BT_FORCE_TT
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsZorderIdxSubPart]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsZorderIdxSubPart]];
#if BARC_BTFORCE
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
      {
#endif
        xDeblockCU(pcCU, uiAbsZorderIdxSubPart, uiDepth, uiWidthSubPart, uiHeight, edgeDir);
#if BT_FORCE_TT
      }
#endif
    }
    return;
  }
#endif
#endif

  xSetLoopfilterParam( pcCU, uiAbsZorderIdx );
#if JVET_C0024_QTBT
#if BARC_DBK
  UInt isBARC = pcCU->getBARCFlag(pcCU->getTextType());
  xSetEdgefilterCU   ( pcCU, uiAbsZorderIdx*(1+3*isBARC), uiWidth<<isBARC, uiHeight<<isBARC);
#else
  xSetEdgefilterCU   ( pcCU, uiAbsZorderIdx, uiWidth, uiHeight);
#endif
#else
  TComTURecurse tuRecurse(pcCU, uiAbsZorderIdx);
  xSetEdgefilterTU   ( tuRecurse );
  xSetEdgefilterPU   ( pcCU, uiAbsZorderIdx );
#endif

#if JVET_C0024_QTBT
 // const UInt uiPelsInPart = sps.getCTUSize() >> sps.getMaxTotalCUDepth();
#if BARC_DBK
  for (UInt j = 0; j < (edgeDir == EDGE_VER ? (uiHeight << isBARC) : 1); j += uiMinCUH)
  {
    for (UInt i = 0; i < (edgeDir == EDGE_VER ? 1 : (uiWidth << isBARC)); i += uiMinCUW)
    {
      UInt uiRaster = j / uiMinCUH * pcCU->getPic()->getNumPartInCtuWidth() + i / uiMinCUW + g_auiZscanToRaster[uiAbsZorderIdx*(1 + 3 * isBARC)];
      UInt uiPartIdx = g_auiRasterToZscan[uiRaster];
      UInt uiBSCheck = 1;
#else
  for (UInt j=0; j<(edgeDir==EDGE_VER ? uiHeight: 1); j+=uiMinCUH)
  {
    for (UInt i=0; i<(edgeDir==EDGE_VER ? 1: uiWidth); i+=uiMinCUW)
    {
      UInt uiRaster = j/uiMinCUH * pcCU->getPic()->getNumPartInCtuWidth() + i/uiMinCUW + g_auiZscanToRaster[uiAbsZorderIdx];
      UInt uiPartIdx = g_auiRasterToZscan[uiRaster];
      UInt uiBSCheck = 1;
#endif
#else
  const UInt uiPelsInPart = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth();

  for( UInt uiPartIdx = uiAbsZorderIdx; uiPartIdx < uiAbsZorderIdx + uiCurNumParts; uiPartIdx++ )
  {
    UInt uiBSCheck;
    if( uiPelsInPart == 4 )
    {
      uiBSCheck = (edgeDir == EDGE_VER && uiPartIdx%2 == 0) || (edgeDir == EDGE_HOR && (uiPartIdx-((uiPartIdx>>2)<<2))/2 == 0);
    }
    else
    {
      uiBSCheck = 1;
    }

#endif
    if ( m_aapbEdgeFilter[edgeDir][uiPartIdx] && uiBSCheck )
    {
#if BARC_DBK
      xGetBoundaryStrengthSingle ( pcCU, edgeDir, uiPartIdx/(1+3*isBARC) );
#else
      xGetBoundaryStrengthSingle ( pcCU, edgeDir, uiPartIdx );
#endif
    }
#if JVET_C0024_QTBT
    }
  }

  UInt uiSizeInPU = 1;
  UInt PartIdxIncr = 1;
#else
  }

  UInt PartIdxIncr = DEBLOCK_SMALLEST_BLOCK / uiPelsInPart ? DEBLOCK_SMALLEST_BLOCK / uiPelsInPart : 1 ;

  UInt uiSizeInPU = pcPic->getNumPartInCtuWidth()>>(uiDepth);
  const ChromaFormat chFmt=pcPic->getChromaFormat();
  const UInt shiftFactor  = edgeDir == EDGE_VER ? pcPic->getComponentScaleX(COMPONENT_Cb) : pcPic->getComponentScaleY(COMPONENT_Cb);
  const Bool bAlwaysDoChroma=chFmt==CHROMA_444;
#endif

  for ( Int iEdge = 0; iEdge < uiSizeInPU ; iEdge+=PartIdxIncr)
  {
#if JVET_C0024_QTBT
    if (isLuma(pcCU->getTextType()))
    {
#endif
    xEdgeFilterLuma     ( pcCU, uiAbsZorderIdx, uiDepth, edgeDir, iEdge );
#if JVET_C0024_QTBT
    }
    if ((isChroma(pcCU->getTextType()) || !pcCU->getSlice()->isIntra()) && pcPic->getChromaFormat()!=CHROMA_400)
    {
#else
    if ( chFmt!=CHROMA_400 && (bAlwaysDoChroma ||
                               (uiPelsInPart>DEBLOCK_SMALLEST_BLOCK) ||
                               (iEdge % ( (DEBLOCK_SMALLEST_BLOCK<<shiftFactor)/uiPelsInPart ) ) == 0
                              )
       )
#endif
    {
      xEdgeFilterChroma   ( pcCU, uiAbsZorderIdx, uiDepth, edgeDir, iEdge );
    }
#if JVET_C0024_QTBT
    }
#endif
  }
}

Void TComLoopFilter::xSetEdgefilterMultiple( TComDataCU*    pcCU,
                                             UInt           uiAbsZorderIdx,
                                             UInt           uiDepth,
                                             DeblockEdgeDir edgeDir,
                                             Int            iEdgeIdx,
                                             Bool           bValue,
                                             UInt           uiWidthInBaseUnits,
                                             UInt           uiHeightInBaseUnits,
                                             const TComRectangle *rect)
{
  if ( uiWidthInBaseUnits == 0 )
  {
    uiWidthInBaseUnits  = pcCU->getPic()->getNumPartInCtuWidth () >> uiDepth;
  }
  if ( uiHeightInBaseUnits == 0 )
  {
    uiHeightInBaseUnits = pcCU->getPic()->getNumPartInCtuHeight() >> uiDepth;
  }
  const UInt uiNumElem = edgeDir == EDGE_VER ? uiHeightInBaseUnits : uiWidthInBaseUnits;
  assert( uiNumElem > 0 );
  assert( uiWidthInBaseUnits > 0 );
  assert( uiHeightInBaseUnits > 0 );
  for( UInt ui = 0; ui < uiNumElem; ui++ )
  {
#if BARC_DBK
    UInt isBARC = pcCU->getBARCFlag(pcCU->getTextType());
    const UInt uiBsIdx = (xCalcBsIdx(pcCU, uiAbsZorderIdx, edgeDir, iEdgeIdx, ui, rect));
#else
    const UInt uiBsIdx = xCalcBsIdx( pcCU, uiAbsZorderIdx, edgeDir, iEdgeIdx, ui, rect );
#endif
    m_aapbEdgeFilter[edgeDir][uiBsIdx] = bValue;
    if (iEdgeIdx == 0)
    {
      m_aapucBS[edgeDir][uiBsIdx] = bValue;
    }
  }
}

#if JVET_C0024_QTBT
Void TComLoopFilter::xSetEdgefilterCU( TComDataCU* pcCU, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight )
{
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());

  UInt uiWidthInBaseUnits  = uiWidth / (sps.getCTUSize()  >> sps.getMaxTotalCUDepth());
  UInt uiHeightInBaseUnits = uiHeight / (sps.getCTUSize()  >> sps.getMaxTotalCUDepth());

  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, 0, EDGE_VER, 0, m_stLFCUParam.bInternalEdge, uiWidthInBaseUnits, uiHeightInBaseUnits );
  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, 0, EDGE_HOR, 0, m_stLFCUParam.bInternalEdge, uiWidthInBaseUnits, uiHeightInBaseUnits );

  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, 0, EDGE_VER, 0, m_stLFCUParam.bLeftEdge, uiWidthInBaseUnits, uiHeightInBaseUnits );
  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, 0, EDGE_HOR, 0, m_stLFCUParam.bTopEdge, uiWidthInBaseUnits, uiHeightInBaseUnits );
}
#else

Void TComLoopFilter::xSetEdgefilterTU(  TComTU &rTu )
{
  TComDataCU* pcCU  = rTu.getCU();
  UInt uiTransDepthTotal = rTu.GetTransformDepthTotal();
#if JVET_C0024_QTBT
  if( pcCU->getDepth( rTu.GetAbsPartIdxTU()) > uiTransDepthTotal )
#else
  if( pcCU->getTransformIdx( rTu.GetAbsPartIdxTU() ) + pcCU->getDepth( rTu.GetAbsPartIdxTU()) > uiTransDepthTotal )
#endif
  {
    TComTURecurse tuChild(rTu, false);
    do
    {
      xSetEdgefilterTU( tuChild );
    } while (tuChild.nextSection(rTu));
    return;
  }

  const TComRectangle &rect = rTu.getRect(COMPONENT_Y);
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());

#if JVET_C0024_QTBT
  const UInt uiWidthInBaseUnits  = rect.width  / (sps.getCTUSize()  >> sps.getMaxTotalCUDepth());
  const UInt uiHeightInBaseUnits = rect.height / (sps.getCTUSize() >> sps.getMaxTotalCUDepth());
#else
  const UInt uiWidthInBaseUnits  = rect.width  / (sps.getMaxCUWidth()  >> sps.getMaxTotalCUDepth());
  const UInt uiHeightInBaseUnits = rect.height / (sps.getMaxCUHeight() >> sps.getMaxTotalCUDepth());
#endif

  xSetEdgefilterMultiple( pcCU, rTu.GetAbsPartIdxCU(), uiTransDepthTotal, EDGE_VER, 0, m_stLFCUParam.bInternalEdge, uiWidthInBaseUnits, uiHeightInBaseUnits, &rect );
  xSetEdgefilterMultiple( pcCU, rTu.GetAbsPartIdxCU(), uiTransDepthTotal, EDGE_HOR, 0, m_stLFCUParam.bInternalEdge, uiWidthInBaseUnits, uiHeightInBaseUnits, &rect );
}

Void TComLoopFilter::xSetEdgefilterPU( TComDataCU* pcCU, UInt uiAbsZorderIdx )
{
  const UInt uiDepth = pcCU->getDepth( uiAbsZorderIdx );
  const UInt uiWidthInBaseUnits  = pcCU->getPic()->getNumPartInCtuWidth () >> uiDepth;
#if !JVET_C0024_QTBT
  const UInt uiHeightInBaseUnits = pcCU->getPic()->getNumPartInCtuHeight() >> uiDepth;
  const UInt uiHWidthInBaseUnits  = uiWidthInBaseUnits  >> 1;
  const UInt uiHHeightInBaseUnits = uiHeightInBaseUnits >> 1;
  const UInt uiQWidthInBaseUnits  = uiWidthInBaseUnits  >> 2;
  const UInt uiQHeightInBaseUnits = uiHeightInBaseUnits >> 2;
#endif

  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, 0, m_stLFCUParam.bLeftEdge );
  xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, 0, m_stLFCUParam.bTopEdge );

#if JVET_C0024_QTBT
  if( pcCU->getFRUCMgrMode( uiAbsZorderIdx ) )
  {
    assert( MIN_PU_SIZE == 4 && FRUC_MERGE_REFINE_MINBLKSIZE == 4 );
    const Int nUnits = FRUC_MERGE_REFINE_MINBLKSIZE >> 2;
    for( UInt nEdgeIdx = nUnits ; nEdgeIdx < uiWidthInBaseUnits ; nEdgeIdx += nUnits )
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, nEdgeIdx, m_stLFCUParam.bInternalEdge );
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, nEdgeIdx, m_stLFCUParam.bInternalEdge );
    }
  }
#else
  switch ( pcCU->getPartitionSize( uiAbsZorderIdx ) )
  {
    case SIZE_2Nx2N:
    {
#if VCEG_AZ07_FRUC_MERGE
      if( pcCU->getFRUCMgrMode( uiAbsZorderIdx ) )
      {
        assert( MIN_PU_SIZE == 4 && FRUC_MERGE_REFINE_MINBLKSIZE == 4 );
        const Int nUnits = FRUC_MERGE_REFINE_MINBLKSIZE >> 2;
        for( UInt nEdgeIdx = nUnits ; nEdgeIdx < uiWidthInBaseUnits ; nEdgeIdx += nUnits )
        {
          xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, nEdgeIdx, m_stLFCUParam.bInternalEdge );
          xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, nEdgeIdx, m_stLFCUParam.bInternalEdge );
        }
      }
#endif
#if JVET_B0038_AFFINE_HARMONIZATION
      if ( pcCU->getAffineFlag(uiAbsZorderIdx) )
      {
        assert( MIN_PU_SIZE == 4 && AFFINE_MIN_BLOCK_SIZE == 4 );
        const Int nUnits = AFFINE_MIN_BLOCK_SIZE >> 2;
        for( UInt nEdgeIdx = nUnits ; nEdgeIdx < uiWidthInBaseUnits ; nEdgeIdx += nUnits )
        {
          xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, nEdgeIdx, m_stLFCUParam.bInternalEdge );
          xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, nEdgeIdx, m_stLFCUParam.bInternalEdge );
        }
      }
#endif
      break;
    }
    case SIZE_2NxN:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, uiHHeightInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_Nx2N:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, uiHWidthInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_NxN:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, uiHWidthInBaseUnits, m_stLFCUParam.bInternalEdge );
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, uiHHeightInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_2NxnU:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, uiQHeightInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_2NxnD:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_HOR, uiHeightInBaseUnits - uiQHeightInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_nLx2N:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, uiQWidthInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    case SIZE_nRx2N:
    {
      xSetEdgefilterMultiple( pcCU, uiAbsZorderIdx, uiDepth, EDGE_VER, uiWidthInBaseUnits - uiQWidthInBaseUnits, m_stLFCUParam.bInternalEdge );
      break;
    }
    default:
    {
      break;
    }
  }
#endif
}


#endif
Void TComLoopFilter::xSetLoopfilterParam( TComDataCU* pcCU, UInt uiAbsZorderIdx )
{
  UInt uiX           = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[ uiAbsZorderIdx ] ];
  UInt uiY           = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[ uiAbsZorderIdx ] ];

  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;

  m_stLFCUParam.bInternalEdge = ! pcCU->getSlice()->getDeblockingFilterDisable();

  if ( (uiX == 0) || pcCU->getSlice()->getDeblockingFilterDisable() )
  {
    m_stLFCUParam.bLeftEdge = false;
  }
  else
  {
    m_stLFCUParam.bLeftEdge = true;
  }
  if ( m_stLFCUParam.bLeftEdge )
  {
    pcTempCU = pcCU->getPULeft( uiTempPartIdx, uiAbsZorderIdx, !pcCU->getSlice()->getLFCrossSliceBoundaryFlag(), !m_bLFCrossTileBoundary);

    if ( pcTempCU != NULL )
    {
      m_stLFCUParam.bLeftEdge = true;
    }
    else
    {
      m_stLFCUParam.bLeftEdge = false;
    }
  }

  if ( (uiY == 0 ) || pcCU->getSlice()->getDeblockingFilterDisable() )
  {
    m_stLFCUParam.bTopEdge = false;
  }
  else
  {
    m_stLFCUParam.bTopEdge = true;
  }
  if ( m_stLFCUParam.bTopEdge )
  {
    pcTempCU = pcCU->getPUAbove( uiTempPartIdx, uiAbsZorderIdx, !pcCU->getSlice()->getLFCrossSliceBoundaryFlag(), false, !m_bLFCrossTileBoundary);

    if ( pcTempCU != NULL )
    {
      m_stLFCUParam.bTopEdge = true;
    }
    else
    {
      m_stLFCUParam.bTopEdge = false;
    }
  }
}

Void TComLoopFilter::xGetBoundaryStrengthSingle ( TComDataCU* pCtu, DeblockEdgeDir edgeDir, UInt uiAbsPartIdx4x4BlockWithinCtu )
{
  TComSlice * const pcSlice = pCtu->getSlice();

  const Bool lfCrossSliceBoundaryFlag=pCtu->getSlice()->getLFCrossSliceBoundaryFlag();

  const UInt uiPartQ = uiAbsPartIdx4x4BlockWithinCtu;
  TComDataCU* const pcCUQ = pCtu;

  UInt uiPartP;
  TComDataCU* pcCUP;
  UInt uiBs = 0;

  //-- Calculate Block Index
  if (edgeDir == EDGE_VER)
  {
    pcCUP = pcCUQ->getPULeft(uiPartP, uiPartQ, !lfCrossSliceBoundaryFlag, !m_bLFCrossTileBoundary);
  }
  else  // (edgeDir == EDGE_HOR)
  {
    pcCUP = pcCUQ->getPUAbove(uiPartP, uiPartQ, !pCtu->getSlice()->getLFCrossSliceBoundaryFlag(), false, !m_bLFCrossTileBoundary);
  }

  //-- Set BS for Intra MB : BS = 4 or 3
  if ( pcCUP->isIntra(uiPartP) || pcCUQ->isIntra(uiPartQ) )
  {
    uiBs = 2;
  }

  //-- Set BS for not Intra MB : BS = 2 or 1 or 0
#if JVET_C0024_QTBT
  if ( isLuma(pCtu->getTextType()) && (!pcCUP->isIntra(uiPartP) && !pcCUQ->isIntra(uiPartQ)) )
#else
  if ( !pcCUP->isIntra(uiPartP) && !pcCUQ->isIntra(uiPartQ) )
#endif
  {
    UInt nsPartQ = uiPartQ;
    UInt nsPartP = uiPartP;
#if JVET_C0024_QTBT
    if ( m_aapucBS[edgeDir][uiAbsPartIdx4x4BlockWithinCtu] && (pcCUQ->getCbf( nsPartQ, COMPONENT_Y, 0) != 0 || pcCUP->getCbf( nsPartP, COMPONENT_Y, 0 ) != 0) )
#else
    if ( m_aapucBS[edgeDir][uiAbsPartIdx4x4BlockWithinCtu] && (pcCUQ->getCbf( nsPartQ, COMPONENT_Y, pcCUQ->getTransformIdx(nsPartQ)) != 0 || pcCUP->getCbf( nsPartP, COMPONENT_Y, pcCUP->getTransformIdx(nsPartP) ) != 0) )
#endif
    {
      uiBs = 1;
    }
    else
    {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      Int nThreshold = 4 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
      if (pcSlice->isInterB() || pcCUP->getSlice()->isInterB())
      {
        Int iRefIdx;
        TComPic *piRefP0, *piRefP1, *piRefQ0, *piRefQ1;
        iRefIdx = pcCUP->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartP);
        piRefP0 = (iRefIdx < 0) ? NULL : pcCUP->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx);
        iRefIdx = pcCUP->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiPartP);
        piRefP1 = (iRefIdx < 0) ? NULL : pcCUP->getSlice()->getRefPic(REF_PIC_LIST_1, iRefIdx);
        iRefIdx = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartQ);
        piRefQ0 = (iRefIdx < 0) ? NULL : pcSlice->getRefPic(REF_PIC_LIST_0, iRefIdx);
        iRefIdx = pcCUQ->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiPartQ);
        piRefQ1 = (iRefIdx < 0) ? NULL : pcSlice->getRefPic(REF_PIC_LIST_1, iRefIdx);

        TComMv pcMvP0 = pcCUP->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartP);
        TComMv pcMvP1 = pcCUP->getCUMvField(REF_PIC_LIST_1)->getMv(uiPartP);
        TComMv pcMvQ0 = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartQ);
        TComMv pcMvQ1 = pcCUQ->getCUMvField(REF_PIC_LIST_1)->getMv(uiPartQ);

        if (piRefP0 == NULL)
        {
          pcMvP0.setZero();
        }
        if (piRefP1 == NULL)
        {
          pcMvP1.setZero();
        }
        if (piRefQ0 == NULL)
        {
          pcMvQ0.setZero();
        }
        if (piRefQ1 == NULL)
        {
          pcMvQ1.setZero();
        }

        if ( ((piRefP0==piRefQ0)&&(piRefP1==piRefQ1)) || ((piRefP0==piRefQ1)&&(piRefP1==piRefQ0)) )
        {
          if ( piRefP0 != piRefP1 )   // Different L0 & L1
          {
            if ( piRefP0 == piRefQ0 )
            {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
              uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= nThreshold) ||
                       (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= nThreshold) ||
                       (abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= nThreshold) ||
                       (abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= nThreshold)) ? 1 : 0;
#else
              uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= 4) ||
                       (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= 4) ||
                       (abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= 4) ||
                       (abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
#endif
            }
            else
            {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
              uiBs  = ((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= nThreshold) ||
                       (abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= nThreshold) ||
                       (abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= nThreshold) ||
                       (abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= nThreshold)) ? 1 : 0;
#else
              uiBs  = ((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= 4) ||
                       (abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= 4) ||
                       (abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= 4) ||
                       (abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
#endif
            }
          }
          else    // Same L0 & L1
          {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
            uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= nThreshold) ||
                     (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= nThreshold) ||
                     (abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= nThreshold) ||
                     (abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= nThreshold)) &&
                    ((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= nThreshold) ||
                     (abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= nThreshold) ||
                     (abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= nThreshold) ||
                     (abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= nThreshold)) ? 1 : 0;
#else
            uiBs  = ((abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= 4) ||
                     (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= 4) ||
                     (abs(pcMvQ1.getHor() - pcMvP1.getHor()) >= 4) ||
                     (abs(pcMvQ1.getVer() - pcMvP1.getVer()) >= 4)) &&
                    ((abs(pcMvQ1.getHor() - pcMvP0.getHor()) >= 4) ||
                     (abs(pcMvQ1.getVer() - pcMvP0.getVer()) >= 4) ||
                     (abs(pcMvQ0.getHor() - pcMvP1.getHor()) >= 4) ||
                     (abs(pcMvQ0.getVer() - pcMvP1.getVer()) >= 4)) ? 1 : 0;
#endif
          }
        }
        else // for all different Ref_Idx
        {
          uiBs = 1;
        }
      }
      else  // pcSlice->isInterP()
      {
        Int iRefIdx;
        TComPic *piRefP0, *piRefQ0;
        iRefIdx = pcCUP->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartP);
        piRefP0 = (iRefIdx < 0) ? NULL : pcCUP->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx);
        iRefIdx = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartQ);
        piRefQ0 = (iRefIdx < 0) ? NULL : pcSlice->getRefPic(REF_PIC_LIST_0, iRefIdx);
        TComMv pcMvP0 = pcCUP->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartP);
        TComMv pcMvQ0 = pcCUQ->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartQ);

        if (piRefP0 == NULL)
        {
          pcMvP0.setZero();
        }
        if (piRefQ0 == NULL)
        {
          pcMvQ0.setZero();
        }

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        uiBs  = ((piRefP0 != piRefQ0) ||
                 (abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= nThreshold) ||
                 (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= nThreshold)) ? 1 : 0;
#else
        uiBs  = ((piRefP0 != piRefQ0) ||
                 (abs(pcMvQ0.getHor() - pcMvP0.getHor()) >= 4) ||
                 (abs(pcMvQ0.getVer() - pcMvP0.getVer()) >= 4)) ? 1 : 0;
#endif
      }
    }   // enf of "if( one of BCBP == 0 )"
  }   // enf of "if( not Intra )"
#if BARC_DBK
  UInt isBARC = pCtu->getBARCFlag(pCtu->getTextType());
  m_aapucBS[edgeDir][uiAbsPartIdx4x4BlockWithinCtu*(1+3*isBARC)] = uiBs;
#else
  m_aapucBS[edgeDir][uiAbsPartIdx4x4BlockWithinCtu] = uiBs;
#endif
}


Void TComLoopFilter::xEdgeFilterLuma( TComDataCU* const pcCU, const UInt uiAbsZorderIdx, const UInt uiDepth, const DeblockEdgeDir edgeDir, const Int iEdge  )
{
        TComPicYuv *pcPicYuvRec                   = pcCU->getPic()->getPicYuvRec();
#if BARC_DBK
        UInt       isBARC                         = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
        Pel        *piSrc                         = pcPicYuvRec->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiAbsZorderIdx*(1 + 3 * isBARC));
        Pel        *piTmpSrc                      = piSrc;
#else
        Pel        *piSrc                         = pcPicYuvRec->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiAbsZorderIdx );
        Pel        *piTmpSrc                      = piSrc;
#endif
  const TComSPS    &sps                           = *(pcCU->getSlice()->getSPS());
  const Bool        ppsTransquantBypassEnableFlag = pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag();
  const Int         bitDepthLuma                  = sps.getBitDepth(CHANNEL_TYPE_LUMA);
  const Bool        lfCrossSliceBoundaryFlag      = pcCU->getSlice()->getLFCrossSliceBoundaryFlag();

  Int  iStride = pcPicYuvRec->getStride(COMPONENT_Y);
  Int iQP = 0;
  Int iQP_P = 0;
  Int iQP_Q = 0;
#if JVET_C0024_QTBT 
  UInt uiNumParts;
  if (edgeDir == EDGE_VER)
  {
    uiNumParts = pcCU->getHeight(uiAbsZorderIdx) / pcCU->getPic()->getMinCUHeight();
  }
  else
  {
    uiNumParts = pcCU->getWidth(uiAbsZorderIdx) / pcCU->getPic()->getMinCUWidth();
  }
#if BARC_DBK
  uiNumParts <<= isBARC;
#endif
#else
  UInt uiNumParts = pcCU->getPic()->getNumPartInCtuWidth()>>uiDepth;
#endif

#if JVET_C0024_QTBT
  UInt  uiPelsInPart = sps.getCTUSize() >> sps.getMaxTotalCUDepth();
#else
  UInt  uiPelsInPart = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth();
#endif
  UInt  uiBsAbsIdx = 0, uiBs = 0;
  Int   iOffset, iSrcStep;

#if JVET_C0024_DF_MODIFY
  const UInt uiLCUWidthInBaseUnits = pcCU->getPic()->getNumPartInCtuWidth();
  UInt uiEdgeNumInLCUVert = g_auiZscanToRaster[uiAbsZorderIdx]%uiLCUWidthInBaseUnits + iEdge;
  UInt uiEdgeNumInLCUHor = g_auiZscanToRaster[uiAbsZorderIdx]/uiLCUWidthInBaseUnits + iEdge;
  
  if ( (uiPelsInPart < DEBLOCK_SMALLEST_BLOCK) && (( (uiEdgeNumInLCUVert%(DEBLOCK_SMALLEST_BLOCK/uiPelsInPart))&&(edgeDir==EDGE_VER) ) || ( (uiEdgeNumInLCUHor%(DEBLOCK_SMALLEST_BLOCK/uiPelsInPart))&& (edgeDir==EDGE_HOR) ) ))
  {
    return;
  }
#endif

  Bool  bPCMFilter = (sps.getUsePCM() && sps.getPCMFilterDisableFlag())? true : false;
  Bool  bPartPNoFilter = false;
  Bool  bPartQNoFilter = false;
  UInt  uiPartPIdx = 0;
  UInt  uiPartQIdx = 0;
  TComDataCU* pcCUP = pcCU;
  TComDataCU* pcCUQ = pcCU;
  Int  betaOffsetDiv2 = pcCUQ->getSlice()->getDeblockingFilterBetaOffsetDiv2();
  Int  tcOffsetDiv2 = pcCUQ->getSlice()->getDeblockingFilterTcOffsetDiv2();

  if (edgeDir == EDGE_VER)
  {
    iOffset = 1;
    iSrcStep = iStride;
    piTmpSrc += iEdge*uiPelsInPart;
  }
  else  // (edgeDir == EDGE_HOR)
  {
    iOffset = iStride;
    iSrcStep = 1;
    piTmpSrc += iEdge*uiPelsInPart*iStride;
  }

  const Int iBitdepthScale = 1 << (bitDepthLuma-8);

  for ( UInt iIdx = 0; iIdx < uiNumParts; iIdx++ )
  {
#if BARC_DBK
    uiBsAbsIdx = (xCalcBsIdx( pcCU, uiAbsZorderIdx*(1+3*isBARC), edgeDir, iEdge, iIdx));
#else
    uiBsAbsIdx = xCalcBsIdx( pcCU, uiAbsZorderIdx, edgeDir, iEdge, iIdx);
#endif
    uiBs = m_aapucBS[edgeDir][uiBsAbsIdx];
#if BARC_DBK
    uiBsAbsIdx /= (1 + 3 * isBARC);
#endif
    if ( uiBs )
    {
      iQP_Q = pcCU->getQP( uiBsAbsIdx );
#if BARC_DBK && BARC_AdpQP
      if (isBARC)
        iQP_Q += BARC_AdpQPSTEP;
#endif
      uiPartQIdx = uiBsAbsIdx;
      // Derive neighboring PU index
      if (edgeDir == EDGE_VER)
      {
        pcCUP = pcCUQ->getPULeft (uiPartPIdx, uiPartQIdx,!lfCrossSliceBoundaryFlag, !m_bLFCrossTileBoundary);
      }
      else  // (iDir == EDGE_HOR)
      {
        pcCUP = pcCUQ->getPUAbove(uiPartPIdx, uiPartQIdx,!lfCrossSliceBoundaryFlag, false, !m_bLFCrossTileBoundary);
      }

      iQP_P = pcCUP->getQP(uiPartPIdx);
#if BARC_DBK && BARC_AdpQP
      if (pcCUP->getBARCFlag(CHANNEL_TYPE_LUMA))
        iQP_P += BARC_AdpQPSTEP;
#endif
      iQP = (iQP_P + iQP_Q + 1) >> 1;

      Int iIndexTC = Clip3(0, MAX_QP+DEFAULT_INTRA_TC_OFFSET, Int(iQP + DEFAULT_INTRA_TC_OFFSET*(uiBs-1) + (tcOffsetDiv2 << 1)));
      Int iIndexB = Clip3(0, MAX_QP, iQP + (betaOffsetDiv2 << 1));

      Int iTc =  sm_tcTable[iIndexTC]*iBitdepthScale;
      Int iBeta = sm_betaTable[iIndexB]*iBitdepthScale;
      Int iSideThreshold = (iBeta+(iBeta>>1))>>3;
      Int iThrCut = iTc*10;


#if JVET_C0024_DF_MODIFY
      UInt  uiBlocksInPart = 1;
#else
      UInt  uiBlocksInPart = uiPelsInPart / 4 ? uiPelsInPart / 4 : 1;
#endif
      for (UInt iBlkIdx = 0; iBlkIdx<uiBlocksInPart; iBlkIdx ++)
      {
        Int dp0 = xCalcDP( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+0), iOffset);
        Int dq0 = xCalcDQ( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+0), iOffset);
#if JVET_C0024_DF_MODIFY
        Int dp3 = xCalcDP( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+uiPelsInPart-1), iOffset);
        Int dq3 = xCalcDQ( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+uiPelsInPart-1), iOffset);
#else
        Int dp3 = xCalcDP( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+3), iOffset);
        Int dq3 = xCalcDQ( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+3), iOffset);
#endif
        Int d0 = dp0 + dq0;
        Int d3 = dp3 + dq3;

        Int dp = dp0 + dp3;
        Int dq = dq0 + dq3;
        Int d =  d0 + d3;

        if (bPCMFilter || ppsTransquantBypassEnableFlag)
        {
          // Check if each of PUs is I_PCM with LF disabling
          bPartPNoFilter = (bPCMFilter && pcCUP->getIPCMFlag(uiPartPIdx));
          bPartQNoFilter = (bPCMFilter && pcCUQ->getIPCMFlag(uiPartQIdx));

          // check if each of PUs is lossless coded
          bPartPNoFilter = bPartPNoFilter || (pcCUP->isLosslessCoded(uiPartPIdx) );
          bPartQNoFilter = bPartQNoFilter || (pcCUQ->isLosslessCoded(uiPartQIdx) );
        }

        if (d < iBeta)
        {
          Bool bFilterP = (dp < iSideThreshold);
          Bool bFilterQ = (dq < iSideThreshold);

          Bool sw =  xUseStrongFiltering( iOffset, 2*d0, iBeta, iTc, piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+0))
#if JVET_C0024_DF_MODIFY
          && xUseStrongFiltering( iOffset, 2*d3, iBeta, iTc, piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+uiPelsInPart-1));
#else
          && xUseStrongFiltering( iOffset, 2*d3, iBeta, iTc, piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+3));
#endif

#if JVET_C0024_DF_MODIFY
          for ( Int i = 0; i < uiPelsInPart; i++)
#else
          for ( Int i = 0; i < DEBLOCK_SMALLEST_BLOCK/2; i++)
#endif
          {
            xPelFilterLuma( piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+i), iOffset, iTc, sw, bPartPNoFilter, bPartQNoFilter, iThrCut, bFilterP, bFilterQ, bitDepthLuma);
          }
        }
      }
    }
  }
}


Void TComLoopFilter::xEdgeFilterChroma( TComDataCU* const pcCU, const UInt uiAbsZorderIdx, const UInt uiDepth, const DeblockEdgeDir edgeDir, const Int iEdge )
{
        TComPicYuv *pcPicYuvRec    = pcCU->getPic()->getPicYuvRec();
        Int         iStride        = pcPicYuvRec->getStride(COMPONENT_Cb);
        Pel        *piSrcCb        = pcPicYuvRec->getAddr( COMPONENT_Cb, pcCU->getCtuRsAddr(), uiAbsZorderIdx );
        Pel        *piSrcCr        = pcPicYuvRec->getAddr( COMPONENT_Cr, pcCU->getCtuRsAddr(), uiAbsZorderIdx );
  const TComSPS    &sps            = *(pcCU->getSlice()->getSPS());
  const Int         bitDepthChroma = sps.getBitDepth(CHANNEL_TYPE_CHROMA);

#if JVET_C0024_QTBT
  const UInt  uiPelsInPartChromaH = sps.getCTUSize() >> (sps.getMaxTotalCUDepth()+pcPicYuvRec->getComponentScaleX(COMPONENT_Cb));
  const UInt  uiPelsInPartChromaV = sps.getCTUSize() >> (sps.getMaxTotalCUDepth()+pcPicYuvRec->getComponentScaleY(COMPONENT_Cb)); //seems Cr ? a bug in HM16.6?
#else
  const UInt  uiPelsInPartChromaH = sps.getMaxCUWidth() >> (sps.getMaxTotalCUDepth()+pcPicYuvRec->getComponentScaleX(COMPONENT_Cb));
  const UInt  uiPelsInPartChromaV = sps.getMaxCUHeight() >> (sps.getMaxTotalCUDepth()+pcPicYuvRec->getComponentScaleY(COMPONENT_Cb));
#endif

  Int iQP = 0;
  Int iQP_P = 0;
  Int iQP_Q = 0;

  Int   iOffset, iSrcStep;
  UInt  uiLoopLength;

  const UInt uiCtuWidthInBaseUnits = pcCU->getPic()->getNumPartInCtuWidth();

  Bool  bPCMFilter = (pcCU->getSlice()->getSPS()->getUsePCM() && pcCU->getSlice()->getSPS()->getPCMFilterDisableFlag())? true : false;
  Bool  bPartPNoFilter = false;
  Bool  bPartQNoFilter = false;
  TComDataCU* pcCUQ = pcCU;
  Int tcOffsetDiv2 = pcCU->getSlice()->getDeblockingFilterTcOffsetDiv2();

  // Vertical Position
  UInt uiEdgeNumInCtuVert = g_auiZscanToRaster[uiAbsZorderIdx]%uiCtuWidthInBaseUnits + iEdge;
  UInt uiEdgeNumInCtuHor = g_auiZscanToRaster[uiAbsZorderIdx]/uiCtuWidthInBaseUnits + iEdge;

  if ( (uiPelsInPartChromaH < DEBLOCK_SMALLEST_BLOCK) && (uiPelsInPartChromaV < DEBLOCK_SMALLEST_BLOCK) &&
       (
         ( (uiEdgeNumInCtuVert%(DEBLOCK_SMALLEST_BLOCK/uiPelsInPartChromaH)) && (edgeDir==EDGE_VER) ) ||
         ( (uiEdgeNumInCtuHor %(DEBLOCK_SMALLEST_BLOCK/uiPelsInPartChromaV)) && (edgeDir==EDGE_HOR) )
       )
     )
  {
    return;
  }


  const Bool lfCrossSliceBoundaryFlag=pcCU->getSlice()->getLFCrossSliceBoundaryFlag();

#if JVET_C0024_QTBT
  UInt uiNumParts;
  if (edgeDir == EDGE_VER)
  {
    uiNumParts = pcCU->getHeight(uiAbsZorderIdx) / pcCU->getPic()->getMinCUHeight();
  }
  else
  {
    uiNumParts = pcCU->getWidth(uiAbsZorderIdx) / pcCU->getPic()->getMinCUWidth();
  }
#else
  UInt  uiNumParts = pcCU->getPic()->getNumPartInCtuWidth()>>uiDepth;
#endif

  UInt  uiBsAbsIdx;
  UChar ucBs;

  Pel* piTmpSrcCb = piSrcCb;
  Pel* piTmpSrcCr = piSrcCr;

  if (edgeDir == EDGE_VER)
  {
    iOffset   = 1;
    iSrcStep  = iStride;
    piTmpSrcCb += iEdge*uiPelsInPartChromaH;
    piTmpSrcCr += iEdge*uiPelsInPartChromaH;
    uiLoopLength=uiPelsInPartChromaV;
  }
  else  // (edgeDir == EDGE_HOR)
  {
    iOffset   = iStride;
    iSrcStep  = 1;
    piTmpSrcCb += iEdge*iStride*uiPelsInPartChromaV;
    piTmpSrcCr += iEdge*iStride*uiPelsInPartChromaV;
    uiLoopLength=uiPelsInPartChromaH;
  }

  const Int iBitdepthScale = 1 << (pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_CHROMA)-8);

  for ( UInt iIdx = 0; iIdx < uiNumParts; iIdx++ )
  {
    uiBsAbsIdx = xCalcBsIdx( pcCU, uiAbsZorderIdx, edgeDir, iEdge, iIdx);
    ucBs = m_aapucBS[edgeDir][uiBsAbsIdx];

    if ( ucBs > 1)
    {
      iQP_Q = pcCU->getQP( uiBsAbsIdx );
      UInt  uiPartQIdx = uiBsAbsIdx;
      // Derive neighboring PU index
      TComDataCU* pcCUP;
      UInt  uiPartPIdx;

      if (edgeDir == EDGE_VER)
      {
        pcCUP = pcCUQ->getPULeft (uiPartPIdx, uiPartQIdx,!lfCrossSliceBoundaryFlag, !m_bLFCrossTileBoundary);
      }
      else  // (edgeDir == EDGE_HOR)
      {
        pcCUP = pcCUQ->getPUAbove(uiPartPIdx, uiPartQIdx,!lfCrossSliceBoundaryFlag, false, !m_bLFCrossTileBoundary);
      }

      iQP_P = pcCUP->getQP(uiPartPIdx);

      if (bPCMFilter || pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
      {
        // Check if each of PUs is I_PCM with LF disabling
        bPartPNoFilter = (bPCMFilter && pcCUP->getIPCMFlag(uiPartPIdx));
        bPartQNoFilter = (bPCMFilter && pcCUQ->getIPCMFlag(uiPartQIdx));

        // check if each of PUs is lossless coded
        bPartPNoFilter = bPartPNoFilter || (pcCUP->isLosslessCoded(uiPartPIdx));
        bPartQNoFilter = bPartQNoFilter || (pcCUQ->isLosslessCoded(uiPartQIdx));
      }

      for ( UInt chromaIdx = 0; chromaIdx < 2; chromaIdx++ )
      {
        Int chromaQPOffset  = pcCU->getSlice()->getPPS()->getQpOffset(ComponentID(chromaIdx + 1));
        Pel* piTmpSrcChroma = (chromaIdx == 0) ? piTmpSrcCb : piTmpSrcCr;

        iQP = ((iQP_P + iQP_Q + 1) >> 1) + chromaQPOffset;
        if (iQP >= chromaQPMappingTableSize)
        {
          if (pcPicYuvRec->getChromaFormat()==CHROMA_420)
          {
            iQP -=6;
          }
          else if (iQP>51)
          {
            iQP=51;
          }
        }
        else if (iQP >= 0 )
        {
          iQP = getScaledChromaQP(iQP, pcPicYuvRec->getChromaFormat());
        }

        Int iIndexTC = Clip3(0, MAX_QP+DEFAULT_INTRA_TC_OFFSET, iQP + DEFAULT_INTRA_TC_OFFSET*(ucBs - 1) + (tcOffsetDiv2 << 1));
        Int iTc =  sm_tcTable[iIndexTC]*iBitdepthScale;

        for ( UInt uiStep = 0; uiStep < uiLoopLength; uiStep++ )
        {
                    xPelFilterChroma( piTmpSrcChroma + iSrcStep*(uiStep+iIdx*uiLoopLength), iOffset, iTc , bPartPNoFilter, bPartQNoFilter, bitDepthChroma
                  #if JVET_D0033_ADAPTIVE_CLIPPING
                                      , ComponentID(chromaIdx + 1)
                  #endif
                                      );
        }
      }
    }
  }
}

/**
 - Deblocking for the luminance component with strong or weak filter
 .
 \param piSrc           pointer to picture data
 \param iOffset         offset value for picture data
 \param tc              tc value
 \param sw              decision strong/weak filter
 \param bPartPNoFilter  indicator to disable filtering on partP
 \param bPartQNoFilter  indicator to disable filtering on partQ
 \param iThrCut         threshold value for weak filter decision
 \param bFilterSecondP  decision weak filter/no filter for partP
 \param bFilterSecondQ  decision weak filter/no filter for partQ
 \param bitDepthLuma    luma bit depth
*/
__inline Void TComLoopFilter::xPelFilterLuma( Pel* piSrc, Int iOffset, Int tc, Bool sw, Bool bPartPNoFilter, Bool bPartQNoFilter, Int iThrCut, Bool bFilterSecondP, Bool bFilterSecondQ, const Int bitDepthLuma)
{
  Int delta;

  Pel m4  = piSrc[0];
  Pel m3  = piSrc[-iOffset];
  Pel m5  = piSrc[ iOffset];
  Pel m2  = piSrc[-iOffset*2];
  Pel m6  = piSrc[ iOffset*2];
  Pel m1  = piSrc[-iOffset*3];
  Pel m7  = piSrc[ iOffset*3];
  Pel m0  = piSrc[-iOffset*4];

  if (sw)
  {
    piSrc[-iOffset]   = Clip3(m3-2*tc, m3+2*tc, ((m1 + 2*m2 + 2*m3 + 2*m4 + m5 + 4) >> 3));
    piSrc[0]          = Clip3(m4-2*tc, m4+2*tc, ((m2 + 2*m3 + 2*m4 + 2*m5 + m6 + 4) >> 3));
    piSrc[-iOffset*2] = Clip3(m2-2*tc, m2+2*tc, ((m1 + m2 + m3 + m4 + 2)>>2));
    piSrc[ iOffset]   = Clip3(m5-2*tc, m5+2*tc, ((m3 + m4 + m5 + m6 + 2)>>2));
    piSrc[-iOffset*3] = Clip3(m1-2*tc, m1+2*tc, ((2*m0 + 3*m1 + m2 + m3 + m4 + 4 )>>3));
    piSrc[ iOffset*2] = Clip3(m6-2*tc, m6+2*tc, ((m3 + m4 + m5 + 3*m6 + 2*m7 +4 )>>3));
  }
  else
  {
    /* Weak filter */
    delta = (9*(m4-m3) -3*(m5-m2) + 8)>>4 ;

    if ( abs(delta) < iThrCut )
    {
      delta = Clip3(-tc, tc, delta);
#if JVET_D0033_ADAPTIVE_CLIPPING
      piSrc[-iOffset] = ClipA((m3+delta), COMPONENT_Y);
      piSrc[0] = ClipA((m4-delta),  COMPONENT_Y);
#else
      piSrc[-iOffset] = ClipBD((m3+delta), bitDepthLuma);
      piSrc[0] = ClipBD((m4-delta), bitDepthLuma);
#endif

      Int tc2 = tc>>1;
      if(bFilterSecondP)
      {
        Int delta1 = Clip3(-tc2, tc2, (( ((m1+m3+1)>>1)- m2+delta)>>1));
#if JVET_D0033_ADAPTIVE_CLIPPING
          piSrc[-iOffset*2] = ClipA((m2+delta1),  COMPONENT_Y);
#else
        piSrc[-iOffset*2] = ClipBD((m2+delta1), bitDepthLuma);
#endif
      }
      if(bFilterSecondQ)
      {
        Int delta2 = Clip3(-tc2, tc2, (( ((m6+m4+1)>>1)- m5-delta)>>1));
#if JVET_D0033_ADAPTIVE_CLIPPING
          piSrc[ iOffset] = ClipA((m5+delta2), COMPONENT_Y);
#else
        piSrc[ iOffset] = ClipBD((m5+delta2), bitDepthLuma);
#endif
      }
    }
  }

  if(bPartPNoFilter)
  {
    piSrc[-iOffset] = m3;
    piSrc[-iOffset*2] = m2;
    piSrc[-iOffset*3] = m1;
  }
  if(bPartQNoFilter)
  {
    piSrc[0] = m4;
    piSrc[ iOffset] = m5;
    piSrc[ iOffset*2] = m6;
  }
}

/**
 - Deblocking of one line/column for the chrominance component
 .
 \param piSrc           pointer to picture data
 \param iOffset         offset value for picture data
 \param tc              tc value
 \param bPartPNoFilter  indicator to disable filtering on partP
 \param bPartQNoFilter  indicator to disable filtering on partQ
 \param bitDepthChroma  chroma bit depth
 */
__inline Void TComLoopFilter::xPelFilterChroma( Pel* piSrc, Int iOffset, Int tc, Bool bPartPNoFilter, Bool bPartQNoFilter, const Int bitDepthChroma
                                                #if JVET_D0033_ADAPTIVE_CLIPPING
                                                , ComponentID compID
                                                #endif
                                                )
{
  Int delta;

  Pel m4  = piSrc[0];
  Pel m3  = piSrc[-iOffset];
  Pel m5  = piSrc[ iOffset];
  Pel m2  = piSrc[-iOffset*2];

  delta = Clip3(-tc,tc, (((( m4 - m3 ) << 2 ) + m2 - m5 + 4 ) >> 3) );
#if JVET_D0033_ADAPTIVE_CLIPPING
    piSrc[-iOffset] = ClipA((m3+delta),  compID);
    piSrc[0] = ClipA((m4-delta),  compID);
#else
  piSrc[-iOffset] = ClipBD((m3+delta), bitDepthChroma);
  piSrc[0] = ClipBD((m4-delta), bitDepthChroma);
#endif

  if(bPartPNoFilter)
  {
    piSrc[-iOffset] = m3;
  }
  if(bPartQNoFilter)
  {
    piSrc[0] = m4;
  }
}

/**
 - Decision between strong and weak filter
 .
 \param offset         offset value for picture data
 \param d               d value
 \param beta            beta value
 \param tc              tc value
 \param piSrc           pointer to picture data
 */
__inline Bool TComLoopFilter::xUseStrongFiltering( Int offset, Int d, Int beta, Int tc, Pel* piSrc)
{
  Pel m4  = piSrc[0];
  Pel m3  = piSrc[-offset];
  Pel m7  = piSrc[ offset*3];
  Pel m0  = piSrc[-offset*4];

  Int d_strong = abs(m0-m3) + abs(m7-m4);

  return ( (d_strong < (beta>>3)) && (d<(beta>>2)) && ( abs(m3-m4) < ((tc*5+1)>>1)) );
}

__inline Int TComLoopFilter::xCalcDP( Pel* piSrc, Int iOffset)
{
  return abs( piSrc[-iOffset*3] - 2*piSrc[-iOffset*2] + piSrc[-iOffset] ) ;
}

__inline Int TComLoopFilter::xCalcDQ( Pel* piSrc, Int iOffset)
{
  return abs( piSrc[0] - 2*piSrc[iOffset] + piSrc[iOffset*2] );
}
//! \}
