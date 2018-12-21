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

/** \file     TDecSlice.cpp
    \brief    slice decoder class
*/

#include "TDecSlice.h"

//! \ingroup TLibDecoder
//! \{

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TDecSlice::TDecSlice()
{
}

TDecSlice::~TDecSlice()
{
}

Void TDecSlice::create()
{
}

Void TDecSlice::destroy()
{
}

Void TDecSlice::init(TDecEntropy* pcEntropyDecoder, TDecCu* pcCuDecoder)
{
  m_pcEntropyDecoder  = pcEntropyDecoder;
  m_pcCuDecoder       = pcCuDecoder;
}

Void TDecSlice::decompressSlice(TComInputBitstream** ppcSubstreams, TComPic* pcPic, TDecSbac* pcSbacDecoder
#if ALF_HM3_REFACTOR
  , ALFParam & alfParam
#endif
  )
{
  TComSlice* pcSlice                 = pcPic->getSlice(pcPic->getCurrSliceIdx());

  const Int  startCtuTsAddr          = pcSlice->getSliceSegmentCurStartCtuTsAddr();
  const Int  startCtuRsAddr          = pcPic->getPicSym()->getCtuTsToRsAddrMap(startCtuTsAddr);
  const UInt numCtusInFrame          = pcPic->getNumberOfCtusInFrame();

  const UInt frameWidthInCtus        = pcPic->getPicSym()->getFrameWidthInCtus();
  const Bool depSliceSegmentsEnabled = pcSlice->getPPS()->getDependentSliceSegmentsEnabledFlag();
  const Bool wavefrontsEnabled       = pcSlice->getPPS()->getEntropyCodingSyncEnabledFlag();

#if BARC_REF_UP_SAMPLE_CNN
#if BARC_INTRA_LUMA_SPS
  if (pcSlice->getSPS()->getUseBarcIntraLuma())
  {
#endif
	  if (pcSlice->getSliceType() == I_SLICE)
	  {
		  pcPic->initCNNModel();
	  }
#if BARC_INTRA_LUMA_SPS
  }
#endif
#endif
  m_pcEntropyDecoder->setEntropyDecoder ( pcSbacDecoder  );
  m_pcEntropyDecoder->setBitstream      ( ppcSubstreams[0] );
  m_pcEntropyDecoder->resetEntropy      (pcSlice);

  // decoder doesn't need prediction & residual frame buffer
  pcPic->setPicYuvPred( 0 );
  pcPic->setPicYuvResi( 0 );

#if ENC_DEC_TRACE
  g_bJustDoIt = g_bEncDecTraceEnable;
#endif
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tPOC: " );
  DTRACE_CABAC_V( pcPic->getPOC() );
  DTRACE_CABAC_T( "\n" );

#if ENC_DEC_TRACE
  g_bJustDoIt = g_bEncDecTraceDisable;
#endif
#if VCEG_AZ08_KLT_COMMON
  pcPic->getPicYuvRec()->fillPicRecBoundary(pcSlice->getSPS()->getBitDepths());
#endif
  // The first CTU of the slice is the first coded substream, but the global substream number, as calculated by getSubstreamForCtuAddr may be higher.
  // This calculates the common offset for all substreams in this slice.
  const UInt subStreamOffset=pcPic->getSubstreamForCtuAddr(startCtuRsAddr, true, pcSlice);


  if (depSliceSegmentsEnabled)
  {
    // modify initial contexts with previous slice segment if this is a dependent slice.
    const UInt startTileIdx=pcPic->getPicSym()->getTileIdxMap(startCtuRsAddr);
    const TComTile *pCurrentTile=pcPic->getPicSym()->getTComTile(startTileIdx);
    const UInt firstCtuRsAddrOfTile = pCurrentTile->getFirstCtuRsAddr();

    if( pcSlice->getDependentSliceSegmentFlag() && startCtuRsAddr != firstCtuRsAddrOfTile)
    {
      if ( pCurrentTile->getTileWidthInCtus() >= 2 || !wavefrontsEnabled)
      {
        pcSbacDecoder->loadContexts(&m_lastSliceSegmentEndContextState);
      }
    }
  }

  // for every CTU in the slice segment...
#if BLOCK_BASED_BACKGROUND_REFERENCE && !DISABLE_DQP_SIGNALING
	Int iOrgSliceQP;
	if (pcSlice->getSPS()->getBackgroundRefsPresent() && pcSlice->getSliceType() != I_SLICE)
	{
		iOrgSliceQP = pcSlice->getSliceQp();
	}
#endif

#if DIMD_INTRA_PRED
  for (UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < numCtusInFrame; ctuTsAddr++)
  {
    const UInt  ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    TComDataCU* pCtu = pcPic->getCtu(ctuRsAddr);
    pCtu->initCtu(pcPic, ctuRsAddr);
  }
#endif

  Bool isLastCtuOfSliceSegment = false;
  for( UInt ctuTsAddr = startCtuTsAddr; !isLastCtuOfSliceSegment && ctuTsAddr < numCtusInFrame; ctuTsAddr++)
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    const TComTile &currentTile = *(pcPic->getPicSym()->getTComTile(pcPic->getPicSym()->getTileIdxMap(ctuRsAddr)));
    const UInt firstCtuRsAddrOfTile = currentTile.getFirstCtuRsAddr();
    const UInt tileXPosInCtus = firstCtuRsAddrOfTile % frameWidthInCtus;
    const UInt tileYPosInCtus = firstCtuRsAddrOfTile / frameWidthInCtus;
    const UInt ctuXPosInCtus  = ctuRsAddr % frameWidthInCtus;
    const UInt ctuYPosInCtus  = ctuRsAddr / frameWidthInCtus;
    const UInt uiSubStrm=pcPic->getSubstreamForCtuAddr(ctuRsAddr, true, pcSlice)-subStreamOffset;
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );
#if !DIMD_INTRA_PRED
    pCtu->initCtu(pcPic, ctuRsAddr);
#endif

    m_pcEntropyDecoder->setBitstream( ppcSubstreams[uiSubStrm] );

#if VCEG_AZ07_INIT_PREVFRAME
    if( pcSlice->getSliceType() != I_SLICE && ctuTsAddr == 0 )
    {
      pcSbacDecoder->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdx(), true, pcSlice->getCtxMapQPIdxforStore(), (pcSlice->getPOC() > pcSlice->getStatsHandle()->m_uiLastIPOC)  ); 
    }
#endif

#if ALF_HM3_REFACTOR
    if ( pcSlice->getSPS()->getUseALF() && ctuRsAddr == 0 )
    {
      m_pcEntropyDecoder->decodeAlfParam(&alfParam, pcSlice->getSPS()->getMaxTotalCUDepth()
#if FIX_TICKET12
        , pcSlice
#endif
        );
    }
#endif

    // set up CABAC contexts' state for this CTU
    if (ctuRsAddr == firstCtuRsAddrOfTile)
    {
      if (ctuTsAddr != startCtuTsAddr) // if it is the first CTU, then the entropy coder has already been reset
      {
        m_pcEntropyDecoder->resetEntropy(pcSlice);
      }
    }
    else if (ctuXPosInCtus == tileXPosInCtus && wavefrontsEnabled)
    {
      // Synchronize cabac probabilities with upper-right CTU if it's available and at the start of a line.
      if (ctuTsAddr != startCtuTsAddr) // if it is the first CTU, then the entropy coder has already been reset
      {
        m_pcEntropyDecoder->resetEntropy(pcSlice);
      }
      TComDataCU *pCtuUp = pCtu->getCtuAbove();
      if ( pCtuUp && ((ctuRsAddr%frameWidthInCtus+1) < frameWidthInCtus)  )
      {
        TComDataCU *pCtuTR = pcPic->getCtu( ctuRsAddr - frameWidthInCtus + 1 );
        if ( pCtu->CUIsFromSameSliceAndTile(pCtuTR) )
        {
          // Top-right is available, so use it.
          pcSbacDecoder->loadContexts( &m_entropyCodingSyncContextState );
        }
      }
    }

#if ENC_DEC_TRACE
    g_bJustDoIt = g_bEncDecTraceEnable;
#endif
#if CNN_BASED_INLOOP_FILTER
		if (pcSlice->getSPS()->getUseCNNLoopFilter())
		{
			CNNLFBlkParam* convNetBlkParam = pcPic->getPicSym()->getConvNetBlkParam();
			for (Int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
			{
				ComponentID compId = ComponentID(comp);
				convNetBlkParam[ctuRsAddr][compId] = false;
			}
			pcSbacDecoder->parseConvNetBlkParam(convNetBlkParam, pCtu);
		}
#endif
    if ( pcSlice->getSPS()->getUseSAO() )
    {
      SAOBlkParam& saoblkParam = (pcPic->getPicSym()->getSAOBlkParam())[ctuRsAddr];
      Bool bIsSAOSliceEnabled = false;
      Bool sliceEnabled[MAX_NUM_COMPONENT];
      for(Int comp=0; comp < MAX_NUM_COMPONENT; comp++)
      {
        ComponentID compId=ComponentID(comp);
        sliceEnabled[compId] = pcSlice->getSaoEnabledFlag(toChannelType(compId)) && (comp < pcPic->getNumberValidComponents());
        if (sliceEnabled[compId])
        {
          bIsSAOSliceEnabled=true;
        }
        saoblkParam[compId].modeIdc = SAO_MODE_OFF;
      }
      if (bIsSAOSliceEnabled)
      {
        Bool leftMergeAvail = false;
        Bool aboveMergeAvail= false;

        //merge left condition
        Int rx = (ctuRsAddr % frameWidthInCtus);
        if(rx > 0)
        {
          leftMergeAvail = pcPic->getSAOMergeAvailability(ctuRsAddr, ctuRsAddr-1);
        }
        //merge up condition
        Int ry = (ctuRsAddr / frameWidthInCtus);
        if(ry > 0)
        {
          aboveMergeAvail = pcPic->getSAOMergeAvailability(ctuRsAddr, ctuRsAddr-frameWidthInCtus);
        }

        pcSbacDecoder->parseSAOBlkParam( saoblkParam, sliceEnabled, leftMergeAvail, aboveMergeAvail, pcSlice->getSPS()->getBitDepths());
      }
    }

#if RRC_NS
		if (pcSlice->getSPS()->getUseNS() && ctuTsAddr == 0)
		{
			bool isAdaptiveNs = m_pcEntropyDecoder->decodeAdaptiveNsFlag();
			pcPic->getPicSym()->setAdaptiveNs(isAdaptiveNs);
			if (isAdaptiveNs)
			{
				const int ns_block_size = m_pcEntropyDecoder->decodeNsBlockSize();
				pcPic->getPicSym()->setNsBlockSize(ns_block_size);
#if JVET_C0024_QTBT
				const int block_width = pcPic->getSlice(0)->getSPS()->getCTUSize() >> ns_block_size;
				const int block_height = pcPic->getSlice(0)->getSPS()->getCTUSize() >> ns_block_size;
#else
				const int block_width = pcPic->getSlice(0)->getSPS()->getMaxCUWidth() >> ns_block_size;
				const int block_height = pcPic->getSlice(0)->getSPS()->getMaxCUHeight() >> ns_block_size;
#endif
				const int width_in_blocks = pcPic->getPicYuvRec()->getWidth(COMPONENT_Y) / block_width + (pcPic->getPicYuvRec()->getWidth(COMPONENT_Y) % block_width != 0);
				const int height_in_blocks = pcPic->getPicYuvRec()->getHeight(COMPONENT_Y) / block_height + (pcPic->getPicYuvRec()->getHeight(COMPONENT_Y) % block_height != 0);

				UInt* nsFlag = pcPic->getPicSym()->getNsFlag();
				for (int i = 0; i < width_in_blocks * height_in_blocks; ++i)
					nsFlag[i] = m_pcEntropyDecoder->decodeNsFlag();
			}
			else
			{
				const bool isNsWholeFrame = m_pcEntropyDecoder->decodeNsFWholeFrameFlag();
				pcPic->getPicSym()->setNsWholeFrame(isNsWholeFrame);
			}
		}
#endif

#if BT_FORCE
#if !BT_FORCE_INTRA_ON
    if (!pcSlice->isIntra())
#endif
      pcSlice->xSetMaxFQTLevel();
#endif
	m_pcCuDecoder->decodeCtu     ( pCtu, isLastCtuOfSliceSegment );
	m_pcCuDecoder->decompressCtu(pCtu);

#if BLOCK_BASED_BACKGROUND_REFERENCE
	//reset QP for filter followed
	if (pcSlice->getSPS()->getBackgroundRefsPresent() && pcSlice->getSliceType() != I_SLICE)
	{
#if !QTBT_DQP_FIX
		if (pCtu->getCodedQP() != pcSlice->getSliceQp() && pcSlice->getSPS()->getPicBg()->getSpliceIdx(iCUAddr) == 0)
		{
			if (pCtu->getCodedQP() > 0)
			{
				Int partStartIdx = 0 - pcPic->getPicSym()->getInverseCUOrderMap(iCUAddr) * pcPic->getNumPartInCU();
				Int firstElement = max<Int>(partStartIdx, 0);
				Int  numElements = pcPic->getNumPartInCU() - firstElement;
				memset(pCtu->getQP() + firstElement, pCtu->getCodedQP(), numElements * sizeof(*pCtu->getQP()));
			}
		}
#endif
		if (pCtu->getLCUBgFlag())
		{
			assert(pcSlice->getSPS()->getPicBg()->getSpliceIdx(ctuRsAddr) == 0);
			pcSlice->getSPS()->getPicBg()->setSpliceIdx(ctuRsAddr, pcSlice->getPOC());
		}
#if QTBT_DQP_FIX
		Int partStartIdx = 0 - pcPic->getPicSym()->getCtuRsToTsAddrMap(ctuRsAddr) * pcPic->getNumPartitionsInCtu();
		Int firstElement = max<Int>(partStartIdx, 0);
		Int  numElements = pcPic->getNumPartitionsInCtu() - firstElement;
		if (pCtu->getLCUBgFlag())
		{
			memset(pCtu->getQP() + firstElement, pcSlice->getSliceBgQp(), numElements * sizeof(*pCtu->getQP()));
		}
		else
		{
			memset(pCtu->getQP() + firstElement, pcSlice->getSliceQp(), numElements * sizeof(*pCtu->getQP()));
		}
#endif
#if !DISABLE_DQP_SIGNALING
		pcSlice->setSliceQp(iOrgSliceQP);
#endif
	}
#endif

#if  TRACE_FLAGS 
	gNumCtu++;
#endif

#if ENC_DEC_TRACE
    g_bJustDoIt = g_bEncDecTraceDisable;
#endif

    //Store probabilities of second CTU in line into buffer
    if ( ctuXPosInCtus == tileXPosInCtus+1 && wavefrontsEnabled)
    {
      m_entropyCodingSyncContextState.loadContexts( pcSbacDecoder );
    }

    if (isLastCtuOfSliceSegment)
    {
#if DECODER_CHECK_SUBSTREAM_AND_SLICE_TRAILING_BYTES
      pcSbacDecoder->parseRemainingBytes(false);
#endif
      if(!pcSlice->getDependentSliceSegmentFlag())
      {
        pcSlice->setSliceCurEndCtuTsAddr( ctuTsAddr+1 );
      }
      pcSlice->setSliceSegmentCurEndCtuTsAddr( ctuTsAddr+1 );
    }
    else if (  ctuXPosInCtus + 1 == tileXPosInCtus + currentTile.getTileWidthInCtus() &&
             ( ctuYPosInCtus + 1 == tileYPosInCtus + currentTile.getTileHeightInCtus() || wavefrontsEnabled)
            )
    {
      // The sub-stream/stream should be terminated after this CTU.
      // (end of slice-segment, end of tile, end of wavefront-CTU-row)
      UInt binVal;
      pcSbacDecoder->parseTerminatingBit( binVal );
      assert( binVal );
#if DECODER_CHECK_SUBSTREAM_AND_SLICE_TRAILING_BYTES
      pcSbacDecoder->parseRemainingBytes(true);
#endif
    }
#if VCEG_AZ07_INIT_PREVFRAME
    if( pcSlice->getSliceType() != I_SLICE )
    {
      UInt uiTargetCUAddr = pcPic->getFrameWidthInCtus()/2 + pcPic->getNumberOfCtusInFrame()/2;
      if( uiTargetCUAddr >= pcPic->getNumberOfCtusInFrame() )
      {
        uiTargetCUAddr = pcPic->getNumberOfCtusInFrame() - 1;
      }
      if( ctuTsAddr == uiTargetCUAddr)
      {        
        pcSbacDecoder->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdxforStore(), false ); 
      }
    }
#endif
  }
#if BARCPostUpSample
  for (UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < numCtusInFrame; ctuTsAddr++)
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    TComDataCU* pCtu = pcPic->getCtu(ctuRsAddr);

    if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
      xUpSampleRec(pCtu, CHANNEL_TYPE_LUMA);
  }
#endif
  assert(isLastCtuOfSliceSegment == true);


  if( depSliceSegmentsEnabled )
  {
    m_lastSliceSegmentEndContextState.loadContexts( pcSbacDecoder );//ctx end of dep.slice
  }
#if TRACE_YUV
  DUMP_ALL_YUV();
#endif
#if REFINE_BLOCK_BOUNDARY_NET
#if BARC_INTRA_LUMA_SPS
  if (pcSlice->getSPS()->getUseBarcIntraLuma())
  {
#endif
#if BARC_ONLY_FOR_INTRA
	  if (pcSlice->getSliceType() == I_SLICE)
	  {
		  pcPic->upSampleDScaleRecon();
	  }
#endif
#if BARC_INTRA_LUMA_SPS
  }
#endif
#endif
}

#if VCEG_AZ08_INTER_KLT
Void TDecSlice::InterpolatePic(TComPic* pcPic)
{
    m_pcCuDecoder->getPointerPrediction()->interpolatePic(pcPic);
}
#endif
#if BARCPostUpSample
Void TDecSlice::xUpSampleRec(TComDataCU*& pCtu, ChannelType chType)
{
  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

  UInt upWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt upHeight = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiWidth = upWidth >> BARC_DScaleRatio;
  UInt uiHeight = upHeight >> BARC_DScaleRatio;
  UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
  UInt starIdx, endIdx;

  Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false : true;
  Bool isLeftCUAvail = pCtu->getCUPelX() == 0 ? false : true;
  Bool isBottomAvail = pCtu->getCUPelY() + upHeight >= picHeight ? false : true;
  Bool isRightAvail = pCtu->getCUPelX() + upWidth >= picWidth ? false : true;

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
  for (UInt comp = starIdx; comp < endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
    Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);

    pCtu->upScaleComponent2x2Post(pRefScale, pSrc, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, isBottomAvail, isRightAvail, compID);
  }
}
#endif
//! \}
