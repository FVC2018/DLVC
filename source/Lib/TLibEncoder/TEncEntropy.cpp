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

/** \file     TEncEntropy.cpp
    \brief    entropy encoder class
*/

#include "TEncEntropy.h"
#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComSampleAdaptiveOffset.h"
#include "TLibCommon/TComTU.h"
#if ALF_HM3_REFACTOR
#include "TLibCommon/TComAdaptiveLoopFilter.h"
#include <cmath>
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
#include "../TLibCommon/Debug.h"
static const Bool bDebugPredEnabled = DebugOptionList::DebugPred.getInt()!=0;
#endif

//! \ingroup TLibEncoder
//! \{
#if VCEG_AZ07_BAC_ADAPT_WDOW
Void TEncEntropy::encodeCtxUpdateInfo( TComSlice* pcSlice,  TComStats* apcStats )
{
  m_pcEntropyCoderIf->codeCtxUpdateInfo( pcSlice, apcStats );
  return;
}
#endif

Void TEncEntropy::setEntropyCoder ( TEncEntropyIf* e )
{
  m_pcEntropyCoderIf = e;
}

Void TEncEntropy::encodeSliceHeader ( TComSlice* pcSlice )
{
  m_pcEntropyCoderIf->codeSliceHeader( pcSlice );
  return;
}

Void  TEncEntropy::encodeTilesWPPEntryPoint( TComSlice* pSlice )
{
  m_pcEntropyCoderIf->codeTilesWPPEntryPoint( pSlice );
}

Void TEncEntropy::encodeTerminatingBit      ( UInt uiIsLast )
{
  m_pcEntropyCoderIf->codeTerminatingBit( uiIsLast );

  return;
}

Void TEncEntropy::encodeSliceFinish()
{
  m_pcEntropyCoderIf->codeSliceFinish();
}

Void TEncEntropy::encodePPS( const TComPPS* pcPPS )
{
  m_pcEntropyCoderIf->codePPS( pcPPS );
  return;
}

Void TEncEntropy::encodeSPS( const TComSPS* pcSPS )
{
  m_pcEntropyCoderIf->codeSPS( pcSPS );
  return;
}

Void TEncEntropy::encodeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }
  m_pcEntropyCoderIf->codeCUTransquantBypassFlag( pcCU, uiAbsPartIdx );
}

Void TEncEntropy::encodeVPS( const TComVPS* pcVPS )
{
  m_pcEntropyCoderIf->codeVPS( pcVPS );
  return;
}

Void TEncEntropy::encodeSkipFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if ( pcCU->getSlice()->isIntra() )
  {
    return;
  }
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }
  m_pcEntropyCoderIf->codeSkipFlag( pcCU, uiAbsPartIdx );
}

#if BARC_UP_SAMPLE_RDO
Void TEncEntropy::encodeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD)
{
	if (bRD)
	{
		uiAbsPartIdx = 0;
	}
	m_pcEntropyCoderIf->codeUpCNNFlag(pcCU, uiAbsPartIdx);
}
#endif

#if BARC
Void TEncEntropy::encodeBARCFlag( TComDataCU* pcCU, ChannelType eCtype )
{
  m_pcEntropyCoderIf->codeBARCFlag( pcCU, eCtype );
}
#endif

#if VCEG_AZ07_IMV
Void TEncEntropy::encodeiMVFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  assert( pcCU->getSlice()->getSPS()->getIMV() );
  if ( pcCU->getSlice()->isIntra() )
  {
    return;
  }
#if JVET_C0024_QTBT
  else if( pcCU->getMergeFlag( uiAbsPartIdx ) )
#else
  else if( pcCU->getMergeFlag( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
  {
    assert( pcCU->getiMVFlag( uiAbsPartIdx ) == 0 );
    return;
  }
#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
  else if( pcCU->getAffineFlag( uiAbsPartIdx ))
#else
  else if( pcCU->getAffineFlag( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
  {
    assert( pcCU->getiMVFlag( uiAbsPartIdx ) == 0 );
    return;
  }
#endif

  m_pcEntropyCoderIf->codeiMVFlag( pcCU, uiAbsPartIdx );
}
#endif

#if COM16_C806_OBMC
Void TEncEntropy::encodeOBMCFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  if ( !pcCU->getSlice()->getSPS()->getOBMC() || !pcCU->isOBMCFlagCoded( uiAbsPartIdx ) )
  {
    return;
  }

  m_pcEntropyCoderIf->codeOBMCFlag( pcCU, uiAbsPartIdx );
}
#endif

#if VCEG_AZ06_IC
Void TEncEntropy::encodeICFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if( pcCU->isICFlagCoded( uiAbsPartIdx ) )
  {
    m_pcEntropyCoderIf->codeICFlag( pcCU, uiAbsPartIdx );
  }
}
#endif

#if VCEG_AZ05_INTRA_MPI
Void TEncEntropy::encodeMPIIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD)
{
  if (bRD)
  {
    uiAbsPartIdx = 0;
  }
  // at least one merge candidate existsput return when not encoded
  m_pcEntropyCoderIf->codeMPIIdx(pcCU, uiAbsPartIdx);
}
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
Void TEncEntropy::encodePDPCIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD)
{
  if (bRD)
  {
    uiAbsPartIdx = 0;
  }
  // at least one merge candidate existsput return when not encoded
  m_pcEntropyCoderIf->codePDPCIdx(pcCU, uiAbsPartIdx);
}
#endif

#if DIMD_INTRA_PRED
Void TEncEntropy::encodeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, Bool bRD)
{
  if (bRD)
  {
    uiAbsPartIdx = 0;
  }
  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
    return;
  }
  m_pcEntropyCoderIf->codeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth);
}
Void TEncEntropy::encodeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, Bool bRD)
{
  if (bRD)
  {
    uiAbsPartIdx = 0;
  }
  assert(pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx));
  m_pcEntropyCoderIf->codeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth);
}
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
Void TEncEntropy::encodeROTIdx( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth, Bool bRD )
{ 

if( bRD )
  {
    uiAbsPartIdx = 0;
  }
  // at least one merge candidate existsput return when not encoded
 m_pcEntropyCoderIf->codeROTIdx( pcCU, uiAbsPartIdx, uiDepth );
}

#if JVET_C0024_QTBT
Void TEncEntropy::encodeROTIdxChroma( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth, Bool bRD )
{ 

  if( bRD )
  {
    uiAbsPartIdx = 0;
  }
  // at least one merge candidate existsput return when not encoded
 m_pcEntropyCoderIf->codeROTIdxChroma( pcCU, uiAbsPartIdx, uiDepth );
}
#endif
#endif
//! encode merge flag
Void TEncEntropy::encodeMergeFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  // at least one merge candidate exists
  m_pcEntropyCoderIf->codeMergeFlag( pcCU, uiAbsPartIdx );
}

#if  BLOCK_BASED_BACKGROUND_REFERENCE
Void TEncEntropy::encodeLCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD)
{
	if (pcCU->getSlice()->isIntra())
	{
		return;
	}
	if (bRD)
	{
		uiAbsPartIdx = 0;
	}
	m_pcEntropyCoderIf->codeLCUBgFlag(pcCU, uiAbsPartIdx);
}
#endif

//! encode merge index
Void TEncEntropy::encodeMergeIndex( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
#if !JVET_C0024_QTBT
    assert( pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N );
#endif
  }
  m_pcEntropyCoderIf->codeMergeIndex( pcCU, uiAbsPartIdx );
}

#if VCEG_AZ07_FRUC_MERGE
Void TEncEntropy::encodeFRUCMgrMode( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx )
{ 
  m_pcEntropyCoderIf->codeFRUCMgrMode( pcCU, uiAbsPartIdx , uiPUIdx );
}
#endif

//! encode prediction mode
Void TEncEntropy::encodePredMode( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  if ( pcCU->getSlice()->isIntra() )
  {
    return;
  }

  m_pcEntropyCoderIf->codePredMode( pcCU, uiAbsPartIdx );
}

//! encode split flag
Void TEncEntropy::encodeSplitFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  m_pcEntropyCoderIf->codeSplitFlag( pcCU, uiAbsPartIdx, uiDepth );
}

#if JVET_C0024_QTBT
Void TEncEntropy::encodeBTSplitMode( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, Bool bRD )
{
  if (bRD)
  {
    uiAbsPartIdx=0;
  }
  m_pcEntropyCoderIf->codeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight);
}
#endif

#if COM16_C806_EMT
Void TEncEntropy::encodeEmtCuFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bCodeCuFlag )
{
  m_pcEntropyCoderIf->codeEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, bCodeCuFlag );
}
#endif

#if !JVET_C0024_QTBT
//! encode partition size
Void TEncEntropy::encodePartSize( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  m_pcEntropyCoderIf->codePartSize( pcCU, uiAbsPartIdx, uiDepth );
}
#endif


/** Encode I_PCM information.
 * \param pcCU          pointer to CU
 * \param uiAbsPartIdx  CU index
 * \param bRD           flag indicating estimation or encoding
 */
Void TEncEntropy::encodeIPCMInfo( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if(!pcCU->getSlice()->getSPS()->getUsePCM()
    || pcCU->getWidth(uiAbsPartIdx) > (1<<pcCU->getSlice()->getSPS()->getPCMLog2MaxSize())
    || pcCU->getWidth(uiAbsPartIdx) < (1<<pcCU->getSlice()->getSPS()->getPCMLog2MinSize()))
  {
    return;
  }

  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  m_pcEntropyCoderIf->codeIPCMInfo ( pcCU, uiAbsPartIdx );

}

#if JVET_C0024_QTBT
Void TEncEntropy::xEncodeTransform( Bool& bCodeDQP, Bool& codeChromaQpAdj, TComTU &rTu, ComponentID compID)
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  const Bool bChroma = isChromaEnabled(pcCU->getPic()->getChromaFormat());
#if COM16_C806_EMT
  const UInt uiDepth = rTu.GetTransformDepthTotal();
#endif
  assert(rTu.GetTransformDepthRel()==0);


  if( !pcCU->isIntra(uiAbsPartIdx) && isLuma(compID) && (!bChroma || (!pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cb, 0 ) && !pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cr, 0 ) ) ) )
  {
    assert( pcCU->getCbf( uiAbsPartIdx, COMPONENT_Y, 0 ) );
    //      printf( "saved one bin! " );
  }
  else
  {
    m_pcEntropyCoderIf->codeQtCbf( rTu, compID, true ); //luma CBF is always at the lowest level
  }

  UChar cbf = pcCU->getCbf( uiAbsPartIdx, compID , 0 );

  if (cbf)
  {
#if COM16_C806_EMT
    if (compID == COMPONENT_Y)
    {
      m_pcEntropyCoderIf->codeEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, true );
    }
#endif
#if  VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    Int dummyCbf=0;
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    Int dummyNzTs=0;
#endif
    m_pcEntropyCoderIf->codeCoeffNxN( rTu, (pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID)), compID, dummyCbf
#if JVET_C0045_C0053_NO_NSST_FOR_TS
      , dummyNzTs
#endif
      );
#else
    m_pcEntropyCoderIf->codeCoeffNxN( rTu, (pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID)), compID);
#endif
  }

  //need to add qp coding, crosscomponentpred,... JCA
}
#else

Void TEncEntropy::xEncodeTransform( Bool& bCodeDQP, Bool& codeChromaQpAdj, TComTU &rTu
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  , Int& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  , Int& iNonZeroCoeffNonTs
#endif
  )
{
//pcCU, absPartIdxCU, uiAbsPartIdx, uiDepth+1, uiTrIdx+1, quadrant,
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  const UInt numValidComponent = pcCU->getPic()->getNumberValidComponents();
  const Bool bChroma = isChromaEnabled(pcCU->getPic()->getChromaFormat());
  const UInt uiTrIdx = rTu.GetTransformDepthRel();
  const UInt uiDepth = rTu.GetTransformDepthTotal();
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  const Bool bDebugRQT=pcCU->getSlice()->getFinalized() && DebugOptionList::DebugRQT.getInt()!=0;
  if (bDebugRQT)
  {
    printf("x..codeTransform: offsetLuma=%d offsetChroma=%d absPartIdx=%d, uiDepth=%d\n width=%d, height=%d, uiTrIdx=%d, uiInnerQuadIdx=%d\n",
           rTu.getCoefficientOffset(COMPONENT_Y), rTu.getCoefficientOffset(COMPONENT_Cb), uiAbsPartIdx, uiDepth, rTu.getRect(COMPONENT_Y).width, rTu.getRect(COMPONENT_Y).height, rTu.GetTransformDepthRel(), rTu.GetSectionNumber());
  }
#endif
  const UInt uiSubdiv = pcCU->getTransformIdx( uiAbsPartIdx ) > uiTrIdx;// + pcCU->getDepth( uiAbsPartIdx ) > uiDepth;
  const UInt uiLog2TrafoSize = rTu.GetLog2LumaTrSize();


  UInt cbf[MAX_NUM_COMPONENT] = {0,0,0};
  Bool bHaveACodedBlock       = false;
  Bool bHaveACodedChromaBlock = false;

  for(UInt ch=0; ch<numValidComponent; ch++)
  {
    const ComponentID compID = ComponentID(ch);

    cbf[compID] = pcCU->getCbf( uiAbsPartIdx, compID , uiTrIdx );
    
    if (cbf[ch] != 0)
    {
      bHaveACodedBlock = true;
      if (isChroma(compID))
      {
        bHaveACodedChromaBlock = true;
      }
    }
  }

  if( pcCU->isIntra(uiAbsPartIdx) && pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_NxN && uiDepth == pcCU->getDepth(uiAbsPartIdx) )
  {
    assert( uiSubdiv );
  }
  else if( pcCU->isInter(uiAbsPartIdx) && (pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N) && uiDepth == pcCU->getDepth(uiAbsPartIdx) &&  (pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) )
  {
    if ( uiLog2TrafoSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) )
    {
      assert( uiSubdiv );
    }
    else
    {
      assert(!uiSubdiv );
    }
  }
  else if( uiLog2TrafoSize > pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() )
  {
    assert( uiSubdiv );
  }
  else if( uiLog2TrafoSize == pcCU->getSlice()->getSPS()->getQuadtreeTULog2MinSize() )
  {
    assert( !uiSubdiv );
  }
  else if( uiLog2TrafoSize == pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) )
  {
    assert( !uiSubdiv );
  }
  else
  {
    assert( uiLog2TrafoSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) );
#if COM16_C806_T64
    m_pcEntropyCoderIf->codeTransformSubdivFlag( uiSubdiv, pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrafoSize );
#else
    m_pcEntropyCoderIf->codeTransformSubdivFlag( uiSubdiv, 5 - uiLog2TrafoSize );
#endif
  }

  const UInt uiTrDepthCurr = uiDepth - pcCU->getDepth( uiAbsPartIdx );
  const Bool bFirstCbfOfCU = uiTrDepthCurr == 0;

  for(UInt ch=COMPONENT_Cb; ch<numValidComponent; ch++)
  {
    const ComponentID compID=ComponentID(ch);
    if( bFirstCbfOfCU || rTu.ProcessingAllQuadrants(compID) )
    {
      if( bFirstCbfOfCU || pcCU->getCbf( uiAbsPartIdx, compID, uiTrDepthCurr - 1 ) )
      {
        m_pcEntropyCoderIf->codeQtCbf( rTu, compID, (uiSubdiv == 0) );
      }
    }
    else
    {
      assert( pcCU->getCbf( uiAbsPartIdx, compID, uiTrDepthCurr ) == pcCU->getCbf( uiAbsPartIdx, compID, uiTrDepthCurr - 1 ) );
    }
  }

  if( uiSubdiv )
  {
    TComTURecurse tuRecurseChild(rTu, true);
#if COM16_C806_EMT
    if( bFirstCbfOfCU )
    {
      m_pcEntropyCoderIf->codeEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, true );
    }
#endif
    do
    {
      xEncodeTransform( bCodeDQP, codeChromaQpAdj, tuRecurseChild 
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
        , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
        , iNonZeroCoeffNonTs
#endif
        );
    } while (tuRecurseChild.nextSection(rTu));
  }
  else
  {
    {
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tTrIdx: abspart=" );
      DTRACE_CABAC_V( uiAbsPartIdx );
      DTRACE_CABAC_T( "\tdepth=" );
      DTRACE_CABAC_V( uiDepth );
      DTRACE_CABAC_T( "\ttrdepth=" );
      DTRACE_CABAC_V( pcCU->getTransformIdx( uiAbsPartIdx ) );
      DTRACE_CABAC_T( "\n" );
    }

    if( !pcCU->isIntra(uiAbsPartIdx) && uiDepth == pcCU->getDepth( uiAbsPartIdx ) && (!bChroma || (!pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cb, 0 ) && !pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cr, 0 ) ) ) )
    {
      assert( pcCU->getCbf( uiAbsPartIdx, COMPONENT_Y, 0 ) );
      //      printf( "saved one bin! " );
    }
    else
    {
      m_pcEntropyCoderIf->codeQtCbf( rTu, COMPONENT_Y, true ); //luma CBF is always at the lowest level
    }

#if COM16_C806_EMT
    if( bFirstCbfOfCU )
    {
      m_pcEntropyCoderIf->codeEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, cbf[COMPONENT_Y] ? true : false );
    }
#endif

#if TRACE_FLAGS  && !JVET_C0024_QTBT
	if (gIsEncodingNow)
	{
		if (pcCU->isIntra(uiAbsPartIdx))
			DUMP_TU_FLAGS_INFO("encoder", pcCU, uiAbsPartIdx, uiAbsPartIdx, rTu, rTu.GetAbsPartIdxCU());
		else
			DUMP_TU_FLAGS_INFO("encoder_tu", pcCU, uiAbsPartIdx, uiAbsPartIdx, rTu, rTu.GetAbsPartIdxCU());
	}
#endif

    if ( bHaveACodedBlock )
    {
      // dQP: only for CTU once
      if ( pcCU->getSlice()->getPPS()->getUseDQP() )
      {
        if ( bCodeDQP )
        {
          encodeQP( pcCU, rTu.GetAbsPartIdxCU() );
          bCodeDQP = false;
        }
      }

      if ( pcCU->getSlice()->getUseChromaQpAdj() )
      {
        if ( bHaveACodedChromaBlock && codeChromaQpAdj && !pcCU->getCUTransquantBypass(rTu.GetAbsPartIdxCU()) )
        {
          encodeChromaQpAdjustment( pcCU, rTu.GetAbsPartIdxCU() );
          codeChromaQpAdj = false;
        }
      }

      const UInt numValidComp=pcCU->getPic()->getNumberValidComponents();

      for(UInt ch=COMPONENT_Y; ch<numValidComp; ch++)
      {
        const ComponentID compID=ComponentID(ch);

        if (rTu.ProcessComponentSection(compID))
        {
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
          if (bDebugRQT)
          {
            printf("Call NxN for chan %d width=%d height=%d cbf=%d\n", compID, rTu.getRect(compID).width, rTu.getRect(compID).height, 1);
          }
#endif

          if (rTu.getRect(compID).width != rTu.getRect(compID).height)
          {
            //code two sub-TUs
            TComTURecurse subTUIterator(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

            do
            {
              const UChar subTUCBF = pcCU->getCbf(subTUIterator.GetAbsPartIdxTU(compID), compID, (uiTrIdx + 1));

              if (subTUCBF != 0)
              {
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
                if (bDebugRQT)
                {
                  printf("Call NxN for chan %d width=%d height=%d cbf=%d\n", compID, subTUIterator.getRect(compID).width, subTUIterator.getRect(compID).height, 1);
                }
#endif
                m_pcEntropyCoderIf->codeCoeffNxN( subTUIterator, (pcCU->getCoeff(compID) + subTUIterator.getCoefficientOffset(compID)), compID 
#if  VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                  , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
                  , iNonZeroCoeffNonTs
#endif
                  );
              }
#if COM16_C806_EMT
              else if ( isLuma(compID) && cbf[COMPONENT_Y] != 0 )
              {
                pcCU->setEmtTuIdxSubParts( DCT2_EMT, uiAbsPartIdx, uiDepth );
              }
#endif
            }
            while (subTUIterator.nextSection(rTu));
          }
          else
          {
            if (isChroma(compID) && (cbf[COMPONENT_Y] != 0))
            {
              m_pcEntropyCoderIf->codeCrossComponentPrediction( rTu, compID );
            }

            if (cbf[compID] != 0)
            {
              m_pcEntropyCoderIf->codeCoeffNxN( rTu, (pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID)), compID
#if  VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
                , iNonZeroCoeffNonTs
#endif
                );
            }
#if COM16_C806_EMT
            else if ( isLuma(compID) && cbf[COMPONENT_Y] != 0 )
            {
              pcCU->setEmtTuIdxSubParts( DCT2_EMT, uiAbsPartIdx, uiDepth );
            }
#endif
          }
        }
      }
    }
  }
}


#endif //#if JVET_C0024_QTBT
//! encode intra direction for luma
Void TEncEntropy::encodeIntraDirModeLuma  ( TComDataCU* pcCU, UInt absPartIdx, Bool isMultiplePU 
#if VCEG_AZ07_INTRA_65ANG_MODES
                                           , Int* piModes, Int  iAboveLeftCase
#endif
                                           )
{
#if DIMD_INTRA_PRED
  if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, absPartIdx))
  {
    return;
  }
#endif
  m_pcEntropyCoderIf->codeIntraDirLumaAng( pcCU, absPartIdx , isMultiplePU
#if VCEG_AZ07_INTRA_65ANG_MODES
    , piModes, iAboveLeftCase
#endif
    );
}


//! encode intra direction for chroma
Void TEncEntropy::encodeIntraDirModeChroma( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
#if DIMD_INTRA_PRED
	if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
	{
		return;
	}
#endif
  m_pcEntropyCoderIf->codeIntraDirChroma( pcCU, uiAbsPartIdx );  

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  if (bDebugPredEnabled && pcCU->getSlice()->getFinalized())
  {
    UInt cdir=pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiAbsPartIdx);
    if (cdir==36)
    {
      cdir=pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
    }
    printf("coding chroma Intra dir: %d, uiAbsPartIdx: %d, luma dir: %d\n", cdir, uiAbsPartIdx, pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx));
  }
#endif
}


Void TEncEntropy::encodePredInfo( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if( pcCU->isIntra( uiAbsPartIdx ) )                                 // If it is Intra mode, encode intra prediction mode.
  {
#if JVET_C0024_QTBT
    if (isLuma(pcCU->getTextType()))
    {
#endif
    encodeIntraDirModeLuma  ( pcCU, uiAbsPartIdx,true );
#if JVET_C0024_QTBT
    }
    if (pcCU->getPic()->getChromaFormat()!=CHROMA_400
      && (isChroma(pcCU->getTextType()) || !pcCU->getSlice()->isIntra()))
#else
    if (pcCU->getPic()->getChromaFormat()!=CHROMA_400)
#endif
    {
      encodeIntraDirModeChroma( pcCU, uiAbsPartIdx );

#if !JVET_C0024_QTBT
      if (enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) && pcCU->getPartitionSize( uiAbsPartIdx )==SIZE_NxN)
      {
        UInt uiPartOffset = ( pcCU->getPic()->getNumPartitionsInCtu() >> ( pcCU->getDepth(uiAbsPartIdx) << 1 ) ) >> 2;
        encodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset   );
        encodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset*2 );
        encodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset*3 );
      }
#endif
    }
#if RRC_INTRA_SHARPENING
#if JVET_C0024_QTBT
    if (isLuma(pcCU->getTextType()))
#endif
    {
      encodeIntraSharpFlag(pcCU, uiAbsPartIdx);
    }
#endif
  }
  else                                                                // if it is Inter mode, encode motion vector and reference index
  {
    encodePUWise( pcCU, uiAbsPartIdx );
  }
}

Void TEncEntropy::encodeCrossComponentPrediction( TComTU &rTu, ComponentID compID )
{
  m_pcEntropyCoderIf->codeCrossComponentPrediction( rTu, compID );
}

#if JVET_D0123_ME_CTX_LUT_BITS
Void TEncEntropy::encodePuMotionInfo(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPartIdx){ 
#if JVET_C0024_QTBT
#if COM16_C1016_AFFINE
  PartSize ePartSize = SIZE_2Nx2N;
#endif
#else
  PartSize ePartSize = pcCU->getPartitionSize( uiAbsPartIdx );
#endif
#if VCEG_AZ07_IMV
  Bool bNonZeroMvd = false;
#endif
    encodeMergeFlag( pcCU, uiAbsPartIdx );
    if ( pcCU->getMergeFlag( uiAbsPartIdx ) )
    {
#if VCEG_AZ07_FRUC_MERGE
      encodeFRUCMgrMode( pcCU, uiAbsPartIdx , uiPartIdx );
      if( !pcCU->getFRUCMgrMode( uiAbsPartIdx ) )
#endif
#if COM16_C1016_AFFINE
      {
        if ( ePartSize == SIZE_2Nx2N && pcCU->isAffineMrgFlagCoded(uiAbsPartIdx, uiPartIdx) )
        {
          encodeAffineFlag( pcCU, uiAbsPartIdx, uiPartIdx );
        }

        if ( !pcCU->isAffine( uiAbsPartIdx ) )
        {
          encodeMergeIndex( pcCU, uiAbsPartIdx );
        }
      }
#else
      encodeMergeIndex( pcCU, uiAbsPartIdx );
#endif
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
      if (bDebugPred)
      {
        std::cout << "Coded merge flag, CU absPartIdx: " << uiAbsPartIdx << " PU(" << uiPartIdx << ") absPartIdx: " << uiSubPartIdx;
        std::cout << " merge index: " << (UInt)pcCU->getMergeIndex(uiSubPartIdx) << std::endl;
      }
#endif
    }
    else
    {
      encodeInterDirPU( pcCU, uiAbsPartIdx );
#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
      if ( pcCU->getWidth(uiAbsPartIdx) > 8 && pcCU->getHeight(uiAbsPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#else
      if ( pcCU->getWidth(uiAbsPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#endif
      {
        encodeAffineFlag( pcCU, uiAbsPartIdx, uiPartIdx );
      }
#endif

      for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
      {
        if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
        {
          encodeRefFrmIdxPU ( pcCU, uiAbsPartIdx, RefPicList( uiRefListIdx ) );
          encodeMvdPU       ( pcCU, uiAbsPartIdx, RefPicList( uiRefListIdx ) );

#if VCEG_AZ07_IMV
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiAbsPartIdx ).getHor() != 0 );
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiAbsPartIdx ).getVer() != 0 );
#endif
          encodeMVPIdxPU    ( pcCU, uiAbsPartIdx, RefPicList( uiRefListIdx ) );
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
          if (bDebugPred)
          {
            std::cout << "refListIdx: " << uiRefListIdx << std::endl;
            std::cout << "MVD horizontal: " << pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd( uiAbsPartIdx ).getHor() << std::endl;
            std::cout << "MVD vertical:   " << pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd( uiAbsPartIdx ).getVer() << std::endl;
            std::cout << "MVPIdxPU: " << pcCU->getMVPIdx(RefPicList( uiRefListIdx ), uiSubPartIdx) << std::endl;
            std::cout << "InterDir: " << (UInt)pcCU->getInterDir(uiSubPartIdx) << std::endl;
          }
#endif
        }
      }
    } 

#if VCEG_AZ07_IMV
  if( bNonZeroMvd && pcCU->getSlice()->getSPS()->getIMV() )
  {
    encodeiMVFlag( pcCU , uiAbsPartIdx );
  }
  if( !bNonZeroMvd )
  {
    assert( pcCU->getiMVFlag( uiAbsPartIdx ) == 0 );
  }
#endif
  return;
}
#endif
//! encode motion information for every PU block
Void TEncEntropy::encodePUWise( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  const Bool bDebugPred = bDebugPredEnabled && pcCU->getSlice()->getFinalized();
#endif

#if JVET_C0024_QTBT
#if COM16_C1016_AFFINE
  PartSize ePartSize = SIZE_2Nx2N;
#endif
  UInt uiNumPU = 1;
  UInt uiPUOffset = 0;
#else
  PartSize ePartSize = pcCU->getPartitionSize( uiAbsPartIdx );
  UInt uiNumPU = ( ePartSize == SIZE_2Nx2N ? 1 : ( ePartSize == SIZE_NxN ? 4 : 2 ) );
  UInt uiDepth = pcCU->getDepth( uiAbsPartIdx );
  UInt uiPUOffset = ( g_auiPUOffset[UInt( ePartSize )] << ( ( pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth ) << 1 ) ) >> 4;
#endif
#if VCEG_AZ07_IMV
  Bool bNonZeroMvd = false;
#endif

  for ( UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset )
  {
    encodeMergeFlag( pcCU, uiSubPartIdx );
    if ( pcCU->getMergeFlag( uiSubPartIdx ) )
    {
#if VCEG_AZ07_FRUC_MERGE
      encodeFRUCMgrMode( pcCU, uiSubPartIdx , uiPartIdx );
      if( !pcCU->getFRUCMgrMode( uiSubPartIdx ) )
#endif
#if MRG_CONSIDER_COMPLEX
      {
        if ((pcCU->getWidth(uiSubPartIdx) >= 16) || (pcCU->getHeight(uiSubPartIdx) >= 16))
        {
          encodeMergeComplexFlag(pcCU, uiSubPartIdx);
          if (pcCU->getMergeComplexFlag(uiSubPartIdx))
          {
            encodeMergeComplexIndex(pcCU, uiSubPartIdx);
          }
        }

        if (!pcCU->getMergeComplexFlag(uiSubPartIdx))
#endif
#if COM16_C1016_AFFINE
      {
        if ( ePartSize == SIZE_2Nx2N && pcCU->isAffineMrgFlagCoded(uiSubPartIdx, uiPartIdx) )
        {
          encodeAffineFlag( pcCU, uiSubPartIdx, uiPartIdx );
        }

        if ( !pcCU->isAffine( uiSubPartIdx ) )
        {
          encodeMergeIndex( pcCU, uiSubPartIdx );
        }
      }
#else
      encodeMergeIndex( pcCU, uiSubPartIdx );
#endif
#if MRG_CONSIDER_COMPLEX
      }
#endif
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
      if (bDebugPred)
      {
        std::cout << "Coded merge flag, CU absPartIdx: " << uiAbsPartIdx << " PU(" << uiPartIdx << ") absPartIdx: " << uiSubPartIdx;
        std::cout << " merge index: " << (UInt)pcCU->getMergeIndex(uiSubPartIdx) << std::endl;
      }
#endif
    }
    else
    {
      encodeInterDirPU( pcCU, uiSubPartIdx );
#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
      if ( pcCU->getWidth(uiSubPartIdx) > 8 && pcCU->getHeight(uiSubPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#else
      if ( pcCU->getWidth(uiSubPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#endif
      {
        encodeAffineFlag( pcCU, uiSubPartIdx, uiPartIdx );
#if JVECT_C0062_AFFINE_SIX_PARAM
        if (pcCU->isAffine(uiSubPartIdx))
          encodeAffineParamFlag(pcCU, uiSubPartIdx, uiPartIdx);
#endif
      }
#endif

      for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
      {
        if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
        {
          encodeRefFrmIdxPU ( pcCU, uiSubPartIdx, RefPicList( uiRefListIdx ) );
          encodeMvdPU       ( pcCU, uiSubPartIdx, RefPicList( uiRefListIdx ) );

#if VCEG_AZ07_IMV
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiSubPartIdx ).getHor() != 0 );
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiSubPartIdx ).getVer() != 0 );
#endif
          encodeMVPIdxPU    ( pcCU, uiSubPartIdx, RefPicList( uiRefListIdx ) );
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
          if (bDebugPred)
          {
            std::cout << "refListIdx: " << uiRefListIdx << std::endl;
            std::cout << "MVD horizontal: " << pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd( uiAbsPartIdx ).getHor() << std::endl;
            std::cout << "MVD vertical:   " << pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd( uiAbsPartIdx ).getVer() << std::endl;
            std::cout << "MVPIdxPU: " << pcCU->getMVPIdx(RefPicList( uiRefListIdx ), uiSubPartIdx) << std::endl;
            std::cout << "InterDir: " << (UInt)pcCU->getInterDir(uiSubPartIdx) << std::endl;
          }
#endif
        }
      }
    }
  }

#if VCEG_AZ07_IMV
  if( bNonZeroMvd && pcCU->getSlice()->getSPS()->getIMV() )
  {
    encodeiMVFlag( pcCU , uiAbsPartIdx );
  }
  if( !bNonZeroMvd )
  {
    assert( pcCU->getiMVFlag( uiAbsPartIdx ) == 0 );
  }
#endif
  return;
}

Void TEncEntropy::encodeInterDirPU( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if ( !pcCU->getSlice()->isInterB() )
  {
    return;
  }

  m_pcEntropyCoderIf->codeInterDir( pcCU, uiAbsPartIdx );

  return;
}

//! encode reference frame index for a PU block
Void TEncEntropy::encodeRefFrmIdxPU( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  assert( pcCU->isInter( uiAbsPartIdx ) );

  if ( ( pcCU->getSlice()->getNumRefIdx( eRefList ) == 1 ) )
  {
    return;
  }

  if ( pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList ) )
  {
    m_pcEntropyCoderIf->codeRefFrmIdx( pcCU, uiAbsPartIdx, eRefList );
  }

  return;
}

//! encode motion vector difference for a PU block
Void TEncEntropy::encodeMvdPU( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  assert( pcCU->isInter( uiAbsPartIdx ) );

#if COM16_C1016_AFFINE
#if JVECT_C0062_AFFINE_SIX_PARAM
  Int width = pcCU->getWidth(uiAbsPartIdx);
  Int height = pcCU->getHeight(uiAbsPartIdx);
  Bool bRecShape = (height / width >= 2);
#endif

  if ( pcCU->isAffine(uiAbsPartIdx) )
  {
    if ( pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList ) )
    {
      // derive LT, RT
      UInt uiPartIdxLT, uiPartIdxRT, uiAbsIndexInLCU;
      uiAbsIndexInLCU = pcCU->getZorderIdxInCtu();
      pcCU->deriveLeftRightTopIdxGeneral( uiAbsPartIdx, 0, uiPartIdxLT, uiPartIdxRT );

#if JVECT_C0062_AFFINE_SIX_PARAM
      UInt uiPartIdxLB;
      pcCU->deriveLeftBottomIdxGeneral(uiAbsPartIdx, 0, uiPartIdxLB);

      if (pcCU->getAffineParamFlag(uiAbsPartIdx))
      {
#if TT_TREE_AFFINE // deriveZIdxInCU
        m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLT), eRefList );
        m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxRT), eRefList );
        m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLB), eRefList );
#else
        m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxLT - uiAbsIndexInLCU, eRefList);
        m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxRT - uiAbsIndexInLCU, eRefList);
        m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxLB - uiAbsIndexInLCU, eRefList);
#endif
      }
      else
      {
        if (bRecShape)
        {
#if TT_TREE_AFFINE // deriveZIdxInCU
          m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLT), eRefList );
          m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLB), eRefList );
#else
          m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxLT - uiAbsIndexInLCU, eRefList);
          m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxLB - uiAbsIndexInLCU, eRefList);
#endif
        }
        else
        {
#if TT_TREE_AFFINE // deriveZIdxInCU
          m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLT), eRefList );
          m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxRT), eRefList );
#else
          m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxLT - uiAbsIndexInLCU, eRefList);
          m_pcEntropyCoderIf->codeMvd(pcCU, uiPartIdxRT - uiAbsIndexInLCU, eRefList);
#endif
        }
      }
#else
#if TT_TREE_AFFINE // deriveZIdxInCU
      m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxLT), eRefList );
      m_pcEntropyCoderIf->codeMvd( pcCU, pcCU->deriveZIdxInCU(uiAbsIndexInLCU, uiPartIdxRT), eRefList );
#else
      m_pcEntropyCoderIf->codeMvd( pcCU, uiPartIdxLT - uiAbsIndexInLCU, eRefList );
      m_pcEntropyCoderIf->codeMvd( pcCU, uiPartIdxRT - uiAbsIndexInLCU, eRefList );
#endif
#endif
    }
    return;
  }
#endif

  if ( pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList ) )
  {
    m_pcEntropyCoderIf->codeMvd( pcCU, uiAbsPartIdx, eRefList );
  }
  return;
}

Void TEncEntropy::encodeMVPIdxPU( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  if ( (pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList )) )
  {
    m_pcEntropyCoderIf->codeMVPIdx( pcCU, uiAbsPartIdx, eRefList );
  }

  return;
}

Void TEncEntropy::encodeQtCbf( TComTU &rTu, const ComponentID compID, const Bool lowestLevel )
{
  m_pcEntropyCoderIf->codeQtCbf( rTu, compID, lowestLevel );
}

Void TEncEntropy::encodeTransformSubdivFlag( UInt uiSymbol, UInt uiCtx )
{
  m_pcEntropyCoderIf->codeTransformSubdivFlag( uiSymbol, uiCtx );
}

Void TEncEntropy::encodeQtRootCbf( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  m_pcEntropyCoderIf->codeQtRootCbf( pcCU, uiAbsPartIdx );
}

Void TEncEntropy::encodeQtCbfZero( TComTU &rTu, const ChannelType chType )
{
  m_pcEntropyCoderIf->codeQtCbfZero( rTu, chType );
}

Void TEncEntropy::encodeQtRootCbfZero( )
{
  m_pcEntropyCoderIf->codeQtRootCbfZero( );
}

// dQP
Void TEncEntropy::encodeQP( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
  {
    uiAbsPartIdx = 0;
  }

  if ( pcCU->getSlice()->getPPS()->getUseDQP() )
  {
#if BLOCK_BASED_BACKGROUND_REFERENCE && !DISABLE_DQP_SIGNALING
		if (!pcCU->getSlice()->getSPS()->getBackgroundRefsPresent())
		{
			m_pcEntropyCoderIf->codeDeltaQP(pcCU, uiAbsPartIdx);
		}
#else
    m_pcEntropyCoderIf->codeDeltaQP( pcCU, uiAbsPartIdx );
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
		pcCU->setCodedQP(pcCU->getQP(uiAbsPartIdx));
#endif
  }
}

//! encode chroma qp adjustment
Void TEncEntropy::encodeChromaQpAdjustment( TComDataCU* cu, UInt absPartIdx, Bool inRd )
{
  if( inRd )
  {
    absPartIdx = 0;
  }

  m_pcEntropyCoderIf->codeChromaQpAdjustment( cu, absPartIdx );
}

// texture

//! encode coefficients
Void TEncEntropy::encodeCoeff( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool& bCodeDQP, Bool& codeChromaQpAdj
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  , Int& bNonZeroCoeff 
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  , Int& iNonZeroCoeffNonTs
#endif
  )
{
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  const Bool bDebugRQT=pcCU->getSlice()->getFinalized() && DebugOptionList::DebugRQT.getInt()!=0;
#endif

  if( pcCU->isIntra(uiAbsPartIdx) )
  {
    if (false)
    {
      DTRACE_CABAC_VL( g_nSymbolCounter++ )
      DTRACE_CABAC_T( "\tdecodeTransformIdx()\tCUDepth=" )
      DTRACE_CABAC_V( uiDepth )
      DTRACE_CABAC_T( "\n" )
    }
  }
  else
  {
#if JVET_C0024_QTBT
    if( !pcCU->getMergeFlag( uiAbsPartIdx) )
#else
    if( !(pcCU->getMergeFlag( uiAbsPartIdx ) && pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N ) )
#endif
    {
      m_pcEntropyCoderIf->codeQtRootCbf( pcCU, uiAbsPartIdx );
    }
    if ( !pcCU->getQtRootCbf( uiAbsPartIdx ) )
    {
      return;
    }
  }

  TComTURecurse tuRecurse(pcCU, uiAbsPartIdx, uiDepth);
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  if (bDebugRQT)
  {
    printf("..codeCoeff: uiAbsPartIdx=%d, PU format=%d, 2Nx2N=%d, NxN=%d\n", uiAbsPartIdx, pcCU->getPartitionSize(uiAbsPartIdx), SIZE_2Nx2N, SIZE_NxN);
  }
#endif

#if JVET_C0024_QTBT
#if JVET_C0024_DELTA_QP_FIX
  Bool validCbf       = false;
  Bool validChromaCbf = false;
#endif
  if (isChroma(pcCU->getTextType()) || !pcCU->getSlice()->isIntra())
  {
    const UInt numValidComponent = pcCU->getPic()->getNumberValidComponents();
    for(UInt ch=COMPONENT_Cb; ch<numValidComponent; ch++)
    {
      const ComponentID compID=ComponentID(ch);
      xEncodeTransform( bCodeDQP, codeChromaQpAdj, tuRecurse, compID);
#if JVET_C0024_DELTA_QP_FIX
      if( pcCU->getCbf( uiAbsPartIdx, compID ) )
      {
        validCbf = true;
        validChromaCbf = true;
      }
#endif
    }
  }
  if (isLuma(pcCU->getTextType()))
  {
    xEncodeTransform( bCodeDQP, codeChromaQpAdj, tuRecurse, COMPONENT_Y);
#if JVET_C0024_DELTA_QP_FIX
    if( pcCU->getCbf( uiAbsPartIdx, COMPONENT_Y ) )
    {
      validCbf = true;
    }
#endif
  }

#if JVET_C0024_DELTA_QP_FIX
  if( pcCU->getSlice()->getPPS()->getUseDQP() && validCbf && bCodeDQP )
  {
    encodeQP( pcCU, uiAbsPartIdx );
    bCodeDQP = false;
  }

  if ( pcCU->getSlice()->getUseChromaQpAdj() )
  {
    if ( validChromaCbf && codeChromaQpAdj && !pcCU->getCUTransquantBypass(uiAbsPartIdx) )
    {
      encodeChromaQpAdjustment( pcCU, uiAbsPartIdx );
      codeChromaQpAdj = false;
    }
  }
#endif

#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
#if !QTBT_NSST
  Int  bCbfCU = bNonZeroCoeff = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y)
    : (pcCU->getCbf(uiAbsPartIdx, COMPONENT_Cb) || pcCU->getCbf(uiAbsPartIdx, COMPONENT_Cr) )): pcCU->getQtRootCbf(uiAbsPartIdx);
#else
  bNonZeroCoeff = 0;
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  iNonZeroCoeffNonTs = 0;
#endif
  const UInt uiFirstComp = isLuma(pcCU->getTextType()) ? COMPONENT_Y : COMPONENT_Cb;
  const UInt uiLastComp  = isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() ? COMPONENT_Y : pcCU->getPic()->getNumberValidComponents()-1;
  const Int iNonZeroCoeffThr = isLuma(pcCU->getTextType()) ? NSST_SIG_NZ_LUMA + (pcCU->getSlice()->isIntra() ? 0 : NSST_SIG_NZ_CHROMA) : NSST_SIG_NZ_CHROMA;
  for(UInt ch= uiFirstComp; ch<=uiLastComp; ch++)
  {
    const ComponentID compID   = ComponentID(ch);
    TCoeff *pcCoef             = pcCU->getCoeff(compID) + tuRecurse.getCoefficientOffset(compID);
    const TComRectangle &rRect = tuRecurse.getRect(compID);
    const UInt uiWidth         = rRect.width;
    const UInt uiHeight        = rRect.height;
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    Int iNumNonZeros = countNonZeroCoeffs(pcCoef, uiWidth * uiHeight);
    bNonZeroCoeff += iNumNonZeros;
    if( !pcCU->getTransformSkip( uiAbsPartIdx, compID) )
    {
      iNonZeroCoeffNonTs += iNumNonZeros;
    }
#else
    bNonZeroCoeff += countNonZeroCoeffs(pcCoef, uiWidth * uiHeight);
#endif
  }
#endif
#endif
#else
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  Int  bCbfCU = false;
#endif   
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  iNonZeroCoeffNonTs = 0;
#endif

  xEncodeTransform( bCodeDQP, codeChromaQpAdj, tuRecurse
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  bNonZeroCoeff = bCbfCU;
#endif
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  /// put conditions if sometimes flag is not encoded
#if JVET_C0024_QTBT
  if (
#if QTBT_NSST
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    iNonZeroCoeffNonTs > iNonZeroCoeffThr
#else
    bNonZeroCoeff > iNonZeroCoeffThr
#endif
#else
    bCbfCU 
#endif
    && isLuma(pcCU->getTextType()))
  {
    encodeROTIdx( pcCU, uiAbsPartIdx, uiDepth );
  }
  else if(
#if QTBT_NSST
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    iNonZeroCoeffNonTs > iNonZeroCoeffThr
#else
    bNonZeroCoeff > iNonZeroCoeffThr
#endif
#else
    bCbfCU 
#endif
     && isChroma(pcCU->getTextType()) && pcCU->getWidth(uiAbsPartIdx) >= 8 && pcCU->getHeight(uiAbsPartIdx) >= 8 )
  {
    assert( pcCU->getSlice()->isIntra() );
    encodeROTIdxChroma( pcCU, uiAbsPartIdx, uiDepth );
  }
#else
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  if ( iNonZeroCoeffNonTs  )
#else
  if (bCbfCU  )
#endif
    encodeROTIdx( pcCU, uiAbsPartIdx, uiDepth );
#endif 
#endif 
}

Void TEncEntropy::encodeCoeffNxN( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID)
{
  TComDataCU *pcCU = rTu.getCU();

#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  Int  bCbfCU = false;
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  Int iNonZeroCoeffNonTs = 0;
#endif

  if (pcCU->getCbf(rTu.GetAbsPartIdxTU(), compID, rTu.GetTransformDepthRel()) != 0)
  {
#if !JVET_C0024_QTBT
    if (rTu.getRect(compID).width != rTu.getRect(compID).height)
    {
      //code two sub-TUs
      TComTURecurse subTUIterator(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

      const UInt subTUSize = subTUIterator.getRect(compID).width * subTUIterator.getRect(compID).height;

      do
      {
        const UChar subTUCBF = pcCU->getCbf(subTUIterator.GetAbsPartIdxTU(compID), compID, (subTUIterator.GetTransformDepthRel() + 1));

        if (subTUCBF != 0)
        {
          m_pcEntropyCoderIf->codeCoeffNxN( subTUIterator, (pcCoef + (subTUIterator.GetSectionNumber() * subTUSize)), compID
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
            , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
            , iNonZeroCoeffNonTs
#endif
            );
        }
      }
      while (subTUIterator.nextSection(rTu));
    }
    else
#endif
    {
      m_pcEntropyCoderIf->codeCoeffNxN(rTu, pcCoef, compID
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
        , bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
        , iNonZeroCoeffNonTs
#endif
        );
    }
  }
}

#if JVET_D0123_ME_CTX_LUT_BITS
Void TEncEntropy::estimatePuMeBit (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac )
{
  m_pcEntropyCoderIf->estPuMeBit ( pcEstPuMeBitsSbac );
}
#endif

Void TEncEntropy::estimateBit (estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, const ChannelType chType
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
  , UInt uiScanIdx
#endif
  )
{
#if JVET_C0024_QTBT
  if (width==2 || height==2)
  {
      return;   //don't use RDOQ for 2xn;
  }

  m_pcEntropyCoderIf->estBit ( pcEstBitsSbac, width, height, chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , uiScanIdx
#endif
    );
#else
  const UInt heightAtEntropyCoding = (width != height) ? (height >> 1) : height;

  m_pcEntropyCoderIf->estBit ( pcEstBitsSbac, width, heightAtEntropyCoding, chType );
#endif
}

Int TEncEntropy::countNonZeroCoeffs( TCoeff* pcCoef, UInt uiSize )
{
  Int count = 0;

  for ( Int i = 0; i < uiSize; i++ )
  {
    count += pcCoef[i] != 0;
  }

  return count;
}

#if ALF_HM3_REFACTOR
Void TEncEntropy::codeFiltCountBit(ALFParam* pAlfParam, Int64* ruiRate, const TComSlice * pSlice)
{
  resetEntropy(pSlice);
  resetBits();
  codeFilt(pAlfParam);
  *ruiRate = getNumberOfWrittenBits();
  resetEntropy(pSlice);
  resetBits();
}

Void TEncEntropy::codeAuxCountBit(ALFParam* pAlfParam, Int64* ruiRate, const TComSlice * pSlice)
{
  resetEntropy(pSlice);
  resetBits();
  codeAux(pAlfParam);
  *ruiRate = getNumberOfWrittenBits();
  resetEntropy(pSlice);
  resetBits();
}

#if JVET_C0038_GALF
Void TEncEntropy::codeAux(ALFParam* pAlfParam)
{
  //number of total filters
  {
    Int iNoVarBins = TComAdaptiveLoopFilter::m_NO_VAR_BINS;
    m_pcEntropyCoderIf->xWriteTruncBinCode((UInt)(pAlfParam->noFilters-1), (UInt)iNoVarBins);
  }

  // Filter tap 
  {
  Int FiltTab[3] = {9, 7, 5};
  Int Tab = FiltTab[pAlfParam->realfiltNo];
    m_pcEntropyCoderIf->codeAlfUvlc((Tab-5)/2); 
  }
  Int iNoVarBins = TComAdaptiveLoopFilter::m_NO_VAR_BINS;

  if(pAlfParam->noFilters > 1) //Filter used for each class
  {
    for (Int i=0; i< iNoVarBins; i++)
    {
      m_pcEntropyCoderIf->xWriteTruncBinCode((UInt)pAlfParam->filterPattern[i], (UInt)pAlfParam->noFilters);
    }
  }
  {
    Int i; 
    UChar predPattern, codetab_pred[3] = {1, 0, 2};
    UInt   availableFilters  = (UInt) pAlfParam->iAvailableFilters;

    // frame that filter is predicted from
    if (availableFilters > 0 )
    { 
      // prediction pattern
      predPattern = pAlfParam->PrevFiltIdx[0] > 0 ? 1 : 0;
      for (i=1; i<iNoVarBins; i++)
      {
        Int iCurrPredPattern = pAlfParam->PrevFiltIdx[i] > 0 ? 1: 0;
        if (iCurrPredPattern != predPattern )
        {
          predPattern=2;
          break;
        }
      }
      pAlfParam->iPredPattern = predPattern;
      m_pcEntropyCoderIf->codeALFPrevFiltType(codetab_pred[predPattern]); //0: all zero, no pred from pre-defined filters; 1: all are predicted but could be different values; 2: some predicted and some not

      if (predPattern==2)
      {
        for (i=0; i<iNoVarBins; i++)
        {    
          m_pcEntropyCoderIf->codeALFPrevFiltFlag(pAlfParam->PrevFiltIdx[i]);
        }
      }
      if (predPattern>0 && availableFilters>1)
      {
        for (i=0; i< iNoVarBins; i++)
        {
          if (pAlfParam->PrevFiltIdx[i] > 0)
          {
            m_pcEntropyCoderIf->xWriteTruncBinCode(pAlfParam->PrevFiltIdx[i]-1, availableFilters);
          }
        }
      }
    }
  }
}
#else
Void TEncEntropy::codeAux(ALFParam* pAlfParam)
{
  Int FiltTab[3] = {9, 7, 5};
  Int Tab = FiltTab[pAlfParam->realfiltNo];
  //  m_pcEntropyCoderIf->codeAlfUvlc(pAlfParam->realfiltNo); 

  m_pcEntropyCoderIf->codeAlfUvlc((Tab-5)/2); 

  if (pAlfParam->filtNo>=0)
  {
    if(pAlfParam->realfiltNo >= 0)
    {
      // filters_per_fr
      m_pcEntropyCoderIf->codeAlfUvlc(pAlfParam->noFilters);

      if(pAlfParam->noFilters == 1)
      {
        m_pcEntropyCoderIf->codeAlfUvlc(pAlfParam->startSecondFilter);
      }
      else if (pAlfParam->noFilters == 2)
      {
        for (int i=1; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++) m_pcEntropyCoderIf->codeAlfFlag (pAlfParam->filterPattern[i]);
      }
    }
  }
}
#endif


Int TEncEntropy::lengthGolomb(int coeffVal, int k)
{
  int m = 2 << (k - 1);
  int q = coeffVal / m;
  if(coeffVal != 0)
    return(q + 2 + k);
  else
    return(q + 1 + k);
}
#if JVET_C0038_GALF
Void TEncEntropy::codeFilterCoeffForce0(ALFParam* ALFp)
{
  int filters_per_group = ALFp->filters_per_group;
  int sqrFiltLength = ALFp->num_coeff;
  int filtNo = ALFp->realfiltNo;
  int flTab[]={9/2, 7/2, 5/2};
  int fl = flTab[filtNo];
  int i, k, kMin, kStart, minBits, ind, scanPos, maxScanVal, coeffVal, 
    kMinTab[TComAdaptiveLoopFilter::m_MAX_SQR_FILT_LENGTH], bitsCoeffScan[TComAdaptiveLoopFilter::m_MAX_SCAN_VAL][TComAdaptiveLoopFilter::m_MAX_EXP_GOLOMB],
    minKStart, minBitsKStart, bitsKStart;

  const Int * pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl-2];
  sqrFiltLength=sqrFiltLength-1;

  maxScanVal = 0;
  for(i = 0; i < sqrFiltLength; i++)
{
    maxScanVal = max(maxScanVal, pDepthInt[i]);
  }

  // vlc for all
  memset(bitsCoeffScan, 0, TComAdaptiveLoopFilter::m_MAX_SCAN_VAL * TComAdaptiveLoopFilter::m_MAX_EXP_GOLOMB * sizeof(int));
  for(ind=0; ind<filters_per_group; ++ind)
  {
    if (ALFp->codedVarBins[ind])
    {
      for(i = 0; i < sqrFiltLength; i++)
      {
        scanPos=pDepthInt[i]-1;
        coeffVal=abs(ALFp->coeffmulti[ind][i]);
        for (k=1; k<15; k++)
        {
          bitsCoeffScan[scanPos][k]+=lengthGolomb(coeffVal, k);
        }
      }
    }
  }

  minBitsKStart = 0;
  minKStart = -1;
  for(k = 1; k < 8; k++)
  { 
    bitsKStart = 0; 
    kStart = k;
    for(scanPos = 0; scanPos < maxScanVal; scanPos++)
    {
      kMin = kStart; 
      minBits = bitsCoeffScan[scanPos][kMin];

      if(bitsCoeffScan[scanPos][kStart+1] < minBits)
      {
        kMin = kStart + 1; 
        minBits = bitsCoeffScan[scanPos][kMin];
      }
      kStart = kMin;
      bitsKStart += minBits;
    }
    if((bitsKStart < minBitsKStart) || (k == 1))
    {
      minBitsKStart = bitsKStart;
      minKStart = k;
    }
  }

  kStart = minKStart; 
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    kMin = kStart; 
    minBits = bitsCoeffScan[scanPos][kMin];

    if(bitsCoeffScan[scanPos][kStart+1] < minBits)
    {
      kMin = kStart + 1; 
      minBits = bitsCoeffScan[scanPos][kMin];
    }

    kMinTab[scanPos] = kMin;
    kStart = kMin;
  }

  // Coding parameters
  ALFp->minKStart = minKStart;
  ALFp->maxScanVal = maxScanVal;
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    ALFp->kMinTab[scanPos] = kMinTab[scanPos];
  }
  writeFilterCodingParams(minKStart, maxScanVal, kMinTab, ALFp->forceCoeff0, ALFp->filters_per_group, ALFp->codedVarBins);

  // Filter coefficients
  writeFilterCoeffs(sqrFiltLength, filters_per_group, pDepthInt, ALFp->coeffmulti, kMinTab, ALFp->codedVarBins);  
}
#endif
Int TEncEntropy::codeFilterCoeff(ALFParam* ALFp
#if JVET_C0038_GALF
    , Bool bChroma 
#endif
  )
{
#if JVET_C0038_GALF
  int filters_per_group = bChroma ? 1: ALFp->filters_per_group_diff;
  int sqrFiltLength = bChroma ? ALFp->num_coeff_chroma : ALFp->num_coeff;
  int filtNo = bChroma ? ( ALFp->tap_chroma==5?2: (ALFp->tap_chroma==7? 1: 0) ) : ALFp->realfiltNo;    
#else
  int filters_per_group = ALFp->filters_per_group_diff;
  int sqrFiltLength = ALFp->num_coeff;
  int filtNo = ALFp->realfiltNo;
#endif

  int flTab[]={9/2, 7/2, 5/2};
  int fl = flTab[filtNo];
  int i, k, kMin, kStart, minBits, ind, scanPos, maxScanVal, coeffVal, len = 0,
    kMinTab[TComAdaptiveLoopFilter::m_MAX_SQR_FILT_LENGTH], bitsCoeffScan[TComAdaptiveLoopFilter::m_MAX_SCAN_VAL][TComAdaptiveLoopFilter::m_MAX_EXP_GOLOMB],
    minKStart, minBitsKStart, bitsKStart;

  const Int * pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl-2];
#if JVET_C0038_GALF
  sqrFiltLength=sqrFiltLength-1;
#endif

  maxScanVal = 0;
  for(i = 0; i < sqrFiltLength; i++)
    maxScanVal = max(maxScanVal, pDepthInt[i]);

  // vlc for all
  memset(bitsCoeffScan, 0, TComAdaptiveLoopFilter::m_MAX_SCAN_VAL * TComAdaptiveLoopFilter::m_MAX_EXP_GOLOMB * sizeof(int));
  for(ind=0; ind<filters_per_group; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos=pDepthInt[i]-1;
#if JVET_C0038_GALF
      coeffVal=abs(bChroma? ALFp->coeff_chroma[i]: ALFp->coeffmulti[ind][i]);
#else
      coeffVal=abs(ALFp->coeffmulti[ind][i]);
#endif
      for (k=1; k<15; k++)
      {
        bitsCoeffScan[scanPos][k]+=lengthGolomb(coeffVal, k);
      }
    }
  }

  minBitsKStart = 0;
  minKStart = -1;
  for(k = 1; k < 8; k++)
  { 
    bitsKStart = 0; 
    kStart = k;
    for(scanPos = 0; scanPos < maxScanVal; scanPos++)
    {
      kMin = kStart; 
      minBits = bitsCoeffScan[scanPos][kMin];

      if(bitsCoeffScan[scanPos][kStart+1] < minBits)
      {
        kMin = kStart + 1; 
        minBits = bitsCoeffScan[scanPos][kMin];
      }
      kStart = kMin;
      bitsKStart += minBits;
    }
    if((bitsKStart < minBitsKStart) || (k == 1))
    {
      minBitsKStart = bitsKStart;
      minKStart = k;
    }
  }

  kStart = minKStart; 
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    kMin = kStart; 
    minBits = bitsCoeffScan[scanPos][kMin];

    if(bitsCoeffScan[scanPos][kStart+1] < minBits)
    {
      kMin = kStart + 1; 
      minBits = bitsCoeffScan[scanPos][kMin];
    }

    kMinTab[scanPos] = kMin;
    kStart = kMin;
  }

  // Coding parameters
  ALFp->minKStart = minKStart;
  ALFp->maxScanVal = maxScanVal;
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    ALFp->kMinTab[scanPos] = kMinTab[scanPos];
  }
#if JVET_C0038_GALF
  if(bChroma)
  {
    len += writeFilterCodingParams(minKStart, maxScanVal, kMinTab, 0, ALFp->filters_per_group, ALFp->codedVarBins);
  }
  else
  {    
    len += writeFilterCodingParams(minKStart, maxScanVal, kMinTab, ALFp->forceCoeff0, ALFp->filters_per_group, ALFp->codedVarBins);
  }
  // Filter coefficients
  if(bChroma)
  {
    len += writeFilterCoeffsForChroma(sqrFiltLength, pDepthInt, ALFp->coeff_chroma, kMinTab);
  }
  else
  {
    len += writeFilterCoeffs(sqrFiltLength, filters_per_group, pDepthInt, ALFp->coeffmulti, kMinTab, ALFp->codedVarBins);
  }
#else
  len += writeFilterCodingParams(minKStart, maxScanVal, kMinTab);

  // Filter coefficients
  len += writeFilterCoeffs(sqrFiltLength, filters_per_group, pDepthInt, ALFp->coeffmulti, kMinTab);
#endif
  return len;
}

#if JVET_C0038_GALF
Int TEncEntropy::writeFilterCodingParams(Int minKStart, Int maxScanVal, Int kMinTab[], Bool forceCoeff0, Int filters_per_group, Bool codedVarBins[])
#else
Int TEncEntropy::writeFilterCodingParams(int minKStart, int maxScanVal, int kMinTab[])
#endif
{
  int scanPos;
  int golombIndexBit;
  int kMin;

  // Golomb parameters
  m_pcEntropyCoderIf->codeAlfUvlc(minKStart - 1);

  kMin = minKStart; 
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    golombIndexBit = (kMinTab[scanPos] != kMin)? 1: 0;

    assert(kMinTab[scanPos] <= kMin + 1);

    m_pcEntropyCoderIf->codeAlfFlag(golombIndexBit);
    kMin = kMinTab[scanPos];
  }    
#if JVET_C0038_GALF  
  if(forceCoeff0)
  {
    for(Int ind=0; ind<filters_per_group; ++ind)
    {
      m_pcEntropyCoderIf->codeALFPrevFiltFlag( codedVarBins[ind] ? 1 : 0 );
    }
  }
#endif
  return 0;
}

Int TEncEntropy::writeFilterCoeffs(int sqrFiltLength, int filters_per_group, const int pDepthInt[], 
  int **FilterCoeff, int kMinTab[]
#if JVET_C0038_GALF
 ,Bool codedVarBins[]
#endif
)
{
  int ind, scanPos, i;

  for(ind = 0; ind < filters_per_group; ++ind)
  {
#if JVET_C0038_GALF
    if(codedVarBins[ind]==false)
    {
      continue;
    }
#endif
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos = pDepthInt[i] - 1;
      golombEncode(FilterCoeff[ind][i], kMinTab[scanPos]);
    }
  }
  return 0;
}
#if JVET_C0038_GALF
Int TEncEntropy::writeFilterCoeffsForChroma(Int sqrFiltLength, const Int pDepthInt[], Int *FilterCoeff, Int kMinTab[])
{
  Int scanPos, i;  
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos = pDepthInt[i] - 1;
      golombEncode(FilterCoeff[i], kMinTab[scanPos]);
    }
  }
  return 0;
}
#endif
Int TEncEntropy::golombEncode(int coeff, int k)
{
  int q, i, m;
  int symbol = abs(coeff);

  m = (int)pow(2.0, k);
  q = symbol / m;

  for (i = 0; i < q; i++)
    m_pcEntropyCoderIf->codeAlfFlag(1);
  m_pcEntropyCoderIf->codeAlfFlag(0);
  // write one zero

  for(i = 0; i < k; i++)
  {
    m_pcEntropyCoderIf->codeAlfFlag(symbol & 0x01);
    symbol >>= 1;
  }

  if(coeff != 0)
  {
    int sign = (coeff > 0)? 1: 0;
    m_pcEntropyCoderIf->codeAlfFlag(sign);
  }
  return 0;
}

Void TEncEntropy::codeFilt(ALFParam* pAlfParam)
{
#if JVET_C0038_GALF
  m_pcEntropyCoderIf->codeALFPrevFiltFlag(pAlfParam->forceCoeff0);
  if( pAlfParam->forceCoeff0 )
  {
    codeFilterCoeffForce0(pAlfParam);
  }
  else
  {
#endif
  if(pAlfParam->filters_per_group > 1)
  {
    m_pcEntropyCoderIf->codeAlfFlag (pAlfParam->predMethod);
  }
  codeFilterCoeff (pAlfParam);
#if JVET_C0038_GALF
}
#endif
}
#if !JVET_C0038_GALF
Void  print(ALFParam* pAlfParam)
{
  Int i=0;
  Int ind=0;
  Int FiltLengthTab[] = {22, 14, 8}; //0:9tap
  Int FiltLength = FiltLengthTab[pAlfParam->realfiltNo];

  printf("set of params\n");
  printf("realfiltNo:%d\n", pAlfParam->realfiltNo);
  printf("filtNo:%d\n", pAlfParam->filtNo);
  printf("filterPattern:");
  for (i=0; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++) printf("%d ", pAlfParam->filterPattern[i]);
  printf("\n");

  printf("startSecondFilter:%d\n", pAlfParam->startSecondFilter);
  printf("noFilters:%d\n", pAlfParam->noFilters);
  printf("varIndTab:");
  for (i=0; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++) printf("%d ", pAlfParam->varIndTab[i]);
  printf("\n");
  printf("filters_per_group_diff:%d\n", pAlfParam->filters_per_group_diff);
  printf("filters_per_group:%d\n", pAlfParam->filters_per_group);
  printf("codedVarBins:");
  for (i=0; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++) printf("%d ", pAlfParam->codedVarBins[i]);
  printf("\n");
  printf("forceCoeff0:%d\n", pAlfParam->forceCoeff0);
  printf("predMethod:%d\n", pAlfParam->predMethod);

  for (ind=0; ind<pAlfParam->filters_per_group_diff; ind++)
  {
    printf("coeffmulti(%d):", ind);
    for (i=0; i<FiltLength; i++) printf("%d ", pAlfParam->coeffmulti[ind][i]);
    printf("\n");
  }

  printf("minKStart:%d\n", pAlfParam->minKStart);  
  printf("maxScanVal:%d\n", pAlfParam->maxScanVal);  
  printf("kMinTab:");
  for(Int scanPos = 0; scanPos < pAlfParam->maxScanVal; scanPos++)
  {
    printf("%d ", pAlfParam->kMinTab[scanPos]);
  }
  printf("\n");

  printf("chroma_idc:%d\n", pAlfParam->chroma_idc);  
  printf("tap_chroma:%d\n", pAlfParam->tap_chroma);  
  printf("chroma_coeff:");
  for(Int scanPos = 0; scanPos < pAlfParam->num_coeff_chroma; scanPos++)
  {
    printf("%d ", pAlfParam->coeff_chroma[scanPos]);
  }
  printf("\n");
}
#endif

Void TEncEntropy::encodeAlfParam(ALFParam* pAlfParam , UInt uiMaxTotalCUDepth
#if FIX_TICKET12
        ,const TComSlice * pSlice
#endif
  )
{
  m_pcEntropyCoderIf->codeAlfFlag(pAlfParam->alf_flag);
  if (!pAlfParam->alf_flag)
    return;
#if !JVET_C0038_GALF 
  Int pos;
#endif

#if COM16_C806_ALF_TEMPPRED_NUM
  //encode temporal prediction flag and index
#if FIX_TICKET12
  if( pSlice->getSliceType() == I_SLICE)
  {
    assert(pAlfParam->temproalPredFlag == false);
  }
  else
#endif
  m_pcEntropyCoderIf->codeAlfFlag( pAlfParam->temproalPredFlag ? 1 : 0 );
  if( pAlfParam->temproalPredFlag )
  {
    m_pcEntropyCoderIf->codeAlfUvlc( pAlfParam->prevIdx );
  }
  else
  {
#endif
    codeAux(pAlfParam);
    codeFilt(pAlfParam);
#if COM16_C806_ALF_TEMPPRED_NUM
  }
#endif
  // filter parameters for chroma
  m_pcEntropyCoderIf->codeAlfUvlc(pAlfParam->chroma_idc);

#if COM16_C806_ALF_TEMPPRED_NUM
  if( !pAlfParam->temproalPredFlag && pAlfParam->chroma_idc )
#else
  if(pAlfParam->chroma_idc)
#endif
  {
#if !JVET_C0038_GALF 
    m_pcEntropyCoderIf->codeAlfUvlc((pAlfParam->tap_chroma-5)/2);
#endif
    // filter coefficients for chroma
#if JVET_C0038_GALF 
    codeFilterCoeff(pAlfParam, true);
#else
    for(pos=0; pos<pAlfParam->num_coeff_chroma; pos++)
    {
      m_pcEntropyCoderIf->codeAlfSvlc(pAlfParam->coeff_chroma[pos]);
    }
#endif
  }
  // region control parameters for luma
  m_pcEntropyCoderIf->codeAlfFlag(pAlfParam->cu_control_flag);
  if (pAlfParam->cu_control_flag)
  {
    assert( (pAlfParam->cu_control_flag && m_pcEntropyCoderIf->getAlfCtrl()) || (!pAlfParam->cu_control_flag && !m_pcEntropyCoderIf->getAlfCtrl()));
    m_pcEntropyCoderIf->codeAlfCtrlDepth( uiMaxTotalCUDepth );
  }
}

#if !JVET_C0024_QTBT
Void TEncEntropy::encodeAlfCtrlFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD )
{
  if( bRD )
    uiAbsPartIdx = 0;

  m_pcEntropyCoderIf->codeAlfCtrlFlag( pcCU, uiAbsPartIdx );
}
#endif

Void TEncEntropy::encodeAlfCtrlParam( ALFParam* pAlfParam )
{
  m_pcEntropyCoderIf->codeAlfFlagNum( pAlfParam->num_alf_cu_flag, pAlfParam->num_cus_in_frame );

  for(UInt i=0; i<pAlfParam->num_alf_cu_flag; i++)
  {
    m_pcEntropyCoderIf->codeAlfCtrlFlag( pAlfParam->alf_cu_flag[i] );
  }
}
#endif

#if COM16_C1016_AFFINE
Void TEncEntropy::encodeAffineFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPuIdx )
{
  if ( pcCU->getSlice()->isIntra() )
  {
    return;
  }

  if ( !pcCU->getSlice()->getSPS()->getUseAffine() )
  {
    return;
  }

  m_pcEntropyCoderIf->codeAffineFlag( pcCU, uiAbsPartIdx );
}
#endif

#if MRG_CONSIDER_COMPLEX
Void TEncEntropy::encodeMergeComplexFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  if (pcCU->getSlice()->isIntra())
  {
    return;
  }

  if (!pcCU->getSlice()->getSPS()->getUseMrgComplex())
  {
    return;
  }

  assert((pcCU->getWidth(uiAbsPartIdx) >= 16) || (pcCU->getHeight(uiAbsPartIdx) >= 16));
  m_pcEntropyCoderIf->codeMergeComplexFlag(pcCU, uiAbsPartIdx);
}

Void TEncEntropy::encodeMergeComplexIndex(TComDataCU* pcCU, UInt uiAbsPartIdx, Bool bRD)
{
  if (bRD)
  {
    uiAbsPartIdx = 0;
  }
  m_pcEntropyCoderIf->codeMergeComplexIndex(pcCU, uiAbsPartIdx);
}
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
Void TEncEntropy::encodeAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPuIdx)
{
  if (pcCU->getSlice()->isIntra())
  {
    return;
  }

  if (!pcCU->getSlice()->getSPS()->getUseAffine())
  {
    return;
  }

  m_pcEntropyCoderIf->codeAffineParamFlag(pcCU, uiAbsPartIdx);
}

#endif

#if RRC_INTRA_SHARPENING
Void TEncEntropy::encodeIntraSharpFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  if (pcCU->isIntraSharpFlagCoded(uiAbsPartIdx))
    m_pcEntropyCoderIf->codeIntraSharpFlag(pcCU->getIntraSharpFlag(uiAbsPartIdx));
}
#endif

#if RRC_SELECTABLE_IF
Void TEncEntropy::encodeAdaptiveIfIndex(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
#if SIF_CU_LEVEL
  if (pcCU->isSIFFlagCoded(uiAbsPartIdx))
#endif
  {
    m_pcEntropyCoderIf->codeAdaptiveIfIndex(pcCU->getSIFFlag(uiAbsPartIdx));
  }
}
#endif
//! \}
