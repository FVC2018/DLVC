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

/** \file     TDecEntropy.cpp
    \brief    entropy decoder class
*/

#include "TDecEntropy.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/TComPrediction.h"
#if ALF_HM3_REFACTOR
#include "TLibCommon/TComAdaptiveLoopFilter.h"
#include <cmath>
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
#include "../TLibCommon/Debug.h"
static const Bool bDebugRQT = DebugOptionList::DebugRQT.getInt()!=0;
static const Bool bDebugPredEnabled = DebugOptionList::DebugPred.getInt()!=0;
#endif

//! \ingroup TLibDecoder
//! \{

Void TDecEntropy::setEntropyDecoder         ( TDecEntropyIf* p )
{
  m_pcEntropyDecoderIf = p;
}

#include "TLibCommon/TComSampleAdaptiveOffset.h"

Void TDecEntropy::decodeSkipFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parseSkipFlag( pcCU, uiAbsPartIdx, uiDepth );
}

#if BARC_UP_SAMPLE_RDO
Void TDecEntropy::decodeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
	m_pcEntropyDecoderIf->parseUpCNNFlag(pcCU, uiAbsPartIdx);
}
#endif

#if BARC
Void TDecEntropy::decodeBARCFlag( TComDataCU* pcCU, ChannelType eCtype )
{
  m_pcEntropyDecoderIf->parseBARCFlag( pcCU,eCtype );
}
#endif
#if VCEG_AZ07_IMV
Void TDecEntropy::decodeiMVFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  assert( pcCU->getSlice()->getSPS()->getIMV() );
  m_pcEntropyDecoderIf->parseiMVFlag( pcCU, uiAbsPartIdx, uiDepth );
}
#endif

#if COM16_C806_OBMC
Void TDecEntropy::decodeOBMCFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parseOBMCFlag( pcCU, uiAbsPartIdx, uiDepth );
}
#endif

#if VCEG_AZ06_IC
Void TDecEntropy::decodeICFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  if( pcCU->isICFlagCoded( uiAbsPartIdx ) )
  {
    m_pcEntropyDecoderIf->parseICFlag( pcCU, uiAbsPartIdx, uiDepth );
  }
}
#endif

Void TDecEntropy::decodeCUTransquantBypassFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parseCUTransquantBypassFlag( pcCU, uiAbsPartIdx, uiDepth );
}

#if VCEG_AZ05_INTRA_MPI
Void TDecEntropy::decodeMPIIdx(TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  m_pcEntropyDecoderIf->parseMPIIdx(pcSubCU, uiAbsPartIdx, uiDepth);
}
#endif

#if DIMD_INTRA_PRED
Void TDecEntropy::decodeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
  m_pcEntropyDecoderIf->parseDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
}
Void TDecEntropy::decodeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
  m_pcEntropyDecoderIf->parseDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
}
#endif

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
Void TDecEntropy::decodePDPCIdx(TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  m_pcEntropyDecoderIf->parsePDPCIdx(pcSubCU, uiAbsPartIdx, uiDepth);
}
#endif


#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
Void TDecEntropy::decodeROTIdx( TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth )
{ 
  m_pcEntropyDecoderIf->parseROTIdx( pcSubCU, uiAbsPartIdx, uiDepth );
#if JVET_C0024_QTBT
  if( !pcSubCU->getSlice()->isIntra() )
  {
    pcSubCU->setROTIdxSubParts( CHANNEL_TYPE_CHROMA, pcSubCU->getROTIdx( CHANNEL_TYPE_LUMA, uiAbsPartIdx ), uiAbsPartIdx,  uiDepth ); 
}
#endif
}

#if JVET_C0024_QTBT
Void TDecEntropy::decodeROTIdxChroma( TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth )
{ 
  m_pcEntropyDecoderIf->parseROTIdxChroma( pcSubCU, uiAbsPartIdx, uiDepth );
}
#endif
#endif
/** decode merge flag
 * \param pcSubCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \param uiPUIdx
 * \returns Void
 */
Void TDecEntropy::decodeMergeFlag( TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx )
{
  // at least one merge candidate exists
  m_pcEntropyDecoderIf->parseMergeFlag( pcSubCU, uiAbsPartIdx, uiDepth, uiPUIdx );
}

/** decode merge index
 * \param pcCU
 * \param uiPartIdx
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TDecEntropy::decodeMergeIndex( TComDataCU* pcCU, UInt uiPartIdx, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiMergeIndex = 0;
  m_pcEntropyDecoderIf->parseMergeIndex( pcCU, uiMergeIndex );
  pcCU->setMergeIndexSubParts( uiMergeIndex, uiAbsPartIdx, uiPartIdx, uiDepth );
}

#if VCEG_AZ07_FRUC_MERGE
Void TDecEntropy::decodeFRUCMgrMode( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx )
{ 
  // at least one merge candidate exists
  m_pcEntropyDecoderIf->parseFRUCMgrMode( pcCU, uiAbsPartIdx, uiDepth, uiPUIdx );
}
#endif

#if JVET_C0024_QTBT
#if TT_TREE_DS1
Void TDecEntropy::decodeBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, UInt uiBTDepth )
{
  m_pcEntropyDecoderIf->parseBTSplitMode( pcCU, uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
}
#else
Void TDecEntropy::decodeBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight )
{
  m_pcEntropyDecoderIf->parseBTSplitMode( pcCU, uiAbsPartIdx, uiWidth, uiHeight );
}
#endif
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
Void TDecEntropy::decodeCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
	m_pcEntropyDecoderIf->parseCUBgFlag(pcCU, uiAbsPartIdx, uiDepth);
}
#endif

Void TDecEntropy::decodeSplitFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parseSplitFlag( pcCU, uiAbsPartIdx, uiDepth );
}

Void TDecEntropy::decodePredMode( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parsePredMode( pcCU, uiAbsPartIdx, uiDepth );
}

#if !JVET_C0024_QTBT
Void TDecEntropy::decodePartSize( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  m_pcEntropyDecoderIf->parsePartSize( pcCU, uiAbsPartIdx, uiDepth );
}
#endif

Void TDecEntropy::decodePredInfo    ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, TComDataCU* pcSubCU )
{
  if( pcCU->isIntra( uiAbsPartIdx ) )                                 // If it is Intra mode, encode intra prediction mode.
  {
#if JVET_C0024_QTBT
    if (isLuma(pcCU->getTextType()))
    {
#endif
// #if PRINTF_INFO
//         printf("CurrLuma CU x is %d, y is %d, absidx is %d\n", pcCU->getCUPelX(), pcCU->getCUPelX(), uiAbsPartIdx);
//         if ((pcCU->getCUPelX()==128)&&(pcCU->getCUPelX()==128)&&(uiAbsPartIdx==991))
//         {
//             int LumaStopHere = 1;
//         }
// #endif
    decodeIntraDirModeLuma  ( pcCU, uiAbsPartIdx, uiDepth );
#if JVET_C0024_QTBT
    }
    if (pcCU->getPic()->getChromaFormat()!=CHROMA_400
      && (isChroma(pcCU->getTextType()) || !pcCU->getSlice()->isIntra()))
#else
    if (pcCU->getPic()->getChromaFormat()!=CHROMA_400)
#endif
    {
// #if PRINTF_INFO
//         printf("CurrChroma CU x is %d, y is %d, absidx is %d\n", pcCU->getCUPelX(), pcCU->getCUPelX(), uiAbsPartIdx);
//         if ((pcCU->getCUPelX() == 128) && (pcCU->getCUPelX() == 128) && (uiAbsPartIdx == 714))
//         {
//             int ChromaStopHere = 1;
//   }
// #endif
      decodeIntraDirModeChroma( pcCU, uiAbsPartIdx, uiDepth );
#if !JVET_C0024_QTBT
      if (enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) && pcCU->getPartitionSize( uiAbsPartIdx )==SIZE_NxN)
      {
        UInt uiPartOffset = ( pcCU->getPic()->getNumPartitionsInCtu() >> ( pcCU->getDepth(uiAbsPartIdx) << 1 ) ) >> 2;
        decodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset,   uiDepth+1 );
        decodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset*2, uiDepth+1 );
        decodeIntraDirModeChroma( pcCU, uiAbsPartIdx + uiPartOffset*3, uiDepth+1 );
      }
#endif
    }
#if RRC_INTRA_SHARPENING
#if JVET_C0024_QTBT
    if (isLuma(pcCU->getTextType()))
#endif
    {
      decodeIntraSharpFlag(pcCU, uiAbsPartIdx, uiDepth);
    }
#endif
  }
  else                                                                // if it is Inter mode, encode motion vector and reference index
  {
    decodePUWise( pcCU, uiAbsPartIdx, uiDepth, pcSubCU );
  }
}

/** Parse I_PCM information.
 * \param pcCU  pointer to CUpointer to CU
 * \param uiAbsPartIdx CU index
 * \param uiDepth CU depth
 * \returns Void
 */
Void TDecEntropy::decodeIPCMInfo( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  if(!pcCU->getSlice()->getSPS()->getUsePCM()
    || pcCU->getWidth(uiAbsPartIdx) > (1<<pcCU->getSlice()->getSPS()->getPCMLog2MaxSize())
    || pcCU->getWidth(uiAbsPartIdx) < (1<<pcCU->getSlice()->getSPS()->getPCMLog2MinSize()) )
  {
    return;
  }

  m_pcEntropyDecoderIf->parseIPCMInfo( pcCU, uiAbsPartIdx, uiDepth );
}

Void TDecEntropy::decodeIntraDirModeLuma  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
#if DIMD_INTRA_PRED
  if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
  {
    return;
  }
#endif
  m_pcEntropyDecoderIf->parseIntraDirLumaAng( pcCU, uiAbsPartIdx, uiDepth );
}

Void TDecEntropy::decodeIntraDirModeChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
#if DIMD_INTRA_PRED
  if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx))
  {
    return;
  }
#endif
  m_pcEntropyDecoderIf->parseIntraDirChroma( pcCU, uiAbsPartIdx, uiDepth );

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  if (bDebugPredEnabled)
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


/** decode motion information for every PU block.
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \param pcSubCU
 * \returns Void
 */
Void TDecEntropy::decodePUWise( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, TComDataCU* pcSubCU )
{
#if JVET_C0024_QTBT
  PartSize ePartSize = SIZE_2Nx2N;
  UInt uiNumPU = 1;
  UInt uiPUOffset = 0;
#else
  PartSize ePartSize = pcCU->getPartitionSize( uiAbsPartIdx );
  UInt uiNumPU = ( ePartSize == SIZE_2Nx2N ? 1 : ( ePartSize == SIZE_NxN ? 4 : 2 ) );
  UInt uiPUOffset = ( g_auiPUOffset[UInt( ePartSize )] << ( ( pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth ) << 1 ) ) >> 4;
#endif
#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT 
  TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UChar    eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
  memset( eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, sizeof(UChar)*MRG_MAX_NUM_CANDS );
#endif
#if VCEG_AZ06_IC
  Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif

  for ( UInt ui = 0; ui < pcCU->getSlice()->getMaxNumMergeCand(); ui++ )
  {
    uhInterDirNeighbours[ui] = 0;
  }
  Int numValidMergeCand = 0;
  Bool hasMergedCandList = false;
#endif

#if VCEG_AZ07_FRUC_MERGE && JVECT_C0062_AFFINE_SIX_PARAM && JVET_C0024_QTBT
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0;
  for (UInt ui = 0; ui < pcCU->getSlice()->getMaxNumMergeCand(); ui++)
  {
    uhInterDirNeighbours[ui] = 0;
  }
#endif

#if VCEG_AZ07_IMV
  Bool bNonZeroMvd = false;
#endif

  pcSubCU->copyInterPredInfoFrom( pcCU, uiAbsPartIdx, REF_PIC_LIST_0 );
  pcSubCU->copyInterPredInfoFrom( pcCU, uiAbsPartIdx, REF_PIC_LIST_1 );
  for ( UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset )
  {
    decodeMergeFlag( pcCU, uiSubPartIdx, uiDepth, uiPartIdx );
    if ( pcCU->getMergeFlag( uiSubPartIdx ) )
    {
#if VCEG_AZ07_FRUC_MERGE
      decodeFRUCMgrMode( pcCU , uiSubPartIdx , uiDepth , uiPartIdx );
      if( !pcCU->getFRUCMgrMode( uiSubPartIdx ) )
      {
#endif
#if COM16_C1016_AFFINE
#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
        TComMvField  cAffineMvField[2][3];
#endif
#if MRG_CONSIDER_COMPLEX
        if ((pcCU->getWidth(uiSubPartIdx) >= 16) || (pcCU->getHeight(uiSubPartIdx) >= 16))
        {
          decodeMergeComplexFlag(pcCU, uiSubPartIdx, uiDepth, uiPartIdx);
          if (pcCU->getMergeComplexFlag(uiSubPartIdx))
          {
            decodeMergeComplexIndex(pcCU, uiPartIdx, uiSubPartIdx, uiDepth);
          }
        }

        if (!pcCU->getMergeComplexFlag(uiSubPartIdx))
        {
#endif
        if ( ePartSize == SIZE_2Nx2N && pcCU->isAffineMrgFlagCoded(uiSubPartIdx, uiPartIdx ) )
        {
          decodeAffineFlag( pcCU, uiSubPartIdx, uiDepth, uiPartIdx );
        }

        if ( pcCU->isAffine(uiSubPartIdx) )
        {
#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT
#if JVECT_C0062_AFFINE_SIX_PARAM
          UInt neighborAffineParam = 0;
          pcSubCU->getAffineMergeCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
          pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiSubPartIdx, uiPartIdx, uiDepth);
#else
          pcSubCU->getAffineMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand );
#endif
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM && VCEG_AZ07_FRUC_MERGE
          UInt neighborAffineParam = 0;
          pcSubCU->getAffineMergeCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
          pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiSubPartIdx, uiPartIdx, uiDepth);
#endif

          pcCU->setMergeIndexSubParts( 0, uiSubPartIdx, uiPartIdx, uiDepth );
        }
        else
        {
          decodeMergeIndex( pcCU, uiPartIdx, uiSubPartIdx, uiDepth );
        }
#else
      decodeMergeIndex( pcCU, uiPartIdx, uiSubPartIdx, uiDepth );
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
      if (bDebugPredEnabled)
      {
        std::cout << "Coded merge flag, CU absPartIdx: " << uiAbsPartIdx << " PU(" << uiPartIdx << ") absPartIdx: " << uiSubPartIdx;
        std::cout << " merge index: " << (UInt)pcCU->getMergeIndex(uiSubPartIdx) << std::endl;
      }
#endif

#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT
      UInt uiMergeIndex = pcCU->getMergeIndex(uiSubPartIdx);
#if COM16_C1016_AFFINE
      if ( pcCU->isAffine(uiSubPartIdx) )
      {
#if JVECT_C0062_AFFINE_SIX_PARAM
        UInt neighborAffineParam = 0;
        pcSubCU->getAffineMergeCandidates(uiSubPartIdx - uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
        pcCU->setAffineParamFlagSubParts(neighborAffineParam, uiSubPartIdx, uiPartIdx, uiDepth);
#else
        pcSubCU->getAffineMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx, cAffineMvField, uhInterDirNeighbours, numValidMergeCand );
#endif
        pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth );
        TComMv cTmpMv( 0, 0 );
#if JVECT_C0062_AFFINE_SIX_PARAM
        Bool b6Param = neighborAffineParam;//ZF
#endif
        for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
        {
          if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
          {
            TComMvField* pcMvField = cAffineMvField[ 2 * uiMergeIndex + uiRefListIdx ];

            pcCU->setMVPIdxSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx, uiPartIdx, uiDepth);
            pcCU->setMVPNumSubParts( 0, RefPicList( uiRefListIdx ), uiSubPartIdx, uiPartIdx, uiDepth);
            pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cTmpMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );
#if JVECT_C0062_AFFINE_SIX_PARAM
            pcCU->setAllAffineMvField(uiSubPartIdx, uiPartIdx, pcMvField, RefPicList(uiRefListIdx), uiDepth, b6Param);
#else
            pcCU->setAllAffineMvField( uiSubPartIdx, uiPartIdx, pcMvField, RefPicList(uiRefListIdx), uiDepth );
#endif
          }
        }
      }
      else
      {
#endif
      if ( pcCU->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() && ePartSize != SIZE_2Nx2N && pcSubCU->getWidth( 0 ) <= 8 )
      {
        if ( !hasMergedCandList )
        {
#if !JVET_C0024_QTBT
          pcSubCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth ); // temporarily set.
#endif
          pcSubCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand 
#if VCEG_AZ06_IC
          , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          , eMergeCandTypeNieghors
          , m_pMvFieldSP
          , m_phInterDirSP
          , uiSubPartIdx
          , pcCU
#endif
            );
#if !JVET_C0024_QTBT
          pcSubCU->setPartSizeSubParts( ePartSize, 0, uiDepth ); // restore.
#endif
          hasMergedCandList = true;
        }
      }
      else
      {
        uiMergeIndex = pcCU->getMergeIndex(uiSubPartIdx);
        pcSubCU->getInterMergeCandidates( uiSubPartIdx-uiAbsPartIdx, uiPartIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#if VCEG_AZ06_IC
          , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          , eMergeCandTypeNieghors
          , m_pMvFieldSP
          , m_phInterDirSP
          , uiSubPartIdx
          , pcCU
#endif
          , uiMergeIndex
          );
      }

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
      pcCU->setMergeTypeSubParts( eMergeCandTypeNieghors[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth ); 
      if( eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_DEFAULT_N )
#endif
      pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth );
#if VCEG_AZ06_IC
#if !JVET_C0024_QTBT
      if( pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
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
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          if( eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_DEFAULT_N )
#endif
          pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvField( cMvFieldNeighbours[ 2*uiMergeIndex + uiRefListIdx ], ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );

        }
      }
#if COM16_C1016_AFFINE
      }
#endif
#endif
#if MRG_CONSIDER_COMPLEX
      }
#endif
#if VCEG_AZ07_FRUC_MERGE
      }
#endif
    }
    else
    {
      decodeInterDirPU( pcCU, uiSubPartIdx, uiDepth, uiPartIdx );
#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
      if ( pcCU->getWidth(uiSubPartIdx) > 8 && pcCU->getHeight(uiSubPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#else
      if ( pcCU->getWidth(uiSubPartIdx) > 8 && ePartSize == SIZE_2Nx2N )
#endif
      {
        decodeAffineFlag ( pcCU, uiSubPartIdx, uiDepth, uiPartIdx );
#if JVECT_C0062_AFFINE_SIX_PARAM
        if (pcCU->isAffine(uiSubPartIdx))
          decodeAffineParamFlag(pcCU, uiSubPartIdx, uiDepth, uiPartIdx);
#endif
      }
#endif

      for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
      {
        if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
        {
          decodeRefFrmIdxPU( pcCU,    uiSubPartIdx,              uiDepth, uiPartIdx, RefPicList( uiRefListIdx ) );
          decodeMvdPU      ( pcCU,    uiSubPartIdx,              uiDepth, uiPartIdx, RefPicList( uiRefListIdx ) );
#if VCEG_AZ07_IMV
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiSubPartIdx ).getHor() != 0 );
          bNonZeroMvd |= ( pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMvd( uiSubPartIdx ).getVer() != 0 );
#endif
#if TT_TREE_DS3 //important fix when TT is ON
          decodeMVPIdxPU   ( pcCU, uiSubPartIdx, uiDepth, uiPartIdx, RefPicList( uiRefListIdx ) );
#else
          decodeMVPIdxPU   ( pcSubCU, uiSubPartIdx-uiAbsPartIdx, uiDepth, uiPartIdx, RefPicList( uiRefListIdx ) );
#endif
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
          if (bDebugPredEnabled)
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

    if ( (pcCU->getInterDir(uiSubPartIdx) == 3) && pcSubCU->isBipredRestriction(uiPartIdx) )
    {
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMv( TComMv(0,0), ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllRefIdx( -1, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->setInterDirSubParts( 1, uiSubPartIdx, uiPartIdx, uiDepth);
    }
  }

#if VCEG_AZ07_IMV
  if( bNonZeroMvd && pcCU->getSlice()->getSPS()->getIMV() )
  {
    decodeiMVFlag( pcCU , uiAbsPartIdx , uiDepth );
  }

#if !VCEG_AZ07_FRUC_MERGE && !JVET_C0024_QTBT
  for ( UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset )
  {
    if( !pcCU->getMergeFlag( uiSubPartIdx ) && pcCU->getiMVFlag( uiSubPartIdx ) )
    {
      for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
      {
        RefPicList eRefList = RefPicList( uiRefListIdx );
        if ( pcCU->getSlice()->getNumRefIdx( eRefList ) > 0 && ( pcCU->getInterDir( uiSubPartIdx ) & ( 1 << uiRefListIdx ) ) )
        {
          Short iMvdHor = pcCU->getCUMvField( eRefList )->getMvd( uiSubPartIdx ).getHor();
          Short iMvdVer = pcCU->getCUMvField( eRefList )->getMvd( uiSubPartIdx ).getVer();
          iMvdHor <<= 2;
          iMvdVer <<= 2;
#if  JVET_E0076_MULTI_PEL_MVD
          if( pcCU->getiMVFlag( uiSubPartIdx ) == 2)
          {
            iMvdHor <<= MULTI_PEL_MVD_BITS;
            iMvdVer <<= MULTI_PEL_MVD_BITS;
          }
#endif
          TComMv cMv( iMvdHor, iMvdVer );

          pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cMv, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx );

          AMVPInfo* pAMVPInfo = pcSubCU->getCUMvField( eRefList )->getAMVPInfo();
          pcSubCU->fillMvpCand(uiPartIdx, uiSubPartIdx - uiAbsPartIdx, eRefList, pcSubCU->getCUMvField( eRefList )->getRefIdx( uiSubPartIdx - uiAbsPartIdx ), pAMVPInfo);
          m_pcPrediction->getMvPredAMVP( pcSubCU, uiPartIdx, uiSubPartIdx - uiAbsPartIdx, RefPicList( uiRefListIdx ), cMv );
          cMv += TComMv( iMvdHor, iMvdVer );
          pcSubCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMv( cMv, ePartSize, uiSubPartIdx - uiAbsPartIdx, 0, uiPartIdx);
        }
      }
    }
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL && VCEG_AZ07_IMV
    else if( pcCU->getMergeFlag( uiSubPartIdx ) && pcCU->getiMVFlag( uiSubPartIdx ) )
    {
      UInt uiMergeIndex = pcCU->getMergeIndex(uiSubPartIdx);
      if ( pcCU->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() && ePartSize != SIZE_2Nx2N && pcSubCU->getWidth( 0 ) <= 8 ) 
      {
#if !JVET_C0024_QTBT
        pcSubCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth );
#endif
        if ( !hasMergedCandList )
        {
          //          pcSubCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand );
          pcSubCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand 
#if VCEG_AZ06_IC
            , abICFlag
#endif
            , eMergeCandTypeNieghors
            , m_pMvFieldSP
            , m_phInterDirSP
            );
          hasMergedCandList = true;
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
          , eMergeCandTypeNieghors
          , m_pMvFieldSP
          , m_phInterDirSP
          );
      }
      pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiSubPartIdx, uiPartIdx, uiDepth );
#if VCEG_AZ06_IC
#if !JVET_C0024_QTBT
      if( pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
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
    }

    if ( (pcCU->getInterDir(uiSubPartIdx) == 3) && pcSubCU->isBipredRestriction(uiPartIdx) )
    {
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMv( TComMv(0,0), ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllRefIdx( -1, ePartSize, uiSubPartIdx, uiDepth, uiPartIdx);
      pcCU->setInterDirSubParts( 1, uiSubPartIdx, uiPartIdx, uiDepth);
    }
#endif
  }
#endif
#endif
  return;
}

/** decode inter direction for a PU block
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \param uiPartIdx
 * \returns Void
 */
Void TDecEntropy::decodeInterDirPU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx )
{
  UInt uiInterDir;

  if ( pcCU->getSlice()->isInterP() )
  {
    uiInterDir = 1;
  }
  else
  {
    m_pcEntropyDecoderIf->parseInterDir( pcCU, uiInterDir, uiAbsPartIdx );
  }

  pcCU->setInterDirSubParts( uiInterDir, uiAbsPartIdx, uiPartIdx, uiDepth );
}

Void TDecEntropy::decodeRefFrmIdxPU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList )
{
  Int iRefFrmIdx = 0;
  Int iParseRefFrmIdx = pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList );

  if ( pcCU->getSlice()->getNumRefIdx( eRefList ) > 1 && iParseRefFrmIdx )
  {
    m_pcEntropyDecoderIf->parseRefFrmIdx( pcCU, iRefFrmIdx, eRefList );
  }
  else if ( !iParseRefFrmIdx )
  {
    iRefFrmIdx = NOT_VALID;
  }
  else
  {
    iRefFrmIdx = 0;
  }

#if JVET_C0024_QTBT
  pcCU->getCUMvField( eRefList )->setAllRefIdx( iRefFrmIdx, SIZE_2Nx2N, uiAbsPartIdx, uiDepth, uiPartIdx );
#else
  PartSize ePartSize = pcCU->getPartitionSize( uiAbsPartIdx );
  pcCU->getCUMvField( eRefList )->setAllRefIdx( iRefFrmIdx, ePartSize, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

/** decode motion vector difference for a PU block
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \param uiPartIdx
 * \param eRefList
 * \returns Void
 */
Void TDecEntropy::decodeMvdPU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList )
{
  if ( pcCU->getInterDir( uiAbsPartIdx ) & ( 1 << eRefList ) )
  {
#if COM16_C1016_AFFINE
    if ( pcCU->isAffine(uiAbsPartIdx) )
    {
      m_pcEntropyDecoderIf->parseAffineMvd( pcCU, uiAbsPartIdx, uiPartIdx, uiDepth, eRefList );
    }
    else
#endif
    m_pcEntropyDecoderIf->parseMvd( pcCU, uiAbsPartIdx, uiPartIdx, uiDepth, eRefList );
  }
}

Void TDecEntropy::decodeMVPIdxPU( TComDataCU* pcSubCU, UInt uiPartAddr, UInt uiDepth, UInt uiPartIdx, RefPicList eRefList )
{
  Int iMVPIdx = -1;

  TComMv cZeroMv( 0, 0 );
  TComMv cMv     = cZeroMv;
  Int    iRefIdx = -1;

  TComCUMvField* pcSubCUMvField = pcSubCU->getCUMvField( eRefList );
  AMVPInfo* pAMVPInfo = pcSubCUMvField->getAMVPInfo();

  iRefIdx = pcSubCUMvField->getRefIdx(uiPartAddr);
  cMv = cZeroMv;

#if COM16_C1016_AFFINE && !VCEG_AZ07_FRUC_MERGE
  if ( pcSubCU->isAffine(uiPartAddr) )
  {
    TComMv acMv[3];
    memset( acMv, 0, sizeof(TComMv) * 3 );

#if JVECT_C0062_AFFINE_SIX_PARAM
    AffineAMVPInfo* pAffineAMVPInfo = pcSubCUMvField->getAffineAMVPInfo(b6Param);
    pcSubCU->fillAffineMvpCand(uiPartIdx, uiPartAddr, eRefList, iRefIdx, pAffineAMVPInfo, b6Param);
#else
    AffineAMVPInfo* pAffineAMVPInfo = pcSubCUMvField->getAffineAMVPInfo();
    pcSubCU->fillAffineMvpCand( uiPartIdx, uiPartAddr, eRefList, iRefIdx, pAffineAMVPInfo );
#endif
    pcSubCU->setMVPNumSubParts( pAffineAMVPInfo->iN, eRefList, uiPartAddr, uiPartIdx, uiDepth );

    if ( (pcSubCU->getInterDir(uiPartAddr) & ( 1 << eRefList )) )
    {
      m_pcEntropyDecoderIf->parseMVPIdx( iMVPIdx );
    }

    pcSubCU->setMVPIdxSubParts( iMVPIdx, eRefList, uiPartAddr, uiPartIdx, uiDepth );

    if ( iRefIdx >= 0 )
    {
#if JVECT_C0062_AFFINE_SIX_PARAM
      m_pcPrediction->getMvPredAffineAMVP(pcSubCU, uiPartIdx, uiPartAddr, eRefList, acMv, b6Param);
#else
      m_pcPrediction->getMvPredAffineAMVP( pcSubCU, uiPartIdx, uiPartAddr, eRefList, acMv );
#endif

      UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB, uiAbsIndexInLCU;
      uiAbsIndexInLCU = pcSubCU->getZorderIdxInCtu();
      pcSubCU->deriveLeftRightTopIdx( uiPartIdx, uiPartIdxLT, uiPartIdxRT );
      pcSubCU->deriveLeftBottomIdx( uiPartIdx, uiPartIdxLB );

      acMv[0] += pcSubCUMvField->getMvd( uiPartIdxLT - uiAbsIndexInLCU );
      acMv[1] += pcSubCUMvField->getMvd( uiPartIdxRT - uiAbsIndexInLCU );

#if JVECT_C0062_AFFINE_SIX_PARAM
      if (b6Param)
      {
        acMv[2] += pcSubCUMvField->getMvd(uiPartIdxLB - uiAbsIndexInLCU);
      }
      else
      {
        Int iWidth = pcSubCU->getWidth(uiPartAddr);
        Int iHeight = pcSubCU->getHeight(uiPartAddr);
        Int vx2 = -(acMv[1].getVer() - acMv[0].getVer()) * iHeight / iWidth + acMv[0].getHor();
        Int vy2 = (acMv[1].getHor() - acMv[0].getHor()) * iHeight / iWidth + acMv[0].getVer();
        acMv[2].set(vx2, vy2);
      }
#else
      Int iWidth = pcSubCU->getWidth(uiPartAddr);
      Int iHeight = pcSubCU->getHeight(uiPartAddr);
      Int vx2 =  - ( acMv[1].getVer() - acMv[0].getVer() ) * iHeight / iWidth + acMv[0].getHor();
      Int vy2 =    ( acMv[1].getHor() - acMv[0].getHor() ) * iHeight / iWidth + acMv[0].getVer();
      acMv[2].set( vx2, vy2 );
#endif

#if !AFFINE_CLIP_REMOVE
      pcSubCU->clipMv(acMv[0]);
      pcSubCU->clipMv(acMv[1]);
      pcSubCU->clipMv(acMv[2]);
#endif
    }

#if JVECT_C0062_AFFINE_SIX_PARAM
    pcSubCU->setAllAffineMv(uiPartAddr, uiPartIdx, acMv, eRefList, 0, b6Param);
#else
    pcSubCU->setAllAffineMv( uiPartAddr, uiPartIdx, acMv, eRefList, 0 );
#endif

    return;
  }
#endif

  if ( (pcSubCU->getInterDir(uiPartAddr) & ( 1 << eRefList )) )
  {
    m_pcEntropyDecoderIf->parseMVPIdx( iMVPIdx );
  }
#if VCEG_AZ07_FRUC_MERGE
  pAMVPInfo->iN = 2;
#else
  pcSubCU->fillMvpCand(uiPartIdx, uiPartAddr, eRefList, iRefIdx, pAMVPInfo);
#endif
  pcSubCU->setMVPNumSubParts(pAMVPInfo->iN, eRefList, uiPartAddr, uiPartIdx, uiDepth);
  pcSubCU->setMVPIdxSubParts( iMVPIdx, eRefList, uiPartAddr, uiPartIdx, uiDepth );
  if ( iRefIdx >= 0 )
  {
    m_pcPrediction->getMvPredAMVP( pcSubCU, uiPartIdx, uiPartAddr, eRefList, cMv);
    cMv += pcSubCUMvField->getMvd( uiPartAddr );
  }

#if JVET_C0024_QTBT
  pcSubCU->getCUMvField( eRefList )->setAllMv(cMv, SIZE_2Nx2N, uiPartAddr, 0, uiPartIdx);
#else
  PartSize ePartSize = pcSubCU->getPartitionSize( uiPartAddr );
  pcSubCU->getCUMvField( eRefList )->setAllMv(cMv, ePartSize, uiPartAddr, 0, uiPartIdx);
#endif
}

#if JVET_C0024_QTBT
Void TDecEntropy::xDecodeTransform        ( Bool& bCodeDQP, Bool& isChromaQpAdjCoded, TComTU &rTu, ComponentID compID)
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  const Bool bChroma = isChromaEnabled(pcCU->getPic()->getChromaFormat());
  const UInt uiDepth=rTu.GetTransformDepthTotal();
  assert( 0 == rTu.GetTransformDepthRel());

  if( (!pcCU->isIntra(uiAbsPartIdx)) && isLuma(compID) && ((!bChroma) || (!pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cb, 0 ) && !pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cr, 0 )) ))
  {
    pcCU->setCbfSubParts( 1, COMPONENT_Y, uiAbsPartIdx, uiDepth );
  }
  else
  {
    m_pcEntropyDecoderIf->parseQtCbf( rTu, compID, true );
  }

  UChar cbf = pcCU->getCbf( uiAbsPartIdx, compID , 0 );

  if (cbf)
  {
#if COM16_C806_EMT
    if (isLuma(compID))
    {
      m_pcEntropyDecoderIf->parseEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, true );
    }
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    Bool dummyCbf;
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    Int  dummyNzTs;
#endif
    m_pcEntropyDecoderIf->parseCoeffNxN( rTu, compID, dummyCbf
#if JVET_C0045_C0053_NO_NSST_FOR_TS
      , dummyNzTs
#endif
      );
#else
    m_pcEntropyDecoderIf->parseCoeffNxN( rTu, compID);
#endif
  }
  //need to add qp coding, crosscomponentpred,... JCA
}
#else

Void TDecEntropy::xDecodeTransform        ( Bool& bCodeDQP, Bool& isChromaQpAdjCoded, TComTU &rTu, const Int quadtreeTULog2MinSizeInCU 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    )
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  const UInt uiDepth=rTu.GetTransformDepthTotal();
  const UInt uiTrDepth = rTu.GetTransformDepthRel();

  UInt uiSubdiv;
  const UInt numValidComponent = pcCU->getPic()->getNumberValidComponents();
  const Bool bChroma = isChromaEnabled(pcCU->getPic()->getChromaFormat());

  const UInt uiLog2TrafoSize = rTu.GetLog2LumaTrSize();
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  if (bDebugRQT)
  {
    printf("x..codeTransform: offsetLuma=%d offsetChroma=%d absPartIdx=%d, uiDepth=%d\n width=%d, height=%d, uiTrIdx=%d, uiInnerQuadIdx=%d\n",
        rTu.getCoefficientOffset(COMPONENT_Y), rTu.getCoefficientOffset(COMPONENT_Cb), uiAbsPartIdx, uiDepth, rTu.getRect(COMPONENT_Y).width, rTu.getRect(COMPONENT_Y).height, rTu.GetTransformDepthRel(), rTu.GetSectionNumber());
    fflush(stdout);
  }
#endif

  if( pcCU->isIntra(uiAbsPartIdx) && pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_NxN && uiDepth == pcCU->getDepth(uiAbsPartIdx) )
  {
    uiSubdiv = 1;
  }
  else if( (pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) && (pcCU->isInter(uiAbsPartIdx)) && ( pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N ) && (uiDepth == pcCU->getDepth(uiAbsPartIdx)) )
  {
    uiSubdiv = (uiLog2TrafoSize >quadtreeTULog2MinSizeInCU);
  }
  else if( uiLog2TrafoSize > pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() )
  {
    uiSubdiv = 1;
  }
  else if( uiLog2TrafoSize == pcCU->getSlice()->getSPS()->getQuadtreeTULog2MinSize() )
  {
    uiSubdiv = 0;
  }
  else if( uiLog2TrafoSize == quadtreeTULog2MinSizeInCU )
  {
    uiSubdiv = 0;
  }
  else
  {
    assert( uiLog2TrafoSize > quadtreeTULog2MinSizeInCU );
#if COM16_C806_T64
    m_pcEntropyDecoderIf->parseTransformSubdivFlag( uiSubdiv, pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrafoSize );
#else
    m_pcEntropyDecoderIf->parseTransformSubdivFlag( uiSubdiv, 5 - uiLog2TrafoSize );
#endif
  }

  for(Int chan=COMPONENT_Cb; chan<numValidComponent; chan++)
  {
    const ComponentID compID=ComponentID(chan);
    const UInt trDepthTotalAdj=rTu.GetTransformDepthTotalAdj(compID);

    const Bool bFirstCbfOfCU = uiTrDepth == 0;

    if( bFirstCbfOfCU )
    {
      pcCU->setCbfSubParts( 0, compID, rTu.GetAbsPartIdxTU(compID), trDepthTotalAdj);
    }
    if( bFirstCbfOfCU || rTu.ProcessingAllQuadrants(compID) )
    {
      if( bFirstCbfOfCU || pcCU->getCbf( uiAbsPartIdx, compID, uiTrDepth - 1 ) )
      {
        m_pcEntropyDecoderIf->parseQtCbf( rTu, compID, (uiSubdiv == 0) );
      }
    }
  }

  if( uiSubdiv )
  {
    const UInt uiQPartNum = pcCU->getPic()->getNumPartitionsInCtu() >> ((uiDepth+1) << 1);
    UInt uiYUVCbf[MAX_NUM_COMPONENT] = {0,0,0};

    TComTURecurse tuRecurseChild(rTu, true);

#if COM16_C806_EMT
    if( 0==uiTrDepth )
    {
      m_pcEntropyDecoderIf->parseEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, true );
    }
#endif

    do
    {
      xDecodeTransform( bCodeDQP, isChromaQpAdjCoded, tuRecurseChild, quadtreeTULog2MinSizeInCU 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    ,  bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
      UInt childTUAbsPartIdx=tuRecurseChild.GetAbsPartIdxTU();
      for(UInt ch=0; ch<numValidComponent; ch++)
      {
        uiYUVCbf[ch] |= pcCU->getCbf(childTUAbsPartIdx , ComponentID(ch),  uiTrDepth+1 );
      }
    } while (tuRecurseChild.nextSection(rTu) );

    for(UInt ch=0; ch<numValidComponent; ch++)
    {
      UChar *pBase = pcCU->getCbf( ComponentID(ch) ) + uiAbsPartIdx;
      const UChar flag = uiYUVCbf[ch] << uiTrDepth;

      for( UInt ui = 0; ui < 4 * uiQPartNum; ++ui )
      {
        pBase[ui] |= flag;
      }
    }
  }
  else
  {
    assert( uiDepth >= pcCU->getDepth( uiAbsPartIdx ) );
    pcCU->setTrIdxSubParts( uiTrDepth, uiAbsPartIdx, uiDepth );

    {
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tTrIdx: abspart=" );
      DTRACE_CABAC_V( uiAbsPartIdx );
      DTRACE_CABAC_T( "\tdepth=" );
      DTRACE_CABAC_V( uiDepth );
      DTRACE_CABAC_T( "\ttrdepth=" );
      DTRACE_CABAC_V( uiTrDepth );
      DTRACE_CABAC_T( "\n" );
    }

    pcCU->setCbfSubParts ( 0, COMPONENT_Y, uiAbsPartIdx, uiDepth );

    if( (!pcCU->isIntra(uiAbsPartIdx)) && uiDepth == pcCU->getDepth( uiAbsPartIdx ) && ((!bChroma) || (!pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cb, 0 ) && !pcCU->getCbf( uiAbsPartIdx, COMPONENT_Cr, 0 )) ))
    {
      pcCU->setCbfSubParts( 1 << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, uiDepth );
    }
    else
    {
      m_pcEntropyDecoderIf->parseQtCbf( rTu, COMPONENT_Y, true );
    }


    // transform_unit begin
    UInt cbf[MAX_NUM_COMPONENT]={0,0,0};
    Bool validCbf       = false;
    Bool validChromaCbf = false;
    const UInt uiTrIdx = rTu.GetTransformDepthRel();

    for(UInt ch=0; ch<pcCU->getPic()->getNumberValidComponents(); ch++)
    {
      const ComponentID compID = ComponentID(ch);

      cbf[compID] = pcCU->getCbf( uiAbsPartIdx, compID, uiTrIdx );

      if (cbf[compID] != 0)
      {
        validCbf = true;
        if (isChroma(compID))
        {
          validChromaCbf = true;
        }
      }
    }

#if COM16_C806_EMT
    if( 0==uiTrDepth )
    {
      m_pcEntropyDecoderIf->parseEmtCuFlag( pcCU, uiAbsPartIdx, uiDepth, cbf[COMPONENT_Y] ? true : false );
    }
#endif

    if ( validCbf )
    {

      // dQP: only for CTU
      if ( pcCU->getSlice()->getPPS()->getUseDQP() )
      {
        if ( bCodeDQP )
        {
          const UInt uiAbsPartIdxCU=rTu.GetAbsPartIdxCU();
#if BLOCK_BASED_BACKGROUND_REFERENCE && !DISABLE_DQP_SIGNALING
					if (pcCU->getLCUBgFlag())
					{
						pcCU->setQPSubParts(pcCU->getSlice()->getSliceBgQp(), uiAbsPartIdxCU, uiDepth);
						pcCU->setCodedQP(pcCU->getSlice()->getSliceBgQp());
					}
					else
					{
						pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), uiAbsPartIdxCU, uiDepth);
						pcCU->setCodedQP(pcCU->getSlice()->getSliceQp());
					}
#else
          decodeQP( pcCU, uiAbsPartIdxCU);
#endif
          bCodeDQP = false;
        }
      }

      if ( pcCU->getSlice()->getUseChromaQpAdj() )
      {
        if ( validChromaCbf && isChromaQpAdjCoded && !pcCU->getCUTransquantBypass(rTu.GetAbsPartIdxCU()) )
        {
          decodeChromaQpAdjustment( pcCU, rTu.GetAbsPartIdxCU() );
          isChromaQpAdjCoded = false;
        }
      }

      const UInt numValidComp=pcCU->getPic()->getNumberValidComponents();

      for(UInt ch=COMPONENT_Y; ch<numValidComp; ch++)
      {
        const ComponentID compID=ComponentID(ch);

        if( rTu.ProcessComponentSection(compID) )
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
              const UInt subTUCBF = pcCU->getCbf(subTUIterator.GetAbsPartIdxTU(), compID, (uiTrIdx + 1));

              if (subTUCBF != 0)
              {
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
                if (bDebugRQT)
                {
                  printf("Call NxN for chan %d width=%d height=%d cbf=%d\n", compID, subTUIterator.getRect(compID).width, subTUIterator.getRect(compID).height, 1);
                }
#endif
                m_pcEntropyDecoderIf->parseCoeffNxN( subTUIterator, compID 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    ,  bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
              }
            } while (subTUIterator.nextSection(rTu));
          }
          else
          {
            if(isChroma(compID) && (cbf[COMPONENT_Y] != 0))
            {
              m_pcEntropyDecoderIf->parseCrossComponentPrediction( rTu, compID );
            }

            if(cbf[compID] != 0)
            {
              m_pcEntropyDecoderIf->parseCoeffNxN( rTu, compID 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    ,  bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
            }
          }
        }
      }
    }
    // transform_unit end
  }
}

#endif //#if !JVET_C0024_QTBT
Void TDecEntropy::decodeQP          ( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if ( pcCU->getSlice()->getPPS()->getUseDQP() )
  {
    m_pcEntropyDecoderIf->parseDeltaQP( pcCU, uiAbsPartIdx, pcCU->getDepth( uiAbsPartIdx ) );
  }
}

Void TDecEntropy::decodeChromaQpAdjustment( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if ( pcCU->getSlice()->getUseChromaQpAdj() )
  {
    m_pcEntropyDecoderIf->parseChromaQpAdjustment( pcCU, uiAbsPartIdx, pcCU->getDepth( uiAbsPartIdx ) );
  }
}


//! decode coefficients
Void TDecEntropy::decodeCoeff( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool& bCodeDQP, Bool& isChromaQpAdjCoded )
{
  if( pcCU->isIntra(uiAbsPartIdx) )
  {
  }
  else
  {
    UInt uiQtRootCbf = 1;
#if JVET_C0024_QTBT
    if( !pcCU->getMergeFlag( uiAbsPartIdx ) ) 
#else
    if( !( pcCU->getPartitionSize( uiAbsPartIdx) == SIZE_2Nx2N && pcCU->getMergeFlag( uiAbsPartIdx ) ) )
#endif
    {
      m_pcEntropyDecoderIf->parseQtRootCbf( uiAbsPartIdx, uiQtRootCbf );
    }
    if ( !uiQtRootCbf )
    {
      static const UInt cbfZero[MAX_NUM_COMPONENT]={0,0,0};
      pcCU->setCbfSubParts( cbfZero, uiAbsPartIdx, uiDepth );
#if !JVET_C0024_QTBT
      pcCU->setTrIdxSubParts( 0 , uiAbsPartIdx, uiDepth );
#endif
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
      xDecodeTransform( bCodeDQP, isChromaQpAdjCoded, tuRecurse, compID);
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
    xDecodeTransform( bCodeDQP, isChromaQpAdjCoded, tuRecurse, COMPONENT_Y);
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
#if BLOCK_BASED_BACKGROUND_REFERENCE && !DISABLE_DQP_SIGNALING
		if (pcCU->getLCUBgFlag())
		{
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceBgQp(), uiAbsPartIdx, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx));
			pcCU->setCodedQP(pcCU->getSlice()->getSliceBgQp());
		}
		else
		{
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), uiAbsPartIdx, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx));
			pcCU->setCodedQP(pcCU->getSlice()->getSliceQp());
		}
#else
    decodeQP( pcCU, uiAbsPartIdx );
#endif
    bCodeDQP = false;
  }

  if ( pcCU->getSlice()->getUseChromaQpAdj() )
  {
    if ( validChromaCbf && isChromaQpAdjCoded && !pcCU->getCUTransquantBypass(uiAbsPartIdx) )
    {
      decodeChromaQpAdjustment( pcCU, uiAbsPartIdx );
      isChromaQpAdjCoded = false;
    }
  }
#endif

#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
#if !QTBT_NSST
  Int bCbfCU = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType()) ? pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y)
    : (pcCU->getCbf(uiAbsPartIdx, COMPONENT_Cb) || pcCU->getCbf(uiAbsPartIdx, COMPONENT_Cr) )): pcCU->getQtRootCbf(uiAbsPartIdx);
#else
  Int iNonZeroCoeff = 0;
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
    if( !pcCU->getTransformSkip( uiAbsPartIdx, compID) )
#endif
    iNonZeroCoeff += countNonZeroCoeffs(pcCoef, uiWidth * uiHeight);
  }
#endif
#endif
#else
  Int quadtreeTULog2MinSizeInCU = pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx);
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  Bool bCbfCU = 0;
#endif  
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  Int iNonZeroCoeffNonTs = 0;
#endif
  xDecodeTransform( bCodeDQP, isChromaQpAdjCoded, tuRecurse, quadtreeTULog2MinSizeInCU 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    ,  bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
  if (
#if QTBT_NSST
    iNonZeroCoeff > iNonZeroCoeffThr
#else
    bCbfCU 
#endif
    && isLuma(pcCU->getTextType()))
  {
    decodeROTIdx( pcCU, uiAbsPartIdx, uiDepth );
  }
  else if(
#if QTBT_NSST
    iNonZeroCoeff > iNonZeroCoeffThr
#else
    bCbfCU 
#endif
    && isChroma(pcCU->getTextType()) && pcCU->getWidth(uiAbsPartIdx) >= 8 && pcCU->getHeight(uiAbsPartIdx) >= 8 )
  {
    assert(pcCU->getSlice()->isIntra());
    decodeROTIdxChroma( pcCU, uiAbsPartIdx, uiDepth );
  }
#else
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  if ( iNonZeroCoeffNonTs )
#else
  if (bCbfCU )
#endif
    decodeROTIdx( pcCU, uiAbsPartIdx, uiDepth );
#endif
#endif
#if RRC_INTRA_SHARPENING
    if (pcCU->getIntraSharpFlag(uiAbsPartIdx))
    {
      pcCU->setLumaIntraFilter(uiAbsPartIdx, false);
    }
#endif
}

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
TDecEntropy::TDecEntropy()
{
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
}

TDecEntropy::~TDecEntropy()
{
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
}
#endif

#if QTBT_NSST
Int TDecEntropy::countNonZeroCoeffs( TCoeff* pcCoef, UInt uiSize )
{
  Int count = 0;

  for ( Int i = 0; i < uiSize; i++ )
  {
    count += pcCoef[i] != 0;
  }

  return count;
}
#endif

#if ALF_HM3_REFACTOR
#if JVET_C0038_GALF
Void TDecEntropy::decodeAux(ALFParam* pAlfParam)
{
  UInt uiSymbol;
  Int sqrFiltLengthTab[3] = {TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_9SYM, TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_7SYM, TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_5SYM};
  Int FiltTab[3] = {9, 7, 5};
  Int i; 
  UChar decodetab_pred[3]={1, 0, 2}; 

  pAlfParam->filters_per_group = 0;

  memset (pAlfParam->filterPattern, 0 , sizeof(Int)*TComAdaptiveLoopFilter::m_NO_VAR_BINS);
  memset (pAlfParam->PrevFiltIdx,   0 , sizeof(Char)*TComAdaptiveLoopFilter::m_NO_VAR_BINS);
  pAlfParam->filtNo = 1; //nonZeroCoeffs

  //number of total filters
  UInt iNoVarBins = TComAdaptiveLoopFilter::m_NO_VAR_BINS;
  m_pcEntropyDecoderIf->parseALFTruncBinVal(uiSymbol, iNoVarBins);
  pAlfParam->noFilters = uiSymbol + 1;
  pAlfParam->filters_per_group = pAlfParam->noFilters;


  //Filter tap
 // if (pAlfParam->noFilters >1 || pAlfParam->forceCoeff0==0) //to be checked
  {
  m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
  Int TabIdx = uiSymbol;
  pAlfParam->realfiltNo = 2-TabIdx;
  pAlfParam->tap = FiltTab[pAlfParam->realfiltNo];
  pAlfParam->num_coeff = sqrFiltLengthTab[pAlfParam->realfiltNo];
  }

  //filter set index for each class
  if(pAlfParam->noFilters > 1)
  {
    for (i=0; i< iNoVarBins; i++)
    {
      m_pcEntropyDecoderIf->parseALFTruncBinVal(uiSymbol, (UInt)pAlfParam->noFilters);
      pAlfParam->filterPattern[i] = (Int)uiSymbol;
    }
  }
  else
  {
    memset( pAlfParam->filterPattern, 0, iNoVarBins* sizeof(Int) );
  }
  memcpy(pAlfParam->varIndTab, pAlfParam->filterPattern, TComAdaptiveLoopFilter::m_NO_VAR_BINS * sizeof(int));

  //prediction from fixed filters
  //0: no pred; 1: all same index; 2: diff index for each variance index

  // frame that filter is predicted from
  if (pAlfParam->iAvailableFilters > 0 )
  { 
    // prediction pattern
    m_pcEntropyDecoderIf->parseALFPrevFiltType(uiSymbol); //0: all zero, no pred from pre-defined filters; 1: all are predicted but could be different values; 2: some predicted and some not
    pAlfParam->iPredPattern  = decodetab_pred[uiSymbol];

    if (pAlfParam->iPredPattern == 0)
    {
      memset(pAlfParam->PrevFiltIdx, 0, sizeof(Char)*TComAdaptiveLoopFilter::m_NO_VAR_BINS);
    }
    else 
    {
      if (pAlfParam->iPredPattern == 2)
      {
        //on/off flags
        for (i=0; i<iNoVarBins; i++)
        {    
          m_pcEntropyDecoderIf->parseALFPrevFiltFlag(uiSymbol);
          pAlfParam->PrevFiltIdx[i] = (Char)uiSymbol;
        }
      }
      else         
      {
        assert (pAlfParam->iPredPattern == 1);
        for (i=0; i<iNoVarBins; i++)
        {            
          pAlfParam->PrevFiltIdx[i] = 1;
        }
      }
      if (pAlfParam->iAvailableFilters > 1)
      {
        assert(pAlfParam->iPredPattern > 0 );
        for (i=0; i< iNoVarBins; i++)
        {
          if (pAlfParam->PrevFiltIdx[i] > 0)
          {
            m_pcEntropyDecoderIf->parseALFTruncBinVal(uiSymbol, (UInt)pAlfParam->iAvailableFilters);
            pAlfParam->PrevFiltIdx[i] = (Char)(uiSymbol + 1);
          }
        }
      }
    }
  }
}
#else
Void TDecEntropy::decodeAux(ALFParam* pAlfParam)
{
  UInt uiSymbol;
  Int sqrFiltLengthTab[3] = {TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_9SYM, TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_7SYM, TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_5SYM};
  Int FiltTab[3] = {9, 7, 5};

  pAlfParam->filters_per_group = 0;

  memset (pAlfParam->filterPattern, 0 , sizeof(Int)*TComAdaptiveLoopFilter::m_NO_VAR_BINS);
  pAlfParam->filtNo = 1; //nonZeroCoeffs

  m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
  Int TabIdx = uiSymbol;
  pAlfParam->realfiltNo = 2-TabIdx;
  pAlfParam->tap = FiltTab[pAlfParam->realfiltNo];
  pAlfParam->num_coeff = sqrFiltLengthTab[pAlfParam->realfiltNo];

  if (pAlfParam->filtNo>=0)
  {
    if(pAlfParam->realfiltNo >= 0)
    {
      // filters_per_fr
      m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
      pAlfParam->noFilters = uiSymbol + 1;
      pAlfParam->filters_per_group = pAlfParam->noFilters; 

      if(pAlfParam->noFilters == 2)
      {
        m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
        pAlfParam->startSecondFilter = uiSymbol;
        pAlfParam->filterPattern [uiSymbol] = 1;
      }
      else if (pAlfParam->noFilters > 2)
      {
        pAlfParam->filters_per_group = 1;
        for (int i=1; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++) 
        {
          m_pcEntropyDecoderIf->parseAlfFlag (uiSymbol);
          pAlfParam->filterPattern[i] = uiSymbol;
          pAlfParam->filters_per_group += uiSymbol;
        }
      }
    }
  }
  else
  {
    memset (pAlfParam->filterPattern, 0, TComAdaptiveLoopFilter::m_NO_VAR_BINS*sizeof(Int));
  }
  // Reconstruct varIndTab[]
  memset(pAlfParam->varIndTab, 0, TComAdaptiveLoopFilter::m_NO_VAR_BINS * sizeof(int));
  if(pAlfParam->filtNo>=0)
  {
    for(Int i = 1; i < TComAdaptiveLoopFilter::m_NO_VAR_BINS; ++i)
    {
      if(pAlfParam->filterPattern[i])
        pAlfParam->varIndTab[i] = pAlfParam->varIndTab[i-1] + 1;
      else
        pAlfParam->varIndTab[i] = pAlfParam->varIndTab[i-1];
    }
  }
}
#endif

Void TDecEntropy::readFilterCodingParams(ALFParam* pAlfParam
#if JVET_C0038_GALF
    , Bool bChroma
#endif
  )
{
  UInt uiSymbol;
  int ind, scanPos;
  int golombIndexBit;
  int kMin;
  int maxScanVal;
  const int *pDepthInt;
  int fl;
#if JVET_C0038_GALF
  if( bChroma )
  {
    maxScanVal = 0;
    pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[pAlfParam->tap_chroma==5? 0: (pAlfParam->tap_chroma==7? 1: 2)];
    Int sqrFiltLength= pAlfParam->num_coeff_chroma -1;
    for(ind = 0; ind < sqrFiltLength; ind++)
    {
      maxScanVal = max(maxScanVal, pDepthInt[ind]);
    }
  }
  else
  {
#endif
  // Determine fl
  if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_9SYM)
    fl = 4;
  else if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_7SYM)
    fl = 3;
  else
    fl = 2;

  // Determine maxScanVal
  maxScanVal = 0;
  pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl - 2];
#if JVET_C0038_GALF
  Int sqrFiltLength= pAlfParam->num_coeff -1;
  for(ind = 0; ind < sqrFiltLength; ind++)
#else
  for(ind = 0; ind < pAlfParam->num_coeff; ind++)
#endif
    maxScanVal = max(maxScanVal, pDepthInt[ind]);
#if JVET_C0038_GALF
  }
#endif
  // Golomb parameters
  m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
  pAlfParam->minKStart = 1 + uiSymbol;

  kMin = pAlfParam->minKStart;
  for(scanPos = 0; scanPos < maxScanVal; scanPos++)
  {
    m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
    golombIndexBit = uiSymbol;
    if(golombIndexBit)
      pAlfParam->kMinTab[scanPos] = kMin + 1;
    else
      pAlfParam->kMinTab[scanPos] = kMin;
    kMin = pAlfParam->kMinTab[scanPos];
  }
#if JVET_C0038_GALF
  if(!bChroma)
  {
    if(pAlfParam->forceCoeff0 )
    {
      for(ind=0; ind<pAlfParam->filters_per_group; ++ind)
      {
        m_pcEntropyDecoderIf->parseALFPrevFiltFlag(uiSymbol);
        pAlfParam->codedVarBins[ind] = (uiSymbol ? true: false);
}
    }
    else
    {
      for(ind = 0; ind < pAlfParam->filters_per_group; ind++)
      {
        pAlfParam->codedVarBins[ind] = true;
      }
    }
  }
#endif
}

Int TDecEntropy::golombDecode(Int k)
{
  UInt uiSymbol;
  Int q = -1;
  Int nr = 0;
  Int m = (Int)pow(2.0, k);
  Int a;

  uiSymbol = 1;
  while (uiSymbol)
  {
    m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
    q++;
  }
  for(a = 0; a < k; ++a)          // read out the sequential log2(M) bits
  {
    m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
    if(uiSymbol)
      nr += 1 << a;
  }
  nr += q * m;                    // add the bits and the multiple of M
  if(nr != 0)
  {
    m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
    nr = (uiSymbol)? nr: -nr;
  }
  return nr;
}


#if JVET_C0038_GALF
Void TDecEntropy::readFilterCoeffs(ALFParam* pAlfParam, Bool bChroma )
{
  int ind, scanPos, i;
  const int *pDepthInt;
  int fl;
  if( bChroma )
  {
    fl = pAlfParam->tap_chroma==9 ? 4 :(pAlfParam->tap_chroma==7 ? 3: 2);
    pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl - 2];
    Int iNumCoeffMinus1 =  pAlfParam->num_coeff_chroma - 1 ;
    for(i = 0; i < iNumCoeffMinus1; i++)
    {
      scanPos = pDepthInt[i] - 1;
      pAlfParam->coeff_chroma [i] = golombDecode(pAlfParam->kMinTab[scanPos]);
    }
  }
  else
  {
    if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_9SYM)
      fl = 4;
    else if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_7SYM)
      fl = 3;
    else
      fl = 2;
    pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl - 2];
    for(ind = 0; ind < pAlfParam->filters_per_group_diff; ++ind)
    {
      if(!pAlfParam->codedVarBins[ind])
      {
        continue;
      }
      Int iNumCoeffMinus1 =  pAlfParam->num_coeff - 1 ;
      for(i = 0; i < iNumCoeffMinus1; i++)
      {
        scanPos = pDepthInt[i] - 1;
        pAlfParam->coeffmulti[ind][i] = golombDecode(pAlfParam->kMinTab[scanPos]);
      }
    }
   }
}
#else
Void TDecEntropy::readFilterCoeffs(ALFParam* pAlfParam)
{
  int ind, scanPos, i;
  const int *pDepthInt;
  int fl;

  if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_9SYM)
    fl = 4;
  else if(pAlfParam->num_coeff == TComAdaptiveLoopFilter::m_SQR_FILT_LENGTH_7SYM)
    fl = 3;
  else
    fl = 2;

  pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl - 2];

  for(ind = 0; ind < pAlfParam->filters_per_group_diff; ++ind)
  {
    for(i = 0; i < pAlfParam->num_coeff; i++)
    {
      scanPos = pDepthInt[i] - 1;
      pAlfParam->coeffmulti[ind][i] = golombDecode(pAlfParam->kMinTab[scanPos]);
    }
  }
}
#endif

#if JVET_C0038_GALF
Void TDecEntropy::decodeFilterCoeff (ALFParam* pAlfParam, Bool bChroma)
{
  readFilterCodingParams (pAlfParam, bChroma);
  readFilterCoeffs (pAlfParam, bChroma);
}
#else
Void TDecEntropy::decodeFilterCoeff (ALFParam* pAlfParam)
{
  readFilterCodingParams (pAlfParam);
  readFilterCoeffs (pAlfParam);
}
#endif

Void TDecEntropy::decodeFilt(ALFParam* pAlfParam)
{
  UInt uiSymbol;

  if (pAlfParam->filtNo >= 0)
  {
    pAlfParam->filters_per_group_diff = pAlfParam->filters_per_group;
#if JVET_C0038_GALF      
    m_pcEntropyDecoderIf->parseALFPrevFiltFlag(uiSymbol); //force0 flag
    pAlfParam->forceCoeff0 = uiSymbol ? true : false;
    if( !pAlfParam->forceCoeff0 && pAlfParam->filters_per_group > 1 )  
    {
      m_pcEntropyDecoderIf->parseAlfFlag (uiSymbol);
      pAlfParam->predMethod = uiSymbol;
    }
    else
    {      
      pAlfParam->predMethod = 0;
    }
#else
    if (pAlfParam->filters_per_group > 1)
    {
      pAlfParam->forceCoeff0 = 0;
      {
        for (int i=0; i<TComAdaptiveLoopFilter::m_NO_VAR_BINS; i++)
          pAlfParam->codedVarBins[i] = 1;

      }
      m_pcEntropyDecoderIf->parseAlfFlag (uiSymbol);
      pAlfParam->predMethod = uiSymbol;
    }
    else
    {
      pAlfParam->forceCoeff0 = 0;
      pAlfParam->predMethod = 0;
    }
#endif

    decodeFilterCoeff (pAlfParam);
  }
}

Void TDecEntropy::decodeAlfParam(ALFParam* pAlfParam, UInt uiMaxTotalCUDepth
#if FIX_TICKET12
        , TComSlice * pSlice
#endif
  )
{
  UInt uiSymbol;
#if !JVET_C0038_GALF 
  Int iSymbol;
#endif
  m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
  pAlfParam->alf_flag = uiSymbol;

  if (!pAlfParam->alf_flag)
  {
    return;
  }
#if COM16_C806_ALF_TEMPPRED_NUM
  //encode temporal reuse flag and index
#if FIX_TICKET12
  if( pSlice->getSliceType() == I_SLICE)
  {
    uiSymbol = 0;
  }
  else
#endif
  m_pcEntropyDecoderIf->parseAlfFlag( uiSymbol );
  pAlfParam->temproalPredFlag = uiSymbol;
  if( uiSymbol )
  {
    m_pcEntropyDecoderIf->parseAlfUvlc( uiSymbol );
    pAlfParam->prevIdx = uiSymbol;
  }
#endif
#if !JVET_C0038_GALF 
  Int pos;
#endif

#if COM16_C806_ALF_TEMPPRED_NUM
  if( !pAlfParam->temproalPredFlag )
  {
#endif
    decodeAux(pAlfParam);
    decodeFilt(pAlfParam);
#if COM16_C806_ALF_TEMPPRED_NUM
  }
#endif
  // filter parameters for chroma
  m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
  pAlfParam->chroma_idc = uiSymbol;
#if COM16_C806_ALF_TEMPPRED_NUM
  if(!pAlfParam->temproalPredFlag && pAlfParam->chroma_idc)
#else
  if(pAlfParam->chroma_idc)
#endif
  {
#if JVET_C0038_GALF 
    pAlfParam->num_coeff_chroma = ((pAlfParam->tap_chroma*pAlfParam->tap_chroma)>>2) + 1;
    decodeFilterCoeff(pAlfParam, true);
#else
    m_pcEntropyDecoderIf->parseAlfUvlc(uiSymbol);
    pAlfParam->tap_chroma = (uiSymbol<<1) + 5;
    pAlfParam->num_coeff_chroma = ((pAlfParam->tap_chroma*pAlfParam->tap_chroma+1)>>1) + 1;

    // filter coefficients for chroma
    for(pos=0; pos<pAlfParam->num_coeff_chroma; pos++)
    {
      m_pcEntropyDecoderIf->parseAlfSvlc(iSymbol);
      pAlfParam->coeff_chroma[pos] = iSymbol;
    }
#endif
  }

  // region control parameters for luma
  m_pcEntropyDecoderIf->parseAlfFlag(uiSymbol);
  pAlfParam->cu_control_flag = uiSymbol;
  if (pAlfParam->cu_control_flag)
  {
    m_pcEntropyDecoderIf->parseAlfCtrlDepth(uiSymbol,uiMaxTotalCUDepth);
    pAlfParam->alf_max_depth = uiSymbol;
    decodeAlfCtrlParam(pAlfParam);
  }
}

Void TDecEntropy::decodeAlfCtrlParam( ALFParam* pAlfParam )
{
  UInt uiSymbol;
  m_pcEntropyDecoderIf->parseAlfFlagNum( uiSymbol, pAlfParam->num_cus_in_frame, pAlfParam->alf_max_depth );
  pAlfParam->num_alf_cu_flag = uiSymbol;

  for(UInt i=0; i<pAlfParam->num_alf_cu_flag; i++)
  {
    m_pcEntropyDecoderIf->parseAlfCtrlFlag( pAlfParam->alf_cu_flag[i] );
  }
}

#endif

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
Void TDecEntropy::updateStates   ( SliceType uiSliceType, UInt uiSliceQP, TComStats*  apcStats)
{
  Int k;
  for (k = 0; k < NUM_QP_PROB; k++)
  {
    if (apcStats-> aaQPUsed[uiSliceType][k].used ==true && apcStats-> aaQPUsed[uiSliceType][k].QP == uiSliceQP)
    {
      apcStats-> aaQPUsed[uiSliceType][k].firstUsed = false;
      break;
    }
    else if (apcStats-> aaQPUsed[uiSliceType][k].used ==false)
    {
      apcStats-> aaQPUsed[uiSliceType][k].used = true;
      apcStats-> aaQPUsed[uiSliceType][k].QP = uiSliceQP;
      apcStats-> aaQPUsed[uiSliceType][k].firstUsed = true;
      break;
    }
  }
}
#endif

#if RRC_NS
UInt TDecEntropy::decodeNsFlag()
{
	return m_pcEntropyDecoderIf->parseNsFlag();
}

UInt TDecEntropy::decodeAdaptiveNsFlag()
{
	return m_pcEntropyDecoderIf->parseAdaptiveNsFlag();
}

UInt TDecEntropy::decodeNsFWholeFrameFlag()
{
	return m_pcEntropyDecoderIf->parseNsFWholeFrameFlag();
}

Int TDecEntropy::decodeNsBlockSize()
{
	return m_pcEntropyDecoderIf->parseNsBlockSize();
}
#endif

#if COM16_C1016_AFFINE
Void TDecEntropy::decodeAffineFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx )
{
  if ( !pcCU->getSlice()->getSPS()->getUseAffine() )
  {
    return;
  }

  m_pcEntropyDecoderIf->parseAffineFlag( pcCU, uiAbsPartIdx, uiDepth, uiPuIdx );
}
#endif

#if MRG_CONSIDER_COMPLEX
Void TDecEntropy::decodeMergeComplexFlag(TComDataCU* pcSubCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx)
{
  if (!pcSubCU->getSlice()->getSPS()->getUseMrgComplex())
  {
    return;
  }

  m_pcEntropyDecoderIf->parseMergeComplexFlag(pcSubCU, uiAbsPartIdx, uiDepth, uiPUIdx);
}
Void TDecEntropy::decodeMergeComplexIndex(TComDataCU* pcCU, UInt uiPartIdx, UInt uiAbsPartIdx, UInt uiDepth)
{
  UInt uiMergeComplexIndex = 0;
  m_pcEntropyDecoderIf->parseMergeComplexIndex(pcCU, uiAbsPartIdx, uiMergeComplexIndex);
#if JVET_C0024_QTBT
  pcCU->setMergeComplexIndexSubParts(uiMergeComplexIndex, uiAbsPartIdx);
#else
  pcCU->setMergeComplexIndexSubParts(uiMergeComplexIndex, uiAbsPartIdx, uiPartIdx, uiDepth);
#endif
}
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
Void TDecEntropy::decodeAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx)
{
  if (!pcCU->getSlice()->getSPS()->getUseAffine())
  {
    return;
  }

  m_pcEntropyDecoderIf->parseAffineParamFlag(pcCU, uiAbsPartIdx, uiDepth, uiPuIdx);
}

#endif

#if RRC_INTRA_SHARPENING
Void TDecEntropy::decodeIntraSharpFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  UInt sharpFlag = 0;
  if (pcCU->isIntraSharpFlagCoded(uiAbsPartIdx))
    sharpFlag = m_pcEntropyDecoderIf->parseIntraSharpFlag();
  pcCU->setIntraSharpFlagSubParts(sharpFlag, uiAbsPartIdx, uiDepth);
}
#endif


#if RRC_SELECTABLE_IF
Void TDecEntropy::decodeAdaptiveIfIndex(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  Int idx = 0;
#if SIF_CU_LEVEL
  if (pcCU->isSIFFlagCoded(uiAbsPartIdx))
#endif
  {
    idx = m_pcEntropyDecoderIf->parseAdaptiveIfIndex();
  }
  pcCU->setSIFFlagSubParts(idx, uiAbsPartIdx, uiDepth);
}
#endif

//! \}
