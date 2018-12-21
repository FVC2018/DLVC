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

/** \file     TEncGOP.cpp
    \brief    GOP encoder class
*/

#include <list>
#include <algorithm>
#include <functional>

#include "TEncTop.h"
#include "TEncGOP.h"
#include "TEncAnalyze.h"
#include "libmd5/MD5.h"
#include "TLibCommon/SEI.h"
#include "TLibCommon/NAL.h"
#include "NALwrite.h"
#if RRC_SELECTABLE_IF
#include "TLibCommon/TComSharpeningFilter.h"
#endif
#include <time.h>
#include <math.h>

#include <deque>
#if RRC_NS
#include "TLibCommon/TComNoiseSuppressor.h"
#endif
using namespace std;

#if JVET_D0033_ADAPTIVE_CLIPPING
namespace {
// simulate a coding decoding of the Intra(intrinsic) prms: to be sync with real encoding
Bound encodeDecodeBoundIntra(Bound b,Int n,Int bd) {

    b.m>>=ClipParam::cquantiz;
    b.m=std::min((1<<(n-ClipParam::cquantiz))-1,b.m);
    b.m<<=ClipParam::cquantiz;
    b.m=Clip3(0,(1<<bd)-1,b.m);

    b.M>>=ClipParam::cquantiz;
    b.M=std::min((1<<(n-ClipParam::cquantiz))-1,b.M);
    b.M<<=ClipParam::cquantiz;
    b.M=Clip3(0,(1<<bd)-1,b.M);

    return b;
}

ClipParam encodeDecodeClipPrmIntra(const ClipParam &prm) {
    ClipParam prm2;

    prm2.isActive=prm.isActive;
    prm2.isChromaActive=prm.isChromaActive;
    prm2.Y()=encodeDecodeBoundIntra(prm.Y(),ClipParam::nbBitsY,ClipParam::ibdLuma);
    prm2.U()=encodeDecodeBoundIntra(prm.U(),ClipParam::nbBitsUV,ClipParam::ibdChroma);
    prm2.V()=encodeDecodeBoundIntra(prm.V(),ClipParam::nbBitsUV,ClipParam::ibdChroma);

    return prm2;

}

ClipParam codingChoice(const ClipParam &prm_tocode,
                          Int delta_disto_luma,Int delta_disto_chroma,
                          Double lambda_luma,Double lambda_chroma,int tbd
                          )
{
    ClipParam prm;

    // intra coding
    prm=encodeDecodeClipPrmIntra(prm_tocode);
#if JVET_D0033_ADAPTIVE_CLIPPING_ENC_METHOD    
    Int dR_luma_only;   // delta cost with off
    Int dR_chroma; // delta cost with off

    dR_luma_only=(ClipParam::nbBitsY-ClipParam::cquantiz)*2     +1 /* chroma off */;
    dR_chroma=(ClipParam::nbBitsUV-ClipParam::cquantiz)*2    *2 /* chroma on */;

    double lambda_coef=2;

    // chroma are set for L2 norm -> use sqrt to correct for L0/L1
    lambda_luma=lambda_coef*sqrt(lambda_luma);
    lambda_chroma=lambda_coef*sqrt(lambda_chroma);
    // then rdo
    // D_on+lamda*R_on < D_off+lamda*R_off -> (D_on-D_off)+lambda*dR<0 -> -dD + lambda * dR <0

    // only luma 0
    Bool chroma_on=false;
    Double best_rd=-delta_disto_luma  +lambda_luma   * dR_luma_only;

    // add chroma
    Double rd=-delta_disto_luma  +lambda_luma   *dR_luma_only-delta_disto_chroma+lambda_chroma*dR_chroma;
    if (rd<best_rd) {
        best_rd=rd;
        chroma_on=true;
    }

    if (best_rd>0.) {
        prm.isActive=false;
    } else {
        prm.isChromaActive=chroma_on;
        if (!chroma_on) { // replace chroma by default values
            prm.U().m=0;
            prm.U().M=(1<<ClipParam::ibdChroma)-1;
            prm.V()=prm.U();
        }
    }
#else
    prm.isActive=true;
    prm.isChromaActive=true;
#endif
    return prm;

}


}
#endif
//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / initialization / destroy
// ====================================================================================================================
Int getLSB(Int poc, Int maxLSB)
{
  if (poc >= 0)
  {
    return poc % maxLSB;
  }
  else
  {
    return (maxLSB - ((-poc) % maxLSB)) % maxLSB;
  }
}

TEncGOP::TEncGOP()
{
  m_iLastIDR            = 0;
  m_iGopSize            = 0;
  m_iNumPicCoded        = 0; //Niko
  m_bFirst              = true;
  m_iLastRecoveryPicPOC = 0;

  m_pcCfg               = NULL;
  m_pcSliceEncoder      = NULL;
  m_pcListPic           = NULL;

  m_pcEntropyCoder      = NULL;
  m_pcCavlcCoder        = NULL;
  m_pcSbacCoder         = NULL;
  m_pcBinCABAC          = NULL;

  m_bSeqFirst           = true;

  m_bRefreshPending     = 0;
  m_pocCRA            = 0;
  m_numLongTermRefPicSPS = 0;
  ::memset(m_ltRefPicPocLsbSps, 0, sizeof(m_ltRefPicPocLsbSps));
  ::memset(m_ltRefPicUsedByCurrPicFlag, 0, sizeof(m_ltRefPicUsedByCurrPicFlag));
  m_lastBPSEI         = 0;
  m_bufferingPeriodSEIPresentInAU = false;
  m_associatedIRAPType = NAL_UNIT_CODED_SLICE_IDR_N_LP;
  m_associatedIRAPPOC  = 0;
  return;
}

TEncGOP::~TEncGOP()
{
}

#if COM16_C806_ALF_TEMPPRED_NUM
#if JVET_E0104_ALF_TEMP_SCALABILITY
Int TEncGOP::m_iStoredAlfParaNum[JVET_E0104_ALF_MAX_TEMPLAYERID] = { 0, 0, 0, 0, 0 };
#else
Int TEncGOP::m_iStoredAlfParaNum = 0;
#endif
#endif

/** Create list to contain pointers to CTU start addresses of slice.
 */
Void  TEncGOP::create()
{
  m_bLongtermTestPictureHasBeenCoded = 0;
  m_bLongtermTestPictureHasBeenCoded2 = 0;
}

Void  TEncGOP::destroy()
{
#if COM16_C806_ALF_TEMPPRED_NUM
  if( m_pcCfg->getUseALF() )
  {
    for( Int i = 0; i < COM16_C806_ALF_TEMPPRED_NUM; i++ )
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      for (Int j = 0; j < JVET_E0104_ALF_MAX_TEMPLAYERID; j++)
      {
        m_pcAdaptiveLoopFilter->freeALFParam(&m_acStoredAlfPara[j][i]);
      }
#else
      m_pcAdaptiveLoopFilter->freeALFParam( &m_acStoredAlfPara[i] );
#endif
    }
  }
#endif
}

Void TEncGOP::init ( TEncTop* pcTEncTop )
{
  m_pcEncTop     = pcTEncTop;
  m_pcCfg                = pcTEncTop;
  m_seiEncoder.init(m_pcCfg, pcTEncTop, this);
  m_pcSliceEncoder       = pcTEncTop->getSliceEncoder();
  m_pcListPic            = pcTEncTop->getListPic();

  m_pcEntropyCoder       = pcTEncTop->getEntropyCoder();
  m_pcCavlcCoder         = pcTEncTop->getCavlcCoder();
  m_pcSbacCoder          = pcTEncTop->getSbacCoder();
  m_pcBinCABAC           = pcTEncTop->getBinCABAC();
  m_pcLoopFilter         = pcTEncTop->getLoopFilter();

#if CNN_BASED_INLOOP_FILTER
	m_pcCNNLoopFilter      = pcTEncTop->getCNNLoopFilter();
#endif
  m_pcSAO                = pcTEncTop->getSAO();
  m_pcRateCtrl           = pcTEncTop->getRateCtrl();
  m_lastBPSEI          = 0;
  m_totalCoded         = 0;

#if ALF_HM3_REFACTOR
  if( m_pcCfg->getUseALF() )
  {
    m_pcAdaptiveLoopFilter = pcTEncTop->getAdaptiveLoopFilter();
#if COM16_C806_ALF_TEMPPRED_NUM
#if JVET_C0024_QTBT
    UInt uiMaxCUWidth = m_pcCfg->getCTUSize();
    UInt uiMaxCUHeight = m_pcCfg->getCTUSize();
#else
    UInt uiMaxCUWidth = m_pcCfg->getMaxCUWidth();
    UInt uiMaxCUHeight = m_pcCfg->getMaxCUHeight();
#endif
    UInt uiPicWidth = m_pcCfg->getSourceWidth();
    UInt uiPicHeight = m_pcCfg->getSourceHeight();
    UInt uiWidthInCU       = ( uiPicWidth %uiMaxCUWidth  ) ? uiPicWidth /uiMaxCUWidth  + 1 : uiPicWidth /uiMaxCUWidth;
    UInt uiHeightInCU      = ( uiPicHeight%uiMaxCUHeight ) ? uiPicHeight/uiMaxCUHeight + 1 : uiPicHeight/uiMaxCUHeight;
    UInt uiNumCUsInFrame   = uiWidthInCU * uiHeightInCU;
    m_pcAdaptiveLoopFilter->setNumCUsInFrame( uiNumCUsInFrame );
    for( Int i = 0; i < COM16_C806_ALF_TEMPPRED_NUM; i++ )
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      for (Int j = 0; j < JVET_E0104_ALF_MAX_TEMPLAYERID; j++)
      {
        m_pcAdaptiveLoopFilter->allocALFParam(&m_acStoredAlfPara[j][i]);
      }
#else
      m_pcAdaptiveLoopFilter->allocALFParam( &m_acStoredAlfPara[i] );
#endif
    }
#endif
  }
#endif
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    m_pcEncTop->getRdCost()->initLumaLevelToWeightTable(m_pcCfg->getIsSDR());
  }
#endif
}

Int TEncGOP::xWriteVPS (AccessUnit &accessUnit, const TComVPS *vps)
{
  OutputNALUnit nalu(NAL_UNIT_VPS);
  m_pcEntropyCoder->setBitstream(&nalu.m_Bitstream);
  m_pcEntropyCoder->encodeVPS(vps);
  accessUnit.push_back(new NALUnitEBSP(nalu));
  return (Int)(accessUnit.back()->m_nalUnitData.str().size()) * 8;
}

Int TEncGOP::xWriteSPS (AccessUnit &accessUnit, const TComSPS *sps)
{
  OutputNALUnit nalu(NAL_UNIT_SPS);
  m_pcEntropyCoder->setBitstream(&nalu.m_Bitstream);
  m_pcEntropyCoder->encodeSPS(sps);
  accessUnit.push_back(new NALUnitEBSP(nalu));
  return (Int)(accessUnit.back()->m_nalUnitData.str().size()) * 8;

}

Int TEncGOP::xWritePPS (AccessUnit &accessUnit, const TComPPS *pps)
{
  OutputNALUnit nalu(NAL_UNIT_PPS);
  m_pcEntropyCoder->setBitstream(&nalu.m_Bitstream);
  m_pcEntropyCoder->encodePPS(pps);
  accessUnit.push_back(new NALUnitEBSP(nalu));
  return (Int)(accessUnit.back()->m_nalUnitData.str().size()) * 8;
}


Int TEncGOP::xWriteParameterSets (AccessUnit &accessUnit, TComSlice *slice)
{
  Int actualTotalBits = 0;

  actualTotalBits += xWriteVPS(accessUnit, m_pcEncTop->getVPS());
  actualTotalBits += xWriteSPS(accessUnit, slice->getSPS());
  actualTotalBits += xWritePPS(accessUnit, slice->getPPS());

  return actualTotalBits;
}

// write SEI list into one NAL unit and add it to the Access unit at auPos
Void TEncGOP::xWriteSEI (NalUnitType naluType, SEIMessages& seiMessages, AccessUnit &accessUnit, AccessUnit::iterator &auPos, Int temporalId, const TComSPS *sps)
{
  // don't do anything, if we get an empty list
  if (seiMessages.empty())
  {
    return;
  }
  OutputNALUnit nalu(naluType, temporalId);
  m_seiWriter.writeSEImessages(nalu.m_Bitstream, seiMessages, sps, false);
  auPos = accessUnit.insert(auPos, new NALUnitEBSP(nalu));
  auPos++;
}

Void TEncGOP::xWriteSEISeparately (NalUnitType naluType, SEIMessages& seiMessages, AccessUnit &accessUnit, AccessUnit::iterator &auPos, Int temporalId, const TComSPS *sps)
{
  // don't do anything, if we get an empty list
  if (seiMessages.empty())
  {
    return;
  }
  for (SEIMessages::const_iterator sei = seiMessages.begin(); sei!=seiMessages.end(); sei++ )
  {
    SEIMessages tmpMessages;
    tmpMessages.push_back(*sei);
    OutputNALUnit nalu(naluType, temporalId);
    m_seiWriter.writeSEImessages(nalu.m_Bitstream, tmpMessages, sps, false);
    auPos = accessUnit.insert(auPos, new NALUnitEBSP(nalu));
    auPos++;
  }
}

Void TEncGOP::xClearSEIs(SEIMessages& seiMessages, Bool deleteMessages)
{
  if (deleteMessages)
  {
    deleteSEIs(seiMessages);
  }
  else
  {
    seiMessages.clear();
  }
}

// write SEI messages as separate NAL units ordered
Void TEncGOP::xWriteLeadingSEIOrdered (SEIMessages& seiMessages, SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, Bool testWrite)
{
  AccessUnit::iterator itNalu = accessUnit.begin();

  while ( (itNalu!=accessUnit.end())&&
    ( (*itNalu)->m_nalUnitType==NAL_UNIT_ACCESS_UNIT_DELIMITER 
    || (*itNalu)->m_nalUnitType==NAL_UNIT_VPS
    || (*itNalu)->m_nalUnitType==NAL_UNIT_SPS
    || (*itNalu)->m_nalUnitType==NAL_UNIT_PPS
    ))
  {
    itNalu++;
  }

  SEIMessages localMessages = seiMessages;
  SEIMessages currentMessages;
  
#if ENC_DEC_TRACE
  g_HLSTraceEnable = !testWrite;
#endif
  // The case that a specific SEI is not present is handled in xWriteSEI (empty list)

  // Active parameter sets SEI must always be the first SEI
  currentMessages = extractSeisByType(localMessages, SEI::ACTIVE_PARAMETER_SETS);
  assert (currentMessages.size() <= 1);
  xWriteSEI(NAL_UNIT_PREFIX_SEI, currentMessages, accessUnit, itNalu, temporalId, sps);
  xClearSEIs(currentMessages, !testWrite);
  
  // Buffering period SEI must always be following active parameter sets
  currentMessages = extractSeisByType(localMessages, SEI::BUFFERING_PERIOD);
  assert (currentMessages.size() <= 1);
  xWriteSEI(NAL_UNIT_PREFIX_SEI, currentMessages, accessUnit, itNalu, temporalId, sps);
  xClearSEIs(currentMessages, !testWrite);

  // Picture timing SEI must always be following buffering period
  currentMessages = extractSeisByType(localMessages, SEI::PICTURE_TIMING);
  assert (currentMessages.size() <= 1);
  xWriteSEI(NAL_UNIT_PREFIX_SEI, currentMessages, accessUnit, itNalu, temporalId, sps);
  xClearSEIs(currentMessages, !testWrite);

  // Decoding unit info SEI must always be following picture timing
  if (!duInfoSeiMessages.empty())
  {
    currentMessages.push_back(duInfoSeiMessages.front());
    if (!testWrite)
    {
      duInfoSeiMessages.pop_front();
    }
    xWriteSEI(NAL_UNIT_PREFIX_SEI, currentMessages, accessUnit, itNalu, temporalId, sps);
    xClearSEIs(currentMessages, !testWrite);
  }

  // Scalable nesting SEI must always be the following DU info
  currentMessages = extractSeisByType(localMessages, SEI::SCALABLE_NESTING);
  xWriteSEISeparately(NAL_UNIT_PREFIX_SEI, currentMessages, accessUnit, itNalu, temporalId, sps);
  xClearSEIs(currentMessages, !testWrite);

  // And finally everything else one by one
  xWriteSEISeparately(NAL_UNIT_PREFIX_SEI, localMessages, accessUnit, itNalu, temporalId, sps);
  xClearSEIs(localMessages, !testWrite);

  if (!testWrite)
  {
    seiMessages.clear();
  }
}


Void TEncGOP::xWriteLeadingSEIMessages (SEIMessages& seiMessages, SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, std::deque<DUData> &duData)
{
  AccessUnit testAU;
  SEIMessages picTimingSEIs = getSeisByType(seiMessages, SEI::PICTURE_TIMING);
  assert (picTimingSEIs.size() < 2);
  SEIPictureTiming * picTiming = picTimingSEIs.empty() ? NULL : (SEIPictureTiming*) picTimingSEIs.front();

  // test writing
  xWriteLeadingSEIOrdered(seiMessages, duInfoSeiMessages, testAU, temporalId, sps, true);
  // update Timing and DU info SEI
  xUpdateDuData(testAU, duData);
  xUpdateTimingSEI(picTiming, duData, sps);
  xUpdateDuInfoSEI(duInfoSeiMessages, picTiming);
  // actual writing
  xWriteLeadingSEIOrdered(seiMessages, duInfoSeiMessages, accessUnit, temporalId, sps, false);

  // testAU will automatically be cleaned up when losing scope
}

Void TEncGOP::xWriteTrailingSEIMessages (SEIMessages& seiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps)
{
  // Note: using accessUnit.end() works only as long as this function is called after slice coding and before EOS/EOB NAL units
  AccessUnit::iterator pos = accessUnit.end();
  xWriteSEISeparately(NAL_UNIT_SUFFIX_SEI, seiMessages, accessUnit, pos, temporalId, sps);
  deleteSEIs(seiMessages);
}

Void TEncGOP::xWriteDuSEIMessages (SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, std::deque<DUData> &duData)
{
  const TComHRD *hrd = sps->getVuiParameters()->getHrdParameters();

  if( m_pcCfg->getDecodingUnitInfoSEIEnabled() && hrd->getSubPicCpbParamsPresentFlag() )
  {
    Int naluIdx = 0;
    AccessUnit::iterator nalu = accessUnit.begin();

    // skip over first DU, we have a DU info SEI there already
    while (naluIdx < duData[0].accumNalsDU && nalu!=accessUnit.end())
    {
      naluIdx++;
      nalu++;
    }

    SEIMessages::iterator duSEI = duInfoSeiMessages.begin();
    // loop over remaining DUs
    for (Int duIdx = 1; duIdx < duData.size(); duIdx++)
    {
      if (duSEI == duInfoSeiMessages.end())
      {
        // if the number of generated SEIs matches the number of DUs, this should not happen
        assert (false);
        return;
      }
      // write the next SEI
      SEIMessages tmpSEI;
      tmpSEI.push_back(*duSEI);
      xWriteSEI(NAL_UNIT_PREFIX_SEI, tmpSEI, accessUnit, nalu, temporalId, sps);
      // nalu points to the position after the SEI, so we have to increase the index as well
      naluIdx++;
      while ((naluIdx < duData[duIdx].accumNalsDU) && nalu!=accessUnit.end())
      {
        naluIdx++;
        nalu++;
      }
      duSEI++;
    }
  }
  deleteSEIs(duInfoSeiMessages);
}


Void TEncGOP::xCreateIRAPLeadingSEIMessages (SEIMessages& seiMessages, const TComSPS *sps, const TComPPS *pps)
{
  OutputNALUnit nalu(NAL_UNIT_PREFIX_SEI);

  if(m_pcCfg->getActiveParameterSetsSEIEnabled())
  {
    SEIActiveParameterSets *sei = new SEIActiveParameterSets;
    m_seiEncoder.initSEIActiveParameterSets (sei, m_pcCfg->getVPS(), sps);
    seiMessages.push_back(sei);
  }

  if(m_pcCfg->getFramePackingArrangementSEIEnabled())
  {
    SEIFramePacking *sei = new SEIFramePacking;
    m_seiEncoder.initSEIFramePacking (sei, m_iNumPicCoded);
    seiMessages.push_back(sei);
  }

  if(m_pcCfg->getSegmentedRectFramePackingArrangementSEIEnabled())
  {
    SEISegmentedRectFramePacking *sei = new SEISegmentedRectFramePacking;
    m_seiEncoder.initSEISegmentedRectFramePacking(sei);
    seiMessages.push_back(sei);
  }

  if (m_pcCfg->getDisplayOrientationSEIAngle())
  {
    SEIDisplayOrientation *sei = new SEIDisplayOrientation;
    m_seiEncoder.initSEIDisplayOrientation(sei);
    seiMessages.push_back(sei);
  }

  if(m_pcCfg->getToneMappingInfoSEIEnabled())
  {
    SEIToneMappingInfo *sei = new SEIToneMappingInfo;
    m_seiEncoder.initSEIToneMappingInfo (sei);
    seiMessages.push_back(sei);
  }

  if(m_pcCfg->getTMCTSSEIEnabled())
  {
    SEITempMotionConstrainedTileSets *sei = new SEITempMotionConstrainedTileSets;
    m_seiEncoder.initSEITempMotionConstrainedTileSets(sei, pps);
    seiMessages.push_back(sei);
  }

  if(m_pcCfg->getTimeCodeSEIEnabled())
  {
    SEITimeCode *seiTimeCode = new SEITimeCode;
    m_seiEncoder.initSEITimeCode(seiTimeCode);
    seiMessages.push_back(seiTimeCode);
  }

  if(m_pcCfg->getKneeSEIEnabled())
  {
    SEIKneeFunctionInfo *sei = new SEIKneeFunctionInfo;
    m_seiEncoder.initSEIKneeFunctionInfo(sei);
    seiMessages.push_back(sei);
  }
    
  if(m_pcCfg->getMasteringDisplaySEI().colourVolumeSEIEnabled)
  {
    const TComSEIMasteringDisplay &seiCfg=m_pcCfg->getMasteringDisplaySEI();
    SEIMasteringDisplayColourVolume *sei = new SEIMasteringDisplayColourVolume;
    sei->values = seiCfg;
    seiMessages.push_back(sei);
  }
}

Void TEncGOP::xCreatePerPictureSEIMessages (Int picInGOP, SEIMessages& seiMessages, SEIMessages& nestedSeiMessages, TComSlice *slice)
{
  if( ( m_pcCfg->getBufferingPeriodSEIEnabled() ) && ( slice->getSliceType() == I_SLICE ) &&
    ( slice->getSPS()->getVuiParametersPresentFlag() ) &&
    ( ( slice->getSPS()->getVuiParameters()->getHrdParameters()->getNalHrdParametersPresentFlag() )
    || ( slice->getSPS()->getVuiParameters()->getHrdParameters()->getVclHrdParametersPresentFlag() ) ) )
  {
    SEIBufferingPeriod *bufferingPeriodSEI = new SEIBufferingPeriod();
    m_seiEncoder.initSEIBufferingPeriod(bufferingPeriodSEI, slice);
    seiMessages.push_back(bufferingPeriodSEI);
    m_bufferingPeriodSEIPresentInAU = true;

    if (m_pcCfg->getScalableNestingSEIEnabled())
    {
      SEIBufferingPeriod *bufferingPeriodSEIcopy = new SEIBufferingPeriod();
      bufferingPeriodSEI->copyTo(*bufferingPeriodSEIcopy);
      nestedSeiMessages.push_back(bufferingPeriodSEIcopy);
    }
  }

  if (picInGOP ==0 && m_pcCfg->getSOPDescriptionSEIEnabled() ) // write SOP description SEI (if enabled) at the beginning of GOP
  {
    SEISOPDescription* sopDescriptionSEI = new SEISOPDescription();
    m_seiEncoder.initSEISOPDescription(sopDescriptionSEI, slice, picInGOP, m_iLastIDR, m_iGopSize);
    seiMessages.push_back(sopDescriptionSEI);
  }

  if( ( m_pcEncTop->getRecoveryPointSEIEnabled() ) && ( slice->getSliceType() == I_SLICE ) )
  {
    if( m_pcEncTop->getGradualDecodingRefreshInfoEnabled() && !slice->getRapPicFlag() )
    {
      // Gradual decoding refresh SEI
      SEIGradualDecodingRefreshInfo *gradualDecodingRefreshInfoSEI = new SEIGradualDecodingRefreshInfo();
      gradualDecodingRefreshInfoSEI->m_gdrForegroundFlag = true; // Indicating all "foreground"
      seiMessages.push_back(gradualDecodingRefreshInfoSEI);
    }
    // Recovery point SEI
    SEIRecoveryPoint *recoveryPointSEI = new SEIRecoveryPoint();
    m_seiEncoder.initSEIRecoveryPoint(recoveryPointSEI, slice);
    seiMessages.push_back(recoveryPointSEI);
  }
  if (m_pcCfg->getTemporalLevel0IndexSEIEnabled())
  {
    SEITemporalLevel0Index *temporalLevel0IndexSEI = new SEITemporalLevel0Index();
    m_seiEncoder.initTemporalLevel0IndexSEI(temporalLevel0IndexSEI, slice);
    seiMessages.push_back(temporalLevel0IndexSEI);
  }

  if(slice->getSPS()->getVuiParametersPresentFlag() && m_pcCfg->getChromaSamplingFilterHintEnabled() && ( slice->getSliceType() == I_SLICE ))
  {
    SEIChromaSamplingFilterHint *seiChromaSamplingFilterHint = new SEIChromaSamplingFilterHint;
    m_seiEncoder.initSEIChromaSamplingFilterHint(seiChromaSamplingFilterHint, m_pcCfg->getChromaSamplingHorFilterIdc(), m_pcCfg->getChromaSamplingVerFilterIdc());
    seiMessages.push_back(seiChromaSamplingFilterHint);
  }

  if( m_pcEncTop->getNoDisplaySEITLayer() && ( slice->getTLayer() >= m_pcEncTop->getNoDisplaySEITLayer() ) )
  {
    SEINoDisplay *seiNoDisplay = new SEINoDisplay;
    seiNoDisplay->m_noDisplay = true;
    seiMessages.push_back(seiNoDisplay);
  }
}

Void TEncGOP::xCreateScalableNestingSEI (SEIMessages& seiMessages, SEIMessages& nestedSeiMessages)
{
  SEIMessages tmpMessages;
  while (!nestedSeiMessages.empty())
  {
    SEI* sei=nestedSeiMessages.front();
    nestedSeiMessages.pop_front();
    tmpMessages.push_back(sei);
    SEIScalableNesting *nestingSEI = new SEIScalableNesting();
    m_seiEncoder.initSEIScalableNesting(nestingSEI, tmpMessages);
    seiMessages.push_back(nestingSEI);
    tmpMessages.clear();
  }
}

Void TEncGOP::xCreatePictureTimingSEI  (Int IRAPGOPid, SEIMessages& seiMessages, SEIMessages& nestedSeiMessages, SEIMessages& duInfoSeiMessages, TComSlice *slice, Bool isField, std::deque<DUData> &duData)
{
  Int picSptDpbOutputDuDelay = 0;
  SEIPictureTiming *pictureTimingSEI = new SEIPictureTiming();

  const TComVUI *vui = slice->getSPS()->getVuiParameters();
  const TComHRD *hrd = vui->getHrdParameters();

  // update decoding unit parameters
  if( ( m_pcCfg->getPictureTimingSEIEnabled() || m_pcCfg->getDecodingUnitInfoSEIEnabled() ) &&
    ( slice->getSPS()->getVuiParametersPresentFlag() ) &&
    (  hrd->getNalHrdParametersPresentFlag() || hrd->getVclHrdParametersPresentFlag() ) )
  {
    // DU parameters
    if( hrd->getSubPicCpbParamsPresentFlag() )
    {
      UInt numDU = (UInt) duData.size();
      pictureTimingSEI->m_numDecodingUnitsMinus1     = ( numDU - 1 );
      pictureTimingSEI->m_duCommonCpbRemovalDelayFlag = false;
      pictureTimingSEI->m_numNalusInDuMinus1.resize( numDU );
      pictureTimingSEI->m_duCpbRemovalDelayMinus1.resize( numDU );
    }
    pictureTimingSEI->m_auCpbRemovalDelay = std::min<Int>(std::max<Int>(1, m_totalCoded - m_lastBPSEI), static_cast<Int>(pow(2, static_cast<Double>(hrd->getCpbRemovalDelayLengthMinus1()+1)))); // Syntax element signalled as minus, hence the .
    pictureTimingSEI->m_picDpbOutputDelay = slice->getSPS()->getNumReorderPics(slice->getSPS()->getMaxTLayers()-1) + slice->getPOC() - m_totalCoded;
    if(m_pcCfg->getEfficientFieldIRAPEnabled() && IRAPGOPid > 0 && IRAPGOPid < m_iGopSize)
    {
      // if pictures have been swapped there is likely one more picture delay on their tid. Very rough approximation
      pictureTimingSEI->m_picDpbOutputDelay ++;
    }
    Int factor = hrd->getTickDivisorMinus2() + 2;
    pictureTimingSEI->m_picDpbOutputDuDelay = factor * pictureTimingSEI->m_picDpbOutputDelay;
    if( m_pcCfg->getDecodingUnitInfoSEIEnabled() )
    {
      picSptDpbOutputDuDelay = factor * pictureTimingSEI->m_picDpbOutputDelay;
    }
    if (m_bufferingPeriodSEIPresentInAU)
    {
      m_lastBPSEI = m_totalCoded;
    }

    if( hrd->getSubPicCpbParamsPresentFlag() )
    {
      Int i;
      UInt64 ui64Tmp;
      UInt uiPrev = 0;
      UInt numDU = ( pictureTimingSEI->m_numDecodingUnitsMinus1 + 1 );
      std::vector<UInt> &rDuCpbRemovalDelayMinus1 = pictureTimingSEI->m_duCpbRemovalDelayMinus1;
      UInt maxDiff = ( hrd->getTickDivisorMinus2() + 2 ) - 1;

      for( i = 0; i < numDU; i ++ )
      {
        pictureTimingSEI->m_numNalusInDuMinus1[ i ]       = ( i == 0 ) ? ( duData[i].accumNalsDU - 1 ) : ( duData[i].accumNalsDU- duData[i-1].accumNalsDU - 1 );
      }

      if( numDU == 1 )
      {
        rDuCpbRemovalDelayMinus1[ 0 ] = 0; /* don't care */
      }
      else
      {
        rDuCpbRemovalDelayMinus1[ numDU - 1 ] = 0;/* by definition */
        UInt tmp = 0;
        UInt accum = 0;

        for( i = ( numDU - 2 ); i >= 0; i -- )
        {
          ui64Tmp = ( ( ( duData[numDU - 1].accumBitsDU  - duData[i].accumBitsDU ) * ( vui->getTimingInfo()->getTimeScale() / vui->getTimingInfo()->getNumUnitsInTick() ) * ( hrd->getTickDivisorMinus2() + 2 ) ) / ( m_pcCfg->getTargetBitrate() ) );
          if( (UInt)ui64Tmp > maxDiff )
          {
            tmp ++;
          }
        }
        uiPrev = 0;

        UInt flag = 0;
        for( i = ( numDU - 2 ); i >= 0; i -- )
        {
          flag = 0;
          ui64Tmp = ( ( ( duData[numDU - 1].accumBitsDU  - duData[i].accumBitsDU ) * ( vui->getTimingInfo()->getTimeScale() / vui->getTimingInfo()->getNumUnitsInTick() ) * ( hrd->getTickDivisorMinus2() + 2 ) ) / ( m_pcCfg->getTargetBitrate() ) );

          if( (UInt)ui64Tmp > maxDiff )
          {
            if(uiPrev >= maxDiff - tmp)
            {
              ui64Tmp = uiPrev + 1;
              flag = 1;
            }
            else                            ui64Tmp = maxDiff - tmp + 1;
          }
          rDuCpbRemovalDelayMinus1[ i ] = (UInt)ui64Tmp - uiPrev - 1;
          if( (Int)rDuCpbRemovalDelayMinus1[ i ] < 0 )
          {
            rDuCpbRemovalDelayMinus1[ i ] = 0;
          }
          else if (tmp > 0 && flag == 1)
          {
            tmp --;
          }
          accum += rDuCpbRemovalDelayMinus1[ i ] + 1;
          uiPrev = accum;
        }
      }
    }
    
    if( m_pcCfg->getPictureTimingSEIEnabled() )
    {
      pictureTimingSEI->m_picStruct = (isField && slice->getPic()->isTopField())? 1 : isField? 2 : 0;
      seiMessages.push_back(pictureTimingSEI);

      if ( m_pcCfg->getScalableNestingSEIEnabled() ) // put picture timing SEI into scalable nesting SEI
      {
        SEIPictureTiming *pictureTimingSEIcopy = new SEIPictureTiming();
        pictureTimingSEI->copyTo(*pictureTimingSEIcopy);
        nestedSeiMessages.push_back(pictureTimingSEIcopy);
      }
    }

    if( m_pcCfg->getDecodingUnitInfoSEIEnabled() && hrd->getSubPicCpbParamsPresentFlag() )
    {
      for( Int i = 0; i < ( pictureTimingSEI->m_numDecodingUnitsMinus1 + 1 ); i ++ )
      {
        SEIDecodingUnitInfo *duInfoSEI = new SEIDecodingUnitInfo();
        duInfoSEI->m_decodingUnitIdx = i;
        duInfoSEI->m_duSptCpbRemovalDelay = pictureTimingSEI->m_duCpbRemovalDelayMinus1[i] + 1;
        duInfoSEI->m_dpbOutputDuDelayPresentFlag = false;
        duInfoSEI->m_picSptDpbOutputDuDelay = picSptDpbOutputDuDelay;

        duInfoSeiMessages.push_back(duInfoSEI);
      }
    }
  }
}

Void TEncGOP::xUpdateDuData(AccessUnit &testAU, std::deque<DUData> &duData)
{
  if (duData.empty())
  {
    return;
  }
  // fix first 
  UInt numNalUnits = (UInt)testAU.size();
  UInt numRBSPBytes = 0;
  for (AccessUnit::const_iterator it = testAU.begin(); it != testAU.end(); it++)
  {
    numRBSPBytes += UInt((*it)->m_nalUnitData.str().size());
  }
  duData[0].accumBitsDU += ( numRBSPBytes << 3 );
  duData[0].accumNalsDU += numNalUnits;

  // adapt cumulative sums for all following DUs
  // and add one DU info SEI, if enabled
  for (Int i=1; i<duData.size(); i++)
  {
    if (m_pcCfg->getDecodingUnitInfoSEIEnabled())
    {
      numNalUnits  += 1;
      numRBSPBytes += ( 5 << 3 );
    }
    duData[i].accumBitsDU += numRBSPBytes; // probably around 5 bytes
    duData[i].accumNalsDU += numNalUnits;
  }

  // The last DU may have a trailing SEI
  if (m_pcCfg->getDecodedPictureHashSEIEnabled())
  {
    duData.back().accumBitsDU += ( 20 << 3 ); // probably around 20 bytes - should be further adjusted, e.g. by type
    duData.back().accumNalsDU += 1;
  }

}
Void TEncGOP::xUpdateTimingSEI(SEIPictureTiming *pictureTimingSEI, std::deque<DUData> &duData, const TComSPS *sps)
{
  if (!pictureTimingSEI)
  {
    return;
  }
  const TComVUI *vui = sps->getVuiParameters();
  const TComHRD *hrd = vui->getHrdParameters();
  if( hrd->getSubPicCpbParamsPresentFlag() )
  {
    Int i;
    UInt64 ui64Tmp;
    UInt uiPrev = 0;
    UInt numDU = ( pictureTimingSEI->m_numDecodingUnitsMinus1 + 1 );
    std::vector<UInt> &rDuCpbRemovalDelayMinus1 = pictureTimingSEI->m_duCpbRemovalDelayMinus1;
    UInt maxDiff = ( hrd->getTickDivisorMinus2() + 2 ) - 1;

    for( i = 0; i < numDU; i ++ )
    {
      pictureTimingSEI->m_numNalusInDuMinus1[ i ]       = ( i == 0 ) ? ( duData[i].accumNalsDU - 1 ) : ( duData[i].accumNalsDU- duData[i-1].accumNalsDU - 1 );
    }

    if( numDU == 1 )
    {
      rDuCpbRemovalDelayMinus1[ 0 ] = 0; /* don't care */
    }
    else
    {
      rDuCpbRemovalDelayMinus1[ numDU - 1 ] = 0;/* by definition */
      UInt tmp = 0;
      UInt accum = 0;

      for( i = ( numDU - 2 ); i >= 0; i -- )
      {
        ui64Tmp = ( ( ( duData[numDU - 1].accumBitsDU  - duData[i].accumBitsDU ) * ( vui->getTimingInfo()->getTimeScale() / vui->getTimingInfo()->getNumUnitsInTick() ) * ( hrd->getTickDivisorMinus2() + 2 ) ) / ( m_pcCfg->getTargetBitrate() ) );
        if( (UInt)ui64Tmp > maxDiff )
        {
          tmp ++;
        }
      }
      uiPrev = 0;

      UInt flag = 0;
      for( i = ( numDU - 2 ); i >= 0; i -- )
      {
        flag = 0;
        ui64Tmp = ( ( ( duData[numDU - 1].accumBitsDU  - duData[i].accumBitsDU ) * ( vui->getTimingInfo()->getTimeScale() / vui->getTimingInfo()->getNumUnitsInTick() ) * ( hrd->getTickDivisorMinus2() + 2 ) ) / ( m_pcCfg->getTargetBitrate() ) );

        if( (UInt)ui64Tmp > maxDiff )
        {
          if(uiPrev >= maxDiff - tmp)
          {
            ui64Tmp = uiPrev + 1;
            flag = 1;
          }
          else                            ui64Tmp = maxDiff - tmp + 1;
        }
        rDuCpbRemovalDelayMinus1[ i ] = (UInt)ui64Tmp - uiPrev - 1;
        if( (Int)rDuCpbRemovalDelayMinus1[ i ] < 0 )
        {
          rDuCpbRemovalDelayMinus1[ i ] = 0;
        }
        else if (tmp > 0 && flag == 1)
        {
          tmp --;
        }
        accum += rDuCpbRemovalDelayMinus1[ i ] + 1;
        uiPrev = accum;
      }
    }
  }
}
Void TEncGOP::xUpdateDuInfoSEI(SEIMessages &duInfoSeiMessages, SEIPictureTiming *pictureTimingSEI)
{
  if (duInfoSeiMessages.empty() || (pictureTimingSEI == NULL))
  {
    return;
  }

  Int i=0;

  for (SEIMessages::iterator du = duInfoSeiMessages.begin(); du!= duInfoSeiMessages.end(); du++)
  {
    SEIDecodingUnitInfo *duInfoSEI = (SEIDecodingUnitInfo*) (*du);
    duInfoSEI->m_decodingUnitIdx = i;
    duInfoSEI->m_duSptCpbRemovalDelay = pictureTimingSEI->m_duCpbRemovalDelayMinus1[i] + 1;
    duInfoSEI->m_dpbOutputDuDelayPresentFlag = false;
    i++;
  }
}

static Void
cabac_zero_word_padding(TComSlice *const pcSlice, TComPic *const pcPic, const std::size_t binCountsInNalUnits, const std::size_t numBytesInVclNalUnits, std::ostringstream &nalUnitData, const Bool cabacZeroWordPaddingEnabled)
{
  const TComSPS &sps=*(pcSlice->getSPS());
  const Int log2subWidthCxsubHeightC = (pcPic->getComponentScaleX(COMPONENT_Cb)+pcPic->getComponentScaleY(COMPONENT_Cb));
  const Int minCuWidth  = pcPic->getMinCUWidth();
  const Int minCuHeight = pcPic->getMinCUHeight();
  const Int paddedWidth = ((sps.getPicWidthInLumaSamples()  + minCuWidth  - 1) / minCuWidth) * minCuWidth;
  const Int paddedHeight= ((sps.getPicHeightInLumaSamples() + minCuHeight - 1) / minCuHeight) * minCuHeight;
  const Int rawBits = paddedWidth * paddedHeight *
                         (sps.getBitDepth(CHANNEL_TYPE_LUMA) + 2*(sps.getBitDepth(CHANNEL_TYPE_CHROMA)>>log2subWidthCxsubHeightC));
  const std::size_t threshold = (32/3)*numBytesInVclNalUnits + (rawBits/32);
  if (binCountsInNalUnits >= threshold)
  {
    // need to add additional cabac zero words (each one accounts for 3 bytes (=00 00 03)) to increase numBytesInVclNalUnits
    const std::size_t targetNumBytesInVclNalUnits = ((binCountsInNalUnits - (rawBits/32))*3+31)/32;

    if (targetNumBytesInVclNalUnits>numBytesInVclNalUnits) // It should be!
    {
      const std::size_t numberOfAdditionalBytesNeeded=targetNumBytesInVclNalUnits - numBytesInVclNalUnits;
      const std::size_t numberOfAdditionalCabacZeroWords=(numberOfAdditionalBytesNeeded+2)/3;
      const std::size_t numberOfAdditionalCabacZeroBytes=numberOfAdditionalCabacZeroWords*3;
      if (cabacZeroWordPaddingEnabled)
      {
        std::vector<Char> zeroBytesPadding(numberOfAdditionalCabacZeroBytes, Char(0));
        for(std::size_t i=0; i<numberOfAdditionalCabacZeroWords; i++)
        {
          zeroBytesPadding[i*3+2]=3;  // 00 00 03
        }
        nalUnitData.write(&(zeroBytesPadding[0]), numberOfAdditionalCabacZeroBytes);
        printf("Adding %d bytes of padding\n", UInt(numberOfAdditionalCabacZeroWords*3));
      }
      else
      {
        printf("Standard would normally require adding %d bytes of padding\n", UInt(numberOfAdditionalCabacZeroWords*3));
      }
    }
  }
}

class EfficientFieldIRAPMapping
{
  private:
    Int  IRAPGOPid;
    Bool IRAPtoReorder;
    Bool swapIRAPForward;

  public:
    EfficientFieldIRAPMapping() :
      IRAPGOPid(-1),
      IRAPtoReorder(false),
      swapIRAPForward(false)
    { }

    Void initialize(const Bool isField, const Int gopSize, const Int POCLast, const Int numPicRcvd, const Int lastIDR, TEncGOP *pEncGop, TEncCfg *pCfg);

    Int adjustGOPid(const Int gopID);
    Int restoreGOPid(const Int gopID);
    Int GetIRAPGOPid() const { return IRAPGOPid; }
};

Void EfficientFieldIRAPMapping::initialize(const Bool isField, const Int gopSize, const Int POCLast, const Int numPicRcvd, const Int lastIDR, TEncGOP *pEncGop, TEncCfg *pCfg )
{
  if(isField)
  {
    Int pocCurr;
    for ( Int iGOPid=0; iGOPid < gopSize; iGOPid++ )
    {
      // determine actual POC
      if(POCLast == 0) //case first frame or first top field
      {
        pocCurr=0;
      }
      else if(POCLast == 1 && isField) //case first bottom field, just like the first frame, the poc computation is not right anymore, we set the right value
      {
        pocCurr = 1;
      }
      else
      {
        pocCurr = POCLast - numPicRcvd + pCfg->getGOPEntry(iGOPid).m_POC - isField;
      }

      // check if POC corresponds to IRAP
      NalUnitType tmpUnitType = pEncGop->getNalUnitType(pocCurr, lastIDR, isField);
      if(tmpUnitType >= NAL_UNIT_CODED_SLICE_BLA_W_LP && tmpUnitType <= NAL_UNIT_CODED_SLICE_CRA) // if picture is an IRAP
      {
        if(pocCurr%2 == 0 && iGOPid < gopSize-1 && pCfg->getGOPEntry(iGOPid).m_POC == pCfg->getGOPEntry(iGOPid+1).m_POC-1)
        { // if top field and following picture in enc order is associated bottom field
          IRAPGOPid = iGOPid;
          IRAPtoReorder = true;
          swapIRAPForward = true; 
          break;
        }
        if(pocCurr%2 != 0 && iGOPid > 0 && pCfg->getGOPEntry(iGOPid).m_POC == pCfg->getGOPEntry(iGOPid-1).m_POC+1)
        {
          // if picture is an IRAP remember to process it first
          IRAPGOPid = iGOPid;
          IRAPtoReorder = true;
          swapIRAPForward = false; 
          break;
        }
      }
    }
  }
}

Int EfficientFieldIRAPMapping::adjustGOPid(const Int GOPid)
{
  if(IRAPtoReorder)
  {
    if(swapIRAPForward)
    {
      if(GOPid == IRAPGOPid)
      {
        return IRAPGOPid +1;
      }
      else if(GOPid == IRAPGOPid +1)
      {
        return IRAPGOPid;
      }
    }
    else
    {
      if(GOPid == IRAPGOPid -1)
      {
        return IRAPGOPid;
      }
      else if(GOPid == IRAPGOPid)
      {
        return IRAPGOPid -1;
      }
    }
  }
  return GOPid;
}

Int EfficientFieldIRAPMapping::restoreGOPid(const Int GOPid)
{
  if(IRAPtoReorder)
  {
    if(swapIRAPForward)
    {
      if(GOPid == IRAPGOPid)
      {
        IRAPtoReorder = false;
        return IRAPGOPid +1;
      }
      else if(GOPid == IRAPGOPid +1)
      {
        return GOPid -1;
      }
    }
    else
    {
      if(GOPid == IRAPGOPid)
      {
        return IRAPGOPid -1;
      }
      else if(GOPid == IRAPGOPid -1)
      {
        IRAPtoReorder = false;
        return IRAPGOPid;
      }
    }
  }
  return GOPid;
}

#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
static UInt calculateCollocatedFromL0Flag(TComSlice *pSlice)
{
  Int refIdx = 0; // Zero always assumed
  TComPic *refPicL0 = pSlice->getRefPic(REF_PIC_LIST_0, refIdx);
  TComPic *refPicL1 = pSlice->getRefPic(REF_PIC_LIST_1, refIdx);
  return refPicL0->getSlice(0)->getSliceQp() > refPicL1->getSlice(0)->getSliceQp();
}
#else
static UInt calculateCollocatedFromL1Flag(TEncCfg *pCfg, const Int GOPid, const Int gopSize)
{
  Int iCloseLeft=1, iCloseRight=-1;
  for(Int i = 0; i<pCfg->getGOPEntry(GOPid).m_numRefPics; i++)
  {
    Int iRef = pCfg->getGOPEntry(GOPid).m_referencePics[i];
    if(iRef>0&&(iRef<iCloseRight||iCloseRight==-1))
    {
      iCloseRight=iRef;
    }
    else if(iRef<0&&(iRef>iCloseLeft||iCloseLeft==1))
    {
      iCloseLeft=iRef;
    }
  }
  if(iCloseRight>-1)
  {
    iCloseRight=iCloseRight+pCfg->getGOPEntry(GOPid).m_POC-1;
  }
  if(iCloseLeft<1)
  {
    iCloseLeft=iCloseLeft+pCfg->getGOPEntry(GOPid).m_POC-1;
    while(iCloseLeft<0)
    {
      iCloseLeft+=gopSize;
    }
  }
  Int iLeftQP=0, iRightQP=0;
  for(Int i=0; i<gopSize; i++)
  {
    if(pCfg->getGOPEntry(i).m_POC==(iCloseLeft%gopSize)+1)
    {
      iLeftQP= pCfg->getGOPEntry(i).m_QPOffset;
    }
    if (pCfg->getGOPEntry(i).m_POC==(iCloseRight%gopSize)+1)
    {
      iRightQP=pCfg->getGOPEntry(i).m_QPOffset;
    }
  }
  if(iCloseRight>-1&&iRightQP<iLeftQP)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
#endif

// ====================================================================================================================
// Public member functions
// ====================================================================================================================
Void TEncGOP::compressGOP( Int iPOCLast, Int iNumPicRcvd, TComList<TComPic*>& rcListPic,
                           TComList<TComPicYuv*>& rcListPicYuvRecOut, std::list<AccessUnit>& accessUnitsInGOP,
                           Bool isField, Bool isTff, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE
#if JVET_F0064_MSSSIM
                           , const Bool printMSSSIM
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
                         , TComStats* m_apcStats
#endif
                           )
{
  // TODO: Split this function up.

  TComPic*        pcPic = NULL;
  TComPicYuv*     pcPicYuvRecOut;
  TComSlice*      pcSlice;
  TComOutputBitstream  *pcBitstreamRedirect;
  pcBitstreamRedirect = new TComOutputBitstream;
  AccessUnit::iterator  itLocationToPushSliceHeaderNALU; // used to store location where NALU containing slice header is to be inserted

  xInitGOP( iPOCLast, iNumPicRcvd, isField );

#if ALF_HM3_REFACTOR
  ALFParam cAlfParam;
  Bool bInitAlfParam = false;
#endif
  m_iNumPicCoded = 0;
  SEIMessages leadingSeiMessages;
  SEIMessages nestedSeiMessages;
  SEIMessages duInfoSeiMessages;
  SEIMessages trailingSeiMessages;
  std::deque<DUData> duData;
  SEIDecodingUnitInfo decodingUnitInfoSEI;

  EfficientFieldIRAPMapping effFieldIRAPMap;
  if (m_pcCfg->getEfficientFieldIRAPEnabled())
  {
    effFieldIRAPMap.initialize(isField, m_iGopSize, iPOCLast, iNumPicRcvd, m_iLastIDR, this, m_pcCfg);
  }

  // reset flag indicating whether pictures have been encoded
  for ( Int iGOPid=0; iGOPid < m_iGopSize; iGOPid++ )
  {
    m_pcCfg->setEncodedFlag(iGOPid, false);
  }

  for ( Int iGOPid=0; iGOPid < m_iGopSize; iGOPid++ )
  {
    if (m_pcCfg->getEfficientFieldIRAPEnabled())
    {
      iGOPid=effFieldIRAPMap.adjustGOPid(iGOPid);
    }

    //-- For time output for each slice
    clock_t iBeforeTime = clock();

#if !JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    UInt uiColDir = calculateCollocatedFromL1Flag(m_pcCfg, iGOPid, m_iGopSize);
#endif
    /////////////////////////////////////////////////////////////////////////////////////////////////// Initial to start encoding
    Int iTimeOffset;
    Int pocCurr;

    if(iPOCLast == 0) //case first frame or first top field
    {
      pocCurr=0;
      iTimeOffset = 1;
    }
    else if(iPOCLast == 1 && isField) //case first bottom field, just like the first frame, the poc computation is not right anymore, we set the right value
    {
      pocCurr = 1;
      iTimeOffset = 1;
    }
    else
    {
      pocCurr = iPOCLast - iNumPicRcvd + m_pcCfg->getGOPEntry(iGOPid).m_POC - ((isField && m_iGopSize>1) ? 1:0);
      iTimeOffset = m_pcCfg->getGOPEntry(iGOPid).m_POC;
    }

    if(pocCurr>=m_pcCfg->getFramesToBeEncoded())
    {
      if (m_pcCfg->getEfficientFieldIRAPEnabled())
      {
        iGOPid=effFieldIRAPMap.restoreGOPid(iGOPid);
      }
      continue;
    }

    if( getNalUnitType(pocCurr, m_iLastIDR, isField) == NAL_UNIT_CODED_SLICE_IDR_W_RADL || getNalUnitType(pocCurr, m_iLastIDR, isField) == NAL_UNIT_CODED_SLICE_IDR_N_LP )
    {
      m_iLastIDR = pocCurr;
    }
    // start a new access unit: create an entry in the list of output access units
    accessUnitsInGOP.push_back(AccessUnit());
    AccessUnit& accessUnit = accessUnitsInGOP.back();
    xGetBuffer( rcListPic, rcListPicYuvRecOut, iNumPicRcvd, iTimeOffset, pcPic, pcPicYuvRecOut, pocCurr, isField );

    //  Slice data initialization
    pcPic->clearSliceBuffer();
    pcPic->allocateNewSlice();
    m_pcSliceEncoder->setSliceIdx(0);
    pcPic->setCurrSliceIdx(0);

    m_pcSliceEncoder->initEncSlice ( pcPic, iPOCLast, pocCurr, iGOPid, pcSlice, isField );

    //Set Frame/Field coding
    pcSlice->getPic()->setField(isField);

    pcSlice->setLastIDR(m_iLastIDR);
    pcSlice->setSliceIdx(0);
    //set default slice level flag to the same as SPS level flag
    pcSlice->setLFCrossSliceBoundaryFlag(  pcSlice->getPPS()->getLoopFilterAcrossSlicesEnabledFlag()  );

    if(pcSlice->getSliceType()==B_SLICE&&m_pcCfg->getGOPEntry(iGOPid).m_sliceType=='P')
    {
      pcSlice->setSliceType(P_SLICE);
    }
    if(pcSlice->getSliceType()==B_SLICE&&m_pcCfg->getGOPEntry(iGOPid).m_sliceType=='I')
    {
      pcSlice->setSliceType(I_SLICE);
    }
    
    // Set the nal unit type
    pcSlice->setNalUnitType(getNalUnitType(pocCurr, m_iLastIDR, isField));
    if(pcSlice->getTemporalLayerNonReferenceFlag())
    {
      if (pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_TRAIL_R &&
          !(m_iGopSize == 1 && pcSlice->getSliceType() == I_SLICE))
        // Add this condition to avoid POC issues with encoder_intra_main.cfg configuration (see #1127 in bug tracker)
      {
        pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_TRAIL_N);
      }
      if(pcSlice->getNalUnitType()==NAL_UNIT_CODED_SLICE_RADL_R)
      {
        pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_RADL_N);
      }
      if(pcSlice->getNalUnitType()==NAL_UNIT_CODED_SLICE_RASL_R)
      {
        pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_RASL_N);
      }
    }

    if (m_pcCfg->getEfficientFieldIRAPEnabled())
    {
      if ( pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_W_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_W_RADL
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_N_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_IDR_W_RADL
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_IDR_N_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_CRA )  // IRAP picture
      {
        m_associatedIRAPType = pcSlice->getNalUnitType();
        m_associatedIRAPPOC = pocCurr;
      }
      pcSlice->setAssociatedIRAPType(m_associatedIRAPType);
      pcSlice->setAssociatedIRAPPOC(m_associatedIRAPPOC);
    }
    // Do decoding refresh marking if any
#if FIX_TICKET12 
    if( pcSlice->getSPS()->getUseALF() )
    {
      if( m_pcAdaptiveLoopFilter->refreshAlfTempPred( pcSlice->getNalUnitType() , pcSlice->getPOC() ) )
      {
#if JVET_E0104_ALF_TEMP_SCALABILITY
        memset(m_iStoredAlfParaNum, 0, sizeof(Int)*JVET_E0104_ALF_MAX_TEMPLAYERID);
#else
        m_iStoredAlfParaNum = 0;
#endif
      }
    }
#endif
    pcSlice->decodingRefreshMarking(m_pocCRA, m_bRefreshPending, rcListPic, m_pcCfg->getEfficientFieldIRAPEnabled());
    m_pcEncTop->selectReferencePictureSet(pcSlice, pocCurr, iGOPid);
    if (!m_pcCfg->getEfficientFieldIRAPEnabled())
    {
      if ( pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_W_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_W_RADL
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_BLA_N_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_IDR_W_RADL
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_IDR_N_LP
        || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_CRA )  // IRAP picture
      {
        m_associatedIRAPType = pcSlice->getNalUnitType();
        m_associatedIRAPPOC = pocCurr;
      }
      pcSlice->setAssociatedIRAPType(m_associatedIRAPType);
      pcSlice->setAssociatedIRAPPOC(m_associatedIRAPPOC);
    }

    if ((pcSlice->checkThatAllRefPicsAreAvailable(rcListPic, pcSlice->getRPS(), false, m_iLastRecoveryPicPOC, m_pcCfg->getDecodingRefreshType() == 3) != 0) || (pcSlice->isIRAP()) 
      || (m_pcCfg->getEfficientFieldIRAPEnabled() && isField && pcSlice->getAssociatedIRAPType() >= NAL_UNIT_CODED_SLICE_BLA_W_LP && pcSlice->getAssociatedIRAPType() <= NAL_UNIT_CODED_SLICE_CRA && pcSlice->getAssociatedIRAPPOC() == pcSlice->getPOC()+1)
      )
    {
      pcSlice->createExplicitReferencePictureSetFromReference(rcListPic, pcSlice->getRPS(), pcSlice->isIRAP(), m_iLastRecoveryPicPOC, m_pcCfg->getDecodingRefreshType() == 3, m_pcCfg->getEfficientFieldIRAPEnabled());
    }

    pcSlice->applyReferencePictureSet(rcListPic, pcSlice->getRPS());

    if(pcSlice->getTLayer() > 0 
      &&  !( pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_RADL_N     // Check if not a leading picture
          || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_RADL_R
          || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_RASL_N
          || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_RASL_R )
        )
    {
      if(pcSlice->isTemporalLayerSwitchingPoint(rcListPic) || pcSlice->getSPS()->getTemporalIdNestingFlag())
      {
        if(pcSlice->getTemporalLayerNonReferenceFlag())
        {
          pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_TSA_N);
        }
        else
        {
          pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_TSA_R);
        }
      }
      else if(pcSlice->isStepwiseTemporalLayerSwitchingPointCandidate(rcListPic))
      {
        Bool isSTSA=true;
        for(Int ii=iGOPid+1;(ii<m_pcCfg->getGOPSize() && isSTSA==true);ii++)
        {
          Int lTid= m_pcCfg->getGOPEntry(ii).m_temporalId;
          if(lTid==pcSlice->getTLayer())
          {
            const TComReferencePictureSet* nRPS = pcSlice->getSPS()->getRPSList()->getReferencePictureSet(ii);
            for(Int jj=0;jj<nRPS->getNumberOfPictures();jj++)
            {
              if(nRPS->getUsed(jj))
              {
                Int tPoc=m_pcCfg->getGOPEntry(ii).m_POC+nRPS->getDeltaPOC(jj);
                Int kk=0;
                for(kk=0;kk<m_pcCfg->getGOPSize();kk++)
                {
                  if(m_pcCfg->getGOPEntry(kk).m_POC==tPoc)
                  {
                    break;
                  }
                }
                Int tTid=m_pcCfg->getGOPEntry(kk).m_temporalId;
                if(tTid >= pcSlice->getTLayer())
                {
                  isSTSA=false;
                  break;
                }
              }
            }
          }
        }
        if(isSTSA==true)
        {
          if(pcSlice->getTemporalLayerNonReferenceFlag())
          {
            pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_STSA_N);
          }
          else
          {
            pcSlice->setNalUnitType(NAL_UNIT_CODED_SLICE_STSA_R);
          }
        }
      }
    }
    arrangeLongtermPicturesInRPS(pcSlice, rcListPic);
    TComRefPicListModification* refPicListModification = pcSlice->getRefPicListModification();
    refPicListModification->setRefPicListModificationFlagL0(0);
    refPicListModification->setRefPicListModificationFlagL1(0);
    pcSlice->setNumRefIdx(REF_PIC_LIST_0,min(m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive,pcSlice->getRPS()->getNumberOfPictures()));
    pcSlice->setNumRefIdx(REF_PIC_LIST_1,min(m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive,pcSlice->getRPS()->getNumberOfPictures()));

#if BLOCK_BASED_BACKGROUND_REFERENCE
		if (pcSlice->getSPS()->getBackgroundRefsPresent() && pcSlice->getSliceType() != I_SLICE)
		{
			ArrangeBackgroundReference(pcSlice, rcListPic, pocCurr);
			pcSlice->setNumRefIdx(REF_PIC_LIST_0, min(m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive, pcSlice->getRPS()->getNumberOfPictures() + 1));
			pcSlice->setNumRefIdx(REF_PIC_LIST_1, min(m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive, pcSlice->getRPS()->getNumberOfPictures() + 1));
		}
#endif

    //  Set reference list
    pcSlice->setRefPicList ( rcListPic );

#if JVET_C0024_AMAX_BT
    if (!pcSlice->isIntra() ) 
    {
      Int refLayer=pcSlice->getDepth();
      if( refLayer>9) refLayer=9; // Max layer is 10  
#if JVET_C0024_AMAX_BT_FIX
      if( g_bInitAMaxBT && pcSlice->getPOC() > g_uiPrevISlicePOC )
      {
        ::memset( g_uiBlkSize, 0, sizeof(g_uiBlkSize) );
        ::memset( g_uiNumBlk, 0, sizeof(g_uiNumBlk) );
        g_bInitAMaxBT = false;
      }
#endif
      if (refLayer >= 0 && g_uiNumBlk[refLayer] != 0) 
      {
        Double dBlkSize = sqrt((Double)g_uiBlkSize[refLayer]/g_uiNumBlk[refLayer]);
        if (dBlkSize < AMAXBT_TH32)
        {
          pcSlice->setMaxBTSize(32>MAX_BT_SIZE_INTER ? MAX_BT_SIZE_INTER: 32);
        }
        else if (dBlkSize < AMAXBT_TH64)
        {
          pcSlice->setMaxBTSize(64>MAX_BT_SIZE_INTER ? MAX_BT_SIZE_INTER: 64);
        }
        else
        {
          pcSlice->setMaxBTSize(128>MAX_BT_SIZE_INTER? MAX_BT_SIZE_INTER: 128);
        }
#if !JVET_C0024_AMAX_BT_FIX
        printf("\n previous layer=%d, avg blk size = %3.2f, current max BT set to %d\n", refLayer, dBlkSize, pcSlice->getMaxBTSize());
#endif

        g_uiBlkSize[refLayer] = 0;
        g_uiNumBlk[refLayer] = 0;
      }
    }
#if JVET_C0024_AMAX_BT_FIX
    else
    {
#if JVET_C0024_AMAX_BT_FIX_TICKET23
      if( g_bInitAMaxBT  )
      {
        ::memset( g_uiBlkSize, 0, sizeof(g_uiBlkSize) );
        ::memset( g_uiNumBlk, 0, sizeof(g_uiNumBlk) );
      }
#endif
      g_uiPrevISlicePOC = pcSlice->getPOC();
      g_bInitAMaxBT = true;
    }
#endif
#endif


    //  Slice info. refinement
    if ( (pcSlice->getSliceType() == B_SLICE) && (pcSlice->getNumRefIdx(REF_PIC_LIST_1) == 0) )
    {
      pcSlice->setSliceType ( P_SLICE );
    }
#if PARALLEL_ENCODING_RAS_CABAC_INIT_PRESENT  
    // Prevent encoder from using cross RAP information
    if ((pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_TRAIL_N || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_TRAIL_R || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_TSA_N || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_TSA_R || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_STSA_R || pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_STSA_N) && // current slice belongs to trailing picture 
      (m_pcSliceEncoder->getLastNALUType() == NAL_UNIT_CODED_SLICE_RADL_N || m_pcSliceEncoder->getLastNALUType() == NAL_UNIT_CODED_SLICE_RADL_R || m_pcSliceEncoder->getLastNALUType() == NAL_UNIT_CODED_SLICE_RASL_N || m_pcSliceEncoder->getLastNALUType() == NAL_UNIT_CODED_SLICE_RASL_R))      // previous slice belongs to leading picture
    {
      m_pcSliceEncoder->setEncCABACTableIdx(pcSlice->getSliceType());
    }
#endif
    pcSlice->setEncCABACTableIdx(m_pcSliceEncoder->getEncCABACTableIdx());


    if (pcSlice->getSliceType() == B_SLICE)
    {
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
      const UInt uiColFromL0 = calculateCollocatedFromL0Flag(pcSlice);
      pcSlice->setColFromL0Flag(uiColFromL0);
#else
      pcSlice->setColFromL0Flag(1-uiColDir);
#endif
      Bool bLowDelay = true;
      Int  iCurrPOC  = pcSlice->getPOC();
      Int iRefIdx = 0;

      for (iRefIdx = 0; iRefIdx < pcSlice->getNumRefIdx(REF_PIC_LIST_0) && bLowDelay; iRefIdx++)
      {
        if ( pcSlice->getRefPic(REF_PIC_LIST_0, iRefIdx)->getPOC() > iCurrPOC )
        {
          bLowDelay = false;
        }
      }
      for (iRefIdx = 0; iRefIdx < pcSlice->getNumRefIdx(REF_PIC_LIST_1) && bLowDelay; iRefIdx++)
      {
        if ( pcSlice->getRefPic(REF_PIC_LIST_1, iRefIdx)->getPOC() > iCurrPOC )
        {
          bLowDelay = false;
        }
      }

      pcSlice->setCheckLDC(bLowDelay);
    }
    else
    {
      pcSlice->setCheckLDC(true);
    }

#if JVET_E0023_FAST_ENCODING_SETTING
    UInt dist = MAX_UINT;
    if ( !pcSlice->isIntra() )
    {
      Int currPOC  = pcSlice->getPOC();
      for (Int refIdx = 0; refIdx < pcSlice->getNumRefIdx(REF_PIC_LIST_0); refIdx++)
      {
        UInt tmp = abs(currPOC - pcSlice->getRefPic(REF_PIC_LIST_0, refIdx)->getPOC());
        if ( tmp < dist )
        {
          dist = tmp;
        }
      }
      if ( pcSlice->getSliceType() == B_SLICE )
      {
        for (Int refIdx = 0; refIdx < pcSlice->getNumRefIdx(REF_PIC_LIST_1); refIdx++)
        {
          UInt tmp = abs(currPOC - pcSlice->getRefPic(REF_PIC_LIST_1, refIdx)->getPOC());
          if ( tmp < dist )
          {
            dist = tmp;
          }
        }
      }
    }
    pcSlice->setPictureDistance(dist);
#endif

#if !JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    uiColDir = 1-uiColDir;
#endif
    //-------------------------------------------------------------
    pcSlice->setRefPOCList();
#if JVET_C0027_BIO
  if (pcSlice->getSliceType() != B_SLICE) // there is no bi-pred
  {
   pcSlice->setBioLDBPossible(false);   
  }
  else if (pcSlice->getNumRefIdx(REF_PIC_LIST_0)<=1 && pcSlice->getNumRefIdx(REF_PIC_LIST_1)<=1)
  {
    Int  iPOCcurrent = pcSlice->getPOC();
    Int  iPOCL0 =  pcSlice->getRefPOC(REF_PIC_LIST_0,0); 
    Int  iPOCL1 =  pcSlice->getRefPOC(REF_PIC_LIST_1,0); 
    if (abs (iPOCL0-iPOCcurrent) ==1 && abs(iPOCL1-iPOCcurrent)==1)
    {
       pcSlice->setBioLDBPossible(true);  // this is LDB
    }
    else
    {
       pcSlice->setBioLDBPossible(false); // this is RA
    }
  }

  else   // check wheather bi-pred from different time direction is possible
  {
    pcSlice->setBioLDBPossible(pcSlice->getCheckLDC()); // this is LDB
  }
#endif
    pcSlice->setList1IdxToList0Idx();

    if (m_pcEncTop->getTMVPModeId() == 2)
    {
      if (iGOPid == 0) // first picture in SOP (i.e. forward B)
      {
        pcSlice->setEnableTMVPFlag(0);
      }
      else
      {
        // Note: pcSlice->getColFromL0Flag() is assumed to be always 0 and getcolRefIdx() is always 0.
        pcSlice->setEnableTMVPFlag(1);
      }
    }
    else if (m_pcEncTop->getTMVPModeId() == 1)
    {
      pcSlice->setEnableTMVPFlag(1);
    }
    else
    {
      pcSlice->setEnableTMVPFlag(0);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////// Compress a slice
    //  Slice compression
    if (m_pcCfg->getUseASR() && pcSlice->getSliceType() != I_SLICE)
    {
      m_pcSliceEncoder->setSearchRange(pcSlice);
    }

    Bool bGPBcheck=false;
    if ( pcSlice->getSliceType() == B_SLICE)
    {
      if ( pcSlice->getNumRefIdx(RefPicList( 0 ) ) == pcSlice->getNumRefIdx(RefPicList( 1 ) ) )
      {
        bGPBcheck=true;
        Int i;
        for ( i=0; i < pcSlice->getNumRefIdx(RefPicList( 1 ) ); i++ )
        {
          if ( pcSlice->getRefPOC(RefPicList(1), i) != pcSlice->getRefPOC(RefPicList(0), i) )
          {
            bGPBcheck=false;
            break;
          }
        }
      }
    }
    if(bGPBcheck)
    {
      pcSlice->setMvdL1ZeroFlag(true);
    }
    else
    {
      pcSlice->setMvdL1ZeroFlag(false);
    }
    pcPic->getSlice(pcSlice->getSliceIdx())->setMvdL1ZeroFlag(pcSlice->getMvdL1ZeroFlag());


    Double lambda            = 0.0;
    Int actualHeadBits       = 0;
    Int actualTotalBits      = 0;
    Int estimatedBits        = 0;
    Int tmpBitsBeforeWriting = 0;
    if ( m_pcCfg->getUseRateCtrl() ) // TODO: does this work with multiple slices and slice-segments?
    {
      Int frameLevel = m_pcRateCtrl->getRCSeq()->getGOPID2Level( iGOPid );
      if ( pcPic->getSlice(0)->getSliceType() == I_SLICE )
      {
        frameLevel = 0;
      }
      m_pcRateCtrl->initRCPic( frameLevel );
      estimatedBits = m_pcRateCtrl->getRCPic()->getTargetBits();

      Int sliceQP = m_pcCfg->getInitialQP();
      if ( ( pcSlice->getPOC() == 0 && m_pcCfg->getInitialQP() > 0 ) || ( frameLevel == 0 && m_pcCfg->getForceIntraQP() ) ) // QP is specified
      {
        Int    NumberBFrames = ( m_pcCfg->getGOPSize() - 1 );
        Double dLambda_scale = 1.0 - Clip3( 0.0, 0.5, 0.05*(Double)NumberBFrames );
        Double dQPFactor     = 0.57*dLambda_scale;
        Int    SHIFT_QP      = 12;
        Int    bitdepth_luma_qp_scale = 0;
        Double qp_temp = (Double) sliceQP + bitdepth_luma_qp_scale - SHIFT_QP;
        lambda = dQPFactor*pow( 2.0, qp_temp/3.0 );
      }
      else if ( frameLevel == 0 )   // intra case, but use the model
      {
        m_pcSliceEncoder->calCostSliceI(pcPic); // TODO: This only analyses the first slice segment - what about the others?

        if ( m_pcCfg->getIntraPeriod() != 1 )   // do not refine allocated bits for all intra case
        {
          Int bits = m_pcRateCtrl->getRCSeq()->getLeftAverageBits();
          bits = m_pcRateCtrl->getRCPic()->getRefineBitsForIntra( bits );
          if ( bits < 200 )
          {
            bits = 200;
          }
          m_pcRateCtrl->getRCPic()->setTargetBits( bits );
        }

        list<TEncRCPic*> listPreviousPicture = m_pcRateCtrl->getPicList();
        m_pcRateCtrl->getRCPic()->getLCUInitTargetBits();
        lambda  = m_pcRateCtrl->getRCPic()->estimatePicLambda( listPreviousPicture, pcSlice->getSliceType());
        sliceQP = m_pcRateCtrl->getRCPic()->estimatePicQP( lambda, listPreviousPicture );
      }
      else    // normal case
      {
        list<TEncRCPic*> listPreviousPicture = m_pcRateCtrl->getPicList();
        lambda  = m_pcRateCtrl->getRCPic()->estimatePicLambda( listPreviousPicture, pcSlice->getSliceType());
        sliceQP = m_pcRateCtrl->getRCPic()->estimatePicQP( lambda, listPreviousPicture );
      }

      sliceQP = Clip3( -pcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, sliceQP );
      m_pcRateCtrl->getRCPic()->setPicEstQP( sliceQP );

      m_pcSliceEncoder->resetQP( pcPic, sliceQP, lambda );
    }

    UInt uiNumSliceSegments = 1;

    // Allocate some coders, now the number of tiles are known.
    const Int numSubstreamsColumns = (pcSlice->getPPS()->getNumTileColumnsMinus1() + 1);
    const Int numSubstreamRows     = pcSlice->getPPS()->getEntropyCodingSyncEnabledFlag() ? pcPic->getFrameHeightInCtus() : (pcSlice->getPPS()->getNumTileRowsMinus1() + 1);
    const Int numSubstreams        = numSubstreamRows * numSubstreamsColumns;
    std::vector<TComOutputBitstream> substreamsOut(numSubstreams);

#if VCEG_AZ07_FRUC_MERGE
    if( pcSlice->getSPS()->getUseFRUCMgrMode() && !pcSlice->isIntra() )
    {
      pcPic->initFRUCMVP();
    }
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
        // set adaptive clipping bounds for current slice
        if (m_pcCfg->getTchClipParam().isActive ) {

            Int tbd;
            if (pcSlice->getSliceType() == I_SLICE) tbd=-1;
            else                                    tbd=pcSlice->getDepth();
            Int delta_disto_luma,delta_disto_chroma;
            ClipParam prm=pcPic->computeTchClipParam(delta_disto_luma,delta_disto_chroma);

            prm=codingChoice(prm,delta_disto_luma,delta_disto_chroma,pcPic->getSlice(0)->getLambdas()[0],pcPic->getSlice(0)->getLambdas()[1],tbd);
            if (prm.isActive) {
                pcPic->m_aclip_prm     = prm;
            } else {
                setOff(pcPic->m_aclip_prm);
            }
        } else {
            setOff(pcPic->m_aclip_prm); // OFF with defaults val
        }
        g_ClipParam =pcPic->m_aclip_prm; // set the global for access from clipBD

#if RRC_DWDIP_MLR_IMPR
        UInt uiCurTuLengthLuma = 1 << (MIN_CU_LOG2 - 1);
        UInt uiCurTuLengthChroma = 1 << (MIN_CU_LOG2 - 2);
        for (UInt uiCnt = 0; uiCnt < BLK_LENGTH_VALUES; ++uiCnt)
        {
          g_ClipParam.Y_.dMaxAngle[uiCnt] = CLIP_FACTOR * TComPrediction::getMaxTangent(g_ClipParam.Y().m, g_ClipParam.Y().M, (uiCurTuLengthLuma <<= 1));
          g_ClipParam.U_.dMaxAngle[uiCnt] = CLIP_FACTOR * TComPrediction::getMaxTangent(g_ClipParam.U().m, g_ClipParam.U().M, (uiCurTuLengthChroma <<= 1));
          g_ClipParam.V_.dMaxAngle[uiCnt] = CLIP_FACTOR * TComPrediction::getMaxTangent(g_ClipParam.V().m, g_ClipParam.V().M, uiCurTuLengthChroma);
        }
#endif
#endif
    // now compress (trial encode) the various slice segments (slices, and dependent slices)
    {
      const UInt numberOfCtusInFrame=pcPic->getPicSym()->getNumberOfCtusInFrame();
      pcSlice->setSliceCurStartCtuTsAddr( 0 );
      pcSlice->setSliceSegmentCurStartCtuTsAddr( 0 );

      for(UInt nextCtuTsAddr = 0; nextCtuTsAddr < numberOfCtusInFrame; )
      {
        m_pcSliceEncoder->precompressSlice( pcPic );
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
        pcSlice->setStatsHandle( m_apcStats );
        pcSlice->initStatsGlobal( );              
#endif  
#if BARC_CORE 
#if BARC_ONLY_FOR_INTRA
		if (pcSlice->getSliceType() == I_SLICE)
		{
#endif
			m_pcSliceEncoder->xDownScalePic(pcPic->getPicYuvOrg(), pcPic->getPicYuvDScaleOrg());
#if BARC_ONLY_FOR_INTRA
		}
#endif
#endif
        m_pcSliceEncoder->compressSlice   ( pcPic, false, false );

        const UInt curSliceSegmentEnd = pcSlice->getSliceSegmentCurEndCtuTsAddr();
        if (curSliceSegmentEnd < numberOfCtusInFrame)
        {
          const Bool bNextSegmentIsDependentSlice=curSliceSegmentEnd<pcSlice->getSliceCurEndCtuTsAddr();
          const UInt sliceBits=pcSlice->getSliceBits();
          pcPic->allocateNewSlice();
          // prepare for next slice
          pcPic->setCurrSliceIdx                    ( uiNumSliceSegments );
          m_pcSliceEncoder->setSliceIdx             ( uiNumSliceSegments   );
          pcSlice = pcPic->getSlice                 ( uiNumSliceSegments   );
          assert(pcSlice->getPPS()!=0);
          pcSlice->copySliceInfo                    ( pcPic->getSlice(uiNumSliceSegments-1)  );
          pcSlice->setSliceIdx                      ( uiNumSliceSegments   );
          if (bNextSegmentIsDependentSlice)
          {
            pcSlice->setSliceBits(sliceBits);
          }
          else
          {
            pcSlice->setSliceCurStartCtuTsAddr      ( curSliceSegmentEnd );
            pcSlice->setSliceBits(0);
          }
          pcSlice->setDependentSliceSegmentFlag(bNextSegmentIsDependentSlice);
          pcSlice->setSliceSegmentCurStartCtuTsAddr ( curSliceSegmentEnd );
          // TODO: optimise cabac_init during compress slice to improve multi-slice operation
          // pcSlice->setEncCABACTableIdx(m_pcSliceEncoder->getEncCABACTableIdx());
          uiNumSliceSegments ++;
        }
        nextCtuTsAddr = curSliceSegmentEnd;
      }
    }

    duData.clear();
    pcSlice = pcPic->getSlice(0);

#if TRACE_YUV
    TRACE_PIC(pcPic->getPicYuvRec(),  "t_rec.yuv");
    TRACE_PIC(pcPic->getPicYuvPred(), "t_predicted_top.yuv");
#endif
    // SAO parameter estimation using non-deblocked pixels for CTU bottom and right boundary areas
    if( pcSlice->getSPS()->getUseSAO() && m_pcCfg->getSaoCtuBoundary() )
    {
      m_pcSAO->getPreDBFStatistics(pcPic);
    }

    //-- Loop filter
    Bool bLFCrossTileBoundary = pcSlice->getPPS()->getLoopFilterAcrossTilesEnabledFlag();
    m_pcLoopFilter->setCfg(bLFCrossTileBoundary);
    if ( m_pcCfg->getDeblockingFilterMetric() )
    {
      applyDeblockingFilterMetric(pcPic, uiNumSliceSegments);
    }
    m_pcLoopFilter->loopFilterPic( pcPic );

#if BLOCK_BASED_BACKGROUND_REFERENCE
		if (pcSlice->getQPList() != NULL)
		{
			delete[] pcSlice->getQPList();
		}
#endif

    /////////////////////////////////////////////////////////////////////////////////////////////////// File writing
    // Set entropy coder
    m_pcEntropyCoder->setEntropyCoder   ( m_pcCavlcCoder );

#if JVET_D0135_PARAMS
    //if ( m_bSeqFirst || (m_pcCfg->getWriteParamSetsIDRFlag() && pcSlice->getIdrPicFlag()) )
    if ( m_bSeqFirst || (m_pcCfg->getReWriteParamSetsFlag() && ( pcPic->getSlice(0)->getSliceType() == I_SLICE )))
#else
      if ( m_bSeqFirst )
#endif
    {
      // write various parameter sets
      actualTotalBits += xWriteParameterSets(accessUnit, pcSlice);

      // create prefix SEI messages at the beginning of the sequence
      assert(leadingSeiMessages.empty());
      xCreateIRAPLeadingSEIMessages(leadingSeiMessages, pcSlice->getSPS(), pcSlice->getPPS());

      m_bSeqFirst = false;
    }

    // reset presence of BP SEI indication
    m_bufferingPeriodSEIPresentInAU = false;
    // create prefix SEI associated with a picture
    xCreatePerPictureSEIMessages(iGOPid, leadingSeiMessages, nestedSeiMessages, pcSlice);

    /* use the main bitstream buffer for storing the marshalled picture */
    m_pcEntropyCoder->setBitstream(NULL);

    pcSlice = pcPic->getSlice(0);
#if CNN_BASED_INLOOP_FILTER 
#if PRINT_INFO
		Double endTime;
		clock_t startClock = clock();
		printf("POC %d cnn loop filter start.\n", pcPic->getPOC());
#endif
		if (pcSlice->getSPS()->getUseCNNLoopFilter())
		{
			m_pcCNNLoopFilter->convNetProcess(pcPic);
		}
#if PRINT_INFO
		endTime = (Double)(clock() - startClock) / CLOCKS_PER_SEC;
		printf("POC %d enc cnn done! Time: %.3f sec.\n", pcPic->getPOC(), endTime);
#endif
#endif
    if (pcSlice->getSPS()->getUseSAO())
    {
      Bool sliceEnabled[MAX_NUM_COMPONENT];
      TComBitCounter tempBitCounter;
      tempBitCounter.resetBits();
      m_pcEncTop->getRDGoOnSbacCoder()->setBitstream(&tempBitCounter);
#if VCEG_AZ07_BAC_ADAPT_WDOW
      m_pcSAO->setEntropyCoder(m_pcEntropyCoder);
#endif
      m_pcSAO->initRDOCabacCoder(m_pcEncTop->getRDGoOnSbacCoder(), pcSlice);
      m_pcSAO->SAOProcess(pcPic, sliceEnabled, pcPic->getSlice(0)->getLambdas(), m_pcCfg->getTestSAODisableAtPictureLevel(), m_pcCfg->getSaoEncodingRate(), m_pcCfg->getSaoEncodingRateChroma(), m_pcCfg->getSaoCtuBoundary());
      m_pcSAO->PCMLFDisableProcess(pcPic);
      m_pcEncTop->getRDGoOnSbacCoder()->setBitstream(NULL);

      //assign SAO slice header
      for(Int s=0; s< uiNumSliceSegments; s++)
      {
        pcPic->getSlice(s)->setSaoEnabledFlag(CHANNEL_TYPE_LUMA, sliceEnabled[COMPONENT_Y]);
        assert(sliceEnabled[COMPONENT_Cb] == sliceEnabled[COMPONENT_Cr]);
        pcPic->getSlice(s)->setSaoEnabledFlag(CHANNEL_TYPE_CHROMA, sliceEnabled[COMPONENT_Cb]);
      }
    }

#if RRC_NS
		if (pcSlice->getSPS()->getUseNS() && pcPic->getTLayer() < RRC_NS_TLAYER)
			filterNs(pcPic, pcSlice);
#endif

#if ALF_HM3_REFACTOR
    if( pcSlice->getSPS()->getUseALF() )
    {
      m_pcAdaptiveLoopFilter->setNumCUsInFrame(pcPic);
      TComBitCounter tempBitCounter;
      tempBitCounter.resetBits();
      m_pcEncTop->getRDGoOnSbacCoder()->setBitstream(&tempBitCounter);
      m_pcEntropyCoder->setEntropyCoder ( m_pcEncTop->getRDGoOnSbacCoder() );
      m_pcAdaptiveLoopFilter->startALFEnc(pcPic, m_pcEntropyCoder );
      UInt64 uiDist, uiBits;
      if( !bInitAlfParam )
      {
        m_pcAdaptiveLoopFilter->setNumCUsInFrame(pcPic);
        m_pcAdaptiveLoopFilter->allocALFParam(&cAlfParam);
        bInitAlfParam = true;
      }
      m_pcAdaptiveLoopFilter->resetALFParam( &cAlfParam );

      /*TComPicYuv* rec = pcPic->getPicYuvRec();
      ofstream file("rec_416x240.yuv", ios::app | ios::binary);
      for (ComponentID comp : {COMPONENT_Y, COMPONENT_Cb, COMPONENT_Cr})
        for (int i = 0; i < rec->getHeight(comp); ++i)
          for (int j = 0; j < rec->getWidth(comp); ++j)
            file << static_cast<char>(rec->getAddr(comp)[j + i * rec->getStride(comp)]);*/

      m_pcAdaptiveLoopFilter->ALFProcess( &cAlfParam, pcPic->getSlice(0)->getLambdas()[0], pcPic->getSlice(0)->getLambdas()[1], uiDist, uiBits, cAlfParam.alf_max_depth
#if COM16_C806_ALF_TEMPPRED_NUM
#if FIX_TICKET12
#if JVET_E0104_ALF_TEMP_SCALABILITY
        , (pcSlice->getSliceType()== I_SLICE ? NULL: m_acStoredAlfPara[pcSlice->getTLayer()]), m_iStoredAlfParaNum[pcSlice->getTLayer()]
#else
        , (pcSlice->getSliceType()== I_SLICE ? NULL: m_acStoredAlfPara), m_iStoredAlfParaNum
#endif
#else
        , m_acStoredAlfPara, m_iStoredAlfParaNum
#endif
#endif
        );

      //ofstream file2("filt_416x240.yuv", ios::app | ios::binary);
      //for (ComponentID comp : {COMPONENT_Y, COMPONENT_Cb, COMPONENT_Cr})
      //  for (int i = 0; i < rec->getHeight(comp); ++i)
      //    for (int j = 0; j < rec->getWidth(comp); ++j)
      //      file2 << static_cast<char>(rec->getAddr(comp)[j + i * rec->getStride(comp)]);

#if COM16_C806_ALF_TEMPPRED_NUM
      if( cAlfParam.alf_flag && !cAlfParam.temproalPredFlag && cAlfParam.filtNo >= 0 )
      {
#if JVET_E0104_ALF_TEMP_SCALABILITY
        Int iCurrTempIdx = pcSlice->getTLayer();
        assert(iCurrTempIdx < JVET_E0104_ALF_MAX_TEMPLAYERID);
        for (Int iLoopedTempIdx = iCurrTempIdx; iLoopedTempIdx < JVET_E0104_ALF_MAX_TEMPLAYERID; iLoopedTempIdx++)
        {
          Int iIdx = m_iStoredAlfParaNum[iLoopedTempIdx] % COM16_C806_ALF_TEMPPRED_NUM;
          m_iStoredAlfParaNum[iLoopedTempIdx] ++;
          m_acStoredAlfPara[iLoopedTempIdx][iIdx].temproalPredFlag = false;
          m_pcAdaptiveLoopFilter->copyALFParam(&m_acStoredAlfPara[iLoopedTempIdx][iIdx], &cAlfParam);
#if JVET_C0038_GALF
          m_pcAdaptiveLoopFilter->resetALFPredParam(&m_acStoredAlfPara[iLoopedTempIdx][iIdx], (pcSlice->getSliceType() == I_SLICE ? true : false));
#endif
        }
#else
        Int iIdx = m_iStoredAlfParaNum % COM16_C806_ALF_TEMPPRED_NUM;
        m_iStoredAlfParaNum++;
        m_acStoredAlfPara[iIdx].temproalPredFlag = false;
        m_pcAdaptiveLoopFilter->copyALFParam( &m_acStoredAlfPara[iIdx], &cAlfParam );
#if JVET_C0038_GALF
        m_pcAdaptiveLoopFilter->resetALFPredParam(&m_acStoredAlfPara[iIdx], (pcSlice->getSliceType()== I_SLICE? true: false));
#endif
#endif
      }
#endif     
      m_pcAdaptiveLoopFilter->endALFEnc();
    }
#endif

    // pcSlice is currently slice 0.
    std::size_t binCountsInNalUnits   = 0; // For implementation of cabac_zero_word stuffing (section 7.4.3.10)
    std::size_t numBytesInVclNalUnits = 0; // For implementation of cabac_zero_word stuffing (section 7.4.3.10)

    for( UInt sliceSegmentStartCtuTsAddr = 0, sliceIdxCount=0; sliceSegmentStartCtuTsAddr < pcPic->getPicSym()->getNumberOfCtusInFrame(); sliceIdxCount++, sliceSegmentStartCtuTsAddr=pcSlice->getSliceSegmentCurEndCtuTsAddr() )
    {
      pcSlice = pcPic->getSlice(sliceIdxCount);
      if(sliceIdxCount > 0 && pcSlice->getSliceType()!= I_SLICE)
      {
        pcSlice->checkColRefIdx(sliceIdxCount, pcPic);
      }
      pcPic->setCurrSliceIdx(sliceIdxCount);
      m_pcSliceEncoder->setSliceIdx(sliceIdxCount);

      pcSlice->setRPS(pcPic->getSlice(0)->getRPS());
      pcSlice->setRPSidx(pcPic->getSlice(0)->getRPSidx());

      for ( UInt ui = 0 ; ui < numSubstreams; ui++ )
      {
        substreamsOut[ui].clear();
      }

      m_pcEntropyCoder->setEntropyCoder   ( m_pcCavlcCoder );
      m_pcEntropyCoder->resetEntropy      ( pcSlice );
      /* start slice NALunit */
      OutputNALUnit nalu( pcSlice->getNalUnitType(), pcSlice->getTLayer() );
      m_pcEntropyCoder->setBitstream(&nalu.m_Bitstream);

      pcSlice->setNoRaslOutputFlag(false);
      if (pcSlice->isIRAP())
      {
        if (pcSlice->getNalUnitType() >= NAL_UNIT_CODED_SLICE_BLA_W_LP && pcSlice->getNalUnitType() <= NAL_UNIT_CODED_SLICE_IDR_N_LP)
        {
          pcSlice->setNoRaslOutputFlag(true);
        }
        //the inference for NoOutputPriorPicsFlag
        // KJS: This cannot happen at the encoder
        if (!m_bFirst && pcSlice->isIRAP() && pcSlice->getNoRaslOutputFlag())
        {
          if (pcSlice->getNalUnitType() == NAL_UNIT_CODED_SLICE_CRA)
          {
            pcSlice->setNoOutputPriorPicsFlag(true);
          }
        }
      }

      pcSlice->setEncCABACTableIdx(m_pcSliceEncoder->getEncCABACTableIdx());

      tmpBitsBeforeWriting = m_pcEntropyCoder->getNumberOfWrittenBits();
      m_pcEntropyCoder->encodeSliceHeader(pcSlice);
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
#if VCEG_AZ07_BAC_ADAPT_WDOW 
      m_pcEntropyCoder->setStatsHandle( m_apcStats );       
      m_pcEntropyCoder->encodeCtxUpdateInfo( pcSlice, m_apcStats );
#endif
      Int iQPIdx = xUpdateTStates (pcSlice->getSliceType(), pcSlice->getSliceQp(), m_apcStats);
      pcSlice->setQPIdx(iQPIdx);
#endif

      actualHeadBits += ( m_pcEntropyCoder->getNumberOfWrittenBits() - tmpBitsBeforeWriting );

      pcSlice->setFinalized(true);

      pcSlice->clearSubstreamSizes(  );

      {
        UInt numBinsCoded = 0;
        m_pcSliceEncoder->encodeSlice(pcPic, &(substreamsOut[0]), numBinsCoded
#if ALF_HM3_REFACTOR
          , cAlfParam
#endif
          );
        binCountsInNalUnits+=numBinsCoded;
      }

      {
        // Construct the final bitstream by concatenating substreams.
        // The final bitstream is either nalu.m_Bitstream or pcBitstreamRedirect;
        // Complete the slice header info.
        m_pcEntropyCoder->setEntropyCoder   ( m_pcCavlcCoder );
        m_pcEntropyCoder->setBitstream(&nalu.m_Bitstream);
        m_pcEntropyCoder->encodeTilesWPPEntryPoint( pcSlice );

        // Append substreams...
        TComOutputBitstream *pcOut = pcBitstreamRedirect;
        const Int numZeroSubstreamsAtStartOfSlice  = pcPic->getSubstreamForCtuAddr(pcSlice->getSliceSegmentCurStartCtuTsAddr(), false, pcSlice);
        const Int numSubstreamsToCode  = pcSlice->getNumberOfSubstreamSizes()+1;
        for ( UInt ui = 0 ; ui < numSubstreamsToCode; ui++ )
        {
          pcOut->addSubstream(&(substreamsOut[ui+numZeroSubstreamsAtStartOfSlice]));
        }
      }

      // If current NALU is the first NALU of slice (containing slice header) and more NALUs exist (due to multiple dependent slices) then buffer it.
      // If current NALU is the last NALU of slice and a NALU was buffered, then (a) Write current NALU (b) Update an write buffered NALU at approproate location in NALU list.
      Bool bNALUAlignedWrittenToList    = false; // used to ensure current NALU is not written more than once to the NALU list.
      xAttachSliceDataToNalUnit(nalu, pcBitstreamRedirect);
      accessUnit.push_back(new NALUnitEBSP(nalu));
      actualTotalBits += UInt(accessUnit.back()->m_nalUnitData.str().size()) * 8;
      numBytesInVclNalUnits += (std::size_t)(accessUnit.back()->m_nalUnitData.str().size());
      bNALUAlignedWrittenToList = true;

      if (!bNALUAlignedWrittenToList)
      {
        nalu.m_Bitstream.writeAlignZero();
        accessUnit.push_back(new NALUnitEBSP(nalu));
      }

      if( ( m_pcCfg->getPictureTimingSEIEnabled() || m_pcCfg->getDecodingUnitInfoSEIEnabled() ) &&
          ( pcSlice->getSPS()->getVuiParametersPresentFlag() ) &&
          ( ( pcSlice->getSPS()->getVuiParameters()->getHrdParameters()->getNalHrdParametersPresentFlag() )
         || ( pcSlice->getSPS()->getVuiParameters()->getHrdParameters()->getVclHrdParametersPresentFlag() ) ) &&
          ( pcSlice->getSPS()->getVuiParameters()->getHrdParameters()->getSubPicCpbParamsPresentFlag() ) )
      {
          UInt numNalus = 0;
        UInt numRBSPBytes = 0;
        for (AccessUnit::const_iterator it = accessUnit.begin(); it != accessUnit.end(); it++)
        {
          numRBSPBytes += UInt((*it)->m_nalUnitData.str().size());
          numNalus ++;
        }
        duData.push_back(DUData());
        duData.back().accumBitsDU = ( numRBSPBytes << 3 );
        duData.back().accumNalsDU = numNalus;
      }
    } // end iteration over slices

    // cabac_zero_words processing
    cabac_zero_word_padding(pcSlice, pcPic, binCountsInNalUnits, numBytesInVclNalUnits, accessUnit.back()->m_nalUnitData, m_pcCfg->getCabacZeroWordPaddingEnabled());
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
    if ( !pcSlice->getSPS()->getAtmvpEnableFlag())
    {
      pcPic->compressMotion();
    }
#else  
    pcPic->compressMotion();
#endif
    //-- For time output for each slice
    Double dEncTime = (Double)(clock()-iBeforeTime) / CLOCKS_PER_SEC;

    std::string digestStr;
    if (m_pcCfg->getDecodedPictureHashSEIEnabled())
    {
      SEIDecodedPictureHash *decodedPictureHashSei = new SEIDecodedPictureHash();
      m_seiEncoder.initDecodedPictureHashSEI(decodedPictureHashSei, pcPic, digestStr, pcSlice->getSPS()->getBitDepths());
      trailingSeiMessages.push_back(decodedPictureHashSei);
    }
    xWriteTrailingSEIMessages(trailingSeiMessages, accessUnit, pcSlice->getTLayer(), pcSlice->getSPS());

    m_pcCfg->setEncodedFlag(iGOPid, true);

#if JVET_F0064_MSSSIM
    xCalculateAddPSNRs( isField, isTff, iGOPid, pcPic, accessUnit, rcListPic, dEncTime, snr_conversion, printFrameMSE, printMSSSIM );
#else
    xCalculateAddPSNRs( isField, isTff, iGOPid, pcPic, accessUnit, rcListPic, dEncTime, snr_conversion, printFrameMSE );
#endif

    if (!digestStr.empty())
    {
      if(m_pcCfg->getDecodedPictureHashSEIEnabled() == 1)
      {
        printf(" [MD5:%s]", digestStr.c_str());
      }
      else if(m_pcCfg->getDecodedPictureHashSEIEnabled() == 2)
      {
        printf(" [CRC:%s]", digestStr.c_str());
      }
      else if(m_pcCfg->getDecodedPictureHashSEIEnabled() == 3)
      {
        printf(" [Checksum:%s]", digestStr.c_str());
      }
    }

    if ( m_pcCfg->getUseRateCtrl() )
    {
      Double avgQP     = m_pcRateCtrl->getRCPic()->calAverageQP();
      Double avgLambda = m_pcRateCtrl->getRCPic()->calAverageLambda();
      if ( avgLambda < 0.0 )
      {
        avgLambda = lambda;
      }

      m_pcRateCtrl->getRCPic()->updateAfterPicture( actualHeadBits, actualTotalBits, avgQP, avgLambda, pcSlice->getSliceType());
      m_pcRateCtrl->getRCPic()->addToPictureLsit( m_pcRateCtrl->getPicList() );

      m_pcRateCtrl->getRCSeq()->updateAfterPic( actualTotalBits );
      if ( pcSlice->getSliceType() != I_SLICE )
      {
        m_pcRateCtrl->getRCGOP()->updateAfterPicture( actualTotalBits );
      }
      else    // for intra picture, the estimated bits are used to update the current status in the GOP
      {
        m_pcRateCtrl->getRCGOP()->updateAfterPicture( estimatedBits );
      }
    }

    xCreatePictureTimingSEI(m_pcCfg->getEfficientFieldIRAPEnabled()?effFieldIRAPMap.GetIRAPGOPid():0, leadingSeiMessages, nestedSeiMessages, duInfoSeiMessages, pcSlice, isField, duData);
    if (m_pcCfg->getScalableNestingSEIEnabled())
    {
      xCreateScalableNestingSEI (leadingSeiMessages, nestedSeiMessages);
    }
    xWriteLeadingSEIMessages(leadingSeiMessages, duInfoSeiMessages, accessUnit, pcSlice->getTLayer(), pcSlice->getSPS(), duData);
    xWriteDuSEIMessages(duInfoSeiMessages, accessUnit, pcSlice->getTLayer(), pcSlice->getSPS(), duData);

#if RRC_SELECTABLE_IF & !SIF_GAUSS_ONLY
    pcPic->getSharpeningFilter() = new TComSharpeningFilter(pcPic->getPicYuvRec(), pcSlice->getSPS());
    pcPic->getSharpeningFilter()->buildEdgeMap(pcPic->getPicYuvRec());
#endif
    pcPic->getPicYuvRec()->copyToPic(pcPicYuvRecOut);

    pcPic->setReconMark   ( true );
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
    if (pcSlice->getSPS()->getUseInterKLT())
    {
#endif
        m_pcSliceEncoder->InterpolatePic(pcPic);
#if VCEG_AZ08_USE_KLT
    }
#endif
#endif
    m_bFirst = false;
    m_iNumPicCoded++;
    m_totalCoded ++;
    /* logging: insert a newline at end of picture period */
    printf("\n");
    fflush(stdout);

    if (m_pcCfg->getEfficientFieldIRAPEnabled())
    {
      iGOPid=effFieldIRAPMap.restoreGOPid(iGOPid);
	}
#if TRACE_FLAGS
	gNumFrame++;
#endif
  } // iGOPid-loop

  delete pcBitstreamRedirect;
#if ALF_HM3_REFACTOR
  if( pcSlice->getSPS()->getUseALF() )
  {
    m_pcAdaptiveLoopFilter->freeALFParam(&cAlfParam);
  }
#endif
  assert ( (m_iNumPicCoded == iNumPicRcvd) );
}

#if BLOCK_BASED_BACKGROUND_REFERENCE
// Function for updating the background reference with high quality CTU
// and detecting the background patches to be encoded with high quality in the current picture.
Void TEncGOP::ArrangeBackgroundReference(TComSlice* pcSlice, TComList<TComPic*>& rcListPic, Int pocCurr)
{
	if (m_pcCfg->getFramesToBeEncoded() - pcSlice->getPOC() < 48)
	{
		return;
	}
	TComPic* rpcPicBg = pcSlice->getSPS()->getPicBg();
	TComPic* rpcPrevPic;
	TComPic* rpcCurPic;
	TComList<TComPic*>::iterator  iterPic = rcListPic.begin();
	while (iterPic != rcListPic.end())
	{
		rpcPrevPic = *(iterPic++);
		if (rpcPrevPic->getPOC() == pocCurr - 1)
		{
			break;
		}
	}
	assert(rpcPrevPic->getPOC() == pocCurr - 1);
	iterPic = rcListPic.begin();
	while (iterPic != rcListPic.end())
	{
		rpcCurPic = *(iterPic++);
		if (rpcCurPic->getPOC() == pocCurr)
		{
			break;
		}
	}
	assert(rpcCurPic->getPOC() == pocCurr);

	Int iWidth = rpcCurPic->getPicYuvOrg()->getWidth(COMPONENT_Y);
	Int iHeight = rpcCurPic->getPicYuvOrg()->getHeight(COMPONENT_Y);
	Int iStride = rpcCurPic->getPicYuvOrg()->getStride(COMPONENT_Y);
	Int iCStride = rpcCurPic->getPicYuvOrg()->getStride(COMPONENT_Cb);
	Pel* pBgLumaAddr = rpcPicBg->getPicYuvRec()->getAddr(COMPONENT_Y);
	Pel* pBgCbAddr = rpcPicBg->getPicYuvRec()->getAddr(COMPONENT_Cb);
	Pel* pBgCrAddr = rpcPicBg->getPicYuvRec()->getAddr(COMPONENT_Cr);
	Pel* pPrevLumaAddr = rpcPrevPic->getPicYuvRec()->getAddr(COMPONENT_Y);
	Pel* pPrevCbAddr = rpcPrevPic->getPicYuvRec()->getAddr(COMPONENT_Cb);
	Pel* pPrevCrAddr = rpcPrevPic->getPicYuvRec()->getAddr(COMPONENT_Cr);

	// Update background reference
	if (pocCurr == 1)
	{
		TComList<TComPic*>::iterator  iterPic = rcListPic.begin();
		TComPic* rpcPictmp;
		while (iterPic != rcListPic.end())
		{
			rpcPictmp = *(iterPic++);
			if (rpcPictmp->getPOC() == 0)
			{
				rpcPictmp->getPicYuvRec()->extendPicBorder();
				rpcPictmp->getPicYuvRec()->setBorderExtension(true);
				rpcPictmp->getPicYuvRec()->copyToPic(rpcPicBg->getPicYuvRec());

				for (Int idx = 0; idx < rpcPicBg->getNumberOfCtusInFrame(); idx++)
				{
					rpcPicBg->getCtu(idx)->initCtu(rpcPicBg, idx);
					rpcPicBg->getCtu(idx)->getCUMvField(REF_PIC_LIST_0)->clearMvField();
					rpcPicBg->getCtu(idx)->getCUMvField(REF_PIC_LIST_1)->clearMvField();
					memset(rpcPicBg->getCtu(idx)->getPredictionMode(), MODE_INTRA, rpcPicBg->getNumPartitionsInCtu());
				}
				rpcPictmp->getPicYuvOrg()->extendPicBorder();
				rpcPictmp->getPicYuvOrg()->setBorderExtension(true);
				rpcPictmp->getPicYuvOrg()->copyToPic(rpcPicBg->getPicYuvOrg());
				break;
			}
		}
		assert(rpcPicBg->getPOC() == 0);
	}
	else
	{
		for (Int y = 0; y < iHeight; y += g_uiMaxCUHeight)
		{
			for (Int x = 0; x < iWidth; x += g_uiMaxCUWidth)
			{
				if (rpcPicBg->getSpliceIdx((y / g_uiMaxCUHeight)*rpcPicBg->getFrameWidthInCtus() + x / g_uiMaxCUWidth) == pocCurr - 1)
				{
					for (Int tmpy = 0; tmpy < g_uiMaxCUHeight; tmpy++)
					{
						if (y + tmpy >= iHeight)
						{
							break;
						}
						for (Int tmpx = 0; tmpx < g_uiMaxCUWidth; tmpx++)
						{
							if (x + tmpx >= iWidth)
							{
								break;
							}
							pBgLumaAddr[(y + tmpy)*iStride + x + tmpx] = pPrevLumaAddr[(y + tmpy)*iStride + x + tmpx];
							if (tmpy % 2 == 0 && tmpx % 2 == 0)
							{
								pBgCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] = pPrevCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2];
								pBgCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] = pPrevCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2];
							}
						}
					}
				}
			}
		}
		rpcPicBg->getPicYuvRec()->setBorderExtension(false);
		rpcPicBg->getPicYuvRec()->extendPicBorder();
		rpcPicBg->getPicYuvRec()->setBorderExtension(true);

		Bool IsCompleted = true;
		for (Int idx = 0; idx < rpcPicBg->getNumberOfCtusInFrame(); idx++)
		{
			if (rpcPicBg->getSpliceIdx(idx) == 0)
			{
				IsCompleted = false;
				break;
			}
		}
		if (IsCompleted)
		{
			return;
		}
		if (m_pcCfg->getFramesToBeEncoded() - pcSlice->getPOC() <= 48)
		{
			return;
		}
	}

	//Identify new background patches for later updating
	pPrevLumaAddr = rpcPrevPic->getPicYuvOrg()->getAddr(COMPONENT_Y);
	pPrevCbAddr = rpcPrevPic->getPicYuvOrg()->getAddr(COMPONENT_Cb);
	pPrevCrAddr = rpcPrevPic->getPicYuvOrg()->getAddr(COMPONENT_Cr);
	Pel* pCurLumaAddr = rpcCurPic->getPicYuvOrg()->getAddr(COMPONENT_Y);
	Pel* pCurCbAddr = rpcCurPic->getPicYuvOrg()->getAddr(COMPONENT_Cb);
	Pel* pCurCrAddr = rpcCurPic->getPicYuvOrg()->getAddr(COMPONENT_Cr);
	Pel* pBgOrgLumaAddr = rpcPicBg->getPicYuvOrg()->getAddr(COMPONENT_Y);
	Pel* pBgOrgCbAddr = rpcPicBg->getPicYuvOrg()->getAddr(COMPONENT_Cb);
	Pel* pBgOrgCrAddr = rpcPicBg->getPicYuvOrg()->getAddr(COMPONENT_Cr);
	Int iMaxReplace = rpcPicBg->getNumberOfCtusInFrame() / 10;
	iMaxReplace = iMaxReplace < 1 ? 1 : iMaxReplace;
	typedef struct tagCostStr { Double dCost; Int iCtuIdx; }  CostStr;
	CostStr* pcMinCtuCost = new CostStr[iMaxReplace];
	for (Int i = 0; i < iMaxReplace; i++)
	{
		pcMinCtuCost[i].dCost = 1e10;
		pcMinCtuCost[i].iCtuIdx = -1;
	}

#if BLOCK_BASED_BACKGROUND_REFERENCE
	TComPic* rpcPrevMPicTmp;
	iterPic = rcListPic.begin();
	Int PrevMPOC = pocCurr > 6 ? pocCurr - 6 : 0;
	while (iterPic != rcListPic.end())
	{
		rpcPrevMPicTmp = *(iterPic++);
		if (rpcPrevMPicTmp->getPOC() == PrevMPOC)
		{
			break;
		}
	}
	Pel* pPrevMLumaAddrTmp = rpcPrevMPicTmp->getPicYuvOrg()->getAddr(COMPONENT_Y);
	Pel* pPrevMCbAddrTmp = rpcPrevMPicTmp->getPicYuvOrg()->getAddr(COMPONENT_Cb);
	Pel* pPrevMCrAddrTmp = rpcPrevMPicTmp->getPicYuvOrg()->getAddr(COMPONENT_Cr);
	int iBlockSize = 32;
	int iWidthInBlock = (iWidth + iBlockSize - 1) / (iBlockSize);
	int iHeightInBlock = (iHeight + iBlockSize - 1) / (iBlockSize);
	int iEdgeBlock = 2 * iWidthInBlock + 2 * iHeightInBlock - 4;
	int iTotalBlock = iWidthInBlock*iHeightInBlock;
	int iEdgeBgBlock = 0;
	int iBgBlcok = 0;

	for (Int y = 0; y < iHeight; y += iBlockSize)
	{
		for (Int x = 0; x < iWidth; x += iBlockSize)
		{
			Double dLCUdist = 0.0;
			Double dLCUdistcb = 0.0;
			Double dLCUdistcr = 0.0;
			Int    iLargeDistCnt = 0;
			Int    iRealPixelCnt = 0;
			Double dMean = 0.0;

			for (Int tmpy = 0; tmpy < iBlockSize; tmpy++)
			{
				if (y + tmpy >= iHeight)
				{
					break;
				}
				for (Int tmpx = 0; tmpx < iBlockSize; tmpx++)
				{
					if (x + tmpx >= iWidth)
					{
						break;
					}

					iRealPixelCnt++;
					dLCUdist += abs(pCurLumaAddr[(y + tmpy)*iStride + x + tmpx] - pPrevMLumaAddrTmp[(y + tmpy)*iStride + x + tmpx]);
					dMean += pCurLumaAddr[(y + tmpy)*iStride + x + tmpx];

					if (tmpy % 2 == 0 && tmpx % 2 == 0)
					{
						dLCUdistcb += abs(pCurCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevMCbAddrTmp[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
						dLCUdistcr += abs(pCurCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevMCrAddrTmp[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
					}

					if (abs(pCurLumaAddr[(y + tmpy)*iStride + x + tmpx] - pPrevMLumaAddrTmp[(y + tmpy)*iStride + x + tmpx]) > 16 ||
						abs(pCurCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevMCbAddrTmp[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]) > 8 ||
						abs(pCurCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevMCrAddrTmp[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]) > 8)
					{
						iLargeDistCnt++;
					}
				}
			}

			if (iLargeDistCnt / (Double)iRealPixelCnt < 0.01 && dLCUdist / iRealPixelCnt < 4.0 &&  dLCUdistcb / iRealPixelCnt < 0.5 &&  dLCUdistcr / iRealPixelCnt < 0.5)
			{
				dMean /= iRealPixelCnt;

				double var = 0.0;
				for (int tmpy = 0; tmpy < iBlockSize; tmpy++)
				{
					if (y + tmpy >= iHeight)
					{
						break;
					}
					for (int tmpx = 0; tmpx < iBlockSize; tmpx++)
					{
						if (x + tmpx >= iWidth)
						{
							break;
						}

						var += (pCurLumaAddr[(y + tmpy)*iStride + x + tmpx] - dMean)*(pCurLumaAddr[(y + tmpy)*iStride + x + tmpx] - dMean);
					}
				}

				if (var / iRealPixelCnt > 16.0)
				{
					iBgBlcok++;
					if (x == 0 || y == 0 || x == (iWidthInBlock - 1)*(iBlockSize) || y == (iHeightInBlock - 1)*(iBlockSize))
					{
						iEdgeBgBlock++;
					}
				}
			}
		}
	}
	Bool isStaticPic = false;
	if (!(iBgBlcok > iTotalBlock / 10 || (iBgBlcok > iTotalBlock / 10 && iEdgeBgBlock > iEdgeBlock / 10)))
	{
		return;
	}
#endif

	for (Int y = 0; y < iHeight; y += g_uiMaxCUHeight)
	{
		for (Int x = 0; x < iWidth; x += g_uiMaxCUWidth)
		{
			Double dLCUdist = 0.0;
			Double dLCUdistcb = 0.0;
			Double dLCUdistcr = 0.0;
			Double dDefaultDist = 0.0;
			Double dDefaultDistcb = 0.0;
			Double dDefaultDistcr = 0.0;
			Int    iRealPixelCnt = 0;
			Double dLCUCost = 1e10;
			Bool IsNormalEdge = true;

			for (Int tmpy = 0; tmpy < g_uiMaxCUHeight; tmpy++)
			{
				if (y + tmpy >= iHeight)
				{
					break;
				}
				for (Int tmpx = 0; tmpx < g_uiMaxCUWidth; tmpx++)
				{
					if (x + tmpx >= iWidth)
					{
						break;
					}

					iRealPixelCnt++;
					dLCUdist += abs(pCurLumaAddr[(y + tmpy)*iStride + x + tmpx] - pPrevLumaAddr[(y + tmpy)*iStride + x + tmpx]);
					if (tmpy % 2 == 0 && tmpx % 2 == 0)
					{
						dLCUdistcb += abs(pCurCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
						dLCUdistcr += abs(pCurCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
					}
					dDefaultDist += abs(pBgOrgLumaAddr[(y + tmpy)*iStride + x + tmpx] - pPrevLumaAddr[(y + tmpy)*iStride + x + tmpx]);
					if (tmpy % 2 == 0 && tmpx % 2 == 0)
					{
						dDefaultDistcb += abs(pBgOrgCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevCbAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
						dDefaultDistcr += abs(pBgOrgCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2] - pPrevCrAddr[(y + tmpy) / 2 * iCStride + (x + tmpx) / 2]);
					}
				}
			}

			//check whether default background reference is better
			if (rpcPrevPic->getSlice(0)->getSliceType() != I_SLICE && dLCUdist + dLCUdistcb + dLCUdistcr > dDefaultDist + dDefaultDistcb + dDefaultDistcr)
			{
				continue;
			}

			//Test the vertical or horizontal edge for background patches candidates
			Int iYinLCU = y / g_uiMaxCUHeight;
			Int iXinLCU = x / g_uiMaxCUWidth;
			Int iLCUIdx = iYinLCU *rpcPicBg->getFrameWidthInCtus() + iXinLCU;
			Int iEdgeThreshold = 8;
			Pel* pcEdgeAddrBg;
			Pel* pcEdgeAddrCur;

			for (Int dx = -1; dx <= 1; dx += 2)
			{
				if (iXinLCU + dx < 0 || iXinLCU + dx >= rpcPicBg->getFrameWidthInCtus())
				{
					continue;
				}
				if (rpcPicBg->getSpliceIdx(iYinLCU *rpcPicBg->getFrameWidthInCtus() + iXinLCU + dx) != 0)
				{
					Int iTmpDif = 0;
					Int iTmpDifBg = 0;

					if (dx == -1)
					{
						pcEdgeAddrBg = &pBgLumaAddr[y*iStride + x - 2];
						pcEdgeAddrCur = &pCurLumaAddr[y*iStride + x];
					}
					if (dx == 1)
					{
						pcEdgeAddrBg = &pBgLumaAddr[y*iStride + x + g_uiMaxCUWidth + 1];
						pcEdgeAddrCur = &pCurLumaAddr[y*iStride + x + g_uiMaxCUWidth - 1];
					}

					Int iRealHeight = 0;
					for (Int tmpy = 0; tmpy < g_uiMaxCUHeight - 2; tmpy++)
					{
						if (y + tmpy >= iHeight - 2)
						{
							break;
						}
						iRealHeight++;

						Int dif = abs(pcEdgeAddrBg[0] + 2 * pcEdgeAddrBg[iStride] + pcEdgeAddrBg[2 * iStride] - pcEdgeAddrCur[0] - 2 * pcEdgeAddrCur[iStride] - pcEdgeAddrCur[2 * iStride]) >> 2;
						if (dif > iEdgeThreshold)
						{
							iTmpDif++;
						}
						dif = abs(pcEdgeAddrBg[-1] + 2 * pcEdgeAddrBg[-1 + iStride] + pcEdgeAddrBg[-1 + 2 * iStride] - (pcEdgeAddrBg[1] + 2 * pcEdgeAddrBg[1 + iStride] + pcEdgeAddrBg[1 + 2 * iStride])) >> 2;
						if (dif > iEdgeThreshold)
						{
							iTmpDifBg++;
						}

						pcEdgeAddrBg += iStride;
						pcEdgeAddrCur += iStride;
					}

					if (IsNormalEdge && (iTmpDif <= iRealHeight / 3 || (iTmpDif > iRealHeight / 3 && iTmpDifBg > iRealHeight / 6)))
					{
						IsNormalEdge = true;
					}
					else
					{
						IsNormalEdge = false;
					}
				}
			}
			for (Int dy = -1; dy <= 1; dy += 2)
			{
				if (iYinLCU + dy < 0 || iYinLCU + dy >= rpcPicBg->getFrameHeightInCtus())
				{
					continue;
				}

				if (rpcPicBg->getSpliceIdx((iYinLCU + dy)*rpcPicBg->getFrameWidthInCtus() + iXinLCU) != 0)
				{
					Int iTmpDif = 0;
					Int iTmpDifBg = 0;

					if (dy == -1)
					{
						pcEdgeAddrBg = &pBgLumaAddr[(y - 2)*iStride + x];
						pcEdgeAddrCur = &pCurLumaAddr[y*iStride + x];
					}

					if (dy == 1)
					{
						pcEdgeAddrCur = &pCurLumaAddr[(y + g_uiMaxCUHeight - 1)*iStride + x];
						pcEdgeAddrBg = &pBgLumaAddr[(y + g_uiMaxCUHeight + 1)*iStride + x];
					}

					Int iRealWidth = 0;
					for (Int tmpx = 0; tmpx < g_uiMaxCUWidth - 2; tmpx++)
					{
						if (x + tmpx >= iWidth - 2)
						{
							break;
						}
						iRealWidth++;

						Int dif = abs(pcEdgeAddrBg[0] + 2 * pcEdgeAddrBg[1] + pcEdgeAddrBg[2] - pcEdgeAddrCur[0] - 2 * pcEdgeAddrCur[1] - pcEdgeAddrCur[2]) >> 2;
						if (dif > iEdgeThreshold)
						{
							iTmpDif++;
						}

						dif = abs(pcEdgeAddrBg[-iStride] + 2 * pcEdgeAddrBg[1 - iStride] + pcEdgeAddrBg[2 - iStride] - (pcEdgeAddrBg[iStride] + 2 * pcEdgeAddrBg[1 + iStride] + pcEdgeAddrBg[2 + iStride])) >> 2;
						if (dif > iEdgeThreshold)
						{
							iTmpDifBg++;
						}

						pcEdgeAddrCur++;
						pcEdgeAddrBg++;
					}

					if (IsNormalEdge && (iTmpDif <= iRealWidth / 3 || (iTmpDif > iRealWidth / 3 && iTmpDifBg > iRealWidth / 6)))
					{
						IsNormalEdge = true;
					}
					else
					{
						IsNormalEdge = false;
					}
				}
			}

			if (IsNormalEdge && dLCUdist / iRealPixelCnt < 3.5 &&  dLCUdistcb / iRealPixelCnt < 0.5 &&  dLCUdistcr / iRealPixelCnt < 0.5 && rpcPicBg->getSpliceIdx(iLCUIdx) == 0)
			{
				dLCUCost = dLCUdist / iRealPixelCnt + dLCUdistcb / iRealPixelCnt + dLCUdistcr / iRealPixelCnt;

				//obtain the iMaxReplace smallest cost
				//1) find the largest cost in the iMaxReplace candidates
				for (Int i = 0; i < iMaxReplace - 1; i++)
				{
					if (pcMinCtuCost[i].dCost > pcMinCtuCost[i + 1].dCost)
					{
						swap(pcMinCtuCost[i].dCost, pcMinCtuCost[i + 1].dCost);
						swap(pcMinCtuCost[i].iCtuIdx, pcMinCtuCost[i + 1].iCtuIdx);
					}
				}
				// 2) compare the current cost with the largest cost
				if (dLCUCost < pcMinCtuCost[iMaxReplace - 1].dCost)
				{
					pcMinCtuCost[iMaxReplace - 1].dCost = dLCUCost;
					pcMinCtuCost[iMaxReplace - 1].iCtuIdx = iLCUIdx;
				}
			}
		}
	}

	// modify QP for background CTU
	if (m_pcCfg->getUseBBRwithDQP())
	{
		Int diffQP = Clip3(0, 15, (Int)(3.5*log(Double(m_pcCfg->getFramesToBeEncoded() - pcSlice->getPOC())) - 8.5));
		Int iNewQP = Clip3(0, m_pcCfg->getQP(), m_pcCfg->getQP() - diffQP);
		Int* QPList = new Int[rpcPicBg->getNumberOfCtusInFrame()];
		for (Int i = 0; i < rpcPicBg->getNumberOfCtusInFrame(); i++)
		{
			QPList[i] = pcSlice->getSliceQp();
		}

		for (Int i = 0; i < iMaxReplace; i++)
		{
			if (pcMinCtuCost[i].iCtuIdx != -1)
			{
				rpcPicBg->setSpliceIdx(pcMinCtuCost[i].iCtuIdx, pocCurr);
				QPList[pcMinCtuCost[i].iCtuIdx] = iNewQP;
			}
		}
		pcSlice->setQPList(QPList);
		pcSlice->setSliceBgQp(iNewQP);
	}
	else
	{
		for (Int i = 0; i < iMaxReplace; i++)
		{
			if (pcMinCtuCost[i].iCtuIdx != -1)
			{
				rpcPicBg->setSpliceIdx(pcMinCtuCost[i].iCtuIdx, pocCurr);
			}
		}
	}
}
#endif

#if JVET_D0186_PRECISEPSNR
Void  TEncGOP::printOutPreciseSummary(const Char * filename, Bool isField, const Bool trueBitdepthPSNR, const BitDepths &bitDepths)
{
  const Int rateMultiplier = (isField ? 2 : 1);
  m_gcAnalyzeAll.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
  m_gcAnalyzeNoFirst.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());

  FILE *fp = fopen(filename, "w");
  fprintf(fp, "\n\nPrecise All --------------------------------------------------------\n");

  const ChromaFormat chFmt = m_pcCfg->getChromaFormatIdc();

#if JVET_D0134_PSNR
#define MAYBE_trueBitdepthPSNR , trueBitdepthPSNR
#else
#define MAYBE_trueBitdepthPSNR
#endif

  m_gcAnalyzeAll.printOutPrecise(fp, 'x', chFmt, false, false MAYBE_trueBitdepthPSNR, bitDepths);
  fprintf(fp, "\n");

  fprintf(fp, "\n\nPrecise Excluding First --------------------------------------------------------\n");
  m_gcAnalyzeNoFirst.printOutPrecise(fp, 'x', chFmt, false, false MAYBE_trueBitdepthPSNR, bitDepths);
  fprintf(fp, "\n");
}
#endif

#if RRC_NS
Void TEncGOP::filterNs(TComPic * pcPic, TComSlice * pcSlice)
{
	const unsigned int channels = 1;
	TComPicYuv* rec = pcPic->getPicYuvRec();
	const unsigned int width = rec->getWidth(COMPONENT_Y);
	const unsigned int height = rec->getHeight(COMPONENT_Y);
	const unsigned int stride = rec->getStride(COMPONENT_Y);
	const int bit_depth = pcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);

	const unsigned int width_chroma = rec->getWidth(COMPONENT_Cb);
	const unsigned int height_chroma = rec->getHeight(COMPONENT_Cb);
	const unsigned int stride_chroma = rec->getStride(COMPONENT_Cb);
	const int bit_depth_chroma = pcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_CHROMA);

	const unsigned int vec_margin = 0;
	const unsigned int vec_stride = width + 2 * vec_margin;
	const unsigned int vec_height = height + 2 * vec_margin;

	const unsigned int vec_margin_chroma = 0;
	const unsigned int vec_stride_chroma = width_chroma + 2 * vec_margin_chroma;
	const unsigned int vec_height_chroma = height_chroma + 2 * vec_margin_chroma;

	//  NsImage src_ns(vec_stride, vec_height);
	NsImage src_ns(vec_stride, vec_height, vec_stride_chroma, vec_height_chroma);

	float sigma = powf(2.0f, (1.0f + 0.141f*(pcSlice->getSliceQp() - 27))); // luma

	if (pcSlice->getSliceType() != I_SLICE)
		sigma *= 1.1f * 0.5f;

	sigma *= (1 << (bit_depth - 8));

	Pel* dst = &src_ns.m_Y[0] + vec_margin + vec_margin * vec_stride;
	for (int i = 0; i < height; ++i, dst += vec_stride)
		copy(&rec->getAddr(COMPONENT_Y)[i * stride], &rec->getAddr(COMPONENT_Y)[i * stride] + width, dst);

	dst = &src_ns.m_U[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
	for (int i = 0; i < height_chroma; ++i, dst += vec_stride_chroma)
		copy(&rec->getAddr(COMPONENT_Cr)[i * stride_chroma], &rec->getAddr(COMPONENT_Cr)[i * stride_chroma] + width_chroma, dst);

	dst = &src_ns.m_V[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
	for (int i = 0; i < height_chroma; ++i, dst += vec_stride_chroma)
		copy(&rec->getAddr(COMPONENT_Cb)[i * stride_chroma], &rec->getAddr(COMPONENT_Cb)[i * stride_chroma] + width_chroma, dst);

	NsImage dst2_ns(src_ns);

	run_ns(sigma, src_ns, dst2_ns, bit_depth);


	//  for (NsImage::iterator it = dst2_ns.begin(); it != dst2_ns.end(); ++it)
	//    *it = (RRC_NS_PEL_TYPE) Clip3(0.f, pow(2.f, bit_depth) - 1, floor(*it + 0.5f));

	enum NsFilteringType { RECONSTRUCTED, FILTERED, ADAPTIVE_DEPTH_0, ADAPTIVE_DEPTH_1, ADAPTIVE_DEPTH_2, ADAPTIVE_DEPTH_3, COUNT };
	long long distortions[COUNT] = { 0 };
	int rates[COUNT] = { 0 };
	vector<vector<int> > ns_application_maps[COUNT - 2];

	TComPicYuv* org = pcPic->getPicYuvOrg();
	const unsigned int org_stride = rec->getStride(COMPONENT_Y);

	for (int i = 0; i < COUNT; ++i)
	{
		TComBitCounter tempBitCounter;
		tempBitCounter.resetBits();
		m_pcEncTop->getRDGoOnSbacCoder()->setBitstream(&tempBitCounter);
		m_pcEntropyCoder->setEntropyCoder(m_pcEncTop->getRDGoOnSbacCoder());
		m_pcEntropyCoder->resetEntropy(pcSlice);
		m_pcEncTop->getRDGoOnSbacCoder()->resetBits();

#if JVET_C0024_QTBT
		const int block_width = i <= FILTERED ? width : pcPic->getSlice(0)->getSPS()->getCTUSize() >> (i - 2);
		const int block_height = i <= FILTERED ? height : pcPic->getSlice(0)->getSPS()->getCTUSize() >> (i - 2);
#else
		const int block_width = i <= FILTERED ? width : pcPic->getSlice(0)->getSPS()->getMaxCUWidth() >> (i - 2);
		const int block_height = i <= FILTERED ? height : pcPic->getSlice(0)->getSPS()->getMaxCUHeight() >> (i - 2);
#endif
		const int width_in_blocks = width / block_width + (width % block_width != 0);
		const int height_in_blocks = height / block_height + (height % block_height != 0);

		if (i > FILTERED)
		{
			m_pcEntropyCoder->encodeNsBlockSize(i - 2);
			ns_application_maps[i - 2].resize(height_in_blocks);
			for (vector<vector<int> >::iterator it = ns_application_maps[i - 2].begin(); it != ns_application_maps[i - 2].end(); ++it)
				it->resize(width_in_blocks);
		}

		for (int block_index_y = 0; block_index_y < height_in_blocks; ++block_index_y)
		{
			const int block_pos_y = block_index_y * block_height;
			const int block_size_y = block_pos_y + block_height <= height ? block_height : height - block_pos_y;
			for (int block_index_x = 0; block_index_x < width_in_blocks; ++block_index_x)
			{
				const int block_pos_x = block_index_x * block_width;
				const int block_size_x = block_pos_x + block_width <= width ? block_width : width - block_pos_x;

				long long block_rec_mse = 0;
				long long block_filt_mse = 0;
				for (int y = block_pos_y; y < block_pos_y + block_size_y; ++y)
				{
					for (int x = block_pos_x; x < block_pos_x + block_size_x; ++x)
					{
						const int org_pix = static_cast<int>(org->getAddr(COMPONENT_Y)[y * org_stride + x]);
						const int rec_pix = static_cast<int>(rec->getAddr(COMPONENT_Y)[y * stride + x]);
						const int filt_pix = static_cast<int>(dst2_ns.m_Y[(vec_margin + x) + (vec_margin + y) * vec_stride]);

						block_rec_mse += (org_pix - rec_pix) * (org_pix - rec_pix);
						block_filt_mse += (org_pix - filt_pix) * (org_pix - filt_pix);
					}
				}

				long long block_mse;

				switch (i)
				{
				case RECONSTRUCTED:
					block_mse = block_rec_mse;
					m_pcEntropyCoder->encodeNsFWholeFrameFlag(false);
					break;
				case FILTERED:
					block_mse = block_filt_mse;
					m_pcEntropyCoder->encodeNsFWholeFrameFlag(true);
					break;
				default:
					block_mse = min(block_rec_mse, block_filt_mse);
					bool ns_application_flag = block_filt_mse < block_rec_mse;
					m_pcEntropyCoder->encodeNsFlag(ns_application_flag);
					ns_application_maps[i - 2][block_index_y][block_index_x] = ns_application_flag;
					assert(i < COUNT);
					break;
				}

				const int bit_increment = bit_depth - 8;
				const UInt shift = bit_increment << 1;

				distortions[i] += block_mse >> shift;
			}
		}
		rates[i] = m_pcEntropyCoder->getNumberOfWrittenBits();
	}

	double costs[COUNT];
	for (int i = 0; i < COUNT; ++i)
	{
		const double lambda = pcPic->getSlice(0)->getLambdas()[0];
		costs[i] = distortions[i] + lambda * rates[i];
	}

#if RRC_NS_FORCE_FILTERING_WHOLE_FRAME
	const NsFilteringType min_cost_pos = NsFilteringType::FILTERED;
#else
	const NsFilteringType min_cost_pos = static_cast<NsFilteringType>(distance(costs, min_element(costs, costs + COUNT)));
#endif

	if (min_cost_pos <= FILTERED)
	{
		pcPic->getPicSym()->setAdaptiveNs(false);

		if (min_cost_pos == FILTERED)
		{
			Pel* src = &dst2_ns.m_Y[0] + vec_margin + vec_margin * vec_stride;
			for (int i = 0; i < height; ++i, src += vec_stride)
				copy(src, src + width, &rec->getAddr(COMPONENT_Y)[i * stride]);

			src = &dst2_ns.m_U[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
			for (int i = 0; i < height_chroma; ++i, src += vec_stride_chroma)
				copy(src, src + width_chroma, &rec->getAddr(COMPONENT_Cr)[i * stride_chroma]);

			src = &dst2_ns.m_V[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
			for (int i = 0; i < height_chroma; ++i, src += vec_stride_chroma)
				copy(src, src + width_chroma, &rec->getAddr(COMPONENT_Cb)[i * stride_chroma]);


			pcPic->getPicSym()->setNsWholeFrame(true);
		}
		else
		{
			pcPic->getPicSym()->setNsWholeFrame(false);
		}
	}
	else
	{
		pcPic->getPicSym()->setAdaptiveNs(true);
		pcPic->getPicSym()->setNsBlockSize(min_cost_pos - 2);

#if JVET_C0024_QTBT
		const int block_width = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (min_cost_pos - 2);
		const int block_height = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (min_cost_pos - 2);
#else
		const int block_width = pcPic->getSlice(0)->getSPS()->getMaxCUWidth() >> (min_cost_pos - 2);
		const int block_height = pcPic->getSlice(0)->getSPS()->getMaxCUHeight() >> (min_cost_pos - 2);
#endif
		const int width_in_blocks = width / block_width + (width % block_width != 0);
		const int height_in_blocks = height / block_height + (height % block_height != 0);

		const int block_width_chroma = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (min_cost_pos - 2 + 1);
		const int block_height_chroma = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (min_cost_pos - 2 + 1);

		const int width_in_blocks_chroma = width_chroma / block_width_chroma + (width_chroma % block_width_chroma != 0);
		const int height_in_blocks_chroma = height_chroma / block_height_chroma + (height_chroma % block_height_chroma != 0);

		UInt* nsFlags = pcPic->getPicSym()->getNsFlag();
		for (int block_index_y = 0; block_index_y < height_in_blocks; ++block_index_y)
		{
			const int block_pos_y = block_index_y * block_height;
			const int block_size_y = block_pos_y + block_height <= height ? block_height : height - block_pos_y;

			for (int block_index_x = 0; block_index_x < width_in_blocks; ++block_index_x)
			{
				const int block_pos_x = block_index_x * block_width;
				const int block_size_x = block_pos_x + block_width <= width ? block_width : width - block_pos_x;

				nsFlags[block_index_y * width_in_blocks + block_index_x] = ns_application_maps[min_cost_pos - 2][block_index_y][block_index_x];
				if (ns_application_maps[min_cost_pos - 2][block_index_y][block_index_x])
				{
					Pel* src = &dst2_ns.m_Y[0] + vec_margin + vec_margin * vec_stride;
					src += block_pos_x + block_pos_y * vec_stride;
					for (int y = block_pos_y; y < block_pos_y + block_size_y; ++y, src += vec_stride)
						copy(src, src + block_size_x, &rec->getAddr(COMPONENT_Y)[y * stride + block_pos_x]);
				}
			}
		}

		for (int block_index_y_chroma = 0; block_index_y_chroma < height_in_blocks_chroma; ++block_index_y_chroma)
		{
			const int block_pos_y_chroma = block_index_y_chroma * block_height_chroma;
			const int block_size_y_chroma = block_pos_y_chroma + block_height_chroma <= height_chroma ? block_height_chroma : height_chroma - block_pos_y_chroma;

			for (int block_index_x_chroma = 0; block_index_x_chroma < width_in_blocks_chroma; ++block_index_x_chroma)
			{
				const int block_pos_x_chroma = block_index_x_chroma * block_width_chroma;
				const int block_size_x_chroma = block_pos_x_chroma + block_width_chroma <= width_chroma ? block_width_chroma : width_chroma - block_pos_x_chroma;

				if (ns_application_maps[min_cost_pos - 2][block_index_y_chroma][block_index_x_chroma])
				{
					Pel* src = &dst2_ns.m_U[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
					src += block_pos_x_chroma + block_pos_y_chroma * vec_stride_chroma;
					for (int y = block_pos_y_chroma; y < block_pos_y_chroma + block_size_y_chroma; ++y, src += vec_stride_chroma)
						copy(src, src + block_size_x_chroma, &rec->getAddr(COMPONENT_Cr)[y * stride_chroma + block_pos_x_chroma]);

					src = &dst2_ns.m_V[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
					src += block_pos_x_chroma + block_pos_y_chroma * vec_stride_chroma;
					for (int y = block_pos_y_chroma; y < block_pos_y_chroma + block_size_y_chroma; ++y, src += vec_stride_chroma)
						copy(src, src + block_size_x_chroma, &rec->getAddr(COMPONENT_Cb)[y * stride_chroma + block_pos_x_chroma]);
				}
			}
		}
	}
}
#endif

#if JVET_D0134_PSNR
#if JVET_F0064_MSSSIM
Void TEncGOP::printOutSummary(UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool printMSSSIM, const Bool trueBitdepthPSNR, const BitDepths &bitDepths)
#else
Void TEncGOP::printOutSummary(UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool trueBitdepthPSNR, const BitDepths &bitDepths)
#endif
#else
#if JVET_F0064_MSSSIM
Void TEncGOP::printOutSummary(UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool printMSSSIM, const BitDepths &bitDepths)
#else
Void TEncGOP::printOutSummary(UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const BitDepths &bitDepths)
#endif
#endif
{
  assert(uiNumAllPicCoded == m_gcAnalyzeAll.getNumPic());

  //--CFG_KDY
  const Int rateMultiplier = (isField ? 2 : 1);
  m_gcAnalyzeAll.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
  m_gcAnalyzeI.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
  m_gcAnalyzeP.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
  m_gcAnalyzeB.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    m_gcAnalyzeWPSNR.setFrmRate(m_pcCfg->getFrameRate()*rateMultiplier / (Double)m_pcCfg->getTemporalSubsampleRatio());
  }
#endif
  const ChromaFormat chFmt = m_pcCfg->getChromaFormatIdc();

#if JVET_D0134_PSNR
  //-- all
#if JVET_F0064_MSSSIM
  printf("\n\nSUMMARY --------------------------------------------------------\n");
  m_gcAnalyzeAll.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, trueBitdepthPSNR, bitDepths);

  printf("\n\nI Slices--------------------------------------------------------\n");
  m_gcAnalyzeI.printOut('i', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, trueBitdepthPSNR, bitDepths);

  printf("\n\nP Slices--------------------------------------------------------\n");
  m_gcAnalyzeP.printOut('p', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, trueBitdepthPSNR, bitDepths);

  printf("\n\nB Slices--------------------------------------------------------\n");
  m_gcAnalyzeB.printOut('b', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, trueBitdepthPSNR, bitDepths);
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    printf("\n\nWPSNR SUMMARY---------------------------------------------------\n");
    m_gcAnalyzeWPSNR.printOut('w', chFmt, printMSEBasedSNR, false, false, trueBitdepthPSNR, bitDepths);
  }
#endif
#else
  printf("\n\nSUMMARY --------------------------------------------------------\n");
  m_gcAnalyzeAll.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, bitDepths);

  printf("\n\nI Slices--------------------------------------------------------\n");
  m_gcAnalyzeI.printOut('i', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, bitDepths);

  printf("\n\nP Slices--------------------------------------------------------\n");
  m_gcAnalyzeP.printOut('p', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, bitDepths);

  printf("\n\nB Slices--------------------------------------------------------\n");
  m_gcAnalyzeB.printOut('b', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, bitDepths);
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    printf("\n\nWPSNR SUMMARY---------------------------------------------------\n");
    m_gcAnalyzeWPSNR.printOut('w', chFmt, printMSEBasedSNR, false, trueBitdepthPSNR, bitDepths);
  }
#endif
#endif


  if (!m_pcCfg->getSummaryOutFilename().empty())
  {
    m_gcAnalyzeAll.printSummary(chFmt, printSequenceMSE, trueBitdepthPSNR, bitDepths, m_pcCfg->getSummaryOutFilename());
  }

  if (!m_pcCfg->getSummaryPicFilenameBase().empty())
  {
    m_gcAnalyzeI.printSummary(chFmt, printSequenceMSE, trueBitdepthPSNR, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "I.txt");
    m_gcAnalyzeP.printSummary(chFmt, printSequenceMSE, trueBitdepthPSNR, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "P.txt");
    m_gcAnalyzeB.printSummary(chFmt, printSequenceMSE, trueBitdepthPSNR, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "B.txt");
  }

  if (isField)
  {
    //-- interlaced summary
    m_gcAnalyzeAll_in.setFrmRate(m_pcCfg->getFrameRate() / (Double)m_pcCfg->getTemporalSubsampleRatio());
    m_gcAnalyzeAll_in.setBits(m_gcAnalyzeAll.getBits());
    // prior to the above statement, the interlace analyser does not contain the correct total number of bits.

    printf("\n\nSUMMARY INTERLACED ---------------------------------------------\n");
#if JVET_F0064_MSSSIM
    m_gcAnalyzeAll_in.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, printMSSSIM, bitDepths);
#else
    m_gcAnalyzeAll_in.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, trueBitdepthPSNR, bitDepths);
#endif

    if (!m_pcCfg->getSummaryOutFilename().empty())
    {
      m_gcAnalyzeAll_in.printSummary(chFmt, printSequenceMSE, trueBitdepthPSNR, bitDepths, m_pcCfg->getSummaryOutFilename());
    }
  }
#else
  //-- all
#if JVET_F0064_MSSSIM
  printf("\n\nSUMMARY --------------------------------------------------------\n");
  m_gcAnalyzeAll.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, bitDepths);

  printf("\n\nI Slices--------------------------------------------------------\n");
  m_gcAnalyzeI.printOut('i', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, bitDepths);

  printf("\n\nP Slices--------------------------------------------------------\n");
  m_gcAnalyzeP.printOut('p', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, bitDepths);

  printf("\n\nB Slices--------------------------------------------------------\n");
  m_gcAnalyzeB.printOut('b', chFmt, printMSEBasedSNR, printSequenceMSE, printMSSSIM, bitDepths);
#else
  printf("\n\nSUMMARY --------------------------------------------------------\n");
  m_gcAnalyzeAll.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, bitDepths);

  printf("\n\nI Slices--------------------------------------------------------\n");
  m_gcAnalyzeI.printOut('i', chFmt, printMSEBasedSNR, printSequenceMSE, bitDepths);

  printf("\n\nP Slices--------------------------------------------------------\n");
  m_gcAnalyzeP.printOut('p', chFmt, printMSEBasedSNR, printSequenceMSE, bitDepths);

  printf("\n\nB Slices--------------------------------------------------------\n");
  m_gcAnalyzeB.printOut('b', chFmt, printMSEBasedSNR, printSequenceMSE, bitDepths);
#endif
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    printf("\n\nWPSNR SUMMARY---------------------------------------------------\n");
    m_gcAnalyzeWPSNR.printOut('w', chFmt, printMSEBasedSNR, false, bitDepths);
  }
#endif

  if (!m_pcCfg->getSummaryOutFilename().empty())
  {
    m_gcAnalyzeAll.printSummary(chFmt, printSequenceMSE, bitDepths, m_pcCfg->getSummaryOutFilename());
  }

  if (!m_pcCfg->getSummaryPicFilenameBase().empty())
  {
    m_gcAnalyzeI.printSummary(chFmt, printSequenceMSE, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "I.txt");
    m_gcAnalyzeP.printSummary(chFmt, printSequenceMSE, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "P.txt");
    m_gcAnalyzeB.printSummary(chFmt, printSequenceMSE, bitDepths, m_pcCfg->getSummaryPicFilenameBase() + "B.txt");
  }

  if (isField)
  {
    //-- interlaced summary
    m_gcAnalyzeAll_in.setFrmRate(m_pcCfg->getFrameRate() / (Double)m_pcCfg->getTemporalSubsampleRatio());
    m_gcAnalyzeAll_in.setBits(m_gcAnalyzeAll.getBits());
    // prior to the above statement, the interlace analyser does not contain the correct total number of bits.

    printf("\n\nSUMMARY INTERLACED ---------------------------------------------\n");
    m_gcAnalyzeAll_in.printOut('a', chFmt, printMSEBasedSNR, printSequenceMSE, bitDepths);

    if (!m_pcCfg->getSummaryOutFilename().empty())
    {
      m_gcAnalyzeAll_in.printSummary(chFmt, printSequenceMSE, bitDepths, m_pcCfg->getSummaryOutFilename());
    }
  }
#endif

  printf("\nRVM: %.3lf\n", xCalculateRVM());
}

#if !JVET_C0038_GALF
Void TEncGOP::preLoopFilterPicAll( TComPic* pcPic, UInt64& ruiDist )
{
  Bool bCalcDist = false;
  m_pcLoopFilter->setCfg(m_pcCfg->getLFCrossTileBoundaryFlag());
  m_pcLoopFilter->loopFilterPic( pcPic );

#if ALF_HM3_REFACTOR
  // Adaptive Loop filter
  if( pcPic->getSlice(0)->getSPS()->getUseALF() )
  {
    ALFParam cAlfParam;
    UInt64 dummy = 0;
    m_pcAdaptiveLoopFilter->setNumCUsInFrame(pcPic);
    m_pcAdaptiveLoopFilter->allocALFParam(&cAlfParam);
    m_pcAdaptiveLoopFilter->startALFEnc(pcPic, m_pcEntropyCoder);
    UInt uiMaxAlfCtrlDepth;
    m_pcAdaptiveLoopFilter->ALFProcess(&cAlfParam, pcPic->getSlice(0)->getLambdas()[0], pcPic->getSlice(0)->getLambdas()[1], ruiDist, dummy, uiMaxAlfCtrlDepth
#if COM16_C806_ALF_TEMPPRED_NUM
      , NULL, 0
#endif
      );
    m_pcAdaptiveLoopFilter->endALFEnc();
  }
#endif

  if (!bCalcDist)
  {
    ruiDist = xFindDistortionFrame(pcPic->getPicYuvOrg(), pcPic->getPicYuvRec(), pcPic->getPicSym()->getSPS().getBitDepths());
  }
}
#endif


// ====================================================================================================================
// Protected member functions
// ====================================================================================================================


Void TEncGOP::xInitGOP( Int iPOCLast, Int iNumPicRcvd, Bool isField )
{
  assert( iNumPicRcvd > 0 );
  //  Exception for the first frames
  if ( ( isField && (iPOCLast == 0 || iPOCLast == 1) ) || (!isField  && (iPOCLast == 0))  )
  {
    m_iGopSize    = 1;
  }
  else
  {
    m_iGopSize    = m_pcCfg->getGOPSize();
  }
  assert (m_iGopSize > 0);

  return;
}


Void TEncGOP::xGetBuffer( TComList<TComPic*>&      rcListPic,
                         TComList<TComPicYuv*>&    rcListPicYuvRecOut,
                         Int                       iNumPicRcvd,
                         Int                       iTimeOffset,
                         TComPic*&                 rpcPic,
                         TComPicYuv*&              rpcPicYuvRecOut,
                         Int                       pocCurr,
                         Bool                      isField)
{
  Int i;
  //  Rec. output
  TComList<TComPicYuv*>::iterator     iterPicYuvRec = rcListPicYuvRecOut.end();

  if (isField && pocCurr > 1 && m_iGopSize!=1)
  {
    iTimeOffset--;
  }

  for ( i = 0; i < (iNumPicRcvd - iTimeOffset + 1); i++ )
  {
    iterPicYuvRec--;
  }

  rpcPicYuvRecOut = *(iterPicYuvRec);

  //  Current pic.
  TComList<TComPic*>::iterator        iterPic       = rcListPic.begin();
  while (iterPic != rcListPic.end())
  {
    rpcPic = *(iterPic);
    rpcPic->setCurrSliceIdx(0);
    if (rpcPic->getPOC() == pocCurr)
    {
      break;
    }
    iterPic++;
  }

  assert (rpcPic != NULL);
  assert (rpcPic->getPOC() == pocCurr);

  return;
}

UInt64 TEncGOP::xFindDistortionFrame (TComPicYuv* pcPic0, TComPicYuv* pcPic1, const BitDepths &bitDepths)
{
  UInt64  uiTotalDiff = 0;

  for(Int chan=0; chan<pcPic0 ->getNumberValidComponents(); chan++)
  {
    const ComponentID ch=ComponentID(chan);
    Pel*  pSrc0   = pcPic0 ->getAddr(ch);
    Pel*  pSrc1   = pcPic1 ->getAddr(ch);
    UInt  uiShift     = 2 * DISTORTION_PRECISION_ADJUSTMENT(bitDepths.recon[toChannelType(ch)]-8);

    const Int   iStride = pcPic0->getStride(ch);
    const Int   iWidth  = pcPic0->getWidth(ch);
    const Int   iHeight = pcPic0->getHeight(ch);

    for(Int y = 0; y < iHeight; y++ )
    {
      for(Int x = 0; x < iWidth; x++ )
      {
        Intermediate_Int iTemp = pSrc0[x] - pSrc1[x];
        uiTotalDiff += UInt64((iTemp*iTemp) >> uiShift);
      }
      pSrc0 += iStride;
      pSrc1 += iStride;
    }
  }

  return uiTotalDiff;
}

#if JVET_F0064_MSSSIM
Void TEncGOP::xCalculateAddPSNRs( const Bool isField, const Bool isFieldTopFieldFirst, const Int iGOPid, TComPic* pcPic, const AccessUnit&accessUnit, TComList<TComPic*> &rcListPic, const Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE, const Bool printMSSSIM )
{
  xCalculateAddPSNR( pcPic, pcPic->getPicYuvRec(), accessUnit, dEncTime, snr_conversion, printFrameMSE, printMSSSIM );
#else
Void TEncGOP::xCalculateAddPSNRs( const Bool isField, const Bool isFieldTopFieldFirst, const Int iGOPid, TComPic* pcPic, const AccessUnit&accessUnit, TComList<TComPic*> &rcListPic, const Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE )
{
  xCalculateAddPSNR( pcPic, pcPic->getPicYuvRec(), accessUnit, dEncTime, snr_conversion, printFrameMSE);
#endif

  //In case of field coding, compute the interlaced PSNR for both fields
  if(isField)
  {
    Bool bothFieldsAreEncoded = false;
    Int correspondingFieldPOC = pcPic->getPOC();
    Int currentPicGOPPoc = m_pcCfg->getGOPEntry(iGOPid).m_POC;
    if(pcPic->getPOC() == 0)
    {
      // particular case for POC 0 and 1.
      // If they are not encoded first and separately from other pictures, we need to change this
      // POC 0 is always encoded first then POC 1 is encoded
      bothFieldsAreEncoded = false;
    }
    else if(pcPic->getPOC() == 1)
    {
      // if we are at POC 1, POC 0 has been encoded for sure
      correspondingFieldPOC = 0;
      bothFieldsAreEncoded = true;
    }
    else
    {
      if(pcPic->getPOC()%2 == 1)
      {
        correspondingFieldPOC -= 1; // all odd POC are associated with the preceding even POC (e.g poc 1 is associated to poc 0)
        currentPicGOPPoc      -= 1;
      }
      else
      {
        correspondingFieldPOC += 1; // all even POC are associated with the following odd POC (e.g poc 0 is associated to poc 1)
        currentPicGOPPoc      += 1;
      }
      for(Int i = 0; i < m_iGopSize; i ++)
      {
        if(m_pcCfg->getGOPEntry(i).m_POC == currentPicGOPPoc)
        {
          bothFieldsAreEncoded = m_pcCfg->getGOPEntry(i).m_isEncoded;
          break;
        }
      }
    }

    if(bothFieldsAreEncoded)
    {
      //get complementary top field
      TComList<TComPic*>::iterator   iterPic = rcListPic.begin();
      while ((*iterPic)->getPOC() != correspondingFieldPOC)
      {
        iterPic ++;
      }
      TComPic* correspondingFieldPic = *(iterPic);

      if( (pcPic->isTopField() && isFieldTopFieldFirst) || (!pcPic->isTopField() && !isFieldTopFieldFirst))
      {
#if JVET_F0064_MSSSIM
        xCalculateInterlacedAddPSNR(pcPic, correspondingFieldPic, pcPic->getPicYuvRec(), correspondingFieldPic->getPicYuvRec(), snr_conversion, printFrameMSE, printMSSSIM );
#else
        xCalculateInterlacedAddPSNR(pcPic, correspondingFieldPic, pcPic->getPicYuvRec(), correspondingFieldPic->getPicYuvRec(), snr_conversion, printFrameMSE);
#endif
      }
      else
      {
#if JVET_F0064_MSSSIM
        xCalculateInterlacedAddPSNR(correspondingFieldPic, pcPic, correspondingFieldPic->getPicYuvRec(), pcPic->getPicYuvRec(), snr_conversion, printFrameMSE, printMSSSIM );
#else
        xCalculateInterlacedAddPSNR(correspondingFieldPic, pcPic, correspondingFieldPic->getPicYuvRec(), pcPic->getPicYuvRec(), snr_conversion, printFrameMSE);
#endif
      }
    }
  }
}

#if JVET_F0064_MSSSIM
Void TEncGOP::xCalculateAddPSNR( TComPic* pcPic, TComPicYuv* pcPicD, const AccessUnit& accessUnit, Double dEncTime, const InputColourSpaceConversion conversion, const Bool printFrameMSE, const Bool printMSSSIM )
#else
Void TEncGOP::xCalculateAddPSNR( TComPic* pcPic, TComPicYuv* pcPicD, const AccessUnit& accessUnit, Double dEncTime, const InputColourSpaceConversion conversion, const Bool printFrameMSE)
#endif
{
  Double  dPSNR[MAX_NUM_COMPONENT];

  for(Int i=0; i<MAX_NUM_COMPONENT; i++)
  {
    dPSNR[i]=0.0;
  }

  TComPicYuv cscd;
  if (conversion!=IPCOLOURSPACE_UNCHANGED)
  {
    cscd.create(pcPicD->getWidth(COMPONENT_Y), pcPicD->getHeight(COMPONENT_Y), pcPicD->getChromaFormat(), pcPicD->getWidth(COMPONENT_Y), pcPicD->getHeight(COMPONENT_Y), 0, false);
    TVideoIOYuv::ColourSpaceConvert(*pcPicD, cscd, conversion, false);
  }
  TComPicYuv &picd=(conversion==IPCOLOURSPACE_UNCHANGED)?*pcPicD : cscd;

  //===== calculate PSNR =====
  Double MSEyuvframe[MAX_NUM_COMPONENT] = {0, 0, 0};
#if WCG_LUMA_DQP_CM_SCALE
  Double dPSNRWeighted[MAX_NUM_COMPONENT];
  Double MSEyuvframeWeighted[MAX_NUM_COMPONENT] = {0, 0, 0};
#endif
  for(Int chan=0; chan<pcPicD->getNumberValidComponents(); chan++)
  {
    const ComponentID ch=ComponentID(chan);
    const TComPicYuv *pOrgPicYuv =(conversion!=IPCOLOURSPACE_UNCHANGED) ? pcPic ->getPicYuvTrueOrg() : pcPic ->getPicYuvOrg();
    const Pel*  pOrg       = pOrgPicYuv->getAddr(ch);
    const Int   iOrgStride = pOrgPicYuv->getStride(ch);
#if WCG_LUMA_DQP_CM_SCALE
    const Pel* pOrgLuma = pOrgPicYuv->getAddr(COMPONENT_Y);
    Int  iOrgStrideLuma = pOrgPicYuv->getStride(COMPONENT_Y);
    Double dSSDtempWeighted = 0;
#endif
    Pel*  pRec             = picd.getAddr(ch);
    const Int   iRecStride = picd.getStride(ch);
    const Int   iWidth  = pcPicD->getWidth (ch) - (m_pcEncTop->getPad(0) >> pcPic->getComponentScaleX(ch));
    const Int   iHeight = pcPicD->getHeight(ch) - ((m_pcEncTop->getPad(1) >> (pcPic->isField()?1:0)) >> pcPic->getComponentScaleY(ch));

    Int   iSize   = iWidth*iHeight;

    UInt64 uiSSDtemp=0;
    for(Int y = 0; y < iHeight; y++ )
    {
      for(Int x = 0; x < iWidth; x++ )
      {
        Intermediate_Int iDiff = (Intermediate_Int)( pOrg[x] - pRec[x] );
        uiSSDtemp   += iDiff * iDiff;
#if WCG_LUMA_DQP_CM_SCALE
        if (m_pcCfg->getUseLumaDeltaQp() > 0)
        {
          Double weight = m_pcEncTop->getRdCost()->getWPSNRLumaLevelWeight(pOrgLuma[x << pcPic->getComponentScaleX(ch)]);
          dSSDtempWeighted += weight*(Double)iDiff*(Double)iDiff;
        }
#endif
      }
      pOrg += iOrgStride;
      pRec += iRecStride;
#if WCG_LUMA_DQP_CM_SCALE
      if (m_pcCfg->getUseLumaDeltaQp() > 0)
      {
        pOrgLuma += iOrgStrideLuma << pcPic->getComponentScaleY(ch);
      }
#endif
    }
    
#if JVET_D0134_PSNR
    const Int maxval = ( m_pcCfg->getTrueBitdepthPSNR() == true )?  (1 << pcPic->getPicSym()->getSPS().getBitDepth(toChannelType(ch))) - 1 : 255 << (pcPic->getPicSym()->getSPS().getBitDepth(toChannelType(ch)) - 8);
#else
    const Int maxval = 255 << (pcPic->getPicSym()->getSPS().getBitDepth(toChannelType(ch)) - 8);
#endif
    const Double fRefValue = (Double) maxval * maxval * iSize;
    dPSNR[ch]         = ( uiSSDtemp ? 10.0 * log10( fRefValue / (Double)uiSSDtemp ) : 999.99 );
    MSEyuvframe[ch]   = (Double)uiSSDtemp/(iSize);
#if WCG_LUMA_DQP_CM_SCALE
    if (m_pcCfg->getUseLumaDeltaQp() > 0)
    {
      dPSNRWeighted[ch]         = ( dSSDtempWeighted ? 10.0 * log10( fRefValue / (Double)dSSDtempWeighted ) : 999.99 );
      MSEyuvframeWeighted[ch] = dSSDtempWeighted / (iSize);
    }
#endif
  }
#if EXTENSION_360_VIDEO
  m_ext360.calculatePSNRs(pcPic);
#endif

#if JVET_F0064_MSSSIM
  //===== calculate MS-SSIM =====
  Double  MSSSIM[MAX_NUM_COMPONENT] = {0,0,0};
  if (printMSSSIM) 
  {
    for(Int chan=0; chan<pcPicD->getNumberValidComponents(); chan++)
    {
      const ComponentID ch=ComponentID(chan);
      const TComPicYuv *pOrgPicYuv =(conversion!=IPCOLOURSPACE_UNCHANGED) ? pcPic ->getPicYuvTrueOrg() : pcPic ->getPicYuvOrg();
      const Pel*  pOrg       = pOrgPicYuv->getAddr(ch);
      const Int   iOrgStride = pOrgPicYuv->getStride(ch);
      Pel*  pRec             = picd.getAddr(ch);
      const Int   iRecStride = picd.getStride(ch);
      const Int   iWidth  = pcPicD->getWidth (ch) - (m_pcEncTop->getPad(0) >> pcPic->getComponentScaleX(ch));
      const Int   iHeight = pcPicD->getHeight(ch) - ((m_pcEncTop->getPad(1) >> (pcPic->isField()?1:0)) >> pcPic->getComponentScaleY(ch));
      const UInt uiBitDepth = m_pcCfg->getBitDepth((const ChannelType)(chan > 0?1:0));

      MSSSIM[ch] = xCalculateMSSSIM (pOrg, iOrgStride, pRec, iRecStride, iWidth, iHeight, uiBitDepth);
    }
  }
#endif

  /* calculate the size of the access unit, excluding:
   *  - SEI NAL units
   */
  UInt numRBSPBytes = 0;
  for (AccessUnit::const_iterator it = accessUnit.begin(); it != accessUnit.end(); it++)
  {
    UInt numRBSPBytes_nal = UInt((*it)->m_nalUnitData.str().size());
    if (m_pcCfg->getSummaryVerboseness() > 0)
    {
      printf("*** %6s numBytesInNALunit: %u\n", nalUnitTypeToString((*it)->m_nalUnitType), numRBSPBytes_nal);
    }
    if ((*it)->m_nalUnitType != NAL_UNIT_PREFIX_SEI && (*it)->m_nalUnitType != NAL_UNIT_SUFFIX_SEI)
    {
      numRBSPBytes += numRBSPBytes_nal;
      // add start code bytes (Annex B)
      if (it == accessUnit.begin() || (*it)->m_nalUnitType == NAL_UNIT_VPS || (*it)->m_nalUnitType == NAL_UNIT_SPS || (*it)->m_nalUnitType == NAL_UNIT_PPS)
      {
        numRBSPBytes += 4;
      }
      else
      {
        numRBSPBytes += 3;
      }
    }
  }

  UInt uibits = numRBSPBytes * 8;
  m_vRVM_RP.push_back( uibits );

  //===== add PSNR =====
#if JVET_F0064_MSSSIM
  m_gcAnalyzeAll.addResult (dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    m_gcAnalyzeWPSNR.addResult(dPSNRWeighted, (Double)uibits, MSEyuvframeWeighted, MSSSIM);
  }
#endif
#else
  m_gcAnalyzeAll.addResult (dPSNR, (Double)uibits, MSEyuvframe);
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    m_gcAnalyzeWPSNR.addResult (dPSNRWeighted, (Double)uibits, MSEyuvframeWeighted);
  }
#endif
#endif

#if EXTENSION_360_VIDEO
  m_ext360.addResult(m_gcAnalyzeAll);
#endif

  TComSlice*  pcSlice = pcPic->getSlice(0);
  if (pcSlice->isIntra())
  {
#if JVET_F0064_MSSSIM
    m_gcAnalyzeI.addResult (dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#else
    m_gcAnalyzeI.addResult (dPSNR, (Double)uibits, MSEyuvframe);
#endif
#if EXTENSION_360_VIDEO
    m_ext360.addResult(m_gcAnalyzeI);
#endif
  }
  if (pcSlice->isInterP())
  {
#if JVET_F0064_MSSSIM
    m_gcAnalyzeP.addResult (dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#else
    m_gcAnalyzeP.addResult (dPSNR, (Double)uibits, MSEyuvframe);
#endif
#if EXTENSION_360_VIDEO
    m_ext360.addResult(m_gcAnalyzeP);
#endif
  }
  if (pcSlice->isInterB())
  {
#if JVET_F0064_MSSSIM
    m_gcAnalyzeB.addResult (dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#else
    m_gcAnalyzeB.addResult (dPSNR, (Double)uibits, MSEyuvframe);
#endif
#if EXTENSION_360_VIDEO
    m_ext360.addResult(m_gcAnalyzeB);
#endif
  }
#if JVET_D0186_PRECISEPSNR
  if (pcSlice->getPOC() > 0)
  {
#if JVET_F0064_MSSSIM
    m_gcAnalyzeNoFirst.addResult(dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#else
    m_gcAnalyzeNoFirst.addResult(dPSNR, (Double)uibits, MSEyuvframe);
#endif
  }
#endif

  Char c = (pcSlice->isIntra() ? 'I' : pcSlice->isInterP() ? 'P' : 'B');
  if (!pcSlice->isReferenced())
  {
    c += 32;
  }

#if ADAPTIVE_QP_SELECTION
  printf("POC %4d TId: %1d ( %c-SLICE, nQP %d QP %d ) %10d bits",
         pcSlice->getPOC(),
         pcSlice->getTLayer(),
         c,
         pcSlice->getSliceQpBase(),
         pcSlice->getSliceQp(),
         uibits );
#else
  printf("POC %4d TId: %1d ( %c-SLICE, QP %d ) %10d bits",
         pcSlice->getPOC()-pcSlice->getLastIDR(),
         pcSlice->getTLayer(),
         c,
         pcSlice->getSliceQp(),
         uibits );
#endif

  printf(" [Y %6.4lf dB    U %6.4lf dB    V %6.4lf dB]", dPSNR[COMPONENT_Y], dPSNR[COMPONENT_Cb], dPSNR[COMPONENT_Cr] );
#if JVET_F0064_MSSSIM
  if (printMSSSIM)
  {
    printf(" [MS-SSIM Y %1.6lf    U %1.6lf    V %1.6lf]", MSSSIM[COMPONENT_Y], MSSSIM[COMPONENT_Cb], MSSSIM[COMPONENT_Cr] );
  }
#endif
  if (printFrameMSE)
  {
    printf(" [Y MSE %6.4lf  U MSE %6.4lf  V MSE %6.4lf]", MSEyuvframe[COMPONENT_Y], MSEyuvframe[COMPONENT_Cb], MSEyuvframe[COMPONENT_Cr] );
  }
#if WCG_LUMA_DQP_CM_SCALE
  if (m_pcCfg->getUseLumaDeltaQp() > 0)
  {
    printf(" [WY %6.4lf dB    WU %6.4lf dB    WV %6.4lf dB]", dPSNRWeighted[COMPONENT_Y], dPSNRWeighted[COMPONENT_Cb], dPSNRWeighted[COMPONENT_Cr]);
  }
#endif
#if EXTENSION_360_VIDEO
  m_ext360.printPerPOCInfo();
#endif
  printf(" [ET %5.0f ]", dEncTime );

  for (Int iRefList = 0; iRefList < 2; iRefList++)
  {
    printf(" [L%d ", iRefList);
    for (Int iRefIndex = 0; iRefIndex < pcSlice->getNumRefIdx(RefPicList(iRefList)); iRefIndex++)
    {
      printf ("%d ", pcSlice->getRefPOC(RefPicList(iRefList), iRefIndex)-pcSlice->getLastIDR());
    }
    printf("]");
  }

  cscd.destroy();
}

#if JVET_F0064_MSSSIM
Double TEncGOP::xCalculateMSSSIM (const Pel *pOrg, const Int iOrgStride, const Pel* pRec, const Int iRecStride, const Int iWidth, const Int iHeight, const UInt uiBitDepth)
{
  Int maxScale = 5; 

  // For low resolution videos determine number of scales 
  if (iWidth < 22 || iHeight < 22)
  {
    maxScale = 1; 
  }
  else if (iWidth < 44 || iHeight < 44)
  {
    maxScale = 2; 
  }
  else if (iWidth < 88 || iHeight < 88) 
  {
    maxScale = 3; 
  }
  else if (iWidth < 176 || iHeight < 176)
  {
    maxScale = 4; 
  } 

  //Normalized gaussian mask design, 11*11, s.d. 1.5
  Double weights[11][11];
  Double coeffSum=0.0;
  for(Int y=0;y<11;y++)
  {
    for(Int x=0;x<11;x++)
    {
      weights[y][x]=exp(-((y-5.0)*(y-5.0)+(x-5.0)*(x-5.0))/4.5);
      coeffSum +=weights[y][x];
    }
  }

  for(Int y=0;y<11;y++)
  {
    for(Int x=0;x<11;x++)
    {
      weights[y][x] /=coeffSum;
    }
  }

  //Resolution based weights
  const Double exponentWeights[5][5] = {{1.0, 0, 0, 0, 0},
                                        {0.1356, 0.8644, 0, 0, 0},
                                        {0.0711, 0.4530, 0.4760, 0, 0},
                                        {0.0517, 0.3295, 0.3462, 0.2726, 0},
                                        {0.0448, 0.2856, 0.3001, 0.2363, 0.1333}};

  //Downsampling luma
  Double* original[5];
  Double* recon[5];

  Int tmpHeight = iHeight;
  Int tmpWidth = iWidth;

  for(Int i=0;i<maxScale;i++)
  {
    original[i] =new Double[tmpHeight*tmpWidth];
    recon[i] = new Double[tmpHeight*tmpWidth];
    memset(original[i],0,tmpHeight*tmpWidth*sizeof(Double));
    memset(recon[i],0,tmpHeight*tmpWidth*sizeof(Double));
    tmpHeight >>=1;
    tmpWidth >>=1;
  }

  for(Int y=0;y<iHeight;y++)
  {
    for(Int x=0;x<iWidth;x++)
    {
      original[0][y*iWidth+x] = pOrg[y*iOrgStride+x];
      recon[0][y*iWidth+x] = pRec[y*iRecStride+x];
    }
  }

  tmpHeight=iHeight>>1;
  tmpWidth=iWidth>>1;

  for(Int i=1;i<maxScale;i++)
  {
    for(Int y=0;y<tmpHeight;y++)
    {
      for(Int x=0;x<tmpWidth;x++)
      {
        original[i][y*tmpWidth+x]= (original[i-1][2*y*(2*tmpWidth)+2*x]+original[i-1][2*y*(2*tmpWidth)+2*x+1]+original[i-1][(2*y+1)*(2*tmpWidth)+2*x]+original[i-1][(2*y+1)*(2*tmpWidth)+2*x+1])/4.0;
        recon[i][y*tmpWidth+x]= (recon[i-1][2*y*(2*tmpWidth)+2*x]+recon[i-1][2*y*(2*tmpWidth)+2*x+1]+recon[i-1][(2*y+1)*(2*tmpWidth)+2*x]+recon[i-1][(2*y+1)*(2*tmpWidth)+2*x+1])/4.0;
      }
    }
    tmpHeight >>=1;
    tmpWidth >>=1;
  }

  tmpHeight = iHeight;
  tmpWidth = iWidth;
  
  // Calculate MS-SSIM
  UInt maxValue  = (1<<uiBitDepth)-1;
  const Double c1 = (0.01*maxValue)*(0.01*maxValue);
  const Double c2 = (0.03*maxValue)*(0.03*maxValue);
  
  Double orgPel;
  Double recPel;
  Double muOrg          =0.0;
  Double muRec          =0.0;
  Double muOrigSqr      =0.0;
  Double muRecSqr       =0.0;
  Double muOrigMultRec  =0.0;
  Double sigmaSqrOrig   =0.0;
  Double sigmaSqrRec    =0.0;
  Double sigmaOrigRec   =0.0;

  Double** blockSSIMluma = new Double*[maxScale];
  Double* meanSSIM = new Double[maxScale];
  tmpWidth = iWidth;
  tmpHeight = iHeight;
 
  for(Int i=0;i<maxScale;i++)
  {
    Int blocks_per_row = tmpWidth-10;
    Int blocks_per_column = tmpHeight-10;
    Int totalblocks= blocks_per_row*blocks_per_column;
    blockSSIMluma[i] =new Double[blocks_per_row*blocks_per_column];

    for(Int blockindex_Y=0; blockindex_Y<blocks_per_column; blockindex_Y++)
    {
      for(Int blockindex_X=0; blockindex_X<blocks_per_row; blockindex_X++)
      {
        for(Int y=0; y<11; y++)
        {
          for(Int x=0;x<11; x++)
          {
            orgPel =original[i][(blockindex_Y+y)*tmpWidth+(blockindex_X+x)];
            recPel=recon[i][(blockindex_Y+y)*tmpWidth+(blockindex_X+x)];

            muOrg+= orgPel*weights[y][x];
            muRec+= recPel*weights[y][x];
            muOrigSqr+=orgPel*orgPel*weights[y][x];
            muRecSqr+=recPel*recPel*weights[y][x];
            muOrigMultRec+=orgPel*recPel*weights[y][x];
          }
        }
      
        sigmaSqrOrig = muOrigSqr-(muOrg*muOrg);
        sigmaSqrRec  = muRecSqr-(muRec*muRec);
        sigmaOrigRec = muOrigMultRec-(muOrg*muRec);
        
        blockSSIMluma[i][blockindex_Y*blocks_per_row+blockindex_X]=((2.0*sigmaOrigRec + c2)/(sigmaSqrOrig+sigmaSqrRec + c2));

        if(i==(maxScale-1))
        {
          blockSSIMluma[i][blockindex_Y*blocks_per_row+blockindex_X]*=(2.0*muOrg*muRec + c1)/(muOrg*muOrg+muRec*muRec + c1);
        }

        muOrg         =0.0;
        muRec         =0.0;
        muOrigSqr     =0.0;
        muRecSqr      =0.0;
        muOrigMultRec =0.0;
      }
    }
    meanSSIM[i]= 0.0;

    for(Int index=0; index<totalblocks; index++)
    {
      meanSSIM[i] += blockSSIMluma[i][index];
    }

    meanSSIM[i] /=totalblocks;
    tmpHeight >>=1;
    tmpWidth >>=1;
  }

  Double finalMSSSIM = 1.0;

  for(Int i=0;i<maxScale;i++)
  {
    finalMSSSIM *= pow(meanSSIM[i], exponentWeights[maxScale-1][i]);
  }

  for(Int i=0;i<maxScale;i++)
  {
    delete blockSSIMluma[i];
    delete original[i];
    delete recon[i];
  }
  return finalMSSSIM;
}
#endif

Void TEncGOP::xCalculateInterlacedAddPSNR( TComPic* pcPicOrgFirstField, TComPic* pcPicOrgSecondField,
                                           TComPicYuv* pcPicRecFirstField, TComPicYuv* pcPicRecSecondField,
#if JVET_F0064_MSSSIM
                                           const InputColourSpaceConversion conversion, const Bool printFrameMSE, const Bool printMSSSIM )
#else                      
                                           const InputColourSpaceConversion conversion, const Bool printFrameMSE )                                           
#endif
{
  const TComSPS &sps=pcPicOrgFirstField->getPicSym()->getSPS();
  Double  dPSNR[MAX_NUM_COMPONENT] = {};
  TComPic    *apcPicOrgFields[2]={pcPicOrgFirstField, pcPicOrgSecondField};
  TComPicYuv *apcPicRecFields[2]={pcPicRecFirstField, pcPicRecSecondField};

  for(Int i=0; i<MAX_NUM_COMPONENT; i++)
  {
    dPSNR[i]=0.0;
  }

  TComPicYuv cscd[2 /* first/second field */];
  if (conversion!=IPCOLOURSPACE_UNCHANGED)
  {
    for(UInt fieldNum=0; fieldNum<2; fieldNum++)
    {
      TComPicYuv &reconField=*(apcPicRecFields[fieldNum]);
      cscd[fieldNum].create(reconField.getWidth(COMPONENT_Y), reconField.getHeight(COMPONENT_Y), reconField.getChromaFormat(), reconField.getWidth(COMPONENT_Y), reconField.getHeight(COMPONENT_Y), 0, false);
      TVideoIOYuv::ColourSpaceConvert(reconField, cscd[fieldNum], conversion, false);
      apcPicRecFields[fieldNum]=cscd+fieldNum;
    }
  }

  //===== calculate PSNR =====
  Double MSEyuvframe[MAX_NUM_COMPONENT] = {0, 0, 0};

  assert(apcPicRecFields[0]->getChromaFormat()==apcPicRecFields[1]->getChromaFormat());
  const UInt numValidComponents=apcPicRecFields[0]->getNumberValidComponents();

  for(Int chan=0; chan<numValidComponents; chan++)
  {
    const ComponentID ch=ComponentID(chan);
    assert(apcPicRecFields[0]->getWidth(ch)==apcPicRecFields[1]->getWidth(ch));
    assert(apcPicRecFields[0]->getHeight(ch)==apcPicRecFields[1]->getHeight(ch));

    UInt64 uiSSDtemp=0;
    const Int   iWidth  = apcPicRecFields[0]->getWidth (ch) - (m_pcEncTop->getPad(0) >> apcPicRecFields[0]->getComponentScaleX(ch));
    const Int   iHeight = apcPicRecFields[0]->getHeight(ch) - ((m_pcEncTop->getPad(1) >> 1) >> apcPicRecFields[0]->getComponentScaleY(ch));

    Int   iSize   = iWidth*iHeight;

    for(UInt fieldNum=0; fieldNum<2; fieldNum++)
    {
      TComPic *pcPic=apcPicOrgFields[fieldNum];
      TComPicYuv *pcPicD=apcPicRecFields[fieldNum];

      const Pel*  pOrg    = (conversion!=IPCOLOURSPACE_UNCHANGED) ? pcPic ->getPicYuvTrueOrg()->getAddr(ch) : pcPic ->getPicYuvOrg()->getAddr(ch);
      Pel*  pRec    = pcPicD->getAddr(ch);
      const Int   iStride = pcPicD->getStride(ch);


      for(Int y = 0; y < iHeight; y++ )
      {
        for(Int x = 0; x < iWidth; x++ )
        {
          Intermediate_Int iDiff = (Intermediate_Int)( pOrg[x] - pRec[x] );
          uiSSDtemp   += iDiff * iDiff;
        }
        pOrg += iStride;
        pRec += iStride;
      }
    }

#if JVET_D0134_PSNR
    const Int maxval = ( m_pcCfg->getTrueBitdepthPSNR() == true )?  (1 << sps.getBitDepth(toChannelType(ch))) - 1 : 255 << (sps.getBitDepth(toChannelType(ch)) - 8);
#else
    const Int maxval = 255 << (sps.getBitDepth(toChannelType(ch)) - 8);
#endif
    const Double fRefValue = (Double) maxval * maxval * iSize*2;
    dPSNR[ch]         = ( uiSSDtemp ? 10.0 * log10( fRefValue / (Double)uiSSDtemp ) : 999.99 );
    MSEyuvframe[ch]   = (Double)uiSSDtemp/(iSize*2);
  }

#if JVET_F0064_MSSSIM
  //===== calculate MS-SSIM =====
  Double  MSSSIM[MAX_NUM_COMPONENT] = {0,0,0};
  if (printMSSSIM)
  {
    for(Int chan=0; chan<numValidComponents; chan++)
    {
         const ComponentID ch=ComponentID(chan);
      assert(apcPicRecFields[0]->getWidth(ch)==apcPicRecFields[1]->getWidth(ch));
      assert(apcPicRecFields[0]->getHeight(ch)==apcPicRecFields[1]->getHeight(ch));

      Double uiMSSSIMtemp = 0.0;
      const Int   iWidth  = apcPicRecFields[0]->getWidth (ch) - (m_pcEncTop->getPad(0) >> apcPicRecFields[0]->getComponentScaleX(ch));
      const Int   iHeight = apcPicRecFields[0]->getHeight(ch) - ((m_pcEncTop->getPad(1) >> 1) >> apcPicRecFields[0]->getComponentScaleY(ch));

      for(UInt fieldNum=0; fieldNum<2; fieldNum++)
      {
        TComPic *pcPic=apcPicOrgFields[fieldNum];
        TComPicYuv *pcPicD=apcPicRecFields[fieldNum];

        const Pel*  pOrg    = (conversion!=IPCOLOURSPACE_UNCHANGED) ? pcPic ->getPicYuvTrueOrg()->getAddr(ch) : pcPic ->getPicYuvOrg()->getAddr(ch);
        Pel*  pRec    = pcPicD->getAddr(ch);
        const Int   iStride = pcPicD->getStride(ch);
        const UInt uiBitDepth = m_pcCfg->getBitDepth((const ChannelType)(chan > 0?1:0));

        uiMSSSIMtemp += xCalculateMSSSIM (pOrg, iStride, pRec, iStride, iWidth, iHeight, uiBitDepth);
      }

      MSSSIM[ch] = uiMSSSIMtemp/2;
    }
  }
#endif

  UInt uibits = 0; // the number of bits for the pair is not calculated here - instead the overall total is used elsewhere.

  //===== add PSNR =====
#if JVET_F0064_MSSSIM
  m_gcAnalyzeAll_in.addResult (dPSNR, (Double)uibits, MSEyuvframe, MSSSIM);
#else
  m_gcAnalyzeAll_in.addResult (dPSNR, (Double)uibits, MSEyuvframe);
#endif

  printf("\n                                      Interlaced frame %d: [Y %6.4lf dB    U %6.4lf dB    V %6.4lf dB]", pcPicOrgSecondField->getPOC()/2 , dPSNR[COMPONENT_Y], dPSNR[COMPONENT_Cb], dPSNR[COMPONENT_Cr] );
#if JVET_F0064_MSSSIM
  if (printMSSSIM)
  {
    printf(" [MS-SSIM Y %1.6lf    U %1.6lf    V %1.6lf]", MSSSIM[COMPONENT_Y], MSSSIM[COMPONENT_Cb], MSSSIM[COMPONENT_Cr] );
  }
#endif
  if (printFrameMSE)
  {
    printf(" [Y MSE %6.4lf  U MSE %6.4lf  V MSE %6.4lf]", MSEyuvframe[COMPONENT_Y], MSEyuvframe[COMPONENT_Cb], MSEyuvframe[COMPONENT_Cr] );
  }

  for(UInt fieldNum=0; fieldNum<2; fieldNum++)
  {
    cscd[fieldNum].destroy();
  }
}

/** Function for deciding the nal_unit_type.
 * \param pocCurr POC of the current picture
 * \param lastIDR  POC of the last IDR picture
 * \param isField  true to indicate field coding
 * \returns the NAL unit type of the picture
 * This function checks the configuration and returns the appropriate nal_unit_type for the picture.
 */
NalUnitType TEncGOP::getNalUnitType(Int pocCurr, Int lastIDR, Bool isField)
{
  if (pocCurr == 0)
  {
    return NAL_UNIT_CODED_SLICE_IDR_W_RADL;
  }

  if(m_pcCfg->getEfficientFieldIRAPEnabled() && isField && pocCurr == 1)
  {
    // to avoid the picture becoming an IRAP
    return NAL_UNIT_CODED_SLICE_TRAIL_R;
  }

  if(m_pcCfg->getDecodingRefreshType() != 3 && (pocCurr - isField) % m_pcCfg->getIntraPeriod() == 0)
  {
    if (m_pcCfg->getDecodingRefreshType() == 1)
    {
      return NAL_UNIT_CODED_SLICE_CRA;
    }
    else if (m_pcCfg->getDecodingRefreshType() == 2)
    {
      return NAL_UNIT_CODED_SLICE_IDR_W_RADL;
    }
  }
  if(m_pocCRA>0)
  {
    if(pocCurr<m_pocCRA)
    {
      // All leading pictures are being marked as TFD pictures here since current encoder uses all
      // reference pictures while encoding leading pictures. An encoder can ensure that a leading
      // picture can be still decodable when random accessing to a CRA/CRANT/BLA/BLANT picture by
      // controlling the reference pictures used for encoding that leading picture. Such a leading
      // picture need not be marked as a TFD picture.
      return NAL_UNIT_CODED_SLICE_RASL_R;
    }
  }
  if (lastIDR>0)
  {
    if (pocCurr < lastIDR)
    {
      return NAL_UNIT_CODED_SLICE_RADL_R;
    }
  }
  return NAL_UNIT_CODED_SLICE_TRAIL_R;
}

Double TEncGOP::xCalculateRVM()
{
  Double dRVM = 0;

  if( m_pcCfg->getGOPSize() == 1 && m_pcCfg->getIntraPeriod() != 1 && m_pcCfg->getFramesToBeEncoded() > RVM_VCEGAM10_M * 2 )
  {
    // calculate RVM only for lowdelay configurations
    std::vector<Double> vRL , vB;
    size_t N = m_vRVM_RP.size();
    vRL.resize( N );
    vB.resize( N );

    Int i;
    Double dRavg = 0 , dBavg = 0;
    vB[RVM_VCEGAM10_M] = 0;
    for( i = RVM_VCEGAM10_M + 1 ; i < N - RVM_VCEGAM10_M + 1 ; i++ )
    {
      vRL[i] = 0;
      for( Int j = i - RVM_VCEGAM10_M ; j <= i + RVM_VCEGAM10_M - 1 ; j++ )
      {
        vRL[i] += m_vRVM_RP[j];
      }
      vRL[i] /= ( 2 * RVM_VCEGAM10_M );
      vB[i] = vB[i-1] + m_vRVM_RP[i] - vRL[i];
      dRavg += m_vRVM_RP[i];
      dBavg += vB[i];
    }

    dRavg /= ( N - 2 * RVM_VCEGAM10_M );
    dBavg /= ( N - 2 * RVM_VCEGAM10_M );

    Double dSigamB = 0;
    for( i = RVM_VCEGAM10_M + 1 ; i < N - RVM_VCEGAM10_M + 1 ; i++ )
    {
      Double tmp = vB[i] - dBavg;
      dSigamB += tmp * tmp;
    }
    dSigamB = sqrt( dSigamB / ( N - 2 * RVM_VCEGAM10_M ) );

    Double f = sqrt( 12.0 * ( RVM_VCEGAM10_M - 1 ) / ( RVM_VCEGAM10_M + 1 ) );

    dRVM = dSigamB / dRavg * f;
  }

  return( dRVM );
}

/** Attaches the input bitstream to the stream in the output NAL unit
    Updates rNalu to contain concatenated bitstream. rpcBitstreamRedirect is cleared at the end of this function call.
 *  \param codedSliceData contains the coded slice data (bitstream) to be concatenated to rNalu
 *  \param rNalu          target NAL unit
 */
Void TEncGOP::xAttachSliceDataToNalUnit (OutputNALUnit& rNalu, TComOutputBitstream* codedSliceData)
{
  // Byte-align
  rNalu.m_Bitstream.writeByteAlignment();   // Slice header byte-alignment

  // Perform bitstream concatenation
  if (codedSliceData->getNumberOfWrittenBits() > 0)
  {
    rNalu.m_Bitstream.addSubstream(codedSliceData);
  }

  m_pcEntropyCoder->setBitstream(&rNalu.m_Bitstream);

  codedSliceData->clear();
}

// Function will arrange the long-term pictures in the decreasing order of poc_lsb_lt,
// and among the pictures with the same lsb, it arranges them in increasing delta_poc_msb_cycle_lt value
Void TEncGOP::arrangeLongtermPicturesInRPS(TComSlice *pcSlice, TComList<TComPic*>& rcListPic)
{
  if(pcSlice->getRPS()->getNumberOfLongtermPictures() == 0)
  {
    return;
  }
  // we can only modify the local RPS!
  assert (pcSlice->getRPSidx()==-1);
  TComReferencePictureSet *rps = pcSlice->getLocalRPS();

  // Arrange long-term reference pictures in the correct order of LSB and MSB,
  // and assign values for pocLSBLT and MSB present flag
  Int longtermPicsPoc[MAX_NUM_REF_PICS], longtermPicsLSB[MAX_NUM_REF_PICS], indices[MAX_NUM_REF_PICS];
  Int longtermPicsMSB[MAX_NUM_REF_PICS];
  Bool mSBPresentFlag[MAX_NUM_REF_PICS];
  ::memset(longtermPicsPoc, 0, sizeof(longtermPicsPoc));    // Store POC values of LTRP
  ::memset(longtermPicsLSB, 0, sizeof(longtermPicsLSB));    // Store POC LSB values of LTRP
  ::memset(longtermPicsMSB, 0, sizeof(longtermPicsMSB));    // Store POC LSB values of LTRP
  ::memset(indices        , 0, sizeof(indices));            // Indices to aid in tracking sorted LTRPs
  ::memset(mSBPresentFlag , 0, sizeof(mSBPresentFlag));     // Indicate if MSB needs to be present

  // Get the long-term reference pictures
  Int offset = rps->getNumberOfNegativePictures() + rps->getNumberOfPositivePictures();
  Int i, ctr = 0;
  Int maxPicOrderCntLSB = 1 << pcSlice->getSPS()->getBitsForPOC();
  for(i = rps->getNumberOfPictures() - 1; i >= offset; i--, ctr++)
  {
    longtermPicsPoc[ctr] = rps->getPOC(i);                                  // LTRP POC
    longtermPicsLSB[ctr] = getLSB(longtermPicsPoc[ctr], maxPicOrderCntLSB); // LTRP POC LSB
    indices[ctr]      = i;
    longtermPicsMSB[ctr] = longtermPicsPoc[ctr] - longtermPicsLSB[ctr];
  }
  Int numLongPics = rps->getNumberOfLongtermPictures();
  assert(ctr == numLongPics);

  // Arrange pictures in decreasing order of MSB;
  for(i = 0; i < numLongPics; i++)
  {
    for(Int j = 0; j < numLongPics - 1; j++)
    {
      if(longtermPicsMSB[j] < longtermPicsMSB[j+1])
      {
        std::swap(longtermPicsPoc[j], longtermPicsPoc[j+1]);
        std::swap(longtermPicsLSB[j], longtermPicsLSB[j+1]);
        std::swap(longtermPicsMSB[j], longtermPicsMSB[j+1]);
        std::swap(indices[j]        , indices[j+1]        );
      }
    }
  }

  for(i = 0; i < numLongPics; i++)
  {
    // Check if MSB present flag should be enabled.
    // Check if the buffer contains any pictures that have the same LSB.
    TComList<TComPic*>::iterator  iterPic = rcListPic.begin();
    TComPic*                      pcPic;
    while ( iterPic != rcListPic.end() )
    {
      pcPic = *iterPic;
      if( (getLSB(pcPic->getPOC(), maxPicOrderCntLSB) == longtermPicsLSB[i])   &&     // Same LSB
                                      (pcPic->getSlice(0)->isReferenced())     &&    // Reference picture
                                        (pcPic->getPOC() != longtermPicsPoc[i])    )  // Not the LTRP itself
      {
        mSBPresentFlag[i] = true;
        break;
      }
      iterPic++;
    }
  }

  // tempArray for usedByCurr flag
  Bool tempArray[MAX_NUM_REF_PICS]; ::memset(tempArray, 0, sizeof(tempArray));
  for(i = 0; i < numLongPics; i++)
  {
    tempArray[i] = rps->getUsed(indices[i]);
  }
  // Now write the final values;
  ctr = 0;
  Int currMSB = 0, currLSB = 0;
  // currPicPoc = currMSB + currLSB
  currLSB = getLSB(pcSlice->getPOC(), maxPicOrderCntLSB);
  currMSB = pcSlice->getPOC() - currLSB;

  for(i = rps->getNumberOfPictures() - 1; i >= offset; i--, ctr++)
  {
    rps->setPOC                   (i, longtermPicsPoc[ctr]);
    rps->setDeltaPOC              (i, - pcSlice->getPOC() + longtermPicsPoc[ctr]);
    rps->setUsed                  (i, tempArray[ctr]);
    rps->setPocLSBLT              (i, longtermPicsLSB[ctr]);
    rps->setDeltaPocMSBCycleLT    (i, (currMSB - (longtermPicsPoc[ctr] - longtermPicsLSB[ctr])) / maxPicOrderCntLSB);
    rps->setDeltaPocMSBPresentFlag(i, mSBPresentFlag[ctr]);

    assert(rps->getDeltaPocMSBCycleLT(i) >= 0);   // Non-negative value
  }
  for(i = rps->getNumberOfPictures() - 1, ctr = 1; i >= offset; i--, ctr++)
  {
    for(Int j = rps->getNumberOfPictures() - 1 - ctr; j >= offset; j--)
    {
      // Here at the encoder we know that we have set the full POC value for the LTRPs, hence we
      // don't have to check the MSB present flag values for this constraint.
      assert( rps->getPOC(i) != rps->getPOC(j) ); // If assert fails, LTRP entry repeated in RPS!!!
    }
  }
}

Void TEncGOP::applyDeblockingFilterMetric( TComPic* pcPic, UInt uiNumSlices )
{
  TComPicYuv* pcPicYuvRec = pcPic->getPicYuvRec();
  Pel* Rec    = pcPicYuvRec->getAddr(COMPONENT_Y);
  Pel* tempRec = Rec;
  Int  stride = pcPicYuvRec->getStride(COMPONENT_Y);
  UInt log2maxTB = pcPic->getSlice(0)->getSPS()->getQuadtreeTULog2MaxSize();
  UInt maxTBsize = (1<<log2maxTB);
  const UInt minBlockArtSize = 8;
  const UInt picWidth = pcPicYuvRec->getWidth(COMPONENT_Y);
  const UInt picHeight = pcPicYuvRec->getHeight(COMPONENT_Y);
  const UInt noCol = (picWidth>>log2maxTB);
  const UInt noRows = (picHeight>>log2maxTB);
  assert(noCol > 1);
  assert(noRows > 1);
  UInt64 *colSAD = (UInt64*)malloc(noCol*sizeof(UInt64));
  UInt64 *rowSAD = (UInt64*)malloc(noRows*sizeof(UInt64));
  UInt colIdx = 0;
  UInt rowIdx = 0;
  Pel p0, p1, p2, q0, q1, q2;

  Int qp = pcPic->getSlice(0)->getSliceQp();
  const Int bitDepthLuma=pcPic->getSlice(0)->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);
  Int bitdepthScale = 1 << (bitDepthLuma-8);
  Int beta = TComLoopFilter::getBeta( qp ) * bitdepthScale;
  const Int thr2 = (beta>>2);
  const Int thr1 = 2*bitdepthScale;
  UInt a = 0;

  memset(colSAD, 0, noCol*sizeof(UInt64));
  memset(rowSAD, 0, noRows*sizeof(UInt64));

  if (maxTBsize > minBlockArtSize)
  {
    // Analyze vertical artifact edges
    for(Int c = maxTBsize; c < picWidth; c += maxTBsize)
    {
      for(Int r = 0; r < picHeight; r++)
      {
        p2 = Rec[c-3];
        p1 = Rec[c-2];
        p0 = Rec[c-1];
        q0 = Rec[c];
        q1 = Rec[c+1];
        q2 = Rec[c+2];
        a = ((abs(p2-(p1<<1)+p0)+abs(q0-(q1<<1)+q2))<<1);
        if ( thr1 < a && a < thr2)
        {
          colSAD[colIdx] += abs(p0 - q0);
        }
        Rec += stride;
      }
      colIdx++;
      Rec = tempRec;
    }

    // Analyze horizontal artifact edges
    for(Int r = maxTBsize; r < picHeight; r += maxTBsize)
    {
      for(Int c = 0; c < picWidth; c++)
      {
        p2 = Rec[c + (r-3)*stride];
        p1 = Rec[c + (r-2)*stride];
        p0 = Rec[c + (r-1)*stride];
        q0 = Rec[c + r*stride];
        q1 = Rec[c + (r+1)*stride];
        q2 = Rec[c + (r+2)*stride];
        a = ((abs(p2-(p1<<1)+p0)+abs(q0-(q1<<1)+q2))<<1);
        if (thr1 < a && a < thr2)
        {
          rowSAD[rowIdx] += abs(p0 - q0);
        }
      }
      rowIdx++;
    }
  }

  UInt64 colSADsum = 0;
  UInt64 rowSADsum = 0;
  for(Int c = 0; c < noCol-1; c++)
  {
    colSADsum += colSAD[c];
  }
  for(Int r = 0; r < noRows-1; r++)
  {
    rowSADsum += rowSAD[r];
  }

  colSADsum <<= 10;
  rowSADsum <<= 10;
  colSADsum /= (noCol-1);
  colSADsum /= picHeight;
  rowSADsum /= (noRows-1);
  rowSADsum /= picWidth;

  UInt64 avgSAD = ((colSADsum + rowSADsum)>>1);
  avgSAD >>= (bitDepthLuma-8);

  if ( avgSAD > 2048 )
  {
    avgSAD >>= 9;
    Int offset = Clip3(2,6,(Int)avgSAD);
    for (Int i=0; i<uiNumSlices; i++)
    {
      pcPic->getSlice(i)->setDeblockingFilterOverrideFlag(true);
      pcPic->getSlice(i)->setDeblockingFilterDisable(false);
      pcPic->getSlice(i)->setDeblockingFilterBetaOffsetDiv2( offset );
      pcPic->getSlice(i)->setDeblockingFilterTcOffsetDiv2( offset );
    }
  }
  else
  {
    for (Int i=0; i<uiNumSlices; i++)
    {
      pcPic->getSlice(i)->setDeblockingFilterOverrideFlag(false);
      pcPic->getSlice(i)->setDeblockingFilterDisable(        pcPic->getSlice(i)->getPPS()->getPicDisableDeblockingFilterFlag() );
      pcPic->getSlice(i)->setDeblockingFilterBetaOffsetDiv2( pcPic->getSlice(i)->getPPS()->getDeblockingFilterBetaOffsetDiv2() );
      pcPic->getSlice(i)->setDeblockingFilterTcOffsetDiv2(   pcPic->getSlice(i)->getPPS()->getDeblockingFilterTcOffsetDiv2()   );
    }
  }

  free(colSAD);
  free(rowSAD);
}

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
Int TEncGOP::xUpdateTStates (UInt uiSliceType, UInt uiSliceQP, TComStats* apcStats )
{ 
  Int iQP = 0, k;
  for (k = 0; k < NUM_QP_PROB; k++)
  {
    if (apcStats-> aaQPUsed[uiSliceType][k].used ==true && apcStats-> aaQPUsed[uiSliceType][k].QP == uiSliceQP)
    {
      iQP  = k;
      apcStats-> aaQPUsed[uiSliceType][k].firstUsed = false;
      break;
    }
    else if (apcStats-> aaQPUsed[uiSliceType][k].used ==false)
    {
      iQP = k;
      apcStats-> aaQPUsed[uiSliceType][k].used = true;
      apcStats-> aaQPUsed[uiSliceType][k].QP = uiSliceQP;
      apcStats-> aaQPUsed[uiSliceType][k].firstUsed = true;
#if VCEG_AZ07_BAC_ADAPT_WDOW
      for (Int index=0; index < MAX_NUM_CTX_MOD; index++) 
      {
        apcStats->m_uiCtxMAP[uiSliceType][iQP][index]     = 0;
        apcStats->m_uiCtxCodeIdx[uiSliceType][iQP][index] = ALPHA0;
      }
#endif
      break;
    }
  }

  iQP = -1;
  for (k = 0; k < NUM_QP_PROB; k++)
  {
    if (apcStats-> aaQPUsed[uiSliceType][k].firstUsed ==true && apcStats-> aaQPUsed[uiSliceType][k].used ==true && apcStats-> aaQPUsed[uiSliceType][k].QP == uiSliceQP)
    {
      iQP  = k;
      break;
    }
  }
  return iQP;
}
#endif

//! \}
