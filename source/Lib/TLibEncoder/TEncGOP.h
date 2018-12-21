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

/** \file     TEncGOP.h
    \brief    GOP encoder class (header)
*/

#ifndef __TENCGOP__
#define __TENCGOP__

#include <list>

#include <stdlib.h>

#include "TLibCommon/TComList.h"
#include "TLibCommon/TComPic.h"
#include "TLibCommon/TComBitCounter.h"
#include "TLibCommon/TComLoopFilter.h"
#include "TLibCommon/AccessUnit.h"
#include "TEncSampleAdaptiveOffset.h"
#include "TEncSlice.h"
#include "TEncEntropy.h"
#include "TEncCavlc.h"
#include "TEncSbac.h"
#include "SEIwrite.h"
#include "SEIEncoder.h"
#if EXTENSION_360_VIDEO
#include "TAppEncHelper360/TExt360EncGop.h"
#endif

#include "TEncAnalyze.h"
#include "TEncRateCtrl.h"
#include <vector>
#if ALF_HM3_REFACTOR
#include "TEncAdaptiveLoopFilter.h"
#endif
#if CNN_BASED_INLOOP_FILTER
#include "TEncCNNLoopFilter.h"
#endif
//! \ingroup TLibEncoder
//! \{

class TEncTop;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

class TEncGOP
{
  class DUData
  {
  public:
    DUData()
    :accumBitsDU(0)
    ,accumNalsDU(0) {};

    Int accumBitsDU;
    Int accumNalsDU;
  };

private:

  TEncAnalyze             m_gcAnalyzeAll;
#if JVET_D0186_PRECISEPSNR
  TEncAnalyze             m_gcAnalyzeNoFirst;
#endif
  TEncAnalyze             m_gcAnalyzeI;
  TEncAnalyze             m_gcAnalyzeP;
  TEncAnalyze             m_gcAnalyzeB;
#if WCG_LUMA_DQP_CM_SCALE
  TEncAnalyze             m_gcAnalyzeWPSNR;
#endif
  TEncAnalyze             m_gcAnalyzeAll_in;

#if EXTENSION_360_VIDEO
  TExt360EncGop           m_ext360;
public:
  TExt360EncGop &getExt360Data() { return m_ext360; }
private:
#endif

  //  Data
  Bool                    m_bLongtermTestPictureHasBeenCoded;
  Bool                    m_bLongtermTestPictureHasBeenCoded2;
  UInt                    m_numLongTermRefPicSPS;
  UInt                    m_ltRefPicPocLsbSps[MAX_NUM_LONG_TERM_REF_PICS];
  Bool                    m_ltRefPicUsedByCurrPicFlag[MAX_NUM_LONG_TERM_REF_PICS];
  Int                     m_iLastIDR;
  Int                     m_iGopSize;
  Int                     m_iNumPicCoded;
  Bool                    m_bFirst;
  Int                     m_iLastRecoveryPicPOC;

  //  Access channel
  TEncTop*                m_pcEncTop;
  TEncCfg*                m_pcCfg;
  TEncSlice*              m_pcSliceEncoder;
  TComList<TComPic*>*     m_pcListPic;

  TEncEntropy*            m_pcEntropyCoder;
  TEncCavlc*              m_pcCavlcCoder;
  TEncSbac*               m_pcSbacCoder;
  TEncBinCABAC*           m_pcBinCABAC;
  TComLoopFilter*         m_pcLoopFilter;

  SEIWriter               m_seiWriter;
#if CNN_BASED_INLOOP_FILTER
  TEncCNNLoopFilter*      m_pcCNNLoopFilter;
#endif
  //--Adaptive Loop filter
  TEncSampleAdaptiveOffset*  m_pcSAO;
  TEncRateCtrl*           m_pcRateCtrl;
  // indicate sequence first
  Bool                    m_bSeqFirst;

  // clean decoding refresh
  Bool                    m_bRefreshPending;
  Int                     m_pocCRA;
  NalUnitType             m_associatedIRAPType;
  Int                     m_associatedIRAPPOC;

  std::vector<Int> m_vRVM_RP;
  UInt                    m_lastBPSEI;
  UInt                    m_totalCoded;
  Bool                    m_bufferingPeriodSEIPresentInAU;
  SEIEncoder              m_seiEncoder;

#if ALF_HM3_REFACTOR
  // Adaptive Loop filter
  TEncAdaptiveLoopFilter* m_pcAdaptiveLoopFilter;
#if COM16_C806_ALF_TEMPPRED_NUM
#if JVET_E0104_ALF_TEMP_SCALABILITY
  static Int           m_iStoredAlfParaNum[JVET_E0104_ALF_MAX_TEMPLAYERID];
  ALFParam             m_acStoredAlfPara[JVET_E0104_ALF_MAX_TEMPLAYERID][COM16_C806_ALF_TEMPPRED_NUM];
#else
  static Int           m_iStoredAlfParaNum;
  ALFParam             m_acStoredAlfPara[COM16_C806_ALF_TEMPPRED_NUM];
#endif
#endif
#endif

#if RRC_NS
	Void filterNs(TComPic * pcPic, TComSlice * pcSlice);
#endif

public:
  TEncGOP();
  virtual ~TEncGOP();

  Void  create      ();
  Void  destroy     ();

#if PVC_AQS // interface
  TEncTop* getEncTop() { return m_pcEncTop; }
#endif

  Void  init        ( TEncTop* pcTEncTop );
  Void  compressGOP ( Int iPOCLast, Int iNumPicRcvd, TComList<TComPic*>& rcListPic, TComList<TComPicYuv*>& rcListPicYuvRec,
                      std::list<AccessUnit>& accessUnitsInGOP, Bool isField, Bool isTff, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE 
#if JVET_F0064_MSSSIM
                      , const Bool printMSSSIM 
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
                      ,TComStats* m_apcStats
#endif
    );

  Void  xAttachSliceDataToNalUnit (OutputNALUnit& rNalu, TComOutputBitstream* pcBitstreamRedirect);


  Int   getGOPSize()          { return  m_iGopSize;  }

  TComList<TComPic*>*   getListPic()      { return m_pcListPic; }

#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void ArrangeBackgroundReference(TComSlice* pcSlice, TComList<TComPic*>& rcListPic, Int pocCurr);
#endif

#if JVET_D0134_PSNR
#if JVET_F0064_MSSSIM
  Void  printOutSummary      ( UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool printMSSSIM, const Bool trueBitdepthPSNR, const BitDepths &bitDepths );
#else
  Void  printOutSummary      ( UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool trueBitdepthPSNR, const BitDepths &bitDepths );
#endif
#else
#if JVET_F0064_MSSSIM
  Void  printOutSummary      ( UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const Bool printMSSSIM, const BitDepths &bitDepths );
#else
  Void  printOutSummary      ( UInt uiNumAllPicCoded, Bool isField, const Bool printMSEBasedSNR, const Bool printSequenceMSE, const BitDepths &bitDepths );
#endif
#endif
#if JVET_D0186_PRECISEPSNR
  Void  printOutPreciseSummary(const Char * filename, Bool isField, const Bool trueBitdepthPSNR, const BitDepths &bitDepths);
#endif

#if !JVET_C0038_GALF
  Void  preLoopFilterPicAll  ( TComPic* pcPic, UInt64& ruiDist );
#endif
  TEncSlice*  getSliceEncoder()   { return m_pcSliceEncoder; }
  NalUnitType getNalUnitType( Int pocCurr, Int lastIdr, Bool isField );
  Void arrangeLongtermPicturesInRPS(TComSlice *, TComList<TComPic*>& );
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  Int xUpdateTStates ( UInt uiSliceType, UInt uiSliceQP,  TComStats* apcStats );  
#endif

#if EXTENSION_360_VIDEO
  TEncAnalyze& getAnalyzeAllData() { return m_gcAnalyzeAll; }
  TEncAnalyze& getAnalyzeIData()   { return m_gcAnalyzeI; }
  TEncAnalyze& getAnalyzePData()   { return m_gcAnalyzeP; }
  TEncAnalyze& getAnalyzeBData()   { return m_gcAnalyzeB; }
#endif

protected:
  TEncRateCtrl* getRateCtrl()       { return m_pcRateCtrl;  }

protected:

  Void  xInitGOP          ( Int iPOCLast, Int iNumPicRcvd, Bool isField );
  Void  xGetBuffer        ( TComList<TComPic*>& rcListPic, TComList<TComPicYuv*>& rcListPicYuvRecOut, Int iNumPicRcvd, Int iTimeOffset, TComPic*& rpcPic, TComPicYuv*& rpcPicYuvRecOut, Int pocCurr, Bool isField );

#if JVET_F0064_MSSSIM
  Void  xCalculateAddPSNRs         ( const Bool isField, const Bool isFieldTopFieldFirst, const Int iGOPid, TComPic* pcPic, const AccessUnit&accessUnit, TComList<TComPic*> &rcListPic, Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE, const Bool printMSSSIM );
  Void  xCalculateAddPSNR          ( TComPic* pcPic, TComPicYuv* pcPicD, const AccessUnit&, Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE, const Bool printMSSSIM );
  Void  xCalculateInterlacedAddPSNR( TComPic* pcPicOrgFirstField, TComPic* pcPicOrgSecondField,
                                     TComPicYuv* pcPicRecFirstField, TComPicYuv* pcPicRecSecondField,
                                     const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE, const Bool printMSSSIM );

  Double xCalculateMSSSIM (const Pel *pOrg, const Int iOrgStride, const Pel* pRec, const Int iRecStride, const Int iWidth, const Int iHeight, const UInt uiBitDepth);
#else
  Void  xCalculateAddPSNRs         ( const Bool isField, const Bool isFieldTopFieldFirst, const Int iGOPid, TComPic* pcPic, const AccessUnit&accessUnit, TComList<TComPic*> &rcListPic, Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE);
  Void  xCalculateAddPSNR          ( TComPic* pcPic, TComPicYuv* pcPicD, const AccessUnit&, Double dEncTime, const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE);
  Void  xCalculateInterlacedAddPSNR( TComPic* pcPicOrgFirstField, TComPic* pcPicOrgSecondField,
                                     TComPicYuv* pcPicRecFirstField, TComPicYuv* pcPicRecSecondField,
                                     const InputColourSpaceConversion snr_conversion, const Bool printFrameMSE);
#endif

  UInt64 xFindDistortionFrame (TComPicYuv* pcPic0, TComPicYuv* pcPic1, const BitDepths &bitDepths);

  Double xCalculateRVM();

  Void xCreateIRAPLeadingSEIMessages (SEIMessages& seiMessages, const TComSPS *sps, const TComPPS *pps);
  Void xCreatePerPictureSEIMessages (Int picInGOP, SEIMessages& seiMessages, SEIMessages& nestedSeiMessages, TComSlice *slice);
  Void xCreatePictureTimingSEI  (Int IRAPGOPid, SEIMessages& seiMessages, SEIMessages& nestedSeiMessages, SEIMessages& duInfoSeiMessages, TComSlice *slice, Bool isField, std::deque<DUData> &duData);
  Void xUpdateDuData(AccessUnit &testAU, std::deque<DUData> &duData);
  Void xUpdateTimingSEI(SEIPictureTiming *pictureTimingSEI, std::deque<DUData> &duData, const TComSPS *sps);
  Void xUpdateDuInfoSEI(SEIMessages &duInfoSeiMessages, SEIPictureTiming *pictureTimingSEI);

  Void xCreateScalableNestingSEI (SEIMessages& seiMessages, SEIMessages& nestedSeiMessages);
  Void xWriteSEI (NalUnitType naluType, SEIMessages& seiMessages, AccessUnit &accessUnit, AccessUnit::iterator &auPos, Int temporalId, const TComSPS *sps);
  Void xWriteSEISeparately (NalUnitType naluType, SEIMessages& seiMessages, AccessUnit &accessUnit, AccessUnit::iterator &auPos, Int temporalId, const TComSPS *sps);
  Void xClearSEIs(SEIMessages& seiMessages, Bool deleteMessages);
  Void xWriteLeadingSEIOrdered (SEIMessages& seiMessages, SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, Bool testWrite);
  Void xWriteLeadingSEIMessages  (SEIMessages& seiMessages, SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, std::deque<DUData> &duData);
  Void xWriteTrailingSEIMessages (SEIMessages& seiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps);
  Void xWriteDuSEIMessages       (SEIMessages& duInfoSeiMessages, AccessUnit &accessUnit, Int temporalId, const TComSPS *sps, std::deque<DUData> &duData);

  Int xWriteVPS (AccessUnit &accessUnit, const TComVPS *vps);
  Int xWriteSPS (AccessUnit &accessUnit, const TComSPS *sps);
  Int xWritePPS (AccessUnit &accessUnit, const TComPPS *pps);
  Int xWriteParameterSets (AccessUnit &accessUnit, TComSlice *slice);

  Void applyDeblockingFilterMetric( TComPic* pcPic, UInt uiNumSlices );
};// END CLASS DEFINITION TEncGOP

//! \}

#endif // __TENCGOP__

