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

/** \file     TEncTop.h
    \brief    encoder class (header)
*/

#ifndef __TENCTOP__
#define __TENCTOP__

// Include files
#include "TLibCommon/TComList.h"
#include "TLibCommon/TComPrediction.h"
#include "TLibCommon/TComTrQuant.h"
#include "TLibCommon/TComLoopFilter.h"
#include "TLibCommon/AccessUnit.h"

#include "TLibVideoIO/TVideoIOYuv.h"

#include "TEncCfg.h"
#include "TEncGOP.h"
#include "TEncSlice.h"
#include "TEncEntropy.h"
#include "TEncCavlc.h"
#include "TEncSbac.h"
#include "TEncSearch.h"
#include "TEncSampleAdaptiveOffset.h"
#include "TEncPreanalyzer.h"
#include "TEncRateCtrl.h"
#if CNN_BASED_INLOOP_FILTER
#include "TEncCNNLoopFilter.h"
#endif

//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// encoder class
class TEncTop : public TEncCfg
{
private:
  // picture
  Int                     m_iPOCLast;                     ///< time index (POC)
  Int                     m_iNumPicRcvd;                  ///< number of received pictures
  UInt                    m_uiNumAllPicCoded;             ///< number of coded pictures
#if BARC_REF_UP_SAMPLE_CNN
  TComConvNet             m_cCNNUpSampleImg;
#if COMPACT_RESOLUTION_CNN
  TComConvNet             m_cCNNDownSampleImg;
#endif
#endif

  TComList<TComPic*>      m_cListPic;                     ///< dynamic list of pictures

  // encoder search
  TEncSearch              m_cSearch;                      ///< encoder search class
  //TEncEntropy*            m_pcEntropyCoder;                     ///< entropy encoder
  TEncCavlc*              m_pcCavlcCoder;                       ///< CAVLC encoder
  // coding tool
  TComTrQuant             m_cTrQuant;                     ///< transform & quantization class
  TComLoopFilter          m_cLoopFilter;                  ///< deblocking filter class
  TEncSampleAdaptiveOffset m_cEncSAO;                     ///< sample adaptive offset class
#if CNN_BASED_INLOOP_FILTER
	TEncCNNLoopFilter       m_cEncCNNLoopFilter;
#endif
  TEncEntropy             m_cEntropyCoder;                ///< entropy encoder
  TEncCavlc               m_cCavlcCoder;                  ///< CAVLC encoder
  TEncSbac                m_cSbacCoder;                   ///< SBAC encoder
  TEncBinCABAC            m_cBinCoderCABAC;               ///< bin coder CABAC

  // processing unit
  TEncGOP                 m_cGOPEncoder;                  ///< GOP encoder
  TEncSlice               m_cSliceEncoder;                ///< slice encoder
  TEncCu                  m_cCuEncoder;                   ///< CU encoder
  // SPS
  TComSPS                 m_cSPS;                         ///< SPS. This is the base value. This is copied to TComPicSym
  TComPPS                 m_cPPS;                         ///< PPS. This is the base value. This is copied to TComPicSym
  // RD cost computation
  TComRdCost              m_cRdCost;                      ///< RD cost computation class
#if JVET_C0024_QTBT
  TEncSbac****            m_ppppcRDSbacCoder;              ///< temporal storage for RD computation
#else
  TEncSbac***             m_pppcRDSbacCoder;              ///< temporal storage for RD computation
#endif
  TEncSbac                m_cRDGoOnSbacCoder;             ///< going on SBAC model for RD stage
#if FAST_BIT_EST
#if JVET_C0024_QTBT
  TEncBinCABACCounter**** m_ppppcBinCoderCABAC;            ///< temporal CABAC state storage for RD computation
#else
  TEncBinCABACCounter***  m_pppcBinCoderCABAC;            ///< temporal CABAC state storage for RD computation
#endif
  TEncBinCABACCounter     m_cRDGoOnBinCoderCABAC;         ///< going on bin coder CABAC for RD stage
#else
  TEncBinCABAC***         m_pppcBinCoderCABAC;            ///< temporal CABAC state storage for RD computation
  TEncBinCABAC            m_cRDGoOnBinCoderCABAC;         ///< going on bin coder CABAC for RD stage
#endif

  // quality control
  TEncPreanalyzer         m_cPreanalyzer;                 ///< image characteristics analyzer for TM5-step3-like adaptive QP

  TEncRateCtrl            m_cRateCtrl;                    ///< Rate control class

#if ALF_HM3_REFACTOR
  TEncAdaptiveLoopFilter  m_cAdaptiveLoopFilter;          ///< adaptive loop filter class
#endif

protected:
  Void  xGetNewPicBuffer  ( TComPic*& rpcPic );           ///< get picture buffer which will be processed
  Void  xInitVPS          ();                             ///< initialize VPS from encoder options
  Void  xInitSPS          ();                             ///< initialize SPS from encoder options
  Void  xInitPPS          ();                             ///< initialize PPS from encoder options
  Void  xInitScalingLists ();                             ///< initialize scaling lists
  Void  xInitHrdParameters();                             ///< initialize HRD parameters

  Void  xInitPPSforTiles  ();
  Void  xInitRPS          (Bool isFieldCoding);           ///< initialize PPS from encoder options

public:
  TEncTop();
  virtual ~TEncTop();

  Void      create          ();
  Void      destroy         ();
  Void      init            (Bool isFieldCoding);
  Void      deletePicBuffer ();

  // -------------------------------------------------------------------------------------------------------------------
  // member access functions
  // -------------------------------------------------------------------------------------------------------------------

  TComList<TComPic*>*     getListPic            () { return  &m_cListPic;             }
  TEncSearch*             getPredSearch         () { return  &m_cSearch;              }

  TComTrQuant*            getTrQuant            () { return  &m_cTrQuant;             }
  TComLoopFilter*         getLoopFilter         () { return  &m_cLoopFilter;          }
  TEncSampleAdaptiveOffset* getSAO              () { return  &m_cEncSAO;              }
#if CNN_BASED_INLOOP_FILTER
	TEncCNNLoopFilter*      getCNNLoopFilter      () { return  &m_cEncCNNLoopFilter;    }
#endif
  TEncGOP*                getGOPEncoder         () { return  &m_cGOPEncoder;          }
  TEncSlice*              getSliceEncoder       () { return  &m_cSliceEncoder;        }
  TEncCu*                 getCuEncoder          () { return  &m_cCuEncoder;           }
  TEncEntropy*            getEntropyCoder       () { return  &m_cEntropyCoder;        }
  TEncCavlc*              getCavlcCoder         () { return  &m_cCavlcCoder;          }
  TEncSbac*               getSbacCoder          () { return  &m_cSbacCoder;           }
  TEncBinCABAC*           getBinCABAC           () { return  &m_cBinCoderCABAC;       }

  TComRdCost*             getRdCost             () { return  &m_cRdCost;              }
#if JVET_C0024_QTBT
  TEncSbac****            getRDSbacCoder        () { return  m_ppppcRDSbacCoder;       }
#else
  TEncSbac***             getRDSbacCoder        () { return  m_pppcRDSbacCoder;       }
#endif
  TEncSbac*               getRDGoOnSbacCoder    () { return  &m_cRDGoOnSbacCoder;     }
  TEncRateCtrl*           getRateCtrl           () { return &m_cRateCtrl;             }
  Void selectReferencePictureSet(TComSlice* slice, Int POCCurr, Int GOPid );
  Int getReferencePictureSetIdxForSOP(Int POCCurr, Int GOPid );
  // -------------------------------------------------------------------------------------------------------------------
  // encoder function
  // -------------------------------------------------------------------------------------------------------------------

  /// encode several number of pictures until end-of-sequence
  Void encode( Bool bEos,
               TComPicYuv* pcPicYuvOrg,
               TComPicYuv* pcPicYuvTrueOrg, const InputColourSpaceConversion snrCSC, // used for SNR calculations. Picture in original colour space.
               TComList<TComPicYuv*>& rcListPicYuvRecOut,
               std::list<AccessUnit>& accessUnitsOut, Int& iNumEncoded 
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
              ,TComStats* m_apcStats
#endif
               );

  /// encode several number of pictures until end-of-sequence
  Void encode( Bool bEos, TComPicYuv* pcPicYuvOrg,
               TComPicYuv* pcPicYuvTrueOrg, const InputColourSpaceConversion snrCSC, // used for SNR calculations. Picture in original colour space.
               TComList<TComPicYuv*>& rcListPicYuvRecOut,
               std::list<AccessUnit>& accessUnitsOut, Int& iNumEncoded, Bool isTff
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
              ,TComStats* m_apcStats
#endif
               );

#if JVET_D0134_PSNR
#if JVET_F0064_MSSSIM
  Void printSummary(Bool isField) { m_cGOPEncoder.printOutSummary(m_uiNumAllPicCoded, isField, m_printMSEBasedSequencePSNR, m_printSequenceMSE, m_printMSSSIM, m_trueBitdepthPSNR, m_cSPS.getBitDepths()); }
#else
  Void printSummary(Bool isField) { m_cGOPEncoder.printOutSummary(m_uiNumAllPicCoded, isField, m_printMSEBasedSequencePSNR, m_printSequenceMSE, m_trueBidepthPSNR, m_cSPS.getBitDepths()); }
#endif
#else
#if JVET_F0064_MSSSIM
  Void printSummary(Bool isField) { m_cGOPEncoder.printOutSummary(m_uiNumAllPicCoded, isField, m_printMSEBasedSequencePSNR, m_printSequenceMSE, m_printMSSSIM, m_cSPS.getBitDepths()); }
#else
  Void printSummary(Bool isField) { m_cGOPEncoder.printOutSummary(m_uiNumAllPicCoded, isField, m_printMSEBasedSequencePSNR, m_printSequenceMSE, m_cSPS.getBitDepths()); }
#endif
#endif

#if JVET_D0186_PRECISEPSNR
  Void printPreciseSummary(const Char * filename, Bool isField) { 
    m_cGOPEncoder.printOutPreciseSummary(filename, isField
#if JVET_D0134_PSNR
      , m_trueBitdepthPSNR
#endif
      , m_cSPS.getBitDepths()
    );
  }
#endif

#if ALF_HM3_REFACTOR
  TEncAdaptiveLoopFilter* getAdaptiveLoopFilter () { return  &m_cAdaptiveLoopFilter;  }
#endif

#if PVC_AQS
  //analyze the original input picture to obtain a error sensitivity map for each pixel
  ESM*   m_pesmErrorSensitivityMap;
public:
  ESM*   getESMap()   { return m_pesmErrorSensitivityMap; }
#endif
};

//! \}

#endif // __TENCTOP__

