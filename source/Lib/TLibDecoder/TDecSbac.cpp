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

/** \file     TDecSbac.cpp
    \brief    Context-adaptive entropy decoder class
*/

#include "TDecSbac.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/TComTrQuant.h"

#if RExt__DECODER_DEBUG_BIT_STATISTICS
#include "TLibCommon/TComCodingStatistics.h"
//
#define RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(a) , a
#else
#define RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(a)
#endif

//! \ingroup TLibDecoder
//! \{

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
#include "../TLibCommon/Debug.h"
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TDecSbac::TDecSbac()
// new structure here
: m_pcBitstream                              ( 0 )
, m_pcTDecBinIf                              ( NULL )
, m_numContextModels                         ( 0 )
, m_cCUSplitFlagSCModel                      ( 1,             1,                      NUM_SPLIT_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#if JVET_C0024_QTBT
, m_cBTSplitFlagSCModel                      ( 1,             1,                      NUM_BTSPLIT_MODE_CTX                 , m_contextModels + m_numContextModels, m_numContextModels )
#endif
, m_cCUSkipFlagSCModel                       ( 1,             1,                      NUM_SKIP_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#if BARC_UP_SAMPLE_RDO
, m_cCUUpCNNFlagSCModel                      ( 1,             1,                      NUM_BARC_UP_FILT_RDO_FLAG_CTX        , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if BARC
, m_cCUBARCFlagSCModel                       ( 1,             1,                      NUM_BARC_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ05_INTRA_MPI
, m_cMPIIdxSCModel                           ( 1,             1,                      NUM_MPI_CTX                          , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if COM16_C1046_PDPC_INTRA
, m_cPDPCIdxSCModel                          ( 1,             1,                      NUM_PDPC_CTX                         , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if DIMD_INTRA_PRED
, m_cDIMDEnabledSCModel(1, 1, NUM_DIMD_CTX, m_contextModels + m_numContextModels, m_numContextModels)
, m_cDIMDNoBTFlagSCModel(1, 1, NUM_DIMD_NO_BT_CTX, m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_DWDIP
, m_cDWDIPIdxSCModel                         ( 1,             1,                      NUM_DWDIP_CTX                        , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_DWDIP_MLR
, m_cDWDIP_GPIdxSCModel                         ( 1,             1,                       NUM_DWDIP_GP_CTX                         , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
, m_cROTidxSCModel                           ( 1,             1,                      NUM_ROT_TR_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if  FAST_EMT_v1
, m_cEMTIntraidxSCModel                      (1,              1,                     NUM_EMT_TR_CTX                ,m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cCUMergeFlagExtSCModel                   ( 1,             1,                      NUM_MERGE_FLAG_EXT_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMergeIdxExtSCModel                    ( 1,             1,                      NUM_MERGE_IDX_EXT_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#if VCEG_AZ07_FRUC_MERGE
, m_cCUFRUCMgrModeSCModel                    ( 1,             1,                      NUM_FRUCMGRMODE_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUFRUCMESCModel                         ( 1,             1,                      NUM_FRUCME_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
#if THIRD_FRUC_MODE
, m_cCUFRUCMESCModel_uni                     ( 1,             1,                      NUM_FRUCME_UNI_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#endif
, m_cCUPartSizeSCModel                       ( 1,             1,                      NUM_PART_SIZE_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUPredModeSCModel                       ( 1,             1,                      NUM_PRED_MODE_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUIntraPredSCModel                      ( 1,             1,                      NUM_INTRA_PREDICT_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#if MLIP_EXTEND_REF_LINE
, m_cCUExternRefSCModel                     (1,               1,                    NUM_EXTERN_REF_CTX,                    m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cCUChromaPredSCModel                     ( 1,             1,                      NUM_CHROMA_PRED_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUDeltaQpSCModel                        ( 1,             1,                      NUM_DELTA_QP_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUInterDirSCModel                       ( 1,             1,                      NUM_INTER_DIR_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCURefPicSCModel                         ( 1,             1,                      NUM_REF_NO_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMvdSCModel                            ( 1,             1,                      NUM_MV_RES_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUQtCbfSCModel                          ( 1,             NUM_QT_CBF_CTX_SETS,    NUM_QT_CBF_CTX_PER_SET               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUTransSubdivFlagSCModel                ( 1,             1,                      NUM_TRANS_SUBDIV_FLAG_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUQtRootCbfSCModel                      ( 1,             1,                      NUM_QT_ROOT_CBF_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUSigCoeffGroupSCModel                  ( 1,             2,                      NUM_SIG_CG_FLAG_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUSigSCModel                            ( 1,             1,                      NUM_SIG_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCuCtxLastX                              ( 1,             NUM_CTX_LAST_FLAG_SETS, NUM_CTX_LAST_FLAG_XY                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCuCtxLastY                              ( 1,             NUM_CTX_LAST_FLAG_SETS, NUM_CTX_LAST_FLAG_XY                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUOneSCModel                            ( 1,             1,                      NUM_ONE_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#if !VCEG_AZ07_CTX_RESIDUALCODING
, m_cCUAbsSCModel                            ( 1,             1,                      NUM_ABS_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cMVPIdxSCModel                           ( 1,             1,                      NUM_MVP_IDX_CTX                      , m_contextModels + m_numContextModels, m_numContextModels)
#if CNN_BASED_INLOOP_FILTER
, m_cCUCNNLoopFilterFlagSCModel              ( 1,             1,                      NUM_CNNLF_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cSaoMergeSCModel                         ( 1,             1,                      NUM_SAO_MERGE_FLAG_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cSaoTypeIdxSCModel                       ( 1,             1,                      NUM_SAO_TYPE_IDX_CTX                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cTransformSkipSCModel                    ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_TRANSFORMSKIP_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
#if VCEG_AZ08_KLT_COMMON
, m_cKLTFlagSCModel                          ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_KLT_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_CUTransquantBypassFlagSCModel            ( 1,             1,                      NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX    , m_contextModels + m_numContextModels, m_numContextModels)
, m_explicitRdpcmFlagSCModel                 ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_EXPLICIT_RDPCM_FLAG_CTX          , m_contextModels + m_numContextModels, m_numContextModels)
, m_explicitRdpcmDirSCModel                  ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_EXPLICIT_RDPCM_DIR_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCrossComponentPredictionSCModel         ( 1,             1,                      NUM_CROSS_COMPONENT_PREDICTION_CTX   , m_contextModels + m_numContextModels, m_numContextModels)
, m_ChromaQpAdjFlagSCModel                   ( 1,             1,                      NUM_CHROMA_QP_ADJ_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
, m_ChromaQpAdjIdcSCModel                    ( 1,             1,                      NUM_CHROMA_QP_ADJ_IDC_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
#if COM16_C806_OBMC
, m_cCUOBMCFlagSCModel                       ( 1,             1,                      NUM_OBMC_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ07_IMV 
, m_cCUiMVFlagSCModel                        ( 1,             1,                      NUM_IMV_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ06_IC
, m_cCUICFlagSCModel                         ( 1,             1,                      NUM_IC_FLAG_CTX                      , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if ALF_HM3_REFACTOR
, m_cCUAlfCtrlFlagSCModel                    ( 1,             1,                      NUM_ALF_CTRL_FLAG_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#if !JVET_C0038_GALF
, m_cALFFlagSCModel                          ( 1,             1,                      NUM_ALF_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cALFUvlcSCModel                          ( 1,             1,                      NUM_ALF_UVLC_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#if !JVET_C0038_GALF
, m_cALFSvlcSCModel                          ( 1,             1,                      NUM_ALF_SVLC_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#endif
#if COM16_C806_EMT
, m_cEmtTuIdxSCModel                         ( 1,             1,                      NUM_EMT_TU_IDX_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
, m_cEmtCuFlagSCModel                        ( 1,             1,                      NUM_EMT_CU_FLAG_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if COM16_C1016_AFFINE
, m_cCUAffineFlagSCModel                     ( 1,             1,                      NUM_AFFINE_FLAG_CTX                  , m_contextModels + m_numContextModels,          m_numContextModels)
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
, m_CUBgFlagSCModel                          ( 1,             1,                      NUM_CU_BG_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if MRG_CONSIDER_COMPLEX
, m_cCUMergeComplexFlagExtSCModel            ( 1,             1,               NUM_MERGE_COMPLEX_FLAG_EXT_CTX       , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMergeComplexIdxExtSCModel             ( 1,             1,               NUM_MERGE_COMPLEX_IDX_EXT_CTX        , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
, m_cCUAffineParamFlagSCModel                ( 1,             1,               NUM_AFFINE_FLAG_CTX                  , m_contextModels + m_numContextModels,        m_numContextModels)
#endif

#if RRC_NS
, m_cNsFlagSCModel(1, 1, NUM_NS_FLAG_CTX, m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_INTRA_SHARPENING
, m_cSharpFlagSCModel                        ( 1,             1,                      NUM_SHARP_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_SELECTABLE_IF
, m_cCUAdaptiveIfIndex                       ( 1,             1,                      NUM_SELECTABLE_IF_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#endif
{
  assert( m_numContextModels <= MAX_NUM_CTX_MOD );
}

TDecSbac::~TDecSbac()
{
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================
#if VCEG_AZ07_BAC_ADAPT_WDOW
Void TDecSbac::xUpdateWindowSize (SliceType eSliceType, Int uiQPIdx, TComStats* apcStats)
{
  if( uiQPIdx == -1 )
  {
    return;
  }
  Int iCtxNr = getCtxNumber();
  for(UInt i=0; i< iCtxNr; i++)
  {
    m_contextModels[i].setWindowSize(apcStats->m_uiCtxCodeIdx[eSliceType][uiQPIdx][i]);
  }
}
#endif

Void TDecSbac::resetEntropy(TComSlice* pSlice)
{
  SliceType sliceType  = pSlice->getSliceType();
  Int       qp         = pSlice->getSliceQp();

  if (pSlice->getPPS()->getCabacInitPresentFlag() && pSlice->getCabacInitFlag())
  {
    switch (sliceType)
    {
    case P_SLICE:           // change initialization table to B_SLICE initialization
      sliceType = B_SLICE;
      break;
    case B_SLICE:           // change initialization table to P_SLICE initialization
      sliceType = P_SLICE;
      break;
    default     :           // should not occur
      assert(0);
      break;
    }
  }

  m_cCUSplitFlagSCModel.initBuffer                ( sliceType, qp, (UChar*)INIT_SPLIT_FLAG );
#if JVET_C0024_QTBT
  m_cBTSplitFlagSCModel.initBuffer                ( sliceType, qp, (UChar*)INIT_BTSPLIT_MODE );
#endif
  m_cCUSkipFlagSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_SKIP_FLAG );
#if BARC_UP_SAMPLE_RDO
  m_cCUUpCNNFlagSCModel.initBuffer                ( sliceType, qp, (UChar*)INIT_UP_CNN_RDO_FLAG);
#endif
#if BARC
  m_cCUBARCFlagSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_BARC_FLAG);
#endif
#if VCEG_AZ05_INTRA_MPI
  m_cMPIIdxSCModel.initBuffer                     ( sliceType, qp, (UChar*)INIT_MPIIdx_FLAG );
#endif
#if COM16_C1046_PDPC_INTRA
  m_cPDPCIdxSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_PDPCIdx_FLAG);
#endif
#if DIMD_INTRA_PRED
  m_cDIMDEnabledSCModel.initBuffer(sliceType, qp, (UChar*)INIT_DIMD_FLAG);
  m_cDIMDNoBTFlagSCModel.initBuffer(sliceType, qp, (UChar*)INIT_DIMD_NO_BT_FLAG);
#endif
#if RRC_DWDIP
  m_cDWDIPIdxSCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_DWDIPIdx_FLAG);
#endif
#if RRC_DWDIP_MLR
  m_cDWDIP_GPIdxSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_DWDIP_GPIdx_FLAG);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  m_cROTidxSCModel.initBuffer        ( sliceType, qp, (UChar*)INIT_ROT_TR_IDX );
#endif
#if FAST_EMT_v1
  m_cEMTIntraidxSCModel.initBuffer   (sliceType, qp, (UChar*)INIT_EMT_TR_IDX);
#endif
  m_cCUMergeFlagExtSCModel.initBuffer             ( sliceType, qp, (UChar*)INIT_MERGE_FLAG_EXT );
  m_cCUMergeIdxExtSCModel.initBuffer              ( sliceType, qp, (UChar*)INIT_MERGE_IDX_EXT );
#if VCEG_AZ07_FRUC_MERGE
  m_cCUFRUCMgrModeSCModel.initBuffer              ( sliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN1 );
  m_cCUFRUCMESCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN2 );
#if THIRD_FRUC_MODE
  m_cCUFRUCMESCModel_uni.initBuffer               ( sliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN3 );
#endif
#endif
  m_cCUPartSizeSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_PART_SIZE );
  m_cCUPredModeSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_PRED_MODE );
  m_cCUIntraPredSCModel.initBuffer                ( sliceType, qp, (UChar*)INIT_INTRA_PRED_MODE );
#if MLIP_EXTEND_REF_LINE
  m_cCUExternRefSCModel.initBuffer              (sliceType, qp, (UChar*)INIT_EXTERN_REF);
#endif
  m_cCUChromaPredSCModel.initBuffer               ( sliceType, qp, (UChar*)INIT_CHROMA_PRED_MODE );
  m_cCUInterDirSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_INTER_DIR );
  m_cCUMvdSCModel.initBuffer                      ( sliceType, qp, (UChar*)INIT_MVD );
  m_cCURefPicSCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_REF_PIC );
  m_cCUDeltaQpSCModel.initBuffer                  ( sliceType, qp, (UChar*)INIT_DQP );
  m_cCUQtCbfSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_QT_CBF );
  m_cCUQtRootCbfSCModel.initBuffer                ( sliceType, qp, (UChar*)INIT_QT_ROOT_CBF );
  m_cCUSigCoeffGroupSCModel.initBuffer            ( sliceType, qp, (UChar*)INIT_SIG_CG_FLAG );
  m_cCUSigSCModel.initBuffer                      ( sliceType, qp, (UChar*)INIT_SIG_FLAG );
  m_cCuCtxLastX.initBuffer                        ( sliceType, qp, (UChar*)INIT_LAST );
  m_cCuCtxLastY.initBuffer                        ( sliceType, qp, (UChar*)INIT_LAST );
  m_cCUOneSCModel.initBuffer                      ( sliceType, qp, (UChar*)INIT_ONE_FLAG );
#if !VCEG_AZ07_CTX_RESIDUALCODING
  m_cCUAbsSCModel.initBuffer                      ( sliceType, qp, (UChar*)INIT_ABS_FLAG );
#endif
  m_cMVPIdxSCModel.initBuffer                     ( sliceType, qp, (UChar*)INIT_MVP_IDX );
#if CNN_BASED_INLOOP_FILTER
	m_cCUCNNLoopFilterFlagSCModel.initBuffer        ( sliceType, qp, (UChar*)INIT_CNNLF_FLAG);
#endif
  m_cSaoMergeSCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_SAO_MERGE_FLAG );
  m_cSaoTypeIdxSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_SAO_TYPE_IDX );
  m_cCUTransSubdivFlagSCModel.initBuffer          ( sliceType, qp, (UChar*)INIT_TRANS_SUBDIV_FLAG );
  m_cTransformSkipSCModel.initBuffer              ( sliceType, qp, (UChar*)INIT_TRANSFORMSKIP_FLAG );
#if VCEG_AZ08_KLT_COMMON
  m_cKLTFlagSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_KLT_FLAG );
#endif
  m_CUTransquantBypassFlagSCModel.initBuffer      ( sliceType, qp, (UChar*)INIT_CU_TRANSQUANT_BYPASS_FLAG );
  m_explicitRdpcmFlagSCModel.initBuffer           ( sliceType, qp, (UChar*)INIT_EXPLICIT_RDPCM_FLAG);
  m_explicitRdpcmDirSCModel.initBuffer            ( sliceType, qp, (UChar*)INIT_EXPLICIT_RDPCM_DIR);
  m_cCrossComponentPredictionSCModel.initBuffer   ( sliceType, qp, (UChar*)INIT_CROSS_COMPONENT_PREDICTION );
  m_ChromaQpAdjFlagSCModel.initBuffer             ( sliceType, qp, (UChar*)INIT_CHROMA_QP_ADJ_FLAG );
  m_ChromaQpAdjIdcSCModel.initBuffer              ( sliceType, qp, (UChar*)INIT_CHROMA_QP_ADJ_IDC );
#if VCEG_AZ07_IMV
  m_cCUiMVFlagSCModel.initBuffer                  ( sliceType, qp, (UChar*)INIT_IMV_FLAG );
#endif
#if COM16_C806_OBMC
  m_cCUOBMCFlagSCModel.initBuffer                 ( sliceType, qp, (UChar*)INIT_OBMC_FLAG );
#endif
#if VCEG_AZ06_IC
  m_cCUICFlagSCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_IC_FLAG );
#endif  
#if ALF_HM3_REFACTOR
  m_cCUAlfCtrlFlagSCModel.initBuffer              ( sliceType, qp, (UChar*)INIT_ALF_CTRL_FLAG );
#if !JVET_C0038_GALF
  m_cALFFlagSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_ALF_FLAG );
#endif
  m_cALFUvlcSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_ALF_UVLC );
#if !JVET_C0038_GALF
  m_cALFSvlcSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_ALF_SVLC );
#endif
#endif
#if COM16_C806_EMT
  m_cEmtTuIdxSCModel.initBuffer                   ( sliceType, qp, (UChar*)INIT_EMT_TU_IDX );
  m_cEmtCuFlagSCModel.initBuffer                  ( sliceType, qp, (UChar*)INIT_EMT_CU_FLAG );
#endif
#if COM16_C1016_AFFINE
  m_cCUAffineFlagSCModel.initBuffer               ( sliceType, qp, (UChar*)INIT_AFFINE_FLAG );
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
	m_CUBgFlagSCModel.initBuffer                    ( sliceType, qp, (UChar*)INIT_CU_BG_FLAG );
#endif
#if MRG_CONSIDER_COMPLEX
  m_cCUMergeComplexFlagExtSCModel.initBuffer      ( sliceType, qp, (UChar*)INIT_MERGE_COMPLEX_FLAG_EXT );
  m_cCUMergeComplexIdxExtSCModel.initBuffer       ( sliceType, qp, (UChar*)INIT_MERGE_COMPLEX_IDX_EXT );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  m_cCUAffineParamFlagSCModel.initBuffer          (sliceType, qp, (UChar*)INIT_AFFINE_FLAG );
#endif

#if RRC_NS
	m_cNsFlagSCModel.initBuffer(sliceType, qp, (UChar*)INIT_NS_FLAG);
#endif
#if RRC_INTRA_SHARPENING
  m_cSharpFlagSCModel.initBuffer                  ( sliceType, qp, (UChar*)INIT_SHARP_FLAG );
#endif
#if RRC_SELECTABLE_IF
  m_cCUAdaptiveIfIndex.initBuffer                 ( sliceType, qp, (UChar*)INIT_ADAPTIVE_IF_INDEX );
#endif

  for (UInt statisticIndex = 0; statisticIndex < RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS ; statisticIndex++)
  {
    m_golombRiceAdaptationStatistics[statisticIndex] = 0;
  }

  m_pcTDecBinIf->start();
#if VCEG_AZ07_BAC_ADAPT_WDOW
  xUpdateWindowSize (pSlice->getSliceType(), pSlice->getCtxMapQPIdx(), pSlice->getStatsHandle()); 
#endif
}

Void TDecSbac::parseTerminatingBit( UInt& ruiBit )
{
  m_pcTDecBinIf->decodeBinTrm( ruiBit );
  if ( ruiBit == 1 )
  {
    m_pcTDecBinIf->finish();

#if RExt__DECODER_DEBUG_BIT_STATISTICS
    TComCodingStatistics::IncrementStatisticEP(STATS__TRAILING_BITS, m_pcBitstream->readOutTrailingBits(),0);
#else
    m_pcBitstream->readOutTrailingBits();
#endif
  }
}

Void TDecSbac::parseRemainingBytes( Bool noTrailingBytesExpected )
{
  if (noTrailingBytesExpected)
  {
    const UInt numberOfRemainingSubstreamBytes=m_pcBitstream->getNumBitsLeft();
    assert (numberOfRemainingSubstreamBytes == 0);
  }
  else
  {
    while (m_pcBitstream->getNumBitsLeft())
    {
      UInt trailingNullByte=m_pcBitstream->readByte();
      if (trailingNullByte!=0)
      {
        printf("Trailing byte should be 0, but has value %02x\n", trailingNullByte);
        assert(trailingNullByte==0);
      }
    }
  }
}

#if RExt__DECODER_DEBUG_BIT_STATISTICS
Void TDecSbac::xReadUnaryMaxSymbol( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol, const class TComCodingStatisticsClassType &whichStat )
#else
Void TDecSbac::xReadUnaryMaxSymbol( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol )
#endif
{
  if (uiMaxSymbol == 0)
  {
    ruiSymbol = 0;
    return;
  }

  m_pcTDecBinIf->decodeBin( ruiSymbol, pcSCModel[0] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );

  if( ruiSymbol == 0 || uiMaxSymbol == 1 )
  {
    return;
  }

  UInt uiSymbol = 0;
  UInt uiCont;

  do
  {
    m_pcTDecBinIf->decodeBin( uiCont, pcSCModel[ iOffset ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );
    uiSymbol++;
  } while( uiCont && ( uiSymbol < uiMaxSymbol - 1 ) );

  if( uiCont && ( uiSymbol == uiMaxSymbol - 1 ) )
  {
    uiSymbol++;
  }

  ruiSymbol = uiSymbol;
}

#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
#if RExt__DECODER_DEBUG_BIT_STATISTICS
Void TDecSbac::xReadTruncBinCode(UInt& ruiSymbol, UInt uiMaxSymbol, const class TComCodingStatisticsClassType &whichStat)
#else
Void TDecSbac::xReadTruncBinCode(UInt& ruiSymbol, UInt uiMaxSymbol)
#endif
{
  UInt uiThresh;
  if (uiMaxSymbol > 256)
  {
    UInt uiThreshVal = 1 << 8;
    uiThresh = 8;
    while (uiThreshVal <= uiMaxSymbol)
    {
      uiThresh++;
      uiThreshVal <<= 1;
    }
    uiThresh--;
  }
  else
  {
    uiThresh = g_NonMPM[uiMaxSymbol];
  }

  UInt uiVal = 1 << uiThresh;
  UInt b = uiMaxSymbol - uiVal;
  m_pcTDecBinIf->decodeBinsEP(ruiSymbol, uiThresh RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
  if (ruiSymbol >= uiVal - b)
  {
    UInt uiSymbol;
    m_pcTDecBinIf->decodeBinEP(uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
    ruiSymbol <<= 1;
    ruiSymbol += uiSymbol;
    ruiSymbol -= (uiVal - b);
  }
}
#endif

#if RExt__DECODER_DEBUG_BIT_STATISTICS
Void TDecSbac::xReadEpExGolomb( UInt& ruiSymbol, UInt uiCount, const class TComCodingStatisticsClassType &whichStat )
#else
Void TDecSbac::xReadEpExGolomb( UInt& ruiSymbol, UInt uiCount )
#endif
{
  UInt uiSymbol = 0;
  UInt uiBit = 1;

  while( uiBit )
  {
    m_pcTDecBinIf->decodeBinEP( uiBit RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );
    uiSymbol += uiBit << uiCount++;
  }

  if ( --uiCount )
  {
    UInt bins;
    m_pcTDecBinIf->decodeBinsEP( bins, uiCount RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );
    uiSymbol += bins;
  }

  ruiSymbol = uiSymbol;
}

#if RExt__DECODER_DEBUG_BIT_STATISTICS
Void TDecSbac::xReadUnarySymbol( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset, const class TComCodingStatisticsClassType &whichStat )
#else
Void TDecSbac::xReadUnarySymbol( UInt& ruiSymbol, ContextModel* pcSCModel, Int iOffset )
#endif
{
  m_pcTDecBinIf->decodeBin( ruiSymbol, pcSCModel[0] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));

  if( !ruiSymbol )
  {
    return;
  }

  UInt uiSymbol = 0;
  UInt uiCont;

  do
  {
    m_pcTDecBinIf->decodeBin( uiCont, pcSCModel[ iOffset ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
    uiSymbol++;
  } while( uiCont );

  ruiSymbol = uiSymbol;
}

/** Parsing of coeff_abs_level_remaing
 * \param rSymbol                 reference to coeff_abs_level_remaing
 * \param rParam                  reference to parameter
 * \param useLimitedPrefixLength
 * \param maxLog2TrDynamicRange
 */
#if RExt__DECODER_DEBUG_BIT_STATISTICS
Void TDecSbac::xReadCoefRemainExGolomb ( UInt &rSymbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange, const class TComCodingStatisticsClassType &whichStat )
#else
Void TDecSbac::xReadCoefRemainExGolomb ( UInt &rSymbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange )
#endif
{
  UInt prefix   = 0;
  UInt codeWord = 0;

  if (useLimitedPrefixLength)
  {
    const UInt longestPossiblePrefix = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange)) + COEF_REMAIN_BIN_REDUCTION;

    do
    {
      prefix++;
      m_pcTDecBinIf->decodeBinEP( codeWord RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );
    } while((codeWord != 0) && (prefix < longestPossiblePrefix));
  }
  else
  {
    do
    {
      prefix++;
      m_pcTDecBinIf->decodeBinEP( codeWord RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat) );
    } while( codeWord);
  }

  codeWord  = 1 - codeWord;
  prefix -= codeWord;
  codeWord=0;
#if VCEG_AZ07_CTX_RESIDUALCODING
  if ( prefix < g_auiGoRiceRange[rParam] )
#else
  if (prefix < COEF_REMAIN_BIN_REDUCTION )
#endif
  {
    m_pcTDecBinIf->decodeBinsEP(codeWord,rParam RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
    rSymbol = (prefix<<rParam) + codeWord;
  }
  else if (useLimitedPrefixLength)
  {
    const UInt maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange));

    const UInt prefixLength = prefix - COEF_REMAIN_BIN_REDUCTION;
    const UInt suffixLength = (prefixLength == maximumPrefixLength) ? (maxLog2TrDynamicRange - rParam) : prefixLength;

    m_pcTDecBinIf->decodeBinsEP(codeWord, (suffixLength + rParam) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));

    rSymbol = codeWord + ((((1 << prefixLength) - 1) + COEF_REMAIN_BIN_REDUCTION) << rParam);
  }
  else
  {
#if VCEG_AZ07_CTX_RESIDUALCODING
    m_pcTDecBinIf->decodeBinsEP(codeWord,prefix-g_auiGoRiceRange[rParam]+rParam RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
    rSymbol = ( ((1<<(prefix-g_auiGoRiceRange[rParam])) + g_auiGoRiceRange[rParam] -1 ) << rParam ) + codeWord;
#else
    m_pcTDecBinIf->decodeBinsEP(codeWord,prefix-COEF_REMAIN_BIN_REDUCTION+rParam RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(whichStat));
    rSymbol = (((1<<(prefix-COEF_REMAIN_BIN_REDUCTION))+COEF_REMAIN_BIN_REDUCTION-1)<<rParam)+codeWord;
#endif
  }
}

/** Parse I_PCM information.
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 *
 * If I_PCM flag indicates that the CU is I_PCM, parse its PCM alignment bits and codes.
 */
Void TDecSbac::parseIPCMInfo ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiSymbol;

  m_pcTDecBinIf->decodeBinTrm(uiSymbol);

  if (uiSymbol == 1)
  {
    Bool bIpcmFlag = true;
#if !JVET_C0024_QTBT
    const TComSPS &sps=*(pcCU->getSlice()->getSPS());
    pcCU->setPartSizeSubParts  ( SIZE_2Nx2N, uiAbsPartIdx, uiDepth );
    pcCU->setSizeSubParts      ( sps.getMaxCUWidth()>>uiDepth, sps.getMaxCUHeight()>>uiDepth, uiAbsPartIdx, uiDepth );
    pcCU->setTrIdxSubParts     ( 0, uiAbsPartIdx, uiDepth );
#endif
    pcCU->setIPCMFlagSubParts  ( bIpcmFlag, uiAbsPartIdx, uiDepth );

    const UInt minCoeffSizeY = pcCU->getPic()->getMinCUWidth() * pcCU->getPic()->getMinCUHeight();
    const UInt offsetY       = minCoeffSizeY * uiAbsPartIdx;
    for (UInt ch=0; ch < pcCU->getPic()->getNumberValidComponents(); ch++)
    {
      const ComponentID compID = ComponentID(ch);
      const UInt offset = offsetY >> (pcCU->getPic()->getComponentScaleX(compID) + pcCU->getPic()->getComponentScaleY(compID));
      Pel * pPCMSample  = pcCU->getPCMSample(compID) + offset;
      const UInt width  = pcCU->getWidth (uiAbsPartIdx) >> pcCU->getPic()->getComponentScaleX(compID);
      const UInt height = pcCU->getHeight(uiAbsPartIdx) >> pcCU->getPic()->getComponentScaleY(compID);
      const UInt sampleBits = pcCU->getSlice()->getSPS()->getPCMBitDepth(toChannelType(compID));
      for (UInt y=0; y<height; y++)
      {
        for (UInt x=0; x<width; x++)
        {
          UInt sample;
          m_pcTDecBinIf->xReadPCMCode(sampleBits, sample);
          pPCMSample[x] = sample;
        }
        pPCMSample += width;
      }
    }

    m_pcTDecBinIf->start();
  }
}

Void TDecSbac::parseCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiSymbol;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_CUTransquantBypassFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__TQ_BYPASS_FLAG) );
  pcCU->setCUTransquantBypassSubParts(uiSymbol ? true : false, uiAbsPartIdx, uiDepth);
}

/** parse skip flag
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parseSkipFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  if( pcCU->getSlice()->isIntra() )
  {
    return;
  }

  UInt uiSymbol = 0;
  UInt uiCtxSkip = pcCU->getCtxSkipFlag( uiAbsPartIdx );
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUSkipFlagSCModel.get( 0, 0, uiCtxSkip ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SKIP_FLAG) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tSkipFlag" );
  DTRACE_CABAC_T( "\tuiCtxSkip: ");
  DTRACE_CABAC_V( uiCtxSkip );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

  if( uiSymbol )
  {
    pcCU->setSkipFlagSubParts( true,        uiAbsPartIdx, uiDepth );
    pcCU->setPredModeSubParts( MODE_INTER,  uiAbsPartIdx, uiDepth );
#if !JVET_C0024_QTBT
    pcCU->setPartSizeSubParts( SIZE_2Nx2N, uiAbsPartIdx, uiDepth );
    pcCU->setSizeSubParts( pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth, pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth, uiAbsPartIdx, uiDepth );
#endif
    pcCU->setMergeFlagSubParts( true , uiAbsPartIdx, 0, uiDepth );
  }
}

#if BARC_UP_SAMPLE_RDO
Void TDecSbac::parseUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
	UInt uiSymbol = 0;
	UInt uiCtxUpCNN = pcCU->getCtxUpCNNFlag(uiAbsPartIdx);
	m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUUpCNNFlagSCModel.get(0, 0, uiCtxUpCNN) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__BARC_FLAG));

	pcCU->setUpCNNFlag(uiSymbol);
}
#endif

#if BARC
/** parse BARC flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parseBARCFlag( TComDataCU* pcCU, ChannelType eCtype )
{
  UInt uiSymbol = 0;
  UInt uiCtxBARC = pcCU->getCtxBARCFlag( 0 );
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUBARCFlagSCModel.get( 0, 0, uiCtxBARC ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__BARC_FLAG));
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tBARCFlag" );
  DTRACE_CABAC_T( "\tuiCtxBARC: ");
  DTRACE_CABAC_V( uiCtx );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

	Bool isBarc = uiSymbol;
	pcCU->setBARCFlag(eCtype, isBarc);

}
#endif

#if VCEG_AZ07_IMV
Void TDecSbac::parseiMVFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  if( pcCU->getSlice()->isIntra() )
  {
    return;
  }
#if JVET_C0024_QTBT
  else if( pcCU->getMergeFlag( uiAbsPartIdx ))
#else
  else if( pcCU->getMergeFlag( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
  {
    return;
  }
#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
  else if( pcCU->getAffineFlag( uiAbsPartIdx ) )
#else
  else if( pcCU->getAffineFlag( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
  {
    return;
  }
#endif

  UInt uiSymbol = 0;
  UInt uiCtxiMV = pcCU->getCtxiMVFlag( uiAbsPartIdx );
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUiMVFlagSCModel.get( 0, 0, uiCtxiMV ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__IMV_FLAG) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tiMVFlag" );
  DTRACE_CABAC_T( "\tuiCtxiMV: ");
  DTRACE_CABAC_V( uiCtxiMV );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

  if( uiSymbol )
  {
#if  JVET_E0076_MULTI_PEL_MVD
    m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUiMVFlagSCModel.get( 0, 0, 3 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__IMV_FLAG) );
    uiSymbol ++;
    pcCU->setiMVFlagSubParts( uiSymbol,        uiAbsPartIdx, uiDepth );
#else
    pcCU->setiMVFlagSubParts( true,        uiAbsPartIdx, uiDepth );
#endif

  }
}
#endif

#if COM16_C806_OBMC
Void TDecSbac::parseOBMCFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  pcCU->setOBMCFlagSubParts( true, uiAbsPartIdx, uiDepth );

  if ( !pcCU->getSlice()->getSPS()->getOBMC() || !pcCU->isOBMCFlagCoded( uiAbsPartIdx ) )
  {
    return;
  }

  UInt uiSymbol = 0;

  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUOBMCFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__OBMC_FLAG) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tOBMCFlag" );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

  pcCU->setOBMCFlagSubParts( uiSymbol ? true : false, uiAbsPartIdx, uiDepth );
}
#endif

#if VCEG_AZ06_IC
/** parse illumination compensation flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parseICFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{ 
  UInt uiSymbol = 0;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUICFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__IC_FLAG) );
 
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tICFlag" );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

  pcCU->setICFlagSubParts( uiSymbol ? true : false , uiAbsPartIdx, uiDepth );
}
#endif

/** parse merge flag
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \param uiPUIdx
 * \returns Void
 */
Void TDecSbac::parseMergeFlag ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx )
{
  UInt uiSymbol;
  m_pcTDecBinIf->decodeBin( uiSymbol, *m_cCUMergeFlagExtSCModel.get( 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MERGE_FLAG) );
  pcCU->setMergeFlagSubParts( uiSymbol ? true : false, uiAbsPartIdx, uiPUIdx, uiDepth );

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tMergeFlag: " );
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\tAddress: " );
  DTRACE_CABAC_V( pcCU->getCtuRsAddr() );
  DTRACE_CABAC_T( "\tuiAbsPartIdx: " );
  DTRACE_CABAC_V( uiAbsPartIdx );
  DTRACE_CABAC_T( "\n" );
}
#if VCEG_AZ05_INTRA_MPI
Void TDecSbac::parseMPIIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  if (!pcCU->getSlice()->getSPS()->getUseMPI()) 
  {
    pcCU->setMPIIdxSubParts(0, uiAbsPartIdx, uiDepth);
    return;
  }
  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
    pcCU->setMPIIdxSubParts(0, uiAbsPartIdx, uiDepth);
    return;
  }
  Int iNumberOfPassesMPI = 1;
  if (pcCU->getSlice()->getSliceType() == I_SLICE) iNumberOfPassesMPI = 2;
  else iNumberOfPassesMPI = 2;
  if (iNumberOfPassesMPI>1) // for only 1 pass no signaling is needed 
  {
    if (iNumberOfPassesMPI>2)  // 3 or 4
    {
      UInt uiSymbol0 = 0;
      UInt uiSymbol1 = 0;
      m_pcTDecBinIf->decodeBin(uiSymbol0, m_cMPIIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MPI_FLAG));
      m_pcTDecBinIf->decodeBin(uiSymbol1, m_cMPIIdxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MPI_FLAG));
      pcCU->setMPIIdxSubParts((uiSymbol0 << 1) + uiSymbol1, uiAbsPartIdx, uiDepth);
    }
    else  //iNumberOfPassesMPI==2
    {
      UInt uiSymbol = 0;
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cMPIIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MPI_FLAG));
      pcCU->setMPIIdxSubParts(uiSymbol, uiAbsPartIdx, uiDepth);
    }
  }
  else
  {
    pcCU->setMPIIdxSubParts(0, uiAbsPartIdx, uiDepth);
  }
}
#endif

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
Void TDecSbac::parsePDPCIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{

  if (!pcCU->getSlice()->getSPS()->getUsePDPC())
  {
    pcCU->setPDPCIdxSubParts(0, uiAbsPartIdx, uiDepth);
    return;
  }
  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
    pcCU->setPDPCIdxSubParts(0, uiAbsPartIdx, uiDepth);
    return;
  }
  Int iNumberOfPassesPDPC = 1;
  if (pcCU->getSlice()->getSliceType() == I_SLICE) iNumberOfPassesPDPC = 2;
  else iNumberOfPassesPDPC = 2;
  if (iNumberOfPassesPDPC > 1) // for only 1 pass no signaling is needed 
  {
    if (iNumberOfPassesPDPC > 2)  // 3 or 4
    {
      UInt uiSymbol0 = 0;
      UInt uiSymbol1 = 0;
      m_pcTDecBinIf->decodeBin(uiSymbol0, m_cPDPCIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__PDPC_FLAG));
      m_pcTDecBinIf->decodeBin(uiSymbol1, m_cPDPCIdxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__PDPC_FLAG));
      pcCU->setPDPCIdxSubParts((uiSymbol0 << 1) + uiSymbol1, uiAbsPartIdx, uiDepth);
     }
    else  //iNumberOfPassesMPI==2
    {
      UInt uiSymbol = 0;
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cPDPCIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__PDPC_FLAG));
      pcCU->setPDPCIdxSubParts(uiSymbol, uiAbsPartIdx, uiDepth);
     }
   }
  else
  {
    pcCU->setPDPCIdxSubParts(0, uiAbsPartIdx, uiDepth);
  }
}
#endif

#if DIMD_INTRA_PRED
Void TDecSbac::parseDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
#if TT_TREE_DIMD // deriveZIdxInCTU	
  UInt uiRaster = g_auiZscanToRaster[pcCU->deriveZIdxInCTU(pcCU->getZorderIdxInCtu(), uiAbsPartIdx)];
#else
  UInt uiShort, uiLong, uiStride;
  UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiAbsPartIdx];
  //UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiAbsPartIdx];
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = pcCU->getPic()->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  UInt uiShortDepth = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = (pcCU->getPic()->getNumPartitionsInCtu() >> (uiShortDepth << 1));
  UInt uiNumPartInShort = (pcCU->getPic()->getNumPartInCtuWidth() >> uiShortDepth);
#endif
  //assert(isLuma(pcCU->getTextType()));

  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
#if TT_TREE_DIMD		
    pcCU->setDIMDEnabledFlagSubParts_TT(CHANNEL_TYPE_LUMA, 0, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
#else
    for (UInt i = 0; i < uiLong; i += uiShort)
    {
      memset(pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA) + g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu(), 0, uiCurrPartNumb);
      uiRaster += uiNumPartInShort * uiStride;
    }

#endif
    return;
  }

  UInt uiSymbol = 0;
  m_pcTDecBinIf->decodeBin(uiSymbol, m_cDIMDEnabledSCModel.get(0, 0, uiDepth) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DIMD_NO_BT_FLAG));

#if TT_TREE_DIMD	
  pcCU->setDIMDEnabledFlagSubParts_TT(CHANNEL_TYPE_LUMA, uiSymbol, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
#else
  for (UInt i = 0; i < uiLong; i += uiShort)
  {
    memset(pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA) + g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu(), uiSymbol, uiCurrPartNumb);
    uiRaster += uiNumPartInShort * uiStride;
  }
#endif
}
Void TDecSbac::parseDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
  assert(pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx));

#if TT_TREE_DIMD // deriveZIdxInCTU	
  UInt uiRaster = g_auiZscanToRaster[pcCU->deriveZIdxInCTU(pcCU->getZorderIdxInCtu(), uiAbsPartIdx)];
#else
  UInt uiShort, uiLong, uiStride;
  UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiAbsPartIdx];
  //UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiAbsPartIdx];
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = pcCU->getPic()->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  UInt uiShortDepth = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = (pcCU->getPic()->getNumPartitionsInCtu() >> (uiShortDepth << 1));
  UInt uiNumPartInShort = (pcCU->getPic()->getNumPartInCtuWidth() >> uiShortDepth);
#endif

  UInt uiSymbol = 0;
  m_pcTDecBinIf->decodeBin(uiSymbol, m_cDIMDNoBTFlagSCModel.get(0, 0, uiDepth) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DIMD_FLAG));

#if TT_TREE_DIMD	
  pcCU->setDIMDNoBTLevelFlagSubParts_TT(CHANNEL_TYPE_LUMA, uiSymbol, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
#else

  for (UInt i = 0; i < uiLong; i += uiShort)
  {
    memset(pcCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA) + g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu(), uiSymbol, uiCurrPartNumb);
    uiRaster += uiNumPartInShort * uiStride;
  }
#endif
}
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
Void TDecSbac::parseROTIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{ 
#if COM16_C1044_NSST
  if (!pcCU->getSlice()->getSPS()->getUseNSST()) 
#else
  if (!pcCU->getSlice()->getSPS()->getUseROT()) 
#endif
  {
#if JVET_C0024_QTBT
    pcCU->setROTIdxSubParts(CHANNEL_TYPE_LUMA, 0, uiAbsPartIdx, uiDepth);
#else
    pcCU->setROTIdxSubParts(0, uiAbsPartIdx, uiDepth);
#endif
    return;
  }
  Int iNumberOfPassesROT = 1;
  if( pcCU->isIntra(uiAbsPartIdx)
#if VCEG_AZ05_INTRA_MPI
    && pcCU->getMPIIdx(uiAbsPartIdx) ==0
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    && pcCU->getPDPCIdx(uiAbsPartIdx) == 0
#endif
    && !pcCU->getCUTransquantBypass(uiAbsPartIdx)
#if MLIP_EXTEND_REF_LINE
      && (pcCU->getExternRef(CHANNEL_TYPE_LUMA, uiAbsPartIdx) == 0)
#endif
    )  iNumberOfPassesROT = 4;

#if COM16_C1044_NSST
#if JVET_C0024_QTBT
#if RRC_DWDIP_MLR
  if (pcCU->getDWDIP_GPIdx(uiAbsPartIdx)==0)
#endif
  if( iNumberOfPassesROT==4)
#else
  if( iNumberOfPassesROT==4 && pcCU->getPartitionSize(uiAbsPartIdx)==SIZE_2Nx2N )
#endif
  {
#if DIMD_INTRA_PRED
    iNumberOfPassesROT = pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx) ? 4 : (
#if FAST_NSST_v1
		pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) <= DC_IDX ? 2 : 4
#else
		pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) <= DC_IDX ? 3 : 4
#endif
		);
#else
#if FAST_NSST_v1
	  iNumberOfPassesROT = pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) <= DC_IDX ? 2 : 4;
#else
    iNumberOfPassesROT = pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) <= DC_IDX ? 3 : 4;
#endif
#endif
  }
#if FAST_NSST_v1
  if (iNumberOfPassesROT == 2)
  {
	  UInt uiSymbol = 0;
	  m_pcTDecBinIf->decodeBin(uiSymbol, m_cROTidxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
	  pcCU->setROTIdxSubParts(CHANNEL_TYPE_LUMA, uiSymbol, uiAbsPartIdx, uiDepth);
	  return;
  }
#endif
  if( iNumberOfPassesROT==3 )
  {
    UInt uiSymbol = 0;
#if JVET_C0042_UNIFIED_BINARIZATION
    m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
    if( uiSymbol )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 3 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
      uiSymbol ++;
    }
#else
    m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
    if( uiSymbol )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
      uiSymbol ++;
    }
#endif
#if JVET_C0024_QTBT
    pcCU->setROTIdxSubParts( CHANNEL_TYPE_LUMA, uiSymbol, uiAbsPartIdx,  uiDepth ); 
#else
    pcCU->setROTIdxSubParts( uiSymbol, uiAbsPartIdx,  uiDepth ); 
#endif
  }
  else
#endif
  if (iNumberOfPassesROT>1) // for only 1 pass no signaling is needed 
  {
#if JVET_C0042_UNIFIED_BINARIZATION
 UInt uiNSST = 0, idx_ROT=0;
    m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
   idx_ROT +=uiNSST;
    if( uiNSST )
    {
          m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
          idx_ROT +=uiNSST;

           if(uiNSST )
            {
              m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 4) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
              idx_ROT +=uiNSST;
            }
    }
#if JVET_C0024_QTBT
    pcCU->setROTIdxSubParts( CHANNEL_TYPE_LUMA, idx_ROT, uiAbsPartIdx,  uiDepth ); 
#else
    pcCU->setROTIdxSubParts( idx_ROT, uiAbsPartIdx,  uiDepth ); 
#endif
#else
      UInt uiSymbol0 = 0;
      UInt uiSymbol1 = 0;
      m_pcTDecBinIf->decodeBin( uiSymbol0, m_cROTidxSCModel.get(0,0, uiDepth ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
      m_pcTDecBinIf->decodeBin( uiSymbol1, m_cROTidxSCModel.get(0,0, uiDepth ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
#if JVET_C0024_QTBT
      pcCU->setROTIdxSubParts( CHANNEL_TYPE_LUMA, (uiSymbol0<<1) +uiSymbol1, uiAbsPartIdx,  uiDepth ); //printf ("%d ",(uiSymbol0<<1) +uiSymbol1);
#else
      pcCU->setROTIdxSubParts( (uiSymbol0<<1) +uiSymbol1, uiAbsPartIdx,  uiDepth ); //printf ("%d ",(uiSymbol0<<1) +uiSymbol1);
#endif
#endif
  }
  else
  {
#if JVET_C0024_QTBT
    pcCU->setROTIdxSubParts( CHANNEL_TYPE_LUMA, 0, uiAbsPartIdx,  uiDepth ); 
#else
       pcCU->setROTIdxSubParts( 0, uiAbsPartIdx,  uiDepth ); 
#endif
  }


}

#if JVET_C0024_QTBT
Void TDecSbac::parseROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{ 
#if DIMD_INTRA_PRED
  assert(!pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
#endif
#if COM16_C1044_NSST
  if (!pcCU->getSlice()->getSPS()->getUseNSST()) 
#else
  if (!pcCU->getSlice()->getSPS()->getUseROT()) 
#endif
  {
    pcCU->setROTIdxSubParts(CHANNEL_TYPE_CHROMA, 0, uiAbsPartIdx, uiDepth);
    return;
  }
  Int iNumberOfPassesROT = 1;
  if( pcCU->isIntra(uiAbsPartIdx)
    && !pcCU->getCUTransquantBypass(uiAbsPartIdx)
    )  iNumberOfPassesROT = 4;

#if COM16_C1044_NSST
  if( iNumberOfPassesROT==4)
  {
    UInt uiIntraMode = pcCU->getIntraDir( CHANNEL_TYPE_CHROMA, uiAbsPartIdx );
#if JVET_E0062_MULTI_DMS && COM16_C806_LMCHROMA
    if( uiIntraMode == LM_CHROMA_IDX )
    {
      uiIntraMode = PLANAR_IDX;
    }
#if JVET_E0077_ENHANCED_LM
    else if (IsLMMode(uiIntraMode))
    {
        uiIntraMode = PLANAR_IDX;
    }
#endif
#else
    if( uiIntraMode == DM_CHROMA_IDX )
    {
#if DIMD_INTRA_PRED
      if (pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx))
      {
        uiIntraMode = DC_IDX + 1;
      }
      else
      {
#endif
        uiIntraMode = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
#if DIMD_INTRA_PRED
      }
#endif
    }
#if COM16_C806_LMCHROMA
    else if( uiIntraMode == LM_CHROMA_IDX )
    {
      uiIntraMode = PLANAR_IDX;
    }
#endif

#if JVET_E0077_ENHANCED_LM
    else if (IsLMMode(uiIntraMode))
    {
        uiIntraMode = PLANAR_IDX;
    }
#endif

#endif
#if FAST_NSST_v1
	iNumberOfPassesROT = uiIntraMode <= DC_IDX ? 2 : 4;
#else
    iNumberOfPassesROT = uiIntraMode <= DC_IDX ? 3 : 4;
#endif
  }
#if FAST_NSST_v1 
  if (iNumberOfPassesROT == 2)
  {
	  UInt uiSymbol = 0;
	  m_pcTDecBinIf->decodeBin(uiSymbol, m_cROTidxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
	  pcCU->setROTIdxSubParts(CHANNEL_TYPE_CHROMA, uiSymbol, uiAbsPartIdx, uiDepth);
	  return;
  }
#endif
  if( iNumberOfPassesROT==3 )
  {
    UInt uiSymbol = 0;
#if JVET_C0042_UNIFIED_BINARIZATION
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
    if( uiSymbol )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 3 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
      uiSymbol ++;
    }
#else
    m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
    if( uiSymbol )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cROTidxSCModel.get(0,0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
      uiSymbol ++;
    }
#endif
    pcCU->setROTIdxSubParts( CHANNEL_TYPE_CHROMA, uiSymbol, uiAbsPartIdx,  uiDepth ); 
  }
  else
#endif
  if (iNumberOfPassesROT>1) // for only 1 pass no signaling is needed 
  {
#if JVET_C0042_UNIFIED_BINARIZATION
UInt uiNSST = 0, idx_ROT=0;
    m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
   idx_ROT +=uiNSST;
    if( uiNSST )
    {
          m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
          idx_ROT +=uiNSST;

           if(uiNSST )
            {
              m_pcTDecBinIf->decodeBin( uiNSST, m_cROTidxSCModel.get(0,0, 4) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG) );
              idx_ROT +=uiNSST;
            }
    }
      pcCU->setROTIdxSubParts( CHANNEL_TYPE_CHROMA, idx_ROT, uiAbsPartIdx,  uiDepth ); 
#else
      UInt uiSymbol0 = 0;
      UInt uiSymbol1 = 0;
      m_pcTDecBinIf->decodeBin( uiSymbol0, m_cROTidxSCModel.get(0,0, uiDepth ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
      m_pcTDecBinIf->decodeBin( uiSymbol1, m_cROTidxSCModel.get(0,0, uiDepth ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ROT_FLAG));
      pcCU->setROTIdxSubParts( CHANNEL_TYPE_CHROMA, (uiSymbol0<<1) +uiSymbol1, uiAbsPartIdx,  uiDepth ); //printf ("%d ",(uiSymbol0<<1) +uiSymbol1);
#endif
  }
  else
  {
    pcCU->setROTIdxSubParts( CHANNEL_TYPE_CHROMA, 0, uiAbsPartIdx,  uiDepth ); 
  }
}
#endif
#endif
Void TDecSbac::parseMergeIndex ( TComDataCU* pcCU, UInt& ruiMergeIndex )
{
  UInt uiUnaryIdx = 0;
  UInt uiNumCand = pcCU->getSlice()->getMaxNumMergeCand();
  if ( uiNumCand > 1 )
  {
    for( ; uiUnaryIdx < uiNumCand - 1; ++uiUnaryIdx )
    {
      UInt uiSymbol = 0;
#if COM16_C806_GEN_MRG_IMPROVEMENT
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, (uiUnaryIdx>NUM_MERGE_IDX_EXT_CTX-1? NUM_MERGE_IDX_EXT_CTX-1:uiUnaryIdx) ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MERGE_INDEX));
#else
      if ( uiUnaryIdx==0 )
      {
        m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MERGE_INDEX) );
      }
      else
      {
        m_pcTDecBinIf->decodeBinEP( uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MERGE_INDEX) );
      }
#endif
      if( uiSymbol == 0 )
      {
        break;
      }
    }
  }
  ruiMergeIndex = uiUnaryIdx;

  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseMergeIndex()" )
  DTRACE_CABAC_T( "\tuiMRGIdx= " )
  DTRACE_CABAC_V( ruiMergeIndex )
  DTRACE_CABAC_T( "\n" )
}

#if VCEG_AZ07_FRUC_MERGE
Void TDecSbac::parseFRUCMgrMode ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx )
{
  if( !pcCU->getSlice()->getSPS()->getUseFRUCMgrMode() )
    return;

  UChar uhFRUCMode = FRUC_MERGE_OFF;
  UInt uiFirstBin = 0;
  UInt uiSecondBin = 0;
#if THIRD_FRUC_MODE
  UInt uiThirdBin = 0;
#endif
#if THIRD_FRUC_MODE_DISABLE_BILATERALFORTID0
  if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0) > 1 && pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 1 )
  {
#endif
#if THIRD_FRUC_MODE || FIX_DECODER_STATISTICS
      m_pcTDecBinIf->decodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_IDX));
#else
      m_pcTDecBinIf->decodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_MERGE_IDX));
#endif
      if (uiFirstBin)
      {
          if (pcCU->getSlice()->isInterP())
          {
              uhFRUCMode = FRUC_MERGE_TEMPLATE;
          }
          else
          {
              m_pcTDecBinIf->decodeBin(uiSecondBin, m_cCUFRUCMESCModel.get(0, 0, pcCU->getCtxFRUCME(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_MERGE_IDX));
              uhFRUCMode = uiSecondBin ? FRUC_MERGE_BILATERALMV : FRUC_MERGE_TEMPLATE;
#if THIRD_FRUC_MODE 
              if (uhFRUCMode == FRUC_MERGE_TEMPLATE)
              {
                  m_pcTDecBinIf->decodeBin(uiThirdBin, m_cCUFRUCMESCModel_uni.get(0, 0, pcCU->getCtxFRUCME_UNI(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_UNI_IDX));
                  uhFRUCMode = uiThirdBin ? FRUC_MERGE_TEMPLATE : FRUC_MERGE_TEMPLATE_UNI;
              }
#endif
          }
      }
#if THIRD_FRUC_MODE_DISABLE_BILATERALFORTID0
  }
  else
  {
#if THIRD_FRUC_MODE || FIX_DECODER_STATISTICS
      m_pcTDecBinIf->decodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_IDX));
      uhFRUCMode = uiFirstBin ? FRUC_MERGE_TEMPLATE : FRUC_MERGE_OFF;
#else
      m_pcTDecBinIf->decodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_MERGE_IDX));
#endif
      if (uiFirstBin)
      {

#if THIRD_FRUC_MODE 
            m_pcTDecBinIf->decodeBin(uiThirdBin, m_cCUFRUCMESCModel_uni.get(0, 0, pcCU->getCtxFRUCME_UNI(uiAbsPartIdx)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__FRUC_UNI_IDX));
            uhFRUCMode = uiThirdBin ? FRUC_MERGE_TEMPLATE : FRUC_MERGE_TEMPLATE_UNI;
#endif
      }
  }
#endif
  pcCU->setFRUCMgrModeSubParts( uhFRUCMode , uiAbsPartIdx, uiPUIdx, uiDepth );
}
#endif

Void TDecSbac::parseMVPIdx      ( Int& riMVPIdx )
{
  UInt uiSymbol;
  xReadUnaryMaxSymbol(uiSymbol, m_cMVPIdxSCModel.get(0), 1, AMVP_MAX_NUM_CANDS-1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVP_IDX) );
  riMVPIdx = uiSymbol;
}

#if JVET_C0024_QTBT
Void TDecSbac::parseBTSplitMode   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
  , UInt uiBTDepth
#endif
  )
{
  UInt uiSymbol;
#if TT_TREE_DS1
  UInt uiSymbol2 = 0;
#endif

#if RExt__DECODER_DEBUG_BIT_STATISTICS
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__BTSPLIT_FLAG, g_aucConvertToBit[uiWidth]+2);
#endif
#if TT_SYNTAX_NEW
  char cAllowSplit = pcCU->deriveAllowedSplitType(uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
  int iHBT = (cAllowSplit>>3) & 1; 
  int iVBT = (cAllowSplit>>2) & 1; 
  int iHTT = (cAllowSplit>>1) & 1;
  int iVTT = (cAllowSplit   ) & 1;
  int sum  = iHBT + iVBT + iHTT + iVTT;

  UInt uiCtx;
  if(sum!=0)
  {
    uiCtx = pcCU->getCtxBTSplitFlag(uiAbsPartIdx, uiWidth, uiHeight);
    m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, uiCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
  }
  else
    uiSymbol = 0;
#else
  UInt uiCtx = pcCU->getCtxBTSplitFlag(uiAbsPartIdx, uiWidth, uiHeight);

  m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, uiCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
#endif
  UInt uiMinBTSize = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType())?MIN_BT_SIZE:MIN_BT_SIZE_C): MIN_BT_SIZE_INTER;

#if TT_SYNTAX //dir ctx
  UInt uiBTSCtx = uiWidth==uiHeight ? 0: (uiWidth>uiHeight? 1: 2);
  uiBTSCtx = (uiWidth>(uiHeight<<1))? 3 : ((uiHeight>(uiWidth<<1))? 4 : uiBTSCtx);
#if !BT_FORCE_BF
  assert(uiBTDepth<=2);
#endif
#if BT_FORCE_FBTD
  UInt uiFBTDDiff = uiBTDepth - pcCU->getFBTDepth(uiAbsPartIdx);
  UInt uiBTSCtx2 = uiFBTDDiff >2 ? 2 : uiFBTDDiff;
#elif BT_FORCE_CTX
  UInt uiBTSCtx2 = uiBTDepth > 2 ? 2 : uiBTDepth;
#elif BT_FORCE_CTX0
  UInt uiBTSCtx2 = uiBTDepth;
  if (pcCU->getBoundaryCtuLocation())
    uiBTSCtx2 = 0;
#else
  UInt uiBTSCtx2 = uiBTDepth;
#endif
  Bool bTreeTypeNeed = true;
#endif
#if TT_LIMIT_SYNTAX
  UInt uiParentType = (uiBTDepth==0)? 0 : pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth-1);
#endif

#if TT_SYNTAX_NEW
  if(uiSymbol != 0)
  {   
    if(sum==4)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      uiSymbol++;
    }
    else if(sum==3)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      if(iHBT==0 || iHTT==0)
      {
        if(uiSymbol==1)
          m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        else
          uiSymbol2 = (iHBT==0);
      }
      else if(iVBT==0 || iVTT==0)
      {
        if(uiSymbol==0)
          m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        else
          uiSymbol2 = (iVBT==0);
      }
      uiSymbol++;
    }
    else if(sum==2)
    {
      if(iHTT==0 && iVTT==0)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        uiSymbol2 = 0;
      }
      else if(iHBT==0 && iHTT==0)
      {
        uiSymbol = 1;
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      }
      else if(iVBT==0 && iVTT==0)
      {
        uiSymbol = 0;
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      }
      else if (iVBT==0 && iHTT==0)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        uiSymbol2 = (uiSymbol==0)? 0 : 1;
      }
      else if (iHBT==0 && iVTT==0)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        uiSymbol2 = (uiSymbol==1)? 0 : 1;
      }
      else
        assert(0);
      uiSymbol++;
    }
    else//i.e., if(sum==1)
    {
      uiSymbol = (iVBT==1 || iVTT==1);
      uiSymbol2= (iHTT==1 || iVTT==1);
      uiSymbol++;
    }
  }
#else
  if (uiWidth==uiMinBTSize && uiSymbol!=0)
  {
#if TT_SYNTAX
    if(uiHeight>=(uiMinBTSize<<2))
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_HORI_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_HORI_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-W" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSymbol2 );
        DTRACE_CABAC_T( "\n");
      }
    }
    else
      uiSymbol2 = 0;
#endif
    uiSymbol = 1;
  }
  else if (uiHeight==uiMinBTSize && uiSymbol!=0)
  {
#if TT_SYNTAX
    if(uiWidth>=(uiMinBTSize<<2))
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_VERT_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_VERT_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-H" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSymbol2 );
        DTRACE_CABAC_T( "\n");
      }
    }
    else
      uiSymbol2 = 0;
#endif
    uiSymbol = 2;
  }
#if TT_SYNTAX
  else if (uiWidth==(uiMinBTSize<<1) && uiHeight==(uiMinBTSize<<1) && uiSymbol!=0)
  {
    uiSymbol2 = 0; //can only use binary split
    m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
    DTRACE_CABAC_VL( g_nSymbolCounter++ );
    DTRACE_CABAC_T( "\tSTDirFlag-W2H2" );
    DTRACE_CABAC_T( "\tuiBTSCtx: ");
    DTRACE_CABAC_V( uiBTSCtx );
    DTRACE_CABAC_T( "\tuiSymbol: ");
    DTRACE_CABAC_V( uiSymbol );
    DTRACE_CABAC_T( "\n");
    uiSymbol ++;
  }
  else if (uiWidth==(uiMinBTSize<<1) && uiSymbol!=0)
  {
    assert(uiHeight>(uiMinBTSize<<1));
    m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
    DTRACE_CABAC_VL( g_nSymbolCounter++ );
    DTRACE_CABAC_T( "\tSTDirFlag-W2H3" );
    DTRACE_CABAC_T( "\tuiBTSCtx: ");
    DTRACE_CABAC_V( uiBTSCtx );
    DTRACE_CABAC_T( "\tuiSymbol: ");
    DTRACE_CABAC_V( uiSymbol );
    DTRACE_CABAC_T( "\n");
    uiSymbol ++;

    if(uiSymbol==2) //if vert split
      uiSymbol2 = 0;
    else
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_HORI_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_HORI_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-W2H3" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSymbol2 );
        DTRACE_CABAC_T( "\n");
      }
    }
  }
  else if (uiHeight==(uiMinBTSize<<1) && uiSymbol!=0)
  {
    assert(uiWidth>(uiMinBTSize<<1));
    m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
    DTRACE_CABAC_VL( g_nSymbolCounter++ );
    DTRACE_CABAC_T( "\tSTDirFlag-W3H2" );
    DTRACE_CABAC_T( "\tuiBTSCtx: ");
    DTRACE_CABAC_V( uiBTSCtx );
    DTRACE_CABAC_T( "\tuiSymbol: ");
    DTRACE_CABAC_V( uiSymbol );
    DTRACE_CABAC_T( "\n");
    uiSymbol ++;

    if(uiSymbol==1) //if hori split
      uiSymbol2 = 0;
    else
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_VERT_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_VERT_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-W3H2" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSymbol2 );
        DTRACE_CABAC_T( "\n");
      }
    }
  }
#endif
#if (TT_LIMIT_AR>0)
  else if (uiWidth>=uiHeight*TT_LIMIT_AR && uiSymbol!=0)
  {
    uiSymbol = 2;
#if TT_LIMIT_SYNTAX
    if(TT_LIMIT_SMALL && uiParentType==ST_VERT_TT_S)
      bTreeTypeNeed = false;
    if(TT_LIMIT_LARGE && uiParentType==ST_VERT_TT_L)
      bTreeTypeNeed = false;
#endif
    if(bTreeTypeNeed)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTModeFlag-ARBTW" );
      DTRACE_CABAC_T( "\tuiBTSCtx2: ");
      DTRACE_CABAC_V( uiBTSCtx2 );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol2 );
      DTRACE_CABAC_T( "\n");
    }
  }
  else if (uiHeight>=uiWidth*TT_LIMIT_AR && uiSymbol!=0)
  {
    uiSymbol = 1;
#if TT_LIMIT_SYNTAX
    if(TT_LIMIT_SMALL && uiParentType==ST_HORI_TT_S)
      bTreeTypeNeed = false;
    if(TT_LIMIT_LARGE && uiParentType==ST_HORI_TT_L)
      bTreeTypeNeed = false;
#endif
    if(bTreeTypeNeed)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTModeFlag-ARBTH" );
      DTRACE_CABAC_T( "\tuiBTSCtx2: ");
      DTRACE_CABAC_V( uiBTSCtx2 );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol2 );
      DTRACE_CABAC_T( "\n");
    }
  }
#endif
#if JVET_C0024_BT_RMV_REDUNDANT
  else if ( pcCU->getSplitConstrain() == 1 && uiSymbol!=0 ) 
  {
#if TT_SYNTAX
    if(bTreeTypeNeed)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTModeFlag-C1" );
      DTRACE_CABAC_T( "\tuiBTSCtx2: ");
      DTRACE_CABAC_V( uiBTSCtx2 );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol2 );
      DTRACE_CABAC_T( "\n");
    }

    if(uiSymbol2)
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
        uiSymbol = (uiParentType==ST_HORI_TT_S)? 2 : 1;
      else if(TT_LIMIT_LARGE && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
        uiSymbol = (uiParentType==ST_HORI_TT_L)? 2 : 1;
      else
      {
#endif
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag-C1" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol );
      DTRACE_CABAC_T( "\n");
      uiSymbol ++;
#if TT_LIMIT_SYNTAX
      }
#endif
    }
    else
#endif
    uiSymbol = 2;
  }
  else if ( pcCU->getSplitConstrain() == 2 && uiSymbol!=0 )
  {
#if TT_SYNTAX
    if(bTreeTypeNeed)
    {
      m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTModeFlag-C2" );
      DTRACE_CABAC_T( "\tuiBTSCtx2: ");
      DTRACE_CABAC_V( uiBTSCtx2 );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol2 );
      DTRACE_CABAC_T( "\n");
    }
    
    if(uiSymbol2)
    {
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
        uiSymbol = (uiParentType==ST_HORI_TT_S)? 2 : 1;
      else if(TT_LIMIT_LARGE && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
        uiSymbol = (uiParentType==ST_HORI_TT_L)? 2 : 1;
      else
      {
#endif
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag-C2" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol );
      DTRACE_CABAC_T( "\n");
      uiSymbol ++;
    }
#if TT_LIMIT_SYNTAX
    }
#endif
    else
#endif
    uiSymbol = 1;
  }
#endif
  else
  {
    if (uiSymbol != 0)
    {
#if TT_SYNTAX
      if(bTreeTypeNeed)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol2, m_cBTSplitFlagSCModel.get( 0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSymbol2 );
        DTRACE_CABAC_T( "\n");
      }
#endif
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiSymbol2==1 && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
        uiSymbol = (uiParentType==ST_HORI_TT_S)? 2 : 1;
      else if(TT_LIMIT_LARGE && uiSymbol2==1 && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
        uiSymbol = (uiParentType==ST_HORI_TT_L)? 2 : 1;
      else
      {
#endif
      UInt uiBTSCtx = uiWidth==uiHeight ? 0: (uiWidth>uiHeight? 1: 2);
#if TT_SYNTAX //dir ctx
      uiBTSCtx = (uiWidth>(uiHeight<<1))? 3 : ((uiHeight>(uiWidth<<1))? 4 : uiBTSCtx);
#endif
      m_pcTDecBinIf->decodeBin(uiSymbol, m_cBTSplitFlagSCModel.get( 0, 0, 3+uiBTSCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSymbol );
      DTRACE_CABAC_T( "\n");
      uiSymbol ++;
#if TT_LIMIT_SYNTAX
      }
#endif
    }
  }
#endif

#if TT_TREE_DS1
  pcCU->setBTSplitModeSubParts((uiSymbol<<2)+(uiSymbol2<<1), uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth);
#else
  pcCU->setBTSplitModeSubParts(uiSymbol, uiAbsPartIdx, uiWidth, uiHeight);
#endif
}
#endif

Void TDecSbac::parseSplitFlag     ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  if ( uiDepth == g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] 
  - g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMinQTSize(pcCU->getSlice()->getSliceType(), pcCU->getTextType())])
#else
  if( uiDepth == pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() )
#endif
  {
    pcCU->setDepthSubParts( uiDepth, uiAbsPartIdx );
    return;
  }
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__SPLIT_FLAG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth]+2);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__SPLIT_FLAG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth]+2);
#endif
#endif

  UInt uiSymbol;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUSplitFlagSCModel.get( 0, 0, pcCU->getCtxSplitFlag( uiAbsPartIdx, uiDepth ) ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tSplitFlag\n" )
  pcCU->setDepthSubParts( uiDepth + uiSymbol, uiAbsPartIdx );

  return;
}

#if !JVET_C0024_QTBT
/** parse partition size
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parsePartSize( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiSymbol, uiMode = 0;
  PartSize eMode;
#if COM16_C806_LARGE_CTU
  const UShort cuWidth =UShort(pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth);
  const UShort cuHeight=UShort(pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth);
#else
  const UChar cuWidth =UChar(pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth);
  const UChar cuHeight=UChar(pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth);
#endif
  const Int log2DiffMaxMinCodingBlockSize = pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize();

#if RExt__DECODER_DEBUG_BIT_STATISTICS
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__PART_SIZE, g_aucConvertToBit[cuWidth]+2);
#endif

  assert ( pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() == log2DiffMaxMinCodingBlockSize);
  if ( pcCU->isIntra( uiAbsPartIdx ) )
  {
    uiSymbol = 1;
    if( uiDepth == log2DiffMaxMinCodingBlockSize )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUPartSizeSCModel.get( 0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }
    eMode = uiSymbol ? SIZE_2Nx2N : SIZE_NxN;
    UInt uiTrLevel = 0;
    UInt uiWidthInBit  = g_aucConvertToBit[pcCU->getWidth(uiAbsPartIdx)]+2;
    UInt uiTrSizeInBit = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxTrSize()]+2;
    uiTrLevel          = uiWidthInBit >= uiTrSizeInBit ? uiWidthInBit - uiTrSizeInBit : 0;
    if( eMode == SIZE_NxN )
    {
      pcCU->setTrIdxSubParts( 1+uiTrLevel, uiAbsPartIdx, uiDepth );
    }
    else
    {
      pcCU->setTrIdxSubParts( uiTrLevel, uiAbsPartIdx, uiDepth );
    }
  }
  else
  {
    UInt uiMaxNumBits = 2;
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL && !COM16_C806_DISABLE_4X4_PU
    if ( (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag() && (uiDepth == log2DiffMaxMinCodingBlockSize)) ||
        (!pcCU->getSlice()->getSPS()->getAtmvpEnableFlag() && (uiDepth == log2DiffMaxMinCodingBlockSize && !( (g_uiMaxCUWidth>>uiDepth) == 8 && (g_uiMaxCUHeight>>uiDepth) == 8 ) ) ))
#else
    if( uiDepth == log2DiffMaxMinCodingBlockSize && !( cuWidth == 8 && cuHeight == 8 ) )
#endif
    {
      uiMaxNumBits ++;
    }

    for ( UInt ui = 0; ui < uiMaxNumBits; ui++ )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUPartSizeSCModel.get( 0, 0, ui) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
      if ( uiSymbol )
      {
        break;
      }
      uiMode++;
    }
    eMode = (PartSize) uiMode;
    if ( pcCU->getSlice()->getSPS()->getUseAMP() && uiDepth < log2DiffMaxMinCodingBlockSize )
    {
      if (eMode == SIZE_2NxN)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUPartSizeSCModel.get( 0, 0, 3 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
        if (uiSymbol == 0)
        {
          m_pcTDecBinIf->decodeBinEP(uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
          eMode = (uiSymbol == 0? SIZE_2NxnU : SIZE_2NxnD);
        }
      }
      else if (eMode == SIZE_Nx2N)
      {
        m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUPartSizeSCModel.get( 0, 0, 3 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
        if (uiSymbol == 0)
        {
          m_pcTDecBinIf->decodeBinEP(uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
          eMode = (uiSymbol == 0? SIZE_nLx2N : SIZE_nRx2N);
        }
      }
    }
  }
  pcCU->setPartSizeSubParts( eMode, uiAbsPartIdx, uiDepth );
  pcCU->setSizeSubParts( cuWidth, cuHeight, uiAbsPartIdx, uiDepth );
}
#endif


/** parse prediction mode
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parsePredMode( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  if( pcCU->getSlice()->isIntra() )
  {
    pcCU->setPredModeSubParts( MODE_INTRA, uiAbsPartIdx, uiDepth );
    return;
  }

  UInt uiSymbol;
  Int  iPredMode = MODE_INTER;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUPredModeSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__PRED_MODE) );
  iPredMode += uiSymbol;
  pcCU->setPredModeSubParts( (PredMode)iPredMode, uiAbsPartIdx, uiDepth );
}


Void TDecSbac::parseIntraDirLumaAng  ( TComDataCU* pcCU, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  UInt partNum = 1;
#else
  PartSize mode = pcCU->getPartitionSize( absPartIdx );
  UInt partNum = mode==SIZE_NxN?4:1;
#endif
  UInt partOffset = ( pcCU->getPic()->getNumPartitionsInCtu() >> ( pcCU->getDepth(absPartIdx) << 1 ) ) >> 2;
  UInt mpmPred[4],symbol;
  Int j,intraPredMode;

#if JVET_C0055_INTRA_MPM
  static const UInt mpmContext[NUM_INTRA_MODE] = { 1, 1, 
#if VCEG_AZ07_INTRA_65ANG_MODES
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,   // 2-34
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3       // 35-67
#else
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,    // 2-18
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3     // 19-35
#endif
  };
#elif VCEG_AZ07_INTRA_65ANG_MODES
  const UInt uiContextMPM0[4] = { 2, 3, 1, 2 };
  const UInt uiContextMPM1[4] = { 4, 5, 5, 6 };
  const UInt uiContextMPM2[4] = { 7, 7, 8, 7 };
#endif


#if MLIP_EXTEND_REF_LINE

  Bool bSPSUseMLIP = pcCU->getSlice()->getSPS()->getUseMLIP();// sps MLIP flag

#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize() >> depth] + 2, CHANNEL_TYPE_LUMA);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth() >> depth] + 2, CHANNEL_TYPE_LUMA);
#endif
#endif
  UInt refflag = 0;
  if(bSPSUseMLIP)
  {

#if MLIP_IMPROVED_INDEX_CODING
      Int iCurNumRef = MLIP_REF_LINE_NUM;
      UInt tempSymbol;
#if MLIP_REF_LINE_NUM > 2
      m_pcTDecBinIf->decodeBin(tempSymbol, m_cCUExternRefSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      if (tempSymbol) {
          m_pcTDecBinIf->decodeBin(tempSymbol, m_cCUExternRefSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
#if MLIP_REF_LINE_NUM > 3
          if (tempSymbol) {
              m_pcTDecBinIf->decodeBin(tempSymbol, m_cCUExternRefSCModel.get(0, 0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
              tempSymbol++;
          }
#endif
          tempSymbol++;
      }
      //MX: swap the coding code of tempSymbol 3 and tempSymbol 1;need to be modified
      if (tempSymbol == MLIP_MODI_INDEX)
      {
          tempSymbol = 1;
      }
      else if (tempSymbol == 1)
      {
          tempSymbol = MLIP_MODI_INDEX;
      }
      refflag = tempSymbol;
      assert(iCurNumRef == 4);//MX: if not NumLine ==4, need to be modified
#else
      m_pcTDecBinIf->decodeBin(symbol, m_cCUExternRefSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      refflag = symbol;
#endif

#else

#if MLIP_REF_LINE_NUM > 2
      m_pcTDecBinIf->decodeBin(symbol, m_cCUExternRefSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      if (symbol) {
          m_pcTDecBinIf->decodeBin(symbol, m_cCUExternRefSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
#if MLIP_REF_LINE_NUM > 3
          if (symbol) {
              m_pcTDecBinIf->decodeBin(symbol, m_cCUExternRefSCModel.get(0, 0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
              symbol++;
          }
#endif
          symbol++;
      }
      refflag = symbol;
#else
      m_pcTDecBinIf->decodeBin(symbol, m_cCUExternRefSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      refflag = symbol;
#endif

#endif 



  }
  pcCU->setExternRefSubParts(CHANNEL_TYPE_LUMA, refflag, absPartIdx, depth);

#endif



#if !JVET_C0024_QTBT
  if (mode==SIZE_NxN)
  {
    depth++;
  }
#endif
#if !MLIP_EXTEND_REF_LINE
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()>>depth]+2, CHANNEL_TYPE_LUMA);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth()>>depth]+2, CHANNEL_TYPE_LUMA);
#endif
#endif
#endif
  for (j=0;j<partNum;j++)
  {
    m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    mpmPred[j] = symbol;
  }
  for (j=0;j<partNum;j++)
  {
#if VCEG_AZ07_INTRA_65ANG_MODES
    Int preds[6] = {-1, -1, -1, -1, -1, -1};
#if !JVET_C0055_INTRA_MPM
    Int iLeftAboveCase=0;
#endif
#else
    Int preds[NUM_MOST_PROBABLE_MODES] = {-1, -1, -1};
#endif
    pcCU->getIntraDirPredictor(absPartIdx+partOffset*j, preds, COMPONENT_Y
#if VCEG_AZ07_INTRA_65ANG_MODES && !JVET_C0055_INTRA_MPM
      , iLeftAboveCase
#endif
#if MLIP_EXTEND_REF_LINE
        , refflag
#endif
      );

    if (mpmPred[j])
    {
#if JVET_C0055_INTRA_MPM
      m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[0]]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#elif VCEG_AZ07_INTRA_65ANG_MODES
      m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM0[iLeftAboveCase]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#else
      m_pcTDecBinIf->decodeBinEP( symbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif

      if (symbol)
      {
#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_C0055_INTRA_MPM
        m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[1]]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
        if( symbol )
        {
          m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[2]]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#else
        m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM1[iLeftAboveCase]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
        if( symbol )
        {
          m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM2[iLeftAboveCase]) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif

          if( symbol )
          {
            m_pcTDecBinIf->decodeBinEP( symbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
            if( symbol )
            {
              m_pcTDecBinIf->decodeBinEP( symbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
              symbol++;
            }
            symbol++;
          }
          symbol++;
        }
#else
        m_pcTDecBinIf->decodeBinEP( symbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif
        symbol++;
      }
      intraPredMode = preds[symbol];
    }
    else
    {
#if JVET_B0051_NON_MPM_MODE
#if JVET_C0024_QTBT
#if JVET_C0024_BT_FIX_TICKET22
      m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, 9) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#else
      m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, 9/3) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif
#else
      m_pcTDecBinIf->decodeBin( symbol, m_cCUIntraPredSCModel.get( 0, 0, 9+mode/3) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif

      if(!symbol) //Non-selected mode
      {
#if MLIP_EXTEND_REF_LINE
          xReadTruncBinCode(symbol, (45 - (refflag > 0)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP));
#else
        xReadTruncBinCode(symbol, 45 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );     
#endif
        symbol  += (symbol/3) ;
        symbol ++;       
      }
      else // selected mode
      {
#if MLIP_EXTEND_REF_LINE
          xReadTruncBinCode(symbol, (16 - (refflag > 0)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP));
#else
        m_pcTDecBinIf->decodeBinsEP( symbol, 4  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif
        symbol <<=2;
      }

      intraPredMode = symbol;

      //postponed sorting of MPMs (only in remaining branch)
#if VCEG_AZ07_INTRA_65ANG_MODES
      std::sort(preds, preds+6);
#else
      if (preds[0] > preds[1])
      {
        std::swap(preds[0], preds[1]);
      }
      if (preds[0] > preds[2])
      {
        std::swap(preds[0], preds[2]);
      }
      if (preds[1] > preds[2])
      {
        std::swap(preds[1], preds[2]);
      }
#endif
#if MLIP_EXTEND_REF_LINE
      if (refflag)
          intraPredMode += 2;
#endif
      for ( UInt i = 0; i < NUM_MOST_PROBABLE_MODES; i++ )
      {
        intraPredMode += ( intraPredMode >= preds[i] );
      }
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
      m_pcTDecBinIf->decodeBinsEP( symbol, 4 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
      symbol <<= 2;
      if( symbol<(NUM_INTRA_MODE-8) )
      {
        UInt symbol0;
        m_pcTDecBinIf->decodeBinsEP( symbol0, 2 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
        symbol += symbol0;
      }
#else
      m_pcTDecBinIf->decodeBinsEP( symbol, 5 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
#endif
      intraPredMode = symbol;

      //postponed sorting of MPMs (only in remaining branch)
#if VCEG_AZ07_INTRA_65ANG_MODES
      std::sort(preds, preds+6);
#else
      if (preds[0] > preds[1])
      {
        std::swap(preds[0], preds[1]);
      }
      if (preds[0] > preds[2])
      {
        std::swap(preds[0], preds[2]);
      }
      if (preds[1] > preds[2])
      {
        std::swap(preds[1], preds[2]);
      }
#endif
      for ( UInt i = 0; i < NUM_MOST_PROBABLE_MODES; i++ )
      {
        intraPredMode += ( intraPredMode >= preds[i] );
      }
#endif
    }
    pcCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, (UChar)intraPredMode, absPartIdx+partOffset*j, depth );
#if DIMD_INTRA_PRED
    pcCU->setSecIntraDirSubParts((UChar)intraPredMode, absPartIdx + partOffset*j, depth);
#endif

  }



#if RRC_DWDIP
  const Bool bUseDWDIP = pcCU->getSlice()->getSPS()->getUseDWDIP();

  UInt symbol0 = 0;
  if ( bUseDWDIP &&
      isEnabledForDWDIP(CHANNEL_TYPE_LUMA, pcCU->getWidth(absPartIdx) , pcCU->getHeight(absPartIdx), pcCU->getIntraDir(CHANNEL_TYPE_LUMA, absPartIdx))  
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC 
    && pcCU->getPDPCIdx(absPartIdx) == 0
#endif
#if !RRC_DWDIP_2Nx2N_CHECK_DISABLE
    && mode == SIZE_2Nx2N
#endif
#if MLIP_EXTEND_REF_LINE
      && refflag == 0
#endif
    )
  {
    UInt uiCtxIdx = 0;
    m_pcTDecBinIf->decodeBin(symbol0, m_cDWDIPIdxSCModel.get(0, 0, uiCtxIdx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DWDIP_FLAG));
  }
  pcCU->setDWDIPIdxSubParts(symbol0, absPartIdx, depth);
  
#endif

#if RRC_DWDIP_MLR
  {
    UInt symbol0 = 0;
    if (bUseDWDIP &&
        isEnabledForDWDIPGP(CHANNEL_TYPE_LUMA, pcCU->getWidth(absPartIdx), pcCU->getHeight(absPartIdx), pcCU->getIntraDir(CHANNEL_TYPE_LUMA, absPartIdx))
#if MLIP_EXTEND_REF_LINE
        && refflag == 0
#endif
      )
    {     
      UInt uiCtxIdx = 0;
      m_pcTDecBinIf->decodeBin(symbol0, m_cDWDIP_GPIdxSCModel.get(0, 0, uiCtxIdx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DWDIP_MLR_FLAG));
    }
    pcCU->setDWDIP_GPIdxSubParts(symbol0, absPartIdx, depth);
  }
#endif
}

#if JVET_E0062_MULTI_DMS
Void TDecSbac::parseIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiSymbol;
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth]+2, CHANNEL_TYPE_CHROMA);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth]+2, CHANNEL_TYPE_CHROMA);
#endif
#endif
#if MLIP_EXTEND_REF_LINE

  UInt refflag = 0;
  pcCU->setExternRefSubParts(CHANNEL_TYPE_CHROMA, refflag, uiAbsPartIdx, uiDepth);

#endif
#if COM16_C806_LMCHROMA
  Int iStartIdx = 1;
#else
  Int iStartIdx = 0;
#endif

#if COM16_C806_LMCHROMA
#if JVET_E0077_ENHANCED_LM
  const UInt csx = getComponentScaleX(COMPONENT_Cb, pcCU->getSlice()->getSPS()->getChromaFormatIdc());
  const UInt csy = getComponentScaleY(COMPONENT_Cb, pcCU->getSlice()->getSPS()->getChromaFormatIdc());

  Int iBlockSize = (pcCU->getHeight(uiAbsPartIdx) >> csy) + (pcCU->getWidth(uiAbsPartIdx) >> csx);
#if JVET_E0077_MMLM
    if (iBlockSize >= g_aiMMLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
    {
        iStartIdx += JVET_E0077_MMLM;
    }
#endif
#if JVET_E0077_LM_MF
    if (iBlockSize >= g_aiMFLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
    {
        iStartIdx += LM_FILTER_NUM;
    }
#endif

#endif
  if( pcCU->getSlice()->getSPS()->getUseLMChroma() )
  {
    m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUChromaPredSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
  }
  else
  {
    uiSymbol = 1;
  }
  if( uiSymbol == 0 )
  {
    uiSymbol = LM_CHROMA_IDX;
#if JVET_E0077_ENHANCED_LM
        Int iCtx = 6;
#if JVET_E0077_MMLM
        if (iBlockSize >= g_aiMMLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
        {
            m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUChromaPredSCModel.get(0, 0, iCtx++) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
            if (uiSymbol == 0)
            {
                uiSymbol = LM_CHROMA_IDX;
            }
            else
            {
                uiSymbol = MMLM_CHROMA_IDX;
            }
        }
#endif
#if JVET_E0077_LM_MF
        if (iBlockSize >= g_aiMFLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
        {
            if (uiSymbol == LM_CHROMA_IDX)
            {
                uiSymbol = 0;

                m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUChromaPredSCModel.get(0, 0, iCtx++) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
                if (uiSymbol == 1)
                {
                    uiSymbol = LM_CHROMA_IDX;
                }
                if (uiSymbol == 0)
                {
                    Int iLabel = 0;
                    m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUChromaPredSCModel.get(0, 0, iCtx++) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
                    iLabel = uiSymbol << 1;
                    m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUChromaPredSCModel.get(0, 0, iCtx++) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
                    iLabel += uiSymbol;
                    uiSymbol = LM_CHROMA_F1_IDX + iLabel;
                }
            }
        }
#endif

#endif
  } 
  else
  {
#endif
    UInt uiAllowedChromaDir[ NUM_CHROMA_MODE ];
    pcCU->getAllowedChromaDir( uiAbsPartIdx, uiAllowedChromaDir );
    {   
      UInt ictxIdx = 1;
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUChromaPredSCModel.get( 0, 0, ictxIdx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
      UInt uiMaxSymbol = NUM_DM_MODES;
      if( uiSymbol == 1 && uiMaxSymbol!= 1)
      {
        uiSymbol = 0;
        UInt uiCont;
        do
        {
          ictxIdx ++;
          m_pcTDecBinIf->decodeBin( uiCont, m_cCUChromaPredSCModel.get( 0, 0, ictxIdx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
          uiSymbol++;
        } while( uiCont && ( uiSymbol < (uiMaxSymbol-1) ) );       
      }      
      uiSymbol = uiAllowedChromaDir[iStartIdx + uiSymbol];
    }
#if COM16_C806_LMCHROMA
  }
#endif
  pcCU->setIntraDirSubParts( CHANNEL_TYPE_CHROMA, uiSymbol, uiAbsPartIdx, uiDepth );
}
#else
Void TDecSbac::parseIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiSymbol;
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth]+2, CHANNEL_TYPE_CHROMA);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth]+2, CHANNEL_TYPE_CHROMA);
#endif
#endif
#if MLIP_EXTEND_REF_LINE

  UInt refflag = 0;
  pcCU->setExternRefSubParts(CHANNEL_TYPE_CHROMA, refflag, uiAbsPartIdx, uiDepth);

#endif
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUChromaPredSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
  if( uiSymbol == 0 )
  {
    uiSymbol = DM_CHROMA_IDX;
  }
  else
  {
#if COM16_C806_LMCHROMA
#if JVET_E0077_ENHANCED_LM
      uiSymbol = -1;
      if (pcCU->getSlice()->getSPS()->getUseLMChroma())
      {
          uiSymbol = parseLMMode(pcCU, uiAbsPartIdx, uiDepth);

      }
      if (uiSymbol != -1)
      {
          //Do nothing
      }
#else
    if( pcCU->getSlice()->getSPS()->getUseLMChroma() )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUChromaPredSCModel.get( 0, 0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }
    else
    {
      uiSymbol = 1;
    }

    if( uiSymbol == 0 )
    {
      uiSymbol = LM_CHROMA_IDX;
    } 
#endif
    else
    {
#endif
      UInt uiIPredMode;
#if MLIP_EXTEND_REF_LINE
      UInt ExternRef = pcCU->getExternRef(CHANNEL_TYPE_CHROMA, uiAbsPartIdx);
#endif
      m_pcTDecBinIf->decodeBinsEP( uiIPredMode, 2 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
      UInt uiAllowedChromaDir[ NUM_CHROMA_MODE ];
      pcCU->getAllowedChromaDir( uiAbsPartIdx, uiAllowedChromaDir );
      uiSymbol = uiAllowedChromaDir[ uiIPredMode ];
#if COM16_C806_LMCHROMA
    }
#endif
  }
  pcCU->setIntraDirSubParts( CHANNEL_TYPE_CHROMA, uiSymbol, uiAbsPartIdx, uiDepth );
}
#endif

Void TDecSbac::parseInterDir( TComDataCU* pcCU, UInt& ruiInterDir, UInt uiAbsPartIdx )
{
  UInt uiSymbol;
  const UInt uiCtx = pcCU->getCtxInterDir( uiAbsPartIdx );
  ContextModel *pCtx = m_cCUInterDirSCModel.get( 0 );

  uiSymbol = 0;
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
#if JVET_C0024_QTBT
  if (1 )
#else
  if (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag() || pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8 )
#endif
#else
#if JVET_C0024_QTBT
  if (1 )
#else
  if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8 )
#endif
#endif
  {
    m_pcTDecBinIf->decodeBin( uiSymbol, *( pCtx + uiCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__INTER_DIR) );
  }

  if( uiSymbol )
  {
    uiSymbol = 2;
  }
  else
  {
    m_pcTDecBinIf->decodeBin( uiSymbol, *( pCtx + 4 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__INTER_DIR) );
    assert(uiSymbol == 0 || uiSymbol == 1);
  }

  uiSymbol++;
  ruiInterDir = uiSymbol;
  return;
}

Void TDecSbac::parseRefFrmIdx( TComDataCU* pcCU, Int& riRefFrmIdx, RefPicList eRefList )
{
  UInt uiSymbol;

  ContextModel *pCtx = m_cCURefPicSCModel.get( 0 );
  m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__REF_FRM_IDX) );

  if( uiSymbol )
  {
    UInt uiRefNum = pcCU->getSlice()->getNumRefIdx( eRefList ) - 2;
    pCtx++;
    UInt ui;
    for( ui = 0; ui < uiRefNum; ++ui )
    {
      if( ui == 0 )
      {
        m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__REF_FRM_IDX) );
      }
      else
      {
        m_pcTDecBinIf->decodeBinEP( uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__REF_FRM_IDX) );
      }
      if( uiSymbol == 0 )
      {
        break;
      }
    }
    uiSymbol = ui + 1;
  }
  riRefFrmIdx = uiSymbol;

  return;
}

Void TDecSbac::parseMvd( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList )
{
  UInt uiSymbol;
  UInt uiHorAbs;
  UInt uiVerAbs;
  UInt uiHorSign = 0;
  UInt uiVerSign = 0;
  ContextModel *pCtx = m_cCUMvdSCModel.get( 0 );

  if(pcCU->getSlice()->getMvdL1ZeroFlag() && eRefList == REF_PIC_LIST_1 && pcCU->getInterDir(uiAbsPartIdx)==3)
  {
    uiHorAbs=0;
    uiVerAbs=0;
  }
  else
  {
    m_pcTDecBinIf->decodeBin( uiHorAbs, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
    m_pcTDecBinIf->decodeBin( uiVerAbs, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );

    const Bool bHorAbsGr0 = uiHorAbs != 0;
    const Bool bVerAbsGr0 = uiVerAbs != 0;
    pCtx++;

    if( bHorAbsGr0 )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
      uiHorAbs += uiSymbol;
    }

    if( bVerAbsGr0 )
    {
      m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
      uiVerAbs += uiSymbol;
    }

    if( bHorAbsGr0 )
    {
      if( 2 == uiHorAbs )
      {
        xReadEpExGolomb( uiSymbol, 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
        uiHorAbs += uiSymbol;
      }

      m_pcTDecBinIf->decodeBinEP( uiHorSign RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
    }

    if( bVerAbsGr0 )
    {
      if( 2 == uiVerAbs )
      {
        xReadEpExGolomb( uiSymbol, 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
        uiVerAbs += uiSymbol;
      }

      m_pcTDecBinIf->decodeBinEP( uiVerSign RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
    }

  }

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  TComMv cMv( uiHorSign ? -Int( uiHorAbs ): uiHorAbs, uiVerSign ? -Int( uiVerAbs ) : uiVerAbs );
  cMv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#else
  const TComMv cMv( uiHorSign ? -Int( uiHorAbs ): uiHorAbs, uiVerSign ? -Int( uiVerAbs ) : uiVerAbs );
#endif
#if JVET_C0024_QTBT
  pcCU->getCUMvField( eRefList )->setAllMvd( cMv, SIZE_2Nx2N, uiAbsPartIdx, uiDepth, uiPartIdx );
#else
  pcCU->getCUMvField( eRefList )->setAllMvd( cMv, pcCU->getPartitionSize( uiAbsPartIdx ), uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
  return;
}

Void TDecSbac::parseCrossComponentPrediction( TComTU &rTu, ComponentID compID )
{
  TComDataCU *pcCU = rTu.getCU();

  if( isLuma(compID) || !pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag() )
  {
    return;
  }

  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#if JVET_E0062_MULTI_DMS
#if JVET_E0077_ENHANCED_LM
  if (!pcCU->isIntra(uiAbsPartIdx) || !IsLMMode(pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiAbsPartIdx)))
#else
  if (!pcCU->isIntra(uiAbsPartIdx) || (pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiAbsPartIdx) != LM_CHROMA_IDX))
#endif
#else
  if (!pcCU->isIntra(uiAbsPartIdx) || (pcCU->getIntraDir( CHANNEL_TYPE_CHROMA, uiAbsPartIdx ) == DM_CHROMA_IDX))
#endif
  {
    Char alpha  = 0;
    UInt symbol = 0;

    DTRACE_CABAC_VL( g_nSymbolCounter++ )
    DTRACE_CABAC_T("\tparseCrossComponentPrediction()")
    DTRACE_CABAC_T( "\tAddr=" )
    DTRACE_CABAC_V( compID )
    DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
    DTRACE_CABAC_V( uiAbsPartIdx )
#if RExt__DECODER_DEBUG_BIT_STATISTICS
    TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__CROSS_COMPONENT_PREDICTION, (g_aucConvertToBit[rTu.getRect(compID).width] + 2), compID);
#endif
    ContextModel *pCtx = m_cCrossComponentPredictionSCModel.get(0, 0) + ((compID == COMPONENT_Cr) ? (NUM_CROSS_COMPONENT_PREDICTION_CTX >> 1) : 0);
    m_pcTDecBinIf->decodeBin( symbol, pCtx[0] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

    if(symbol != 0)
    {
      // Cross-component prediction alpha is non-zero.
      UInt sign = 0;
      m_pcTDecBinIf->decodeBin( symbol, pCtx[1] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

      if (symbol != 0)
      {
        // alpha is 2 (symbol=1), 4(symbol=2) or 8(symbol=3).
        // Read up to two more bits
        xReadUnaryMaxSymbol( symbol, (pCtx + 2), 1, 2 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
        symbol += 1;
      }
      m_pcTDecBinIf->decodeBin( sign, pCtx[4] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

      alpha = (sign != 0) ? -(1 << symbol) : (1 << symbol);
    }
    DTRACE_CABAC_T( "\tAlpha=" )
    DTRACE_CABAC_V( alpha )
    DTRACE_CABAC_T( "\n" )

    pcCU->setCrossComponentPredictionAlphaPartRange( alpha, compID, uiAbsPartIdx, rTu.GetAbsPartIdxNumParts( compID ) );
  }
}

Void TDecSbac::parseTransformSubdivFlag( UInt& ruiSubdivFlag, UInt uiLog2TransformBlockSize )
{
  m_pcTDecBinIf->decodeBin( ruiSubdivFlag, m_cCUTransSubdivFlagSCModel.get( 0, 0, uiLog2TransformBlockSize )
      RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(TComCodingStatisticsClassType(STATS__CABAC_BITS__TRANSFORM_SUBDIV_FLAG, 5-uiLog2TransformBlockSize))
                          );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseTransformSubdivFlag()" )
  DTRACE_CABAC_T( "\tsymbol=" )
  DTRACE_CABAC_V( ruiSubdivFlag )
  DTRACE_CABAC_T( "\tctx=" )
  DTRACE_CABAC_V( uiLog2TransformBlockSize )
  DTRACE_CABAC_T( "\n" )
}

Void TDecSbac::parseQtRootCbf( UInt uiAbsPartIdx, UInt& uiQtRootCbf )
{
  UInt uiSymbol;
  const UInt uiCtx = 0;
  m_pcTDecBinIf->decodeBin( uiSymbol , m_cCUQtRootCbfSCModel.get( 0, 0, uiCtx ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__QT_ROOT_CBF) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseQtRootCbf()" )
  DTRACE_CABAC_T( "\tsymbol=" )
  DTRACE_CABAC_V( uiSymbol )
  DTRACE_CABAC_T( "\tctx=" )
  DTRACE_CABAC_V( uiCtx )
  DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
  DTRACE_CABAC_V( uiAbsPartIdx )
  DTRACE_CABAC_T( "\n" )

  uiQtRootCbf = uiSymbol;
}

Void TDecSbac::parseDeltaQP( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  Int qp;
  UInt uiDQp;
  Int  iDQp;

  UInt uiSymbol;

  xReadUnaryMaxSymbol (uiDQp,  &m_cCUDeltaQpSCModel.get( 0, 0, 0 ), 1, CU_DQP_TU_CMAX RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DELTA_QP_EP) );

  if( uiDQp >= CU_DQP_TU_CMAX)
  {
    xReadEpExGolomb( uiSymbol, CU_DQP_EG_k RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DELTA_QP_EP));
    uiDQp+=uiSymbol;
  }

  if ( uiDQp > 0 )
  {
    UInt uiSign;
    Int qpBdOffsetY = pcCU->getSlice()->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA);
    m_pcTDecBinIf->decodeBinEP(uiSign RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__DELTA_QP_EP));
    iDQp = uiDQp;
    if(uiSign)
    {
      iDQp = -iDQp;
    }
    qp = (((Int) pcCU->getRefQP( uiAbsPartIdx ) + iDQp + 52 + 2*qpBdOffsetY )%(52+qpBdOffsetY)) - qpBdOffsetY;
  }
  else
  {
    qp = pcCU->getRefQP(uiAbsPartIdx);
  }

#if JVET_C0024_DELTA_QP_FIX
  pcCU->setQPSubParts(qp, uiAbsPartIdx, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx));
#else
  pcCU->setQPSubParts(qp, uiAbsPartIdx, uiDepth);
#endif
  pcCU->setCodedQP(qp);
}

/** parse chroma qp adjustment, converting to the internal table representation.
 * \returns Void
 */
Void TDecSbac::parseChromaQpAdjustment( TComDataCU* cu, UInt absPartIdx, UInt depth )
{
  UInt symbol;
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__CHROMA_QP_ADJUSTMENT, g_aucConvertToBit[cu->getSlice()->getSPS()->getCTUSize()>>depth]+2, CHANNEL_TYPE_CHROMA);
#else
  const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__CHROMA_QP_ADJUSTMENT, g_aucConvertToBit[cu->getSlice()->getSPS()->getMaxCUWidth()>>depth]+2, CHANNEL_TYPE_CHROMA);
#endif
#endif

  Int chromaQpOffsetListLen = cu->getSlice()->getPPS()->getPpsRangeExtension().getChromaQpOffsetListLen();

  // cu_chroma_qp_offset_flag
  m_pcTDecBinIf->decodeBin( symbol, m_ChromaQpAdjFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

  if (symbol && chromaQpOffsetListLen > 1)
  {
    // cu_chroma_qp_offset_idx
    xReadUnaryMaxSymbol( symbol,  &m_ChromaQpAdjIdcSCModel.get( 0, 0, 0 ), 0, chromaQpOffsetListLen - 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    symbol++;
  }
  /* NB, symbol = 0 if outer flag is not set,
   *              1 if outer flag is set and there is no inner flag
   *              1+ otherwise */
  cu->setChromaQpAdjSubParts( symbol, absPartIdx, depth );
  cu->setCodedChromaQpAdj(symbol);
}

Void TDecSbac::parseQtCbf( TComTU &rTu, const ComponentID compID, const Bool lowestLevel )
{
  TComDataCU* pcCU = rTu.getCU();

  const UInt absPartIdx       = rTu.GetAbsPartIdxTU(compID);
#if !JVET_C0024_QTBT
  const UInt TUDepth          = rTu.GetTransformDepthRel();
#endif
  const UInt uiCtx            = pcCU->getCtxQtCbf( rTu, toChannelType(compID) );
  const UInt contextSet       = toChannelType(compID);

#if !JVET_C0024_QTBT
  const UInt width            = rTu.getRect(compID).width;
  const UInt height           = rTu.getRect(compID).height;
  const Bool canQuadSplit     = (width >= (MIN_TU_SIZE * 2)) && (height >= (MIN_TU_SIZE * 2));
  const UInt coveredPartIdxes = rTu.GetAbsPartIdxNumParts(compID);

  //             Since the CBF for chroma is coded at the highest level possible, if sub-TUs are
  //             to be coded for a 4x8 chroma TU, their CBFs must be coded at the highest 4x8 level
  //             (i.e. where luma TUs are 8x8 rather than 4x4)
  //    ___ ___
  //   |   |   | <- 4 x (8x8 luma + 4x8 4:2:2 chroma)
  //   |___|___|    each quadrant has its own chroma CBF
  //   |   |   | _ _ _ _
  //   |___|___|        |
  //   <--16--->        V
  //                   _ _
  //                  |_|_| <- 4 x 4x4 luma + 1 x 4x8 4:2:2 chroma
  //                  |_|_|    no chroma CBF is coded - instead the parent CBF is inherited
  //                  <-8->    if sub-TUs are present, their CBFs had to be coded at the parent level

  const UInt lowestTUDepth = TUDepth + ((!lowestLevel && !canQuadSplit) ? 1 : 0); //unsplittable TUs inherit their parent's CBF
        UInt lowestTUCBF   = 0;

  if ((width != height) && (lowestLevel || !canQuadSplit)) //if sub-TUs are present
  {
    const UInt subTUDepth        = lowestTUDepth + 1;
    const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

    UInt combinedSubTUCBF = 0;

    for (UInt subTU = 0; subTU < 2; subTU++)
    {
      UInt uiCbf = MAX_UINT;
      m_pcTDecBinIf->decodeBin(uiCbf, m_cCUQtCbfSCModel.get(0, contextSet, uiCtx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(TComCodingStatisticsClassType(STATS__CABAC_BITS__QT_CBF, g_aucConvertToBit[rTu.getRect(compID).width]+2, compID)));

      const UInt subTUAbsPartIdx = absPartIdx + (subTU * partIdxesPerSubTU);
      pcCU->setCbfPartRange((uiCbf << subTUDepth), compID, subTUAbsPartIdx, partIdxesPerSubTU);
      combinedSubTUCBF |= uiCbf;

      DTRACE_CABAC_VL( g_nSymbolCounter++ )
      DTRACE_CABAC_T( "\tparseQtCbf()" )
      DTRACE_CABAC_T( "\tsub-TU=" )
      DTRACE_CABAC_V( subTU )
      DTRACE_CABAC_T( "\tsymbol=" )
      DTRACE_CABAC_V( uiCbf )
      DTRACE_CABAC_T( "\tctx=" )
      DTRACE_CABAC_V( uiCtx )
      DTRACE_CABAC_T( "\tetype=" )
      DTRACE_CABAC_V( compID )
      DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
      DTRACE_CABAC_V( subTUAbsPartIdx )
      DTRACE_CABAC_T( "\n" )
    }

    //propagate the sub-TU CBF up to the lowest TU level
    if (combinedSubTUCBF != 0)
    {
      pcCU->bitwiseOrCbfPartRange((combinedSubTUCBF << lowestTUDepth), compID, absPartIdx, coveredPartIdxes);
      lowestTUCBF = combinedSubTUCBF;
    }
  }
  else
  {
#endif
    UInt uiCbf = MAX_UINT;
    m_pcTDecBinIf->decodeBin(uiCbf, m_cCUQtCbfSCModel.get(0, contextSet, uiCtx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(TComCodingStatisticsClassType(STATS__CABAC_BITS__QT_CBF, g_aucConvertToBit[rTu.getRect(compID).width]+2, compID)));

#if JVET_C0024_QTBT
    pcCU->setCbfSubParts(uiCbf, compID, absPartIdx, 0);
#else
    pcCU->setCbfSubParts((uiCbf << lowestTUDepth), compID, absPartIdx, rTu.GetTransformDepthTotalAdj(compID));
#endif

    DTRACE_CABAC_VL( g_nSymbolCounter++ )
    DTRACE_CABAC_T( "\tparseQtCbf()" )
    DTRACE_CABAC_T( "\tsymbol=" )
    DTRACE_CABAC_V( uiCbf )
    DTRACE_CABAC_T( "\tctx=" )
    DTRACE_CABAC_V( uiCtx )
    DTRACE_CABAC_T( "\tetype=" )
    DTRACE_CABAC_V( compID )
    DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
    DTRACE_CABAC_V( rTu.GetAbsPartIdxTU(compID) )
    DTRACE_CABAC_T( "\n" )

#if !JVET_C0024_QTBT
    lowestTUCBF = uiCbf;
  }

  //propagate the lowest level CBF up to the current level
  if (lowestTUCBF != 0)
  {
    for (UInt depth = TUDepth; depth < lowestTUDepth; depth++)
    {
      pcCU->bitwiseOrCbfPartRange((lowestTUCBF << depth), compID, absPartIdx, coveredPartIdxes);
    }
  }
#endif
}


Void TDecSbac::parseTransformSkipFlags (TComTU &rTu, ComponentID component)
{
  TComDataCU* pcCU=rTu.getCU();
  UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU(component);

  if (pcCU->getCUTransquantBypass(uiAbsPartIdx))
  {
    return;
  }

#if JVET_C0024_QTBT
  if (!TUCompRectHasAssociatedTransformSkipFlag(pcCU->getSlice()->isIntra(), rTu.getRect(component), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize()))
#else
  if (!TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(component), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize()))
#endif
  {
    return;
  }
#if JVET_F0031_RMV_REDUNDANT_TRSKIP && COM16_C806_EMT
  if (isLuma(component) && pcCU->getEmtCuFlag(uiAbsPartIdx))
  {
    return;
  }
#endif

  UInt useTransformSkip;

  m_pcTDecBinIf->decodeBin( useTransformSkip , m_cTransformSkipSCModel.get( 0, toChannelType(component), 0 )
      RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(TComCodingStatisticsClassType(STATS__CABAC_BITS__TRANSFORM_SKIP_FLAGS, component))
                          );

  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T("\tparseTransformSkip()");
  DTRACE_CABAC_T( "\tsymbol=" )
  DTRACE_CABAC_V( useTransformSkip )
  DTRACE_CABAC_T( "\tAddr=" )
  DTRACE_CABAC_V( pcCU->getCtuRsAddr() )
  DTRACE_CABAC_T( "\tetype=" )
  DTRACE_CABAC_V( component )
  DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
  DTRACE_CABAC_V( rTu.GetAbsPartIdxTU() )
  DTRACE_CABAC_T( "\n" )

#if JVET_C0024_QTBT
  pcCU->setTransformSkipPartRange( useTransformSkip, component, uiAbsPartIdx, 0);
#else
  pcCU->setTransformSkipPartRange( useTransformSkip, component, uiAbsPartIdx, rTu.GetAbsPartIdxNumParts(component));
#endif
}

#if VCEG_AZ08_KLT_COMMON
Void TDecSbac::parseKLTFlags(TComTU &rTu, ComponentID component)
{
    TComDataCU* pcCU = rTu.getCU();
    UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU(component);

    if (pcCU->getCUTransquantBypass(uiAbsPartIdx))
    {
        return;
    }
    UInt useKLTFlag = 0;

    m_pcTDecBinIf->decodeBin(useKLTFlag, m_cKLTFlagSCModel.get(0, toChannelType(component), 0)
        RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(TComCodingStatisticsClassType(STATS__CABAC_BITS__TRANSFORM_SKIP_FLAGS, component))
        );

    DTRACE_CABAC_VL(g_nSymbolCounter++)
    DTRACE_CABAC_T("\tparseKLTFlag()");
    DTRACE_CABAC_T("\tsymbol=")
    DTRACE_CABAC_V(useKLTFlag)
    DTRACE_CABAC_T("\tAddr=")
    DTRACE_CABAC_V(pcCU->getCtuRsAddr())
    DTRACE_CABAC_T("\tetype=")
    DTRACE_CABAC_V(component)
    DTRACE_CABAC_T("\tuiAbsPartIdx=")
    DTRACE_CABAC_V(rTu.GetAbsPartIdxTU())
    DTRACE_CABAC_T("\n")

    pcCU->setKLTPartRange(useKLTFlag, component, uiAbsPartIdx, rTu.GetAbsPartIdxNumParts(component));
}
#endif

/** Parse (X,Y) position of the last significant coefficient
 * \param uiPosLastX reference to X component of last coefficient
 * \param uiPosLastY reference to Y component of last coefficient
 * \param width  Block width
 * \param height Block height
 * \param component chroma compinent ID
 * \param uiScanIdx scan type (zig-zag, hor, ver)
 *
 * This method decodes the X and Y component within a block of the last significant coefficient.
 */
Void TDecSbac::parseLastSignificantXY( UInt& uiPosLastX, UInt& uiPosLastY, Int width, Int height, ComponentID component, UInt uiScanIdx )
{
  UInt uiLast;

  ContextModel *pCtxX = m_cCuCtxLastX.get( 0, toChannelType(component) );
  ContextModel *pCtxY = m_cCuCtxLastY.get( 0, toChannelType(component) );

#if RExt__DECODER_DEBUG_BIT_STATISTICS
  TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__LAST_SIG_X_Y, g_aucConvertToBit[width]+2, component);
#endif


  if ( uiScanIdx == SCAN_VER )
  {
    swap( width, height );
  }

#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
// posX
  Int widthCtx = component ? 4 : width;
  const UInt *puiCtxIdxX = g_uiLastCtx + ( g_aucConvertToBit[ widthCtx ] * ( g_aucConvertToBit[ widthCtx ] + 3 ) );

  for( uiPosLastX = 0; uiPosLastX < g_uiGroupIdx[ width - 1 ]; uiPosLastX++ )
  {
    if ( component  )
    {
      m_pcTDecBinIf->decodeBin( uiLast, *( pCtxX + (uiPosLastX>>g_aucConvertToBit[ width ])  ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }
    else
    {
      m_pcTDecBinIf->decodeBin( uiLast, *( pCtxX + puiCtxIdxX[ uiPosLastX ] ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }

    if( !uiLast )
    {
      break;
    }
  }

  // posY

  Int heightCtx = component ? 4 : height;
  const UInt *puiCtxIdxY = g_uiLastCtx + ( g_aucConvertToBit[ heightCtx ] * ( g_aucConvertToBit[ heightCtx ] + 3 ) );

  for( uiPosLastY = 0; uiPosLastY < g_uiGroupIdx[ height - 1 ]; uiPosLastY++ )
  {
    if (component)
    {
      m_pcTDecBinIf->decodeBin( uiLast, *( pCtxY + (uiPosLastY>>g_aucConvertToBit[ height ]) ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }
    else
    {
      m_pcTDecBinIf->decodeBin( uiLast, *( pCtxY + puiCtxIdxY[ uiPosLastY ] ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
    }

    if( !uiLast )
    {
      break;
    }
  }
#else
  Int blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY;
  getLastSignificantContextParameters(component, width, height, blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY);

  //------------------

  // posX

  for( uiPosLastX = 0; uiPosLastX < g_uiGroupIdx[ width - 1 ]; uiPosLastX++ )
  {
    m_pcTDecBinIf->decodeBin( uiLast, *( pCtxX + blkSizeOffsetX + (uiPosLastX >>shiftX) ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

    if( !uiLast )
    {
      break;
    }
  }

  // posY

  for( uiPosLastY = 0; uiPosLastY < g_uiGroupIdx[ height - 1 ]; uiPosLastY++ )
  {
    m_pcTDecBinIf->decodeBin( uiLast, *( pCtxY + blkSizeOffsetY + (uiPosLastY >>shiftY)) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );

    if( !uiLast )
    {
      break;
    }
  }
#endif
  // EP-coded part

  if ( uiPosLastX > 3 )
  {
    UInt uiTemp  = 0;
    UInt uiCount = ( uiPosLastX - 2 ) >> 1;
    for ( Int i = uiCount - 1; i >= 0; i-- )
    {
      m_pcTDecBinIf->decodeBinEP( uiLast RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
      uiTemp += uiLast << i;
    }
    uiPosLastX = g_uiMinInGroup[ uiPosLastX ] + uiTemp;
  }
  if ( uiPosLastY > 3 )
  {
    UInt uiTemp  = 0;
    UInt uiCount = ( uiPosLastY - 2 ) >> 1;
    for ( Int i = uiCount - 1; i >= 0; i-- )
    {
      m_pcTDecBinIf->decodeBinEP( uiLast RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype) );
      uiTemp += uiLast << i;
    }
    uiPosLastY = g_uiMinInGroup[ uiPosLastY ] + uiTemp;
  }

  if( uiScanIdx == SCAN_VER )
  {
    swap( uiPosLastX, uiPosLastY );
  }
}

Void TDecSbac::parseCoeffNxN(  TComTU &rTu, ComponentID compID 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    , Bool& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , Int& iNonZeroCoeffNonTs
#endif
    )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU(compID);
  const TComRectangle &rRect=rTu.getRect(compID);
  const UInt uiWidth=rRect.width;
  const UInt uiHeight=rRect.height;
  TCoeff* pcCoef=(pcCU->getCoeff(compID)+rTu.getCoefficientOffset(compID));
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());
#if COM16_C806_EMT
  UInt uiNumSig=0;
#endif

  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseCoeffNxN()\teType=" )
  DTRACE_CABAC_V( compID )
  DTRACE_CABAC_T( "\twidth=" )
  DTRACE_CABAC_V( uiWidth )
  DTRACE_CABAC_T( "\theight=" )
  DTRACE_CABAC_V( uiHeight )
  DTRACE_CABAC_T( "\tdepth=" )
//  DTRACE_CABAC_V( rTu.GetTransformDepthTotalAdj(compID) )
  DTRACE_CABAC_V( rTu.GetTransformDepthTotal() )
  DTRACE_CABAC_T( "\tabspartidx=" )
//  DTRACE_CABAC_V( uiAbsPartIdx )
  DTRACE_CABAC_V( rTu.GetAbsPartIdxTU(compID) )
  DTRACE_CABAC_T( "\ttoCU-X=" )
  DTRACE_CABAC_V( pcCU->getCUPelX() )
  DTRACE_CABAC_T( "\ttoCU-Y=" )
  DTRACE_CABAC_V( pcCU->getCUPelY() )
  DTRACE_CABAC_T( "\tCU-addr=" )
  DTRACE_CABAC_V(  pcCU->getCtuRsAddr() )
  DTRACE_CABAC_T( "\tinCU-X=" )
//  DTRACE_CABAC_V( g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ] )
  DTRACE_CABAC_V( g_auiRasterToPelX[ g_auiZscanToRaster[rTu.GetAbsPartIdxTU(compID)] ] )
  DTRACE_CABAC_T( "\tinCU-Y=" )
// DTRACE_CABAC_V( g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ] )
  DTRACE_CABAC_V( g_auiRasterToPelY[ g_auiZscanToRaster[rTu.GetAbsPartIdxTU(compID)] ] )
  DTRACE_CABAC_T( "\tpredmode=" )
  DTRACE_CABAC_V(  pcCU->getPredictionMode( uiAbsPartIdx ) )
  DTRACE_CABAC_T( "\n" )

  //--------------------------------------------------------------------------------------------------

#if JVET_C0024_QTBT
  if( uiWidth > sps.getCTUSize() )
#else
  if( uiWidth > sps.getMaxTrSize() )
#endif
  {
    std::cerr << "ERROR: parseCoeffNxN was passed a TU with dimensions larger than the maximum allowed size" << std::endl;
    assert(false);
    exit(1);
  }

  //--------------------------------------------------------------------------------------------------

  //set parameters

  const ChannelType  chType            = toChannelType(compID);
#if JVET_C0024_QTBT
  const UInt         uiLog2BlockWidth  = g_aucConvertToBit[ uiWidth  ] + MIN_CU_LOG2;
  const UInt         uiLog2BlockHeight  = g_aucConvertToBit[ uiHeight  ] + MIN_CU_LOG2;
#if VCEG_AZ07_CTX_RESIDUALCODING
  const UInt         uiLog2BlockSize = (uiLog2BlockWidth+uiLog2BlockHeight)>>1;
#endif
#else
  const UInt         uiLog2BlockWidth  = g_aucConvertToBit[ uiWidth  ] + 2;
//#if !VCEG_AZ07_CTX_RESIDUALCODING || JVET_C0046_ZO_ASSERT ticket#24
#if !VCEG_AZ07_CTX_RESIDUALCODING
  const UInt         uiLog2BlockHeight = g_aucConvertToBit[ uiHeight ] + 2;
#endif
#endif
  const UInt         uiMaxNumCoeff     = uiWidth * uiHeight;
  const UInt         uiMaxNumCoeffM1   = uiMaxNumCoeff - 1;
  const Int          log2CoeffGroupSize = 2 - ((uiWidth & 0x03) != 0 || (uiHeight & 0x03) != 0);


  const ChannelType  channelType       = toChannelType(compID);
  const Bool         extendedPrecision = sps.getSpsRangeExtension().getExtendedPrecisionProcessingFlag();

  const Bool         alignCABACBeforeBypass = sps.getSpsRangeExtension().getCabacBypassAlignmentEnabledFlag();
  const Int          maxLog2TrDynamicRange  = sps.getMaxLog2TrDynamicRange(channelType);

#if RExt__DECODER_DEBUG_BIT_STATISTICS
  TComCodingStatisticsClassType ctype_group(STATS__CABAC_BITS__SIG_COEFF_GROUP_FLAG, uiLog2BlockWidth, compID);
  TComCodingStatisticsClassType ctype_map(STATS__CABAC_BITS__SIG_COEFF_MAP_FLAG, uiLog2BlockWidth, compID);
  TComCodingStatisticsClassType ctype_gt1(STATS__CABAC_BITS__GT1_FLAG, uiLog2BlockWidth, compID);
  TComCodingStatisticsClassType ctype_gt2(STATS__CABAC_BITS__GT2_FLAG, uiLog2BlockWidth, compID);
#endif

  Bool beValid;
  if (pcCU->getCUTransquantBypass(uiAbsPartIdx))
  {
    beValid = false;
    if((!pcCU->isIntra(uiAbsPartIdx)) && pcCU->isRDPCMEnabled(uiAbsPartIdx))
    {
      parseExplicitRdpcmMode(rTu, compID);
    }
  }
  else
  {
    beValid = pcCU->getSlice()->getPPS()->getSignHideFlag();
  }

  UInt absSum = 0;

  //--------------------------------------------------------------------------------------------------

  if(pcCU->getSlice()->getPPS()->getUseTransformSkip())
  {
    parseTransformSkipFlags(rTu, compID);
    //  This TU has coefficients and is transform skipped. Check whether is inter coded and if yes decode the explicit RDPCM mode
    if(pcCU->getTransformSkip(uiAbsPartIdx, compID) && (!pcCU->isIntra(uiAbsPartIdx)) && pcCU->isRDPCMEnabled(uiAbsPartIdx) )
    {
      parseExplicitRdpcmMode(rTu, compID);
      if(pcCU->getExplicitRdpcmMode(compID, uiAbsPartIdx) != RDPCM_OFF)
      {
        //  Sign data hiding is avoided for horizontal and vertical RDPCM modes
        beValid = false;
      }
    }
  }

#if VCEG_AZ08_KLT_COMMON
#if VCEG_AZ08_USE_KLT
  if (pcCU->getSlice()->getSPS()->getUseKLT())
  {
#endif
      UInt uiMaxTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
      UInt uiMinTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MIN - 1];
      Bool bCheckKLTFlag = (toChannelType(compID) == CHANNEL_TYPE_LUMA) && (uiWidth == uiHeight) && (uiWidth <= uiMaxTrWidth) && (uiWidth >= uiMinTrWidth);
      if (bCheckKLTFlag && pcCU->getSlice()->getPPS()->getUseTransformSkip())
      {
          UInt useTransformSkip = pcCU->getTransformSkip(uiAbsPartIdx, compID);
          bCheckKLTFlag &= !useTransformSkip;
      }

#if VCEG_AZ08_USE_KLT
      if (pcCU->getSlice()->getSPS()->getUseInterKLT() && !pcCU->getSlice()->getSPS()->getUseIntraKLT()) //only inter
      {
          bCheckKLTFlag &= (!pcCU->isIntra(uiAbsPartIdx));
      }
      else if(!pcCU->getSlice()->getSPS()->getUseInterKLT() && pcCU->getSlice()->getSPS()->getUseIntraKLT()) //only intra
      {
          bCheckKLTFlag &= (pcCU->isIntra(uiAbsPartIdx));
      }
      else if ((!pcCU->getSlice()->getSPS()->getUseInterKLT()) && (!pcCU->getSlice()->getSPS()->getUseIntraKLT())) //neither
      {
          bCheckKLTFlag = false;
      }
#else
#if VCEG_AZ08_INTER_KLT && !VCEG_AZ08_INTRA_KLT //only inter
      bCheckKLTFlag &= (!pcCU->isIntra(uiAbsPartIdx));
#endif
#if !VCEG_AZ08_INTER_KLT && VCEG_AZ08_INTRA_KLT //only intra
      bCheckKLTFlag &= (pcCU->isIntra(uiAbsPartIdx));
#endif
#if !VCEG_AZ08_INTER_KLT && !VCEG_AZ08_INTRA_KLT //none
      bCheckKLTFlag = false;
#endif
#endif

      if (bCheckKLTFlag)
      {
          parseKLTFlags(rTu, compID);
      }
#if VCEG_AZ08_USE_KLT
  }
#endif
#endif

  Int uiIntraMode = -1;
  const Bool       bIsLuma = isLuma(compID);
  Int isIntra = pcCU->isIntra(uiAbsPartIdx) ? 1 : 0;
  if ( isIntra && pcCU->isRDPCMEnabled(uiAbsPartIdx) )
  {
#if !JVET_C0024_QTBT
    const UInt partsPerMinCU = 1<<(2*(sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
#endif
    uiIntraMode = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx );
#if JVET_C0024_QTBT
#if !JVET_E0062_MULTI_DMS
    uiIntraMode = (uiIntraMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) : uiIntraMode;
#endif
#else
    uiIntraMode = (uiIntraMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, rTu.GetChromaFormat(), partsPerMinCU)) : uiIntraMode;
#endif
    uiIntraMode = ((rTu.GetChromaFormat() == CHROMA_422) && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiIntraMode] : uiIntraMode;
    Bool transformSkip = pcCU->getTransformSkip( uiAbsPartIdx,compID);
    Bool rdpcm_lossy = ( transformSkip /*&& isIntra*/ && ( (uiIntraMode == HOR_IDX) || (uiIntraMode == VER_IDX) ) );
    if ( rdpcm_lossy )
    {
      beValid = false;
    }
  }

  //--------------------------------------------------------------------------------------------------
  const Bool  bUseGolombRiceParameterAdaptation = sps.getSpsRangeExtension().getPersistentRiceAdaptationEnabledFlag();
        UInt &currentGolombRiceStatistic        = m_golombRiceAdaptationStatistics[rTu.getGolombRiceStatisticsIndex(compID)];

  //select scans
  TUEntropyCodingParameters codingParameters;
  getTUEntropyCodingParameters(codingParameters, rTu, compID);

  //===== decode last significant =====
  UInt uiPosLastX, uiPosLastY;
  parseLastSignificantXY( uiPosLastX, uiPosLastY, uiWidth, uiHeight, compID, codingParameters.scanType );
  UInt uiBlkPosLast      = uiPosLastX + (uiPosLastY<<uiLog2BlockWidth);
  pcCoef[ uiBlkPosLast ] = 1;

#if JVET_C0046_ZO_ASSERT && JVET_C0046_ZO_ASSERT_LAST_COEF
  if ( ((uiWidth > ZERO_OUT_TH) || (uiHeight > ZERO_OUT_TH)) &&
       (!pcCU->getTransformSkip(compID) && !pcCU->getCUTransquantBypass(uiAbsPartIdx)))
  {
     // last coeff shall be in the low freqecy domain
     assert((uiPosLastX < ZERO_OUT_TH) && (uiPosLastY < ZERO_OUT_TH));
  }
#endif


  //===== decode significance flags =====
  UInt uiScanPosLast;
  for( uiScanPosLast = 0; uiScanPosLast < uiMaxNumCoeffM1; uiScanPosLast++ )
  {
    UInt uiBlkPos = codingParameters.scan[ uiScanPosLast ];
    if( uiBlkPosLast == uiBlkPos )
    {
      break;
    }
  }
#if JVET_C0024_QTBT
  const UInt uiCGSizeLog2 = (uiWidth==2 || uiHeight==2) ? MLS_CG_SIZE-2: MLS_CG_SIZE;
#endif

  ContextModel * const baseCoeffGroupCtx = m_cCUSigCoeffGroupSCModel.get( 0, isChroma(chType) );
#if VCEG_AZ07_CTX_RESIDUALCODING
  UInt uiPrevGRParam = 0;
#if JVET_C0024_QTBT
  UInt iOffsetonTU = uiLog2BlockSize<=2 ? 0: NUM_SIG_FLAG_CTX_LUMA_TU << ( min(1, (Int)(uiLog2BlockSize - 3)) );
#else
  UInt iOffsetonTU = (g_aucConvertToBit[ uiWidth ] > 2 ? 2:g_aucConvertToBit[ uiWidth ])*NUM_SIG_FLAG_CTX_LUMA_TU;
#endif
  ContextModel * const baseCtx = (chType == CHANNEL_TYPE_LUMA) ? m_cCUSigSCModel.get( 0, 0 ) + iOffsetonTU : m_cCUSigSCModel.get( 0, 0 ) + NUM_SIG_FLAG_CTX_LUMA;
  ContextModel * const greXCtx = (chType == CHANNEL_TYPE_LUMA) ? m_cCUOneSCModel.get( 0, 0 ) : m_cCUOneSCModel.get( 0, 0 ) + NUM_ONE_FLAG_CTX_LUMA;
#else
  ContextModel * const baseCtx = m_cCUSigSCModel.get( 0, 0 ) + getSignificanceMapContextOffset(compID);
#endif

#if JVET_C0024_QTBT
  const Int  iLastScanSet  = uiScanPosLast >> uiCGSizeLog2;
#else
  const Int  iLastScanSet  = uiScanPosLast >> MLS_CG_SIZE;
#endif
#if !VCEG_AZ07_CTX_RESIDUALCODING
  UInt c1                  = 1;
#endif
  UInt uiGoRiceParam       = 0;


  UInt uiSigCoeffGroupFlag[ MLS_GRP_NUM ];
  memset( uiSigCoeffGroupFlag, 0, sizeof(UInt) * MLS_GRP_NUM );

  Int  iScanPosSig             = (Int) uiScanPosLast;
#if VCEG_AZ07_CTX_RESIDUALCODING
  Bool bHor8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_HOR;
  Bool bVer8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_VER;
  Bool bNonZig8x8 = bHor8x8 || bVer8x8; 
#endif
  for( Int iSubSet = iLastScanSet; iSubSet >= 0; iSubSet-- )
  {
#if JVET_C0024_QTBT
    Int  iSubPos   = iSubSet << uiCGSizeLog2;
#else
    Int  iSubPos   = iSubSet << MLS_CG_SIZE;
#endif
#if !VCEG_AZ07_CTX_RESIDUALCODING
    uiGoRiceParam  = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;
#endif
    Bool updateGolombRiceStatistics = bUseGolombRiceParameterAdaptation; //leave the statistics at 0 when not using the adaptation system
    Int numNonZero = 0;
#if VCEG_AZ07_CTX_RESIDUALCODING
    UInt ctxG1 = 0;
    UInt ctxG2 = 0;
#endif

    Int lastNZPosInCG  = -1;
#if JVET_C0024_QTBT
    Int firstNZPosInCG = 1 << uiCGSizeLog2;
#else
    Int firstNZPosInCG = 1 << MLS_CG_SIZE;
#endif

    Bool escapeDataPresentInGroup = false;

    Int pos[1 << MLS_CG_SIZE];

    if( iScanPosSig == (Int) uiScanPosLast )
    {
      lastNZPosInCG  = iScanPosSig;
      firstNZPosInCG = iScanPosSig;
      iScanPosSig--;
      pos[ numNonZero ] = uiBlkPosLast;
      numNonZero = 1;
    }

    // decode significant_coeffgroup_flag
    Int iCGBlkPos = codingParameters.scanCG[ iSubSet ];
    Int iCGPosY   = iCGBlkPos / codingParameters.widthInGroups;
    Int iCGPosX   = iCGBlkPos - (iCGPosY * codingParameters.widthInGroups);
#if VCEG_AZ07_CTX_RESIDUALCODING
    if(bNonZig8x8)
    {
      iCGPosY = (bHor8x8 ? iCGBlkPos : 0);
      iCGPosX = (bVer8x8 ? iCGBlkPos : 0);
    }
#endif

    if( iSubSet == iLastScanSet || iSubSet == 0)
    {
      uiSigCoeffGroupFlag[ iCGBlkPos ] = 1;
    }
#if COM16_C806_T64 && !JVET_C0024_QTBT && !JVET_C0046_ZO_ASSERT
    else if (iCGPosY >= (ZERO_OUT_TH >> log2CoeffGroupSize) || iCGPosX >= (ZERO_OUT_TH >> log2CoeffGroupSize))
    {
      uiSigCoeffGroupFlag[ iCGBlkPos ] = 0;
    }
#endif
    else
    {
      UInt uiSigCoeffGroup;
#if VCEG_AZ07_CTX_RESIDUALCODING
      UInt uiCtxSig  = TComTrQuant::getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups, codingParameters.scanType);
#else
      UInt uiCtxSig  = TComTrQuant::getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups );
#endif
      m_pcTDecBinIf->decodeBin( uiSigCoeffGroup, baseCoeffGroupCtx[ uiCtxSig ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_group) );
      uiSigCoeffGroupFlag[ iCGBlkPos ] = uiSigCoeffGroup;
    }

#if JVET_C0046_ZO_ASSERT && JVET_C0046_ZO_ASSERT_CODED_SBK_FLAG
   if ( ((uiWidth > ZERO_OUT_TH) || (uiHeight > ZERO_OUT_TH)) 
     && (!pcCU->getTransformSkip(compID) && !pcCU->getCUTransquantBypass(uiAbsPartIdx)) )
   {
     if (iCGPosY >= (ZERO_OUT_TH >> log2CoeffGroupSize) || iCGPosX >= (ZERO_OUT_TH >> log2CoeffGroupSize))
      {
         //coded_sbk_flag(iCGX,iCGY) shall be equal to 0
         assert(0 == uiSigCoeffGroupFlag[iCGBlkPos]);
      }
   }
#endif 
    // decode significant_coeff_flag
#if !VCEG_AZ07_CTX_RESIDUALCODING
    const Int patternSigCtx = TComTrQuant::calcPatternSigCtx(uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups);
#endif

    UInt uiBlkPos, uiSig, uiCtxSig;
    for( ; iScanPosSig >= iSubPos; iScanPosSig-- )
    {
      uiBlkPos  = codingParameters.scan[ iScanPosSig ];
      uiSig     = 0;

      if( uiSigCoeffGroupFlag[ iCGBlkPos ] )
      {
        if( iScanPosSig > iSubPos || iSubSet == 0  || numNonZero )
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          uiCtxSig  = TComTrQuant::getGrtZeroCtxInc( pcCoef, uiBlkPos, uiWidth, uiHeight, chType );
#else
          uiCtxSig  = TComTrQuant::getSigCtxInc( patternSigCtx, codingParameters, iScanPosSig, uiLog2BlockWidth, uiLog2BlockHeight, chType );
#endif
          m_pcTDecBinIf->decodeBin( uiSig, baseCtx[ uiCtxSig ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_map) );
        }
        else
        {
          uiSig = 1;
        }
      }
      pcCoef[ uiBlkPos ] = uiSig;
      if( uiSig )
      {
        pos[ numNonZero ] = uiBlkPos;
        numNonZero ++;
        if( lastNZPosInCG == -1 )
        {
          lastNZPosInCG = iScanPosSig;
        }
        firstNZPosInCG = iScanPosSig;
      }
    }

    if( numNonZero > 0 )
    {
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  bCbfCU = true;
#endif
     Bool signHidden = ( lastNZPosInCG - firstNZPosInCG >= SBH_THRESHOLD );
      absSum = 0;
#if VCEG_AZ07_CTX_RESIDUALCODING
      UInt uiBin;
#else
      const UInt uiCtxSet = getContextSetIndex(compID, iSubSet, (c1 == 0));
      c1 = 1;
      UInt uiBin;

      ContextModel *baseCtxMod = m_cCUOneSCModel.get( 0, 0 ) + (NUM_ONE_FLAG_CTX_PER_SET * uiCtxSet);
#endif
      Int absCoeff[1 << MLS_CG_SIZE];

      for ( Int i = 0; i < numNonZero; i++)
      {
        absCoeff[i] = 1;
      }
      Int numC1Flag = min(numNonZero, C1FLAG_NUMBER);
      Int firstC2FlagIdx = -1;

      for( Int idx = 0; idx < numC1Flag; idx++ )
      {
#if VCEG_AZ07_CTX_RESIDUALCODING
        if( idx || iSubSet != iLastScanSet)
        {
          ctxG1 = TComTrQuant::getGrtOneCtxInc( pcCoef, pos[idx], uiWidth, uiHeight, chType );
        }
        m_pcTDecBinIf->decodeBin( uiBin, greXCtx[ ctxG1 ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_gt1) );
#else
        m_pcTDecBinIf->decodeBin( uiBin, baseCtxMod[c1] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_gt1) );
#endif
        if( uiBin == 1 )
        {
#if !VCEG_AZ07_CTX_RESIDUALCODING
          c1 = 0;
#endif
          if (firstC2FlagIdx == -1)
          {
            firstC2FlagIdx = idx;
          }
          else //if a greater-than-one has been encountered already this group
          {
            escapeDataPresentInGroup = true;
          }
        }
#if VCEG_AZ07_CTX_RESIDUALCODING
        pcCoef[ pos[idx] ] += uiBin ;
#else
        else if( (c1 < 3) && (c1 > 0) )
        {
          c1++;
        }
#endif
        absCoeff[ idx ] = uiBin + 1;
      }

#if !VCEG_AZ07_CTX_RESIDUALCODING
      if (c1 == 0)
      {
        baseCtxMod = m_cCUAbsSCModel.get( 0, 0 ) + (NUM_ABS_FLAG_CTX_PER_SET * uiCtxSet);
#endif
        if ( firstC2FlagIdx != -1)
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          if( firstC2FlagIdx || iSubSet != iLastScanSet )
          {
            ctxG2 = TComTrQuant::getGrtTwoCtxInc( pcCoef, pos[firstC2FlagIdx], uiWidth, uiHeight, chType );
          }
          m_pcTDecBinIf->decodeBin( uiBin, greXCtx[ ctxG2 ] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_gt2) );
          pcCoef[ pos[firstC2FlagIdx] ] += uiBin ;
#else
          m_pcTDecBinIf->decodeBin( uiBin, baseCtxMod[0] RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_gt2) );
#endif
          absCoeff[ firstC2FlagIdx ] = uiBin + 2;
          if (uiBin != 0)
          {
            escapeDataPresentInGroup = true;
          }
        }
#if !VCEG_AZ07_CTX_RESIDUALCODING
      }
#endif

      escapeDataPresentInGroup = escapeDataPresentInGroup || (numNonZero > C1FLAG_NUMBER);

      const Bool alignGroup = escapeDataPresentInGroup && alignCABACBeforeBypass;

#if RExt__DECODER_DEBUG_BIT_STATISTICS
      TComCodingStatisticsClassType ctype_signs((alignGroup ? STATS__CABAC_BITS__ALIGNED_SIGN_BIT    : STATS__CABAC_BITS__SIGN_BIT   ), uiLog2BlockWidth, compID);
      TComCodingStatisticsClassType ctype_escs ((alignGroup ? STATS__CABAC_BITS__ALIGNED_ESCAPE_BITS : STATS__CABAC_BITS__ESCAPE_BITS), uiLog2BlockWidth, compID);
#endif

      if (alignGroup)
      {
        m_pcTDecBinIf->align();
      }
#if !VCEG_AZ07_CTX_RESIDUALCODING
      UInt coeffSigns;
      if ( signHidden && beValid )
      {
        m_pcTDecBinIf->decodeBinsEP( coeffSigns, numNonZero-1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_signs) );
        coeffSigns <<= 32 - (numNonZero-1);
      }
      else
      {
        m_pcTDecBinIf->decodeBinsEP( coeffSigns, numNonZero RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_signs) );
        coeffSigns <<= 32 - numNonZero;
      }
#endif
      Int iFirstCoeff2 = 1;
      if (escapeDataPresentInGroup)
      {
        for( Int idx = 0; idx < numNonZero; idx++ )
        {
          UInt baseLevel  = (idx < C1FLAG_NUMBER)? (2 + iFirstCoeff2) : 1;

          if( absCoeff[ idx ] == baseLevel)
          {
            UInt uiLevel;
#if VCEG_AZ07_CTX_RESIDUALCODING
            if( updateGolombRiceStatistics && iSubSet == iLastScanSet )
            {
              uiGoRiceParam = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;
            }
            else
            {
              uiGoRiceParam = TComTrQuant::getRemainCoeffCtxInc( pcCoef, pos[idx], uiWidth, uiHeight );
              if(bUseGolombRiceParameterAdaptation)
              {
                uiGoRiceParam = max(uiPrevGRParam, uiGoRiceParam);
              }
            }
            uiPrevGRParam = max(0, (Int)uiGoRiceParam - 1);
            xReadCoefRemainExGolomb( uiLevel, uiGoRiceParam, extendedPrecision, maxLog2TrDynamicRange RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_escs) );
            pcCoef[ pos[idx] ] = uiLevel + baseLevel;
#else
            xReadCoefRemainExGolomb( uiLevel, uiGoRiceParam, extendedPrecision, maxLog2TrDynamicRange RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_escs) );
#endif
            absCoeff[ idx ] = uiLevel + baseLevel;

#if !VCEG_AZ07_CTX_RESIDUALCODING
            if (absCoeff[idx] > (3 << uiGoRiceParam))
            {
              uiGoRiceParam = bUseGolombRiceParameterAdaptation ? (uiGoRiceParam + 1) : (std::min<UInt>((uiGoRiceParam + 1), 4));
            }
#endif
            if (updateGolombRiceStatistics)
            {
              const UInt initialGolombRiceParameter = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;

              if (uiLevel >= (3 << initialGolombRiceParameter))
              {
                currentGolombRiceStatistic++;
              }
              else if (((uiLevel * 2) < (1 << initialGolombRiceParameter)) && (currentGolombRiceStatistic > 0))
              {
                currentGolombRiceStatistic--;
              }

              updateGolombRiceStatistics = false;
            }

          }

          if(absCoeff[ idx ] >= 2)
          {
            iFirstCoeff2 = 0;
          }
        }
      }

#if VCEG_AZ07_CTX_RESIDUALCODING
      UInt coeffSigns;
#if JVET_C0024_QTBT
      if ( signHidden && beValid && uiWidth>=4 && uiHeight>=4)
#else
      if ( signHidden && beValid )
#endif
      {
        m_pcTDecBinIf->decodeBinsEP( coeffSigns, numNonZero-1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_signs) );
        coeffSigns <<= 32 - (numNonZero-1);
      }
      else
      {
        m_pcTDecBinIf->decodeBinsEP( coeffSigns, numNonZero RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype_signs) );
        coeffSigns <<= 32 - numNonZero;
      }
#endif
      for( Int idx = 0; idx < numNonZero; idx++ )
      {
        Int blkPos = pos[ idx ];
        // Signs applied later.
        pcCoef[ blkPos ] = absCoeff[ idx ];
        absSum += absCoeff[ idx ];

        if ( idx == numNonZero-1 && signHidden && beValid )
        {
          // Infer sign of 1st element.
          if (absSum&0x1)
          {
            pcCoef[ blkPos ] = -pcCoef[ blkPos ];
          }
        }
        else
        {
          Int sign = static_cast<Int>( coeffSigns ) >> 31;
          pcCoef[ blkPos ] = ( pcCoef[ blkPos ] ^ sign ) - sign;
          coeffSigns <<= 1;
        }
      }
    }
#if COM16_C806_EMT
    uiNumSig += numNonZero;
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    if( !pcCU->getTransformSkip( uiAbsPartIdx,compID) )
    {
      iNonZeroCoeffNonTs += numNonZero;
    }
#endif
  }

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  printSBACCoeffData(uiPosLastX, uiPosLastY, uiWidth, uiHeight, compID, uiAbsPartIdx, codingParameters.scanType, pcCoef);
#endif

#if COM16_C806_EMT
  if (!pcCU->getTransformSkip( uiAbsPartIdx, compID) && compID == COMPONENT_Y )
  {
    if ( pcCU->getEmtCuFlag( uiAbsPartIdx ) && pcCU->isIntra( uiAbsPartIdx ) )
    {
      if( uiNumSig > g_iEmtSigNumThr )
      {
        parseEmtTuIdx( pcCU, uiAbsPartIdx, rTu.GetTransformDepthTotal() ); 
      }
      else
      {
        pcCU->setEmtTuIdxSubParts( 0, uiAbsPartIdx, rTu.GetTransformDepthTotal() );
      }
    }
    if ( pcCU->getEmtCuFlag( uiAbsPartIdx ) && !pcCU->isIntra( uiAbsPartIdx ) )
    {
      parseEmtTuIdx( pcCU, uiAbsPartIdx, rTu.GetTransformDepthTotal() ); 
    }
  }
#endif

#if JVET_C0024_ITSKIP
#if VCEG_AZ08_USE_KLT
  if( pcCU->getKLTFlag(uiAbsPartIdx, compID) )
  {
    return;
  }
#endif
#if !ITH_C2_CNT_ZR_ZC
  Int maxX = uiPosLastX;
  Int maxY = uiPosLastY;
  for (Int i=0; i<uiWidth*uiHeight; i++)
  {
    UInt uiBlkPos = codingParameters.scan[ i ];
    if( uiBlkPosLast == uiBlkPos )
    {
      break;
    }

    Int x = uiBlkPos%uiWidth;
    Int y = uiBlkPos>>uiLog2BlockWidth;

    if (pcCoef[uiBlkPos]==0)
    {
      continue;
    }
    else 
    {
      if (x > maxX )
      {
        maxX = x;
      }
      if (y > maxY)
      {
        maxY = y;
      }
    }
  }

  pcCU->getTUSkipWidth(compID, uiAbsPartIdx) = uiWidth-1-maxX;
  pcCU->getTUSkipHeight(compID, uiAbsPartIdx) = uiHeight-1-maxY;
#else
  countZeroRowsZeroCols(uiWidth, uiHeight, pcCoef, uiWidth, &pcCU->getTUSkipWidth(compID, uiAbsPartIdx), &pcCU->getTUSkipHeight(compID, uiAbsPartIdx));
#endif
#endif

  return;
}
#if CNN_BASED_INLOOP_FILTER
Void TDecSbac::parseCNNLoopFilterTypeIdx(TComDataCU *pcCU, UInt&  ruiVal, CNNLFBlkParam* &ctuParam, ComponentID compIdx)
{
	UInt uiSymbol = 0;
	UInt uiCtx = pcCU->getCNNLoopFilterCtxFlag(compIdx, ctuParam);
	m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUCNNLoopFilterFlagSCModel.get(0, 0, uiCtx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__CNNLF_FLAG));
	ruiVal = (Int)uiSymbol;
}
Void TDecSbac::parseConvNetBlkParam(CNNLFBlkParam*& convNetBlkParam, TComDataCU *pcCU)
{
	UInt uiSymbol;
	for (Int compId = COMPONENT_Y; compId < MAX_NUM_COMPONENT; compId++)
	{
		const ComponentID compIdx = ComponentID(compId);
		const ComponentID firstCompOfChType = getFirstComponentOfChannel(toChannelType(compIdx));
		Bool& ctbParam = convNetBlkParam[pcCU->getCtuRsAddr()][compIdx];

		if (compIdx == firstCompOfChType)
		{
			parseCNNLoopFilterTypeIdx(pcCU, uiSymbol, convNetBlkParam, (ComponentID)compId);
			if (uiSymbol == 0) //OFF
			{
				ctbParam = false;
			}
			else if (uiSymbol == 1)
			{
				ctbParam = true;
			}
		}
		else //Cr, follow Cb  type
		{
			ctbParam = convNetBlkParam[pcCU->getCtuRsAddr()][COMPONENT_Cb];
		}
	}
}
#endif
Void TDecSbac::parseSaoMaxUvlc ( UInt& val, UInt maxSymbol )
{
  if (maxSymbol == 0)
  {
    val = 0;
    return;
  }

  UInt code;
  Int  i;
  m_pcTDecBinIf->decodeBinEP( code RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
  if ( code == 0 )
  {
    val = 0;
    return;
  }

  i=1;
  while (1)
  {
    m_pcTDecBinIf->decodeBinEP( code RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
    if ( code == 0 )
    {
      break;
    }
    i++;
    if (i == maxSymbol)
    {
      break;
    }
  }

  val = i;
}

Void TDecSbac::parseSaoUflc (UInt uiLength, UInt&  riVal)
{
  m_pcTDecBinIf->decodeBinsEP ( riVal, uiLength RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
}

Void TDecSbac::parseSaoMerge (UInt&  ruiVal)
{
  UInt uiCode;
  m_pcTDecBinIf->decodeBin( uiCode, m_cSaoMergeSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
  ruiVal = (Int)uiCode;
}

Void TDecSbac::parseSaoTypeIdx (UInt&  ruiVal)
{
  UInt uiCode;
  m_pcTDecBinIf->decodeBin( uiCode, m_cSaoTypeIdxSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
  if (uiCode == 0)
  {
    ruiVal = 0;
  }
  else
  {
    m_pcTDecBinIf->decodeBinEP( uiCode RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
    if (uiCode == 0)
    {
      ruiVal = 1;
    }
    else
    {
      ruiVal = 2;
    }
  }
}

Void TDecSbac::parseSaoSign(UInt& val)
{
  m_pcTDecBinIf->decodeBinEP ( val RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SAO) );
}

Void TDecSbac::parseSAOBlkParam (SAOBlkParam& saoBlkParam
                                , Bool* sliceEnabled
                                , Bool leftMergeAvail
                                , Bool aboveMergeAvail
                                , const BitDepths &bitDepths
                                )
{
  UInt uiSymbol;

  Bool isLeftMerge = false;
  Bool isAboveMerge= false;

  if(leftMergeAvail)
  {
    parseSaoMerge(uiSymbol); //sao_merge_left_flag
    isLeftMerge = (uiSymbol?true:false);
  }

  if( aboveMergeAvail && !isLeftMerge)
  {
    parseSaoMerge(uiSymbol); //sao_merge_up_flag
    isAboveMerge = (uiSymbol?true:false);
  }

  if(isLeftMerge || isAboveMerge) //merge mode
  {
    for (UInt componentIndex = 0; componentIndex < MAX_NUM_COMPONENT; componentIndex++)
    {
      saoBlkParam[componentIndex].modeIdc = (sliceEnabled[componentIndex]) ? SAO_MODE_MERGE : SAO_MODE_OFF;
      saoBlkParam[componentIndex].typeIdc = (isLeftMerge)?SAO_MERGE_LEFT:SAO_MERGE_ABOVE;
    }
  }
  else //new or off mode
  {
    for(Int compId=COMPONENT_Y; compId < MAX_NUM_COMPONENT; compId++)
    {
      const ComponentID compIdx=ComponentID(compId);
      const ComponentID firstCompOfChType = getFirstComponentOfChannel(toChannelType(compIdx));
      SAOOffset& ctbParam = saoBlkParam[compIdx];
#if O0043_BEST_EFFORT_DECODING
      const Int bitDepthOrig = bitDepths.stream[toChannelType(compIdx)];
      const Int forceBitDepthAdjust = bitDepthOrig - bitDepths.recon[toChannelType(compIdx)];
#else
      const Int bitDepthOrig = bitDepths.recon[toChannelType(compIdx)];
#endif
      const Int maxOffsetQVal=TComSampleAdaptiveOffset::getMaxOffsetQVal(bitDepthOrig);
      if(!sliceEnabled[compIdx])
      {
        //off
        ctbParam.modeIdc = SAO_MODE_OFF;
        continue;
      }

      //type
      if(compIdx == firstCompOfChType)
      {
        parseSaoTypeIdx(uiSymbol); //sao_type_idx_luma or sao_type_idx_chroma

        assert(uiSymbol ==0 || uiSymbol ==1 || uiSymbol ==2);

        if(uiSymbol ==0) //OFF
        {
          ctbParam.modeIdc = SAO_MODE_OFF;
        }
        else if(uiSymbol == 1) //BO
        {
          ctbParam.modeIdc = SAO_MODE_NEW;
          ctbParam.typeIdc = SAO_TYPE_START_BO;
        }
        else //2, EO
        {
          ctbParam.modeIdc = SAO_MODE_NEW;
          ctbParam.typeIdc = SAO_TYPE_START_EO;
        }

      }
      else //Cr, follow Cb SAO type
      {
        ctbParam.modeIdc = saoBlkParam[COMPONENT_Cb].modeIdc;
        ctbParam.typeIdc = saoBlkParam[COMPONENT_Cb].typeIdc;
      }

      if(ctbParam.modeIdc == SAO_MODE_NEW)
      {
        Int offset[4];
        for(Int i=0; i< 4; i++)
        {
          parseSaoMaxUvlc(uiSymbol, maxOffsetQVal ); //sao_offset_abs
          offset[i] = (Int)uiSymbol;
        }

        if(ctbParam.typeIdc == SAO_TYPE_START_BO)
        {
          for(Int i=0; i< 4; i++)
          {
            if(offset[i] != 0)
            {
              parseSaoSign(uiSymbol); //sao_offset_sign
              if(uiSymbol)
              {
#if O0043_BEST_EFFORT_DECODING
                offset[i] >>= forceBitDepthAdjust;
#endif
                offset[i] = -offset[i];
              }
            }
          }
          parseSaoUflc(NUM_SAO_BO_CLASSES_LOG2, uiSymbol ); //sao_band_position
          ctbParam.typeAuxInfo = uiSymbol;

          for(Int i=0; i<4; i++)
          {
            ctbParam.offset[(ctbParam.typeAuxInfo+i)%MAX_NUM_SAO_CLASSES] = offset[i];
          }

        }
        else //EO
        {
          ctbParam.typeAuxInfo = 0;

          if(firstCompOfChType == compIdx)
          {
            parseSaoUflc(NUM_SAO_EO_TYPES_LOG2, uiSymbol ); //sao_eo_class_luma or sao_eo_class_chroma
            ctbParam.typeIdc += uiSymbol;
          }
          else
          {
            ctbParam.typeIdc = saoBlkParam[firstCompOfChType].typeIdc;
          }
          ctbParam.offset[SAO_CLASS_EO_FULL_VALLEY] = offset[0];
          ctbParam.offset[SAO_CLASS_EO_HALF_VALLEY] = offset[1];
          ctbParam.offset[SAO_CLASS_EO_PLAIN      ] = 0;
          ctbParam.offset[SAO_CLASS_EO_HALF_PEAK  ] = -offset[2];
          ctbParam.offset[SAO_CLASS_EO_FULL_PEAK  ] = -offset[3];
        }
      }
    }
  }
}

/**
 - Initialize our contexts from the nominated source.
 .
 \param pSrc Contexts to be copied.
 */
Void TDecSbac::xCopyContextsFrom( const TDecSbac* pSrc )
{
  memcpy(m_contextModels, pSrc->m_contextModels, m_numContextModels*sizeof(m_contextModels[0]));
  memcpy(m_golombRiceAdaptationStatistics, pSrc->m_golombRiceAdaptationStatistics, (sizeof(UInt) * RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS));
}

Void TDecSbac::xCopyFrom( const TDecSbac* pSrc )
{
  m_pcTDecBinIf->copyState( pSrc->m_pcTDecBinIf );
  xCopyContextsFrom( pSrc );
}

Void TDecSbac::load ( const TDecSbac* pSrc )
{
  xCopyFrom(pSrc);
}

Void TDecSbac::loadContexts ( const TDecSbac* pSrc )
{
  xCopyContextsFrom(pSrc);
}

/** Performs CABAC decoding of the explicit RDPCM mode
 * \param rTu current TU data structure
 * \param compID component identifier
 */
Void TDecSbac::parseExplicitRdpcmMode( TComTU &rTu, ComponentID compID )
{
  TComDataCU* cu = rTu.getCU();
  const UInt absPartIdx=rTu.GetAbsPartIdxTU(compID);
#if !JVET_C0046_OMIT_ASSERT_ERDPCM
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt tuHeight = g_aucConvertToBit[rect.height];
  const UInt tuWidth  = g_aucConvertToBit[rect.width];
#endif  
  UInt code = 0;

#if !JVET_C0046_OMIT_ASSERT_ERDPCM
  assert(tuHeight == tuWidth);
#endif 

#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
  const TComCodingStatisticsClassType ctype(STATS__EXPLICIT_RDPCM_BITS, g_aucConvertToBit[cu->getSlice()->getSPS()->getCTUSize()>>rTu.GetTransformDepthTotal()]+2);
#else
  const TComCodingStatisticsClassType ctype(STATS__EXPLICIT_RDPCM_BITS, g_aucConvertToBit[cu->getSlice()->getSPS()->getMaxCUWidth()>>rTu.GetTransformDepthTotal()]+2);
#endif
#endif

  m_pcTDecBinIf->decodeBin(code, m_explicitRdpcmFlagSCModel.get (0, toChannelType(compID), 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));

  if(code == 0)
  {
    cu->setExplicitRdpcmModePartRange( RDPCM_OFF, compID, absPartIdx, rTu.GetAbsPartIdxNumParts(compID));
  }
  else
  {
    m_pcTDecBinIf->decodeBin(code, m_explicitRdpcmDirSCModel.get (0, toChannelType(compID), 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(ctype));
    if(code == 0)
    {
      cu->setExplicitRdpcmModePartRange( RDPCM_HOR, compID, absPartIdx, rTu.GetAbsPartIdxNumParts(compID));
    }
    else
    {
      cu->setExplicitRdpcmModePartRange( RDPCM_VER, compID, absPartIdx, rTu.GetAbsPartIdxNumParts(compID));
    }
  }
}

#if ALF_HM3_REFACTOR
Void TDecSbac::parseAlfCtrlDepth( UInt& ruiAlfCtrlDepth , UInt uiMaxTotalCUDepth )
{
  UInt uiSymbol;
  xReadUnaryMaxSymbol( uiSymbol, m_cALFUvlcSCModel.get( 0 ), 1, uiMaxTotalCUDepth - 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_CTRL_DEPTH));
  ruiAlfCtrlDepth = uiSymbol;
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlDepth\n" )
}

#if !JVET_C0024_QTBT
Void TDecSbac::parseAlfCtrlFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiMaxAlfCtrlDepth )
{
  if( uiDepth > uiMaxAlfCtrlDepth && !pcCU->isFirstAbsZorderIdxInDepth( uiAbsPartIdx, uiMaxAlfCtrlDepth ) )
  {
    return;
  }

  UInt uiSymbol;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUAlfCtrlFlagSCModel.get( 0, 0, pcCU->getCtxAlfCtrlFlag( uiAbsPartIdx ) ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_CTRL_FLAG) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlFlag\n" )

  if( uiDepth > uiMaxAlfCtrlDepth )
  {
    pcCU->setAlfCtrlFlagSubParts( uiSymbol, uiAbsPartIdx, uiMaxAlfCtrlDepth );
  }
  else
  {
    pcCU->setAlfCtrlFlagSubParts( uiSymbol, uiAbsPartIdx, uiDepth );
  }
}
#endif
#if JVET_C0038_GALF  
Void TDecSbac::parseALFTruncBinVal( UInt& ruiSymbol, UInt uiMaxSymbol )
{
  xReadTruncBinCode(ruiSymbol, uiMaxSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP));
}
//0: no pred; 1: all same index; 2: diff index for each variance index
Void TDecSbac::parseALFPrevFiltType( UInt& uiCode)
{
  xReadEpExGolomb (uiCode, 0 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_FLAG) );
}
Void TDecSbac::parseALFPrevFiltFlag( UInt& uiCode)
{
   m_pcTDecBinIf->decodeBinEP( uiCode RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_FLAG) );
}
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
Void TDecSbac::parseCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
	UInt uiSymbol;
	m_pcTDecBinIf->decodeBin(uiSymbol, m_CUBgFlagSCModel.get(0, 0, 0)RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__BG_FLAG));
	pcCU->setLCUBgFlag(uiSymbol ? true : false);
}
#endif

Void TDecSbac::parseAlfFlag (UInt& ruiVal)
{
  UInt uiSymbol;
#if JVET_C0038_GALF
  m_pcTDecBinIf->decodeBinEP( uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_FLAG));
#else
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cALFFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_FLAG));
#endif
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfFlag\n" )

  ruiVal = uiSymbol;
}

Void TDecSbac::parseAlfFlagNum( UInt& ruiVal, UInt minValue, UInt depth )
{
  UInt uiLength = 0;
  UInt maxValue = (minValue << (depth*2));
  UInt temp = maxValue - minValue;
  for(UInt i=0; i<32; i++)
  {
    if(temp&0x1)
    {
      uiLength = i+1;
    }
    temp = (temp >> 1);
  }
  ruiVal = 0;
  UInt uiBit;
  if(uiLength)
  {
    while( uiLength-- )
    {
      m_pcTDecBinIf->decodeBinEP( uiBit RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_FLAGNUM));
      ruiVal += uiBit << uiLength;
    }
  }
  else
  {
    ruiVal = 0;
  }
  ruiVal += minValue;
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfFlagNum\n" )
}

Void TDecSbac::parseAlfCtrlFlag( UInt &ruiAlfCtrlFlag )
{
  UInt uiSymbol;
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUAlfCtrlFlagSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_CTRL_FLAG));
  ruiAlfCtrlFlag = uiSymbol;
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlFlag\n" )
}

Void TDecSbac::parseAlfUvlc (UInt& ruiVal)
{
  UInt uiCode;
  Int  i;
#if JVET_C0038_GALF
  m_pcTDecBinIf->decodeBinEP( uiCode  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_UVLC));
  if ( uiCode == 0 )
  {
    ruiVal = 0;
    return;
  }

  i=1;
  while (1)
  {
    m_pcTDecBinIf->decodeBinEP( uiCode  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_UVLC));
    if ( uiCode == 0 ) break;
    i++;
  }
#else
  m_pcTDecBinIf->decodeBin( uiCode, m_cALFUvlcSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_UVLC));
  if ( uiCode == 0 )
  {
    ruiVal = 0;
    return;
  }

  i=1;
  while (1)
  {
    m_pcTDecBinIf->decodeBin( uiCode, m_cALFUvlcSCModel.get( 0, 0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_UVLC));
    if ( uiCode == 0 ) break;
    i++;
  }
#endif
  ruiVal = i;
}

Void TDecSbac::parseAlfSvlc (Int&  riVal)
{
  UInt uiCode;
  Int  iSign;
  Int  i;
#if JVET_C0038_GALF
  m_pcTDecBinIf->decodeBinEP( uiCode  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));
  if ( uiCode == 0 )
  {
    riVal = 0;
    return;
  }

  // read sign
  m_pcTDecBinIf->decodeBinEP( uiCode RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));

  if ( uiCode == 0 ) iSign =  1;
  else               iSign = -1;

  // read magnitude
  i=1;
  while (1)
  {
    m_pcTDecBinIf->decodeBinEP( uiCode  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));
    if ( uiCode == 0 ) break;
    i++;
  }
#else
  m_pcTDecBinIf->decodeBin( uiCode, m_cALFSvlcSCModel.get( 0, 0, 0 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));

  if ( uiCode == 0 )
  {
    riVal = 0;
    return;
  }

  // read sign
  m_pcTDecBinIf->decodeBin( uiCode, m_cALFSvlcSCModel.get( 0, 0, 1 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));

  if ( uiCode == 0 ) iSign =  1;
  else               iSign = -1;

  // read magnitude
  i=1;
  while (1)
  {
    m_pcTDecBinIf->decodeBin( uiCode, m_cALFSvlcSCModel.get( 0, 0, 2 ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__ALF_SVLC));
    if ( uiCode == 0 ) break;
    i++;
  }
#endif
  riVal = i*iSign;
}
#endif

#if COM16_C806_EMT
#if FAST_EMT_v1
Void TDecSbac::parseEmtTuIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
	UChar trIdx = 0;

#if JVET_C0024_QTBT
	if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTRA_MAX_CU)
#else
	if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU)
#endif
	{
		//UInt uiSymbol1 = 0, uiSymbol2 = 0;

		UInt uiSymbol = 0;

		m_pcTDecBinIf->decodeBin(uiSymbol, m_cEMTIntraidxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
		if (uiSymbol)
		{
			m_pcTDecBinIf->decodeBin(uiSymbol, m_cEMTIntraidxSCModel.get(0, 0, 3) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
			uiSymbol++;
		}

		trIdx = uiSymbol;
		//	m_pcTDecBinIf->decodeBin( uiSymbol1, m_cEmtTuIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
		//	m_pcTDecBinIf->decodeBin( uiSymbol2, m_cEmtTuIdxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
		//	trIdx = (uiSymbol2 << 1) | uiSymbol1; 
		//		printf("%d\t%d\t%d\n", pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx),trIdx);
	}
#if JVET_C0024_QTBT
	if (!pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTER_MAX_CU)
#else
	if (!pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU)
#endif
	{
#if FAST_EMT_v1 //EMT_SIGNAL_INTER
		UInt uiSymbol = 0;
		m_pcTDecBinIf->decodeBin(uiSymbol, m_cEMTIntraidxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
		if (uiSymbol)
		{
			m_pcTDecBinIf->decodeBin(uiSymbol, m_cEMTIntraidxSCModel.get(0, 0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
			uiSymbol++;
		}
		trIdx = uiSymbol;
#else
		UInt uiSymbol1 = 0;
		UInt uiSymbol2 = 0;
		m_pcTDecBinIf->decodeBin(uiSymbol1, m_cEmtTuIdxSCModel.get(0, 0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
		m_pcTDecBinIf->decodeBin(uiSymbol2, m_cEmtTuIdxSCModel.get(0, 0, 3) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX));
		trIdx = (uiSymbol2 << 1) | uiSymbol1;
#endif
	}

	pcCU->setEmtTuIdxSubParts(trIdx, uiAbsPartIdx, uiDepth);
}
#else
Void TDecSbac::parseEmtTuIdx( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  UChar trIdx = 0;

#if JVET_C0024_QTBT
  if ( pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx)<= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTRA_MAX_CU  )
#else
  if ( pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU  )
#endif
  {
    UInt uiSymbol1 = 0, uiSymbol2 = 0;
    m_pcTDecBinIf->decodeBin( uiSymbol1, m_cEmtTuIdxSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
    m_pcTDecBinIf->decodeBin( uiSymbol2, m_cEmtTuIdxSCModel.get(0, 0, 1) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
    trIdx = (uiSymbol2 << 1) | uiSymbol1; 
  }
#if JVET_C0024_QTBT
  if( !pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx)<= EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTER_MAX_CU  )
#else
  if( !pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU )
#endif
  {
    UInt uiSymbol1 = 0, uiSymbol2 = 0;
    m_pcTDecBinIf->decodeBin( uiSymbol1, m_cEmtTuIdxSCModel.get(0, 0, 2) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
    m_pcTDecBinIf->decodeBin( uiSymbol2, m_cEmtTuIdxSCModel.get(0, 0, 3) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_TU_INDX) );
    trIdx = (uiSymbol2 << 1) | uiSymbol1; 
  }

  pcCU->setEmtTuIdxSubParts( trIdx, uiAbsPartIdx, uiDepth );
}
#endif
Void TDecSbac::parseEmtCuFlag  ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bParseCuFlag )
{
#if JVET_C0024_QTBT
  if (uiDepth >= NUM_EMT_CU_FLAG_CTX)
  {
    uiDepth = NUM_EMT_CU_FLAG_CTX - 1;
  }
#else
  assert( uiDepth < NUM_EMT_CU_FLAG_CTX );
#endif

  pcCU->setEmtCuFlagSubParts( 0, uiAbsPartIdx, uiDepth );
  pcCU->setEmtTuIdxSubParts( DCT2_EMT, uiAbsPartIdx, uiDepth );

#if JVET_C0024_QTBT
  if ( pcCU->isIntra( uiAbsPartIdx ) && bParseCuFlag && pcCU->getWidth(uiAbsPartIdx)<= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTRA_MAX_CU 
      && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
#else
  if ( pcCU->isIntra( uiAbsPartIdx ) && bParseCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
#endif
  {
    UInt uiCuFlag = 0;
    m_pcTDecBinIf->decodeBin( uiCuFlag, m_cEmtCuFlagSCModel.get(0, 0, uiDepth) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_CU_FLAG) );
    pcCU->setEmtCuFlagSubParts( uiCuFlag, uiAbsPartIdx, uiDepth );
  }
#if JVET_C0024_QTBT
  if( !pcCU->isIntra( uiAbsPartIdx ) && bParseCuFlag && pcCU->getWidth(uiAbsPartIdx)<= EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTER_MAX_CU
      && pcCU->getSlice()->getSPS()->getUseInterEMT() )
#else
  if( !pcCU->isIntra( uiAbsPartIdx ) && bParseCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU && pcCU->getSlice()->getSPS()->getUseInterEMT() )
#endif
  {
    UInt uiCuFlag = 0;
    m_pcTDecBinIf->decodeBin( uiCuFlag, m_cEmtCuFlagSCModel.get(0, 0, uiDepth) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__EMT_CU_FLAG));
    pcCU->setEmtCuFlagSubParts( uiCuFlag, uiAbsPartIdx, uiDepth );
  }
}
#endif

#if VCEG_AZ07_INIT_PREVFRAME
Void TDecSbac::loadContextsFromPrev (TComStats* apcStats, SliceType eSliceType, Int iQPIdx, Bool bFromGloble, Int iQPIdxRst, Bool bAfterLastISlice )
{
  if(bFromGloble)
  {
#if VCEG_AZ07_INIT_PREVFRAME_FIX
    if(iQPIdx==-1 || !apcStats->aaQPUsed[eSliceType][iQPIdxRst].resetInit)
#else
    if(iQPIdx==-1 || (bAfterLastISlice && !apcStats->aaQPUsed[eSliceType][iQPIdxRst].resetInit))
#endif
    {
      return;
    }
    Int iCtxNr = getCtxNumber();
    for(UInt i = 0; i < iCtxNr; i++)
    {
#if VCEG_AZ07_BAC_ADAPT_WDOW  || VCEG_AZ05_MULTI_PARAM_CABAC
      m_contextModels[i].setState(apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i]);
#else
      m_contextModels[i].setState( (UChar) apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i] );
#endif
    }
  }
  else
  {
    Int iCtxNr = getCtxNumber();
    for(UInt i = 0; i < iCtxNr; i++)
    {
#if VCEG_AZ07_BAC_ADAPT_WDOW  || VCEG_AZ05_MULTI_PARAM_CABAC
      apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i] = m_contextModels[i].getState();
#else
      apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i] = (UShort) m_contextModels[i].getOrigState();
#endif
    }
  }
}
#endif

#if COM16_C1016_AFFINE
/** parse affine flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \param uiDepth
 * \returns Void
 */
Void TDecSbac::parseAffineFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx )
{
  if( pcCU->getSlice()->isIntra() )
  {
    return;
  }

  UInt uiSymbol = 0;
  UInt uiCtxAffine = pcCU->getCtxAffineFlag( uiAbsPartIdx );
  m_pcTDecBinIf->decodeBin( uiSymbol, m_cCUAffineFlagSCModel.get( 0, 0, uiCtxAffine ) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__AFFINE_FLAG) );

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tAffineFlag" );
  DTRACE_CABAC_T( "\tuiCtxAffine: ");
  DTRACE_CABAC_V( uiCtxAffine );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");

  pcCU->setAffineFlagSubParts( uiSymbol ? true : false, uiAbsPartIdx, uiPuIdx, uiDepth );
#if JVECT_C0062_AFFINE_SIX_PARAM
  pcCU->setAffineParamFlagSubParts(false, uiAbsPartIdx, uiPuIdx, uiDepth);
#endif
}

#if JVECT_C0062_AFFINE_SIX_PARAM
Void TDecSbac::parseAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPuIdx)
{
  if (pcCU->getSlice()->isIntra())
  {
    return;
  }

  UInt uiSymbol = 0;
  UInt uiCtxAffineParam = pcCU->getCtxAffineParamFlag(uiAbsPartIdx);
 
  m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUAffineParamFlagSCModel.get(0, 0, uiCtxAffineParam) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__AFFINE_FLAG));

  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tAffineFlag");
  DTRACE_CABAC_T("\tuiCtxAffine: ");
  DTRACE_CABAC_V(uiCtxAffine);
  DTRACE_CABAC_T("\tuiSymbol: ");
  DTRACE_CABAC_V(uiSymbol);
  DTRACE_CABAC_T("\n");

  pcCU->setAffineParamFlagSubParts(uiSymbol ? true : false, uiAbsPartIdx, uiPuIdx, uiDepth);
}
#endif

Void TDecSbac::parseAffineMvd( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth, RefPicList eRefList )
{
  TComMv acMvd[3];

#if JVECT_C0062_AFFINE_SIX_PARAM
  Int width = pcCU->getWidth(uiAbsPartIdx);
  Int height = pcCU->getHeight(uiAbsPartIdx);
  Bool bRecShape = (height / width >= 2);

  Bool b6Param = pcCU->getAffineParamFlag(uiAbsPartIdx);
  Int iVerIdx = b6Param ? 3 : 2;
  for (Int i = 0; i < iVerIdx; i++)
#else
  for ( Int i=0; i<2; i++ )
#endif
  {
    UInt uiSymbol;
    UInt uiHorAbs;
    UInt uiVerAbs;
    UInt uiHorSign = 0;
    UInt uiVerSign = 0;
    ContextModel *pCtx = m_cCUMvdSCModel.get( 0 );

    if(pcCU->getSlice()->getMvdL1ZeroFlag() && eRefList == REF_PIC_LIST_1 && pcCU->getInterDir(uiAbsPartIdx)==3)
    {
      uiHorAbs=0;
      uiVerAbs=0;
    }
    else
    {
      m_pcTDecBinIf->decodeBin( uiHorAbs, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
      m_pcTDecBinIf->decodeBin( uiVerAbs, *pCtx  RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );

      const Bool bHorAbsGr0 = uiHorAbs != 0;
      const Bool bVerAbsGr0 = uiVerAbs != 0;
      pCtx++;

      if( bHorAbsGr0 )
      {
        m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
        uiHorAbs += uiSymbol;
      }

      if( bVerAbsGr0 )
      {
        m_pcTDecBinIf->decodeBin( uiSymbol, *pCtx RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD) );
        uiVerAbs += uiSymbol;
      }

      if( bHorAbsGr0 )
      {
        if( 2 == uiHorAbs )
        {
          xReadEpExGolomb( uiSymbol, 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
          uiHorAbs += uiSymbol;
        }

        m_pcTDecBinIf->decodeBinEP( uiHorSign RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
      }

      if( bVerAbsGr0 )
      {
        if( 2 == uiVerAbs )
        {
          xReadEpExGolomb( uiSymbol, 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
          uiVerAbs += uiSymbol;
        }

        m_pcTDecBinIf->decodeBinEP( uiVerSign RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MVD_EP) );
      }
    }

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    TComMv cMv( uiHorSign ? -Int( uiHorAbs ): uiHorAbs, uiVerSign ? -Int( uiVerAbs ) : uiVerAbs );
    cMv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#else
    const TComMv cMv( uiHorSign ? -Int( uiHorAbs ): uiHorAbs, uiVerSign ? -Int( uiVerAbs ) : uiVerAbs );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    acMvd[(bRecShape && i == 1 && b6Param == 0) ? (i + 1) : i] = cMv;
#else
    acMvd[i] = cMv;
#endif
  }

  pcCU->setAllAffineMvd( uiAbsPartIdx, uiPartIdx, acMvd, eRefList, uiDepth );
  return;
}
#endif

#if MRG_CONSIDER_COMPLEX
Void TDecSbac::parseMergeComplexFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiPUIdx)
{
  UInt uiSymbol = 0;
  UInt uiCtx = pcCU->getCtxMergeComplexFlag(uiAbsPartIdx);
  m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUMergeComplexFlagExtSCModel.get(0, 0, uiCtx) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MRG_COMPLEX_FLAG));
#if JVET_C0024_QTBT
  pcCU->setMergeComplexFlagSubParts(uiSymbol ? true : false, uiAbsPartIdx);
#else 
  pcCU->setMergeComplexFlagSubParts(uiSymbol ? true : false, uiAbsPartIdx, uiPUIdx, uiDepth);
#endif
  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tMergeComplexFlag: ");
  DTRACE_CABAC_V(uiSymbol);
  DTRACE_CABAC_T("\n");
}

Void TDecSbac::parseMergeComplexIndex(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt& ruiMergeComplexIndex)
{
  UInt uiUnaryIdx = 0;
  UChar numCand = MRG_COMPLEX_NUM_CANDS_REAL;
  Int n = 0;
  Int numBits = 0;
  if (numCand == 2)
  {
    numBits = 1;
  }
  else if ((numCand > 2) && (numCand <= 4))
  {
    numBits = 2;
  }
  else if ((numCand > 4) && (numCand <= 8))
  {
    numBits = 3;
  }
  else if ((numCand > 8) && (numCand <= 16))
  {
    numBits = 4;
  }

  if (numCand > 1)
  {
    UInt uiCtx = 0;
    UInt dec_min_val = 0;
    UInt max_val = numCand - 1;
    n = numBits - 1;

    while (n >= 0)
    {
      Int diff = max_val - dec_min_val;

      if (diff <= 0)
      {
        n = min(n, -1);
      }
      else if (diff <= 1)
      {
        n = min(n, 0);
      }
      else if (diff <= 3)
      {
        n = min(n, 1);
      }
      else if (diff <= 7)
      {
        n = min(n, 2);
      }
      else
      {
        n = min(n, 3);
      }

      if (n >= 0)
      {
        UInt uiSymbol = 0;

        if (n == numBits - 1)
          m_pcTDecBinIf->decodeBin(uiSymbol, m_cCUMergeComplexIdxExtSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MRG_COMPLEX_INDX));
        else
          m_pcTDecBinIf->decodeBinEP(uiSymbol RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__MRG_COMPLEX_INDX));
        dec_min_val = dec_min_val | (uiSymbol << n);
      }
      n--;
    }
    uiUnaryIdx = dec_min_val;
  }
  ruiMergeComplexIndex = uiUnaryIdx;

  DTRACE_CABAC_VL(g_nSymbolCounter++)
    DTRACE_CABAC_T("\tparseMergeComplexIndex()")
    DTRACE_CABAC_T("\tuiMRGComplexIdx= ")
    DTRACE_CABAC_V(ruiMergeComplexIndex)
    DTRACE_CABAC_T("\n")
}
#endif

#if JVET_E0077_ENHANCED_LM
Int TDecSbac::parseLMMode(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
#if RExt__DECODER_DEBUG_BIT_STATISTICS
#if JVET_C0024_QTBT
    const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize() >> uiDepth] + 2, CHANNEL_TYPE_CHROMA);
#else
    const TComCodingStatisticsClassType ctype(STATS__CABAC_BITS__INTRA_DIR_ANG, g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth() >> uiDepth] + 2, CHANNEL_TYPE_CHROMA);
#endif
#endif

    UInt uiSymbol = -1;
    Int aiLMModeList[10];
    Int iSymbolNum = pcCU->getLMSymbolList(aiLMModeList, uiAbsPartIdx);

    xReadUnaryMaxSymbol(uiSymbol, &m_cCUChromaPredSCModel.get(0, 0, 1), 1, iSymbolNum - 1 RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__INTRA_DIR_ANG));
    return aiLMModeList[uiSymbol];
}
#endif

#if RRC_NS
UInt TDecSbac::parseNsFlag()
{
	UInt uiCode;
	m_pcTDecBinIf->decodeBin(uiCode, m_cNsFlagSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__NS_FLAG));
	return uiCode;
}

UInt TDecSbac::parseAdaptiveNsFlag()
{
	UInt uiCode;
	m_pcTDecBinIf->decodeBinEP(uiCode RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__NS_ADAPT_FLAG));
	return uiCode;
}

UInt TDecSbac::parseNsFWholeFrameFlag()
{
	UInt uiCode;
	m_pcTDecBinIf->decodeBinEP(uiCode RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__NS_WHOLE_FRAME_FLAG));
	return uiCode;
}

Int TDecSbac::parseNsBlockSize()
{
	UInt first_bit, second_bit;
	m_pcTDecBinIf->decodeBinEP(first_bit RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__NS_BLOCK_SIZE));
	m_pcTDecBinIf->decodeBinEP(second_bit RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__NS_BLOCK_SIZE));

	return (first_bit << 1) | second_bit;
}
#endif

#if RRC_SELECTABLE_IF
Int TDecSbac::parseAdaptiveIfIndex()
{
  UInt uiCode = 0;
  m_pcTDecBinIf->decodeBin(uiCode, m_cCUAdaptiveIfIndex.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__SIF));
  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tSIF index: ");
  DTRACE_CABAC_V(uiCode);
  DTRACE_CABAC_T("\n");
  return uiCode;
}
#endif

#if RRC_INTRA_SHARPENING
UInt TDecSbac::parseIntraSharpFlag()
{
  UInt uiCode = 0;
  m_pcTDecBinIf->decodeBin(uiCode, m_cSharpFlagSCModel.get(0, 0, 0) RExt__DECODER_DEBUG_BIT_STATISTICS_PASS_OPT_ARG(STATS__CABAC_BITS__INTRA_SHARP_FLAG));

  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tSharpFlag: ");
  DTRACE_CABAC_V(uiCode);
  DTRACE_CABAC_T("\n");
  return uiCode;
}
#endif

//! \}
