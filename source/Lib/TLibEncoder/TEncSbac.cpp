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

/** \file     TEncSbac.cpp
    \brief    SBAC encoder class
*/

#include "TEncTop.h"
#include "TEncSbac.h"
#include "TLibCommon/TComTU.h"

#include <map>
#include <algorithm>
#if VCEG_AZ07_CTX_RESIDUALCODING
#include <math.h>
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
#include "../TLibCommon/Debug.h"
#endif


//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TEncSbac::TEncSbac()
// new structure here
: m_pcBitIf                            ( NULL )
, m_pcBinIf                            ( NULL )
, m_numContextModels                   ( 0 )
, m_cCUSplitFlagSCModel                ( 1,             1,                      NUM_SPLIT_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#if JVET_C0024_QTBT
, m_cBTSplitFlagSCModel                ( 1,             1,                      NUM_BTSPLIT_MODE_CTX                 , m_contextModels + m_numContextModels, m_numContextModels )
#endif
, m_cCUSkipFlagSCModel                 ( 1,             1,                      NUM_SKIP_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#if BARC_UP_SAMPLE_RDO
, m_cCUUpCNNFlagSCModel                ( 1,             1,                      NUM_BARC_UP_FILT_RDO_FLAG_CTX        , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if BARC
, m_cCUBARCFlagSCModel                 ( 1,             1,                      NUM_BARC_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ05_INTRA_MPI
, m_cMPIIdxSCModel                     ( 1,             1,                      NUM_MPI_CTX                          , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if COM16_C1046_PDPC_INTRA
 , m_cPDPCIdxSCModel                   ( 1,             1,                      NUM_PDPC_CTX                         , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if DIMD_INTRA_PRED
  , m_cDIMDEnabledSCModel(1, 1, NUM_DIMD_CTX, m_contextModels + m_numContextModels, m_numContextModels)
  , m_cDIMDNoBTFlagSCModel(1, 1, NUM_DIMD_NO_BT_CTX, m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_DWDIP
 , m_cDWDIPIdxSCModel                  ( 1,             1,                      NUM_DWDIP_CTX                        , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_DWDIP_MLR
 , m_cDWDIP_GPIdxSCModel               ( 1,             1,                      NUM_DWDIP_GP_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
, m_cROTidxSCModel                     ( 1,             1,                      NUM_ROT_TR_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if  FAST_EMT_v1
, m_cEMTIntraidxSCModel                (1,             1,                        NUM_EMT_TR_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cCUMergeFlagExtSCModel             ( 1,             1,                      NUM_MERGE_FLAG_EXT_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMergeIdxExtSCModel              ( 1,             1,                      NUM_MERGE_IDX_EXT_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#if VCEG_AZ07_FRUC_MERGE
, m_cCUFRUCMgrModeSCModel              ( 1,             1,                      NUM_FRUCMGRMODE_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUFRUCMESCModel                   ( 1,             1,                      NUM_FRUCME_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
#if THIRD_FRUC_MODE
, m_cCUFRUCMESCModel_uni               ( 1,             1,                      NUM_FRUCME_UNI_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#endif
, m_cCUPartSizeSCModel                 ( 1,             1,                      NUM_PART_SIZE_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUPredModeSCModel                 ( 1,             1,                      NUM_PRED_MODE_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUIntraPredSCModel                ( 1,             1,                      NUM_INTRA_PREDICT_CTX                , m_contextModels + m_numContextModels, m_numContextModels)
#if MLIP_EXTEND_REF_LINE
, m_cCUExternRefSCModel(1, 1, NUM_EXTERN_REF_CTX, m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cCUChromaPredSCModel               ( 1,             1,                      NUM_CHROMA_PRED_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUDeltaQpSCModel                  ( 1,             1,                      NUM_DELTA_QP_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUInterDirSCModel                 ( 1,             1,                      NUM_INTER_DIR_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCURefPicSCModel                   ( 1,             1,                      NUM_REF_NO_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMvdSCModel                      ( 1,             1,                      NUM_MV_RES_CTX                       , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUQtCbfSCModel                    ( 1,             NUM_QT_CBF_CTX_SETS,    NUM_QT_CBF_CTX_PER_SET               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUTransSubdivFlagSCModel          ( 1,             1,                      NUM_TRANS_SUBDIV_FLAG_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUQtRootCbfSCModel                ( 1,             1,                      NUM_QT_ROOT_CBF_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUSigCoeffGroupSCModel            ( 1,             2,                      NUM_SIG_CG_FLAG_CTX                  , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUSigSCModel                      ( 1,             1,                      NUM_SIG_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCuCtxLastX                        ( 1,             NUM_CTX_LAST_FLAG_SETS, NUM_CTX_LAST_FLAG_XY                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCuCtxLastY                        ( 1,             NUM_CTX_LAST_FLAG_SETS, NUM_CTX_LAST_FLAG_XY                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUOneSCModel                      ( 1,             1,                      NUM_ONE_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#if !VCEG_AZ07_CTX_RESIDUALCODING
, m_cCUAbsSCModel                      ( 1,             1,                      NUM_ABS_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cMVPIdxSCModel                     ( 1,             1,                      NUM_MVP_IDX_CTX                      , m_contextModels + m_numContextModels, m_numContextModels)
#if CNN_BASED_INLOOP_FILTER
, m_cCUCNNLoopFilterFlagSCModel        ( 1,             1,                      NUM_CNNLF_FLAG_CTX                   , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cSaoMergeSCModel                   ( 1,             1,                      NUM_SAO_MERGE_FLAG_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
, m_cSaoTypeIdxSCModel                 ( 1,             1,                      NUM_SAO_TYPE_IDX_CTX                 , m_contextModels + m_numContextModels, m_numContextModels)
, m_cTransformSkipSCModel              ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_TRANSFORMSKIP_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
#if VCEG_AZ08_KLT_COMMON
, m_cKLTFlagSCModel                    ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_KLT_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_CUTransquantBypassFlagSCModel      ( 1,             1,                      NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX    , m_contextModels + m_numContextModels, m_numContextModels)
, m_explicitRdpcmFlagSCModel           ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_EXPLICIT_RDPCM_FLAG_CTX          , m_contextModels + m_numContextModels, m_numContextModels)
, m_explicitRdpcmDirSCModel            ( 1,             MAX_NUM_CHANNEL_TYPE,   NUM_EXPLICIT_RDPCM_DIR_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCrossComponentPredictionSCModel   ( 1,             1,                      NUM_CROSS_COMPONENT_PREDICTION_CTX   , m_contextModels + m_numContextModels, m_numContextModels)
, m_ChromaQpAdjFlagSCModel             ( 1,             1,                      NUM_CHROMA_QP_ADJ_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
, m_ChromaQpAdjIdcSCModel              ( 1,             1,                      NUM_CHROMA_QP_ADJ_IDC_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
#if COM16_C806_OBMC
, m_cCUOBMCFlagSCModel                 ( 1,             1,                      NUM_OBMC_FLAG_CTX                    , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ07_IMV
, m_cCUiMVFlagSCModel                  ( 1,             1,                      NUM_IMV_FLAG_CTX                     , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if VCEG_AZ06_IC
, m_cCUICFlagSCModel                   ( 1,             1,                      NUM_IC_FLAG_CTX                      , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if ALF_HM3_REFACTOR
, m_bAlfCtrl                           ( false )
, m_uiMaxAlfCtrlDepth                  ( 0 )
, m_cCUAlfCtrlFlagSCModel              ( 1,             1,               NUM_ALF_CTRL_FLAG_CTX         , m_contextModels + m_numContextModels, m_numContextModels)
#if !JVET_C0038_GALF
, m_cALFFlagSCModel                    ( 1,             1,               NUM_ALF_FLAG_CTX              , m_contextModels + m_numContextModels, m_numContextModels)
#endif
, m_cALFUvlcSCModel                    ( 1,             1,               NUM_ALF_UVLC_CTX              , m_contextModels + m_numContextModels, m_numContextModels)
#if !JVET_C0038_GALF
, m_cALFSvlcSCModel                    ( 1,             1,               NUM_ALF_SVLC_CTX              , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#endif
#if COM16_C806_EMT
, m_cEmtTuIdxSCModel                   ( 1,             1,               NUM_EMT_TU_IDX_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
, m_cEmtCuFlagSCModel                  ( 1,             1,               NUM_EMT_CU_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if COM16_C1016_AFFINE
, m_cCUAffineFlagSCModel               ( 1,             1,               NUM_AFFINE_FLAG_CTX           , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
, m_CUBgFlagSCModel                    ( 1,             1,               NUM_CU_BG_FLAG_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if MRG_CONSIDER_COMPLEX
, m_cCUMergeComplexFlagExtSCModel      ( 1,             1,               NUM_MERGE_COMPLEX_FLAG_EXT_CTX , m_contextModels + m_numContextModels, m_numContextModels)
, m_cCUMergeComplexIdxExtSCModel       ( 1,             1,               NUM_MERGE_COMPLEX_IDX_EXT_CTX  , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
, m_cCUAffineParamFlagSCModel          ( 1,             1,               NUM_AFFINE_FLAG_CTX           , m_contextModels + m_numContextModels,                   m_numContextModels)
#endif
#if RRC_NS
, m_cNsFlagSCModel                     ( 1,             1,               NUM_NS_FLAG_CTX               , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_INTRA_SHARPENING
, m_cSharpFlagSCModel                  ( 1,             1,               NUM_SHARP_FLAG_CTX            , m_contextModels + m_numContextModels, m_numContextModels)
#endif
#if RRC_SELECTABLE_IF
, m_cCUAdaptiveIfIndex                 ( 1,             1,               NUM_SELECTABLE_IF_CTX         , m_contextModels + m_numContextModels, m_numContextModels)
#endif
{
  assert( m_numContextModels <= MAX_NUM_CTX_MOD );
}

TEncSbac::~TEncSbac()
{
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

Void TEncSbac::resetEntropy           (const TComSlice *pSlice)
{
  Int  iQp              = pSlice->getSliceQp();
  SliceType eSliceType  = pSlice->getSliceType();

  SliceType encCABACTableIdx = pSlice->getEncCABACTableIdx();
  if (!pSlice->isIntra() && (encCABACTableIdx==B_SLICE || encCABACTableIdx==P_SLICE) && pSlice->getPPS()->getCabacInitPresentFlag())
  {
    eSliceType = encCABACTableIdx;
  }

  m_cCUSplitFlagSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_SPLIT_FLAG );
#if JVET_C0024_QTBT
  m_cBTSplitFlagSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_BTSPLIT_MODE );
#endif
  m_cCUSkipFlagSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_SKIP_FLAG );
#if BARC_UP_SAMPLE_RDO
  m_cCUUpCNNFlagSCModel.initBuffer                (eSliceType, iQp, (UChar*)INIT_UP_CNN_RDO_FLAG);
#endif
#if BARC
  m_cCUBARCFlagSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_BARC_FLAG);
#endif
#if VCEG_AZ05_INTRA_MPI
  m_cMPIIdxSCModel.initBuffer                     ( eSliceType, iQp, (UChar*)INIT_MPIIdx_FLAG );
#endif 
#if COM16_C1046_PDPC_INTRA
  m_cPDPCIdxSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_PDPCIdx_FLAG);
#endif
#if DIMD_INTRA_PRED
  m_cDIMDEnabledSCModel.initBuffer(eSliceType, iQp, (UChar*)INIT_DIMD_FLAG);
  m_cDIMDNoBTFlagSCModel.initBuffer(eSliceType, iQp, (UChar*)INIT_DIMD_NO_BT_FLAG);
#endif
#if RRC_DWDIP
  m_cDWDIPIdxSCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_DWDIPIdx_FLAG);
#endif
#if RRC_DWDIP_MLR
  m_cDWDIP_GPIdxSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_DWDIP_GPIdx_FLAG);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  m_cROTidxSCModel.initBuffer                     ( eSliceType, iQp, (UChar*)INIT_ROT_TR_IDX );
#endif
#if FAST_EMT_v1
  m_cEMTIntraidxSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_EMT_TR_IDX);
#endif
  m_cCUMergeFlagExtSCModel.initBuffer             ( eSliceType, iQp, (UChar*)INIT_MERGE_FLAG_EXT);
  m_cCUMergeIdxExtSCModel.initBuffer              ( eSliceType, iQp, (UChar*)INIT_MERGE_IDX_EXT);
#if VCEG_AZ07_FRUC_MERGE
  m_cCUFRUCMgrModeSCModel.initBuffer              ( eSliceType, iQp, (UChar*)INIT_FRUCMGRMODEBIN1 );
  m_cCUFRUCMESCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_FRUCMGRMODEBIN2 );
#if THIRD_FRUC_MODE
  m_cCUFRUCMESCModel_uni.initBuffer               ( eSliceType, iQp, (UChar*)INIT_FRUCMGRMODEBIN3 );
#endif
#endif
  m_cCUPartSizeSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_PART_SIZE );
  m_cCUPredModeSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_PRED_MODE );
  m_cCUIntraPredSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_INTRA_PRED_MODE );
#if MLIP_EXTEND_REF_LINE
  m_cCUExternRefSCModel.initBuffer(eSliceType, iQp, (UChar*)INIT_EXTERN_REF);
#endif
  m_cCUChromaPredSCModel.initBuffer               ( eSliceType, iQp, (UChar*)INIT_CHROMA_PRED_MODE );
  m_cCUInterDirSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_INTER_DIR );
  m_cCUMvdSCModel.initBuffer                      ( eSliceType, iQp, (UChar*)INIT_MVD );
  m_cCURefPicSCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_REF_PIC );
  m_cCUDeltaQpSCModel.initBuffer                  ( eSliceType, iQp, (UChar*)INIT_DQP );
  m_cCUQtCbfSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_QT_CBF );
  m_cCUQtRootCbfSCModel.initBuffer                ( eSliceType, iQp, (UChar*)INIT_QT_ROOT_CBF );
  m_cCUSigCoeffGroupSCModel.initBuffer            ( eSliceType, iQp, (UChar*)INIT_SIG_CG_FLAG );
  m_cCUSigSCModel.initBuffer                      ( eSliceType, iQp, (UChar*)INIT_SIG_FLAG );
  m_cCuCtxLastX.initBuffer                        ( eSliceType, iQp, (UChar*)INIT_LAST );
  m_cCuCtxLastY.initBuffer                        ( eSliceType, iQp, (UChar*)INIT_LAST );
  m_cCUOneSCModel.initBuffer                      ( eSliceType, iQp, (UChar*)INIT_ONE_FLAG );
#if !VCEG_AZ07_CTX_RESIDUALCODING
  m_cCUAbsSCModel.initBuffer                      ( eSliceType, iQp, (UChar*)INIT_ABS_FLAG );
#endif
  m_cMVPIdxSCModel.initBuffer                     ( eSliceType, iQp, (UChar*)INIT_MVP_IDX );
  m_cCUTransSubdivFlagSCModel.initBuffer          ( eSliceType, iQp, (UChar*)INIT_TRANS_SUBDIV_FLAG );
#if CNN_BASED_INLOOP_FILTER
	m_cCUCNNLoopFilterFlagSCModel.initBuffer        ( eSliceType, iQp, (UChar*)INIT_CNNLF_FLAG);
#endif
  m_cSaoMergeSCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_SAO_MERGE_FLAG );
  m_cSaoTypeIdxSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_SAO_TYPE_IDX );
  m_cTransformSkipSCModel.initBuffer              ( eSliceType, iQp, (UChar*)INIT_TRANSFORMSKIP_FLAG );
#if VCEG_AZ08_KLT_COMMON
  m_cKLTFlagSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_KLT_FLAG);
#endif
  m_CUTransquantBypassFlagSCModel.initBuffer      ( eSliceType, iQp, (UChar*)INIT_CU_TRANSQUANT_BYPASS_FLAG );
  m_explicitRdpcmFlagSCModel.initBuffer           ( eSliceType, iQp, (UChar*)INIT_EXPLICIT_RDPCM_FLAG);
  m_explicitRdpcmDirSCModel.initBuffer            ( eSliceType, iQp, (UChar*)INIT_EXPLICIT_RDPCM_DIR);
  m_cCrossComponentPredictionSCModel.initBuffer   ( eSliceType, iQp, (UChar*)INIT_CROSS_COMPONENT_PREDICTION  );
  m_ChromaQpAdjFlagSCModel.initBuffer             ( eSliceType, iQp, (UChar*)INIT_CHROMA_QP_ADJ_FLAG );
  m_ChromaQpAdjIdcSCModel.initBuffer              ( eSliceType, iQp, (UChar*)INIT_CHROMA_QP_ADJ_IDC );
#if COM16_C806_OBMC
  m_cCUOBMCFlagSCModel.initBuffer                 ( eSliceType, iQp, (UChar*)INIT_OBMC_FLAG );
#endif
#if VCEG_AZ07_IMV
  m_cCUiMVFlagSCModel.initBuffer                  ( eSliceType, iQp, (UChar*)INIT_IMV_FLAG );
#endif
#if VCEG_AZ06_IC
  m_cCUICFlagSCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_IC_FLAG );
#endif
#if ALF_HM3_REFACTOR
  m_cCUAlfCtrlFlagSCModel.initBuffer              ( eSliceType, iQp, (UChar*)INIT_ALF_CTRL_FLAG );
#if !JVET_C0038_GALF
  m_cALFFlagSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_ALF_FLAG );
#endif
  m_cALFUvlcSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_ALF_UVLC );
#if !JVET_C0038_GALF
  m_cALFSvlcSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_ALF_SVLC );
#endif
#endif
#if COM16_C806_EMT
  m_cEmtTuIdxSCModel.initBuffer                   ( eSliceType, iQp, (UChar*)INIT_EMT_TU_IDX );
  m_cEmtCuFlagSCModel.initBuffer                  ( eSliceType, iQp, (UChar*)INIT_EMT_CU_FLAG );
#endif 
#if COM16_C1016_AFFINE
  m_cCUAffineFlagSCModel.initBuffer               ( eSliceType, iQp, (UChar*)INIT_AFFINE_FLAG );
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
	m_CUBgFlagSCModel.initBuffer                    ( eSliceType, iQp, (UChar*)INIT_CU_BG_FLAG );
#endif
#if MRG_CONSIDER_COMPLEX
  m_cCUMergeComplexFlagExtSCModel.initBuffer      ( eSliceType, iQp, (UChar*)INIT_MERGE_COMPLEX_FLAG_EXT );
  m_cCUMergeComplexIdxExtSCModel.initBuffer       ( eSliceType, iQp, (UChar*)INIT_MERGE_COMPLEX_IDX_EXT );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
  m_cCUAffineParamFlagSCModel.initBuffer          ( eSliceType, iQp, (UChar*)INIT_AFFINE_FLAG );
#endif
#if RRC_NS
  m_cNsFlagSCModel.initBuffer                   (eSliceType, iQp, (UChar*) INIT_NS_FLAG);
#endif
#if RRC_INTRA_SHARPENING
  m_cSharpFlagSCModel.initBuffer                  ( eSliceType, iQp, (UChar*) INIT_SHARP_FLAG );
#endif
#if RRC_SELECTABLE_IF
  m_cCUAdaptiveIfIndex.initBuffer                 (eSliceType, iQp, (UChar*) INIT_ADAPTIVE_IF_INDEX);
#endif

  for (UInt statisticIndex = 0; statisticIndex < RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS ; statisticIndex++)
  {
    m_golombRiceAdaptationStatistics[statisticIndex] = 0;
  }

  m_pcBinIf->start();
#if VCEG_AZ07_BAC_ADAPT_WDOW 
  xUpdateWindowSize (pSlice->getSliceType(), pSlice->getCtxMapQPIdx(), pSlice->getStatsHandle());
#endif
  return;
}

/** The function does the following:
 * If current slice type is P/B then it determines the distance of initialisation type 1 and 2 from the current CABAC states and
 * stores the index of the closest table.  This index is used for the next P/B slice when cabac_init_present_flag is true.
 */
SliceType TEncSbac::determineCabacInitIdx(const TComSlice *pSlice)
{
  Int  qp              = pSlice->getSliceQp();

  if (!pSlice->isIntra())
  {
    SliceType aSliceTypeChoices[] = {B_SLICE, P_SLICE};

    UInt bestCost             = MAX_UINT;
    SliceType bestSliceType   = aSliceTypeChoices[0];
    for (UInt idx=0; idx<2; idx++)
    {
      UInt curCost          = 0;
      SliceType curSliceType  = aSliceTypeChoices[idx];

      curCost  = m_cCUSplitFlagSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_SPLIT_FLAG );
#if JVET_C0024_QTBT
      curCost += m_cBTSplitFlagSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_BTSPLIT_MODE );
#endif
      curCost += m_cCUSkipFlagSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_SKIP_FLAG );
#if BARC_UP_SAMPLE_RDO
	  curCost += m_cCUUpCNNFlagSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_UP_CNN_RDO_FLAG);
#endif
#if BARC
      curCost += m_cCUBARCFlagSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_BARC_FLAG);
#endif
#if VCEG_AZ05_INTRA_MPI
      curCost += m_cMPIIdxSCModel.calcCost                     ( curSliceType, qp, (UChar*)INIT_MPIIdx_FLAG );
#endif
#if COM16_C1046_PDPC_INTRA
      curCost += m_cPDPCIdxSCModel.calcCost                    ( curSliceType, qp, (UChar*)INIT_PDPCIdx_FLAG);
#endif
#if DIMD_INTRA_PRED
      curCost += m_cDIMDEnabledSCModel.calcCost(curSliceType, qp, (UChar*)INIT_DIMD_FLAG);
      curCost += m_cDIMDNoBTFlagSCModel.calcCost(curSliceType, qp, (UChar*)INIT_DIMD_NO_BT_FLAG);
#endif
#if RRC_DWDIP  
      curCost += m_cDWDIPIdxSCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_DWDIPIdx_FLAG);
#endif
#if RRC_DWDIP_MLR
      curCost += m_cDWDIP_GPIdxSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_DWDIP_GPIdx_FLAG);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
      curCost += m_cROTidxSCModel.calcCost        ( curSliceType, qp, (UChar*)INIT_ROT_TR_IDX );
#endif
#if FAST_EMT_v1
	  curCost += m_cEMTIntraidxSCModel.calcCost    (curSliceType, qp, (UChar*)INIT_EMT_TR_IDX);
#endif
      curCost += m_cCUMergeFlagExtSCModel.calcCost             ( curSliceType, qp, (UChar*)INIT_MERGE_FLAG_EXT);
      curCost += m_cCUMergeIdxExtSCModel.calcCost              ( curSliceType, qp, (UChar*)INIT_MERGE_IDX_EXT);
#if VCEG_AZ07_FRUC_MERGE
      curCost += m_cCUFRUCMgrModeSCModel.calcCost              ( curSliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN1);
      curCost += m_cCUFRUCMESCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN2);
#if THIRD_FRUC_MODE
      curCost += m_cCUFRUCMESCModel_uni.calcCost               ( curSliceType, qp, (UChar*)INIT_FRUCMGRMODEBIN3);
#endif
#endif
      curCost += m_cCUPartSizeSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_PART_SIZE );
      curCost += m_cCUPredModeSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_PRED_MODE );
      curCost += m_cCUIntraPredSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_INTRA_PRED_MODE );
#if MLIP_EXTEND_REF_LINE
      curCost += m_cCUExternRefSCModel.calcCost(curSliceType, qp, (UChar*)INIT_EXTERN_REF);
#endif
      curCost += m_cCUChromaPredSCModel.calcCost               ( curSliceType, qp, (UChar*)INIT_CHROMA_PRED_MODE );
      curCost += m_cCUInterDirSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_INTER_DIR );
      curCost += m_cCUMvdSCModel.calcCost                      ( curSliceType, qp, (UChar*)INIT_MVD );
      curCost += m_cCURefPicSCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_REF_PIC );
      curCost += m_cCUDeltaQpSCModel.calcCost                  ( curSliceType, qp, (UChar*)INIT_DQP );
      curCost += m_cCUQtCbfSCModel.calcCost                    ( curSliceType, qp, (UChar*)INIT_QT_CBF );
      curCost += m_cCUQtRootCbfSCModel.calcCost                ( curSliceType, qp, (UChar*)INIT_QT_ROOT_CBF );
      curCost += m_cCUSigCoeffGroupSCModel.calcCost            ( curSliceType, qp, (UChar*)INIT_SIG_CG_FLAG );
      curCost += m_cCUSigSCModel.calcCost                      ( curSliceType, qp, (UChar*)INIT_SIG_FLAG );
      curCost += m_cCuCtxLastX.calcCost                        ( curSliceType, qp, (UChar*)INIT_LAST );
      curCost += m_cCuCtxLastY.calcCost                        ( curSliceType, qp, (UChar*)INIT_LAST );
      curCost += m_cCUOneSCModel.calcCost                      ( curSliceType, qp, (UChar*)INIT_ONE_FLAG );
#if !VCEG_AZ07_CTX_RESIDUALCODING
      curCost += m_cCUAbsSCModel.calcCost                      ( curSliceType, qp, (UChar*)INIT_ABS_FLAG );
#endif
      curCost += m_cMVPIdxSCModel.calcCost                     ( curSliceType, qp, (UChar*)INIT_MVP_IDX );
      curCost += m_cCUTransSubdivFlagSCModel.calcCost          ( curSliceType, qp, (UChar*)INIT_TRANS_SUBDIV_FLAG );
#if CNN_BASED_INLOOP_FILTER
      curCost += m_cCUCNNLoopFilterFlagSCModel.calcCost        ( curSliceType, qp, (UChar*)INIT_CNNLF_FLAG);
#endif
      curCost += m_cSaoMergeSCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_SAO_MERGE_FLAG );
      curCost += m_cSaoTypeIdxSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_SAO_TYPE_IDX );
      curCost += m_cTransformSkipSCModel.calcCost              ( curSliceType, qp, (UChar*)INIT_TRANSFORMSKIP_FLAG );
#if VCEG_AZ08_KLT_COMMON
      curCost += m_cKLTFlagSCModel.calcCost                    ( curSliceType, qp, (UChar*)INIT_KLT_FLAG);
#endif
      curCost += m_CUTransquantBypassFlagSCModel.calcCost      ( curSliceType, qp, (UChar*)INIT_CU_TRANSQUANT_BYPASS_FLAG );
      curCost += m_explicitRdpcmFlagSCModel.calcCost           ( curSliceType, qp, (UChar*)INIT_EXPLICIT_RDPCM_FLAG);
      curCost += m_explicitRdpcmDirSCModel.calcCost            ( curSliceType, qp, (UChar*)INIT_EXPLICIT_RDPCM_DIR);
      curCost += m_cCrossComponentPredictionSCModel.calcCost   ( curSliceType, qp, (UChar*)INIT_CROSS_COMPONENT_PREDICTION );
      curCost += m_ChromaQpAdjFlagSCModel.calcCost             ( curSliceType, qp, (UChar*)INIT_CHROMA_QP_ADJ_FLAG );
      curCost += m_ChromaQpAdjIdcSCModel.calcCost              ( curSliceType, qp, (UChar*)INIT_CHROMA_QP_ADJ_IDC );
#if COM16_C806_OBMC
      curCost += m_cCUOBMCFlagSCModel.calcCost                 ( curSliceType, qp, (UChar*)INIT_OBMC_FLAG );
#endif
#if COM16_C806_EMT
      curCost += m_cEmtTuIdxSCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_EMT_TU_IDX );
      curCost += m_cEmtCuFlagSCModel.calcCost                  ( curSliceType, qp, (UChar*)INIT_EMT_CU_FLAG );
#endif
#if VCEG_AZ07_IMV
      curCost += m_cCUiMVFlagSCModel.calcCost                  ( curSliceType, qp, (UChar*)INIT_IMV_FLAG );
#endif
#if VCEG_AZ06_IC
      curCost += m_cCUICFlagSCModel.calcCost                   ( curSliceType, qp, (UChar*)INIT_IC_FLAG );
#endif
#if COM16_C1016_AFFINE
      curCost += m_cCUAffineFlagSCModel.calcCost               ( curSliceType, qp, (UChar*)INIT_AFFINE_FLAG );
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
			curCost += m_CUBgFlagSCModel.calcCost                    ( curSliceType, qp, (UChar*)INIT_CU_BG_FLAG );
#endif
#if MRG_CONSIDER_COMPLEX
      curCost += m_cCUMergeComplexFlagExtSCModel.calcCost      ( curSliceType, qp, (UChar*)INIT_MERGE_COMPLEX_FLAG_EXT);
      curCost += m_cCUMergeComplexIdxExtSCModel.calcCost       ( curSliceType, qp, (UChar*)INIT_MERGE_COMPLEX_IDX_EXT);
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
      curCost += m_cCUAffineParamFlagSCModel.calcCost          ( curSliceType, qp, (UChar*)INIT_AFFINE_FLAG );
#endif
#if RRC_NS
      curCost += m_cNsFlagSCModel.calcCost                   (curSliceType, qp, (UChar*) INIT_NS_FLAG);
#endif
#if RRC_INTRA_SHARPENING
      curCost += m_cSharpFlagSCModel.calcCost                  ( curSliceType, qp, (UChar*)INIT_SHARP_FLAG );
#endif
#if RRC_SELECTABLE_IF
      curCost += m_cCUAdaptiveIfIndex.calcCost                 (curSliceType, qp, (UChar*) INIT_ADAPTIVE_IF_INDEX);
#endif

      if (curCost < bestCost)
      {
        bestSliceType = curSliceType;
        bestCost      = curCost;
      }
    }
    return bestSliceType;
  }
  else
  {
    return I_SLICE;
  }
}

Void TEncSbac::codeVPS( const TComVPS* /*pcVPS*/ )
{
  assert (0);
  return;
}

Void TEncSbac::codeSPS( const TComSPS* /*pcSPS*/ )
{
  assert (0);
  return;
}

Void TEncSbac::codePPS( const TComPPS* /*pcPPS*/ )
{
  assert (0);
  return;
}

Void TEncSbac::codeSliceHeader( TComSlice* /*pcSlice*/ )
{
  assert (0);
  return;
}

Void TEncSbac::codeTilesWPPEntryPoint( TComSlice* /*pSlice*/ )
{
  assert (0);
  return;
}

Void TEncSbac::codeTerminatingBit( UInt uilsLast )
{
  m_pcBinIf->encodeBinTrm( uilsLast );
}

Void TEncSbac::codeSliceFinish()
{
  m_pcBinIf->finish();
}

Void TEncSbac::xWriteUnarySymbol( UInt uiSymbol, ContextModel* pcSCModel, Int iOffset )
{
  m_pcBinIf->encodeBin( uiSymbol ? 1 : 0, pcSCModel[0] );

  if( 0 == uiSymbol)
  {
    return;
  }

  while( uiSymbol-- )
  {
    m_pcBinIf->encodeBin( uiSymbol ? 1 : 0, pcSCModel[ iOffset ] );
  }

  return;
}

Void TEncSbac::xWriteUnaryMaxSymbol( UInt uiSymbol, ContextModel* pcSCModel, Int iOffset, UInt uiMaxSymbol )
{
  if (uiMaxSymbol == 0)
  {
    return;
  }

  m_pcBinIf->encodeBin( uiSymbol ? 1 : 0, pcSCModel[ 0 ] );

  if ( uiSymbol == 0 )
  {
    return;
  }

  Bool bCodeLast = ( uiMaxSymbol > uiSymbol );

  while( --uiSymbol )
  {
    m_pcBinIf->encodeBin( 1, pcSCModel[ iOffset ] );
  }
  if( bCodeLast )
  {
    m_pcBinIf->encodeBin( 0, pcSCModel[ iOffset ] );
  }

  return;
}

#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
Void TEncSbac::xWriteTruncBinCode(UInt uiSymbol, UInt uiMaxSymbol)
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
  assert(uiVal <= uiMaxSymbol);
  assert((uiVal << 1) > uiMaxSymbol);
  assert(uiSymbol < uiMaxSymbol);
  UInt b = uiMaxSymbol - uiVal;
  assert(b < uiVal);
  if (uiSymbol < uiVal - b)
  {
    m_pcBinIf->encodeBinsEP(uiSymbol, uiThresh);
  }
  else
  {
    uiSymbol += uiVal - b;
    assert(uiSymbol < (uiVal << 1));
    assert((uiSymbol >> 1) >= uiVal - b);
    m_pcBinIf->encodeBinsEP(uiSymbol, uiThresh + 1);
  }
}
#endif

Void TEncSbac::xWriteEpExGolomb( UInt uiSymbol, UInt uiCount )
{
  UInt bins = 0;
  Int numBins = 0;

  while( uiSymbol >= (UInt)(1<<uiCount) )
  {
    bins = 2 * bins + 1;
    numBins++;
    uiSymbol -= 1 << uiCount;
    uiCount  ++;
  }
  bins = 2 * bins + 0;
  numBins++;

  bins = (bins << uiCount) | uiSymbol;
  numBins += uiCount;

  assert( numBins <= 32 );
  m_pcBinIf->encodeBinsEP( bins, numBins );
}


/** Coding of coeff_abs_level_minus3
 * \param symbol                  value of coeff_abs_level_minus3
 * \param rParam                  reference to Rice parameter
 * \param useLimitedPrefixLength
 * \param maxLog2TrDynamicRange 
 */

Void TEncSbac::xWriteCoefRemainExGolomb ( UInt symbol, UInt &rParam, const Bool useLimitedPrefixLength, const Int maxLog2TrDynamicRange )
{
  Int codeNumber  = (Int)symbol;
  UInt length;
#if VCEG_AZ07_CTX_RESIDUALCODING
  if (codeNumber < (g_auiGoRiceRange[rParam] << rParam))
#else
  if (codeNumber < (COEF_REMAIN_BIN_REDUCTION << rParam))
#endif
  {
    length = codeNumber>>rParam;
    m_pcBinIf->encodeBinsEP( (1<<(length+1))-2 , length+1);
    m_pcBinIf->encodeBinsEP((codeNumber%(1<<rParam)),rParam);
  }
  else if (useLimitedPrefixLength)
  {
    const UInt maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange));

    UInt prefixLength = 0;
    UInt suffixLength = MAX_UINT;
    UInt codeValue    = (symbol >> rParam) - COEF_REMAIN_BIN_REDUCTION;

    if (codeValue >= ((1 << maximumPrefixLength) - 1))
    {
      prefixLength = maximumPrefixLength;
      suffixLength = maxLog2TrDynamicRange - rParam;
    }
    else
    {
      while (codeValue > ((2 << prefixLength) - 2))
      {
        prefixLength++;
      }

      suffixLength = prefixLength + 1; //+1 for the separator bit
    }

    const UInt suffix = codeValue - ((1 << prefixLength) - 1);

    const UInt totalPrefixLength = prefixLength + COEF_REMAIN_BIN_REDUCTION;
    const UInt prefix            = (1 << totalPrefixLength) - 1;
    const UInt rParamBitMask     = (1 << rParam) - 1;

    m_pcBinIf->encodeBinsEP(  prefix,                                        totalPrefixLength      ); //prefix
    m_pcBinIf->encodeBinsEP(((suffix << rParam) | (symbol & rParamBitMask)), (suffixLength + rParam)); //separator, suffix, and rParam bits
  }
  else
  {
    length = rParam;
#if VCEG_AZ07_CTX_RESIDUALCODING
    codeNumber  = codeNumber - ( g_auiGoRiceRange[rParam] << rParam);
#else
    codeNumber  = codeNumber - ( COEF_REMAIN_BIN_REDUCTION << rParam);
#endif

    while (codeNumber >= (1<<length))
    {
      codeNumber -=  (1<<(length++));
    }
#if VCEG_AZ07_CTX_RESIDUALCODING
    m_pcBinIf->encodeBinsEP((1<<(g_auiGoRiceRange[rParam] + length + 1 - rParam))-2, g_auiGoRiceRange[rParam] + length + 1 - rParam);
#else
    m_pcBinIf->encodeBinsEP((1<<(COEF_REMAIN_BIN_REDUCTION+length+1-rParam))-2,COEF_REMAIN_BIN_REDUCTION+length+1-rParam);
#endif
    m_pcBinIf->encodeBinsEP(codeNumber,length);
  }
}

// SBAC RD
Void  TEncSbac::load ( const TEncSbac* pSrc)
{
  this->xCopyFrom(pSrc);
}

Void  TEncSbac::loadIntraDirMode( const TEncSbac* pSrc, const ChannelType chType )
{
  m_pcBinIf->copyState( pSrc->m_pcBinIf );
  if (isLuma(chType))
  {
    this->m_cCUIntraPredSCModel      .copyFrom( &pSrc->m_cCUIntraPredSCModel       );
#if RRC_DWDIP
    this->m_cDWDIPIdxSCModel         .copyFrom( &pSrc->m_cDWDIPIdxSCModel          );
#endif
#if RRC_DWDIP_MLR
    this->m_cDWDIP_GPIdxSCModel      .copyFrom( &pSrc->m_cDWDIP_GPIdxSCModel       );
#endif
  }
  else
  {
    this->m_cCUChromaPredSCModel     .copyFrom( &pSrc->m_cCUChromaPredSCModel      );
  }
}


Void  TEncSbac::store( TEncSbac* pDest) const
{
  pDest->xCopyFrom( this );
}


Void TEncSbac::xCopyFrom( const TEncSbac* pSrc )
{
  m_pcBinIf->copyState( pSrc->m_pcBinIf );
  xCopyContextsFrom(pSrc);
}

Void TEncSbac::codeMVPIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  Int iSymbol = pcCU->getMVPIdx(eRefList, uiAbsPartIdx);
  Int iNum = AMVP_MAX_NUM_CANDS;

  xWriteUnaryMaxSymbol(iSymbol, m_cMVPIdxSCModel.get(0), 1, iNum-1);
}

#if !JVET_C0024_QTBT
Void TEncSbac::codePartSize( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  PartSize eSize         = pcCU->getPartitionSize( uiAbsPartIdx );
  const UInt log2DiffMaxMinCodingBlockSize = pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize();

  if ( pcCU->isIntra( uiAbsPartIdx ) )
  {
    if( uiDepth == log2DiffMaxMinCodingBlockSize )
    {
      m_pcBinIf->encodeBin( eSize == SIZE_2Nx2N? 1 : 0, m_cCUPartSizeSCModel.get( 0, 0, 0 ) );
    }
    return;
  }

  switch(eSize)
  {
    case SIZE_2Nx2N:
    {
      m_pcBinIf->encodeBin( 1, m_cCUPartSizeSCModel.get( 0, 0, 0) );
      break;
    }
    case SIZE_2NxN:
    case SIZE_2NxnU:
    case SIZE_2NxnD:
    {
      m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 0) );
      m_pcBinIf->encodeBin( 1, m_cCUPartSizeSCModel.get( 0, 0, 1) );
      if ( pcCU->getSlice()->getSPS()->getUseAMP() && uiDepth < log2DiffMaxMinCodingBlockSize )
      {
        if (eSize == SIZE_2NxN)
        {
          m_pcBinIf->encodeBin(1, m_cCUPartSizeSCModel.get( 0, 0, 3 ));
        }
        else
        {
          m_pcBinIf->encodeBin(0, m_cCUPartSizeSCModel.get( 0, 0, 3 ));
          m_pcBinIf->encodeBinEP((eSize == SIZE_2NxnU? 0: 1));
        }
      }
      break;
    }
    case SIZE_Nx2N:
    case SIZE_nLx2N:
    case SIZE_nRx2N:
    {
      m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 0) );
      m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 1) );
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL && !COM16_C806_DISABLE_4X4_PU
      if( uiDepth == log2DiffMaxMinCodingBlockSize && (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag()||!( pcCU->getWidth(uiAbsPartIdx) == 8 && pcCU->getHeight(uiAbsPartIdx) == 8 )))
#else
      if( uiDepth == log2DiffMaxMinCodingBlockSize && !( pcCU->getWidth(uiAbsPartIdx) == 8 && pcCU->getHeight(uiAbsPartIdx) == 8 ) )
#endif
      {
        m_pcBinIf->encodeBin( 1, m_cCUPartSizeSCModel.get( 0, 0, 2) );
      }

      if ( pcCU->getSlice()->getSPS()->getUseAMP() && uiDepth < log2DiffMaxMinCodingBlockSize )
      {
        if (eSize == SIZE_Nx2N)
        {
          m_pcBinIf->encodeBin(1, m_cCUPartSizeSCModel.get( 0, 0, 3 ));
        }
        else
        {
          m_pcBinIf->encodeBin(0, m_cCUPartSizeSCModel.get( 0, 0, 3 ));
          m_pcBinIf->encodeBinEP((eSize == SIZE_nLx2N? 0: 1));
        }
      }
      break;
    }
    case SIZE_NxN:
    {
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL && !COM16_C806_DISABLE_4X4_PU
      if( uiDepth == log2DiffMaxMinCodingBlockSize && 
        (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag()||!( pcCU->getWidth(uiAbsPartIdx) == 8 && pcCU->getHeight(uiAbsPartIdx) == 8 )))
#else
      if( uiDepth == log2DiffMaxMinCodingBlockSize && !( pcCU->getWidth(uiAbsPartIdx) == 8 && pcCU->getHeight(uiAbsPartIdx) == 8 ) )
#endif
      {
        m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 0) );
        m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 1) );
        m_pcBinIf->encodeBin( 0, m_cCUPartSizeSCModel.get( 0, 0, 2) );
      }
      break;
    }
    default:
    {
      assert(0);
      break;
    }
  }
}
#endif


/** code prediction mode
 * \param pcCU
 * \param uiAbsPartIdx
 * \returns Void
 */
Void TEncSbac::codePredMode( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  // get context function is here
  m_pcBinIf->encodeBin( pcCU->isIntra( uiAbsPartIdx ) ? 1 : 0, m_cCUPredModeSCModel.get( 0, 0, 0 ) );
}

Void TEncSbac::codeCUTransquantBypassFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiSymbol = pcCU->getCUTransquantBypass(uiAbsPartIdx);
  m_pcBinIf->encodeBin( uiSymbol, m_CUTransquantBypassFlagSCModel.get( 0, 0, 0 ) );
}

#if BLOCK_BASED_BACKGROUND_REFERENCE
Void TEncSbac::codeLCUBgFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
	UInt uiSymbol = pcCU->getLCUBgFlag() ? 1 : 0;
	m_pcBinIf->encodeBin(uiSymbol, m_CUBgFlagSCModel.get(0, 0, 0));
}
#endif

/** code skip flag
 * \param pcCU
 * \param uiAbsPartIdx
 * \returns Void
 */
Void TEncSbac::codeSkipFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  // get context function is here
  UInt uiSymbol = pcCU->isSkipped( uiAbsPartIdx ) ? 1 : 0;
  UInt uiCtxSkip = pcCU->getCtxSkipFlag( uiAbsPartIdx ) ;
  m_pcBinIf->encodeBin( uiSymbol, m_cCUSkipFlagSCModel.get( 0, 0, uiCtxSkip ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tSkipFlag" );
  DTRACE_CABAC_T( "\tuiCtxSkip: ");
  DTRACE_CABAC_V( uiCtxSkip );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}

#if COM16_C806_OBMC
Void TEncSbac::codeOBMCFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiSymbol = pcCU->getOBMCFlag( uiAbsPartIdx ) ? 1 : 0;

  m_pcBinIf->encodeBin( uiSymbol, m_cCUOBMCFlagSCModel.get( 0, 0, 0 ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tOBMCFlag" );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}
#endif

#if VCEG_AZ06_IC
/** code Illumination Compensation flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \returns Void
 */
Void TEncSbac::codeICFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiSymbol = pcCU->getICFlag( uiAbsPartIdx ) ? 1 : 0;
  m_pcBinIf->encodeBin( uiSymbol, m_cCUICFlagSCModel.get( 0, 0, 0 ) );

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tICFlag" );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}
#endif

#if VCEG_AZ07_IMV
Void TEncSbac::codeiMVFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  // get context function is here
  UInt uiSymbol = pcCU->getiMVFlag( uiAbsPartIdx ) ? 1 : 0;
  UInt uiCtxiMV = pcCU->getCtxiMVFlag( uiAbsPartIdx ) ;
  m_pcBinIf->encodeBin( uiSymbol, m_cCUiMVFlagSCModel.get( 0, 0, uiCtxiMV ) );

#if  JVET_E0076_MULTI_PEL_MVD
  if (uiSymbol)
  {
    uiSymbol = pcCU->getiMVFlag( uiAbsPartIdx ) > 1 ? 1 : 0;  
    m_pcBinIf->encodeBin( uiSymbol, m_cCUiMVFlagSCModel.get( 0, 0, 3 ) );
  }
#endif

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tiMVFlag" );
  DTRACE_CABAC_T( "\tuiCtxiMV: ");
  DTRACE_CABAC_V( uiCtxiMV );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}
#endif


#if VCEG_AZ05_INTRA_MPI
Void TEncSbac::codeMPIIdx(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  if (!pcCU->getSlice()->getSPS()->getUseMPI()) return;
  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
    return;
  }
  Int iNumberOfPassesMPI = 1;
  if (pcCU->getSlice()->getSliceType() == I_SLICE) iNumberOfPassesMPI = 2;
  else iNumberOfPassesMPI = 2;
  if (iNumberOfPassesMPI>1) // for only 1 pass no signaling is needed 
  {
    if (iNumberOfPassesMPI>2)  // 3 or 4
    {
      Int idxMPI = pcCU->getMPIIdx(uiAbsPartIdx);
      const UInt uiSymbol0 = (idxMPI >> 1);
      const UInt uiSymbol1 = (idxMPI % 2);
      m_pcBinIf->encodeBin(uiSymbol0, m_cMPIIdxSCModel.get(0, 0, 0));
      m_pcBinIf->encodeBin(uiSymbol1, m_cMPIIdxSCModel.get(0, 0, 1));
    }
    else //iNumberOfPassesMPI==2
    {
      Int idxMPI = pcCU->getMPIIdx(uiAbsPartIdx);
      const UInt uiSymbol = idxMPI;
      m_pcBinIf->encodeBin(uiSymbol, m_cMPIIdxSCModel.get(0, 0, 0));
    }
  }
}
#endif

#if DIMD_INTRA_PRED
Void TEncSbac::codeDIMDFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  UChar ucDIMDFlag = pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
  m_pcBinIf->encodeBin(ucDIMDFlag, m_cDIMDEnabledSCModel.get(0, 0, uiDepth));
}

Void TEncSbac::codeDIMDNoBTFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  UChar ucDIMDNoBTFlag = pcCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
  m_pcBinIf->encodeBin(ucDIMDNoBTFlag, m_cDIMDNoBTFlagSCModel.get(0, 0, uiDepth));
}
#endif

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
 Void TEncSbac::codePDPCIdx(TComDataCU* pcCU, UInt uiAbsPartIdx)
 {

  if (!pcCU->getSlice()->getSPS()->getUsePDPC()) return;
  if (pcCU->getPredictionMode(uiAbsPartIdx) == MODE_INTER)
  {
    return;
  }
  Int iNumberOfPassesPDPC = 1;
  if (pcCU->getSlice()->getSliceType() == I_SLICE) iNumberOfPassesPDPC = 2;
  else iNumberOfPassesPDPC = 2;
  if (iNumberOfPassesPDPC > 1) // for only 1 pass no signaling is needed 
  {
    if (iNumberOfPassesPDPC > 2)  // 3 or 4
    {
      Int idxPDPC = pcCU->getPDPCIdx(uiAbsPartIdx);
      const UInt uiSymbol0 = (idxPDPC >> 1);
      const UInt uiSymbol1 = (idxPDPC % 2);
      m_pcBinIf->encodeBin(uiSymbol0, m_cPDPCIdxSCModel.get(0, 0, 0));
      m_pcBinIf->encodeBin(uiSymbol1, m_cPDPCIdxSCModel.get(0, 0, 1));
    }
    else //iNumberOfPassesMPI==2
    {
      Int idxPDPC = pcCU->getPDPCIdx(uiAbsPartIdx);
      const UInt uiSymbol = idxPDPC;
      m_pcBinIf->encodeBin(uiSymbol, m_cPDPCIdxSCModel.get(0, 0, 0));
    }
   }
 }
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
Void TEncSbac::codeROTIdx ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth  )
{


#if COM16_C1044_NSST
  if (!pcCU->getSlice()->getSPS()->getUseNSST()) return;
#else
  if (!pcCU->getSlice()->getSPS()->getUseROT()) return;
#endif
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
  if( iNumberOfPassesROT==4 )
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
	  Int idxROT = pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
	  m_pcBinIf->encodeBin(idxROT ? 1 : 0, m_cROTidxSCModel.get(0, 0, 1));
	  return;
  }
#endif
  if( iNumberOfPassesROT==3 )
  {
#if JVET_C0024_QTBT
    Int idxROT = pcCU->getROTIdx( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
#else
    Int idxROT = pcCU->getROTIdx( uiAbsPartIdx );
#endif

#if RRC_DWDIP_MLR
    if (pcCU->getDWDIP_GPIdx(uiAbsPartIdx)==0)
#endif
    assert(idxROT<3);

#if JVET_C0042_UNIFIED_BINARIZATION
            m_pcBinIf->encodeBin(  idxROT ? 1 : 0 , m_cROTidxSCModel.get(0,0, 1) );
            if( idxROT )
            {
                   if(idxROT ==1) m_pcBinIf->encodeBin( 0 , m_cROTidxSCModel.get(0,0, 3) );
                   else m_pcBinIf->encodeBin( 1 , m_cROTidxSCModel.get(0,0, 3) );
            } 
#else //#if JVET_C0042_UNIFIED_BINARIZATION

    m_pcBinIf->encodeBin( idxROT ? 1 : 0, m_cROTidxSCModel.get(0,0, 0 ) );
    if( idxROT )
    {
      m_pcBinIf->encodeBin( (idxROT-1) ? 1 : 0, m_cROTidxSCModel.get(0,0, 1 ) );
    }
 #endif //#if JVET_C0042_UNIFIED_BINARIZATION
  }
  else
#endif
  if (iNumberOfPassesROT>1) // for only 1 pass no signaling is needed 
  {
#if JVET_C0024_QTBT
      Int idxROT = pcCU->getROTIdx( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
#else
      Int idxROT = pcCU->getROTIdx( uiAbsPartIdx );
#endif

#if JVET_C0042_UNIFIED_BINARIZATION
            m_pcBinIf->encodeBin(  idxROT ? 1 : 0 , m_cROTidxSCModel.get(0,0, 0) );
            if( idxROT )
            {
                     m_pcBinIf->encodeBin( (idxROT-1) ? 1 : 0 , m_cROTidxSCModel.get(0,0, 2) );

                    if(idxROT >1 )
                    {
                        m_pcBinIf->encodeBin( (idxROT-2) ? 1 : 0, m_cROTidxSCModel.get(0,0, 4) );
                    }
               
            } 
#else
      const UInt uiSymbol0 = (idxROT >>1);
      const UInt uiSymbol1 = (idxROT %2 );
      m_pcBinIf->encodeBin( uiSymbol0, m_cROTidxSCModel.get(0,0, uiDepth ) );
      m_pcBinIf->encodeBin( uiSymbol1, m_cROTidxSCModel.get(0,0, uiDepth ) );
#endif
  }
}

#if JVET_C0024_QTBT
Void TEncSbac::codeROTIdxChroma ( TComDataCU* pcCU, UInt uiAbsPartIdx,UInt uiDepth  )
{
#if DIMD_INTRA_PRED
  assert(!pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, uiAbsPartIdx));
#endif

#if COM16_C1044_NSST
  if (!pcCU->getSlice()->getSPS()->getUseNSST()) return;
#else
  if (!pcCU->getSlice()->getSPS()->getUseROT()) return;
#endif
  Int iNumberOfPassesROT = 1;
  if( pcCU->isIntra(uiAbsPartIdx)
    && !pcCU->getCUTransquantBypass(uiAbsPartIdx)
    )  iNumberOfPassesROT = 4;

#if COM16_C1044_NSST
  if( iNumberOfPassesROT==4 )
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
      uiIntraMode = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu()+uiAbsPartIdx);
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
	  Int idxROT = pcCU->getROTIdx(CHANNEL_TYPE_CHROMA, uiAbsPartIdx);
	  m_pcBinIf->encodeBin(idxROT ? 1 : 0, m_cROTidxSCModel.get(0, 0, 1));
	  return;
  }
#endif
  if( iNumberOfPassesROT==3 )
  {
    Int idxROT = pcCU->getROTIdx( CHANNEL_TYPE_CHROMA, uiAbsPartIdx );
    assert(idxROT<3);
#if JVET_C0042_UNIFIED_BINARIZATION
      m_pcBinIf->encodeBin(  idxROT ? 1 : 0 , m_cROTidxSCModel.get(0,0, 1) );
            if( idxROT )
            {
                   if(idxROT ==1) m_pcBinIf->encodeBin( 0 , m_cROTidxSCModel.get(0,0, 3) );
                   else m_pcBinIf->encodeBin( 1 , m_cROTidxSCModel.get(0,0, 3) );
            } 
#else
    m_pcBinIf->encodeBin( idxROT ? 1 : 0, m_cROTidxSCModel.get(0,0, 0 ) );
    if( idxROT )
    {
      m_pcBinIf->encodeBin( (idxROT-1) ? 1 : 0, m_cROTidxSCModel.get(0,0, 1 ) );
    }
 #endif //#if JVET_C0042_UNIFIED_BINARIZATION
  }
  else
#endif
  if (iNumberOfPassesROT>1) // for only 1 pass no signaling is needed 
  {
      Int idxROT = pcCU->getROTIdx( CHANNEL_TYPE_CHROMA, uiAbsPartIdx );
#if JVET_C0042_UNIFIED_BINARIZATION
      m_pcBinIf->encodeBin(  idxROT ? 1 : 0 , m_cROTidxSCModel.get(0,0, 0) );
            if( idxROT )
            {
                     m_pcBinIf->encodeBin( (idxROT-1) ? 1 : 0 , m_cROTidxSCModel.get(0,0, 2) );

                    if(idxROT >1 )
                    {
                        m_pcBinIf->encodeBin( (idxROT-2) ? 1 : 0, m_cROTidxSCModel.get(0,0, 4) );
                    }
               
            } 
#else
      const UInt uiSymbol0 = (idxROT >>1);
      const UInt uiSymbol1 = (idxROT %2 );
      m_pcBinIf->encodeBin( uiSymbol0, m_cROTidxSCModel.get(0,0, uiDepth ) );
      m_pcBinIf->encodeBin( uiSymbol1, m_cROTidxSCModel.get(0,0, uiDepth ) );
#endif
  }
}
#endif
#endif
/** code merge flag
 * \param pcCU
 * \param uiAbsPartIdx
 * \returns Void
 */
Void TEncSbac::codeMergeFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  const UInt uiSymbol = pcCU->getMergeFlag( uiAbsPartIdx ) ? 1 : 0;
  m_pcBinIf->encodeBin( uiSymbol, *m_cCUMergeFlagExtSCModel.get( 0 ) );

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tMergeFlag: " );
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\tAddress: " );
  DTRACE_CABAC_V( pcCU->getCtuRsAddr() );
  DTRACE_CABAC_T( "\tuiAbsPartIdx: " );
  DTRACE_CABAC_V( uiAbsPartIdx );
  DTRACE_CABAC_T( "\n" );
}

/** code merge index
 * \param pcCU
 * \param uiAbsPartIdx
 * \returns Void
 */
Void TEncSbac::codeMergeIndex( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiUnaryIdx = pcCU->getMergeIndex( uiAbsPartIdx );
  UInt uiNumCand = pcCU->getSlice()->getMaxNumMergeCand();
  if ( uiNumCand > 1 )
  {
    for( UInt ui = 0; ui < uiNumCand - 1; ++ui )
    {
      const UInt uiSymbol = ui == uiUnaryIdx ? 0 : 1;
#if COM16_C806_GEN_MRG_IMPROVEMENT
      m_pcBinIf->encodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, (ui > NUM_MERGE_IDX_EXT_CTX-1? NUM_MERGE_IDX_EXT_CTX-1:ui) ) );
#else
      if ( ui==0 )
      {
        m_pcBinIf->encodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, 0 ) );
      }
      else
      {
        m_pcBinIf->encodeBinEP( uiSymbol );
      }
#endif
      if( uiSymbol == 0 )
      {
        break;
      }
    }
  }
  
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tparseMergeIndex()" );
  DTRACE_CABAC_T( "\tuiMRGIdx= " );
  DTRACE_CABAC_V( pcCU->getMergeIndex( uiAbsPartIdx ) );
  DTRACE_CABAC_T( "\n" );
}

#if VCEG_AZ07_FRUC_MERGE
Void TEncSbac::codeFRUCMgrMode( TComDataCU* pcCU, UInt uiAbsPartIdx , UInt uiPUIdx )
{
  if( !pcCU->getSlice()->getSPS()->getUseFRUCMgrMode() )
    return;

  UInt uiFirstBin = pcCU->getFRUCMgrMode( uiAbsPartIdx ) != FRUC_MERGE_OFF;
#if THIRD_FRUC_MODE_DISABLE_BILATERALFORTID0 //ToDo: low delay'i de dikkate al.

  if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0) > 1 && pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 1)
  {
#endif
      m_pcBinIf->encodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)));
      if (uiFirstBin)
      {
          if (pcCU->getSlice()->isInterP())
          {
              assert(pcCU->getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_TEMPLATE);
          }
          else
          {
              UInt uiSecondBin = pcCU->getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_BILATERALMV;
              m_pcBinIf->encodeBin(uiSecondBin, m_cCUFRUCMESCModel.get(0, 0, pcCU->getCtxFRUCME(uiAbsPartIdx)));
#if THIRD_FRUC_MODE
              if (uiSecondBin == 0  )
              {
                  UInt uiThirdBin = pcCU->getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_TEMPLATE;  //FRUC_MERGE_TEMPLATE = 1 , FRUC_MERGE_TEMPLATE_UNI = 0
                  m_pcBinIf->encodeBin(uiThirdBin, m_cCUFRUCMESCModel_uni.get(0, 0, pcCU->getCtxFRUCME_UNI(uiAbsPartIdx)));
              }
#endif
          }
      }
#if THIRD_FRUC_MODE_DISABLE_BILATERALFORTID0 //ToDo: low delay'i de dikkate al.
  }
  else
  {
      m_pcBinIf->encodeBin(uiFirstBin, m_cCUFRUCMgrModeSCModel.get(0, 0, pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)));

      if (uiFirstBin )
      {
#if THIRD_FRUC_MODE
        UInt uiThirdBin = pcCU->getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_TEMPLATE;  //FRUC_MERGE_TEMPLATE = 1 , FRUC_MERGE_TEMPLATE_UNI = 0
        m_pcBinIf->encodeBin(uiThirdBin, m_cCUFRUCMESCModel_uni.get(0, 0, pcCU->getCtxFRUCME_UNI(uiAbsPartIdx)));
#endif
      }
  }
#endif
}
#endif

#if JVET_C0024_QTBT
Void TEncSbac::codeBTSplitMode ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight)
{
#if TT_LIMIT_SIZEIDX
#if JVET_C0024_SPS_MAX_BT_DEPTH
  TComPic* pcPic        = pcCU->getPic();
  const TComSPS &sps    = pcPic->getPicSym()->getSPS();
  TComSlice* pcSlice    = pcCU->getSlice();
  UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(pcCU->getTextType())?sps.getMaxBTDepthISliceL():sps.getMaxBTDepthISliceC()): sps.getMaxBTDepth();
#if BT_FORCE_EXT
  Int iMaxELv = 
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;
  if (
#if !BT_FORCE_INTRA_ON
    !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx) &&
#endif
    pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
#if BT_FORCE_QTN
    Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(pcCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
    Int iMaxQtLv = 
#if BT_FORCE_INTRA_ON
      pcCU->getSlice()->isIntra() ? MAX_FQT_LV_INTRA :
#endif
      MAX_FQT_LV;
#endif
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
#else
#if BT_FORCE_EXT
  Int iMaxELv = 
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;
  if (
#if !BT_FORCE_INTRA_ON
    !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx) &&
#endif
    pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
#if BT_FORCE_QTN
    Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(pcCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
    Int iMaxQtLv = 
#if BT_FORCE_INTRA_ON
      pcCU->getSlice()->isIntra() ? MAX_FQT_LV_INTRA :
#endif
      MAX_FQT_LV;
#endif
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
  UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(pcCU->getTextType())?MAX_BT_DEPTH:MAX_BT_DEPTH_C): MAX_BT_DEPTH_INTER;
#endif
  UInt uiDepth   = pcCU->getDepth(uiAbsPartIdx);
  UInt uiQTWidth = pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth;
  UInt uiQTHeight = pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth;
  UInt uiDeltaSizeIdx = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight];
  if(uiDeltaSizeIdx>=(uiMaxBTD+1))
    return;
#endif
#if TT_LEAF_NODE
  UInt uiBTDepthCurr = pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);
  UInt uiBTSplitModeLastBtDepth = (uiBTDepthCurr==0)? 0 : pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepthCurr-1);
  if(uiBTSplitModeLastBtDepth & 2) //if tt
      return;
#endif
  UInt uiBTSplitMode = pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight));  //0: no split; 1: hor; 2: ver
#if TT_SYNTAX_NEW
  UInt uiBTDepth1 = pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);
  char cAllowSplit = pcCU->deriveAllowedSplitType(uiAbsPartIdx, uiWidth, uiHeight, uiBTDepth1);
  int iHBT = (cAllowSplit>>3) & 1; 
  int iVBT = (cAllowSplit>>2) & 1; 
  int iHTT = (cAllowSplit>>1) & 1;
  int iVTT = (cAllowSplit   ) & 1;
  int sum  = iHBT + iVBT + iHTT + iVTT;
  UInt uiCtx;
  if(sum!=0)
  {
    uiCtx = pcCU->getCtxBTSplitFlag(uiAbsPartIdx, uiWidth, uiHeight);
    m_pcBinIf->encodeBin(uiBTSplitMode!=0, m_cBTSplitFlagSCModel.get(0, 0, uiCtx));
  }
  else
  {
    assert(uiBTSplitMode==0);
  }
#else
  UInt uiCtx = pcCU->getCtxBTSplitFlag(uiAbsPartIdx, uiWidth, uiHeight);

  //signal bits: no split: 0; hor split: 10; ver split: 11.
  m_pcBinIf->encodeBin(uiBTSplitMode!=0, m_cBTSplitFlagSCModel.get(0, 0, uiCtx));
#endif
  UInt uiMinBTSize = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType())?MIN_BT_SIZE:MIN_BT_SIZE_C): MIN_BT_SIZE_INTER;

#if TT_SYNTAX //dir ctx
  UInt uiBTSCtx = uiWidth==uiHeight ? 0: (uiWidth>uiHeight? 1: 2);
  uiBTSCtx = (uiWidth>(uiHeight<<1))? 3 : ((uiHeight>(uiWidth<<1))? 4 : uiBTSCtx);
  UInt uiBTDepth = pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);
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
  UInt uiSTMode = (uiBTSplitMode>>1) & 1;  //0: BT,   1: TT
  UInt uiSTDir  = (uiBTSplitMode>>2) != 1; //0: hori, 1: vert
#if !BT_FORCE_BF
  assert(uiBTDepth<3);
#endif
  Bool bTreeTypeNeed = true;
#endif

#if TT_LIMIT_SYNTAX
  UInt uiParentType = (uiBTDepth==0)? 0 : pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth-1);
#endif

#if TT_SYNTAX_NEW
  if(uiBTSplitMode != 0)
  {
    if(sum==4)
    {
      m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
      m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
    }
    else if(sum==3)
    {
      m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
      if(iHBT==0 || iHTT==0)
      {
        if(uiSTDir==1)
          m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        else
          assert( uiSTMode == (iHBT==0) );
      }
      else if(iVBT==0 || iVTT==0)
      {
        if(uiSTDir==0)
          m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        else
          assert( uiSTMode == (iVBT==0) );
      }
    }
    else if(sum==2)
    {
      if(iHTT==0 && iVTT==0)
      {
        m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
        assert(uiSTMode == 0);
      }
      else if(iHBT==0 && iHTT==0)
      {
        assert(uiSTDir == 1);
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
      }
      else if(iVBT==0 && iVTT==0)
      {
        assert(uiSTDir == 0);
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
      }
      else if(iVBT==0 && iHTT==0)
      {
        m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
        assert(uiSTMode == (uiSTDir==1));
      }
      else if(iHBT==0 && iVTT==0)
      {
        m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
        assert(uiSTMode == (uiSTDir==0));
      }
      else
        assert(0);
    }
    else//i.e., if(sum==1)
    {
      assert(uiSTDir == (iVBT==1 || iVTT==1));
      assert(uiSTMode== (iHTT==1 || iVTT==1));
    }
  }
#else
  if (uiWidth==uiMinBTSize || uiHeight==uiMinBTSize)
  {
#if TT_SYNTAX
    if(uiBTSplitMode!=0)
    {
      if(uiHeight>=(uiMinBTSize<<2) || uiWidth>=(uiMinBTSize<<2))
      {
#if TT_LIMIT_SYNTAX && TT_LIMIT_SMALL
        if((uiParentType==ST_HORI_TT_S && uiWidth==uiMinBTSize) || (uiParentType==ST_VERT_TT_S && uiHeight==uiMinBTSize))
        {
          bTreeTypeNeed = false;
          assert(uiSTMode==0);
        }
#endif
#if TT_LIMIT_SYNTAX && TT_LIMIT_LARGE
        if((uiParentType==ST_HORI_TT_L && uiWidth==uiMinBTSize) || (uiParentType==ST_VERT_TT_L && uiHeight==uiMinBTSize))
        {
          bTreeTypeNeed = false;
          assert(uiSTMode==0);
        }
#endif
        if(bTreeTypeNeed)
        {
          m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
          DTRACE_CABAC_VL( g_nSymbolCounter++ );
          if(uiWidth==uiMinBTSize){
            DTRACE_CABAC_T( "\tSTModeFlag-W" );}
          else{
            DTRACE_CABAC_T( "\tSTModeFlag-H" );}
          DTRACE_CABAC_T( "\tuiBTSCtx2: ");
          DTRACE_CABAC_V( uiBTSCtx2 );
          DTRACE_CABAC_T( "\tuiSymbol: ");
          DTRACE_CABAC_V( uiSTMode );
          DTRACE_CABAC_T( "\n");
        }
      }
    }
#endif
    return;
  }

#if TT_SYNTAX
  if (uiWidth==(uiMinBTSize<<1) && uiHeight==(uiMinBTSize<<1))
  {
    if(uiBTSplitMode!=0)
    {
      m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag-W2H2" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSTDir );
      DTRACE_CABAC_T( "\n");

      assert( uiBTSplitMode == ST_HORI_BT || uiBTSplitMode == ST_VERT_BT);
    }
    return;
  }

  if (uiWidth==(uiMinBTSize<<1))
  {
    assert(uiHeight>(uiMinBTSize<<1));
    if(uiBTSplitMode!=0)
    {
      m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag-W2H3" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSTDir );
      DTRACE_CABAC_T( "\n");

      if(uiSTDir==1)
        assert( uiBTSplitMode == ST_VERT_BT);
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
          m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
          DTRACE_CABAC_VL( g_nSymbolCounter++ );
          DTRACE_CABAC_T( "\tSTModeFlag-W2H3" );
          DTRACE_CABAC_T( "\tuiBTSCtx2: ");
          DTRACE_CABAC_V( uiBTSCtx2 );
          DTRACE_CABAC_T( "\tuiSymbol: ");
          DTRACE_CABAC_V( uiSTMode );
          DTRACE_CABAC_T( "\n");
        }
      }
    }
    return;
  }

  if (uiHeight==(uiMinBTSize<<1))
  {
    assert(uiWidth>(uiMinBTSize<<1));
    if(uiBTSplitMode!=0)
    {
      m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTDirFlag-W3H2" );
      DTRACE_CABAC_T( "\tuiBTSCtx: ");
      DTRACE_CABAC_V( uiBTSCtx );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSTDir );
      DTRACE_CABAC_T( "\n");

      if(uiSTDir==0)
        assert( uiBTSplitMode == ST_HORI_BT);
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
          m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
          DTRACE_CABAC_VL( g_nSymbolCounter++ );
          DTRACE_CABAC_T( "\tSTModeFlag-W3H2" );
          DTRACE_CABAC_T( "\tuiBTSCtx2: ");
          DTRACE_CABAC_V( uiBTSCtx2 );
          DTRACE_CABAC_T( "\tuiSymbol: ");
          DTRACE_CABAC_V( uiSTMode );
          DTRACE_CABAC_T( "\n");
        }
      }
    }
    return;
  }
#endif

#if (TT_LIMIT_AR>0)
  if (uiWidth>=(uiHeight*TT_LIMIT_AR))
  {
    if(uiBTSplitMode!=0)
    {
      assert(uiSTDir==1);
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_VERT_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_VERT_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-ARBTW" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTMode );
        DTRACE_CABAC_T( "\n");
      }
    }
    return;
  }

  if (uiHeight>=(uiWidth*TT_LIMIT_AR))
  {
    if(uiBTSplitMode!=0)
    {
      assert(uiSTDir==0);
#if TT_LIMIT_SYNTAX
      if(TT_LIMIT_SMALL && uiParentType==ST_HORI_TT_S)
        bTreeTypeNeed = false;
      if(TT_LIMIT_LARGE && uiParentType==ST_HORI_TT_L)
        bTreeTypeNeed = false;
#endif
      if(bTreeTypeNeed)
      {
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-ARBTH" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTMode );
        DTRACE_CABAC_T( "\n");
      }
    }
    return;
  }
#endif

#if JVET_C0024_BT_RMV_REDUNDANT
  if ( pcCU->getSplitConstrain() == 1 )
  {
#if TT_SYNTAX
    if(uiBTSplitMode!=0)
    {
      if(bTreeTypeNeed)
      {
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-C1" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTMode );
        DTRACE_CABAC_T( "\n");
      }

      if(uiSTMode)
      {
#if TT_LIMIT_SYNTAX
        if(TT_LIMIT_SMALL && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
          assert( ((uiParentType==ST_HORI_TT_S)? 2 : 1) == uiSTDir+1);
        else if(TT_LIMIT_LARGE && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
          assert( ((uiParentType==ST_HORI_TT_L)? 2 : 1) == uiSTDir+1);
        else
        {
#endif
        m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTDirFlag-C1" );
        DTRACE_CABAC_T( "\tuiBTSCtx: ");
        DTRACE_CABAC_V( uiBTSCtx );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTDir );
        DTRACE_CABAC_T( "\n");
#if TT_LIMIT_SYNTAX
        }
#endif
      }
      else
        assert( uiBTSplitMode != ST_HORI_BT );
    }
#else
#if TT_TREE_DS1
    assert( uiBTSplitMode != ST_HORI_BT );
#else
    assert( uiBTSplitMode != 1 );
#endif
#endif
    return;
  }
  if ( pcCU->getSplitConstrain() == 2 )
  {
#if TT_SYNTAX
    if(uiBTSplitMode!=0)
    {
      if(bTreeTypeNeed)
      {
        m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTModeFlag-C2" );
        DTRACE_CABAC_T( "\tuiBTSCtx2: ");
        DTRACE_CABAC_V( uiBTSCtx2 );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTMode );
        DTRACE_CABAC_T( "\n");
      }

      if(uiSTMode)
      {
#if TT_LIMIT_SYNTAX
        if(TT_LIMIT_SMALL && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
          assert( ((uiParentType==ST_HORI_TT_S)? 2 : 1) == uiSTDir+1);
        else if(TT_LIMIT_LARGE && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
          assert( ((uiParentType==ST_HORI_TT_L)? 2 : 1) == uiSTDir+1);
        else
        {
#endif
        m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
        DTRACE_CABAC_VL( g_nSymbolCounter++ );
        DTRACE_CABAC_T( "\tSTDirFlag-C2" );
        DTRACE_CABAC_T( "\tuiBTSCtx: ");
        DTRACE_CABAC_V( uiBTSCtx );
        DTRACE_CABAC_T( "\tuiSymbol: ");
        DTRACE_CABAC_V( uiSTDir );
        DTRACE_CABAC_T( "\n");
#if TT_LIMIT_SYNTAX
        }
#endif
      }
      else
        assert( uiBTSplitMode != ST_VERT_BT );
    }
#else
#if TT_TREE_DS1
    assert( uiBTSplitMode != ST_VERT_BT );
#else
    assert( uiBTSplitMode != 2 );
#endif
#endif
    return;
  }
#endif

  if (uiBTSplitMode != 0)
  {
#if TT_SYNTAX
    if(bTreeTypeNeed)
    {
      m_pcBinIf->encodeBin(uiSTMode, m_cBTSplitFlagSCModel.get(0, 0, BT_TYPE_CTX_OFFSET+uiBTSCtx2));
      DTRACE_CABAC_VL( g_nSymbolCounter++ );
      DTRACE_CABAC_T( "\tSTModeFlag" );
      DTRACE_CABAC_T( "\tuiBTSCtx2: ");
      DTRACE_CABAC_V( uiBTSCtx2 );
      DTRACE_CABAC_T( "\tuiSymbol: ");
      DTRACE_CABAC_V( uiSTMode );
      DTRACE_CABAC_T( "\n");
    }
#endif
#if TT_LIMIT_SYNTAX
    if(TT_LIMIT_SMALL && uiSTMode==1 && (uiParentType==ST_HORI_TT_S || uiParentType==ST_VERT_TT_S))
      assert( ((uiParentType==ST_HORI_TT_S)? 2 : 1) == uiSTDir+1);
    else if(TT_LIMIT_LARGE && uiSTMode==1 && (uiParentType==ST_HORI_TT_L || uiParentType==ST_VERT_TT_L))
      assert( ((uiParentType==ST_HORI_TT_L)? 2 : 1) == uiSTDir+1);
    else
    {
#endif
    UInt uiBTSCtx = uiWidth==uiHeight ? 0: (uiWidth>uiHeight? 1: 2);
#if TT_SYNTAX //dir ctx
    uiBTSCtx = (uiWidth>(uiHeight<<1))? 3 : ((uiHeight>(uiWidth<<1))? 4 : uiBTSCtx);
#endif
#if TT_TREE_DS1
    m_pcBinIf->encodeBin(uiSTDir, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
    DTRACE_CABAC_VL( g_nSymbolCounter++ );
    DTRACE_CABAC_T( "\tSTDirFlag" );
    DTRACE_CABAC_T( "\tuiBTSCtx: ");
    DTRACE_CABAC_V( uiBTSCtx );
    DTRACE_CABAC_T( "\tuiSymbol: ");
    DTRACE_CABAC_V( uiSTDir );
    DTRACE_CABAC_T( "\n");
#else
    m_pcBinIf->encodeBin(uiBTSplitMode!=1, m_cBTSplitFlagSCModel.get(0, 0, 3+uiBTSCtx));
#endif
#if TT_LIMIT_SYNTAX
    }
#endif
  }
#endif
}
#endif

Void TEncSbac::codeSplitFlag   ( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  if ( uiDepth == g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] 
  - g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMinQTSize(pcCU->getSlice()->getSliceType(), pcCU->getTextType())])
#else
  if( uiDepth == pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() )
#endif
  {
    return;
  }

  UInt uiCtx           = pcCU->getCtxSplitFlag( uiAbsPartIdx, uiDepth );
  UInt uiCurrSplitFlag = ( pcCU->getDepth( uiAbsPartIdx ) > uiDepth ) ? 1 : 0;

#if COM16_C806_LARGE_CTU
  assert( uiCtx < NUM_SPLIT_FLAG_CTX );
#else
  assert( uiCtx < 3 );
#endif
  m_pcBinIf->encodeBin( uiCurrSplitFlag, m_cCUSplitFlagSCModel.get( 0, 0, uiCtx ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tSplitFlag\n" )
  return;
}

Void TEncSbac::codeTransformSubdivFlag( UInt uiSymbol, UInt uiCtx )
{
  m_pcBinIf->encodeBin( uiSymbol, m_cCUTransSubdivFlagSCModel.get( 0, 0, uiCtx ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseTransformSubdivFlag()" )
  DTRACE_CABAC_T( "\tsymbol=" )
  DTRACE_CABAC_V( uiSymbol )
  DTRACE_CABAC_T( "\tctx=" )
  DTRACE_CABAC_V( uiCtx )
  DTRACE_CABAC_T( "\n" )
}


Void TEncSbac::codeIntraDirLumaAng( TComDataCU* pcCU, UInt absPartIdx, Bool isMultiple
#if VCEG_AZ07_INTRA_65ANG_MODES
                                   , Int* piModes, Int iAboveLeftCase
#endif
                                   )
{
  UInt dir[4],j;
#if VCEG_AZ07_INTRA_65ANG_MODES
  Int preds[4][NUM_MOST_PROBABLE_MODES] = {{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1}};
#else
  Int preds[4][NUM_MOST_PROBABLE_MODES] = {{-1, -1, -1},{-1, -1, -1},{-1, -1, -1},{-1, -1, -1}};
#endif
  Int predIdx[4] ={ -1,-1,-1,-1};
#if JVET_C0024_QTBT
  UInt partNum = 1;
#else
  PartSize mode = pcCU->getPartitionSize( absPartIdx );
  UInt partNum = isMultiple?(mode==SIZE_NxN?4:1):1;
#endif
  UInt partOffset = ( pcCU->getPic()->getNumPartitionsInCtu() >> ( pcCU->getDepth(absPartIdx) << 1 ) ) >> 2;
  
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
  Int aiCase[4]={0,0,0,0};
#endif


#if MLIP_EXTEND_REF_LINE
  Bool bSPSUseMLIP = pcCU->getSlice()->getSPS()->getUseMLIP();//sps flag,
  UInt ExternRef = pcCU->getExternRef(CHANNEL_TYPE_LUMA, absPartIdx);
  if (bSPSUseMLIP)
  {
#if MLIP_IMPROVED_INDEX_CODING
      Int iCurNumRef = MLIP_REF_LINE_NUM;
      UInt tempExternRef = ExternRef;

      if (iCurNumRef == 4)
      {
          //MX: swap the coding code of ExternRef 3 and ExternRef 1; need to be modified
          if (tempExternRef == MLIP_MODI_INDEX)
          {
              tempExternRef = 1;
          }
          else if (tempExternRef == 1)
          {
              tempExternRef = MLIP_MODI_INDEX;
          }
      }
      assert(iCurNumRef == 4);//MX: if not NumLine ==4; need to be modified
      m_pcBinIf->encodeBin((tempExternRef > 0), m_cCUExternRefSCModel.get(0, 0, 0));
#if MLIP_REF_LINE_NUM > 2
      if (tempExternRef) {
          m_pcBinIf->encodeBin((tempExternRef - 1 > 0), m_cCUExternRefSCModel.get(0, 0, 1));
#if MLIP_REF_LINE_NUM > 3
          if (tempExternRef - 1)
              m_pcBinIf->encodeBin(tempExternRef - 2, m_cCUExternRefSCModel.get(0, 0, 2));
#endif
      }
#endif
#else
      m_pcBinIf->encodeBin((ExternRef > 0), m_cCUExternRefSCModel.get(0, 0, 0));
#if MLIP_REF_LINE_NUM > 2
      if (ExternRef) {
          m_pcBinIf->encodeBin((ExternRef - 1 > 0), m_cCUExternRefSCModel.get(0, 0, 1));
#if MLIP_REF_LINE_NUM > 3
          if (ExternRef - 1)
              m_pcBinIf->encodeBin(ExternRef - 2, m_cCUExternRefSCModel.get(0, 0, 2));
#endif
      }
#endif
#endif
  }


#endif


  for (j=0;j<partNum;j++)
  {
    dir[j] = pcCU->getIntraDir( CHANNEL_TYPE_LUMA, absPartIdx+partOffset*j );
#if VCEG_AZ07_INTRA_65ANG_MODES
    if( piModes )
    {
      assert( !isMultiple );
      memcpy( preds[j], piModes, 6*sizeof(Int) );
#if !JVET_C0055_INTRA_MPM
      aiCase[j] = iAboveLeftCase;
#endif
    }
    else
#endif
    pcCU->getIntraDirPredictor(absPartIdx+partOffset*j, preds[j], COMPONENT_Y
#if VCEG_AZ07_INTRA_65ANG_MODES && !JVET_C0055_INTRA_MPM
    , aiCase[j]
#endif
#if MLIP_EXTEND_REF_LINE
        , ExternRef
#endif
    );
    for(UInt i = 0; i < NUM_MOST_PROBABLE_MODES; i++)
    {
      if(dir[j] == preds[j][i])
      {
        predIdx[j] = i;
      }
    }
    m_pcBinIf->encodeBin((predIdx[j] != -1)? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, 0 ) );
  }
  for (j=0;j<partNum;j++)
  {
    if(predIdx[j] != -1)
    {
#if JVET_C0055_INTRA_MPM
      m_pcBinIf->encodeBin( predIdx[j] ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[j][0]] ) );
#elif VCEG_AZ07_INTRA_65ANG_MODES
      m_pcBinIf->encodeBin( predIdx[j] ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM0[aiCase[j]] ) );
#else
      m_pcBinIf->encodeBinEP( predIdx[j] ? 1 : 0 );
#endif
      if (predIdx[j])
      {
#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_C0055_INTRA_MPM
        m_pcBinIf->encodeBin( (predIdx[j]-1) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[j][1]] ) );
        if ( (predIdx[j]-1) )
        {
          m_pcBinIf->encodeBin( (predIdx[j]-2) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, mpmContext[preds[j][2]] ) );
#else
        m_pcBinIf->encodeBin( (predIdx[j]-1) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM1[aiCase[j]] ) );
        if ( (predIdx[j]-1) )
        {
          m_pcBinIf->encodeBin( (predIdx[j]-2) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, uiContextMPM2[aiCase[j]] ) );
#endif

          if (predIdx[j]-2)
          {
            m_pcBinIf->encodeBinEP( (predIdx[j]-3) ? 1 : 0 );
            if (predIdx[j]-3)
            {
              m_pcBinIf->encodeBinEP( (predIdx[j]-4) ? 1 : 0 );
            }
          }
        }
#else
        m_pcBinIf->encodeBinEP( predIdx[j]-1 );
#endif
      }
    }
    else
    {
#if VCEG_AZ07_INTRA_65ANG_MODES
      std::sort(preds[j], preds[j]+6);
#else
      if (preds[j][0] > preds[j][1])
      {
        std::swap(preds[j][0], preds[j][1]);
      }
      if (preds[j][0] > preds[j][2])
      {
        std::swap(preds[j][0], preds[j][2]);
      }
      if (preds[j][1] > preds[j][2])
      {
        std::swap(preds[j][1], preds[j][2]);
      }
#endif
      for(Int i = (Int(NUM_MOST_PROBABLE_MODES) - 1); i >= 0; i--)
      {
        dir[j] = dir[j] > preds[j][i] ? dir[j] - 1 : dir[j];
      }
#if MLIP_EXTEND_REF_LINE
      if (ExternRef)
          dir[j] -= 2;
#endif
#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_E0077_ENHANCED_LM
      assert(dir[j]<(NUM_INTRA_MODE - NUM_INTRA_MODE_NON_ANG - 5));
#else
      assert( dir[j]<(NUM_INTRA_MODE-7) );
#endif

#if JVET_B0051_NON_MPM_MODE
#if JVET_C0024_QTBT
#if JVET_C0024_BT_FIX_TICKET22
      m_pcBinIf->encodeBin(( (dir[j]%4) ==0 ) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, 9 ) ); // flag to indicate if it is selected mode or non-selected mode
#else
      m_pcBinIf->encodeBin(( (dir[j]%4) ==0 ) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, 9/3 ) ); // flag to indicate if it is selected mode or non-selected mode
#endif
#else
      m_pcBinIf->encodeBin(( (dir[j]%4) ==0 ) ? 1 : 0, m_cCUIntraPredSCModel.get( 0, 0, 9+mode/3 ) ); // flag to indicate if it is selected mode or non-selected mode
#endif
      if(dir[j] %4 ==0) 
      {
#if MLIP_EXTEND_REF_LINE
          xWriteTruncBinCode(dir[j] >> 2, 16 - (ExternRef > 0));
#else
        m_pcBinIf->encodeBinsEP( dir[j]>>2, 4 );  // selected mode is 4-bit FLC coded
#endif
      }
      else
      {
        dir[j] -= dir[j]>>2 ;
        dir[j] --;     
#if MLIP_EXTEND_REF_LINE
        xWriteTruncBinCode(dir[j], 45 - (ExternRef > 0));
#else
        xWriteTruncBinCode(dir[j] , 45);  // Non-selected mode is truncated binary coded
#endif
      }
#else
      if( dir[j]>=(NUM_INTRA_MODE-8) )
        m_pcBinIf->encodeBinsEP( dir[j]>>2, 4 );
      else
        m_pcBinIf->encodeBinsEP( dir[j], 6 );
#endif
#else
      m_pcBinIf->encodeBinsEP( dir[j], 5 );
#endif
    }
  }


#if RRC_DWDIP
  const Bool bUseDWDIP = pcCU->getSlice()->getSPS()->getUseDWDIP();

  if (bUseDWDIP &&
    isEnabledForDWDIP(CHANNEL_TYPE_LUMA, pcCU->getWidth(absPartIdx), pcCU->getHeight(absPartIdx), pcCU->getIntraDir(CHANNEL_TYPE_LUMA, absPartIdx))
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    && pcCU->getPDPCIdx(absPartIdx) == 0
#endif
#if !RRC_DWDIP_2Nx2N_CHECK_DISABLE
    && pcCU->getPartitionSize(absPartIdx) == SIZE_2Nx2N
#endif    
#if MLIP_EXTEND_REF_LINE
      && pcCU->getExternRef(CHANNEL_TYPE_LUMA, absPartIdx) == 0
#endif
    )
  {
    UInt uiCtxIdx = 0;
    m_pcBinIf->encodeBin(pcCU->getDWDIPIdx(absPartIdx), m_cDWDIPIdxSCModel.get(0, 0, uiCtxIdx));
  }
#endif

#if RRC_DWDIP_MLR
  if (bUseDWDIP &&
      isEnabledForDWDIPGP(CHANNEL_TYPE_LUMA, pcCU->getWidth(absPartIdx), pcCU->getHeight(absPartIdx), pcCU->getIntraDir(CHANNEL_TYPE_LUMA, absPartIdx))
#if MLIP_EXTEND_REF_LINE
      && pcCU->getExternRef(CHANNEL_TYPE_LUMA, absPartIdx) == 0
#endif
    )
  {
    UInt uiCtxIdx = 0;
    m_pcBinIf->encodeBin(pcCU->getDWDIP_GPIdx(absPartIdx), m_cDWDIP_GPIdxSCModel.get(0, 0, uiCtxIdx));
  }
#endif

  return;
}

#if JVET_E0077_ENHANCED_LM
Void TEncSbac::codeLMModes(TComDataCU* pcCU, UInt uiAbsPartIdx, Int iMode)
{
    Int symbol = -1;
    Int k = 0;
    Int aiLMModeList[10];
    Int iSymbolNum = pcCU->getLMSymbolList(aiLMModeList, uiAbsPartIdx);
    for (k = 0; k < LM_SYMBOL_NUM; k++)
    {
        if (aiLMModeList[k] == iMode || (aiLMModeList[k] == -1 && iMode < LM_CHROMA_IDX))
        {
            symbol = k;
            break;
        }
    }
    assert(symbol >= 0);


    xWriteUnaryMaxSymbol(symbol, &m_cCUChromaPredSCModel.get(0, 0, 1), 1, iSymbolNum - 1);
    return;
}
#endif

Void TEncSbac::codeIntraDirChroma( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiIntraDirChroma = pcCU->getIntraDir( CHANNEL_TYPE_CHROMA, uiAbsPartIdx );
#if JVET_E0062_MULTI_DMS
  Int i = 0;
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
  if (IsLMMode(uiIntraDirChroma) && pcCU->getSlice()->getSPS()->getUseLMChroma())
#else
  if( uiIntraDirChroma == LM_CHROMA_IDX && pcCU->getSlice()->getSPS()->getUseLMChroma() )
#endif
  {
    m_pcBinIf->encodeBin( 0, m_cCUChromaPredSCModel.get( 0, 0, 0 ) );   
#if JVET_E0077_ENHANCED_LM
        Int iCtx = 6;
#if JVET_E0077_MMLM
        if (iBlockSize >= g_aiMMLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
        {
            UInt uiFlag = uiIntraDirChroma == MMLM_CHROMA_IDX;
            m_pcBinIf->encodeBin(uiFlag, m_cCUChromaPredSCModel.get(0, 0, iCtx++));
        }
#endif

#if JVET_E0077_LM_MF
        if (iBlockSize >= g_aiMFLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
        {
            if (((uiIntraDirChroma == LM_CHROMA_IDX)
                || (uiIntraDirChroma >= LM_CHROMA_F1_IDX && uiIntraDirChroma < LM_CHROMA_F1_IDX + LM_FILTER_NUM)))
            {
                m_pcBinIf->encodeBin(uiIntraDirChroma == LM_CHROMA_IDX, m_cCUChromaPredSCModel.get(0, 0, iCtx++));

                if (uiIntraDirChroma >= LM_CHROMA_F1_IDX && uiIntraDirChroma < LM_CHROMA_F1_IDX + LM_FILTER_NUM)
                {
                    Int iLable = uiIntraDirChroma - LM_CHROMA_F1_IDX;
                    m_pcBinIf->encodeBin((iLable >> 1) & 1, m_cCUChromaPredSCModel.get(0, 0, iCtx++));
                    m_pcBinIf->encodeBin(iLable & 1, m_cCUChromaPredSCModel.get(0, 0, iCtx++));
                }
            }
        }

#endif
#endif
  }
  else
#endif
  {
    UInt iDMIdx = 0, uiAllowedChromaDir[NUM_CHROMA_MODE];
    pcCU->getAllowedChromaDir( uiAbsPartIdx, uiAllowedChromaDir );

    for(i = iStartIdx; i < (NUM_DM_MODES + iStartIdx); i++)
    {
      if(uiIntraDirChroma == uiAllowedChromaDir[i])
      {
        iDMIdx = i - iStartIdx;
        break;
      }
    }

#if COM16_C806_LMCHROMA
    if( pcCU->getSlice()->getSPS()->getUseLMChroma() )
    {
      m_pcBinIf->encodeBin( 1, m_cCUChromaPredSCModel.get( 0, 0, 0 ) );
    }
#endif
            
    UInt ictxIdx = 1;
    m_pcBinIf->encodeBin(iDMIdx ? 1 : 0, m_cCUChromaPredSCModel.get(0, 0, ictxIdx));
    UInt uiMaxSymbol = NUM_DM_MODES;
    if (iDMIdx)
    {
      Bool bCodeLast = (uiMaxSymbol > iDMIdx);
      while (--iDMIdx)
      {
        ictxIdx++;
        m_pcBinIf->encodeBin(1, m_cCUChromaPredSCModel.get(0, 0, ictxIdx));
      }
      if (bCodeLast)
      {
        ictxIdx++;
        m_pcBinIf->encodeBin(0, m_cCUChromaPredSCModel.get(0, 0, ictxIdx));
      }
    }
      
  }
#else
  if( uiIntraDirChroma == DM_CHROMA_IDX )
  {
    m_pcBinIf->encodeBin( 0, m_cCUChromaPredSCModel.get( 0, 0, 0 ) );
  }

#if COM16_C806_LMCHROMA
#if JVET_E0077_ENHANCED_LM
  else
  {
      m_pcBinIf->encodeBin(1, m_cCUChromaPredSCModel.get(0, 0, 0));
      if (pcCU->getSlice()->getSPS()->getUseLMChroma())
      {
          codeLMModes(pcCU, uiAbsPartIdx, uiIntraDirChroma);
      }
  }
  if (uiIntraDirChroma == DM_CHROMA_IDX || IsLMMode(uiIntraDirChroma))
  {
      //Do nothing
  }
#else
  else if( uiIntraDirChroma == LM_CHROMA_IDX && pcCU->getSlice()->getSPS()->getUseLMChroma() )
  {
    m_pcBinIf->encodeBin( 1, m_cCUChromaPredSCModel.get( 0, 0, 0 ) );
    m_pcBinIf->encodeBin( 0, m_cCUChromaPredSCModel.get( 0, 0, 1 ) );
  }
#endif

#endif
  else
  {
#if !JVET_E0077_ENHANCED_LM
    m_pcBinIf->encodeBin( 1, m_cCUChromaPredSCModel.get( 0, 0, 0 ) );
#if COM16_C806_LMCHROMA
    if (pcCU->getSlice()->getSPS()->getUseLMChroma())
    {
      m_pcBinIf->encodeBin( 1, m_cCUChromaPredSCModel.get( 0, 0, 1 ));
    }
#endif
#endif
    UInt uiAllowedChromaDir[ NUM_CHROMA_MODE ];
    pcCU->getAllowedChromaDir( uiAbsPartIdx, uiAllowedChromaDir );

#if COM16_C806_LMCHROMA
    for( Int i = 0; i < NUM_CHROMA_MODE - 2; i++ )
#else
    for( Int i = 0; i < NUM_CHROMA_MODE - 1; i++ )
#endif
    {
      if( uiIntraDirChroma == uiAllowedChromaDir[i] )
      {
        uiIntraDirChroma = i;
        break;
      }
    }

    m_pcBinIf->encodeBinsEP( uiIntraDirChroma, 2 );
  }
#endif
  return;
}


Void TEncSbac::codeInterDir( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  const UInt uiInterDir = pcCU->getInterDir( uiAbsPartIdx ) - 1;
  const UInt uiCtx      = pcCU->getCtxInterDir( uiAbsPartIdx );
  ContextModel *pCtx    = m_cCUInterDirSCModel.get( 0 );
#if COM16_C806_LARGE_CTU
  assert( uiCtx < 4 );  // uiCtx=4 is only for the last bin when uiInterDir < 2
#endif

#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
#if JVET_C0024_QTBT
  if (1 )
#else
  if (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag() || pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8 )
#endif
#else
#if JVET_C0024_QTBT
  if(1)
#else
  if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8 )
#endif
#endif
  {
    m_pcBinIf->encodeBin( uiInterDir == 2 ? 1 : 0, *( pCtx + uiCtx ) );
  }

  if (uiInterDir < 2)
  {
    m_pcBinIf->encodeBin( uiInterDir, *( pCtx + 4 ) );
  }

  return;
}

Void TEncSbac::codeRefFrmIdx( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  Int iRefFrame = pcCU->getCUMvField( eRefList )->getRefIdx( uiAbsPartIdx );
  ContextModel *pCtx = m_cCURefPicSCModel.get( 0 );
  m_pcBinIf->encodeBin( ( iRefFrame == 0 ? 0 : 1 ), *pCtx );

  if( iRefFrame > 0 )
  {
    UInt uiRefNum = pcCU->getSlice()->getNumRefIdx( eRefList ) - 2;
    pCtx++;
    iRefFrame--;
    for( UInt ui = 0; ui < uiRefNum; ++ui )
    {
      const UInt uiSymbol = ui == iRefFrame ? 0 : 1;
      if( ui == 0 )
      {
        m_pcBinIf->encodeBin( uiSymbol, *pCtx );
      }
      else
      {
        m_pcBinIf->encodeBinEP( uiSymbol );
      }
      if( uiSymbol == 0 )
      {
        break;
      }
    }
  }
  return;
}

Void TEncSbac::codeMvd( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefList )
{
  if(pcCU->getSlice()->getMvdL1ZeroFlag() && eRefList == REF_PIC_LIST_1 && pcCU->getInterDir(uiAbsPartIdx)==3)
  {
    return;
  }

  const TComCUMvField* pcCUMvField = pcCU->getCUMvField( eRefList );
#if VCEG_AZ07_IMV || VCEG_AZ07_FRUC_MERGE || VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  Int iHor = pcCUMvField->getMvd( uiAbsPartIdx ).getHor();
  Int iVer = pcCUMvField->getMvd( uiAbsPartIdx ).getVer();
#else
  const Int iHor = pcCUMvField->getMvd( uiAbsPartIdx ).getHor();
  const Int iVer = pcCUMvField->getMvd( uiAbsPartIdx ).getVer();
#endif
  ContextModel* pCtx = m_cCUMvdSCModel.get( 0 );

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  assert( iHor == ( iHor >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE ) );
  assert( iVer == ( iVer >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE ) );
  iHor >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
  iVer >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif

#if VCEG_AZ07_IMV
  if( pcCU->getiMVFlag( uiAbsPartIdx ) && pcCU->getSlice()->getSPS()->getIMV() )
  {
    assert( ( iHor & 0x03 ) == 0 && ( iVer & 0x03 ) == 0 );
    iHor >>= 2;
    iVer >>= 2;
  }
#if  JVET_E0076_MULTI_PEL_MVD
    if (pcCU->getiMVFlag( uiAbsPartIdx ) == 2)
    {
      assert( ( iHor % (1 << MULTI_PEL_MVD_BITS)) == 0 && ( iVer % (1 << MULTI_PEL_MVD_BITS) ) == 0 );
      iHor >>= MULTI_PEL_MVD_BITS;
      iVer >>= MULTI_PEL_MVD_BITS;
    }
#endif
#endif

  m_pcBinIf->encodeBin( iHor != 0 ? 1 : 0, *pCtx );
  m_pcBinIf->encodeBin( iVer != 0 ? 1 : 0, *pCtx );

  const Bool bHorAbsGr0 = iHor != 0;
  const Bool bVerAbsGr0 = iVer != 0;
  const UInt uiHorAbs   = 0 > iHor ? -iHor : iHor;
  const UInt uiVerAbs   = 0 > iVer ? -iVer : iVer;

  pCtx++;

  if( bHorAbsGr0 )
  {
    m_pcBinIf->encodeBin( uiHorAbs > 1 ? 1 : 0, *pCtx );
  }

  if( bVerAbsGr0 )
  {
    m_pcBinIf->encodeBin( uiVerAbs > 1 ? 1 : 0, *pCtx );
  }

  if( bHorAbsGr0 )
  {
    if( uiHorAbs > 1 )
    {
      xWriteEpExGolomb( uiHorAbs-2, 1 );
    }

    m_pcBinIf->encodeBinEP( 0 > iHor ? 1 : 0 );
  }

  if( bVerAbsGr0 )
  {
    if( uiVerAbs > 1 )
    {
      xWriteEpExGolomb( uiVerAbs-2, 1 );
    }

    m_pcBinIf->encodeBinEP( 0 > iVer ? 1 : 0 );
  }

  return;
}

Void TEncSbac::codeCrossComponentPrediction( TComTU &rTu, ComponentID compID )
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
    DTRACE_CABAC_VL( g_nSymbolCounter++ )
    DTRACE_CABAC_T("\tparseCrossComponentPrediction()")
    DTRACE_CABAC_T( "\tAddr=" )
    DTRACE_CABAC_V( compID )
    DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
    DTRACE_CABAC_V( uiAbsPartIdx )

    Int alpha = pcCU->getCrossComponentPredictionAlpha( uiAbsPartIdx, compID );
    ContextModel *pCtx = m_cCrossComponentPredictionSCModel.get(0, 0) + ((compID == COMPONENT_Cr) ? (NUM_CROSS_COMPONENT_PREDICTION_CTX >> 1) : 0);
    m_pcBinIf->encodeBin(((alpha != 0) ? 1 : 0), pCtx[0]);

    if (alpha != 0)
    {
      static const Int log2AbsAlphaMinus1Table[8] = { 0, 1, 1, 2, 2, 2, 3, 3 };
      assert(abs(alpha) <= 8);

      if (abs(alpha)>1)
      {
        m_pcBinIf->encodeBin(1, pCtx[1]);
        xWriteUnaryMaxSymbol( log2AbsAlphaMinus1Table[abs(alpha) - 1] - 1, (pCtx + 2), 1, 2 );
      }
      else
      {
        m_pcBinIf->encodeBin(0, pCtx[1]);
      }
      m_pcBinIf->encodeBin( ((alpha < 0) ? 1 : 0), pCtx[4] );
    }
    DTRACE_CABAC_T( "\tAlpha=" )
    DTRACE_CABAC_V( pcCU->getCrossComponentPredictionAlpha( uiAbsPartIdx, compID ) )
    DTRACE_CABAC_T( "\n" )
  }
}

Void TEncSbac::codeDeltaQP( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  Int iDQp  = pcCU->getQP( uiAbsPartIdx ) - pcCU->getRefQP( uiAbsPartIdx );

  Int qpBdOffsetY =  pcCU->getSlice()->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA);
  iDQp = (iDQp + 78 + qpBdOffsetY + (qpBdOffsetY/2)) % (52 + qpBdOffsetY) - 26 - (qpBdOffsetY/2);

  UInt uiAbsDQp = (UInt)((iDQp > 0)? iDQp  : (-iDQp));
  UInt TUValue = min((Int)uiAbsDQp, CU_DQP_TU_CMAX);
  xWriteUnaryMaxSymbol( TUValue, &m_cCUDeltaQpSCModel.get( 0, 0, 0 ), 1, CU_DQP_TU_CMAX);
  if( uiAbsDQp >= CU_DQP_TU_CMAX )
  {
    xWriteEpExGolomb( uiAbsDQp - CU_DQP_TU_CMAX, CU_DQP_EG_k );
  }

  if ( uiAbsDQp > 0)
  {
    UInt uiSign = (iDQp > 0 ? 0 : 1);
    m_pcBinIf->encodeBinEP(uiSign);
  }

  return;
}

/** code chroma qp adjustment, converting from the internal table representation
 * \returns Void
 */
Void TEncSbac::codeChromaQpAdjustment( TComDataCU* cu, UInt absPartIdx )
{
  Int internalIdc = cu->getChromaQpAdj( absPartIdx );
  Int chromaQpOffsetListLen = cu->getSlice()->getPPS()->getPpsRangeExtension().getChromaQpOffsetListLen();
  /* internal_idc == 0 => flag = 0
   * internal_idc > 1 => code idc value (if table size warrents) */
  m_pcBinIf->encodeBin( internalIdc > 0, m_ChromaQpAdjFlagSCModel.get( 0, 0, 0 ) );

  if (internalIdc > 0 && chromaQpOffsetListLen > 1)
  {
    xWriteUnaryMaxSymbol( internalIdc - 1, &m_ChromaQpAdjIdcSCModel.get( 0, 0, 0 ), 0, chromaQpOffsetListLen - 1 );
  }
}

Void TEncSbac::codeQtCbf( TComTU &rTu, const ComponentID compID, const Bool lowestLevel )
{
  TComDataCU* pcCU = rTu.getCU();

  const UInt absPartIdx   = rTu.GetAbsPartIdxTU(compID);
#if !JVET_C0024_QTBT
  const UInt TUDepth      = rTu.GetTransformDepthRel();
#endif
        UInt uiCtx        = pcCU->getCtxQtCbf( rTu, toChannelType(compID) );
  const UInt contextSet   = toChannelType(compID);

#if !JVET_C0024_QTBT
  const UInt width        = rTu.getRect(compID).width;
  const UInt height       = rTu.getRect(compID).height;
  const Bool canQuadSplit = (width >= (MIN_TU_SIZE * 2)) && (height >= (MIN_TU_SIZE * 2));

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

  if ((width != height) && (lowestLevel || !canQuadSplit)) //if sub-TUs are present
  {
    const UInt subTUDepth        = lowestTUDepth + 1;                      //if this is the lowest level of the TU-tree, the sub-TUs are directly below. Otherwise, this must be the level above the lowest level (as specified above)
    const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

    for (UInt subTU = 0; subTU < 2; subTU++)
    {
      const UInt subTUAbsPartIdx = absPartIdx + (subTU * partIdxesPerSubTU);
      const UInt uiCbf           = pcCU->getCbf(subTUAbsPartIdx, compID, subTUDepth);

      m_pcBinIf->encodeBin(uiCbf, m_cCUQtCbfSCModel.get(0, contextSet, uiCtx));

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
  }
  else
#endif
  {
#if JVET_C0024_QTBT
    const UInt uiCbf = pcCU->getCbf( absPartIdx, compID, 0 );
#else
    const UInt uiCbf = pcCU->getCbf( absPartIdx, compID, lowestTUDepth );
#endif
    m_pcBinIf->encodeBin( uiCbf , m_cCUQtCbfSCModel.get( 0, contextSet, uiCtx ) );


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
  }
}

#if VCEG_AZ08_KLT_COMMON
Void TEncSbac::codeKLTFlags(TComTU &rTu, ComponentID component)
{
    TComDataCU* pcCU = rTu.getCU();
    const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();

    if (pcCU->getCUTransquantBypass(uiAbsPartIdx))
    {
        return;
    }

    const TComRectangle &rect = rTu.getRect(component);
    const UInt uiWidth = rect.width;
    const UInt uiHeight = rect.height;
    UInt uiMaxTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
    UInt uiMinTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MIN - 1];
    Bool checkKLTY = ((uiWidth == uiHeight) && (uiWidth <= uiMaxTrWidth) && (uiWidth >= uiMinTrWidth) && (toChannelType(component) == CHANNEL_TYPE_LUMA));
    if (checkKLTY == false)
    {
        return;
    }

    UInt useKLTFlag = pcCU->getKLTFlag(uiAbsPartIdx, component);
    m_pcBinIf->encodeBin(useKLTFlag, m_cKLTFlagSCModel.get(0, toChannelType(component), 0));

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
}
#endif

Void TEncSbac::codeTransformSkipFlags (TComTU &rTu, ComponentID component )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

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
  UInt useTransformSkip = pcCU->getTransformSkip( uiAbsPartIdx,component);
  m_pcBinIf->encodeBin( useTransformSkip, m_cTransformSkipSCModel.get( 0, toChannelType(component), 0 ) );

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
}


/** Code I_PCM information.
 * \param pcCU pointer to CU
 * \param uiAbsPartIdx CU index
 * \returns Void
 */
Void TEncSbac::codeIPCMInfo( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiIPCM = (pcCU->getIPCMFlag(uiAbsPartIdx) == true)? 1 : 0;

  Bool writePCMSampleFlag = pcCU->getIPCMFlag(uiAbsPartIdx);

  m_pcBinIf->encodeBinTrm (uiIPCM);

  if (writePCMSampleFlag)
  {
    m_pcBinIf->encodePCMAlignBits();

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
          UInt sample = pPCMSample[x];
          m_pcBinIf->xWritePCMCode(sample, sampleBits);
        }
        pPCMSample += width;
      }
    }

    m_pcBinIf->resetBac();
  }
}

Void TEncSbac::codeQtRootCbf( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiCbf = pcCU->getQtRootCbf( uiAbsPartIdx );
  UInt uiCtx = 0;
  m_pcBinIf->encodeBin( uiCbf , m_cCUQtRootCbfSCModel.get( 0, 0, uiCtx ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tparseQtRootCbf()" )
  DTRACE_CABAC_T( "\tsymbol=" )
  DTRACE_CABAC_V( uiCbf )
  DTRACE_CABAC_T( "\tctx=" )
  DTRACE_CABAC_V( uiCtx )
  DTRACE_CABAC_T( "\tuiAbsPartIdx=" )
  DTRACE_CABAC_V( uiAbsPartIdx )
  DTRACE_CABAC_T( "\n" )
}

Void TEncSbac::codeQtCbfZero( TComTU & rTu, const ChannelType chType )
{
  // this function is only used to estimate the bits when cbf is 0
  // and will never be called when writing the bistream. do not need to write log
  UInt uiCbf = 0;
  UInt uiCtx = rTu.getCU()->getCtxQtCbf( rTu, chType );

  m_pcBinIf->encodeBin( uiCbf , m_cCUQtCbfSCModel.get( 0, chType, uiCtx ) );
}

Void TEncSbac::codeQtRootCbfZero( )
{
  // this function is only used to estimate the bits when cbf is 0
  // and will never be called when writing the bistream. do not need to write log
  UInt uiCbf = 0;
  UInt uiCtx = 0;
  m_pcBinIf->encodeBin( uiCbf , m_cCUQtRootCbfSCModel.get( 0, 0, uiCtx ) );
}

/** Encode (X,Y) position of the last significant coefficient
 * \param uiPosX     X component of last coefficient
 * \param uiPosY     Y component of last coefficient
 * \param width      Block width
 * \param height     Block height
 * \param component  chroma component ID
 * \param uiScanIdx  scan type (zig-zag, hor, ver)
 * This method encodes the X and Y component within a block of the last significant coefficient.
 */
Void TEncSbac::codeLastSignificantXY( UInt uiPosX, UInt uiPosY, Int width, Int height, ComponentID component, UInt uiScanIdx )
{
  // swap
  if( uiScanIdx == SCAN_VER )
  {
    swap( uiPosX, uiPosY );
    swap( width,  height );
  }

  UInt uiCtxLast;
  UInt uiGroupIdxX    = g_uiGroupIdx[ uiPosX ];
  UInt uiGroupIdxY    = g_uiGroupIdx[ uiPosY ];

  ContextModel *pCtxX = m_cCuCtxLastX.get( 0, toChannelType(component) );
  ContextModel *pCtxY = m_cCuCtxLastY.get( 0, toChannelType(component) );

#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
  Int widthCtx = component ? 4: width;
  const UInt *puiCtxIdxX = g_uiLastCtx + ( g_aucConvertToBit[ widthCtx ] * ( g_aucConvertToBit[ widthCtx ] + 3 ) );

  for( uiCtxLast = 0; uiCtxLast < uiGroupIdxX; uiCtxLast++ )
  {
    if (component)
    {
      m_pcBinIf->encodeBin( 1, *( pCtxX + (uiCtxLast >> g_aucConvertToBit[ width ]) ) );
    }
    else
    {
      m_pcBinIf->encodeBin( 1, *( pCtxX + puiCtxIdxX[ uiCtxLast ] ) );
    }
  }
  if( uiGroupIdxX < g_uiGroupIdx[ width - 1 ])
  {
    if ( component )
    {
      m_pcBinIf->encodeBin( 0, *( pCtxX + (uiCtxLast >> g_aucConvertToBit[ width ]) ) );
    }
    else
    {
      m_pcBinIf->encodeBin( 0, *( pCtxX + puiCtxIdxX[ uiCtxLast ] ) );
    }
  }

  // posY
  Int heightCtx = component ? 4: height;
  const UInt *puiCtxIdxY = g_uiLastCtx + ( g_aucConvertToBit[ heightCtx ] * ( g_aucConvertToBit[ heightCtx ] + 3 ) );
  for( uiCtxLast = 0; uiCtxLast < uiGroupIdxY; uiCtxLast++ )
  {
    if (component)
    {
      m_pcBinIf->encodeBin( 1, *( pCtxY + (uiCtxLast >>  g_aucConvertToBit[ height ])));
    }
    else
    {
      m_pcBinIf->encodeBin( 1, *( pCtxY + puiCtxIdxY[ uiCtxLast ] ) );
    }
  }
  if( uiGroupIdxY < g_uiGroupIdx[ height - 1 ])
  {
    if (component)
    {
      m_pcBinIf->encodeBin( 0, *( pCtxY + (uiCtxLast >> g_aucConvertToBit[ height ]) ) );
    }
    else
    {
      m_pcBinIf->encodeBin( 0, *( pCtxY + puiCtxIdxY[ uiCtxLast ] ) );
    }
   }
#else
  Int blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY;
  getLastSignificantContextParameters(component, width, height, blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY);

  //------------------

  // posX

  for( uiCtxLast = 0; uiCtxLast < uiGroupIdxX; uiCtxLast++ )
  {
    m_pcBinIf->encodeBin( 1, *( pCtxX + blkSizeOffsetX + (uiCtxLast >>shiftX) ) );
  }
  if( uiGroupIdxX < g_uiGroupIdx[ width - 1 ])
  {
    m_pcBinIf->encodeBin( 0, *( pCtxX + blkSizeOffsetX + (uiCtxLast >>shiftX) ) );
  }

  // posY

  for( uiCtxLast = 0; uiCtxLast < uiGroupIdxY; uiCtxLast++ )
  {
    m_pcBinIf->encodeBin( 1, *( pCtxY + blkSizeOffsetY + (uiCtxLast >>shiftY) ) );
  }
  if( uiGroupIdxY < g_uiGroupIdx[ height - 1 ])
  {
    m_pcBinIf->encodeBin( 0, *( pCtxY + blkSizeOffsetY + (uiCtxLast >>shiftY) ) );
  }
#endif
  // EP-coded part

  if ( uiGroupIdxX > 3 )
  {
    UInt uiCount = ( uiGroupIdxX - 2 ) >> 1;
    uiPosX       = uiPosX - g_uiMinInGroup[ uiGroupIdxX ];
    for (Int i = uiCount - 1 ; i >= 0; i-- )
    {
      m_pcBinIf->encodeBinEP( ( uiPosX >> i ) & 1 );
    }
  }
  if ( uiGroupIdxY > 3 )
  {
    UInt uiCount = ( uiGroupIdxY - 2 ) >> 1;
    uiPosY       = uiPosY - g_uiMinInGroup[ uiGroupIdxY ];
    for ( Int i = uiCount - 1 ; i >= 0; i-- )
    {
      m_pcBinIf->encodeBinEP( ( uiPosY >> i ) & 1 );
    }
  }
}

Void TEncSbac::codeCoeffNxN( TComTU &rTu, TCoeff* pcCoef, const ComponentID compID 
#if VCEG_AZ05_ROT_TR    || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  , Int& bCbfCU
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  , Int& iNonZeroCoeffNonTs
#endif
  )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU(compID);
  const TComRectangle &tuRect=rTu.getRect(compID);
  const UInt uiWidth=tuRect.width;
  const UInt uiHeight=tuRect.height;
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());

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
  DTRACE_CABAC_V( uiAbsPartIdx )
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
    std::cerr << "ERROR: codeCoeffNxN was passed a TU with dimensions larger than the maximum allowed size" << std::endl;
    assert(false);
    exit(1);
  }

  // compute number of significant coefficients
  UInt uiNumSig = TEncEntropy::countNonZeroCoeffs(pcCoef, uiWidth * uiHeight);

#if COM16_C806_EMT
  UInt uiTuNumSig = uiNumSig;
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  if( !pcCU->getTransformSkip( uiAbsPartIdx,compID) )
  {
    iNonZeroCoeffNonTs += uiNumSig;
  }
#endif

  if ( uiNumSig == 0 )
  {
    std::cerr << "ERROR: codeCoeffNxN called for empty TU!" << std::endl;
    assert(false);
    exit(1);
  }

  //--------------------------------------------------------------------------------------------------

  //set parameters

  const ChannelType  chType            = toChannelType(compID);
#if JVET_C0024_QTBT
  const UInt         uiLog2BlockWidth  = g_aucConvertToBit[ uiWidth  ] + MIN_CU_LOG2;
  const UInt         uiLog2BlockHeight = g_aucConvertToBit[ uiHeight ] + MIN_CU_LOG2;
#if VCEG_AZ07_CTX_RESIDUALCODING
  const UInt         uiLog2BlockSize = (uiLog2BlockWidth + uiLog2BlockHeight)>>1;
#endif
#else
  const UInt         uiLog2BlockWidth  = g_aucConvertToBit[ uiWidth  ] + 2;
#if !VCEG_AZ07_CTX_RESIDUALCODING || JVET_C0046_ZO_ASSERT
  const UInt         uiLog2BlockHeight = g_aucConvertToBit[ uiHeight ] + 2;
#endif
#endif

  const ChannelType  channelType       = toChannelType(compID);
  const Bool         extendedPrecision = sps.getSpsRangeExtension().getExtendedPrecisionProcessingFlag();


  const Bool         alignCABACBeforeBypass = sps.getSpsRangeExtension().getCabacBypassAlignmentEnabledFlag();
  const Int          maxLog2TrDynamicRange  = sps.getMaxLog2TrDynamicRange(channelType);

  const Int          log2CoeffGroupSize = 2 - ((uiWidth & 0x03) != 0 || (uiHeight & 0x03) != 0);

  Bool beValid;

  {
    Int uiIntraMode = -1;
    const Bool       bIsLuma = isLuma(compID);
    Int isIntra = pcCU->isIntra(uiAbsPartIdx) ? 1 : 0;
    if ( isIntra )
    {
      uiIntraMode = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx );

#if JVET_C0024_QTBT
#if !JVET_E0062_MULTI_DMS
      uiIntraMode = (uiIntraMode==DM_CHROMA_IDX && !bIsLuma) 
        ? (pcCU->getSlice()->isIntra()? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu()+uiAbsPartIdx)
        :pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx)) : uiIntraMode;
#endif
#else
      const UInt partsPerMinCU = 1<<(2*(sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
      uiIntraMode = (uiIntraMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, rTu.GetChromaFormat(), partsPerMinCU)) : uiIntraMode;
#endif
      uiIntraMode = ((rTu.GetChromaFormat() == CHROMA_422) && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiIntraMode] : uiIntraMode;
    }

    Int transformSkip = pcCU->getTransformSkip( uiAbsPartIdx,compID) ? 1 : 0;
    Bool rdpcm_lossy = ( transformSkip && isIntra && ( (uiIntraMode == HOR_IDX) || (uiIntraMode == VER_IDX) ) ) && pcCU->isRDPCMEnabled(uiAbsPartIdx);

    if ( (pcCU->getCUTransquantBypass(uiAbsPartIdx)) || rdpcm_lossy )
    {
      beValid = false;
      if ( (!pcCU->isIntra(uiAbsPartIdx)) && pcCU->isRDPCMEnabled(uiAbsPartIdx))
      {
        codeExplicitRdpcmMode( rTu, compID);
      }
    }
    else
    {
      beValid = pcCU->getSlice()->getPPS()->getSignHideFlag();
    }
  }

  //--------------------------------------------------------------------------------------------------

  if(pcCU->getSlice()->getPPS()->getUseTransformSkip())
  {
    codeTransformSkipFlags(rTu, compID);
    if(pcCU->getTransformSkip(uiAbsPartIdx, compID) && !pcCU->isIntra(uiAbsPartIdx) && pcCU->isRDPCMEnabled(uiAbsPartIdx))
    {
      //  This TU has coefficients and is transform skipped. Check whether is inter coded and if yes encode the explicit RDPCM mode
      codeExplicitRdpcmMode( rTu, compID);

      if(pcCU->getExplicitRdpcmMode(compID, uiAbsPartIdx) != RDPCM_OFF)
      {
        //  Sign data hiding is avoided for horizontal and vertical explicit RDPCM modes
        beValid = false;
      }
    }
  }

  //--------------------------------------------------------------------------------------------------
  const Bool  bUseGolombRiceParameterAdaptation = sps.getSpsRangeExtension().getPersistentRiceAdaptationEnabledFlag();
        UInt &currentGolombRiceStatistic        = m_golombRiceAdaptationStatistics[rTu.getGolombRiceStatisticsIndex(compID)];

  //select scans
  TUEntropyCodingParameters codingParameters;
  getTUEntropyCodingParameters(codingParameters, rTu, compID);

#if VCEG_AZ07_CTX_RESIDUALCODING
  Bool bHor8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_HOR;
  Bool bVer8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_VER;
  Bool bNonZig8x8 = bHor8x8 || bVer8x8; 
#endif

#if VCEG_AZ08_USE_KLT
  if (pcCU->getSlice()->getSPS()->getUseKLT())
  {
#endif
#if VCEG_AZ08_KLT_COMMON
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
      else if (!pcCU->getSlice()->getSPS()->getUseInterKLT() && pcCU->getSlice()->getSPS()->getUseIntraKLT()) //only intra
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
          codeKLTFlags(rTu, compID);
      }
#endif
#if VCEG_AZ08_USE_KLT
  }
#endif
  //----- encode significance map -----

  // Find position of last coefficient
  Int scanPosLast = -1;
  Int posLast;

#if JVET_C0024_QTBT
  const UInt uiCGLog2 = (uiWidth==2 || uiHeight==2) ? MLS_CG_LOG2_HEIGHT-1: MLS_CG_LOG2_HEIGHT;
  const UInt uiCGSizeLog2 = (uiWidth==2 || uiHeight==2) ? MLS_CG_SIZE-2: MLS_CG_SIZE;
#endif

  UInt uiSigCoeffGroupFlag[ MLS_GRP_NUM ];

  memset( uiSigCoeffGroupFlag, 0, sizeof(UInt) * MLS_GRP_NUM );
  do
  {
    posLast = codingParameters.scan[ ++scanPosLast ];

    if( pcCoef[ posLast ] != 0 )
    {
      // get L1 sig map
      UInt uiPosY   = posLast >> uiLog2BlockWidth;
      UInt uiPosX   = posLast - ( uiPosY << uiLog2BlockWidth );

#if JVET_C0024_QTBT
      UInt uiBlkIdx = (codingParameters.widthInGroups * (uiPosY >> uiCGLog2)) + (uiPosX >> uiCGLog2);
#else
      UInt uiBlkIdx = (codingParameters.widthInGroups * (uiPosY >> MLS_CG_LOG2_HEIGHT)) + (uiPosX >> MLS_CG_LOG2_WIDTH);
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING 
      if( bHor8x8 )
      {
        uiBlkIdx = uiPosY >> 1;
      }
      else if( bVer8x8)
      {
        uiBlkIdx = uiPosX >> 1;
      }
#endif
      uiSigCoeffGroupFlag[ uiBlkIdx ] = 1;

      uiNumSig--;
    }
  } while ( uiNumSig > 0 );

  // Code position of last coefficient
  Int posLastY = posLast >> uiLog2BlockWidth;
  Int posLastX = posLast - ( posLastY << uiLog2BlockWidth );
  codeLastSignificantXY(posLastX, posLastY, uiWidth, uiHeight, compID, codingParameters.scanType);

 
#if JVET_C0046_ZO_ASSERT && JVET_C0046_ZO_ASSERT_LAST_COEF
  if ( ((uiWidth > ZERO_OUT_TH) || (uiHeight > ZERO_OUT_TH)) 
    && (!pcCU->getTransformSkip(compID) && !pcCU->getCUTransquantBypass(uiAbsPartIdx)))
  {
     // last coeff shall be in the low freqecy domain
     assert((posLastX < ZERO_OUT_TH) && (posLastY < ZERO_OUT_TH));
  }
#endif


  //===== code significance flag =====
#if VCEG_AZ07_CTX_RESIDUALCODING
  UInt uiPrevGRParam = 0; 
  ContextModel * const baseCoeffGroupCtx = m_cCUSigCoeffGroupSCModel.get( 0, chType );
#if JVET_C0024_QTBT
  UInt uiOffsetonTU = uiLog2BlockSize<=2 ? 0: NUM_SIG_FLAG_CTX_LUMA_TU << ( min(1, (Int)(uiLog2BlockSize - 3)) );
#else
  UInt uiOffsetonTU = uiLog2BlockWidth==2 ? 0: NUM_SIG_FLAG_CTX_LUMA_TU << ( min(1, (Int)(uiLog2BlockWidth - 3)) );
#endif
  ContextModel * const baseCtx = (chType == CHANNEL_TYPE_LUMA) ? m_cCUSigSCModel.get( 0, 0 ) + uiOffsetonTU : m_cCUSigSCModel.get( 0, 0 ) + NUM_SIG_FLAG_CTX_LUMA;
  ContextModel * const greXCtx = (chType == CHANNEL_TYPE_LUMA) ? m_cCUOneSCModel.get( 0, 0 ) : m_cCUOneSCModel.get( 0, 0 ) + NUM_ONE_FLAG_CTX_LUMA;
#else
  ContextModel * const baseCoeffGroupCtx = m_cCUSigCoeffGroupSCModel.get( 0, chType );
  ContextModel * const baseCtx = m_cCUSigSCModel.get( 0, 0 ) + getSignificanceMapContextOffset(compID);
#endif

#if JVET_C0024_QTBT
  const Int  iLastScanSet  = scanPosLast >> uiCGSizeLog2;
#else
  const Int  iLastScanSet  = scanPosLast >> MLS_CG_SIZE;
#endif

#if !VCEG_AZ07_CTX_RESIDUALCODING
  UInt c1                  = 1;
#endif
  UInt uiGoRiceParam       = 0;
  Int  iScanPosSig         = scanPosLast;

  for( Int iSubSet = iLastScanSet; iSubSet >= 0; iSubSet-- )
  {
    Int numNonZero = 0;
#if JVET_C0024_QTBT
    Int  iSubPos   = iSubSet << uiCGSizeLog2;
#else
    Int  iSubPos   = iSubSet << MLS_CG_SIZE;
#endif
#if !VCEG_AZ07_CTX_RESIDUALCODING
    uiGoRiceParam  = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;
#endif
    Bool updateGolombRiceStatistics = bUseGolombRiceParameterAdaptation; //leave the statistics at 0 when not using the adaptation system
    UInt coeffSigns = 0;

    Int absCoeff[1 << MLS_CG_SIZE];
#if VCEG_AZ07_CTX_RESIDUALCODING
    Int pos     [1 << MLS_CG_SIZE];
    UInt ctxG1     = 0;
    UInt ctxG2     = 0;
#endif

    Int lastNZPosInCG  = -1;
#if JVET_C0024_QTBT
    Int firstNZPosInCG = 1 << uiCGSizeLog2;
#else
    Int firstNZPosInCG = 1 << MLS_CG_SIZE;
#endif

    Bool escapeDataPresentInGroup = false;

    if( iScanPosSig == scanPosLast )
    {
      absCoeff[ 0 ] = Int(abs( pcCoef[ posLast ] ));
#if VCEG_AZ07_CTX_RESIDUALCODING
      pos[numNonZero]= posLastX + (posLastY<< uiLog2BlockWidth);
#endif
      coeffSigns    = ( pcCoef[ posLast ] < 0 );
      numNonZero    = 1;
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
      bCbfCU += abs(pcCoef[posLast]);
#endif
      lastNZPosInCG  = iScanPosSig;
      firstNZPosInCG = iScanPosSig;
      iScanPosSig--;
    }

    // encode significant_coeffgroup_flag
    Int iCGBlkPos = codingParameters.scanCG[ iSubSet ];
    Int iCGPosY   = iCGBlkPos / codingParameters.widthInGroups;
    Int iCGPosX   = iCGBlkPos - (iCGPosY * codingParameters.widthInGroups);

#if VCEG_AZ07_CTX_RESIDUALCODING
    if(bNonZig8x8 )
    {
      iCGPosY = (bHor8x8 ? iCGBlkPos : 0);
      iCGPosX = (bVer8x8 ? iCGBlkPos : 0);
    }
#endif

    if( iSubSet == iLastScanSet || iSubSet == 0)
    {
      uiSigCoeffGroupFlag[ iCGBlkPos ] = 1;
    }
#if COM16_C806_T64  && !JVET_C0024_QTBT && !JVET_C0046_ZO_ASSERT
    else if( iCGPosY>=(ZERO_OUT_TH>>log2CoeffGroupSize) || iCGPosX>=(ZERO_OUT_TH>> log2CoeffGroupSize) )
    {
      assert( 0 == uiSigCoeffGroupFlag[ iCGBlkPos ] );
    }
#endif
#if JVET_C0046_ZO_ASSERT && JVET_C0046_ZO_ASSERT_CODED_SBK_FLAG
    else if ((iCGPosY >= (ZERO_OUT_TH >> log2CoeffGroupSize) || iCGPosX >= (ZERO_OUT_TH >> log2CoeffGroupSize))
      && (!pcCU->getTransformSkip(compID) && !pcCU->getCUTransquantBypass(uiAbsPartIdx)))
    {
      // coded_sbk_flag(iCGX,iCGY) shall be equal to 0
      assert(0 == uiSigCoeffGroupFlag[iCGBlkPos]);
    }
#endif
    else
    {
      UInt uiSigCoeffGroup   = (uiSigCoeffGroupFlag[ iCGBlkPos ] != 0);
#if VCEG_AZ07_CTX_RESIDUALCODING
      UInt uiCtxSig  = TComTrQuant::getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups, codingParameters.scanType );
#else
      UInt uiCtxSig  = TComTrQuant::getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups );
#endif
      m_pcBinIf->encodeBin( uiSigCoeffGroup, baseCoeffGroupCtx[ uiCtxSig ] );
    }

    // encode significant_coeff_flag
    if( uiSigCoeffGroupFlag[ iCGBlkPos ] )
    {
#if !VCEG_AZ07_CTX_RESIDUALCODING
      const Int patternSigCtx = TComTrQuant::calcPatternSigCtx(uiSigCoeffGroupFlag, iCGPosX, iCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups);
#endif

      UInt uiBlkPos, uiSig, uiCtxSig;
      for( ; iScanPosSig >= iSubPos; iScanPosSig-- )
      {
        uiBlkPos  = codingParameters.scan[ iScanPosSig ];
        uiSig     = (pcCoef[ uiBlkPos ] != 0);
        if( iScanPosSig > iSubPos || iSubSet == 0 || numNonZero )
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          uiCtxSig  = TComTrQuant::getGrtZeroCtxInc( pcCoef, codingParameters.scan[iScanPosSig], uiWidth, uiHeight, chType );
#else
          uiCtxSig  = TComTrQuant::getSigCtxInc( patternSigCtx, codingParameters, iScanPosSig, uiLog2BlockWidth, uiLog2BlockHeight, chType );
#endif
          m_pcBinIf->encodeBin( uiSig, baseCtx[ uiCtxSig ] );
        }
        if( uiSig )
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          pos[ numNonZero ]      = uiBlkPos;  
#endif
          absCoeff[ numNonZero ] = Int(abs( pcCoef[ uiBlkPos ] ));
          coeffSigns = 2 * coeffSigns + ( pcCoef[ uiBlkPos ] < 0 );
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
          bCbfCU += absCoeff[numNonZero];
#endif
          numNonZero++;
          if( lastNZPosInCG == -1 )
          {
            lastNZPosInCG = iScanPosSig;
          }
          firstNZPosInCG = iScanPosSig;
        }
      }
    }
    else
    {
      iScanPosSig = iSubPos - 1;
    }

    if( numNonZero > 0 )
    {
      Bool signHidden = ( lastNZPosInCG - firstNZPosInCG >= SBH_THRESHOLD );

#if !VCEG_AZ07_CTX_RESIDUALCODING
      const UInt uiCtxSet = getContextSetIndex(compID, iSubSet, (c1 == 0));
      c1 = 1;

      ContextModel *baseCtxMod = m_cCUOneSCModel.get( 0, 0 ) + (NUM_ONE_FLAG_CTX_PER_SET * uiCtxSet);
#endif

      Int numC1Flag = min(numNonZero, C1FLAG_NUMBER);
      Int firstC2FlagIdx = -1;
      for( Int idx = 0; idx < numC1Flag; idx++ )
      {
        UInt uiSymbol = absCoeff[ idx ] > 1;
#if VCEG_AZ07_CTX_RESIDUALCODING 
        if(idx || iSubSet != iLastScanSet)
        {
          ctxG1 = TComTrQuant::getGrtOneCtxInc( pcCoef, pos[idx], uiWidth, uiHeight, chType );
        }
        m_pcBinIf->encodeBin( uiSymbol, greXCtx[ ctxG1 ] );
#else
        m_pcBinIf->encodeBin( uiSymbol, baseCtxMod[c1] );
#endif
        if( uiSymbol )
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
#if !VCEG_AZ07_CTX_RESIDUALCODING
        else if( (c1 < 3) && (c1 > 0) )
        {
          c1++;
        }
#endif
      }
#if !VCEG_AZ07_CTX_RESIDUALCODING
      if (c1 == 0)
      {
        baseCtxMod = m_cCUAbsSCModel.get( 0, 0 ) + (NUM_ABS_FLAG_CTX_PER_SET * uiCtxSet);
#endif
        if ( firstC2FlagIdx != -1)
        {
          UInt symbol = absCoeff[ firstC2FlagIdx ] > 2;
#if VCEG_AZ07_CTX_RESIDUALCODING
          if( firstC2FlagIdx || iSubSet != iLastScanSet )
          {
            ctxG2 = TComTrQuant::getGrtTwoCtxInc( pcCoef, pos[firstC2FlagIdx], uiWidth, uiHeight, chType );
          }
          m_pcBinIf->encodeBin( symbol, greXCtx[ ctxG2 ] );
#else
          m_pcBinIf->encodeBin( symbol, baseCtxMod[0] );
#endif
          if (symbol != 0)
          {
            escapeDataPresentInGroup = true;
          }
        }
#if !VCEG_AZ07_CTX_RESIDUALCODING
      }
#endif
      escapeDataPresentInGroup = escapeDataPresentInGroup || (numNonZero > C1FLAG_NUMBER);

      if (escapeDataPresentInGroup && alignCABACBeforeBypass)
      {
        m_pcBinIf->align();
      }
#if !VCEG_AZ07_CTX_RESIDUALCODING
      if( beValid && signHidden )
      {
        m_pcBinIf->encodeBinsEP( (coeffSigns >> 1), numNonZero-1 );
      }
      else
      {
        m_pcBinIf->encodeBinsEP( coeffSigns, numNonZero );
      }
#endif
      Int iFirstCoeff2 = 1;
      if (escapeDataPresentInGroup)
      {
        for ( Int idx = 0; idx < numNonZero; idx++ )
        {
          UInt baseLevel  = (idx < C1FLAG_NUMBER)? (2 + iFirstCoeff2 ) : 1;

          if( absCoeff[ idx ] >= baseLevel)
          {
            const UInt escapeCodeValue = absCoeff[idx] - baseLevel;
#if VCEG_AZ07_CTX_RESIDUALCODING
            if( updateGolombRiceStatistics && iSubSet == iLastScanSet )
            {
              uiGoRiceParam = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;
            }
            else
            {
              uiGoRiceParam = TComTrQuant::getRemainCoeffCtxInc( pcCoef, pos[idx], uiWidth, uiHeight);
              if(bUseGolombRiceParameterAdaptation)
              {
                uiGoRiceParam = max(uiPrevGRParam , uiGoRiceParam);
              }
            }
            uiPrevGRParam = max(0, (Int)uiGoRiceParam - 1);
#endif
            xWriteCoefRemainExGolomb( escapeCodeValue, uiGoRiceParam, extendedPrecision, maxLog2TrDynamicRange );
#if !VCEG_AZ07_CTX_RESIDUALCODING 
            if (absCoeff[idx] > (3 << uiGoRiceParam))
            {
              uiGoRiceParam = bUseGolombRiceParameterAdaptation ? (uiGoRiceParam + 1) : (std::min<UInt>((uiGoRiceParam + 1), 4));
            }
#endif
            if (updateGolombRiceStatistics)
            {
              const UInt initialGolombRiceParameter = currentGolombRiceStatistic / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;

              if (escapeCodeValue >= (3 << initialGolombRiceParameter))
              {
                currentGolombRiceStatistic++;
              }
              else if (((escapeCodeValue * 2) < (1 << initialGolombRiceParameter)) && (currentGolombRiceStatistic > 0))
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
#if JVET_C0024_QTBT
      if( beValid && signHidden && uiWidth>=4 && uiHeight>=4)
#else
      if( beValid && signHidden )
#endif
      {
        m_pcBinIf->encodeBinsEP( (coeffSigns >> 1), numNonZero-1 );
      }
      else
      {
        m_pcBinIf->encodeBinsEP( coeffSigns, numNonZero );
      }
#endif
    }
  }
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  printSBACCoeffData(posLastX, posLastY, uiWidth, uiHeight, compID, uiAbsPartIdx, codingParameters.scanType, pcCoef, pcCU->getSlice()->getFinalized());
#endif

#if COM16_C806_EMT
  if ( !pcCU->getTransformSkip( uiAbsPartIdx, compID) && compID == COMPONENT_Y )
  {
    if( pcCU->getEmtCuFlag( uiAbsPartIdx ) && pcCU->isIntra( uiAbsPartIdx ) )
    {
      if ( uiTuNumSig>g_iEmtSigNumThr )
      {
        codeEmtTuIdx( pcCU, uiAbsPartIdx, rTu.GetTransformDepthTotal() ); 
      }
      else
      {
        assert( pcCU->getEmtTuIdx( uiAbsPartIdx )==0 );
      }
    }
    if( pcCU->getEmtCuFlag( uiAbsPartIdx ) && !pcCU->isIntra( uiAbsPartIdx ) )
    {
      codeEmtTuIdx( pcCU, uiAbsPartIdx, rTu.GetTransformDepthTotal() ); 
    }
  }
#endif
  return;
}
#if CNN_BASED_INLOOP_FILTER
Void TEncSbac::codeCNNLFBlkParam(ComponentID compIdx, CNNLFBlkParam* &ctuParam, TComDataCU *pcCU)
{
	const Bool bIsFirstCompOfChType = (getFirstComponentOfChannel(toChannelType(compIdx)) == compIdx);
	UInt uiCtxConvNet = pcCU->getCNNLoopFilterCtxFlag(compIdx, ctuParam);
	if (bIsFirstCompOfChType)
	{
		m_pcBinIf->encodeBin(UInt(ctuParam[pcCU->getCtuRsAddr()][compIdx]), m_cCUCNNLoopFilterFlagSCModel.get(0, 0, uiCtxConvNet));
	}
}
#endif
/** code SAO offset sign
 * \param code sign value
 */
Void TEncSbac::codeSAOSign( UInt code )
{
  m_pcBinIf->encodeBinEP( code );
}

Void TEncSbac::codeSaoMaxUvlc    ( UInt code, UInt maxSymbol )
{
  if (maxSymbol == 0)
  {
    return;
  }

  Int i;
  Bool bCodeLast = ( maxSymbol > code );

  if ( code == 0 )
  {
    m_pcBinIf->encodeBinEP( 0 );
  }
  else
  {
    m_pcBinIf->encodeBinEP( 1 );
    for ( i=0; i<code-1; i++ )
    {
      m_pcBinIf->encodeBinEP( 1 );
    }
    if( bCodeLast )
    {
      m_pcBinIf->encodeBinEP( 0 );
    }
  }
}

/** Code SAO EO class or BO band position
 */
Void TEncSbac::codeSaoUflc       ( UInt uiLength, UInt uiCode )
{
  m_pcBinIf->encodeBinsEP ( uiCode, uiLength );
}

/** Code SAO merge flags
 */
Void TEncSbac::codeSaoMerge       ( UInt uiCode )
{
  m_pcBinIf->encodeBin(((uiCode == 0) ? 0 : 1),  m_cSaoMergeSCModel.get( 0, 0, 0 ));
}

/** Code SAO type index
 */
Void TEncSbac::codeSaoTypeIdx       ( UInt uiCode)
{
  if (uiCode == 0)
  {
    m_pcBinIf->encodeBin( 0, m_cSaoTypeIdxSCModel.get( 0, 0, 0 ) );
  }
  else
  {
    m_pcBinIf->encodeBin( 1, m_cSaoTypeIdxSCModel.get( 0, 0, 0 ) );
    m_pcBinIf->encodeBinEP( uiCode == 1 ? 0 : 1 );
  }
}

Void TEncSbac::codeSAOOffsetParam(ComponentID compIdx, SAOOffset& ctbParam, Bool sliceEnabled, const Int channelBitDepth)
{
  UInt uiSymbol;
  if(!sliceEnabled)
  {
    assert(ctbParam.modeIdc == SAO_MODE_OFF);
    return;
  }
  const Bool bIsFirstCompOfChType = (getFirstComponentOfChannel(toChannelType(compIdx)) == compIdx);

  //type
  if(bIsFirstCompOfChType)
  {
    //sao_type_idx_luma or sao_type_idx_chroma
    if(ctbParam.modeIdc == SAO_MODE_OFF)
    {
      uiSymbol =0;
    }
    else if(ctbParam.typeIdc == SAO_TYPE_BO) //BO
    {
      uiSymbol = 1;
    }
    else
    {
      assert(ctbParam.typeIdc < SAO_TYPE_START_BO); //EO
      uiSymbol = 2;
    }
    codeSaoTypeIdx(uiSymbol);
  }

  if(ctbParam.modeIdc == SAO_MODE_NEW)
  {
    Int numClasses = (ctbParam.typeIdc == SAO_TYPE_BO)?4:NUM_SAO_EO_CLASSES;
    Int offset[4];
    Int k=0;
    for(Int i=0; i< numClasses; i++)
    {
      if(ctbParam.typeIdc != SAO_TYPE_BO && i == SAO_CLASS_EO_PLAIN)
      {
        continue;
      }
      Int classIdx = (ctbParam.typeIdc == SAO_TYPE_BO)?(  (ctbParam.typeAuxInfo+i)% NUM_SAO_BO_CLASSES   ):i;
      offset[k] = ctbParam.offset[classIdx];
      k++;
    }

    const Int  maxOffsetQVal = TComSampleAdaptiveOffset::getMaxOffsetQVal(channelBitDepth);
    for(Int i=0; i< 4; i++)
    {
      codeSaoMaxUvlc((offset[i]<0)?(-offset[i]):(offset[i]),  maxOffsetQVal ); //sao_offset_abs
    }


    if(ctbParam.typeIdc == SAO_TYPE_BO)
    {
      for(Int i=0; i< 4; i++)
      {
        if(offset[i] != 0)
        {
          codeSAOSign((offset[i]< 0)?1:0);
        }
      }

      codeSaoUflc(NUM_SAO_BO_CLASSES_LOG2, ctbParam.typeAuxInfo ); //sao_band_position
    }
    else //EO
    {
      if(bIsFirstCompOfChType)
      {
        assert(ctbParam.typeIdc - SAO_TYPE_START_EO >=0);
        codeSaoUflc(NUM_SAO_EO_TYPES_LOG2, ctbParam.typeIdc - SAO_TYPE_START_EO ); //sao_eo_class_luma or sao_eo_class_chroma
      }
    }

  }
}


Void TEncSbac::codeSAOBlkParam(SAOBlkParam& saoBlkParam, const BitDepths &bitDepths
                              , Bool* sliceEnabled
                              , Bool leftMergeAvail
                              , Bool aboveMergeAvail
                              , Bool onlyEstMergeInfo // = false
                              )
{

  Bool isLeftMerge = false;
  Bool isAboveMerge= false;

  if(leftMergeAvail)
  {
    isLeftMerge = ((saoBlkParam[COMPONENT_Y].modeIdc == SAO_MODE_MERGE) && (saoBlkParam[COMPONENT_Y].typeIdc == SAO_MERGE_LEFT));
    codeSaoMerge( isLeftMerge?1:0  ); //sao_merge_left_flag
  }

  if( aboveMergeAvail && !isLeftMerge)
  {
    isAboveMerge = ((saoBlkParam[COMPONENT_Y].modeIdc == SAO_MODE_MERGE) && (saoBlkParam[COMPONENT_Y].typeIdc == SAO_MERGE_ABOVE));
    codeSaoMerge( isAboveMerge?1:0  ); //sao_merge_left_flag
  }

  if(onlyEstMergeInfo)
  {
    return; //only for RDO
  }

  if(!isLeftMerge && !isAboveMerge) //not merge mode
  {
    for(Int compIdx=0; compIdx < MAX_NUM_COMPONENT; compIdx++)
    {
      codeSAOOffsetParam(ComponentID(compIdx), saoBlkParam[compIdx], sliceEnabled[compIdx], bitDepths.recon[toChannelType(ComponentID(compIdx))]);
    }
  }
}

#if JVET_D0123_ME_CTX_LUT_BITS
Void TEncSbac::estPuMeBit  (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  estMvdBit            (pcEstPuMeBitsSbac);
  estMvpIdxBit         (pcEstPuMeBitsSbac);
  estRefIdxBit         (pcEstPuMeBitsSbac);
  estMrgFlagBit        (pcEstPuMeBitsSbac);
  estMrgIdxBit         (pcEstPuMeBitsSbac);  
  estInterDirBit       (pcEstPuMeBitsSbac);  
#if VCEG_AZ07_FRUC_MERGE
  estFrucModeBit     (pcEstPuMeBitsSbac);  
#endif
#if COM16_C1016_AFFINE
  estAffineFlagBit     (pcEstPuMeBitsSbac);
#endif
#if MRG_CONSIDER_COMPLEX
  estMrgComplexFlagBit(pcEstPuMeBitsSbac);
  estMrgComplexIdxBit(pcEstPuMeBitsSbac);
#endif
}

Void TEncSbac::estMvdBit            (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUMvdSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MV_RES_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mvdBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->mvdBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}

Void TEncSbac::estMvpIdxBit         (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac){
  ContextModel *pCtx = m_cMVPIdxSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MVP_IDX_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mvpIdxBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->mvpIdxBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}

Void TEncSbac::estRefIdxBit          (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCURefPicSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_REF_NO_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->refIdxBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->refIdxBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}
#if VCEG_AZ07_FRUC_MERGE
Void TEncSbac::estFrucModeBit      (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUFRUCMgrModeSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_FRUCMGRMODE_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->frucMrgBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->frucMrgBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
  pCtx = m_cCUFRUCMESCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_FRUCME_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->frucMeBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->frucMeBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}
#endif
Void TEncSbac::estMrgFlagBit        (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUMergeFlagExtSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MERGE_FLAG_EXT_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mrgFlagBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->mrgFlagBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}
Void TEncSbac::estMrgIdxBit         (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUMergeIdxExtSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MERGE_IDX_EXT_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mrgIdxBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->mrgIdxBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}
#if COM16_C1016_AFFINE
Void TEncSbac::estAffineFlagBit     (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUAffineFlagSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_AFFINE_FLAG_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->affineFlagBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits(0);
    pcEstPuMeBitsSbac->affineFlagBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits(1);
  }
}
#endif
#if MRG_CONSIDER_COMPLEX
Void TEncSbac::estMrgComplexFlagBit(estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUMergeFlagExtSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MERGE_COMPLEX_FLAG_EXT_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mrgComplexFlagBits[uiCtxInc][0] = pCtx[uiCtxInc].getEntropyBits(0);
    pcEstPuMeBitsSbac->mrgComplexFlagBits[uiCtxInc][1] = pCtx[uiCtxInc].getEntropyBits(1);
  }
}

Void TEncSbac::estMrgComplexIdxBit(estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUMergeIdxExtSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_MERGE_COMPLEX_IDX_EXT_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->mrgComplexIdxBits[uiCtxInc][0] = pCtx[uiCtxInc].getEntropyBits(0);
    pcEstPuMeBitsSbac->mrgComplexIdxBits[uiCtxInc][1] = pCtx[uiCtxInc].getEntropyBits(1);
  }
}
#endif
#if VCEG_AZ07_IMV
Void TEncSbac::estIMVFlagBit        (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUiMVFlagSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_IMV_FLAG_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->iMVFlagBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits(0);
    pcEstPuMeBitsSbac->iMVFlagBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits(1);
  }
}
#endif
Void TEncSbac::estInterDirBit     (estPuMeBitsSbacStruct* pcEstPuMeBitsSbac)
{
  ContextModel *pCtx = m_cCUInterDirSCModel.get(0);
  for (UInt uiCtxInc = 0; uiCtxInc < NUM_INTER_DIR_CTX; uiCtxInc++)
  {
    pcEstPuMeBitsSbac->interDirBits[uiCtxInc][0] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstPuMeBitsSbac->interDirBits[uiCtxInc][1] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}
#endif
/*!
 ****************************************************************************
 * \brief
 *   estimate bit cost for CBP, significant map and significant coefficients
 ****************************************************************************
 */
Void TEncSbac::estBit( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
  , UInt uiScanIdx
#endif
  )
{
  estCBFBit( pcEstBitsSbac );

  estSignificantCoeffGroupMapBit( pcEstBitsSbac, chType );

  // encode significance map
  estSignificantMapBit( pcEstBitsSbac, width, height, chType );

  // encode last significant position
  estLastSignificantPositionBit( pcEstBitsSbac, width, height, chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
    , uiScanIdx
#endif
    );

  // encode significant coefficients
  estSignificantCoefficientsBit( pcEstBitsSbac, chType );

  memcpy(pcEstBitsSbac->golombRiceAdaptationStatistics, m_golombRiceAdaptationStatistics, (sizeof(UInt) * RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS));
}

/*!
 ****************************************************************************
 * \brief
 *    estimate bit cost for each CBP bit
 ****************************************************************************
 */
Void TEncSbac::estCBFBit( estBitsSbacStruct* pcEstBitsSbac )
{
  ContextModel *pCtx = m_cCUQtCbfSCModel.get( 0 );

  for( UInt uiCtxInc = 0; uiCtxInc < (NUM_QT_CBF_CTX_SETS * NUM_QT_CBF_CTX_PER_SET); uiCtxInc++ )
  {
    pcEstBitsSbac->blockCbpBits[ uiCtxInc ][ 0 ] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstBitsSbac->blockCbpBits[ uiCtxInc ][ 1 ] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }

  pCtx = m_cCUQtRootCbfSCModel.get( 0 );

  for( UInt uiCtxInc = 0; uiCtxInc < 4; uiCtxInc++ )
  {
    pcEstBitsSbac->blockRootCbpBits[ uiCtxInc ][ 0 ] = pCtx[ uiCtxInc ].getEntropyBits( 0 );
    pcEstBitsSbac->blockRootCbpBits[ uiCtxInc ][ 1 ] = pCtx[ uiCtxInc ].getEntropyBits( 1 );
  }
}


/*!
 ****************************************************************************
 * \brief
 *    estimate SAMBAC bit cost for significant coefficient group map
 ****************************************************************************
 */
Void TEncSbac::estSignificantCoeffGroupMapBit( estBitsSbacStruct* pcEstBitsSbac, ChannelType chType )
{
  Int firstCtx = 0, numCtx = NUM_SIG_CG_FLAG_CTX;

  for ( Int ctxIdx = firstCtx; ctxIdx < firstCtx + numCtx; ctxIdx++ )
  {
    for( UInt uiBin = 0; uiBin < 2; uiBin++ )
    {
      pcEstBitsSbac->significantCoeffGroupBits[ ctxIdx ][ uiBin ] = m_cCUSigCoeffGroupSCModel.get(  0, chType, ctxIdx ).getEntropyBits( uiBin );
    }
  }
}


/*!
 ****************************************************************************
 * \brief
 *    estimate SAMBAC bit cost for significant coefficient map
 ****************************************************************************
 */
Void TEncSbac::estSignificantMapBit( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType )
{
  //--------------------------------------------------------------------------------------------------

  //set up the number of channels and context variables

  const UInt firstComponent = ((isLuma(chType)) ? (COMPONENT_Y) : (COMPONENT_Cb));
  const UInt lastComponent  = ((isLuma(chType)) ? (COMPONENT_Y) : (COMPONENT_Cb));

  //----------------------------------------------------------
#if !VCEG_AZ07_CTX_RESIDUALCODING
  Int firstCtx = MAX_INT;
  Int numCtx   = MAX_INT;

  if      ((width == 4) && (height == 4))
  {
    firstCtx = significanceMapContextSetStart[chType][CONTEXT_TYPE_4x4];
    numCtx   = significanceMapContextSetSize [chType][CONTEXT_TYPE_4x4];
  }
  else if ((width == 8) && (height == 8))
  {
    firstCtx = significanceMapContextSetStart[chType][CONTEXT_TYPE_8x8];
    numCtx   = significanceMapContextSetSize [chType][CONTEXT_TYPE_8x8];
  }
  else
  {
    firstCtx = significanceMapContextSetStart[chType][CONTEXT_TYPE_NxN];
    numCtx   = significanceMapContextSetSize [chType][CONTEXT_TYPE_NxN];
  }
#endif
  //--------------------------------------------------------------------------------------------------

  //fill the data for the significace map

  for (UInt component = firstComponent; component <= lastComponent; component++)
  {
#if VCEG_AZ07_CTX_RESIDUALCODING
#if JVET_C0024_QTBT
    Int log2Size = ((g_aucConvertToBit[ width ] + g_aucConvertToBit[ height ]) >>1) + MIN_CU_LOG2;
    Int firstCtx   = !component  ? (log2Size>4? 2: (log2Size<=2 ? 0: log2Size-2)) * NUM_SIG_FLAG_CTX_LUMA_TU : 0;
#else
    Int firstCtx   = !component  ? (g_aucConvertToBit[ width ] > 2 ? 2 : g_aucConvertToBit[width]) * NUM_SIG_FLAG_CTX_LUMA_TU : 0;
#endif
    Int numCtx     = !component  ? NUM_SIG_FLAG_CTX_LUMA_TU : NUM_SIG_FLAG_CTX_CHROMA;
    Int iCtxOffset = !component  ? 0 : NUM_SIG_FLAG_CTX_LUMA;
#else
    const UInt contextOffset = getSignificanceMapContextOffset(ComponentID(component));

    if (firstCtx > 0)
    {
      for( UInt bin = 0; bin < 2; bin++ ) //always get the DC
      {
        pcEstBitsSbac->significantBits[ contextOffset ][ bin ] = m_cCUSigSCModel.get( 0, 0, contextOffset ).getEntropyBits( bin );
      }
    }
#endif
    // This could be made optional, but would require this function to have knowledge of whether the
    // TU is transform-skipped or transquant-bypassed and whether the SPS flag is set
#if !VCEG_AZ07_CTX_RESIDUALCODING
    for( UInt bin = 0; bin < 2; bin++ )
    {
      const Int ctxIdx = significanceMapContextSetStart[chType][CONTEXT_TYPE_SINGLE];
      pcEstBitsSbac->significantBits[ contextOffset + ctxIdx ][ bin ] = m_cCUSigSCModel.get( 0, 0, (contextOffset + ctxIdx) ).getEntropyBits( bin );
    }
#endif

    for ( Int ctxIdx = firstCtx; ctxIdx < firstCtx + numCtx; ctxIdx++ )
    {
      for( UInt uiBin = 0; uiBin < 2; uiBin++ )
      {
#if VCEG_AZ07_CTX_RESIDUALCODING
        pcEstBitsSbac->significantBits[ ctxIdx ][ uiBin ] = m_cCUSigSCModel.get(  0, 0, iCtxOffset + ctxIdx ).getEntropyBits( uiBin );
#else
        pcEstBitsSbac->significantBits[ contextOffset + ctxIdx ][ uiBin ] = m_cCUSigSCModel.get(  0, 0, (contextOffset + ctxIdx) ).getEntropyBits( uiBin );
#endif
      }
    }
  }

  //--------------------------------------------------------------------------------------------------
}


/*!
 ****************************************************************************
 * \brief
 *    estimate bit cost of significant coefficient
 ****************************************************************************
 */

Void TEncSbac::estLastSignificantPositionBit( estBitsSbacStruct* pcEstBitsSbac, Int width, Int height, ChannelType chType 
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
  , UInt uiScanIdx
#endif
  )
{
  //--------------------------------------------------------------------------------------------------.
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
  // swap
  if (uiScanIdx == SCAN_VER)
  {
    swap(width, height);
  }
#endif

  //set up the number of channels

  const UInt firstComponent = ((isLuma(chType)) ? (COMPONENT_Y) : (COMPONENT_Cb));
  const UInt lastComponent  = ((isLuma(chType)) ? (COMPONENT_Y) : (COMPONENT_Cb));

#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
  const UInt uiLog2BlockWidthIdx = g_aucConvertToBit[ width ];
  const UInt  ctxWidth   = isLuma(chType) ? width : 4;
  const UInt *puiCtxXIdx = g_uiLastCtx + (g_aucConvertToBit[ ctxWidth ] * (g_aucConvertToBit[ ctxWidth ] + 3 ) );   
#endif
  //--------------------------------------------------------------------------------------------------

  //fill the data for the last-significant-coefficient position

  for (UInt componentIndex = firstComponent; componentIndex <= lastComponent; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

    Int iBitsX = 0, iBitsY = 0;
#if !VCEG_AZ07_CTX_RESIDUALCODING || COM16_C806_T64
    Int blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY;
    getLastSignificantContextParameters(ComponentID(component), width, height, blkSizeOffsetX, blkSizeOffsetY, shiftX, shiftY);
#endif

    Int ctx;

    const ChannelType channelType = toChannelType(ComponentID(component));

    ContextModel *const pCtxX = m_cCuCtxLastX.get( 0, channelType );
    ContextModel *const pCtxY = m_cCuCtxLastY.get( 0, channelType );
    Int          *const lastXBitsArray = pcEstBitsSbac->lastXBits[channelType];
    Int          *const lastYBitsArray = pcEstBitsSbac->lastYBits[channelType];

    //------------------------------------------------

    //X-coordinate

    for (ctx = 0; ctx < g_uiGroupIdx[ width - 1 ]; ctx++)
    {
#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
      Int ctxOffset = isLuma(chType)? puiCtxXIdx[ ctx ]: (ctx >> uiLog2BlockWidthIdx );
#else
      Int ctxOffset = blkSizeOffsetX + (ctx >>shiftX);
#endif
      lastXBitsArray[ ctx ] = iBitsX + pCtxX[ ctxOffset ].getEntropyBits( 0 );
      iBitsX += pCtxX[ ctxOffset ].getEntropyBits( 1 );
    }

    lastXBitsArray[ctx] = iBitsX;

    //------------------------------------------------

    //Y-coordinate
#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
    const UInt uiLog2BlockheightIdx = g_aucConvertToBit[ height ];
    const UInt  ctxHeight  = isLuma(chType) ? height : 4;
    const UInt* puiCtxYIdx = g_uiLastCtx + (g_aucConvertToBit[ ctxHeight ] * (g_aucConvertToBit[ ctxHeight ] + 3 ) );   
#endif
    for (ctx = 0; ctx < g_uiGroupIdx[ height - 1 ]; ctx++)
    {
#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
      Int ctxOffset = isLuma(chType)? puiCtxYIdx[ ctx ]: (ctx >> uiLog2BlockheightIdx );
#else
      Int ctxOffset = blkSizeOffsetY + (ctx >>shiftY);
#endif
      lastYBitsArray[ ctx ] = iBitsY + pCtxY[ ctxOffset ].getEntropyBits( 0 );
      iBitsY += pCtxY[ ctxOffset ].getEntropyBits( 1 );
    }

    lastYBitsArray[ctx] = iBitsY;

  } //end of component loop

  //--------------------------------------------------------------------------------------------------
}


/*!
 ****************************************************************************
 * \brief
 *    estimate bit cost of significant coefficient
 ****************************************************************************
 */
Void TEncSbac::estSignificantCoefficientsBit( estBitsSbacStruct* pcEstBitsSbac, ChannelType chType )
{
#if VCEG_AZ07_CTX_RESIDUALCODING
  const UInt oneStopIndex  = ((isLuma(chType)) ? (NUM_ONE_FLAG_CTX_LUMA) : (NUM_ONE_FLAG_CTX_CHROMA));
  ContextModel *ctxOne = m_cCUOneSCModel.get(0, 0) + ((isLuma(chType)) ? (0) : (NUM_ONE_FLAG_CTX_LUMA));
  for (Int ctxIdx = 0; ctxIdx < oneStopIndex; ctxIdx++)
#else
  ContextModel *ctxOne = m_cCUOneSCModel.get(0, 0);
  ContextModel *ctxAbs = m_cCUAbsSCModel.get(0, 0);

  const UInt oneStartIndex = ((isLuma(chType)) ? (0)                     : (NUM_ONE_FLAG_CTX_LUMA));
  const UInt oneStopIndex  = ((isLuma(chType)) ? (NUM_ONE_FLAG_CTX_LUMA) : (NUM_ONE_FLAG_CTX));

  const UInt absStartIndex = ((isLuma(chType)) ? (0)                     : (NUM_ABS_FLAG_CTX_LUMA));
  const UInt absStopIndex  = ((isLuma(chType)) ? (NUM_ABS_FLAG_CTX_LUMA) : (NUM_ABS_FLAG_CTX));


  for (Int ctxIdx = oneStartIndex; ctxIdx < oneStopIndex; ctxIdx++)
#endif
  {
    pcEstBitsSbac->m_greaterOneBits[ ctxIdx ][ 0 ] = ctxOne[ ctxIdx ].getEntropyBits( 0 );
    pcEstBitsSbac->m_greaterOneBits[ ctxIdx ][ 1 ] = ctxOne[ ctxIdx ].getEntropyBits( 1 );
  }
#if !VCEG_AZ07_CTX_RESIDUALCODING
  for (Int ctxIdx = absStartIndex; ctxIdx < absStopIndex; ctxIdx++)
  {
    pcEstBitsSbac->m_levelAbsBits[ ctxIdx ][ 0 ] = ctxAbs[ ctxIdx ].getEntropyBits( 0 );
    pcEstBitsSbac->m_levelAbsBits[ ctxIdx ][ 1 ] = ctxAbs[ ctxIdx ].getEntropyBits( 1 );
  }
#endif
}

/**
 - Initialize our context information from the nominated source.
 .
 \param pSrc From where to copy context information.
 */
Void TEncSbac::xCopyContextsFrom( const TEncSbac* pSrc )
{
  memcpy(m_contextModels, pSrc->m_contextModels, m_numContextModels*sizeof(m_contextModels[0]));
  memcpy(m_golombRiceAdaptationStatistics, pSrc->m_golombRiceAdaptationStatistics, (sizeof(UInt) * RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS));
}

Void  TEncSbac::loadContexts ( const TEncSbac* pSrc)
{
  xCopyContextsFrom(pSrc);
}

#if VCEG_AZ07_INIT_PREVFRAME
Void  TEncSbac::loadContextsFromPrev (TComStats* apcStats, SliceType eSliceType, Int iQPIdx, Bool bFromGloble, Int iQPIdxRst, Bool bAfterLastISlice )
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
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ05_MULTI_PARAM_CABAC
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
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ05_MULTI_PARAM_CABAC
      apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i] = m_contextModels[i].getState();
#else
      apcStats->m_uiCtxProbIdx[eSliceType][iQPIdx][0][i] = (UShort)m_contextModels[i].getOrigState();
#endif
    }
  }
}
#endif
/** Performs CABAC encoding of the explicit RDPCM mode
 * \param rTu current TU data structure
 * \param compID component identifier
 */
Void TEncSbac::codeExplicitRdpcmMode( TComTU &rTu, const ComponentID compID )
{
  TComDataCU *cu = rTu.getCU();
#if !JVET_C0046_OMIT_ASSERT_ERDPCM
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt tuHeight = g_aucConvertToBit[rect.height];
  const UInt tuWidth  = g_aucConvertToBit[rect.width];
#endif
  const UInt absPartIdx   = rTu.GetAbsPartIdxTU(compID);

#if !JVET_C0046_OMIT_ASSERT_ERDPCM
  assert(tuHeight == tuWidth);
  assert(tuHeight < 4);
#endif

  UInt explicitRdpcmMode = cu->getExplicitRdpcmMode(compID, absPartIdx);

  if( explicitRdpcmMode == RDPCM_OFF )
  {
    m_pcBinIf->encodeBin (0, m_explicitRdpcmFlagSCModel.get (0, toChannelType(compID), 0));
  }
  else if( explicitRdpcmMode == RDPCM_HOR || explicitRdpcmMode == RDPCM_VER )
  {
    m_pcBinIf->encodeBin (1, m_explicitRdpcmFlagSCModel.get (0, toChannelType(compID), 0));
    if(explicitRdpcmMode == RDPCM_HOR)
    {
      m_pcBinIf->encodeBin ( 0, m_explicitRdpcmDirSCModel.get(0, toChannelType(compID), 0));
    }
    else
    {
      m_pcBinIf->encodeBin ( 1, m_explicitRdpcmDirSCModel.get(0, toChannelType(compID), 0));
    }
  }
  else
  {
    assert(0);
  }
}

#if ALF_HM3_REFACTOR
#if !JVET_C0024_QTBT
Void TEncSbac::codeAlfCtrlFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  if (!m_bAlfCtrl)
    return;

  if( pcCU->getDepth(uiAbsPartIdx) > m_uiMaxAlfCtrlDepth && !pcCU->isFirstAbsZorderIdxInDepth(uiAbsPartIdx, m_uiMaxAlfCtrlDepth))
  {
    return;
  }

  // get context function is here
  UInt uiSymbol = pcCU->getAlfCtrlFlag( uiAbsPartIdx ) ? 1 : 0;

  m_pcBinIf->encodeBin( uiSymbol, m_cCUAlfCtrlFlagSCModel.get( 0, 0, pcCU->getCtxAlfCtrlFlag( uiAbsPartIdx) ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlFlag\n" )
}
#endif

Void TEncSbac::codeAlfCtrlDepth( UInt uiMaxTotalCUDepth )
{
  UInt uiDepth = m_uiMaxAlfCtrlDepth;
  xWriteUnaryMaxSymbol(uiDepth, m_cALFUvlcSCModel.get(0), 1, uiMaxTotalCUDepth-1);
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlDepth\n" )
}
#if JVET_C0038_GALF  
//0: no pred; 1: all same index; 2: diff index for each variance index
Void TEncSbac::codeALFPrevFiltType( UInt uiCode)
{
  xWriteEpExGolomb( uiCode, 0);
}
Void TEncSbac::codeALFPrevFiltFlag( Int uiCode)
{
   m_pcBinIf->encodeBinEP( uiCode > 0 ? 1: 0 );
}
#endif
Void TEncSbac::codeAlfFlag       ( UInt uiCode )
{
  UInt uiSymbol = ( ( uiCode == 0 ) ? 0 : 1 );
#if JVET_C0038_GALF
  m_pcBinIf->encodeBinEP( uiSymbol);
#else
  m_pcBinIf->encodeBin( uiSymbol, m_cALFFlagSCModel.get( 0, 0, 0 ) );
#endif
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfFlag\n" )
}

Void TEncSbac::codeAlfFlagNum( UInt uiCode, UInt minValue )
{
  UInt uiLength = 0;
  UInt maxValue = (minValue << (this->getMaxAlfCtrlDepth()*2));
  assert((uiCode>=minValue)&&(uiCode<=maxValue));
  UInt temp = maxValue - minValue;
  for(UInt i=0; i<32; i++)
  {
    if(temp&0x1)
    {
      uiLength = i+1;
    }
    temp = (temp >> 1);
  }
  UInt uiSymbol = uiCode - minValue;
  if(uiLength)
  {
    while( uiLength-- )
    {
      m_pcBinIf->encodeBinEP( (uiSymbol>>uiLength) & 0x1 );
    }
  }
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfFlagNum\n" )
}

Void TEncSbac::codeAlfCtrlFlag( UInt uiSymbol )
{
  m_pcBinIf->encodeBin( uiSymbol, m_cCUAlfCtrlFlagSCModel.get( 0, 0, 0) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ )
  DTRACE_CABAC_T( "\tAlfCtrlFlag\n" )
}

Void TEncSbac::codeAlfUvlc       ( UInt uiCode )
{
  Int i;
#if JVET_C0038_GALF
  if ( uiCode == 0 )
  {
    m_pcBinIf->encodeBinEP(0);
  }
  else
  {
    m_pcBinIf->encodeBinEP(1);
    for ( i=0; i<uiCode-1; i++ )
    {
        m_pcBinIf->encodeBinEP(1);
    }
    m_pcBinIf->encodeBinEP(0);
  }
#else
  if ( uiCode == 0 )
  {
    m_pcBinIf->encodeBin( 0, m_cALFUvlcSCModel.get( 0, 0, 0 ) );
  }
  else
  {
    m_pcBinIf->encodeBin( 1, m_cALFUvlcSCModel.get( 0, 0, 0 ) );
    for ( i=0; i<uiCode-1; i++ )
    {
      m_pcBinIf->encodeBin( 1, m_cALFUvlcSCModel.get( 0, 0, 1 ) );
    }
    m_pcBinIf->encodeBin( 0, m_cALFUvlcSCModel.get( 0, 0, 1 ) );
  }
#endif
}

Void TEncSbac::codeAlfSvlc       ( Int iCode )
{
  Int i;
#if JVET_C0038_GALF
  if ( iCode == 0 )
  {
    m_pcBinIf->encodeBinEP(0);
  }
  else
  {
    m_pcBinIf->encodeBinEP(1);

    // write sign
    if ( iCode > 0 )
    {
      m_pcBinIf->encodeBinEP(0);
    }
    else
    {
     m_pcBinIf->encodeBinEP(1);
      iCode = -iCode;
    }

    // write magnitude
    for ( i=0; i<iCode-1; i++ )
    {
      m_pcBinIf->encodeBinEP(1);
    }
    m_pcBinIf->encodeBinEP(0);
  }
#else
  if ( iCode == 0 )
  {
    m_pcBinIf->encodeBin( 0, m_cALFSvlcSCModel.get( 0, 0, 0 ) );
  }
  else
  {
    m_pcBinIf->encodeBin( 1, m_cALFSvlcSCModel.get( 0, 0, 0 ) );

    // write sign
    if ( iCode > 0 )
    {
      m_pcBinIf->encodeBin( 0, m_cALFSvlcSCModel.get( 0, 0, 1 ) );
    }
    else
    {
      m_pcBinIf->encodeBin( 1, m_cALFSvlcSCModel.get( 0, 0, 1 ) );
      iCode = -iCode;
    }

    // write magnitude
    for ( i=0; i<iCode-1; i++ )
    {
      m_pcBinIf->encodeBin( 1, m_cALFSvlcSCModel.get( 0, 0, 2 ) );
    }
    m_pcBinIf->encodeBin( 0, m_cALFSvlcSCModel.get( 0, 0, 2 ) );
  }
#endif
}
#endif

#if COM16_C806_EMT
#if FAST_EMT_v1
Void TEncSbac::codeEmtTuIdx(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
	if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTRA_MAX_CU)
#else
	if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU)
#endif
	{


		UChar ucTrIdx = pcCU->getEmtTuIdx(uiAbsPartIdx);

		m_pcBinIf->encodeBin(ucTrIdx ? 1 : 0, m_cEMTIntraidxSCModel.get(0, 0, 1));
		if (ucTrIdx)
		{
			if (ucTrIdx == 1) m_pcBinIf->encodeBin(0, m_cEMTIntraidxSCModel.get(0, 0, 3));
			else m_pcBinIf->encodeBin(1, m_cEMTIntraidxSCModel.get(0, 0, 3));
		}

		//	m_pcBinIf->encodeBin( ( ucTrIdx & 1 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 0));
		//		m_pcBinIf->encodeBin( ( ucTrIdx / 2 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 1));
	}
#if JVET_C0024_QTBT
	if (!pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx) <= EMT_INTER_MAX_CU)
#else
	if (!pcCU->isIntra(uiAbsPartIdx) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU)
#endif
	{
		UChar ucTrIdx = pcCU->getEmtTuIdx(uiAbsPartIdx);
#if FAST_EMT_v1 //EMT_SIGNAL_INTER
		m_pcBinIf->encodeBin(ucTrIdx ? 1 : 0, m_cEMTIntraidxSCModel.get(0, 0, 0));
		if (ucTrIdx)
		{
			if (ucTrIdx == 1) m_pcBinIf->encodeBin(0, m_cEMTIntraidxSCModel.get(0, 0, 2));
			else m_pcBinIf->encodeBin(1, m_cEMTIntraidxSCModel.get(0, 0, 2));
		}
#else
		m_pcBinIf->encodeBin((ucTrIdx & 1) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 2));
		m_pcBinIf->encodeBin((ucTrIdx / 2) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 3));
#endif
	}
}
#else
Void TEncSbac::codeEmtTuIdx( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  if ( pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx)<= EMT_INTRA_MAX_CU )
#else
  if ( pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU )
#endif
  {
    UChar ucTrIdx = pcCU->getEmtTuIdx( uiAbsPartIdx );
    m_pcBinIf->encodeBin( ( ucTrIdx & 1 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 0));
    m_pcBinIf->encodeBin( ( ucTrIdx / 2 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 1));
  }
#if JVET_C0024_QTBT
  if ( !pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <=EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx)<= EMT_INTER_MAX_CU  )
#else
  if ( !pcCU->isIntra( uiAbsPartIdx ) && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU )
#endif
  {
    UChar ucTrIdx = pcCU->getEmtTuIdx( uiAbsPartIdx );
    m_pcBinIf->encodeBin( ( ucTrIdx & 1 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 2));
    m_pcBinIf->encodeBin( ( ucTrIdx / 2 ) ? 1 : 0, m_cEmtTuIdxSCModel.get(0, 0, 3));
  }
}
#endif
Void TEncSbac::codeEmtCuFlag( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool bCodeCuFlag )
{
#if JVET_C0024_QTBT
  if (uiDepth >= NUM_EMT_CU_FLAG_CTX)
  {
    uiDepth = NUM_EMT_CU_FLAG_CTX - 1;
  }
#else
  assert( uiDepth < NUM_EMT_CU_FLAG_CTX );
#endif
  
#if JVET_C0024_QTBT
  if ( pcCU->isIntra( uiAbsPartIdx ) && bCodeCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getHeight(uiAbsPartIdx)<= EMT_INTRA_MAX_CU
      && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
#else
  if ( pcCU->isIntra( uiAbsPartIdx ) && bCodeCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTRA_MAX_CU && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
#endif
  {
    UChar ucCuFlag = pcCU->getEmtCuFlag( uiAbsPartIdx );
    m_pcBinIf->encodeBin( ucCuFlag, m_cEmtCuFlagSCModel.get(0, 0, uiDepth));
  }
#if JVET_C0024_QTBT
  if( !pcCU->isIntra( uiAbsPartIdx ) && bCodeCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU && pcCU->getHeight(uiAbsPartIdx)<= EMT_INTER_MAX_CU
      && pcCU->getSlice()->getSPS()->getUseInterEMT() )
#else
  if( !pcCU->isIntra( uiAbsPartIdx ) && bCodeCuFlag && pcCU->getWidth(uiAbsPartIdx) <= EMT_INTER_MAX_CU && pcCU->getSlice()->getSPS()->getUseInterEMT() )
#endif
  {
    UChar ucCuFlag = pcCU->getEmtCuFlag( uiAbsPartIdx );
    m_pcBinIf->encodeBin( ucCuFlag, m_cEmtCuFlagSCModel.get(0, 0, uiDepth));
  }
}
#endif

#if COM16_C1016_AFFINE
/** code affine flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \returns Void
 */
Void TEncSbac::codeAffineFlag( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  // get context function is here
  UInt uiSymbol = pcCU->isAffine( uiAbsPartIdx ) ? 1 : 0;
  UInt uiCtxAffine = pcCU->getCtxAffineFlag( uiAbsPartIdx );
  m_pcBinIf->encodeBin( uiSymbol, m_cCUAffineFlagSCModel.get( 0, 0, uiCtxAffine ) );

  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tAffineFlag" );
  DTRACE_CABAC_T( "\tuiCtxAffine: ");
  DTRACE_CABAC_V( uiCtxAffine );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}
#endif

#if MRG_CONSIDER_COMPLEX
Void TEncSbac::codeMergeComplexFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  UInt uiCtx = 0;
  uiCtx = pcCU->getCtxMergeComplexFlag(uiAbsPartIdx);
  const UInt uiSymbol = pcCU->getMergeComplexFlag(uiAbsPartIdx) ? 1 : 0;
  m_pcBinIf->encodeBin(uiSymbol, m_cCUMergeComplexFlagExtSCModel.get(0, 0, uiCtx));
  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tMergeComplexFlag: ");
  DTRACE_CABAC_V(uiSymbol);
  DTRACE_CABAC_T("\n");
}

Void TEncSbac::codeMergeComplexIndex(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  UInt uiUnaryIdx = pcCU->getMergeComplexIndex(uiAbsPartIdx);
  UChar numCand = MRG_COMPLEX_NUM_CANDS_REAL;
  Int n = 0;
  UChar bin[4];
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

  while (uiUnaryIdx > 0)
  {
    bin[n++] = uiUnaryIdx % 2;
    uiUnaryIdx /= 2;
  }

  while (n < numBits)
  {
    bin[n++] = 0;
  }

  UInt uiCtx = 0;
  UInt dec_min_val = 0;
  UInt max_val = numCand - 1;

  n = n - 1; // Point to the MSB(Most significant bit) , it is like bin_idx
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
      const UInt uiSymbol = bin[n];

      if (n == numBits - 1)
        m_pcBinIf->encodeBin(uiSymbol, m_cCUMergeComplexIdxExtSCModel.get(0, 0, 0));
      else
        m_pcBinIf->encodeBinEP(uiSymbol);

      dec_min_val = dec_min_val | (uiSymbol << n);
    }
    n--;
  }

  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tparseMergeComplexIndex()");
  DTRACE_CABAC_T("\tuiMRGComplexIdx= ");
  DTRACE_CABAC_V(pcCU->getMergeComplexIndex(uiAbsPartIdx));
  DTRACE_CABAC_T("\n");
}
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
Void TEncSbac::codeAffineParamFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  // get context function is here
  UInt uiSymbol = pcCU->getAffineParamFlag(uiAbsPartIdx) ? 1 : 0;
  assert(pcCU->isAffine(uiAbsPartIdx));


  UInt uiCtxAffine = pcCU->getCtxAffineParamFlag(uiAbsPartIdx);
  m_pcBinIf->encodeBin(uiSymbol, m_cCUAffineParamFlagSCModel.get(0, 0, uiCtxAffine));


  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tAffineParamFlag");
  DTRACE_CABAC_T("\tuiCtxAffineParam: ");
  DTRACE_CABAC_V(uiCtxAffine);
  DTRACE_CABAC_T("\tuiSymbol: ");
  DTRACE_CABAC_V(uiSymbol);
  DTRACE_CABAC_T("\n");
}
#endif

#if VCEG_AZ07_BAC_ADAPT_WDOW
Void TEncSbac::xUpdateWindowSize ( SliceType eSliceType, Int uiQPIdx, TComStats* apcStats )
{
  if( uiQPIdx == -1 )
  {
    return;
  }
  Int iCtxNr = getCtxNumber();
  for(UInt i=0; i<iCtxNr; i++)
  {
    m_contextModels[i].setWindowSize(apcStats->m_uiCtxCodeIdx[eSliceType][uiQPIdx][i]);
  }
}
Void TEncSbac::codeCtxUpdateInfo  ( TComSlice* pcSlice,  TComStats* apcStats )
{
  assert(0);
  return;
}
#endif

#if RRC_NS
Void TEncSbac::codeNsBlockSize(UInt uiCode)
{
	assert(uiCode < 4);
	m_pcBinIf->encodeBinEP(uiCode >> 1);
	m_pcBinIf->encodeBinEP(uiCode & 1);
}

Void TEncSbac::codeNsFlag(UInt uiCode)
{
	m_pcBinIf->encodeBin(!!uiCode, m_cNsFlagSCModel.get(0, 0, 0));
}

Void TEncSbac::codeAdaptiveNsFlag(UInt uiCode)
{
	m_pcBinIf->encodeBinEP(uiCode);
}

Void TEncSbac::codeNsFWholeFrameFlag(UInt uiCode)
{
	m_pcBinIf->encodeBinEP(uiCode);
}

#endif

#if RRC_INTRA_SHARPENING
Void TEncSbac::codeIntraSharpFlag(UInt uiCode)
{
  m_pcBinIf->encodeBin(uiCode, m_cSharpFlagSCModel.get(0, 0, 0));

  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tSharp flag: ");
  DTRACE_CABAC_V(symbol);
  DTRACE_CABAC_T("\n");
}
#endif

#if RRC_SELECTABLE_IF
Void TEncSbac::codeAdaptiveIfIndex(Int idx)
{
  m_pcBinIf->encodeBin(idx, m_cCUAdaptiveIfIndex.get(0, 0, 0));
  DTRACE_CABAC_VL(g_nSymbolCounter++);
  DTRACE_CABAC_T("\tSIF index: ");
  DTRACE_CABAC_V(symbol);
  DTRACE_CABAC_T("\n");
}
#endif

#if BARC_UP_SAMPLE_RDO
Void TEncSbac::codeUpCNNFlag(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
	// get context function is here
	UInt uiSymbol = pcCU->getUpCNNFlag() ? 1 : 0;
	UInt uiCtxNetRDO = pcCU->getCtxUpCNNFlag(uiAbsPartIdx);
	m_pcBinIf->encodeBin(uiSymbol, m_cCUUpCNNFlagSCModel.get(0, 0, uiCtxNetRDO));
}
#endif

#if BARC
/** code BARC flag
 * \param pcCU
 * \param uiAbsPartIdx 
 * \returns Void
 */
Void TEncSbac::codeBARCFlag( TComDataCU* pcCU, ChannelType eCtype )
{
  // get context function is here
  UInt uiSymbol = pcCU->getBARCFlag(eCtype) ? 1 : 0;
  UInt uiCtxBARC = pcCU->getCtxBARCFlag( 0 ) ;
  m_pcBinIf->encodeBin( uiSymbol, m_cCUBARCFlagSCModel.get( 0, 0, uiCtxBARC ) );
  DTRACE_CABAC_VL( g_nSymbolCounter++ );
  DTRACE_CABAC_T( "\tBARCFlag" );
  DTRACE_CABAC_T( "\tuiCtxBARC: ");
  DTRACE_CABAC_V( uiCtxBARC );
  DTRACE_CABAC_T( "\tuiSymbol: ");
  DTRACE_CABAC_V( uiSymbol );
  DTRACE_CABAC_T( "\n");
}
#endif
//! \}
