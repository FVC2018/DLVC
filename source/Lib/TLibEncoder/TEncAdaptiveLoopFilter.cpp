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

/** \file     TEncAdaptiveLoopFilter.cpp
 \brief    estimation part of adaptive loop filter class
 */
#include "TEncAdaptiveLoopFilter.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// ====================================================================================================================
// Constants
// ====================================================================================================================

#if ALF_HM3_REFACTOR

#define ALF_NUM_OF_REDESIGN 1

// ====================================================================================================================
// Tables
// ====================================================================================================================
const Int TEncAdaptiveLoopFilter::m_aiSymmetricArray9x9[81] =
{
   0,  1,  2,  3,  4,  5,  6,  7,  8,
   9, 10, 11, 12, 13, 14, 15, 16, 17,
  18, 19, 20, 21, 22, 23, 24, 25, 26,
  27, 28, 29, 30, 31, 32, 33, 34, 35,
  36, 37, 38, 39, 40, 39, 38, 37, 36,
  35, 34, 33, 32, 31, 30, 29, 28, 27,
  26, 25, 24, 23, 22, 21, 20, 19, 18,
  17, 16, 15, 14, 13, 12, 11, 10,  9,
   8,  7,  6,  5,  4,  3,  2,  1,  0
};

const Int TEncAdaptiveLoopFilter::m_aiSymmetricArray7x7[49] =
{
  0,  1,  2,  3,  4,  5,  6,
  7,  8,  9, 10, 11, 12, 13,
  14, 15, 16, 17, 18, 19, 20,
  21, 22, 23, 24, 23, 22, 21,
  20, 19, 18, 17, 16, 15, 14,
  13, 12, 11, 10,  9,  8,  7,
  6,  5,  4,  3,  2,  1,  0,
};
#if JVET_C0038_GALF 
const Int TEncAdaptiveLoopFilter::m_aiSymmetricArray5x5[25] =
{
  0,  0,  1,  0,  0,
  0,  2,  3,  4,  0,
  5,  6,  7,  6,  5,
  0,  4,  3,  2,  0,
  0,  0,  1,  0,  0,
};
#else
const Int TEncAdaptiveLoopFilter::m_aiSymmetricArray5x5[25] =
{
  0,  1,  2,  3,  4,
  5,  6,  7,  8,  9,
  10, 11, 12, 11, 10,
  9,  8,  7,  6,  5,
  4,  3,  2,  1,  0,
};
#endif
const Int TEncAdaptiveLoopFilter::m_aiSymmetricArray9x7[63] =
{
   0,  1,  2,  3,  4,  5,  6,  7,  8,
   9, 10, 11, 12, 13, 14, 15, 16, 17,
  18, 19, 20, 21, 22, 23, 24, 25, 26,
  27, 28, 29, 30, 31, 30, 29, 28, 27,
  26, 25, 24, 23, 22, 21, 20, 19, 18,
  17, 16, 15, 14, 13, 12, 11, 10,  9,
   8,  7,  6,  5,  4,  3,  2,  1,  0
};

Int TEncAdaptiveLoopFilter::m_aiTapPos9x9_In9x9Sym[21] =
{
                  0,  1,  2,
              3,  4,  5,  6,  7,
          8,  9, 10, 11, 12, 13, 14,
     15, 16, 17, 18, 19, 20
};

Int TEncAdaptiveLoopFilter::m_aiTapPos7x7_In9x9Sym[14] =
{                 
                  1,    
              4,  5,  6,    
          9, 10, 11, 12, 13,    
     16, 17, 18, 19, 20
};

Int TEncAdaptiveLoopFilter::m_aiTapPos5x5_In9x9Sym[8]  =
{

            5,
       10, 11, 12, 
   17, 18, 19, 20
};

Int* TEncAdaptiveLoopFilter::m_iTapPosTabIn9x9Sym[m_NO_TEST_FILT] =
{
  m_aiTapPos9x9_In9x9Sym, m_aiTapPos7x7_In9x9Sym, m_aiTapPos5x5_In9x9Sym
};

// ====================================================================================================================
// Constructor / destructor
// ====================================================================================================================

TEncAdaptiveLoopFilter::TEncAdaptiveLoopFilter()
{
  m_ppdAlfCorr = NULL;
  m_pdDoubleAlfCoeff = NULL;
  m_pcPic = NULL;
  m_pcEntropyCoder = NULL;
  m_pcBestAlfParam = NULL;
  m_pcTempAlfParam = NULL;
  m_pcPicYuvBest = NULL;
  m_pcPicYuvTmp = NULL;
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/**
 \param pcPic           picture (TComPic) pointer
 \param pcEntropyCoder  entropy coder class
 */
Void TEncAdaptiveLoopFilter::startALFEnc( TComPic* pcPic, TEncEntropy* pcEntropyCoder )
{
  m_pcPic = pcPic;
  m_pcEntropyCoder = pcEntropyCoder;
  
  m_eSliceType = pcPic->getSlice(0)->getSliceType();
  m_iPicNalReferenceIdc = (pcPic->getSlice(0)->isReferenced() ? 1 :0);
  
  m_uiNumSCUInCU = m_pcPic->getNumPartitionsInCtu();
 
  xInitParam();
  xCreateTmpAlfCtrlFlags();
  
  Int iWidth = pcPic->getPicYuvOrg()->getWidth(COMPONENT_Y);
  Int iHeight = pcPic->getPicYuvOrg()->getHeight(COMPONENT_Y);
  
  m_pcPicYuvTmp = new TComPicYuv();
#if JVET_C0024_QTBT
  m_pcPicYuvTmp->create(iWidth, iHeight, pcPic->getChromaFormat(), pcPic->getSlice(0)->getSPS()->getCTUSize(), pcPic->getSlice(0)->getSPS()->getCTUSize(), pcPic->getSlice(0)->getSPS()->getMaxTotalCUDepth(), true);
#else
  m_pcPicYuvTmp->create(iWidth, iHeight, pcPic->getChromaFormat(), pcPic->getSlice(0)->getSPS()->getMaxCUWidth(), pcPic->getSlice(0)->getSPS()->getMaxCUHeight(), pcPic->getSlice(0)->getSPS()->getMaxTotalCUDepth(), true);
#endif
  m_pcPicYuvBest = pcPic->getPicYuvPred();
  
  m_pcBestAlfParam = new ALFParam;
  m_pcTempAlfParam = new ALFParam;
  allocALFParam(m_pcBestAlfParam);
  allocALFParam(m_pcTempAlfParam);

  m_im_width = iWidth;
  m_im_height = iHeight;
  
  // init qc_filter
  initMatrix4D_double(&m_EGlobalSym, m_NO_TEST_FILT,  m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);
  initMatrix3D_double(&m_yGlobalSym, m_NO_TEST_FILT, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH); 
  initMatrix_int(&m_filterCoeffSymQuant, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH); 
  
  m_pixAcc = (double *) calloc(m_NO_VAR_BINS, sizeof(double));
  get_mem2Dpel(&m_maskImg, m_im_height, m_im_width);
#if JVET_C0038_GALF 
  get_mem2Dpel(&m_maskBestpImg, m_im_height, m_im_width);
#endif

  initMatrix_double(&m_E_temp, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);//
  m_y_temp = (double *) calloc(m_MAX_SQR_FILT_LENGTH, sizeof(double));//
  initMatrix3D_double(&m_E_merged, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);//
  initMatrix_double(&m_y_merged, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH); //
  m_pixAcc_merged = (double *) calloc(m_NO_VAR_BINS, sizeof(double));//
  
  m_filterCoeffQuantMod = (int *) calloc(m_MAX_SQR_FILT_LENGTH, sizeof(int));//
  m_filterCoeff = (double *) calloc(m_MAX_SQR_FILT_LENGTH, sizeof(double));//
  m_filterCoeffQuant = (int *) calloc(m_MAX_SQR_FILT_LENGTH, sizeof(int));//
  initMatrix_int(&m_diffFilterCoeffQuant, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);//
  initMatrix_int(&m_FilterCoeffQuantTemp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);//
  
  m_tempALFp = new ALFParam;
  allocALFParam(m_tempALFp);
  m_pcDummyEntropyCoder = m_pcEntropyCoder;
#if JVET_C0038_GALF
  initMatrix_int(&m_imgY_preFilter, iHeight, iWidth);
#endif
}

Void TEncAdaptiveLoopFilter::endALFEnc()
{
  xUninitParam();
  xDestroyTmpAlfCtrlFlags();
  
  m_pcPicYuvTmp->destroy();
  delete m_pcPicYuvTmp;
  m_pcPicYuvTmp = NULL;
  m_pcPic = NULL;
  m_pcEntropyCoder = NULL;
  
  freeALFParam(m_pcBestAlfParam);
  freeALFParam(m_pcTempAlfParam);
  delete m_pcBestAlfParam;
  delete m_pcTempAlfParam;

  // delete qc filters
  destroyMatrix4D_double(m_EGlobalSym, m_NO_TEST_FILT,  m_NO_VAR_BINS);
  destroyMatrix3D_double(m_yGlobalSym, m_NO_TEST_FILT);
  destroyMatrix_int(m_filterCoeffSymQuant);
  
  free(m_pixAcc);
  free_mem2Dpel(m_maskImg);
#if JVET_C0038_GALF
  free_mem2Dpel(m_maskBestpImg);
#endif 
  destroyMatrix3D_double(m_E_merged, m_NO_VAR_BINS);
  destroyMatrix_double(m_y_merged);
  destroyMatrix_double(m_E_temp);
  free(m_pixAcc_merged);
  
  free(m_filterCoeffQuantMod);
  free(m_y_temp);
  
  free(m_filterCoeff);
  free(m_filterCoeffQuant);
  destroyMatrix_int(m_diffFilterCoeffQuant);
  destroyMatrix_int(m_FilterCoeffQuantTemp);
#if JVET_C0038_GALF
  destroyMatrix_int(m_imgY_preFilter);
#endif

  freeALFParam(m_tempALFp);
  delete m_tempALFp;
}

/**
 \param pcAlfParam           ALF parameter
 \param dLambda              lambda value for RD cost computation
 \retval ruiDist             distortion
 \retval ruiBits             required bits
 \retval ruiMaxAlfCtrlDepth  optimal partition depth
 */
Void TEncAdaptiveLoopFilter::ALFProcess( ALFParam* pcAlfParam, Double dLambdaLuma , Double dLambdaChroma, UInt64& ruiDist, UInt64& ruiBits, UInt& ruiMaxAlfCtrlDepth
#if COM16_C806_ALF_TEMPPRED_NUM
  , ALFParam *pcStoredAlfPara, Int iStoredAlfParaNum
#endif
  )
{
  Int tap, num_coef;
  
#if JVET_C0024_QTBT
  m_pcPic->getSlice(0)->setTextType(CHANNEL_TYPE_LUMA);  //cu level on off only for luma
#endif
  // set global variables
  tap         = m_ALF_MAX_NUM_TAP;
  Int tapV = TComAdaptiveLoopFilter::ALFTapHToTapV(tap);
  num_coef = (tap * tapV + 1) >> 1;
#if !JVET_C0038_GALF
  num_coef    = num_coef + 1; // DC offset
#endif

  resetALFParam( m_pcBestAlfParam );
  resetALFParam( m_pcTempAlfParam );

#if JVET_C0038_GALF  
  xInitFixedFilters();
#endif
  // set lambda
  m_dLambdaLuma   = dLambdaLuma;
  m_dLambdaChroma = dLambdaChroma;
  
  TComPicYuv* pcPicOrg = m_pcPic->getPicYuvOrg();
  
  // extend image for filtering
  TComPicYuv* pcPicYuvRec    = m_pcPic->getPicYuvRec();
  TComPicYuv* pcPicYuvExtRec = m_pcTempPicYuv;
  
  pcPicYuvRec->copyToPic(pcPicYuvExtRec);
  pcPicYuvExtRec->setBorderExtension( false );
#if JVET_C0038_GALF
  pcPicYuvExtRec->extendPicBorder   (m_PADDING_W_ALF);
  pcPicYuvRec->copyToPic(m_pcPicYuvTmp, COMPONENT_Y);
#else
  pcPicYuvExtRec->extendPicBorder   (m_FILTER_LENGTH >> 1);
#endif 
  // set min cost
  UInt64 uiMinRate = MAX_INT;
  UInt64 uiMinDist = MAX_INT;
  Double dMinCost  = MAX_DOUBLE;

  UInt64  uiOrigRate;
  UInt64  uiOrigDist;
  Double  dOrigCost;
 
  // calc original cost
  xCalcRDCost( pcPicOrg, pcPicYuvRec, NULL, uiOrigRate, uiOrigDist, dOrigCost , m_pcPic->getSlice(0));
  m_pcBestAlfParam->alf_flag = 0;
  m_pcBestAlfParam->cu_control_flag = 0;

#if JVET_C0038_GALF  
  TComAdaptiveLoopFilter::m_max_NO_VAR_BINS = m_NO_VAR_BINS;
  TComAdaptiveLoopFilter::m_max_NO_FILTERS  = m_NO_FILTERS ;    
#endif
  // initialize temp_alfps
  m_pcTempAlfParam->alf_flag        = 1;
  m_pcTempAlfParam->tap             = tap;
  m_pcTempAlfParam->tapV            = tapV;
  m_pcTempAlfParam->num_coeff       = num_coef;
  m_pcTempAlfParam->chroma_idc      = 0;
  m_pcTempAlfParam->cu_control_flag = 0;
  
  // adaptive in-loop wiener filtering
  xEncALFLuma_qc( pcPicOrg, pcPicYuvExtRec, pcPicYuvRec, uiMinRate, uiMinDist, dMinCost , m_pcPic->getSlice(0) );
  
  // cu-based filter on/off control
#if JVET_C0038_GALF
  if(!m_pcPic->getSlice(0)->isIntra())
#endif
  xCUAdaptiveControl_qc( pcPicOrg, pcPicYuvExtRec, pcPicYuvRec, uiMinRate, uiMinDist, dMinCost , m_pcPic->getSlice(0) ); 


  // adaptive tap-length
  xFilterTapDecision_qc( pcPicOrg, pcPicYuvExtRec, pcPicYuvRec, uiMinRate, uiMinDist, dMinCost, m_pcPic->getSlice(0));
  
  // compute RD cost
  xCalcRDCost( pcPicOrg, pcPicYuvRec, m_pcBestAlfParam, uiMinRate, uiMinDist, dMinCost, m_pcPic->getSlice(0) );

  // compare RD cost to non-ALF case
  if( dMinCost < dOrigCost )
  {
    m_pcBestAlfParam->alf_flag = 1;
    
    ruiBits = uiMinRate;
    ruiDist = uiMinDist;
  } 
  else
  {
    m_pcBestAlfParam->alf_flag        = 0;
    m_pcBestAlfParam->cu_control_flag = 0;
    
    uiMinRate = uiOrigRate;
    uiMinDist = uiOrigDist;
    dMinCost = dOrigCost;
    m_pcEntropyCoder->setAlfCtrl(false);
    pcPicYuvExtRec->copyToPic(pcPicYuvRec, COMPONENT_Y);
    
    ruiBits = uiOrigRate;
    ruiDist = uiOrigDist;
  }
  
  // if ALF works
  if( m_pcBestAlfParam->alf_flag )
  {
#if !JVET_C0038_GALF
    // predict ALF coefficients
    predictALFCoeff( m_pcBestAlfParam );
#endif
    // do additional ALF process for chroma
    xEncALFChroma( uiMinRate, pcPicOrg, pcPicYuvExtRec, pcPicYuvRec, ruiDist, ruiBits , m_pcPic->getSlice(0));
  }
  
  // copy to best storage
  copyALFParam(pcAlfParam, m_pcBestAlfParam);
  
  // store best depth
  ruiMaxAlfCtrlDepth = m_pcEntropyCoder->getMaxAlfCtrlDepth();
  
#if COM16_C806_ALF_TEMPPRED_NUM
  if( pcStoredAlfPara != NULL && iStoredAlfParaNum > 0 )
  {
    pcPicYuvRec->copyToPic( m_pcPicYuvBest );
    //test stored ALF coefficients of both luma and chroma.
    for( UInt iAlfIdx = 0; iAlfIdx < iStoredAlfParaNum && iAlfIdx < COM16_C806_ALF_TEMPPRED_NUM; iAlfIdx++ )
    {
      pcPicYuvExtRec->copyToPic( pcPicYuvRec );
      ALFParam pcStoredAlf = pcStoredAlfPara[iAlfIdx];

      if( xFilteringLumaChroma( &pcStoredAlf, pcPicOrg, pcPicYuvExtRec, pcPicYuvRec, ruiBits, ruiDist, dMinCost, iAlfIdx , m_pcPic->getSlice( 0 )) )
      {
        pcAlfParam->temproalPredFlag = false;
        copyALFParam(pcAlfParam, m_pcBestAlfParam);
        pcAlfParam->temproalPredFlag = true;
        pcAlfParam->prevIdx = iAlfIdx;
        if( m_pcBestAlfParam->cu_control_flag )
        {
          pcAlfParam->alf_max_depth = m_pcBestAlfParam->alf_max_depth;
        }
      }
    }
    m_pcPicYuvBest->copyToPic( pcPicYuvRec );
  }
#endif
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

Void TEncAdaptiveLoopFilter::xEncALFChroma( UInt64 uiLumaRate, TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiDist, UInt64& ruiBits, const TComSlice *pSlice )
{
  // restriction for non-referenced B-slice
  if (m_eSliceType == B_SLICE && m_iPicNalReferenceIdc == 0)
  {
    return;
  }
  
  Int tap, num_coef;
  
  // set global variables
  tap         = m_ALF_MAX_NUM_TAP_C;
#if JVET_C0038_GALF 
  num_coef    = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(tap);
#else
  num_coef    = (tap*tap+1)>>1;
  num_coef    = num_coef + 1; // DC offset
#endif

  // set min cost
  UInt64 uiMinRate = uiLumaRate;
  UInt64 uiMinDist = MAX_INT;
  Double dMinCost  = MAX_DOUBLE;

  // calc original cost
  copyALFParam(m_pcTempAlfParam, m_pcBestAlfParam);
  xCalcRDCostChroma(pcPicOrg, pcPicRest, m_pcTempAlfParam, uiMinRate, uiMinDist, dMinCost, pSlice);
  
  // initialize temp_alfps
  m_pcTempAlfParam->chroma_idc = 3;
  m_pcTempAlfParam->tap_chroma       = tap;
  m_pcTempAlfParam->num_coeff_chroma = num_coef;
  
  // Adaptive in-loop wiener filtering for chroma
  xFilteringFrameChroma(pcPicOrg, pcPicDec, pcPicRest);
  
  // filter on/off decision for chroma
  Int iCWidth = pcPicOrg->getWidth(COMPONENT_Cb);
  Int iCHeight = pcPicOrg->getHeight(COMPONENT_Cb);
  Int iCStride = pcPicOrg->getStride(COMPONENT_Cb);
  UInt64 uiFiltDistCb = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cb), pcPicRest->getAddr(COMPONENT_Cb), iCWidth, iCHeight, iCStride);
  UInt64 uiFiltDistCr = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cr), pcPicRest->getAddr(COMPONENT_Cr), iCWidth, iCHeight, iCStride);
  UInt64 uiOrgDistCb = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cb), pcPicDec->getAddr(COMPONENT_Cb), iCWidth, iCHeight, iCStride);
  UInt64 uiOrgDistCr = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cr), pcPicDec->getAddr(COMPONENT_Cr), iCWidth, iCHeight, iCStride);
  
  m_pcTempAlfParam->chroma_idc = 0;
  if(uiOrgDistCb > uiFiltDistCb)
    m_pcTempAlfParam->chroma_idc += 2;
  if(uiOrgDistCr  > uiFiltDistCr )
    m_pcTempAlfParam->chroma_idc += 1;
  
  if(m_pcTempAlfParam->chroma_idc)
  {
    if(m_pcTempAlfParam->chroma_idc!=3)
    {
      // chroma filter re-design
      xFilteringFrameChroma(pcPicOrg, pcPicDec, pcPicRest);
    }
    
    UInt64 uiRate, uiDist;
    Double dCost;
    xCalcRDCostChroma(pcPicOrg, pcPicRest, m_pcTempAlfParam, uiRate, uiDist, dCost, pSlice);
    
    if( dCost < dMinCost )
    {
      copyALFParam(m_pcBestAlfParam, m_pcTempAlfParam);
#if !JVET_C0038_GALF
      predictALFCoeffChroma(m_pcBestAlfParam);
#endif
      ruiBits += uiRate;
      ruiDist += uiDist;
    }
    else
    {
      m_pcBestAlfParam->chroma_idc = 0;
      
      if((m_pcTempAlfParam->chroma_idc>>1)&0x01)
        pcPicDec->copyToPic(pcPicRest,COMPONENT_Cb);
      if(m_pcTempAlfParam->chroma_idc&0x01)
        pcPicDec->copyToPic(pcPicRest,COMPONENT_Cr);
      
      ruiBits += uiMinRate;
      ruiDist += uiMinDist;
    }
  }
  else
  {
    m_pcBestAlfParam->chroma_idc = 0;
    
    ruiBits += uiMinRate;
    ruiDist += uiMinDist;
    
    pcPicDec->copyToPic(pcPicRest,COMPONENT_Cb);
    pcPicDec->copyToPic(pcPicRest,COMPONENT_Cr);
  }
}

// ====================================================================================================================
// Private member functions
// ====================================================================================================================

Void TEncAdaptiveLoopFilter::xInitParam()
{
  Int i, j;
  
  if (m_ppdAlfCorr != NULL)
  {
    for (i = 0; i < m_ALF_MAX_NUM_COEF; i++)
    {
      for (j = 0; j < m_ALF_MAX_NUM_COEF+1; j++)
      {
        m_ppdAlfCorr[i][j] = 0;
      }
    }
  }
  else
  {
    m_ppdAlfCorr = new Double*[m_ALF_MAX_NUM_COEF];
    for (i = 0; i < m_ALF_MAX_NUM_COEF; i++)
    {
      m_ppdAlfCorr[i] = new Double[m_ALF_MAX_NUM_COEF+1];
      for (j = 0; j < m_ALF_MAX_NUM_COEF+1; j++)
      {
        m_ppdAlfCorr[i][j] = 0;
      }
    }
  }
  
  if (m_pdDoubleAlfCoeff != NULL)
  {
    for (i = 0; i < m_ALF_MAX_NUM_COEF; i++)
    {
      m_pdDoubleAlfCoeff[i] = 0;
    }
  }
  else
  {
    m_pdDoubleAlfCoeff = new Double[m_ALF_MAX_NUM_COEF];
    for (i = 0; i < m_ALF_MAX_NUM_COEF; i++)
    {
      m_pdDoubleAlfCoeff[i] = 0;
    }
  }
}

Void TEncAdaptiveLoopFilter::xUninitParam()
{
  Int i;
  
  if (m_ppdAlfCorr != NULL)
  {
    for (i = 0; i < m_ALF_MAX_NUM_COEF; i++)
    {
      delete[] m_ppdAlfCorr[i];
      m_ppdAlfCorr[i] = NULL;
    }
    delete[] m_ppdAlfCorr;
    m_ppdAlfCorr = NULL;
  }
  
  if (m_pdDoubleAlfCoeff != NULL)
  {
    delete[] m_pdDoubleAlfCoeff;
    m_pdDoubleAlfCoeff = NULL;
  }
}

Void TEncAdaptiveLoopFilter::xCreateTmpAlfCtrlFlags()
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
    pcCU->createTmpAlfCtrlFlag();
  }
}

Void TEncAdaptiveLoopFilter::xDestroyTmpAlfCtrlFlags()
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
    pcCU->destroyTmpAlfCtrlFlag();
  }
}

Void TEncAdaptiveLoopFilter::xCopyTmpAlfCtrlFlagsTo()
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
    pcCU->copyAlfCtrlFlagFromTmp();
  }
}

Void TEncAdaptiveLoopFilter::xCopyTmpAlfCtrlFlagsFrom()
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
    pcCU->copyAlfCtrlFlagToTmp();
  }
}

#if !JVET_C0024_QTBT
Void TEncAdaptiveLoopFilter::xEncodeCUAlfCtrlFlags()
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
    xEncodeCUAlfCtrlFlag(pcCU, 0, 0);
  }
}

Void TEncAdaptiveLoopFilter::xEncodeCUAlfCtrlFlag(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
  Bool bBoundary = false;
  UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiRPelX   = uiLPelX + (pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth)  - 1;
#else
  UInt uiRPelX   = uiLPelX + (pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth)  - 1;
#endif
  UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiBPelY   = uiTPelY + (pcCU->getSlice()->getSPS()->getCTUSize()>>uiDepth) - 1;
#else
  UInt uiBPelY   = uiTPelY + (pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth) - 1;
#endif
  
  if( ( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) || ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
  {
    bBoundary = true;
  }
  
#if JVET_C0024_QTBT
  if( ( uiDepth < pcCU->getDepth( uiAbsPartIdx )  ) || bBoundary )
#else
  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
#endif
  {
    UInt uiQNumParts = ( m_pcPic->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
      
      if( ( uiLPelX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
        xEncodeCUAlfCtrlFlag(pcCU, uiAbsPartIdx, uiDepth+1);
    }
    return;
  }
  
  m_pcEntropyCoder->encodeAlfCtrlFlag(pcCU, uiAbsPartIdx);
}

#endif
#if JVET_C0038_GALF 
Void TEncAdaptiveLoopFilter::xCalcCorrelationFunc(Pel* pOrg, Pel* pCmp, Int iTap, Int iWidth, Int iHeight, Int iOrgStride, Int iCmpStride, Int N)
#else
Void TEncAdaptiveLoopFilter::xCalcCorrelationFunc(Pel* pOrg, Pel* pCmp, Int iTap, Int iWidth, Int iHeight, Int iOrgStride, Int iCmpStride)
#endif
{ 
  Int iTapV   = TComAdaptiveLoopFilter::ALFTapHToTapV(iTap);
#if !JVET_C0038_GALF 
  Int N       = (iTap * iTapV + 1) >> 1;
#endif
  Int offsetV = iTapV >> 1;
  Int offset = iTap>>1;

#if JVET_C0038_GALF 
  const Int* pFiltPos = iTap==5 ? m_aiSymmetricArray5x5 : (iTap==7 ? m_aiSymmetricArray7x7:m_aiSymmetricArray9x9);  
#else
  const Int* pFiltPos;
  
  switch(iTap)
  {
    case 5:
      pFiltPos = m_aiSymmetricArray5x5;
      break;
    case 7:
      pFiltPos = m_aiSymmetricArray7x7;
      break;
    case 9:
      pFiltPos = m_aiSymmetricArray9x7;
      break;
    default:
      pFiltPos = m_aiSymmetricArray9x7;
      assert(0);
      break;
  }
#endif

  Pel* pTerm = new Pel[N];
  
  Int i, j;
  for (Int y = 0; y < iHeight; y++)
  {
    for (Int x = 0; x < iWidth; x++)
    {
      i = 0;
      ::memset(pTerm, 0, sizeof(Pel)*N);
      for (Int yy = y - offsetV; yy <= y + offsetV; yy++)
      {
        for(Int xx=x-offset; xx<=x+offset; xx++)
        {
#if JVET_C0038_GALF 
          if(pFiltPos[i])
          {
            pTerm[pFiltPos[i]-1] += pCmp[xx + yy*iCmpStride];
          }
#else
          pTerm[pFiltPos[i]] += pCmp[xx + yy*iCmpStride];
#endif
          i++;
        }
      }
      
      for(j=0; j<N; j++)
      {
        m_ppdAlfCorr[j][j] += pTerm[j]*pTerm[j];   //auto  m_EGlobalSym
        for(i=j+1; i<N; i++)
          m_ppdAlfCorr[j][i] += pTerm[j]*pTerm[i]; //cross corr
  
#if JVET_C0038_GALF 
        m_ppdAlfCorr[j][N] += pOrg[x+y*iOrgStride]*pTerm[j]; //m_yGlobalSym
#else
        // DC offset
        m_ppdAlfCorr[j][N]   += pTerm[j]; //ELocal[p_pattern[k++]]
        m_ppdAlfCorr[j][N+1] += pOrg[x+y*iOrgStride]*pTerm[j]; //m_yGlobalSym
#endif
      }

#if !JVET_C0038_GALF 
      // DC offset
      for(i=0; i<N; i++)
        m_ppdAlfCorr[N][i] += pTerm[i];
      m_ppdAlfCorr[N][N]   += 1;
      m_ppdAlfCorr[N][N+1] += pOrg[x+y*iOrgStride];
#endif
    }
  }
  for(j=0; j<N-1; j++)
  {
    for(i=j+1; i<N; i++)
      m_ppdAlfCorr[i][j] = m_ppdAlfCorr[j][i];
  }

  delete[] pTerm;
  pTerm = NULL;
}

UInt64 TEncAdaptiveLoopFilter::xCalcSSD(Pel* pOrg, Pel* pCmp, Int iWidth, Int iHeight, Int iStride )
{
  UInt64 uiSSD = 0;
  Int x, y;
  
  UInt uiShift = m_nBitIncrement<<1;
  Int iTemp;
  
  for( y = 0; y < iHeight; y++ )
  {
    for( x = 0; x < iWidth; x++ )
    {
      iTemp = pOrg[x] - pCmp[x]; uiSSD += ( iTemp * iTemp ) >> uiShift;
    }
    pOrg += iStride;
    pCmp += iStride;
  }
  
  return uiSSD;;
}

Int TEncAdaptiveLoopFilter::xGauss(Double **a, Int N)
{
  Int i, j, k;
  Double t;
  
  for(k=0; k<N; k++)
  {
    if (a[k][k] <0.000001)
      return 1;
  }
  
  for(k=0; k<N-1; k++)
  {
    for(i=k+1;i<N; i++)
    {
      t=a[i][k]/a[k][k];
      for(j=k+1; j<=N; j++)
      {
        a[i][j] -= t * a[k][j];
        if(i==j && fabs(a[i][j])<0.000001) return 1;
      }
    }
  }
  for(i=N-1; i>=0; i--)
  {
    t = a[i][N];
    for(j=i+1; j<N; j++)
      t -= a[i][j] * a[j][N];
    a[i][N] = t / a[i][i];
  }
  return 0;
}

Void TEncAdaptiveLoopFilter::xFilterCoefQuickSort( Double *coef_data, Int *coef_num, Int upper, Int lower )
{
  Double mid, tmp_data;
  Int i, j, tmp_num;
  
  i = upper;
  j = lower;
  mid = coef_data[(lower+upper)>>1];
  do
  {
    while( coef_data[i] < mid ) i++;
    while( mid < coef_data[j] ) j--;
    if( i <= j )
    {
      tmp_data = coef_data[i];
      tmp_num  = coef_num[i];
      coef_data[i] = coef_data[j];
      coef_num[i]  = coef_num[j];
      coef_data[j] = tmp_data;
      coef_num[j]  = tmp_num;
      i++;
      j--;
    }
  } while( i <= j );
  if ( upper < j ) xFilterCoefQuickSort(coef_data, coef_num, upper, j);
  if ( i < lower ) xFilterCoefQuickSort(coef_data, coef_num, i, lower);
}

Void TEncAdaptiveLoopFilter::xQuantFilterCoef(Double* h, Int* qh, Int tap, int bit_depth)
{
  Int i, N;
  Int max_value, min_value;
  Double dbl_total_gain;
  Int total_gain, q_total_gain;
  Int upper, lower;
  Double *dh;
  Int    *nc;

#if JVET_C0038_GALF 
  Int factor     = (1 << m_ALF_NUM_BIT_SHIFT); 
  Int filtNo     = (tap==9) ? 0 : (tap==7 ? 1: 2);
  Int * pFiltMag = m_weightsTab[filtNo]; 
#else
  const Int    *pFiltMag;

  switch(tap)
  {
    case 5:
      pFiltMag = m_aiSymmetricMag5x5;
      break;
    case 7:
      pFiltMag = m_aiSymmetricMag7x7;
      break;
    case 9:
      pFiltMag = m_aiSymmetricMag9x7;
      break;
    default:
      pFiltMag = m_aiSymmetricMag9x7;
      assert(0);
      break;
  }
#endif  
#if JVET_C0038_GALF 
  N = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(tap);
#else
  Int tapV = TComAdaptiveLoopFilter::ALFTapHToTapV(tap);
  N = (tap * tapV + 1) >> 1;
#endif

  dh = new Double[N];
  nc = new Int[N];
  
  max_value =   (1<<(1+m_ALF_NUM_BIT_SHIFT))-1;
  min_value = 0-(1<<(1+m_ALF_NUM_BIT_SHIFT));
  
  dbl_total_gain=0.0;
  q_total_gain=0;
#if JVET_C0038_GALF 
  Int quantCoeffSum=0;
#endif

  for(i=0; i<N; i++)
  {
#if JVET_C0038_GALF 
    if(i == (N-1))
    {
      qh[i]=factor-quantCoeffSum;
    }
    else
    {
#endif
    if(h[i]>=0.0)
      qh[i] =  (Int)( h[i]*(1<<m_ALF_NUM_BIT_SHIFT)+0.5);
    else
      qh[i] = -(Int)(-h[i]*(1<<m_ALF_NUM_BIT_SHIFT)+0.5);
#if JVET_C0038_GALF 
      quantCoeffSum += pFiltMag[i]*qh[i];
    }
#endif
    dh[i] = (Double)qh[i]/(Double)(1<<m_ALF_NUM_BIT_SHIFT) - h[i];
    dh[i]*=pFiltMag[i];
    dbl_total_gain += h[i]*pFiltMag[i];
    q_total_gain   += qh[i]*pFiltMag[i];
    nc[i] = i;
  }
  
  // modification of quantized filter coefficients
  total_gain = (Int)(dbl_total_gain*(1<<m_ALF_NUM_BIT_SHIFT)+0.5);
  
  if( q_total_gain != total_gain )
  {
    xFilterCoefQuickSort(dh, nc, 0, N-1);
    if( q_total_gain > total_gain )
    {
      upper = N-1;
      while( q_total_gain > total_gain+1 )
      {
        i = nc[upper%N];
        qh[i]--;
        q_total_gain -= pFiltMag[i];
        upper--;
      }
      if( q_total_gain == total_gain+1 )
      {
        if(dh[N-1]>0)
          qh[N-1]--;
        else
        {
          i=nc[upper%N];
          qh[i]--;
          qh[N-1]++;
        }
      }
    }
    else if( q_total_gain < total_gain )
    {
      lower = 0;
      while( q_total_gain < total_gain-1 )
      {
        i=nc[lower%N];
        qh[i]++;
        q_total_gain += pFiltMag[i];
        lower++;
      }
      if( q_total_gain == total_gain-1 )
      {
        if(dh[N-1]<0)
          qh[N-1]++;
        else
        {
          i=nc[lower%N];
          qh[i]++;
          qh[N-1]--;
        }
      }
    }
  }
  
  // set of filter coefficients
#if JVET_C0038_GALF 
  quantCoeffSum=0;
  for (i=0; i< N-1; i++)
  {
    quantCoeffSum += pFiltMag[i]*qh[i];
  }
  qh[N-1]= factor - quantCoeffSum;
#endif

  for(i=0; i<N; i++)
  {
    qh[i] = std::max(min_value,std::min(max_value, qh[i]));
  }
#if !JVET_C0038_GALF 
  // DC offset
  //  max_value = Min(  (1<<(3+Max(img_bitdepth_luma,img_bitdepth_chroma)))-1, (1<<14)-1);
  //  min_value = Max( -(1<<(3+Max(img_bitdepth_luma,img_bitdepth_chroma))),  -(1<<14)  );
  max_value = std::min(  (1<<(3+m_nInputBitDepth + m_nBitIncrement))-1, (1<<14)-1);
  min_value = std::max( -(1<<(3+m_nInputBitDepth + m_nBitIncrement)),  -(1<<14)  );

  qh[N] =  (h[N]>=0.0)? (Int)( h[N]*(1<<(m_ALF_NUM_BIT_SHIFT-bit_depth+8)) + 0.5) : -(Int)(-h[N]*(1<<(m_ALF_NUM_BIT_SHIFT-bit_depth+8)) + 0.5);
  qh[N] = std::max(min_value,std::min(max_value, qh[N]));
#endif

  delete[] dh;
  dh = NULL;
  
  delete[] nc;
  nc = NULL;
}

Void TEncAdaptiveLoopFilter::xClearFilterCoefInt(Int* qh, Int N)
{
  // clear
  memset( qh, 0, sizeof( Int ) * N );
  
  // center pos
#if JVET_C0038_GALF
  qh[N-1]  = 1<<m_ALF_NUM_BIT_SHIFT;
#else
  qh[N-2]  = 1<<m_ALF_NUM_BIT_SHIFT;
#endif
}

Void TEncAdaptiveLoopFilter::xCalcRDCost(ALFParam* pAlfParam, UInt64& ruiRate, UInt64 uiDist, Double& rdCost, const TComSlice * pSlice)
{
  if(pAlfParam != NULL)
  {
#if !JVET_C0038_GALF
    Int* piTmpCoef;
    piTmpCoef = new Int[m_ALF_MAX_NUM_COEF];
    
    memcpy(piTmpCoef, pAlfParam->coeff, sizeof(Int)*pAlfParam->num_coeff);
#if COM16_C806_ALF_TEMPPRED_NUM
    if( !pAlfParam->temproalPredFlag )
#endif
    predictALFCoeff(pAlfParam);
#endif

    m_pcEntropyCoder->resetEntropy(pSlice);
    m_pcEntropyCoder->resetBits();

    m_pcEntropyCoder->encodeAlfParam(pAlfParam, pSlice->getSPS()->getMaxTotalCUDepth()
#if FIX_TICKET12
        , pSlice
#endif
      );
    
    if(pAlfParam->cu_control_flag)
    {
      m_pcEntropyCoder->encodeAlfCtrlParam(pAlfParam);
    }

    ruiRate = m_pcEntropyCoder->getNumberOfWrittenBits();

#if !JVET_C0038_GALF
    memcpy(pAlfParam->coeff, piTmpCoef, sizeof(int)*pAlfParam->num_coeff);
    delete[] piTmpCoef;
    piTmpCoef = NULL;
#endif
  }
  else
  {
    ruiRate = 1;
  }
  
  rdCost      = (Double)(ruiRate) * m_dLambdaLuma + (Double)(uiDist);
}

Void TEncAdaptiveLoopFilter::xCalcRDCost(TComPicYuv* pcPicOrg, TComPicYuv* pcPicCmp, ALFParam* pAlfParam, UInt64& ruiRate, UInt64& ruiDist, Double& rdCost, const TComSlice * pSlice)
{
  if(pAlfParam != NULL)
  {
#if !JVET_C0038_GALF
    Int* piTmpCoef;
    piTmpCoef = new Int[m_ALF_MAX_NUM_COEF];
    
    memcpy(piTmpCoef, pAlfParam->coeff, sizeof(Int)*pAlfParam->num_coeff);
#if COM16_C806_ALF_TEMPPRED_NUM
    if( !pAlfParam->temproalPredFlag )
#endif
    predictALFCoeff(pAlfParam);
#endif

    m_pcEntropyCoder->resetEntropy(pSlice);
    m_pcEntropyCoder->resetBits();

    m_pcEntropyCoder->encodeAlfParam(pAlfParam, pSlice->getSPS()->getMaxTotalCUDepth()
#if FIX_TICKET12
        , pSlice
#endif
      );
    
    if(pAlfParam->cu_control_flag)
    {
      m_pcEntropyCoder->encodeAlfCtrlParam(pAlfParam);
    }

    ruiRate = m_pcEntropyCoder->getNumberOfWrittenBits();

#if !JVET_C0038_GALF
    memcpy(pAlfParam->coeff, piTmpCoef, sizeof(int)*pAlfParam->num_coeff);
    delete[] piTmpCoef;
    piTmpCoef = NULL;
#endif
  }
  else
  {
    ruiRate = 1;
  }
  
  ruiDist     = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Y), pcPicCmp->getAddr(COMPONENT_Y), pcPicOrg->getWidth(COMPONENT_Y), pcPicOrg->getHeight(COMPONENT_Y), pcPicOrg->getStride(COMPONENT_Y));
  rdCost      = (Double)(ruiRate) * m_dLambdaLuma + (Double)(ruiDist);
}

Void TEncAdaptiveLoopFilter::xCalcRDCostChroma(TComPicYuv* pcPicOrg, TComPicYuv* pcPicCmp, ALFParam* pAlfParam, UInt64& ruiRate, UInt64& ruiDist, Double& rdCost, const TComSlice * pSlice)
{
  if(pAlfParam->chroma_idc)
  {
#if !JVET_C0038_GALF
    Int* piTmpCoef;
    piTmpCoef = new Int[m_ALF_MAX_NUM_COEF_C];
    
    memcpy(piTmpCoef, pAlfParam->coeff_chroma, sizeof(Int)*pAlfParam->num_coeff_chroma);
    
    predictALFCoeffChroma(pAlfParam);
#endif
    m_pcEntropyCoder->resetEntropy(pSlice);
    m_pcEntropyCoder->resetBits();

    m_pcEntropyCoder->encodeAlfParam(pAlfParam, pSlice->getSPS()->getMaxTotalCUDepth()
#if FIX_TICKET12
        , pSlice
#endif
      );
    
    if(pAlfParam->cu_control_flag)
    {
      m_pcEntropyCoder->encodeAlfCtrlParam(pAlfParam);
    }

    ruiRate = m_pcEntropyCoder->getNumberOfWrittenBits();

#if !JVET_C0038_GALF
    memcpy(pAlfParam->coeff_chroma, piTmpCoef, sizeof(int)*pAlfParam->num_coeff_chroma);
    delete[] piTmpCoef;
    piTmpCoef = NULL;
#endif
  }
  ruiDist = 0;
  ruiDist += xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cb), pcPicCmp->getAddr(COMPONENT_Cb), pcPicOrg->getWidth(COMPONENT_Cb), pcPicOrg->getHeight(COMPONENT_Cb), pcPicOrg->getStride(COMPONENT_Cb));
  ruiDist += xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cr), pcPicCmp->getAddr(COMPONENT_Cr), pcPicOrg->getWidth(COMPONENT_Cr), pcPicOrg->getHeight(COMPONENT_Cr), pcPicOrg->getStride(COMPONENT_Cr));
  rdCost  = (Double)(ruiRate) * m_dLambdaChroma + (Double)(ruiDist);
}

Void TEncAdaptiveLoopFilter::xFilteringFrameChroma(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest)
{
#if JVET_C0038_GALF
  Int tap;
#else
  Int    i, tap, N, err_code;
#endif
  Int* qh;
  
  tap  = m_pcTempAlfParam->tap_chroma;
#if !JVET_C0038_GALF
  N    = m_pcTempAlfParam->num_coeff_chroma;
#endif
  qh   = m_pcTempAlfParam->coeff_chroma;

#if JVET_C0038_GALF
  if(m_pcTempAlfParam->chroma_idc)
  {
    Int filtNo = 2;
    if (tap==9)
      filtNo =0;
    else if (tap==7)
      filtNo =1;

    xstoreInBlockMatrixForChroma(pcPicOrg, pcPicDec, tap, m_pcTempAlfParam->chroma_idc);
    {
      //static double pixAcc_temp;
      Int sqrFiltLength = m_sqrFiltLengthTab[filtNo];   
      Int* weights      = m_weightsTab      [filtNo];  
      Int  bit_depth    = m_NUM_BITS;
      // Find coeffcients
      QuantizeIntegerFilterPP(m_filterCoeff, qh, m_EGlobalSym[filtNo][0], m_yGlobalSym[filtNo][0], sqrFiltLength, weights, bit_depth, true); 
#if COM16_C806_ALF_TEMPPRED_NUM
      memcpy( m_pcTempAlfParam->alfCoeffChroma, qh, sizeof(Int)*sqrFiltLength );
#endif
    }
  }
#else
  // initialize correlation
  for(i=0; i<N; i++)
    memset(m_ppdAlfCorr[i], 0, sizeof(Double)*(N+1));
  
  if ((m_pcTempAlfParam->chroma_idc>>1)&0x01)
  {
    Pel* pOrg = pcPicOrg->getAddr(COMPONENT_Cb);
    Pel* pCmp = pcPicDec->getAddr(COMPONENT_Cb);
    xCalcCorrelationFunc(pOrg, pCmp, tap, pcPicOrg->getWidth(COMPONENT_Cb), pcPicOrg->getHeight(COMPONENT_Cb), pcPicOrg->getStride(COMPONENT_Cb), pcPicDec->getStride(COMPONENT_Cb));
  }
  if ((m_pcTempAlfParam->chroma_idc)&0x01)
  {
    Pel* pOrg = pcPicOrg->getAddr(COMPONENT_Cr);
    Pel* pCmp = pcPicDec->getAddr(COMPONENT_Cr);
    xCalcCorrelationFunc(pOrg, pCmp, tap, pcPicOrg->getWidth(COMPONENT_Cr), pcPicOrg->getHeight(COMPONENT_Cr), pcPicOrg->getStride(COMPONENT_Cr), pcPicDec->getStride(COMPONENT_Cr));
  }
  
  err_code = xGauss(m_ppdAlfCorr, N);
  
  if(err_code)
  {
    xClearFilterCoefInt(qh, N);
  }
  else
  {
    for(i=0; i<N; i++)
      m_pdDoubleAlfCoeff[i] = m_ppdAlfCorr[i][N];
    
    xQuantFilterCoef(m_pdDoubleAlfCoeff, qh, tap, m_nInputBitDepth + m_nBitIncrement);
  }

#if COM16_C806_ALF_TEMPPRED_NUM 
  memcpy( m_pcTempAlfParam->alfCoeffChroma, qh, sizeof(Int)*m_ALF_MAX_NUM_COEF_C );
#endif
#endif

#if JVET_C0038_GALF
  if(m_pcTempAlfParam->chroma_idc)
  {
    initVarForChroma(m_pcTempAlfParam);
  }
  if ((m_pcTempAlfParam->chroma_idc>>1)&0x01)
  {
    xFrameChroma( m_pcTempAlfParam, pcPicDec, pcPicRest, 0 );
  }
  if ((m_pcTempAlfParam->chroma_idc)&0x01)
  {
    xFrameChroma( m_pcTempAlfParam, pcPicDec, pcPicRest, 1 );
  }
#else
  if ((m_pcTempAlfParam->chroma_idc>>1)&0x01)
  {
    xFrameChroma(pcPicDec, pcPicRest, qh, tap, 0);
  }
  if ((m_pcTempAlfParam->chroma_idc)&0x01)
  {
    xFrameChroma(pcPicDec, pcPicRest, qh, tap, 1);
  }
#endif 
  if(m_pcTempAlfParam->chroma_idc<3)
  {
    if(m_pcTempAlfParam->chroma_idc==1)
    {
      pcPicDec->copyToPic(pcPicRest, COMPONENT_Cb);
    }
    if(m_pcTempAlfParam->chroma_idc==2)
    {
      pcPicDec->copyToPic(pcPicRest, COMPONENT_Cr);
    }
  }
  
}

Void TEncAdaptiveLoopFilter::xCopyDecToRestCUs(TComPicYuv* pcPicDec, TComPicYuv* pcPicRest)
{
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
#if JVET_C0024_QTBT
    xCopyDecToRestCU(pcCU, 0, 0, pcCU->getSlice()->getSPS()->getCTUSize(), pcCU->getSlice()->getSPS()->getCTUSize(), pcPicDec, pcPicRest);
#else
    xCopyDecToRestCU(pcCU, 0, 0, pcPicDec, pcPicRest);
#endif
  }
}

#if JVET_C0024_QTBT
Void TEncAdaptiveLoopFilter::xCopyDecToRestCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest)
#else
Void TEncAdaptiveLoopFilter::xCopyDecToRestCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest)
#endif
{
  Bool bBoundary = false;
  UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiRPelX   = uiLPelX + uiWidth  - 1;
#else
  UInt uiRPelX   = uiLPelX + (pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth)  - 1;
#endif
  UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiBPelY   = uiTPelY + uiHeight - 1;
#else
  UInt uiBPelY   = uiTPelY + (pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth) - 1;
#endif
  
#if BARC_Boundary
  UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX =  (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
  UInt deltaY =  (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());

  if( ( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - deltaX) || ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - deltaY) )
#else
  if( ( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) || ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
#endif
  {
    bBoundary = true;
  }
#if JVET_C0024_QTBT
  if( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) || bBoundary )
#else
  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < (pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize()) ) ) || bBoundary )
#endif
  {
    UInt uiQNumParts = ( m_pcPic->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
      
#if BARC_Boundary
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if( ( uiLPelX + deltaX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY + deltaY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
#else
      if( ( uiLPelX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )      
#endif
#if JVET_C0024_QTBT
        xCopyDecToRestCU(pcCU, uiAbsPartIdx, uiDepth+1, uiWidth>>1, uiHeight>>1, pcPicDec, pcPicRest);
#else
        xCopyDecToRestCU(pcCU, uiAbsPartIdx, uiDepth+1, pcPicDec, pcPicRest);
#endif
    }
    return;
  }
  
  if (!pcCU->getAlfCtrlFlag(uiAbsPartIdx))
  {
    UInt uiCUAddr = pcCU->getCtuRsAddr();
    
#if JVET_C0024_QTBT
    Int iWidth = uiWidth;
    Int iHeight = uiHeight;
#else
    Int iWidth = pcCU->getWidth(uiAbsPartIdx);
    Int iHeight = pcCU->getHeight(uiAbsPartIdx);
#endif
    
    Pel* pRec = pcPicDec->getAddr(COMPONENT_Y, uiCUAddr, uiAbsPartIdx);
    Pel* pFilt = pcPicRest->getAddr(COMPONENT_Y, uiCUAddr, uiAbsPartIdx);
    
    Int iRecStride = pcPicDec->getStride(COMPONENT_Y);
    Int iFiltStride = pcPicRest->getStride(COMPONENT_Y);
    
    for (Int y = 0; y < iHeight; y++)
    {
      for (Int x = 0; x < iWidth; x++)
      {
        pFilt[x] = pRec[x];
      }
      pRec += iRecStride;
      pFilt += iFiltStride;
    }
  }
}

Void TEncAdaptiveLoopFilter::xcollectStatCodeFilterCoeffForce0(int **pDiffQFilterCoeffIntPP, int fl, int sqrFiltLength, 
                                                               int filters_per_group, int bitsVarBin[])
{
  int i, k, kMin, kStart, minBits, ind, scanPos, maxScanVal, coeffVal, 
    kMinTab[m_MAX_SQR_FILT_LENGTH], bitsCoeffScan[m_MAX_SCAN_VAL][m_MAX_EXP_GOLOMB],
    minKStart, minBitsKStart, bitsKStart;
  
  const Int * pDepthInt=TComAdaptiveLoopFilter::m_pDepthIntTab[fl-2];

#if JVET_C0038_GALF
  sqrFiltLength=sqrFiltLength-1;
#endif

  maxScanVal=0;
  for (i=0; i<sqrFiltLength; i++)
  {
    maxScanVal=max(maxScanVal, pDepthInt[i]);
  }
  
  // vlc for all
  memset(bitsCoeffScan, 0, m_MAX_SCAN_VAL * m_MAX_EXP_GOLOMB * sizeof(int));
  for(ind=0; ind<filters_per_group; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos=pDepthInt[i]-1;
      coeffVal=abs(pDiffQFilterCoeffIntPP[ind][i]);
      for (k=1; k<15; k++)
      {
        bitsCoeffScan[scanPos][k] += lengthGolomb(coeffVal, k);
      }
    }
  }
  
  minBitsKStart=0;
  minKStart = -1;
  for (k=1; k<8; k++)
  { 
    bitsKStart=0; kStart=k;
    for (scanPos=0; scanPos<maxScanVal; scanPos++)
    {
      kMin=kStart; minBits=bitsCoeffScan[scanPos][kMin];
      
      if (bitsCoeffScan[scanPos][kStart+1]<minBits)
      {
        kMin=kStart+1; minBits=bitsCoeffScan[scanPos][kMin];
      }
      kStart=kMin;
      bitsKStart+=minBits;
    }
    if (bitsKStart<minBitsKStart || k==1)
    {
      minBitsKStart=bitsKStart;
      minKStart=k;
    }
  }
  
  kStart = minKStart; 
  for (scanPos=0; scanPos<maxScanVal; scanPos++)
  {
    kMin=kStart; minBits=bitsCoeffScan[scanPos][kMin];
    
    if (bitsCoeffScan[scanPos][kStart+1]<minBits)
    {
      kMin = kStart+1; 
      minBits = bitsCoeffScan[scanPos][kMin];
    }
    
    kMinTab[scanPos] = kMin;
    kStart = kMin;
  }
  
  for(ind=0; ind<filters_per_group; ++ind)
  {
    bitsVarBin[ind]=0;
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos=pDepthInt[i]-1;
      bitsVarBin[ind] += lengthGolomb(abs(pDiffQFilterCoeffIntPP[ind][i]), kMinTab[scanPos]);
    }
  }
}
#if JVET_C0038_GALF  
Void TEncAdaptiveLoopFilter::xdecideCoeffForce0(Bool codedVarBins[m_NO_VAR_BINS], Double errorForce0Coeff[], Double errorForce0CoeffTab[m_NO_VAR_BINS][2], Int bitsVarBin[m_NO_VAR_BINS], Double lambda, Int filters_per_fr)
#else
Void TEncAdaptiveLoopFilter::xdecideCoeffForce0(int codedVarBins[m_NO_VAR_BINS], double errorForce0Coeff[], double errorForce0CoeffTab[m_NO_VAR_BINS][2], int bitsVarBin[m_NO_VAR_BINS], double lambda, int filters_per_fr)
#endif
{
  int filtNo;
  
  int ind;
#if JVET_C0038_GALF  
  errorForce0Coeff[0] = 0;
  for (ind=0; ind< m_NO_VAR_BINS; ind++)
  {
    codedVarBins[ind] = false;
  }
#else
  errorForce0Coeff[0]=errorForce0Coeff[1]=0;
  for (ind=0; ind<16; ind++) codedVarBins[ind]=0;
#endif
  
  for(filtNo=0; filtNo<filters_per_fr; filtNo++)
  {
#if JVET_C0038_GALF
    Double lagrangianDiff=errorForce0CoeffTab[filtNo][0]-(errorForce0CoeffTab[filtNo][1]+lambda*bitsVarBin[filtNo]);
    codedVarBins[filtNo] = (lagrangianDiff>0) ? true : false;
    errorForce0Coeff[0] += errorForce0CoeffTab[filtNo][codedVarBins[filtNo]?1:0];
#else
    // No coeffcient prediction bits used
    //lagrangianDiff=errorForce0CoeffTab[filtNo][0]-(errorForce0CoeffTab[filtNo][1]+lambda*bitsVarBin[filtNo]);
    codedVarBins[filtNo]= 1;
    errorForce0Coeff[0]+=errorForce0CoeffTab[filtNo][codedVarBins[filtNo]];
    errorForce0Coeff[1]+=errorForce0CoeffTab[filtNo][1];
#endif
  }   
}
#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::xInitFixedFilters()
{
  Int factor = (1<<(m_NUM_BITS-1));
  Int maxFilterLength = m_MAX_SQR_FILT_LENGTH/2 + 1;
  for(Int i = 0; i < maxFilterLength; i++)
  {
    for (Int j=0; j<m_NO_FILTERS*JVET_C0038_NO_PREV_FILTERS; j++)
    {
      m_filterCoeffPrev[j][i]=(Double)m_ALFfilterCoeffFixed[j][i]/(Double)factor;
    }
  }
  memset( m_filterCoeffDefault, 0, (maxFilterLength-1)*sizeof(Double) );
  m_filterCoeffDefault[maxFilterLength - 1] = 1.0;
}
Double TEncAdaptiveLoopFilter::xfindBestCoeffCodMethod(Bool codedVarBins[m_NO_VAR_BINS], Bool *forceCoeff0, Int fl, Int sqrFiltLength, 
                                                       Int filters_per_fr, Double errorForce0CoeffTab[m_NO_VAR_BINS][2], 
                                                       Double *errorQuant, Double lambda, const TComSlice * pSlice)
{
  Int bitsVarBin[m_NO_VAR_BINS], createBistream, coeffBits, coeffBitsForce0;
  Double errorForce0Coeff[2], lagrangianForce0, lagrangian;

  xcollectStatCodeFilterCoeffForce0(m_filterCoeffSymQuant, fl, sqrFiltLength,  
    filters_per_fr, bitsVarBin);

  xdecideCoeffForce0(codedVarBins, errorForce0Coeff, errorForce0CoeffTab, bitsVarBin, lambda, filters_per_fr);

  coeffBitsForce0 = xsendAllFiltersPPPredForce0(m_filterCoeffSymQuant, fl, sqrFiltLength, 
    filters_per_fr, codedVarBins, createBistream=0, m_tempALFp, pSlice);

  coeffBits = xsendAllFiltersPPPred(m_filterCoeffSymQuant, fl, sqrFiltLength, filters_per_fr, 
    createBistream=0, m_tempALFp, pSlice);

  lagrangianForce0 = (errorForce0Coeff[0]+lambda*coeffBitsForce0);
  lagrangian= (*errorQuant) + lambda * coeffBits;

  if (lagrangianForce0 < lagrangian)
  {
    *forceCoeff0 = true;
    return(lagrangianForce0);
  }
  else
  {
    *forceCoeff0 = false;
    return(lagrangian);
  }
}
#else
double TEncAdaptiveLoopFilter::xfindBestCoeffCodMethod(int codedVarBins[m_NO_VAR_BINS], int *forceCoeff0, 
                                                       int **filterCoeffSymQuant, int fl, int sqrFiltLength, 
                                                       int filters_per_fr, double errorForce0CoeffTab[m_NO_VAR_BINS][2], 
                                                       double *errorQuant, double lambda, const TComSlice * pSlice)

{
  int bitsVarBin[m_NO_VAR_BINS], createBistream, coeffBits, coeffBitsForce0;
  double errorForce0Coeff[2], lagrangianForce0, lagrangian;
  
  xcollectStatCodeFilterCoeffForce0(filterCoeffSymQuant, fl, sqrFiltLength,  
                                    filters_per_fr, bitsVarBin);
  
  xdecideCoeffForce0(codedVarBins, errorForce0Coeff, errorForce0CoeffTab, bitsVarBin, lambda, filters_per_fr);
  
  coeffBitsForce0 = xsendAllFiltersPPPredForce0(filterCoeffSymQuant, fl, sqrFiltLength, 
                                                filters_per_fr, codedVarBins, createBistream=0, m_tempALFp, pSlice);
  
  coeffBits = xsendAllFiltersPPPred(filterCoeffSymQuant, fl, sqrFiltLength, filters_per_fr, 
                                    createBistream=0, m_tempALFp, pSlice);
  
  lagrangianForce0=errorForce0Coeff[0]+lambda*coeffBitsForce0;
  lagrangian=errorForce0Coeff[1]+lambda*coeffBits;
  if (lagrangianForce0<lagrangian)
  {
    *errorQuant=errorForce0Coeff[0];
    *forceCoeff0=1;
    return(lagrangianForce0);
  }
  else
  {
    *errorQuant=errorForce0Coeff[1];
    *forceCoeff0=0;
    return(lagrangian);
  }
}
#endif
Int TEncAdaptiveLoopFilter::xsendAllFiltersPPPred(int **FilterCoeffQuant, int fl, int sqrFiltLength, 
                                                  int filters_per_group, int createBistream, ALFParam* ALFp, const TComSlice * pSlice)
{
#if JVET_C0038_GALF
  int ind, bit_ct0 = 0, i;
#else
  int ind, bit_ct = 0, bit_ct0 = 0, i;
#endif
  int predMethod = 0;
#if !JVET_C0038_GALF
  int force0 = 0;
#endif
  Int64 Newbit_ct;
  
  bit_ct0 = xcodeFilterCoeff(FilterCoeffQuant, fl, sqrFiltLength, filters_per_group, 0);
  
  for(ind = 0; ind < filters_per_group; ++ind)
  {
    if(ind == 0)
    {
      for(i = 0; i < sqrFiltLength; i++)
        m_diffFilterCoeffQuant[ind][i] = FilterCoeffQuant[ind][i];
    }
    else
    {
      for(i = 0; i < sqrFiltLength; i++)
        m_diffFilterCoeffQuant[ind][i] = FilterCoeffQuant[ind][i] - FilterCoeffQuant[ind-1][i];
    }
  }

#if JVET_C0038_GALF
  predMethod = (xcodeFilterCoeff(m_diffFilterCoeffQuant, fl, sqrFiltLength, filters_per_group, 0) >= bit_ct0) ? 0: 1;
  ALFp->forceCoeff0 = false;
#else
  if(xcodeFilterCoeff(m_diffFilterCoeffQuant, fl, sqrFiltLength, filters_per_group, 0) >= bit_ct0)
  {
    predMethod = 0;  
    if(filters_per_group > 1)
      bit_ct += lengthPredFlags(force0, predMethod, NULL, 0, createBistream);
    bit_ct += xcodeFilterCoeff(FilterCoeffQuant, fl, sqrFiltLength, filters_per_group, createBistream);
  }
  else
  {
    predMethod = 1;
    if(filters_per_group > 1)
      bit_ct += lengthPredFlags(force0, predMethod, NULL, 0, createBistream);
    bit_ct += xcodeFilterCoeff(m_diffFilterCoeffQuant, fl, sqrFiltLength, filters_per_group, createBistream);
  }


  ALFp->forceCoeff0 = 0;
#endif
  ALFp->filters_per_group_diff = filters_per_group;
  ALFp->filters_per_group = filters_per_group;
  ALFp->predMethod = predMethod;
  ALFp->num_coeff = sqrFiltLength;
  if (ALFp->num_coeff == m_SQR_FILT_LENGTH_5SYM)
    ALFp->realfiltNo=2;
  else if (ALFp->num_coeff == m_SQR_FILT_LENGTH_7SYM)
    ALFp->realfiltNo=1;
  else
    ALFp->realfiltNo=0;
  
  for(ind = 0; ind < filters_per_group; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      if (predMethod) ALFp->coeffmulti[ind][i] = m_diffFilterCoeffQuant[ind][i];
      else ALFp->coeffmulti[ind][i] = FilterCoeffQuant[ind][i];
    }
  }
#if JVET_C0038_GALF
  for(ind = 0; ind < ALF_NO_VAR_BIN; ind++)
  {
    ALFp->codedVarBins[ind] = true;
  }
  if(createBistream==0) //Li: use estimated bits
  {
    Newbit_ct = 1;
    Newbit_ct += (Int64)xcodeFilterCoeff(ALFp->coeffmulti, fl, sqrFiltLength, filters_per_group, createBistream);    
  }
  else  
#endif
  m_pcDummyEntropyCoder->codeFiltCountBit(ALFp, &Newbit_ct, pSlice);
  
  
  //  return(bit_ct);
  return ((Int)Newbit_ct);
}

#if JVET_C0038_GALF  
Int TEncAdaptiveLoopFilter::xsendAllFiltersPPPredForce0(Int **FilterCoeffQuant, Int fl, Int sqrFiltLength, Int filters_per_group, 
                                                        Bool codedVarBins[m_NO_VAR_BINS], Int createBistream, ALFParam* ALFp, const TComSlice * pSlice)
#else
Int TEncAdaptiveLoopFilter::xsendAllFiltersPPPredForce0(int **FilterCoeffQuant, int fl, int sqrFiltLength, int filters_per_group, 
                                                        int codedVarBins[m_NO_VAR_BINS], int createBistream, ALFParam* ALFp, const TComSlice * pSlice)
#endif
{
#if JVET_C0038_GALF
  int ind, i, j, len = 0 ;
#else
  int ind, bit_ct=0, bit_ct0, i, j;
  int filters_per_group_temp, filters_per_group_diff;
  int chosenPred = 0;
  int force0 = 1;
  Int64 Newbit_ct;
#endif
  i = 0;
  for(ind = 0; ind < filters_per_group; ind++)
  {
#if JVET_C0038_GALF
    if(codedVarBins[ind] == true)
#else
    if(codedVarBins[ind] == 1)
#endif
    {
      for(j = 0; j < sqrFiltLength; j++)
#if JVET_C0038_GALF
        m_FilterCoeffQuantTemp[ind][j]=FilterCoeffQuant[ind][j];
#else
        m_FilterCoeffQuantTemp[i][j]=FilterCoeffQuant[ind][j];
#endif
      i++;
    }
#if JVET_C0038_GALF
    else
    {
      for(j = 0; j < sqrFiltLength; j++)
      {
        m_FilterCoeffQuantTemp[ind][j] = 0;
      }
    }
#endif
  }
#if JVET_C0038_GALF
  if(createBistream==0) //Li: use estimated bits
  {
    len = 1;
    len += xcodeFilterCoeffForce0(m_FilterCoeffQuantTemp, fl, sqrFiltLength, filters_per_group, codedVarBins, createBistream);
    return (len);
  }
  else
  {
    Int64 Newbit_ct = 0;
    ALFp->forceCoeff0 = true;

    ALFp->filters_per_group_diff = filters_per_group;
    ALFp->filters_per_group = filters_per_group;
    ALFp->predMethod = 0;
    ALFp->num_coeff = sqrFiltLength;
    if (ALFp->num_coeff == m_SQR_FILT_LENGTH_5SYM)
      ALFp->realfiltNo=2;
    else if (ALFp->num_coeff == m_SQR_FILT_LENGTH_7SYM)
      ALFp->realfiltNo=1;
    else
      ALFp->realfiltNo=0;

    for(ind = 0; ind < filters_per_group; ++ind)
    {
      for(i = 0; i < sqrFiltLength; i++)
      {
        ALFp->coeffmulti[ind][i] = FilterCoeffQuant[ind][i];
      }
    }
    memcpy ( ALFp->codedVarBins, codedVarBins, sizeof(Bool)*ALF_NO_VAR_BIN );
    m_pcDummyEntropyCoder->codeFiltCountBit(ALFp, &Newbit_ct, pSlice);
    return ((Int)Newbit_ct);
  } 
#else
  filters_per_group_diff = filters_per_group_temp = i;
  
  for(ind = 0; ind < filters_per_group; ++ind)
  {
    if(ind == 0)
    {
      for(i = 0; i < sqrFiltLength; i++)
        m_diffFilterCoeffQuant[ind][i] = m_FilterCoeffQuantTemp[ind][i];
    }
    else
    {
      for(i = 0; i < sqrFiltLength; i++)
        m_diffFilterCoeffQuant[ind][i] = m_FilterCoeffQuantTemp[ind][i] - m_FilterCoeffQuantTemp[ind-1][i];
    }
  }
  
  if(!((filters_per_group_temp == 0) && (filters_per_group == 1)))
  {
    bit_ct0 = xcodeFilterCoeff(m_FilterCoeffQuantTemp, fl, sqrFiltLength, filters_per_group_temp, 0);
    
    if(xcodeFilterCoeff(m_diffFilterCoeffQuant, fl, sqrFiltLength, filters_per_group_diff, 0) >= bit_ct0)
    {
      chosenPred = 0;
      bit_ct += lengthPredFlags(force0, chosenPred, codedVarBins, filters_per_group, createBistream);
      bit_ct += xcodeFilterCoeff(m_FilterCoeffQuantTemp, fl, sqrFiltLength, filters_per_group_temp, createBistream);
    }
    else
    {
      chosenPred = 1;
      bit_ct += lengthPredFlags(force0, chosenPred, codedVarBins, filters_per_group, createBistream);
      bit_ct += xcodeFilterCoeff(m_diffFilterCoeffQuant, fl, sqrFiltLength, filters_per_group_temp, createBistream);
    }
  }
  ALFp->forceCoeff0 = 1;
  ALFp->predMethod = chosenPred;
  ALFp->filters_per_group_diff = filters_per_group_diff;
  ALFp->filters_per_group = filters_per_group;
  ALFp->num_coeff = sqrFiltLength;
  if (ALFp->num_coeff == m_SQR_FILT_LENGTH_5SYM)
    ALFp->realfiltNo=2;
  else if (ALFp->num_coeff == m_SQR_FILT_LENGTH_7SYM)
    ALFp->realfiltNo=1;
  else
    ALFp->realfiltNo=0;
  
  for(ind = 0; ind < filters_per_group; ++ind)
  {
    ALFp->codedVarBins[ind] = codedVarBins[ind];
  }
  for(ind = 0; ind < filters_per_group_diff; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      if (chosenPred) ALFp->coeffmulti[ind][i] = m_diffFilterCoeffQuant[ind][i];
      else ALFp->coeffmulti[ind][i] = m_FilterCoeffQuantTemp[ind][i];
    }
  }
  m_pcDummyEntropyCoder->codeFiltCountBit(ALFp, &Newbit_ct, pSlice);
  
  return ((Int)Newbit_ct);
#endif
}

//filtNo==-1/realfiltNo, noFilters=filters_per_frames, realfiltNo=filtNo
//Li's note: variable frNo is useless
#if JVET_C0038_GALF
Int TEncAdaptiveLoopFilter::xcodeAuxInfo(Int filtNo, Int noFilters, Int varIndTab[m_NO_VAR_BINS], Int createBitstream, ALFParam* ALFp, const TComSlice * pSlice, Char usePrevFiltBest[m_NO_VAR_BINS], UChar availableFilters, Int noVarBins, Bool forceCoeff0)
#else
Int TEncAdaptiveLoopFilter::xcodeAuxInfo(int filtNo, int noFilters, int varIndTab[m_NO_VAR_BINS], int frNo, int createBitstream,int realfiltNo, ALFParam* ALFp, const TComSlice * pSlice)
#endif
{
#if JVET_C0038_GALF
  Int64 NewbitCt;
  ALFp->realfiltNo = filtNo;
  ALFp->filtNo     = filtNo;
  ALFp->noFilters  = noFilters;
  assert(filtNo >= 0);

  ALFp->forceCoeff0 = forceCoeff0;
  
  memcpy (ALFp->filterPattern, varIndTab, m_max_NO_VAR_BINS * sizeof(Int));     
  memcpy (ALFp->PrevFiltIdx, usePrevFiltBest, m_max_NO_VAR_BINS * sizeof(Char));
  
  m_pcDummyEntropyCoder->codeAuxCountBit(ALFp, &NewbitCt, pSlice);

  return((Int)NewbitCt);
#else
  int i, filterPattern[m_NO_VAR_BINS], startSecondFilter=0, bitCt=0;
  Int64 NewbitCt;
  
  
  
  //send realfiltNo (tap related)
  ALFp->realfiltNo = realfiltNo;
  ALFp->filtNo = filtNo;
  
  if(filtNo >= 0)
  {
    // decide startSecondFilter and filterPattern
    if(noFilters > 1)
    {
      memset(filterPattern, 0, m_NO_VAR_BINS * sizeof(int)); 
      for(i = 1; i < m_NO_VAR_BINS; ++i)
      {
        if(varIndTab[i] != varIndTab[i-1])
        {
          filterPattern[i] = 1;
          startSecondFilter = i;
        }
      }
      memcpy (ALFp->filterPattern, filterPattern, m_NO_VAR_BINS * sizeof(int));
      ALFp->startSecondFilter = startSecondFilter;
    }
    
    //send noFilters (filters_per_frame)
    //0: filters_per_frame = 1
    //1: filters_per_frame = 2
    //2: filters_per_frame > 2 (exact number from filterPattern)

    ALFp->noFilters = min(noFilters-1,2);
    if (noFilters<=0) printf("error\n");
  }
  m_pcDummyEntropyCoder->codeAuxCountBit(ALFp, &NewbitCt, pSlice);
  bitCt = (int) NewbitCt;
  return(bitCt);
#endif
}
#if JVET_C0038_GALF
Int   TEncAdaptiveLoopFilter::xcodeFilterCoeffForce0(Int **pDiffQFilterCoeffIntPP, Int fl, Int sqrFiltLength, Int filters_per_group, Bool codedVarBins[], Int createBitstream)
{
  int i, k, kMin, kStart, minBits, ind, scanPos, maxScanVal, coeffVal, len = 0,
    kMinTab[m_MAX_SQR_FILT_LENGTH], bitsCoeffScan[m_MAX_SCAN_VAL][m_MAX_EXP_GOLOMB],
    minKStart, minBitsKStart, bitsKStart;
  
  const Int * pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl-2];
  sqrFiltLength=sqrFiltLength-1;
  maxScanVal = 0;
  for(i = 0; i < sqrFiltLength; i++)
  {
    maxScanVal = max(maxScanVal, pDepthInt[i]);
  }
  
  // vlc for all
  memset(bitsCoeffScan, 0, m_MAX_SCAN_VAL * m_MAX_EXP_GOLOMB * sizeof(int));
  for(ind=0; ind<filters_per_group; ++ind)
  {
    if(!codedVarBins[ind])
    {
      continue;
    }
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos=pDepthInt[i]-1;
      coeffVal=abs(pDiffQFilterCoeffIntPP[ind][i]);
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
  assert(createBitstream == 0);
  len += (3 + maxScanVal);
  len += filters_per_group;
  // Filter coefficients
  for(ind=0; ind<filters_per_group; ++ind)
  {
    if (codedVarBins[ind]==1)
    {
      for(i = 0; i < sqrFiltLength; i++)
      {
        scanPos=pDepthInt[i]-1;
        len += lengthGolomb(abs(pDiffQFilterCoeffIntPP[ind][i]), kMinTab[scanPos]);
      }
    }
  }

  return len;
}
#endif
Int   TEncAdaptiveLoopFilter::xcodeFilterCoeff(int **pDiffQFilterCoeffIntPP, int fl, int sqrFiltLength, 
                                               int filters_per_group, int createBitstream)
{
  int i, k, kMin, kStart, minBits, ind, scanPos, maxScanVal, coeffVal, len = 0,
    kMinTab[m_MAX_SQR_FILT_LENGTH], bitsCoeffScan[m_MAX_SCAN_VAL][m_MAX_EXP_GOLOMB],
    minKStart, minBitsKStart, bitsKStart;
  
  const Int * pDepthInt = TComAdaptiveLoopFilter::m_pDepthIntTab[fl-2];
#if JVET_C0038_GALF
  sqrFiltLength=sqrFiltLength-1;
#endif

  maxScanVal = 0;
  for(i = 0; i < sqrFiltLength; i++)
  {
    maxScanVal = max(maxScanVal, pDepthInt[i]);
  }
  
  // vlc for all
  memset(bitsCoeffScan, 0, m_MAX_SCAN_VAL * m_MAX_EXP_GOLOMB * sizeof(int));
  for(ind=0; ind<filters_per_group; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos=pDepthInt[i]-1;
      coeffVal=abs(pDiffQFilterCoeffIntPP[ind][i]);
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
  //  len += lengthFilterCodingParams(minKStart, maxScanVal, kMinTab, createBitstream);
  len += (3 + maxScanVal);
  
  // Filter coefficients
  len += lengthFilterCoeffs(sqrFiltLength, filters_per_group, pDepthInt, pDiffQFilterCoeffIntPP, 
                            kMinTab, createBitstream);
  
  return len;
}

Int TEncAdaptiveLoopFilter::lengthGolomb(int coeffVal, int k)
{
  int m = 2 << (k - 1);
  int q = coeffVal / m;
  if(coeffVal != 0)
    return(q + 2 + k);
  else
    return(q + 1 + k);
}
#if !JVET_C0038_GALF
Int TEncAdaptiveLoopFilter::lengthPredFlags(int force0, int predMethod, int codedVarBins[m_NO_VAR_BINS], 
                                            int filters_per_group, int createBitstream)
{
  int bit_cnt = 0;
  
  if(force0)
    bit_cnt = 2 + filters_per_group;
  else
    bit_cnt = 2;
  return bit_cnt;
  
}
#endif
//important
Int TEncAdaptiveLoopFilter::lengthFilterCoeffs(int sqrFiltLength, int filters_per_group, const int pDepthInt[], 
                                               int **FilterCoeff, int kMinTab[], int createBitstream)
{
  int ind, scanPos, i;
  int bit_cnt = 0;
  
  for(ind = 0; ind < filters_per_group; ++ind)
  {
    for(i = 0; i < sqrFiltLength; i++)
    {
      scanPos = pDepthInt[i] - 1;
      bit_cnt += lengthGolomb(abs(FilterCoeff[ind][i]), kMinTab[scanPos]);
    }
  }
  return bit_cnt;
}

Void   TEncAdaptiveLoopFilter::xEncALFLuma_qc ( TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiMinRate, UInt64& ruiMinDist, Double& rdMinCost, const TComSlice * pSlice )
{
  //pcPicDec: extended decoded
  //pcPicRest: original decoded: filtered signal will be stored
  
  UInt64  uiRate;
  UInt64  uiDist;
  Double dCost;
  Int    LumaStride = pcPicOrg->getStride(COMPONENT_Y);
  imgpel* pOrg = (imgpel*) pcPicOrg->getAddr(COMPONENT_Y);
  imgpel* pDec = (imgpel*) pcPicDec->getAddr(COMPONENT_Y);

#if JVET_C0038_GALF
  Int tap               = m_ALF_MAX_NUM_TAP;
  bTapDecision          = true; 
  bUpdateMatrix         = true;
  bFindBestFixedFilter  = false; //for initilization of selecting best fixed filter
#else
  Int tap               = m_ALF_MIN_NUM_TAP;
#endif

  m_pcTempAlfParam->tap = tap;
  m_pcTempAlfParam->tapV = TComAdaptiveLoopFilter::ALFTapHToTapV(m_pcTempAlfParam->tap);
  m_pcTempAlfParam->num_coeff = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(m_pcTempAlfParam->tap);

  Int64       iDist;
  Int64       iMinMethodDist = MAX_INT;
  UInt64      uiMinMethodRate = MAX_INT;
  Double      dMinMethodCost = MAX_DOUBLE;
  Double    adExtraCostReduction;
  ALFParam  cFrmAlfParam        ;
  ALFParam* pcAlfParam = NULL;

  pcAlfParam = &(cFrmAlfParam);
  allocALFParam(pcAlfParam);
  resetALFParam(pcAlfParam);

  pcAlfParam->alf_flag        = 1;
  pcAlfParam->chroma_idc      = 0;
  pcAlfParam->cu_control_flag = 0;
  pcAlfParam->tap = tap;
  pcAlfParam->tapV = TComAdaptiveLoopFilter::ALFTapHToTapV(pcAlfParam->tap);
  pcAlfParam->num_coeff = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(pcAlfParam->tap);

  adExtraCostReduction = 0.0;

  pcAlfParam       = &(cFrmAlfParam);
  m_varImg         = m_varImgMethods;

#if JVET_C0038_GALF
  setInitialMask(pcPicOrg, pcPicDec,pcAlfParam);
#else
  setInitialMask(pcPicOrg, pcPicDec);
#endif

  xFirstFilteringFrameLuma(pOrg, pDec, (imgpel*)m_pcPicYuvTmp->getAddr(COMPONENT_Y), pcAlfParam, pcAlfParam->tap, LumaStride, pSlice); 
  xCalcRDCost(pcPicOrg, m_pcPicYuvTmp, pcAlfParam, uiRate, uiDist, dCost, pSlice); 
  iDist = (Int64)uiDist;

  dCost -= adExtraCostReduction;

  if(dCost < dMinMethodCost)
  {
    dMinMethodCost = dCost;
    uiMinMethodRate= uiRate;
    iMinMethodDist = iDist;
    m_pcPicYuvTmp->copyToPic(pcPicRest, COMPONENT_Y);
  }  

  dMinMethodCost += adExtraCostReduction;


  m_varImg= m_varImgMethods;

  copyALFParam(m_pcBestAlfParam, &cFrmAlfParam); 

  ruiMinRate = uiMinMethodRate;
  ruiMinDist = (UInt64)iMinMethodDist;
  rdMinCost =  dMinMethodCost;


  freeALFParam(&cFrmAlfParam);
}

Void   TEncAdaptiveLoopFilter::xFirstFilteringFrameLuma(imgpel* ImgOrg, imgpel* ImgDec, imgpel* ImgRest, ALFParam* ALFp, Int tap, Int Stride, const TComSlice * pSlice)
{
  xstoreInBlockMatrix(ImgOrg, ImgDec, tap, Stride);
  xFilteringFrameLuma_qc(ImgOrg, ImgDec, ImgRest, ALFp, tap, Stride, pSlice);
}


#if JVET_C0038_GALF
Void   TEncAdaptiveLoopFilter::xstoreInBlockMatrixForChroma(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, Int tap, Int chroma_idc)
{
  imgpel* ImgOrg;
  imgpel* ImgDec;
  Int i,j, k,l, varInd = 0, ii, jj;
  Int x, y, yLocal;
  Int fl =tap/2;
  Int flV = TComAdaptiveLoopFilter::ALFFlHToFlV(fl);
  Int sqrFiltLength = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(tap);
  Int fl2=9/2; //extended size at each side of the frame
  Int ELocal[m_MAX_SQR_FILT_LENGTH];
  Int *p_pattern;
  Int filtNo = 2; //for chroma
  Int iImgHeight = pcPicDec->getHeight(COMPONENT_Cb);
  Int iImgWidth  = pcPicDec->getWidth (COMPONENT_Cb);
  Double **E,*yy;
 
  p_pattern= m_patternTab[filtNo];
  
  memset( m_pixAcc, 0,sizeof(Double)*m_NO_VAR_BINS);  
  {
    memset(m_yGlobalSym[filtNo][varInd], 0, sizeof(Double)*m_MAX_SQR_FILT_LENGTH);
    for (k=0; k<sqrFiltLength; k++)
    {
      memset(m_EGlobalSym[filtNo][varInd][k], 0, sizeof(Double)*m_MAX_SQR_FILT_LENGTH);
    }
  }
  for(Int iColorIdx = 0; iColorIdx < 2; iColorIdx ++)
  {
    if( (iColorIdx==0 && chroma_idc < 2) || (iColorIdx==1 && (chroma_idc & 0x01) == 0))
    {
      continue;
    }
    ImgOrg = (imgpel*) pcPicOrg->getAddr(ComponentID(iColorIdx+1));
    ImgDec = (imgpel*) pcPicDec->getAddr(ComponentID(iColorIdx+1));
    Int Stride = pcPicDec->getStride(ComponentID(iColorIdx+1));

    for (i = 0, y = fl2; i < iImgHeight; i ++, y ++)
    {
      for (j = 0, x = fl2; j < iImgWidth; j ++, x ++)
      {
        varInd = 0;
        k=0; 
        memset(ELocal, 0, sqrFiltLength*sizeof(Int));
        for (ii = -flV; ii < 0; ii++)
        {
          for (jj=-fl-ii; jj<=fl+ii; jj++)
          {  
            ELocal[p_pattern[k++]] += (ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
          }
        }
        Int iOffset = i * Stride + j;
        for (jj=-fl; jj<0; jj++)
        {
          ELocal[p_pattern[k++]] += (ImgDec[iOffset + jj]+ImgDec[iOffset - jj]);
        }
        ELocal[p_pattern[k++]] += ImgDec[iOffset];

        yLocal = ImgOrg[iOffset];
        m_pixAcc[varInd] += (yLocal*yLocal);
        E  = m_EGlobalSym[filtNo][varInd];
        yy = m_yGlobalSym[filtNo][varInd];

        for (k=0; k<sqrFiltLength; k++)
        {
          for (l=k; l<sqrFiltLength; l++)
          {
            E[k][l]+=(Double)(ELocal[k]*ELocal[l]);
          }
          yy[k]+=(Double)(ELocal[k]*yLocal);
        }
      }      
    }
  }

  // Matrix EGlobalSeq is symmetric, only part of it is calculated
  Double **pE = m_EGlobalSym[filtNo][varInd];
  for (k = 1; k < sqrFiltLength; k ++)
  {
    for (l = 0; l < k; l ++)
    {
      pE[k][l] = pE[l][k];
    }
  }
}
#endif

#if JVET_C0038_GALF
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////        Prediction from Fixed Filter     ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
Double TEncAdaptiveLoopFilter::xTestFixedFilterFast(Double ***A, Double **b, Double *pixAcc, Double *filterCoeffSym, Double *filterCoeffDefault, int varInd)
{

  Int j;
  Double seFilt=0, seOrig=0, error, filterCoeffTemp[m_MAX_SQR_FILT_LENGTH];
  Int sqrFiltLength= (m_MAX_SQR_FILT_LENGTH/2+1);

  for (j = 0; j < sqrFiltLength; j ++)
  {
    filterCoeffTemp[j] = filterCoeffSym[j] - filterCoeffDefault[j];
  }
  
  seOrig = pixAcc[varInd];
  seFilt = pixAcc[varInd] + calculateErrorCoeffProvided(A[varInd], b[varInd], filterCoeffTemp, sqrFiltLength);

  error = 0;
  if (seFilt < seOrig)
  {
    error += seFilt;
  }
  else
  {
    error += seOrig;
  } 
  return(error);
}
                                            
Double TEncAdaptiveLoopFilter::xTestFixedFilter(imgpel *imgY_Rec, imgpel *imgY_org, imgpel* imgY_append, Int usePrevFilt[], Int noVarBins, Int Stride, Int filtNo)
{

  Int i, j, varInd, pixelInt, temp;
  Double seFilt[m_NO_VAR_BINS]={0}, seOrig[m_NO_VAR_BINS]={0}, error;
  Int fl = m_FILTER_LENGTH/2; 
  Int offset = (1<<(m_NUM_BITS - 2));

  for (i = fl; i < m_img_height + fl; i++)
  {
    for (j = fl; j < m_img_width + fl; j++)
    {
      if(m_maskImg[i-fl][j-fl])
      {
        varInd = m_varImg[(i-fl)][(j-fl)];
        pixelInt = xFilterPixel (imgY_append, &varInd,  m_filterCoeffFinal, NULL, i, j, fl, Stride, filtNo);
        pixelInt= ((pixelInt+offset) >> (m_NUM_BITS - 1));
#if JVET_D0033_ADAPTIVE_CLIPPING
        pixelInt = ClipA(pixelInt, COMPONENT_Y); // always on luma
#else
        pixelInt = Clip3(0, m_nIBDIMax, pixelInt);

#endif
        Int iOffset     = (i-fl)*Stride + (j-fl);
        temp            = pixelInt-imgY_org[iOffset];
        seFilt[varInd] += temp*temp;

        temp            = imgY_Rec[iOffset]-imgY_org[iOffset];
        seOrig[varInd] += temp*temp;
      }
    }
  }

  error=0;
  for (i = 0; i< noVarBins; i ++)
  {
    if (seFilt[i] < seOrig[i])
    {
      usePrevFilt[i] = 1;
      error += seFilt[i];
    }
    else
    {
      usePrevFilt[i] = 0;
      error += seOrig[i];
    } 
  }
  return(error);
}

void TEncAdaptiveLoopFilter::xPreFilterFr(Int** imgY_preFilter, imgpel* imgY_rec, imgpel * imgY_org, imgpel* imgY_append, Int usePrevFilt[], Int Stride, Int filtNo)
{

  Int i, j;
  Int fl = m_FILTER_LENGTH/2;

  Int temp = 0, pixelInt = 0, offset = (1<<(m_NUM_BITS - 2));

  for (i = fl; i < m_img_height + fl; i++)
  {
    for (j = fl; j < m_img_width + fl; j++)
    {      
      Int varInd = m_varImg[(i-fl)][(j-fl)];
      Int varIndAfterMapping = selectTransposeVarInd(varInd, &temp);
      if (m_maskImg[i-fl][j-fl] && usePrevFilt[varIndAfterMapping] > 0)
      {
        pixelInt = xFilterPixel(imgY_append, &varInd, m_filterCoeffFinal, NULL, i, j, fl, Stride, filtNo);
        pixelInt= ((pixelInt+offset) >> (m_NUM_BITS - 1));        
#if JVET_D0033_ADAPTIVE_CLIPPING
        imgY_preFilter[(i-fl)][(j-fl)] = ClipA(pixelInt,COMPONENT_Y) ; // always luma
#else
        imgY_preFilter[(i-fl)][(j-fl)] = Clip3(0, m_nIBDIMax, pixelInt) ;
#endif
      }
      else
      {
        imgY_preFilter[(i-fl)][(j-fl)] = imgY_rec[(i-fl)*Stride+ (j-fl)];
      }
    }
  }
}
Void TEncAdaptiveLoopFilter::calcMatrixE(int *ELocal, imgpel *ImgDec, int *p_pattern, int i, int j, int flV, int fl, int transpose, int Stride)
{
  int ii, jj, k=0;

  if (transpose==0)
  {
    for (ii = -flV; ii < 0; ii++)
    {
      for (jj=-fl-ii; jj<=fl+ii; jj++)
      {  
        ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
      }
    }
    for (jj=-fl; jj<0; jj++)
    {
      ELocal[p_pattern[k++]]+=(ImgDec[(i)*Stride + (j+jj)]+ImgDec[(i)*Stride + (j-jj)]);
    }
    ELocal[p_pattern[k++]]+=ImgDec[(i)*Stride + (j)];
  }
  else if (transpose==1)
  {
    for (jj = -flV; jj < 0; jj++)
    {
      for (ii=-fl-jj; ii<=fl+jj; ii++)
      {  
        ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
      }
    }
    for (ii=-fl; ii<0; ii++)
    {
      ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + j]+ImgDec[(i-ii)*Stride + j]);
    }
    ELocal[p_pattern[k++]]+=ImgDec[(i)*Stride + (j)];
  }
  else if (transpose==2)
  {
    for (ii = -flV; ii < 0; ii++)
    {
      for (jj=fl+ii; jj>=-fl-ii; jj--)
      {  
        ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
      }
    }
    for (jj=-fl; jj<0; jj++)
    {
      ELocal[p_pattern[k++]]+=(ImgDec[(i)*Stride + (j+jj)]+ImgDec[(i)*Stride + (j-jj)]);
    }
    ELocal[p_pattern[k++]]+=ImgDec[(i)*Stride + (j)];
  }
  else
  {
    for (jj = -flV; jj < 0; jj++)
    {
      for (ii=fl+jj; ii>=-fl-jj; ii--)
      {  
        ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
      }
    }
    for (ii=-fl; ii<0; ii++)
    {
      ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + j]+ImgDec[(i-ii)*Stride + j]);
    }
    ELocal[p_pattern[k++]]+=ImgDec[(i)*Stride + (j)];
  }
}
Int TEncAdaptiveLoopFilter::xFilterPixel(imgpel *ImgDec, Int* varIndBeforeMapping,  Int **filterCoeffSym, Int *pattern, Int i, Int j, Int fl, Int Stride, Int filtNo)
{
  Int varInd;  
  Int transpose;
  imgpel *imgY_rec = ImgDec;
  imgpel *p_imgY_pad, *p_imgY_pad0;

  varInd = selectTransposeVarInd((*varIndBeforeMapping), &transpose);
  (*varIndBeforeMapping) = varInd;
  Int *coef = filterCoeffSym == NULL? m_filterCoeffPrevSelected[varInd]: filterCoeffSym[varInd];  
  Int pixelInt= 0;

  if(filterCoeffSym == NULL)
  {
    if(transpose==1)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[38]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[30]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[39]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[32]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[22]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[31]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[40]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[37]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[29]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[21]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[13]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);
      }
      else
      {
        pixelInt += coef[36]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[28]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[37]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[34]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[20]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[24]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[12]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);
      }
    }
    else if(transpose==3)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[38]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[32]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[39]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[30]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[22]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[31]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[40]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[37]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[33]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[23]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[13]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
      else
      {
        pixelInt += coef[36]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[34]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[37]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[28]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[24]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[14]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[4]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
    }
    else if(transpose==2)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[22]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[32]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[31]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[30]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[38]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[39]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[40]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[13]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[23]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[33]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
      else
      {
        pixelInt += coef[4]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[14]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[24]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[34]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[28]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[36]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);
      }
    }
    else
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[22]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[30]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[31]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[32]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[38]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[39]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[40]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[13]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[21]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[29]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
      else
      {
        pixelInt += coef[4]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[12]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[20]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[24]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[28]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[34]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[36]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
    }
  }
  else
  {
    //test prev filter
    if(transpose==1)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[4 ]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[1 ]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[5 ]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[3 ]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[0 ]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[2 ]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[6 ]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[9 ]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[0 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]);

      }
      else
      {
        pixelInt += coef[16]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[17]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[15]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[18]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[19]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[0 ]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl]);
      }
    }
    else if(transpose==3)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[4 ]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[3 ]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[5 ]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[1 ]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[0 ]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[2 ]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[6 ]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[9 ]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[0 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]);
      }
      else
      {
        pixelInt += coef[16]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[15]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[17]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[18]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[19]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[0 ]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl]);
      }
    }
    else if(transpose==2)
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[0 ]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[3 ]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[2 ]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[1 ]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[4 ]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[5 ]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[6 ]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[0 ]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[ 3]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[ 2]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[ 1]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]);

      }
      else
      {
        pixelInt += coef[0 ]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[15]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[16]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[17]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[18]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[19]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl]);
      }
    }
    else
    {
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[0]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[1]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[2]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[3]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[4]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[5]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[6]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[0]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[1]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[2]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[3]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]);
      }
      else
      {
        pixelInt += coef[0 ]* (imgY_rec[(i-fl+4)*Stride + j-fl]+imgY_rec[(i-fl-4)*Stride + j-fl]);
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[1 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[2 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[3 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[4 ]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[5 ]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[6 ]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[7 ]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[8 ]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[9 ]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[10]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[11]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[12]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[15]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[16]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[17]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[18]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[19]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[20]* (p_imgY_pad[j-fl]);
      }
    }
  }
  return(pixelInt);
}
#endif

Void   TEncAdaptiveLoopFilter::xstoreInBlockMatrix(imgpel* ImgOrg, imgpel* ImgDec, Int tap, Int Stride)
{
#if JVET_C0038_GALF
  if(bUpdateMatrix)
  {
    Int i,j,k,l,varInd;
#else
  Int var_step_size_w = m_ALF_VAR_SIZE_W;
  Int var_step_size_h = m_ALF_VAR_SIZE_H;
  Int i,j,k,l,varInd,ii,jj;
#endif
  Int x, y;
  Int fl =tap/2;
  Int flV = TComAdaptiveLoopFilter::ALFFlHToFlV(fl);
  Int sqrFiltLength = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(tap);
  Int fl2=9/2; //extended size at each side of the frame
  Int ELocal[m_MAX_SQR_FILT_LENGTH];
  Int yLocal;
  Int *p_pattern;
  Int filtNo =2; 
  double **E,*yy;
  Int count_valid=0;
  if (tap==9)
    filtNo =0;
  else if (tap==7)
    filtNo =1;
  
  p_pattern= m_patternTab[filtNo];
  
  memset( m_pixAcc, 0,sizeof(double)*m_NO_VAR_BINS);
#if JVET_C0038_GALF
  for (varInd=0; varInd<m_max_NO_VAR_BINS; varInd++)
#else
  for (varInd=0; varInd<m_NO_VAR_BINS; varInd++)
#endif
  {
    memset(m_yGlobalSym[filtNo][varInd], 0, sizeof(double)*m_MAX_SQR_FILT_LENGTH);
    for (k=0; k<sqrFiltLength; k++)
    {
      memset(m_EGlobalSym[filtNo][varInd][k], 0, sizeof(double)*m_MAX_SQR_FILT_LENGTH);
    }
  }
  for (i = fl2; i < m_im_height+fl2; i++)
  {
    for (j = fl2; j < m_im_width+fl2; j++)
    {
      if (m_maskImg[i-fl2][j-fl2] == 1)
      {
        count_valid++;
      }
    }
  }

  {
    for (i=0,y=fl2; i<m_im_height; i++,y++)
    {
      for (j=0,x=fl2; j<m_im_width; j++,x++)
      {
        Int condition = (m_maskImg[i][j] == 0 && count_valid > 0);
        if(!condition)
        {
#if JVET_C0038_GALF
          varInd = m_varImg[i][j];
#else
          varInd = m_varImg[i/var_step_size_h][j/var_step_size_w];
#endif
          k=0; 
          memset(ELocal, 0, sqrFiltLength*sizeof(int));
#if JVET_C0038_GALF
          Int transpose = 0;
          Int varIndMod = selectTransposeVarInd(varInd, &transpose);
          yLocal=ImgOrg[(i)*Stride + (j)]- ImgDec[(i)*Stride + (j)];
          calcMatrixE(ELocal, ImgDec, p_pattern, i, j, flV, fl, transpose, Stride);
          E= m_EGlobalSym[filtNo][varIndMod];
          yy= m_yGlobalSym[filtNo][varIndMod];

          for (k=0; k<sqrFiltLength; k++)
          {
            for (l=k; l<sqrFiltLength; l++)
            {
              E[k][l]+=(double)(ELocal[k]*ELocal[l]);
            }
            yy[k]+=(double)(ELocal[k]*yLocal);
          }
          m_pixAcc[varIndMod]+=(yLocal*yLocal);
#else
          for (ii = -flV; ii < 0; ii++)
          {
            for (jj=-fl-ii; jj<=fl+ii; jj++)
            {  
              ELocal[p_pattern[k++]]+=(ImgDec[(i+ii)*Stride + (j+jj)]+ImgDec[(i-ii)*Stride + (j-jj)]);
            }
          }
          for (jj=-fl; jj<0; jj++)
          {
            ELocal[p_pattern[k++]]+=(ImgDec[(i)*Stride + (j+jj)]+ImgDec[(i)*Stride + (j-jj)]);
          }
          ELocal[p_pattern[k++]]+=ImgDec[(i)*Stride + (j)];
          ELocal[sqrFiltLength-1]=1;
          yLocal=ImgOrg[(i)*Stride + (j)];

          m_pixAcc[varInd]+=(yLocal*yLocal);
          E= m_EGlobalSym[filtNo][varInd];
          yy= m_yGlobalSym[filtNo][varInd];

          for (k=0; k<sqrFiltLength; k++)
          {
            for (l=k; l<sqrFiltLength; l++)
              E[k][l]+=(double)(ELocal[k]*ELocal[l]);
            yy[k]+=(double)(ELocal[k]*yLocal);
          }
#endif
        }
      }
    }
  }

  // Matrix EGlobalSeq is symmetric, only part of it is calculated
  for (varInd=0; varInd<m_NO_VAR_BINS; varInd++)
  {
    double **pE = m_EGlobalSym[filtNo][varInd];
    for (k=1; k<sqrFiltLength; k++)
    {
      for (l=0; l<k; l++)
      {
        pE[k][l]=pE[l][k];
      }
    }
  }
#if JVET_C0038_GALF 
  }
  else
  {
    assert(tap != 9);
    Int filtNo = (tap == 7 ? 1: 2);
    for (Int varInd=0; varInd < m_max_NO_VAR_BINS; varInd ++)
    {
      xDeriveGlobalEyFromLgrTapFilter(m_EGlobalSym[0][varInd], m_yGlobalSym[0][varInd], m_EGlobalSym[filtNo][varInd], m_yGlobalSym[filtNo][varInd], m_patternMapTab[0], m_patternMapTab[filtNo]);
    }
  }
#endif
}

Void   TEncAdaptiveLoopFilter::xFilteringFrameLuma_qc(imgpel* ImgOrg, imgpel* imgY_pad, imgpel* ImgFilt, ALFParam* ALFp, Int tap, Int Stride, const TComSlice * pSlice)
{
#if JVET_C0038_GALF
  static Double **y_temp, **y_temp9x9;  
  static Int first = 0;
#endif

  int  filtNo,filters_per_fr;
  static double **ySym, ***ESym;
  int lambda_val = (Int) m_dLambdaLuma;
  lambda_val = lambda_val * (1<<(2*m_nBitIncrement));
  if (tap==9)
    filtNo =0;
  else if (tap==7)
    filtNo =1;
  else
    filtNo=2;
  
  ESym=m_EGlobalSym[filtNo];  
  ySym=m_yGlobalSym[filtNo];

#if JVET_C0038_GALF
  if (first==0)
  {
    initMatrix_double  (&y_temp,    m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);
    initMatrix_double  (&y_temp9x9, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);     
    first=1;
  }
  else if( !bUpdateMatrix )
  {    
    for(Int varInd = 0; varInd < m_NO_VAR_BINS; varInd ++)
    {
      xDeriveLocalEyFromLgrTapFilter(y_temp9x9[varInd], y_temp[varInd], m_patternMapTab[0], m_patternMapTab[filtNo]);
    }
  }
  xfindBestFilterVarPred(ySym, ESym, filtNo, &filters_per_fr, ImgFilt, imgY_pad, lambda_val, pSlice, ImgOrg, Stride, y_temp);

#else
  xfindBestFilterVarPred(ySym, ESym, m_pixAcc, m_filterCoeffSym, m_filterCoeffSymQuant, filtNo, &filters_per_fr,
                         m_varIndTab, NULL, m_varImg, m_maskImg, NULL, lambda_val, pSlice);
#endif
  
#if JVET_C0038_GALF
  if( bUpdateMatrix )
  {
    for(Int i = 0; i< m_NO_VAR_BINS; i ++)
    {
      for(Int j = 0; j< m_MAX_SQR_FILT_LENGTH; j++)
      {
        y_temp9x9[i][j] = y_temp[i][j];
      }
    }
  }
  xcalcPredFilterCoeff(filtNo);
#else
  // g_filterCoeffPrevSelected = g_filterCoeffSym
  xcalcPredFilterCoeff(filtNo);
#endif

#if JVET_C0038_GALF
  Int iActualFiltNo = filtNo;
  filtNo = 0; 
#endif
#if COM16_C806_ALF_TEMPPRED_NUM || JVET_C0038_GALF
  Int* patternMap=m_patternMapTab[filtNo];
  for(Int varInd=0; varInd<m_NO_VAR_BINS; ++varInd)
  {
    Int k=0;
    for(Int i = 0; i < m_MAX_SQR_FILT_LENGTH; i++)
    {
      if( patternMap[i] > 0 )
      {
        ALFp->alfCoeffLuma[varInd][k] = m_filterCoeffPrevSelected[varInd][i];
        k++;
      }
    }
  }        
#endif
  //filter the frame with g_filterCoeffPrevSelected
#if JVET_C0038_GALF
  xfilterFrame_en(imgY_pad, ImgFilt, filtNo, Stride);
  filtNo = iActualFiltNo;
  xcodeFiltCoeff(m_filterCoeffSymQuant, filtNo, m_varIndTab, filters_per_fr, 0,  ALFp, pSlice, m_predIdxTab, m_focedCoeff0, m_codedVarBins);
#else
  xfilterFrame_en(imgY_pad, ImgFilt, filtNo, Stride);
  xcodeFiltCoeff(m_filterCoeffSymQuant, filtNo, m_varIndTab, filters_per_fr,0, ALFp, pSlice);
#endif
}

Void TEncAdaptiveLoopFilter::xfilterFrame_en(imgpel* ImgDec, imgpel* ImgRest,int filtNo, int Stride)
{
  imgpel *imgY_rec = ImgDec;
#if !JVET_C0038_GALF
  imgpel *p_imgY_pad, *p_imgY_pad0;
  int var_step_size_w = m_ALF_VAR_SIZE_W;
  int var_step_size_h = m_ALF_VAR_SIZE_H;
#endif

  int i,j,y,x;
  int fl;
  int pixelInt;
  int offset = (1<<(m_NUM_BITS - 2));
#if !JVET_C0038_GALF
  int sqrFiltLength;
  sqrFiltLength=m_MAX_SQR_FILT_LENGTH;  
#endif
  fl=m_FILTER_LENGTH/2;
  
  for (y=0, i = fl; i < m_im_height+fl; i++, y++)
  {
    for (x=0, j = fl; j < m_im_width+fl; j++, x++)
    {
#if JVET_C0038_GALF
      Int varInd = m_varImg[y][x];
      pixelInt   = xFilterPixel(imgY_rec, &varInd,  NULL, NULL, i, j, fl, Stride, filtNo);
#else
      int varInd=m_varImg[y/var_step_size_h][x/var_step_size_w];
      int *coef = m_filterCoeffPrevSelected[varInd];
      pixelInt= m_filterCoeffPrevSelected[varInd][sqrFiltLength-1]; 
      if (filtNo == 2) //5x5
      {
        pixelInt += coef[22]* (imgY_rec[(i-fl+2)*Stride + j-fl]+imgY_rec[(i-fl-2)*Stride + j-fl]);

        pixelInt += coef[30]* (imgY_rec[(i-fl+1)*Stride + j-fl+1]+imgY_rec[(i-fl-1)*Stride + j-fl-1]);
        pixelInt += coef[31]* (imgY_rec[(i-fl+1)*Stride + j-fl]  +imgY_rec[(i-fl-1)*Stride + j-fl]);
        pixelInt += coef[32]* (imgY_rec[(i-fl+1)*Stride + j-fl-1]+imgY_rec[(i-fl-1)*Stride + j-fl+1]);

        pixelInt += coef[38]* (imgY_rec[(i-fl)*Stride + j-fl-2]+imgY_rec[(i-fl)*Stride + j-fl+2]);
        pixelInt += coef[39]* (imgY_rec[(i-fl)*Stride + j-fl-1]+imgY_rec[(i-fl)*Stride + j-fl+1]);
        pixelInt += coef[40]* (imgY_rec[(i-fl)*Stride + j-fl]);
      }
      else if (filtNo == 1) //7x7
      {
        pixelInt += coef[13]* (imgY_rec[(i-fl+3)*Stride + j-fl]+imgY_rec[(i-fl-3)*Stride + j-fl]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[21]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[29]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
      else
      {
        p_imgY_pad = imgY_rec + (i-fl+3)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-3)*Stride;
        pixelInt += coef[12]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[13]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[14]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);

        p_imgY_pad = imgY_rec + (i-fl+2)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-2)*Stride;
        pixelInt += coef[20]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[21]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[22]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[23]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[24]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);

        p_imgY_pad = imgY_rec + (i-fl+1)*Stride;
        p_imgY_pad0 = imgY_rec + (i-fl-1)*Stride;
        pixelInt += coef[28]* (p_imgY_pad[j-fl+3]+p_imgY_pad0[j-fl-3]);
        pixelInt += coef[29]* (p_imgY_pad[j-fl+2]+p_imgY_pad0[j-fl-2]);
        pixelInt += coef[30]* (p_imgY_pad[j-fl+1]+p_imgY_pad0[j-fl-1]);
        pixelInt += coef[31]* (p_imgY_pad[j-fl]+p_imgY_pad0[j-fl]);
        pixelInt += coef[32]* (p_imgY_pad[j-fl-1]+p_imgY_pad0[j-fl+1]);
        pixelInt += coef[33]* (p_imgY_pad[j-fl-2]+p_imgY_pad0[j-fl+2]);
        pixelInt += coef[34]* (p_imgY_pad[j-fl-3]+p_imgY_pad0[j-fl+3]);

        p_imgY_pad = imgY_rec + (i-fl)*Stride;
        pixelInt += coef[36]* (p_imgY_pad[j-fl+4]+p_imgY_pad[j-fl-4]);
        pixelInt += coef[37]* (p_imgY_pad[j-fl+3]+p_imgY_pad[j-fl-3]);
        pixelInt += coef[38]* (p_imgY_pad[j-fl+2]+p_imgY_pad[j-fl-2]);
        pixelInt += coef[39]* (p_imgY_pad[j-fl+1]+p_imgY_pad[j-fl-1]);
        pixelInt += coef[40]* (p_imgY_pad[j-fl]);

      }
#endif
      pixelInt=(int)((pixelInt+offset) >> (m_NUM_BITS - 1));
#if JVET_D0033_ADAPTIVE_CLIPPING
      ImgRest[y*Stride + x] = ClipA(pixelInt, COMPONENT_Y); // always luma here
#else
      ImgRest[y*Stride + x] = Clip3(0, m_nIBDIMax, pixelInt);
#endif
    }
  }
}
#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::xfindBestFilterPredictor(Double ***E_temp, Double**y_temp, Double *pixAcc_temp, Int filtNo, const TComSlice * pSlice
                                     ,imgpel* ImgOrg, imgpel* ImgDec, Int Stride, Bool* forceCoeff0, Double errorForce0CoeffTab[m_NO_VAR_BINS][2], Char* usePrevFiltBest, Bool*  codedVarBins, Int sqrFiltLength, Int fl )
{
  Int    varInd, i, j, k, yLocal, filterNo;
  Int    ELocal[m_MAX_SQR_FILT_LENGTH];
  Int    usePrevFilt[m_NO_VAR_BINS];  
  
  Double errorMin, error;

  Int noVarBins = m_max_NO_VAR_BINS;
  Int factor = (1<<(TComAdaptiveLoopFilter::m_NUM_BITS-1));
  Char iFixedFilters = m_tempALFp->iAvailableFilters;

  memset(ELocal, 0, sqrFiltLength*sizeof(Int));

  if(bTapDecision == false)
  {
    memset(usePrevFiltBest, 0, m_NO_VAR_BINS*sizeof(Char));
    for(varInd=0; varInd<noVarBins; ++varInd)
    {
      for(i = 0; i < (m_MAX_SQR_FILT_LENGTH/2+1); i++)
      {
        m_filterCoeffFinal[varInd][i] = (Int)(m_filterCoeffDefault[i]*factor);
      }
    }
  }
  else
  {
    if( bFindBestFixedFilter == false || bUpdateMatrix) 
    {

      for(varInd=0; varInd<noVarBins; ++varInd)
      {
        errorMin=0;
        for (i=0; i<iFixedFilters; i++)
        {
          filterNo=varInd*iFixedFilters + i;
          error = xTestFixedFilterFast(E_temp, y_temp, pixAcc_temp, m_filterCoeffPrev[filterNo], m_filterCoeffDefault, varInd);

          if (error<errorMin || i==0)
          {
            errorMin=error;
            usePrevFiltBest[varInd] = (Char)i;
          }
        }
      }

      for(varInd=0; varInd<noVarBins; ++varInd)
      {
        for(i = 0; i < (m_MAX_SQR_FILT_LENGTH/2+1); i++)
        {
          filterNo = (Int)(varInd * iFixedFilters + usePrevFiltBest[varInd]);
          m_filterCoeffFinal[varInd][i] = (Int) (m_filterCoeffPrev[filterNo][i]* factor);
        } 
      }

      error = xTestFixedFilter(ImgDec, ImgOrg, ImgDec, usePrevFilt, noVarBins, Stride, filtNo);
      xPreFilterFr(m_imgY_preFilter, ImgDec, ImgOrg, ImgDec, usePrevFilt, Stride, filtNo);
      
      for(varInd=0; varInd<noVarBins; ++varInd)
      {
        if (usePrevFilt[varInd]==1)
        {
          usePrevFiltBest[varInd]++;
        }
        else
        {
          usePrevFiltBest[varInd]=0;
          for(i = 0; i < TComAdaptiveLoopFilter::m_MAX_SQT_FILT_SYM_LENGTH; i++)
          {
            m_filterCoeffFinal[varInd][i] = (Int)(m_filterCoeffDefault[i]*factor);
          }
        }
      }
      bFindBestFixedFilter = true;
    }
  }
  
  // If frNo>0 pixAcc and yGlobalSym have to be calculated again
  if(bTapDecision == true && m_tempALFp->iAvailableFilters)
  {
    if( bUpdateMatrix )
    {
      for (varInd=0; varInd<noVarBins; varInd++)
      {
        pixAcc_temp[varInd]=0;
        for (k=0; k<sqrFiltLength; k++)
        {
          y_temp[varInd][k]=0;
        }
      }
      Int transpose;
      Int flV = TComAdaptiveLoopFilter::ALFFlHToFlV(fl);
      for (i = fl; i < m_im_height+fl; i++)
      {
        for (j = fl; j < m_im_width+fl; j++)
        {
          if(m_maskImg[i-fl][j-fl])
          {
            memset(ELocal, 0, sqrFiltLength*sizeof(Int));
            varInd  = selectTransposeVarInd(m_varImg[i-fl][j-fl], &transpose);
            Int pos = (i-fl)*Stride + (j-fl);
            yLocal=ImgOrg[pos]- m_imgY_preFilter[i-fl][j-fl];
            calcMatrixE (ELocal, ImgDec, m_patternTab[filtNo], i-fl, j-fl, flV, fl, transpose, Stride);

            for (k=0; k<sqrFiltLength; k++)
            {
              y_temp[varInd][k]+=(Double)(ELocal[k]*yLocal);
            }
            pixAcc_temp[varInd]+=(yLocal*yLocal);
          }
        }
      }
    }
  }
}

Void TEncAdaptiveLoopFilter::xfindBestFilterVarPred(Double **ySym, Double ***ESym, 
                              Int filtNo, Int *filters_per_fr_best, imgpel *imgY_rec, 
                              imgpel *ImgDec, Double lambda_val, const TComSlice * pSlice, imgpel* ImgOrg, Int Stride, Double **y_temp)
{
  Int filters_per_fr, firstFilt, varIndTab[m_NO_VAR_BINS][m_NO_VAR_BINS], sqrFiltLength, numBits, coeffBits, createBistream;
  Int i, k, varInd;
  Int *weights;
 
  Double  error, lambda, lagrangian, lagrangianMin; 
    
  static Int first = 0;
  static Double ***E_temp, *pixAcc_temp;
  static Int **FilterCoeffQuantTemp;
 
  lambda = lambda_val;
  sqrFiltLength=m_MAX_SQR_FILT_LENGTH;
  
  if (first == 0)
  {
    initMatrix3D_double(&E_temp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);
    pixAcc_temp = (Double *) calloc(m_NO_VAR_BINS, sizeof(Double));
    initMatrix_int(&FilterCoeffQuantTemp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);
    first=1;
  }
  
  sqrFiltLength = m_sqrFiltLengthTab[filtNo];   
  Int fl  = m_flTab[filtNo];

  weights = m_weightsTab[filtNo];                 
  memcpy(pixAcc_temp, m_pixAcc, sizeof(Double)*m_NO_VAR_BINS);

  for (varInd=0; varInd<m_max_NO_VAR_BINS; varInd++)
  {
    if(bUpdateMatrix)    
    {
      memcpy(y_temp[varInd],ySym[varInd],sizeof(Double)*sqrFiltLength);
    }
    for (k=0; k<sqrFiltLength; k++)
    {
      memcpy(E_temp[varInd][k],ESym[varInd][k],sizeof(Double)*sqrFiltLength);
    }
  }
 
  // zero all variables 
  memset(varIndTab, 0, sizeof(Int)*m_NO_VAR_BINS);

  Bool forceCoeff0, codedVarBins[m_NO_VAR_BINS];
  Char iFixedFilters = m_tempALFp->iAvailableFilters;
  static Char usePrevFiltBest[m_NO_VAR_BINS];
  Double errorForce0CoeffTab[m_NO_VAR_BINS][2];

  xfindBestFilterPredictor(E_temp, y_temp, pixAcc_temp, filtNo, pSlice, ImgOrg, ImgDec, Stride, &forceCoeff0,
                           errorForce0CoeffTab, usePrevFiltBest, codedVarBins, sqrFiltLength, fl);

  for(i = 0; i < m_max_NO_VAR_BINS; i++)
  {
    memset(m_filterCoeffSym[i], 0, sizeof(Int)*m_MAX_SQR_FILT_LENGTH);
    memset(m_filterCoeffSymQuant[i], 0, sizeof(Int)*m_MAX_SQR_FILT_LENGTH);
  }
  
  firstFilt = 1;  
  lagrangianMin = 0;
  filters_per_fr = m_max_NO_FILTERS;

  xFindFilterGroupingError(E_temp, y_temp, pixAcc_temp, varIndTab, sqrFiltLength, m_max_NO_VAR_BINS, m_max_NO_VAR_BINS);

  while(filters_per_fr>=1)
  {
    error = findFilterCoeff(E_temp, y_temp, pixAcc_temp, 
                    varIndTab[filters_per_fr-1], m_max_NO_VAR_BINS, sqrFiltLength, filters_per_fr, weights, numBits=m_NUM_BITS,errorForce0CoeffTab);
    lagrangian = xfindBestCoeffCodMethod(codedVarBins, &forceCoeff0, fl, sqrFiltLength, filters_per_fr, errorForce0CoeffTab, &error, lambda, pSlice);

    if (lagrangian<lagrangianMin || firstFilt==1)
    {
      firstFilt = 0;
      lagrangianMin=lagrangian;
      (*filters_per_fr_best) = filters_per_fr;
      memcpy(m_varIndTab, varIndTab[filters_per_fr - 1], m_NO_VAR_BINS*sizeof(Int));
      m_focedCoeff0 = forceCoeff0;
      memcpy(m_codedVarBins, codedVarBins, m_NO_VAR_BINS*sizeof(Bool));
    }

    filters_per_fr--;

  }


  error = findFilterCoeff(E_temp, y_temp, pixAcc_temp, m_varIndTab,
                  m_max_NO_VAR_BINS, sqrFiltLength, (*filters_per_fr_best), weights, numBits=m_NUM_BITS
                   ,errorForce0CoeffTab );  

  xfindBestCoeffCodMethod(m_codedVarBins, &forceCoeff0, fl, sqrFiltLength, 
    (*filters_per_fr_best), errorForce0CoeffTab, &error, lambda, pSlice);

  coeffBits = xcodeAuxInfo(filtNo, (*filters_per_fr_best), m_varIndTab, createBistream = 0,  m_tempALFp, pSlice, usePrevFiltBest, iFixedFilters - 1, m_max_NO_VAR_BINS, forceCoeff0);

  if (forceCoeff0 == false)
  {
    coeffBits += xsendAllFiltersPPPred(m_filterCoeffSymQuant, fl, sqrFiltLength, (*filters_per_fr_best), createBistream = 0, m_tempALFp, pSlice);
  }
  else
  {
    coeffBits += xsendAllFiltersPPPredForce0(m_filterCoeffSymQuant, fl, sqrFiltLength, 
      (*filters_per_fr_best), m_codedVarBins, createBistream=0, m_tempALFp, pSlice);
    for(varInd=0; varInd < (*filters_per_fr_best); varInd++)
    {
      if (m_codedVarBins[varInd] == false)
      {
        memset(m_filterCoeffSym[varInd], 0, sizeof(Int)*m_MAX_SQR_FILT_LENGTH);
        memset(m_filterCoeffSymQuant[varInd], 0, sizeof(Int)*m_MAX_SQR_FILT_LENGTH);
      }
    }
  }
  m_focedCoeff0 = forceCoeff0;
  memcpy( m_predIdxTab, usePrevFiltBest, sizeof(Char)*m_NO_VAR_BINS);  
}
#else
Void TEncAdaptiveLoopFilter::xfindBestFilterVarPred(double **ySym, double ***ESym, double *pixAcc, int **filterCoeffSym, int **filterCoeffSymQuant, int filtNo, int *filters_per_fr_best, int varIndTab[], imgpel **imgY_rec, imgpel **varImg, imgpel **maskImg, imgpel **imgY_pad, double lambda_val, const TComSlice * pSlice)
{
  int filters_per_fr, firstFilt, forceCoeff0,
  interval[m_NO_VAR_BINS][2], intervalBest[m_NO_VAR_BINS][2];
  int i, k, varInd;
  static double ***E_temp, **y_temp, *pixAcc_temp;
  static int **FilterCoeffQuantTemp;
  double  error, lambda, lagrangian, lagrangianMin;
  
  int sqrFiltLength;
  int *weights;
  int numBits, coeffBits;
  double errorForce0CoeffTab[m_NO_VAR_BINS][2];
  int  codedVarBins[m_NO_VAR_BINS], createBistream /*, forceCoeff0 */;
  //int  usePrevFilt[m_NO_VAR_BINS];
  static int first=0;
  
  //for (i = 0; i < m_NO_VAR_BINS; i++)
  //  usePrevFiltDefault[i]=usePrevFilt[i]=1;
  lambda = lambda_val;
  sqrFiltLength=m_MAX_SQR_FILT_LENGTH;
  
  if (first==0)
  {
    initMatrix3D_double(&E_temp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);
    initMatrix_double(&y_temp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH); 
    pixAcc_temp = (double *) calloc(m_NO_VAR_BINS, sizeof(double));
    initMatrix_int(&FilterCoeffQuantTemp, m_NO_VAR_BINS, m_MAX_SQR_FILT_LENGTH);
    first=1;
  }
  
  sqrFiltLength=m_sqrFiltLengthTab[filtNo];   
  Int fl = m_flTab[filtNo];
  weights=m_weightsTab[filtNo];               
  //patternMap=m_patternMapTab[filtNo];  
  //pattern=m_patternTab[filtNo];
  
  memcpy(pixAcc_temp,pixAcc,sizeof(double)*m_NO_VAR_BINS);

  for (varInd=0; varInd<m_NO_VAR_BINS; varInd++)
  {
    memcpy(y_temp[varInd],ySym[varInd],sizeof(double)*sqrFiltLength);
    for (k=0; k<sqrFiltLength; k++)
      memcpy(E_temp[varInd][k],ESym[varInd][k],sizeof(double)*sqrFiltLength);
  }
  
  // zero all variables 
  memset(varIndTab,0,sizeof(int)*m_NO_VAR_BINS);
  
  for(i = 0; i < m_NO_VAR_BINS; i++)
  {
    memset(filterCoeffSym[i],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
    memset(filterCoeffSymQuant[i],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
  }
  
  firstFilt=1;  lagrangianMin=0;
  filters_per_fr=m_NO_FILTERS;

  while(filters_per_fr>=1)
  {
    findFilterGroupingError(E_temp, y_temp, pixAcc_temp, interval, sqrFiltLength, filters_per_fr);
    findFilterCoeff(E_temp, y_temp, pixAcc_temp, filterCoeffSym, filterCoeffSymQuant, interval,
                    varIndTab, sqrFiltLength, filters_per_fr, weights, numBits=m_NUM_BITS,  errorForce0CoeffTab);
    lagrangian=xfindBestCoeffCodMethod(codedVarBins, &forceCoeff0, filterCoeffSymQuant, fl, 
                                       sqrFiltLength, filters_per_fr, errorForce0CoeffTab, &error, lambda, pSlice);

    if (lagrangian<lagrangianMin || firstFilt==1)
    {
      firstFilt=0;
      lagrangianMin=lagrangian;

      (*filters_per_fr_best)=filters_per_fr;
      memcpy(intervalBest, interval, m_NO_VAR_BINS*2*sizeof(int));
    }
    filters_per_fr--;
  }


  findFilterCoeff(E_temp, y_temp, pixAcc_temp, filterCoeffSym, filterCoeffSymQuant, intervalBest,
                  varIndTab, sqrFiltLength, (*filters_per_fr_best), weights, numBits=m_NUM_BITS, errorForce0CoeffTab);
  
  xfindBestCoeffCodMethod(codedVarBins, &forceCoeff0, filterCoeffSymQuant, fl, sqrFiltLength, 
                          (*filters_per_fr_best), errorForce0CoeffTab, &error, lambda, pSlice);

  
  //coded=1;
  if (forceCoeff0==1 && (*filters_per_fr_best)==1)
  {
    //coded=0;
    coeffBits = xcodeAuxInfo(-1, (*filters_per_fr_best), varIndTab, 0, createBistream=0,filtNo, m_tempALFp, pSlice);
  }
  else
  {
    coeffBits = xcodeAuxInfo(filtNo, (*filters_per_fr_best), varIndTab, 0, createBistream=0,filtNo, m_tempALFp, pSlice);
  }

  if (forceCoeff0==0)
  {
    coeffBits += xsendAllFiltersPPPred(filterCoeffSymQuant, fl, sqrFiltLength, 
      (*filters_per_fr_best), createBistream=0, m_tempALFp, pSlice);
  }
  else
  {
    if ((*filters_per_fr_best)==1)
    {
      for(varInd=0; varInd<(*filters_per_fr_best); varInd++)
      {
        memset(filterCoeffSym[varInd],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
        memset(filterCoeffSymQuant[varInd],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
      }
    }
    else
    {
      coeffBits += xsendAllFiltersPPPredForce0(filterCoeffSymQuant, fl, sqrFiltLength, 
        (*filters_per_fr_best), codedVarBins, createBistream=0, m_tempALFp, pSlice);

      for(varInd=0; varInd<(*filters_per_fr_best); varInd++)
      {
        if (codedVarBins[varInd]==0)
        {
          memset(filterCoeffSym[varInd],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
          memset(filterCoeffSymQuant[varInd],0,sizeof(int)*m_MAX_SQR_FILT_LENGTH);
        }
      }
    }
  }
}
#endif

Void TEncAdaptiveLoopFilter::xcalcPredFilterCoeff(int filtNo)
{
  int *patternMap, varInd, i, k;

  patternMap=m_patternMapTab[filtNo];

#if JVET_C0038_GALF
  Int quantCoeffSum = 0;
  Int * weights=TComAdaptiveLoopFilter::m_weightsTab[0]; 
  Int factor = (1<<(TComAdaptiveLoopFilter::m_NUM_BITS-1)); 
  for(varInd=0; varInd<m_max_NO_VAR_BINS; ++varInd)
#else
  for(varInd=0; varInd<m_NO_VAR_BINS; ++varInd)
#endif
  {
    k=0;
#if JVET_C0038_GALF
    quantCoeffSum = 0;
#endif
    for(i = 0; i < m_MAX_SQR_FILT_LENGTH; i++)
    {
#if JVET_C0038_GALF
     if (m_pattern9x9Sym_Quart[i]>0 && patternMap[i]>0)
     {
       m_filterCoeffPrevSelected[varInd][i] = (m_filterCoeffFinal[varInd][m_pattern9x9Sym_Quart[i]-1]+m_filterCoeffSym[m_varIndTab[varInd]][patternMap[i]-1]);

       if(i!= (m_MAX_SQR_FILT_LENGTH-1))
       {
         quantCoeffSum += (weights[m_pattern9x9Sym_Quart[i]-1]* m_filterCoeffPrevSelected[varInd][i]);
       }
       k++;
     }
     else if( m_pattern9x9Sym_Quart[i] > 0 )
     {
       m_filterCoeffPrevSelected[varInd][i] = m_filterCoeffFinal[varInd][m_pattern9x9Sym_Quart[i]-1];
       if(i!= (m_MAX_SQR_FILT_LENGTH-1))
       {
         quantCoeffSum += (weights[m_pattern9x9Sym_Quart[i]-1]* m_filterCoeffPrevSelected[varInd][i]);
       }
     }
     else
     {
       m_filterCoeffPrevSelected[varInd][i]= 0;
     }
#else
      if (patternMap[i]>0)
      {
        m_filterCoeffPrevSelected[varInd][i]=m_filterCoeffSym[m_varIndTab[varInd]][k];
        k++;
      }
      else
      {
        m_filterCoeffPrevSelected[varInd][i]=0;
      }
#endif
    }
#if JVET_C0038_GALF
    m_filterCoeffPrevSelected[varInd][m_MAX_SQR_FILT_LENGTH - 1] = factor-quantCoeffSum;
#endif
  }
}

#if COM16_C806_ALF_TEMPPRED_NUM
Void TEncAdaptiveLoopFilter::xcopyFilterCoeff(int filtNo, int **filterCoeff)
{
  int *patternMap, varInd, i, k;

  patternMap=m_patternMapTab[filtNo];
  for(varInd=0; varInd<m_NO_VAR_BINS; ++varInd)
  {
    k=0;
    for(i = 0; i < m_MAX_SQR_FILT_LENGTH; i++)
    {
      if (patternMap[i]>0)
      {
        m_filterCoeffPrevSelected[varInd][i]=filterCoeff[varInd][k];
        k++;
      }
      else
      {
        m_filterCoeffPrevSelected[varInd][i]=0;
      }
    }
  }
}
#endif

UInt TEncAdaptiveLoopFilter::xcodeFiltCoeff(int **filterCoeffSymQuant, int filtNo, int varIndTab[], int filters_per_fr_best, int frNo, ALFParam* ALFp, const TComSlice * pSlice
#if JVET_C0038_GALF
  , Char PrevFiltIdx[]
  , Bool forcedCoeff0
  , Bool codedVarBins[]
#endif
  )
{
#if JVET_C0038_GALF 
  Int coeffBits, createBistream;
  Int sqrFiltLength = m_sqrFiltLengthTab[filtNo], fl=m_flTab[filtNo];
  
  ALFp->filters_per_group_diff = filters_per_fr_best;
  ALFp->filters_per_group = filters_per_fr_best;  

  ALFp->forceCoeff0 = forcedCoeff0;
  if(forcedCoeff0)
  {
    memcpy(ALFp->codedVarBins, codedVarBins, sizeof(Bool)*ALF_NO_VAR_BIN);
  }
  else
  {
    for(Int ind = 0; ind < ALF_NO_VAR_BIN; ind++)
    {
      ALFp->codedVarBins[ind] = true;
    }
  }
  coeffBits = xcodeAuxInfo(filtNo, filters_per_fr_best, varIndTab, createBistream=1, ALFp, pSlice, PrevFiltIdx, ALFp->iAvailableFilters - 1, m_max_NO_VAR_BINS, forcedCoeff0);
  ALFp->predMethod = 0;
  ALFp->num_coeff = sqrFiltLength;
  ALFp->realfiltNo=filtNo;
  if (filters_per_fr_best <= 1)
  {
    ALFp->predMethod = 0;
  } 

  if ( !forcedCoeff0 )
  {
    coeffBits += xsendAllFiltersPPPred(filterCoeffSymQuant, fl, sqrFiltLength, filters_per_fr_best, createBistream=1, ALFp, pSlice);
  }
  else
  {
    coeffBits += xsendAllFiltersPPPredForce0(filterCoeffSymQuant, fl, sqrFiltLength, filters_per_fr_best, codedVarBins, createBistream = 1, ALFp, pSlice);
  }
  return (UInt)coeffBits;
#else
  int varInd, forceCoeff0, codedVarBins[m_NO_VAR_BINS], coeffBits, createBistream,   sqrFiltLength=m_sqrFiltLengthTab[filtNo], 
  fl=m_flTab[filtNo];
  
  ALFp->filters_per_group_diff = filters_per_fr_best;
  ALFp->filters_per_group = filters_per_fr_best;
  
  for(varInd=0; varInd<filters_per_fr_best; varInd++)
  {
    codedVarBins[varInd] = 1;
  }
  memcpy (ALFp->codedVarBins, codedVarBins, sizeof(int)*m_NO_VAR_BINS);
  forceCoeff0=0;
  for(varInd=0; varInd<filters_per_fr_best; varInd++)
  {
    if (codedVarBins[varInd] == 0)
    {
      forceCoeff0=1;
      break;
    }
  }

  //coded=1;
  if (forceCoeff0==1 && filters_per_fr_best==1)
  {
    //coded=0;
    coeffBits = xcodeAuxInfo(-1, filters_per_fr_best, varIndTab, frNo, createBistream=1,filtNo, ALFp, pSlice);
  }
  else
  {
    coeffBits = xcodeAuxInfo(filtNo, filters_per_fr_best, varIndTab, frNo, createBistream=1,filtNo, ALFp, pSlice);
  }
  ALFp->forceCoeff0 = forceCoeff0;
  ALFp->predMethod = 0;
  ALFp->num_coeff = sqrFiltLength;
  ALFp->realfiltNo=filtNo;
  if (filters_per_fr_best <= 1)
  {
    ALFp->forceCoeff0 = 0;
    ALFp->predMethod = 0;
  }
  if (forceCoeff0==0) 
  {
    coeffBits += xsendAllFiltersPPPred(filterCoeffSymQuant, fl, sqrFiltLength, 
                                       filters_per_fr_best, createBistream=1, ALFp, pSlice);
  }
  else if (filters_per_fr_best>1)
  {
    coeffBits += xsendAllFiltersPPPredForce0(filterCoeffSymQuant, fl, sqrFiltLength, 
                                             filters_per_fr_best, codedVarBins, createBistream=1, ALFp, pSlice);
  }
  return (UInt)coeffBits;
#endif
}



Void TEncAdaptiveLoopFilter::xSetCUAlfCtrlFlags_qc(UInt uiAlfCtrlDepth, TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiDist, ALFParam *pAlfParam)
{
  ruiDist = 0;
  pAlfParam->num_alf_cu_flag = 0;
  for( UInt uiCUAddr = 0; uiCUAddr < m_pcPic->getNumberOfCtusInFrame() ; uiCUAddr++ )
  {
    TComDataCU* pcCU = m_pcPic->getCtu( uiCUAddr );
#if JVET_C0024_QTBT
    xSetCUAlfCtrlFlag_qc(pcCU, 0, 0, pcCU->getSlice()->getSPS()->getCTUSize(), pcCU->getSlice()->getSPS()->getCTUSize(), uiAlfCtrlDepth, pcPicOrg, pcPicDec, pcPicRest, ruiDist, pAlfParam);
#else
    xSetCUAlfCtrlFlag_qc(pcCU, 0, 0, uiAlfCtrlDepth, pcPicOrg, pcPicDec, pcPicRest, ruiDist, pAlfParam);
#endif
  }
}

#if JVET_C0024_QTBT
Void TEncAdaptiveLoopFilter::xSetCUAlfCtrlFlag_qc(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiAlfCtrlDepth, TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiDist, ALFParam *pAlfParam)
#else
Void TEncAdaptiveLoopFilter::xSetCUAlfCtrlFlag_qc(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiAlfCtrlDepth, TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiDist, ALFParam *pAlfParam)
#endif
{
#if !JVET_C0024_QTBT
  Bool bBoundary = false;
#endif
  UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiRPelX   = uiLPelX + uiWidth  - 1;
#else
  UInt uiRPelX   = uiLPelX + (pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiDepth)  - 1;
#endif
  UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  UInt uiBPelY   = uiTPelY + uiHeight - 1;
#else
  UInt uiBPelY   = uiTPelY + (pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiDepth) - 1;
#endif
  
#if !JVET_C0024_QTBT
  if( ( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) || ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
  {
    bBoundary = true;
  }
#endif

#if JVET_C0024_QTBT
#if JVET_C0024_QTBT //alf, add BT structure?  JCA 
  if( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) && uiDepth < uiAlfCtrlDepth))
#else
  if( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) || bBoundary )
#endif
#else
  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < pcCU->getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
#endif
  {
    UInt uiQNumParts = ( m_pcPic->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
      
#if BARC_Boundary
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if( ( uiLPelX + deltaX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY + deltaY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
#else
      if( ( uiLPelX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
#endif
#if JVET_C0024_QTBT
        xSetCUAlfCtrlFlag_qc(pcCU, uiAbsPartIdx, uiDepth+1, uiWidth>>1, uiHeight>>1, uiAlfCtrlDepth, pcPicOrg, pcPicDec, pcPicRest, ruiDist, pAlfParam);
#else
        xSetCUAlfCtrlFlag_qc(pcCU, uiAbsPartIdx, uiDepth+1, uiAlfCtrlDepth, pcPicOrg, pcPicDec, pcPicRest, ruiDist, pAlfParam);
#endif
    }
    return;
  }
  
#if JVET_C0024_QTBT
  UInt uiWidthPel = uiWidth;
  UInt uiHeightPel = uiHeight;

#if BARC_Boundary
  UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX =  (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
  UInt deltaY =  (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());

  if( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - deltaX)
  {
    //bBoundary = true;
    uiWidthPel = pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - deltaX - uiLPelX;
    uiRPelX = pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - deltaX - 1;
  }
  if ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - deltaY )
  {
    //bBoundary = true;
    uiHeightPel = pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - deltaY - uiTPelY;
    uiBPelY = pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - deltaY - 1;
  }
#else
  if( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() )
  {
    //bBoundary = true;
    uiWidthPel = pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - uiLPelX;
    uiRPelX = pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - 1;
  }
  if ( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() )
  {
      //bBoundary = true;
      uiHeightPel = pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - uiTPelY;
      uiBPelY = pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - 1;
  }
#endif
#else
  if( uiDepth > uiAlfCtrlDepth && !pcCU->isFirstAbsZorderIdxInDepth(uiAbsPartIdx, uiAlfCtrlDepth))
  {
    return;
  }
#endif
  
  UInt uiCUAddr = pcCU->getCtuRsAddr();
  UInt64 uiRecSSD = 0;
  UInt64 uiFiltSSD = 0;
  
#if !JVET_C0024_QTBT
  Int iWidth;
  Int iHeight;
  UInt uiSetDepth;
  
  if (uiDepth > uiAlfCtrlDepth && pcCU->isFirstAbsZorderIdxInDepth(uiAbsPartIdx, uiAlfCtrlDepth))
  {
#if JVET_C0024_QTBT
    iWidth = pcCU->getSlice()->getSPS()->getCTUSize() >> uiAlfCtrlDepth;
    iHeight = pcCU->getSlice()->getSPS()->getCTUSize() >> uiAlfCtrlDepth;
#else
    iWidth = pcCU->getSlice()->getSPS()->getMaxCUWidth() >> uiAlfCtrlDepth;
    iHeight = pcCU->getSlice()->getSPS()->getMaxCUHeight() >> uiAlfCtrlDepth;
#endif
    
    uiRPelX   = uiLPelX + iWidth  - 1;
    uiBPelY   = uiTPelY + iHeight - 1;

    if( uiRPelX >= pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() )
    {
      iWidth = pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() - uiLPelX;
    }
    
    if( uiBPelY >= pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() )
    {
      iHeight = pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() - uiTPelY;
    }
    
    uiSetDepth = uiAlfCtrlDepth;
  }
  else
  {
    iWidth = pcCU->getWidth(uiAbsPartIdx);
    iHeight = pcCU->getHeight(uiAbsPartIdx);
    uiSetDepth = uiDepth;
  }
#endif
  
  Pel* pOrg = pcPicOrg->getAddr(COMPONENT_Y, uiCUAddr, uiAbsPartIdx);
  Pel* pRec = pcPicDec->getAddr(COMPONENT_Y, uiCUAddr, uiAbsPartIdx);
  Pel* pFilt = pcPicRest->getAddr(COMPONENT_Y, uiCUAddr, uiAbsPartIdx);
  
#if JVET_C0024_QTBT
  uiRecSSD  += xCalcSSD( pOrg, pRec,  uiWidthPel, uiHeightPel, pcPicOrg->getStride(COMPONENT_Y) );
  uiFiltSSD += xCalcSSD( pOrg, pFilt, uiWidthPel, uiHeightPel, pcPicOrg->getStride(COMPONENT_Y) );
#else
  uiRecSSD  += xCalcSSD( pOrg, pRec,  iWidth, iHeight, pcPicOrg->getStride(COMPONENT_Y) );
  uiFiltSSD += xCalcSSD( pOrg, pFilt, iWidth, iHeight, pcPicOrg->getStride(COMPONENT_Y) );
#endif
  
  if (uiFiltSSD < uiRecSSD)
  {
    ruiDist += uiFiltSSD;
#if JVET_C0024_QTBT
    pcCU->setAlfCtrlFlagSubParts(1, uiAbsPartIdx, uiWidth, uiHeight);
#else
    pcCU->setAlfCtrlFlagSubParts(1, uiAbsPartIdx, uiSetDepth);
#endif
    pAlfParam->alf_cu_flag[pAlfParam->num_alf_cu_flag]=1;
    for (int i=uiTPelY ;i<=min(uiBPelY,(unsigned int)(pcPicOrg->getHeight(COMPONENT_Y)-1))  ;i++)
    {
      for (int j=uiLPelX ;j<=min(uiRPelX,(unsigned int)(pcPicOrg->getWidth(COMPONENT_Y)-1)) ;j++)
      { 
        m_maskImg[i][j]=1;
      }
    }
  }
  else
  {
    ruiDist += uiRecSSD;
#if JVET_C0024_QTBT
    pcCU->setAlfCtrlFlagSubParts(0, uiAbsPartIdx, uiWidth, uiHeight);
#else
    pcCU->setAlfCtrlFlagSubParts(0, uiAbsPartIdx, uiSetDepth);
#endif
    pAlfParam->alf_cu_flag[pAlfParam->num_alf_cu_flag]=0;
    for (int i=uiTPelY ;i<=min(uiBPelY,(unsigned int)(pcPicOrg->getHeight(COMPONENT_Y)-1))  ;i++)
    {
      for (int j=uiLPelX ;j<=min(uiRPelX,(unsigned int)(pcPicOrg->getWidth(COMPONENT_Y)-1)) ;j++)
      { 
        m_maskImg[i][j]=0;
      }
    }
  }
  pAlfParam->num_alf_cu_flag++;
}

Void TEncAdaptiveLoopFilter::xReDesignFilterCoeff_qc(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Bool bReadCorr, const TComSlice * pSlice)
{
  
  Int tap = m_pcTempAlfParam->tap;
  Int    LumaStride = pcPicOrg->getStride(COMPONENT_Y);
  imgpel* pOrg = (imgpel*)pcPicOrg->getAddr(COMPONENT_Y);
  imgpel* pDec = (imgpel*)pcPicDec->getAddr(COMPONENT_Y);
  imgpel* pRest = (imgpel*)pcPicRest->getAddr(COMPONENT_Y);
  xFirstFilteringFrameLuma(pOrg, pDec, pRest, m_pcTempAlfParam, tap, LumaStride, pSlice); 
  
}
Void TEncAdaptiveLoopFilter::xCUAdaptiveControl_qc(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiMinRate, UInt64& ruiMinDist, Double& rdMinCost, const TComSlice * pSlice)
{
  m_pcEntropyCoder->setAlfCtrl(true);
  
  UInt uiBestDepth = 0;
  
  ALFParam cFrmAlfParam;
  allocALFParam(&cFrmAlfParam);
  copyALFParam(&cFrmAlfParam, m_pcBestAlfParam);
#if JVET_C0038_GALF
  Int    Height = pcPicOrg->getHeight(COMPONENT_Y);
  Int    Width = pcPicOrg->getWidth(COMPONENT_Y);
#endif  

  for (UInt uiDepth = 0; uiDepth < pSlice->getSPS()->getMaxTotalCUDepth(); uiDepth++)
  {
#if JVET_C0024_QTBT
    Int nBlkSize = ( pSlice->getSPS()->getCTUSize() * pSlice->getSPS()->getCTUSize() ) >> ( uiDepth << 1 );
#else
    Int nBlkSize = ( pSlice->getSPS()->getMaxCUHeight() * pSlice->getSPS()->getMaxCUWidth() ) >> ( uiDepth << 1 );
#endif
    Int nPicSize = pcPicOrg->getWidth(COMPONENT_Y) * pcPicOrg->getHeight(COMPONENT_Y);
    if( ( nBlkSize << 4 ) > nPicSize )
      // block is too large
      continue;
    m_pcEntropyCoder->setMaxAlfCtrlDepth(uiDepth);
    pcPicRest->copyToPic(m_pcPicYuvTmp, COMPONENT_Y);
    copyALFParam(m_pcTempAlfParam, &cFrmAlfParam);
    m_pcTempAlfParam->cu_control_flag = 1;   
    for (UInt uiRD = 0; uiRD <= ALF_NUM_OF_REDESIGN; uiRD++)
    {
      if (uiRD)
      {
        // re-design filter coefficients
        xReDesignFilterCoeff_qc(pcPicOrg, pcPicDec, m_pcPicYuvTmp, true, pSlice); //use filtering of mine
      }
      
      UInt64 uiRate, uiDist;
      Double dCost;
     //m_pcPicYuvTmp: filtered signal, pcPicDec: orig reconst
      xSetCUAlfCtrlFlags_qc(uiDepth, pcPicOrg, pcPicDec, m_pcPicYuvTmp, uiDist, m_pcTempAlfParam); //set up varImg here
      
      xCalcRDCost(m_pcTempAlfParam, uiRate, uiDist, dCost, pSlice);
      
      if (dCost < rdMinCost)
      {
        uiBestDepth = uiDepth;
        rdMinCost = dCost;
        ruiMinDist = uiDist;
        ruiMinRate = uiRate;
        m_pcPicYuvTmp->copyToPic(m_pcPicYuvBest,COMPONENT_Y);
        copyALFParam(m_pcBestAlfParam, m_pcTempAlfParam);
        //save maskImg
        xCopyTmpAlfCtrlFlagsFrom();
#if JVET_C0038_GALF
        for (Int i=0; i<Height; i++)
        {
          for (Int j=0; j<Width; j++)
          {
            m_maskBestpImg[i][j] = m_maskImg[i][j];
          }
        }
#endif    
      }
    }
  }
  
  if (m_pcBestAlfParam->cu_control_flag)
  {
    m_pcEntropyCoder->setAlfCtrl(true);
    m_pcEntropyCoder->setMaxAlfCtrlDepth(uiBestDepth);
    xCopyTmpAlfCtrlFlagsTo();
    m_pcPicYuvBest->copyToPic(pcPicRest, COMPONENT_Y);//copy m_pcPicYuvBest to pcPicRest
    xCopyDecToRestCUs(pcPicDec, pcPicRest); //pcPicRest = pcPicDec
#if JVET_C0038_GALF    
    for (Int i=0; i<Height; i++)
    {
      for (Int j=0; j<Width; j++)
      {
        m_maskImg[i][j] = m_maskBestpImg[i][j];
      }
    }    
#endif 
  }
  else
  {
    m_pcEntropyCoder->setAlfCtrl(false);
    m_pcEntropyCoder->setMaxAlfCtrlDepth(0);
  }
  freeALFParam(&cFrmAlfParam);
}


Void TEncAdaptiveLoopFilter::xFilterTapDecision_qc(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiMinRate, UInt64& ruiMinDist, Double& rdMinCost, const TComSlice * pSlice)
{
  // restriction for non-referenced B-slice
  if (m_eSliceType == B_SLICE && m_iPicNalReferenceIdc == 0)
  {
    return;
  }
  
  UInt64 uiRate, uiDist;
  Double dCost;
  
  if (m_pcBestAlfParam->cu_control_flag)
  {
    xCopyTmpAlfCtrlFlagsFrom();
  }
  
  Bool bChanged = false;
#if JVET_C0038_GALF
  for (Int iTap = m_ALF_MAX_NUM_TAP; iTap >= m_ALF_MIN_NUM_TAP; iTap -= 2)
  {
    if(m_eSliceType==I_SLICE && iTap == m_ALF_MAX_NUM_TAP)
    {
      continue;
    }
    if(iTap != m_ALF_MAX_NUM_TAP && !m_pcTempAlfParam->cu_control_flag)
    {
      bUpdateMatrix = false;
    }
#else
  for (Int iTap = m_ALF_MIN_NUM_TAP; iTap <= m_ALF_MAX_NUM_TAP; iTap += 2)
  {
#endif
    copyALFParam(m_pcTempAlfParam, m_pcBestAlfParam);
    m_pcTempAlfParam->tap = iTap;
    m_pcTempAlfParam->tapV = TComAdaptiveLoopFilter::ALFTapHToTapV(m_pcTempAlfParam->tap);
    m_pcTempAlfParam->num_coeff = TComAdaptiveLoopFilter::ALFTapHToNumCoeff(m_pcTempAlfParam->tap); 
    
    if (m_pcTempAlfParam->cu_control_flag)
    {
      xReDesignFilterCoeff_qc(pcPicOrg, pcPicDec, m_pcPicYuvTmp, false, pSlice);
      xSetCUAlfCtrlFlags_qc(m_pcEntropyCoder->getMaxAlfCtrlDepth(), pcPicOrg, pcPicDec, m_pcPicYuvTmp, uiDist, m_pcTempAlfParam);
      xCalcRDCost(m_pcTempAlfParam, uiRate, uiDist, dCost, pSlice);
    }

    else
    {
      Int    Height = pcPicOrg->getHeight(COMPONENT_Y);
      Int    Width = pcPicOrg->getWidth(COMPONENT_Y);

      for (Int i=0; i<Height; i++)
      {
        for (Int j=0; j<Width; j++)
        {
          m_maskImg[i][j] = 1;
        }
      }
      xReDesignFilterCoeff_qc(pcPicOrg, pcPicDec, m_pcPicYuvTmp, false, pSlice);

      xCalcRDCost(pcPicOrg, m_pcPicYuvTmp, m_pcTempAlfParam, uiRate, uiDist, dCost, pSlice);
    }

    if (dCost < rdMinCost)
    {
      rdMinCost = dCost;
      ruiMinDist = uiDist;
      ruiMinRate = uiRate;
      m_pcPicYuvTmp->copyToPic(m_pcPicYuvBest, COMPONENT_Y);
      copyALFParam(m_pcBestAlfParam, m_pcTempAlfParam);
      bChanged = true;
      if (m_pcTempAlfParam->cu_control_flag)
      {
        xCopyTmpAlfCtrlFlagsFrom();
      }
    }
  }
  
  if (m_pcBestAlfParam->cu_control_flag)
  {
    xCopyTmpAlfCtrlFlagsTo();
    if (bChanged)
    {
      m_pcPicYuvBest->copyToPic(pcPicRest, COMPONENT_Y);
      xCopyDecToRestCUs(pcPicDec, pcPicRest);
    }
  }
  else if (m_pcBestAlfParam->tap >= m_ALF_MIN_NUM_TAP && bChanged ) // fixed a bug with HM-3
  {
    m_pcPicYuvBest->copyToPic(pcPicRest, COMPONENT_Y);
  }
  
  copyALFParam(m_pcTempAlfParam, m_pcBestAlfParam);
}


#define ROUND(a)  (((a) < 0)? (int)((a) - 0.5) : (int)((a) + 0.5))
#define REG              0.0001
#define REG_SQR          0.0000001

//Find filter coeff related
Int TEncAdaptiveLoopFilter::gnsCholeskyDec(double **inpMatr, double outMatr[m_MAX_SQR_FILT_LENGTH][m_MAX_SQR_FILT_LENGTH], int noEq)
{ 
  int 
  i, j, k;     /* Looping Variables */
  double 
  scale;       /* scaling factor for each row */
  double 
  invDiag[m_MAX_SQR_FILT_LENGTH];  /* Vector of the inverse of diagonal entries of outMatr */
  
  
  /*
   *  Cholesky decomposition starts
   */
  
  for(i = 0; i < noEq; i++)
  {
    for(j = i; j < noEq; j++)
    {
      /* Compute the scaling factor */
      scale=inpMatr[i][j];
      if ( i > 0) for( k = i - 1 ; k >= 0 ; k--)
        scale -= outMatr[k][j] * outMatr[k][i];
      
      /* Compute i'th row of outMatr */
      if(i==j)
      {
        if(scale <= REG_SQR ) // if(scale <= 0 )  /* If inpMatr is singular */
        {
          return(0);
        }
        else              /* Normal operation */
          invDiag[i] =  1.0/(outMatr[i][i]=sqrt(scale));
      }
      else
      {
        outMatr[i][j] = scale*invDiag[i]; /* Upper triangular part          */
        outMatr[j][i] = 0.0;              /* Lower triangular part set to 0 */
      }                    
    }
  }
  return(1); /* Signal that Cholesky factorization is successfully performed */
}


Void TEncAdaptiveLoopFilter::gnsTransposeBacksubstitution(double U[m_MAX_SQR_FILT_LENGTH][m_MAX_SQR_FILT_LENGTH], double rhs[], double x[], int order)
{
  int 
  i,j;              /* Looping variables */
  double 
  sum;              /* Holds backsubstitution from already handled rows */
  
  /* Backsubstitution starts */

  x[0] = U[0][0] != 0 ? rhs[0]/U[0][0] : 0;      /* First row of U'                   */
  for (i = 1; i < order; i++)
  {         /* For the rows 1..order-1           */
    
    for (j = 0, sum = 0.0; j < i; j++) /* Backsubst already solved unknowns */
      sum += x[j]*U[j][i];
    
    x[i]=(rhs[i] - sum)/U[i][i];       /* i'th component of solution vect.  */
  }
}



Void  TEncAdaptiveLoopFilter::gnsBacksubstitution(double R[m_MAX_SQR_FILT_LENGTH][m_MAX_SQR_FILT_LENGTH], double z[m_MAX_SQR_FILT_LENGTH], int R_size, double A[m_MAX_SQR_FILT_LENGTH])
{
  int i, j;
  double sum;
  
  R_size--;
  
  A[R_size] = z[R_size] / R[R_size][R_size];
  
  for (i = R_size-1; i >= 0; i--)
  {
    for (j = i+1, sum = 0.0; j <= R_size; j++)
      sum += R[i][j] * A[j];
    
    A[i] = (z[i] - sum) / R[i][i];
  }
}


Int TEncAdaptiveLoopFilter::gnsSolveByChol(double **LHS, double *rhs, double *x, int noEq)
{ 
  double aux[m_MAX_SQR_FILT_LENGTH];     /* Auxiliary vector */
  double U[m_MAX_SQR_FILT_LENGTH][m_MAX_SQR_FILT_LENGTH];    /* Upper triangular Cholesky factor of LHS */
  int  i, singular;          /* Looping variable */
  
  /* The equation to be solved is LHSx = rhs */
  
  /* Compute upper triangular U such that U'*U = LHS */
  if(gnsCholeskyDec(LHS, U, noEq)) /* If Cholesky decomposition has been successful */
  {
    singular=1;
    /* Now, the equation is  U'*U*x = rhs, where U is upper triangular
     * Solve U'*aux = rhs for aux
     */
    gnsTransposeBacksubstitution(U, rhs, aux, noEq);         
    
    /* The equation is now U*x = aux, solve it for x (new motion coefficients) */
    gnsBacksubstitution(U, aux, noEq, x);   
    
  }
  else /* LHS was singular */ 
  {
    singular=0;
    
    /* Regularize LHS */
    for(i=0; i<noEq; i++)
      LHS[i][i] += REG;
    /* Compute upper triangular U such that U'*U = regularized LHS */
    singular = gnsCholeskyDec(LHS, U, noEq);
    /* Solve  U'*aux = rhs for aux */  
    gnsTransposeBacksubstitution(U, rhs, aux, noEq);   
    
    /* Solve U*x = aux for x */
    gnsBacksubstitution(U, aux, noEq, x);
  }  
  return(singular);
}


//////////////////////////////////////////////////////////////////////////////////////////

#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::add_A(Double **Amerged, Double ***A, Int varIndTab[], Int filtNo, Int noVarBins, Int size)
{ 
  Int i, j, ind;          

  for (i = 0; i < size; i ++)
  {
    for (j = 0; j < size; j ++)
    {
      Amerged[i][j] = 0;
      for (ind = 0; ind < noVarBins; ind ++)
      {
        if (varIndTab[ind] == filtNo)
        {
          Amerged[i][j] += A[ind][i][j];
        }
      }
    }
  }
}

Void TEncAdaptiveLoopFilter::add_b(Double *bmerged, Double **b, Int varIndTab[], Int filtNo, Int noVarBins, Int size)
{ 
  Int i, ind;        

  for (i = 0; i < size; i ++)
  {
    bmerged[i] = 0;
    for (ind = 0; ind < noVarBins; ind ++)
    {
      if (varIndTab[ind] == filtNo)
      {
        bmerged[i] += b[ind][i];
      }
    }
  }
}
#else
Void TEncAdaptiveLoopFilter::add_A(double **Amerged, double ***A, int start, int stop, int size)
{ 
  int
  i, j, ind;          /* Looping variable */
  
  for (i=0; i<size; i++)
  {
    for (j=0; j<size; j++)
    {
      Amerged[i][j]=0;
      for (ind=start; ind<=stop; ind++)
      {
        Amerged[i][j]+=A[ind][i][j];
      }
    }
  }
}

Void TEncAdaptiveLoopFilter::add_b(double *bmerged, double **b, int start, int stop, int size)
{ 
  int
  i, ind;          /* Looping variable */
  
  for (i=0; i<size; i++)
  {
    bmerged[i]=0;
    for (ind=start; ind<=stop; ind++)
    {
      bmerged[i]+=b[ind][i];
    }
  }
}
#endif

double TEncAdaptiveLoopFilter::calculateErrorCoeffProvided(double **A, double *b, double *c, int size)
{
  int i, j;
  double error, sum=0;
  
  error=0;
  for (i=0; i<size; i++)   //diagonal
  {
    sum=0;
    for (j=i+1; j<size; j++)
      sum+=(A[j][i]+A[i][j])*c[j];
    error+=(A[i][i]*c[i]+sum-2*b[i])*c[i];
  }
  
  return(error);
}

double TEncAdaptiveLoopFilter::calculateErrorAbs(double **A, double *b, double y, int size)
{
  int i;
  double error, sum;
  double c[m_MAX_SQR_FILT_LENGTH];
  
  gnsSolveByChol(A, b, c, size);
  
  sum=0;
  for (i=0; i<size; i++)
  {
    sum+=c[i]*b[i];
  }
  error=y-sum;
  
  return(error);
}
#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::xMergeFiltersGreedy(Double **yGlobalSeq, Double ***EGlobalSeq, Double *pixAccGlobalSeq, Int indexListBest[m_NO_VAR_BINS][m_NO_VAR_BINS], Int sqrFiltLength, Int noVarBins, Int maxNoIntervals)
{
  Int first, ind, ind1, ind2, noRemaining, i, j, exist, indexList[m_NO_VAR_BINS], indexListTemp[m_NO_VAR_BINS], available[m_NO_VAR_BINS], bestToMerge[2];
  Double error, error1, error2, errorMin;
  static Double *y_temp, **E_temp, pixAcc_temp;
  static Int init = 0;
  if (init == 0)
  {
    initMatrix_double(&E_temp, m_MAX_SQR_FILT_LENGTH, m_MAX_SQR_FILT_LENGTH);
    y_temp = new Double[m_MAX_SQR_FILT_LENGTH];
    init = 1;
  }

  noRemaining = noVarBins;  
  for (ind=0; ind<m_max_NO_VAR_BINS; ind++)
  {
    indexListBest[noRemaining-1][ind] = ind;
    indexList[ind]                    = ind; 
    available[ind]                    = 1;
    m_pixAcc_merged[ind] = pixAccGlobalSeq[ind];

    memcpy(m_y_merged[ind], yGlobalSeq[ind], sizeof(Double)*sqrFiltLength);
    for (i = 0; i < sqrFiltLength; i ++)
    {
      memcpy(m_E_merged[ind][i], EGlobalSeq[ind][i], sizeof(Double)*sqrFiltLength);
    }
  }
  for (ind=0; ind<m_max_NO_VAR_BINS; ind++)
  {
    indexListBest[0][ind] = 0;
  }

  // Try merging different matrices

  while (noRemaining > 2)
  {
    errorMin = 0; first = 1; bestToMerge[0] = 0; bestToMerge[1] = 1;
    for (ind1=0; ind1<noVarBins-1; ind1++)
    {
      if (available[ind1])
      {
        for (ind2=ind1+1; ind2<noVarBins; ind2++)
        {
          if (available[ind2])
          {
            error1=calculateErrorAbs(m_E_merged[ind1], m_y_merged[ind1], m_pixAcc_merged[ind1], sqrFiltLength);
            error2=calculateErrorAbs(m_E_merged[ind2], m_y_merged[ind2], m_pixAcc_merged[ind2], sqrFiltLength);

            pixAcc_temp=m_pixAcc_merged[ind1]+m_pixAcc_merged[ind2];
            for (i=0; i<sqrFiltLength; i++){
              y_temp[i]=m_y_merged[ind1][i]+m_y_merged[ind2][i];
              for (j=0; j<sqrFiltLength; j++){
                E_temp[i][j]=m_E_merged[ind1][i][j]+m_E_merged[ind2][i][j];
              }
            }
            error = calculateErrorAbs(E_temp, y_temp, pixAcc_temp, sqrFiltLength)-error1-error2;

            if (error<errorMin || first==1)
            {
              errorMin=error;
              bestToMerge[0]=ind1;
              bestToMerge[1]=ind2;
              first=0;
            }
          }
        }
      }
    }

    ind1 = bestToMerge[0];
    ind2 = bestToMerge[1];

    m_pixAcc_merged[ind1] += m_pixAcc_merged[ind2];
    for (i=0; i<sqrFiltLength; i++)
    {
      m_y_merged[ind1][i]+=m_y_merged[ind2][i];
      for (j=0; j<sqrFiltLength; j++)
      {
        m_E_merged[ind1][i][j]+=m_E_merged[ind2][i][j];
      }
    }

    available[ind2]=0;

    for (i=0; i<noVarBins; i++)
    {
      if (indexList[i]==ind2)
      {
        indexList[i]=ind1;
      }
    }
    noRemaining--;
    if (noRemaining<=maxNoIntervals)
    {
      for (i=0; i<noVarBins; i++)
      {
        indexListTemp[i]=indexList[i];
      }

      exist=0; ind=0;
      for (j=0; j<noVarBins; j++)
      {
        exist=0;
        for (i=0; i<noVarBins; i++)
        {
          if (indexListTemp[i]==j)
          {
            exist=1;
            break;
          }
        }

        if (exist)
        {
          for (i=0; i<noVarBins; i++)
          {
            if (indexListTemp[i]==j)
            {
              indexListBest[noRemaining-1][i]=ind;
              indexListTemp[i]=-1;
            }
          }
          ind++;
        }
      }
    }
  }
}
#else
double TEncAdaptiveLoopFilter::mergeFiltersGreedy(double **yGlobalSeq, double ***EGlobalSeq, double *pixAccGlobalSeq, int intervalBest[m_NO_VAR_BINS][2], int sqrFiltLength, int noIntervals)
{
  int first, ind, ind1, ind2, i, j, bestToMerge ;
  double error, error1, error2, errorMin;
  static double pixAcc_temp, error_tab[m_NO_VAR_BINS],error_comb_tab[m_NO_VAR_BINS];
  static int indexList[m_NO_VAR_BINS], available[m_NO_VAR_BINS], noRemaining;
  if (noIntervals == m_NO_FILTERS)
  {
    noRemaining=m_NO_VAR_BINS;
    for (ind=0; ind<m_NO_VAR_BINS; ind++)
    {
      indexList[ind]=ind; 
      available[ind]=1;
      m_pixAcc_merged[ind]=pixAccGlobalSeq[ind];
      memcpy(m_y_merged[ind],yGlobalSeq[ind],sizeof(double)*sqrFiltLength);
      for (i=0; i<sqrFiltLength; i++)
      {
        memcpy(m_E_merged[ind][i],EGlobalSeq[ind][i],sizeof(double)*sqrFiltLength);
      }
    }
  }
  // Try merging different matrices
  if (noIntervals == m_NO_FILTERS)
  {
    for (ind=0; ind<m_NO_VAR_BINS; ind++)
    {
      error_tab[ind]=calculateErrorAbs(m_E_merged[ind], m_y_merged[ind], m_pixAcc_merged[ind], sqrFiltLength);
    }
    for (ind=0; ind<m_NO_VAR_BINS-1; ind++)
    {
      ind1=indexList[ind];
      ind2=indexList[ind+1];
      
      error1=error_tab[ind1];
      error2=error_tab[ind2];
      
      pixAcc_temp=m_pixAcc_merged[ind1]+m_pixAcc_merged[ind2];
      for (i=0; i<sqrFiltLength; i++)
      {
        m_y_temp[i]=m_y_merged[ind1][i]+m_y_merged[ind2][i];
        for (j=0; j<sqrFiltLength; j++)
        {
          m_E_temp[i][j]=m_E_merged[ind1][i][j]+m_E_merged[ind2][i][j];
        }
      }
      error_comb_tab[ind1]=calculateErrorAbs(m_E_temp, m_y_temp, pixAcc_temp, sqrFiltLength)-error1-error2;
    }
  }
  while (noRemaining>noIntervals)
  {
    errorMin=0; first=1;
    bestToMerge = 0;
    for (ind=0; ind<noRemaining-1; ind++)
    {
      error = error_comb_tab[indexList[ind]];
      if ((error<errorMin || first==1))
      {
        errorMin=error;
        bestToMerge=ind;
        first=0;
      }
    }
    ind1=indexList[bestToMerge];
    ind2=indexList[bestToMerge+1];
    m_pixAcc_merged[ind1]+=m_pixAcc_merged[ind2];
    for (i=0; i<sqrFiltLength; i++)
    {
      m_y_merged[ind1][i]+=m_y_merged[ind2][i];
      for (j=0; j<sqrFiltLength; j++)
      {
        m_E_merged[ind1][i][j]+=m_E_merged[ind2][i][j];
      }
    }
    available[ind2]=0;
    
    //update error tables
    error_tab[ind1]=error_comb_tab[ind1]+error_tab[ind1]+error_tab[ind2];
    if (indexList[bestToMerge] > 0)
    {
      ind1=indexList[bestToMerge-1];
      ind2=indexList[bestToMerge];
      error1=error_tab[ind1];
      error2=error_tab[ind2];
      pixAcc_temp=m_pixAcc_merged[ind1]+m_pixAcc_merged[ind2];
      for (i=0; i<sqrFiltLength; i++)
      {
        m_y_temp[i]=m_y_merged[ind1][i]+m_y_merged[ind2][i];
        for (j=0; j<sqrFiltLength; j++)
        {
          m_E_temp[i][j]=m_E_merged[ind1][i][j]+m_E_merged[ind2][i][j];
        }
      }
      error_comb_tab[ind1]=calculateErrorAbs(m_E_temp, m_y_temp, pixAcc_temp, sqrFiltLength)-error1-error2;
    }
    if (indexList[bestToMerge+1] < m_NO_VAR_BINS-1)
    {
      ind1=indexList[bestToMerge];
      ind2=indexList[bestToMerge+2];
      error1=error_tab[ind1];
      error2=error_tab[ind2];
      pixAcc_temp=m_pixAcc_merged[ind1]+m_pixAcc_merged[ind2];
      for (i=0; i<sqrFiltLength; i++)
      {
        m_y_temp[i]=m_y_merged[ind1][i]+m_y_merged[ind2][i];
        for (j=0; j<sqrFiltLength; j++)
        {
          m_E_temp[i][j]=m_E_merged[ind1][i][j]+m_E_merged[ind2][i][j];
        }
      }
      error_comb_tab[ind1]=calculateErrorAbs(m_E_temp, m_y_temp, pixAcc_temp, sqrFiltLength)-error1-error2;
    }
    
    ind=0;
    for (i=0; i<m_NO_VAR_BINS; i++)
    {
      if (available[i]==1)
      {
        indexList[ind]=i;
        ind++;
      }
    }
    noRemaining--;
  }
  
  
  errorMin=0;
  for (ind=0; ind<noIntervals; ind++)
  {
    errorMin+=error_tab[indexList[ind]];
  }
  
  for (ind=0; ind<noIntervals-1; ind++)
  {
    intervalBest[ind][0]=indexList[ind]; intervalBest[ind][1]=indexList[ind+1]-1;
  }
  
  intervalBest[noIntervals-1][0]=indexList[noIntervals-1]; 
  intervalBest[noIntervals-1][1]=m_NO_VAR_BINS-1;

  
  return(errorMin);
}
#endif

#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::xFindFilterGroupingError(Double ***EGlobalSeq, Double **yGlobalSeq, Double *pixAccGlobalSeq, Int intervalBest[m_NO_VAR_BINS][m_NO_VAR_BINS], Int sqrFiltLength, Int filters_per_fr, Int noVarBins)
{
  return (xMergeFiltersGreedy(yGlobalSeq, EGlobalSeq, pixAccGlobalSeq, intervalBest, sqrFiltLength, noVarBins, filters_per_fr));
}
#else
double TEncAdaptiveLoopFilter::findFilterGroupingError(double ***EGlobalSeq, double **yGlobalSeq, double *pixAccGlobalSeq, int intervalBest[m_NO_VAR_BINS][2], int sqrFiltLength, int filters_per_fr)
{
  double error;
  
  // find best filters for each frame group
  error = 0;
  error += mergeFiltersGreedy(yGlobalSeq, EGlobalSeq, pixAccGlobalSeq, intervalBest, sqrFiltLength, filters_per_fr);
  
  return(error);
}
#endif


Void TEncAdaptiveLoopFilter::roundFiltCoeff(int *FilterCoeffQuan, double *FilterCoeff, int sqrFiltLength, int factor)
{
  int i;
  double diff; 
  int diffInt, sign; 
  
  for(i = 0; i < sqrFiltLength; i++)
  {
    sign               = (FilterCoeff[i]>0) ?  1: -1; 
    diff               = FilterCoeff[i]*sign; 
    diffInt            = (int)(diff*(double)factor+0.5); 
    FilterCoeffQuan[i] = diffInt*sign;
  }
}

Double TEncAdaptiveLoopFilter::QuantizeIntegerFilterPP(double *filterCoeff, int *filterCoeffQuant, double **E, double *y, int sqrFiltLength, int *weights, int bit_depth
#if JVET_C0038_GALF
  ,Bool bChroma 
#endif
  )
{
  double error;
  
  int factor = (1<<(bit_depth-1)), i; 
  int quantCoeffSum, minInd, targetCoeffSumInt, k, diff;
#if !JVET_C0038_GALF
  double targetCoeffSum, errMin;
#else
  double errMin;
#endif
  
  gnsSolveByChol(E, y, filterCoeff, sqrFiltLength);

#if JVET_C0038_GALF
   roundFiltCoeff(filterCoeffQuant, filterCoeff, sqrFiltLength, factor);
   if(bChroma)
   {
     targetCoeffSumInt = factor;
   }
   else
   targetCoeffSumInt=0;
#else
  targetCoeffSum=0;
  for (i=0; i<sqrFiltLength; i++)
  {
    targetCoeffSum+=(weights[i]*filterCoeff[i]*factor);
  }
  targetCoeffSumInt=ROUND(targetCoeffSum);
  roundFiltCoeff(filterCoeffQuant, filterCoeff, sqrFiltLength, factor);
#endif
  quantCoeffSum=0;
  for (i=0; i<sqrFiltLength; i++)
  {
    quantCoeffSum+=weights[i]*filterCoeffQuant[i];
  }
  
  int count=0;
  while(quantCoeffSum!=targetCoeffSumInt && count < 10)
  {
    if (quantCoeffSum>targetCoeffSumInt)
    {
      diff=quantCoeffSum-targetCoeffSumInt;
      errMin=0; minInd=-1;
      for (k=0; k<sqrFiltLength; k++)
      {
        if (weights[k]<=diff)
        {
          for (i=0; i<sqrFiltLength; i++)
          {
            m_filterCoeffQuantMod[i]=filterCoeffQuant[i];
          }
          m_filterCoeffQuantMod[k]--;
          for (i=0; i<sqrFiltLength; i++)
          {
            filterCoeff[i]=(double)m_filterCoeffQuantMod[i]/(double)factor;
          }
          error=calculateErrorCoeffProvided(E, y, filterCoeff, sqrFiltLength);
          if (error<errMin || minInd==-1)
          {
            errMin=error;
            minInd=k;
          }
        } // if (weights(k)<=diff){
      } // for (k=0; k<sqrFiltLength; k++){
      filterCoeffQuant[minInd]--;
    }
    else
    {
      diff=targetCoeffSumInt-quantCoeffSum;
      errMin=0; minInd=-1;
      for (k=0; k<sqrFiltLength; k++)
      {
        if (weights[k]<=diff)
        {
          for (i=0; i<sqrFiltLength; i++)
          {
            m_filterCoeffQuantMod[i]=filterCoeffQuant[i];
          }
          m_filterCoeffQuantMod[k]++;
          for (i=0; i<sqrFiltLength; i++)
          {
            filterCoeff[i]=(double)m_filterCoeffQuantMod[i]/(double)factor;
          }
          error=calculateErrorCoeffProvided(E, y, filterCoeff, sqrFiltLength);
          if (error<errMin || minInd==-1)
          {
            errMin=error;
            minInd=k;
          }
        } // if (weights(k)<=diff){
      } // for (k=0; k<sqrFiltLength; k++){
      filterCoeffQuant[minInd]++;
    }
    
    quantCoeffSum=0;
    for (i=0; i<sqrFiltLength; i++)
    {
      quantCoeffSum+=weights[i]*filterCoeffQuant[i];
    }
  }
  if( count == 10 )
  {
    for (i=0; i<sqrFiltLength; i++)
    {
      filterCoeffQuant[i] = 0;
    }
  }
 
#if JVET_C0038_GALF
  quantCoeffSum=0;
  for (i=0; i<sqrFiltLength-1; i++){
    quantCoeffSum+=weights[i]*filterCoeffQuant[i];
  }
  if(bChroma)
  {
    filterCoeffQuant[sqrFiltLength-1]=factor-quantCoeffSum;
  }
  else
  {
    filterCoeffQuant[sqrFiltLength-1]=-quantCoeffSum;
  }
#endif


  Int max_value = std::min(  (1<<(3+m_nInputBitDepth + m_nBitIncrement))-1, (1<<14)-1);
  Int min_value = std::max( -(1<<(3+m_nInputBitDepth + m_nBitIncrement)),  -(1<<14)  );
  for (i=0; i<sqrFiltLength; i++)
  {
    filterCoeffQuant[i] = min( max_value , max( min_value , filterCoeffQuant[i] ) );
    filterCoeff[i]=(double)filterCoeffQuant[i]/(double)factor;
  }
 
#if !JVET_C0038_GALF 
  // Encoder-side restriction on the ALF coefficients for limiting the Clipping LUT size
  Int sumCoef[2] = {0, 0};
  Int maxPxlVal = m_nIBDIMax;
  Int maxSampleValue, minSampleValue;
  Int numBitsMinus1= m_NUM_BITS-1;
  Int offset = (1<<(m_NUM_BITS-2));
  Int lastCoef = sqrFiltLength-1;
  Int centerCoef = sqrFiltLength-2;
  Int *coef = filterCoeffQuant;
  Int clipTableMax = ( ( m_ALF_HM3_QC_CLIP_RANGE - m_ALF_HM3_QC_CLIP_OFFSET - 1 ) << m_nBitIncrement );
  Int clipTableMin = ( ( - m_ALF_HM3_QC_CLIP_OFFSET ) << m_nBitIncrement );

  for (i=0; i<centerCoef; i++)
  {
    sumCoef[coef[i]>0?0:1] += (coef[i]<<1);
  }
  sumCoef[ coef[centerCoef] > 0 ? 0 : 1 ] += coef[centerCoef];
  maxSampleValue = ( maxPxlVal * sumCoef[0] + coef[lastCoef] + offset ) >> numBitsMinus1;
  minSampleValue = ( maxPxlVal * sumCoef[1] + coef[lastCoef] + offset ) >> numBitsMinus1;

  if ( maxSampleValue > clipTableMax || minSampleValue < clipTableMin )
  {
    memset( coef, 0, sizeof(Int) * sqrFiltLength );
    coef[centerCoef] = (1<<numBitsMinus1);
    for (Int j=0; j<sqrFiltLength; j++)
    {
      assert( filterCoeffQuant[j]>=min_value && filterCoeffQuant[j]<=max_value );
      filterCoeff[j]=(double)filterCoeffQuant[j]/(double)factor;
    }
  }
#endif
  error=calculateErrorCoeffProvided(E, y, filterCoeff, sqrFiltLength);
  return(error);
}
#if JVET_C0038_GALF
Double TEncAdaptiveLoopFilter::findFilterCoeff(Double ***EGlobalSeq, Double **yGlobalSeq, Double *pixAccGlobalSeq, Int varIndTab[m_NO_VAR_BINS], Int noVarBins, Int sqrFiltLength, Int filters_per_fr, Int *weights, Int bit_depth, Double errorTabForce0Coeff[m_NO_VAR_BINS][2])
#else
Double TEncAdaptiveLoopFilter::findFilterCoeff(double ***EGlobalSeq, double **yGlobalSeq, double *pixAccGlobalSeq, int **filterCoeffSeq, int **filterCoeffQuantSeq, int intervalBest[m_NO_VAR_BINS][2], int varIndTab[m_NO_VAR_BINS], int sqrFiltLength, int filters_per_fr, int *weights, int bit_depth, double errorTabForce0Coeff[m_NO_VAR_BINS][2])
#endif
{
  static double pixAcc_temp;
  double error;
  int k, filtNo;
  
  error = 0;
  for(filtNo = 0; filtNo < filters_per_fr; filtNo++)
  {
#if JVET_C0038_GALF
    add_A(m_E_temp, EGlobalSeq, varIndTab, filtNo, noVarBins, sqrFiltLength);
    add_b(m_y_temp, yGlobalSeq, varIndTab, filtNo, noVarBins, sqrFiltLength);

    pixAcc_temp=0;    
    for (k=0; k<noVarBins; k++)
    {
      if (varIndTab[k]==filtNo)
      {
        pixAcc_temp+=pixAccGlobalSeq[k];
      }
    }
#else
    add_A(m_E_temp, EGlobalSeq, intervalBest[filtNo][0], intervalBest[filtNo][1], sqrFiltLength);
    add_b(m_y_temp, yGlobalSeq, intervalBest[filtNo][0], intervalBest[filtNo][1], sqrFiltLength);
    
    pixAcc_temp = 0;    
    for(k = intervalBest[filtNo][0]; k <= intervalBest[filtNo][1]; k++)
      pixAcc_temp += pixAccGlobalSeq[k];
#endif  
    // Find coeffcients
    errorTabForce0Coeff[filtNo][1] = pixAcc_temp + QuantizeIntegerFilterPP(m_filterCoeff, m_filterCoeffQuant, m_E_temp, m_y_temp, sqrFiltLength, weights, bit_depth);
    errorTabForce0Coeff[filtNo][0] = pixAcc_temp;
    error += errorTabForce0Coeff[filtNo][1];
    for(k = 0; k < sqrFiltLength; k++)
    {
#if JVET_C0038_GALF
      m_filterCoeffSym[filtNo][k] = m_filterCoeffQuant[k];
      m_filterCoeffSymQuant[filtNo][k] = m_filterCoeffQuant[k];
#else
      filterCoeffSeq[filtNo][k] = m_filterCoeffQuant[k];
      filterCoeffQuantSeq[filtNo][k] = m_filterCoeffQuant[k];
#endif
    }
  }
#if !JVET_C0038_GALF 
  for(filtNo = 0; filtNo < filters_per_fr; filtNo++)
  {
    for(k = intervalBest[filtNo][0]; k <= intervalBest[filtNo][1]; k++)
      varIndTab[k] = filtNo;
  }
#endif
  return(error);
}

#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::setInitialMask(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, ALFParam* pcAlfParam )
#else
Void TEncAdaptiveLoopFilter::setInitialMask(TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec)
#endif
{
  Int Height = pcPicOrg->getHeight(COMPONENT_Y);
  Int Width = pcPicOrg->getWidth(COMPONENT_Y);
  Int LumaStride = pcPicOrg->getStride(COMPONENT_Y);
  imgpel* pDec = (imgpel*)pcPicDec->getAddr(COMPONENT_Y);



#if JVET_C0038_GALF
  calcVar(m_varImg, pDec, 9/2, JVET_C0038_SHIFT_VAL_HALFW, Height, Width, LumaStride, 0, 0);
#else
  calcVar(m_varImg, pDec, 9/2, m_VAR_SIZE, Height, Width, LumaStride);
#endif


  for(Int y=0; y<Height; y++)
  {
    for(Int x=0; x<Width; x++)
    {
      m_maskImg[y][x] = 1;
    }
  }

}


#if COM16_C806_ALF_TEMPPRED_NUM
Bool TEncAdaptiveLoopFilter::xFilteringLumaChroma( ALFParam *pAlfParam, TComPicYuv* pcPicOrg, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, UInt64& ruiMinRate, UInt64& ruiMinDist, Double& rdMinCost, Int uiIndex, const TComSlice * pSlice )
{
  UInt64 uiRate, uiDist = 0;
  Double dCost;

  imgpel* pDec = (imgpel*)pcPicDec->getAddr(COMPONENT_Y);
  imgpel* pRest = (imgpel*)pcPicRest->getAddr(COMPONENT_Y);
  Int    LumaStride = pcPicOrg->getStride(COMPONENT_Y);
  UInt   uiTmpMaxDepth = m_pcEntropyCoder->getMaxAlfCtrlDepth();
  UInt   uiTmpAlfCtrlFlag = m_pcEntropyCoder->getAlfCtrl();
  m_pcTempAlfParam->temproalPredFlag = true;
  copyALFParam(m_pcTempAlfParam, pAlfParam);
  xcopyFilterCoeff( m_pcTempAlfParam->filtNo, m_pcTempAlfParam->alfCoeffLuma );
  m_pcTempAlfParam->cu_control_flag = 0;
  m_pcTempAlfParam->prevIdx = uiIndex;
  m_pcTempAlfParam->alf_flag = 1;
  m_pcTempAlfParam->chroma_idc = 0;
  m_pcBestAlfParam->temproalPredFlag = false;
  m_varImg = m_varImgMethods;

  xfilterFrame_en(pDec, pRest, m_pcTempAlfParam->filtNo, LumaStride );
  uiDist     = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Y), pcPicRest->getAddr(COMPONENT_Y), pcPicOrg->getWidth(COMPONENT_Y), pcPicOrg->getHeight(COMPONENT_Y), pcPicOrg->getStride(COMPONENT_Y));
  xCalcRDCost(m_pcTempAlfParam, uiRate, uiDist, dCost, pSlice);

  if (dCost < rdMinCost)
  {
    rdMinCost = dCost;
    ruiMinDist = uiDist;
    ruiMinRate = uiRate;
    m_pcBestAlfParam->temproalPredFlag = true;
    copyALFParam(m_pcBestAlfParam, m_pcTempAlfParam);
    m_pcBestAlfParam->prevIdx = uiIndex;
    m_pcBestAlfParam->alf_flag = 1;
    m_pcBestAlfParam->cu_control_flag = 0;

    if (m_pcBestAlfParam->tap >= m_ALF_MIN_NUM_TAP) // fixed a bug with HM-3
    {
      pcPicRest->copyToPic(m_pcPicYuvBest, COMPONENT_Y);
    }
    m_pcEntropyCoder->setAlfCtrl( false );
    m_pcEntropyCoder->setMaxAlfCtrlDepth( 0 );
  }

  m_pcEntropyCoder->setAlfCtrl( true );
  m_pcTempAlfParam->cu_control_flag = 1;
  UInt uiBestDepth = 0;
  Bool bChanged = false;

  for (UInt uiDepth = 0; uiDepth < pSlice->getSPS()->getMaxTotalCUDepth(); uiDepth++)
  {
    m_pcEntropyCoder->setMaxAlfCtrlDepth(uiDepth);
    pcPicRest->copyToPic(m_pcPicYuvTmp, COMPONENT_Y);

    UInt64 uiTmpRate, uiTmpDist;
    Double dTmpCost;
    //m_pcPicYuvTmp: filtered signal, pcPicDec: orig reconst
    xSetCUAlfCtrlFlags_qc(uiDepth, pcPicOrg, pcPicDec, m_pcPicYuvTmp, uiTmpDist, m_pcTempAlfParam); //set up varImg here

    xCalcRDCost(m_pcTempAlfParam, uiTmpRate, uiTmpDist, dTmpCost, pSlice);
    if (dTmpCost < rdMinCost)
    {
      bChanged = true;
      uiBestDepth = uiDepth;
      rdMinCost = dTmpCost;
      ruiMinDist = uiTmpDist;
      ruiMinRate = uiTmpRate;
      m_pcPicYuvTmp->copyToPic(m_pcPicYuvBest, COMPONENT_Y);
      m_pcBestAlfParam->temproalPredFlag = false;
      copyALFParam(m_pcBestAlfParam, m_pcTempAlfParam);
      m_pcBestAlfParam->temproalPredFlag = true;
      //save maskImg
      xCopyTmpAlfCtrlFlagsFrom();
    }
  }

  if (bChanged)
  {
    m_pcBestAlfParam->alf_flag = true;
    m_pcBestAlfParam->prevIdx = uiIndex;
    m_pcBestAlfParam->cu_control_flag = true;
    m_pcBestAlfParam->alf_max_depth = uiBestDepth;
    m_pcEntropyCoder->setMaxAlfCtrlDepth(uiBestDepth);
    //m_pcPicYuvBest->copyToPicLuma(pcPicRest);//copy m_pcPicYuvBest to pcPicRest
    xCopyTmpAlfCtrlFlagsTo();
    xCopyDecToRestCUs(pcPicDec, pcPicRest); //pcPicRest = pcPicDec
    pcPicRest->copyToPic( m_pcPicYuvBest, COMPONENT_Y );
  }
  else
  {
    m_pcEntropyCoder->setAlfCtrl( uiTmpAlfCtrlFlag );
    m_pcEntropyCoder->setMaxAlfCtrlDepth( uiTmpMaxDepth );
  }


  if( !m_pcBestAlfParam->temproalPredFlag )
  {
    return false;
  }

  m_pcBestAlfParam->chroma_idc = 0;

  if( pAlfParam->chroma_idc != 0 )
  {
#if JVET_C0038_GALF
    memcpy( pAlfParam->coeff_chroma, pAlfParam->alfCoeffChroma, sizeof(Int)*m_ALF_MAX_NUM_COEF_C );
    initVarForChroma(pAlfParam);
    xFrameChroma(pAlfParam, pcPicDec, pcPicRest, 0);
    xFrameChroma(pAlfParam, pcPicDec, pcPicRest, 1);
#else
    Int *qh = new Int[m_ALF_MAX_NUM_COEF_C];

    memcpy( qh, m_pcTempAlfParam->alfCoeffChroma, sizeof(Int)*m_ALF_MAX_NUM_COEF_C );

    //filter chroma
    xFrameChroma(pcPicDec, pcPicRest, qh, m_pcTempAlfParam->tap_chroma, 0);
    xFrameChroma(pcPicDec, pcPicRest, qh, m_pcTempAlfParam->tap_chroma, 1);
#endif
    UInt64 uiDistOrg;

    uiDist    = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cb), pcPicRest->getAddr(COMPONENT_Cb), pcPicOrg->getWidth(COMPONENT_Cb), pcPicOrg->getHeight(COMPONENT_Cb), pcPicOrg->getStride(COMPONENT_Cb));
    uiDistOrg = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cb), pcPicDec->getAddr(COMPONENT_Cb), pcPicOrg->getWidth(COMPONENT_Cb), pcPicOrg->getHeight(COMPONENT_Cb), pcPicOrg->getStride(COMPONENT_Cb));

    if( uiDist < uiDistOrg )
    {
      m_pcBestAlfParam->chroma_idc |= 2;
      pcPicRest->copyToPic(m_pcPicYuvBest,COMPONENT_Cb);
    }
    else
    {
      pcPicDec->copyToPic(m_pcPicYuvBest,COMPONENT_Cb);
    }

    uiDist    = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cr), pcPicRest->getAddr(COMPONENT_Cr), pcPicOrg->getWidth(COMPONENT_Cr), pcPicOrg->getHeight(COMPONENT_Cr), pcPicOrg->getStride(COMPONENT_Cr));
    uiDistOrg = xCalcSSD(pcPicOrg->getAddr(COMPONENT_Cr), pcPicDec->getAddr(COMPONENT_Cr), pcPicOrg->getWidth(COMPONENT_Cr), pcPicOrg->getHeight(COMPONENT_Cr), pcPicOrg->getStride(COMPONENT_Cr));

    if( uiDist < uiDistOrg )
    {
      m_pcBestAlfParam->chroma_idc |= 1;
      pcPicRest->copyToPic(m_pcPicYuvBest, COMPONENT_Cr);
    }
    else
    {
      pcPicDec->copyToPic(m_pcPicYuvBest, COMPONENT_Cr);
    }
#if !JVET_C0038_GALF 
    delete [] qh;
#endif
  }
  else
  {
    pcPicDec->copyToPic(m_pcPicYuvBest, COMPONENT_Cb);
    pcPicDec->copyToPic(m_pcPicYuvBest, COMPONENT_Cr);
  }

  return true;
}
#endif
#if JVET_C0038_GALF
Void TEncAdaptiveLoopFilter::xDeriveGlobalEyFromLgrTapFilter( Double **E0, Double *y0, Double **E1, Double *y1, Int *pattern0, Int *pattern1 )
{
  Int i, j, l, k = 0;
  for (i = 0; i < m_MAX_SQR_FILT_LENGTH; i ++)
  {
    if ( pattern0[i] > 0 )
    {
      if ( pattern1[i] > 0 )
      {
        l = 0;
        y1[pattern1[i]-1] = y0[k]; 
        for (j = 0; j < m_MAX_SQR_FILT_LENGTH; j ++)
        {
          if (pattern0[j] > 0)
          {
            if (pattern1[j] > 0)
            {
              E1[pattern1[i]-1][pattern1[j]-1] = E0[k][l];
            }
            l++;
          }
        }
      }
      k++;
    }
  }
}
Void TEncAdaptiveLoopFilter::xDeriveLocalEyFromLgrTapFilter ( Double *y0,  Double *y1, Int *pattern0, Int *pattern1 )
{
  Int i, k = 0;
  for (i = 0; i < m_MAX_SQR_FILT_LENGTH; i ++)
  {
    if (pattern0[i] > 0)
    {
      if (pattern1[i] > 0)
      {
        y1[pattern1[i]-1] = y0[k]; 
      }
      k++;
    }
  }
}

#endif

#endif
