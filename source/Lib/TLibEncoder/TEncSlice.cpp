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

/** \file     TEncSlice.cpp
    \brief    slice encoder class
*/

#include "TEncTop.h"
#include "TEncSlice.h"
#include <math.h>

//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TEncSlice::TEncSlice()
 : m_encCABACTableIdx(I_SLICE)
{
  m_apcPicYuvPred = NULL;
  m_apcPicYuvResi = NULL;

  m_pdRdPicLambda = NULL;
  m_pdRdPicQp     = NULL;
  m_piRdPicQp     = NULL;
}

TEncSlice::~TEncSlice()
{
}

Void TEncSlice::create( Int iWidth, Int iHeight, ChromaFormat chromaFormat, UInt iMaxCUWidth, UInt iMaxCUHeight, UChar uhTotalDepth )
{
  // create prediction picture
  if ( m_apcPicYuvPred == NULL )
  {
    m_apcPicYuvPred  = new TComPicYuv;
    m_apcPicYuvPred->create( iWidth, iHeight, chromaFormat, iMaxCUWidth, iMaxCUHeight, uhTotalDepth, true );
  }

  // create residual picture
  if( m_apcPicYuvResi == NULL )
  {
    m_apcPicYuvResi  = new TComPicYuv;
    m_apcPicYuvResi->create( iWidth, iHeight, chromaFormat, iMaxCUWidth, iMaxCUHeight, uhTotalDepth, true );
  }
}

Void TEncSlice::destroy()
{
  // destroy prediction picture
  if ( m_apcPicYuvPred )
  {
    m_apcPicYuvPred->destroy();
    delete m_apcPicYuvPred;
    m_apcPicYuvPred  = NULL;
  }

  // destroy residual picture
  if ( m_apcPicYuvResi )
  {
    m_apcPicYuvResi->destroy();
    delete m_apcPicYuvResi;
    m_apcPicYuvResi  = NULL;
  }

  // free lambda and QP arrays
  if ( m_pdRdPicLambda )
  {
    xFree( m_pdRdPicLambda );
    m_pdRdPicLambda = NULL;
  }
  if ( m_pdRdPicQp )
  {
    xFree( m_pdRdPicQp );
    m_pdRdPicQp = NULL;
  }
  if ( m_piRdPicQp )
  {
    xFree( m_piRdPicQp );
    m_piRdPicQp = NULL;
  }
}

Void TEncSlice::init( TEncTop* pcEncTop )
{
  m_pcCfg             = pcEncTop;
  m_pcListPic         = pcEncTop->getListPic();

  m_pcGOPEncoder      = pcEncTop->getGOPEncoder();
  m_pcCuEncoder       = pcEncTop->getCuEncoder();
  m_pcPredSearch      = pcEncTop->getPredSearch();

  m_pcEntropyCoder    = pcEncTop->getEntropyCoder();
  m_pcSbacCoder       = pcEncTop->getSbacCoder();
  m_pcBinCABAC        = pcEncTop->getBinCABAC();
  m_pcTrQuant         = pcEncTop->getTrQuant();

  m_pcRdCost          = pcEncTop->getRdCost();
#if JVET_C0024_QTBT
  m_ppppcRDSbacCoder  = pcEncTop->getRDSbacCoder();
#else
  m_pppcRDSbacCoder   = pcEncTop->getRDSbacCoder();
#endif
  m_pcRDGoOnSbacCoder = pcEncTop->getRDGoOnSbacCoder();

  // create lambda and QP arrays
  m_pdRdPicLambda     = (Double*)xMalloc( Double, m_pcCfg->getDeltaQpRD() * 2 + 1 );
  m_pdRdPicQp         = (Double*)xMalloc( Double, m_pcCfg->getDeltaQpRD() * 2 + 1 );
  m_piRdPicQp         = (Int*   )xMalloc( Int,    m_pcCfg->getDeltaQpRD() * 2 + 1 );
  m_pcRateCtrl        = pcEncTop->getRateCtrl();
}

#if WCG_LUMA_DQP_CM_SCALE//  copied from initEncSlice
Void TEncSlice::updateLambda(TComSlice* rpcSlice, Int dQP) {
    Int    NumberBFrames = ( m_pcCfg->getGOPSize() - 1 );
    Int    SHIFT_QP = 12;    
    Double dLambda_scale = 1.0 - Clip3( 0.0, 0.5, 0.05*(Double)(rpcSlice->getPic()->isField() ? NumberBFrames/2 : NumberBFrames) );
#if FULL_NBIT
    Int    bitdepth_luma_qp_scale = 6 * (rpcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8);
#else
    Int    bitdepth_luma_qp_scale = 0;
#endif
#if JVET_B0039_QP_FIX  // dQP may have been adjusted by lambda, get the QP without the lambda offset adjustment
  Int lambdaQPOffset = rpcSlice->getSliceQpLambdaOffset();
  dQP -= lambdaQPOffset;
#endif
    Double qp_temp = (Double) dQP + bitdepth_luma_qp_scale - SHIFT_QP;
#if FULL_NBIT
    Double qp_temp_orig = (Double) dQP - SHIFT_QP;
#endif
 enum SliceType eSliceType = rpcSlice->getSliceType();
    Double dQPFactor = m_pcCfg->getGOPEntry(m_GopID).m_QPFactor;
   if ( eSliceType==I_SLICE )
   {
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
     if (m_pcCfg->getIntraQpFactor() >= 0.0 && m_pcCfg->getGOPEntry(m_GopID).m_sliceType != I_SLICE)
     {
       dQPFactor = m_pcCfg->getIntraQpFactor();
     }
     else
     {
       if (m_pcCfg->getLambdaFromQPEnable())
       {
         dQPFactor = 0.57;
       }
       else
       {
         dQPFactor = 0.57*dLambda_scale;
       }
     }
#else
     dQPFactor=0.57*dLambda_scale;
#endif
   }
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
   else if (m_pcCfg->getLambdaFromQPEnable())
   {
     dQPFactor = 0.57*dQPFactor;
   }
#endif
    Double dLambda = dQPFactor*pow( 2.0, qp_temp/3.0 );
    Int depth = rpcSlice->getDepth();
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    if (!m_pcCfg->getLambdaFromQPEnable() && depth>0)
#else
    if ( depth>0 )
#endif
    {  
#if FULL_NBIT
      qp_temp_orig = rpcSlice->getSliceQp() - SHIFT_QP; // avoid lambda  over adjustment,  use slice_qp here
      dLambda *= Clip3( 2.00, 4.00, (qp_temp_orig / 6.0) ); // (j == B_SLICE && p_cur_frm->layer != 0 )
#else
      Int qp_temp_slice = rpcSlice->getSliceQp() + bitdepth_luma_qp_scale - SHIFT_QP; // avoid lambda  over adjustment,  use slice_qp here
      dLambda *= Clip3( 2.00, 4.00, (qp_temp_slice / 6.0) ); // (j == B_SLICE && p_cur_frm->layer != 0 )
#endif
    }
    if ( !m_pcCfg->getUseHADME() && rpcSlice->getSliceType( ) != I_SLICE )
    {
      dLambda *= 0.95;
    }

#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    const Int temporalId = m_pcCfg->getGOPEntry(m_GopID).m_temporalId;
    const std::vector<Double> &intraLambdaModifiers = m_pcCfg->getIntraLambdaModifier();
    Double lambdaModifier;
    if( rpcSlice->getSliceType( ) != I_SLICE || intraLambdaModifiers.empty())
    {
      lambdaModifier = m_pcCfg->getLambdaModifier(temporalId);
    }
    else
    {
      lambdaModifier = intraLambdaModifiers[(temporalId < intraLambdaModifiers.size()) ? temporalId : (intraLambdaModifiers.size() - 1)];
    }
    dLambda *= lambdaModifier;
#endif

#if JVET_B0039_QP_FIX
    Double lambdaRef = 0.57*pow(2.0, qp_temp/3.0);
    Double qpOffset = floor((3.0*log(dLambda/lambdaRef)/log(2.0)) +0.5);
    dQP += qpOffset;
#endif
    Int qpBDoffset = rpcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA);
    Int iQP = max( -qpBDoffset, min( MAX_QP, (Int) floor( dQP + 0.5 ) ) ); 
    setUpLambda(rpcSlice, dLambda, iQP);
}
#endif


Void
TEncSlice::setUpLambda(TComSlice* slice, const Double dLambda, Int iQP)
{
  // store lambda
  m_pcRdCost ->setLambda( dLambda, slice->getSPS()->getBitDepths() );

  // for RDO
  // in RdCost there is only one lambda because the luma and chroma bits are not separated, instead we weight the distortion of chroma.
  Double dLambdas[MAX_NUM_COMPONENT] = { dLambda };
  for(UInt compIdx=1; compIdx<MAX_NUM_COMPONENT; compIdx++)
  {
    const ComponentID compID=ComponentID(compIdx);
    Int chromaQPOffset = slice->getPPS()->getQpOffset(compID) + slice->getSliceChromaQpDelta(compID);
    Int qpc=(iQP + chromaQPOffset < 0) ? iQP : getScaledChromaQP(iQP + chromaQPOffset, m_pcCfg->getChromaFormatIdc());
    Double tmpWeight = pow( 2.0, (iQP-qpc)/3.0 );  // takes into account of the chroma qp mapping and chroma qp Offset
#if RRC_FRAC_CHROMA_QP_OFFSET
    Int qpc_2 = (iQP + chromaQPOffset < 0) ? iQP : getScaledChromaQP(iQP + chromaQPOffset + 1, m_pcCfg->getChromaFormatIdc());
    Double qpcD = (qpc*(1. - RRC_FRAC_CHROMA_QP_OFFSET_ADD) + qpc_2*RRC_FRAC_CHROMA_QP_OFFSET_ADD);
    tmpWeight = pow(2.0, (iQP - qpcD) / 3.0);  // takes into account of the chroma qp mapping and chroma qp Offset
#endif
    m_pcRdCost->setDistortionWeight(compID, tmpWeight);
    dLambdas[compIdx]=dLambda/tmpWeight;
  }

#if RDOQ_CHROMA_LAMBDA
// for RDOQ
  m_pcTrQuant->setLambdas( dLambdas );
#else
  m_pcTrQuant->setLambda( dLambda );
#endif

// For SAO
  slice   ->setLambdas( dLambdas );
}



/**
 - non-referenced frame marking
 - QP computation based on temporal structure
 - lambda computation based on QP
 - set temporal layer ID and the parameter sets
 .
 \param pcPic         picture class
 \param pocLast       POC of last picture
 \param pocCurr       current POC
 \param iNumPicRcvd   number of received pictures
 \param iGOPid        POC offset for hierarchical structure
 \param rpcSlice      slice header class
 \param isField       true for field coding
 */

Void TEncSlice::initEncSlice( TComPic* pcPic, Int pocLast, Int pocCurr, Int iGOPid, TComSlice*& rpcSlice, Bool isField )
{
  Double dQP;
  Double dLambda;

  rpcSlice = pcPic->getSlice(0);
  rpcSlice->setSliceBits(0);
  rpcSlice->setPic( pcPic );
  rpcSlice->initSlice();
  rpcSlice->setPicOutputFlag( true );
  rpcSlice->setPOC( pocCurr );
#if VCEG_AZ06_IC
  rpcSlice->setApplyIC( false );
#endif  
#if BT_FORCE
  for (Int i = 0; i < 2; i++)
    rpcSlice->setMaxFQTLevel(0, i);
#endif
#if WCG_LUMA_DQP_CM_SCALE
  m_GopID = iGOPid;
#endif  
  // depth computation based on GOP size
  Int depth;
  {
    Int poc = rpcSlice->getPOC();
    if(isField)
    {
      poc = (poc/2) % (m_pcCfg->getGOPSize()/2);
    }
    else
    {
      poc = poc % m_pcCfg->getGOPSize();   
    }

    if ( poc == 0 )
    {
      depth = 0;
    }
    else
    {
      Int step = m_pcCfg->getGOPSize();
      depth    = 0;
      for( Int i=step>>1; i>=1; i>>=1 )
      {
        for ( Int j=i; j<m_pcCfg->getGOPSize(); j+=step )
        {
          if ( j == poc )
          {
            i=0;
            break;
          }
        }
        step >>= 1;
        depth++;
      }
    }

    if(m_pcCfg->getHarmonizeGopFirstFieldCoupleEnabled() && poc != 0)
    {
      if (isField && ((rpcSlice->getPOC() % 2) == 1))
      {
        depth ++;
      }
    }
  }

  // slice type
  SliceType eSliceType;

  eSliceType=B_SLICE;
  if(!(isField && pocLast == 1) || !m_pcCfg->getEfficientFieldIRAPEnabled())
  {
    if(m_pcCfg->getDecodingRefreshType() == 3)
    {
      eSliceType = (pocLast == 0 || pocCurr % m_pcCfg->getIntraPeriod() == 0             || m_pcGOPEncoder->getGOPSize() == 0) ? I_SLICE : eSliceType;
    }
    else
    {
      eSliceType = (pocLast == 0 || (pocCurr - (isField ? 1 : 0)) % m_pcCfg->getIntraPeriod() == 0 || m_pcGOPEncoder->getGOPSize() == 0) ? I_SLICE : eSliceType;
    }
  }

  rpcSlice->setSliceType    ( eSliceType );

  // ------------------------------------------------------------------------------------------------------------------
  // Non-referenced frame marking
  // ------------------------------------------------------------------------------------------------------------------

  if(pocLast == 0)
  {
    rpcSlice->setTemporalLayerNonReferenceFlag(false);
  }
  else
  {
    rpcSlice->setTemporalLayerNonReferenceFlag(!m_pcCfg->getGOPEntry(iGOPid).m_refPic);
  }
  rpcSlice->setReferenced(true);

  // ------------------------------------------------------------------------------------------------------------------
  // QP setting
  // ------------------------------------------------------------------------------------------------------------------

#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
  dQP = m_pcCfg->getQPForPicture(iGOPid, rpcSlice);
#else
  dQP = m_pcCfg->getQP();
  if(eSliceType!=I_SLICE)
  {
#if WCG_LUMA_DQP_CM_SCALE
 if (!(( m_pcCfg->getMaxDeltaQP() == 0) && (m_pcCfg->getUseLumaDeltaQp() == 0) && (dQP == -rpcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA) ) && (rpcSlice->getPPS()->getTransquantBypassEnableFlag())))
#else
    if (!(( m_pcCfg->getMaxDeltaQP() == 0 ) && (dQP == -rpcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA) ) && (rpcSlice->getPPS()->getTransquantBypassEnableFlag())))
#endif
    {
      dQP += m_pcCfg->getGOPEntry(iGOPid).m_QPOffset;
    }
  }

  // modify QP
  Int* pdQPs = m_pcCfg->getdQPs();
  if ( pdQPs )
  {
    dQP += pdQPs[ rpcSlice->getPOC() ];
  }

  if (m_pcCfg->getCostMode()==COST_LOSSLESS_CODING)
  {
    dQP=LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP;
    m_pcCfg->setDeltaQpRD(0);
  }
#endif

  // ------------------------------------------------------------------------------------------------------------------
  // Lambda computation
  // ------------------------------------------------------------------------------------------------------------------

#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
  const Int temporalId=m_pcCfg->getGOPEntry(iGOPid).m_temporalId;
  const std::vector<Double> &intraLambdaModifiers=m_pcCfg->getIntraLambdaModifier();
#endif
  Int iQP;
  Double dOrigQP = dQP;

  // pre-compute lambda and QP values for all possible QP candidates
  for ( Int iDQpIdx = 0; iDQpIdx < 2 * m_pcCfg->getDeltaQpRD() + 1; iDQpIdx++ )
  {
    // compute QP value
    dQP = dOrigQP + ((iDQpIdx+1)>>1)*(iDQpIdx%2 ? -1 : 1);

    // compute lambda value
    Int    NumberBFrames = ( m_pcCfg->getGOPSize() - 1 );
    Int    SHIFT_QP = 12;

    Double dLambda_scale = 1.0 - Clip3( 0.0, 0.5, 0.05*(Double)(isField ? NumberBFrames/2 : NumberBFrames) );

#if FULL_NBIT
    Int    bitdepth_luma_qp_scale = 6 * (rpcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8);
#else
    Int    bitdepth_luma_qp_scale = 0;
#endif
    Double qp_temp = (Double) dQP + bitdepth_luma_qp_scale - SHIFT_QP;
#if FULL_NBIT
    Double qp_temp_orig = (Double) dQP - SHIFT_QP;
#endif
    // Case #1: I or P-slices (key-frame)
    Double dQPFactor = m_pcCfg->getGOPEntry(iGOPid).m_QPFactor;
    if ( eSliceType==I_SLICE )
    {
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
      if (m_pcCfg->getIntraQpFactor()>=0.0 && m_pcCfg->getGOPEntry(iGOPid).m_sliceType != I_SLICE)
      {
        dQPFactor=m_pcCfg->getIntraQpFactor();
      }
      else
      {
        if(m_pcCfg->getLambdaFromQPEnable())
        {
          dQPFactor=0.57;
        }
        else
        {
          dQPFactor=0.57*dLambda_scale;
        }
      }
#else
      dQPFactor=0.57*dLambda_scale;
#endif
    }
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    else if( m_pcCfg->getLambdaFromQPEnable() )
    {
      dQPFactor=0.57*dQPFactor;
    }
#endif
    dLambda = dQPFactor*pow( 2.0, qp_temp/3.0 );

#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    if(!m_pcCfg->getLambdaFromQPEnable() && depth>0)
#else
    if ( depth>0 )
#endif
    {
#if FULL_NBIT
        dLambda *= Clip3( 2.00, 4.00, (qp_temp_orig / 6.0) ); // (j == B_SLICE && p_cur_frm->layer != 0 )
#else
        dLambda *= Clip3( 2.00, 4.00, (qp_temp / 6.0) ); // (j == B_SLICE && p_cur_frm->layer != 0 )
#endif
    }

    // if hadamard is used in ME process
    if ( !m_pcCfg->getUseHADME() && rpcSlice->getSliceType( ) != I_SLICE )
    {
      dLambda *= 0.95;
    }
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
    Double lambdaModifier;
    if( rpcSlice->getSliceType( ) != I_SLICE || intraLambdaModifiers.empty())
    {
      lambdaModifier = m_pcCfg->getLambdaModifier( temporalId );
    }
    else
    {
      lambdaModifier = intraLambdaModifiers[ (temporalId < intraLambdaModifiers.size()) ? temporalId : (intraLambdaModifiers.size()-1) ];
    }
    dLambda *= lambdaModifier;
#endif

#if JVET_B0039_QP_FIX
    Double lambdaRef = 0.57*pow(2.0, qp_temp/3.0);
    // qp correction to get HM lambda
    Double qpOffset = floor((3.0*log(dLambda/lambdaRef)/log(2.0)) +0.5);
    dQP += qpOffset;
#if WCG_LUMA_DQP_CM_SCALE
    if (iDQpIdx==0)
      rpcSlice->setSliceQpLambdaOffset( qpOffset );        // save the qpOffset for updateLambda to get original slice qp
#endif
#endif

    iQP = max( -rpcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), min( MAX_QP, (Int) floor( dQP + 0.5 ) ) );

    m_pdRdPicLambda[iDQpIdx] = dLambda;
    m_pdRdPicQp    [iDQpIdx] = dQP;
    m_piRdPicQp    [iDQpIdx] = iQP;
  }

  // obtain dQP = 0 case
  dLambda = m_pdRdPicLambda[0];
  dQP     = m_pdRdPicQp    [0];
  iQP     = m_piRdPicQp    [0];

#if !JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
  const Int temporalId=m_pcCfg->getGOPEntry(iGOPid).m_temporalId;
  const std::vector<Double> &intraLambdaModifiers=m_pcCfg->getIntraLambdaModifier();
#endif

#if !JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
  Double lambdaModifier;
  if( rpcSlice->getSliceType( ) != I_SLICE || intraLambdaModifiers.empty())
  {
    lambdaModifier = m_pcCfg->getLambdaModifier( temporalId );
  }
  else
  {
    lambdaModifier = intraLambdaModifiers[ (temporalId < intraLambdaModifiers.size()) ? temporalId : (intraLambdaModifiers.size()-1) ];
  }

  dLambda *= lambdaModifier;
#endif

  setUpLambda(rpcSlice, dLambda, iQP);
#if WCG_LUMA_DQP_CM_SCALE
  // cost = Distortion + Lambda*R, 
  // when QP is adjusted by luma, distortion is changed, so we have to adjust lambda to match the distortion, then the cost function becomes
  // costA = Distortion + AdjustedLambda * R          -- currently, costA is still used when calculating intermediate cost of using SAD, HAD, resisual etc. 
  // an alternative way is to weight the distortion to before the luma QP adjustment, then the cost function becomes
  // costB = weightedDistortion + Lambda * R          -- currently, costB is used to calculat final cost, and when DF_FUNC is DF_DEFAULT 
  m_pcRdCost->saveUnadjustedLambda();
#endif

  if (m_pcCfg->getFastMEForGenBLowDelayEnabled())
  {
    // restore original slice type

    if(!(isField && pocLast == 1) || !m_pcCfg->getEfficientFieldIRAPEnabled())
    {
      if(m_pcCfg->getDecodingRefreshType() == 3)
      {
        eSliceType = (pocLast == 0 || (pocCurr)                     % m_pcCfg->getIntraPeriod() == 0 || m_pcGOPEncoder->getGOPSize() == 0) ? I_SLICE : eSliceType;
      }
      else
      {
        eSliceType = (pocLast == 0 || (pocCurr - (isField ? 1 : 0)) % m_pcCfg->getIntraPeriod() == 0 || m_pcGOPEncoder->getGOPSize() == 0) ? I_SLICE : eSliceType;
      }
    }

    rpcSlice->setSliceType        ( eSliceType );
  }

  if (m_pcCfg->getUseRecalculateQPAccordingToLambda())
  {
    dQP = xGetQPValueAccordingToLambda( dLambda );
    iQP = max( -rpcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), min( MAX_QP, (Int) floor( dQP + 0.5 ) ) );
  }

  rpcSlice->setSliceQp           ( iQP );
#if !DISABLE_DQP_SIGNALING
	rpcSlice->setSliceBgQp         ( iQP );
#endif
#if ADAPTIVE_QP_SELECTION
  rpcSlice->setSliceQpBase       ( iQP );
#endif
  rpcSlice->setSliceQpDelta      ( 0 );
  rpcSlice->setSliceChromaQpDelta( COMPONENT_Cb, 0 );
  rpcSlice->setSliceChromaQpDelta( COMPONENT_Cr, 0 );
  rpcSlice->setUseChromaQpAdj( rpcSlice->getPPS()->getPpsRangeExtension().getChromaQpOffsetListEnabledFlag() );
  rpcSlice->setNumRefIdx(REF_PIC_LIST_0,m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive);
  rpcSlice->setNumRefIdx(REF_PIC_LIST_1,m_pcCfg->getGOPEntry(iGOPid).m_numRefPicsActive);

  if ( m_pcCfg->getDeblockingFilterMetric() )
  {
    rpcSlice->setDeblockingFilterOverrideFlag(true);
    rpcSlice->setDeblockingFilterDisable(false);
    rpcSlice->setDeblockingFilterBetaOffsetDiv2( 0 );
    rpcSlice->setDeblockingFilterTcOffsetDiv2( 0 );
  }
  else if (rpcSlice->getPPS()->getDeblockingFilterControlPresentFlag())
  {
    rpcSlice->setDeblockingFilterOverrideFlag( rpcSlice->getPPS()->getDeblockingFilterOverrideEnabledFlag() );
    rpcSlice->setDeblockingFilterDisable( rpcSlice->getPPS()->getPicDisableDeblockingFilterFlag() );
    if ( !rpcSlice->getDeblockingFilterDisable())
    {
      if ( rpcSlice->getDeblockingFilterOverrideFlag() && eSliceType!=I_SLICE)
      {
        rpcSlice->setDeblockingFilterBetaOffsetDiv2( m_pcCfg->getGOPEntry(iGOPid).m_betaOffsetDiv2 + m_pcCfg->getLoopFilterBetaOffset()  );
        rpcSlice->setDeblockingFilterTcOffsetDiv2( m_pcCfg->getGOPEntry(iGOPid).m_tcOffsetDiv2 + m_pcCfg->getLoopFilterTcOffset() );
      }
      else
      {
        rpcSlice->setDeblockingFilterBetaOffsetDiv2( m_pcCfg->getLoopFilterBetaOffset() );
        rpcSlice->setDeblockingFilterTcOffsetDiv2( m_pcCfg->getLoopFilterTcOffset() );
      }
    }
  }
  else
  {
    rpcSlice->setDeblockingFilterOverrideFlag( false );
    rpcSlice->setDeblockingFilterDisable( false );
    rpcSlice->setDeblockingFilterBetaOffsetDiv2( 0 );
    rpcSlice->setDeblockingFilterTcOffsetDiv2( 0 );
  }

  rpcSlice->setDepth            ( depth );

  pcPic->setTLayer( m_pcCfg->getGOPEntry(iGOPid).m_temporalId );
  if(eSliceType==I_SLICE)
  {
    pcPic->setTLayer(0);
  }
  rpcSlice->setTLayer( pcPic->getTLayer() );

  assert( m_apcPicYuvPred );
  assert( m_apcPicYuvResi );

  pcPic->setPicYuvPred( m_apcPicYuvPred );
  pcPic->setPicYuvResi( m_apcPicYuvResi );
  rpcSlice->setSliceMode            ( m_pcCfg->getSliceMode()            );
  rpcSlice->setSliceArgument        ( m_pcCfg->getSliceArgument()        );
  rpcSlice->setSliceSegmentMode     ( m_pcCfg->getSliceSegmentMode()     );
  rpcSlice->setSliceSegmentArgument ( m_pcCfg->getSliceSegmentArgument() );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  rpcSlice->setMaxNumMergeCand        ( m_pcCfg->getMaxNumMergeCand() + (m_pcCfg->getAtmvp()? 2:0)        );
#else
  rpcSlice->setMaxNumMergeCand        ( m_pcCfg->getMaxNumMergeCand()        );
#endif
#if JVET_C0024_QTBT
  rpcSlice->setMaxBTSize(rpcSlice->isIntra() ? MAX_BT_SIZE: MAX_BT_SIZE_INTER);
#endif

#if PVC_AQS
  m_pesmErrorSensitivityMap = m_pcGOPEncoder->getEncTop()->getESMap();
  initErrorSensitivityMap(pcPic);
#endif
}

Void TEncSlice::resetQP( TComPic* pic, Int sliceQP, Double lambda )
{
  TComSlice* slice = pic->getSlice(0);

  // store lambda
  slice->setSliceQp( sliceQP );
#if ADAPTIVE_QP_SELECTION
  slice->setSliceQpBase ( sliceQP );
#endif
  setUpLambda(slice, lambda, sliceQP);
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

Void TEncSlice::setSearchRange( TComSlice* pcSlice )
{
  Int iCurrPOC = pcSlice->getPOC();
  Int iRefPOC;
  Int iGOPSize = m_pcCfg->getGOPSize();
  Int iOffset = (iGOPSize >> 1);
  Int iMaxSR = m_pcCfg->getSearchRange();
  Int iNumPredDir = pcSlice->isInterP() ? 1 : 2;

  for (Int iDir = 0; iDir < iNumPredDir; iDir++)
  {
    //RefPicList e = (RefPicList)iDir;
    RefPicList  e = ( iDir ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );
    for (Int iRefIdx = 0; iRefIdx < pcSlice->getNumRefIdx(e); iRefIdx++)
    {
      iRefPOC = pcSlice->getRefPic(e, iRefIdx)->getPOC();
      Int newSearchRange = Clip3(m_pcCfg->getMinSearchWindow(), iMaxSR, (iMaxSR*ADAPT_SR_SCALE*abs(iCurrPOC - iRefPOC) + iOffset) / iGOPSize);
      m_pcPredSearch->setAdaptiveSearchRange(iDir, iRefIdx, newSearchRange);
    }
  }
}

/**
 Multi-loop slice encoding for different slice QP

 \param pcPic    picture class
 */
Void TEncSlice::precompressSlice( TComPic* pcPic )
{
  // if deltaQP RD is not used, simply return
  if ( m_pcCfg->getDeltaQpRD() == 0 )
  {
    return;
  }

  if ( m_pcCfg->getUseRateCtrl() )
  {
    printf( "\nMultiple QP optimization is not allowed when rate control is enabled." );
    assert(0);
    return;
  }

  TComSlice* pcSlice        = pcPic->getSlice(getSliceIdx());

  if (pcSlice->getDependentSliceSegmentFlag())
  {
    // if this is a dependent slice segment, then it was optimised
    // when analysing the entire slice.
    return;
  }

  if (pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES)
  {
    // TODO: investigate use of average cost per CTU so that this Slice Mode can be used.
    printf( "\nUnable to optimise Slice-level QP if Slice Mode is set to FIXED_NUMBER_OF_BYTES\n" );
    assert(0);
    return;
  }

  Double     dPicRdCostBest = MAX_DOUBLE;
  UInt       uiQpIdxBest = 0;

  Double dFrameLambda;
#if FULL_NBIT
  Int    SHIFT_QP = 12 + 6 * (pcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8);
#else
  Int    SHIFT_QP = 12;
#endif

  // set frame lambda
  if (m_pcCfg->getGOPSize() > 1)
  {
    dFrameLambda = 0.68 * pow (2, (m_piRdPicQp[0]  - SHIFT_QP) / 3.0) * (pcSlice->isInterB()? 2 : 1);
  }
  else
  {
    dFrameLambda = 0.68 * pow (2, (m_piRdPicQp[0] - SHIFT_QP) / 3.0);
  }
  m_pcRdCost      ->setFrameLambda(dFrameLambda);

  // for each QP candidate
  for ( UInt uiQpIdx = 0; uiQpIdx < 2 * m_pcCfg->getDeltaQpRD() + 1; uiQpIdx++ )
  {
    pcSlice       ->setSliceQp             ( m_piRdPicQp    [uiQpIdx] );
#if ADAPTIVE_QP_SELECTION
    pcSlice       ->setSliceQpBase         ( m_piRdPicQp    [uiQpIdx] );
#endif
    setUpLambda(pcSlice, m_pdRdPicLambda[uiQpIdx], m_piRdPicQp    [uiQpIdx]);

    // try compress
    compressSlice   ( pcPic, true, m_pcCfg->getFastDeltaQp());

    UInt64 uiPicDist        = m_uiPicDist; // Distortion, as calculated by compressSlice.
    // NOTE: This distortion is the chroma-weighted SSE distortion for the slice.
    //       Previously a standard SSE distortion was calculated (for the entire frame).
    //       Which is correct?

    // TODO: Update loop filter, SAO and distortion calculation to work on one slice only.
    // m_pcGOPEncoder->preLoopFilterPicAll( pcPic, uiPicDist );

    // compute RD cost and choose the best
    Double dPicRdCost = m_pcRdCost->calcRdCost64( m_uiPicTotalBits, uiPicDist, true, DF_SSE_FRAME); // NOTE: Is the 'true' parameter really necessary?

    if ( dPicRdCost < dPicRdCostBest )
    {
      uiQpIdxBest    = uiQpIdx;
      dPicRdCostBest = dPicRdCost;
    }
  }

  // set best values
  pcSlice       ->setSliceQp             ( m_piRdPicQp    [uiQpIdxBest] );
#if ADAPTIVE_QP_SELECTION
  pcSlice       ->setSliceQpBase         ( m_piRdPicQp    [uiQpIdxBest] );
#endif
  setUpLambda(pcSlice, m_pdRdPicLambda[uiQpIdxBest], m_piRdPicQp    [uiQpIdxBest]);
}

Void TEncSlice::calCostSliceI(TComPic* pcPic) // TODO: this only analyses the first slice segment. What about the others?
{
  Double            iSumHadSlice      = 0;
  TComSlice * const pcSlice           = pcPic->getSlice(getSliceIdx());
  const TComSPS    &sps               = *(pcSlice->getSPS());
  const Int         shift             = sps.getBitDepth(CHANNEL_TYPE_LUMA)-8;
  const Int         offset            = (shift>0)?(1<<(shift-1)):0;

  pcSlice->setSliceSegmentBits(0);

  UInt startCtuTsAddr, boundingCtuTsAddr;
  xDetermineStartAndBoundingCtuTsAddr ( startCtuTsAddr, boundingCtuTsAddr, pcPic );

  for( UInt ctuTsAddr = startCtuTsAddr, ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap( startCtuTsAddr);
       ctuTsAddr < boundingCtuTsAddr;
       ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(++ctuTsAddr) )
  {
    // initialize CU encoder
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );
    pCtu->initCtu( pcPic, ctuRsAddr );

#if JVET_C0024_QTBT
    Int height  = min( sps.getCTUSize(),sps.getPicHeightInLumaSamples() - ctuRsAddr / pcPic->getFrameWidthInCtus() * sps.getCTUSize() );
    Int width   = min( sps.getCTUSize(), sps.getPicWidthInLumaSamples()  - ctuRsAddr % pcPic->getFrameWidthInCtus() * sps.getCTUSize() );
#else
    Int height  = min( sps.getMaxCUHeight(),sps.getPicHeightInLumaSamples() - ctuRsAddr / pcPic->getFrameWidthInCtus() * sps.getMaxCUHeight() );
    Int width   = min( sps.getMaxCUWidth(), sps.getPicWidthInLumaSamples()  - ctuRsAddr % pcPic->getFrameWidthInCtus() * sps.getMaxCUWidth() );
#endif

    Int iSumHad = m_pcCuEncoder->updateCtuDataISlice(pCtu, width, height);

    (m_pcRateCtrl->getRCPic()->getLCU(ctuRsAddr)).m_costIntra=(iSumHad+offset)>>shift;
    iSumHadSlice += (m_pcRateCtrl->getRCPic()->getLCU(ctuRsAddr)).m_costIntra;

  }
  m_pcRateCtrl->getRCPic()->setTotalIntraCost(iSumHadSlice);
}

/** \param pcPic   picture class
 */
Void TEncSlice::compressSlice( TComPic* pcPic, const Bool bCompressEntireSlice, const Bool bFastDeltaQP )
{
  // if bCompressEntireSlice is true, then the entire slice (not slice segment) is compressed,
  //   effectively disabling the slice-segment-mode.

  UInt   startCtuTsAddr;
  UInt   boundingCtuTsAddr;
  TComSlice* const pcSlice            = pcPic->getSlice(getSliceIdx());
  pcSlice->setSliceSegmentBits(0);
  xDetermineStartAndBoundingCtuTsAddr ( startCtuTsAddr, boundingCtuTsAddr, pcPic );
#if VCEG_AZ08_KLT_COMMON
  pcPic->getPicYuvRec()->fillPicRecBoundary(pcSlice->getSPS()->getBitDepths());
#endif
  if (bCompressEntireSlice)
  {
    boundingCtuTsAddr = pcSlice->getSliceCurEndCtuTsAddr();
    pcSlice->setSliceSegmentCurEndCtuTsAddr(boundingCtuTsAddr);
  }

  // initialize cost values - these are used by precompressSlice (they should be parameters).
  m_uiPicTotalBits  = 0;
  m_dPicRdCost      = 0; // NOTE: This is a write-only variable!
  m_uiPicDist       = 0;
#if BARC_REF_UP_SAMPLE_CNN
#if BARC_ONLY_FOR_INTRA
  if (pcSlice->getSliceType() == I_SLICE)
  {

#endif
	  pcPic->initCNNModel();
#if BARC_ONLY_FOR_INTRA
  }
#endif
#if COMPACT_RESOLUTION_CNN
#if BARC_ONLY_FOR_INTRA
  if (pcSlice->getSliceType() == I_SLICE)
  {

#endif
#if PRINT_INFO
	  Double endTime;
	  clock_t startClock = clock();
	  printf("POC %d downScaleOrgCRCNN start.\n", pcPic->getPOC());
#endif
	  pcPic->downScaleOrgCRCNN();
#if PRINT_INFO
	  endTime = (Double)(clock() - startClock) / CLOCKS_PER_SEC;
	  printf("POC %d enc downScaleOrgCRCNN done! Time: %.3f sec.\n", pcPic->getPOC(), endTime);
#endif
#if BARC_ONLY_FOR_INTRA
  }
#endif
#endif
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
	Int iOrgSliceQP;
	if (pcSlice->getSPS()->getBackgroundRefsPresent() && pcSlice->getSliceType() != I_SLICE)
	{
		iOrgSliceQP = pcSlice->getSliceQp();
	}
#endif

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  m_pcEntropyCoder->setStatsHandle ( pcSlice->getStatsHandle() );
#endif  
#if JVET_C0024_QTBT
  UInt uiMaxWIdx = g_aucConvertToBit[pcSlice->getSPS()->getCTUSize()];
  UInt uiMaxHIdx = g_aucConvertToBit[pcSlice->getSPS()->getCTUSize()];
  m_pcEntropyCoder->setEntropyCoder   ( m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST] );
#else
  m_pcEntropyCoder->setEntropyCoder   ( m_pppcRDSbacCoder[0][CI_CURR_BEST] );
#endif
  m_pcEntropyCoder->resetEntropy      ( pcSlice );

#if JVET_C0024_QTBT
  TEncBinCABAC* pRDSbacCoder = (TEncBinCABAC *) m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->getEncBinIf();
#else
  TEncBinCABAC* pRDSbacCoder = (TEncBinCABAC *) m_pppcRDSbacCoder[0][CI_CURR_BEST]->getEncBinIf();
#endif
  pRDSbacCoder->setBinCountingEnableFlag( false );
  pRDSbacCoder->setBinsCoded( 0 );

  TComBitCounter  tempBitCounter;
  const UInt      frameWidthInCtus = pcPic->getPicSym()->getFrameWidthInCtus();
  
  m_pcCuEncoder->setFastDeltaQp(bFastDeltaQP);

#if ALF_HM3_REFACTOR
  // initialize ALF parameters
  m_pcEntropyCoder->setAlfCtrl(false);
  m_pcEntropyCoder->setMaxAlfCtrlDepth(0); //unnecessary
#endif

  //------------------------------------------------------------------------------
  //  Weighted Prediction parameters estimation.
  //------------------------------------------------------------------------------
  // calculate AC/DC values for current picture
  if( pcSlice->getPPS()->getUseWP() || pcSlice->getPPS()->getWPBiPred() )
  {
    xCalcACDCParamSlice(pcSlice);
  }

  const Bool bWp_explicit = (pcSlice->getSliceType()==P_SLICE && pcSlice->getPPS()->getUseWP()) || (pcSlice->getSliceType()==B_SLICE && pcSlice->getPPS()->getWPBiPred());

  if ( bWp_explicit )
  {
    //------------------------------------------------------------------------------
    //  Weighted Prediction implemented at Slice level. SliceMode=2 is not supported yet.
    //------------------------------------------------------------------------------
    if ( pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES || pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES )
    {
      printf("Weighted Prediction is not supported with slice mode determined by max number of bins.\n"); exit(0);
    }

    xEstimateWPParamSlice( pcSlice );
    pcSlice->initWpScaling(pcSlice->getSPS());

    // check WP on/off
    xCheckWPEnable( pcSlice );
  }

#if ADAPTIVE_QP_SELECTION
  if( m_pcCfg->getUseAdaptQpSelect() && !(pcSlice->getDependentSliceSegmentFlag()))
  {
    // TODO: this won't work with dependent slices: they do not have their own QP. Check fix to mask clause execution with && !(pcSlice->getDependentSliceSegmentFlag())
    m_pcTrQuant->clearSliceARLCnt(); // TODO: this looks wrong for multiple slices - the results of all but the last slice will be cleared before they are used (all slices compressed, and then all slices encoded)
    if(pcSlice->getSliceType()!=I_SLICE)
    {
      Int qpBase = pcSlice->getSliceQpBase();
      pcSlice->setSliceQp(qpBase + m_pcTrQuant->getQpDelta(qpBase));
    }
  }
#endif

#if VCEG_AZ06_IC
  if ( m_pcCfg->getUseIC() )
  {
#if VCEG_AZ06_IC_SPEEDUP || JVET_C0024_QTBT
    pcSlice->xSetApplyIC();
#else
    pcSlice->setApplyIC( pcSlice->isIntra() ? false : true );
#endif
  }
#endif
#if BT_FORCE
#if !BT_FORCE_INTRA_ON
  if (!pcSlice->isIntra())
#endif
  {
    pcSlice->xSetMaxFQTLevel();
  }
#endif

  // Adjust initial state if this is the start of a dependent slice.
  {
    const UInt      ctuRsAddr               = pcPic->getPicSym()->getCtuTsToRsAddrMap( startCtuTsAddr);
    const UInt      currentTileIdx          = pcPic->getPicSym()->getTileIdxMap(ctuRsAddr);
    const TComTile *pCurrentTile            = pcPic->getPicSym()->getTComTile(currentTileIdx);
    const UInt      firstCtuRsAddrOfTile    = pCurrentTile->getFirstCtuRsAddr();
    if( pcSlice->getDependentSliceSegmentFlag() && ctuRsAddr != firstCtuRsAddrOfTile )
    {
      // This will only occur if dependent slice-segments (m_entropyCodingSyncContextState=true) are being used.
      if( pCurrentTile->getTileWidthInCtus() >= 2 || !m_pcCfg->getWaveFrontsynchro() )
      {
#if JVET_C0024_QTBT
        m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->loadContexts( &m_lastSliceSegmentEndContextState );
#else
        m_pppcRDSbacCoder[0][CI_CURR_BEST]->loadContexts( &m_lastSliceSegmentEndContextState );
#endif
      }
    }
  }

  // for every CTU in the slice segment (may terminate sooner if there is a byte limit on the slice-segment)
#if DIMD_INTRA_PRED
  for (UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < boundingCtuTsAddr; ++ctuTsAddr)
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    TComDataCU* pCtu = pcPic->getCtu(ctuRsAddr);
    pCtu->initCtu(pcPic, ctuRsAddr);
  }
#endif

  for( UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < boundingCtuTsAddr; ++ctuTsAddr )
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    // initialize CTU encoder
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );

#if !DIMD_INTRA_PRED
    pCtu->initCtu(pcPic, ctuRsAddr);
#endif

    // update CABAC state
    const UInt firstCtuRsAddrOfTile = pcPic->getPicSym()->getTComTile(pcPic->getPicSym()->getTileIdxMap(ctuRsAddr))->getFirstCtuRsAddr();
    const UInt tileXPosInCtus = firstCtuRsAddrOfTile % frameWidthInCtus;
    const UInt ctuXPosInCtus  = ctuRsAddr % frameWidthInCtus;
    
    if (ctuRsAddr == firstCtuRsAddrOfTile)
    {
#if JVET_C0024_QTBT
      m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->resetEntropy(pcSlice);
#else
      m_pppcRDSbacCoder[0][CI_CURR_BEST]->resetEntropy(pcSlice);
#endif
    }
    else if ( ctuXPosInCtus == tileXPosInCtus && m_pcCfg->getWaveFrontsynchro())
    {
      // reset and then update contexts to the state at the end of the top-right CTU (if within current slice and tile).
#if JVET_C0024_QTBT
      m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->resetEntropy(pcSlice);
#else
      m_pppcRDSbacCoder[0][CI_CURR_BEST]->resetEntropy(pcSlice);
#endif
      // Sync if the Top-Right is available.
      TComDataCU *pCtuUp = pCtu->getCtuAbove();
      if ( pCtuUp && ((ctuRsAddr%frameWidthInCtus+1) < frameWidthInCtus)  )
      {
        TComDataCU *pCtuTR = pcPic->getCtu( ctuRsAddr - frameWidthInCtus + 1 );
        if ( pCtu->CUIsFromSameSliceAndTile(pCtuTR) )
        {
          // Top-Right is available, we use it.
#if JVET_C0024_QTBT
          m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->loadContexts( &m_entropyCodingSyncContextState );
#else
          m_pppcRDSbacCoder[0][CI_CURR_BEST]->loadContexts( &m_entropyCodingSyncContextState );
#endif
        }
      }
    }

    // set go-on entropy coder (used for all trial encodings - the cu encoder and encoder search also have a copy of the same pointer)
    m_pcEntropyCoder->setEntropyCoder ( m_pcRDGoOnSbacCoder );
    m_pcEntropyCoder->setBitstream( &tempBitCounter );
    tempBitCounter.resetBits();

#if VCEG_AZ07_INIT_PREVFRAME
    if( pcSlice->getSliceType() != I_SLICE && ctuTsAddr == 0 )
    {
#if JVET_C0024_QTBT
      m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdx(), true, pcSlice->getCtxMapQPIdxforStore(), (pcSlice->getPOC() >  pcSlice->getStatsHandle()->m_uiLastIPOC) ); 
#else
      m_pppcRDSbacCoder[0][CI_CURR_BEST]->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdx(), true, pcSlice->getCtxMapQPIdxforStore(), (pcSlice->getPOC() >  pcSlice->getStatsHandle()->m_uiLastIPOC) ); 
#endif
    }
#endif

#if JVET_C0024_QTBT
    m_pcRDGoOnSbacCoder->load( m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST] ); // this copy is not strictly necessary here, but indicates that the GoOnSbacCoder
#else
    m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[0][CI_CURR_BEST] ); // this copy is not strictly necessary here, but indicates that the GoOnSbacCoder
#endif
                                                                     // is reset to a known state before every decision process.

    ((TEncBinCABAC*)m_pcRDGoOnSbacCoder->getEncBinIf())->setBinCountingEnableFlag(true);

    Double oldLambda = m_pcRdCost->getLambda();
    if ( m_pcCfg->getUseRateCtrl() )
    {
      Int estQP        = pcSlice->getSliceQp();
      Double estLambda = -1.0;
      Double bpp       = -1.0;

      if ( ( pcPic->getSlice( 0 )->getSliceType() == I_SLICE && m_pcCfg->getForceIntraQP() ) || !m_pcCfg->getLCULevelRC() )
      {
        estQP = pcSlice->getSliceQp();
      }
      else
      {
        bpp = m_pcRateCtrl->getRCPic()->getLCUTargetBpp(pcSlice->getSliceType());
        if ( pcPic->getSlice( 0 )->getSliceType() == I_SLICE)
        {
          estLambda = m_pcRateCtrl->getRCPic()->getLCUEstLambdaAndQP(bpp, pcSlice->getSliceQp(), &estQP);
        }
        else
        {
          estLambda = m_pcRateCtrl->getRCPic()->getLCUEstLambda( bpp );
          estQP     = m_pcRateCtrl->getRCPic()->getLCUEstQP    ( estLambda, pcSlice->getSliceQp() );
        }

        estQP     = Clip3( -pcSlice->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, estQP );

        m_pcRdCost->setLambda(estLambda, pcSlice->getSPS()->getBitDepths());

#if RDOQ_CHROMA_LAMBDA
        // set lambda for RDOQ
        const Double chromaLambda = estLambda / m_pcRdCost->getChromaWeight();
        const Double lambdaArray[MAX_NUM_COMPONENT] = { estLambda, chromaLambda, chromaLambda };
        m_pcTrQuant->setLambdas( lambdaArray );
#else
        m_pcTrQuant->setLambda( estLambda );
#endif
      }

      m_pcRateCtrl->setRCQP( estQP );
    }

#if BLOCK_BASED_BACKGROUND_REFERENCE
		if (pcSlice->getSPS()->getBackgroundRefsPresent())
		{
			if (pcSlice->getSliceType() != I_SLICE && pCtu->getSlice()->getQPList() != NULL)
			{
				Int iQPdif;
				if (ctuRsAddr == 0)
				{
					iQPdif = pCtu->getSlice()->getSliceQp() - pCtu->getSlice()->getQPList(0);
				}
				else
				{
					iQPdif = pCtu->getSlice()->getQPList(ctuRsAddr - 1) - pCtu->getSlice()->getQPList(ctuRsAddr);
				}
				Double dNewLambda = pCtu->getSlice()->getLambdas()[0] / pow(2.0, iQPdif / 3.0);
				resetQP(pcPic, pCtu->getSlice()->getQPList(ctuRsAddr), dNewLambda);
			}
		}
#endif

    // run CTU trial encoder
    m_pcCuEncoder->compressCtu( pCtu );

#if  BLOCK_BASED_BACKGROUND_REFERENCE
		if (pcSlice->getSPS()->getBackgroundRefsPresent())
		{
			if (pcSlice->getSliceType() != I_SLICE)
			{
				if (pCtu->getSlice()->getQPList() != NULL)
				{
					pCtu->setCodedQP(0);
				}
				if (pCtu->getSlice()->getSPS()->getPicBg()->getSpliceIdx(ctuRsAddr) == pCtu->getSlice()->getPOC())
				{
					pCtu->setLCUBgFlag(true);
				}
				else
				{
					pCtu->setLCUBgFlag(false);
				}
			}
		}
#endif

    // All CTU decisions have now been made. Restore entropy coder to an initial stage, ready to make a true encode,
    // which will result in the state of the contexts being correct. It will also count up the number of bits coded,
    // which is used if there is a limit of the number of bytes per slice-segment.

#if JVET_C0024_QTBT
    m_pcEntropyCoder->setEntropyCoder ( m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST] );
#else
    m_pcEntropyCoder->setEntropyCoder ( m_pppcRDSbacCoder[0][CI_CURR_BEST] );
#endif
    m_pcEntropyCoder->setBitstream( &tempBitCounter );
    pRDSbacCoder->setBinCountingEnableFlag( true );
#if JVET_C0024_QTBT
    m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->resetBits();
#else
    m_pppcRDSbacCoder[0][CI_CURR_BEST]->resetBits();
#endif
    pRDSbacCoder->setBinsCoded( 0 );

    // encode CTU and calculate the true bit counters.
    m_pcCuEncoder->encodeCtu( pCtu );


    pRDSbacCoder->setBinCountingEnableFlag( false );

    const Int numberOfWrittenBits = m_pcEntropyCoder->getNumberOfWrittenBits();

    // Calculate if this CTU puts us over slice bit size.
    // cannot terminate if current slice/slice-segment would be 0 Ctu in size,
    const UInt validEndOfSliceCtuTsAddr = ctuTsAddr + (ctuTsAddr == startCtuTsAddr ? 1 : 0);
    // Set slice end parameter
    if(pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES && pcSlice->getSliceBits()+numberOfWrittenBits > (pcSlice->getSliceArgument()<<3))
    {
      pcSlice->setSliceSegmentCurEndCtuTsAddr(validEndOfSliceCtuTsAddr);
      pcSlice->setSliceCurEndCtuTsAddr(validEndOfSliceCtuTsAddr);
      boundingCtuTsAddr=validEndOfSliceCtuTsAddr;
    }
    else if((!bCompressEntireSlice) && pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES && pcSlice->getSliceSegmentBits()+numberOfWrittenBits > (pcSlice->getSliceSegmentArgument()<<3))
    {
      pcSlice->setSliceSegmentCurEndCtuTsAddr(validEndOfSliceCtuTsAddr);
      boundingCtuTsAddr=validEndOfSliceCtuTsAddr;
    }

    if (boundingCtuTsAddr <= ctuTsAddr)
    {
      break;
    }

    pcSlice->setSliceBits( (UInt)(pcSlice->getSliceBits() + numberOfWrittenBits) );
    pcSlice->setSliceSegmentBits(pcSlice->getSliceSegmentBits()+numberOfWrittenBits);

    // Store probabilities of second CTU in line into buffer - used only if wavefront-parallel-processing is enabled.
    if ( ctuXPosInCtus == tileXPosInCtus+1 && m_pcCfg->getWaveFrontsynchro())
    {
#if JVET_C0024_QTBT
      m_entropyCodingSyncContextState.loadContexts(m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]);
#else
      m_entropyCodingSyncContextState.loadContexts(m_pppcRDSbacCoder[0][CI_CURR_BEST]);
#endif
    }


    if ( m_pcCfg->getUseRateCtrl() )
    {
      Int actualQP        = g_RCInvalidQPValue;
      Double actualLambda = m_pcRdCost->getLambda();
      Int actualBits      = pCtu->getTotalBits();
      Int numberOfEffectivePixels    = 0;
      for ( Int idx = 0; idx < pcPic->getNumPartitionsInCtu(); idx++ )
      {
        if ( pCtu->getPredictionMode( idx ) != NUMBER_OF_PREDICTION_MODES && ( !pCtu->isSkipped( idx ) ) )
        {
          numberOfEffectivePixels = numberOfEffectivePixels + 16;
          break;
        }
      }

      if ( numberOfEffectivePixels == 0 )
      {
        actualQP = g_RCInvalidQPValue;
      }
      else
      {
        actualQP = pCtu->getQP( 0 );
      }
      m_pcRdCost->setLambda(oldLambda, pcSlice->getSPS()->getBitDepths());
      m_pcRateCtrl->getRCPic()->updateAfterCTU( m_pcRateCtrl->getRCPic()->getLCUCoded(), actualBits, actualQP, actualLambda,
                                                pCtu->getSlice()->getSliceType() == I_SLICE ? 0 : m_pcCfg->getLCULevelRC() );
    }

    m_uiPicTotalBits += pCtu->getTotalBits();
    m_dPicRdCost     += pCtu->getTotalCost();
    m_uiPicDist      += pCtu->getTotalDistortion();
#if SHOW_MODE
#if BARC_REF_UP_SAMPLE_CNN
#if BARC_ONLY_FOR_INTRA
	if (pcSlice->getSliceType() == I_SLICE)
	{
		cout << '\n' << "ctuRsAddr: " << ctuRsAddr << ' ' << "BARC Flag: " << pCtu->getBARCFlag(CHANNEL_TYPE_LUMA) << " " << "CNN Up Flag: " << pCtu->getUpCNNFlag();
	}
#endif
#endif
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
		//reset QP for filter followed
		if (pcSlice->getSPS()->getBackgroundRefsPresent())
		{
			if (pcSlice->getSliceType() != I_SLICE && pcSlice->getQPList() != NULL)
			{
#if !QTBT_DQP_FIX
				if (pCtu->getCodedQP() == 0)
				{
					if (pcSlice->getSPS()->getPicBg()->getSpliceIdx(ctuRsAddr) == pcSlice->getPOC())
					{
						if (ctuRsAddr == 0)
						{
							Int partStartIdx = 0 - pcPic->getPicSym()->getCtuRsToTsAddrMap(ctuRsAddr) * pcPic->getNumPartitionsInCtu();
							Int firstElement = max<Int>(partStartIdx, 0);
							Int numElements = pcPic->getNumPartitionsInCtu() - firstElement;
							memset(pCtu->getQP() + firstElement, iOrgSliceQP, numElements * sizeof(*pCtu->getQP()));
						}
					}
				}
				else if (pcSlice->getSliceQp() > 0 && pcSlice->getSliceQp() != iOrgSliceQP)
				{
					Int partStartIdx = 0 - pcPic->getPicSym()->getCtuRsToTsAddrMap(ctuRsAddr) * pcPic->getNumPartitionsInCtu();
					Int firstElement = max<Int>(partStartIdx, 0);
					Int numElements = pcPic->getNumPartitionsInCtu() - firstElement;
					memset(pCtu->getQP() + firstElement, pcSlice->getSliceQp(), numElements * sizeof(*pCtu->getQP()));
					assert(pCtu->getLCUBgFlag());
				}

#else
				Int partStartIdx = 0 - pcPic->getPicSym()->getCtuRsToTsAddrMap(ctuRsAddr) * pcPic->getNumPartitionsInCtu();
				Int firstElement = max<Int>(partStartIdx, 0);
				Int  numElements = pcPic->getNumPartitionsInCtu() - firstElement;
				if (pCtu->getLCUBgFlag())
				{
					memset(pCtu->getQP() + firstElement, pcSlice->getSliceBgQp(), numElements * sizeof(*pCtu->getQP()));
				}
				else
				{
					memset(pCtu->getQP() + firstElement, iOrgSliceQP, numElements * sizeof(*pCtu->getQP()));
				}
#endif
			}
		}
#endif
  }

#if BLOCK_BASED_BACKGROUND_REFERENCE
	if (pcSlice->getSPS()->getBackgroundRefsPresent())
	{
		if (pcSlice->getSliceType() != I_SLICE && pcSlice->getQPList() != NULL)
		{
			if (pcSlice->getSliceQp() != iOrgSliceQP)
			{
				Double dOrgLambda = pcSlice->getLambdas()[0] / pow(2.0, (pcSlice->getSliceQp() - iOrgSliceQP) / 3.0);
				resetQP(pcPic, iOrgSliceQP, dOrgLambda);
			}
		}
	}
#endif

#if BARCPostUpSample
  for (UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < boundingCtuTsAddr; ++ctuTsAddr)
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    // initialize CTU encoder
    TComDataCU* pCtu = pcPic->getCtu(ctuRsAddr);

    if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
      xUpSampleRec(pCtu, CHANNEL_TYPE_LUMA);
  }
#endif
#if REFINE_BLOCK_BOUNDARY_NET
  if (pcSlice->getSliceType() == I_SLICE)
  {
	  pcPic->upSampleDScaleRecon();
  }
#endif
  // store context state at the end of this slice-segment, in case the next slice is a dependent slice and continues using the CABAC contexts.
  if( pcSlice->getPPS()->getDependentSliceSegmentsEnabledFlag() )
  {
#if JVET_C0024_QTBT
    m_lastSliceSegmentEndContextState.loadContexts( m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST] );//ctx end of dep.slice
#else
    m_lastSliceSegmentEndContextState.loadContexts( m_pppcRDSbacCoder[0][CI_CURR_BEST] );//ctx end of dep.slice
#endif
  }

  // stop use of temporary bit counter object.
#if JVET_C0024_QTBT
  m_ppppcRDSbacCoder[uiMaxWIdx][uiMaxHIdx][CI_CURR_BEST]->setBitstream(NULL);
#else
  m_pppcRDSbacCoder[0][CI_CURR_BEST]->setBitstream(NULL);
#endif
  m_pcRDGoOnSbacCoder->setBitstream(NULL); // stop use of tempBitCounter.

  // TODO: optimise cabac_init during compress slice to improve multi-slice operation
  //if (pcSlice->getPPS()->getCabacInitPresentFlag() && !pcSlice->getPPS()->getDependentSliceSegmentsEnabledFlag())
  //{
  //  m_encCABACTableIdx = m_pcEntropyCoder->determineCabacInitIdx();
  //}
  //else
  //{
  //  m_encCABACTableIdx = pcSlice->getSliceType();
  //}
}

Void TEncSlice::encodeSlice   ( TComPic* pcPic, TComOutputBitstream* pcSubstreams, UInt &numBinsCoded 
#if ALF_HM3_REFACTOR
  , ALFParam & alfParam
#endif
  )
{
#if TRACE_YUV || TRACE_FLAGS
  gIsEncodingNow = 1;
#endif
  TComSlice *const pcSlice           = pcPic->getSlice(getSliceIdx());

  const UInt startCtuTsAddr          = pcSlice->getSliceSegmentCurStartCtuTsAddr();
  const UInt boundingCtuTsAddr       = pcSlice->getSliceSegmentCurEndCtuTsAddr();

  const UInt frameWidthInCtus        = pcPic->getPicSym()->getFrameWidthInCtus();
  const Bool depSliceSegmentsEnabled = pcSlice->getPPS()->getDependentSliceSegmentsEnabledFlag();
  const Bool wavefrontsEnabled       = pcSlice->getPPS()->getEntropyCodingSyncEnabledFlag();

  // initialise entropy coder for the slice
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  TComStats*  pcStats = pcSlice->getStatsHandle();
  m_pcEntropyCoder->setStatsHandle ( pcStats );
#endif 
  m_pcSbacCoder->init( (TEncBinIf*)m_pcBinCABAC );
  m_pcEntropyCoder->setEntropyCoder ( m_pcSbacCoder );
  m_pcEntropyCoder->resetEntropy    ( pcSlice );

  numBinsCoded = 0;
  m_pcBinCABAC->setBinCountingEnableFlag( true );
  m_pcBinCABAC->setBinsCoded(0);

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

#if VCEG_AZ07_BAC_ADAPT_WDOW 
  Int iQPIdx = pcSlice->getQPIdx();
  TEncBinCABAC* pEncBin = m_pcEntropyCoder->getCABACCoder()->getBinIf()->getTEncBinCABAC();
  if(iQPIdx != -1)
  {
    pEncBin->allocateMemoryforBinStrings();
  }
#endif

  if (depSliceSegmentsEnabled)
  {
    // modify initial contexts with previous slice segment if this is a dependent slice.
    const UInt ctuRsAddr        = pcPic->getPicSym()->getCtuTsToRsAddrMap( startCtuTsAddr );
    const UInt currentTileIdx=pcPic->getPicSym()->getTileIdxMap(ctuRsAddr);
    const TComTile *pCurrentTile=pcPic->getPicSym()->getTComTile(currentTileIdx);
    const UInt firstCtuRsAddrOfTile = pCurrentTile->getFirstCtuRsAddr();

    if( pcSlice->getDependentSliceSegmentFlag() && ctuRsAddr != firstCtuRsAddrOfTile )
    {
      if( pCurrentTile->getTileWidthInCtus() >= 2 || !wavefrontsEnabled )
      {
        m_pcSbacCoder->loadContexts(&m_lastSliceSegmentEndContextState);
      }
    }
  }

  // for every CTU in the slice segment...

  for( UInt ctuTsAddr = startCtuTsAddr; ctuTsAddr < boundingCtuTsAddr; ++ctuTsAddr )
  {
    const UInt ctuRsAddr = pcPic->getPicSym()->getCtuTsToRsAddrMap(ctuTsAddr);
    const TComTile &currentTile = *(pcPic->getPicSym()->getTComTile(pcPic->getPicSym()->getTileIdxMap(ctuRsAddr)));
    const UInt firstCtuRsAddrOfTile = currentTile.getFirstCtuRsAddr();
    const UInt tileXPosInCtus       = firstCtuRsAddrOfTile % frameWidthInCtus;
    const UInt tileYPosInCtus       = firstCtuRsAddrOfTile / frameWidthInCtus;
    const UInt ctuXPosInCtus        = ctuRsAddr % frameWidthInCtus;
    const UInt ctuYPosInCtus        = ctuRsAddr / frameWidthInCtus;
    const UInt uiSubStrm=pcPic->getSubstreamForCtuAddr(ctuRsAddr, true, pcSlice);
    TComDataCU* pCtu = pcPic->getCtu( ctuRsAddr );

    m_pcEntropyCoder->setBitstream( &pcSubstreams[uiSubStrm] );

    // set up CABAC contexts' state for this CTU
    if (ctuRsAddr == firstCtuRsAddrOfTile)
    {
      if (ctuTsAddr != startCtuTsAddr) // if it is the first CTU, then the entropy coder has already been reset
      {
        m_pcEntropyCoder->resetEntropy(pcSlice);
      }
    }
    else if (ctuXPosInCtus == tileXPosInCtus && wavefrontsEnabled)
    {
      // Synchronize cabac probabilities with upper-right CTU if it's available and at the start of a line.
      if (ctuTsAddr != startCtuTsAddr) // if it is the first CTU, then the entropy coder has already been reset
      {
        m_pcEntropyCoder->resetEntropy(pcSlice);
      }
      TComDataCU *pCtuUp = pCtu->getCtuAbove();
      if ( pCtuUp && ((ctuRsAddr%frameWidthInCtus+1) < frameWidthInCtus)  )
      {
        TComDataCU *pCtuTR = pcPic->getCtu( ctuRsAddr - frameWidthInCtus + 1 );
        if ( pCtu->CUIsFromSameSliceAndTile(pCtuTR) )
        {
          // Top-right is available, so use it.
          m_pcSbacCoder->loadContexts( &m_entropyCodingSyncContextState );
        }
      }
    }
#if VCEG_AZ07_INIT_PREVFRAME
    if( pcSlice->getSliceType() != I_SLICE && ctuRsAddr == 0 )
    {
      m_pcSbacCoder->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdx(), true, pcSlice->getCtxMapQPIdxforStore(), (pcSlice->getPOC() > pcSlice->getStatsHandle()->m_uiLastIPOC) ); 
    }
#endif

#if ALF_HM3_REFACTOR
    if( pcSlice->getSPS()->getUseALF() && ctuRsAddr == 0 )
    {
      if (alfParam.cu_control_flag)
      {
        m_pcEntropyCoder->setAlfCtrl( true );
        m_pcEntropyCoder->setMaxAlfCtrlDepth(alfParam.alf_max_depth);
      }
      else
      {
        m_pcEntropyCoder->setAlfCtrl(false);
      }

      m_pcEntropyCoder->encodeAlfParam(&alfParam,pcSlice->getSPS()->getMaxTotalCUDepth()
#if FIX_TICKET12
        , pcSlice
#endif
        );
      if(alfParam.cu_control_flag)
      {
        m_pcEntropyCoder->encodeAlfCtrlParam(&alfParam);
      }
    }
#endif
#if CNN_BASED_INLOOP_FILTER
		if (pcSlice->getSPS()->getUseCNNLoopFilter())
		{
			CNNLFBlkParam* convNetBlkParam = (pcPic->getPicSym()->getConvNetBlkParam());
			for (Int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
			{
				ComponentID compId = ComponentID(comp);
				m_pcEntropyCoder->encodeCNNLFBlkParam(compId, convNetBlkParam, pCtu);
			}
		}
#endif
    if ( pcSlice->getSPS()->getUseSAO() )
    {
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
      }
      if (bIsSAOSliceEnabled)
      {
        SAOBlkParam& saoblkParam = (pcPic->getPicSym()->getSAOBlkParam())[ctuRsAddr];

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

        m_pcEntropyCoder->encodeSAOBlkParam(saoblkParam, pcPic->getPicSym()->getSPS().getBitDepths(), sliceEnabled, leftMergeAvail, aboveMergeAvail);
      }
    }

#if RRC_NS && !RRC_NS_DISABLE_SIGNALING
		if (pcSlice->getSPS()->getUseNS() && ctuTsAddr == 0)
		{
			const bool isAdaptiveNs = pcSlice->getPic()->getPicSym()->getAdaptiveNs();
			m_pcEntropyCoder->encodeAdaptiveNsFlag(isAdaptiveNs);
			if (isAdaptiveNs)
			{
				const int ns_block_size = pcPic->getPicSym()->getNsBlockSize();
#if JVET_C0024_QTBT
				const int block_width = pcPic->getSlice(0)->getSPS()->getCTUSize() >> ns_block_size;
				const int block_height = pcPic->getSlice(0)->getSPS()->getCTUSize() >> ns_block_size;
#else
				const int block_width = pcPic->getSlice(0)->getSPS()->getMaxCUWidth() >> ns_block_size;
				const int block_height = pcPic->getSlice(0)->getSPS()->getMaxCUHeight() >> ns_block_size;
#endif
				const int width_in_blocks = pcPic->getPicYuvRec()->getWidth(COMPONENT_Y) / block_width + (pcPic->getPicYuvRec()->getWidth(COMPONENT_Y) % block_width != 0);
				const int height_in_blocks = pcPic->getPicYuvRec()->getHeight(COMPONENT_Y) / block_height + (pcPic->getPicYuvRec()->getHeight(COMPONENT_Y) % block_height != 0);

				m_pcEntropyCoder->encodeNsBlockSize(ns_block_size);
				UInt* nsMap = pcSlice->getPic()->getPicSym()->getNsFlag();
				for (int i = 0; i < width_in_blocks * height_in_blocks; ++i)
					m_pcEntropyCoder->encodeNsFlag(nsMap[i]);
			}
			else
			{
				const bool isNsWholeFrame = pcSlice->getPic()->getPicSym()->getNsWholeFrame();
				m_pcEntropyCoder->encodeNsFWholeFrameFlag(isNsWholeFrame);
			}
		}
#endif

#if ENC_DEC_TRACE
    g_bJustDoIt = g_bEncDecTraceEnable;
#endif
      m_pcCuEncoder->encodeCtu( pCtu );
#if ENC_DEC_TRACE
    g_bJustDoIt = g_bEncDecTraceDisable;
#endif

    //Store probabilities of second CTU in line into buffer
    if ( ctuXPosInCtus == tileXPosInCtus+1 && wavefrontsEnabled)
    {
      m_entropyCodingSyncContextState.loadContexts( m_pcSbacCoder );
    }

#if VCEG_AZ07_INIT_PREVFRAME
    UInt uiTargetCUAddr = pcPic->getFrameWidthInCtus()/2 + pcPic->getNumberOfCtusInFrame()/2;
    if( uiTargetCUAddr >= pcPic->getNumberOfCtusInFrame() )
    {
      uiTargetCUAddr = pcPic->getNumberOfCtusInFrame() - 1;
    }
    if( pcSlice->getSliceType() != I_SLICE && ctuTsAddr ==  uiTargetCUAddr )
    {
      m_pcSbacCoder->loadContextsFromPrev( pcSlice->getStatsHandle(), pcSlice->getSliceType(), pcSlice->getCtxMapQPIdxforStore(), false ); 
    }
#endif

    // terminate the sub-stream, if required (end of slice-segment, end of tile, end of wavefront-CTU-row):
    if (ctuTsAddr+1 == boundingCtuTsAddr ||
         (  ctuXPosInCtus + 1 == tileXPosInCtus + currentTile.getTileWidthInCtus() &&
          ( ctuYPosInCtus + 1 == tileYPosInCtus + currentTile.getTileHeightInCtus() || wavefrontsEnabled)
         )
       )
    {
      m_pcEntropyCoder->encodeTerminatingBit(1);
      m_pcEntropyCoder->encodeSliceFinish();
      // Byte-alignment in slice_data() when new tile
      pcSubstreams[uiSubStrm].writeByteAlignment();

      // write sub-stream size
      if (ctuTsAddr+1 != boundingCtuTsAddr)
      {
        pcSlice->addSubstreamSize( (pcSubstreams[uiSubStrm].getNumberOfWrittenBits() >> 3) + pcSubstreams[uiSubStrm].countStartCodeEmulations() );
      }
    }
  } // CTU-loop

  if( depSliceSegmentsEnabled )
  {
    m_lastSliceSegmentEndContextState.loadContexts( m_pcSbacCoder );//ctx end of dep.slice
  }

#if ADAPTIVE_QP_SELECTION
  if( m_pcCfg->getUseAdaptQpSelect() )
  {
    m_pcTrQuant->storeSliceQpNext(pcSlice); // TODO: this will only be storing the adaptive QP state of the very last slice-segment that is not dependent in the frame... Perhaps this should be moved to the compress slice loop.
  }
#endif

#if PARALLEL_ENCODING_RAS_CABAC_INIT_PRESENT  
  m_eLastNALUType = pcSlice->getNalUnitType();
#endif

  if (pcSlice->getPPS()->getCabacInitPresentFlag() && !pcSlice->getPPS()->getDependentSliceSegmentsEnabledFlag())
  {
    m_encCABACTableIdx = m_pcEntropyCoder->determineCabacInitIdx(pcSlice);
  }
  else
  {
    m_encCABACTableIdx = pcSlice->getSliceType();
  }
  
  numBinsCoded = m_pcBinCABAC->getBinsCoded();

#if VCEG_AZ07_BAC_ADAPT_WDOW 
  if(iQPIdx != -1)
  {
    xGenUpdateMap    (pcSlice->getSliceType(), pcSlice->getQPIdx(), pcStats);
    pEncBin->freeMemoryforBinStrings();
  }
#endif
#if TRACE_YUV
  DUMP_ALL_YUV();
#endif
#if TRACE_YUV || TRACE_FLAGS
  gIsEncodingNow = 0;
#endif
#if VCEG_AZ07_INIT_PREVFRAME_FIX
  pcSlice->updateStatsGlobal();
#endif
}

Void TEncSlice::calculateBoundingCtuTsAddrForSlice(UInt &startCtuTSAddrSlice, UInt &boundingCtuTSAddrSlice, Bool &haveReachedTileBoundary,
                                                   TComPic* pcPic, const Int sliceMode, const Int sliceArgument)
{
  TComSlice* pcSlice = pcPic->getSlice(getSliceIdx());
  const UInt numberOfCtusInFrame = pcPic->getNumberOfCtusInFrame();
  const TComPPS &pps=*(pcSlice->getPPS());
  boundingCtuTSAddrSlice=0;
  haveReachedTileBoundary=false;

  switch (sliceMode)
  {
    case FIXED_NUMBER_OF_CTU:
      {
        UInt ctuAddrIncrement    = sliceArgument;
        boundingCtuTSAddrSlice  = ((startCtuTSAddrSlice + ctuAddrIncrement) < numberOfCtusInFrame) ? (startCtuTSAddrSlice + ctuAddrIncrement) : numberOfCtusInFrame;
      }
      break;
    case FIXED_NUMBER_OF_BYTES:
      boundingCtuTSAddrSlice  = numberOfCtusInFrame; // This will be adjusted later if required.
      break;
    case FIXED_NUMBER_OF_TILES:
      {
        const UInt tileIdx        = pcPic->getPicSym()->getTileIdxMap( pcPic->getPicSym()->getCtuTsToRsAddrMap(startCtuTSAddrSlice) );
        const UInt tileTotalCount = (pcPic->getPicSym()->getNumTileColumnsMinus1()+1) * (pcPic->getPicSym()->getNumTileRowsMinus1()+1);
        UInt ctuAddrIncrement   = 0;

        for(UInt tileIdxIncrement = 0; tileIdxIncrement < sliceArgument; tileIdxIncrement++)
        {
          if((tileIdx + tileIdxIncrement) < tileTotalCount)
          {
            UInt tileWidthInCtus   = pcPic->getPicSym()->getTComTile(tileIdx + tileIdxIncrement)->getTileWidthInCtus();
            UInt tileHeightInCtus  = pcPic->getPicSym()->getTComTile(tileIdx + tileIdxIncrement)->getTileHeightInCtus();
            ctuAddrIncrement    += (tileWidthInCtus * tileHeightInCtus);
          }
        }

        boundingCtuTSAddrSlice  = ((startCtuTSAddrSlice + ctuAddrIncrement) < numberOfCtusInFrame) ? (startCtuTSAddrSlice + ctuAddrIncrement) : numberOfCtusInFrame;
      }
      break;
    default:
      boundingCtuTSAddrSlice    = numberOfCtusInFrame;
      break;
  }

  // Adjust for tiles and wavefronts.
  const Bool wavefrontsAreEnabled = pps.getEntropyCodingSyncEnabledFlag();

  if ((sliceMode == FIXED_NUMBER_OF_CTU || sliceMode == FIXED_NUMBER_OF_BYTES) &&
      (pps.getNumTileRowsMinus1() > 0 || pps.getNumTileColumnsMinus1() > 0))
  {
    const UInt ctuRSAddr                  = pcPic->getPicSym()->getCtuTsToRsAddrMap(startCtuTSAddrSlice);
    const UInt startTileIdx               = pcPic->getPicSym()->getTileIdxMap(ctuRSAddr);

    const TComTile *pStartingTile         = pcPic->getPicSym()->getTComTile(startTileIdx);
    const UInt tileStartTsAddr            = pcPic->getPicSym()->getCtuRsToTsAddrMap(pStartingTile->getFirstCtuRsAddr());
    const UInt tileStartWidth             = pStartingTile->getTileWidthInCtus();
    const UInt tileStartHeight            = pStartingTile->getTileHeightInCtus();
    const UInt tileLastTsAddr_excl        = tileStartTsAddr + tileStartWidth*tileStartHeight;
    const UInt tileBoundingCtuTsAddrSlice = tileLastTsAddr_excl;

    const UInt ctuColumnOfStartingTile    = ((startCtuTSAddrSlice-tileStartTsAddr)%tileStartWidth);
    if (wavefrontsAreEnabled && ctuColumnOfStartingTile!=0)
    {
      // WPP: if a slice does not start at the beginning of a CTB row, it must end within the same CTB row
      const UInt numberOfCTUsToEndOfRow            = tileStartWidth - ctuColumnOfStartingTile;
      const UInt wavefrontTileBoundingCtuAddrSlice = startCtuTSAddrSlice + numberOfCTUsToEndOfRow;
      if (wavefrontTileBoundingCtuAddrSlice < boundingCtuTSAddrSlice)
      {
        boundingCtuTSAddrSlice = wavefrontTileBoundingCtuAddrSlice;
      }
    }

    if (tileBoundingCtuTsAddrSlice < boundingCtuTSAddrSlice)
    {
      boundingCtuTSAddrSlice = tileBoundingCtuTsAddrSlice;
      haveReachedTileBoundary = true;
    }
  }
  else if ((sliceMode == FIXED_NUMBER_OF_CTU || sliceMode == FIXED_NUMBER_OF_BYTES) && wavefrontsAreEnabled && ((startCtuTSAddrSlice % pcPic->getFrameWidthInCtus()) != 0))
  {
    // Adjust for wavefronts (no tiles).
    // WPP: if a slice does not start at the beginning of a CTB row, it must end within the same CTB row
    boundingCtuTSAddrSlice = min(boundingCtuTSAddrSlice, startCtuTSAddrSlice - (startCtuTSAddrSlice % pcPic->getFrameWidthInCtus()) + (pcPic->getFrameWidthInCtus()));
  }
}

/** Determines the starting and bounding CTU address of current slice / dependent slice
 * \param [out] startCtuTsAddr
 * \param [out] boundingCtuTsAddr
 * \param [in]  pcPic

 * Updates startCtuTsAddr, boundingCtuTsAddr with appropriate CTU address
 */
Void TEncSlice::xDetermineStartAndBoundingCtuTsAddr  ( UInt& startCtuTsAddr, UInt& boundingCtuTsAddr, TComPic* pcPic )
{
  TComSlice* pcSlice                 = pcPic->getSlice(getSliceIdx());

  // Non-dependent slice
  UInt startCtuTsAddrSlice           = pcSlice->getSliceCurStartCtuTsAddr();
  Bool haveReachedTileBoundarySlice  = false;
  UInt boundingCtuTsAddrSlice;
  calculateBoundingCtuTsAddrForSlice(startCtuTsAddrSlice, boundingCtuTsAddrSlice, haveReachedTileBoundarySlice, pcPic,
                                     m_pcCfg->getSliceMode(), m_pcCfg->getSliceArgument());
  pcSlice->setSliceCurEndCtuTsAddr(   boundingCtuTsAddrSlice );
  pcSlice->setSliceCurStartCtuTsAddr( startCtuTsAddrSlice    );

  // Dependent slice
  UInt startCtuTsAddrSliceSegment          = pcSlice->getSliceSegmentCurStartCtuTsAddr();
  Bool haveReachedTileBoundarySliceSegment = false;
  UInt boundingCtuTsAddrSliceSegment;
  calculateBoundingCtuTsAddrForSlice(startCtuTsAddrSliceSegment, boundingCtuTsAddrSliceSegment, haveReachedTileBoundarySliceSegment, pcPic,
                                     m_pcCfg->getSliceSegmentMode(), m_pcCfg->getSliceSegmentArgument());
  if (boundingCtuTsAddrSliceSegment>boundingCtuTsAddrSlice)
  {
    boundingCtuTsAddrSliceSegment = boundingCtuTsAddrSlice;
  }
  pcSlice->setSliceSegmentCurEndCtuTsAddr( boundingCtuTsAddrSliceSegment );
  pcSlice->setSliceSegmentCurStartCtuTsAddr(startCtuTsAddrSliceSegment);

  // Make a joint decision based on reconstruction and dependent slice bounds
  startCtuTsAddr    = max(startCtuTsAddrSlice   , startCtuTsAddrSliceSegment   );
  boundingCtuTsAddr = boundingCtuTsAddrSliceSegment;
}

Double TEncSlice::xGetQPValueAccordingToLambda ( Double lambda )
{
  return 4.2005*log(lambda) + 13.7122;
}

#if VCEG_AZ07_BAC_ADAPT_WDOW 
#if VCEG_AZ07_INIT_PREVFRAME
Void TEncSlice::xContextWdowSizeUpdateDecision(TEncSbac* pTestEncSbac, UInt &uiCtxStartPos, ContextModel* pSliceCtx, Bool *uiCtxMap, UChar *uiCtxCodeIdx, Bool** pCodedBinStr, Int* pCounter, UShort* uiCTX)
#else
Void TEncSlice::xContextWdowSizeUpdateDecision(TEncSbac* pTestEncSbac, UInt &uiCtxStartPos, ContextModel* pSliceCtx, Bool *uiCtxMap, UChar *uiCtxCodeIdx, Bool** pCodedBinStr, Int* pCounter)
#endif
{
  //derive the best window size
  UInt uiBestW = 0, currBits = 0, minBits = MAX_UINT;
  Bool* pBinStr  = pCodedBinStr[uiCtxStartPos];
  UInt  uiStrlen = pCounter[uiCtxStartPos];
  if(uiStrlen==0)
  {
    uiCtxMap[uiCtxStartPos] = 0;
    uiCtxStartPos += 1;
    return;
  }


  for(UInt uiW = 0; uiW < NUM_WDOW; uiW++)
  {
   //initilization        
    pTestEncSbac->getEncBinIf()->start(); 
    currBits = pTestEncSbac->getNumberOfWrittenBits();

    ContextModel cCurrCtx = pSliceCtx[uiCtxStartPos];  
#if VCEG_AZ07_INIT_PREVFRAME
    if(uiCTX != NULL)
    {
      cCurrCtx.setState(uiCTX[uiCtxStartPos]);
    }
#endif

    cCurrCtx.setWindowSize(uiW + 4);    


    for(UInt k = 0; k < uiStrlen; k ++)
    {         
      UInt val = pBinStr[k];
      pTestEncSbac->getEncBinIf()->encodeBin(val, cCurrCtx);
    }
    currBits = pTestEncSbac->getNumberOfWrittenBits()-currBits;

    if(currBits < minBits)
    {
      minBits = currBits;
      uiBestW = uiW;
    }   
  }

  if( uiBestW != (ALPHA0 - 4) )
  {
    uiCtxMap[uiCtxStartPos] = 1;
    uiCtxCodeIdx[uiCtxStartPos] =(UChar) (uiBestW + 4); //the best window size, to be signalled
  }
  else
  {
    uiCtxMap[uiCtxStartPos] = 0;
    uiCtxCodeIdx[uiCtxStartPos] = ALPHA0;
  }
  uiCtxStartPos += 1;
}

Void TEncSlice::xGenUpdateMap (UInt uiSliceType, Int iQP,  TComStats* apcStats)
{
  UInt uiCtxStartPos = 0;
  Bool bUpdate       = false;

  TEncSbac*     pSbacCoder   =  m_pcEntropyCoder->getCABACCoder();  
  TEncBinCABAC* pEncBinCABAC = pSbacCoder->getBinIf()->getTEncBinCABAC();
  Bool** pCodedBinStr        = pEncBinCABAC->getCodedBinStr();
  Int*   pCounter            = pEncBinCABAC->getCodedNumBins();
  Int    iCtxNr              = pSbacCoder->getCtxNumber();

  ContextModel* pCtx         =  pSbacCoder->getContextModel();
  pEncBinCABAC->setUpdateStr(false);

  TEncSbac*     pTestEncSbac        = new TEncSbac;  
  TEncBinCABAC* pcTestBinCoderCABAC = new TEncBinCABAC;
  TComOutputBitstream* pBitIf       = new TComOutputBitstream;

  pTestEncSbac->init(pcTestBinCoderCABAC);
  pTestEncSbac->setBitstream(pBitIf);

  for(UInt i=0; i<iCtxNr; i++)
  {
#if VCEG_AZ07_INIT_PREVFRAME
    xContextWdowSizeUpdateDecision(pTestEncSbac, uiCtxStartPos, pCtx, apcStats->m_uiCtxMAP[uiSliceType][iQP], apcStats->m_uiCtxCodeIdx[uiSliceType][iQP], pCodedBinStr, pCounter, (uiSliceType != I_SLICE ?apcStats->m_uiCtxProbIdx[uiSliceType][iQP][0]: NULL));
#else
    xContextWdowSizeUpdateDecision(pTestEncSbac, uiCtxStartPos, pCtx, apcStats->m_uiCtxMAP[uiSliceType][iQP], apcStats->m_uiCtxCodeIdx[uiSliceType][iQP], pCodedBinStr, pCounter);
#endif

    if(bUpdate==false && apcStats->m_uiCtxMAP[uiSliceType][iQP][i])
    {
      bUpdate = true;
    }
  }
  apcStats->m_uiNumCtx[uiSliceType][iQP] = uiCtxStartPos; 
  assert( uiCtxStartPos == iCtxNr );

  if(bUpdate==false)
  {
    apcStats->m_uiNumCtx[uiSliceType][iQP]= 0; 
  }

  if(pTestEncSbac)        {delete pTestEncSbac;             pTestEncSbac        = NULL;      }
  if(pcTestBinCoderCABAC) {delete pcTestBinCoderCABAC;      pcTestBinCoderCABAC = NULL;      }
  if(pBitIf)              {delete pBitIf;                   pBitIf              = NULL;      }

}
#endif

#if PVC_AQS
//analyze the original input picture to obtain a error sensitivity map for each pixel
void TEncSlice::initErrorSensitivityMap(TComPic*  pcPic)
{
  Int i, j, k, ind;
  Pel* pPel_pic = pcPic->getPicYuvOrg()->getAddr(COMPONENT_Y);
  Int width = pcPic->getPicYuvOrg()->getWidth(COMPONENT_Y);
  Int height = pcPic->getPicYuvOrg()->getHeight(COMPONENT_Y);
  Int stride = pcPic->getPicYuvOrg()->getStride(COMPONENT_Y);
  Int s = width*height;
  UInt fac_lumi = 1; //luminance
  UInt fac_cont = 1; //contrast

  Int blk_w = 8;//PVC_ENCODER_ESMAP_BLK_SIZE;
  Int blk_h = 8;//PVC_ENCODER_ESMAP_BLK_SIZE;
  Int bg_w = 8;//PVC_ENCODER_ESMAP_BG_SIZE;
  Int bg_h = 8;//PVC_ENCODER_ESMAP_BG_SIZE;

  Int vd = 3;
  Double bg_angle = 0.5; //0.5 visual degree
  Double pi = 3.1415926535;
  UChar* puc_bg = NULL;
  UChar* puc_blk = NULL;
  Int blk_size = 0, bg_size = 0;
  Int i_offset = (bg_w - blk_w) / 2;
  Int j_offset = (bg_h - blk_h) / 2;
  Int a, b;
  Int avg_bg, std_bg;
  ESM fac_all_mul = 0;
  ESM pesm_weight_blk;
  Int sum = 0;

  // normalizing the ES value into [0.5, 2.0]
  ESM up = 1 << (ESM_SHIFT + 1); //MAX_ESM
  ESM low = 1 << (ESM_SHIFT - 1); //MIN_ESM

  //--------------- modify Upper bound for high bitrate case (QP<28) - in ESMap
  Int iQP = pcPic->getSlice(0)->getSliceQp();
  Int iNormFactor = 0;
  Int iFactor = 0;

  Int iOffsetQPIslice = 28;
  Double dScalerIsliceStep = 0.10;
  Double dScalerIsliceMax = 1.00;
  Double dScalerIsliceMin = 0.25;

  Int iOffsetQPPBslice = 28;
  Double dScalerPBsliceStep = 0.10;
  Double dScalerPBsliceMax = 1.00;
  Double dScalerPBsliceMin = 0.25;

  Double dMaxScaleFactor = 0;

  if (pcPic->getSlice(0)->getSliceType() == I_SLICE)
  {
    dMaxScaleFactor = Clip3(dScalerIsliceMin, dScalerIsliceMax, (iQP - iOffsetQPIslice) * dScalerIsliceStep + 1.0) * (2.0 - 1.0) + 1.0;
    up = (Int)((1 << ESM_SHIFT) * dMaxScaleFactor + 0.5);
  }
  else
  {
    dMaxScaleFactor = Clip3(dScalerPBsliceMin, dScalerPBsliceMax, (iQP - iOffsetQPPBslice) * dScalerPBsliceStep + 1.0) * (2.0 - 1.0) + 1.0;
    up = (Int)((1 << ESM_SHIFT) * dMaxScaleFactor + 0.5);
  }

  //initial look-up tables for factors
  initFactorValue();

  //allocate memory
  blk_size = blk_w * blk_h;
  bg_size = bg_w  * bg_h;
  puc_blk = (UChar*)malloc(blk_size*sizeof(UChar));
  puc_bg = (UChar*)malloc(bg_size *sizeof(UChar));

  for (j = 0; j<height; j = j + blk_h)
  {
    k = j*width;
    for (i = 0; i<width; i = i + blk_w)
    {
      ind = k + i;
      // prepare blk and its background
      fetch_block(pPel_pic, puc_blk, width, height, i, j, blk_w, blk_h, stride);
      if (blk_w == bg_w)
        fetch_block(pPel_pic, puc_bg, width, height, i, j, bg_w, bg_h, stride);
      else if (bg_w > blk_w)
      {
        a = i - i_offset;
        b = j - j_offset;
        if (a >= 0 && b >= 0)
          fetch_block(pPel_pic, puc_bg, width, height, a, b, bg_w, bg_h, stride);
        else //approximation for blk at image margin
          fetch_block(pPel_pic, puc_bg, width, height, (i / bg_w)*bg_w, (j / bg_h)*bg_h, bg_w, bg_h, stride);
      }

      // get average luminance and variance of the blk and the background
      // note:avg and var can be calcuated in one loop, since var = E(x^2) - avg^2;
      avg_bg = get_block_avg(puc_bg, bg_w, bg_h);
      std_bg = get_block_std_simp(puc_bg, bg_w, bg_h, avg_bg);

      // shift according to bit-depth
      avg_bg >>= (pcPic->getSlice(0)->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8);
      std_bg >>= (pcPic->getSlice(0)->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8);

      //luminance
      fac_lumi = g_aiFactorLuminance[avg_bg];
      //contrast
      fac_cont = g_aiFactorContrast[std_bg];
      //perceptual ESM weight of the blk
      pesm_weight_blk = Clip3(low, up, (fac_lumi * fac_cont + (1 << (ESM_SHIFT - 2))) >> (ESM_SHIFT - 1));

      assign_block(m_pesmErrorSensitivityMap, pesm_weight_blk, width, height, i, j, blk_w, blk_h);
    }
  }

  // modify ESMap by considering overall seq feature for bitrate alignment
  // it is of low efficiency now, needs optimization
  Int64 ESMAvg = 0;
  Int   ESM_offset = 0;
  Int Num = ((height + (blk_h - 1)) / blk_h) * ((width + (blk_w - 1)) / blk_w);
  Int ESM_mul_factor = 0;
  ESM val = 0;

  Double ESMGeoMean = 0;
  for (j = 0; j<height; j = j + blk_h)
  {
    k = j*width;
    for (i = 0; i<width; i = i + blk_w)
    {
      ESMGeoMean += log((Double)m_pesmErrorSensitivityMap[k + i]);
    }
  }
  ESMGeoMean = exp(ESMGeoMean / (Double)Num);
  ESMAvg = (ESM)ESMGeoMean;

  if (pcPic->getSlice(0)->getSliceType() == I_SLICE)
  {
    ESMAvg = (ESMAvg * 5 + (1 << ESM_SHIFT) * 3) >> 3;
  }
  else
  {
    ESMAvg = (ESMAvg * 3 + (1 << ESM_SHIFT) * 5) >> 3;
  }

  ESM_offset = (Int)(ESMAvg - (1 << ESM_SHIFT));

  //***************** quantize, since this value will be transmitted but it does not need high precision ***********
  ESM tabel_th_pos[8] = { 279, 304, 331, 362, 394, 430, 469, 512 };
  ESM tabel_th_neg[8] = { 234, 215, 197, 181, 165, 152, 139, 128 };
  ESM table_factor_pos[8] = { 256, 234, 215, 197, 181, 165, 152, 139 }; //also 8-bit 
  ESM table_factor_neg[8] = { 256, 279, 304, 331, 362, 394, 430, 469 };
  int clipVal = 6;

  // get factor index
  if (ESM_offset >= 0)
  {
    g_ESMFactorIdx = clipVal;
    for (Int i = 0; i <= clipVal; i++)
    {
      if (ESMAvg < tabel_th_pos[i])
      {
        g_ESMFactorIdx = i;
        break;
      }
    }
  }
  else
  {
    g_ESMFactorIdx = -clipVal;
    for (Int i = 0; i <= clipVal; i++)
    {
      if (ESMAvg > tabel_th_neg[i])
      {
        g_ESMFactorIdx = -i;
        break;
      }
    }
  }

  // get multiplication factor
  if (g_ESMFactorIdx >= 0)
    g_ESMFactor = table_factor_pos[g_ESMFactorIdx];
  else
    g_ESMFactor = table_factor_neg[-g_ESMFactorIdx];

  // modify ESM values
  for (j = 0; j<height; j++)
  {
    k = j*width;
    for (i = 0; i<width; i++)
    {
      val = m_pesmErrorSensitivityMap[k + i];
      m_pesmErrorSensitivityMap[k + i] = Clip3(low, up, (val * g_ESMFactor + 128) >> 8);
    }
  }

  if (puc_bg)
    free(puc_bg);
  if (puc_blk)
    free(puc_blk);
  puc_bg = puc_blk = NULL;
}
#endif

#if VCEG_AZ08_INTER_KLT
Void TEncSlice::InterpolatePic(TComPic* pcPic)
{
  m_pcPredSearch->interpolatePic(pcPic);
}
#endif
#if BARC_FUNC
Void TEncSlice::xDownScalePic(TComPicYuv* pcYuvSrc, TComPicYuv* pcYuvDest)
{
  pcYuvSrc->setBorderExtension(false);
  pcYuvSrc->extendPicBorder();
  downScaleComponent2x2ZeroPhase( pcYuvSrc->getAddr(COMPONENT_Y), pcYuvDest->getAddr(COMPONENT_Y), pcYuvSrc->getStride(COMPONENT_Y), pcYuvDest->getStride(COMPONENT_Y), pcYuvSrc->getWidth(COMPONENT_Y), pcYuvSrc->getHeight(COMPONENT_Y), COMPONENT_Y);
  pcYuvDest->setBorderExtension(false);
  pcYuvDest->extendPicBorder();
}

static int dn_filt_taps0[4][12] = { {36, 0, -4, 0, 0, 0}, /* 3-lobe Lanczos*/
{40, 0, -12, 0, 4, 0}, /*h11*/
{26, 19, 5, -3, -4, 0 ,2},// wrong
{39, -19, -10, 4, 1}};
static int dn_filt_half_lens0[4] = {3, 5, 7,6};
static int dn_centerTap[4] = {64, 64, 64, 128};
static int dn_filt_shifts[4] = {128, 128, 128, 158};

/** Downscale a component array by 2/2 with samples at phase 0 (co-sited) with respect to the high-res grid
 * \returns void
 */
Void  TEncSlice:: downScaleComponent2x2ZeroPhase( const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight, ComponentID compID)
{
  const int iDestWidth = iSrcWidth/2;
  const int iDestHeight = iSrcHeight/2;

  const int HALF_NUM_DOWN_TAPS = dn_filt_half_lens0[BARC_REF_DN_FILT];
  const int* iTaps = dn_filt_taps0[BARC_REF_DN_FILT];
  const int iFilterShift = dn_filt_shifts[BARC_REF_DN_FILT];
  const int iCentreTap = dn_centerTap[BARC_REF_DN_FILT];
  Pel pDstTemp[7];
  Int iTemp;

  // Vertical
//#if DownBottom
  const Pel* pDataSrc = pSrc + iSrcStride; 
//#else
 // const Pel* pDataSrc = pSrc; 
//#endif
  Pel* pDataDest = pDest;

  for (int j=0; j<iDestHeight; ++j) 
  {
    for (int i=-HALF_NUM_DOWN_TAPS; i<iSrcWidth+HALF_NUM_DOWN_TAPS; ++i)
    {
      iTemp = pDataSrc[i]*iCentreTap;
      for (int k=0; k<HALF_NUM_DOWN_TAPS; ++k)
      {
          iTemp += iTaps[k]*(pDataSrc[i+(k+1)*iSrcStride] + pDataSrc[i-(k+1)*iSrcStride]);
      }
      iTemp += (iFilterShift>>1); 
      pDataDest[i+HALF_NUM_DOWN_TAPS] = ClipA(iTemp/iFilterShift, compID);
    }
   // std::fill(pDataDest - HALF_NUM_DOWN_TAPS, pDataDest, pDataDest[0] );
   // std::fill(pDataDest + iSrcWidth, pDataDest+iSrcWidth+HALF_NUM_DOWN_TAPS, pDataDest[iSrcWidth-1] );
    pDataSrc += 2*iSrcStride;
    pDataDest += iDestStride;
  }

  // Horizontal
  Pel* pDataSrcTmp = pDest+HALF_NUM_DOWN_TAPS;
  pDataDest = pDstTemp;
  Int tmpIdx = 0;

  for (int j = 0; j < iDestHeight; ++j)
  {
    tmpIdx = 0;
    //#if DownRight
    for (int i = 0, ii = 1; i < iDestWidth; ++i, ii += 2)
      //#else
      //    for (int i=0,ii=0; i<iDestWidth; ++i,ii+=2) 
      //#endif
    {
      iTemp = pDataSrcTmp[ii] * iCentreTap;
      for (int k = 0; k < HALF_NUM_DOWN_TAPS; ++k)
      {
        iTemp += iTaps[k] * (pDataSrcTmp[ii + k + 1] + pDataSrcTmp[ii - k - 1]);
      }
      iTemp += (iFilterShift >> 1);
      if (tmpIdx < HALF_NUM_DOWN_TAPS)
      {
        pDataDest[tmpIdx++] = ClipA(iTemp / iFilterShift, compID);
      }
      else
      {
        for (Int tmpNum = 0; tmpNum < HALF_NUM_DOWN_TAPS; tmpNum++)
        {
          pDataSrcTmp[i - 2*HALF_NUM_DOWN_TAPS + tmpNum] = pDataDest[tmpNum];
        }
        tmpIdx = 0;
        pDataDest[tmpIdx++] = ClipA(iTemp / iFilterShift, compID);
      }
    }
    for (Int tmpNum = 0; tmpNum < tmpIdx; tmpNum++)
    {
      pDataSrcTmp[iDestWidth - HALF_NUM_DOWN_TAPS - tmpIdx + tmpNum] = pDataDest[tmpNum];
    }
    memset(pDataSrcTmp - HALF_NUM_DOWN_TAPS + iDestWidth, 0, sizeof(Pel)*(iSrcWidth - iDestWidth));
    pDataSrcTmp += iDestStride;
  }
} 
#endif
#if BARCPostUpSample
static int up_filt_taps0[4][4] = { { 40, -11, 4,  -1 },            /*f7*/
{ 36, -4, 0, 0 },           /*Linear*/
{ 20, -5, 1, 0 },
{ 20, -4, 0, 0 } };
static int up_filt_half_lens0[4] = { 4, 2, 3, 2 };
static int up_filt_shifts0[4] = { 6, 6, 5, 5 };
Void TEncSlice::xUpSampleRec(TComDataCU*& pCtu, ChannelType chType)
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
