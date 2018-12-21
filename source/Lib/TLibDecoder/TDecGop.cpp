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

/** \file     TDecGop.cpp
    \brief    GOP decoder class
*/

#include "TDecGop.h"
#include "TDecCAVLC.h"
#include "TDecSbac.h"
#include "TDecBinCoder.h"
#include "TDecBinCoderCABAC.h"
#include "libmd5/MD5.h"
#include "TLibCommon/SEI.h"
#if RRC_SELECTABLE_IF
#include "TLibCommon/TComSharpeningFilter.h"
#endif

#include <time.h>
#if RRC_NS
#include "TLibCommon/TComNoiseSuppressor.h"
#endif

//! \ingroup TLibDecoder
//! \{
static Void calcAndPrintHashStatus(TComPicYuv& pic, const SEIDecodedPictureHash* pictureHashSEI, const BitDepths &bitDepths, UInt &numChecksumErrors);
// ====================================================================================================================
// Constructor / destructor / initialization / destroy
// ====================================================================================================================

TDecGop::TDecGop()
 : m_numberOfChecksumErrorsDetected(0)
{
  m_dDecTime = 0;
}

TDecGop::~TDecGop()
{

}

Void TDecGop::create()
{

}


Void TDecGop::destroy()
{
#if COM16_C806_ALF_TEMPPRED_NUM
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
#endif
}

#if COM16_C806_ALF_TEMPPRED_NUM
#if JVET_E0104_ALF_TEMP_SCALABILITY
Int TDecGop::m_iStoredAlfParaNum[JVET_E0104_ALF_MAX_TEMPLAYERID] = { 0, 0, 0, 0, 0 };
#else
Int TDecGop::m_iStoredAlfParaNum = 0;
#endif
#endif

Void TDecGop::init( TDecEntropy*            pcEntropyDecoder,
                   TDecSbac*               pcSbacDecoder,
                   TDecBinCABAC*           pcBinCABAC,
                   TDecCavlc*              pcCavlcDecoder,
                   TDecSlice*              pcSliceDecoder,
                   TComLoopFilter*         pcLoopFilter,
#if ALF_HM3_REFACTOR
                   TComAdaptiveLoopFilter* pcAdaptiveLoopFilter, 
#endif
                   TComSampleAdaptiveOffset* pcSAO
#if CNN_BASED_INLOOP_FILTER
	                 , TComCNNLoopFilter*      pcCNNLoopFilter
#endif
                   )
{
  m_pcEntropyDecoder      = pcEntropyDecoder;
  m_pcSbacDecoder         = pcSbacDecoder;
  m_pcBinCABAC            = pcBinCABAC;
  m_pcCavlcDecoder        = pcCavlcDecoder;
  m_pcSliceDecoder        = pcSliceDecoder;
  m_pcLoopFilter          = pcLoopFilter;
  m_pcSAO                 = pcSAO;
  m_numberOfChecksumErrorsDetected = 0;
#if ALF_HM3_REFACTOR
  m_pcAdaptiveLoopFilter  = pcAdaptiveLoopFilter;
#endif
#if CNN_BASED_INLOOP_FILTER
	m_pcCNNLoopFilter       = pcCNNLoopFilter;
#endif
}


// ====================================================================================================================
// Private member functions
// ====================================================================================================================
// ====================================================================================================================
// Public member functions
// ====================================================================================================================
Void TDecGop::decompressSlice(TComInputBitstream* pcBitstream, TComPic* pcPic
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
                            , TComStats*  m_apcStats
#endif
  )
{
  TComSlice*  pcSlice = pcPic->getSlice(pcPic->getCurrSliceIdx());
  // Table of extracted substreams.
  // These must be deallocated AND their internal fifos, too.
  TComInputBitstream **ppcSubstreams = NULL;

  //-- For time output for each slice
  clock_t iBeforeTime = clock();
  m_pcSbacDecoder->init( (TDecBinIf*)m_pcBinCABAC );
  m_pcEntropyDecoder->setEntropyDecoder (m_pcSbacDecoder);
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  pcSlice->setStatsHandle( m_apcStats );
  m_pcEntropyDecoder->setStatsHandle ( m_apcStats );
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
    pcPic->m_aclip_prm=pcSlice->m_clip_decoded;
    g_ClipParam=pcPic->m_aclip_prm;
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
  const UInt uiNumSubstreams = pcSlice->getNumberOfSubstreamSizes()+1;

  // init each couple {EntropyDecoder, Substream}
  ppcSubstreams    = new TComInputBitstream*[uiNumSubstreams];
  for ( UInt ui = 0 ; ui < uiNumSubstreams ; ui++ )
  {
    ppcSubstreams[ui] = pcBitstream->extractSubstream(ui+1 < uiNumSubstreams ? (pcSlice->getSubstreamSize(ui)<<3) : pcBitstream->getNumBitsLeft());
  }

#if ALF_HM3_REFACTOR
  if ( pcSlice->getSPS()->getUseALF() )
  {
    m_pcAdaptiveLoopFilter->setNumCUsInFrame(pcPic);
    m_pcAdaptiveLoopFilter->allocALFParam(&m_cAlfParam);
    m_pcAdaptiveLoopFilter->resetALFParam(&m_cAlfParam);
#if COM16_C806_ALF_TEMPPRED_NUM
    static int iFirstLoop = 0;
    for( Int i = 0; i < COM16_C806_ALF_TEMPPRED_NUM && !iFirstLoop; i++ )
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
    iFirstLoop++;
#endif
  }
#endif
#if COM16_C806_ALF_TEMPPRED_NUM
  else
  {
    for( Int i = 0; i < COM16_C806_ALF_TEMPPRED_NUM; i++ )
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      for (Int j = 0; j < JVET_E0104_ALF_MAX_TEMPLAYERID; j++)
      {
#if !JVET_C0038_GALF
        m_acStoredAlfPara[j][i].coeff = NULL;
#endif
        m_acStoredAlfPara[j][i].coeff_chroma = NULL;
        m_acStoredAlfPara[j][i].coeffmulti = NULL;
        m_acStoredAlfPara[j][i].alf_cu_flag = NULL;
        m_acStoredAlfPara[j][i].alfCoeffLuma = NULL;
        m_acStoredAlfPara[j][i].alfCoeffChroma = NULL;
      }
#else
#if !JVET_C0038_GALF
      m_acStoredAlfPara[i].coeff        = NULL;
#endif
      m_acStoredAlfPara[i].coeff_chroma = NULL;
      m_acStoredAlfPara[i].coeffmulti = NULL;
      m_acStoredAlfPara[i].alf_cu_flag = NULL;
      m_acStoredAlfPara[i].alfCoeffLuma   = NULL;
      m_acStoredAlfPara[i].alfCoeffChroma = NULL;
#endif
    }
  }
#endif

#if VCEG_AZ07_INIT_PREVFRAME
  pcSlice->setStatsHandle(m_apcStats);
  pcSlice->initStatsGlobal();
#endif

  m_pcSliceDecoder->decompressSlice( ppcSubstreams, pcPic, m_pcSbacDecoder
#if ALF_HM3_REFACTOR
    , m_cAlfParam
#endif
    );

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
  m_pcEntropyDecoder->updateStates (pcSlice->getSliceType(), pcSlice->getSliceQp(), m_apcStats);
#if VCEG_AZ07_INIT_PREVFRAME_FIX
  pcSlice->updateStatsGlobal();
#endif
#endif

  // deallocate all created substreams, including internal buffers.
  for (UInt ui = 0; ui < uiNumSubstreams; ui++)
  {
    delete ppcSubstreams[ui];
  }
  delete[] ppcSubstreams;

  m_dDecTime += (Double)(clock()-iBeforeTime) / CLOCKS_PER_SEC;
}

Void TDecGop::filterPicture(TComPic* pcPic)
{
  TComSlice*  pcSlice = pcPic->getSlice(pcPic->getCurrSliceIdx());

  //-- For time output for each slice
  clock_t iBeforeTime = clock();

  // deblocking filter
  Bool bLFCrossTileBoundary = pcSlice->getPPS()->getLoopFilterAcrossTilesEnabledFlag();
  m_pcLoopFilter->setCfg(bLFCrossTileBoundary);
  m_pcLoopFilter->loopFilterPic(pcPic);
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
	printf("Total dec cnn done! Time: %.3f sec.\n", endTime);
#endif
#endif
  if (pcSlice->getSPS()->getUseSAO())
  {
    m_pcSAO->reconstructBlkSAOParams(pcPic, pcPic->getPicSym()->getSAOBlkParam());
    m_pcSAO->SAOProcess(pcPic);
    m_pcSAO->PCMLFDisableProcess(pcPic);
  }

#if RRC_NS
	if (pcSlice->getSPS()->getUseNS() && pcPic->getTLayer() < RRC_NS_TLAYER)
		filterNs(pcPic, pcSlice);
#endif

#if ALF_HM3_REFACTOR
  // adaptive loop filter
  if (pcSlice->getSPS()->getUseALF())
  {
#if COM16_C806_ALF_TEMPPRED_NUM
#if FIX_TICKET12
    if( m_pcAdaptiveLoopFilter->refreshAlfTempPred( pcSlice->getNalUnitType() , pcSlice->getPOC() ) )
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      memset(m_iStoredAlfParaNum, 0, sizeof(Int)*JVET_E0104_ALF_MAX_TEMPLAYERID);
#else
      m_iStoredAlfParaNum = 0;
#endif
      assert( m_cAlfParam.temproalPredFlag == false );
    }
#endif
    if (m_cAlfParam.temproalPredFlag)
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      m_pcAdaptiveLoopFilter->copyALFParam(&m_cAlfParam, &m_acStoredAlfPara[pcSlice->getTLayer()][m_cAlfParam.prevIdx]);
#else
      m_pcAdaptiveLoopFilter->copyALFParam(&m_cAlfParam, &m_acStoredAlfPara[m_cAlfParam.prevIdx]);
#endif
    }
#endif
    m_pcAdaptiveLoopFilter->ALFProcess(pcPic, &m_cAlfParam);
#if COM16_C806_ALF_TEMPPRED_NUM
    if (m_cAlfParam.alf_flag && !m_cAlfParam.temproalPredFlag && m_cAlfParam.filtNo >= 0)
    {
#if JVET_E0104_ALF_TEMP_SCALABILITY
      Int iCurrTempIdx = pcSlice->getTLayer();
      assert(iCurrTempIdx < JVET_E0104_ALF_MAX_TEMPLAYERID);
      for (Int iLoopedTempIdx = iCurrTempIdx; iLoopedTempIdx < JVET_E0104_ALF_MAX_TEMPLAYERID; iLoopedTempIdx++)
      {
        Int iIdx = m_iStoredAlfParaNum[iLoopedTempIdx] % COM16_C806_ALF_TEMPPRED_NUM;
        m_iStoredAlfParaNum[iLoopedTempIdx] ++;
        m_acStoredAlfPara[iLoopedTempIdx][iIdx].temproalPredFlag = false;
        m_pcAdaptiveLoopFilter->copyALFParam(&m_acStoredAlfPara[iLoopedTempIdx][iIdx], &m_cAlfParam);
#if JVET_C0038_GALF
        m_pcAdaptiveLoopFilter->resetALFPredParam(&m_acStoredAlfPara[iLoopedTempIdx][iIdx], (pcSlice->getSliceType() == I_SLICE ? true : false));
#endif
      }
#else
      Int iIdx = m_iStoredAlfParaNum % COM16_C806_ALF_TEMPPRED_NUM;
      m_iStoredAlfParaNum++;
      m_acStoredAlfPara[iIdx].temproalPredFlag = false;
      m_pcAdaptiveLoopFilter->copyALFParam(&m_acStoredAlfPara[iIdx], &m_cAlfParam);
#if JVET_C0038_GALF
      m_pcAdaptiveLoopFilter->resetALFPredParam(&m_acStoredAlfPara[iIdx], (pcSlice->getSliceType()== I_SLICE? true: false));
#endif
#endif
    }
#endif
    m_pcAdaptiveLoopFilter->freeALFParam(&m_cAlfParam);
  }
#endif

#if RRC_SELECTABLE_IF && !SIF_GAUSS_ONLY
  pcPic->getSharpeningFilter() = new TComSharpeningFilter(pcPic->getPicYuvRec(), pcSlice->getSPS());
#if !RRC_BLOCK_SHARPENING_DECODER
  pcPic->getSharpeningFilter()->buildEdgeMap(pcPic->getPicYuvRec());
#endif
#endif

#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
  if (!pcSlice->getSPS()->getAtmvpEnableFlag())
  {
    pcPic->compressMotion();
  }
#else
  pcPic->compressMotion();
#endif
  Char c = (pcSlice->isIntra() ? 'I' : pcSlice->isInterP() ? 'P' : 'B');
  if (!pcSlice->isReferenced())
  {
    c += 32;
  }

  //-- For time output for each slice
  printf("POC %4d TId: %1d ( %c-SLICE, QP%3d ) ", pcSlice->getPOC(),
    pcSlice->getTLayer(),
    c,
    pcSlice->getSliceQp());

  m_dDecTime += (Double) (clock() - iBeforeTime) / CLOCKS_PER_SEC;
  printf("[DT %6.3f] ", m_dDecTime);
  m_dDecTime = 0;

  for (Int iRefList = 0; iRefList < 2; iRefList++)
  {
    printf("[L%d ", iRefList);
    for (Int iRefIndex = 0; iRefIndex < pcSlice->getNumRefIdx(RefPicList(iRefList)); iRefIndex++)
    {
      printf("%d ", pcSlice->getRefPOC(RefPicList(iRefList), iRefIndex));
    }
    printf("] ");
  }
  if (m_decodedPictureHashSEIEnabled)
  {
    SEIMessages pictureHashes = getSeisByType(pcPic->getSEIs(), SEI::DECODED_PICTURE_HASH);
    const SEIDecodedPictureHash *hash = (pictureHashes.size() > 0) ? (SEIDecodedPictureHash*) *(pictureHashes.begin()) : NULL;
    if (pictureHashes.size() > 1)
    {
      printf("Warning: Got multiple decoded picture hash SEI messages. Using first.");
    }
    calcAndPrintHashStatus(*(pcPic->getPicYuvRec()), hash, pcSlice->getSPS()->getBitDepths(), m_numberOfChecksumErrorsDetected);
  }

  printf("\n");

  pcPic->setOutputMark(pcPic->getSlice(0)->getPicOutputFlag() ? true : false);
  pcPic->setReconMark(true);
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
  if (pcSlice->getSPS()->getUseInterKLT())
  {
#endif
      m_pcSliceDecoder->InterpolatePic(pcPic);
#if VCEG_AZ08_USE_KLT
}
#endif
#endif
}

#if RRC_NS
Void TDecGop::filterNs(TComPic * pcPic, TComSlice * pcSlice)
{

	if (pcPic->getPicSym()->getAdaptiveNs() || pcPic->getPicSym()->getNsWholeFrame())
	{
		const unsigned int channels = 1;
		TComPicYuv* rec = pcPic->getPicYuvRec();
		const unsigned int width = rec->getWidth(COMPONENT_Y);
		const unsigned int height = rec->getHeight(COMPONENT_Y);
		const unsigned int stride = rec->getStride(COMPONENT_Y);
		const int bit_depth = pcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);

		const unsigned int vec_margin = 0;
		const unsigned int vec_stride = width + 2 * vec_margin;
		const unsigned int vec_height = height + 2 * vec_margin;

		const unsigned int width_chroma = rec->getWidth(COMPONENT_Cb);
		const unsigned int height_chroma = rec->getHeight(COMPONENT_Cb);
		const unsigned int stride_chroma = rec->getStride(COMPONENT_Cb);
		const int bit_depth_chroma = pcSlice->getSPS()->getBitDepth(CHANNEL_TYPE_CHROMA);

		const unsigned int vec_margin_chroma = 0;
		const unsigned int vec_stride_chroma = width_chroma + 2 * vec_margin_chroma;
		const unsigned int vec_height_chroma = height_chroma + 2 * vec_margin_chroma;

		//    NsImage src_ns(vec_stride, vec_height);
		NsImage src_ns(vec_stride, vec_height, vec_stride_chroma, vec_height_chroma);

		float sigma = powf(2.0f, (1.0f + 0.141f * (pcSlice->getSliceQp() - 27))); // luma

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

		//    NsImage dst2_ns(src_ns);
		//    run_ns(sigma, src_ns, dst2_ns, bit_depth);

		run_ns(sigma, src_ns, src_ns, bit_depth);

		if (pcPic->getPicSym()->getAdaptiveNs())
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
			UInt *nsFlags = pcPic->getPicSym()->getNsFlag();

			const int block_width_chroma = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (ns_block_size + 1);
			const int block_height_chroma = pcPic->getSlice(0)->getSPS()->getCTUSize() >> (ns_block_size + 1);

			const int width_in_blocks_chroma = width_chroma / block_width_chroma + (width_chroma % block_width_chroma != 0);
			const int height_in_blocks_chroma = height_chroma / block_height_chroma + (height_chroma % block_height_chroma != 0);

			for (int block_index_y = 0; block_index_y < height_in_blocks; ++block_index_y)
			{
				const int block_pos_y = block_index_y * block_height;
				const int block_size_y = block_pos_y + block_height <= height ? block_height : height - block_pos_y;

				for (int block_index_x = 0; block_index_x < width_in_blocks; ++block_index_x)
				{
					const int block_pos_x = block_index_x * block_width;
					const int block_size_x = block_pos_x + block_width <= width ? block_width : width - block_pos_x;

					if (nsFlags[block_index_y * width_in_blocks + block_index_x])
					{
						Pel* src = &src_ns.m_Y[0] + vec_margin + vec_margin * vec_stride;
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

					if (nsFlags[block_index_y_chroma * width_in_blocks_chroma + block_index_x_chroma])
					{
						Pel* src = &src_ns.m_U[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
						src += block_pos_x_chroma + block_pos_y_chroma * vec_stride_chroma;
						for (int y = block_pos_y_chroma; y < block_pos_y_chroma + block_size_y_chroma; ++y, src += vec_stride_chroma)
							copy(src, src + block_size_x_chroma, &rec->getAddr(COMPONENT_Cr)[y * stride_chroma + block_pos_x_chroma]);
						src = &src_ns.m_V[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
						src += block_pos_x_chroma + block_pos_y_chroma * vec_stride_chroma;
						for (int y = block_pos_y_chroma; y < block_pos_y_chroma + block_size_y_chroma; ++y, src += vec_stride_chroma)
							copy(src, src + block_size_x_chroma, &rec->getAddr(COMPONENT_Cb)[y * stride_chroma + block_pos_x_chroma]);
					}
				}
			}
		}
		else
		{
			Pel* src = &src_ns.m_Y[0] + vec_margin + vec_margin * vec_stride;
			for (int i = 0; i < height; ++i, src += vec_stride)
				copy(src, src + width, &rec->getAddr(COMPONENT_Y)[i * stride]);

			src = &src_ns.m_U[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
			for (int i = 0; i < height_chroma; ++i, src += vec_stride_chroma)
				copy(src, src + width_chroma, &rec->getAddr(COMPONENT_Cr)[i * stride_chroma]);

			src = &src_ns.m_V[0] + vec_margin_chroma + vec_margin_chroma * vec_stride_chroma;
			for (int i = 0; i < height_chroma; ++i, src += vec_stride_chroma)
				copy(src, src + width_chroma, &rec->getAddr(COMPONENT_Cb)[i * stride_chroma]);
		}
	}//if (rpcPic->getPicSym()->getAdaptiveNs() || rpcPic->getPicSym()->getNsWholeFrame())
}
#endif

/**
 * Calculate and print hash for pic, compare to picture_digest SEI if
 * present in seis.  seis may be NULL.  Hash is printed to stdout, in
 * a manner suitable for the status line. Theformat is:
 *  [Hash_type:xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx,(yyy)]
 * Where, x..x is the hash
 *        yyy has the following meanings:
 *            OK          - calculated hash matches the SEI message
 *            ***ERROR*** - calculated hash does not match the SEI message
 *            unk         - no SEI message was available for comparison
 */
static Void calcAndPrintHashStatus(TComPicYuv& pic, const SEIDecodedPictureHash* pictureHashSEI, const BitDepths &bitDepths, UInt &numChecksumErrors)
{
  /* calculate MD5sum for entire reconstructed picture */
  TComPictureHash recon_digest;
  Int numChar=0;
  const Char* hashType = "\0";

  if (pictureHashSEI)
  {
    switch (pictureHashSEI->method)
    {
      case SEIDecodedPictureHash::MD5:
        {
          hashType = "MD5";
          numChar = calcMD5(pic, recon_digest, bitDepths);
          break;
        }
      case SEIDecodedPictureHash::CRC:
        {
          hashType = "CRC";
          numChar = calcCRC(pic, recon_digest, bitDepths);
          break;
        }
      case SEIDecodedPictureHash::CHECKSUM:
        {
          hashType = "Checksum";
          numChar = calcChecksum(pic, recon_digest, bitDepths);
          break;
        }
      default:
        {
          assert (!"unknown hash type");
          break;
        }
    }
  }

  /* compare digest against received version */
  const Char* ok = "(unk)";
  Bool mismatch = false;

  if (pictureHashSEI)
  {
    ok = "(OK)";
    if (recon_digest != pictureHashSEI->m_pictureHash)
    {
      ok = "(***ERROR***)";
      mismatch = true;
    }
  }

  printf("[%s:%s,%s] ", hashType, hashToString(recon_digest, numChar).c_str(), ok);

  if (mismatch)
  {
    numChecksumErrors++;
    printf("[rx%s:%s] ", hashType, hashToString(pictureHashSEI->m_pictureHash, numChar).c_str());
  }
}
//! \}
