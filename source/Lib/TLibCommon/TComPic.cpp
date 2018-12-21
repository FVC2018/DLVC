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

/** \file     TComPic.cpp
    \brief    picture class
*/

#include "TComPic.h"
#include "SEI.h"
#if RRC_SELECTABLE_IF
#include "TComSharpeningFilter.h"
#endif
#if BARC_REF_UP_SAMPLE_CNN
#include <string>
using std::string;
#endif
//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TComPic::TComPic()
: m_uiTLayer                              (0)
, m_bUsedByCurr                           (false)
, m_bIsLongTerm                           (false)
, m_pcPicYuvPred                          (NULL)
, m_pcPicYuvResi                          (NULL)
, m_bReconstructed                        (false)
, m_bNeededForOutput                      (false)
, m_uiCurrSliceIdx                        (0)
, m_bCheckLTMSB                           (false)
#if RRC_SELECTABLE_IF
, m_pSharpeningFilter                     (NULL)
#endif
{
  for(UInt i=0; i<NUM_PIC_YUV; i++)
  {
    m_apcPicYuv[i]      = NULL;
  }
#if VCEG_AZ08_INTER_KLT
  m_apcQuaPicYuv[0][0] = NULL; m_apcQuaPicYuv[0][1] = NULL; m_apcQuaPicYuv[0][2] = NULL; m_apcQuaPicYuv[0][3] = NULL;
  m_apcQuaPicYuv[1][0] = NULL; m_apcQuaPicYuv[1][1] = NULL; m_apcQuaPicYuv[1][2] = NULL; m_apcQuaPicYuv[1][3] = NULL;
  m_apcQuaPicYuv[2][0] = NULL; m_apcQuaPicYuv[2][1] = NULL; m_apcQuaPicYuv[2][2] = NULL; m_apcQuaPicYuv[2][3] = NULL;
  m_apcQuaPicYuv[3][0] = NULL; m_apcQuaPicYuv[3][1] = NULL; m_apcQuaPicYuv[3][2] = NULL; m_apcQuaPicYuv[3][3] = NULL;
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
	m_SpliceIdx = NULL;
#endif
}

TComPic::~TComPic()
{
}

#if BLOCK_BASED_BACKGROUND_REFERENCE
Void TComPic::createSpliceIdx()
{
	m_SpliceIdx = new Int[getFrameWidthInCtus()*getFrameHeightInCtus()];
	memset(m_SpliceIdx, 0, getFrameWidthInCtus()*getFrameHeightInCtus()*sizeof(Int));
}
#endif

Void TComPic::create( const TComSPS &sps, const TComPPS &pps, const Bool bIsVirtual
#if BARC_REF_UP_SAMPLE_CNN
  , TComConvNet*           cCNNUpSampleImg
#if COMPACT_RESOLUTION_CNN
  , TComConvNet*           cCNNDownSampleImg
#endif
#endif
)
{
#if BARC_REF_UP_SAMPLE_CNN
  m_cCNNUpSampleImg = cCNNUpSampleImg;
#if COMPACT_RESOLUTION_CNN
  m_cCNNDownSampleImg = cCNNDownSampleImg;
#endif
#endif
  const ChromaFormat chromaFormatIDC = sps.getChromaFormatIdc();
  const Int          iWidth          = sps.getPicWidthInLumaSamples();
  const Int          iHeight         = sps.getPicHeightInLumaSamples();
#if JVET_C0024_QTBT
  const UInt         uiMaxCuWidth    = sps.getCTUSize();
  const UInt         uiMaxCuHeight   = sps.getCTUSize();
#else
  const UInt         uiMaxCuWidth    = sps.getMaxCUWidth();
  const UInt         uiMaxCuHeight   = sps.getMaxCUHeight();
#endif
  const UInt         uiMaxDepth      = sps.getMaxTotalCUDepth();

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  m_iNumCuInWidth                    = iWidth / uiMaxCuWidth;
  m_iNumCuInWidth                   += ( iWidth % uiMaxCuWidth ) ? 1 : 0;
  m_iBaseUnitWidth                   = uiMaxCuWidth  >> uiMaxDepth;
  m_iBaseUnitHeight                  = uiMaxCuHeight >> uiMaxDepth;
#endif

  m_picSym.create( sps, pps, uiMaxDepth );
  if (!bIsVirtual)
  {
    m_apcPicYuv[PIC_YUV_ORG    ]   = new TComPicYuv;  m_apcPicYuv[PIC_YUV_ORG     ]->create( iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true );
    m_apcPicYuv[PIC_YUV_TRUE_ORG]  = new TComPicYuv;  m_apcPicYuv[PIC_YUV_TRUE_ORG]->create( iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true );
  }
  m_apcPicYuv[PIC_YUV_REC]  = new TComPicYuv;  m_apcPicYuv[PIC_YUV_REC]->create( iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true );
#if BARC
  m_apcPicYuv[PIC_YUV_Dscale]      = new TComPicYuv;  m_apcPicYuv[PIC_YUV_Dscale]->create( iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true );
  m_apcPicYuv[PIC_YUV_Dscale_REF]  = new TComPicYuv;  m_apcPicYuv[PIC_YUV_Dscale_REF]->create( iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true );
#if COMPACT_RESOLUTION_CNN
  m_apcPicYuv[PIC_YUV_CRCNN] = new TComPicYuv;  m_apcPicYuv[PIC_YUV_CRCNN]->create(iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true);
#endif
#endif
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
  if (sps.getUseInterKLT())
  {
#endif
      for (UInt uiRow = 0; uiRow < 4; uiRow++)
      {
          for (UInt uiCol = 0; uiCol < 4; uiCol++)
          {
              if (uiRow == 0 && uiCol == 0)
              {
                  m_apcQuaPicYuv[uiRow][uiCol] = m_apcPicYuv[PIC_YUV_REC];
              }
              else
              {
                  m_apcQuaPicYuv[uiRow][uiCol] = new TComPicYuv;
                  m_apcQuaPicYuv[uiRow][uiCol]->create(iWidth, iHeight, chromaFormatIDC, uiMaxCuWidth, uiMaxCuHeight, uiMaxDepth, true);
              }
          }
      }
#if VCEG_AZ08_USE_KLT
  }
#endif
#endif 

  // there are no SEI messages associated with this picture initially
  if (m_SEIs.size() > 0)
  {
    deleteSEIs (m_SEIs);
  }
  m_bUsedByCurr = false;
}

Void TComPic::destroy()
{
  m_picSym.destroy();

  for(UInt i=0; i<NUM_PIC_YUV; i++)
  {
    if (m_apcPicYuv[i])
    {
      m_apcPicYuv[i]->destroy();
      delete m_apcPicYuv[i];
      m_apcPicYuv[i]  = NULL;
    }
  }

#if VCEG_AZ08_INTER_KLT
  for (UInt uiRow = 0; uiRow < 4; uiRow++)
  {
      for (UInt uiCol = 0; uiCol < 4; uiCol++)
      {
          if (uiRow == 0 && uiCol == 0)
          {
              continue;
          }
          if (m_apcQuaPicYuv[uiRow][uiCol])
          {
              m_apcQuaPicYuv[uiRow][uiCol]->destroy();
              delete m_apcQuaPicYuv[uiRow][uiCol];
              m_apcQuaPicYuv[uiRow][uiCol] = NULL;
          }
      }
  }
#endif
#if BLOCK_BASED_BACKGROUND_REFERENCE
	if (m_SpliceIdx)
	{
		delete[] m_SpliceIdx;
		m_SpliceIdx = NULL;
}
#endif
  deleteSEIs(m_SEIs);
#if RRC_SELECTABLE_IF
  delete m_pSharpeningFilter; // TODO: check why it is not created inside this class
#endif
}
#if REFINE_BLOCK_BOUNDARY_NET
Void TComPic::upSampleDScaleRecon()
{
	Pel* piSrc = getPicYuvDScaleRef()->getAddr(COMPONENT_Y);
	Int  iSrcStride = getPicYuvDScaleRef()->getStride(COMPONENT_Y);

	Int  iDstWidth = getPicYuvRec()->getWidth(COMPONENT_Y);
	Int  iDstHeight = getPicYuvRec()->getHeight(COMPONENT_Y);

	Double* piDst = new Double[iDstWidth * iDstHeight];
	Int  iDstStride = iDstWidth;

	Int  iSrcWidth = iDstWidth / 2;
	Int  iSrcHeight = iDstHeight / 2;
	Int  iMaxBlockSize = 128;
	Int  iWidthInMaxBlock = iSrcWidth / iMaxBlockSize;
	Int  iHeightInMaxBlock = iSrcHeight / iMaxBlockSize;
#if LIYE_DEBUG
	cout << "\n" << "DstWidth: " << iDstWidth << " DstHeidht: " << iDstHeight << endl;
	cout << "\n" << "SrcWidth: " << iSrcWidth << " SrcHeidht: " << iSrcHeight << endl;
#endif
	/****Temp CTU pointer****/
	TComDataCU* pCtu = getCtu(0);
	for (Int i = 0; i < iHeightInMaxBlock; i++)
	{
		for (Int j = 0; j < iWidthInMaxBlock; j++)
		{
			Pel* srcY = piSrc + i * iMaxBlockSize * iSrcStride + j * iMaxBlockSize;
			Double* dstY = piDst + i * iMaxBlockSize * 2 * iDstStride + j * iMaxBlockSize * 2;
			CNNUpScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iMaxBlockSize, iMaxBlockSize);

		}
	}
	Int iWidthRemain = iSrcWidth - iMaxBlockSize * iWidthInMaxBlock;
	Int iHeightRemain = iSrcHeight - iMaxBlockSize * iHeightInMaxBlock;

	if (iWidthRemain > 0)
	{
		for (Int j = 0; j < iHeightInMaxBlock; j++)
		{
			Pel* srcY = piSrc + j * iMaxBlockSize * iSrcStride + iWidthInMaxBlock * iMaxBlockSize;
			Double* dstY = piDst + j * iMaxBlockSize * 2 * iDstStride + iWidthInMaxBlock * iMaxBlockSize * 2;
			CNNUpScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iWidthRemain, iMaxBlockSize);

		}
	}
	if (iHeightRemain > 0)
	{
		for (Int i = 0; i < iWidthInMaxBlock; i++)
		{
			Pel* srcY = piSrc + iHeightInMaxBlock * iMaxBlockSize * iSrcStride + i * iMaxBlockSize;
			Double* dstY = piDst + iHeightInMaxBlock * iMaxBlockSize * 2 * iDstStride + i * iMaxBlockSize * 2;
			CNNUpScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iMaxBlockSize, iHeightRemain);
		}
	}
	if (iHeightRemain > 0 && iWidthRemain > 0)
	{
		Pel* srcY = piSrc + iHeightInMaxBlock * iMaxBlockSize * iSrcStride + iWidthInMaxBlock * iMaxBlockSize;
		Double* dstY = piDst + iHeightInMaxBlock * iMaxBlockSize * 2 * iDstStride + iWidthInMaxBlock * iMaxBlockSize * 2;
		CNNUpScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iWidthRemain, iHeightRemain);
	}
	//DCTIF up-sampling
	Pel* piDstDctif = new Pel[iDstWidth * iDstHeight];
	getPicSym()->getCtu(0)->upScaleComponent2x2(piSrc, piDstDctif, iSrcStride, iDstStride, iSrcWidth, iSrcHeight, false, false, COMPONENT_Y);
	// copy up-sampled recon to real recon where the mode is barc && the up-sampling method is learning

	iSrcStride = iDstWidth;
	iDstStride = getPicYuvRec()->getStride(COMPONENT_Y);

	Int iWidthInCU = getPicSym()->getFrameWidthInCtus();
	Int iHeightInCU = getPicSym()->getFrameHeightInCtus();
	Int iTotalCU = getPicSym()->getNumberOfCtusInFrame();
	for (Int iCU = 0; iCU < iTotalCU; iCU++)
	{
		TComDataCU* pcCU = getPicSym()->getCtu(iCU);
		if (!pcCU->getBARCFlag(CHANNEL_TYPE_LUMA))
		{
			continue;
		}
		Int iCUPelXLeft = pcCU->getCUPelX();
		Int iCUPelYTop = pcCU->getCUPelY();
		Int iCUPelXRight = iCUPelXLeft + pcCU->getSlice()->getSPS()->getCTUSize() - 1;
		Int iCUPelYBottom = iCUPelYTop + pcCU->getSlice()->getSPS()->getCTUSize() - 1;
		if (iCUPelXRight > iDstWidth - 1)
		{
			iCUPelXRight = iDstWidth - 1;
		}
		if (iCUPelYBottom > iDstHeight - 1)
		{
			iCUPelYBottom = iDstHeight - 1;
		}
		Double* pSrc1 = piDst + iCUPelYTop * iSrcStride;
		Pel* pSrc2 = piDstDctif + iCUPelYTop * iSrcStride;
		Pel* piDstRecon = getPicYuvRec()->getAddr(COMPONENT_Y) + iCUPelYTop * iDstStride;
#if BARC_UP_SAMPLE_RDO
		if (pcCU->getUpCNNFlag())
#else
		if (1)
#endif
		{
			for (Int i = iCUPelYTop; i <= iCUPelYBottom; i++)
			{
				for (Int j = iCUPelXLeft; j <= iCUPelXRight; j++)
				{
#if COMPACT_RESOLUTION_CNN

					piDstRecon[j] = ClipA(int(pSrc1[j] * PixelRangeInt + 0.5), COMPONENT_Y);
#else
					piDstRecon[j] = ClipA(int(pSrc2[j] + pSrc1[j] * PixelRangeInt + 0.5), COMPONENT_Y);
#endif

				}
				pSrc1 += iSrcStride;
				pSrc2 += iSrcStride;
				piDstRecon += iDstStride;
			}
		}
		else
		{
			for (Int i = iCUPelYTop; i <= iCUPelYBottom; i++)
			{
				for (Int j = iCUPelXLeft; j <= iCUPelXRight; j++)
				{
					piDstRecon[j] = pSrc2[j];
				}
				pSrc1 += iSrcStride;
				pSrc2 += iSrcStride;
				piDstRecon += iDstStride;
			}
		}
	}
	if (piDst)
	{
		delete[] piDst;
		piDst = NULL;
	}
	if (piDstDctif)
	{
		delete[] piDstDctif;
		piDstDctif = NULL;
	}
}
#endif

#if BARC_REF_UP_SAMPLE_CNN
Void TComPic::initCNNModel()
{
  //init CNN Up-sample model
#if COMPACT_RESOLUTION_CNN
  string caffeModelFile = DEFAULT_CNN_UP_SAMPLE_MODEL_FILE_CRCNN;
  string caffeNetFile = DEFAULT_CNN_UP_SAMPLE_NET_FILE_CRCNN;
  char*  queryBlobName = "conv18";
#else
  string caffeModelFile = DEFAULT_CNN_UP_SAMPLE_MODEL_FILE_BICUBIC;
  string caffeNetFile = DEFAULT_CNN_UP_SAMPLE_NET_FILE_BICUBIC;
  char*  queryBlobName = "conv4";
#endif
#if frameworkCombine
  string platform = getSlice(getCurrSliceIdx())->getSPS()->getDLframework();
  if (platform == "Tensorflow")
  {
    platform = DEFAULT_DL_FRAMEWORK;
  }
#endif
  Int iGPUid = getSlice(getCurrSliceIdx())->getSPS()->getGPUid();
  //cout << iGPUid;
  //char*  mode = DEFAULT_WORKING_MODE;
  string WorkingMode = getSlice(getCurrSliceIdx())->getSPS()->getWorkingMode();
  Char* mode = new Char[WorkingMode.length() + 1];
  WorkingMode.copy(mode, WorkingMode.length(), 0);
  mode[WorkingMode.length()] = '\0';
  Int iQP = getSlice(getCurrSliceIdx())->getSliceQp();
  Int qp[] = { 27, 32, 37, 42 };
  Int iDiff, iBestQP;
  Int iMin = 100;
  for (Int i = 0; i < 4; i++)
  {
    iDiff = abs(iQP - qp[i]);
    if (iMin > iDiff)
    {
      iBestQP = qp[i];
      iMin = iDiff;
    }
  }
  stringstream ss;
  string iTmpQP;
  ss << iBestQP;
  ss >> iTmpQP;
  caffeModelFile.replace(caffeModelFile.find("?"), 2, iTmpQP);
#if frameworkCombine 
  m_cCNNUpSampleImg->create(caffeNetFile, caffeModelFile, queryBlobName, platform, mode, iGPUid, cv::Size(0, 0), cv::Size(0, 0));
#else
  m_cCNNUpSampleImg->create(caffeNetFile, caffeModelFile, queryBlobName, mode, iGPUid, cv::Size(0, 0));
#endif
#if COMPACT_RESOLUTION_CNN
  // init CNN Down Sample Model
  caffeNetFile = DEFAULT_CNN_DOWN_SAMPLE_NET_FILE_CRCNN;
  caffeModelFile = DEFAULT_CNN_DOWN_SAMPLE_MODEL_FILE_CRCNN;
#if frameworkCombine
  if (strcmp(platform.c_str(), "Opencv") == 0)
  {
    char* query_blob_name_crcnn = "sum_ds_global";
    m_cCNNDownSampleImg->create(caffeNetFile, caffeModelFile, query_blob_name_crcnn, platform, mode, iGPUid, cv::Size(0, 0), cv::Size(0, 0));
  }
  else if (strcmp(platform.c_str(), "Caffe") == 0)
  {
    char* query_blob_name_crcnn = "ds_conv10_global";
    m_cCNNDownSampleImg->create(caffeNetFile, caffeModelFile, query_blob_name_crcnn, platform, mode, iGPUid, cv::Size(0, 0), cv::Size(0, 0));
  }
  else
  {
    assert(0);
  }
#else
  char* query_blob_name_crcnn = "ds_conv10_global";
  m_cCNNDownSampleImg->create(caffeNetFile, caffeModelFile, query_blob_name_crcnn, mode, iGPUid, cv::Size(0, 0));
#endif

#endif
  }

Void TComPic::CNNUpScaleLumaComponent2x2(Bool uiFlag, const Pel* pSrc, Double* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight)
{
  Int M = iSrcHeight;//行 
  Int N = iSrcWidth;//列
  Int istrideREC = iSrcStride;
  const Pel* srcYUV = pSrc;
#if BitDepth10
  cv::Mat input_img(M, N, CV_16UC1);
#else
  cv::Mat input_img(M, N, CV_8UC1);
#endif

  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
#if BitDepth10
      input_img.at<short>(i, j) = (short)srcYUV[j];
#else
      input_img.at<uchar>(i, j) = (uchar)srcYUV[j];
#endif
    }
    srcYUV += istrideREC;
  }
  CHECK(!input_img.empty()) << "Unable to reconstruct image ";
  getConvNetUpSample()->setInputBlockSize(input_img.cols, input_img.rows);
#if frameworkCombine
  getConvNetUpSample()->setOutputBlockSize(2 * N, 2 * M);
#endif
  getConvNetUpSample()->ConvNetForward(input_img);

  UInt frameworkMode = 0;//0=Caffe,1=Opencv,2=Tensorflow
  std::vector<float> result = getConvNetUpSample()->getConvNetOutput();
  istrideREC = iDestStride;
  for (int i = 0; i < 2 * M; i++)
  {
    for (int j = 0; j < 2 * N; j++)
    {
#if COMPACT_RESOLUTION_CNN
      pDest[j] = result[i * 2 * N + j] + 0.5;
#else
      pDest[j] = result[i * 2 * N + j];
#endif
    }
    pDest += istrideREC;
  }
}


#if COMPACT_RESOLUTION_CNN
Void TComPic::downScaleOrgCRCNN()
{
  Pel* piSrc = getPicYuvOrg()->getAddr(COMPONENT_Y);
  Int  iSrcStride = getPicYuvOrg()->getStride(COMPONENT_Y);

  Int  iDstWidth = getPicYuvOrg()->getWidth(COMPONENT_Y) / 2;
  Int  iDstHeight = getPicYuvOrg()->getHeight(COMPONENT_Y) / 2;

  Pel* piDst = getPicYuvCRCNNOrg()->getAddr(COMPONENT_Y);
  Int  iDstStride = iSrcStride;

  Int  iSrcWidth = iDstWidth * 2;
  Int  iSrcHeight = iDstHeight * 2;
  Int  iMaxBlockSize = 512;

  Int  iWidthInMaxBlock = iSrcWidth / iMaxBlockSize;
  Int  iHeightInMaxBlock = iSrcHeight / iMaxBlockSize;
  for (Int i = 0; i < iHeightInMaxBlock; i++)
  {
    for (Int j = 0; j < iWidthInMaxBlock; j++)
    {
      Pel* srcY = piSrc + i * iMaxBlockSize * iSrcStride + j * iMaxBlockSize;
      Pel* dstY = piDst + i * iMaxBlockSize / 2 * iDstStride + j * iMaxBlockSize / 2;
      CNNDownScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iMaxBlockSize, iMaxBlockSize);
    }
  }

  Int iWidthRemain = iSrcWidth - iMaxBlockSize * iWidthInMaxBlock;
  Int iHeightRemain = iSrcHeight - iMaxBlockSize * iHeightInMaxBlock;

  if (iWidthRemain > 0)
  {
    for (Int j = 0; j < iHeightInMaxBlock; j++)
    {
      Pel* srcY = piSrc + j * iMaxBlockSize * iSrcStride + iWidthInMaxBlock * iMaxBlockSize;
      Pel* dstY = piDst + j * iMaxBlockSize / 2 * iDstStride + iWidthInMaxBlock * iMaxBlockSize / 2;
      CNNDownScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iWidthRemain, iMaxBlockSize);
    }
  }
  if (iHeightRemain > 0)
  {
    for (Int i = 0; i < iWidthInMaxBlock; i++)
    {
      Pel* srcY = piSrc + iHeightInMaxBlock * iMaxBlockSize * iSrcStride + i * iMaxBlockSize;
      Pel* dstY = piDst + iHeightInMaxBlock * iMaxBlockSize / 2 * iDstStride + i * iMaxBlockSize / 2;
      CNNDownScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iMaxBlockSize, iHeightRemain);
    }
  }
  if (iHeightRemain > 0 && iWidthRemain > 0)
  {
    Pel* srcY = piSrc + iHeightInMaxBlock * iMaxBlockSize * iSrcStride + iWidthInMaxBlock * iMaxBlockSize;
    Pel* dstY = piDst + iHeightInMaxBlock * iMaxBlockSize / 2 * iDstStride + iWidthInMaxBlock * iMaxBlockSize / 2;
    CNNDownScaleLumaComponent2x2(false, srcY, dstY, iSrcStride, iDstStride, iWidthRemain, iHeightRemain);
  }
}
Void  TComPic::CNNDownScaleLumaComponent2x2(Bool uiFlag, const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight)
{

  Int M = iSrcHeight;//行 
  Int N = iSrcWidth;//列
  Int istrideREC = iSrcStride;
  const Pel* srcYUV = pSrc;
  cv::Mat input_img(M, N, CV_16UC1);

  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      input_img.at<ushort>(i, j) = (ushort)srcYUV[j];
    }
    srcYUV += istrideREC;
  }

  CHECK(!input_img.empty()) << "Unable to reconstruct image ";
  getConvNetDownSample()->setInputBlockSize(input_img.cols, input_img.rows);
#if frameworkCombine
  getConvNetDownSample()->setOutputBlockSize(N / 2, M / 2);
#endif
  getConvNetDownSample()->ConvNetForward(input_img);
  std::vector<float> result = getConvNetDownSample()->getConvNetOutput();
  //std::vector<float> result = m_pcPic->getConvNet()->getConvNetOutput();
  istrideREC = iDestStride;

  //std::vector<float> result_img = m_netCompactImgLuma.recon(input_img, M, N);
  /*for (int i = 0; i < result_img.size(); i++)
  {
  std::cout << result_img[i] << " ";
  }*/
  //	float tmp_val;
  //istrideREC = iDestStride;
  for (int i = 0; i < M / 2; i++)
  {
    for (int j = 0; j < N / 2; j++)
    {
      pDest[j] = ClipA(int(result[i * N / 2 + j] * PixelRangeDouble + 0.5), COMPONENT_Y);
    }
    pDest += istrideREC;
  }

}
#endif
#endif
Void TComPic::compressMotion()
{
  TComPicSym* pPicSym = getPicSym();
  for ( UInt uiCUAddr = 0; uiCUAddr < pPicSym->getNumberOfCtusInFrame(); uiCUAddr++ )
  {
    TComDataCU* pCtu = pPicSym->getCtu(uiCUAddr);
    pCtu->compressMV();
  }
}

Bool  TComPic::getSAOMergeAvailability(Int currAddr, Int mergeAddr)
{
  Bool mergeCtbInSliceSeg = (mergeAddr >= getPicSym()->getCtuTsToRsAddrMap(getCtu(currAddr)->getSlice()->getSliceCurStartCtuTsAddr()));
  Bool mergeCtbInTile     = (getPicSym()->getTileIdxMap(mergeAddr) == getPicSym()->getTileIdxMap(currAddr));
  return (mergeCtbInSliceSeg && mergeCtbInTile);
}

UInt TComPic::getSubstreamForCtuAddr(const UInt ctuAddr, const Bool bAddressInRaster, TComSlice *pcSlice)
{
  UInt subStrm;
  const bool bWPPEnabled=pcSlice->getPPS()->getEntropyCodingSyncEnabledFlag();
  const TComPicSym &picSym            = *(getPicSym());

  if ((bWPPEnabled && picSym.getFrameHeightInCtus()>1) || (picSym.getNumTiles()>1)) // wavefronts, and possibly tiles being used.
  {
    if (bWPPEnabled)
    {
      const UInt ctuRsAddr                = bAddressInRaster?ctuAddr : picSym.getCtuTsToRsAddrMap(ctuAddr);
      const UInt frameWidthInCtus         = picSym.getFrameWidthInCtus();
      const UInt tileIndex                = picSym.getTileIdxMap(ctuRsAddr);
      const UInt numTileColumns           = (picSym.getNumTileColumnsMinus1()+1);
      const TComTile *pTile               = picSym.getTComTile(tileIndex);
      const UInt firstCtuRsAddrOfTile     = pTile->getFirstCtuRsAddr();
      const UInt tileYInCtus              = firstCtuRsAddrOfTile / frameWidthInCtus;
      // independent tiles => substreams are "per tile"
      const UInt ctuLine                  = ctuRsAddr / frameWidthInCtus;
      const UInt startingSubstreamForTile =(tileYInCtus*numTileColumns) + (pTile->getTileHeightInCtus()*(tileIndex%numTileColumns));
      subStrm = startingSubstreamForTile + (ctuLine - tileYInCtus);
    }
    else
    {
      const UInt ctuRsAddr                = bAddressInRaster?ctuAddr : picSym.getCtuTsToRsAddrMap(ctuAddr);
      const UInt tileIndex                = picSym.getTileIdxMap(ctuRsAddr);
      subStrm=tileIndex;
    }
  }
  else
  {
    // dependent tiles => substreams are "per frame".
    subStrm = 0;
  }
  return subStrm;
}
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
Void TComPic::getCUAddrAndPartIdx( Int iX, Int iY, Int& riCuAddr, Int& riAbsZorderIdx )
{
#if JVET_C0024_QTBT
  Int iMaxCUWidth   = (Int) ( getPicSym()->getSPS().getCTUSize()  );
  Int iMaxCuHeight  = (Int) ( getPicSym()->getSPS().getCTUSize() );
#else
  Int iMaxCUWidth   = (Int) ( getPicSym()->getSPS().getMaxCUWidth()  );
  Int iMaxCuHeight  = (Int) ( getPicSym()->getSPS().getMaxCUHeight() );
#endif

  Int iCuX            = iX / iMaxCUWidth;
  Int iCuY            = iY / iMaxCuHeight;
  Int iBaseX          = ( iX - iCuX * iMaxCUWidth  ) / m_iBaseUnitWidth;
  Int iBaseY          = ( iY - iCuY * iMaxCuHeight ) / m_iBaseUnitHeight;
  Int iCuSizeInBases  = iMaxCUWidth                  / m_iBaseUnitWidth;
  riCuAddr            = iCuY   * m_iNumCuInWidth + iCuX;
  Int iRastPartIdx    = iBaseY * iCuSizeInBases  + iBaseX;
  riAbsZorderIdx      = g_auiRasterToZscan[ iRastPartIdx ];
}
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
namespace {
Bound computeBoundsComp(const Pel *p,Int height,Int width,Int stride)
{
    Bound b;
    b.m=b.M=*p;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
        {
            const int x=p[i*stride + j];
            if (x> b.M) b.M = x;
            if (x< b.m) b.m = x;
        }
    return b;
}

int count(const Pel *p,Int height,Int width,Int stride,Int m,Int M) {
    Int s=0;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
        {
            const int x=p[i*stride + j];
            if (x>=m&&x<=M) ++s;
        }
    return s;
}
}

ClipParam TComPic::computeTchClipParam(Int &delta_disto_luma,Int &delta_disto_chroma) const
{
    ClipParam prm;
    prm.isActive=true;
    prm.isChromaActive=true;
    delta_disto_luma=delta_disto_chroma=0;

    const TComPicYuv &picorg = *m_apcPicYuv[PIC_YUV_ORG];

    const Pel *pY = picorg.getAddr(COMPONENT_Y);
    Int strideY   = picorg.getStride(COMPONENT_Y);
    Int widthY    = picorg.getWidth(COMPONENT_Y);
    Int heightY   = picorg.getHeight(COMPONENT_Y);

    Bound Y=computeBoundsComp(pY,heightY,widthY,strideY);
    prm.Y().m = Y.m;
    prm.Y().M = Y.M;
    const int kMargin=8; // margin (in pixel value) to consider near the bound
    {
        if (Y.m>0)                            delta_disto_luma+=count(pY,heightY,widthY,strideY,Y.m        ,Y.m+kMargin);
        if (Y.M<(1<<ClipParam::ibdLuma)-1) delta_disto_luma+=count(pY,heightY,widthY,strideY,Y.M-kMargin,Y.M        );
    }
    const Pel *pCb = picorg.getAddr(COMPONENT_Cb);
    Int strideC    = picorg.getStride(COMPONENT_Cb);
    Int widthC     = picorg.getWidth(COMPONENT_Cb);
    Int heightC    = picorg.getHeight(COMPONENT_Cb);

    Bound U=computeBoundsComp(pCb,heightC,widthC,strideC);

    const Pel *pCr = picorg.getAddr(COMPONENT_Cr);
    Bound V=computeBoundsComp(pCr,heightC,widthC,strideC);
    {
        if (U.m>0)                              delta_disto_chroma+=count(pCb,heightC,widthC,strideC,U.m        ,U.m+kMargin);
        if (U.M<(1<<ClipParam::ibdChroma)-1) delta_disto_chroma+=count(pCb,heightC,widthC,strideC,U.M-kMargin,U.M        );
        if (V.m>0)                              delta_disto_chroma+=count(pCr,heightC,widthC,strideC,V.m        ,V.m+kMargin);
        if (V.M<(1<<ClipParam::ibdChroma)-1) delta_disto_chroma+=count(pCr,heightC,widthC,strideC,V.M-kMargin,V.M        );
    }
    prm.U()=U;
    prm.V()=V;

    return prm;
}
#endif
#if JVET_C0024_QTBT
Void TComPic::setCodedBlkInCTU(Bool bCoded, UInt uiBlkX, UInt uiBlkY, UInt uiWidth, UInt uiHeight)
{
  assert(sizeof(**m_bCodedBlkInCTU)==1);
  for (UInt i=uiBlkY; i<uiBlkY+uiHeight; i++)
  {
    memset(&m_bCodedBlkInCTU[i][uiBlkX], bCoded, uiWidth);
  }
}
Int TComPic::getCodedAreaInCTU()
{
  return m_iCodedArea;
}

Void TComPic::setCodedAreaInCTU(Int iArea)
{
  m_iCodedArea = iArea;
}

Void TComPic::addCodedAreaInCTU(Int iArea)
{
  m_iCodedArea += iArea;
  assert(m_iCodedArea>=0);
}

Void  TComPic::setSkiped(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bSkiped)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  m_bSkiped[uiZorder][uiWIdx][uiHIdx] = bSkiped;  
}
Bool  TComPic::getSkiped(UInt uiZorder, UInt uiWidth, UInt uiHeight)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  return m_bSkiped[uiZorder][uiWIdx][uiHIdx];
}
Void  TComPic::clearAllSkiped()
{
  memset(m_bSkiped, 0, (1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1))*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*sizeof(Bool));
}

Void  TComPic::setInter(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bInter)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  m_bInter[uiZorder][uiWIdx][uiHIdx] = bInter; 
}
Bool  TComPic::getInter(UInt uiZorder, UInt uiWidth, UInt uiHeight)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  return m_bInter[uiZorder][uiWIdx][uiHIdx];
}
Void  TComPic::clearAllInter()
{
  memset(m_bInter, 0, (1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1))*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*sizeof(Bool));
}

Void  TComPic::setIntra(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bIntra)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  m_bIntra[uiZorder][uiWIdx][uiHIdx] = bIntra; 
}
Bool  TComPic::getIntra(UInt uiZorder, UInt uiWidth, UInt uiHeight)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  return m_bIntra[uiZorder][uiWIdx][uiHIdx];
}
Void  TComPic::clearAllIntra()
{
  memset(m_bIntra, 0, (1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1))*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*sizeof(Bool));
}

Void  TComPic::setIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx, TComMv cMv)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  m_cIntMv[uiZorder][uiWIdx][uiHIdx][(UInt)eRefList][uiRefIdx] = cMv;
  m_bSetIntMv[uiZorder][uiWIdx][uiHIdx][(UInt)eRefList][uiRefIdx] = true; 
}

TComMv  TComPic::getIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  return m_cIntMv[uiZorder][uiWIdx][uiHIdx][(UInt)eRefList][uiRefIdx];
}

Bool  TComPic::IsSetIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx)
{
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  return m_bSetIntMv[uiZorder][uiWIdx][uiHIdx][(UInt)eRefList][uiRefIdx];
}

Void  TComPic::clearAllIntMv()
{
  memset(m_bSetIntMv, 0, (1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1))*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*(MAX_CU_DEPTH-MIN_CU_LOG2+1)*2*5*sizeof(Bool));
}
#endif
//! \}
