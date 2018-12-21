/* The copyright in this software is being made available under the BSD
* License, included below. This software may be subject to other third party
* and contributor rights, including patent rights, and no such rights are
* granted under this license.
*
* Copyright (c) 2010-2013, ITU/ISO/IEC
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


/** \file     TComSplitImg.cpp
*/
#include "TComSplitImg.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include "TComConvNet.h"
#include <string>

#if USE_TENSORFLOW
#include "TensorFlowHeader.h"
#include "TComTF.h"
#endif

#if CNN_BASED_INLOOP_FILTER

TComSplitImg::TComSplitImg(Int maxWidth, Int maxHeight)
{
  m_maxBlockWidth = maxWidth;
  m_maxBlockHeight = maxHeight;
  m_blockNumInVer = 1;
  m_blockNumInHor = 1;
}
TComSplitImg::~TComSplitImg()
{
}
Void TComSplitImg::picYuv2Mat(TComPic* pcPic, Int compIdx, cv::Mat &inputImg, TComPicYuv* pcPicYuvRec)
{
  Int istride = pcPicYuvRec->getStride(ComponentID(compIdx));
  Int imgHeight = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicHeightInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
  Int imgWidth = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicWidthInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
  Pel* srcYUV = pcPicYuvRec->getAddr(ComponentID(compIdx), 0, 0);

  for (Int i = 0; i < imgHeight; i++)
  {
    for (Int j = 0; j < imgWidth; j++)
    {
      inputImg.at<ushort>(i, j) = srcYUV[j];
    }
    srcYUV += istride;
  }

}
Void TComSplitImg::mat2PicYuv(TComPic* pcPic, Int compIdx, cv::Mat &recImg, TComPicYuv* pcPicYuvRec)
{
  Int istride = pcPicYuvRec->getStride(ComponentID(compIdx));
  Int imgHeight = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicHeightInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
  Int imgWidth = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicWidthInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
  Pel* srcYUV = pcPicYuvRec->getAddr(ComponentID(compIdx), 0, 0);

  for (Int i = 0; i < imgHeight; i++)
  {
    for (Int j = 0; j < imgWidth; j++)
    {
      srcYUV[j] = Clip3(0, 1023, (Int)recImg.at<ushort>(i, j));
    }
    srcYUV += istride;
  }
}

Void TComSplitImg::splitImg(TComPic* pcPic, std::vector<BlockInfo>&subImg, Int compIdx, TComPicYuv* pcPicYuvRec)
{
  Int imgHeight = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicHeightInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
  Int imgWidth = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicWidthInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));

  cv::Mat inputImg;
  inputImg.create(imgHeight, imgWidth, CV_16UC1);

  picYuv2Mat(pcPic, compIdx, inputImg, pcPicYuvRec);

  Bool isSplit = false;

  m_blockNumInVer = 1;
  m_blockNumInHor = 1;
  if (imgWidth>m_maxBlockWidth || imgHeight > m_maxBlockHeight)
  {
    isSplit = true;
    m_blockNumInVer = (Int)ceil(imgHeight / float(m_maxBlockHeight));
    m_blockNumInHor = (Int)ceil(imgWidth / float(m_maxBlockWidth));
  }

  BlockInfo tempImg;

  for (Int row = 0; row < m_blockNumInVer; row++)
    for (Int col = 0; col < m_blockNumInHor; col++)
    {
			if (isSplit)
			{
				tempImg.blockX = col*(m_maxBlockWidth);
				tempImg.blockY = row*(m_maxBlockHeight);
				tempImg.blockHeight = (tempImg.blockY + m_maxBlockHeight>imgHeight) ? (imgHeight - tempImg.blockY) : m_maxBlockHeight;
				tempImg.blockWidth = (tempImg.blockX + m_maxBlockWidth>imgWidth) ? (imgWidth - tempImg.blockX) : m_maxBlockWidth;

				tempImg.block = inputImg(cv::Rect(col*(m_maxBlockWidth), row*(m_maxBlockHeight), tempImg.blockWidth, tempImg.blockHeight));
			}
			else
			{
				tempImg.block = inputImg(cv::Rect(0, 0, imgWidth, imgHeight));
				tempImg.blockX = 0;
				tempImg.blockY = 0;
				tempImg.blockHeight = imgHeight;
				tempImg.blockWidth = imgWidth;
			}
			subImg.push_back(tempImg);
    }
}
Void TComSplitImg::blockCombine(TComPic* pcPic, cv::Mat &recImg, BlockInfo subImg, std::vector<float>result)
{
  for (Int i = 0; i < subImg.blockHeight; i++)
  {
    for (Int j = 0; j < subImg.blockWidth; j++)
    {
      recImg.at<ushort>(subImg.blockY + i, subImg.blockX + j) = Clip3(0, 1023, (Int)((result[i*subImg.blockWidth + j]) * 1023));
    }
  }
  result.erase(result.begin(), result.end());
}

Void TComSplitImg::postProcess(TComPic* pcPic, TComPicYuv *pcPicYuvRec
#if frameworkCombine
  , string platform
#endif
)
{ 
	if (pcPicYuvRec == NULL)
		pcPicYuvRec = pcPic->getPicYuvRec();
	Int iQP = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSliceQp();
#if USE_NEW_MODEL
	Int qp[] = { 28,34,40,46 };
#else
	Int qp[] = {22,27,32,37 };
#endif
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

  //double endTime;
  //clock_t startClock = clock();
  //printf("caffe cnn loop filter start.\n");

  string caffeModelPath = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getCaffeModelPath();
  string caffeModelFile = caffeModelPath + DEFAULT_LOOP_FILTER_CAFFEMODEL_FILENAME_TEMPLATE;
  string caffeNetFile = caffeModelPath + DEFAULT_LOOP_FILTER_NET_PROTO_FILENAME;
  Char *queryBlobName = DEFAULT_CAFFEMODEL_OUTPUT_LAYER_NAME;
  string WorkingMode = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getWorkingMode();
  Char* mode = new Char[WorkingMode.length() + 1];
  WorkingMode.copy(mode, WorkingMode.length(), 0);
  mode[WorkingMode.length()] = '\0';
  Int GPUid = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getGPUid();
  caffeModelFile.replace(caffeModelFile.find("?"), 2, iTmpQP);

  TComConvNet CNNReconImg;
#if frameworkCombine
  CNNReconImg.create(caffeNetFile, caffeModelFile, queryBlobName, platform, mode, GPUid, cv::Size(0, 0), cv::Size(0, 0));
#else
  CNNReconImg.create(caffeNetFile, caffeModelFile, queryBlobName, mode, GPUid, cv::Size(0, 0));
#endif

  for (Int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    Int imgHeight = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicHeightInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
    Int imgWidth = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicWidthInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));

    cv::Mat recImg(imgHeight, imgWidth, CV_16UC1);
    vector<BlockInfo>subImg;

    splitImg(pcPic, subImg, compIdx, pcPicYuvRec);

    for (Int i = 0; i < subImg.size(); i++)
    {
      //printf("block %d,height %d, width %d\n", i, subImg[i].block.cols, subImg[i].block.rows);
      CHECK(!subImg[i].block.empty()) << "Unable to reconstruct image ";
      
      CNNReconImg.setInputBlockSize(subImg[i].block.cols, subImg[i].block.rows);
#if frameworkCombine
      CNNReconImg.setOutputBlockSize(subImg[i].block.cols, subImg[i].block.rows);
#endif
      //printf("CNN forword begin\n");
      //system("pause");
      CNNReconImg.ConvNetForward(subImg[i].block);
      //printf("CNN forword end\n");
      //system("pause");
      std::vector<float> result = CNNReconImg.getConvNetOutput();
      blockCombine(pcPic, recImg, subImg[i], result);
    }
    subImg.erase(subImg.begin(), subImg.end());
    mat2PicYuv(pcPic, compIdx, recImg, pcPicYuvRec);
  }
  if (platform == "Caffe")
  {
    CNNReconImg.destroy(caffeNetFile, caffeModelFile, mode, platform);
  }
  //endTime = (double)(clock() - startClock) / CLOCKS_PER_SEC;
  //printf("caffe cnn loop filter done! Time: %.3f sec.\n", endTime);
}
#if USE_TENSORFLOW
Void TComSplitImg::postProcess_tf(TComPic* pcPic, TComPicYuv *pcPicYuvRec)
{
  if (pcPicYuvRec == NULL)
    pcPicYuvRec = pcPic->getPicYuvRec();
  Int iQP = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSliceQp();
#if USE_NEW_MODEL
  Int qp[] = { 28,34,40,46 };
#else
  Int qp[] = { 22,27,32,37 };
#endif
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

  //double endTime;
  //clock_t startClock = clock();
  //printf("tensorflow cnn loop filter start.\n");

  string graph_path = "./TensorflowModel/frozen_model.pb";
  Int GPUid = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getGPUid();
  stringstream ss1;
  string sGPUid;
  ss1 << GPUid;
  ss1 >> sGPUid;
  Char* cGPUid = new Char[sGPUid.length() + 1];
  sGPUid.copy(cGPUid, sGPUid.length(), 0);
  cGPUid[sGPUid.length()] = '\0';

  for (Int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    Int imgHeight = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicHeightInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
    Int imgWidth = pcPic->getSlice(pcPic->getCurrSliceIdx())->getSPS()->getPicWidthInLumaSamples() >> pcPicYuvRec->getComponentScaleX(ComponentID(compIdx));
    Int istride = pcPicYuvRec->getStride(ComponentID(compIdx));
    Pel* srcYUV = pcPicYuvRec->getAddr(ComponentID(compIdx), 0, 0);

    cv::Mat recImg(imgHeight, imgWidth, CV_16UC1);
    vector<BlockInfo>subImg;

    splitImg(pcPic, subImg, compIdx, pcPicYuvRec);

    for (Int i = 0; i < subImg.size(); i++)
    {
      //printf("block %d,height %d, width %d\n", i, subImg[i].block.cols, subImg[i].block.rows);
      CHECK(!subImg[i].block.empty()) << "Unable to reconstruct image ";

      Pel* piSubImg= srcYUV + subImg[i].blockY*istride + subImg[i].blockX;

      TComTF TFTest;
      TFTest.TFNetForward(graph_path, cGPUid, piSubImg, subImg[i].block.rows, subImg[i].block.cols, istride, "images", "output/pred/add");
    }
    subImg.erase(subImg.begin(), subImg.end());
  }
  //endTime = (double)(clock() - startClock) / CLOCKS_PER_SEC;
  //printf("tensorflow cnn loop filter done! Time: %.3f sec.\n", endTime);
}
#endif
#endif