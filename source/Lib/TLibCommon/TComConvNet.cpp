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


/** \file     TComConvNet.cpp
*/

#include "TComCaffeHeader.h"
#include "TComConvNet.h"
#include "TComCaffe.h"
#include "cassert"
#include <caffe/caffe.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "TComPicYuv.h"
#include "TComYuv.h"
#include "TComMv.h"
#include "TComDataCU.h"
#include "TComPic.h"
#include "TComInterpolationFilter.h"

using std::string;

TComConvNet::TComConvNet()
{
}

TComConvNet::~TComConvNet()
{

}

Void TComConvNet::create(const string& netFile, const string& modelFile, Char* query_blob_name,
#if frameworkCombine
  string platform,
#endif
  Char* mode, Int iGPUid, cv::Size inputSize
#if frameworkCombine
  , cv::Size outputSize
#endif
)
{
#if frameworkCombine
  m_dlPlatform = platform;
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {
    m_pcOpencv.create(netFile, modelFile, query_blob_name, inputSize, outputSize);
  }
  else if (strcmp(getPlatform().c_str(), "Caffe") == 0)
  {
    m_pcCaffe.create(netFile, modelFile, query_blob_name, mode, iGPUid, inputSize);
  }
  else
  {
    assert(0);
  }
#else
  m_pcCaffe.create(netFile, modelFile, query_blob_name, mode, iGPUid, inputSize);
#endif
}

Void TComConvNet::destroy(const string& netFile, const string& modelFile, Char* mode
#if frameworkCombine
  ,string platform
#endif
)
{
#if frameworkCombine
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {

  }
  else if (strcmp(getPlatform().c_str(), "Caffe") == 0)
  {
    m_pcCaffe.destroy(netFile, modelFile, mode);
  }
  else
  {
    assert(0);
  }
#else
  m_pcCaffe.destroy(netFile, modelFile, mode);
#endif
  
}
Void TComConvNet::setInputBlockSize(Int inputWidth, Int inputHeight)
{
#if frameworkCombine
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {
    m_pcOpencv.setInputSize(cv::Size(inputWidth, inputHeight));
  }
  else if (strcmp(getPlatform().c_str(), "Caffe") == 0)
  {
    m_pcCaffe.setInputSize(cv::Size(inputWidth, inputHeight));
  }
  else
  {
    assert(0);
  }
#else
  m_pcCaffe.setInputSize(cv::Size(inputWidth, inputHeight));
#endif
}
#if frameworkCombine
Void TComConvNet::setOutputBlockSize(Int outputWidth, Int outputHeight)
{
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {
    m_pcOpencv.setOutputSize(cv::Size(outputWidth, outputHeight));
  }
}
#endif
Void TComConvNet::setInputBlock(Pel* src, Int Stride, cv::Mat& inputMat, Int inputHeight, Int inputWidth)
{
#if frameworkCombine
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {
    m_pcOpencv.setInputSize(cv::Size(inputWidth, inputHeight));
  }
  else if (strcmp(getPlatform().c_str(), "Caffe") == 0)
  {
    m_pcCaffe.setInputSize(cv::Size(inputWidth, inputHeight));
  }
  else
  {
    assert(0);
  }
#else
  m_pcCaffe.setInputSize(cv::Size(inputWidth, inputHeight));
#endif
  m_iInputHeight = inputHeight;
  m_iInputWidth = inputWidth;

  for (Int y = 0; y < inputHeight; y++)
  {
    for (Int x = 0; x < inputWidth; x++)
    {
      inputMat.at<UChar>(y, x) = (UChar)src[x];
    }
    src += Stride;
  }
}

Void TComConvNet::ConvNetForward(cv::Mat& inputMat)
{
#if frameworkCombine
  if (strcmp(getPlatform().c_str(), "Opencv") == 0)
  {
    m_pcOpencv.netForward(inputMat);
    m_networkOutput = m_pcOpencv.getNetOutput();
  }
  else if (strcmp(getPlatform().c_str(), "Caffe") == 0)
  {
    m_pcCaffe.netForward(inputMat);
    m_networkOutput = m_pcCaffe.getNetOutput();
  }
  else
  {
    assert(0);
  }

#else
  m_pcCaffe.netForward(inputMat);
  m_networkOutput = m_pcCaffe.getNetOutput();
#endif
}


std::vector<Float> TComConvNet::getConvNetOutput()
{
  return m_networkOutput;
}

std::vector<Float> TComConvNet::getConvNetOutputCenter(Int iMargin)
{
  Int iheight = getInputHeight();
  Int iwidth = getInputWidth();
  Int iheightDst = iheight - 2 * iMargin;
  Int iwidthDst = iwidth - 2 * iMargin;

  std::vector<Float> netOut = m_networkOutput;
  std::vector<Float>::iterator it = netOut.begin();
  it += (iMargin*iwidth + iMargin);

  vector<Float> dst;
  for (Int y = 0; y < iheightDst; y++)
  {
    for (Int x = 0; x < iwidthDst; x++)
    {
      dst.push_back(it[y*iwidth + x]);
    }
  }
  return dst;
}
