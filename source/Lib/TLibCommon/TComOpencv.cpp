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

/** \file     TComOpencv.cpp
\brief    Caffe interface c lass
*/

#include "TComOpencv.h"

#if CNN_BASED_FRAC_INTERPOLATION
//#include "TComInterpolationFilter.h"
#include "CommonDef.h"
#endif

#if frameworkCombine

TComOpencv::TComOpencv()
{

}

TComOpencv::~TComOpencv()
{

}

static bool init_glog = false;

Void TComOpencv::create(const string& netFile, const string& modelFile, Char* query_blob_name, cv::Size inputSize, cv::Size outputSize)
{
  m_blobName = query_blob_name;

  m_net = cv::dnn::readNetFromCaffe(netFile, modelFile);
  m_net.setPreferableBackend(0);
  m_net.setPreferableTarget(0);

  //without checking the channel of input layer
  m_inputSize = inputSize;
  m_outputSize = outputSize;
}

Void TComOpencv::netForward(cv::Mat img)
{
	CHECK(!img.empty()) << "Unable to reconstruct image ";
	//!!!set number of channel is equal to the channel of input data
	m_numChannel = img.channels();
	preProcess(img);

	auto output_raw = m_net.forward(m_blobName);

	std::vector<int> siz;
	siz.push_back(m_outputSize.height);
	siz.push_back(m_outputSize.width);
	auto output = output_raw.reshape(1, siz);

	cv::Mat outputMat;
	if (m_numChannel == 3)
	{
		outputMat = cv::Mat(m_outputSize.height, m_outputSize.width, CV_32FC3);
		output.convertTo(outputMat, CV_32FC3);
	}
	else
	{
		outputMat = cv::Mat(m_outputSize.height, m_outputSize.width, CV_32FC1);
		output.convertTo(outputMat, CV_32FC1);
	}

	const Float* begin = (Float*)outputMat.data;
	const Float* end = begin + outputMat.total();

	m_NetOutput = std::vector<Float>(begin, end);

}

Void TComOpencv::preProcess( const cv::Mat& img)
{
	cv::Mat sample;
	sample = cv::Mat(m_inputSize.height, m_inputSize.width, CV_32FC1);
	for (int y = 0; y < m_inputSize.height; ++y)
	{
		for (int x = 0; x < m_inputSize.width; ++x)
		{
			sample.at<float>(y, x) = (float)img.at<ushort>(y, x);
		}
	}

  std::vector<int> siz;
  siz.push_back(1);
  siz.push_back(m_numChannel);
  siz.push_back(m_inputSize.height);
  siz.push_back(m_inputSize.width);

  cv::Mat sample_float;
  if (m_numChannel == 3)
  {
	  sample.reshape(1, siz).convertTo(sample_float, CV_32FC3, 1.0 / PixelRangeDouble);
  }
  else
  {
	  sample.reshape(1, siz).convertTo(sample_float, CV_32FC1, 1.0 / PixelRangeDouble);
  }

  m_net.setInput(sample_float);
}

#endif
//! \}