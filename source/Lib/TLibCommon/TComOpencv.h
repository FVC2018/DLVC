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

/** \file     TComOpencv.h
\brief    Caffe interface class (header)
*/

#ifndef __TCOMOPENCV__
#define __TCOMOPENCV__
#include "CommonDef.h"

#if frameworkCombine
#include "opencv2/dnn/dnn.hpp"
#include "opencv2/dnn/all_layers.hpp"
#include <caffe/caffe.hpp>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

//using namespace caffe;
using std::string;

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// basic Opencv interface class

class TComOpencv
{
public:
	TComOpencv();
	~TComOpencv();

	Void create(const string& netFile, const string& modelFile, Char* query_blob_name, cv::Size inputSize, cv::Size outputSize);
	Void netForward(cv::Mat img);
	Void setInputSize(cv::Size inputSize) { m_inputSize = inputSize; }
	Void setOutputSize(cv::Size outputSize) { m_outputSize = outputSize; }

	cv::dnn::Net& getNet() { return m_net; };

	std::vector<Float> getNetOutput() { return m_NetOutput; }

private:
	Void preProcess(const cv::Mat& img);

private:
	cv::dnn::Net					 m_net;
	cv::Size                         m_inputSize;
	cv::Size		                 m_outputSize;
	Int                              m_numChannel;
	Char*                            m_blobName;
	std::vector<Float>			     m_NetOutput;
};// END CLASS DEFINITION TCOMOPENCV
#endif
  //! \}

#endif // __TCOMOPENCV__
