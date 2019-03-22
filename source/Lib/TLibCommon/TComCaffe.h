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

/** \file     TComCaffe.h
\brief    Caffe interface class (header)
*/

#ifndef __TCOMCAFFE__
#define __TCOMCAFFE__
#include "CommonDef.h"
#ifdef __GNUC__
#include <openblas/cblas.h>
#endif
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
#include "caffe/common.hpp"
#include "caffe/layers/input_layer.hpp"
#include "caffe/layers/inner_product_layer.hpp"
#include "caffe/layers/dropout_layer.hpp"
#include "caffe/layers/conv_layer.hpp"
#include "caffe/layers/relu_layer.hpp"

#include "caffe/layers/pooling_layer.hpp"
#include "caffe/layers/lrn_layer.hpp"
#include "caffe/layers/softmax_layer.hpp"
#include "caffe/layers/concat_layer.hpp"
#include "caffe/layers/euclidean_loss_layer.hpp"

using namespace caffe;
using std::string;

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// basic Caffe interface class

class TComCaffe
{
public:
  TComCaffe();
  ~TComCaffe();

  Void create( const string& netFile, const string& modelFile, Char* query_blob_name, Char* mode, Int iGPUid, cv::Size inputSize );

  Void destroy(const string& netFile, const string& modelFile, Char* mode);

  Void setInputSize(cv::Size inputSize){ m_inputSize = inputSize; }

  Void netForward(cv::Mat img);

  boost::shared_ptr<Net<Float> >& getNet(){ return m_net; };

  Int getBlobIndex(boost::shared_ptr< Net<Float> >& net, Char *query_blob_name);

  std::vector<Float> getNetOutput(){ return m_NetOutput; }
  std::vector<Float> getNetOutputCenter(Int iMargin);

private:
  Void wrapInputLayer( std::vector<cv::Mat>* input_channels );
  Void preProcess( const cv::Mat& img, std::vector<cv::Mat>* input_channels );

private:
  boost::shared_ptr< Net<Float> >  m_net;
  cv::Size                         m_inputSize;
  Int                              m_numChannel;
  Char*                            m_blobName;
  std::vector<Float>               m_NetOutput;
};// END CLASS DEFINITION TCOMCAFFE

//! \}

#endif // __TCOMCAFFE__
