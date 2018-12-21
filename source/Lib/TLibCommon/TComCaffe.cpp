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

/** \file     TComCaffe.cpp
\brief    Caffe interface class
*/

#include "TComCaffe.h"

#if CNN_BASED_FRAC_INTERPOLATION
//#include "TComInterpolationFilter.h"
#include "CommonDef.h"
#endif
//using namespace caffe;
//! \ingroup TLibCommon
//! \{

TComCaffe::TComCaffe()
{

}

TComCaffe::~TComCaffe()
{

}

static bool init_glog = false;

Void TComCaffe::create(const string& netFile, const string& modelFile, Char* query_blob_name, Char* mode, Int iGPUid, cv::Size inputSize)
{
  if (!init_glog) {
    google::InitGoogleLogging("XXX");
    google::SetCommandLineOption("GLOG_minloglevel", "2");
    init_glog = true;
		openblas_set_num_threads(8);
  }
  Phase phase = TEST;
  if (strcmp(mode, "GPU") == 0)
  {
    Caffe::set_mode(Caffe::GPU);
    Caffe::SetDevice(iGPUid);
  }
  else
  {
    Caffe::set_mode(Caffe::CPU);
  }
  m_blobName = query_blob_name;
  m_net.reset(new Net<Float>(netFile, phase));
  m_net->CopyTrainedLayersFrom(modelFile);

  Blob<Float>* input_layer = m_net->input_blobs()[0];
  m_numChannel = input_layer->channels();
  CHECK(m_numChannel == 3 || m_numChannel == 1) << "Input layer should have 1 or 3 channels.";

  m_inputSize = inputSize;
}

Void TComCaffe::destroy(const string& netFile, const string& modelFile, Char* mode)
{
  Phase phase = TEST;
  if (strcmp(mode, "GPU") == 0)
  {
    Caffe::set_mode(Caffe::GPU);
  }
  else
  {
    Caffe::set_mode(Caffe::CPU);
  }
  m_net.reset(new Net<Float>(netFile, phase));
  m_net->CopyTrainedLayersFrom(modelFile);
}


Void TComCaffe::wrapInputLayer(std::vector<cv::Mat>* input_channels)
{
  Blob<Float>* input_layer = m_net->input_blobs()[0];

  Int width  = input_layer->width();
  Int height = input_layer->height();
  Float* input_data = input_layer->mutable_cpu_data();

  for (Int i = 0; i < input_layer->channels(); ++i) 
  {
    cv::Mat channel(height, width, CV_32FC1, input_data);
    input_channels->push_back(channel);
    input_data += width * height;
  }
}


Void TComCaffe::preProcess( const cv::Mat& img, std::vector<cv::Mat>* input_channels )
{
  /* Convert the input image to the input image format of the network. */
  cv::Mat sample;
  if (img.channels() == 3 && m_numChannel == 1)
  {
    cv::cvtColor(img, sample, cv::COLOR_BGR2GRAY);
  }
  else if (img.channels() == 4 && m_numChannel == 1)
  {
    cv::cvtColor(img, sample, cv::COLOR_BGRA2GRAY);
  }
  else if (img.channels() == 4 && m_numChannel == 3)
  {
    cv::cvtColor(img, sample, cv::COLOR_BGRA2BGR);
  }
  else if (img.channels() == 1 && m_numChannel == 3)
  {
    cv::cvtColor(img, sample, cv::COLOR_GRAY2BGR);
  }
  else
  {
    sample = img;
  }

  cv::Mat sample_resized;
  if (sample.size() != m_inputSize)
  {
    cv::resize(sample, sample_resized, m_inputSize);
  }
  else
  {
    sample_resized = sample;
  }

  cv::Mat sample_float;
  if (m_numChannel == 3)
  {
    sample_resized.convertTo(sample_float, CV_32FC3, 1.0 / PixelRangeDouble);
  }
  else
  {
    sample_resized.convertTo(sample_float, CV_32FC1, 1.0 / PixelRangeDouble);
  }

  /* This operation will write the separate BGR planes directly to the
  * input layer of the network because it is wrapped by the cv::Mat
  * objects in input_channels. */
  cv::split(sample_float, *input_channels);

  CHECK(reinterpret_cast<Float*>(input_channels->at(0).data) == m_net->input_blobs()[0]->cpu_data())
    << "Input channels are not wrapping the input layer of the network.";
}


Int TComCaffe::getBlobIndex( boost::shared_ptr< Net<Float> > & net, Char *query_blob_name )
{
  std::string str_query(query_blob_name);
  vector< string > const & blob_names = net->blob_names();

  for (UInt i = 0; i != blob_names.size(); ++i)
  {
    if (str_query == blob_names[i])
    {
      return i;
    }
  }
  LOG(FATAL) << "Unknown blob name: " << str_query;
  return -1;
}

Void TComCaffe::netForward(cv::Mat img)
{
  CHECK(!img.empty()) << "Unable to reconstruct image ";
  Blob<Float>* input_layer = m_net->input_blobs()[0];

  input_layer->Reshape( 1, m_numChannel, m_inputSize.height, m_inputSize.width );
  //printf("CNN Reshape begin\n");
  //system("pause");
  m_net->Reshape( );
  //printf("CNN Reshape end\n");
  //system("pause");
  std::vector<cv::Mat> input_channels;
  wrapInputLayer( &input_channels );
  preProcess(img, &input_channels);
  m_net->Forward( );

  Int uiBlobId = getBlobIndex(m_net, m_blobName);
  boost::shared_ptr<Blob<Float> > blob = getNet()->blobs()[uiBlobId];
  UInt uiNumData = blob->count();

  const Float* begin = blob->cpu_data();
  const Float* end   = begin + uiNumData;
  m_NetOutput = std::vector<Float>(begin, end);
}

std::vector<Float> TComCaffe::getNetOutputCenter( Int iMargin )
{
  Int iheight = m_inputSize.height;
  Int iwidth  = m_inputSize.width;
  Int iheightDst = iheight - 2 * iMargin;
  Int iwidthDst  = iwidth  - 2 * iMargin;

  std::vector<float>::iterator it = m_NetOutput.begin();
  it += (iMargin*iwidth + iMargin);

  vector<float> dst;
  for (Int y = 0; y < iheightDst; y++)
  {
    for (Int x = 0; x < iwidthDst; x++)
    {
      dst.push_back(it[y*iwidth + x]);
    }
  }
  return dst;
}
//! \}