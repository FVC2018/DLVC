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

/** \file     TComCaffeHeader.h
*/

#ifndef __TCOMCAFFEHEADER__
#define __TCOMCAFFEHEADER__


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
#include "caffe/layers/scale_layer.hpp"
#include "caffe/layers/bias_layer.hpp"
#include "caffe/layers/deconv_layer.hpp"
namespace caffe
{
  extern INSTANTIATE_CLASS(InputLayer);
  extern INSTANTIATE_CLASS(InnerProductLayer);
  extern INSTANTIATE_CLASS(DropoutLayer);
  extern INSTANTIATE_CLASS(ConvolutionLayer);
  //REGISTER_LAYER_CLASS(Convolution);
  extern INSTANTIATE_CLASS(ReLULayer);
  //REGISTER_LAYER_CLASS(ReLU);
  extern INSTANTIATE_CLASS(PoolingLayer);
  //REGISTER_LAYER_CLASS(Pooling);
  extern INSTANTIATE_CLASS(LRNLayer);
  //REGISTER_LAYER_CLASS(LRN);
  extern INSTANTIATE_CLASS(SoftmaxLayer);
  //REGISTER_LAYER_CLASS(Softmax);
  extern INSTANTIATE_CLASS(ConcatLayer);
  extern INSTANTIATE_CLASS(EuclideanLossLayer);
  extern INSTANTIATE_CLASS(ScaleLayer);
  extern INSTANTIATE_CLASS(BiasLayer);
  extern INSTANTIATE_CLASS(DeconvolutionLayer);
}

#endif