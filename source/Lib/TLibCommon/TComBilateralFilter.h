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
#ifndef TComBilateralFilter_h
#define TComBilateralFilter_h

#include "TComPattern.h"
#include "CommonDef.h"

#define BILATERAL_FILTER_MAX_DENOMINATOR_PLUS_ONE  325
#define BITS_PER_DIV_LUT_ENTRY 14
#define SIGN_IF_NEG(x) (-((x)<0))                             // Returns -1 if x is negative, otherwise 0

class TComBilateralFilter
{
private:
  static TComBilateralFilter* m_bilateralFilterInstance;
  TComBilateralFilter();

public:
  static const Int SpatialSigmaValue;
  static const Int spatialSigmaBlockLengthOffsets[5];
  ~TComBilateralFilter();
  UShort** m_bilateralFilterTable;
  Int m_bilateralCenterWeightTable[5];
#if JVET_F0096_BILATERAL_FILTER
	bool m_initFlag;
	bool getInitFlag() const { return m_initFlag; };
	Void setInitFlag(bool b) { m_initFlag = b; }
#endif
  UInt divToMulOneOverN[BILATERAL_FILTER_MAX_DENOMINATOR_PLUS_ONE];
  UChar divToMulShift[BILATERAL_FILTER_MAX_DENOMINATOR_PLUS_ONE];
  Void createdivToMulLUTs();

  static TComBilateralFilter* instance();
  Void createBilateralFilterTable(Int qp);
  Void smoothBlockBilateralFilter(TComDataCU* pcCU, UInt uiWidth, UInt uiHeight, Short block[], Int length, Int isInterBlock, Int qp);
  Void smoothBlockBilateralFilter(TComDataCU* pcCU, UInt uiWidth, UInt uiHeight, Short dstRecon[], Int length, Int isInterBlock, Int qp, Pel *piReco, UInt uiSrcStride, UInt uiDstStride);
  Void bilateralFilterInter(TComDataCU *pcCU, UInt uiWidth, UInt uiHeight, Pel *piResi, UInt uiStrideRes, Pel *piPred, UInt uiPredStride, Pel *piReco, UInt uiRecStride, Int clipbd, Int qp);
  Void bilateralFilterIntra(TComDataCU *pcCU, UInt uiWidth, UInt uiHeight, Pel *piReco, UInt uiStride, Int qp);
};

#endif /* TComBilateralFilter_h */
