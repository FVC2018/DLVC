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

#ifndef __TCOMCNNLOOPFILTER__
#define __TCOMCNNLOOPFILTER__

#include "CommonDef.h"
#include "TComPic.h"

#if CNN_BASED_INLOOP_FILTER

class TComCNNLoopFilter
{
public:
  TComCNNLoopFilter();
  virtual ~TComCNNLoopFilter();
  Void convNetProcess(TComPic* pDecPic);
  Void create(Int picWidth, Int picHeight, ChromaFormat format, UInt maxCUWidth, UInt maxCUHeight, UInt maxCUDepth, UInt lumaBitShift, UInt chromaBitShift);
  Void destroy();
protected:
  Void dealCTU(Int ctu, TComPicYuv* srcYuv, TComPicYuv* resYuv, CNNLFBlkParam* netblkParam, TComPic* pPic);
  Void reconstructCTU(Pel* srcBlk, Pel* resBlk, Int srcStride, Int resStride, Int width, Int height );
protected:

  TComPicYuv*   m_tempPicYuv; //temporary buffer
  Int m_picWidth;
  Int m_picHeight;
  Int m_maxCUWidth;
  Int m_maxCUHeight;
  Int m_numCTUInWidth;
  Int m_numCTUInHeight;
  Int m_numCTUsPic;
  ChromaFormat m_chromaFormatIDC;
};
#endif
//! \}
#endif

