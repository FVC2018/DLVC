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

#ifndef __TENCCNNLOOPFILTER__
#define __TENCCNNLOOPFILTER__

#include "TLibCommon/TComCNNLoopFilter.h"
#include "TLibCommon/TComPic.h"

#include "TEncEntropy.h"
#include "TEncSbac.h"
#include "TLibCommon/TComBitCounter.h"

#if CNN_BASED_INLOOP_FILTER

class TEncCNNLoopFilter : public TComCNNLoopFilter
{
public:
  TEncCNNLoopFilter();
  virtual ~TEncCNNLoopFilter();
  Void createEncData();
  Void destroyEncData();
  Void initRDOCabacCoder(TEncSbac* pcRDGoOnSbacCoder, TComSlice* pcSlice);
  Void convNetProcess(TComPic* pPic);
private: //methods
 // Void decideBlkParams(TComPic* pic, TComPicYuv* srcYuv, TComPicYuv* resYuv, SAOBlkParam* reconParams, SAOBlkParam* codedParams,TComPicYuv *netRecYuv);
 
  Int64 getDistortion(Pel *srcBlk, Pel *orgBlk, Pel *netRecBlk, Int srcStride, Int orgStride, Int netRecStride, Int width, Int height);
  Void  decideBlkState(TComPic* pPic, TComPicYuv* orgYuv, TComPicYuv* srcYuv, TComPicYuv* netRecYuv, CNNLFBlkParam *& blkStats );
private: //members
  //for RDO
  TEncSbac**             m_pppcRDSbacCoder;
  TEncSbac*              m_pcRDGoOnSbacCoder;
  TEncBinCABACCounter**  m_pppcBinCoderCABAC;

  Double                 m_convNetDisabledRate[MAX_NUM_COMPONENT];
};
#endif

//! \}

#endif
