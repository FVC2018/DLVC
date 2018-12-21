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

#include "TEncCNNLoopFilter.h"
#include "TLibCommon/TComSplitImg.h"

#if CNN_BASED_INLOOP_FILTER

TEncCNNLoopFilter::TEncCNNLoopFilter()
{
}

TEncCNNLoopFilter::~TEncCNNLoopFilter()
{
}

Void TEncCNNLoopFilter::createEncData()
{
  ::memset(m_convNetDisabledRate, 0, sizeof(m_convNetDisabledRate));
 }

Void TEncCNNLoopFilter::convNetProcess(TComPic* pPic)
{
  string DLframework = pPic->getSlice(pPic->getCurrSliceIdx())->getSPS()->getDLframework();
  TComPicYuv* orgYuv = pPic->getPicYuvOrg();
  TComPicYuv* resYuv = pPic->getPicYuvRec();
  Int maxBlkWidth = DEFAULT_INPUT_BLOCK_SIZE;
  Int maxBlkHeight = DEFAULT_INPUT_BLOCK_SIZE;
#if USE_TENSORFLOW
  if (DLframework == "Tensorflow")
  {
    maxBlkWidth = 512;
    maxBlkHeight = 512;
  }
#endif
  TComPicYuv* netRecYuv = m_tempPicYuv;
  resYuv->copyToPic(netRecYuv);
  netRecYuv->setBorderExtension(false);
  netRecYuv->extendPicBorder();
  TComSplitImg splitImg(maxBlkWidth, maxBlkHeight);
  
  if ((DLframework == "Caffe") || (DLframework == "Opencv"))
  {
    splitImg.postProcess(pPic, netRecYuv, DLframework);
  }
#if USE_TENSORFLOW
  else if(DLframework == "Tensorflow")
  {
    splitImg.postProcess_tf(pPic, netRecYuv);
  }
#endif

	CNNLFBlkParam* reconParams = pPic->getPicSym()->getConvNetBlkParam();
  decideBlkState(pPic, orgYuv, resYuv, netRecYuv, reconParams);
}

Void TEncCNNLoopFilter::decideBlkState(TComPic* pPic, TComPicYuv* orgYuv, TComPicYuv* srcYuv, TComPicYuv* netRecYuv, CNNLFBlkParam*& blkStats )
{
  const Int numberOfComponents = getNumberValidComponents(m_chromaFormatIDC);
  Int64 cost = 0;
 
  for (Int ctu = 0; ctu < m_numCTUsPic; ctu++)
  {
    Int yPos = (ctu / m_numCTUInWidth)*m_maxCUHeight;
    Int xPos = (ctu % m_numCTUInWidth)*m_maxCUWidth;
    Int height = (yPos + m_maxCUHeight > m_picHeight) ? (m_picHeight - yPos) : m_maxCUHeight;
    Int width = (xPos + m_maxCUWidth > m_picWidth) ? (m_picWidth - xPos) : m_maxCUWidth;

    //luma
    ComponentID component = COMPONENT_Y;

    const UInt componentScaleX = getComponentScaleX(component, pPic->getChromaFormat());
    const UInt componentScaleY = getComponentScaleY(component, pPic->getChromaFormat());

    Int  srcStride = srcYuv->getStride(component);
    Pel* srcBlk = srcYuv->getAddr(component) + ((yPos >> componentScaleY) * srcStride) + (xPos >> componentScaleX);

    Int  orgStride = orgYuv->getStride(component);
    Pel* orgBlk = orgYuv->getAddr(component) + ((yPos >> componentScaleY) * orgStride) + (xPos >> componentScaleX);

    Int  netRecStride = netRecYuv->getStride(component);
    Pel* netRecBlk = netRecYuv->getAddr(component) + ((yPos >> componentScaleY) * netRecStride) + (xPos >> componentScaleX);

    cost= getDistortion(srcBlk, orgBlk, netRecBlk, srcStride, orgStride, netRecStride, (width >> componentScaleX), (height >> componentScaleY));
    if (cost > 0)
    {
      blkStats[ctu][component] = true;
    }
    else
    {
      blkStats[ctu][component] = false;
    }
	
    if (blkStats[ctu][component])
    {
      reconstructCTU(srcBlk, netRecBlk, srcStride, netRecStride, (width >> componentScaleX), (height >> componentScaleY));
    }
     //chroma
	  cost = 0;
    for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberOfComponents; componentIndex++)
    {
      ComponentID component = ComponentID(componentIndex);
      const UInt componentScaleX = getComponentScaleX(component, pPic->getChromaFormat());
      const UInt componentScaleY = getComponentScaleY(component, pPic->getChromaFormat());

      Int  srcStride = srcYuv->getStride(component);
      Pel* srcBlk = srcYuv->getAddr(component) + ((yPos >> componentScaleY) * srcStride) + (xPos >> componentScaleX);

      Int  orgStride = orgYuv->getStride(component);
      Pel* orgBlk = orgYuv->getAddr(component) + ((yPos >> componentScaleY) * orgStride) + (xPos >> componentScaleX);

      Int  netRecStride = netRecYuv->getStride(component);
      Pel* netRecBlk = netRecYuv->getAddr(component) + ((yPos >> componentScaleY) * netRecStride) + (xPos >> componentScaleX);

      cost += getDistortion(srcBlk, orgBlk, netRecBlk, srcStride, orgStride, netRecStride, (width >> componentScaleX), (height >> componentScaleY));
    }
   
    if (cost > 0)
    {
      blkStats[ctu][COMPONENT_Cb] = true;
      blkStats[ctu][COMPONENT_Cr] = true;
    }
    else
    {
      blkStats[ctu][COMPONENT_Cb] = false;
      blkStats[ctu][COMPONENT_Cr] = false;
    }
       
    for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberOfComponents; componentIndex++)
    {
      ComponentID component = ComponentID(componentIndex);
      const UInt componentScaleX = getComponentScaleX(component, pPic->getChromaFormat());
      const UInt componentScaleY = getComponentScaleY(component, pPic->getChromaFormat());

      Int  srcStride = srcYuv->getStride(component);
      Pel* srcBlk = srcYuv->getAddr(component) + ((yPos >> componentScaleY) * srcStride) + (xPos >> componentScaleX);

      Int  orgStride = orgYuv->getStride(component);
      Pel* orgBlk = orgYuv->getAddr(component) + ((yPos >> componentScaleY) * orgStride) + (xPos >> componentScaleX);

      Int  netRecStride = netRecYuv->getStride(component);
      Pel* netRecBlk = netRecYuv->getAddr(component) + ((yPos >> componentScaleY) * netRecStride) + (xPos >> componentScaleX);

      if (blkStats[ctu][component])
      {
        reconstructCTU(srcBlk, netRecBlk, srcStride, netRecStride, (width >> componentScaleX), (height >> componentScaleY));
      }
    }
  }

  Int numLcusForConvNetOff[MAX_NUM_COMPONENT];
  for (Int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    numLcusForConvNetOff[compIdx] = 0;
    for (Int ctu = 0; ctu< m_numCTUsPic; ctu++)
    {
      if (blkStats[ctu][compIdx] == false)
      {
        numLcusForConvNetOff[compIdx]++;
      }
    }
  }
  for (Int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    m_convNetDisabledRate[compIdx] = (Double)numLcusForConvNetOff[compIdx] / (Double)m_numCTUsPic;
  }
}

Int64 TEncCNNLoopFilter::getDistortion(Pel *srcBlk, Pel *orgBlk, Pel *netRecBlk, Int srcStride, Int orgStride, Int netRecStride, Int width, Int height)
{
  Int x, y;
  Int64 diff[2] = {0,0};
  for (y = 0; y < height; y++)
  {
    for (x = 0; x < width; x++)
    {
      diff[0] += (orgBlk[x] - srcBlk[x])*(orgBlk[x] - srcBlk[x]);
      diff[1] += (orgBlk[x] - netRecBlk[x])*(orgBlk[x] - netRecBlk[x]);
    }
    orgBlk += orgStride;
    srcBlk += srcStride;
    netRecBlk += netRecStride;
  }
  return (diff[0] - diff[1]);
}

#endif
//! \}


