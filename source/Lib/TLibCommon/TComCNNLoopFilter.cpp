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

#include "TComCNNLoopFilter.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "TComSplitImg.h"

#if CNN_BASED_INLOOP_FILTER

TComCNNLoopFilter::TComCNNLoopFilter()
{
  m_tempPicYuv = NULL;
}

Void TComCNNLoopFilter::destroy()
{
  if (m_tempPicYuv)
  {
    m_tempPicYuv->destroy();
    delete m_tempPicYuv;
    m_tempPicYuv = NULL;
  }
}

TComCNNLoopFilter::~TComCNNLoopFilter()
{
  destroy();
}

Void TComCNNLoopFilter::create(Int picWidth, Int picHeight, ChromaFormat format, UInt maxCUWidth, UInt maxCUHeight, UInt maxCUDepth, UInt lumaBitShift, UInt chromaBitShift)
{
  destroy();

  m_picWidth = picWidth;
  m_picHeight = picHeight;
  m_chromaFormatIDC = format;
  m_maxCUWidth = maxCUWidth;
  m_maxCUHeight = maxCUHeight;

  m_numCTUInWidth = (m_picWidth / m_maxCUWidth) + ((m_picWidth % m_maxCUWidth) ? 1 : 0);
  m_numCTUInHeight = (m_picHeight / m_maxCUHeight) + ((m_picHeight % m_maxCUHeight) ? 1 : 0);
  m_numCTUsPic = m_numCTUInHeight*m_numCTUInWidth;

  //temporary picture buffer
  if (!m_tempPicYuv)
  {
    m_tempPicYuv = new TComPicYuv;
    m_tempPicYuv->create(m_picWidth, m_picHeight, m_chromaFormatIDC, m_maxCUWidth, m_maxCUWidth, maxCUDepth, true);
  }
}

Void TComCNNLoopFilter::reconstructCTU(Pel* srcBlk, Pel* netBlk, Int srcStride, Int netStride, Int width, Int height)
{
  Int x, y;
  Pel* srcLine = srcBlk;
  Pel* netLine = netBlk;

  for (y = 0; y< height; y++)
  {
    for (x = 0; x< width; x++)
    {
      srcLine[x] = netLine[x];
    }
    srcLine += srcStride;
    netLine += netStride;
  }
}

Void TComCNNLoopFilter::dealCTU(Int ctu, TComPicYuv* srcYuv, TComPicYuv* netYuv, CNNLFBlkParam* netblkParam, TComPic* pPic)
{
  const Int numberOfComponents = getNumberValidComponents(m_chromaFormatIDC);
  Int yPos = (ctu / m_numCTUInWidth)*m_maxCUHeight;
  Int xPos = (ctu % m_numCTUInWidth)*m_maxCUWidth;
  Int height = (yPos + m_maxCUHeight > m_picHeight) ? (m_picHeight - yPos) : m_maxCUHeight;
  Int width = (xPos + m_maxCUWidth  > m_picWidth) ? (m_picWidth - xPos) : m_maxCUWidth;
	for (Int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    const ComponentID component = ComponentID(compIdx);

    if (netblkParam[ctu][compIdx])
    {
      const UInt componentScaleX = getComponentScaleX(component, pPic->getChromaFormat());
      const UInt componentScaleY = getComponentScaleY(component, pPic->getChromaFormat());

      Int  blkWidth = (width >> componentScaleX);
      Int  blkHeight = (height >> componentScaleY);
      Int  blkXPos = (xPos >> componentScaleX);
      Int  blkYPos = (yPos >> componentScaleY);

      Int  srcStride = srcYuv->getStride(component);
      Pel* srcBlk = srcYuv->getAddr(component) + blkYPos*srcStride + blkXPos;

      Int  netStride = netYuv->getStride(component);
      Pel* netBlk = netYuv->getAddr(component) + blkYPos*netStride + blkXPos;
      
      reconstructCTU(srcBlk, netBlk, srcStride, netStride, blkWidth, blkHeight );
    }
  } //compIdx
}


Void TComCNNLoopFilter::convNetProcess(TComPic* pDecPic)
{
  string DLframework = pDecPic->getSlice(pDecPic->getCurrSliceIdx())->getSPS()->getDLframework();
  TComPicYuv* resYuv = pDecPic->getPicYuvRec();
  Int maxWidth = DEFAULT_INPUT_BLOCK_SIZE;
  Int maxHeight = DEFAULT_INPUT_BLOCK_SIZE;
#if USE_TENSORFLOW
  if (DLframework == "Tensorflow")
  {
    maxWidth = 512;
    maxHeight = 512;
  }
#endif

  TComPicYuv* netRecYuv = m_tempPicYuv;
  resYuv->copyToPic(netRecYuv);
  netRecYuv->setBorderExtension(false);
  netRecYuv->extendPicBorder();
  TComSplitImg splitImg(maxWidth, maxHeight);

  if ((DLframework == "Caffe") || (DLframework == "Opencv"))
  {
    splitImg.postProcess(pDecPic, netRecYuv, DLframework);
  }
#if USE_TENSORFLOW
  else if (DLframework == "Tensorflow")
  {
    splitImg.postProcess_tf(pDecPic, netRecYuv);
  }
#endif
  Int numLcusForConvNetOff[MAX_NUM_COMPONENT] = { 0, 0, 0 };

  for (Int ctu = 0; ctu < m_numCTUsPic; ctu++)
  {
		CNNLFBlkParam* netblkParam = pDecPic->getPicSym()->getConvNetBlkParam();

    dealCTU(ctu, resYuv, netRecYuv, netblkParam, pDecPic );
    for (Int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
    {
      if (netblkParam[ctu][compIdx] == false)
      {
        numLcusForConvNetOff[compIdx]++;
      }
    }
  } //ctu

 /* for (Int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  { 
    cout << "component " << compIdx << " disable ratio:" << (Double)numLcusForConvNetOff[compIdx] / (Double)m_numCTUsPic<< endl;
  }*/
}

#endif
