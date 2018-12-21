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

/** \file     TComPicYuv.cpp
    \brief    picture YUV buffer class
*/

#include <cstdlib>
#include <assert.h>
#include <memory.h>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include "TComPicYuv.h"
#include "TLibVideoIO/TVideoIOYuv.h"

//! \ingroup TLibCommon
//! \{

TComPicYuv::TComPicYuv()
{
  for(UInt i=0; i<MAX_NUM_COMPONENT; i++)
  {
    m_apiPicBuf[i]    = NULL;   // Buffer (including margin)
    m_piPicOrg[i]     = NULL;    // m_apiPicBufY + m_iMarginLuma*getStride() + m_iMarginLuma
  }

  for(UInt i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
  {
    m_ctuOffsetInBuffer[i]=0;
    m_subCuOffsetInBuffer[i]=0;
  }

#if BARC
  for(UInt i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
  {
    m_ctuOffsetInBufferDscale[i]=0;
    m_subCuOffsetInBufferDscale[i]=0;
  }
#endif
  m_bIsBorderExtended = false;
}




TComPicYuv::~TComPicYuv()
{
}

#if BLOCK_BASED_BACKGROUND_REFERENCE
Void TComPicYuv::clear()
{
	Int chan = 0;
	for (; chan < getNumberValidComponents(); chan++)
	{
		const ComponentID ch = ComponentID(chan);
		memset(m_apiPicBuf[chan], 0, sizeof(Pel)*getStride(ch)       * getTotalHeight(ch));
	}
}
#endif


Void TComPicYuv::create ( const Int iPicWidth,                ///< picture width
                          const Int iPicHeight,               ///< picture height
                          const ChromaFormat chromaFormatIDC, ///< chroma format
                          const UInt uiMaxCUWidth,            ///< used for generating offsets to CUs. Can use iPicWidth if no offsets are required
                          const UInt uiMaxCUHeight,           ///< used for generating offsets to CUs. Can use iPicHeight if no offsets are required
                          const UInt uiMaxCUDepth,            ///< used for generating offsets to CUs. Can use 0 if no offsets are required
                          const Bool bUseMargin)              ///< if true, then a margin of uiMaxCUWidth+16 and uiMaxCUHeight+16 is created around the image.

{
  m_iPicWidth         = iPicWidth;
  m_iPicHeight        = iPicHeight;
  m_chromaFormatIDC   = chromaFormatIDC;
  m_iMarginX          = (bUseMargin?uiMaxCUWidth:0) + 16;   // for 16-byte alignment
  m_iMarginY          = (bUseMargin?uiMaxCUHeight:0) + 16;  // margin for 8-tap filter and infinite padding
  m_bIsBorderExtended = false;

  // assign the picture arrays and set up the ptr to the top left of the original picture
  {
    Int chan=0;
    for(; chan<getNumberValidComponents(); chan++)
    {
      const ComponentID ch=ComponentID(chan);
      m_apiPicBuf[chan] = (Pel*)xMalloc( Pel, getStride(ch)       * getTotalHeight(ch));
      m_piPicOrg[chan]  = m_apiPicBuf[chan] + (m_iMarginY >> getComponentScaleY(ch))   * getStride(ch)       + (m_iMarginX >> getComponentScaleX(ch));
    }
    for(;chan<MAX_NUM_COMPONENT; chan++)
    {
      m_apiPicBuf[chan] = NULL;
      m_piPicOrg[chan]  = NULL;
    }
  }


  const Int numCuInWidth  = m_iPicWidth  / uiMaxCUWidth  + (m_iPicWidth  % uiMaxCUWidth  != 0);
  const Int numCuInHeight = m_iPicHeight / uiMaxCUHeight + (m_iPicHeight % uiMaxCUHeight != 0);
  for(Int chan=0; chan<2; chan++)
  {
    const ComponentID ch=ComponentID(chan);
    const Int ctuHeight=uiMaxCUHeight>>getComponentScaleY(ch);
    const Int ctuWidth=uiMaxCUWidth>>getComponentScaleX(ch);
    const Int stride = getStride(ch);

    m_ctuOffsetInBuffer[chan] = new Int[numCuInWidth * numCuInHeight];

    for (Int cuRow = 0; cuRow < numCuInHeight; cuRow++)
    {
      for (Int cuCol = 0; cuCol < numCuInWidth; cuCol++)
      {
        m_ctuOffsetInBuffer[chan][cuRow * numCuInWidth + cuCol] = stride * cuRow * ctuHeight + cuCol * ctuWidth;
      }
    }
#if BARC
    m_ctuOffsetInBufferDscale[chan] = new Int[numCuInWidth * numCuInHeight];

    for (Int cuRow = 0; cuRow < numCuInHeight; cuRow++)
    {
      for (Int cuCol = 0; cuCol < numCuInWidth; cuCol++)
      {
        m_ctuOffsetInBufferDscale[chan][cuRow * numCuInWidth + cuCol] = stride * cuRow * ctuHeight / 2 + cuCol * ctuWidth /2;
      }
    }
#endif
    m_subCuOffsetInBuffer[chan] = new Int[(size_t)1 << (2 * uiMaxCUDepth)];

    const Int numSubBlockPartitions=(1<<uiMaxCUDepth);
    const Int minSubBlockHeight    =(ctuHeight >> uiMaxCUDepth);
    const Int minSubBlockWidth     =(ctuWidth  >> uiMaxCUDepth);

    for (Int buRow = 0; buRow < numSubBlockPartitions; buRow++)
    {
      for (Int buCol = 0; buCol < numSubBlockPartitions; buCol++)
      {
        m_subCuOffsetInBuffer[chan][(buRow << uiMaxCUDepth) + buCol] = stride  * buRow * minSubBlockHeight + buCol * minSubBlockWidth;
      }
    }
#if BARC
    m_subCuOffsetInBufferDscale[chan] = new Int[(size_t)1 << (2 * uiMaxCUDepth)];

    for (Int buRow = 0; buRow < numSubBlockPartitions; buRow++)
    {
      for (Int buCol = 0; buCol < numSubBlockPartitions; buCol++)
      {
        m_subCuOffsetInBufferDscale[chan][(buRow << uiMaxCUDepth) + buCol] = stride  * buRow * minSubBlockHeight / 2 + buCol * minSubBlockWidth / 2;
      }
    }
#endif
  }
  return;
}



Void TComPicYuv::destroy()
{
  for(Int chan=0; chan<MAX_NUM_COMPONENT; chan++)
  {
    m_piPicOrg[chan] = NULL;

    if( m_apiPicBuf[chan] )
    {
      xFree( m_apiPicBuf[chan] );
      m_apiPicBuf[chan] = NULL;
    }
  }

  for(UInt chan=0; chan<MAX_NUM_CHANNEL_TYPE; chan++)
  {
    if (m_ctuOffsetInBuffer[chan])
    {
      delete[] m_ctuOffsetInBuffer[chan];
      m_ctuOffsetInBuffer[chan] = NULL;
    }
    if (m_subCuOffsetInBuffer[chan])
    {
      delete[] m_subCuOffsetInBuffer[chan];
      m_subCuOffsetInBuffer[chan] = NULL;
    }
  }
#if BARC
  for(UInt chan=0; chan<MAX_NUM_CHANNEL_TYPE; chan++)
  {
    if (m_ctuOffsetInBufferDscale[chan])
    {
      delete[] m_ctuOffsetInBufferDscale[chan];
      m_ctuOffsetInBufferDscale[chan] = NULL;
    }
    if (m_subCuOffsetInBufferDscale[chan])
    {
      delete[] m_subCuOffsetInBufferDscale[chan];
      m_subCuOffsetInBufferDscale[chan] = NULL;
    }
  }
#endif
}



Void  TComPicYuv::copyToPic (TComPicYuv*  pcPicYuvDst
#if ALF_HM3_REFACTOR
  , ComponentID chDst
  , Bool bMarginIncluded 
#endif
  ) const
{
  assert( m_iPicWidth  == pcPicYuvDst->getWidth(COMPONENT_Y)  );
  assert( m_iPicHeight == pcPicYuvDst->getHeight(COMPONENT_Y) );
  assert( m_chromaFormatIDC == pcPicYuvDst->getChromaFormat() );

#if ALF_HM3_REFACTOR
  Int chStart = 0 , chEnd = getNumberValidComponents();
  if( chDst != MAX_NUM_COMPONENT )
  {
    chStart = chDst;
    chEnd = chStart + 1;
  }
  for(Int chan=chStart ; chan < chEnd ; chan++ )
#else
  for(Int chan=0; chan<getNumberValidComponents(); chan++)
#endif
  {
    const ComponentID ch=ComponentID(chan);
#if ALF_HM3_REFACTOR
    if( !bMarginIncluded )
    {
      const Pel * pSrc = getAddr( ch );
      Int nSrcStride = getStride( ch );
      Pel * pDst = pcPicYuvDst->getAddr( ch );
      Int nDstStride = pcPicYuvDst->getStride( ch );
      Int nSize = getWidth( ch ) * sizeof( Pel );

	  twoDimBufCopier
		  (
			  (UChar *)pDst,
			  nDstStride * sizeof(Pel),
			  (UChar *)pSrc,
			  nSrcStride * sizeof(Pel),
			  nSize,
			  getHeight(ch)
			  );
      //for( Int n = 0 ; n < getHeight( ch ) ; n++ , pSrc += nSrcStride , pDst += nDstStride )
      //{
      //  ::memcpy( pDst , pSrc , nSize );
      //}
    }
    else
#endif
    ::memcpy ( pcPicYuvDst->getBuf(ch), m_apiPicBuf[ch], sizeof (Pel) * getStride(ch) * getTotalHeight(ch));
  }
  return;
}


Void TComPicYuv::extendPicBorder (
#if ALF_HM3_REFACTOR
  Int nMargin 
#endif
  )
{
  if ( m_bIsBorderExtended )
  {
    return;
  }

  for(Int chan=0; chan<getNumberValidComponents(); chan++)
  {
    const ComponentID ch=ComponentID(chan);
    Pel *piTxt=getAddr(ch); // piTxt = point to (0,0) of image within bigger picture.
    const Int iStride=getStride(ch);
    const Int iWidth=getWidth(ch);
    const Int iHeight=getHeight(ch);
#if ALF_HM3_REFACTOR
    const Int iMarginX=nMargin>0?nMargin:getMarginX(ch);
    const Int iMarginY=nMargin>0?nMargin:getMarginY(ch);
#else
    const Int iMarginX=getMarginX(ch);
    const Int iMarginY=getMarginY(ch);
#endif

    Pel*  pi = piTxt;
    // do left and right margins
    for (Int y = 0; y < iHeight; y++)
    {
      for (Int x = 0; x < iMarginX; x++ )
      {
        pi[ -iMarginX + x ] = pi[0];
        pi[    iWidth + x ] = pi[iWidth-1];
      }
      pi += iStride;
    }

    // pi is now the (0,height) (bottom left of image within bigger picture
    pi -= (iStride + iMarginX);
    // pi is now the (-marginX, height-1)
    for (Int y = 0; y < iMarginY; y++ )
    {
      ::memcpy( pi + (y+1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX<<1)) );
    }

    // pi is still (-marginX, height-1)
    pi -= ((iHeight-1) * iStride);
    // pi is now (-marginX, 0)
    for (Int y = 0; y < iMarginY; y++ )
    {
      ::memcpy( pi - (y+1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX<<1)) );
    }
  }

  m_bIsBorderExtended = true;
}

#if BARC_FUNC
Void TComPicYuv::extendPicBorderBARC (
#if ALF_HM3_REFACTOR
  Int nMargin 
#endif
  )
{
  if ( m_bIsBorderExtended )
  {
    return;
  }

  for(Int chan=0; chan<getNumberValidComponents(); chan++)
  {
    const ComponentID ch=ComponentID(chan);
    Pel *piTxt=getAddr(ch); // piTxt = point to (0,0) of image within bigger picture.
    const Int iStride=getStride(ch);
    const Int iWidth=getWidth(ch);
    const Int iHeight=getHeight(ch);
#if ALF_HM3_REFACTOR
    const Int iMarginX=nMargin>0?nMargin:getMarginX(ch);
    const Int iMarginY=nMargin>0?nMargin:getMarginY(ch);
#else
    const Int iMarginX=getMarginX(ch);
    const Int iMarginY=getMarginY(ch);
#endif

    Pel*  pi = piTxt;
    // do left and right margins
    for (Int y = 0; y < iHeight; y++)
    {
      for (Int x = 0; x < iMarginX; x++ )
      {
        pi[ -iMarginX + x ] = pi[0];
        pi[    iWidth + x ] = pi[iWidth-1];
      }
      pi += iStride;
    }

    // pi is now the (0,height) (bottom left of image within bigger picture
    pi -= (iStride + iMarginX);
    // pi is now the (-marginX, height-1)
    for (Int y = 0; y < iMarginY; y++ )
    {
      ::memcpy( pi + (y+1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX<<1)) );
    }

    // pi is still (-marginX, height-1)
    pi -= ((iHeight-1) * iStride);
    // pi is now (-marginX, 0)
    for (Int y = 0; y < iMarginY; y++ )
    {
      ::memcpy( pi - (y+1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX<<1)) );
    }
  }

  m_bIsBorderExtended = true;
}

#endif

// NOTE: This function is never called, but may be useful for developers.
Void TComPicYuv::dump (const Char* pFileName, const BitDepths &bitDepths, Bool bAdd) const
{
  FILE* pFile;
  if (!bAdd)
  {
    pFile = fopen (pFileName, "wb");
  }
  else
  {
    pFile = fopen (pFileName, "ab");
  }


  for(Int chan = 0; chan < getNumberValidComponents(); chan++)
  {
    const ComponentID  ch     = ComponentID(chan);
    const Int          shift  = bitDepths.recon[toChannelType(ch)] - 8;
    const Int          offset = (shift>0)?(1<<(shift-1)):0;
    const Pel         *pi     = getAddr(ch);
    const Int          stride = getStride(ch);
    const Int          height = getHeight(ch);
    const Int          width  = getWidth(ch);

    for (Int y = 0; y < height; y++ )
    {
      for (Int x = 0; x < width; x++ )
      {
        UChar uc = (UChar)Clip3<Pel>(0, 255, (pi[x]+offset)>>shift);
        fwrite( &uc, sizeof(UChar), 1, pFile );
      }
      pi += stride;
    }
  }

  fclose(pFile);
}

#if EXTENSION_360_VIDEO
Void TComPicYuv::createWithoutCUInfo ( const Int picWidth,                 ///< picture width
                                       const Int picHeight,                ///< picture height
                                       const ChromaFormat chromaFormatIDC, ///< chroma format
                                       const Bool bUseMargin,              ///< if true, then a margin of uiMaxCUWidth+16 and uiMaxCUHeight+16 is created around the image.
                                       const UInt maxCUWidth,              ///< used for margin only
                                       const UInt maxCUHeight)             ///< used for margin only
{
  destroy();

  m_iPicWidth          = picWidth;
  m_iPicHeight         = picHeight;
  m_chromaFormatIDC   = chromaFormatIDC;
  m_iMarginX          = (bUseMargin?maxCUWidth:0) + 16;   // for 16-byte alignment
  m_iMarginY          = (bUseMargin?maxCUHeight:0) + 16;  // margin for 8-tap filter and infinite padding
  m_bIsBorderExtended = false;

  // assign the picture arrays and set up the ptr to the top left of the original picture
  for(UInt comp=0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID ch=ComponentID(comp);
    m_apiPicBuf[comp] = (Pel*)xMalloc( Pel, getStride(ch) * getTotalHeight(ch));
    m_piPicOrg[comp]  = m_apiPicBuf[comp] + (m_iMarginY >> getComponentScaleY(ch)) * getStride(ch) + (m_iMarginX >> getComponentScaleX(ch));
  }
  // initialize pointers for unused components to NULL
  for(UInt comp=getNumberValidComponents();comp<MAX_NUM_COMPONENT; comp++)
  {
    m_apiPicBuf[comp] = NULL;
    m_piPicOrg[comp]  = NULL;
  }

  for(Int chan=0; chan<MAX_NUM_CHANNEL_TYPE; chan++)
  {
    m_ctuOffsetInBuffer[chan]   = NULL;
    m_subCuOffsetInBuffer[chan] = NULL;
  }
}
#endif

//! \}
