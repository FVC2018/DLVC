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

/** \file     TComYuv.cpp
    \brief    general YUV buffer class
    \todo     this should be merged with TComPicYuv
*/

#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

#include "CommonDef.h"
#include "TComYuv.h"
#include "TComInterpolationFilter.h"

//! \ingroup TLibCommon
//! \{

TComYuv::TComYuv()
{
  for(Int comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    m_apiBuf[comp] = NULL;
  }
}

TComYuv::~TComYuv()
{
}

Void TComYuv::create( UInt iWidth, UInt iHeight, ChromaFormat chromaFormatIDC )
{
  // set width and height
  m_iWidth   = iWidth;
  m_iHeight  = iHeight;
  m_chromaFormatIDC = chromaFormatIDC;

  for(Int comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    // memory allocation
    m_apiBuf[comp]  = (Pel*)xMalloc( Pel, getWidth(ComponentID(comp))*getHeight(ComponentID(comp)) );
  }
}

Void TComYuv::destroy()
{
  // memory free
  for(Int comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    if (m_apiBuf[comp]!=NULL)
    {
      xFree( m_apiBuf[comp] );
      m_apiBuf[comp] = NULL;
    }
  }
}

Void TComYuv::clear()
{
  for(Int comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    if (m_apiBuf[comp]!=NULL)
    {
      ::memset( m_apiBuf[comp], 0, ( getWidth(ComponentID(comp)) * getHeight(ComponentID(comp))  )*sizeof(Pel) );
    }
  }
}




Void TComYuv::copyToPicYuv   ( TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth, const UInt uiPartIdx
#if BARC
  ,const Bool isBARC
#endif
) const
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    copyToPicComponent  ( ComponentID(comp), pcPicYuvDst, ctuRsAddr, uiAbsZorderIdx, uiPartDepth, uiPartIdx
#if BARC
      ,isBARC
#endif
    );
  }
}

Void TComYuv::copyToPicComponent  ( const ComponentID compID, TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth, const UInt uiPartIdx
#if BARC
  , const Bool isBARC
#endif
) const
{
  const Int iWidth  = getWidth(compID) >>uiPartDepth;
  const Int iHeight = getHeight(compID)>>uiPartDepth;

  const Pel* pSrc     = getAddr(compID, uiPartIdx, iWidth);
        Pel* pDst     = NULL;
#if BARC
        if (isBARC)
          pDst = pcPicYuvDst->getAddrDscale(compID, ctuRsAddr, uiAbsZorderIdx);
        else
#endif
          pDst = pcPicYuvDst->getAddr(compID, ctuRsAddr, uiAbsZorderIdx);
  const UInt  iSrcStride  = getStride(compID);
  const UInt  iDstStride  = pcPicYuvDst->getStride(compID);

  //for ( Int y = iHeight; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, sizeof(Pel)*iWidth);
  //  pDst += iDstStride;
  //  pSrc += iSrcStride;
  //}
  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  iWidth * sizeof(Pel),
		  iHeight
		  );
}




Void TComYuv::copyFromPicYuv   ( const TComPicYuv* pcPicYuvSrc, const UInt ctuRsAddr, const UInt uiAbsZorderIdx 
#if BARC
  , const Bool isBARC
#endif
)
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    copyFromPicComponent  ( ComponentID(comp), pcPicYuvSrc, ctuRsAddr, uiAbsZorderIdx
#if BARC
      , isBARC
#endif
    );
  }
}

Void TComYuv::copyFromPicComponent  ( const ComponentID compID, const TComPicYuv* pcPicYuvSrc, const UInt ctuRsAddr, const UInt uiAbsZorderIdx 
#if BARC
  , const Bool isBARC
#endif
)
{
        Pel* pDst     = getAddr(compID);
  const Pel* pSrc     = pcPicYuvSrc->getAddr ( compID, ctuRsAddr, uiAbsZorderIdx );
#if BARC
  const Pel* pSrcDscale = pcPicYuvSrc->getAddrDscale(compID, ctuRsAddr, uiAbsZorderIdx);
#endif

  const UInt iDstStride  = getStride(compID);
  const UInt iSrcStride  = pcPicYuvSrc->getStride(compID);
  const Int  iWidth=getWidth(compID);
  const Int  iHeight=getHeight(compID);
#if BARC
  if (isBARC)
	  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrcDscale,
		  iSrcStride * sizeof(Pel),
		  iWidth * sizeof(Pel),
		  iHeight
		  );
  else
#endif
  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  iWidth * sizeof(Pel),
		  iHeight
		  );

//  for (Int y = iHeight; y != 0; y-- )
//  {
//#if BARC
//    if (isBARC)
//      ::memcpy(pDst, pSrcDscale, sizeof(Pel)*iWidth);
//    else
//#endif
//    ::memcpy( pDst, pSrc, sizeof(Pel)*iWidth);
//    pDst += iDstStride;
//#if BARC
//    if (isBARC)
//      pSrcDscale += iSrcStride;
//    else
//#endif
//    pSrc += iSrcStride;
//  }
}



#if JVET_D0033_ADAPTIVE_CLIPPING
// encoder only for inter skip
Void TComYuv::clipToPartYuv(TComYuv* pcYuvDst, const UInt uiDstPartIdx, const BitDepths &bitDepths) const
{

    for (Int comp = 0; comp<getNumberValidComponents(); comp++)
    {
        clipToPartComponent(ComponentID(comp), pcYuvDst, uiDstPartIdx);
    }
}

Void TComYuv::clipToPartComponent(const ComponentID compID, TComYuv* pcYuvDst, const UInt uiDstPartIdx) const
{
    const Pel* pSrc = getAddr(compID);
    Pel* pDst = pcYuvDst->getAddr(compID, uiDstPartIdx);


    const UInt iSrcStride = getStride(compID);
    const UInt iDstStride = pcYuvDst->getStride(compID);
    const Int  iWidth = getWidth(compID);
    const Int  iHeight = getHeight(compID);

    for (Int y = iHeight; y != 0; y--)
    {
        for (Int x = 0; x<iWidth; x ++)
        {
            pDst[x] = ClipA(pSrc[x], compID);
        }
        pDst += iDstStride;
        pSrc += iSrcStride;
    }
}

// QT BT
Void TComYuv::clipToPicYuv   ( TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth, const UInt uiPartIdx ) const
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    clipToPicComponent  ( ComponentID(comp), pcPicYuvDst, ctuRsAddr, uiAbsZorderIdx, uiPartDepth, uiPartIdx );
  }
}

Void TComYuv::clipToPicComponent  ( const ComponentID compID, TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth, const UInt uiPartIdx ) const
{
    const Int iWidth  = getWidth(compID) >>uiPartDepth;
    const Int iHeight = getHeight(compID)>>uiPartDepth;

    const Pel* pSrc     = getAddr(compID, uiPartIdx, iWidth);
    Pel* pDst     = pcPicYuvDst->getAddr ( compID, ctuRsAddr, uiAbsZorderIdx );

    const UInt  iSrcStride  = getStride(compID);
    const UInt  iDstStride  = pcPicYuvDst->getStride(compID);

    for ( Int y = iHeight; y != 0; y-- )
    {
        for (Int x = 0; x<iWidth; x ++)
        {
            pDst[x] = ClipA(pSrc[x],  compID);
        }
        pDst += iDstStride;
        pSrc += iSrcStride;
    }
}


// decoder only for no residual inter
Void TComYuv::clipPartToPartYuv(TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidth, const UInt iHeight) const
{
    for (Int comp = 0; comp<getNumberValidComponents(); comp++)
    {
        clipPartToPartComponent(ComponentID(comp),
                                pcYuvDst,
                                uiPartIdx,
                                iWidth >> getComponentScaleX(ComponentID(comp)),
                                iHeight >> getComponentScaleY(ComponentID(comp)));
    }
}

Void TComYuv::clipPartToPartComponent(const ComponentID compID, TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidthComponent, const UInt iHeightComponent) const
{
    const Pel* pSrc = getAddr(compID, uiPartIdx);
    Pel* pDst = pcYuvDst->getAddr(compID, uiPartIdx);
    const UInt  iSrcStride = getStride(compID);
    const UInt  iDstStride = pcYuvDst->getStride(compID);
    for (UInt y = iHeightComponent; y != 0; y--)
    {
        for (UInt x = 0; x<iWidthComponent; x++)
        {
            pDst[x] = ClipA(pSrc[x],  compID);
        }
        pSrc += iSrcStride;
        pDst += iDstStride;
    }
}
#endif

Void TComYuv::copyToPartYuv( TComYuv* pcYuvDst, const UInt uiDstPartIdx ) const
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    copyToPartComponent  ( ComponentID(comp), pcYuvDst, uiDstPartIdx );
  }
}

Void TComYuv::copyToPartComponent( const ComponentID compID, TComYuv* pcYuvDst, const UInt uiDstPartIdx ) const
{
  const Pel* pSrc     = getAddr(compID);
        Pel* pDst     = pcYuvDst->getAddr( compID, uiDstPartIdx );

  const UInt iSrcStride  = getStride(compID);
  const UInt iDstStride  = pcYuvDst->getStride(compID);
  const Int  iWidth=getWidth(compID);
  const Int  iHeight=getHeight(compID);

  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  iWidth * sizeof(Pel),
		  iHeight
		  );
  //for (Int y = iHeight; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, sizeof(Pel)*iWidth);
  //  pDst += iDstStride;
  //  pSrc += iSrcStride;
  //}
}

#if JVET_E0052_DMVR
Void TComYuv::copyToPartXYComponent( const ComponentID compID, const UInt uiAbsZorderIdx, TComYuv* pcYuvDst, const UInt uiPelX, const UInt uiPelY, const Int iWidth, const Int iHeight ) const
{
  const Pel* pSrc     = getAddr(compID, uiAbsZorderIdx);
  Pel* pDst     = pcYuvDst->getAddrPix( compID, uiPelX, uiPelY );

  const UInt iSrcStride  = getStride(compID);
  const UInt iDstStride  = pcYuvDst->getStride(compID);

  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  iWidth * sizeof(Pel),
		  iHeight
		  );
  //for (Int y = iHeight; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, sizeof(Pel)*iWidth);
  //  pDst += iDstStride;
  //  pSrc += iSrcStride;
  //}
}

#if ITH_C2_MC_OPT
Void TComYuv::copyFromPartXYComponent(const ComponentID compID, const UInt uiAbsZorderIdx, TComYuv* pcYuvDst, const UInt uiSrcPelX, const UInt uiSrcPelY, const Int iWidth, const Int iHeight) const
{
  const Pel* pSrc = getAddrPix(compID, uiSrcPelX, uiSrcPelY);
  Pel* pDst = pcYuvDst->getAddr(compID, uiAbsZorderIdx);

  const UInt iSrcStride = getStride(compID);
  const UInt iDstStride = pcYuvDst->getStride(compID);

  for (Int y = iHeight; y != 0; y--)
  {
    ::memcpy(pDst, pSrc, sizeof(Pel)*iWidth);
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}
#endif
#endif

Void TComYuv::copyPartToYuv( TComYuv* pcYuvDst, const UInt uiSrcPartIdx ) const
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    copyPartToComponent  ( ComponentID(comp), pcYuvDst, uiSrcPartIdx );
  }
}

Void TComYuv::copyPartToComponent( const ComponentID compID, TComYuv* pcYuvDst, const UInt uiSrcPartIdx ) const
{
  const Pel* pSrc     = getAddr(compID, uiSrcPartIdx);
        Pel* pDst     = pcYuvDst->getAddr(compID, 0 );

  const UInt  iSrcStride  = getStride(compID);
  const UInt  iDstStride  = pcYuvDst->getStride(compID);

  const UInt uiHeight = pcYuvDst->getHeight(compID);
  const UInt uiWidth = pcYuvDst->getWidth(compID);

  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  uiWidth * sizeof(Pel),
		  uiHeight
		  );
  //for ( UInt y = uiHeight; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, sizeof(Pel)*uiWidth);
  //  pDst += iDstStride;
  //  pSrc += iSrcStride;
  //}
}




Void TComYuv::copyPartToPartYuv   ( TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidth, const UInt iHeight ) const
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    copyPartToPartComponent   (ComponentID(comp), pcYuvDst, uiPartIdx, iWidth>>getComponentScaleX(ComponentID(comp)), iHeight>>getComponentScaleY(ComponentID(comp)) );
  }
}

Void TComYuv::copyPartToPartComponent  ( const ComponentID compID, TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidthComponent, const UInt iHeightComponent ) const
{
  const Pel* pSrc =           getAddr(compID, uiPartIdx);
        Pel* pDst = pcYuvDst->getAddr(compID, uiPartIdx);
  if( pSrc == pDst )
  {
    //th not a good idea
    //th best would be to fix the caller
    return ;
  }

  const UInt  iSrcStride = getStride(compID);
  const UInt  iDstStride = pcYuvDst->getStride(compID);

  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  iWidthComponent * sizeof(Pel),
		  (int)iHeightComponent
		  );
  //for ( UInt y = iHeightComponent; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, iWidthComponent * sizeof(Pel) );
  //  pSrc += iSrcStride;
  //  pDst += iDstStride;
  //}
}




Void TComYuv::copyPartToPartComponentMxN  ( const ComponentID compID, TComYuv* pcYuvDst, const TComRectangle &rect) const
{
  const Pel* pSrc =           getAddrPix( compID, rect.x0, rect.y0 );
        Pel* pDst = pcYuvDst->getAddrPix( compID, rect.x0, rect.y0 );
  if( pSrc == pDst )
  {
    //th not a good idea
    //th best would be to fix the caller
    return ;
  }

  const UInt  iSrcStride = getStride(compID);
  const UInt  iDstStride = pcYuvDst->getStride(compID);
  const UInt uiHeightComponent=rect.height;
  const UInt uiWidthComponent=rect.width;

  twoDimBufCopier
	  (
		  (UChar *)pDst,
		  iDstStride * sizeof(Pel),
		  (UChar *)pSrc,
		  iSrcStride * sizeof(Pel),
		  uiWidthComponent * sizeof(Pel),
		  uiHeightComponent
		  );
  //for ( UInt y = uiHeightComponent; y != 0; y-- )
  //{
  //  ::memcpy( pDst, pSrc, uiWidthComponent * sizeof( Pel ) );
  //  pSrc += iSrcStride;
  //  pDst += iDstStride;
  //}
}

#if RRC_SELECTABLE_IF
Void TComYuv::copyPartToPartComponentDst(const ComponentID compID, TComYuv* pcYuvDst, UInt uiPartIdx, UInt uiPartIdxDst, UInt iWidth, UInt iHeight) const
{
  const Pel* pSrc =           getAddr(compID, uiPartIdx);
        Pel* pDst = pcYuvDst->getAddr(compID, uiPartIdxDst);
  if (pSrc == pDst)
  {
    return;
  }
  UInt  iSrcStride = getStride(compID);
  UInt  iDstStride = pcYuvDst->getStride(compID);
  for (UInt y = iHeight; y != 0; y--)
  {
    ::memcpy(pDst, pSrc, iWidth * sizeof(Pel));
    pSrc += iSrcStride;
    pDst += iDstStride;
  }
}
#endif
#if JVET_C0024_QTBT
Void TComYuv::addClip( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiWidth, const UInt uiHeight, const BitDepths &clipBitDepths )
#else
Void TComYuv::addClip( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiPartSize, const BitDepths &clipBitDepths )
#endif
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID compID=ComponentID(comp);
#if JVET_C0024_QTBT
    const Int uiPartWidth =uiWidth>>getComponentScaleX(compID);
    const Int uiPartHeight=uiHeight>>getComponentScaleY(compID);
#else
    const Int uiPartWidth =uiPartSize>>getComponentScaleX(compID);
    const Int uiPartHeight=uiPartSize>>getComponentScaleY(compID);
#endif

    const Pel* pSrc0 = pcYuvSrc0->getAddr(compID, uiTrUnitIdx, uiPartWidth );
    const Pel* pSrc1 = pcYuvSrc1->getAddr(compID, uiTrUnitIdx, uiPartWidth );
          Pel* pDst  = getAddr(compID, uiTrUnitIdx, uiPartWidth );

    const UInt iSrc0Stride = pcYuvSrc0->getStride(compID);
    const UInt iSrc1Stride = pcYuvSrc1->getStride(compID);
    const UInt iDstStride  = getStride(compID);
#if !JVET_D0033_ADAPTIVE_CLIPPING
    const Int clipbd = clipBitDepths.recon[toChannelType(compID)];
#endif
#if O0043_BEST_EFFORT_DECODING
    const Int bitDepthDelta = clipBitDepths.stream[toChannelType(compID)] - clipbd;
#endif
#if !O0043_BEST_EFFORT_DECODING
	reconComputer
		(
			const_cast<Pel *>(pSrc0),
			iSrc0Stride,
			const_cast<Pel *>(pSrc1),
			iSrc1Stride,
			pDst,
			iDstStride,
#if JVET_D0033_ADAPTIVE_CLIPPING
			g_ClipParam.min(compID),
			g_ClipParam.max(compID),
#else
			0,
			(1 << clipbd) - 1,
#endif
			uiPartWidth,
			uiPartHeight
			);
#else
    for ( Int y = uiPartHeight-1; y >= 0; y-- )
    {
      for ( Int x = uiPartWidth-1; x >= 0; x-- )
      {
#if O0043_BEST_EFFORT_DECODING
        pDst[x] = Pel(ClipBD<Int>( Int(pSrc0[x]) + rightShiftEvenRounding<Pel>(pSrc1[x], bitDepthDelta), clipbd));
#else
#if JVET_D0033_ADAPTIVE_CLIPPING // enc and dec for inter pred+res
        pDst[x] = Pel(ClipA<Int>( Int(pSrc0[x]) + Int(pSrc1[x]), compID));
#else
        pDst[x] = Pel(ClipBD<Int>( Int(pSrc0[x]) + Int(pSrc1[x]), clipbd));
#endif
#endif
      }
      pSrc0 += iSrc0Stride;
      pSrc1 += iSrc1Stride;
      pDst  += iDstStride;
    }
#endif
  }
}




#if JVET_C0024_QTBT
Void TComYuv::subtract( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiWidth, const UInt uiHeight )
#else
Void TComYuv::subtract( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiPartSize )
#endif
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID compID=ComponentID(comp);
#if JVET_C0024_QTBT
    const Int uiPartWidth =uiWidth>>getComponentScaleX(compID);
    const Int uiPartHeight=uiHeight>>getComponentScaleY(compID);
#else
    const Int uiPartWidth =uiPartSize>>getComponentScaleX(compID);
    const Int uiPartHeight=uiPartSize>>getComponentScaleY(compID);
#endif

    const Pel* pSrc0 = pcYuvSrc0->getAddr( compID, uiTrUnitIdx, uiPartWidth );
    const Pel* pSrc1 = pcYuvSrc1->getAddr( compID, uiTrUnitIdx, uiPartWidth );
          Pel* pDst  = getAddr( compID, uiTrUnitIdx, uiPartWidth );

    const Int  iSrc0Stride = pcYuvSrc0->getStride(compID);
    const Int  iSrc1Stride = pcYuvSrc1->getStride(compID);
    const Int  iDstStride  = getStride(compID);

	residueComputer
		(
			const_cast<Pel *>(pSrc0),
			iSrc0Stride,
			const_cast<Pel *>(pSrc1),
			iSrc1Stride,
			pDst,
			iDstStride,
			uiPartWidth,
			uiPartHeight
			);
    //for (Int y = uiPartHeight-1; y >= 0; y-- )
    //{
    //  for (Int x = uiPartWidth-1; x >= 0; x-- )
    //  {
    //    pDst[x] = pSrc0[x] - pSrc1[x];
    //  }
    //  pSrc0 += iSrc0Stride;
    //  pSrc1 += iSrc1Stride;
    //  pDst  += iDstStride;
    //}
  }
}

#if COM16_C806_OBMC
UInt TComYuv::sadLuma( TComYuv* pcYuvSrc0 )
{
  assert( pcYuvSrc0->getWidth( COMPONENT_Y ) == getWidth( COMPONENT_Y ) );
  Pel* pSrc0 = pcYuvSrc0->getAddr( COMPONENT_Y );
  Pel* pSrc1 = getAddr( COMPONENT_Y );

  Int  iSrc0Stride = pcYuvSrc0->getStride( COMPONENT_Y );
  Int  iSrc1Stride = getStride( COMPONENT_Y );
  UInt uiSAD = 0;
  for ( Int y = pcYuvSrc0->getHeight( COMPONENT_Y ) - 1 ; y >= 0; y-- )
  {
    for ( Int x = pcYuvSrc0->getWidth( COMPONENT_Y ) - 1 ; x >= 0; x-- )
    {
      uiSAD += abs( pSrc0[x] - pSrc1[x] );
    }
    pSrc0 += iSrc0Stride;
    pSrc1 += iSrc1Stride;
  }

  return( uiSAD );
}
#endif

#if JVET_E0052_DMVR
Void TComYuv::toLast( const UInt iPartUnitIdx, const UInt uiWidth, const UInt uiHeight, const BitDepths &clipBitDepths )
{
  Pel* p  = getAddr( COMPONENT_Y, iPartUnitIdx );
  const UInt  iSrcStride = getStride(COMPONENT_Y);
  const Int   clipbd      = clipBitDepths.recon[toChannelType(COMPONENT_Y)];
  const Int   shiftNum    = std::max<Int>(2, (IF_INTERNAL_PREC - clipbd));
  const Int   offset      = ( 1 << ( shiftNum - 1 ) ) + IF_INTERNAL_OFFS;

  const Int   iWidth      = uiWidth  >> getComponentScaleX(COMPONENT_Y);
  const Int   iHeight     = uiHeight >> getComponentScaleY(COMPONENT_Y);

  if (iWidth&1)
  {
    assert(0);
    exit(-1);
  }
  else if (iWidth&2)
  {
    for ( Int y = 0; y < iHeight; y++ )
    {
      for (Int x=0 ; x < iWidth; x+=2 )
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
         p[ x + 0 ] = ClipA( rightShift(( p[ x + 0 ] + offset ), shiftNum),  COMPONENT_Y );
         p[ x + 1 ] = ClipA( rightShift(( p[ x + 1 ] + offset ), shiftNum) , COMPONENT_Y);
#else
         p[ x + 0 ] = ClipBD( rightShift(( p[ x + 0 ] + offset ), shiftNum), clipbd );
         p[ x + 1 ] = ClipBD( rightShift(( p[ x + 1 ] + offset ), shiftNum), clipbd );
#endif
      }
       p += iSrcStride;
    }
  }
  else
  {
    for ( Int y = 0; y < iHeight; y++ )
    {
      for (Int x=0 ; x < iWidth; x+=4 )
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
         p[ x + 0 ] = ClipA( rightShift(( p[ x + 0 ] + offset ), shiftNum),  COMPONENT_Y );
         p[ x + 1 ] = ClipA( rightShift(( p[ x + 1 ] + offset ), shiftNum) , COMPONENT_Y);
         p[ x + 2 ] = ClipA( rightShift(( p[ x + 2 ] + offset ), shiftNum),  COMPONENT_Y );
         p[ x + 3 ] = ClipA( rightShift(( p[ x + 3 ] + offset ), shiftNum) , COMPONENT_Y);
#else
         p[ x + 0 ] = ClipBD( rightShift(( p[ x + 0 ] + offset ), shiftNum), clipbd );
         p[ x + 1 ] = ClipBD( rightShift(( p[ x + 1 ] + offset ), shiftNum), clipbd );
         p[ x + 2 ] = ClipBD( rightShift(( p[ x + 2 ] + offset ), shiftNum), clipbd );
         p[ x + 3 ] = ClipBD( rightShift(( p[ x + 3 ] + offset ), shiftNum), clipbd );
#endif
       }
       p += iSrcStride;
     }
  }
}
#endif

Void TComYuv::addAvg( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt iPartUnitIdx, const UInt uiWidth, const UInt uiHeight, const BitDepths &clipBitDepths 
#if VCEG_AZ05_BIO                  
  ,bool bBIOapplied
#endif
#if JVET_E0052_DMVR
  ,Bool bOnlyLuma
#endif
)
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID compID=ComponentID(comp);
#if JVET_E0052_DMVR
    if (bOnlyLuma)
    {
      if (compID!=COMPONENT_Y)
      {
        continue;
      }
    }
#endif
#if VCEG_AZ05_BIO
    if(!bBIOapplied || compID!=COMPONENT_Y)
    {
#endif 
    const Pel* pSrc0  = pcYuvSrc0->getAddr( compID, iPartUnitIdx );
    const Pel* pSrc1  = pcYuvSrc1->getAddr( compID, iPartUnitIdx );
    Pel* pDst   = getAddr( compID, iPartUnitIdx );

    const UInt  iSrc0Stride = pcYuvSrc0->getStride(compID);
    const UInt  iSrc1Stride = pcYuvSrc1->getStride(compID);
    const UInt  iDstStride  = getStride(compID);
    const Int   clipbd      = clipBitDepths.recon[toChannelType(compID)];
    const Int   shiftNum    = std::max<Int>(2, (IF_INTERNAL_PREC - clipbd)) + 1;
    const Int   offset      = ( 1 << ( shiftNum - 1 ) ) + 2 * IF_INTERNAL_OFFS;

    const Int   iWidth      = uiWidth  >> getComponentScaleX(compID);
    const Int   iHeight     = uiHeight >> getComponentScaleY(compID);

    if (iWidth&1)
    {
      assert(0);
      exit(-1);
    }
    else if (iWidth&2)
    {
      for ( Int y = 0; y < iHeight; y++ )
      {
        for (Int x=0 ; x < iWidth; x+=2 )
        {
#if JVET_D0033_ADAPTIVE_CLIPPING
            pDst[ x + 0 ] = ClipA( rightShift(( pSrc0[ x + 0 ] + pSrc1[ x + 0 ] + offset ), shiftNum),  compID );
            pDst[ x + 1 ] = ClipA( rightShift(( pSrc0[ x + 1 ] + pSrc1[ x + 1 ] + offset ), shiftNum) , compID);

#else
          pDst[ x + 0 ] = ClipBD( rightShift(( pSrc0[ x + 0 ] + pSrc1[ x + 0 ] + offset ), shiftNum), clipbd );
          pDst[ x + 1 ] = ClipBD( rightShift(( pSrc0[ x + 1 ] + pSrc1[ x + 1 ] + offset ), shiftNum), clipbd );
#endif
        }
        pSrc0 += iSrc0Stride;
        pSrc1 += iSrc1Stride;
        pDst  += iDstStride;
      }
    }
    else
    {
      for ( Int y = 0; y < iHeight; y++ )
      {
        for (Int x=0 ; x < iWidth; x+=4 )
        {
#if JVET_D0033_ADAPTIVE_CLIPPING
                pDst[ x + 0 ] = ClipA( rightShift(( pSrc0[ x + 0 ] + pSrc1[ x + 0 ] + offset ), shiftNum),  compID );
                pDst[ x + 1 ] = ClipA( rightShift(( pSrc0[ x + 1 ] + pSrc1[ x + 1 ] + offset ), shiftNum),  compID );
                pDst[ x + 2 ] = ClipA( rightShift(( pSrc0[ x + 2 ] + pSrc1[ x + 2 ] + offset ), shiftNum),  compID );
                pDst[ x + 3 ] = ClipA( rightShift(( pSrc0[ x + 3 ] + pSrc1[ x + 3 ] + offset ), shiftNum),  compID );
#else
          pDst[ x + 0 ] = ClipBD( rightShift(( pSrc0[ x + 0 ] + pSrc1[ x + 0 ] + offset ), shiftNum), clipbd );
          pDst[ x + 1 ] = ClipBD( rightShift(( pSrc0[ x + 1 ] + pSrc1[ x + 1 ] + offset ), shiftNum), clipbd );
          pDst[ x + 2 ] = ClipBD( rightShift(( pSrc0[ x + 2 ] + pSrc1[ x + 2 ] + offset ), shiftNum), clipbd );
          pDst[ x + 3 ] = ClipBD( rightShift(( pSrc0[ x + 3 ] + pSrc1[ x + 3 ] + offset ), shiftNum), clipbd );
#endif
        }
        pSrc0 += iSrc0Stride;
        pSrc1 += iSrc1Stride;
        pDst  += iDstStride;
      }
    }
#if VCEG_AZ05_BIO
    }
#endif 
  }
}

#if ITH_C2_UNCLIPPED_DMVR
Void averagerWithRoundingUnclipped
(
  Pel *pSrc0,
  Int iSrc0Stride,
  Pel *pSrc1,
  Int iSrc1Stride,
  Pel *pDst,
  Int iDstStride,
  const Int iOffset,
  const Int iShift,
  const Int iMin,
  const Int iMax,
  const Int iWidth,
  const Int iHeight
)
{
  Int iIntermediary;

  for (Int j = 0; j < iHeight; j++)
  {
    for (Int i = 0; i < iWidth; i++)
    {
      /*iIntermediary = pSrc0[i] + pSrc1[i] + iOffset;
      pDst[i + j * iDstStride] =(iIntermediary >> iShift);*/

      iIntermediary = pSrc0[i] + pSrc1[i] + 1;
      pDst[i + j * iDstStride] = (iIntermediary >> 1);
    }

    pSrc0 += iSrc0Stride;
    pSrc1 += iSrc1Stride;
  }
}

Void TComYuv::addAvgUnclipped(const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt iPartUnitIdx, const UInt uiWidth, const UInt uiHeight, const BitDepths &clipBitDepths
#if VCEG_AZ05_BIO
  , bool bBIOapplied
#endif
#if JVET_E0052_DMVR
  , Bool bOnlyLuma
#endif
)
{
  for (Int comp = 0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID compID = ComponentID(comp);
#if JVET_E0052_DMVR
    if (bOnlyLuma)
    {
      if (compID != COMPONENT_Y)
      {
        continue;
      }
    }
#endif
#if VCEG_AZ05_BIO
    if (!bBIOapplied || compID != COMPONENT_Y)
    {
#endif 
      const Pel* pSrc0 = pcYuvSrc0->getAddr(compID, iPartUnitIdx);
      const Pel* pSrc1 = pcYuvSrc1->getAddr(compID, iPartUnitIdx);
      Pel* pDst = getAddr(compID, iPartUnitIdx);

      const UInt  iSrc0Stride = pcYuvSrc0->getStride(compID);
      const UInt  iSrc1Stride = pcYuvSrc1->getStride(compID);
      const UInt  iDstStride = getStride(compID);
      const Int   clipbd = clipBitDepths.recon[toChannelType(compID)];
      const Int   shiftNum = std::max<Int>(2, (IF_INTERNAL_PREC - clipbd)) + 1;
      const Int   offset = (1 << (shiftNum - 1)) + 2 * IF_INTERNAL_OFFS;

      const Int   iWidth = uiWidth >> getComponentScaleX(compID);
      const Int   iHeight = uiHeight >> getComponentScaleY(compID);

      if (iWidth & 1)
      {
        assert(0);
        exit(-1);
      }

      {
        Int iMin, iMax;

#if JVET_D0033_ADAPTIVE_CLIPPING
        iMin = g_ClipParam.min(compID);
        iMax = g_ClipParam.max(compID);
#else
        iMin = 0;
        iMax = (1 << clipbd) - 1;
#endif

        averagerWithRoundingUnclipped
        (
          const_cast<Pel *>(pSrc0),
          iSrc0Stride,
          const_cast<Pel *>(pSrc1),
          iSrc1Stride,
          pDst,
          iDstStride,
          offset,
          shiftNum,
          iMin,
          iMax,
          iWidth,
          iHeight
        );
      }
#if VCEG_AZ05_BIO
    }
#endif 
  }
}
#endif

Void TComYuv::removeHighFreq( const TComYuv* pcYuvSrc,
                              const UInt uiPartIdx,
                              const UInt uiWidth,
                              const UInt uiHeight,
                              const Int bitDepths[MAX_NUM_CHANNEL_TYPE],
                              const Bool bClipToBitDepths
                              )
{
  for(Int comp=0; comp<getNumberValidComponents(); comp++)
  {
    const ComponentID compID=ComponentID(comp);
    const Pel* pSrc  = pcYuvSrc->getAddr(compID, uiPartIdx);
    Pel* pDst  = getAddr(compID, uiPartIdx);

    const Int iSrcStride = pcYuvSrc->getStride(compID);
    const Int iDstStride = getStride(compID);
    const Int iWidth  = uiWidth >>getComponentScaleX(compID);
    const Int iHeight = uiHeight>>getComponentScaleY(compID);
    if (bClipToBitDepths)
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      const Int clipBd=bitDepths[toChannelType(compID)];
#endif
      for ( Int y = iHeight-1; y >= 0; y-- )
      {
        for ( Int x = iWidth-1; x >= 0; x-- )
        {
#if JVET_D0033_ADAPTIVE_CLIPPING
          pDst[x ] = ClipA((2 * pDst[x]) - pSrc[x],  compID);
#else
          pDst[x ] = ClipBD((2 * pDst[x]) - pSrc[x], clipBd);
#endif
        }
        pSrc += iSrcStride;
        pDst += iDstStride;
      }
    }
    else
    {
      for ( Int y = iHeight-1; y >= 0; y-- )
      {
        for ( Int x = iWidth-1; x >= 0; x-- )
        {
          pDst[x ] = (2 * pDst[x]) - pSrc[x];
        }
        pSrc += iSrcStride;
        pDst += iDstStride;
      }
    }
  }
}

#if VCEG_AZ08_KLT_COMMON
Void TComPicYuv::fillPicRecBoundary(const BitDepths bitDepths)
{
  const ComponentID compID = COMPONENT_Y;
  Pel* piTxt = getAddr(compID);
  Int iStride = getStride(compID);
  Int iWidth = getWidth(compID);
  Int iHeight = getHeight(compID);
  Int iMarginX = getMarginX(compID);
  Int iMarginY = getMarginY(compID);
  const Int channelBitDepth = bitDepths.recon[toChannelType(compID)];
  Int   x, y;
  Pel*  pi;
  Pel fillValue = 1 << (channelBitDepth - 1);
  pi = piTxt;
  for (y = 0; y < iHeight; y++)
  {
    for (x = 0; x < iMarginX; x++)
    {
      pi[-iMarginX + x] = fillValue;
      pi[iWidth + x] = fillValue;
    }
    pi += iStride;
  }
  
  pi = piTxt - (iMarginY*iStride + iMarginX);
  for (y = 0; y < iMarginY; y++)
  {
      for (x = 0; x < iStride; x++)
      {
          pi[x] = fillValue;
      }
      pi += iStride;
  }
}
#endif
//! \}

