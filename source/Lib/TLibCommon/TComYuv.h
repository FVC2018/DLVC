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

/** \file     TComYuv.h
    \brief    general YUV buffer class (header)
    \todo     this should be merged with TComPicYuv \n
              check usage of removeHighFreq function
*/

#ifndef __TCOMYUV__
#define __TCOMYUV__
#include "CommonDef.h"
#include "TComPicYuv.h"
#include "TComRectangle.h"

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// general YUV buffer class
class TComYuv
{
private:

  // ------------------------------------------------------------------------------------------------------------------
  //  YUV buffer
  // ------------------------------------------------------------------------------------------------------------------

  Pel*    m_apiBuf[MAX_NUM_COMPONENT];

  // ------------------------------------------------------------------------------------------------------------------
  //  Parameter for general YUV buffer usage
  // ------------------------------------------------------------------------------------------------------------------

  UInt     m_iWidth;
  UInt     m_iHeight;
  ChromaFormat m_chromaFormatIDC; ////< Chroma Format

  // dims 16x16
  // blkSize=4x4

  // these functions assume a square CU, of size width*width, split into square TUs each of size blkSize*blkSize.
  // iTransUnitIdx is the raster-scanned index of the sub-block (TU) in question.
  // eg for a 16x16 CU, with 4x4 TUs:
  //  0  1  2  3
  //  4  5  6  7
  //  8  9 10 11
  // 12 13 14 15

  // So, for iTransUnitIdx=14, 14*4 & 15 =8=X offset.
  //                           14*4 / 16 =3=Y block offset
  //                                      3*4*16 = Y offset within buffer


public:

               TComYuv                    ();
  virtual     ~TComYuv                    ();

  // ------------------------------------------------------------------------------------------------------------------
  //  Memory management
  // ------------------------------------------------------------------------------------------------------------------

  Void         create                     ( const UInt iWidth, const UInt iHeight, const ChromaFormat chromaFormatIDC );  ///< Create  YUV buffer
  Void         destroy                    ();                             ///< Destroy YUV buffer
  Void         clear                      ();                             ///< clear   YUV buffer

  // ------------------------------------------------------------------------------------------------------------------
  //  Copy, load, store YUV buffer
  // ------------------------------------------------------------------------------------------------------------------

  //  Copy YUV buffer to picture buffer
  Void         copyToPicYuv               ( TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth = 0, const UInt uiPartIdx = 0
#if BARC
    ,const Bool isBARC = false
#endif
  ) const ;
  Void         copyToPicComponent         ( const ComponentID id, TComPicYuv* pcPicYuvDst, const UInt iCtuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth = 0, const UInt uiPartIdx = 0 
#if BARC
    , const Bool isBARC = false
#endif
  ) const ;

  //  Copy YUV buffer from picture buffer
  Void         copyFromPicYuv(const TComPicYuv* pcPicYuvSrc, const  UInt ctuRsAddr, const UInt uiAbsZorderIdx
#if BARC
    ,const Bool isBARC = false
#endif
  );
  Void         copyFromPicComponent       ( const ComponentID id, const TComPicYuv* pcPicYuvSrc, const UInt iCtuRsAddr, const UInt uiAbsZorderIdx 
#if BARC
    , const Bool isBARC = false
#endif
  );

#if JVET_D0033_ADAPTIVE_CLIPPING
  Void clipToPartYuv(TComYuv* pcYuvDst, const UInt uiDstPartIdx, const BitDepths &bitDepths) const;
  Void clipToPartComponent(const ComponentID compID, TComYuv* pcYuvDst, const UInt uiDstPartIdx) const;
  Void clipToPicYuv   ( TComPicYuv* pcPicYuvDst, const UInt ctuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth = 0, const UInt uiPartIdx = 0 ) const ;
  Void clipToPicComponent   ( const ComponentID id, TComPicYuv* pcPicYuvDst, const UInt iCtuRsAddr, const UInt uiAbsZorderIdx, const UInt uiPartDepth = 0, const UInt uiPartIdx = 0 ) const ;
  Void clipPartToPartYuv(TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidth, const UInt iHeight) const;
  Void clipPartToPartComponent(const ComponentID compID, TComYuv* pcYuvDst, const UInt uiPartIdx, const UInt iWidthComponent, const UInt iHeightComponent) const;
#endif
  //  Copy Small YUV buffer to the part of other Big YUV buffer
  Void         copyToPartYuv              ( TComYuv*    pcYuvDst,    const UInt uiDstPartIdx ) const ;
  Void         copyToPartComponent        ( const ComponentID id, TComYuv*    pcYuvDst,    const UInt uiDstPartIdx ) const ;
#if JVET_E0052_DMVR
  Void         copyToPartXYComponent      ( const ComponentID compID, const UInt uiAbsZorderIdx, TComYuv* pcYuvDst, const UInt uiPelX, const UInt uiPelY, const Int iWidth, const Int iHeight ) const;

#if ITH_C2_MC_OPT
  Void         copyFromPartXYComponent    ( const ComponentID compID, const UInt uiAbsZorderIdx, TComYuv* pcYuvDst, const UInt uiSrcPelX, const UInt uiSrcPelY, const Int iWidth, const Int iHeight ) const;
#endif
#endif
  //  Copy the part of Big YUV buffer to other Small YUV buffer
  Void         copyPartToYuv              ( TComYuv*    pcYuvDst,   const UInt uiSrcPartIdx ) const;
  Void         copyPartToComponent        ( const ComponentID id, TComYuv*    pcYuvDst,    const UInt uiSrcPartIdx ) const;

  //  Copy YUV partition buffer to other YUV partition buffer
  Void         copyPartToPartYuv          ( TComYuv*    pcYuvDst, const UInt uiPartIdx, const UInt uiWidth, const UInt uiHeight ) const;
  Void         copyPartToPartComponent    ( const ComponentID id, TComYuv*    pcYuvDst, const UInt uiPartIdx, const UInt uiWidthComponent, const UInt uiHeightComponent ) const;

 // Copy YUV partition buffer to other YUV partition buffer for non-square blocks
  Void         copyPartToPartComponentMxN ( const ComponentID id, TComYuv*    pcYuvDst, const TComRectangle &rect ) const;

#if RRC_SELECTABLE_IF
  Void    copyPartToPartComponentDst      ( const ComponentID id, TComYuv*    pcYuvDst, UInt uiPartIdx, UInt uiPartIdxDst, UInt uiWidth, UInt uiHeight ) const;
#endif
  // ------------------------------------------------------------------------------------------------------------------
  //  Algebraic operation for YUV buffer
  // ------------------------------------------------------------------------------------------------------------------

  //  Clip(pcYuvSrc0 + pcYuvSrc1) -> m_apiBuf
#if JVET_C0024_QTBT
  Void         addClip                    ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiWidth, const UInt uiHeight, const BitDepths &clipBitDepths );
#else
  Void         addClip                    ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiPartSize, const BitDepths &clipBitDepths );
#endif

  //  pcYuvSrc0 - pcYuvSrc1 -> m_apiBuf
#if JVET_C0024_QTBT
  Void         subtract                   ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiWidth, const UInt uiHeight );
#else
  Void         subtract                   ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt uiTrUnitIdx, const UInt uiPartSize );
#endif
#if COM16_C806_OBMC
  UInt         sadLuma                    ( TComYuv* pcYuvSrc0 );
#endif
  //  (pcYuvSrc0 + pcYuvSrc1)/2 for YUV partition
  Void         addAvg                     ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt iPartUnitIdx, const UInt iWidth, const UInt iHeight, const BitDepths &clipBitDepths 
#if VCEG_AZ05_BIO                  
    , bool bBIOapplied
#endif
#if JVET_E0052_DMVR
   , Bool bOnlyLuma=false
#endif
);

#if ITH_C2_UNCLIPPED_DMVR
  Void         addAvgUnclipped            ( const TComYuv* pcYuvSrc0, const TComYuv* pcYuvSrc1, const UInt iPartUnitIdx, const UInt iWidth, const UInt iHeight, const BitDepths &clipBitDepths 
#if VCEG_AZ05_BIO
    , bool bBIOapplied
#endif
#if JVET_E0052_DMVR
    , Bool bOnlyLuma=false
#endif
);
#endif

#if JVET_E0052_DMVR
  Void         toLast                     ( const UInt iPartUnitIdx, const UInt iWidth, const UInt iHeight, const BitDepths &clipBitDepths );
#endif
  Void         removeHighFreq             ( const TComYuv* pcYuvSrc, const UInt uiPartIdx, const UInt uiWidth, const UInt uiHeight
                                          , const Int bitDepths[MAX_NUM_CHANNEL_TYPE], const Bool bClipToBitDepths
                                          );

  // ------------------------------------------------------------------------------------------------------------------
  //  Access function for YUV buffer
  // ------------------------------------------------------------------------------------------------------------------

  //  Access starting position of YUV buffer
  Pel*         getAddr                    (const ComponentID id)                    { return m_apiBuf[id]; }
  const Pel*   getAddr                    (const ComponentID id) const              { return m_apiBuf[id]; }

  //  Access starting position of YUV partition unit buffer
  Pel*         getAddr                    (const ComponentID id, const UInt uiPartUnitIdx)
                                              {
                                                  Int blkX = g_auiRasterToPelX[ g_auiZscanToRaster[ uiPartUnitIdx ] ] >> getComponentScaleX(id);
                                                  Int blkY = g_auiRasterToPelY[ g_auiZscanToRaster[ uiPartUnitIdx ] ] >> getComponentScaleY(id);
                                                  assert((blkX<getWidth(id) && blkY<getHeight(id)));
                                                  return m_apiBuf[id] + blkX + blkY * getStride(id);
                                              }
  const Pel*   getAddr                    (const ComponentID id, const UInt uiPartUnitIdx) const
                                              {
                                                  Int blkX = g_auiRasterToPelX[ g_auiZscanToRaster[ uiPartUnitIdx ] ] >> getComponentScaleX(id);
                                                  Int blkY = g_auiRasterToPelY[ g_auiZscanToRaster[ uiPartUnitIdx ] ] >> getComponentScaleY(id);
                                                  assert((blkX<getWidth(id) && blkY<getHeight(id)));
                                                  return m_apiBuf[id] + blkX + blkY * getStride(id);
                                              }

  //  Access starting position of YUV transform unit buffer
  Pel*         getAddr                    (const ComponentID id, const UInt iTransUnitIdx, const UInt iBlkSizeForComponent)
                                              {
                                                UInt width=getWidth(id);
                                                Int blkX = ( iTransUnitIdx * iBlkSizeForComponent ) &  ( width - 1 );
                                                Int blkY = ( iTransUnitIdx * iBlkSizeForComponent ) &~ ( width - 1 );
                                                if (m_chromaFormatIDC==CHROMA_422 && id!=COMPONENT_Y)
                                                {
                                                  blkY<<=1;
                                                }
                                                return m_apiBuf[id] + blkX + blkY * iBlkSizeForComponent;
                                              }

  const Pel*   getAddr                    (const ComponentID id, const UInt iTransUnitIdx, const UInt iBlkSizeForComponent) const
                                              {
                                                UInt width=getWidth(id);
                                                Int blkX = ( iTransUnitIdx * iBlkSizeForComponent ) &  ( width - 1 );
                                                Int blkY = ( iTransUnitIdx * iBlkSizeForComponent ) &~ ( width - 1 );
                                                if (m_chromaFormatIDC==CHROMA_422 && id!=COMPONENT_Y)
                                                {
                                                  blkY<<=1;
                                                }
                                                return m_apiBuf[id] + blkX + blkY * iBlkSizeForComponent;
                                              }

  // Access starting position of YUV transform unit buffer by pix offset for square & non-square blocks
  Pel*         getAddrPix                 (const ComponentID id, const UInt iPixX, const UInt iPixY )       { return m_apiBuf[id] + iPixY * getStride(id) + iPixX; }
  const Pel*   getAddrPix                 (const ComponentID id, const UInt iPixX, const UInt iPixY ) const { return m_apiBuf[id] + iPixY * getStride(id) + iPixX; }

  //  Get stride value of YUV buffer
  UInt         getStride                  (const ComponentID id) const { return m_iWidth >> getComponentScaleX(id);   }
  UInt         getHeight                  (const ComponentID id) const { return m_iHeight >> getComponentScaleY(id);  }
  UInt         getWidth                   (const ComponentID id) const { return m_iWidth >> getComponentScaleX(id);   }
  ChromaFormat getChromaFormat            ()                     const { return m_chromaFormatIDC; }
  UInt         getNumberValidComponents   ()                     const { return ::getNumberValidComponents(m_chromaFormatIDC); }
  UInt         getComponentScaleX         (const ComponentID id) const { return ::getComponentScaleX(id, m_chromaFormatIDC); }
  UInt         getComponentScaleY         (const ComponentID id) const { return ::getComponentScaleY(id, m_chromaFormatIDC); }
#if JVET_C0024_FAST_MRG
  Void         setWidth                   (UInt uiWidth) { m_iWidth = uiWidth;}
  Void         setHeight                  (UInt uiHeight) { m_iHeight = uiHeight;}
#endif
};// END CLASS DEFINITION TComYuv

//! \}

#endif // __TCOMYUV__
