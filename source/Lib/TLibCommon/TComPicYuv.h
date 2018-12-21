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

/** \file     TComPicYuv.h
    \brief    picture YUV buffer class (header)
*/

#ifndef __TCOMPICYUV__
#define __TCOMPICYUV__

#include <stdio.h>
#include "CommonDef.h"
#include "TComRom.h"
#include "TComChromaFormat.h"
#include "SEI.h"

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// picture YUV buffer class
class TComPicYuv
{
private:

  // ------------------------------------------------------------------------------------------------
  //  YUV buffer
  // ------------------------------------------------------------------------------------------------

  Pel*  m_apiPicBuf[MAX_NUM_COMPONENT];             ///< Buffer (including margin)

  Pel*  m_piPicOrg[MAX_NUM_COMPONENT];              ///< m_apiPicBufY + m_iMarginLuma*getStride() + m_iMarginLuma

  // ------------------------------------------------------------------------------------------------
  //  Parameter for general YUV buffer usage
  // ------------------------------------------------------------------------------------------------

  Int   m_iPicWidth;                                ///< Width of picture in pixels
  Int   m_iPicHeight;                               ///< Height of picture in pixels
  ChromaFormat m_chromaFormatIDC;                   ///< Chroma Format

  Int*  m_ctuOffsetInBuffer[MAX_NUM_CHANNEL_TYPE];  ///< Gives an offset in the buffer for a given CTU (and channel)
  Int*  m_subCuOffsetInBuffer[MAX_NUM_CHANNEL_TYPE];///< Gives an offset in the buffer for a given sub-CU (and channel), relative to start of CTU
#if BARC
  Int*  m_ctuOffsetInBufferDscale[MAX_NUM_CHANNEL_TYPE];  ///< Gives an offset in the buffer for a given CTU (and channel)
  Int*  m_subCuOffsetInBufferDscale[MAX_NUM_CHANNEL_TYPE];///< Gives an offset in the buffer for a given sub-CU (and channel), relative to start of CTU
#endif
  Int   m_iMarginX;                                 ///< margin of Luma channel (chroma's may be smaller, depending on ratio)
  Int   m_iMarginY;                                 ///< margin of Luma channel (chroma's may be smaller, depending on ratio)

  Bool  m_bIsBorderExtended;

public:
               TComPicYuv         ();
  virtual     ~TComPicYuv         ();


  // ------------------------------------------------------------------------------------------------
  //  Memory management
  // ------------------------------------------------------------------------------------------------

  Void          create            (const Int iPicWidth,
                                   const Int iPicHeight,
                                   const ChromaFormat chromaFormatIDC,
                                   const UInt uiMaxCUWidth,  ///< used for generating offsets to CUs. Can use iPicWidth if no offsets are required
                                   const UInt uiMaxCUHeight, ///< used for generating offsets to CUs. Can use iPicHeight if no offsets are required
                                   const UInt uiMaxCUDepth,  ///< used for generating offsets to CUs. Can use 0 if no offsets are required
                                   const Bool bUseMargin);   ///< if true, then a margin of uiMaxCUWidth+16 and uiMaxCUHeight+16 is created around the image.

  Void          destroy           ();

#if BLOCK_BASED_BACKGROUND_REFERENCE
	Void  clear();
#endif

  // The following have been removed - Use CHROMA_400 in the above function call.
  //Void  createLuma  ( Int iPicWidth, Int iPicHeight, UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uhMaxCUDepth );
  //Void  destroyLuma ();

  // ------------------------------------------------------------------------------------------------
  //  Get information of picture
  // ------------------------------------------------------------------------------------------------

  Int           getWidth          (const ComponentID id) const { return  m_iPicWidth >> getComponentScaleX(id);   }
  Int           getHeight         (const ComponentID id) const { return  m_iPicHeight >> getComponentScaleY(id);  }
  ChromaFormat  getChromaFormat   ()                     const { return m_chromaFormatIDC; }
  UInt          getNumberValidComponents() const { return ::getNumberValidComponents(m_chromaFormatIDC); }

  Int           getStride         (const ComponentID id) const { return ((m_iPicWidth     ) + (m_iMarginX  <<1)) >> getComponentScaleX(id); }
  Int           getTotalHeight    (const ComponentID id) const { return ((m_iPicHeight    ) + (m_iMarginY  <<1)) >> getComponentScaleY(id); }

  Int           getMarginX        (const ComponentID id) const { return m_iMarginX >> getComponentScaleX(id);  }
  Int           getMarginY        (const ComponentID id) const { return m_iMarginY >> getComponentScaleY(id);  }

  // ------------------------------------------------------------------------------------------------
  //  Access function for picture buffer
  // ------------------------------------------------------------------------------------------------

  //  Access starting position of picture buffer with margin
  Pel*          getBuf            (const ComponentID ch)       { return  m_apiPicBuf[ch];   }

  //  Access starting position of original picture
  Pel*          getAddr           (const ComponentID ch)       { return  m_piPicOrg[ch];   }
  const Pel*    getAddr           (const ComponentID ch) const { return  m_piPicOrg[ch];   }

  //  Access starting position of original picture for specific coding unit (CU) or partition unit (PU)
  Pel*          getAddr           (const ComponentID ch, const Int ctuRSAddr )       { return m_piPicOrg[ch] + m_ctuOffsetInBuffer[ch==0?0:1][ ctuRSAddr ]; }
  const Pel*    getAddr           (const ComponentID ch, const Int ctuRSAddr ) const { return m_piPicOrg[ch] + m_ctuOffsetInBuffer[ch==0?0:1][ ctuRSAddr ]; }
  Pel*          getAddr           (const ComponentID ch, const Int ctuRSAddr, const Int uiAbsZorderIdx )
                                     { return m_piPicOrg[ch] + m_ctuOffsetInBuffer[ch==0?0:1][ctuRSAddr] + m_subCuOffsetInBuffer[ch==0?0:1][g_auiZscanToRaster[uiAbsZorderIdx]]; }
  const Pel*    getAddr           (const ComponentID ch, const Int ctuRSAddr, const Int uiAbsZorderIdx ) const
                                     { return m_piPicOrg[ch] + m_ctuOffsetInBuffer[ch==0?0:1][ctuRSAddr] + m_subCuOffsetInBuffer[ch==0?0:1][g_auiZscanToRaster[uiAbsZorderIdx]]; }

#if BARC_FUNC
  //  Access starting position of original picture for specific coding unit (CU) or partition unit (PU)
  Pel*          getAddrDscale           (const ComponentID ch, const Int ctuRSAddr, const Int uiAbsZorderIdx )
                                     { return m_piPicOrg[ch] + m_ctuOffsetInBufferDscale[ch==0?0:1][ctuRSAddr] + m_subCuOffsetInBuffer[ch==0?0:1][g_auiZscanToRaster[uiAbsZorderIdx]]; }
  const Pel*    getAddrDscale           (const ComponentID ch, const Int ctuRSAddr, const Int uiAbsZorderIdx ) const
                                     { return m_piPicOrg[ch] + m_ctuOffsetInBufferDscale[ch==0?0:1][ctuRSAddr] + m_subCuOffsetInBuffer[ch==0?0:1][g_auiZscanToRaster[uiAbsZorderIdx]]; }
#endif
  UInt          getComponentScaleX(const ComponentID id) const { return ::getComponentScaleX(id, m_chromaFormatIDC); }
  UInt          getComponentScaleY(const ComponentID id) const { return ::getComponentScaleY(id, m_chromaFormatIDC); }

  // ------------------------------------------------------------------------------------------------
  //  Miscellaneous
  // ------------------------------------------------------------------------------------------------

  //  Copy function to picture
  Void          copyToPic         ( TComPicYuv*  pcPicYuvDst 
#if ALF_HM3_REFACTOR
    , ComponentID chDst = MAX_NUM_COMPONENT
    , Bool bMarginIncluded = true
#endif
    ) const ;

  //  Extend function of picture buffer
  Void          extendPicBorder   (
#if ALF_HM3_REFACTOR
    Int nMargin = -1      // use default margin
#endif
    );
#if BARC_FUNC
  Void  extendPicBorderBARC  (
#if ALF_HM3_REFACTOR
    Int nMargin = -1      // use default margin
#endif
);
#endif
#if VCEG_AZ08_KLT_COMMON
  Void          fillPicRecBoundary(const BitDepths bitDepths);
#endif
  //  Dump picture
  Void          dump              (const Char* pFileName, const BitDepths &bitDepths, Bool bAdd = false) const ;

  // Set border extension flag
  Void          setBorderExtension(Bool b) { m_bIsBorderExtended = b; }

#if EXTENSION_360_VIDEO
  Void          createWithoutCUInfo(const Int picWidth,
                                    const Int picHeight,
                                    const ChromaFormat chromaFormatIDC,
                                    const Bool bUseMargin=false, ///< if true, then a margin of uiMaxCUWidth+16 and uiMaxCUHeight+16 is created around the image.
                                    const UInt maxCUWidth=0,   ///< used for margin only
                                    const UInt maxCUHeight=0); ///< used for margin only
#endif
};// END CLASS DEFINITION TComPicYuv


// These functions now return the length of the digest strings.
UInt calcChecksum(const TComPicYuv& pic, TComPictureHash &digest, const BitDepths &bitDepths);
UInt calcCRC     (const TComPicYuv& pic, TComPictureHash &digest, const BitDepths &bitDepths);
UInt calcMD5     (const TComPicYuv& pic, TComPictureHash &digest, const BitDepths &bitDepths);
std::string hashToString(const TComPictureHash &digest, Int numChar);
//! \}

#endif // __TCOMPICYUV__