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

/** \file     TComPattern.cpp
    \brief    neighbouring pixel access classes
*/

#include "TComPic.h"
#include "TComPattern.h"
#include "TComDataCU.h"
#include "TComTU.h"
#include "Debug.h"
#include "TComPrediction.h"


//! \ingroup TLibCommon
//! \{

// Forward declarations

/// padding of unavailable reference samples for intra prediction
Void fillReferenceSamples( const Int bitDepth, 
#if O0043_BEST_EFFORT_DECODING
                           const Int bitDepthDelta, 
#endif
                           const Pel* piRoiOrigin, 
                                 Pel* piIntraTemp,
                           const Bool* bNeighborFlags,
                           const Int iNumIntraNeighbor, 
                           const Int unitWidth, 
                           const Int unitHeight, 
                           const Int iAboveUnits, 
                           const Int iLeftUnits,
                           const UInt uiWidth, 
                           const UInt uiHeight, 
                           const Int iPicStride 
#if MLIP_EXTEND_REF_LINE
    , UInt ExternRefNum
    , Pel* piSecondIntraTemp
#endif
);

/// constrained intra prediction
Bool  isAboveLeftAvailable  ( TComDataCU* pcCU, UInt uiPartIdxLT );
Int   isAboveAvailable      ( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxRT, Bool* bValidFlags );
Int   isLeftAvailable       ( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxLB, Bool* bValidFlags );
#if JVET_C0024_QTBT
Int   isAboveRightAvailable ( TComDataCU* pcCU, UInt uiPartIdxRT, UInt uiNumUnitsInPU, Bool *bValidFlags );
Int   isBelowLeftAvailable  ( TComDataCU* pcCU, UInt uiPartIdxLB, UInt uiNumUnitsInPU, Bool* bValidFlags );
#else
Int   isAboveRightAvailable ( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxRT, Bool* bValidFlags );
Int   isBelowLeftAvailable  ( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxLB, Bool* bValidFlags );
#endif


// ====================================================================================================================
// Public member functions (TComPatternParam)
// ====================================================================================================================

/** 
 \param  piTexture     pixel data
 \param  iRoiWidth     pattern width
 \param  iRoiHeight    pattern height
 \param  iStride       buffer stride
 \param  bitDepth      bit depth
 */
Void TComPatternParam::setPatternParamPel ( Pel* piTexture,
                                           Int iRoiWidth,
                                           Int iRoiHeight,
                                           Int iStride,
                                           Int bitDepth
                                           )
{
  m_piROIOrigin    = piTexture;
  m_iROIWidth      = iRoiWidth;
  m_iROIHeight     = iRoiHeight;
  m_iPatternStride = iStride;
  m_bitDepth       = bitDepth;
}


// ====================================================================================================================
// Initialize static members (TComPattern)
// ====================================================================================================================


// ====================================================================================================================
// Public member functions (TComPattern)
// ====================================================================================================================

Void TComPattern::initPattern (Pel* piY,
                               Int iRoiWidth,
                               Int iRoiHeight,
                               Int iStride,
                               Int bitDepthLuma)
{
  m_cPatternY. setPatternParamPel( piY,  iRoiWidth, iRoiHeight, iStride, bitDepthLuma);
}



// TODO: move this function to TComPrediction.cpp.
Void TComPrediction::initIntraPatternChType( TComTU &rTu, const ComponentID compID, const Bool bFilterRefSamples 
#if RRC_INTRA_SHARPENING
   , Bool bSharpRefSamples
#endif
#if BARC
  , Bool isCTUBoundary
#endif
   DEBUG_STRING_FN_DECLARE(sDebug)
  )
{
  const ChannelType chType    = toChannelType(compID);

  TComDataCU *pcCU=rTu.getCU();
  const TComSPS &sps = *(pcCU->getSlice()->getSPS());
  const UInt uiZorderIdxInPart=rTu.GetAbsPartIdxTU();
  const UInt uiTuWidth        = rTu.getRect(compID).width;
  const UInt uiTuHeight       = rTu.getRect(compID).height;
#if JVET_C0024_QTBT
  const UInt uiTuWidth2       = uiTuWidth  + uiTuHeight;
  const UInt uiTuHeight2      = uiTuHeight + uiTuWidth;
#else
  const UInt uiTuWidth2       = uiTuWidth  << 1;
  const UInt uiTuHeight2      = uiTuHeight << 1;
#endif
#if MLIP_EXTEND_REF_LINE
  UInt ExternRefNum = pcCU->getExternRef(chType, uiZorderIdxInPart);
  if (ExternRefNum > 1 && compID != COMPONENT_Y)
      ExternRefNum = 1;
#endif
#if JVET_C0024_QTBT
  assert(uiZorderIdxInPart==0);
  const Int  iBaseUnitSize    = sps.getCTUSize() >> sps.getMaxTotalCUDepth();
#else
  const Int  iBaseUnitSize    = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth();
#endif
  const Int  iUnitWidth       = iBaseUnitSize  >> pcCU->getPic()->getPicYuvRec()->getComponentScaleX(compID);
  const Int  iUnitHeight      = iBaseUnitSize  >> pcCU->getPic()->getPicYuvRec()->getComponentScaleY(compID);
  const Int  iTUWidthInUnits  = uiTuWidth  / iUnitWidth;
  const Int  iTUHeightInUnits = uiTuHeight / iUnitHeight;
#if JVET_C0024_QTBT
  const Int  iAboveUnits      = iTUWidthInUnits  + iTUHeightInUnits;
  const Int  iLeftUnits       = iTUHeightInUnits + iTUWidthInUnits;
#else
  const Int  iAboveUnits      = iTUWidthInUnits  << 1;
  const Int  iLeftUnits       = iTUHeightInUnits << 1;
#endif
  const Int  bitDepthForChannel = sps.getBitDepth(chType);

  assert(iTUHeightInUnits > 0 && iTUWidthInUnits > 0);

  const Int  iPartIdxStride   = pcCU->getPic()->getNumPartInCtuWidth();
  const UInt uiPartIdxLT      = pcCU->getZorderIdxInCtu() + uiZorderIdxInPart;
  const UInt uiPartIdxRT      = g_auiRasterToZscan[ g_auiZscanToRaster[ uiPartIdxLT ] +   iTUWidthInUnits  - 1                   ];
  const UInt uiPartIdxLB      = g_auiRasterToZscan[ g_auiZscanToRaster[ uiPartIdxLT ] + ((iTUHeightInUnits - 1) * iPartIdxStride)];

  Int   iPicStride = pcCU->getPic()->getStride(compID);
  Bool  bNeighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  Int   iNumIntraNeighbor = 0;

  bNeighborFlags[iLeftUnits] = isAboveLeftAvailable( pcCU, uiPartIdxLT );
  iNumIntraNeighbor += bNeighborFlags[iLeftUnits] ? 1 : 0;
  iNumIntraNeighbor  += isAboveAvailable     ( pcCU, uiPartIdxLT, uiPartIdxRT, (bNeighborFlags + iLeftUnits + 1)                    );
#if JVET_C0024_QTBT
  iNumIntraNeighbor  += isAboveRightAvailable( pcCU, uiPartIdxRT, iTUHeightInUnits, (bNeighborFlags + iLeftUnits + 1 + iTUWidthInUnits ) );
#else
  iNumIntraNeighbor  += isAboveRightAvailable( pcCU, uiPartIdxLT, uiPartIdxRT, (bNeighborFlags + iLeftUnits + 1 + iTUWidthInUnits ) );
#endif
  iNumIntraNeighbor  += isLeftAvailable      ( pcCU, uiPartIdxLT, uiPartIdxLB, (bNeighborFlags + iLeftUnits - 1)                    );
#if JVET_C0024_QTBT
  iNumIntraNeighbor  += isBelowLeftAvailable ( pcCU, uiPartIdxLB, iTUWidthInUnits, (bNeighborFlags + iLeftUnits - 1 - iTUHeightInUnits) );
#else
  iNumIntraNeighbor  += isBelowLeftAvailable ( pcCU, uiPartIdxLT, uiPartIdxLB, (bNeighborFlags + iLeftUnits - 1 - iTUHeightInUnits) );
#endif

  const UInt         uiROIWidth  = uiTuWidth2+1;
  const UInt         uiROIHeight = uiTuHeight2+1;

  assert(uiROIWidth*uiROIHeight <= m_iYuvExtSize);

#if DEBUG_STRING
  std::stringstream ss(stringstream::out);
#endif

  {
    Pel *piIntraTemp   = m_piYuvExt[compID][PRED_BUF_UNFILTERED];
#if MLIP_EXTEND_REF_LINE
    Pel *piSecondIntraTemp = m_piSecondYuvExt[compID][PRED_BUF_UNFILTERED];
#endif
    Pel *piRoiOrigin = NULL;
#if BARC
    if (isCTUBoundary)
      piRoiOrigin = pcCU->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiZorderIdxInPart);
    else
#endif
      piRoiOrigin = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiZorderIdxInPart);
#if O0043_BEST_EFFORT_DECODING
    const Int  bitDepthForChannelInStream = sps.getStreamBitDepth(chType);
    fillReferenceSamples (bitDepthForChannelInStream, bitDepthForChannelInStream - bitDepthForChannel,
#else
    fillReferenceSamples (bitDepthForChannel,
#endif
                          piRoiOrigin, piIntraTemp, bNeighborFlags, iNumIntraNeighbor,  iUnitWidth, iUnitHeight, iAboveUnits, iLeftUnits,
                          uiROIWidth, uiROIHeight, iPicStride
#if MLIP_EXTEND_REF_LINE
        , ExternRefNum
        , piSecondIntraTemp
#endif
    );

#if DEBUG_STRING
    if (DebugOptionList::DebugString_Pred.getInt()&DebugStringGetPredModeMask(MODE_INTRA))
    {
      ss << "###: generating Ref Samples for channel " << compID << " and " << rTu.getRect(compID).width << " x " << rTu.getRect(compID).height << "\n";
      for (UInt y=0; y<uiROIHeight; y++)
      {
        ss << "###: - ";
        for (UInt x=0; x<uiROIWidth; x++)
        {
          if (x==0 || y==0)
          {
            ss << piIntraTemp[y*uiROIWidth + x] << ", ";
//          if (x%16==15) ss << "\nPart size: ~ ";
          }
        }
        ss << "\n";
      }
    }
#endif

#if RRC_INTRA_SHARPENING
    if (!bSharpRefSamples)
#endif
    if (bFilterRefSamples)
    {
      // generate filtered intra prediction samples

            Int          stride    = uiROIWidth;
      const Pel         *piSrcPtr  = piIntraTemp                           + (stride * uiTuHeight2); // bottom left
            Pel         *piDestPtr = m_piYuvExt[compID][PRED_BUF_FILTERED] + (stride * uiTuHeight2); // bottom left
#if MLIP_EXTEND_REF_LINE
            Int         Secondstride = uiROIWidth + ExternRefNum;
            const Pel   *piSecondSrcPtr = piSecondIntraTemp + (Secondstride * (uiTuHeight2 + ExternRefNum));
            Pel         *piSecondDestPtr = m_piSecondYuvExt[compID][PRED_BUF_FILTERED] + (Secondstride * (uiTuHeight2 + ExternRefNum));
#endif
      //------------------------------------------------

      Bool useStrongIntraSmoothing = isLuma(chType) && sps.getUseStrongIntraSmoothing();
      const Pel bottomLeft = piIntraTemp[stride * uiTuHeight2];
      const Pel topLeft    = piIntraTemp[0];
      const Pel topRight   = piIntraTemp[uiTuWidth2];
#if MLIP_EXTEND_REF_LINE
      const Pel SecondbottomLeft = piSecondIntraTemp[Secondstride * (uiTuHeight2 + ExternRefNum)];
      const Pel SecondtopLeft = piSecondIntraTemp[0];
      const Pel SecondtopRight = piSecondIntraTemp[uiTuWidth2 + ExternRefNum];
      Bool SeconduseStrongIntraSmoothing = useStrongIntraSmoothing;
#if !JVET_G0104_PLANAR_PDPC
      const Bool bIsSecondWeakSmoothing = bIsWeakSmoothing;
#endif
#endif
      if (useStrongIntraSmoothing)
      {
#if O0043_BEST_EFFORT_DECODING
        const Int  threshold     = 1 << (bitDepthForChannelInStream - 5);
#else
        const Int  threshold     = 1 << (bitDepthForChannel - 5);
#endif
        const Bool bilinearLeft  = abs((bottomLeft + topLeft ) - (2 * piIntraTemp[stride * uiTuHeight])) < threshold; //difference between the
        const Bool bilinearAbove = abs((topLeft    + topRight) - (2 * piIntraTemp[         uiTuWidth ])) < threshold; //ends and the middle
        if ((uiTuWidth < 32) || (!bilinearLeft) || (!bilinearAbove))
        {
          useStrongIntraSmoothing = false;
#if MLIP_EXTEND_REF_LINE
          const Bool SecondbilinearLeft = abs((SecondbottomLeft + SecondtopLeft) - (2 * piSecondIntraTemp[Secondstride * uiTuHeight])) < threshold; //difference between the
          const Bool SecondbilinearAbove = abs((SecondtopLeft + SecondtopRight) - (2 * piSecondIntraTemp[uiTuWidth])) < threshold; //ends and the middle
          if ((uiTuWidth < 32) || (!SecondbilinearLeft) || (!SecondbilinearAbove))
          {
              SeconduseStrongIntraSmoothing = false;
          }
#endif
        }
      }

      *piDestPtr = *piSrcPtr; // bottom left is not filtered
      piDestPtr -= stride;
      piSrcPtr  -= stride;
#if MLIP_EXTEND_REF_LINE
      *piSecondDestPtr = *piSecondSrcPtr; // bottom left is not filtered
      piSecondDestPtr -= Secondstride;
      piSecondSrcPtr -= Secondstride;
#endif
      //------------------------------------------------

      //left column (bottom to top)

      if (useStrongIntraSmoothing)
      {
#if JVET_C0024_QTBT
        const Int shift = g_aucConvertToBit[uiTuHeight] + MIN_CU_LOG2 + 1; //log2(uiTuHeight2)    //it is a bug for non-square PU for strong filter, JCA
#else
        const Int shift = g_aucConvertToBit[uiTuHeight] + 3; //log2(uiTuHeight2)    //it is a bug for non-square PU for strong filter, 
#endif

        for(UInt i=1; i<uiTuHeight2; i++, piDestPtr-=stride)
        {
          *piDestPtr = (((uiTuHeight2 - i) * bottomLeft) + (i * topLeft) + uiTuHeight) >> shift;
        }

        piSrcPtr -= stride * (uiTuHeight2 - 1);
      }
      else
      {
        for(UInt i=1; i<uiTuHeight2; i++, piDestPtr-=stride, piSrcPtr-=stride)
        {
          *piDestPtr = ( piSrcPtr[stride] + 2*piSrcPtr[0] + piSrcPtr[-stride] + 2 ) >> 2;
        }
      }
#if MLIP_EXTEND_REF_LINE
      if (SeconduseStrongIntraSmoothing)
      {
          const Int shift = g_aucConvertToBit[uiTuHeight] + 3; //log2(uiTuHeight2)

          for (UInt i = 1; i < uiTuHeight2 + ExternRefNum; i++, piSecondDestPtr -= Secondstride)
          {
              *piSecondDestPtr = (((uiTuHeight2 + ExternRefNum - i) * SecondbottomLeft) + (i * SecondtopLeft) + uiTuHeight) >> shift;
          }

          piSecondSrcPtr -= Secondstride * (uiTuHeight2 - 1 + ExternRefNum);
      }
      else
      {
          for (UInt i = 1; i < uiTuHeight2; i++, piDestPtr -= stride, piSrcPtr -= stride)
          {
              *piSecondDestPtr = (piSecondSrcPtr[Secondstride] + 2 * piSecondSrcPtr[0] + piSecondSrcPtr[-Secondstride] + 2) >> 2;
          }
      }
#endif
      //------------------------------------------------

      //top-left

      if (useStrongIntraSmoothing)
      {
        *piDestPtr = piSrcPtr[0];
      }
      else
      {
        *piDestPtr = ( piSrcPtr[stride] + 2*piSrcPtr[0] + piSrcPtr[1] + 2 ) >> 2;
      }
      piDestPtr += 1;
      piSrcPtr  += 1;
#if MLIP_EXTEND_REF_LINE
      if (SeconduseStrongIntraSmoothing)
      {
          *piSecondDestPtr = piSecondSrcPtr[0];
      }
      else
      {
          *piSecondDestPtr = (piSecondSrcPtr[Secondstride] + 2 * piSecondSrcPtr[0] + piSecondSrcPtr[1] + 2) >> 2;
      }
      piSecondDestPtr += 1;
      piSecondSrcPtr += 1;
#endif
      //------------------------------------------------

      //top row (left-to-right)

      if (useStrongIntraSmoothing)
      {
#if JVET_C0024_QTBT
        const Int shift = g_aucConvertToBit[uiTuWidth] + MIN_CU_LOG2 + 1; //log2(uiTuWidth2)
#else
        const Int shift = g_aucConvertToBit[uiTuWidth] + 3; //log2(uiTuWidth2)
#endif

        for(UInt i=1; i<uiTuWidth2; i++, piDestPtr++)
        {
          *piDestPtr = (((uiTuWidth2 - i) * topLeft) + (i * topRight) + uiTuWidth) >> shift;
        }

        piSrcPtr += uiTuWidth2 - 1;
      }
      else
      {
        for(UInt i=1; i<uiTuWidth2; i++, piDestPtr++, piSrcPtr++)
        {
          *piDestPtr = ( piSrcPtr[1] + 2*piSrcPtr[0] + piSrcPtr[-1] + 2 ) >> 2;
        }
      }
#if MLIP_EXTEND_REF_LINE
      if (SeconduseStrongIntraSmoothing)
      {
          const Int shift = g_aucConvertToBit[uiTuWidth] + 3; //log2(uiTuWidth2)

          for (UInt i = 1; i < uiTuWidth2 + ExternRefNum; i++, piSecondDestPtr++)
          {
              *piSecondDestPtr = (((uiTuWidth2 + ExternRefNum - i) * SecondtopLeft) + (i * SecondtopRight) + uiTuWidth) >> shift;
          }
          piSecondSrcPtr += uiTuWidth2 - 1 + ExternRefNum;
      }
      else
      {
          for (UInt i = 1; i < uiTuWidth2; i++, piDestPtr++, piSrcPtr++)
          {
              *piSecondDestPtr = (piSecondSrcPtr[1] + 2 * piSecondSrcPtr[0] + piSecondSrcPtr[-1] + 2) >> 2;
          }
      }
#endif

      //------------------------------------------------

      *piDestPtr=*piSrcPtr; // far right is not filtered
#if MLIP_EXTEND_REF_LINE
      *piSecondDestPtr = *piSecondSrcPtr;
#endif

#if DEBUG_STRING
    if (DebugOptionList::DebugString_Pred.getInt()&DebugStringGetPredModeMask(MODE_INTRA))
    {
      ss << "###: filtered result for channel " << compID <<"\n";
      for (UInt y=0; y<uiROIHeight; y++)
      {
        ss << "###: - ";
        for (UInt x=0; x<uiROIWidth; x++)
        {
          if (x==0 || y==0)
          {
            ss << m_piYuvExt[compID][PRED_BUF_FILTERED][y*uiROIWidth + x] << ", ";
//          if (x%16==15) ss << "\nPart size: ~ ";
          }
        }
        ss << "\n";
      }
    }
#endif
    }
#if RRC_INTRA_SHARPENING
    if (bSharpRefSamples)
    {
      const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();

      Int             stride = uiROIWidth;
      const Pel   *piSrcPtr = piIntraTemp + (stride * uiTuHeight2); // bottom left
      Pel         *piDestPtr = m_piYuvExt[compID][PRED_BUF_FILTERED] + (stride * uiTuHeight2); // bottom left

      const Pel bottomLeft = piIntraTemp[stride * uiTuHeight2];
      const Pel topLeft = piIntraTemp[0];
      const Pel topRight = piIntraTemp[uiTuWidth2];

      Bool ret = initSharpPattern(piIntraTemp, m_piYuvExt[compID][PRED_BUF_FILTERED], uiTuWidth2 + 1, uiTuWidth2 + 1, uiTuHeight2 + 1,
        rTu,
        pcCU, uiAbsPartIdx,
        m_piYuvSharpExt, m_iYuvExtSize,
        m_piYuvSharpExt2, m_iYuvExtSize,
        m_piYuvSharpExt3, m_iYuvExtSize,
        bitDepthForChannel
        );
    }
#endif
  }
  DEBUG_STRING_APPEND(sDebug, ss.str())
}

#if DIMD_INTRA_PRED
Void TComPrediction::initDIMDLumaFlexibleIntraPattern(TComDataCU* pcCU, UInt uiZorderIdxInPart, Int iRefX, Int iRefY, UInt uiRefWidth, UInt uiRefHeight, Bool bFilterRefSamples
#if BARC_DIMD
  , Bool isBARC, Bool isInCurCtu, UInt uiRefZorderIdxPart
#endif
  )
{
  const TComSPS &sps = *(pcCU->getSlice()->getSPS());
  const UInt uiTuWidth = uiRefWidth;
  const UInt uiTuHeight = uiRefHeight;

  const UInt uiTuWidth2 = uiTuWidth + uiTuHeight;
  const UInt uiTuHeight2 = uiTuHeight + uiTuWidth;

  const UInt uiROIWidth = uiTuWidth2 + 1;
  const UInt uiROIHeight = uiTuHeight2 + 1;

  const Int  iUnitSize = sps.getCTUSize() >> sps.getMaxTotalCUDepth();  // 4
                                      //assert(isLuma(pcCU->getTextType()));

  const Int  iTotalSamples = uiTuWidth2 + uiTuHeight2 + 1;
  const Int  bitDepthForChannel = sps.getBitDepth(CHANNEL_TYPE_LUMA);

  //assert(pcCU->getZorderIdxInCtu() == 0);
  Pel *piRoiOrigin = pcCU->getPic()->getPicYuvRec()->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiZorderIdxInPart);
#if BARC_DIMD
  if (isBARC && !isInCurCtu)
    piRoiOrigin = pcCU->getPic()->getPicYuvDScaleRef()->getAddrDscale(COMPONENT_Y, pcCU->getCtuRsAddr(), uiRefZorderIdxPart);
#endif
  Int  iPicStride = pcCU->getPic()->getStride(COMPONENT_Y);
  Int  iOffsetX = iRefX - (pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiZorderIdxInPart]]);
  Int  iOffsetY = iRefY - (pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiZorderIdxInPart]]);
  piRoiOrigin += iOffsetY * iPicStride + iOffsetX;

  Pel  piAdiLine[5 * MAX_CU_SIZE];
  Bool bNeighSampleAvail[5 * MAX_CU_SIZE];
  Int  iDCValue = 1 << (bitDepthForChannel - 1);

  ///< initialization
  for (Int i = 0; i < iTotalSamples; i++)
  {
    piAdiLine[i] = iDCValue;
    bNeighSampleAvail[i] = false;
  }

  TComDataCU* pcRefCU = NULL;
  UInt        uiRefPartIdx = 0;
  Pel*        piRoiTemp = piRoiOrigin - iPicStride - 1;

  Int iOrigX = iRefX - 1;
  Int iOrigY = iRefY - 1;
  Int iEndX = iOrigX + uiTuWidth2;
  Int iEndY = iOrigY + uiTuHeight2;
  Int iNumAvailNeighborSamples = 0;

  Pel*  piAdiLineTemp = piAdiLine + uiTuHeight2;
  Bool* pbNeighSampleAvail = bNeighSampleAvail + uiTuHeight2;

  // above-left reference pixel
  pcRefCU = pcCU->getOrigRefReg(uiRefPartIdx, iOrigX, iOrigY);
  if (pcRefCU != NULL && pcCU->isReferenceAvailable(pcRefCU, uiRefPartIdx))
  {
    piAdiLineTemp[0] = piRoiTemp[0];
    pbNeighSampleAvail[0] = true;
    iNumAvailNeighborSamples++;
  }

  ///< above reference pixels, width + height
  Int iScanX = iOrigX + 1;
  while (iScanX <= iEndX)
  {
    Int iSegLength = iUnitSize - (iScanX%iUnitSize);
    if (iScanX + iSegLength > iEndX)
    {
      iSegLength = iEndX - iScanX + 1;
    }
    pcRefCU = pcCU->getOrigRefReg(uiRefPartIdx, iScanX, iOrigY
#if BARC_DIMD
      ,isBARC
#endif
    );
    if (pcRefCU && pcCU->isReferenceAvailable(pcRefCU, uiRefPartIdx))
    {
      for (Int iIdx = 0; iIdx < iSegLength; iIdx++)
      {
        Int iTempIdx = iScanX + iIdx - iOrigX;
        piAdiLineTemp[iTempIdx] = piRoiTemp[iTempIdx];
        pbNeighSampleAvail[iTempIdx] = true;
      }
      iNumAvailNeighborSamples += iSegLength;
    }
    iScanX += iSegLength;
  }
  //assert(iScanX == iEndX+1);

  ///< left reference pixels, width + height
  Int iScanY = iOrigY + 1;
  while (iScanY <= iEndY)
  {
    Int iSegLength = iUnitSize - (iScanY%iUnitSize);
    if (iScanY + iSegLength > iEndY)
    {
      iSegLength = iEndY - iScanY + 1;
    }
    pcRefCU = pcCU->getOrigRefReg(uiRefPartIdx, iOrigX, iScanY
#if BARC_DIMD
      ,isBARC
#endif
    );
    if (pcRefCU && pcCU->isReferenceAvailable(pcRefCU, uiRefPartIdx))
    {
      for (Int iIdx = 0; iIdx < iSegLength; iIdx++)
      {
        Int iTempIdx = iScanY + iIdx - iOrigY;
        piAdiLineTemp[-iTempIdx] = piRoiTemp[iTempIdx*iPicStride];
        pbNeighSampleAvail[-iTempIdx] = true;
      }
      iNumAvailNeighborSamples += iSegLength;
    }
    iScanY += iSegLength;
  }
  //assert(iScanY = iEndY+1);

  //assert(iNumAvailNeighborSamples <= uiTuWidth2 + uiTuHeight2 + 1);
  if (iNumAvailNeighborSamples != 0 && iNumAvailNeighborSamples < uiTuWidth2 + uiTuHeight2 + 1)
  {
    Int iCurr = 0, iNext = 1;
    Pel piRef = 0;
    while (iCurr < iTotalSamples)
    {
      if (!bNeighSampleAvail[iCurr])
      {
        if (iCurr == 0)
        {
          while (iNext < iTotalSamples && !bNeighSampleAvail[iNext])
          {
            iNext++;
          }
          //assert(iNext < iTotalSamples);
          piRef = piAdiLine[iNext];

          while (iCurr < iNext)
          {
            piAdiLine[iCurr] = piRef;
            iCurr++;
          }
        }
        else
        {
          piAdiLine[iCurr] = piAdiLine[iCurr - 1];
          iCurr++;
        }
      }
      else
      {
        iCurr++;
      }
    }
  }

  //assert(uiROIWidth*uiROIHeight <= m_iYuvExtSize);

  Pel *piIntraTemp = m_piYuvExt[COMPONENT_Y][PRED_BUF_UNFILTERED];
  piAdiLineTemp = piAdiLine + uiTuHeight2;
  for (Int i = 0; i < uiROIWidth; i++)
  {
    piIntraTemp[i] = piAdiLineTemp[i];
  }
  for (Int i = 1; i < uiROIHeight; i++)
  {
    piIntraTemp[i*uiROIWidth] = piAdiLineTemp[-i];
  }

  if (bFilterRefSamples)
  {
    // generate filtered intra prediction samples

    Int          stride = uiROIWidth;
    const Pel    *piSrcPtr = piIntraTemp + (stride * uiTuHeight2); // bottom left
    Pel          *piDestPtr = m_piYuvExt[COMPONENT_Y][PRED_BUF_FILTERED] + (stride * uiTuHeight2); // bottom left

                                                     //------------------------------------------------

    Bool useStrongIntraSmoothing = isLuma(COMPONENT_Y) && sps.getUseStrongIntraSmoothing();
    const Pel bottomLeft = piIntraTemp[stride * uiTuHeight2];
    const Pel topLeft = piIntraTemp[0];
    const Pel topRight = piIntraTemp[uiTuWidth2];

    if (useStrongIntraSmoothing)
    {
      const Int  threshold = 1 << (bitDepthForChannel - 5);
      const Bool bilinearLeft = abs((bottomLeft + topLeft) - (2 * piIntraTemp[stride * uiTuHeight])) < threshold; //difference between the
      const Bool bilinearAbove = abs((topLeft + topRight) - (2 * piIntraTemp[uiTuWidth])) < threshold; //ends and the middle
      if ((uiTuWidth < 32) || (!bilinearLeft) || (!bilinearAbove))
      {
        useStrongIntraSmoothing = false;
      }
    }

    *piDestPtr = *piSrcPtr; // bottom left is not filtered
    piDestPtr -= stride;
    piSrcPtr -= stride;

    //------------------------------------------------
    //left column (bottom to top)

    if (useStrongIntraSmoothing)
    {
      const Int shift = g_aucConvertToBit[uiTuHeight] + MIN_CU_LOG2 + 1; //log2(uiTuHeight2)    //it is a bug for non-square PU for strong filter, JCA

      for (UInt i = 1; i < uiTuHeight2; i++, piDestPtr -= stride)
      {
        *piDestPtr = (((uiTuHeight2 - i) * bottomLeft) + (i * topLeft) + uiTuHeight) >> shift;
      }

      piSrcPtr -= stride * (uiTuHeight2 - 1);
    }
    else
    {
      for (UInt i = 1; i < uiTuHeight2; i++, piDestPtr -= stride, piSrcPtr -= stride)
      {
        *piDestPtr = (piSrcPtr[stride] + 2 * piSrcPtr[0] + piSrcPtr[-stride] + 2) >> 2;
      }
    }

    //------------------------------------------------

    //top-left

    if (useStrongIntraSmoothing)
    {
      *piDestPtr = piSrcPtr[0];
    }
    else
    {
      *piDestPtr = (piSrcPtr[stride] + 2 * piSrcPtr[0] + piSrcPtr[1] + 2) >> 2;
    }
    piDestPtr += 1;
    piSrcPtr += 1;

    //------------------------------------------------

    //top row (left-to-right)

    if (useStrongIntraSmoothing)
    {
      const Int shift = g_aucConvertToBit[uiTuWidth] + MIN_CU_LOG2 + 1; //log2(uiTuWidth2)

      for (UInt i = 1; i < uiTuWidth2; i++, piDestPtr++)
      {
        *piDestPtr = (((uiTuWidth2 - i) * topLeft) + (i * topRight) + uiTuWidth) >> shift;
      }

      piSrcPtr += uiTuWidth2 - 1;
    }
    else
    {
      for (UInt i = 1; i < uiTuWidth2; i++, piDestPtr++, piSrcPtr++)
      {
        *piDestPtr = (piSrcPtr[1] + 2 * piSrcPtr[0] + piSrcPtr[-1] + 2) >> 2;
      }
    }

    //------------------------------------------------

    *piDestPtr = *piSrcPtr; // far right is not filtered
  }
}

Bool TComPrediction::filteringDIMDIntraReferenceSamples(const ComponentID compID, UInt uiDirMode, UInt uiTuChWidth, UInt uiTuChHeight, const ChromaFormat chFmt, const Bool intraReferenceSmoothingDisabled )
{
  Bool bFilter;

  if (!filterIntraReferenceSamples(toChannelType(compID), chFmt, intraReferenceSmoothingDisabled))
  {
    bFilter = false;
  }
  else
  {
    //assert(uiTuChWidth >= MIN_PU_SIZE && uiTuChHeight >= MIN_PU_SIZE && uiTuChWidth <= MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE && uiTuChHeight <= MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE);

    if (uiDirMode == DC_IDX)
    {
      bFilter = false;
    }
    else
    {
#if DIMD_NUM_INTRA_DIR_INC
      assert(uiDirMode >= PLANAR_IDX && uiDirMode <= EXT_VDIA_IDX);
      Int diff = min<Int>(abs((Int)uiDirMode - EXT_HOR_IDX), abs((Int)uiDirMode - EXT_VER_IDX));
#else
      Int diff = min<Int>(abs((Int)uiDirMode - HOR_IDX), abs((Int)uiDirMode - VER_IDX));
#endif

      UInt uiLog2BlkWidthSize = 0, uiLog2BlkHeightSize = 0;
      while ((uiTuChWidth >> uiLog2BlkWidthSize) > 1) uiLog2BlkWidthSize++;
      while ((uiTuChHeight >> uiLog2BlkHeightSize) > 1) uiLog2BlkHeightSize++;

      UInt log2Size = (uiLog2BlkWidthSize + uiLog2BlkHeightSize) >> 1;
      UInt sizeIndex = log2Size - 1;
      //assert(sizeIndex < MAX_INTRA_FILTER_DEPTHS);

#if DIMD_NUM_INTRA_DIR_INC
      bFilter = diff > m_aucExtIntraFilter[toChannelType(compID)][sizeIndex];
#else
      bFilter = diff > m_aucIntraFilter[toChannelType(compID)][sizeIndex];
#endif
    }
  }

  return bFilter;
}

Void TComPrediction::predDIMDIntraLumaAng(TComDataCU* pcCU, UInt uiDirMode, Pel* piPred, UInt uiStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  const TComSPS &sps = *(pcCU->getSlice()->getSPS());
  const ChromaFormat   chFmt = pcCU->getPic()->getChromaFormat();
  //assert(isLuma(pcCU->getTextType()));

  Bool bUseFilteredPredictions = TComPrediction::filteringDIMDIntraReferenceSamples(COMPONENT_Y, uiDirMode, iWidth, iHeight, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag());

  const Pel *ptrSrc = getPredictorPtr(COMPONENT_Y, bUseFilteredPredictions);
  const Int sw = iWidth + iHeight + 1;
  Pel *pDst = piPred;

  if (uiDirMode == PLANAR_IDX)
  {
    xPredDIMDIntraPlanar(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight);
  }
  else
  {
    const ChannelType    channelType = toChannelType(COMPONENT_Y);
    const Bool           enableEdgeFilters = true;
    const Int   channelsBitDepthForPrediction = pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);
#if VCEG_AZ07_INTRA_4TAP_FILTER
    const Bool  enable4TapFilter = pcCU->getSlice()->getSPS()->getUseIntra4TapFilter();
#endif
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
    const Bool  enableBoundaryFilter = pcCU->getSlice()->getSPS()->getUseIntraBoundaryFilter();
#endif

    xPredDIMDIntraAng(channelsBitDepthForPrediction, ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, channelType, uiDirMode, enableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
      , enable4TapFilter
#endif
      , eTempType, iTemplateWidth, iTemplateHeight);

    if (uiDirMode == DC_IDX)
    {
      xDCPredFiltering(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, channelType);
    }
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
    else if (enableBoundaryFilter)
    {
#if DIMD_NUM_INTRA_DIR_INC
      if (uiDirMode == EXT_VDIA_IDX)
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
      if (uiDirMode == VDIA_IDX)
#else
      if (uiDirMode == 34)
#endif
#endif
      {
        xDIMDIntraPredFilteringMode34(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight);
      }
      else  if (uiDirMode == 2)
      {
        xDIMDIntraPredFilteringMode02(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight);
      }
#if DIMD_NUM_INTRA_DIR_INC
      else if ((uiDirMode <= 18 && uiDirMode > 2) || (uiDirMode >= (EXT_VDIA_IDX - 16) && uiDirMode < EXT_VDIA_IDX))
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
      else if ((uiDirMode <= 10 && uiDirMode > 2) || (uiDirMode >= (VDIA_IDX - 8) && uiDirMode < VDIA_IDX))
#else
      else if ((uiDirMode <= 6 && uiDirMode>2) || (uiDirMode >= 30 && uiDirMode < 34))
#endif
#endif
      {
        xDIMDIntraPredFilteringModeDGL(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, uiDirMode, eTempType, iTemplateWidth, iTemplateHeight);
      }
    }
#endif
  }
}

Void TComPrediction::xPredDIMDIntraPlanar(const Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  //assert(width <= MAX_CU_SIZE+MAX_TEMP_SIZE && height <= MAX_CU_SIZE+MAX_TEMP_SIZE);
  Int leftColumn[MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE + 1], topRow[MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE + 1], bottomRow[MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE], rightColumn[MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE];
  UInt uiDiv1Dhor = width;
  UInt uiDiv1Dver = height;
  UInt delt = width * height;

  /*
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight > 0);
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight == 0);
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth == 0 && iTemplateHeight > 0);
  }
  */

  if (eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    for (Int k = 0; k < width; k++)
    {
      topRow[k] = pSrc[k - srcStride];
    }

    for (Int k = 0; k < height; k++)
    {
      leftColumn[k] = pSrc[k*srcStride - 1];
    }

    Int bottomLeft = pSrc[Int(height)*srcStride - 1];
    Int topRight = pSrc[Int(width) - srcStride];

    for (Int k = 0; k < width; k++)
    {
      bottomRow[k] = bottomLeft - topRow[k];
      topRow[k] *= uiDiv1Dver;
    }

    for (Int k = 0; k < height; k++)
    {
      rightColumn[k] = topRight - leftColumn[k];
      leftColumn[k] *= uiDiv1Dhor;
    }

    const UInt topRowShift = 0;

    for (Int y = 0; y < height; y++)
    {
      Int horPred = leftColumn[y];

      if (y < iTemplateHeight)
      {
        for (Int x = 0; x < width; x++)
        {
          horPred += rightColumn[y];
          topRow[x] += bottomRow[x];

          Int vertPred = ((topRow[x] + topRowShift) >> topRowShift);
          if (x >= iTemplateWidth)
          {
            rpDst[y*dstStride + x] = ((horPred*uiDiv1Dver) + (vertPred*uiDiv1Dhor) + delt) / (uiDiv1Dhor*uiDiv1Dver * 2);
          }
        }
      }
      else
      {
        for (Int x = 0; x < iTemplateWidth; x++)
        {
          horPred += rightColumn[y];
          topRow[x] += bottomRow[x];

          Int vertPred = ((topRow[x] + topRowShift) >> topRowShift);
          rpDst[y*dstStride + x] = ((horPred*uiDiv1Dver) + (vertPred*uiDiv1Dhor) + delt) / (uiDiv1Dhor*uiDiv1Dver * 2);
        }
      }
    }
  }
  else if (eTempType == LEFT_NEIGHBOR)
  {
    for (Int k = 0; k < height; k++)
    {
      leftColumn[k] = pSrc[k*srcStride - 1];
    }
    for (Int k = 0; k < iTemplateWidth; k++)
    {
      topRow[k] = pSrc[k - srcStride];
    }

    Int bottomLeft = pSrc[Int(height)*srcStride - 1];
    Int topRight = pSrc[Int(width) - srcStride];

    for (Int k = 0; k < height; k++)
    {
      rightColumn[k] = topRight - leftColumn[k];
      leftColumn[k] *= uiDiv1Dhor;
    }

    for (Int k = 0; k < iTemplateWidth; k++)
    {
      bottomRow[k] = bottomLeft - topRow[k];
      topRow[k] *= uiDiv1Dver;
    }

    const UInt topRowShift = 0;

    for (Int y = 0; y < height; y++)
    {
      Int horPred = leftColumn[y];
      for (Int x = 0; x < iTemplateWidth; x++)
      {
        horPred += rightColumn[y];
        topRow[x] += bottomRow[x];

        Int vertPred = ((topRow[x] + topRowShift) >> topRowShift);
        rpDst[y*dstStride + x] = ((horPred*uiDiv1Dver) + (vertPred*uiDiv1Dhor) + delt) / (uiDiv1Dhor*uiDiv1Dver * 2);
      }
    }
  }
  else if (eTempType == ABOVE_NEIGHBOR)
  {
    for (Int k = 0; k < width; k++)
    {
      topRow[k] = pSrc[k - srcStride];
    }

    for (Int k = 0; k < iTemplateHeight; k++)
    {
      leftColumn[k] = pSrc[k*srcStride - 1];
    }

    Int bottomLeft = pSrc[Int(height)*srcStride - 1];
    Int topRight = pSrc[Int(width) - srcStride];

    for (Int k = 0; k < width; k++)
    {
      bottomRow[k] = bottomLeft - topRow[k];
      topRow[k] *= uiDiv1Dver;
    }

    for (Int k = 0; k < iTemplateHeight; k++)
    {
      rightColumn[k] = topRight - leftColumn[k];
      leftColumn[k] *= uiDiv1Dhor;
    }

    const UInt topRowShift = 0;

    for (Int y = 0; y < iTemplateHeight; y++)
    {
      Int horPred = leftColumn[y];
      for (Int x = 0; x < width; x++)
      {
        horPred += rightColumn[y];
        topRow[x] += bottomRow[x];

        Int vertPred = ((topRow[x] + topRowShift) >> topRowShift);
        rpDst[y*dstStride + x] = ((horPred*uiDiv1Dver) + (vertPred*uiDiv1Dhor) + delt) / (uiDiv1Dhor*uiDiv1Dver * 2);
      }
    }
  }
  else
  {
    assert(0);
  }
}

Void TComPrediction::xPredDIMDIntraAng(Int bitDepth, const Pel* pSrc, Int srcStride, Pel* pTrueDst, Int dstStrideTrue, UInt uiWidth, UInt uiHeight, ChannelType channelType, UInt dirMode, const Bool bEnableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
  , Bool enable4TapFilter
#endif 
  , TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  Int width = Int(uiWidth);
  Int height = Int(uiHeight);

  // Map the mode index to main prediction direction and angle
  //assert( dirMode != PLANAR_IDX ); //no planar
  const Bool modeDC = dirMode == DC_IDX;

  /*
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight > 0);
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight == 0);
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth == 0 && iTemplateHeight > 0);
  }
  */

  // Do the DC prediction
  if (modeDC)
  {
    const Pel dcval = predIntraGetPredValDC(pSrc, srcStride, width, height);

    if (eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      for (Int y = 0; y < height; y++, pTrueDst += dstStrideTrue)
      {
        if (y < iTemplateHeight)
        {
          for (Int x = iTemplateWidth; x < width;)
          {
            pTrueDst[x++] = dcval;
          }
        }
        else
        {
          for (Int x = 0; x < iTemplateWidth;)
          {
            pTrueDst[x++] = dcval;
          }
        }
      }
    }
    else if (eTempType == LEFT_NEIGHBOR)
    {
      for (Int y = 0; y < height; y++, pTrueDst += dstStrideTrue)
      {
        for (Int x = 0; x < iTemplateWidth;)
        {
          pTrueDst[x++] = dcval;
        }
      }
    }
    else if (eTempType == ABOVE_NEIGHBOR)
    {
      for (Int y = 0; y < iTemplateHeight; y++, pTrueDst += dstStrideTrue)
      {
        for (Int x = 0; x < width;)
        {
          pTrueDst[x++] = dcval;
        }
      }
    }
    else
    {
      assert(0);
    }
  }
  else // Do angular predictions
  {
#if DIMD_NUM_INTRA_DIR_INC
    const Bool       bIsModeVer = (dirMode >= EXT_DIA_IDX);
    const Int        intraPredAngleMode = (bIsModeVer) ? (Int)dirMode - EXT_VER_IDX : -((Int)dirMode - EXT_HOR_IDX);
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
    const Bool       bIsModeVer = (dirMode >= DIA_IDX);
#else
    const Bool       bIsModeVer = (dirMode >= 18);
#endif
    const Int        intraPredAngleMode = (bIsModeVer) ? (Int)dirMode - VER_IDX : -((Int)dirMode - HOR_IDX);
#endif
    const Int        absAngMode = abs(intraPredAngleMode);
    const Int        signAng = intraPredAngleMode < 0 ? -1 : 1;
    const Bool       edgeFilter = bEnableEdgeFilters && isLuma(channelType) && (width <= MAXIMUM_INTRA_FILTERED_WIDTH) && (height <= MAXIMUM_INTRA_FILTERED_HEIGHT);

    // Set bitshifts and scale the angle parameter to block size
#if DIMD_NUM_INTRA_DIR_INC
    static const Int extAngTable[33] = { 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 49, 52, 55, 58, 61, 64 };   ///< module 64
    static const Int extInvAngTable[33] = { 0, 16384, 8192, 5461, 4096, 3277, 2731, 2048, 1638, 1365, 1170, 1024, 910, 819, 745, 683, 630, 585, 546, 512, 482, 455, 431, 410, 390, 372, 356, 334, 315, 298, 282, 269, 256 };   ///< (256 * 64) / Angle
    Int invAngle = extInvAngTable[absAngMode];
    Int absAng = extAngTable[absAngMode];
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
    static const Int angTable[17] = { 0,    1,    2,    3,    5,    7,    9,   11,   13,   15,   17,   19,   21,   23,   26,   29,   32 };
    static const Int invAngTable[17] = { 0, 8192, 4096, 2731, 1638, 1170,  910,  745,  630,  546,  482,  431,  390,  356,  315,  282,  256 }; // (256 * 32) / Angle
#else
    static const Int angTable[9] = { 0,    2,    5,   9,  13,  17,  21,  26,  32 };
    static const Int invAngTable[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / Angle
#endif
    Int invAngle = invAngTable[absAngMode];
    Int absAng = angTable[absAngMode];
#endif
    Int intraPredAngle = signAng * absAng;

    Pel* refMain;
    Pel* refSide;

    //assert(width+height <= 2*MAX_CU_SIZE);
    Pel  refAbove[2 * MAX_CU_SIZE + 1];  ///< buffer may not be enough
    Pel  refLeft[2 * MAX_CU_SIZE + 1];   ///< buffer may not be enough

                       // Initialize the Main and Left reference array.
    if (intraPredAngle < 0)
    {
      const Int refMainOffsetPreScale = (bIsModeVer ? height : width) - 1;
      for (Int x = 0; x < width + 1; x++)
      {
        refAbove[x + height - 1] = pSrc[x - srcStride - 1];
      }
      for (Int y = 0; y<height + 1; y++)
      {
        refLeft[y + width - 1] = pSrc[(y - 1)*srcStride - 1];
      }
      refMain = (bIsModeVer ? refAbove + height : refLeft + width) - 1;
      refSide = (bIsModeVer ? refLeft + width : refAbove + height) - 1;

      // Extend the Main reference to the left.
      Int invAngleSum = 128;       // rounding for (shift by 8)
#if DIMD_NUM_INTRA_DIR_INC
      for (Int k = -1; k>(refMainOffsetPreScale + 1)*intraPredAngle >> 6; k--)
#else
      for (Int k = -1; k > (refMainOffsetPreScale + 1)*intraPredAngle >> 5; k--)
#endif
      {
        invAngleSum += invAngle;
        refMain[k] = refSide[invAngleSum >> 8];
      }
    }
    else
    {
      for (Int x = 0; x < width + height + 1; x++)
      {
        refAbove[x] = pSrc[x - srcStride - 1];
        refLeft[x] = pSrc[(x - 1)*srcStride - 1];
      }
      refMain = bIsModeVer ? refAbove : refLeft;
      refSide = bIsModeVer ? refLeft : refAbove;
    }

    // swap width/height if we are doing a horizontal mode:
    //assert(width <= MAX_CU_SIZE+MAX_TEMP_SIZE && height <= MAX_CU_SIZE+MAX_TEMP_SIZE);
    Pel tempArray[(MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE)*(MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE)];  ///< buffer size may not be big enough
    const Int dstStride = bIsModeVer ? dstStrideTrue : (MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE);
    Pel *pDst = bIsModeVer ? pTrueDst : tempArray;
    if (!bIsModeVer)
    {
      std::swap(width, height);
      std::swap(iTemplateWidth, iTemplateHeight);
    }

    if (intraPredAngle == 0)  // pure vertical or pure horizontal
    {
      if (eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        for (Int y = 0; y < height; y++)
        {
          if (y < iTemplateHeight)
          {
            for (Int x = iTemplateWidth; x < width; x++)
            {
              pDst[y*dstStride + x] = refMain[x + 1];
            }
          }
          else
          {
            for (Int x = 0; x < iTemplateWidth; x++)
            {
              pDst[y*dstStride + x] = refMain[x + 1];
            }
          }
        }

        if (edgeFilter)
        {
          for (Int y = iTemplateHeight; y < height; y++)
          {
            pDst[y*dstStride] = Clip3(0, ((1 << bitDepth) - 1), pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 1));
          }
        }
      }
      else if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
      {
        if ((eTempType == LEFT_NEIGHBOR && bIsModeVer) || (eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
        {
          for (Int y = 0; y < height; y++)
          {
            for (Int x = 0; x < iTemplateWidth; x++)
            {
              pDst[y*dstStride + x] = refMain[x + 1];
            }
          }

          if (edgeFilter)
          {
            for (Int y = 0; y < height; y++)
            {
              pDst[y*dstStride] = Clip3(0, ((1 << bitDepth) - 1), pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 1));
            }
          }
        }
        else
        {
          for (Int y = 0; y < iTemplateHeight; y++)
          {
            for (Int x = 0; x < width; x++)
            {
              pDst[y*dstStride + x] = refMain[x + 1];
            }
          }

          if (edgeFilter)
          {
            for (Int y = 0; y < iTemplateHeight; y++)
            {
              pDst[y*dstStride] = Clip3(0, ((1 << bitDepth) - 1), pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 1));
            }
          }
        }
      }
      else
      {
        assert(0);
      }
    }
    else
    {
      Pel *pDsty = pDst;

      if (eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        for (Int y = 0, deltaPos = intraPredAngle; y < height; y++, deltaPos += intraPredAngle, pDsty += dstStride)
        {
#if DIMD_NUM_INTRA_DIR_INC
          const Int deltaInt = deltaPos >> 6;
          const Int deltaFract = deltaPos & (64 - 1);
#else
          const Int deltaInt = deltaPos >> 5;
          const Int deltaFract = deltaPos & (32 - 1);
#endif

          Int iStartIdx, iEndIdx;
          if (y < iTemplateHeight)
          {
            iStartIdx = iTemplateWidth;
            iEndIdx = width - 1;
          }
          else
          {
            iStartIdx = 0;
            iEndIdx = iTemplateWidth - 1;
          }

          if (deltaFract)
          {
            // Do linear filtering
#if VCEG_AZ07_INTRA_4TAP_FILTER
            if (enable4TapFilter)
            {
              Int p[4], x, refMainIndex;
              const Pel nMin = 0, nMax = (1 << bitDepth) - 1;
#if DIMD_NUM_INTRA_DIR_INC
              Int *f = (width < 16) ? g_aiExtIntraCubicFilter[deltaFract] : g_aiExtIntraGaussFilter[deltaFract];
#else
              Int *f = (width < 16) ? g_aiIntraCubicFilter[deltaFract] : g_aiIntraGaussFilter[deltaFract];
#endif

              for (x = iStartIdx; x <= iEndIdx; x++)
              {
                refMainIndex = x + deltaInt + 1;

                p[1] = refMain[refMainIndex];
                p[2] = refMain[refMainIndex + 1];

                p[0] = x == 0 ? p[1] : refMain[refMainIndex - 1];
                p[3] = x == (width - 1) ? p[2] : refMain[refMainIndex + 2];

                pDst[y*dstStride + x] = (Pel)((f[0] * p[0] + f[1] * p[1] + f[2] * p[2] + f[3] * p[3] + 128) >> 8);

                if (width < 16) // for blocks larger than 8x8, Gaussian interpolation filter with positive coefficients is used, no Clipping is necessary
                {
                  pDst[y*dstStride + x] = Clip3(nMin, nMax, pDst[y*dstStride + x]);
                }
              }
            }
            else
            {
#endif
              const Pel *pRM = refMain + deltaInt + 1 + iStartIdx;
              Int lastRefMainPel = *pRM++;
              for (Int x = iStartIdx; x <= iEndIdx; pRM++, x++)
              {
                Int thisRefMainPel = *pRM;
#if DIMD_NUM_INTRA_DIR_INC
                pDsty[x + 0] = (Pel)(((64 - deltaFract)*lastRefMainPel + deltaFract*thisRefMainPel + 32) >> 6);
#else
                pDsty[x + 0] = (Pel)(((32 - deltaFract)*lastRefMainPel + deltaFract*thisRefMainPel + 16) >> 5);
#endif
                lastRefMainPel = thisRefMainPel;
              }
#if VCEG_AZ07_INTRA_4TAP_FILTER
            }
#endif
          }
          else
          {
            // Just copy the integer samples
            for (Int x = iStartIdx; x <= iEndIdx; x++)
            {
              pDsty[x] = refMain[x + deltaInt + 1];
            }
          }
        }
      }
      else if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
      {
        Int iRegionWidth, iRegionHeight;
        if ((eTempType == LEFT_NEIGHBOR && bIsModeVer) || (eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
        {
          iRegionWidth = iTemplateWidth;
          iRegionHeight = height;
        }
        else
        {
          iRegionWidth = width;
          iRegionHeight = iTemplateHeight;
        }

        for (Int y = 0, deltaPos = intraPredAngle; y < iRegionHeight; y++, deltaPos += intraPredAngle, pDsty += dstStride)
        {
#if DIMD_NUM_INTRA_DIR_INC
          const Int deltaInt = deltaPos >> 6;
          const Int deltaFract = deltaPos & (64 - 1);
#else
          const Int deltaInt = deltaPos >> 5;
          const Int deltaFract = deltaPos & (32 - 1);
#endif

          if (deltaFract)
          {
            // Do linear filtering
#if VCEG_AZ07_INTRA_4TAP_FILTER
            if (enable4TapFilter)
            {
              Int p[4], x, refMainIndex;
              const Pel nMin = 0, nMax = (1 << bitDepth) - 1;
#if DIMD_NUM_INTRA_DIR_INC
              Int *f = (width < 16) ? g_aiExtIntraCubicFilter[deltaFract] : g_aiExtIntraGaussFilter[deltaFract];
#else
              Int *f = (width < 16) ? g_aiIntraCubicFilter[deltaFract] : g_aiIntraGaussFilter[deltaFract];
#endif

              for (x = 0; x < iRegionWidth; x++)
              {
                refMainIndex = x + deltaInt + 1;

                p[1] = refMain[refMainIndex];
                p[2] = refMain[refMainIndex + 1];

                p[0] = x == 0 ? p[1] : refMain[refMainIndex - 1];
                p[3] = x == (width - 1) ? p[2] : refMain[refMainIndex + 2];

                pDst[y*dstStride + x] = (Pel)((f[0] * p[0] + f[1] * p[1] + f[2] * p[2] + f[3] * p[3] + 128) >> 8);

                if (width < 16) // for blocks larger than 8x8, Gaussian interpolation filter with positive coefficients is used, no Clipping is necessary
                {
                  pDst[y*dstStride + x] = Clip3(nMin, nMax, pDst[y*dstStride + x]);
                }
              }
            }
            else
            {
#endif
              const Pel *pRM = refMain + deltaInt + 1;
              Int lastRefMainPel = *pRM++;
              for (Int x = 0; x < iRegionWidth; pRM++, x++)
              {
                Int thisRefMainPel = *pRM;
#if DIMD_NUM_INTRA_DIR_INC
                pDsty[x + 0] = (Pel)(((64 - deltaFract)*lastRefMainPel + deltaFract*thisRefMainPel + 32) >> 6);
#else
                pDsty[x + 0] = (Pel)(((32 - deltaFract)*lastRefMainPel + deltaFract*thisRefMainPel + 16) >> 5);
#endif
                lastRefMainPel = thisRefMainPel;
              }
#if VCEG_AZ07_INTRA_4TAP_FILTER
            }
#endif
          }
          else
          {
            // Just copy the integer samples
            for (Int x = 0; x < iRegionWidth; x++)
            {
              pDsty[x] = refMain[x + deltaInt + 1];
            }
          }
        }
      }
      else
      {
        assert(0);
      }

#if VCEG_AZ07_INTRA_65ANG_MODES
#if DIMD_NUM_INTRA_DIR_INC
      if (edgeFilter && absAng <= 2)
#else
      if (edgeFilter && absAng <= 1)
#endif
      {
        if (eTempType == LEFT_ABOVE_NEIGHBOR)
        {
          for (Int y = iTemplateHeight; y < height; y++)
          {
            pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
          }
        }
        else if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
        {
          if ((eTempType == LEFT_NEIGHBOR && bIsModeVer) || (eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
          {
            for (Int y = 0; y < height; y++)
            {
              pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
            }
          }
          else
          {
            for (Int y = 0; y < iTemplateHeight; y++)
            {
              pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
            }
          }
        }
        else
        {
          assert(0);
        }
      }
#elif DIMD_NUM_INTRA_DIR_INC
      if (edgeFilter && absAng <= 2)
      {
        if (eTempType == LEFT_ABOVE_NEIGHBOR)
        {
          for (Int y = iTemplateHeight; y < height; y++)
          {
            pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
          }
        }
        else if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
        {
          if ((eTempType == LEFT_NEIGHBOR && bIsModeVer) || (eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
          {
            for (Int y = 0; y < height; y++)
            {
              pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
            }
          }
          else
          {
            for (Int y = 0; y < iTemplateHeight; y++)
            {
              pDst[y*dstStride] = Clip3(0, (1 << bitDepth) - 1, pDst[y*dstStride] + ((refSide[y + 1] - refSide[0]) >> 2));
            }
          }
        }
        else
        {
          assert(0);
        }
      }
#endif
    }

    // Flip the block if this is the horizontal mode
    if (!bIsModeVer)
    {
      if (eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        for (Int y = 0; y < height; y++)
        {
          Int iStartIdx, iEndIdx;
          if (y < iTemplateHeight)
          {
            iStartIdx = iTemplateWidth;
            iEndIdx = width - 1;
          }
          else
          {
            iStartIdx = 0;
            iEndIdx = iTemplateWidth - 1;
          }

          for (Int x = iStartIdx; x <= iEndIdx; x++)
          {
            pTrueDst[x*dstStrideTrue + y] = pDst[y*dstStride + x];
          }
        }
      }
      else if (eTempType == LEFT_NEIGHBOR)
      {
        for (Int y = 0; y < iTemplateHeight; y++)
        {
          for (Int x = 0; x < width; x++)
          {
            pTrueDst[x*dstStrideTrue + y] = pDst[y*dstStride + x];
          }
        }
      }
      else if (eTempType == ABOVE_NEIGHBOR)
      {
        for (Int y = 0; y < height; y++)
        {
          for (Int x = 0; x < iTemplateWidth; x++)
          {
            pTrueDst[x*dstStrideTrue + y] = pDst[y*dstStride + x];
          }
        }
      }
      else
      {
        assert(0);
      }
    }
  }
}

Void TComPrediction::xDIMDIntraPredFilteringMode34(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  Pel* pDst = rpDst;
  /*
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight > 0);
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight == 0);
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth == 0 && iTemplateHeight > 0);
  }
  */

  if (eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    for (Int y = 0, iDstStride2 = 0, iSrcStride2 = -1; y < iHeight; y++, iDstStride2 += iDstStride, iSrcStride2 += iSrcStride)
    {
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      if (y < iTemplateHeight)
      {
        if (iTemplateWidth == 1)
        {
          pDst[iDstStride2 + 1] = (12 * pDst[iDstStride2 + 1] + 4 * pSrc[iSrcStride2 + iSrcStride * 2] + 8) >> 4;
          pDst[iDstStride2 + 2] = (14 * pDst[iDstStride2 + 2] + 2 * pSrc[iSrcStride2 + iSrcStride * 3] + 8) >> 4;
          pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
        }
        else if (iTemplateWidth == 2)
        {
          pDst[iDstStride2 + 2] = (14 * pDst[iDstStride2 + 2] + 2 * pSrc[iSrcStride2 + iSrcStride * 3] + 8) >> 4;
          pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
        }
        else if (iTemplateWidth == 3)
        {
          pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
        }
      }
      else
#endif
        if (y >= iTemplateHeight)
        {
          pDst[iDstStride2] = (8 * pDst[iDstStride2] + 8 * pSrc[iSrcStride2 + iSrcStride] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
          if (iTemplateWidth > 1)
          {
            pDst[iDstStride2 + 1] = (12 * pDst[iDstStride2 + 1] + 4 * pSrc[iSrcStride2 + iSrcStride * 2] + 8) >> 4;
          }
          if (iTemplateWidth > 2)
          {
            pDst[iDstStride2 + 2] = (14 * pDst[iDstStride2 + 2] + 2 * pSrc[iSrcStride2 + iSrcStride * 3] + 8) >> 4;
          }
          if (iTemplateWidth > 3)
          {
            pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
          }
#endif
        }
    }
  }
  else if (eTempType == LEFT_NEIGHBOR)
  {
    for (Int y = 0, iDstStride2 = 0, iSrcStride2 = -1; y < iHeight; y++, iDstStride2 += iDstStride, iSrcStride2 += iSrcStride)
    {
      pDst[iDstStride2] = (8 * pDst[iDstStride2] + 8 * pSrc[iSrcStride2 + iSrcStride] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      if (iTemplateWidth > 1)
      {
        pDst[iDstStride2 + 1] = (12 * pDst[iDstStride2 + 1] + 4 * pSrc[iSrcStride2 + iSrcStride * 2] + 8) >> 4;
      }
      if (iTemplateWidth > 2)
      {
        pDst[iDstStride2 + 2] = (14 * pDst[iDstStride2 + 2] + 2 * pSrc[iSrcStride2 + iSrcStride * 3] + 8) >> 4;
      }
      if (iTemplateWidth > 3)
      {
        pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
      }
#endif
    }
  }
  else if (eTempType == ABOVE_NEIGHBOR)
  {
    //assert(iWidth >= 4);    
    for (Int y = 0, iDstStride2 = 0, iSrcStride2 = -1; y < iTemplateHeight; y++, iDstStride2 += iDstStride, iSrcStride2 += iSrcStride)
    {
      pDst[iDstStride2] = (8 * pDst[iDstStride2] + 8 * pSrc[iSrcStride2 + iSrcStride] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      pDst[iDstStride2 + 1] = (12 * pDst[iDstStride2 + 1] + 4 * pSrc[iSrcStride2 + iSrcStride * 2] + 8) >> 4;
      pDst[iDstStride2 + 2] = (14 * pDst[iDstStride2 + 2] + 2 * pSrc[iSrcStride2 + iSrcStride * 3] + 8) >> 4;
      pDst[iDstStride2 + 3] = (15 * pDst[iDstStride2 + 3] + pSrc[iSrcStride2 + iSrcStride * 4] + 8) >> 4;
#endif
    }
  }
  else
  {
    assert(0);
  }

  return;
}

Void TComPrediction::xDIMDIntraPredFilteringMode02(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  Pel* pDst = rpDst;
  /*
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight > 0);
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight == 0);
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth == 0 && iTemplateHeight > 0);
  }
  */

  if (eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    for (Int x = 0; x < iWidth; x++)
    {
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      if (x < iTemplateWidth)
      {
        if (iTemplateHeight == 1)
        {
          pDst[x + iDstStride] = (12 * pDst[x + iDstStride] + 4 * pSrc[x - iSrcStride + 2] + 8) >> 4;
          pDst[x + iDstStride * 2] = (14 * pDst[x + iDstStride * 2] + 2 * pSrc[x - iSrcStride + 3] + 8) >> 4;
          pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
        }
        else if (iTemplateHeight == 2)
        {
          pDst[x + iDstStride * 2] = (14 * pDst[x + iDstStride * 2] + 2 * pSrc[x - iSrcStride + 3] + 8) >> 4;
          pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
        }
        else if (iTemplateHeight == 3)
        {
          pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
        }
      }
      else
#endif
        if (x >= iTemplateWidth)
        {
          pDst[x] = (8 * pDst[x] + 8 * pSrc[x - iSrcStride + 1] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
          if (iTemplateHeight > 1)
          {
            pDst[x + iDstStride] = (12 * pDst[x + iDstStride] + 4 * pSrc[x - iSrcStride + 2] + 8) >> 4;
          }
          if (iTemplateHeight > 2)
          {
            pDst[x + iDstStride * 2] = (14 * pDst[x + iDstStride * 2] + 2 * pSrc[x - iSrcStride + 3] + 8) >> 4;
          }
          if (iTemplateHeight > 3)
          {
            pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
          }
#endif
        }
    }
  }
  else if (eTempType == LEFT_NEIGHBOR)
  {
    //assert(iHeight >= 4);
    for (Int x = 0; x < iTemplateWidth; x++)
    {
      pDst[x] = (8 * pDst[x] + 8 * pSrc[x - iSrcStride + 1] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      pDst[x + iDstStride] = (12 * pDst[x + iDstStride] + 4 * pSrc[x - iSrcStride + 2] + 8) >> 4;
      pDst[x + iDstStride * 2] = (14 * pDst[x + iDstStride * 2] + 2 * pSrc[x - iSrcStride + 3] + 8) >> 4;
      pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
#endif
    }
  }
  else if (eTempType == ABOVE_NEIGHBOR)
  {
    for (Int x = 0; x < iWidth; x++)
    {
      pDst[x] = (8 * pDst[x] + 8 * pSrc[x - iSrcStride + 1] + 8) >> 4;
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE
      if (iTemplateHeight > 1)
      {
        pDst[x + iDstStride] = (12 * pDst[x + iDstStride] + 4 * pSrc[x - iSrcStride + 2] + 8) >> 4;
      }
      if (iTemplateHeight > 2)
      {
        pDst[x + iDstStride * 2] = (14 * pDst[x + iDstStride * 2] + 2 * pSrc[x - iSrcStride + 3] + 8) >> 4;
      }
      if (iTemplateHeight > 3)
      {
        pDst[x + iDstStride * 3] = (15 * pDst[x + iDstStride * 3] + pSrc[x - iSrcStride + 4] + 8) >> 4;
      }
#endif
    }
  }
  else
  {
    assert(0);
  }

  return;
}

Void TComPrediction::xDIMDIntraPredFilteringModeDGL(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, UInt uiMode, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  Pel* pDst = rpDst;

#if DIMD_NUM_INTRA_DIR_INC
  const Int aucExtAngPredFilterCoef[16][3] = {
    { 12, 3, 1 }, ///< 18, 114
    { 12, 3, 1 }, ///< 17, 115
    { 12, 3, 1 }, ///< 16, 116
    { 12, 3, 1 }, ///< 15, 117
    { 12, 1, 3 }, ///< 14, 118
    { 12, 1, 3 }, ///< 13, 119
    { 12, 2, 2 }, ///< 12, 120
    { 12, 2, 2 }, ///< 11, 121
    { 12, 2, 2 }, ///< 10, 122
    { 12, 2, 2 }, ///<  9, 123 
    { 12, 3, 1 }, ///<  8, 124
    { 12, 3, 1 }, ///<  7, 125
    { 8, 6, 2 }, ///<  6, 126
    { 8, 6, 2 }, ///<  5, 127 
    { 8, 7, 1 }, ///<  4, 128 
    { 8, 7, 1 }, ///<  3, 129
  };
  const Int aucExtAngPredPosiOffset[16][2] = {
    { 2, 3 }, ///< 18, 114
    { 2, 3 }, ///< 17, 115 
    { 2, 3 }, ///< 16, 116
    { 2, 3 }, ///< 15, 117 
    { 1, 2 }, ///< 14, 118 
    { 1, 2 }, ///< 13, 119
    { 1, 2 }, ///< 12, 120
    { 1, 2 }, ///< 11, 121
    { 1, 2 }, ///< 10, 122
    { 1, 2 }, ///<  9, 123
    { 1, 2 }, ///<  8, 124
    { 1, 2 }, ///<  7, 125
    { 1, 2 }, ///<  6, 126
    { 1, 2 }, ///<  5, 127
    { 1, 2 }, ///<  4, 128
    { 1, 2 }, ///<  3, 129
  };
  assert((uiMode >= (EXT_VDIA_IDX - 16) && uiMode < EXT_VDIA_IDX) || (uiMode>2 && uiMode <= (2 + 16)));

  Bool bHorz = (uiMode < EXT_DIA_IDX);
  UInt deltaAng = bHorz ? ((2 + 16) - uiMode) : (uiMode - (EXT_VDIA_IDX - 16));
  const Int *offset = aucExtAngPredPosiOffset[deltaAng];
  const Int *filter = aucExtAngPredFilterCoef[deltaAng];
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
  const Int aucAngPredFilterCoef[8][3] = {
    { 12, 3, 1 },
    { 12, 3, 1 },
    { 12, 1, 3 },
    { 12, 2, 2 },
    { 12, 2, 2 },
    { 12, 3, 1 },
    { 8, 6, 2 },
    { 8, 7, 1 },
  };
  const Int aucAngPredPosiOffset[8][2] = {
    { 2, 3 },
    { 2, 3 },
    { 1, 2 },
    { 1, 2 },
    { 1, 2 },
    { 1, 2 },
    { 1, 2 },
    { 1, 2 },
  };
  //assert( ( uiMode>=(VDIA_IDX-8) && uiMode<VDIA_IDX ) || ( uiMode>2 && uiMode<=(2+8) ) );
#else
  const Int aucAngPredFilterCoef[4][3] = {
    { 12, 3, 1 },
    { 12, 1, 3 },
    { 12, 2, 2 },
    { 8, 6, 2 },
  };
  const Int aucAngPredPosiOffset[4][2] = {
    { 2, 3 },
    { 1, 2 },
    { 1, 2 },
    { 1, 2 },
  };
  assert((uiMode >= 30 && uiMode < 34) || (uiMode>2 && uiMode <= 6));
#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
  Bool bHorz = (uiMode < DIA_IDX);
  UInt deltaAng = bHorz ? ((2 + 8) - uiMode) : (uiMode - (VDIA_IDX - 8));
#else
  Bool bHorz = (uiMode < 18);
  UInt deltaAng = bHorz ? (6 - uiMode) : (uiMode - 30);
#endif
  const Int *offset = aucAngPredPosiOffset[deltaAng];
  const Int *filter = aucAngPredFilterCoef[deltaAng];
#endif

  if (bHorz)
  {
    Int iStartIdx, iEndIdx;
    if (eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      iStartIdx = iTemplateWidth;
      iEndIdx = iWidth - 1;
    }
    else if (eTempType == LEFT_NEIGHBOR)
    {
      iStartIdx = 0;
      iEndIdx = iTemplateWidth - 1;
    }
    else if (eTempType == ABOVE_NEIGHBOR)
    {
      iStartIdx = 0;
      iEndIdx = iWidth - 1;
    }
    else
    {
      assert(0);
    }
    for (Int x = iStartIdx; x <= iEndIdx; x++)
    {
      pDst[x] = (filter[0] * pDst[x]
        + filter[1] * pSrc[x - iSrcStride + offset[0]]
        + filter[2] * pSrc[x - iSrcStride + offset[1]] + 8) >> 4;
    }
  }
  else
  {
    Int iStartIdx, iEndIdx;
    if (eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      iStartIdx = iTemplateHeight;
      iEndIdx = iHeight - 1;
    }
    else if (eTempType == LEFT_NEIGHBOR)
    {
      iStartIdx = 0;
      iEndIdx = iHeight - 1;
    }
    else if (eTempType == ABOVE_NEIGHBOR)
    {
      iStartIdx = 0;
      iEndIdx = iTemplateHeight - 1;
    }
    else
    {
      assert(0);
    }
    for (Int y = iStartIdx; y <= iEndIdx; y++)
    {
      pDst[y * iDstStride] = (filter[0] * pDst[y * iDstStride]
        + filter[1] * pSrc[(y + offset[0]) * iSrcStride - 1]
        + filter[2] * pSrc[(y + offset[1]) * iSrcStride - 1] + 8) >> 4;
    }
  }

  return;
}

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
Void TComPrediction::predDIMDIntraLumaPDPC(TComDataCU* pcCU, UInt uiDirMode, Pel* piPred, UInt uiStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight)
{
  /*
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight > 0);
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
  assert(iTemplateWidth > 0 && iTemplateHeight == 0);
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
  assert(iTemplateWidth == 0 && iTemplateHeight > 0);
  }
  assert(iWidth >= MIN_PU_SIZE && iHeight >= MIN_PU_SIZE && iWidth <= MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE && iHeight <= MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE);
  assert(isLuma(pcCU->getTextType()));
  */

  Int iLog2Width = 0, iLog2Height = 0;
  while ((iWidth >> iLog2Width) > 1) iLog2Width++;
  while ((iHeight >> iLog2Height) > 1) iLog2Height++;

  Int blkSizeGroup[2] = { std::min(4, 1 + iLog2Width - 2), std::min(4, 1 + iLog2Height - 2) };

  Pel *ptrSrc = getPredictorPtr(COMPONENT_Y, false);
  Pel *pDst = piPred;
  const Int channelsBitDepthForPrediction = pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);
  const Int iSrcStride = iWidth + iHeight + 1;
  const Int sw = iSrcStride;
  const Int iDoubleSize = iWidth + iHeight;

#if DIMD_NUM_INTRA_DIR_INC
  Int iSelMode = (uiDirMode > 1 ? 18 + ((Int(uiDirMode) - 66) >> 2) : uiDirMode);
  const Int *pdpcParam[2] = { g_pdpc_pred_param[blkSizeGroup[0]][iSelMode], g_pdpc_pred_param[blkSizeGroup[1]][iSelMode] };
  const Int *pPdpcPar = pdpcParam[iWidth < iHeight];
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
  Int   iSelMode = (uiDirMode > 1 ? 18 + ((Int(uiDirMode) - 34) >> 1) : uiDirMode);
  const Int *pdpcParam[2] = { g_pdpc_pred_param[blkSizeGroup[0]][iSelMode], g_pdpc_pred_param[blkSizeGroup[1]][iSelMode] };
  const Int *pPdpcPar = pdpcParam[iWidth < iHeight];
#else
  const Int *pdpcParam[2] = { g_pdpc_pred_param[blkSizeGroup[0]][uiDirMode], g_pdpc_pred_param[blkSizeGroup[1]][uiDirMode] };
  const Int *pPdpcPar = pdpcParam[iWidth < iHeight];
#endif
#endif

  Int * piRefVector = piTempRef + iDoubleSize;
  Int * piLowpRefer = piFiltRef + iDoubleSize;

  //unfiltered reference
  for (Int j = 0; j <= iDoubleSize; j++)
  {
    piRefVector[j] = ptrSrc[j];
  }

  for (Int i = 1; i <= iDoubleSize; i++)
  {
    piRefVector[-i] = ptrSrc[i*iSrcStride];
  }

  if (pPdpcPar[5] != 0)
  {
    // filter reference samples
    xReferenceFilter(iDoubleSize, pPdpcPar[4], pPdpcPar[5], piRefVector, piLowpRefer);
    for (Int j = 0; j <= iDoubleSize; j++)
    {
      ptrSrc[j] = piLowpRefer[j];
    }
    for (Int i = 1; i <= iDoubleSize; i++)
    {
      ptrSrc[i*iSrcStride] = piLowpRefer[-i];
    }
  }

  if (uiDirMode == PLANAR_IDX)
  {
    xPredDIMDIntraPlanar(ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight);
  }
  else
  {
    const ChannelType     channelType = toChannelType(COMPONENT_Y);
    const Bool            enableEdgeFilters = true;

#if VCEG_AZ07_INTRA_4TAP_FILTER
    const Bool            enable4TapFilter = pcCU->getSlice()->getSPS()->getUseIntra4TapFilter();
#endif

    xPredDIMDIntraAng(channelsBitDepthForPrediction, ptrSrc + sw + 1, sw, pDst, uiStride, iWidth, iHeight, channelType, uiDirMode, enableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
      , enable4TapFilter
#endif
      , eTempType, iTemplateWidth, iTemplateHeight);
  }

  // use unfiltered reference sample for weighted prediction
  if (pPdpcPar[5] != 0)
  {
    for (int j = 0; j <= iDoubleSize; j++)
    {
      ptrSrc[j] = piRefVector[j];
    }

    for (int i = 1; i <= iDoubleSize; i++)
    {
      ptrSrc[i*iSrcStride] = piRefVector[-i];
    }
  }

  Int scale = ((iLog2Width + iLog2Height) < 10) ? 0 : 1;
  Int ParShift = 6;
  Int ParScale = 1 << ParShift;
  Int ParOffset = 1 << (ParShift - 1);

  if (eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    for (Int row = 0; row < iHeight; row++)
    {
      Int iStartCol, iEndCol;
      if (row < iTemplateHeight)
      {
        iStartCol = iTemplateWidth;
        iEndCol = iWidth - 1;
      }
      else
      {
        iStartCol = 0;
        iEndCol = iTemplateWidth - 1;
      }

      Int pos = row * uiStride;
      Int shiftRow = row >> scale;
      Int Coeff_Top = pdpcParam[1][2] >> shiftRow;
      Int Coeff_offset = pdpcParam[1][3] >> shiftRow;

      for (Int col = iStartCol; col <= iEndCol; col++, pos++)
      {
        Int shiftCol = col >> scale;
        Int Coeff_Left = pdpcParam[0][0] >> shiftCol;
        Int Coeff_TopLeft = (pdpcParam[0][1] >> shiftCol) + Coeff_offset;
        Int Coeff_Cur = ParScale - Coeff_Left - Coeff_Top + Coeff_TopLeft;

        Int sampleVal = (Coeff_Left* piRefVector[-row - 1] + Coeff_Top * piRefVector[col + 1] - Coeff_TopLeft * piRefVector[0] + Coeff_Cur * pDst[pos] + ParOffset) >> ParShift;
        pDst[pos] = Clip3(0, ((1 << channelsBitDepthForPrediction) - 1), sampleVal);
      }
    }
  }
  else if (eTempType == ABOVE_NEIGHBOR)
  {
    for (Int row = 0; row < iTemplateHeight; row++)
    {
      Int pos = row * uiStride;
      Int shiftRow = row >> scale;
      Int Coeff_Top = pdpcParam[1][2] >> shiftRow;
      Int Coeff_offset = pdpcParam[1][3] >> shiftRow;

      for (Int col = 0; col < iWidth; col++, pos++)
      {
        Int shiftCol = col >> scale;
        Int Coeff_Left = pdpcParam[0][0] >> shiftCol;
        Int Coeff_TopLeft = (pdpcParam[0][1] >> shiftCol) + Coeff_offset;
        Int Coeff_Cur = ParScale - Coeff_Left - Coeff_Top + Coeff_TopLeft;

        Int sampleVal = (Coeff_Left* piRefVector[-row - 1] + Coeff_Top * piRefVector[col + 1] - Coeff_TopLeft * piRefVector[0] + Coeff_Cur * pDst[pos] + ParOffset) >> ParShift;
        pDst[pos] = Clip3(0, ((1 << channelsBitDepthForPrediction) - 1), sampleVal);
      }
    }
  }
  else if (eTempType == LEFT_NEIGHBOR)
  {
    for (Int row = 0; row < iHeight; row++)
    {
      Int pos = row * uiStride;
      Int shiftRow = row >> scale;
      Int Coeff_Top = pdpcParam[1][2] >> shiftRow;
      Int Coeff_offset = pdpcParam[1][3] >> shiftRow;

      for (Int col = 0; col < iTemplateWidth; col++, pos++)
      {
        Int shiftCol = col >> scale;
        Int Coeff_Left = pdpcParam[0][0] >> shiftCol;
        Int Coeff_TopLeft = (pdpcParam[0][1] >> shiftCol) + Coeff_offset;
        Int Coeff_Cur = ParScale - Coeff_Left - Coeff_Top + Coeff_TopLeft;

        Int sampleVal = (Coeff_Left* piRefVector[-row - 1] + Coeff_Top * piRefVector[col + 1] - Coeff_TopLeft * piRefVector[0] + Coeff_Cur * pDst[pos] + ParOffset) >> ParShift;
        pDst[pos] = Clip3(0, ((1 << channelsBitDepthForPrediction) - 1), sampleVal);
      }
    }
  }
  else
  {
    assert(0);
  }
}
#endif

UInt TComPrediction::calcTemplateSAD(Int bitDepth, Pel* pi0, Int iStride0, Pel* pi1, Int iStride1, Int iWidth, Int iHeight, Bool bBDClip)
{
  Pel* piOrg = pi0;
  Pel* piCur = pi1;
  UInt uiSum = 0;

  for (UInt iRow = 0; iRow < iHeight; iRow++)
  {
    for (UInt iCol = 0; iCol < iWidth; iCol++)
    {
      uiSum += abs(piOrg[iCol] - piCur[iCol]);
    }
    piOrg += iStride0;
    piCur += iStride1;
  }

  if (!bBDClip)
  {
    return uiSum;
  }
  else
  {
    return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(bitDepth - 8);
  }
}
#endif


Void fillReferenceSamples( const Int bitDepth, 
#if O0043_BEST_EFFORT_DECODING
                           const Int bitDepthDelta, 
#endif
                           const Pel* piRoiOrigin, 
                                 Pel* piIntraTemp,
                           const Bool* bNeighborFlags,
                           const Int iNumIntraNeighbor, 
                           const Int unitWidth, 
                           const Int unitHeight, 
                           const Int iAboveUnits, 
                           const Int iLeftUnits,
                           const UInt uiWidth, 
                           const UInt uiHeight, 
                           const Int iPicStride
#if MLIP_EXTEND_REF_LINE
    , UInt ExternRefNum
    , Pel* piSecondIntraTemp
#endif
)
{
  const Pel* piRoiTemp;

#if MLIP_EXTEND_REF_LINE
  const Pel* piSecondRoiTemp;
  UInt ExternUnit = ((2 * ExternRefNum + 1) / unitWidth) + 1;
  const Int iTotalSecondUnits = iAboveUnits + iLeftUnits + ExternUnit;
  Bool bSecondNeighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 10];
  for (Int i = 0; i < iLeftUnits; i++)
      bSecondNeighborFlags[i] = bNeighborFlags[i];
  for (Int i = iLeftUnits; i < iLeftUnits + ExternUnit; i++)
      bSecondNeighborFlags[i] = bNeighborFlags[iLeftUnits];
  for (Int i = iLeftUnits + ExternUnit; i < iTotalSecondUnits; i++)
      bSecondNeighborFlags[i] = bNeighborFlags[i - ExternUnit + 1];
#endif


  Int  i, j;
  Int  iDCValue = 1 << (bitDepth - 1);
  const Int iTotalUnits = iAboveUnits + iLeftUnits + 1; //+1 for top-left

  if (iNumIntraNeighbor == 0)
  {
    // Fill border with DC value
    for (i=0; i<uiWidth; i++)
    {
      piIntraTemp[i] = iDCValue;
    }
    for (i=1; i<uiHeight; i++)
    {
      piIntraTemp[i*uiWidth] = iDCValue;
    }
#if MLIP_EXTEND_REF_LINE
    for (i = 0; i < uiWidth + ExternRefNum; i++)
    {
        piSecondIntraTemp[i] = iDCValue;
    }
    for (i = 1; i < uiHeight + ExternRefNum; i++)
    {
        piSecondIntraTemp[i*(uiWidth + ExternRefNum)] = iDCValue;
    }
#endif
  }
  else if (iNumIntraNeighbor == iTotalUnits)
  {
    // Fill top-left border and top and top right with rec. samples
    piRoiTemp = piRoiOrigin - iPicStride - 1;

    for (i=0; i<uiWidth; i++)
    {
#if O0043_BEST_EFFORT_DECODING
      piIntraTemp[i] = piRoiTemp[i] << bitDepthDelta;
#else
      piIntraTemp[i] = piRoiTemp[i];
#endif
    }

    // Fill left and below left border with rec. samples
    piRoiTemp = piRoiOrigin - 1;

    for (i=1; i<uiHeight; i++)
    {
#if O0043_BEST_EFFORT_DECODING
      piIntraTemp[i*uiWidth] = (*(piRoiTemp)) << bitDepthDelta;
#else
      piIntraTemp[i*uiWidth] = *(piRoiTemp);
#endif
      piRoiTemp += iPicStride;
    }
#if MLIP_EXTEND_REF_LINE
    // Fill top-left border and top and top right with rec. samples
    piRoiTemp = piRoiOrigin - (ExternRefNum + 1) * iPicStride - (ExternRefNum + 1);

    for (i = 0; i < uiWidth + ExternRefNum; i++)
    {
#if O0043_BEST_EFFORT_DECODING
        piSecondIntraTemp[i] = piRoiTemp[i] << bitDepthDelta;
#else
        piSecondIntraTemp[i] = piRoiTemp[i];
#endif
    }

    // Fill left and below left border with rec. samples
    piRoiTemp = piRoiOrigin - (ExternRefNum + 1) - ExternRefNum * iPicStride;

    for (i = 1; i < uiHeight + ExternRefNum; i++)
    {
#if O0043_BEST_EFFORT_DECODING
        piSecondIntraTemp[i*uiWidth] = (*(piRoiTemp)) << bitDepthDelta;
#else
        piSecondIntraTemp[i*(uiWidth + ExternRefNum)] = *(piRoiTemp);
#endif
        piRoiTemp += iPicStride;
    }
#endif

  }
  else // reference samples are partially available
  {
    // all above units have "unitWidth" samples each, all left/below-left units have "unitHeight" samples each
    const Int  iTotalSamples = (iLeftUnits * unitHeight) + ((iAboveUnits + 1) * unitWidth);
    Pel  piIntraLine[5 * MAX_CU_SIZE];
    Pel  *piIntraLineTemp;
    const Bool *pbNeighborFlags;
#if MLIP_EXTEND_REF_LINE
    const Int iTotalSecondSamples = (iLeftUnits * unitHeight) + ((iAboveUnits + ExternUnit) * unitWidth);
    Pel  piIntraSecondLine[5 * MAX_CU_SIZE];
    Pel  *piIntraSecondLineTemp;
#endif

    // Initialize
    for (i=0; i<iTotalSamples; i++)
    {
      piIntraLine[i] = iDCValue;
    }
#if MLIP_EXTEND_REF_LINE
    for (i = 0; i < iTotalSecondSamples; i++)
    {
        piIntraSecondLine[i] = iDCValue;
    }
#endif
    // Fill top-left sample
    piRoiTemp = piRoiOrigin - iPicStride - 1;
    piIntraLineTemp = piIntraLine + (iLeftUnits * unitHeight);
    pbNeighborFlags = bNeighborFlags + iLeftUnits;
#if MLIP_EXTEND_REF_LINE
    piSecondRoiTemp = piRoiOrigin - (ExternRefNum + 1) * iPicStride - (ExternRefNum + 1);
    piIntraSecondLineTemp = piIntraSecondLine + (iLeftUnits * unitHeight);
#endif

    if (*pbNeighborFlags)
    {
#if O0043_BEST_EFFORT_DECODING
      Pel topLeftVal=piRoiTemp[0] << bitDepthDelta;
#else
      Pel topLeftVal=piRoiTemp[0];
#endif
      for (i=0; i<unitWidth; i++)
      {
        piIntraLineTemp[i] = topLeftVal;
      }
#if MLIP_EXTEND_REF_LINE
      for (i = 0; i < ExternRefNum; i++)
          piIntraSecondLineTemp[i] = piSecondRoiTemp[iPicStride * (ExternRefNum - i)];
      piIntraSecondLineTemp[ExternRefNum] = piSecondRoiTemp[0];
      for (i = 1; i < ExternRefNum; i++)
          piIntraSecondLineTemp[ExternRefNum + i] = piSecondRoiTemp[i];
      Pel topLeftSecondVal = piSecondRoiTemp[ExternRefNum];
      for (i = 2 * ExternRefNum; i < ExternUnit*unitWidth; i++)
      {
          piIntraSecondLineTemp[i] = topLeftSecondVal;
      }
#endif
    }

    // Fill left & below-left samples (downwards)
    piRoiTemp += iPicStride;
    piIntraLineTemp--;
    pbNeighborFlags--;
#if MLIP_EXTEND_REF_LINE
    piSecondRoiTemp += (ExternRefNum + 1) * iPicStride;
    piIntraSecondLineTemp--;
#endif

    for (j=0; j<iLeftUnits; j++)
    {
      if (*pbNeighborFlags)
      {
        for (i=0; i<unitHeight; i++)
        {
#if O0043_BEST_EFFORT_DECODING
          piIntraLineTemp[-i] = piRoiTemp[i*iPicStride] << bitDepthDelta;
#else
          piIntraLineTemp[-i] = piRoiTemp[i*iPicStride];
#endif
#if MLIP_EXTEND_REF_LINE
          piIntraSecondLineTemp[-i] = piSecondRoiTemp[i*iPicStride];
#endif

        }
      }
      piRoiTemp += unitHeight*iPicStride;
      piIntraLineTemp -= unitHeight;
#if MLIP_EXTEND_REF_LINE
      piSecondRoiTemp += unitHeight*iPicStride;
      piIntraSecondLineTemp -= unitHeight;
#endif
      pbNeighborFlags--;
    }

    // Fill above & above-right samples (left-to-right) (each unit has "unitWidth" samples)
    piRoiTemp = piRoiOrigin - iPicStride;
    // offset line buffer by iNumUints2*unitHeight (for left/below-left) + unitWidth (for above-left)
    piIntraLineTemp = piIntraLine + (iLeftUnits * unitHeight) + unitWidth;

#if MLIP_EXTEND_REF_LINE
    piSecondRoiTemp = piRoiOrigin - (ExternRefNum + 1) * iPicStride;
    piIntraSecondLineTemp = piIntraSecondLine + (iLeftUnits * unitHeight) + ExternUnit*unitWidth;
#endif

    pbNeighborFlags = bNeighborFlags + iLeftUnits + 1;
    for (j=0; j<iAboveUnits; j++)
    {
      if (*pbNeighborFlags)
      {
        for (i=0; i<unitWidth; i++)
        {
#if O0043_BEST_EFFORT_DECODING
          piIntraLineTemp[i] = piRoiTemp[i] << bitDepthDelta;
#else
          piIntraLineTemp[i] = piRoiTemp[i];
#if MLIP_EXTEND_REF_LINE
          piIntraSecondLineTemp[i] = piSecondRoiTemp[i];
#endif
#endif
        }
      }
      piRoiTemp += unitWidth;
      piIntraLineTemp += unitWidth;
#if MLIP_EXTEND_REF_LINE
      piSecondRoiTemp += unitWidth;
      piIntraSecondLineTemp += unitWidth;
#endif
      pbNeighborFlags++;
    }

    // Pad reference samples when necessary
    Int iCurrJnit = 0;
    Pel  *piIntraLineCur   = piIntraLine;
    const UInt piIntraLineTopRowOffset = iLeftUnits * (unitHeight - unitWidth);

#if MLIP_EXTEND_REF_LINE
    Int iSecondCurrJnit = 0;
    Pel  *piIntraSecondLineCur = piIntraSecondLine;
#endif


    if (!bNeighborFlags[0])
    {
      // very bottom unit of bottom-left; at least one unit will be valid.
      {
        Int   iNext = 1;
        while (iNext < iTotalUnits && !bNeighborFlags[iNext])
        {
          iNext++;
        }
        Pel *piIntraLineNext = piIntraLine + ((iNext < iLeftUnits) ? (iNext * unitHeight) : (piIntraLineTopRowOffset + (iNext * unitWidth)));
        const Pel refSample = *piIntraLineNext;
        // Pad unavailable samples with new value
        Int iNextOrTop = std::min<Int>(iNext, iLeftUnits);
        // fill left column
        while (iCurrJnit < iNextOrTop)
        {
          for (i=0; i<unitHeight; i++)
          {
            piIntraLineCur[i] = refSample;
          }
          piIntraLineCur += unitHeight;
          iCurrJnit++;
        }
        // fill top row
        while (iCurrJnit < iNext)
        {
          for (i=0; i<unitWidth; i++)
          {
            piIntraLineCur[i] = refSample;
          }
          piIntraLineCur += unitWidth;
          iCurrJnit++;
        }
      }
#if MLIP_EXTEND_REF_LINE
      {
          Int iSecondNext = 1;
          while (iSecondNext < iTotalSecondUnits && !bSecondNeighborFlags[iSecondNext])
          {
              iSecondNext++;
          }
          Pel *piIntraSecondLineNext = piIntraSecondLine + ((iSecondNext < iLeftUnits) ? (iSecondNext * unitHeight) : (piIntraLineTopRowOffset + (iSecondNext * unitWidth)));
          const Pel refSecondSample = *piIntraSecondLineNext;

          // Pad unavailable samples with new value
          Int iNextOrTop = std::min<Int>(iSecondNext, iLeftUnits);
          // fill left column
          while (iSecondCurrJnit < iNextOrTop)
          {
              for (i = 0; i < unitHeight; i++)
              {
                  piIntraSecondLineCur[i] = refSecondSample;
              }
              piIntraSecondLineCur += unitHeight;
              iSecondCurrJnit++;
          }
          // fill top row
          while (iSecondCurrJnit < iSecondNext)
          {
              for (i = 0; i < unitWidth; i++)
              {
                  piIntraSecondLineCur[i] = refSecondSample;
              }
              piIntraSecondLineCur += unitWidth;
              iSecondCurrJnit++;
          }
      }
#endif

    }

    // pad all other reference samples.
    while (iCurrJnit < iTotalUnits)
    {
      if (!bNeighborFlags[iCurrJnit]) // samples not available
      {
        {
          const Int numSamplesInCurrUnit = (iCurrJnit >= iLeftUnits) ? unitWidth : unitHeight;
          const Pel refSample = *(piIntraLineCur-1);
          for (i=0; i<numSamplesInCurrUnit; i++)
          {
            piIntraLineCur[i] = refSample;
          }
          piIntraLineCur += numSamplesInCurrUnit;
          iCurrJnit++;
        }
      }
      else
      {
        piIntraLineCur += (iCurrJnit >= iLeftUnits) ? unitWidth : unitHeight;
        iCurrJnit++;
      }
    }
#if MLIP_EXTEND_REF_LINE
    while (iSecondCurrJnit < iTotalSecondUnits)
    {
        if (!bSecondNeighborFlags[iSecondCurrJnit]) // samples not available
        {
            {
                const Int numSamplesInCurrUnit = (iSecondCurrJnit >= iLeftUnits) ? unitWidth : unitHeight;
                const Pel refSecondSample = *(piIntraSecondLineCur - 1);

                for (i = 0; i < numSamplesInCurrUnit; i++)
                {
                    piIntraSecondLineCur[i] = refSecondSample;
                }
                piIntraSecondLineCur += numSamplesInCurrUnit;
                iSecondCurrJnit++;
            }
        }
        else
        {
            piIntraSecondLineCur += (iSecondCurrJnit >= iLeftUnits) ? unitWidth : unitHeight;
            iSecondCurrJnit++;
        }
    }
#endif
    // Copy processed samples

    piIntraLineTemp = piIntraLine + uiHeight + unitWidth - 2;
    // top left, top and top right samples
    for (i=0; i<uiWidth; i++)
    {
      piIntraTemp[i] = piIntraLineTemp[i];
    }

    piIntraLineTemp = piIntraLine + uiHeight - 1;
    for (i=1; i<uiHeight; i++)
    {
      piIntraTemp[i*uiWidth] = piIntraLineTemp[-i];
    }
#if MLIP_EXTEND_REF_LINE
    piIntraSecondLineTemp = piIntraSecondLine + (iLeftUnits * unitHeight);
    for (i = 0; i < ExternRefNum; i++)
        piSecondIntraTemp[(ExternRefNum - i)*(uiWidth + ExternRefNum)] = piIntraSecondLineTemp[i];
    piSecondIntraTemp[0] = piIntraSecondLineTemp[ExternRefNum];
    for (i = 0; i < ExternRefNum; i++)
        piSecondIntraTemp[1 + i] = piIntraSecondLineTemp[ExternRefNum + i + 1];

    piIntraSecondLineTemp = piIntraSecondLine + (iLeftUnits * unitHeight) + ExternUnit * unitWidth - 1;
    // top left, top and top right samples
    for (i = 1; i < uiWidth; i++)
    {
        piSecondIntraTemp[i + ExternRefNum] = piIntraSecondLineTemp[i];
    }

    piIntraSecondLineTemp = piIntraSecondLine + uiHeight - 1;
    for (i = 1; i < uiHeight; i++)
    {
        piSecondIntraTemp[(i + ExternRefNum)*(uiWidth + ExternRefNum)] = piIntraSecondLineTemp[-i];
    }
#endif

  }
}

Bool TComPrediction::filteringIntraReferenceSamples(const ComponentID compID, UInt uiDirMode, UInt uiTuChWidth, UInt uiTuChHeight, const ChromaFormat chFmt, const Bool intraReferenceSmoothingDisabled
#if DIMD_NUM_INTRA_DIR_INC
  , TComDataCU* pcCU, UInt uiAbsPartIdx
#endif 
                                                   )
{
  Bool bFilter;

  if (!filterIntraReferenceSamples(toChannelType(compID), chFmt, intraReferenceSmoothingDisabled))
  {
    bFilter=false;
  }
  else
  {
#if COM16_C806_LARGE_CTU || JVET_C0024_QTBT
    assert(uiTuChWidth>=MIN_PU_SIZE && uiTuChHeight>=MIN_PU_SIZE && uiTuChWidth<=MAX_CU_SIZE && uiTuChHeight<=MAX_CU_SIZE);
#else
    assert(uiTuChWidth>=4 && uiTuChHeight>=4 && uiTuChWidth<128 && uiTuChHeight<128);
#endif

#if DIMD_NUM_INTRA_DIR_INC
    if (pcCU->getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx))
    {
      uiDirMode = pcCU->getExtIntraDir(toChannelType(compID), uiAbsPartIdx);
      assert(uiDirMode >= PLANAR_IDX && uiDirMode <= EXT_VDIA_IDX);
    }
    else
    {
      assert(uiDirMode >= PLANAR_IDX && uiDirMode < NUM_INTRA_MODE);
    }
#endif

    if (uiDirMode == DC_IDX)
    {
      bFilter=false; //no smoothing for DC or LM chroma
    }
    else
    {
#if DIMD_NUM_INTRA_DIR_INC
      if (pcCU->getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx))
      {
        Int diff = min<Int>(abs((Int)uiDirMode - EXT_HOR_IDX), abs((Int)uiDirMode - EXT_VER_IDX));
        UInt log2Size = ((g_aucConvertToBit[uiTuChWidth] + g_aucConvertToBit[uiTuChHeight]) >> 1) + MIN_CU_LOG2;
        UInt sizeIndex = log2Size - 1;
        assert(sizeIndex < MAX_INTRA_FILTER_DEPTHS);
        bFilter = diff > m_aucExtIntraFilter[toChannelType(compID)][sizeIndex];
      }
      else
      {
#endif
      Int diff = min<Int>(abs((Int) uiDirMode - HOR_IDX), abs((Int)uiDirMode - VER_IDX));
#if JVET_C0024_QTBT
      UInt log2Size=((g_aucConvertToBit[uiTuChWidth]+g_aucConvertToBit[uiTuChHeight])>>1) + MIN_CU_LOG2;
      UInt sizeIndex = log2Size - 1;
#else
      UInt sizeIndex=g_aucConvertToBit[uiTuChWidth];
#endif
      assert(sizeIndex < MAX_INTRA_FILTER_DEPTHS);
      bFilter = diff > m_aucIntraFilter[toChannelType(compID)][sizeIndex];
#if DIMD_NUM_INTRA_DIR_INC
    }
#endif
    }
  }
  return bFilter;
}

Bool isAboveLeftAvailable( TComDataCU* pcCU, UInt uiPartIdxLT )
{
  Bool bAboveLeftFlag;
  UInt uiPartAboveLeft;
  TComDataCU* pcCUAboveLeft = pcCU->getPUAboveLeft( uiPartAboveLeft, uiPartIdxLT );
  if(pcCU->getSlice()->getPPS()->getConstrainedIntraPred())
  {
    bAboveLeftFlag = ( pcCUAboveLeft && pcCUAboveLeft->isIntra( uiPartAboveLeft ) );
  }
  else
  {
    bAboveLeftFlag = (pcCUAboveLeft ? true : false);
  }
  return bAboveLeftFlag;
}

Int isAboveAvailable( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxRT, Bool *bValidFlags )
{
  const UInt uiRasterPartBegin = g_auiZscanToRaster[uiPartIdxLT];
  const UInt uiRasterPartEnd = g_auiZscanToRaster[uiPartIdxRT]+1;
  const UInt uiIdxStep = 1;
  Bool *pbValidFlags = bValidFlags;
  Int iNumIntra = 0;

  for ( UInt uiRasterPart = uiRasterPartBegin; uiRasterPart < uiRasterPartEnd; uiRasterPart += uiIdxStep )
  {
    UInt uiPartAbove;
    TComDataCU* pcCUAbove = pcCU->getPUAbove( uiPartAbove, g_auiRasterToZscan[uiRasterPart] );
    if(pcCU->getSlice()->getPPS()->getConstrainedIntraPred())
    {
      if ( pcCUAbove && pcCUAbove->isIntra( uiPartAbove ) )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    else
    {
      if (pcCUAbove)
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    pbValidFlags++;
  }
  return iNumIntra;
}

Int isLeftAvailable( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxLB, Bool *bValidFlags )
{
  const UInt uiRasterPartBegin = g_auiZscanToRaster[uiPartIdxLT];
  const UInt uiRasterPartEnd = g_auiZscanToRaster[uiPartIdxLB]+1;
  const UInt uiIdxStep = pcCU->getPic()->getNumPartInCtuWidth();
  Bool *pbValidFlags = bValidFlags;
  Int iNumIntra = 0;

  for ( UInt uiRasterPart = uiRasterPartBegin; uiRasterPart < uiRasterPartEnd; uiRasterPart += uiIdxStep )
  {
    UInt uiPartLeft;
    TComDataCU* pcCULeft = pcCU->getPULeft( uiPartLeft, g_auiRasterToZscan[uiRasterPart] );
    if(pcCU->getSlice()->getPPS()->getConstrainedIntraPred())
    {
      if ( pcCULeft && pcCULeft->isIntra( uiPartLeft ) )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    else
    {
      if ( pcCULeft )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    pbValidFlags--; // opposite direction
  }

  return iNumIntra;
}

#if JVET_C0024_QTBT
Int isAboveRightAvailable( TComDataCU* pcCU, UInt uiPartIdxRT, UInt uiNumUnitsInPU, Bool *bValidFlags )
#else
Int isAboveRightAvailable( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxRT, Bool *bValidFlags )
#endif
{
#if !JVET_C0024_QTBT
  const UInt uiNumUnitsInPU = g_auiZscanToRaster[uiPartIdxRT] - g_auiZscanToRaster[uiPartIdxLT] + 1;
#endif
  Bool *pbValidFlags = bValidFlags;
  Int iNumIntra = 0;

  for ( UInt uiOffset = 1; uiOffset <= uiNumUnitsInPU; uiOffset++ )
  {
    UInt uiPartAboveRight;
    TComDataCU* pcCUAboveRight = pcCU->getPUAboveRight( uiPartAboveRight, uiPartIdxRT, uiOffset );
    if(pcCU->getSlice()->getPPS()->getConstrainedIntraPred())
    {
      if ( pcCUAboveRight && pcCUAboveRight->isIntra( uiPartAboveRight ) )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    else
    {
      if ( pcCUAboveRight )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    pbValidFlags++;
  }

  return iNumIntra;
}

#if JVET_C0024_QTBT
Int isBelowLeftAvailable( TComDataCU* pcCU, UInt uiPartIdxLB, UInt uiNumUnitsInPU, Bool *bValidFlags )
#else
Int isBelowLeftAvailable( TComDataCU* pcCU, UInt uiPartIdxLT, UInt uiPartIdxLB, Bool *bValidFlags )
#endif
{
#if !JVET_C0024_QTBT
  const UInt uiNumUnitsInPU = (g_auiZscanToRaster[uiPartIdxLB] - g_auiZscanToRaster[uiPartIdxLT]) / pcCU->getPic()->getNumPartInCtuWidth() + 1;
#endif
  Bool *pbValidFlags = bValidFlags;
  Int iNumIntra = 0;

  for ( UInt uiOffset = 1; uiOffset <= uiNumUnitsInPU; uiOffset++ )
  {
    UInt uiPartBelowLeft;
    TComDataCU* pcCUBelowLeft = pcCU->getPUBelowLeft( uiPartBelowLeft, uiPartIdxLB, uiOffset );
    if(pcCU->getSlice()->getPPS()->getConstrainedIntraPred())
    {
      if ( pcCUBelowLeft && pcCUBelowLeft->isIntra( uiPartBelowLeft ) )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    else
    {
      if ( pcCUBelowLeft )
      {
        iNumIntra++;
        *pbValidFlags = true;
      }
      else
      {
        *pbValidFlags = false;
      }
    }
    pbValidFlags--; // opposite direction
  }

  return iNumIntra;
}
#if RRC_INTRA_SHARPENING
#define FLT_MARGIN               (1 + 2 + 1)
#define IS_FP_SHIFT                20
Bool TComPrediction::initSharpPattern(Pel* piAdiBuf, Pel* piAdiBufSharp, Int iStride, Int iWidth, Int iHeight,
  TComTU &rTu, TComDataCU *pcCU, UInt uiAbsPartIdx,
  PelMap *tmpArrSource, Int tmpArrSourceSize,
  PelMap *tmpArrMap, Int tmpArrMapSize,
  PelMap *tmpArrMapTemp, Int tmpArrMap2Size,
  Int bitDepth
  )
{
#if !JVET_C0024_QTBT
  if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N)
  {
    assert(pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) != PLANAR_IDX);
    assert(pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx) != DC_IDX);
  }
#else
  assert(pcCU->getIntraDir(pcCU->getTextType(), uiAbsPartIdx) < LM_CHROMA_IDX);
  assert(pcCU->getIntraDir(pcCU->getTextType(), uiAbsPartIdx) != PLANAR_IDX);
  assert(pcCU->getIntraDir(pcCU->getTextType(), uiAbsPartIdx) != DC_IDX);
#endif

  {
    assert(tmpArrSourceSize >= iHeight + 2 * FLT_MARGIN && tmpArrSourceSize >= iWidth + 2 * FLT_MARGIN);
    assert(tmpArrMap2Size >= iHeight && tmpArrMap2Size >= iWidth);
    assert(tmpArrMapSize >= iHeight && tmpArrMapSize >= iWidth);
    assert(((MAX_CU_SIZE + 8) << 4) > iStride + 2 * FLT_MARGIN);
    assert(((MAX_CU_SIZE + 2) << 4) > iHeight + 2 * FLT_MARGIN);

    // fill linear patterns
    PelMap *pDst = tmpArrSource + FLT_MARGIN;
    for (int i = 0; i < iWidth; i++)
    {
      pDst[i] = piAdiBuf[i];
    }
    for (int i = 0; i < FLT_MARGIN; i++)
    {
      tmpArrSource[i] = pDst[0];
      pDst[iWidth + i] = pDst[iWidth - 1];
    }

    Int k = (Int)(-0.18f * (1 << IS_FP_SHIFT));

    calcEdgePattern(tmpArrSource, tmpArrMap, tmpArrMapTemp, iWidth, bitDepth, pcCU->getQP(uiAbsPartIdx));
    displaceVector(piAdiBuf, piAdiBufSharp, 1, tmpArrMap, 0, bitDepth, k, iWidth, true, 0);

    pDst = tmpArrSource + FLT_MARGIN;
    for (int i = 0; i < iHeight; i++)
    {
      pDst[i] = piAdiBuf[i*iStride];
    }
    for (int i = 0; i < FLT_MARGIN; i++)
    {
      tmpArrSource[i] = pDst[0];
      pDst[iHeight + i] = pDst[iHeight - 1];
    }
    calcEdgePattern(tmpArrSource, tmpArrMap, tmpArrMapTemp, iHeight, bitDepth, pcCU->getQP(uiAbsPartIdx));
    displaceVector(piAdiBuf, piAdiBufSharp, iStride, tmpArrMap, 0, bitDepth, k, iHeight, true, 0);
  }
  return true;
}

Void TComPrediction::calcEdgePattern(PelMap *pSamples, PelMap *pEdgeMap, PelMap *pTmpArr, int iLen, Int bitDepth, Int qp)
{
  int dx_half_len = 1;
  PelMap *pSrc = pSamples + dx_half_len;
  PelMap* pDst = pEdgeMap + dx_half_len;
  PelMap *pDstTmp = pTmpArr + dx_half_len;
  int len = iLen + 2 * FLT_MARGIN - 2 * dx_half_len;
  assert(dx_half_len <= FLT_MARGIN);

  int offset = dx_half_len;
  for (Int x = 0; x < len; x++)
  {
    PelMap hor = abs((pSrc[x + 1] >> 0) - (pSrc[x - 1] >> 0));
    pDst[x] = hor;
  }

  Int saturation_limit = 14;
  Int scale = bitDepth - 8;
  saturation_limit <<= scale;
  for (Int x = 0; x < len; x++)
  {
    pDst[x] = pDst[x] < saturation_limit ? pDst[x] : saturation_limit;
  }

  // blur
  int g_half_len = 2;
  pDst += g_half_len; pDstTmp += g_half_len;
  len -= 2 * g_half_len;
  for (int i = 0; i < g_half_len + 1; i++)
  {
    pDst[0 - i] = pDst[1];
    pDst[len - 1 + i] = pDst[len - 2];
  }

  offset = g_half_len;
  for (Int x = 0; x < len; x++)
  {
    pDstTmp[x]
      = ((pDst[x - 2] + pDst[x + 2])
      + ((pDst[x - 1] + pDst[x + 1]) << 2)
      + ((pDst[x + 0]) << 3) - pDst[x + 0]);
  }

  // 2nd derivations // TODO: merge with filtering
  pDst += dx_half_len; pDstTmp += dx_half_len;
  len -= 2 * dx_half_len;
  assert(len == iLen);
  offset = dx_half_len;
  for (Int x = 0; x < len; x++)
  {
    pEdgeMap[x] = pDstTmp[x + 1] - pDstTmp[x - 1];
  }
}

void TComPrediction::displaceVector(PelMap* pcPicSrc, Pel* pcPicDst, Int iStride, PelMap *pcMapX,
  int map_shift, int bitDepth, Int k, int block_width, int clipOut, int margin)
{
  Int widthSrc = block_width;
  PelMap* src = pcPicSrc;
  Int widthDst = block_width;
  Pel* dst = pcPicDst;
  PelMap* mapX = pcMapX;

  int maxVal = (1 << bitDepth) - 1;
  int midVal = 1 << (bitDepth - 1);

  int x_start = 0;
  if (x_start < 0) x_start = 0;

  Int scale = 12 * (1 << (bitDepth - 8));

  k = k / scale;

  for (int px = x_start; px < x_start + block_width; px++)
  {
    Int x = px * (1 << IS_FP_SHIFT);
    Int dx = mapX[px];

    x += k * dx;

    // check block borders

    if (x < (x_start - margin) * (1 << IS_FP_SHIFT)) 
      x = (x_start - margin) * (1 << IS_FP_SHIFT);
    if (x > (x_start + block_width - 1 + margin) * (1 << IS_FP_SHIFT))
      x = (x_start + block_width - 1 + margin) * (1 << IS_FP_SHIFT);

    // linear interpolation
    {
      Int pelVal = 0;
      Int intX = x >> IS_FP_SHIFT;
      Int fracDx = x & ((1 << IS_FP_SHIFT) - 1);

      {
        pelVal = (short)((((1 << IS_FP_SHIFT) - fracDx) * src[intX * iStride] + fracDx * src[(intX + 1) * iStride] + (1 << (IS_FP_SHIFT-1))) >> IS_FP_SHIFT);
      }

      if (clipOut)
      {
        pelVal = (pelVal < 0) ? 0 : pelVal;
        pelVal = (pelVal > maxVal) ? maxVal : pelVal;
      }
      dst[px*iStride] = pelVal;
    }
  }
}

#endif

//! \}
