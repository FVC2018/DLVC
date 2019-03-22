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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include "TComBilateralFilter.h"

const int TComBilateralFilter::SpatialSigmaValue = 62;

const int TComBilateralFilter::spatialSigmaBlockLengthOffsets[] = {20, 10, -10, 0, -10};
TComBilateralFilter* TComBilateralFilter::m_bilateralFilterInstance = NULL;

const UShort maxPosList[34] = {6, 12, 18, 23, 29, 35, 41, 46, 52, 58, 64, 69, 75, 81, 87, 92, 98, 104, 110, 115, 121, 127, 133, 138, 144, 150, 156, 161, 167, 173, 179, 184, 190, 196};

const Int maxDivToMulShift = 9;

TComBilateralFilter::TComBilateralFilter()
{
  int numQP = MAX_QP-18+1;
  // allocation
  m_bilateralFilterTable = new UShort*[numQP];
  for(int i = 0; i < numQP; i++)
  {
    m_bilateralFilterTable[i] = new UShort[maxPosList[i]+1];
  }

  // initialization
  for(int i = 0; i < numQP; i++)
  {
    for(int k = 0; k < (maxPosList[i]+1); k++)
    {
      m_bilateralFilterTable[i][k] = 0;
    }
  }
}

TComBilateralFilter::~TComBilateralFilter()
{
  int numQP = MAX_QP-18+1;
  for(int i = 0; i < numQP; ++i)
  {
    delete [] m_bilateralFilterTable[i];
  }
  delete [] m_bilateralFilterTable;
}

TComBilateralFilter* TComBilateralFilter::instance()
{
  if (m_bilateralFilterInstance == NULL)
    m_bilateralFilterInstance = new TComBilateralFilter();
  return m_bilateralFilterInstance;
}

Void TComBilateralFilter::createdivToMulLUTs()
{
  UInt one = 1 << BITS_PER_DIV_LUT_ENTRY; // 1 is represented by 2^14 (not 2^14 -1)
#if ITH_C2_DIV2MUL_LUT_MOD
  divToMulOneOverN[0] = one << maxDivToMulShift; // We can never divide by zero since the centerweight is non-zero, so we can set this value to something arbitrary.
  divToMulShift[0] = maxDivToMulShift;
#else
  divToMulOneOverN[0] = one; // We can never divide by zero since the centerweight is non-zero, so we can set this value to something arbitrary.
  divToMulShift[0] = 0;
#endif
  for (UInt n=1; n<BILATERAL_FILTER_MAX_DENOMINATOR_PLUS_ONE; n++)
  {
    UInt tryLUT = one / n;
    
    UInt tryShift = 0;
    // Make sure the LUT entry stored does not start with (binary) zeros.
    while(tryLUT <= one)
    {
      // This value of tryLUT
      divToMulOneOverN[n] = tryLUT;
      divToMulShift[n] = tryShift;
      
      tryShift++;
      tryLUT = (one << tryShift) / n;
    }
    
    // We may need to add 1 to the LUT entry in order to make 3/3, 4/4, 5/5, ... come out right.
    UInt adiv = divToMulOneOverN[n] * n / (one << divToMulShift[n]);
    if(adiv != 1)
      divToMulOneOverN[n]++;

#if ITH_C2_DIV2MUL_LUT_MOD
	divToMulOneOverN[n] = divToMulOneOverN[n] << (maxDivToMulShift - divToMulShift[n]);
	divToMulShift[n] = maxDivToMulShift;
#endif
  }
}

Void TComBilateralFilter::createBilateralFilterTable(Int qp)
{
  Int spatialSigmaValue;
  Int intensitySigmaValue = (qp - 17) * 50;
  Int sqrtSpatialSigmaMulTwo;
  Int sqrtIntensitySigmaMulTwo = 2 * intensitySigmaValue * intensitySigmaValue;
  int centerWeightTableSize = 5;

  spatialSigmaValue = SpatialSigmaValue;;
  for (Int i = 0; i < centerWeightTableSize; i++)
  {
    sqrtSpatialSigmaMulTwo = 2 * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[i]) * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[i]);

    // Calculate the multiplication factor that we will use to convert the first table (with the strongest filter) to one of the
    // tables that gives weaker filtering (such as when TU = 8 or 16 or when we have inter filtering). 
    Int sqrtSpatialSigmaMulTwoStrongestFiltering = 2 * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[0]) * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[0]);

    // multiplication factor equals exp(-1/stronger)/exp(-1/weaker)
    double centerWeightMultiplier = exp(-(10000.0 / sqrtSpatialSigmaMulTwoStrongestFiltering))/exp(-(10000.0 / sqrtSpatialSigmaMulTwo));
    m_bilateralCenterWeightTable[i] = (Int)(centerWeightMultiplier*65 + 0.5);
  }
  Int i = 0;
  sqrtSpatialSigmaMulTwo = 2 * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[i]) * (spatialSigmaValue + spatialSigmaBlockLengthOffsets[i]);
  for (Int j = 0; j < (maxPosList[qp-18]+1); j++)
  {        
    Int temp = j * 25;
    m_bilateralFilterTable[qp-18][j] = UShort(exp(-(10000.0 / sqrtSpatialSigmaMulTwo) - (temp * temp / (sqrtIntensitySigmaMulTwo * 1.0))) * 65 + 0.5);
  }
}
#if ITH_C3_SMT_BL_FILT==C_FXN
// UInt, UInt, Short*, Int, Int, Int, Pel*, UInt, UInt
Void TComBilateralFilter::smoothBlockBilateralFilter(TComDataCU* pcCU, UInt uiWidth, UInt uiHeight, Short block[], Int length, Int isInterBlock, Int qp)
{
  Int rightPixel, centerPixel;
  Int rightWeight, bottomWeight, centerWeight;
  Int sumWeights[MAX_CU_SIZE];
  Int sumDelta[MAX_CU_SIZE];
  Int blockLengthIndex;
  
  Int dIB, dIR;
    
  switch (length)
  {
    case 4:
      blockLengthIndex = 0;
      break;
    case 8:
      blockLengthIndex = 1;
      break;
    default:
      blockLengthIndex = 2;
      break;
  }
  
  UShort *lookupTablePtr;
  
  centerWeight = m_bilateralCenterWeightTable[blockLengthIndex + 3 * isInterBlock];
  
  Int theMaxPos = maxPosList[qp-18];
  lookupTablePtr = m_bilateralFilterTable[qp-18];

  // for each pixel in block
  
  // These are the types of pixels:
  //
  // A BB C
  //
  // D EE F
  // D EE F
  //
  // G HH I
  //
  // If the block is larger than 4x4, the E-part is larger.
  //
  // Filter types:
  //
  // AA  BBB   CC
  // A    B     C
  //
  // D    E     F
  // DD  EEE   FF
  // D    E     F
  //
  // G    H     I
  // GG  HHH   II
  // C uses a filter of type x
  Int currentPixelDeltaSum;
  Int currentPixelSumWeights;
  Int rightPixelDeltaSum;
  Int rightPixelSumWeights;
  Int rightWeightTimesdIR;
  Int bottomWeightTimesdIB;
  
  Int mySignIfNeg;
  Int mySign;
  
  Short *blockCurrentPixelPtr = block;
  Short *blockRightPixelPtr = blockCurrentPixelPtr+1;
  Short *blockNextLinePixelPtr = blockCurrentPixelPtr + uiWidth;
  Int *sumWeightsPtr = sumWeights;
  Int *sumDeltaPtr = sumDelta;
  
  // A pixel. uses filter type xx
  //                           x
  //
  // No information from previous row
  // No information from previous pixel
  // Storing information to next row
  // Storing information to next pixel
  
  // top left pixel; i = 0, j = 0;
  
  centerPixel = *(blockCurrentPixelPtr);
  rightPixel = *(blockRightPixelPtr++);
  dIR = rightPixel - centerPixel;
  dIB = *(blockNextLinePixelPtr++) - centerPixel;
  
  rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
  bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
  
  rightWeightTimesdIR = rightWeight*dIR;
  bottomWeightTimesdIB = bottomWeight*dIB;
  
  currentPixelSumWeights = centerWeight + rightWeight + bottomWeight;
  currentPixelDeltaSum = rightWeightTimesdIR + bottomWeightTimesdIB;
  
  rightPixelSumWeights = rightWeight; //next pixel to the right
  rightPixelDeltaSum = rightWeightTimesdIR;
  
  *(sumWeightsPtr++) = bottomWeight; //next pixel to the bottom
  *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
  
  mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
  mySign = 1 | mySignIfNeg;

  *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
  
  for (Int i = 1; i < (uiWidth-1); i++)
  {
    // B pixel. uses filter type xxx
    //                            x
    //
    // No information from previous row
    // Information reused from previous pixel
    // Storing information to next row
    // Storing information to next pixel
    
    centerPixel = rightPixel;
    rightPixel = *(blockRightPixelPtr++);
    dIR = rightPixel - centerPixel;
    dIB = *(blockNextLinePixelPtr++) - centerPixel;
    
    rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
    bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
    
    rightWeightTimesdIR = rightWeight*dIR;
    bottomWeightTimesdIB = bottomWeight*dIB;
    
    currentPixelSumWeights = centerWeight + rightPixelSumWeights + rightWeight + bottomWeight;
    currentPixelDeltaSum = rightWeightTimesdIR + bottomWeightTimesdIB - rightPixelDeltaSum;
    
    rightPixelSumWeights = rightWeight; //next pixel to the right
    rightPixelDeltaSum = rightWeightTimesdIR; //next pixel to the right
    
    *(sumWeightsPtr++) = bottomWeight; //next pixel to the bottom
    *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
    
    mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
    mySign = 1 | mySignIfNeg;

    *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
    
  }
  
  // C pixel. uses filter type xx
  //                            x
  //
  // No information from previous row
  // Information reused from previous pixel
  // Storing information to next row
  // No information to store to next pixel
  
  centerPixel = rightPixel;
  blockRightPixelPtr++;
  dIB = *(blockNextLinePixelPtr++) - centerPixel;
  
  bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
  bottomWeightTimesdIB = bottomWeight*dIB;
  
  currentPixelSumWeights = centerWeight + rightPixelSumWeights + bottomWeight;
  currentPixelDeltaSum = bottomWeightTimesdIB - rightPixelDeltaSum;
  
  *(sumWeightsPtr) = bottomWeight; //next pixel to the bottom
  *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
  
  mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
  mySign = 1 | mySignIfNeg;

  *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
  
  for (Int j = 1; j < (uiHeight-1); j++)
  {
    sumWeightsPtr = sumWeights;
    sumDeltaPtr = sumDelta;
    
    //                           x
    // D pixel. uses filter type xx
    //                           x
    //
    // Uses information from previous row
    // No information from previous pixel
    // Storing information to next row
    // Storing information to next pixel
    
    centerPixel = *(blockCurrentPixelPtr);
    rightPixel = *(blockRightPixelPtr++);
    dIR = rightPixel - centerPixel;
    dIB = *(blockNextLinePixelPtr++) - centerPixel;
    
    rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
    bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
    
    rightWeightTimesdIR = rightWeight*dIR;
    bottomWeightTimesdIB = bottomWeight*dIB;
    
    currentPixelSumWeights = centerWeight + *(sumWeightsPtr) + rightWeight + bottomWeight;
    currentPixelDeltaSum = rightWeightTimesdIR + bottomWeightTimesdIB - *(sumDeltaPtr);
    
    rightPixelSumWeights = rightWeight; //next pixel to the right
    rightPixelDeltaSum = rightWeightTimesdIR;
    
    *(sumWeightsPtr++) = bottomWeight; //next pixel to the bottom
    *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
    
    mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
    mySign = 1 | mySignIfNeg;

    *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
    
    for (Int i = 1; i < (uiWidth-1); i++)
    {
      //                            x
      // E pixel. uses filter type xxx
      //                            x
      //
      // Uses information from previous row
      // Uses information from previous pixel
      // Storing information to next row
      // No information to store to next pixel
      
      centerPixel = rightPixel;
      rightPixel = *(blockRightPixelPtr++);
      dIR = rightPixel - centerPixel;
      dIB = *(blockNextLinePixelPtr++) - centerPixel;
      
      rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
      bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
      
      rightWeightTimesdIR = rightWeight*dIR;
      bottomWeightTimesdIB = bottomWeight*dIB;
      
      currentPixelSumWeights = centerWeight + *(sumWeightsPtr) + rightPixelSumWeights + rightWeight + bottomWeight;
      currentPixelDeltaSum = rightWeightTimesdIR + bottomWeightTimesdIB - rightPixelDeltaSum - *(sumDeltaPtr);
      
      rightPixelSumWeights = rightWeight; //next pixel to the right
      rightPixelDeltaSum = rightWeightTimesdIR;
      
      *(sumWeightsPtr++) = bottomWeight; //next pixel to the bottom
      *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
      
      mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
      mySign = 1 | mySignIfNeg;

      *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
    }
    
    //                            x
    // F pixel. uses filter type xx
    //                            x
    //
    // Uses information from previous row
    // Uses information from previous pixel
    // Storing information to next row
    // Storing information to next pixel
    
    centerPixel = rightPixel;
    blockRightPixelPtr++;
    dIB = *(blockNextLinePixelPtr++) - centerPixel;
    
    bottomWeight = lookupTablePtr[std::min(theMaxPos, abs(dIB))];
    bottomWeightTimesdIB = bottomWeight*dIB;
    
    currentPixelSumWeights = centerWeight + *(sumWeightsPtr) + rightPixelSumWeights + bottomWeight;
    currentPixelDeltaSum = bottomWeightTimesdIB - rightPixelDeltaSum - *(sumDeltaPtr);
    
    *(sumWeightsPtr) = bottomWeight; //next pixel to the bottom
    *(sumDeltaPtr++) = bottomWeightTimesdIB; //next pixel to the bottom
    
    mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
    mySign = 1 | mySignIfNeg;

    *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
    
  }
  
  sumWeightsPtr = sumWeights;
  sumDeltaPtr = sumDelta;
  
  //                           x
  // G pixel. uses filter type xx
  //
  // Uses information from previous row
  // No information from previous pixel
  // No information to store to next row
  // Storing information to next pixel
  
  
  centerPixel = *(blockCurrentPixelPtr);
  rightPixel = *(blockRightPixelPtr++);
  dIR = rightPixel - centerPixel;
  
  rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
  rightWeightTimesdIR = rightWeight*dIR;
  
  currentPixelSumWeights = centerWeight + *(sumWeightsPtr++) + rightWeight;
  currentPixelDeltaSum = rightWeightTimesdIR - *(sumDeltaPtr++);
  
  rightPixelSumWeights = rightWeight; //next pixel to the right
  rightPixelDeltaSum = rightWeightTimesdIR;
  
  mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
  mySign = 1 | mySignIfNeg;

  *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
  
  for (Int i = 1; i < (uiWidth-1); i++)
  {
    //                            x
    // H pixel. uses filter type xxx
    //
    // Uses information from previous row
    // Uses information from previous pixel
    // No information to store to next row
    // Storing information to next pixel
    
    centerPixel = rightPixel;
    rightPixel = *(blockRightPixelPtr++);
    dIR = rightPixel - centerPixel;
    
    rightWeight = lookupTablePtr[std::min(theMaxPos, abs(dIR))];
    rightWeightTimesdIR = rightWeight*dIR;
    
    currentPixelSumWeights = centerWeight + *(sumWeightsPtr++) + rightWeight + rightPixelSumWeights;
    currentPixelDeltaSum = rightWeightTimesdIR - rightPixelDeltaSum - *(sumDeltaPtr++);
    
    rightPixelSumWeights = rightWeight; //next pixel to the right
    rightPixelDeltaSum = rightWeightTimesdIR;
    
    mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
    mySign = 1 | mySignIfNeg;

    *(blockCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
    
  }
  
  //                            x
  // I pixel. uses filter type xx
  //
  // Uses information from previous row
  // Uses information from previous pixel
  // No information to store to next row
  // No information to store to nex pixel
  
  centerPixel = rightPixel;
  
  currentPixelSumWeights = centerWeight + *(sumWeightsPtr) + rightPixelSumWeights;
  currentPixelDeltaSum = - rightPixelDeltaSum - *(sumDeltaPtr);
  
  mySignIfNeg = SIGN_IF_NEG(currentPixelDeltaSum);
  mySign = 1 | mySignIfNeg;

  *(blockCurrentPixelPtr) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights+mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights]) >> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));
  
}
#elif ITH_C3_SMT_BL_FILT==SIMD_FXN
Void TComBilateralFilter::smoothBlockBilateralFilter(TComDataCU* pcCU, UInt uiWidth, UInt uiHeight, Short dstRecon[], Int length, Int isInterBlock, Int qp, Pel *piReco, UInt uiSrcStride, UInt uiDstStride)
//Void TComBilateralFilter::smoothBlockBilateralFilter(UInt uiWidth, UInt uiHeight, Pel *piReco, UInt uiSrcStride, Short dstRecon[], UInt uiDstStride, Int length, Int isInterBlock, Int qp)
{
	// Filter types:
	//
	// AA  BBB   CC
	// A    B     C
	//
	// D    E     F
	// DD  EEE   FF
	// D    E     F
	//
	// G    H     I
	// GG  HHH   II

	//Pointers
	Pel *recon = piReco;
	Pel *recon_col_start = piReco;
	Short *dstRecon_col_start = dstRecon;

	//Registers
	__m128i PixelCenter, PixelBottom, PixelLeft, PixelRight;
	__m128i WeightsCenter, WeightsTop, WeightsBottom, WeightsLeft, WeightsRight;
	__m128i dLeft, dRight, dTop, dBottom;
	__m128i sumWeights, temp_sumWeights;
	__m128i maxPos_16x8;
	__m128i LeftWeightPointer, RightWeightPointer, BottomWeightPointer;
	__m128i zero_128, minus1_sign;
	__m128i WeightsTL_lo, WeightsTL_hi, dTL_lo, dTL_hi, WeightedTL_lo, WeightedTL_hi;
	__m128i abs_sumWeightedDelta_lo, abs_sumWeightedDelta_hi;
	__m128i WeightsBR_lo, WeightsBR_hi, dBR_lo, dBR_hi, WeightedBR_lo, WeightedBR_hi;
	__m128i sumWeightedDelta_lo, sumWeightedDelta_hi;
	__m128i sumWeights_lo, sumWeights_hi, mySignIfNeg_32bx4_lo, mySignIfNeg_32bx4_hi;
	__m128i sumWeights_with_neg_lo, sumWeights_with_neg_hi, sumWeights_with_neg_RightShift2_lo, sumWeights_with_neg_RightShift2_hi;
	__m128i sumWeightedDelta_sumWeights_lo, sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_lo, divToMulOneOverN_32bx4_hi;
	__m128i WeightsR_lo, WeightsR_hi, dR_lo, dR_hi, WeightedR_lo, WeightedR_hi;
	__m128i FilterWeight_lo, FilterWeight_hi, FilterResult;
	Int BlockLengthIndex;
	Short WeightCenter;
	UShort *lookupTablePtr;
	UShort BottomWeightsPointerBuffer[8];
	UShort RightWeightsPointerBuffer[8];
	UShort LeftWeightsPointerBuffer[8];
	UInt divToMulBuffer_lo[4];
	UInt divToMulBuffer_hi[4];
	Int theMaxPos = maxPosList[qp - 18];

	switch (length)
	{
	case 4:
		BlockLengthIndex = 0;
		break;
	case 8:
		BlockLengthIndex = 1;
		break;
	default:
		BlockLengthIndex = 2;
		break;
	}

	WeightCenter = m_bilateralCenterWeightTable[BlockLengthIndex + 3 * isInterBlock];

	lookupTablePtr = m_bilateralFilterTable[qp - 18];
	maxPos_16x8 = _mm_set1_epi16(theMaxPos);

	WeightsCenter = _mm_set1_epi16(WeightCenter);

	WeightsTop = _mm_setzero_si128();
	dTop = _mm_setzero_si128();
	zero_128 = _mm_setzero_si128();
	minus1_sign = _mm_set1_epi16(-1);

	PixelCenter = _mm_loadu_si128((__m128i *) (recon)); //C0 C1 C2 C3 C4 C5 C6 C7

	if (uiWidth == 4)
	{
		__m128i PixelCenter2, PixelBottom2, PixelLeft2, PixelRight2;

		//A, D pixel
		for (int row = 0; row < ((uiHeight / 2) - 1); row++)
		{
			//Loading of Registers
			PixelRight = _mm_srli_si128(PixelCenter, 2);						//R0 R1 R2 0
			PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 0 0 0 0 
			PixelLeft = _mm_slli_si128(PixelCenter, 2);							//0 L1 L2 L3

			PixelCenter2 = PixelBottom;

			PixelRight2 = _mm_srli_si128(PixelCenter2, 2);							//R0 R1 R2 0
			PixelBottom2 = _mm_loadu_si128((__m128i *) (recon + 2 * uiSrcStride));	//B0 B1 B2 B3 0 0 0 0 
			PixelLeft2 = _mm_slli_si128(PixelCenter2, 2);							//0 L1 L2 L3

			PixelCenter = _mm_unpacklo_epi64(PixelCenter, PixelCenter2);
			PixelLeft = _mm_unpacklo_epi64(PixelLeft, PixelLeft2);
			PixelRight = _mm_unpacklo_epi64(PixelRight, PixelRight2);
			PixelBottom = _mm_unpacklo_epi64(PixelBottom, PixelBottom2);

			dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//0; ......; L3-C3;
			dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; 0;
			dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B3-C3;
			dTop = _mm_unpacklo_epi64(dTop, _mm_sign_epi16(dBottom, minus1_sign));

			//Using Max position for finding minimum Weight
			RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));
			BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

			_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
			_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

			//Find Weights
			WeightsBottom = _mm_set_epi16(lookupTablePtr[BottomWeightsPointerBuffer[7]],
				lookupTablePtr[BottomWeightsPointerBuffer[6]],
				lookupTablePtr[BottomWeightsPointerBuffer[5]],
				lookupTablePtr[BottomWeightsPointerBuffer[4]],
				lookupTablePtr[BottomWeightsPointerBuffer[3]],
				lookupTablePtr[BottomWeightsPointerBuffer[2]],
				lookupTablePtr[BottomWeightsPointerBuffer[1]],
				lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

			WeightsRight = _mm_set_epi16(0,
				lookupTablePtr[RightWeightsPointerBuffer[6]],
				lookupTablePtr[RightWeightsPointerBuffer[5]],
				lookupTablePtr[RightWeightsPointerBuffer[4]],
				0,
				lookupTablePtr[RightWeightsPointerBuffer[2]],
				lookupTablePtr[RightWeightsPointerBuffer[1]],
				lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

															   // for n = (0 to 6): wLn+1 = wRn

			WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL2 wL3 wL4 wL5 wL6 wL7
			WeightsTop = _mm_unpacklo_epi64(WeightsTop, WeightsBottom);

			sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
			temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
			sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
			sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

			WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
			WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

			dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
			dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

			WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
			WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

			WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
			WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

			dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
			dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

			WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
			WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

			sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
			sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																			   //Result calculation
																			   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																			   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

			abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
			abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

			sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
			sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

			mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
			mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

			sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
			sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

			sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
			sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

			sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
			sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

			_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
			_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

			divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
				divToMulOneOverN[divToMulBuffer_lo[2]],
				divToMulOneOverN[divToMulBuffer_lo[1]],
				divToMulOneOverN[divToMulBuffer_lo[0]]);

			divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
				divToMulOneOverN[divToMulBuffer_hi[2]],
				divToMulOneOverN[divToMulBuffer_hi[1]],
				divToMulOneOverN[divToMulBuffer_hi[0]]);

			FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
			FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

			FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
			FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

			FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
			FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

			FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

			FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

			//store
			_mm_storel_epi64((__m128i *)dstRecon, FilterResult);
			dstRecon += uiDstStride;
			_mm_storel_epi64((__m128i *)dstRecon, _mm_srli_si128(FilterResult, 8));
			dstRecon += uiDstStride;

			dTop = _mm_sign_epi16(_mm_srli_si128(dBottom, 8), minus1_sign);
			WeightsTop = _mm_srli_si128(WeightsBottom, 8);
			PixelCenter = PixelBottom2;

			recon += 2 * uiSrcStride;
		}

		//Loop unrolling for last two rows
		//Loading of Registers
		PixelRight = _mm_srli_si128(PixelCenter, 2);						//R0 R1 R2 0
		PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 0 0 0 0 
		PixelLeft = _mm_slli_si128(PixelCenter, 2);							//0 L1 L2 L3

		PixelCenter2 = PixelBottom;

		PixelRight2 = _mm_srli_si128(PixelCenter2, 2);							//R0 R1 R2 0
		PixelLeft2 = _mm_slli_si128(PixelCenter2, 2);							//0 L1 L2 L3

		PixelCenter = _mm_unpacklo_epi64(PixelCenter, PixelCenter2);
		PixelLeft = _mm_unpacklo_epi64(PixelLeft, PixelLeft2);
		PixelRight = _mm_unpacklo_epi64(PixelRight, PixelRight2);

		dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//0; ......; L3-C3;
		dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; 0;
		dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B3-C3;
		dTop = _mm_unpacklo_epi64(dTop, _mm_sign_epi16(dBottom, minus1_sign));

		//Using Max position for finding minimum Weight
		RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));
		BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

		_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
		_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

		//Find Weights
		WeightsBottom = _mm_set_epi16(0,
			0,
			0,
			0,
			lookupTablePtr[BottomWeightsPointerBuffer[3]],
			lookupTablePtr[BottomWeightsPointerBuffer[2]],
			lookupTablePtr[BottomWeightsPointerBuffer[1]],
			lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

		WeightsRight = _mm_set_epi16(0,
			lookupTablePtr[RightWeightsPointerBuffer[6]],
			lookupTablePtr[RightWeightsPointerBuffer[5]],
			lookupTablePtr[RightWeightsPointerBuffer[4]],
			0,
			lookupTablePtr[RightWeightsPointerBuffer[2]],
			lookupTablePtr[RightWeightsPointerBuffer[1]],
			lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

														   // for n = (0 to 6): wLn+1 = wRn

		WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL2 wL3 wL4 wL5 wL6 wL7
		WeightsTop = _mm_unpacklo_epi64(WeightsTop, WeightsBottom);

		sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
		temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
		sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
		sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

		WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
		WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

		dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
		dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

		WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
		WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

		WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
		WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

		dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
		dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

		WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
		WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

		sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
		sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																		   //Result calculation
																		   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																		   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

		abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
		abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

		sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
		sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

		mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
		mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

		sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
		sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

		sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
		sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

		sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
		sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

		_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
		_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

		divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
			divToMulOneOverN[divToMulBuffer_lo[2]],
			divToMulOneOverN[divToMulBuffer_lo[1]],
			divToMulOneOverN[divToMulBuffer_lo[0]]);

		divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
			divToMulOneOverN[divToMulBuffer_hi[2]],
			divToMulOneOverN[divToMulBuffer_hi[1]],
			divToMulOneOverN[divToMulBuffer_hi[0]]);

		FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
		FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

		FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
		FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

		FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
		FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

		FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

		FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

		//store
		_mm_storel_epi64((__m128i *)dstRecon, FilterResult);
		dstRecon += uiDstStride;
		_mm_storel_epi64((__m128i *)dstRecon, _mm_srli_si128(FilterResult, 8));
		dstRecon += uiDstStride;

		return;
	}

	if (uiWidth == 8)
	{
		//A, D pixel
		for (int row = 0; row < (uiHeight - 1); row++)
		{
			//Loading of Registers
			PixelRight = _mm_srli_si128(PixelCenter, 2);						//R0 R1 R2 R3 R4 R5 R6 0
			PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 B4 B5 B6 B7
			PixelLeft = _mm_slli_si128(PixelCenter, 2);							//0 L1 L2 L3 L4 L5 L6 L7

			dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
			dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;
			dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B7-C7;

																//Using Max position for finding minimum Weight
			RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));
			BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

			_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
			_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

			//Find Weights
			WeightsBottom = _mm_set_epi16(lookupTablePtr[BottomWeightsPointerBuffer[7]],
				lookupTablePtr[BottomWeightsPointerBuffer[6]],
				lookupTablePtr[BottomWeightsPointerBuffer[5]],
				lookupTablePtr[BottomWeightsPointerBuffer[4]],
				lookupTablePtr[BottomWeightsPointerBuffer[3]],
				lookupTablePtr[BottomWeightsPointerBuffer[2]],
				lookupTablePtr[BottomWeightsPointerBuffer[1]],
				lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

			WeightsRight = _mm_set_epi16(0,
				lookupTablePtr[RightWeightsPointerBuffer[6]],
				lookupTablePtr[RightWeightsPointerBuffer[5]],
				lookupTablePtr[RightWeightsPointerBuffer[4]],
				lookupTablePtr[RightWeightsPointerBuffer[3]],
				lookupTablePtr[RightWeightsPointerBuffer[2]],
				lookupTablePtr[RightWeightsPointerBuffer[1]],
				lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

															   // for n = (0 to 6): wLn+1 = wRn

			WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL2 wL3 wL4 wL5 wL6 wL7

			sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
			temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
			sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
			sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

			WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
			WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

			dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
			dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

			WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
			WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

			WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
			WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

			dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
			dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

			WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
			WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

			sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
			sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																			   //Result calculation
																			   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																			   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

			abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
			abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

			sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
			sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

			mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
			mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

			sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
			sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

			sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
			sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

			sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
			sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

			_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
			_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

			divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
				divToMulOneOverN[divToMulBuffer_lo[2]],
				divToMulOneOverN[divToMulBuffer_lo[1]],
				divToMulOneOverN[divToMulBuffer_lo[0]]);

			divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
				divToMulOneOverN[divToMulBuffer_hi[2]],
				divToMulOneOverN[divToMulBuffer_hi[1]],
				divToMulOneOverN[divToMulBuffer_hi[0]]);

			FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
			FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

			FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
			FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

			FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
			FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

			FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

			FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

			//store
			_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
			dstRecon += uiDstStride;

			dTop = _mm_sign_epi16(dBottom, minus1_sign);
			WeightsTop = WeightsBottom;
			PixelCenter = PixelBottom;


			recon += uiSrcStride;
		}

		//loop unroll for G Pixel
		PixelRight = _mm_srli_si128(PixelCenter, 2);
		PixelLeft = _mm_slli_si128(PixelCenter, 2);

		dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);
		dRight = _mm_sub_epi16(PixelRight, PixelCenter);

		RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));

		_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

		//Find Weights
		WeightsRight = _mm_set_epi16(0,
			lookupTablePtr[RightWeightsPointerBuffer[6]],
			lookupTablePtr[RightWeightsPointerBuffer[5]],
			lookupTablePtr[RightWeightsPointerBuffer[4]],
			lookupTablePtr[RightWeightsPointerBuffer[3]],
			lookupTablePtr[RightWeightsPointerBuffer[2]],
			lookupTablePtr[RightWeightsPointerBuffer[1]],
			lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

		WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL1 wL2 wL3 wL4 wL5 wL6 wL7

		sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);
		temp_sumWeights = _mm_add_epi16(WeightsCenter, WeightsTop);
		sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);

		WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
		WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

		dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
		dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

		WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
		WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

		WeightsR_lo = _mm_unpacklo_epi16(zero_128, WeightsRight);	//0 wR0; 0 wR1; 0 wR2; 0 wR3;
		WeightsR_hi = _mm_unpackhi_epi16(zero_128, WeightsRight);	//0 wR4; 0 wR5; 0 wR6; 0 wR7;

		dR_lo = _mm_unpacklo_epi16(zero_128, dRight); //0 R0; 0 R1; 0 R2; 0 R3;
		dR_hi = _mm_unpackhi_epi16(zero_128, dRight); //0 R4; 0 R5; 0 R6; 0 R7;

		WeightedR_lo = _mm_madd_epi16(WeightsR_lo, dR_lo); //wR0*R0; wR1*R1; wR2*R2; wR7*R7;
		WeightedR_hi = _mm_madd_epi16(WeightsR_hi, dR_hi); //wR4*R4; wR5*R5; wR6*R6; wR7*R7;

		sumWeightedDelta_lo = _mm_add_epi32(WeightedR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
		sumWeightedDelta_hi = _mm_add_epi32(WeightedR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																		  //Result calculation
		abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
		abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

		sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
		sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

		mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
		mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

		sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
		sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

		sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
		sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

		sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
		sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

		_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
		_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

		divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
			divToMulOneOverN[divToMulBuffer_lo[2]],
			divToMulOneOverN[divToMulBuffer_lo[1]],
			divToMulOneOverN[divToMulBuffer_lo[0]]);

		divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
			divToMulOneOverN[divToMulBuffer_hi[2]],
			divToMulOneOverN[divToMulBuffer_hi[1]],
			divToMulOneOverN[divToMulBuffer_hi[0]]);


		FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
		FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

		FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
		FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

		FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
		FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

		FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

		FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

		//store
		_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
		return;
	}

	if (uiWidth > 8)
	{
		//A, D pixel
		for (int row = 0; row < (uiHeight - 1); row++)
		{
			//Loading of Registers
			PixelRight = _mm_loadu_si128((__m128i *) (recon + 1));				//R0 R1 R2 R3 R4 R5 R6 R7
			PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 B4 B5 B6 B7
			PixelLeft = _mm_slli_si128(PixelCenter, 2);							//L0 L1 L2 L3 L4 L5 L6 L7

			dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
			dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;
			dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B7-C7;

																//Using Max position for finding minimum Weight
			RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));
			BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

			_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
			_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

			//Find Weights
			WeightsBottom = _mm_set_epi16(lookupTablePtr[BottomWeightsPointerBuffer[7]],
				lookupTablePtr[BottomWeightsPointerBuffer[6]],
				lookupTablePtr[BottomWeightsPointerBuffer[5]],
				lookupTablePtr[BottomWeightsPointerBuffer[4]],
				lookupTablePtr[BottomWeightsPointerBuffer[3]],
				lookupTablePtr[BottomWeightsPointerBuffer[2]],
				lookupTablePtr[BottomWeightsPointerBuffer[1]],
				lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

			WeightsRight = _mm_set_epi16(lookupTablePtr[RightWeightsPointerBuffer[7]],
				lookupTablePtr[RightWeightsPointerBuffer[6]],
				lookupTablePtr[RightWeightsPointerBuffer[5]],
				lookupTablePtr[RightWeightsPointerBuffer[4]],
				lookupTablePtr[RightWeightsPointerBuffer[3]],
				lookupTablePtr[RightWeightsPointerBuffer[2]],
				lookupTablePtr[RightWeightsPointerBuffer[1]],
				lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

															   // for n = (0 to 6): wLn+1 = wRn

			WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL1 wL2 wL3 wL4 wL5 wL6 wL7

			sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
			temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
			sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
			sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

			WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
			WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

			dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
			dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

			WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
			WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

			WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
			WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

			dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
			dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

			WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
			WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

			sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
			sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																			   //Result calculation
																			   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																			   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

			abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
			abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

			sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
			sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

			mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
			mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

			sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
			sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

			sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
			sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

			sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
			sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

			_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
			_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

			divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
				divToMulOneOverN[divToMulBuffer_lo[2]],
				divToMulOneOverN[divToMulBuffer_lo[1]],
				divToMulOneOverN[divToMulBuffer_lo[0]]);

			divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
				divToMulOneOverN[divToMulBuffer_hi[2]],
				divToMulOneOverN[divToMulBuffer_hi[1]],
				divToMulOneOverN[divToMulBuffer_hi[0]]);

			FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
			FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

			FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
			FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

			FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
			FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

			FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

			FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

			//store
			_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
			dstRecon += uiDstStride;

			dTop = _mm_sign_epi16(dBottom, minus1_sign);
			WeightsTop = WeightsBottom;
			PixelCenter = PixelBottom;


			recon += uiSrcStride;
		}

		//loop unroll for G Pixel
		PixelRight = _mm_loadu_si128((__m128i *) (recon + 1));
		PixelLeft = _mm_slli_si128(PixelCenter, 2);

		dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);
		dRight = _mm_sub_epi16(PixelRight, PixelCenter);

		RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));

		_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

		//Find Weights
		WeightsRight = _mm_set_epi16(lookupTablePtr[RightWeightsPointerBuffer[7]],
			lookupTablePtr[RightWeightsPointerBuffer[6]],
			lookupTablePtr[RightWeightsPointerBuffer[5]],
			lookupTablePtr[RightWeightsPointerBuffer[4]],
			lookupTablePtr[RightWeightsPointerBuffer[3]],
			lookupTablePtr[RightWeightsPointerBuffer[2]],
			lookupTablePtr[RightWeightsPointerBuffer[1]],
			lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

		WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL1 wL2 wL3 wL4 wL5 wL6 wL7

		sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);
		temp_sumWeights = _mm_add_epi16(WeightsCenter, WeightsTop);
		sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);

		WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
		WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

		dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
		dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

		WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
		WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

		WeightsR_lo = _mm_unpacklo_epi16(zero_128, WeightsRight);	//0 wR0; 0 wR1; 0 wR2; 0 wR3;
		WeightsR_hi = _mm_unpackhi_epi16(zero_128, WeightsRight);	//0 wR4; 0 wR5; 0 wR6; 0 wR7;

		dR_lo = _mm_unpacklo_epi16(zero_128, dRight); //0 R0; 0 R1; 0 R2; 0 R3;
		dR_hi = _mm_unpackhi_epi16(zero_128, dRight); //0 R4; 0 R5; 0 R6; 0 R7;

		WeightedR_lo = _mm_madd_epi16(WeightsR_lo, dR_lo); //wR0*R0; wR1*R1; wR2*R2; wR7*R7;
		WeightedR_hi = _mm_madd_epi16(WeightsR_hi, dR_hi); //wR4*R4; wR5*R5; wR6*R6; wR7*R7;

		sumWeightedDelta_lo = _mm_add_epi32(WeightedR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
		sumWeightedDelta_hi = _mm_add_epi32(WeightedR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																		  //Result calculation
		abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
		abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

		sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
		sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

		mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
		mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

		sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
		sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

		sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
		sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

		sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
		sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

		_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
		_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

		divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
			divToMulOneOverN[divToMulBuffer_lo[2]],
			divToMulOneOverN[divToMulBuffer_lo[1]],
			divToMulOneOverN[divToMulBuffer_lo[0]]);

		divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
			divToMulOneOverN[divToMulBuffer_hi[2]],
			divToMulOneOverN[divToMulBuffer_hi[1]],
			divToMulOneOverN[divToMulBuffer_hi[0]]);

		FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
		FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

		FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
		FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

		FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
		FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

		FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

		FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

		//store
		_mm_storeu_si128((__m128i *)dstRecon, FilterResult);

		//B E H Pixel
		for (int col = 0; col < ((uiWidth / 8) - 2); col++)
		{
			WeightsTop = _mm_setzero_si128();
			dTop = _mm_setzero_si128();

			recon_col_start += 8;
			dstRecon_col_start += 8;
			recon = recon_col_start;
			dstRecon = dstRecon_col_start;

			PixelCenter = _mm_loadu_si128((__m128i *) (recon));

			for (int row = 0; row < (uiHeight - 1); row++)
			{
				//Loading of Registers
				PixelRight = _mm_loadu_si128((__m128i *) (recon + 1));				//R0 R1 R2 R3 R4 R5 R6 R7
				PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 B4 B5 B6 B7
				PixelLeft = _mm_loadu_si128((__m128i *) (recon - 1));				//L0 L1 L2 L3 L4 L5 L6 L7

				dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
				dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;
				dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B7-C7;

																	//Using Max position for finding minimum Weight
				LeftWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dLeft));
				RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));
				BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

				_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
				_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

				//Find Weights
				WeightsBottom = _mm_set_epi16(lookupTablePtr[BottomWeightsPointerBuffer[7]],
					lookupTablePtr[BottomWeightsPointerBuffer[6]],
					lookupTablePtr[BottomWeightsPointerBuffer[5]],
					lookupTablePtr[BottomWeightsPointerBuffer[4]],
					lookupTablePtr[BottomWeightsPointerBuffer[3]],
					lookupTablePtr[BottomWeightsPointerBuffer[2]],
					lookupTablePtr[BottomWeightsPointerBuffer[1]],
					lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

				WeightsRight = _mm_set_epi16(lookupTablePtr[RightWeightsPointerBuffer[7]],
					lookupTablePtr[RightWeightsPointerBuffer[6]],
					lookupTablePtr[RightWeightsPointerBuffer[5]],
					lookupTablePtr[RightWeightsPointerBuffer[4]],
					lookupTablePtr[RightWeightsPointerBuffer[3]],
					lookupTablePtr[RightWeightsPointerBuffer[2]],
					lookupTablePtr[RightWeightsPointerBuffer[1]],
					lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

																   // for n = (0 to 6): wLn+1 = wRn

				WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL1 wL2 wL3 wL4 wL5 wL6 wL7
				WeightsLeft = _mm_insert_epi16(WeightsLeft, lookupTablePtr[_mm_extract_epi16(LeftWeightPointer, 0)], 0); // wL0 wL1 wL2 wL3 wL4 wL5 wL6 wL7

				sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
				temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
				sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
				sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

				WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
				WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

				dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
				dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

				WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
				WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

				WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
				WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

				dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
				dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

				WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
				WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

				sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
				sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																				   //Result calculation
																				   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																				   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

				abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
				abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

				sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
				sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

				mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
				mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

				sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
				sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

				sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
				sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

				sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
				sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

				_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
				_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

				divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
					divToMulOneOverN[divToMulBuffer_lo[2]],
					divToMulOneOverN[divToMulBuffer_lo[1]],
					divToMulOneOverN[divToMulBuffer_lo[0]]);

				divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
					divToMulOneOverN[divToMulBuffer_hi[2]],
					divToMulOneOverN[divToMulBuffer_hi[1]],
					divToMulOneOverN[divToMulBuffer_hi[0]]);

				FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
				FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

				FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
				FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

				FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
				FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

				FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

				FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

				//store
				_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
				dstRecon += uiDstStride;

				dTop = _mm_sign_epi16(dBottom, minus1_sign);
				WeightsTop = WeightsBottom;
				PixelCenter = PixelBottom;

				recon += uiSrcStride;
			}

			//loop unroll for H Pixel
			PixelRight = _mm_loadu_si128((__m128i *) (recon + 1));	//R0 R1 R2 R3 R4 R5 R6 R7
			PixelLeft = _mm_loadu_si128((__m128i *) (recon - 1));	//L0 L1 L2 L3 L4 L5 L6 L7

			dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
			dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;

			LeftWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dLeft));
			RightWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dRight));

			_mm_storeu_si128((__m128i *) RightWeightsPointerBuffer, RightWeightPointer);

			//Find Weights
			WeightsRight = _mm_set_epi16(lookupTablePtr[RightWeightsPointerBuffer[7]],
				lookupTablePtr[RightWeightsPointerBuffer[6]],
				lookupTablePtr[RightWeightsPointerBuffer[5]],
				lookupTablePtr[RightWeightsPointerBuffer[4]],
				lookupTablePtr[RightWeightsPointerBuffer[3]],
				lookupTablePtr[RightWeightsPointerBuffer[2]],
				lookupTablePtr[RightWeightsPointerBuffer[1]],
				lookupTablePtr[RightWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

			WeightsLeft = _mm_slli_si128(WeightsRight, 2); // 0 wL1 wL2 wL3 wL4 wL5 wL6 wL7
			WeightsLeft = _mm_insert_epi16(WeightsLeft, lookupTablePtr[_mm_extract_epi16(LeftWeightPointer, 0)], 0); // wL0 wL1 wL2 wL3 wL4 wL5 wL6 wL7

			sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);
			temp_sumWeights = _mm_add_epi16(WeightsCenter, WeightsTop);
			sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);

			WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
			WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

			dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
			dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

			WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
			WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

			WeightsR_lo = _mm_unpacklo_epi16(zero_128, WeightsRight);	//0 wR0; 0 wR1; 0 wR2; 0 wR3;
			WeightsR_hi = _mm_unpackhi_epi16(zero_128, WeightsRight);	//0 wR4; 0 wR5; 0 wR6; 0 wR7;

			dR_lo = _mm_unpacklo_epi16(zero_128, dRight); //0 R0; 0 R1; 0 R2; 0 R3;
			dR_hi = _mm_unpackhi_epi16(zero_128, dRight); //0 R4; 0 R5; 0 R6; 0 R7;

			WeightedR_lo = _mm_madd_epi16(WeightsR_lo, dR_lo); //wR0*R0; wR1*R1; wR2*R2; wR7*R7;
			WeightedR_hi = _mm_madd_epi16(WeightsR_hi, dR_hi); //wR4*R4; wR5*R5; wR6*R6; wR7*R7;

			sumWeightedDelta_lo = _mm_add_epi32(WeightedR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
			sumWeightedDelta_hi = _mm_add_epi32(WeightedR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																			  //Result calculation
			abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
			abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

			sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
			sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

			mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
			mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

			sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
			sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

			sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
			sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

			sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
			sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

			_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
			_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

			divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
				divToMulOneOverN[divToMulBuffer_lo[2]],
				divToMulOneOverN[divToMulBuffer_lo[1]],
				divToMulOneOverN[divToMulBuffer_lo[0]]);

			divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
				divToMulOneOverN[divToMulBuffer_hi[2]],
				divToMulOneOverN[divToMulBuffer_hi[1]],
				divToMulOneOverN[divToMulBuffer_hi[0]]);

			FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
			FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

			FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
			FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

			FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
			FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

			FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

			FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

			//store
			_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
		}

		WeightsTop = _mm_setzero_si128();
		dTop = _mm_setzero_si128();

		recon_col_start += 8;
		dstRecon_col_start += 8;
		recon = recon_col_start;
		dstRecon = dstRecon_col_start;

		PixelCenter = _mm_loadu_si128((__m128i *) (recon));

		//C, F pixel
		for (int row = 0; row < (uiHeight - 1); row++)
		{

			//Loading of Registers
			PixelBottom = _mm_loadu_si128((__m128i *) (recon + uiSrcStride));	//B0 B1 B2 B3 B4 B5 B6 B7
			PixelLeft = _mm_loadu_si128((__m128i *) (recon - 1));				//L0 L1 L2 L3 L4 L5 L6 L7
			PixelRight = _mm_srli_si128(PixelCenter, 2);						//R0 R1 R2 R3 R4 R5 R6 0

			dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
			dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;
			dBottom = _mm_sub_epi16(PixelBottom, PixelCenter);	//B0-C0; ......; B7-C7;

																//Using Max position for finding minimum Weight
			LeftWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dLeft));
			BottomWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dBottom));

			_mm_storeu_si128((__m128i *) BottomWeightsPointerBuffer, BottomWeightPointer);
			_mm_storeu_si128((__m128i *) LeftWeightsPointerBuffer, LeftWeightPointer);

			//Find Weights
			WeightsBottom = _mm_set_epi16(lookupTablePtr[BottomWeightsPointerBuffer[7]],
				lookupTablePtr[BottomWeightsPointerBuffer[6]],
				lookupTablePtr[BottomWeightsPointerBuffer[5]],
				lookupTablePtr[BottomWeightsPointerBuffer[4]],
				lookupTablePtr[BottomWeightsPointerBuffer[3]],
				lookupTablePtr[BottomWeightsPointerBuffer[2]],
				lookupTablePtr[BottomWeightsPointerBuffer[1]],
				lookupTablePtr[BottomWeightsPointerBuffer[0]]); //wB0 wB1 wB2 wB3 wB4 wB5 wB6 wB7 

			WeightsLeft = _mm_set_epi16(lookupTablePtr[LeftWeightsPointerBuffer[7]],
				lookupTablePtr[LeftWeightsPointerBuffer[6]],
				lookupTablePtr[LeftWeightsPointerBuffer[5]],
				lookupTablePtr[LeftWeightsPointerBuffer[4]],
				lookupTablePtr[LeftWeightsPointerBuffer[3]],
				lookupTablePtr[LeftWeightsPointerBuffer[2]],
				lookupTablePtr[LeftWeightsPointerBuffer[1]],
				lookupTablePtr[LeftWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

															  // for n = (0 to 6): wLn+1 = wRn

			WeightsRight = _mm_srli_si128(WeightsLeft, 2); // wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

			sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);		// wL0+wR0; wL1+wR1; wL2+wR2; wL3+wR3; wL4+wR4; wL5+wR5; wL6+wR6; wL7+wR7;
			temp_sumWeights = _mm_add_epi16(WeightsTop, WeightsBottom); // wT0+wB0; wT1+wB1; wT2+wB2; wT3+wB3; wT4+wB4; wT5+wB5; wT6+wB6; wT7+wB7;
			sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);	// sum(wL0,wR0,wT0,wB0);....;sum(wL7,wR7,wT7,wB7);
			sumWeights = _mm_add_epi16(sumWeights, WeightsCenter);		// sum(wC0,wL0,wR0,wT0,wB0);....;sum(wC0,wL7,wR7,wT7,wB7);

			WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
			WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

			dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
			dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

			WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
			WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

			WeightsBR_lo = _mm_unpacklo_epi16(WeightsBottom, WeightsRight);	//wB0 wR0; wB1 wR1; wB2 wR2; wB3 wR3;
			WeightsBR_hi = _mm_unpackhi_epi16(WeightsBottom, WeightsRight);	//wB4 wR4; wB5 wR5; wB6 wR6; wB7 wR7;

			dBR_lo = _mm_unpacklo_epi16(dBottom, dRight); //B0 R0; B1 R1; B2 R2; B3 R3;
			dBR_hi = _mm_unpackhi_epi16(dBottom, dRight); //B4 R4; B5 R5; B6 R6; B7 R7;

			WeightedBR_lo = _mm_madd_epi16(WeightsBR_lo, dBR_lo); //wB0*B0 + wR0*R0; wB1*B1 + wR1*R1; wB2*B2 + wR2*R2; wB3*B3 + wR7*R7;
			WeightedBR_hi = _mm_madd_epi16(WeightsBR_hi, dBR_hi); //wB4*B4 + wR4*R4; wB5*B5 + wR5*R5; wB6*B6 + wR6*R6; wB3*B3 + wR7*R7;

			sumWeightedDelta_lo = _mm_add_epi32(WeightedBR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
			sumWeightedDelta_hi = _mm_add_epi32(WeightedBR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																			   //Result calculation
																			   //*(dstReconCurrentPixelPtr++) = centerPixel + mySign*((((mySign*currentPixelDeltaSum + ((currentPixelSumWeights + mySignIfNeg) >> 1))*divToMulOneOverN[currentPixelSumWeights])
																			   //>> (BITS_PER_DIV_LUT_ENTRY + divToMulShift[currentPixelSumWeights])));

			abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
			abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

			sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
			sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

			mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
			mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

			sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
			sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

			sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
			sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

			sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
			sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

			_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
			_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

			divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
				divToMulOneOverN[divToMulBuffer_lo[2]],
				divToMulOneOverN[divToMulBuffer_lo[1]],
				divToMulOneOverN[divToMulBuffer_lo[0]]);

			divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
				divToMulOneOverN[divToMulBuffer_hi[2]],
				divToMulOneOverN[divToMulBuffer_hi[1]],
				divToMulOneOverN[divToMulBuffer_hi[0]]);

			FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
			FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

			FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
			FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

			FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
			FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

			FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

			FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

			//store
			_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
			dstRecon += uiDstStride;

			dTop = _mm_sign_epi16(dBottom, minus1_sign);
			WeightsTop = WeightsBottom;
			PixelCenter = PixelBottom;

			recon += uiSrcStride;
		}

		//loop unroll for I Pixel
		PixelRight = _mm_srli_si128(PixelCenter, 2);			//R0 R1 R2 R3 R4 R5 R6 0
		PixelLeft = _mm_loadu_si128((__m128i *) (recon - 1));	//L0 L1 L2 L3 L4 L5 L6 L7

		dLeft = _mm_sub_epi16(PixelLeft, PixelCenter);		//L0-C0; ......; L7-C7;
		dRight = _mm_sub_epi16(PixelRight, PixelCenter);	//R0-C0; ......; R7-C7;

		LeftWeightPointer = _mm_min_epu16(maxPos_16x8, _mm_abs_epi16(dLeft));

		_mm_storeu_si128((__m128i *) LeftWeightsPointerBuffer, LeftWeightPointer);

		//Find Weights
		WeightsLeft = _mm_set_epi16(lookupTablePtr[LeftWeightsPointerBuffer[7]],
			lookupTablePtr[LeftWeightsPointerBuffer[6]],
			lookupTablePtr[LeftWeightsPointerBuffer[5]],
			lookupTablePtr[LeftWeightsPointerBuffer[4]],
			lookupTablePtr[LeftWeightsPointerBuffer[3]],
			lookupTablePtr[LeftWeightsPointerBuffer[2]],
			lookupTablePtr[LeftWeightsPointerBuffer[1]],
			lookupTablePtr[LeftWeightsPointerBuffer[0]]); //wR0 wR1 wR2 wR3 wR4 wR5 wR6 wR7

														  // for n = (0 to 6): wLn+1 = wRn

		WeightsRight = _mm_srli_si128(WeightsLeft, 2); // wR0 wR1 wR2 wR3 wR4 wR5 wR6 0

		sumWeights = _mm_add_epi16(WeightsLeft, WeightsRight);
		temp_sumWeights = _mm_add_epi16(WeightsCenter, WeightsTop);
		sumWeights = _mm_add_epi16(sumWeights, temp_sumWeights);

		WeightsTL_lo = _mm_unpacklo_epi16(WeightsTop, WeightsLeft);	//wT0 wL0 ; wT1 wL1; wT2 wL2; wT3 wL3;
		WeightsTL_hi = _mm_unpackhi_epi16(WeightsTop, WeightsLeft);	//wT4 wL4 ; wT5 wL5; wT6 wL6; wT7 wL7;

		dTL_lo = _mm_unpacklo_epi16(dTop, dLeft);	//T0 L0; T1 L1; T2 L2; T3 L3;
		dTL_hi = _mm_unpackhi_epi16(dTop, dLeft);	//T4 L4; T5 L5; T6 L6; T7 L7;

		WeightedTL_lo = _mm_madd_epi16(WeightsTL_lo, dTL_lo); //wT0*T0 + wL0*L0; wT1*T1 + wL1*L1; wT2*T2 + wL2*L2; wT3*T3 + wL7*L7;
		WeightedTL_hi = _mm_madd_epi16(WeightsTL_hi, dTL_hi); //wT4*T4 + wL4*L4; wT5*T5 + wL5*L5; wT6*T6 + wL6*L6; wT3*T3 + wL7*L7;

		WeightsR_lo = _mm_unpacklo_epi16(zero_128, WeightsRight);	//0 wR0; 0 wR1; 0 wR2; 0 wR3;
		WeightsR_hi = _mm_unpackhi_epi16(zero_128, WeightsRight);	//0 wR4; 0 wR5; 0 wR6; 0 wR7;

		dR_lo = _mm_unpacklo_epi16(zero_128, dRight); //0 R0; 0 R1; 0 R2; 0 R3;
		dR_hi = _mm_unpackhi_epi16(zero_128, dRight); //0 R4; 0 R5; 0 R6; 0 R7;

		WeightedR_lo = _mm_madd_epi16(WeightsR_lo, dR_lo); //wR0*R0; wR1*R1; wR2*R2; wR7*R7;
		WeightedR_hi = _mm_madd_epi16(WeightsR_hi, dR_hi); //wR4*R4; wR5*R5; wR6*R6; wR7*R7;

		sumWeightedDelta_lo = _mm_add_epi32(WeightedR_lo, WeightedTL_lo); //sum(wB0*B0,wR0*R0,wT0*T0,wL0*L0);.......;sum(wB3*B3,wR3*R3,wT3*T3,wL3*L3);
		sumWeightedDelta_hi = _mm_add_epi32(WeightedR_hi, WeightedTL_hi); //sum(wB4*B4,wR4*R4,wT4*T4,wL4*L4);.......;sum(wB7*B7,wR7*R7,wT7*T7,wL7*L7);

																		  //Result calculation
		abs_sumWeightedDelta_lo = _mm_abs_epi32(sumWeightedDelta_lo); //abs(SumdW0); abs(SumdW1);abs(SumdW2);abs(SumdW3); 
		abs_sumWeightedDelta_hi = _mm_abs_epi32(sumWeightedDelta_hi); //abs(SumdW4); abs(SumdW5);abs(SumdW6);abs(SumdW7);

		sumWeights_lo = _mm_unpacklo_epi16(sumWeights, zero_128); // SumW0; SumW1; SumW2; SumW3;
		sumWeights_hi = _mm_unpackhi_epi16(sumWeights, zero_128); // SumW4; SumW5; SumW6; SumW7;

		mySignIfNeg_32bx4_lo = _mm_cmplt_epi32(sumWeightedDelta_lo, zero_128); //SIGN_IF_NEG(currentPixelDeltaSum)
		mySignIfNeg_32bx4_hi = _mm_cmplt_epi32(sumWeightedDelta_hi, zero_128);	//SIGN_IF_NEG(currentPixelDeltaSum)

		sumWeights_with_neg_lo = _mm_add_epi32(sumWeights_lo, mySignIfNeg_32bx4_lo); //(currentPixelSumWeights + mySignIfNeg)
		sumWeights_with_neg_hi = _mm_add_epi32(sumWeights_hi, mySignIfNeg_32bx4_hi); //(currentPixelSumWeights + mySignIfNeg)

		sumWeights_with_neg_RightShift2_lo = _mm_srli_epi32(sumWeights_with_neg_lo, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)
		sumWeights_with_neg_RightShift2_hi = _mm_srli_epi32(sumWeights_with_neg_hi, 1); //((currentPixelSumWeights + mySignIfNeg) >> 1)

		sumWeightedDelta_sumWeights_lo = _mm_add_epi32(abs_sumWeightedDelta_lo, sumWeights_with_neg_RightShift2_lo); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))
		sumWeightedDelta_sumWeights_hi = _mm_add_epi32(abs_sumWeightedDelta_hi, sumWeights_with_neg_RightShift2_hi); //(abs(currentPixelDeltaSum) + ((currentPixelSumWeights + mySignIfNeg) >> 1))

		_mm_storeu_si128((__m128i *)divToMulBuffer_lo, sumWeights_lo);
		_mm_storeu_si128((__m128i *)divToMulBuffer_hi, sumWeights_hi);

		divToMulOneOverN_32bx4_lo = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_lo[3]],
			divToMulOneOverN[divToMulBuffer_lo[2]],
			divToMulOneOverN[divToMulBuffer_lo[1]],
			divToMulOneOverN[divToMulBuffer_lo[0]]);

		divToMulOneOverN_32bx4_hi = _mm_set_epi32(divToMulOneOverN[divToMulBuffer_hi[3]],
			divToMulOneOverN[divToMulBuffer_hi[2]],
			divToMulOneOverN[divToMulBuffer_hi[1]],
			divToMulOneOverN[divToMulBuffer_hi[0]]);

		FilterWeight_lo = _mm_mullo_epi32(sumWeightedDelta_sumWeights_lo, divToMulOneOverN_32bx4_lo);
		FilterWeight_hi = _mm_mullo_epi32(sumWeightedDelta_sumWeights_hi, divToMulOneOverN_32bx4_hi);

		FilterWeight_lo = _mm_srli_epi32(FilterWeight_lo, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);
		FilterWeight_hi = _mm_srli_epi32(FilterWeight_hi, maxDivToMulShift + BITS_PER_DIV_LUT_ENTRY);

		FilterWeight_lo = _mm_sign_epi32(FilterWeight_lo, sumWeightedDelta_lo);
		FilterWeight_hi = _mm_sign_epi32(FilterWeight_hi, sumWeightedDelta_hi);

		FilterResult = _mm_packs_epi32(FilterWeight_lo, FilterWeight_hi);

		FilterResult = _mm_add_epi16(PixelCenter, FilterResult);

		//store
		_mm_storeu_si128((__m128i *)dstRecon, FilterResult);
	}
}
#endif
Void TComBilateralFilter::bilateralFilterIntra(TComDataCU *pcCU, UInt uiWidth, UInt uiHeight, Pel *piReco, UInt uiStride, Int qp)
{
  UInt uiMinSize = std::min(uiWidth, uiHeight);
  Short *tempblock = new Short[ uiWidth * uiHeight ];
  
  twoDimBufCopier
	  (
		  (UChar *)tempblock,
		  uiWidth * sizeof(Short),
		  (UChar *)piReco,
		  uiStride * sizeof(Short),
		  uiWidth * sizeof(Short),
		  uiHeight
		  );

  //for (UInt j = 0; j < uiHeight; j++)   
  //{
  //  memcpy(tempblock + j * uiWidth, piReco + j * uiStride, uiWidth * sizeof(Short));
  //}
  smoothBlockBilateralFilter(pcCU, uiWidth, uiHeight, tempblock, uiMinSize, 0, qp, piReco, uiStride, uiWidth);
  twoDimBufCopier
	  (
		  (UChar *)piReco,
		  uiStride * sizeof(Short),
		  (UChar *)tempblock,
		  uiWidth * sizeof(Short),
		  uiWidth * sizeof(Short),
		  uiHeight
		  );
  //for (UInt j = 0; j < uiHeight; j++)
  //{
  //  memcpy(piReco + j * uiStride, tempblock + j * uiWidth, uiWidth * sizeof(Short));
  //}
  delete[] tempblock;
}

Void TComBilateralFilter::bilateralFilterInter(TComDataCU *pcCU, UInt uiWidth, UInt uiHeight, Pel *piResi, UInt uiStrideRes, Pel *piPred, UInt uiPredStride, Pel *piReco, UInt uiRecStride, Int clipbd, Int qp)
{
  UInt uiMinSize = std::min(uiWidth, uiHeight);
  Short *tempblock = new Short[ uiWidth * uiHeight ];
  Pel *piPredTemp = piPred;
  Pel *piResiTemp = piResi;
  Pel *piRecoTemp = piReco;
  // Reco = Pred + Resi
  reconComputer
	  (
		  piResi,
		  uiStrideRes,
		  piPred,
		  uiPredStride,
		  piReco,
		  uiRecStride,
#if JVET_D0033_ADAPTIVE_CLIPPING
		  g_ClipParam.min(COMPONENT_Y),
		  g_ClipParam.max(COMPONENT_Y),
#else
		  0,
		  (1 << clipbd) - 1,
#endif
		  uiWidth,
		  uiHeight
		  );
//  for (UInt uiY = 0; uiY < uiHeight; ++uiY)
//  {
//    for (UInt uiX = 0; uiX < uiWidth; ++uiX)
//    {
//#if JVET_D0033_ADAPTIVE_CLIPPING
//      piReco[uiX] = ClipA(piPred[uiX] + piResi[uiX], COMPONENT_Y);
//#else
//      piReco[uiX] = ClipBD(piPred[uiX] + piResi[uiX], clipbd);
//#endif
//    }
//    piPred += uiPredStride;
//    piResi += uiStrideRes;
//    piReco += uiRecStride;
//  }

  piPred = piPredTemp;
  piResi = piResiTemp;
  piReco = piRecoTemp;

  // Reco' = filter(Reco)
  twoDimBufCopier
	  (
		  (UChar *)tempblock,
		  uiWidth * sizeof(Short),
		  (UChar *)piReco,
		  uiRecStride * sizeof(Short),
		  uiWidth * sizeof(Short),
		  uiHeight
		  );
  //for (UInt j = 0; j < uiHeight; j++)
  //{
  //  memcpy(tempblock + j * uiWidth, piReco + j * uiRecStride, uiWidth * sizeof(Short));
  //}
  smoothBlockBilateralFilter(pcCU, uiWidth, uiHeight, tempblock, uiMinSize, 1, qp, piReco, uiRecStride, uiWidth);

  twoDimBufCopier
	  (
		  (UChar *)piReco,
		  uiRecStride * sizeof(Short),
		  (UChar *)tempblock,
		  uiWidth * sizeof(Short),
		  uiWidth * sizeof(Short),
		  uiHeight
		  );

  //for (UInt j = 0; j < uiHeight; j++)
  //{
  //  memcpy(piReco + j * uiRecStride, tempblock + j * uiWidth, uiWidth * sizeof(Short));
  //}
  delete[] tempblock;

  // need to be performed if residual  is used
  // Resi' = Reco' - Pred
  residueComputer
	  (
		  piReco,
		  uiRecStride,
		  piPred,
		  uiPredStride,
		  piResi,
		  uiStrideRes,
		  uiWidth,
		  uiHeight
		  );
  //for (UInt uiY = 0; uiY < uiHeight; ++uiY)
  //{
  //  for (UInt uiX = 0; uiX < uiWidth; ++uiX)
  //  {
  //    piResi[uiX] = piReco[uiX] - piPred[uiX];
  //  }
  //  piPred += uiPredStride;
  //  piResi += uiStrideRes;
  //  piReco += uiRecStride;
  //}
}

