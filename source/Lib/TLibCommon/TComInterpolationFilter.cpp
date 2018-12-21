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

/**
 * \file
 * \brief Implementation of TComInterpolationFilter class
 */

// ====================================================================================================================
// Includes
// ====================================================================================================================

#include "TComRom.h"
#include "TComInterpolationFilter.h"
#include <assert.h>

#include "TComChromaFormat.h"

#if COM16_C806_SIMD_OPT
#include <emmintrin.h>  
#endif


//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Tables
// ====================================================================================================================

#if JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC
const TFilterCoeff TComInterpolationFilter::m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_LUMA] =
{
  {  0, 0,   0, 64,  0,   0,  0,  0 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  {  0, 1,  -3, 63,  4,  -2,  1,  0 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -1, 2,  -5, 62,  8,  -3,  1,  0 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 3,  -8, 60, 13,  -4,  1,  0 },
#endif
  { -1, 4, -10, 58, 17,  -5,  1,  0 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 4, -11, 52, 26,  -8,  3, -1 }, 
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1
  { -1, 3,  -9, 47, 31, -10,  4, -1 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 4, -11, 45, 34, -10,  4, -1 },
#endif
  { -1, 4, -11, 40, 40, -11,  4, -1 }, 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 4, -10, 34, 45, -11,  4, -1 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -1, 4, -10, 31, 47,  -9,  3, -1 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 3,  -8, 26, 52, -11,  4, -1 }, 
#endif
  {  0, 1,  -5, 17, 58, -10,  4, -1 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  {  0, 1,  -4, 13, 60,  -8,  3, -1 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  {  0, 1,  -3,  8, 62,  -5,  2, -1 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  {  0, 1,  -2,  4, 63,  -3,  1,  0 }
#endif
};

const TFilterCoeff TComInterpolationFilter::m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_CHROMA] =
{
  {  0, 64,  0,  0 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1, 63,  2,  0 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -2, 62,  4,  0 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -2, 60,  7,  -1 },
#endif
  { -2, 58, 10, -2 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -3, 57, 12, -2 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -4, 56, 14, -2 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -4, 55, 15, -2 },
#endif
  { -4, 54, 16, -2 }, 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -5, 53, 18, -2 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -6, 52, 20, -2 }, 
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -6, 49, 24, -3 },
#endif
  { -6, 46, 28, -4 }, 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -5, 44, 29, -4 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -4, 42, 30, -4 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -4, 39, 33, -4 },
#endif
  { -4, 36, 36, -4 }, 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -4, 33, 39, -4 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -4, 30, 42, -4 }, 
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -4, 29, 44, -5 },
#endif
  { -4, 28, 46, -6 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -3, 24, 49, -6 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -2, 20, 52, -6 }, 
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -2, 18, 53, -5 },
#endif
  { -2, 16, 54, -4 },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -2, 15, 55, -4 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  { -2, 14, 56, -4 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -2, 12, 57, -3 },
#endif
  { -2, 10, 58, -2 }, 
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { -1,  7, 60, -2 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 1 
  {  0,  4, 62, -2 },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >= 2 
  { 0,  2, 63, -1 },
#endif
};
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
// from SHVC upsampling filter
const Short TComInterpolationFilter::m_lumaFilter[8][NTAPS_LUMA] =
{
  {  0, 0,   0, 64,  0,   0,  0,  0 },
  { -1, 2,  -5, 62,  8,  -3,  1,  0 },
  { -1, 4, -10, 58, 17,  -5,  1,  0 },
  { -1, 3,  -9, 47, 31, -10,  4, -1 },
  { -1, 4, -11, 40, 40, -11,  4, -1 }, 
  { -1, 4, -10, 31, 47,  -9,  3, -1 },
  {  0, 1,  -5, 17, 58, -10,  4, -1 },
  {  0, 1,  -3,  8, 62,  -5,  2, -1 },
};

const Short TComInterpolationFilter::m_chromaFilter[16][NTAPS_CHROMA] =
{
  {  0, 64,  0,  0 },
  { -2, 62,  4,  0 },
  { -2, 58, 10, -2 },
  { -4, 56, 14, -2 },
  { -4, 54, 16, -2 }, 
  { -6, 52, 20, -2 }, 
  { -6, 46, 28, -4 }, 
  { -4, 42, 30, -4 },
  { -4, 36, 36, -4 }, 
  { -4, 30, 42, -4 }, 
  { -4, 28, 46, -6 },
  { -2, 20, 52, -6 }, 
  { -2, 16, 54, -4 },
  { -2, 14, 56, -4 },
  { -2, 10, 58, -2 }, 
  {  0,  4, 62, -2 }  
};
#else
const TFilterCoeff TComInterpolationFilter::m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA] =
{
  {  0, 0,   0, 64,  0,   0, 0,  0 },
  { -1, 4, -10, 58, 17,  -5, 1,  0 },
  { -1, 4, -11, 40, 40, -11, 4, -1 },
  {  0, 1,  -5, 17, 58, -10, 4, -1 }
};

const TFilterCoeff TComInterpolationFilter::m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA] =
{
  {  0, 64,  0,  0 },
  { -2, 58, 10, -2 },
  { -4, 54, 16, -2 },
  { -6, 46, 28, -4 },
  { -4, 36, 36, -4 },
  { -4, 28, 46, -6 },
  { -2, 16, 54, -4 },
  { -2, 10, 58, -2 }
};
#endif
#endif

#if VCEG_AZ07_FRUC_MERGE
#if JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC
const Short TComInterpolationFilter::m_lumaFilterBilinear[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_LUMA_FRUC] =
{
  { 64,  0, },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 60,  4, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=1 
  { 56,  8, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 52, 12, },
#endif
  { 48, 16, },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 44, 20, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=1 
  { 40, 24, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 36, 28, },
#endif
  { 32, 32, },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 28, 36, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=1 
  { 24, 40, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 20, 44, },
#endif
  { 16, 48, },
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  { 12, 52, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=1 
  {  8, 56, },
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE >=2 
  {  4, 60, },
#endif
};
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
const Short TComInterpolationFilter::m_lumaFilterBilinear[8][NTAPS_LUMA_FRUC] =
{
  { 64,  0, },
  { 56,  8, },
  { 48, 16, },
  { 40, 24, },
  { 32, 32, },
  { 24, 40, },
  { 16, 48, },
  {  8, 56, },
};
#else
const Short TComInterpolationFilter::m_lumaFilterBilinear[4][NTAPS_LUMA_FRUC] =
{
  { 64,  0, },
  { 48, 16, },
  { 32, 32, },
  { 16, 48, },
};
#endif
#endif
#endif

#if RRC_SELECTABLE_IF && SIF_GAUSS
#if !VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
const Short TComInterpolationFilter::m_lumaFilterGauss[4][NTAPS_LUMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  {  0,  0,  0, 64,  0,  0,  0,  0 },
#else
  {  0,  4, 15, 26, 15,  4,  0,  0 },
#endif
//{  0,  2, 11, 26, 20,  5,  0,  0 }, // old
  {  0,  2, 12, 25, 19,  6,  0,  0 }, // classic
#if BIF_BLUR_HALF
  {  0,  1,  8, 23, 23,  8,  1,  0 },// 1/2, Gauss
#else
  { -1, 4, -11, 40, 40, -11, 4, -1 }, // 1/2, original
#endif
//{  0,  0,  5, 20, 26, 11,  2,  0 },// old
  {  0,  0,  6, 19, 25, 12,  2,  0 }, // classic
};
const TFilterCoeff TComInterpolationFilter::m_chromaFilterGauss[8][NTAPS_CHROMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  {  0, 64,  0,  0 },
#else
//{ 17, 27, 17,  3 }, // non symmetric :(
  { 18, 28, 18,  0 },
#endif
  { 14, 27, 18,  5 },
  { 12, 26, 20,  6 },
  { 10, 25, 22,  7 },
#if BIF_BLUR_HALF
  { 9, 23, 23, 9 }, // 1/2, Gauss
#else
  { -4, 36, 36, -4 }, // 1/2, original
#endif
  {  7, 22, 25, 10 },
  {  6, 20, 26, 12 },
  {  5, 18, 27, 14 },

  /* old
  {  0, 64,  0,  0 },
  { 14, 28, 18,  4 }, // Gauss
  { -4, 54, 16, -2 }, // original
  { 10, 25, 22,  7 }, // Gauss
#if BIF_BLUR_HALF
  {  9, 23, 23,  9 }, // 1/2, Gauss
#else
  { -4, 36, 36, -4 }, // 1/2, original
#endif
  {  7, 22, 25, 10 }, // Gauss
  { -2, 16, 54, -4 }, // original
  {  4, 18, 28, 14 }  // Gauss
  */
};
#elif VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
const Short TComInterpolationFilter::m_lumaFilterGauss[8][NTAPS_LUMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  {  0,  0,  0, 64,  0,  0,  0,  0 },
#else
  {  0,  4, 15, 26, 15,  4,  0,  0 },
#endif
  {  0,  3, 14, 26, 17,  4,  0,  0 },
  {  0,  2, 12, 25, 19,  6,  0,  0 },
  {  0,  2, 10, 24, 21,  7,  0,  0 },
  {  0,  1,  8, 23, 23,  8,  1,  0 },
  {  0,  0,  7, 21, 24, 10,  2,  0 },
  {  0,  0,  6, 19, 25, 12,  2,  0 },
  {  0,  0,  4, 17, 26, 14,  3,  0 },

};
const TFilterCoeff TComInterpolationFilter::m_chromaFilterGauss[16][NTAPS_CHROMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  {  0, 64,  0,  0 },
#else
  { 18, 28, 18,  0 },
#endif
  { 15, 27, 18,  4 },
  { 14, 27, 18,  5 },
  { 13, 26, 20,  5 },
  { 12, 26, 20,  6 },
  { 11, 25, 21,  7 },
  { 10, 25, 22,  7 },
  {  9, 24, 23,  8 },
  {  9, 23, 23,  9 },
  {  8, 23, 24,  9 },
  {  7, 22, 25, 10 },
  {  7, 21, 25, 11 },
  {  6, 20, 26, 12 },
  {  5, 20, 26, 13 },
  {  5, 18, 27, 14 },
  {  4, 18, 27, 15 },

};
#elif VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 2
const Short TComInterpolationFilter::m_lumaFilterGauss[16][NTAPS_LUMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  { 0, 0,  0, 64,  0,  0, 0, 0 },  // 0
#else
  { 0, 4, 15, 26, 15,  4, 0, 0 },  // 0
#endif
  { 0, 3, 15, 26, 16,  4, 0, 0 },
  { 0, 3, 14, 26, 17,  4, 0, 0 },
  { 0, 2, 13, 26, 18,  5, 0, 0 },
  { 0, 2, 12, 25, 19,  6, 0, 0 },  // 4
  { 0, 2, 11, 24, 20,  6, 1, 0 },
  { 0, 2, 10, 24, 21,  7, 0, 0 },
  { 0, 1,  9, 23, 22,  8, 1, 0 },
  { 0, 1,  8, 23, 23,  8, 1, 0 },   // 8
  { 0, 1,  8, 22, 23,  9, 1, 0 },
  { 0, 0,  7, 21, 24, 10, 2, 0 },
  { 0, 1,  6, 20, 24, 11, 2, 0 },
  { 0, 0,  6, 19, 25, 12, 2, 0 },  // 12
  { 0, 0,  5, 18, 26, 13, 2, 0 },
  { 0, 0,  4, 17, 26, 14, 3, 0 },
  { 0, 0,  4, 16, 26, 15, 3, 0 },
};

const TFilterCoeff TComInterpolationFilter::m_chromaFilterGauss[32][NTAPS_CHROMA] =
{
#if !BIF_FILTER_ZERO_FRAC
  {  0, 64,  0,  0 },
#else
  { 18, 28, 18,  0 },
#endif
  { 16, 27, 17,  4 },
  { 15, 27, 18,  4 },
  { 15, 27, 18,  4 },
  { 14, 27, 18,  5 },
  { 14, 27, 18,  5 },
  { 13, 26, 20,  5 },
  { 13, 26, 20,  5 },
  { 12, 26, 20,  6 },
  { 12, 26, 20,  6 },
  { 11, 25, 21,  7 },
  { 11, 25, 21,  7 },
  { 10, 25, 22,  7 },
  { 10, 25, 22,  7 },
  { 9,  24, 23,  8 },
  { 9,  24, 23,  8 },
  { 9,  23, 23,  9 },
  { 8,  23, 24,  9 },
  { 8,  23, 24,  9 },
  { 7,  22, 25, 10 },
  { 7,  22, 25, 10 },
  { 7,  21, 25, 11 },
  { 7,  21, 25, 11 },
  { 6,  20, 26, 12 },
  { 6,  20, 26, 12 },
  { 5,  20, 26, 13 },
  { 5,  20, 26, 13 },
  { 5,  18, 27, 14 },
  { 5,  18, 27, 14 },
  { 4,  18, 27, 15 },
  { 4,  18, 27, 15 },
  { 4,  17, 27, 16 },
};
#else
#error not supported
#endif
#endif
#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
const Short TComInterpolationFilter::m_lumaFilterAffine[(NFRACS_LUMA_AFFINE)*NTAPS_LUMA] =
{
   0,  0,   0, 256,   0,   0,  0,  0,
   0,  1,  -3, 256,   4,  -2,  0,  0,
   0,  2,  -7, 255,   8,  -3,  1,  0,
  -1,  3, -10, 255,  12,  -4,  1,  0,
  -1,  4, -13, 254,  16,  -5,  2, -1,
  -1,  5, -16, 253,  20,  -7,  2,  0,
  -1,  6, -18, 251,  25,  -9,  3, -1,
  -2,  7, -21, 250,  29, -10,  4, -1,
  -2,  8, -23, 248,  34, -12,  4, -1,
  -2,  8, -25, 246,  38, -13,  5, -1,
  -2,  9, -27, 244,  43, -15,  5, -1,
  -2, 10, -30, 242,  48, -16,  6, -2,
  -2, 10, -31, 239,  52, -17,  5,  0,
  -2, 10, -32, 237,  57, -18,  6, -2,
  -2, 11, -34, 234,  63, -21,  7, -2,
  -2, 11, -35, 231,  68, -21,  6, -2,
  -3, 13, -38, 228,  74, -24,  9, -3,
  -2, 12, -38, 224,  78, -24,  7, -1,
  -3, 14, -40, 221,  84, -27, 10, -3,
  -2, 12, -39, 217,  88, -27,  8, -1,
  -3, 13, -40, 213,  94, -28,  9, -2,
  -3, 15, -43, 210, 100, -31, 11, -3,
  -3, 13, -41, 205, 104, -30,  9, -1,
  -3, 12, -41, 201, 110, -31,  9, -1,
  -3, 15, -43, 197, 116, -35, 12, -3,
  -3, 14, -43, 192, 121, -35, 12, -2,
  -2, 13, -42, 187, 126, -35, 10, -1,
  -3, 14, -43, 183, 132, -37, 12, -2,
  -2, 13, -42, 178, 137, -38, 12, -2,
  -3, 14, -42, 173, 143, -39, 12, -2,
  -3, 15, -43, 169, 148, -41, 14, -3,
  -3, 13, -41, 163, 153, -40, 13, -2,
  -3, 13, -40, 158, 158, -40, 13, -3,
  -2, 13, -40, 153, 163, -41, 13, -3,
  -3, 14, -41, 148, 169, -43, 15, -3,
  -2, 12, -39, 143, 173, -42, 14, -3,
  -2, 12, -38, 137, 178, -42, 13, -2,
  -2, 12, -37, 132, 183, -43, 14, -3,
  -1, 10, -35, 126, 187, -42, 13, -2,
  -2, 12, -35, 121, 192, -43, 14, -3,
  -3, 12, -35, 116, 197, -43, 15, -3,
  -1,  9, -31, 110, 201, -41, 12, -3,
  -1,  9, -30, 104, 205, -41, 13, -3,
  -3, 11, -31, 100, 210, -43, 15, -3,
  -2,  9, -28,  94, 213, -40, 13, -3,
  -1,  8, -27,  88, 217, -39, 12, -2,
  -3, 10, -27,  84, 221, -40, 14, -3,
  -1,  7, -24,  78, 224, -38, 12, -2,
  -3,  9, -24,  74, 228, -38, 13, -3,
  -2,  6, -21,  68, 231, -35, 11, -2,
  -2,  7, -21,  63, 234, -34, 11, -2,
  -2,  6, -18,  57, 237, -32, 10, -2,
   0,  5, -17,  52, 239, -31, 10, -2,
  -2,  6, -16,  48, 242, -30, 10, -2,
  -1,  5, -15,  43, 244, -27,  9, -2,
  -1,  5, -13,  38, 246, -25,  8, -2,
  -1,  4, -12,  34, 248, -23,  8, -2,
  -1,  4, -10,  29, 250, -21,  7, -2,
  -1,  3,  -9,  25, 251, -18,  6, -1,
   0,  2,  -7,  20, 253, -16,  5, -1,
  -1,  2,  -5,  16, 254, -13,  4, -1,
   0,  1,  -4,  12, 255, -10,  3, -1,
   0,  1,  -3,   8, 255,  -7,  2,  0,
   0,  0,  -2,   4, 256,  -3,  1,  0
};

const Short TComInterpolationFilter::m_chromaFilterAffine[(NFRACS_CHROMA_AFFINE) * NTAPS_CHROMA] = 
{
    0, 256,   0,   0,
   -2, 255,   4,  -1,
   11, 225,  23,  -3,
    8, 225,  26,  -3,
   -6, 248,  18,  -4,
   -6, 244,  24,  -6,
   -8, 243,  27,  -6,
  -11, 242,  32,  -7,
  -12, 240,  36,  -8,
  -13, 236,  42,  -9,
  -15, 235,  45,  -9,
  -15, 230,  51, -10,
  -15, 226,  56, -11,
  -15, 221,  62, -12,
  -17, 221,  65, -13,
  -19, 219,  69, -13,
  -20, 216,  74, -14,
  -17, 208,  80, -15,
  -19, 207,  84, -16,
  -19, 201,  90, -16,
  -19, 198,  94, -17,
  -20, 195,  98, -17,
  -19, 190, 103, -18,
  -19, 185, 108, -18,
  -20, 182, 112, -18,
  -21, 178, 117, -18,
  -23, 176, 122, -19,
  -21, 170, 126, -19,
  -22, 168, 131, -21,
  -21, 162, 135, -20,
  -19, 156, 139, -20,
  -20, 153, 144, -21,
  -20, 148, 148, -20,
  -21, 144, 153, -20,
  -20, 139, 156, -19,
  -20, 135, 162, -21,
  -21, 131, 168, -22,
  -19, 126, 170, -21,
  -19, 122, 176, -23,
  -18, 117, 178, -21,
  -18, 112, 182, -20,
  -18, 108, 185, -19,
  -18, 103, 190, -19,
  -17,  98, 195, -20,
  -17,  94, 198, -19,
  -16,  90, 201, -19,
  -16,  84, 207, -19,
  -15,  80, 208, -17,
  -14,  74, 216, -20,
  -13,  69, 219, -19,
  -13,  65, 221, -17,
  -12,  62, 221, -15,
  -11,  56, 226, -15,
  -10,  51, 230, -15,
   -9,  45, 235, -15,
   -9,  42, 236, -13,
   -8,  36, 240, -12,
   -7,  32, 242, -11,
   -6,  27, 243,  -8,
   -6,  24, 244,  -6,
   -4,  18, 248,  -6,
   -3,  26, 225,   8,
   -3,  23, 225,  11,
   -1,   4, 255,  -2
};
#endif
#endif


#if COM16_C806_SIMD_OPT

static const UChar shuffle_2Tap[16] = { 0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7, 8, 9 };
static const UChar shuffle_4Tap[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 8, 9 };

inline __m128i simdInterpolateLuma4( Short const *src , Int srcStride , __m128i *mmCoeff , const __m128i & mmOffset , Int shift )
{
  __m128i sumHi = _mm_setzero_si128();
  __m128i sumLo = _mm_setzero_si128();
  for( Int n = 0 ; n < 8 ; n++ )
  {
    __m128i mmPix = _mm_loadl_epi64( ( __m128i* )src );
    __m128i hi = _mm_mulhi_epi16( mmPix , mmCoeff[n] );
    __m128i lo = _mm_mullo_epi16( mmPix , mmCoeff[n] );
    sumHi = _mm_add_epi32( sumHi , _mm_unpackhi_epi16( lo , hi ) );
    sumLo = _mm_add_epi32( sumLo , _mm_unpacklo_epi16( lo , hi ) );
    src += srcStride;
  }
  sumHi = _mm_srai_epi32( _mm_add_epi32( sumHi , mmOffset ) , shift );
  sumLo = _mm_srai_epi32( _mm_add_epi32( sumLo , mmOffset ) , shift );
  return( _mm_packs_epi32( sumLo , sumHi ) );
}

inline __m128i simdInterpolateChroma4( Short const *src , Int srcStride , __m128i *mmCoeff , const __m128i & mmOffset , Int shift )
{
  __m128i sumHi = _mm_setzero_si128();
  __m128i sumLo = _mm_setzero_si128();
  for( Int n = 0 ; n < 4 ; n++ )
  {
    __m128i mmPix = _mm_loadl_epi64( ( __m128i* )src );
    __m128i hi = _mm_mulhi_epi16( mmPix , mmCoeff[n] );
    __m128i lo = _mm_mullo_epi16( mmPix , mmCoeff[n] );
    sumHi = _mm_add_epi32( sumHi , _mm_unpackhi_epi16( lo , hi ) );
    sumLo = _mm_add_epi32( sumLo , _mm_unpacklo_epi16( lo , hi ) );
    src += srcStride;
  }
  sumHi = _mm_srai_epi32( _mm_add_epi32( sumHi , mmOffset ) , shift );
  sumLo = _mm_srai_epi32( _mm_add_epi32( sumLo , mmOffset ) , shift );
  return( _mm_packs_epi32( sumLo , sumHi ) );
}

inline __m128i simdInterpolateLuma8( Short const *src , Int srcStride , __m128i *mmCoeff , const __m128i & mmOffset , Int shift )
{
  __m128i sumHi = _mm_setzero_si128();
  __m128i sumLo = _mm_setzero_si128();
  for( Int n = 0 ; n < 8 ; n++ )
  {
    __m128i mmPix = _mm_loadu_si128( ( __m128i* )src );
    __m128i hi = _mm_mulhi_epi16( mmPix , mmCoeff[n] );
    __m128i lo = _mm_mullo_epi16( mmPix , mmCoeff[n] );
    sumHi = _mm_add_epi32( sumHi , _mm_unpackhi_epi16( lo , hi ) );
    sumLo = _mm_add_epi32( sumLo , _mm_unpacklo_epi16( lo , hi ) );
    src += srcStride;
  }
  sumHi = _mm_srai_epi32( _mm_add_epi32( sumHi , mmOffset ) , shift );
  sumLo = _mm_srai_epi32( _mm_add_epi32( sumLo , mmOffset ) , shift );
  return( _mm_packs_epi32( sumLo , sumHi ) );
}

inline __m128i simdInterpolateLuma2P8( Short const *src , Int srcStride , __m128i *mmCoeff , const __m128i & mmOffset , Int shift )
{
  __m128i sumHi = _mm_setzero_si128();
  __m128i sumLo = _mm_setzero_si128();
  for( Int n = 0 ; n < 2 ; n++ )
  {
    __m128i mmPix = _mm_loadu_si128( ( __m128i* )src );
    __m128i hi = _mm_mulhi_epi16( mmPix , mmCoeff[n] );
    __m128i lo = _mm_mullo_epi16( mmPix , mmCoeff[n] );
    sumHi = _mm_add_epi32( sumHi , _mm_unpackhi_epi16( lo , hi ) );
    sumLo = _mm_add_epi32( sumLo , _mm_unpacklo_epi16( lo , hi ) );
    src += srcStride;
  }
  sumHi = _mm_srai_epi32( _mm_add_epi32( sumHi , mmOffset ) , shift );
  sumLo = _mm_srai_epi32( _mm_add_epi32( sumLo , mmOffset ) , shift );
  return( _mm_packs_epi32( sumLo , sumHi ) );
}

inline __m128i simdInterpolateLuma2P4( Short const *src , Int srcStride , __m128i *mmCoeff , const __m128i & mmOffset , Int shift )
{
  __m128i sumHi = _mm_setzero_si128();
  __m128i sumLo = _mm_setzero_si128();
  for( Int n = 0 ; n < 2 ; n++ )
  {
    __m128i mmPix = _mm_loadl_epi64( ( __m128i* )src );
    __m128i hi = _mm_mulhi_epi16( mmPix , mmCoeff[n] );
    __m128i lo = _mm_mullo_epi16( mmPix , mmCoeff[n] );
    sumHi = _mm_add_epi32( sumHi , _mm_unpackhi_epi16( lo , hi ) );
    sumLo = _mm_add_epi32( sumLo , _mm_unpacklo_epi16( lo , hi ) );
    src += srcStride;
  }
  sumHi = _mm_srai_epi32( _mm_add_epi32( sumHi , mmOffset ) , shift );
  sumLo = _mm_srai_epi32( _mm_add_epi32( sumLo , mmOffset ) , shift );
  return( _mm_packs_epi32( sumLo , sumHi ) );
}

inline __m128i simdClip3( __m128i mmMin , __m128i mmMax , __m128i mmPix )
{
  __m128i mmMask = _mm_cmpgt_epi16( mmPix , mmMin );
  mmPix = _mm_or_si128( _mm_and_si128( mmMask , mmPix ) , _mm_andnot_si128( mmMask , mmMin ) );
  mmMask = _mm_cmplt_epi16( mmPix , mmMax );
  mmPix = _mm_or_si128( _mm_and_si128( mmMask , mmPix ) , _mm_andnot_si128( mmMask , mmMax ) );
  return( mmPix );
}
#endif

#if ITH_C3_INTERP_FILT==COARSE_SIMD_FXN
template < Bool isLast >
Void vertInterpolate8Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = srcStride;
	Int N = 8;

	if (!(iWidth & 0x07))
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 8; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 8)
			{
				__m128i mmFiltered = simdInterpolateLuma8(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storeu_si128((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}
	else if (!(iWidth & 0x03))
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 8; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateLuma4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
				//val = (val < iMinVal) ? iMinVal:(val>iMaxVal?iMaxVal:val);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}


template < Bool isLast >
Void horzInterpolate4Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = 1;
	Int N = 4;

	if (!(iWidth & 0x03))  /*if iWidth multiples of 4*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 4; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateChroma4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}

template < Bool isLast >
Void vertInterpolate2Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = srcStride;
	Int N = 2;

	if (!(iWidth & 0x07))  /*if iWidth multiples of 8*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[2];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 2; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 8)
			{
				__m128i mmFiltered = simdInterpolateLuma2P8(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storeu_si128((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}
	else if (!(iWidth & 0x03))  /*if iWidth multiples of 4*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 2; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateLuma2P4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}
		src += srcStride;
		dst += dstStride;
	}
}

template <Bool isLast>
Void horzInterpolate8Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = 1;
	Int N = 8;

	if (!(iWidth & 0x07))
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 8; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 8)
			{
				__m128i mmFiltered = simdInterpolateLuma8(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storeu_si128((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}
	else if (!(iWidth & 0x03))
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 8; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateLuma4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				// val = ClipA(val,compID);
				val = (val < iMinVal) ? iMinVal : (val>iMaxVal ? iMaxVal : val);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}

template < Bool isLast >
Void horzInterpolate2Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = 1;
	Int N = 2;

	if (!(iWidth & 0x07))  /*if iWidth multiples of 8*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[2];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 2; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 8)
			{
				__m128i mmFiltered = simdInterpolateLuma2P8(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storeu_si128((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}
	else if (!(iWidth & 0x03))  /*if iWidth multiples of 4*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 2; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateLuma2P4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}
		src += srcStride;
		dst += dstStride;
	}
}


template < Bool isLast >
Void vertInterpolate4Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = srcStride;
	Int N = 4;

	if (!(iWidth & 0x03))  /*if iWidth multiples of 4*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 4; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateChroma4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}

#elif ITH_C3_INTERP_FILT==SIMD_FXN
template < Bool isLast >
Void horzInterpolate8Tap
(
	Pel const *pSrc,
	Int iSrcStride,
	Pel *pDst,
	Int iDstStride,
	Short *pCoeff,
	Int iWidth,
	Int iHeight,
	Int iMinVal,
	Int iMaxVal,
	Int iOffset,
	Int iShift,
	ComponentID compID)
{
	Int iRow, iCol, iRemWd;
	Pel const *pSrcTmp;
	Pel const *inp_copy;
	Pel *dst_copy;

	/* all 128 bit registers are named with a suffix mxnb, where m is the */
	/* number of n bits packed in the register                            */

	__m128i offset_8x16b = _mm_set1_epi32(iOffset);
	__m128i	mmMin = _mm_set1_epi16(iMinVal);
	__m128i mmMax = _mm_set1_epi16(iMaxVal);
	__m128i src_temp1_16x8b, src_temp2_16x8b, src_temp3_16x8b, src_temp4_16x8b, src_temp5_16x8b, src_temp6_16x8b;
	__m128i src_temp7_16x8b, src_temp8_16x8b, src_temp9_16x8b, src_temp0_16x8b;
	__m128i src_temp11_16x8b, src_temp12_16x8b, src_temp13_16x8b, src_temp14_16x8b, src_temp15_16x8b, src_temp16_16x8b;
	__m128i res_temp1_8x16b, res_temp2_8x16b, res_temp3_8x16b, res_temp4_8x16b, res_temp5_8x16b, res_temp6_8x16b, res_temp7_8x16b, res_temp8_8x16b;
	__m128i res_temp9_8x16b, res_temp0_8x16b;
	__m128i res_temp11_8x16b, res_temp12_8x16b, res_temp13_8x16b, res_temp14_8x16b, res_temp15_8x16b, res_temp16_8x16b;
	__m128i coeff0_1_8x16b, coeff2_3_8x16b, coeff4_5_8x16b, coeff6_7_8x16b, mmMask;

	pSrcTmp = pSrc;
	iRemWd = iWidth;
	inp_copy = pSrcTmp;
	dst_copy = pDst;

	/* load 8 8-bit coefficients and convert 8-bit into 16-bit  */
	coeff0_1_8x16b = _mm_loadu_si128((__m128i*)pCoeff);

	coeff2_3_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x55);
	coeff4_5_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0xaa);
	coeff6_7_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0xff);
	coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);

	if (!(iHeight & 1))	/*even height*/
	{
		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < iHeight; iRow += 1)
			{
				int iCount = 0;
				for (iCol = iWidth; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));

					src_temp3_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp7_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp1_8x16b = _mm_madd_epi16(src_temp3_16x8b, coeff0_1_8x16b);
					res_temp7_8x16b = _mm_madd_epi16(src_temp7_16x8b, coeff0_1_8x16b);

					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));

					src_temp4_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp8_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp2_8x16b = _mm_madd_epi16(src_temp4_16x8b, coeff2_3_8x16b);
					res_temp8_8x16b = _mm_madd_epi16(src_temp8_16x8b, coeff2_3_8x16b);

					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 4));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 5));

					src_temp5_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp9_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp3_8x16b = _mm_madd_epi16(src_temp5_16x8b, coeff4_5_8x16b);
					res_temp9_8x16b = _mm_madd_epi16(src_temp9_16x8b, coeff4_5_8x16b);

					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 6));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 7));

					src_temp6_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp0_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp4_8x16b = _mm_madd_epi16(src_temp6_16x8b, coeff6_7_8x16b);
					res_temp0_8x16b = _mm_madd_epi16(src_temp0_16x8b, coeff6_7_8x16b);

					res_temp5_8x16b = _mm_add_epi32(res_temp1_8x16b, res_temp2_8x16b);
					res_temp6_8x16b = _mm_add_epi32(res_temp3_8x16b, res_temp4_8x16b);
					res_temp5_8x16b = _mm_add_epi32(res_temp5_8x16b, res_temp6_8x16b);

					res_temp6_8x16b = _mm_add_epi32(res_temp7_8x16b, res_temp8_8x16b);
					res_temp7_8x16b = _mm_add_epi32(res_temp9_8x16b, res_temp0_8x16b);
					res_temp8_8x16b = _mm_add_epi32(res_temp6_8x16b, res_temp7_8x16b);

					res_temp6_8x16b = _mm_add_epi32(res_temp5_8x16b, offset_8x16b);
					res_temp7_8x16b = _mm_add_epi32(res_temp8_8x16b, offset_8x16b);
					res_temp6_8x16b = _mm_srai_epi32(res_temp6_8x16b, iShift);
					res_temp7_8x16b = _mm_srai_epi32(res_temp7_8x16b, iShift);
					res_temp5_8x16b = _mm_packs_epi32(res_temp6_8x16b, res_temp7_8x16b);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res_temp5_8x16b, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res_temp5_8x16b, mmMax);
						res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res_temp5_8x16b);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += iSrcStride; /* pointer updates*/
				dst_copy += iDstStride; /* pointer updates*/
			}
		}

		iRemWd &= 0x7;
		inp_copy = pSrcTmp + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < iHeight; iRow += 2)
			{
				/*load 8 pixel values */
				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy));				/* row = 0 */
				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/* row = 1 */

				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 1));
				src_temp3_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp1_8x16b = _mm_madd_epi16(src_temp3_16x8b, coeff0_1_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 2));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 3));

				src_temp4_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp2_8x16b = _mm_madd_epi16(src_temp4_16x8b, coeff2_3_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 4));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 5));

				src_temp5_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp3_8x16b = _mm_madd_epi16(src_temp5_16x8b, coeff4_5_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 6));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 7));

				src_temp6_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp4_8x16b = _mm_madd_epi16(src_temp6_16x8b, coeff6_7_8x16b);

				res_temp5_8x16b = _mm_add_epi32(res_temp1_8x16b, res_temp2_8x16b);
				res_temp6_8x16b = _mm_add_epi32(res_temp3_8x16b, res_temp4_8x16b);
				res_temp5_8x16b = _mm_add_epi32(res_temp5_8x16b, res_temp6_8x16b);

				res_temp6_8x16b = _mm_add_epi32(res_temp5_8x16b, offset_8x16b);
				res_temp6_8x16b = _mm_srai_epi32(res_temp6_8x16b, iShift);
				res_temp5_8x16b = _mm_packs_epi32(res_temp6_8x16b, res_temp6_8x16b);

				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 1));

				src_temp13_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp11_8x16b = _mm_madd_epi16(src_temp13_16x8b, coeff0_1_8x16b);

				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 2));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 3));

				src_temp14_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp12_8x16b = _mm_madd_epi16(src_temp14_16x8b, coeff2_3_8x16b);

				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 4));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 5));

				src_temp15_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp13_8x16b = _mm_madd_epi16(src_temp15_16x8b, coeff4_5_8x16b);

				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 6));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 7));

				src_temp16_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp14_8x16b = _mm_madd_epi16(src_temp16_16x8b, coeff6_7_8x16b);

				res_temp15_8x16b = _mm_add_epi32(res_temp11_8x16b, res_temp12_8x16b);
				res_temp16_8x16b = _mm_add_epi32(res_temp13_8x16b, res_temp14_8x16b);
				res_temp15_8x16b = _mm_add_epi32(res_temp15_8x16b, res_temp16_8x16b);

				res_temp16_8x16b = _mm_add_epi32(res_temp15_8x16b, offset_8x16b);
				res_temp16_8x16b = _mm_srai_epi32(res_temp16_8x16b, iShift);
				res_temp15_8x16b = _mm_packs_epi32(res_temp16_8x16b, res_temp16_8x16b);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res_temp5_8x16b, mmMin);
					res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res_temp5_8x16b, mmMax);
					res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res_temp15_8x16b, mmMin);
					res_temp15_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp15_8x16b), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res_temp15_8x16b, mmMax);
					res_temp15_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp15_8x16b), _mm_andnot_si128(mmMask, mmMax));
				}

				/* to store the 1st 4 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res_temp5_8x16b);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), res_temp15_8x16b);
				inp_copy += (iSrcStride << 1);  /* Pointer update */
				dst_copy += (iDstStride << 1);  /* Pointer update */
			}
		}

		iRemWd &= 0x3;
		inp_copy = pSrcTmp + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			__m128i filt_coef;
			Pel sum, sum1;

			filt_coef = _mm_loadu_si128((__m128i*)pCoeff);   //w0 w1 w2 w3 w4 w5 w6 w7

			for (iRow = 0; iRow < iHeight; iRow += 2)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCol));
					src_temp5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCol));

					src_temp1_16x8b = _mm_madd_epi16(src_temp1_16x8b, filt_coef);
					src_temp5_16x8b = _mm_madd_epi16(src_temp5_16x8b, filt_coef);

					src_temp1_16x8b = _mm_hadd_epi32(src_temp1_16x8b, src_temp5_16x8b);
					src_temp1_16x8b = _mm_hadd_epi32(src_temp1_16x8b, src_temp1_16x8b);

					src_temp1_16x8b = _mm_add_epi32(src_temp1_16x8b, offset_8x16b);
					src_temp1_16x8b = _mm_srai_epi32(src_temp1_16x8b, iShift);
					src_temp1_16x8b = _mm_packs_epi32(src_temp1_16x8b, src_temp1_16x8b);

					sum = _mm_extract_epi16(src_temp1_16x8b, 0);
					sum1 = _mm_extract_epi16(src_temp1_16x8b, 1);

					if (isLast)
					{
						sum = (sum < iMinVal) ? iMinVal : (sum>iMaxVal ? iMaxVal : sum);
						sum1 = (sum1 < iMinVal) ? iMinVal : (sum1>iMaxVal ? iMaxVal : sum1);
					}
					dst_copy[iCol] = (sum);
					dst_copy[iCol + iDstStride] = (sum1);
				}
				inp_copy += (iSrcStride << 1);
				dst_copy += (iDstStride << 1);
			}
		}
	}
	else
	{
		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < iHeight; iRow += 1)
			{
				int iCount = 0;
				for (iCol = iWidth; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));

					src_temp3_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp7_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp1_8x16b = _mm_madd_epi16(src_temp3_16x8b, coeff0_1_8x16b);
					res_temp7_8x16b = _mm_madd_epi16(src_temp7_16x8b, coeff0_1_8x16b);
					/* row = 0 */
					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));

					src_temp4_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp8_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp2_8x16b = _mm_madd_epi16(src_temp4_16x8b, coeff2_3_8x16b);
					res_temp8_8x16b = _mm_madd_epi16(src_temp8_16x8b, coeff2_3_8x16b);

					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 4));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 5));

					src_temp5_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp9_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp3_8x16b = _mm_madd_epi16(src_temp5_16x8b, coeff4_5_8x16b);
					res_temp9_8x16b = _mm_madd_epi16(src_temp9_16x8b, coeff4_5_8x16b);

					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 6));
					src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 7));

					src_temp6_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
					src_temp0_16x8b = _mm_unpackhi_epi16(src_temp1_16x8b, src_temp2_16x8b);
					res_temp4_8x16b = _mm_madd_epi16(src_temp6_16x8b, coeff6_7_8x16b);
					res_temp0_8x16b = _mm_madd_epi16(src_temp0_16x8b, coeff6_7_8x16b);

					res_temp5_8x16b = _mm_add_epi32(res_temp1_8x16b, res_temp2_8x16b);
					res_temp6_8x16b = _mm_add_epi32(res_temp3_8x16b, res_temp4_8x16b);
					res_temp5_8x16b = _mm_add_epi32(res_temp5_8x16b, res_temp6_8x16b);

					res_temp6_8x16b = _mm_add_epi32(res_temp7_8x16b, res_temp8_8x16b);
					res_temp7_8x16b = _mm_add_epi32(res_temp9_8x16b, res_temp0_8x16b);
					res_temp8_8x16b = _mm_add_epi32(res_temp6_8x16b, res_temp7_8x16b);

					res_temp6_8x16b = _mm_add_epi32(res_temp5_8x16b, offset_8x16b);
					res_temp7_8x16b = _mm_add_epi32(res_temp8_8x16b, offset_8x16b);
					res_temp6_8x16b = _mm_srai_epi32(res_temp6_8x16b, iShift);
					res_temp7_8x16b = _mm_srai_epi32(res_temp7_8x16b, iShift);
					res_temp5_8x16b = _mm_packs_epi32(res_temp6_8x16b, res_temp7_8x16b);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res_temp5_8x16b, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res_temp5_8x16b, mmMax);
						res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res_temp5_8x16b);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += iSrcStride; /* pointer updates*/
				dst_copy += iDstStride; /* pointer updates*/
			}
		}

		iRemWd &= 0x7;
		inp_copy = pSrcTmp + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < (iHeight - 1); iRow += 2)
			{

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy));
				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));

				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 1));

				src_temp3_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp1_8x16b = _mm_madd_epi16(src_temp3_16x8b, coeff0_1_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 2));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 3));

				src_temp4_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp2_8x16b = _mm_madd_epi16(src_temp4_16x8b, coeff2_3_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 4));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 5));

				src_temp5_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp3_8x16b = _mm_madd_epi16(src_temp5_16x8b, coeff4_5_8x16b);

				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 6));
				src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 7));

				src_temp6_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
				res_temp4_8x16b = _mm_madd_epi16(src_temp6_16x8b, coeff6_7_8x16b);

				res_temp5_8x16b = _mm_add_epi32(res_temp1_8x16b, res_temp2_8x16b);
				res_temp6_8x16b = _mm_add_epi32(res_temp3_8x16b, res_temp4_8x16b);
				res_temp5_8x16b = _mm_add_epi32(res_temp5_8x16b, res_temp6_8x16b);

				res_temp6_8x16b = _mm_add_epi32(res_temp5_8x16b, offset_8x16b);
				res_temp6_8x16b = _mm_srai_epi32(res_temp6_8x16b, iShift);
				res_temp5_8x16b = _mm_packs_epi32(res_temp6_8x16b, res_temp6_8x16b);

				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 1));

				src_temp13_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp11_8x16b = _mm_madd_epi16(src_temp13_16x8b, coeff0_1_8x16b);
				/* row =1 */
				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 2));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 3));

				src_temp14_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp12_8x16b = _mm_madd_epi16(src_temp14_16x8b, coeff2_3_8x16b);

				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 4));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 5));

				src_temp15_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp13_8x16b = _mm_madd_epi16(src_temp15_16x8b, coeff4_5_8x16b);

				src_temp11_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 6));
				src_temp12_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + 7));

				src_temp16_16x8b = _mm_unpacklo_epi16(src_temp11_16x8b, src_temp12_16x8b);
				res_temp14_8x16b = _mm_madd_epi16(src_temp16_16x8b, coeff6_7_8x16b);

				res_temp15_8x16b = _mm_add_epi32(res_temp11_8x16b, res_temp12_8x16b);
				res_temp16_8x16b = _mm_add_epi32(res_temp13_8x16b, res_temp14_8x16b);
				res_temp15_8x16b = _mm_add_epi32(res_temp15_8x16b, res_temp16_8x16b);

				res_temp16_8x16b = _mm_add_epi32(res_temp15_8x16b, offset_8x16b);
				res_temp16_8x16b = _mm_srai_epi32(res_temp16_8x16b, iShift);
				res_temp15_8x16b = _mm_packs_epi32(res_temp16_8x16b, res_temp16_8x16b);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res_temp5_8x16b, mmMin);
					res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res_temp5_8x16b, mmMax);
					res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res_temp15_8x16b, mmMin);
					res_temp15_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp15_8x16b), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res_temp15_8x16b, mmMax);
					res_temp15_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp15_8x16b), _mm_andnot_si128(mmMask, mmMax));
				}

				/* to store the 1st 4 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res_temp5_8x16b);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), res_temp15_8x16b);
				inp_copy += (iSrcStride << 1);  /* Pointer update */
				dst_copy += (iDstStride << 1);  /* Pointer update */
			}

			/*extra one height to be done*/
			src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy));

			src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 1));

			src_temp3_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
			res_temp1_8x16b = _mm_madd_epi16(src_temp3_16x8b, coeff0_1_8x16b);

			src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 2));
			src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 3));

			src_temp4_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
			res_temp2_8x16b = _mm_madd_epi16(src_temp4_16x8b, coeff2_3_8x16b);

			src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 4));
			src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 5));

			src_temp5_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
			res_temp3_8x16b = _mm_madd_epi16(src_temp5_16x8b, coeff4_5_8x16b);

			src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 6));
			src_temp2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + 7));

			src_temp6_16x8b = _mm_unpacklo_epi16(src_temp1_16x8b, src_temp2_16x8b);
			res_temp4_8x16b = _mm_madd_epi16(src_temp6_16x8b, coeff6_7_8x16b);

			res_temp5_8x16b = _mm_add_epi32(res_temp1_8x16b, res_temp2_8x16b);
			res_temp6_8x16b = _mm_add_epi32(res_temp3_8x16b, res_temp4_8x16b);
			res_temp5_8x16b = _mm_add_epi32(res_temp5_8x16b, res_temp6_8x16b);

			res_temp6_8x16b = _mm_add_epi32(res_temp5_8x16b, offset_8x16b);
			res_temp6_8x16b = _mm_srai_epi32(res_temp6_8x16b, iShift);
			res_temp5_8x16b = _mm_packs_epi32(res_temp6_8x16b, res_temp6_8x16b);

			if (isLast)
			{
				mmMask = _mm_cmpgt_epi16(res_temp5_8x16b, mmMin);
				res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMin));
				mmMask = _mm_cmplt_epi16(res_temp5_8x16b, mmMax);
				res_temp5_8x16b = _mm_or_si128(_mm_and_si128(mmMask, res_temp5_8x16b), _mm_andnot_si128(mmMask, mmMax));
			}

			_mm_storel_epi64((__m128i *)(dst_copy), res_temp5_8x16b);

		}

		iRemWd &= 0x3;
		inp_copy = pSrcTmp + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			__m128i filt_coef;
			Pel sum, sum1;

			filt_coef = _mm_loadu_si128((__m128i*)pCoeff);   //w0 w1 w2 w3 w4 w5 w6 w7

			for (iRow = 0; iRow < (iHeight - 1); iRow += 2)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCol));
					src_temp5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCol));

					src_temp1_16x8b = _mm_madd_epi16(src_temp1_16x8b, filt_coef);
					src_temp5_16x8b = _mm_madd_epi16(src_temp5_16x8b, filt_coef);

					src_temp1_16x8b = _mm_hadd_epi32(src_temp1_16x8b, src_temp5_16x8b);
					src_temp1_16x8b = _mm_hadd_epi32(src_temp1_16x8b, src_temp1_16x8b);

					src_temp1_16x8b = _mm_add_epi32(src_temp1_16x8b, offset_8x16b);
					src_temp1_16x8b = _mm_srai_epi32(src_temp1_16x8b, iShift);
					src_temp1_16x8b = _mm_packs_epi32(src_temp1_16x8b, src_temp1_16x8b);

					sum = _mm_extract_epi16(src_temp1_16x8b, 0);
					sum1 = _mm_extract_epi16(src_temp1_16x8b, 1);

					if (isLast)
					{
						sum = (sum < iMinVal) ? iMinVal : (sum>iMaxVal ? iMaxVal : sum);
						sum1 = (sum1 < iMinVal) ? iMinVal : (sum1>iMaxVal ? iMaxVal : sum1);
					}
					dst_copy[iCol] = (sum);
					dst_copy[iCol + iDstStride] = (sum1);
				}
				inp_copy += (iSrcStride << 1);
				dst_copy += (iDstStride << 1);
			}

			for (iCol = 0; iCol < iRemWd; iCol++)
			{
				src_temp1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCol));

				src_temp1_16x8b = _mm_madd_epi16(src_temp1_16x8b, filt_coef);

				src_temp1_16x8b = _mm_hadd_epi32(src_temp1_16x8b, src_temp1_16x8b);
				src_temp2_16x8b = _mm_srli_si128(src_temp1_16x8b, 4);
				src_temp1_16x8b = _mm_add_epi32(src_temp1_16x8b, src_temp2_16x8b);

				src_temp1_16x8b = _mm_add_epi32(src_temp1_16x8b, offset_8x16b);
				src_temp1_16x8b = _mm_srai_epi32(src_temp1_16x8b, iShift);
				src_temp1_16x8b = _mm_packs_epi32(src_temp1_16x8b, src_temp1_16x8b);

				sum = (Pel)_mm_extract_epi16(src_temp1_16x8b, 0);

				if (isLast)
				{
					sum = (sum < iMinVal) ? iMinVal : (sum>iMaxVal ? iMaxVal : sum);
				}
				dst_copy[iCol] = (sum);
			}

		}
	}
}


template < Bool isLast >
Void horzInterpolate2Tap
(
	Pel const *pSrc,
	Int iSrcStride,
	Pel *pDst,
	Int iDstStride,
	Short *pCoeff,
	Int iWidth,
	Int iHeight,
	Int iMinVal,
	Int iMaxVal,
	Int iOffset,
	Int iShift,
	ComponentID compID)
{

	assert(iHeight >= 4);	/*height has to be >= 4*/

	Int iRow, iCol, iRemWd;
	Int iSrcStride2, iSrcStride3;
	Pel const *inp_copy;
	Pel *dst_copy;

	__m128i offset_4x32b = _mm_set1_epi32(iOffset);
	__m128i	mmMin = _mm_set1_epi16(iMinVal);
	__m128i mmMax = _mm_set1_epi16(iMaxVal);

	__m128i row1, row11, row2, row22, row3, row33, row4, row44;
	__m128i res0, res1, res2, res3;
	__m128i coeff0_1_8x16b, shuffle, mmMask;

	iRemWd = iWidth;
	inp_copy = pSrc;
	dst_copy = pDst;

	iSrcStride2 = (iSrcStride << 1);
	iSrcStride3 = (iSrcStride * 3);

	/* load 8 8-bit coefficients and convert 8-bit into 16-bit  */
	coeff0_1_8x16b = _mm_loadl_epi64((__m128i*)pCoeff);			/*w0 w1 x x x x x x*/
	coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);		/*w0 w1 w0 w1 w0 w1 w0 w1*/

	shuffle = _mm_loadu_si128((__m128i*)shuffle_2Tap);

	if (!(iHeight & 3))	/*height multiples of 4*/
	{
		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < iHeight; iRow += 4)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));				/*a1 a2 a3 a4 a5 a6 a7 a8*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));		/*b0 b1 b2 b3 b4 b5 b6 b7*/
					row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));	/*b1 b2 b3 b4 b5 b6 b7 b8*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row33 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 1));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
					row44 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 1));

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/
					row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
					row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);
					row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row11 = _mm_add_epi32(row11, offset_4x32b);
					row2 = _mm_add_epi32(row2, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row3 = _mm_add_epi32(row3, offset_4x32b);
					row33 = _mm_add_epi32(row33, offset_4x32b);
					row4 = _mm_add_epi32(row4, offset_4x32b);
					row44 = _mm_add_epi32(row44, offset_4x32b);

					row1 = _mm_srai_epi32(row1, iShift);
					row11 = _mm_srai_epi32(row11, iShift);
					row2 = _mm_srai_epi32(row2, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row3 = _mm_srai_epi32(row3, iShift);
					row33 = _mm_srai_epi32(row33, iShift);
					row4 = _mm_srai_epi32(row4, iShift);
					row44 = _mm_srai_epi32(row44, iShift);

					row1 = _mm_packs_epi32(row1, row2);
					row11 = _mm_packs_epi32(row11, row22);
					row3 = _mm_packs_epi32(row3, row4);
					row33 = _mm_packs_epi32(row33, row44);

					res0 = _mm_unpacklo_epi16(row1, row11);
					res1 = _mm_unpackhi_epi16(row1, row11);
					res2 = _mm_unpacklo_epi16(row3, row33);
					res3 = _mm_unpackhi_epi16(row3, row33);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res2, mmMin);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res3, mmMin);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res3, mmMax);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 2 + iCount), res2);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 3 + iCount), res3);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}
		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < iHeight; iRow += 4)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 a4 a5 a6 a7*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/*a1 a2 a3 a4 a5 a6 a7 a8*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
				row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));

				row1 = _mm_shuffle_epi8(row1, shuffle);	/*a0 a1 a1 a2 a2 a3 a3 a4 */
				row2 = _mm_shuffle_epi8(row2, shuffle);
				row3 = _mm_shuffle_epi8(row3, shuffle);
				row4 = _mm_shuffle_epi8(row4, shuffle);

				row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
				row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
				row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);

				row1 = _mm_add_epi32(row1, offset_4x32b);
				row2 = _mm_add_epi32(row2, offset_4x32b);
				row3 = _mm_add_epi32(row3, offset_4x32b);
				row4 = _mm_add_epi32(row4, offset_4x32b);

				row1 = _mm_srai_epi32(row1, iShift);
				row2 = _mm_srai_epi32(row2, iShift);
				row3 = _mm_srai_epi32(row3, iShift);
				row4 = _mm_srai_epi32(row4, iShift);

				res0 = _mm_packs_epi32(row1, row2);
				res1 = _mm_packs_epi32(row3, row4);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
				}

				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 2), res1);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 3), _mm_unpackhi_epi64(res1, res1));

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}
		}

		iRemWd &= 0x3;
		inp_copy = pSrc + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			Int iSum, iSum1;
			for (iRow = 0; iRow < (iHeight); iRow += 4)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCol));				/*a0 a1 x x x x x x*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCol));	/*b0 b1 x x x x x x*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCol));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCol));

					row1 = _mm_unpacklo_epi32(row1, row2);	/*a0 a1 b0 b1*/
					row3 = _mm_unpacklo_epi32(row3, row4); /*c0 c1 d0 d1*/
					row1 = _mm_unpacklo_epi64(row1, row3);	/*a0 a1 b0 b1 c0 c1 d0 d1*/

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 b0+b1 c0+c1 d0+d1*/

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row1 = _mm_srai_epi32(row1, iShift);
					res0 = _mm_packs_epi32(row1, row1);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));
					}

					/*extract 32 bit integer form register and store it in dst_copy*/
					iSum = _mm_extract_epi32(res0, 0);
					iSum1 = _mm_extract_epi32(res0, 1);

					dst_copy[iCol] = (Pel)iSum;
					dst_copy[iCol + iDstStride] = (Pel)(iSum >> 16);
					dst_copy[iCol + (iDstStride << 1)] = (Pel)iSum1;
					dst_copy[iCol + (iDstStride * 3)] = (Pel)(iSum1 >> 16);
				}
				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

		}
	}
	else
	{
		Int iRemHt = (iHeight % 4);

		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));				/*a1 a2 a3 a4 a5 a6 a7 a8*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));		/*b0 b1 b2 b3 b4 b5 b6 b7*/
					row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));	/*b1 b2 b3 b4 b5 b6 b7 b8*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row33 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 1));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
					row44 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 1));

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/
					row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
					row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);
					row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row11 = _mm_add_epi32(row11, offset_4x32b);
					row2 = _mm_add_epi32(row2, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row3 = _mm_add_epi32(row3, offset_4x32b);
					row33 = _mm_add_epi32(row33, offset_4x32b);
					row4 = _mm_add_epi32(row4, offset_4x32b);
					row44 = _mm_add_epi32(row44, offset_4x32b);

					row1 = _mm_srai_epi32(row1, iShift);
					row11 = _mm_srai_epi32(row11, iShift);
					row2 = _mm_srai_epi32(row2, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row3 = _mm_srai_epi32(row3, iShift);
					row33 = _mm_srai_epi32(row33, iShift);
					row4 = _mm_srai_epi32(row4, iShift);
					row44 = _mm_srai_epi32(row44, iShift);

					row1 = _mm_packs_epi32(row1, row2);
					row11 = _mm_packs_epi32(row11, row22);
					row3 = _mm_packs_epi32(row3, row4);
					row33 = _mm_packs_epi32(row33, row44);

					res0 = _mm_unpacklo_epi16(row1, row11);
					res1 = _mm_unpackhi_epi16(row1, row11);
					res2 = _mm_unpacklo_epi16(row3, row33);
					res3 = _mm_unpackhi_epi16(row3, row33);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res2, mmMin);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res3, mmMin);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res3, mmMax);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 2 + iCount), res2);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 3 + iCount), res3);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			/*extra height to be done --- one row at a time*/
			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));				/*a1 a2 a3 a4 a5 a6 a7 a8*/

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row11 = _mm_add_epi32(row11, offset_4x32b);

					row1 = _mm_srai_epi32(row1, iShift);
					row11 = _mm_srai_epi32(row11, iShift);

					row1 = _mm_packs_epi32(row1, row11);	/*a0 a2 a4 a6 a1 a3 a5 a7*/

					res0 = _mm_unpackhi_epi64(row1, row1);	/*a1 a3 a5 a7*/
					res1 = _mm_unpacklo_epi16(row1, res0);	/*a0 a1 a2 a3 a4 a5 a6 a7*/

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res1);

					iCount += 8;
				}

				inp_copy += (iSrcStride << 1);
				dst_copy += (iDstStride << 1);
			}

		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 a4 a5 a6 a7*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/*a1 a2 a3 a4 a5 a6 a7 a8*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
				row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));

				row1 = _mm_shuffle_epi8(row1, shuffle);	/*a0 a1 a1 a2 a2 a3 a3 a4 */
				row2 = _mm_shuffle_epi8(row2, shuffle);
				row3 = _mm_shuffle_epi8(row3, shuffle);
				row4 = _mm_shuffle_epi8(row4, shuffle);

				row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
				row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
				row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);

				row1 = _mm_add_epi32(row1, offset_4x32b);
				row2 = _mm_add_epi32(row2, offset_4x32b);
				row3 = _mm_add_epi32(row3, offset_4x32b);
				row4 = _mm_add_epi32(row4, offset_4x32b);

				row1 = _mm_srai_epi32(row1, iShift);
				row2 = _mm_srai_epi32(row2, iShift);
				row3 = _mm_srai_epi32(row3, iShift);
				row4 = _mm_srai_epi32(row4, iShift);

				res0 = _mm_packs_epi32(row1, row2);
				res1 = _mm_packs_epi32(row3, row4);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
				}

				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 2), res1);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 3), _mm_unpackhi_epi64(res1, res1));

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 a4 a5 a6 a7*/

				res0 = _mm_shuffle_epi8(row1, shuffle);	/*a0 a1 a1 a2 a2 a3 a3 a4 */
				res0 = _mm_madd_epi16(res0, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				res0 = _mm_add_epi32(res0, offset_4x32b);
				res0 = _mm_srai_epi32(res0, iShift);
				res0 = _mm_packs_epi32(res0, res0);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

				}

				_mm_storel_epi64((__m128i *)(dst_copy), res0);

				inp_copy += (iSrcStride);
				dst_copy += (iDstStride);
			}
		}

		iRemWd &= 0x3;
		inp_copy = pSrc + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			Int iSum, iSum1;
			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCol));				/*a0 a1 x x x x x x*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCol));	/*b0 b1 x x x x x x*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCol));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCol));

					row1 = _mm_unpacklo_epi32(row1, row2);	/*a0 a1 b0 b1*/
					row3 = _mm_unpacklo_epi32(row3, row4); /*c0 c1 d0 d1*/
					row1 = _mm_unpacklo_epi64(row1, row3);	/*a0 a1 b0 b1 c0 c1 d0 d1*/

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 b0+b1 c0+c1 d0+d1*/

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row1 = _mm_srai_epi32(row1, iShift);
					res0 = _mm_packs_epi32(row1, row1);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					}

					/*extract 32 bit integer form register and store it in dst_copy*/
					iSum = _mm_extract_epi32(res0, 0);
					iSum1 = _mm_extract_epi32(res0, 1);

					dst_copy[iCol] = (Pel)iSum;
					dst_copy[iCol + iDstStride] = (Pel)(iSum >> 16);
					dst_copy[iCol + (iDstStride << 1)] = (Pel)iSum1;
					dst_copy[iCol + (iDstStride * 3)] = (Pel)(iSum1 >> 16);
				}
				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					Int sum;

					sum = inp_copy[iCol + 0] * pCoeff[0];
					sum += inp_copy[iCol + 1] * pCoeff[1];

					Pel val = (sum + iOffset) >> iShift;
					if (isLast)
					{
						val = (val < iMinVal) ? iMinVal : (val>iMaxVal ? iMaxVal : val);
					}
					dst_copy[iCol] = val;
				}
				inp_copy += iSrcStride;
				dst_copy += iDstStride;
			}
		}
	}
}


template < Bool isLast >
Void vertInterpolate2Tap
(
	Pel const *pSrc,
	Int iSrcStride,
	Pel *pDst,
	Int iDstStride,
	Short *pCoeff,
	Int iWidth,
	Int iHeight,
	Int iMinVal,
	Int iMaxVal,
	Int iOffset,
	Int iShift,
	ComponentID compID)
{
	assert(iHeight >= 4);	/*height has to be >= 4*/

	Int iRow, iCol, iRemWd;
	Int iSrcStride2, iSrcStride3, iSrcStride4;
	Pel const *inp_copy;
	Pel *dst_copy;

	__m128i offset_4x32b = _mm_set1_epi32(iOffset);
	__m128i	mmMin = _mm_set1_epi16(iMinVal);
	__m128i mmMax = _mm_set1_epi16(iMaxVal);

	__m128i row1, row11, row2, row22, row3, row33, row4, row44, row5;
	__m128i res0, res1, res2, res3, mmMask;
	__m128i coeff0_1_8x16b, mmMask0, mmMask1, mmMask2, mmMask3;

	iRemWd = iWidth;
	inp_copy = pSrc;
	dst_copy = pDst;

	iSrcStride2 = (iSrcStride << 1);
	iSrcStride3 = (iSrcStride * 3);
	iSrcStride4 = (iSrcStride << 2);

	coeff0_1_8x16b = _mm_loadl_epi64((__m128i*)pCoeff);			/*w0 w1 x x x x x x*/
	coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);		/*w0 w1 w0 w1 w0 w1 w0 w1*/

	if (!(iHeight & 3))	/*height multiples of 4*/
	{
		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < iHeight - 1; iRow += 4)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));	/*b0 b1 b2 b3 b4 b5 b6 b7*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
					row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4 + iCount));

					row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
					row1 = _mm_unpackhi_epi16(row1, row2);		/*a4 b4 a5 b5 a6 b6 a7 b7*/
					row22 = _mm_unpacklo_epi16(row2, row3);
					row2 = _mm_unpackhi_epi16(row2, row3);
					row33 = _mm_unpacklo_epi16(row3, row4);
					row3 = _mm_unpackhi_epi16(row3, row4);
					row44 = _mm_unpacklo_epi16(row4, row5);
					row4 = _mm_unpackhi_epi16(row4, row5);

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
					row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
					row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);
					row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);

					row11 = _mm_add_epi32(row11, offset_4x32b);
					row1 = _mm_add_epi32(row1, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row2 = _mm_add_epi32(row2, offset_4x32b);
					row33 = _mm_add_epi32(row33, offset_4x32b);
					row3 = _mm_add_epi32(row3, offset_4x32b);
					row44 = _mm_add_epi32(row44, offset_4x32b);
					row4 = _mm_add_epi32(row4, offset_4x32b);

					row11 = _mm_srai_epi32(row11, iShift);
					row1 = _mm_srai_epi32(row1, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row2 = _mm_srai_epi32(row2, iShift);
					row33 = _mm_srai_epi32(row33, iShift);
					row3 = _mm_srai_epi32(row3, iShift);
					row44 = _mm_srai_epi32(row44, iShift);
					row4 = _mm_srai_epi32(row4, iShift);

					res0 = _mm_packs_epi32(row11, row1);
					res1 = _mm_packs_epi32(row22, row2);
					res2 = _mm_packs_epi32(row33, row3);
					res3 = _mm_packs_epi32(row44, row4);

					if (isLast)
					{
						mmMask0 = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask0, res0), _mm_andnot_si128(mmMask0, mmMin));
						mmMask0 = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask0, res0), _mm_andnot_si128(mmMask0, mmMax));

						mmMask1 = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask1, res1), _mm_andnot_si128(mmMask1, mmMin));
						mmMask1 = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask1, res1), _mm_andnot_si128(mmMask1, mmMax));

						mmMask2 = _mm_cmpgt_epi16(res2, mmMin);
						res2 = _mm_or_si128(_mm_and_si128(mmMask2, res2), _mm_andnot_si128(mmMask2, mmMin));
						mmMask2 = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask2, res2), _mm_andnot_si128(mmMask2, mmMax));

						mmMask3 = _mm_cmpgt_epi16(res3, mmMin);
						res3 = _mm_or_si128(_mm_and_si128(mmMask3, res3), _mm_andnot_si128(mmMask3, mmMin));
						mmMask3 = _mm_cmplt_epi16(res3, mmMax);
						res3 = _mm_or_si128(_mm_and_si128(mmMask3, res3), _mm_andnot_si128(mmMask3, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 2 + iCount), res2);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 3 + iCount), res3);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}
		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < iHeight - 1; iRow += 4)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 x x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/*b0 b1 b2 b3 x x x x*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
				row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));
				row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4));

				row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
				row22 = _mm_unpacklo_epi16(row2, row3);
				row33 = _mm_unpacklo_epi16(row3, row4);
				row44 = _mm_unpacklo_epi16(row4, row5);

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
				row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
				row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row22 = _mm_add_epi32(row22, offset_4x32b);
				row33 = _mm_add_epi32(row33, offset_4x32b);
				row44 = _mm_add_epi32(row44, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row22 = _mm_srai_epi32(row22, iShift);
				row33 = _mm_srai_epi32(row33, iShift);
				row44 = _mm_srai_epi32(row44, iShift);

				res0 = _mm_packs_epi32(row11, row22);
				res1 = _mm_packs_epi32(row33, row44);

				if (isLast)
				{
					mmMask0 = _mm_cmpgt_epi16(res0, mmMin);
					res0 = _mm_or_si128(_mm_and_si128(mmMask0, res0), _mm_andnot_si128(mmMask0, mmMin));
					mmMask0 = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask0, res0), _mm_andnot_si128(mmMask0, mmMax));

					mmMask1 = _mm_cmpgt_epi16(res1, mmMin);
					res1 = _mm_or_si128(_mm_and_si128(mmMask1, res1), _mm_andnot_si128(mmMask1, mmMin));
					mmMask1 = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask1, res1), _mm_andnot_si128(mmMask1, mmMax));
				}

				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 2), res1);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 3), _mm_unpackhi_epi64(res1, res1));

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}
		}

		iRemWd &= 0x3;
		inp_copy = pSrc + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			Int iSum, iSum1;

			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{

					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy));					/*a0 x x x x x x x*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));		/*b0 x x x x x x x*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));
					row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4));

					row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
					row22 = _mm_unpacklo_epi16(row2, row3);
					row33 = _mm_unpacklo_epi16(row3, row4);
					row44 = _mm_unpacklo_epi16(row4, row5);

					row11 = _mm_unpacklo_epi32(row11, row22);
					row33 = _mm_unpacklo_epi32(row33, row44);

					row11 = _mm_unpacklo_epi64(row11, row33);

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
					row11 = _mm_add_epi32(row11, offset_4x32b);
					row11 = _mm_srai_epi32(row11, iShift);
					row11 = _mm_packs_epi32(row11, row11);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(row11, mmMin);
						row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(row11, mmMax);
						row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));
					}

					iSum = (Int)_mm_cvtsi128_si32(row11);
					iSum1 = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 8));

					dst_copy[iCol] = (iSum & 0xffff);
					dst_copy[iCol + iDstStride] = (iSum >> 16) & 0xffff;
					dst_copy[iCol + iDstStride * 2] = (iSum1 & 0xffff);
					dst_copy[iCol + iDstStride * 3] = (iSum1 >> 16) & 0xffff;

					inp_copy += (iSrcStride << 2);
					dst_copy += (iDstStride << 2);
				}
			}
		}

	}
	else
	{
		Int iRemHt = (iHeight % 4);

		if (iRemWd > 7)
		{
			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));	/*b0 b1 b2 b3 b4 b5 b6 b7*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
					row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4 + iCount));

					row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
					row1 = _mm_unpackhi_epi16(row1, row2);		/*a4 b4 a5 b5 a6 b6 a7 b7*/
					row22 = _mm_unpacklo_epi16(row2, row3);
					row2 = _mm_unpackhi_epi16(row2, row3);
					row33 = _mm_unpacklo_epi16(row3, row4);
					row3 = _mm_unpackhi_epi16(row3, row4);
					row44 = _mm_unpacklo_epi16(row4, row5);
					row4 = _mm_unpackhi_epi16(row4, row5);

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row2 = _mm_madd_epi16(row2, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
					row3 = _mm_madd_epi16(row3, coeff0_1_8x16b);
					row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);
					row4 = _mm_madd_epi16(row4, coeff0_1_8x16b);

					row11 = _mm_add_epi32(row11, offset_4x32b);
					row1 = _mm_add_epi32(row1, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row2 = _mm_add_epi32(row2, offset_4x32b);
					row33 = _mm_add_epi32(row33, offset_4x32b);
					row3 = _mm_add_epi32(row3, offset_4x32b);
					row44 = _mm_add_epi32(row44, offset_4x32b);
					row4 = _mm_add_epi32(row4, offset_4x32b);

					row11 = _mm_srai_epi32(row11, iShift);
					row1 = _mm_srai_epi32(row1, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row2 = _mm_srai_epi32(row2, iShift);
					row33 = _mm_srai_epi32(row33, iShift);
					row3 = _mm_srai_epi32(row3, iShift);
					row44 = _mm_srai_epi32(row44, iShift);
					row4 = _mm_srai_epi32(row4, iShift);

					res0 = _mm_packs_epi32(row11, row1);
					res1 = _mm_packs_epi32(row22, row2);
					res2 = _mm_packs_epi32(row33, row3);
					res3 = _mm_packs_epi32(row44, row4);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res2, mmMin);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res3, mmMin);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res3, mmMax);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 2 + iCount), res2);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride * 3 + iCount), res3);

					iCount += 8; /* To pointer updates*/
				}

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			/*extra height to be done --- one row at a time*/
			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				int iCount = 0;
				for (iCol = iRemWd; iCol > 7; iCol -= 8)
				{
					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));					/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));	/*b0 b1 b2 b3 b4 b5 b6 b7*/

					row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
					row1 = _mm_unpackhi_epi16(row1, row2);		/*a4 b4 a5 b5 a6 b6 a7 b7*/

					row1 = _mm_madd_epi16(row1, coeff0_1_8x16b);			/*a0+a1 a2+a3 a4+a5 a6+a7*/
					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a1+a2 a3+a4 a5+a6 a7+a8*/

					row1 = _mm_add_epi32(row1, offset_4x32b);
					row11 = _mm_add_epi32(row11, offset_4x32b);

					row1 = _mm_srai_epi32(row1, iShift);
					row11 = _mm_srai_epi32(row11, iShift);

					res1 = _mm_packs_epi32(row11, row1);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res1, mmMin);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res1);

					iCount += 8;
				}

				inp_copy += (iSrcStride << 1);
				dst_copy += (iDstStride << 1);
			}

		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 x x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/*b0 b1 b2 b3 x x x x*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
				row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));
				row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4));

				row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
				row22 = _mm_unpacklo_epi16(row2, row3);
				row33 = _mm_unpacklo_epi16(row3, row4);
				row44 = _mm_unpacklo_epi16(row4, row5);

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
				row33 = _mm_madd_epi16(row33, coeff0_1_8x16b);
				row44 = _mm_madd_epi16(row44, coeff0_1_8x16b);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row22 = _mm_add_epi32(row22, offset_4x32b);
				row33 = _mm_add_epi32(row33, offset_4x32b);
				row44 = _mm_add_epi32(row44, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row22 = _mm_srai_epi32(row22, iShift);
				row33 = _mm_srai_epi32(row33, iShift);
				row44 = _mm_srai_epi32(row44, iShift);

				res0 = _mm_packs_epi32(row11, row22);
				res1 = _mm_packs_epi32(row33, row44);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
				}

				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 2), res1);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride * 3), _mm_unpackhi_epi64(res1, res1));

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				/*load 8 pixel values from row 0*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));				/*a0 a1 a2 a3 x x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));	/*b0 b1 b2 b3 x x x x*/

				row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
				row11 = _mm_add_epi32(row11, offset_4x32b);
				row11 = _mm_srai_epi32(row11, iShift);
				row11 = _mm_packs_epi32(row11, row11);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(row11, mmMin);
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(row11, mmMax);
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));

				}

				_mm_storel_epi64((__m128i *)(dst_copy), row11);

				inp_copy += (iSrcStride);
				dst_copy += (iDstStride);
			}
		}

		iRemWd &= 0x3;
		inp_copy = pSrc + ((iWidth / 4) * 4);
		dst_copy = pDst + ((iWidth / 4) * 4);

		if (iRemWd)
		{
			Int iSum, iSum1;

			for (iRow = 0; iRow < (iHeight - 1); iRow += 4)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{

					/*load 8 pixel values from row 0*/
					row1 = _mm_loadu_si128((__m128i*)(inp_copy));					/*a0 x x x x x x x*/
					row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));		/*b0 x x x x x x x*/
					row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
					row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));
					row5 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride4));

					row11 = _mm_unpacklo_epi16(row1, row2);		/*a0 b0 a1 b1 a2 b2 a3 b3*/
					row22 = _mm_unpacklo_epi16(row2, row3);
					row33 = _mm_unpacklo_epi16(row3, row4);
					row44 = _mm_unpacklo_epi16(row4, row5);

					row11 = _mm_unpacklo_epi32(row11, row22);
					row33 = _mm_unpacklo_epi32(row33, row44);

					row11 = _mm_unpacklo_epi64(row11, row33);

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);			/*a0+a1 a1+a2 a2+a3 a3+a4*/
					row11 = _mm_add_epi32(row11, offset_4x32b);
					row11 = _mm_srai_epi32(row11, iShift);
					row11 = _mm_packs_epi32(row11, row11);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(row11, mmMin);
						row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(row11, mmMax);
						row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));
					}

					iSum = (Int)_mm_cvtsi128_si32(row11);
					iSum1 = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 8));

					dst_copy[iCol] = (iSum & 0xffff);
					dst_copy[iCol + iDstStride] = (iSum >> 16) & 0xffff;
					dst_copy[iCol + iDstStride * 2] = (iSum1 & 0xffff);
					dst_copy[iCol + iDstStride * 3] = (iSum1 >> 16) & 0xffff;

					inp_copy += (iSrcStride << 2);
					dst_copy += (iDstStride << 2);
				}
			}

			for (iRow = 0; iRow < iRemHt; iRow++)
			{
				for (iCol = 0; iCol < iRemWd; iCol++)
				{
					Int sum;

					sum = inp_copy[iCol + 0 * iSrcStride] * pCoeff[0];
					sum += inp_copy[iCol + 1 * iSrcStride] * pCoeff[1];

					Pel val = (sum + iOffset) >> iShift;
					if (isLast)
					{
#if JVET_D0033_ADAPTIVE_CLIPPING
						// val = ClipA(val,compID);
						val = (val < iMinVal) ? iMinVal : (val>iMaxVal ? iMaxVal : val);
#else
						val = (val < 0) ? 0 : val;
						val = (val > maxVal) ? maxVal : val;
#endif
					}
					dst_copy[iCol] = val;
				}

				inp_copy += iSrcStride;
				dst_copy += iDstStride;
			}

		}
	}
}


template <Bool isLast>
Void vertInterpolate8Tap
(
	Pel const *pSrc,
	Int iSrcStride,
	Pel *pDst,
	Int iDstStride,
	Short *pCoeff,
	Int iWidth,
	Int iHeight,
	Int iMinVal,
	Int iMaxVal,
	Int iOffset,
	Int iShift,
	ComponentID compID)
{
	Int iRow, iCol, iRemWd;
	Pel const *pSrcTmp;
	Pel const *inp_copy;
	Pel *dst_copy;

	__m128i coeff0_1_8x16b, coeff2_3_8x16b, coeff4_5_8x16b, coeff6_7_8x16b, mmMask;
	__m128i s0_8x16b, s1_8x16b, s2_8x16b, s3_8x16b, s4_8x16b, s5_8x16b, s6_8x16b, s7_8x16b, s8_8x16b, s9_8x16b;
	__m128i s2_0_16x8b, s2_1_16x8b, s2_2_16x8b, s2_3_16x8b, s2_4_16x8b, s2_5_16x8b, s2_6_16x8b, s2_7_16x8b;
	__m128i s3_0_16x8b, s3_1_16x8b, s3_2_16x8b, s3_3_16x8b, s3_4_16x8b, s3_5_16x8b, s3_6_16x8b, s3_7_16x8b;

	__m128i mmMin = _mm_set1_epi16(iMinVal);
	__m128i mmMax = _mm_set1_epi16(iMaxVal);
	__m128i offset_8x16b = _mm_set1_epi32(iOffset); /* for offset addition */

	pSrcTmp = pSrc;
	iRemWd = iWidth;
	inp_copy = pSrc;
	dst_copy = pDst;

	/* load 8 8-bit coefficients and convert 8-bit into 16-bit  */
	coeff0_1_8x16b = _mm_loadu_si128((__m128i*)pCoeff);

	coeff2_3_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x55);
	coeff4_5_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0xaa);
	coeff6_7_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0xff);
	coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);

	if (iRemWd > 7)
	{
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			Int iCount = 0;
			for (iCol = iWidth; iCol > 7; iCol -= 8)
			{
				/*load 8 pixel values.*/
				s2_0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount));

				/*load 8 pixel values*/
				s2_1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));

				s3_0_16x8b = _mm_unpacklo_epi16(s2_0_16x8b, s2_1_16x8b);
				s3_4_16x8b = _mm_unpackhi_epi16(s2_0_16x8b, s2_1_16x8b);

				s0_8x16b = _mm_madd_epi16(s3_0_16x8b, coeff0_1_8x16b);
				s4_8x16b = _mm_madd_epi16(s3_4_16x8b, coeff0_1_8x16b);
				/*load 8 pixel values*/
				s2_2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 1) + iCount));

				/*load 8 pixel values*/
				s2_3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride * 3) + iCount));

				s3_1_16x8b = _mm_unpacklo_epi16(s2_2_16x8b, s2_3_16x8b);
				s3_5_16x8b = _mm_unpackhi_epi16(s2_2_16x8b, s2_3_16x8b);

				s1_8x16b = _mm_madd_epi16(s3_1_16x8b, coeff2_3_8x16b);
				s5_8x16b = _mm_madd_epi16(s3_5_16x8b, coeff2_3_8x16b);

				/*load 8 pixel values*/
				s2_4_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 2) + iCount));

				/*load 8 pixel values*/
				s2_5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride * 5) + iCount));

				s3_2_16x8b = _mm_unpacklo_epi16(s2_4_16x8b, s2_5_16x8b);
				s3_6_16x8b = _mm_unpackhi_epi16(s2_4_16x8b, s2_5_16x8b);

				s2_8x16b = _mm_madd_epi16(s3_2_16x8b, coeff4_5_8x16b);
				s6_8x16b = _mm_madd_epi16(s3_6_16x8b, coeff4_5_8x16b);

				/*load 8 pixel values*/
				s2_6_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride * 6) + iCount));

				/*load 8 pixel values*/
				s2_7_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride * 7) + iCount));

				s3_3_16x8b = _mm_unpacklo_epi16(s2_6_16x8b, s2_7_16x8b);
				s3_7_16x8b = _mm_unpackhi_epi16(s2_6_16x8b, s2_7_16x8b);

				s3_8x16b = _mm_madd_epi16(s3_3_16x8b, coeff6_7_8x16b);
				s7_8x16b = _mm_madd_epi16(s3_7_16x8b, coeff6_7_8x16b);

				s0_8x16b = _mm_add_epi32(s0_8x16b, s1_8x16b);
				s2_8x16b = _mm_add_epi32(s2_8x16b, s3_8x16b);
				s4_8x16b = _mm_add_epi32(s4_8x16b, s5_8x16b);
				s6_8x16b = _mm_add_epi32(s6_8x16b, s7_8x16b);
				s0_8x16b = _mm_add_epi32(s0_8x16b, s2_8x16b);
				s4_8x16b = _mm_add_epi32(s4_8x16b, s6_8x16b);

				s0_8x16b = _mm_add_epi32(s0_8x16b, offset_8x16b);
				s4_8x16b = _mm_add_epi32(s4_8x16b, offset_8x16b);

				s7_8x16b = _mm_srai_epi32(s0_8x16b, iShift);
				s8_8x16b = _mm_srai_epi32(s4_8x16b, iShift);

				/* i2_tmp = CLIP_U8(i2_tmp);*/
				s9_8x16b = _mm_packs_epi32(s7_8x16b, s8_8x16b);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(s9_8x16b, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					s9_8x16b = _mm_or_si128(_mm_and_si128(mmMask, s9_8x16b), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(s9_8x16b, mmMax);
					s9_8x16b = _mm_or_si128(_mm_and_si128(mmMask, s9_8x16b), _mm_andnot_si128(mmMask, mmMax));
				}

				_mm_storeu_si128((__m128i*)(dst_copy + iCount), s9_8x16b);

				iCount += 8;
			}
			inp_copy += (iSrcStride);
			dst_copy += (iDstStride);
		}
	}

	iRemWd &= 0x7;
	inp_copy = pSrcTmp + ((iWidth / 8) * 8);
	dst_copy = pDst + ((iWidth / 8) * 8);

	if (iRemWd > 3)
	{
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			/*load 8 pixel values */
			s2_0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy));

			/*load 8 pixel values */
			s2_1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (iSrcStride)));

			s3_0_16x8b = _mm_unpacklo_epi16(s2_0_16x8b, s2_1_16x8b);

			s0_8x16b = _mm_madd_epi16(s3_0_16x8b, coeff0_1_8x16b);
			/*load 8 pixel values*/
			s2_2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (2 * iSrcStride)));

			/*load 8 pixel values*/
			s2_3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (3 * iSrcStride)));

			s3_1_16x8b = _mm_unpacklo_epi16(s2_2_16x8b, s2_3_16x8b);

			s1_8x16b = _mm_madd_epi16(s3_1_16x8b, coeff2_3_8x16b);
			/*load 8 pixel values*/
			s2_4_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (4 * iSrcStride)));

			/*load 8 pixel values*/
			s2_5_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (5 * iSrcStride)));

			s3_2_16x8b = _mm_unpacklo_epi16(s2_4_16x8b, s2_5_16x8b);

			s2_8x16b = _mm_madd_epi16(s3_2_16x8b, coeff4_5_8x16b);
			/*load 8 pixel values*/
			s2_6_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (6 * iSrcStride)));

			/*load 8 pixel values*/
			s2_7_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + (7 * iSrcStride)));

			s3_3_16x8b = _mm_unpacklo_epi16(s2_6_16x8b, s2_7_16x8b);

			s3_8x16b = _mm_madd_epi16(s3_3_16x8b, coeff6_7_8x16b);

			s4_8x16b = _mm_add_epi32(s0_8x16b, s1_8x16b);
			s5_8x16b = _mm_add_epi32(s2_8x16b, s3_8x16b);
			s6_8x16b = _mm_add_epi32(s4_8x16b, s5_8x16b);

			s7_8x16b = _mm_add_epi32(s6_8x16b, offset_8x16b);

			/*(i2_tmp + OFFSET_14_MINUS_BIT_DEPTH) >> SHIFT_14_MINUS_BIT_DEPTH */
			s8_8x16b = _mm_srai_epi32(s7_8x16b, iShift);

			/* i2_tmp = CLIP_U8(i2_tmp);*/
			s9_8x16b = _mm_packs_epi32(s8_8x16b, s8_8x16b);

			if (isLast)
			{
				mmMask = _mm_cmpgt_epi16(s9_8x16b, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
				s9_8x16b = _mm_or_si128(_mm_and_si128(mmMask, s9_8x16b), _mm_andnot_si128(mmMask, mmMin));
				mmMask = _mm_cmplt_epi16(s9_8x16b, mmMax);
				s9_8x16b = _mm_or_si128(_mm_and_si128(mmMask, s9_8x16b), _mm_andnot_si128(mmMask, mmMax));
			}
			_mm_storel_epi64((__m128i*)(dst_copy), s9_8x16b);

			inp_copy += (iSrcStride);
			dst_copy += (iDstStride);
		}
	}

	iRemWd &= 0x3;
	inp_copy = pSrcTmp + ((iWidth / 4) * 4);
	dst_copy = pDst + ((iWidth / 4) * 4);

	if (iRemWd)
	{
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iRemWd; iCol++)
			{
				Int sum;

				sum = inp_copy[iCol + 0 * iSrcStride] * pCoeff[0];
				sum += inp_copy[iCol + 1 * iSrcStride] * pCoeff[1];
				sum += inp_copy[iCol + 2 * iSrcStride] * pCoeff[2];
				sum += inp_copy[iCol + 3 * iSrcStride] * pCoeff[3];
				sum += inp_copy[iCol + 4 * iSrcStride] * pCoeff[4];
				sum += inp_copy[iCol + 5 * iSrcStride] * pCoeff[5];
				sum += inp_copy[iCol + 6 * iSrcStride] * pCoeff[6];
				sum += inp_copy[iCol + 7 * iSrcStride] * pCoeff[7];

				Pel val = (sum + iOffset) >> iShift;
				if (isLast)
				{
					val = ClipA(val, compID);
				}
				dst_copy[iCol] = val;
			}

			inp_copy += iSrcStride;
			dst_copy += iDstStride;
		}
	}
}

template < Bool isLast >
Void horzInterpolate4Tap
(
	Pel const *pSrc,
	Int iSrcStride,
	Pel *pDst,
	Int iDstStride,
	Short *pCoeff,
	Int iWidth,
	Int iHeight,
	Int iMinVal,
	Int iMaxVal,
	Int iOffset,
	Int iShift,
	ComponentID compID)
{
	assert(!(iWidth & 1));	/*width has to be even*/

	Int iRow, iCol, iRemWd, iRemHt, iCount;
	Int iSrcStride2, iSrcStride3;
	Pel const *inp_copy;
	Pel *dst_copy;

	/* all 128 bit registers are named with a suffix mxnb, where m is the */
	/* number of n bits packed in the register                            */

	__m128i offset_4x32b = _mm_set1_epi32(iOffset);
	__m128i	mmMin = _mm_set1_epi16(iMinVal);
	__m128i mmMax = _mm_set1_epi16(iMaxVal);
	__m128i coeff0_1_8x16b, coeff2_3_8x16b, coeff_8x16b, mmMask, shuffle;

	__m128i row1, row2, row3, row4;
	__m128i res0, res1, res2, res3;
	__m128i row11, row12, row13, row14, row21, row22, row23, row24;
	__m128i row31, row32, row33, row34, row41, row42, row43, row44;

	iSrcStride2 = (iSrcStride << 1);
	iSrcStride3 = (iSrcStride * 3);

	iRemWd = iWidth;
	inp_copy = pSrc;
	dst_copy = pDst;

	/* load 8 8-bit coefficients and convert 8-bit into 16-bit  */
	coeff0_1_8x16b = _mm_loadu_si128((__m128i*)pCoeff);
	shuffle = _mm_loadu_si128((__m128i*)shuffle_4Tap);

	if ((2 == iWidth) || (2 == iHeight))		/*height can be either 2 or 5 (2+filtersize-1)*/
	{
		coeff_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x44);  /*w0 w1 w2 w3 w0 w1 w2 w3*/

		if (iHeight >= 4)	/*height greater than or eaual to 4*/
		{
			iRemHt = iHeight % 4;

			for (iRow = 0; iRow < (iHeight - 3); iRow += 4)
			{

				/*load pixel values from row 1*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));		/*a0 a1 a2 a3 a4 x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));		/*b0 b1 b2 b3 b4 x x x*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));
				row4 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3));

				row11 = _mm_shuffle_epi8(row1, shuffle);			/*a0 a1 a2 a3 a1 a2 a3 a4*/
				row22 = _mm_shuffle_epi8(row2, shuffle);			/*b0 b1 b2 b3 b1 b2 b3 b4*/
				row33 = _mm_shuffle_epi8(row3, shuffle);
				row44 = _mm_shuffle_epi8(row4, shuffle);

				row11 = _mm_madd_epi16(row11, coeff_8x16b);		/*a0+a1 a2+a3 a1+a2 a3+a4*/
				row22 = _mm_madd_epi16(row22, coeff_8x16b);
				row33 = _mm_madd_epi16(row33, coeff_8x16b);
				row44 = _mm_madd_epi16(row44, coeff_8x16b);

				row11 = _mm_hadd_epi32(row11, row22);	/*a0+a1+a2+a3 a1+a2+a3+a4  b0+b1+b2+b3 b1+b2+b3+b4*/
				row33 = _mm_hadd_epi32(row33, row44);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row33 = _mm_add_epi32(row33, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row33 = _mm_srai_epi32(row33, iShift);

				row11 = _mm_packs_epi32(row11, row33);	/*1st 2nd 3rd 4th rows*/

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(row11, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(row11, mmMax);
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));

				}

				*((Int *)(dst_copy + 0 * iDstStride)) = (Int)_mm_cvtsi128_si32(row11);
				*((Int *)(dst_copy + iDstStride)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 4));
				*((Int *)(dst_copy + iDstStride * 2)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 8));
				*((Int *)(dst_copy + iDstStride * 3)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 12));

				inp_copy += (iSrcStride << 2);
				dst_copy += (iDstStride << 2);
			}

			if (1 == iRemHt)
			{
				/*load pixel values from row 1*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));	/*a0 a1 a2 a3 a4 x x x*/
				row11 = _mm_shuffle_epi8(row1, shuffle);		/*a0 a1 a2 a3 a1 a2 a3 a4*/
				row11 = _mm_madd_epi16(row11, coeff_8x16b);		/*a0+a1 a2+a3 a1+a2 a3+a4*/
				row11 = _mm_hadd_epi32(row11, row11);			/*a0+a1+a2+a3 a1+a2+a3+a4*/
				row11 = _mm_add_epi32(row11, offset_4x32b);
				row11 = _mm_srai_epi32(row11, iShift);
				row11 = _mm_packs_epi32(row11, row11);	/*1st*/

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(row11, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(row11, mmMax);
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));
				}

				*((Int *)(dst_copy)) = (Int)_mm_cvtsi128_si32(row11);
			}

			if (3 == iRemHt)
			{
				/*load pixel values from row 1*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));		/*a0 a1 a2 a3 a4 x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));		/*b0 b1 b2 b3 b4 x x x*/
				row3 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2));

				row11 = _mm_shuffle_epi8(row1, shuffle);			/*a0 a1 a2 a3 a1 a2 a3 a4*/
				row22 = _mm_shuffle_epi8(row2, shuffle);			/*b0 b1 b2 b3 b1 b2 b3 b4*/
				row33 = _mm_shuffle_epi8(row3, shuffle);

				row11 = _mm_madd_epi16(row11, coeff_8x16b);		/*a0+a1 a2+a3 a1+a2 a3+a4*/
				row22 = _mm_madd_epi16(row22, coeff_8x16b);
				row33 = _mm_madd_epi16(row33, coeff_8x16b);

				row11 = _mm_hadd_epi32(row11, row22);	/*a0+a1+a2+a3 a1+a2+a3+a4  b0+b1+b2+b3 b1+b2+b3+b4*/
				row33 = _mm_hadd_epi32(row33, row33);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row33 = _mm_add_epi32(row33, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row33 = _mm_srai_epi32(row33, iShift);

				row11 = _mm_packs_epi32(row11, row33);	/*1st 2nd 3rd 4th rows*/

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(row11, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(row11, mmMax);
					row11 = _mm_or_si128(_mm_and_si128(mmMask, row11), _mm_andnot_si128(mmMask, mmMax));
				}

				*((Int *)(dst_copy)) = (Int)_mm_cvtsi128_si32(row11);
				*((Int *)(dst_copy + iDstStride)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 4));
				*((Int *)(dst_copy + iDstStride * 2)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row11, 8));
			}
		}
		else
		{

			coeff2_3_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x55);  /*w2 w3 w2 w3 w2 w3 w2 w3*/
			coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);		/*w0 w1 w0 w1 w0 w1 w0 w1*/

			if (iRemWd > 7)
			{
				iCount = 0;
				for (iCol = iWidth; iCol > 7; iCol -= 8)
				{
					/*load pixel values from row 1*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row12 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
					row13 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
					row14 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																					/*load pixel values from row 2*/
					row21 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));
					row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));
					row23 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 2));
					row24 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 3));

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);    /*a0+a1 a2+a3 a4+a5 a6+a7*/
					row12 = _mm_madd_epi16(row12, coeff0_1_8x16b);	   /*a1+a2 a3+a4 a5+a6 a7+a8*/
					row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);	   /*a2+a3 a4+a5 a6+a7 a8+a9*/
					row14 = _mm_madd_epi16(row14, coeff2_3_8x16b);	   /*a3+a4 a5+a6 a7+a8 a9+a10*/
					row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
					row24 = _mm_madd_epi16(row24, coeff2_3_8x16b);

					row11 = _mm_add_epi32(row11, row13);
					row12 = _mm_add_epi32(row12, row14);
					row21 = _mm_add_epi32(row21, row23);
					row22 = _mm_add_epi32(row22, row24);

					row11 = _mm_add_epi32(row11, offset_4x32b);
					row12 = _mm_add_epi32(row12, offset_4x32b);
					row21 = _mm_add_epi32(row21, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);

					row11 = _mm_srai_epi32(row11, iShift);
					row12 = _mm_srai_epi32(row12, iShift);
					row21 = _mm_srai_epi32(row21, iShift);
					row22 = _mm_srai_epi32(row22, iShift);

					row11 = _mm_packs_epi32(row11, row21);
					row12 = _mm_packs_epi32(row12, row22);

					res0 = _mm_unpacklo_epi16(row11, row12);
					res1 = _mm_unpackhi_epi16(row11, row12);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
					}
					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);

					iCount += 8;
				}
			}

			iRemWd &= 0x7;
			inp_copy = pSrc + ((iWidth / 8) * 8);
			dst_copy = pDst + ((iWidth / 8) * 8);

			if (iRemWd > 3)
			{
				/*load pixel values from row 1*/
				row11 = _mm_loadl_epi64((__m128i*)(inp_copy));			/*a0 a1 a2 a3 */
				row12 = _mm_loadl_epi64((__m128i*)(inp_copy + 1));		/*a1 a2 a3 a4 */
				row13 = _mm_loadl_epi64((__m128i*)(inp_copy + 2));       /*a2 a3 a4 a5 */
				row14 = _mm_loadl_epi64((__m128i*)(inp_copy + 3));		/*a3 a4 a5 a6 */
																		/*load pixel values from row 2*/
				row21 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride));
				row22 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 1));
				row23 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 2));
				row24 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 3));

				row11 = _mm_unpacklo_epi32(row11, row12);
				row13 = _mm_unpacklo_epi32(row13, row14);
				row21 = _mm_unpacklo_epi32(row21, row22);
				row23 = _mm_unpacklo_epi32(row23, row24);

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);
				row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);
				row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
				row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);

				row11 = _mm_add_epi32(row11, row13);
				row21 = _mm_add_epi32(row21, row23);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row21 = _mm_add_epi32(row21, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row21 = _mm_srai_epi32(row21, iShift);

				res0 = _mm_packs_epi32(row11, row21);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));
				}
				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
			}

			iRemWd &= 0x3;
			inp_copy = pSrc + ((iWidth / 4) * 4);
			dst_copy = pDst + ((iWidth / 4) * 4);

			if (iRemWd > 1)	/*width = 2 case*/
			{

				/*load pixel values from row 1*/
				row1 = _mm_loadu_si128((__m128i*)(inp_copy));		/*a0 a1 a2 a3 a4 x x x*/
				row2 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride));		/*b0 b1 b2 b3 b4 x x x*/

				row11 = _mm_shuffle_epi8(row1, shuffle);			/*a0 a1 a2 a3 a1 a2 a3 a4*/
				row22 = _mm_shuffle_epi8(row2, shuffle);			/*b0 b1 b2 b3 b1 b2 b3 b4*/

				row11 = _mm_madd_epi16(row11, coeff_8x16b);		/*a0+a1 a2+a3 a1+a2 a3+a4*/
				row22 = _mm_madd_epi16(row22, coeff_8x16b);

				row1 = _mm_hadd_epi32(row11, row22);	/*a0+a1+a2+a3 a1+a2+a3+a4  b0+b1+b2+b3 b1+b2+b3+b4*/

				row2 = _mm_add_epi32(row1, offset_4x32b);
				row3 = _mm_srai_epi32(row2, iShift);
				row4 = _mm_packs_epi32(row3, row3);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(row4, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					row4 = _mm_or_si128(_mm_and_si128(mmMask, row4), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(row4, mmMax);
					row4 = _mm_or_si128(_mm_and_si128(mmMask, row4), _mm_andnot_si128(mmMask, mmMax));
				}

				*((Int *)(dst_copy)) = (Int)_mm_cvtsi128_si32(row4);
				*((Int *)(dst_copy + iDstStride)) = (Int)_mm_cvtsi128_si32(_mm_srli_si128(row4, 4));
			}
		}
	}
	else if (5 == iHeight)
	{
		__m128i row54, row51, row52, row53, res4;
		iRemWd = iWidth;
		coeff2_3_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x55);  /*w2 w3 w2 w3 w2 w3 w2 w3*/
		coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);		/*w0 w1 w0 w1 w0 w1 w0 w1*/

		iCount = 0;
		if (iRemWd > 7)
		{

			for (iCol = iWidth; iCol > 7; iCol -= 8)
			{
				/*load pixel values from row 1*/
				row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
				row12 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
				row13 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
				row14 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																				/*load pixel values from row 2*/
				row21 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));
				row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));
				row23 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 2));
				row24 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 3));

				/*load pixel values from row 3*/
				row31 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
				row32 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 1));
				row33 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 2));
				row34 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 3));
				/*load pixel values from row 4*/
				row41 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
				row42 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 1));
				row43 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 2));
				row44 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 3));
				/*load pixel values from row 5*/
				row51 = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 2) + iCount));
				row52 = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 2) + iCount + 1));
				row53 = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 2) + iCount + 2));
				row54 = _mm_loadu_si128((__m128i*)(inp_copy + (iSrcStride << 2) + iCount + 3));

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);    /*a0+a1 a2+a3 a4+a5 a6+a7*/
				row12 = _mm_madd_epi16(row12, coeff0_1_8x16b);	   /*a1+a2 a3+a4 a5+a6 a7+a8*/
				row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);	   /*a2+a3 a4+a5 a6+a7 a8+a9*/
				row14 = _mm_madd_epi16(row14, coeff2_3_8x16b);	   /*a3+a4 a5+a6 a7+a8 a9+a10*/
				row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
				row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
				row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
				row24 = _mm_madd_epi16(row24, coeff2_3_8x16b);
				row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);
				row32 = _mm_madd_epi16(row32, coeff0_1_8x16b);
				row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);
				row34 = _mm_madd_epi16(row34, coeff2_3_8x16b);
				row41 = _mm_madd_epi16(row41, coeff0_1_8x16b);
				row42 = _mm_madd_epi16(row42, coeff0_1_8x16b);
				row43 = _mm_madd_epi16(row43, coeff2_3_8x16b);
				row44 = _mm_madd_epi16(row44, coeff2_3_8x16b);
				row51 = _mm_madd_epi16(row51, coeff0_1_8x16b);
				row52 = _mm_madd_epi16(row52, coeff0_1_8x16b);
				row53 = _mm_madd_epi16(row53, coeff2_3_8x16b);
				row54 = _mm_madd_epi16(row54, coeff2_3_8x16b);

				row11 = _mm_add_epi32(row11, row13);
				row12 = _mm_add_epi32(row12, row14);
				row21 = _mm_add_epi32(row21, row23);
				row22 = _mm_add_epi32(row22, row24);
				row31 = _mm_add_epi32(row31, row33);
				row32 = _mm_add_epi32(row32, row34);
				row41 = _mm_add_epi32(row41, row43);
				row42 = _mm_add_epi32(row42, row44);
				row51 = _mm_add_epi32(row51, row53);
				row52 = _mm_add_epi32(row52, row54);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row12 = _mm_add_epi32(row12, offset_4x32b);
				row21 = _mm_add_epi32(row21, offset_4x32b);
				row22 = _mm_add_epi32(row22, offset_4x32b);
				row31 = _mm_add_epi32(row31, offset_4x32b);
				row32 = _mm_add_epi32(row32, offset_4x32b);
				row41 = _mm_add_epi32(row41, offset_4x32b);
				row42 = _mm_add_epi32(row42, offset_4x32b);
				row51 = _mm_add_epi32(row51, offset_4x32b);
				row52 = _mm_add_epi32(row52, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row12 = _mm_srai_epi32(row12, iShift);
				row21 = _mm_srai_epi32(row21, iShift);
				row22 = _mm_srai_epi32(row22, iShift);
				row31 = _mm_srai_epi32(row31, iShift);
				row32 = _mm_srai_epi32(row32, iShift);
				row41 = _mm_srai_epi32(row41, iShift);
				row42 = _mm_srai_epi32(row42, iShift);
				row51 = _mm_srai_epi32(row51, iShift);	/*r0 r2 r4 r6*/
				row52 = _mm_srai_epi32(row52, iShift); /*r1 r3 r5 r7*/

				row11 = _mm_packs_epi32(row11, row21);
				row12 = _mm_packs_epi32(row12, row22);
				row31 = _mm_packs_epi32(row31, row41);
				row32 = _mm_packs_epi32(row32, row42);
				row51 = _mm_packs_epi32(row51, row51);	/*r0 r2 r4 r6*/
				row52 = _mm_packs_epi32(row52, row52);	/*r1 r3 r5 r7*/

				res0 = _mm_unpacklo_epi16(row11, row12);
				res1 = _mm_unpackhi_epi16(row11, row12);
				res2 = _mm_unpacklo_epi16(row31, row32);
				res3 = _mm_unpackhi_epi16(row31, row32);
				res4 = _mm_unpacklo_epi16(row51, row52);	/*r0 r1 r2 r3 r4 r5 r6 r7*/

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res2, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res2, mmMax);
					res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res3, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res3, mmMax);
					res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res4, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res4 = _mm_or_si128(_mm_and_si128(mmMask, res4), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res4, mmMax);
					res4 = _mm_or_si128(_mm_and_si128(mmMask, res4), _mm_andnot_si128(mmMask, mmMax));
				}
				/* to store the 8 pixels res. */
				_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
				_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
				_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride << 1) + iCount), res2);
				_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride * 3) + iCount), res3);
				_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride << 2) + iCount), res4);

				iCount += 8;
			}
		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			/*load pixel values from row 1*/
			row11 = _mm_loadl_epi64((__m128i*)(inp_copy));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
			row12 = _mm_loadl_epi64((__m128i*)(inp_copy + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
			row13 = _mm_loadl_epi64((__m128i*)(inp_copy + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
			row14 = _mm_loadl_epi64((__m128i*)(inp_copy + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																	/*load pixel values from row 2*/
			row21 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride));
			row22 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 1));
			row23 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 2));
			row24 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 3));

			/*load pixel values from row 3*/
			row31 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2));
			row32 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 1));
			row33 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 2));
			row34 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 3));
			/*load pixel values from row 4*/
			row41 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3));
			row42 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 1));
			row43 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 2));
			row44 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 3));
			/*load pixel values from row 5*/
			row51 = _mm_loadl_epi64((__m128i*)(inp_copy + (iSrcStride << 2)));
			row52 = _mm_loadl_epi64((__m128i*)(inp_copy + (iSrcStride << 2) + 1));
			row53 = _mm_loadl_epi64((__m128i*)(inp_copy + (iSrcStride << 2) + 2));
			row54 = _mm_loadl_epi64((__m128i*)(inp_copy + (iSrcStride << 2) + 3));

			row11 = _mm_unpacklo_epi32(row11, row12);
			row13 = _mm_unpacklo_epi32(row13, row14);
			row21 = _mm_unpacklo_epi32(row21, row22);
			row23 = _mm_unpacklo_epi32(row23, row24);
			row31 = _mm_unpacklo_epi32(row31, row32);
			row33 = _mm_unpacklo_epi32(row33, row34);
			row41 = _mm_unpacklo_epi32(row41, row42);
			row43 = _mm_unpacklo_epi32(row43, row44);
			row51 = _mm_unpacklo_epi32(row51, row52);
			row53 = _mm_unpacklo_epi32(row53, row54);

			row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);
			row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);
			row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
			row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
			row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);
			row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);
			row41 = _mm_madd_epi16(row41, coeff0_1_8x16b);
			row43 = _mm_madd_epi16(row43, coeff2_3_8x16b);
			row51 = _mm_madd_epi16(row51, coeff0_1_8x16b);
			row53 = _mm_madd_epi16(row53, coeff2_3_8x16b);

			row11 = _mm_add_epi32(row11, row13);
			row21 = _mm_add_epi32(row21, row23);
			row31 = _mm_add_epi32(row31, row33);
			row41 = _mm_add_epi32(row41, row43);
			row51 = _mm_add_epi32(row51, row53);

			row11 = _mm_add_epi32(row11, offset_4x32b);
			row21 = _mm_add_epi32(row21, offset_4x32b);
			row31 = _mm_add_epi32(row31, offset_4x32b);
			row41 = _mm_add_epi32(row41, offset_4x32b);
			row51 = _mm_add_epi32(row51, offset_4x32b);

			row11 = _mm_srai_epi32(row11, iShift);
			row21 = _mm_srai_epi32(row21, iShift);
			row31 = _mm_srai_epi32(row31, iShift);
			row41 = _mm_srai_epi32(row41, iShift);
			row51 = _mm_srai_epi32(row51, iShift);

			res0 = _mm_packs_epi32(row11, row21);
			res1 = _mm_packs_epi32(row31, row41);
			res2 = _mm_packs_epi32(row51, row51);

			if (isLast)
			{
				mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
				res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
				mmMask = _mm_cmplt_epi16(res0, mmMax);
				res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

				mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
				res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
				mmMask = _mm_cmplt_epi16(res1, mmMax);
				res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

				mmMask = _mm_cmpgt_epi16(res2, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
				res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
				mmMask = _mm_cmplt_epi16(res2, mmMax);
				res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));
			}
			/* to store the 8 pixels res. */
			_mm_storel_epi64((__m128i *)(dst_copy), res0);
			_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
			_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride << 1)), res1);
			_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride * 3)), _mm_unpackhi_epi64(res1, res1));
			_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride << 2)), res2);
		}
	}
	else	/*width 4 or M4 ---- height 4 or 4+3*/
	{
		assert(!(iHeight == 5));

		iRemHt = iHeight % 4;
		iRemWd = iWidth;

		coeff2_3_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0x55);  /*w2 w3 w2 w3 w2 w3 w2 w3*/
		coeff0_1_8x16b = _mm_shuffle_epi32(coeff0_1_8x16b, 0);		/*w0 w1 w0 w1 w0 w1 w0 w1*/

		if (iRemWd > 7)
		{
			iCount = 0;
			for (iRow = 0; iRow < (iHeight - 3); iRow += 4)
			{
				for (iCol = iWidth; iCol > 7; iCol -= 8)
				{
					/*load pixel values from row 1*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row12 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
					row13 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
					row14 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																					/*load pixel values from row 2*/
					row21 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));
					row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));
					row23 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 2));
					row24 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 3));

					/*load pixel values from row 3*/
					row31 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row32 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 1));
					row33 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 2));
					row34 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 3));
					/*load pixel values from row 4*/
					row41 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount));
					row42 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 1));
					row43 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 2));
					row44 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride3 + iCount + 3));

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);    /*a0+a1 a2+a3 a4+a5 a6+a7*/
					row12 = _mm_madd_epi16(row12, coeff0_1_8x16b);	   /*a1+a2 a3+a4 a5+a6 a7+a8*/
					row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);	   /*a2+a3 a4+a5 a6+a7 a8+a9*/
					row14 = _mm_madd_epi16(row14, coeff2_3_8x16b);	   /*a3+a4 a5+a6 a7+a8 a9+a10*/
					row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
					row24 = _mm_madd_epi16(row24, coeff2_3_8x16b);
					row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);
					row32 = _mm_madd_epi16(row32, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);
					row34 = _mm_madd_epi16(row34, coeff2_3_8x16b);
					row41 = _mm_madd_epi16(row41, coeff0_1_8x16b);
					row42 = _mm_madd_epi16(row42, coeff0_1_8x16b);
					row43 = _mm_madd_epi16(row43, coeff2_3_8x16b);
					row44 = _mm_madd_epi16(row44, coeff2_3_8x16b);

					row11 = _mm_add_epi32(row11, row13);
					row12 = _mm_add_epi32(row12, row14);
					row21 = _mm_add_epi32(row21, row23);
					row22 = _mm_add_epi32(row22, row24);
					row31 = _mm_add_epi32(row31, row33);
					row32 = _mm_add_epi32(row32, row34);
					row41 = _mm_add_epi32(row41, row43);
					row42 = _mm_add_epi32(row42, row44);

					row11 = _mm_add_epi32(row11, offset_4x32b);
					row12 = _mm_add_epi32(row12, offset_4x32b);
					row21 = _mm_add_epi32(row21, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row31 = _mm_add_epi32(row31, offset_4x32b);
					row32 = _mm_add_epi32(row32, offset_4x32b);
					row41 = _mm_add_epi32(row41, offset_4x32b);
					row42 = _mm_add_epi32(row42, offset_4x32b);

					row11 = _mm_srai_epi32(row11, iShift);
					row12 = _mm_srai_epi32(row12, iShift);
					row21 = _mm_srai_epi32(row21, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row31 = _mm_srai_epi32(row31, iShift);
					row32 = _mm_srai_epi32(row32, iShift);
					row41 = _mm_srai_epi32(row41, iShift);
					row42 = _mm_srai_epi32(row42, iShift);

					row11 = _mm_packs_epi32(row11, row21);
					row12 = _mm_packs_epi32(row12, row22);
					row31 = _mm_packs_epi32(row31, row41);
					row32 = _mm_packs_epi32(row32, row42);

					res0 = _mm_unpacklo_epi16(row11, row12);
					res1 = _mm_unpackhi_epi16(row11, row12);
					res2 = _mm_unpacklo_epi16(row31, row32);
					res3 = _mm_unpackhi_epi16(row31, row32);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res2, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res3, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res3, mmMax);
						res3 = _mm_or_si128(_mm_and_si128(mmMask, res3), _mm_andnot_si128(mmMask, mmMax));
					}
					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride << 1) + iCount), res2);
					_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride * 3) + iCount), res3);

					iCount += 8;
				}

				iCount = 0;
				inp_copy += (iSrcStride << 2); /* pointer updates*/
				dst_copy += (iDstStride << 2); /* pointer updates*/
			}

			/*remaining 3 rows at a time*/
			if (3 == iRemHt)
			{
				iCount = 0;
				for (iCol = iWidth; iCol > 7; iCol -= 8)
				{
					/*load pixel values from row 1*/
					row11 = _mm_loadu_si128((__m128i*)(inp_copy + iCount));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
					row12 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
					row13 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
					row14 = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																					/*load pixel values from row 2*/
					row21 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount));
					row22 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 1));
					row23 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 2));
					row24 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride + iCount + 3));

					/*load pixel values from row 3*/
					row31 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount));
					row32 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 1));
					row33 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 2));
					row34 = _mm_loadu_si128((__m128i*)(inp_copy + iSrcStride2 + iCount + 3));

					row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);    /*a0+a1 a2+a3 a4+a5 a6+a7*/
					row12 = _mm_madd_epi16(row12, coeff0_1_8x16b);	   /*a1+a2 a3+a4 a5+a6 a7+a8*/
					row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);	   /*a2+a3 a4+a5 a6+a7 a8+a9*/
					row14 = _mm_madd_epi16(row14, coeff2_3_8x16b);	   /*a3+a4 a5+a6 a7+a8 a9+a10*/
					row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
					row22 = _mm_madd_epi16(row22, coeff0_1_8x16b);
					row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
					row24 = _mm_madd_epi16(row24, coeff2_3_8x16b);
					row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);
					row32 = _mm_madd_epi16(row32, coeff0_1_8x16b);
					row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);
					row34 = _mm_madd_epi16(row34, coeff2_3_8x16b);

					row11 = _mm_add_epi32(row11, row13);
					row12 = _mm_add_epi32(row12, row14);
					row21 = _mm_add_epi32(row21, row23);
					row22 = _mm_add_epi32(row22, row24);
					row31 = _mm_add_epi32(row31, row33);
					row32 = _mm_add_epi32(row32, row34);

					row11 = _mm_add_epi32(row11, offset_4x32b);
					row12 = _mm_add_epi32(row12, offset_4x32b);
					row21 = _mm_add_epi32(row21, offset_4x32b);
					row22 = _mm_add_epi32(row22, offset_4x32b);
					row31 = _mm_add_epi32(row31, offset_4x32b);
					row32 = _mm_add_epi32(row32, offset_4x32b);

					row11 = _mm_srai_epi32(row11, iShift);
					row12 = _mm_srai_epi32(row12, iShift);
					row21 = _mm_srai_epi32(row21, iShift);
					row22 = _mm_srai_epi32(row22, iShift);
					row31 = _mm_srai_epi32(row31, iShift);
					row32 = _mm_srai_epi32(row32, iShift);

					row11 = _mm_packs_epi32(row11, row21);
					row12 = _mm_packs_epi32(row12, row22);
					row31 = _mm_packs_epi32(row31, row32);

					res0 = _mm_unpacklo_epi16(row11, row12);
					res1 = _mm_unpackhi_epi16(row11, row12);
					res3 = _mm_unpackhi_epi64(row31, row31);
					res2 = _mm_unpacklo_epi16(row31, res3);

					if (isLast)
					{
						mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res0, mmMax);
						res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

						mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res1, mmMax);
						res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));


						mmMask = _mm_cmpgt_epi16(res2, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMin));
						mmMask = _mm_cmplt_epi16(res2, mmMax);
						res2 = _mm_or_si128(_mm_and_si128(mmMask, res2), _mm_andnot_si128(mmMask, mmMax));
					}

					/* to store the 8 pixels res. */
					_mm_storeu_si128((__m128i *)(dst_copy + iCount), res0);
					_mm_storeu_si128((__m128i *)(dst_copy + iDstStride + iCount), res1);
					_mm_storeu_si128((__m128i *)(dst_copy + (iDstStride << 1) + iCount), res2);

					iCount += 8;
				}
			}

		}

		iRemWd &= 0x7;
		inp_copy = pSrc + ((iWidth / 8) * 8);
		dst_copy = pDst + ((iWidth / 8) * 8);

		if (iRemWd > 3)
		{
			for (iRow = 0; iRow < (iHeight - 3); iRow += 4)
			{
				/*load pixel values from row 1*/
				row11 = _mm_loadl_epi64((__m128i*)(inp_copy));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
				row12 = _mm_loadl_epi64((__m128i*)(inp_copy + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
				row13 = _mm_loadl_epi64((__m128i*)(inp_copy + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
				row14 = _mm_loadl_epi64((__m128i*)(inp_copy + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																		/*load pixel values from row 2*/
				row21 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride));
				row22 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 1));
				row23 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 2));
				row24 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 3));

				/*load pixel values from row 3*/
				row31 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2));
				row32 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 1));
				row33 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 2));
				row34 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 3));
				/*load pixel values from row 4*/
				row41 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3));
				row42 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 1));
				row43 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 2));
				row44 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride3 + 3));

				row11 = _mm_unpacklo_epi32(row11, row12);
				row13 = _mm_unpacklo_epi32(row13, row14);
				row21 = _mm_unpacklo_epi32(row21, row22);
				row23 = _mm_unpacklo_epi32(row23, row24);
				row31 = _mm_unpacklo_epi32(row31, row32);
				row33 = _mm_unpacklo_epi32(row33, row34);
				row41 = _mm_unpacklo_epi32(row41, row42);
				row43 = _mm_unpacklo_epi32(row43, row44);

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);
				row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);
				row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);
				row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);
				row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);
				row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);
				row41 = _mm_madd_epi16(row41, coeff0_1_8x16b);
				row43 = _mm_madd_epi16(row43, coeff2_3_8x16b);

				row11 = _mm_add_epi32(row11, row13);
				row21 = _mm_add_epi32(row21, row23);
				row31 = _mm_add_epi32(row31, row33);
				row41 = _mm_add_epi32(row41, row43);

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row21 = _mm_add_epi32(row21, offset_4x32b);
				row31 = _mm_add_epi32(row31, offset_4x32b);
				row41 = _mm_add_epi32(row41, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row21 = _mm_srai_epi32(row21, iShift);
				row31 = _mm_srai_epi32(row31, iShift);
				row41 = _mm_srai_epi32(row41, iShift);

				res0 = _mm_packs_epi32(row11, row21);
				res1 = _mm_packs_epi32(row31, row41);

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
				}
				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride << 1)), res1);
				_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride * 3)), _mm_unpackhi_epi64(res1, res1));

				inp_copy += (iSrcStride << 2); /* pointer updates*/
				dst_copy += (iDstStride << 2); /* pointer updates*/
			}

			for (iRow = 0; iRow < iRemHt; iRow += 3)
			{
				/*load pixel values from row 1*/
				row11 = _mm_loadl_epi64((__m128i*)(inp_copy));			/*a0 a1 a2 a3 a4 a5 a6 a7*/
				row12 = _mm_loadl_epi64((__m128i*)(inp_copy + 1));		/*a1 a2 a3 a4 a5 a6 a7 a8*/
				row13 = _mm_loadl_epi64((__m128i*)(inp_copy + 2));       /*a2 a3 a4 a5 a6 a7 a8 a9*/
				row14 = _mm_loadl_epi64((__m128i*)(inp_copy + 3));		/*a3 a4 a5 a6 a7 a8 a9 a10*/
																		/*load pixel values from row 2*/
				row21 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride));
				row22 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 1));
				row23 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 2));
				row24 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride + 3));

				/*load pixel values from row 3*/
				row31 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2));		/*c0 c1 c2 c3*/
				row32 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 1));		/*c1 c2 c3 c4*/
				row33 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 2));		/*c2 c3 c4 c5*/
				row34 = _mm_loadl_epi64((__m128i*)(inp_copy + iSrcStride2 + 3));		/*c3 c4 c5 c6*/

				row11 = _mm_unpacklo_epi32(row11, row12);		/*a0 a1 a1 a2 a2 a3 a3 a4*/
				row13 = _mm_unpacklo_epi32(row13, row14);		/*a2 a3 a3 a4 a4 a5 a5 a6*/
				row21 = _mm_unpacklo_epi32(row21, row22);		/*b0 b1 b1 b2 b2 b3 b3 b4*/
				row23 = _mm_unpacklo_epi32(row23, row24);		/*b2 b3 b3 b4 b4 b5 b5 b6*/
				row31 = _mm_unpacklo_epi32(row31, row32);		/*c0 c1 c1 c2 c2 c3 c3 c4*/
				row33 = _mm_unpacklo_epi32(row33, row34);		/*c2 c3 c3 c4 c4 c5 c5 c6*/

				row11 = _mm_madd_epi16(row11, coeff0_1_8x16b);    /*a0+a1 a1+a2 a2+a3 a3+a4*/
				row13 = _mm_madd_epi16(row13, coeff2_3_8x16b);	   /*a2+a3 a3+a4 a4+a5 a5+a6*/
				row21 = _mm_madd_epi16(row21, coeff0_1_8x16b);	   /*b0+b1 b1+b2 b2+b3 b3+b4*/
				row23 = _mm_madd_epi16(row23, coeff2_3_8x16b);	   /*b2+b3 b3+b4 b4+b5 b5+b6*/
				row31 = _mm_madd_epi16(row31, coeff0_1_8x16b);	   /*c0+c1 c1+c2 c2+c3 c3+c4*/
				row33 = _mm_madd_epi16(row33, coeff2_3_8x16b);	   /*c2+c3 c3+c4 c4+c5 c5+c6*/

				row11 = _mm_add_epi32(row11, row13);	/*r00 r01  r02  r03*/
				row21 = _mm_add_epi32(row21, row23);	/*r10 r11  r12  r13*/
				row31 = _mm_add_epi32(row31, row33);   /*r20 r21  r22  r23*/

				row11 = _mm_add_epi32(row11, offset_4x32b);
				row21 = _mm_add_epi32(row21, offset_4x32b);
				row31 = _mm_add_epi32(row31, offset_4x32b);

				row11 = _mm_srai_epi32(row11, iShift);
				row21 = _mm_srai_epi32(row21, iShift);
				row31 = _mm_srai_epi32(row31, iShift);

				res0 = _mm_packs_epi32(row11, row21);	/*r00 r01  r02  r03	r10 r11  r12  r13*/
				res1 = _mm_packs_epi32(row31, row31);	/*r20 r21  r22  r23*/

				if (isLast)
				{
					mmMask = _mm_cmpgt_epi16(res0, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res0, mmMax);
					res0 = _mm_or_si128(_mm_and_si128(mmMask, res0), _mm_andnot_si128(mmMask, mmMax));

					mmMask = _mm_cmpgt_epi16(res1, mmMin);  /*if gt = -1...  -1 -1 0 0 -1 */
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMin));
					mmMask = _mm_cmplt_epi16(res1, mmMax);
					res1 = _mm_or_si128(_mm_and_si128(mmMask, res1), _mm_andnot_si128(mmMask, mmMax));
				}
				/* to store the 8 pixels res. */
				_mm_storel_epi64((__m128i *)(dst_copy), res0);
				_mm_storel_epi64((__m128i *)(dst_copy + iDstStride), _mm_unpackhi_epi64(res0, res0));
				_mm_storel_epi64((__m128i *)(dst_copy + (iDstStride << 1)), res1);
			}
		}
	}
}

template < Bool isLast >
Void vertInterpolate4Tap(Pel const *src, Int srcStride, Pel *dst,
	Int dstStride, Short *pCoeff, Int iWidth, Int iHeight,
	Int iMinVal, Int iMaxVal,
	Int iOffset, Int iShift, ComponentID compID)
{
	Int iRow, iCol;
	Int cStride = srcStride;
	Int N = 4;

	if (!(iWidth & 0x03))  /*if iWidth multiples of 4*/
	{
#if !JVET_D0033_ADAPTIVE_CLIPPING
		Short iMinVal = 0;
#endif
		__m128i mmOffset = _mm_set1_epi32(iOffset);
		__m128i mmCoeff[8];
		__m128i mmMin = _mm_set1_epi16(iMinVal);
		__m128i mmMax = _mm_set1_epi16(iMaxVal);
		for (Int n = 0; n < 4; n++)
			mmCoeff[n] = _mm_set1_epi16(pCoeff[n]);
		for (iRow = 0; iRow < iHeight; iRow++)
		{
			for (iCol = 0; iCol < iWidth; iCol += 4)
			{
				__m128i mmFiltered = simdInterpolateChroma4(src + iCol, cStride, mmCoeff, mmOffset, iShift);
				if (isLast)
				{
					mmFiltered = simdClip3(mmMin, mmMax, mmFiltered);
				}
				_mm_storel_epi64((__m128i *)(dst + iCol), mmFiltered);
			}
			src += srcStride;
			dst += dstStride;
		}
		return;
	}

	for (iRow = 0; iRow < iHeight; iRow++)
	{
		for (iCol = 0; iCol < iWidth; iCol++)
		{
			Int sum;

			sum = src[iCol + 0 * cStride] * pCoeff[0];
			sum += src[iCol + 1 * cStride] * pCoeff[1];
			if (N >= 4)
			{
				sum += src[iCol + 2 * cStride] * pCoeff[2];
				sum += src[iCol + 3 * cStride] * pCoeff[3];
			}
			if (N >= 6)
			{
				sum += src[iCol + 4 * cStride] * pCoeff[4];
				sum += src[iCol + 5 * cStride] * pCoeff[5];
			}
			if (N == 8)
			{
				sum += src[iCol + 6 * cStride] * pCoeff[6];
				sum += src[iCol + 7 * cStride] * pCoeff[7];
			}

			Pel val = (sum + iOffset) >> iShift;
			if (isLast)
			{
#if JVET_D0033_ADAPTIVE_CLIPPING
				val = ClipA(val, compID);
#else
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
#endif
			}
			dst[iCol] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}
#endif

// ====================================================================================================================
// Private member functions
// ====================================================================================================================

/**
 * \brief Apply unit FIR filter to a block of samples
 *
 * \param bitDepth   bitDepth of samples
 * \param src        Pointer to source samples
 * \param srcStride  Stride of source samples
 * \param dst        Pointer to destination samples
 * \param dstStride  Stride of destination samples
 * \param width      Width of block
 * \param height     Height of block
 * \param isFirst    Flag indicating whether it is the first filtering operation
 * \param isLast     Flag indicating whether it is the last filtering operation
 */
#if JVET_D0033_ADAPTIVE_CLIPPING
Void TComInterpolationFilter::filterCopy(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, ComponentID compID)
#else
Void TComInterpolationFilter::filterCopy(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast)
#endif
{
  Int row, col;

  if ( isFirst == isLast )
  {
	  twoDimBufCopier
		  (
			  (UChar *)dst,
			  dstStride * sizeof(Pel),
			  (UChar *)src,
			  srcStride * sizeof(Pel),
			  width * sizeof(Pel),
			  height
			  );
    //for (row = 0; row < height; row++)
    //{
    //  for (col = 0; col < width; col++)
    //  {
    //    dst[col] = src[col];
    //  }

    //  src += srcStride;
    //  dst += dstStride;
    //}
  }
  else if ( isFirst )
  {
    const Int shift = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));

	leftShiftSub
		(
			src,
			srcStride,
			dst,
			dstStride,
			width,
			height,
			shift,
			IF_INTERNAL_OFFS
			);
    //for (row = 0; row < height; row++)
    //{
    //  for (col = 0; col < width; col++)
    //  {
    //    Pel val = leftShift_round(src[col], shift);
    //    dst[col] = val - (Pel)IF_INTERNAL_OFFS;
    //  }

    //  src += srcStride;
    //  dst += dstStride;
    //}
  }
  else
  {
    const Int shift = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));

#if !JVET_D0033_ADAPTIVE_CLIPPING
    Pel maxVal = (1 << bitDepth) - 1;
    Pel minVal = 0;
#else
	Pel maxVal = g_ClipParam.max(compID);
	Pel minVal = g_ClipParam.min(compID);
#endif
	rightShiftClip
		(
			src,
			srcStride,
			dst,
			dstStride,
			width,
			height,
			shift,
			IF_INTERNAL_OFFS,
			minVal,
			maxVal,
			compID
			);
//    for (row = 0; row < height; row++)
//    {
//      for (col = 0; col < width; col++)
//      {
//        Pel val = src[ col ];
//        val = rightShift_round((val + IF_INTERNAL_OFFS), shift);
//#if JVET_D0033_ADAPTIVE_CLIPPING
//        val=ClipA(val,compID);
//#else
//        if (val < minVal)
//        {
//          val = minVal;
//        }
//        if (val > maxVal)
//        {
//          val = maxVal;
//        }
//#endif
//        dst[col] = val;
//      }
//
//      src += srcStride;
//      dst += dstStride;
//    }
  }
}

/**
 * \brief Apply FIR filter to a block of samples
 *
 * \tparam N          Number of taps
 * \tparam isVertical Flag indicating filtering along vertical direction
 * \tparam isFirst    Flag indicating whether it is the first filtering operation
 * \tparam isLast     Flag indicating whether it is the last filtering operation
 * \param  bitDepth   Bit depth of samples
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  coeff      Pointer to filter taps
 */
#if JVET_D0033_ADAPTIVE_CLIPPING
template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
Void TComInterpolationFilter::filter(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, TFilterCoeff const *coeff, ComponentID compID)
#else
template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
Void TComInterpolationFilter::filter(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, TFilterCoeff const *coeff)
#endif
{
  Int row, col;

  Pel c[8];
  c[0] = coeff[0];
  c[1] = coeff[1];
  if ( N >= 4 )
  {
    c[2] = coeff[2];
    c[3] = coeff[3];
  }
  if ( N >= 6 )
  {
    c[4] = coeff[4];
    c[5] = coeff[5];
  }
  if ( N == 8 )
  {
    c[6] = coeff[6];
    c[7] = coeff[7];
  }

  Int cStride = ( isVertical ) ? srcStride : 1;
  src -= ( N/2 - 1 ) * cStride;

  Int offset;
  Pel maxVal;
  Int headRoom = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));
  Int shift    = IF_FILTER_PREC;
  // with the current settings (IF_INTERNAL_PREC = 14 and IF_FILTER_PREC = 6), though headroom can be
  // negative for bit depths greater than 14, shift will remain non-negative for bit depths of 8->20
  assert(shift >= 0);

#if JVET_D0033_ADAPTIVE_CLIPPING
  Short minVal=0;
#endif
  if ( isLast )
  {
    shift += (isFirst) ? 0 : headRoom;
    offset = 1 << (shift - 1);
    offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
#if JVET_D0033_ADAPTIVE_CLIPPING
    minVal = g_ClipParam.min(compID);
    maxVal = g_ClipParam.max(compID);
#else
    maxVal = (1 << bitDepth) - 1;
#endif
  }
  else
  {
    shift -= (isFirst) ? headRoom : 0;
    offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    maxVal = 0;
  }

//#if COM16_C806_SIMD_OPT
//  if( bitDepth <= 10 )
//  {
//    if( N == 8 && !( width & 0x07 ) )
//    {
//#if !JVET_D0033_ADAPTIVE_CLIPPING
//      Short minVal = 0;
//#endif
//      __m128i mmOffset = _mm_set1_epi32( offset );
//      __m128i mmCoeff[8];
//      __m128i mmMin = _mm_set1_epi16( minVal );
//      __m128i mmMax = _mm_set1_epi16( maxVal );
//      for( Int n = 0 ; n < 8 ; n++ )
//        mmCoeff[n] = _mm_set1_epi16( c[n] );
//      for( row = 0 ; row < height ; row++ )
//      {
//        for( col = 0 ; col < width ; col += 8 )
//        {
//          __m128i mmFiltered = simdInterpolateLuma8( src + col , cStride , mmCoeff , mmOffset , shift );
//          if( isLast )
//          {
//            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
//          }
//          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
//        }
//        src += srcStride;
//        dst += dstStride;
//      }
//      return;
//    }
//    else if( N == 8 && !( width & 0x03 ) )
//    {
//#if !JVET_D0033_ADAPTIVE_CLIPPING
//      Short minVal = 0;
//#endif
//      __m128i mmOffset = _mm_set1_epi32( offset );
//      __m128i mmCoeff[8];
//      __m128i mmMin = _mm_set1_epi16( minVal );
//      __m128i mmMax = _mm_set1_epi16( maxVal );
//      for( Int n = 0 ; n < 8 ; n++ )
//        mmCoeff[n] = _mm_set1_epi16( c[n] );
//      for( row = 0 ; row < height ; row++ )
//      {
//        for( col = 0 ; col < width ; col += 4 )
//        {
//          __m128i mmFiltered = simdInterpolateLuma4( src + col , cStride , mmCoeff , mmOffset , shift );
//          if( isLast )
//          {
//            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
//          }
//          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
//        }
//        src += srcStride;
//        dst += dstStride;
//      }
//      return;
//    }
//    else if( N == 4 && !( width & 0x03 ) )
//    {
//#if !JVET_D0033_ADAPTIVE_CLIPPING
//      Short minVal = 0;
//#endif
//      __m128i mmOffset = _mm_set1_epi32( offset );
//      __m128i mmCoeff[8];
//      __m128i mmMin = _mm_set1_epi16( minVal );
//      __m128i mmMax = _mm_set1_epi16( maxVal );
//      for( Int n = 0 ; n < 4 ; n++ )
//        mmCoeff[n] = _mm_set1_epi16( c[n] );
//      for( row = 0 ; row < height ; row++ )
//      {
//        for( col = 0 ; col < width ; col += 4 )
//        {
//          __m128i mmFiltered = simdInterpolateChroma4( src + col , cStride , mmCoeff , mmOffset , shift );
//          if( isLast )
//          {
//            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
//          }
//          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
//        }
//        src += srcStride;
//        dst += dstStride;
//      }
//      return;
//    }
//    else if( N == 2 && !( width & 0x07 ) )
//    {
//#if !JVET_D0033_ADAPTIVE_CLIPPING
//      Short minVal = 0;
//#endif
//      __m128i mmOffset = _mm_set1_epi32( offset );
//      __m128i mmCoeff[2];
//      __m128i mmMin = _mm_set1_epi16( minVal );
//      __m128i mmMax = _mm_set1_epi16( maxVal );
//      for( Int n = 0 ; n < 2 ; n++ )
//        mmCoeff[n] = _mm_set1_epi16( c[n] );
//      for( row = 0 ; row < height ; row++ )
//      {
//        for( col = 0 ; col < width ; col += 8 )
//        {
//          __m128i mmFiltered = simdInterpolateLuma2P8( src + col , cStride , mmCoeff , mmOffset , shift );
//          if( isLast )
//          {
//            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
//          }
//          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
//        }
//        src += srcStride;
//        dst += dstStride;
//      }
//      return;
//    }
//    else if( N == 2 && !( width & 0x03 ) )
//    {
//#if !JVET_D0033_ADAPTIVE_CLIPPING
//      Short minVal = 0;
//#endif
//      __m128i mmOffset = _mm_set1_epi32( offset );
//      __m128i mmCoeff[8];
//      __m128i mmMin = _mm_set1_epi16( minVal );
//      __m128i mmMax = _mm_set1_epi16( maxVal );
//      for( Int n = 0 ; n < 2 ; n++ )
//        mmCoeff[n] = _mm_set1_epi16( c[n] );
//      for( row = 0 ; row < height ; row++ )
//      {
//        for( col = 0 ; col < width ; col += 4 )
//        {
//          __m128i mmFiltered = simdInterpolateLuma2P4( src + col , cStride , mmCoeff , mmOffset , shift );
//          if( isLast )
//          {
//            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
//          }
//          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
//        }
//        src += srcStride;
//        dst += dstStride;
//      }
//      return;
//    }
//  }
//#endif
#if (ITH_C3_INTERP_FILT==SIMD_FXN) || (ITH_C3_INTERP_FILT==COARSE_SIMD_FXN)
  if (bitDepth <= 10)
  {
	  if (N == 8)
	  {
		  if (isVertical)
		  {
			  vertInterpolate8Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
		  else
		  {
			  horzInterpolate8Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
	  }
	  else if (N == 4)
	  {
		  if (isVertical)
		  {
			  vertInterpolate4Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
		  else
		  {
			  horzInterpolate4Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
	  }
	  else if (N == 2)
	  {
		  if (isVertical)
		  {
			  vertInterpolate2Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
		  else
		  {
			  horzInterpolate2Tap< isLast >(src, srcStride, dst, dstStride, c, width, height,
				  minVal, maxVal, offset, shift, compID);
		  }
	  }
  }
  else
  {
	  for (Int row = 0; row < height; row++)
	  {
		  for (Int col = 0; col < width; col++)
		  {
			  Int sum;

			  sum = src[col + 0 * cStride] * c[0];
			  sum += src[col + 1 * cStride] * c[1];
			  if (N >= 4)
			  {
				  sum += src[col + 2 * cStride] * c[2];
				  sum += src[col + 3 * cStride] * c[3];
			  }
			  if (N >= 6)
			  {
				  sum += src[col + 4 * cStride] * c[4];
				  sum += src[col + 5 * cStride] * c[5];
			  }
			  if (N == 8)
			  {
				  sum += src[col + 6 * cStride] * c[6];
				  sum += src[col + 7 * cStride] * c[7];
			  }

			  Pel val = (sum + offset) >> shift;
			  if (isLast)
			  {
#if JVET_D0033_ADAPTIVE_CLIPPING
				  val = ClipA(val, compID);
#else
				  val = (val < 0) ? 0 : val;
				  val = (val > maxVal) ? maxVal : val;
#endif
			  }
			  dst[col] = val;
		  }

		  src += srcStride;
		  dst += dstStride;
	  }
  }
#elif ITH_C3_INTERP_FILT==C_FXN
  for (row = 0; row < height; row++)
  {
    for (col = 0; col < width; col++)
    {
      Int sum;

      sum  = src[ col + 0 * cStride] * c[0];
      sum += src[ col + 1 * cStride] * c[1];
      if ( N >= 4 )
      {
        sum += src[ col + 2 * cStride] * c[2];
        sum += src[ col + 3 * cStride] * c[3];
      }
      if ( N >= 6 )
      {
        sum += src[ col + 4 * cStride] * c[4];
        sum += src[ col + 5 * cStride] * c[5];
      }
      if ( N == 8 )
      {
        sum += src[ col + 6 * cStride] * c[6];
        sum += src[ col + 7 * cStride] * c[7];
      }

      Pel val = ( sum + offset ) >> shift;
      if ( isLast )
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
        val = ClipA(val,compID);
#else
        val = ( val < 0 ) ? 0 : val;
        val = ( val > maxVal ) ? maxVal : val;
#endif
      }
      dst[col] = val;
    }

    src += srcStride;
    dst += dstStride;
  }
#endif
}

/**
 * \brief Filter a block of samples (horizontal)
 *
 * \tparam N          Number of taps
 * \param  bitDepth   Bit depth of samples
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  coeff      Pointer to filter taps
 */
#if JVET_D0033_ADAPTIVE_CLIPPING
template<Int N>
Void TComInterpolationFilter::filterHor(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isLast, TFilterCoeff const *coeff, ComponentID compID)
#else
template<Int N>
Void TComInterpolationFilter::filterHor(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isLast, TFilterCoeff const *coeff)
#endif
{
  if ( isLast )
  {
#if JVET_D0033_ADAPTIVE_CLIPPING
      filter<N, false, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
#else
    filter<N, false, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
#endif
  }
  else
  {
#if JVET_D0033_ADAPTIVE_CLIPPING
      filter<N, false, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
#else
    filter<N, false, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
#endif
  }
}

/**
 * \brief Filter a block of samples (vertical)
 *
 * \tparam N          Number of taps
 * \param  bitDepth   Bit depth
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  isFirst    Flag indicating whether it is the first filtering operation
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  coeff      Pointer to filter taps
 */
#if JVET_D0033_ADAPTIVE_CLIPPING
template<Int N>
Void TComInterpolationFilter::filterVer(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, TFilterCoeff const *coeff, ComponentID compID)
{
    if (isFirst && isLast)
    {
        filter<N, true, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
    }
    else if (isFirst && !isLast)
    {
        filter<N, true, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
    }
    else if (!isFirst && isLast)
    {
        filter<N, true, false, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
    }
    else
    {
        filter<N, true, false, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff, compID);
    }
}
#else
template<Int N>
Void TComInterpolationFilter::filterVer(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, TFilterCoeff const *coeff)
{
  if ( isFirst && isLast )
  {
    filter<N, true, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else if ( isFirst && !isLast )
  {
    filter<N, true, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else if ( !isFirst && isLast )
  {
    filter<N, true, false, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else
  {
    filter<N, true, false, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
}
#endif

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
/**
 * \brief Apply FIR filter to a block of samples
 *
 * \tparam N          Number of taps
 * \tparam isVertical Flag indicating filtering along vertical direction
 * \tparam isFirst    Flag indicating whether it is the first filtering operation
 * \tparam isLast     Flag indicating whether it is the last filtering operation
 * \param  bitDepth   Bit depth of samples
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  coeff      Pointer to filter taps
 */
template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
Void TComInterpolationFilter::filterAffine(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, TFilterCoeff const *coeff
                                           #if JVET_D0033_ADAPTIVE_CLIPPING
                                                       , ComponentID compID
                                           #endif
                                           )
{
  Int row, col;

  Pel c[8];
  c[0] = coeff[0];
  c[1] = coeff[1];
  if ( N >= 4 )
  {
    c[2] = coeff[2];
    c[3] = coeff[3];
  }
  if ( N >= 6 )
  {
    c[4] = coeff[4];
    c[5] = coeff[5];
  }
  if ( N == 8 )
  {
    c[6] = coeff[6];
    c[7] = coeff[7];
  }

  Int cStride = ( isVertical ) ? srcStride : 1;
  src -= ( N/2 - 1 ) * cStride;

  Int offset;
  Pel maxVal;
  Int headRoom = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));
  Int shift    = IF_FILTER_PREC_AFFINE;
  // with the current settings (IF_INTERNAL_PREC = 14 and IF_FILTER_PREC = 6), though headroom can be
  // negative for bit depths greater than 14, shift will remain non-negative for bit depths of 8->20
  assert(shift >= 0);

#if JVET_D0033_ADAPTIVE_CLIPPING
    Short minVal=0;
#endif
  if ( isLast )
  {
    shift += (isFirst) ? 0 : headRoom;
    offset = 1 << (shift - 1);
    offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC_AFFINE;
#if JVET_D0033_ADAPTIVE_CLIPPING
        minVal = g_ClipParam.min(compID);
        maxVal = g_ClipParam.max(compID);
#else
    maxVal = (1 << bitDepth) - 1;
#endif
  }
  else
  {
    shift -= (isFirst) ? headRoom : 0;
    offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    maxVal = 0;
  }

#if COM16_C806_SIMD_OPT
  if( bitDepth <= 10 )
  {
    if( N == 8 && !( width & 0x07 ) )
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      Short minVal = 0;
#endif
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 8 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 8 )
        {
          __m128i mmFiltered = simdInterpolateLuma8( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 8 && !( width & 0x03 ) )
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      Short minVal = 0;
#endif
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 8 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateLuma4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 4 && !( width & 0x03 ) )
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      Short minVal = 0;
#endif
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 4 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateChroma4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 2 && !( width & 0x07 ) )
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      Short minVal = 0;
#endif
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[2];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 2 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 8 )
        {
          __m128i mmFiltered = simdInterpolateLuma2P8( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 2 && !( width & 0x03 ) )
    {
#if !JVET_D0033_ADAPTIVE_CLIPPING
      Short minVal = 0;
#endif
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 2 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateLuma2P4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
  }
#endif

  for (row = 0; row < height; row++)
  {
    for (col = 0; col < width; col++)
    {
      Int sum;

      sum  = src[ col + 0 * cStride] * c[0];
      sum += src[ col + 1 * cStride] * c[1];
      if ( N >= 4 )
      {
        sum += src[ col + 2 * cStride] * c[2];
        sum += src[ col + 3 * cStride] * c[3];
      }
      if ( N >= 6 )
      {
        sum += src[ col + 4 * cStride] * c[4];
        sum += src[ col + 5 * cStride] * c[5];
      }
      if ( N == 8 )
      {
        sum += src[ col + 6 * cStride] * c[6];
        sum += src[ col + 7 * cStride] * c[7];
      }

      Pel val = ( sum + offset ) >> shift;
      if ( isLast )
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
                val = ClipA(val,compID);
#else
        val = ( val < 0 ) ? 0 : val;
        val = ( val > maxVal ) ? maxVal : val;
#endif
      }
      dst[col] = val;
    }

    src += srcStride;
    dst += dstStride;
  }
}

/**
 * \brief Filter a block of samples (horizontal)
 *
 * \tparam N          Number of taps
 * \param  bitDepth   Bit depth of samples
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  coeff      Pointer to filter taps
 */
template<Int N>
Void TComInterpolationFilter::filterHorAffine(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isLast, TFilterCoeff const *coeff
                                              #if JVET_D0033_ADAPTIVE_CLIPPING
                                                          ,ComponentID compID
                                              #endif
                                              )
{
  if ( isLast )
  {
        filterAffine<N, false, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                   #if JVET_D0033_ADAPTIVE_CLIPPING
                                               , compID
                                   #endif
                                           );
  }
  else
  {
        filterAffine<N, false, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                    #if JVET_D0033_ADAPTIVE_CLIPPING
                                                , compID
                                    #endif
                                            );
  }
}

/**
 * \brief Filter a block of samples (vertical)
 *
 * \tparam N          Number of taps
 * \param  bitDepth   Bit depth
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  isFirst    Flag indicating whether it is the first filtering operation
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  coeff      Pointer to filter taps
 */
template<Int N>
Void TComInterpolationFilter::filterVerAffine(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, TFilterCoeff const *coeff
                                              #if JVET_D0033_ADAPTIVE_CLIPPING
                                              ,ComponentID compID
                                              #endif
                                              )
{
  if ( isFirst && isLast )
  {
        filterAffine<N, true, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                  #if JVET_D0033_ADAPTIVE_CLIPPING
                                              , compID
                                  #endif
                                          );
  }
  else if ( isFirst && !isLast )
  {
        filterAffine<N, true, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                   #if JVET_D0033_ADAPTIVE_CLIPPING
                                               , compID
                                   #endif
                                           );
  }
  else if ( !isFirst && isLast )
  {
        filterAffine<N, true, false, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                   #if JVET_D0033_ADAPTIVE_CLIPPING
                                               , compID
                                   #endif
                                           );
  }
  else
  {
        filterAffine<N, true, false, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff
                                    #if JVET_D0033_ADAPTIVE_CLIPPING
                                                , compID
                                    #endif
                                            );
  }
}
#endif
#endif

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/**
 * \brief Filter a block of Luma/Chroma samples (horizontal)
 *
 * \param  compID     Chroma component ID
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  frac       Fractional sample offset
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  fmt        Chroma format
 * \param  bitDepth   Bit depth
 */
Void TComInterpolationFilter::filterHor(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isLast, const ChromaFormat fmt, const Int bitDepth 
#if VCEG_AZ07_FRUC_MERGE
  , Int nFilterIdx
#endif
#if SIF_GAUSS
  , Int nFilterGauss
#endif
  )
{
  if ( frac == 0 
#if BIF_FILTER_ZERO_FRAC
    && nFilterGauss == 0
#endif
     )
  {
#if JVET_D0033_ADAPTIVE_CLIPPING
      filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, true, isLast, compID);
#else
    filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, true, isLast );
#endif
  }
  else 
    if (isLuma(compID))
  {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    assert(frac >= 0 && frac < (LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
    assert(frac >= 0 && frac < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS);
#endif
#if VCEG_AZ07_FRUC_MERGE
    if( nFilterIdx == 1 )
      filterHor<NTAPS_LUMA_FRUC>(bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_lumaFilterBilinear[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
      , compID
#endif
      );
    else
#endif
#if SIF_GAUSS
    if ( nFilterGauss == 1 )
    {
      filterHor<NTAPS_LUMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_lumaFilterGauss[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
    }
    else
#endif
    filterHor<NTAPS_LUMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_lumaFilter[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
    , compID
#endif
    );
  }
  else
  {
    const UInt csx = getComponentScaleX(compID, fmt);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    assert(frac >=0 && csx<2 && (frac<<(1-csx)) < (CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
    assert(frac >=0 && csx<2 && (frac<<(1-csx)) < CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS);
#endif
#if SIF_GAUSS
    if (nFilterGauss == 1)
    {
      filterHor<NTAPS_CHROMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_chromaFilterGauss[frac<<(1-csx)]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
    }
    else
#endif
    filterHor<NTAPS_CHROMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_chromaFilter[frac<<(1-csx)]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
  }
}


/**
 * \brief Filter a block of Luma/Chroma samples (vertical)
 *
 * \param  compID     Colour component ID
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  frac       Fractional sample offset
 * \param  isFirst    Flag indicating whether it is the first filtering operation
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  fmt        Chroma format
 * \param  bitDepth   Bit depth
 */
Void TComInterpolationFilter::filterVer(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth 
#if VCEG_AZ07_FRUC_MERGE
  , Int nFilterIdx
#endif
#if SIF_GAUSS
  , Int nFilterGauss
#endif
  )
{
  if ( frac == 0 
#if BIF_FILTER_ZERO_FRAC
    && nFilterGauss == 0
#endif
    )
  {
#if JVET_D0033_ADAPTIVE_CLIPPING
      filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, compID);
#else
    filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast );
#endif
  }
  else 
    if (isLuma(compID))
  {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    assert(frac >= 0 && frac < (LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
    assert(frac >= 0 && frac < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS);
#endif
#if VCEG_AZ07_FRUC_MERGE
    if( nFilterIdx == 1 )
      filterVer<NTAPS_LUMA_FRUC>(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_lumaFilterBilinear[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
      , compID
#endif
      );
    else
#endif
#if SIF_GAUSS
    if( nFilterGauss == 1 )
    {
      filterVer<NTAPS_LUMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_lumaFilterGauss[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
    }
    else
#endif
    filterVer<NTAPS_LUMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_lumaFilter[frac]
#if JVET_D0033_ADAPTIVE_CLIPPING
    , compID
#endif
    );
  }
  else
  {
    const UInt csy = getComponentScaleY(compID, fmt);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    assert(frac >=0 && csy<2 && (frac<<(1-csy)) < (CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
    assert(frac >=0 && csy<2 && (frac<<(1-csy)) < CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS);
#endif
#if SIF_GAUSS
    if (nFilterGauss == 1)
    {
      filterVer<NTAPS_CHROMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_chromaFilterGauss[frac<<(1-csy)]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
    }
    else
#endif
    filterVer<NTAPS_CHROMA>(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_chromaFilter[frac<<(1-csy)]
#if JVET_D0033_ADAPTIVE_CLIPPING
        , compID
#endif
        );
  }
}

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
/**
 * \brief Filter a block of Luma/Chroma samples (horizontal)
 *
 * \param  compID     Color component ID
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  frac       Fractional sample offset
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  fmt        Chroma format
 * \param  bitDepth   Bit depth
 */
Void TComInterpolationFilter::filterHorAffine(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isLast, const ChromaFormat fmt, const Int bitDepth )
{
  if ( frac == 0 )
  {
        filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, true, isLast
           #if JVET_D0033_ADAPTIVE_CLIPPING
                   ,compID
           #endif
                   );
  }
  else if (isLuma(compID))
  {
    assert(frac >= 0 && frac < NFRACS_LUMA_AFFINE);
        filterHorAffine<NTAPS_LUMA>( bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_lumaFilterAffine+frac*NTAPS_LUMA
                             #if JVET_D0033_ADAPTIVE_CLIPPING
                                     ,compID
                             #endif
                                     );
  }
  else
  {
    assert(frac >= 0 && frac < NFRACS_CHROMA_AFFINE);
        filterHorAffine<NTAPS_CHROMA>( bitDepth, src, srcStride, dst, dstStride, width, height, isLast, m_chromaFilterAffine+frac*NTAPS_CHROMA
                               #if JVET_D0033_ADAPTIVE_CLIPPING
                                       ,compID
                               #endif
                                       );
  }
}


/**
 * \brief Filter a block of Luma/Chroma samples (vertical)
 *
 * \param  compID     Color component ID
 * \param  src        Pointer to source samples
 * \param  srcStride  Stride of source samples
 * \param  dst        Pointer to destination samples
 * \param  dstStride  Stride of destination samples
 * \param  width      Width of block
 * \param  height     Height of block
 * \param  frac       Fractional sample offset
 * \param  isFirst    Flag indicating whether it is the first filtering operation
 * \param  isLast     Flag indicating whether it is the last filtering operation
 * \param  fmt        Chroma format
 * \param  bitDepth   Bit depth
 */
Void TComInterpolationFilter::filterVerAffine(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth )
{
  if ( frac == 0 )
  {
        filterCopy(bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast
           #if JVET_D0033_ADAPTIVE_CLIPPING
                   ,compID
           #endif
                   );
  }
  else if (isLuma(compID))
  {
    assert(frac >= 0 && frac < NFRACS_LUMA_AFFINE);
        filterVerAffine<NTAPS_LUMA>( bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_lumaFilterAffine+frac*NTAPS_LUMA
                             #if JVET_D0033_ADAPTIVE_CLIPPING
                                     ,compID
                             #endif
                                     );
  }
  else
  {
    assert(frac >= 0 && frac < NFRACS_CHROMA_AFFINE);
    filterVerAffine<NTAPS_CHROMA>( bitDepth, src, srcStride, dst, dstStride, width, height, isFirst, isLast, m_chromaFilterAffine+frac*NTAPS_CHROMA );
  }
}
#endif
#endif

#if RRC_WIF
Void TComInterpolationFilter::filterDownSample( Pel *src, Int srcStride, Pel *dst, Int dstStride, Pel *tmp, Int tmpStride, Int width, Int height, Int iScale, const Int bitDepth, Bool bBi, Bool bHor, Bool bVer )
{
  Int iTaps  = NTAPS_LUMA_DS;
  Int iLtaps = iTaps / 2;

  Int offset;
  Int shift;
  Int headRoom = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));
  Pel maxVal;

  // horizontal
  shift  = IF_FILTER_PREC_DS - headRoom;
  offset = -IF_INTERNAL_OFFS << shift;
  Pel* _tmp = tmp - iLtaps * tmpStride;
  src      -= iLtaps * srcStride;
  if ( bHor )
  {
    for ( Int h = 0; h < height*iScale+iTaps; h++ )
    {
      for ( Int w = 0; w < width; w++ )
      {
        Int val = 0;
        Int sw = w * iScale;
        for ( Int i = 0; i < iTaps; i++ )
        {
          val += src[sw-iLtaps+i] * g_iFilterCoeff[i];
        }
        _tmp[w] = ( val + offset ) >> shift;
      }
      src  += srcStride;
      _tmp += tmpStride;
    }
  }
  else
  {
    for ( Int h = 0; h < height*iScale+iTaps; h++ )
    {
      for ( Int w = 0; w < width; w++ )
      {
        _tmp[w] = ( src[w*iScale] << headRoom ) - IF_INTERNAL_OFFS;
      }
      src  += srcStride;
      _tmp += tmpStride;
    }
  }

  // vertical
  if ( bBi ) // not last precision
  {
    shift = IF_FILTER_PREC_DS;
    offset = 0;
    maxVal = 0;
  }
  else
  {
    shift = IF_FILTER_PREC_DS + headRoom;
    offset = 1 << (shift - 1);
    offset += IF_INTERNAL_OFFS << IF_FILTER_PREC_DS;
    maxVal = (1<<bitDepth) - 1;
  }

  _tmp = tmp - iLtaps * tmpStride;
  if ( bVer )
  {
    for ( Int w = 0; w < width; w++ )
    {
      for ( Int h = 0; h < height; h++ )
      {
        Int val = 0;
        Int sh = h * iScale;
        for ( Int i = 0; i < iTaps; i++ )
        {
          val += _tmp[(sh+i)*tmpStride] * g_iFilterCoeff[i];
        }
        val = ( val + offset ) >> shift;
        if ( !bBi )
        {
          val = ( val < 0 ) ? 0 : val;
          val = ( val > maxVal ) ? maxVal : val;
        }
        dst[h*dstStride] = val;
      }
      dst++;
      _tmp++;
    }
  }
  else if ( bBi )
  {
    for ( Int w = 0; w < width; w++ )
    {
      for ( Int h = 0; h < height; h++ )
      {
        Int sh = h * iScale + iLtaps;
        dst[h*dstStride] = _tmp[sh*tmpStride];
      }
      dst++;
      _tmp++;
    }
  }
  else
  {
    shift = headRoom;
    for ( Int w = 0; w < width; w++ )
    {
      for ( Int h = 0; h < height; h++ )
      {
        Int sh = h*iScale + iLtaps;
        //Pel val = _tmp[sh*tmpStride];
        Pel val = rightShift_round( (_tmp[sh*tmpStride]+IF_INTERNAL_OFFS), shift );
        val = ( val < 0 ) ? 0 : val;
        val = ( val > maxVal ) ? maxVal : val;
        dst[h*dstStride] = val;
      }
      dst++;
      _tmp++;
    }
  }

  return;
}

Void TComInterpolationFilter::filterHorDs(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Int scale, Bool isLast, const ChromaFormat fmt, const Int bitDepth )
{
  if ( frac == 0 )
  {
    filterCopyDs(bitDepth, src, srcStride, dst, dstStride, width, height, scale, true, isLast );
  }
  else if (isLuma(compID))
  {
    filterHorDs<NTAPS_LUMA_DS>( bitDepth, src, srcStride, dst, dstStride, width, height, scale, isLast, m_lumaFilterAffine+frac*NTAPS_LUMA );
  }
  else
  {
    filterHorDs<NTAPS_CHROMA_DS>( bitDepth, src, srcStride, dst, dstStride, width, height, scale, isLast, m_chromaFilterAffine+frac*NTAPS_CHROMA );
  }
}

Void TComInterpolationFilter::filterVerDs(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Int scale, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth )
{
  if ( frac == 0 )
  {
    filterCopyDs(bitDepth, src, srcStride, dst, dstStride, width, height, scale, isFirst, isLast );
  }
  else if (isLuma(compID))
  {
    filterVerDs<NTAPS_LUMA_DS>( bitDepth, src, srcStride, dst, dstStride, width, height, scale, isFirst, isLast, m_lumaFilterAffine+frac*NTAPS_LUMA );
  }
  else
  {
    filterVerDs<NTAPS_CHROMA_DS>( bitDepth, src, srcStride, dst, dstStride, width, height, scale, isFirst, isLast, m_chromaFilterAffine+frac*NTAPS_CHROMA );
  }
}


Void TComInterpolationFilter::filterCopyDs(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int scale, Bool isFirst, Bool isLast)
{
  Int row, col;

  if ( isFirst == isLast )
  {
    for (row = 0; row < height; row++)
    {
      for (col = 0; col < width; col++)
      {
        dst[col] = src[col];
      }

      src += srcStride;
      dst += dstStride;
    }
  }
  else if ( isFirst )
  {
    const Int shift = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));

    for (row = 0; row < height; row++)
    {
      for (col = 0; col < width; col++)
      {
        Pel val = leftShift_round(src[col], shift);
        dst[col] = val - (Pel)IF_INTERNAL_OFFS;
      }

      src += srcStride;
      dst += dstStride;
    }
  }
  else
  {
    const Int shift = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));

    Pel maxVal = (1 << bitDepth) - 1;
    Pel minVal = 0;
    for (row = 0; row < height; row++)
    {
      for (col = 0; col < width; col++)
      {
        Pel val = src[ col ];
        val = rightShift_round((val + IF_INTERNAL_OFFS), shift);
        if (val < minVal)
        {
          val = minVal;
        }
        if (val > maxVal)
        {
          val = maxVal;
        }
        dst[col] = val;
      }

      src += srcStride;
      dst += dstStride;
    }
  }
}

template<Int N>
Void TComInterpolationFilter::filterHorDs(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int scale, Bool isLast, TFilterCoeff const *coeff)
{
  if ( isLast )
  {
    filter<N, false, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else
  {
    filter<N, false, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
}

template<Int N>
Void TComInterpolationFilter::filterVerDs(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int scale, Bool isFirst, Bool isLast, TFilterCoeff const *coeff)
{
  if ( isFirst && isLast )
  {
    filter<N, true, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else if ( isFirst && !isLast )
  {
    filter<N, true, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else if ( !isFirst && isLast )
  {
    filter<N, true, false, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
  else
  {
    filter<N, true, false, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
  }
}

template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
Void TComInterpolationFilter::filterDs(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int scale, TFilterCoeff const *coeff)
{
  Int row, col;

  Pel c[8];
  c[0] = coeff[0];
  c[1] = coeff[1];
  if ( N >= 4 )
  {
    c[2] = coeff[2];
    c[3] = coeff[3];
  }
  if ( N >= 6 )
  {
    c[4] = coeff[4];
    c[5] = coeff[5];
  }
  if ( N == 8 )
  {
    c[6] = coeff[6];
    c[7] = coeff[7];
  }

  Int cStride = ( isVertical ) ? srcStride : 1;
  src -= ( N/2 - 1 ) * cStride;

  Int offset;
  Pel maxVal;
  Int headRoom = std::max<Int>(2, (IF_INTERNAL_PREC - bitDepth));
  Int shift    = IF_FILTER_PREC_AFFINE;
  // with the current settings (IF_INTERNAL_PREC = 14 and IF_FILTER_PREC = 6), though headroom can be
  // negative for bit depths greater than 14, shift will remain non-negative for bit depths of 8->20
  assert(shift >= 0);

  if ( isLast )
  {
    shift += (isFirst) ? 0 : headRoom;
    offset = 1 << (shift - 1);
    offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC_AFFINE;
    maxVal = (1 << bitDepth) - 1;
  }
  else
  {
    shift -= (isFirst) ? headRoom : 0;
    offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    maxVal = 0;
  }

#if COM16_C806_SIMD_OPT
  if( bitDepth <= 10 )
  {
    if( N == 8 && !( width & 0x07 ) )
    {
      Short minVal = 0;
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 8 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 8 )
        {
          __m128i mmFiltered = simdInterpolateLuma8( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 8 && !( width & 0x03 ) )
    {
      Short minVal = 0;
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 8 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateLuma4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 4 && !( width & 0x03 ) )
    {
      Short minVal = 0;
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 4 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateChroma4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 2 && !( width & 0x07 ) )
    {
      Short minVal = 0;
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[2];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 2 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 8 )
        {
          __m128i mmFiltered = simdInterpolateLuma2P8( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storeu_si128( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
    else if( N == 2 && !( width & 0x03 ) )
    {
      Short minVal = 0;
      __m128i mmOffset = _mm_set1_epi32( offset );
      __m128i mmCoeff[8];
      __m128i mmMin = _mm_set1_epi16( minVal );
      __m128i mmMax = _mm_set1_epi16( maxVal );
      for( Int n = 0 ; n < 2 ; n++ )
        mmCoeff[n] = _mm_set1_epi16( c[n] );
      for( row = 0 ; row < height ; row++ )
      {
        for( col = 0 ; col < width ; col += 4 )
        {
          __m128i mmFiltered = simdInterpolateLuma2P4( src + col , cStride , mmCoeff , mmOffset , shift );
          if( isLast )
          {
            mmFiltered = simdClip3( mmMin , mmMax , mmFiltered );
          }
          _mm_storel_epi64( ( __m128i * )( dst + col ) , mmFiltered );
        }
        src += srcStride;
        dst += dstStride;
      }
      return;
    }
  }
#endif

  for (row = 0; row < height; row++)
  {
    for (col = 0; col < width; col++)
    {
      Int sum;

      sum  = src[ col + 0 * cStride] * c[0];
      sum += src[ col + 1 * cStride] * c[1];
      if ( N >= 4 )
      {
        sum += src[ col + 2 * cStride] * c[2];
        sum += src[ col + 3 * cStride] * c[3];
      }
      if ( N >= 6 )
      {
        sum += src[ col + 4 * cStride] * c[4];
        sum += src[ col + 5 * cStride] * c[5];
      }
      if ( N == 8 )
      {
        sum += src[ col + 6 * cStride] * c[6];
        sum += src[ col + 7 * cStride] * c[7];
      }

      Pel val = ( sum + offset ) >> shift;
      if ( isLast )
      {
        val = ( val < 0 ) ? 0 : val;
        val = ( val > maxVal ) ? maxVal : val;
      }
      dst[col] = val;
    }

    src += srcStride;
    dst += dstStride;
  }
}
#endif

//! \}
