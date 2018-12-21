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
 * \brief Declaration of TComInterpolationFilter class
 */

#ifndef __TCOMINTERPOLATIONFILTER__
#define __TCOMINTERPOLATIONFILTER__

#include "CommonDef.h"

//! \ingroup TLibCommon
//! \{

#define NTAPS_LUMA        8 ///< Number of taps for luma
#define NTAPS_CHROMA      4 ///< Number of taps for chroma
#if VCEG_AZ07_FRUC_MERGE
#define NTAPS_LUMA_FRUC   2
#endif
#define IF_INTERNAL_PREC 14 ///< Number of bits for internal precision
#define IF_FILTER_PREC    6 ///< Log2 of sum of filter taps
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1)) ///< Offset used internally

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#define IF_FILTER_PREC_AFFINE    8  ///< Log2 of sum of affine filter taps
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
#define NFRACS_LUMA_AFFINE       64 ///< Number of fraction positions for luma affine MCP
#define NFRACS_CHROMA_AFFINE     64 ///< Number of fraction positions for chroma affine MCP
#endif
#endif

#if RRC_WIF
#define IF_FILTER_PREC_DS        6
#define NTAPS_LUMA_DS            7
#define NTAPS_CHROMA_DS          7
#endif

/**
 * \brief Interpolation filter class
 */
class TComInterpolationFilter
{
#if JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC
  static const TFilterCoeff m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_LUMA];     ///< Luma filter taps
  static const TFilterCoeff m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_CHROMA]; ///< Chroma filter taps
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
  static const Short m_lumaFilter[8][NTAPS_LUMA];     ///< Luma filter taps
  static const Short m_chromaFilter[16][NTAPS_CHROMA]; ///< Chroma filter taps
#else
  static const TFilterCoeff m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA];     ///< Luma filter taps
  static const TFilterCoeff m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps
#endif
#endif
#if VCEG_AZ07_FRUC_MERGE
#if JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC
  static const Short m_lumaFilterBilinear[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_LUMA_FRUC];     ///< Luma filter taps
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
  static const Short m_lumaFilterBilinear[8][NTAPS_LUMA_FRUC];     ///< Luma filter taps
#else
  static const Short m_lumaFilterBilinear[4][NTAPS_LUMA_FRUC];     ///< Luma filter taps
#endif
#endif
#endif

#if SIF_GAUSS
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  static const Short m_lumaFilterGauss  [LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE  ][NTAPS_LUMA  ];  ///< Luma filter taps
  static const Short m_chromaFilterGauss[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE][NTAPS_CHROMA];  ///< Chroma filter taps
#else
  static const Short m_lumaFilterGauss  [LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS  ][NTAPS_LUMA  ];  ///< Luma filter taps
  static const Short m_chromaFilterGauss[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA];  ///< Chroma filter taps
#endif
#endif
#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
  static const Short m_lumaFilterAffine[(NFRACS_LUMA_AFFINE)*NTAPS_LUMA];
  static const Short m_chromaFilterAffine[(NFRACS_CHROMA_AFFINE)*NTAPS_CHROMA];
#endif
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
  static Void filterCopy(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, ComponentID compID);
  template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
  static Void filter(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, TFilterCoeff const *coeff, ComponentID compID);
#else
  static Void filterCopy(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast);

  template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
  static Void filter(Int bitDepth, Pel const *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, TFilterCoeff const *coeff);
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
  template<Int N>
  static Void filterHor(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height,               Bool isLast, TFilterCoeff const *coeff, ComponentID compID);
  template<Int N>
  static Void filterVer(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, TFilterCoeff const *coeff, ComponentID compID);
#else
  template<Int N>
  static Void filterHor(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height,               Bool isLast, TFilterCoeff const *coeff);
  template<Int N>
  static Void filterVer(Int bitDepth, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, TFilterCoeff const *coeff);
#endif

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
#if JVET_D0033_ADAPTIVE_CLIPPING
  template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
  static Void filterAffine(Int bitDepth, Pel const *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Short const *coeff,ComponentID compID);

  template<Int N>
  static Void filterHorAffine(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height,               Bool isLast, Short const *coeff, ComponentID compID);
  template<Int N>
  static Void filterVerAffine(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, Short const *coeff, ComponentID compID);
#else
  template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
  static Void filterAffine(Int bitDepth, Pel const *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Short const *coeff);

  template<Int N>
  static Void filterHorAffine(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height,               Bool isLast, Short const *coeff);
  template<Int N>
  static Void filterVerAffine(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Bool isFirst, Bool isLast, Short const *coeff);
#endif
#endif
#endif

#if RRC_WIF
  static Void filterCopyDs(Int bitDepth, const Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int scale, Bool isFirst, Bool isLast);

  template<Int N, Bool isVertical, Bool isFirst, Bool isLast>
  static Void filterDs(Int bitDepth, Pel const *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Int scale, Short const *coeff);

  template<Int N>
  static Void filterHorDs(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Int scale,                Bool isLast, Short const *coeff);
  template<Int N>
  static Void filterVerDs(Int bitDepth, Pel *src, Int srcStride, Short *dst, Int dstStride, Int width, Int height, Int scale,  Bool isFirst, Bool isLast, Short const *coeff);
#endif

public:
  TComInterpolationFilter() {}
  ~TComInterpolationFilter() {}

  Void filterHor(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac,               Bool isLast, const ChromaFormat fmt, const Int bitDepth 
#if VCEG_AZ07_FRUC_MERGE
    , Int nFilterIdx = 0
#endif
#if SIF_GAUSS
    , Int nFilterGauss = 0
#endif
    );
  Void filterVer(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth 
#if VCEG_AZ07_FRUC_MERGE
    , Int nFilterIdx = 0
#endif
#if SIF_GAUSS
    , Int nFilterGauss = 0
#endif
    );

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
#if !JVET_C0025_AFFINE_FILTER_SIMPLIFICATION
  Void filterHorAffine(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac,               Bool isLast, const ChromaFormat fmt, const Int bitDepth );
  Void filterVerAffine(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth );
#endif


#if RRC_WIF
  Void filterDownSample( Pel *src, Int srcStride, Pel *dst, Int dstStride, Pel *tmp, Int tmpStride, Int width, Int height, Int iScale, const Int bitDepth, Bool bBi, Bool bHor, Bool bVer );

  Void filterHorDs(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Int scale,               Bool isLast, const ChromaFormat fmt, const Int bitDepth );
  Void filterVerDs(const ComponentID compID, Pel *src, Int srcStride, Pel *dst, Int dstStride, Int width, Int height, Int frac, Int scale, Bool isFirst, Bool isLast, const ChromaFormat fmt, const Int bitDepth );
#endif

#endif
};

//! \}

#endif
