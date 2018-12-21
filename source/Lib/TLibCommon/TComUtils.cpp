/*!
*******************************************************************************
*                                JEM_ITT Encoder
*                      ITTIAM SYSTEMS PVT LTD, BANGALORE
*                               COPYRIGHT(C) 2017
*
*  This program is proprietary to Ittiam Systems Pvt. Ltd. and is protected
*  under Indian Copyright Act as an unpublished work. Its use and disclosure is
*  limited by the terms and conditions of a license agreement. It may not be
*  copied or otherwise reproduced or disclosed to persons outside the licensee's
*  organization except in accordance with the terms and conditions of such an
*  agreement. All copies and reproductions shall be the property of Ittiam
*******************************************************************************
*/

/** \file     TComUtils.cpp
\brief    Common utility funcions
*/

#include <assert.h>

#include "CommonDef.h"
#include "TComInterpolationFilter.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <stdlib.h>
#include <cstring>

#if ITH_C3_TWODIMBUFCOPIER==SIMD_FXN
void twoDimBufCopier
(
UChar *pDst,
const Int iDstStride,
const UChar *pSrc,
const Int iSrcStride,
const Int iBlkWd,
const Int iBlkHt
)
{
  Int iRows, iCols;
  const UChar *inp_copy;
  UChar *dst_copy;
  Int iRemWd = iBlkWd;

  __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;
  __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;

  inp_copy = pSrc;
  dst_copy = pDst;

  if (iRemWd >= 16)
  {
    Int iRemHt = iBlkHt % 8;
    Int iCount = 0;
    iRemWd &= 0xf;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)	/*loop over height of 8*/
    {
      for (iCols = iBlkWd; iCols > 15; iCols -= 16)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride)); // i = 0
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride)); // i = 1
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride)); // i = 2
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride)); // i = 3
        src4_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 4 * iSrcStride)); // i = 0
        src5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 5 * iSrcStride)); // i = 1
        src6_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 6 * iSrcStride)); // i = 2
        src7_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 7 * iSrcStride)); // i = 3

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 4 * iDstStride), src4_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 5 * iDstStride), src5_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 6 * iDstStride), src6_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 7 * iDstStride), src7_16x8b);

        iCount += 16;
      }

      iCount = 0;
      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);

    }

    if (iRemHt > 3) /*loop over height of 4*/
    {
      iRemHt = iRemHt % 4;

      for (iCols = iBlkWd; iCols > 15; iCols -= 16)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride)); // i = 0
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride)); // i = 1
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride)); // i = 2
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride)); // i = 3

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);

        iCount += 16;
      }
      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)	/*loop over height of 1*/
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {
        iCount = 0;
        for (iCols = iBlkWd; iCols > 15; iCols -= 16)
        {
          src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount)); // i = 3
          _mm_storeu_si128((__m128i*)(dst_copy + iCount), src0_16x8b);
          iCount += 16;
        }
        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }

    }
  }

  inp_copy = pSrc + ((iBlkWd / 16) * 16);
  dst_copy = pDst + ((iBlkWd / 16) * 16);

  if (iRemWd > 7)
  {
    Int iRemHt = iBlkHt % 8;
    Int iCount = 0;

    for (Int j = 0; j < (iBlkHt - 7); j += 8)
    {
      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 0 * iSrcStride)); // i = 0
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 1 * iSrcStride)); // i = 1
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 2 * iSrcStride)); // i = 2
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 3 * iSrcStride)); // i = 3
      src4_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 4 * iSrcStride)); // i = 0
      src5_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 5 * iSrcStride)); // i = 1
      src6_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 6 * iSrcStride)); // i = 2
      src7_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + iCount + 7 * iSrcStride)); // i = 3

      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 4 * iDstStride), src4_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 5 * iDstStride), src5_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 6 * iDstStride), src6_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 7 * iDstStride), src7_16x8b);

      iCount = 0;
      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);
    }

    if (iRemHt > 3)
    {
      iRemHt = iRemHt % 4;

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride)); // i = 0
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride)); // i = 1
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride)); // i = 2
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride)); // i = 3

      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);

      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {
        iCount = 0;

        src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy)); // i = 3
        _mm_storel_epi64((__m128i*)(dst_copy), src0_16x8b);

        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }
    }
  }

  iRemWd &= 0x7;
  inp_copy = pSrc + ((iBlkWd / 8) * 8);
  dst_copy = pDst + ((iBlkWd / 8) * 8);

  if (iRemWd > 3)
  {
    Int iRemHt = iBlkHt % 8;
    iRemWd &= 0x3;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)
    {

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride)); // i = 0
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride)); // i = 1
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride)); // i = 2
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride)); // i = 3
      src4_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 4 * iSrcStride)); // i = 0
      src5_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 5 * iSrcStride)); // i = 1
      src6_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 6 * iSrcStride)); // i = 2
      src7_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 7 * iSrcStride)); // i = 3

      *(Int *)(&dst_copy[0 * iDstStride]) = _mm_cvtsi128_si32(src0_16x8b);
      *(Int *)(&dst_copy[1 * iDstStride]) = _mm_cvtsi128_si32(src1_16x8b);
      *(Int *)(&dst_copy[2 * iDstStride]) = _mm_cvtsi128_si32(src2_16x8b);
      *(Int *)(&dst_copy[3 * iDstStride]) = _mm_cvtsi128_si32(src3_16x8b);
      *(Int *)(&dst_copy[4 * iDstStride]) = _mm_cvtsi128_si32(src4_16x8b);
      *(Int *)(&dst_copy[5 * iDstStride]) = _mm_cvtsi128_si32(src5_16x8b);
      *(Int *)(&dst_copy[6 * iDstStride]) = _mm_cvtsi128_si32(src6_16x8b);
      *(Int *)(&dst_copy[7 * iDstStride]) = _mm_cvtsi128_si32(src7_16x8b);

      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);
    }

    if (iRemHt > 3)
    {
      iRemHt = iRemHt % 4;

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride)); // i = 0
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride)); // i = 1
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride)); // i = 2
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride)); // i = 3

      *(Int *)(&dst_copy[0 * iDstStride]) = _mm_cvtsi128_si32(src0_16x8b);
      *(Int *)(&dst_copy[1 * iDstStride]) = _mm_cvtsi128_si32(src1_16x8b);
      *(Int *)(&dst_copy[2 * iDstStride]) = _mm_cvtsi128_si32(src2_16x8b);
      *(Int *)(&dst_copy[3 * iDstStride]) = _mm_cvtsi128_si32(src3_16x8b);

      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {

        src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy));
        *(Int *)(&dst_copy[0]) = _mm_cvtsi128_si32(src0_16x8b);

        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }
    }
  }

  inp_copy = pSrc + ((iBlkWd / 4) * 4);
  dst_copy = pDst + ((iBlkWd / 4) * 4);

  if (iRemWd)
  {
    for (Int i = 0; i < iBlkHt; i++)
    {
      memcpy(dst_copy, inp_copy, iRemWd*sizeof(UChar));

      dst_copy += iDstStride;
      inp_copy += iSrcStride;
    }
  }
}
#elif ITH_C3_TWODIMBUFCOPIER==C_FXN
void twoDimBufCopier
(
UChar *pDst,
const Int iDstStride,
const UChar *pSrc,
const Int iSrcStride,
const Int iBlkWd,
const Int iBlkHt
)
{
  for (Int i = 0; i < iBlkHt; i++)
  {
    ::memcpy
      (
      pDst + i * iDstStride,
      pSrc + i * iSrcStride,
      iBlkWd
      );
  }
}
#endif

#if ITH_C3_RESIDUECOMPUTER==SIMD_FXN
Void residueComputer
(
Pel *pSrc,
const Int iSrcStride,
Pel *pPred,
const Int iPredStride,
Pel *pResidue,
const Int iResidueStride,
const Int iWidth,
const Int iHeight
)
{
  __m128i mmSrc[4];
  __m128i mmPred[4];
  __m128i mmRes[4];

  Int iNumRowsLeft = iHeight;
  Pel *const pSrcOrig = pSrc;
  Pel *const pPredOrig = pPred;
  Pel *const pResidueOrig = pResidue;

  if (!(iWidth % 8))
  {
    /* Guard is 'y < (iHeight - 2)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 2); y += 4, pSrc += 4 * iSrcStride, pPred += 4 * iPredStride, pResidue += 4 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pSrc += 8, pPred += 8, pResidue += 8)
      {
        mmSrc[0] = _mm_loadu_si128((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadu_si128((const __m128i *) (pSrc + 1 * iSrcStride));
        mmSrc[2] = _mm_loadu_si128((const __m128i *) (pSrc + 2 * iSrcStride));
        mmSrc[3] = _mm_loadu_si128((const __m128i *) (pSrc + 3 * iSrcStride));

        mmPred[0] = _mm_loadu_si128((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadu_si128((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadu_si128((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadu_si128((const __m128i *) (pPred + 3 * iPredStride));

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[1] = _mm_sub_epi16(mmSrc[1], mmPred[1]);
        mmRes[2] = _mm_sub_epi16(mmSrc[2], mmPred[2]);
        mmRes[3] = _mm_sub_epi16(mmSrc[3], mmPred[3]);

        _mm_storeu_si128((__m128i *) (pResidue + 0 * iResidueStride), mmRes[0]);
        _mm_storeu_si128((__m128i *) (pResidue + 1 * iResidueStride), mmRes[1]);
        _mm_storeu_si128((__m128i *) (pResidue + 2 * iResidueStride), mmRes[2]);
        _mm_storeu_si128((__m128i *) (pResidue + 3 * iResidueStride), mmRes[3]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc = pSrcOrig + (iHeight - iNumRowsLeft) * iSrcStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc += 2 * iSrcStride, pPred += 2 * iPredStride, pResidue += 2 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pSrc += 8, pPred += 8, pResidue += 8)
      {
        mmSrc[0] = _mm_loadu_si128((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadu_si128((const __m128i *) (pSrc + 1 * iSrcStride));

        mmPred[0] = _mm_loadu_si128((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadu_si128((const __m128i *) (pPred + 1 * iPredStride));

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[1] = _mm_sub_epi16(mmSrc[1], mmPred[1]);

        _mm_storeu_si128((__m128i *) (pResidue + 0 * iResidueStride), mmRes[0]);
        _mm_storeu_si128((__m128i *) (pResidue + 1 * iResidueStride), mmRes[1]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;
    }
  }
  else if (!(iWidth % 4))
  {
    /* Guard is 'y < (iHeight - 2)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 2); y += 4, pSrc += 4 * iSrcStride, pPred += 4 * iPredStride, pResidue += 4 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pSrc += 4, pPred += 4, pResidue += 4)
      {
        mmSrc[0] = _mm_loadl_epi64((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadl_epi64((const __m128i *) (pSrc + 1 * iSrcStride));
        mmSrc[2] = _mm_loadl_epi64((const __m128i *) (pSrc + 2 * iSrcStride));
        mmSrc[3] = _mm_loadl_epi64((const __m128i *) (pSrc + 3 * iSrcStride));

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[1] = _mm_sub_epi16(mmSrc[1], mmPred[1]);
        mmRes[2] = _mm_sub_epi16(mmSrc[2], mmPred[2]);
        mmRes[3] = _mm_sub_epi16(mmSrc[3], mmPred[3]);

        _mm_storel_epi64((__m128i *) (pResidue + 0 * iResidueStride), mmRes[0]);
        _mm_storel_epi64((__m128i *) (pResidue + 1 * iResidueStride), mmRes[1]);
        _mm_storel_epi64((__m128i *) (pResidue + 2 * iResidueStride), mmRes[2]);
        _mm_storel_epi64((__m128i *) (pResidue + 3 * iResidueStride), mmRes[3]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc = pSrcOrig + (iHeight - iNumRowsLeft) * iSrcStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc += 2 * iSrcStride, pPred += 2 * iPredStride, pResidue += 2 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pSrc += 4, pPred += 4, pResidue += 4)
      {
        mmSrc[0] = _mm_loadl_epi64((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadl_epi64((const __m128i *) (pSrc + 1 * iSrcStride));

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[1] = _mm_sub_epi16(mmSrc[1], mmPred[1]);

        _mm_storel_epi64((__m128i *) (pResidue + 0 * iResidueStride), mmRes[0]);
        _mm_storel_epi64((__m128i *) (pResidue + 1 * iResidueStride), mmRes[1]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;
    }
  }
  else
  {
    /* Guard is 'y < (iHeight - 2)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 2); y += 4, pSrc += 4 * iSrcStride, pPred += 4 * iPredStride, pResidue += 4 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pSrc += 2, pPred += 2, pResidue += 2)
      {
        mmSrc[0] = _mm_loadl_epi64((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadl_epi64((const __m128i *) (pSrc + 1 * iSrcStride));
        mmSrc[2] = _mm_loadl_epi64((const __m128i *) (pSrc + 2 * iSrcStride));
        mmSrc[3] = _mm_loadl_epi64((const __m128i *) (pSrc + 3 * iSrcStride));

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));

        mmSrc[0] = _mm_unpacklo_epi32(mmSrc[0], mmSrc[1]);
        mmSrc[1] = _mm_unpacklo_epi32(mmSrc[2], mmSrc[3]);
        mmPred[0] = _mm_unpacklo_epi32(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi32(mmPred[2], mmPred[3]);

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[2] = _mm_sub_epi16(mmSrc[1], mmPred[1]);
        mmRes[1] = _mm_srli_si128(mmRes[0], 4);
        mmRes[3] = _mm_srli_si128(mmRes[2], 4);

        ((Int *)(pResidue + 0 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[0]);
        ((Int *)(pResidue + 1 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[1]);
        ((Int *)(pResidue + 2 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[2]);
        ((Int *)(pResidue + 3 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[3]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc = pSrcOrig + (iHeight - iNumRowsLeft) * iSrcStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc += 2 * iSrcStride, pPred += 2 * iPredStride, pResidue += 2 * iResidueStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pSrc += 2, pPred += 2, pResidue += 2)
      {
        mmSrc[0] = _mm_loadl_epi64((const __m128i *) (pSrc + 0 * iSrcStride));
        mmSrc[1] = _mm_loadl_epi64((const __m128i *) (pSrc + 1 * iSrcStride));

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));

        mmSrc[0] = _mm_unpacklo_epi32(mmSrc[0], mmSrc[1]);
        mmPred[0] = _mm_unpacklo_epi32(mmPred[0], mmPred[1]);

        mmRes[0] = _mm_sub_epi16(mmSrc[0], mmPred[0]);
        mmRes[1] = _mm_srli_si128(mmRes[0], 4);

        ((Int *)(pResidue + 0 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[0]);
        ((Int *)(pResidue + 1 * iResidueStride))[0] = _mm_cvtsi128_si32(mmRes[1]);
      }

      pSrc -= iWidth;
      pPred -= iWidth;
      pResidue -= iWidth;
    }
  }
}
#elif ITH_C3_RESIDUECOMPUTER==C_FXN
Void residueComputer
(
Pel *pSrc,
const Int iSrcStride,
Pel *pPred,
const Int iPredStride,
Pel *pResidue,
const Int iResidueStride,
const Int iWidth,
const Int iHeight
)
{
  for (Int j = 0; j < iHeight; j++)
  {
    for (Int i = 0; i < iWidth; i++)
    {
      pResidue[i + j * iResidueStride] = pSrc[i] - pPred[i];
    }

    pSrc += iSrcStride;
    pPred += iPredStride;
  }
}
#endif

Void LongLong2DoubleConverter
(
Int64 *pi64Src,
Double *pdDst,
Int iSize
)
{
  for (Int i = 0; i < iSize; i++)
  {
    pdDst[i] = (Double)pi64Src[i];
  }
}

#if ITH_C3_RECONCOMPUTER==C_FXN
Void reconComputer
(
Pel *pResidue,
Int iResidueStride,
Pel *pPred,
Int iPredStride,
Pel *pRecon,
Int iReconStride,
Int iMin,
Int iMax,
Int iWidth,
Int iHeight
)
{
  Int iIntermediary;

  for (Int j = 0; j < iHeight; j++)
  {
    for (Int i = 0; i < iWidth; i++)
    {
      iIntermediary = pResidue[i] + pPred[i];
      pRecon[i + j * iReconStride] = Clip3<Int>(iMin, iMax, iIntermediary);
    }

    pResidue += iResidueStride;
    pPred += iPredStride;
  }
}
#elif ITH_C3_RECONCOMPUTER==SIMD_FXN
Void reconComputer
(
Pel *pResidue,
Int iResidueStride,
Pel *pPred,
Int iPredStride,
Pel *pRecon,
Int iReconStride,
Int iMin,
Int iMax,
Int iWidth,
Int iHeight
)
{
  __m128i mmResidue[8];
  __m128i mmPred[8];
  __m128i mmRecon[8];
  __m128i mmIntermediary[4];

  const __m128i mmMin = _mm_set1_epi16(iMin);
  const __m128i mmMax = _mm_set1_epi16(iMax);
  Int iNumRowsLeft = iHeight;
  Pel *const pResidueOrig = pResidue;
  Pel *const pPredOrig = pPred;
  Pel *const pReconOrig = pRecon;

  if (!(iWidth % 8))
  {
    /* Guard is 'y < (iHeight - 2)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 2); y += 4, pResidue += 4 * iResidueStride, pPred += 4 * iPredStride, pRecon += 4 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pResidue += 8, pPred += 8, pRecon += 8)
      {
        mmResidue[0] = _mm_loadu_si128((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadu_si128((const __m128i *) (pResidue + 1 * iResidueStride));
        mmResidue[2] = _mm_loadu_si128((const __m128i *) (pResidue + 2 * iResidueStride));
        mmResidue[3] = _mm_loadu_si128((const __m128i *) (pResidue + 3 * iResidueStride));

        mmPred[0] = _mm_loadu_si128((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadu_si128((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadu_si128((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadu_si128((const __m128i *) (pPred + 3 * iPredStride));

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);
        mmRecon[1] = _mm_add_epi16(mmResidue[1], mmPred[1]);
        mmRecon[2] = _mm_add_epi16(mmResidue[2], mmPred[2]);
        mmRecon[3] = _mm_add_epi16(mmResidue[3], mmPred[3]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmRecon[1], mmMax);
        mmIntermediary[2] = _mm_cmpgt_epi16(mmRecon[2], mmMax);
        mmIntermediary[3] = _mm_cmpgt_epi16(mmRecon[3], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMax, mmIntermediary[1]);
        mmRecon[2] = _mm_blendv_epi8(mmRecon[2], mmMax, mmIntermediary[2]);
        mmRecon[3] = _mm_blendv_epi8(mmRecon[3], mmMax, mmIntermediary[3]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmRecon[1], mmMin);
        mmIntermediary[2] = _mm_cmplt_epi16(mmRecon[2], mmMin);
        mmIntermediary[3] = _mm_cmplt_epi16(mmRecon[3], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMin, mmIntermediary[1]);
        mmRecon[2] = _mm_blendv_epi8(mmRecon[2], mmMin, mmIntermediary[2]);
        mmRecon[3] = _mm_blendv_epi8(mmRecon[3], mmMin, mmIntermediary[3]);

        _mm_storeu_si128((__m128i *) (pRecon + 0 * iReconStride), mmRecon[0]);
        _mm_storeu_si128((__m128i *) (pRecon + 1 * iReconStride), mmRecon[1]);
        _mm_storeu_si128((__m128i *) (pRecon + 2 * iReconStride), mmRecon[2]);
        _mm_storeu_si128((__m128i *) (pRecon + 3 * iReconStride), mmRecon[3]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;

      iNumRowsLeft -= 4;
    }

    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pRecon = pReconOrig + (iHeight - iNumRowsLeft) * iReconStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pResidue += 2 * iResidueStride, pPred += 2 * iPredStride, pRecon += 2 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pResidue += 8, pPred += 8, pRecon += 8)
      {
        mmResidue[0] = _mm_loadu_si128((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadu_si128((const __m128i *) (pResidue + 1 * iResidueStride));

        mmPred[0] = _mm_loadu_si128((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadu_si128((const __m128i *) (pPred + 1 * iPredStride));

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);
        mmRecon[1] = _mm_add_epi16(mmResidue[1], mmPred[1]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmRecon[1], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmRecon[1], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMin, mmIntermediary[1]);

        _mm_storeu_si128((__m128i *) (pRecon + 0 * iReconStride), mmRecon[0]);
        _mm_storeu_si128((__m128i *) (pRecon + 1 * iReconStride), mmRecon[1]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;
    }
  }
  else if (!(iWidth % 4))
  {
    /* Guard is 'y < (iHeight - 6)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 6); y += 8, pResidue += 8 * iResidueStride, pPred += 8 * iPredStride, pRecon += 8 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pResidue += 4, pPred += 4, pRecon += 4)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));
        mmResidue[2] = _mm_loadl_epi64((const __m128i *) (pResidue + 2 * iResidueStride));
        mmResidue[3] = _mm_loadl_epi64((const __m128i *) (pResidue + 3 * iResidueStride));
        mmResidue[4] = _mm_loadl_epi64((const __m128i *) (pResidue + 4 * iResidueStride));
        mmResidue[5] = _mm_loadl_epi64((const __m128i *) (pResidue + 5 * iResidueStride));
        mmResidue[6] = _mm_loadl_epi64((const __m128i *) (pResidue + 6 * iResidueStride));
        mmResidue[7] = _mm_loadl_epi64((const __m128i *) (pResidue + 7 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi64(mmResidue[0], mmResidue[1]);
        mmResidue[1] = _mm_unpacklo_epi64(mmResidue[2], mmResidue[3]);
        mmResidue[2] = _mm_unpacklo_epi64(mmResidue[4], mmResidue[5]);
        mmResidue[3] = _mm_unpacklo_epi64(mmResidue[6], mmResidue[7]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));
        mmPred[4] = _mm_loadl_epi64((const __m128i *) (pPred + 4 * iPredStride));
        mmPred[5] = _mm_loadl_epi64((const __m128i *) (pPred + 5 * iPredStride));
        mmPred[6] = _mm_loadl_epi64((const __m128i *) (pPred + 6 * iPredStride));
        mmPred[7] = _mm_loadl_epi64((const __m128i *) (pPred + 7 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi64(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi64(mmPred[2], mmPred[3]);
        mmPred[2] = _mm_unpacklo_epi64(mmPred[4], mmPred[5]);
        mmPred[3] = _mm_unpacklo_epi64(mmPred[6], mmPred[7]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);
        mmRecon[1] = _mm_add_epi16(mmResidue[1], mmPred[1]);
        mmRecon[2] = _mm_add_epi16(mmResidue[2], mmPred[2]);
        mmRecon[3] = _mm_add_epi16(mmResidue[3], mmPred[3]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmRecon[1], mmMax);
        mmIntermediary[2] = _mm_cmpgt_epi16(mmRecon[2], mmMax);
        mmIntermediary[3] = _mm_cmpgt_epi16(mmRecon[3], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMax, mmIntermediary[1]);
        mmRecon[2] = _mm_blendv_epi8(mmRecon[2], mmMax, mmIntermediary[2]);
        mmRecon[3] = _mm_blendv_epi8(mmRecon[3], mmMax, mmIntermediary[3]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmRecon[1], mmMin);
        mmIntermediary[2] = _mm_cmplt_epi16(mmRecon[2], mmMin);
        mmIntermediary[3] = _mm_cmplt_epi16(mmRecon[3], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMin, mmIntermediary[1]);
        mmRecon[2] = _mm_blendv_epi8(mmRecon[2], mmMin, mmIntermediary[2]);
        mmRecon[3] = _mm_blendv_epi8(mmRecon[3], mmMin, mmIntermediary[3]);

        mmRecon[4] = _mm_srli_si128(mmRecon[0], 8);
        mmRecon[5] = _mm_srli_si128(mmRecon[1], 8);
        mmRecon[6] = _mm_srli_si128(mmRecon[2], 8);
        mmRecon[7] = _mm_srli_si128(mmRecon[3], 8);

        _mm_storel_epi64((__m128i *) (pRecon + 0 * iReconStride), mmRecon[0]);
        _mm_storel_epi64((__m128i *) (pRecon + 1 * iReconStride), mmRecon[4]);
        _mm_storel_epi64((__m128i *) (pRecon + 2 * iReconStride), mmRecon[1]);
        _mm_storel_epi64((__m128i *) (pRecon + 3 * iReconStride), mmRecon[5]);
        _mm_storel_epi64((__m128i *) (pRecon + 4 * iReconStride), mmRecon[2]);
        _mm_storel_epi64((__m128i *) (pRecon + 5 * iReconStride), mmRecon[6]);
        _mm_storel_epi64((__m128i *) (pRecon + 6 * iReconStride), mmRecon[3]);
        _mm_storel_epi64((__m128i *) (pRecon + 7 * iReconStride), mmRecon[7]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;

      iNumRowsLeft -= 8;
    }

    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pRecon = pReconOrig + (iHeight - iNumRowsLeft) * iReconStride;

    for (Int y = iHeight - iNumRowsLeft; y < (iHeight - 2); y += 4, pResidue += 4 * iResidueStride, pPred += 4 * iPredStride, pRecon += 4 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pResidue += 4, pPred += 4, pRecon += 4)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));
        mmResidue[2] = _mm_loadl_epi64((const __m128i *) (pResidue + 2 * iResidueStride));
        mmResidue[3] = _mm_loadl_epi64((const __m128i *) (pResidue + 3 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi64(mmResidue[0], mmResidue[1]);
        mmResidue[1] = _mm_unpacklo_epi64(mmResidue[2], mmResidue[3]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi64(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi64(mmPred[2], mmPred[3]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);
        mmRecon[1] = _mm_add_epi16(mmResidue[1], mmPred[1]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmRecon[1], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmRecon[1], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);
        mmRecon[2] = _mm_blendv_epi8(mmRecon[1], mmMin, mmIntermediary[1]);

        mmRecon[1] = _mm_srli_si128(mmRecon[0], 8);
        mmRecon[3] = _mm_srli_si128(mmRecon[2], 8);

        _mm_storel_epi64((__m128i *) (pRecon + 0 * iReconStride), mmRecon[0]);
        _mm_storel_epi64((__m128i *) (pRecon + 2 * iReconStride), mmRecon[2]);
        _mm_storel_epi64((__m128i *) (pRecon + 1 * iReconStride), mmRecon[1]);
        _mm_storel_epi64((__m128i *) (pRecon + 3 * iReconStride), mmRecon[3]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;

      iNumRowsLeft -= 4;
    }

    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pRecon = pReconOrig + (iHeight - iNumRowsLeft) * iReconStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pResidue += 2 * iResidueStride, pPred += 2 * iPredStride, pRecon += 2 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pResidue += 4, pPred += 4, pRecon += 4)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi64(mmResidue[0], mmResidue[1]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi64(mmPred[0], mmPred[1]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);

        mmRecon[1] = _mm_srli_si128(mmRecon[0], 8);

        _mm_storel_epi64((__m128i *) (pRecon + 0 * iReconStride), mmRecon[0]);
        _mm_storel_epi64((__m128i *) (pRecon + 1 * iReconStride), mmRecon[1]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;
    }
  }
  else
  {
    /* Guard is 'y < (iHeight - 6)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 6); y += 8, pResidue += 8 * iResidueStride, pPred += 8 * iPredStride, pRecon += 8 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pResidue += 2, pPred += 2, pRecon += 2)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));
        mmResidue[2] = _mm_loadl_epi64((const __m128i *) (pResidue + 2 * iResidueStride));
        mmResidue[3] = _mm_loadl_epi64((const __m128i *) (pResidue + 3 * iResidueStride));
        mmResidue[4] = _mm_loadl_epi64((const __m128i *) (pResidue + 4 * iResidueStride));
        mmResidue[5] = _mm_loadl_epi64((const __m128i *) (pResidue + 5 * iResidueStride));
        mmResidue[6] = _mm_loadl_epi64((const __m128i *) (pResidue + 6 * iResidueStride));
        mmResidue[7] = _mm_loadl_epi64((const __m128i *) (pResidue + 7 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi32(mmResidue[0], mmResidue[1]);
        mmResidue[1] = _mm_unpacklo_epi32(mmResidue[2], mmResidue[3]);
        mmResidue[2] = _mm_unpacklo_epi32(mmResidue[4], mmResidue[5]);
        mmResidue[3] = _mm_unpacklo_epi32(mmResidue[6], mmResidue[7]);

        mmResidue[0] = _mm_unpacklo_epi64(mmResidue[0], mmResidue[1]);
        mmResidue[1] = _mm_unpacklo_epi64(mmResidue[2], mmResidue[3]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));
        mmPred[4] = _mm_loadl_epi64((const __m128i *) (pPred + 4 * iPredStride));
        mmPred[5] = _mm_loadl_epi64((const __m128i *) (pPred + 5 * iPredStride));
        mmPred[6] = _mm_loadl_epi64((const __m128i *) (pPred + 6 * iPredStride));
        mmPred[7] = _mm_loadl_epi64((const __m128i *) (pPred + 7 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi32(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi32(mmPred[2], mmPred[3]);
        mmPred[2] = _mm_unpacklo_epi32(mmPred[4], mmPred[5]);
        mmPred[3] = _mm_unpacklo_epi32(mmPred[6], mmPred[7]);

        mmPred[0] = _mm_unpacklo_epi64(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi64(mmPred[2], mmPred[3]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);
        mmRecon[1] = _mm_add_epi16(mmResidue[1], mmPred[1]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmRecon[1], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);
        mmRecon[1] = _mm_blendv_epi8(mmRecon[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmRecon[1], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);
        mmRecon[4] = _mm_blendv_epi8(mmRecon[1], mmMin, mmIntermediary[1]);

        mmRecon[1] = _mm_srli_si128(mmRecon[0], 4);
        mmRecon[2] = _mm_srli_si128(mmRecon[0], 8);
        mmRecon[3] = _mm_srli_si128(mmRecon[0], 12);
        mmRecon[5] = _mm_srli_si128(mmRecon[4], 4);
        mmRecon[6] = _mm_srli_si128(mmRecon[4], 8);
        mmRecon[7] = _mm_srli_si128(mmRecon[4], 12);

        ((Int *)(pRecon + 0 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[0]);
        ((Int *)(pRecon + 4 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[4]);
        ((Int *)(pRecon + 1 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[1]);
        ((Int *)(pRecon + 2 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[2]);
        ((Int *)(pRecon + 3 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[3]);
        ((Int *)(pRecon + 5 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[5]);
        ((Int *)(pRecon + 6 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[6]);
        ((Int *)(pRecon + 7 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[7]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;

      iNumRowsLeft -= 8;
    }

    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pRecon = pReconOrig + (iHeight - iNumRowsLeft) * iReconStride;

    for (Int y = iHeight - iNumRowsLeft; y < (iHeight - 2); y += 4, pResidue += 4 * iResidueStride, pPred += 4 * iPredStride, pRecon += 4 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pResidue += 2, pPred += 2, pRecon += 2)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));
        mmResidue[2] = _mm_loadl_epi64((const __m128i *) (pResidue + 2 * iResidueStride));
        mmResidue[3] = _mm_loadl_epi64((const __m128i *) (pResidue + 3 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi32(mmResidue[0], mmResidue[1]);
        mmResidue[1] = _mm_unpacklo_epi32(mmResidue[2], mmResidue[3]);

        mmResidue[0] = _mm_unpacklo_epi64(mmResidue[0], mmResidue[1]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));
        mmPred[2] = _mm_loadl_epi64((const __m128i *) (pPred + 2 * iPredStride));
        mmPred[3] = _mm_loadl_epi64((const __m128i *) (pPred + 3 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi32(mmPred[0], mmPred[1]);
        mmPred[1] = _mm_unpacklo_epi32(mmPred[2], mmPred[3]);

        mmPred[0] = _mm_unpacklo_epi64(mmPred[0], mmPred[1]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);

        mmRecon[1] = _mm_srli_si128(mmRecon[0], 4);
        mmRecon[2] = _mm_srli_si128(mmRecon[0], 8);
        mmRecon[3] = _mm_srli_si128(mmRecon[0], 12);

        ((Int *)(pRecon + 0 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[0]);
        ((Int *)(pRecon + 1 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[1]);
        ((Int *)(pRecon + 2 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[2]);
        ((Int *)(pRecon + 3 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[3]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;

      iNumRowsLeft -= 4;
    }

    pResidue = pResidueOrig + (iHeight - iNumRowsLeft) * iResidueStride;
    pPred = pPredOrig + (iHeight - iNumRowsLeft) * iPredStride;
    pRecon = pReconOrig + (iHeight - iNumRowsLeft) * iReconStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pResidue += 2 * iResidueStride, pPred += 2 * iPredStride, pRecon += 2 * iReconStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pResidue += 2, pPred += 2, pRecon += 2)
      {
        mmResidue[0] = _mm_loadl_epi64((const __m128i *) (pResidue + 0 * iResidueStride));
        mmResidue[1] = _mm_loadl_epi64((const __m128i *) (pResidue + 1 * iResidueStride));

        mmResidue[0] = _mm_unpacklo_epi32(mmResidue[0], mmResidue[1]);

        mmPred[0] = _mm_loadl_epi64((const __m128i *) (pPred + 0 * iPredStride));
        mmPred[1] = _mm_loadl_epi64((const __m128i *) (pPred + 1 * iPredStride));

        mmPred[0] = _mm_unpacklo_epi32(mmPred[0], mmPred[1]);

        mmRecon[0] = _mm_add_epi16(mmResidue[0], mmPred[0]);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmRecon[0], mmMax);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmRecon[0], mmMin);

        mmRecon[0] = _mm_blendv_epi8(mmRecon[0], mmMin, mmIntermediary[0]);

        mmRecon[1] = _mm_srli_si128(mmRecon[0], 4);

        ((Int *)(pRecon + 0 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[0]);
        ((Int *)(pRecon + 1 * iReconStride))[0] = _mm_cvtsi128_si32(mmRecon[1]);
      }

      pResidue -= iWidth;
      pPred -= iWidth;
      pRecon -= iWidth;
    }
  }
}
#endif

#if 1//!COM16_C806_SIMD_OPT
Void averagerWithRounding
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
      iIntermediary = pSrc0[i] + pSrc1[i] + iOffset;
      pDst[i + j * iDstStride] = Clip3<Int>(iMin, iMax, iIntermediary >> iShift);
    }

    pSrc0 += iSrc0Stride;
    pSrc1 += iSrc1Stride;
  }
}
#else
Void averagerWithRounding
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
  __m128i mmSrc0[8];
  __m128i mmSrc1[8];
  __m128i mmDst[8];
  __m128i mmIntermediary[4];

  const __m128i mmMin = _mm_set1_epi16(iMin);
  const __m128i mmMax = _mm_set1_epi16(iMax);
  const __m128i mmOffset = _mm_set1_epi16(iOffset);
  Int iNumRowsLeft = iHeight;
  Pel *pSrc0Orig = pSrc0;
  Pel *pSrc1Orig = pSrc1;
  Pel *pDstOrig = pDst;

  if (!(iWidth % 8))
  {
    /* Guard is 'y < (iHeight - 2)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 2); y += 4, pSrc0 += 4 * iSrc0Stride, pSrc1 += 4 * iSrc1Stride, pDst += 4 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pSrc0 += 8, pSrc1 += 8, pDst += 8)
      {
        mmSrc0[0] = _mm_loadu_si128((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadu_si128((const __m128i *) (pSrc0 + 1 * iSrc0Stride));
        mmSrc0[2] = _mm_loadu_si128((const __m128i *) (pSrc0 + 2 * iSrc0Stride));
        mmSrc0[3] = _mm_loadu_si128((const __m128i *) (pSrc0 + 3 * iSrc0Stride));

        mmSrc1[0] = _mm_loadu_si128((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadu_si128((const __m128i *) (pSrc1 + 1 * iSrc1Stride));
        mmSrc1[2] = _mm_loadu_si128((const __m128i *) (pSrc1 + 2 * iSrc1Stride));
        mmSrc1[3] = _mm_loadu_si128((const __m128i *) (pSrc1 + 3 * iSrc1Stride));

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);
        mmDst[1] = _mm_add_epi16(mmSrc0[1], mmSrc1[1]);
        mmDst[2] = _mm_add_epi16(mmSrc0[2], mmSrc1[2]);
        mmDst[3] = _mm_add_epi16(mmSrc0[3], mmSrc1[3]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);
        mmDst[1] = _mm_add_epi16(mmDst[1], mmOffset);
        mmDst[2] = _mm_add_epi16(mmDst[2], mmOffset);
        mmDst[3] = _mm_add_epi16(mmDst[3], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);
        mmDst[1] = _mm_srai_epi16(mmDst[1], iShift);
        mmDst[2] = _mm_srai_epi16(mmDst[2], iShift);
        mmDst[3] = _mm_srai_epi16(mmDst[3], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmDst[1], mmMax);
        mmIntermediary[2] = _mm_cmpgt_epi16(mmDst[2], mmMax);
        mmIntermediary[3] = _mm_cmpgt_epi16(mmDst[3], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMax, mmIntermediary[1]);
        mmDst[2] = _mm_blendv_epi8(mmDst[2], mmMax, mmIntermediary[2]);
        mmDst[3] = _mm_blendv_epi8(mmDst[3], mmMax, mmIntermediary[3]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmDst[1], mmMin);
        mmIntermediary[2] = _mm_cmplt_epi16(mmDst[2], mmMin);
        mmIntermediary[3] = _mm_cmplt_epi16(mmDst[3], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMin, mmIntermediary[1]);
        mmDst[2] = _mm_blendv_epi8(mmDst[2], mmMin, mmIntermediary[2]);
        mmDst[3] = _mm_blendv_epi8(mmDst[3], mmMin, mmIntermediary[3]);

        _mm_storeu_si128((__m128i *) (pDst + 0 * iDstStride), mmDst[0]);
        _mm_storeu_si128((__m128i *) (pDst + 1 * iDstStride), mmDst[1]);
        _mm_storeu_si128((__m128i *) (pDst + 2 * iDstStride), mmDst[2]);
        _mm_storeu_si128((__m128i *) (pDst + 3 * iDstStride), mmDst[3]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc0 = pSrc0Orig + (iHeight - iNumRowsLeft) * iSrc0Stride;
    pSrc1 = pSrc1Orig + (iHeight - iNumRowsLeft) * iSrc1Stride;
    pDst = pDstOrig + (iHeight - iNumRowsLeft) * iDstStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc0 += 2 * iSrc0Stride, pSrc1 += 2 * iSrc1Stride, pDst += 2 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 8, pSrc0 += 8, pSrc1 += 8, pDst += 8)
      {
        mmSrc0[0] = _mm_loadu_si128((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadu_si128((const __m128i *) (pSrc0 + 1 * iSrc0Stride));

        mmSrc1[0] = _mm_loadu_si128((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadu_si128((const __m128i *) (pSrc1 + 1 * iSrc1Stride));

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);
        mmDst[1] = _mm_add_epi16(mmSrc0[1], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);
        mmDst[1] = _mm_add_epi16(mmDst[1], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);
        mmDst[1] = _mm_srai_epi16(mmDst[1], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmDst[1], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmDst[1], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMin, mmIntermediary[1]);

        _mm_storeu_si128((__m128i *) (pDst + 0 * iDstStride), mmDst[0]);
        _mm_storeu_si128((__m128i *) (pDst + 1 * iDstStride), mmDst[1]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;
    }
  }
  else if (!(iWidth % 4))
  {
    /* Guard is 'y < (iHeight - 6)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 6); y += 8, pSrc0 += 8 * iSrc0Stride, pSrc1 += 8 * iSrc1Stride, pDst += 8 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pSrc0 += 4, pSrc1 += 4, pDst += 4)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));
        mmSrc0[2] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 2 * iSrc0Stride));
        mmSrc0[3] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 3 * iSrc0Stride));
        mmSrc0[4] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 4 * iSrc0Stride));
        mmSrc0[5] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 5 * iSrc0Stride));
        mmSrc0[6] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 6 * iSrc0Stride));
        mmSrc0[7] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 7 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi64(mmSrc0[0], mmSrc0[1]);
        mmSrc0[1] = _mm_unpacklo_epi64(mmSrc0[2], mmSrc0[3]);
        mmSrc0[2] = _mm_unpacklo_epi64(mmSrc0[4], mmSrc0[5]);
        mmSrc0[3] = _mm_unpacklo_epi64(mmSrc0[6], mmSrc0[7]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));
        mmSrc1[2] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 2 * iSrc1Stride));
        mmSrc1[3] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 3 * iSrc1Stride));
        mmSrc1[4] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 4 * iSrc1Stride));
        mmSrc1[5] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 5 * iSrc1Stride));
        mmSrc1[6] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 6 * iSrc1Stride));
        mmSrc1[7] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 7 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi64(mmSrc1[0], mmSrc1[1]);
        mmSrc1[1] = _mm_unpacklo_epi64(mmSrc1[2], mmSrc1[3]);
        mmSrc1[2] = _mm_unpacklo_epi64(mmSrc1[4], mmSrc1[5]);
        mmSrc1[3] = _mm_unpacklo_epi64(mmSrc1[6], mmSrc1[7]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);
        mmDst[1] = _mm_add_epi16(mmSrc0[1], mmSrc1[1]);
        mmDst[2] = _mm_add_epi16(mmSrc0[2], mmSrc1[2]);
        mmDst[3] = _mm_add_epi16(mmSrc0[3], mmSrc1[3]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);
        mmDst[1] = _mm_add_epi16(mmDst[1], mmOffset);
        mmDst[2] = _mm_add_epi16(mmDst[2], mmOffset);
        mmDst[3] = _mm_add_epi16(mmDst[3], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);
        mmDst[1] = _mm_srai_epi16(mmDst[1], iShift);
        mmDst[2] = _mm_srai_epi16(mmDst[2], iShift);
        mmDst[3] = _mm_srai_epi16(mmDst[3], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmDst[1], mmMax);
        mmIntermediary[2] = _mm_cmpgt_epi16(mmDst[2], mmMax);
        mmIntermediary[3] = _mm_cmpgt_epi16(mmDst[3], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMax, mmIntermediary[1]);
        mmDst[2] = _mm_blendv_epi8(mmDst[2], mmMax, mmIntermediary[2]);
        mmDst[3] = _mm_blendv_epi8(mmDst[3], mmMax, mmIntermediary[3]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmDst[1], mmMin);
        mmIntermediary[2] = _mm_cmplt_epi16(mmDst[2], mmMin);
        mmIntermediary[3] = _mm_cmplt_epi16(mmDst[3], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMin, mmIntermediary[1]);
        mmDst[2] = _mm_blendv_epi8(mmDst[2], mmMin, mmIntermediary[2]);
        mmDst[3] = _mm_blendv_epi8(mmDst[3], mmMin, mmIntermediary[3]);

        mmDst[4] = _mm_srli_si128(mmDst[0], 8);
        mmDst[5] = _mm_srli_si128(mmDst[1], 8);
        mmDst[6] = _mm_srli_si128(mmDst[2], 8);
        mmDst[7] = _mm_srli_si128(mmDst[3], 8);

        _mm_storel_epi64((__m128i *) (pDst + 0 * iDstStride), mmDst[0]);
        _mm_storel_epi64((__m128i *) (pDst + 1 * iDstStride), mmDst[4]);
        _mm_storel_epi64((__m128i *) (pDst + 2 * iDstStride), mmDst[1]);
        _mm_storel_epi64((__m128i *) (pDst + 3 * iDstStride), mmDst[5]);
        _mm_storel_epi64((__m128i *) (pDst + 4 * iDstStride), mmDst[2]);
        _mm_storel_epi64((__m128i *) (pDst + 5 * iDstStride), mmDst[6]);
        _mm_storel_epi64((__m128i *) (pDst + 6 * iDstStride), mmDst[3]);
        _mm_storel_epi64((__m128i *) (pDst + 7 * iDstStride), mmDst[7]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;

      iNumRowsLeft -= 8;
    }

    pSrc0 = pSrc0Orig + (iHeight - iNumRowsLeft) * iSrc0Stride;
    pSrc1 = pSrc1Orig + (iHeight - iNumRowsLeft) * iSrc1Stride;
    pDst = pDstOrig + (iHeight - iNumRowsLeft) * iDstStride;

    for (Int y = iHeight - iNumRowsLeft; y < (iHeight - 2); y += 4, pSrc0 += 4 * iSrc0Stride, pSrc1 += 4 * iSrc1Stride, pDst += 4 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pSrc0 += 4, pSrc1 += 4, pDst += 4)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));
        mmSrc0[2] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 2 * iSrc0Stride));
        mmSrc0[3] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 3 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi64(mmSrc0[0], mmSrc0[1]);
        mmSrc0[1] = _mm_unpacklo_epi64(mmSrc0[2], mmSrc0[3]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));
        mmSrc1[2] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 2 * iSrc1Stride));
        mmSrc1[3] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 3 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi64(mmSrc1[0], mmSrc1[1]);
        mmSrc1[1] = _mm_unpacklo_epi64(mmSrc1[2], mmSrc1[3]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);
        mmDst[1] = _mm_add_epi16(mmSrc0[1], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);
        mmDst[1] = _mm_add_epi16(mmDst[1], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);
        mmDst[1] = _mm_srai_epi16(mmDst[1], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmDst[1], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmDst[1], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);
        mmDst[2] = _mm_blendv_epi8(mmDst[1], mmMin, mmIntermediary[1]);

        mmDst[1] = _mm_srli_si128(mmDst[0], 8);
        mmDst[3] = _mm_srli_si128(mmDst[2], 8);

        _mm_storel_epi64((__m128i *) (pDst + 0 * iDstStride), mmDst[0]);
        _mm_storel_epi64((__m128i *) (pDst + 2 * iDstStride), mmDst[2]);
        _mm_storel_epi64((__m128i *) (pDst + 1 * iDstStride), mmDst[1]);
        _mm_storel_epi64((__m128i *) (pDst + 3 * iDstStride), mmDst[3]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc0 = pSrc0Orig + (iHeight - iNumRowsLeft) * iSrc0Stride;
    pSrc1 = pSrc1Orig + (iHeight - iNumRowsLeft) * iSrc1Stride;
    pDst = pDstOrig + (iHeight - iNumRowsLeft) * iDstStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc0 += 2 * iSrc0Stride, pSrc1 += 2 * iSrc1Stride, pDst += 2 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 4, pSrc0 += 4, pSrc1 += 4, pDst += 4)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi64(mmSrc0[0], mmSrc0[1]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi64(mmSrc1[0], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);

        mmDst[1] = _mm_srli_si128(mmDst[0], 8);

        _mm_storel_epi64((__m128i *) (pDst + 0 * iDstStride), mmDst[0]);
        _mm_storel_epi64((__m128i *) (pDst + 1 * iDstStride), mmDst[1]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;
    }
  }
  else
  {
    /* Guard is 'y < (iHeight - 6)' to ensure that only blocks with iNumRowsLeft >= 4 enter */
    for (Int y = 0; y < (iHeight - 6); y += 8, pSrc0 += 8 * iSrc0Stride, pSrc1 += 8 * iSrc1Stride, pDst += 8 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pSrc0 += 2, pSrc1 += 2, pDst += 2)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));
        mmSrc0[2] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 2 * iSrc0Stride));
        mmSrc0[3] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 3 * iSrc0Stride));
        mmSrc0[4] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 4 * iSrc0Stride));
        mmSrc0[5] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 5 * iSrc0Stride));
        mmSrc0[6] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 6 * iSrc0Stride));
        mmSrc0[7] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 7 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi32(mmSrc0[0], mmSrc0[1]);
        mmSrc0[1] = _mm_unpacklo_epi32(mmSrc0[2], mmSrc0[3]);
        mmSrc0[2] = _mm_unpacklo_epi32(mmSrc0[4], mmSrc0[5]);
        mmSrc0[3] = _mm_unpacklo_epi32(mmSrc0[6], mmSrc0[7]);

        mmSrc0[0] = _mm_unpacklo_epi64(mmSrc0[0], mmSrc0[1]);
        mmSrc0[1] = _mm_unpacklo_epi64(mmSrc0[2], mmSrc0[3]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));
        mmSrc1[2] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 2 * iSrc1Stride));
        mmSrc1[3] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 3 * iSrc1Stride));
        mmSrc1[4] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 4 * iSrc1Stride));
        mmSrc1[5] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 5 * iSrc1Stride));
        mmSrc1[6] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 6 * iSrc1Stride));
        mmSrc1[7] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 7 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi32(mmSrc1[0], mmSrc1[1]);
        mmSrc1[1] = _mm_unpacklo_epi32(mmSrc1[2], mmSrc1[3]);
        mmSrc1[2] = _mm_unpacklo_epi32(mmSrc1[4], mmSrc1[5]);
        mmSrc1[3] = _mm_unpacklo_epi32(mmSrc1[6], mmSrc1[7]);

        mmSrc1[0] = _mm_unpacklo_epi64(mmSrc1[0], mmSrc1[1]);
        mmSrc1[1] = _mm_unpacklo_epi64(mmSrc1[2], mmSrc1[3]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);
        mmDst[1] = _mm_add_epi16(mmSrc0[1], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);
        mmDst[1] = _mm_add_epi16(mmDst[1], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);
        mmDst[1] = _mm_srai_epi16(mmDst[1], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);
        mmIntermediary[1] = _mm_cmpgt_epi16(mmDst[1], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);
        mmDst[1] = _mm_blendv_epi8(mmDst[1], mmMax, mmIntermediary[1]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);
        mmIntermediary[1] = _mm_cmplt_epi16(mmDst[1], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);
        mmDst[4] = _mm_blendv_epi8(mmDst[1], mmMin, mmIntermediary[1]);

        mmDst[1] = _mm_srli_si128(mmDst[0], 4);
        mmDst[2] = _mm_srli_si128(mmDst[0], 8);
        mmDst[3] = _mm_srli_si128(mmDst[0], 12);
        mmDst[5] = _mm_srli_si128(mmDst[4], 4);
        mmDst[6] = _mm_srli_si128(mmDst[4], 8);
        mmDst[7] = _mm_srli_si128(mmDst[4], 12);

        ((Int *)(pDst + 0 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[0]);
        ((Int *)(pDst + 4 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[4]);
        ((Int *)(pDst + 1 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[1]);
        ((Int *)(pDst + 2 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[2]);
        ((Int *)(pDst + 3 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[3]);
        ((Int *)(pDst + 5 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[5]);
        ((Int *)(pDst + 6 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[6]);
        ((Int *)(pDst + 7 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[7]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;

      iNumRowsLeft -= 8;
    }

    pSrc0 = pSrc0Orig + (iHeight - iNumRowsLeft) * iSrc0Stride;
    pSrc1 = pSrc1Orig + (iHeight - iNumRowsLeft) * iSrc1Stride;
    pDst = pDstOrig + (iHeight - iNumRowsLeft) * iDstStride;

    for (Int y = iHeight - iNumRowsLeft; y < (iHeight - 2); y += 4, pSrc0 += 4 * iSrc0Stride, pSrc1 += 4 * iSrc1Stride, pDst += 4 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pSrc0 += 2, pSrc1 += 2, pDst += 2)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));
        mmSrc0[2] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 2 * iSrc0Stride));
        mmSrc0[3] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 3 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi32(mmSrc0[0], mmSrc0[1]);
        mmSrc0[1] = _mm_unpacklo_epi32(mmSrc0[2], mmSrc0[3]);

        mmSrc0[0] = _mm_unpacklo_epi64(mmSrc0[0], mmSrc0[1]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));
        mmSrc1[2] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 2 * iSrc1Stride));
        mmSrc1[3] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 3 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi32(mmSrc1[0], mmSrc1[1]);
        mmSrc1[1] = _mm_unpacklo_epi32(mmSrc1[2], mmSrc1[3]);

        mmSrc1[0] = _mm_unpacklo_epi64(mmSrc1[0], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);

        mmDst[1] = _mm_srli_si128(mmDst[0], 4);
        mmDst[2] = _mm_srli_si128(mmDst[0], 8);
        mmDst[3] = _mm_srli_si128(mmDst[0], 12);

        ((Int *)(pDst + 0 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[0]);
        ((Int *)(pDst + 1 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[1]);
        ((Int *)(pDst + 2 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[2]);
        ((Int *)(pDst + 3 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[3]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;

      iNumRowsLeft -= 4;
    }

    pSrc0 = pSrc0Orig + (iHeight - iNumRowsLeft) * iSrc0Stride;
    pSrc1 = pSrc1Orig + (iHeight - iNumRowsLeft) * iSrc1Stride;
    pDst = pDstOrig + (iHeight - iNumRowsLeft) * iDstStride;

    for (Int y = iHeight - iNumRowsLeft; y < iHeight; y += 2, pSrc0 += 2 * iSrc0Stride, pSrc1 += 2 * iSrc1Stride, pDst += 2 * iDstStride)
    {
      for (Int x = 0; x < iWidth; x += 2, pSrc0 += 2, pSrc1 += 2, pDst += 2)
      {
        mmSrc0[0] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 0 * iSrc0Stride));
        mmSrc0[1] = _mm_loadl_epi64((const __m128i *) (pSrc0 + 1 * iSrc0Stride));

        mmSrc0[0] = _mm_unpacklo_epi32(mmSrc0[0], mmSrc0[1]);

        mmSrc1[0] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 0 * iSrc1Stride));
        mmSrc1[1] = _mm_loadl_epi64((const __m128i *) (pSrc1 + 1 * iSrc1Stride));

        mmSrc1[0] = _mm_unpacklo_epi32(mmSrc1[0], mmSrc1[1]);

        mmDst[0] = _mm_add_epi16(mmSrc0[0], mmSrc1[0]);

        mmDst[0] = _mm_add_epi16(mmDst[0], mmOffset);

        mmDst[0] = _mm_srai_epi16(mmDst[0], iShift);

        mmIntermediary[0] = _mm_cmpgt_epi16(mmDst[0], mmMax);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMax, mmIntermediary[0]);

        mmIntermediary[0] = _mm_cmplt_epi16(mmDst[0], mmMin);

        mmDst[0] = _mm_blendv_epi8(mmDst[0], mmMin, mmIntermediary[0]);

        mmDst[1] = _mm_srli_si128(mmDst[0], 4);

        ((Int *)(pDst + 0 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[0]);
        ((Int *)(pDst + 1 * iDstStride))[0] = _mm_cvtsi128_si32(mmDst[1]);
      }

      pSrc0 -= iWidth;
      pSrc1 -= iWidth;
      pDst -= iWidth;
    }
  }
}
#endif

#if ITH_C3_SCALEADDCLIP==C_FXN
void scaleAddClip
(
Pel *piSrc,
Int iSrcStride,
Pel *piDst,
Int iDstStride,
UInt cxWidth,
UInt cxHeight,
Int a,
const Int iShift,
Int b,
const ComponentID compID,
const Int bitDepth
)
{
  int i, j;

  for (i = 0; i < cxHeight; i++)
  {
    for (j = 0; j < cxWidth; j++)
    {
#if JVET_D0033_ADAPTIVE_CLIPPING
      piDst[j] = ClipA(((a*piSrc[j]) >> iShift) + b, compID);
#endif
    }
    piDst += iDstStride;
    piSrc += iSrcStride;
  }
}
#elif ITH_C3_SCALEADDCLIP==SIMD_FXN
void scaleAddClip
(
Pel *piSrc,
Int iSrcStride,
Pel *piDst,
Int iDstStride,
UInt cxWidth,
UInt cxHeight,
Int a,
const Int iShift,
Int b,
const ComponentID compID,
const Int bitDepth
)
{
  if (bitDepth <= 10)
  {
    Int row, col, width, height;
    Pel m, M;

    /* all 128 bit registers are named with a suffix mxnb, where m is the */
    /* number of n bits packed in the register                            */
    __m128i src_temp0_4x32b, src_temp1_4x32b, src_temp2_4x32b, src_temp3_4x32b;
    __m128i wgt0_8x16b, off0_8x16b, off0_8x16b_shift, madd_const_8x16b;
    __m128i clip_min_8x16b, clip_max_8x16b;
    __m128i zero_8x16b, one_8x16b;

    width = cxWidth;
    height = cxHeight;

    assert(width % 2 == 0);
    assert(height % 2 == 0);

    m = g_ClipParam.min(compID);
    M = g_ClipParam.max(compID);

    // seting values in register
    wgt0_8x16b = _mm_set1_epi16(a);
    off0_8x16b = _mm_set1_epi16(b);
    clip_min_8x16b = _mm_set1_epi16(m);
    clip_max_8x16b = _mm_set1_epi16(M);
    zero_8x16b = _mm_setzero_si128();
    one_8x16b = _mm_set1_epi16(1);

    off0_8x16b_shift = _mm_slli_epi16(off0_8x16b, iShift);
    madd_const_8x16b = _mm_unpacklo_epi16(wgt0_8x16b, off0_8x16b_shift);

    if ((0 == (width & 7)) && (0 == (height & 3))) /* width multiple of 8 case */
    {
      __m128i src_temp4_4x32b, src_temp5_4x32b, src_temp6_4x32b, src_temp7_4x32b;
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        for (col = 0; col < width; col += 8)
        {
          /* for row =0 ,1,2,3*/
          /* row = 0 */
          src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
          /* row = 1 */
          src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));
          /* row = 2 */
          src_temp2_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride));
          /* row = 3 */
          src_temp3_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride));

          /* row = 0 */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 4));
          /* row = 1 */
          src_temp5_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride + 4));
          /* row = 2 */
          src_temp6_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride + 4));
          /* row = 3 */
          src_temp7_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride + 4));

          /* First 4 pixels */
          src_temp0_4x32b = _mm_unpacklo_epi16(src_temp0_4x32b, one_8x16b);
          src_temp1_4x32b = _mm_unpacklo_epi16(src_temp1_4x32b, one_8x16b);
          src_temp2_4x32b = _mm_unpacklo_epi16(src_temp2_4x32b, one_8x16b);
          src_temp3_4x32b = _mm_unpacklo_epi16(src_temp3_4x32b, one_8x16b);

          src_temp0_4x32b = _mm_madd_epi16(src_temp0_4x32b, madd_const_8x16b);
          src_temp1_4x32b = _mm_madd_epi16(src_temp1_4x32b, madd_const_8x16b);
          src_temp2_4x32b = _mm_madd_epi16(src_temp2_4x32b, madd_const_8x16b);
          src_temp3_4x32b = _mm_madd_epi16(src_temp3_4x32b, madd_const_8x16b);

          /* Last 4 pixels */
          src_temp4_4x32b = _mm_unpacklo_epi16(src_temp4_4x32b, one_8x16b);
          src_temp5_4x32b = _mm_unpacklo_epi16(src_temp5_4x32b, one_8x16b);
          src_temp6_4x32b = _mm_unpacklo_epi16(src_temp6_4x32b, one_8x16b);
          src_temp7_4x32b = _mm_unpacklo_epi16(src_temp7_4x32b, one_8x16b);

          /* Last 4 pixels */
          src_temp4_4x32b = _mm_madd_epi16(src_temp4_4x32b, madd_const_8x16b);
          src_temp5_4x32b = _mm_madd_epi16(src_temp5_4x32b, madd_const_8x16b);
          src_temp6_4x32b = _mm_madd_epi16(src_temp6_4x32b, madd_const_8x16b);
          src_temp7_4x32b = _mm_madd_epi16(src_temp7_4x32b, madd_const_8x16b);

          /* (i4_tmp = i4_tmp >> shift) */ /* First 4 pixels */
          src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
          src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);
          src_temp2_4x32b = _mm_srai_epi32(src_temp2_4x32b, iShift);
          src_temp3_4x32b = _mm_srai_epi32(src_temp3_4x32b, iShift);

          /* (i4_tmp = i4_tmp >> shift) */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_srai_epi32(src_temp4_4x32b, iShift);
          src_temp5_4x32b = _mm_srai_epi32(src_temp5_4x32b, iShift);
          src_temp6_4x32b = _mm_srai_epi32(src_temp6_4x32b, iShift);
          src_temp7_4x32b = _mm_srai_epi32(src_temp7_4x32b, iShift);

          src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp4_4x32b);
          src_temp1_4x32b = _mm_packs_epi32(src_temp1_4x32b, src_temp5_4x32b);
          src_temp2_4x32b = _mm_packs_epi32(src_temp2_4x32b, src_temp6_4x32b);
          src_temp3_4x32b = _mm_packs_epi32(src_temp3_4x32b, src_temp7_4x32b);

          // Clipping
          src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
          src_temp1_4x32b = _mm_min_epi16(src_temp1_4x32b, clip_max_8x16b);
          src_temp2_4x32b = _mm_min_epi16(src_temp2_4x32b, clip_max_8x16b);
          src_temp3_4x32b = _mm_min_epi16(src_temp3_4x32b, clip_max_8x16b);

          src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);
          src_temp1_4x32b = _mm_max_epi16(src_temp1_4x32b, clip_min_8x16b);
          src_temp2_4x32b = _mm_max_epi16(src_temp2_4x32b, clip_min_8x16b);
          src_temp3_4x32b = _mm_max_epi16(src_temp3_4x32b, clip_min_8x16b);

          /* store eight 16-bit output values  */
          _mm_storeu_si128((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/
          _mm_storeu_si128((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/
          _mm_storeu_si128((__m128i*)(piDst + 2 * iDstStride), src_temp2_4x32b); /* row = 2*/
          _mm_storeu_si128((__m128i*)(piDst + 3 * iDstStride), src_temp3_4x32b); /* row = 3*/

          /* Update poInter */
          piSrc += 8;
          piDst += 8;
        } /* inner loop ends here */
        /* PoInter update */
        piSrc = piSrc - width + 4 * iSrcStride;
        piDst = piDst - width + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 3)) && (0 == (height & 3)))
    {
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        for (col = 0; col < width; col += 4)
        {
          /* for row =0 ,1,2,3*/
          /* row = 0 */ /*load 8 pixel values from 7:0 pos. relative to cur. pos.*/
          src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
          /* row = 1 */
          src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));
          /* row = 2 */
          src_temp2_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride));
          /* row = 3 */
          src_temp3_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride));

          /* First 4 pixels */
          src_temp0_4x32b = _mm_unpacklo_epi16(src_temp0_4x32b, one_8x16b);
          src_temp1_4x32b = _mm_unpacklo_epi16(src_temp1_4x32b, one_8x16b);
          src_temp2_4x32b = _mm_unpacklo_epi16(src_temp2_4x32b, one_8x16b);
          src_temp3_4x32b = _mm_unpacklo_epi16(src_temp3_4x32b, one_8x16b);

          src_temp0_4x32b = _mm_madd_epi16(src_temp0_4x32b, madd_const_8x16b);
          src_temp1_4x32b = _mm_madd_epi16(src_temp1_4x32b, madd_const_8x16b);
          src_temp2_4x32b = _mm_madd_epi16(src_temp2_4x32b, madd_const_8x16b);
          src_temp3_4x32b = _mm_madd_epi16(src_temp3_4x32b, madd_const_8x16b);

          /* (i4_tmp = i4_tmp >> shift) */
          src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
          src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);
          src_temp2_4x32b = _mm_srai_epi32(src_temp2_4x32b, iShift);
          src_temp3_4x32b = _mm_srai_epi32(src_temp3_4x32b, iShift);

          src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);
          src_temp2_4x32b = _mm_packs_epi32(src_temp2_4x32b, src_temp3_4x32b);

          // Clipping
          src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
          src_temp2_4x32b = _mm_min_epi16(src_temp2_4x32b, clip_max_8x16b);

          src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);
          src_temp2_4x32b = _mm_max_epi16(src_temp2_4x32b, clip_min_8x16b);

          _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/
          _mm_storel_epi64((__m128i*)(piDst + 2 * iDstStride), src_temp2_4x32b); /* row = 2*/

          src_temp1_4x32b = _mm_srli_si128(src_temp0_4x32b, 8);
          src_temp3_4x32b = _mm_srli_si128(src_temp2_4x32b, 8);

          _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/
          _mm_storel_epi64((__m128i*)(piDst + 3 * iDstStride), src_temp3_4x32b); /* row = 3*/

          /* Update poInter */
          piSrc += 4;
          piDst += 4;
        } /* inner loop ends here(4-output values in single iteration) */
        /* PoInter update */
        piSrc = piSrc - width + 4 * iDstStride;
        piDst = piDst - width + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 1)) && (0 == (height & 3)))
    {
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        /* for row =0 ,1,2,3*/
        /* row = 0 */ /*load 4 pixel values from 3:0 pos. relative to cur. pos.*/
        src_temp0_4x32b = _mm_loadl_epi64((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));
        /* row = 2 */
        src_temp2_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + 2 * iSrcStride));
        /* row = 3 */
        src_temp3_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + 3 * iSrcStride));

        src_temp0_4x32b = _mm_unpacklo_epi32(src_temp0_4x32b, src_temp1_4x32b);
        src_temp1_4x32b = _mm_unpacklo_epi32(src_temp2_4x32b, src_temp3_4x32b);

        src_temp0_4x32b = _mm_unpacklo_epi16(src_temp0_4x32b, one_8x16b);
        src_temp1_4x32b = _mm_unpacklo_epi16(src_temp1_4x32b, one_8x16b);

        src_temp0_4x32b = _mm_madd_epi16(src_temp0_4x32b, madd_const_8x16b);
        src_temp1_4x32b = _mm_madd_epi16(src_temp1_4x32b, madd_const_8x16b);

        src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
        src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);

        src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);

        src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
        src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

        *(Int *)piDst = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + 2 * iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + 3 * iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        /* PoInter update */
        piSrc = piSrc + 4 * iDstStride;
        piDst = piDst + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 3)) && (0 == (height & 1)))
    {
      for (col = 0; col < width; col += 4)
      {
        /* for row =0 ,1*/
        /* row = 0 */ /*load 8 pixel values from 7:0 pos. relative to cur. pos.*/
        src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));

        /* First 4 pixels */
        src_temp0_4x32b = _mm_unpacklo_epi16(src_temp0_4x32b, one_8x16b);
        src_temp1_4x32b = _mm_unpacklo_epi16(src_temp1_4x32b, one_8x16b);

        src_temp0_4x32b = _mm_madd_epi16(src_temp0_4x32b, madd_const_8x16b);
        src_temp1_4x32b = _mm_madd_epi16(src_temp1_4x32b, madd_const_8x16b);

        /* (i4_tmp = i4_tmp >> shift) */
        src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
        src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);

        src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);

        // Clipping
        src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
        src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

        _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/

        src_temp1_4x32b = _mm_srli_si128(src_temp0_4x32b, 8);
        _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/

        /* Update poInter */
        piSrc += 4;
        piDst += 4;
      } /* inner loop ends here(4-output values in single iteration) */
    }
    else if ((0 == (width & 1)) && (0 == (height & 1)))
    {

      /* for row =0 ,1*/
      /* row = 0 */ /*load 4 pixel values from 3:0 pos. relative to cur. pos.*/
      src_temp0_4x32b = _mm_loadl_epi64((__m128i*)(piSrc));
      /* row = 1 */
      src_temp1_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));

      src_temp0_4x32b = _mm_unpacklo_epi32(src_temp0_4x32b, src_temp1_4x32b);
      src_temp0_4x32b = _mm_unpacklo_epi16(src_temp0_4x32b, one_8x16b);

      src_temp0_4x32b = _mm_madd_epi16(src_temp0_4x32b, madd_const_8x16b);

      src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);

      src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, zero_8x16b);

      src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
      src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

      *(Int *)piDst = _mm_cvtsi128_si32(src_temp0_4x32b);
      src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);

      *(Int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);
    }
  }
  else // bit depth greater than 10
  {
    Int row, col, width, height;
    Pel m, M;

    /* all 128 bit registers are named with a suffix mxnb, where m is the */
    /* number of n bits packed in the register                            */
    __m128i src_temp0_4x32b, src_temp1_4x32b, src_temp2_4x32b, src_temp3_4x32b;
    __m128i wgt0_4x32b, off0_4x32b;
    __m128i clip_min_8x16b, clip_max_8x16b;
    __m128i zero_8x16b;

    width = cxWidth;
    height = cxHeight;

    assert(width % 2 == 0);
    assert(height % 2 == 0);

    m = g_ClipParam.min(compID);
    M = g_ClipParam.max(compID);

    // seting values in register
    wgt0_4x32b = _mm_set1_epi32(a);
    off0_4x32b = _mm_set1_epi32(b);
    clip_min_8x16b = _mm_set1_epi16(m);
    clip_max_8x16b = _mm_set1_epi16(M);
    zero_8x16b = _mm_setzero_si128();

    if ((0 == (width & 7)) && (0 == (height & 3))) /* width multiple of 8 case */
    {
      __m128i src_temp4_4x32b, src_temp5_4x32b, src_temp6_4x32b, src_temp7_4x32b;
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        for (col = 0; col < width; col += 8)
        {
          /* for row =0 ,1,2,3*/
          /* row = 0 */
          src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
          /* row = 1 */
          src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));
          /* row = 2 */
          src_temp2_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride));
          /* row = 3 */
          src_temp3_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride));

          /* row = 0 */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 4));
          /* row = 1 */
          src_temp5_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride + 4));
          /* row = 2 */
          src_temp6_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride + 4));
          /* row = 3 */
          src_temp7_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride + 4));

          /* considering pix. 4:0 by converting 16-Into 32 bit */ /* First 4 pixels */
          src_temp0_4x32b = _mm_cvtepi16_epi32(src_temp0_4x32b);
          src_temp1_4x32b = _mm_cvtepi16_epi32(src_temp1_4x32b);
          src_temp2_4x32b = _mm_cvtepi16_epi32(src_temp2_4x32b);
          src_temp3_4x32b = _mm_cvtepi16_epi32(src_temp3_4x32b);

          /* (i4_tmp = i4_tmp * wgt0) */ /* First 4 pixels */
          src_temp0_4x32b = _mm_mullo_epi32(src_temp0_4x32b, wgt0_4x32b);
          src_temp1_4x32b = _mm_mullo_epi32(src_temp1_4x32b, wgt0_4x32b);
          src_temp2_4x32b = _mm_mullo_epi32(src_temp2_4x32b, wgt0_4x32b);
          src_temp3_4x32b = _mm_mullo_epi32(src_temp3_4x32b, wgt0_4x32b);

          /* considering pix. 4:0 by converting 16-Into 32 bit */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_cvtepi16_epi32(src_temp4_4x32b);
          src_temp5_4x32b = _mm_cvtepi16_epi32(src_temp5_4x32b);
          src_temp6_4x32b = _mm_cvtepi16_epi32(src_temp6_4x32b);
          src_temp7_4x32b = _mm_cvtepi16_epi32(src_temp7_4x32b);

          /* (i4_tmp = i4_tmp * wgt0) */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_mullo_epi32(src_temp4_4x32b, wgt0_4x32b);
          src_temp5_4x32b = _mm_mullo_epi32(src_temp5_4x32b, wgt0_4x32b);
          src_temp6_4x32b = _mm_mullo_epi32(src_temp6_4x32b, wgt0_4x32b);
          src_temp7_4x32b = _mm_mullo_epi32(src_temp7_4x32b, wgt0_4x32b);

          /* (i4_tmp = i4_tmp >> shift) */ /* First 4 pixels */
          src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
          src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);
          src_temp2_4x32b = _mm_srai_epi32(src_temp2_4x32b, iShift);
          src_temp3_4x32b = _mm_srai_epi32(src_temp3_4x32b, iShift);

          /* (i4_tmp = i4_tmp >> shift) */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_srai_epi32(src_temp4_4x32b, iShift);
          src_temp5_4x32b = _mm_srai_epi32(src_temp5_4x32b, iShift);
          src_temp6_4x32b = _mm_srai_epi32(src_temp6_4x32b, iShift);
          src_temp7_4x32b = _mm_srai_epi32(src_temp7_4x32b, iShift);

          /* (i4_tmp = i4_tmp + off0); */ /* First 4 pixels */
          src_temp0_4x32b = _mm_add_epi32(src_temp0_4x32b, off0_4x32b);
          src_temp1_4x32b = _mm_add_epi32(src_temp1_4x32b, off0_4x32b);
          src_temp2_4x32b = _mm_add_epi32(src_temp2_4x32b, off0_4x32b);
          src_temp3_4x32b = _mm_add_epi32(src_temp3_4x32b, off0_4x32b);

          /* (i4_tmp = i4_tmp + off0); */ /* Last 4 pixels */
          src_temp4_4x32b = _mm_add_epi32(src_temp4_4x32b, off0_4x32b);
          src_temp5_4x32b = _mm_add_epi32(src_temp5_4x32b, off0_4x32b);
          src_temp6_4x32b = _mm_add_epi32(src_temp6_4x32b, off0_4x32b);
          src_temp7_4x32b = _mm_add_epi32(src_temp7_4x32b, off0_4x32b);

          src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp4_4x32b);
          src_temp1_4x32b = _mm_packs_epi32(src_temp1_4x32b, src_temp5_4x32b);
          src_temp2_4x32b = _mm_packs_epi32(src_temp2_4x32b, src_temp6_4x32b);
          src_temp3_4x32b = _mm_packs_epi32(src_temp3_4x32b, src_temp7_4x32b);

          // Clipping
          src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
          src_temp1_4x32b = _mm_min_epi16(src_temp1_4x32b, clip_max_8x16b);
          src_temp2_4x32b = _mm_min_epi16(src_temp2_4x32b, clip_max_8x16b);
          src_temp3_4x32b = _mm_min_epi16(src_temp3_4x32b, clip_max_8x16b);

          src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);
          src_temp1_4x32b = _mm_max_epi16(src_temp1_4x32b, clip_min_8x16b);
          src_temp2_4x32b = _mm_max_epi16(src_temp2_4x32b, clip_min_8x16b);
          src_temp3_4x32b = _mm_max_epi16(src_temp3_4x32b, clip_min_8x16b);

          /* store eight 16-bit output values  */
          _mm_storeu_si128((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/
          _mm_storeu_si128((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/
          _mm_storeu_si128((__m128i*)(piDst + 2 * iDstStride), src_temp2_4x32b); /* row = 2*/
          _mm_storeu_si128((__m128i*)(piDst + 3 * iDstStride), src_temp3_4x32b); /* row = 3*/

          /* Update poInter */
          piSrc += 8;
          piDst += 8;
        }
        /* PoInter update */
        piSrc = piSrc - width + 4 * iSrcStride;
        piDst = piDst - width + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 3)) && (0 == (height & 3)))
    {
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        for (col = 0; col < width; col += 4)
        {
          /* for row =0 ,1,2,3*/
          /* row = 0 */ /*load 8 pixel values from 7:0 pos. relative to cur. pos.*/
          src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
          /* row = 1 */
          src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));
          /* row = 2 */
          src_temp2_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride));
          /* row = 3 */
          src_temp3_4x32b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride));

          /* considering pix. 4:0 by converting 16-Into 32 bit */
          src_temp0_4x32b = _mm_cvtepi16_epi32(src_temp0_4x32b);
          src_temp1_4x32b = _mm_cvtepi16_epi32(src_temp1_4x32b);
          src_temp2_4x32b = _mm_cvtepi16_epi32(src_temp2_4x32b);
          src_temp3_4x32b = _mm_cvtepi16_epi32(src_temp3_4x32b);

          src_temp0_4x32b = _mm_mullo_epi32(src_temp0_4x32b, wgt0_4x32b);
          src_temp1_4x32b = _mm_mullo_epi32(src_temp1_4x32b, wgt0_4x32b);
          src_temp2_4x32b = _mm_mullo_epi32(src_temp2_4x32b, wgt0_4x32b);
          src_temp3_4x32b = _mm_mullo_epi32(src_temp3_4x32b, wgt0_4x32b);

          /* (i4_tmp = i4_tmp >> shift) */
          src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
          src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);
          src_temp2_4x32b = _mm_srai_epi32(src_temp2_4x32b, iShift);
          src_temp3_4x32b = _mm_srai_epi32(src_temp3_4x32b, iShift);

          /* (i4_tmp = i4_tmp + off0) */
          src_temp0_4x32b = _mm_add_epi32(src_temp0_4x32b, off0_4x32b);
          src_temp1_4x32b = _mm_add_epi32(src_temp1_4x32b, off0_4x32b);
          src_temp2_4x32b = _mm_add_epi32(src_temp2_4x32b, off0_4x32b);
          src_temp3_4x32b = _mm_add_epi32(src_temp3_4x32b, off0_4x32b);

          src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);
          src_temp2_4x32b = _mm_packs_epi32(src_temp2_4x32b, src_temp3_4x32b);

          // Clipping
          src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
          src_temp2_4x32b = _mm_min_epi16(src_temp2_4x32b, clip_max_8x16b);

          src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);
          src_temp2_4x32b = _mm_max_epi16(src_temp2_4x32b, clip_min_8x16b);

          _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/
          _mm_storel_epi64((__m128i*)(piDst + 2 * iDstStride), src_temp2_4x32b); /* row = 2*/

          src_temp1_4x32b = _mm_srli_si128(src_temp0_4x32b, 8);
          src_temp3_4x32b = _mm_srli_si128(src_temp2_4x32b, 8);

          _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/
          _mm_storel_epi64((__m128i*)(piDst + 3 * iDstStride), src_temp3_4x32b); /* row = 3*/

          /* Update poInter */
          piSrc += 4;
          piDst += 4;
        }
        /* PoInter update */
        piSrc = piSrc - width + 4 * iDstStride;
        piDst = piDst - width + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 1)) && (0 == (height & 3)))
    {
      /*  outer for loop starts from here */
      for (row = 0; row < height; row += 4)
      {
        /* for row =0 ,1,2,3*/
        /* row = 0 */ /*load 4 pixel values from 3:0 pos. relative to cur. pos.*/
        src_temp0_4x32b = _mm_loadl_epi64((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));
        /* row = 2 */
        src_temp2_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + 2 * iSrcStride));
        /* row = 3 */
        src_temp3_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + 3 * iSrcStride));

        /* considering pix. 4:0 by converting 16-Into 32 bit */
        src_temp0_4x32b = _mm_cvtepi16_epi32(src_temp0_4x32b);
        src_temp1_4x32b = _mm_cvtepi16_epi32(src_temp1_4x32b);
        src_temp2_4x32b = _mm_cvtepi16_epi32(src_temp2_4x32b);
        src_temp3_4x32b = _mm_cvtepi16_epi32(src_temp3_4x32b);

        src_temp0_4x32b = _mm_unpacklo_epi64(src_temp0_4x32b, src_temp1_4x32b);
        src_temp1_4x32b = _mm_unpacklo_epi64(src_temp2_4x32b, src_temp3_4x32b);

        src_temp0_4x32b = _mm_mullo_epi32(src_temp0_4x32b, wgt0_4x32b);
        src_temp1_4x32b = _mm_mullo_epi32(src_temp1_4x32b, wgt0_4x32b);

        src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
        src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);

        src_temp0_4x32b = _mm_add_epi32(src_temp0_4x32b, off0_4x32b);
        src_temp1_4x32b = _mm_add_epi32(src_temp1_4x32b, off0_4x32b);

        src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);

        src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
        src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

        *(Int *)piDst = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + 2 * iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
        *(Int *)(piDst + 3 * iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);

        /* PoInter update */
        piSrc = piSrc + 4 * iDstStride;
        piDst = piDst + 4 * iDstStride;
      }
    }
    else if ((0 == (width & 3)) && (0 == (height & 1)))
    {
      for (col = 0; col < width; col += 4)
      {
        /* for row =0 ,1*/
        /* row = 0 */ /*load 8 pixel values from 7:0 pos. relative to cur. pos.*/
        src_temp0_4x32b = _mm_loadu_si128((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_4x32b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));

        /* considering pix. 4:0 by converting 16-Into 32 bit */
        src_temp0_4x32b = _mm_cvtepi16_epi32(src_temp0_4x32b);
        src_temp1_4x32b = _mm_cvtepi16_epi32(src_temp1_4x32b);

        src_temp0_4x32b = _mm_mullo_epi32(src_temp0_4x32b, wgt0_4x32b);
        src_temp1_4x32b = _mm_mullo_epi32(src_temp1_4x32b, wgt0_4x32b);

        /* (i4_tmp = i4_tmp >> shift) */
        src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);
        src_temp1_4x32b = _mm_srai_epi32(src_temp1_4x32b, iShift);

        /* (i4_tmp = i4_tmp + off0) */
        src_temp0_4x32b = _mm_add_epi32(src_temp0_4x32b, off0_4x32b);
        src_temp1_4x32b = _mm_add_epi32(src_temp1_4x32b, off0_4x32b);

        src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, src_temp1_4x32b);

        // Clipping
        src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
        src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

        _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_4x32b); /* row = 0*/

        src_temp1_4x32b = _mm_srli_si128(src_temp0_4x32b, 8);

        _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_4x32b); /* row = 1*/
        /* Update poInter */
        piSrc += 4;
        piDst += 4;
      }
    }
    else if ((0 == (width & 1)) && (0 == (height & 1)))
    {
      /* for row =0 ,1*/
      /* row = 0 */ /*load 4 pixel values from 7:0 pos. relative to cur. pos.*/
      src_temp0_4x32b = _mm_loadl_epi64((__m128i*)(piSrc));
      /* row = 1 */
      src_temp1_4x32b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));

      /* considering pix. 4:0 by converting 16-Into 32 bit */
      src_temp0_4x32b = _mm_cvtepi16_epi32(src_temp0_4x32b);
      src_temp1_4x32b = _mm_cvtepi16_epi32(src_temp1_4x32b);

      src_temp0_4x32b = _mm_unpacklo_epi64(src_temp0_4x32b, src_temp1_4x32b);

      src_temp0_4x32b = _mm_mullo_epi32(src_temp0_4x32b, wgt0_4x32b);

      src_temp0_4x32b = _mm_srai_epi32(src_temp0_4x32b, iShift);

      src_temp0_4x32b = _mm_add_epi32(src_temp0_4x32b, off0_4x32b);

      src_temp0_4x32b = _mm_packs_epi32(src_temp0_4x32b, zero_8x16b);

      src_temp0_4x32b = _mm_min_epi16(src_temp0_4x32b, clip_max_8x16b);
      src_temp0_4x32b = _mm_max_epi16(src_temp0_4x32b, clip_min_8x16b);

      *(Int *)piDst = _mm_cvtsi128_si32(src_temp0_4x32b);

      src_temp0_4x32b = _mm_srli_si128(src_temp0_4x32b, 4);
      *(Int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp0_4x32b);
    }
  }
}
#endif

#if ITH_C3_SCALESUB==C_FXN
void scaleSub
(
Pel *piSrc,
Int iSrcStride,
Pel *piDst,
Int iDstStride,
UInt cxWidth,
UInt cxHeight,
const Int iShift,
Int iOffset
)
{
  int i, j;
  for (i = 0; i < cxHeight; i++)
  {
    for (j = 0; j < cxWidth; j++)
    {
      Short val = piSrc[j] << iShift;
      piDst[j] = val - iOffset;
    }
    piDst += iDstStride;
    piSrc += iSrcStride;
  }
}
#elif ITH_C3_SCALESUB==SIMD_FXN
void scaleSub
(
Pel *piSrc,
Int iSrcStride,
Pel *piDst,
Int iDstStride,
UInt cxWidth,
UInt cxHeight,
const Int iShift,
Int iOffset
)
{
  int row, col, width, height;

  /* all 128 bit registers are named with a suffix mxnb, where m is the */
  /* number of n bits packed in the register                            */
  __m128i src_temp0_8x16b, src_temp1_8x16b, src_temp2_8x16b, src_temp3_8x16b;
  __m128i off0_8x16b, shift_4x32b;

  width = cxWidth;
  height = cxHeight;

  assert(width % 2 == 0);
  assert(height % 2 == 0);

  // seting values in register
  off0_8x16b = _mm_set1_epi16(iOffset);
  shift_4x32b = _mm_set_epi32(0, 0, 0, iShift);//Deepa

  if ((0 == (width & 7)) && (0 == (height & 3)))
  {
    /*  outer for loop starts from here */
    for (row = 0; row < height; row += 4)
    {
      for (col = 0; col < width; col += 8)
      {
        /* for row =0 ,1,2,3*/
        /* row = 0 */
        src_temp0_8x16b = _mm_loadu_si128((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_8x16b = _mm_loadu_si128((__m128i*)(piSrc + iSrcStride));
        /* row = 2 */
        src_temp2_8x16b = _mm_loadu_si128((__m128i*)(piSrc + 2 * iSrcStride));
        /* row = 3 */
        src_temp3_8x16b = _mm_loadu_si128((__m128i*)(piSrc + 3 * iSrcStride));

        /* (i4_tmp = i4_tmp << shift) */ /* First 8 pixels */
        src_temp0_8x16b = _mm_sll_epi16(src_temp0_8x16b, shift_4x32b);
        src_temp1_8x16b = _mm_sll_epi16(src_temp1_8x16b, shift_4x32b);
        src_temp2_8x16b = _mm_sll_epi16(src_temp2_8x16b, shift_4x32b);
        src_temp3_8x16b = _mm_sll_epi16(src_temp3_8x16b, shift_4x32b);

        /* (i4_tmp = i4_tmp - off0); */ /* First 8 pixels */
        src_temp0_8x16b = _mm_sub_epi16(src_temp0_8x16b, off0_8x16b);
        src_temp1_8x16b = _mm_sub_epi16(src_temp1_8x16b, off0_8x16b);
        src_temp2_8x16b = _mm_sub_epi16(src_temp2_8x16b, off0_8x16b);
        src_temp3_8x16b = _mm_sub_epi16(src_temp3_8x16b, off0_8x16b);

        /* store eight 16-bit output values  */
        _mm_storeu_si128((__m128i*)(piDst + 0 * iDstStride), src_temp0_8x16b); /* row = 0*/
        _mm_storeu_si128((__m128i*)(piDst + 1 * iDstStride), src_temp1_8x16b); /* row = 1*/
        _mm_storeu_si128((__m128i*)(piDst + 2 * iDstStride), src_temp2_8x16b); /* row = 2*/
        _mm_storeu_si128((__m128i*)(piDst + 3 * iDstStride), src_temp3_8x16b); /* row = 3*/

        /* Update pointer */
        piSrc += 8;
        piDst += 8;
      } /* inner loop ends here */
      /* Pointer update */
      piSrc = piSrc - width + 4 * iSrcStride;
      piDst = piDst - width + 4 * iDstStride;
    }
  }
  else if ((0 == (width & 3)) && (0 == (height & 3)))
  {
    /*  outer for loop starts from here */
    for (row = 0; row < height; row += 4)
    {
      for (col = 0; col < width; col += 4)
      {
        /* for row =0 ,1,2,3*/
        /* row = 0 */
        src_temp0_8x16b = _mm_loadl_epi64((__m128i*)(piSrc));
        /* row = 1 */
        src_temp1_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));
        /* row = 2 */
        src_temp2_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + 2 * iSrcStride));
        /* row = 3 */
        src_temp3_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + 3 * iSrcStride));

        /* (i4_tmp = i4_tmp << shift) */ /* First 8 pixels */
        src_temp0_8x16b = _mm_sll_epi16(src_temp0_8x16b, shift_4x32b);
        src_temp1_8x16b = _mm_sll_epi16(src_temp1_8x16b, shift_4x32b);
        src_temp2_8x16b = _mm_sll_epi16(src_temp2_8x16b, shift_4x32b);
        src_temp3_8x16b = _mm_sll_epi16(src_temp3_8x16b, shift_4x32b);

        /* (i4_tmp = i4_tmp - off0); */ /* First 8 pixels */
        src_temp0_8x16b = _mm_sub_epi16(src_temp0_8x16b, off0_8x16b);
        src_temp1_8x16b = _mm_sub_epi16(src_temp1_8x16b, off0_8x16b);
        src_temp2_8x16b = _mm_sub_epi16(src_temp2_8x16b, off0_8x16b);
        src_temp3_8x16b = _mm_sub_epi16(src_temp3_8x16b, off0_8x16b);

        /* store four 16-bit output values  */
        _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_8x16b); /* row = 0*/
        _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_8x16b); /* row = 1*/
        _mm_storel_epi64((__m128i*)(piDst + 2 * iDstStride), src_temp2_8x16b); /* row = 2*/
        _mm_storel_epi64((__m128i*)(piDst + 3 * iDstStride), src_temp3_8x16b); /* row = 3*/

        /* Update pointer */
        piSrc += 4;
        piDst += 4;
      }
      /* Pointer update */
      piSrc = piSrc - width + 4 * iSrcStride;
      piDst = piDst - width + 4 * iDstStride;
    }
  }
  else if ((0 == (width & 1)) && (0 == (height & 3)))
  {
    /*  outer for loop starts from here */
    for (row = 0; row < height; row += 4)
    {
      /* for row =0 ,1,2,3*/
      /* row = 0 */
      src_temp0_8x16b = _mm_loadl_epi64((__m128i*)(piSrc));
      /* row = 1 */
      src_temp1_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));
      /* row = 2 */
      src_temp2_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + 2 * iSrcStride));
      /* row = 3 */
      src_temp3_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + 3 * iSrcStride));

      /* (i4_tmp = i4_tmp << shift) */ /* First 8 pixels */
      src_temp0_8x16b = _mm_sll_epi16(src_temp0_8x16b, shift_4x32b);
      src_temp1_8x16b = _mm_sll_epi16(src_temp1_8x16b, shift_4x32b);
      src_temp2_8x16b = _mm_sll_epi16(src_temp2_8x16b, shift_4x32b);
      src_temp3_8x16b = _mm_sll_epi16(src_temp3_8x16b, shift_4x32b);

      /* (i4_tmp = i4_tmp - off0); */ /* First 8 pixels */
      src_temp0_8x16b = _mm_sub_epi16(src_temp0_8x16b, off0_8x16b);
      src_temp1_8x16b = _mm_sub_epi16(src_temp1_8x16b, off0_8x16b);
      src_temp2_8x16b = _mm_sub_epi16(src_temp2_8x16b, off0_8x16b);
      src_temp3_8x16b = _mm_sub_epi16(src_temp3_8x16b, off0_8x16b);

      *(int *)piDst = _mm_cvtsi128_si32(src_temp0_8x16b);
      *(int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp1_8x16b);
      *(int *)(piDst + 2 * iDstStride) = _mm_cvtsi128_si32(src_temp2_8x16b);
      *(int *)(piDst + 3 * iDstStride) = _mm_cvtsi128_si32(src_temp3_8x16b);

      /* Pointer update */
      piSrc = piSrc + 4 * iSrcStride;
      piDst = piDst + 4 * iDstStride;
    }
  }
  else if ((0 == (width & 3)) && (0 == (height & 1)))
  {
    for (col = 0; col < width; col += 4)
    {
      /* for row =0 ,1,2,3*/
      /* row = 0 */
      src_temp0_8x16b = _mm_loadl_epi64((__m128i*)(piSrc));
      /* row = 1 */
      src_temp1_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));

      /* (i4_tmp = i4_tmp << shift) */ /* First 8 pixels */
      src_temp0_8x16b = _mm_sll_epi16(src_temp0_8x16b, shift_4x32b);
      src_temp1_8x16b = _mm_sll_epi16(src_temp1_8x16b, shift_4x32b);

      /* (i4_tmp = i4_tmp - off0); */ /* First 8 pixels */
      src_temp0_8x16b = _mm_sub_epi16(src_temp0_8x16b, off0_8x16b);
      src_temp1_8x16b = _mm_sub_epi16(src_temp1_8x16b, off0_8x16b);

      /* store eight 16-bit output values  */
      _mm_storel_epi64((__m128i*)(piDst + 0 * iDstStride), src_temp0_8x16b); /* row = 0*/
      _mm_storel_epi64((__m128i*)(piDst + 1 * iDstStride), src_temp1_8x16b); /* row = 1*/

      /* Update pointer */
      piSrc += 4;
      piDst += 4;
    } /* inner loop ends here */
  }
  else if ((0 == (width & 1)) && (0 == (height & 1)))
  {
    /* for row =0 ,1,2,3*/
    /* row = 0 */
    src_temp0_8x16b = _mm_loadl_epi64((__m128i*)(piSrc));
    /* row = 1 */
    src_temp1_8x16b = _mm_loadl_epi64((__m128i*)(piSrc + iSrcStride));

    /* (i4_tmp = i4_tmp << shift) */ /* First 8 pixels */
    src_temp0_8x16b = _mm_sll_epi16(src_temp0_8x16b, shift_4x32b);
    src_temp1_8x16b = _mm_sll_epi16(src_temp1_8x16b, shift_4x32b);

    /* (i4_tmp = i4_tmp - off0); */ /* First 8 pixels */
    src_temp0_8x16b = _mm_sub_epi16(src_temp0_8x16b, off0_8x16b);
    src_temp1_8x16b = _mm_sub_epi16(src_temp1_8x16b, off0_8x16b);

    *(int *)piDst = _mm_cvtsi128_si32(src_temp0_8x16b);
    *(int *)(piDst + iDstStride) = _mm_cvtsi128_si32(src_temp1_8x16b);
  }
}
#endif

void countZeroRowsZeroCols
(
UInt uiWidth,
UInt uiHeight,
TCoeff *pRef,
UInt uiSrcStride,
UInt *puiZerosWidth,
UInt *puiZerosHeight
)
{
  TCoeff *pSrc = pRef + ((uiHeight - 1) * uiSrcStride) + uiWidth - 1;
  TCoeff *pSrc_start = pSrc;
  Int iRow, iCol;
  UInt zerosWidth = 0;
  UInt zerosHeight = 0;
  Bool bBreakFlag = false;

  for (iRow = 0; (iRow < uiHeight) && !bBreakFlag; iRow++)
  {
    for (iCol = 0; iCol < uiWidth; iCol++)
    {
      if (*(pSrc--) != 0)
      {
        bBreakFlag = true;
        zerosWidth = iCol;
        zerosHeight = iRow;
        break;
      }
    }
    pSrc_start -= uiSrcStride;
    pSrc = pSrc_start;
  }

  bBreakFlag = false;
  pSrc = pRef + ((uiHeight - zerosHeight - 1) * uiSrcStride) + uiWidth - 1;
  pSrc_start = pSrc;

  for (iCol = 0; (iCol < zerosWidth) && !bBreakFlag; iCol++)
  {
    for (iRow = 0; iRow < (uiHeight - zerosHeight); iRow++)
    {
      if (*(pSrc) != 0)
      {
        bBreakFlag = true;
        zerosWidth = iCol;
        break;
      }
      pSrc -= uiSrcStride;
    }
    pSrc_start--;
    pSrc = pSrc_start;
  }

  *(puiZerosWidth) = zerosWidth;
  *(puiZerosHeight) = zerosHeight;
}

#if ITH_C3_LEFTSHIFTSUB==C_FXN
void leftShiftSub
(
const Pel *src,
const Int srcStride,
Pel *dst,
const Int dstStride,
const Int width,
const Int height,
const Int shift,
const Pel offset
)
{
  int row, col;
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
#elif  ITH_C3_LEFTSHIFTSUB==SIMD_FXN
void leftShiftSub
(
const Pel *pSrc,
const Int iSrcStride,
Pel *pDst,
const Int iDstStride,
const Int iBlkWd,
const Int iBlkHt,
const Int shift,
const Pel offset
)
{
  Int iRows, iCols;
  const Pel *inp_copy;
  Pel *dst_copy;
  Int iRemWd = iBlkWd;

  __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;
  __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;
  __m128i offset_16x8b, shift_16x8b;

  offset_16x8b = _mm_set1_epi16(offset);
  shift_16x8b = _mm_set_epi32(0, 0, 0, shift);

  inp_copy = pSrc;
  dst_copy = pDst;

  assert(shift >= 0);

  if (iRemWd >= 8)
  {
    Int iRemHt = iBlkHt % 8;
    Int iCount = 0;
    iRemWd &= 0x7;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)	/*loop over height of 8*/
    {
      for (iCols = iBlkWd; iCols > 7; iCols -= 8)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride));
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride));
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride));
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride));
        src4_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 4 * iSrcStride));
        src5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 5 * iSrcStride));
        src6_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 6 * iSrcStride));
        src7_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 7 * iSrcStride));

        src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
        src1_16x8b = _mm_sll_epi16(src1_16x8b, shift_16x8b);
        src2_16x8b = _mm_sll_epi16(src2_16x8b, shift_16x8b);
        src3_16x8b = _mm_sll_epi16(src3_16x8b, shift_16x8b);
        src4_16x8b = _mm_sll_epi16(src4_16x8b, shift_16x8b);
        src5_16x8b = _mm_sll_epi16(src5_16x8b, shift_16x8b);
        src6_16x8b = _mm_sll_epi16(src6_16x8b, shift_16x8b);
        src7_16x8b = _mm_sll_epi16(src7_16x8b, shift_16x8b);

        src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
        src1_16x8b = _mm_sub_epi16(src1_16x8b, offset_16x8b);
        src2_16x8b = _mm_sub_epi16(src2_16x8b, offset_16x8b);
        src3_16x8b = _mm_sub_epi16(src3_16x8b, offset_16x8b);
        src4_16x8b = _mm_sub_epi16(src4_16x8b, offset_16x8b);
        src5_16x8b = _mm_sub_epi16(src5_16x8b, offset_16x8b);
        src6_16x8b = _mm_sub_epi16(src6_16x8b, offset_16x8b);
        src7_16x8b = _mm_sub_epi16(src7_16x8b, offset_16x8b);

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 4 * iDstStride), src4_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 5 * iDstStride), src5_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 6 * iDstStride), src6_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 7 * iDstStride), src7_16x8b);

        iCount += 8;
      }

      iCount = 0;
      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);

    }

    if (iRemHt > 3) /*loop over height of 4*/
    {
      iRemHt = iRemHt % 4;

      for (iCols = iBlkWd; iCols > 7; iCols -= 8)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride));
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride));
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride));
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride));

        src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
        src1_16x8b = _mm_sll_epi16(src1_16x8b, shift_16x8b);
        src2_16x8b = _mm_sll_epi16(src2_16x8b, shift_16x8b);
        src3_16x8b = _mm_sll_epi16(src3_16x8b, shift_16x8b);

        src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
        src1_16x8b = _mm_sub_epi16(src1_16x8b, offset_16x8b);
        src2_16x8b = _mm_sub_epi16(src2_16x8b, offset_16x8b);
        src3_16x8b = _mm_sub_epi16(src3_16x8b, offset_16x8b);

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);

        iCount += 8;
      }
      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)	/*loop over height of 1*/
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {
        iCount = 0;
        for (iCols = iBlkWd; iCols > 7; iCols -= 8)
        {
          src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount));
          src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
          src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
          _mm_storeu_si128((__m128i*)(dst_copy + iCount), src0_16x8b);
          iCount += 8;
        }
        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }

    }
  }

  inp_copy = pSrc + ((iBlkWd / 8) * 8);
  dst_copy = pDst + ((iBlkWd / 8) * 8);

  if (iRemWd > 3)
  {
    Int iRemHt = iBlkHt % 8;
    iRemWd &= 0x3;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)
    {

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride));
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride));
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride));
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride));
      src4_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 4 * iSrcStride));
      src5_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 5 * iSrcStride));
      src6_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 6 * iSrcStride));
      src7_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 7 * iSrcStride));

      src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
      src1_16x8b = _mm_sll_epi16(src1_16x8b, shift_16x8b);
      src2_16x8b = _mm_sll_epi16(src2_16x8b, shift_16x8b);
      src3_16x8b = _mm_sll_epi16(src3_16x8b, shift_16x8b);
      src4_16x8b = _mm_sll_epi16(src4_16x8b, shift_16x8b);
      src5_16x8b = _mm_sll_epi16(src5_16x8b, shift_16x8b);
      src6_16x8b = _mm_sll_epi16(src6_16x8b, shift_16x8b);
      src7_16x8b = _mm_sll_epi16(src7_16x8b, shift_16x8b);

      src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
      src1_16x8b = _mm_sub_epi16(src1_16x8b, offset_16x8b);
      src2_16x8b = _mm_sub_epi16(src2_16x8b, offset_16x8b);
      src3_16x8b = _mm_sub_epi16(src3_16x8b, offset_16x8b);
      src4_16x8b = _mm_sub_epi16(src4_16x8b, offset_16x8b);
      src5_16x8b = _mm_sub_epi16(src5_16x8b, offset_16x8b);
      src6_16x8b = _mm_sub_epi16(src6_16x8b, offset_16x8b);
      src7_16x8b = _mm_sub_epi16(src7_16x8b, offset_16x8b);

      _mm_storel_epi64((__m128i *)(dst_copy + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 3 * iDstStride), src3_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 4 * iDstStride), src4_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 5 * iDstStride), src5_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 6 * iDstStride), src6_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 7 * iDstStride), src7_16x8b);

      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);
    }

    if (iRemHt > 3)
    {
      iRemHt = iRemHt % 4;

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride));
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride));
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride));
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride));

      src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
      src1_16x8b = _mm_sll_epi16(src1_16x8b, shift_16x8b);
      src2_16x8b = _mm_sll_epi16(src2_16x8b, shift_16x8b);
      src3_16x8b = _mm_sll_epi16(src3_16x8b, shift_16x8b);

      src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
      src1_16x8b = _mm_sub_epi16(src1_16x8b, offset_16x8b);
      src2_16x8b = _mm_sub_epi16(src2_16x8b, offset_16x8b);
      src3_16x8b = _mm_sub_epi16(src3_16x8b, offset_16x8b);

      _mm_storel_epi64((__m128i *)(dst_copy + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 3 * iDstStride), src3_16x8b);

      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {

        src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy));
        src0_16x8b = _mm_sll_epi16(src0_16x8b, shift_16x8b);
        src0_16x8b = _mm_sub_epi16(src0_16x8b, offset_16x8b);
        _mm_storel_epi64((__m128i *)dst_copy, src0_16x8b);

        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }
    }
  }

  inp_copy = pSrc + ((iBlkWd / 4) * 4);
  dst_copy = pDst + ((iBlkWd / 4) * 4);

  if (iRemWd)
  {
    for (Int i = 0; i < iBlkHt; i++)
    {
      for (Int j = 0; j < iRemWd; j++)
      {
        *(dst_copy + j) = ((*(inp_copy + j)) << shift) - offset;
      }

      dst_copy += iDstStride;
      inp_copy += iSrcStride;
    }
  }
}
#endif

#if ITH_C3_RIGHTSHIFTCLIP==C_FXN
void rightShiftClip
(
const Pel *src,
const Int srcStride,
Pel *dst,
const Int dstStride,
const Int width,
const Int height,
const Int shift,
const Int offset,
const Pel minVal,
const Pel maxVal,
ComponentID compID
)
{
  int row, col;
  for (row = 0; row < height; row++)
  {
    for (col = 0; col < width; col++)
    {
      Pel val = src[col];
      val = rightShift_round((val + offset), shift);
      dst[col] = Clip3<Pel>(minVal, maxVal, val);
    }

    src += srcStride;
    dst += dstStride;
  }
}
#elif ITH_C3_RIGHTSHIFTCLIP==SIMD_FXN
void rightShiftClip
(
const Pel *pSrc,
const Int iSrcStride,
Pel *pDst,
const Int iDstStride,
const Int iBlkWd,
const Int iBlkHt,
const Int shift,
const Int offset,
const Pel minVal,
const Pel maxVal,
ComponentID compID
)
{
  Int iRows, iCols;
  const Pel *inp_copy;
  Pel *dst_copy;
  Int iRemWd = iBlkWd;
  Pel m = minVal;
  Pel M = maxVal;

  __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;
  __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;
  __m128i src8_16x8b, src9_16x8b, src10_16x8b, src11_16x8b;
  __m128i src12_16x8b, src13_16x8b, src14_16x8b, src15_16x8b;
  __m128i shift_32x4b, add_32x4b, clip_max_8x16b, clip_min_8x16b;

  assert(shift >= 0);

  shift_32x4b = _mm_set_epi32(0, 0, 0, shift);
  add_32x4b = _mm_set_epi32((1 << (shift - 1)) + offset,
                            (1 << (shift - 1)) + offset,
                            (1 << (shift - 1)) + offset,
                            (1 << (shift - 1)) + offset);

  clip_max_8x16b = _mm_set1_epi16(M);
  clip_min_8x16b = _mm_set1_epi16(m);

  inp_copy = pSrc;
  dst_copy = pDst;

  if (iRemWd >= 8)
  {
    Int iRemHt = iBlkHt % 8;
    Int iCount = 0;
    iRemWd &= 0x7;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)	/*loop over height of 8*/
    {
      for (iCols = iBlkWd; iCols > 7; iCols -= 8)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride));
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride));
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride));
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride));
        src4_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 4 * iSrcStride));
        src5_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 5 * iSrcStride));
        src6_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 6 * iSrcStride));
        src7_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 7 * iSrcStride));

        src8_16x8b = _mm_srli_si128(src0_16x8b, 8);
        src9_16x8b = _mm_srli_si128(src1_16x8b, 8);
        src10_16x8b = _mm_srli_si128(src2_16x8b, 8);
        src11_16x8b = _mm_srli_si128(src3_16x8b, 8);
        src12_16x8b = _mm_srli_si128(src4_16x8b, 8);
        src13_16x8b = _mm_srli_si128(src5_16x8b, 8);
        src14_16x8b = _mm_srli_si128(src6_16x8b, 8);
        src15_16x8b = _mm_srli_si128(src7_16x8b, 8);

        src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
        src1_16x8b = _mm_cvtepi16_epi32(src1_16x8b);
        src2_16x8b = _mm_cvtepi16_epi32(src2_16x8b);
        src3_16x8b = _mm_cvtepi16_epi32(src3_16x8b);
        src4_16x8b = _mm_cvtepi16_epi32(src4_16x8b);
        src5_16x8b = _mm_cvtepi16_epi32(src5_16x8b);
        src6_16x8b = _mm_cvtepi16_epi32(src6_16x8b);
        src7_16x8b = _mm_cvtepi16_epi32(src7_16x8b);

        src8_16x8b = _mm_cvtepi16_epi32(src8_16x8b);
        src9_16x8b = _mm_cvtepi16_epi32(src9_16x8b);
        src10_16x8b = _mm_cvtepi16_epi32(src10_16x8b);
        src11_16x8b = _mm_cvtepi16_epi32(src11_16x8b);
        src12_16x8b = _mm_cvtepi16_epi32(src12_16x8b);
        src13_16x8b = _mm_cvtepi16_epi32(src13_16x8b);
        src14_16x8b = _mm_cvtepi16_epi32(src14_16x8b);
        src15_16x8b = _mm_cvtepi16_epi32(src15_16x8b);

        src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
        src1_16x8b = _mm_add_epi32(src1_16x8b, add_32x4b);
        src2_16x8b = _mm_add_epi32(src2_16x8b, add_32x4b);
        src3_16x8b = _mm_add_epi32(src3_16x8b, add_32x4b);
        src4_16x8b = _mm_add_epi32(src4_16x8b, add_32x4b);
        src5_16x8b = _mm_add_epi32(src5_16x8b, add_32x4b);
        src6_16x8b = _mm_add_epi32(src6_16x8b, add_32x4b);
        src7_16x8b = _mm_add_epi32(src7_16x8b, add_32x4b);

        src8_16x8b = _mm_add_epi32(src8_16x8b, add_32x4b);
        src9_16x8b = _mm_add_epi32(src9_16x8b, add_32x4b);
        src10_16x8b = _mm_add_epi32(src10_16x8b, add_32x4b);
        src11_16x8b = _mm_add_epi32(src11_16x8b, add_32x4b);
        src12_16x8b = _mm_add_epi32(src12_16x8b, add_32x4b);
        src13_16x8b = _mm_add_epi32(src13_16x8b, add_32x4b);
        src14_16x8b = _mm_add_epi32(src14_16x8b, add_32x4b);
        src15_16x8b = _mm_add_epi32(src15_16x8b, add_32x4b);

        src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
        src1_16x8b = _mm_sra_epi32(src1_16x8b, shift_32x4b);
        src2_16x8b = _mm_sra_epi32(src2_16x8b, shift_32x4b);
        src3_16x8b = _mm_sra_epi32(src3_16x8b, shift_32x4b);
        src4_16x8b = _mm_sra_epi32(src4_16x8b, shift_32x4b);
        src5_16x8b = _mm_sra_epi32(src5_16x8b, shift_32x4b);
        src6_16x8b = _mm_sra_epi32(src6_16x8b, shift_32x4b);
        src7_16x8b = _mm_sra_epi32(src7_16x8b, shift_32x4b);

        src8_16x8b = _mm_sra_epi32(src8_16x8b, shift_32x4b);
        src9_16x8b = _mm_sra_epi32(src9_16x8b, shift_32x4b);
        src10_16x8b = _mm_sra_epi32(src10_16x8b, shift_32x4b);
        src11_16x8b = _mm_sra_epi32(src11_16x8b, shift_32x4b);
        src12_16x8b = _mm_sra_epi32(src12_16x8b, shift_32x4b);
        src13_16x8b = _mm_sra_epi32(src13_16x8b, shift_32x4b);
        src14_16x8b = _mm_sra_epi32(src14_16x8b, shift_32x4b);
        src15_16x8b = _mm_sra_epi32(src15_16x8b, shift_32x4b);

        src0_16x8b = _mm_packs_epi32(src0_16x8b, src8_16x8b);
        src1_16x8b = _mm_packs_epi32(src1_16x8b, src9_16x8b);
        src2_16x8b = _mm_packs_epi32(src2_16x8b, src10_16x8b);
        src3_16x8b = _mm_packs_epi32(src3_16x8b, src11_16x8b);
        src4_16x8b = _mm_packs_epi32(src4_16x8b, src12_16x8b);
        src5_16x8b = _mm_packs_epi32(src5_16x8b, src13_16x8b);
        src6_16x8b = _mm_packs_epi32(src6_16x8b, src14_16x8b);
        src7_16x8b = _mm_packs_epi32(src7_16x8b, src15_16x8b);

        src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
        src1_16x8b = _mm_min_epi16(src1_16x8b, clip_max_8x16b);
        src2_16x8b = _mm_min_epi16(src2_16x8b, clip_max_8x16b);
        src3_16x8b = _mm_min_epi16(src3_16x8b, clip_max_8x16b);
        src4_16x8b = _mm_min_epi16(src4_16x8b, clip_max_8x16b);
        src5_16x8b = _mm_min_epi16(src5_16x8b, clip_max_8x16b);
        src6_16x8b = _mm_min_epi16(src6_16x8b, clip_max_8x16b);
        src7_16x8b = _mm_min_epi16(src7_16x8b, clip_max_8x16b);

        src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);
        src1_16x8b = _mm_max_epi16(src1_16x8b, clip_min_8x16b);
        src2_16x8b = _mm_max_epi16(src2_16x8b, clip_min_8x16b);
        src3_16x8b = _mm_max_epi16(src3_16x8b, clip_min_8x16b);
        src4_16x8b = _mm_max_epi16(src4_16x8b, clip_min_8x16b);
        src5_16x8b = _mm_max_epi16(src5_16x8b, clip_min_8x16b);
        src6_16x8b = _mm_max_epi16(src6_16x8b, clip_min_8x16b);
        src7_16x8b = _mm_max_epi16(src7_16x8b, clip_min_8x16b);

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 4 * iDstStride), src4_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 5 * iDstStride), src5_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 6 * iDstStride), src6_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 7 * iDstStride), src7_16x8b);

        iCount += 8;
      }

      iCount = 0;
      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);

    }

    if (iRemHt > 3) /*loop over height of 4*/
    {
      iRemHt = iRemHt % 4;

      for (iCols = iBlkWd; iCols > 7; iCols -= 8)
      {
        src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 0 * iSrcStride));
        src1_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 1 * iSrcStride));
        src2_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 2 * iSrcStride));
        src3_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount + 3 * iSrcStride));

        src4_16x8b = _mm_srli_si128(src0_16x8b, 8);
        src5_16x8b = _mm_srli_si128(src1_16x8b, 8);
        src6_16x8b = _mm_srli_si128(src2_16x8b, 8);
        src7_16x8b = _mm_srli_si128(src3_16x8b, 8);

        src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
        src1_16x8b = _mm_cvtepi16_epi32(src1_16x8b);
        src2_16x8b = _mm_cvtepi16_epi32(src2_16x8b);
        src3_16x8b = _mm_cvtepi16_epi32(src3_16x8b);
        src4_16x8b = _mm_cvtepi16_epi32(src4_16x8b);
        src5_16x8b = _mm_cvtepi16_epi32(src5_16x8b);
        src6_16x8b = _mm_cvtepi16_epi32(src6_16x8b);
        src7_16x8b = _mm_cvtepi16_epi32(src7_16x8b);

        src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
        src1_16x8b = _mm_add_epi32(src1_16x8b, add_32x4b);
        src2_16x8b = _mm_add_epi32(src2_16x8b, add_32x4b);
        src3_16x8b = _mm_add_epi32(src3_16x8b, add_32x4b);
        src4_16x8b = _mm_add_epi32(src4_16x8b, add_32x4b);
        src5_16x8b = _mm_add_epi32(src5_16x8b, add_32x4b);
        src6_16x8b = _mm_add_epi32(src6_16x8b, add_32x4b);
        src7_16x8b = _mm_add_epi32(src7_16x8b, add_32x4b);

        src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
        src1_16x8b = _mm_sra_epi32(src1_16x8b, shift_32x4b);
        src2_16x8b = _mm_sra_epi32(src2_16x8b, shift_32x4b);
        src3_16x8b = _mm_sra_epi32(src3_16x8b, shift_32x4b);
        src4_16x8b = _mm_sra_epi32(src4_16x8b, shift_32x4b);
        src5_16x8b = _mm_sra_epi32(src5_16x8b, shift_32x4b);
        src6_16x8b = _mm_sra_epi32(src6_16x8b, shift_32x4b);
        src7_16x8b = _mm_sra_epi32(src7_16x8b, shift_32x4b);

        src0_16x8b = _mm_packs_epi32(src0_16x8b, src4_16x8b);
        src1_16x8b = _mm_packs_epi32(src1_16x8b, src5_16x8b);
        src2_16x8b = _mm_packs_epi32(src2_16x8b, src6_16x8b);
        src3_16x8b = _mm_packs_epi32(src3_16x8b, src7_16x8b);

        src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
        src1_16x8b = _mm_min_epi16(src1_16x8b, clip_max_8x16b);
        src2_16x8b = _mm_min_epi16(src2_16x8b, clip_max_8x16b);
        src3_16x8b = _mm_min_epi16(src3_16x8b, clip_max_8x16b);

        src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);
        src1_16x8b = _mm_max_epi16(src1_16x8b, clip_min_8x16b);
        src2_16x8b = _mm_max_epi16(src2_16x8b, clip_min_8x16b);
        src3_16x8b = _mm_max_epi16(src3_16x8b, clip_min_8x16b);

        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 0 * iDstStride), src0_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 1 * iDstStride), src1_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 2 * iDstStride), src2_16x8b);
        _mm_storeu_si128((__m128i*)(dst_copy + iCount + 3 * iDstStride), src3_16x8b);

        iCount += 8;
      }
      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)	/*loop over height of 1*/
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {
        iCount = 0;
        for (iCols = iBlkWd; iCols > 7; iCols -= 8)
        {
          src0_16x8b = _mm_loadu_si128((__m128i*)(inp_copy + iCount));
          src1_16x8b = _mm_srli_si128(src0_16x8b, 8);
          src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
          src1_16x8b = _mm_cvtepi16_epi32(src1_16x8b);
          src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
          src1_16x8b = _mm_add_epi32(src1_16x8b, add_32x4b);
          src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
          src1_16x8b = _mm_sra_epi32(src1_16x8b, shift_32x4b);
          src0_16x8b = _mm_packs_epi32(src0_16x8b, src1_16x8b);
          src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
          src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);
          _mm_storeu_si128((__m128i*)(dst_copy + iCount), src0_16x8b);
          iCount += 8;
        }
        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }

    }
  }

  inp_copy = pSrc + ((iBlkWd / 8) * 8);
  dst_copy = pDst + ((iBlkWd / 8) * 8);

  if (iRemWd > 3)
  {
    Int iRemHt = iBlkHt % 8;
    iRemWd &= 0x3;
    for (Int j = 0; j < (iBlkHt - 7); j += 8)
    {

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride));
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride));
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride));
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride));
      src4_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 4 * iSrcStride));
      src5_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 5 * iSrcStride));
      src6_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 6 * iSrcStride));
      src7_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 7 * iSrcStride));

      src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
      src1_16x8b = _mm_cvtepi16_epi32(src1_16x8b);
      src2_16x8b = _mm_cvtepi16_epi32(src2_16x8b);
      src3_16x8b = _mm_cvtepi16_epi32(src3_16x8b);
      src4_16x8b = _mm_cvtepi16_epi32(src4_16x8b);
      src5_16x8b = _mm_cvtepi16_epi32(src5_16x8b);
      src6_16x8b = _mm_cvtepi16_epi32(src6_16x8b);
      src7_16x8b = _mm_cvtepi16_epi32(src7_16x8b);

      src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
      src1_16x8b = _mm_add_epi32(src1_16x8b, add_32x4b);
      src2_16x8b = _mm_add_epi32(src2_16x8b, add_32x4b);
      src3_16x8b = _mm_add_epi32(src3_16x8b, add_32x4b);
      src4_16x8b = _mm_add_epi32(src4_16x8b, add_32x4b);
      src5_16x8b = _mm_add_epi32(src5_16x8b, add_32x4b);
      src6_16x8b = _mm_add_epi32(src6_16x8b, add_32x4b);
      src7_16x8b = _mm_add_epi32(src7_16x8b, add_32x4b);

      src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
      src1_16x8b = _mm_sra_epi32(src1_16x8b, shift_32x4b);
      src2_16x8b = _mm_sra_epi32(src2_16x8b, shift_32x4b);
      src3_16x8b = _mm_sra_epi32(src3_16x8b, shift_32x4b);
      src4_16x8b = _mm_sra_epi32(src4_16x8b, shift_32x4b);
      src5_16x8b = _mm_sra_epi32(src5_16x8b, shift_32x4b);
      src6_16x8b = _mm_sra_epi32(src6_16x8b, shift_32x4b);
      src7_16x8b = _mm_sra_epi32(src7_16x8b, shift_32x4b);

      src0_16x8b = _mm_packs_epi32(src0_16x8b, src0_16x8b);
      src1_16x8b = _mm_packs_epi32(src1_16x8b, src1_16x8b);
      src2_16x8b = _mm_packs_epi32(src2_16x8b, src2_16x8b);
      src3_16x8b = _mm_packs_epi32(src3_16x8b, src3_16x8b);
      src4_16x8b = _mm_packs_epi32(src4_16x8b, src4_16x8b);
      src5_16x8b = _mm_packs_epi32(src5_16x8b, src5_16x8b);
      src6_16x8b = _mm_packs_epi32(src6_16x8b, src6_16x8b);
      src7_16x8b = _mm_packs_epi32(src7_16x8b, src7_16x8b);

      src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
      src1_16x8b = _mm_min_epi16(src1_16x8b, clip_max_8x16b);
      src2_16x8b = _mm_min_epi16(src2_16x8b, clip_max_8x16b);
      src3_16x8b = _mm_min_epi16(src3_16x8b, clip_max_8x16b);
      src4_16x8b = _mm_min_epi16(src4_16x8b, clip_max_8x16b);
      src5_16x8b = _mm_min_epi16(src5_16x8b, clip_max_8x16b);
      src6_16x8b = _mm_min_epi16(src6_16x8b, clip_max_8x16b);
      src7_16x8b = _mm_min_epi16(src7_16x8b, clip_max_8x16b);

      src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);
      src1_16x8b = _mm_max_epi16(src1_16x8b, clip_min_8x16b);
      src2_16x8b = _mm_max_epi16(src2_16x8b, clip_min_8x16b);
      src3_16x8b = _mm_max_epi16(src3_16x8b, clip_min_8x16b);
      src4_16x8b = _mm_max_epi16(src4_16x8b, clip_min_8x16b);
      src5_16x8b = _mm_max_epi16(src5_16x8b, clip_min_8x16b);
      src6_16x8b = _mm_max_epi16(src6_16x8b, clip_min_8x16b);
      src7_16x8b = _mm_max_epi16(src7_16x8b, clip_min_8x16b);

      _mm_storel_epi64((__m128i *)(dst_copy + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 3 * iDstStride), src3_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 4 * iDstStride), src4_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 5 * iDstStride), src5_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 6 * iDstStride), src6_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 7 * iDstStride), src7_16x8b);

      inp_copy += (8 * iSrcStride);
      dst_copy += (8 * iDstStride);
    }

    if (iRemHt > 3)
    {
      iRemHt = iRemHt % 4;

      src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 0 * iSrcStride));
      src1_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 1 * iSrcStride));
      src2_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 2 * iSrcStride));
      src3_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy + 3 * iSrcStride));

      src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
      src1_16x8b = _mm_cvtepi16_epi32(src1_16x8b);
      src2_16x8b = _mm_cvtepi16_epi32(src2_16x8b);
      src3_16x8b = _mm_cvtepi16_epi32(src3_16x8b);

      src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
      src1_16x8b = _mm_add_epi32(src1_16x8b, add_32x4b);
      src2_16x8b = _mm_add_epi32(src2_16x8b, add_32x4b);
      src3_16x8b = _mm_add_epi32(src3_16x8b, add_32x4b);

      src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
      src1_16x8b = _mm_sra_epi32(src1_16x8b, shift_32x4b);
      src2_16x8b = _mm_sra_epi32(src2_16x8b, shift_32x4b);
      src3_16x8b = _mm_sra_epi32(src3_16x8b, shift_32x4b);

      src0_16x8b = _mm_packs_epi32(src0_16x8b, src0_16x8b);
      src1_16x8b = _mm_packs_epi32(src1_16x8b, src1_16x8b);
      src2_16x8b = _mm_packs_epi32(src2_16x8b, src2_16x8b);
      src3_16x8b = _mm_packs_epi32(src3_16x8b, src3_16x8b);

      src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
      src1_16x8b = _mm_min_epi16(src1_16x8b, clip_max_8x16b);
      src2_16x8b = _mm_min_epi16(src2_16x8b, clip_max_8x16b);
      src3_16x8b = _mm_min_epi16(src3_16x8b, clip_max_8x16b);

      src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);
      src1_16x8b = _mm_max_epi16(src1_16x8b, clip_min_8x16b);
      src2_16x8b = _mm_max_epi16(src2_16x8b, clip_min_8x16b);
      src3_16x8b = _mm_max_epi16(src3_16x8b, clip_min_8x16b);

      _mm_storel_epi64((__m128i *)(dst_copy + 0 * iDstStride), src0_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 1 * iDstStride), src1_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 2 * iDstStride), src2_16x8b);
      _mm_storel_epi64((__m128i *)(dst_copy + 3 * iDstStride), src3_16x8b);

      inp_copy += (4 * iSrcStride);
      dst_copy += (4 * iDstStride);
    }

    if (iRemHt)
    {
      for (iRows = 0; iRows < iRemHt; iRows++)
      {

        src0_16x8b = _mm_loadl_epi64((__m128i*)(inp_copy));

        src0_16x8b = _mm_cvtepi16_epi32(src0_16x8b);
        src0_16x8b = _mm_add_epi32(src0_16x8b, add_32x4b);
        src0_16x8b = _mm_sra_epi32(src0_16x8b, shift_32x4b);
        src0_16x8b = _mm_packs_epi32(src0_16x8b, src0_16x8b);
        src0_16x8b = _mm_min_epi16(src0_16x8b, clip_max_8x16b);
        src0_16x8b = _mm_max_epi16(src0_16x8b, clip_min_8x16b);

        _mm_storel_epi64((__m128i *)dst_copy, src0_16x8b);

        inp_copy += (iSrcStride);
        dst_copy += (iDstStride);
      }
    }
  }

  inp_copy = pSrc + ((iBlkWd / 4) * 4);
  dst_copy = pDst + ((iBlkWd / 4) * 4);

  if (iRemWd)
  {
    for (Int i = 0; i < iBlkHt; i++)
    {
      for (Int j = 0; j < iRemWd; j++)
      {
        Pel val = *(inp_copy + j);
        val = rightShift_round((val + offset), shift);
#if JVET_D0033_ADAPTIVE_CLIPPING
        val = ClipA(val, compID);
#else
        if (val < minVal)
        {
          val = minVal;
        }
        if (val > maxVal)
        {
          val = maxVal;
        }
#endif
        *(dst_copy + j) = val;
      }
      dst_copy += iDstStride;
      inp_copy += iSrcStride;
    }
  }
}
#endif

