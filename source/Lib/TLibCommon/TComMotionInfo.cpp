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

/** \file     TComMotionInfo.cpp
    \brief    motion information handling classes
*/

#include <memory.h>
#include "TComMotionInfo.h"
#include "assert.h"
#include <stdlib.h>
#include <cmath>
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP || VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
#include "TComDataCU.h"
#include "TComPic.h"
#endif
//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

// --------------------------------------------------------------------------------------------------------------------
// Create / destroy
// --------------------------------------------------------------------------------------------------------------------

#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
Void TComCUMvField::create( UInt uiNumPartition,  UInt uiWidth, UInt uiHeight )
#else
Void TComCUMvField::create( UInt uiNumPartition )
#endif
{
  assert(m_pcMv     == NULL);
  assert(m_pcMvd    == NULL);
  assert(m_piRefIdx == NULL);

  m_pcMv     = new TComMv[ uiNumPartition ];
  m_pcMvd    = new TComMv[ uiNumPartition ];
  m_piRefIdx = new Char  [ uiNumPartition ];

  m_uiNumPartition = uiNumPartition;
#if MRG_CONSIDER_COMPLEX
  // Initialization for motion vector, motion vector difference and reference index
  for ( Int i = 0; i < uiNumPartition; i++ )
  {
    m_pcMv [i].setZero();
    m_pcMvd[i].setZero();
  }
  memset( m_piRefIdx, NOT_VALID, uiNumPartition * sizeof( *m_piRefIdx ) );
#if !JVET_C0024_QTBT
  m_uiWidth  = uiWidth;
  m_uiHeight = uiHeight;
#endif
#endif
}

Void TComCUMvField::destroy()
{
  assert(m_pcMv     != NULL);
  assert(m_pcMvd    != NULL);
  assert(m_piRefIdx != NULL);

  delete[] m_pcMv;
  delete[] m_pcMvd;
  delete[] m_piRefIdx;

  m_pcMv     = NULL;
  m_pcMvd    = NULL;
  m_piRefIdx = NULL;

  m_uiNumPartition = 0;
}

// --------------------------------------------------------------------------------------------------------------------
// Clear / copy
// --------------------------------------------------------------------------------------------------------------------

#if JVET_C0024_QTBT //only for FRUC, which clears the MvField before initCtu.
Void TComCUMvField::clearCtuMvField()
{
  for ( Int i = 0; i < m_uiNumPartition; i++ )
  {
    m_pcMv [ i ].setZero();
    m_pcMvd[ i ].setZero();
  }
  assert( sizeof( *m_piRefIdx ) == 1 );
  memset( m_piRefIdx, NOT_VALID, m_uiNumPartition * sizeof( *m_piRefIdx ) );
}
#endif

Void TComCUMvField::clearMvField()
{
#if JVET_C0024_QTBT
  UInt uiRaster = g_auiZscanToRaster[m_pcCU->getZorderIdxInCtu()];
  UInt uiShort, uiLong;
  UInt uiStride;

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiWidth   = m_pcCU->getWidth(0);
  UInt uiHeight  = m_pcCU->getHeight(0);
  UInt uiNumXDiv = (uiWidth  > m_pcCU->getPic()->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiHeight > m_pcCU->getPic()->getMinCUHeight())? 2 : 1;
  uiWidth  = (uiNumXDiv==1)? uiWidth  : (uiWidth>>1);
  uiHeight = (uiNumYDiv==1)? uiHeight : (uiHeight>>1);
  UInt uiLenX    = uiWidth  / m_pcCU->getPic()->getMinCUWidth();
  UInt uiLenY    = uiHeight / m_pcCU->getPic()->getMinCUHeight();

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      assert(m_pcCU->getUseAbsoluteZIdxToCtu()==false);
      uiRaster = (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcCU->getPic()->getNumPartInCtuWidth();

      if (uiHeight > uiWidth)
      {
        uiShort = uiWidth;
        uiLong  = uiHeight;
        uiStride = m_pcCU->getPic()->getNumPartInCtuWidth();
      }
      else
      {
        uiShort = uiHeight;
        uiLong  = uiWidth;
        uiStride = 1;
      }
#else
  if (m_pcCU->getHeight(0) > m_pcCU->getWidth(0))
  {
    uiShort = m_pcCU->getWidth(0);
    uiLong = m_pcCU->getHeight(0);
    uiStride = m_pcCU->getPic()->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = m_pcCU->getHeight(0);
    uiLong = m_pcCU->getWidth(0);
    uiStride = 1;
  }
#endif

  UInt uiDepth = g_aucConvertToBit[m_pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcCU->getPic()->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcCU->getPic()->getNumPartInCtuWidth() >> uiDepth;

  for (UInt i=0; i<uiLong; i+=uiShort)
  {
#if TT_TREE_DS1 //only called at encoder
    assert(m_pcCU->getUseAbsoluteZIdxToCtu()==false);
    memset(m_piRefIdx + g_auiRasterToZscan[uiRaster], NOT_VALID, uiCurrPartNumb );
    TComMv* pMv = m_pcMv + g_auiRasterToZscan[uiRaster];
    TComMv* pMvd = m_pcMvd + g_auiRasterToZscan[uiRaster];
#else
    memset(m_piRefIdx + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu(), NOT_VALID, uiCurrPartNumb );
    TComMv* pMv = m_pcMv + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu();
    TComMv* pMvd = m_pcMvd + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu();
#endif

    for (UInt j=0; j<uiCurrPartNumb; j++)
    {
      pMv [ j ].setZero();
      pMvd[ j ].setZero();  
    }

    uiRaster += uiNumPartInShort * uiStride;
  }
#if TT_TREE_DS1 //data copy end
  }}
#endif
#else
  for ( Int i = 0; i < m_uiNumPartition; i++ )
  {
    m_pcMv [ i ].setZero();
    m_pcMvd[ i ].setZero();
  }
  assert( sizeof( *m_piRefIdx ) == 1 );
  memset( m_piRefIdx, NOT_VALID, m_uiNumPartition * sizeof( *m_piRefIdx ) );
#endif
}

#if JVET_C0024_QTBT
Void TComCUMvField::copyFromTo( TComCUMvField const * pcCUMvFieldSrc, Int iNumPart, Int iPartAddrSrc, Int iPartAddrDst )
{
  Int iSizeInTComMv = sizeof( TComMv ) * iNumPart;
  
  memcpy( m_pcMv     + iPartAddrDst, pcCUMvFieldSrc->m_pcMv + iPartAddrSrc,     iSizeInTComMv );
  memcpy( m_pcMvd    + iPartAddrDst, pcCUMvFieldSrc->m_pcMvd + iPartAddrSrc,    iSizeInTComMv );
  memcpy( m_piRefIdx + iPartAddrDst, pcCUMvFieldSrc->m_piRefIdx + iPartAddrSrc, sizeof( *m_piRefIdx ) * iNumPart );
}
#endif

Void TComCUMvField::copyFrom( TComCUMvField const * pcCUMvFieldSrc, Int iNumPartSrc, Int iPartAddrDst )
{
  Int iSizeInTComMv = sizeof( TComMv ) * iNumPartSrc;

  memcpy( m_pcMv     + iPartAddrDst, pcCUMvFieldSrc->m_pcMv,     iSizeInTComMv );
  memcpy( m_pcMvd    + iPartAddrDst, pcCUMvFieldSrc->m_pcMvd,    iSizeInTComMv );
  memcpy( m_piRefIdx + iPartAddrDst, pcCUMvFieldSrc->m_piRefIdx, sizeof( *m_piRefIdx ) * iNumPartSrc );
}

Void TComCUMvField::copyTo( TComCUMvField* pcCUMvFieldDst, Int iPartAddrDst ) const
{
  copyTo( pcCUMvFieldDst, iPartAddrDst, 0, m_uiNumPartition );
}

Void TComCUMvField::copyTo( TComCUMvField* pcCUMvFieldDst, Int iPartAddrDst, UInt uiOffset, UInt uiNumPart ) const
{
  Int iSizeInTComMv = sizeof( TComMv ) * uiNumPart;
  Int iOffset = uiOffset + iPartAddrDst;

  memcpy( pcCUMvFieldDst->m_pcMv     + iOffset, m_pcMv     + uiOffset, iSizeInTComMv );
  memcpy( pcCUMvFieldDst->m_pcMvd    + iOffset, m_pcMvd    + uiOffset, iSizeInTComMv );
  memcpy( pcCUMvFieldDst->m_piRefIdx + iOffset, m_piRefIdx + uiOffset, sizeof( *m_piRefIdx ) * uiNumPart );
}

// --------------------------------------------------------------------------------------------------------------------
// Set
// --------------------------------------------------------------------------------------------------------------------

template <typename T>
Void TComCUMvField::setAll( T *p, T const & val, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx  )
{
#if JVET_C0024_QTBT
  assert(eCUMode == SIZE_2Nx2N); 
  UInt uiRaster = g_auiZscanToRaster[m_pcCU->getZorderIdxInCtu()+iPartAddr];
  UInt uiShort, uiLong;
  UInt uiStride;
  UInt uiWidth = m_pcCU->getWidth(iPartAddr);
  UInt uiHeight = m_pcCU->getHeight(iPartAddr);

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiNumXDiv = (uiWidth  > m_pcCU->getPic()->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiHeight > m_pcCU->getPic()->getMinCUHeight())? 2 : 1;
  uiWidth  = (uiNumXDiv==1)? uiWidth  : (uiWidth>>1);
  uiHeight = (uiNumYDiv==1)? uiHeight : (uiHeight>>1);
  UInt uiLenX    = uiWidth  / m_pcCU->getPic()->getMinCUWidth();
  UInt uiLenY    = uiHeight / m_pcCU->getPic()->getMinCUHeight();

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      //assert(m_pcCU->getUseAbsoluteZIdxToCtu()==false);
      if(m_pcCU->getUseAbsoluteZIdxToCtu()==false)
        uiRaster = g_auiZscanToRaster[iPartAddr] + (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcCU->getPic()->getNumPartInCtuWidth();
      else
        uiRaster = g_auiZscanToRaster[m_pcCU->getZorderIdxInCtu()+iPartAddr] + (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcCU->getPic()->getNumPartInCtuWidth();
#endif
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = m_pcCU->getPic()->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  uiDepth = g_aucConvertToBit[m_pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcCU->getPic()->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcCU->getPic()->getNumPartInCtuWidth() >> uiDepth;

  for (UInt i=0; i<uiLong; i+=uiShort)
  {
#if TT_TREE_DS1
    //assert(m_pcCU->getUseAbsoluteZIdxToCtu()==false);
    T* p1;
    if(m_pcCU->getUseAbsoluteZIdxToCtu()==false)
      p1 = p + g_auiRasterToZscan[uiRaster];
    else
      p1 = p + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu();
#else
    T* p1 = p + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu();
#endif
    for (UInt j=0; j<uiCurrPartNumb; j++)
    {
      p1[ j ] = val;
    }
    uiRaster += uiNumPartInShort * uiStride;
  }
#if TT_TREE_DS1 //data copy end
  }}
#endif
#else
  Int i;
  p += iPartAddr;
  Int numElements = m_uiNumPartition >> ( 2 * uiDepth );

  switch( eCUMode )
  {
    case SIZE_2Nx2N:
      for ( i = 0; i < numElements; i++ )
      {
        p[ i ] = val;
      }
      break;

    case SIZE_2NxN:
      numElements >>= 1;
      for ( i = 0; i < numElements; i++ )
      {
        p[ i ] = val;
      }
      break;

    case SIZE_Nx2N:
      numElements >>= 2;
      for ( i = 0; i < numElements; i++ )
      {
        p[ i                   ] = val;
        p[ i + 2 * numElements ] = val;
      }
      break;

    case SIZE_NxN:
      numElements >>= 2;
      for ( i = 0; i < numElements; i++)
      {
        p[ i ] = val;
      }
      break;
    case SIZE_2NxnU:
    {
      Int iCurrPartNumQ = numElements>>2;
      if( iPartIdx == 0 )
      {
        T *pT  = p;
        T *pT2 = p + iCurrPartNumQ;
        for (i = 0; i < (iCurrPartNumQ>>1); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }
      }
      else
      {
        T *pT  = p;
        for (i = 0; i < (iCurrPartNumQ>>1); i++)
        {
          pT[i] = val;
        }

        pT = p + iCurrPartNumQ;
        for (i = 0; i < ( (iCurrPartNumQ>>1) + (iCurrPartNumQ<<1) ); i++)
        {
          pT[i] = val;
        }
      }
      break;
    }
  case SIZE_2NxnD:
    {
      Int iCurrPartNumQ = numElements>>2;
      if( iPartIdx == 0 )
      {
        T *pT  = p;
        for (i = 0; i < ( (iCurrPartNumQ>>1) + (iCurrPartNumQ<<1) ); i++)
        {
          pT[i] = val;
        }
        pT = p + ( numElements - iCurrPartNumQ );
        for (i = 0; i < (iCurrPartNumQ>>1); i++)
        {
          pT[i] = val;
        }
      }
      else
      {
        T *pT  = p;
        T *pT2 = p + iCurrPartNumQ;
        for (i = 0; i < (iCurrPartNumQ>>1); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }
      }
      break;
    }
  case SIZE_nLx2N:
    {
      Int iCurrPartNumQ = numElements>>2;
      if( iPartIdx == 0 )
      {
        T *pT  = p;
        T *pT2 = p + (iCurrPartNumQ<<1);
        T *pT3 = p + (iCurrPartNumQ>>1);
        T *pT4 = p + (iCurrPartNumQ<<1) + (iCurrPartNumQ>>1);

        for (i = 0; i < (iCurrPartNumQ>>2); i++)
        {
          pT [i] = val;
          pT2[i] = val;
          pT3[i] = val;
          pT4[i] = val;
        }
      }
      else
      {
        T *pT  = p;
        T *pT2 = p + (iCurrPartNumQ<<1);
        for (i = 0; i < (iCurrPartNumQ>>2); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }

        pT  = p + (iCurrPartNumQ>>1);
        pT2 = p + (iCurrPartNumQ<<1) + (iCurrPartNumQ>>1);
        for (i = 0; i < ( (iCurrPartNumQ>>2) + iCurrPartNumQ ); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }
      }
      break;
    }
  case SIZE_nRx2N:
    {
      Int iCurrPartNumQ = numElements>>2;
      if( iPartIdx == 0 )
      {
        T *pT  = p;
        T *pT2 = p + (iCurrPartNumQ<<1);
        for (i = 0; i < ( (iCurrPartNumQ>>2) + iCurrPartNumQ ); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }

        pT  = p + iCurrPartNumQ + (iCurrPartNumQ>>1);
        pT2 = p + numElements - iCurrPartNumQ + (iCurrPartNumQ>>1);
        for (i = 0; i < (iCurrPartNumQ>>2); i++)
        {
          pT [i] = val;
          pT2[i] = val;
        }
      }
      else
      {
        T *pT  = p;
        T *pT2 = p + (iCurrPartNumQ>>1);
        T *pT3 = p + (iCurrPartNumQ<<1);
        T *pT4 = p + (iCurrPartNumQ<<1) + (iCurrPartNumQ>>1);
        for (i = 0; i < (iCurrPartNumQ>>2); i++)
        {
          pT [i] = val;
          pT2[i] = val;
          pT3[i] = val;
          pT4[i] = val;
        }
      }
      break;
    }
    default:
      assert(0);
      break;
  }
#endif
}

Void TComCUMvField::setAllMv( TComMv const & mv, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  setAll(m_pcMv, mv, eCUMode, iPartAddr, uiDepth, iPartIdx);
}

Void TComCUMvField::setAllMvd( TComMv const & mvd, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  setAll(m_pcMvd, mvd, eCUMode, iPartAddr, uiDepth, iPartIdx);
}

Void TComCUMvField::setAllRefIdx ( Int iRefIdx, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  setAll(m_piRefIdx, static_cast<Char>(iRefIdx), eCUMode, iPartAddr, uiDepth, iPartIdx);
}

Void TComCUMvField::setAllMvField( TComMvField const & mvField, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  setAllMv    ( mvField.getMv(),     eCUMode, iPartAddr, uiDepth, iPartIdx );
  setAllRefIdx( mvField.getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
}

/**Subsampling of the stored prediction mode, reference index and motion vector
 * \param pePredMode Pointer to prediction modes
 * \param scale      Factor by which to subsample motion information
 */
Void TComCUMvField::compress(Char* pePredMode, Int scale)
{
  Int N = scale * scale;
  assert( N > 0 && N <= m_uiNumPartition);

  for ( Int uiPartIdx = 0; uiPartIdx < m_uiNumPartition; uiPartIdx += N )
  {
    TComMv cMv(0,0);
    Int iRefIdx = 0;

    cMv = m_pcMv[ uiPartIdx ];
    PredMode predMode = static_cast<PredMode>( pePredMode[ uiPartIdx ] );
    iRefIdx = m_piRefIdx[ uiPartIdx ];
    for ( Int i = 0; i < N; i++ )
    {
      m_pcMv[ uiPartIdx + i ] = cMv;
      pePredMode[ uiPartIdx + i ] = predMode;
      m_piRefIdx[ uiPartIdx + i ] = iRefIdx;
    }
  }
}

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP || VCEG_AZ07_FRUC_MERGE
Void TComCUMvField::setMvFieldSP( TComDataCU* pcCU, UInt uiAbsPartIdx, TComMvField cMvField, Int iWidth, Int iHeight  )
{
#if TT_TREE_ATMVP // deriveZIdxInCTU
  uiAbsPartIdx  = pcCU->deriveZIdxInCTU(pcCU->getZorderIdxInCtu(), uiAbsPartIdx);
#else
  uiAbsPartIdx += pcCU->getZorderIdxInCtu();
#endif

  Int iMinCUW  = pcCU->getPic()->getMinCUHeight();
  Int iStartPelX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]]/iMinCUW;
  Int iStartPelY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]]/iMinCUW;
  Int iEndPelX = iStartPelX + iWidth/iMinCUW;
  Int iEndPelY = iStartPelY + iHeight/iMinCUW;

  Int iNumPart = pcCU->getPic()->getNumPartInCtuWidth();
  Int iCurrRaster, uiPartAddr;

  for (Int i=iStartPelY; i < iEndPelY; i ++)
  {
    for (Int j=iStartPelX; j < iEndPelX; j ++)
    {
      iCurrRaster = i * iNumPart + j;
      uiPartAddr  = g_auiRasterToZscan[iCurrRaster];
#if TT_TREE_ATMVP // deriveZIdxInCU
      uiPartAddr  = pcCU->deriveZIdxInCU(pcCU->getZorderIdxInCtu(), uiPartAddr); 
#else
      uiPartAddr -= pcCU->getZorderIdxInCtu();  
#endif

      m_pcMv[uiPartAddr]     = cMvField.getMv();
      m_piRefIdx[uiPartAddr] = cMvField.getRefIdx();
    }
  }
}
#endif

#if MRG_CONSIDER_COMPLEX
Int TComCUMvField::agaus(Double *a,Double *b,Int n)
{
  Int *js, l, k, i, j, is, p, q;
  Double d, t;
  js = new Int[n];
  l = 1;

  for (k=0;k<=n-2;k++)
  {
    d=0.0;
    for (i=k;i<=n-1;i++)
    {
      for (j=k;j<=n-1;j++)
      {
        t=fabs(a[i*n+j]);
        if (t>d)
        {
          d=t;
          js[k]=j;
          is=i;
        }
      }
    }
    if (d+1.0==1.0)
    {
      l=0;
    }
    else
    {
      if (js[k]!=k)
        for (i=0;i<=n-1;i++)
        {
          p=i*n+k;
          q=i*n+js[k];
          t=a[p];
          a[p]=a[q];
          a[q]=t;
        }
        if (is!=k)
        {
          for (j=k;j<=n-1;j++)
          {
            p=k*n+j;
            q=is*n+j;
            t=a[p];
            a[p]=a[q];
            a[q]=t;
          }
          t=b[k];
          b[k]=b[is];
          b[is]=t;
        }
    }
    if (l==0)
    {
      delete[] js;
      return(0);
    }
    d=a[k*n+k];
    for (j=k+1;j<=n-1;j++)
    {
      p=k*n+j;
      a[p]=a[p]/d;
    }
    b[k]=b[k]/d;
    for (i=k+1;i<=n-1;i++)
    {
      for (j=k+1;j<=n-1;j++)
      {
        p=i*n+j;
        a[p]=a[p]-a[i*n+k]*a[k*n+j];
      }
      b[i]=b[i]-a[i*n+k]*b[k];
    }
  } //k

  d=a[(n-1)*n+n-1];
  if (fabs(d)+1.0==1.0)
  {
    delete[] js;
    return(0);
  }
  b[n-1]=b[n-1]/d;

  for (i=n-2;i>=0;i--)
  {
    t=0.0;
    for (j=i+1;j<=n-1;j++)
    {
      t=t+a[i*n+j]*b[j];
    }
    b[i]=b[i]-t;
  }
  js[n-1]=n-1;

  for (k=n-1;k>=0;k--)
  {
    if (js[k]!=k)
    {
      t=b[k];
      b[k]=b[js[k]];
      b[js[k]]=t;
    }
  }
  delete[] js;
  return(1);
}

#if JVET_C0024_QTBT
Void TComCUMvField::setAll( TComMv* pMv, Char* pRefIdx, NodeInfo NodeMotionInfo, UInt MergeComplexIndex,  PartSize eCUMode, RefPicList eRefPicList, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  Int i, j;
  Int l;
  Int iRow;
  Int iCol;

  assert(eCUMode == SIZE_2Nx2N);
#if TT_TREE_AFFINE
  assert(m_pcCU->getZorderIdxInCtu()==0 || iPartAddr==0);
#endif
  UInt uiRaster = g_auiZscanToRaster[m_pcCU->getZorderIdxInCtu()+iPartAddr];
  UInt uiShort, uiLong;
  UInt uiStride;
  Int iPW = m_pcCU->getWidth(iPartAddr);
  Int iPH = m_pcCU->getHeight(iPartAddr);
  if (iPH > iPW)
  {
    uiShort = iPW;
    uiLong  = iPH;
    uiStride = m_pcCU->getPic()->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = iPH;
    uiLong  = iPW;
    uiStride = 1;
  }

  UInt numSubS= uiShort >> 2;
  uiDepth = g_aucConvertToBit[m_pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcCU->getPic()->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcCU->getPic()->getNumPartInCtuWidth() >> uiDepth;

  TComMvField  cMvFieldNeighbourNodes[MRG_NODE_NUM <<1];
  UChar InterDirNeighbourNodes[MRG_NODE_NUM ];
  Int numValidMergeNode[MRG_COMPLEX_NUM_CANDS <<1];
  Int RealComplexIndex[MRG_COMPLEX_NUM_CANDS]={0};

  UInt uiNumPartInWidth = m_pcCU->getPic()->getNumPartInCtuWidth();
  UInt shiftRZ = g_aucConvertToBit[ uiNumPartInWidth ] + 2;

  for (i=0; i<(MRG_NODE_NUM<<1); i++)
  {
    cMvFieldNeighbourNodes[i] = NodeMotionInfo.MvFieldNbNodes[i];
  }
  for (i=0; i < MRG_NODE_NUM; i++)
  {
    InterDirNeighbourNodes[i] = NodeMotionInfo.InterDirNbNodes[i];
  }
  for (i=0; i<(MRG_COMPLEX_NUM_CANDS<<1); i++)
  {
    numValidMergeNode[i] = NodeMotionInfo.numValidNodes[i];
  }
  for (i=0; i<MRG_COMPLEX_NUM_CANDS;i++)
  {
    RealComplexIndex[i] = NodeMotionInfo.RealComplexIdx[i];
  }

  l = (eRefPicList == REF_PIC_LIST_0 ? 0 : 1);

  // mergeComplexIndex > 0, affine model
  Int NodeNum;
  if ( RealComplexIndex[MergeComplexIndex] == 0 )
  {
    NodeNum = 4;
  }
  else if ( RealComplexIndex[MergeComplexIndex] < 5 )
  {
    NodeNum = 3;
  }
  else
  {
    NodeNum = 2;
  }

  Int countMV=1;
  for ( i = 1; i < NodeNum; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if (cMvFieldNeighbourNodes[(i<<1) + l].getMv()==cMvFieldNeighbourNodes[(j<<1) + l].getMv())
        break;
    }
    if ( j == i )
      countMV++;
  }

  if ( countMV == 1 )
  {
    setAllMv( cMvFieldNeighbourNodes[l].getMv(), eCUMode, iPartAddr, uiDepth, iPartIdx );
    setAllRefIdx( cMvFieldNeighbourNodes[l].getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
  }
  else
  {
    // assignment of Reference  Index and motion vector for complex motion,
    setAllRefIdx( cMvFieldNeighbourNodes[l].getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
    
    // when a MVP is derived from a MV pointing to a different reference picture, the MV is scaled to the target reference picture as the final MVP
    TComMv sMV[4];
    for ( j = 0; j < NodeNum; j++ )
    {
      sMV[j] = cMvFieldNeighbourNodes[(j<<1) + l].getMv();
    }

    // calculate the affine coefficients
    Double a[36] = {0};
    Double m[6]  = {0};
    Double b[16] = {0};
    Double n[4]  = {0};
    switch ( RealComplexIndex[MergeComplexIndex] )
    {
    case 1:
      a[0] = (-2);
      a[1] = (-2);
      a[2] = 1;

      a[9]  = (-2);
      a[10] = (-2);
      a[11] = 1;

      a[12] = (iPW + 2);
      a[13] = (-2);
      a[14] = 1;

      a[21] = (-2);
      a[22] = (iPW + 2);
      a[23] = 1;

      a[24] = (-2);
      a[25] = (iPH + 2);
      a[26] = 1;

      a[33] = (iPH + 2);
      a[34] = (-2);
      a[35] = 1;

      m[0] = sMV[0].getHor();
      m[1] = sMV[0].getVer();
      m[2] = sMV[1].getHor();
      m[3] = sMV[1].getVer();
      m[4] = sMV[2].getHor();
      m[5] = sMV[2].getVer();
      agaus(a,m,6);
      break;
    case 2:
      a[0] = (-2);
      a[1] = (-2);
      a[2] = 1;

      a[9]  = (-2);
      a[10] = (-2);
      a[11] = 1;

      a[12] = (iPW + 2);
      a[13] = (-2);
      a[14] = 1;

      a[21] = (-2);
      a[22] = (iPW + 2);
      a[23] = 1;

      a[24] = (iPW + 2);
      a[25] = (iPH + 2);
      a[26] = 1;

      a[33] = (iPH + 2);
      a[34] = (iPW + 2);
      a[35] = 1;

      m[0] = sMV[0].getHor();
      m[1] = sMV[0].getVer();
      m[2] = sMV[1].getHor();
      m[3] = sMV[1].getVer();
      m[4] = sMV[2].getHor();
      m[5] = sMV[2].getVer();
      agaus(a,m,6);
      break;
    case 3:
      a[0] = (-2);
      a[1] = (-2);
      a[2] = 1;

      a[9]  = (-2);
      a[10] = (-2);
      a[11] = 1;

      a[12] = (-2);
      a[13] = (iPH + 2);
      a[14] = 1;

      a[21] = (iPH + 2);
      a[22] = (-2);
      a[23] = 1;

      a[24] = (iPW + 2);
      a[25] = (iPH + 2);
      a[26] = 1;

      a[33] = (iPH + 2);
      a[34] = (iPW + 2);
      a[35] = 1;

      m[0] = sMV[0].getHor();
      m[1] = sMV[0].getVer();
      m[2] = sMV[1].getHor();
      m[3] = sMV[1].getVer();
      m[4] = sMV[2].getHor();
      m[5] = sMV[2].getVer();
      agaus(a,m,6);
      break;
    case 4:
      a[0] = (iPW + 2);
      a[1] = (-2);
      a[2] = 1;

      a[9]  = (-2);
      a[10] = (iPW + 2);
      a[11] = 1;

      a[12] = (-2);
      a[13] = (iPH + 2);
      a[14] = 1;

      a[21] = (iPH + 2);
      a[22] = (-2);
      a[23] = 1;

      a[24] = (iPW + 2);
      a[25] = (iPH + 2);
      a[26] = 1;

      a[33] = (iPH + 2);
      a[34] = (iPW + 2);
      a[35] = 1;

      m[0] = sMV[0].getHor();
      m[1] = sMV[0].getVer();
      m[2] = sMV[1].getHor();
      m[3] = sMV[1].getVer();
      m[4] = sMV[2].getHor();
      m[5] = sMV[2].getVer();
      agaus(a,m,6);
      break;
    case 5:
      b[0] = (-2);
      b[1] = (-2);
      b[2] = 1;

      b[4] = (-2);
      b[5] = 2;
      b[7] = 1;

      b[8]  = (iPW + 2);
      b[9]  = (-2);
      b[10] = 1;

      b[12] = (-2);
      b[13] = -((iPW + 2));
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    case 6:
      b[0] = (-2);
      b[1] = (-2);
      b[2] = 1;

      b[4] = (-2);
      b[5] = 2;
      b[7] = 1;

      b[8]  = (-2);
      b[9]  = (iPH + 2);
      b[10] = 1;

      b[12] = (iPH + 2);
      b[13] =  2 ;
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    case 7:
      b[0] = (iPW + 2);
      b[1] = (-2);
      b[2] = 1;

      b[4] = (-2);
      b[5] = -((iPW + 2));
      b[7] = 1;

      b[8]  = (-2);
      b[9]  = (iPH + 2);
      b[10] = 1;

      b[12] = (iPH + 2);
      b[13] =  2 ;
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    case 8:
      b[0] = (-2);
      b[1] = (-2);
      b[2] = 1;

      b[4] = (-2);
      b[5] = 2;
      b[7] = 1;

      b[8]  = (iPW + 2);
      b[9]  = (iPH + 2);
      b[10] = 1;

      b[12] = (iPH + 2);
      b[13] = -((iPW + 2));
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    case 9:
      b[0] = (iPW + 2);
      b[1] = (-2);
      b[2] = 1;

      b[4] = (-2);
      b[5] = -((iPW + 2));
      b[7] = 1;

      b[8]  = (iPW + 2);
      b[9]  = (iPH + 2);
      b[10] = 1;

      b[12] = (iPH + 2);
      b[13] = -((iPW + 2));
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    case 10:
      b[0] = (-2);
      b[1] = (iPH + 2);
      b[2] = 1;

      b[4] = (iPH + 2);
      b[5] = 2;
      b[7] = 1;

      b[8]  = (iPW + 2);
      b[9]  = (iPH + 2);
      b[10] = 1;

      b[12] = (iPH + 2);
      b[13] = -((iPW + 2));
      b[15] = 1;

      n[0] = sMV[0].getHor();
      n[1] = sMV[0].getVer();
      n[2] = sMV[1].getHor();
      n[3] = sMV[1].getVer();
      agaus(b,n,4);
      break;
    default:
      break;
    }

    Double dHor,dVer;
#if !TT_TREE_AFFINE
    TComMv* pMvT;
#endif
    for (i = 0; i < uiLong; i += uiShort)
    {
#if !TT_TREE_AFFINE
      pMvT = pMv + g_auiRasterToZscan[uiRaster] - m_pcCU->getZorderIdxInCtu();
#endif
      for (UInt j = 0; j < uiCurrPartNumb; j++)
      {
        iRow = (uiShort==iPW) ? (j / numSubS) * 4 + 2 + i : (j / numSubS) * 4 + 2; // y
        iCol = (uiShort==iPH) ? (j % numSubS) * 4 + 2 + i : (j % numSubS) * 4 + 2; // x
        if ( RealComplexIndex[MergeComplexIndex] == 0 ) // NodeNum == 4
        {
          Double weight[4]; // Weighting factor
          weight[0] = (Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4) * (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4) * (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4) * (iPH + 4));
          weight[3] = (Double) ((iRow + 2) * (iCol + 2 ))/((iPW + 4) * (iPH + 4));

          dHor = ( weight[0] * sMV[0].getHor() +weight[1] * sMV[1].getHor() + weight[2] * sMV[2].getHor() +weight[3] * sMV[3].getHor() );
          dVer = ( weight[0] * sMV[0].getVer() +weight[1] * sMV[1].getVer() + weight[2] * sMV[2].getVer() +weight[3] * sMV[3].getVer() );
        }
        else
        {
          if ( NodeNum == 3 ) // 6 parameter affine model
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else // 4 parameter affine model
          {
            dHor = (Double)(  n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }

        Int iSign = dHor < 0 ? -1 : 1;
#if TT_TREE_AFFINE
        assert(m_pcCU->getUseAbsoluteZIdxToCtu()==false);
        UInt uiSubPartRaster = ((j / numSubS) * m_pcCU->getPic()->getNumPartInCtuWidth()) + (j % numSubS) + uiRaster; // this is CTUZIdx
        UInt uiSubPartZIdx   = g_auiRasterToZscan[uiSubPartRaster];
        //Note: when m_pcCU->getZorderIdxInCtu()!=0, it is called in encoder and m_pcCU is a CU not LCU
        UInt uiOffset = m_pcCU->getZorderIdxInCtu()==0? uiSubPartZIdx : m_pcCU->deriveZIdxInCUNoCheck(m_pcCU->getZorderIdxInCtu(),uiSubPartZIdx);
        pMv[ uiOffset ].setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );
#else
        pMvT[ g_auiRasterToZscan[ ((j / numSubS)<<shiftRZ) + (j % numSubS) ] ].setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );
#endif

        iSign = dVer < 0 ? -1 : 1;
#if TT_TREE_AFFINE
        pMv[ uiOffset ].setVer( iSign * (Short)floor( abs(dVer) + 0.5 ) );
#else
        pMvT[ g_auiRasterToZscan[ ((j / numSubS)<<shiftRZ) + (j % numSubS) ] ].setVer( iSign * (Short)floor( abs(dVer) + 0.5 ) );
#endif
      }
      uiRaster += uiNumPartInShort * uiStride;
    }
  }
}
#else
Void TComCUMvField::setAll( TComMv* pMv, Char* pRefIdx, NodeInfo NodeMotionInfo, UInt MergeComplexIndex, PartSize eCUMode, RefPicList eRefPicList, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  Int i, j;
  Int l;
  Int iRow;
  Int iCol;
  pMv += iPartAddr;
  pRefIdx += iPartAddr;
  Int numElements = m_uiNumPartition >> (  uiDepth<<1 );
  Int iWidth = (m_uiWidth>>uiDepth)>>2; //CU Width mesured With SCU
  Int iHeight =(m_uiHeight>>uiDepth)>>2;//CU Height mesured With SCU
  Int iPW = (m_uiWidth>>uiDepth); //CU Width mesured With Pixel
  Int iPH = (m_uiHeight>>uiDepth);//CU Height mesured With Pixel
  Double weight[4];//Weighting factor
  TComMvField  cMvFieldNeighbourNodes[MRG_NODE_NUM <<1];
  UChar InterDirNeighbourNodes[MRG_NODE_NUM ];
  Int numValidMergeNode[MRG_COMPLEX_NUM_CANDS <<1];
  Int RealComplexIndex[MRG_COMPLEX_NUM_CANDS]={0};
  UInt  uiMinCUWidth  = g_uiMaxCUWidth  >> g_uiMaxCUDepth;
  UInt  uiNumPartInWidth  = g_uiMaxCUWidth  / uiMinCUWidth;
  UInt shiftRZ = g_aucConvertToBit[ uiNumPartInWidth ] + 2;

  for(Int i=0; i<(MRG_NODE_NUM <<1); i++)
  {
    cMvFieldNeighbourNodes[i] = NodeMotionInfo.MvFieldNbNodes[i];
  }
  for(Int i=0; i < MRG_NODE_NUM; i++)
  {
    InterDirNeighbourNodes[i] = NodeMotionInfo.InterDirNbNodes[i];
  }
  for(Int i=0; i<(MRG_COMPLEX_NUM_CANDS <<1); i++)
  {
    numValidMergeNode[i] = NodeMotionInfo.numValidNodes[i];
  }
  for(Int i=0; i<MRG_COMPLEX_NUM_CANDS;i++)
  {
    RealComplexIndex[i] = NodeMotionInfo.RealComplexIdx[i];
  }
  //assert ( eCUMode == SIZE_2Nx2N );
  //assinment of motion vector for complex motion
  l = (eRefPicList == REF_PIC_LIST_0 ? 0 : 1);

  // mergeComplexIndex > 0, affine model
  Int NodeNum;
  if(RealComplexIndex[MergeComplexIndex]==0)
  {
    NodeNum = 4;
  }
  else if(RealComplexIndex[MergeComplexIndex] < 5)
  {
    NodeNum = 3;
  }
  else
  {
    NodeNum = 2;
  }

  if((eCUMode==SIZE_2NxnU) &&(iPartIdx ==1))
  {
    iHeight = (iHeight>>2) + (iHeight>>1);
    iPH =  (iPH>>2) + (iPH>>1);
    Int iCurrPartNumQ = numElements>>2;
    Int countMV=1;
    for(i = 1; i < NodeNum; i ++)
    {
      for(j = 0; j < i; j ++)
      {
        if(cMvFieldNeighbourNodes[(i<<1)  + l]. getMv()==cMvFieldNeighbourNodes[(j<<1)  + l]. getMv())
          break;
      }
      if(j==i)
        countMV++;
    }
    if(countMV==1)
    {
      setAllMv( cMvFieldNeighbourNodes[ l ]. getMv(), eCUMode, iPartAddr, uiDepth, iPartIdx );
      setAllRefIdx( cMvFieldNeighbourNodes[ l ]. getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
    }
    else
    {

      setAllRefIdx( cMvFieldNeighbourNodes[l ]. getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
      //when a MVP is derived from a MV pointing to a different reference picture, the MV is scaled to the target reference picture as the final MVP
      TComMv sMV[4];

      for(j = 0; j < NodeNum; j++)
      {
        sMV[ j ]=  cMvFieldNeighbourNodes[(j<<1)  + l].getMv();
      }

      //calculate the affine coefficients
      Double a[36] = {0};
      Double m[6] = {0};
      Double b[16] = {0};
      Double n[4] = {0};
      switch(RealComplexIndex[MergeComplexIndex])
      {
      case 1:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(-2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (-2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 2:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 3:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 4:
        a[0] = (iPW + 2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (iPW + 2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 5:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (-2);
        b[10] = 1;
        b[12] = (-2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 6:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 7:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 8:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 9:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 10:
        b[0] = (-2);
        b[1] = (iPH + 2);
        b[2] = 1;
        b[4] = (iPH + 2);
        b[5] =  2 ;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      default:
        break;

      }
      Double dHor,dVer;
      for ( i = 0; i < iCurrPartNumQ; i++ )
      {
        iRow = (i / iWidth)*4 + 2;// pixel row
        iCol  = (i % iWidth)*4 + 2; //pixel col
        if(RealComplexIndex[MergeComplexIndex]==0)
        {
          //weighting factor
          weight[0] =(Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[3] =(Double)  ((iRow + 2) * (iCol + 2 ))/((iPW + 4)* (iPH + 4));

          //interpolate the motion vector
          dHor=(weight[0] * sMV[ 0 ].getHor() +weight[1] * sMV[ 1 ].getHor() + weight[2] *sMV[ 2 ].getHor() +weight[3] * sMV[ 3 ].getHor() );
          dVer = (weight[0] * sMV[ 0 ].getVer() +weight[1] * sMV[ 1 ].getVer() + weight[2] * sMV[ 2 ].getVer() +weight[3] * sMV[ 3 ].getVer());
        }
        else
        {
          //derive the motion factor
          if(NodeNum==3)
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else
          {
            dHor = (Double)(n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }
        Int iSign = dHor < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );
        iSign = dVer < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setVer(iSign * (Short)floor( abs(dVer) + 0.5 ));
      }

      pMv +=  iCurrPartNumQ + (iCurrPartNumQ>>1);
      for( i= 0; i < (iCurrPartNumQ<<1); i++)
      {
        iRow = ((i + iCurrPartNumQ)/ iWidth)*4 + 2; //pixel row
        iCol  =((i + iCurrPartNumQ) % iWidth)*4 + 2;//pixel height
        if(RealComplexIndex[MergeComplexIndex]==0)
        {
          //weighting factor
          weight[0] =(Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[3] =(Double)  ((iRow + 2) * (iCol + 2 ))/((iPW + 4)* (iPH + 4));

          //interpolate the motion vector
          dHor=(weight[0] * sMV[ 0 ].getHor() +weight[1] * sMV[ 1 ].getHor() + weight[2] *sMV[ 2 ].getHor() +weight[3] * sMV[ 3 ].getHor() );
          dVer = (weight[0] * sMV[ 0 ].getVer() +weight[1] * sMV[ 1 ].getVer() + weight[2] * sMV[ 2 ].getVer() +weight[3] * sMV[ 3 ].getVer());
        }
        else
        {
          if(NodeNum==3)
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else
          {
            dHor = (Double)(n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }
        Int iSign = dHor < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );

        iSign = dVer < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setVer(iSign * (Short)floor( abs(dVer) + 0.5 ));
      }

    }
  }
  else if((eCUMode ==SIZE_nLx2N)&&(iPartIdx ==1))
  {
    Int iCurrPartNumQ = numElements>>2;
    Int iWidth1 = iWidth>>2;
    Int iWidth2 = iWidth>>1;
    iWidth = iWidth1 + iWidth2;
    iPW = (iPW>>2) + (iPW>>1);
    Int countMV=1;
    for(i = 1; i < NodeNum; i ++)
    {
      for(j = 0; j < i; j ++)
      {
        if(cMvFieldNeighbourNodes[(i<<1)  + l]. getMv()==cMvFieldNeighbourNodes[(j<<1)  + l]. getMv())
          break;
      }
      if(j==i)
        countMV++;
    }
    if(countMV==1)
    {
      setAllMv( cMvFieldNeighbourNodes[l ].getMv(), eCUMode, iPartAddr, uiDepth, iPartIdx );
      setAllRefIdx( cMvFieldNeighbourNodes[l ].getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
    }
    else
    {
      //assinment of Reference  Index and motion vector for complex motion,
      setAllRefIdx( cMvFieldNeighbourNodes[l].getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );

      //when a MVP is derived from a MV pointing to a different reference picture, the MV is scaled to the target reference picture as the final MVP
      TComMv sMV[4];

      for(j = 0; j < NodeNum; j++)
      {
        sMV[ j ]=  cMvFieldNeighbourNodes[(j<<1)  + l].getMv();
      }

      //calculate the projection coefficients

      //calculate the affine coefficients
      Double a[36] = {0};
      Double m[6] = {0};
      Double b[16] = {0};
      Double n[4] = {0};
      switch(RealComplexIndex[MergeComplexIndex])
      {
      case 1:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(-2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (-2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 2:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 3:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 4:
        a[0] = (iPW + 2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (iPW + 2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 5:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (-2);
        b[10] = 1;
        b[12] = (-2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 6:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 7:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 8:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 9:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 10:
        b[0] = (-2);
        b[1] = (iPH + 2);
        b[2] = 1;
        b[4] = (iPH + 2);
        b[5] =  2 ;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      default:
        break;

      }
      Double dHor,dVer;
      for ( i = 0; i < iCurrPartNumQ; i++ )
      {
        iRow = (i / iWidth1)*4 + 2;
        iCol  = (i % iWidth1)*4 + 2;
        if(RealComplexIndex[MergeComplexIndex]==0)
        {
          weight[0] =(Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[3] =(Double)  ((iRow + 2) * (iCol + 2 ))/((iPW + 4)* (iPH + 4));

          dHor=(weight[0] * sMV[ 0 ].getHor() +weight[1] * sMV[ 1 ].getHor() + weight[2] *sMV[ 2 ].getHor() +weight[3] * sMV[ 3 ].getHor() );
          dVer = (weight[0] * sMV[ 0 ].getVer() +weight[1] * sMV[ 1 ].getVer() + weight[2] * sMV[ 2 ].getVer() +weight[3] * sMV[ 3 ].getVer());
        }
        else
        {

          //derive the motion factor

          if(NodeNum==3)
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else
          {
            dHor = (Double)(n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }
        Int iSign = dHor < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth1)<<shiftRZ) + (i % iWidth1)]] .setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );

        iSign = dVer < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth1)<<shiftRZ) + (i % iWidth1)]] .setVer(iSign * (Short)floor( abs(dVer) + 0.5 ));
      }

      pMv +=  (iCurrPartNumQ>>2) + (iCurrPartNumQ>>1);

      for ( i = 0; i < (iCurrPartNumQ<<1); i++ )
      {
        iRow = (i / iWidth2)*4 + 2;
        iCol  = ((i % iWidth2)+ iWidth1)*4 + 2;
        if(RealComplexIndex[MergeComplexIndex]==0)
        {

          weight[0] =(Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[3] =(Double)  ((iRow + 2) * (iCol + 2 ))/((iPW + 4)* (iPH + 4));

          dHor=(weight[0] * sMV[ 0 ].getHor() +weight[1] * sMV[ 1 ].getHor() + weight[2] *sMV[ 2 ].getHor() +weight[3] * sMV[ 3 ].getHor() );
          dVer = (weight[0] * sMV[ 0 ].getVer() +weight[1] * sMV[ 1 ].getVer() + weight[2] * sMV[ 2 ].getVer() +weight[3] * sMV[ 3 ].getVer());
        }
        else
        {
          //derive the motion factor
          if(NodeNum==3)
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else
          {
            dHor = (Double)(n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }
        Int iSign = dHor < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth2)<<shiftRZ) + (i % iWidth2)]] .setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );

        iSign = dVer < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth2)<<shiftRZ) + (i % iWidth2)]] .setVer(iSign * (Short)floor( abs(dVer) + 0.5 ));
      }
    }
  }
  else
  {
    switch ( eCUMode )
    {
    case SIZE_2NxN:
      iHeight  >>= 1; iPH  >>= 1;     numElements >>= 1;
      break;
    case SIZE_Nx2N:
      iWidth  >>= 1;  iPW  >>= 1;   numElements >>= 1;
      break;
    case SIZE_NxN:
      iHeight  >>= 1;       iWidth  >>= 1;  iPH  >>= 1;       iPW  >>= 1;    numElements >>= 2;
      break;
    case SIZE_2NxnU:
      iHeight = ( iPartIdx == 0 ) ?  iHeight >> 2 : ( iHeight >> 2 ) + ( iHeight >> 1 );
      iPH = ( iPartIdx == 0 ) ?  iPH >> 2 : ( iPH >> 2 ) + ( iPH >> 1 );
      numElements = ( iPartIdx == 0 ) ? numElements >> 2: ( numElements >> 2 ) + ( numElements >> 1 );
      break;
    case SIZE_2NxnD:
      iHeight = ( iPartIdx == 0 ) ?  ( iHeight >> 2 ) + ( iHeight >> 1 ) : iHeight >> 2;
      iPH = ( iPartIdx == 0 ) ?  ( iPH >> 2 ) + ( iPH >> 1 ) : iPH >> 2;
      numElements = ( iPartIdx == 0 ) ? ( numElements >> 2 ) + ( numElements >> 1 ) : numElements >> 2;
      break;
    case SIZE_nLx2N:
      iWidth = ( iPartIdx == 0 ) ? iWidth >> 2 : ( iWidth >> 2 ) + ( iWidth >> 1 );
      iPW = ( iPartIdx == 0 ) ? iPW >> 2 : ( iPW >> 2 ) + ( iPW >> 1 );
      numElements = ( iPartIdx == 0 ) ? numElements >> 2: ( numElements >> 2 ) + ( numElements >> 1 );
      break;
    case SIZE_nRx2N:
      iWidth     = ( iPartIdx == 0 ) ? ( iWidth >> 2 ) + ( iWidth >> 1 ) : iWidth >> 2;
      iPW     = ( iPartIdx == 0 ) ? ( iPW >> 2 ) + ( iPW >> 1 ) : iPW >> 2;
      numElements = ( iPartIdx == 0 ) ? ( numElements >> 2 ) + ( numElements >> 1 ) : numElements >> 2;
      break;
    default:
      assert (  eCUMode == SIZE_2Nx2N );
      break;
    }
    Int countMV=1;
    for(i = 1; i < NodeNum; i ++)
    {
      for(j = 0; j < i; j ++)
      {
        if(cMvFieldNeighbourNodes[(i<<1)  + l]. getMv()==cMvFieldNeighbourNodes[(j<<1)  + l]. getMv())
          break;
      }
      if(j==i)
        countMV++;
    }
    if(countMV==1)
    {
      setAllMv( cMvFieldNeighbourNodes[l]. getMv(), eCUMode, iPartAddr, uiDepth, iPartIdx );
      setAllRefIdx( cMvFieldNeighbourNodes[l]. getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
    }
    else
    {
      //assinment of Reference  Index and motion vector for complex motion,

      setAllRefIdx( cMvFieldNeighbourNodes[l ]. getRefIdx(), eCUMode, iPartAddr, uiDepth, iPartIdx );
      //when a MVP is derived from a MV pointing to a different reference picture, the MV is scaled to the target reference picture as the final MVP
      TComMv sMV[4];

      for(j = 0; j < NodeNum; j++)
      {
        sMV[ j ]=  cMvFieldNeighbourNodes[(j<<1)  + l].getMv();
      }

      //calculate the affine coefficients
      Double a[36] = {0};
      Double m[6] = {0};
      Double b[16] = {0};
      Double n[4] = {0};
      switch(RealComplexIndex[MergeComplexIndex])
      {
      case 1:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(-2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (-2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 2:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (iPW + 2);
        a[13] = (-2);
        a[14] = 1;
        a[21]=(-2);
        a[22]=(iPW + 2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 3:
        a[0] = (-2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (-2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 4:
        a[0] = (iPW + 2);
        a[1] = (-2);
        a[2] = 1;
        a[9] = (-2);
        a[10] = (iPW + 2);
        a[11] = 1;
        a[12] = (-2);
        a[13] = (iPH + 2);
        a[14] = 1;
        a[21]=(iPH + 2);
        a[22]=(-2);
        a[23]=1;
        a[24]=(iPW + 2);
        a[25]=(iPH + 2);
        a[26]=1;
        a[33] = (iPH + 2);
        a[34] = (iPW + 2);
        a[35] = 1;
        m[0]=sMV[ 0 ].getHor();
        m[1]=sMV[ 0 ].getVer();
        m[2]=sMV[ 1 ].getHor();
        m[3]=sMV[ 1 ].getVer();
        m[4]=sMV[ 2 ].getHor();
        m[5]=sMV[ 2 ].getVer();
        agaus(a,m,6);
        break;
      case 5:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (-2);
        b[10] = 1;
        b[12] = (-2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 6:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 7:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (-2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] =  2 ;
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 8:
        b[0] = (-2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = 2;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 9:
        b[0] = (iPW + 2);
        b[1] = (-2);
        b[2] = 1;
        b[4] = (-2);
        b[5] = -((iPW + 2));
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      case 10:
        b[0] = (-2);
        b[1] = (iPH + 2);
        b[2] = 1;
        b[4] = (iPH + 2);
        b[5] =  2 ;
        b[7] = 1;
        b[8] = (iPW + 2);
        b[9] = (iPH + 2);
        b[10] = 1;
        b[12] = (iPH + 2);
        b[13] = -((iPW + 2));
        b[15] = 1;
        n[0]=sMV[ 0 ].getHor();
        n[1]=sMV[ 0 ].getVer();
        n[2]=sMV[ 1 ].getHor();
        n[3]=sMV[ 1 ].getVer();
        agaus(b,n,4);
        break;
      default:
        break;

      }
      Double dHor,dVer;
      for ( i = 0; i < numElements; i++ )
      {
        iRow = (i / iWidth)*4 + 2;
        iCol  = (i % iWidth)*4 + 2;
        if(RealComplexIndex[MergeComplexIndex]==0)
        {
          weight[0] =(Double) ((iPH - iRow + 2 ) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[1] = (Double) ((iPH - iRow + 2 ) * (iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[2] = (Double) ((iRow + 2) * (iPW - iCol + 2))/((iPW + 4)* (iPH + 4));
          weight[3] =(Double)  ((iRow + 2) * (iCol + 2 ))/((iPW + 4)* (iPH + 4));

          dHor=(weight[0] * sMV[ 0 ].getHor() +weight[1] * sMV[ 1 ].getHor() + weight[2] *sMV[ 2 ].getHor() +weight[3] * sMV[ 3 ].getHor() );
          dVer = (weight[0] * sMV[ 0 ].getVer() +weight[1] * sMV[ 1 ].getVer() + weight[2] * sMV[ 2 ].getVer() +weight[3] * sMV[ 3 ].getVer());
        }
        else
        {
          if(NodeNum==3)
          {
            dHor = (Double)(m[0]*iCol + m[1]*iRow + m[2]);
            dVer = (Double)(m[3]*iRow + m[4]*iCol + m[5]);
          }
          else
          {
            dHor = (Double)(n[0]*iCol + n[1]*iRow + n[2]);
            dVer = (Double)(-(n[1]*iCol)+ n[0]*iRow + n[3]);
          }
        }

        Int iSign = dHor < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setHor( iSign * (Short)floor( abs(dHor) + 0.5 ) );

        iSign = dVer < 0 ? -1 : 1;
        pMv[ g_auiRasterToZscan [((i / iWidth)<<shiftRZ) + (i % iWidth)]] .setVer(iSign * (Short)floor( abs(dVer) + 0.5 ));
      }
    }
  }
}
#endif

Void TComCUMvField::setAllMvFieldNode( NodeInfo NodeMotionInfo,UInt MergeComplexIndex, PartSize eCUMode, RefPicList eRefPicList, Int iPartAddr, UInt uiDepth, Int iPartIdx )
{
  setAll(m_pcMv, m_piRefIdx, NodeMotionInfo, MergeComplexIndex, eCUMode, eRefPicList, iPartAddr, uiDepth, iPartIdx);
}
#endif

//! \}
