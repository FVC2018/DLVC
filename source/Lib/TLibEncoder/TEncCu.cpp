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

/** \file     TEncCu.cpp
    \brief    Coding Unit (CU) encoder class
*/

#include <stdio.h>
#include "TEncTop.h"
#include "TEncCu.h"
#include "TEncAnalyze.h"
#include "TLibCommon/Debug.h"

#include <cmath>
#include <algorithm>
using namespace std;

//! \ingroup TLibEncoder
//! \{
#if VCEG_AZ08_INTER_KLT
extern Bool g_bEnableCheck;
#endif
// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

/**
 \param    uhTotalDepth  total number of allowable depth
 \param    uiMaxWidth    largest CU width
 \param    uiMaxHeight   largest CU height
 \param    chromaFormat  chroma format
 */
Void TEncCu::create(UChar uhTotalDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormat)
{
#if !JVET_C0024_QTBT
  Int i;
#endif

  m_uhTotalDepth   = uhTotalDepth + 1;
#if JVET_C0024_QTBT
  UInt uiNumPartitions = 1<<( ( m_uhTotalDepth - 1 )<<1 );
  assert(uiNumPartitions == 1<<(g_aucConvertToBit[uiMaxWidth] + g_aucConvertToBit[uiMaxHeight]));
  assert(uiMaxWidth>>(m_uhTotalDepth-1) == (1<<MIN_CU_LOG2));

  UInt uiNumWidthIdx = g_aucConvertToBit[uiMaxWidth] + 1;
  UInt uiNumHeightIdx = g_aucConvertToBit[uiMaxHeight] + 1;

  m_pppcBestCU      = new TComDataCU**[uiNumWidthIdx];
  m_pppcTempCU      = new TComDataCU**[uiNumWidthIdx];
#if COM16_C806_OBMC
  m_pppcTempCUWoOBMC  = new TComDataCU**[uiNumWidthIdx];
  m_pppcTmpYuv1       = new TComYuv**[uiNumWidthIdx];
  m_pppcTmpYuv2       = new TComYuv**[uiNumWidthIdx];
  m_pppcPredYuvWoOBMC = new TComYuv**[uiNumWidthIdx];
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_pppcFRUCBufferCU  = new TComDataCU**[uiNumWidthIdx];
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    m_pppcTempCUIMVCache[size] = new TComDataCU**[uiNumWidthIdx];
  }
#endif
  m_pppcPredYuvBest = new TComYuv**[uiNumWidthIdx];
  m_pppcResiYuvBest = new TComYuv**[uiNumWidthIdx];
  m_pppcRecoYuvBest = new TComYuv**[uiNumWidthIdx];
  m_pppcPredYuvTemp = new TComYuv**[uiNumWidthIdx];
  m_pppcResiYuvTemp = new TComYuv**[uiNumWidthIdx];
  m_pppcRecoYuvTemp = new TComYuv**[uiNumWidthIdx];
  m_pppcOrigYuv     = new TComYuv**[uiNumWidthIdx];

#if COM16_C806_LARGE_CTU
  for( UInt i = 0 ; i < NUMBER_OF_STORED_RESIDUAL_TYPES ; i++ )
  {
    m_resiBuffer[i] = new Pel [MAX_CU_SIZE * MAX_CU_SIZE];
  }
#endif

#if BARC_FUNC
#if BARC_ONLY_FOR_LUMA
  m_reconLumaBufferHM = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#if COMPACT_RESOLUTION_CNN
  m_reconBufferBicubic = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
  m_DScaleBufferBicubic = new Pel[(MAX_CU_SIZE >> BARC_DScaleRatio)*(MAX_CU_SIZE >> BARC_DScaleRatio)];
#endif
#else
  for (UInt i = 0; i < MAX_NUM_COMPONENT; i++)
  {
    m_reconBuffer[i] = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#if COMPACT_RESOLUTION_CNN
	m_reconBufferBicubic[i] = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#endif
  }
  m_reconLumaBuffer = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#endif
#endif

  for (UInt wIdx=0; wIdx<uiNumWidthIdx; wIdx++)
  {
    m_pppcBestCU[wIdx] = new TComDataCU*[uiNumHeightIdx];
    m_pppcTempCU[wIdx] = new TComDataCU*[uiNumHeightIdx];

#if COM16_C806_OBMC
    m_pppcTempCUWoOBMC[wIdx]  = new TComDataCU*[uiNumHeightIdx];
    m_pppcTmpYuv1[wIdx]       = new TComYuv*[uiNumHeightIdx];
    m_pppcTmpYuv2[wIdx]       = new TComYuv*[uiNumHeightIdx];
    m_pppcPredYuvWoOBMC[wIdx] = new TComYuv*[uiNumHeightIdx];
#endif
#if VCEG_AZ07_FRUC_MERGE
    m_pppcFRUCBufferCU[wIdx]  = new TComDataCU*[uiNumHeightIdx];
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
    for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
    {
      m_pppcTempCUIMVCache[size][wIdx] = new TComDataCU*[uiNumHeightIdx];
    }
#endif
    m_pppcPredYuvBest[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcResiYuvBest[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcRecoYuvBest[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcPredYuvTemp[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcResiYuvTemp[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcRecoYuvTemp[wIdx] = new TComYuv*[uiNumHeightIdx];
    m_pppcOrigYuv[wIdx]     = new TComYuv*[uiNumHeightIdx];

    for (UInt hIdx=0; hIdx<uiNumHeightIdx; hIdx++)
    {
      UInt uiWidth = 1<<(wIdx+MIN_CU_LOG2);
      UInt uiHeight = 1<<(hIdx+MIN_CU_LOG2);

      m_pppcBestCU[wIdx][hIdx] = new TComDataCU; m_pppcBestCU[wIdx][hIdx]->create(chromaFormat, uiNumPartitions, uiMaxWidth, uiMaxHeight, false, 1<<MIN_CU_LOG2, uiWidth, uiHeight);
      m_pppcTempCU[wIdx][hIdx] = new TComDataCU; m_pppcTempCU[wIdx][hIdx]->create(chromaFormat, uiNumPartitions, uiMaxWidth, uiMaxHeight, false, 1<<MIN_CU_LOG2, uiWidth, uiHeight);

#if COM16_C806_OBMC
      m_pppcTempCUWoOBMC[wIdx][hIdx]  = new TComDataCU; m_pppcTempCUWoOBMC[wIdx][hIdx]->create(chromaFormat, uiNumPartitions, uiMaxWidth, uiMaxHeight, false, 1<<MIN_CU_LOG2, uiWidth, uiHeight);
      m_pppcTmpYuv1[wIdx][hIdx]       = new TComYuv; m_pppcTmpYuv1[wIdx][hIdx]->create(uiWidth+4, uiHeight+4, chromaFormat );
      m_pppcTmpYuv2[wIdx][hIdx]       = new TComYuv; m_pppcTmpYuv2[wIdx][hIdx]->create(uiWidth+4, uiHeight+4, chromaFormat );
      m_pppcPredYuvWoOBMC[wIdx][hIdx] = new TComYuv; m_pppcPredYuvWoOBMC[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
#endif
#if VCEG_AZ07_FRUC_MERGE
      m_pppcFRUCBufferCU[wIdx][hIdx]  = new TComDataCU; m_pppcFRUCBufferCU[wIdx][hIdx]->create(chromaFormat, uiNumPartitions, uiMaxWidth, uiMaxHeight, false, 1<<MIN_CU_LOG2, uiWidth, uiHeight);
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
      for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
      {
        m_pppcTempCUIMVCache[size][wIdx][hIdx] = new TComDataCU; m_pppcTempCUIMVCache[size][wIdx][hIdx]->create(chromaFormat, uiNumPartitions, uiMaxWidth, uiMaxHeight, false, 1<<MIN_CU_LOG2, uiWidth, uiHeight);
      }
#endif
      m_pppcPredYuvBest[wIdx][hIdx] = new TComYuv; m_pppcPredYuvBest[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcResiYuvBest[wIdx][hIdx] = new TComYuv; m_pppcResiYuvBest[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcRecoYuvBest[wIdx][hIdx] = new TComYuv; m_pppcRecoYuvBest[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcPredYuvTemp[wIdx][hIdx] = new TComYuv; m_pppcPredYuvTemp[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcResiYuvTemp[wIdx][hIdx] = new TComYuv; m_pppcResiYuvTemp[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcRecoYuvTemp[wIdx][hIdx] = new TComYuv; m_pppcRecoYuvTemp[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
      m_pppcOrigYuv[wIdx][hIdx]     = new TComYuv; m_pppcOrigYuv[wIdx][hIdx]->create( uiWidth, uiHeight, chromaFormat );
    }
  }
#if JVET_C0024_FAST_MRG
  for (UInt i=0; i<MRG_MAX_NUM_CANDS; i++)
  {
    m_pcMrgPredTempYuv[i] = new TComYuv;  m_pcMrgPredTempYuv[i]->create(uiMaxWidth, uiMaxHeight, chromaFormat);
  }
#endif
#else
  m_ppcBestCU      = new TComDataCU*[m_uhTotalDepth-1];
  m_ppcTempCU      = new TComDataCU*[m_uhTotalDepth-1];
#if COM16_C806_OBMC
  m_ppcTempCUWoOBMC  = new TComDataCU*[m_uhTotalDepth-1];
  m_ppcTmpYuv1       = new TComYuv*[m_uhTotalDepth-1];
  m_ppcTmpYuv2       = new TComYuv*[m_uhTotalDepth-1];
  m_ppcPredYuvWoOBMC = new TComYuv*[m_uhTotalDepth-1];
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_ppcFRUCBufferCU  = new TComDataCU*[m_uhTotalDepth-1];
#endif
#if VCEG_AZ07_IMV
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    m_ppcTempCUIMVCache[size] = new TComDataCU*[m_uhTotalDepth-1];
  }
#endif
m_ppcPredYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcPredYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcOrigYuv     = new TComYuv*[m_uhTotalDepth-1];

#if COM16_C806_LARGE_CTU
  for( i = 0 ; i < NUMBER_OF_STORED_RESIDUAL_TYPES ; i++ )
  {
    m_resiBuffer[i] = new Pel [MAX_CU_SIZE * MAX_CU_SIZE];
  }
#endif

  UInt uiNumPartitions;
#if JVET_C0024_QTBT
  for( i=0 ; i<m_uhTotalDepth ; i++)
#else
  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
#endif
  {
    uiNumPartitions = 1<<( ( m_uhTotalDepth - i - 1 )<<1 );
    UInt uiWidth  = uiMaxWidth  >> i;
    UInt uiHeight = uiMaxHeight >> i;

    m_ppcBestCU[i] = new TComDataCU; m_ppcBestCU[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
    m_ppcTempCU[i] = new TComDataCU; m_ppcTempCU[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );

    m_ppcPredYuvBest[i] = new TComYuv; m_ppcPredYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcResiYuvBest[i] = new TComYuv; m_ppcResiYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcRecoYuvBest[i] = new TComYuv; m_ppcRecoYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);

    m_ppcPredYuvTemp[i] = new TComYuv; m_ppcPredYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcResiYuvTemp[i] = new TComYuv; m_ppcResiYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcRecoYuvTemp[i] = new TComYuv; m_ppcRecoYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);

    m_ppcOrigYuv    [i] = new TComYuv; m_ppcOrigYuv    [i]->create(uiWidth, uiHeight, chromaFormat);

#if VCEG_AZ07_FRUC_MERGE
    m_ppcFRUCBufferCU[i] = new TComDataCU; m_ppcFRUCBufferCU[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
#endif
#if VCEG_AZ07_IMV
    for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
    {
      m_ppcTempCUIMVCache[size][i] = new TComDataCU; m_ppcTempCUIMVCache[size][i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
    }
#endif
#if COM16_C806_OBMC
    m_ppcTempCUWoOBMC [i] = new TComDataCU; m_ppcTempCUWoOBMC[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
#if VCEG_AZ05_BIO
    m_ppcTmpYuv1      [i] = new TComYuv;    m_ppcTmpYuv1      [i]->create( uiWidth+4, uiHeight+4, chromaFormat );
    m_ppcTmpYuv2      [i] = new TComYuv;    m_ppcTmpYuv2      [i]->create( uiWidth+4, uiHeight+4, chromaFormat );
#else
    m_ppcTmpYuv1      [i] = new TComYuv;    m_ppcTmpYuv1      [i]->create( uiWidth, uiHeight, chromaFormat );
    m_ppcTmpYuv2      [i] = new TComYuv;    m_ppcTmpYuv2      [i]->create( uiWidth, uiHeight, chromaFormat );
#endif
    m_ppcPredYuvWoOBMC[i] = new TComYuv;    m_ppcPredYuvWoOBMC[i]->create( uiWidth, uiHeight, chromaFormat );
#endif
  }
#endif  //JVET_C0024_QTBT

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
  for (Int ii=0; ii< NUM_MGR_TYPE; ii++)
  {
    m_pMvFieldSP[ii] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
    m_phInterDirSP[ii] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
    assert( m_pMvFieldSP[ii] != NULL && m_phInterDirSP[ii] != NULL );
  }
#else
  m_pMvFieldSP[0] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
  m_pMvFieldSP[1] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH*2];
  m_phInterDirSP[0] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
  m_phInterDirSP[1] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
  assert( m_pMvFieldSP[0] != NULL && m_phInterDirSP[0] != NULL );
  assert( m_pMvFieldSP[1] != NULL && m_phInterDirSP[1] != NULL );
#endif
#endif

  m_bEncodeDQP                     = false;
  m_stillToCodeChromaQpOffsetFlag  = false;
  m_cuChromaQpOffsetIdxPlus1       = 0;
  m_bFastDeltaQP                   = false;

  // initialize partition order.
  UInt* piTmp = &g_auiZscanToRaster[0];
  initZscanToRaster( m_uhTotalDepth, 1, 0, piTmp);
  initRasterToZscan( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );

  // initialize conversion matrix from partition index to pel
  initRasterToPelXY( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );

#if PVC_AQS
  initFactorValue();
#endif
}

Void TEncCu::destroy()
{
  Int i;

#if JVET_C0024_QTBT
  UInt uiNumWidthIdx = m_uhTotalDepth;
  UInt uiNumHeightIdx = m_uhTotalDepth;

  for (UInt wIdx = 0; wIdx < uiNumWidthIdx; wIdx++)
  {
    for (UInt hIdx=0; hIdx<uiNumHeightIdx; hIdx++)
    {
      if(m_pppcBestCU[wIdx][hIdx])
      {
        m_pppcBestCU[wIdx][hIdx]->destroy();      delete m_pppcBestCU[wIdx][hIdx];      m_pppcBestCU[wIdx][hIdx] = NULL;
      }
      if(m_pppcTempCU[wIdx][hIdx])
      {
        m_pppcTempCU[wIdx][hIdx]->destroy();      delete m_pppcTempCU[wIdx][hIdx];      m_pppcTempCU[wIdx][hIdx] = NULL;
      }
      if(m_pppcPredYuvBest[wIdx][hIdx])
      {
        m_pppcPredYuvBest[wIdx][hIdx]->destroy(); delete m_pppcPredYuvBest[wIdx][hIdx]; m_pppcPredYuvBest[wIdx][hIdx] = NULL;
      }
      if(m_pppcResiYuvBest[wIdx][hIdx])
      {
        m_pppcResiYuvBest[wIdx][hIdx]->destroy(); delete m_pppcResiYuvBest[wIdx][hIdx]; m_pppcResiYuvBest[wIdx][hIdx] = NULL;
      }
      if(m_pppcRecoYuvBest[wIdx][hIdx])
      {
        m_pppcRecoYuvBest[wIdx][hIdx]->destroy(); delete m_pppcRecoYuvBest[wIdx][hIdx]; m_pppcRecoYuvBest[wIdx][hIdx] = NULL;
      }
      if(m_pppcPredYuvTemp[wIdx][hIdx])
      {
        m_pppcPredYuvTemp[wIdx][hIdx]->destroy(); delete m_pppcPredYuvTemp[wIdx][hIdx]; m_pppcPredYuvTemp[wIdx][hIdx] = NULL;
      }
      if(m_pppcResiYuvTemp[wIdx][hIdx])
      {
        m_pppcResiYuvTemp[wIdx][hIdx]->destroy(); delete m_pppcResiYuvTemp[wIdx][hIdx]; m_pppcResiYuvTemp[wIdx][hIdx] = NULL;
      }
      if(m_pppcRecoYuvTemp[wIdx][hIdx])
      {
        m_pppcRecoYuvTemp[wIdx][hIdx]->destroy(); delete m_pppcRecoYuvTemp[wIdx][hIdx]; m_pppcRecoYuvTemp[wIdx][hIdx] = NULL;
      }
      if(m_pppcOrigYuv[wIdx][hIdx])
      {
        m_pppcOrigYuv[wIdx][hIdx]->destroy();     delete m_pppcOrigYuv[wIdx][hIdx];     m_pppcOrigYuv[wIdx][hIdx] = NULL;
      }
#if VCEG_AZ07_FRUC_MERGE
      if(m_pppcFRUCBufferCU[wIdx][hIdx])
      {
        m_pppcFRUCBufferCU[wIdx][hIdx]->destroy();  delete m_pppcFRUCBufferCU[wIdx][hIdx];      m_pppcFRUCBufferCU[wIdx][hIdx] = NULL;
      }
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
      for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
      {
        if(m_pppcTempCUIMVCache[size][wIdx][hIdx])
        {
          m_pppcTempCUIMVCache[size][wIdx][hIdx]->destroy(); delete m_pppcTempCUIMVCache[size][wIdx][hIdx]; m_pppcTempCUIMVCache[size][wIdx][hIdx] = NULL;
        }
      }
#endif
#if COM16_C806_OBMC
      if(m_pppcTempCUWoOBMC[wIdx][hIdx])
      {
        m_pppcTempCUWoOBMC[wIdx][hIdx]->destroy();  delete m_pppcTempCUWoOBMC[wIdx][hIdx];      m_pppcTempCUWoOBMC[wIdx][hIdx] = NULL;
      }

      if(m_pppcTmpYuv1[wIdx][hIdx])
      {
        m_pppcTmpYuv1[wIdx][hIdx]->destroy();       delete m_pppcTmpYuv1[wIdx][hIdx]; m_pppcTmpYuv1[wIdx][hIdx] = NULL;
      }
      if(m_pppcTmpYuv2[wIdx][hIdx])
      {
        m_pppcTmpYuv2[wIdx][hIdx]->destroy();       delete m_pppcTmpYuv2[wIdx][hIdx]; m_pppcTmpYuv2[wIdx][hIdx] = NULL;
      }
      if(m_pppcPredYuvWoOBMC[wIdx][hIdx])
      {
        m_pppcPredYuvWoOBMC[wIdx][hIdx]->destroy(); delete m_pppcPredYuvWoOBMC[wIdx][hIdx]; m_pppcPredYuvWoOBMC[wIdx][hIdx] = NULL;
      }
#endif
    }
    if(m_pppcBestCU[wIdx])
    {
      delete [] m_pppcBestCU[wIdx];
      m_pppcBestCU[wIdx] = NULL;
    }
    if(m_pppcTempCU[wIdx])
    {
      delete [] m_pppcTempCU[wIdx];
      m_pppcTempCU[wIdx] = NULL;
    }

    if(m_pppcPredYuvBest[wIdx])
    {
      delete [] m_pppcPredYuvBest[wIdx];
      m_pppcPredYuvBest[wIdx] = NULL;
    }
    if(m_pppcResiYuvBest[wIdx])
    {
      delete [] m_pppcResiYuvBest[wIdx];
      m_pppcResiYuvBest[wIdx] = NULL;
    }
    if(m_pppcRecoYuvBest[wIdx])
    {
      delete [] m_pppcRecoYuvBest[wIdx];
      m_pppcRecoYuvBest[wIdx] = NULL;
    }
    if(m_pppcPredYuvTemp[wIdx])
    {
      delete [] m_pppcPredYuvTemp[wIdx];
      m_pppcPredYuvTemp[wIdx] = NULL;
    }
    if(m_pppcResiYuvTemp[wIdx])
    {
      delete [] m_pppcResiYuvTemp[wIdx];
      m_pppcResiYuvTemp[wIdx] = NULL;
    }
    if(m_pppcRecoYuvTemp[wIdx])
    {
      delete [] m_pppcRecoYuvTemp[wIdx];
      m_pppcRecoYuvTemp[wIdx] = NULL;
    }
    if(m_pppcOrigYuv[wIdx])
    {
      delete [] m_pppcOrigYuv[wIdx];
      m_pppcOrigYuv[wIdx] = NULL;
    }
#if VCEG_AZ07_FRUC_MERGE
    if(m_pppcFRUCBufferCU[wIdx])
    {
      delete [] m_pppcFRUCBufferCU[wIdx];
      m_pppcFRUCBufferCU[wIdx] = NULL;
    }
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
    for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
    {
      delete [] m_pppcTempCUIMVCache[size][wIdx];
      m_pppcTempCUIMVCache[size][wIdx] = NULL;
    }
#endif
#if COM16_C806_OBMC
    if(m_pppcTempCUWoOBMC[wIdx])
    {
      delete [] m_pppcTempCUWoOBMC[wIdx];
      m_pppcTempCUWoOBMC[wIdx] = NULL;
    }

    if(m_pppcTmpYuv1[wIdx])
    {
      delete [] m_pppcTmpYuv1[wIdx];
      m_pppcTmpYuv1[wIdx] = NULL;
    }
    if(m_pppcTmpYuv2[wIdx])
    {
      delete [] m_pppcTmpYuv2[wIdx];
      m_pppcTmpYuv2[wIdx] = NULL;
    }
    if(m_pppcPredYuvWoOBMC[wIdx])
    {
      delete [] m_pppcPredYuvWoOBMC[wIdx];
      m_pppcPredYuvWoOBMC[wIdx] = NULL;
    }
#endif
  }

  if(m_pppcBestCU)
  {
    delete [] m_pppcBestCU;
    m_pppcBestCU = NULL;
  }
  if(m_pppcTempCU)
  {
    delete [] m_pppcTempCU;
    m_pppcTempCU = NULL;
  }

#if JVET_C0024_FAST_MRG
  for (UInt idx=0; idx<MRG_MAX_NUM_CANDS; idx++)
  {
    if(m_pcMrgPredTempYuv[idx])
    {
      m_pcMrgPredTempYuv[idx]->destroy();
      delete m_pcMrgPredTempYuv[idx];
      m_pcMrgPredTempYuv[idx] = NULL;
    }
  }

#endif
  if(m_pppcPredYuvBest)
  {
    delete [] m_pppcPredYuvBest;
    m_pppcPredYuvBest = NULL;
  }
  if(m_pppcResiYuvBest)
  {
    delete [] m_pppcResiYuvBest;
    m_pppcResiYuvBest = NULL;
  }
  if(m_pppcRecoYuvBest)
  {
    delete [] m_pppcRecoYuvBest;
    m_pppcRecoYuvBest = NULL;
  }
  if(m_pppcPredYuvTemp)
  {
    delete [] m_pppcPredYuvTemp;
    m_pppcPredYuvTemp = NULL;
  }
  if(m_pppcResiYuvTemp)
  {
    delete [] m_pppcResiYuvTemp;
    m_pppcResiYuvTemp = NULL;
  }
  if(m_pppcRecoYuvTemp)
  {
    delete [] m_pppcRecoYuvTemp;
    m_pppcRecoYuvTemp = NULL;
  }
  if(m_pppcOrigYuv)
  {
    delete [] m_pppcOrigYuv;
    m_pppcOrigYuv = NULL;
  }
#if VCEG_AZ07_FRUC_MERGE
  if(m_pppcFRUCBufferCU)
  {
    delete [] m_pppcFRUCBufferCU;
    m_pppcFRUCBufferCU = NULL;
  }
#endif
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    delete [] m_pppcTempCUIMVCache[size];
    m_pppcTempCUIMVCache[size] = NULL;
  }
#endif
#if COM16_C806_OBMC
  if(m_pppcTempCUWoOBMC)
  {
    delete [] m_pppcTempCUWoOBMC;
    m_pppcTempCUWoOBMC = NULL;
  }

  if(m_pppcTmpYuv1)
  {
    delete [] m_pppcTmpYuv1;
    m_pppcTmpYuv1 = NULL;
  }
  if(m_pppcTmpYuv2)
  {
    delete [] m_pppcTmpYuv2;
    m_pppcTmpYuv2 = NULL;
  }
  if(m_pppcPredYuvWoOBMC)
  {
    delete [] m_pppcPredYuvWoOBMC;
    m_pppcPredYuvWoOBMC = NULL;
  }
#endif

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
  for (UInt ui=0;ui<NUM_MGR_TYPE;ui++) 
#else
  for (UInt ui=0;ui<2;ui++)
#endif
  {
    if( m_pMvFieldSP[ui] != NULL )
    {
      delete [] m_pMvFieldSP[ui];
      m_pMvFieldSP[ui] = NULL;
    }
    if( m_phInterDirSP[ui] != NULL )
    {
      delete [] m_phInterDirSP[ui];
      m_phInterDirSP[ui] = NULL;
    }
  }
#endif

#else
  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
  {
    if(m_ppcBestCU[i])
    {
      m_ppcBestCU[i]->destroy();      delete m_ppcBestCU[i];      m_ppcBestCU[i] = NULL;
    }
    if(m_ppcTempCU[i])
    {
      m_ppcTempCU[i]->destroy();      delete m_ppcTempCU[i];      m_ppcTempCU[i] = NULL;
    }
    if(m_ppcPredYuvBest[i])
    {
      m_ppcPredYuvBest[i]->destroy(); delete m_ppcPredYuvBest[i]; m_ppcPredYuvBest[i] = NULL;
    }
    if(m_ppcResiYuvBest[i])
    {
      m_ppcResiYuvBest[i]->destroy(); delete m_ppcResiYuvBest[i]; m_ppcResiYuvBest[i] = NULL;
    }
    if(m_ppcRecoYuvBest[i])
    {
      m_ppcRecoYuvBest[i]->destroy(); delete m_ppcRecoYuvBest[i]; m_ppcRecoYuvBest[i] = NULL;
    }
    if(m_ppcPredYuvTemp[i])
    {
      m_ppcPredYuvTemp[i]->destroy(); delete m_ppcPredYuvTemp[i]; m_ppcPredYuvTemp[i] = NULL;
    }
    if(m_ppcResiYuvTemp[i])
    {
      m_ppcResiYuvTemp[i]->destroy(); delete m_ppcResiYuvTemp[i]; m_ppcResiYuvTemp[i] = NULL;
    }
    if(m_ppcRecoYuvTemp[i])
    {
      m_ppcRecoYuvTemp[i]->destroy(); delete m_ppcRecoYuvTemp[i]; m_ppcRecoYuvTemp[i] = NULL;
    }
    if(m_ppcOrigYuv[i])
    {
      m_ppcOrigYuv[i]->destroy();     delete m_ppcOrigYuv[i];     m_ppcOrigYuv[i] = NULL;
    }
#if VCEG_AZ07_FRUC_MERGE
    if(m_ppcFRUCBufferCU[i])
    {
      m_ppcFRUCBufferCU[i]->destroy();  delete m_ppcFRUCBufferCU[i];      m_ppcFRUCBufferCU[i] = NULL;
    }
#endif
#if VCEG_AZ07_IMV
    for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
    {
      if(m_ppcTempCUIMVCache[size][i])
      {
        m_ppcTempCUIMVCache[size][i]->destroy(); delete m_ppcTempCUIMVCache[size][i]; m_ppcTempCUIMVCache[size][i] = NULL;
      }
    }
#endif
#if COM16_C806_OBMC
    if(m_ppcTempCUWoOBMC[i])
    {
      m_ppcTempCUWoOBMC[i]->destroy();  delete m_ppcTempCUWoOBMC[i];      m_ppcTempCUWoOBMC[i] = NULL;
    }

    if(m_ppcTmpYuv1[i])
    {
      m_ppcTmpYuv1[i]->destroy();       delete m_ppcTmpYuv1[i]; m_ppcTmpYuv1[i] = NULL;
    }
    if(m_ppcTmpYuv2[i])
    {
      m_ppcTmpYuv2[i]->destroy();       delete m_ppcTmpYuv2[i]; m_ppcTmpYuv2[i] = NULL;
    }
    if(m_ppcPredYuvWoOBMC[i])
    {
      m_ppcPredYuvWoOBMC[i]->destroy(); delete m_ppcPredYuvWoOBMC[i]; m_ppcPredYuvWoOBMC[i] = NULL;
    }
#endif
  }
  if(m_ppcBestCU)
  {
    delete [] m_ppcBestCU;
    m_ppcBestCU = NULL;
  }
  if(m_ppcTempCU)
  {
    delete [] m_ppcTempCU;
    m_ppcTempCU = NULL;
  }

  if(m_ppcPredYuvBest)
  {
    delete [] m_ppcPredYuvBest;
    m_ppcPredYuvBest = NULL;
  }
  if(m_ppcResiYuvBest)
  {
    delete [] m_ppcResiYuvBest;
    m_ppcResiYuvBest = NULL;
  }
  if(m_ppcRecoYuvBest)
  {
    delete [] m_ppcRecoYuvBest;
    m_ppcRecoYuvBest = NULL;
  }
  if(m_ppcPredYuvTemp)
  {
    delete [] m_ppcPredYuvTemp;
    m_ppcPredYuvTemp = NULL;
  }
  if(m_ppcResiYuvTemp)
  {
    delete [] m_ppcResiYuvTemp;
    m_ppcResiYuvTemp = NULL;
  }
  if(m_ppcRecoYuvTemp)
  {
    delete [] m_ppcRecoYuvTemp;
    m_ppcRecoYuvTemp = NULL;
  }
  if(m_ppcOrigYuv)
  {
    delete [] m_ppcOrigYuv;
    m_ppcOrigYuv = NULL;
  }
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  for (UInt ui=0;ui<2;ui++)
  {
    if( m_pMvFieldSP[ui] != NULL )
    {
      delete [] m_pMvFieldSP[ui];
      m_pMvFieldSP[ui] = NULL;
    }
    if( m_phInterDirSP[ui] != NULL )
    {
      delete [] m_phInterDirSP[ui];
      m_phInterDirSP[ui] = NULL;
    }
  }
#endif

#if VCEG_AZ07_FRUC_MERGE
  if(m_ppcFRUCBufferCU)
  {
    delete [] m_ppcFRUCBufferCU;
    m_ppcFRUCBufferCU = NULL;
  }
#endif
#if VCEG_AZ07_IMV
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    delete [] m_ppcTempCUIMVCache[size];
    m_ppcTempCUIMVCache[size] = NULL;
  }
#endif
#if COM16_C806_OBMC
  if(m_ppcTempCUWoOBMC)
  {
    delete [] m_ppcTempCUWoOBMC;
    m_ppcTempCUWoOBMC = NULL;
  }

  if(m_ppcTmpYuv1)
  {
    delete [] m_ppcTmpYuv1;
    m_ppcTmpYuv1 = NULL;
  }
  if(m_ppcTmpYuv2)
  {
    delete [] m_ppcTmpYuv2;
    m_ppcTmpYuv2 = NULL;
  }
  if(m_ppcPredYuvWoOBMC)
  {
    delete [] m_ppcPredYuvWoOBMC;
    m_ppcPredYuvWoOBMC = NULL;
  }
#endif
#endif //JVET_C0024_QTBT

#if COM16_C806_LARGE_CTU
  for( i = 0 ; i < NUMBER_OF_STORED_RESIDUAL_TYPES ; i++ )
  {
    if( m_resiBuffer[i] )
    {
      delete [] m_resiBuffer[i];
      m_resiBuffer[i] = NULL;
    }
  }
#endif
#if BARC_FUNC
#if BARC_ONLY_FOR_LUMA
  if (m_reconLumaBufferHM)
  {
	  delete[] m_reconLumaBufferHM;
	  m_reconLumaBufferHM = NULL;
  }
#if COMPACT_RESOLUTION_CNN
  if (m_reconBufferBicubic)
  {
	  delete[] m_reconBufferBicubic;
	  m_reconBufferBicubic = NULL;
  }
  if (m_DScaleBufferBicubic)
  {
	  delete[] m_DScaleBufferBicubic;
	  m_DScaleBufferBicubic = NULL;
  }
#endif
#else
  for (i = 0; i < MAX_NUM_COMPONENT; i++)
  {
    if (m_reconBuffer[i])
    {
      delete[] m_reconBuffer[i];
      m_reconBuffer[i] = NULL;
    }
#if COMPACT_RESOLUTION_CNN
	if (m_reconBufferBicubic[i])
	{
		delete[] m_reconBufferBicubic[i];
		m_reconBufferBicubic[i] = NULL;
	}
#endif
  }
  delete[] m_reconLumaBuffer;
  m_reconLumaBuffer = NULL;
#endif
#endif
}

/** \param    pcEncTop      pointer of encoder class
 */
Void TEncCu::init( TEncTop* pcEncTop )
{
  m_pcEncCfg           = pcEncTop;
  m_pcPredSearch       = pcEncTop->getPredSearch();
  m_pcTrQuant          = pcEncTop->getTrQuant();
  m_pcRdCost           = pcEncTop->getRdCost();

  m_pcEntropyCoder     = pcEncTop->getEntropyCoder();
  m_pcBinCABAC         = pcEncTop->getBinCABAC();

#if JVET_C0024_QTBT
  m_ppppcRDSbacCoder   = pcEncTop->getRDSbacCoder();
#else
  m_pppcRDSbacCoder    = pcEncTop->getRDSbacCoder();
#endif
  m_pcRDGoOnSbacCoder  = pcEncTop->getRDGoOnSbacCoder();

  m_pcRateCtrl         = pcEncTop->getRateCtrl();
#if WCG_LUMA_DQP_CM_SCALE
  m_LumaQPOffset=0;
  if (m_pcEncCfg->getUseLumaDeltaQp())
  {
    initLumaDeltaQpLUT();
  }
#endif
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/** 
 \param  pCtu pointer of CU data class
 */
Void TEncCu::compressCtu( TComDataCU* pCtu )
{
#if JVET_C0024_QTBT
  pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
  if (!pCtu->getSlice()->isIntra())
  {
    pCtu->getPic()->clearAllIntMv();
    pCtu->getPic()->clearAllSkiped();
    pCtu->getPic()->clearAllInter();
    pCtu->getPic()->clearAllIntra();
  }
#endif

  // initialize CU data
#if JVET_C0024_QTBT
  pCtu->getPic()->setCodedAreaInCTU(0);
  UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiWidthIdx = g_aucConvertToBit[uiCTUSize];
  UInt uiHeightIdx = g_aucConvertToBit[uiCTUSize];

  m_pppcBestCU[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
  m_pppcTempCU[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#if JVET_C0024_DELTA_QP_FIX
  if (pCtu->getSlice()->getPPS()->getUseDQP())
  {
    setdQPFlag(true);   // JVET_C0024_DELTA_QP_FIX2 additional bug fix

    Char qp = pCtu->getCtuLastCodedQP();
    m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setCodedQP( qp );
    m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setCodedQP( qp );
    m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setQuLastCodedQP( qp );
    m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setQuLastCodedQP( qp );
    m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setQuPartIdx(0); // init current quantization unit partition index
    m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setQuPartIdx(0); // init current quantization unit partition index
  }
#endif

#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    m_pppcTempCUIMVCache[size][uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic() , pCtu->getCtuRsAddr() );
  }
#endif
#if COM16_C806_OBMC
  m_pppcTempCUWoOBMC[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_pppcFRUCBufferCU[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#endif
#if  SIF_SPLIT_ONLY
  m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setSIFToCheck(0);
  m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setSIFToCheck(0);
#endif
  // analysis of CU
  DEBUG_STRING_NEW(sDebug)

  xCompressCU( m_pppcBestCU[uiWidthIdx][uiHeightIdx], m_pppcTempCU[uiWidthIdx][uiHeightIdx], 0, uiCTUSize, uiCTUSize, 0 DEBUG_STRING_PASS_INTO(sDebug) );     
#if BARC
#if COMPACT_RESOLUTION_CNN
  if (pCtu->getSlice()->getSliceType() == I_SLICE)
  {
#if  BARC_DIMD
	  const TComSPS &sps = *(pCtu->getSlice()->getSPS());
	  UInt remainX = sps.getPicWidthInLumaSamples() - pCtu->getCUPelX();
	  UInt remainY = sps.getPicHeightInLumaSamples() - pCtu->getCUPelY();
	  if (!(remainX < uiCTUSize && remainX % 16 != 0) && !(remainY < uiCTUSize && remainY % 16 != 0))
	  {
#endif

		  copyRecToTmp(pCtu, true);//copy picYuvRec  to m_reconBuffer[]
		  TComDataCU* pCtuLuma = pCtu->getPic()->getCtuTmp();
		  pCtuLuma->copyCtu(pCtu);

		  Double hmCost = pCtu->getTotalCost();
		  UInt   hmDistortion = pCtu->getTotalDistortion();

		  pCtu->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  pCtu->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  pCtu->getPic()->setCodedAreaInCTU(0);
		  pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize >> MIN_CU_LOG2, uiCTUSize >> MIN_CU_LOG2); //only used for affine merge code or not

		  UInt uiBARCCTUSize = (uiCTUSize >> BARC_DScaleRatio);
		  UInt uiBARCWidthIdx = g_aucConvertToBit[uiBARCCTUSize];
		  UInt uiBARCHeightIdx = g_aucConvertToBit[uiBARCCTUSize];

		  m_ppppcRDSbacCoder[uiBARCWidthIdx][uiBARCHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);

#if BARC_CORE
		  TComSlice*  pcSlice = pCtu->getSlice();
		  initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, true);
#endif
		  m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  xCompressCU(m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx], m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx], BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize, 0 DEBUG_STRING_PASS_INTO(sDebug));
		  upSampleRec(pCtu, CHANNEL_TYPE_LUMA);
#if BARC_CORE
		  initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, false);
#endif
		  // save rec yuv and coding information for filter-up and filter-down mode.
		  copyBicubicRecToTmp(pCtu, true);//copy picYuvRec  to m_reconBuffer[]
		  reComputeBARCDistortion(pCtu, CHANNEL_TYPE_LUMA);
		  TComDataCU* pCtuLumaBicubic = pCtu->getPic()->getCtuTmpBicubic();
		  pCtuLumaBicubic->copyCtu(pCtu);


		  Double bicubicDSCost = pCtu->getTotalCost();//down-sample is implemented through Bicubic
		  TComPicYuv* pcPicYuvDScaleOrgBackup = pCtu->getPic()->getPicYuvDScaleOrg();
		  pCtu->getPic()->setPicYuvDScaleOrg(pCtu->getPic()->getPicYuvCRCNNOrg());//set PicYuvDScaleOrg to PicYuvCRCNNOrg

		  pCtu->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  pCtu->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  pCtu->getPic()->setCodedAreaInCTU(0);
		  m_ppppcRDSbacCoder[uiBARCWidthIdx][uiBARCHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);
#if BARC_CORE
		  initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, true);
#endif
		  m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
		  xCompressCU(m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx], m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx], BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize, 0 DEBUG_STRING_PASS_INTO(sDebug));
		  //upSampleRec(pCtu, CHANNEL_TYPE_LUMA);
#if BARC_CORE
		  initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, false);
#endif
		  upSampleRecDeep(pCtu, CHANNEL_TYPE_LUMA);
		  reComputeBARCDistortion(pCtu, CHANNEL_TYPE_LUMA);

		  Double crcnnCost = pCtu->getTotalCost();
		  pCtu->getPic()->setPicYuvDScaleOrg(pcPicYuvDScaleOrgBackup);
		  //#if LIYE_DEBUG
		  //		  cout << "\n" << "CtuAddr: "<<pCtu->getCtuRsAddr<<" "<<hmCost << " " << bicubicDSCost << " " << crcnnCost << "\n";
		  //#endif
		  if (bicubicDSCost > hmCost && crcnnCost > hmCost)
		  {
			  pCtu->copyCtu(pCtuLuma);
			  copyRecToTmp(pCtu, false);
			  downSampleRec(pCtu, CHANNEL_TYPE_LUMA);
		  }
		  else if (bicubicDSCost < crcnnCost)
		  {
			  pCtu->copyCtu(pCtuLumaBicubic);
			  copyBicubicRecToTmp(pCtu, false);
			  reComputeBARCDistortion(pCtu, CHANNEL_TYPE_LUMA);
		  }
#if BARC_DIMD
	  }
	  else
	  {
		  downSampleRec(pCtu, CHANNEL_TYPE_LUMA);
	  }
#endif
  }
#else
  if (pCtu->getSlice()->getSliceType() == I_SLICE)
  {
#if  BARC_DIMD
    const TComSPS &sps = *(pCtu->getSlice()->getSPS());
    UInt remainX = sps.getPicWidthInLumaSamples() - pCtu->getCUPelX();
    UInt remainY = sps.getPicHeightInLumaSamples() - pCtu->getCUPelY();
    if (!(remainX < uiCTUSize && remainX % 16 != 0) && !(remainY < uiCTUSize && remainY % 16 != 0))
    {
#endif
#if BARC 
      copyRecToTmp(pCtu, true);
      TComDataCU* pCtuLuma = pCtu->getPic()->getCtuTmp();
      pCtuLuma->copyCtu(pCtu);
#endif
      Double hmCost = pCtu->getTotalCost();
      UInt   hmDistortion = pCtu->getTotalDistortion();

      pCtu->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
      pCtu->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
      pCtu->getPic()->setCodedAreaInCTU(0);
      pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize >> MIN_CU_LOG2, uiCTUSize >> MIN_CU_LOG2); //only used for affine merge code or not

      UInt uiBARCCTUSize = (uiCTUSize >> BARC_DScaleRatio);
      UInt uiBARCWidthIdx = g_aucConvertToBit[uiBARCCTUSize];
      UInt uiBARCHeightIdx = g_aucConvertToBit[uiBARCCTUSize];

      m_ppppcRDSbacCoder[uiBARCWidthIdx][uiBARCHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);

#if BARC_CORE
      TComSlice*  pcSlice = pCtu->getSlice();
      initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, true);
#endif
      m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
      m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
      m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
      m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx]->initDscaleCtu(CHANNEL_TYPE_LUMA, pCtu->getPic(), pCtu->getCtuRsAddr());
      xCompressCU(m_pppcBestCU[uiBARCWidthIdx][uiBARCHeightIdx], m_pppcTempCU[uiBARCWidthIdx][uiBARCHeightIdx], BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize, 0 DEBUG_STRING_PASS_INTO(sDebug));
      upSampleRec(pCtu, CHANNEL_TYPE_LUMA);
#if BARC_CORE
      initEncSliceBARC(CHANNEL_TYPE_LUMA, pCtu->getPic(), 0, pcSlice, false);
#endif
      reComputeBARCDistortion(pCtu, CHANNEL_TYPE_LUMA);
      if (pCtu->getTotalCost() > hmCost)
      {
        pCtu->copyCtu(pCtuLuma);
        copyRecToTmp(pCtu, false);
        downSampleRec(pCtu, CHANNEL_TYPE_LUMA);
      }
#if BARC_DIMD
    }
    else
    {
      downSampleRec(pCtu, CHANNEL_TYPE_LUMA);
    }
#endif
  }
#endif
#endif
#else
  m_ppcBestCU[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
  m_ppcTempCU[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#if VCEG_AZ07_IMV
  for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
  {
    m_ppcTempCUIMVCache[size][0]->initCtu( pCtu->getPic() , pCtu->getCtuRsAddr() );
  }
#endif
#if COM16_C806_OBMC
  m_ppcTempCUWoOBMC[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_ppcFRUCBufferCU[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#endif
#if  SIF_SPLIT_ONLY
  m_ppcBestCU[0]->setSIFToCheck(0);
  m_ppcTempCU[0]->setSIFToCheck(0);
#endif
  // analysis of CU
  DEBUG_STRING_NEW(sDebug)

  xCompressCU( m_ppcBestCU[0], m_ppcTempCU[0], 0 DEBUG_STRING_PASS_INTO(sDebug) );
#endif
  DEBUG_STRING_OUTPUT(std::cout, sDebug)

#if JVET_C0024_QTBT
    if (pCtu->getSlice()->isIntra())
    {
      pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
      // initialize CU data
      pCtu->getPic()->setCodedAreaInCTU(0);
      m_pppcBestCU[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
      m_pppcTempCU[uiWidthIdx][uiHeightIdx]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
#if JVET_C0024_DELTA_QP_FIX
      if (pCtu->getSlice()->getPPS()->getUseDQP())
      {
        Char qp = pCtu->getCtuLastCodedQP();
        m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setCodedQP( qp );
        m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setCodedQP( qp );
        m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setQuLastCodedQP( qp );
        m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setQuLastCodedQP( qp );
        m_pppcBestCU[uiWidthIdx][uiHeightIdx]->setQuPartIdx(0); // init current quantization unit partition index
        m_pppcTempCU[uiWidthIdx][uiHeightIdx]->setQuPartIdx(0); // init current quantization unit partition index
      }
#endif

      Double lumaCTBTotalCost = pCtu->getTotalCost();
      Distortion lumaCTBTotalDistortion = pCtu->getTotalDistortion();
      UInt lumaCTBTotalBits = pCtu->getTotalBits();
      UInt lumaCTBTotalBins = pCtu->getTotalBins();
      xCompressCU(m_pppcBestCU[uiWidthIdx][uiHeightIdx], m_pppcTempCU[uiWidthIdx][uiHeightIdx], 0, uiCTUSize, uiCTUSize, 0 DEBUG_STRING_PASS_INTO(sDebug));
      pCtu->getTotalCost() += lumaCTBTotalCost;
      pCtu->getTotalDistortion() += lumaCTBTotalDistortion;
      pCtu->getTotalBits() += lumaCTBTotalBits;
      pCtu->getTotalBins() += lumaCTBTotalBins;
    }
#endif

#if ADAPTIVE_QP_SELECTION
  if( m_pcEncCfg->getUseAdaptQpSelect() )
  {
    if(pCtu->getSlice()->getSliceType()!=I_SLICE) //IIII
    {
      xCtuCollectARLStats( pCtu );
    }
  }
#endif
}
/** \param  pCtu  pointer of CU data class
 */
Void TEncCu::encodeCtu ( TComDataCU* pCtu )
{
#if JVET_C0024_QTBT
  pCtu->getSlice()->setTextType(CHANNEL_TYPE_LUMA);
#endif
  if ( pCtu->getSlice()->getPPS()->getUseDQP() )
  {
    setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
    Char qp = pCtu->getCtuLastCodedQP();
    pCtu->setQuPartIdx(0); 
    pCtu->setQuLastCodedQP( qp );
    pCtu->setCodedQP( qp );
#endif
  }

  if ( pCtu->getSlice()->getUseChromaQpAdj() )
  {
    setCodeChromaQpAdjFlag(true);
  }

  // Encode CU data
#if JVET_C0024_QTBT
  UInt uiCTUSize = pCtu->getSlice()->getSPS()->getCTUSize();
  pCtu->getPic()->setCodedAreaInCTU(0);
  pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2); //only used for affine merge code or not

#if BARC_CORE
  if (pCtu->getSlice()->isIntra())
  {
	  m_pcEntropyCoder->encodeBARCFlag(pCtu, CHANNEL_TYPE_LUMA);
#if BARC_UP_SAMPLE_RDO
	  if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
	  {
		  m_pcEntropyCoder->encodeUpCNNFlag(pCtu, 0);
	  }
#endif
  }
  UInt uiBARCCTUSize = uiCTUSize;
  if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
  {
      uiBARCCTUSize >>= BARC_DScaleRatio;
      xEncodeCU( pCtu, 0, BARC_DScaleRatio, uiBARCCTUSize, uiBARCCTUSize );
  }
  else
#endif
  xEncodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize );  
#else
  xEncodeCU( pCtu, 0, 0 );
#endif

#if JVET_C0024_QTBT
  if (pCtu->getSlice()->isIntra())
  {
    pCtu->getSlice()->setTextType(CHANNEL_TYPE_CHROMA);
    // Encode CU data
    pCtu->getPic()->setCodedAreaInCTU(0);
#if JVET_E0077_ENHANCED_LM
    pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize >> MIN_CU_LOG2, uiCTUSize >> MIN_CU_LOG2);
#endif
#if JVET_C0024_DELTA_QP_FIX
    if (pCtu->getSlice()->getPPS()->getUseDQP())
    {
      setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
      Char qp = pCtu->getCtuLastCodedQP();
      pCtu->setQuPartIdx(0); 
      pCtu->setQuLastCodedQP( qp );
      pCtu->setCodedQP( qp );
#endif
    }
#endif
#if JVET_E0062_MULTI_DMS
    pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2); 
#endif

    xEncodeCU( pCtu, 0, 0, uiCTUSize, uiCTUSize );  
}
#if JVET_E0062_MULTI_DMS
  pCtu->getPic()->setCodedBlkInCTU(false, 0, 0, uiCTUSize>>MIN_CU_LOG2, uiCTUSize>>MIN_CU_LOG2); 
#endif
  // --- write terminating bit ---
  finishCU(pCtu,0);
#endif
}

#if !JVET_C0024_QTBT
// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
//! Derive small set of test modes for AMP encoder speed-up
#if AMP_ENC_SPEEDUP
#if AMP_MRG
Void TEncCu::deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver, Bool &bTestMergeAMP_Hor, Bool &bTestMergeAMP_Ver
  )
#else
Void TEncCu::deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver
  )
#endif
{
  if ( pcBestCU->getPartitionSize(0) == SIZE_2NxN 
      )
  {
    bTestAMP_Hor = true;
  }
  else if ( pcBestCU->getPartitionSize(0) == SIZE_Nx2N 
          )
  {
    bTestAMP_Ver = true;
  }
  else if ( pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->getMergeFlag(0) == false && pcBestCU->isSkipped(0) == false 
          )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  }

#if AMP_MRG
  //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( eParentPartSize == NUMBER_OF_PART_SIZES ) //! if parent is intra
  {
    if ( pcBestCU->getPartitionSize(0) == SIZE_2NxN 
        )
    {
      bTestMergeAMP_Hor = true;
    }
    else if ( pcBestCU->getPartitionSize(0) == SIZE_Nx2N 
            )
    {
      bTestMergeAMP_Ver = true;
    }
  }

  if ( pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->isSkipped(0) == false )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( pcBestCU->getWidth(0) == 64 )
  {
    bTestAMP_Hor = false;
    bTestAMP_Ver = false;
  }
#else
  //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  }

  if ( eParentPartSize == SIZE_2Nx2N )
  {
    bTestAMP_Hor = false;
    bTestAMP_Ver = false;
  }
#endif
}
#endif
#endif


// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
#if WCG_LUMA_DQP_CM_SCALE
Void TEncCu::initLumaDeltaQpLUT() {
  const LumaLevelToDeltaQPMapping &mapping = m_pcEncCfg->getLumaLevelToDeltaQPMapping();

  if (!mapping.isEnabled())
  {
    return;
  }

  // map the sparse LumaLevelToDeltaQPMapping.mapping to a fully populated linear table.
  Int         lastDeltaQPValue = 0;
  std::size_t nextSparseIndex = 0;
  for (Int index = 0; index<LUMA_LEVEL_TO_DQP_LUT_MAXSIZE; index++)
  {
    while (nextSparseIndex < mapping.mapping.size() && index >= mapping.mapping[nextSparseIndex].first)
    {
      lastDeltaQPValue = mapping.mapping[nextSparseIndex].second;
      nextSparseIndex++;
    }
    m_lumaLevelToDeltaQPLUT[index] = lastDeltaQPValue;
  }
}

// Get QP offset derived from Luma activity
Int TEncCu::calculateLumaDQP(TComDataCU *pcCU, const UInt uiAbsPartIdx, const TComYuv * pOrgYuv)
{
  const Int      stride  = pOrgYuv->getStride(COMPONENT_Y);  
  Int      width   = pcCU->getWidth(uiAbsPartIdx);
  Int      height  = pcCU->getHeight(uiAbsPartIdx);

  // limit the block by picture size
  const TComSPS* pSPS = pcCU->getSlice()->getSPS();
  if ( pcCU->getCUPelX() + width > pSPS->getPicWidthInLumaSamples())
    width = pSPS->getPicWidthInLumaSamples() - pcCU->getCUPelX();

  if ( pcCU->getCUPelY() + height > pSPS->getPicHeightInLumaSamples())
    height = pSPS->getPicHeightInLumaSamples() - pcCU->getCUPelY();

  // Get Luma
  Int Sum = 0;
  Double avg=0;

  const Pel *pY = pOrgYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);

  for (Int y = 0; y < height; y++)
  {
    for (Int x = 0; x < width; x++)
    {
      Sum += pY[x];
    }
    pY += stride;
  }
  avg = (double)Sum/(width*height); 

  Int lumaIdx = min(Int(avg+0.5), LUMA_LEVEL_TO_DQP_LUT_MAXSIZE-1);
  Int QP = m_lumaLevelToDeltaQPLUT[lumaIdx];

  return QP;
}
#endif

/** Compress a CU block recursively with enabling sub-CTU-level delta QP
 *  - for loop of QP value to compress the current CU with all possible QP
*/
#if JVET_C0024_QTBT
Void TEncCu::CompressBT( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug_), Int* piParam, Double* pdParam, UInt uiSplitConstrain )
{
  TComPic* pcPic = rpcBestCU->getPic();
  DEBUG_STRING_NEW(sDebug)
    const TComPPS &pps=*(rpcTempCU->getSlice()->getPPS());
  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  TComSlice * pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());
#if BT_FORCE
  const UInt uiLPelX = rpcBestCU->getCUPelX();
  const UInt uiRPelX = uiLPelX + rpcBestCU->getWidth(0) - 1;
  const UInt uiTPelY = rpcBestCU->getCUPelY();
  const UInt uiBPelY = uiTPelY + rpcBestCU->getHeight(0) - 1;
#if BARC_BTFORCE
  UInt uiPelXInCTUTmp = uiLPelX - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTUTmp = uiTPelY - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTUTmp + uiWidth)*rpcBestCU->getBARCFlag(rpcBestCU->getTextType());
  UInt deltaY = (uiPelYInCTUTmp + uiHeight)*rpcBestCU->getBARCFlag(rpcBestCU->getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#else
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#endif
#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(rpcBestCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
#endif
#if JVET_C0024_DELTA_QP_FIX
  UInt uiQTWidth  = sps.getCTUSize()>>uiDepth;
  UInt uiQTHeight = sps.getCTUSize()>>uiDepth;
#if TT_TREE_DS1
  UInt uiBTDepth  = rpcTempCU->getBTDepthFromBTSplitMode(uiBTSplitMode);
  const UInt uiQTBTDepth = (uiDepth<<1) + (g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight]);
#else
  UInt uiBTDepth  = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight];
  const UInt uiQTBTDepth = (uiDepth<<1) + uiBTDepth;
#endif
  const UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif
#if BARC
  Bool isBARCFlag = rpcBestCU->getBARCFlag(rpcBestCU->getTextType());
#endif
  Int iC1x = uiWidth >> (MIN_CU_LOG2+2);
  Int iC1y = 0;
  Int iC2x = (uiWidth >> (MIN_CU_LOG2+1)) + (uiWidth >> (MIN_CU_LOG2+2)) - 1;
  Int iC2y = (uiHeight >> MIN_CU_LOG2) - 1;
  iC2x = iC2x<0? 0 : iC2x;
  iC2y = iC2y<0? 0 : iC2y;
  UInt uiZorderTL = 0; //top left
  UInt uiZorderC1 = g_auiRasterToZscan[iC1y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC1x]; //center, closer to TL
  UInt uiZorderC2 = g_auiRasterToZscan[iC2y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC2x]; //center, closer to BR
  UInt uiZorderBR = g_auiRasterToZscan[((uiHeight>> MIN_CU_LOG2)-1) * (sps.getCTUSize()>> MIN_CU_LOG2) + (uiWidth>> MIN_CU_LOG2)-1];   //bottom-right part.
  UInt uiPelXInCTU = rpcBestCU->getCUPelX() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = rpcBestCU->getCUPelY() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelY();

  //parameter reading
  Int  iSplitDir       = piParam[0];
  Int  iQP             = piParam[1];
  Int  lastCodedQP     = piParam[2];
  Bool bBTHorRmvEnable = piParam[3];
  Bool bBTVerRmvEnable = piParam[4];
  piParam[6] = 0;

#if TT_BUGFIX_ZIDX
  if(iSplitDir==0)
  {
    iC1x = 0;
    iC1y = uiHeight >> (MIN_CU_LOG2+2);
    iC2x = (uiWidth >> MIN_CU_LOG2) - 1;
    iC2y = (uiHeight >> (MIN_CU_LOG2+1)) + (uiHeight >> (MIN_CU_LOG2+2)) - 1;
    iC2x = iC2x<0? 0 : iC2x;
    iC2y = iC2y<0? 0 : iC2y;
    uiZorderC1 = g_auiRasterToZscan[iC1y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC1x]; //center, closer to TL
    uiZorderC2 = g_auiRasterToZscan[iC2y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC2x]; //center, closer to BR
  }
#endif

  //initial parameter related to split dir
  UInt uiWidthSubPart   = iSplitDir==1? (uiWidth>>1)  : uiWidth;
  UInt uiHeightSubPart  = iSplitDir==0? (uiHeight>>1) : uiHeight;
#if TT_TREE_DS1
  UInt uiSTSplitMode    = iSplitDir==0? ST_HORI_BT    : ST_VERT_BT;  //bt split mode of the sub-cu
  UInt uiSTSplitMode2   = iSplitDir==0? ST_VERT_BT    : ST_HORI_BT;  //the other bt split mode of the sub-cu, used in redundancy removal
#else
  UInt uiSTSplitMode    = iSplitDir==0? 1 : 2;
  UInt uiSTSplitMode2   = iSplitDir==0? 2 : 1;
#endif
  Bool bBTRmvEnable     = iSplitDir==0? bBTHorRmvEnable : bBTVerRmvEnable;
  UInt uiSTSplitMethod  = iSplitDir==0? 1 : 2; //bt split method idx of the sub-cu (0: non-split, 1: HBT, 2: VBT, 3: HTT, 4: VTT)
  UInt uiSTSplitMethod2 = iSplitDir==0? 2 : 1; //the other bt split method idx of the sub-cu, used in redundancy removal

#if DIMD_INTRA_PRED
  UInt uiStartDIMDId = 0, uiEndStartDIMDId = 0;
  Bool bIsBTDepthSigDIMDFlag = false;
  if (isLuma(rpcTempCU->getTextType()))
  {
#if BT_FORCE_DIMD_HARM
    if (rpcTempCU->getBTDepth(0) < uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
    if (rpcTempCU->getBTDepth(0) < DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
    {
      uiStartDIMDId = 0;
      uiEndStartDIMDId = 0;
    }
#if BT_FORCE_DIMD_HARM
    else if (rpcTempCU->getBTDepth(0) == uiDIMDMaxBtDepht && rpcTempCU->getSlice()->isIntra())
#else
    else if (rpcTempCU->getBTDepth(0) == DIMD_MAX_BT_DEPTH_SIG_AI && rpcTempCU->getSlice()->isIntra())
#endif
    {
      uiStartDIMDId = 0;
      uiEndStartDIMDId = 1;
      bIsBTDepthSigDIMDFlag = true;
    }
    else
    {
      assert(rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0));
      uiStartDIMDId = uiEndStartDIMDId = rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0);
    }
  }
  for (UInt uiDIMDModeId = uiStartDIMDId; uiDIMDModeId <= uiEndStartDIMDId; uiDIMDModeId++)
  {
    UInt uiStartDIMDNoBTId = 0, uiEndDIMDNoBTId = 0;
    Bool bIsBTDepthSigDIMDNoBTLevelFlag = false;
    if (isLuma(rpcTempCU->getTextType()) && uiDIMDModeId == 1)
    {
#if BT_FORCE_DIMD_HARM
      assert(rpcTempCU->getBTDepth(0) >= uiDIMDMaxBtDepht && rpcTempCU->getSlice()->isIntra());
      if (rpcTempCU->getBTDepth(0) == uiDIMDMaxBtDepht)
#else
      assert(rpcTempCU->getBTDepth(0) >= DIMD_MAX_BT_DEPTH_SIG_AI && rpcTempCU->getSlice()->isIntra());
      if (rpcTempCU->getBTDepth(0) == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        uiStartDIMDNoBTId = 0;
        uiEndDIMDNoBTId = 1;
        bIsBTDepthSigDIMDNoBTLevelFlag = true;
      }
      else
      {
        assert(rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0));
        uiStartDIMDNoBTId = uiEndDIMDNoBTId = rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0);
      }
    }
    for (UInt uiDIMDNoBTId = uiStartDIMDNoBTId; uiDIMDNoBTId <= uiEndDIMDNoBTId; uiDIMDNoBTId++)
    {
#endif
      rpcTempCU->initEstData( uiDepth, iQP, false );

      UChar uhNextDepth  = uiDepth;
      UInt uiWidthIdx    = g_aucConvertToBit[uiWidth];   
      UInt uiHeightIdx   = g_aucConvertToBit[uiHeight];
      UInt uiBTWidthIdx  = g_aucConvertToBit[uiWidthSubPart];   
      UInt uiBTHeightIdx = g_aucConvertToBit[uiHeightSubPart];
      TComDataCU* pcSubBestPartCU = m_pppcBestCU[uiBTWidthIdx][uiBTHeightIdx];
      TComDataCU* pcSubTempPartCU = m_pppcTempCU[uiBTWidthIdx][uiBTHeightIdx];
#if TT_TREE_DS1
      rpcTempCU->setBTSplitModeSubParts(uiSTSplitMode, 0, uiWidth, uiHeight, uiBTDepth);
#else
      rpcTempCU->setBTSplitModeSubParts(uiSTSplitMode, 0, uiWidth, uiHeight);
#endif
#if JVET_C0024_BT_RMV_REDUNDANT
      uiSplitConstrain = 0;
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( pps.getUseDQP() ) // inherit quantization group info. from parent CU.
      {
        pcSubBestPartCU->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP );           
        pcSubTempPartCU->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP );           // clear sub partition datas or init.
        pcSubBestPartCU->setCodedQP( lastCodedQP );
        pcSubTempPartCU->setCodedQP( lastCodedQP );
        pcSubBestPartCU->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubTempPartCU->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubBestPartCU->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );
        pcSubTempPartCU->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );
      }
#endif

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
      {
        pcSubBestPartCU->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP );           
        pcSubTempPartCU->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP );           // clear sub partition datas or init.

#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
        for( Int size = 0 ; size < NUMBER_OF_PART_SIZES ; size++ )
          m_pppcTempCUIMVCache[size][uiBTWidthIdx][uiBTHeightIdx]->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP ); 
#endif
#if COM16_C806_OBMC
        m_pppcTempCUWoOBMC[uiBTWidthIdx][uiBTHeightIdx]->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP );  // clear sub partition datas or init.
#endif
#if VCEG_AZ07_FRUC_MERGE
        m_pppcFRUCBufferCU[uiBTWidthIdx][uiBTHeightIdx]->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiSTSplitMethod, iQP ); 
#endif
#if BARC_Boundary
        UInt uiPelXInCTUtmp = pcSubBestPartCU->getCUPelX() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUtmp = pcSubBestPartCU->getCUPelY() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUtmp*isBARCFlag;
        UInt deltaY = uiPelYInCTUtmp*isBARCFlag;

        if ((pcSubBestPartCU->getCUPelX() + deltaX < pcSlice->getSPS()->getPicWidthInLumaSamples()) && (pcSubBestPartCU->getCUPelY() + deltaY < pcSlice->getSPS()->getPicHeightInLumaSamples()))
#else
        if(( pcSubBestPartCU->getCUPelX() < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( pcSubBestPartCU->getCUPelY() < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
#endif
        {
          if ( 0 == uiPartUnitIdx) //initialize RD with previous depth buffer
            m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);
          else
            m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_NEXT_BEST]);

#if DIMD_INTRA_PRED
          if (isLuma(rpcTempCU->getTextType()))
          {
#if BT_FORCE_DIMD_HARM
            if (pcSubTempPartCU->getBTDepth(0) > uiDIMDMaxBtDepht && pcSubTempPartCU->getSlice()->isIntra())
#else
            if (pcSubTempPartCU->getBTDepth(0) > DIMD_MAX_BT_DEPTH_SIG_AI && pcSubTempPartCU->getSlice()->isIntra())
#endif
            {
              pcSubBestPartCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDModeId, 0, uiDepth);
              pcSubTempPartCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDModeId, 0, uiDepth);
              pcSubBestPartCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDNoBTId, 0, uiDepth);
              pcSubTempPartCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDNoBTId, 0, uiDepth);
            }
          }
#endif
#if JVET_C0024_BT_RMV_REDUNDANT
          xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uiDepth, uiWidthSubPart, uiHeightSubPart, pcSubBestPartCU->getBTSplitMode(0), uiSplitConstrain );

          if( uiPartUnitIdx == 0 && pcSubBestPartCU->getBTSplitModeForBTDepth(0, uiBTDepth+1) == uiSTSplitMode2 && bBTRmvEnable )
            uiSplitConstrain = uiSTSplitMethod2;
#if TT_RMV_PART1_BTBT && !TT_LEAF_NODE
          if( uiPartUnitIdx == 0 && pcSubBestPartCU->getBTSplitModeForBTDepth(0, uiBTDepth+1) == uiSTSplitMode )
            uiSplitConstrain = uiSTSplitMethod;
#endif
#else
          xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uiDepth, uiWidthSubPart, uiHeightSubPart, pcSubBestPartCU->getBTSplitMode(0) );
#endif
          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidthSubPart, uiHeightSubPart );         // Keep best part data to current temporary data.
#if JVET_C0024_DELTA_QP_FIX
          if( pps.getUseDQP() ) // update coded QP
          { 
            rpcTempCU->setCodedQP( pcSubBestPartCU->getCodedQP() ); 
            pcSubBestPartCU->setCodedQP( rpcTempCU->getCodedQP() );
            pcSubTempPartCU->setCodedQP( rpcTempCU->getCodedQP() );
          }
#endif
#if TT_TREE_DS1
          xCopyYuv2Tmp( pcSubBestPartCU->getZorderIdxInCtu(), rpcTempCU->getZorderIdxInCtu(), uiWidth, uiHeight, uiSTSplitMethod );
#else
          xCopyYuv2Tmp( pcSubBestPartCU->getZorderIdxInCtu()-rpcTempCU->getZorderIdxInCtu(), uiWidth, uiHeight, uiSTSplitMethod );
#endif
        }
#if BT_FORCE
        else
        {
          pcSubBestPartCU->copyToPic(uhNextDepth, uiWidthSubPart, uiHeightSubPart);
          rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidthSubPart, uiHeightSubPart);
          rpcBestCU->getPic()->addCodedAreaInCTU(uiWidthSubPart*uiHeightSubPart);
        }
#endif
      }
      m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_NEXT_BEST]);
        m_pcEntropyCoder->resetBits();
#if BT_FORCE
        if (
#if !BT_FORCE_INTRA_ON
          rpcTempCU->getSlice()->isIntra() ||
#endif
          (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
        {
#endif
        if (uiBTSplitMode == 0)
          m_pcEntropyCoder->encodeSplitFlag(rpcTempCU, 0, uiDepth, true);
        m_pcEntropyCoder->encodeBTSplitMode(rpcTempCU, 0, uiWidth, uiHeight, true);
#if BT_FORCE
        }
#endif
#if DIMD_INTRA_PRED
          if (isLuma(rpcTempCU->getTextType()) && bIsBTDepthSigDIMDFlag)
            m_pcEntropyCoder->encodeDIMDFlag(rpcTempCU, 0, uiDepth, uiWidth, uiHeight, true);
          if (isLuma(rpcTempCU->getTextType()) && bIsBTDepthSigDIMDNoBTLevelFlag)
            m_pcEntropyCoder->encodeDIMDNoBTFlag(rpcTempCU, 0, uiDepth, uiWidth, uiHeight, true);
#endif

        rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
        rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL && SIF_QTBT_FIX
        if (modeSIFnum == 2)
        {
          m_pcEntropyCoder->resetBits();
          m_pcEntropyCoder->encodeAdaptiveIfIndex(rpcTempCU, 0);
          rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits();
          rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
        }
#endif//RRC_SELECTABLE_IF

      rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        Bool foundNonZeroCbf = false;
        UInt uiFirstNonZeroPartIdx = 0;
        rpcTempCU->setQPSubCUs( rpcTempCU->getRefQP( 0 ), 0, uiDepth, uiWidth, uiHeight, uiFirstNonZeroPartIdx, foundNonZeroCbf );
        if ( foundNonZeroCbf )
        {
          m_pcEntropyCoder->resetBits();
          m_pcEntropyCoder->encodeQP( rpcTempCU, uiFirstNonZeroPartIdx, false );
          rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
          rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
          rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
        }
        else
        {
          rpcTempCU->setQPSubParts( rpcTempCU->getRefQP( 0 ), 0, uiWidth, uiHeight ); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          if( pps.getUseDQP() ) // update coded QP
          { 
            rpcTempCU->setCodedQP( rpcTempCU->getQP( 0 ) ); 
          }
#endif
        }
      }
#endif
      m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_TEMP_BEST]);

#if TT_ENC_OPT //store rd cost
      pdParam[0] = rpcTempCU->getTotalCost();
      pdParam[1] = rpcTempCU->getTotalDistortion();
      piParam[5] = rpcTempCU->getTotalBits();
#endif
#if TT_FAST_SMALL_ES || TT_FAST_LARGE_ES
      if(uiBTDepth==0)
      {
        rpcTempCU->setSTSplitModeSave( iSplitDir, rpcTempCU->getBTSplitMode(uiZorderTL), rpcTempCU->getBTSplitMode(uiZorderBR) );
        rpcBestCU->setSTSplitModeSave( iSplitDir, rpcTempCU->getBTSplitMode(uiZorderTL), rpcTempCU->getBTSplitMode(uiZorderBR) );
        rpcTempCU->setSTSplitModeSaveL( iSplitDir, rpcTempCU->getBTSplitMode(uiZorderC1), rpcTempCU->getBTSplitMode(uiZorderC2) );
        rpcBestCU->setSTSplitModeSaveL( iSplitDir, rpcTempCU->getBTSplitMode(uiZorderC1), rpcTempCU->getBTSplitMode(uiZorderC2) );
      }
#endif
#if ENC_FAST_SKIP_RDC_EQUAL //decide second BT
      if (piParam[6] == 0)
      {
        if (rpcTempCU->getTotalDistortion() == rpcBestCU->getTotalDistortion())
          piParam[6] = 1; //true
        else
          piParam[6] = 0; //false
      }
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
      if( rpcTempCU->getTotalCost() < pdParam[2] )
        pdParam[2] = rpcTempCU->getTotalCost();
#endif

      if (rpcBestCU->getTotalCost()!=MAX_DOUBLE)
      {
        const Bool isEndOfSlice        =    pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES
          && ((pcSlice->getSliceBits()+rpcBestCU->getTotalBits())>pcSlice->getSliceArgument()<<3)
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceCurStartCtuTsAddr())
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
        const Bool isEndOfSliceSegment =    pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES
          && ((pcSlice->getSliceSegmentBits()+rpcBestCU->getTotalBits()) > pcSlice->getSliceSegmentArgument()<<3)
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
        // Do not need to check slice condition for slice-segment since a slice-segment is a subset of a slice.
        if(isEndOfSlice||isEndOfSliceSegment)
        {
          rpcBestCU->getTotalCost()=MAX_DOUBLE;
        }
      }

      xCheckBestMode( rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight);
      rpcTempCU->initEstData( uiDepth, iQP, false, uiWidth, uiHeight, uiBTSplitMode );

      rpcBestCU->getPic()->setCodedBlkInCTU(false, uiPelXInCTU>> MIN_CU_LOG2, uiPelYInCTU>> MIN_CU_LOG2, uiWidth>> MIN_CU_LOG2, uiHeight>> MIN_CU_LOG2 );  
      rpcBestCU->getPic()->addCodedAreaInCTU(-(Int)uiWidth*uiHeight);
#if DIMD_INTRA_PRED
    }
  }
#endif
}

#if TT_ENC_OPT
Void TEncCu::CompressTT( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug_), Int* piParam, Double* pdParam, UInt uiSplitConstrain )
{
  TComPic* pcPic = rpcBestCU->getPic();
  DEBUG_STRING_NEW(sDebug)
    const TComPPS &pps=*(rpcTempCU->getSlice()->getPPS());
  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  TComSlice * pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());
#if BT_FORCE_TT
  const UInt uiLPelX = rpcBestCU->getCUPelX();
  const UInt uiRPelX = uiLPelX + rpcBestCU->getWidth(0) - 1;
  const UInt uiTPelY = rpcBestCU->getCUPelY();
  const UInt uiBPelY = uiTPelY + rpcBestCU->getHeight(0) - 1;
#if BARC_BTFORCE
  UInt uiPelXInCTUTmp = uiLPelX - rpcTempCU->getPic()->getCtu(rpcTempCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTUTmp = uiTPelY - rpcTempCU->getPic()->getCtu(rpcTempCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTUTmp + uiWidth)*rpcTempCU->getBARCFlag(rpcTempCU->getTextType());
  UInt deltaY = (uiPelYInCTUTmp + uiHeight)*rpcTempCU->getBARCFlag(rpcTempCU->getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#else
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#endif
#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(rpcBestCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
#endif
#if JVET_C0024_DELTA_QP_FIX
  UInt uiQTWidth  = sps.getCTUSize()>>uiDepth;
  UInt uiQTHeight = sps.getCTUSize()>>uiDepth;
#if TT_TREE_DS1
  UInt uiBTDepth  = rpcTempCU->getBTDepthFromBTSplitMode(uiBTSplitMode);
  const UInt uiQTBTDepth = (uiDepth<<1) + (g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight]);
#else
  UInt uiBTDepth  = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight];
  const UInt uiQTBTDepth = (uiDepth<<1) + uiBTDepth;
#endif
  const UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif
#if BARC
  Bool isBARCFlag = rpcBestCU->getBARCFlag(rpcBestCU->getTextType());
#endif
  Int iC1x = uiWidth >> (MIN_CU_LOG2+2);
  Int iC1y = 0;
  Int iC2x = (uiWidth >> (MIN_CU_LOG2+1)) + (uiWidth >> (MIN_CU_LOG2+2)) - 1;
  Int iC2y = (uiHeight >> MIN_CU_LOG2) - 1;
  iC2x = iC2x<0? 0 : iC2x;
  iC2y = iC2y<0? 0 : iC2y;
  UInt uiZorderTL = 0; //top left
  UInt uiZorderC1 = g_auiRasterToZscan[iC1y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC1x]; //center, closer to TL
  UInt uiZorderC2 = g_auiRasterToZscan[iC2y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC2x]; //center, closer to BR
  UInt uiZorderBR = g_auiRasterToZscan[((uiHeight>> MIN_CU_LOG2)-1) * (sps.getCTUSize()>> MIN_CU_LOG2) + (uiWidth>> MIN_CU_LOG2)-1];   //bottom-right part.
  UInt uiPelXInCTU = rpcBestCU->getCUPelX() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = rpcBestCU->getCUPelY() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelY();

  //parameter reading
  Int  iSplitDir       = piParam[0];
  Int  iQP             = piParam[1];
  Int  lastCodedQP     = piParam[2];
  piParam[6] = 0;

#if TT_BUGFIX_ZIDX
  if(iSplitDir==0)
  {
    iC1x = 0;
    iC1y = uiHeight >> (MIN_CU_LOG2+2);
    iC2x = (uiWidth >> MIN_CU_LOG2) - 1;
    iC2y = (uiHeight >> (MIN_CU_LOG2+1)) + (uiHeight >> (MIN_CU_LOG2+2)) - 1;
    iC2x = iC2x<0? 0 : iC2x;
    iC2y = iC2y<0? 0 : iC2y;
    uiZorderC1 = g_auiRasterToZscan[iC1y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC1x]; //center, closer to TL
    uiZorderC2 = g_auiRasterToZscan[iC2y * (sps.getCTUSize()>> MIN_CU_LOG2) + iC2x]; //center, closer to BR
  }
#endif

  //initial parameter related to split dir
  UInt uiWidthSubPartL  = iSplitDir==1? (uiWidth>>1)  : uiWidth;
  UInt uiWidthSubPartS  = iSplitDir==1? (uiWidth>>2)  : uiWidth;
  UInt uiHeightSubPartL = iSplitDir==0? (uiHeight>>1) : uiHeight;
  UInt uiHeightSubPartS = iSplitDir==0? (uiHeight>>2) : uiHeight;
  UInt uiSTSplitMode    = iSplitDir==0? ST_HORI_TT_L  : ST_VERT_TT_L;  //tt split mode of the large sub-cu
  UInt uiSTSplitMethodL = iSplitDir==0? 3             : 4;             //tt split mode idx of the large sub-cu
  UInt uiSTSplitMethodS = iSplitDir==0? 5             : 6;             //tt split mode idx of the small sub-cu

#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(rpcBestCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
  const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
  UInt uiStartDIMDId = 0, uiEndStartDIMDId = 0;
  Bool bIsBTDepthSigDIMDFlag = false;
  if (isLuma(rpcTempCU->getTextType()))
  {
#if BT_FORCE_DIMD_HARM
    if (rpcTempCU->getBTDepth(0) < uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
    if (rpcTempCU->getBTDepth(0) < DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
    {
      uiStartDIMDId = 0;
      uiEndStartDIMDId = 0;
    }
#if BT_FORCE_DIMD_HARM
    else if (rpcTempCU->getBTDepth(0) == uiDIMDMaxBtDepht && rpcTempCU->getSlice()->isIntra())
#else
    else if (rpcTempCU->getBTDepth(0) == DIMD_MAX_BT_DEPTH_SIG_AI && rpcTempCU->getSlice()->isIntra())
#endif
    {
      uiStartDIMDId = 0;
      uiEndStartDIMDId = 1;
      bIsBTDepthSigDIMDFlag = true;
    }
    else
    {
      assert(rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0));
      uiStartDIMDId = uiEndStartDIMDId = rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0);
    }
  }
  for (UInt uiDIMDModeId = uiStartDIMDId; uiDIMDModeId <= uiEndStartDIMDId; uiDIMDModeId++)
  {
    UInt uiStartDIMDNoBTId = 0, uiEndDIMDNoBTId = 0;
    Bool bIsBTDepthSigDIMDNoBTLevelFlag = false;
    if (isLuma(rpcTempCU->getTextType()) && uiDIMDModeId == 1)
    {
#if BT_FORCE_DIMD_HARM
      assert(rpcTempCU->getBTDepth(0) >= uiDIMDMaxBtDepht && rpcTempCU->getSlice()->isIntra());
      if (rpcTempCU->getBTDepth(0) == uiDIMDMaxBtDepht)
#else
      assert(rpcTempCU->getBTDepth(0) >= DIMD_MAX_BT_DEPTH_SIG_AI && rpcTempCU->getSlice()->isIntra());
      if (rpcTempCU->getBTDepth(0) == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        uiStartDIMDNoBTId = 0;
        uiEndDIMDNoBTId = 1;
        bIsBTDepthSigDIMDNoBTLevelFlag = true;
      }
      else
      {
        assert(rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0));
        uiStartDIMDNoBTId = uiEndDIMDNoBTId = rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0);
      }
    }
    for (UInt uiDIMDNoBTId = uiStartDIMDNoBTId; uiDIMDNoBTId <= uiEndDIMDNoBTId; uiDIMDNoBTId++)
    {
#endif
      rpcTempCU->initEstData( uiDepth, iQP, false );

      UChar uhNextDepth  = uiDepth;
      UInt uiWidthIdx    = g_aucConvertToBit[uiWidth];   
      UInt uiHeightIdx   = g_aucConvertToBit[uiHeight];
      UInt uiBTWidthIdxL = g_aucConvertToBit[uiWidthSubPartL];   
      UInt uiBTWidthIdxS = g_aucConvertToBit[uiWidthSubPartS];   
      UInt uiBTHeightIdxL = g_aucConvertToBit[uiHeightSubPartL];
      UInt uiBTHeightIdxS = g_aucConvertToBit[uiHeightSubPartS];
      TComDataCU* pcSubBestPartCUL     = m_pppcBestCU[uiBTWidthIdxL][uiBTHeightIdxL];
      TComDataCU* pcSubTempPartCUL     = m_pppcTempCU[uiBTWidthIdxL][uiBTHeightIdxL];
      TComDataCU* pcSubBestPartCUS     = m_pppcBestCU[uiBTWidthIdxS][uiBTHeightIdxS];
      TComDataCU* pcSubTempPartCUS     = m_pppcTempCU[uiBTWidthIdxS][uiBTHeightIdxS];
      UInt uiBTWidthIdx, uiBTHeightIdx, uiBTSPlitModeSubPart;
      TComDataCU* pcSubBestPartCU;
      TComDataCU* pcSubTempPartCU;
      UInt uiWidthSubPart, uiHeightSubPart;

      rpcTempCU->setBTSplitModeSubParts(uiSTSplitMode, 0, uiWidth, uiHeight, uiBTDepth); //either ST_HORI_TT_L or ST_HORI_TT_S is OK
#if JVET_C0024_DELTA_QP_FIX
      if( pps.getUseDQP() ) // inherit quantization group info. from parent CU.
      {
        pcSubBestPartCUL->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPartL, uiHeightSubPartL, uiSTSplitMethodL, iQP );           
        pcSubTempPartCUL->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPartL, uiHeightSubPartL, uiSTSplitMethodL, iQP );           // clear sub partition datas or init.
        pcSubBestPartCUL->setCodedQP( lastCodedQP );
        pcSubTempPartCUL->setCodedQP( lastCodedQP );
        pcSubBestPartCUL->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubTempPartCUL->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubBestPartCUL->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );
        pcSubTempPartCUL->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );

        pcSubBestPartCUS->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPartS, uiHeightSubPartS, uiSTSplitMethodS, iQP );           
        pcSubTempPartCUS->initSubBT( rpcTempCU, 0, uiDepth, uiWidthSubPartS, uiHeightSubPartS, uiSTSplitMethodS, iQP );           // clear sub partition datas or init.
        pcSubBestPartCUS->setCodedQP( lastCodedQP );
        pcSubTempPartCUS->setCodedQP( lastCodedQP );
        pcSubBestPartCUS->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubTempPartCUS->setQuPartIdx( rpcTempCU->getQuPartIdx() );
        pcSubBestPartCUS->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );
        pcSubTempPartCUS->setQuLastCodedQP( rpcTempCU->getQuLastCodedQP() );
      }
#endif

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
      {
        pcSubBestPartCU = uiPartUnitIdx==1? pcSubBestPartCUL : pcSubBestPartCUS;
        pcSubTempPartCU = uiPartUnitIdx==1? pcSubTempPartCUL : pcSubTempPartCUS;
        uiWidthSubPart  = uiPartUnitIdx==1? uiWidthSubPartL  : uiWidthSubPartS;
        uiHeightSubPart = uiPartUnitIdx==1? uiHeightSubPartL : uiHeightSubPartS;
        uiBTWidthIdx    = uiPartUnitIdx==1? uiBTWidthIdxL    : uiBTWidthIdxS;
        uiBTHeightIdx   = uiPartUnitIdx==1? uiBTHeightIdxL   : uiBTHeightIdxS;
        uiBTSPlitModeSubPart = uiPartUnitIdx==1? uiSTSplitMethodL : uiSTSplitMethodS;

        pcSubBestPartCU->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiBTSPlitModeSubPart, iQP );           
        pcSubTempPartCU->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiBTSPlitModeSubPart, iQP );           // clear sub partition datas or init.
#if COM16_C806_OBMC
        m_pppcTempCUWoOBMC[uiBTWidthIdx][uiBTHeightIdx]->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiBTSPlitModeSubPart, iQP );  // clear sub partition datas or init.
#endif
#if VCEG_AZ07_FRUC_MERGE
        m_pppcFRUCBufferCU[uiBTWidthIdx][uiBTHeightIdx]->initSubBT( rpcTempCU, uiPartUnitIdx, uiDepth, uiWidthSubPart, uiHeightSubPart, uiBTSPlitModeSubPart, iQP ); 
#endif
#if BARC_Boundary
        UInt uiPelXInCTUtmp = pcSubBestPartCU->getCUPelX() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUtmp = pcSubBestPartCU->getCUPelY() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUtmp*isBARCFlag;
        UInt deltaY = uiPelYInCTUtmp*isBARCFlag;

        if ((pcSubBestPartCU->getCUPelX() + deltaX < pcSlice->getSPS()->getPicWidthInLumaSamples()) && (pcSubBestPartCU->getCUPelY() + deltaY < pcSlice->getSPS()->getPicHeightInLumaSamples()))
#else
        if(( pcSubBestPartCU->getCUPelX() < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( pcSubBestPartCU->getCUPelY() < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
#endif
        {
          if ( 0 == uiPartUnitIdx) //initialize RD with previous depth buffer
            m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);
          else if ( 1 == uiPartUnitIdx) //copy from part 0
            m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiBTWidthIdxS][uiBTHeightIdxS][CI_NEXT_BEST]);
          else //copy from part 1
            m_ppppcRDSbacCoder[uiBTWidthIdx][uiBTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiBTWidthIdxL][uiBTHeightIdxL][CI_NEXT_BEST]);

#if DIMD_INTRA_PRED
          if (isLuma(rpcTempCU->getTextType()))
          {
#if BT_FORCE_DIMD_HARM
            if (pcSubTempPartCU->getBTDepth(0) > uiDIMDMaxBtDepht && pcSubTempPartCU->getSlice()->isIntra())
#else
            if (pcSubTempPartCU->getBTDepth(0) > DIMD_MAX_BT_DEPTH_SIG_AI && pcSubTempPartCU->getSlice()->isIntra())
#endif
            {
              pcSubBestPartCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDModeId, 0, uiDepth);
              pcSubTempPartCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDModeId, 0, uiDepth);
              pcSubBestPartCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDNoBTId, 0, uiDepth);
              pcSubTempPartCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_LUMA, uiDIMDNoBTId, 0, uiDepth);
            }
          }
#endif

          xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uiDepth, uiWidthSubPart, uiHeightSubPart, pcSubBestPartCU->getBTSplitMode(0) );

          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidthSubPart, uiHeightSubPart, true );         // Keep best part data to current temporary data.
#if JVET_C0024_DELTA_QP_FIX
          if( pps.getUseDQP() ) // update coded QP
          { 
            rpcTempCU->setCodedQP( pcSubBestPartCU->getCodedQP() ); 
            pcSubBestPartCU->setCodedQP( rpcTempCU->getCodedQP() );
            pcSubTempPartCU->setCodedQP( rpcTempCU->getCodedQP() );
          }
#endif
#if TT_TREE_DS1
          xCopyYuv2Tmp( pcSubBestPartCU->getZorderIdxInCtu(), rpcTempCU->getZorderIdxInCtu(), uiWidth, uiHeight, uiBTSPlitModeSubPart );
#else
          xCopyYuv2Tmp( pcSubBestPartCU->getZorderIdxInCtu()-rpcTempCU->getZorderIdxInCtu(), uiWidth, uiHeight, uiBTSPlitModeSubPart );
#endif
        }
#if BT_FORCE_TT
        else
        {
          pcSubBestPartCU->copyToPic(uhNextDepth, uiWidthSubPart, uiHeightSubPart);
          rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidthSubPart, uiHeightSubPart, true);
          rpcBestCU->getPic()->addCodedAreaInCTU(uiWidthSubPart*uiHeightSubPart);
        }
#endif
      }
      m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiBTWidthIdxS][uiBTHeightIdxS][CI_NEXT_BEST]);
      m_pcEntropyCoder->resetBits();
#if BT_FORCE_TT
      if (
#if !BT_FORCE_INTRA_ON
        rpcTempCU->getSlice()->isIntra() ||
#endif
        ( !(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
      {
#endif
        if (uiBTSplitMode == 0)
        {
          m_pcEntropyCoder->encodeSplitFlag(rpcTempCU, 0, uiDepth, true);
        }
        m_pcEntropyCoder->encodeBTSplitMode(rpcTempCU, 0, uiWidth, uiHeight, true);
#if BT_FORCE_TT
      }
#endif
#if DIMD_INTRA_PRED
      if (isLuma(rpcTempCU->getTextType()) && bIsBTDepthSigDIMDFlag)
        m_pcEntropyCoder->encodeDIMDFlag(rpcTempCU, 0, uiDepth, uiWidth, uiHeight, true);
      if (isLuma(rpcTempCU->getTextType()) && bIsBTDepthSigDIMDNoBTLevelFlag)
        m_pcEntropyCoder->encodeDIMDNoBTFlag(rpcTempCU, 0, uiDepth, uiWidth, uiHeight, true);
#endif
      rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
      rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        Bool foundNonZeroCbf = false;
        UInt uiFirstNonZeroPartIdx = 0;
        rpcTempCU->setQPSubCUs( rpcTempCU->getRefQP( 0 ), 0, uiDepth, uiWidth, uiHeight, uiFirstNonZeroPartIdx, foundNonZeroCbf );
        if ( foundNonZeroCbf )
        {
          m_pcEntropyCoder->resetBits();
          m_pcEntropyCoder->encodeQP( rpcTempCU, uiFirstNonZeroPartIdx, false );
          rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
          rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
          rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
        }
        else
        {
          rpcTempCU->setQPSubParts( rpcTempCU->getRefQP( 0 ), 0, uiWidth, uiHeight ); // set QP to default QP
#if JVET_C0024_DELTA_QP_FIX
          if( pps.getUseDQP() ) // update coded QP
          { 
            rpcTempCU->setCodedQP( rpcTempCU->getQP( 0 ) ); 
          }
#endif
        }
      }
#endif
      m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_TEMP_BEST]);

#if TT_ENC_OPT //store rd cost
      pdParam[0] = rpcTempCU->getTotalCost();
      pdParam[1] = rpcTempCU->getTotalDistortion();
      piParam[5] = rpcTempCU->getTotalBits();
#endif
#if ENC_FAST_SKIP_RDC_EQUAL //decide second TT
      if (piParam[6]==0)

      {
        if (rpcTempCU->getTotalDistortion() == rpcBestCU->getTotalDistortion())
          piParam[6] = 1; //true
        else
          piParam[6] = 0; //false
      }
#endif

      if (rpcBestCU->getTotalCost()!=MAX_DOUBLE)
      {
        const Bool isEndOfSlice        =    pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES
          && ((pcSlice->getSliceBits()+rpcBestCU->getTotalBits())>pcSlice->getSliceArgument()<<3)
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceCurStartCtuTsAddr())
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
        const Bool isEndOfSliceSegment =    pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES
          && ((pcSlice->getSliceSegmentBits()+rpcBestCU->getTotalBits()) > pcSlice->getSliceSegmentArgument()<<3)
          && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
        // Do not need to check slice condition for slice-segment since a slice-segment is a subset of a slice.
        if(isEndOfSlice||isEndOfSliceSegment)
        {
          rpcBestCU->getTotalCost()=MAX_DOUBLE;
        }
      }

      xCheckBestMode( rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight);
      rpcTempCU->initEstData( uiDepth, iQP, false, uiWidth, uiHeight, uiBTSplitMode );

      rpcBestCU->getPic()->setCodedBlkInCTU(false, uiPelXInCTU>> MIN_CU_LOG2, uiPelYInCTU>> MIN_CU_LOG2, uiWidth>> MIN_CU_LOG2, uiHeight>> MIN_CU_LOG2 );  
      rpcBestCU->getPic()->addCodedAreaInCTU(-(Int)uiWidth*uiHeight);
#if DIMD_INTRA_PRED
    }
  }
#endif
}
#endif
#endif


#if AMP_ENC_SPEEDUP
#if JVET_C0024_QTBT
Void TEncCu::xCompressCU(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiBTSplitMode DEBUG_STRING_FN_DECLARE(sDebug_), UInt uiSplitConstrain)
#else
Void TEncCu::xCompressCU(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth DEBUG_STRING_FN_DECLARE(sDebug_), PartSize eParentPartSize)
#endif
#else
Void TEncCu::xCompressCU(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth)
#endif
{
#if JVET_C0024_BT_RMV_REDUNDANT
  rpcBestCU->setSplitConstrain(uiSplitConstrain);
  rpcTempCU->setSplitConstrain(uiSplitConstrain);

  Bool bQTreeValid = false;
#endif
  TComPic* pcPic = rpcBestCU->getPic();
  DEBUG_STRING_NEW(sDebug)
    const TComPPS &pps = *(rpcTempCU->getSlice()->getPPS());
  const TComSPS &sps = *(rpcTempCU->getSlice()->getSPS());

#if JVET_C0024_QTBT
  // These are only used if getFastDeltaQp() is true
  const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMinQTSize(rpcBestCU->getSlice()->getSliceType(), rpcBestCU->getTextType())
    , sps.getCTUSize(), 32u);
#else
  // These are only used if getFastDeltaQp() is true
  const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight() >> sps.getLog2DiffMaxMinCodingBlockSize(), sps.getMaxCUHeight(), 32u);
#endif

  // get Original YUV data from picture
#if JVET_C0024_QTBT
  UInt uiWidthIdx = g_aucConvertToBit[rpcTempCU->getWidth(0)];
  UInt uiHeightIdx = g_aucConvertToBit[rpcTempCU->getHeight(0)];

#if BARC
  Bool isBARCFlag = rpcBestCU->getBARCFlag(rpcBestCU->getTextType());
  if (isBARCFlag)
    m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->copyFromPicYuv( pcPic->getPicYuvDScaleOrg(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu(), true);
  else
#endif
  m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->copyFromPicYuv(pcPic->getPicYuvOrg(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu());

  UInt uiPelXInCTU = rpcBestCU->getCUPelX() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = rpcBestCU->getCUPelY() - rpcBestCU->getPic()->getCtu(rpcBestCU->getCtuRsAddr())->getCUPelY();
  rpcBestCU->getPic()->setCodedBlkInCTU(false, uiPelXInCTU >> MIN_CU_LOG2, uiPelYInCTU >> MIN_CU_LOG2, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
#else
  m_ppcOrigYuv[uiDepth]->copyFromPicYuv(pcPic->getPicYuvOrg(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu());
#endif

  // variable for Cbf fast mode PU decision
#if !JVET_C0024_QTBT
  Bool    doNotBlockPu = true;
#endif
  Bool    earlyDetectionSkipMode = false;

  const UInt uiLPelX = rpcBestCU->getCUPelX();
  const UInt uiRPelX = uiLPelX + rpcBestCU->getWidth(0) - 1;
  const UInt uiTPelY = rpcBestCU->getCUPelY();
  const UInt uiBPelY = uiTPelY + rpcBestCU->getHeight(0) - 1;
#if !JVET_C0024_QTBT
  const UInt uiWidth = rpcBestCU->getWidth(0);
#endif
#if JVET_C0024_DELTA_QP_FIX
  UInt uiQTWidth = sps.getCTUSize() >> uiDepth;
  UInt uiQTHeight = sps.getCTUSize() >> uiDepth;
#if TT_TREE_DS1
  UInt uiBTDepth = rpcTempCU->getBTDepthFromBTSplitMode(uiBTSplitMode);
  const UInt uiQTBTDepth = (uiDepth << 1) + (g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight]);
#else
  UInt uiBTDepth = g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight];
  const UInt uiQTBTDepth = (uiDepth << 1) + uiBTDepth;
#endif
  const UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  Bool bUseSaveLoad = m_pcEncCfg->getUseSaveLoadEncInfo();
#else
  Bool bUseSaveLoad = m_pcEncCfg->getUseSaveLoadEncInfo() && uiWidthIdx > 0 && uiHeightIdx > 0;
#endif
  Bool bUseSaveLoadSplitDecision = bUseSaveLoad && m_pcEncCfg->getUseSaveLoadSplitDecision();
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  ChannelType eChannelType = rpcBestCU->getTextType();
#endif
  UInt uiZorderIdx = rpcBestCU->getZorderIdxInCtu();
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT //search
  UInt uiQTLeafWidth = (rpcBestCU->getSlice()->getSPS()->getCTUSize() >> (uiDepth + MIN_CU_LOG2));
  UInt uiCTUSize = rpcBestCU->getSlice()->getSPS()->getCTUSize();
  UInt uiNumPartInQTLeaf = uiQTLeafWidth * uiQTLeafWidth;

  if (uiBTDepth == 0) //init at the QTLeafNode
  {
    for (UInt uiWIdx = 0; uiWIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiWIdx++)
    {
      for (UInt uiHIdx = 0; uiHIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiHIdx++)
      {
        for (UInt uiPartIdx = uiZorderIdx; uiPartIdx < uiZorderIdx + uiNumPartInQTLeaf; uiPartIdx++)
          m_pcPredSearch->cleanSaveLoadInfo(uiPartIdx, uiWIdx, uiHIdx);
      }
    }
  }
  //match CU Route Info in the database
  UInt uiSaveLoadIdx = m_pcPredSearch->searchSaveLoadIdx(uiZorderIdx, uiWidthIdx, uiHeightIdx, uiBTDepth, rpcBestCU->getCURouteInfo());
  rpcBestCU->setCURouteIdx(uiSaveLoadIdx);
  rpcTempCU->setCURouteIdx(uiSaveLoadIdx);
#endif
#if ENC_FAST_SKIP_RDC_EQUAL //store rd cost
  double dRDCost[6] = { MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE }; //rd cost for non-split, HBT, VBT, HTT, VTT, QT
  double dRDDist[6] = { MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE,MAX_DOUBLE }; //distortion for non-split, HBT, VBT, HTT, VTT, QT
  int    iRDBits[6] = { MAX_INT,   MAX_INT,   MAX_INT,   MAX_INT,   MAX_INT,  MAX_INT };    //bits for non-split, HBT, VBT, HTT, VTT, QT
#endif


  Int iBaseQP = xComputeQP(rpcBestCU, uiDepth);
  Int iMinQP;
  Int iMaxQP;
  Bool isAddLowestQP = false;

#if COM16_C806_EMT
  Double dBestInterCost = MAX_DOUBLE;
  Bool   bEarlySkipIntra = false;
#if !JVET_C0024_QTBT
  Bool   bAllIntra = (m_pcEncCfg->getIntraPeriod() == 1);
  Double dIntra2Nx2NCost = MAX_DOUBLE;
#endif
#if !JVET_C0024_QTBT
  Double dIntraNxNCost = MAX_DOUBLE;
#endif
#endif

#if JVET_C0024_DELTA_QP_FIX
  const Char lastCodedQP = rpcBestCU->getCodedQP();
#endif

#if JVET_E0076_MULTI_PEL_MVD
  m_dBestMvDPelCost[0] = m_dBestMvDPelCost[1] = m_dBestMvDPelCost[2] = MAX_DOUBLE / 2;
#endif

  const UInt numberValidComponents = rpcBestCU->getPic()->getNumberValidComponents();
#if VCEG_AZ08_INTER_KLT
  g_bEnableCheck = false;
#endif
#if JVET_C0024_DELTA_QP_FIX
  if (uiQTBTDepth <= uiMaxDQPDepthQTBT)
#else
  if (uiDepth <= pps.getMaxCuDQPDepth())
#endif
  {
    Int idQP = m_pcEncCfg->getMaxDeltaQP();
    iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - idQP);
    iMaxQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP + idQP);
#if JVET_C0024_DELTA_QP_FIX
    if (pps.getUseDQP())
    {
      UInt uiCurrPartIdxInCtu = rpcBestCU->getZorderIdxInCtu();
      rpcBestCU->setQuPartIdx(uiCurrPartIdxInCtu);
      rpcTempCU->setQuPartIdx(uiCurrPartIdxInCtu);
      rpcBestCU->setQuLastCodedQP(lastCodedQP);
      rpcTempCU->setQuLastCodedQP(lastCodedQP);
    }
#endif
  }
  else
  {
    iMinQP = rpcTempCU->getQP(0);
    iMaxQP = rpcTempCU->getQP(0);
  }
#if WCG_LUMA_DQP_CM_SCALE
  Int targetQP = iBaseQP;
  if (m_pcEncCfg->getUseLumaDeltaQp())
  {
    if (uiQTBTDepth <= uiMaxDQPDepthQTBT)

      m_LumaQPOffset = calculateLumaDQP(rpcTempCU, 0, m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]);  // keep using the same m_QP_LUMA_OFFSET in the same LCU
		targetQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - m_LumaQPOffset);  // targetQP is used for control lambda
		{
      iMinQP = targetQP;
      iMaxQP = iMinQP;   // make it same as MinQP to force encode choose the modified QP
    }
  }

#endif
  if (m_pcEncCfg->getUseRateCtrl())
  {
    iMinQP = m_pcRateCtrl->getRCQP();
    iMaxQP = m_pcRateCtrl->getRCQP();
  }

  // transquant-bypass (TQB) processing loop variable initialisation ---

  const Int lowestQP = iMinQP; // For TQB, use this QP which is the lowest non TQB QP tested (rather than QP'=0) - that way delta QPs are smaller, and TQB can be tested at all CU levels.

  if ((pps.getTransquantBypassEnableFlag()))
  {
    isAddLowestQP = true; // mark that the first iteration is to cost TQB mode.
    iMinQP = iMinQP - 1;  // increase loop variable range by 1, to allow testing of TQB mode along with other QPs
    if (m_pcEncCfg->getCUTransquantBypassFlagForceValue())
    {
      iMaxQP = iMinQP;
    }
  }

#if VCEG_AZ06_IC
  Bool bICEnabled = rpcTempCU->getSlice()->getApplyIC();
#if JVET_C0024_QTBT
#if VCEG_AZ06_IC_SPEEDUP
  if (uiWidth * uiHeight <= 32)
  {
    bICEnabled = false;
  }
#endif
#endif
#endif

  TComSlice * pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());

#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
  UChar ucMinDepth = 0, ucMaxDepth = g_aucConvertToBit[pcSlice->getSPS()->getCTUSize()]
    - g_aucConvertToBit[pcSlice->getSPS()->getMinQTSize(pcSlice->getSliceType(), rpcTempCU->getTextType())];
#else
  UChar ucMinDepth = 0, ucMaxDepth = (UChar)pcSlice->getSPS()->getLog2DiffMaxMinCodingBlockSize();
#endif
  if (m_pcEncCfg->getUseFastLCTU())
  {
    rpcTempCU->getMaxMinCUDepth(ucMinDepth, ucMaxDepth, rpcTempCU->getZorderIdxInCtu());
  }
#endif


#if JVET_C0024_QTBT
  Bool bBoundary = !(uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples());
#if BARC_Boundary 
  if (isBARCFlag)
  {
    UInt deltaX =  uiPelXInCTU + rpcBestCU->getWidth(0);
    UInt deltaY =  uiPelYInCTU + rpcBestCU->getHeight(0);

    bBoundary = !( uiRPelX < sps.getPicWidthInLumaSamples() - deltaX && uiBPelY < sps.getPicHeightInLumaSamples()-deltaY );
  }
#endif
#if BT_FORCE
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#if BARC_BTFORCE
  if (isBARCFlag)
  {
     UInt deltaX =  uiPelXInCTU + rpcBestCU->getWidth(0);
     UInt deltaY =  uiPelYInCTU + rpcBestCU->getHeight(0);

     bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
     bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
  }
#endif
#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(rpcTempCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
#endif
#else
  const Bool bBoundary = !(uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples());
#endif
#if JVET_C0024_QTBT
  Bool bForceQT = uiWidth > MAX_TU_SIZE;
  if (bForceQT)
  {
    assert(uiWidth == uiHeight);
  }
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  UChar saveLoadTag = m_pcPredSearch->getSaveLoadTag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx);
  UChar saveLoadSplit = ((bUseSaveLoadSplitDecision && saveLoadTag == LOAD_ENC_INFO) ? m_pcPredSearch->getSaveLoadSplit(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx) : 0);
#else
  UChar saveLoadTag = m_pcPredSearch->getSaveLoadTag(uiZorderIdx, uiWidthIdx, uiHeightIdx);
  UChar saveLoadSplit = ((bUseSaveLoadSplitDecision && saveLoadTag == LOAD_ENC_INFO) ? m_pcPredSearch->getSaveLoadSplit(uiWidthIdx, uiHeightIdx) : 0);
#endif
  Double dCostTempBest = MAX_DOUBLE;
  Double dNonSplitCost = MAX_DOUBLE;
  Double dHorSplitCost = MAX_DOUBLE;
  Double dVerSplitCost = MAX_DOUBLE;
#endif
  if (!bBoundary
#if COM16_C806_LARGE_CTU
    && ucMinDepth <= uiDepth
#endif
#if JVET_C0024_QTBT
    && !bForceQT
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
    && !(saveLoadSplit & 0x01)
#endif
    )
  {
#if JVET_D0077_FAST_EXT
    Bool bPrevSameBlockIsIntra = rpcBestCU->getPic()->getIntra(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight);
    Bool bPrevSameBlockIsSkip = rpcBestCU->getPic()->getSkiped(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight);
#endif
#if BARC_CORE & BARC_AdpQP
    if (isBARCFlag && (iMinQP == iBaseQP))
    {
      if (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA)
      {
        iMinQP -= BARC_AdpQPSTEP;
        iMaxQP -= BARC_AdpQPSTEP;
      }
    }
#endif
    for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
    {
#if RRC_SELECTABLE_IF
      Int modeSIFnum = 1;
      if (rpcTempCU->isCodeBIF(uiWidth, uiDepth
#if JVET_C0024_QTBT && SIF_QTBT_FIX
        + uiBTDepth // check total depth
#endif
        )
#if SIF_SPLIT_ONLY
        && 0
#endif
#if SIF_CU_LEVEL
        || 1
#endif
        )
      {
        modeSIFnum = 2;
      }
      for (int modeSIF = 0; modeSIF < modeSIFnum; modeSIF++)
      {
        if (modeSIFnum == 2)
        {
          rpcTempCU->setSIFToCheck(modeSIF);
          rpcBestCU->setSIFToCheck(modeSIF);
        }
#endif
        const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

        if (bIsLosslessMode)
        {
          iQP = lowestQP;
        }
#if WCG_LUMA_DQP_CM_SCALE
        if (m_pcEncCfg->getUseLumaDeltaQp() && (uiQTBTDepth <= uiMaxDQPDepthQTBT))
        {
          getSliceEncoder()->updateLambda(pcSlice, targetQP);
        }
#endif

        m_cuChromaQpOffsetIdxPlus1 = 0;
        if (pcSlice->getUseChromaQpAdj())
        {
          /* Pre-estimation of chroma QP based on input block activity may be performed
           * here, using for example m_ppcOrigYuv[uiDepth] */
           /* To exercise the current code, the index used for adjustment is based on
            * block position
            */
#if JVET_C0024_QTBT
          Int lgMinCuSize = g_aucConvertToBit[sps.getMinQTSize(pcSlice->getSliceType(), rpcBestCU->getTextType())] + MIN_CU_LOG2 +
            std::max<Int>(0, g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[sps.getMinQTSize(pcSlice->getSliceType(), rpcBestCU->getTextType())] - Int(pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth()));
#else
          Int lgMinCuSize = sps.getLog2MinCodingBlockSize() +
            std::max<Int>(0, sps.getLog2DiffMaxMinCodingBlockSize() - Int(pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth()));
#endif
          m_cuChromaQpOffsetIdxPlus1 = ((uiLPelX >> lgMinCuSize) + (uiTPelY >> lgMinCuSize)) % (pps.getPpsRangeExtension().getChromaQpOffsetListLen() + 1);
        }

        rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
        m_setInterCand.clear();
        for (Int n = 0; n < NUMBER_OF_PART_SIZES; n++)
        {
#if JVET_C0024_QTBT
          m_pppcTempCUIMVCache[n][uiWidthIdx][uiHeightIdx]->initEstData(uiDepth, iQP, bIsLosslessMode);
#else
          m_ppcTempCUIMVCache[n][uiDepth]->initEstData(uiDepth, iQP, bIsLosslessMode);
#endif
        }
#endif

        // do inter modes, SKIP and 2Nx2N
        if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
        {
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT 
#if VCEG_AZ07_FRUC_MERGE
          Bool bFastSkipFruc = (saveLoadTag == LOAD_ENC_INFO && 0 == m_pcPredSearch->getSaveLoadFrucMode(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx));
#endif
#if VCEG_AZ07_IMV
          Bool bFastSkipIMV = (saveLoadTag == LOAD_ENC_INFO && (!m_pcPredSearch->getSaveLoadIMVFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx) || m_pcPredSearch->getSaveLoadMergeFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx)));
#endif
          Bool bFastSkipInter = (saveLoadTag == LOAD_ENC_INFO && m_pcPredSearch->getSaveLoadMergeFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx));
#if COM16_C1016_AFFINE
          Bool bFastSkipAffine = (saveLoadTag == LOAD_ENC_INFO && !m_pcPredSearch->getSaveLoadAffineFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx));
#endif
#else
#if VCEG_AZ07_FRUC_MERGE
          Bool bFastSkipFruc = (saveLoadTag == LOAD_ENC_INFO && 0 == m_pcPredSearch->getSaveLoadFrucMode(uiWidthIdx, uiHeightIdx));
#endif
#if VCEG_AZ07_IMV
          Bool bFastSkipIMV = (saveLoadTag == LOAD_ENC_INFO && (!m_pcPredSearch->getSaveLoadIMVFlag(uiWidthIdx, uiHeightIdx) || m_pcPredSearch->getSaveLoadMergeFlag(uiWidthIdx, uiHeightIdx)));
#endif
          Bool bFastSkipInter = (saveLoadTag == LOAD_ENC_INFO && m_pcPredSearch->getSaveLoadMergeFlag(uiWidthIdx, uiHeightIdx));
#if COM16_C1016_AFFINE
          Bool bFastSkipAffine = (saveLoadTag == LOAD_ENC_INFO && !m_pcPredSearch->getSaveLoadAffineFlag(uiWidthIdx, uiHeightIdx));
#endif
#endif
#endif
#if VCEG_AZ06_IC
          for (UInt uiICId = 0; uiICId < (bICEnabled ? 2 : 1); uiICId++)
          {
            Bool bICFlag = uiICId ? true : false;
#endif
            // 2Nx2N
            if (m_pcEncCfg->getUseEarlySkipDetection())
            {
#if VCEG_AZ06_IC
              rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
              xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);//by Competition for inter_2Nx2N
            }
            // SKIP
#if VCEG_AZ06_IC
            if (!bICFlag)
            {
#endif
#if SIF_AFFINE_HARMONIZATION
              if (!modeSIF)
              {
#endif
#if COM16_C1016_AFFINE
#if JVET_D0077_FAST_EXT
                if (rpcTempCU->getSlice()->getSPS()->getUseAffine() && !bPrevSameBlockIsIntra
#if JVET_D0077_SAVE_LOAD_ENC_INFO
                  && !bFastSkipAffine
#endif
                  )
#else
                if (rpcTempCU->getSlice()->getSPS()->getUseAffine())
#endif
                {
                  xCheckRDCostAffineMerge2Nx2N(rpcBestCU, rpcTempCU);
                  rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                }
#endif
#if SIF_AFFINE_HARMONIZATION
              }
#endif

              xCheckRDCostMerge2Nx2N(rpcBestCU, rpcTempCU DEBUG_STRING_PASS_INTO(sDebug), &earlyDetectionSkipMode);//by Merge for inter_2Nx2N
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if VCEG_AZ06_IC
            }
#endif

#if VCEG_AZ07_FRUC_MERGE
#if JVET_D0077_FAST_EXT
            if (rpcTempCU->getSlice()->getSPS()->getUseFRUCMgrMode() && !bPrevSameBlockIsIntra
#if JVET_D0077_SAVE_LOAD_ENC_INFO
              && !bFastSkipFruc
#endif
              )
#else
            if (rpcTempCU->getSlice()->getSPS()->getUseFRUCMgrMode())
#endif
            {
#if VCEG_AZ06_IC
              rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
              xCheckRDCostMerge2Nx2NFRUC(rpcBestCU, rpcTempCU, &earlyDetectionSkipMode);
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
            }
#endif

#if JVET_C0024_QTBT
#if JVET_D0077_FAST_EXT
            if (!m_pcEncCfg->getUseEarlySkipDetection() && !bPrevSameBlockIsSkip && !bPrevSameBlockIsIntra
#if JVET_D0077_SAVE_LOAD_ENC_INFO
              && !bFastSkipInter
#if VCEG_AZ06_IC
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
              && !(saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
              && !(saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiWidthIdx, uiHeightIdx))
#endif
#endif
#endif
#else
            if (!m_pcEncCfg->getUseEarlySkipDetection() && !rpcBestCU->getPic()->getSkiped(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight)
              && !rpcBestCU->getPic()->getIntra(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight)
#endif
              )
#else
            if (!m_pcEncCfg->getUseEarlySkipDetection())
#endif
            {
              // 2Nx2N, NxN
#if VCEG_AZ06_IC
              rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
              xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if !JVET_C0024_QTBT
              if (m_pcEncCfg->getUseCbfFastMode())
              {
                doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
              }
#endif
            }
#if VCEG_AZ06_IC
          }
#endif
#if VCEG_AZ07_IMV && JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_FAST_EXT
          if (m_pcEncCfg->getIMV() && !rpcBestCU->getSlice()->isIntra() && !bPrevSameBlockIsIntra && !bPrevSameBlockIsSkip
#if JVET_D0077_SAVE_LOAD_ENC_INFO
            && !bFastSkipIMV
#endif
            )
#else
          if (m_pcEncCfg->getIMV() && !rpcBestCU->getSlice()->isIntra())
#endif
          {
#if JVET_E0076_MULTI_PEL_MVD
            for (UChar iMv = 1; iMv < 3; iMv++)
            {
              if (iMv == 2 && m_dBestMvDPelCost[0] * 1.06 < m_dBestMvDPelCost[1])
              {
                break;
              }
#if VCEG_AZ06_IC
              for (UInt uiICId = 0; uiICId < (bICEnabled ? 2 : 1); uiICId++)
              {
                Bool bICFlag = uiICId ? true : false;
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
                if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
                if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiWidthIdx, uiHeightIdx))
#endif
                {
                  continue;
                }
#endif
                rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N, false, iMv);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if VCEG_AZ06_IC
              }
#endif
            }
#else
            // always check SIZE_2Nx2N
#if VCEG_AZ06_IC
            for (UInt uiICId = 0; uiICId < (bICEnabled ? 2 : 1); uiICId++)
            {
              Bool bICFlag = uiICId ? true : false;
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
              if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
              if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiWidthIdx, uiHeightIdx))
#endif
              {
                continue;
              }
#endif
              rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
              xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N, false, true);
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if VCEG_AZ06_IC
            }
#endif
#endif
          }
#endif
        }

        if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
        {
          iQP = iMinQP;
        }
#if RRC_SELECTABLE_IF
      }
#endif
    }
#if VCEG_AZ06_IC
#if !JVET_C0024_QTBT
    Bool bTestICNon2Nx2N = bICEnabled && rpcBestCU->getICFlag(0);
#if VCEG_AZ06_IC_SPEEDUP
    bTestICNon2Nx2N = false;
#endif
#endif
#endif
    if (!earlyDetectionSkipMode)
    {
      for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
      {
        const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP); // If lossless, then iQP is irrelevant for subsequent modules.
#if RRC_SELECTABLE_IF
        int modeSIFnum = 1;
        if (rpcTempCU->isCodeBIF(uiWidth, uiDepth
#if JVET_C0024_QTBT && SIF_QTBT_FIX
          + uiBTDepth // check total depth
#endif
          )
#if SIF_SPLIT_ONLY
          && 0
#endif
#if SIF_CU_LEVEL
          || 1
#endif
          )
        {
          modeSIFnum = 2;
        }
        for (int modeSIF = 0; modeSIF < modeSIFnum; modeSIF++)
        {
          if (modeSIFnum == 2)
          {
            rpcTempCU->setSIFToCheck(modeSIF);
            rpcBestCU->setSIFToCheck(modeSIF);
          }
#endif

          if (bIsLosslessMode)
          {
            iQP = lowestQP;
          }

          rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

#if !JVET_C0024_QTBT
          // do inter modes, NxN, 2NxN, and Nx2N
          if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
          {
#if VCEG_AZ06_IC
            for (UInt uiICId = 0; uiICId < (bTestICNon2Nx2N ? 2 : 1); uiICId++)
            {
              Bool bICFlag = uiICId ? true : false;
#endif
              // 2Nx2N, NxN
#if !COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL || COM16_C806_DISABLE_4X4_PU
              if (!((rpcBestCU->getWidth(0) == 8) && (rpcBestCU->getHeight(0) == 8)))
#else
              if (rpcBestCU->getSlice()->getSPS()->getAtmvpEnableFlag() || !((rpcBestCU->getWidth(0) == 8) && (rpcBestCU->getHeight(0) == 8)))
#endif
              {
                if (uiDepth == sps.getLog2DiffMaxMinCodingBlockSize() && doNotBlockPu)
                {
#if VCEG_AZ06_IC
                  rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                  xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug));
                  rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                }
              }

              if (doNotBlockPu)
              {
#if VCEG_AZ06_IC
                rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_Nx2N DEBUG_STRING_PASS_INTO(sDebug));
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_Nx2N)
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                }
              }
              if (doNotBlockPu)
              {
#if VCEG_AZ06_IC
                rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxN DEBUG_STRING_PASS_INTO(sDebug));
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxN)
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                }
              }
#if VCEG_AZ06_IC
            }
            bTestICNon2Nx2N &= rpcBestCU->getICFlag(0);
#endif
            //! Try AMP (SIZE_2NxnU, SIZE_2NxnD, SIZE_nLx2N, SIZE_nRx2N)
            if (sps.getUseAMP() && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize())
            {
#if VCEG_AZ06_IC
              for (UInt uiICId = 0; uiICId < (bTestICNon2Nx2N ? 2 : 1); uiICId++)
              {
                Bool bICFlag = uiICId ? true : false;
#endif
#if AMP_ENC_SPEEDUP
                Bool bTestAMP_Hor = false, bTestAMP_Ver = false;

#if AMP_MRG
                Bool bTestMergeAMP_Hor = false, bTestMergeAMP_Ver = false;

                deriveTestModeAMP(rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver, bTestMergeAMP_Hor, bTestMergeAMP_Ver
                  );
#else
                deriveTestModeAMP(rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver
                  );
#endif

                //! Do horizontal AMP
                if (bTestAMP_Hor)
                {
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                }
#if AMP_MRG
                else if (bTestMergeAMP_Hor)
                {
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug), true);
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug), true);
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                }
#endif

                //! Do horizontal AMP
                if (bTestAMP_Ver)
                {
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                  }
                }
#if AMP_MRG
                else if (bTestMergeAMP_Ver)
                {
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug), true);
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                    if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N)
                    {
                      doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
                    }
                  }
                  if (doNotBlockPu)
                  {
#if VCEG_AZ06_IC
                    rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                    xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug), true);
                    rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                  }
                }
#endif

#else
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

#endif
#if VCEG_AZ06_IC
              }
#endif
            }
          }
#endif //#if !JVET_C0024_QTBT

#if VCEG_AZ07_IMV && !JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_FAST_EXT
          if (m_pcEncCfg->getIMV() && !rpcBestCU->getSlice()->isIntra() && !bPrevSameBlockIsIntra && !bPrevSameBlockIsSkip
#if JVET_D0077_SAVE_LOAD_ENC_INFO
            && !bFastSkipIMV
#endif
#if SIF_HARMO_IMV
            && !modeSIF
#endif
            )
#else
          if (m_pcEncCfg->getIMV() && !rpcBestCU->getSlice()->isIntra()
#if SIF_HARMO_IMV
            && !modeSIF
#endif
            )
#endif
          {
            // always check SIZE_2Nx2N
#if JVET_E0076_MULTI_PEL_MVD
            for (UChar iMv = 1; iMv < 3; iMv++)
            {
              if (iMv == 2 && m_dBestMvDPelCost[0] * 1.06 < m_dBestMvDPelCost[1])
              {
                break;
              }
#if VCEG_AZ06_IC
              for (UInt uiICId = 0; uiICId < (bICEnabled ? 2 : 1); uiICId++)
              {
                Bool bICFlag = uiICId ? true : false;
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
                if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
                if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiWidthIdx, uiHeightIdx))
#endif
                {
                  continue;
                }
#endif
                rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
                xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N, false, iMv);
                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if VCEG_AZ06_IC
              }
            }
#endif
#else
#if VCEG_AZ06_IC
            for (UInt uiICId = 0; uiICId < (bICEnabled ? 2 : 1); uiICId++)
            {
              Bool bICFlag = uiICId ? true : false;
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
              if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
              if (saveLoadTag == LOAD_ENC_INFO && bICFlag != m_pcPredSearch->getSaveLoadICFlag(uiWidthIdx, uiHeightIdx))
#endif
              {
                continue;
              }
#endif
              rpcTempCU->setICFlagSubParts(bICFlag, 0, uiDepth);
#endif
              xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N, false, true);
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#if VCEG_AZ06_IC
            }
#endif
#endif
#if !JVET_C0024_QTBT
            // check other candidates
            const Char cMinCandNum = 1, cMaxCandNum = m_pcEncCfg->getIMVMaxCand();
            Char cMaxNeighboriMVCandNum = max(rpcTempCU->getMaxNeighboriMVCandNum(0), cMinCandNum);
            std::set<SModeCand, cmpModeCand>::iterator pos;
            Char n, nBestCand = 0;
            for (n = 1, pos = m_setInterCand.begin(); n <= cMaxNeighboriMVCandNum && pos != m_setInterCand.end(); n++, pos++)
            {
              Double dCurrentRD = rpcBestCU->getTotalCost();

              if (pos->eInterPartSize == SIZE_2Nx2N)
              {
                continue;
              }

              xCheckRDCostInter(rpcBestCU, rpcTempCU, pos->eInterPartSize, false, true, pos->bUseMrg ? NULL : pos->pcCUMode);

              if (rpcTempCU->getTotalCost() == MAX_DOUBLE)
              {
                cMaxNeighboriMVCandNum++;
              }

              if (rpcBestCU->getTotalCost() < dCurrentRD)
              {
                nBestCand = n;
                rpcBestCU->setiMVCandNumSubParts(nBestCand, 0, uiDepth);
                cMaxNeighboriMVCandNum += 1;
                cMaxNeighboriMVCandNum = min(cMaxCandNum, cMaxNeighboriMVCandNum);
              }
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
            }
            for (; n <= cMaxNeighboriMVCandNum + 1 && pos != m_setInterCand.end(); n++, pos++)
            {
              Double dCurrentRD = rpcBestCU->getTotalCost();

              if (pos->eInterPartSize == SIZE_2Nx2N)
              {
                continue;
              }

              xCheckRDCostInter(rpcBestCU, rpcTempCU, pos->eInterPartSize, false, true, pos->bUseMrg ? NULL : pos->pcCUMode);

              if (rpcTempCU->getTotalCost() == MAX_DOUBLE)
              {
                cMaxNeighboriMVCandNum++;
              }

              if (rpcBestCU->getTotalCost() < dCurrentRD)
              {
                nBestCand = n;
                rpcBestCU->setiMVCandNumSubParts(nBestCand, 0, uiDepth);
              }
              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
            }
#endif //#if !JVET_C0024_QTBT
          }
#endif
#if RRC_SELECTABLE_IF
        }// for modeSIF
#endif

        // do normal intra modes
        // speedup for inter frames
        Double intraCost = 0.0;

#if COM16_C806_EMT
        if (rpcBestCU->getSlice()->getSliceType() != I_SLICE && m_pcEncCfg->getUseFastInterEMT())
        {
          dBestInterCost = rpcBestCU->getTotalCost();
        }
#endif

        if ((rpcBestCU->getSlice()->getSliceType() == I_SLICE) ||
#if JVET_C0024_QTBT
          ((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && (!rpcBestCU->getPic()->getInter(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight)) && (
#else
            ((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && (
#endif
              (rpcBestCU->getCbf(0, COMPONENT_Y) != 0) ||
              ((rpcBestCU->getCbf(0, COMPONENT_Cb) != 0) && (numberValidComponents > COMPONENT_Cb)) ||
              ((rpcBestCU->getCbf(0, COMPONENT_Cr) != 0) && (numberValidComponents > COMPONENT_Cr))  // avoid very complex intra if it is unlikely
              )))
        {
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
          Int bNonZeroCoeff = 0;
#endif
#if VCEG_AZ05_INTRA_MPI
          Char iMPIidx = 0;  Char iNumberOfPassesMPI = 2;
          if (rpcTempCU->getSlice()->getSliceType() != I_SLICE)  iNumberOfPassesMPI = 2;
          if (!rpcTempCU->getSlice()->getSPS()->getUseMPI()) iNumberOfPassesMPI = 1;
#if JVET_C0024_QTBT
          if (isChroma(rpcBestCU->getTextType()))
          {
            iNumberOfPassesMPI = 1;
          }
#endif
#endif
#if JVET_C0024_PBINTRA_FAST
          rpcTempCU->getInterHAD() = MAX_UINT;
          if (rpcBestCU->getPredictionMode(0) == MODE_INTER && !rpcBestCU->getSlice()->isIntra())
          {
#if JVET_C0024_PBINTRA_FAST_FIX          
            // redundant MC process to make sure that m_pppcPredYuvBest has the correct moiton compensation prediction data
#if JVET_E0052_DMVR
            m_pcPredSearch->motionCompensation(rpcBestCU, m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx], false);
#if COM16_C806_OBMC
            m_pcPredSearch->subBlockOBMC(rpcBestCU, 0, m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv1[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv2[uiWidthIdx][uiHeightIdx], false);
#endif
#else
            m_pcPredSearch->motionCompensation(rpcBestCU, m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx]);
#if COM16_C806_OBMC
            m_pcPredSearch->subBlockOBMC(rpcBestCU, 0, m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv1[uiWidthIdx][uiHeightIdx], m_pppcTmpYuv2[uiWidthIdx][uiHeightIdx]);
#endif
#endif
            DistParam distParam;
            const Bool bUseHadamard = rpcTempCU->getCUTransquantBypass(0) == 0;
            m_pcRdCost->setDistParam(distParam, rpcTempCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->getAddr(COMPONENT_Y)
              , m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->getStride(COMPONENT_Y)
              , m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx]->getAddr(COMPONENT_Y), m_pppcPredYuvBest[uiWidthIdx][uiHeightIdx]->getStride(COMPONENT_Y)
              , rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), bUseHadamard);
#else
            m_pcPredSearch->motionCompensation(rpcTempCU, m_pppcPredYuvTemp[uiWidthIdx][uiHeightIdx]);
            // use hadamard transform here
            DistParam distParam;
            const Bool bUseHadamard = rpcTempCU->getCUTransquantBypass(0) == 0;
            m_pcRdCost->setDistParam(distParam, rpcTempCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->getAddr(COMPONENT_Y)
              , m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]->getStride(COMPONENT_Y)
              , m_pppcPredYuvTemp[uiWidthIdx][uiHeightIdx]->getAddr(COMPONENT_Y), m_pppcPredYuvTemp[uiWidthIdx][uiHeightIdx]->getStride(COMPONENT_Y)
              , rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), bUseHadamard);
#endif
            distParam.bApplyWeight = false;

            UInt uiSad = distParam.DistFunc(&distParam);

            rpcTempCU->getInterHAD() = uiSad;
          }
#endif


#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
          const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
          const UInt uiDIMDMaxBtDepht = (rpcBestCU->getBoundaryCtuLocation() && rpcBestCU->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
          UInt uiDIMDIntraMode = DC_IDX;
#if COM16_C1046_PDPC_INTRA
          UInt uiDIMDIntraModePDPC = DC_IDX;
#endif
          if (isLuma(rpcBestCU->getTextType()))
          {
            if (!m_pcEncCfg->getUseFastLCTU() || rpcTempCU->getWidth(0)*rpcTempCU->getHeight(0) <= 4096)
            {
#if BT_FORCE_DIMD_HARM
              if (rpcTempCU->getBTDepth(0) <= uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
              if (rpcTempCU->getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
              {
                uiDIMDIntraMode = m_pcPredSearch->deriveNeighborIntraDirs(rpcTempCU, 0);
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                if (rpcTempCU->getSlice()->getSPS()->getUsePDPC())
                {
                  rpcTempCU->setPDPCIdxSubParts(1, 0, uiDepth);
                  uiDIMDIntraModePDPC = m_pcPredSearch->deriveNeighborIntraDirs(rpcTempCU, 0);
                }
#endif
#if BT_FORCE_DIMD_HARM
                if (rpcTempCU->getBTDepth(0) == uiDIMDMaxBtDepht && rpcTempCU->getSlice()->isIntra())
#else
                if (rpcTempCU->getBTDepth(0) == DIMD_MAX_BT_DEPTH_SIG_AI && rpcTempCU->getSlice()->isIntra())
#endif
                {
                  m_pcPredSearch->setParentDIMDIntraMode(uiDIMDIntraMode);
                }
              }
              else
              {
                assert(rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0));
                if (rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0))
                {
                  uiDIMDIntraMode = m_pcPredSearch->getParentDIMDIntraMode();
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                  if (rpcTempCU->getSlice()->getSPS()->getUsePDPC())
                  {
                    uiDIMDIntraModePDPC = m_pcPredSearch->getParentDIMDIntraMode();
                  }
#endif
                }
                else
                {
                  uiDIMDIntraMode = m_pcPredSearch->deriveNeighborIntraDirs(rpcTempCU, 0);
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                  if (rpcTempCU->getSlice()->getSPS()->getUsePDPC())
                  {
                    rpcTempCU->setPDPCIdxSubParts(1, 0, uiDepth);
                    uiDIMDIntraModePDPC = m_pcPredSearch->deriveNeighborIntraDirs(rpcTempCU, 0);
                  }
#endif
                }
              }
            }
          }

          Bool bTestNonDIMDIntra = true;
          if (isLuma(rpcTempCU->getTextType()))
          {
#if BT_FORCE_DIMD_HARM
            if (rpcTempCU->getBTDepth(0) <= uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
            if (rpcTempCU->getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
            {
              bTestNonDIMDIntra = true;
            }
            else
            {
              assert(rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0));
              bTestNonDIMDIntra = !rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0);
            }
          }

          Bool bTestDIMDIntra = true;
          if (isChroma(rpcBestCU->getTextType()))
          {
            bTestDIMDIntra = false;
            assert(bTestNonDIMDIntra);
          }
          if (isLuma(rpcTempCU->getTextType()))
          {
#if BT_FORCE_DIMD_HARM
            if (rpcTempCU->getBTDepth(0) <= uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
            if (rpcTempCU->getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
            {
              bTestDIMDIntra = true;
            }
            else
            {
              assert(rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0) == rpcBestCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0));
              bTestDIMDIntra = rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0);
            }
          }

          Bool btestDIMD[2];
          btestDIMD[0] = bTestNonDIMDIntra;
          btestDIMD[1] = bTestDIMDIntra;

          for (int dimdIdx = 0; dimdIdx < 2; dimdIdx++)
          {
#endif
#if INTRA_LOOPS_ORDER_FIX
#if RRC_INTRA_SHARPENING
            int iNumberOfPassesIntraSharp = rpcTempCU->getSlice()->getSPS()->getUseIntraSharp() ? 2 : 1;
            if (isChroma(rpcBestCU->getTextType()))
            {
              iNumberOfPassesIntraSharp = 1;
            }

#if DIMD_INTRA_SHARP_HARM
            if ((dimdIdx == 1) && btestDIMD[dimdIdx])
            { // TestDIMDIntra settings

              iNumberOfPassesIntraSharp = 1;
            }
#endif
            for (Int iSharpIdx = 0; iSharpIdx < iNumberOfPassesIntraSharp; iSharpIdx++)
            {
#if JVET_D0077_SAVE_LOAD_ENC_INFO && INTRA_SHARP_D0077
              if (CHANNEL_TYPE_LUMA == eChannelType && saveLoadTag == LOAD_ENC_INFO && iSharpIdx != m_pcPredSearch->getSaveLoadSharpIdx(uiWidthIdx, uiHeightIdx))
              {
                continue;
              }
#endif
#endif

#if RRC_DWDIP
              const Bool bUseDWDIP = rpcBestCU->getSlice()->getSPS()->getUseDWDIP();
#endif

#if RRC_DWDIP_MLR
              Int iNumberOfPassesDWDIP_GP = 2;

              if (isChroma(rpcBestCU->getSlice()->getTextType()) || !bUseDWDIP
#if RRC_INTRA_SHARPENING
                || iSharpIdx != 0
#endif
                )
                iNumberOfPassesDWDIP_GP = 1;

#if DIMD_DWDIP_HARM 
              if ((dimdIdx == 1) && btestDIMD[dimdIdx])
              { // TestDIMDIntra settings
                iNumberOfPassesDWDIP_GP = 1;
              }
#endif

              for (Int iDWDIP_GP = 0; iDWDIP_GP < iNumberOfPassesDWDIP_GP; ++iDWDIP_GP)
              {
                rpcTempCU->setDWDIP_GPIdxSubParts(iDWDIP_GP, 0, uiDepth);
#endif

#if RRC_DWDIP
                int iDWDIPPasses = (!bUseDWDIP ||
                  !isEnabledForDWDIP(rpcBestCU->getTextType(), rpcTempCU->getWidth(0), rpcTempCU->getHeight(0))
                  //#if COM16_C1046_PDPC_INTRA 
                  //           || iPDPCidx != 0
                  //#endif
#if RRC_INTRA_SHARPENING
                  || iSharpIdx != 0
#endif
                  ) ?
                  1 :
                  2
                  ;

#if DIMD_DWDIP_HARM 
                if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                { // TestDIMDIntra settings           
                  iDWDIPPasses = 1;
                } // end of TestDIMDIntra settings
#endif
                for (Int iDWDIPIdx = 0; iDWDIPIdx < iDWDIPPasses; iDWDIPIdx++)
                {
                  rpcTempCU->setDWDIPIdxSubParts(iDWDIPIdx, 0, uiDepth);
#endif
#endif //INTRA_LOOPS_ORDER_FIX

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                  Char iPDPCidx = 0;  Char iNumberOfPassesPDPC = 2;
                  if (rpcTempCU->getSlice()->getSliceType() != I_SLICE)  iNumberOfPassesPDPC = 2;
                  if (!rpcTempCU->getSlice()->getSPS()->getUsePDPC()) iNumberOfPassesPDPC = 1;
#if JVET_C0024_QTBT
                  if (isChroma(rpcBestCU->getTextType()))
                  {
                    iNumberOfPassesPDPC = 1;
                  }
#endif
#endif


#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                  Char iROTidx = 0; Char iNumberOfPassesROT = 4;
#if JVET_C0024_QTBT
                  if (isChroma(rpcBestCU->getTextType()) && !(uiWidth >= 8 && uiHeight >= 8))
                  {
                    iNumberOfPassesROT = 1;
                  }
#endif
#if COM16_C1044_NSST
                  if (!rpcTempCU->getSlice()->getSPS()->getUseNSST()) iNumberOfPassesROT = 1;
#else
                  if (!rpcTempCU->getSlice()->getSPS()->getUseROT()) iNumberOfPassesROT = 1;
#endif

                  for (iROTidx = 0; iROTidx < iNumberOfPassesROT; iROTidx++)
                  {
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                    if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                    { // TestNonDIMDIntra settings
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
                      if (saveLoadTag == LOAD_ENC_INFO && iROTidx && iROTidx != m_pcPredSearch->getSaveLoadRotIdx(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
                      if (saveLoadTag == LOAD_ENC_INFO && iROTidx && iROTidx != m_pcPredSearch->getSaveLoadRotIdx(uiWidthIdx, uiHeightIdx))
#endif
                      {
                        continue;
                      }
#endif
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                    } // end of TestNonDIMDIntra settings
#endif
#endif

#if !INTRA_LOOPS_ORDER_FIX
#if RRC_INTRA_SHARPENING
                    int iNumberOfPassesIntraSharp = rpcTempCU->getSlice()->getSPS()->getUseIntraSharp() ? 2 : 1;
                    if (isChroma(rpcBestCU->getTextType()))
                    {
                      iNumberOfPassesIntraSharp = 1;
                    }

#if DIMD_INTRA_SHARP_HARM
                    if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                    { // TestDIMDIntra settings

                      iNumberOfPassesIntraSharp = 1;
                    }
#endif
                    for (Int iSharpIdx = 0; iSharpIdx < iNumberOfPassesIntraSharp; iSharpIdx++)
                    {
#if JVET_D0077_SAVE_LOAD_ENC_INFO && INTRA_SHARP_D0077
                      if (CHANNEL_TYPE_LUMA == eChannelType && saveLoadTag == LOAD_ENC_INFO && iSharpIdx != m_pcPredSearch->getSaveLoadSharpIdx(uiWidthIdx, uiHeightIdx))
                      {
                        continue;
                      }
#endif
#endif
#endif //#if !INTRA_LOOPS_ORDER_FIX

#if VCEG_AZ05_INTRA_MPI
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                      if (iROTidx) iNumberOfPassesMPI = 1;
#endif
                      for (iMPIidx = 0; iMPIidx < iNumberOfPassesMPI; iMPIidx++)
                      {
#endif

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                        if (iROTidx) iNumberOfPassesPDPC = 1;
#endif
                        for (iPDPCidx = 0; iPDPCidx < iNumberOfPassesPDPC; iPDPCidx++)
                        {
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                          if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                          { // TestNonDIMDIntra settings
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
                            if (CHANNEL_TYPE_LUMA == eChannelType && saveLoadTag == LOAD_ENC_INFO && iPDPCidx != m_pcPredSearch->getSaveLoadPdpcIdx(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
                            if (CHANNEL_TYPE_LUMA == eChannelType && saveLoadTag == LOAD_ENC_INFO && iPDPCidx != m_pcPredSearch->getSaveLoadPdpcIdx(uiWidthIdx, uiHeightIdx))
#endif
                            {
                              continue;
                            }
#endif
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                          } // end of TestNonDIMDIntra settings
#endif
#endif

#if !INTRA_LOOPS_ORDER_FIX
#if RRC_DWDIP
                          const Bool bUseDWDIP = rpcBestCU->getSlice()->getSPS()->getUseDWDIP();
#endif

#if RRC_DWDIP_MLR
                          Int iNumberOfPassesDWDIP_GP = 2;

                          if (isChroma(rpcBestCU->getSlice()->getTextType()) || !bUseDWDIP
#if RRC_INTRA_SHARPENING
                            || iSharpIdx != 0
#endif
                            )
                            iNumberOfPassesDWDIP_GP = 1;

#if DIMD_DWDIP_HARM 
                          if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                          { // TestDIMDIntra settings

                            iNumberOfPassesDWDIP_GP = 1;
                          }
#endif

                          for (Int iDWDIP_GP = 0; iDWDIP_GP < iNumberOfPassesDWDIP_GP; ++iDWDIP_GP)
                          {
                            rpcTempCU->setDWDIP_GPIdxSubParts(iDWDIP_GP, 0, uiDepth);
#endif

#if RRC_DWDIP
                            int iDWDIPPasses = (!bUseDWDIP ||
                              !isEnabledForDWDIP(rpcBestCU->getTextType(), rpcTempCU->getWidth(0), rpcTempCU->getHeight(0))
#if COM16_C1046_PDPC_INTRA 
                              || iPDPCidx != 0
#endif
#if RRC_INTRA_SHARPENING
                              || iSharpIdx != 0
#endif
                              ) ?
                              1 :
                              2
                              ;

#if DIMD_DWDIP_HARM 
                            if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                            { // TestDIMDIntra settings           
                              iDWDIPPasses = 1;
                            } // end of TestDIMDIntra settings
#endif
                            for (Int iDWDIPIdx = 0; iDWDIPIdx < iDWDIPPasses; iDWDIPIdx++)
                            {
                              rpcTempCU->setDWDIPIdxSubParts(iDWDIPIdx, 0, uiDepth);
#endif
#endif //!INTRA_LOOPS_ORDER_FIX

#if COM16_C806_EMT
#if JVET_C0024_QTBT
                              UChar ucEmtUsage = (uiWidth > EMT_INTRA_MAX_CU || uiHeight > EMT_INTRA_MAX_CU || (rpcTempCU->getSlice()->getSPS()->getUseIntraEMT() == 0) || isChroma(rpcBestCU->getTextType())) ? 1 : 2;
#else
                              UChar ucEmtUsage = ((rpcTempCU->getWidth(0) > EMT_INTRA_MAX_CU) || (rpcTempCU->getSlice()->getSPS()->getUseIntraEMT() == 0)) ? 1 : 2;
#endif
                              for (UChar ucCuFlag = 0; ucCuFlag < ucEmtUsage; ucCuFlag++)
                              {
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                                { // TestNonDIMDIntra settings             
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
                                  if (saveLoadTag == LOAD_ENC_INFO && /*ucCuFlag && */ucCuFlag != m_pcPredSearch->getSaveLoadEmtFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx))
#else
                                  if (saveLoadTag == LOAD_ENC_INFO && /*ucCuFlag && */ucCuFlag != m_pcPredSearch->getSaveLoadEmtFlag(uiWidthIdx, uiHeightIdx))
#endif
                                  {
                                    continue;
                                  }
#endif
                                  if (ucCuFlag && bEarlySkipIntra && m_pcEncCfg->getUseFastInterEMT())
                                  {
                                    continue;
                                  }

#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                } // end of TestNonDIMDIntra settings
#endif

                                rpcTempCU->setEmtCuFlagSubParts(ucCuFlag, 0, uiDepth);
#if VCEG_AZ05_INTRA_MPI
                                rpcTempCU->setMPIIdxSubParts(iMPIidx, 0, uiDepth);
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                                rpcTempCU->setPDPCIdxSubParts(iPDPCidx, 0, uiDepth);
#endif
#if MLIP_EXTEND_REF_LINE
                                rpcTempCU->setExternRefSubParts(rpcTempCU->getTextType(), 0, 0, uiDepth);
#endif
#if BARC_CORE
                                rpcTempCU->setBARCFlag(rpcTempCU->getTextType(), isBARCFlag );
#endif
#if RRC_DWDIP
                                rpcTempCU->setDWDIPIdxSubParts(iDWDIPIdx, 0, uiDepth);
#endif

#if RRC_DWDIP_MLR
                                rpcTempCU->setDWDIP_GPIdxSubParts(iDWDIP_GP, 0, uiDepth);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
                                rpcTempCU->setROTIdxSubParts(rpcTempCU->getTextType(), iROTidx, 0, uiDepth);
                                if (!rpcTempCU->getSlice()->isIntra())
                                {
                                  rpcTempCU->setROTIdxSubParts(CHANNEL_TYPE_CHROMA, iROTidx, 0, uiDepth);
                                }
#else
                                rpcTempCU->setROTIdxSubParts(iROTidx, 0, uiDepth);
#endif
#endif

#if RRC_INTRA_SHARPENING         
                                rpcTempCU->setIntraSharpFlagSubParts(iSharpIdx, 0, uiDepth);
#endif
#if INTRA_SHARP_PDPC_HARM
                                if (iSharpIdx && iPDPCidx)
                                  continue;
#endif
#if INTRA_LOOPS_ORDER_FIX && RRC_DWDIP
#if JVET_G0104_PLANAR_PDPC
                                if (iDWDIPIdx)
#else
																if (iDWDIPIdx && iPDPCidx)
#endif
                                  continue;
#endif
#if INTRA_SHARP_MAX_ASPECT_RATIO
                                int AspectRatio = std::max(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0)) / std::min(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0));
                                if (iSharpIdx && AspectRatio > INTRA_SHARP_MAX_ASPECT_RATIO)
                                  continue;
#endif
#if INTRA_SHARP_MAX_SIZE
                                if (iSharpIdx && std::max(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0)) > INTRA_SHARP_MAX_SIZE)
                                  continue;
#endif

#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                                { // TestNonDIMDIntra settings             
#endif
#if JVET_C0024_PBINTRA_FAST
                                  if (rpcBestCU->getPredictionMode(0) == MODE_INTER && !rpcBestCU->getSlice()->isIntra())
                                  {
                                    if (rpcTempCU->getInterHAD() == 0)
                                    {
                                      continue;  //has calculated the best intra HAD much larger than that of inter; so skip intra mode RDO
                                    }
                                  }
#endif
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                } // end of TestNonDIMDIntra settings
#endif

             /* DIMD processing start */
#if DIMD_INTRA_PRED
                                if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                                { // TestNonDIMDIntra settings
#if DIMD_NUM_INTRA_DIR_INC
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                                  UInt uiDIMDDir = iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode;
#else
                                  UInt uiDIMDDir = uiDIMDIntraMode;
#endif
#if VCEG_AZ07_INTRA_65ANG_MODES
                                  m_pcPredSearch->setDIMDIntraMode(MAP131TO67(uiDIMDDir));
#else
                                  m_pcPredSearch->setDIMDIntraMode(MAP131TO35(uiDIMDDir));
#endif
#else
                                  m_pcPredSearch->setDIMDIntraMode(iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode);
#endif

                                }
                                else if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                                { // TestDIMDIntra settings
                                  UInt uiIntraMode = DC_IDX;
                                  if (isLuma(rpcTempCU->getTextType()))  // luma of I-slice, and luma and chroma of P/B-slice
                                  {
#if JVET_G0104_PLANAR_PDPC
																		uiIntraMode = uiDIMDIntraMode;
#else
                                    uiIntraMode = iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode;
#endif
                                  }
                                  else                                  // chroma 
                                  {
                                    assert(0);
                                  }

                                  if (isLuma(rpcTempCU->getTextType()))
                                  {
                                    UChar ucSecIntraMode = rpcTempCU->generateSecIntraDir(0, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0));
                                    rpcTempCU->setSecIntraDirSubParts(ucSecIntraMode, 0, uiDepth);
                                  }

                                  if (rpcTempCU->getSlice()->isIntra())
                                  {
#if DIMD_NUM_INTRA_DIR_INC
                                    rpcTempCU->setExtIntraDirSubParts(rpcTempCU->getTextType(), uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO67(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#else
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO35(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#endif
#else
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? uiIntraMode : DM_CHROMA_IDX, 0, uiDepth);
#endif
                                  }
#if	DIMD_DISABLE_CHROMA_P_B_SLICE_JB
                                  else
                                  {
#if DIMD_NUM_INTRA_DIR_INC
                                    rpcTempCU->setExtIntraDirSubParts(rpcTempCU->getTextType(), uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES								
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO67(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#else
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO35(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#endif
#else
                                    rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? uiIntraMode : DM_CHROMA_IDX, 0, uiDepth);
#endif

                                  }
#else
                                  else
                                  {
#if DIMD_NUM_INTRA_DIR_INC
                                    rpcTempCU->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiIntraMode, 0, uiDepth);
                                    rpcTempCU->setExtIntraDirSubParts(CHANNEL_TYPE_CHROMA, uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
                                    rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiIntraMode), 0, uiDepth);
#else
                                    rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiIntraMode), 0, uiDepth);
#endif
#else
                                    rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiIntraMode, 0, uiDepth);
#endif
#if	DIMD_DISABLE_CHROMA_P_B_SLICE_JB

#else
                                    rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_CHROMA, DM_CHROMA_IDX, 0, uiDepth);
#endif
                                  }
#endif
                                }

#endif
                                /* DIMD processing end */
#if DIMD_INTRA_PRED
                                if (btestDIMD[dimdIdx])
                                {
#endif
                                  xCheckRDCostIntra(rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug)
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                                    , bNonZeroCoeff
#endif
#if DIMD_INTRA_PRED
                                    , dimdIdx
#endif
// #if DIMD_MLIP_HARM
//                                       , ((dimdIdx == 1) && btestDIMD[dimdIdx])
// #endif
                                    );
#if DIMD_INTRA_PRED
                                }
#endif

#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                                { // TestNonDIMDIntra settings             
#endif
                                  if (!ucCuFlag && !rpcBestCU->isIntra(0) && m_pcEncCfg->getUseFastInterEMT()
#if JVET_D0077_SAVE_LOAD_ENC_INFO
                                    && dBestInterCost < MAX_DOUBLE
#endif
                                    )
                                  {
                                    static const Double thEmtInterFastSkipIntra = 1.4; // Skip checking Intra if "2Nx2N using DCT2" is worse than best Inter mode
                                    if (rpcTempCU->getTotalCost() > thEmtInterFastSkipIntra*dBestInterCost)
                                    {
                                      bEarlySkipIntra = true;
                                    }
                                  }
                                  if (!ucCuFlag && m_pcEncCfg->getUseFastIntraEMT())
                                  {
#if JVET_C0024_QTBT
                                    //dIntra2Nx2NCost = rpcBestCU->isIntra(0) ? rpcBestCU->getTotalCost() : rpcTempCU->getTotalCost();
#else
                                    dIntra2Nx2NCost = (rpcBestCU->isIntra(0) && rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N) ? rpcBestCU->getTotalCost() : rpcTempCU->getTotalCost();
#endif
                                  }
#if HEM_2_0_MIMICS && DIMD_INTRA_PRED
                                } // end of TestNonDIMDIntra settings
#endif
                                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                              }
#else
#if VCEG_AZ05_INTRA_MPI
                              rpcTempCU->setMPIIdxSubParts(iMPIidx, 0, uiDepth);
#endif
#if COM16_C1046_PDPC_INTRA
                              rpcTempCU->setPDPCIdxSubParts(iPDPCidx, 0, uiDepth);
#endif
#if MLIP_EXTEND_REF_LINE
                              rpcTempCU->setExternRefSubParts(rpcTempCU->getTextType(), 0, 0, uiDepth);
#endif
#if RRC_DWDIP
                              rpcTempCU->setDWDIPIdxSubParts(iDWDIPIdx, 0, uiDepth);
#endif

#if RRC_INTRA_SHARPENING
                              rpcTempCU->setIntraSharpFlagSubParts(iSharpIdx, 0, uiDepth);
#endif


#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
                              rpcTempCU->setROTIdxSubParts(rpcTempCU->getTextType(), iROTidx, 0, uiDepth);
                              if (!rpcTempCU->getSlice()->isIntra())
                              {
                                rpcTempCU->setROTIdxSubParts(CHANNEL_TYPE_CHROMA, iROTidx, 0, uiDepth);
                              }
#else
                              rpcTempCU->setROTIdxSubParts(iROTidx, 0, uiDepth);
#endif
#endif

                              /* DIMD processing start */
#if DIMD_INTRA_PRED
                              if ((dimdIdx == 0) && btestDIMD[dimdIdx])
                              { // TestNonDIMDIntra settings
#if DIMD_NUM_INTRA_DIR_INC
#if COM16_C1046_PDPC_INTRA
                                UInt uiDIMDDir = iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode;
#else
                                UInt uiDIMDDir = uiDIMDIntraMode;
#endif
#if VCEG_AZ07_INTRA_65ANG_MODES
                                m_pcPredSearch->setDIMDIntraMode(MAP131TO67(uiDIMDDir));
#else
                                m_pcPredSearch->setDIMDIntraMode(MAP131TO35(uiDIMDDir));
#endif
#else
                                m_pcPredSearch->setDIMDIntraMode(iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode);
#endif

                              }
                              else if ((dimdIdx == 1) && btestDIMD[dimdIdx])
                              { // TestDIMDIntra settings
                                UInt uiIntraMode = DC_IDX;
                                if (isLuma(rpcTempCU->getTextType()))  // luma of I-slice, and luma and chroma of P/B-slice
                                {
                                  uiIntraMode = iPDPCidx ? uiDIMDIntraModePDPC : uiDIMDIntraMode;
                                }
                                else                                  // chroma 
                                {
                                  assert(0);
                                }

                                if (isLuma(rpcTempCU->getTextType()))
                                {
                                  UChar ucSecIntraMode = rpcTempCU->generateSecIntraDir(0, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0));
                                  rpcTempCU->setSecIntraDirSubParts(ucSecIntraMode, 0, uiDepth);
                                }

                                if (rpcTempCU->getSlice()->isIntra())
                                {
#if DIMD_NUM_INTRA_DIR_INC
                                  rpcTempCU->setExtIntraDirSubParts(rpcTempCU->getTextType(), uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO67(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#else
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO35(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#endif
#else
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? uiIntraMode : DM_CHROMA_IDX, 0, uiDepth);
#endif
                                }
#if	DIMD_DISABLE_CHROMA_P_B_SLICE_JB
                                else
                                {
#if DIMD_NUM_INTRA_DIR_INC
                                  rpcTempCU->setExtIntraDirSubParts(rpcTempCU->getTextType(), uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES								
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO67(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#else
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? MAP131TO35(uiIntraMode) : DM_CHROMA_IDX, 0, uiDepth);
#endif
#else
                                  rpcTempCU->setIntraDirSubParts(rpcTempCU->getTextType(), (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA) ? uiIntraMode : DM_CHROMA_IDX, 0, uiDepth);
#endif

                                }


#else
                                else
                                {
#if DIMD_NUM_INTRA_DIR_INC
                                  rpcTempCU->setExtIntraDirSubParts(CHANNEL_TYPE_LUMA, uiIntraMode, 0, uiDepth);
                                  rpcTempCU->setExtIntraDirSubParts(CHANNEL_TYPE_CHROMA, uiIntraMode, 0, uiDepth);
#if VCEG_AZ07_INTRA_65ANG_MODES
                                  rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO67(uiIntraMode), 0, uiDepth);
#else
                                  rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, MAP131TO35(uiIntraMode), 0, uiDepth);
#endif
#else
                                  rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiIntraMode, 0, uiDepth);
#endif
#if	DIMD_DISABLE_CHROMA_P_B_SLICE_JB

#else
                                  rpcTempCU->setIntraDirSubParts(CHANNEL_TYPE_CHROMA, DM_CHROMA_IDX, 0, uiDepth);
#endif
                                }
#endif
                              }

#endif
#if DIMD_INTRA_PRED
                              /* DIMD processing end */
                              if (btestDIMD[dimdIdx])
#endif
                              {
                                xCheckRDCostIntra(rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug)
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                                  , bNonZeroCoeff
#endif
#if DIMD_INTRA_PRED
                                  , dimdIdx
#endif
// #if DIMD_MLIP_HARM
//                                     , ((dimdIdx == 1) && btestDIMD[dimdIdx])
// #endif
                                  );
                              }
                              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#endif

#if RRC_DWDIP
                        }
#endif

#if VCEG_AZ05_INTRA_MPI
                            if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                          }
#endif
#if RRC_DWDIP_MLR
                  } // for (iDWDIP_GP = 0; ...)
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
                        if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                }
#endif
#if RRC_INTRA_SHARPENING
                      // if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break; // TODO: check this
              }
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                    if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                    //if (bNonZeroCoeff>rpcTempCU->getWidth(0)*rpcTempCU->getHeight(0) && rpcTempCU->getSlice()->getSliceType() == I_SLICE) break;
            }
#endif

#if DIMD_INTRA_PRED
          }
#endif

#if !JVET_C0024_QTBT
                if (uiDepth == sps.getLog2DiffMaxMinCodingBlockSize()
#if COM16_C806_EMT
                  && !bEarlySkipIntra
#endif
                  )
                {
#if COM16_C806_EMT && (VCEG_AZ05_ROT_TR ||VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST|| COM16_C1046_PDPC_INTRA)
                  UChar ucEmtUsage = ((rpcTempCU->getWidth(0) > EMT_INTRA_MAX_CU) || (rpcTempCU->getSlice()->getSPS()->getUseIntraEMT() == 0)) ? 1 : 2;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                  for (iROTidx = 0; iROTidx < iNumberOfPassesROT; iROTidx++)
                  {
#endif
#if RRC_INTRA_SHARPENING
                    for (Int iSharpIdx = 0; iSharpIdx < 2; iSharpIdx++)
                    {
#endif
#if VCEG_AZ05_INTRA_MPI
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                      iMPIidx = 0;   iNumberOfPassesMPI = 2;
                      if (rpcTempCU->getSlice()->getSliceType() != I_SLICE)  iNumberOfPassesMPI = 2;
                      if (!rpcTempCU->getSlice()->getSPS()->getUseMPI()) iNumberOfPassesMPI = 1;
                      if (iROTidx) iNumberOfPassesMPI = 1;
#endif
                      for (iMPIidx = 0; iMPIidx < iNumberOfPassesMPI; iMPIidx++)
                      {
#endif 

#if COM16_C1046_PDPC_INTRA
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                        iPDPCidx = 0;   iNumberOfPassesPDPC = 2;
                        if (rpcTempCU->getSlice()->getSliceType() != I_SLICE)  iNumberOfPassesPDPC = 2;
                        if (!rpcTempCU->getSlice()->getSPS()->getUsePDPC()) iNumberOfPassesPDPC = 1;
                        if (iROTidx) iNumberOfPassesPDPC = 1;
#endif
                        for (iPDPCidx = 0; iPDPCidx<iNumberOfPassesPDPC; iPDPCidx++)
                        {
#endif 

                          if (rpcTempCU->getWidth(0) >(1 << sps.getQuadtreeTULog2MinSize()))
                          {
                            Double tmpIntraCost;
#if COM16_C806_EMT
                            for (UChar ucCuFlag = 0; ucCuFlag < ucEmtUsage; ucCuFlag++)
                            {
                              if (ucCuFlag && bAllIntra && m_pcEncCfg->getUseFastIntraEMT())
                              {
                                static const Double thEmtIntraFastSkipNxN = 1.2; // Skip checking "NxN using EMT" if "NxN using DCT2" is worse than "2Nx2N using DCT2"
                                if (dIntraNxNCost > thEmtIntraFastSkipNxN*dIntra2Nx2NCost)
                                {
                                  break;
                                }
                              }
                              rpcTempCU->setEmtCuFlagSubParts(ucCuFlag, 0, uiDepth);
#if VCEG_AZ05_INTRA_MPI
                              rpcTempCU->setMPIIdxSubParts(iMPIidx, 0, uiDepth);
#endif
#if COM16_C1046_PDPC_INTRA
                              rpcTempCU->setPDPCIdxSubParts(iPDPCidx, 0, uiDepth);
#endif
#if MLIP_EXTEND_REF_LINE
                              rpcTempCU->setExternRefSubParts(rpcTempCU->getTextType(), 0, 0, uiDepth);
#endif
#if RRC_DWDIP
                              rpcTempCU->setDWDIPIdxSubParts(0, 0, uiDepth);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                              rpcTempCU->setROTIdxSubParts(iROTidx, 0, uiDepth);
#endif
#if RRC_INTRA_SHARPENING
                              rpcTempCU->setIntraSharpFlagSubParts(iSharpIdx, 0, uiDepth);
#endif
#if INTRA_SHARP_PDPC_HARM
                              if (iSharpIdx && iPDPCidx)
                                continue;
#endif
#if INTRA_SHARP_MAX_ASPECT_RATIO
                              int AspectRatio = std::max(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0)) / std::min(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0));
                              if (iSharpIdx && AspectRatio > INTRA_SHARP_MAX_ASPECT_RATIO)
                                continue;
#endif
#if INTRA_SHARP_MAX_SIZE
                              if (iSharpIdx && std::max(rpcTempCU->getWidth(0), rpcTempCU->getHeight(0)) > INTRA_SHARP_MAX_SIZE)
                                continue;
#endif
                              xCheckRDCostIntra(rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug)
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                                , bNonZeroCoeff
#endif
                                );
                              if (!ucCuFlag && m_pcEncCfg->getUseFastIntraEMT())
                              {
                                dIntraNxNCost = (rpcBestCU->isIntra(0) && rpcBestCU->getPartitionSize(0) == SIZE_NxN) ? rpcBestCU->getTotalCost() : rpcTempCU->getTotalCost();
                              }
                              intraCost = std::min(intraCost, tmpIntraCost);
                              rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                            }
#else
#if VCEG_AZ05_INTRA_MPI
                            rpcTempCU->setMPIIdxSubParts(iMPIidx, 0, uiDepth);
#endif
#if COM16_C1046_PDPC_INTRA
                            rpcTempCU->setPDPCIdxSubParts(iPDPCidx, 0, uiDepth);
#endif
#if MLIP_EXTEND_REF_LINE
                            rpcTempCU->setExternRefSubParts(rpcTempCU->getTextType(), 0, 0, uiDepth);
#endif
#if RRC_DWDIP
                            rpcTempCU->setDWDIPIdxSubParts(0, 0, uiDepth);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
                            rpcTempCU->setROTIdxSubParts(iROTidx, 0, uiDepth);
#endif
#if RRC_INTRA_SHARPENING
                            rpcTempCU->setIntraSharpFlagSubParts(iSharpIdx, 0, uiDepth);
#endif
                            xCheckRDCostIntra(rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug)
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                              , bNonZeroCoeff
#endif
                              );
                            intraCost = std::min(intraCost, tmpIntraCost);
                            rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
#endif
                          }
#if VCEG_AZ05_INTRA_MPI
                          if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                        }
#endif
#if COM16_C1046_PDPC_INTRA
                        if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                      }
#endif
#if RRC_INTRA_SHARPENING
                      //if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break; // TODO: check this
                    }
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
                    if (rpcBestCU->isIntra(0) && !bNonZeroCoeff) break;
                  }
#endif
                }
#endif //#if !JVET_C0024_QTBT
        }

              // test PCM
              if (sps.getUsePCM()
                && rpcTempCU->getWidth(0) <= (1 << sps.getPCMLog2MaxSize())
                && rpcTempCU->getWidth(0) >= (1 << sps.getPCMLog2MinSize())
#if BARC
                && !isBARCFlag
#endif
 )
              {
                UInt uiRawBits = getTotalBits(rpcBestCU->getWidth(0), rpcBestCU->getHeight(0), rpcBestCU->getPic()->getChromaFormat(), sps.getBitDepths().recon);
                UInt uiBestBits = rpcBestCU->getTotalBits();
                if ((uiBestBits > uiRawBits) || (rpcBestCU->getTotalCost() > m_pcRdCost->calcRdCost(uiRawBits, 0)))
                {
                  xCheckIntraPCM(rpcBestCU, rpcTempCU);
                  rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
                }
              }

              if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
              {
                iQP = iMinQP;
              }
      }
    }

#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
          if (sps.getUseInterKLT())
          {
#endif
            if (!rpcBestCU->isIntra(0) && rpcBestCU->getQtRootCbf(0) != 0)
            {
              //Only check from the best modes for speeding up
              g_bEnableCheck = true;
              Int iQP = rpcBestCU->getQP(0);
#if JVET_C0024_QTBT
              xCheckRDCostInterKLT(rpcBestCU, rpcTempCU, SIZE_2Nx2N);
#else
              PartSize eSize = rpcBestCU->getPartitionSize(0);
              xCheckRDCostInterKLT(rpcBestCU, rpcTempCU, eSize);
#endif
              rpcTempCU->initEstData(uiDepth, iQP, false);
            }
#if VCEG_AZ08_USE_KLT
          }
#endif
#endif
#if JVET_C0024_QTBT
          if (rpcBestCU->getTotalCost() < MAX_DOUBLE - 1)
#else
          if (rpcBestCU->getTotalCost() != MAX_DOUBLE)
#endif
          {
#if JVET_C0024_QTBT
            m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_NEXT_BEST]);
#else
            m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
#endif
            m_pcEntropyCoder->resetBits();
#if JVET_C0024_QTBT
            if (uiBTSplitMode == 0)
            {
#endif
              m_pcEntropyCoder->encodeSplitFlag(rpcBestCU, 0, uiDepth, true);
#if JVET_C0024_QTBT
            }
#if JVET_C0024_SPS_MAX_BT_DEPTH
            UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? sps.getMaxBTDepthISliceL() : sps.getMaxBTDepthISliceC()) : sps.getMaxBTDepth();
#if BT_FORCE_EXT
            Int iMaxELv = 
#if BT_FORCE_INTRA_ON
              pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
              FORCE_BTD_EXT;
            if (
#if !BT_FORCE_INTRA_ON
              !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
              !rpcBestCU->isFromForcedQT(0) &&
#endif
              rpcBestCU->getBoundaryCtuLocation() && iMaxELv)
            {
              UInt uiMaxBTDTmp = rpcBestCU->getSlice()->getMaxFQTLevel(rpcBestCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
              uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
            }
#endif 
#else
            UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? MAX_BT_DEPTH : MAX_BT_DEPTH_C) : MAX_BT_DEPTH_INTER;
#if BT_FORCE_EXT
            Int iMaxELv = 
#if BT_FORCE_INTRA_ON
              pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
              FORCE_BTD_EXT;
            if (
#if !BT_FORCE_INTRA_ON
              !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
              !rpcBestCU->isFromForcedQT(0) &&
#endif
              rpcBestCU->getBoundaryCtuLocation() && iMaxELv)
            {
              UInt uiMaxBTDTmp = rpcBestCU->getSlice()->getMaxFQTLevel(rpcBestCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
              uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
            }
#endif 
#endif
#if JVET_C0024_SPS_MAX_BT_SIZE
            UInt uiMaxBTSize = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? sps.getMaxBTSizeISliceL() : sps.getMaxBTSizeISliceC()) : sps.getMaxBTSize();
#else
            UInt uiMaxBTSize = isLuma(rpcTempCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C;
#endif
            UInt uiMinBTSize = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? MIN_BT_SIZE : MIN_BT_SIZE_C) : MIN_BT_SIZE_INTER;

            if (rpcBestCU->getWidth(0) <= uiMaxBTSize && rpcBestCU->getHeight(0) <= uiMaxBTSize
              && (rpcBestCU->getWidth(0) > uiMinBTSize || rpcBestCU->getHeight(0) > uiMinBTSize)
              && rpcBestCU->getBTDepth(0) < uiMaxBTD)
            {
              m_pcEntropyCoder->encodeBTSplitMode(rpcBestCU, 0, rpcBestCU->getWidth(0), rpcBestCU->getHeight(0), true);
            }
#endif
            rpcBestCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
            rpcBestCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
            rpcBestCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcBestCU->getTotalBits(), rpcBestCU->getTotalDistortion());
#if TT_ENC_OPT //store rd cost
            dRDCost[0] = rpcBestCU->getTotalCost();
            dRDDist[0] = rpcBestCU->getTotalDistortion();
            iRDBits[0] = rpcBestCU->getTotalBits();
#endif
#if JVET_C0024_QTBT
            m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_NEXT_BEST]);
#else
            m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
            if (saveLoadTag == SAVE_ENC_INFO)
            {
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
#if COM16_C806_EMT
              m_pcPredSearch->setSaveLoadEmtFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getEmtCuFlag(0));
              m_pcPredSearch->setSaveLoadEmtIdx(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getEmtTuIdx(0));
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
              m_pcPredSearch->setSaveLoadRotIdx(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getROTIdx(rpcBestCU->getTextType(), 0));
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
              m_pcPredSearch->setSaveLoadPdpcIdx(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getPDPCIdx(0));
#endif
              if (!rpcBestCU->isIntra(0))
              {
#if VCEG_AZ07_FRUC_MERGE
                m_pcPredSearch->setSaveLoadFrucMode(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getFRUCMgrMode(0));
#endif
#if VCEG_AZ07_IMV
                m_pcPredSearch->setSaveLoadIMVFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getiMVFlag(0));
#endif
#if VCEG_AZ06_IC
                m_pcPredSearch->setSaveLoadICFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getICFlag(0));
#endif
                m_pcPredSearch->setSaveLoadMergeFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getMergeFlag(0));
#if COM16_C1016_AFFINE
                m_pcPredSearch->setSaveLoadAffineFlag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getAffineFlag(0));
#endif
                m_pcPredSearch->setSaveLoadInterDir(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, rpcBestCU->getInterDir(0));
              }
              dNonSplitCost = dCostTempBest = rpcBestCU->getTotalCost();
              m_pcPredSearch->setSaveLoadTag(uiZorderIdx, rpcBestCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx, LOAD_ENC_INFO);
#else
#if COM16_C806_EMT
              m_pcPredSearch->setSaveLoadEmtFlag(uiWidthIdx, uiHeightIdx, rpcBestCU->getEmtCuFlag(0));
              m_pcPredSearch->setSaveLoadEmtIdx(uiWidthIdx, uiHeightIdx, rpcBestCU->getEmtTuIdx(0));
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
              m_pcPredSearch->setSaveLoadRotIdx(uiWidthIdx, uiHeightIdx, rpcBestCU->getROTIdx(rpcBestCU->getTextType(), 0));
#endif
#if RRC_INTRA_SHARPENING && INTRA_SHARP_D0077
              m_pcPredSearch->setSaveLoadSharpIdx(uiWidthIdx, uiHeightIdx, rpcBestCU->getIntraSharpFlag(0));
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
              m_pcPredSearch->setSaveLoadPdpcIdx(uiWidthIdx, uiHeightIdx, rpcBestCU->getPDPCIdx(0));
#endif
              if (!rpcBestCU->isIntra(0))
              {
#if VCEG_AZ07_FRUC_MERGE
                m_pcPredSearch->setSaveLoadFrucMode(uiWidthIdx, uiHeightIdx, rpcBestCU->getFRUCMgrMode(0));
#endif
#if VCEG_AZ07_IMV
                m_pcPredSearch->setSaveLoadIMVFlag(uiWidthIdx, uiHeightIdx, rpcBestCU->getiMVFlag(0));
#endif
#if VCEG_AZ06_IC
                m_pcPredSearch->setSaveLoadICFlag(uiWidthIdx, uiHeightIdx, rpcBestCU->getICFlag(0));
#endif
                m_pcPredSearch->setSaveLoadMergeFlag(uiWidthIdx, uiHeightIdx, rpcBestCU->getMergeFlag(0));
#if COM16_C1016_AFFINE
                m_pcPredSearch->setSaveLoadAffineFlag(uiWidthIdx, uiHeightIdx, rpcBestCU->getAffineFlag(0));
#endif
                m_pcPredSearch->setSaveLoadInterDir(uiWidthIdx, uiHeightIdx, rpcBestCU->getInterDir(0));
              }
              dNonSplitCost = dCostTempBest = rpcBestCU->getTotalCost();
              m_pcPredSearch->setSaveLoadTag(uiZorderIdx, uiWidthIdx, uiHeightIdx, LOAD_ENC_INFO);
#endif
            }
#endif
          }
#if JVET_C0024_DELTA_QP_FIX
          if (pps.getUseDQP())
          {
            rpcBestCU->setCodedQP(rpcBestCU->getQP(0));
          }
#endif
  }

#if JVET_C0024_QTBT && COM16_C806_LARGE_CTU
        if (ucMinDepth > uiDepth)
        {
          bBoundary = true; //to force not to try BT split. JCA
        }
#endif
        // copy original YUV samples to PCM buffer
        if (rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isLosslessCoded(0) && (rpcBestCU->getIPCMFlag(0) == false))
        {
#if JVET_C0024_QTBT
          xFillPCMBuffer(rpcBestCU, m_pppcOrigYuv[uiWidthIdx][uiHeightIdx]);
#else
          xFillPCMBuffer(rpcBestCU, m_ppcOrigYuv[uiDepth]);
#endif
        }

#if JVET_C0024_DELTA_QP_FIX
        if (uiQTBTDepth == uiMaxDQPDepthQTBT)
#else
        if (uiDepth == pps.getMaxCuDQPDepth())
#endif
        {
          Int idQP = m_pcEncCfg->getMaxDeltaQP();
          iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - idQP);
          iMaxQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP + idQP);
        }
#if JVET_C0024_DELTA_QP_FIX
        else if (uiQTBTDepth < uiMaxDQPDepthQTBT)
#else
        else if (uiDepth < pps.getMaxCuDQPDepth())
#endif
        {
          iMinQP = iBaseQP;
          iMaxQP = iBaseQP;
        }
        else
        {
          const Int iStartQP = rpcTempCU->getQP(0);
          iMinQP = iStartQP;
          iMaxQP = iStartQP;
        }

        if (m_pcEncCfg->getUseRateCtrl())
        {
          iMinQP = m_pcRateCtrl->getRCQP();
          iMaxQP = m_pcRateCtrl->getRCQP();
        }

        if (m_pcEncCfg->getCUTransquantBypassFlagForceValue())
        {
          iMaxQP = iMinQP; // If all TUs are forced into using transquant bypass, do not loop here.
        }

        Bool bSubBranch = bBoundary ||
#if COM16_C806_LARGE_CTU
          (!(m_pcEncCfg->getUseEarlyCU() && rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isSkipped(0)) && (!m_pcEncCfg->getUseFastLCTU() || uiDepth < ucMaxDepth));
#else
          !(m_pcEncCfg->getUseEarlyCU() && rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isSkipped(0));
#endif

#if JVET_C0024_BT_RMV_REDUNDANT
        bQTreeValid = false;
        if (uiBTSplitMode == 0)
        {
          UInt uiMinQTSize = sps.getMinQTSize(rpcBestCU->getSlice()->getSliceType(), rpcBestCU->getTextType());
          bQTreeValid = (uiWidth >= 2 * uiMinQTSize);
        }

        Bool bBTHorRmvEnable = false;
        Bool bBTVerRmvEnable = false;
        if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
        {
          bBTHorRmvEnable = true;
          bBTVerRmvEnable = bQTreeValid;
        }
#endif

#if JVET_C0024_QTBT
        Bool bQTSplit = bSubBranch && uiBTSplitMode == 0 && (!getFastDeltaQp() || uiWidth > fastDeltaQPCuMaxSize || bBoundary);
        if (bQTSplit)
        {
          assert(uiWidth == uiHeight);
        }
        bQTSplit = bQTSplit && (uiWidth > sps.getMinQTSize(pcSlice->getSliceType(), rpcBestCU->getTextType()) || bBoundary);

#if BT_FORCE
#if !BT_FORCE_INTRA_ON
        if (!pcSlice->isIntra())
#endif
          bQTSplit = bQTSplit && (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv));
#endif

        if (!bBoundary && rpcBestCU->isSkipped(0))
        {
          rpcBestCU->getPic()->setSkiped(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight, true);
        }
        if (!bBoundary && rpcBestCU->getPredictionMode(0) == MODE_INTER)
        {
          rpcBestCU->getPic()->setInter(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight, true);
        }
        else if (!bBoundary && rpcBestCU->getPredictionMode(0) == MODE_INTRA)
        {
          rpcBestCU->getPic()->setIntra(rpcBestCU->getZorderIdxInCtu(), uiWidth, uiHeight, true);
        }

#if !JVET_C0024_DELTA_QP_FIX
        UInt uiQTWidth = sps.getCTUSize() >> uiDepth;
        UInt uiQTHeight = sps.getCTUSize() >> uiDepth;
        UInt uiBTDepth = g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight];
#endif

#if JVET_C0024_SPS_MAX_BT_DEPTH
        UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? sps.getMaxBTDepthISliceL() : sps.getMaxBTDepthISliceC()) : sps.getMaxBTDepth();
#if BT_FORCE_EXT
        Int iMaxELv = 
#if BT_FORCE_INTRA_ON
          pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
          FORCE_BTD_EXT;
        if (
#if !BT_FORCE_INTRA_ON
          !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
          !rpcTempCU->isFromForcedQT(0) &&
#endif
          rpcTempCU->getBoundaryCtuLocation() && iMaxELv)
        {
          UInt uiMaxBTDTmp = rpcTempCU->getSlice()->getMaxFQTLevel(rpcTempCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
          uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
        }
#endif 
#else
#if BT_FORCE_EXT
        Int iMaxELv = 
#if BT_FORCE_INTRA_ON
          pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
          FORCE_BTD_EXT;
        if (
#if !BT_FORCE_INTRA_ON
          !pcSlice->isIntra() &&
#endif
#if BT_FORCE_INTRA_SPEED_UP
          !rpcTempCU->isFromForcedQT(0) &&
#endif
          rpcTempCU->getBoundaryCtuLocation() && iMaxELv)
        {
          UInt uiMaxBTDTmp = rpcTempCU->getSlice()->getMaxFQTLevel(rpcTempCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
          uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
        }
#endif 
        UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? MAX_BT_DEPTH : MAX_BT_DEPTH_C) : MAX_BT_DEPTH_INTER;
#endif
#if JVET_C0024_SPS_MAX_BT_SIZE
        UInt uiMaxBTSize = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? sps.getMaxBTSizeISliceL() : sps.getMaxBTSizeISliceC()) : sps.getMaxBTSize();
#else
        UInt uiMaxBTSize = isLuma(rpcTempCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C;
#endif
        UInt uiMinBTSize = pcSlice->isIntra() ? (isLuma(rpcTempCU->getTextType()) ? MIN_BT_SIZE : MIN_BT_SIZE_C) : MIN_BT_SIZE_INTER;

        Bool bTestHorSplit = (!bBoundary && uiHeight > uiMinBTSize
          && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD
#if JVET_C0024_QTBT
          && !bForceQT
#endif
          );

        Bool bTestVerSplit = (!bBoundary && uiWidth > uiMinBTSize
          && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD
#if JVET_C0024_QTBT
          && !bForceQT
#endif
          );

#if TT_ENC_OPT
        Bool bTestHorSplitTT = (!bBoundary && uiHeight > (uiMinBTSize << 1)
          && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD
          && !bForceQT);

        Bool bTestVerSplitTT = (!bBoundary && uiWidth > (uiMinBTSize << 1)
          && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD
          && !bForceQT);

        UInt uiParentType = (uiBTDepth == 0) ? 0 : rpcBestCU->getBTSplitModeForBTDepth(0, uiBTDepth - 1);
#endif

        /*-- TT or BT off switch here ---*/


#if (TT_LIMIT_AR>0) //normative
        if (uiWidth >= uiHeight*TT_LIMIT_AR)
          bTestHorSplit = bTestHorSplitTT = false;
        if (uiHeight >= uiWidth*TT_LIMIT_AR)
          bTestVerSplit = bTestVerSplitTT = false;
#endif
#if TT_LIMIT_SIZEIDX //normative
        UInt uiDeltaSizeIdx = g_aucConvertToBit[uiQTWidth] - g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight] - g_aucConvertToBit[uiHeight];
        if (uiDeltaSizeIdx >= (uiMaxBTD + 1))
          bTestHorSplit = bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;
#endif
#if TT_LEAF_NODE //normative
        UInt uiBTSplitModeLastBtDepth = (uiBTDepth == 0) ? 0 : ((uiBTSplitMode >> ((uiBTDepth - 1) << 2)) & 15);
        if (uiBTSplitModeLastBtDepth & 2) //if tt
          bTestHorSplit = bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;
#endif
#if TT_LIMIT_SMALL || TT_LIMIT_LARGE //normative
        if (TT_LIMIT_SMALL)
        {
          if (uiParentType == ST_HORI_TT_S)
            bTestHorSplitTT = false;
          else if (uiParentType == ST_VERT_TT_S)
            bTestVerSplitTT = false;
        }
        if (TT_LIMIT_LARGE)
        {
          if (uiParentType == ST_HORI_TT_L)
            bTestHorSplitTT = false;
          else if (uiParentType == ST_VERT_TT_L)
            bTestVerSplitTT = false;
        }
#endif

#if ENC_FAST_SKIP_RDC_EQUAL
        Bool bSkipSecondBT = false;
        Bool bSkipSecondTT = false;
#endif
#if TT_FAST_USE_ONE_TT
        UInt uiSkipTTbyRD = 0;
#endif
#if TT_FAST_SMALL_ES //decide
        if (uiBTDepth == 0)
        {
          rpcTempCU->setSTSplitModeSave(0, 0, 0);
          rpcBestCU->setSTSplitModeSave(0, 0, 0);
          rpcTempCU->setSTSplitModeSave(1, 0, 0);
          rpcBestCU->setSTSplitModeSave(1, 0, 0);
        }

        if (uiBTDepth == 1 && (uiParentType == ST_HORI_TT_S || uiParentType == ST_VERT_TT_S) && rpcBestCU->getSlice()->getSliceType() != I_SLICE)
        {
          Bool bTryCopyMode = false;
          UInt uiSplitModeD0, uiSplitModeD1, uiSplitModeD2;
          UInt uiBTTLMode = uiParentType == ST_HORI_TT_S ? rpcBestCU->getSTSplitModeSave(0, 0) : rpcBestCU->getSTSplitModeSave(1, 0);
          UInt uiBTBRMode = uiParentType == ST_HORI_TT_S ? rpcBestCU->getSTSplitModeSave(0, 1) : rpcBestCU->getSTSplitModeSave(1, 1);

          if (rpcBestCU->getSubCUIdx() == 0 && uiBTTLMode != 0)
          {
            uiSplitModeD0 = uiBTTLMode & 15;
            uiSplitModeD1 = (uiBTTLMode >> 4) & 15;
            uiSplitModeD2 = (uiBTTLMode >> 8) & 15;
            bTryCopyMode = true;
          }
          else if (rpcBestCU->getSubCUIdx() == 2 && uiBTBRMode != 0)
          {
            uiSplitModeD0 = uiBTBRMode & 15;
            uiSplitModeD1 = (uiBTBRMode >> 4) & 15;
            uiSplitModeD2 = (uiBTBRMode >> 8) & 15;
            bTryCopyMode = true;
          }

          if (bTryCopyMode)
          {
            //part 1 : if BT CU not split, TT CU not split
            if (uiSplitModeD1 == ST_NO_SPLIT)
              bTestHorSplit = bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;

            //part 2 : if BT CU are same size as the TT CU, copy the split mode
            if ((uiSplitModeD0 == ST_HORI_BT && uiSplitModeD1 == ST_HORI_BT) || (uiSplitModeD0 == ST_VERT_BT && uiSplitModeD1 == ST_VERT_BT))
            {
              //force the same split
              if (uiSplitModeD2 == ST_NO_SPLIT)
                bTestHorSplit = bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;
              else if (uiSplitModeD2 == ST_HORI_BT)
                bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;
              else if (uiSplitModeD2 == ST_VERT_BT)
                bTestHorSplit = bTestHorSplitTT = bTestVerSplitTT = false;
              else if (uiSplitModeD2 == ST_HORI_TT_S)
                bTestHorSplit = bTestVerSplit = bTestVerSplitTT = false;
              else if (uiSplitModeD2 == ST_VERT_TT_S)
                bTestHorSplit = bTestVerSplit = bTestHorSplitTT = false;
              else
                assert(0);
            }
          }
        }
#endif
#if TT_FAST_LARGE_ES //decide
        if (uiBTDepth == 0)
        {
          rpcTempCU->setSTSplitModeSaveL(0, 0, 0);
          rpcBestCU->setSTSplitModeSaveL(0, 0, 0);
          rpcTempCU->setSTSplitModeSaveL(1, 0, 0);
          rpcBestCU->setSTSplitModeSaveL(1, 0, 0);
        }

        if (uiBTDepth == 1 && (uiParentType == ST_HORI_TT_L || uiParentType == ST_VERT_TT_L) && rpcBestCU->getSlice()->getSliceType() != I_SLICE)
        {
          Bool bTryCopyMode = false;
          UInt uiSplitModeC1D0, uiSplitModeC1D1, uiSplitModeC1D2;
          UInt uiSplitModeC2D0, uiSplitModeC2D1, uiSplitModeC2D2;
          UInt uiBTTLMode = uiParentType == ST_HORI_TT_L ? rpcBestCU->getSTSplitModeSave(0, 2) : rpcBestCU->getSTSplitModeSave(1, 2);
          UInt uiBTBRMode = uiParentType == ST_HORI_TT_L ? rpcBestCU->getSTSplitModeSave(0, 3) : rpcBestCU->getSTSplitModeSave(1, 3);

          if (rpcBestCU->getSubCUIdx() == 1 && uiBTTLMode != 0)
          {
            uiSplitModeC1D0 = uiBTTLMode & 15;
            uiSplitModeC1D1 = (uiBTTLMode >> 4) & 15;
            uiSplitModeC1D2 = (uiBTTLMode >> 8) & 15;
            uiSplitModeC2D0 = uiBTBRMode & 15;
            uiSplitModeC2D1 = (uiBTBRMode >> 4) & 15;
            uiSplitModeC2D2 = (uiBTBRMode >> 8) & 15;
            bTryCopyMode = true;
          }

          if (bTryCopyMode)
          {
            //part 1 : if two BT CUs not split, TT CU try one BT split 
            if (uiSplitModeC1D1 == ST_NO_SPLIT && uiSplitModeC2D1 == ST_NO_SPLIT)
            {
              if (uiSplitModeC1D0 == ST_HORI_BT) //try HBT
                bTestVerSplit = bTestHorSplitTT = bTestVerSplitTT = false;
              else if (uiSplitModeC1D0 == ST_VERT_BT) //try VBT
                bTestHorSplit = bTestHorSplitTT = bTestVerSplitTT = false;
            }
          }
        }
#endif

        //for encoder speedup
#if JVET_E0023_FAST_ENCODING_SETTING
#if TT_ENC_OPT
        if (rpcBestCU->getSkipFlag(0) && (bTestHorSplit || bTestVerSplit || bTestHorSplitTT || bTestVerSplitTT) && uiBTDepth >= ((pcSlice->getPictureDistance() <= PICTURE_DISTANCE_TH) ? FAST_SKIP_DEPTH : SKIP_DEPTH))
#else
        if (rpcBestCU->getSkipFlag(0) && (bTestHorSplit || bTestVerSplit) && uiBTDepth >= ((pcSlice->getPictureDistance() <= PICTURE_DISTANCE_TH) ? FAST_SKIP_DEPTH : SKIP_DEPTH))
#endif
#else
        if (rpcBestCU->getSkipFlag(0) && (bTestHorSplit || bTestVerSplit) && uiBTDepth >= SKIP_DEPTH)
#endif
        {
          bTestHorSplit = bTestVerSplit = bQTSplit = false;
#if TT_ENC_OPT
          bTestHorSplitTT = bTestVerSplitTT = false;
#endif
        }

#if JVET_C0024_BT_RMV_REDUNDANT
        if (uiSplitConstrain == 1)
        {
          bTestHorSplit = false;
        }
        if (uiSplitConstrain == 2)
        {
          bTestVerSplit = false;
        }
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO && !JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
        if (bUseSaveLoad)
        {
          m_pcPredSearch->setSaveLoadTag(uiZorderIdx, uiWidthIdx - 1, uiHeightIdx - 1, SAVE_ENC_INFO);
        }
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO && !JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
        if (bUseSaveLoadSplitDecision && saveLoadTag == LOAD_ENC_INFO)
        {
          if (bTestHorSplit && (saveLoadSplit & 0x02))
          {
            bTestHorSplit = false;
          }
          if (bTestVerSplit && (saveLoadSplit & 0x04))
          {
            bTestVerSplit = false;
          }
        }
#endif

#if TT_ENC_OPT
#if TT_BUGFIX_ZIDX
        Int iC1x = uiWidth >> (MIN_CU_LOG2+1);
        Int iC1y = uiHeight>> (MIN_CU_LOG2+1);
#else
        Int iC1x = uiWidth >> (MIN_CU_LOG2 + 2);
        Int iC1y = 0;
#endif
        Int iC2x = (uiWidth >> (MIN_CU_LOG2 + 1)) + (uiWidth >> (MIN_CU_LOG2 + 2)) - 1;
        Int iC2y = (uiHeight >> MIN_CU_LOG2) - 1;
        iC2x = iC2x<0 ? 0 : iC2x;
        iC2y = iC2y<0 ? 0 : iC2y;
        UInt uiZorderTL = 0; //top left
        UInt uiZorderC1 = g_auiRasterToZscan[iC1y * (sps.getCTUSize() >> MIN_CU_LOG2) + iC1x]; //center, closer to TL
        UInt uiZorderC2 = g_auiRasterToZscan[iC2y * (sps.getCTUSize() >> MIN_CU_LOG2) + iC2x]; //center, closer to BR
        UInt uiZorderBR = g_auiRasterToZscan[((uiHeight >> MIN_CU_LOG2) - 1) * (sps.getCTUSize() >> MIN_CU_LOG2) + (uiWidth >> MIN_CU_LOG2) - 1];   //bottom-right part.
        Bool bIsPBSlice = pcSlice->getSliceType() != I_SLICE;
#endif
        Int aiParam[7];
        Double adParam[3];
        Bool bTestFirstBTSplit, bTestSecondBTSplit;
        UInt iFirstBTDir, iSecondBTDir;
        UInt uiFirstBTSplitIdx, uiSecondBTSplitIdx;

#if BT_FORCE
          bTestHorSplit = bTestHorSplit || (
#if !BT_FORCE_INTRA_ON
            !rpcBestCU->getSlice()->isIntra() &&
#endif
            (
#if BT_FORCE_TT
              rpcTempCU->isForcedBT(0, uiWidth, uiHeight)&&
#endif
              bBBoundary && (uiDepth >= iMaxQtLv)));
        bTestVerSplit = bTestVerSplit || (
#if !BT_FORCE_INTRA_ON
          !rpcBestCU->getSlice()->isIntra() && 
#endif
          (
#if BT_FORCE_TT
            rpcTempCU->isForcedBT(0, uiWidth, uiHeight) &&
#endif
            bRBoundary && (uiDepth >= iMaxQtLv)));
#endif
#if (BT_FORCE_AR > 0)
        if (bBBoundary && bTestHorSplit)
        {
          if (uiWidth >= uiHeight*BT_FORCE_AR)
          {
            bTestHorSplit = false;
            bTestVerSplit = true;
          }
        }
        if (bRBoundary && bTestVerSplit)
        {
          if (uiHeight >= uiWidth*BT_FORCE_AR)
          {
            bTestVerSplit = false;
            bTestHorSplit = true;
          }
        }
#endif


#if ADAPTIVE_RD_ORDER
        if (uiWidth < uiHeight) {
#endif
          bTestFirstBTSplit = bTestHorSplit;
          bTestSecondBTSplit = bTestVerSplit;
          iFirstBTDir = 0;
          iSecondBTDir = 1;
          uiFirstBTSplitIdx = 1;
          uiSecondBTSplitIdx = 2;
#if ADAPTIVE_RD_ORDER
  }else{
          bTestFirstBTSplit = bTestVerSplit;
          bTestSecondBTSplit = bTestHorSplit;
          iFirstBTDir = 1;
          iSecondBTDir = 0;
          uiFirstBTSplitIdx = 2;
          uiSecondBTSplitIdx = 1;
        }
#endif

        if (bTestFirstBTSplit)
        {
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL && SIF_QTBT_FIX
#if JVET_C0024_QTBT
          UInt uiQTWidthIdx = g_aucConvertToBit[uiWidth >> 1];
          UInt uiQTHeightIdx = g_aucConvertToBit[uiHeight >> 1];
#else
          PartSize parentSize = rpcBestCU->getPartitionSize(0);
#endif
          int modeSIFnum = 1;
          if (rpcTempCU->isCodeBIF(uiWidth, /*uiDepth*/uiQTBTDepth
#if SCU_SIF
            , m_pcEntropyCoder->getEntropyCoder()->getSigDepth(0)
#endif        
            ))
          {
            modeSIFnum = 2;
          }
          for (int modeSIF = 0; modeSIF < modeSIFnum; modeSIF++)
          {
            if (modeSIFnum == 2)
            {
              rpcTempCU->setSIFToCheck(modeSIF);
              rpcBestCU->setSIFToCheck(modeSIF);
            }
#endif //RRC_SELECTABLE_IF
            // further split
            for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
            {
              //prepare parameters
              aiParam[0] = iFirstBTDir;
              aiParam[1] = iQP;
              aiParam[2] = lastCodedQP;
              aiParam[3] = bBTHorRmvEnable;
              aiParam[4] = bBTVerRmvEnable;
#if JVET_D0077_SAVE_LOAD_ENC_INFO
              if (iFirstBTDir == 0)
                adParam[2] = dHorSplitCost;
              else
                adParam[2] = dVerSplitCost;
#endif

              CompressBT(rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight, uiBTSplitMode, aiParam, adParam, uiSplitConstrain);
#if TT_ENC_OPT
              dRDCost[uiFirstBTSplitIdx] = adParam[0];
              dRDDist[uiFirstBTSplitIdx] = adParam[1];
              iRDBits[uiFirstBTSplitIdx] = aiParam[5];
#endif
#if ENC_FAST_SKIP_RDC_EQUAL
              bSkipSecondBT = aiParam[6];
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
              if (iFirstBTDir == 0)
                dHorSplitCost = adParam[2];
              else
                dVerSplitCost = adParam[2];
#endif
            }
#if JVET_D0077_SAVE_LOAD_ENC_INFO
            if (dCostTempBest > adParam[2])
              dCostTempBest = adParam[2];
#endif
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL && SIF_QTBT_FIX
          }
#endif
        }

        //for encoder speedup
        if (bTestFirstBTSplit && rpcBestCU->isSkipped(0) && rpcBestCU->getBTDepth(0) == uiBTDepth && uiBTDepth >= SKIPHORNOVERQT_DEPTH_TH)
        {
          bTestSecondBTSplit = bQTSplit = false;
#if TT_ENC_OPT
          bTestHorSplitTT = bTestVerSplitTT = false;
#endif
        }

#if ENC_FAST_SKIP_RDC_EQUAL //action VBT
        if (bTestSecondBTSplit && !bSkipSecondBT)
#else
        if (bTestSecondBTSplit)
#endif
        {
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL && SIF_QTBT_FIX
#if JVET_C0024_QTBT
          UInt uiQTWidthIdx = g_aucConvertToBit[uiWidth >> 1];
          UInt uiQTHeightIdx = g_aucConvertToBit[uiHeight >> 1];
#else
          PartSize parentSize = rpcBestCU->getPartitionSize(0);
#endif
          int modeSIFnum = 1;
          if (rpcTempCU->isCodeBIF(uiWidth, /*uiDepth*/uiQTBTDepth
#if SCU_SIF
            , m_pcEntropyCoder->getEntropyCoder()->getSigDepth(0)
#endif        
            ))
          {
            modeSIFnum = 2;
          }
          for (int modeSIF = 0; modeSIF < modeSIFnum; modeSIF++)
          {
            if (modeSIFnum == 2)
            {
              rpcTempCU->setSIFToCheck(modeSIF);
              rpcBestCU->setSIFToCheck(modeSIF);
            }
#endif //RRC_SELECTABLE_IF
            // further split
            for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
            {
              //prepare parameters
              aiParam[0] = iSecondBTDir;
              aiParam[1] = iQP;
              aiParam[2] = lastCodedQP;
              aiParam[3] = bBTHorRmvEnable;
              aiParam[4] = bBTVerRmvEnable;

              CompressBT(rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight, uiBTSplitMode, aiParam, adParam, uiSplitConstrain);
#if TT_ENC_OPT
              dRDCost[uiSecondBTSplitIdx] = adParam[0];
              dRDDist[uiSecondBTSplitIdx] = adParam[1];
              iRDBits[uiSecondBTSplitIdx] = aiParam[5];
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
              if (iSecondBTDir == 0)
                dHorSplitCost = adParam[2];
              else
                dVerSplitCost = adParam[2];
#endif
            }
#if JVET_D0077_SAVE_LOAD_ENC_INFO
            if (dCostTempBest > adParam[2])
              dCostTempBest = adParam[2];
#endif
        }

#if TT_FAST_BAN_TT //ban tt from RDcost
          if (bTestFirstBTSplit || bTestSecondBTSplit)
          {
            if (pcSlice->getSliceType() != I_SLICE)
            {
              double dTH2 = (uiBTDepth <= 1) ? 0.95 : 0.92;
              bool bFullIsBest = true;
              if (bTestFirstBTSplit)
              {
                if (dRDCost[0] >= (dRDCost[uiFirstBTSplitIdx] * dTH2))
                  bFullIsBest = false;
              }
              if (bTestSecondBTSplit)
              {
                if (dRDCost[0] >= (dRDCost[uiSecondBTSplitIdx] * dTH2))
                  bFullIsBest = false;
              }
              if (!bTestFirstBTSplit && !bTestSecondBTSplit)
                bFullIsBest = false;

              if (bFullIsBest && uiBTDepth <= TT_FAST_BAN_TT_TH)
                bTestHorSplitTT = bTestVerSplitTT = bQTSplit = false;

#if TT_FAST_USE_ONE_TT //decide
              if (bTestFirstBTSplit && bTestSecondBTSplit)
              {
                if (uiBTDepth <= 1)
                {
                  if ((dRDCost[1] > dRDCost[2] * 1.125) || (dRDDist[1] > (dRDDist[2] * 1.03) && iRDBits[1] >= iRDBits[2] - 2))
                    uiSkipTTbyRD = 1;

                  if ((dRDCost[2] > dRDCost[1] * 1.125) || (dRDDist[2] > (dRDDist[1] * 1.03) && iRDBits[2] >= iRDBits[1] - 2))
                    uiSkipTTbyRD = 2;
                }
                else
                {
                  if ((dRDCost[1] > dRDCost[2] * 1.125) || (dRDDist[1] > (dRDDist[2] * 1.05) && iRDBits[1] >= iRDBits[2] - 2) && dRDCost[0] > dRDCost[2])
                    uiSkipTTbyRD = 1;

                  if ((dRDCost[2] > dRDCost[1] * 1.125) || (dRDDist[2] > (dRDDist[1] * 1.05) && iRDBits[2] >= iRDBits[1] - 2) && dRDCost[0] > dRDCost[1])
                    uiSkipTTbyRD = 2;
                }
                if (uiSkipTTbyRD == 1)
                  bTestHorSplitTT = false;
                if (uiSkipTTbyRD == 2)
                  bTestVerSplitTT = false;
              }
#endif
            }
          }
#endif

#if TT_ENC_OPT //xCompressCU, HORI_TT
          //for encoder speedup
          if ((bTestFirstBTSplit || bTestSecondBTSplit) && rpcBestCU->isSkipped(0) && rpcBestCU->getBTDepth(0) == uiBTDepth && uiBTDepth >= SKIPHORNOVERQT_DEPTH_TH)
            bTestHorSplitTT = bTestVerSplitTT = bQTSplit = false;

#if BT_FORCE_TT
          bTestHorSplitTT = bTestHorSplitTT || (
#if !BT_FORCE_INTRA_ON
            !rpcBestCU->getSlice()->isIntra() &&
#endif
            ((!rpcTempCU->isForcedBT(0, uiWidth, uiHeight)) && bBBoundary && (uiDepth >= iMaxQtLv)));
          bTestVerSplitTT = bTestVerSplitTT || (
#if !BT_FORCE_INTRA_ON
            !rpcBestCU->getSlice()->isIntra() &&
#endif
            ((!rpcTempCU->isForcedBT(0, uiWidth, uiHeight)) && bRBoundary && (uiDepth >= iMaxQtLv)));
#if (BT_FORCE_AR_TT>0)
          if (bBBoundary && bTestHorSplitTT)
          {
            if (uiWidth >= uiHeight*BT_FORCE_AR_TT)
            {
              bTestHorSplitTT = false;
              bTestVerSplitTT = true;
            }
          }
          if (bRBoundary && bTestVerSplitTT)
          {
            if (uiHeight >= uiWidth*BT_FORCE_AR_TT)
            {
              bTestVerSplitTT = false;
              bTestHorSplitTT = true;
            }
          }
#endif
#endif

          Bool bTestFirstTTSplit, bTestSecondTTSplit;
          UInt iFirstTTDir, iSecondTTDir;
          UInt uiFirstTTSplitIdx, uiSecondTTSplitIdx;
#if ADAPTIVE_RD_ORDER
          if (uiWidth < uiHeight) {
#endif
            bTestFirstTTSplit = bTestHorSplitTT;
            bTestSecondTTSplit = bTestVerSplitTT;
            iFirstTTDir = 0;
            iSecondTTDir = 1;
            uiFirstTTSplitIdx = 3;
            uiSecondTTSplitIdx = 4;
#if ADAPTIVE_RD_ORDER
  }else{
            bTestFirstTTSplit = bTestVerSplitTT;
            bTestSecondTTSplit = bTestHorSplitTT;
            iFirstTTDir = 1;
            iSecondTTDir = 0;
            uiFirstTTSplitIdx = 4;
            uiSecondTTSplitIdx = 3;
          }
#endif

          if (bTestFirstTTSplit)
          {
            // further split
            for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
            {
              //prepare parameters
              aiParam[0] = iFirstTTDir;
              aiParam[1] = iQP;
              aiParam[2] = lastCodedQP;

              CompressTT(rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight, uiBTSplitMode, aiParam, adParam, uiSplitConstrain);
#if TT_ENC_OPT
              dRDCost[uiFirstTTSplitIdx] = adParam[0];
              dRDDist[uiFirstTTSplitIdx] = adParam[1];
              iRDBits[uiFirstTTSplitIdx] = aiParam[5];
#endif
#if ENC_FAST_SKIP_RDC_EQUAL
              bSkipSecondTT = aiParam[6];
#endif
            }
          }
#endif

#if TT_ENC_OPT //xCompressCU, VERT_TT
          //for encoder speedup
          if ((bTestFirstBTSplit || bTestSecondBTSplit || bTestFirstTTSplit) && rpcBestCU->isSkipped(0) && rpcBestCU->getBTDepth(0) == uiBTDepth && uiBTDepth >= SKIPHORNOVERQT_DEPTH_TH)
            bTestSecondTTSplit = bQTSplit = false;

#if ENC_FAST_SKIP_RDC_EQUAL //action VTT
          if (bTestSecondTTSplit && !bSkipSecondTT)
#else
          if (bTestSecondTTSplit)
#endif
          {
            // further split
            for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
            {
              //prepare parameters
              aiParam[0] = iSecondTTDir;
              aiParam[1] = iQP;
              aiParam[2] = lastCodedQP;

              CompressTT(rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight, uiBTSplitMode, aiParam, adParam, uiSplitConstrain);
#if TT_ENC_OPT
              dRDCost[uiSecondTTSplitIdx] = adParam[0];
              dRDDist[uiSecondTTSplitIdx] = adParam[1];
              iRDBits[uiSecondTTSplitIdx] = aiParam[5];
#endif
            }
          }
#endif

#if !TT_ENC_OPT
          UInt uiZorderBR = g_auiRasterToZscan[((uiHeight >> MIN_CU_LOG2) - 1) * (sps.getCTUSize() >> MIN_CU_LOG2) + (uiWidth >> MIN_CU_LOG2) - 1];  //bottom-right part.
#endif

#if TT_ENC_OPT
          if (bQTSplit && ((rpcBestCU->getBTDepth(0) == 0 && uiMaxBTD >= (rpcBestCU->getSlice()->isIntra() ? 3 : 2))
            || (rpcBestCU->getBTDepth(0) == 1 && rpcBestCU->getBTDepth(uiZorderBR) == 1 && rpcBestCU->getBTDepth(uiZorderC1) == 1 && uiMaxBTD >= (rpcBestCU->getSlice()->isIntra() ? 4 : 3)))
            && (bTestFirstBTSplit || bTestFirstTTSplit) && (bTestSecondBTSplit || bTestSecondTTSplit))
          {
            bQTSplit = false;
          }
#else
          if (bQTSplit && ((rpcBestCU->getBTDepth(0) == 0 && uiMaxBTD >= (rpcBestCU->getSlice()->isIntra() ? 3 : 2))
            || (rpcBestCU->getBTDepth(0) == 1 && rpcBestCU->getBTDepth(uiZorderBR) == 1 && uiMaxBTD >= (rpcBestCU->getSlice()->isIntra() ? 4 : 3)))
            && bTestHorSplit && bTestVerSplit)
          {
            bQTSplit = false;
          }
#endif
#endif

#if JVET_C0024_QTBT
          if (bQTSplit)
#else
          if (bSubBranch && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() && (!getFastDeltaQp() || uiWidth > fastDeltaQPCuMaxSize || bBoundary))
#endif
          {
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL
#if JVET_C0024_QTBT
            UInt uiQTWidthIdx = g_aucConvertToBit[uiWidth >> 1];
            UInt uiQTHeightIdx = g_aucConvertToBit[uiHeight >> 1];
#else
            PartSize parentSize = rpcBestCU->getPartitionSize(0);
#endif
            int modeSIFnum = 1;
            if (rpcTempCU->isCodeBIF(uiWidth, uiDepth))
            {
              modeSIFnum = 2;
            }
            for (int modeSIF = 0; modeSIF < modeSIFnum; modeSIF++)
            {
              if (modeSIFnum == 2)
              {
                rpcTempCU->setSIFToCheck(modeSIF);
                rpcBestCU->setSIFToCheck(modeSIF);
              }
#endif
              // further split
              for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
              {
                const Bool bIsLosslessMode = false; // False at this level. Next level down may set it to true.

                rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

                UChar       uhNextDepth = uiDepth + 1;
#if JVET_C0024_QTBT
                UInt uiQTWidthIdx = g_aucConvertToBit[uiWidth >> 1];
                UInt uiQTHeightIdx = g_aucConvertToBit[uiHeight >> 1];
                TComDataCU* pcSubBestPartCU = m_pppcBestCU[uiQTWidthIdx][uiQTHeightIdx];
                TComDataCU* pcSubTempPartCU = m_pppcTempCU[uiQTWidthIdx][uiQTHeightIdx];
#else
                TComDataCU* pcSubBestPartCU = m_ppcBestCU[uhNextDepth];
                TComDataCU* pcSubTempPartCU = m_ppcTempCU[uhNextDepth];
#endif
                DEBUG_STRING_NEW(sTempDebug)
#if JVET_C0024_DELTA_QP_FIX // inherit quantization group info. from parent CU.
                  if (pps.getUseDQP())
                  {
                    pcSubBestPartCU->initSubCU(rpcTempCU, 0, uhNextDepth, iQP);           // clear sub partition datas or init.
                    pcSubTempPartCU->initSubCU(rpcTempCU, 0, uhNextDepth, iQP);           // clear sub partition datas or init.
                    pcSubBestPartCU->setCodedQP(lastCodedQP);
                    pcSubTempPartCU->setCodedQP(lastCodedQP);
                    pcSubBestPartCU->setQuPartIdx(rpcTempCU->getQuPartIdx());
                    pcSubTempPartCU->setQuPartIdx(rpcTempCU->getQuPartIdx());
                    pcSubBestPartCU->setQuLastCodedQP(rpcTempCU->getQuLastCodedQP());
                    pcSubTempPartCU->setQuLastCodedQP(rpcTempCU->getQuLastCodedQP());
                  }
#endif

                for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++)
                {
                  pcSubBestPartCU->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);           // clear sub partition datas or init.
                  pcSubTempPartCU->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);           // clear sub partition datas or init.
#if VCEG_AZ07_IMV && !JVET_C0024_QTBT
                  for (Int size = 0; size < NUMBER_OF_PART_SIZES; size++)
                  {
#if JVET_C0024_QTBT
                    m_pppcTempCUIMVCache[size][uiQTWidthIdx][uiQTHeightIdx]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);
#else
                    m_ppcTempCUIMVCache[size][uhNextDepth]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);
#endif
                  }
#endif
#if COM16_C806_OBMC
#if JVET_C0024_QTBT
                  m_pppcTempCUWoOBMC[uiQTWidthIdx][uiQTHeightIdx]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP); // clear sub partition datas or init.
#else
                  m_ppcTempCUWoOBMC[uhNextDepth]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP); // clear sub partition datas or init.
#endif
#endif
#if VCEG_AZ07_FRUC_MERGE
#if JVET_C0024_QTBT
                  m_pppcFRUCBufferCU[uiQTWidthIdx][uiQTHeightIdx]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);
#else
                  m_ppcFRUCBufferCU[uhNextDepth]->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP);
#endif
#endif
#if BARC_Boundary
                  UInt uiPelXInCTUtmp = pcSubBestPartCU->getCUPelX() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelX();
                  UInt uiPelYInCTUtmp = pcSubBestPartCU->getCUPelY() - pcSubBestPartCU->getPic()->getCtu(pcSubBestPartCU->getCtuRsAddr())->getCUPelY();

                  UInt deltaX = uiPelXInCTUtmp*isBARCFlag;
                  UInt deltaY = uiPelYInCTUtmp*isBARCFlag;

                  if(( pcSubBestPartCU->getCUPelX() + deltaX < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( pcSubBestPartCU->getCUPelY() + deltaY < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
#else
                  if ((pcSubBestPartCU->getCUPelX() < sps.getPicWidthInLumaSamples()) && (pcSubBestPartCU->getCUPelY() < sps.getPicHeightInLumaSamples()))
#endif
                  {
#if JVET_C0024_QTBT
                    if (0 == uiPartUnitIdx) //initialize RD with previous depth buffer
                    {
                      m_ppppcRDSbacCoder[uiQTWidthIdx][uiQTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_CURR_BEST]);
                    }
                    else
                    {
                      m_ppppcRDSbacCoder[uiQTWidthIdx][uiQTHeightIdx][CI_CURR_BEST]->load(m_ppppcRDSbacCoder[uiQTWidthIdx][uiQTHeightIdx][CI_NEXT_BEST]);
                    }
#else
                    if (0 == uiPartUnitIdx) //initialize RD with previous depth buffer
                    {
                      m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
                    }
                    else
                    {
                      m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
                    }
#endif

#if AMP_ENC_SPEEDUP
                    DEBUG_STRING_NEW(sChild)
#if JVET_C0024_QTBT
                      xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth, uiWidth >> 1, uiHeight >> 1 DEBUG_STRING_PASS_INTO(sChild), SIZE_2Nx2N);
#else
                      if (!(rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isInter(0)))
                      {
                        xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild), NUMBER_OF_PART_SIZES
                          );
                      }
                      else
                      {

                        xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild),
#if !RRC_SELECTABLE_IF || SIF_CU_LEVEL
                          rpcBestCU->getPartitionSize(0)
#else
                          parentSize
#endif 
                          );
                      }
#endif
                    DEBUG_STRING_APPEND(sTempDebug, sChild)
#else
                    xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth);
#endif

#if JVET_C0024_QTBT
                    rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidth >> 1, uiHeight >> 1);         // Keep best part data to current temporary data.
#if JVET_C0024_DELTA_QP_FIX
                    if (pps.getUseDQP())  // update coded QP
                    {
                      rpcTempCU->setCodedQP(pcSubBestPartCU->getCodedQP());
                      pcSubBestPartCU->setCodedQP(rpcTempCU->getCodedQP());
                      pcSubTempPartCU->setCodedQP(rpcTempCU->getCodedQP());
                    }
#endif
                    assert(pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx == pcSubBestPartCU->getZorderIdxInCtu() - rpcTempCU->getZorderIdxInCtu());
#if TT_TREE_DS1
                    xCopyYuv2Tmp(pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx, 0, uiWidth, uiHeight);
#else
                    xCopyYuv2Tmp(pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx, uiWidth, uiHeight);
#endif
#else
                    rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth);         // Keep best part data to current temporary data.
                    xCopyYuv2Tmp(pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx, uhNextDepth);
#endif
                  }
                  else
                  {
#if JVET_C0024_QTBT
                    pcSubBestPartCU->copyToPic(uhNextDepth, uiWidth >> 1, uiHeight >> 1);
                    rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth, uiWidth >> 1, uiHeight >> 1);

                    rpcBestCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 2);
#else
                    pcSubBestPartCU->copyToPic(uhNextDepth);
                    rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth);
#endif
                  }
                }

#if JVET_C0024_QTBT
                m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiQTWidthIdx][uiQTHeightIdx][CI_NEXT_BEST]);
#else
                m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
#endif
                if (!bBoundary)
                {
                  m_pcEntropyCoder->resetBits();
#if JVET_C0024_QTBT
                  if (!bForceQT)
#endif
                    m_pcEntropyCoder->encodeSplitFlag(rpcTempCU, 0, uiDepth, true);

                  rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
                  rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
                }
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL
                if (modeSIFnum == 2)
                {
                  m_pcEntropyCoder->resetBits();
                  m_pcEntropyCoder->encodeAdaptiveIfIndex(rpcTempCU, 0);
                  rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits();
                  rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
                }
#endif
                rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());
#if TT_ENC_OPT //store rd cost
                dRDCost[5] = rpcTempCU->getTotalCost();
                dRDDist[5] = rpcTempCU->getTotalDistortion();
                iRDBits[5] = rpcTempCU->getTotalBits();
#endif

#if JVET_C0024_DELTA_QP_FIX
                if (uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
                if (uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
                {
                  Bool hasResidual = false;
                  for (UInt uiBlkIdx = 0; uiBlkIdx < rpcTempCU->getTotalNumPart(); uiBlkIdx++)
                  {
#if JVET_C0024_DELTA_QP_FIX
                    if (rpcTempCU->getSlice()->isIntra())
                    {
                      if (rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA)
                      {
                        if (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Y))
                        {
                          hasResidual = true;
                        }
                      }
                      else
                      {
                        if ((rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cb) && (numberValidComponents > COMPONENT_Cb)) ||
                          (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cr) && (numberValidComponents > COMPONENT_Cr)))
                        {
                          hasResidual = true;
                        }
                      }
                    }
                    else
#endif
                      if ((rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Y)
                        || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cb) && (numberValidComponents > COMPONENT_Cb))
                        || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cr) && (numberValidComponents > COMPONENT_Cr))))
                      {
                        hasResidual = true;
                        break;
                      }
                  }

                  if (hasResidual)
                  {
#if JVET_C0024_DELTA_QP_FIX
                    Bool foundNonZeroCbf = false;
                    UInt uiFirstNonZeroPartIdx = 0;

                    rpcTempCU->setQPSubCUs(rpcTempCU->getRefQP(0), 0, uiDepth, uiWidth, uiHeight, uiFirstNonZeroPartIdx, foundNonZeroCbf);

                    m_pcEntropyCoder->resetBits();
                    m_pcEntropyCoder->encodeQP(rpcTempCU, uiFirstNonZeroPartIdx, false);
                    rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
                    rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
                    rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());
#else
                    m_pcEntropyCoder->resetBits();
                    m_pcEntropyCoder->encodeQP(rpcTempCU, 0, false);
                    rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
                    rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
                    rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

                    Bool foundNonZeroCbf = false;
                    rpcTempCU->setQPSubCUs(rpcTempCU->getRefQP(0), 0, uiDepth, foundNonZeroCbf);
#endif
                    assert(foundNonZeroCbf);
                  }
                  else
                  {
#if JVET_C0024_DELTA_QP_FIX
                    rpcTempCU->setQPSubParts(rpcTempCU->getRefQP(0), 0, uiWidth, uiHeight); // set QP to default QP
#else
                    rpcTempCU->setQPSubParts(rpcTempCU->getRefQP(0), 0, uiDepth); // set QP to default QP
#endif
#if JVET_C0024_DELTA_QP_FIX
                    if (pps.getUseDQP()) // update coded QP
                    {
                      rpcTempCU->setCodedQP(rpcTempCU->getQP(0));
                    }
#endif
                  }
                }

#if JVET_C0024_QTBT
                m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWidthIdx][uiHeightIdx][CI_TEMP_BEST]);
#else
                m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);
#endif

                // If the configuration being tested exceeds the maximum number of bytes for a slice / slice-segment, then
                // a proper RD evaluation cannot be performed. Therefore, termination of the
                // slice/slice-segment must be made prior to this CTU.
                // This can be achieved by forcing the decision to be that of the rpcTempCU.
                // The exception is each slice / slice-segment must have at least one CTU.
                if (rpcBestCU->getTotalCost() != MAX_DOUBLE)
                {
                  const Bool isEndOfSlice = pcSlice->getSliceMode() == FIXED_NUMBER_OF_BYTES
                    && ((pcSlice->getSliceBits() + rpcBestCU->getTotalBits()) > pcSlice->getSliceArgument() << 3)
                    && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceCurStartCtuTsAddr())
                    && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
                  const Bool isEndOfSliceSegment = pcSlice->getSliceSegmentMode() == FIXED_NUMBER_OF_BYTES
                    && ((pcSlice->getSliceSegmentBits() + rpcBestCU->getTotalBits()) > pcSlice->getSliceSegmentArgument() << 3)
                    && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
                  // Do not need to check slice condition for slice-segment since a slice-segment is a subset of a slice.
                  if (isEndOfSlice || isEndOfSliceSegment)
                  {
                    rpcBestCU->getTotalCost() = MAX_DOUBLE;
                  }
                }

#if JVET_C0024_QTBT
                xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth, uiWidth, uiHeight DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTempDebug) DEBUG_STRING_PASS_INTO(false)); // RD compare current larger prediction
                rpcBestCU->getPic()->setCodedBlkInCTU(false, uiPelXInCTU >> MIN_CU_LOG2, uiPelYInCTU >> MIN_CU_LOG2, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
                rpcBestCU->getPic()->addCodedAreaInCTU(-(Int)uiWidth*uiHeight);
#else
                xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTempDebug) DEBUG_STRING_PASS_INTO(false)); // RD compare current larger prediction
#endif
                                                                                                                                                       // with sub partitioned prediction.
              }
#if RRC_SELECTABLE_IF && !SIF_CU_LEVEL
            }
#endif
          }

          DEBUG_STRING_APPEND(sDebug_, sDebug);

#if JVET_C0024_QTBT
          rpcBestCU->copyToPic(uiDepth, uiWidth, uiHeight);                                                     // Copy Best data to Picture for next partition prediction.
          rpcBestCU->getPic()->setCodedBlkInCTU(true, uiPelXInCTU >> MIN_CU_LOG2, uiPelYInCTU >> MIN_CU_LOG2, uiWidth >> MIN_CU_LOG2, uiHeight >> MIN_CU_LOG2);
          rpcBestCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight);
          xCopyYuv2Pic( rpcBestCU->getPic(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu(), uiDepth, uiDepth, uiWidth, uiHeight
#if BARC
            ,isBARCFlag
#endif
    );   // Copy Yuv data to picture Yuv
#else
          rpcBestCU->copyToPic(uiDepth);                                                     // Copy Best data to Picture for next partition prediction.

          xCopyYuv2Pic(rpcBestCU->getPic(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu(), uiDepth, uiDepth);   // Copy Yuv data to picture Yuv
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO && !JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
          if (bUseSaveLoad)
          {
            m_pcPredSearch->setSaveLoadTag(MAX_UINT, uiWidthIdx - 1, uiHeightIdx - 1, SAVE_LOAD_INIT);
          }
          if (bUseSaveLoadSplitDecision && saveLoadTag == SAVE_ENC_INFO)
          {
            UChar c = 0;
            Double TH = JVET_D0077_SPLIT_DECISION_COST_SCALE * dCostTempBest;
            if (dNonSplitCost > TH)
            {
              c = c | 0x01;
            }
            if (dHorSplitCost > TH)
            {
              c = c | 0x02;
            }
            if (dVerSplitCost > TH)
            {
              c = c | 0x04;
            }
            m_pcPredSearch->setSaveLoadSplit(uiWidthIdx, uiHeightIdx, c);
          }
#endif

          if (bBoundary)
          {
            return;
          }

          // Assert if Best prediction mode is NONE
          // Selected mode's RD-cost must be not MAX_DOUBLE.
#if !JVET_C0024_QTBT
          assert(rpcBestCU->getPartitionSize(0) != NUMBER_OF_PART_SIZES);
          assert(rpcBestCU->getPredictionMode(0) != NUMBER_OF_PREDICTION_MODES);
          assert(rpcBestCU->getTotalCost() != MAX_DOUBLE);
#endif
}

/** finish encoding a cu and handle end-of-slice conditions
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TEncCu::finishCU( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
#if (TRACE_YUV || TRACE_FLAGS) && !JVET_C0024_QTBT
  if (gIsEncodingNow)
  {
#if TRACE_YUV
    MARK_PARTIONING(pcCU, uiAbsPartIdx, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx), "t_predicted.yuv");
    MARK_PARTIONING(pcCU, uiAbsPartIdx, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx), "t_rec_from_best.yuv");
#endif
#if TRACE_FLAGS
	DUMP_PU_FLAGS_INTER("encoder", pcCU, uiAbsPartIdx,-1);
#endif
  }
#endif
  TComPic* pcPic = pcCU->getPic();
  TComSlice * pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());

  //Calculate end address
  const Int  currentCTUTsAddr = pcPic->getPicSym()->getCtuRsToTsAddrMap(pcCU->getCtuRsAddr());
#if !JVET_C0024_QTBT
  const Bool isLastSubCUOfCtu = pcCU->isLastSubCUOfCtu(uiAbsPartIdx);
  if ( isLastSubCUOfCtu )
#endif
  {
    // The 1-terminating bit is added to all streams, so don't add it here when it's 1.
    // i.e. when the slice segment CurEnd CTU address is the current CTU address+1.
    if (pcSlice->getSliceSegmentCurEndCtuTsAddr() != currentCTUTsAddr+1)
    {
      m_pcEntropyCoder->encodeTerminatingBit( 0 );
    }
  }
}

/** Compute QP for each CU
 * \param pcCU Target CU
 * \param uiDepth CU depth
 * \returns quantization parameter
 */
Int TEncCu::xComputeQP( TComDataCU* pcCU, UInt uiDepth )
{
  Int iBaseQp = pcCU->getSlice()->getSliceQp();
  Int iQpOffset = 0;
  if ( m_pcEncCfg->getUseAdaptiveQP() )
  {
    TEncPic* pcEPic = dynamic_cast<TEncPic*>( pcCU->getPic() );
    UInt uiAQDepth = min( uiDepth, pcEPic->getMaxAQDepth()-1 );
    TEncPicQPAdaptationLayer* pcAQLayer = pcEPic->getAQLayer( uiAQDepth );
    UInt uiAQUPosX = pcCU->getCUPelX() / pcAQLayer->getAQPartWidth();
    UInt uiAQUPosY = pcCU->getCUPelY() / pcAQLayer->getAQPartHeight();
    UInt uiAQUStride = pcAQLayer->getAQPartStride();
    TEncQPAdaptationUnit* acAQU = pcAQLayer->getQPAdaptationUnit();

    Double dMaxQScale = pow(2.0, m_pcEncCfg->getQPAdaptationRange()/6.0);
    Double dAvgAct = pcAQLayer->getAvgActivity();
    Double dCUAct = acAQU[uiAQUPosY * uiAQUStride + uiAQUPosX].getActivity();
    Double dNormAct = (dMaxQScale*dCUAct + dAvgAct) / (dCUAct + dMaxQScale*dAvgAct);
    Double dQpOffset = log(dNormAct) / log(2.0) * 6.0;
    iQpOffset = Int(floor( dQpOffset + 0.49999 ));
  }

  return Clip3(-pcCU->getSlice()->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQp+iQpOffset );
}

/** encode a CU block recursively
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
#if JVET_C0024_QTBT
Void TEncCu::xEncodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt uiSplitConstrain )
#else
Void TEncCu::xEncodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
#endif
{
#if JVET_C0024_BT_RMV_REDUNDANT
  pcCU->setSplitConstrain( uiSplitConstrain );
  Bool bQTreeValid = false;
#endif

        TComPic   *const pcPic   = pcCU->getPic();
        TComSlice *const pcSlice = pcCU->getSlice();
  const TComSPS   &sps =*(pcSlice->getSPS());
  const TComPPS   &pps =*(pcSlice->getPPS());

#if !JVET_C0024_QTBT
  const UInt maxCUWidth  = sps.getMaxCUWidth();
  const UInt maxCUHeight = sps.getMaxCUHeight();

        Bool bBoundary = false;
#endif
        UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  const UInt uiRPelX   = uiLPelX + uiWidth  - 1;
#else
  const UInt uiRPelX   = uiLPelX + (maxCUWidth>>uiDepth)  - 1;
#endif
        UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if JVET_C0024_QTBT
  const UInt uiBPelY   = uiTPelY + uiHeight - 1;
#else
  const UInt uiBPelY   = uiTPelY + (maxCUHeight>>uiDepth) - 1;
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
	if (pcSlice->getSPS()->getBackgroundRefsPresent() && pcCU->getSlice()->getSliceType() != I_SLICE && uiDepth == 0 && uiAbsPartIdx == 0)
	{
		m_pcEntropyCoder->encodeLCUBgFlag(pcCU, 0);
	}
#endif

#if BT_FORCE
#if BARC_BTFORCE
  UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
  UInt deltaY = (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#else
  Bool bBBoundary = uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples();
  Bool bRBoundary = uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples();
#endif

#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(pcCU->getTextType()) ? pcSlice->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  const UInt uiDIMDMaxBtDepht = (pcCU->getBoundaryCtuLocation() && pcSlice->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (pcCU->getBoundaryCtuLocation() && pcSlice->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
#endif

#if JVET_C0024_QTBT
  Bool bForceQT = uiWidth > MAX_TU_SIZE;
  if( bForceQT )
  {
    assert(uiWidth == uiHeight);
  }
  UInt uiCTUSize = pcCU->getSlice()->getSPS()->getCTUSize();
#if JVET_C0024_DELTA_QP_FIX
  UInt uiQTWidth = uiCTUSize>>uiDepth;
  UInt uiQTHeight = uiCTUSize>>uiDepth;
  const UInt uiQTBTDepth = (uiDepth<<1) + (g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiHeight]);
  const UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif
  if (uiCTUSize>>uiDepth == uiWidth && uiWidth==uiHeight)
  {
#endif
#if BARC_Boundary
    UInt uiPelXInCTU = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
    UInt uiPelYInCTU = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

    UInt deltaX =  (uiPelXInCTU + uiWidth)*pcCU->getBARCFlag(pcCU->getTextType());
    UInt deltaY =  (uiPelYInCTU + uiHeight)*pcCU->getBARCFlag(pcCU->getTextType());

    if( ( uiRPelX < sps.getPicWidthInLumaSamples() - deltaX) && ( uiBPelY < sps.getPicHeightInLumaSamples() - deltaY) )
#else
  if( ( uiRPelX < sps.getPicWidthInLumaSamples() ) && ( uiBPelY < sps.getPicHeightInLumaSamples() ) )
#endif
  {
#if JVET_C0024_QTBT
    if( bForceQT )
    {
      assert(uiDepth < pcCU->getDepth( uiAbsPartIdx ) ); 
    }
    else
#endif
    m_pcEntropyCoder->encodeSplitFlag( pcCU, uiAbsPartIdx, uiDepth );
  }
  else
  {
#if JVET_C0024_QTBT
#if BT_FORCE
    if (
#if !BT_FORCE_INTRA_ON 
      pcSlice->isIntra() || 
#endif
      (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
      assert(uiDepth < pcCU->getDepth(uiAbsPartIdx));
#else
    bBoundary = true;
#endif
  }

#if JVET_C0024_BT_RMV_REDUNDANT
  bQTreeValid = true;
  UInt uiMinQTSize = sps.getMinQTSize(pcCU->getSlice()->getSliceType(), pcCU->getTextType());
  if ((uiCTUSize>>uiDepth) <= uiMinQTSize)
  {
    bQTreeValid = false;
  }
#endif

#if RRC_SELECTABLE_IF
#if !SIF_SPLIT_ONLY && !SIF_CU_LEVEL
  if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiDepth))
  {
    m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
  }
#endif
#endif

#if JVET_C0024_QTBT
  if ((uiDepth < pcCU->getDepth(uiAbsPartIdx))
#if BT_FORCE
    && (
#if !BT_FORCE_INTRA_ON
      pcSlice->isIntra() || 
#endif
      (!(bBBoundary || bRBoundary) || (uiDepth < iMaxQtLv)))
#endif
    )
#else
  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
#endif
  {
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL
    if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiDepth))
    {
      m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
    }
#endif
    UInt uiQNumParts = ( pcPic->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
#if JVET_C0024_DELTA_QP_FIX
    if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
    if( uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
    {
      setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setQuPartIdx( uiAbsPartIdx );
      pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
    }

    if( uiDepth == pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
    {
      setCodeChromaQpAdjFlag(true);
    }

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if BARC_Boundary
      UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
      UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

      UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
      if((  uiLPelX + deltaX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY + deltaY < sps.getPicHeightInLumaSamples() ) )
#else
      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
#endif
      {
#if JVET_C0024_QTBT
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth+1, uiWidth>>1, uiHeight>>1 );
#else
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth+1 );
#endif
      }
#if JVET_C0024_QTBT
      else
      {
        pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight>>2);
    }
#endif
    }
    return;
  }

#if JVET_C0024_QTBT
  }
#if JVET_C0024_BT_RMV_REDUNDANT
  Bool bBTHorRmvEnable = false;
  Bool bBTVerRmvEnable = false;
  if (pcCU->getSlice()->getSliceType() != I_SLICE)
  {
    bBTHorRmvEnable = true;
    bBTVerRmvEnable = bQTreeValid;
  }
#endif
#if JVET_C0024_SPS_MAX_BT_DEPTH
  UInt uiMaxBTD = pcSlice->isIntra() ? (isLuma(pcCU->getTextType())?sps.getMaxBTDepthISliceL():sps.getMaxBTDepthISliceC()): sps.getMaxBTDepth();
#if BT_FORCE_EXT
  Int iMaxELv =
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;

  if (
#if !BT_FORCE_INTRA_ON
    !pcSlice->isIntra() && 
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx) &&
#endif
    pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
#else
#if BT_FORCE_EXT
  Int iMaxELv = 
#if BT_FORCE_INTRA_ON
    pcSlice->isIntra() ? FORCE_BTD_EXT_INTRA :
#endif
    FORCE_BTD_EXT;
  if (
#if BT_FORCE_INTRA_ON
    !pcSlice->isIntra() &&
#endif
#if BT_FORCE_INTRA_SPEED_UP
    !pcCU->isFromForcedQT(uiAbsPartIdx) &&
#endif
     pcCU->getBoundaryCtuLocation() && iMaxELv)
  {
    UInt uiMaxBTDTmp = pcCU->getSlice()->getMaxFQTLevel(pcCU->getBoundaryCtuLocation() - 1) - iMaxQtLv + iMaxELv;
    uiMaxBTD = (uiMaxBTDTmp > uiMaxBTD) ? uiMaxBTDTmp : uiMaxBTD;
  }
#endif 
  UInt uiMaxBTD = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType())?MAX_BT_DEPTH:MAX_BT_DEPTH_C): MAX_BT_DEPTH_INTER;
#endif
#if JVET_C0024_SPS_MAX_BT_SIZE
  UInt uiMaxBTSize = pcSlice->isIntra() ? (isLuma(pcCU->getTextType())?sps.getMaxBTSizeISliceL():sps.getMaxBTSizeISliceC()): sps.getMaxBTSize();
#else
  UInt uiMaxBTSize = isLuma(pcCU->getTextType()) ? pcCU->getSlice()->getMaxBTSize(): MAX_BT_SIZE_C;
#endif
  UInt uiMinBTSize = pcCU->getSlice()->isIntra() ? (isLuma(pcCU->getTextType())?MIN_BT_SIZE:MIN_BT_SIZE_C): MIN_BT_SIZE_INTER;
  UInt uiBTDepth = pcCU->getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);

  if (((uiHeight > uiMinBTSize || uiWidth > uiMinBTSize)
    && uiWidth <= uiMaxBTSize && uiHeight <= uiMaxBTSize && uiBTDepth < uiMaxBTD)
#if BT_FORCE
    || (
#if !BT_FORCE_INTRA_ON
      !pcSlice->isIntra() && 
#endif
      ((bBBoundary || bRBoundary) && (uiDepth >= iMaxQtLv)))
#endif
    )
  {
#if JVET_C0024_BT_RMV_REDUNDANT
    uiSplitConstrain = 0;
#endif
#if BT_FORCE
#if !BT_FORCE_INTRA_ON
    if (pcSlice->isIntra())
    {
      m_pcEntropyCoder->encodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight);
    }
    else
#endif
      if (!(bBBoundary || bRBoundary))
      {
        m_pcEntropyCoder->encodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight);
      }
#else
    m_pcEntropyCoder->encodeBTSplitMode(pcCU, uiAbsPartIdx, uiWidth, uiHeight);
#endif
#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
    if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, uiQTBTDepth))
    {
      m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
    }
#endif
#if TT_TREE_DS1
    if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==ST_HORI_BT)
#else
    if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==1)
#endif
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyCoder->encodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyCoder->encodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx( uiAbsPartIdx );
        pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
      }
#endif
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
      if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, /*uiDepth*/uiQTBTDepth))
      {
        m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
      }
#endif
      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
      {
        if (uiPartUnitIdx==1)
        {
          uiAbsPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] 
          + (uiHeight>>1)/pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
        }
#if JVET_C0024_BT_RMV_REDUNDANT
#if BT_FORCE
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if BARC_BTFORCE
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xEncodeCU(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight >> 1, uiSplitConstrain);
#if TT_TREE_DS1
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == ST_VERT_BT && bBTHorRmvEnable && uiPartUnitIdx == 0)
#else
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == 2 && bBTHorRmvEnable && uiPartUnitIdx == 0)
#endif
          {
            uiSplitConstrain = 2;
          }
#if TT_RMV_PART1_BTBT && !TT_LEAF_NODE
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == ST_HORI_BT && uiPartUnitIdx == 0)
            uiSplitConstrain = 1;
#endif
#if BT_FORCE
        }
        else
        {
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
        }
#endif
#else
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight>>1 );
#endif
      }
      return;
    }
#if TT_TREE_DS1
    else if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==ST_VERT_BT)
#else
    else if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth)==2)
#endif
    {
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyCoder->encodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyCoder->encodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx( uiAbsPartIdx );
        pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
      }
#endif
#if RRC_SELECTABLE_IF && SIF_SPLIT_ONLY && !SIF_CU_LEVEL && SIF_QTBT_FIX
      if (pcCU->isCodeBIF(uiRPelX - uiLPelX + 1, /*uiDepth*/uiQTBTDepth))
      {
        m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
      }
#endif
      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
      {
        if (uiPartUnitIdx==1)
        {
          uiAbsPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] 
          + (uiWidth>>1)/pcCU->getPic()->getMinCUWidth()];
        }
#if JVET_C0024_BT_RMV_REDUNDANT
#if BT_FORCE
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
#if BARC_BTFORCE
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xEncodeCU(pcCU, uiAbsPartIdx, uiDepth, uiWidth >> 1, uiHeight, uiSplitConstrain);
#if TT_TREE_DS1
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == ST_HORI_BT && bBTVerRmvEnable && uiPartUnitIdx == 0)
#else
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == 1 && bBTVerRmvEnable && uiPartUnitIdx == 0)
#endif
          {
            uiSplitConstrain = 1;
          }
#if TT_RMV_PART1_BTBT && !TT_LEAF_NODE
          if (pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth + 1) == ST_VERT_BT && uiPartUnitIdx == 0)
            uiSplitConstrain = 2;
#endif
#if BT_FORCE
        }
        else
        {
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight >> 1);
        }
#endif
#else
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth, uiWidth>>1, uiHeight );
#endif
      }
      return;
    }

#if TT_ENC_OPT //xEncodeCU
    UInt uiHeightSubPart, uiWidthSubPart, uiAbstPartIdxSubPart;
    UInt uiSTSplitMode = pcCU->getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth);
    if (uiSTSplitMode==ST_HORI_TT_L || uiSTSplitMode==ST_HORI_TT_S )
    {
      assert(uiHeight>(uiMinBTSize<<1));
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyCoder->encodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyCoder->encodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx( uiAbsPartIdx );
        pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
      }
#endif

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
      {
        if (uiPartUnitIdx==1)
        {
          uiAbstPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + (uiHeight>>2)/pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
          uiHeightSubPart = uiHeight>>1;
        }
        else
        {
          if(uiPartUnitIdx==2)
            uiAbstPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + ((uiHeight>>2)+(uiHeight>>1))/pcCU->getPic()->getMinCUHeight()*pcCU->getPic()->getNumPartInCtuWidth()];
          else
            uiAbstPartIdxSubPart = uiAbsPartIdx;
          uiHeightSubPart = uiHeight>>2;
        }
#if BT_FORCE_TT
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbstPartIdxSubPart]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbstPartIdxSubPart]];
#if BARC_Boundary
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xEncodeCU(pcCU, uiAbstPartIdxSubPart, uiDepth, uiWidth, uiHeightSubPart, 0);
#if BT_FORCE_TT
        }
        else
        {
          pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeightSubPart);
        }
#endif
      }
      return;
    }

    if(uiSTSplitMode == ST_VERT_TT_L || uiSTSplitMode == ST_VERT_TT_S )
    {
      assert(uiWidth>(uiMinBTSize<<1));
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == uiDIMDMaxBtDepht)
#else
      if (isLuma(pcCU->getTextType()) && pcCU->getSlice()->isIntra() && uiBTDepth == DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_pcEntropyCoder->encodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        if (pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
        {
          m_pcEntropyCoder->encodeDIMDNoBTFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
        }
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( uiQTBTDepth == uiMaxDQPDepthQTBT && pps.getUseDQP())
      {
        setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
        pcCU->setQuPartIdx( uiAbsPartIdx );
        pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
      }
#endif

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
      {
        if (uiPartUnitIdx==1)
        {
          uiAbstPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + (uiWidth>>2)/pcCU->getPic()->getMinCUWidth()];
          uiWidthSubPart = uiWidth>>1;
        }
        else
        {
          if(uiPartUnitIdx==2)
            uiAbstPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx] + ((uiWidth>>2)+(uiWidth>>1))/pcCU->getPic()->getMinCUWidth()];
          else
            uiAbstPartIdxSubPart = uiAbsPartIdx;
          uiWidthSubPart = uiWidth>>2;
        }
#if BT_FORCE_TT
        uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbstPartIdxSubPart]];
        uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbstPartIdxSubPart]];
#if BARC_Boundary
        UInt uiPelXInCTUTmp = uiLPelX - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelX();
        UInt uiPelYInCTUTmp = uiTPelY - pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getCUPelY();

        UInt deltaX = uiPelXInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        UInt deltaY = uiPelYInCTUTmp*pcCU->getBARCFlag(pcCU->getTextType());
        if ((uiLPelX + deltaX < sps.getPicWidthInLumaSamples()) && (uiTPelY + deltaY < sps.getPicHeightInLumaSamples()))
#else
        if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
#endif
        {
#endif
          xEncodeCU(pcCU, uiAbstPartIdxSubPart, uiDepth, uiWidthSubPart, uiHeight, 0);
#if BT_FORCE_TT
        }
        else
        {
          pcCU->getPic()->addCodedAreaInCTU(uiWidthSubPart*uiHeight);
        }
#endif
      }
      return;
    }
#endif
  }

  pcCU->getPic()->addCodedAreaInCTU(uiWidth*uiHeight);
  UInt uiBlkX = g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ] >> MIN_CU_LOG2;
  UInt uiBlkY = g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ] >> MIN_CU_LOG2;
  pcCU->getPic()->setCodedBlkInCTU(true, uiBlkX, uiBlkY, uiWidth>> MIN_CU_LOG2, uiHeight>> MIN_CU_LOG2);

#endif

#if JVET_C0024_AMAX_BT
  if (!pcCU->getSlice()->isIntra())
  {
    g_uiBlkSize[pcCU->getSlice()->getDepth()] += uiWidth*uiHeight;
    g_uiNumBlk[pcCU->getSlice()->getDepth()]++;
  }
#endif
#if JVET_C0024_DELTA_QP_FIX
  if( uiQTBTDepth <= uiMaxDQPDepthQTBT && pps.getUseDQP())
#else
  if( uiDepth <= pps.getMaxCuDQPDepth() && pps.getUseDQP())
#endif
  {
    setdQPFlag(true);
#if JVET_C0024_DELTA_QP_FIX
    pcCU->setQuPartIdx( uiAbsPartIdx );
    pcCU->setQuLastCodedQP( pcCU->getCodedQP() );
#endif
  }

  if( uiDepth <= pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
  {
    setCodeChromaQpAdjFlag(true);
  }

  if (pps.getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( pcCU, uiAbsPartIdx );
  }

  if( !pcSlice->isIntra() )
  {
    m_pcEntropyCoder->encodeSkipFlag( pcCU, uiAbsPartIdx );
  }

  if( pcCU->isSkipped( uiAbsPartIdx ) )
  {
#if VCEG_AZ07_FRUC_MERGE
    m_pcEntropyCoder->encodeFRUCMgrMode( pcCU , uiAbsPartIdx , 0 );
    if( !pcCU->getFRUCMgrMode( uiAbsPartIdx ) )
#endif

#if MRG_CONSIDER_COMPLEX
    {
      if ((pcCU->getWidth(uiAbsPartIdx) >= 16) || (pcCU->getHeight(uiAbsPartIdx) >= 16))
      {
        m_pcEntropyCoder->encodeMergeComplexFlag(pcCU, uiAbsPartIdx);
        if (pcCU->getMergeComplexFlag(uiAbsPartIdx))
        {
          m_pcEntropyCoder->encodeMergeComplexIndex(pcCU, uiAbsPartIdx);
        }
      }

      if (!pcCU->getMergeComplexFlag(uiAbsPartIdx))
#endif
#if COM16_C1016_AFFINE
    {
      if ( pcCU->isAffineMrgFlagCoded(uiAbsPartIdx, 0) )
      {
        m_pcEntropyCoder->encodeAffineFlag( pcCU, uiAbsPartIdx, 0 );
      }

      if ( !pcCU->isAffine(uiAbsPartIdx) )
      {
        m_pcEntropyCoder->encodeMergeIndex( pcCU, uiAbsPartIdx );
      }
    }
#else
    m_pcEntropyCoder->encodeMergeIndex( pcCU, uiAbsPartIdx );
#endif
#if MRG_CONSIDER_COMPLEX
  }
#endif
#if VCEG_AZ06_IC
    m_pcEntropyCoder->encodeICFlag  ( pcCU, uiAbsPartIdx );
#endif
#if SIF_CU_LEVEL
    m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
#endif
#if !JVET_C0024_QTBT
    finishCU(pcCU,uiAbsPartIdx);
#endif
#if JVET_C0024_DELTA_QP_FIX
    if( pps.getUseDQP() )
    { 
      pcCU->setCodedQP( pcCU->getQP(uiAbsPartIdx) );
    }
#endif
#if TRACE_YUV || TRACE_FLAGS
    if (gIsEncodingNow)
    {
#if JVET_C0024_QTBT && TRACE_YUV
      MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_predicted.yuv");
      MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_rec_from_best.yuv");
#endif
#if  TRACE_FLAGS
#if !JVET_C0024_QTBT
	  if (pcCU->isIntra(uiAbsPartIdx))
#endif
	  {
		  DUMP_CU_FLAGS_INFO("encoder", pcCU, uiAbsPartIdx, uiAbsPartIdx);
	  }
#endif
    }
#endif
    return;
  }

  m_pcEntropyCoder->encodePredMode( pcCU, uiAbsPartIdx );
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
  if (isLuma(pcCU->getTextType()) && (uiBTDepth <= uiDIMDMaxBtDepht || !pcCU->getSlice()->isIntra()))
#else
  if (isLuma(pcCU->getTextType()) && (uiBTDepth <= DIMD_MAX_BT_DEPTH_SIG_AI || !pcCU->getSlice()->isIntra()))
#endif
  {
    m_pcEntropyCoder->encodeDIMDFlag(pcCU, uiAbsPartIdx, uiDepth, uiWidth, uiHeight);
  }
#endif
#if JVET_C0024_QTBT
  if (isLuma(pcCU->getTextType()))
  {

#endif
#if VCEG_AZ05_INTRA_MPI
  m_pcEntropyCoder->encodeMPIIdx(pcCU, uiAbsPartIdx);
#endif   
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  m_pcEntropyCoder->encodePDPCIdx(pcCU, uiAbsPartIdx);
#endif  
#if JVET_C0024_QTBT
  }
#else
  m_pcEntropyCoder->encodePartSize( pcCU, uiAbsPartIdx, uiDepth );
#endif

#if JVET_C0024_QTBT
  if (pcCU->isIntra( uiAbsPartIdx ) )
#else
  if (pcCU->isIntra( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
  {
    m_pcEntropyCoder->encodeIPCMInfo( pcCU, uiAbsPartIdx );

    if(pcCU->getIPCMFlag(uiAbsPartIdx))
    {
#if !JVET_C0024_QTBT
      // Encode slice finish
      finishCU(pcCU,uiAbsPartIdx);
#endif
#if JVET_C0024_DELTA_QP_FIX
      if( pps.getUseDQP() )
      { 
        pcCU->setCodedQP( pcCU->getQP(uiAbsPartIdx) );
      }
#endif
#if TRACE_YUV || TRACE_FLAGS
      if (gIsEncodingNow)
      {
#if JVET_C0024_QTBT && TRACE_YUV
        MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_predicted.yuv");
        MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_rec_from_best.yuv");
#endif
#if TRACE_FLAGS
		DUMP_CU_FLAGS_INFO("encoder", pcCU, uiAbsPartIdx, uiAbsPartIdx);
#endif
      }
#endif
      return;
    }
  }

  // prediction Info ( Intra : direction mode, Inter : Mv, reference idx )
  m_pcEntropyCoder->encodePredInfo( pcCU, uiAbsPartIdx );
#if COM16_C806_OBMC
  m_pcEntropyCoder->encodeOBMCFlag( pcCU, uiAbsPartIdx );
#endif
#if SIF_CU_LEVEL
  m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, uiAbsPartIdx);
#endif
#if VCEG_AZ06_IC
  m_pcEntropyCoder->encodeICFlag  ( pcCU, uiAbsPartIdx );
#endif
  // Encode Coefficients
  Bool bCodeDQP = getdQPFlag();
  Bool codeChromaQpAdj = getCodeChromaQpAdjFlag();
#if  VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
  Int bNonZeroCoeff = false;
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  Int iNonZeroCoeffNonTs;
#endif
  m_pcEntropyCoder->encodeCoeff( pcCU, uiAbsPartIdx, uiDepth, bCodeDQP, codeChromaQpAdj
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , bNonZeroCoeff
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
  setCodeChromaQpAdjFlag( codeChromaQpAdj );
  setdQPFlag( bCodeDQP );
#if JVET_C0024_DELTA_QP_FIX
  if( pps.getUseDQP() )
  { 
    pcCU->setCodedQP( pcCU->getQP(uiAbsPartIdx) );
  }
#endif

#if !JVET_C0024_QTBT
  // --- write terminating bit ---
  finishCU(pcCU,uiAbsPartIdx);
#endif
#if (TRACE_FLAGS || TRACE_YUV) && JVET_C0024_QTBT
  if (gIsEncodingNow)
  {
  
#if TRACE_YUV
    MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_predicted.yuv");
    MARK_PARTIONING(pcCU, uiAbsPartIdx, uiWidth, uiHeight, "t_rec_from_best.yuv");
#endif
#if  TRACE_FLAGS
	DUMP_CU_FLAGS_INFO("encoder", pcCU, uiAbsPartIdx, uiAbsPartIdx);
#endif
  }
#endif
}

Int xCalcHADs8x8_ISlice(Pel *piOrg, Int iStrideOrg)
{
  Int k, i, j, jj;
  Int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

  for( k = 0; k < 64; k += 8 )
  {
    diff[k+0] = piOrg[0] ;
    diff[k+1] = piOrg[1] ;
    diff[k+2] = piOrg[2] ;
    diff[k+3] = piOrg[3] ;
    diff[k+4] = piOrg[4] ;
    diff[k+5] = piOrg[5] ;
    diff[k+6] = piOrg[6] ;
    diff[k+7] = piOrg[7] ;

    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 3;
    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i=0; i < 8; i++)
  {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      iSumHad += abs(m2[i][j]);
    }
  }
  iSumHad -= abs(m2[0][0]);
  iSumHad =(iSumHad+2)>>2;
  return(iSumHad);
}

Int  TEncCu::updateCtuDataISlice(TComDataCU* pCtu, Int width, Int height)
{
  Int  xBl, yBl;
  const Int iBlkSize = 8;

  Pel* pOrgInit   = pCtu->getPic()->getPicYuvOrg()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);
  Int  iStrideOrig = pCtu->getPic()->getPicYuvOrg()->getStride(COMPONENT_Y);
  Pel  *pOrg;

  Int iSumHad = 0;
  for ( yBl=0; (yBl+iBlkSize)<=height; yBl+= iBlkSize)
  {
    for ( xBl=0; (xBl+iBlkSize)<=width; xBl+= iBlkSize)
    {
      pOrg = pOrgInit + iStrideOrig*yBl + xBl;
      iSumHad += xCalcHADs8x8_ISlice(pOrg, iStrideOrig);
    }
  }
  return(iSumHad);
}

/** check RD costs for a CU block encoded with merge
 * \param rpcBestCU
 * \param rpcTempCU
 * \param earlyDetectionSkipMode
 */
Void TEncCu::xCheckRDCostMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU DEBUG_STRING_FN_DECLARE(sDebug), Bool *earlyDetectionSkipMode )
{
#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[rpcBestCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcBestCU->getHeight(0)];
#else
#if COM16_C806_LARGE_CTU
  if( m_pcEncCfg->getUseFastLCTU() && rpcTempCU->getHeight( 0 ) * 2 > rpcTempCU->getSlice()->getSPS()->getPicHeightInLumaSamples() )
  {
    rpcTempCU->getTotalCost() = MAX_DOUBLE / 4;
    rpcTempCU->getTotalDistortion() = MAX_INT;
    xCheckBestMode(rpcBestCU, rpcTempCU, rpcTempCU->getDepth( 0 ));
    return;
  }
#endif
#endif

  assert( rpcTempCU->getSlice()->getSliceType() != I_SLICE );
  if(getFastDeltaQp())
  {
    return;   // never check merge in fast deltaqp mode
  }
  TComMvField  cMvFieldNeighbours[2 * MRG_MAX_NUM_CANDS]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0;
  const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UChar  eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
  memset ( eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, sizeof(UChar)*MRG_MAX_NUM_CANDS );
#endif
#if VCEG_AZ06_IC
  Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif

  for( UInt ui = 0; ui < rpcTempCU->getSlice()->getMaxNumMergeCand(); ++ui )
  {
    uhInterDirNeighbours[ui] = 0;
  }
  UChar uhDepth = rpcTempCU->getDepth( 0 );
#if !JVET_C0024_QTBT
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to CTU level
#endif


#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if !JVET_C0024_QTBT  
  for (Int i=0 , i2 = 0 ; i< rpcTempCU->getTotalNumPart(); i++ , i2 += 2)
  {
#if JVET_C0035_ATMVP_SIMPLIFICATION
    for (Int j=0; j < NUM_MGR_TYPE ; j++)
    {
      m_phInterDirSP[j][i] = 0;
      m_pMvFieldSP[j][i2].setRefIdx(-1);
      m_pMvFieldSP[j][i2+1].setRefIdx(-1);
    }
#else
    m_phInterDirSP[1][i] = 0;
    m_pMvFieldSP[1][i2].setRefIdx(-1);
    m_pMvFieldSP[1][i2+1].setRefIdx(-1);
#endif
  }
#endif
#endif

  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,uhInterDirNeighbours, numValidMergeCand
#if VCEG_AZ06_IC
    , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    , eMergeCandTypeNieghors
    , m_pMvFieldSP
    , m_phInterDirSP
#endif
    );

  Int mergeCandBuffer[MRG_MAX_NUM_CANDS];
  for( UInt ui = 0; ui < numValidMergeCand; ++ui )
  {
    mergeCandBuffer[ui] = 0;
  }

#if JVET_C0024_FAST_MRG
  Bool bestIsSkip = rpcBestCU->getPic()->getSkiped(rpcBestCU->getZorderIdxInCtu(), rpcBestCU->getWidth(0), rpcBestCU->getHeight(0));
  UInt uiNumMrgSATDCand = NUM_MRG_SATD_CAND;
  UInt uiRdModeList[MRG_MAX_NUM_CANDS];
  Double CandCostList[MRG_MAX_NUM_CANDS];
  for (UInt i=0; i<MRG_MAX_NUM_CANDS; i++)
  {
    uiRdModeList[i] = i;
    CandCostList[i] = MAX_DOUBLE;
  }

  Bool bMrgTempBufSet = false;
  if (!bestIsSkip)
  {
#if JVET_D0123_ME_CTX_LUT_BITS
    UInt uiMrgIdxBits[MRG_MAX_NUM_CANDS];
    m_pcPredSearch->getMrgCandBits(rpcBestCU, 0, uiMrgIdxBits);
#endif
    bMrgTempBufSet = true;
    for( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
    {
      // set MC parameters
      rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth ); // interprets depth relative to CTU level
      rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag, 0, uhDepth );
      rpcTempCU->setChromaQpAdjSubParts( bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );

      rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth ); // interprets depth relative to CTU level
      rpcTempCU->setMergeIndexSubParts( uiMergeCand, 0, 0, uhDepth ); // interprets depth relative to CTU level
#if COM16_C806_OBMC
      rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
      rpcTempCU->setAffineFlagSubParts(false, 0, 0, uhDepth);
#endif
#if MRG_CONSIDER_COMPLEX
      rpcTempCU->setMergeComplexFlagSubParts(false, 0);
#endif
#if VCEG_AZ06_IC
      rpcTempCU->setICFlagSubParts( rpcTempCU->getSlice()->getApplyIC() ? abICFlag[uiMergeCand] : 0, 0, uhDepth );
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
      rpcTempCU->setMergeTypeSubParts(eMergeCandTypeNieghors[uiMergeCand] , 0, 0, uhDepth ); 
      if( eMergeCandTypeNieghors[uiMergeCand])
      {
        UInt uiSPAddr;
        Int iWidth = rpcTempCU->getWidth(0);
        Int iHeight = rpcTempCU->getHeight(0);

        Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
        UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand];
#else
        UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand] == MGR_TYPE_SUBPU_TMVP?0:1;
#endif
        rpcTempCU->getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

        for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
        {
          rpcTempCU->getSPAbsPartIdx(0, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
          rpcTempCU->setInterDirSP(m_phInterDirSP[uiSPListIndex][iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
          rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(rpcTempCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx], iSPWidth, iSPHeight);
          rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(rpcTempCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx + 1], iSPWidth, iSPHeight);
        }
      }
      else
      {
#endif
        rpcTempCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth ); // interprets depth relative to CTU level
        rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMvFieldNeighbours[0 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
        rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMvFieldNeighbours[1 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
      }
#endif

      m_pcMrgPredTempYuv[uiMergeCand]->setWidth(rpcBestCU->getWidth(0));
      m_pcMrgPredTempYuv[uiMergeCand]->setHeight(rpcBestCU->getHeight(0));

      m_pcPredSearch->motionCompensation ( rpcTempCU, m_pcMrgPredTempYuv[uiMergeCand] );
#if COM16_C806_OBMC
      m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pcMrgPredTempYuv[uiMergeCand], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
#if JVET_E0052_DMVR
// the mv may be changed in the mv refinement during the MC
  if (uhInterDirNeighbours[uiMergeCand]==3 && eMergeCandTypeNieghors[uiMergeCand]== MGR_TYPE_DEFAULT_N)
  {
      cMvFieldNeighbours[0 + 2*uiMergeCand].getMv() = rpcTempCU->getCUMvField(REF_PIC_LIST_0)->getMv(0);
      cMvFieldNeighbours[1 + 2*uiMergeCand].getMv() = rpcTempCU->getCUMvField(REF_PIC_LIST_1)->getMv(0);
  }
  else if (eMergeCandTypeNieghors[uiMergeCand] != MGR_TYPE_DEFAULT_N)
  {
      UInt uiSPAddr;
      Int iWidth = rpcTempCU->getWidth(0);
      Int iHeight = rpcTempCU->getHeight(0);

      Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
      UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand];
#else
      UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand] == MGR_TYPE_SUBPU_TMVP ? 0 : 1;
#endif
      rpcTempCU->getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

    for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
    {
      if (m_phInterDirSP[uiSPListIndex][iPartitionIdx] == 3)
      {
        rpcTempCU->getSPAbsPartIdx(0, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
        m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx].getMv() = rpcTempCU->getCUMvField(REF_PIC_LIST_0)->getMv(uiSPAddr);
        m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx + 1].getMv() = rpcTempCU->getCUMvField(REF_PIC_LIST_1)->getMv(uiSPAddr);
      }
    }
  }   
#endif
      // use hadamard transform here
      DistParam distParam;
      const Bool bUseHadamard=rpcTempCU->getCUTransquantBypass(0) == 0;
      m_pcRdCost->setDistParam(distParam, rpcTempCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), m_pppcOrigYuv[uiWIdx][uiHIdx]->getAddr(COMPONENT_Y)
        , m_pppcOrigYuv[uiWIdx][uiHIdx]->getStride(COMPONENT_Y)
        , m_pcMrgPredTempYuv[uiMergeCand]->getAddr(COMPONENT_Y), m_pcMrgPredTempYuv[uiMergeCand]->getStride(COMPONENT_Y)
        , rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), bUseHadamard);
      distParam.bApplyWeight = false;

      UInt uiSad = distParam.DistFunc(&distParam);
#if JVET_D0123_ME_CTX_LUT_BITS
      UInt uiBitsCand = uiMrgIdxBits[uiMergeCand];
      Double cost  = (Double)uiSad + (Double)uiBitsCand/((Double)EPBIT) * m_pcRdCost->getSqrtLambda();
#else
      UInt uiBitsCand = uiMergeCand + 1;                                         
      if (uiMergeCand == rpcTempCU->getSlice()->getMaxNumMergeCand() -1)
      {
        uiBitsCand--;
      }   
      Double cost      = (Double)uiSad + (Double)uiBitsCand * m_pcRdCost->getSqrtLambda();
#endif
      TEncSearch::updateCandList( uiMergeCand, cost, uiNumMrgSATDCand, uiRdModeList, CandCostList );
    }
    for (UInt i=1; i<uiNumMrgSATDCand; i++)
    {
      if (CandCostList[i] > MRG_FAST_RATIO*CandCostList[0])
      {
        uiNumMrgSATDCand = i;
        break;
      }
    }
  }
#else
  Bool bestIsSkip = false;
#endif

  UInt iteration;
  if ( rpcTempCU->isLosslessCoded(0))
  {
    iteration = 1;
  }
  else
  {
    iteration = 2;
  }
  DEBUG_STRING_NEW(bestStr)

  for( UInt uiNoResidual = 0; uiNoResidual < iteration; ++uiNoResidual )
  {
#if JVET_C0024_FAST_MRG
      for( UInt uiMrgHADIdx = 0; uiMrgHADIdx < uiNumMrgSATDCand; ++uiMrgHADIdx )
      {
        UInt uiMergeCand = uiRdModeList[uiMrgHADIdx];
#else
    for( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
    {
#endif
      if(!(uiNoResidual==1 && mergeCandBuffer[uiMergeCand]==1))
      {
        if( !(bestIsSkip && uiNoResidual == 0) )
        {
          DEBUG_STRING_NEW(tmpStr)
          // set MC parameters
          rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag, 0, uhDepth );
          rpcTempCU->setChromaQpAdjSubParts( bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );
#if !JVET_C0024_QTBT
          rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to CTU level
#endif
          rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setMergeIndexSubParts( uiMergeCand, 0, 0, uhDepth ); // interprets depth relative to CTU level

#if COM16_C806_OBMC
          rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
#endif
#if MRG_CONSIDER_COMPLEX
          rpcTempCU->setMergeComplexFlagSubParts(false, 0);
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
          rpcTempCU->setAffineFlagSubParts(false, 0, 0, uhDepth);
#endif
#if VCEG_AZ06_IC
          rpcTempCU->setICFlagSubParts( rpcTempCU->getSlice()->getApplyIC() ? abICFlag[uiMergeCand] : 0, 0, uhDepth );
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          rpcTempCU->setMergeTypeSubParts(eMergeCandTypeNieghors[uiMergeCand] , 0, 0, uhDepth ); 
          if( eMergeCandTypeNieghors[uiMergeCand])
          {
            UInt uiSPAddr;
            Int iWidth = rpcTempCU->getWidth(0);
            Int iHeight = rpcTempCU->getHeight(0);

            Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
            UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand];
#else
            UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeCand] == MGR_TYPE_SUBPU_TMVP?0:1;
#endif
            rpcTempCU->getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

            for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
            {
              rpcTempCU->getSPAbsPartIdx(0, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
              rpcTempCU->setInterDirSP(m_phInterDirSP[uiSPListIndex][iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
              rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(rpcTempCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx], iSPWidth, iSPHeight);
              rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(rpcTempCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2*iPartitionIdx + 1], iSPWidth, iSPHeight);
            }
          }
          else
          {
#endif
          rpcTempCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMvFieldNeighbours[0 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMvFieldNeighbours[1 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          }
#endif
          // do MC
#if JVET_C0024_QTBT
#if JVET_C0024_FAST_MRG
            if (bMrgTempBufSet)
            {
              m_pcMrgPredTempYuv[uiMergeCand]->copyToPartYuv(m_pppcPredYuvTemp[uiWIdx][uiHIdx], 0);
            }
            else
            {
              m_pcPredSearch->motionCompensation ( rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx] );
#if COM16_C806_OBMC
              m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
            }
#else
            m_pcPredSearch->motionCompensation ( rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx] );
#if COM16_C806_OBMC
            m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
#endif
            // estimate residual and encode everything
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
              m_pppcOrigYuv    [uiWIdx][uiHIdx],
              m_pppcPredYuvTemp[uiWIdx][uiHIdx],
              m_pppcResiYuvTemp[uiWIdx][uiHIdx],
              m_pppcResiYuvBest[uiWIdx][uiHIdx],
              m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
              (uiNoResidual != 0) 
#if COM16_C806_EMT
              , rpcBestCU->getTotalCost()
#endif
              DEBUG_STRING_PASS_INTO(tmpStr) );

#if DEBUG_STRING
            DebugInterPredResiReco(tmpStr, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif
#else
          m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[uhDepth] );
#if COM16_C806_OBMC
          m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth] );
#endif
          // estimate residual and encode everything
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                     m_ppcOrigYuv    [uhDepth],
                                                     m_ppcPredYuvTemp[uhDepth],
                                                     m_ppcResiYuvTemp[uhDepth],
                                                     m_ppcResiYuvBest[uhDepth],
                                                     m_ppcRecoYuvTemp[uhDepth],
                                                     (uiNoResidual != 0) 
#if COM16_C806_EMT
                                                     , rpcBestCU->getTotalCost()
#endif
                                                     DEBUG_STRING_PASS_INTO(tmpStr) );

#if DEBUG_STRING
          DebugInterPredResiReco(tmpStr, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif
#endif
          if ((uiNoResidual == 0) && (rpcTempCU->getQtRootCbf(0) == 0))
          {
            // If no residual when allowing for one, then set mark to not try case where residual is forced to 0
            mergeCandBuffer[uiMergeCand] = 1;
          }

          Int orgQP = rpcTempCU->getQP( 0 );
          xCheckDQP( rpcTempCU );
#if BIF_NO_SKIP_ONLY
          if (rpcTempCU->getSIFFlag(0) && rpcTempCU->isSkipped(0))
          {
            rpcTempCU->getTotalCost() = DBL_MAX;
          }
#endif
          xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr));

          rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag );

          if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
          {
            bestIsSkip = rpcBestCU->getQtRootCbf(0) == 0;
          }
        }
      }
    }

    if(uiNoResidual == 0 && m_pcEncCfg->getUseEarlySkipDetection())
    {
      if(rpcBestCU->getQtRootCbf( 0 ) == 0)
      {
        if( rpcBestCU->getMergeFlag( 0 ))
        {
          *earlyDetectionSkipMode = true;
        }
        else if(m_pcEncCfg->getFastSearch() != SELECTIVE)
        {
          Int absoulte_MV=0;
          for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if ( rpcBestCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              TComCUMvField* pcCUMvField = rpcBestCU->getCUMvField(RefPicList( uiRefListIdx ));
              Int iHor = pcCUMvField->getMvd( 0 ).getAbsHor();
              Int iVer = pcCUMvField->getMvd( 0 ).getAbsVer();
              absoulte_MV+=iHor+iVer;
            }
          }

          if(absoulte_MV == 0)
          {
            *earlyDetectionSkipMode = true;
          }
        }
      }
    }
  }
  DEBUG_STRING_APPEND(sDebug, bestStr)

#if MRG_CONSIDER_COMPLEX
    if (rpcTempCU->getSlice()->getSPS()->getUseMrgComplex() && ((rpcTempCU->getWidth(0) >= 16) || (rpcTempCU->getHeight(0) >= 16)))
    {
      TComSlice*    pSlice = rpcTempCU->getSlice();
      TComMvField  cMvFieldNeighbourNodes[MRG_NODE_NUM << 1];
      UChar uhInterDirNeighbourNodes[MRG_NODE_NUM];
      Int numExistMergeNode[2] = { 0 }; //numExistMergeNode[i] represents number of Exist merge nodes  of list i (i=0,1) without removing duplication
      Int numValidMergeNode[MRG_COMPLEX_NUM_CANDS << 1] = { 0 };//numValidMergeNode[i] represents number of valid merge nodes of list i  (i=0,1) after removing duplication
      UChar numValidMergeComplexCand = 0;
      Int RealComplexIndex[11] = { 0 };
      for (UInt ui = 0; ui < MRG_NODE_NUM; ++ui)
      {
        uhInterDirNeighbourNodes[ui] = 0;
      }

      rpcTempCU->getInterMergeComplexCandidates(0, 0, numValidMergeComplexCand, cMvFieldNeighbourNodes, uhInterDirNeighbourNodes, numExistMergeNode, numValidMergeNode, RealComplexIndex);

      if (numValidMergeComplexCand > 0)
      {
        UInt uiCostCand = MAX_UINT;
        UInt uiMergeComplexIndex = 0;
        UInt uiComplexInterDir = 0;

        NodeInfo NodeMotionInfo;
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          NodeMotionInfo.MvFieldNbNodes[i] = cMvFieldNeighbourNodes[i];
        }
        for (Int i = 0; i < MRG_NODE_NUM; i++)
        {
          NodeMotionInfo.InterDirNbNodes[i] = uhInterDirNeighbourNodes[i];
        }
        for (Int i = 0; i < (MRG_COMPLEX_NUM_CANDS << 1); i++)
        {
          NodeMotionInfo.numValidNodes[i] = numValidMergeNode[i];
        }
        for (Int i = 0; i < MRG_COMPLEX_NUM_CANDS;i++)
        {
          NodeMotionInfo.RealComplexIdx[i] = RealComplexIndex[i];
        }

#if JVET_C0024_QTBT
        rpcTempCU->setNodeMotionInfoSubParts(NodeMotionInfo, 0);
        rpcTempCU->setNumMergeComplexCandSubParts(numValidMergeComplexCand, 0);
        m_pcPredSearch->xMergeComplexEstimation(rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], 0, uiComplexInterDir, NodeMotionInfo, uiMergeComplexIndex, uiCostCand);
#else
        rpcTempCU->setNodeMotionInfoSubParts(NodeMotionInfo, 0, 0, uhDepth);
        rpcTempCU->setNumMergeComplexCandSubParts(numValidMergeComplexCand, 0, 0, uhDepth);
        m_pcPredSearch->xMergeComplexEstimation(rpcTempCU, m_ppcOrigYuv[uhDepth], 0, uiComplexInterDir, pacComplexCUMvField, NodeMotionInfo, uiMergeComplexIndex, uiCostCand);
#endif

        TComYuv* pcPredYuvTemp = NULL;
        UInt iteration;
        if (rpcTempCU->isLosslessCoded(0))
        {
          iteration = 1;
        }
        else
        {
          iteration = 2;
        }

        for (UInt uiNoResidual = 0; uiNoResidual < iteration; ++uiNoResidual)
        {
          DEBUG_STRING_NEW(tmpStr)
            // set MC parameters
            rpcTempCU->setPredModeSubParts(MODE_INTER, 0, uhDepth); // interprets depth relative to CTU level
          rpcTempCU->setCUTransquantBypassSubParts(bTransquantBypassFlag, 0, uhDepth);
          rpcTempCU->setChromaQpAdjSubParts(bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth);
#if !JVET_C0024_QTBT
          rpcTempCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uhDepth); // interprets depth relative to CTU level
#endif
          rpcTempCU->setMergeFlagSubParts(true, 0, 0, uhDepth); // interprets depth relative to CTU level
#if COM16_C806_OBMC
          rpcTempCU->setOBMCFlagSubParts(true, 0, uhDepth);
#endif
#if VCEG_AZ06_IC
          rpcTempCU->setICFlagSubParts(false, 0, uhDepth);
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          rpcTempCU->setMergeTypeSubParts(MGR_TYPE_DEFAULT_N, 0, 0, uhDepth);
#endif

#if JVET_C0024_QTBT
          rpcTempCU->setMergeComplexFlagSubParts(true, 0);
          rpcTempCU->setMergeComplexIndexSubParts(uiMergeComplexIndex, 0); // interprets depth relative to LCU level
          rpcTempCU->setNodeMotionInfoSubParts(NodeMotionInfo, 0);
          rpcTempCU->setNumMergeComplexCandSubParts(numValidMergeComplexCand, 0);
          rpcTempCU->setNodeMotionInfoSubParts(NodeMotionInfo, 0);
          rpcTempCU->setMotionInfoField(uiMergeComplexIndex, uiComplexInterDir, SIZE_2Nx2N, 0, 0);
#else
          rpcTempCU->setMergeComplexFlagSubParts(true, 0, 0, uhDepth);
          rpcTempCU->setMergeComplexIndexSubParts(uiMergeComplexIndex, 0, 0, uhDepth); // interprets depth relative to LCU level
          rpcTempCU->setNodeMotionInfoSubParts(NodeMotionInfo, 0, 0, uhDepth);
          rpcTempCU->setNumMergeComplexCandSubParts(numValidMergeComplexCand, 0, 0, uhDepth);
#endif
          rpcTempCU->setInterDirSubParts(uiComplexInterDir, 0, 0, 0); // interprets depth relative to LCU level

#if JVET_C0024_QTBT
                                                                      // do MC
          if (uiNoResidual == 0)
          {
            m_pcPredSearch->motionCompensation(rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx]);
#if COM16_C806_OBMC
            m_pcPredSearch->subBlockOBMC(rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx]);
#endif
            // save pred address
            pcPredYuvTemp = m_pppcPredYuvTemp[uiWIdx][uiHIdx];
          }
          else
          {
            if (pcPredYuvTemp != m_pppcPredYuvTemp[uiWIdx][uiHIdx])
            {
              //address changes take best (old temp)
              pcPredYuvTemp = m_pppcPredYuvBest[uiWIdx][uiHIdx];
            }
          }

          // estimate residual and encode everything
          m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU,
            m_pppcOrigYuv[uiWIdx][uiHIdx],
            m_pppcPredYuvTemp[uiWIdx][uiHIdx],
            m_pppcResiYuvTemp[uiWIdx][uiHIdx],
            m_pppcResiYuvBest[uiWIdx][uiHIdx],
            m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
            (uiNoResidual != 0)
#if COM16_C806_EMT
            , rpcBestCU->getTotalCost()
#endif
            DEBUG_STRING_PASS_INTO(tmpStr));
          rpcTempCU->setSkipFlagSubParts(rpcTempCU->getQtRootCbf(0) == 0, 0, 0);
#else
                                                                      // do MC
          if (uiNoResidual == 0)
          {
            m_pcPredSearch->motionCompensation(rpcTempCU, m_ppcPredYuvTemp[uhDepth]);
#if COM16_C806_OBMC
            m_pcPredSearch->subBlockOBMC(rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth]);
#endif
            // save pred address
            pcPredYuvTemp = m_ppcPredYuvTemp[uhDepth];
          }
          else
          {
            if (pcPredYuvTemp != m_ppcPredYuvTemp[uhDepth])
            {
              //address changes take best (old temp)
              pcPredYuvTemp = m_ppcPredYuvBest[uhDepth];
            }
          }

          // estimate residual and encode everything
          m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU,
            m_ppcOrigYuv[uhDepth],
            m_ppcPredYuvTemp[uhDepth],
            m_ppcResiYuvTemp[uhDepth],
            m_ppcResiYuvBest[uhDepth],
            m_ppcRecoYuvTemp[uhDepth],
            (uiNoResidual != 0)
#if COM16_C806_EMT
            , rpcBestCU->getTotalCost()
#endif
            DEBUG_STRING_PASS_INTO(tmpStr));
          rpcTempCU->setSkipFlagSubParts(rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth);
#endif
          Int orgQP = rpcTempCU->getQP(0);
          xCheckDQP(rpcTempCU);
          xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr));

          rpcTempCU->initEstData(uhDepth, orgQP, bTransquantBypassFlag);
          Bool bQtRootCbf = rpcTempCU->getQtRootCbf(0) == 1;
          if (!bQtRootCbf)
            break;
        }
      }
    }
#endif
}

#if VCEG_AZ07_FRUC_MERGE
Void TEncCu::xCheckRDCostMerge2Nx2NFRUC( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU , Bool *earlyDetectionSkipMode )
{
#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[rpcBestCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcBestCU->getHeight(0)];
#endif
  UChar uhDepth = rpcTempCU->getDepth( 0 );
#if THIRD_FRUC_MODE
  const UChar uhFRUCME[3] = { FRUC_MERGE_BILATERALMV , FRUC_MERGE_TEMPLATE,FRUC_MERGE_TEMPLATE_UNI };
#else
  const UChar uhFRUCME[2] = { FRUC_MERGE_BILATERALMV , FRUC_MERGE_TEMPLATE };
#endif
#if VCEG_AZ06_IC
  Bool bICFlag = rpcTempCU->getICFlag( 0 );
#endif
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
  Bool bFrucTemplateSearch = false;
#endif
#if THIRD_FRUC_MODE
  for (Int nME = 0; nME < 3; nME++)
#else
  for( Int nME = 0 ; nME < 2 ; nME++ )
#endif
  {
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
    if( m_pcPredSearch->getSaveLoadTag(rpcBestCU->getZorderIdxInCtu(), rpcBestCU->getCURouteIdx(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && 
      uhFRUCME[nME] != m_pcPredSearch->getSaveLoadFrucMode( rpcBestCU->getZorderIdxInCtu(), rpcBestCU->getCURouteIdx(), uiWIdx, uiHIdx ) )
#else
    if( m_pcPredSearch->getSaveLoadTag(rpcBestCU->getZorderIdxInCtu(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && 
      uhFRUCME[nME] != m_pcPredSearch->getSaveLoadFrucMode( uiWIdx, uiHIdx ) )
#endif
    {
      continue;
    }
#endif
#if !JVET_C0024_QTBT
    rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); 
#endif
    rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth ); // interprets depth relative to LCU level
    rpcTempCU->setCUTransquantBypassSubParts( false,     0, uhDepth );
    rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth ); // interprets depth relative to LCU level
    rpcTempCU->setMergeIndexSubParts( 0, 0, 0, uhDepth ); // interprets depth relative to LCU level
    rpcTempCU->setFRUCMgrModeSubParts( uhFRUCME[nME] , 0 , 0 , uhDepth );
#if VCEG_AZ06_IC
    rpcTempCU->setICFlagSubParts( bICFlag, 0, uhDepth );
#endif
#if JVECT_C0062_AFFINE_SIX_PARAM
    rpcTempCU->setAffineFlagSubParts(false, 0, 0, uhDepth);
    rpcTempCU->setAffineParamFlagSubParts(0, 0, 0, uhDepth);
#endif
    Bool bAvailable = m_pcPredSearch->deriveFRUCMV( rpcTempCU , uhDepth , 0 , 0
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
        ,-1,REF_PIC_LIST_0,NULL,bFrucTemplateSearch
#endif
        ); 
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
    if (nME == 1)
        bFrucTemplateSearch = bAvailable;
#endif
#if JVET_C0024_QTBT
    m_pppcFRUCBufferCU[uiWIdx][uiHIdx]->copyPartFrom( rpcTempCU , 0 , uhDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0) );
#else
    m_ppcFRUCBufferCU[uhDepth]->copyPartFrom( rpcTempCU , 0 , uhDepth );
#endif
    if( bAvailable )
    {
      UInt iteration = 1 + !rpcTempCU->isLosslessCoded(0) ;
      for( UInt uiNoResidual = 0; uiNoResidual < iteration; uiNoResidual++ )
      {
        if( uiNoResidual > 0 )
        {
#if JVET_C0024_QTBT
            rpcTempCU->copyPartFrom( m_pppcFRUCBufferCU[uiWIdx][uiHIdx] , 0 , uhDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0) );
#else
          rpcTempCU->copyPartFrom( m_ppcFRUCBufferCU[uhDepth] , 0 , uhDepth );
#endif
        }
        // do MC
#if JVET_C0024_QTBT
        m_pcPredSearch->motionCompensation ( rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx] );
#if COM16_C806_OBMC
        m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
        // estimate residual and encode everything
        m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
          m_pppcOrigYuv    [uiWIdx][uiHIdx],
          m_pppcPredYuvTemp[uiWIdx][uiHIdx],
          m_pppcResiYuvTemp[uiWIdx][uiHIdx],
          m_pppcResiYuvBest[uiWIdx][uiHIdx],
          m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
          (uiNoResidual? true:false)
#if COM16_C806_EMT
          , rpcBestCU->getTotalCost()
#endif
          );
#else
        m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[uhDepth] );
#if COM16_C806_OBMC
        m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth] );
#endif
        // estimate residual and encode everything
        m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
          m_ppcOrigYuv    [uhDepth],
          m_ppcPredYuvTemp[uhDepth],
          m_ppcResiYuvTemp[uhDepth],
          m_ppcResiYuvBest[uhDepth],
          m_ppcRecoYuvTemp[uhDepth],
          (uiNoResidual? true:false)
#if COM16_C806_EMT
          , rpcBestCU->getTotalCost()
#endif
          );
#endif
        if ( uiNoResidual == 0 && rpcTempCU->getQtRootCbf(0) == 0 )
        {
          uiNoResidual++;
        }
        rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth );
        Int orgQP = rpcTempCU->getQP( 0 );
        xCheckDQP( rpcTempCU );
#if BIF_NO_SKIP_ONLY
        if (rpcTempCU->getSIFFlag(0) && rpcTempCU->isSkipped(0))
        {
          rpcTempCU->getTotalCost() = DBL_MAX;
        }
#endif
        xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
        rpcTempCU->initEstData( uhDepth, orgQP, false );
      }
    }
  }
}
#endif

#if AMP_MRG
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseMRG
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
  , UChar bIMV , TComDataCU * pcCUInfo2Reuse 
#else
  , Bool bIMV , TComDataCU * pcCUInfo2Reuse 
#endif
#endif
  )
#else
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize 
#if VCEG_AZ07_IMV
  , Bool bIMV , TComDataCU * pcCUInfo2Reuse 
#endif
  )
#endif
{
#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[rpcBestCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcBestCU->getHeight(0)];
#endif

  DEBUG_STRING_NEW(sTest)

  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
#if JVET_C0024_QTBT
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMinQTSize(rpcBestCU->getSlice()->getSliceType(), rpcBestCU->getTextType()), sps.getCTUSize(), 32u);
#else
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight()>>(sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
#endif
    if(ePartSize != SIZE_2Nx2N || rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Inter in fast deltaqp mode
    }
  }

  // prior to this, rpcTempCU will have just been reset using rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
  UChar uhDepth = rpcTempCU->getDepth( 0 );

#if COM16_C806_LARGE_CTU && !JVET_C0024_QTBT
  if( m_pcEncCfg->getUseFastLCTU() )
  {
    if( ePartSize != SIZE_2Nx2N && rpcTempCU->getWidth( 0 ) > 64 )
    {
      rpcTempCU->getTotalCost() = MAX_DOUBLE / 4;
      rpcTempCU->getTotalDistortion() = MAX_INT;
      xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
      return;
    }
  }
#endif

#if !JVET_C0024_QTBT
  rpcTempCU->setPartSizeSubParts  ( ePartSize,  0, uhDepth );
#endif
  rpcTempCU->setPredModeSubParts  ( MODE_INTER, 0, uhDepth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );
#if COM16_C806_OBMC
  rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
#endif
#if VCEG_AZ06_IC
  Bool bICFlag = rpcTempCU->getICFlag( 0 );
#endif
#if VCEG_AZ07_IMV
  rpcTempCU->setiMVFlagSubParts( bIMV,  0, uhDepth );
  if( bIMV && pcCUInfo2Reuse != NULL )
  {
    // reuse the motion info from pcCUInfo2Reuse
#if JVET_C0024_QTBT
      rpcTempCU->copyPartFrom( pcCUInfo2Reuse , 0 , uhDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0) );
#else
    assert( pcCUInfo2Reuse->getPartitionSize( 0 ) == ePartSize );
    rpcTempCU->copyPartFrom( pcCUInfo2Reuse , 0 , uhDepth );
#endif
    rpcTempCU->setiMVFlagSubParts( bIMV , 0 , uhDepth );
    rpcTempCU->resetMVDandMV2Int( true 
#if VCEG_AZ07_FRUC_MERGE
      , m_pcPredSearch
#endif
      );
#if VCEG_AZ06_IC
    bICFlag = rpcTempCU->getICFlag( 0 );
#endif
    if( !rpcTempCU->hasSubCUNonZeroMVd() )
    {
      return;
    }
    else
    {
#if JVET_C0024_QTBT
      m_pcPredSearch->motionCompensation( rpcTempCU , m_pppcPredYuvTemp[uiWIdx][uiHIdx] );
#if COM16_C806_OBMC
      m_pppcPredYuvTemp[uiWIdx][uiHIdx]->copyToPartYuv( m_pppcPredYuvWoOBMC[uiWIdx][uiHIdx] , 0 );
      rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
      m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
#else
      m_pcPredSearch->motionCompensation( rpcTempCU , m_ppcPredYuvTemp[uhDepth] );
#if COM16_C806_OBMC
      m_ppcPredYuvTemp[uhDepth]->copyToPartYuv( m_ppcPredYuvWoOBMC[uhDepth] , 0 );
      rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
      m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth] );
#endif
#endif
    }
  }
  else
  {
#endif
#if AMP_MRG
#if !JVET_C0024_QTBT
  rpcTempCU->setMergeAMP (true);
#endif
#if COM16_C806_OBMC
#if JVET_C0024_QTBT
#if AFF_EARLY_TERMINATION
	m_pcPredSearch->predInterSearch(rpcTempCU, rpcBestCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], m_pppcPredYuvWoOBMC[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx], false, bUseMRG);
#else
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], m_pppcPredYuvWoOBMC[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx], false, bUseMRG );
#endif
#else
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth], m_ppcPredYuvWoOBMC[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth], false, bUseMRG );
#endif
#else
#if JVET_C0024_QTBT
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx] , false, bUseMRG );
#else
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] DEBUG_STRING_PASS_INTO(sTest), false, bUseMRG );
#endif
#endif
#else
#if COM16_C806_OBMC
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth], m_ppcPredYuvWoOBMC[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth] );
#else
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] );
#endif
#endif

#if AMP_MRG && !JVET_C0024_QTBT
  if ( !rpcTempCU->getMergeAMP() )
  {
    return;
  }
#endif

#if VCEG_AZ07_IMV
  if( bIMV )
  {
    if( !rpcTempCU->hasSubCUNonZeroMVd() )
    {
      return;
    }
  }
  }
#endif

#if COM16_C806_OBMC
#if JVET_C0024_QTBT
  m_pppcTempCUWoOBMC[uiWIdx][uiHIdx]->initEstData( uhDepth, rpcTempCU->getQP( 0 ), rpcTempCU->getCUTransquantBypass( 0 ) );
  m_pppcTempCUWoOBMC[uiWIdx][uiHIdx]->copyPartFrom( rpcTempCU, 0, uhDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0) );
#else
  m_ppcTempCUWoOBMC[uhDepth]->initEstData( uhDepth, rpcTempCU->getQP( 0 ), rpcTempCU->getCUTransquantBypass( 0 ) );
  m_ppcTempCUWoOBMC[uhDepth]->copyPartFrom( rpcTempCU, 0, uhDepth );
#endif
  Bool bCheckOBMC[2] = { true , true };
  bCheckOBMC[0] = rpcTempCU->isOBMCFlagCoded( 0 ); // check OBMC off only when the flag is to be coded
  if( !rpcTempCU->getSlice()->getSPS()->getOBMC() )
  {
    bCheckOBMC[1] = false; bCheckOBMC[0] = true;
  }
  else if( rpcTempCU->isOBMCFlagCoded( 0 ) )
  {
    const Double dOBMCThOn  = 0.0;
    const Double dOBMCThOff = 1.0;
#if JVET_C0024_QTBT
    UInt uiSADOBMCOff = m_pppcOrigYuv[uiWIdx][uiHIdx]->sadLuma( m_pppcPredYuvWoOBMC[uiWIdx][uiHIdx] );
    UInt uiSADOBMCOn  = m_pppcOrigYuv[uiWIdx][uiHIdx]->sadLuma( m_pppcPredYuvTemp[uiWIdx][uiHIdx] );
#else
    UInt uiSADOBMCOff = m_ppcOrigYuv[uhDepth]->sadLuma( m_ppcPredYuvWoOBMC[uhDepth] );
    UInt uiSADOBMCOn  = m_ppcOrigYuv[uhDepth]->sadLuma( m_ppcPredYuvTemp[uhDepth] );
#endif
    // OBMC off
    bCheckOBMC[0] = uiSADOBMCOff * dOBMCThOff < uiSADOBMCOn;
    // OBMC on
    bCheckOBMC[1] = !bCheckOBMC[0] || uiSADOBMCOn * dOBMCThOn < uiSADOBMCOff;
  }
  for( Int nOBMC = 1 ; nOBMC >= 0 ; nOBMC-- )
  {
    if( !bCheckOBMC[nOBMC] )
    {
      continue;
    }
#if JVET_C0024_QTBT
    rpcTempCU->copyPartFrom( m_pppcTempCUWoOBMC[uiWIdx][uiHIdx], 0, uhDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0) );
    if (nOBMC == 0)
      m_pppcPredYuvWoOBMC[uiWIdx][uiHIdx]->copyToPartYuv(m_pppcPredYuvTemp[uiWIdx][uiHIdx], 0);
#else
    rpcTempCU->copyPartFrom( m_ppcTempCUWoOBMC[uhDepth], 0, uhDepth );
    if (nOBMC == 0)
      m_ppcPredYuvWoOBMC[uhDepth]->copyToPartYuv(m_ppcPredYuvTemp[uhDepth], 0);
#endif
    rpcTempCU->setOBMCFlagSubParts( ( Bool )nOBMC , 0 , uhDepth );
#endif
#if VCEG_AZ06_IC
    rpcTempCU->setICFlagSubParts( bICFlag, 0, uhDepth );
#endif
#if JVET_C0024_QTBT
    m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], 

    m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvBest[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], false 
#else
    m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], 

    m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false 
#endif
#if COM16_C806_EMT
    , rpcBestCU->getTotalCost() 
#endif
    DEBUG_STRING_PASS_INTO(sTest) 
  );
  rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
#if VCEG_AZ07_IMV
#if !JVET_C0024_QTBT
  if( rpcTempCU->getiMVFlag( 0 ) == 0 )
  {
    if( rpcTempCU->getTotalCost() < m_ppcTempCUIMVCache[ePartSize][uhDepth]->getTotalCost() )
    {
      SModeCand tmpCand;
      tmpCand.eInterPartSize = ePartSize;
      tmpCand.bUseMrg = bUseMRG;
      tmpCand.dRDCost = rpcTempCU->getTotalCost();
      m_ppcTempCUIMVCache[ePartSize][uhDepth]->copyPartFrom( rpcTempCU , 0 , uhDepth );
      m_ppcTempCUIMVCache[ePartSize][uhDepth]->getTotalCost() = tmpCand.dRDCost;
      tmpCand.pcCUMode = m_ppcTempCUIMVCache[ePartSize][uhDepth];
      m_setInterCand.insert( tmpCand );
    }
  }
#endif
#endif

#if JVET_E0076_MULTI_PEL_MVD
  if (rpcTempCU->getTotalCost() < m_dBestMvDPelCost[rpcTempCU->getiMVFlag(0)])
  {
    m_dBestMvDPelCost[rpcTempCU->getiMVFlag(0)] = rpcTempCU->getTotalCost();
  }
#endif

#if DEBUG_STRING
  DebugInterPredResiReco(sTest, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

#if JVET_C0024_DELTA_QP_FIX && COM16_C806_OBMC
  Int orgQP = rpcTempCU->getQP( 0 );
#endif
  xCheckDQP( rpcTempCU );
  xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
#if COM16_C806_OBMC
#if JVET_C0024_DELTA_QP_FIX
  rpcTempCU->initEstData( uhDepth, orgQP, rpcTempCU->getCUTransquantBypass( 0 ) );
#else
  rpcTempCU->initEstData( uhDepth, rpcTempCU->getQP( 0 ), rpcTempCU->getCUTransquantBypass( 0 ) );
#endif
  }
#endif
}

#if VCEG_AZ08_INTER_KLT
Void TEncCu::xCheckRDCostInterKLT(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize)
{
    DEBUG_STRING_NEW(sTest)

#if JVET_C0024_QTBT
    if( rpcTempCU->getWidth(0) != rpcTempCU->getHeight(0) )
    {
      return;
    }
#endif

    if(getFastDeltaQp())
    {
        const TComSPS &sps=*( rpcTempCU->getSlice()->getSPS());
#if JVET_C0024_QTBT
        const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMinQTSize(rpcBestCU->getSlice()->getSliceType(), rpcBestCU->getTextType()), sps.getCTUSize(), 32u);
#else
        const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight()>>(sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
#endif
        if(ePartSize != SIZE_2Nx2N || rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxSize)
        {
            return; // only check necessary 2Nx2N Inter in fast deltaqp mode
        }
    }

    // prior to this, rpcTempCU will have just been reset using rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
    UChar uhDepth = rpcTempCU->getDepth(0);

#if COM16_C806_LARGE_CTU
    if (m_pcEncCfg->getUseFastLCTU())
    {
        if (ePartSize != SIZE_2Nx2N && rpcTempCU->getWidth(0) > 64)
        {
            rpcTempCU->getTotalCost() = MAX_DOUBLE / 4;
            rpcTempCU->getTotalDistortion() = MAX_INT;
            xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
            return;
        }
    }
#endif

    Bool bSkipPossible = false;
    rpcTempCU->copySameSizeCUFrom(rpcBestCU, 0, uhDepth);
    UInt uiWidth = rpcTempCU->getWidth(0);
    UInt uiHeigh = rpcTempCU->getHeight(0);
#if JVET_C0024_QTBT
    UInt uiWIdx = g_aucConvertToBit[uiWidth];
    UInt uiHIdx = g_aucConvertToBit[uiHeigh];
    Pel *pYSrc = m_pppcPredYuvBest[uiWIdx][uiHIdx]->getAddr(COMPONENT_Y);
    Pel *pYDst = m_pppcPredYuvTemp[uiWIdx][uiHIdx]->getAddr(COMPONENT_Y);
#else
    Pel *pYSrc = m_ppcPredYuvBest[uhDepth]->getAddr(COMPONENT_Y);
    Pel *pYDst = m_ppcPredYuvTemp[uhDepth]->getAddr(COMPONENT_Y);
#endif
    memcpy(pYDst, pYSrc, sizeof(Pel)*uiWidth*uiHeigh);
    const UInt componentShiftCb = rpcTempCU->getPic()->getComponentScaleX(COMPONENT_Cb) + rpcTempCU->getPic()->getComponentScaleY(COMPONENT_Cb);
#if JVET_C0024_QTBT
    pYSrc = m_pppcPredYuvBest[uiWIdx][uiHIdx]->getAddr(COMPONENT_Cb);
    pYDst = m_pppcPredYuvTemp[uiWIdx][uiHIdx]->getAddr(COMPONENT_Cb);
#else
    pYSrc = m_ppcPredYuvBest[uhDepth]->getAddr(COMPONENT_Cb);
    pYDst = m_ppcPredYuvTemp[uhDepth]->getAddr(COMPONENT_Cb);
#endif
    memcpy(pYDst, pYSrc, sizeof(Pel)* uiWidth * uiHeigh >> componentShiftCb);
    const UInt componentShiftCr = rpcTempCU->getPic()->getComponentScaleX(COMPONENT_Cb) + rpcTempCU->getPic()->getComponentScaleY(COMPONENT_Cb);
#if JVET_C0024_QTBT
    pYSrc = m_pppcPredYuvBest[uiWIdx][uiHIdx]->getAddr(COMPONENT_Cr);
    pYDst = m_pppcPredYuvTemp[uiWIdx][uiHIdx]->getAddr(COMPONENT_Cb);
#else
    pYSrc = m_ppcPredYuvBest[uhDepth]->getAddr(COMPONENT_Cr);
    pYDst = m_ppcPredYuvTemp[uhDepth]->getAddr(COMPONENT_Cb);
#endif
    memcpy(pYDst, pYSrc, sizeof(Pel)*uiWidth*uiHeigh >> componentShiftCr);
    bSkipPossible = rpcBestCU->getSkipFlag(0);
#if AMP_MRG && !JVET_C0024_QTBT
    if (!rpcTempCU->getMergeAMP())
    {
        return;
    }
#endif

#if COM16_C806_OBMC //QC_OBMC
#if JVET_C0024_QTBT
    m_pcPredSearch->motionCompensation(rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx]);
    rpcTempCU->setOBMCFlagSubParts(true, 0, uhDepth);
    m_pcPredSearch->subBlockOBMC(rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx]);
#else
    m_pcPredSearch->motionCompensation(rpcTempCU, m_ppcPredYuvTemp[uhDepth]);
    rpcTempCU->setOBMCFlagSubParts(true, 0, uhDepth);
    m_pcPredSearch->subBlockOBMC(rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth]);
#endif
#endif

#if JVET_C0024_QTBT
    m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], 
        m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvBest[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], false
#else
    m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], 
        m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false
#endif
#if COM16_C806_EMT
        , rpcBestCU->getTotalCost()
#endif
        DEBUG_STRING_PASS_INTO(sTest)
        );
    if (bSkipPossible)
    {
        rpcTempCU->setSkipFlagSubParts(rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth);
    }
    xCheckDQP(rpcTempCU);
    xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
}
#endif

Void TEncCu::xCheckRDCostIntra( TComDataCU *&rpcBestCU,
                                TComDataCU *&rpcTempCU,
                                Double      &cost,
                                PartSize     eSize
                                DEBUG_STRING_FN_DECLARE(sDebug)
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
                                , Int& bNonZeroCoeff
#endif
#if DIMD_INTRA_PRED
  , UChar ucDIMDIntraRD
#endif
// #if DIMD_MLIP_HARM
//     , Bool bUseDIMD
// #endif
                                )
{
  DEBUG_STRING_NEW(sTest)

  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
#if JVET_C0024_QTBT
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMinQTSize(rpcBestCU->getSlice()->getSliceType(), rpcBestCU->getTextType()), sps.getCTUSize(), 32u);
#else
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight()>>(sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
#endif
    if(rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Intra in fast deltaqp mode
    }
  }

#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[rpcBestCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcBestCU->getHeight(0)];
#endif
  UInt uiDepth = rpcTempCU->getDepth( 0 );
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
  Int iMaxQtLv = g_aucConvertToBit[rpcTempCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[isLuma(rpcTempCU->getTextType()) ? rpcTempCU->getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
  const UInt uiDIMDMaxBtDepht = (rpcTempCU->getBoundaryCtuLocation() && rpcTempCU->getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
  const UInt uiDIMDMaxBtDepht = (rpcTempCU->getBoundaryCtuLocation() && rpcTempCU->getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif

#if COM16_C806_LARGE_CTU
  if( m_pcEncCfg->getUseFastLCTU() ) 
  {
#if JVET_C0024_QTBT
    if( rpcTempCU->getWidth( 0 )*rpcTempCU->getHeight(0) > 4096 )
#else
    if( rpcTempCU->getWidth( 0 ) > 64 )
#endif
    {
      return;
    }
  }
#endif

#if DIMD_INTRA_PRED
  rpcTempCU->setDIMDEnabledFlagSubParts(rpcTempCU->getTextType(), ucDIMDIntraRD, 0, uiDepth);
  if (!rpcTempCU->getSlice()->isIntra())
  {
    //assert(rpcTempCU->getTextType() == CHANNEL_TYPE_LUMA);
#if DIMD_DISABLE_CHROMA_P_B_SLICE_JB	 
    rpcTempCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, false, 0, uiDepth);
#else
    rpcTempCU->setDIMDEnabledFlagSubParts(CHANNEL_TYPE_CHROMA, ucDIMDIntraRD, 0, uiDepth);
#endif
  }
  if (isLuma(rpcTempCU->getTextType()))
  {
#if BT_FORCE_DIMD_HARM
    if (rpcTempCU->getBTDepth(0) <= uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra())
#else
    if (rpcTempCU->getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra())
#endif
    {
      assert(!rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0));
    }
    if (!rpcTempCU->getSlice()->isIntra())
    {
#if DIMD_DISABLE_CHROMA_P_B_SLICE_JB
      rpcTempCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, false, 0, uiDepth);
#else
      rpcTempCU->setDIMDNoBTLevelFlagSubParts(CHANNEL_TYPE_CHROMA, rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_LUMA, 0), 0, uiDepth);
#endif
    }
  }
  else
  {
    assert(!rpcTempCU->getDIMDNoBTLevelFlag(CHANNEL_TYPE_CHROMA, 0));
  }
#endif

  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );
#if COM16_C806_OBMC
  rpcTempCU->setOBMCFlagSubParts( false, 0, uiDepth );
#endif
#if !JVET_C0024_QTBT
  rpcTempCU->setPartSizeSubParts( eSize, 0, uiDepth );
#endif
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth );

#if !COM16_C806_LARGE_CTU
  Pel resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
#endif

#if JVET_C0024_QTBT
  if (isLuma(rpcBestCU->getTextType()))
  {
#if DIMD_INTRA_PRED
    if (rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0))
    {
      m_pcPredSearch->estDIMDIntraPredLumaQT(rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
#if COM16_C806_LARGE_CTU
        m_resiBuffer
#else
        resiLuma
#endif
        );
    }
    else
    {
#endif
  m_pcPredSearch->estIntraPredLumaQT( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], 
// #if DIMD_MLIP_HARM
//       bUseDIMD,
// #endif
#else
  m_pcPredSearch->estIntraPredLumaQT( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], 
#endif
#if COM16_C806_LARGE_CTU
    m_resiBuffer 
#else
    resiLuma 
#endif
    DEBUG_STRING_PASS_INTO(sTest) );
#if JVET_C0024_PBINTRA_FAST
  if (rpcTempCU->getTotalDistortion()==MAX_UINT)
  {
    assert(rpcTempCU->getInterHAD()==0);
    return;
  }
#endif
#if DIMD_INTRA_PRED
    }
#endif

#if JVET_C0024_QTBT
  if (!rpcBestCU->getSlice()->isIntra())
  {
  m_pppcRecoYuvTemp[uiWIdx][uiHIdx]->copyToPicComponent(COMPONENT_Y, rpcTempCU->getPic()->getPicYuvRec(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu() );
#else
  m_ppcRecoYuvTemp[uiDepth]->copyToPicComponent(COMPONENT_Y, rpcTempCU->getPic()->getPicYuvRec(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu() );
#endif
#if JVET_C0024_QTBT
  }
  }
  if (rpcBestCU->getPic()->getChromaFormat()!=CHROMA_400 
    && (isChroma(rpcBestCU->getTextType()) || !rpcBestCU->getSlice()->isIntra()))
  {
#if DIMD_INTRA_PRED
    if (rpcTempCU->getDIMDEnabledFlag(CHANNEL_TYPE_CHROMA, 0))
    {
      m_pcPredSearch->estDIMDIntraPredChromaQT(rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
#if COM16_C806_LARGE_CTU
        m_resiBuffer
#else
        resiLuma
#endif
        );
    }
    else
    {
#endif
    m_pcPredSearch->estIntraPredChromaQT( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx], 
#else
  if (rpcBestCU->getPic()->getChromaFormat()!=CHROMA_400)
  {
    m_pcPredSearch->estIntraPredChromaQT( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], 
#endif
#if COM16_C806_LARGE_CTU
      m_resiBuffer 
#else
      resiLuma 
#endif
      DEBUG_STRING_PASS_INTO(sTest) 
      );
#if DIMD_INTRA_PRED
    }
#endif
  }

  m_pcEntropyCoder->resetBits();

  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }

  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode( rpcTempCU, 0,          true );
#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
  if (isLuma(rpcTempCU->getTextType()) && (rpcTempCU->getBTDepth(0) <= uiDIMDMaxBtDepht || !rpcTempCU->getSlice()->isIntra()))
#else
  if (isLuma(rpcTempCU->getTextType()) && (rpcTempCU->getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI || !rpcTempCU->getSlice()->isIntra()))
#endif
  {
    m_pcEntropyCoder->encodeDIMDFlag(rpcTempCU, 0, uiDepth, rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), true);
  }
#endif
#if JVET_C0024_QTBT
  if (isLuma(rpcBestCU->getTextType()))
  {
#endif
#if VCEG_AZ05_INTRA_MPI
  m_pcEntropyCoder->encodeMPIIdx(rpcTempCU, 0, true);
#endif 

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  m_pcEntropyCoder->encodePDPCIdx(rpcTempCU, 0, true);
#endif
#if JVET_C0024_QTBT
    }
#else
  m_pcEntropyCoder->encodePartSize( rpcTempCU, 0, uiDepth, true );
#endif
  m_pcEntropyCoder->encodePredInfo( rpcTempCU, 0 );
  m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true );

  // Encode Coefficients
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
#if JVET_C0024_QTBT
  bNonZeroCoeff = 0;
#else
  bNonZeroCoeff = false;
#endif
#endif

#if QTBT_NSST
  UChar ucNsstIdx = rpcTempCU->getROTIdx( rpcTempCU->getTextType(), 0 );
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  Int   iNonZeroCoeffNonTs;
#endif

  Bool bCodeDQP = getdQPFlag();
  Bool codeChromaQpAdjFlag = getCodeChromaQpAdjFlag();
  m_pcEntropyCoder->encodeCoeff( rpcTempCU, 0, uiDepth, bCodeDQP, codeChromaQpAdjFlag
#if VCEG_AZ05_ROT_TR || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
    , bNonZeroCoeff
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
    , iNonZeroCoeffNonTs
#endif
    );
  setCodeChromaQpAdjFlag( codeChromaQpAdjFlag );
  setdQPFlag( bCodeDQP );

#if JVET_C0024_QTBT
  m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);
#endif

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

#if JVET_C0024_DELTA_QP_FIX
  xCheckDQP( rpcTempCU );
#endif
#if QTBT_NSST
  const Int iNonZeroCoeffThr = isLuma(rpcTempCU->getTextType()) ? NSST_SIG_NZ_LUMA + (rpcTempCU->getSlice()->isIntra() ? 0 : NSST_SIG_NZ_CHROMA) : NSST_SIG_NZ_CHROMA;

#if JVET_G0104_PLANAR_PDPC
	if (ucNsstIdx && iNonZeroCoeffNonTs <= iNonZeroCoeffThr)
	{
			Bool isMDIS = false;
			if (isLuma(rpcTempCU->getTextType()))
			{
#if DIMD_NUM_INTRA_DIR_INC
					isMDIS = TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, rpcTempCU->getIntraDir(CHANNEL_TYPE_LUMA, 0), rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), rpcTempCU->getSlice()->getSPS()->getChromaFormatIdc(), rpcTempCU->getSlice()->getSPS()->getSpsRangeExtension().getIntraSmoothingDisabledFlag(), rpcTempCU, 0);
#else
					isMDIS = TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, rpcTempCU->getIntraDir(CHANNEL_TYPE_LUMA, 0), rpcTempCU->getWidth(0), rpcTempCU->getHeight(0), rpcTempCU->getSlice()->getSPS()->getChromaFormatIdc(), rpcTempCU->getSlice()->getSPS()->getSpsRangeExtension().getIntraSmoothingDisabledFlag());
#endif
			}

			if (iNonZeroCoeffNonTs > 0 || isMDIS)
			{
					rpcTempCU->getTotalCost() = MAX_DOUBLE;
			}
	}
#else
#if JVET_C0045_C0053_NO_NSST_FOR_TS
  if ( ucNsstIdx && iNonZeroCoeffNonTs <= iNonZeroCoeffThr && iNonZeroCoeffNonTs>0 )
#else
  if ( ucNsstIdx && bNonZeroCoeff <= iNonZeroCoeffThr && bNonZeroCoeff>0 )
#endif
  {
    rpcTempCU->getTotalCost() = MAX_DOUBLE;
  }
#endif
#endif

#if !JVET_C0024_DELTA_QP_FIX
  xCheckDQP( rpcTempCU );
#endif

  cost = rpcTempCU->getTotalCost();

#if JVET_C0024_QTBT
  UInt uiHAD = rpcTempCU->getInterHAD();
#endif
  xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
#if JVET_C0024_QTBT
  rpcTempCU->getInterHAD() = uiHAD;
#endif
}


/** Check R-D costs for a CU with PCM mode.
 * \param rpcBestCU pointer to best mode CU data structure
 * \param rpcTempCU pointer to testing mode CU data structure
 * \returns Void
 *
 * \note Current PCM implementation encodes sample values in a lossless way. The distortion of PCM mode CUs are zero. PCM mode is selected if the best mode yields bits greater than that of PCM mode.
 */
Void TEncCu::xCheckIntraPCM( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU )
{
  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
    const UInt fastDeltaQPCuMaxPCMSize = Clip3((UInt)1<<sps.getPCMLog2MinSize(), (UInt)1<<sps.getPCMLog2MaxSize(), 32u);
    if (rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxPCMSize)
    {
      return;   // only check necessary PCM in fast deltaqp mode
    }
  }
  
  UInt uiDepth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );
#if COM16_C806_OBMC
  rpcTempCU->setOBMCFlagSubParts( false, 0, uiDepth );
#endif
  rpcTempCU->setIPCMFlag(0, true);
  rpcTempCU->setIPCMFlagSubParts (true, 0, rpcTempCU->getDepth(0));
#if !JVET_C0024_QTBT
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth );
#endif
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
#if !JVET_C0024_QTBT
  rpcTempCU->setTrIdxSubParts ( 0, 0, uiDepth );
#endif
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth );

#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[rpcTempCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcTempCU->getHeight(0)];
  m_pcPredSearch->IPCMSearch( rpcTempCU, m_pppcOrigYuv[uiWIdx][uiHIdx], m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcResiYuvTemp[uiWIdx][uiHIdx], m_pppcRecoYuvTemp[uiWIdx][uiHIdx]);

  m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
  m_pcPredSearch->IPCMSearch( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth]);

  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
#endif

  m_pcEntropyCoder->resetBits();

  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }

  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode ( rpcTempCU, 0,          true );
#if !JVET_C0024_QTBT
  m_pcEntropyCoder->encodePartSize ( rpcTempCU, 0, uiDepth, true );
#endif
  m_pcEntropyCoder->encodeIPCMInfo ( rpcTempCU, 0, true );

#if JVET_C0024_QTBT
  m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);
#endif

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );
  DEBUG_STRING_NEW(a)
  DEBUG_STRING_NEW(b)
  xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(a) DEBUG_STRING_PASS_INTO(b));
}

/** check whether current try is the best with identifying the depth of current try
 * \param rpcBestCU
 * \param rpcTempCU
 * \param uiDepth
 */
#if JVET_C0024_QTBT
Void TEncCu::xCheckBestMode( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth, UInt uiWidth, UInt uiHeight DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo) )
#else
Void TEncCu::xCheckBestMode( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo) )
#endif
{
  if( rpcTempCU->getTotalCost() < rpcBestCU->getTotalCost() )
  {
    TComYuv* pcYuv;
    // Change Information data
    TComDataCU* pcCU = rpcBestCU;
    rpcBestCU = rpcTempCU;
    rpcTempCU = pcCU;

#if JVET_C0024_QTBT
    if (uiWidth==0 || uiHeight==0)
    {
      uiWidth = rpcTempCU->getWidth(0);
      uiHeight = rpcTempCU->getHeight(0);
    }
    UInt uiWIdx = g_aucConvertToBit[uiWidth]; 
    UInt uiHIdx = g_aucConvertToBit[uiHeight];
    // Change Prediction data
    pcYuv = m_pppcPredYuvBest[uiWIdx][uiHIdx];
    m_pppcPredYuvBest[uiWIdx][uiHIdx] = m_pppcPredYuvTemp[uiWIdx][uiHIdx];
    m_pppcPredYuvTemp[uiWIdx][uiHIdx] = pcYuv;

    // Change Reconstruction data
    pcYuv = m_pppcRecoYuvBest[uiWIdx][uiHIdx];
    m_pppcRecoYuvBest[uiWIdx][uiHIdx] = m_pppcRecoYuvTemp[uiWIdx][uiHIdx];
    m_pppcRecoYuvTemp[uiWIdx][uiHIdx] = pcYuv;
#else
    // Change Prediction data
    pcYuv = m_ppcPredYuvBest[uiDepth];
    m_ppcPredYuvBest[uiDepth] = m_ppcPredYuvTemp[uiDepth];
    m_ppcPredYuvTemp[uiDepth] = pcYuv;

    // Change Reconstruction data
    pcYuv = m_ppcRecoYuvBest[uiDepth];
    m_ppcRecoYuvBest[uiDepth] = m_ppcRecoYuvTemp[uiDepth];
    m_ppcRecoYuvTemp[uiDepth] = pcYuv;
#endif

    pcYuv = NULL;
    pcCU  = NULL;

    // store temp best CI for next CU coding
#if JVET_C0024_QTBT
    m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_NEXT_BEST]);
#else
    m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
#endif


#if DEBUG_STRING
    DEBUG_STRING_SWAP(sParent, sTest)
    const PredMode predMode=rpcBestCU->getPredictionMode(0);
    if ((DebugOptionList::DebugString_Structure.getInt()&DebugStringGetPredModeMask(predMode)) && bAddSizeInfo)
    {
      std::stringstream ss(stringstream::out);
      ss <<"###: " << (predMode==MODE_INTRA?"Intra   ":"Inter   ") << partSizeToString[rpcBestCU->getPartitionSize(0)] << " CU at " << rpcBestCU->getCUPelX() << ", " << rpcBestCU->getCUPelY() << " width=" << UInt(rpcBestCU->getWidth(0)) << std::endl;
      sParent+=ss.str();
    }
#endif
  }

#if TRACE_YUV
#if JVET_C0024_QTBT
  if (uiWidth == 0 || uiHeight == 0)
  {
    uiWidth = rpcTempCU->getWidth(0);
    uiHeight = rpcTempCU->getHeight(0);
  }
  UInt uiWIdx = g_aucConvertToBit[uiWidth]; 
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  TRACE_CU_PART(m_pppcOrigYuv[uiWIdx][uiHIdx],     rpcBestCU, "t_orig.yuv");
  TRACE_CU_PART(m_pppcPredYuvBest[uiWIdx][uiHIdx], rpcBestCU, "t_predicted.yuv");
  TRACE_CU_PART(m_pppcRecoYuvBest[uiWIdx][uiHIdx], rpcBestCU, "t_rec_from_best.yuv");
#else
  TRACE_CU_PART(m_ppcOrigYuv[uiDepth],     rpcBestCU, "t_orig.yuv");
  TRACE_CU_PART(m_ppcPredYuvBest[uiDepth], rpcBestCU, "t_predicted.yuv");
  TRACE_CU_PART(m_ppcRecoYuvBest[uiDepth], rpcBestCU, "t_rec_from_best.yuv");
#endif
#endif
}

Void TEncCu::xCheckDQP( TComDataCU* pcCU )
{
  UInt uiDepth = pcCU->getDepth( 0 );

  const TComPPS &pps = *(pcCU->getSlice()->getPPS());
#if JVET_C0024_DELTA_QP_FIX
  UInt uiBTDepth = pcCU->getBTDepth( 0 );
#if TT_TREE_DS1
  UInt uiQTBTDepth = (uiDepth<<1) + pcCU->getDeltaSizeIdx( 0 );
#else
  UInt uiQTBTDepth = (uiDepth<<1) + uiBTDepth;
#endif
  UInt uiMaxDQPDepthQTBT = pps.getMaxCuDQPDepth() << 1;
#endif
#if JVET_C0024_DELTA_QP_FIX
  if ( pps.getUseDQP() && uiQTBTDepth <= uiMaxDQPDepthQTBT )
#else
  if ( pps.getUseDQP() && uiDepth <= pps.getMaxCuDQPDepth() )
#endif
  {
#if JVET_C0024_DELTA_QP_FIX
    const UInt numValidComp = pcCU->getPic()->getNumberValidComponents();
    if( pcCU->getSlice()->isIntra() )
    {
      if ( ( pcCU->getTextType() == CHANNEL_TYPE_LUMA && pcCU->getCbf(0, COMPONENT_Y) )
        || ( pcCU->getTextType() == CHANNEL_TYPE_CHROMA && ( ( numValidComp>COMPONENT_Cb && pcCU->getCbf( 0, COMPONENT_Cb ) ) || ( numValidComp>COMPONENT_Cr && pcCU->getCbf( 0, COMPONENT_Cr) ) ) ) 
        )
      {
        m_pcEntropyCoder->resetBits();
        m_pcEntropyCoder->encodeQP( pcCU, 0, false );
        pcCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
        pcCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
        pcCU->getTotalCost() = m_pcRdCost->calcRdCost( pcCU->getTotalBits(), pcCU->getTotalDistortion() );
      }
      else
      {
        pcCU->setQPSubParts( pcCU->getRefQP( 0 ), 0, pcCU->getWidth(0), pcCU->getHeight(0) ); // set QP to default QP
      }
    }
    else
    {
#endif
    if ( pcCU->getQtRootCbf( 0) )
    {
      m_pcEntropyCoder->resetBits();
      m_pcEntropyCoder->encodeQP( pcCU, 0, false );
      pcCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
      pcCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      pcCU->getTotalCost() = m_pcRdCost->calcRdCost( pcCU->getTotalBits(), pcCU->getTotalDistortion() );
    }
    else
    {
#if JVET_C0024_DELTA_QP_FIX
      pcCU->setQPSubParts( pcCU->getRefQP( 0 ), 0, pcCU->getWidth(0), pcCU->getHeight(0) ); // set QP to default QP
#else
      pcCU->setQPSubParts( pcCU->getRefQP( 0 ), 0, uiDepth ); // set QP to default QP
#endif
    }
#if JVET_C0024_DELTA_QP_FIX
    }
#endif
  }
}

Void TEncCu::xCopyAMVPInfo (AMVPInfo* pSrc, AMVPInfo* pDst)
{
  pDst->iN = pSrc->iN;
  for (Int i = 0; i < pSrc->iN; i++)
  {
    pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
  }
}
#if JVET_C0024_QTBT
Void TEncCu::xCopyYuv2Pic(TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth, UInt uiWidth, UInt uiHeight
#if BARC
  ,Bool isBARCFlag
#endif
  )
#else
Void TEncCu::xCopyYuv2Pic(TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth )
#endif
{
  UInt uiAbsPartIdxInRaster = g_auiZscanToRaster[uiAbsPartIdx];
  UInt uiSrcBlkWidth = rpcPic->getNumPartInCtuWidth() >> (uiSrcDepth);
  UInt uiBlkWidth    = rpcPic->getNumPartInCtuWidth() >> (uiDepth);
  UInt uiPartIdxX = ( ( uiAbsPartIdxInRaster % rpcPic->getNumPartInCtuWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
  UInt uiPartIdxY = ( ( uiAbsPartIdxInRaster / rpcPic->getNumPartInCtuWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
  UInt uiPartIdx = uiPartIdxY * ( uiSrcBlkWidth / uiBlkWidth ) + uiPartIdxX;
#if JVET_C0024_QTBT
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
  if (isLuma(rpcPic->getCtu(uiCUAddr)->getSlice()->getTextType()))
  {
    m_pppcRecoYuvBest[uiWIdx][uiHIdx]->copyToPicYuv( rpcPic->getPicYuvRec (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);

    m_pppcPredYuvBest[uiWIdx][uiHIdx]->copyToPicYuv( rpcPic->getPicYuvPred (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
#if BARC
  if (isBARCFlag)
    m_pppcRecoYuvBest[uiWIdx][uiHIdx]->copyToPicYuv( rpcPic->getPicYuvDScaleRef(), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx, true);
#endif
  }
  else
  {
    assert(isChroma(rpcPic->getCtu(uiCUAddr)->getSlice()->getTextType()));
    for(Int comp=1; comp<m_pppcRecoYuvBest[uiWIdx][uiHIdx]->getNumberValidComponents(); comp++)
    {
      m_pppcRecoYuvBest[uiWIdx][uiHIdx]->copyToPicComponent  ( ComponentID(comp), rpcPic->getPicYuvRec (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx );
      m_pppcPredYuvBest[uiWIdx][uiHIdx]->copyToPicComponent  ( ComponentID(comp), rpcPic->getPicYuvPred (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx );
    }
  }
#else
  m_ppcRecoYuvBest[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvRec (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);

  m_ppcPredYuvBest[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvPred (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
#endif
}

#if JVET_C0024_QTBT //uiSplitMethod: 0: quadtree; 1: hor; 2: ver
#if TT_TREE_DS1
Void TEncCu::xCopyYuv2Tmp( UInt uiPartUnitIdxSrc, UInt uiPartUnitIdxDst, UInt uiWidth, UInt uiHeight, UInt uiSplitMethod )
#else
Void TEncCu::xCopyYuv2Tmp( UInt uiPartUnitIdx, UInt uiWidth, UInt uiHeight, UInt uiSplitMethod )
#endif
#else
Void TEncCu::xCopyYuv2Tmp( UInt uiPartUnitIdx, UInt uiNextDepth )
#endif
{
#if JVET_C0024_QTBT 
  UInt uiWIdx = g_aucConvertToBit[uiWidth];
  UInt uiHIdx = g_aucConvertToBit[uiHeight];
#if TT_TREE_DS1 //xCopyYuv2Tmp
  UInt uiNextWIdx,uiNextHIdx;
  if(uiSplitMethod<=2)
  {
    uiNextWIdx = uiWIdx - ((uiSplitMethod & 1)==0 ? 1: 0);
    uiNextHIdx = uiHIdx - ((uiSplitMethod & 2)==0 ? 1: 0);
  }
  else if(uiSplitMethod==3 || uiSplitMethod == 5)
  {
    uiNextWIdx = uiWIdx;
    uiNextHIdx = uiHIdx - (uiSplitMethod==3 ? 1: 2);
  }
  else if(uiSplitMethod==4 || uiSplitMethod == 6)
  {
    uiNextWIdx = uiWIdx - (uiSplitMethod==4 ? 1: 2);
    uiNextHIdx = uiHIdx;
  }
  
  UInt uiRasterSrc = g_auiZscanToRaster[uiPartUnitIdxSrc];
  UInt uiRasterDst = g_auiZscanToRaster[uiPartUnitIdxDst];
  assert(uiRasterSrc>=uiRasterDst);
  UInt uiStrideInPart = m_pppcBestCU[uiNextWIdx][uiNextHIdx]->getPic()->getNumPartInCtuWidth();
  UInt uiDeltaY       = (uiRasterSrc - uiRasterDst) / uiStrideInPart;
  UInt uiDeltaX       = (uiRasterSrc - uiRasterDst) % uiStrideInPart;
  UInt uiPartUnitIdx  = g_auiRasterToZscan[uiDeltaY * uiStrideInPart + uiDeltaX];
#else
  UInt uiNextWIdx = uiWIdx - ((uiSplitMethod & 1)==0 ? 1: 0);
  UInt uiNextHIdx = uiHIdx - ((uiSplitMethod & 2)==0 ? 1: 0);
#endif
  m_pppcRecoYuvBest[uiNextWIdx][uiNextHIdx]->copyToPartYuv( m_pppcRecoYuvTemp[uiWIdx][uiHIdx], uiPartUnitIdx );
  m_pppcPredYuvBest[uiNextWIdx][uiNextHIdx]->copyToPartYuv( m_pppcPredYuvTemp[uiWIdx][uiHIdx], uiPartUnitIdx );
#else
  UInt uiCurrDepth = uiNextDepth - 1;
  m_ppcRecoYuvBest[uiNextDepth]->copyToPartYuv(m_ppcRecoYuvTemp[uiCurrDepth], uiPartUnitIdx);
#if TRACE_YUV
  m_ppcPredYuvBest[uiNextDepth]->copyToPartYuv(m_ppcPredYuvTemp[uiCurrDepth], uiPartUnitIdx);
#else
  m_ppcPredYuvBest[uiNextDepth]->copyToPartYuv(m_ppcPredYuvBest[uiCurrDepth], uiPartUnitIdx);
#endif
#endif
}

/** Function for filling the PCM buffer of a CU using its original sample array
 * \param pCU pointer to current CU
 * \param pOrgYuv pointer to original sample array
 */
Void TEncCu::xFillPCMBuffer     ( TComDataCU* pCU, TComYuv* pOrgYuv )
{
  const ChromaFormat format = pCU->getPic()->getChromaFormat();
  const UInt numberValidComponents = getNumberValidComponents(format);
  for (UInt componentIndex = 0; componentIndex < numberValidComponents; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

    const UInt width  = pCU->getWidth(0)  >> getComponentScaleX(component, format);
    const UInt height = pCU->getHeight(0) >> getComponentScaleY(component, format);

    Pel *source      = pOrgYuv->getAddr(component, 0, width);
    Pel *destination = pCU->getPCMSample(component);

    const UInt sourceStride = pOrgYuv->getStride(component);

    for (Int line = 0; line < height; line++)
    {
      for (Int column = 0; column < width; column++)
      {
        destination[column] = source[column];
      }

      source      += sourceStride;
      destination += width;
    }
  }
}

#if ADAPTIVE_QP_SELECTION
/** Collect ARL statistics from one block
  */
Int TEncCu::xTuCollectARLStats(TCoeff* rpcCoeff, TCoeff* rpcArlCoeff, Int NumCoeffInCU, Double* cSum, UInt* numSamples )
{
  for( Int n = 0; n < NumCoeffInCU; n++ )
  {
    TCoeff u = abs( rpcCoeff[ n ] );
    TCoeff absc = rpcArlCoeff[ n ];

    if( u != 0 )
    {
      if( u < LEVEL_RANGE )
      {
        cSum[ u ] += ( Double )absc;
        numSamples[ u ]++;
      }
      else
      {
        cSum[ LEVEL_RANGE ] += ( Double )absc - ( Double )( u << ARL_C_PRECISION );
        numSamples[ LEVEL_RANGE ]++;
      }
    }
  }

  return 0;
}

//! Collect ARL statistics from one CTU
Void TEncCu::xCtuCollectARLStats(TComDataCU* pCtu )
{
  Double cSum[ LEVEL_RANGE + 1 ];     //: the sum of DCT coefficients corresponding to data type and quantization output
  UInt numSamples[ LEVEL_RANGE + 1 ]; //: the number of coefficients corresponding to data type and quantization output

  TCoeff* pCoeffY = pCtu->getCoeff(COMPONENT_Y);
  TCoeff* pArlCoeffY = pCtu->getArlCoeff(COMPONENT_Y);
  const TComSPS &sps = *(pCtu->getSlice()->getSPS());

#if JVET_C0024_QTBT
  const UInt uiMinCUWidth = sps.getCTUSize() >> sps.getMaxTotalCUDepth(); // NOTE: ed - this is not the minimum CU width. It is the square-root of the number of coefficients per part.
#else
  const UInt uiMinCUWidth = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth(); // NOTE: ed - this is not the minimum CU width. It is the square-root of the number of coefficients per part.
#endif
  const UInt uiMinNumCoeffInCU = 1 << uiMinCUWidth;                          // NOTE: ed - what is this?

  memset( cSum, 0, sizeof( Double )*(LEVEL_RANGE+1) );
  memset( numSamples, 0, sizeof( UInt )*(LEVEL_RANGE+1) );

  // Collect stats to cSum[][] and numSamples[][]
  for(Int i = 0; i < pCtu->getTotalNumPart(); i ++ )
  {
#if JVET_C0024_QTBT
    UInt uiTrIdx = 0;
#else
    UInt uiTrIdx = pCtu->getTransformIdx(i);
#endif

    if(pCtu->isInter(i) && pCtu->getCbf( i, COMPONENT_Y, uiTrIdx ) )
    {
      xTuCollectARLStats(pCoeffY, pArlCoeffY, uiMinNumCoeffInCU, cSum, numSamples);
    }//Note that only InterY is processed. QP rounding is based on InterY data only.

    pCoeffY  += uiMinNumCoeffInCU;
    pArlCoeffY  += uiMinNumCoeffInCU;
  }

  for(Int u=1; u<LEVEL_RANGE;u++)
  {
    m_pcTrQuant->getSliceSumC()[u] += cSum[ u ] ;
    m_pcTrQuant->getSliceNSamples()[u] += numSamples[ u ] ;
  }
  m_pcTrQuant->getSliceSumC()[LEVEL_RANGE] += cSum[ LEVEL_RANGE ] ;
  m_pcTrQuant->getSliceNSamples()[LEVEL_RANGE] += numSamples[ LEVEL_RANGE ] ;
}
#endif

#if COM16_C1016_AFFINE
Void TEncCu::xCheckRDCostAffineMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU )
{
  assert( rpcTempCU->getSlice()->getSliceType() != I_SLICE );
  if ( getFastDeltaQp() )
  {
    return;
  }

  TComMvField cAffineMvField[2][3];
  UChar uhInterDirNeighbours[1] = {0};
  Int numValidMergeCand = 0;
  const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);

  UChar uhDepth = rpcTempCU->getDepth( 0 );
#if JVET_C0024_QTBT
  if (rpcBestCU->getWidth(0) * rpcBestCU->getHeight(0) < 64)
  {
    return;
  }
  UInt uiWIdx = g_aucConvertToBit[rpcBestCU->getWidth(0)];
  UInt uiHIdx = g_aucConvertToBit[rpcBestCU->getHeight(0)];
#else
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  UInt neighborAffineParam = 0;
  rpcTempCU->getAffineMergeCandidates(0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand, neighborAffineParam);
  assert(neighborAffineParam == 0 || neighborAffineParam == 1);
#else
  rpcTempCU->getAffineMergeCandidates( 0, 0, cAffineMvField, uhInterDirNeighbours, numValidMergeCand );
#endif
  if ( numValidMergeCand == -1 )
  {
    return;
  }
  Int mergeCandBuffer[1] = {0};

  UInt iteration;
  if ( rpcTempCU->isLosslessCoded(0))
  {
    iteration = 1;
  }
  else
  {
    iteration = 2;
  }

  for ( UInt uiNoResidual = 0; uiNoResidual < iteration; ++uiNoResidual )
  {
    for ( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
    {
      if ( !( uiNoResidual==1 && mergeCandBuffer[uiMergeCand]==1 ) )
      {
        // set MC parameters
        rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth );
        rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag, 0, uhDepth );
        rpcTempCU->setChromaQpAdjSubParts( bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );
#if !JVET_C0024_QTBT
        rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth );
#endif
        rpcTempCU->setAffineFlagSubParts( true, 0, 0, uhDepth );
#if JVECT_C0062_AFFINE_SIX_PARAM
        rpcTempCU->setAffineParamFlagSubParts(neighborAffineParam, 0, 0, uhDepth);
#endif
        rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth );
        rpcTempCU->setMergeIndexSubParts( uiMergeCand, 0, 0, uhDepth );
        rpcTempCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth );
#if JVECT_C0062_AFFINE_SIX_PARAM
        Bool b6Param = neighborAffineParam;
        rpcTempCU->setAllAffineMvField(0, 0, cAffineMvField[0 + 2 * uiMergeCand], REF_PIC_LIST_0, 0, b6Param);
        rpcTempCU->setAllAffineMvField(0, 0, cAffineMvField[1 + 2 * uiMergeCand], REF_PIC_LIST_1, 0, b6Param);
#else
        rpcTempCU->setAllAffineMvField( 0, 0, cAffineMvField[0 + 2*uiMergeCand], REF_PIC_LIST_0, 0 );
        rpcTempCU->setAllAffineMvField( 0, 0, cAffineMvField[1 + 2*uiMergeCand], REF_PIC_LIST_1, 0 );
#endif

#if COM16_C806_OBMC
        rpcTempCU->setOBMCFlagSubParts( true, 0, uhDepth );
#endif

        // do MC
#if JVET_C0024_QTBT
        m_pcPredSearch->motionCompensation ( rpcTempCU, m_pppcPredYuvTemp[uiWIdx][uiHIdx] );

#if COM16_C806_OBMC
        m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_pppcPredYuvTemp[uiWIdx][uiHIdx], m_pppcTmpYuv1[uiWIdx][uiHIdx], m_pppcTmpYuv2[uiWIdx][uiHIdx] );
#endif
        // estimate residual and encode everything
        m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                   m_pppcOrigYuv    [uiWIdx][uiHIdx],
                                                   m_pppcPredYuvTemp[uiWIdx][uiHIdx],
                                                   m_pppcResiYuvTemp[uiWIdx][uiHIdx],
                                                   m_pppcResiYuvBest[uiWIdx][uiHIdx],
                                                   m_pppcRecoYuvTemp[uiWIdx][uiHIdx],
                                                   (uiNoResidual != 0) 
#if COM16_C806_EMT
                                                   , rpcBestCU->getTotalCost()
#endif
          );
#else
        m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[uhDepth] );

#if COM16_C806_OBMC
        m_pcPredSearch->subBlockOBMC( rpcTempCU, 0, m_ppcPredYuvTemp[uhDepth], m_ppcTmpYuv1[uhDepth], m_ppcTmpYuv2[uhDepth] );
#endif
        // estimate residual and encode everything
        m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                   m_ppcOrigYuv    [uhDepth],
                                                   m_ppcPredYuvTemp[uhDepth],
                                                   m_ppcResiYuvTemp[uhDepth],
                                                   m_ppcResiYuvBest[uhDepth],
                                                   m_ppcRecoYuvTemp[uhDepth],
                                                   (uiNoResidual != 0) 
#if COM16_C806_EMT
                                                   , rpcBestCU->getTotalCost()
#endif
          );
#endif
        if ( uiNoResidual == 0 && rpcTempCU->getQtRootCbf(0) == 0 )
        {
          // If no residual when allowing for one, then set mark to not try case where residual is forced to 0
          mergeCandBuffer[uiMergeCand] = 1;
        }

        rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth );
        Int orgQP = rpcTempCU->getQP( 0 );
        xCheckDQP( rpcTempCU );
        xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
#if BIF_NO_SKIP_ONLY
        if (rpcTempCU->getSIFFlag(0) && rpcTempCU->isSkipped(0))
        {
          rpcTempCU->getTotalCost() = DBL_MAX;
        }
#endif
        rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag );
      }
    }
  }
}
#endif
#if BARC_FUNC
Void TEncCu::reComputeBARCDistortion( TComDataCU* pCtu, ChannelType chType)
{
  UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();
  
  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();
  Distortion distortion = 0;
  const Int channelBitDepth = pCtu->getSlice()->getSPS()->getBitDepth(chType);
  UInt startIdx, endIdx;

  if (m_uiCUPelX + uiWidth > picWidth)
  {
    uiWidth = picWidth - m_uiCUPelX;
  }
  if (m_uiCUPelY + uiHeight > picHeight)
  {
    uiHeight = picHeight - m_uiCUPelY;
  }

  if (chType == CHANNEL_TYPE_LUMA)
  {
    startIdx = 0;
    endIdx = 1;
  }
  else
  {
    startIdx = 1;
    endIdx = 3;
  }
  for (UInt comp=startIdx; comp<endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    if (comp == 1)
    {
      uiWidth >>= 1;
      uiHeight >>= 1;
    }

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    Pel* pSrcOrg = pCtu->getPic()->getPicYuvOrg()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);

    distortion += m_pcRdCost->getDistPart(channelBitDepth, pSrc, iSrcStride, pSrcOrg, iSrcStride, uiWidth, uiHeight, compID);

    if (comp == endIdx -1 )
    {
      pCtu->getTotalCost() = m_pcRdCost->calcRdCost( pCtu->getTotalBits(), distortion);
      pCtu->getTotalDistortion() = distortion;
    }
  }
}
#endif
#if COMPACT_RESOLUTION_CNN
Void TEncCu::upSampleRecDeep(TComDataCU*& pCtu, ChannelType chType)
{
	Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false : true;
	Bool isLeftCUAvail = pCtu->getCUPelX() == 0 ? false : true;

	UInt m_uiCUPelX = pCtu->getCUPelX();
	UInt m_uiCUPelY = pCtu->getCUPelY();
	UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
	UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

	UInt upWidth = pCtu->getSlice()->getSPS()->getCTUSize();
	UInt upHeight = pCtu->getSlice()->getSPS()->getCTUSize();
	UInt uiWidth = upWidth >> BARC_DScaleRatio;
	UInt uiHeight = upHeight >> BARC_DScaleRatio;
	UInt starIdx, endIdx;

	if (m_uiCUPelX + upWidth > picWidth)
	{
		upWidth = picWidth - m_uiCUPelX;
		uiWidth = upWidth >> BARC_DScaleRatio;
	}
	if (m_uiCUPelY + upHeight > picHeight)
	{
		upHeight = picHeight - m_uiCUPelY;
		uiHeight = upHeight >> BARC_DScaleRatio;
	}

	if (chType == CHANNEL_TYPE_LUMA)
	{
		starIdx = 0;
		endIdx = 1;
	}
	else
	{
		starIdx = 1;
		endIdx = 3;
		upWidth >>= 1;
		upHeight >>= 1;
		uiWidth >>= 1;
		uiHeight >>= 1;
	}
	for (UInt comp = starIdx; comp < endIdx; comp++)
	{
		const ComponentID compID = ComponentID(comp);

		Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
		Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pCtu->getCtuRsAddr(), 0);
		UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);

		for (Int i = 0; i < uiHeight; i++)
			memcpy(pRefScale + i*iSrcStride, pSrc + i*iSrcStride, uiWidth*sizeof(Pel));

		pCtu->upScaleComponent2x2(pRefScale, pSrc, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, compID);
	}

#if BARC_REF_UP_SAMPLE_CNN
	if (chType == CHANNEL_TYPE_LUMA)
	{
		Double* pDstNet = new Double[MAX_CU_SIZE*MAX_CU_SIZE];
		Pel*       pDst = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];

		Pel* pSrc = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);
		Pel* pRec = pCtu->getPic()->getPicYuvRec()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);

		UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
		UInt iDstStride = MAX_CU_SIZE;

		pCtu->getPic()->CNNUpScaleLumaComponent2x2(true, pSrc, pDstNet, iSrcStride, iDstStride, uiWidth, uiHeight);

		Double* tpDstNet = pDstNet;
		Pel*    tpSrcRDO = pRec;
		Pel*    tpRec = pDst;
		Pel*    pSrcOrgRDO = pCtu->getPic()->getPicYuvOrg()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);

		for (Int i = 0; i < upHeight; i++)
		{
			for (Int j = 0; j < upWidth; j++)
			{
				tpRec[j] = ClipA(int(tpDstNet[j] * PixelRangeInt + 0.5), COMPONENT_Y);
			}
			tpSrcRDO += iSrcStride;
			tpRec += iDstStride;
			tpDstNet += iDstStride;
		}
		Double distortionYRDONet = getDist(pDst, iDstStride, pSrcOrgRDO, iSrcStride, upWidth, upHeight);
		Double distortionYRDO = getDist(pRec, iSrcStride, pSrcOrgRDO, iSrcStride, upWidth, upHeight);
		{
			pCtu->setUpCNNFlag(true);
			for (Int i = 0; i < upHeight; i++)
				memcpy(pRec + i*iSrcStride, pDst + i*iDstStride, upWidth*sizeof(Pel));
		}

		if (pDst)
		{
			delete[] pDst;
			pDst = NULL;
		}
		if (pDstNet)
		{
			delete[] pDstNet;
			pDstNet = NULL;
		}
	}
#endif	
}
#endif

#if BARC_FUNC
Void TEncCu::upSampleRec ( TComDataCU*& pCtu, ChannelType chType)
{
  Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false : true;
  Bool isLeftCUAvail = pCtu->getCUPelX() == 0 ? false : true;

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

  UInt upWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt upHeight = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiWidth = upWidth >> BARC_DScaleRatio;
  UInt uiHeight = upHeight >> BARC_DScaleRatio;
  UInt starIdx, endIdx;

  if (m_uiCUPelX + upWidth > picWidth)
  {
    upWidth = picWidth - m_uiCUPelX;
    uiWidth = upWidth >> BARC_DScaleRatio;
  }
  if (m_uiCUPelY + upHeight > picHeight)
  {
    upHeight = picHeight - m_uiCUPelY;
    uiHeight = upHeight >> BARC_DScaleRatio;
  }

  if (chType == CHANNEL_TYPE_LUMA)
  {
    starIdx = 0;
    endIdx = 1;
  }
  else
  {
    starIdx = 1;
    endIdx = 3;
    upWidth >>= 1;
    upHeight >>= 1;
    uiWidth >>= 1;
    uiHeight >>= 1;
  }
  for (UInt comp=starIdx; comp<endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID,pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);

#if COMPACT_RESOLUTION_CNN
	Pel* pDScaleRec = m_DScaleBufferBicubic;
	UInt iDstStride = MAX_CU_SIZE >> BARC_DScaleRatio;
	// copy rec to DScaleBufferBicubic
	for (Int i = 0; i < uiHeight; i++)
		memcpy(pDScaleRec + i*iDstStride, pSrc + i*iSrcStride, uiWidth*sizeof(Pel));
#endif
    for (Int i = 0; i < uiHeight; i++)
      memcpy(pRefScale+i*iSrcStride, pSrc+i*iSrcStride, uiWidth*sizeof(Pel));

    pCtu->upScaleComponent2x2( pRefScale, pSrc, iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail,isLeftCUAvail, compID);
  }
#if BARC_REF_UP_SAMPLE_CNN
#if COMPACT_RESOLUTION_CNN
  if (chType == CHANNEL_TYPE_LUMA)
  {
	  pCtu->setUpCNNFlag(false);
  }
#else
  if (chType = CHANNEL_TYPE_LUMA)
  {
	  Double*  pCNNRes = new Double[MAX_CU_SIZE*MAX_CU_SIZE];  // The reconstructed residual from CNN output
	  Pel*     pLRRec = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(COMPONENT_Y, pCtu->getCtuRsAddr(), 0); //Low resolution rec
	  Pel*     pSRRec = pCtu->getPic()->getPicYuvRec()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);             //Super resolution rec from DCTIF
	  UInt     iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
	  UInt     iDstStride = MAX_CU_SIZE;
	  pCtu->getPic()->CNNUpScaleLumaComponent2x2(true, pLRRec, pCNNRes, iSrcStride, iDstStride, uiWidth, uiHeight);

	  Pel*     pSrcOrg = pCtu->getPic()->getPicYuvOrg()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);  //The original pixel value
	  Pel*     pCNNRec = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];                                               // CNNRec = CNNRes + SRRec 
	  Pel*     pTempCNNRec = pCNNRec;
	  Pel*     pTempSRRec = pSRRec;
	  Double*  pTempCNNRes = pCNNRes;
	  for (Int i = 0; i < upHeight; i++)
	  {
		  for (Int j = 0; j < upWidth; j++)
		  {
			  pTempCNNRec[j] = ClipA(int(pTempSRRec[j] + pTempCNNRes[j] * PixelRangeInt + 0.5), COMPONENT_Y);
		  }
		  pTempCNNRec += iSrcStride;
		  pTempSRRec += iSrcStride;
		  pTempCNNRes += iDstStride;
	  }
	  // RDO to determine best Up-sample mode
#if BARC_UP_SAMPLE_RDO
	  Double filtLumaDistortion = getDist(pSRRec, iSrcStride, pSrcOrg, iSrcStride, upWidth, upHeight);
	  Double CNNLumaDistortion = getDist(pCNNRec, iDstStride, pSrcOrg, iSrcStride, upWidth, upHeight);

	  if (CNNLumaDistortion < filtLumaDistortion)
	  {
		  pCtu->setUpCNNFlag(true);
		  for (Int i = 0; i < upHeight; i++)
			  memcpy(pSRRec + i*iSrcStride, pCNNRec + iSrcStride, upWidth*sizeof(Pel));
	  }
	  else
	  {
		  pCtu->setUpCNNFlag(false);
	  }
#else
	  for (Int i = 0; i < upHeight; i++)
		  memcpy(pSRRec + i*iSrcStride, pCNNRec + iSrcStride, upWidth*sizeof(Pel));
#endif
	  if (pCNNRes)
	  {
		  delete[] pCNNRes;
		  pCNNRes = NULL;
	  }
	  if (pCNNRec)
	  {
		  delete[] pCNNRec;
		  pCNNRec = NULL;
	  }
  }
#endif
#endif
}
#endif

#if BARC_FUNC
UInt TEncCu::getDist(Pel* piCur, Int iCurStride, Pel* piOrg, Int iOrgStride, UInt uiBlkWidth, UInt uiBlkHeight, UInt bitDepth)
{
	UInt dist = 0;
	UInt uiShift = DISTORTION_PRECISION_ADJUSTMENT((bitDepth - 8) << 1);
	Intermediate_Int  iTemp;

	for (Int i = 0; i < uiBlkHeight; i++)
	{
		for (Int j = 0; j < uiBlkWidth; j++)
		{
			iTemp = (piOrg[i*iOrgStride + j] - piCur[i*iCurStride + j]);
			dist += Distortion((iTemp * iTemp) >> uiShift);
		}
	}
	return dist;
}
#endif

#if BARC_FUNC
Void TEncCu::downSampleRec ( TComDataCU*& pCtu, ChannelType chType)
{
  UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

  Bool isAboveCUAvail = pCtu->getCUPelY() == 0 ? false : true;
  Bool isLeftCUAvail = pCtu->getCUPelX() == 0 ? false : true;

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();
  Double distortion = 0;
  UInt startIdx, endIdx;		

  if (m_uiCUPelX + uiWidth > picWidth)
  {
    uiWidth = picWidth - m_uiCUPelX;
  }
  if (m_uiCUPelY + uiHeight > picHeight)
  {
    uiHeight = picHeight - m_uiCUPelY;
  }

  if (chType == CHANNEL_TYPE_LUMA)
  {
    startIdx = 0;
    endIdx = 1;
  }
  else
  {
    startIdx = 1;
    endIdx = 3;
    uiWidth >>= 1;
    uiHeight >>= 1;
  }

  for (UInt comp=startIdx; comp<endIdx; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pSrc = pCtu->getPic()->getPicYuvRec()->getAddr(compID,pCtu->getCtuRsAddr(), 0);
    Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID,pCtu->getCtuRsAddr(), 0);
    UInt iSrcStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
#if COMPACT_RESOLUTION_CNN
	pCtu->getPic()->CNNDownScaleLumaComponent2x2(false, pSrc, pRefScale, iSrcStride, iSrcStride, uiWidth, uiHeight);
#else
    pCtu->downScaleComponent2x2ZeroPhaseBoundary(pSrc, pRefScale,iSrcStride, iSrcStride, uiWidth, uiHeight, isAboveCUAvail, isLeftCUAvail, compID);
#endif
  }
}
#endif

#if BARC_FUNC
Void TEncCu::initEncSliceBARC(ChannelType chType, TComPic* pcPic, Int iGOPid, TComSlice*& rpcSlice, Bool isBARC )
{
//  TEncSlice* rEncSlice;
  rpcSlice = pcPic->getSlice(0);

  Double dQP,dLambda;
  dQP = rpcSlice->getSliceQp();

  if (isBARC && chType == CHANNEL_TYPE_LUMA) dQP -= BARC_AdpQPSTEP;

  Int SHIFT_QP = 12;
  Double qp_temp = (Double) dQP  - SHIFT_QP;
  Double dQPFactor = 0.57;

  dLambda = dQPFactor*pow( 2.0, qp_temp/3.0 );
  dQP = max( -rpcSlice->getSPS()->getQpBDOffset(chType), min( MAX_QP, (Int) floor( dQP + 0.5 ) ) );

  getSliceEncoder()->setUpLambda(rpcSlice, dLambda, (Int)dQP);
#if WCG_LUMA_DQP_CM_SCALE
  // cost = Distortion + Lambda*R, 
  // when QP is adjusted by luma, distortion is changed, so we have to adjust lambda to match the distortion, then the cost function becomes
  // costA = Distortion + AdjustedLambda * R          -- currently, costA is still used when calculating intermediate cost of using SAD, HAD, resisual etc. 
  // an alternative way is to weight the distortion to before the luma QP adjustment, then the cost function becomes
  // costB = weightedDistortion + Lambda * R          -- currently, costB is used to calculat final cost, and when DF_FUNC is DF_DEFAULT 
  m_pcRdCost->saveUnadjustedLambda();
#endif
}
#endif
#if BARC_FUNC
Void TEncCu::copyLumaDscaletoRec(TComDataCU* pCtu, Bool isCopyDscaletoRec)
{
  UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

  if (m_uiCUPelX + uiWidth > picWidth)
  {
    uiWidth = picWidth - m_uiCUPelX;
  }
  if (m_uiCUPelY + uiHeight > picHeight)
  {
    uiHeight = picHeight - m_uiCUPelY;
  }

  for (UInt comp = 0; comp < 1; comp++)
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pRec = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
    Pel* pRecScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(compID, pCtu->getCtuRsAddr(), 0);
#if BARC_ONLY_FOR_INTRA
	Pel* pTmp = m_reconLumaBufferHM;
#else
    Pel* pTmp = m_reconLumaBuffer;
#endif
    UInt iRecStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
    UInt iTmpStride = MAX_CU_SIZE;

    if (isCopyDscaletoRec)
    {
      for (Int y = 0; y < uiHeight; y++)
        memcpy(pTmp + y*iTmpStride, pRec + y*iRecStride, sizeof(Pel)*uiWidth);
      for (Int y = 0; y < uiHeight / 2; y++)
        memcpy(pRec + y*iRecStride, pRecScale + y*iRecStride, sizeof(Pel)*uiWidth / 2);
    }
    else
    {
      for (Int y = 0; y < uiHeight; y++)
        memcpy(pRec + y*iRecStride, pTmp + y*iTmpStride, sizeof(Pel)*uiWidth);
    }
  }
}
#endif
#if BARC_FUNC
Void TEncCu::copyRecToTmp(TComDataCU* pCtu, Bool isCopyRecToTmp)
{
  UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
  UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

  UInt m_uiCUPelX = pCtu->getCUPelX();
  UInt m_uiCUPelY = pCtu->getCUPelY();
  UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

  if (m_uiCUPelX + uiWidth > picWidth)
  {
    uiWidth = picWidth - m_uiCUPelX;
  }
  if (m_uiCUPelY + uiHeight > picHeight)
  {
    uiHeight = picHeight - m_uiCUPelY;
  }
#if BARC_ONLY_FOR_LUMA
  for (UInt comp = 0; comp < 1; comp++)
#else
  for (UInt comp = 0; comp < 3; comp++)
#endif
  {
    const ComponentID compID = ComponentID(comp);

    Pel* pRec     = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
#if BARC_ONLY_FOR_LUMA
	Pel* pTmp = m_reconLumaBufferHM;
#else
    Pel* pTmp     = m_reconBuffer[comp];
#endif
    UInt iRecStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
    UInt iTmpStride = MAX_CU_SIZE;

    if (comp == 1)
    {
      uiHeight >>= 1;
      uiWidth >>= 1;
    }
    if (isCopyRecToTmp)
    {
      for (Int y = 0; y < uiHeight; y++)
        memcpy(pTmp + y*iTmpStride, pRec + y*iRecStride, sizeof(Pel)*uiWidth);
    }
    else
    {
      for (Int y = 0; y < uiHeight; y++)
        memcpy(pRec + y*iRecStride, pTmp + y*iTmpStride, sizeof(Pel)*uiWidth);
    }
  }
}
#endif
#if COMPACT_RESOLUTION_CNN
Void TEncCu::copyBicubicRecToTmp(TComDataCU* pCtu, Bool isCopyRecToTmp)
{
	UInt uiWidth = pCtu->getSlice()->getSPS()->getCTUSize();
	UInt uiHeight = pCtu->getSlice()->getSPS()->getCTUSize();

	UInt m_uiCUPelX = pCtu->getCUPelX();
	UInt m_uiCUPelY = pCtu->getCUPelY();
	UInt picWidth = pCtu->getSlice()->getSPS()->getPicWidthInLumaSamples();
	UInt picHeight = pCtu->getSlice()->getSPS()->getPicHeightInLumaSamples();

	if (m_uiCUPelX + uiWidth > picWidth)
	{
		uiWidth = picWidth - m_uiCUPelX;
	}
	if (m_uiCUPelY + uiHeight > picHeight)
	{
		uiHeight = picHeight - m_uiCUPelY;
	}

#if BARC_ONLY_FOR_LUMA
	for (UInt comp = 0; comp < 1; comp++)
#else
	for (UInt comp = 0; comp < 3; comp++)
#endif
	{
		const ComponentID compID = ComponentID(comp);

		Pel* pRec = pCtu->getPic()->getPicYuvRec()->getAddr(compID, pCtu->getCtuRsAddr(), 0);
		Pel* pTmp = m_reconBufferBicubic;
		UInt iRecStride = pCtu->getPic()->getPicYuvRec()->getStride(compID);
		UInt iTmpStride = MAX_CU_SIZE;

		if (comp == 1)
		{
			uiHeight >>= 1;
			uiWidth >>= 1;
		}
		if (isCopyRecToTmp)
		{
			for (Int y = 0; y < uiHeight; y++)
				memcpy(pTmp + y*iTmpStride, pRec + y*iRecStride, sizeof(Pel)*uiWidth);
		}
		else
		{
			for (Int y = 0; y < uiHeight; y++)
				memcpy(pRec + y*iRecStride, pTmp + y*iTmpStride, sizeof(Pel)*uiWidth);
		}
	}
	if (!isCopyRecToTmp)
	{
		Pel* pRefScale = pCtu->getPic()->getPicYuvDScaleRef()->getAddrDscale(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);
		UInt iDstStride = pCtu->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
		Pel* pSrc = m_DScaleBufferBicubic;
		UInt iSrcStride = MAX_CU_SIZE >> BARC_DScaleRatio;
		UInt DScaleWidth = uiWidth >> BARC_DScaleRatio;
		UInt DScaleHeidht = uiHeight >> BARC_DScaleRatio;

		for (Int i = 0; i < DScaleHeidht; i++)
			memcpy(pRefScale + i*iDstStride, pSrc + i*iSrcStride, DScaleWidth*sizeof(Pel));
	}
}
#endif


//! \}