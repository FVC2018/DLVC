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

/** \file     TComTrQuant.cpp
    \brief    transform and quantization class
*/

#include <stdlib.h>
#include <math.h>
#include <limits>
#include <memory.h>
#include "TComTrQuant.h"
#include "TComPic.h"
#include "ContextTables.h"
#include "TComTU.h"
#include "Debug.h"

#if VCEG_AZ08_KLT_COMMON
#include <algorithm>
#endif

#if VCEG_AZ08_USE_SSE_SPEEDUP
#include <intrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <smmintrin.h>
typedef unsigned short U16;
typedef unsigned int UInt;
typedef short I16;
typedef int   Int;

float InnerProduct_SSE_FLOATXSHORT(float *pa, short *pb, int m);
int InnerProduct_SSE_SHORT(short *pa, short *pb, int m);
void scaleMatrix(float **ppx, short **ppout, float scale, int rows, int cols);
void scaleMatrix(float **ppx, short **ppout, float scale, int rows, int cols);
int AbsSumOfVector(short *pa, short *pb, int m);
int AbsSumOfVectorLesseqthan8(short *pa, short *pb, int m);

UInt GetSAD4x4_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD);
UInt GetSAD8x8_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD);
UInt GetSAD16x16_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD);
UInt GetSAD32x32_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD);
#endif

#if VCEG_AZ08_KLT_COMMON //only support 4x4-32x32 now
UInt g_uiDepth2MaxCandiNum[5] = { MAX_CANDI_NUM, MAX_CANDI_NUM, MAX_CANDI_NUM, MAX_CANDI_NUM, MAX_CANDI_NUM };
UInt g_uiDepth2MinCandiNum[5] = { 8, 8, 8, 8, 8 };
Int g_maxValueThre = MAX_INT;

UInt g_uiDepth2Width[5] = { 4, 8, 16, 32, 64 };
//template size ===========
#if VCEG_AZ08_INTER_KLT
UInt g_uiDepth2InterTempSize[5] = { 3, 3, 3, 3, 3 };
#endif
#if VCEG_AZ08_INTRA_KLT
UInt g_uiDepth2IntraTempSize[5] = { 3, 3, 3, 3, 3 };
#endif

#define USE_EIGLIBDOUBLE                    1 ///<(default 1) If defined, will use double type for deriving eigen vectors
#define USE_FLOAT_COV                       1 ///<(default 1) If defined, will use the float rather than double for covariance
#define KLT_MODE                        65533 ///< Mark the mode as KLT mode

#include <iostream>

#ifdef __GNUC__
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wshadow"
#if __GNUC__>=6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include "../../extlib/Eigen/Dense"
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic pop
#else
#pragma GCC diagnostic warning "-Wshadow"
#endif
#if __GNUC__>=6
#pragma GCC diagnostic warning "-Wignored-attributes"
#pragma GCC diagnostic warning "-Wmisleading-indentation"
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#endif
#else
#include "../../extlib/Eigen/Dense"
#endif

using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;

#if USE_FLOAT_COV
#if USE_EIGLIBDOUBLE
typedef MatrixXd matrixTypeDefined;
typedef VectorXd vectorType;
#else
typedef MatrixXf matrixTypeDefined; //MatrixXd using double will be a little better -4.6% vs. -4.3% (bitrate) for classD with time complexity increased 7% due to the derivation of eigenvectors
typedef VectorXf vectorType; //VectorXd
#endif
#else
typedef MatrixXd matrixTypeDefined; //MatrixXd
typedef VectorXd vectorType; //VectorXd
#endif
void xKLTr(Int bitDepth, TCoeff *block, TCoeff *coeff, UInt uiTrSize);
void xIKLTr(Int bitDepth, TCoeff *coeff, TCoeff *block, UInt uiTrSize);
#endif

typedef struct
{
  Int    iNNZbeforePos0;
  Double d64CodedLevelandDist; // distortion and level cost only
  Double d64UncodedDist;    // all zero coded block distortion
  Double d64SigCost;
  Double d64SigCost_0;
} coeffGroupRDStats;

#if COM16_C806_EMT
#if JVET_C0024_QTBT
FwdTrans *fastFwdTrans[16][7] = 
#else
FwdTrans *fastFwdTrans[16][5] = 
#endif
{
#if JVET_C0024_QTBT
  {fastForwardDCT2_B2, fastForwardDCT2_B4, fastForwardDCT2_B8, fastForwardDCT2_B16, fastForwardDCT2_B32, fastForwardDCT2_B64, fastForwardDCT2_B128},
  {NULL,               fastForwardDCT5_B4, fastForwardDCT5_B8, fastForwardDCT5_B16, fastForwardDCT5_B32, fastForwardDCT5_B64, fastForwardDCT5_B128},
  {NULL,               fastForwardDCT8_B4, fastForwardDCT8_B8, fastForwardDCT8_B16, fastForwardDCT8_B32, fastForwardDCT8_B64, fastForwardDCT8_B128},
  {NULL,               fastForwardDST1_B4, fastForwardDST1_B8, fastForwardDST1_B16, fastForwardDST1_B32, fastForwardDST1_B64, fastForwardDST1_B128},
  {NULL,               fastForwardDST7_B4, fastForwardDST7_B8, fastForwardDST7_B16, fastForwardDST7_B32, fastForwardDST7_B64, fastForwardDST7_B128},
#else
  {fastForwardDCT2_B4, fastForwardDCT2_B8, fastForwardDCT2_B16, fastForwardDCT2_B32, fastForwardDCT2_B64},
  {fastForwardDCT5_B4, fastForwardDCT5_B8, fastForwardDCT5_B16, fastForwardDCT5_B32, NULL               },
  {fastForwardDCT8_B4, fastForwardDCT8_B8, fastForwardDCT8_B16, fastForwardDCT8_B32, NULL               },
  {fastForwardDST1_B4, fastForwardDST1_B8, fastForwardDST1_B16, fastForwardDST1_B32, NULL               },
  {fastForwardDST7_B4, fastForwardDST7_B8, fastForwardDST7_B16, fastForwardDST7_B32, NULL               },
#endif
};

#if JVET_C0024_QTBT
InvTrans *fastInvTrans[16][7] = 
{
  {fastInverseDCT2_B2, fastInverseDCT2_B4, fastInverseDCT2_B8, fastInverseDCT2_B16, fastInverseDCT2_B32, fastInverseDCT2_B64, fastInverseDCT2_B128},
  {NULL,               fastInverseDCT5_B4, fastInverseDCT5_B8, fastInverseDCT5_B16, fastInverseDCT5_B32, fastInverseDCT5_B64, fastInverseDCT5_B128},
  {NULL,               fastInverseDCT8_B4, fastInverseDCT8_B8, fastInverseDCT8_B16, fastInverseDCT8_B32, fastInverseDCT8_B64, fastInverseDCT8_B128},
  {NULL,               fastInverseDST1_B4, fastInverseDST1_B8, fastInverseDST1_B16, fastInverseDST1_B32, fastInverseDST1_B64, fastInverseDST1_B128},
  {NULL,               fastInverseDST7_B4, fastInverseDST7_B8, fastInverseDST7_B16, fastInverseDST7_B32, fastInverseDST7_B64, fastInverseDST7_B128},
};
#else
InvTrans *fastInvTrans[16][5] = 
{
  {fastInverseDCT2_B4, fastInverseDCT2_B8, fastInverseDCT2_B16, fastInverseDCT2_B32, fastInverseDCT2_B64},
  {fastInverseDCT5_B4, fastInverseDCT5_B8, fastInverseDCT5_B16, fastInverseDCT5_B32, NULL               },
  {fastInverseDCT8_B4, fastInverseDCT8_B8, fastInverseDCT8_B16, fastInverseDCT8_B32, NULL               },
  {fastInverseDST1_B4, fastInverseDST1_B8, fastInverseDST1_B16, fastInverseDST1_B32, NULL               },
  {fastInverseDST7_B4, fastInverseDST7_B8, fastInverseDST7_B16, fastInverseDST7_B32, NULL               },
};
#endif
#endif

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constants
// ====================================================================================================================

#define RDOQ_CHROMA                 1           ///< use of RDOQ in chroma


// ====================================================================================================================
// QpParam constructor
// ====================================================================================================================

QpParam::QpParam(const Int           qpy,
                 const ChannelType   chType,
                 const Int           qpBdOffset,
                 const Int           chromaQPOffset,
                 const ChromaFormat  chFmt )
{
  Int baseQp;

  if(isLuma(chType))
  {
    baseQp = qpy + qpBdOffset;
  }
  else
  {
    baseQp = Clip3( -qpBdOffset, (chromaQPMappingTableSize - 1), qpy + chromaQPOffset );

    if(baseQp < 0)
    {
      baseQp = baseQp + qpBdOffset;
    }
    else
    {
      baseQp = getScaledChromaQP(baseQp, chFmt) + qpBdOffset;
    }
  }

  Qp =baseQp;
  per=baseQp/6;
  rem=baseQp%6;
}

QpParam::QpParam(const TComDataCU &cu, const ComponentID compID)
{
  Int chromaQpOffset = 0;

  if (isChroma(compID))
  {
    chromaQpOffset += cu.getSlice()->getPPS()->getQpOffset(compID);
    chromaQpOffset += cu.getSlice()->getSliceChromaQpDelta(compID);

    chromaQpOffset += cu.getSlice()->getPPS()->getPpsRangeExtension().getChromaQpOffsetListEntry(cu.getChromaQpAdj(0)).u.offset[Int(compID)-1];
  }

  *this = QpParam(cu.getQP( 0 ),
                  toChannelType(compID),
                  cu.getSlice()->getSPS()->getQpBDOffset(toChannelType(compID)),
                  chromaQpOffset,
                  cu.getPic()->getChromaFormat());
}


// ====================================================================================================================
// TComTrQuant class member functions
// ====================================================================================================================

TComTrQuant::TComTrQuant()
{
  // allocate temporary buffers
  m_plTempCoeff  = new TCoeff[ MAX_CU_SIZE*MAX_CU_SIZE ];

#if JVET_C0024_QTBT && ITH_C1_RDOQ_MEM_OPT
  m_piRateIncUp = new  Int[MAX_TU_SIZE * MAX_TU_SIZE];
  m_piRateIncDown = new  Int[MAX_TU_SIZE * MAX_TU_SIZE];
  m_piSigRateDelta = new  Int[MAX_TU_SIZE * MAX_TU_SIZE];
  m_pdCostCoeff = new  Double[MAX_TU_SIZE * MAX_TU_SIZE];
  m_pdCostSig = new  Double[MAX_TU_SIZE * MAX_TU_SIZE];
  m_pdCostCoeff0 = new  Double[MAX_TU_SIZE * MAX_TU_SIZE];
  m_pcDeltaU = new  TCoeff[MAX_TU_SIZE * MAX_TU_SIZE];
#endif

  // allocate bit estimation class  (for RDOQ)
  m_pcEstBitsSbac = new estBitsSbacStruct;
  initScalingList();
#if PVC_AQS
  m_iEstQSFactor = 1 << ESM_SHIFT;
#endif
#if VCEG_AZ08_KLT_COMMON
  memset( m_pData , 0 , sizeof( m_pData ) );
  m_pppTarPatch = NULL;
  m_pCovMatrix = NULL;
  m_pppdTmpEigenVector = NULL;
  m_pppdEigenVector = NULL;
  m_pppsEigenVector = NULL;
#endif
}

TComTrQuant::~TComTrQuant()
{
  // delete temporary buffers
  if ( m_plTempCoeff )
  {
    delete [] m_plTempCoeff;
    m_plTempCoeff = NULL;
  }

#if ITH_C1_RDOQ_MEM_OPT
  if (m_piRateIncUp)
  {
	  delete[] m_piRateIncUp;
	  m_piRateIncUp = NULL;
  }

  if (m_piRateIncDown)
  {
	  delete[] m_piRateIncDown;
	  m_piRateIncDown = NULL;
  }

  if (m_piSigRateDelta)
  {
	  delete[] m_piSigRateDelta;
	  m_piSigRateDelta = NULL;
  }

  if (m_pdCostCoeff)
  {
	  delete[] m_pdCostCoeff;
	  m_pdCostCoeff = NULL;
  }

  if (m_pdCostSig)
  {
	  delete[] m_pdCostSig;
	  m_pdCostSig = NULL;
  }

  if (m_pdCostCoeff0)
  {
	  delete[] m_pdCostCoeff0;
	  m_pdCostCoeff0 = NULL;
  }

  if (m_pcDeltaU)
  {
	  delete[] m_pcDeltaU;
	  m_pcDeltaU = NULL;
  }
#endif

  // delete bit estimation class
  if ( m_pcEstBitsSbac )
  {
    delete m_pcEstBitsSbac;
  }
  destroyScalingList();
#if VCEG_AZ08_USE_KLT
  if (m_useKLT)
  {
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
      Int useIntraKLT = m_useKLT & 1;
      Int useInterKLT = (m_useKLT >> 1) & 1;
#endif
#endif
#if VCEG_AZ08_KLT_COMMON
      for (UInt i = 0; i < MAX_CANDI_NUM; i++)
      {
        if( m_pData[i] != NULL )
        {
          delete[] m_pData[i];
          m_pData[i] = NULL;
}
      }

#if VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY
      for (UInt i = 0; i < MAX_1DTRANS_LEN; i++)
      {
          delete[]m_pDataT[i];
          m_pDataT[i] = NULL;
      }
#endif
      if (m_pppTarPatch != NULL)
      {
          for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
          {
              UInt blkSize = g_uiDepth2Width[uiDepth];
#if VCEG_AZ08_USE_KLT
              UInt tempSize = 0;
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
              if(useIntraKLT && useInterKLT)
              {
                  tempSize = max(g_uiDepth2IntraTempSize[uiDepth], g_uiDepth2InterTempSize[uiDepth]);
              }
              else if (useIntraKLT && !useInterKLT)
              {
                  tempSize = g_uiDepth2IntraTempSize[uiDepth];
              }
              else if (!useIntraKLT && useInterKLT)
              {
                  tempSize = g_uiDepth2InterTempSize[uiDepth];
              }
#endif
#if VCEG_AZ08_INTRA_KLT && !VCEG_AZ08_INTER_KLT
              tempSize = g_uiDepth2IntraTempSize[uiDepth];
#endif
#if !VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
              tempSize = g_uiDepth2InterTempSize[uiDepth];
#endif
#else
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
              UInt tempSize = max(g_uiDepth2IntraTempSize[uiDepth], g_uiDepth2InterTempSize[uiDepth]);
#endif
#if VCEG_AZ08_INTRA_KLT && !VCEG_AZ08_INTER_KLT
              UInt tempSize = g_uiDepth2IntraTempSize[uiDepth];
#endif
#if !VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
              UInt tempSize = g_uiDepth2InterTempSize[uiDepth];
#endif
#endif
              UInt patchSize = blkSize + tempSize;
              for (UInt uiRow = 0; uiRow < patchSize; uiRow++)
              {
                  if (m_pppTarPatch[uiDepth][uiRow] != NULL)
                  {
                      delete[]m_pppTarPatch[uiDepth][uiRow]; m_pppTarPatch[uiDepth][uiRow] = NULL;
                  }
              }
              if (m_pppTarPatch[uiDepth] != NULL)
              {
                  delete[]m_pppTarPatch[uiDepth]; m_pppTarPatch[uiDepth] = NULL;
              }
          }
          delete [] m_pppTarPatch;
          m_pppTarPatch = NULL;
      }

      if (m_pCovMatrix != NULL)
      {
          for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
          {
              if (m_pCovMatrix[uiDepth] != NULL)
              {
                  delete[]m_pCovMatrix[uiDepth]; m_pCovMatrix[uiDepth] = NULL;
              }
          }
          delete[]m_pCovMatrix; m_pCovMatrix = NULL;
      }
#if VCEG_AZ08_FAST_DERIVE_KLT
      if (m_pppdTmpEigenVector != NULL)
      {
          UInt blkSize = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX-1];
          UInt uiDim = blkSize*blkSize;
          for (UInt k = 0; k < uiDim; k++)
          {
              delete[]m_pppdTmpEigenVector[k];
              m_pppdTmpEigenVector[k] = NULL;
          }
          delete[] m_pppdTmpEigenVector;
          m_pppdTmpEigenVector = NULL;
      }
#endif
      if (m_pppdEigenVector != NULL)
      {
          for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
          {
              UInt blkSize = g_uiDepth2Width[uiDepth];
              UInt uiDim = blkSize*blkSize;
              for (UInt k = 0; k < uiDim; k++)
              {
                  if (m_pppdEigenVector[uiDepth][k] != NULL)
                  {
                      delete[]m_pppdEigenVector[uiDepth][k]; m_pppdEigenVector[uiDepth][k] = NULL;
                  }
              }
              delete [] m_pppdEigenVector[uiDepth];
          }
          delete[]m_pppdEigenVector; m_pppdEigenVector = NULL;
      }
      if (m_pppsEigenVector != NULL)
      {
          for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
          {
              UInt blkSize = g_uiDepth2Width[uiDepth];
              UInt uiDim = blkSize*blkSize;
              for (UInt k = 0; k < uiDim; k++)
              {
                  if (m_pppsEigenVector[uiDepth][k] != NULL)
                  {
                      delete[]m_pppsEigenVector[uiDepth][k]; m_pppsEigenVector[uiDepth][k] = NULL;
                  }
              }
          }
          delete[]m_pppsEigenVector; m_pppsEigenVector = NULL;
      }
#endif
#if VCEG_AZ08_USE_KLT
  }
#endif
}

#if ADAPTIVE_QP_SELECTION
Void TComTrQuant::storeSliceQpNext(TComSlice* pcSlice)
{
  // NOTE: does this work with negative QPs or when some blocks are transquant-bypass enabled?

  Int qpBase = pcSlice->getSliceQpBase();
  Int sliceQpused = pcSlice->getSliceQp();
  Int sliceQpnext;
  Double alpha = qpBase < 17 ? 0.5 : 1;

  Int cnt=0;
  for(Int u=1; u<=LEVEL_RANGE; u++)
  {
    cnt += m_sliceNsamples[u] ;
  }

  if( !m_useRDOQ )
  {
    sliceQpused = qpBase;
    alpha = 0.5;
  }

  if( cnt > 120 )
  {
    Double sum = 0;
    Int k = 0;
    for(Int u=1; u<LEVEL_RANGE; u++)
    {
      sum += u*m_sliceSumC[u];
      k += u*u*m_sliceNsamples[u];
    }

    Int v;
    Double q[MAX_QP+1] ;
    for(v=0; v<=MAX_QP; v++)
    {
      q[v] = (Double)(g_invQuantScales[v%6] * (1<<(v/6)))/64 ;
    }

    Double qnext = sum/k * q[sliceQpused] / (1<<ARL_C_PRECISION);

    for(v=0; v<MAX_QP; v++)
    {
      if(qnext < alpha * q[v] + (1 - alpha) * q[v+1] )
      {
        break;
      }
    }
    sliceQpnext = Clip3(sliceQpused - 3, sliceQpused + 3, v);
  }
  else
  {
    sliceQpnext = sliceQpused;
  }

  m_qpDelta[qpBase] = sliceQpnext - qpBase;
}

Void TComTrQuant::initSliceQpDelta()
{
  for(Int qp=0; qp<=MAX_QP; qp++)
  {
    m_qpDelta[qp] = qp < 17 ? 0 : 1;
  }
}

Void TComTrQuant::clearSliceARLCnt()
{
  memset(m_sliceSumC, 0, sizeof(Double)*(LEVEL_RANGE+1));
  memset(m_sliceNsamples, 0, sizeof(Int)*(LEVEL_RANGE+1));
}
#endif



#if MATRIX_MULT
/** NxN forward transform (2D) using brute force matrix multiplication (3 nested loops)
 *  \param block pointer to input data (residual)
 *  \param coeff pointer to output data (transform coefficients)
 *  \param uiStride stride of input data
 *  \param uiTrSize transform size (uiTrSize x uiTrSize)
 *  \param uiMode is Intra Prediction mode used in Mode-Dependent DCT/DST only
 */
#if COM16_C806_EMT
Void xTr_EMT(Int bitDepth, Pel *block, TCoeff *coeff, UInt uiStride, UInt uiTrSize, Bool useDST, const Int maxLog2TrDynamicRange, UChar ucMode, UChar ucTrIdx)
{
  UInt i,j,k;
  TCoeff iSum;
  TCoeff tmp[MAX_TU_SIZE * MAX_TU_SIZE];
  const TMatrixCoeff *iTH, *iTV;
  UInt uiLog2TrSize = g_aucConvertToBit[ uiTrSize ] + 2;

#if COM16_C806_T64
  Int zo = ( ( ucMode == INTER_MODE_IDX && uiTrSize>=32 && ucTrIdx != DCT2_EMT ) || uiTrSize==64 ) ? 1 : 0;
#else
  Int zo = ( ucMode == INTER_MODE_IDX && uiTrSize>=32 && ucTrIdx != DCT2_EMT ) ? 1 : 0;
#endif

  UInt  nTrIdxHor = DCT2, nTrIdxVer = DCT2;
  if( ucMode != INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    UInt  nTrSubsetHor = g_aucTrSetHorz[ucMode];
    UInt  nTrSubsetVer = g_aucTrSetVert[ucMode];
    nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor][ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer][ucTrIdx>>1];
  }
  if ( ucMode == INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    nTrIdxHor = g_aiTrSubsetInter[ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetInter[ucTrIdx>>1];
  }

  if (uiTrSize==4)
  {
    iTH  = g_aiTr4[nTrIdxHor][0];
    iTV  = g_aiTr4[nTrIdxVer][0];
  }
  else if (uiTrSize==8)
  {
    iTH  = g_aiTr8[nTrIdxHor][0];
    iTV  = g_aiTr8[nTrIdxVer][0];
  }
  else if (uiTrSize==16)
  {
    iTH  = g_aiTr16[nTrIdxHor][0];
    iTV  = g_aiTr16[nTrIdxVer][0];
  }
  else if (uiTrSize==32)
  {
    iTH  = g_aiTr32[nTrIdxHor][0];
    iTV  = g_aiTr32[nTrIdxVer][0];
  }
#if COM16_C806_T64
  else if (uiTrSize==64)
  {
    assert( ucTrIdx == DCT2_EMT );
    iTH = iTV = g_aiTr64[DCT2][0];
  }
#endif
  else
  {
    assert(0);
  }

  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];
  
  Int shift_1st = (uiLog2TrSize +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
  Int shift_2nd = uiLog2TrSize + TRANSFORM_MATRIX_SHIFT + COM16_C806_TRANS_PREC;
  const Int add_1st = (shift_1st>0) ? (1<<(shift_1st-1)) : 0;
  const Int add_2nd = 1<<(shift_2nd-1);

  /* Horizontal transform */
  for (i=0; i<(uiTrSize>>zo); i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += iTH[i*uiTrSize+k]*block[j*uiStride+k];
      }
      tmp[i*uiTrSize+j] = (iSum + add_1st)>>shift_1st;
    }
  }

  /* Vertical transform */
  for (i=0; i<(uiTrSize>>zo); i++)
  {
    for (j=0; j<(uiTrSize>>zo); j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += iTV[i*uiTrSize+k]*tmp[j*uiTrSize+k];
      }
      coeff[i*uiTrSize+j] = (iSum + add_2nd)>>shift_2nd;
    }
  }
  if(zo)
  {
    memset( coeff + uiTrSize*uiTrSize/2, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
    coeff += uiTrSize/2;
    for (j=0; j<uiTrSize/2; j++)
    {
      memset( coeff, 0, sizeof(TCoeff)*uiTrSize/2 );
      coeff += uiTrSize;
    }
  }
}

Void xITr_EMT(Int bitDepth, TCoeff *coeff, Pel *block, UInt uiStride, UInt uiTrSize, Bool useDST, const Int maxLog2TrDynamicRange, UChar ucMode, UChar ucTrIdx)
{
  UInt i,j,k;
  TCoeff iSum;
  TCoeff tmp[MAX_TU_SIZE * MAX_TU_SIZE];
  const TMatrixCoeff *iTH, *iTV;

#if COM16_C806_T64
  Int zo = ( ( ucMode == INTER_MODE_IDX && uiTrSize>=32 && ucTrIdx != DCT2_EMT ) || uiTrSize==64 ) ? 1 : 0;
#else
  Int zo = ( ucMode == INTER_MODE_IDX && uiTrSize>=32 && ucTrIdx != DCT2_EMT ) ? 1 : 0;
#endif

  UInt  nTrIdxHor = DCT2, nTrIdxVer = DCT2;
  if( ucMode != INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    UInt  nTrSubsetHor = g_aucTrSetHorz[ucMode];
    UInt  nTrSubsetVer = g_aucTrSetVert[ucMode];
    nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor][ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer][ucTrIdx>>1];
  }
  if ( ucMode == INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    nTrIdxHor = g_aiTrSubsetInter[ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetInter[ucTrIdx>>1];
  }

  if (uiTrSize==4)
  {
    iTH  = g_aiTr4[nTrIdxHor][0];
    iTV  = g_aiTr4[nTrIdxVer][0];
  }
  else if (uiTrSize==8)
  {
    iTH  = g_aiTr8[nTrIdxHor][0];
    iTV  = g_aiTr8[nTrIdxVer][0];
  }
  else if (uiTrSize==16)
  {
    iTH  = g_aiTr16[nTrIdxHor][0];
    iTV  = g_aiTr16[nTrIdxVer][0];
  }
  else if (uiTrSize==32)
  {
    iTH  = g_aiTr32[nTrIdxHor][0];
    iTV  = g_aiTr32[nTrIdxVer][0];
  }
#if COM16_C806_T64
  else if (uiTrSize==64)
  {
    assert( ucTrIdx == DCT2_EMT );
    iTH = iTV = g_aiTr64[DCT2][0];
  }
#endif
  else
  {
    assert(0);
  }

  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_INVERSE];


  Int shift_1st = TRANSFORM_MATRIX_SHIFT + 1 + COM16_C806_TRANS_PREC; //1 has been added to shift_1st at the expense of shift_2nd
  Int shift_2nd = (TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1) - bitDepth + COM16_C806_TRANS_PREC;
  const TCoeff clipMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff clipMaximum =  (1 << maxLog2TrDynamicRange) - 1;
  assert(shift_2nd>=0);
  const Int add_1st = 1<<(shift_1st-1);
  const Int add_2nd = (shift_2nd>0) ? (1<<(shift_2nd-1)) : 0;

  /* Horizontal transform */
  for (i=0; i<uiTrSize; i++)
  {
    for (j=0; j<uiTrSize>>zo; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize>>zo; k++)
      {
        iSum += iTV[k*uiTrSize+i]*coeff[k*uiTrSize+j];
      }

      // Clipping here is not in the standard, but is used to protect the "Pel" data type into which the inverse-transformed samples will be copied
      tmp[i*uiTrSize+j] = Clip3<TCoeff>(clipMinimum, clipMaximum, (iSum + add_1st)>>shift_1st);
    }
  }

  /* Vertical transform */
  for (i=0; i<uiTrSize; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize>>zo; k++)
      {
        iSum += iTH[k*uiTrSize+j]*tmp[i*uiTrSize+k];
      }

      block[i*uiStride+j] = Clip3<TCoeff>(std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max(), (iSum + add_2nd)>>shift_2nd);
    }
  }
}
#endif

Void xTr(Int bitDepth, Pel *block, TCoeff *coeff, UInt uiStride, UInt uiTrSize, Bool useDST, const Int maxLog2TrDynamicRange)
{
  UInt i,j,k;
  TCoeff iSum;
  TCoeff tmp[MAX_TU_SIZE * MAX_TU_SIZE];
  const TMatrixCoeff *iT;
  UInt uiLog2TrSize = g_aucConvertToBit[ uiTrSize ] + 2;

#if COM16_C806_T64
  Int zo = 0;
#endif

  if (uiTrSize==4)
  {
    iT  = (useDST ? g_as_DST_MAT_4[TRANSFORM_FORWARD][0] : g_aiT4[TRANSFORM_FORWARD][0]);
  }
  else if (uiTrSize==8)
  {
    iT = g_aiT8[TRANSFORM_FORWARD][0];
  }
  else if (uiTrSize==16)
  {
    iT = g_aiT16[TRANSFORM_FORWARD][0];
  }
  else if (uiTrSize==32)
  {
    iT = g_aiT32[TRANSFORM_FORWARD][0];
  }
#if COM16_C806_T64
  else if (uiTrSize==64)
  {
    iT = g_aiTr64[DCT2][0];
    zo = 1;
  }
#endif
  else
  {
    assert(0);
  }

  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];

#if !COM16_C806_T64
  const 
#endif
    Int shift_1st = (uiLog2TrSize +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange;
#if !COM16_C806_T64  
  const 
#endif
    Int shift_2nd = uiLog2TrSize + TRANSFORM_MATRIX_SHIFT;
#if COM16_C806_T64
  if( uiTrSize==64 )
  {
    shift_1st += COM16_C806_TRANS_PREC;
    shift_2nd += COM16_C806_TRANS_PREC;
  }
#endif
  const Int add_1st = (shift_1st>0) ? (1<<(shift_1st-1)) : 0;
  const Int add_2nd = 1<<(shift_2nd-1);

  /* Horizontal transform */

#if COM16_C806_T64
  for (i=0; i<(uiTrSize>>zo); i++)
#else
  for (i=0; i<uiTrSize; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += iT[i*uiTrSize+k]*block[j*uiStride+k];
      }
      tmp[i*uiTrSize+j] = (iSum + add_1st)>>shift_1st;
    }
  }

  /* Vertical transform */
#if COM16_C806_T64
  for (i=0; i<(uiTrSize>>zo); i++)
#else
  for (i=0; i<uiTrSize; i++)
#endif
  {
#if COM16_C806_T64
    for (j=0; j<(uiTrSize>>zo); j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += iT[i*uiTrSize+k]*tmp[j*uiTrSize+k];
      }
      coeff[i*uiTrSize+j] = (iSum + add_2nd)>>shift_2nd;
    }
  }
#if COM16_C806_T64
  if(zo)
  {
    memset( coeff + uiTrSize*uiTrSize/2, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
    coeff += uiTrSize/2;
    for (j=0; j<uiTrSize/2; j++)
    {
      memset( coeff, 0, sizeof(TCoeff)*uiTrSize/2 );
      coeff += uiTrSize;
    }
  }
#endif
}

/** NxN inverse transform (2D) using brute force matrix multiplication (3 nested loops)
 *  \param coeff pointer to input data (transform coefficients)
 *  \param block pointer to output data (residual)
 *  \param uiStride stride of output data
 *  \param uiTrSize transform size (uiTrSize x uiTrSize)
 *  \param uiMode is Intra Prediction mode used in Mode-Dependent DCT/DST only
 */
Void xITr(Int bitDepth, TCoeff *coeff, Pel *block, UInt uiStride, UInt uiTrSize, Bool useDST, const Int maxLog2TrDynamicRange)
{
  UInt i,j,k;
  TCoeff iSum;
  TCoeff tmp[MAX_TU_SIZE * MAX_TU_SIZE];
  const TMatrixCoeff *iT;
#if COM16_C806_T64
  Int zo = 0;
#endif

  if (uiTrSize==4)
  {
    iT  = (useDST ? g_as_DST_MAT_4[TRANSFORM_INVERSE][0] : g_aiT4[TRANSFORM_INVERSE][0]);
  }
  else if (uiTrSize==8)
  {
    iT = g_aiT8[TRANSFORM_INVERSE][0];
  }
  else if (uiTrSize==16)
  {
    iT = g_aiT16[TRANSFORM_INVERSE][0];
  }
  else if (uiTrSize==32)
  {
    iT = g_aiT32[TRANSFORM_INVERSE][0];
  }
#if COM16_C806_T64
  else if (uiTrSize==64)
  {
    iT = g_aiTr64[DCT2][0];
    zo = 1;
  }
#endif
  else
  {
    assert(0);
  }

  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_INVERSE];

#if !COM16_C806_T64
  const 
#endif
    Int shift_1st = TRANSFORM_MATRIX_SHIFT + 1; //1 has been added to shift_1st at the expense of shift_2nd
#if !COM16_C806_T64
  const 
#endif
    Int shift_2nd = (TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1) - bitDepth;
#if COM16_C806_T64
  if( uiTrSize==64 )
  {
    shift_1st += COM16_C806_TRANS_PREC;
    shift_2nd += COM16_C806_TRANS_PREC;
  }
#endif
  const TCoeff clipMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff clipMaximum =  (1 << maxLog2TrDynamicRange) - 1;
  assert(shift_2nd>=0);
  const Int add_1st = 1<<(shift_1st-1);
  const Int add_2nd = (shift_2nd>0) ? (1<<(shift_2nd-1)) : 0;

  /* Horizontal transform */
  for (i=0; i<uiTrSize; i++)
  {
#if COM16_C806_T64
    for (j=0; j<uiTrSize>>zo; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
#if COM16_C806_T64
      for (k=0; k<uiTrSize>>zo; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += iT[k*uiTrSize+i]*coeff[k*uiTrSize+j];
      }

      // Clipping here is not in the standard, but is used to protect the "Pel" data type into which the inverse-transformed samples will be copied
      tmp[i*uiTrSize+j] = Clip3<TCoeff>(clipMinimum, clipMaximum, (iSum + add_1st)>>shift_1st);
    }
  }

  /* Vertical transform */
  for (i=0; i<uiTrSize; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if COM16_C806_T64
      for (k=0; k<uiTrSize>>zo; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += iT[k*uiTrSize+j]*tmp[i*uiTrSize+k];
      }

      block[i*uiStride+j] = Clip3<TCoeff>(std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max(), (iSum + add_2nd)>>shift_2nd);
    }
  }
}

#endif //MATRIX_MULT


/** 4x4 forward transform implemented using partial butterfly structure (1D)
 *  \param src   input data (residual)
 *  \param dst   output data (transform coefficients)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterfly4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine)
#else
Void partialButterfly4(TCoeff *src, TCoeff *dst, Int shift, Int line)
#endif
{
  Int j;
  TCoeff E[2],O[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;

#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* E and O */
    E[0] = src[0] + src[3];
    O[0] = src[0] - src[3];
    E[1] = src[1] + src[2];
    O[1] = src[1] - src[2];

    dst[0]      = (g_aiT4[TRANSFORM_FORWARD][0][0]*E[0] + g_aiT4[TRANSFORM_FORWARD][0][1]*E[1] + add)>>shift;
    dst[2*line] = (g_aiT4[TRANSFORM_FORWARD][2][0]*E[0] + g_aiT4[TRANSFORM_FORWARD][2][1]*E[1] + add)>>shift;
    dst[line]   = (g_aiT4[TRANSFORM_FORWARD][1][0]*O[0] + g_aiT4[TRANSFORM_FORWARD][1][1]*O[1] + add)>>shift;
    dst[3*line] = (g_aiT4[TRANSFORM_FORWARD][3][0]*O[0] + g_aiT4[TRANSFORM_FORWARD][3][1]*O[1] + add)>>shift;

    src += 4;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    Int size = sizeof(TCoeff)*iSkipLine;
    for( j = 0; j < 4; j++ )
    { 
      memset( dst, 0, size );
      dst += line;
    }
  }
#endif
}

// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm
// give identical results
Void fastForwardDst(TCoeff *block, TCoeff *coeff, Int shift)  // input block, output coeff
{
  Int i;
  TCoeff c[4];
  TCoeff rnd_factor = (shift > 0) ? (1<<(shift-1)) : 0;
  for (i=0; i<4; i++)
  {
    // Intermediate Variables
    c[0] = block[4*i+0];
    c[1] = block[4*i+1];
    c[2] = block[4*i+2];
    c[3] = block[4*i+3];

    for (Int row = 0; row < 4; row++)
    {
      TCoeff result = 0;
      for (Int column = 0; column < 4; column++)
      {
        result += c[column] * g_as_DST_MAT_4[TRANSFORM_FORWARD][row][column]; // use the defined matrix, rather than hard-wired numbers
      }

      coeff[(row * 4) + i] = rightShift((result + rnd_factor), shift);
    }
  }
}

Void fastInverseDst(TCoeff *tmp, TCoeff *block, Int shift, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
{
  Int i;
  TCoeff c[4];
  TCoeff rnd_factor = (shift > 0) ? (1<<(shift-1)) : 0;
  for (i=0; i<4; i++)
  {
    // Intermediate Variables
    c[0] = tmp[   i];
    c[1] = tmp[4 +i];
    c[2] = tmp[8 +i];
    c[3] = tmp[12+i];

    for (Int column = 0; column < 4; column++)
    {
      TCoeff &result = block[(i * 4) + column];

      result = 0;
      for (Int row = 0; row < 4; row++)
      {
        result += c[row] * g_as_DST_MAT_4[TRANSFORM_INVERSE][row][column]; // use the defined matrix, rather than hard-wired numbers
      }

      result = Clip3( outputMinimum, outputMaximum, rightShift((result + rnd_factor), shift));
    }
  }
}

/** 4x4 inverse transform implemented using partial butterfly structure (1D)
 *  \param src   input data (transform coefficients)
 *  \param dst   output data (residual)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 *  \param outputMinimum  minimum for clipping
 *  \param outputMaximum  maximum for clipping
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterflyInverse4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
Void partialButterflyInverse4(TCoeff *src, TCoeff *dst, Int shift, Int line, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j;
  TCoeff E[2],O[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;

#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    O[0] = g_aiT4[TRANSFORM_INVERSE][1][0]*src[line] + g_aiT4[TRANSFORM_INVERSE][3][0]*src[3*line];
    O[1] = g_aiT4[TRANSFORM_INVERSE][1][1]*src[line] + g_aiT4[TRANSFORM_INVERSE][3][1]*src[3*line];
    E[0] = g_aiT4[TRANSFORM_INVERSE][0][0]*src[0]    + g_aiT4[TRANSFORM_INVERSE][2][0]*src[2*line];
    E[1] = g_aiT4[TRANSFORM_INVERSE][0][1]*src[0]    + g_aiT4[TRANSFORM_INVERSE][2][1]*src[2*line];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    dst[0] = Clip3( outputMinimum, outputMaximum, (E[0] + O[0] + add)>>shift );
    dst[1] = Clip3( outputMinimum, outputMaximum, (E[1] + O[1] + add)>>shift );
    dst[2] = Clip3( outputMinimum, outputMaximum, (E[1] - O[1] + add)>>shift );
    dst[3] = Clip3( outputMinimum, outputMaximum, (E[0] - O[0] + add)>>shift );

    src   ++;
    dst += 4;
  }
#if JVET_D0077_TRANSFORM_OPT
  memset(dst, 0, (iSkipLine<<2)*sizeof(TCoeff));
#endif
}

/** 8x8 forward transform implemented using partial butterfly structure (1D)
 *  \param src   input data (residual)
 *  \param dst   output data (transform coefficients)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterfly8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine)
#else
Void partialButterfly8(TCoeff *src, TCoeff *dst, Int shift, Int line)
#endif
{
  Int j,k;
  TCoeff E[4],O[4];
  TCoeff EE[2],EO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;

#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* E and O*/
    for (k=0;k<4;k++)
    {
      E[k] = src[k] + src[7-k];
      O[k] = src[k] - src[7-k];
    }
    /* EE and EO */
    EE[0] = E[0] + E[3];
    EO[0] = E[0] - E[3];
    EE[1] = E[1] + E[2];
    EO[1] = E[1] - E[2];

    dst[0]      = (g_aiT8[TRANSFORM_FORWARD][0][0]*EE[0] + g_aiT8[TRANSFORM_FORWARD][0][1]*EE[1] + add)>>shift;
    dst[4*line] = (g_aiT8[TRANSFORM_FORWARD][4][0]*EE[0] + g_aiT8[TRANSFORM_FORWARD][4][1]*EE[1] + add)>>shift;
    dst[2*line] = (g_aiT8[TRANSFORM_FORWARD][2][0]*EO[0] + g_aiT8[TRANSFORM_FORWARD][2][1]*EO[1] + add)>>shift;
    dst[6*line] = (g_aiT8[TRANSFORM_FORWARD][6][0]*EO[0] + g_aiT8[TRANSFORM_FORWARD][6][1]*EO[1] + add)>>shift;

    dst[line]   = (g_aiT8[TRANSFORM_FORWARD][1][0]*O[0] + g_aiT8[TRANSFORM_FORWARD][1][1]*O[1] + g_aiT8[TRANSFORM_FORWARD][1][2]*O[2] + g_aiT8[TRANSFORM_FORWARD][1][3]*O[3] + add)>>shift;
    dst[3*line] = (g_aiT8[TRANSFORM_FORWARD][3][0]*O[0] + g_aiT8[TRANSFORM_FORWARD][3][1]*O[1] + g_aiT8[TRANSFORM_FORWARD][3][2]*O[2] + g_aiT8[TRANSFORM_FORWARD][3][3]*O[3] + add)>>shift;
    dst[5*line] = (g_aiT8[TRANSFORM_FORWARD][5][0]*O[0] + g_aiT8[TRANSFORM_FORWARD][5][1]*O[1] + g_aiT8[TRANSFORM_FORWARD][5][2]*O[2] + g_aiT8[TRANSFORM_FORWARD][5][3]*O[3] + add)>>shift;
    dst[7*line] = (g_aiT8[TRANSFORM_FORWARD][7][0]*O[0] + g_aiT8[TRANSFORM_FORWARD][7][1]*O[1] + g_aiT8[TRANSFORM_FORWARD][7][2]*O[2] + g_aiT8[TRANSFORM_FORWARD][7][3]*O[3] + add)>>shift;

    src += 8;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    Int size = sizeof(TCoeff)*iSkipLine;
    for( j = 0; j < 8; j++ )
    { 
      memset( dst, 0, size );
      dst += line;
    }
  }
#endif
}

/** 8x8 inverse transform implemented using partial butterfly structure (1D)
 *  \param src   input data (transform coefficients)
 *  \param dst   output data (residual)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 *  \param outputMinimum  minimum for clipping
 *  \param outputMaximum  maximum for clipping
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterflyInverse8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
Void partialButterflyInverse8(TCoeff *src, TCoeff *dst, Int shift, Int line, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  TCoeff E[4];
  TCoeff EE[2];
#if _WIN64 && _MSC_VER == 1900
  // Workaround for avoid bug in Visual Studio 2015 Update 3 optimizer 
  // https://hevc.hhi.fraunhofer.de/trac/hevc/ticket/1463
  volatile
#endif
  TCoeff O[4], EO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;
#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<4;k++)
    {
      O[k] = g_aiT8[TRANSFORM_INVERSE][ 1][k]*src[line]   + g_aiT8[TRANSFORM_INVERSE][ 3][k]*src[3*line] +
             g_aiT8[TRANSFORM_INVERSE][ 5][k]*src[5*line] + g_aiT8[TRANSFORM_INVERSE][ 7][k]*src[7*line];
    }

    EO[0] = g_aiT8[TRANSFORM_INVERSE][2][0]*src[ 2*line ] + g_aiT8[TRANSFORM_INVERSE][6][0]*src[ 6*line ];
    EO[1] = g_aiT8[TRANSFORM_INVERSE][2][1]*src[ 2*line ] + g_aiT8[TRANSFORM_INVERSE][6][1]*src[ 6*line ];
    EE[0] = g_aiT8[TRANSFORM_INVERSE][0][0]*src[ 0      ] + g_aiT8[TRANSFORM_INVERSE][4][0]*src[ 4*line ];
    EE[1] = g_aiT8[TRANSFORM_INVERSE][0][1]*src[ 0      ] + g_aiT8[TRANSFORM_INVERSE][4][1]*src[ 4*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    E[0] = EE[0] + EO[0];
    E[3] = EE[0] - EO[0];
    E[1] = EE[1] + EO[1];
    E[2] = EE[1] - EO[1];
    for (k=0;k<4;k++)
    {
      dst[ k   ] = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[ k+4 ] = Clip3( outputMinimum, outputMaximum, (E[3-k] - O[3-k] + add)>>shift );
    }
    src ++;
    dst += 8;
  }
#if JVET_D0077_TRANSFORM_OPT
  memset(dst, 0, (iSkipLine<<3)*sizeof(TCoeff));
#endif
}

/** 16x16 forward transform implemented using partial butterfly structure (1D)
 *  \param src   input data (residual)
 *  \param dst   output data (transform coefficients)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterfly16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine)
#else
Void partialButterfly16(TCoeff *src, TCoeff *dst, Int shift, Int line)
#endif
{
  Int j,k;
  TCoeff E[8],O[8];
  TCoeff EE[4],EO[4];
  TCoeff EEE[2],EEO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;
#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* E and O*/
    for (k=0;k<8;k++)
    {
      E[k] = src[k] + src[15-k];
      O[k] = src[k] - src[15-k];
    }
    /* EE and EO */
    for (k=0;k<4;k++)
    {
      EE[k] = E[k] + E[7-k];
      EO[k] = E[k] - E[7-k];
    }
    /* EEE and EEO */
    EEE[0] = EE[0] + EE[3];
    EEO[0] = EE[0] - EE[3];
    EEE[1] = EE[1] + EE[2];
    EEO[1] = EE[1] - EE[2];

    dst[ 0      ] = (g_aiT16[TRANSFORM_FORWARD][ 0][0]*EEE[0] + g_aiT16[TRANSFORM_FORWARD][ 0][1]*EEE[1] + add)>>shift;
    dst[ 8*line ] = (g_aiT16[TRANSFORM_FORWARD][ 8][0]*EEE[0] + g_aiT16[TRANSFORM_FORWARD][ 8][1]*EEE[1] + add)>>shift;
    dst[ 4*line ] = (g_aiT16[TRANSFORM_FORWARD][ 4][0]*EEO[0] + g_aiT16[TRANSFORM_FORWARD][ 4][1]*EEO[1] + add)>>shift;
    dst[ 12*line] = (g_aiT16[TRANSFORM_FORWARD][12][0]*EEO[0] + g_aiT16[TRANSFORM_FORWARD][12][1]*EEO[1] + add)>>shift;

    for (k=2;k<16;k+=4)
    {
      dst[ k*line ] = (g_aiT16[TRANSFORM_FORWARD][k][0]*EO[0] + g_aiT16[TRANSFORM_FORWARD][k][1]*EO[1] +
                       g_aiT16[TRANSFORM_FORWARD][k][2]*EO[2] + g_aiT16[TRANSFORM_FORWARD][k][3]*EO[3] + add)>>shift;
    }

    for (k=1;k<16;k+=2)
    {
      dst[ k*line ] = (g_aiT16[TRANSFORM_FORWARD][k][0]*O[0] + g_aiT16[TRANSFORM_FORWARD][k][1]*O[1] +
                       g_aiT16[TRANSFORM_FORWARD][k][2]*O[2] + g_aiT16[TRANSFORM_FORWARD][k][3]*O[3] +
                       g_aiT16[TRANSFORM_FORWARD][k][4]*O[4] + g_aiT16[TRANSFORM_FORWARD][k][5]*O[5] +
                       g_aiT16[TRANSFORM_FORWARD][k][6]*O[6] + g_aiT16[TRANSFORM_FORWARD][k][7]*O[7] + add)>>shift;
    }

    src += 16;
    dst ++;

  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    Int size = sizeof(TCoeff)*iSkipLine;
    for( j = 0; j < 16; j++ )
    { 
      memset( dst, 0, size );
      dst += line;
    }
  }
#endif
}

/** 16x16 inverse transform implemented using partial butterfly structure (1D)
 *  \param src            input data (transform coefficients)
 *  \param dst            output data (residual)
 *  \param shift          specifies right shift after 1D transform
 *  \param line
 *  \param outputMinimum  minimum for clipping
 *  \param outputMaximum  maximum for clipping
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterflyInverse16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
Void partialButterflyInverse16(TCoeff *src, TCoeff *dst, Int shift, Int line, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  TCoeff E[8];
  TCoeff EE[4];
  TCoeff EEE[2];
#if _WIN64 && _MSC_VER == 1900
  // Workaround for avoid bug in Visual Studio 2015 Update 3 optimizer 
  // https://hevc.hhi.fraunhofer.de/trac/hevc/ticket/1463
  volatile
#endif
  TCoeff O[8], EO[4], EEO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;
#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<8;k++)
    {
      O[k] = g_aiT16[TRANSFORM_INVERSE][ 1][k]*src[ line]   + g_aiT16[TRANSFORM_INVERSE][ 3][k]*src[ 3*line] +
             g_aiT16[TRANSFORM_INVERSE][ 5][k]*src[ 5*line] + g_aiT16[TRANSFORM_INVERSE][ 7][k]*src[ 7*line] +
             g_aiT16[TRANSFORM_INVERSE][ 9][k]*src[ 9*line] + g_aiT16[TRANSFORM_INVERSE][11][k]*src[11*line] +
             g_aiT16[TRANSFORM_INVERSE][13][k]*src[13*line] + g_aiT16[TRANSFORM_INVERSE][15][k]*src[15*line];
    }
    for (k=0;k<4;k++)
    {
      EO[k] = g_aiT16[TRANSFORM_INVERSE][ 2][k]*src[ 2*line] + g_aiT16[TRANSFORM_INVERSE][ 6][k]*src[ 6*line] +
              g_aiT16[TRANSFORM_INVERSE][10][k]*src[10*line] + g_aiT16[TRANSFORM_INVERSE][14][k]*src[14*line];
    }
    EEO[0] = g_aiT16[TRANSFORM_INVERSE][4][0]*src[ 4*line ] + g_aiT16[TRANSFORM_INVERSE][12][0]*src[ 12*line ];
    EEE[0] = g_aiT16[TRANSFORM_INVERSE][0][0]*src[ 0      ] + g_aiT16[TRANSFORM_INVERSE][ 8][0]*src[ 8*line  ];
    EEO[1] = g_aiT16[TRANSFORM_INVERSE][4][1]*src[ 4*line ] + g_aiT16[TRANSFORM_INVERSE][12][1]*src[ 12*line ];
    EEE[1] = g_aiT16[TRANSFORM_INVERSE][0][1]*src[ 0      ] + g_aiT16[TRANSFORM_INVERSE][ 8][1]*src[ 8*line  ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    for (k=0;k<2;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }
    for (k=0;k<4;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }
    for (k=0;k<8;k++)
    {
      dst[k]   = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[k+8] = Clip3( outputMinimum, outputMaximum, (E[7-k] - O[7-k] + add)>>shift );
    }
    src ++;
    dst += 16;
  }
#if JVET_D0077_TRANSFORM_OPT
  memset(dst, 0, (iSkipLine<<4)*sizeof(TCoeff));
#endif
}

/** 32x32 forward transform implemented using partial butterfly structure (1D)
 *  \param src   input data (residual)
 *  \param dst   output data (transform coefficients)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterfly32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine)
#else
Void partialButterfly32(TCoeff *src, TCoeff *dst, Int shift, Int line)
#endif
{
  Int j,k;
  TCoeff E[16],O[16];
  TCoeff EE[8],EO[8];
  TCoeff EEE[4],EEO[4];
  TCoeff EEEE[2],EEEO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;
#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* E and O*/
    for (k=0;k<16;k++)
    {
      E[k] = src[k] + src[31-k];
      O[k] = src[k] - src[31-k];
    }
    /* EE and EO */
    for (k=0;k<8;k++)
    {
      EE[k] = E[k] + E[15-k];
      EO[k] = E[k] - E[15-k];
    }
    /* EEE and EEO */
    for (k=0;k<4;k++)
    {
      EEE[k] = EE[k] + EE[7-k];
      EEO[k] = EE[k] - EE[7-k];
    }
    /* EEEE and EEEO */
    EEEE[0] = EEE[0] + EEE[3];
    EEEO[0] = EEE[0] - EEE[3];
    EEEE[1] = EEE[1] + EEE[2];
    EEEO[1] = EEE[1] - EEE[2];

    dst[ 0       ] = (g_aiT32[TRANSFORM_FORWARD][ 0][0]*EEEE[0] + g_aiT32[TRANSFORM_FORWARD][ 0][1]*EEEE[1] + add)>>shift;
    dst[ 16*line ] = (g_aiT32[TRANSFORM_FORWARD][16][0]*EEEE[0] + g_aiT32[TRANSFORM_FORWARD][16][1]*EEEE[1] + add)>>shift;
    dst[ 8*line  ] = (g_aiT32[TRANSFORM_FORWARD][ 8][0]*EEEO[0] + g_aiT32[TRANSFORM_FORWARD][ 8][1]*EEEO[1] + add)>>shift;
    dst[ 24*line ] = (g_aiT32[TRANSFORM_FORWARD][24][0]*EEEO[0] + g_aiT32[TRANSFORM_FORWARD][24][1]*EEEO[1] + add)>>shift;
    for (k=4;k<32;k+=8)
    {
      dst[ k*line ] = (g_aiT32[TRANSFORM_FORWARD][k][0]*EEO[0] + g_aiT32[TRANSFORM_FORWARD][k][1]*EEO[1] +
                       g_aiT32[TRANSFORM_FORWARD][k][2]*EEO[2] + g_aiT32[TRANSFORM_FORWARD][k][3]*EEO[3] + add)>>shift;
    }
    for (k=2;k<32;k+=4)
    {
      dst[ k*line ] = (g_aiT32[TRANSFORM_FORWARD][k][0]*EO[0] + g_aiT32[TRANSFORM_FORWARD][k][1]*EO[1] +
                       g_aiT32[TRANSFORM_FORWARD][k][2]*EO[2] + g_aiT32[TRANSFORM_FORWARD][k][3]*EO[3] +
                       g_aiT32[TRANSFORM_FORWARD][k][4]*EO[4] + g_aiT32[TRANSFORM_FORWARD][k][5]*EO[5] +
                       g_aiT32[TRANSFORM_FORWARD][k][6]*EO[6] + g_aiT32[TRANSFORM_FORWARD][k][7]*EO[7] + add)>>shift;
    }
    for (k=1;k<32;k+=2)
    {
      dst[ k*line ] = (g_aiT32[TRANSFORM_FORWARD][k][ 0]*O[ 0] + g_aiT32[TRANSFORM_FORWARD][k][ 1]*O[ 1] +
                       g_aiT32[TRANSFORM_FORWARD][k][ 2]*O[ 2] + g_aiT32[TRANSFORM_FORWARD][k][ 3]*O[ 3] +
                       g_aiT32[TRANSFORM_FORWARD][k][ 4]*O[ 4] + g_aiT32[TRANSFORM_FORWARD][k][ 5]*O[ 5] +
                       g_aiT32[TRANSFORM_FORWARD][k][ 6]*O[ 6] + g_aiT32[TRANSFORM_FORWARD][k][ 7]*O[ 7] +
                       g_aiT32[TRANSFORM_FORWARD][k][ 8]*O[ 8] + g_aiT32[TRANSFORM_FORWARD][k][ 9]*O[ 9] +
                       g_aiT32[TRANSFORM_FORWARD][k][10]*O[10] + g_aiT32[TRANSFORM_FORWARD][k][11]*O[11] +
                       g_aiT32[TRANSFORM_FORWARD][k][12]*O[12] + g_aiT32[TRANSFORM_FORWARD][k][13]*O[13] +
                       g_aiT32[TRANSFORM_FORWARD][k][14]*O[14] + g_aiT32[TRANSFORM_FORWARD][k][15]*O[15] + add)>>shift;
    }

    src += 32;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    Int size = sizeof(TCoeff)*iSkipLine;
    for( j = 0; j < 32; j++ )
    { 
      memset( dst, 0, size );
      dst += line;
    }
  }
#endif
}

/** 32x32 inverse transform implemented using partial butterfly structure (1D)
 *  \param src   input data (transform coefficients)
 *  \param dst   output data (residual)
 *  \param shift specifies right shift after 1D transform
 *  \param line
 *  \param outputMinimum  minimum for clipping
 *  \param outputMaximum  maximum for clipping
 */
#if JVET_D0077_TRANSFORM_OPT
Void partialButterflyInverse32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
Void partialButterflyInverse32(TCoeff *src, TCoeff *dst, Int shift, Int line, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  TCoeff E[16];
  TCoeff EE[8];
  TCoeff EEE[4];
  TCoeff EEEE[2];
#if _WIN64 && _MSC_VER == 1900
  // Workaround for avoid bug in Visual Studio 2015 Update 3 optimizer 
  // https://hevc.hhi.fraunhofer.de/trac/hevc/ticket/1463
  volatile
#endif
  TCoeff O[16], EO[8], EEO[4], EEEO[2];
  TCoeff add = (shift > 0) ? (1<<(shift-1)) : 0;
#if JVET_D0077_TRANSFORM_OPT
  for (j=0; j<line-iSkipLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<16;k++)
    {
      O[k] = g_aiT32[TRANSFORM_INVERSE][ 1][k]*src[ line    ] + g_aiT32[TRANSFORM_INVERSE][ 3][k]*src[ 3*line  ] +
             g_aiT32[TRANSFORM_INVERSE][ 5][k]*src[ 5*line  ] + g_aiT32[TRANSFORM_INVERSE][ 7][k]*src[ 7*line  ] +
             g_aiT32[TRANSFORM_INVERSE][ 9][k]*src[ 9*line  ] + g_aiT32[TRANSFORM_INVERSE][11][k]*src[ 11*line ] +
             g_aiT32[TRANSFORM_INVERSE][13][k]*src[ 13*line ] + g_aiT32[TRANSFORM_INVERSE][15][k]*src[ 15*line ] +
             g_aiT32[TRANSFORM_INVERSE][17][k]*src[ 17*line ] + g_aiT32[TRANSFORM_INVERSE][19][k]*src[ 19*line ] +
             g_aiT32[TRANSFORM_INVERSE][21][k]*src[ 21*line ] + g_aiT32[TRANSFORM_INVERSE][23][k]*src[ 23*line ] +
             g_aiT32[TRANSFORM_INVERSE][25][k]*src[ 25*line ] + g_aiT32[TRANSFORM_INVERSE][27][k]*src[ 27*line ] +
             g_aiT32[TRANSFORM_INVERSE][29][k]*src[ 29*line ] + g_aiT32[TRANSFORM_INVERSE][31][k]*src[ 31*line ];
    }
    for (k=0;k<8;k++)
    {
      EO[k] = g_aiT32[TRANSFORM_INVERSE][ 2][k]*src[ 2*line  ] + g_aiT32[TRANSFORM_INVERSE][ 6][k]*src[ 6*line  ] +
              g_aiT32[TRANSFORM_INVERSE][10][k]*src[ 10*line ] + g_aiT32[TRANSFORM_INVERSE][14][k]*src[ 14*line ] +
              g_aiT32[TRANSFORM_INVERSE][18][k]*src[ 18*line ] + g_aiT32[TRANSFORM_INVERSE][22][k]*src[ 22*line ] +
              g_aiT32[TRANSFORM_INVERSE][26][k]*src[ 26*line ] + g_aiT32[TRANSFORM_INVERSE][30][k]*src[ 30*line ];
    }
    for (k=0;k<4;k++)
    {
      EEO[k] = g_aiT32[TRANSFORM_INVERSE][ 4][k]*src[  4*line ] + g_aiT32[TRANSFORM_INVERSE][12][k]*src[ 12*line ] +
               g_aiT32[TRANSFORM_INVERSE][20][k]*src[ 20*line ] + g_aiT32[TRANSFORM_INVERSE][28][k]*src[ 28*line ];
    }
    EEEO[0] = g_aiT32[TRANSFORM_INVERSE][8][0]*src[ 8*line ] + g_aiT32[TRANSFORM_INVERSE][24][0]*src[ 24*line ];
    EEEO[1] = g_aiT32[TRANSFORM_INVERSE][8][1]*src[ 8*line ] + g_aiT32[TRANSFORM_INVERSE][24][1]*src[ 24*line ];
    EEEE[0] = g_aiT32[TRANSFORM_INVERSE][0][0]*src[ 0      ] + g_aiT32[TRANSFORM_INVERSE][16][0]*src[ 16*line ];
    EEEE[1] = g_aiT32[TRANSFORM_INVERSE][0][1]*src[ 0      ] + g_aiT32[TRANSFORM_INVERSE][16][1]*src[ 16*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    EEE[0] = EEEE[0] + EEEO[0];
    EEE[3] = EEEE[0] - EEEO[0];
    EEE[1] = EEEE[1] + EEEO[1];
    EEE[2] = EEEE[1] - EEEO[1];
    for (k=0;k<4;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+4] = EEE[3-k] - EEO[3-k];
    }
    for (k=0;k<8;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+8] = EE[7-k] - EO[7-k];
    }
    for (k=0;k<16;k++)
    {
      dst[k]    = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[k+16] = Clip3( outputMinimum, outputMaximum, (E[15-k] - O[15-k] + add)>>shift );
    }
    src ++;
    dst += 32;
  }
#if JVET_D0077_TRANSFORM_OPT
  memset(dst, 0, (iSkipLine<<5)*sizeof(TCoeff));
#endif
}

#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT2_B2(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)
#else
void fastForwardDCT2_B2(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use)
#endif
{
  Int j;
  Int E,O;
  Int add = 1<<(shift-1);
  const TMatrixCoeff *iT = g_aiTr2[DCT2][0];

#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoef = dst;
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* E and O */
    E = src[0] + src[1];
    O = src[0] - src[1];

    dst[0] = (iT[0]*E + add)>>shift;
    dst[line] = (iT[2]*O + add)>>shift;

    //assert(iT[2]==256 && iT[0]==256); 

    src += 2;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    dst = pCoef + reducedLine;
    for (j=0; j<2; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT2_B2(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
void fastInverseDCT2_B2(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j;
  Int E,O;
  Int add = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr2[DCT2][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */    
    E = iT[0]*(src[0] + src[line]);
    O = iT[2]*(src[0] - src[line]);

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    dst[0] = Clip3( outputMinimum, outputMaximum, (E + add)>>shift );
    dst[1] = Clip3( outputMinimum, outputMaximum, (O + add)>>shift );

    src   ++;
    dst += 2;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(dst, 0, (iSkipLine<<1)*sizeof(TCoeff));
  }
#endif
}
#endif

/** 4x4 forward transform implemented using partial butterfly structure (1D)
 *  \param src   input data (residual)
 *  \param dst   output data (transform coefficients)
 *  \param shift specifies right shift after 1D transform
 */
// ********************************** DCT-II **********************************
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT2_B4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)
#else
void fastForwardDCT2_B4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use)
#endif
{
  Int j;
  Int E[2],O[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr4[DCT2][0] : g_aiT4[0][0];
#else
  const Short *iT = g_aiT4[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoef = dst;
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* E and O */
    E[0] = src[0] + src[3];
    O[0] = src[0] - src[3];
    E[1] = src[1] + src[2];
    O[1] = src[1] - src[2];

    dst[0] = (iT[0]*E[0] + iT[1]*E[1] + add)>>shift;
    dst[2*line] = (iT[8]*E[0] + iT[9]*E[1] + add)>>shift;
    dst[line] = (iT[4]*O[0] + iT[5]*O[1] + add)>>shift;
    dst[3*line] = (iT[12]*O[0] + iT[13]*O[1] + add)>>shift;

    src += 4;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    dst = pCoef + reducedLine;
    for (j=0; j<4; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT2_B4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
void fastInverseDCT2_B4(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j;
  Int E[2],O[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr4[DCT2][0] : g_aiT4[0][0];
#else
  const Short *iT = g_aiT4[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */    
    O[0] = iT[1*4+0]*src[line] + iT[3*4+0]*src[3*line];
    O[1] = iT[1*4+1]*src[line] + iT[3*4+1]*src[3*line];
    E[0] = iT[0*4+0]*src[0] + iT[2*4+0]*src[2*line];
    E[1] = iT[0*4+1]*src[0] + iT[2*4+1]*src[2*line];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    dst[0] = Clip3( outputMinimum, outputMaximum, (E[0] + O[0] + add)>>shift );
    dst[1] = Clip3( outputMinimum, outputMaximum, (E[1] + O[1] + add)>>shift );
    dst[2] = Clip3( outputMinimum, outputMaximum, (E[1] - O[1] + add)>>shift );
    dst[3] = Clip3( outputMinimum, outputMaximum, (E[0] - O[0] + add)>>shift );
            
    src   ++;
    dst += 4;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(dst, 0, (iSkipLine<<2)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT2_B8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)
#else
void fastForwardDCT2_B8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use)
#endif
{
  Int j,k;
  Int E[4],O[4];
  Int EE[2],EO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr8[DCT2][0] : g_aiT8[0][0];
#else
  const Short *iT = g_aiT8[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoef = dst;
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {  
    /* E and O*/
    for (k=0;k<4;k++)
    {
      E[k] = src[k] + src[7-k];
      O[k] = src[k] - src[7-k];
    }    
    /* EE and EO */
    EE[0] = E[0] + E[3];    
    EO[0] = E[0] - E[3];
    EE[1] = E[1] + E[2];
    EO[1] = E[1] - E[2];

    dst[0] = (iT[0]*EE[0] + iT[1]*EE[1] + add)>>shift;
    dst[4*line] = (iT[32]*EE[0] + iT[33]*EE[1] + add)>>shift; 
    dst[2*line] = (iT[16]*EO[0] + iT[17]*EO[1] + add)>>shift;
    dst[6*line] = (iT[48]*EO[0] + iT[49]*EO[1] + add)>>shift; 

    dst[line] = (iT[8]*O[0] + iT[9]*O[1] + iT[10]*O[2] + iT[11]*O[3] + add)>>shift;
    dst[3*line] = (iT[24]*O[0] + iT[25]*O[1] + iT[26]*O[2] + iT[27]*O[3] + add)>>shift;
    dst[5*line] = (iT[40]*O[0] + iT[41]*O[1] + iT[42]*O[2] + iT[43]*O[3] + add)>>shift;
    dst[7*line] = (iT[56]*O[0] + iT[57]*O[1] + iT[58]*O[2] + iT[59]*O[3] + add)>>shift;

    src += 8;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    dst = pCoef + reducedLine;
    for (j=0; j<8; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT2_B8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
void fastInverseDCT2_B8(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  Int E[4],O[4];
  Int EE[2],EO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr8[DCT2][0] : g_aiT8[0][0];
#else
  const Short *iT = g_aiT8[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++) 
#endif
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<4;k++)
    {
      O[k] = iT[ 1*8+k]*src[line] + iT[ 3*8+k]*src[3*line] + iT[ 5*8+k]*src[5*line] + iT[ 7*8+k]*src[7*line];
    }

    EO[0] = iT[2*8+0]*src[ 2*line ] + iT[6*8+0]*src[ 6*line ];
    EO[1] = iT[2*8+1]*src[ 2*line ] + iT[6*8+1]*src[ 6*line ];
    EE[0] = iT[0*8+0]*src[ 0      ] + iT[4*8+0]*src[ 4*line ];
    EE[1] = iT[0*8+1]*src[ 0      ] + iT[4*8+1]*src[ 4*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    E[0] = EE[0] + EO[0];
    E[3] = EE[0] - EO[0];
    E[1] = EE[1] + EO[1];
    E[2] = EE[1] - EO[1];
    for (k=0;k<4;k++)
    {
      dst[ k   ] = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[ k+4 ] = Clip3( outputMinimum, outputMaximum, (E[3-k] - O[3-k] + add)>>shift );
    }   
    src ++;
    dst += 8;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(dst, 0, (iSkipLine<<3)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT2_B16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)
#else
void fastForwardDCT2_B16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use)
#endif
{
  Int j,k;
  Int E[8],O[8];
  Int EE[4],EO[4];
  Int EEE[2],EEO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr16[DCT2][0] : g_aiT16[0][0];
#else
  const Short *iT = g_aiT16[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoef = dst;
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++) 
#endif
  {    
    /* E and O*/
    for (k=0;k<8;k++)
    {
      E[k] = src[k] + src[15-k];
      O[k] = src[k] - src[15-k];
    } 
    /* EE and EO */
    for (k=0;k<4;k++)
    {
      EE[k] = E[k] + E[7-k];
      EO[k] = E[k] - E[7-k];
    }
    /* EEE and EEO */
    EEE[0] = EE[0] + EE[3];    
    EEO[0] = EE[0] - EE[3];
    EEE[1] = EE[1] + EE[2];
    EEO[1] = EE[1] - EE[2];

    dst[ 0      ] = (iT[0]*EEE[0] + iT[1]*EEE[1] + add)>>shift;        
    dst[ 8*line ] = (iT[8*16]*EEE[0] + iT[8*16+1]*EEE[1] + add)>>shift;    
    dst[ 4*line ] = (iT[4*16]*EEO[0] + iT[4*16+1]*EEO[1] + add)>>shift;        
    dst[ 12*line] = (iT[12*16]*EEO[0] + iT[12*16+1]*EEO[1] + add)>>shift;

    for (k=2;k<16;k+=4)
    {
      dst[ k*line ] = (iT[k*16]*EO[0] + iT[k*16+1]*EO[1] + iT[k*16+2]*EO[2] + iT[k*16+3]*EO[3] + add)>>shift;      
    }

    for (k=1;k<16;k+=2)
    {
      dst[ k*line ] = (iT[k*16]*O[0] + iT[k*16+1]*O[1] + iT[k*16+2]*O[2] + iT[k*16+3]*O[3] + 
        iT[k*16+4]*O[4] + iT[k*16+5]*O[5] + iT[k*16+6]*O[6] + iT[k*16+7]*O[7] + add)>>shift;
    }

    src += 16;
    dst ++; 

  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    dst = pCoef + reducedLine;
    for (j=0; j<16; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT2_B16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
void fastInverseDCT2_B16(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  Int E[8],O[8];
  Int EE[4],EO[4];
  Int EEE[2],EEO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr16[DCT2][0] : g_aiT16[0][0];
#else
  const Short *iT = g_aiT16[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<8;k++)
    {
      O[k] = iT[ 1*16+k]*src[ line] + iT[ 3*16+k]*src[ 3*line] + iT[ 5*16+k]*src[ 5*line] + iT[ 7*16+k]*src[ 7*line] + 
        iT[ 9*16+k]*src[ 9*line] + iT[11*16+k]*src[11*line] + iT[13*16+k]*src[13*line] + iT[15*16+k]*src[15*line];
    }
    for (k=0;k<4;k++)
    {
      EO[k] = iT[ 2*16+k]*src[ 2*line] + iT[ 6*16+k]*src[ 6*line] + iT[10*16+k]*src[10*line] + iT[14*16+k]*src[14*line];
    }
    EEO[0] = iT[4*16]*src[ 4*line ] + iT[12*16]*src[ 12*line ];
    EEE[0] = iT[0]*src[ 0      ] + iT[ 8*16]*src[ 8*line  ];
    EEO[1] = iT[4*16+1]*src[ 4*line ] + iT[12*16+1]*src[ 12*line ];
    EEE[1] = iT[0*16+1]*src[ 0      ] + iT[ 8*16+1]*src[ 8*line  ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    for (k=0;k<2;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }    
    for (k=0;k<4;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      dst[k]   = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[k+8] = Clip3( outputMinimum, outputMaximum, (E[7-k] - O[7-k] + add)>>shift );
    }   
    src ++; 
    dst += 16;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(dst, 0, (iSkipLine<<4)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT2_B32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)
#else
void fastForwardDCT2_B32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use)
#endif
{
  Int j,k;
  Int E[16],O[16];
  Int EE[8],EO[8];
  Int EEE[4],EEO[4];
  Int EEEE[2],EEEO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr32[DCT2][0] : g_aiT32[0][0];
#else
  const Short *iT = g_aiT32[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoef = dst;
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* E and O*/
    for (k=0;k<16;k++)
    {
      E[k] = src[k] + src[31-k];
      O[k] = src[k] - src[31-k];
    } 
    /* EE and EO */
    for (k=0;k<8;k++)
    {
      EE[k] = E[k] + E[15-k];
      EO[k] = E[k] - E[15-k];
    }
    /* EEE and EEO */
    for (k=0;k<4;k++)
    {
      EEE[k] = EE[k] + EE[7-k];
      EEO[k] = EE[k] - EE[7-k];
    }
    /* EEEE and EEEO */
    EEEE[0] = EEE[0] + EEE[3];    
    EEEO[0] = EEE[0] - EEE[3];
    EEEE[1] = EEE[1] + EEE[2];
    EEEO[1] = EEE[1] - EEE[2];

    dst[ 0       ] = (iT[ 0*32+0]*EEEE[0] + iT[ 0*32+1]*EEEE[1] + add)>>shift;
    dst[ 16*line ] = (iT[16*32+0]*EEEE[0] + iT[16*32+1]*EEEE[1] + add)>>shift;
    dst[ 8*line  ] = (iT[ 8*32+0]*EEEO[0] + iT[ 8*32+1]*EEEO[1] + add)>>shift; 
    dst[ 24*line ] = (iT[24*32+0]*EEEO[0] + iT[24*32+1]*EEEO[1] + add)>>shift;
    for (k=4;k<32;k+=8)
    {
      dst[ k*line ] = (iT[k*32+0]*EEO[0] + iT[k*32+1]*EEO[1] + iT[k*32+2]*EEO[2] + iT[k*32+3]*EEO[3] + add)>>shift;
    }       
    for (k=2;k<32;k+=4)
    {
      dst[ k*line ] = (iT[k*32+0]*EO[0] + iT[k*32+1]*EO[1] + iT[k*32+2]*EO[2] + iT[k*32+3]*EO[3] + 
        iT[k*32+4]*EO[4] + iT[k*32+5]*EO[5] + iT[k*32+6]*EO[6] + iT[k*32+7]*EO[7] + add)>>shift;
    }       
    for (k=1;k<32;k+=2)
    {
      dst[ k*line ] = (iT[k*32+ 0]*O[ 0] + iT[k*32+ 1]*O[ 1] + iT[k*32+ 2]*O[ 2] + iT[k*32+ 3]*O[ 3] + 
        iT[k*32+ 4]*O[ 4] + iT[k*32+ 5]*O[ 5] + iT[k*32+ 6]*O[ 6] + iT[k*32+ 7]*O[ 7] +
        iT[k*32+ 8]*O[ 8] + iT[k*32+ 9]*O[ 9] + iT[k*32+10]*O[10] + iT[k*32+11]*O[11] + 
        iT[k*32+12]*O[12] + iT[k*32+13]*O[13] + iT[k*32+14]*O[14] + iT[k*32+15]*O[15] + add)>>shift;
    }
    src += 32;
    dst ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    dst = pCoef + reducedLine;
    for (j=0; j<32; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT2_B32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#else
void fastInverseDCT2_B32(TCoeff *src, TCoeff *dst, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)
#endif
{
  Int j,k;
  Int E[16],O[16];
  Int EE[8],EO[8];
  Int EEE[4],EEO[4];
  Int EEEE[2],EEEO[2];
  Int add = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr32[DCT2][0] : g_aiT32[0][0];
#else
  const Short *iT = g_aiT32[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (j=0; j<reducedLine; j++)
#else
  for (j=0; j<line; j++)
#endif
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<16;k++)
    {
      O[k] = iT[ 1*32+k]*src[ line  ] + iT[ 3*32+k]*src[ 3*line  ] + iT[ 5*32+k]*src[ 5*line  ] + iT[ 7*32+k]*src[ 7*line  ] + 
        iT[ 9*32+k]*src[ 9*line  ] + iT[11*32+k]*src[ 11*line ] + iT[13*32+k]*src[ 13*line ] + iT[15*32+k]*src[ 15*line ] + 
        iT[17*32+k]*src[ 17*line ] + iT[19*32+k]*src[ 19*line ] + iT[21*32+k]*src[ 21*line ] + iT[23*32+k]*src[ 23*line ] + 
        iT[25*32+k]*src[ 25*line ] + iT[27*32+k]*src[ 27*line ] + iT[29*32+k]*src[ 29*line ] + iT[31*32+k]*src[ 31*line ];
    }
    for (k=0;k<8;k++)
    {
      EO[k] = iT[ 2*32+k]*src[ 2*line  ] + iT[ 6*32+k]*src[ 6*line  ] + iT[10*32+k]*src[ 10*line ] + iT[14*32+k]*src[ 14*line ] + 
        iT[18*32+k]*src[ 18*line ] + iT[22*32+k]*src[ 22*line ] + iT[26*32+k]*src[ 26*line ] + iT[30*32+k]*src[ 30*line ];
    }
    for (k=0;k<4;k++)
    {
      EEO[k] = iT[4*32+k]*src[ 4*line ] + iT[12*32+k]*src[ 12*line ] + iT[20*32+k]*src[ 20*line ] + iT[28*32+k]*src[ 28*line ];
    }
    EEEO[0] = iT[8*32+0]*src[ 8*line ] + iT[24*32+0]*src[ 24*line ];
    EEEO[1] = iT[8*32+1]*src[ 8*line ] + iT[24*32+1]*src[ 24*line ];
    EEEE[0] = iT[0*32+0]*src[ 0      ] + iT[16*32+0]*src[ 16*line ];    
    EEEE[1] = iT[0*32+1]*src[ 0      ] + iT[16*32+1]*src[ 16*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    EEE[0] = EEEE[0] + EEEO[0];
    EEE[3] = EEEE[0] - EEEO[0];
    EEE[1] = EEEE[1] + EEEO[1];
    EEE[2] = EEEE[1] - EEEO[1];    
    for (k=0;k<4;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+4] = EEE[3-k] - EEO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+8] = EE[7-k] - EO[7-k];
    }    
    for (k=0;k<16;k++)
    {
      dst[k]    = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + add)>>shift );
      dst[k+16] = Clip3( outputMinimum, outputMaximum, (E[15-k] - O[15-k] + add)>>shift );
    }
    src ++;
    dst += 32;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(dst, 0, (iSkipLine<<5)*sizeof(TCoeff));
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastForwardDCT2_B64(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                         , Int use
#endif
                         )
{
  Int rnd_factor = 1<<(shift-1);
  const Int uiTrSize = 64;
#if COM16_C806_T64
  const TMatrixCoeff *iT = g_aiTr64[DCT2][0];
#else
  const TMatrixCoeff *iT = NULL;
  assert(0);
#endif

  Int j, k;
  Int E[32],O[32];
  Int EE[16],EO[16];
  Int EEE[8],EEO[8];
  Int EEEE[4],EEEO[4];
  Int EEEEE[2],EEEEO[2];
  TCoeff *tmp = dst;

  Bool zo = iSkipLine2 >= 32;
  for (j=0; j<line-iSkipLine; j++)
  {    
    /* E and O*/
    for (k=0;k<32;k++)
    {
      E[k] = src[k] + src[63-k];
      O[k] = src[k] - src[63-k];
    } 
    /* EE and EO */
    for (k=0;k<16;k++)
    {
      EE[k] = E[k] + E[31-k];
      EO[k] = E[k] - E[31-k];
    }
    /* EEE and EEO */
    for (k=0;k<8;k++)
    {
      EEE[k] = EE[k] + EE[15-k];
      EEO[k] = EE[k] - EE[15-k];
    }
    /* EEEE and EEEO */
    for (k=0;k<4;k++)
    {
      EEEE[k] = EEE[k] + EEE[7-k];
      EEEO[k] = EEE[k] - EEE[7-k];
    }
    /* EEEEE and EEEEO */
    EEEEE[0] = EEEE[0] + EEEE[3];    
    EEEEO[0] = EEEE[0] - EEEE[3];
    EEEEE[1] = EEEE[1] + EEEE[2];
    EEEEO[1] = EEEE[1] - EEEE[2];

    dst[ 0       ] = (iT[ 0*64+0]*EEEEE[0] + iT[ 0*64+1]*EEEEE[1] + rnd_factor)>>shift;
    dst[ 16*line ] = (iT[16*64+0]*EEEEO[0] + iT[16*64+1]*EEEEO[1] + rnd_factor)>>shift; 
    
    if( !zo )
    {
      dst[ 32*line ] = (iT[32*64+0]*EEEEE[0] + iT[32*64+1]*EEEEE[1] + rnd_factor)>>shift;
      dst[ 48*line ] = (iT[48*64+0]*EEEEO[0] + iT[48*64+1]*EEEEO[1] + rnd_factor)>>shift;
    }
    for (k=8;k<(zo?32:64);k+=16)
    {
      dst[ k*line ] = (iT[k*64+0]*EEEO[0] + iT[k*64+1]*EEEO[1] + iT[k*64+2]*EEEO[2] + iT[k*64+3]*EEEO[3] + rnd_factor)>>shift;
    }       
    for (k=4;k<(zo?32:64);k+=8)
    {
      dst[ k*line ] = (iT[k*64+0]*EEO[0] + iT[k*64+1]*EEO[1] + iT[k*64+2]*EEO[2] + iT[k*64+3]*EEO[3] + 
        iT[k*64+4]*EEO[4] + iT[k*64+5]*EEO[5] + iT[k*64+6]*EEO[6] + iT[k*64+7]*EEO[7] + rnd_factor)>>shift;
    }       
    for (k=2;k<(zo?32:64);k+=4)
    {
      dst[ k*line ] = (iT[k*64+ 0]*EO[ 0] + iT[k*64+ 1]*EO[ 1] + iT[k*64+ 2]*EO[ 2] + iT[k*64+ 3]*EO[ 3] + 
        iT[k*64+ 4]*EO[ 4] + iT[k*64+ 5]*EO[ 5] + iT[k*64+ 6]*EO[ 6] + iT[k*64+ 7]*EO[ 7] +
        iT[k*64+ 8]*EO[ 8] + iT[k*64+ 9]*EO[ 9] + iT[k*64+10]*EO[10] + iT[k*64+11]*EO[11] + 
        iT[k*64+12]*EO[12] + iT[k*64+13]*EO[13] + iT[k*64+14]*EO[14] + iT[k*64+15]*EO[15] + rnd_factor)>>shift;
    }
    for (k=1;k<(zo?32:64);k+=2)
    {
      dst[ k*line ] = (iT[k*64+ 0]*O[ 0] + iT[k*64+ 1]*O[ 1] + iT[k*64+ 2]*O[ 2] + iT[k*64+ 3]*O[ 3] + 
        iT[k*64+ 4]*O[ 4] + iT[k*64+ 5]*O[ 5] + iT[k*64+ 6]*O[ 6] + iT[k*64+ 7]*O[ 7] +
        iT[k*64+ 8]*O[ 8] + iT[k*64+ 9]*O[ 9] + iT[k*64+10]*O[10] + iT[k*64+11]*O[11] + 
        iT[k*64+12]*O[12] + iT[k*64+13]*O[13] + iT[k*64+14]*O[14] + iT[k*64+15]*O[15] + 
        iT[k*64+16]*O[16] + iT[k*64+17]*O[17] + iT[k*64+18]*O[18] + iT[k*64+19]*O[19] + 
        iT[k*64+20]*O[20] + iT[k*64+21]*O[21] + iT[k*64+22]*O[22] + iT[k*64+23]*O[23] + 
        iT[k*64+24]*O[24] + iT[k*64+25]*O[25] + iT[k*64+26]*O[26] + iT[k*64+27]*O[27] + 
        iT[k*64+28]*O[28] + iT[k*64+29]*O[29] + iT[k*64+30]*O[30] + iT[k*64+31]*O[31] + rnd_factor)>>shift;
    }
    src += uiTrSize;
    dst ++;
  }

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  if( iSkipLine )
  {
    dst= tmp + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
    }
  }
  if( iSkipLine2 )
  {
    dst = tmp + line*cutoff;
    memset( dst, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}


#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastInverseDCT2_B64(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2,  
#if JVET_D0077_TRANSFORM_OPT
                         Int use,
#endif
                         const TCoeff outputMinimum, const TCoeff outputMaximum)
{
  Int rnd_factor = 1<<(shift-1);
  const Int uiTrSize = 64;
#if COM16_C806_T64
  const TMatrixCoeff *iT = g_aiTr64[DCT2][0];
#else
  const TMatrixCoeff *iT = NULL;
  assert(0);
#endif

  Int j, k;
  Int E[32],O[32];
  Int EE[16],EO[16];
  Int EEE[8],EEO[8];
  Int EEEE[4],EEEO[4];
  Int EEEEE[2],EEEEO[2];
  Bool zo = iSkipLine2 >= 32;
  for (j=0; j<line-iSkipLine; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<32;k++)
    {
      O[k] = iT[ 1*64+k]*coeff[ line  ] + iT[ 3*64+k]*coeff[ 3*line  ] + iT[ 5*64+k]*coeff[ 5*line  ] + iT[ 7*64+k]*coeff[ 7*line  ] + 
        iT[ 9*64+k]*coeff[ 9*line  ] + iT[11*64+k]*coeff[ 11*line ] + iT[13*64+k]*coeff[ 13*line ] + iT[15*64+k]*coeff[ 15*line ] + 
        iT[17*64+k]*coeff[ 17*line ] + iT[19*64+k]*coeff[ 19*line ] + iT[21*64+k]*coeff[ 21*line ] + iT[23*64+k]*coeff[ 23*line ] + 
        iT[25*64+k]*coeff[ 25*line ] + iT[27*64+k]*coeff[ 27*line ] + iT[29*64+k]*coeff[ 29*line ] + iT[31*64+k]*coeff[ 31*line ] +
        ( zo ? 0 : (
        iT[33*64+k]*coeff[ 33*line ] + iT[35*64+k]*coeff[ 35*line ] + iT[37*64+k]*coeff[ 37*line ] + iT[39*64+k]*coeff[ 39*line ] +
        iT[41*64+k]*coeff[ 41*line ] + iT[43*64+k]*coeff[ 43*line ] + iT[45*64+k]*coeff[ 45*line ] + iT[47*64+k]*coeff[ 47*line ] +
        iT[49*64+k]*coeff[ 49*line ] + iT[51*64+k]*coeff[ 51*line ] + iT[53*64+k]*coeff[ 53*line ] + iT[55*64+k]*coeff[ 55*line ] +
        iT[57*64+k]*coeff[ 57*line ] + iT[59*64+k]*coeff[ 59*line ] + iT[61*64+k]*coeff[ 61*line ] + iT[63*64+k]*coeff[ 63*line ] ) );
    }
    for (k=0;k<16;k++)
    {
      EO[k] = iT[ 2*64+k]*coeff[ 2*line  ] + iT[ 6*64+k]*coeff[ 6*line  ] + iT[10*64+k]*coeff[ 10*line ] + iT[14*64+k]*coeff[ 14*line ] + 
        iT[18*64+k]*coeff[ 18*line ] + iT[22*64+k]*coeff[ 22*line ] + iT[26*64+k]*coeff[ 26*line ] + iT[30*64+k]*coeff[ 30*line ] + 
        ( zo ? 0 : (
        iT[34*64+k]*coeff[ 34*line ] + iT[38*64+k]*coeff[ 38*line ] + iT[42*64+k]*coeff[ 42*line ] + iT[46*64+k]*coeff[ 46*line ] +
        iT[50*64+k]*coeff[ 50*line ] + iT[54*64+k]*coeff[ 54*line ] + iT[58*64+k]*coeff[ 58*line ] + iT[62*64+k]*coeff[ 62*line ] ) );
    }
    for (k=0;k<8;k++)
    {
      EEO[k] = iT[4*64+k]*coeff[ 4*line ] + iT[12*64+k]*coeff[ 12*line ] + iT[20*64+k]*coeff[ 20*line ] + iT[28*64+k]*coeff[ 28*line ] +
        ( zo ? 0 : (
        iT[36*64+k]*coeff[ 36*line ] + iT[44*64+k]*coeff[ 44*line ] + iT[52*64+k]*coeff[ 52*line ] + iT[60*64+k]*coeff[ 60*line ] ) );
    }
    for (k=0;k<4;k++)
    {
      EEEO[k] = iT[8*64+k]*coeff[ 8*line ] + iT[24*64+k]*coeff[ 24*line ] + ( zo ? 0 : ( iT[40*64+k]*coeff[ 40*line ] + iT[56*64+k]*coeff[ 56*line ] ) );
    }
    EEEEO[0] = iT[16*64+0]*coeff[ 16*line ] + ( zo ? 0 : iT[48*64+0]*coeff[ 48*line ] );
    EEEEO[1] = iT[16*64+1]*coeff[ 16*line ] + ( zo ? 0 : iT[48*64+1]*coeff[ 48*line ] );
    EEEEE[0] = iT[ 0*64+0]*coeff[  0      ] + ( zo ? 0 : iT[32*64+0]*coeff[ 32*line ] );    
    EEEEE[1] = iT[ 0*64+1]*coeff[  0      ] + ( zo ? 0 : iT[32*64+1]*coeff[ 32*line ] );

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    for (k=0;k<2;k++)
    {
      EEEE[k] = EEEEE[k] + EEEEO[k];
      EEEE[k+2] = EEEEE[1-k] - EEEEO[1-k];
    } 
    for (k=0;k<4;k++)
    {
      EEE[k] = EEEE[k] + EEEO[k];
      EEE[k+4] = EEEE[3-k] - EEEO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+8] = EEE[7-k] - EEO[7-k];
    }   
    for (k=0;k<16;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+16] = EE[15-k] - EO[15-k];
    }    
    for (k=0;k<32;k++)
    {
      block[k]    = Clip3( outputMinimum, outputMaximum, (E[k] + O[k] + rnd_factor)>>shift );
      block[k+32] = Clip3( outputMinimum, outputMaximum, (E[31-k] - O[31-k] + rnd_factor)>>shift );
    }
    coeff ++;
    block += uiTrSize;
  }
  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
}



#if JVET_C0024_QTBT
void fastForwardDCT2_B128(TCoeff *src, TCoeff *dst, Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                          , Int use
#endif
                          )
{
  Int j,k;
  Int E[64],O[64];
  Int EE[32],EO[32];
  Int EEE[16],EEO[16];
  Int EEEE[8],EEEO[8];
  Int EEEEE[4], EEEEO[4];
  Int EEEEEE[2], EEEEEO[2];
  Int add = 1<<(shift-1);
  const TMatrixCoeff (*iT)[128] = g_aiTr128[DCT2];

  TCoeff* tmp = dst;
  for (j=0; j<line-iSkipLine; j++)
  {    
    /* E and O*/
    for (k=0;k< 64;k++)
    {
      E[k] = src[k] + src[127-k];
      O[k] = src[k] - src[127-k];
    } 
    /* EE and EO */
    for (k=0;k< 32;k++)
    {
      EE[k] = E[k] + E[63-k];
      EO[k] = E[k] - E[63-k];
    }

    /* EEE and EEO */
    for (k=0;k< 16;k++)
    {
      EEE[k] = EE[k] + EE[31-k];
      EEO[k] = EE[k] - EE[31-k];
    }

    /* EEEE and EEEO */
    for( k=0; k< 8; k++)
    {
      EEEE[k]= EEE[k]+ EEE[15-k];
      EEEO[k]= EEE[k]- EEE[15-k];
    }

    for( k=0; k< 4; k++)
    {
      EEEEE[k] = EEEE[k]+ EEEE[7-k];
      EEEEO[k] = EEEE[k]- EEEE[7-k];
    }

    for( k=0; k< 2; k++)
    {
      EEEEEE[k] = EEEEE[k]+ EEEEE[3-k];
      EEEEEO[k] = EEEEE[k]- EEEEE[3-k];
    }

    //0
    dst[ 0       ] = (iT[0][0]*EEEEEE[0] 
    + iT[0][ 1]*EEEEEE[1] 
    + add)>>shift;
    dst[ 64*line ] = (iT[64][0]*EEEEEE[0] 
    + iT[64][1]*EEEEEE[1] 
    + add)>>shift;

    //2
    for (k=32;k<128;k+=64)
    {
      dst[ k*line ] = ( iT[k][0]*EEEEEO[0] 
      + iT[k][1]*EEEEEO[1] 
      + add)>>shift;
    }

    //4
    for (k=16;k<128;k+=32) 
    {
      dst[ k*line ] = ( iT[k][0]*EEEEO[0] 
      + iT[k][1]*EEEEO[1] 
      + iT[k][2]*EEEEO[2] 
      + iT[k][3]*EEEEO[3] 
      + add)>>shift;
    }       

    //8
    for (k=8;k<128;k+=16)
    {
      dst[ k*line ] = (iT[k][0]*EEEO[0] 
      + iT[k][1]*EEEO[1] 
      + iT[k][2]*EEEO[2] 
      + iT[k][3]*EEEO[3] 
      + iT[k][4]*EEEO[4] 
      + iT[k][5]*EEEO[5] 
      + iT[k][6]*EEEO[6] 
      + iT[k][7]*EEEO[7] 
      + add)>>shift;
    }  

    //16
    for (k=4;k<128;k+=8)
    {
      dst[ k*line ] = (iT[k][0]*EEO[ 0] 
      + iT[k][1]*EEO[ 1] 
      + iT[k][2]*EEO[ 2] 
      + iT[k][3]*EEO[ 3] 
      + iT[k][4]*EEO[ 4] 
      + iT[k][5]*EEO[ 5] 
      + iT[k][6]*EEO[ 6] 
      + iT[k][7]*EEO[ 7] 
      + iT[k][8]*EEO[ 8] 
      + iT[k][9]*EEO[ 9] 
      + iT[k][10]*EEO[10] 
      + iT[k][11]*EEO[11] 
      + iT[k][12]*EEO[12] 
      + iT[k][13]*EEO[13] 
      + iT[k][14]*EEO[14] 
      + iT[k][15]*EEO[15] 
      + add)>>shift;
    }


    //32
    for (k=2;k<128;k+=4)
    {
      dst[ k*line ] = (iT[k][0]*EO[ 0] 
      + iT[k][1]*EO[ 1] 
      + iT[k][2]*EO[ 2] 
      + iT[k][3]*EO[ 3] 
      + iT[k][4]*EO[ 4] 
      + iT[k][5]*EO[ 5] 
      + iT[k][6]*EO[ 6] 
      + iT[k][7]*EO[ 7] 
      + iT[k][8]*EO[ 8] 
      + iT[k][9]*EO[ 9] 
      + iT[k][10]*EO[10] 
      + iT[k][11]*EO[11] 
      + iT[k][12]*EO[12] 
      + iT[k][13]*EO[13] 
      + iT[k][14]*EO[14] 
      + iT[k][15]*EO[15] 
      + iT[k][16]*EO[16] 
      + iT[k][17]*EO[17] 
      + iT[k][18]*EO[18] 
      + iT[k][19]*EO[19] 
      + iT[k][20]*EO[20] 
      + iT[k][21]*EO[21] 
      + iT[k][22]*EO[22] 
      + iT[k][23]*EO[23] 
      + iT[k][24]*EO[24] 
      + iT[k][25]*EO[25] 
      + iT[k][26]*EO[26] 
      + iT[k][27]*EO[27] 
      + iT[k][28]*EO[28] 
      + iT[k][29]*EO[29] 
      + iT[k][30]*EO[30] 
      + iT[k][31]*EO[31] 
      + add)>>shift;
    }

    //64
    for (k=1;k<128;k+=2)
    {
      dst[ k*line ] = (iT[k][ 0]*O[ 0] 
      + iT[k][1]*O[ 1] 
      + iT[k][2]*O[ 2] 
      + iT[k][3]*O[ 3] 
      + iT[k][4]*O[ 4] 
      + iT[k][5]*O[ 5] 
      + iT[k][6]*O[ 6] 
      + iT[k][7]*O[ 7] 
      + iT[k][8]*O[ 8] 
      + iT[k][9]*O[ 9] 
      + iT[k][10]*O[10] 
      + iT[k][11]*O[11] 
      + iT[k][12]*O[12] 
      + iT[k][13]*O[13] 
      + iT[k][14]*O[14] 
      + iT[k][15]*O[15] 
      + iT[k][16]*O[16] 
      + iT[k][17]*O[17] 
      + iT[k][18]*O[18] 
      + iT[k][19]*O[19] 
      + iT[k][20]*O[20] 
      + iT[k][21]*O[21] 
      + iT[k][22]*O[22] 
      + iT[k][23]*O[23] 
      + iT[k][24]*O[24] 
      + iT[k][25]*O[25] 
      + iT[k][26]*O[26] 
      + iT[k][27]*O[27] 
      + iT[k][28]*O[28] 
      + iT[k][29]*O[29] 
      + iT[k][30]*O[30] 
      + iT[k][31]*O[31] 

      + iT[k][32]*O[32] 
      + iT[k][33]*O[33] 
      + iT[k][34]*O[34] 
      + iT[k][35]*O[35] 
      + iT[k][36]*O[36] 
      + iT[k][37]*O[37] 
      + iT[k][38]*O[38] 
      + iT[k][39]*O[39] 
      + iT[k][40]*O[40] 
      + iT[k][41]*O[41] 
      + iT[k][42]*O[42] 
      + iT[k][43]*O[43] 
      + iT[k][44]*O[44] 
      + iT[k][45]*O[45] 
      + iT[k][46]*O[46] 
      + iT[k][47]*O[47] 
      + iT[k][48]*O[48] 
      + iT[k][49]*O[49] 
      + iT[k][50]*O[50] 
      + iT[k][51]*O[51] 
      + iT[k][52]*O[52] 
      + iT[k][53]*O[53] 
      + iT[k][54]*O[54] 
      + iT[k][55]*O[55] 
      + iT[k][56]*O[56] 
      + iT[k][57]*O[57] 
      + iT[k][58]*O[58] 
      + iT[k][59]*O[59] 
      + iT[k][60]*O[60] 
      + iT[k][61]*O[61] 
      + iT[k][62]*O[62] 
      + iT[k][63]*O[63] 
      + add)>>shift;
    }
    src += 128;
    dst ++;
  }
  const UInt uiTrSize = 128;
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  if( iSkipLine )
  {
    dst= tmp + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( dst, 0, sizeof(TCoeff)*iSkipLine );
      dst += line;
}
  }
  if( iSkipLine2 )
  {
    dst = tmp + line*cutoff;
    memset( dst, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}
void fastInverseDCT2_B128(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2,  
#if JVET_D0077_TRANSFORM_OPT 
                          Int use,
#endif
                          const TCoeff outputMinimum, const TCoeff outputMaximum)
{
  Int j,k;
  Int E[64],O[64];
  Int EE[32],EO[32];
  Int EEE[16],EEO[16];
  Int EEEE[8],EEEO[8];
  Int EEEEE[4],EEEEO[4];
  Int EEEEEE[2],EEEEEO[2];
  Int add = 1<<(shift-1);
  const TMatrixCoeff (*iT)[128] = g_aiTr128[DCT2];

  Bool c1 = iSkipLine2 >= 96 ;
  Bool c2 = iSkipLine2 >= 64 ;
  Bool c3 = iSkipLine2 >= 32 ;

  for (j=0; j<line - iSkipLine; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    if (c1)
    {
      for (k=0;k<64;k++) //+2
      {
        O[k] = iT[ 1][k]*coeff[ line  ] 
        + iT[ 3][k]*coeff[ 3*line  ] 
        + iT[ 5][k]*coeff[ 5*line  ] 
        + iT[ 7][k]*coeff[ 7*line  ] 
        + iT[ 9][k]*coeff[ 9*line  ] 
        + iT[11][k]*coeff[ 11*line ] 
        + iT[13][k]*coeff[ 13*line ] 
        + iT[15][k]*coeff[ 15*line ] 
        + iT[17][k]*coeff[ 17*line ] 
        + iT[19][k]*coeff[ 19*line ] 
        + iT[21][k]*coeff[ 21*line ]
        + iT[23][k]*coeff[ 23*line ] 
        + iT[25][k]*coeff[ 25*line ] 
        + iT[27][k]*coeff[ 27*line ] 
        + iT[29][k]*coeff[ 29*line ] 
        + iT[31][k]*coeff[ 31*line ]
        ;
      }

      for (k=0;k<32;k++) //+4
      {
        EO[k] = iT[ 2][k]*coeff[ 2*line  ] 
        + iT[ 6][k]*coeff[ 6*line  ] 
        + iT[10][k]*coeff[ 10*line ] 
        + iT[14][k]*coeff[ 14*line ] 
        + iT[18][k]*coeff[ 18*line ] 
        + iT[22][k]*coeff[ 22*line ] 
        + iT[26][k]*coeff[ 26*line ] 
        + iT[30][k]*coeff[ 30*line ]
        ;
      }
    }
    else if (c2)
    {
      for (k=0;k<64;k++) //+2
      {
        O[k] = iT[ 1][k]*coeff[ line  ] 
        + iT[ 3][k]*coeff[ 3*line  ] 
        + iT[ 5][k]*coeff[ 5*line  ] 
        + iT[ 7][k]*coeff[ 7*line  ] 
        + iT[ 9][k]*coeff[ 9*line  ] 
        + iT[11][k]*coeff[ 11*line ] 
        + iT[13][k]*coeff[ 13*line ] 
        + iT[15][k]*coeff[ 15*line ] 
        + iT[17][k]*coeff[ 17*line ] 
        + iT[19][k]*coeff[ 19*line ] 
        + iT[21][k]*coeff[ 21*line ]
        + iT[23][k]*coeff[ 23*line ] 
        + iT[25][k]*coeff[ 25*line ] 
        + iT[27][k]*coeff[ 27*line ] 
        + iT[29][k]*coeff[ 29*line ] 
        + iT[31][k]*coeff[ 31*line ]
        + iT[33][k]*coeff[ 33*line ] 
        + iT[35][k]*coeff[ 35*line ] 
        + iT[37][k]*coeff[ 37*line ] 
        + iT[39][k]*coeff[ 39*line ]
        + iT[41][k]*coeff[ 41*line ]
        + iT[43][k]*coeff[ 43*line ] 
        + iT[45][k]*coeff[ 45*line ] 
        + iT[47][k]*coeff[ 47*line ] 
        + iT[49][k]*coeff[ 49*line ] 
        + iT[51][k]*coeff[ 51*line ]
        + iT[53][k]*coeff[ 53*line ] 
        + iT[55][k]*coeff[ 55*line ] 
        + iT[57][k]*coeff[ 57*line ] 
        + iT[59][k]*coeff[ 59*line ] 
        + iT[61][k]*coeff[ 61*line ]
        + iT[63][k]*coeff[ 63*line ] 
        ;
      }

      for (k=0;k<32;k++) //+4
      {
        EO[k] = iT[ 2][k]*coeff[ 2*line  ] 
        + iT[ 6][k]*coeff[ 6*line  ] 
        + iT[10][k]*coeff[ 10*line ] 
        + iT[14][k]*coeff[ 14*line ] 
        + iT[18][k]*coeff[ 18*line ] 
        + iT[22][k]*coeff[ 22*line ] 
        + iT[26][k]*coeff[ 26*line ] 
        + iT[30][k]*coeff[ 30*line ]
        + iT[34][k]*coeff[ 34*line ] 
        + iT[38][k]*coeff[ 38*line ] 
        + iT[42][k]*coeff[ 42*line ] 
        + iT[46][k]*coeff[ 46*line ] 
        + iT[50][k]*coeff[ 50*line ]
        + iT[54][k]*coeff[ 54*line ] 
        + iT[58][k]*coeff[ 58*line ] 
        + iT[62][k]*coeff[ 62*line ] 
        ;
      }
    }
    else if (c3)
    {
      for (k=0;k<64;k++) //+2
      {
        O[k] = iT[ 1][k]*coeff[ line  ] 
        + iT[ 3][k]*coeff[ 3*line  ] 
        + iT[ 5][k]*coeff[ 5*line  ] 
        + iT[ 7][k]*coeff[ 7*line  ] 
        + iT[ 9][k]*coeff[ 9*line  ] 
        + iT[11][k]*coeff[ 11*line ] 
        + iT[13][k]*coeff[ 13*line ] 
        + iT[15][k]*coeff[ 15*line ] 
        + iT[17][k]*coeff[ 17*line ] 
        + iT[19][k]*coeff[ 19*line ] 
        + iT[21][k]*coeff[ 21*line ]
        + iT[23][k]*coeff[ 23*line ] 
        + iT[25][k]*coeff[ 25*line ] 
        + iT[27][k]*coeff[ 27*line ] 
        + iT[29][k]*coeff[ 29*line ] 
        + iT[31][k]*coeff[ 31*line ]
        + iT[33][k]*coeff[ 33*line ] 
        + iT[35][k]*coeff[ 35*line ] 
        + iT[37][k]*coeff[ 37*line ] 
        + iT[39][k]*coeff[ 39*line ]
        + iT[41][k]*coeff[ 41*line ]
        + iT[43][k]*coeff[ 43*line ] 
        + iT[45][k]*coeff[ 45*line ] 
        + iT[47][k]*coeff[ 47*line ] 
        + iT[49][k]*coeff[ 49*line ] 
        + iT[51][k]*coeff[ 51*line ]
        + iT[53][k]*coeff[ 53*line ] 
        + iT[55][k]*coeff[ 55*line ] 
        + iT[57][k]*coeff[ 57*line ] 
        + iT[59][k]*coeff[ 59*line ] 
        + iT[61][k]*coeff[ 61*line ]
        + iT[63][k]*coeff[ 63*line ] 
        + iT[65][k]*coeff[ 65*line ] 
        + iT[67][k]*coeff[ 67*line ] 
        + iT[69][k]*coeff[ 69*line ] 
        + iT[71][k]*coeff[ 71*line ] 
        + iT[73][k]*coeff[ 73*line ] 
        + iT[75][k]*coeff[ 75*line ] 
        + iT[77][k]*coeff[ 77*line ] 
        + iT[79][k]*coeff[ 79*line ] 
        + iT[81][k]*coeff[ 81*line ] 
        + iT[83][k]*coeff[ 83*line ] 
        + iT[85][k]*coeff[ 85*line ]
        + iT[87][k]*coeff[ 87*line ] 
        + iT[89][k]*coeff[ 89*line ] 
        + iT[91][k]*coeff[ 91*line ] 
        + iT[93][k]*coeff[ 93*line ] 
        + iT[95][k]*coeff[ 95*line ]
        ;
      }

      for (k=0;k<32;k++) //+4
      {
        EO[k] = iT[ 2][k]*coeff[ 2*line  ] 
        + iT[ 6][k]*coeff[ 6*line  ] 
        + iT[10][k]*coeff[ 10*line ] 
        + iT[14][k]*coeff[ 14*line ] 
        + iT[18][k]*coeff[ 18*line ] 
        + iT[22][k]*coeff[ 22*line ] 
        + iT[26][k]*coeff[ 26*line ] 
        + iT[30][k]*coeff[ 30*line ]
        + iT[34][k]*coeff[ 34*line ] 
        + iT[38][k]*coeff[ 38*line ] 
        + iT[42][k]*coeff[ 42*line ] 
        + iT[46][k]*coeff[ 46*line ] 
        + iT[50][k]*coeff[ 50*line ]
        + iT[54][k]*coeff[ 54*line ] 
        + iT[58][k]*coeff[ 58*line ] 
        + iT[62][k]*coeff[ 62*line ] 
        + iT[66][k]*coeff[ 66*line  ] 
        + iT[70][k]*coeff[ 70*line  ] 
        + iT[74][k]*coeff[ 74*line ] 
        + iT[78][k]*coeff[ 78*line ] 
        + iT[82][k]*coeff[ 82*line ] 
        + iT[86][k]*coeff[ 86*line ] 
        + iT[90][k]*coeff[ 90*line ] 
        + iT[94][k]*coeff[ 94*line ]
        ;
      }
    }
    else
    {
      for (k=0;k<64;k++) //+2
      {
        O[k] = iT[ 1][k]*coeff[ line  ] 
        + iT[ 3][k]*coeff[ 3*line  ] 
        + iT[ 5][k]*coeff[ 5*line  ] 
        + iT[ 7][k]*coeff[ 7*line  ] 
        + iT[ 9][k]*coeff[ 9*line  ] 
        + iT[11][k]*coeff[ 11*line ] 
        + iT[13][k]*coeff[ 13*line ] 
        + iT[15][k]*coeff[ 15*line ] 
        + iT[17][k]*coeff[ 17*line ] 
        + iT[19][k]*coeff[ 19*line ] 
        + iT[21][k]*coeff[ 21*line ]
        + iT[23][k]*coeff[ 23*line ] 
        + iT[25][k]*coeff[ 25*line ] 
        + iT[27][k]*coeff[ 27*line ] 
        + iT[29][k]*coeff[ 29*line ] 
        + iT[31][k]*coeff[ 31*line ]
        + iT[33][k]*coeff[ 33*line ] 
        + iT[35][k]*coeff[ 35*line ] 
        + iT[37][k]*coeff[ 37*line ] 
        + iT[39][k]*coeff[ 39*line ]
        + iT[41][k]*coeff[ 41*line ]
        + iT[43][k]*coeff[ 43*line ] 
        + iT[45][k]*coeff[ 45*line ] 
        + iT[47][k]*coeff[ 47*line ] 
        + iT[49][k]*coeff[ 49*line ] 
        + iT[51][k]*coeff[ 51*line ]
        + iT[53][k]*coeff[ 53*line ] 
        + iT[55][k]*coeff[ 55*line ] 
        + iT[57][k]*coeff[ 57*line ] 
        + iT[59][k]*coeff[ 59*line ] 
        + iT[61][k]*coeff[ 61*line ]
        + iT[63][k]*coeff[ 63*line ] 
        + iT[65][k]*coeff[ 65*line ] 
        + iT[67][k]*coeff[ 67*line ] 
        + iT[69][k]*coeff[ 69*line ] 
        + iT[71][k]*coeff[ 71*line ] 
        + iT[73][k]*coeff[ 73*line ] 
        + iT[75][k]*coeff[ 75*line ] 
        + iT[77][k]*coeff[ 77*line ] 
        + iT[79][k]*coeff[ 79*line ] 
        + iT[81][k]*coeff[ 81*line ] 
        + iT[83][k]*coeff[ 83*line ] 
        + iT[85][k]*coeff[ 85*line ]
        + iT[87][k]*coeff[ 87*line ] 
        + iT[89][k]*coeff[ 89*line ] 
        + iT[91][k]*coeff[ 91*line ] 
        + iT[93][k]*coeff[ 93*line ] 
        + iT[95][k]*coeff[ 95*line ]
        + iT[97][k]*coeff[ 97*line ] 
        + iT[99][k]*coeff[ 99*line ] 
        + iT[101][k]*coeff[ 101*line ] 
        + iT[103][k]*coeff[ 103*line ]
        + iT[105][k]*coeff[ 105*line ]
        + iT[107][k]*coeff[ 107*line ] 
        + iT[109][k]*coeff[ 109*line ] 
        + iT[111][k]*coeff[ 111*line ] 
        + iT[113][k]*coeff[ 113*line ] 
        + iT[115][k]*coeff[ 115*line ]
        + iT[117][k]*coeff[ 117*line ] 
        + iT[119][k]*coeff[ 119*line ] 
        + iT[121][k]*coeff[ 121*line ] 
        + iT[123][k]*coeff[ 123*line ] 
        + iT[125][k]*coeff[ 125*line ]
        + iT[127][k]*coeff[ 127*line ] 
        ;
      }

      for (k=0;k<32;k++) //+4
      {
        EO[k] = iT[ 2][k]*coeff[ 2*line  ] 
        + iT[ 6][k]*coeff[ 6*line  ] 
        + iT[10][k]*coeff[ 10*line ] 
        + iT[14][k]*coeff[ 14*line ] 
        + iT[18][k]*coeff[ 18*line ] 
        + iT[22][k]*coeff[ 22*line ] 
        + iT[26][k]*coeff[ 26*line ] 
        + iT[30][k]*coeff[ 30*line ]
        + iT[34][k]*coeff[ 34*line ] 
        + iT[38][k]*coeff[ 38*line ] 
        + iT[42][k]*coeff[ 42*line ] 
        + iT[46][k]*coeff[ 46*line ] 
        + iT[50][k]*coeff[ 50*line ]
        + iT[54][k]*coeff[ 54*line ] 
        + iT[58][k]*coeff[ 58*line ] 
        + iT[62][k]*coeff[ 62*line ] 
        + iT[66][k]*coeff[ 66*line  ] 
        + iT[70][k]*coeff[ 70*line  ] 
        + iT[74][k]*coeff[ 74*line ] 
        + iT[78][k]*coeff[ 78*line ] 
        + iT[82][k]*coeff[ 82*line ] 
        + iT[86][k]*coeff[ 86*line ] 
        + iT[90][k]*coeff[ 90*line ] 
        + iT[94][k]*coeff[ 94*line ]
        + iT[98][k]*coeff[ 98*line ] 
        + iT[102][k]*coeff[ 102*line ] 
        + iT[106][k]*coeff[ 106*line ] 
        + iT[110][k]*coeff[ 110*line ] 
        + iT[114][k]*coeff[ 114*line ]
        + iT[118][k]*coeff[ 118*line ] 
        + iT[122][k]*coeff[ 122*line ] 
        + iT[126][k]*coeff[ 126*line ] 
        ;
      }
    }

    for (k=0;k<16;k++) //+8
    {
      EEO[k] = iT[4 ][k]*coeff[ 4*line ] 
      + iT[12][k]*coeff[ 12*line ] 
      + iT[20][k]*coeff[ 20*line ] 
      + iT[28][k]*coeff[ 28*line ]
      + iT[36][k]*coeff[ 36*line ]
      + iT[44][k]*coeff[ 44*line ]
      + iT[52][k]*coeff[ 52*line ] 
      + iT[60][k]*coeff[ 60*line ]
      + iT[68][k]*coeff[ 68*line ] 
      + iT[76][k]*coeff[ 76*line ] 
      + iT[84][k]*coeff[ 84*line ] 
      + iT[92][k]*coeff[ 92*line ]
      + iT[100][k]*coeff[ 100*line ]
      + iT[108][k]*coeff[ 108*line ]
      + iT[116][k]*coeff[ 116*line ] 
      + iT[124][k]*coeff[ 124*line ]
      ;
    }

    for (k=0;k<8;k++) //+16
    {
      EEEO[k] = iT[8 ][k]*coeff[ 8*line ] 
      + iT[24][k]*coeff[ 24*line ] 
      + iT[40][k]*coeff[ 40*line ] 
      + iT[56][k]*coeff[ 56*line ] 
      + iT[72][k]*coeff[ 72*line ] 
      + iT[88][k]*coeff[ 88*line ] 
      + iT[104][k]*coeff[ 104*line ] 
      + iT[120][k]*coeff[ 120*line ] 
      ;
    }


    for(k=0; k< 4; k++) //+32
    {
      EEEEO[k] = iT[16][k]*coeff[ 16*line ] 
      + iT[48][k]*coeff[ 48*line ]
      + iT[80][k]*coeff[ 80*line ] 
      + iT[112][k]*coeff[ 112*line ]
      ;
    }

    for(k=0; k< 2; k++) //+64
    {
      EEEEEO[k] = iT[32][k]*coeff[ 32*line ] 
      + iT[96][k]*coeff[ 96*line ]
      ;
    }

    EEEEEE[0] = iT[0][0]*coeff[ 0      ] + iT[64][0]*coeff[ 64*line ];    
    EEEEEE[1] = iT[0][1]*coeff[ 0      ] + iT[64][1]*coeff[ 64*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    for (k=0;k<2;k++)
    {
      EEEEE[k  ] = EEEEEE[k  ] + EEEEEO[k  ];
      EEEEE[k+2] = EEEEEE[1-k] - EEEEEO[1-k];
    }    

    for (k=0;k<4;k++)
    {
      EEEE[k  ] = EEEEE[k  ] + EEEEO[k  ];
      EEEE[k+4] = EEEEE[3-k] - EEEEO[3-k];
    }    

    for (k=0;k<8;k++)
    {
      EEE[k  ] = EEEE[k  ] + EEEO[k  ];
      EEE[k+8] = EEEE[7-k] - EEEO[7-k];
    }    

    for (k=0;k<16;k++)
    {
      EE[k   ] = EEE[k   ] + EEO[k   ];
      EE[k+16] = EEE[15-k] - EEO[15-k];
    }    

    for (k=0;k<32;k++)
    {
      E[k   ] = EE[k   ] + EO[k   ];
      E[k+32] = EE[31-k] - EO[31-k];
    }   

    for (k=0;k<64;k++)
    {
      block[k]    = Clip3( outputMinimum, outputMaximum, (E[k   ] + O[k   ] + add)>>shift );
      block[k+64] = Clip3( outputMinimum, outputMaximum, (E[63-k] - O[63-k] + add)>>shift );
    }
    coeff ++;
    block += 128;
  }
  memset(block, 0, 128*iSkipLine*sizeof(TCoeff));
}
#endif
// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm 
// give identical results
// ********************************** DST-VII **********************************
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST7_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST7_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i;
  Int rnd_factor = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr4[DST7][0] : g_as_DST_MAT_4[0][0];
#else
  const Short *iT = g_as_DST_MAT_4[0][0];
#endif

  Int c[4];
#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoeff = coeff;
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    // Intermediate Variables
    c[0] = block[0] + block[3];
    c[1] = block[1] + block[3];
    c[2] = block[0] - block[1];
    c[3] = iT[2]* block[2];

#if JVET_C0024_QTBT
    coeff[ 0*line ] =  ( iT[0] * c[0] + iT[1] * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 1*line ] =  ( iT[2] * (block[0]+ block[1] - block[3])   + rnd_factor ) >> shift;
    coeff[ 2*line ] =  ( iT[0] * c[2] + iT[1] * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[ 3*line ] =  ( iT[1] * c[2] - iT[0] * c[1]         + c[3]               + rnd_factor ) >> shift;
#else
    coeff[ 0] =  ( iT[0] * c[0] + iT[1] * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 4] =  ( iT[2] * (block[0]+ block[1] - block[3])   + rnd_factor ) >> shift;
    coeff[ 8] =  ( iT[0] * c[2] + iT[1] * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[12] =  ( iT[1] * c[2] - iT[0] * c[1]         + c[3]               + rnd_factor ) >> shift;
#endif

    block+=4;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    coeff = pCoeff + reducedLine;
    for (i=0; i<4; i++)
    {
      memset( coeff, 0, sizeof(TCoeff)*iSkipLine );
      coeff += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST7_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#else
void fastInverseDST7_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#endif
{
  Int i, c[4];
  Int rnd_factor = 1<<(shift-1);

#if COM16_C806_EMT
  const TMatrixCoeff *iT = use ? g_aiTr4[DST7][0] : g_as_DST_MAT_4[0][0];
#else
  const Short *iT = g_as_DST_MAT_4[0][0];
#endif

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {  
    // Intermediate Variables
#if JVET_C0024_QTBT
    c[0] = coeff[0*line] + coeff[2*line];
    c[1] = coeff[2*line] + coeff[3*line];
    c[2] = coeff[0*line] - coeff[3*line];
    c[3] = iT[2]* coeff[1*line];
#else
    c[0] = coeff[0] + coeff[ 8];
    c[1] = coeff[8] + coeff[12];
    c[2] = coeff[0] - coeff[12];
    c[3] = iT[2]* coeff[4];
#endif

    block[0] = Clip3( outputMinimum, outputMaximum, ( iT[0] * c[0] + iT[1] * c[1]     + c[3]               + rnd_factor ) >> shift );
    block[1] = Clip3( outputMinimum, outputMaximum, ( iT[1] * c[2] - iT[0] * c[1]     + c[3]               + rnd_factor ) >> shift );
#if JVET_C0024_QTBT
    block[2] = Clip3( outputMinimum, outputMaximum, ( iT[2] * (coeff[0*line] - coeff[2*line]  + coeff[3*line])      + rnd_factor ) >> shift );
#else
    block[2] = Clip3( outputMinimum, outputMaximum, ( iT[2] * (coeff[0] - coeff[8]  + coeff[12])      + rnd_factor ) >> shift );
#endif
    block[3] = Clip3( outputMinimum, outputMaximum, ( iT[1] * c[0] + iT[0] * c[2]     - c[3]               + rnd_factor ) >> shift );

    block+=4;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<2)*sizeof(TCoeff));
  }
#endif
}

#if COM16_C806_EMT
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST7_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST7_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr8[DST7][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST7_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST7_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT = g_aiTr8[DST7][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<3)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST7_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST7_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr16[DST7][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST7_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST7_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT = g_aiTr16[DST7][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<4)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST7_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST7_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr32[DST7][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#else
  if ( zo )
  {
    TCoeff *tmp = coeff;
    for (i=0; i<(line>>(zo-1)); i++)
    {
      pCoef = coeff;
      iT = g_aiTr32[DST7][0];
      for (j=0; j<uiTrSize/2; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += block[k]*iT[k];
        }
        pCoef[i] = (iSum + rnd_factor)>>shift;
        iT += uiTrSize;
        pCoef += line;
      }
      block+=uiTrSize;
    }

    coeff += (line>>(zo-1));
    if( zo==2 )
    {
      for (j=0; j<uiTrSize/2; j++)
      {
        memset( coeff, 0, sizeof(TCoeff)*uiTrSize/2 );
        coeff += uiTrSize;
      }
    }
    coeff = tmp + uiTrSize*uiTrSize/2;
    memset( coeff, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
  }
  else
  {
    for (i=0; i<line; i++)
    {
      pCoef = coeff;
      iT = g_aiTr32[DST7][0];
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += block[k]*iT[k];
        }
        pCoef[i] = (iSum + rnd_factor)>>shift;
        pCoef += line;
        iT += uiTrSize;
      }
      block += uiTrSize;
    }
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST7_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST7_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT = g_aiTr32[DST7][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<5)*sizeof(TCoeff));
  }
#else
  if ( zo )
  {
    for (i=0; i<(line>>(zo-1)); i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
    /*if( zo==2 )
    {
      memset( block, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
    }*/
  }
  else
  {
    for (i=0; i<line; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
#endif
}





#if JVET_C0024_QTBT
void fastForwardDST7_B64(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                         , Int use
#endif
                         )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr64[DST7][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastInverseDST7_B64(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                         Int use,
#endif
                         const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT = g_aiTr64[DST7][0];
#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<6)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 32;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }

  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}

void fastForwardDST7_B128(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                          , Int use
#endif
                          )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr128[DST7][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastInverseDST7_B128(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                          Int use,
#endif
                          const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT = g_aiTr128[DST7][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<7)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 64;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}
#endif




// ********************************** DCT-VIII **********************************
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT8_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT8_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i;
  Int rnd_factor = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr4[DCT8][0];

  Int c[4];
#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoeff = coeff;
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    // Intermediate Variables
    c[0] = block[0] + block[3];
    c[1] = block[2] + block[0];
    c[2] = block[3] - block[2];
    c[3] = iT[1]* block[1];

#if JVET_C0024_QTBT
    coeff[ 0*line ] =  ( iT[3] * c[0] + iT[2] * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 1*line ] =  ( iT[1] * (block[0] - block[2] - block[3])   + rnd_factor ) >> shift;
    coeff[ 2*line ] =  ( iT[3] * c[2] + iT[2] * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[ 3*line ] =  ( iT[3] * c[1] - iT[2] * c[2]         - c[3]               + rnd_factor ) >> shift;
#else
    coeff[ 0] =  ( iT[3] * c[0] + iT[2] * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 4] =  ( iT[1] * (block[0] - block[2] - block[3])   + rnd_factor ) >> shift;
    coeff[ 8] =  ( iT[3] * c[2] + iT[2] * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[12] =  ( iT[3] * c[1] - iT[2] * c[2]         - c[3]               + rnd_factor ) >> shift;
#endif

    block+=4;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    coeff = pCoeff + reducedLine;
    for (i=0; i<4; i++)
    {
      memset( coeff, 0, sizeof(TCoeff)*iSkipLine );
      coeff += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT8_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#else
void fastInverseDCT8_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#endif
{
  Int i;
  Int rnd_factor = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr4[DCT8][0];
  
  Int c[4];
#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    // Intermediate Variables
#if JVET_C0024_QTBT
    c[0] = coeff[0*line] + coeff[3*line];
    c[1] = coeff[2*line] + coeff[0*line];
    c[2] = coeff[3*line] - coeff[2*line];
    c[3] = iT[1]* coeff[1*line];
#else
    c[0] = coeff[ 0] + coeff[12];
    c[1] = coeff[ 8] + coeff[ 0];
    c[2] = coeff[12] - coeff[ 8];
    c[3] = iT[1]* coeff[4];
#endif

    block[0] =  Clip3(outputMinimum, outputMaximum, ( iT[3] * c[0] + iT[2] * c[1]         + c[3]               + rnd_factor ) >> shift);
#if JVET_C0024_QTBT
    block[1] =  Clip3(outputMinimum, outputMaximum, ( iT[1] * (coeff[0*line] - coeff[2*line] - coeff[3*line])   + rnd_factor ) >> shift);
#else
    block[1] =  Clip3(outputMinimum, outputMaximum, ( iT[1] * (coeff[0] - coeff[8] - coeff[12])   + rnd_factor ) >> shift);
#endif
    block[2] =  Clip3(outputMinimum, outputMaximum, ( iT[3] * c[2] + iT[2] * c[0]         - c[3]               + rnd_factor ) >> shift);
    block[3] =  Clip3(outputMinimum, outputMaximum, ( iT[3] * c[1] - iT[2] * c[2]         - c[3]               + rnd_factor ) >> shift);

    block+=4;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<2)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT8_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT8_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr8[DCT8][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT8_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT8_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT = g_aiTr8[DCT8][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = 8 - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<3)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT8_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT8_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr16[DCT8][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT8_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT8_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT = g_aiTr16[DCT8][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<4)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT8_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT8_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr32[DCT8][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#else
  if ( zo )
  {
    TCoeff *tmp = coeff;
    for (i=0; i<(line>>(zo-1)); i++)
    {
      pCoef = coeff;
      iT = g_aiTr32[DCT8][0];
      for (j=0; j<uiTrSize/2; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += block[k]*iT[k];
        }
        pCoef[i] = (iSum + rnd_factor)>>shift;
        iT += uiTrSize;
        pCoef += line;
      }
      block+=uiTrSize;
    }

    coeff += (line>>(zo-1));
    if( zo==2 )
    {
      for (j=0; j<uiTrSize/2; j++)
      {
        memset( coeff, 0, sizeof(TCoeff)*uiTrSize/2 );
        coeff += uiTrSize;
      }
    }
    coeff = tmp + uiTrSize*uiTrSize/2;
    memset( coeff, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
  }
  else
  {
    for (i=0; i<line; i++)
    {
      pCoef = coeff;
      iT = g_aiTr32[DCT8][0];
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += block[k]*iT[k];
        }
        pCoef[i] = (iSum + rnd_factor)>>shift;
        pCoef += line;
        iT += uiTrSize;
      }
      block += uiTrSize;
    }
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT8_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT8_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT = g_aiTr32[DCT8][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<5)*sizeof(TCoeff));
  }
#else
  if ( zo )
  {
    for (i=0; i<(line>>(zo-1)); i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
    if( zo==2 )
    {
      memset( block, 0, sizeof(TCoeff)*uiTrSize*uiTrSize/2 );
    }
  }
  else
  {
    for (i=0; i<line; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
#endif
}

#if JVET_C0024_QTBT
void fastForwardDCT8_B64(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                         , Int use
#endif
                         )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr64[DCT8][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastInverseDCT8_B64(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                         Int use,
#endif
                         const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT = g_aiTr64[DCT8][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<6)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 32;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}

void fastForwardDCT8_B128(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                          , Int use
#endif
                          )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr128[DCT8][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
void fastInverseDCT8_B128(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                          Int use,
#endif
                          const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT = g_aiTr128[DCT8][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<7)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 64;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }

  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}
#endif




// ********************************** DCT-VIII **********************************
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT5_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT5_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 4;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr4[DCT5][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}


#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT5_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#else
void fastInverseDCT5_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr4[DCT5][0];
  const Int uiTrSize = 4;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<2)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT5_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT5_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr8[DCT5][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT5_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT5_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT = g_aiTr8[DCT5][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<3)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT5_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT5_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr16[DCT5][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT5_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT5_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT = g_aiTr16[DCT5][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<4)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDCT5_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDCT5_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr32[DCT5][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDCT5_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDCT5_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT = g_aiTr32[DCT5][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<5)*sizeof(TCoeff));
  }
#endif
}

#if JVET_C0024_QTBT
void fastForwardDCT5_B64(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                         , Int use
#endif
                         )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr64[DCT5][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

void fastInverseDCT5_B64(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                         Int use,
#endif
                         const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT = g_aiTr64[DCT5][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<6)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 32;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }

  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}

void fastForwardDCT5_B128(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                          , Int use
#endif
                          )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr128[DCT5][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

void fastInverseDCT5_B128(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                          Int use,
#endif
                          const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT = g_aiTr128[DCT5][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<7)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 64;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }

  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}
#endif



// ********************************** DST-I **********************************
#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST1_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST1_B4(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i;
  Int rnd_factor = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr4[DST1][0];

  Int E[2],O[2];
#if JVET_D0077_TRANSFORM_OPT
  TCoeff *pCoeff = coeff;
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {    
    /* E and O */
    E[0] = block[0] + block[3];
    O[0] = block[0] - block[3];
    E[1] = block[1] + block[2];
    O[1] = block[1] - block[2];

    coeff[0     ] = (E[0]*iT[0] + E[1]*iT[1] + rnd_factor)>>shift;
    coeff[line  ] = (O[0]*iT[1] + O[1]*iT[0] + rnd_factor)>>shift;
    coeff[2*line] = (E[0]*iT[1] - E[1]*iT[0] + rnd_factor)>>shift;
    coeff[3*line] = (O[0]*iT[0] - O[1]*iT[1] + rnd_factor)>>shift;

    block += 4;
    coeff ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    coeff = pCoeff + reducedLine;
    for (i=0; i<4; i++)
    {
      memset( coeff, 0, sizeof(TCoeff)*iSkipLine );
      coeff += line;
    }
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST1_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#else
void fastInverseDST1_B4(TCoeff *coeff, TCoeff *block, Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input tmp, output block
#endif
{
  Int i;
  Int rnd_factor = 1<<(shift-1);

  const TMatrixCoeff *iT = g_aiTr4[DST1][0];

  Int E[2],O[2];
#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {    
    /* E and O */
#if JVET_C0024_QTBT
    E[0] = coeff[0*line] + coeff[3*line];
    O[0] = coeff[0*line] - coeff[3*line];
    E[1] = coeff[1*line] + coeff[2*line];
    O[1] = coeff[1*line] - coeff[2*line];
#else
    E[0] = coeff[0*4] + coeff[3*4];
    O[0] = coeff[0*4] - coeff[3*4];
    E[1] = coeff[1*4] + coeff[2*4];
    O[1] = coeff[1*4] - coeff[2*4];
#endif

    block[0] = Clip3(outputMinimum, outputMaximum, (E[0]*iT[0] + E[1]*iT[1] + rnd_factor)>>shift);
    block[1] = Clip3(outputMinimum, outputMaximum, (O[0]*iT[1] + O[1]*iT[0] + rnd_factor)>>shift);
    block[2] = Clip3(outputMinimum, outputMaximum, (E[0]*iT[1] - E[1]*iT[0] + rnd_factor)>>shift);
    block[3] = Clip3(outputMinimum, outputMaximum, (O[0]*iT[0] - O[1]*iT[1] + rnd_factor)>>shift);

    block += 4;
    coeff ++;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<2)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST1_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST1_B8(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr8[DST1][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST1_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST1_B8(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 8;
  const TMatrixCoeff *iT = g_aiTr8[DST1][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line+i]*iT[k*uiTrSize+j];
      }
      block[i*uiTrSize+j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block+(reducedLine<<3), 0, (iSkipLine<<3)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST1_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST1_B16(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr16[DST1][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600 && GCC_VERSION < 40700
__attribute__((optimize("no-tree-vrp")))
#endif
#endif
#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST1_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST1_B16(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 16;
  const TMatrixCoeff *iT = g_aiTr16[DST1][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line+i]*iT[k*uiTrSize+j];
      }
      block[i*uiTrSize+j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block+(reducedLine<<4), 0, (iSkipLine<<4)*sizeof(TCoeff));
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastForwardDST1_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use)  // input block, output coeff
#else
void fastForwardDST1_B32(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int zo, Int use)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    pCoef = coeff;
    iT = g_aiTr32[DST1][0];
#if JVET_D0077_TRANSFORM_OPT
    for (j=0; j<cutoff; j++)
#else
    for (j=0; j<uiTrSize; j++)
#endif
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      pCoef += line;
      iT += uiTrSize;
    }
    block += uiTrSize;
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
#endif
}

#if JVET_D0077_TRANSFORM_OPT
void fastInverseDST1_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#else
void fastInverseDST1_B32(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int zo, Int use, const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
#endif
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 32;
  const TMatrixCoeff *iT = g_aiTr32[DST1][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
#else
  for (i=0; i<line; i++)
#endif
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
#if JVET_D0077_TRANSFORM_OPT
      for (k=0; k<cutoff; k++)
#else
      for (k=0; k<uiTrSize; k++)
#endif
      {
        iSum += coeff[k*line+i]*iT[k*uiTrSize+j];
      }
      block[i*uiTrSize+j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
  }
#if JVET_D0077_TRANSFORM_OPT
  if( iSkipLine )
  {
    memset(block+(reducedLine<<5), 0, (iSkipLine<<5)*sizeof(TCoeff));
  }
#endif
}
#if JVET_C0024_QTBT
void fastForwardDST1_B64(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                         , Int use
#endif
                         )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr64[DST1][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

void fastInverseDST1_B64(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                         Int use,
#endif
                         const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 64;
  const TMatrixCoeff *iT = g_aiTr64[DST1][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<6)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 32;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }

  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}

void fastForwardDST1_B128(TCoeff *block, TCoeff *coeff,Int shift, Int line, Int iSkipLine, Int iSkipLine2
#if JVET_D0077_TRANSFORM_OPT
                          , Int use
#endif
                          )  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT;
  TCoeff *pCoef;

  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    pCoef = coeff;
    iT = g_aiTr128[DST1][0];
    for (j=0; j<cutoff; j++)
    {
      iSum = 0;
      for (k=0; k<uiTrSize; k++)
      {
        iSum += block[k]*iT[k];
      }
      pCoef[i] = (iSum + rnd_factor)>>shift;
      iT += uiTrSize;
      pCoef += line;
    }
    block+=uiTrSize;
  }

  if( iSkipLine )
  {
    pCoef = coeff + reducedLine;
    for (j=0; j<cutoff; j++)
    {
      memset( pCoef, 0, sizeof(TCoeff)*iSkipLine );
      pCoef += line;
    }
  }
  if( iSkipLine2 )
  {
    pCoef = coeff + line*cutoff;
    memset( pCoef, 0, sizeof(TCoeff)*line*iSkipLine2 );
  }
}

void fastInverseDST1_B128(TCoeff *coeff, TCoeff *block,Int shift, Int line, Int iSkipLine, Int iSkipLine2, 
#if JVET_D0077_TRANSFORM_OPT
                          Int use,
#endif
                          const TCoeff outputMinimum, const TCoeff outputMaximum)  // input block, output coeff
{
  Int i, j, k, iSum;
  Int rnd_factor = 1<<(shift-1);

  const Int uiTrSize = 128;
  const TMatrixCoeff *iT = g_aiTr128[DST1][0];

#if JVET_D0077_TRANSFORM_OPT
  const Int  reducedLine = line - iSkipLine;
  const Int  cutoff = uiTrSize - iSkipLine2;
  for (i=0; i<reducedLine; i++)
  {
    for (j=0; j<uiTrSize; j++)
    {
      iSum = 0;
      for (k=0; k<cutoff; k++)
      {
        iSum += coeff[k*line]*iT[k*uiTrSize+j];
      }
      block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
    }
    block+=uiTrSize;
    coeff++;
  }
  if( iSkipLine )
  {
    memset(block, 0, (iSkipLine<<7)*sizeof(TCoeff));
  }
#else
  Bool zo = iSkipLine2 >= 64;

  if ( zo )
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize/2; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  else
  {
    for (i=0; i<line-iSkipLine; i++)
    {
      for (j=0; j<uiTrSize; j++)
      {
        iSum = 0;
        for (k=0; k<uiTrSize; k++)
        {
          iSum += coeff[k*line]*iT[k*uiTrSize+j];
        }
        block[j] = Clip3(outputMinimum, outputMaximum, (Int)(iSum + rnd_factor)>>shift);
      }
      block+=uiTrSize;
      coeff++;
    }
  }
  memset(block, 0, uiTrSize*iSkipLine*sizeof(TCoeff));
#endif
}
#endif

#endif

/** MxN forward transform (2D)
*  \param bitDepth              [in]  bit depth
*  \param block                 [in]  residual block
*  \param coeff                 [out] transform coefficients
*  \param iWidth                [in]  width of transform
*  \param iHeight               [in]  height of transform
*  \param useDST                [in]
*  \param maxLog2TrDynamicRange [in]

*/

#if COM16_C806_EMT
void xTrMxN_EMT(Int bitDepth, TCoeff *block, TCoeff *coeff, Int iWidth, Int iHeight, Bool useDST, const Int maxLog2TrDynamicRange, UChar ucMode, UChar ucTrIdx )
{
  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];

#if JVET_C0024_QTBT
  const Int shift_1st        = ((g_aucConvertToBit[iWidth] + MIN_CU_LOG2) +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
  const Int shift_2nd        = (g_aucConvertToBit[iHeight] + MIN_CU_LOG2) + TRANSFORM_MATRIX_SHIFT + COM16_C806_TRANS_PREC;
  const UInt nLog2WidthMinus1 = g_aucConvertToBit[iWidth] + MIN_CU_LOG2 - 1;  //nLog2WidthMinus1, since transform start from 2-point
  const UInt nLog2HeightMinus1 = g_aucConvertToBit[iHeight] + MIN_CU_LOG2 - 1;  //nLog2HeightMinus1, since transform start from 2-point
  Int iSkipWidth = (iWidth > ZERO_OUT_TH ? iWidth-ZERO_OUT_TH : 0);
  Int iSkipHeight = (iHeight > ZERO_OUT_TH ? iHeight-ZERO_OUT_TH : 0);
#else
  const Int shift_1st        = ((g_aucConvertToBit[iWidth] + 2) +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
  const Int shift_2nd        = (g_aucConvertToBit[iHeight] + 2) + TRANSFORM_MATRIX_SHIFT + COM16_C806_TRANS_PREC;
  const UInt nLog2SizeMinus2 = g_aucConvertToBit[iWidth];

  const Bool bZeroOut        = ( ucMode == INTER_MODE_IDX
#if COM16_C806_T64
    || iWidth==64
#endif
    );
#endif

  TCoeff tmp[ MAX_TU_SIZE * MAX_TU_SIZE ];

  UInt  nTrIdxHor = DCT2, nTrIdxVer = DCT2;
  if( ucMode != INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    UInt  nTrSubsetHor = g_aucTrSetHorz[ucMode];
    UInt  nTrSubsetVer = g_aucTrSetVert[ucMode];
#if FAST_EMT_v1
	if (ucTrIdx == 0)
	{
		nTrIdxHor = DST7;
		nTrIdxVer = DST7;
	}
	else if (ucTrIdx == 1)
	{
		nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor];
		nTrIdxVer = DST7;
	}
	else if (ucTrIdx == 2)
	{
		nTrIdxHor = DST7;
		nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer];
	}
#else
    nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor][ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer][ucTrIdx>>1];
#endif
  }
  if ( ucMode == INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    nTrIdxHor = g_aiTrSubsetInter[ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetInter[ucTrIdx>>1];
  }

#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
  fastFwdTrans[nTrIdxHor][nLog2WidthMinus1]( block, tmp, shift_1st, iHeight, 0, iSkipWidth, 1 );
  fastFwdTrans[nTrIdxVer][nLog2HeightMinus1]( tmp, coeff, shift_2nd,  iWidth, iSkipWidth, iSkipHeight, 1 );
#else
  if( iWidth >= 64 )
  {
    fastFwdTrans[nTrIdxHor][nLog2WidthMinus1]( block, tmp, shift_1st, iHeight, 0, iSkipWidth );
  }
  else
  {
    fastFwdTrans[nTrIdxHor][nLog2WidthMinus1]( block, tmp, shift_1st, iHeight, 0,1 );
  }

  if( iHeight >= 64 )
  {
    fastFwdTrans[nTrIdxVer][nLog2HeightMinus1]( tmp, coeff, shift_2nd,  iWidth, iSkipWidth, iSkipHeight );
  }
  else
  {
    fastFwdTrans[nTrIdxVer][nLog2HeightMinus1]( tmp, coeff, shift_2nd,  iWidth, 0,1 );
  }
#endif
#else
  fastFwdTrans[nTrIdxHor][nLog2SizeMinus2]( block, tmp, shift_1st, iHeight, bZeroOut?1:0, 1 );
  fastFwdTrans[nTrIdxVer][nLog2SizeMinus2]( tmp, coeff, shift_2nd,  iWidth, bZeroOut?2:0, 1 );
#endif
}
#endif

Void xTrMxN(Int bitDepth, TCoeff *block, TCoeff *coeff, Int iWidth, Int iHeight, Bool useDST, const Int maxLog2TrDynamicRange
#if VCEG_AZ08_KLT_COMMON
    , Bool useKLT
#endif
    )
{
#if VCEG_AZ08_KLT_COMMON
    if (useKLT == true)
    {
      xKLTr(bitDepth, block, coeff, iWidth);
      return;
    }
#endif
  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];

#if !COM16_C806_T64
  const
#endif
#if JVET_C0024_QTBT
    Int shift_1st = ((g_aucConvertToBit[iWidth] + MIN_CU_LOG2) +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange;
#else
    Int shift_1st = ((g_aucConvertToBit[iWidth] + 2) +  bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange;
#endif
#if !COM16_C806_T64
  const 
#endif
#if JVET_C0024_QTBT
    Int shift_2nd = (g_aucConvertToBit[iHeight] + MIN_CU_LOG2) + TRANSFORM_MATRIX_SHIFT;
#else
    Int shift_2nd = (g_aucConvertToBit[iHeight] + 2) + TRANSFORM_MATRIX_SHIFT;
#endif

#if COM16_C806_T64
#if JVET_C0024_QTBT
  if( iWidth>=64 || iWidth==2)
  {
      shift_1st += COM16_C806_TRANS_PREC;
  }
  if( iHeight>=64 || iHeight==2)
  {
      shift_2nd += COM16_C806_TRANS_PREC;
  }
#else
  if( iWidth==64 )
  {
    shift_1st += COM16_C806_TRANS_PREC;
    shift_2nd += COM16_C806_TRANS_PREC;
  }
#endif
#endif

  assert(shift_1st >= 0);
  assert(shift_2nd >= 0);

  TCoeff tmp[ MAX_TU_SIZE * MAX_TU_SIZE ];
  Int iSkipWidth = (iWidth > ZERO_OUT_TH ? iWidth-ZERO_OUT_TH : 0);
  Int iSkipHeight = (iHeight > ZERO_OUT_TH ? iHeight-ZERO_OUT_TH : 0);
  switch (iWidth)
  {
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 2:    fastForwardDCT2_B2( block, tmp, shift_1st, iHeight, 0, iSkipWidth, 0 );  break;
#else
    case 2:    fastForwardDCT2_B2( block, tmp, shift_1st, iHeight, 0, 0 );  break;
#endif
#endif
    case 4:
      {
        if ((iHeight == 4) && useDST)    // Check for DCT or DST
        {
           fastForwardDst( block, tmp, shift_1st );
        }
        else
        {
#if JVET_D0077_TRANSFORM_OPT
          partialButterfly4 ( block, tmp, shift_1st, iHeight, iSkipHeight );
#else
          partialButterfly4 ( block, tmp, shift_1st, iHeight );
#endif
        }
      }
      break;

#if JVET_D0077_TRANSFORM_OPT
    case 8:     partialButterfly8 ( block, tmp, shift_1st, iHeight, iSkipHeight );  break;
    case 16:    partialButterfly16( block, tmp, shift_1st, iHeight, iSkipHeight );  break;
    case 32:    partialButterfly32( block, tmp, shift_1st, iHeight, iSkipHeight );  break;
#else
    case 8:     partialButterfly8 ( block, tmp, shift_1st, iHeight );  break;
    case 16:    partialButterfly16( block, tmp, shift_1st, iHeight );  break;
    case 32:    partialButterfly32( block, tmp, shift_1st, iHeight );  break;
#endif
#if COM16_C806_T64
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 64:    fastForwardDCT2_B64( block, tmp, shift_1st, iHeight, 0, iSkipWidth, 0 );  break;
    case 128:   fastForwardDCT2_B128( block, tmp, shift_1st, iHeight, 0, iSkipWidth, 0 );  break;
#else
    case 64:    fastForwardDCT2_B64( block, tmp, shift_1st, iHeight, 0, iSkipWidth );  break;
    case 128:   fastForwardDCT2_B128( block, tmp, shift_1st, iHeight, 0, iSkipWidth );  break;
#endif
#else
    case 64:    fastForwardDCT2_B64( block, tmp, shift_1st, iHeight, 1, 0 );  break;
#endif
#endif
    default:
      assert(0); exit (1); break;
  }

  switch (iHeight)
  {
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 2:    fastForwardDCT2_B2( tmp, coeff, shift_2nd, iWidth, iSkipWidth, iSkipHeight, 0 );  break;
#else
    case 2:    fastForwardDCT2_B2( tmp, coeff, shift_2nd, iWidth, 0, 0 );  break;
#endif
#endif
    case 4:
      {
        if ((iWidth == 4) && useDST)    // Check for DCT or DST
        {
          fastForwardDst( tmp, coeff, shift_2nd );
        }
        else
        {
#if JVET_D0077_TRANSFORM_OPT
          partialButterfly4 ( tmp, coeff, shift_2nd, iWidth, iSkipWidth );
#else
          partialButterfly4 ( tmp, coeff, shift_2nd, iWidth );
#endif
        }
      }
      break;

#if JVET_D0077_TRANSFORM_OPT
    case 8:     partialButterfly8 ( tmp, coeff, shift_2nd, iWidth, iSkipWidth );    break;
    case 16:    partialButterfly16( tmp, coeff, shift_2nd, iWidth, iSkipWidth );    break;
    case 32:    partialButterfly32( tmp, coeff, shift_2nd, iWidth, iSkipWidth );    break;
#else
    case 8:     partialButterfly8 ( tmp, coeff, shift_2nd, iWidth );    break;
    case 16:    partialButterfly16( tmp, coeff, shift_2nd, iWidth );    break;
    case 32:    partialButterfly32( tmp, coeff, shift_2nd, iWidth );    break;
#endif
#if COM16_C806_T64
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 64:    fastForwardDCT2_B64( tmp, coeff, shift_2nd, iWidth, iSkipWidth, iSkipHeight, 0 );  break;
    case 128:   fastForwardDCT2_B128( tmp, coeff, shift_2nd, iWidth, iSkipWidth, iSkipHeight, 0 );  break;
#else
    case 64:    fastForwardDCT2_B64( tmp, coeff, shift_2nd, iWidth, iSkipWidth, iSkipHeight );  break;
    case 128:   fastForwardDCT2_B128( tmp, coeff, shift_2nd, iWidth, iSkipWidth, iSkipHeight );  break;
#endif
#else
    case 64:    fastForwardDCT2_B64( tmp, coeff, shift_2nd, iWidth, 2, 0 );  break;
#endif
#endif
    default:
      assert(0); exit (1); break;
  }
}


/** MxN inverse transform (2D)
*  \param bitDepth              [in]  bit depth
*  \param coeff                 [in]  transform coefficients
*  \param block                 [out] residual block
*  \param iWidth                [in]  width of transform
*  \param iHeight               [in]  height of transform
*  \param useDST                [in]
*  \param maxLog2TrDynamicRange [in]
*/

#if COM16_C806_EMT
#if JVET_C0024_ITSKIP
void xITrMxN_EMT(Int bitDepth, TCoeff *coeff, TCoeff *block, Int iWidth, Int iHeight, UInt uiSkipWidth, UInt uiSkipHeight, Bool useDST, const Int maxLog2TrDynamicRange, UChar ucMode, UChar ucTrIdx )
#else
void xITrMxN_EMT(Int bitDepth, TCoeff *coeff, TCoeff *block, Int iWidth, Int iHeight, Bool useDST, const Int maxLog2TrDynamicRange, UChar ucMode, UChar ucTrIdx )
#endif
{
  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_INVERSE];

  const Int shift_1st        = TRANSFORM_MATRIX_SHIFT + 1 + COM16_C806_TRANS_PREC; //1 has been added to shift_1st at the expense of shift_2nd
  const Int shift_2nd        = (TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1) - bitDepth + COM16_C806_TRANS_PREC;
  const TCoeff clipMinimum   = -(1 << maxLog2TrDynamicRange);
  const TCoeff clipMaximum   =  (1 << maxLog2TrDynamicRange) - 1;
#if JVET_C0024_QTBT
  const UInt nLog2WidthMinus1 = g_aucConvertToBit[iWidth] + MIN_CU_LOG2 - 1;  //nLog2WidthMinus1, since transform start from 2-point
  const UInt nLog2HeightMinus1 = g_aucConvertToBit[iHeight] + MIN_CU_LOG2 - 1;  //nLog2HeightMinus1, since transform start from 2-point
#else
  const UInt nLog2SizeMinus2 = g_aucConvertToBit[iWidth];
  const Bool bZeroOut        = ( ucMode == INTER_MODE_IDX
#if COM16_C806_T64
    || iWidth==64
#endif
    );
#endif

  TCoeff tmp[ MAX_TU_SIZE * MAX_TU_SIZE ];

  UInt  nTrIdxHor = DCT2, nTrIdxVer = DCT2;
  if ( ucMode != INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    UInt  nTrSubsetHor = g_aucTrSetHorz[ucMode];
    UInt  nTrSubsetVer = g_aucTrSetVert[ucMode];
#if FAST_EMT_v1
	if (ucTrIdx == 0)
	{
		nTrIdxHor = DST7;
		nTrIdxVer = DST7;
	}
	else if (ucTrIdx == 1)
	{
		nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor];
		nTrIdxVer = DST7;
	}
	else if (ucTrIdx == 2)
	{
		nTrIdxHor = DST7;
		nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer];
	}
#else
    nTrIdxHor = g_aiTrSubsetIntra[nTrSubsetHor][ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetIntra[nTrSubsetVer][ucTrIdx>>1];
#endif
  }
  if ( ucMode == INTER_MODE_IDX && ucTrIdx != DCT2_EMT )
  {
    nTrIdxHor = g_aiTrSubsetInter[ucTrIdx &1];
    nTrIdxVer = g_aiTrSubsetInter[ucTrIdx>>1];
  }

#if JVET_C0024_QTBT
#if JVET_C0024_ITSKIP
  // uiSkipWidth and uiSkipHeight (initilized as 0) may not be set at encoder side, need to be calculated  
  if( iWidth > ZERO_OUT_TH && uiSkipWidth == 0 )
  {
    uiSkipWidth = iWidth - ZERO_OUT_TH;
  }
  if( iHeight > ZERO_OUT_TH && uiSkipHeight == 0 )
  {
    uiSkipHeight = iHeight - ZERO_OUT_TH;
  }
#else
  UInt uiSkipWidth = (iWidth > ZERO_OUT_TH ? iWidth - ZERO_OUT_TH : 0);
  UInt uiSkipHeight = (iHeight > ZERO_OUT_TH ? iHeight - ZERO_OUT_TH : 0);
#endif

#if JVET_D0077_TRANSFORM_OPT
  fastInvTrans[nTrIdxVer][nLog2HeightMinus1]( coeff, tmp, shift_1st,  iWidth, uiSkipWidth, uiSkipHeight, 1, clipMinimum, clipMaximum );
  fastInvTrans[nTrIdxHor][nLog2WidthMinus1]( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, 1, clipMinimum, clipMaximum );
#else
  if (nLog2HeightMinus1 + 1 >= 6)
  {
#if JVET_C0024_ITSKIP && (ZERO_OUT_FIX || JVET_D0077_TRANSFORM_OPT)
    assert( uiSkipHeight >= iHeight - ZERO_OUT_TH );
#endif
    fastInvTrans[nTrIdxVer][nLog2HeightMinus1]( coeff, tmp, shift_1st,  iWidth, uiSkipWidth, uiSkipHeight, clipMinimum, clipMaximum );
  }
  else
  {
    fastInvTrans[nTrIdxVer][nLog2HeightMinus1]( coeff, tmp, shift_1st,  iWidth, 0, 1, clipMinimum, clipMaximum );
  }

  if (nLog2WidthMinus1 + 1>=6)
  {
#if JVET_C0024_ITSKIP && JVET_D0077_TRANSFORM_OPT
    assert( uiSkipWidth >= iWidth - ZERO_OUT_TH );
#endif
    fastInvTrans[nTrIdxHor][nLog2WidthMinus1]( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, clipMinimum, clipMaximum );
  }
  else
  {
    fastInvTrans[nTrIdxHor][nLog2WidthMinus1]( tmp, block, shift_2nd, iHeight, 0, 1, clipMinimum, clipMaximum );
  }
#endif
#else
  fastInvTrans[nTrIdxVer][nLog2SizeMinus2]( coeff, tmp, shift_1st,  iWidth, bZeroOut?2:0, 1, clipMinimum, clipMaximum );
  fastInvTrans[nTrIdxHor][nLog2SizeMinus2]( tmp, block, shift_2nd, iHeight, bZeroOut?1:0, 1, clipMinimum, clipMaximum );
#endif
}
#endif

#if JVET_C0024_ITSKIP
Void xITrMxN(Int bitDepth, TCoeff *coeff, TCoeff *block, Int iWidth, Int iHeight, UInt uiSkipWidth, UInt uiSkipHeight, Bool useDST, const Int maxLog2TrDynamicRange
#if VCEG_AZ08_KLT_COMMON
    , Bool useKLT
#endif
)
#else
Void xITrMxN(Int bitDepth, TCoeff *coeff, TCoeff *block, Int iWidth, Int iHeight, Bool useDST, const Int maxLog2TrDynamicRange
#if VCEG_AZ08_KLT_COMMON
    , Bool useKLT
#endif
    )
#endif
{
#if VCEG_AZ08_KLT_COMMON
    if (useKLT == true)
    {
      assert(iWidth == iHeight);
      xIKLTr(bitDepth, coeff, block, iHeight);
      return;
    }
#endif
  const Int TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_INVERSE];

  Int shift_1st = TRANSFORM_MATRIX_SHIFT + 1; //1 has been added to shift_1st at the expense of shift_2nd
  Int shift_2nd = (TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1) - bitDepth;
  const TCoeff clipMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff clipMaximum =  (1 << maxLog2TrDynamicRange) - 1;

#if COM16_C806_T64
#if JVET_C0024_QTBT
  if( iWidth>=64 || iWidth==2)
  {
      shift_2nd += COM16_C806_TRANS_PREC;
  }
  if( iHeight>=64 || iHeight==2)
  {
      shift_1st += COM16_C806_TRANS_PREC;
  }
#else
  if( iWidth==64 )
  {
    shift_1st += COM16_C806_TRANS_PREC;
    shift_2nd += COM16_C806_TRANS_PREC;
  }
#endif
#endif

  assert(shift_1st >= 0);
  assert(shift_2nd >= 0);

  TCoeff tmp[MAX_TU_SIZE * MAX_TU_SIZE];
#if JVET_D0077_TRANSFORM_OPT && !JVET_C0024_ITSKIP
  UInt uiSkipWidth = (iWidth > ZERO_OUT_TH ? iWidth - ZERO_OUT_TH : 0);
  UInt uiSkipHeight = (iHeight > ZERO_OUT_TH ? iHeight - ZERO_OUT_TH : 0);
#endif

  switch (iHeight)
  {
#if JVET_C0024_QTBT
    case 2:
#if JVET_D0077_TRANSFORM_OPT
        fastInverseDCT2_B2(coeff, tmp, shift_1st, iWidth, uiSkipWidth, uiSkipHeight, 0, clipMinimum, clipMaximum); break;
#else
        fastInverseDCT2_B2(coeff, tmp, shift_1st, iWidth, 0, 0, clipMinimum, clipMaximum); break;
#endif
#endif
    case 4:
      {
        if ((iWidth == 4) && useDST)    // Check for DCT or DST
        {
          fastInverseDst( coeff, tmp, shift_1st, clipMinimum, clipMaximum);
        }
        else
        {
#if JVET_D0077_TRANSFORM_OPT
          partialButterflyInverse4 ( coeff, tmp, shift_1st, iWidth, uiSkipWidth, clipMinimum, clipMaximum);
#else
          partialButterflyInverse4 ( coeff, tmp, shift_1st, iWidth, clipMinimum, clipMaximum);
#endif
        }
      }
      break;

#if JVET_D0077_TRANSFORM_OPT
    case  8: partialButterflyInverse8 ( coeff, tmp, shift_1st, iWidth, uiSkipWidth, clipMinimum, clipMaximum); break;
    case 16: partialButterflyInverse16( coeff, tmp, shift_1st, iWidth, uiSkipWidth, clipMinimum, clipMaximum); break;
    case 32: partialButterflyInverse32( coeff, tmp, shift_1st, iWidth, uiSkipWidth, clipMinimum, clipMaximum); break;
#else
    case  8: partialButterflyInverse8 ( coeff, tmp, shift_1st, iWidth, clipMinimum, clipMaximum); break;
    case 16: partialButterflyInverse16( coeff, tmp, shift_1st, iWidth, clipMinimum, clipMaximum); break;
    case 32: partialButterflyInverse32( coeff, tmp, shift_1st, iWidth, clipMinimum, clipMaximum); break;
#endif
#if COM16_C806_T64
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 64: fastInverseDCT2_B64( coeff, tmp, shift_1st, iWidth, uiSkipWidth, uiSkipHeight, 0, clipMinimum, clipMaximum); break;
    case 128: fastInverseDCT2_B128( coeff, tmp, shift_1st, iWidth, uiSkipWidth, uiSkipHeight, 0, clipMinimum, clipMaximum); break;
#else
    case 64: fastInverseDCT2_B64( coeff, tmp, shift_1st, iWidth, uiSkipWidth, uiSkipHeight, clipMinimum, clipMaximum); break;
    case 128: fastInverseDCT2_B128( coeff, tmp, shift_1st, iWidth, uiSkipWidth, uiSkipHeight, clipMinimum, clipMaximum); break;
#endif
#else
    case 64: fastInverseDCT2_B64( coeff, tmp, shift_1st, iWidth, 2, 0, clipMinimum, clipMaximum); break;
#endif
#endif
    default:
      assert(0); exit (1); break;
  }

  switch (iWidth)
  {
    // Clipping here is not in the standard, but is used to protect the "Pel" data type into which the inverse-transformed samples will be copied
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 2: fastInverseDCT2_B2( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#else
    case 2: fastInverseDCT2_B2( tmp, block, shift_2nd, iHeight, 0, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#endif
#endif
    case 4:
      {
        if ((iHeight == 4) && useDST)    // Check for DCT or DST
        {
          fastInverseDst( tmp, block, shift_2nd, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max() );
        }
        else
        {
#if JVET_D0077_TRANSFORM_OPT
          partialButterflyInverse4 ( tmp, block, shift_2nd, iHeight, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max());
#else
          partialButterflyInverse4 ( tmp, block, shift_2nd, iHeight, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max());
#endif
        }
      }
      break;

#if JVET_D0077_TRANSFORM_OPT
    case  8: partialButterflyInverse8 ( tmp, block, shift_2nd, iHeight, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 16: partialButterflyInverse16( tmp, block, shift_2nd, iHeight, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 32: partialButterflyInverse32( tmp, block, shift_2nd, iHeight, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#else
    case  8: partialButterflyInverse8 ( tmp, block, shift_2nd, iHeight, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 16: partialButterflyInverse16( tmp, block, shift_2nd, iHeight, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 32: partialButterflyInverse32( tmp, block, shift_2nd, iHeight, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#endif
#if COM16_C806_T64
#if JVET_C0024_QTBT
#if JVET_D0077_TRANSFORM_OPT
    case 64: fastInverseDCT2_B64( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 128: fastInverseDCT2_B128( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#else
    case 64: fastInverseDCT2_B64( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
    case 128: fastInverseDCT2_B128( tmp, block, shift_2nd, iHeight, 0, uiSkipWidth, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#endif
#else
    case 64: fastInverseDCT2_B64( tmp, block, shift_2nd, iHeight, 1, 0, std::numeric_limits<Pel>::min(), std::numeric_limits<Pel>::max()); break;
#endif
#endif
    default:
      assert(0); exit (1); break;
  }
}


// To minimize the distortion only. No rate is considered.
Void TComTrQuant::signBitHidingHDQ( TCoeff* pQCoef, TCoeff* pCoef, TCoeff* deltaU, const TUEntropyCodingParameters &codingParameters, const Int maxLog2TrDynamicRange )
{
  const UInt width     = codingParameters.widthInGroups  << MLS_CG_LOG2_WIDTH;
  const UInt height    = codingParameters.heightInGroups << MLS_CG_LOG2_HEIGHT;
  const UInt groupSize = 1 << MLS_CG_SIZE;

  const TCoeff entropyCodingMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff entropyCodingMaximum =  (1 << maxLog2TrDynamicRange) - 1;

  Int lastCG = -1;
  Int absSum = 0 ;
  Int n ;

  for( Int subSet = (width*height-1) >> MLS_CG_SIZE; subSet >= 0; subSet-- )
  {
    Int  subPos = subSet << MLS_CG_SIZE;
    Int  firstNZPosInCG=groupSize , lastNZPosInCG=-1 ;
    absSum = 0 ;

    for(n = groupSize-1; n >= 0; --n )
    {
      if( pQCoef[ codingParameters.scan[ n + subPos ]] )
      {
        lastNZPosInCG = n;
        break;
      }
    }

    for(n = 0; n <groupSize; n++ )
    {
      if( pQCoef[ codingParameters.scan[ n + subPos ]] )
      {
        firstNZPosInCG = n;
        break;
      }
    }

    for(n = firstNZPosInCG; n <=lastNZPosInCG; n++ )
    {
      absSum += Int(pQCoef[ codingParameters.scan[ n + subPos ]]);
    }

    if(lastNZPosInCG>=0 && lastCG==-1)
    {
      lastCG = 1 ;
    }

    if( lastNZPosInCG-firstNZPosInCG>=SBH_THRESHOLD )
    {
      UInt signbit = (pQCoef[codingParameters.scan[subPos+firstNZPosInCG]]>0?0:1) ;
      if( signbit!=(absSum&0x1) )  //compare signbit with sum_parity
      {
        TCoeff curCost    = std::numeric_limits<TCoeff>::max();
        TCoeff minCostInc = std::numeric_limits<TCoeff>::max();
        Int minPos =-1, finalChange=0, curChange=0;

        for( n = (lastCG==1?lastNZPosInCG:groupSize-1) ; n >= 0; --n )
        {
          UInt blkPos   = codingParameters.scan[ n+subPos ];
          if(pQCoef[ blkPos ] != 0 )
          {
            if(deltaU[blkPos]>0)
            {
              curCost = - deltaU[blkPos];
              curChange=1 ;
            }
            else
            {
              //curChange =-1;
              if(n==firstNZPosInCG && abs(pQCoef[blkPos])==1)
              {
                curCost = std::numeric_limits<TCoeff>::max();
              }
              else
              {
                curCost = deltaU[blkPos];
                curChange =-1;
              }
            }
          }
          else
          {
            if(n<firstNZPosInCG)
            {
              UInt thisSignBit = (pCoef[blkPos]>=0?0:1);
              if(thisSignBit != signbit )
              {
                curCost = std::numeric_limits<TCoeff>::max();
              }
              else
              {
                curCost = - (deltaU[blkPos])  ;
                curChange = 1 ;
              }
            }
            else
            {
              curCost = - (deltaU[blkPos])  ;
              curChange = 1 ;
            }
          }

          if( curCost<minCostInc)
          {
            minCostInc = curCost ;
            finalChange = curChange ;
            minPos = blkPos ;
          }
        } //CG loop

        if(pQCoef[minPos] == entropyCodingMaximum || pQCoef[minPos] == entropyCodingMinimum)
        {
          finalChange = -1;
        }

        if(pCoef[minPos]>=0)
        {
          pQCoef[minPos] += finalChange ;
        }
        else
        {
          pQCoef[minPos] -= finalChange ;
        }
      } // Hide
    }
    if(lastCG==1)
    {
      lastCG=0 ;
    }
  } // TU loop

  return;
}


Void TComTrQuant::xQuant(       TComTU       &rTu,
                                TCoeff      * pSrc,
                                TCoeff      * pDes,
#if ADAPTIVE_QP_SELECTION
                                TCoeff      *pArlDes,
#endif
                                TCoeff       &uiAbsSum,
                          const ComponentID   compID,
                          const QpParam      &cQP 
                          )
{
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt uiWidth        = rect.width;
  const UInt uiHeight       = rect.height;
  TComDataCU* pcCU          = rTu.getCU();
  const UInt uiAbsPartIdx   = rTu.GetAbsPartIdxTU();
  const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));

  TCoeff* piCoef    = pSrc;
  TCoeff* piQCoef   = pDes;
#if ADAPTIVE_QP_SELECTION
  TCoeff* piArlCCoef = pArlDes;
#endif

  const Bool useTransformSkip      = pcCU->getTransformSkip(uiAbsPartIdx, compID);
  const Int  maxLog2TrDynamicRange = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));

  Bool useRDOQ = useTransformSkip ? m_useRDOQTS : m_useRDOQ;
#if JVET_C0024_QTBT
  if (uiWidth<4 || uiHeight<4)
  {
    useRDOQ = false;
  }
#endif
  if ( useRDOQ && (isLuma(compID) || RDOQ_CHROMA) )
  {
#if T0196_SELECTIVE_RDOQ
    if ( !m_useSelectiveRDOQ || xNeedRDOQ( rTu, piCoef, compID, cQP ) )
    {
#endif
#if ADAPTIVE_QP_SELECTION
      xRateDistOptQuant( rTu, piCoef, pDes, pArlDes, uiAbsSum, compID, cQP );
#else
      xRateDistOptQuant( rTu, piCoef, pDes, uiAbsSum, compID, cQP );
#endif
#if T0196_SELECTIVE_RDOQ
    }
    else
    {
      memset( pDes, 0, sizeof( TCoeff ) * uiWidth *uiHeight );
      uiAbsSum = 0;
    }
#endif
  }
  else
  {
    TUEntropyCodingParameters codingParameters;
    getTUEntropyCodingParameters(codingParameters, rTu, compID);

    const TCoeff entropyCodingMinimum = -(1 << maxLog2TrDynamicRange);
    const TCoeff entropyCodingMaximum =  (1 << maxLog2TrDynamicRange) - 1;

    TCoeff deltaU[MAX_TU_SIZE * MAX_TU_SIZE];

    const UInt uiLog2TrSize = rTu.GetEquivalentLog2TrSize(compID);

    Int scalingListType = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
    assert(scalingListType < SCALING_LIST_NUM);
#if JVET_C0024_QTBT
    const UInt uiLog2TrWidth = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
    const UInt uiLog2TrHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
    Int *piQuantCoeff = getQuantCoeff(scalingListType, cQP.rem, uiLog2TrWidth-1, uiLog2TrHeight-1);
#else
    Int *piQuantCoeff = getQuantCoeff(scalingListType, cQP.rem, uiLog2TrSize-2);
#endif

    const Bool enableScalingLists             = getUseScalingList(uiWidth, uiHeight, (pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0));
    const Int  defaultQuantisationCoefficient = g_quantScales[cQP.rem];

    /* for 422 chroma blocks, the effective scaling applied during transformation is not a power of 2, hence it cannot be
     * implemented as a bit-shift (the quantised result will be sqrt(2) * larger than required). Alternatively, adjust the
     * uiLog2TrSize applied in iTransformShift, such that the result is 1/sqrt(2) the required result (i.e. smaller)
     * Then a QP+3 (sqrt(2)) or QP-3 (1/sqrt(2)) method could be used to get the required result
     */

    // Represents scaling through forward transform
    Int iTransformShift = getTransformShift(channelBitDepth, uiLog2TrSize, maxLog2TrDynamicRange);
    if (useTransformSkip && pcCU->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag())
    {
      iTransformShift = std::max<Int>(0, iTransformShift);
    }

#if JVET_C0024_QTBT
  Int iWHScale = 1;
  if ((g_aucConvertToBit[ uiWidth ] + g_aucConvertToBit[ uiHeight ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iTransformShift += 7;
    iWHScale = 181;
  }
#endif

    const Int iQBits = QUANT_SHIFT + cQP.per + iTransformShift;
    // QBits will be OK for any internal bit depth as the reduction in transform shift is balanced by an increase in Qp_per due to QpBDOffset

#if ADAPTIVE_QP_SELECTION
    Int iQBitsC = MAX_INT;
    Int iAddC   = MAX_INT;

    if (m_bUseAdaptQpSelect)
    {
      iQBitsC = iQBits - ARL_C_PRECISION;
      iAddC   = 1 << (iQBitsC-1);
    }
#endif

#if JVET_C0024_ENCODER_OVERFLOW_FIX
    const Int64 iAdd   = (Int64)(pcCU->getSlice()->getSliceType()==I_SLICE ? 171 : 85) << ((Int64)iQBits-9);
#else
    const Int iAdd   = (pcCU->getSlice()->getSliceType()==I_SLICE ? 171 : 85) << (iQBits-9);
#endif
    const Int qBits8 = iQBits - 8;

    for( Int uiBlockPos = 0; uiBlockPos < uiWidth*uiHeight; uiBlockPos++ )
    {
      const TCoeff iLevel   = piCoef[uiBlockPos];
      const TCoeff iSign    = (iLevel < 0 ? -1: 1);

      const Int64  tmpLevel = (Int64)abs(iLevel) * (enableScalingLists ? piQuantCoeff[uiBlockPos] : defaultQuantisationCoefficient);

#if ADAPTIVE_QP_SELECTION
      if( m_bUseAdaptQpSelect )
      {
#if JVET_C0024_QTBT
      piArlCCoef[uiBlockPos] = (TCoeff)((tmpLevel * iWHScale + iAddC ) >> iQBitsC);
#else
        piArlCCoef[uiBlockPos] = (TCoeff)((tmpLevel + iAddC ) >> iQBitsC);
#endif
      }
#endif

#if JVET_C0024_QTBT
      const TCoeff quantisedMagnitude = TCoeff((tmpLevel * iWHScale + iAdd ) >> iQBits);
#if JVET_C0024_ENCODER_OVERFLOW_FIX
    deltaU[uiBlockPos] = (TCoeff)((tmpLevel * iWHScale - ((Int64)quantisedMagnitude<<iQBits) )>> qBits8);
#else
    deltaU[uiBlockPos] = (TCoeff)((tmpLevel * iWHScale - (quantisedMagnitude<<iQBits) )>> qBits8);
#endif
#else
#if PVC_AQS
      TCoeff quantisedMagnitude;
      if (compID == COMPONENT_Y)
      {
        Int iEstQSFactorShift14 = m_iEstQSFactor << (14 - ESM_SHIFT);
        Int64 tmpLevelShift = tmpLevel * iEstQSFactorShift14;
        quantisedMagnitude = TCoeff((tmpLevelShift + ((Int64) iAdd << 14)) >> (iQBits + 14));
        deltaU[uiBlockPos] = (TCoeff) ((tmpLevelShift - (quantisedMagnitude << (iQBits + 14))) >> (qBits8 + 14));
      }
      else
      {
        quantisedMagnitude = TCoeff((tmpLevel + iAdd) >> iQBits);
        deltaU[uiBlockPos] = (TCoeff) ((tmpLevel - (quantisedMagnitude << iQBits)) >> qBits8);
      }
#else
      const TCoeff quantisedMagnitude = TCoeff((tmpLevel + iAdd) >> iQBits);
      deltaU[uiBlockPos] = (TCoeff) ((tmpLevel - (quantisedMagnitude << iQBits)) >> qBits8);
#endif
#endif

      uiAbsSum += quantisedMagnitude;
      const TCoeff quantisedCoefficient = quantisedMagnitude * iSign;

      piQCoef[uiBlockPos] = Clip3<TCoeff>( entropyCodingMinimum, entropyCodingMaximum, quantisedCoefficient );
    } // for n

#if JVET_C0024_QTBT
    if( pcCU->getSlice()->getPPS()->getSignHideFlag() && uiWidth>=4 && uiHeight>=4)
#else
    if( pcCU->getSlice()->getPPS()->getSignHideFlag() )
#endif
    {
      if(uiAbsSum >= 2) //this prevents TUs with only one coefficient of value 1 from being tested
      {
        signBitHidingHDQ( piQCoef, piCoef, deltaU, codingParameters, maxLog2TrDynamicRange ) ;
      }
    }
  } //if RDOQ
  //return;
}

#if T0196_SELECTIVE_RDOQ
Bool TComTrQuant::xNeedRDOQ( TComTU &rTu, TCoeff * pSrc, const ComponentID compID, const QpParam &cQP )
{
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt uiWidth        = rect.width;
  const UInt uiHeight       = rect.height;
  TComDataCU* pcCU          = rTu.getCU();
  const UInt uiAbsPartIdx   = rTu.GetAbsPartIdxTU();
  const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));

  TCoeff* piCoef    = pSrc;

  const Bool useTransformSkip      = pcCU->getTransformSkip(uiAbsPartIdx, compID);
  const Int  maxLog2TrDynamicRange = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));

  const UInt uiLog2TrSize = rTu.GetEquivalentLog2TrSize(compID);

  Int scalingListType = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
  assert(scalingListType < SCALING_LIST_NUM);
#if JVET_C0024_QTBT
  const UInt uiLog2TrWidth = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
  const UInt uiLog2TrHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
  Int *piQuantCoeff = getQuantCoeff(scalingListType, cQP.rem, uiLog2TrWidth-1, uiLog2TrHeight-1);
#else
  Int *piQuantCoeff = getQuantCoeff(scalingListType, cQP.rem, uiLog2TrSize-2);
#endif

  const Bool enableScalingLists             = getUseScalingList(uiWidth, uiHeight, (pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0));
  const Int  defaultQuantisationCoefficient = g_quantScales[cQP.rem];

  /* for 422 chroma blocks, the effective scaling applied during transformation is not a power of 2, hence it cannot be
    * implemented as a bit-shift (the quantised result will be sqrt(2) * larger than required). Alternatively, adjust the
    * uiLog2TrSize applied in iTransformShift, such that the result is 1/sqrt(2) the required result (i.e. smaller)
    * Then a QP+3 (sqrt(2)) or QP-3 (1/sqrt(2)) method could be used to get the required result
    */

  // Represents scaling through forward transform
  Int iTransformShift = getTransformShift(channelBitDepth, uiLog2TrSize, maxLog2TrDynamicRange);
  if (useTransformSkip && pcCU->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag())
  {
    iTransformShift = std::max<Int>(0, iTransformShift);
  }

#if JVET_C0024_QTBT
  Int iWHScale = 1;
  if ((g_aucConvertToBit[ uiWidth ] + g_aucConvertToBit[ uiHeight ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iTransformShift += 7;
    iWHScale = 181;
  }
#endif
  const Int iQBits = QUANT_SHIFT + cQP.per + iTransformShift;
  // QBits will be OK for any internal bit depth as the reduction in transform shift is balanced by an increase in Qp_per due to QpBDOffset

  // iAdd is different from the iAdd used in normal quantization
#if JVET_C0024_QTBT
#if JVET_C0024_ENCODER_OVERFLOW_FIX
  const Int64 iAdd   = (Int64)(compID == COMPONENT_Y ? 171 : 256) << (iQBits-9);
#else
  const Int64 iAdd   = (compID == COMPONENT_Y ? 171 : 256) << (iQBits-9);
#endif
#else
  const Int iAdd   = (compID == COMPONENT_Y ? 171 : 256) << (iQBits-9);
#endif

  for( Int uiBlockPos = 0; uiBlockPos < uiWidth*uiHeight; uiBlockPos++ )
  {
    const TCoeff iLevel   = piCoef[uiBlockPos];
    const Int64  tmpLevel = (Int64)abs(iLevel) * (enableScalingLists ? piQuantCoeff[uiBlockPos] : defaultQuantisationCoefficient);
#if JVET_C0024_QTBT
    const TCoeff quantisedMagnitude = TCoeff((tmpLevel * iWHScale + iAdd ) >> iQBits);
#else
    const TCoeff quantisedMagnitude = TCoeff((tmpLevel + iAdd ) >> iQBits);
#endif

    if ( quantisedMagnitude != 0 )
    {
      return true;
    }
  } // for n
  return false;
}
#endif

#if ITH_C3_DEQUANT==C_FXN
Void rightShiftWithScaleFactorAndDequantCoef(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int *piDequantCoef, Int rightShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	const Intermediate_Int iAdd = 1 << (rightShift - 1);

#if ITH_C2_ZR_ZC_SL
	Int n = 0;
	for (UInt j = 0; j < uiHeight - uiTUSkipHeight; j++)
	{
		n = j*uiWidth;
		for (UInt i = 0; i < uiWidth - uiTUSkipWidth; i++)
		{
#else
	Int numSamplesInBlock = uiWidth * uiHeight;
	for (Int n = 0; n < numSamplesInBlock; n++)
	{
#endif
		const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
#if JVET_C0024_QTBT
		const Intermediate_Int iCoeffQ = ((Intermediate_Int(clipQCoef) * piDequantCoef[n] * scale) + iAdd) >> rightShift;
#else
		const Intermediate_Int iCoeffQ = ((Intermediate_Int(clipQCoef) * piDequantCoef[n]) + iAdd) >> rightShift;
#endif

		piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ));
#if ITH_C2_ZR_ZC_SL
		n++;
	}
	memset(piCoef + n, 0, uiTUSkipWidth*sizeof(TCoeff));
		}
memset(piCoef + (uiHeight - uiTUSkipHeight)*uiWidth, 0, uiTUSkipHeight*uiWidth * sizeof(TCoeff));
#else
	}
#endif
}

Void leftShiftWithScaleFactorAndDequantCoef(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int *piDequantCoef, Int leftShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
#if ITH_C2_ZR_ZC_SL
	Int n = 0;
	for (UInt j = 0; j < uiHeight - uiTUSkipHeight; j++)
	{
		n = j*uiWidth;
		for (UInt i = 0; i < uiWidth - uiTUSkipWidth; i++)
		{
#else
	Int numSamplesInBlock = uiWidth * uiHeight;
	for (Int n = 0; n < numSamplesInBlock; n++)
	{
#endif
		const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
#if JVET_C0024_QTBT
		const Intermediate_Int iCoeffQ = (Intermediate_Int(clipQCoef) * piDequantCoef[n] * scale) << leftShift;
#else
		const Intermediate_Int iCoeffQ = (Intermediate_Int(clipQCoef) * piDequantCoef[n]) << leftShift;
#endif

		piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ));
#if ITH_C2_ZR_ZC_SL
		n++;
	}
	memset(piCoef + n, 0, uiTUSkipWidth*sizeof(TCoeff));
		}
memset(piCoef + (uiHeight - uiTUSkipHeight)*uiWidth, 0, uiTUSkipHeight*uiWidth * sizeof(TCoeff));
#else
	}
#endif
}

Void rightShiftWithScaleFactor(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int rightShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	const Intermediate_Int iAdd = 1 << (rightShift - 1);

#if JVET_C0024_ITSKIP
	Int n = 0;
	for (UInt j = 0; j < uiHeight - uiTUSkipHeight; j++)
	{
		n = j*uiWidth;
		for (UInt i = 0; i < uiWidth - uiTUSkipWidth; i++)
		{
#else
	for (Int n = 0; n < numSamplesInBlock; n++)
	{
#endif
		const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
#if JVET_C0024_QTBT
		const Intermediate_Int iCoeffQ = ((Intermediate_Int(clipQCoef) * scale) + iAdd) >> rightShift;
#else
		const Intermediate_Int iCoeffQ = ((Intermediate_Int(clipQCoef)) + iAdd) >> rightShift;
#endif

		piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ));
#if JVET_C0024_ITSKIP
		n++;
	}
	memset(piCoef + n, 0, uiTUSkipWidth * sizeof(TCoeff));
		}
memset(piCoef + (uiHeight - uiTUSkipHeight)*uiWidth, 0, uiTUSkipHeight * uiWidth * sizeof(TCoeff));
#else
	}
#endif
}

Void leftShiftWithScaleFactor(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int leftShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
#if JVET_C0024_ITSKIP
	Int n = 0;
	for (UInt j = 0; j < uiHeight - uiTUSkipHeight; j++)
	{
		n = j*uiWidth;
		for (UInt i = 0; i < uiWidth - uiTUSkipWidth; i++)
		{
#else
	for (Int n = 0; n < numSamplesInBlock; n++)
	{
#endif
		const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
#if JVET_C0024_QTBT
		const Intermediate_Int iCoeffQ = (Intermediate_Int(clipQCoef) * scale) << leftShift;
#else
		const Intermediate_Int iCoeffQ = (Intermediate_Int(clipQCoef)) << leftShift;
#endif

		piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ));
#if JVET_C0024_ITSKIP
		n++;
	}
	memset(piCoef + n, 0, uiTUSkipWidth * sizeof(TCoeff));
		}
memset(piCoef + (uiHeight - uiTUSkipHeight)*uiWidth, 0, uiTUSkipHeight * uiWidth * sizeof(TCoeff));
#else
	}
#endif
}

#elif ITH_C3_DEQUANT==SIMD_FXN
Void rightShiftWithScaleFactor(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int rightShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	TCoeff *srcQCoeff = (TCoeff *)piQCoef;
	TCoeff *dstCoeff = piCoef;

	TCoeff *srcQCoeff_rowPtr = (TCoeff *)piQCoef;
	TCoeff *dstCoeff_rowPtr = piCoef;

	TCoeff *dstCoeff_ref = piCoef;

	__m128i QCoeff_r0_0, QCoeff_r0_1, QCoeff_r1_0, QCoeff_r1_1, QCoeff_r2_0, QCoeff_r2_1, QCoeff_r3_0, QCoeff_r3_1;
	__m128i iCoeffQ_r0_0, iCoeffQ_r0_1, iCoeffQ_r1_0, iCoeffQ_r1_1, iCoeffQ_r2_0, iCoeffQ_r2_1, iCoeffQ_r3_0, iCoeffQ_r3_1;

	__m128i iWHScale_32bx4 = _mm_set1_epi32(scale);
	__m128i inputMinimum_32bx4 = _mm_set1_epi32(inputMinimum);
	__m128i inputMaximum_32bx4 = _mm_set1_epi32(inputMaximum);
	__m128i transformMaximum_32bx4 = _mm_set1_epi32(transformMaximum);
	__m128i transformMinimum_32bx4 = _mm_set1_epi32(transformMinimum);
	__m128i zero_128 = _mm_setzero_si128();



	UInt non_zero_cols = uiWidth - uiTUSkipWidth;
	UInt non_zero_rows = uiHeight - uiTUSkipHeight;

	const Intermediate_Int iAdd = 1 << (rightShift - 1);

	__m128i iAdd_32bx4 = _mm_set1_epi32(iAdd);

	for (Int row = 0; row < non_zero_rows / 4; row++)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));
			QCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride + 4));


			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);

			QCoeff_r3_1 = _mm_min_epi32(QCoeff_r3_1, inputMaximum_32bx4);
			QCoeff_r3_1 = _mm_max_epi32(QCoeff_r3_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);
			iCoeffQ_r3_1 = _mm_mullo_epi32(QCoeff_r3_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r2_1 = _mm_add_epi32(iCoeffQ_r2_1, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);
			iCoeffQ_r3_1 = _mm_add_epi32(iCoeffQ_r3_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r2_1 = _mm_srai_epi32(iCoeffQ_r2_1, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);
			iCoeffQ_r3_1 = _mm_srai_epi32(iCoeffQ_r3_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			iCoeffQ_r3_1 = _mm_min_epi32(iCoeffQ_r3_1, transformMaximum_32bx4);
			iCoeffQ_r3_1 = _mm_max_epi32(iCoeffQ_r3_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride + 4), iCoeffQ_r3_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);
			iCoeffQ_r3_0 = _mm_blend_epi16(iCoeffQ_r3_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		srcQCoeff_rowPtr += 4 * uiSrcStride;
		dstCoeff_rowPtr += 4 * uiDstStride;

		srcQCoeff = srcQCoeff_rowPtr;
		dstCoeff = dstCoeff_rowPtr;
	}

	if ((non_zero_rows & 0x3) == 3)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r2_1 = _mm_add_epi32(iCoeffQ_r2_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r2_1 = _mm_srai_epi32(iCoeffQ_r2_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}
	}

	if ((non_zero_rows & 0x3) == 2)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}

		if ((non_zero_cols % 4) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}
	}

	if ((non_zero_rows & 0x3) == 1)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}
	}

	if (uiTUSkipWidth != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + non_zero_cols + (non_zero_cols & 0x1);
		dstCoeff = dstCoeff_rowPtr;

		//Columns of zeros
		for (Int row = 0; row < ((non_zero_rows) / 2); row++)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}
		}
	}

	if (uiTUSkipHeight != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + (uiHeight - uiTUSkipHeight) * uiDstStride;
		dstCoeff = dstCoeff_rowPtr;

		//Full rows of zero
		for (Int row = 0; row < uiTUSkipHeight / 2; row++)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
			}
		}
	}
}
Void leftShiftWithScaleFactor(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int leftShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	TCoeff *srcQCoeff = (TCoeff *)piQCoef;
	TCoeff *dstCoeff = piCoef;

	TCoeff *srcQCoeff_rowPtr = (TCoeff *)piQCoef;
	TCoeff *dstCoeff_rowPtr = piCoef;

	TCoeff *dstCoeff_ref = piCoef;

	__m128i QCoeff_r0_0, QCoeff_r0_1, QCoeff_r1_0, QCoeff_r1_1, QCoeff_r2_0, QCoeff_r2_1, QCoeff_r3_0, QCoeff_r3_1;
	__m128i iCoeffQ_r0_0, iCoeffQ_r0_1, iCoeffQ_r1_0, iCoeffQ_r1_1, iCoeffQ_r2_0, iCoeffQ_r2_1, iCoeffQ_r3_0, iCoeffQ_r3_1;

	__m128i iWHScale_32bx4 = _mm_set1_epi32(scale);
	__m128i inputMinimum_32bx4 = _mm_set1_epi32(inputMinimum);
	__m128i inputMaximum_32bx4 = _mm_set1_epi32(inputMaximum);
	__m128i transformMaximum_32bx4 = _mm_set1_epi32(transformMaximum);
	__m128i transformMinimum_32bx4 = _mm_set1_epi32(transformMinimum);
	__m128i zero_128 = _mm_setzero_si128();



	UInt non_zero_cols = uiWidth - uiTUSkipWidth;
	UInt non_zero_rows = uiHeight - uiTUSkipHeight;


	for (Int row = 0; row < non_zero_rows / 4; row++)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));
			QCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride + 4));


			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);

			QCoeff_r3_1 = _mm_min_epi32(QCoeff_r3_1, inputMaximum_32bx4);
			QCoeff_r3_1 = _mm_max_epi32(QCoeff_r3_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);
			iCoeffQ_r3_1 = _mm_mullo_epi32(QCoeff_r3_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r2_1 = _mm_slli_epi32(iCoeffQ_r2_1, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);
			iCoeffQ_r3_1 = _mm_slli_epi32(iCoeffQ_r3_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			iCoeffQ_r3_1 = _mm_min_epi32(iCoeffQ_r3_1, transformMaximum_32bx4);
			iCoeffQ_r3_1 = _mm_max_epi32(iCoeffQ_r3_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride + 4), iCoeffQ_r3_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);
			iCoeffQ_r3_0 = _mm_blend_epi16(iCoeffQ_r3_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		srcQCoeff_rowPtr += 4 * uiSrcStride;
		dstCoeff_rowPtr += 4 * uiDstStride;

		srcQCoeff = srcQCoeff_rowPtr;
		dstCoeff = dstCoeff_rowPtr;
	}

	if ((non_zero_rows & 0x3) == 3)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r2_1 = _mm_slli_epi32(iCoeffQ_r2_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}
	}

	if ((non_zero_rows & 0x3) == 2)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}

		if ((non_zero_cols % 4) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}
	}

	if ((non_zero_rows & 0x3) == 1)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);

			//Pointer Update
			srcQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);

			//Pointer Update
			srcQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}
	}

	if (uiTUSkipWidth != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + non_zero_cols + (non_zero_cols & 0x1);
		dstCoeff = dstCoeff_rowPtr;

		//Columns of zeros

		for (Int row = 0; row < ((non_zero_rows) / 2); row++)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}
		}
	}

	if (uiTUSkipHeight != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + (uiHeight - uiTUSkipHeight) * uiDstStride;
		dstCoeff = dstCoeff_rowPtr;

		//Full rows of zero
		for (Int row = 0; row < uiTUSkipHeight / 2; row++)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
			}
		}
	}
}
Void rightShiftWithScaleFactorAndDequantCoef(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int *piDequantCoef, Int rightShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	TCoeff *srcQCoeff = (TCoeff *)piQCoef;
	Int *srcdeQCoeff = piDequantCoef;
	TCoeff *dstCoeff = piCoef;

	TCoeff *srcQCoeff_rowPtr = (TCoeff *)piQCoef;
	Int *srcdeQCoeff_rowPtr = piDequantCoef;
	TCoeff *dstCoeff_rowPtr = piCoef;

	TCoeff *dstCoeff_ref = piCoef;

	__m128i QCoeff_r0_0, QCoeff_r0_1, QCoeff_r1_0, QCoeff_r1_1, QCoeff_r2_0, QCoeff_r2_1, QCoeff_r3_0, QCoeff_r3_1;
	__m128i iCoeffQ_r0_0, iCoeffQ_r0_1, iCoeffQ_r1_0, iCoeffQ_r1_1, iCoeffQ_r2_0, iCoeffQ_r2_1, iCoeffQ_r3_0, iCoeffQ_r3_1;
	__m128i deQCoeff_r0_0, deQCoeff_r0_1, deQCoeff_r1_0, deQCoeff_r1_1, deQCoeff_r2_0, deQCoeff_r2_1, deQCoeff_r3_0, deQCoeff_r3_1;

	__m128i iWHScale_32bx4 = _mm_set1_epi32(scale);
	__m128i inputMinimum_32bx4 = _mm_set1_epi32(inputMinimum);
	__m128i inputMaximum_32bx4 = _mm_set1_epi32(inputMaximum);
	__m128i transformMaximum_32bx4 = _mm_set1_epi32(transformMaximum);
	__m128i transformMinimum_32bx4 = _mm_set1_epi32(transformMinimum);
	__m128i zero_128 = _mm_setzero_si128();



	UInt non_zero_cols = uiWidth - uiTUSkipWidth;
	UInt non_zero_rows = uiHeight - uiTUSkipHeight;

	const Intermediate_Int iAdd = 1 << (rightShift - 1);

	__m128i iAdd_32bx4 = _mm_set1_epi32(iAdd);

	for (Int row = 0; row < non_zero_rows / 4; row++)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));
			QCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride + 4));

			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));
			deQCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride + 4));


			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);

			QCoeff_r3_1 = _mm_min_epi32(QCoeff_r3_1, inputMaximum_32bx4);
			QCoeff_r3_1 = _mm_max_epi32(QCoeff_r3_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);
			iCoeffQ_r3_1 = _mm_mullo_epi32(QCoeff_r3_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r2_1 = _mm_mullo_epi32(iCoeffQ_r2_1, deQCoeff_r2_1);

			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);
			iCoeffQ_r3_1 = _mm_mullo_epi32(iCoeffQ_r3_1, deQCoeff_r3_1);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r2_1 = _mm_add_epi32(iCoeffQ_r2_1, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);
			iCoeffQ_r3_1 = _mm_add_epi32(iCoeffQ_r3_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r2_1 = _mm_srai_epi32(iCoeffQ_r2_1, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);
			iCoeffQ_r3_1 = _mm_srai_epi32(iCoeffQ_r3_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			iCoeffQ_r3_1 = _mm_min_epi32(iCoeffQ_r3_1, transformMaximum_32bx4);
			iCoeffQ_r3_1 = _mm_max_epi32(iCoeffQ_r3_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride + 4), iCoeffQ_r3_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r3_0 = _mm_add_epi32(iCoeffQ_r3_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r3_0 = _mm_srai_epi32(iCoeffQ_r3_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);
			iCoeffQ_r3_0 = _mm_blend_epi16(iCoeffQ_r3_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);
			deQCoeff_r2_0 = _mm_unpacklo_epi64(deQCoeff_r2_0, deQCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);
			deQCoeff_r2_0 = _mm_unpacklo_epi64(deQCoeff_r2_0, deQCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		srcQCoeff_rowPtr += 4 * uiSrcStride;
		srcdeQCoeff_rowPtr += 4 * uiSrcStride;
		dstCoeff_rowPtr += 4 * uiDstStride;

		srcQCoeff = srcQCoeff_rowPtr;
		srcdeQCoeff = srcdeQCoeff_rowPtr;
		dstCoeff = dstCoeff_rowPtr;
	}

	if ((non_zero_rows & 0x3) == 3)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r2_1 = _mm_mullo_epi32(iCoeffQ_r2_1, deQCoeff_r2_1);



			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);
			iCoeffQ_r2_1 = _mm_add_epi32(iCoeffQ_r2_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);
			iCoeffQ_r2_1 = _mm_srai_epi32(iCoeffQ_r2_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r2_0 = _mm_add_epi32(iCoeffQ_r2_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r2_0 = _mm_srai_epi32(iCoeffQ_r2_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}
	}

	if ((non_zero_rows & 0x3) == 2)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);
			iCoeffQ_r1_1 = _mm_add_epi32(iCoeffQ_r1_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);
			iCoeffQ_r1_1 = _mm_srai_epi32(iCoeffQ_r1_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r1_0 = _mm_add_epi32(iCoeffQ_r1_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r1_0 = _mm_srai_epi32(iCoeffQ_r1_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}

		if ((non_zero_cols % 4) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}
	}

	if ((non_zero_rows & 0x3) == 1)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);
			iCoeffQ_r0_1 = _mm_add_epi32(iCoeffQ_r0_1, iAdd_32bx4);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);
			iCoeffQ_r0_1 = _mm_srai_epi32(iCoeffQ_r0_1, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//Add Rounding Factor iAdd
			iCoeffQ_r0_0 = _mm_add_epi32(iCoeffQ_r0_0, iAdd_32bx4);

			//Right shift
			iCoeffQ_r0_0 = _mm_srai_epi32(iCoeffQ_r0_0, rightShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}
	}

	if (uiTUSkipWidth != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + non_zero_cols + (non_zero_cols & 0x1);
		dstCoeff = dstCoeff_rowPtr;

		//Columns of zeros
		for (Int row = 0; row < ((non_zero_rows) / 2); row++)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}
		}
	}

	if (uiTUSkipHeight != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + (uiHeight - uiTUSkipHeight) * uiDstStride;
		dstCoeff = dstCoeff_rowPtr;

		//Full rows of zero
		for (Int row = 0; row < uiTUSkipHeight / 2; row++)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
			}
		}
	}
}
Void leftShiftWithScaleFactorAndDequantCoef(UInt uiWidth, UInt uiHeight, UInt uiTUSkipWidth, UInt uiTUSkipHeight, const TCoeff *piQCoef, UInt uiSrcStride, TCoeff *piCoef, UInt uiDstStride, Int *piDequantCoef, Int leftShift, Int scale, Intermediate_Int inputMinimum, Intermediate_Int inputMaximum, TCoeff transformMinimum, TCoeff transformMaximum)
{
	TCoeff *srcQCoeff = (TCoeff *)piQCoef;
	Int *srcdeQCoeff = piDequantCoef;
	TCoeff *dstCoeff = piCoef;

	TCoeff *srcQCoeff_rowPtr = (TCoeff *)piQCoef;
	Int *srcdeQCoeff_rowPtr = piDequantCoef;
	TCoeff *dstCoeff_rowPtr = piCoef;

	TCoeff *dstCoeff_ref = piCoef;

	__m128i QCoeff_r0_0, QCoeff_r0_1, QCoeff_r1_0, QCoeff_r1_1, QCoeff_r2_0, QCoeff_r2_1, QCoeff_r3_0, QCoeff_r3_1;
	__m128i iCoeffQ_r0_0, iCoeffQ_r0_1, iCoeffQ_r1_0, iCoeffQ_r1_1, iCoeffQ_r2_0, iCoeffQ_r2_1, iCoeffQ_r3_0, iCoeffQ_r3_1;
	__m128i deQCoeff_r0_0, deQCoeff_r0_1, deQCoeff_r1_0, deQCoeff_r1_1, deQCoeff_r2_0, deQCoeff_r2_1, deQCoeff_r3_0, deQCoeff_r3_1;

	__m128i iWHScale_32bx4 = _mm_set1_epi32(scale);
	__m128i inputMinimum_32bx4 = _mm_set1_epi32(inputMinimum);
	__m128i inputMaximum_32bx4 = _mm_set1_epi32(inputMaximum);
	__m128i transformMaximum_32bx4 = _mm_set1_epi32(transformMaximum);
	__m128i transformMinimum_32bx4 = _mm_set1_epi32(transformMinimum);
	__m128i zero_128 = _mm_setzero_si128();

	UInt non_zero_cols = uiWidth - uiTUSkipWidth;
	UInt non_zero_rows = uiHeight - uiTUSkipHeight;

	for (Int row = 0; row < non_zero_rows / 4; row++)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));
			QCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride + 4));

			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));
			deQCoeff_r3_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride + 4));


			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);

			QCoeff_r3_1 = _mm_min_epi32(QCoeff_r3_1, inputMaximum_32bx4);
			QCoeff_r3_1 = _mm_max_epi32(QCoeff_r3_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);
			iCoeffQ_r3_1 = _mm_mullo_epi32(QCoeff_r3_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r2_1 = _mm_mullo_epi32(iCoeffQ_r2_1, deQCoeff_r2_1);

			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);
			iCoeffQ_r3_1 = _mm_mullo_epi32(iCoeffQ_r3_1, deQCoeff_r3_1);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r2_1 = _mm_slli_epi32(iCoeffQ_r2_1, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);
			iCoeffQ_r3_1 = _mm_slli_epi32(iCoeffQ_r3_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			iCoeffQ_r3_1 = _mm_min_epi32(iCoeffQ_r3_1, transformMaximum_32bx4);
			iCoeffQ_r3_1 = _mm_max_epi32(iCoeffQ_r3_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride + 4), iCoeffQ_r3_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r3_0 = _mm_min_epi32(QCoeff_r3_0, inputMaximum_32bx4);
			QCoeff_r3_0 = _mm_max_epi32(QCoeff_r3_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r3_0 = _mm_mullo_epi32(QCoeff_r3_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r3_0 = _mm_mullo_epi32(iCoeffQ_r3_0, deQCoeff_r3_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r3_0 = _mm_slli_epi32(iCoeffQ_r3_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r3_0 = _mm_min_epi32(iCoeffQ_r3_0, transformMaximum_32bx4);
			iCoeffQ_r3_0 = _mm_max_epi32(iCoeffQ_r3_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);
			iCoeffQ_r3_0 = _mm_blend_epi16(iCoeffQ_r3_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 3 * uiDstStride), iCoeffQ_r3_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);
			deQCoeff_r2_0 = _mm_unpacklo_epi64(deQCoeff_r2_0, deQCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 3 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);
			QCoeff_r2_0 = _mm_unpacklo_epi64(QCoeff_r2_0, QCoeff_r3_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r3_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 3 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);
			deQCoeff_r2_0 = _mm_unpacklo_epi64(deQCoeff_r2_0, deQCoeff_r3_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + 3 * uiDstStride), _mm_srli_si128(iCoeffQ_r2_0, 8));
		}

		srcQCoeff_rowPtr += 4 * uiSrcStride;
		srcdeQCoeff_rowPtr += 4 * uiSrcStride;
		dstCoeff_rowPtr += 4 * uiDstStride;

		srcQCoeff = srcQCoeff_rowPtr;
		srcdeQCoeff = srcdeQCoeff_rowPtr;
		dstCoeff = dstCoeff_rowPtr;
	}

	if ((non_zero_rows & 0x3) == 3)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));
			QCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));
			deQCoeff_r2_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			QCoeff_r2_1 = _mm_min_epi32(QCoeff_r2_1, inputMaximum_32bx4);
			QCoeff_r2_1 = _mm_max_epi32(QCoeff_r2_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);
			iCoeffQ_r2_1 = _mm_mullo_epi32(QCoeff_r2_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);
			iCoeffQ_r2_1 = _mm_mullo_epi32(iCoeffQ_r2_1, deQCoeff_r2_1);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);
			iCoeffQ_r2_1 = _mm_slli_epi32(iCoeffQ_r2_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r2_1 = _mm_min_epi32(iCoeffQ_r2_1, transformMaximum_32bx4);
			iCoeffQ_r2_1 = _mm_max_epi32(iCoeffQ_r2_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride + 4), iCoeffQ_r2_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + 2 * uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r2_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 2 * uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r2_0 = _mm_min_epi32(QCoeff_r2_0, inputMaximum_32bx4);
			QCoeff_r2_0 = _mm_max_epi32(QCoeff_r2_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r2_0 = _mm_mullo_epi32(QCoeff_r2_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r2_0 = _mm_mullo_epi32(iCoeffQ_r2_0, deQCoeff_r2_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r2_0 = _mm_slli_epi32(iCoeffQ_r2_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r2_0 = _mm_min_epi32(iCoeffQ_r2_0, transformMaximum_32bx4);
			iCoeffQ_r2_0 = _mm_max_epi32(iCoeffQ_r2_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);
			iCoeffQ_r2_0 = _mm_blend_epi16(iCoeffQ_r2_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
			_mm_storel_epi64((__m128i *) (dstCoeff + 2 * uiDstStride), iCoeffQ_r2_0);
		}
	}

	if ((non_zero_rows & 0x3) == 2)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));
			QCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + 4));

			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));
			deQCoeff_r1_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride + 4));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			QCoeff_r1_1 = _mm_min_epi32(QCoeff_r1_1, inputMaximum_32bx4);
			QCoeff_r1_1 = _mm_max_epi32(QCoeff_r1_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);
			iCoeffQ_r1_1 = _mm_mullo_epi32(QCoeff_r1_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);
			iCoeffQ_r1_1 = _mm_mullo_epi32(iCoeffQ_r1_1, deQCoeff_r1_1);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);
			iCoeffQ_r1_1 = _mm_slli_epi32(iCoeffQ_r1_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			iCoeffQ_r1_1 = _mm_min_epi32(iCoeffQ_r1_1, transformMaximum_32bx4);
			iCoeffQ_r1_1 = _mm_max_epi32(iCoeffQ_r1_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride + 4), iCoeffQ_r1_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}

		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r1_0 = _mm_min_epi32(QCoeff_r1_0, inputMaximum_32bx4);
			QCoeff_r1_0 = _mm_max_epi32(QCoeff_r1_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r1_0 = _mm_mullo_epi32(QCoeff_r1_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r1_0 = _mm_mullo_epi32(iCoeffQ_r1_0, deQCoeff_r1_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r1_0 = _mm_slli_epi32(iCoeffQ_r1_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r1_0 = _mm_min_epi32(iCoeffQ_r1_0, transformMaximum_32bx4);
			iCoeffQ_r1_0 = _mm_max_epi32(iCoeffQ_r1_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);
			iCoeffQ_r1_0 = _mm_blend_epi16(iCoeffQ_r1_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), iCoeffQ_r1_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}

		if ((non_zero_cols % 4) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcQCoeff + uiSrcStride));

			QCoeff_r0_0 = _mm_unpacklo_epi64(QCoeff_r0_0, QCoeff_r1_0);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r1_0 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			deQCoeff_r0_0 = _mm_unpacklo_epi64(deQCoeff_r0_0, deQCoeff_r1_0);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
			_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), _mm_srli_si128(iCoeffQ_r0_0, 8));
		}
	}

	if ((non_zero_rows & 0x3) == 1)
	{
		for (Int col = 0; col < (non_zero_cols / 8); col++)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);
			QCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcQCoeff + 4));

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);
			deQCoeff_r0_1 = _mm_loadu_si128((__m128i *) (srcdeQCoeff + uiSrcStride));

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			QCoeff_r0_1 = _mm_min_epi32(QCoeff_r0_1, inputMaximum_32bx4);
			QCoeff_r0_1 = _mm_max_epi32(QCoeff_r0_1, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_1 = _mm_mullo_epi32(QCoeff_r0_1, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);
			iCoeffQ_r0_1 = _mm_mullo_epi32(iCoeffQ_r0_1, deQCoeff_r0_1);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);
			iCoeffQ_r0_1 = _mm_slli_epi32(iCoeffQ_r0_1, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_1 = _mm_min_epi32(iCoeffQ_r0_1, transformMaximum_32bx4);
			iCoeffQ_r0_1 = _mm_max_epi32(iCoeffQ_r0_1, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) (dstCoeff), iCoeffQ_r0_0);
			_mm_storeu_si128((__m128i *) (dstCoeff + 4), iCoeffQ_r0_1);

			//Pointer Update
			srcQCoeff += 8;
			srcdeQCoeff += 8;
			dstCoeff += 8;
		}


		//If leftover is greater than 3
		if ((non_zero_cols & 0x7) > 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);

			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);

			//Pointer Update
			srcQCoeff += 4;
			srcdeQCoeff += 4;
			dstCoeff += 4;
		}

		if ((non_zero_cols & 0x3) == 3)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xc0);

			_mm_storeu_si128((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 2)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);


			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}

		if ((non_zero_cols & 0x3) == 1)
		{
			//Loading QCoeff values
			QCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcQCoeff);

			//Loading deQCoeff values
			deQCoeff_r0_0 = _mm_loadu_si128((__m128i *) srcdeQCoeff);

			//Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n])
			QCoeff_r0_0 = _mm_min_epi32(QCoeff_r0_0, inputMaximum_32bx4);
			QCoeff_r0_0 = _mm_max_epi32(QCoeff_r0_0, inputMinimum_32bx4);

			//Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale
			iCoeffQ_r0_0 = _mm_mullo_epi32(QCoeff_r0_0, iWHScale_32bx4);
			iCoeffQ_r0_0 = _mm_mullo_epi32(iCoeffQ_r0_0, deQCoeff_r0_0);

			//leftShift
			iCoeffQ_r0_0 = _mm_slli_epi32(iCoeffQ_r0_0, leftShift);

			//Clip3<Intermediate_Int>(transformMinimum, transformMaximum, iCoeffQ)
			iCoeffQ_r0_0 = _mm_min_epi32(iCoeffQ_r0_0, transformMaximum_32bx4);
			iCoeffQ_r0_0 = _mm_max_epi32(iCoeffQ_r0_0, transformMinimum_32bx4);

			iCoeffQ_r0_0 = _mm_blend_epi16(iCoeffQ_r0_0, zero_128, 0xcc);

			//store
			_mm_storel_epi64((__m128i *) dstCoeff, iCoeffQ_r0_0);
		}
	}

	if (uiTUSkipWidth != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + non_zero_cols + (non_zero_cols & 0x1);
		dstCoeff = dstCoeff_rowPtr;

		//Columns of zeros
		for (Int row = 0; row < ((non_zero_rows) / 2); row++)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			Int num_zero_cols = (uiTUSkipWidth & 0xfffe);
			for (Int col = (num_zero_cols / 4); col > 0; col--)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiSrcStride), zero_128);
				dstCoeff += 4;
			}

			if ((num_zero_cols & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiSrcStride), zero_128);
			}
		}
	}

	if (uiTUSkipHeight != 0)
	{
		dstCoeff_rowPtr = dstCoeff_ref + (uiHeight - uiTUSkipHeight) * uiDstStride;
		dstCoeff = dstCoeff_rowPtr;

		//Full rows of zero
		for (Int row = 0; row < uiTUSkipHeight / 2; row++)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				_mm_storeu_si128((__m128i *) (dstCoeff + uiDstStride), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
				_mm_storel_epi64((__m128i *) (dstCoeff + uiDstStride), zero_128);
			}

			dstCoeff_rowPtr += 2 * uiDstStride;
			dstCoeff = dstCoeff_rowPtr;
		}

		if ((uiTUSkipHeight & 0x1) == 1)
		{
			for (Int col = 0; col < uiWidth; col += 4)
			{
				_mm_storeu_si128((__m128i *) (dstCoeff), zero_128);
				dstCoeff += 4;
			}

			if ((uiWidth & 0x3) == 2)
			{
				_mm_storel_epi64((__m128i *) (dstCoeff), zero_128);
			}
		}
	}
}
#endif

Void TComTrQuant::xDeQuant(       TComTU        &rTu,
                            const TCoeff       * pSrc,
                                  TCoeff       * pDes,
                            const ComponentID    compID,
                            const QpParam       &cQP )
{
  assert(compID<MAX_NUM_COMPONENT);

        TComDataCU          *pcCU               = rTu.getCU();
  const UInt                 uiAbsPartIdx       = rTu.GetAbsPartIdxTU();
  const TComRectangle       &rect               = rTu.getRect(compID);
  const UInt                 uiWidth            = rect.width;
  const UInt                 uiHeight           = rect.height;
  const TCoeff        *const piQCoef            = pSrc;
        TCoeff        *const piCoef             = pDes;
  const UInt                 uiLog2TrSize       = rTu.GetEquivalentLog2TrSize(compID);
  const UInt                 numSamplesInBlock  = uiWidth*uiHeight;
  const Int                  maxLog2TrDynamicRange  = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));
  const TCoeff               transformMinimum   = -(1 << maxLog2TrDynamicRange);
  const TCoeff               transformMaximum   =  (1 << maxLog2TrDynamicRange) - 1;
  const Bool                 enableScalingLists = getUseScalingList(uiWidth, uiHeight, (pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0));
  const Int                  scalingListType    = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
#if O0043_BEST_EFFORT_DECODING
  const Int                  channelBitDepth    = pcCU->getSlice()->getSPS()->getStreamBitDepth(toChannelType(compID));
#else
  const Int                  channelBitDepth    = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#endif

  assert (scalingListType < SCALING_LIST_NUM);
  assert ( uiWidth <= m_uiMaxTrSize );

  // Represents scaling through forward transform
  const Bool bClipTransformShiftTo0 = (pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0) && pcCU->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag();
  const Int  originalTransformShift = getTransformShift(channelBitDepth, uiLog2TrSize, maxLog2TrDynamicRange);
  const Int  iTransformShift        = bClipTransformShiftTo0 ? std::max<Int>(0, originalTransformShift) : originalTransformShift;

  const Int QP_per = cQP.per;
  const Int QP_rem = cQP.rem;

#if JVET_C0024_QTBT
  Int rightShift = (IQUANT_SHIFT - (iTransformShift + QP_per)) + (enableScalingLists ? LOG2_SCALING_LIST_NEUTRAL_VALUE : 0);
  Int iWHScale = 1;
  Int iWHShift = 0;
  if ((g_aucConvertToBit[ uiWidth ] + g_aucConvertToBit[ uiHeight ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iWHShift += 8;
    iWHScale = 181;
  }
#else
#if PVC_AQS
  // Note: now, enableScalingLists is false, so we did not add code to this part
  if (enableScalingLists)
  {
    fprintf(stderr, "\n Error: AQS dequant is not fully implemented when using scaling list");
    exit(-1);
  }

  const Int rightShift = (IQUANT_SHIFT - (iTransformShift + QP_per)) + (enableScalingLists ? LOG2_SCALING_LIST_NEUTRAL_VALUE : 0) + (compID == COMPONENT_Y ? 14 : 0);
#else
  const Int rightShift = (IQUANT_SHIFT - (iTransformShift + QP_per)) + (enableScalingLists ? LOG2_SCALING_LIST_NEUTRAL_VALUE : 0);
#endif
#endif

  if(enableScalingLists)
  {
    //from the dequantisation equation:
    //iCoeffQ                         = ((Intermediate_Int(clipQCoef) * piDequantCoef[deQuantIdx]) + iAdd ) >> rightShift
    //(sizeof(Intermediate_Int) * 8)  =              inputBitDepth    +    dequantCoefBits                   - rightShift
    const UInt             dequantCoefBits     = 1 + IQUANT_SHIFT + SCALING_LIST_BITS;
    const UInt             targetInputBitDepth = std::min<UInt>((maxLog2TrDynamicRange + 1), (((sizeof(Intermediate_Int) * 8) + rightShift) - dequantCoefBits));

    const Intermediate_Int inputMinimum        = -(1 << (targetInputBitDepth - 1));
    const Intermediate_Int inputMaximum        =  (1 << (targetInputBitDepth - 1)) - 1;

#if JVET_C0024_QTBT
    const UInt uiLog2TrWidth = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
    const UInt uiLog2TrHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
    Int *piDequantCoef = getDequantCoeff(scalingListType,QP_rem,uiLog2TrWidth-1, uiLog2TrHeight-1);
#else
    Int *piDequantCoef = getDequantCoeff(scalingListType,QP_rem,uiLog2TrSize-2);
#endif

#if JVET_C0024_QTBT
  rightShift += iWHShift;
#endif
    if(rightShift > 0)
    {
#if ITH_C2_ZR_ZC_SL
#if JVET_C0024_QTBT
		rightShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, pcCU->getTUSkipWidth(compID, uiAbsPartIdx), pcCU->getTUSkipHeight(compID, uiAbsPartIdx), piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, rightShift, iWHScale, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#else
		rightShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, pcCU->getTUSkipWidth(compID, uiAbsPartIdx), pcCU->getTUSkipHeight(compID, uiAbsPartIdx), piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, rightShift, 1, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#endif
#else
#if JVET_C0024_QTBT
		rightShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, 0, 0, piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, rightShift, iWHScale, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#else
		rightShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, 0, 0, piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, rightShift, 1, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#endif
#endif
//      const Intermediate_Int iAdd = 1 << (rightShift - 1);
//
//      for( Int n = 0; n < numSamplesInBlock; n++ )
//      {
//        const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
//#if JVET_C0024_QTBT
//        const Intermediate_Int iCoeffQ   = ((Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale) + iAdd ) >> rightShift;
//#else
//        const Intermediate_Int iCoeffQ   = ((Intermediate_Int(clipQCoef) * piDequantCoef[n]) + iAdd ) >> rightShift;
//#endif
//
//        piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
//      }
    }
    else
    {
#if ITH_C2_ZR_ZC_SL
#if JVET_C0024_QTBT
		leftShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, pcCU->getTUSkipWidth(compID, uiAbsPartIdx), pcCU->getTUSkipHeight(compID, uiAbsPartIdx), piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, -rightShift, iWHScale, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#else
		leftShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, pcCU->getTUSkipWidth(compID, uiAbsPartIdx), pcCU->getTUSkipHeight(compID, uiAbsPartIdx), piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, -rightShift, 1, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#endif
#else
#if JVET_C0024_QTBT
		leftShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, 0, 0, piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, -rightShift, iWHScale, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#else
		leftShiftWithScaleFactorAndDequantCoef(uiWidth, uiHeight, 0, 0, piQCoef, uiWidth, piCoef, uiWidth, piDequantCoef, -rightShift, 1, inputMinimum, inputMaximum, transformMinimum, transformMaximum);
#endif
#endif
//      const Int leftShift = -rightShift;
//
//      for( Int n = 0; n < numSamplesInBlock; n++ )
//      {
//        const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
//#if JVET_C0024_QTBT
//        const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * piDequantCoef[n] * iWHScale) << leftShift;
//#else
//        const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * piDequantCoef[n]) << leftShift;
//#endif
//
//        piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
//      }
    }
  }
  else
  {
#if JVET_C0024_QTBT
    const Int scale     =  g_invQuantScales[QP_rem] * iWHScale;
#else
    const Int scale     =  g_invQuantScales[QP_rem];
#endif
    const Int scaleBits =     (IQUANT_SHIFT + 1)   ;

    //from the dequantisation equation:
    //iCoeffQ                         = Intermediate_Int((Int64(clipQCoef) * scale + iAdd) >> rightShift);
    //(sizeof(Intermediate_Int) * 8)  =                    inputBitDepth   + scaleBits      - rightShift
    const UInt             targetInputBitDepth = std::min<UInt>((maxLog2TrDynamicRange + 1), (((sizeof(Intermediate_Int) * 8) + rightShift) - scaleBits));
    const Intermediate_Int inputMinimum        = -(1 << (targetInputBitDepth - 1));
    const Intermediate_Int inputMaximum        =  (1 << (targetInputBitDepth - 1)) - 1;

#if JVET_C0024_QTBT
    rightShift += iWHShift;
#endif
    if (rightShift > 0)
    {
      const Intermediate_Int iAdd = 1 << (rightShift - 1);

#if JVET_C0024_ITSKIP
      Int n=0;
      for (UInt j=0; j<uiHeight-pcCU->getTUSkipHeight(compID, uiAbsPartIdx); j++)
      {
        n = j*uiWidth;
        for (UInt i=0; i<uiWidth-pcCU->getTUSkipWidth(compID, uiAbsPartIdx); i++)
        {
#else
      for( Int n = 0; n < numSamplesInBlock; n++ )
      {
#endif
        const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));

#if PVC_AQS  // de-quantization
        Intermediate_Int iCoeffQ;
        if (compID == COMPONENT_Y)
        {
          Int iEstQSFactorInvShift14 = ((1 << (14 + ESM_SHIFT)) + (m_iEstQSFactor >> 1)) / m_iEstQSFactor;
          iCoeffQ = ((Intermediate_Int) clipQCoef * scale * (Int64) iEstQSFactorInvShift14 + iAdd) >> rightShift;
        }
        else
        {
          iCoeffQ = (Intermediate_Int(clipQCoef) * scale + iAdd) >> rightShift;
        }
#else
        const Intermediate_Int iCoeffQ = (Intermediate_Int(clipQCoef) * scale + iAdd) >> rightShift;
#endif

        piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
#if JVET_C0024_ITSKIP
        n++;
      }
      memset(piCoef+n, 0, pcCU->getTUSkipWidth(compID, uiAbsPartIdx)*sizeof(TCoeff));
    }
    memset(piCoef+(uiHeight-pcCU->getTUSkipHeight(compID, uiAbsPartIdx))*uiWidth, 0, pcCU->getTUSkipHeight(compID, uiAbsPartIdx)*uiWidth * sizeof(TCoeff));
#else
      }
#endif
    }
    else
    {
      const Int leftShift = -rightShift;

#if JVET_C0024_ITSKIP
      Int n=0;
      for (UInt j=0; j<uiHeight-pcCU->getTUSkipHeight(compID, uiAbsPartIdx); j++)
      {
        n = j*uiWidth;
        for (UInt i=0; i<uiWidth-pcCU->getTUSkipWidth(compID, uiAbsPartIdx); i++)
        {
#else
      for( Int n = 0; n < numSamplesInBlock; n++ )
      {
#endif
        const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, piQCoef[n]));
        const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * scale) << leftShift;

        piCoef[n] = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
#if JVET_C0024_ITSKIP
        n++;
      }
        memset(piCoef+n, 0, pcCU->getTUSkipWidth(compID, uiAbsPartIdx)*sizeof(TCoeff));
    }
      memset(piCoef+(uiHeight-pcCU->getTUSkipHeight(compID, uiAbsPartIdx))*uiWidth, 0, pcCU->getTUSkipHeight(compID, uiAbsPartIdx)*uiWidth * sizeof(TCoeff));
#else
  }
#endif
}
  }
}


Void TComTrQuant::init(   UInt  uiMaxTrSize,
#if VCEG_AZ08_USE_KLT
                          UInt  uiUseKLT,
#endif
                          Bool  bUseRDOQ,
                          Bool  bUseRDOQTS,
#if T0196_SELECTIVE_RDOQ
                          Bool  useSelectiveRDOQ,
#endif
                          Bool  bEnc,
                          Bool  useTransformSkipFast
#if ADAPTIVE_QP_SELECTION
                        , Bool bUseAdaptQpSelect
#endif
                       )
{
  m_uiMaxTrSize  = uiMaxTrSize;

  m_bEnc         = bEnc;
  m_useRDOQ      = bUseRDOQ;
  m_useRDOQTS    = bUseRDOQTS;
#if T0196_SELECTIVE_RDOQ
  m_useSelectiveRDOQ = useSelectiveRDOQ;
#endif
#if ADAPTIVE_QP_SELECTION
  m_bUseAdaptQpSelect = bUseAdaptQpSelect;
#endif
  m_useTransformSkipFast = useTransformSkipFast;
#if VCEG_AZ08_USE_KLT
  m_useKLT = uiUseKLT;
  if (uiUseKLT)
  {
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
      Int useIntraKLT = m_useKLT & 1;
      Int useInterKLT = (m_useKLT >> 1) & 1;
#endif
#endif
#if VCEG_AZ08_KLT_COMMON
      UInt blkSize;
      UInt uiDim;
      m_tempLibFast.init(MAX_CANDI_NUM);
      for (UInt i = 0; i < MAX_CANDI_NUM; i++)
      {
        if( m_pData[i] == NULL )
        {
          m_pData[i] = new TrainDataType[MAX_1DTRANS_LEN];
}
      }
#if VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY
      for (UInt i = 0; i < MAX_1DTRANS_LEN; i++)
      {
          m_pDataT[i] = new TrainDataType[MAX_CANDI_NUM];
      }
#endif
      if( m_pppTarPatch == NULL )
      {
        m_pppTarPatch = new Pel**[USE_MORE_BLOCKSIZE_DEPTH_MAX];
        for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
        {
          blkSize = g_uiDepth2Width[uiDepth];
#if VCEG_AZ08_USE_KLT
          UInt tempSize = 0;
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
          if(useIntraKLT && useInterKLT)
          {
            tempSize = max(g_uiDepth2IntraTempSize[uiDepth], g_uiDepth2InterTempSize[uiDepth]);
          }
          else if(useIntraKLT && !useInterKLT)
          {
            tempSize = g_uiDepth2IntraTempSize[uiDepth];
          }
          else if (!useIntraKLT && useInterKLT)
          {
            tempSize = g_uiDepth2InterTempSize[uiDepth];
          }
#endif
#if VCEG_AZ08_INTRA_KLT && !VCEG_AZ08_INTER_KLT
          tempSize = g_uiDepth2IntraTempSize[uiDepth];
#endif
#if !VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
          tempSize = g_uiDepth2InterTempSize[uiDepth];
#endif
#else
#if VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
          UInt tempSize = max(g_uiDepth2IntraTempSize[uiDepth], g_uiDepth2InterTempSize[uiDepth]);
#endif
#if VCEG_AZ08_INTRA_KLT && !VCEG_AZ08_INTER_KLT
          UInt tempSize = g_uiDepth2IntraTempSize[uiDepth];
#endif
#if !VCEG_AZ08_INTRA_KLT && VCEG_AZ08_INTER_KLT
          UInt tempSize = g_uiDepth2InterTempSize[uiDepth];
#endif
#endif
          UInt patchSize = blkSize + tempSize;
          m_pppTarPatch[uiDepth] = new Pel *[patchSize];
          for (UInt uiRow = 0; uiRow < patchSize; uiRow++)
          {
            m_pppTarPatch[uiDepth][uiRow] = new Pel[patchSize];
          }
        }
      }

      if( m_pCovMatrix == NULL )
      {
        m_pCovMatrix = new covMatrixType*[USE_MORE_BLOCKSIZE_DEPTH_MAX];
        for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
        {
          blkSize = g_uiDepth2Width[uiDepth];
          uiDim = blkSize*blkSize;
          UInt uiMatrixDim = uiDim*uiDim;
          m_pCovMatrix[uiDepth] = new covMatrixType[uiMatrixDim];
        }
      }
#if VCEG_AZ08_FAST_DERIVE_KLT
      if( m_pppdTmpEigenVector == NULL )
      {
        blkSize = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
        uiDim = blkSize*blkSize;
        m_pppdTmpEigenVector = new EigenType *[uiDim];
        for (UInt k = 0; k < uiDim; k++)
        {
          m_pppdTmpEigenVector[k] = new EigenType[uiDim];
        }
      }
#endif
      if( m_pppdEigenVector == NULL )
      {
        m_pppdEigenVector = new EigenType**[USE_MORE_BLOCKSIZE_DEPTH_MAX];
        for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
        {
          blkSize = g_uiDepth2Width[uiDepth];
          uiDim = blkSize*blkSize;
          m_pppdEigenVector[uiDepth] = new EigenType *[uiDim];
          for (UInt k = 0; k < uiDim; k++)
          {
            m_pppdEigenVector[uiDepth][k] = new EigenType[uiDim];
          }
        }
      }
      if( m_pppsEigenVector == NULL )
      {
        m_pppsEigenVector = new Short**[USE_MORE_BLOCKSIZE_DEPTH_MAX];
        for (UInt uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
        {
          blkSize = g_uiDepth2Width[uiDepth];
          uiDim = blkSize*blkSize;
          m_pppsEigenVector[uiDepth] = new Short *[uiDim];
          for (UInt k = 0; k < uiDim; k++)
          {
            m_pppsEigenVector[uiDepth][k] = new Short[uiDim];
          }
        }
      }
#endif
#if VCEG_AZ08_USE_KLT
  }
#endif
}
#if VCEG_AZ05_ROT_TR
// ====================================================================================================================
// VCEG_AZ05_ROT_TR
// ====================================================================================================================
const short g_ROT_U[4][6] =
{//----------horizontal-----    
  { -7,  -7,   5,  -4,  -2,  -2, },
  {  4, -22,   9,  -7,  -4,   5, },
  { -3,  -5,  -3,   2, -29,   5, },
  {-10,  -6,  -4,  -1, -29,   6, },
};    
const short g_ROT_V[4][6] =
{//----------horizontal-----      
  { 14,  13, -10,   8,   3,   4, },
  { -8,  30, -17,  13,   7, -10, },
  {  5,   9,   6,  -5,  32, -10, },
  { 18,  11,   8,   1,  32, -12, },
};    
#define LIFT_PRECISION   5
///4x4
Void TComTrQuant::InvRotTransform4I(  Int* matrix, UChar index )
{
 Int temp[16];
  Int n=0;
  // Rot process
  Int liftPrecision = LIFT_PRECISION;
  
  Int dum;
  Int A0, B1,  C0,  D1,  F1,  G0;
  
  for (n=0; n<4; n++)
  { // vertical ROT
   temp[n+12] = matrix[n+12];
     
    G0 = matrix[n   ] - xMult(g_ROT_U[index][2] * matrix[n+ 4], liftPrecision);
    F1 = matrix[n+4 ] - xMult(g_ROT_V[index][2]*G0, liftPrecision);
    C0 = G0        - xMult(g_ROT_U[index][2]*F1, liftPrecision);

    D1      = F1        - xMult(g_ROT_U[index][1] * matrix[n+8], liftPrecision);
    temp[n+8] = matrix[n+8] - xMult(g_ROT_V[index][1] * D1, liftPrecision);
    B1      = D1        - xMult(g_ROT_U[index][1] * temp[n+8], liftPrecision);

    if (index == 1)
    {
      dum = B1;
      B1 = temp[n+8];
      temp[n+8] = -dum;
    } // if index == 1

    A0      = C0 - xMult(g_ROT_U[index][0] * B1, liftPrecision);
    temp[n+ 4] = B1 - xMult(g_ROT_V[index][0] * A0, liftPrecision);
    temp[n]    = A0 - xMult(g_ROT_U[index][0] * temp[n+4], liftPrecision);
  } // for n, vertical ROT

  for (n=0; n<16; n+=4)
  { // horizontal ROT

    matrix[n+3] = temp[n+3];
    G0 = temp[n  ] - xMult(g_ROT_U[index][5] * temp[n+1], liftPrecision);
    F1 = temp[n+1] - xMult(g_ROT_V[index][5] * G0, liftPrecision);
    C0 = G0   - xMult(g_ROT_U[index][5] * F1, liftPrecision);

    D1 = F1       - xMult(g_ROT_U[index][4] * temp[n+2], liftPrecision);
    matrix[n+2] = temp[n+2] - xMult(g_ROT_V[index][4] * D1, liftPrecision);
    B1 = D1       - xMult(g_ROT_U[index][4] * matrix[n+2], liftPrecision);

     if ((index == 2) || (index == 3))
    {
      dum = B1; 
      B1 = matrix[n+2];
      matrix[n+2] = -dum;
    }// if index = 2 or 3

    A0     = C0 - xMult(g_ROT_U[index][3] * B1, liftPrecision);
    matrix[n+1] = B1 - xMult(g_ROT_V[index][3] * A0, liftPrecision);
    matrix[n  ] = A0 - xMult(g_ROT_U[index][3] * matrix[n+1], liftPrecision);
  } // for n, horizontal ROT
}

Void TComTrQuant::RotTransform4I( Int* matrix, UChar index )
{

  Int temp[16];
  Int n = 0;

  Int liftPrecision = LIFT_PRECISION;
  Int  dum;
  Int A0,  B1,  C0,  D1,  F1,  G0;


  for (n=0; n<16; n+=4)
  { // horizontal ROT
    
    A0 = matrix[n  ] + xMult(g_ROT_U[index][3] * matrix[n+1], liftPrecision);
      B1 = xMult(g_ROT_V[index][3]*A0, liftPrecision) + matrix[n+1];
    C0 = A0     + xMult(g_ROT_U[index][3] * B1, liftPrecision);
    if ((index == 2) || (index == 3))
    {
      dum = B1; 
      B1 = -matrix[n+2];
      matrix[n+2] = dum;
    }// if index = 2 or 3
    D1   = B1     + xMult(g_ROT_U[index][4] * matrix[n+2], liftPrecision);
    temp[n+2] = xMult(g_ROT_V[index][4]*D1, liftPrecision) + matrix[n+2];
    F1   = D1     + xMult(g_ROT_U[index][4] * temp[n+2], liftPrecision);

    G0   = C0     + xMult(g_ROT_U[index][5] * F1, liftPrecision);
    temp[n+1] = xMult(g_ROT_V[index][5]*G0, liftPrecision) + F1;
    temp[n  ] = G0     + xMult(g_ROT_U[index][5] * temp[n+1], liftPrecision);

    temp[n+3] = matrix[n+3];
  } // for n, horizontal ROT

  for (n=0; n<4; n++)
  { // vertical ROT
    A0   = temp[n]       + xMult(g_ROT_U[index][0] * temp[n+4 ], liftPrecision);
    B1   = xMult(g_ROT_V[index][0]*A0, liftPrecision) + temp[n+4];
    C0   = A0         + xMult(g_ROT_U[index][0] * B1, liftPrecision);

    if (index == 1)
    {
      dum = B1;
      B1 = -temp[n+8];
      temp[n+8] = dum;
    } // if index == 1
    D1   = B1         + xMult(g_ROT_U[index][1] * temp[n+8], liftPrecision);
    matrix[n+8] = xMult(g_ROT_V[index][1]*D1, liftPrecision) + temp[n+8];
    F1   = D1         + xMult(g_ROT_U[index][1] * matrix[n+8], liftPrecision);

    G0   = C0         + xMult(g_ROT_U[index][2] * F1, liftPrecision);
    matrix[n+4 ] = xMult(g_ROT_V[index][2]*G0, liftPrecision) + F1;
    matrix[n   ] = G0    + xMult(g_ROT_U[index][2] * matrix[n+ 4], liftPrecision);

    matrix[n+12] = temp[n+12];
  } // for n, vertical ROT
}
#elif COM16_C1044_NSST

#if JVET_D0120_NSST_IMPROV

Void TComTrQuant::FwdNsst4x4( Int* src, UInt uiMode, UChar index )
{
  const Int   rnd = NSST_HYGT_RNDS_4x4;
  const Int   shl = 5;
  const Int * par = g_nsstHyGTPar4x4[uiMode][index];
  const Int   cof = 1 << (shl + 9);

  assert( index<4 );

  for (Int k = 0; k < 16; k++) src[k] <<= shl;

  for (Int r = 0, q = (4 * rnd - 1); r < rnd; r++) 
  {
    for (Int d = 0; d < 4; d++, q--) 
    {
      const Int   s = 1 << d;
      const Int * p = par + ((r * 4 + d) * 8);
      if (q > 0 )
      {
        for (Int i = 0; i < 8; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a - t->s * b + 512) >> 10;
          src[j + s] = (t->c * b + t->s * a + 512) >> 10;
        }
      }
      else
      {
        for (Int i = 0; i < 8; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a - t->s * b + cof) >> (10 + shl);
          src[j + s] = (t->c * b + t->s * a + cof) >> (10 + shl);
        }
      }
    }
  }
}

Void TComTrQuant::InvNsst4x4( Int* src, UInt uiMode, UChar index )
{
  const Int    rnd = NSST_HYGT_RNDS_4x4;
  const Int    shl = 5;
  const Int  * par = g_nsstHyGTPar4x4[uiMode][index];
  const Int    cof = 1 << (shl + 9);

  assert( index<4 );

  for (Int k = 0; k < 16; k++) src[k] <<= shl;

  for (Int r = rnd, q = (4 * rnd - 1); --r >= 0;) 
  {
    for (Int d = 4; --d >= 0; q--) 
    {
      const Int   s = 1 << d;
      const Int * p = par + ((r * 4 + d) * 8);
      if (q > 0 )
      {
        for (Int i = 0; i < 8; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a + t->s * b + 512) >> 10;
          src[j + s] = (t->c * b - t->s * a + 512) >> 10;
        }
      }
      else
      {
        for (Int i = 0; i < 8; i++) {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a + t->s * b + cof) >> (10 + shl);
          src[j + s] = (t->c * b - t->s * a + cof) >> (10 + shl);
        }
      }
    }
  }
}

Void TComTrQuant::FwdNsst8x8( Int* src, UInt uiMode, UChar index )
{
  const Int   rnd = NSST_HYGT_RNDS_8x8;
  const Int   shl = 5;
  const Int * par = g_nsstHyGTPar8x8[uiMode][index];
  const Int   cof = 1 << (shl + 9);

  assert( index<4 );

  for (Int k = 0; k < 64; k++) src[k] <<= shl;

  for (Int r = 0, q = (6 * rnd - 1); r < rnd; r++) 
  {
    for (Int d = 0; d < 6; d++, q--) 
    {
      const Int   s = 1 << d;
      const Int * p = par + ((r * 6 + d) * 32);
      if (q > 0)
      {
        for (Int i = 0; i < 32; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a - t->s * b + 512) >> 10;
          src[j + s] = (t->c * b + t->s * a + 512) >> 10;
        }
      }
      else
      {
        for (int i = 0; i < 32; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a - t->s * b + cof) >> (10 + shl);
          src[j + s] = (t->c * b + t->s * a + cof) >> (10 + shl);
        }
      }
    }
  }
}

Void TComTrQuant::InvNsst8x8( Int* src, UInt uiMode, UChar index )
{
  const Int    rnd = NSST_HYGT_RNDS_8x8;
  const Int    shl = 5;
  const Int  * par = g_nsstHyGTPar8x8[uiMode][index];
  const Int    cof = 1 << (shl + 9);
  
  assert( index<4 );

  for (Int k = 0; k < 64; k++) src[k] <<= shl;

  for (Int r = rnd, q = (6 * rnd - 1); --r >= 0;) 
  {
    for (Int d = 6; --d >= 0; q--) 
    {
      const Int   s = 1 << d;
      const Int * p = par + ((r * 6 + d) * 32);
      if (q > 0)
      {
        for (Int i = 0; i < 32; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a + t->s * b + 512) >> 10;
          src[j + s] = (t->c * b - t->s * a + 512) >> 10;
        }
      }
      else
      {
        for (Int i = 0; i < 32; i++) 
        {
          const tabSinCos * t = &g_tabSinCos[*p++];
          register Int j = i + (i & -s);
          register Int a = src[j];
          register Int b = src[j + s];
          src[j]     = (t->c * a + t->s * b + cof) >> (10 + shl);
          src[j + s] = (t->c * b - t->s * a + cof) >> (10 + shl);
        }
      }
    }
  }
}

#else

Void TComTrQuant::FwdNsst4x4( Int* src, UInt uiMode, UChar index )
{
  const Int *iT = g_aiNsst4x4[uiMode][index][0];
  Int coef, temp[16];

  assert( index<4 );

  for (Int j=0; j<16; j++)
  {
    coef = src[ 0]*iT[ 0] + src[ 1]*iT[ 1] + src[ 2]*iT[ 2] + src[ 3]*iT[ 3] +
           src[ 4]*iT[ 4] + src[ 5]*iT[ 5] + src[ 6]*iT[ 6] + src[ 7]*iT[ 7] +
           src[ 8]*iT[ 8] + src[ 9]*iT[ 9] + src[10]*iT[10] + src[11]*iT[11] +
           src[12]*iT[12] + src[13]*iT[13] + src[14]*iT[14] + src[15]*iT[15];
    temp[j] = ( coef + 128 ) >> 8;
    iT  += 16;
  }

  memcpy( src, temp, 16*sizeof(Int) );
}

Void TComTrQuant::InvNsst4x4( Int* src, UInt uiMode, UChar index )
{
  const Int *iT = g_aiNsst4x4[uiMode][index][0];
  Int  temp[16], resi;

  assert( index<4 );

  for (Int j=0; j<16; j++)
  {
    resi = src[ 0]*iT[ 0*16] + src[ 1]*iT[ 1*16] + src[ 2]*iT[ 2*16] + src[ 3]*iT[ 3*16] +
           src[ 4]*iT[ 4*16] + src[ 5]*iT[ 5*16] + src[ 6]*iT[ 6*16] + src[ 7]*iT[ 7*16] +
           src[ 8]*iT[ 8*16] + src[ 9]*iT[ 9*16] + src[10]*iT[10*16] + src[11]*iT[11*16] +
           src[12]*iT[12*16] + src[13]*iT[13*16] + src[14]*iT[14*16] + src[15]*iT[15*16];
    temp[j] = ( resi + 128 ) >> 8;
    iT  ++;
  }

  memcpy( src, temp, 16*sizeof(Int) );
}

#endif
#endif

#if COM16_C806_EMT
UChar TComTrQuant::getEmtTrIdx( TComTU &rTu, const ComponentID compID )
{
  TComDataCU* pcCU          = rTu.getCU();
  const UInt uiAbsPartIdx   = rTu.GetAbsPartIdxTU();
  UChar ucTrIdx             = DCT2_HEVC;

  if( compID == COMPONENT_Y )
  {
    if (  pcCU->isIntra(uiAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
    {
      ucTrIdx = pcCU->getEmtCuFlag( uiAbsPartIdx ) ? pcCU->getEmtTuIdx( uiAbsPartIdx ) : DCT2_EMT;
    }
    if ( !pcCU->isIntra(uiAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseInterEMT() )
    {
      ucTrIdx = pcCU->getEmtCuFlag( uiAbsPartIdx ) ? pcCU->getEmtTuIdx( uiAbsPartIdx ) : DCT2_EMT;
    }
  }
  else
  {
    if (  pcCU->isIntra(uiAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseIntraEMT() )
    {
      ucTrIdx = DCT2_EMT;
    }
    if ( !pcCU->isIntra(uiAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseInterEMT() )
    {
      ucTrIdx = DCT2_EMT;
    }
  }

  return ucTrIdx;
}

UChar TComTrQuant::getEmtMode( TComTU &rTu, const ComponentID compID )
{
  TComDataCU* pcCU          = rTu.getCU();
  const UInt uiAbsPartIdx   = rTu.GetAbsPartIdxTU();
  UChar ucMode              = 0;

  if( pcCU->isIntra( uiAbsPartIdx ) )
  {
    ucMode = pcCU->getIntraDir( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
  }
  else
  {
    ucMode = INTER_MODE_IDX;
  }
  
  return ucMode;
}
#endif

Void TComTrQuant::transformNxN(       TComTU        & rTu,
                                const ComponentID     compID,
                                      Pel          *  pcResidual,
                                const UInt            uiStride,
                                      TCoeff       *  rpcCoeff,
#if ADAPTIVE_QP_SELECTION
                                      TCoeff       *  pcArlCoeff,
#endif
                                      TCoeff        & uiAbsSum,
                                const QpParam       & cQP 
#if VCEG_AZ08_KLT_COMMON
                                , Bool useKLT
#endif
                                )
{
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt uiWidth        = rect.width;
  const UInt uiHeight       = rect.height;
  TComDataCU* pcCU          = rTu.getCU();
  const UInt uiAbsPartIdx   = rTu.GetAbsPartIdxTU();
  const UInt uiOrgTrDepth   = rTu.GetTransformDepthRel();

  uiAbsSum=0;

  RDPCMMode rdpcmMode = RDPCM_OFF;
  rdpcmNxN( rTu, compID, pcResidual, uiStride, cQP, rpcCoeff, uiAbsSum, rdpcmMode );

  if (rdpcmMode == RDPCM_OFF)
  {
    uiAbsSum = 0;
    //transform and quantise
    if(pcCU->getCUTransquantBypass(uiAbsPartIdx))
    {
      const Bool rotateResidual = rTu.isNonTransformedResidualRotated(compID);
      const UInt uiSizeMinus1   = (uiWidth * uiHeight) - 1;

      for (UInt y = 0, coefficientIndex = 0; y<uiHeight; y++)
      {
        for (UInt x = 0; x<uiWidth; x++, coefficientIndex++)
        {
          const Pel currentSample = pcResidual[(y * uiStride) + x];

          rpcCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] = currentSample;
          uiAbsSum += TCoeff(abs(currentSample));
        }
      }
    }
    else
    {
#if DEBUG_TRANSFORM_AND_QUANTISE
      std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU at input to transform\n";
      printBlock(pcResidual, uiWidth, uiHeight, uiStride);
#endif
#if JVET_C0024_QTBT
      assert( (pcCU->getSlice()->getSPS()->getCTUSize() >= uiWidth) );
#else
      assert( (pcCU->getSlice()->getSPS()->getMaxTrSize() >= uiWidth) );
#endif

      if(pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0)
      {
        xTransformSkip( pcResidual, uiStride, m_plTempCoeff, rTu, compID );
      }
      else
      {
        const Int channelBitDepth=pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
        xT( channelBitDepth, rTu.useDST(compID), pcResidual, uiStride, m_plTempCoeff, uiWidth, uiHeight, pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID)) 
#if COM16_C806_EMT
          , getEmtMode ( rTu, compID )
          , getEmtTrIdx( rTu, compID )
#endif
#if VCEG_AZ08_KLT_COMMON
          , useKLT && (compID == 0)
#endif
          );
      }

#if DEBUG_TRANSFORM_AND_QUANTISE
      std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU between transform and quantiser\n";
      printBlock(m_plTempCoeff, uiWidth, uiHeight, uiWidth);
#endif
#if VCEG_AZ05_ROT_TR
#if JVET_C0024_QTBT
  Char ucROTIdx = pcCU->getROTIdx(toChannelType(compID), uiAbsPartIdx) ;
  if (ucROTIdx && uiWidth>=4 && uiHeight>=4)
#else
  if (pcCU->getROTIdx(uiAbsPartIdx) )
#endif
  {           
            static Int ROT_MATRIX[16];
      Int iSubGroupXMax = Clip3 (1,16,(Int)( (uiWidth>>2)));
      Int iSubGroupYMax = Clip3 (1,16,(Int)( (uiHeight>>2)));

      Int iOffSetX = 0;
      Int iOffSetY = 0;
      Int y = 0;
      Int* piCoeffTemp = m_plTempCoeff;
      Int* piROTTemp = ROT_MATRIX;
      Int iString2CopySize = 4*sizeof(Int);
    for (Int iSubGroupX = 0; iSubGroupX<iSubGroupXMax; iSubGroupX++)
      for (Int iSubGroupY = 0; iSubGroupY<iSubGroupYMax; iSubGroupY++)
      {
        iOffSetX = 4*iSubGroupX;
        iOffSetY = 4*iSubGroupY*uiWidth;
          piROTTemp = ROT_MATRIX;
          piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;
          for(  y = 0; y < 4; y++ )
           {  
             ::memcpy(piROTTemp, piCoeffTemp, iString2CopySize); 
                       piROTTemp +=4;
             piCoeffTemp +=uiWidth;
            }
#if JVET_C0024_QTBT
           RotTransform4I( ROT_MATRIX, ucROTIdx-1 );
#else
           RotTransform4I( ROT_MATRIX, pcCU->getROTIdx(uiAbsPartIdx)-1 );
#endif
          piROTTemp = ROT_MATRIX;
          piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;
          for(  y = 0; y < 4; y++ )
          {    
             ::memcpy(piCoeffTemp,piROTTemp, iString2CopySize);
                       piROTTemp +=4;
             piCoeffTemp +=uiWidth;
          }
             }
  }
#elif COM16_C1044_NSST
#if JVET_C0024_QTBT
    Char ucNsstIdx = pcCU->getROTIdx(toChannelType(compID), uiAbsPartIdx) ;
    if (ucNsstIdx && uiWidth>=4 && uiHeight>=4
#if  JVET_C0045_C0053_NO_NSST_FOR_TS
      && !pcCU->getTransformSkip(uiAbsPartIdx, compID)
#endif
      )
#else
      if (pcCU->getROTIdx(uiAbsPartIdx) 
#if  JVET_C0045_C0053_NO_NSST_FOR_TS
        && !pcCU->getTransformSkip(uiAbsPartIdx, compID)
#endif
        )
#endif
      {           
#if JVET_D0120_NSST_IMPROV
        static Int NSST_MATRIX[64];
        const  Int iLog2SbSize   = (uiWidth > 4 && uiHeight > 4) ? 3 : 2;
        const  Int iSbSize       = (uiWidth > 4 && uiHeight > 4) ? 8 : 4;
#if FAST_NSST_v1 //NSST_TU_44
		const  Int uiSize = min(uiWidth, uiHeight);
		const  Int iSubGroupXMax = Clip3(1, 8, (Int)uiSize) >> iLog2SbSize;
		const  Int iSubGroupYMax = iSubGroupXMax;

#else
        const  Int iSubGroupXMax = Clip3(1, 8, (Int)uiWidth ) >> iLog2SbSize;
        const  Int iSubGroupYMax = Clip3(1, 8, (Int)uiHeight) >> iLog2SbSize;
#endif
#else
        static Int NSST_MATRIX[16];
        Int iSubGroupXMax = Clip3 (1,16,(Int)( (uiWidth>>2)));
        Int iSubGroupYMax = Clip3 (1,16,(Int)( (uiHeight>>2)));
#endif

        Int iOffSetX = 0;
        Int iOffSetY = 0;
        Int y = 0;
        Int* piCoeffTemp = m_plTempCoeff;
        Int* piNsstTemp = NSST_MATRIX;
#if JVET_D0120_NSST_IMPROV
        Int iString2CopySize = iSbSize*sizeof(Int);
#else
        Int iString2CopySize = 4*sizeof(Int);
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
#if !JVET_C0024_QTBT
        const UInt uiLog2BlockSize = g_aucConvertToBit[ uiWidth ] + 2;
#endif
#endif
        const UInt uiScanIdx = pcCU->getCoefScanIdx(uiAbsPartIdx, uiWidth, uiHeight, compID);
#if JVET_C0024_QTBT
        const UInt log2BlockWidth  = g_aucConvertToBit[uiWidth]  + MIN_CU_LOG2;
        const UInt log2BlockHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
#else
        const UInt log2BlockWidth  = g_aucConvertToBit[uiWidth]  + 2;
        const UInt log2BlockHeight = g_aucConvertToBit[uiHeight] + 2;
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
        const UInt *scan = (log2BlockWidth>=3 && log2BlockHeight>=3) ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-3] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
        const UInt *scan = uiLog2BlockSize>3 ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-4] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#else
#if JVET_C0024_QTBT
        const UInt *scan = (log2BlockWidth==3 && log2BlockHeight==3) ? g_auiCoefScanFirstCG8x8[uiScanIdx] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
        const UInt *scan = uiLog2BlockSize==3 ? g_auiCoefScanFirstCG8x8[uiScanIdx] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#endif
#else
#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
        const UInt *scan = (log2BlockWidth>=3 && log2BlockHeight>=3) ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-3] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
        const UInt *scan = uiLog2BlockSize>3 ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-4] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#else
        const UInt *scan = g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#endif
        UInt uiIntraMode = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx);
        if( compID != COMPONENT_Y )
        {
#if JVET_E0062_MULTI_DMS && COM16_C806_LMCHROMA
          if( uiIntraMode == LM_CHROMA_IDX )
          {
            uiIntraMode = PLANAR_IDX;
          }
#if JVET_E0077_ENHANCED_LM
          if (IsLMMode(uiIntraMode))
          {
              uiIntraMode = PLANAR_IDX;
          }
#endif
#else
          if( uiIntraMode == DM_CHROMA_IDX )
          {
#if JVET_C0024_QTBT
            if( pcCU->getSlice()->isIntra() )
            {
#if DIMD_INTRA_PRED
              if (pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx))
              {
                uiIntraMode = DC_IDX + 1;
              }
              else
              {
#endif
                uiIntraMode = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
#if DIMD_INTRA_PRED
              }
#endif
            }
            else
            {
#endif
            uiIntraMode = pcCU->getIntraDir( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
#if JVET_C0024_QTBT
          }
#endif
          }
#if COM16_C806_LMCHROMA
          else if( uiIntraMode == LM_CHROMA_IDX )
          {
            uiIntraMode = PLANAR_IDX;
          }
#endif


#if JVET_E0077_ENHANCED_LM
          if (IsLMMode(uiIntraMode))
          {
              uiIntraMode = PLANAR_IDX;
          }
#endif

#endif

        }

        assert(uiIntraMode < NUM_INTRA_MODE - 1);
#if DIMD_INTRA_PRED
        const Int iNsstCandNum = (pcCU->getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx)) ? 4 : (
#if FAST_NSST_v1
		    (uiIntraMode <= DC_IDX) ? 2 : 4
#else
			(uiIntraMode <= DC_IDX) ? 3 : 4
#endif
			);
#else
#if FAST_NSST_v1
		const Int iNsstCandNum = (uiIntraMode <= DC_IDX) ? 2 : 4;
#else
        const Int iNsstCandNum = (uiIntraMode <= DC_IDX) ? 3 : 4;
#endif
#endif

#if JVET_C0024_QTBT
        if( iNsstCandNum > ucNsstIdx )
#else
        if( iNsstCandNum > pcCU->getROTIdx(uiAbsPartIdx) )
#endif
        {
#if JVET_D0120_NSST_IMPROV 
#if JVET_C0024_QTBT
          const Int * permut = iSbSize>4 ? g_nsstHyGTPermut8x8[g_NsstLut[uiIntraMode]][ucNsstIdx - 1] : g_nsstHyGTPermut4x4[g_NsstLut[uiIntraMode]][ucNsstIdx - 1];
#else
          const Int * permut = iSbSize>4 ? g_nsstHyGTPermut8x8[g_NsstLut[uiIntraMode]][pcCU->getROTIdx(uiAbsPartIdx) - 1] : g_nsstHyGTPermut4x4[g_NsstLut[uiIntraMode]][pcCU->getROTIdx(uiAbsPartIdx) - 1];
#endif
#endif
          for (Int iSubGroupX = 0; iSubGroupX<iSubGroupXMax; iSubGroupX++)
          {
            for (Int iSubGroupY = 0; iSubGroupY<iSubGroupYMax; iSubGroupY++)
            {
#if JVET_D0120_NSST_IMPROV
              iOffSetX = iSbSize*iSubGroupX;
              iOffSetY = iSbSize*iSubGroupY*uiWidth;
#else
              iOffSetX = 4*iSubGroupX;
              iOffSetY = 4*iSubGroupY*uiWidth;
#endif
              piNsstTemp = NSST_MATRIX;
              piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;

#if JVET_D0120_NSST_IMPROV
              for(  y = 0; y < iSbSize; y++ )
#else
              for(  y = 0; y < 4; y++ )
#endif
              {  
                if( uiIntraMode>DIA_IDX )
                {
#if JVET_D0120_NSST_IMPROV
                  if( iSbSize==4 )
                  {
                    piNsstTemp[ 0] = piCoeffTemp[0]; piNsstTemp[ 4] = piCoeffTemp[1];
                    piNsstTemp[ 8] = piCoeffTemp[2]; piNsstTemp[12] = piCoeffTemp[3];
                  }
                  else if( iSbSize==8 )
                  {
                    piNsstTemp[ 0] = piCoeffTemp[0]; piNsstTemp[ 8] = piCoeffTemp[1];
                    piNsstTemp[16] = piCoeffTemp[2]; piNsstTemp[24] = piCoeffTemp[3];
                    piNsstTemp[32] = piCoeffTemp[4]; piNsstTemp[40] = piCoeffTemp[5];
                    piNsstTemp[48] = piCoeffTemp[6]; piNsstTemp[56] = piCoeffTemp[7];
                  }
#else
                  piNsstTemp[ 0] = piCoeffTemp[0]; piNsstTemp[ 4] = piCoeffTemp[1];
                  piNsstTemp[ 8] = piCoeffTemp[2]; piNsstTemp[12] = piCoeffTemp[3];
#endif
                  piNsstTemp ++;
                }
                else
                {
                  ::memcpy(piNsstTemp, piCoeffTemp, iString2CopySize); 
#if JVET_D0120_NSST_IMPROV
                  piNsstTemp += iSbSize;
#else
                  piNsstTemp +=4;
#endif
                }
                piCoeffTemp +=uiWidth;
              }

#if JVET_C0024_QTBT
#if JVET_D0120_NSST_IMPROV
              if ( iSbSize>4 )
              {
                FwdNsst8x8( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
              }
              else
              {
                FwdNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
              }
#else
              FwdNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
#endif
#else
#if JVET_D0120_NSST_IMPROV
              if ( iSbSize>4 )
              {
                FwdNsst8x8( NSST_MATRIX, g_NsstLut[uiIntraMode], pcCU->getROTIdx(uiAbsPartIdx)-1 );
              }
              else 
              {
                FwdNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], pcCU->getROTIdx(uiAbsPartIdx)-1 );
              }
#else
              FwdNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], pcCU->getROTIdx(uiAbsPartIdx)-1 );
#endif
#endif

              piNsstTemp = NSST_MATRIX;
              piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;

#if JVET_D0120_NSST_IMPROV
              for(  y = 0; y < iSbSize*iSbSize; y++ )
#else
              for(  y = 0; y < 16; y++ )
#endif
              {
#if JVET_D0120_NSST_IMPROV 
                piCoeffTemp[scan[y]] = piNsstTemp[permut[y]];
#else
                piCoeffTemp[scan[y]] = piNsstTemp[y];
#endif
              }
            }
          }
        }
      }
#endif
#if VCEG_AZ08_KLT_COMMON
      if (useKLT && (compID == 0))
      {
          TUEntropyCodingParameters codingParameters;
          getTUEntropyCodingParameters(codingParameters, rTu, compID);
          const UInt *scan = codingParameters.scan; 
          reOrderCoeff(m_plTempCoeff, scan, uiWidth, uiHeight);
      }
#endif 
      xQuant( rTu, m_plTempCoeff, rpcCoeff,

#if ADAPTIVE_QP_SELECTION
              pcArlCoeff,
#endif
              uiAbsSum, compID, cQP 
            );

#if DEBUG_TRANSFORM_AND_QUANTISE
      std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU at output of quantiser\n";
      printBlock(rpcCoeff, uiWidth, uiHeight, uiWidth);
#endif
    }
  }

    //set the CBF
#if JVET_C0024_QTBT
  assert(uiOrgTrDepth==0);
  pcCU->setCbfPartRange((uiAbsSum > 0) ? 1 : 0, compID, uiAbsPartIdx, 0);
#else
  pcCU->setCbfPartRange((((uiAbsSum > 0) ? 1 : 0) << uiOrgTrDepth), compID, uiAbsPartIdx, rTu.GetAbsPartIdxNumParts(compID));
#endif
}


Void TComTrQuant::invTransformNxN(      TComTU        &rTu,
                                  const ComponentID    compID,
                                        Pel          *pcResidual,
                                  const UInt           uiStride,
                                        TCoeff       * pcCoeff,
                                  const QpParam       &cQP
#if VCEG_AZ08_KLT_COMMON
                                  , Bool useKLT
#endif
                                        DEBUG_STRING_FN_DECLAREP(psDebug)
                                        )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const TComRectangle &rect = rTu.getRect(compID);
  const UInt uiWidth = rect.width;
  const UInt uiHeight = rect.height;

#if ITH_C2_ZR_ZC_SL
  countZeroRowsZeroCols(uiWidth, uiHeight, pcCoeff, uiWidth, &pcCU->getTUSkipWidth(compID, uiAbsPartIdx), &pcCU->getTUSkipHeight(compID, uiAbsPartIdx));
#endif

#if JVET_C0024_ITSKIP
  UInt skipWidth = pcCU->getTUSkipWidth(compID, uiAbsPartIdx);
  UInt skipHeight = pcCU->getTUSkipHeight(compID, uiAbsPartIdx);
#endif
#if !JVET_C0024_QTBT
  if (uiWidth != uiHeight) //for intra, the TU will have been split above this level, so this condition won't be true, hence this only affects inter
  {
    //------------------------------------------------

    //recurse deeper

    TComTURecurse subTURecurse(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

    do
    {
      //------------------

      const UInt lineOffset = subTURecurse.GetSectionNumber() * subTURecurse.getRect(compID).height;

      Pel    *subTUResidual     = pcResidual + (lineOffset * uiStride);
      TCoeff *subTUCoefficients = pcCoeff     + (lineOffset * subTURecurse.getRect(compID).width);

      invTransformNxN(subTURecurse, compID, subTUResidual, uiStride, subTUCoefficients, cQP
#if VCEG_AZ08_KLT_COMMON
      , useKLT
#endif
      DEBUG_STRING_PASS_INTO(psDebug)
          );

      //------------------

    } while (subTURecurse.nextSection(rTu));

    //------------------------------------------------

    return;
  }
#endif

#if DEBUG_STRING
  if (psDebug)
  {
    std::stringstream ss(stringstream::out);
    printBlockToStream(ss, (compID==0)?"###InvTran ip Ch0: " : ((compID==1)?"###InvTran ip Ch1: ":"###InvTran ip Ch2: "), pcCoeff, uiWidth, uiHeight, uiWidth);
    DEBUG_STRING_APPEND((*psDebug), ss.str())
  }
#endif

  if(pcCU->getCUTransquantBypass(uiAbsPartIdx))
  {
    const Bool rotateResidual = rTu.isNonTransformedResidualRotated(compID);
    const UInt uiSizeMinus1   = (uiWidth * uiHeight) - 1;

    for (UInt y = 0, coefficientIndex = 0; y<uiHeight; y++)
    {
      for (UInt x = 0; x<uiWidth; x++, coefficientIndex++)
      {
        pcResidual[(y * uiStride) + x] = Pel(pcCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex]);
      }
    }
  }
  else
  {
#if DEBUG_TRANSFORM_AND_QUANTISE
    std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU at input to dequantiser\n";
    printBlock(pcCoeff, uiWidth, uiHeight, uiWidth);
#endif
    xDeQuant(rTu, pcCoeff, m_plTempCoeff, compID, cQP);
#if VCEG_AZ05_ROT_TR
#if JVET_C0024_QTBT
  Char ucROTIdx = pcCU->getROTIdx(toChannelType(compID), uiAbsPartIdx) ;
  if (ucROTIdx && uiWidth>=4 && uiHeight>=4)
#else
  if (pcCU->getROTIdx(uiAbsPartIdx))
#endif
  {    
#if !JVET_C0024_QTBT
    Char ucROTIdx = pcCU->getROTIdx(uiAbsPartIdx) ;
#endif
       static Int ROT_MATRIX[16];
      Int iSubGroupXMax = Clip3 (1,16,(Int)( (uiWidth>>2)));
      Int iSubGroupYMax = Clip3 (1,16,(Int)( (uiHeight>>2)));
      Int iOffSetX = 0;
      Int iOffSetY = 0;
      Int y = 0;
      Int* piCoeffTemp = m_plTempCoeff;
      Int* piROTTemp = ROT_MATRIX;
      Int iString2CopySize = 4*sizeof(Int);
      for (Int iSubGroupX = 0; iSubGroupX<iSubGroupXMax; iSubGroupX++)
      for (Int iSubGroupY = 0; iSubGroupY<iSubGroupYMax; iSubGroupY++)
      {
        iOffSetX = 4*iSubGroupX;
        iOffSetY = 4*iSubGroupY*uiWidth;
          piROTTemp = ROT_MATRIX;
          piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;
          for(  y = 0; y < 4; y++ )
          {    
             ::memcpy(piROTTemp, piCoeffTemp, iString2CopySize); 
                       piROTTemp +=4;
             piCoeffTemp +=uiWidth;            
          }
           InvRotTransform4I( ROT_MATRIX, ucROTIdx-1);
          piROTTemp = ROT_MATRIX;
          piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;
          for(  y = 0; y < 4; y++ )
          {    
             ::memcpy( piCoeffTemp, piROTTemp, iString2CopySize); 
                       piROTTemp +=4;
             piCoeffTemp +=uiWidth;          
          }
        }
  }
#elif COM16_C1044_NSST
#if JVET_C0024_QTBT
    Char ucNsstIdx = pcCU->getROTIdx(toChannelType(compID), uiAbsPartIdx) ;
    if (ucNsstIdx && uiWidth>=4 && uiHeight>=4
#if  JVET_C0045_C0053_NO_NSST_FOR_TS
      && !pcCU->getTransformSkip(uiAbsPartIdx, compID)
#endif
      )
#else
    if (pcCU->getROTIdx(uiAbsPartIdx)
#if  JVET_C0045_C0053_NO_NSST_FOR_TS
      && !pcCU->getTransformSkip(uiAbsPartIdx, compID)
#endif
      )
#endif
    {   
#if !JVET_C0024_QTBT
      Char ucNsstIdx = pcCU->getROTIdx(uiAbsPartIdx) ;
#endif
#if JVET_D0120_NSST_IMPROV
      static Int NSST_MATRIX[64];
      const  Int iLog2SbSize   = (uiWidth > 4 && uiHeight > 4) ? 3 : 2;
      const  Int iSbSize       = (uiWidth > 4 && uiHeight > 4) ? 8 : 4;
#if FAST_NSST_v1  //NSST_TU_44
	  const  Int uiSize = min(uiWidth, uiHeight);
	  const  Int iSubGroupXMax = Clip3(1, 8, (Int)uiSize) >> iLog2SbSize;
	  const  Int iSubGroupYMax = iSubGroupXMax;

#else
      const  Int iSubGroupXMax = Clip3(1, 8, (Int)uiWidth ) >> iLog2SbSize;
      const  Int iSubGroupYMax = Clip3(1, 8, (Int)uiHeight) >> iLog2SbSize;
#endif
#else
      static Int NSST_MATRIX[16];
      Int iSubGroupXMax = Clip3 (1,16,(Int)( (uiWidth>>2)));
      Int iSubGroupYMax = Clip3 (1,16,(Int)( (uiHeight>>2)));
#endif

      Int iOffSetX = 0;
      Int iOffSetY = 0;
      Int y = 0;
      Int* piCoeffTemp = m_plTempCoeff;
      Int* piNsstTemp = NSST_MATRIX;
#if JVET_D0120_NSST_IMPROV
      Int iString2CopySize = iSbSize*sizeof(Int);
#else
      Int iString2CopySize = 4*sizeof(Int);
#endif

#if VCEG_AZ07_CTX_RESIDUALCODING
#if !JVET_C0024_QTBT
      const UInt uiLog2BlockSize = g_aucConvertToBit[ uiWidth ] + 2;
#endif
#endif
      const UInt uiScanIdx = pcCU->getCoefScanIdx(uiAbsPartIdx, uiWidth, uiHeight, compID);
#if JVET_C0024_QTBT
      const UInt log2BlockWidth  = g_aucConvertToBit[uiWidth]  + MIN_CU_LOG2;
      const UInt log2BlockHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
#else
      const UInt log2BlockWidth  = g_aucConvertToBit[uiWidth]  + 2;
      const UInt log2BlockHeight = g_aucConvertToBit[uiHeight] + 2;
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
      const UInt *scan = (log2BlockWidth>=3 && log2BlockHeight>=3) ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-3] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
      const UInt *scan = uiLog2BlockSize>3 ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-4] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#else
#if JVET_C0024_QTBT
      const UInt *scan = (log2BlockWidth==3 && log2BlockHeight==3)? g_auiCoefScanFirstCG8x8[uiScanIdx] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
      const UInt *scan = uiLog2BlockSize==3 ? g_auiCoefScanFirstCG8x8[uiScanIdx] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#endif
#else
#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
      const UInt *scan = (log2BlockWidth>=3 && log2BlockHeight>=3) ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-3] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#else
      const UInt *scan = uiLog2BlockSize>3 ? g_auiCoefTopLeftDiagScan8x8[log2BlockWidth-4] : g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#else
      const UInt *scan = g_scanOrder[ SCAN_GROUPED_4x4 ][ uiScanIdx ][ log2BlockWidth    ][ log2BlockHeight    ];
#endif
#endif

      UInt uiIntraMode = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx);
      if( compID != COMPONENT_Y )
      {
#if JVET_E0062_MULTI_DMS && COM16_C806_LMCHROMA
        if( uiIntraMode == LM_CHROMA_IDX )
        {
          uiIntraMode = PLANAR_IDX;
        }
#if JVET_E0077_ENHANCED_LM
        if (IsLMMode(uiIntraMode))
        {
            uiIntraMode = PLANAR_IDX;
        }
#endif
#else
        if( uiIntraMode == DM_CHROMA_IDX )
        {
#if JVET_C0024_QTBT
          if( pcCU->getSlice()->isIntra() )
          {
#if DIMD_INTRA_PRED
            if (pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx))
            {
              uiIntraMode = DC_IDX + 1;
            }
            else
            {
#endif
              uiIntraMode = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
#if DIMD_INTRA_PRED
            }
#endif
          }
          else
          {
#endif
          uiIntraMode = pcCU->getIntraDir( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
#if JVET_C0024_QTBT
        }
#endif
        }
#if COM16_C806_LMCHROMA
        else if( uiIntraMode == LM_CHROMA_IDX )
        {
          uiIntraMode = PLANAR_IDX;
        }
#endif


#if JVET_E0077_ENHANCED_LM
        if (IsLMMode(uiIntraMode))
        {
            uiIntraMode = PLANAR_IDX;
        }
#endif

#endif

      }

      assert( uiIntraMode<NUM_INTRA_MODE-1 );

#if DIMD_INTRA_PRED
      const Int iNsstCandNum = (pcCU->getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx)) ? 4 : (
#if FAST_NSST_v1
	      (uiIntraMode <= DC_IDX) ? 2 : 4
#else
		  (uiIntraMode <= DC_IDX) ? 3 : 4
#endif
		  );
#else
#if FAST_NSST_v1
	  const Int iNsstCandNum = (uiIntraMode <= DC_IDX) ? 2 : 4;
#else
      const Int iNsstCandNum = (uiIntraMode <= DC_IDX) ? 3 : 4;
#endif
#endif

      if( iNsstCandNum > ucNsstIdx )
      {
#if JVET_D0120_NSST_IMPROV 
        const Int * permut = iSbSize>4 ? g_nsstHyGTPermut8x8[g_NsstLut[uiIntraMode]][ucNsstIdx - 1] : g_nsstHyGTPermut4x4[g_NsstLut[uiIntraMode]][ucNsstIdx - 1];
#endif
        for (Int iSubGroupX = 0; iSubGroupX<iSubGroupXMax; iSubGroupX++)
        {
          for (Int iSubGroupY = 0; iSubGroupY<iSubGroupYMax; iSubGroupY++)
          {
#if JVET_D0120_NSST_IMPROV
            iOffSetX = iSbSize*iSubGroupX;
            iOffSetY = iSbSize*iSubGroupY*uiWidth;
#else
            iOffSetX = 4*iSubGroupX;
            iOffSetY = 4*iSubGroupY*uiWidth;
#endif
            piNsstTemp = NSST_MATRIX;
            piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;

#if JVET_D0120_NSST_IMPROV
            for(  y = 0; y < iSbSize*iSbSize; y++ )
#else
            for(  y = 0; y < 16; y++ )
#endif
            {    
#if JVET_D0120_NSST_IMPROV
              piNsstTemp[permut[y]] = piCoeffTemp[scan[y]];
#else
              piNsstTemp[y] = piCoeffTemp[scan[y]];
#endif
            }

#if JVET_D0120_NSST_IMPROV
            if ( iSbSize>4 )
            {
              InvNsst8x8( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
            }
            else
            {
              InvNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
            }
#else
            InvNsst4x4( NSST_MATRIX, g_NsstLut[uiIntraMode], ucNsstIdx-1 );
#endif

            piNsstTemp = NSST_MATRIX;
            piCoeffTemp = m_plTempCoeff+iOffSetX+iOffSetY;
#if JVET_D0120_NSST_IMPROV
            for(  y = 0; y < iSbSize; y++ )
#else
            for(  y = 0; y < 4; y++ )
#endif
            {    
              if( uiIntraMode>DIA_IDX )
              {
#if JVET_D0120_NSST_IMPROV
                if ( iSbSize==4 )
                {
                  piCoeffTemp[0] = piNsstTemp[ 0]; piCoeffTemp[1] = piNsstTemp[ 4];
                  piCoeffTemp[2] = piNsstTemp[ 8]; piCoeffTemp[3] = piNsstTemp[12];
                }
                else if ( iSbSize==8 )
                {
                  piCoeffTemp[0] = piNsstTemp[ 0]; piCoeffTemp[1] = piNsstTemp[ 8];
                  piCoeffTemp[2] = piNsstTemp[16]; piCoeffTemp[3] = piNsstTemp[24];
                  piCoeffTemp[4] = piNsstTemp[32]; piCoeffTemp[5] = piNsstTemp[40];
                  piCoeffTemp[6] = piNsstTemp[48]; piCoeffTemp[7] = piNsstTemp[56];
                }
#else
                piCoeffTemp[0] = piNsstTemp[ 0]; piCoeffTemp[1] = piNsstTemp[ 4];
                piCoeffTemp[2] = piNsstTemp[ 8]; piCoeffTemp[3] = piNsstTemp[12];
#endif
                piNsstTemp ++;
              }
              else
              {
                ::memcpy( piCoeffTemp, piNsstTemp, iString2CopySize); 
#if JVET_D0120_NSST_IMPROV
                piNsstTemp +=iSbSize;
#else
                piNsstTemp +=4;
#endif
              }
              piCoeffTemp +=uiWidth;          
            }
          }
        }
#if JVET_C0024_ITSKIP
#if !ITH_C2_CNT_ZR_ZC
        Int maxX = 0;
        Int maxY = 0;
        piCoeffTemp = m_plTempCoeff;
        for (Int yy=0; yy<uiHeight; yy++)
        {
          for (Int xx=0; xx<uiWidth; xx++)
          {
            if (piCoeffTemp[xx]!=0)
            {
              if (xx > maxX )
              {
                maxX = xx;
      }
              if (yy > maxY)
              {
                maxY = yy;
    }
            }
          }
          piCoeffTemp+=uiWidth;
        }

        skipWidth = uiWidth-1-maxX;
        skipHeight = uiHeight-1-maxY;
#else
		countZeroRowsZeroCols(uiWidth, uiHeight, m_plTempCoeff, uiWidth, &skipWidth, &skipHeight);
#endif
#endif
      }
    }
#endif
#if DEBUG_TRANSFORM_AND_QUANTISE
    std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU between dequantiser and inverse-transform\n";
    printBlock(m_plTempCoeff, uiWidth, uiHeight, uiWidth);
#endif

#if DEBUG_STRING
    if (psDebug)
    {
      std::stringstream ss(stringstream::out);
      printBlockToStream(ss, "###InvTran deq: ", m_plTempCoeff, uiWidth, uiHeight, uiWidth);
      (*psDebug)+=ss.str();
    }
#endif

    if(pcCU->getTransformSkip(uiAbsPartIdx, compID))
    {
      xITransformSkip( m_plTempCoeff, pcResidual, uiStride, rTu, compID );

#if DEBUG_STRING
      if (psDebug)
      {
        std::stringstream ss(stringstream::out);
        printBlockToStream(ss, "###InvTran resi: ", pcResidual, uiWidth, uiHeight, uiStride);
        (*psDebug)+=ss.str();
        (*psDebug)+="(<- was a Transform-skipped block)\n";
      }
#endif
    }
    else
    {
#if O0043_BEST_EFFORT_DECODING
      const Int channelBitDepth = pcCU->getSlice()->getSPS()->getStreamBitDepth(toChannelType(compID));
#else
      const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#endif
      xIT( channelBitDepth, rTu.useDST(compID), m_plTempCoeff, pcResidual, uiStride, uiWidth, uiHeight, pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID)) 
#if COM16_C806_EMT
        , getEmtMode ( rTu, compID )
        , getEmtTrIdx( rTu, compID )
#endif
#if JVET_C0024_ITSKIP
        , skipWidth
        , skipHeight
#endif
#if VCEG_AZ08_KLT_COMMON
        , useKLT
#endif
        );

#if DEBUG_STRING
      if (psDebug)
      {
        std::stringstream ss(stringstream::out);
        printBlockToStream(ss, "###InvTran resi: ", pcResidual, uiWidth, uiHeight, uiStride);
        (*psDebug)+=ss.str();
        (*psDebug)+="(<- was a Transformed block)\n";
      }
#endif
    }

#if DEBUG_TRANSFORM_AND_QUANTISE
    std::cout << g_debugCounter << ": " << uiWidth << "x" << uiHeight << " channel " << compID << " TU at output of inverse-transform\n";
    printBlock(pcResidual, uiWidth, uiHeight, uiStride);
    g_debugCounter++;
#endif
  }

  invRdpcmNxN( rTu, compID, pcResidual, uiStride );
}

Void TComTrQuant::invRecurTransformNxN( const ComponentID compID,
                                        TComYuv *pResidual,
                                        TComTU &rTu
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
                                        , TComYuv *pcPred
#endif
                                        )
{
  if (!rTu.ProcessComponentSection(compID))
  {
    return;
  }

  TComDataCU* pcCU = rTu.getCU();
  UInt absPartIdxTU = rTu.GetAbsPartIdxTU();
  UInt uiTrMode=rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
  assert(uiTrMode==0);
#endif
  if( (pcCU->getCbf(absPartIdxTU, compID, uiTrMode) == 0) && (isLuma(compID) || !pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()) )
  {
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
      if (pcCU->getSlice()->getSPS()->getUseInterKLT())
      {
#endif
          UInt uiDepth = pcCU->getDepth(absPartIdxTU) + uiTrMode;
#if JVET_C0024_QTBT
          UInt uiLog2TrSize = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize() >> uiDepth] + 2;
#else
          UInt uiLog2TrSize = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getMaxCUWidth() >> uiDepth] + 2;
#endif
          Int  tuWidth = rTu.getRect(compID).width;
          Int  tuHeight = rTu.getRect(compID).height;
          if (isChroma(compID) && uiLog2TrSize == 2)
          {
              UInt uiQPDiv = pcCU->getPic()->getNumPartitionsInCtu() >> ((uiDepth - 1) << 1);
              if ((absPartIdxTU % uiQPDiv) != 0)
              {
    return;
  }
              tuWidth <<= 1;
              tuHeight <<= 1;
          }
          //===== reconstruction =====
          UInt    uiZOrder = pcCU->getZorderIdxInCtu() + absPartIdxTU;
          TComPicYuv *picRec = pcCU->getPic()->getPicYuvRec();
          Pel* piReco = picRec->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
          UInt uiRecStride = picRec->getStride(compID);
          Pel* piPred = pcPred->getAddr(compID, absPartIdxTU);
          UInt uiPredStride = pcPred->getStride(compID);

          for (UInt uiY = 0; uiY < tuHeight; uiY++)
          {
              for (UInt uiX = 0; uiX < tuWidth; uiX++)
              {
                  piReco[uiX] = piPred[uiX];
              }
              piPred += uiPredStride;
              piReco += uiRecStride;
          }
#if VCEG_AZ08_USE_KLT
      }
#endif
#endif
    return;
  }

#if !JVET_C0024_QTBT
  if( uiTrMode == pcCU->getTransformIdx( absPartIdxTU ) )
  {
#endif
    const TComRectangle &tuRect      = rTu.getRect(compID);
    const Int            uiStride    = pResidual->getStride( compID );
          Pel           *rpcResidual = pResidual->getAddr( compID );
          UInt           uiAddr      = (tuRect.x0 + uiStride*tuRect.y0);
          Pel           *pResi       = rpcResidual + uiAddr;
          TCoeff        *pcCoeff     = pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID);

    const QpParam cQP(*pcCU, compID);

#if TRACE_FLAGS
#if !JVET_C0024_QTBT
    if (isLuma(compID))
#endif
    {
      DUMP_CU_FLAGS_INFO("decoder_tu", pcCU, absPartIdxTU, pcCU->getZorderIdxInCtu() + absPartIdxTU);
    }
#endif

    if (pcCU->getCbf(absPartIdxTU, compID, uiTrMode) != 0)
    {
      DEBUG_STRING_NEW(sTemp)
#if DEBUG_STRING
      std::string *psDebug=((DebugOptionList::DebugString_InvTran.getInt()&(pcCU->isIntra(absPartIdxTU)?1:(pcCU->isInter(absPartIdxTU)?2:4)))!=0) ? &sTemp : 0;
#endif
#if VCEG_AZ08_INTER_KLT
        Bool useKLT = false;
#if VCEG_AZ08_USE_KLT
        if (pcCU->getSlice()->getSPS()->getUseInterKLT())
        {
#endif
            Int  tuWidth = tuRect.width;
            Int  tuHeight = tuRect.height;
            if (compID == COMPONENT_Y)
            {
                UInt uiMaxTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
                UInt uiMinTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MIN - 1];
                Bool checkKLTY = (tuWidth == tuHeight) && (tuWidth <= uiMaxTrWidth) && (tuWidth >= uiMinTrWidth);
                if (pcCU->getSlice()->getPPS()->getUseTransformSkip())
                {
                    checkKLTY &= (0 == pcCU->getTransformSkip(absPartIdxTU, compID));
                }
                checkKLTY &= (Bool)pcCU->getKLTFlag(absPartIdxTU, compID);
                if (checkKLTY)
                {
                    UInt uiTarDepth = g_aucConvertToBit[tuWidth];
                    UInt uiTempSize = g_uiDepth2InterTempSize[uiTarDepth];
                    getTargetPatch(pcCU, absPartIdxTU, absPartIdxTU, pcPred, tuWidth, uiTempSize);
                    candidateSearch(pcCU, absPartIdxTU, tuWidth, uiTempSize);
                    Bool bKLTAvailable = candidateTrain(tuWidth, uiTempSize);
                    useKLT = bKLTAvailable;
                    if (useKLT)
                    {
                        TUEntropyCodingParameters codingParameters;
                        getTUEntropyCodingParameters(codingParameters, rTu, compID);
                        const UInt *scan = codingParameters.scan;
                        assert(tuWidth == tuHeight);
                        recoverOrderCoeff(pcCoeff, scan, tuWidth, tuHeight);
                    }
                }
            }
            invTransformNxN(rTu, compID, pResi, uiStride, pcCoeff, cQP, useKLT DEBUG_STRING_PASS_INTO(psDebug));
            if (compID == COMPONENT_Y)
            {
#if !JVET_D0033_ADAPTIVE_CLIPPING
                const Int clipbd = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#endif
#if O0043_BEST_EFFORT_DECODING
                const Int bitDepthDelta = pcCU->getSlice()->getSPS()->getStreamBitDepth(toChannelType(compID)) - clipbd;
#endif
                //===== reconstruction =====
                UInt    uiZOrder = pcCU->getZorderIdxInCtu() + absPartIdxTU;
                TComPicYuv *picRec = pcCU->getPic()->getPicYuvRec();
                Pel* piResi = pResi;
                Pel* piReco = picRec->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
                UInt uiRecStride = picRec->getStride(compID);
                Pel* piPred = pcPred->getAddr(compID, absPartIdxTU);
                UInt uiPredStride = pcPred->getStride(compID);
                for (UInt uiY = 0; uiY < tuHeight; uiY++)
                {
                    for (UInt uiX = 0; uiX < tuWidth; uiX++)
                    {
#if O0043_BEST_EFFORT_DECODING
#if JVET_D0033_ADAPTIVE_CLIPPING
                      piReco[uiX] = ClipA(rightShiftEvenRounding<Pel>(piPred[uiX] + piResi[uiX], bitDepthDelta), compID);
#else
                      piReco[uiX] = ClipBD(rightShiftEvenRounding<Pel>(piPred[uiX] + piResi[uiX], bitDepthDelta), clipbd);
#endif
#else
#if JVET_D0033_ADAPTIVE_CLIPPING
                      piReco[uiX] = ClipA(piPred[uiX] + piResi[uiX], compID);
#else
                      piReco[uiX] = ClipBD(piPred[uiX] + piResi[uiX], clipbd);
#endif
#endif

                    }
                    piPred += uiPredStride;
                    piResi += uiStride;
                    piReco += uiRecStride;
                }
            }
#if VCEG_AZ08_USE_KLT
        }
        else
        {
#if PVC_AQS // de-quantization
          setEstQSFactor(pcPred->getAddr(compID) + uiAddr, uiStride, tuRect.width, tuRect.height, pcCU->getSlice()->getSliceType(), pcCU, absPartIdxTU, compID);
#endif
          invTransformNxN(rTu, compID, pResi, uiStride, pcCoeff, cQP, useKLT DEBUG_STRING_PASS_INTO(psDebug));
        }
#endif
#else
#if PVC_AQS // de-quantization
      setEstQSFactor(pcPred->getAddr(compID) + uiAddr, uiStride, tuRect.width, tuRect.height, pcCU->getSlice()->getSliceType(), pcCU, absPartIdxTU, compID);
#endif
      invTransformNxN( rTu, compID, pResi, uiStride, pcCoeff, cQP DEBUG_STRING_PASS_INTO(psDebug) );
#endif
#if DEBUG_STRING
      if (psDebug != 0)
      {
        std::cout << (*psDebug);
      }
#endif
#if JVET_F0096_BILATERAL_FILTER
	  if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
	  {
		  if (isLuma(compID))
		  {
			  UInt minSize = std::min(tuRect.width, tuRect.height);
			  if ((pcCU->getCbf(absPartIdxTU, compID, uiTrMode) != 0) && (pcCU->getQP(absPartIdxTU) > 17) && (minSize < 16))
			  {
				  Pel* piPred = pcPred->getAddr(compID, absPartIdxTU);
				  UInt uiPredStride = pcPred->getStride(compID);
				  Pel* piResi = pResidual->getAddr(compID, absPartIdxTU);
				  UInt uiStrideRes = pResidual->getStride(compID);
				  UInt uiZOrder = pcCU->getZorderIdxInCtu() + absPartIdxTU;
				  TComPicYuv *picRec = pcCU->getPic()->getPicYuvRec();
				  Pel* piReco = picRec->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
				  UInt uiRecStride = picRec->getStride(compID);
				  const Int clipbd = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
				  TComBilateralFilter::instance()->bilateralFilterInter(pcCU, tuRect.width, tuRect.height, piResi, uiStrideRes, piPred, uiPredStride, piReco, uiRecStride, clipbd, pcCU->getQP(absPartIdxTU));
			  }
		  }
	  }
#endif
    }

    if (isChroma(compID) && (pcCU->getCrossComponentPredictionAlpha(absPartIdxTU, compID) != 0))
    {
      const Pel *piResiLuma = pResidual->getAddr( COMPONENT_Y );
      const Int  strideLuma = pResidual->getStride( COMPONENT_Y );
      const Int  tuWidth    = rTu.getRect( compID ).width;
      const Int  tuHeight   = rTu.getRect( compID ).height;

      if(pcCU->getCbf(absPartIdxTU, COMPONENT_Y, uiTrMode) != 0)
      {
        pResi = rpcResidual + uiAddr;
        const Pel *pResiLuma = piResiLuma + uiAddr;

        crossComponentPrediction( rTu, compID, pResiLuma, pResi, pResi, tuWidth, tuHeight, strideLuma, uiStride, uiStride, true );
      }
    }
#if !JVET_C0024_QTBT
  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
        invRecurTransformNxN(compID, pResidual, tuRecurseChild 
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
            , pcPred
#endif
            );
    } while (tuRecurseChild.nextSection(rTu));
  }
#endif
}

Void TComTrQuant::applyForwardRDPCM( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride, const QpParam& cQP, TCoeff* pcCoeff, TCoeff &uiAbsSum, const RDPCMMode mode )
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

  const Bool bLossless      = pcCU->getCUTransquantBypass( uiAbsPartIdx );
  const UInt uiWidth        = rTu.getRect(compID).width;
  const UInt uiHeight       = rTu.getRect(compID).height;
  const Bool rotateResidual = rTu.isNonTransformedResidualRotated(compID);
  const UInt uiSizeMinus1   = (uiWidth * uiHeight) - 1;

  UInt uiX = 0;
  UInt uiY = 0;

        UInt &majorAxis             = (mode == RDPCM_VER) ? uiX      : uiY;
        UInt &minorAxis             = (mode == RDPCM_VER) ? uiY      : uiX;
  const UInt  majorAxisLimit        = (mode == RDPCM_VER) ? uiWidth  : uiHeight;
  const UInt  minorAxisLimit        = (mode == RDPCM_VER) ? uiHeight : uiWidth;

  const Bool bUseHalfRoundingPoint  = (mode != RDPCM_OFF);

  uiAbsSum = 0;

  for ( majorAxis = 0; majorAxis < majorAxisLimit; majorAxis++ )
  {
    TCoeff accumulatorValue = 0; // 32-bit accumulator
    for ( minorAxis = 0; minorAxis < minorAxisLimit; minorAxis++ )
    {
      const UInt sampleIndex      = (uiY * uiWidth) + uiX;
      const UInt coefficientIndex = (rotateResidual ? (uiSizeMinus1-sampleIndex) : sampleIndex);
      const Pel  currentSample    = pcResidual[(uiY * uiStride) + uiX];
      const TCoeff encoderSideDelta = TCoeff(currentSample) - accumulatorValue;

      Pel reconstructedDelta;
      if ( bLossless )
      {
        pcCoeff[coefficientIndex] = encoderSideDelta;
        reconstructedDelta        = (Pel) encoderSideDelta;
      }
      else
      {
        transformSkipQuantOneSample(rTu, compID, encoderSideDelta, pcCoeff, coefficientIndex, cQP, bUseHalfRoundingPoint);
        invTrSkipDeQuantOneSample  (rTu, compID, pcCoeff[coefficientIndex], reconstructedDelta, cQP, coefficientIndex);
      }

      uiAbsSum += abs(pcCoeff[coefficientIndex]);

      if (mode != RDPCM_OFF)
      {
        accumulatorValue += reconstructedDelta;
      }
    }
  }
}

Void TComTrQuant::rdpcmNxN   ( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride, const QpParam& cQP, TCoeff* pcCoeff, TCoeff &uiAbsSum, RDPCMMode& rdpcmMode )
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

  if (!pcCU->isRDPCMEnabled(uiAbsPartIdx) || ((pcCU->getTransformSkip(uiAbsPartIdx, compID) == 0) && !pcCU->getCUTransquantBypass(uiAbsPartIdx)))
  {
    rdpcmMode = RDPCM_OFF;
  }
  else if ( pcCU->isIntra( uiAbsPartIdx ) )
  {
    const ChromaFormat chFmt = pcCU->getPic()->getPicYuvOrg()->getChromaFormat();
    const ChannelType chType = toChannelType(compID);
    const UInt uiChPredMode  = pcCU->getIntraDir( chType, uiAbsPartIdx );
#if JVET_C0024_QTBT
#if JVET_E0062_MULTI_DMS
    const UInt uiChCodedMode = uiChPredMode;
#else
    const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && isChroma(compID)) 
      ? (pcCU->getSlice()->isIntra()? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu()+uiAbsPartIdx)
      :pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx)) : uiChPredMode;
#endif
#else
    const TComSPS *sps=pcCU->getSlice()->getSPS();
    const UInt partsPerMinCU = 1<<(2*(sps->getMaxTotalCUDepth() - sps->getLog2DiffMaxMinCodingBlockSize()));
    const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && isChroma(compID)) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
#endif
    const UInt uiChFinalMode = ((chFmt == CHROMA_422)       && isChroma(compID)) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;
    if (uiChFinalMode == VER_IDX || uiChFinalMode == HOR_IDX)
    {
      rdpcmMode = (uiChFinalMode == VER_IDX) ? RDPCM_VER : RDPCM_HOR;
      applyForwardRDPCM( rTu, compID, pcResidual, uiStride, cQP, pcCoeff, uiAbsSum, rdpcmMode );
    }
    else
    {
      rdpcmMode = RDPCM_OFF;
    }
  }
  else // not intra, need to select the best mode
  {
    const UInt uiWidth  = rTu.getRect(compID).width;
    const UInt uiHeight = rTu.getRect(compID).height;

    RDPCMMode bestMode   = NUMBER_OF_RDPCM_MODES;
    TCoeff    bestAbsSum = std::numeric_limits<TCoeff>::max();
    TCoeff    bestCoefficients[MAX_TU_SIZE * MAX_TU_SIZE];

    for (UInt modeIndex = 0; modeIndex < NUMBER_OF_RDPCM_MODES; modeIndex++)
    {
      const RDPCMMode mode = RDPCMMode(modeIndex);

      TCoeff currAbsSum = 0;

      applyForwardRDPCM( rTu, compID, pcResidual, uiStride, cQP, pcCoeff, currAbsSum, mode );

      if (currAbsSum < bestAbsSum)
      {
        bestMode   = mode;
        bestAbsSum = currAbsSum;
        if (mode != RDPCM_OFF)
        {
          memcpy(bestCoefficients, pcCoeff, (uiWidth * uiHeight * sizeof(TCoeff)));
        }
      }
    }

    rdpcmMode = bestMode;
    uiAbsSum  = bestAbsSum;

    if (rdpcmMode != RDPCM_OFF) //the TU is re-transformed and quantised if DPCM_OFF is returned, so there is no need to preserve it here
    {
      memcpy(pcCoeff, bestCoefficients, (uiWidth * uiHeight * sizeof(TCoeff)));
    }
  }

  pcCU->setExplicitRdpcmModePartRange(rdpcmMode, compID, uiAbsPartIdx, rTu.GetAbsPartIdxNumParts(compID));
}

Void TComTrQuant::invRdpcmNxN( TComTU& rTu, const ComponentID compID, Pel* pcResidual, const UInt uiStride )
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

  if (pcCU->isRDPCMEnabled( uiAbsPartIdx ) && ((pcCU->getTransformSkip(uiAbsPartIdx, compID ) != 0) || pcCU->getCUTransquantBypass(uiAbsPartIdx)))
  {
    const UInt uiWidth  = rTu.getRect(compID).width;
    const UInt uiHeight = rTu.getRect(compID).height;

    RDPCMMode rdpcmMode = RDPCM_OFF;

    if ( pcCU->isIntra( uiAbsPartIdx ) )
    {
      const ChromaFormat chFmt = pcCU->getPic()->getPicYuvRec()->getChromaFormat();
      const ChannelType chType = toChannelType(compID);
      const UInt uiChPredMode  = pcCU->getIntraDir( chType, uiAbsPartIdx );
#if JVET_C0024_QTBT
#if JVET_E0062_MULTI_DMS
      const UInt uiChCodedMode = uiChPredMode;
#else
      const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && isChroma(compID)) 
        ? (pcCU->getSlice()->isIntra()? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu()+uiAbsPartIdx)
        :pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx)) : uiChPredMode;
#endif
#else
      const TComSPS *sps=pcCU->getSlice()->getSPS();
      const UInt partsPerMinCU = 1<<(2*(sps->getMaxTotalCUDepth() - sps->getLog2DiffMaxMinCodingBlockSize()));
      const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && isChroma(compID)) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
#endif
      const UInt uiChFinalMode = ((chFmt == CHROMA_422)       && isChroma(compID)) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;
      if (uiChFinalMode == VER_IDX || uiChFinalMode == HOR_IDX)
      {
        rdpcmMode = (uiChFinalMode == VER_IDX) ? RDPCM_VER : RDPCM_HOR;
      }
    }
    else  // not intra case
    {
      rdpcmMode = RDPCMMode(pcCU->getExplicitRdpcmMode( compID, uiAbsPartIdx ));
    }

    const TCoeff pelMin=(TCoeff) std::numeric_limits<Pel>::min();
    const TCoeff pelMax=(TCoeff) std::numeric_limits<Pel>::max();
    if (rdpcmMode == RDPCM_VER)
    {
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        Pel *pcCurResidual = pcResidual+uiX;
        TCoeff accumulator = *pcCurResidual; // 32-bit accumulator
        pcCurResidual+=uiStride;
        for( UInt uiY = 1; uiY < uiHeight; uiY++, pcCurResidual+=uiStride )
        {
          accumulator += *(pcCurResidual);
          *pcCurResidual = (Pel)Clip3<TCoeff>(pelMin, pelMax, accumulator);
        }
      }
    }
    else if (rdpcmMode == RDPCM_HOR)
    {
      for( UInt uiY = 0; uiY < uiHeight; uiY++ )
      {
        Pel *pcCurResidual = pcResidual+uiY*uiStride;
        TCoeff accumulator = *pcCurResidual;
        pcCurResidual++;
        for( UInt uiX = 1; uiX < uiWidth; uiX++, pcCurResidual++ )
        {
          accumulator += *(pcCurResidual);
          *pcCurResidual = (Pel)Clip3<TCoeff>(pelMin, pelMax, accumulator);
        }
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------
// Logical transform
// ------------------------------------------------------------------------------------------------

/** Wrapper function between HM interface and core NxN forward transform (2D)
 *  \param channelBitDepth bit depth of channel
 *  \param useDST
 *  \param piBlkResi input data (residual)
 *  \param uiStride stride of input residual data
 *  \param psCoeff output data (transform coefficients)
 *  \param iWidth transform width
 *  \param iHeight transform height
 *  \param maxLog2TrDynamicRange
 */
Void TComTrQuant::xT( const Int channelBitDepth, Bool useDST, Pel* piBlkResi, UInt uiStride, TCoeff* psCoeff, Int iWidth, Int iHeight, const Int maxLog2TrDynamicRange 
#if COM16_C806_EMT
                     , UChar ucMode
                     , UChar ucTrIdx
#endif
#if VCEG_AZ08_KLT_COMMON
                     , Bool useKLT
#endif
                     )
{
#if MATRIX_MULT
  if( iWidth == iHeight)
  {
#if COM16_C806_EMT
    if( ucTrIdx!=DCT2_HEVC )
    {
      xTr_EMT(channelBitDepth, piBlkResi, psCoeff, uiStride, (UInt)iWidth, useDST, maxLog2TrDynamicRange, ucMode, ucTrIdx );
    }
    else
#endif
    xTr(channelBitDepth, piBlkResi, psCoeff, uiStride, (UInt)iWidth, useDST, maxLog2TrDynamicRange);
    return;
  }
#endif

  TCoeff block[ MAX_TU_SIZE * MAX_TU_SIZE ];
  TCoeff coeff[ MAX_TU_SIZE * MAX_TU_SIZE ];

  for (Int y = 0; y < iHeight; y++)
  {
    for (Int x = 0; x < iWidth; x++)
    {
      block[(y * iWidth) + x] = piBlkResi[(y * uiStride) + x];
    }
  }
#if COM16_C806_EMT
#if VCEG_AZ08_KLT_COMMON
  if( ucTrIdx!=DCT2_HEVC && useKLT == false)
#else
  if( ucTrIdx!=DCT2_HEVC )
#endif
  {
    xTrMxN_EMT(channelBitDepth, block, coeff, iWidth, iHeight, useDST, maxLog2TrDynamicRange, ucMode, ucTrIdx );
  }
  else
#endif

#if VCEG_AZ08_KLT_COMMON
  xTrMxN( channelBitDepth, block, coeff, iWidth, iHeight, useDST, maxLog2TrDynamicRange, useKLT);
#else
  xTrMxN( channelBitDepth, block, coeff, iWidth, iHeight, useDST, maxLog2TrDynamicRange );
#endif
  memcpy(psCoeff, coeff, (iWidth * iHeight * sizeof(TCoeff)));
}

/** Wrapper function between HM interface and core NxN inverse transform (2D)
 *  \param channelBitDepth bit depth of channel
 *  \param useDST
 *  \param plCoef input data (transform coefficients)
 *  \param pResidual output data (residual)
 *  \param uiStride stride of input residual data
 *  \param iWidth transform width
 *  \param iHeight transform height
 *  \param maxLog2TrDynamicRange
 */
Void TComTrQuant::xIT( const Int channelBitDepth, Bool useDST, TCoeff* plCoef, Pel* pResidual, UInt uiStride, Int iWidth, Int iHeight, const Int maxLog2TrDynamicRange 
#if COM16_C806_EMT
                      , UChar ucMode
                      , UChar ucTrIdx
#endif
#if JVET_C0024_ITSKIP
                      , UInt uiSkipWidth
                      , UInt uiSkipHeight
#endif
#if VCEG_AZ08_KLT_COMMON
                      , Bool useKLT
#endif
                      )
{
#if MATRIX_MULT
  if( iWidth == iHeight )
  {
#if COM16_C806_EMT
    if( ucTrIdx!=DCT2_HEVC )
    {
      xITr_EMT(channelBitDepth, plCoef, pResidual, uiStride, (UInt)iWidth, useDST, maxLog2TrDynamicRange, ucMode, ucTrIdx );
    }
    else
#endif
    xITr(channelBitDepth, plCoef, pResidual, uiStride, (UInt)iWidth, useDST, maxLog2TrDynamicRange);
    return;
  }
#endif

  TCoeff block[ MAX_TU_SIZE * MAX_TU_SIZE ];
  TCoeff coeff[ MAX_TU_SIZE * MAX_TU_SIZE ];

  memcpy(coeff, plCoef, (iWidth * iHeight * sizeof(TCoeff)));

#if COM16_C806_EMT
#if VCEG_AZ08_KLT_COMMON
  if(ucTrIdx!=DCT2_HEVC && useKLT == false)
#else
  if(ucTrIdx!=DCT2_HEVC)
#endif
  {
#if JVET_C0024_ITSKIP
    xITrMxN_EMT( channelBitDepth, coeff, block, iWidth, iHeight, uiSkipWidth, uiSkipHeight, useDST, maxLog2TrDynamicRange, ucMode, ucTrIdx );
#else
    xITrMxN_EMT( channelBitDepth, coeff, block, iWidth, iHeight, useDST, maxLog2TrDynamicRange, ucMode, ucTrIdx );
#endif
  }
  else
#endif
#if JVET_C0024_ITSKIP
  xITrMxN( channelBitDepth, coeff, block, iWidth, iHeight, uiSkipWidth, uiSkipHeight, useDST, maxLog2TrDynamicRange
#if VCEG_AZ08_KLT_COMMON
  , useKLT
#endif
 );
#else
  xITrMxN( channelBitDepth, coeff, block, iWidth, iHeight, useDST, maxLog2TrDynamicRange
#if VCEG_AZ08_KLT_COMMON
  , useKLT
#endif
  );
#endif

  for (Int y = 0; y < iHeight; y++)
  {
    for (Int x = 0; x < iWidth; x++)
    {
      pResidual[(y * uiStride) + x] = Pel(block[(y * iWidth) + x]);
    }
  }
}

/** Wrapper function between HM interface and core 4x4 transform skipping
 *  \param piBlkResi input data (residual)
 *  \param uiStride stride of input residual data
 *  \param psCoeff output data (transform coefficients)
 *  \param rTu reference to transform data
 *  \param component colour component
 */
Void TComTrQuant::xTransformSkip( Pel* piBlkResi, UInt uiStride, TCoeff* psCoeff, TComTU &rTu, const ComponentID component )
{
  const TComRectangle &rect = rTu.getRect(component);
  const Int width           = rect.width;
  const Int height          = rect.height;
  const Int maxLog2TrDynamicRange = rTu.getCU()->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(component));
  const Int channelBitDepth = rTu.getCU()->getSlice()->getSPS()->getBitDepth(toChannelType(component));

  Int iTransformShift = getTransformShift(channelBitDepth, rTu.GetEquivalentLog2TrSize(component), maxLog2TrDynamicRange);
  if (rTu.getCU()->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag())
  {
    iTransformShift = std::max<Int>(0, iTransformShift);
  }

#if JVET_C0024_QTBT
  Int iWHScale = 1;
  if ((g_aucConvertToBit[ rTu.getRect(component).width ] + g_aucConvertToBit[ rTu.getRect(component).height ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iTransformShift -= 8;
    iWHScale = 181;
  }
#endif

  const Bool rotateResidual = rTu.isNonTransformedResidualRotated(component);
  const UInt uiSizeMinus1   = (width * height) - 1;

  if (iTransformShift >= 0)
  {
#if JVET_C0024_QTBT
    assert((g_aucConvertToBit[ rTu.getRect(component).width ] + g_aucConvertToBit[ rTu.getRect(component).height ] + (MIN_CU_LOG2<<1))%2 ==0);
#endif
    for (UInt y = 0, coefficientIndex = 0; y < height; y++)
    {
      for (UInt x = 0; x < width; x++, coefficientIndex++)
      {
#if JVET_C0024_QTBT
        psCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] = (TCoeff(piBlkResi[(y * uiStride) + x]) * iWHScale) << iTransformShift;
#else
        psCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] = TCoeff(piBlkResi[(y * uiStride) + x]) << iTransformShift;
#endif
      }
    }
  }
  else //for very high bit depths
  {
    iTransformShift = -iTransformShift;
    const TCoeff offset = 1 << (iTransformShift - 1);

    for (UInt y = 0, coefficientIndex = 0; y < height; y++)
    {
      for (UInt x = 0; x < width; x++, coefficientIndex++)
      {
#if JVET_C0024_QTBT
        psCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] = (TCoeff(piBlkResi[(y * uiStride) + x]) * iWHScale + offset) >> iTransformShift;
#else
        psCoeff[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] = (TCoeff(piBlkResi[(y * uiStride) + x]) + offset) >> iTransformShift;
#endif
      }
    }
  }
}

/** Wrapper function between HM interface and core NxN transform skipping
 *  \param plCoef input data (coefficients)
 *  \param pResidual output data (residual)
 *  \param uiStride stride of input residual data
 *  \param rTu reference to transform data
 *  \param component colour component ID
 */
Void TComTrQuant::xITransformSkip( TCoeff* plCoef, Pel* pResidual, UInt uiStride, TComTU &rTu, const ComponentID component )
{
  const TComRectangle &rect = rTu.getRect(component);
  const Int width           = rect.width;
  const Int height          = rect.height;
  const Int maxLog2TrDynamicRange = rTu.getCU()->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(component));
#if O0043_BEST_EFFORT_DECODING
  const Int channelBitDepth = rTu.getCU()->getSlice()->getSPS()->getStreamBitDepth(toChannelType(component));
#else
  const Int channelBitDepth = rTu.getCU()->getSlice()->getSPS()->getBitDepth(toChannelType(component));
#endif

  Int iTransformShift = getTransformShift(channelBitDepth, rTu.GetEquivalentLog2TrSize(component), maxLog2TrDynamicRange);
  if (rTu.getCU()->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag())
  {
    iTransformShift = std::max<Int>(0, iTransformShift);
  }
#if JVET_C0024_QTBT
  Int iWHScale = 1;
  if ((g_aucConvertToBit[ rTu.getRect(component).width ] + g_aucConvertToBit[ rTu.getRect(component).height ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iTransformShift += 7;
    iWHScale = 181;
  }
#endif

  const Bool rotateResidual = rTu.isNonTransformedResidualRotated(component);
  const UInt uiSizeMinus1   = (width * height) - 1;

  if (iTransformShift >= 0)
  {
    const TCoeff offset = iTransformShift==0 ? 0 : (1 << (iTransformShift - 1));

    for (UInt y = 0, coefficientIndex = 0; y < height; y++)
    {
      for (UInt x = 0; x < width; x++, coefficientIndex++)
      {
#if JVET_C0024_QTBT
        pResidual[(y * uiStride) + x] =  Pel((plCoef[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] * iWHScale + offset) >> iTransformShift);
#else
        pResidual[(y * uiStride) + x] =  Pel((plCoef[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] + offset) >> iTransformShift);
#endif
      }
    }
  }
  else //for very high bit depths
  {
    iTransformShift = -iTransformShift;

    for (UInt y = 0, coefficientIndex = 0; y < height; y++)
    {
      for (UInt x = 0; x < width; x++, coefficientIndex++)
      {
#if JVET_C0024_QTBT
        pResidual[(y * uiStride) + x] = Pel(plCoef[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] * iWHScale << iTransformShift);
#else
        pResidual[(y * uiStride) + x] = Pel(plCoef[rotateResidual ? (uiSizeMinus1 - coefficientIndex) : coefficientIndex] << iTransformShift);
#endif
      }
    }
  }
}

/** RDOQ with CABAC
 * \param rTu reference to transform data
 * \param plSrcCoeff pointer to input buffer
 * \param piDstCoeff reference to pointer to output buffer
 * \param piArlDstCoeff
 * \param uiAbsSum reference to absolute sum of quantized transform coefficient
 * \param compID colour component ID
 * \param cQP reference to quantization parameters

 * Rate distortion optimized quantization for entropy
 * coding engines using probability models like CABAC
 */
Void TComTrQuant::xRateDistOptQuant                 (       TComTU       &rTu,
                                                            TCoeff      * plSrcCoeff,
                                                            TCoeff      * piDstCoeff,
#if ADAPTIVE_QP_SELECTION
                                                            TCoeff      * piArlDstCoeff,
#endif
                                                            TCoeff       &uiAbsSum,
                                                      const ComponentID   compID,
                                                      const QpParam      &cQP  
                                                    )
{
  const TComRectangle  & rect             = rTu.getRect(compID);
  const UInt             uiWidth          = rect.width;
  const UInt             uiHeight         = rect.height;
        TComDataCU    *  pcCU             = rTu.getCU();
  const UInt             uiAbsPartIdx     = rTu.GetAbsPartIdxTU();
  const ChannelType      channelType      = toChannelType(compID);
  const UInt             uiLog2TrSize     = rTu.GetEquivalentLog2TrSize(compID);
#if JVET_C0024_QTBT
  const UInt             uiLog2BlockWidth    = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
  const UInt             uiLog2BlockHeight   = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
#endif

  const Bool             extendedPrecision = pcCU->getSlice()->getSPS()->getSpsRangeExtension().getExtendedPrecisionProcessingFlag();
  const Int              maxLog2TrDynamicRange = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));
  const Int              channelBitDepth = rTu.getCU()->getSlice()->getSPS()->getBitDepth(channelType);

  /* for 422 chroma blocks, the effective scaling applied during transformation is not a power of 2, hence it cannot be
   * implemented as a bit-shift (the quantised result will be sqrt(2) * larger than required). Alternatively, adjust the
   * uiLog2TrSize applied in iTransformShift, such that the result is 1/sqrt(2) the required result (i.e. smaller)
   * Then a QP+3 (sqrt(2)) or QP-3 (1/sqrt(2)) method could be used to get the required result
   */

  // Represents scaling through forward transform
  Int iTransformShift = getTransformShift(channelBitDepth, uiLog2TrSize, maxLog2TrDynamicRange);
  if ((pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0) && extendedPrecision)
  {
    iTransformShift = std::max<Int>(0, iTransformShift);
  }
#if JVET_C0024_QTBT && !JVET_C0024_QTBT_FIX_QUANT_TICKET25
  Int iWHScale = 1;
  Int iWHShift = 0;
  Int iWHOffset = 0;
  if ((g_aucConvertToBit[ uiWidth ] + g_aucConvertToBit[ uiHeight ] + (MIN_CU_LOG2<<1))%2 !=0)
  {
    iWHOffset = 64; //1<<(iWHShift-1)
    iWHShift = 7;
    iWHScale = 181;
  }
#endif
  const Bool bUseGolombRiceParameterAdaptation = pcCU->getSlice()->getSPS()->getSpsRangeExtension().getPersistentRiceAdaptationEnabledFlag();
#if VCEG_AZ07_CTX_RESIDUALCODING
  UInt uiPrevGRParam = 0;
#endif
  const UInt initialGolombRiceParameter        = m_pcEstBitsSbac->golombRiceAdaptationStatistics[rTu.getGolombRiceStatisticsIndex(compID)] / RExt__GOLOMB_RICE_INCREMENT_DIVISOR;
        UInt uiGoRiceParam                     = initialGolombRiceParameter;

        Double     d64BlockUncodedCost               = 0;
#if !JVET_C0024_QTBT
  const UInt uiLog2BlockWidth                  = g_aucConvertToBit[ uiWidth  ] + 2;
#endif
#if !VCEG_AZ07_CTX_RESIDUALCODING && !JVET_C0024_QTBT
  const UInt uiLog2BlockHeight                 = g_aucConvertToBit[ uiHeight ] + 2;
#endif
  const UInt uiMaxNumCoeff                     = uiWidth * uiHeight;
  assert(compID<MAX_NUM_COMPONENT);

  Int scalingListType = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
  assert(scalingListType < SCALING_LIST_NUM);

#if ADAPTIVE_QP_SELECTION
  memset(piArlDstCoeff, 0, sizeof(TCoeff) *  uiMaxNumCoeff);
#endif

#if JVET_C0024_QTBT //to prevent stack overflow
#if ITH_C1_RDOQ_MEM_OPT
  Double* pdCostCoeff = m_pdCostCoeff;
  Double* pdCostSig = m_pdCostSig;
  Double* pdCostCoeff0 = m_pdCostCoeff0;
#else
  Double* pdCostCoeff = new Double [uiMaxNumCoeff];
  Double* pdCostSig = new Double [uiMaxNumCoeff];
  Double* pdCostCoeff0 = new Double [uiMaxNumCoeff];
#endif
#else
  Double pdCostCoeff [ MAX_TU_SIZE * MAX_TU_SIZE ];
  Double pdCostSig   [ MAX_TU_SIZE * MAX_TU_SIZE ];
  Double pdCostCoeff0[ MAX_TU_SIZE * MAX_TU_SIZE ];
#endif
#if !ITH_C1_RDOQ_MEM_OPT
  memset( pdCostCoeff, 0, sizeof(Double) *  uiMaxNumCoeff );
  memset( pdCostSig,   0, sizeof(Double) *  uiMaxNumCoeff );
#endif
#if JVET_C0024_QTBT
#if ITH_C1_RDOQ_MEM_OPT
  Int* rateIncUp = m_piRateIncUp;
  Int* rateIncDown = m_piRateIncDown;
  Int* sigRateDelta = m_piSigRateDelta;
  TCoeff* deltaU = m_pcDeltaU;
#else
  Int* rateIncUp   = new Int [uiMaxNumCoeff];
  Int* rateIncDown = new Int [uiMaxNumCoeff];
  Int* sigRateDelta= new Int [uiMaxNumCoeff];
  TCoeff* deltaU   = new Int [uiMaxNumCoeff];
#endif
#else
  Int rateIncUp   [ MAX_TU_SIZE * MAX_TU_SIZE ];
  Int rateIncDown [ MAX_TU_SIZE * MAX_TU_SIZE ];
  Int sigRateDelta[ MAX_TU_SIZE * MAX_TU_SIZE ];
  TCoeff deltaU   [ MAX_TU_SIZE * MAX_TU_SIZE ];
#endif
#if !ITH_C1_RDOQ_MEM_OPT
  memset( rateIncUp,    0, sizeof(Int   ) *  uiMaxNumCoeff );
  memset( rateIncDown,  0, sizeof(Int   ) *  uiMaxNumCoeff );
  memset( sigRateDelta, 0, sizeof(Int   ) *  uiMaxNumCoeff );
  memset( deltaU,       0, sizeof(TCoeff) *  uiMaxNumCoeff );
#endif
  const Int iQBits = QUANT_SHIFT + cQP.per + iTransformShift;                   // Right shift of non-RDOQ quantizer;  level = (coeff*uiQ + offset)>>q_bits
#if JVET_C0024_QTBT
  const Double *const pdErrScale = getErrScaleCoeff(scalingListType, (uiLog2BlockWidth-1), (uiLog2BlockHeight-1), cQP.rem);
  const Int    *const piQCoef    = getQuantCoeff(scalingListType, cQP.rem, (uiLog2BlockWidth-1), (uiLog2BlockHeight-1));
#else
  const Double *const pdErrScale = getErrScaleCoeff(scalingListType, (uiLog2TrSize-2), cQP.rem);
  const Int    *const piQCoef    = getQuantCoeff(scalingListType, cQP.rem, (uiLog2TrSize-2));
#endif

  const Bool   enableScalingLists             = getUseScalingList(uiWidth, uiHeight, (pcCU->getTransformSkip(uiAbsPartIdx, compID) != 0));
#if JVET_C0024_QTBT_FIX_QUANT_TICKET25
  const Int    defaultQuantisationCoefficient = ( uiLog2BlockWidth + uiLog2BlockHeight ) & 0x01 ? ( g_quantScales[cQP.rem] * 181 ) >> 7 : g_quantScales[cQP.rem];
#else
  const Int    defaultQuantisationCoefficient = g_quantScales[cQP.rem];
#endif
#if JVET_C0024_QTBT
  const Double defaultErrorScale              = getErrScaleCoeffNoScalingList(scalingListType, (uiLog2BlockWidth-1), (uiLog2BlockHeight-1), cQP.rem);
#else
  const Double defaultErrorScale              = getErrScaleCoeffNoScalingList(scalingListType, (uiLog2TrSize-2), cQP.rem);
#endif

  const TCoeff entropyCodingMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff entropyCodingMaximum =  (1 << maxLog2TrDynamicRange) - 1;

#if ADAPTIVE_QP_SELECTION
  Int iQBitsC = iQBits - ARL_C_PRECISION;
  Int iAddC =  1 << (iQBitsC-1);
#endif

  TUEntropyCodingParameters codingParameters;
  getTUEntropyCodingParameters(codingParameters, rTu, compID);
  const UInt uiCGSize = (1 << MLS_CG_SIZE);
#if VCEG_AZ07_CTX_RESIDUALCODING
  Bool bHor8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_HOR;
  Bool bVer8x8 = uiWidth == 8 && uiHeight == 8 && codingParameters.scanType == SCAN_VER;
  Bool bNonZig8x8 = bHor8x8 || bVer8x8; 
#endif

  Double pdCostCoeffGroupSig[ MLS_GRP_NUM ];
  UInt uiSigCoeffGroupFlag[ MLS_GRP_NUM ];
  Int iCGLastScanPos = -1;
#if VCEG_AZ07_CTX_RESIDUALCODING
#if JVET_C0024_QTBT
  UInt uiOffsetonTU = (compID == COMPONENT_Y) ? (uiLog2TrSize<=2 ? 0: NUM_SIG_FLAG_CTX_LUMA_TU << ( min(1, (Int)(uiLog2TrSize - 3)))) : 0 ; 
#else
  UInt uiOffsetonTU = (compID == COMPONENT_Y) ? (uiLog2BlockWidth==2 ? 0: NUM_SIG_FLAG_CTX_LUMA_TU << ( min(1, (Int)(uiLog2BlockWidth - 3)))) : 0 ; 
#endif
#else
  UInt    uiCtxSet            = 0;
  Int     c1                  = 1;
  Int     c2                  = 0;
#endif
  Double  d64BaseCost         = 0;
  Int     iLastScanPos        = -1;

  UInt    c1Idx     = 0;
  UInt    c2Idx     = 0;
#if !VCEG_AZ07_CTX_RESIDUALCODING
  Int     baseLevel;
#endif

  memset( pdCostCoeffGroupSig,   0, sizeof(Double) * MLS_GRP_NUM );
  memset( uiSigCoeffGroupFlag,   0, sizeof(UInt) * MLS_GRP_NUM );

  UInt uiCGNum = uiWidth * uiHeight >> MLS_CG_SIZE;
  Int iScanPos;
  coeffGroupRDStats rdStats;

#if VCEG_AZ07_CTX_RESIDUALCODING
  const UInt significanceMapContextOffset = 0; 
#else
  const UInt significanceMapContextOffset = getSignificanceMapContextOffset(compID);
#endif

  for (Int iCGScanPos = uiCGNum-1; iCGScanPos >= 0; iCGScanPos--)
  {
    UInt uiCGBlkPos = codingParameters.scanCG[ iCGScanPos ];
    UInt uiCGPosY   = uiCGBlkPos / codingParameters.widthInGroups;
    UInt uiCGPosX   = uiCGBlkPos - (uiCGPosY * codingParameters.widthInGroups);
#if VCEG_AZ07_CTX_RESIDUALCODING
    if( bNonZig8x8 )
    {
      uiCGPosY = (bHor8x8 ? uiCGBlkPos : 0);
      uiCGPosX = (bVer8x8 ? uiCGBlkPos : 0);
    }
#endif
    memset( &rdStats, 0, sizeof (coeffGroupRDStats));
#if !VCEG_AZ07_CTX_RESIDUALCODING
    const Int patternSigCtx = TComTrQuant::calcPatternSigCtx(uiSigCoeffGroupFlag, uiCGPosX, uiCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups);
#endif
    for (Int iScanPosinCG = uiCGSize-1; iScanPosinCG >= 0; iScanPosinCG--)
    {
      iScanPos = iCGScanPos*uiCGSize + iScanPosinCG;
      //===== quantization =====
      UInt    uiBlkPos          = codingParameters.scan[iScanPos];
      // set coeff

      const Int    quantisationCoefficient = (enableScalingLists) ? piQCoef   [uiBlkPos] : defaultQuantisationCoefficient;
      const Double errorScale              = (enableScalingLists) ? pdErrScale[uiBlkPos] : defaultErrorScale;

#if JVET_C0024_QTBT && !JVET_C0024_QTBT_FIX_QUANT_TICKET25
      const Int64  tmpLevel                = (Int64(abs(plSrcCoeff[ uiBlkPos ])) * quantisationCoefficient * iWHScale + iWHOffset)>>iWHShift;
#else
#if PVC_AQS // in RDOQ
      Int64 tmpLevel = Int64(abs(plSrcCoeff[uiBlkPos])) * quantisationCoefficient;
      if (compID == COMPONENT_Y)
      {
        Int iEstQSFactorShift14 = m_iEstQSFactor << (14 - ESM_SHIFT);
        tmpLevel = (tmpLevel * (Int64) iEstQSFactorShift14 + (1 >> 13)) >> 14;
      }
#else
      const Int64  tmpLevel = Int64(abs(plSrcCoeff[uiBlkPos])) * quantisationCoefficient;
#endif
#endif

      const Intermediate_Int lLevelDouble  = (Intermediate_Int)min<Int64>(tmpLevel, std::numeric_limits<Intermediate_Int>::max() - (Intermediate_Int(1) << (iQBits - 1)));

#if ADAPTIVE_QP_SELECTION
      if( m_bUseAdaptQpSelect )
      {
        piArlDstCoeff[uiBlkPos]   = (TCoeff)(( lLevelDouble + iAddC) >> iQBitsC );
      }
#endif
      const UInt uiMaxAbsLevel  = std::min<UInt>(UInt(entropyCodingMaximum), UInt((lLevelDouble + (Intermediate_Int(1) << (iQBits - 1))) >> iQBits));

      const Double dErr         = Double( lLevelDouble );
      pdCostCoeff0[ iScanPos ]  = dErr * dErr * errorScale;
      d64BlockUncodedCost      += pdCostCoeff0[ iScanPos ];
      piDstCoeff[ uiBlkPos ]    = uiMaxAbsLevel;

      if ( uiMaxAbsLevel > 0 && iLastScanPos < 0 )
      {
        iLastScanPos            = iScanPos;
#if !VCEG_AZ07_CTX_RESIDUALCODING
        uiCtxSet                = getContextSetIndex(compID, (iScanPos >> MLS_CG_SIZE), 0);
#endif
        iCGLastScanPos          = iCGScanPos;
      }

      if ( iLastScanPos >= 0 )
      {
        //===== coefficient level estimation =====
        UInt  uiLevel;
#if VCEG_AZ07_CTX_RESIDUALCODING
        UInt  uiOneCtx = 0;
        UInt  uiAbsCtx = 0;
#else
        UInt  uiOneCtx         = (NUM_ONE_FLAG_CTX_PER_SET * uiCtxSet) + c1;
        UInt  uiAbsCtx         = (NUM_ABS_FLAG_CTX_PER_SET * uiCtxSet) + c2;
#endif
        if( iScanPos == iLastScanPos )
        {
          uiLevel              = xGetCodedLevel( pdCostCoeff[ iScanPos ], pdCostCoeff0[ iScanPos ], pdCostSig[ iScanPos ],
                                                  lLevelDouble, uiMaxAbsLevel, significanceMapContextOffset, uiOneCtx, uiAbsCtx, uiGoRiceParam,
                                                  c1Idx, c2Idx, iQBits, errorScale, 1, extendedPrecision, maxLog2TrDynamicRange
                                                  );
#if ITH_C1_RDOQ_MEM_OPT
		  sigRateDelta[uiBlkPos] = 0;
#endif
        }
        else
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          UShort uiCtxSig      = uiOffsetonTU + getSigCtxInc( piDstCoeff, codingParameters.scan[iScanPos], uiWidth, uiHeight, channelType, uiOneCtx, uiAbsCtx, uiGoRiceParam ) ;
          if(bUseGolombRiceParameterAdaptation)
          {
            uiGoRiceParam = max(uiPrevGRParam, uiGoRiceParam);
          }
#else
          UShort uiCtxSig      = significanceMapContextOffset + getSigCtxInc( patternSigCtx, codingParameters, iScanPos, uiLog2BlockWidth, uiLog2BlockHeight, channelType );
#endif
          uiLevel              = xGetCodedLevel( pdCostCoeff[ iScanPos ], pdCostCoeff0[ iScanPos ], pdCostSig[ iScanPos ],
                                                  lLevelDouble, uiMaxAbsLevel, uiCtxSig, uiOneCtx, uiAbsCtx, uiGoRiceParam,
                                                  c1Idx, c2Idx, iQBits, errorScale, 0, extendedPrecision, maxLog2TrDynamicRange
                                                  );

          sigRateDelta[ uiBlkPos ] = m_pcEstBitsSbac->significantBits[ uiCtxSig ][ 1 ] - m_pcEstBitsSbac->significantBits[ uiCtxSig ][ 0 ];
          
        }

#if VCEG_AZ07_CTX_RESIDUALCODING
        uiPrevGRParam = max(0, (Int)uiGoRiceParam - 1);
#endif
        deltaU[ uiBlkPos ]        = TCoeff((lLevelDouble - (Intermediate_Int(uiLevel) << iQBits)) >> (iQBits-8));

        if( uiLevel > 0 )
        {
          Int rateNow = xGetICRate( uiLevel, uiOneCtx, uiAbsCtx, uiGoRiceParam, c1Idx, c2Idx, extendedPrecision, maxLog2TrDynamicRange );
          rateIncUp   [ uiBlkPos ] = xGetICRate( uiLevel+1, uiOneCtx, uiAbsCtx, uiGoRiceParam, c1Idx, c2Idx, extendedPrecision, maxLog2TrDynamicRange ) - rateNow;
          rateIncDown [ uiBlkPos ] = xGetICRate( uiLevel-1, uiOneCtx, uiAbsCtx, uiGoRiceParam, c1Idx, c2Idx, extendedPrecision, maxLog2TrDynamicRange ) - rateNow;
        }
        else // uiLevel == 0
        {
          rateIncUp   [ uiBlkPos ] = m_pcEstBitsSbac->m_greaterOneBits[ uiOneCtx ][ 0 ];
#if ITH_C1_RDOQ_MEM_OPT
		  rateIncDown[uiBlkPos] = 0;
#endif
        }
        piDstCoeff[ uiBlkPos ] = uiLevel;
        d64BaseCost           += pdCostCoeff [ iScanPos ];

#if !VCEG_AZ07_CTX_RESIDUALCODING
        baseLevel = (c1Idx < C1FLAG_NUMBER) ? (2 + (c2Idx < C2FLAG_NUMBER)) : 1;
        if( uiLevel >= baseLevel )
        {
          if (uiLevel > 3*(1<<uiGoRiceParam))
          {
            uiGoRiceParam = bUseGolombRiceParameterAdaptation ? (uiGoRiceParam + 1) : (std::min<UInt>((uiGoRiceParam + 1), 4));
          }
        }
#endif
        if ( uiLevel >= 1)
        {
          c1Idx ++;
#if VCEG_AZ07_CTX_RESIDUALCODING
          if( uiLevel>1)
          {
            c2Idx ++;
          }
#endif
        }

        //===== update bin model =====
#if !VCEG_AZ07_CTX_RESIDUALCODING
        if( uiLevel > 1 )
        {
          c1 = 0;
          c2 += (c2 < 2);
          c2Idx ++;
        }
        else if( (c1 < 3) && (c1 > 0) && uiLevel)
        {
          c1++;
        }
#endif

        //===== context set update =====
        if( ( iScanPos % uiCGSize == 0 ) && ( iScanPos > 0 ) )
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          c1Idx = 0;
          c2Idx = 0; 
#else
          uiCtxSet          = getContextSetIndex(compID, ((iScanPos - 1) >> MLS_CG_SIZE), (c1 == 0)); //(iScanPos - 1) because we do this **before** entering the final group
          c1                = 1;
          c2                = 0;
          c1Idx             = 0;
          c2Idx             = 0;
          uiGoRiceParam     = initialGolombRiceParameter;
#endif
        }
      }
      else
      {
        d64BaseCost    += pdCostCoeff0[ iScanPos ];
      }
      
      rdStats.d64SigCost += pdCostSig[ iScanPos ];
      if (iScanPosinCG == 0 )
      {
        rdStats.d64SigCost_0 = pdCostSig[ iScanPos ];
      }
      if (piDstCoeff[ uiBlkPos ] )
      {
        uiSigCoeffGroupFlag[ uiCGBlkPos ] = 1;
        rdStats.d64CodedLevelandDist += pdCostCoeff[ iScanPos ] - pdCostSig[ iScanPos ];
        rdStats.d64UncodedDist += pdCostCoeff0[ iScanPos ];
        if ( iScanPosinCG != 0 )
        {
          rdStats.iNNZbeforePos0++;
        }
      }
    } //end for (iScanPosinCG)

    if (iCGLastScanPos >= 0)
    {
      if( iCGScanPos )
      {
        if (uiSigCoeffGroupFlag[ uiCGBlkPos ] == 0)
        {
#if VCEG_AZ07_CTX_RESIDUALCODING
          UInt  uiCtxSig = getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, uiCGPosX, uiCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups, codingParameters.scanType );
#else
          UInt  uiCtxSig = getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, uiCGPosX, uiCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups );
#endif
          d64BaseCost += xGetRateSigCoeffGroup(0, uiCtxSig) - rdStats.d64SigCost;;
          pdCostCoeffGroupSig[ iCGScanPos ] = xGetRateSigCoeffGroup(0, uiCtxSig);
        }
        else
        {
          if (iCGScanPos < iCGLastScanPos) //skip the last coefficient group, which will be handled together with last position below.
          {
            if ( rdStats.iNNZbeforePos0 == 0 )
            {
              d64BaseCost -= rdStats.d64SigCost_0;
              rdStats.d64SigCost -= rdStats.d64SigCost_0;
            }
            // rd-cost if SigCoeffGroupFlag = 0, initialization
            Double d64CostZeroCG = d64BaseCost;

            // add SigCoeffGroupFlag cost to total cost
#if VCEG_AZ07_CTX_RESIDUALCODING
            UInt  uiCtxSig = getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, uiCGPosX, uiCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups, codingParameters.scanType );
#else
            UInt  uiCtxSig = getSigCoeffGroupCtxInc( uiSigCoeffGroupFlag, uiCGPosX, uiCGPosY, codingParameters.widthInGroups, codingParameters.heightInGroups );
#endif
            if (iCGScanPos < iCGLastScanPos)
            {
              d64BaseCost  += xGetRateSigCoeffGroup(1, uiCtxSig);
              d64CostZeroCG += xGetRateSigCoeffGroup(0, uiCtxSig);
              pdCostCoeffGroupSig[ iCGScanPos ] = xGetRateSigCoeffGroup(1, uiCtxSig);
            }

            // try to convert the current coeff group from non-zero to all-zero
            d64CostZeroCG += rdStats.d64UncodedDist;  // distortion for resetting non-zero levels to zero levels
            d64CostZeroCG -= rdStats.d64CodedLevelandDist;   // distortion and level cost for keeping all non-zero levels
            d64CostZeroCG -= rdStats.d64SigCost;     // sig cost for all coeffs, including zero levels and non-zerl levels

            // if we can save cost, change this block to all-zero block
            if ( d64CostZeroCG < d64BaseCost )
            {
              uiSigCoeffGroupFlag[ uiCGBlkPos ] = 0;
              d64BaseCost = d64CostZeroCG;
              if (iCGScanPos < iCGLastScanPos)
              {
                pdCostCoeffGroupSig[ iCGScanPos ] = xGetRateSigCoeffGroup(0, uiCtxSig);
              }
              // reset coeffs to 0 in this block
              for (Int iScanPosinCG = uiCGSize-1; iScanPosinCG >= 0; iScanPosinCG--)
              {
                iScanPos      = iCGScanPos*uiCGSize + iScanPosinCG;
                UInt uiBlkPos = codingParameters.scan[ iScanPos ];

                if (piDstCoeff[ uiBlkPos ])
                {
                  piDstCoeff [ uiBlkPos ] = 0;
                  pdCostCoeff[ iScanPos ] = pdCostCoeff0[ iScanPos ];
                  pdCostSig  [ iScanPos ] = 0;
                }
              }
            } // end if ( d64CostAllZeros < d64BaseCost )
          }
        } // end if if (uiSigCoeffGroupFlag[ uiCGBlkPos ] == 0)
      }
      else
      {
        uiSigCoeffGroupFlag[ uiCGBlkPos ] = 1;
      }
    }
  } //end for (iCGScanPos)

  //===== estimate last position =====
  if ( iLastScanPos < 0 )
  {
#if JVET_C0024_QTBT && !ITH_C1_RDOQ_MEM_OPT
    delete [] pdCostCoeff;
    delete [] pdCostSig;
    delete [] pdCostCoeff0;
    delete [] rateIncUp;
    delete [] rateIncDown;
    delete [] sigRateDelta;
    delete [] deltaU;
#endif
    return;
  }

  Double  d64BestCost         = 0;
  Int     ui16CtxCbf          = 0;
  Int     iBestLastIdxP1      = 0;
#if JVET_C0024_QTBT
  if( !pcCU->isIntra( uiAbsPartIdx ) && isLuma(compID))
#else
  if( !pcCU->isIntra( uiAbsPartIdx ) && isLuma(compID) && pcCU->getTransformIdx( uiAbsPartIdx ) == 0 )
#endif
  {
    ui16CtxCbf   = 0;
    d64BestCost  = d64BlockUncodedCost + xGetICost( m_pcEstBitsSbac->blockRootCbpBits[ ui16CtxCbf ][ 0 ] );
    d64BaseCost += xGetICost( m_pcEstBitsSbac->blockRootCbpBits[ ui16CtxCbf ][ 1 ] );
  }
  else
  {
    ui16CtxCbf   = pcCU->getCtxQtCbf( rTu, channelType );
    ui16CtxCbf  += getCBFContextOffset(compID);
    d64BestCost  = d64BlockUncodedCost + xGetICost( m_pcEstBitsSbac->blockCbpBits[ ui16CtxCbf ][ 0 ] );
    d64BaseCost += xGetICost( m_pcEstBitsSbac->blockCbpBits[ ui16CtxCbf ][ 1 ] );
  }


  Bool bFoundLast = false;
  for (Int iCGScanPos = iCGLastScanPos; iCGScanPos >= 0; iCGScanPos--)
  {
    UInt uiCGBlkPos = codingParameters.scanCG[ iCGScanPos ];

    d64BaseCost -= pdCostCoeffGroupSig [ iCGScanPos ];
    if (uiSigCoeffGroupFlag[ uiCGBlkPos ])
    {
      for (Int iScanPosinCG = uiCGSize-1; iScanPosinCG >= 0; iScanPosinCG--)
      {
        iScanPos = iCGScanPos*uiCGSize + iScanPosinCG;

        if (iScanPos > iLastScanPos)
        {
          continue;
        }
        UInt   uiBlkPos     = codingParameters.scan[iScanPos];

        if( piDstCoeff[ uiBlkPos ] )
        {
          UInt   uiPosY       = uiBlkPos >> uiLog2BlockWidth;
          UInt   uiPosX       = uiBlkPos - ( uiPosY << uiLog2BlockWidth );

          Double d64CostLast= codingParameters.scanType == SCAN_VER ? xGetRateLast( uiPosY, uiPosX, compID ) : xGetRateLast( uiPosX, uiPosY, compID );
          Double totalCost = d64BaseCost + d64CostLast - pdCostSig[ iScanPos ];

          if( totalCost < d64BestCost )
          {
            iBestLastIdxP1  = iScanPos + 1;
            d64BestCost     = totalCost;
          }
          if( piDstCoeff[ uiBlkPos ] > 1 )
          {
            bFoundLast = true;
            break;
          }
          d64BaseCost      -= pdCostCoeff[ iScanPos ];
          d64BaseCost      += pdCostCoeff0[ iScanPos ];
        }
        else
        {
          d64BaseCost      -= pdCostSig[ iScanPos ];
        }
      } //end for
      if (bFoundLast)
      {
        break;
      }
    } // end if (uiSigCoeffGroupFlag[ uiCGBlkPos ])
  } // end for

#if JVET_C0024_QTBT && !ITH_C1_RDOQ_MEM_OPT
  delete [] pdCostCoeff;
  delete [] pdCostSig;
  delete [] pdCostCoeff0;
#endif

  for ( Int scanPos = 0; scanPos < iBestLastIdxP1; scanPos++ )
  {
    Int blkPos = codingParameters.scan[ scanPos ];
    TCoeff level = piDstCoeff[ blkPos ];
    uiAbsSum += level;
    piDstCoeff[ blkPos ] = ( plSrcCoeff[ blkPos ] < 0 ) ? -level : level;
  }

  //===== clean uncoded coefficients =====
  for ( Int scanPos = iBestLastIdxP1; scanPos <= iLastScanPos; scanPos++ )
  {
    piDstCoeff[ codingParameters.scan[ scanPos ] ] = 0;
  }


  if( pcCU->getSlice()->getPPS()->getSignHideFlag() && uiAbsSum>=2)
  {
    const Double inverseQuantScale = Double(g_invQuantScales[cQP.rem]);
    Int64 rdFactor = (Int64)(inverseQuantScale * inverseQuantScale * (1 << (2 * cQP.per))
                             / m_dLambda / 16 / (1 << (2 * DISTORTION_PRECISION_ADJUSTMENT(channelBitDepth - 8)))
                             + 0.5);

    Int lastCG = -1;
    Int absSum = 0 ;

    Int n ;

    for( Int subSet = (uiWidth*uiHeight-1) >> MLS_CG_SIZE; subSet >= 0; subSet-- )
    {
      Int  subPos     = subSet << MLS_CG_SIZE;
      Int  firstNZPosInCG=uiCGSize , lastNZPosInCG=-1 ;
      absSum = 0 ;

      for(n = uiCGSize-1; n >= 0; --n )
      {
        if( piDstCoeff[ codingParameters.scan[ n + subPos ]] )
        {
          lastNZPosInCG = n;
          break;
        }
      }

      for(n = 0; n <uiCGSize; n++ )
      {
        if( piDstCoeff[ codingParameters.scan[ n + subPos ]] )
        {
          firstNZPosInCG = n;
          break;
        }
      }

      for(n = firstNZPosInCG; n <=lastNZPosInCG; n++ )
      {
        absSum += Int(piDstCoeff[ codingParameters.scan[ n + subPos ]]);
      }

      if(lastNZPosInCG>=0 && lastCG==-1)
      {
        lastCG = 1;
      }

      if( lastNZPosInCG-firstNZPosInCG>=SBH_THRESHOLD )
      {
        UInt signbit = (piDstCoeff[codingParameters.scan[subPos+firstNZPosInCG]]>0?0:1);
        if( signbit!=(absSum&0x1) )  // hide but need tune
        {
          // calculate the cost
          Int64 minCostInc = std::numeric_limits<Int64>::max(), curCost = std::numeric_limits<Int64>::max();
          Int minPos = -1, finalChange = 0, curChange = 0;

          for( n = (lastCG==1?lastNZPosInCG:uiCGSize-1) ; n >= 0; --n )
          {
            UInt uiBlkPos   = codingParameters.scan[ n + subPos ];
            if(piDstCoeff[ uiBlkPos ] != 0 )
            {
              Int64 costUp   = rdFactor * ( - deltaU[uiBlkPos] ) + rateIncUp[uiBlkPos];
              Int64 costDown = rdFactor * (   deltaU[uiBlkPos] ) + rateIncDown[uiBlkPos]
                               -   ((abs(piDstCoeff[uiBlkPos]) == 1) ? sigRateDelta[uiBlkPos] : 0);

              if(lastCG==1 && lastNZPosInCG==n && abs(piDstCoeff[uiBlkPos])==1)
              {
                costDown -= (4<<15);
              }

              if(costUp<costDown)
              {
                curCost = costUp;
                curChange =  1;
              }
              else
              {
                curChange = -1;
                if(n==firstNZPosInCG && abs(piDstCoeff[uiBlkPos])==1)
                {
                  curCost = std::numeric_limits<Int64>::max();
                }
                else
                {
                  curCost = costDown;
                }
              }
            }
            else
            {
              curCost = rdFactor * ( - (abs(deltaU[uiBlkPos])) ) + (1<<15) + rateIncUp[uiBlkPos] + sigRateDelta[uiBlkPos] ;
              curChange = 1 ;

              if(n<firstNZPosInCG)
              {
                UInt thissignbit = (plSrcCoeff[uiBlkPos]>=0?0:1);
                if(thissignbit != signbit )
                {
                  curCost = std::numeric_limits<Int64>::max();
                }
              }
            }

            if( curCost<minCostInc)
            {
              minCostInc = curCost;
              finalChange = curChange;
              minPos = uiBlkPos;
            }
          }

          if(piDstCoeff[minPos] == entropyCodingMaximum || piDstCoeff[minPos] == entropyCodingMinimum)
          {
            finalChange = -1;
          }

          if(plSrcCoeff[minPos]>=0)
          {
            piDstCoeff[minPos] += finalChange ;
          }
          else
          {
            piDstCoeff[minPos] -= finalChange ;
          }
        }
      }

      if(lastCG==1)
      {
        lastCG=0 ;
      }
    }
  }
#if JVET_C0024_QTBT && !ITH_C1_RDOQ_MEM_OPT
  delete [] rateIncUp;
  delete [] rateIncDown;
  delete [] sigRateDelta;
  delete [] deltaU;
#endif
}

#if VCEG_AZ07_CTX_RESIDUALCODING
Int TComTrQuant::getGrtZeroCtxInc( TCoeff*                         pcCoeff,
                                   const Int                       rasterPosition,
                                   Int                             width,
                                   Int                             height,
                                  const ChannelType                chanType                                  
                                  )
{
#if JVET_C0024_QTBT
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + MIN_CU_LOG2;
#else
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + 2;
#endif
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  const TCoeff *pData  = pcCoeff + posX + posY * width;
  const Int  widthM1   = width  - 1;
  const Int  heightM1  = height - 1;
  const Int  diag      = posX + posY;

  Int numPos = 0;
  if( posX < widthM1 )
  {
    numPos += pData[ 1 ] != 0;
    if( posX < widthM1 - 1 )
    {
      numPos += pData[ 2 ] != 0;
    }
    if( posY < heightM1 )
    {
       numPos += pData[ width + 1 ] != 0;
    }
  }
  if( posY < heightM1 )
  {
    numPos += pData[ width ] != 0;
    if( posY < heightM1 - 1 )
    {
      numPos += pData[ 2 * width ] != 0;
    }
  }

  const Int ctxIdx = min( numPos, 5 );
        Int ctxOfs = diag < 2 ? 6 : 0;

  if( chanType == CHANNEL_TYPE_LUMA )
  {
    ctxOfs += diag < 5 ? 6 : 0;
  }
  return ctxOfs + ctxIdx;
}
Int TComTrQuant::getGrtOneCtxInc ( TCoeff*                         pcCoeff,                                   
                                   const Int                       rasterPosition,
                                   Int                             width,
                                   Int                             height,
                                   const ChannelType               chanType
                                   )
{
#if JVET_C0024_QTBT
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + MIN_CU_LOG2;
#else
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + 2;
#endif
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  const TCoeff *pData  = pcCoeff + posX + (posY << log2BlockWidth);
  const Int  widthM1   = width  - 1;
  const Int  heightM1  = height - 1;
  UInt  diag           = posX + posY;
  Int sumAbsOne        = 0;

  if( posX < widthM1 )
  {
    sumAbsOne += abs( pData[ 1 ] )>1;
    if( posX < widthM1 - 1 )
    {
      sumAbsOne += abs( pData[ 2 ] )>1;
    }
    if( posY < heightM1 )
    {
      sumAbsOne += abs( pData[ width + 1 ] )>1;
    }
  }
  if( posY < heightM1 )
  {
    sumAbsOne += abs( pData[ width ] )>1;
    if( posY < heightM1 - 1 )
    {
      sumAbsOne += abs( pData[ 2 * width ] )>1;
    }
  }

   sumAbsOne       = min<UInt>( sumAbsOne, 4 ) + 1;
   if(chanType == CHANNEL_TYPE_LUMA)
   {
     sumAbsOne   += diag < 3 ? 10 : ( diag < 10 ? 5 : 0 );
   }

  return sumAbsOne;
}
Int TComTrQuant::getGrtTwoCtxInc ( TCoeff*                         pcCoeff,                                  
                                   const Int                       rasterPosition,
                                   Int                             width,
                                   Int                             height,
                                  const ChannelType                chanType
                                  )
{
#if JVET_C0024_QTBT
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + MIN_CU_LOG2;
#else
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + 2;
#endif
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  const TCoeff *pData  = pcCoeff + posX + (posY << log2BlockWidth);
  const Int  widthM1   = width  - 1;
  const Int  heightM1  = height - 1;
  UInt diag            = posX + posY; 

  Int sumAbsTwo = 0;
  if( posX < widthM1 )
  {
    sumAbsTwo += abs( pData[ 1 ] )>2;
    if( posX < widthM1 - 1 )
    {
      sumAbsTwo += abs( pData[ 2 ] )>2;
    }
    if( posY < heightM1 )
    {
      sumAbsTwo += abs( pData[ width + 1 ] )>2;
    }
  }
  if( posY < heightM1 )
  {
    sumAbsTwo += abs( pData[ width ] )>2;
    if( posY < heightM1 - 1 )
    {
      sumAbsTwo += abs( pData[ 2 * width ] )>2;
    }
  }
 
  sumAbsTwo  = min<UInt>( sumAbsTwo, 4 ) + 1;
  if( chanType == CHANNEL_TYPE_LUMA )
  {  
    sumAbsTwo += diag < 3 ? 10 : ( diag < 10 ? 5 : 0 );
  }

  return sumAbsTwo;
}
Int TComTrQuant::getRemainCoeffCtxInc(   TCoeff*                            pcCoeff,
                                   const Int                          rasterPosition,
                                         Int                                   width,
                                         Int                                   height
                                  )
{
#if JVET_C0024_QTBT
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + MIN_CU_LOG2;
#else
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + 2;
#endif
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  const TCoeff *pData  = pcCoeff + posX + (posY << log2BlockWidth);
  const Int  widthM1   = width  - 1;
  const Int  heightM1  = height - 1;
  
  Int numPos = 0;
  Int sumAbsAll = 0;
  if( posX < widthM1 )
  {
    sumAbsAll += abs( pData[ 1 ] );
    numPos += pData[ 1 ] != 0;
    if( posX < widthM1 - 1 )
    {
      sumAbsAll += abs( pData[ 2 ] );
      numPos += pData[ 2 ] != 0;
    }
    if( posY < heightM1 )
    {
      sumAbsAll += abs( pData[ width + 1 ] );
      numPos += pData[ width + 1 ] != 0;
    }
  }
  if( posY < heightM1 )
  {
    sumAbsAll += abs( pData[ width ] );
    numPos += pData[ width ] != 0;
    if( posY < heightM1 - 1 )
    {
      sumAbsAll += abs( pData[ 2 * width ] );
      numPos += pData[ 2 * width ] != 0;
    }
  }

  UInt uiVal = (sumAbsAll - numPos);
  UInt iOrder;
  for( iOrder = 0; iOrder < MAX_GR_ORDER_RESIDUAL; iOrder ++ )
  {
    if( (1 << (iOrder + 3)) > (uiVal + 4))
    {
      break;
    }
  }
  return ( iOrder == MAX_GR_ORDER_RESIDUAL ? (MAX_GR_ORDER_RESIDUAL - 1): iOrder);
}
Int TComTrQuant::getSigCtxInc    ( TCoeff*                         pcCoeff,
                                  const Int                        rasterPosition,
                                   Int                             width,
                                   Int                             height,
                                  const ChannelType                chanType,
                                   UInt&                           sumOne,
                                   UInt&                           sumTwo,
                                   UInt&                           sumAbs
                                  )
{
#if JVET_C0024_QTBT
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + MIN_CU_LOG2;
#else
  const UInt log2BlockWidth = g_aucConvertToBit[ width  ] + 2;
#endif
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  const TCoeff *pData  = pcCoeff + posX + (posY << log2BlockWidth);
  const Int  widthM1   = width  - 1;
  const Int  heightM1  = height - 1;
  const Int  diag      = posX + posY;

  Int sumAbsOne = 0;
  Int sumAbsTwo = 0;
  Int numPos = 0;
  Int sumAbsAll = 0;
  if( posX < widthM1 )
  {
    sumAbsOne += abs( pData[ 1 ] )>1;
    sumAbsTwo += abs( pData[ 1 ] )>2;
    sumAbsAll += abs( pData[ 1 ] );

    numPos += pData[ 1 ] != 0;
    if( posX < widthM1 - 1 )
    {
      sumAbsOne += abs( pData[ 2 ] )>1;
      sumAbsTwo += abs( pData[ 2 ] )>2;
      sumAbsAll += abs( pData[ 2 ] );
      numPos += pData[ 2 ] != 0;
    }
    if( posY < heightM1 )
    {
      sumAbsOne += abs( pData[ width + 1 ] )>1;
      sumAbsTwo += abs( pData[ width + 1 ] )>2;
      sumAbsAll += abs( pData[ width + 1 ] );
      numPos += pData[ width + 1 ] != 0;
    }
  }
  if( posY < heightM1 )
  {
    sumAbsOne += abs( pData[ width ] )>1;
    sumAbsTwo += abs( pData[ width ] )>2;
    sumAbsAll += abs( pData[ width ] );
    numPos += pData[ width ] != 0;
    if( posY < heightM1 - 1 )
    {
      sumAbsOne += abs( pData[ 2 * width ] )>1;
      sumAbsTwo += abs( pData[ 2 * width ] )>2;
      sumAbsAll += abs( pData[ 2 * width ] );
      numPos += pData[ 2 * width ] != 0;
    }
  }

  Int ctxIdx = min( numPos, 5 );
  Int ctxOfs = diag < 2 ? 6 : 0;
  sumOne     = min<UInt>( sumAbsOne, 4 ) + 1;
  sumTwo     = min<UInt>( sumAbsTwo, 4 ) + 1;

  UInt uiVal = (sumAbsAll - numPos);
  UInt iOrder;
  for( iOrder = 0; iOrder < MAX_GR_ORDER_RESIDUAL; iOrder ++ )
  {
    if( (1 << (iOrder + 3)) > (uiVal + 4))
    {
      break;
    }
  }
  sumAbs = ( iOrder == MAX_GR_ORDER_RESIDUAL ? (MAX_GR_ORDER_RESIDUAL - 1): iOrder);

  if( chanType == CHANNEL_TYPE_LUMA )
  {
    Int  iOfs = diag < 3 ? 10 : ( diag < 10 ? 5 : 0 );
    ctxOfs += ((diag < 5 ? 6 : 0));
    sumOne            +=  iOfs;
    sumTwo            +=  iOfs;
  }

  return ctxOfs + ctxIdx;
}
#else
/** Pattern decision for context derivation process of significant_coeff_flag
 * \param sigCoeffGroupFlag pointer to prior coded significant coeff group
 * \param uiCGPosX column of current coefficient group
 * \param uiCGPosY row of current coefficient group
 * \param widthInGroups width of the block
 * \param heightInGroups height of the block
 * \returns pattern for current coefficient group
 */
Int  TComTrQuant::calcPatternSigCtx( const UInt* sigCoeffGroupFlag, UInt uiCGPosX, UInt uiCGPosY, UInt widthInGroups, UInt heightInGroups )
{
  if ((widthInGroups <= 1) && (heightInGroups <= 1))
  {
    return 0;
  }

  const Bool rightAvailable = uiCGPosX < (widthInGroups  - 1);
  const Bool belowAvailable = uiCGPosY < (heightInGroups - 1);

  UInt sigRight = 0;
  UInt sigLower = 0;

  if (rightAvailable)
  {
    sigRight = ((sigCoeffGroupFlag[ (uiCGPosY * widthInGroups) + uiCGPosX + 1 ] != 0) ? 1 : 0);
  }
  if (belowAvailable)
  {
    sigLower = ((sigCoeffGroupFlag[ (uiCGPosY + 1) * widthInGroups + uiCGPosX ] != 0) ? 1 : 0);
  }

  return sigRight + (sigLower << 1);
}


/** Context derivation process of coeff_abs_significant_flag
 * \param patternSigCtx pattern for current coefficient group
 * \param codingParameters coding parameters for the TU (includes the scan)
 * \param scanPosition current position in scan order
 * \param log2BlockWidth log2 width of the block
 * \param log2BlockHeight log2 height of the block
 * \param chanType channel type (CHANNEL_TYPE_LUMA/CHROMA)
 * \returns ctxInc for current scan position
 */
Int TComTrQuant::getSigCtxInc    (       Int                        patternSigCtx,
                                   const TUEntropyCodingParameters &codingParameters,
                                   const Int                        scanPosition,
                                   const Int                        log2BlockWidth,
                                   const Int                        log2BlockHeight,
                                   const ChannelType                chanType)
{
  if (codingParameters.firstSignificanceMapContext == significanceMapContextSetStart[chanType][CONTEXT_TYPE_SINGLE])
  {
    //single context mode
    return significanceMapContextSetStart[chanType][CONTEXT_TYPE_SINGLE];
  }

  const UInt rasterPosition = codingParameters.scan[scanPosition];
  const UInt posY           = rasterPosition >> log2BlockWidth;
  const UInt posX           = rasterPosition - (posY << log2BlockWidth);

  if ((posX + posY) == 0)
  {
    return 0; //special case for the DC context variable
  }

  Int offset = MAX_INT;

  if ((log2BlockWidth == 2) && (log2BlockHeight == 2)) //4x4
  {
    offset = ctxIndMap4x4[ (4 * posY) + posX ];
  }
  else
  {
    Int cnt = 0;

    switch (patternSigCtx)
    {
      //------------------

      case 0: //neither neighbouring group is significant
        {
          const Int posXinSubset     = posX & ((1 << MLS_CG_LOG2_WIDTH)  - 1);
          const Int posYinSubset     = posY & ((1 << MLS_CG_LOG2_HEIGHT) - 1);
          const Int posTotalInSubset = posXinSubset + posYinSubset;

          //first N coefficients in scan order use 2; the next few use 1; the rest use 0.
          const UInt context1Threshold = NEIGHBOURHOOD_00_CONTEXT_1_THRESHOLD_4x4;
          const UInt context2Threshold = NEIGHBOURHOOD_00_CONTEXT_2_THRESHOLD_4x4;

          cnt = (posTotalInSubset >= context1Threshold) ? 0 : ((posTotalInSubset >= context2Threshold) ? 1 : 2);
        }
        break;

      //------------------

      case 1: //right group is significant, below is not
        {
          const Int posYinSubset = posY & ((1 << MLS_CG_LOG2_HEIGHT) - 1);
          const Int groupHeight  = 1 << MLS_CG_LOG2_HEIGHT;

          cnt = (posYinSubset >= (groupHeight >> 1)) ? 0 : ((posYinSubset >= (groupHeight >> 2)) ? 1 : 2); //top quarter uses 2; second-from-top quarter uses 1; bottom half uses 0
        }
        break;

      //------------------

      case 2: //below group is significant, right is not
        {
          const Int posXinSubset = posX & ((1 << MLS_CG_LOG2_WIDTH)  - 1);
          const Int groupWidth   = 1 << MLS_CG_LOG2_WIDTH;

          cnt = (posXinSubset >= (groupWidth >> 1)) ? 0 : ((posXinSubset >= (groupWidth >> 2)) ? 1 : 2); //left quarter uses 2; second-from-left quarter uses 1; right half uses 0
        }
        break;

      //------------------

      case 3: //both neighbouring groups are significant
        {
          cnt = 2;
        }
        break;

      //------------------

      default:
        std::cerr << "ERROR: Invalid patternSigCtx \"" << Int(patternSigCtx) << "\" in getSigCtxInc" << std::endl;
        exit(1);
        break;
    }

    //------------------------------------------------

    const Bool notFirstGroup = ((posX >> MLS_CG_LOG2_WIDTH) + (posY >> MLS_CG_LOG2_HEIGHT)) > 0;

    offset = (notFirstGroup ? notFirstGroupNeighbourhoodContextOffset[chanType] : 0) + cnt;
  }

  return codingParameters.firstSignificanceMapContext + offset;
}

#endif
/** Get the best level in RD sense
 *
 * \returns best quantized transform level for given scan position
 *
 * This method calculates the best quantized transform level for a given scan position.
 */
__inline UInt TComTrQuant::xGetCodedLevel ( Double&          rd64CodedCost,          //< reference to coded cost
                                            Double&          rd64CodedCost0,         //< reference to cost when coefficient is 0
                                            Double&          rd64CodedCostSig,       //< rd64CodedCostSig reference to cost of significant coefficient
                                            Intermediate_Int lLevelDouble,           //< reference to unscaled quantized level
                                            UInt             uiMaxAbsLevel,          //< scaled quantized level
                                            UShort           ui16CtxNumSig,          //< current ctxInc for coeff_abs_significant_flag
                                            UShort           ui16CtxNumOne,          //< current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
                                            UShort           ui16CtxNumAbs,          //< current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
                                            UShort           ui16AbsGoRice,          //< current Rice parameter for coeff_abs_level_minus3
                                            UInt             c1Idx,                  //< 
                                            UInt             c2Idx,                  //< 
                                            Int              iQBits,                 //< quantization step size
                                            Double           errorScale,             //< 
                                            Bool             bLast,                  //< indicates if the coefficient is the last significant
                                            Bool             useLimitedPrefixLength, //< 
                                            const Int        maxLog2TrDynamicRange   //< 
                                            ) const
{
  Double dCurrCostSig   = 0;
  UInt   uiBestAbsLevel = 0;

  if( !bLast && uiMaxAbsLevel < 3 )
  {
    rd64CodedCostSig    = xGetRateSigCoef( 0, ui16CtxNumSig );
    rd64CodedCost       = rd64CodedCost0 + rd64CodedCostSig;
    if( uiMaxAbsLevel == 0 )
    {
      return uiBestAbsLevel;
    }
  }
  else
  {
    rd64CodedCost       = MAX_DOUBLE;
  }

  if( !bLast )
  {
    dCurrCostSig        = xGetRateSigCoef( 1, ui16CtxNumSig );
  }

  UInt uiMinAbsLevel    = ( uiMaxAbsLevel > 1 ? uiMaxAbsLevel - 1 : 1 );
  for( Int uiAbsLevel  = uiMaxAbsLevel; uiAbsLevel >= uiMinAbsLevel ; uiAbsLevel-- )
  {
    Double dErr         = Double( lLevelDouble  - ( Intermediate_Int(uiAbsLevel) << iQBits ) );
    Double dCurrCost    = dErr * dErr * errorScale + xGetICost( xGetICRate( uiAbsLevel, ui16CtxNumOne, ui16CtxNumAbs, ui16AbsGoRice, c1Idx, c2Idx, useLimitedPrefixLength, maxLog2TrDynamicRange ) );
    dCurrCost          += dCurrCostSig;

    if( dCurrCost < rd64CodedCost )
    {
      uiBestAbsLevel    = uiAbsLevel;
      rd64CodedCost     = dCurrCost;
      rd64CodedCostSig  = dCurrCostSig;
    }
  }

  return uiBestAbsLevel;
}

/** Calculates the cost for specific absolute transform level
 * \param uiAbsLevel scaled quantized level
 * \param ui16CtxNumOne current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
 * \param ui16CtxNumAbs current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
 * \param ui16AbsGoRice Rice parameter for coeff_abs_level_minus3
 * \param c1Idx
 * \param c2Idx
 * \param useLimitedPrefixLength
 * \param maxLog2TrDynamicRange
 * \returns cost of given absolute transform level
 */
__inline Int TComTrQuant::xGetICRate         ( const UInt    uiAbsLevel,
                                               const UShort  ui16CtxNumOne,
                                               const UShort  ui16CtxNumAbs,
                                               const UShort  ui16AbsGoRice,
                                               const UInt    c1Idx,
                                               const UInt    c2Idx,
                                               const Bool    useLimitedPrefixLength,
                                               const Int     maxLog2TrDynamicRange
                                               ) const
{
  Int  iRate      = Int(xGetIEPRate()); // cost of sign bit
  UInt baseLevel  = (c1Idx < C1FLAG_NUMBER) ? (2 + (c2Idx < C2FLAG_NUMBER)) : 1;

  if ( uiAbsLevel >= baseLevel )
  {
    UInt symbol     = uiAbsLevel - baseLevel;
    UInt length;
#if VCEG_AZ07_CTX_RESIDUALCODING
    if ( symbol < (g_auiGoRiceRange[ui16AbsGoRice] << ui16AbsGoRice) )
#else
    if (symbol < (COEF_REMAIN_BIN_REDUCTION << ui16AbsGoRice))
#endif
    {
      length = symbol>>ui16AbsGoRice;
      iRate += (length+1+ui16AbsGoRice)<< 15;
    }
    else if (useLimitedPrefixLength)
    {
      const UInt maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange));

      UInt prefixLength = 0;
      UInt suffix       = (symbol >> ui16AbsGoRice) - COEF_REMAIN_BIN_REDUCTION;

      while ((prefixLength < maximumPrefixLength) && (suffix > ((2 << prefixLength) - 2)))
      {
        prefixLength++;
      }

      const UInt suffixLength = (prefixLength == maximumPrefixLength) ? (maxLog2TrDynamicRange - ui16AbsGoRice) : (prefixLength + 1/*separator*/);

      iRate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + ui16AbsGoRice) << 15;
    }
    else
    {
      length = ui16AbsGoRice;
#if VCEG_AZ07_CTX_RESIDUALCODING
      symbol  = symbol - ( g_auiGoRiceRange[ui16AbsGoRice] << ui16AbsGoRice);
#else
      symbol  = symbol - ( COEF_REMAIN_BIN_REDUCTION << ui16AbsGoRice);
#endif
      while (symbol >= (1<<length))
      {
        symbol -=  (1<<(length++));
      }
#if VCEG_AZ07_CTX_RESIDUALCODING
      iRate += (g_auiGoRiceRange[ui16AbsGoRice] + length + 1 - ui16AbsGoRice + length)<< 15;
#else
      iRate += (COEF_REMAIN_BIN_REDUCTION+length+1-ui16AbsGoRice+length)<< 15;
#endif
    }

    if (c1Idx < C1FLAG_NUMBER)
    {
      iRate += m_pcEstBitsSbac->m_greaterOneBits[ ui16CtxNumOne ][ 1 ];

      if (c2Idx < C2FLAG_NUMBER)
      {
#if VCEG_AZ07_CTX_RESIDUALCODING
        iRate += m_pcEstBitsSbac->m_greaterOneBits[ ui16CtxNumAbs ][ 1 ];
#else
        iRate += m_pcEstBitsSbac->m_levelAbsBits[ ui16CtxNumAbs ][ 1 ];
#endif
      }
    }
  }
  else if( uiAbsLevel == 1 )
  {
    iRate += m_pcEstBitsSbac->m_greaterOneBits[ ui16CtxNumOne ][ 0 ];
  }
  else if( uiAbsLevel == 2 )
  {
    iRate += m_pcEstBitsSbac->m_greaterOneBits[ ui16CtxNumOne ][ 1 ];
#if VCEG_AZ07_CTX_RESIDUALCODING
    iRate += m_pcEstBitsSbac->m_greaterOneBits[ ui16CtxNumAbs ][ 0 ];
#else
    iRate += m_pcEstBitsSbac->m_levelAbsBits[ ui16CtxNumAbs ][ 0 ];
#endif
  }
  else
  {
    iRate = 0;
  }

  return  iRate;
}

__inline Double TComTrQuant::xGetRateSigCoeffGroup  ( UShort                    uiSignificanceCoeffGroup,
                                                UShort                          ui16CtxNumSig ) const
{
  return xGetICost( m_pcEstBitsSbac->significantCoeffGroupBits[ ui16CtxNumSig ][ uiSignificanceCoeffGroup ] );
}

/** Calculates the cost of signaling the last significant coefficient in the block
 * \param uiPosX X coordinate of the last significant coefficient
 * \param uiPosY Y coordinate of the last significant coefficient
 * \param component colour component ID
 * \returns cost of last significant coefficient
 */
/*
 * \param uiWidth width of the transform unit (TU)
*/
__inline Double TComTrQuant::xGetRateLast   ( const UInt                      uiPosX,
                                              const UInt                      uiPosY,
                                              const ComponentID               component  ) const
{
  UInt uiCtxX   = g_uiGroupIdx[uiPosX];
  UInt uiCtxY   = g_uiGroupIdx[uiPosY];

  Double uiCost = m_pcEstBitsSbac->lastXBits[toChannelType(component)][ uiCtxX ] + m_pcEstBitsSbac->lastYBits[toChannelType(component)][ uiCtxY ];

  if( uiCtxX > 3 )
  {
    uiCost += xGetIEPRate() * ((uiCtxX-2)>>1);
  }
  if( uiCtxY > 3 )
  {
    uiCost += xGetIEPRate() * ((uiCtxY-2)>>1);
  }
  return xGetICost( uiCost );
}

__inline Double TComTrQuant::xGetRateSigCoef  ( UShort                          uiSignificance,
                                                UShort                          ui16CtxNumSig ) const
{
  return xGetICost( m_pcEstBitsSbac->significantBits[ ui16CtxNumSig ][ uiSignificance ] );
}

/** Get the cost for a specific rate
 * \param dRate rate of a bit
 * \returns cost at the specific rate
 */
__inline Double TComTrQuant::xGetICost        ( Double                          dRate         ) const
{
  return m_dLambda * dRate;
}

/** Get the cost of an equal probable bit
 * \returns cost of equal probable bit
 */
__inline Double TComTrQuant::xGetIEPRate      (                                               ) const
{
  return 32768;
}

/** Context derivation process of coeff_abs_significant_flag
 * \param uiSigCoeffGroupFlag significance map of L1
 * \param uiCGPosX column of current scan position
 * \param uiCGPosY row of current scan position
 * \param widthInGroups width of the block
 * \param heightInGroups height of the block
 * \returns ctxInc for current scan position
 */
UInt TComTrQuant::getSigCoeffGroupCtxInc  (const UInt*  uiSigCoeffGroupFlag,
                                           const UInt   uiCGPosX,
                                           const UInt   uiCGPosY,
#if VCEG_AZ07_CTX_RESIDUALCODING
                                                 UInt   widthInGroups,
                                                 UInt   heightInGroups,
                                           COEFF_SCAN_TYPE   scanIdx
#else
                                           const UInt   widthInGroups,
                                           const UInt   heightInGroups
#endif
                                           )
{
  UInt sigRight = 0;
  UInt sigLower = 0;
#if VCEG_AZ07_CTX_RESIDUALCODING
  if( widthInGroups == 2 && heightInGroups == 2 ) // 8x8
  {
    if( scanIdx == SCAN_HOR )  
    {
      widthInGroups  = 1;
      heightInGroups = 4;
    }
    else if( scanIdx == SCAN_VER )
    {
      widthInGroups  = 4;
      heightInGroups = 1;
    }
  }
#endif
  if (uiCGPosX < (widthInGroups  - 1))
  {
    sigRight = ((uiSigCoeffGroupFlag[ (uiCGPosY * widthInGroups) + uiCGPosX + 1 ] != 0) ? 1 : 0);
  }
  if (uiCGPosY < (heightInGroups - 1))
  {
    sigLower = ((uiSigCoeffGroupFlag[ (uiCGPosY + 1) * widthInGroups + uiCGPosX ] != 0) ? 1 : 0);
  }

  return ((sigRight + sigLower) != 0) ? 1 : 0;
}


/** set quantized matrix coefficient for encode
 * \param scalingList            quantized matrix address
 * \param format                 chroma format
 * \param maxLog2TrDynamicRange
 * \param bitDepths              reference to bit depth array for all channels
 */
Void TComTrQuant::setScalingList(TComScalingList *scalingList, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths)
{
  const Int minimumQp = 0;
  const Int maximumQp = SCALING_LIST_REM_NUM;

  for(UInt size = 0; size < SCALING_LIST_SIZE_NUM; size++)
  {
    for(UInt list = 0; list < SCALING_LIST_NUM; list++)
    {
      for(Int qp = minimumQp; qp < maximumQp; qp++)
      {
        xSetScalingListEnc(scalingList,list,size,qp);
        xSetScalingListDec(*scalingList,list,size,qp);
#if JVET_C0024_QTBT
        assert(0);  //no impact for CTC
        setErrScaleCoeff(list,size, size,qp,maxLog2TrDynamicRange, bitDepths);
#else
        setErrScaleCoeff(list,size,qp,maxLog2TrDynamicRange, bitDepths);
#endif
      }
    }
  }
}
/** set quantized matrix coefficient for decode
 * \param scalingList quantized matrix address
 * \param format      chroma format
 */
Void TComTrQuant::setScalingListDec(const TComScalingList &scalingList)
{
  const Int minimumQp = 0;
  const Int maximumQp = SCALING_LIST_REM_NUM;

  for(UInt size = 0; size < SCALING_LIST_SIZE_NUM; size++)
  {
    for(UInt list = 0; list < SCALING_LIST_NUM; list++)
    {
      for(Int qp = minimumQp; qp < maximumQp; qp++)
      {
        xSetScalingListDec(scalingList,list,size,qp);
      }
    }
  }
}
/** set error scale coefficients
 * \param list                   list ID
 * \param size                   
 * \param qp                     quantization parameter
 * \param maxLog2TrDynamicRange
 * \param bitDepths              reference to bit depth array for all channels
 */
#if JVET_C0024_QTBT
Void TComTrQuant::setErrScaleCoeff(UInt list, UInt w, UInt h, Int qp, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths)
#else
Void TComTrQuant::setErrScaleCoeff(UInt list, UInt size, Int qp, const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths)
#endif
{
#if JVET_C0024_QTBT
  const UInt uiLog2TrSize = (g_aucConvertToBit[ g_scalingListSizeX[w] ] + MIN_CU_LOG2 
      + g_aucConvertToBit[ g_scalingListSizeX[h] ] + MIN_CU_LOG2)>>1;
#else
  const UInt uiLog2TrSize = g_aucConvertToBit[ g_scalingListSizeX[size] ] + 2;
#endif
  const ChannelType channelType = ((list == 0) || (list == MAX_NUM_COMPONENT)) ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

  const Int channelBitDepth    = bitDepths.recon[channelType];
  const Int iTransformShift = getTransformShift(channelBitDepth, uiLog2TrSize, maxLog2TrDynamicRange[channelType]);  // Represents scaling through forward transform

#if JVET_C0024_QTBT
  Double dShift=0;
  if (((w+h)&1)!=0)
  {
      dShift = -0.5;
  }
  UInt i,uiMaxNumCoeff = g_scalingListSizeX[w]*g_scalingListSizeX[h];
  Int *piQuantcoeff;
  Double *pdErrScale;
  piQuantcoeff   = getQuantCoeff(list, qp, w, h);
  pdErrScale     = getErrScaleCoeff(list, w, h, qp);
#else
  UInt i,uiMaxNumCoeff = g_scalingListSize[size];
  Int *piQuantcoeff;
  Double *pdErrScale;
  piQuantcoeff   = getQuantCoeff(list, qp,size);
  pdErrScale     = getErrScaleCoeff(list, size, qp);
#endif

  Double dErrScale = (Double)(1<<SCALE_BITS);                                // Compensate for scaling of bitcount in Lagrange cost function
#if JVET_C0024_QTBT
  dErrScale = dErrScale*pow(2.0,(-2.0*(iTransformShift+dShift)));                     // Compensate for scaling through forward transform
#else
  dErrScale = dErrScale*pow(2.0,(-2.0*iTransformShift));                     // Compensate for scaling through forward transform
#endif

  for(i=0;i<uiMaxNumCoeff;i++)
  {
    pdErrScale[i] =  dErrScale / piQuantcoeff[i] / piQuantcoeff[i] / (1 << DISTORTION_PRECISION_ADJUSTMENT(2 * (bitDepths.recon[channelType] - 8)));
  }

#if JVET_C0024_QTBT
#if JVET_C0024_QTBT_FIX_QUANT_TICKET25
  if( (w+h)&1 )
  {
    Int quant = ( g_quantScales[qp] * 181 ) >> 7;
    getErrScaleCoeffNoScalingList(list, w, h, qp) = dErrScale / quant / quant / (1 << DISTORTION_PRECISION_ADJUSTMENT(2 * (bitDepths.recon[channelType] - 8)));
  }
  else
#endif
  getErrScaleCoeffNoScalingList(list, w, h, qp) = dErrScale / g_quantScales[qp] / g_quantScales[qp] / (1 << DISTORTION_PRECISION_ADJUSTMENT(2 * (bitDepths.recon[channelType] - 8)));
#else
  getErrScaleCoeffNoScalingList(list, size, qp) = dErrScale / g_quantScales[qp] / g_quantScales[qp] / (1 << DISTORTION_PRECISION_ADJUSTMENT(2 * (bitDepths.recon[channelType] - 8)));
#endif
}

/** set quantized matrix coefficient for encode
 * \param scalingList quantized matrix address
 * \param listId List index
 * \param sizeId size index
 * \param qp Quantization parameter
 * \param format chroma format
 */
Void TComTrQuant::xSetScalingListEnc(TComScalingList *scalingList, UInt listId, UInt sizeId, Int qp)
{
#if JVET_C0024_QTBT
  assert(0);  //scaling list, no impact for CTC, need to change later
#endif
  UInt width  = g_scalingListSizeX[sizeId];
  UInt height = g_scalingListSizeX[sizeId];
  UInt ratio  = g_scalingListSizeX[sizeId]/min(MAX_MATRIX_SIZE_NUM,(Int)g_scalingListSizeX[sizeId]);
  Int *quantcoeff;
  Int *coeff  = scalingList->getScalingListAddress(sizeId,listId);
#if JVET_C0024_QTBT 
  quantcoeff  = getQuantCoeff(listId, qp, sizeId, sizeId);
#else
  quantcoeff  = getQuantCoeff(listId, qp, sizeId);
#endif

  Int quantScales = g_quantScales[qp];

  processScalingListEnc(coeff,
                        quantcoeff,
                        (quantScales << LOG2_SCALING_LIST_NEUTRAL_VALUE),
                        height, width, ratio,
                        min(MAX_MATRIX_SIZE_NUM, (Int)g_scalingListSizeX[sizeId]),
                        scalingList->getScalingListDC(sizeId,listId));
}

/** set quantized matrix coefficient for decode
 * \param scalingList quantaized matrix address
 * \param listId List index
 * \param sizeId size index
 * \param qp Quantization parameter
 * \param format chroma format
 */
Void TComTrQuant::xSetScalingListDec(const TComScalingList &scalingList, UInt listId, UInt sizeId, Int qp)
{
  UInt width  = g_scalingListSizeX[sizeId];
  UInt height = g_scalingListSizeX[sizeId];
  UInt ratio  = g_scalingListSizeX[sizeId]/min(MAX_MATRIX_SIZE_NUM,(Int)g_scalingListSizeX[sizeId]);
  Int *dequantcoeff;
  const Int *coeff  = scalingList.getScalingListAddress(sizeId,listId);

#if JVET_C0024_QTBT   //scaling list, no impact for CTC, JCA
  assert(0);
  dequantcoeff = getDequantCoeff(listId, qp, sizeId, sizeId);
#else
  dequantcoeff = getDequantCoeff(listId, qp, sizeId);
#endif

  Int invQuantScale = g_invQuantScales[qp];

  processScalingListDec(coeff,
                        dequantcoeff,
                        invQuantScale,
                        height, width, ratio,
                        min(MAX_MATRIX_SIZE_NUM, (Int)g_scalingListSizeX[sizeId]),
                        scalingList.getScalingListDC(sizeId,listId));
}

/** set flat matrix value to quantized coefficient
 */
Void TComTrQuant::setFlatScalingList(const Int maxLog2TrDynamicRange[MAX_NUM_CHANNEL_TYPE], const BitDepths &bitDepths)
{
  const Int minimumQp = 0;
  const Int maximumQp = SCALING_LIST_REM_NUM;

#if JVET_C0024_QTBT
  for (UInt w=0; w<SCALING_LIST_SIZE_NUM; w++)
  {
  for (UInt h=0; h<SCALING_LIST_SIZE_NUM; h++)
#else
  for(UInt size = 0; size < SCALING_LIST_SIZE_NUM; size++)
#endif
  {
    for(UInt list = 0; list < SCALING_LIST_NUM; list++)
    {
      for(Int qp = minimumQp; qp < maximumQp; qp++)
      {
#if JVET_C0024_QTBT
        xsetFlatScalingList(list,w,h,qp);
        setErrScaleCoeff(list,w,h,qp,maxLog2TrDynamicRange, bitDepths);
#else
        xsetFlatScalingList(list,size,qp);
        setErrScaleCoeff(list,size,qp,maxLog2TrDynamicRange, bitDepths);
#endif
      }
    }
  }
#if JVET_C0024_QTBT
}
#endif
}

/** set flat matrix value to quantized coefficient
 * \param list List ID
 * \param size size index
 * \param qp Quantization parameter
 * \param format chroma format
 */
#if JVET_C0024_QTBT
Void TComTrQuant::xsetFlatScalingList(UInt list, UInt w, UInt h, Int qp)
#else
Void TComTrQuant::xsetFlatScalingList(UInt list, UInt size, Int qp)
#endif
{
#if JVET_C0024_QTBT
  UInt i,num = g_scalingListSizeX[w]*g_scalingListSizeX[h];
#else
  UInt i,num = g_scalingListSize[size];
#endif
  Int *quantcoeff;
  Int *dequantcoeff;

  Int quantScales    = g_quantScales   [qp];
  Int invQuantScales = g_invQuantScales[qp] << 4;

#if JVET_C0024_QTBT
  quantcoeff   = getQuantCoeff(list, qp, w, h);
  dequantcoeff = getDequantCoeff(list, qp, w, h);
#else
  quantcoeff   = getQuantCoeff(list, qp, size);
  dequantcoeff = getDequantCoeff(list, qp, size);
#endif

  for(i=0;i<num;i++)
  {
    *quantcoeff++ = quantScales;
    *dequantcoeff++ = invQuantScales;
  }
}

/** set quantized matrix coefficient for encode
 * \param coeff quantaized matrix address
 * \param quantcoeff quantaized matrix address
 * \param quantScales Q(QP%6)
 * \param height height
 * \param width width
 * \param ratio ratio for upscale
 * \param sizuNum matrix size
 * \param dc dc parameter
 */
Void TComTrQuant::processScalingListEnc( Int *coeff, Int *quantcoeff, Int quantScales, UInt height, UInt width, UInt ratio, Int sizuNum, UInt dc)
{
  for(UInt j=0;j<height;j++)
  {
    for(UInt i=0;i<width;i++)
    {
      quantcoeff[j*width + i] = quantScales / coeff[sizuNum * (j / ratio) + i / ratio];
    }
  }

  if(ratio > 1)
  {
    quantcoeff[0] = quantScales / dc;
  }
}

/** set quantized matrix coefficient for decode
 * \param coeff quantaized matrix address
 * \param dequantcoeff quantaized matrix address
 * \param invQuantScales IQ(QP%6))
 * \param height height
 * \param width width
 * \param ratio ratio for upscale
 * \param sizuNum matrix size
 * \param dc dc parameter
 */
Void TComTrQuant::processScalingListDec( const Int *coeff, Int *dequantcoeff, Int invQuantScales, UInt height, UInt width, UInt ratio, Int sizuNum, UInt dc)
{
  for(UInt j=0;j<height;j++)
  {
    for(UInt i=0;i<width;i++)
    {
      dequantcoeff[j*width + i] = invQuantScales * coeff[sizuNum * (j / ratio) + i / ratio];
    }
  }

  if(ratio > 1)
  {
    dequantcoeff[0] = invQuantScales * dc;
  }
}

/** initialization process of scaling list array
 */
Void TComTrQuant::initScalingList()
{
#if JVET_C0024_QTBT
    for(UInt w = 0; w < SCALING_LIST_SIZE_NUM; w++)
    {
        for (UInt h = 0; h < SCALING_LIST_SIZE_NUM; h++)
        {
#else
  for(UInt sizeId = 0; sizeId < SCALING_LIST_SIZE_NUM; sizeId++)
  {
#endif
    for(UInt qp = 0; qp < SCALING_LIST_REM_NUM; qp++)
    {
      for(UInt listId = 0; listId < SCALING_LIST_NUM; listId++)
      {
#if JVET_C0024_QTBT
                m_quantCoef   [w][h][listId][qp] = new Int    [g_scalingListSizeX[w]*g_scalingListSizeX[h]];
                m_dequantCoef [w][h][listId][qp] = new Int    [g_scalingListSizeX[w]*g_scalingListSizeX[h]];
                m_errScale    [w][h][listId][qp] = new Double [g_scalingListSizeX[w]*g_scalingListSizeX[h]];
#else
        m_quantCoef   [sizeId][listId][qp] = new Int    [g_scalingListSize[sizeId]];
        m_dequantCoef [sizeId][listId][qp] = new Int    [g_scalingListSize[sizeId]];
        m_errScale    [sizeId][listId][qp] = new Double [g_scalingListSize[sizeId]];
#endif
      } // listID loop
    }
  }
#if JVET_C0024_QTBT
}
#endif
}

/** destroy quantization matrix array
 */
Void TComTrQuant::destroyScalingList()
{
#if JVET_C0024_QTBT
  for(UInt w = 0; w < SCALING_LIST_SIZE_NUM; w++)
  {
      for (UInt h = 0; h < SCALING_LIST_SIZE_NUM; h++)
      {
#else
  for(UInt sizeId = 0; sizeId < SCALING_LIST_SIZE_NUM; sizeId++)
  {
#endif
    for(UInt listId = 0; listId < SCALING_LIST_NUM; listId++)
    {
      for(UInt qp = 0; qp < SCALING_LIST_REM_NUM; qp++)
      {
#if JVET_C0024_QTBT
        if(m_quantCoef[w][h][listId][qp])
        {
          delete [] m_quantCoef[w][h][listId][qp];
        }
        if(m_dequantCoef[w][h][listId][qp])
        {
          delete [] m_dequantCoef[w][h][listId][qp];
        }
        if(m_errScale[w][h][listId][qp])
        {
          delete [] m_errScale[w][h][listId][qp];
        }
      }
#else
        if(m_quantCoef[sizeId][listId][qp])
        {
          delete [] m_quantCoef[sizeId][listId][qp];
        }
        if(m_dequantCoef[sizeId][listId][qp])
        {
          delete [] m_dequantCoef[sizeId][listId][qp];
        }
        if(m_errScale[sizeId][listId][qp])
        {
          delete [] m_errScale[sizeId][listId][qp];
        }
#endif
      }
    }
  }
}

Void TComTrQuant::transformSkipQuantOneSample(TComTU &rTu, const ComponentID compID, const TCoeff resiDiff, TCoeff* pcCoeff, const UInt uiPos, const QpParam &cQP, const Bool bUseHalfRoundingPoint)
{
        TComDataCU    *pcCU                           = rTu.getCU();
  const UInt           uiAbsPartIdx                   = rTu.GetAbsPartIdxTU();
  const TComRectangle &rect                           = rTu.getRect(compID);
  const UInt           uiWidth                        = rect.width;
  const UInt           uiHeight                       = rect.height;
  const Int            maxLog2TrDynamicRange          = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));
  const Int            channelBitDepth                = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
  const Int            iTransformShift                = getTransformShift(channelBitDepth, rTu.GetEquivalentLog2TrSize(compID), maxLog2TrDynamicRange);
  const Int            scalingListType                = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
  const Bool           enableScalingLists             = getUseScalingList(uiWidth, uiHeight, true);
  const Int            defaultQuantisationCoefficient = g_quantScales[cQP.rem];

  assert( scalingListType < SCALING_LIST_NUM );
#if JVET_C0024_QTBT
  const UInt uiLog2TrWidth = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
  const UInt uiLog2TrHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
  const Int *const piQuantCoeff = getQuantCoeff( scalingListType, cQP.rem, uiLog2TrWidth-1, uiLog2TrHeight-1 );
#else
  const Int *const piQuantCoeff = getQuantCoeff( scalingListType, cQP.rem, (rTu.GetEquivalentLog2TrSize(compID)-2) );
#endif


  /* for 422 chroma blocks, the effective scaling applied during transformation is not a power of 2, hence it cannot be
  * implemented as a bit-shift (the quantised result will be sqrt(2) * larger than required). Alternatively, adjust the
  * uiLog2TrSize applied in iTransformShift, such that the result is 1/sqrt(2) the required result (i.e. smaller)
  * Then a QP+3 (sqrt(2)) or QP-3 (1/sqrt(2)) method could be used to get the required result
  */

  const Int iQBits = QUANT_SHIFT + cQP.per + iTransformShift;
  // QBits will be OK for any internal bit depth as the reduction in transform shift is balanced by an increase in Qp_per due to QpBDOffset

#if JVET_C0024_ENCODER_OVERFLOW_FIX
  const Int iAdd = (Int64)( bUseHalfRoundingPoint ? 256 : (pcCU->getSlice()->getSliceType() == I_SLICE ? 171 : 85) ) << ((Int64)iQBits - 9);
#else
  const Int iAdd = ( bUseHalfRoundingPoint ? 256 : (pcCU->getSlice()->getSliceType() == I_SLICE ? 171 : 85) ) << (iQBits - 9);
#endif

  TCoeff transformedCoefficient;

  // transform-skip
  if (iTransformShift >= 0)
  {
    transformedCoefficient = resiDiff << iTransformShift;
  }
  else // for very high bit depths
  {
    const Int iTrShiftNeg  = -iTransformShift;
    const Int offset       = 1 << (iTrShiftNeg - 1);
    transformedCoefficient = ( resiDiff + offset ) >> iTrShiftNeg;
  }

  // quantization
  const TCoeff iSign = (transformedCoefficient < 0 ? -1: 1);

  const Int quantisationCoefficient = enableScalingLists ? piQuantCoeff[uiPos] : defaultQuantisationCoefficient;

  const Int64 tmpLevel = (Int64)abs(transformedCoefficient) * quantisationCoefficient;

  const TCoeff quantisedCoefficient = (TCoeff((tmpLevel + iAdd ) >> iQBits)) * iSign;

  const TCoeff entropyCodingMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff entropyCodingMaximum =  (1 << maxLog2TrDynamicRange) - 1;
  pcCoeff[ uiPos ] = Clip3<TCoeff>( entropyCodingMinimum, entropyCodingMaximum, quantisedCoefficient );
}


Void TComTrQuant::invTrSkipDeQuantOneSample( TComTU &rTu, ComponentID compID, TCoeff inSample, Pel &reconSample, const QpParam &cQP, UInt uiPos )
{
        TComDataCU    *pcCU               = rTu.getCU();
  const UInt           uiAbsPartIdx       = rTu.GetAbsPartIdxTU();
  const TComRectangle &rect               = rTu.getRect(compID);
  const UInt           uiWidth            = rect.width;
  const UInt           uiHeight           = rect.height;
  const Int            QP_per             = cQP.per;
  const Int            QP_rem             = cQP.rem;
  const Int            maxLog2TrDynamicRange = pcCU->getSlice()->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID));
#if O0043_BEST_EFFORT_DECODING
  const Int            channelBitDepth    = pcCU->getSlice()->getSPS()->getStreamBitDepth(toChannelType(compID));
#else
  const Int            channelBitDepth    = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#endif
  const Int            iTransformShift    = getTransformShift(channelBitDepth, rTu.GetEquivalentLog2TrSize(compID), maxLog2TrDynamicRange);
  const Int            scalingListType    = getScalingListType(pcCU->getPredictionMode(uiAbsPartIdx), compID);
  const Bool           enableScalingLists = getUseScalingList(uiWidth, uiHeight, true);
#if !JVET_C0024_QTBT
  const UInt           uiLog2TrSize       = rTu.GetEquivalentLog2TrSize(compID);
#endif

  assert( scalingListType < SCALING_LIST_NUM );

  const Int rightShift = (IQUANT_SHIFT - (iTransformShift + QP_per)) + (enableScalingLists ? LOG2_SCALING_LIST_NEUTRAL_VALUE : 0);

  const TCoeff transformMinimum = -(1 << maxLog2TrDynamicRange);
  const TCoeff transformMaximum =  (1 << maxLog2TrDynamicRange) - 1;

  // Dequantisation

  TCoeff dequantisedSample;

  if(enableScalingLists)
  {
    const UInt             dequantCoefBits     = 1 + IQUANT_SHIFT + SCALING_LIST_BITS;
    const UInt             targetInputBitDepth = std::min<UInt>((maxLog2TrDynamicRange + 1), (((sizeof(Intermediate_Int) * 8) + rightShift) - dequantCoefBits));

    const Intermediate_Int inputMinimum        = -(1 << (targetInputBitDepth - 1));
    const Intermediate_Int inputMaximum        =  (1 << (targetInputBitDepth - 1)) - 1;

#if JVET_C0024_QTBT
    const UInt uiLog2TrWidth = g_aucConvertToBit[uiWidth] + MIN_CU_LOG2;
    const UInt uiLog2TrHeight = g_aucConvertToBit[uiHeight] + MIN_CU_LOG2;
    Int *piDequantCoef = getDequantCoeff(scalingListType,QP_rem,uiLog2TrWidth-1, uiLog2TrHeight-1);
#else
    Int *piDequantCoef = getDequantCoeff(scalingListType,QP_rem,uiLog2TrSize-2);
#endif

    if(rightShift > 0)
    {
      const Intermediate_Int iAdd      = 1 << (rightShift - 1);
      const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, inSample));
      const Intermediate_Int iCoeffQ   = ((Intermediate_Int(clipQCoef) * piDequantCoef[uiPos]) + iAdd ) >> rightShift;

      dequantisedSample = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
    }
    else
    {
      const Int              leftShift = -rightShift;
      const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, inSample));
      const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * piDequantCoef[uiPos]) << leftShift;

      dequantisedSample = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
    }
  }
  else
  {
    const Int scale     =  g_invQuantScales[QP_rem];
    const Int scaleBits =     (IQUANT_SHIFT + 1)   ;

    const UInt             targetInputBitDepth = std::min<UInt>((maxLog2TrDynamicRange + 1), (((sizeof(Intermediate_Int) * 8) + rightShift) - scaleBits));
    const Intermediate_Int inputMinimum        = -(1 << (targetInputBitDepth - 1));
    const Intermediate_Int inputMaximum        =  (1 << (targetInputBitDepth - 1)) - 1;

    if (rightShift > 0)
    {
      const Intermediate_Int iAdd      = 1 << (rightShift - 1);
      const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, inSample));
      const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * scale + iAdd) >> rightShift;

      dequantisedSample = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
    }
    else
    {
      const Int              leftShift = -rightShift;
      const TCoeff           clipQCoef = TCoeff(Clip3<Intermediate_Int>(inputMinimum, inputMaximum, inSample));
      const Intermediate_Int iCoeffQ   = (Intermediate_Int(clipQCoef) * scale) << leftShift;

      dequantisedSample = TCoeff(Clip3<Intermediate_Int>(transformMinimum,transformMaximum,iCoeffQ));
    }
  }

  // Inverse transform-skip

  if (iTransformShift >= 0)
  {
    const TCoeff offset = iTransformShift==0 ? 0 : (1 << (iTransformShift - 1));
    reconSample =  Pel(( dequantisedSample + offset ) >> iTransformShift);
  }
  else //for very high bit depths
  {
    const Int iTrShiftNeg = -iTransformShift;
    reconSample = Pel(dequantisedSample << iTrShiftNeg);
  }
}


Void TComTrQuant::crossComponentPrediction(       TComTU      & rTu,
                                            const ComponentID   compID,
                                            const Pel         * piResiL,
                                            const Pel         * piResiC,
                                                  Pel         * piResiT,
                                            const Int           width,
                                            const Int           height,
                                            const Int           strideL,
                                            const Int           strideC,
                                            const Int           strideT,
                                            const Bool          reverse )
{
  const Pel *pResiL = piResiL;
  const Pel *pResiC = piResiC;
        Pel *pResiT = piResiT;

  TComDataCU *pCU = rTu.getCU();
  const Int alpha = pCU->getCrossComponentPredictionAlpha( rTu.GetAbsPartIdxTU( compID ), compID );
  const Int diffBitDepth = pCU->getSlice()->getSPS()->getDifferentialLumaChromaBitDepth();

  for( Int y = 0; y < height; y++ )
  {
    if (reverse)
    {
      // A constraint is to be added to the HEVC Standard to limit the size of pResiL and pResiC at this point.
      // The likely form of the constraint is to either restrict the values to CoeffMin to CoeffMax,
      // or to be representable in a bitDepthY+4 or bitDepthC+4 signed integer.
      //  The result of the constraint is that for 8/10/12bit profiles, the input values
      //  can be represented within a 16-bit Pel-type.
#if RExt__HIGH_BIT_DEPTH_SUPPORT
      for( Int x = 0; x < width; x++ )
      {
        pResiT[x] = pResiC[x] + (( alpha * rightShift( pResiL[x], diffBitDepth) ) >> 3);
      }
#else
      const Int minPel=std::numeric_limits<Pel>::min();
      const Int maxPel=std::numeric_limits<Pel>::max();
      for( Int x = 0; x < width; x++ )
      {
        pResiT[x] = Clip3<Int>(minPel, maxPel, pResiC[x] + (( alpha * rightShift<Int>(Int(pResiL[x]), diffBitDepth) ) >> 3));
      }
#endif
    }
    else
    {
      // Forward does not need clipping. Pel type should always be big enough.
      for( Int x = 0; x < width; x++ )
      {
        pResiT[x] = pResiC[x] - (( alpha * rightShift<Int>(Int(pResiL[x]), diffBitDepth) ) >> 3);
      }
    }

    pResiL += strideL;
    pResiC += strideC;
    pResiT += strideT;
  }
}

#if PVC_AQS // core, decide Scaling factor
Void   TComTrQuant::setEstQSFactor(Pel* pPred, UInt uiStride, UInt uiWidth, UInt uiHeight, SliceType eSliceType, TComDataCU* pcCU, UInt uiAbsPartIdx, ComponentID eText)
{
  if (eText != COMPONENT_Y)
  {
    m_iEstQSFactor = 1 << ESM_SHIFT;
    return;
  }

  // first row and column
  if ((pcCU->getCUPelX() == 0 || pcCU->getCUPelY() == 0) && eSliceType == I_SLICE)
  {
    m_iEstQSFactor = 1 << ESM_SHIFT;
    return;
  }

  // get factor for luma residual coefficients only
  UInt sum;
  UInt temp, avg, st, num;
  Int  iTemp;
  Int  i, j;
  Int  step = 1;
  Pel* pPel = pPred;
  Int  up = (1 << ESM_SHIFT) + (1 << (ESM_SHIFT - 1)); // 1.5  * (1<<ESM_SHIFT)
  Int  low = (1 << ESM_SHIFT) - (1 << (ESM_SHIFT - 2)); // 0.75 * (1<<ESM_SHIFT)
  up = (Int) (sqrt(2.0) * (1 << ESM_SHIFT) + 0.5);
  low = (Int) (sqrt(0.5) * (1 << ESM_SHIFT) + 0.5);

  Int i_start, j_start;
  Int i_curr, j_curr;
  Int blk_w;
  Int cu_w;
  TComPicYuv* pcPic = pcCU->getPic()->getPicYuvRec();
  Pel* pPelPic = pcPic->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);
  UInt uiPicStride = pcPic->getStride(COMPONENT_Y);
  Int TBStride = 16;
  Pel *pSB1, *pSB2;

  Int BG_SIZE = 16;

  //--------------- modify Upper bound for high bitrate case (QP<28) - in AQS
  Int iQP = pcCU->getSlice()->getSliceQp();
  Int iNormFactor = 0;
  Int iFactor = 0;

  Int iOffsetQPIslice = 28;
  Double dScalerIsliceStep = 0.10;
  Double dScalerIsliceMax = 1.00;
  Double dScalerIsliceMin = 0.25;

  Int iOffsetQPPBslice = 28;
  Double dScalerPBsliceStep = 0.10;
  Double dScalerPBsliceMax = 1.00;
  Double dScalerPBsliceMin = 0.25;

  Double dMaxScaleFactor = 0;

  Int iBitShfit = pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) - 8;

  if (pcCU->getSlice()->getSliceType() == I_SLICE)
  {
    dMaxScaleFactor = Clip3(dScalerIsliceMin, dScalerIsliceMax, (iQP - iOffsetQPIslice) * dScalerIsliceStep + 1.0) * (1.414 - 1.0) + 1.0;
    up = (Int) ((1 << ESM_SHIFT) * dMaxScaleFactor + 0.5);
  }
  else
  {
    dMaxScaleFactor = Clip3(dScalerPBsliceMin, dScalerPBsliceMax, (iQP - iOffsetQPPBslice) * dScalerPBsliceStep + 1.0) * (1.414 - 1.0) + 1.0;
    up = (int) ((1 << ESM_SHIFT) * dMaxScaleFactor + 0.5);
  }
  //----------------

  if (uiWidth <= BG_SIZE)
  {
    // get avg and var
    i_start = pcCU->getCUPelX();
    j_start = pcCU->getCUPelY();
    i_curr = i_start + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
    j_curr = j_start + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
    cu_w = pcCU->getWidth(0);
    if (uiWidth == BG_SIZE)
    {
      pPel = pPred;
      blk_w = uiWidth;
    }
    else // if TB is smaller than a 16x16 blk;
    {
      if (pcCU->getPredictionMode(0) == MODE_INTRA)
      {
        i_start = i_curr + uiWidth - BG_SIZE;
        j_start = j_curr + uiHeight - BG_SIZE;
        if (i_start < 0 || j_start < 0)
        {
          pPel = pPred;
          blk_w = uiWidth;
        }
        else
        {
          Pel apTB[256];
          pSB1 = pPelPic - (j_curr - j_start)*uiPicStride - (i_curr - i_start);
          pSB2 = pPred;
          composite_block(apTB, TBStride, pSB1, BG_SIZE, BG_SIZE, uiPicStride, pSB2, uiWidth, uiHeight, uiStride, BG_SIZE - uiWidth, BG_SIZE - uiHeight);
          pPel = apTB;
          blk_w = BG_SIZE;
          uiStride = TBStride;
        }
      }
      else //Inter Slices
      {
        if (uiStride < BG_SIZE) // cu_width < 16, 
        {
          i_start = i_curr + uiWidth - BG_SIZE;
          j_start = j_curr + uiHeight - BG_SIZE;
          if (i_start < 0 || j_start < 0) // NOT all the 16x16 blk is available, then use the current CU prediction sample
          {
            i_start = pcCU->getCUPelX();
            j_start = pcCU->getCUPelY();
            pPel = pPred - (j_curr - j_start) * uiStride - (i_curr - i_start);
            blk_w = uiStride;
          }
          else
          {
            Pel apTB[256];

            pSB1 = pPelPic - (j_curr - j_start)*uiPicStride - (i_curr - i_start); //reconstructed sample outside the current CU
            i_start = pcCU->getCUPelX();
            j_start = pcCU->getCUPelY();
            pSB2 = pPred - (j_curr - j_start)*uiStride - (i_curr - i_start); //prediction sample of the current CU
            composite_block(apTB, TBStride, pSB1, BG_SIZE, BG_SIZE, uiPicStride, pSB2, uiStride, uiStride, uiStride, BG_SIZE - uiStride, BG_SIZE - uiStride);
            pPel = apTB;
            blk_w = BG_SIZE;
            uiStride = TBStride;
          }
        }
        else
        {
          i_start = (i_curr / BG_SIZE) * BG_SIZE;
          j_start = (j_curr / BG_SIZE) * BG_SIZE; //find the 16x16 block that the TB belongs to
          pPel = pPred - (j_curr - j_start) * uiStride - (i_curr - i_start);
          blk_w = BG_SIZE;
        }
      }
    }

    Pel* pPel2 = pPel;

    sum = 0;
    num = 0;
    for (j = 0; j<blk_w; j += step)
    {
      for (i = 0; i<blk_w; i += step)
      {
        temp = *(pPel + i);
        sum += temp;
        num++;
      }
      pPel += uiStride * step;
    }
    avg = (sum + (num >> 1)) / num;

    sum = 0;
    num = 0;
    pPel = pPel2;
    for (j = 0; j<blk_w; j += step)
    {
      for (i = 0; i<blk_w; i += step)
      {
        temp = *(pPel + i);
        sum += abs((int) temp - (int) avg);
        num++;
      }
      pPel += uiStride * step;
    }
    st = (sum + (num >> 1)) / num;

    //get ESM value
    //factor = sqrt(g_aiFactorLuminance[avg] * g_aiFactorContrast[st]/normalizer);
    //simplified version - sqrt(x) is approximated by (x-1)/2 + 1
    iTemp = (g_aiFactorLuminance[avg >> iBitShfit] * g_aiFactorContrast[st >> iBitShfit] + (1 << (ESM_SHIFT - 2))) >> (ESM_SHIFT - 1); //normalizer here is 0.5
    iTemp = (iTemp * g_ESMFactor + 128) >> 8;
    iTemp = (iTemp + (1 << ESM_SHIFT) + 1) >> 1;
    m_iEstQSFactor = Clip3(low, up, iTemp);
  }
  else
  {
    int m, n;
    m_iEstQSFactor = 0;

    for (n = 0; n<uiHeight; n += BG_SIZE)
    {
      for (m = 0; m<uiWidth; m += BG_SIZE)
      {
        //get avg and var
        sum = 0;
        num = 0;
        pPel = pPred + uiStride * n;

        Pel* pPel2 = pPel;

        for (j = n; j<n + BG_SIZE; j += step)
        {
          for (i = m; i<m + BG_SIZE; i += step)
          {
            temp = *(pPel + i);
            sum += temp;
            num++;
          }
          pPel += uiStride * step;
        }
        avg = (sum + (num >> 1)) / num;

        sum = 0;
        num = 0;
        pPel = pPel2;
        for (j = n; j<n + BG_SIZE; j += step)
        {
          for (i = m; i<m + BG_SIZE; i += step)
          {
            temp = *(pPel + i);
            sum += abs((int) temp - (int) avg);
            num++;
          }
          pPel += uiStride * step;
        }
        st = (sum + (num >> 1)) / num;

        //get ESM value
        ///factor = sqrt(g_aiFactorLuminance[avg] * g_aiFactorContrast[st]/normalizer);
        //simplified version - sqrt(x) is approximated by (x-1)/2 + 1
        iTemp = (g_aiFactorLuminance[avg >> iBitShfit] * g_aiFactorContrast[st >> iBitShfit] + (1 << (ESM_SHIFT - 2))) >> (ESM_SHIFT - 1); //normalizer here is 0.5
        iTemp = (iTemp * g_ESMFactor + 128) >> 8;
        iTemp = (iTemp + (1 << ESM_SHIFT) + 1) >> 1;
        iTemp = Clip3(low, up, iTemp);
        m_iEstQSFactor += iTemp;
      }
    }
    num = (uiWidth / BG_SIZE)*(uiHeight / BG_SIZE);
    m_iEstQSFactor = (m_iEstQSFactor + (num >> 1)) / num;
  }
}
#endif



#if VCEG_AZ08_KLT_COMMON
inline int MY_INT(double x)
{
    return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5);
}

TempLibFast::TempLibFast()
{
    m_pX = NULL;
    m_pY = NULL;
    m_pDiff = NULL;
    m_pId = NULL;
    m_pXInteger = NULL;
    m_pYInteger = NULL;
    m_pDiffInteger = NULL;
    m_pIdInteger = NULL;
}

TempLibFast::~TempLibFast()
{
    if (m_pX != NULL)
    {
        delete[]m_pX;
        m_pX = NULL;
    }

    if (m_pY != NULL)
    {
        delete[]m_pY;
        m_pY = NULL;
    }


    if (m_pXInteger != NULL)
    {
        delete[]m_pXInteger;
        m_pXInteger = NULL;
    }
    if (m_pYInteger != NULL)
    {
        delete[]m_pYInteger;
        m_pYInteger = NULL;
    }
    if (m_pDiffInteger != NULL)
    {
        delete[]m_pDiffInteger;
        m_pDiffInteger = NULL;
    }
    if (m_pIdInteger != NULL)
    {
        delete[]m_pIdInteger;
        m_pIdInteger = NULL;
    }

    if (m_pDiff != NULL)
    {
        delete[]m_pDiff;
        m_pDiff = NULL;
    }

    if (m_pId != NULL)
    {
        delete[]m_pId;
        m_pId = NULL;
    }
}

Void TempLibFast::init(UInt iSize)
{
    m_iSize = iSize;
    if( m_pY == NULL )
    {
      m_pX = new Int[iSize];
      m_pY = new Int[iSize];
      m_pDiff = new DistType[iSize];
      m_pId = new Short[iSize];

      m_pXInteger = new Int[iSize];
      m_pYInteger = new Int[iSize];
      m_pDiffInteger = new DistType[iSize];
      m_pIdInteger = new Short[iSize];
    }
}

Void TempLibFast::initDiff(UInt uiPatchSize, Int bitDepth, Int iCandiNumber)
{
#if VCEG_AZ08_USE_SAD_DISTANCE
    DistType maxValue = ((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*(uiPatchSize*uiPatchSize);
#endif
#if VCEG_AZ08_USE_SSD_DISTANCE
    DistType maxValue = ((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*uiPatchSize*uiPatchSize;
#endif
    g_maxValueThre = maxValue;
    m_diffMax = maxValue;
    for (Int i = 0; i < iCandiNumber; i++)
    {
        m_pDiff[i] = maxValue;
    }
}


Void insertNode(DistType diff, Int iXOffset, Int iYOffset, DistType *pDiff, Int *pX, Int *pY, Short *pId, UInt uiLibSizeMinusOne, Int setId)
{
    Int j; //should be Int rather than UInt
    //bin-search insert
    Int iStart = 0, iEnd = uiLibSizeMinusOne, iMiddle = (iStart + iEnd) >> 1;
    while (iStart < iEnd - 1)
    {
        if (pDiff[iMiddle] > diff)
        {
            iEnd = iMiddle;
        }
        else if (pDiff[iMiddle] < diff)
        {
            iStart = iMiddle;
        }
        else
        {
            iStart = iMiddle;
            iEnd = iMiddle;
            break;
        }
        iMiddle = (iStart + iEnd) >> 1;
    }
    if (pDiff[iStart] <= diff)
    {
        j = iStart;
    }
    else
    {
        j = iStart - 1;
    }

    if (!(j >= 0 && diff == pDiff[j]))
    {
        Int iPlacedPos = j + 1;
        //Insert the new node: memory copy and assign values to that node
        if (iPlacedPos < uiLibSizeMinusOne)
        {
            Int copyNum = (uiLibSizeMinusOne - iPlacedPos);
            Int iPlacePosAddOne = iPlacedPos + 1;
            memmove(&pDiff[iPlacePosAddOne], &pDiff[iPlacedPos], copyNum*sizeof(DistType));
            memmove(&pX[iPlacePosAddOne], &pX[iPlacedPos], copyNum*sizeof(Int));
            memmove(&pY[iPlacePosAddOne], &pY[iPlacedPos], copyNum*sizeof(Int));
            memmove(&pId[iPlacePosAddOne], &pId[iPlacedPos], copyNum*sizeof(Short));
        }
        pDiff[iPlacedPos] = diff;
        pX[iPlacedPos] = iXOffset;
        pY[iPlacedPos] = iYOffset;
        pId[iPlacedPos] = setId;
    }
}

/** NxN forward KL-transform (1D) using brute force matrix multiplication
*  \param block pointer to input data (residual)
*  \param coeff pointer to output data (transform coefficients)
*  \param uiTrSize transform size (uiTrSize x uiTrSize)
*/
void xKLTr(Int bitDepth, TCoeff *block, TCoeff *coeff, UInt uiTrSize)  //void xKLTr(Int bitDepth, Pel *block, Short *coeff, UInt uiTrSize)  
{
    Int i, k, iSum;
    Int uiDim = uiTrSize*uiTrSize;
    UInt uiLog2TrSize = g_aucConvertToBit[uiTrSize] + 2;
    Int shift = bitDepth + 2 * uiLog2TrSize + KLTBASIS_SHIFTBIT - 15;
    Int add = 1 << (shift - 1);
    UInt uiTarDepth = g_aucConvertToBit[uiTrSize];
    Short **pTMat = g_ppsEigenVector[uiTarDepth];
    for (i = 0; i< uiDim; i++)
    {
        iSum = 0;
        Short *pT = pTMat[i]; //g_psEigenVectorDim16[i];
        for (k = 0; k<uiDim; k++)
        {
            iSum += pT[k] * block[k];
        }
        coeff[i] = (iSum + add) >> shift;
    }
}

/** NxN inverse KL-transform (1D) using brute force matrix multiplication
*  \param coeff pointer to input data (transform coefficients)
*  \param block pointer to output data (residual)
*  \param uiTrSize transform size (uiTrSize x uiTrSize)
*/
void xIKLTr(Int bitDepth, TCoeff *coeff, TCoeff *block, UInt uiTrSize)  //void xIKLTr(Short *coeff, Pel *block, UInt uiTrSize) 
{
    Int i, k, iSum;
    UInt uiDim = uiTrSize*uiTrSize;
    //Int shift = 7 + KLTBASIS_SHIFTBIT - (g_bitDepthY - 8);
    //ComponentID component = COMPONENT_Y;
    //const Int channelBitDepth = rTu.getCU()->getSlice()->getSPS()->getBitDepth(toChannelType(component));
    Int shift = 7 + KLTBASIS_SHIFTBIT - (bitDepth - 8);
    Int add = 1 << (shift - 1);
    UInt uiTarDepth = g_aucConvertToBit[uiTrSize];
    Short **pTMat = g_ppsEigenVector[uiTarDepth];
    for (i = 0; i < uiDim; i++)
    {
        iSum = 0;
        for (k = 0; k < uiDim; k++)
        {
            iSum += pTMat[k][i] * coeff[k];
        }
        block[i] = (iSum + add) >> shift;
    }
}

Bool TComTrQuant::deriveKLT(UInt uiBlkSize, UInt uiUseCandiNumber)
{
    Bool bSucceedFlag = true;
#if VCEG_AZ08_FAST_DERIVE_KLT
    Int fast_klt_blksize = 8;
    if (FAST_KLT_CANDINUM > 64)
    {
        fast_klt_blksize = 16;
    }
    else if (FAST_KLT_CANDINUM > 256)
    {
        fast_klt_blksize = 32;
    }
    else if (FAST_KLT_CANDINUM > 1024)
    {
        fast_klt_blksize = 64;
    }

    if (uiBlkSize >= fast_klt_blksize)
    {
        bSucceedFlag = derive1DimKLT_Fast(uiBlkSize, uiUseCandiNumber);
    }
    else
    {
        bSucceedFlag = derive1DimKLT(uiBlkSize, uiUseCandiNumber);
    }
#else
    bSucceedFlag = derive1DimKLT(uiBlkSize, iUseCandiNumber);
#endif
    return bSucceedFlag;
}


typedef Double _Ty;
void OrderData(int *IdRecord, _Ty *pData, UInt uiSize)
{
    int k;
    int i, j;
    _Ty temp;
    int tempId;
    assert(pData);
    for (k = 0; k < uiSize; k++)
    {
        IdRecord[k] = k;
        if (pData[k] < 0)
        {
            pData[k] = -pData[k]; // original values changed here
        }
    }
    for (i = 1; i < uiSize; i++)
    {
        for (j = uiSize - 1; j >= i; j--)
        {
            if (pData[j] > pData[j - 1])
            {
                temp = pData[j];
                pData[j] = pData[j - 1];
                pData[j - 1] = temp;

                tempId = IdRecord[j];
                IdRecord[j] = IdRecord[j - 1];
                IdRecord[j - 1] = tempId;
            }
        }
    }
}

Bool TComTrQuant::derive1DimKLT(UInt uiBlkSize, UInt uiUseCandiNumber)
{
    Bool bSucceed = true;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    //calculate covariance matrix
    const UInt uiDim = uiBlkSize*uiBlkSize;
    covMatrixType *covMatrix = m_pCovMatrix[uiTarDepth];
    calcCovMatrix(m_pData, uiUseCandiNumber, covMatrix, uiDim);
    EigenType **pdEigenVector = m_pppdEigenVector[uiTarDepth];
    Short **psEigenVector = m_pppsEigenVector[uiTarDepth];
    g_ppsEigenVector[uiTarDepth] = psEigenVector;

    matrixTypeDefined Cov(uiDim, uiDim);
    UInt i = 0;
    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        for (UInt uiCol = 0; uiCol < uiDim; uiCol++)
        {
            Cov(uiRow, uiCol) = covMatrix[i++];
        }
    }
    SelfAdjointEigenSolver<matrixTypeDefined> es(Cov);
    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        m_pEigenValues[uiRow] = es.eigenvalues()[uiRow];
    }
    OrderData(m_pIDTmp, m_pEigenValues, uiDim);
    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        UInt uiMapedCol = m_pIDTmp[uiRow];
        vectorType v = es.eigenvectors().col(uiMapedCol);
        for (UInt uiCol = 0; uiCol < uiDim; uiCol++)
        {
            pdEigenVector[uiRow][uiCol] = (EigenType)v(uiCol);
        }
    }
    Int scale = uiBlkSize*(1 << KLTBASIS_SHIFTBIT);
#if VCEG_AZ08_USE_SSE_SCLAE
    scaleMatrix(pdEigenVector, psEigenVector, scale, uiDim, uiDim);
#else
    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        for (UInt uiCol = 0; uiCol < uiDim; uiCol++)
        {
            psEigenVector[uiRow][uiCol] = MY_INT(pdEigenVector[uiRow][uiCol] * scale);
        }
    }
#endif
    return bSucceed;
}

Bool TComTrQuant::derive1DimKLT_Fast(UInt uiBlkSize, UInt uiUseCandiNumber)
{
    Bool bSucceed = true;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    const UInt uiDim = uiBlkSize*uiBlkSize;
    covMatrixType *covMatrix = m_pCovMatrix[uiTarDepth];
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    UInt uiSampleNum = min(uiUseCandiNumber, uiTargetCandiNum);
    //calculate covariance matrix
    calcCovMatrixXXt(m_pData, uiSampleNum, covMatrix, uiDim);
    EigenType **pdEigenVector = m_pppdTmpEigenVector;
    EigenType **pdEigenVectorTarget = m_pppdEigenVector[uiTarDepth];
    Short **psEigenVector = m_pppsEigenVector[uiTarDepth];
    g_ppsEigenVector[uiTarDepth] = psEigenVector;

    //depend on eigen libarary 
    matrixTypeDefined Cov(uiSampleNum, uiSampleNum);
    UInt i = 0;
    for (UInt uiRow = 0; uiRow < uiSampleNum; uiRow++)
    {
        for (UInt uiCol = 0; uiCol < uiSampleNum; uiCol++)
        {
            Cov(uiRow, uiCol) = covMatrix[i++];
        }
    }
    SelfAdjointEigenSolver<matrixTypeDefined> es(Cov);
    for (UInt uiRow = 0; uiRow < uiSampleNum; uiRow++)
    {
        m_pEigenValues[uiRow] = es.eigenvalues()[uiRow];
    }
    OrderData(m_pIDTmp, m_pEigenValues, uiSampleNum);
    EigenType *pEigenVectorRow;
    for (UInt uiRow = 0; uiRow < uiSampleNum; uiRow++)
    {
        UInt uiMapedCol = m_pIDTmp[uiRow];
        vectorType v = es.eigenvectors().col(uiMapedCol);
        pEigenVectorRow = pdEigenVector[uiRow];
        for (UInt uiCol = 0; uiCol < uiSampleNum; uiCol++)
        {
            pEigenVectorRow[uiCol] = (EigenType)v(uiCol);
        }
    }

    UInt uiCalcEigNum = uiSampleNum;
#if VCEG_AZ08_FORCE_USE_GIVENNUM_BASIS
    uiCalcEigNum = min(uiSampleNum, (UInt)FORCE_BASIS_NUM);
#endif

    UInt uiFilledDim = uiCalcEigNum;
    Double dValue;

    assert(uiSampleNum <= uiDim);
    double dIgnoreThreshouldFinal = IGNORE_THRESHOULD_OF_LARGEST;

    for (UInt uiRow = 0; uiRow < uiCalcEigNum; uiRow++)
    {
        if (m_pEigenValues[uiRow] < dIgnoreThreshouldFinal)
        {
            uiFilledDim = uiRow;
            break;
        }
        EigenType *pdEigenVectorRow = pdEigenVector[uiRow];
        EigenType *pThisBasisRow = pdEigenVectorTarget[uiRow];
        Double dValueNorm = sqrt(m_pEigenValues[uiRow]);
        for (UInt uiCol = 0; uiCol < uiDim; uiCol++)
        {
#if VCEG_AZ08_USE_FLOATXSHORT_SSE
            dValue = InnerProduct_SSE_FLOATXSHORT(pdEigenVectorRow, m_pDataT[uiCol], uiSampleNum);
#else
            dValue = 0;
            for (UInt k = 0; k < uiSampleNum; k++)
            {
                dValue += pdEigenVectorRow[k] * m_pData[k][uiCol];
            }
#endif
            pThisBasisRow[uiCol] = (EigenType)(dValue / dValueNorm);
        }
    }

    UInt uiBasisNum = uiFilledDim;
    Int iScale = uiBlkSize*(1 << KLTBASIS_SHIFTBIT);
#if VCEG_AZ08_USE_SSE_SCLAE
    scaleMatrix(pdEigenVectorTarget, psEigenVector, iScale, uiBasisNum, uiDim);
#else
    for (UInt uiRow = 0; uiRow < uiBasisNum; uiRow++)
    {
        for (UInt uiCol = 0; uiCol < uiDim; uiCol++)
        {
            psEigenVector[uiRow][uiCol] = MY_INT(pdEigenVectorTarget[uiRow][uiCol] * iScale);
        }
    }
#endif

    for (int row = uiBasisNum; row < uiDim; row++)
    {
        memset(psEigenVector[row], 0, sizeof(Short)*uiDim);
    }
    return bSucceed;
}


#if VCEG_AZ08_FAST_DERIVE_KLT
//calculate XX' ranther than X'X  X:N * Dim
Void TComTrQuant::calcCovMatrixXXt(TrainDataType **pData, UInt uiSampleNum, covMatrixType *pCovMatrix, UInt uiDim)
{
    //Get the covariance matrix
    Int covValue; //should be int; if float, the accuracy will be low.
    TrainDataType *pDataCol;
    for (UInt uiRow = 0; uiRow < uiSampleNum; uiRow++)
    {
        TrainDataType *pDataRow = pData[uiRow];
        Int offset = uiRow*uiSampleNum;
        for (UInt uiCol = 0; uiCol <= uiRow; uiCol++)
        {
            pDataRow = pData[uiRow];
            pDataCol = pData[uiCol];
#if VCEG_AZ08_USE_SHORTXSHORT_SSE
            covValue = InnerProduct_SSE_SHORT(pDataRow, pDataCol, uiDim);
#else
            covValue = 0;
            for (Int i = 0; i < uiDim; i++)
            {
                covValue += (*(pDataRow++)) * (*(pDataCol++));
            }
#endif
            pCovMatrix[offset + uiCol] = (covMatrixType)covValue;
        }
    }

    for (UInt uiRow = 0; uiRow < uiSampleNum; uiRow++)
    {
        for (UInt uiCol = uiRow + 1; uiCol < uiSampleNum; uiCol++)
        {
            pCovMatrix[uiRow*uiSampleNum + uiCol] = pCovMatrix[uiCol*uiSampleNum + uiRow];
        }
    }
}
#endif

Void TComTrQuant::calcCovMatrix(TrainDataType **pData, UInt uiSampleNum, covMatrixType *pCovMatrix, UInt uiDim)
{
    //Get the covariance matrix
    covMatrixType dCovValue;
    Int covValue;
#if !(VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY && VCEG_AZ08_USE_SHORTXSHORT_SSE)
    UInt i;
    TrainDataType *pSample;
#endif
    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        for (UInt uiCol = 0; uiCol <= uiRow; uiCol++)
        {
#if VCEG_AZ08_USE_SHORTXSHORT_SSE
            covValue = InnerProduct_SSE_SHORT(m_pDataT[uiRow], m_pDataT[uiCol], uiSampleNum);
#else
            covValue = 0;
            for (i = 0; i < uiSampleNum; i++)
            {
                pSample = pData[i];
                covValue += pSample[uiRow] * pSample[uiCol];
            }
#endif
            dCovValue = (covMatrixType)covValue / (covMatrixType)(uiSampleNum);
            pCovMatrix[uiRow*uiDim + uiCol] = dCovValue;
        }
    }

    for (UInt uiRow = 0; uiRow < uiDim; uiRow++)
    {
        for (UInt uiCol = uiRow + 1; uiCol < uiDim; uiCol++)
        {
            pCovMatrix[uiRow*uiDim + uiCol] = pCovMatrix[uiCol*uiDim + uiRow];
        }
    }
}
#endif

#if VCEG_AZ08_INTER_KLT 
Void TComTrQuant::getTargetPatch(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absTUPartIdx, TComYuv* pcPred, UInt uiBlkSize, UInt uiTempSize)
{
    ComponentID compID = COMPONENT_Y;
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    Pel **tarPatch = m_pppTarPatch[uiTarDepth];
    UInt  uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
    Pel   *pCurrStart = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
    UInt  uiPicStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);

    Pel   *pPred = pcPred->getAddr(compID, absTUPartIdx);
    UInt  uiPredStride = pcPred->getStride(compID);

    //fill the target block
    UInt uiY, uiX;
    for (uiY = uiTempSize; uiY < uiPatchSize; uiY++)
    {
        for (uiX = uiTempSize; uiX < uiPatchSize; uiX++)
        {
            tarPatch[uiY][uiX] = pPred[uiX - uiTempSize];
        }
        pPred += uiPredStride;
    }

    //fill template
    //up-left & up 
    Pel  *tarTemp;
    Pel  *pCurrTemp = pCurrStart - uiTempSize*uiPicStride - uiTempSize;
    for (uiY = 0; uiY < uiTempSize; uiY++)
    {
        tarTemp = tarPatch[uiY];
        for (uiX = 0; uiX < uiPatchSize; uiX++)
        {
            tarTemp[uiX] = pCurrTemp[uiX];
        }
        pCurrTemp += uiPicStride;
    }
    //left
    for (uiY = uiTempSize; uiY < uiPatchSize; uiY++)
    {
        tarTemp = tarPatch[uiY];
        for (uiX = 0; uiX < uiTempSize; uiX++)
        {
            tarTemp[uiX] = pCurrTemp[uiX];
        }
        pCurrTemp += uiPicStride;
    }
}


Void TComTrQuant::candidateSearch(TComDataCU *pcCU, UInt uiPartAddr, UInt uiBlkSize, UInt uiTempSize)
{
    ComponentID compID = COMPONENT_Y;
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    Pel **tarPatch = m_pppTarPatch[uiTarDepth];
    Int      iRefIdx[2] = { -1, -1 };
#if !VCEG_AZ08_INTER_KLT_MV_BUGFIXED
    TComMv      cMvs[2];
#endif
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
    //Initialize the library for saving the best candidates
    m_tempLibFast.initDiff(uiPatchSize, channelBitDepth, uiTargetCandiNum);
    Short setId = -1; //to record the reference picture.

    Int iCurrPOC = pcCU->getPic()->getPOC();
    Int iRefPOC = 0;
    TComMv cMvRef;

    Int iRefPOCs[2] = { 0, 0 };
    TComMv      cMvRefs[2];
    cMvRefs[0].setZero();
    cMvRefs[1].setZero();
    Int filledNum = 0;
    for (Int iRefList = 0; iRefList < 2; iRefList++)
    {
        RefPicList eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
        iRefIdx[eRefPicList] = pcCU->getCUMvField(eRefPicList)->getRefIdx(uiPartAddr);
        if (iRefIdx[eRefPicList] >= 0)
        {
#if VCEG_AZ08_INTER_KLT_MV_BUGFIXED
            cMvRefs[filledNum] = pcCU->getCUMvField(eRefPicList)->getMv(uiPartAddr);
#else
            cMvRefs[filledNum] = cMvs[eRefPicList];
#endif
            iRefPOCs[filledNum] = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdx[eRefPicList])->getPOC();
            filledNum++;
        }
    }

    //for other reference frames not related with prediction
    for (Int iRefList = 0; iRefList < 2; iRefList++)
    {
        RefPicList eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
        for (Int iRefIdxTemp = 0; iRefIdxTemp < pcCU->getSlice()->getNumRefIdx(eRefPicList); iRefIdxTemp++)
        {
            TComMv   cMv;
            Int iTargetPOC = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxTemp)->getPOC();
            Int minDistancePOC = MAX_INT;
            Int distance;
            //use the reference frame (having mvs) with nearest distacne from the current iTargetPOC 
            for (Int filledId = 0; filledId < filledNum; filledId++)
            {
                distance = abs(iRefPOCs[filledId] - iTargetPOC);
                if (distance < minDistancePOC)
                {
                    minDistancePOC = distance;
                    iRefPOC = iRefPOCs[filledId];
                    cMvRef = cMvRefs[filledId];
                }
            }

            Double dScale = ((Double)(iCurrPOC - iTargetPOC) / (Double)(iCurrPOC - iRefPOC));
            cMv.set((Short)(cMvRef.getHor()*dScale + 0.5), (Short)(cMvRef.getVer()*dScale + 0.5));
            if (eRefPicList == REF_PIC_LIST_1)
            {
                //check whether the reference frames in list1 had appeared in list0, if so, will not use it to find candidates.
                bool bSimilarSearchRegion = false;
                for (Int iRefIdxList0 = 0; iRefIdxList0 < pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_0); iRefIdxList0++)
                {
                    UInt prevPicPOCList0 = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxList0)->getPOC();
                    if (prevPicPOCList0 == iTargetPOC)
                    {
                        bSimilarSearchRegion = true;
                        break; //had appeared in list0 for this reference in list1
                    }
                }
                if (bSimilarSearchRegion)
                {
                    continue;
                }
            }

            TComPic* refPic = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxTemp);
            UInt uiRow = 0;
            UInt uiCol = 0;
            TComPicYuv *refPicRec = refPic->getPicQuaYuvRec(uiRow, uiCol);
            Pel *ref = refPicRec->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);
            setId++;
            setRefPicUsed(setId, ref); //to facilitate the access of each candidate point 
            setRefPicBuf(setId, refPic);
            setStride(refPicRec->getStride(compID));
            searchCandidateFromOnePicInteger(pcCU, uiPartAddr, refPicRec, cMv, tarPatch, uiPatchSize, uiTempSize, setId);
        }
    }

    Short setIdFraStart = setId + 1;
    RecordPosition(uiTargetCandiNum);
    searchCandidateFraBasedOnInteger(pcCU, tarPatch, uiPatchSize, uiTempSize, uiPartAddr, setIdFraStart);
}

Void TComTrQuant::searchCandidateFromOnePicInteger(TComDataCU *pcCU, UInt uiPartAddr, TComPicYuv *refPic, TComMv  cMv, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt setId)
{
    const ComponentID compID = COMPONENT_Y;
    UInt uiBlkSize = uiPatchSize - uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    UInt  uiLibSizeMinusOne = uiTargetCandiNum - 1;
    m_uiPartLibSize = uiTargetCandiNum;

    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    DistType *pDiff = m_tempLibFast.getDiff();
    Short *pId = m_tempLibFast.getId();

    Int  refStride = refPic->getStride(compID);
    Pel  *ref = refPic->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);

    Int      iSrchRng = SEARCHRANGE;
    TComMv  cMvSrchRngLT;
    TComMv  cMvSrchRngRB;
    xSetSearchRange(pcCU, cMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);

    Int mvYMin = cMvSrchRngLT.getVer();
    Int mvYMax = cMvSrchRngRB.getVer();

    Int mvXMin = cMvSrchRngLT.getHor();
    Int mvXMax = cMvSrchRngRB.getHor();

    //search
    Pel *refMove = ref + mvYMin*refStride + mvXMin;
    DistType *pDiffEnd = &pDiff[uiLibSizeMinusOne];

    DistType diff;
    for (Int iYOffset = mvYMin; iYOffset <= mvYMax; iYOffset++)
    {
        Pel *refCurr = refMove;
        for (Int iXOffset = mvXMin; iXOffset <= mvXMax; iXOffset++)
        {
            //The position of the leftup pixel within this block: refCurr = ref + iYOffset*refStride + iXOffset;
            diff = calcPatchDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
            refCurr++;
            if (diff > 0 && diff < (*pDiffEnd))//when residual is zero, may not contribute to the distribution.
            {
                insertNode(diff, iXOffset, iYOffset, pDiff, pX, pY, pId, uiLibSizeMinusOne, setId);
            }
        }
        refMove += refStride;
    }
}


Void TComTrQuant::searchCandidateFraBasedOnInteger(TComDataCU *pcCU, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt uiPartAddr, Short setIdFraStart)
{
    const ComponentID compID = COMPONENT_Y;
    UInt uiBlkSize = uiPatchSize - uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    UInt  uiLibSizeMinusOne = uiTargetCandiNum - 1;
    Int  refStride = getStride();
    Pel *ref;
    UInt setId;
    Int iCandiPosNum = m_uiPartLibSize;

    Int k;
    Int *pXInteger = m_tempLibFast.getXInteger();
    Int *pYInteger = m_tempLibFast.getYInteger();
    Short *pIdInteger = m_tempLibFast.getIdInteger();

    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    DistType *pDiff = m_tempLibFast.getDiff();
    Short *pId = m_tempLibFast.getId();

    DistType diff;
    Pel *refCurr, *refCenter;
    DistType *pDiffEnd = &pDiff[uiLibSizeMinusOne];
    Int iOffsetY, iOffsetX;
    TComPic* refPic;
    UInt uiIdxAddr = pcCU->getZorderIdxInCtu() + uiPartAddr;
    Short setIdFra = setIdFraStart - 1;
    for (k = 0; k < iCandiPosNum; k++)
    {
        setId = pIdInteger[k];
        refPic = getRefPicBuf(setId);

        iOffsetY = pYInteger[k];
        iOffsetX = pXInteger[k];

        for (UInt uiRow = 0; uiRow < 4; uiRow++)
        {
            for (UInt uiCol = 0; uiCol < 4; uiCol++)
            {
                if (uiRow != 0 || uiCol != 0)
                {
                    TComPicYuv *refPicRec = refPic->getPicQuaYuvRec(uiRow, uiCol);
                    ref = refPicRec->getAddr(compID, pcCU->getCtuRsAddr(), uiIdxAddr);

                    setIdFra++;
                    setRefPicUsed(setIdFra, ref);
                    refCenter = ref + iOffsetY*refStride + iOffsetX;
                    //center
                    refCurr = refCenter;
                    diff = calcPatchDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                    if (diff > 0 && diff < (*pDiffEnd))
                    {
                        insertNode(diff, iOffsetX, iOffsetY, pDiff, pX, pY, pId, uiLibSizeMinusOne, setIdFra);
                    }
                    //left
                    refCurr = refCenter - 1;
                    diff = calcPatchDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                    if (diff > 0 && diff < (*pDiffEnd))
                    {
                        insertNode(diff, iOffsetX - 1, iOffsetY, pDiff, pX, pY, pId, uiLibSizeMinusOne, setIdFra);
                    }
                    //up
                    refCurr = refCenter - refStride;
                    diff = calcPatchDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                    if (diff > 0 && diff < (*pDiffEnd))
                    {
                        insertNode(diff, iOffsetX, iOffsetY - 1, pDiff, pX, pY, pId, uiLibSizeMinusOne, setIdFra);
                    }
                    //left-up
                    refCurr = refCenter - refStride - 1;
                    diff = calcPatchDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                    if (diff > 0 && diff < (*pDiffEnd))
                    {
                        insertNode(diff, iOffsetX - 1, iOffsetY - 1, pDiff, pX, pY, pId, uiLibSizeMinusOne, setIdFra);
                    }
                }
            }
        }
    }
}

Void TComTrQuant::xSetSearchRange(TComDataCU* pcCU, TComMv& cMvPred, Int iSrchRng, TComMv& rcMvSrchRngLT, TComMv& rcMvSrchRngRB)
{
    Int  iMvShift = 2;
#if VCEG_AZ08_INTER_KLT_MV_BUGFIXED
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    iMvShift += VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
#endif
    TComMv cTmpMvPred = cMvPred;
    pcCU->clipMv(cTmpMvPred);

    rcMvSrchRngLT.setHor(cTmpMvPred.getHor() - (iSrchRng << iMvShift));
    rcMvSrchRngLT.setVer(cTmpMvPred.getVer() - (iSrchRng << iMvShift));

    rcMvSrchRngRB.setHor(cTmpMvPred.getHor() + (iSrchRng << iMvShift));
    rcMvSrchRngRB.setVer(cTmpMvPred.getVer() + (iSrchRng << iMvShift));
    pcCU->clipMv(rcMvSrchRngLT);
    pcCU->clipMv(rcMvSrchRngRB);

    rcMvSrchRngLT >>= iMvShift;
    rcMvSrchRngRB >>= iMvShift;
}

Void TComTrQuant::RecordPosition(UInt uiTargetCandiNum)
{
    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    DistType *pDiff = m_tempLibFast.getDiff();
    Short *pId = m_tempLibFast.getId();
    Short *pIdInteger = m_tempLibFast.getIdInteger();
    Int *pXInteger = m_tempLibFast.getXInteger();
    Int *pYInteger = m_tempLibFast.getYInteger();

    Int maxDiff = m_tempLibFast.getDiffMax();
    UInt uiInvalidCount = 0;
    Int k;
    for (k = uiTargetCandiNum - 1; k >= 0; k--)
    {
        if (pDiff[k] >= maxDiff)
        {
            uiInvalidCount++;
        }
        else
        {
            break;
        }
    }
    Int uiVaildCandiNum = uiTargetCandiNum - uiInvalidCount;
    memcpy(pXInteger, pX, sizeof(Int)*uiVaildCandiNum);
    memcpy(pYInteger, pY, sizeof(Int)*uiVaildCandiNum);
    memcpy(pIdInteger, pId, sizeof(Short)*uiVaildCandiNum);
    m_uiPartLibSize = uiVaildCandiNum;
}


Bool TComTrQuant::candidateTrain(UInt uiBlkSize, UInt uiTempSize)
{
    Bool bSucceedFlag = true;
    bSucceedFlag = prepareKLTSamplesInter(uiBlkSize, uiTempSize);
    if (bSucceedFlag)
    {
        bSucceedFlag = deriveKLT(uiBlkSize, m_uiVaildCandiNum);
    }
    return bSucceedFlag;
}

Bool TComTrQuant::prepareKLTSamplesInter(UInt uiBlkSize, UInt uiTempSize)
{
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiHeight = uiBlkSize;
    UInt uiWidth = uiHeight;
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    Pel **tarPatch = m_pppTarPatch[uiTarDepth];
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    UInt uiAllowMinCandiNum = g_uiDepth2MinCandiNum[uiTarDepth];
    //count collected candidate number
    Int k;
    Int maxDiff = m_tempLibFast.getDiffMax();
    DistType *pDiff = m_tempLibFast.getDiff();
    UInt uiInvalidCount = 0;
    for (k = uiTargetCandiNum - 1; k >= 0; k--)
    {
        if (pDiff[k] >= maxDiff)
        {
            uiInvalidCount++;
        }
        else
        {
            break;
        }
    }
    m_uiVaildCandiNum = uiTargetCandiNum - uiInvalidCount;
    if (m_uiVaildCandiNum < uiAllowMinCandiNum)
    {
        return false;
    }
    Int picStride = getStride();
    Pel predBlk[MAX_1DTRANS_LEN];
    Int i = 0;

    Pel *tarPatchRow;
    for (UInt uiY = uiTempSize; uiY < uiPatchSize; uiY++)
    {
        tarPatchRow = tarPatch[uiY];
        for (UInt uiX = uiTempSize; uiX < uiPatchSize; uiX++)
        {
            predBlk[i++] = tarPatchRow[uiX];
        }
    }

    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    Short *pId = m_tempLibFast.getId();
    Short setId;
    Pel *ref;
    Int iOffsetY, iOffsetX;
    Pel *refTarget;

    TrainDataType *pData;
    Int iCandiNum = m_uiVaildCandiNum;
    for (k = 0; k < iCandiNum; k++)
    {
        setId = pId[k];
        pData = m_pData[k];
        iOffsetY = pY[k];
        iOffsetX = pX[k];
        ref = getRefPicUsed(setId);
        refTarget = ref + iOffsetY*picStride + iOffsetX;
        i = 0;
        for (UInt uiY = 0; uiY < uiHeight; uiY++)
        {
            for (UInt uiX = 0; uiX < uiWidth; uiX++)
            {
                *pData++ = refTarget[uiX] - predBlk[i++];
            }
            refTarget += picStride;
        }
    }

#if VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY
    Int dim = uiWidth * uiHeight;
    Int d;
    for (k = 0; k < iCandiNum; k++)
    {
        for (d = 0; d < dim; d++)
        {
            m_pDataT[d][k] = m_pData[k][d];
        }
    }
#endif
    return true;
}

#endif



#if VCEG_AZ08_INTRA_KLT
Void TempLibFast::initTemplateDiff(UInt uiPatchSize, UInt uiBlkSize, Int bitDepth, Int iCandiNumber)
{
#if VCEG_AZ08_USE_SAD_DISTANCE
    DistType maxValue = ((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*(uiPatchSize*uiPatchSize - uiBlkSize*uiBlkSize);
#endif
#if VCEG_AZ08_USE_SSD_DISTANCE
    DistType maxValue = ((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS))*(uiPatchSize*uiPatchSize - uiBlkSize*uiBlkSize);
#endif
    m_diffMax = maxValue;
    for (Int i = 0; i < iCandiNumber; i++)
    {
        m_pDiff[i] = maxValue;
    }
}

Void TComTrQuant::getTargetTemplate(TComDataCU *pcCU, UInt uiAbsPartIdx, UInt uiBlkSize, UInt uiTempSize)
{
    const ComponentID compID = COMPONENT_Y;
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    Pel **tarPatch = m_pppTarPatch[uiTarDepth];
    UInt  uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
    Pel   *pCurrStart = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
    UInt  uiPicStride = pcCU->getPic()->getStride(compID);
    UInt uiY, uiX;

    //fill template
    //up-left & up 
    Pel  *tarTemp;
    Pel  *pCurrTemp = pCurrStart - uiTempSize*uiPicStride - uiTempSize;
    for (uiY = 0; uiY < uiTempSize; uiY++)
    {
        tarTemp = tarPatch[uiY];
        for (uiX = 0; uiX < uiPatchSize; uiX++)
        {
            tarTemp[uiX] = pCurrTemp[uiX];
        }
        pCurrTemp += uiPicStride;
    }
    //left
    for (uiY = uiTempSize; uiY < uiPatchSize; uiY++)
    {
        tarTemp = tarPatch[uiY];
        for (uiX = 0; uiX < uiTempSize; uiX++)
        {
            tarTemp[uiX] = pCurrTemp[uiX];
        }
        pCurrTemp += uiPicStride;
    }
}

Void TComTrQuant::candidateSearchIntra(TComDataCU *pcCU, UInt uiPartAddr, UInt uiBlkSize, UInt uiTempSize)
{
    const ComponentID compID = COMPONENT_Y;
    const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    Pel **tarPatch = getTargetPatch(uiTarDepth);
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    TComPic *pPic = pcCU->getPic();
    TComPicYuv *pPicYuv = pPic->getPicYuvRec();
    //Initialize the library for saving the best candidates
    m_tempLibFast.initTemplateDiff(uiPatchSize, uiBlkSize, channelBitDepth, uiTargetCandiNum);
    Short setId = 0; //record the reference picture.
    searchCandidateFromOnePicIntra(pcCU, uiPartAddr, pPic, pPicYuv, tarPatch, uiPatchSize, uiTempSize, setId);
}

Void  TComTrQuant::searchCandidateFromOnePicIntra(TComDataCU *pcCU, UInt uiPartAddr, TComPic* refPicSrc, TComPicYuv *refPic, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, UInt setId)
{
    const ComponentID compID = COMPONENT_Y;
    UInt uiBlkSize = uiPatchSize - uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];
    if (TMPRED0_TMPREDKLT1_ORI2 == 0)
    {
        uiTargetCandiNum = min((UInt)TMPRED_CANDI_NUM, uiTargetCandiNum);
    }
    UInt  uiLibSizeMinusOne = uiTargetCandiNum - 1;

    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    DistType *pDiff = m_tempLibFast.getDiff();
    Short *pId = m_tempLibFast.getId();
    Int     refStride = refPic->getStride(compID);
    Int zOrder = pcCU->getZorderIdxInCtu() + uiPartAddr;
    Pel *ref = refPic->getAddr(compID, pcCU->getCtuRsAddr(), zOrder);
    setRefPicUsed(setId, ref); //facilitate the access of each candidate point 
    setStride(refPic->getStride(compID));

    Int     iSrchRng = SEARCHRANGEINTRA;
    TComMv  cMvSrchRngLT;
    TComMv  cMvSrchRngRB;

    Int  iMvShift = 0;
    TComMv cTmpMvPred;
    cTmpMvPred.setZero();

    UInt uiCUPelY = pcCU->getCUPelY();
    UInt uiCUPelX = pcCU->getCUPelX();
    Int blkX = g_auiRasterToPelX[g_auiZscanToRaster[uiPartAddr]];
    Int blkY = g_auiRasterToPelY[g_auiZscanToRaster[uiPartAddr]];
    Int iCurrY = uiCUPelY + blkY;
    Int iCurrX = uiCUPelX + blkX;
    Int offsetLCUY = g_auiRasterToPelY[g_auiZscanToRaster[zOrder]]; //offset in this LCU
    Int offsetLCUX = g_auiRasterToPelX[g_auiZscanToRaster[zOrder]];

    Int iYOffset, iXOffset;
    DistType diff;
    DistType *pDiffEnd = &pDiff[uiLibSizeMinusOne];
    Pel *refCurr;

#define REGION_NUM 3
    Int mvYMins[REGION_NUM];
    Int mvYMaxs[REGION_NUM];
    Int mvXMins[REGION_NUM];
    Int mvXMaxs[REGION_NUM];
    Int regionNum = REGION_NUM;
    Int regionId = 0;

    //1. check the near pixels within LCU
    //above pixels in LCU
    Int iTemplateSize = uiTempSize;
    Int iBlkSize = uiBlkSize;
    regionId = 0;
#if JVET_C0024_QTBT
    const UInt maxCUWidth = pcCU->getPic()->getPicSym()->getSPS().getCTUSize();
    const UInt maxCUHeight = pcCU->getPic()->getPicSym()->getSPS().getCTUSize();
#else
    const UInt maxCUWidth = pcCU->getPic()->getPicSym()->getSPS().getMaxCUWidth();
    const UInt maxCUHeight = pcCU->getPic()->getPicSym()->getSPS().getMaxCUHeight();
#endif

    Int iVerMin = max(((iTemplateSize) << iMvShift), (iCurrY - offsetLCUY - iBlkSize + 1) << iMvShift);
    Int iVerMax = (iCurrY - iBlkSize) << iMvShift;
    Int iHorMin = max((iTemplateSize) << iMvShift, (iCurrX - offsetLCUX - iBlkSize + 1) << iMvShift);
    Int iHorMax = min((iCurrX - offsetLCUX + maxCUWidth - iBlkSize) << iMvShift, refPicSrc->getSlice(0)->getSPS()->getPicWidthInLumaSamples() - iBlkSize);

    mvXMins[regionId] = iHorMin - iCurrX;
    mvXMaxs[regionId] = iHorMax - iCurrX;
    mvYMins[regionId] = iVerMin - iCurrY;
    mvYMaxs[regionId] = iVerMax - iCurrY;
    //left pixels in LCU
    regionId = 1;
    iVerMin = max(((iTemplateSize) << iMvShift), (iCurrY - iBlkSize + 1) << iMvShift);
    iVerMax = min((iCurrY - offsetLCUY + maxCUHeight - iBlkSize) << iMvShift, refPicSrc->getSlice(0)->getSPS()->getPicHeightInLumaSamples() - iBlkSize);

    iHorMin = max((iTemplateSize) << iMvShift, (iCurrX - offsetLCUX - iBlkSize + 1) << iMvShift);
    iHorMax = (iCurrX - iBlkSize) << iMvShift;
    mvXMins[regionId] = iHorMin - iCurrX;
    mvXMaxs[regionId] = iHorMax - iCurrX;
    mvYMins[regionId] = iVerMin - iCurrY;
    mvYMaxs[regionId] = iVerMax - iCurrY;
    Int combinedX = offsetLCUX + iBlkSize - 1;
    Int combinedY = offsetLCUY + iBlkSize - 1;
#if JVET_C0024_QTBT
    Int NumInRow = pcCU->getSlice()->getSPS()->getCTUSize() >> 2;
#else
    Int NumInRow = pcCU->getSlice()->getSPS()->getMaxCUHeight() >> 2;
#endif

    //check within LCU pixels
    for (regionId = 0; regionId < 2; regionId++)
    {
        Int mvYMin = mvYMins[regionId];
        Int mvYMax = mvYMaxs[regionId];
        Int mvXMin = mvXMins[regionId];
        Int mvXMax = mvXMaxs[regionId];
        if (mvYMax < mvYMin || mvXMax < mvXMin)
        {
            continue;
        }
        for (iYOffset = mvYMax; iYOffset >= mvYMin; iYOffset--)
        {
            for (iXOffset = mvXMax; iXOffset >= mvXMin; iXOffset--)
            {
                refCurr = ref + iYOffset*refStride + iXOffset;
                Int iLCUX = iXOffset + combinedX;
                Int iLCUY = iYOffset + combinedY;
                Int ZorderTmp = getZorder(iLCUX, iLCUY, NumInRow);
                if (ZorderTmp >= zOrder)
                {
                    //Ignore the blocks that have not been coded.
                    continue;
                }
                diff = calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                if (diff < (*pDiffEnd))
                {
                    insertNode(diff, iXOffset, iYOffset, pDiff, pX, pY, pId, uiLibSizeMinusOne, setId);
                }
            }
        }
    }

    //2. check the pixels outside LCU
    for (regionId = 0; regionId < regionNum; regionId++)
    {
        pcCU->clipMvIntraConstraint(regionId, mvXMins[regionId], mvXMaxs[regionId], mvYMins[regionId], mvYMaxs[regionId], iSrchRng, uiTempSize, uiBlkSize, iCurrY, iCurrX, offsetLCUY, offsetLCUX);
    }
    for (regionId = 0; regionId < regionNum; regionId++)
    {
        Int mvYMin = mvYMins[regionId];
        Int mvYMax = mvYMaxs[regionId];
        Int mvXMin = mvXMins[regionId];
        Int mvXMax = mvXMaxs[regionId];
        if (mvXMax < mvXMin)
        {
            continue;
        }
        for (iYOffset = mvYMax; iYOffset >= mvYMin; iYOffset--)
        {
            for (iXOffset = mvXMax; iXOffset >= mvXMin; iXOffset--)
            {
                refCurr = ref + iYOffset*refStride + iXOffset;
                diff = calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchSize, uiTempSize, *pDiffEnd);
                if (diff < (*pDiffEnd))
                {
                    insertNode(diff, iXOffset, iYOffset, pDiff, pX, pY, pId, uiLibSizeMinusOne, setId);
                }
            }
        }
    }
}

Bool TComTrQuant::generateTMPrediction(Pel *piPred, UInt uiStride, UInt uiBlkSize, UInt uiTempSize, Int &foundCandiNum)
{
    Bool bSucceedFlag = true;
    UInt uiPatchSize = uiBlkSize + uiTempSize;
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
     //count collected candidate number
    DistType *pDiff = m_tempLibFast.getDiff();
    DistType maxDiff = m_tempLibFast.getDiffMax();
    Int k;
    UInt uiInvalidCount = 0;
    UInt uiTargetCandiNum = g_uiDepth2MaxCandiNum[uiTarDepth];

    for (k = uiTargetCandiNum - 1; k >= 0; k--)
    {
        if (pDiff[k] >= maxDiff)
        {
            uiInvalidCount++;
        }
        else
        {
            break;
        }
    }
    m_uiVaildCandiNum = uiTargetCandiNum - uiInvalidCount;
    foundCandiNum = m_uiVaildCandiNum;
    Int iCandiNum;
    iCandiNum = min((UInt)TMPRED_CANDI_NUM, m_uiVaildCandiNum);
    if (iCandiNum < 1)
    {
        return false;
    }

    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    Short setId;
    Pel *ref;
    Int picStride = getStride();
    Int iOffsetY, iOffsetX;
    Pel *refTarget;
    UInt uiHeight = uiPatchSize - uiTempSize;
    UInt uiWidth = uiHeight;
    TrainDataType *pData;

    //the data center: we use the prediction block as the center now.
    Pel predBlk[MAX_1DTRANS_LEN] = { 0 };
    UInt i = 0;
    //collect the candidates
    setId = 0;
    ref = getRefPicUsed(setId);

    for (k = 0; k < iCandiNum; k++)
    {
        pData = m_pData[k];
        iOffsetY = pY[k];
        iOffsetX = pX[k];
        refTarget = ref + iOffsetY*picStride + iOffsetX;
        for (UInt uiY = 0; uiY < uiHeight; uiY++)
        {
            for (UInt uiX = 0; uiX < uiWidth; uiX++)
            {
                *pData++ = refTarget[uiX];
            }
            refTarget += picStride;
        }
    }
    //average of the first several candidates as prediction
    Int iSize = uiWidth*uiHeight;
    for (k = 0; k < iCandiNum; k++)
    {
        pData = m_pData[k];
        for (i = 0; i < iSize; i++)
        {
            predBlk[i] += pData[i];
        }
    }
    Int iShift = iCandiNum >> 1;
    for (i = 0; i < iSize; i++)
    {
        predBlk[i] = (predBlk[i] + iShift) / iCandiNum;
    }

    Pel*  pPred = piPred;
    i = 0;
    for (UInt uiY = 0; uiY < uiHeight; uiY++)
    {
        for (UInt uiX = 0; uiX < uiWidth; uiX++)
        {
            pPred[uiX] = predBlk[i++];
        }
        pPred += uiStride;
    }
    return bSucceedFlag;
}

Bool TComTrQuant::calcKLTIntra(Pel *piPred, UInt uiStride, UInt uiBlkSize)
{
    Bool bSucceedFlag = prepareKLTSamplesIntra(piPred, uiStride, uiBlkSize);
    if (bSucceedFlag == false)
    {
        return bSucceedFlag;
    }
    bSucceedFlag = deriveKLT(uiBlkSize, m_uiVaildCandiNum);
    return bSucceedFlag;
}

Bool TComTrQuant::prepareKLTSamplesIntra(Pel *piPred, UInt uiStride, UInt uiBlkSize)
{
    UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
    UInt uiAllowMinCandiNum = g_uiDepth2MinCandiNum[uiTarDepth];
    if (m_uiVaildCandiNum < uiAllowMinCandiNum)
    {
        return false;
    }
    UInt uiHeight = uiBlkSize;
    UInt uiWidth = uiHeight;
    Int k;
    Pel*  pPred = piPred;
    Int i = 0;
    Pel predBlk[MAX_1DTRANS_LEN] = { 0 };
    for (UInt uiY = 0; uiY < uiHeight; uiY++)
    {
        for (UInt uiX = 0; uiX < uiWidth; uiX++)
        {
            predBlk[i++] = pPred[uiX];
        }
        pPred += uiStride;
    }
    Int *pX = m_tempLibFast.getX();
    Int *pY = m_tempLibFast.getY();
    Short setId;
    Pel *ref;
    Int picStride = getStride();
    Int iOffsetY, iOffsetX;
    Pel *refTarget;

    TrainDataType *pData;
    Int iCandiNum = m_uiVaildCandiNum;
    setId = 0;
    ref = getRefPicUsed(setId);
    for (k = 0; k < iCandiNum; k++)
    {
        pData = m_pData[k];
        iOffsetY = pY[k];
        iOffsetX = pX[k];
        refTarget = ref + iOffsetY*picStride + iOffsetX;
        i = 0;
        for (UInt uiY = 0; uiY < uiHeight; uiY++)
        {
            for (UInt uiX = 0; uiX < uiWidth; uiX++)
            {
                *pData++ = refTarget[uiX] - predBlk[i++];

            }
            refTarget += picStride;
        }
    }
#if VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY
    Int dim = uiWidth * uiHeight;
    Int d;
    for (k = 0; k < m_uiVaildCandiNum; k++)
    {
        for (d = 0; d < dim; d++)
        {
            m_pDataT[d][k] = m_pData[k][d];
        }
    }
#endif
    return true;
}

#if VCEG_AZ08_USE_SSE_TMP_SAD
DistType TComTrQuant::calcTemplateDiff(Pel *ref, UInt uiStride, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, DistType iMax)
{
    Int iY, iX;
    Int iDiffSum = 0;
    Pel *refPatchRow = ref - uiTempSize*uiStride - uiTempSize;
    Pel *tarPatchRow;
    static int id = 0;
    id++;
    for (iY = 0; iY < uiTempSize; iY++)
    {
        tarPatchRow = tarPatch[iY];
#if VCEG_AZ08_USE_SSE_TMP_SAD
        iDiffSum += AbsSumOfVector(refPatchRow, tarPatchRow, uiPatchSize);
#else
        for (iX = 0; iX < uiPatchSize; iX++)
        {
            iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
        }
        iDiffSum += tmpDiff2;
#endif
        if (iDiffSum > iMax) //for speeding up
        {
            return iDiffSum;
        }
        refPatchRow += uiStride;
    }

    if (uiTempSize > 2 && uiTempSize <= 8)
    {
        for (iY = uiTempSize; iY < uiPatchSize; iY++)
        {
            tarPatchRow = tarPatch[iY];
#if VCEG_AZ08_USE_SSE_TMP_SAD
            iDiffSum += AbsSumOfVectorLesseqthan8(refPatchRow, tarPatchRow, uiTempSize);
#else
            for (iX = 0; iX < uiTempSize; iX++)
            {
                iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
            }
#endif
            if (iDiffSum > iMax)
            {
                return iDiffSum;
            }
            refPatchRow += uiStride;
        }
    }
    else
    {
        for (iY = uiTempSize; iY < uiPatchSize; iY++)
        {
            tarPatchRow = tarPatch[iY];
            for (iX = 0; iX < uiTempSize; iX++)
            {
                iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
            }
            if (iDiffSum > iMax)
            {
                return iDiffSum;
            }
            refPatchRow += uiStride;
        }
    }
    return iDiffSum;
}
#else
DistType TComTrQuant::calcTemplateDiff(Pel *ref, UInt uiStride, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, DistType iMax)
{
    Int iY, iX;
#if VCEG_AZ08_USE_SSD_DISTANCE
    Int iDiff;
#endif
    DistType iDiffSum = 0;
    Pel *refPatchRow = ref - uiTempSize*uiStride - uiTempSize;
    Pel *tarPatchRow;
    for (iY = 0; iY < uiTempSize; iY++)
    {
        tarPatchRow = tarPatch[iY];
        for (iX = 0; iX < uiPatchSize; iX++)
        {
#if VCEG_AZ08_USE_SAD_DISTANCE
            iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
#endif
#if VCEG_AZ08_USE_SSD_DISTANCE
            iDiff = refPatchRow[iX] - tarPatchRow[iX];
            iDiffSum += iDiff * iDiff;
#endif
        }
        if (iDiffSum > iMax) //for speeding up
        {
            return iDiffSum;
        }
        refPatchRow += uiStride;
    }
    for (iY = uiTempSize; iY < uiPatchSize; iY++)
    {
        tarPatchRow = tarPatch[iY];
        for (iX = 0; iX < uiTempSize; iX++)
        {
#if VCEG_AZ08_USE_SAD_DISTANCE
            iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
#endif
#if VCEG_AZ08_USE_SSD_DISTANCE
            iDiff = refPatchRow[iX] - tarPatchRow[iX];
            iDiffSum += iDiff * iDiff;
#endif
        }
        if (iDiffSum > iMax) //for speeding up
        {
            return iDiffSum;
        }
        refPatchRow += uiStride;
    }
    return iDiffSum;
}
#endif
#endif

#if VCEG_AZ08_INTER_KLT
DistType TComTrQuant::calcPatchDiff(Pel *ref, UInt uiStride, Pel **tarPatch, UInt uiPatchSize, UInt uiTempSize, DistType iMax)
{
#if VCEG_AZ08_USE_SSE_BLK_SAD
    UInt uiBlkSize = uiPatchSize - uiTempSize;
    UInt blkDiffSum = 0;
    switch (uiBlkSize)
    {
    case 4:
    {
              blkDiffSum = GetSAD4x4_SSE_U16(tarPatch, ref, uiStride, uiTempSize, uiTempSize, iMax);
              break;
    }
    case 8:
    {
              blkDiffSum = GetSAD8x8_SSE_U16(tarPatch, ref, uiStride, uiTempSize, uiTempSize, iMax);
              break;
    }
    case 16:
    {
               blkDiffSum = GetSAD16x16_SSE_U16(tarPatch, ref, uiStride, uiTempSize, uiTempSize, iMax);
               break;
    }
    case 32:
    {
               blkDiffSum = GetSAD32x32_SSE_U16(tarPatch, ref, uiStride, uiTempSize, uiTempSize, iMax);
               break;
    }
    default:
    {
               //common calculation rather than using SSE for other cases.
               Int iY, iX;
               DistType iDiffSum = 0;
               Pel *refPatchRow = ref - uiTempSize*uiStride - uiTempSize;
               Pel *tarPatchRow;
               for (iY = 0; iY < uiPatchSize; iY++)
               {
                   tarPatchRow = tarPatch[iY];
                   for (iX = 0; iX < uiPatchSize; iX++)
                   {
                       iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
                   }
                   if (iDiffSum > iMax) //for speeding up
                   {
                       break;
                   }
                   refPatchRow += uiStride;
               }
    }
    }

    if (uiTempSize >0 && blkDiffSum < iMax)
    {
        Int remainMax = iMax - blkDiffSum;
        UInt tempDiffSum = calcTemplateDiff(ref, uiStride, tarPatch, uiPatchSize, uiTempSize, remainMax);
        blkDiffSum += tempDiffSum;
    }
    return blkDiffSum;
#else
    Int iY, iX;
#if VCEG_AZ08_USE_SSD_DISTANCE
    Int iDiff;
#endif
    DistType iDiffSum = 0;
    Pel *refPatchRow = ref - uiTempSize*uiStride - uiTempSize;
    Pel *tarPatchRow;
    for (iY = 0; iY < uiPatchSize; iY++)
    {
        tarPatchRow = tarPatch[iY];
        for (iX = 0; iX < uiPatchSize; iX++)
        {
#if VCEG_AZ08_USE_SAD_DISTANCE
            iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
#endif
#if VCEG_AZ08_USE_SSD_DISTANCE
            iDiff = refPatchRow[iX] - tarPatchRow[iX];
            iDiffSum += iDiff * iDiff;
#endif
        }
        if (iDiffSum > iMax) //for speeding up
        {
            break;
        }
        refPatchRow += uiStride;
    }
    return iDiffSum;
#endif
}
#endif





//For SSE speed up
#if VCEG_AZ08_USE_SSE_SPEEDUP
#define USE_SUM 1 //If defined, faster
#if USE_SUM
inline Int GetSum(__m128i x)
{
    __m128i hi = _mm_srli_si128(x, 8);
    __m128i low64_16bit = _mm_add_epi16(hi, x); //0+4; 1+5; 2+6; 3+7
    __m128i low64_32bit = _mm_cvtepi16_epi32(low64_16bit); //16 bytes

    __m128i hi2 = _mm_srli_si128(low64_32bit, 8);
    low64_32bit = _mm_add_epi32(hi2, low64_32bit); //0+4 + 2+6;  1+5 + 3+7;

    return low64_32bit.m128i_i32[0] + low64_32bit.m128i_i32[1];
}
#endif

float InnerProduct_SSE_FLOATXSHORT(float *pa, short *pb, int m)
{
    const int mask = 0xf1; //11110001
    int m_int = (m >> 2) << 2;
    int m_remain = m % 4;
    int i;
    __m128 Xf, Yf, innerprod;
    __m128i Ys, Y32;
    float sumfinal = 0;
    for (i = 0; i < m_int; i += 4)
    {
        Xf = _mm_castsi128_ps(_mm_loadu_si128((__m128i*) (pa + i))); //4 
        Ys = _mm_loadu_si128((__m128i*) (pb + i));
        Y32 = _mm_cvtepi16_epi32(Ys);
        Yf = _mm_cvtepi32_ps(Y32); //Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        innerprod = _mm_dp_ps(Xf, Yf, mask);
        sumfinal += innerprod.m128_f32[0];
    }
    if (m_remain > 0)
    {
        for (i = m_int; i < m; i++)
        {
            sumfinal += pa[i] * pb[i];
        }
    }
    return sumfinal;
}

int InnerProduct_SSE_SHORT(short *pa, short *pb, int m)
{
    int m_int = (m >> 3) << 3;
    int m_remain = m % 8;
    int i;
    __m128i X, Y, X32, Y32, prod;
    int sumfinal = 0;
    for (i = 0; i < m_int; i += 8)
    {
        X = _mm_loadu_si128((__m128i*) (pa + i));
        Y = _mm_loadu_si128((__m128i*) (pb + i));
        prod = _mm_madd_epi16(X, Y);
        sumfinal += prod.m128i_i32[0] + prod.m128i_i32[1] + prod.m128i_i32[2] + prod.m128i_i32[3];
    }
    if (m_remain > 0)
    {
        int m_remain_remain = m_remain % 4;
        int m_int_int = (m_remain >> 2) << 2;
        if (m_int_int > 0)
        {
            i = m_int;
            X = _mm_loadl_epi64(reinterpret_cast<const __m128i*> (pa + i)); //Load the lower 64 bits of the value pointed to by p into the lower 64 bits of the result, zeroing the upper 64 bits of the result.
            Y = _mm_loadl_epi64(reinterpret_cast<const __m128i*> (pb + i));
            X32 = _mm_cvtepi16_epi32(X); //performs a conversion of signed integers from 16-bit to 32-bit.
            Y32 = _mm_cvtepi16_epi32(Y);
            prod = _mm_mullo_epi32(X32, Y32); //multiplies two sets of 32-bit signed integers
            sumfinal += prod.m128i_i32[0] + prod.m128i_i32[1] + prod.m128i_i32[2] + prod.m128i_i32[3];
        }
        if (m_remain_remain > 0)
        {
            for (i = m_int + m_int_int; i < m; i++)
            {
                sumfinal += pa[i] * pb[i];
            }
        }
    }
    return sumfinal;
}

inline int MY_INT(float x)
{
    return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5);
}

void scaleVector(float *px, __m128 scale, short *pout, int m)
{
    int m_int = (m >> 2) << 2;
    int m_remain = m % 4;
    int i;
    __m128 Xf, Xscaled;
    float scalevalue = scale.m128_f32[0];
    for (i = 0; i < m_int; i += 4)
    {
        Xf = _mm_castsi128_ps(_mm_loadu_si128((__m128i*) (px + i))); //4 
        Xscaled = _mm_mul_ps(Xf, scale);
        pout[i] = MY_INT(Xscaled.m128_f32[0]);
        pout[i + 1] = MY_INT(Xscaled.m128_f32[1]);
        pout[i + 2] = MY_INT(Xscaled.m128_f32[2]);
        pout[i + 3] = MY_INT(Xscaled.m128_f32[3]);
    }
    if (m_remain > 0)
    {
        for (i = m_int; i < m; i++)
        {
            pout[i] = MY_INT(px[i] * scalevalue);
        }
    }
}

void scaleMatrix(float **ppx, short **ppout, float scale, int rows, int cols)
{
    __m128 scaleArray;
    scaleArray.m128_f32[0] = scale;
    scaleArray.m128_f32[1] = scale;
    scaleArray.m128_f32[2] = scale;
    scaleArray.m128_f32[3] = scale;

    for (int iRow = 0; iRow < rows; iRow++)
    {
        scaleVector(ppx[iRow], scaleArray, ppout[iRow], cols);
    }
}

UInt GetSAD4x4_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD)
{
    int blkSize = 4;
    __m128i s0, s1, s2, s3;
    __m128i r0, r1, r2, r3;

    r0 = _mm_cvtsi64_si128(*(__int64*)(pRef)); //first line
    r1 = _mm_cvtsi64_si128(*(__int64*)(pRef + iRefStride)); //second line
    r2 = _mm_or_si128(_mm_slli_si128(r1, 8), r0); //first & second line

    s0 = _mm_cvtsi64_si128(*(__int64*)(pSrc[iYOffset] + iXOffset)); //first line
    s1 = _mm_cvtsi64_si128(*(__int64*)(pSrc[iYOffset + 1] + iXOffset)); //second line
    s2 = _mm_or_si128(_mm_slli_si128(s1, 8), s0); //first & second line
    r2 = _mm_subs_epi16(r2, s2);
    r2 = _mm_abs_epi16(r2);
#if USE_SUM
    UInt sum1 = GetSum(r2);
#else
    UInt sum1 = r2.m128i_i16[0] + r2.m128i_i16[1] + r2.m128i_i16[2] + r2.m128i_i16[3] + r2.m128i_i16[4] + r2.m128i_i16[5] + r2.m128i_i16[6] + r2.m128i_i16[7];
#endif
    if (sum1 > uiBestSAD)
    {
        return sum1;
    }

    r0 = _mm_cvtsi64_si128(*(__int64*)(pRef + 2 * iRefStride)); //3 line
    r1 = _mm_cvtsi64_si128(*(__int64*)(pRef + 3 * iRefStride)); //4 line
    r3 = _mm_or_si128(_mm_slli_si128(r1, 8), r0); //3 & 4 line

    s0 = _mm_cvtsi64_si128(*(__int64*)(pSrc[iYOffset + 2] + iXOffset)); //3 line
    s1 = _mm_cvtsi64_si128(*(__int64*)(pSrc[iYOffset + 3] + iXOffset)); //4 line
    s3 = _mm_or_si128(_mm_slli_si128(s1, 8), s0); //3 & 4 line

    r3 = _mm_subs_epi16(r3, s3);
    r3 = _mm_abs_epi16(r3);
#if USE_SUM
    UInt sum2 = GetSum(r3);
#else
    UInt sum2 = r3.m128i_i16[0] + r3.m128i_i16[1] + r3.m128i_i16[2] + r3.m128i_i16[3] + r3.m128i_i16[4] + r3.m128i_i16[5] + r3.m128i_i16[6] + r3.m128i_i16[7];
#endif
    return sum1 + sum2;
}

UInt GetSAD8x8_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD)
{
    int blkSize = 8;
    __m128i s0, s1;
    __m128i r0, r1;
    __m128i diff, diff0, diff1;
    UInt sumOfRow;
    UInt sum = 0;
    for (int row = 0; row < 8; row += 2)
    {
        r0 = _mm_loadu_si128((__m128i*) (pRef + iRefStride*row));
        s0 = _mm_loadu_si128((__m128i*) (pSrc[iYOffset + row] + iXOffset));
        diff0 = _mm_subs_epi16(r0, s0);
        diff0 = _mm_abs_epi16(diff0);

        r1 = _mm_loadu_si128((__m128i*) (pRef + iRefStride*(1 + row)));
        s1 = _mm_loadu_si128((__m128i*) (pSrc[iYOffset + row + 1] + iXOffset));
        diff1 = _mm_subs_epi16(r1, s1);
        diff1 = _mm_abs_epi16(diff1);

        diff = _mm_adds_epi16(diff0, diff1);
#if USE_SUM
        sumOfRow = GetSum(diff);
#else
        sumOfRow = diff.m128i_i16[0] + diff.m128i_i16[1] + diff.m128i_i16[2] + diff.m128i_i16[3] + diff.m128i_i16[4] + diff.m128i_i16[5] + diff.m128i_i16[6] + diff.m128i_i16[7];
#endif
        sum += sumOfRow;
        if (sum > uiBestSAD)
        {
            return sum;
        }
    }
    return sum;
}

UInt GetSAD16x16_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD)
{
    int blkSize = 16;
    __m128i s0, s1;
    __m128i r0, r1;
    __m128i diff, diff0, diff1;
    UInt sumOfRow;
    UInt sum = 0;
    for (int blkRow = 0; blkRow < 2; blkRow++)
    {
        for (int blkCol = 0; blkCol < 2; blkCol++)
        {
            Int iRelativeOffset = ((blkRow *iRefStride + blkCol) << 3);
            I16 *pRefCurr = pRef + iRelativeOffset;
            Int iYOffsetCurr = iYOffset + (blkRow << 3);
            Int iXOffsetCurr = iXOffset + (blkCol << 3);

            for (int row = 0; row < 8; row += 2)
            {
                r0 = _mm_loadu_si128((__m128i*) (pRefCurr + iRefStride*row));
                s0 = _mm_loadu_si128((__m128i*) (pSrc[iYOffsetCurr + row] + iXOffsetCurr));
                diff0 = _mm_subs_epi16(r0, s0);
                diff0 = _mm_abs_epi16(diff0);

                r1 = _mm_loadu_si128((__m128i*) (pRefCurr + iRefStride*(1 + row)));
                s1 = _mm_loadu_si128((__m128i*) (pSrc[iYOffsetCurr + row + 1] + iXOffsetCurr));
                diff1 = _mm_subs_epi16(r1, s1);
                diff1 = _mm_abs_epi16(diff1);

                diff = _mm_adds_epi16(diff0, diff1);
#if USE_SUM
                sumOfRow = GetSum(diff);
#else
                sumOfRow = diff.m128i_i16[0] + diff.m128i_i16[1] + diff.m128i_i16[2] + diff.m128i_i16[3] + diff.m128i_i16[4] + diff.m128i_i16[5] + diff.m128i_i16[6] + diff.m128i_i16[7];
#endif
                sum += sumOfRow;
                if (sum > uiBestSAD)
                {
                    return sum;
                }
            }
        }
    }
    return sum;
}

UInt GetSAD32x32_SSE_U16(I16 **pSrc, I16 *pRef, Int iRefStride, Int iYOffset, Int iXOffset, UInt uiBestSAD)
{
    int blkSize = 32;
    __m128i s0, s1;
    __m128i r0, r1;
    __m128i diff, diff0, diff1;
    UInt sumOfRow;
    UInt sum = 0;
    for (int blkRow = 0; blkRow < 4; blkRow++)
    {
        for (int blkCol = 0; blkCol < 4; blkCol++)
        {
            Int iRelativeOffset = ((blkRow *iRefStride + blkCol) << 3);
            I16 *pRefCurr = pRef + iRelativeOffset;
            Int iYOffsetCurr = iYOffset + (blkRow << 3);
            Int iXOffsetCurr = iXOffset + (blkCol << 3);
            for (int row = 0; row < 8; row += 2)
            {
                r0 = _mm_loadu_si128((__m128i*) (pRefCurr + iRefStride*row));
                s0 = _mm_loadu_si128((__m128i*) (pSrc[iYOffsetCurr + row] + iXOffsetCurr));
                diff0 = _mm_subs_epi16(r0, s0);
                diff0 = _mm_abs_epi16(diff0);

                r1 = _mm_loadu_si128((__m128i*) (pRefCurr + iRefStride*(1 + row)));
                s1 = _mm_loadu_si128((__m128i*) (pSrc[iYOffsetCurr + row + 1] + iXOffsetCurr));
                diff1 = _mm_subs_epi16(r1, s1);
                diff1 = _mm_abs_epi16(diff1);

                diff = _mm_adds_epi16(diff0, diff1);
#if USE_SUM
                sumOfRow = GetSum(diff);
#else
                sumOfRow = diff.m128i_i16[0] + diff.m128i_i16[1] + diff.m128i_i16[2] + diff.m128i_i16[3] + diff.m128i_i16[4] + diff.m128i_i16[5] + diff.m128i_i16[6] + diff.m128i_i16[7];
#endif
                sum += sumOfRow;
                if (sum > uiBestSAD)
                {
                    return sum;
                }
            }
        }
    }
    return sum;
}

// sum of abs for two vectors of m length
int AbsSumOfVector(short *pa, short *pb, int m)
{
    int sum = 0;
    int m_int = (m >> 3) << 3;
    int m_remain = m % 8;
    int i;
    __m128i a, b, diff;
    for (i = 0; i < m_int; i += 8)
    {
        a = _mm_loadu_si128((__m128i*) (pa + i));
        b = _mm_loadu_si128((__m128i*) (pb + i));
        diff = _mm_subs_epi16(a, b);
        diff = _mm_abs_epi16(diff);
        sum += GetSum(diff);
    }
    if (m_remain > 0)
    {
        i = m_int;
        a = _mm_loadu_si128((__m128i*) (pa + i));
        b = _mm_loadu_si128((__m128i*) (pb + i));
        diff = _mm_subs_epi16(a, b);
        diff = _mm_abs_epi16(diff);
        if (m_remain == 3)
        {
            sum += (diff.m128i_i16[0] + diff.m128i_i16[1] + diff.m128i_i16[2]);
        }
        else
        {
            for (i = 0; i < m_remain; i++)
            {
                sum += diff.m128i_i16[i];
            }
        }
    }
    return sum;
}

// sum of abs for two vectors of m length (m<=8)
int AbsSumOfVectorLesseqthan8(short *pa, short *pb, int m)
{
    //assert(m <= 8);
    int sum = 0;
    int i;
    __m128i a, b, diff;

    a = _mm_loadu_si128((__m128i*) (pa));
    b = _mm_loadu_si128((__m128i*) (pb));
    diff = _mm_subs_epi16(a, b);
    diff = _mm_abs_epi16(diff);
    for (i = 0; i < m; i++)
    {
        sum += diff.m128i_i16[i];
    }
    return sum;
}
#endif
//! \}
