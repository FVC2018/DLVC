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

/** \file     TComRom.h
    \brief    global variables & functions (header)
*/

#ifndef __TCOMROM__
#define __TCOMROM__

#include "CommonDef.h"

#include<stdio.h>
#include<iostream>
#if GuoLei
extern Bool g_new;
#endif
#if VCEG_AZ08_KLT_COMMON
extern short **g_ppsEigenVector[USE_MORE_BLOCKSIZE_DEPTH_MAX];
#endif
//! \ingroup TLibCommon
//! \{

#if VCEG_AZ08_KLT_COMMON
#if VCEG_AZ08_INTER_KLT
extern Bool g_bEnableCheck;
#endif
Void         reOrderCoeff(TCoeff *pcCoef, const UInt *scan, UInt uiWidth, UInt uiHeight);
Void         recoverOrderCoeff(TCoeff *pcCoef, const UInt *scan, UInt uiWidth, UInt uiHeight);
#endif

#if VCEG_AZ08_INTRA_KLT
Int          getZorder(Int iLCUX, Int iLCUY, Int NumInRow);
#endif
// ====================================================================================================================
// Initialize / destroy functions
// ====================================================================================================================

Void         initROM();
Void         destroyROM();

// ====================================================================================================================
// Data structure related table & variable
// ====================================================================================================================
#if JVET_C0024_AMAX_BT
extern UInt g_uiBlkSize[ 10 ];
extern UInt g_uiNumBlk[ 10 ];
#if JVET_C0024_AMAX_BT_FIX
extern UInt g_uiPrevISlicePOC;
extern Bool g_bInitAMaxBT;
#endif
#endif

// flexible conversion from relative to absolute index
extern       UInt   g_auiZscanToRaster[ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ];
extern       UInt   g_auiRasterToZscan[ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ];
#if COM16_C806_T64
extern       UInt*  g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES][ MAX_LOG2_TU_SIZE_PLUS_ONE ][ MAX_LOG2_TU_SIZE_PLUS_ONE ];
#else
extern       UInt*  g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES][ MAX_CU_DEPTH ][ MAX_CU_DEPTH ];
#endif

Void         initZscanToRaster ( Int iMaxDepth, Int iDepth, UInt uiStartVal, UInt*& rpuiCurrIdx );
Void         initRasterToZscan ( UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxDepth         );

// conversion of partition index to picture pel position
extern       UInt   g_auiRasterToPelX[ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ];
extern       UInt   g_auiRasterToPelY[ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ];

Void         initRasterToPelXY ( UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxDepth );

#if !JVET_C0024_QTBT
extern const UInt g_auiPUOffset[NUMBER_OF_PART_SIZES];
#endif

extern const Int g_quantScales[SCALING_LIST_REM_NUM];             // Q(QP%6)
extern const Int g_invQuantScales[SCALING_LIST_REM_NUM];          // IQ(QP%6)

#if RExt__HIGH_PRECISION_FORWARD_TRANSFORM
static const Int g_transformMatrixShift[TRANSFORM_NUMBER_OF_DIRECTIONS] = { 14, 6 };
#else
static const Int g_transformMatrixShift[TRANSFORM_NUMBER_OF_DIRECTIONS] = {  6, 6 };
#endif

extern const TMatrixCoeff g_aiT4 [TRANSFORM_NUMBER_OF_DIRECTIONS][4][4];
extern const TMatrixCoeff g_aiT8 [TRANSFORM_NUMBER_OF_DIRECTIONS][8][8];
extern const TMatrixCoeff g_aiT16[TRANSFORM_NUMBER_OF_DIRECTIONS][16][16];
extern const TMatrixCoeff g_aiT32[TRANSFORM_NUMBER_OF_DIRECTIONS][32][32];

#if PVC_AQS
extern Int                g_aiFactorLuminance[256];
extern Int                g_aiFactorContrast[256];
extern Int                g_ESMFactor;
extern Int                g_ESMFactorIdx;
Void                      initFactorValue();
#endif

// ====================================================================================================================
// Luma QP to Chroma QP mapping
// ====================================================================================================================

static const Int chromaQPMappingTableSize = 58;

extern const UChar  g_aucChromaScale[NUM_CHROMA_FORMAT][chromaQPMappingTableSize];

#if BLOCK_BASED_BACKGROUND_REFERENCE
extern UInt g_uiMaxCUWidth;
extern UInt g_uiMaxCUHeight;
#endif


// ====================================================================================================================
// Scanning order & context mapping table
// ====================================================================================================================

extern const UInt   ctxIndMap4x4[4*4];

extern const UInt   g_uiGroupIdx[ MAX_TU_SIZE ];
extern const UInt   g_uiMinInGroup[ LAST_SIGNIFICANT_GROUPS ];

#if VCEG_AZ07_CTX_RESIDUALCODING
// ====================================================================================================================
// coefficients coding
// ====================================================================================================================
extern const UInt   g_auiGoRiceRange[MAX_GR_ORDER_RESIDUAL];                  //!< maximum value coded with Rice codes
#if !COM16_C806_T64
extern const UInt   g_uiLastCtx[ 28 ];
#endif
#endif

// ====================================================================================================================
// Intra prediction table
// ====================================================================================================================

#if JVET_C0024_QTBT
extern const UChar  g_aucIntraModeNumFast_UseMPM[7-MIN_CU_LOG2+1][7-MIN_CU_LOG2+1];
#else
extern const UChar  g_aucIntraModeNumFast_UseMPM[MAX_CU_DEPTH];
#endif
extern const UChar  g_aucIntraModeNumFast_NotUseMPM[MAX_CU_DEPTH];

extern const UChar  g_chroma422IntraAngleMappingTable[NUM_INTRA_MODE];

// ====================================================================================================================
// Mode-Dependent DST Matrices
// ====================================================================================================================

extern const TMatrixCoeff g_as_DST_MAT_4 [TRANSFORM_NUMBER_OF_DIRECTIONS][4][4];

#if COM16_C806_EMT
#if FAST_EMT_v1
extern Int g_aiTrSubsetIntra[3];
#else
extern Int g_aiTrSubsetIntra[3][2];
#endif
extern Int g_aiTrSubsetInter[4];
#if VCEG_AZ07_INTRA_65ANG_MODES
extern const UChar g_aucTrSetVert[NUM_INTRA_MODE-1];
extern const UChar g_aucTrSetHorz[NUM_INTRA_MODE-1];
#else
extern const UChar g_aucTrSetVert[35];
extern const UChar g_aucTrSetHorz[35];
#endif
extern const UInt g_iEmtSigNumThr;
#endif

#if COM16_C806_EMT || COM16_C806_T64
#if JVET_C0024_QTBT
extern TMatrixCoeff g_aiTr2 [NUM_TRANS_TYPE][ 2][ 2];
extern TMatrixCoeff g_aiTr128 [NUM_TRANS_TYPE][ 128][ 128];
#endif
extern TMatrixCoeff g_aiTr4 [NUM_TRANS_TYPE][ 4][ 4];
extern TMatrixCoeff g_aiTr8 [NUM_TRANS_TYPE][ 8][ 8];
extern TMatrixCoeff g_aiTr16[NUM_TRANS_TYPE][16][16];
extern TMatrixCoeff g_aiTr32[NUM_TRANS_TYPE][32][32];
extern TMatrixCoeff g_aiTr64[NUM_TRANS_TYPE][64][64];
#endif

#if COM16_C1044_NSST
extern const UChar g_NsstLut[NUM_INTRA_MODE-1];
#if JVET_D0120_NSST_IMPROV
struct tabSinCos { Int c, s; };
extern tabSinCos g_tabSinCos     [NSST_HYGT_PTS]; 
extern const Int g_nsstHyGTPermut4x4 [35][3][16];
extern const Int g_nsstHyGTPar4x4    [35][3][64];
extern const Int g_nsstHyGTPermut8x8 [35][3][64];
extern const Int g_nsstHyGTPar8x8    [35][3][768];
#else
extern const Int   g_aiNsst4x4[12][3][16][16];
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
extern const UInt  g_auiCoefScanFirstCG8x8[3][16];
#endif
#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
extern const UInt g_auiCoefTopLeftDiagScan8x8[5][64];
#else
extern const UInt g_auiCoefTopLeftDiagScan8x8[3][64];
#endif
#endif
#endif

#if VCEG_AZ07_INTRA_4TAP_FILTER
extern Int g_aiIntraCubicFilter[32][4];
extern Int g_aiIntraGaussFilter[32][4];
#if DIMD_NUM_INTRA_DIR_INC
extern Int g_aiExtIntraCubicFilter[64][4];
extern Int g_aiExtIntraGaussFilter[64][4];
#endif
#endif

#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
extern        UChar g_NonMPM[257];
#endif
// ====================================================================================================================
// Misc.
// ====================================================================================================================

extern       Char   g_aucConvertToBit  [ MAX_CU_SIZE+1 ];   // from width to log2(width)-2


#if ENC_DEC_TRACE
extern FILE*  g_hTrace;
extern Bool   g_bJustDoIt;
extern const Bool g_bEncDecTraceEnable;
extern const Bool g_bEncDecTraceDisable;
extern Bool   g_HLSTraceEnable;
extern UInt64 g_nSymbolCounter;

#define COUNTER_START    1
#define COUNTER_END      0 //( UInt64(1) << 63 )

#define DTRACE_CABAC_F(x)     if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "%f", x );
#define DTRACE_CABAC_V(x)     if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "%d", x );
#define DTRACE_CABAC_VL(x)    if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "%lld", x );
#define DTRACE_CABAC_T(x)     if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "%s", x );
#define DTRACE_CABAC_X(x)     if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "%x", x );
#define DTRACE_CABAC_R( x,y ) if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, x,    y );
#define DTRACE_CABAC_N        if ( ( g_nSymbolCounter >= COUNTER_START && g_nSymbolCounter <= COUNTER_END )|| g_bJustDoIt ) fprintf( g_hTrace, "\n"    );

#else

#define DTRACE_CABAC_F(x)
#define DTRACE_CABAC_V(x)
#define DTRACE_CABAC_VL(x)
#define DTRACE_CABAC_T(x)
#define DTRACE_CABAC_X(x)
#define DTRACE_CABAC_R( x,y )
#define DTRACE_CABAC_N

#endif

#if RRC_WIF
extern Int g_iFilterCoeff[7];
#endif

const Char* nalUnitTypeToString(NalUnitType type);

extern const Char *MatrixType[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
extern const Char *MatrixType_DC[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];

#if JVET_E0077_MMLM
extern Int g_aiLMDivTable[];
extern Int g_aiLMDivTableLow[];
extern Int g_aiLMDivTableHigh[];
#endif

extern const Int g_quantTSDefault4x4[4*4];
extern const Int g_quantIntraDefault8x8[8*8];
extern const Int g_quantInterDefault8x8[8*8];

extern const UInt g_scalingListSize [SCALING_LIST_SIZE_NUM];
extern const UInt g_scalingListSizeX[SCALING_LIST_SIZE_NUM];

#if COM16_C1046_PDPC_INTRA
#if JVET_G0104_PLANAR_PDPC
extern const Short g_pdpcParam[5][6];
#elif JVET_C0024_QTBT // lossless change, just remove unused entries from the table
extern const Int g_pdpc_pred_param[5][35][6];
#else
extern const Int g_pdpc_pred_param[5][2][35][7];
#endif
#endif

//! \}
#if JVET_E0077_ENHANCED_LM
Bool   IsLMMode(UInt uiIntraMode);
extern const Int g_aiLAP_MinSize[];
extern const Int g_aiMFLM_MinSize[];
extern const Int g_aiMMLM_MinSize[];
extern const Int g_aiNonLMPosThrs[];
#endif


#if RRC_DWDIP

#define DECL_DWDIP_TABLE_SYMBOL(NAME,W,H)  g_DWDIP##NAME##W##x##H##gen
#define DECL_DWDIP_TABLE(NAME,W,H)  DECL_DWDIP_TABLE_SYMBOL(NAME,W,H)[33][W*H]


#define DECL_DWDIP_TABLE_SET(NAME)           \
  Int DECL_DWDIP_TABLE(NAME, 4, 16);         \
  Int DECL_DWDIP_TABLE(NAME, 4, 32);         \
                                             \
  Int DECL_DWDIP_TABLE(NAME, 8, 16);         \
  Int DECL_DWDIP_TABLE(NAME, 8, 32);         \
                                             \
  Int DECL_DWDIP_TABLE(NAME, 16, 4);         \
  Int DECL_DWDIP_TABLE(NAME, 16, 8);         \
  Int DECL_DWDIP_TABLE(NAME, 16, 16);        \
  Int DECL_DWDIP_TABLE(NAME, 16, 32);        \
                                             \
  Int DECL_DWDIP_TABLE(NAME, 32, 4);         \
  Int DECL_DWDIP_TABLE(NAME, 32, 8);         \
  Int DECL_DWDIP_TABLE(NAME, 32, 16);        \
  Int DECL_DWDIP_TABLE(NAME, 32, 32);


#define DECL_SELECTOR_FUN(NAME) \
            Int * const getDWDIP##NAME##Gen( Int iWidth, Int iHeight, Int iDirMode) {                         \
              if (iDirMode < DIA_IDX) swap(iWidth, iHeight);                                                  \
              switch ( (iWidth<<8)|iHeight ) {                                                                \
              case (1040):  return DECL_DWDIP_TABLE_SYMBOL(NAME,4,16)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (1056):  return DECL_DWDIP_TABLE_SYMBOL(NAME,4,32)   [ abs(DIA_IDX - iDirMode) ]; break;   \
                                                                                                              \
              case (2064):  return DECL_DWDIP_TABLE_SYMBOL(NAME,8,16)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (2080):  return DECL_DWDIP_TABLE_SYMBOL(NAME,8,32)   [ abs(DIA_IDX - iDirMode) ]; break;   \
                                                                                                              \
              case (4100):  return DECL_DWDIP_TABLE_SYMBOL(NAME,16,4)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (4104):  return DECL_DWDIP_TABLE_SYMBOL(NAME,16,8)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (4112):  return DECL_DWDIP_TABLE_SYMBOL(NAME,16,16)  [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (4128):  return DECL_DWDIP_TABLE_SYMBOL(NAME,16,32)  [ abs(DIA_IDX - iDirMode) ]; break;   \
                                                                                                              \
              case (8196):  return DECL_DWDIP_TABLE_SYMBOL(NAME,32,4)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (8200):  return DECL_DWDIP_TABLE_SYMBOL(NAME,32,8)   [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (8208):  return DECL_DWDIP_TABLE_SYMBOL(NAME,32,16)  [ abs(DIA_IDX - iDirMode) ]; break;   \
              case (8224):  return DECL_DWDIP_TABLE_SYMBOL(NAME,32,32)  [ abs(DIA_IDX - iDirMode) ]; break;   \
              }                                                                                               \
              assert(0); return NULL;                                                                         \
            } 

extern Int * const getDWDIPweightsGen(Int iWidth, Int iHeight, Int iDirMode);
extern Int * const getDWDIPSideDecisionGen(Int iWidth, Int iHeight, Int iDirMode);

#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
extern Int * const getDWDIPSideDeltaPosGen(Int iWidth, Int iHeight, Int iDirMode);
#endif


#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
extern const Int64 g_iMLRpxlNumLUT   [6*6];
extern const Int64 g_iMLRaverXLUT    [6*6];
extern const Int64 g_iMLRaverYLUT    [6*6];

extern const Int64 g_iMLRfactorYp2LUT[6*6];
extern const Int64 g_iMLRfactorXp2LUT[6*6];
extern const UInt64 g_iMLRfactorXYLUT[6*6];
#endif

#endif // RRC_DWDIP


#endif  //__TCOMROM__
