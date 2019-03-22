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

/** \file     CommonDef.h
    \brief    Defines version information, constants and small in-line functions
*/

#ifndef __COMMONDEF__
#define __COMMONDEF__

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <limits>
#ifdef __GNUC__
#include <tmmintrin.h>
#include <smmintrin.h>
#endif

#if _MSC_VER > 1000
// disable "signed and unsigned mismatch"
#pragma warning( disable : 4018 )
// disable Bool coercion "performance warning"
#pragma warning( disable : 4800 )
#endif // _MSC_VER > 1000
#include "TypeDef.h"
#include "TComFxnVariantSelector.h"

#ifdef _MSC_VER
#if _MSC_VER <= 1500
inline Int64 abs (Int64 x) { return _abs64(x); };
#endif
#if RRC_DWDIP && __cplusplus < 201103L
#define round(a)  (((a) < 0)? (int)((a) - 0.5) : (int)((a) + 0.5))
#endif
#endif

//////////////////////////////////////////////////////////////////////////
// RRC
#include "../trace/trace_yuv.h"
#include "../trace/trace_flags.h"
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Version information
// ====================================================================================================================

#define NV_VERSION        "HM-16.6-JEM-6.0-HEM4.0"                 ///< Current software version

// ====================================================================================================================
// Platform information
// ====================================================================================================================

#ifdef __GNUC__
#define NVM_COMPILEDBY  "[GCC %d.%d.%d]", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__
#ifdef __IA64__
#define NVM_ONARCH    "[on 64-bit] "
#else
#define NVM_ONARCH    "[on 32-bit] "
#endif
#endif

#ifdef __INTEL_COMPILER
#define NVM_COMPILEDBY  "[ICC %d]", __INTEL_COMPILER
#elif  _MSC_VER
#define NVM_COMPILEDBY  "[VS %d]", _MSC_VER
#endif

#ifndef NVM_COMPILEDBY
#define NVM_COMPILEDBY "[Unk-CXX]"
#endif

#ifdef _WIN32
#define NVM_ONOS        "[Windows]"
#elif  __linux
#define NVM_ONOS        "[Linux]"
#elif  __CYGWIN__
#define NVM_ONOS        "[Cygwin]"
#elif __APPLE__
#define NVM_ONOS        "[Mac OS X]"
#else
#define NVM_ONOS "[Unk-OS]"
#endif

#define NVM_BITS          "[%d bit] ", (sizeof(Void*) == 8 ? 64 : 32) ///< used for checking 64-bit O/S

#ifndef NULL
#define NULL              0
#endif

// ====================================================================================================================
// Common constants
// ====================================================================================================================

static const UInt   MAX_UINT =                            0xFFFFFFFFU; ///< max. value of unsigned 32-bit integer
static const Int    MAX_INT =                              2147483647; ///< max. value of signed 32-bit integer
static const Double MAX_DOUBLE =                             1.7e+308; ///< max. value of Double-type value

// ====================================================================================================================
// Coding tool configuration
// ====================================================================================================================
// Most of these should not be changed - they resolve the meaning of otherwise magic numbers.

static const Int MAX_GOP =                                         64; ///< max. value of hierarchical GOP size
static const Int MAX_NUM_REF_PICS =                                16; ///< max. number of pictures used for reference
//intra KLT
#if VCEG_AZ08_INTRA_KLT
static const Int TMPRED0_TMPREDKLT1_ORI2 =                          1; ///< (default 1) 0: Template matching prediction; 1: TM prediction + KLT; 2: Original method
static const Int TMPRED_CANDI_NUM =                                 8; ///< Candidate number for intra prediction (should <= 32)
static const Int SEARCHRANGEINTRA =                                64; ///< Intra search range (-SEARCHRANGE,+SEARCHRANGE)
#endif
//inter KLT
#if VCEG_AZ08_INTER_KLT
static const Int SEARCHRANGE =                                     32; ///< (default 32) Search range for inter coding (-SEARCHRANGE,+SEARCHRANGE)
static const Int SEARCH_SIZE =                                     ((SEARCHRANGE << 1) + 1)*((SEARCHRANGE << 1) + 1);
#endif

#if VCEG_AZ08_KLT_COMMON
static const Int USE_MORE_BLOCKSIZE_DEPTH_MIN =                     1; ///< (default 1) To indicate minimum block size for KLT. 1~4 means 4x4, 8x8, 16x16, 32x32 respectively.
static const Int USE_MORE_BLOCKSIZE_DEPTH_MAX =                     4; ///< (default 4) To indicate maximum block size for KLT. 1~4 means 4x4, 8x8, 16x16, 32x32 respectively.
static const Int MAX_CANDI_NUM =                                  100; ///< Max allowable candidate number. The candidate number for different size blocks can be set respectively.
#if VCEG_AZ08_FAST_DERIVE_KLT
static const Int FAST_KLT_CANDINUM =                    MAX_CANDI_NUM; ///If MAX_CANDI_NUM > blkSize, fast algorithm will be performed.
#endif
static const Int KLTBASIS_SHIFTBIT =                               10; ///< KLT scale factor is BLOCK_SIZE*(1<<KLTBASIS_SHIFTBIT); (log2(width)+KLTBASIS_SHIFTBIT <= 15); If 6, then the first base vector is {64,...,64}. We use 10.
static const Int INIT_THRESHOULD_SHIFTBITS =                        2;  ///< (default 2) Early skip threshold for checking distance.
static const Int MAX_NUM_REF_IDS =                                 ((MAX_NUM_REF_PICS << 1) + (MAX_CANDI_NUM << 4)); 
static const Double IGNORE_THRESHOULD_OF_LARGEST =               1e-6;
#if VCEG_AZ08_FORCE_USE_GIVENNUM_BASIS
static const Int FORCE_BASIS_NUM =                                 32; /// Forced number of basis utilized (for speeding up).
#endif
#endif

static const Int MAX_NUM_REF =                                     16; ///< max. number of entries in picture reference list
static const Int MAX_QP =                                          51;
static const Int NOT_VALID =                                       -1;

static const Int AMVP_MAX_NUM_CANDS =                               2; ///< AMVP: advanced motion vector prediction - max number of final candidates
static const Int AMVP_MAX_NUM_CANDS_MEM =                           3; ///< AMVP: advanced motion vector prediction - max number of candidates
static const Int AMVP_DECIMATION_FACTOR =                           4;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
static const Int MRG_MAX_NUM_CANDS =                                7; ///< MERGE
#else
static const Int MRG_MAX_NUM_CANDS =                                5; ///< MERGE
#endif

static const Int MAX_TLAYER =                                       7; ///< Explicit temporal layer QP offset - max number of temporal layer

static const Int ADAPT_SR_SCALE =                                   1; ///< division factor for adaptive search range

static const Int MAX_NUM_PICS_IN_SOP =                           1024;

static const Int MAX_NESTING_NUM_OPS =                           1024;
static const Int MAX_NESTING_NUM_LAYER =                           64;

static const Int MAX_VPS_NUM_HRD_PARAMETERS =                       1;
static const Int MAX_VPS_OP_SETS_PLUS1 =                         1024;
static const Int MAX_VPS_NUH_RESERVED_ZERO_LAYER_ID_PLUS1 =         1;

static const Int MAXIMUM_INTRA_FILTERED_WIDTH =                    16;
static const Int MAXIMUM_INTRA_FILTERED_HEIGHT =                   16;

static const Int MAX_CPB_CNT =                                     32; ///< Upper bound of (cpb_cnt_minus1 + 1)
static const Int MAX_NUM_LAYER_IDS =                               64;

static const Int COEF_REMAIN_BIN_REDUCTION =                        3; ///< indicates the level at which the VLC transitions from Golomb-Rice to TU+EG(k)

static const Int CU_DQP_TU_CMAX =                                   5; ///< max number bins for truncated unary
static const Int CU_DQP_EG_k =                                      0; ///< expgolomb order

static const Int SBH_THRESHOLD =                                    4; ///< value of the fixed SBH controlling threshold

static const Int C1FLAG_NUMBER =                                    8; // maximum number of largerThan1 flag coded in one chunk:  16 in HM5
static const Int C2FLAG_NUMBER =                                    1; // maximum number of largerThan2 flag coded in one chunk:  16 in HM5

static const Int MAX_NUM_VPS =                                     16;
static const Int MAX_NUM_SPS =                                     16;
static const Int MAX_NUM_PPS =                                     64;

#if QTBT_NSST
static const Int NSST_SIG_NZ_LUMA =                                 1;
static const Int NSST_SIG_NZ_CHROMA =                               1;
#endif

#if JVET_E0077_LM_MF
static const Int LM_FILTER_NUM =                                    4;
#else
static const Int LM_FILTER_NUM =                                    0;
#endif
#if JVET_E0077_ENHANCED_LM
static const Int NUM_INTRA_MODE_NON_ANG = (1 + JVET_E0077_MMLM + LM_FILTER_NUM);
#endif

#if JVET_E0077_MMLM
static const Int MMLM_SAMPLE_NEIGHBOR_LINES =                       2;
#endif

#if JVET_E0077_ENHANCED_LM
static const Int  LM_SYMBOL_NUM = (2 + JVET_E0077_MMLM + LM_FILTER_NUM);
#endif

static const Int ZERO_OUT_TH =                                    32;  ///< threshold of zeor-out transform

#if COM16_C806_T64
#if JVET_C0024_QTBT
static const Int MLS_GRP_NUM =                                   1024; ///< Max number of coefficient groups, max(16, 256)
#else
static const Int MLS_GRP_NUM =                                    256; ///< Max number of coefficient groups, max(16, 256)
#endif
#else
static const Int MLS_GRP_NUM =                                     64; ///< Max number of coefficient groups, max(16, 64)
#endif
static const Int MLS_CG_LOG2_WIDTH =                                2;
static const Int MLS_CG_LOG2_HEIGHT =                               2;
static const Int MLS_CG_SIZE =                                      4; ///< Coefficient group size of 4x4; = MLS_CG_LOG2_WIDTH + MLS_CG_LOG2_HEIGHT

#if JVET_C0046_ZO_ASSERT
static const Int TH_LOG2TBAREASIZE =                               10; ///< Threshold for zeroing out
#endif

#if ADAPTIVE_QP_SELECTION
static const Int ARL_C_PRECISION =                                  7; ///< G382: 7-bit arithmetic precision
static const Int LEVEL_RANGE =                                     30; ///< G382: max coefficient level in statistics collection
#endif

static const Int RVM_VCEGAM10_M =                                   4;

#if VCEG_AZ07_INTRA_65ANG_MODES
static const Int FAST_UDI_MAX_RDMODE_NUM =                         67;
#else
static const Int FAST_UDI_MAX_RDMODE_NUM =                         35; ///< maximum number of RD comparison in fast-UDI estimation loop
#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_E0077_ENHANCED_LM
static const Int NUM_INTRA_MODE = (67 + NUM_INTRA_MODE_NON_ANG);
static const Int NUM_DIR = ((64 >> 2) + 1);
#else
static const Int NUM_INTRA_MODE =                                  68;
static const Int NUM_DIR =                (((NUM_INTRA_MODE-4)>>2)+1);
#endif
static const Int PLANAR_IDX =                                       0;
static const Int DC_IDX =                                           1; ///< index for intra DC mode
static const Int VER_IDX =                          (3*(NUM_DIR-1)+2); ///< index for intra VERTICAL   mode
static const Int HOR_IDX =                          (1*(NUM_DIR-1)+2); ///< index for intra HORIZONTAL mode
static const Int DIA_IDX =                          (2*(NUM_DIR-1)+2); ///< index for intra Diagonal mode
static const Int VDIA_IDX =                         (4*(NUM_DIR-1)+2); ///< index for intra DC mode


#if COM16_C806_LMCHROMA
#if JVET_E0077_ENHANCED_LM
static const Int NUM_CHROMA_MODE = (6 + JVET_E0077_MMLM + LM_FILTER_NUM); ///< total number of chroma modes
#else
static const Int NUM_CHROMA_MODE =                                  6; ///< total number of chroma modes
static const Int LM_CHROMA_IDX =                 (NUM_INTRA_MODE - 1); ///< chroma mode index for derived from LM mode
#endif
#else
static const Int NUM_CHROMA_MODE =                                  5; ///< total number of chroma modes
#endif
static const Int DM_CHROMA_IDX =                       NUM_INTRA_MODE; ///< chroma mode index for derived from luma intra mode
#else
static const Int NUM_INTRA_MODE =                                  36;
static const Int PLANAR_IDX =                                       0;
static const Int VER_IDX =                                         26; ///< index for intra VERTICAL   mode
static const Int HOR_IDX =                                         10; ///< index for intra HORIZONTAL mode
static const Int DC_IDX =                                           1; ///< index for intra DC mode
#if COM16_C1044_NSST
static const Int DIA_IDX =                                         18; ///< index for intra Diagonal mode
#endif
#if COM16_C806_LMCHROMA
static const Int NUM_CHROMA_MODE =                                  6; ///< total number of chroma modes
static const Int LM_CHROMA_IDX =                                   35; ///< chroma mode index for derived from LM mode
#else
static const Int NUM_CHROMA_MODE =                                  5; ///< total number of chroma modes
#endif
static const Int DM_CHROMA_IDX =                                   36; ///< chroma mode index for derived from luma intra mode
#endif

#if DIMD_NUM_INTRA_DIR_INC
static const Int EXT_HOR_IDX = 34;
static const Int EXT_DIA_IDX = 66;
static const Int EXT_VER_IDX = 98;
static const Int EXT_VDIA_IDX = 130;
#if VCEG_AZ07_INTRA_65ANG_MODES
#define MAP131TO67( mode )                 (mode<2?mode:((mode>>1)+1))
#define MAP67TO131( mode )                 (mode<2?mode:((mode<<1)-2))
#else
#define MAP131TO35( mode )       (mode<2?mode:(MAP131TO67(mode)>>1)+1)
#define MAP35TO131( mode )       MAP67TO131( (mode<2?mode:((mode<<1)-2)) )
#endif
#endif

#if RRC_DWDIP
#ifndef COM16_C1044_NSST 
#ifndef VCEG_AZ07_INTRA_65ANG_MODES
static const Int DIA_IDX = (VER_IDX + HOR_IDX) >> 1;
#endif
#endif 
static const Int RIGHT_ANGLE = VER_IDX - HOR_IDX; ///< right angle
static const Int FIRST_ANGULAR_IPM_IDX = DIA_IDX - RIGHT_ANGLE; ///< the 1st angular intra-prediction mode (-45 dgrs)
static const Int LAST_ANGULAR_IPM_IDX = DIA_IDX + RIGHT_ANGLE; ///< the last angular intra-prediction mode (135 dgrs)
#endif

#if JVET_E0062_MULTI_DMS
static const Int NUM_DM_MODES =                                     5; ///< total number of chroma DM modes
static const Int INVALID_CHROMAMODE = NUM_INTRA_MODE + NUM_CHROMA_MODE; ///< initilized value for chroma mode candidate list
#endif

#if COM16_C806_EMT
static const UChar INTER_MODE_IDX =                               255; ///< index for inter modes
#if JVET_C0024_QTBT
static const UInt  EMT_INTRA_MAX_CU =                              64; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64, 128
static const UInt  EMT_INTER_MAX_CU =                              64; ///< Max Inter CU size applying EMT, supported values: 8, 16, 32, 64, 128
#else
static const UInt  EMT_INTRA_MAX_CU =                              32; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32
static const UInt  EMT_INTER_MAX_CU =                              32; ///< Max Inter CU size applying EMT, supported values: 8, 16, 32
#endif
#endif

#if JVET_D0120_NSST_IMPROV
static const Int  NSST_HYGT_RNDS_4x4 =                              2;
static const Int  NSST_HYGT_RNDS_8x8 =                              4;
static const Int  NSST_HYGT_PTS =                            (1 << 8);
#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
static const Int MDCS_ANGLE_LIMIT =                                 9; ///< 0 = Horizontal/vertical only, 1 = Horizontal/vertical +/- 1, 2 = Horizontal/vertical +/- 2 etc...
#else
static const Int MDCS_ANGLE_LIMIT =                                 4; ///< 0 = Horizontal/vertical only, 1 = Horizontal/vertical +/- 1, 2 = Horizontal/vertical +/- 2 etc...
#endif

#if DIMD_NUM_INTRA_DIR_INC
static const Int MDCS_EXT_ANGLE_LIMIT = 19;
#endif

static const Int MDCS_MAXIMUM_WIDTH =                               8; ///< (measured in pixels) TUs with width greater than this can only use diagonal scan
static const Int MDCS_MAXIMUM_HEIGHT =                              8; ///< (measured in pixels) TUs with height greater than this can only use diagonal scan


static const Int LOG2_MAX_NUM_COLUMNS_MINUS1 =                      7;
static const Int LOG2_MAX_NUM_ROWS_MINUS1 =                         7;

static const Int CABAC_INIT_PRESENT_FLAG =                          1;

static const Int LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS =   4;
static const Int CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS = 8;

static const Int MAX_NUM_LONG_TERM_REF_PICS =                      33;
static const Int NUM_LONG_TERM_REF_PIC_SPS =                        0;


static const Int MAX_QP_OFFSET_LIST_SIZE =                          6; ///< Maximum size of QP offset list is 6 entries

// Cost mode support
static const Int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP =      0; ///< QP to use for lossless coding.
static const Int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME =4; ///< QP' to use for mixed_lossy_lossless coding.
#if COM16_C806_LMCHROMA
static const Int CR_FROM_CB_REG_COST_SHIFT =                        9;     
#endif

static const Int RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS =     4;
static const Int RExt__GOLOMB_RICE_INCREMENT_DIVISOR =              4;

static const Int RExt__PREDICTION_WEIGHTING_ANALYSIS_DC_PRECISION = 0; ///< Additional fixed bit precision used during encoder-side weighting prediction analysis. Currently only used when high_precision_prediction_weighting_flag is set, for backwards compatibility reasons.

static const Int MAX_TIMECODE_SEI_SETS =                            3; ///< Maximum number of time sets

#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT && !JVET_C0024_CTU_256
static const Int MAX_CU_DEPTH =                                     7; ///< log2(CTUSize)
static const Int MAX_CU_SIZE =                                    128; ///< = 1<<(MAX_CU_DEPTH)
#else
static const Int MAX_CU_DEPTH =                                     8; ///< log2(CTUSize)
static const Int MAX_CU_SIZE =                                    256; ///< = 1<<(MAX_CU_DEPTH)
#endif
#else
static const Int MAX_CU_DEPTH =                                     6; ///< log2(CTUSize)
static const Int MAX_CU_SIZE =                                     64; ///< = 1<<(MAX_CU_DEPTH)
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
static const Int MAX_PART_IN_CTU =                                  1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1); 
#endif
#if JVET_C0024_QTBT
static const Int MIN_PU_SIZE =                                      1<<MIN_CU_LOG2;
static const Int MIN_TU_SIZE =                                      1<<MIN_CU_LOG2;
#else
static const Int MIN_PU_SIZE =                                      4;
static const Int MIN_TU_SIZE =                                      4;
#endif
#if COM16_C806_T64
#if JVET_C0024_QTBT
static const Int MAX_TU_SIZE =                                      128; 
static const Int MAX_LOG2_TU_SIZE_PLUS_ONE =                        8; ///< log2(MAX_TU_SIZE) + 1 
#else
static const Int MAX_TU_SIZE =                                     64;
static const Int MAX_LOG2_TU_SIZE_PLUS_ONE =                        7; ///< log2(MAX_TU_SIZE) + 1
#endif
#else
static const Int MAX_TU_SIZE =                                     32;
#endif
static const Int MAX_NUM_PART_IDXS_IN_CTU_WIDTH = MAX_CU_SIZE/MIN_PU_SIZE; ///< maximum number of partition indices across the width of a CTU (or height of a CTU)
static const Int SCALING_LIST_REM_NUM =                             6;

static const Int QUANT_SHIFT =                                     14; ///< Q(4) = 2^14
static const Int IQUANT_SHIFT =                                     6;
static const Int SCALE_BITS =                                      15; ///< For fractional bit estimates in RDOQ

static const Int SCALING_LIST_NUM = MAX_NUM_COMPONENT * NUMBER_OF_PREDICTION_MODES; ///< list number for quantization matrix

static const Int SCALING_LIST_START_VALUE =                        8 ; ///< start value for dpcm mode
static const Int MAX_MATRIX_COEF_NUM =                            64 ; ///< max coefficient number for quantization matrix
static const Int MAX_MATRIX_SIZE_NUM =                             8 ; ///< max size number for quantization matrix
static const Int SCALING_LIST_BITS =                               8 ; ///< bit depth of scaling list entries
static const Int LOG2_SCALING_LIST_NEUTRAL_VALUE =                 4 ; ///< log2 of the value that, when used in a scaling list, has no effect on quantisation
static const Int SCALING_LIST_DC =                                16 ; ///< default DC value

static const Int CONTEXT_STATE_BITS =                              6 ;
#if COM16_C806_T64
#if JVET_C0024_QTBT
static const Int LAST_SIGNIFICANT_GROUPS =                        14 ;
#else
static const Int LAST_SIGNIFICANT_GROUPS =                        12 ;
#endif
#else
static const Int LAST_SIGNIFICANT_GROUPS =                        10 ;
#endif
#if VCEG_AZ07_FRUC_MERGE
static const Int FRUC_MERGE_OFF =                                0x0 ;
static const Int FRUC_MERGE_BILATERALMV =                        0x01;
static const Int FRUC_MERGE_TEMPLATE =                           0x02;
#if THIRD_FRUC_MODE
static const Int FRUC_MERGE_TEMPLATE_UNI =                       0x03;
#endif
static const Int FRUC_MERGE_TEMPLATE_SIZE =                        4 ;
static const Int FRUC_MERGE_REFINE_MVWEIGHT =                      4 ;
static const Int FRUC_MERGE_REFINE_MINBLKSIZE =                    4 ;
#endif
#if JVET_E0060_FRUC_CAND
static const Int NB_FRUC_CAND_ADDED =                              2 ; ///< for entire (AMVP and merge) CU, number of added spatial candidates in top, left, top-left, top-right, below-left <0-5>
static const Int NB_FRUC_CAND_ADDED_SUB =                          4 ; ///< for sub-blocks of merge CU, number of added spatial candidates in top, left, top-left, top-right, below-left <0-5>
static const UInt NB_FRUC_CAND_ATMVP =                             4 ; ///< for sub-blocks of merge CU, number of ATMVP candidates
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
static const Int MAX_GR_ORDER_RESIDUAL =                          10 ;
#endif
#if VCEG_AZ07_BAC_ADAPT_WDOW
#if !VCEG_AZ05_MULTI_PARAM_CABAC
static const Int ALPHA0 =                                           6; ///< 2^ALPHA0 is "window size" for probability up-date
#endif
static const Int CABAC_NUM_BINS =                              100000; ///< max number of bins for window calculation
static const Int NUM_WDOW =                                         4; ///< could be 16, 32, 64, 128
#endif
#if VCEG_AZ05_MULTI_PARAM_CABAC
static const Int ALPHA0 =                                           4; ///< 2^ALPHA0 is "window size" for probability up-date
#endif

#if COM16_C1016_AFFINE || MRG_CONSIDER_COMPLEX
static const Int AFFINE_MAX_NUM_V0 =                                3; ///< max number of motion candidates in top-left corner
static const Int AFFINE_MAX_NUM_V1 =                                2; ///< max number of motion candidates in top-right corner
static const Int AFFINE_MAX_NUM_V2 =                                2; ///< max number of motion candidates in left-bottom corner
static const Int AFFINE_MAX_NUM_COMB =                             12; ///< max number of combined motion candidates
static const Int AFFINE_MIN_BLOCK_SIZE =                            4; ///< Minimum affine MC block size
#endif

#if MRG_CONSIDER_COMPLEX
static const Int MRG_NODE_NUM =                                     4;
static const Int MRG_COMPLEX_NUM_CANDS =                           11;
static const Int MRG_COMPLEX_NUM_CANDS_REAL =                       8;
#endif

#if JVET_C0024_QTBT

#if JVET_C0024_AMAX_BT
static const Double AMAXBT_TH32 =                                   15.0;
static const Double AMAXBT_TH64 =                                   30.0;
#endif
 
#if JVET_E0023_FAST_ENCODING_SETTING
static const Int    SKIP_DEPTH =                                    SKIP_DEPTH_VALUE;
static const Int    FAST_SKIP_DEPTH =                               FAST_SKIP_DEPTH_VALUE;
#else
static const Int    SKIP_DEPTH =                                    3;
#endif
static const Int    SKIPHORNOVERQT_DEPTH_TH =                       2;

#if JVET_C0024_FAST_MRG
static const Int    NUM_MRG_SATD_CAND =                             4;
static const Double MRG_FAST_RATIO    =                             1.25;
#endif

#if JVET_C0024_PBINTRA_FAST
static const Double PBINTRA_RATIO     =                             1.1;
#endif

//QTBT high level parameters
//for I slice luma CTB configuration para.
static const Int    MAX_BT_DEPTH  =                                 4;      ///<  <=7
static const Int    MAX_BT_SIZE   =                                 32;     ///<  [1<<MIN_QT_SIZE, 1<<CTU_LOG2]
static const Int    MIN_BT_SIZE   =                                 4;      ///<  can be set down to 1<<MIN_CU_LOG2

//for I slice chroma CTB configuration para. (in luma samples)
static const Int    MAX_BT_DEPTH_C =                                0;      ///< <=7   
static const Int    MAX_BT_SIZE_C  =                               64;      ///< [1<<MIN_QT_SIZE_C, 1<<CTU_LOG2], in luma samples
static const Int    MIN_BT_SIZE_C  =                                4;      ///< can be set down to 4, in luma samples

//for P/B slice CTU config. para.
static const Int    MAX_BT_DEPTH_INTER =                            4;      ///< <=7
static const Int    MAX_BT_SIZE_INTER  =                          128;      ///< for initialization, [1<<MIN_BT_SIZE_INTER, 1<<CTU_LOG2]
static const Int    MIN_BT_SIZE_INTER  =                            4;      ///<

#if BT_FORCE
static const Int    MAX_FQT_LV =                                    1;      ///max force qt level
#if BT_FORCE_EXT
static const Int    FORCE_BTD_EXT =                                 3;      ///boudnary max bt depth extension, for HD maximal value is 3 (to apply 3 MAX_NUM_ROUTE need to extend to 48)
#endif
#if BT_FORCE_INTRA_ON
static const Int    MAX_FQT_LV_INTRA =                              2;      ///max I slice force qt level
#if BT_FORCE_EXT
static const Int    FORCE_BTD_EXT_INTRA =                           3;      ///I slice boudnary max bt depth extension, for HD maximal value is 3 (to apply 3 MAX_NUM_ROUTE need to extend to 48)
#endif
#endif
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO
static const Double JVET_D0077_SPLIT_DECISION_COST_SCALE =       1.05;
#endif
#endif

#if JVET_D0123_ME_CTX_LUT_BITS
static const Int    EPBIT =                                   (1<<15);
#endif

#if PVC_AQS
static const Int    ESM_SHIFT = 8;
#endif

#if JVET_E0052_DMVR
static const Int DMVR_INTME_RANGE =                                 1;
#endif

#if JVET_E0076_MULTI_PEL_MVD                  
static const Int MULTI_PEL_MVD_BITS =                               2;
#endif

#if WCG_LUMA_DQP_CM_SCALE
static const Int LUMA_LEVEL_TO_DQP_LUT_MAXSIZE =                 1024;      ///< max LUT size for QP offset based on luma
#endif

// ====================================================================================================================
// Macro functions
// ====================================================================================================================

#if COM16_C1016_AFFINE
#define SIGN(x)                       ( x >= 0 ? 1 : -1 )
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
// for clipping pixel values, use ClipA NOT Clip3
template <typename T> inline T Clip3(const T minVal, const T maxVal, const T a) { return (a<minVal)?minVal:(a>maxVal?maxVal:a); }  ///< general min/max clip

struct Bound { Int m,M; 
#if RRC_DWDIP_MLR_IMPR
               Double dMaxAngle[BLK_LENGTH_VALUES];
               Void reset() { std::fill_n(dMaxAngle, BLK_LENGTH_VALUES, std::numeric_limits<Double>::max()); }
               Void setOff() { std::fill_n(dMaxAngle, BLK_LENGTH_VALUES, 0.0); }
               Bound() { reset(); };
#endif
};

#if !!ITH_C1_CLIPA
struct ClipParam
{
	static Int nbBitsY;
	static Int nbBitsUV;
	static Int ibdLuma;
	static Int ibdChroma;
	static Int cquantiz;
	Bool isActive;
	Bool isChromaActive;

	Bound componentBounds[MAX_NUM_COMPONENT];

	Int max(ComponentID c) const
	{
		return componentBounds[c].M;
	}

	Int min(ComponentID c) const
	{
		return componentBounds[c].m;
	}

	Bound getComponentBound(ComponentID c) const
	{
		return componentBounds[c];
	}

	Bound Y()  const { return componentBounds[0]; }
	Bound &Y() { return componentBounds[0]; }
	Bound U()  const { return componentBounds[1]; }
	Bound &U() { return componentBounds[1]; }
	Bound V()  const { return componentBounds[2]; }
	Bound &V() { return componentBounds[2]; }

	void setInvalid()
	{
		componentBounds[0].m = componentBounds[0].M = 0;
		componentBounds[1] = componentBounds[2] = componentBounds[0];
	}

	Int operator[](int idx) const
	{
		switch (idx)
		{
		case 0: return Y().m; break;
		case 1: return Y().M; break;
		case 2: return U().m; break;
		case 3: return U().M; break;
		case 4: return V().m; break;
		default: return V().M; break;
		}
	}

	Int &operator[](int idx)
	{
		switch (idx)
		{
		case 0: return Y().m; break;
		case 1: return Y().M; break;
		case 2: return U().m; break;
		case 3: return U().M; break;
		case 4: return V().m; break;
		default: return V().M; break;
		}
	}
};

inline void resetBounds(ClipParam &prm)
{
	prm.Y().m = 0;
	prm.Y().M = (1 << ClipParam::ibdLuma) - 1; // (1<<BD)-1
	prm.U().m = 0;
	prm.U().M = (1 << ClipParam::ibdChroma) - 1;; // (1<<BD)-1
	prm.V() = prm.U();
}

inline void setOff(ClipParam &prm)
{
	resetBounds(prm);
	prm.isActive = false;
	prm.isChromaActive = false;
}

extern ClipParam g_ClipParam;

template <typename T> T ClipA(const T x, const ComponentID compID)
{
	return Clip3((T)g_ClipParam.min(compID), (T)g_ClipParam.max(compID), x);
}
#else
struct ClipParam{
    static Int nbBitsY;
    static Int nbBitsUV;
    static Int ibdLuma;
    static Int ibdChroma;
    static Int cquantiz;
    Bool isActive;
    Bool isChromaActive;
    Int max(ComponentID c) const {
        switch(c) {
        case COMPONENT_Y: return Y().M; break;
        case COMPONENT_Cb: return U().M; break;
        case COMPONENT_Cr: return V().M; break;
        default: assert(false); return 0;
        }
    }
    Int min(ComponentID c) const {
        switch(c) {
        case COMPONENT_Y: return Y().m; break;
        case COMPONENT_Cb: return U().m; break;
        case COMPONENT_Cr: return V().m; break;
        default: assert(false); return 0;
        }
    }
    Bound Y()  const { return Y_;  }
    Bound &Y()       { return Y_; }

    Bound Y_;

    Bound U()  const { return U_; }
    Bound &U()       { return U_; }
    Bound V()  const { return V_; }
    Bound &V()       { return V_; }
    Bound U_;
    Bound V_;


    void setInvalid() { Y_.m=Y_.M=0; U_=V_=Y_; }
    Int operator[](int idx) const {
        switch(idx) {
        case 0: return Y_.m; break;
        case 1: return Y_.M; break;
        case 2: return U_.m; break;
        case 3: return U_.M; break;
        case 4: return V_.m; break;
        default: return V_.M; break;
        }
    }
    Int &operator[](int idx) {
        switch(idx) {
        case 0: return Y_.m; break;
        case 1: return Y_.M; break;
        case 2: return U_.m; break;
        case 3: return U_.M; break;
        case 4: return V_.m; break;
        default: return V_.M; break;
        }
    }
};



inline void resetBounds(ClipParam &prm) {
    prm.Y().m = 0;
    prm.Y().M = (1<<ClipParam::ibdLuma)-1; // (1<<BD)-1
    prm.U().m = 0;
    prm.U().M = (1<<ClipParam::ibdChroma)-1;; // (1<<BD)-1
    prm.V()=prm.U();
#if RRC_DWDIP_MLR_IMPR
    prm.Y_.reset();
    prm.U_.reset();
    prm.V_.reset();
#endif
}

inline void setOff(ClipParam &prm) {
    resetBounds(prm);
    prm.isActive=false;
    prm.isChromaActive=false;
#if RRC_DWDIP_MLR_IMPR
    prm.Y_.setOff();
    prm.U_.setOff();
    prm.V_.setOff();
#endif
}

extern ClipParam g_ClipParam;

template <typename T> T ClipA(const T x, const ComponentID compID)
{
    switch(compID) {
    case COMPONENT_Y:    return Clip3((T)g_ClipParam.Y().m,(T)g_ClipParam.Y().M,x); break;
    case COMPONENT_Cb:   return Clip3((T)g_ClipParam.U().m,(T)g_ClipParam.U().M,x); break;
    case COMPONENT_Cr:   return Clip3((T)g_ClipParam.V().m,(T)g_ClipParam.V().M,x); break;
    default: std::cerr << "ClipA: Invalid compID value " << compID << " . Exiting." << std::endl; assert(false); exit(0);  return 0;
    }
}
#endif
#if EXTENSION_360_VIDEO
template <typename T> inline T ClipBD(const T x, const Int bitDepth)             { return Clip3(T(0), T((1 << bitDepth)-1), x);           }
#endif
#else
template <typename T> inline T Clip3 (const T minVal, const T maxVal, const T a) { return std::min<T> (std::max<T> (minVal, a) , maxVal); }  ///< general min/max clip
template <typename T> inline T ClipBD(const T x, const Int bitDepth)             { return Clip3(T(0), T((1 << bitDepth)-1), x);           }
#endif

template <typename T> inline Void Check3( T minVal, T maxVal, T a)
{
  if ((a > maxVal) || (a < minVal))
  {
    std::cerr << "ERROR: Range check " << minVal << " >= " << a << " <= " << maxVal << " failed" << std::endl;
    assert(false);
    exit(1);
  }
}  ///< general min/max clip

#define DATA_ALIGN                  1                                                                 ///< use 32-bit aligned malloc/free
#if     DATA_ALIGN && _WIN32 && ( _MSC_VER > 1300 )
#define xMalloc( type, len )        _aligned_malloc( sizeof(type)*(len), 32 )
#define xFree( ptr )                _aligned_free  ( ptr )
#else
#define xMalloc( type, len )        malloc   ( sizeof(type)*(len) )
#define xFree( ptr )                free     ( ptr )
#endif

#define FATAL_ERROR_0(MESSAGE, EXITCODE)                      \
{                                                             \
  printf(MESSAGE);                                            \
  exit(EXITCODE);                                             \
}

template <typename ValueType> inline ValueType leftShift       (const ValueType value, const Int shift) { return (shift >= 0) ? ( value                                  << shift) : ( value                                   >> -shift); }
template <typename ValueType> inline ValueType rightShift      (const ValueType value, const Int shift) { return (shift >= 0) ? ( value                                  >> shift) : ( value                                   << -shift); }
template <typename ValueType> inline ValueType leftShift_round (const ValueType value, const Int shift) { return (shift >= 0) ? ( value                                  << shift) : ((value + (ValueType(1) << (-shift - 1))) >> -shift); }
template <typename ValueType> inline ValueType rightShift_round(const ValueType value, const Int shift) { return (shift >= 0) ? ((value + (ValueType(1) << (shift - 1))) >> shift) : ( value                                   << -shift); }
#if O0043_BEST_EFFORT_DECODING
// when shift = 0, returns value
// when shift = 1, (value + 0 + value[1]) >> 1
// when shift = 2, (value + 1 + value[2]) >> 2
// when shift = 3, (value + 3 + value[3]) >> 3
template <typename ValueType> inline ValueType rightShiftEvenRounding(const ValueType value, const UInt shift) { return (shift == 0) ? value : ((value + (1<<(shift-1))-1 + ((value>>shift)&1)) >> shift) ; }
#endif

extern void twoDimBufCopier
(
	UChar *pDst,
	const Int iDstStride,
	const UChar *pSrc,
	const Int iSrcStride,
	const Int iBlkWd,
	const Int iBlkHt
	);

extern Void residueComputer
(
	Pel *pSrc,
	const Int iSrcStride,
	Pel *pPred,
	const Int iPredStride,
	Pel *pResidue,
	const Int iResidueStride,
	const Int iWidth,
	const Int iHeight
	);

extern Void reconComputer
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
	);

extern void scaleAddClip
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
	);

extern void scaleSub
(
	Pel *piSrc,
	Int iSrcStride,
	Pel *piDst,
	Int iDstStride,
	UInt cxWidth,
	UInt cxHeight,
	const Int iShift,
	Int iOffset
	);

extern void leftShiftSub
(
	const Pel *src,
	const Int srcStride,
	Pel *dst,
	const Int dstStride,
	const Int width,
	const Int height,
	const Int shift,
	const Pel offset
	);

extern void rightShiftClip
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
	);

extern void countZeroRowsZeroCols
(
	UInt uiWidth,
	UInt uiHeight,
	TCoeff *pRef,
	UInt uiSrcStride,
	UInt *puiZerosWidth,
	UInt *puiZerosHeight
	);

extern Void LongLong2DoubleConverter
(
	Int64 *pi64Src,
	Double *pdDst,
	Int iSize
	);
//! \}

#endif // end of #ifndef  __COMMONDEF__

