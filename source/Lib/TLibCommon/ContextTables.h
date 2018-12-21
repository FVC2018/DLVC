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

/** \file     ContextTables.h
    \brief    Defines constants and tables for SBAC
    \todo     number of context models is not matched to actual use, should be fixed
*/

#ifndef __CONTEXTTABLES__
#define __CONTEXTTABLES__

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constants
// ====================================================================================================================
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
#define MAX_NUM_CTX_MOD             384       ///< number of context models for last coefficient position
#else
#define MAX_NUM_CTX_MOD             512       ///< maximum number of supported contexts
#endif

#if COM16_C806_LARGE_CTU
#define NUM_SPLIT_FLAG_CTX            5       ///< number of context models for split flag
#else
#define NUM_SPLIT_FLAG_CTX            3       ///< number of context models for split flag
#endif
#if JVET_C0024_QTBT
#if TT_SYNTAX
#define NUM_BTSPLIT_DIR_CTX           5
#define NUM_BTSPLIT_TYPE_CTX          3
#define BT_TYPE_CTX_OFFSET            (3+NUM_BTSPLIT_DIR_CTX)
#define NUM_BTSPLIT_MODE_CTX          (3+NUM_BTSPLIT_DIR_CTX+NUM_BTSPLIT_TYPE_CTX)      // [6,7,8] for tree type
#else
#define NUM_BTSPLIT_MODE_CTX          6
#endif
#endif
#define NUM_SKIP_FLAG_CTX             3       ///< number of context models for skip flag
#if BARC
#define NUM_BARC_FLAG_CTX             3       ///< number of context models for BARC flag
#endif
#if BARC_UP_SAMPLE_RDO
#define NUM_BARC_UP_FILT_RDO_FLAG_CTX 3
#endif
#if VCEG_AZ05_INTRA_MPI
#define NUM_MPI_CTX                   2       /// < number of context models for MPI Idx coding
#endif
#if COM16_C1046_PDPC_INTRA
#define NUM_PDPC_CTX                  2      /// < number of context models for MPI Idx coding
#endif

#if DIMD_INTRA_PRED
#define NUM_DIMD_CTX                  6
#define NUM_DIMD_NO_BT_CTX            6
#endif

#if RRC_DWDIP  
#define NUM_DWDIP_CTX                 1      /// < number of context models for DWDIP Idx coding
#endif
#if RRC_DWDIP_MLR
#define NUM_DWDIP_GP_CTX              1      /// < number of context models for DWDIP GP Idx coding
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
 #define NUM_ROT_TR_CTX               7       /// < number of context models for ROT Idx coding
#endif
#define NUM_MERGE_FLAG_EXT_CTX        1       ///< number of context models for merge flag of merge extended
#if COM16_C806_GEN_MRG_IMPROVEMENT
#define NUM_MERGE_IDX_EXT_CTX         5       ///< number of context models for merge index of merge extended
#else
#define NUM_MERGE_IDX_EXT_CTX         1       ///< number of context models for merge index of merge extended
#endif
#if COM16_C806_OBMC
#define NUM_OBMC_FLAG_CTX             1       ///< number of context models for OBMC flag
#endif

#if COM16_C1016_AFFINE
#define NUM_AFFINE_FLAG_CTX           3       ///< number of context models for affine flag
#endif

#if MRG_CONSIDER_COMPLEX
#define NUM_MERGE_COMPLEX_FLAG_EXT_CTX 3      ///< number of context models for merge flag of merge extended
#define NUM_MERGE_COMPLEX_IDX_EXT_CTX  1
#endif

#if VCEG_AZ07_IMV
#if  JVET_E0076_MULTI_PEL_MVD
#define NUM_IMV_FLAG_CTX              4       ///< number of context models for iMV flag
#else
#define NUM_IMV_FLAG_CTX              3       ///< number of context models for iMV flag
#endif
#endif
#if VCEG_AZ07_FRUC_MERGE
#define NUM_FRUCMGRMODE_CTX           3
#define NUM_FRUCME_CTX                3
#if THIRD_FRUC_MODE
#define NUM_FRUCME_UNI_CTX            2
#endif
#endif
#if VCEG_AZ06_IC
#define NUM_IC_FLAG_CTX               1       ///< number of context models for illumination compensation flag
#endif

#define NUM_PART_SIZE_CTX             4       ///< number of context models for partition size
#define NUM_PRED_MODE_CTX             1       ///< number of context models for prediction mode

#if MLIP_EXTEND_REF_LINE
#define NUM_EXTERN_REF_CTX            3
#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_B0051_NON_MPM_MODE
#define NUM_INTRA_PREDICT_CTX         12       ///< number of context models for intra prediction
#else
#define NUM_INTRA_PREDICT_CTX         9       ///< number of context models for intra prediction
#endif
#else
#define NUM_INTRA_PREDICT_CTX         1       ///< number of context models for intra prediction
#endif


#if JVET_E0077_ENHANCED_LM
#if JVET_E0062_MULTI_DMS
#define NUM_CHROMA_PRED_CTX           12   ///< number of context models for intra prediction (chroma)
#else
#define NUM_CHROMA_PRED_CTX           8     
#endif
#else
#if JVET_E0062_MULTI_DMS
#define NUM_CHROMA_PRED_CTX           6  ///< number of context models for intra prediction (chroma)
#else
#define NUM_CHROMA_PRED_CTX           2       ///< number of context models for intra prediction (chroma)
#endif
#endif

#define NUM_INTER_DIR_CTX             5       ///< number of context models for inter prediction direction
#define NUM_MV_RES_CTX                2       ///< number of context models for motion vector difference
#define NUM_CHROMA_QP_ADJ_FLAG_CTX    1       ///< number of context models for chroma_qp_adjustment_flag
#define NUM_CHROMA_QP_ADJ_IDC_CTX     1       ///< number of context models for chroma_qp_adjustment_idc

#define NUM_REF_NO_CTX                2       ///< number of context models for reference index
#if COM16_C806_T64
#define NUM_TRANS_SUBDIV_FLAG_CTX     4       ///< number of context models for transform subdivision flags
#else
#define NUM_TRANS_SUBDIV_FLAG_CTX     3       ///< number of context models for transform subdivision flags
#endif
#define NUM_QT_ROOT_CBF_CTX           1       ///< number of context models for QT ROOT CBF
#define NUM_DELTA_QP_CTX              3       ///< number of context models for dQP

#define NUM_SIG_CG_FLAG_CTX           2       ///< number of context models for MULTI_LEVEL_SIGNIFICANCE
#define NUM_EXPLICIT_RDPCM_FLAG_CTX   1       ///< number of context models for the flag which specifies whether to use RDPCM on inter coded residues
#define NUM_EXPLICIT_RDPCM_DIR_CTX    1       ///< number of context models for the flag which specifies which RDPCM direction is used on inter coded residues

#if VCEG_AZ08_KLT_COMMON
#define NUM_KLT_FLAG_CTX              1       ///< number of context models for KLT 
#endif
//--------------------------------------------------------------------------------------------------

// context size definitions for significance map
#if VCEG_AZ07_CTX_RESIDUALCODING
#define NUM_SIG_FLAG_CTX_LUMA         54      ///< number of context models for luma sig flag
#define NUM_SIG_FLAG_CTX_CHROMA       12      ///< number of context models for chroma sig flag
#define NUM_SIG_FLAG_CTX_LUMA_TU      18      ///< number of context models for luma sig flag per TU
#else
#define NUM_SIG_FLAG_CTX_LUMA        28       ///< number of context models for luma sig flag
#define NUM_SIG_FLAG_CTX_CHROMA      16       ///< number of context models for chroma sig flag

//                                                                                                           |----Luma-----|  |---Chroma----|
static const UInt significanceMapContextSetStart         [MAX_NUM_CHANNEL_TYPE][CONTEXT_NUMBER_OF_TYPES] = { {0,  9, 21, 27}, {0,  9, 12, 15} }; 
static const UInt significanceMapContextSetSize          [MAX_NUM_CHANNEL_TYPE][CONTEXT_NUMBER_OF_TYPES] = { {9, 12,  6,  1}, {9,  3,  3,  1} };
static const UInt nonDiagonalScan8x8ContextOffset        [MAX_NUM_CHANNEL_TYPE]                          = {  6,               0              };
static const UInt notFirstGroupNeighbourhoodContextOffset[MAX_NUM_CHANNEL_TYPE]                          = {  3,               0              };

//------------------

#define NEIGHBOURHOOD_00_CONTEXT_1_THRESHOLD_4x4  3
#define NEIGHBOURHOOD_00_CONTEXT_2_THRESHOLD_4x4  1

//------------------

#define FIRST_SIG_FLAG_CTX_LUMA                   0
#define FIRST_SIG_FLAG_CTX_CHROMA     (FIRST_SIG_FLAG_CTX_LUMA + NUM_SIG_FLAG_CTX_LUMA)
#endif
#define NUM_SIG_FLAG_CTX              (NUM_SIG_FLAG_CTX_LUMA + NUM_SIG_FLAG_CTX_CHROMA)       ///< number of context models for sig flag

//--------------------------------------------------------------------------------------------------

#if BLOCK_BASED_BACKGROUND_REFERENCE
#define NUM_CU_BG_FLAG_CTX            1
#endif

// context size definitions for last significant coefficient position

#define NUM_CTX_LAST_FLAG_SETS         2

#if COM16_C806_T64
#if JVET_C0024_QTBT
#define NUM_CTX_LAST_FLAG_XY          25      ///< number of context models for last coefficient position
#else
#define NUM_CTX_LAST_FLAG_XY          19      ///< number of context models for last coefficient position
#endif
#else
#define NUM_CTX_LAST_FLAG_XY          15      ///< number of context models for last coefficient position
#endif

//--------------------------------------------------------------------------------------------------

// context size definitions for greater-than-one and greater-than-two maps
#if !VCEG_AZ07_CTX_RESIDUALCODING
#define NUM_ONE_FLAG_CTX_PER_SET       4      ///< number of context models for greater than 1 flag in a set
#define NUM_ABS_FLAG_CTX_PER_SET       1      ///< number of context models for greater than 2 flag in a set

//------------------

#define NUM_CTX_SETS_LUMA              4      ///< number of context model sets for luminance
#define NUM_CTX_SETS_CHROMA            2      ///< number of context model sets for combined chrominance

#define FIRST_CTX_SET_LUMA             0      ///< index of first luminance context set
#endif
//------------------

#if VCEG_AZ07_CTX_RESIDUALCODING
#define NUM_ONE_FLAG_CTX_LUMA          16
#define NUM_ONE_FLAG_CTX_CHROMA        6                                                       ///< number of context models for greater than 1 flag of chroma
#else
#define NUM_ONE_FLAG_CTX_LUMA         (NUM_ONE_FLAG_CTX_PER_SET * NUM_CTX_SETS_LUMA)           ///< number of context models for greater than 1 flag of luma
#define NUM_ONE_FLAG_CTX_CHROMA       (NUM_ONE_FLAG_CTX_PER_SET * NUM_CTX_SETS_CHROMA)         ///< number of context models for greater than 1 flag of chroma

#define NUM_ABS_FLAG_CTX_LUMA         (NUM_ABS_FLAG_CTX_PER_SET * NUM_CTX_SETS_LUMA)           ///< number of context models for greater than 2 flag of luma
#define NUM_ABS_FLAG_CTX_CHROMA       (NUM_ABS_FLAG_CTX_PER_SET * NUM_CTX_SETS_CHROMA)         ///< number of context models for greater than 2 flag of chroma
#endif

#define NUM_ONE_FLAG_CTX              (NUM_ONE_FLAG_CTX_LUMA + NUM_ONE_FLAG_CTX_CHROMA)        ///< number of context models for greater than 1 flag
#if !VCEG_AZ07_CTX_RESIDUALCODING
#define NUM_ABS_FLAG_CTX              (NUM_ABS_FLAG_CTX_LUMA + NUM_ABS_FLAG_CTX_CHROMA)        ///< number of context models for greater than 2 flag


#define FIRST_CTX_SET_CHROMA          (FIRST_CTX_SET_LUMA + NUM_CTX_SETS_LUMA)                 ///< index of first chrominance context set
#endif
//--------------------------------------------------------------------------------------------------

// context size definitions for CBF

#define NUM_QT_CBF_CTX_SETS           2

#define NUM_QT_CBF_CTX_PER_SET        5       ///< number of context models for QT CBF

#define FIRST_CBF_CTX_LUMA            0       ///< index of first luminance CBF context

#define FIRST_CBF_CTX_CHROMA          (FIRST_CBF_CTX_LUMA + NUM_QT_CBF_CTX_PER_SET)  ///< index of first chrominance CBF context


//--------------------------------------------------------------------------------------------------

#define NUM_MVP_IDX_CTX               1       ///< number of context models for MVP index
#if CNN_BASED_INLOOP_FILTER
#define NUM_CNNLF_FLAG_CTX            1       ///< number of context models for convnet flag
#endif
#define NUM_SAO_MERGE_FLAG_CTX        1       ///< number of context models for SAO merge flags
#define NUM_SAO_TYPE_IDX_CTX          1       ///< number of context models for SAO type index

#define NUM_TRANSFORMSKIP_FLAG_CTX    1       ///< number of context models for transform skipping

#define NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX  1

#define NUM_CROSS_COMPONENT_PREDICTION_CTX 10

#if ALF_HM3_REFACTOR
#if JVET_C0038_GALF
#define NUM_ALF_CTRL_FLAG_CTX         1       ///< number of context models for ALF control flag
#else
#define NUM_ALF_CTRL_FLAG_CTX         3       ///< number of context models for ALF control flag
#endif
#define NUM_ALF_FLAG_CTX              1       ///< number of context models for ALF flag
#define NUM_ALF_UVLC_CTX              2       ///< number of context models for ALF UVLC (filter length)
#define NUM_ALF_SVLC_CTX              3       ///< number of context models for ALF SVLC (filter coeff.)
#endif

#if COM16_C806_EMT
#define NUM_EMT_TR_CTX               7       /// < number of context models for EMT Intra Idx coding
#define NUM_EMT_TU_IDX_CTX            4       ///< number of context models for EMT TU-level transform index
#if COM16_C806_LARGE_CTU
#define NUM_EMT_CU_FLAG_CTX           6       ///< number of context models for EMT CU-level flag
#else
#define NUM_EMT_CU_FLAG_CTX           4       ///< number of context models for EMT CU-level flag
#endif
#endif

#if RRC_NS
#define NUM_NS_FLAG_CTX             1       ///< number of context models for NS flags
#endif

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
#define NUM_QP_PROB                  52                //could be set to N (N>5, depending on the allowed QPs in a coded sequence)
#define NUM_CTX_PBSLICE              MAX_NUM_CTX_MOD //could be set to the exact number of used contexts later
#endif
#if RRC_SELECTABLE_IF
#define NUM_SELECTABLE_IF_CTX         1
#endif

#if RRC_INTRA_SHARPENING
#define NUM_SHARP_FLAG_CTX            1
#endif

#define CNU                          154      ///< dummy initialization value for unused context models 'Context model Not Used'


// ====================================================================================================================
// Tables
// ====================================================================================================================

#if BLOCK_BASED_BACKGROUND_REFERENCE
static const UChar
INIT_CU_BG_FLAG[3][NUM_CU_BG_FLAG_CTX] =
{
	{ 154 },
	{ 154 },
	{ 154 },
};
#endif

// initial probability for cu_transquant_bypass flag
static const UChar
INIT_CU_TRANSQUANT_BYPASS_FLAG[NUMBER_OF_SLICE_TYPES][NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX] =
{
  { 154 },
  { 154 },
  { 154 },
};

// initial probability for split flag
static const UChar
INIT_SPLIT_FLAG[NUMBER_OF_SLICE_TYPES][NUM_SPLIT_FLAG_CTX] =  
{
#if COM16_C806_LARGE_CTU
  { 107,  139,  126, 255, 0, },
  { 107,  139,  126, 255, 0, }, 
  { 139,  141,  157, 255, 0, }, 
#else
  { 107,  139,  126, },
  { 107,  139,  126, }, 
  { 139,  141,  157, }, 
#endif
};

#if JVET_C0024_QTBT
static const UChar
INIT_BTSPLIT_MODE[3][NUM_BTSPLIT_MODE_CTX] = 
{
#if TT_SYNTAX
  { 107,  139,  126, 154, 154, 154, 154, 154, 154, 154, 154, },
  { 107,  139,  126, 154, 154, 154, 154, 154, 154, 154, 154, },
  { 139,  141,  157, 154, 154, 154, 154, 154, 154, 154, 154, },
#else
  { 107,  139,  126, 154, 154, 154},
  { 107,  139,  126, 154, 154, 154}, 
  { 139,  141,  157, 154, 154, 154}, 
#endif
};
#endif

static const UChar
INIT_SKIP_FLAG[NUMBER_OF_SLICE_TYPES][NUM_SKIP_FLAG_CTX] =
{
  { 197,  185,  201, },
  { 197,  185,  201, },
  { CNU,  CNU,  CNU, },
};
#if BARC
static const UChar
INIT_BARC_FLAG[NUMBER_OF_SLICE_TYPES][NUM_BARC_FLAG_CTX] =
{
  { 197,  185,  201, },
  { 197,  185,  201, },
  { CNU,  CNU,  CNU, },
};
#endif
#if BARC_UP_SAMPLE_RDO
static const UChar
INIT_UP_CNN_RDO_FLAG[3][NUM_BARC_UP_FILT_RDO_FLAG_CTX] =
{
	{ 197, 185, 201, },
	{ 197, 185, 201, },
	{ CNU, CNU, CNU, },
};
#endif
#if VCEG_AZ05_INTRA_MPI
static const UChar
INIT_MPIIdx_FLAG[NUMBER_OF_SLICE_TYPES][NUM_MPI_CTX] =
{
  { 107, 107 },
  { 107, 107 },
  { 139, 139 },
};
#endif
#if COM16_C1046_PDPC_INTRA
static const UChar
INIT_PDPCIdx_FLAG[NUMBER_OF_SLICE_TYPES][NUM_PDPC_CTX] =
{
  { 107, 107 },
  { 107, 107 },
  { 139, 139 },
};
#endif

#if FAST_EMT_v1
static const UChar
INIT_EMT_TR_IDX[3][NUM_EMT_TR_CTX] =
{
	{ 107,107,107,107,107,107,107 },
	{ 107,107,107,107,107,107,107 },
	{ 139,139,139,139,139,139,139 },
};
#endif

#if DIMD_INTRA_PRED
static const UChar
INIT_DIMD_FLAG[NUMBER_OF_SLICE_TYPES][NUM_DIMD_CTX] =
{
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
};
static const UChar
INIT_DIMD_NO_BT_FLAG[NUMBER_OF_SLICE_TYPES][NUM_DIMD_NO_BT_CTX] =
{
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
};
#endif

#if RRC_DWDIP  
static const UChar
INIT_DWDIPIdx_FLAG[NUMBER_OF_SLICE_TYPES][NUM_DWDIP_CTX] =
{
  { 107 },
  { 107 },
  { 139 },
};
#endif
#if RRC_DWDIP_MLR
static const UChar
INIT_DWDIP_GPIdx_FLAG[NUMBER_OF_SLICE_TYPES][NUM_DWDIP_GP_CTX] =
{
  { 107 },
  { 107 },
  { 139 },
};
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
static const UChar 
INIT_ROT_TR_IDX[3][NUM_ROT_TR_CTX] =  
{
  { 107,107,107,107,107,107,107 }, 
  { 107,107,107,107,107,107,107 }, 
  { 139,139,139,139,139,139,139 }, 
};
#endif
static const UChar
INIT_MERGE_FLAG_EXT[NUMBER_OF_SLICE_TYPES][NUM_MERGE_FLAG_EXT_CTX] =
{
  { 154, },
  { 110, },
  { CNU, },
};

static const UChar
INIT_MERGE_IDX_EXT[NUMBER_OF_SLICE_TYPES][NUM_MERGE_IDX_EXT_CTX] =
{
#if COM16_C806_GEN_MRG_IMPROVEMENT
  { 137, CNU, CNU, CNU, CNU}, 
  { 122, CNU, CNU, CNU, CNU}, 
  { CNU, CNU, CNU, CNU, CNU}, 
#else
  { 137, },
  { 122, },
  { CNU, },
#endif
};

#if VCEG_AZ07_FRUC_MERGE
static const UChar
  INIT_FRUCMGRMODEBIN1[NUMBER_OF_SLICE_TYPES][NUM_FRUCMGRMODE_CTX] = 
{
  { 197,  185,  201 }, 
  { 197,  185,  201 }, 
  { CNU,  CNU,  CNU }, 
};

static const UChar
  INIT_FRUCMGRMODEBIN2[NUMBER_OF_SLICE_TYPES][NUM_FRUCME_CTX] = 
{
  { 197,  185,  201 }, 
  { 197,  185,  201 }, 
  { CNU,  CNU,  CNU }, 
};
#if THIRD_FRUC_MODE
static const UChar
INIT_FRUCMGRMODEBIN3[NUMBER_OF_SLICE_TYPES][NUM_FRUCME_UNI_CTX] =
{
  { 197,  140 },
  { 197,  140 },
  { CNU,  CNU },
};
#endif
#endif
#if VCEG_AZ07_IMV
static const UChar 
  INIT_IMV_FLAG[NUMBER_OF_SLICE_TYPES][NUM_IMV_FLAG_CTX] =
{
#if  JVET_E0076_MULTI_PEL_MVD
  { 197,  185,  201, 185,}, 
  { 197,  185,  201, 185,}, 
  { CNU,  CNU,  CNU, 185,}, 
#else
  { 197,  185,  201, }, 
  { 197,  185,  201, }, 
  { CNU,  CNU,  CNU, }, 
#endif
};
#endif

#if COM16_C806_OBMC
static const UChar 
INIT_OBMC_FLAG[NUMBER_OF_SLICE_TYPES][NUM_OBMC_FLAG_CTX] =  
{
  { 201, }, 
  { 201, }, 
  { CNU, }, 
};
#endif

#if COM16_C1016_AFFINE
static const UChar 
INIT_AFFINE_FLAG[3][NUM_AFFINE_FLAG_CTX] =  
{
  { 197,  185,  201, }, 
  { 197,  185,  201, }, 
  { CNU,  CNU,  CNU, }, 
};
#endif

#if MRG_CONSIDER_COMPLEX
static const Short
INIT_MERGE_COMPLEX_FLAG_EXT[3][NUM_MERGE_COMPLEX_FLAG_EXT_CTX] =
{
  { 196,  131,  233, },
  { 195,  131,  232, },
  { CNU,  CNU,  CNU, },
};

static const UChar
INIT_MERGE_COMPLEX_IDX_EXT[3][NUM_MERGE_COMPLEX_IDX_EXT_CTX] =
{
  { CNU, },
  { CNU, },
  { CNU, },
};
#endif

#if VCEG_AZ06_IC
static const UChar 
INIT_IC_FLAG[NUMBER_OF_SLICE_TYPES][NUM_IC_FLAG_CTX] =  
{
  { 154 },
  { 154 },
  { CNU },
};
#endif

static const UChar
INIT_PART_SIZE[NUMBER_OF_SLICE_TYPES][NUM_PART_SIZE_CTX] =
{
  { 154,  139,  154, 154 },
  { 154,  139,  154, 154 },
  { 184,  CNU,  CNU, CNU },
};

static const UChar
INIT_PRED_MODE[NUMBER_OF_SLICE_TYPES][NUM_PRED_MODE_CTX] =
{
  { 134, },
  { 149, },
  { CNU, },
};

#if MLIP_EXTEND_REF_LINE
static const UChar
INIT_EXTERN_REF[NUMBER_OF_SLICE_TYPES][NUM_EXTERN_REF_CTX] =
{
    { 154, 154, 154 },
    { 154, 154, 154 },
    { CNU, CNU, CNU },
};
#endif

static const UChar
INIT_INTRA_PRED_MODE[NUMBER_OF_SLICE_TYPES][NUM_INTRA_PREDICT_CTX] =
{
#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_B0051_NON_MPM_MODE
  { 183, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,184,184,184 }, 
  { 154, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,184,184,184 }, 
  { 184, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,184,184,184 },
#else
  { 183, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU }, 
  { 154, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU }, 
  { 184, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU },
#endif
#else
  { 183, },
  { 154, },
  { 184, },
#endif
};

#if JVET_E0077_ENHANCED_LM
static const UChar
INIT_CHROMA_PRED_MODE[NUMBER_OF_SLICE_TYPES][NUM_CHROMA_PRED_CTX] =
{
#if JVET_E0062_MULTI_DMS
    //LM, DM0,   DM1,   DM2,    DM3,    DM4,    DM5, DM6 
    { 139, 152, 139, 154, 154, 154, 154, 154, 154, 154, 154 , 154},
    { 139, 152, 139, 154, 154, 154, 154, 154, 154, 154, 154 , 154},
    { 139,  63, 139, 154, 154, 154, 154, 154, 154, 154, 154 , 154},
#else
    { 152, 139, 154, 154, 154, 154, 154, 154},
    { 152, 139, 154, 154, 154, 154, 154, 154},
    {  63, 139, 154, 154, 154, 154, 154, 154},
#endif
};
#else
#if JVET_E0062_MULTI_DMS
static const UChar
INIT_CHROMA_PRED_MODE[NUMBER_OF_SLICE_TYPES][NUM_CHROMA_PRED_CTX] =
{
  //LM,      DMIdx0, DMIdx1, DMIdx2, DMIdx3, DMIdx4
  { 139,      152,   139,    154,    154,    154 },
  { 139,      152,   139,    154,    154,    154 },
  { 139,       63,   139,    154,    154,    154 },
};

#else
static const UChar
INIT_CHROMA_PRED_MODE[NUMBER_OF_SLICE_TYPES][NUM_CHROMA_PRED_CTX] =
{
  { 152,  139, },
  { 152,  139, },
  {  63,  139, },
};
#endif
#endif


static const UChar
INIT_INTER_DIR[NUMBER_OF_SLICE_TYPES][NUM_INTER_DIR_CTX] =
{
  {  95,   79,   63,   31,  31, },
  {  95,   79,   63,   31,  31, },
  { CNU,  CNU,  CNU,  CNU, CNU, },
};

static const UChar
INIT_MVD[NUMBER_OF_SLICE_TYPES][NUM_MV_RES_CTX] =
{
  { 169,  198, },
  { 140,  198, },
  { CNU,  CNU, },
};

static const UChar
INIT_REF_PIC[NUMBER_OF_SLICE_TYPES][NUM_REF_NO_CTX] =
{
  { 153,  153 },
  { 153,  153 },
  { CNU,  CNU },
};

static const UChar
INIT_DQP[NUMBER_OF_SLICE_TYPES][NUM_DELTA_QP_CTX] =
{
  { 154,  154,  154, },
  { 154,  154,  154, },
  { 154,  154,  154, },
};

static const UChar
INIT_CHROMA_QP_ADJ_FLAG[NUMBER_OF_SLICE_TYPES][NUM_CHROMA_QP_ADJ_FLAG_CTX] =
{
  { 154, },
  { 154, },
  { 154, },
};

static const UChar
INIT_CHROMA_QP_ADJ_IDC[NUMBER_OF_SLICE_TYPES][NUM_CHROMA_QP_ADJ_IDC_CTX] =
{
  { 154, },
  { 154, },
  { 154, },
};

//--------------------------------------------------------------------------------------------------

//Initialisation for CBF

//                                 |---------Luminance---------|
#define BSLICE_LUMA_CBF_CONTEXT     153,  111,  CNU,  CNU,  CNU
#define PSLICE_LUMA_CBF_CONTEXT     153,  111,  CNU,  CNU,  CNU
#define ISLICE_LUMA_CBF_CONTEXT     111,  141,  CNU,  CNU,  CNU
//                                 |--------Chrominance--------|
#define BSLICE_CHROMA_CBF_CONTEXT   149,   92,  167,  154,  154
#define PSLICE_CHROMA_CBF_CONTEXT   149,  107,  167,  154,  154
#define ISLICE_CHROMA_CBF_CONTEXT    94,  138,  182,  154,  154


static const UChar
INIT_QT_CBF[NUMBER_OF_SLICE_TYPES][NUM_QT_CBF_CTX_SETS * NUM_QT_CBF_CTX_PER_SET] =
{
  { BSLICE_LUMA_CBF_CONTEXT, BSLICE_CHROMA_CBF_CONTEXT },
  { PSLICE_LUMA_CBF_CONTEXT, PSLICE_CHROMA_CBF_CONTEXT },
  { ISLICE_LUMA_CBF_CONTEXT, ISLICE_CHROMA_CBF_CONTEXT },
};


//--------------------------------------------------------------------------------------------------

static const UChar
INIT_QT_ROOT_CBF[NUMBER_OF_SLICE_TYPES][NUM_QT_ROOT_CBF_CTX] =
{
  {  79, },
  {  79, },
  { CNU, },
};


//--------------------------------------------------------------------------------------------------

//Initialisation for last-significant-position
#if VCEG_AZ07_CTX_RESIDUALCODING && !COM16_C806_T64
//                                           |------------------------------Luminance----------------------------------|
#define BSLICE_LUMA_LAST_POSITION_CONTEXT     110, 110,  94, 110, 140, 140, 111, 126, 126, 125, 126, 127, 143, 126, 125
#define PSLICE_LUMA_LAST_POSITION_CONTEXT     111, 125, 124, 111, 111, 111, 111, 126, 126, 110, 111, 141, 127, 111, 125
#define ISLICE_LUMA_LAST_POSITION_CONTEXT     125,  95, 109, 110, 125, 110, 125, 125, 110, 110, 154, 140, 140, 111, 111
//                                           |------------------------------Chrominance--------------------------------|
#define BSLICE_CHROMA_LAST_POSITION_CONTEXT   108, 108,  62, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#define PSLICE_CHROMA_LAST_POSITION_CONTEXT   109,  94,  63, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#define ISLICE_CHROMA_LAST_POSITION_CONTEXT   123,  93,  77, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#else
//                                           |------------------------------Luminance----------------------------------|
#define BSLICE_LUMA_LAST_POSITION_CONTEXT     125, 110, 124, 110,  95,  94, 125, 111, 111,  79, 125, 126, 111, 111,  79
#define PSLICE_LUMA_LAST_POSITION_CONTEXT     125, 110,  94, 110,  95,  79, 125, 111, 110,  78, 110, 111, 111,  95,  94
#define ISLICE_LUMA_LAST_POSITION_CONTEXT     110, 110, 124, 125, 140, 153, 125, 127, 140, 109, 111, 143, 127, 111,  79
//                                           |------------------------------Chrominance--------------------------------|
#define BSLICE_CHROMA_LAST_POSITION_CONTEXT   108, 123,  93, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#define PSLICE_CHROMA_LAST_POSITION_CONTEXT   108, 123, 108, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#define ISLICE_CHROMA_LAST_POSITION_CONTEXT   108, 123,  63, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU
#endif

#if COM16_C806_T64 
static const UChar 
INIT_LAST[NUMBER_OF_SLICE_TYPES][NUM_CTX_LAST_FLAG_SETS * NUM_CTX_LAST_FLAG_XY] =  
{
#if FIX_TICKET32
  { BSLICE_LUMA_LAST_POSITION_CONTEXT,    126,  111,  111,   79,  CNU, CNU, CNU, CNU, CNU, CNU,
    BSLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,  CNU, CNU, CNU, CNU, CNU, CNU,
  },                                                              
  { PSLICE_LUMA_LAST_POSITION_CONTEXT,    111,  111,   95,   94,  CNU, CNU, CNU, CNU, CNU, CNU,
    PSLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,  CNU, CNU, CNU, CNU, CNU, CNU,
  }, 
  { ISLICE_LUMA_LAST_POSITION_CONTEXT,    143,  127,  111,   79,  CNU, CNU, CNU, CNU, CNU, CNU,
    ISLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,  CNU, CNU, CNU, CNU, CNU, CNU,
  }, 
#else
  { BSLICE_LUMA_LAST_POSITION_CONTEXT,    126,  111,  111,   79,
    BSLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,
  }, 
  { PSLICE_LUMA_LAST_POSITION_CONTEXT,    111,  111,   95,   94,
    PSLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,
  }, 
  { ISLICE_LUMA_LAST_POSITION_CONTEXT,    143,  127,  111,   79, 
    ISLICE_CHROMA_LAST_POSITION_CONTEXT,  CNU,  CNU,  CNU,  CNU,
  }, 
#endif
};
#else
static const UChar
INIT_LAST[NUMBER_OF_SLICE_TYPES][NUM_CTX_LAST_FLAG_SETS * NUM_CTX_LAST_FLAG_XY] =
{
  { BSLICE_LUMA_LAST_POSITION_CONTEXT, BSLICE_CHROMA_LAST_POSITION_CONTEXT },
  { PSLICE_LUMA_LAST_POSITION_CONTEXT, PSLICE_CHROMA_LAST_POSITION_CONTEXT },
  { ISLICE_LUMA_LAST_POSITION_CONTEXT, ISLICE_CHROMA_LAST_POSITION_CONTEXT },
};
#endif

//--------------------------------------------------------------------------------------------------

static const UChar
INIT_SIG_CG_FLAG[NUMBER_OF_SLICE_TYPES][2 * NUM_SIG_CG_FLAG_CTX] =
{
#if VCEG_AZ07_CTX_RESIDUALCODING
  { 122, 143,
    91, 141,
  },
  { 78, 111,
    60, 140,
  },
  { 135, 155,
    104, 139,
  },
#else
  { 121,  140,
    61,  154,
  },
  { 121,  140,
    61,  154,
  },
  {  91,  171,
    134,  141,
  },
#endif
};


//--------------------------------------------------------------------------------------------------

//Initialisation for significance map
#if VCEG_AZ07_CTX_RESIDUALCODING 
//                                          |-------------------------------------------- 4x4 ---------------------------------------|-------------------------------------------- 8x8 ---------------------------------------|-------------------------------------------- 16x16&above---------------------------------|
#define BSLICE_LUMA_SIGNIFICANCE_CONTEXT     107, 139, 154, 140, 140, 141, 108, 154, 125, 155, 126, 127, 139, 155, 155, 141, 156, 143, 107, 139, 154, 140, 140, 141, 108, 154, 125, 155, 126, 127, 139, 155, 155, 141, 156, 143, 107, 139, 154, 140, 140, 141, 108, 154, 125, 155, 126, 127, 139, 155, 155, 141, 156, 143
#define PSLICE_LUMA_SIGNIFICANCE_CONTEXT     121, 167, 153, 139, 154, 140, 137, 168, 139, 154, 169, 155, 167, 169, 169, 184, 199, 156, 121, 167, 153, 139, 154, 140, 137, 168, 139, 154, 169, 155, 167, 169, 169, 184, 199, 156, 121, 167, 153, 139, 154, 140, 137, 168, 139, 154, 169, 155, 167, 169, 169, 184, 199, 156
#define ISLICE_LUMA_SIGNIFICANCE_CONTEXT     152, 139, 154, 154, 169, 155, 182, 154, 169, 184, 155, 141, 168, 214, 199, 170, 170, 171, 152, 139, 154, 154, 169, 155, 182, 154, 169, 184, 155, 141, 168, 214, 199, 170, 170, 171, 152, 139, 154, 154, 169, 155, 182, 154, 169, 184, 155, 141, 168, 214, 199, 170, 170, 171

#define BSLICE_CHROMA_SIGNIFICANCE_CONTEXT   137, 154, 154, 155, 155, 156, 124, 185, 156, 171, 142, 158
#define PSLICE_CHROMA_SIGNIFICANCE_CONTEXT   136, 153, 139, 154, 125, 140, 122, 154, 184, 185, 171, 157
#define ISLICE_CHROMA_SIGNIFICANCE_CONTEXT   167, 154, 169, 140, 155, 141, 153, 171, 185, 156, 171, 172

#else
//                                          |-DC-|  |-----------------4x4------------------|  |------8x8 Diagonal Scan------|  |----8x8 Non-Diagonal Scan----|  |-NxN First group-|  |-NxN Other group-| |-Single context-|
//                                          |    |  |                                      |  |-First Group-| |-Other Group-|  |-First Group-| |-Other Group-|  |                 |  |                 | |                |
#define BSLICE_LUMA_SIGNIFICANCE_CONTEXT     170,    154, 139, 153, 139, 123, 123,  63, 124,   166, 183, 140,  136, 153, 154,   166, 183, 140,  136, 153, 154,   166,   183,   140,   136,   153,   154,        140
#define PSLICE_LUMA_SIGNIFICANCE_CONTEXT     155,    154, 139, 153, 139, 123, 123,  63, 153,   166, 183, 140,  136, 153, 154,   166, 183, 140,  136, 153, 154,   166,   183,   140,   136,   153,   154,        140
#define ISLICE_LUMA_SIGNIFICANCE_CONTEXT     111,    111, 125, 110, 110,  94, 124, 108, 124,   107, 125, 141,  179, 153, 125,   107, 125, 141,  179, 153, 125,   107,   125,   141,   179,   153,   125,        141

//                                          |-DC-|  |-----------------4x4------------------|  |-8x8 Any group-|  |-NxN Any group-| |-Single context-|
#define BSLICE_CHROMA_SIGNIFICANCE_CONTEXT   170,    153, 138, 138, 122, 121, 122, 121, 167,   151,  183,  140,   151,  183,  140,        140
#define PSLICE_CHROMA_SIGNIFICANCE_CONTEXT   170,    153, 123, 123, 107, 121, 107, 121, 167,   151,  183,  140,   151,  183,  140,        140
#define ISLICE_CHROMA_SIGNIFICANCE_CONTEXT   140,    139, 182, 182, 152, 136, 152, 136, 153,   136,  139,  111,   136,  139,  111,        111
#endif
//------------------------------------------------

static const UChar
INIT_SIG_FLAG[NUMBER_OF_SLICE_TYPES][NUM_SIG_FLAG_CTX] =
{
  { BSLICE_LUMA_SIGNIFICANCE_CONTEXT, BSLICE_CHROMA_SIGNIFICANCE_CONTEXT },
  { PSLICE_LUMA_SIGNIFICANCE_CONTEXT, PSLICE_CHROMA_SIGNIFICANCE_CONTEXT },
  { ISLICE_LUMA_SIGNIFICANCE_CONTEXT, ISLICE_CHROMA_SIGNIFICANCE_CONTEXT },
};


//--------------------------------------------------------------------------------------------------

//Initialisation for greater-than-one flags and greater-than-two flags
#if VCEG_AZ07_CTX_RESIDUALCODING          
#define BSLICE_LUMA_ONE_CONTEXT    121, 135, 123, 124, 139, 125,  92, 124, 154, 125, 155, 138, 169, 155, 170, 156
#define PSLICE_LUMA_ONE_CONTEXT    165,  75, 152, 153, 139, 154, 121, 138, 139, 154, 140, 167, 183, 169, 170, 156
#define ISLICE_LUMA_ONE_CONTEXT    196, 105, 152, 153, 139, 154, 136, 138, 139, 169, 140, 196, 183, 169, 170, 171

#define BSLICE_CHROMA_ONE_CONTEXT  166, 152, 140, 170, 171, 157
#define PSLICE_CHROMA_ONE_CONTEXT  193, 181, 169, 170, 171, 172
#define ISLICE_CHROMA_ONE_CONTEXT  195, 181, 169, 170, 156, 157
#else
//                                 |------Set 0-------| |------Set 1-------| |------Set 2-------| |------Set 3-------|
#define BSLICE_LUMA_ONE_CONTEXT     154, 196, 167, 167,  154, 152, 167, 182,  182, 134, 149, 136,  153, 121, 136, 122
#define PSLICE_LUMA_ONE_CONTEXT     154, 196, 196, 167,  154, 152, 167, 182,  182, 134, 149, 136,  153, 121, 136, 137
#define ISLICE_LUMA_ONE_CONTEXT     140,  92, 137, 138,  140, 152, 138, 139,  153,  74, 149,  92,  139, 107, 122, 152


#define BSLICE_LUMA_ABS_CONTEXT     107,                 167,                  91,                 107
#define PSLICE_LUMA_ABS_CONTEXT     107,                 167,                  91,                 122
#define ISLICE_LUMA_ABS_CONTEXT     138,                 153,                 136,                 167

//                                 |------Set 4-------| |------Set 5-------|
#define BSLICE_CHROMA_ONE_CONTEXT   169, 208, 166, 167,  154, 152, 167, 182
#define PSLICE_CHROMA_ONE_CONTEXT   169, 194, 166, 167,  154, 167, 137, 182
#define ISLICE_CHROMA_ONE_CONTEXT   140, 179, 166, 182,  140, 227, 122, 197

#define BSLICE_CHROMA_ABS_CONTEXT   107,                 167
#define PSLICE_CHROMA_ABS_CONTEXT   107,                 167
#define ISLICE_CHROMA_ABS_CONTEXT   152,                 152
#endif

//------------------------------------------------

static const UChar
INIT_ONE_FLAG[NUMBER_OF_SLICE_TYPES][NUM_ONE_FLAG_CTX] =
{
  { BSLICE_LUMA_ONE_CONTEXT, BSLICE_CHROMA_ONE_CONTEXT },
  { PSLICE_LUMA_ONE_CONTEXT, PSLICE_CHROMA_ONE_CONTEXT },
  { ISLICE_LUMA_ONE_CONTEXT, ISLICE_CHROMA_ONE_CONTEXT },
};

#if !VCEG_AZ07_CTX_RESIDUALCODING 
static const UChar
INIT_ABS_FLAG[NUMBER_OF_SLICE_TYPES][NUM_ABS_FLAG_CTX] =
{
  { BSLICE_LUMA_ABS_CONTEXT, BSLICE_CHROMA_ABS_CONTEXT },
  { PSLICE_LUMA_ABS_CONTEXT, PSLICE_CHROMA_ABS_CONTEXT },
  { ISLICE_LUMA_ABS_CONTEXT, ISLICE_CHROMA_ABS_CONTEXT },
};
#endif

//--------------------------------------------------------------------------------------------------

static const UChar
INIT_MVP_IDX[NUMBER_OF_SLICE_TYPES][NUM_MVP_IDX_CTX] =
{
  { 168, },
  { 168, },
  { CNU, },
};

#if CNN_BASED_INLOOP_FILTER
static const UChar
INIT_CNNLF_FLAG[NUMBER_OF_SLICE_TYPES][NUM_CNNLF_FLAG_CTX] =
{
	{ 153, },
	{ 153, },
	{ 153, },
};
#endif

static const UChar
INIT_SAO_MERGE_FLAG[NUMBER_OF_SLICE_TYPES][NUM_SAO_MERGE_FLAG_CTX] =
{
  { 153,  },
  { 153,  },
  { 153,  },
};

static const UChar
INIT_SAO_TYPE_IDX[NUMBER_OF_SLICE_TYPES][NUM_SAO_TYPE_IDX_CTX] =
{
  { 160, },
  { 185, },
  { 200, },
};

static const UChar
INIT_TRANS_SUBDIV_FLAG[NUMBER_OF_SLICE_TYPES][NUM_TRANS_SUBDIV_FLAG_CTX] =
{
  { 224,  167,  122,
#if COM16_C806_T64
  122
#endif
  },
  { 124,  138,   94,
#if COM16_C806_T64
  94
#endif
   },
  { 153,  138,  138,
#if COM16_C806_T64
  138
#endif
   },
};

static const UChar
INIT_TRANSFORMSKIP_FLAG[NUMBER_OF_SLICE_TYPES][2*NUM_TRANSFORMSKIP_FLAG_CTX] =
{
  { 139,  139},
  { 139,  139},
  { 139,  139},
};

static const UChar
INIT_EXPLICIT_RDPCM_FLAG[NUMBER_OF_SLICE_TYPES][2*NUM_EXPLICIT_RDPCM_FLAG_CTX] =
{
  {139, 139},
  {139, 139},
  {CNU, CNU}
};

static const UChar
INIT_EXPLICIT_RDPCM_DIR[NUMBER_OF_SLICE_TYPES][2*NUM_EXPLICIT_RDPCM_DIR_CTX] =
{
  {139, 139},
  {139, 139},
  {CNU, CNU}
};

static const UChar
INIT_CROSS_COMPONENT_PREDICTION[NUMBER_OF_SLICE_TYPES][NUM_CROSS_COMPONENT_PREDICTION_CTX] =
{
  { 154, 154, 154, 154, 154, 154, 154, 154, 154, 154 },
  { 154, 154, 154, 154, 154, 154, 154, 154, 154, 154 },
  { 154, 154, 154, 154, 154, 154, 154, 154, 154, 154 },
};

#if ALF_HM3_REFACTOR
static const UChar
  INIT_ALF_CTRL_FLAG[NUMBER_OF_SLICE_TYPES][NUM_ALF_CTRL_FLAG_CTX] =
{
#if JVET_C0038_GALF 
  { CNU, },
  { CNU, },
  { CNU, },
#else
  { CNU, CNU, CNU },
  { CNU, CNU, CNU },
  { CNU, CNU, CNU },
#endif
};

// initial probability for ALF flag
static const UChar
  INIT_ALF_FLAG[NUMBER_OF_SLICE_TYPES][NUM_ALF_FLAG_CTX] =
{
  {240},
  {224},
  {224}
};

// initial probability for ALF side information (unsigned)
static const UChar
  INIT_ALF_UVLC[NUMBER_OF_SLICE_TYPES][NUM_ALF_UVLC_CTX] =
{
  {154, 140},
  {140, 110},
  {139, 139},
};

// initial probability for ALF side information (signed)
static const UChar
  INIT_ALF_SVLC[NUMBER_OF_SLICE_TYPES][NUM_ALF_SVLC_CTX] =
{
  { 185, 185, CNU },
  { CNU, CNU, CNU },
  { CNU, CNU, CNU },
};
#endif

#if COM16_C806_EMT
static const UChar 
INIT_EMT_TU_IDX[NUMBER_OF_SLICE_TYPES][NUM_EMT_TU_IDX_CTX] =  
{
  { CNU,  CNU, CNU,  CNU }, 
  { CNU,  CNU, CNU,  CNU }, 
  { CNU,  CNU, CNU,  CNU }, 
};

static const UChar 
INIT_EMT_CU_FLAG[NUMBER_OF_SLICE_TYPES][NUM_EMT_CU_FLAG_CTX] = 
{
#if COM16_C806_LARGE_CTU
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU, CNU,  CNU }, 
#else
  { CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU },
  { CNU,  CNU, CNU,  CNU }, 
#endif
};
#endif

#if RRC_NS
static const UChar
INIT_NS_FLAG[3][NUM_NS_FLAG_CTX] =
{
	{ CNU },
	{ CNU },
	{ CNU },
};
#endif

#if RRC_INTRA_SHARPENING
static const UChar
INIT_SHARP_FLAG[3][NUM_SHARP_FLAG_CTX] =
{
  { CNU },
  { CNU },
  { CNU },
};
#endif
#if VCEG_AZ08_KLT_COMMON
static const UChar
INIT_KLT_FLAG[3][2 * NUM_KLT_FLAG_CTX] =
{
  { 139, 139 },
  { 139, 139 },
  { 139, 139 },
};
#endif
#if RRC_SELECTABLE_IF
static const UChar
INIT_ADAPTIVE_IF_INDEX[3][NUM_SELECTABLE_IF_CTX] =
{
  { CNU, },
  { CNU, },
  { CNU, },
};
#endif
//! \}

#endif
