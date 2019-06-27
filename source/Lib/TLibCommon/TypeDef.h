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

/** \file     TypeDef.h
    \brief    Define macros, basic types, new types and enumerations
*/

#ifndef __TYPEDEF__
#define __TYPEDEF__

#ifndef __COMMONDEF__
#error Include CommonDef.h not TypeDef.h
#endif

#include <vector>

//! \ingroup TLibCommon
//! \{
///////////////////////////////////////////////////////////
// KTA tools section start
///////////////////////////////////////////////////////////
//linjinap2
//linjianp

#define JEM_ENABLE                                        1
#if JEM_ENABLE  

#define JVET_F0096_BILATERAL_FILTER                       1   // for bitexact implementation with division see JVET-F0096
#define JVET_F0031_RMV_REDUNDANT_TRSKIP                   1

#define JVET_E0062_MULTI_DMS                              1   ///< Extended chroma multiple DM modes

#define JVET_E0077_ENHANCED_LM                            1   ///< Enhanced LM mode
#if JVET_E0077_ENHANCED_LM
#define JVET_E0077_MMLM                                   1
#define JVET_E0077_LM_MF                                  1
#endif

#define JVET_E0052_DMVR                                   1 //Decoder-side motion vector refinement based on bilateral template matching
#if JVET_E0052_DMVR
#define DMVR_HALF_ME                                      0
#endif

#define JVET_D0134_PSNR                                   1
#define JVET_D0135_PARAMS                                 1
#define JVET_D0186_PRECISEPSNR                            1
#define JVET_F0064_MSSSIM                                 1 // Calculate MS-SSIM scores

#define JVET_D0033_ADAPTIVE_CLIPPING                      1
#if JVET_D0033_ADAPTIVE_CLIPPING
#define JVET_D0033_ADAPTIVE_CLIPPING_ENC_METHOD           1
#endif

#define JVET_D0123_ME_CTX_LUT_BITS                        1

#define JVET_C0024_QTBT                                   1

#if JVET_C0024_QTBT
#define JVET_C0024_BT_FIX_TICKET22                        1
#define JVET_C0024_QTBT_FIX_QUANT_TICKET25                1  // JVET-D0025
#define RDOQ_BIT_ESTIMATE_FIX_TICKET29                    1  // correct RDOQ uninitialized values in case of vertical scan 
#define FIX_TICKET32                                      1  ///< Fix of ticket #32 (Initialization for INIT_LAST)

#define MIN_CU_LOG2                                       2
#if MIN_CU_LOG2==1
#define JVET_C0024_DF_MODIFY                              1 //deblocking modifications
#else
#define JVET_C0024_DF_MODIFY                              0
#endif

#define JVET_C0024_BT_RMV_REDUNDANT                       1  ///< Remove redundant BT structure for B/P slice

#define JVET_C0024_SPS_MAX_BT_SIZE                        0  ///< signal max BT size in SPS
#define JVET_C0024_SPS_MAX_BT_DEPTH                       1  ///< signal max BT depth in SPS 

#define JVET_C0024_CTU_256                                0  ///< support CTU 256 for QTBT, force QT split for CU 256x256 
#define JVET_C0024_ENCODER_OVERFLOW_FIX                   1  ///< fix the encoder overflow in the case of very high QP
#define JVET_C0024_DELTA_QP_FIX                           1  ///< support delta QP signaling in QTBT

// for fast algorithms
#define JVET_C0024_AMAX_BT                                1  ///< slice level adaptive maximum BT size (encoder only)
#define JVET_C0024_AMAX_BT_FIX                            1  ///< Support parallel encoding
#define JVET_C0024_AMAX_BT_FIX_TICKET23                   1
#define JVET_C0024_FAST_MRG                               1
#define JVET_C0024_PBINTRA_FAST                           1
#if JVET_C0024_PBINTRA_FAST
#define JVET_C0024_PBINTRA_FAST_FIX                       1  ///< JVET-D0077
#endif
#define JVET_C0024_ITSKIP                                 1  ///< skip zero row/column in inverse transform (decoder speedup)

#define JVET_D0077_FAST_EXT                               1  ///< extension of exsiting fast algorithm 

#define JVET_D0077_TRANSFORM_OPT                          1  ///< software optimization to take full advantages of zero rows/columns in transform coefficients
#define JVET_D0077_SAVE_LOAD_ENC_INFO                     1  ///< save and load encoder decision for speedup

#define WCG_LUMA_DQP_CM_SCALE                             1  ///< enable luma adaptive QP and chroma QP scale, intended for data in ST-2084 container

#define JVET_E0023_FAST_ENCODING_SETTING                  1
#if JVET_E0023_FAST_ENCODING_SETTING
#define PICTURE_DISTANCE_TH                               1  // If a distance between current picture and reference picture is smaller than or equal to PICTURE_DISTANCE_TH,
                                                             // FAST_SKIP_DEPTH_VALUE is used as a threshold of early CU determination. Otherwise a higher value (SKIP_DEPTH_VALUE) is used.
#define SKIP_DEPTH_VALUE                                  3
#define FAST_SKIP_DEPTH_VALUE                             2
#endif

#endif // end of JVET_C0024_QTBT

#define JVET_C0046_OMIT_ASSERT_ERDPCM                     1  ///< for RExt, omit assertion related to Explict Residual DPCM

#define JVET_C0046_FIX_DECODER_DEBUG_BIT_STATISTICS       1  ///< bugfix on DECODER_DEBUG_BIT_STATISTICS for QTBT

#define VCEG_AZ08_USE_KLT                                 1  ///< KLT (if defined 1, use cfg option of KLT to control the enablement of intra KLT and inter KLT (INTERA_KLT, VCEG_AZ08_INTER_KLT should be set as 1); if 0, use INTERA_KLT, VCEG_AZ08_INTER_KLT to control the enablement.)

#define VCEG_AZ08_INTER_KLT                               1  ///< (default 1) Enable inter KLT
#define VCEG_AZ08_INTRA_KLT                               1  ///< (default 1) Enable intra KLT
#define VCEG_AZ08_KLT_COMMON                              (VCEG_AZ08_INTER_KLT || VCEG_AZ08_INTRA_KLT)

#if VCEG_AZ08_KLT_COMMON
#define VCEG_AZ08_FORCE_USE_GIVENNUM_BASIS                1  /// (default 1) If defined, force to use up to FORCE_BASIS_NUM basis to reduce complexity.
#define VCEG_AZ08_USE_SSD_DISTANCE                        0  ///< (default 0) If defined, use SSD distance.
#define VCEG_AZ08_USE_SAD_DISTANCE                        1  ///< (default 1) If defined, use SAD distance.
//Speed up
#define VCEG_AZ08_FAST_DERIVE_KLT                         1  ///< (default 1) If defined, will use fast algorithm to calculate KLT basis
#define VCEG_AZ08_USE_SSE_SPEEDUP                         0  ///< (default 0) If defined, will use sse for speeding up (Note: should use x64 compile mode)
#if VCEG_AZ08_USE_SSE_SPEEDUP
#define VCEG_AZ08_USE_SSE_SCLAE                           1  ///< (default 1) If defined, will use SSE for calculating the scaling of float to get integer KLT basis
#define VCEG_AZ08_USE_FLOATXSHORT_SSE                     1  ///< (default 1) If defined, will use float*short
#define VCEG_AZ08_USE_SHORTXSHORT_SSE                     1  ///< (default 1) If defined, will use SSE to calculate short x short multiplication
#if VCEG_AZ08_USE_SHORTXSHORT_SSE
#define VCEG_AZ08_USE_TRANSPOSE_CANDDIATEARRAY            1  ///< (default 1) If defined, will use transpose of candidate array to facilitate the vector multiplication
#endif
#if VCEG_AZ08_USE_SAD_DISTANCE
#define VCEG_AZ08_USE_SSE_BLK_SAD                         1  ///< (default 1) If defined, will calculate the block difference use intel SSE 
#define VCEG_AZ08_USE_SSE_TMP_SAD                         1  ///< (default 1) If defined, will calculate the template difference use intel SSE 
#endif
#endif

#define VCEG_AZ08_INTER_KLT_MV_BUGFIXED                   1  /// Fixed the bug related with mv (2016.2.29)
#endif


#define ALF_HM3_REFACTOR                                  1  ///< Adaptive loop filter (ALF)
#if ALF_HM3_REFACTOR
#define COM16_C806_ALF_TEMPPRED_NUM                       6  ///< 0: no temporal prediction
#if COM16_C806_ALF_TEMPPRED_NUM
#define FIX_TICKET12                                      1  ///< fixed ticket #12
#define JVET_E0104_ALF_TEMP_SCALABILITY                   1  ///< ALF temporal prediction with temporal scalability as in JVET-E0104
#if JVET_E0104_ALF_TEMP_SCALABILITY
#define JVET_E0104_ALF_MAX_TEMPLAYERID                    5  ///< maximum number of temporal layers
#endif
#endif

#define JVET_C0038_GALF                                   1 ///<JVET-C0038 GALF
#if JVET_C0038_GALF
#define JVET_C0038_SHIFT_VAL_HALFW                        1  ///<clean up
#define JVET_C0038_NO_PREV_FILTERS                        16 ///<number of fixed filters per class
#endif
#endif


#define COM16_C806_VCEG_AZ10_SUB_PU_TMVP                  1  ///< Sub-PU level motion vector prediction
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP                     
#define COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL         1
#define COM16_C806_DISABLE_4X4_PU                         1
#define COM16_C806_GEN_MRG_IMPROVEMENT                    1
#define JVET_C0035_ATMVP_SIMPLIFICATION                   1 ///< JVET-C0035 ATMVP_SIMPLIFICATION
#endif

#define COM16_C806_OBMC                                   1  ///< Overlapped block motion compensation (OBMC)
#if COM16_C806_OBMC
#define COM16_C806_AOBMC_MAXCUSIZE                        16 ///< Maximum CU size which can apply OBMC adaptively, larger CUs always apply OBMC
#endif

#define COM16_C806_EMT                                    1  ///< Explicit multiple core transform
#if COM16_C806_EMT

// reducing complexity stage-1
#define FAST_EMT_v1                                       1//1123,fast EMT v1
//#define EMT_SIGNAL                                        1 //1121
//#define EMT_SIGNAL_INTER                                  1 //1121
//#define EMT_MOD_TR                                        1 //1121
#endif

#define COM16_C806_T64                                    1  ///< 64x64 transform

#if COM16_C806_T64
#define JVET_C0046_ZO_ASSERT                              1  ///< assertion on last coeff and coded_sbk_flag when zeroing out is used (no TS and no TQBypass and using large transform is satisfied)
#if JVET_C0046_ZO_ASSERT
#define JVET_C0046_ZO_ASSERT_CODED_SBK_FLAG               1  ///< if (iCGX > TH1 || iCGY > TH1) and (no TS && no TQBypass), then coded_sbk_flag(iCGX, iCGY) shall be 0.
#define JVET_C0046_ZO_ASSERT_LAST_COEF                    1  ///< if (posLastX>TH2 || posLastY>TH2) and (no TS && no TQBypass), then last coef (x,y) shall be in the low frequency domain.
#endif
#endif

#if COM16_C806_EMT || COM16_C806_T64
#define COM16_C806_TRANS_PREC                             2  ///< Integer transform matrix precision
#endif

#define COM16_C806_LARGE_CTU                              1  ///< Large CTU up to 256x256

#define COM16_C806_LMCHROMA                               1  ///< Cross-component linear model (CCLM) prediction
#if COM16_C806_LMCHROMA
#define COM16_C806_CR_FROM_CB_LAMBDA_ADJUSTMENT           1
#endif

#define VCEG_AZ07_IMV                                     1  ///< Locally adaptive motion vector resolution (AMVR)

#if VCEG_AZ07_IMV
#define JVET_E0076_MULTI_PEL_MVD                          1
#endif

#define VCEG_AZ07_FRUC_MERGE                              1  ///< Pattern matched motion vector derivation
#if VCEG_AZ07_FRUC_MERGE
#define JVET_E0060_FRUC_CAND                              1  ///< E0060: Add candidates to FRUC lists, and adjust the number of these candidates
#define JVET_F0032_UNI_BI_SELECTION                       1  ///< JVET_F0032: selection between uni-prediction and bi-prediction for FRUC template matching mode
#endif

#define JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC         1

#if JVET_B058_HIGH_PRECISION_MOTION_VECTOR_MC
#define VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE          2   ///< additional precision bit for MV storage
#elif VCEG_AZ07_FRUC_MERGE 
#define VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE          1   ///< additional precision bit for MV storage
#endif

#define VCEG_AZ06_IC                                      1  ///< Local illumination compensation (LIC)
#if VCEG_AZ06_IC
#define VCEG_AZ06_IC_SPEEDUP                              1  ///< speedup of IC
#if JVET_C0024_QTBT
#undef VCEG_AZ06_IC_SPEEDUP                              
#define VCEG_AZ06_IC_SPEEDUP                              0  ///< speedup of IC
#define IC_THRESHOLD                                      0.03
#endif
#endif

#define VCEG_AZ07_INTRA_4TAP_FILTER                       1  ///< 4-tap interpolation filter for intra prediction
#define VCEG_AZ07_INTRA_BOUNDARY_FILTER                   1  ///< Addtional boundary filter or intra prediction
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER                          
#define VCEG_AZ07_INTRA_BOUNDARY_FILTER_MULTI_LINE        1  ///< 0: Filter one boundary line, 1: Filter 4 boundary lines
#endif                                                       

#define VCEG_AZ07_INTRA_65ANG_MODES                       1  ///< 65 intra prediction directions
#if VCEG_AZ07_INTRA_65ANG_MODES
#define JVET_B0051_NON_MPM_MODE                           1  // Use two mode sets for non-MPM mode coding
#endif
#define JVET_C0055_INTRA_MPM                              1  ///< Intra MPM derivation from JVET-C0055

#define VCEG_AZ07_ECABAC                                  1  ///< CABAC improvements
#if VCEG_AZ07_ECABAC                                         
#define VCEG_AZ07_CTX_RESIDUALCODING                      1  ///< Context modeling for transform coefficient levels
#define VCEG_AZ07_BAC_ADAPT_WDOW                          1  ///< Multi-hypothesis probability estimation
#define VCEG_AZ07_INIT_PREVFRAME                          1  ///< Initialization for context models
#if VCEG_AZ07_INIT_PREVFRAME
#define VCEG_AZ07_INIT_PREVFRAME_FIX                      1  ///< Fix for GOP16
#endif
#endif                                                       

#define VCEG_AZ05_MULTI_PARAM_CABAC                       1  ///< CABAC probability estimation with 2 windows 
#define VCEG_AZ05_BIO                                     1  ///< Bi-directional optical flow (BIO)
#define VCEG_AZ05_INTRA_MPI                               0  ///< Multi-parameter Intra prediction
#define VCEG_AZ05_ROT_TR                                  0  

#if VCEG_AZ05_BIO                                            
#define COM16_C1045_BIO_HARMO_IMPROV                      1  ///< Improvement of BIO
#define JVET_C0027_BIO                                    1   /// MV refinement max value up, BIO_LDB check optimization,  BIO  for 1/16 pel MV support
#define JVET_F0028_BIO_NO_BLOCK_EXTENTION                 1
#define JVET_G0082                                        1
#endif                                                       

#define COM16_C1016_AFFINE                                1  ///< Affine motion prediction
#if COM16_C1016_AFFINE
#define JVET_B0038_AFFINE_HARMONIZATION                   1  ///< Harmonization of affine, OBMC and DBF
#define JVET_C0025_AFFINE_FILTER_SIMPLIFICATION           1  ///< Simplification of MC filters for affine
#if JVET_C0025_AFFINE_FILTER_SIMPLIFICATION && VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE != 2
#error JVET_C0025_AFFINE_FILTER_SIMPLIFICATION shall be off if no 1/16 MV accuracy
#endif
#endif

#define JVET_G0104_PLANAR_PDPC                            1  ///< apply PDPC for planar mode

#define COM16_C1044_NSST                                  1  ///< Mode dependent non-separable secondary transforms
#if COM16_C1044_NSST

// reducing complexity stage-1
#define FAST_NSST_v1                                                1 //1123,NSST fast version1
//#define NSST_PLANAR                                                 1 // 1121
//#define NSST_TU_44                                                  1 // replace 2 4x4 NSST of 8x4,4X8 with 1 4x4 NSST on left-top block //1121
#endif
#if COM16_C1044_NSST || VCEG_AZ05_ROT_TR
#define JVET_C0042_UNIFIED_BINARIZATION 1       // unified binarization for NSST index and context modeling based on partition type and intra prediction mode
#endif
#if COM16_C1044_NSST && VCEG_AZ05_ROT_TR                     
#error                                                       
#endif                                                       

#if COM16_C1044_NSST && JVET_C0024_QTBT
#define QTBT_NSST                                         1
#endif

#if COM16_C1044_NSST
#define JVET_C0045_C0053_NO_NSST_FOR_TS                   1  ///< JVET-C0045/C0053: Disable NSST for TS coded blocks 
#define JVET_D0120_NSST_IMPROV                            1  ///< JVET-D0120: NSST improvements using HyGT and 8x8 NSST  
#define JVET_D0127_REDUNDANCY_REMOVAL                     1
#endif

#define COM16_C1046_PDPC_INTRA                            1  ///< Position dependent intra prediction combination

#if COM16_C1046_PDPC_INTRA && VCEG_AZ05_INTRA_MPI
#error
#endif

#define PARALLEL_ENCODING_SAO_FIX                         1  ///< Fix of SAO for parallel encoding proposed in JVET-B0036
#define PARALLEL_ENCODING_RAS_CABAC_INIT_PRESENT          1  ///< Fix of CABAC initialization for parallel encoding proposed in JVET-B0036

// encoder only changes
#define COM16_C806_SIMD_OPT                               1  ///< SIMD optimization, no impact on RD performance

#define JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY             1 ///< This approach derives lambda from QP+QPoffset+QPoffset2. QPoffset2 is derived from QP+QPoffset using a linear model that is clipped between 0 and 3.
                                                            // To use this capability enable config parameter LambdaFromQpEnable
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
#define JVET_B0039_QP_FIX                                 0  ///< Recalcualtes QP to align with a HM lambda (same relation as for all intra coding is used)
#else
#define JVET_B0039_QP_FIX                                 1  ///< Recalcualtes QP to align with a HM lambda (same relation as for all intra coding is used)
#endif

#endif //JEM_ENABLE
///////////////////////////////////////////////////////////
// KTA tools section end
///////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////
// Huawei start section
///////////////////////////////////////////////////////////
#define HEM_TOOLS_ON                                      1
#if HEM_TOOLS_ON
#define GuoLei                                            0
// Deep learning based tools SPS
#define BARC_INTER_SPS                                    1
#define BARC_INTRA_CHROMA_SPS															1
#define INTERPOLATE_FRAME_SPS                             1
// Deep learning based tools SPS

#define BLOCK_BASED_BACKGROUND_REFERENCE                  1   //block based background reference
#if BLOCK_BASED_BACKGROUND_REFERENCE
#define CLOSE_AFFINE_FOR_BG_REFERENCE                     1
#define FAST_ME_FOR_BG_REFERENCE                          1
#define DISABLE_DQP_SIGNALING                             0
#define QTBT_DQP_FIX                                      1
#endif

#define USE_CONVLUTIONAL_NEURAL_NETWORK                   1
#define USE_TENSORFLOW                                    0
#define USE_LIBTORCH                                      1
#define TEST_LIBTORCH                                     0
#define frameworkCombine                                  1
#if USE_CONVLUTIONAL_NEURAL_NETWORK
#define CNN_BASED_INLOOP_FILTER                           1 
#endif
#define USE_NEW_MODEL                                     1
#if USE_CONVLUTIONAL_NEURAL_NETWORK
#define DEFAULT_DL_FRAMEWORK                              "Caffe"
#define DEFAULT_WORKING_MODE                              "CPU"
#define DEFAULT_CAFFE_MODEL_PATH                          "./CaffeModel/"
#define DEFAULT_GPU_ID                                    0
#endif

#if CNN_BASED_INLOOP_FILTER
#define CNN_BASED_INLOOP_FILTER_SPS                        1
#define PRINT_INFO                                         1
#define DEFAULT_LOOP_FILTER_CAFFEMODEL_FILENAME_TEMPLATE  "QP??.caffemodel"
#define DEFAULT_LOOP_FILTER_NET_PROTO_FILENAME            "deploy.prototxt"    
#define DEFAULT_CAFFEMODEL_OUTPUT_LAYER_NAME              "res"
#define DEFAULT_INPUT_BLOCK_SIZE                          140 
#endif

#define AFFINE_IMPR                                       1
#if AFFINE_IMPR
#define AFFINE_IMPR_V1                                    1
#if AFFINE_IMPR_V1
#if COM16_C1016_AFFINE
#define JVECT_C0062_AFFINE_SIX_PARAM                      1 ///< enable 6 parameter affine model
#define AFFINE_MERGE_IMPROVE                              1
#define AFFINE_AMVP_IMPROVE                               1
#endif
#define MRG_CONSIDER_COMPLEX                              1
#endif

#define AFFINE_IMPR_V2                                    1
#if AFFINE_IMPR_V2
#define AFFINE_CLIP_REMOVE                                1 ///< remove clip in affine motion vector
#define AFFINE_UNIFY_MC                                   1 ///< unify affine store and MC motion vector
#define AFFINE_MC_PRECISION_FIX                           1 ///< fix affine MC precision
#endif

// AFFINE ME related
#define AFFINE_INHERIT_FOUR_PARAM                         1
#define AFFINE_REDUCE_ITERATION                           1
#define AFF_SKIP_ME                                       1
#define AFF_EARLY_TERMINATION                             1
#endif

#define RRC_NS                                            1
#if RRC_NS                                          
namespace ns
{
	enum ClassifierType
	{
		BLOCK_VARIANCE
	};
}

#define RRC_NS_USE_SSE                                  1
#define RRC_NS_PEL_TYPE                                 short
#define RRC_NS_TABLE_DISTANCE_TYPE                      int
//#define RRC_NS_BM_ALGO                                  ALG2
#define RRC_NS_BM_ALGO                                  SB_2
#define RRC_NS_LIFTING                                  1
#define RRC_NS_ITER                                     0
#define RRC_NS_BM_SB_SIZE                               12
#define RRC_NS_NHW                                      8
#define RRC_NS_nHW                                      18
#define RRC_NS_USE_SAD                                  0
#define RRC_NS_TLAYER                                   MAX_UINT // Control which temporal layers to filter. Layers starting from RRC_NS_TLAYER are not filtered
#define RRC_NS_CHROMA                                   1 
// Debug modes
// Attention!
// This defines can make the steram undecodable
#define RRC_NS_OPTIMIZE_OFF                             0
#define RRC_NS_FORCE_FILTERING_WHOLE_FRAME              0   
#define RRC_NS_DISABLE_SIGNALING                        0   
#endif // RRC_NS

#define RRC_INTRA_SHARPENING                              0  ///< Enable Intra Sharpening tool
#if RRC_INTRA_SHARPENING
#define INTRA_SHARP_CHROMA                                1  ///< Intra Sharpening for chroma components
#if COM16_C1046_PDPC_INTRA
#define INTRA_SHARP_PDPC_HARM                             1  ///< Harmonization between Intra Sharpening and PDPC
#endif
#if VCEG_AZ07_INTRA_4TAP_FILTER
#define INTRA_SHARP_4TAP_HARM                             1  ///< Harmonization between Intra Sharpening and 4tap interpolation filter
#define INTRA_SHARP_CUBIC_CHROMA                          1  ///< Cubic filter for chroma if sharpened
#endif
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
#define INTRA_SHARP_INTRA_BOUNDARY_HARM                   0  ///< Harmonization between Intra Sharpening and Intra Boundary Filter (not fully implemented)
#endif
#define INTRA_SHARP_MAX_SIZE                              32 ///< Limit max block size
#define INTRA_SHARP_MAX_ASPECT_RATIO                      2  ///< Limit aspect ratio
#define INTRA_SHARP_D0077                                 0
#endif//RRC_INTRA_SHARPENING

#define RRC_DWDIP                                         0  ///< Enable Distance-Weighted Directional Intra Prediction
#if RRC_DWDIP

#define RRC_DWDIP_MIN_PU_SIZE                             16 
#define RRC_DWDIP_INIT_W_GRAD                             (0.75)
#define RRC_DWDIP_2Nx2N_CHECK_DISABLE                     (JVET_C0024_QTBT)
#define RRC_DWDIP_ASPECT_RATIO_THRESHOLD                  8

#define RRC_INTRA_PRED_DATACU_ABSPARTIDX                  1

#define RRC_DWDIP_MLR                                     1
#define RRC_DWDIP_MLR_IMPR                                1
#define RRC_DWDIP_DIVISION_FREE_ADDSHIFT                  1

#if RRC_DWDIP_MLR_IMPR  
#define BLK_LENGTH_VALUES                                 (MAX_LOG2_TU_SIZE_PLUS_ONE - MIN_CU_LOG2)
#define CLIP_FACTOR                                       (0.3)
#endif

#endif //RRC_DWDIP

#if JVET_E0052_DMVR
#define DMVR_IMP_PHASE1                                   1  /// DMVR improvements by Semih - Max number of iteration = 7
#if DMVR_IMP_PHASE1
#undef DMVR_HALF_ME
#define DMVR_HALF_ME                                      1
#define DMVR_IMP_PHASE2                                   1
#if DMVR_IMP_PHASE2
#define MAX_NUM_DMVR_ITERATION                            7
#endif //DMVR_IMP_PHASE2
#endif //DMVR_IMP_PHASE1
#endif //JVET_E0052_DMVR

#define HEM_2_0_MIMICS                                    1 //Temporary macro to ensure HEM 2.1 has same BD-Rate as HEM 2.0
#define INTRA_LOOPS_ORDER_FIX                             1
#define TT_BUGFIX_ZIDX                                    1

#define DIMD_TOOL                                         0
#if DIMD_TOOL
#define DIMD_INTRA_PRED                                   1 // Enable DIMD Intra prediction
#define DIMD_INTRA_SHARP_HARM                             1 // When DIMD is enabled, disable Intra Sharpening
#define DIMD_DWDIP_HARM                                   1 // When DIMD is enabled, disable DWDIP
#define DIMD_SM_TEMP_SIZE                                 2 // Small block template size i.e. (W, H <= 8)
#define DIMD_LG_TEMP_SIZE                                 3 // Large block template size i.e. (W, H > 8)
#define DIMD_NUM_INTRA_DIR_INC                            1 // Increase the number of angular modes for DIMD mode calculation
#define DIMD_MAX_TEMP_SIZE                                4 // Maximum DIMD template size 
#define DIMD_MAX_BT_DEPTH_SIG_AI                          1 // Used only for I Slices. Maximum depth of the tree, where DIMD can be used 
#define DIMD_DISABLE_CHROMA_P_B_SLICE_JB                  1 // Chroma in P,B slices do not use DIMD derived intra mode
#define DIMD_PARSING_FIX_INTRA_JB                         1 // DIMD parsing fix for Intra blocks w.r.t to chomra intra mode derivation
#define DIMD_REDUCE_LINE_MEMORY_JMB                       0 // DIMD reduce line memory at CTU line borders by only using left template
#define DIMD_ADAPTIVE_TS_BLOCK_BASED_JMB                  1 // Make template size (TS) adaptive based on the actual block size (have one more level of distinction for TS). 
                                                            // Also takes care of Line memory access (limited to max. 4 lines as other tools such as MLIP and DBF).
#endif // end DIMD_TOOL

#define MLIP_EXTEND_REF_LINE                              0 /// MLIP
#if MLIP_EXTEND_REF_LINE
#define MLIP_REF_LINE_NUM                                 4
#define MLIP_IMPROVED_INDEX_CODING                        1 ///MX: improveing the method coding ext ref index; 
#if MLIP_IMPROVED_INDEX_CODING
#define MLIP_MODI_INDEX                                   3 /// org: 0 1 2 3; 3: 0 3 2 1;  2: 0 2 1 3
#endif
#define MLIP_SKIP_NO_USE_PROCESS                          1 /// skip some useless process to speedup
#undef JVET_D0127_REDUNDANCY_REMOVAL
#define JVET_D0127_REDUNDANCY_REMOVAL                     0
#define HAMONIZATION_FOR_MLIP_AND_D0127                   1 // if use Extend ref line, use anohter method like D0127; current implementation only work when NSST/ROT not use for further reference line
#define HAMONZATION_MLIP_INTRA4TAPOFF                     1 // harmonization MLIP and intra 4taps turn off. can not turn off.
#endif

//Triple Tree
#define TT_DEV                                            1
#if TT_DEV
#define TT_TREE_DS1                                       1 //data structure of TT & BT, the same to D0117
#define TT_TREE_DS3                                       1 //fix on motion vector field setting at decoder
#define TT_TREE_DS5                                       1 //deriveZIdxInCTU and deriveZIdxInCU for one PU

#define TT_SYNTAX                                         1 //tree syntax encoding & decoding, the same to D0117
#define TT_SYNTAX_NEW                                     1
#define TT_ENC_OPT                                        1 //encoder related
#define TT_DEC_OPT                                        1 //decoder related

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP                     
#define TT_TREE_ATMVP                                     1  //tightly related to ATMVP harmonization with TT
#endif
#if COM16_C806_OBMC
#define TT_TREE_OBMC                                      1  //tightly related to OBMC harmonization with TT
#endif
#if VCEG_AZ07_FRUC_MERGE
#define TT_TREE_FRUC                                      1  //tightly related to FRUC harmonization with TT
#endif
#if COM16_C1016_AFFINE
#define TT_TREE_AFFINE                                    1  //tightly related to affine harmonization with TT
#endif
#if DIMD_TOOL
#define TT_TREE_DIMD                                      1 // Triple Tree harmonization for DIMD 
#endif

//optimization - normative
#define TT_RMV_PART1_BTBT                                 1 //forbid BT+(BT+BT) in the same dir
#define TT_LIMIT_SIZEIDX                                  1 //ban tt based on deltaSizeIdx
#define TT_LIMIT_AR                                       8 //ban tt/bt based on aspect ratio 1:N, N can be 4 or 8

#define TT_LIMIT_SMALL                                    1 //TT small part cannot use same dir
#define TT_LIMIT_LARGE                                    1 //TT large part cannot use same dir
#if TT_LIMIT_SMALL || TT_LIMIT_LARGE
#define TT_LIMIT_SYNTAX                                   1 //1: change syntax, 0: just encoder fast algorithm
#endif
//specific rule on TT
#define TT_LEAF_NODE                                      0 //tt use as the leaf node

//fast algorithm - encoder only
#define TT_FAST_BAN_TT                                    1
#if TT_FAST_BAN_TT
#define TT_FAST_BAN_TT_TH                                 1 //2: 0.1% loss, 12% Enc
#endif
#define TT_FAST_USE_ONE_TT                                1

#define TT_FAST_SMALL_ES                                  1 //for tt small CU, refer to split info of the same directional BT
#define TT_FAST_LARGE_ES                                  1 //for tt large CU, refer to split info of the same directional BT

//new RDO options, can be applied to QTBT also
#define ADAPTIVE_RD_ORDER                                 1 //adaptive RDO order
#define ENC_FAST_SKIP_RDC_EQUAL                           1 //skip one BT/TT if FULL & the other BT/TT have the same distortion

#define JVET_D0077_SAVE_LOAD_ENC_INFO_EXT                 1  ///< extension by storing and retrieving a cu coding mode info; this macro marks te major changes (cannot switch off)
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
#define MAX_NUM_ROUTE                                     24  ///< max number of routes to save, which shall be no less than 16 -> moved to BT_FORCE part - Han
#define FAST_ADDITIONAL_REF_ROUTE                         1  ///< route matching criterion is loosen by not checking left and up neighbors' BT/TT Depth
//todo: save and load split info is not finished

//for debug
#define DEBUG_TOOLS_OFF                                   0 //switch off JEM tools for debug
#endif


/* to do list
1. check deltaQP issue: related to JVET_C0024_DELTA_QP_FIX
2. Now, only support MaxBTDepth <=3
*/

enum STSplitMode
{
  ST_NO_SPLIT = 0 * 4 + 0 * 2 + 0,
  ST_HORI_BT = 1 * 4 + 0 * 2 + 0,
  ST_VERT_BT = 2 * 4 + 0 * 2 + 0,
  ST_HORI_TT_L = 1 * 4 + 1 * 2 + 0,
  ST_HORI_TT_S = 1 * 4 + 1 * 2 + 1,
  ST_VERT_TT_L = 2 * 4 + 1 * 2 + 0,
  ST_VERT_TT_S = 2 * 4 + 1 * 2 + 1,
};

#endif

#define BT_FORCE                                          1 //Enable Forced BT
#if BT_FORCE
#define BT_FORCE_EXT                                      1 //Enable Maixmum BT depth extension for boundary CU
#define BT_FORCE_INTRA_ON                                 1 //Enalbe I frame forced BT partition
#if BT_FORCE_INTRA_ON
#define BT_FORCE_INTRA_SPEED_UP                           1 //speed up intra FBP
#endif
#if TT_TREE_DS1
#define BT_FORCE_FBTD                                     0 //fix Ctx bug by calclulting FBT/FTT depth
#define BT_FORCE_CTX0                                     0 //fix Ctx bug using 0 CTX for boundary CTU
#define BT_FORCE_CTX                                      1 //fix Ctx bug using uiBTdepht and not larger than 2
#define BT_FORCE_BF                                       1 //Extend BT flag buffer
#define BT_FORCE_TT                                       1 //Enable forced TT
#define BT_FORCE_AR_TT                                    0 // Aspect ratio limitation 1:N for forced TT, 0 means off
#endif
#define BT_FORCE_AR                                       0 // Aspect ratio limitation 1:N for forced BT, 0 means off
#if BT_FORCE_INTRA_ON&&DIMD_INTRA_PRED &&!BT_FORCE_TT
#define BT_FORCE_DIMD_HARM                                1 //increse DIMD_MAX_BT_DEPTH_SIG_AI on I slice boundary Cu (might be not necessary for forced TT case)
#endif
#if BT_FORCE_EXT
#define BT_FORCE_QTN                                      1 //calculate QT level from MaxBTSize QTN = log2(CTUSize) - log2 (MaxBTSize);
#undef MAX_NUM_ROUTE
#define MAX_NUM_ROUTE                                     88  ///< max number of routes to save, which shall be no less than 16
#else
#undef MAX_NUM_ROUTE
#define MAX_NUM_ROUTE                                     24  ///< max number of routes to save, which shall be no less than 16
#endif
#endif //end of BT_FORCE

#define RRC_FRAC_CHROMA_QP_OFFSET                         1  ///< Enable Fractional Chroma QP offset
#if RRC_FRAC_CHROMA_QP_OFFSET
#define RRC_FRAC_CHROMA_QP_OFFSET_ADD                 (0.5)  ///< Fractional Chroma QP offset value
#endif

#define THIRD_FRUC_MODE                                   1 // Switch between Uni and Bi prediction for FRUC template matching mode.
#define FIX_DECODER_STATISTICS                            1 
#if THIRD_FRUC_MODE
#define THIRD_FRUC_MODE_CTX_SWITCH_LEVEL                  8 // use 2 Ctx to code Uni-Bi prediction flag. TODO: (SEMIH) - Context initializers seem to be wrong.
#define THIRD_FRUC_MODE_DISABLE_BILATERALFORTID0          1 // do not signal Bilateral Mode if there is only 1 ref pic in L0 and L1 which is identical.
#define THIRD_FRUC_MODE_ENCODER_SPEEDUP                   1 // CU based FRUC estimation speedup
#define THIRD_FRUC_MODE_ENCODER_SPEEDUP2                  0 // Sub-PU level FRUC estimation speedup, not tested throughly, there might be decoder mismatches for some sequences (although I did not observe any)
#define THIRD_FRUC_MODE_TT_OFF_FIX                        1 // Fix compilation error when TT is turned off
#endif

// DMVR improvement
#define ITH_C2_DMVR_IMPR                                  1
#if ITH_C2_DMVR_IMPR
#define ITH_C2_PROPOSAL1                                  1 ///< Unclipped Pred + MC opt
#if ITH_C2_PROPOSAL1
#define ITH_C2_UNCLIPPED_DMVR                             1 ///< Reduce complexity without affecting BD-Rate
#if ITH_C2_UNCLIPPED_DMVR
#define ITH_C2_MC_OPT                                     1 ///< code level optimize for MC
#endif
#endif

#define ITH_C2_PROPOSAL2                                  1 ///< Do half Pel refinement outside the iteration loop and control number of half pel point to search
#if ITH_C2_PROPOSAL2
#undef MAX_NUM_DMVR_ITERATION
#define MAX_NUM_DMVR_ITERATION                           11 ///< Number of DMVR iteration
#define ITH_C2_NUM_HALF_POINT                             4 ///< Number of half pel points to consider for half pel search
#define ITH_NON_NORMATIVE_DISABLE                         0 ///< disable non-normative optimize for DMVR half MC
#endif
#endif

#if !ITH_C2_PROPOSAL2
#define ITH_C2_NUM_HALF_POINT                             0 ///< non-normative optimize for DMVR half MC
#endif

#define DMVR_MRSAD_SPEEDUP                                1
#define DMVR_TEMPLATE_SWITCH                              1  ///< Use MRSAD instead of SAD cost metric.
#define DMVR_ITERATION_TERMINATE                          1  ///< TErminate DMVR iterations earlier

#define BARC                                              1 //Enable BARC 
#if BARC
#define BARC_INTRA_LUMA_SPS                               1
#define BARC_CORE                                         1
#define BARC_Boundary                                     1 //re-define boundary for BARC CTU
#define BARC_FUNC                                         1 
#define BARC_Neighbour                                    1 //get neighbour for BARC block
#define BARC_DBK                                          1 //refine DBK for BARC CTU
#define BARC_AdpQP                                        1 
#define BARC_AdpQPSTEP                                    7 //BARC Luma QP offset
#define BARC_DScaleRatio                                  1 //Scale Ratio, 2:1 for horizontal and vertical 
#define BARCPostUpSample                                  0 //Up sample right and below boundary after right and blow reconstructed samples are available
#define BARC_REF_DN_FILT                                  1 // down sample filter
#define BARC_REF_UP_FILT                                  0 // up sample filter
#define BARC_DIMD                                         1 //Harmonization with DIMD
#define BARC_BTFORCE                                      1 //Harmonization with BTFORCE
#define BARC_MAX_CTU_BUFFER_LINE                          MAX_CU_SIZE + 14
#define BARC_MAX_CTU_BUFFER                              (MAX_CU_SIZE + 14)*(MAX_CU_SIZE + 14)
#define BARC_ONLY_FOR_INTRA                              1
#define BARC_REF_UP_SAMPLE_CNN                           1
#if BARC_REF_UP_SAMPLE_CNN
#define BARC_DEBUG                                        0
#define DEFAULT_CNN_UP_SAMPLE_MODEL_FILE_BICUBIC          "./CaffeModel/Barc_??_5_iter_10000.caffemodel"
#define DEFAULT_CNN_UP_SAMPLE_NET_FILE_BICUBIC            "./CaffeModel/netLuma.prototxt"
#define BARC_CREATE_CAFFE_MODEL_PICTURE_LEVEL             0 // Create TComCaffe object in picture level
#define BARC_UP_SAMPLE_RDO                                1
#if BARC_UP_SAMPLE_RDO
#define BARC_POST_UP_SAMPLE_ONLY_FOR_FILT                 1
#endif
#define REFINE_BLOCK_BOUNDARY_NET                         1
#define COMPACT_RESOLUTION_CNN                            1
#endif
#if COMPACT_RESOLUTION_CNN
#define DEFAULT_CNN_UP_SAMPLE_MODEL_FILE_CRCNN         "./CaffeModel/Barc_EDSR_Baseline_??_3.caffemodel"
#define DEFAULT_CNN_UP_SAMPLE_NET_FILE_CRCNN           "./CaffeModel/net_edsr_baseline.prototxt"
#define DEFAULT_CNN_DOWN_SAMPLE_NET_FILE_CRCNN         "./CaffeModel/net_crcnn.prototxt"
#define DEFAULT_CNN_DOWN_SAMPLE_MODEL_FILE_CRCNN       "./CaffeModel/test54_params_iter_80000.caffemodel"
#endif
#define LIYE_FIX_BUG                                     1
#define LIYE_DEBUG                                       0
#define SHOW_MODE                                        0
#define LIYE_DEBUG_2                                     0
#define MULTI_COMPRESS_CTU                               0
#define BARC_ONLY_FOR_LUMA                               1  
#endif  //BARC_ON
#define BitDepth10                                        1
#if BitDepth10
#define PixelRangeInt                                    1023
#define PixelRangeDouble                                 1023.0
#else
#define PixelRangeInt                                    255
#define PixelRangeDouble                                 255.0
#endif
#endif //HEM_TOOLS_ON

#define RRC_TOOLS_ON                                      0
#if RRC_TOOLS_ON
#define PVC_AQS                                           0 ///< Adaptive Quantization Step

#if COM16_C1016_AFFINE
#define RRC_WIF                                           0  ///< Affine 2 step bilinear filter
#endif

#define RRC_SELECTABLE_IF                                 0 //enable selectable interpolation filter
#if RRC_SELECTABLE_IF
#define RRC_SIF                                           0
#define RRC_BIF                                           0
#define RRC_MIE                                           1

#define SIF_QTBT_FIX                                      0
#define SIF_CU_LEVEL                                      0
#if !SIF_CU_LEVEL
#define SIF_SPLIT_ONLY                                    1
#endif
#define SIF_DEPTH                                         0
#define SIF_SHIFT_MV_SEARCH                               1  // TODO: check this flag, what is the best

#define SIF_ENC_SPEEDUP                                   0

#define SIF_STD                                           0
#define SIF_BLUR                                          1
#define SIF_SHARP                                         2

#if VCEG_AZ07_IMV
#define SIF_HARMO_IMV                                     0
#endif
#if COM16_C1016_AFFINE
#define SIF_AFFINE_HARMONIZATION                          0  ///< Harmonization with affine
#endif
#if RRC_BIF
#define SHARP_ALL_FRAC                                    0
#define SIF_GAUSS                                         1
#define SIF_GAUSS_ONLY                                    1
#define BIF_BLUR_HALF                                     1
#define BIF_FILTER_ZERO_FRAC                              1
#if SIF_CU_LEVEL
#define BIF_NO_SKIP_ONLY                                  1  ///< don't use BIF for Skip Merge 2Nx2N
#endif
#if RRC_SIF || RRC_MIE
#error wrong combination
#endif
#endif//RRC_BIF
#if RRC_SIF
#define SHARP_ALL_FRAC                                    0
#define SIF_GAUSS                                         1
#define SIF_GAUSS_ONLY                                    0
#define BIF_FILTER_ZERO_FRAC                              0
#define BIF_BLUR_HALF                                     1
#if RRC_BIF || RRC_MIE
#error wrong combination
#endif
#endif
#if RRC_MIE
#define SHARP_ALL_FRAC                                    1
#define SIF_GAUSS                                         0
#define SIF_GAUSS_ONLY                                    0
#if SIF_CU_LEVEL
#define BIF_NO_SKIP_ONLY                                  1  ///< don't use MIE for Skip Merge 2Nx2N
#endif
#if RRC_BIF || RRC_SIF
#error wrong combination
#endif
#endif

#endif //RRC_SELECTABLE_IF

#endif //RRC_TOOLS_ON

///////////////////////////////////////////////////////////
//  Huawei stop section
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// ITH tools section start
///////////////////////////////////////////////////////////
#define ITH_TOOLS                                         1   ///< Enables all changes made by Ittiam

#if ITH_TOOLS

#define ITH_OPT_C1_GEN                                    1   ///< Enable all software optimization related to code rewrite in C / control code optimization / removal of malloc / memset / buffer copies
#define ITH_OPT_C2_ALG                                    1   ///< Enable all algorithmic optimizations including normative and non-normative
#define ITH_OPT_C3_SIMD                                   (1 && COM16_C806_SIMD_OPT)   ///< Enable all SIMD related optimizations
#define ITH_FIX_C4                                        1   ///< Enable bux fixes in the reference
#define ITH_NEWTOOL_C5                                    1   ///< Enable new tools proposed
#define ITH_IMPROV_C6                                     1   ///< Enable optimistations for improvements in existing tools that give quality gains

#endif

#if ITH_OPT_C1_GEN //全部移植

#define ITH_C1_SMT_RESIDUE                                1   ///< Remove redundant buffer copies and replace vector objects with C style arrays
#define ITH_C1_CLIPA                                      1   ///< Optimised ClipA
#define ITH_C1_RDOQ_MEM_OPT                               1   ///< Removed multiple on the fly mallocs and unrequired memsets

#endif

#if ITH_OPT_C2_ALG

#define ITH_C2_CNT_ZR_ZC                                  1   ///< Code for finding zero columns and zero rows of coefficients array [Bit Exact; Time Reduction = ~0]
#define ITH_C2_ZR_ZC_SL                                   (1 && (ITH_C2_CNT_ZR_ZC) && (JVET_C0024_ITSKIP)) ///< Zero detection of columns and rows before dequantization of coefficients for speed up [Bit Exact; Time Reduction = To be tested]

#define ITH_C2_FXD_PT_AFF_PARAMS                          1   ///< Integer Arithmetic instead of floating point arithmetic [Bit Exact]   这个需要先移植，然后，再移植AFFINE汇编（ITH_C3_EQUAL_COEFF_COMPUTER，ITH_C3_HORIZ_SOBEL_FILT，ITH_C3_VERT_SOBEL_FILT）
#define ITH_C2_DIV2MUL_LUT_MOD                            1   ///< Converted Variable Q format to Fixed Q format to make it SIMD friendly keeping precision same as Variable Q format [Bit Exact; Time reduction = 0]


#endif
///////////////////////////////////////////////////////////
// ITH tools section end
///////////////////////////////////////////////////////////
// ====================================================================================================================
// Debugging
// ====================================================================================================================

#define DEBUG_STRING                                      0 ///< When enabled, prints out final decision debug info at encoder and decoder
#define DEBUG_ENCODER_SEARCH_BINS                         0 ///< When enabled, prints out each bin as it is coded during encoder search
#define DEBUG_CABAC_BINS                                  0 ///< When enabled, prints out each bin as it is coded during final encode and decode
#define DEBUG_INTRA_SEARCH_COSTS                          0 ///< When enabled, prints out the cost for each mode during encoder search
#define DEBUG_TRANSFORM_AND_QUANTISE                      0 ///< When enabled, prints out each TU as it passes through the transform-quantise-dequantise-inverseTransform process

#define ENVIRONMENT_VARIABLE_DEBUG_AND_TEST               0 ///< When enabled, allows control of debug modifications via environment variables
#define PRINT_MACRO_VALUES                                1 ///< When enabled, the encoder prints out a list of the non-environment-variable controlled macros and their values on startup

// TODO: rename this macro to DECODER_DEBUG_BIT_STATISTICS (may currently cause merge issues with other branches)
// This can be enabled by the makefile
#ifndef RExt__DECODER_DEBUG_BIT_STATISTICS
#define RExt__DECODER_DEBUG_BIT_STATISTICS                0 ///< 0 (default) = decoder reports as normal, 1 = decoder produces bit usage statistics (will impact decoder run time by up to ~10%)
#endif

// This can be enabled by the makefile
#ifndef ENC_DEC_TRACE
#define ENC_DEC_TRACE                                     0
#endif
#define DEC_NUH_TRACE                                     0 ///< When trace enabled, enable tracing of NAL unit headers at the decoder (currently not possible at the encoder)

#define PRINT_RPS_INFO                                    0 ///< Enable/disable the printing of bits used to send the RPS.

// ====================================================================================================================
// Tool Switches - transitory (these macros are likely to be removed in future revisions)
// ====================================================================================================================

#define DECODER_CHECK_SUBSTREAM_AND_SLICE_TRAILING_BYTES  1 ///< TODO: integrate this macro into a broader conformance checking system.
#define T0196_SELECTIVE_RDOQ                              1 ///< selective RDOQ

#ifndef EXTENSION_360_VIDEO
#define EXTENSION_360_VIDEO                               0   ///< extension for 360/spherical video coding support; this macro should be controlled by makefile, as it would be used to control whether the library is built and linked
#endif

// ====================================================================================================================
// Tool Switches
// ====================================================================================================================

#define ADAPTIVE_QP_SELECTION                             1 ///< G382: Adaptive reconstruction levels, non-normative part for adaptive QP selection

#define JVET_E0059_FLOATING_POINT_QP_FIX                  1 ///< Replace floating point QP with a source-file frame number.
#define JVET_G0101_QP_SWITCHING                           1 ///< After switching POC, increase base QP instead of frame level QP.

#define AMP_ENC_SPEEDUP                                   1 ///< encoder only speed-up by AMP mode skipping
#if AMP_ENC_SPEEDUP
#define AMP_MRG                                           1 ///< encoder only force merge for AMP partition (no motion search for AMP)
#endif

#define FAST_BIT_EST                                      1 ///< G763: Table-based bit estimation for CABAC

#define HHI_RQT_INTRA_SPEEDUP                             1 ///< tests one best mode with full rqt
#define HHI_RQT_INTRA_SPEEDUP_MOD                         0 ///< tests two best modes with full rqt

#if HHI_RQT_INTRA_SPEEDUP_MOD && !HHI_RQT_INTRA_SPEEDUP
#error
#endif

#define MATRIX_MULT                                       0 ///< Brute force matrix multiplication instead of partial butterfly

#define O0043_BEST_EFFORT_DECODING                        0 ///< 0 (default) = disable code related to best effort decoding, 1 = enable code relating to best effort decoding [ decode-side only ].

#define RDOQ_CHROMA_LAMBDA                                1 ///< F386: weighting of chroma for RDOQ

// This can be enabled by the makefile
#ifndef RExt__HIGH_BIT_DEPTH_SUPPORT
#define RExt__HIGH_BIT_DEPTH_SUPPORT                      0 ///< 0 (default) use data type definitions for 8-10 bit video, 1 = use larger data types to allow for up to 16-bit video (originally developed as part of N0188)
#endif


// ====================================================================================================================
// Derived macros
// ====================================================================================================================

#if RExt__HIGH_BIT_DEPTH_SUPPORT
#define FULL_NBIT                                         1 ///< When enabled, use distortion measure derived from all bits of source data, otherwise discard (bitDepth - 8) least-significant bits of distortion
#define RExt__HIGH_PRECISION_FORWARD_TRANSFORM            1 ///< 0 use original 6-bit transform matrices for both forward and inverse transform, 1 (default) = use original matrices for inverse transform and high precision matrices for forward transform
#else
#define FULL_NBIT                                         0 ///< When enabled, use distortion measure derived from all bits of source data, otherwise discard (bitDepth - 8) least-significant bits of distortion
#define RExt__HIGH_PRECISION_FORWARD_TRANSFORM            0 ///< 0 (default) use original 6-bit transform matrices for both forward and inverse transform, 1 = use original matrices for inverse transform and high precision matrices for forward transform
#endif

#if FULL_NBIT
# define DISTORTION_PRECISION_ADJUSTMENT(x)  0
#else
# define DISTORTION_PRECISION_ADJUSTMENT(x) (x)
#endif

#if DEBUG_STRING
  #define DEBUG_STRING_PASS_INTO(name) , name
  #define DEBUG_STRING_PASS_INTO_OPTIONAL(name, exp) , (exp==0)?0:name
  #define DEBUG_STRING_FN_DECLARE(name) , std::string &name
  #define DEBUG_STRING_FN_DECLAREP(name) , std::string *name
  #define DEBUG_STRING_NEW(name) std::string name;
  #define DEBUG_STRING_OUTPUT(os, name) os << name;
  #define DEBUG_STRING_APPEND(str1, str2) str1+=str2;
  #define DEBUG_STRING_SWAP(str1, str2) str1.swap(str2);
  #define DEBUG_STRING_CHANNEL_CONDITION(compID) (true)
  #include <sstream>
  #include <iomanip>
#else
  #define DEBUG_STRING_PASS_INTO(name)
  #define DEBUG_STRING_PASS_INTO_OPTIONAL(name, exp)
  #define DEBUG_STRING_FN_DECLARE(name)
  #define DEBUG_STRING_FN_DECLAREP(name)
  #define DEBUG_STRING_NEW(name)
  #define DEBUG_STRING_OUTPUT(os, name)
  #define DEBUG_STRING_APPEND(str1, str2)
  #define DEBUG_STRING_SWAP(srt1, str2)
  #define DEBUG_STRING_CHANNEL_CONDITION(compID)
#endif

// ====================================================================================================================
// Error checks
// ====================================================================================================================

#if ((RExt__HIGH_PRECISION_FORWARD_TRANSFORM != 0) && (RExt__HIGH_BIT_DEPTH_SUPPORT == 0))
#error ERROR: cannot enable RExt__HIGH_PRECISION_FORWARD_TRANSFORM without RExt__HIGH_BIT_DEPTH_SUPPORT
#endif

// ====================================================================================================================
// Basic type redefinition
// ====================================================================================================================

typedef       void                Void;
typedef       bool                Bool;

#ifdef __arm__
typedef       signed char         Char;
#else
typedef       char                Char;
#endif
#if EXTENSION_360_VIDEO
typedef       char                TChar; // Used for text/characters
#endif
typedef       unsigned char       UChar;
typedef       short               Short;
typedef       unsigned short      UShort;
typedef       int                 Int;
typedef       unsigned int        UInt;
typedef       double              Double;
typedef       float               Float;


// ====================================================================================================================
// 64-bit integer type
// ====================================================================================================================

#ifdef _MSC_VER
typedef       __int64             Int64;

#if _MSC_VER <= 1200 // MS VC6
typedef       __int64             UInt64;   // MS VC6 does not support unsigned __int64 to double conversion
#else
typedef       unsigned __int64    UInt64;
#endif

#else

typedef       long long           Int64;
typedef       unsigned long long  UInt64;

#endif

// ====================================================================================================================
// Named numerical types
// ====================================================================================================================

#if RExt__HIGH_BIT_DEPTH_SUPPORT
typedef       Int             Pel;               ///< pixel type
typedef       Int64           TCoeff;            ///< transform coefficient
typedef       Int             TMatrixCoeff;      ///< transform matrix coefficient
typedef       Short           TFilterCoeff;      ///< filter coefficient
typedef       Int64           Intermediate_Int;  ///< used as intermediate value in calculations
typedef       UInt64          Intermediate_UInt; ///< used as intermediate value in calculations
#else
typedef       Short           Pel;               ///< pixel type
typedef       Int             TCoeff;            ///< transform coefficient
typedef       Short           TMatrixCoeff;      ///< transform matrix coefficient
typedef       Short           TFilterCoeff;      ///< filter coefficient
typedef       Int             Intermediate_Int;  ///< used as intermediate value in calculations
typedef       UInt            Intermediate_UInt; ///< used as intermediate value in calculations
#endif

#if FULL_NBIT
typedef       UInt64          Distortion;        ///< distortion measurement
#else
typedef       UInt            Distortion;        ///< distortion measurement
#endif

#if VCEG_AZ08_KLT_COMMON
#if VCEG_AZ08_USE_SSD_DISTANCE || VCEG_AZ08_USE_SAD_DISTANCE
typedef       int             DistType;
#endif
typedef       Float           covMatrixType;
typedef       float           EigenType;
#endif
// ====================================================================================================================
// Enumeration
// ====================================================================================================================
#if JVET_D0077_SAVE_LOAD_ENC_INFO
enum SaveLoadTag
{
  SAVE_LOAD_INIT = 0,
  SAVE_ENC_INFO = 1,
  LOAD_ENC_INFO = 2
};
#endif

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
enum MergeType
{
#if JVET_C0035_ATMVP_SIMPLIFICATION
  MGR_TYPE_DEFAULT_N  = 0, // 0
  MGR_TYPE_SUBPU_ATMVP = 1, // 1
  MGR_TYPE_SUBPU_ATMVP_EXT =2, // 2
  NUM_MGR_TYPE =3,              // 3
#else
  MGR_TYPE_DEFAULT_N  = 0, // 0
  MGR_TYPE_SUBPU_TMVP = 1, // 1
  MGR_TYPE_SUBPU_TMVP_EXT =2, // 2
  NUM_MGR_TYPE =3, 
#endif
};
#endif

enum RDPCMMode
{
  RDPCM_OFF             = 0,
  RDPCM_HOR             = 1,
  RDPCM_VER             = 2,
  NUMBER_OF_RDPCM_MODES = 3
};

enum RDPCMSignallingMode
{
  RDPCM_SIGNAL_IMPLICIT            = 0,
  RDPCM_SIGNAL_EXPLICIT            = 1,
  NUMBER_OF_RDPCM_SIGNALLING_MODES = 2
};

/// supported slice type
enum SliceType
{
  B_SLICE               = 0,
  P_SLICE               = 1,
  I_SLICE               = 2,
  NUMBER_OF_SLICE_TYPES = 3
};

/// chroma formats (according to semantics of chroma_format_idc)
enum ChromaFormat
{
  CHROMA_400        = 0,
  CHROMA_420        = 1,
  CHROMA_422        = 2,
  CHROMA_444        = 3,
  NUM_CHROMA_FORMAT = 4
};

enum ChannelType
{
  CHANNEL_TYPE_LUMA    = 0,
  CHANNEL_TYPE_CHROMA  = 1,
  MAX_NUM_CHANNEL_TYPE = 2
};

enum ComponentID
{
  COMPONENT_Y       = 0,
  COMPONENT_Cb      = 1,
  COMPONENT_Cr      = 2,
  MAX_NUM_COMPONENT = 3
};

enum InputColourSpaceConversion // defined in terms of conversion prior to input of encoder.
{
  IPCOLOURSPACE_UNCHANGED               = 0,
  IPCOLOURSPACE_YCbCrtoYCrCb            = 1, // Mainly used for debug!
  IPCOLOURSPACE_YCbCrtoYYY              = 2, // Mainly used for debug!
  IPCOLOURSPACE_RGBtoGBR                = 3,
  NUMBER_INPUT_COLOUR_SPACE_CONVERSIONS = 4
};

enum DeblockEdgeDir
{
  EDGE_VER     = 0,
  EDGE_HOR     = 1,
  NUM_EDGE_DIR = 2
};

/// supported partition shape
enum PartSize
{
  SIZE_2Nx2N           = 0,           ///< symmetric motion partition,  2Nx2N
#if JVET_C0024_QTBT
  NUMBER_OF_PART_SIZES = 1
#else
  SIZE_2NxN            = 1,           ///< symmetric motion partition,  2Nx N
  SIZE_Nx2N            = 2,           ///< symmetric motion partition,   Nx2N
  SIZE_NxN             = 3,           ///< symmetric motion partition,   Nx N
  SIZE_2NxnU           = 4,           ///< asymmetric motion partition, 2Nx( N/2) + 2Nx(3N/2)
  SIZE_2NxnD           = 5,           ///< asymmetric motion partition, 2Nx(3N/2) + 2Nx( N/2)
  SIZE_nLx2N           = 6,           ///< asymmetric motion partition, ( N/2)x2N + (3N/2)x2N
  SIZE_nRx2N           = 7,           ///< asymmetric motion partition, (3N/2)x2N + ( N/2)x2N
  NUMBER_OF_PART_SIZES = 8
#endif
};

/// supported prediction type
enum PredMode
{
  MODE_INTER                 = 0,     ///< inter-prediction mode
  MODE_INTRA                 = 1,     ///< intra-prediction mode
  NUMBER_OF_PREDICTION_MODES = 2,
};

/// reference list index
enum RefPicList
{
  REF_PIC_LIST_0               = 0,   ///< reference list 0
  REF_PIC_LIST_1               = 1,   ///< reference list 1
  NUM_REF_PIC_LIST_01          = 2,
  REF_PIC_LIST_X               = 100  ///< special mark
};

/// distortion function index
enum DFunc
{
  DF_DEFAULT         = 0,
  DF_SSE             = 1,      ///< general size SSE
  DF_SSE4            = 2,      ///<   4xM SSE
  DF_SSE8            = 3,      ///<   8xM SSE
  DF_SSE16           = 4,      ///<  16xM SSE
  DF_SSE32           = 5,      ///<  32xM SSE
  DF_SSE64           = 6,      ///<  64xM SSE
  DF_SSE16N          = 7,      ///< 16NxM SSE

  DF_SAD             = 8,      ///< general size SAD
  DF_SAD4            = 9,      ///<   4xM SAD
  DF_SAD8            = 10,     ///<   8xM SAD
  DF_SAD16           = 11,     ///<  16xM SAD
  DF_SAD32           = 12,     ///<  32xM SAD
  DF_SAD64           = 13,     ///<  64xM SAD
  DF_SAD16N          = 14,     ///< 16NxM SAD

  DF_SADS            = 15,     ///< general size SAD with step
  DF_SADS4           = 16,     ///<   4xM SAD with step
  DF_SADS8           = 17,     ///<   8xM SAD with step
  DF_SADS16          = 18,     ///<  16xM SAD with step
  DF_SADS32          = 19,     ///<  32xM SAD with step
  DF_SADS64          = 20,     ///<  64xM SAD with step
  DF_SADS16N         = 21,     ///< 16NxM SAD with step

  DF_HADS            = 22,     ///< general size Hadamard with step
  DF_HADS4           = 23,     ///<   4xM HAD with step
  DF_HADS8           = 24,     ///<   8xM HAD with step
  DF_HADS16          = 25,     ///<  16xM HAD with step
  DF_HADS32          = 26,     ///<  32xM HAD with step
  DF_HADS64          = 27,     ///<  64xM HAD with step
  DF_HADS16N         = 28,     ///< 16NxM HAD with step

  DF_SAD12           = 43,
  DF_SAD24           = 44,
  DF_SAD48           = 45,

  DF_SADS12          = 46,
  DF_SADS24          = 47,
  DF_SADS48          = 48,

  DF_SSE_FRAME       = 50,     ///< Frame-based SSE
#if WCG_LUMA_DQP_CM_SCALE         ///< Weighted SSE
  DF_SSE_WTD             = 51,      ///< general size SSE
  DF_SSE4_WTD            = 52,      ///<   4xM SSE
  DF_SSE8_WTD            = 53,      ///<   8xM SSE
  DF_SSE16_WTD           = 54,      ///<  16xM SSE
  DF_SSE32_WTD           = 55,      ///<  32xM SSE
  DF_SSE64_WTD           = 56,      ///<  64xM SSE
  DF_SSE16N_WTD          = 57,      ///< 16NxM SSE
  DF_DEFAULT_ORI         = 58,
#endif
  DF_TOTAL_FUNCTIONS = 64
};

/// index for SBAC based RD optimization
enum CI_IDX
{
  CI_CURR_BEST = 0,     ///< best mode index
  CI_NEXT_BEST,         ///< next best index
  CI_TEMP_BEST,         ///< temporal index
  CI_CHROMA_INTRA,      ///< chroma intra index
  CI_QT_TRAFO_TEST,
  CI_QT_TRAFO_ROOT,
#if JVET_D0123_ME_CTX_LUT_BITS
  CI_PU_NEXT_BEST,
#endif
  CI_NUM,               ///< total number
};

/// motion vector predictor direction used in AMVP
enum MVP_DIR
{
  MD_LEFT = 0,          ///< MVP of left block
  MD_ABOVE,             ///< MVP of above block
  MD_ABOVE_RIGHT,       ///< MVP of above right block
  MD_BELOW_LEFT,        ///< MVP of below left block
  MD_ABOVE_LEFT         ///< MVP of above left block
};

enum StoredResidualType
{
  RESIDUAL_RECONSTRUCTED          = 0,
  RESIDUAL_ENCODER_SIDE           = 1,
  NUMBER_OF_STORED_RESIDUAL_TYPES = 2
};

enum TransformDirection
{
  TRANSFORM_FORWARD              = 0,
  TRANSFORM_INVERSE              = 1,
  TRANSFORM_NUMBER_OF_DIRECTIONS = 2
};

/// supported ME search methods
enum MESearchMethod
{
  FULL_SEARCH                = 0,     ///< Full search
  DIAMOND                    = 1,     ///< Fast search
  SELECTIVE                  = 2      ///< Selective search
};

#if DIMD_INTRA_PRED
enum TEMPLATE_TYPE
{
  NO_NEIGHBOR = 0,
  LEFT_NEIGHBOR,
  ABOVE_NEIGHBOR,
  LEFT_ABOVE_NEIGHBOR
};
#endif

/// coefficient scanning type used in ACS
enum COEFF_SCAN_TYPE
{
  SCAN_DIAG = 0,        ///< up-right diagonal scan
  SCAN_HOR  = 1,        ///< horizontal first scan
  SCAN_VER  = 2,        ///< vertical first scan
  SCAN_NUMBER_OF_TYPES = 3
};

#if COM16_C806_EMT || COM16_C806_T64
enum TRANS_TYPE
{
  DCT2, DCT5, DCT8, DST1, DST7, NUM_TRANS_TYPE,
  DCT2_HEVC, DCT2_EMT
};
#endif

enum COEFF_SCAN_GROUP_TYPE
{
  SCAN_UNGROUPED   = 0,
  SCAN_GROUPED_4x4 = 1,
  SCAN_NUMBER_OF_GROUP_TYPES = 2
};

#if !VCEG_AZ07_CTX_RESIDUALCODING
enum SignificanceMapContextType
{
  CONTEXT_TYPE_4x4    = 0,
  CONTEXT_TYPE_8x8    = 1,
  CONTEXT_TYPE_NxN    = 2,
  CONTEXT_TYPE_SINGLE = 3,
  CONTEXT_NUMBER_OF_TYPES = 4
};
#endif

enum ScalingListMode
{
  SCALING_LIST_OFF,
  SCALING_LIST_DEFAULT,
  SCALING_LIST_FILE_READ
};

enum ScalingListSize
{
#if JVET_C0024_QTBT
  SCALING_LIST_2x2 = 0,
  SCALING_LIST_4x4,
#else
  SCALING_LIST_4x4 = 0,
#endif
  SCALING_LIST_8x8,
  SCALING_LIST_16x16,
  SCALING_LIST_32x32,
#if COM16_C806_T64
  SCALING_LIST_64x64,
#endif
#if JVET_C0024_QTBT_FIX_QUANT_TICKET25
  SCALING_LIST_128x128,
#endif
  SCALING_LIST_SIZE_NUM
};

// Slice / Slice segment encoding modes
enum SliceConstraint
{
  NO_SLICES              = 0,          ///< don't use slices / slice segments
  FIXED_NUMBER_OF_CTU    = 1,          ///< Limit maximum number of largest coding tree units in a slice / slice segments
  FIXED_NUMBER_OF_BYTES  = 2,          ///< Limit maximum number of bytes in a slice / slice segment
  FIXED_NUMBER_OF_TILES  = 3,          ///< slices / slice segments span an integer number of tiles
};

enum SAOMode //mode
{
  SAO_MODE_OFF = 0,
  SAO_MODE_NEW,
  SAO_MODE_MERGE,
  NUM_SAO_MODES
};

enum SAOModeMergeTypes
{
  SAO_MERGE_LEFT =0,
  SAO_MERGE_ABOVE,
  NUM_SAO_MERGE_TYPES
};


enum SAOModeNewTypes
{
  SAO_TYPE_START_EO =0,
  SAO_TYPE_EO_0 = SAO_TYPE_START_EO,
  SAO_TYPE_EO_90,
  SAO_TYPE_EO_135,
  SAO_TYPE_EO_45,

  SAO_TYPE_START_BO,
  SAO_TYPE_BO = SAO_TYPE_START_BO,

  NUM_SAO_NEW_TYPES
};
#define NUM_SAO_EO_TYPES_LOG2 2

enum SAOEOClasses
{
  SAO_CLASS_EO_FULL_VALLEY = 0,
  SAO_CLASS_EO_HALF_VALLEY = 1,
  SAO_CLASS_EO_PLAIN       = 2,
  SAO_CLASS_EO_HALF_PEAK   = 3,
  SAO_CLASS_EO_FULL_PEAK   = 4,
  NUM_SAO_EO_CLASSES,
};

#define NUM_SAO_BO_CLASSES_LOG2  5
#define NUM_SAO_BO_CLASSES       (1<<NUM_SAO_BO_CLASSES_LOG2)

namespace Profile
{
  enum Name
  {
    NONE = 0,
    MAIN = 1,
    MAIN10 = 2,
    MAINSTILLPICTURE = 3,
    MAINREXT = 4,
    HIGHTHROUGHPUTREXT = 5
  };
}

namespace Level
{
  enum Tier
  {
    MAIN = 0,
    HIGH = 1,
  };

  enum Name
  {
    // code = (level * 30)
    NONE     = 0,
    LEVEL1   = 30,
    LEVEL2   = 60,
    LEVEL2_1 = 63,
    LEVEL3   = 90,
    LEVEL3_1 = 93,
    LEVEL4   = 120,
    LEVEL4_1 = 123,
    LEVEL5   = 150,
    LEVEL5_1 = 153,
    LEVEL5_2 = 156,
    LEVEL6   = 180,
    LEVEL6_1 = 183,
    LEVEL6_2 = 186,
    LEVEL8_5 = 255,
  };
}

enum CostMode
{
  COST_STANDARD_LOSSY              = 0,
  COST_SEQUENCE_LEVEL_LOSSLESS     = 1,
  COST_LOSSLESS_CODING             = 2,
  COST_MIXED_LOSSLESS_LOSSY_CODING = 3
};

enum SPSExtensionFlagIndex
{
  SPS_EXT__REXT           = 0,
//SPS_EXT__MVHEVC         = 1, //for use in future versions
//SPS_EXT__SHVC           = 2, //for use in future versions
  NUM_SPS_EXTENSION_FLAGS = 8
};

enum PPSExtensionFlagIndex
{
  PPS_EXT__REXT           = 0,
//PPS_EXT__MVHEVC         = 1, //for use in future versions
//PPS_EXT__SHVC           = 2, //for use in future versions
  NUM_PPS_EXTENSION_FLAGS = 8
};

// TODO: Existing names used for the different NAL unit types can be altered to better reflect the names in the spec.
//       However, the names in the spec are not yet stable at this point. Once the names are stable, a cleanup
//       effort can be done without use of macros to alter the names used to indicate the different NAL unit types.
enum NalUnitType
{
  NAL_UNIT_CODED_SLICE_TRAIL_N = 0, // 0
  NAL_UNIT_CODED_SLICE_TRAIL_R,     // 1

  NAL_UNIT_CODED_SLICE_TSA_N,       // 2
  NAL_UNIT_CODED_SLICE_TSA_R,       // 3

  NAL_UNIT_CODED_SLICE_STSA_N,      // 4
  NAL_UNIT_CODED_SLICE_STSA_R,      // 5

  NAL_UNIT_CODED_SLICE_RADL_N,      // 6
  NAL_UNIT_CODED_SLICE_RADL_R,      // 7

  NAL_UNIT_CODED_SLICE_RASL_N,      // 8
  NAL_UNIT_CODED_SLICE_RASL_R,      // 9

  NAL_UNIT_RESERVED_VCL_N10,
  NAL_UNIT_RESERVED_VCL_R11,
  NAL_UNIT_RESERVED_VCL_N12,
  NAL_UNIT_RESERVED_VCL_R13,
  NAL_UNIT_RESERVED_VCL_N14,
  NAL_UNIT_RESERVED_VCL_R15,

  NAL_UNIT_CODED_SLICE_BLA_W_LP,    // 16
  NAL_UNIT_CODED_SLICE_BLA_W_RADL,  // 17
  NAL_UNIT_CODED_SLICE_BLA_N_LP,    // 18
  NAL_UNIT_CODED_SLICE_IDR_W_RADL,  // 19
  NAL_UNIT_CODED_SLICE_IDR_N_LP,    // 20
  NAL_UNIT_CODED_SLICE_CRA,         // 21
  NAL_UNIT_RESERVED_IRAP_VCL22,
  NAL_UNIT_RESERVED_IRAP_VCL23,

  NAL_UNIT_RESERVED_VCL24,
  NAL_UNIT_RESERVED_VCL25,
  NAL_UNIT_RESERVED_VCL26,
  NAL_UNIT_RESERVED_VCL27,
  NAL_UNIT_RESERVED_VCL28,
  NAL_UNIT_RESERVED_VCL29,
  NAL_UNIT_RESERVED_VCL30,
  NAL_UNIT_RESERVED_VCL31,

  NAL_UNIT_VPS,                     // 32
  NAL_UNIT_SPS,                     // 33
  NAL_UNIT_PPS,                     // 34
  NAL_UNIT_ACCESS_UNIT_DELIMITER,   // 35
  NAL_UNIT_EOS,                     // 36
  NAL_UNIT_EOB,                     // 37
  NAL_UNIT_FILLER_DATA,             // 38
  NAL_UNIT_PREFIX_SEI,              // 39
  NAL_UNIT_SUFFIX_SEI,              // 40

  NAL_UNIT_RESERVED_NVCL41,
  NAL_UNIT_RESERVED_NVCL42,
  NAL_UNIT_RESERVED_NVCL43,
  NAL_UNIT_RESERVED_NVCL44,
  NAL_UNIT_RESERVED_NVCL45,
  NAL_UNIT_RESERVED_NVCL46,
  NAL_UNIT_RESERVED_NVCL47,
  NAL_UNIT_UNSPECIFIED_48,
  NAL_UNIT_UNSPECIFIED_49,
  NAL_UNIT_UNSPECIFIED_50,
  NAL_UNIT_UNSPECIFIED_51,
  NAL_UNIT_UNSPECIFIED_52,
  NAL_UNIT_UNSPECIFIED_53,
  NAL_UNIT_UNSPECIFIED_54,
  NAL_UNIT_UNSPECIFIED_55,
  NAL_UNIT_UNSPECIFIED_56,
  NAL_UNIT_UNSPECIFIED_57,
  NAL_UNIT_UNSPECIFIED_58,
  NAL_UNIT_UNSPECIFIED_59,
  NAL_UNIT_UNSPECIFIED_60,
  NAL_UNIT_UNSPECIFIED_61,
  NAL_UNIT_UNSPECIFIED_62,
  NAL_UNIT_UNSPECIFIED_63,
  NAL_UNIT_INVALID,
};

// ====================================================================================================================
// Type definition
// ====================================================================================================================

/// parameters for adaptive loop filter
class TComPicSym;

#define MAX_NUM_SAO_CLASSES  32  //(NUM_SAO_EO_GROUPS > NUM_SAO_BO_GROUPS)?NUM_SAO_EO_GROUPS:NUM_SAO_BO_GROUPS

struct SAOOffset
{
  SAOMode modeIdc; // NEW, MERGE, OFF
  Int typeIdc;     // union of SAOModeMergeTypes and SAOModeNewTypes, depending on modeIdc.
  Int typeAuxInfo; // BO: starting band index
  Int offset[MAX_NUM_SAO_CLASSES];

  SAOOffset();
  ~SAOOffset();
  Void reset();

  const SAOOffset& operator= (const SAOOffset& src);
};
#if CNN_BASED_INLOOP_FILTER
struct CNNLFBlkParam
{
	Bool modeIdc[MAX_NUM_COMPONENT];
	Bool& operator[](Int compIdx) { return modeIdc[compIdx]; }
};
#endif
struct SAOBlkParam
{

  SAOBlkParam();
  ~SAOBlkParam();
  Void reset();
  const SAOBlkParam& operator= (const SAOBlkParam& src);
  SAOOffset& operator[](Int compIdx){ return offsetParam[compIdx];}
private:
  SAOOffset offsetParam[MAX_NUM_COMPONENT];

};


struct BitDepths
{
#if O0043_BEST_EFFORT_DECODING
  Int recon[MAX_NUM_CHANNEL_TYPE]; ///< the bit depth used for reconstructing the video
  Int stream[MAX_NUM_CHANNEL_TYPE];///< the bit depth used indicated in the SPS
#else
  Int recon[MAX_NUM_CHANNEL_TYPE]; ///< the bit depth as indicated in the SPS
#endif
};

/// parameters for deblocking filter
typedef struct _LFCUParam
{
  Bool bInternalEdge;                     ///< indicates internal edge
  Bool bLeftEdge;                         ///< indicates left edge
  Bool bTopEdge;                          ///< indicates top edge
} LFCUParam;



//TU settings for entropy encoding
struct TUEntropyCodingParameters
{
  const UInt            *scan;
  const UInt            *scanCG;
        COEFF_SCAN_TYPE  scanType;
        UInt             widthInGroups;
        UInt             heightInGroups;
        UInt             firstSignificanceMapContext;
};


struct TComPictureHash
{
  std::vector<UChar> hash;

  Bool operator==(const TComPictureHash &other) const
  {
    if (other.hash.size() != hash.size())
    {
      return false;
    }
    for(UInt i=0; i<UInt(hash.size()); i++)
    {
      if (other.hash[i] != hash[i])
      {
        return false;
      }
    }
    return true;
  }

  Bool operator!=(const TComPictureHash &other) const
  {
    return !(*this == other);
  }
};

struct TComSEITimeSet
{
  TComSEITimeSet() : clockTimeStampFlag(false),
                     numUnitFieldBasedFlag(false),
                     countingType(0),
                     fullTimeStampFlag(false),
                     discontinuityFlag(false),
                     cntDroppedFlag(false),
                     numberOfFrames(0),
                     secondsValue(0),
                     minutesValue(0),
                     hoursValue(0),
                     secondsFlag(false),
                     minutesFlag(false),
                     hoursFlag(false),
                     timeOffsetLength(0),
                     timeOffsetValue(0)
  { }
  Bool clockTimeStampFlag;
  Bool numUnitFieldBasedFlag;
  Int  countingType;
  Bool fullTimeStampFlag;
  Bool discontinuityFlag;
  Bool cntDroppedFlag;
  Int  numberOfFrames;
  Int  secondsValue;
  Int  minutesValue;
  Int  hoursValue;
  Bool secondsFlag;
  Bool minutesFlag;
  Bool hoursFlag;
  Int  timeOffsetLength;
  Int  timeOffsetValue;
};

struct TComSEIMasteringDisplay
{
  Bool      colourVolumeSEIEnabled;
  UInt      maxLuminance;
  UInt      minLuminance;
  UShort    primaries[3][2];
  UShort    whitePoint[2];
};

#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ07_INIT_PREVFRAME
typedef struct _QPFLAG
{
  UInt      QP;       
  Bool      used;      //same QP, same type has appearaed
  Bool      firstUsed; //same QP, same type was firstly signaled
#if VCEG_AZ07_INIT_PREVFRAME
#if VCEG_AZ07_INIT_PREVFRAME_FIX
  Int       resetInit; 
#else
  UInt      resetInit; //for the first B/P frame after intra slice, no init update
#endif
#endif
} QPFlag;
#endif

#if JVET_E0077_ENHANCED_LM
enum ADDITIONAL_CHROMA_MODE
{
  LM_CHROMA_IDX = 67,
#if JVET_E0077_MMLM
  MMLM_CHROMA_IDX,
#endif
#if JVET_E0077_LM_MF
  LM_CHROMA_F1_IDX,
  LM_CHROMA_F2_IDX,
  LM_CHROMA_F3_IDX,
  LM_CHROMA_F4_IDX,
#endif
};
#endif

#if WCG_LUMA_DQP_CM_SCALE
enum LumaLevelToDQPMode
{
  LUMALVL_TO_DQP_DISABLED = 0,
  LUMALVL_TO_DQP_AVG_METHOD = 1, // use average of CTU to determine luma level
  LUMALVL_TO_DQP_NUM_MODES = 2
};

struct LumaLevelToDeltaQPMapping
{
  LumaLevelToDQPMode                 mode;             ///< use deltaQP determined by block luma level
  Bool                               isSDR;            ///< wheter inputis SDR converted to BT2100 contrainer or true HDR content 
  Double                             maxMethodWeight;  ///< weight of max luma value when mode = 2
  std::vector< std::pair<Int, Int> > mapping;          ///< first=luma level, second=delta QP.
  Bool isEnabled() const { return mode != LUMALVL_TO_DQP_DISABLED; }
};

struct WCGChromaQPControl
{
  Bool isEnabled() const { return enabled; }
  Bool   enabled;         ///< Enabled flag (0:default)
  Double chromaCbQpScale; ///< Chroma Cb QP Scale (1.0:default)
  Double chromaCrQpScale; ///< Chroma Cr QP Scale (1.0:default)
  Double chromaQpScale;   ///< Chroma QP Scale (0.0:default)
  Double chromaQpOffset;  ///< Chroma QP Offset (0.0:default)
};
#endif
//! \}

#endif


