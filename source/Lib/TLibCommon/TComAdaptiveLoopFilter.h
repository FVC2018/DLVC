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

/** \file     TComAdaptiveLoopFilter.h
    \brief    adaptive loop filter class (header)
*/

#ifndef __TCOMADAPTIVELOOPFILTER__
#define __TCOMADAPTIVELOOPFILTER__

#include "CommonDef.h"
#include "TComPic.h"

#if ALF_HM3_REFACTOR
#if JVET_C0038_GALF
const static Int ALF_NO_VAR_BIN = 25;
#else
const static Int ALF_NO_VAR_BIN = 16;
#endif


typedef struct _AlfParam
{
  Int alf_flag;                           ///< indicates use of ALF
  Int cu_control_flag;                    ///< coding unit based control flag
  Int chroma_idc;                         ///< indicates use of ALF for chroma
  Int tap;                                ///< number of filter taps - horizontal
  Int tapV;                               ///< number of filter taps - vertical
  Int num_coeff;                          ///< number of filter coefficients
#if !JVET_C0038_GALF
  Int *coeff;                             ///< filter coefficient array
#endif
  Int tap_chroma;                         ///< number of filter taps (chroma)
  Int num_coeff_chroma;                   ///< number of filter coefficients (chroma)
  Int *coeff_chroma;                      ///< filter coefficient array (chroma)

  //CodeAux related
  Int realfiltNo;
  Int filtNo;
  Int filterPattern[ALF_NO_VAR_BIN];
#if JVET_C0038_GALF
  Char  PrevFiltIdx  [ALF_NO_VAR_BIN];
  UChar iAvailableFilters;
  UChar iPredPattern; //whether predicted from fixed filters
#endif
  Int startSecondFilter;
  Int noFilters;
  Int varIndTab[ALF_NO_VAR_BIN];
#if COM16_C806_ALF_TEMPPRED_NUM
  Bool temproalPredFlag; //indicate whether reuse previous ALF coefficients
  Int  prevIdx;          //index of the reused ALF coefficients
  Int  **alfCoeffLuma;    
  Int  *alfCoeffChroma;
#endif
  //Coeff send related
  Int filters_per_group_diff; //this can be updated using codedVarBins
  Int filters_per_group;
#if JVET_C0038_GALF
  Bool forceCoeff0;
  Bool codedVarBins[ALF_NO_VAR_BIN]; 
#else
  Int codedVarBins[ALF_NO_VAR_BIN]; 
  Int forceCoeff0;
#endif
  Int predMethod;
  Int **coeffmulti; //to be written in bitstream, after prediction from fixed filter or previous calss
  Int minKStart;
  Int maxScanVal;
  Int kMinTab[42];
  UInt num_alf_cu_flag;
  UInt num_cus_in_frame;
  UInt alf_max_depth;
  UInt *alf_cu_flag;

}ALFParam;

typedef unsigned short imgpel;

void destroyMatrix_int(int **m2D);
void initMatrix_int(int ***m2D, int d1, int d2);

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// adaptive loop filter class
class TComAdaptiveLoopFilter
{
protected:
#if !JVET_C0038_GALF
  static const Int m_ALF_VAR_SIZE_H        = 4;
  static const Int m_ALF_VAR_SIZE_W        = 4;
#endif
  static const Int m_ALF_WIN_VERSIZE       = 32;
  static const Int m_ALF_WIN_HORSIZE       = 32;

  static const Int m_ALF_MAX_NUM_TAP       = 9;                                     ///< maximum number of filter taps (9x9)
  static const Int m_ALF_MIN_NUM_TAP       = 5;                                     ///< minimum number of filter taps (5x5)
  static const Int m_ALF_MAX_NUM_TAP_C     = 5;                                     ///< number of filter taps for chroma (5x5)
  static const Int m_ALF_MAX_NUM_COEF      = 42;                                    ///< maximum number of filter coefficients
  static const Int m_ALF_MIN_NUM_COEF      = 14;                                    ///< minimum number of filter coefficients
  static const Int m_ALF_MAX_NUM_COEF_C    = 14;                                    ///< number of filter taps for chroma
  static const Int m_ALF_NUM_BIT_SHIFT     = 8;                                     ///< bit shift parameter for quantization of ALF param.
  static const Int m_ALF_ROUND_OFFSET      = ( 1 << ( m_ALF_NUM_BIT_SHIFT - 1 ) );  ///< rounding offset for ALF quantization
  
  static const Int m_VAR_SIZE              = 1;                                     ///< JCTVC-E323+E046

  static const Int m_FILTER_LENGTH         = 9;
#if JVET_C0038_GALF
  static const Int m_NUM_ALF_SETS_PERDIR   = 5;
  static const Int m_PADDING_W_ALF         = ( m_FILTER_LENGTH/2 > JVET_C0038_SHIFT_VAL_HALFW ? m_FILTER_LENGTH/2 : JVET_C0038_SHIFT_VAL_HALFW);
#endif
  static const Int m_ALF_HM3_QC_CLIP_RANGE = 1024;
  static const Int m_ALF_HM3_QC_CLIP_OFFSET= 384;

public:
#if JVET_C0038_GALF 
  static const Int m_NUM_BITS              =  10;
#else
  static const Int m_NUM_BITS              =  9;
#endif
  static const Int m_NO_TEST_FILT          =  3;                                    ///< Filter supports (5/7/9)
#if JVET_C0038_GALF
  #define DIR_TH                           2.0
  #define NO_VALS_LAGR                     5
  #define NO_VALS_LAGR_SHIFT               3
  static const Int m_NO_VAR_BINS           = 25;
  static const Int m_NO_FILTERS            = 25; 
  static const Int m_MAX_SQR_FILT_LENGTH   = ((m_FILTER_LENGTH*m_FILTER_LENGTH) / 2 + 1);
  static const Int m_MAX_SQT_FILT_SYM_LENGTH = ((m_FILTER_LENGTH*m_FILTER_LENGTH) / 4 + 1);
  static const Int m_SQR_FILT_LENGTH_9SYM  = ((9*9) / 4 + 1); 
  static const Int m_SQR_FILT_LENGTH_7SYM  = ((7*7) / 4 + 1); 
  static const Int m_SQR_FILT_LENGTH_5SYM  = ((5*5) / 4 + 1); 
#else
  static const Int m_NO_VAR_BINS           = 16; 
  static const Int m_NO_FILTERS            = 16; 
  static const Int m_MAX_SQR_FILT_LENGTH   = ((m_FILTER_LENGTH*m_FILTER_LENGTH) / 2 + 2);
  static const Int m_SQR_FILT_LENGTH_9SYM  = ((9*9) / 4 + 2 - 1); 
  static const Int m_SQR_FILT_LENGTH_7SYM  = ((7*7) / 4 + 2); 
  static const Int m_SQR_FILT_LENGTH_5SYM  = ((5*5) / 4 + 2); 
#endif
  static const Int m_MAX_SCAN_VAL          = 11;
  static const Int m_MAX_EXP_GOLOMB        = 16;

  // quantized filter coefficients
  static const Int m_aiSymmetricMag9x9[41];                                         ///< quantization scaling factor for 9x9 filter
  static const Int m_aiSymmetricMag7x7[25];                                         ///< quantization scaling factor for 7x7 filter
  static const Int m_aiSymmetricMag5x5[13];                                         ///< quantization scaling factor for 5x5 filter
  static const Int m_aiSymmetricMag9x7[32];                                         ///< quantization scaling factor for 9x7 filter
  
  // temporary picture buffer
  TComPicYuv*   m_pcTempPicYuv;                                                     ///< temporary picture buffer for ALF processing

public:
  static const Int* m_pDepthIntTab[m_NO_TEST_FILT];
#if JVET_C0038_GALF
  Int m_max_NO_VAR_BINS; 
  Int m_max_NO_FILTERS;
#endif

protected:
#if JVET_C0038_GALF
  static const Int m_ALFfilterCoeffFixed[m_NO_FILTERS*JVET_C0038_NO_PREV_FILTERS][21]; /// fixed filters used in ALF.
  static const Int depthInt9x9Cut[21];
  static const Int depthInt7x7Cut[14];
  static const Int depthInt5x5Cut[8];
#else
  static const Int  m_depthInt9x9Sym[21];
  static const Int  m_depthInt7x7Sym[14];
  static const Int  m_depthInt5x5Sym[8];
#endif
  // ------------------------------------------------------------------------------------------------------------------
  // For luma component
  // ------------------------------------------------------------------------------------------------------------------
#if JVET_C0038_GALF 
  static Int m_pattern9x9Sym[41];
  static Int m_weights9x9Sym[22];
#else
  static Int m_pattern9x9Sym[39];
  static Int m_weights9x9Sym[21];
#endif
  static Int m_pattern9x9Sym_Quart[42];
  static Int m_pattern7x7Sym[25];
  static Int m_weights7x7Sym[14];
  static Int m_pattern7x7Sym_Quart[42];
  static Int m_pattern5x5Sym[13];
  static Int m_weights5x5Sym[8];
  static Int m_pattern5x5Sym_Quart[45];
  static Int m_pattern9x9Sym_9[39];
  static Int m_pattern9x9Sym_7[25];
  static Int m_pattern9x9Sym_5[13];
  
  static Int *m_patternTab_filt[m_NO_TEST_FILT];
  static Int m_flTab[m_NO_TEST_FILT];
  static Int *m_patternTab[m_NO_TEST_FILT]; 
  static Int *m_patternMapTab[m_NO_TEST_FILT];
  static Int *m_weightsTab[m_NO_TEST_FILT];
  static Int m_sqrFiltLengthTab[m_NO_TEST_FILT];
  
  Int       m_img_height,m_img_width;
  Int       m_nInputBitDepth;
  Int       m_nInternalBitDepth;
  Int       m_nBitIncrement;
  Int       m_nIBDIMax;
  UInt      m_uiMaxTotalCUDepth;
  
  imgpel**  m_imgY_var;
  Int**     m_imgY_temp;
  
  Int**     m_imgY_ver;
  Int**     m_imgY_hor;
#if JVET_C0038_GALF
  Int**     m_imgY_dig0 ;
  Int**     m_imgY_dig1 ;
#endif
  imgpel**  m_varImgMethods;
  Int **    m_filterCoeffSym;
#if JVET_C0038_GALF
  Int **    m_filterCoeffFinal;
#endif
  Int **    m_filterCoeffPrevSelected;
  Short **  m_filterCoeffShort;
  imgpel *  m_alfClipTable;
  Int       m_alfClipOffset;
  Int **    m_filterCoeffTmp;
  Int **    m_filterCoeffSymTmp;
  UInt      m_uiNumCUsInFrame;
  /// ALF for luma component
  Void xALFLuma_qc( TComPic* pcPic, ALFParam* pcAlfParam, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest);
  
  Void reconstructFilterCoeffs(ALFParam* pcAlfParam,int **pfilterCoeffSym, int bit_depth);
  Void getCurrentFilter(int **filterCoeffSym,ALFParam* pcAlfParam);
  // memory allocation
  Void destroyMatrix_short(short **m2D);
  Void initMatrix_short(short ***m2D, int d1, int d2);
  Void destroyMatrix_imgpel(imgpel **m2D);
  Void destroyMatrix_int(int **m2D);
  Void initMatrix_int(int ***m2D, int d1, int d2);
  Void initMatrix_imgpel(imgpel ***m2D, int d1, int d2);
  Void destroyMatrix4D_double(double ****m4D, int d1, int d2);
  Void destroyMatrix3D_double(double ***m3D, int d1);
  Void destroyMatrix_double(double **m2D);
  Void initMatrix4D_double(double *****m4D, int d1, int d2, int d3, int d4);
  Void initMatrix3D_double(double ****m3D, int d1, int d2, int d3);
  Void initMatrix_double(double ***m2D, int d1, int d2);
  Void free_mem2Dpel(imgpel **array2D);
  Void get_mem2Dpel(imgpel ***array2D, int rows, int columns);
  Void no_mem_exit(const char *where);
  Void xError(const char *text, int code);
  Void calcVar(imgpel **imgY_var, imgpel *imgY_pad, int pad_size, int fl, int img_height, int img_width, int img_stride, int start_width = 0 , int start_height = 0 );
  Void xCalcVar(imgpel **imgY_var, imgpel *imgY_pad, int pad_size, int fl, int img_height, int img_width, int img_stride, int start_width , int start_height );
#if JVET_C0038_GALF
  Void xCalcVarPerPixel(imgpel **imgY_var, imgpel *imgY_pad, Int pad_size, Int fl, Int img_height, Int img_width, Int img_stride, Int start_width , Int start_height);
#endif

  Void DecFilter_qc(imgpel* imgY_rec,ALFParam* pcAlfParam, int Stride);
#if JVET_C0024_QTBT
  Void xSubCUAdaptive_qc(TComDataCU* pcCU, ALFParam* pcAlfParam, imgpel *imgY_rec_post, imgpel *imgY_rec, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, Int Stride);
#else
  Void xSubCUAdaptive_qc(TComDataCU* pcCU, ALFParam* pcAlfParam, imgpel *imgY_rec_post, imgpel *imgY_rec, UInt uiAbsPartIdx, UInt uiDepth, Int Stride);
#endif
  Void xCUAdaptive_qc(TComPic* pcPic, ALFParam* pcAlfParam, imgpel *imgY_rec_post, imgpel *imgY_rec, Int Stride);
#if JVET_C0038_GALF
  Void subfilterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, ALFParam* pcAlfPara, int start_height, int end_height, int start_width, int end_width, int Stride,
                    #if JVET_D0033_ADAPTIVE_CLIPPING
                    ComponentID compid
                    #else
                    Bool bChroma = false
                    #endif
                    );
  Void filterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, ALFParam* pcAlfPara, int Stride);
#else
  Void subfilterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, int filtNo, int start_height, int end_height, int start_width, int end_width, int Stride);
  Void filterFrame(imgpel *imgY_rec_post, imgpel *imgY_rec, int filtNo, int Stride);
#endif
#if JVET_C0024_QTBT
  Void setAlfCtrlFlags (ALFParam *pAlfParam, TComDataCU *pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight, UInt &idx);
#else
  Void setAlfCtrlFlags (ALFParam *pAlfParam, TComDataCU *pcCU, UInt uiAbsPartIdx, UInt uiDepth, UInt &idx);
#endif
#if JVET_C0038_GALF
  Int selectTransposeVarInd(Int varInd, Int *transpose);
#endif  
  // ------------------------------------------------------------------------------------------------------------------
  // For chroma component
  // ------------------------------------------------------------------------------------------------------------------
  
  /// ALF for chroma component
  Void xALFChroma   ( ALFParam* pcAlfParam, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest );
  
  /// sub function: non-adaptive ALF process for chroma
#if JVET_C0038_GALF
  Void xFrameChroma ( ALFParam* pcAlfParam, TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Int iColor );
#else
  Void xFrameChroma ( TComPicYuv* pcPicDec, TComPicYuv* pcPicRest, Int *qh, Int iTap, Int iColor );
#endif

public:
  TComAdaptiveLoopFilter();
  virtual ~TComAdaptiveLoopFilter() {}
  
  // initialize & destory temporary buffer
  Void create  ( Int iPicWidth, Int iPicHeight, ChromaFormat chromaFormatIDC, Int uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxCUDepth , Int nInputBitDepth , Int nInternalBitDepth );
  Void destroy ();
  
  // alloc & free & set functions
  Void allocALFParam  ( ALFParam* pAlfParam );
  Void freeALFParam   ( ALFParam* pAlfParam );
  Void copyALFParam   ( ALFParam* pDesAlfParam, ALFParam* pSrcAlfParam );
  Void setNumCUsInFrame        (TComPic *pcPic);
#if COM16_C806_ALF_TEMPPRED_NUM
  Void setNumCUsInFrame        (UInt uiNumCUsInFrame);
#endif
#if JVET_C0038_GALF  
  Void resetALFPredParam(ALFParam *pAlfParam, Bool bIntra);
#else
  Void printALFParam  ( ALFParam* pAlfParam , Bool bDecoder = false );
#endif
  Void resetALFParam  ( ALFParam* pDesAlfParam);

  // predict filter coefficients
#if JVET_C0038_GALF 
  Void initVarForChroma       ( ALFParam* pcAlfParam, Bool bUpdatedDCCoef = true);
#else
  Void predictALFCoeff        ( ALFParam* pAlfParam );                  ///< prediction of luma ALF coefficients
  Void predictALFCoeffChroma  ( ALFParam* pAlfParam );                  ///< prediction of chroma ALF coefficients
#endif
  // interface function
  Void ALFProcess             ( TComPic* pcPic, ALFParam* pcAlfParam); ///< interface function for ALF process

#if FIX_TICKET12
  Bool refreshAlfTempPred( NalUnitType nalu , Int poc );
#endif
  
  static Int ALFTapHToTapV(Int tapH);
  static Int ALFTapHToNumCoeff(Int tapH);
  static Int ALFFlHToFlV(Int flH);
};
#endif
#endif
