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

/** \file     TComPrediction.h
    \brief    prediction class (header)
*/

#ifndef __TCOMPREDICTION__
#define __TCOMPREDICTION__


// Include files
#include "TComYuv.h"
#include "TComInterpolationFilter.h"
#include "TComWeightPrediction.h"
#if VCEG_AZ07_FRUC_MERGE
#include "TComRdCost.h"
#include <list>
#endif
#if COM16_C1046_PDPC_INTRA
#include "TComRom.h"
#endif
#include "TComPic.h"


// forward declaration
class TComMv;
class TComTU; 
#if VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
class TComMvField;
#endif

#if RRC_DWDIP 
typedef struct
{
  Double gradualWeight; //[0..1]
} RSGIntraParams;
#endif

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// prediction class
typedef enum PRED_BUF_E
{
  PRED_BUF_UNFILTERED=0,
  PRED_BUF_FILTERED=1,
  NUM_PRED_BUF=2
} PRED_BUF;

#if JVET_C0024_QTBT
static const UInt MAX_INTRA_FILTER_DEPTHS=MAX_CU_DEPTH;
#else
#if COM16_C806_LARGE_CTU
static const UInt MAX_INTRA_FILTER_DEPTHS=MAX_CU_DEPTH-1;
#else
static const UInt MAX_INTRA_FILTER_DEPTHS=5;
#endif
#endif

class TComPrediction : public TComWeightPrediction
{
#if TRACE_YUV
  friend class TDecCu;
#endif
private:
  static const UChar m_aucIntraFilter[MAX_NUM_CHANNEL_TYPE][MAX_INTRA_FILTER_DEPTHS];
#if DIMD_NUM_INTRA_DIR_INC
  static const UChar m_aucExtIntraFilter[MAX_NUM_CHANNEL_TYPE][MAX_INTRA_FILTER_DEPTHS];
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UInt*  m_puiW;
  UInt*  m_puiH;
  UInt*  m_puiSPAddr;
#endif

protected:
#if VCEG_AZ05_BIO  
  Pel*   m_pGradX0;
  Pel*   m_pGradY0;
  Pel*   m_pGradX1;
  Pel*   m_pGradY1;
#if !JVET_F0028_BIO_NO_BLOCK_EXTENTION
  Pel*   m_pPred0 ;
  Pel*   m_pPred1 ;
#endif
  Int    iRefListIdx;
#endif

#if COM16_C1046_PDPC_INTRA
  Int* piTempRef;
  Int* piFiltRef;
  Int* piBinBuff;
#endif
  Pel*      m_piYuvExt[MAX_NUM_COMPONENT][NUM_PRED_BUF];
#if RRC_INTRA_SHARPENING
  // TODO: optimize mem
  Pel*      m_piYuvSharpExt;
  Pel*      m_piYuvSharpExt2;
  Pel*      m_piYuvSharpExt3;
#endif
  Int       m_iYuvExtSize;
#if MLIP_EXTEND_REF_LINE
  Pel*      m_piSecondYuvExt[MAX_NUM_COMPONENT][NUM_PRED_BUF];
  Int       m_iSecondYuvExtSize;
#endif
#if BARC
  Pel*      m_piYuvExtDscale[MAX_NUM_COMPONENT][NUM_PRED_BUF];
#endif
  TComYuv   m_acYuvPred[NUM_REF_PIC_LIST_01];
  TComYuv   m_cYuvPredTemp;
  TComYuv m_filteredBlock[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS];
  TComYuv m_filteredBlockTmp[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS];

  TComInterpolationFilter m_if;

  Pel*   m_pLumaRecBuffer;       ///< array for downsampled reconstructed luma sample
  Int    m_iLumaRecStride;       ///< stride of #m_pLumaRecBuffer array

#if JVET_E0077_LM_MF
  Pel*   m_pLumaRecBufferMul[LM_FILTER_NUM];
#endif

#if JVET_E0077_ENHANCED_LM
  Int  m_iCurAngMode;
#endif

#if COM16_C806_LMCHROMA
  UInt m_uiaLMShift[ 32 ];       // Table for multiplication to substitue of division operation
#endif

#if VCEG_AZ06_IC
  UInt   m_uiaICShift[ 64 ];     // Table for multiplication to substitue of division operation
  static const Int m_ICRegCostShift = 7;
  static const Int m_ICConstShift = 5;
  static const Int m_ICShiftDiff = 12;

#if JVET_G0082
	UInt   m_uiaBIOShift[64];
#endif

#if RRC_WIF
  TComYuv         m_cYuvScaleTemp[2];
#endif
#endif

#if VCEG_AZ08_INTER_KLT
  TComPicYuv* m_tempPicYuv;
#endif
  
#if RRC_DWDIP
  inline Int getIntPart(Int pos);

  inline Int getFractPart(Int pos);
  Void xInterpolateRefSamplesAng(Int               iWidth,
    Int               iHeight,
    UInt              dirMode,
    Pel const * const pRefAbove,
    Pel const * const pRefLeft,
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
    Int*              pRefGen
#else
    Double*           pRefGen
#endif
#if JVET_D0033_ADAPTIVE_CLIPPING
    , ComponentID       compID
#endif
#if VCEG_AZ07_INTRA_4TAP_FILTER
    , Bool              bEnable4TapFilter,
    Bool              bIsBicubicInterpUsed = true
#endif
    );
  Int xInterpolateIntraAng(Int               xPos,
    Int               yPos,
    Int               iWidth,
    Int               iHeight,
    Int               intraPredAngle,
    Pel const * const pRef
#if JVET_D0033_ADAPTIVE_CLIPPING
    , ComponentID       compID
#endif
#if VCEG_AZ07_INTRA_4TAP_FILTER
    , Bool              bEnable4TapFilter,
    Bool              bIsBicubicInterpUsed = true
#endif
    );

  Void xPredictRowLUT(Pel * pDst  // pointer to picture start
    , Pel * pDsty // pointer to the current row
    , Int dstStride
    , Int y
    , pair<Int, Int> x_range
    , Pel * refMain
    , Pel * refSide
    , Int deltaPos
    , Int bitDepth
    , Int width
    , Int height
    , const Int * piWeight
    , const Int * pbDecision
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
    , const Int * piSideDeltaPos
#endif 
#if JVET_D0033_ADAPTIVE_CLIPPING
    , ComponentID compID
#else
    , ChannelType channelType
#endif
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
    , Int*    pGenRefSamples
#else
    , Double* pGenRefSamples
#endif
  );

#endif//RRC_DWDIP

    Void xPredIntraAng            ( Int bitDepth, const Pel* pSrc, Int srcStride, Pel* pDst, Int dstStride, UInt width, UInt height,
#if JVET_D0033_ADAPTIVE_CLIPPING
                                    ComponentID compID,
#else
                                    ChannelType channelType,
#endif
                                    UInt dirMode, const Bool bEnableEdgeFilters
#if MLIP_EXTEND_REF_LINE
        , UInt ExternRef
#endif
#if RRC_INTRA_PRED_DATACU_ABSPARTIDX
      , TComDataCU* const pcCU, const UInt uiAbsPartIdx
#endif
#if VCEG_AZ07_INTRA_4TAP_FILTER
    , Bool enable4TapFilter = false
#endif
#if DIMD_NUM_INTRA_DIR_INC
      , Bool bExtIntraDir = false
#endif
#if INTRA_SHARP_4TAP_HARM
    , Bool enableSharp = false
#endif
#if RRC_DWDIP
    , Int iAngModeIdx = 0
#endif
    );


#if RRC_DWDIP_MLR
#if JVET_D0033_ADAPTIVE_CLIPPING && RRC_DWDIP_MLR_IMPR
public:
  static inline Double getMaxTangent(const Int iMinVal, const Int iMaxVal, const Int iBlockSize)
  {
    return(static_cast<Double>(iMaxVal - iMinVal) / static_cast<Double>(iBlockSize));
  }

protected:
  template<typename T>
  static inline T ClipCoeff(const ComponentID compID, const UInt uiBlkSideLength, const T coeff)
  {
    assert(uiBlkSideLength >= 2);
    UInt uiBlkSideLengthIdx;
    T retVal;

    switch (compID)
    {
    case COMPONENT_Y:  uiBlkSideLengthIdx = g_aucConvertToBit[uiBlkSideLength];
      retVal = Clip3(static_cast<T>(-g_ClipParam.Y().dMaxAngle[uiBlkSideLengthIdx]),
        static_cast<T>(g_ClipParam.Y().dMaxAngle[uiBlkSideLengthIdx]), coeff);
      break;
    case COMPONENT_Cb: uiBlkSideLengthIdx = (2 == uiBlkSideLength) ? 0 : (g_aucConvertToBit[uiBlkSideLength] + 1);
      retVal = Clip3(static_cast<T>(-g_ClipParam.U().dMaxAngle[uiBlkSideLengthIdx]),
        static_cast<T>(g_ClipParam.U().dMaxAngle[uiBlkSideLengthIdx]), coeff);
      break;
    case COMPONENT_Cr: uiBlkSideLengthIdx = (2 == uiBlkSideLength) ? 0 : (g_aucConvertToBit[uiBlkSideLength] + 1);
      retVal = Clip3(static_cast<T>(-g_ClipParam.V().dMaxAngle[uiBlkSideLengthIdx]),
        static_cast<T>(g_ClipParam.V().dMaxAngle[uiBlkSideLengthIdx]), coeff);
      break;
    default: std::cerr << "ClipA: Invalid compID value " << compID << " . Exiting." << std::endl; assert(false); exit(0);
    }

    return(retVal);
  }
#endif
  Void xPredIntraGlobalPlane( const  Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height
                            , UInt refWidth, UInt refHeight
                            , const Int bitDepth
#if RRC_DWDIP_MLR_IMPR
                            , const ComponentID compID
#endif
#if RRC_INTRA_PRED_DATACU_ABSPARTIDX
                            , TComDataCU* const pcCU, const UInt uiAbsPartIdx
#endif
                            , Bool bEnableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
                            , Bool enable4TapFilter
#endif
    );
    //alias function
    Void xPredIntraGlobalPlane(const Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height, const Int bitDepth
#if RRC_DWDIP_MLR_IMPR
                             , const ComponentID compID
#endif
#if RRC_INTRA_PRED_DATACU_ABSPARTIDX
                             , TComDataCU* const pcCU, const UInt uiAbsPartIdx
#endif
                             , Bool bEnableEdgeFilters = false 
#if VCEG_AZ07_INTRA_4TAP_FILTER
                             , Bool enable4TapFilter = false
#endif   
                              )
    {
      xPredIntraGlobalPlane(pSrc, srcStride, rpDst, dstStride, width, height, width+1, height+1, bitDepth
#if RRC_DWDIP_MLR_IMPR
                          , compID
#endif
#if RRC_INTRA_PRED_DATACU_ABSPARTIDX
                         , pcCU, uiAbsPartIdx
#endif
                          , bEnableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
                          , enable4TapFilter
#endif   
                           );
    } 
#endif // RRC_DWDIP_MLR

  Void xPredIntraPlanar         ( const Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height );

#if VCEG_AZ07_FRUC_MERGE
  TComRdCost              m_cFRUCRDCost;
  std::list <TComMvField> m_listMVFieldCand[2];
#if JVET_E0060_FRUC_CAND
  RefPicList              m_bilatBestRefPicList;
#endif
  TComYuv                 m_cYuvPredFrucTemplate[2];      // 0: top, 1: left
  Bool                    m_bFrucTemplateAvailabe[2];
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UChar                   m_eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComMvField *           m_cMvFieldSP[NUM_MGR_TYPE];
  UChar *                 m_uhInterDirSP[NUM_MGR_TYPE];
#else
  TComMvField *           m_cMvFieldSP[2];
  UChar *                 m_uhInterDirSP[2];
#endif
#endif
#endif

  // motion compensation functions
#if VCEG_AZ05_BIO  
#define BIO_FILTER_LENGTH                 6
#define BIO_FILTER_LENGTH_MINUS_1         (BIO_FILTER_LENGTH-1)
#define BIO_FILTER_HALF_LENGTH_MINUS_1    ((BIO_FILTER_LENGTH>>1)-1)
#if !JVET_F0028_BIO_NO_BLOCK_EXTENTION
  Void  xPredInterFrac(Pel* ref,Pel* dst,Int dstStride,Int refStride,Int xFrac,Int yFrac,Int width, Int height,Bool bi,ChromaFormat chFmt, const Int bitDepth);
#endif
  Void  xGradFilterX(Pel*  piRefY, Int iRefStride,Pel*  piDstY,Int iDstStride, Int iWidth, Int iHeight,Int iMVyFrac,Int iMVxFrac, const Int bitDepth);
  Void  xGradFilterY(Pel*  piRefY, Int iRefStride,Pel*  piDstY,Int iDstStride, Int iWidth, Int iHeight,Int iMVyFrac,Int iMVxFrac, const Int bitDepth);
  __inline Void gradFilter2DVer (Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMv, const Int iShift);
  __inline Void gradFilter2DHor (Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMV, const Int iShift);
  __inline Void fracFilter2DHor(Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMV, const Int iShift);
  __inline Void fracFilter2DVer(Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMv, const Int iShift);
  __inline Void gradFilter1DHor (Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMV, const Int iShift);
  __inline Void gradFilter1DVer (Pel* piSrc, Int iSrcStride,  Int iWidth, Int iHeight, Int iDstStride,  Pel*& rpiDst, Int iMV, const Int iShift);

#if JVET_G0082
	__inline Int64 divide64(Int64 numer, Int64 denom);
	__inline Void calcBlkGradient(Int sx, Int sy, Int64 *arraysGx2, Int64 *arraysGxGy, Int64 *arraysGxdI, Int64 *arraysGy2, Int64 *arraysGydI,
			Int64 &sGx2, Int64 &sGy2, Int64 &sGxGy, Int64 &sGxdI, Int64 &sGydI, Int iWidth, Int iHeight);
#endif
#endif
  Void xPredInterUni            ( TComDataCU* pcCU,                          UInt uiPartAddr,               Int iWidth, Int iHeight, RefPicList eRefPicList, TComYuv* pcYuvPred
#if JVET_E0052_DMVR
    , Bool bRefineflag
#endif
#if VCEG_AZ05_BIO
    , Bool bBIOApplied =false
#endif
    , Bool bi=false          
#if VCEG_AZ07_FRUC_MERGE || RRC_SELECTABLE_IF
    , Bool bOBMC = false
#endif
    );
  Void xPredInterBi             ( TComDataCU* pcCU,                          UInt uiPartAddr,               Int iWidth, Int iHeight,                         TComYuv* pcYuvPred          
#if JVET_E0052_DMVR    
    , Bool bRefineflag
#endif
#if VCEG_AZ07_FRUC_MERGE || JVET_G0082 || RRC_SELECTABLE_IF
    , Bool bOBMC = false
#endif
    );
#if RRC_SELECTABLE_IF
public:
#endif
  Void xPredInterBlk(const ComponentID compID, TComDataCU *cu, TComPicYuv *refPic, UInt partAddr, TComMv *mv, Int width, Int height, TComYuv *dstPic, Bool bi, const Int bitDepth
#if VCEG_AZ05_BIO                  
    , Bool bBIOapplied =false
#endif
#if VCEG_AZ07_FRUC_MERGE
    , Int nFRUCMode = FRUC_MERGE_OFF
#endif
#if VCEG_AZ06_IC
    , Bool bICFlag      = false
#endif
#if ITH_C2_MC_OPT
    , Bool doPred       = true
#endif
    );
#if JVET_E0052_DMVR
    Void xBIPMVRefine(TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefPicList, Int iWidth, Int iHeight, TComYuv* pOrgYuv, TComYuv* pDstYuv, UInt uiMaxSearchRounds, UInt nSearchStepShift, UInt& uiMinCost
#if DMVR_MRSAD_SPEEDUP
      , DistParam * cDistParam
#endif
#if ITH_C2_MC_OPT
      , TComMv *refineMv
#endif
#if DMVR_IMP_PHASE1
, TComMv  &SecondBestMv, UInt& uiSecondBestCost, TComMv mvDiff = TComMv(0,0), bool bFracPel = 0
#endif
);
  UInt xDirectMCCost(Int iBitDepth, Pel* pRef, UInt uiRefStride, Pel* pOrg, UInt uiOrgStride, Int iWidth, Int iHeight
#if DMVR_MRSAD_SPEEDUP
    , DistParam * cDistParam
#endif
  );
  Void xPredInterLines(TComDataCU *cu, TComPicYuv *refPic, UInt partAddr, TComMv *mv, Int width, Int height, Pel* dstPix, Int dstStride, Bool bi, const Int bitDepth);
  Void xFillPredBorder(TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefPicList, Int iWidth, Int iHeight, TComYuv* pDstYuv);
#if DMVR_HALF_ME
  Void xGenerateFracPixel(TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefPicList, Int iWidth, Int iHeight, UInt nSearchStepShift
#if (ITH_C2_NUM_HALF_POINT == 2) && !ITH_NON_NORMATIVE_DISABLE
    , TComMv mvDiff = TComMv(0, 0)
#endif
  );
#endif
#endif
#if RRC_SELECTABLE_IF
protected:
#endif
  Void xWeightedAverage         ( TComDataCU * pCu, TComYuv* pcYuvSrc0, TComYuv* pcYuvSrc1, Int iRefIdx0, Int iRefIdx1, UInt uiPartAddr, Int iWidth, Int iHeight, TComYuv* pcYuvDst, const BitDepths &clipBitDepths  
#if VCEG_AZ05_BIO                  
    , Bool bBIOapplied 
#endif
#if JVET_E0052_DMVR
    , Bool bRefineflag
#endif
#if JVET_E0052_DMVR || JVET_G0082
    , Bool bOBMC
#endif
    );

  Void xGetLLSPrediction ( const Pel* pSrc0, Int iSrcStride, Pel* pDst0, Int iDstStride, UInt uiWidth, UInt uiHeight, UInt uiExt0, const ChromaFormat chFmt  DEBUG_STRING_FN_DECLARE(sDebug) );
#if COM16_C1046_PDPC_INTRA
  Void xReferenceFilter  (Int iBlkSize, Int iOrigWeight, Int iFilterOrder, Int * piRefVector, Int * piLowPassRef);
#endif
  Void xDCPredFiltering( const Pel* pSrc, Int iSrcStride, Pel* pDst, Int iDstStride, Int iWidth, Int iHeight, ChannelType channelType );

#if VCEG_AZ05_INTRA_MPI
  Void xMPIredFiltering(Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, Int idxMPI);
#endif

#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
  Void xIntraPredFilteringModeDGL(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, UInt uiMode
#if DIMD_NUM_INTRA_DIR_INC
    , Bool bExtIntraDir
#endif
    );
  Void xIntraPredFilteringMode34(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight);
  Void xIntraPredFilteringMode02(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight);
#endif
  Bool xCheckIdenticalMotion(TComDataCU* pcCU, UInt PartAddr);
  Void destroy();

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  Bool xCheckTwoSPMotion ( TComDataCU* pcCU, UInt PartAddr0, UInt PartAddr1 );
  Void xGetSubPUAddrAndMerge(TComDataCU* pcCU, UInt uiPartAddr, Int iSPWidth, Int iSPHeight, Int iNumSPInOneLine, Int iNumSP, UInt* uiMergedSPW, UInt* uiMergedSPH, UInt* uiSPAddr );
#endif
#if COM16_C806_OBMC
  Void xSubblockOBMC ( const ComponentID eComp, TComDataCU* pcCU, Int uiAbsPartIdx, TComYuv* pcYuvPredDst, TComYuv* pcYuvPredSrc, Int iWidth, Int iHeight, Int iDir, Bool bOBMCSimp );
  Void xSubtractOBMC ( TComDataCU* pcCU, Int uiAbsPartIdx, TComYuv* pcYuvPredDst, TComYuv* pcYuvPredSrc, Int iWidth, Int iHeight, Int iDir, Bool bOBMCSimp );
  Void xSubBlockMotionCompensation ( TComDataCU* pcCU, TComYuv* pcYuvPred, Int uiPartAddr, Int iWidth, Int iHeight  
#if JVET_E0052_DMVR
    , Bool bRefineflag
#endif
    );
#endif
#if VCEG_AZ07_FRUC_MERGE
  Bool xFrucFindBlkMv( TComDataCU * pCU , UInt uiPUIdx 
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
      , Bool bTemplateModeSearched = true
#endif
      );
  Bool xFrucFindBlkMv4Pred( TComDataCU * pCU , UInt uiPUIdx , RefPicList eTargetRefPicList , Int nTargetRefIdx
#if JVET_E0060_FRUC_CAND
                          , AMVPInfo* pInfo = NULL
#endif
                          );
  Bool xFrucRefineSubBlkMv( TComDataCU * pCU , UInt uiDepth , UInt uiPUIdx , Bool bTM 
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
      , Bool bTemplateModeSearched = true
#endif
      );

  Void xFrucCollectBlkStartMv( TComDataCU * pCU , UInt uiPUIdx , RefPicList eTargetRefList = REF_PIC_LIST_0 , Int nTargetRefIdx = -1
#if JVET_E0060_FRUC_CAND
                             , AMVPInfo* pInfo = NULL
#endif
                             );
  Void xFrucCollectSubBlkStartMv( TComDataCU * pCU , UInt uiAbsPartIdx , RefPicList eRefPicList , const TComMvField & rMvStart , Int nSubBlkWidth , Int nSubBlkHeight 
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    , UInt uiSubBlkRasterIdx , UInt uiSubBlkRasterStep
#endif
    );

  UInt xFrucFindBestMvFromList( TComMvField * pBestMvField , RefPicList & rBestRefPicList , TComDataCU * pCU , UInt uiAbsPartIdx , const TComMvField & rMvStart , Int nBlkWidth , Int nBlkHeight , Bool bTM , Bool bMvCost );

  UInt xFrucRefineMv( TComMvField * pBestMvField , RefPicList eCurRefPicList , UInt uiMinCost , Int nSearchMethod , TComDataCU * pCU , UInt uiAbsPartIdx , const TComMvField & rMvStart , Int nBlkWidth , Int nBlkHeight , Bool bTM
#if JVET_F0032_UNI_BI_SELECTION
      , Bool bMvCostZero = false
#endif
  );
  template<Int SearchPattern>
  UInt xFrucRefineMvSearch( TComMvField * pBestMvField , RefPicList eCurRefPicList , TComDataCU * pCU , UInt uiAbsPartIdx , TComMvField const & rMvStart , Int nBlkWidth , Int nBlkHeight , UInt uiMinDist , Bool bTM , Int nSearchStepShift , UInt uiMaxSearchRounds = MAX_UINT 
#if JVET_F0032_UNI_BI_SELECTION
      , Bool bMvCostZero = false
#endif
  );

  UInt xFrucGetMvCost( const TComMv & rMvStart , const TComMv & rMvCur , Int nSearchRange , Int nWeighting );
  UInt xFrucGetBilaMatchCost( TComDataCU * pcCU , UInt uiAbsPartIdx , Int nWidth , Int nHeight , RefPicList eCurRefPicList , const TComMvField & rCurMvField , TComMvField & rPairMVField , UInt uiMVCost );
  UInt xFrucGetTempMatchCost( TComDataCU * pcCU , UInt uiAbsPartIdx , Int nWidth , Int nHeight , RefPicList eCurRefPicList , const TComMvField & rCurMvField , UInt uiMVCost );
#if JVET_F0032_UNI_BI_SELECTION
  Void xFrucUpdateTemplate(TComDataCU * pcCU, UInt uiAbsPartIdx, Int nWidth, Int nHeight, RefPicList eCurRefPicList, TComMvField rCurMvField);
#endif

  Void xFrucInsertMv2StartList( const TComMvField & rMvField , std::list<TComMvField> & rList );
  Bool xFrucIsInList( const TComMvField & rMvField , std::list<TComMvField> & rList );

  Bool xFrucGetCurBlkTemplate( TComDataCU * pCU , UInt uiAbsPartIdx , Int nCurBlkWidth , Int nCurBlkHeight );
  Bool xFrucIsTopTempAvailable( TComDataCU * pCU , UInt uiAbsPartIdx );
  Bool xFrucIsLeftTempAvailable( TComDataCU * pCU , UInt uiAbsPartIdx );
  Int  xFrucGetSubBlkSize( TComDataCU * pcCU , UInt uiAbsPartIdx , Int nBlkWidth , Int nBlkHeight );
#endif
#if VCEG_AZ06_IC
  Void xGetLLSICPrediction( TComDataCU* pcCU, TComMv *pMv, TComPicYuv *pRefPic, Int &a, Int &b, const ComponentID eComp, Int nBitDepth );
#endif
public:

#if RRC_DWDIP          
  RSGIntraParams rsgIntraParams; //Reference samples generation parameters for intra-prediction
  RSGIntraParams & rrsgIntraParams; // Effective reference samples generation parameters being used
#endif

  TComPrediction();
  virtual ~TComPrediction();
#if COM16_C806_OBMC
  Void subBlockOBMC ( TComDataCU*  pcCU, UInt uiAbsPartIdx, TComYuv *pcYuvPred, TComYuv *pcYuvTmpPred1, TComYuv *pcYuvTmpPred2
#if JVET_E0052_DMVR
    , Bool bRefineflag = true
#endif
    , Bool bOBMC4ME = false );
#endif
#if COM16_C806_LMCHROMA
  Void    initTempBuff(ChromaFormat chromaFormatIDC, Int bitDepthY
#if VCEG_AZ08_INTER_KLT
    , bool interKLT , const Int iPicWidth, const Int iPicHeight, const UInt uiMaxCUWidth, const UInt uiMaxCUHeight, const UInt uiMaxCUDepth
#endif
    );
#else
  Void    initTempBuff(ChromaFormat chromaFormatIDC
#if VCEG_AZ08_INTER_KLT
    , bool interKLT , const Int iPicWidth, const Int iPicHeight, const UInt uiMaxCUWidth, const UInt uiMaxCUHeight, const UInt uiMaxCUDepth
#endif
    );
#endif

  ChromaFormat getChromaFormat() const { return m_cYuvPredTemp.getChromaFormat(); }

  // inter
  Void motionCompensation         ( TComDataCU*  pcCU, TComYuv* pcYuvPred
#if JVET_E0052_DMVR
    , Bool bRefineflag = true
#endif
    , RefPicList eRefPicList = REF_PIC_LIST_X, Int iPartIdx = -1 );

#if VCEG_AZ07_FRUC_MERGE
  Bool deriveFRUCMV( TComDataCU * pCU , UInt uiDepth , UInt uiAbsPartIdx , UInt uiPUIdx ,   Int nTargetRefIdx = -1 , RefPicList eTargetRefList = REF_PIC_LIST_0
#if JVET_E0060_FRUC_CAND
                   , AMVPInfo* pInfo = NULL
#endif
#if THIRD_FRUC_MODE_ENCODER_SPEEDUP
      ,Bool bTemplateModeSearched = 
#if THIRD_FRUC_MODE_TT_OFF_FIX
      false
#else
      true
#endif
#endif
                   );
#endif

  // motion vector prediction
  Void getMvPredAMVP              ( TComDataCU* pcCU, UInt uiPartIdx, UInt uiPartAddr, RefPicList eRefPicList, TComMv& rcMvPred );

  // Angular Intra
  Void predIntraAng               ( const ComponentID compID, UInt uiDirMode, Pel *piOrg /* Will be null for decoding */, UInt uiOrgStride, Pel* piPred, UInt uiStride, TComTU &rTu, const Bool bUseFilteredPredSamples
                                  , const Bool bUseLosslessDPCM = false );

#if JVET_E0077_MMLM
  struct MMLM_parameter
  {
      Int Inf;  // Inferio boundary
      Int Sup;  // Superior bounday
      Int a;
      Int b;
      Int shift;
  };
  Int xCalcLMParametersGeneralized(Int x, Int y, Int xx, Int xy, Int iCountShift, Int bitDepth, Int &a, Int &b, Int &iShift);
  Int xLMSampleClassifiedTraining(Int count, Int LumaSamples[], Int ChrmSamples[], Int GroupNum, Int bitDepth, MMLM_parameter parameters[]);
  Int xGetMMLMParameters(TComTU& rTu, const ComponentID compID, UInt uiWidth, UInt uiHeight, Int &numClass, MMLM_parameter parameters[]);
#endif



#if COM16_C806_LMCHROMA
#if JVET_E0077_LM_MF
  Void xFilterGroup(Pel* pMulDst[LM_FILTER_NUM], Int i, Pel* piSrc, Int iRecStride, Bool bAboveAvaillable, Bool bLeftAvaillable);
#endif

  Void predLMIntraChroma(TComTU& rTu, ComponentID compID, Pel* pPred, UInt uiPredStride, UInt uiCWidth, UInt uiCHeight
#if JVET_E0077_ENHANCED_LM
      , Int LMtype = LM_CHROMA_IDX
#endif
      );
  Void getLumaRecPixels  ( TComTU& rTu, UInt uiCWidth, UInt uiCHeight );
  Void addCrossColorResi ( TComTU& rTu, ComponentID compID, Pel* piPred, UInt uiPredStride, UInt uiWidth, UInt uiHeight, Pel* piResi, UInt uiResiStride );
  Void xGetLMParameters  ( TComTU& rTu,  ComponentID compID, UInt uiWidth, UInt uiHeight, Int iPredType, Int &a, Int &b, Int &iShift );
  Void xCalcLMParameters ( Int x, Int y, Int xx, Int xy, Int iCountShift, Int iPredType, Int bitDepth, Int &a, Int &b, Int &iShift );
#endif

#if COM16_C1016_AFFINE
  Bool xCheckIdenticalAffineMotion ( TComDataCU* pcCU, UInt PartAddr, Int iWidth, Int iHeight );
  Void xPredAffineBlk              ( const ComponentID compID, TComDataCU *cu, TComPicYuv *refPic, UInt partAddr, TComMv acMv[3], Int width, Int height, TComYuv *dstPic, Bool bi, const Int bitDepth 
#if JVECT_C0062_AFFINE_SIX_PARAM
    , Bool b6Param
#endif
  );

  Void getMvPredAffineAMVP         ( TComDataCU* pcCU, UInt uiPartIdx, UInt uiPartAddr, RefPicList eRefPicList, TComMv acMvPred[3]
#if JVECT_C0062_AFFINE_SIX_PARAM
    , Bool b6Param
#endif
  );
#if RRC_WIF
  Void xPredAffineNewBlk           ( const ComponentID compID, TComDataCU *cu, TComPicYuv *refPic, UInt partAddr, TComMv acMv[3], Int width, Int height, TComYuv *dstPic, Bool bi, const Int bitDepth );
#endif
#endif

  Pel  predIntraGetPredValDC      ( const Pel* pSrc, Int iSrcStride, UInt iWidth, UInt iHeight);

  Pel*  getPredictorPtr           ( const ComponentID compID, const Bool bUseFilteredPredictions )
  {
    return m_piYuvExt[compID][bUseFilteredPredictions ? PRED_BUF_FILTERED : PRED_BUF_UNFILTERED];
  }
#if MLIP_EXTEND_REF_LINE
  Pel*  getSecondPredictorPtr(const ComponentID compID, const Bool bUseFilteredPredictions)
  {
      return m_piSecondYuvExt[compID][bUseFilteredPredictions ? PRED_BUF_FILTERED : PRED_BUF_UNFILTERED];
  }
#endif
#if BARC_FUNC
  Pel*  getDscalePredictorPtr     ( const ComponentID compID, const Bool bUseFilteredPredictions )
  {
    return m_piYuvExtDscale[compID][bUseFilteredPredictions?PRED_BUF_FILTERED:PRED_BUF_UNFILTERED];
  }
#endif
#if DIMD_INTRA_PRED
  UInt deriveNeighborIntraDirs(TComDataCU* pcCU, UInt uiAbsPartIdx);
  Void initDIMDLumaFlexibleIntraPattern(TComDataCU* pcCU,
    UInt        uiZorderIdxInPart,
    Int         iRefX,
    Int         iRefY,
    UInt        uiRefWidth,
    UInt        uiRefHeight,
    Bool        bFilterRefSamples
#if BARC_DIMD
    , Bool isBARC = false, Bool isInCurCtu = false, UInt uiRefZorderIdxInPart = 0
#endif
    );
  static Bool filteringDIMDIntraReferenceSamples(const ComponentID compID, UInt uiDirMode, UInt uiTuChWidth, UInt uiTuChHeight, const ChromaFormat chFmt, const Bool intraReferenceSmoothingDisabled );
  Void predDIMDIntraLumaAng(TComDataCU* pcCU, UInt uiDirMode, Pel* piPred, UInt uiStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
  Void xPredDIMDIntraPlanar(const Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
  Void xPredDIMDIntraAng(Int bitDepth, const Pel* pSrc, Int srcStride, Pel* rpDst, Int dstStride, UInt width, UInt height, ChannelType channelType, UInt dirMode, const Bool bEnableEdgeFilters
#if VCEG_AZ07_INTRA_4TAP_FILTER
    , Bool          enable4TapFilter
#endif
    , TEMPLATE_TYPE eTempType
    , Int           iTemplateWidth
    , Int           iTemplateHeight
    );
#if COM16_C1046_PDPC_INTRA
  Void predDIMDIntraLumaPDPC(TComDataCU* pcCU, UInt uiDirMode, Pel* piPred, UInt uiStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
#endif
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
  Void xDIMDIntraPredFilteringMode34(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
  Void xDIMDIntraPredFilteringMode02(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
  Void xDIMDIntraPredFilteringModeDGL(const Pel* pSrc, Int iSrcStride, Pel*& rpDst, Int iDstStride, Int iWidth, Int iHeight, UInt uiMode, TEMPLATE_TYPE eTempType, Int iTemplateWidth, Int iTemplateHeight);
#endif
  UInt calcTemplateSAD(Int bitDepth, Pel* pi0, Int iStride0, Pel* pi1, Int iStride1, Int iWidth, Int iHeight, Bool bBDClip);
#endif

  // This function is actually still in TComPattern.cpp
  /// set parameters from CU data for accessing intra data
  Void initIntraPatternChType ( TComTU &rTu,
                              const ComponentID compID, const Bool bFilterRefSamples
#if RRC_INTRA_SHARPENING
                              , Bool bSharpRefSamples = false
#endif
#if BARC
                               , Bool isBarc = false
#endif
                              DEBUG_STRING_FN_DECLARE(sDebug)
                              );

  static Bool filteringIntraReferenceSamples(const ComponentID compID, UInt uiDirMode, UInt uiTuChWidth, UInt uiTuChHeight, const ChromaFormat chFmt, const Bool intraReferenceSmoothingDisabled
#if DIMD_NUM_INTRA_DIR_INC
  , TComDataCU* pcCU, UInt uiAbsPartIdx
#endif
                                            );

  static Bool UseDPCMForFirstPassIntraEstimation(TComTU &rTu, const UInt uiDirMode);
#if VCEG_AZ08_INTER_KLT
  Void interpolatePic(TComPic* pcPic);
#endif

#if RRC_INTRA_SHARPENING
  typedef Short PelMap;
  static Bool initSharpPattern(Pel* piAdiBuf, Pel* piAdiBufSharp, Int iStride, Int iWidth, Int iHeight, 
    TComTU &rTu, TComDataCU *pcCU, UInt uiAbsPartIdx,
    PelMap *tmpArrSource, Int tmpArrSourceSize,
    PelMap *tmpArrMap, Int tmpArrMapSize,
    PelMap *tmpArrMap2, Int tmpArrMap2Size,
    Int bitDepth
    );

protected:
  static Void calcEdgePattern(PelMap *pSamples, PelMap *pEdgeMap, PelMap *pTmpArr, int iLen, Int bitDepth, Int qp);

  static void displaceVector(PelMap* pcPicSrc, Pel* pcPicDst, Int iStride, PelMap *pcMapX,
    int map_shift, int bitDepth, Int k, int block_width, int clipOut = 0, int margin = 0);
#endif
};

//! \}

#endif // __TCOMPREDICTION__

