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

/** \file     TComPic.h
    \brief    picture class (header)
*/

#ifndef __TCOMPIC__
#define __TCOMPIC__

// Include files
#include "CommonDef.h"
#include "TComPicSym.h"
#include "TComPicYuv.h"
#include "TComBitStream.h"
#if BARC_REF_UP_SAMPLE_CNN
#include "TComConvNet.h"
#endif
//! \ingroup TLibCommon
//! \{

#if RRC_SELECTABLE_IF
class TComSharpeningFilter;
#endif
// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// picture class (symbol + YUV buffers)

class TComPic
{
public:
  typedef enum { PIC_YUV_ORG=0, PIC_YUV_REC=1, PIC_YUV_TRUE_ORG=2, 
#if BARC
#if COMPACT_RESOLUTION_CNN
	  PIC_YUV_Dscale = 3,
	  PIC_YUV_Dscale_REF = 4,
	  //PIC_YUV_REC_TMP = 5,
	  //PIC_YUV_CRCNN = 6,
	  //NUM_PIC_YUV = 7
	  PIC_YUV_CRCNN = 5,
	  NUM_PIC_YUV = 6
#else
	  PIC_YUV_Dscale = 3,
	  PIC_YUV_Dscale_REF = 4,
	  NUM_PIC_YUV = 5
#endif
#else
    NUM_PIC_YUV=3
#endif
  } PIC_YUV_T;
     // TRUE_ORG is the input file without any pre-encoder colour space conversion (but with possible bit depth increment)
  TComPicYuv*   getPicYuvTrueOrg()        { return  m_apcPicYuv[PIC_YUV_TRUE_ORG]; }

private:
  UInt                  m_uiTLayer;               //  Temporal layer
  Bool                  m_bUsedByCurr;            //  Used by current picture
  Bool                  m_bIsLongTerm;            //  IS long term picture
  TComPicSym            m_picSym;                 //  Symbol
  TComPicYuv*           m_apcPicYuv[NUM_PIC_YUV];

#if BARC_REF_UP_SAMPLE_CNN
  TComConvNet*           m_cCNNUpSampleImg;
#if COMPACT_RESOLUTION_CNN
  TComConvNet*           m_cCNNDownSampleImg;
#endif
#endif
  TComPicYuv*           m_pcPicYuvPred;           //  Prediction
  TComPicYuv*           m_pcPicYuvResi;           //  Residual
  Bool                  m_bReconstructed;
  Bool                  m_bNeededForOutput;
  UInt                  m_uiCurrSliceIdx;         // Index of current slice
  Bool                  m_bCheckLTMSB;

  Bool                  m_isTop;
  Bool                  m_isField;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  Int                   m_iBaseUnitWidth;       ///< Width of Base Unit (with maximum depth or minimum size, m_iCuWidth >> Max. Depth)
  Int                   m_iBaseUnitHeight;      ///< Height of Base Unit (with maximum depth or minimum size, m_iCuHeight >> Max. Depth)
  Int                   m_iNumCuInWidth;
#endif
#if JVET_C0024_QTBT
  //for record codec block info.
  Bool m_bCodedBlkInCTU[MAX_CU_SIZE>>MIN_CU_LOG2][MAX_CU_SIZE>>MIN_CU_LOG2];    //[CTUSize>>MIN_CU_Log2][CTUSize>>MIN_CU_Log2]; [h][w]
  Int  m_iCodedArea;

  //for encoder speedup
  TComMv                m_cIntMv[1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1)][MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1][2][5]; //[zorder][w][h][refList][refIdx]
  Bool                  m_bSetIntMv[1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1)][MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1][2][5]; //[zorder][w][h][refList][refIdx]
  Bool                  m_bSkiped[1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1)][MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; //[zorder][w][h] , if skip mode, not try inter, intra
  Bool                  m_bInter[1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1)][MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; //[zorder][w][h] , if inter mode, not try intra
  Bool                  m_bIntra[1<<((MAX_CU_DEPTH-MIN_CU_LOG2)<<1)][MAX_CU_DEPTH-MIN_CU_LOG2+1][MAX_CU_DEPTH-MIN_CU_LOG2+1]; // if intra mode, not try inter
#endif

  std::vector<std::vector<TComDataCU*> > m_vSliceCUDataLink;

#if BLOCK_BASED_BACKGROUND_REFERENCE
	Int*                   m_SpliceIdx;
#endif

  SEIMessages  m_SEIs; ///< Any SEI messages that have been received.  If !NULL we own the object.

#if RRC_SELECTABLE_IF
  TComSharpeningFilter*  m_pSharpeningFilter;
#endif
public:
  TComPic();
  virtual ~TComPic();

#if JVET_C0024_QTBT
  //to record coded block info.
  Void          setCodedBlkInCTU(Bool bCoded, UInt uiBlkX, UInt uiBlkY, UInt uiWidth, UInt uiHeight);
  Bool          getCodedBlkInCTU(UInt uiBlkX, UInt uiBlkY) {return m_bCodedBlkInCTU[uiBlkY][uiBlkX];}
  Void          setCodedAreaInCTU(Int iArea);
  Void          addCodedAreaInCTU(Int iArea);
  Int           getCodedAreaInCTU();

  //for encoder speed-up
  Void          setSkiped(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bSkip);
  Bool          getSkiped(UInt uiZorder, UInt uiWidth, UInt uiHeight);
  Void          clearAllSkiped();
  Void          setInter(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bInter);
  Bool          getInter(UInt uiZorder, UInt uiWidth, UInt uiHeight);
  Void          clearAllInter();
  Void          setIntra(UInt uiZorder, UInt uiWidth, UInt uiHeight, Bool bIntra);
  Bool          getIntra(UInt uiZorder, UInt uiWidth, UInt uiHeight);
  Void          clearAllIntra();

  Void          setIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx, TComMv cMv);
  TComMv        getIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx);
  Void          clearAllIntMv();
  Bool          IsSetIntMv(UInt uiZorder, UInt uiWidth, UInt uiHeight, RefPicList eRefList, UInt uiRefIdx);
#endif

  Void          create( const TComSPS &sps, const TComPPS &pps, const Bool bIsVirtual /*= false*/
#if BARC_REF_UP_SAMPLE_CNN
    , TComConvNet*           cCNNUpSampleImg = NULL
#if COMPACT_RESOLUTION_CNN
    , TComConvNet*           cCNNDownSampleImg = NULL
#endif
#endif
  );

#if BLOCK_BASED_BACKGROUND_REFERENCE
	Int           getSpliceIdx(UInt idx) { return m_SpliceIdx[idx]; }
	Void          setSpliceIdx(UInt idx, Int poc) { m_SpliceIdx[idx] = poc; }
	Void          createSpliceIdx();
#endif

  virtual Void  destroy();

  UInt          getTLayer() const               { return m_uiTLayer;   }
  Void          setTLayer( UInt uiTLayer ) { m_uiTLayer = uiTLayer; }

  Bool          getUsedByCurr() const            { return m_bUsedByCurr; }
  Void          setUsedByCurr( Bool bUsed ) { m_bUsedByCurr = bUsed; }
  Bool          getIsLongTerm() const            { return m_bIsLongTerm; }
  Void          setIsLongTerm( Bool lt ) { m_bIsLongTerm = lt; }
  Void          setCheckLTMSBPresent     (Bool b ) {m_bCheckLTMSB=b;}
  Bool          getCheckLTMSBPresent     () { return m_bCheckLTMSB;}

#if VCEG_AZ07_FRUC_MERGE
  Void          initFRUCMVP()         { return m_picSym.initFRUCMVP();  }
#endif
  TComPicSym*   getPicSym()           { return  &m_picSym;    }
  TComSlice*    getSlice(Int i)       { return  m_picSym.getSlice(i);  }
  Int           getPOC() const        { return  m_picSym.getSlice(m_uiCurrSliceIdx)->getPOC();  }
  TComDataCU*   getCtu( UInt ctuRsAddr )           { return  m_picSym.getCtu( ctuRsAddr ); }
  const TComDataCU* getCtu( UInt ctuRsAddr ) const { return  m_picSym.getCtu( ctuRsAddr ); }

  TComPicYuv*   getPicYuvOrg()        { return  m_apcPicYuv[PIC_YUV_ORG]; }
  TComPicYuv*   getPicYuvRec()        { return  m_apcPicYuv[PIC_YUV_REC]; }
#if BARC
  TComDataCU*   getCtuTmp()                 { return m_picSym.getCtuTmp();     }
  TComPicYuv*   getPicYuvDScaleOrg()        { return  m_apcPicYuv[PIC_YUV_Dscale]; }
  TComPicYuv*   getPicYuvDScaleRef()        { return  m_apcPicYuv[PIC_YUV_Dscale_REF]; }
#if COMPACT_RESOLUTION_CNN
  TComDataCU*   getCtuTmpBicubic() { return m_picSym.getCtuTmpBicubic(); }
  TComPicYuv*   getPicYuvCRCNNOrg() { return  m_apcPicYuv[PIC_YUV_CRCNN]; }
  Void          setPicYuvDScaleOrg(TComPicYuv* pcPicYuv) { m_apcPicYuv[PIC_YUV_Dscale] = pcPicYuv; }
#endif
#endif
  TComPicYuv*   getPicYuvPred()       { return  m_pcPicYuvPred; }
  TComPicYuv*   getPicYuvResi()       { return  m_pcPicYuvResi; }
  Void          setPicYuvPred( TComPicYuv* pcPicYuv )       { m_pcPicYuvPred = pcPicYuv; }
  Void          setPicYuvResi( TComPicYuv* pcPicYuv )       { m_pcPicYuvResi = pcPicYuv; }

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  Void          getCUAddrAndPartIdx( Int iX, Int iY, Int& riCuAddr, Int& riAbsZorderIdx );
#endif

#if ALF_HM3_REFACTOR
  TComPicYuv*   replacePicYuvRecPointer( TComPicYuv * pNewPicYuvRec )
  {
    TComPicYuv * pOldPicYuvRec = m_apcPicYuv[1];
    m_apcPicYuv[1] = pNewPicYuvRec;
    return ( pOldPicYuvRec );
  }
#endif

  UInt          getNumberOfCtusInFrame() const     { return m_picSym.getNumberOfCtusInFrame(); }
  UInt          getNumPartInCtuWidth() const       { return m_picSym.getNumPartInCtuWidth();   }
  UInt          getNumPartInCtuHeight() const      { return m_picSym.getNumPartInCtuHeight();  }
  UInt          getNumPartitionsInCtu() const      { return m_picSym.getNumPartitionsInCtu();  }
  UInt          getFrameWidthInCtus() const        { return m_picSym.getFrameWidthInCtus();    }
  UInt          getFrameHeightInCtus() const       { return m_picSym.getFrameHeightInCtus();   }
  UInt          getMinCUWidth() const              { return m_picSym.getMinCUWidth();          }
  UInt          getMinCUHeight() const             { return m_picSym.getMinCUHeight();         }

  Int           getStride(const ComponentID id) const          { return m_apcPicYuv[PIC_YUV_REC]->getStride(id); }
  Int           getComponentScaleX(const ComponentID id) const    { return m_apcPicYuv[PIC_YUV_REC]->getComponentScaleX(id); }
  Int           getComponentScaleY(const ComponentID id) const    { return m_apcPicYuv[PIC_YUV_REC]->getComponentScaleY(id); }
  ChromaFormat  getChromaFormat() const                           { return m_apcPicYuv[PIC_YUV_REC]->getChromaFormat(); }
  Int           getNumberValidComponents() const                  { return m_apcPicYuv[PIC_YUV_REC]->getNumberValidComponents(); }

  Void          setReconMark (Bool b) { m_bReconstructed = b;     }
  Bool          getReconMark () const      { return m_bReconstructed;  }
  Void          setOutputMark (Bool b) { m_bNeededForOutput = b;     }
  Bool          getOutputMark () const      { return m_bNeededForOutput;  }

  Void          compressMotion();
  UInt          getCurrSliceIdx() const           { return m_uiCurrSliceIdx;                }
  Void          setCurrSliceIdx(UInt i)      { m_uiCurrSliceIdx = i;                   }
  UInt          getNumAllocatedSlice() const      {return m_picSym.getNumAllocatedSlice();}
  Void          allocateNewSlice()           {m_picSym.allocateNewSlice();         }
  Void          clearSliceBuffer()           {m_picSym.clearSliceBuffer();         }
#if BARC_REF_UP_SAMPLE_CNN
  Void          initCNNModel();
  TComConvNet*  getConvNetUpSample() { return m_cCNNUpSampleImg; }
  Void          CNNUpScaleLumaComponent2x2(Bool uiFlag, const Pel* pSrc, Double* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight);
#if COMPACT_RESOLUTION_CNN
  TComConvNet*  getConvNetDownSample() { return m_cCNNDownSampleImg; }
  Void          downScaleOrgCRCNN();
  Void          CNNDownScaleLumaComponent2x2(Bool uiFlag, const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight);
#endif
#endif
#if REFINE_BLOCK_BOUNDARY_NET
  Void          upSampleDScaleRecon();
#endif
  const Window& getConformanceWindow() const { return m_picSym.getSPS().getConformanceWindow(); }
  Window        getDefDisplayWindow() const  { return m_picSym.getSPS().getVuiParametersPresentFlag() ? m_picSym.getSPS().getVuiParameters()->getDefaultDisplayWindow() : Window(); }

  Bool          getSAOMergeAvailability(Int currAddr, Int mergeAddr);

  UInt          getSubstreamForCtuAddr(const UInt ctuAddr, const Bool bAddressInRaster, TComSlice *pcSlice);

#if RRC_SELECTABLE_IF
  TComSharpeningFilter*& getSharpeningFilter() { return m_pSharpeningFilter; }
#endif
  /* field coding parameters*/

   Void              setTopField(Bool b)                  {m_isTop = b;}
   Bool              isTopField()                         {return m_isTop;}
   Void              setField(Bool b)                     {m_isField = b;}
   Bool              isField()                            {return m_isField;}

  /** transfer ownership of seis to this picture */
  Void setSEIs(SEIMessages& seis) { m_SEIs = seis; }

  /**
   * return the current list of SEI messages associated with this picture.
   * Pointer is valid until this->destroy() is called */
  SEIMessages& getSEIs() { return m_SEIs; }

  /**
   * return the current list of SEI messages associated with this picture.
   * Pointer is valid until this->destroy() is called */
  const SEIMessages& getSEIs() const { return m_SEIs; }

#if VCEG_AZ08_INTER_KLT
  TComPicYuv*   m_apcQuaPicYuv[4][4];     //   quarter pixel reconstructed pictures (fractional pixels); 
  TComPicYuv*   getPicQuaYuvRec(UInt uiRow, UInt uiCol)     { return  m_apcQuaPicYuv[uiRow][uiCol]; }
#endif
#if JVET_D0033_ADAPTIVE_CLIPPING
  ClipParam m_aclip_prm;
  ClipParam computeTchClipParam(Int &delta_disto_luma, Int &delta_disto_chroma) const;
#endif
};// END CLASS DEFINITION TComPic

//! \}

#endif // __TCOMPIC__
