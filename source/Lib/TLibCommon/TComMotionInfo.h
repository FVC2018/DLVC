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

/** \file     TComMotionInfo.h
    \brief    motion information handling classes (header)
    \todo     TComMvField seems to be better to be inherited from TComMv
*/

#ifndef __TCOMMOTIONINFO__
#define __TCOMMOTIONINFO__

#include <memory.h>
#include "CommonDef.h"
#include "TComMv.h"

#if MRG_CONSIDER_COMPLEX
#include "TComRom.h"
#include "math.h"
#include "TComSlice.h"
#endif

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP || VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
class TComDataCU;
class TComPic;
#endif
//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Type definition
// ====================================================================================================================

/// parameters for AMVP
typedef struct _AMVPInfo
{
  TComMv m_acMvCand[ AMVP_MAX_NUM_CANDS_MEM ];  ///< array of motion vector predictor candidates
  Int    iN;                                ///< number of motion vector predictor candidates
} AMVPInfo;

#if COM16_C1016_AFFINE
typedef struct _AffineAMVPInfo
{
  TComMv m_acMvCand[ AMVP_MAX_NUM_CANDS_MEM ][3];  ///< array of affine motion vector predictor candidates
  Int    iN;                                       ///< number of affine motion vector predictor candidates
} AffineAMVPInfo;
#endif

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// class for motion vector with reference index
class TComMvField
{
private:
  TComMv    m_acMv;
  Int       m_iRefIdx;

public:
  TComMvField() : m_iRefIdx( NOT_VALID ) {}

  Void setMvField( TComMv const & cMv, Int iRefIdx )
  {
    m_acMv    = cMv;
    m_iRefIdx = iRefIdx;
  }

  Void setRefIdx( Int refIdx ) { m_iRefIdx = refIdx; }

  TComMv const & getMv() const { return  m_acMv; }
  TComMv       & getMv()       { return  m_acMv; }

  Int getRefIdx() const { return  m_iRefIdx;       }
  Int getHor   () const { return  m_acMv.getHor(); }
  Int getVer   () const { return  m_acMv.getVer(); }
#if COM16_C806_GEN_MRG_IMPROVEMENT || COM16_C806_OBMC || VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
  Bool operator== (const TComMvField& rcMv) const
  {
    return (m_acMv.getHor()== rcMv.getHor() && m_acMv.getVer()== rcMv.getVer() && m_iRefIdx == rcMv.getRefIdx());
  }
#endif
};

#if MRG_CONSIDER_COMPLEX
typedef struct _NodeCand
{
  TComMvField MvFieldNbNodes[MRG_NODE_NUM<<1];//motion information of each node
  UChar InterDirNbNodes[MRG_NODE_NUM ];
  Int numValidNodes[MRG_COMPLEX_NUM_CANDS<<1];//numValidMergeNode[i] represents number of valid merge nodes of list i  (i=0,1) after removing duplication,
  Int RealComplexIdx[MRG_COMPLEX_NUM_CANDS];
} NodeInfo;
#endif

/// class for motion information in one CU
class TComCUMvField
{
private:
  TComMv*   m_pcMv;
  TComMv*   m_pcMvd;
  Char*     m_piRefIdx;
  UInt      m_uiNumPartition;
  AMVPInfo  m_cAMVPInfo;
#if COM16_C1016_AFFINE
#if JVECT_C0062_AFFINE_SIX_PARAM
  AffineAMVPInfo m_cAffineAMVPInfo[2];
#else
  AffineAMVPInfo m_cAffineAMVPInfo;
#endif
#endif
#if JVET_C0024_QTBT
  TComDataCU* m_pcCU;
#endif

#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
  UInt      m_uiWidth;
  UInt      m_uiHeight;
#endif

  template <typename T>
  Void setAll( T *p, T const & val, PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx );

#if MRG_CONSIDER_COMPLEX
  Int  agaus ( Double *a,Double *b,Int n );
  Void setAll( TComMv*  pMv, Char* pRefIdx, NodeInfo NodeMotionInfo, UInt MergeComplexIndex, PartSize eCUMode, RefPicList eRefPicList, Int iPartAddr, UInt uiDepth, Int iPartIdx );
#endif

public:
  TComCUMvField() : m_pcMv(NULL), m_pcMvd(NULL), m_piRefIdx(NULL), m_uiNumPartition(0) {}
  ~TComCUMvField() {}

  // ------------------------------------------------------------------------------------------------------------------
  // create / destroy
  // ------------------------------------------------------------------------------------------------------------------

#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
  Void    create( UInt uiNumPartition,  UInt uiWidth, UInt uiHeight);
#else
  Void    create( UInt uiNumPartition );
#endif
  Void    destroy();

  // ------------------------------------------------------------------------------------------------------------------
  // clear / copy
  // ------------------------------------------------------------------------------------------------------------------

  Void    clearMvField();

#if JVET_C0024_QTBT
  Void    clearCtuMvField();
  Void    copyFromTo( TComCUMvField const * pcCUMvFieldSrc, Int iNumPart, Int iPartAddrSrc, Int iPartAddrDst );
#endif
  Void    copyFrom( TComCUMvField const * pcCUMvFieldSrc, Int iNumPartSrc, Int iPartAddrDst );
  Void    copyTo  ( TComCUMvField* pcCUMvFieldDst, Int iPartAddrDst ) const;
  Void    copyTo  ( TComCUMvField* pcCUMvFieldDst, Int iPartAddrDst, UInt uiOffset, UInt uiNumPart ) const;

  // ------------------------------------------------------------------------------------------------------------------
  // get
  // ------------------------------------------------------------------------------------------------------------------

  TComMv const & getMv    ( Int iIdx ) const { return  m_pcMv    [iIdx]; }
  TComMv const & getMvd   ( Int iIdx ) const { return  m_pcMvd   [iIdx]; }
  Int            getRefIdx( Int iIdx ) const { return  m_piRefIdx[iIdx]; }
#if VCEG_AZ07_FRUC_MERGE || COM16_C1016_AFFINE || JVET_C0024_QTBT
  Void           setMv    ( TComMv  cMv,     Int iIdx ) { m_pcMv    [iIdx] = cMv;     }
  Void           setRefIdx( Int     iRefIdx, Int iIdx ) { m_piRefIdx[iIdx] = iRefIdx; }
#endif

#if COM16_C1016_AFFINE
  Void            setMvd   ( TComMv cMvd, Int iIdx )     { m_pcMvd   [iIdx] = cMvd;    }
#if JVECT_C0062_AFFINE_SIX_PARAM
  AffineAMVPInfo* getAffineAMVPInfo(UInt idx) { return &m_cAffineAMVPInfo[idx]; }
#else
  AffineAMVPInfo* getAffineAMVPInfo () { return &m_cAffineAMVPInfo; }
#endif
#endif

  AMVPInfo* getAMVPInfo () { return &m_cAMVPInfo; }

  // ------------------------------------------------------------------------------------------------------------------
  // set
  // ------------------------------------------------------------------------------------------------------------------

#if JVET_C0024_QTBT
  Void    setCU(TComDataCU* pcCU) { m_pcCU = pcCU;}
#endif
  Void    setAllMv     ( TComMv const & rcMv,         PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx=0 );
  Void    setAllMvd    ( TComMv const & rcMvd,        PartSize eCUMode, Int iPartAddr, UInt uiDepth, Int iPartIdx=0 );
  Void    setAllRefIdx ( Int iRefIdx,                 PartSize eMbMode, Int iPartAddr, UInt uiDepth, Int iPartIdx=0 );
  Void    setAllMvField( TComMvField const & mvField, PartSize eMbMode, Int iPartAddr, UInt uiDepth, Int iPartIdx=0 );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP || VCEG_AZ07_FRUC_MERGE || JVET_C0024_QTBT
  Void    setMvFieldSP ( TComDataCU* pcCU, UInt uiAbsPartIdx, TComMvField cMvField, Int iWidth, Int iHeight  );
#endif

  Void setNumPartition( Int iNumPart )
  {
    m_uiNumPartition = iNumPart;
  }

  Void linkToWithOffset( TComCUMvField const * src, Int offset )
  {
    m_pcMv     = src->m_pcMv     + offset;
    m_pcMvd    = src->m_pcMvd    + offset;
    m_piRefIdx = src->m_piRefIdx + offset;
  }

#if MRG_CONSIDER_COMPLEX
  Void    setAllMvFieldNode( NodeInfo NodeMotionInfo,UInt MergeComplexIndex, PartSize eMbMode, RefPicList eRefPicList, Int iPartAddr, UInt uiDepth, Int iPartIdx=0 );
#endif

#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
  Void    setWidth( Int iWidth )      { m_uiWidth  = iWidth;  }
  Void    setHeight( Int iHeight )    { m_uiHeight = iHeight; }
#endif

  Void compress(Char* pePredMode, Int scale);
};

//! \}

#endif // __TCOMMOTIONINFO__
