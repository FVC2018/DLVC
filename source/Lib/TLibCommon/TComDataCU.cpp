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

/** \file     TComDataCU.cpp
    \brief    CU data structure
    \todo     not all entities are documented
*/

#include "TComDataCU.h"
#include "TComTU.h"
#include "TComPic.h"
#if VCEG_AZ07_FRUC_MERGE
#include "TComPrediction.h"
#endif



//! \ingroup TLibCommon
#if VCEG_AZ08_KLT_COMMON
extern UInt g_uiDepth2Width[5];
#endif
// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TComDataCU::TComDataCU()
{
  m_pcPic              = NULL;
  m_pcSlice            = NULL;
#if JVET_C0024_QTBT
  m_puhDepth[0] = m_puhDepth[1] = NULL;
#else
  m_puhDepth           = NULL;
#endif
#if TT_TREE_AFFINE //we have to store ZIdx of the top-left corner of each CU for affine
  m_puhCUabsZIdxInCtu[0][0] = m_puhCUabsZIdxInCtu[0][1] = NULL;
  m_puhCUabsZIdxInCtu[1][0] = m_puhCUabsZIdxInCtu[1][1] = NULL;
#endif

  m_skipFlag           = NULL;

#if VCEG_AZ05_INTRA_MPI
  m_MPIIdx             = NULL;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  m_PDPCIdx           = NULL;
#endif

#if RRC_DWDIP
  m_DWDIPIdx = NULL;
#endif
#if RRC_DWDIP_MLR
  m_DWDIP_GPIdx = NULL;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
  for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
  {
    m_ROTIdx[i] = NULL;
  }
#else
  m_ROTIdx           = NULL;
#endif
#endif
#if !JVET_C0024_QTBT
  m_pePartSize         = NULL;
#endif
  m_pePredMode         = NULL;
  m_CUTransquantBypass = NULL;
#if JVET_C0024_QTBT
  m_puhWidth[0] = m_puhWidth[1] = NULL;
  m_puhHeight[0] = m_puhHeight[1] = NULL;
#else
  m_puhWidth           = NULL;
  m_puhHeight          = NULL;
#endif
#if JVET_C0024_DELTA_QP_FIX
  for(UInt i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
  {
    m_phQP[i]          = NULL;
  }
#else
  m_phQP               = NULL;
#endif
  m_ChromaQpAdj        = NULL;
  m_pbMergeFlag        = NULL;
  m_puhMergeIndex      = NULL;
#if VCEG_AZ07_FRUC_MERGE
  m_puhFRUCMgrMode     = NULL;
#endif
  for(UInt i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
  {
    m_puhIntraDir[i]     = NULL;
#if MLIP_EXTEND_REF_LINE
    m_puhExternRef[i] = NULL;
#endif
  }
  m_puhInterDir        = NULL;
#if !JVET_C0024_QTBT
  m_puhTrIdx           = NULL;
#endif

  for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    m_puhCbf[comp]                        = NULL;
    m_crossComponentPredictionAlpha[comp] = NULL;
    m_puhTransformSkip[comp]              = NULL;
    m_pcTrCoeff[comp]                     = NULL;
#if ADAPTIVE_QP_SELECTION
    m_pcArlCoeff[comp]                    = NULL;
#endif
    m_pcIPCMSample[comp]                  = NULL;
    m_explicitRdpcmMode[comp]             = NULL;
#if VCEG_AZ08_KLT_COMMON
    m_puhKLTFlag[comp]                    = NULL;
#endif
  }
#if ADAPTIVE_QP_SELECTION
  m_ArlCoeffIsAliasedAllocation = false;
#endif
  m_pbIPCMFlag         = NULL;

  m_pCtuAboveLeft      = NULL;
  m_pCtuAboveRight     = NULL;
  m_pCtuAbove          = NULL;
  m_pCtuLeft           = NULL;

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i]  = NULL;
    m_apiMVPIdx[i]       = NULL;
    m_apiMVPNum[i]       = NULL;
  }

  m_bDecSubCu          = false;

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  m_peMergeType        = NULL;
#endif
#if VCEG_AZ07_IMV
  m_iMVFlag            = NULL;
  m_piMVCandNum        = NULL;
#endif
#if COM16_C806_OBMC
  m_OBMCFlag           = NULL;
#endif
#if VCEG_AZ06_IC
  m_pbICFlag           = NULL;
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
  m_puhAlfCtrlFlag     = NULL;
  m_puhTmpAlfCtrlFlag  = NULL;
#else
  m_puiAlfCtrlFlag     = NULL;
  m_puiTmpAlfCtrlFlag  = NULL;
#endif
#endif

#if DIMD_INTRA_PRED
  for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
  {
    m_ucDIMDEnabledFlag[ch] = NULL;
#if DIMD_NUM_INTRA_DIR_INC
    m_ucExtIntraDir[ch] = NULL;
#endif
    m_ucDIMDNoBTLevelFlag[ch] = NULL;
  }
  m_puhSecIntraDir = NULL;
#endif

#if COM16_C806_EMT
  m_puhEmtTuIdx        = NULL;
  m_puhEmtCuFlag       = NULL;
#endif

#if COM16_C1016_AFFINE
  m_affineFlag         = NULL;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  m_affineParamFlag    = NULL;
#endif

#if MRG_CONSIDER_COMPLEX
  m_pbMergeComplexFlag   = NULL;
  m_pNodeMotionInfo      = NULL;
  m_pNumMergeComplexCand = NULL;
  m_puhMergeComplexIndex = NULL;
#endif

#if RRC_INTRA_SHARPENING
  m_sharpFlag         = NULL;
#endif  
  
#if JVET_C0024_QTBT
  m_acCUMvField[0].setCU(this);
  m_acCUMvField[1].setCU(this);
#if TT_TREE_DS3
  m_bUseAbsoluteZIdxToCtu = false;
  m_pCtuCurr              = NULL;
  m_acCUMvFieldCtu[0].setCU(this);
  m_acCUMvFieldCtu[1].setCU(this);
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_acFRUCUniLateralMVField[0].setCU(this);
  m_acFRUCUniLateralMVField[1].setCU(this);
#endif
#endif

#if BLOCK_BASED_BACKGROUND_REFERENCE
	m_codedQP[CHANNEL_TYPE_LUMA] = 0;
	m_codedQP[CHANNEL_TYPE_CHROMA] = 0;
	m_LCUBgFlag = false;
#endif

#if RRC_SELECTABLE_IF
  m_pSIFFlag      = NULL;
  m_sifToCheck      = 0;
#endif 
}

TComDataCU::~TComDataCU()
{
}

Void TComDataCU::create( ChromaFormat chromaFormatIDC, UInt uiNumPartition, UInt uiWidth, UInt uiHeight, Bool bDecSubCu, Int unitSize
#if JVET_C0024_QTBT
                        , UInt uiCUWidth, UInt uiCUHeight
#endif
#if ADAPTIVE_QP_SELECTION
                        , TCoeff *pParentARLBuffer
#endif
                        )
{
  m_bDecSubCu = bDecSubCu;

  m_pcPic              = NULL;
  m_pcSlice            = NULL;
  m_uiNumPartition     = uiNumPartition;
  m_unitSize = unitSize;
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  m_uiSubCUIdxCompST = 0;
#endif
#if BT_FORCE_FBTD
  m_uiFBTDepth = 0;
#endif

#if DIMD_INTRA_PRED
  m_ctuRsAddr = MAX_UINT;
#endif

  if ( !bDecSubCu )
  {
#if !JVET_C0024_DELTA_QP_FIX
    m_phQP               = (Char*     )xMalloc(Char,     uiNumPartition);
#endif
#if JVET_C0024_QTBT
    for (Int i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
    {
#if JVET_C0024_DELTA_QP_FIX
      m_phQP[i]               = (Char*     )xMalloc( Char,    uiNumPartition);
#endif
      m_puhDepth[i]           = (UChar*    )xMalloc(UChar,    uiNumPartition);
      m_puhWidth[i]           = (UChar*    )xMalloc(UChar,    uiNumPartition);
      m_puhHeight[i]          = (UChar*    )xMalloc(UChar,    uiNumPartition);
#if TT_TREE_AFFINE
      m_puhCUabsZIdxInCtu[i][0] = (UChar*    )xMalloc(UChar,    uiNumPartition);
      m_puhCUabsZIdxInCtu[i][1] = (UChar*    )xMalloc(UChar,    uiNumPartition);
#endif

      m_puhBTSplitMode[i][0]     = new UChar[ uiNumPartition ];
      memset( m_puhBTSplitMode[i][0], 0,uiNumPartition * sizeof( *m_puhBTSplitMode[i][0] ) );
      m_puhBTSplitMode[i][1]     = new UChar[ uiNumPartition ];
      memset( m_puhBTSplitMode[i][1], 0,uiNumPartition * sizeof( *m_puhBTSplitMode[i][1] ) );
#if BT_FORCE_BF
      m_puhBTSplitMode[i][2] = new UChar[uiNumPartition];
      memset(m_puhBTSplitMode[i][2], 0, uiNumPartition * sizeof(*m_puhBTSplitMode[i][2]));
      m_puhBTSplitMode[i][3] = new UChar[uiNumPartition];
      memset(m_puhBTSplitMode[i][3], 0, uiNumPartition * sizeof(*m_puhBTSplitMode[i][3]));
#endif
    }
#else
    m_puhDepth           = (UChar*    )xMalloc(UChar,    uiNumPartition);
#if COM16_C806_LARGE_CTU
    m_puhWidth           = (UShort*   )xMalloc(UShort,   uiNumPartition);
    m_puhHeight          = (UShort*   )xMalloc(UShort,   uiNumPartition);
#else
    m_puhWidth           = (UChar*    )xMalloc(UChar,    uiNumPartition);
    m_puhHeight          = (UChar*    )xMalloc(UChar,    uiNumPartition);
#endif
#endif

    m_ChromaQpAdj        = new UChar[ uiNumPartition ];
    m_skipFlag           = new Bool[ uiNumPartition ];
#if VCEG_AZ05_INTRA_MPI
    m_MPIIdx             = new Char[ uiNumPartition ];
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    m_PDPCIdx            = new Char[ uiNumPartition];
#endif
#if RRC_DWDIP
    m_DWDIPIdx           = new Char[ uiNumPartition];
#endif

#if RRC_DWDIP_MLR
    m_DWDIP_GPIdx            = new Char[ uiNumPartition];
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
    {
      m_ROTIdx[i] = (Char* )xMalloc(Char,  uiNumPartition);
    }
#else
    m_ROTIdx           = new Char[ uiNumPartition ];
#endif
#endif
#if !JVET_C0024_QTBT
    m_pePartSize         = new Char[ uiNumPartition ];
    memset( m_pePartSize, NUMBER_OF_PART_SIZES,uiNumPartition * sizeof( *m_pePartSize ) );
#endif

    m_pePredMode         = new Char[ uiNumPartition ];
    m_CUTransquantBypass = new Bool[ uiNumPartition ];

    m_pbMergeFlag        = (Bool*  )xMalloc(Bool,   uiNumPartition);
    m_puhMergeIndex      = (UChar* )xMalloc(UChar,  uiNumPartition);
#if VCEG_AZ07_FRUC_MERGE
    m_puhFRUCMgrMode     = (UChar* )xMalloc(UChar,  uiNumPartition);
#endif

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    m_peMergeType        = (UChar*  )xMalloc(MergeType,   uiNumPartition);
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
    m_iMVFlag            = new UChar[ uiNumPartition ];
#else
    m_iMVFlag            = new Bool[ uiNumPartition ];
#endif
    m_piMVCandNum        = new Char[ uiNumPartition ];
#endif
#if COM16_C806_OBMC
    m_OBMCFlag           = new Bool[ uiNumPartition ];
#endif
#if VCEG_AZ06_IC
    m_pbICFlag           = new Bool[ uiNumPartition ];
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    m_puhAlfCtrlFlag     = (UChar*  )xMalloc(UChar,   uiNumPartition);
#else
    m_puiAlfCtrlFlag     = (UInt*  )xMalloc(UInt,   uiNumPartition);
#endif
#endif

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
      m_ucDIMDEnabledFlag[ch] = (UChar*)xMalloc(UChar, uiNumPartition);
#if DIMD_NUM_INTRA_DIR_INC
      m_ucExtIntraDir[ch] = (UChar*)xMalloc(UChar, uiNumPartition);
#endif
      m_ucDIMDNoBTLevelFlag[ch] = (UChar*)xMalloc(UChar, uiNumPartition);
    }
    m_puhSecIntraDir = (UChar*)xMalloc(UChar, uiNumPartition);
#endif

#if COM16_C806_EMT
    m_puhEmtTuIdx        = (UChar*)xMalloc(UChar, uiNumPartition);
    m_puhEmtCuFlag       = (UChar*)xMalloc(UChar, uiNumPartition);
#endif

#if COM16_C1016_AFFINE
    m_affineFlag         = (Bool*  )xMalloc(Bool,   uiNumPartition);
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    m_affineParamFlag    = (Bool*)xMalloc(Bool, uiNumPartition);
#endif

#if MRG_CONSIDER_COMPLEX
    m_pbMergeComplexFlag   = (Bool*  )   xMalloc(Bool,     uiNumPartition);
    m_pNodeMotionInfo      = (NodeInfo* )xMalloc(NodeInfo, uiNumPartition);
    m_pNumMergeComplexCand = (UChar* )   xMalloc(UChar,    uiNumPartition);
    m_puhMergeComplexIndex = (UChar* )   xMalloc(UChar,    uiNumPartition);
#endif

#if RRC_INTRA_SHARPENING
    m_sharpFlag         = NULL;
    m_sharpFlag         = new Char[uiNumPartition];
    memset(m_sharpFlag, 0, sizeof(*m_sharpFlag)*uiNumPartition); // TODO: check if really need this
#endif

#if RRC_SELECTABLE_IF
    m_pSIFFlag = (Char*)xMalloc(Char, uiNumPartition); //new Char[uiNumPartition];
    memset(m_pSIFFlag, 0, sizeof(Char)*uiNumPartition);
#endif
    
    for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
    {
      m_puhIntraDir[ch] = (UChar* )xMalloc(UChar,  uiNumPartition);
#if MLIP_EXTEND_REF_LINE
      m_puhExternRef[ch] = (UChar*)xMalloc(UChar, uiNumPartition);
#endif
    }
    m_puhInterDir        = (UChar* )xMalloc(UChar,  uiNumPartition);

#if !JVET_C0024_QTBT
    m_puhTrIdx           = (UChar* )xMalloc(UChar,  uiNumPartition);
#endif

    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      const RefPicList rpl=RefPicList(i);
      m_apiMVPIdx[rpl]       = new Char[ uiNumPartition ];
      m_apiMVPNum[rpl]       = new Char[ uiNumPartition ];
      memset( m_apiMVPIdx[rpl], -1,uiNumPartition * sizeof( Char ) );
    }

    for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      const ComponentID compID = ComponentID(comp);
      const UInt chromaShift = getComponentScaleX(compID, chromaFormatIDC) + getComponentScaleY(compID, chromaFormatIDC);
      const UInt totalSize   = (uiWidth * uiHeight) >> chromaShift;

      m_crossComponentPredictionAlpha[compID] = (Char*  )xMalloc(Char,   uiNumPartition);
      m_puhTransformSkip[compID]              = (UChar* )xMalloc(UChar,  uiNumPartition);
      m_explicitRdpcmMode[compID]             = (UChar* )xMalloc(UChar,  uiNumPartition);
      m_puhCbf[compID]                        = (UChar* )xMalloc(UChar,  uiNumPartition);
      m_pcTrCoeff[compID]                     = (TCoeff*)xMalloc(TCoeff, totalSize);
      memset( m_pcTrCoeff[compID], 0, (totalSize * sizeof( TCoeff )) );
#if VCEG_AZ08_KLT_COMMON
      m_puhKLTFlag[compID]                    = (UChar*)xMalloc(UChar,  uiNumPartition);
#endif
#if ADAPTIVE_QP_SELECTION
      if( pParentARLBuffer != 0 )
      {
        m_pcArlCoeff[compID] = pParentARLBuffer;
        m_ArlCoeffIsAliasedAllocation = true;
        pParentARLBuffer += totalSize;
      }
      else
      {
        m_pcArlCoeff[compID] = (TCoeff*)xMalloc(TCoeff, totalSize);
        m_ArlCoeffIsAliasedAllocation = false;
      }
#endif
      m_pcIPCMSample[compID] = (Pel*   )xMalloc(Pel , totalSize);
#if JVET_C0024_ITSKIP
      m_puiSkipWidth[compID] = (UInt*)xMalloc(UInt,    uiNumPartition);
      memset( m_puiSkipWidth[compID], 0, uiNumPartition * sizeof( UInt ) );
      m_puiSkipHeight[compID] = (UInt*)xMalloc(UInt,    uiNumPartition);
      memset( m_puiSkipHeight[compID], 0, uiNumPartition * sizeof( UInt ) );
#endif
    }

    m_pbIPCMFlag         = (Bool*  )xMalloc(Bool, uiNumPartition);

    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
      m_acCUMvField[i].create(uiNumPartition, uiWidth, uiHeight);
#else
      m_acCUMvField[i].create( uiNumPartition );
#endif
#if JVET_C0024_QTBT
      UInt uiCUNumPart = uiCUWidth*uiCUHeight /(unitSize*unitSize);
      //create big buffer since the BT not follow zorder,but set the buffer size as CU size
      m_acCUMvField[i].setNumPartition( uiCUNumPart);
#endif

#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
      m_acCUMvField[i].setWidth(uiWidth);
      m_acCUMvField[i].setHeight(uiHeight);
#endif
    }

  }
  else
  {
    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
#if JVET_C0024_QTBT
      UInt uiCUNumPart = uiCUWidth*uiCUHeight /(unitSize*unitSize);
      //set the buffer size as CU size
      m_acCUMvField[i].setNumPartition( uiCUNumPart);
#else
      m_acCUMvField[i].setNumPartition(uiNumPartition );
#endif
    }
  }

#if VCEG_AZ07_FRUC_MERGE
#if MRG_CONSIDER_COMPLEX && !JVET_C0024_QTBT
  m_acFRUCUniLateralMVField[0].create(uiNumPartition, uiWidth, uiHeight);
  m_acFRUCUniLateralMVField[1].create(uiNumPartition, uiWidth, uiHeight);
#else
  m_acFRUCUniLateralMVField[0].create( uiNumPartition );
  m_acFRUCUniLateralMVField[1].create( uiNumPartition );
#endif
#if JVET_C0024_QTBT
  UInt uiCUNumPart = uiCUWidth*uiCUHeight /(unitSize*unitSize);
  //create big buffer since the BT not follow zorder,but set the buffer size as CU size
  m_acFRUCUniLateralMVField[0].setNumPartition( uiCUNumPart);
  m_acFRUCUniLateralMVField[1].setNumPartition( uiCUNumPart);
#endif
#endif

  // create motion vector fields

  m_pCtuAboveLeft      = NULL;
  m_pCtuAboveRight     = NULL;
  m_pCtuAbove          = NULL;
  m_pCtuLeft           = NULL;

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i]  = NULL;
  }
#if JVET_C0024_QTBT
  m_uiNumPartition = uiCUWidth*uiCUHeight /(unitSize*unitSize);
#endif
}

Void TComDataCU::destroy()
{
  // encoder-side buffer free
  if ( !m_bDecSubCu )
  {
#if !JVET_C0024_DELTA_QP_FIX
    if ( m_phQP )
    {
      xFree(m_phQP);
      m_phQP = NULL;
    }
#endif
#if JVET_C0024_QTBT
    for (Int i=0; i<MAX_NUM_CHANNEL_TYPE; i++)
    {
#if JVET_C0024_DELTA_QP_FIX
      if ( m_phQP[i]               ) { xFree(m_phQP[i]);                m_phQP[i]              = NULL; }
#endif
      if ( m_puhDepth[i]           ) { xFree(m_puhDepth[i]);            m_puhDepth[i]          = NULL; }
      if ( m_puhWidth[i]           ) { xFree(m_puhWidth[i]);            m_puhWidth[i]          = NULL; }
      if ( m_puhHeight[i]          ) { xFree(m_puhHeight[i]);           m_puhHeight[i]         = NULL; }
#if TT_TREE_AFFINE
      if ( m_puhCUabsZIdxInCtu[i][0]  ) { xFree(m_puhCUabsZIdxInCtu[i][0]);   m_puhCUabsZIdxInCtu[i][0] = NULL; }
      if ( m_puhCUabsZIdxInCtu[i][1]  ) { xFree(m_puhCUabsZIdxInCtu[i][1]);   m_puhCUabsZIdxInCtu[i][1] = NULL; }
#endif

      if ( m_puhBTSplitMode[i][0]     ) { delete[] m_puhBTSplitMode[i][0];    m_puhBTSplitMode[i][0]    = NULL; }
      if ( m_puhBTSplitMode[i][1]     ) { delete[] m_puhBTSplitMode[i][1];    m_puhBTSplitMode[i][1]    = NULL; }
#if BT_FORCE_BF
      if (m_puhBTSplitMode[i][2]) { delete[] m_puhBTSplitMode[i][2];    m_puhBTSplitMode[i][2] = NULL; }
      if (m_puhBTSplitMode[i][3]) { delete[] m_puhBTSplitMode[i][3];    m_puhBTSplitMode[i][3] = NULL; }
#endif
    }
#else
    if ( m_puhDepth )
    {
      xFree(m_puhDepth);
      m_puhDepth = NULL;
    }
    if ( m_puhWidth )
    {
      xFree(m_puhWidth);
      m_puhWidth = NULL;
    }
    if ( m_puhHeight )
    {
      xFree(m_puhHeight);
      m_puhHeight = NULL;
    }
#endif

    if ( m_skipFlag )
    {
      delete[] m_skipFlag;
      m_skipFlag = NULL;
    }
#if VCEG_AZ05_INTRA_MPI
    if ( m_MPIIdx ) 
    { 
      delete[] m_MPIIdx;          
      m_MPIIdx = NULL; 
    }
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    if (m_PDPCIdx)
    {
      delete[] m_PDPCIdx;
      m_PDPCIdx = NULL;
    }
#endif
#if RRC_DWDIP
    if (m_DWDIPIdx)
    {
      delete[] m_DWDIPIdx;
      m_DWDIPIdx = NULL;
    }
#endif
#if RRC_DWDIP_MLR
    if (m_DWDIP_GPIdx)
    {
      delete[] m_DWDIP_GPIdx;
      m_DWDIP_GPIdx = NULL;
    }
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
    {
      xFree( m_ROTIdx[i] );          
      m_ROTIdx[i]     = NULL; 
    }
#else
    if ( m_ROTIdx           ) 
    { 
      delete[] m_ROTIdx;          
      m_ROTIdx          = NULL; 
    }
#endif
#endif
#if !JVET_C0024_QTBT
    if ( m_pePartSize )
    {
      delete[] m_pePartSize;
      m_pePartSize = NULL;
    }
#endif
    if ( m_pePredMode )
    {
      delete[] m_pePredMode;
      m_pePredMode = NULL;
    }
    if ( m_ChromaQpAdj )
    {
      delete[] m_ChromaQpAdj;
      m_ChromaQpAdj = NULL;
    }
    if ( m_CUTransquantBypass )
    {
      delete[] m_CUTransquantBypass;
      m_CUTransquantBypass = NULL;
    }
    if ( m_puhInterDir )
    {
      xFree(m_puhInterDir);
      m_puhInterDir = NULL;
    }
    if ( m_pbMergeFlag )
    {
      xFree(m_pbMergeFlag);
      m_pbMergeFlag = NULL;
    }
    if ( m_puhMergeIndex )
    {
      xFree(m_puhMergeIndex);
      m_puhMergeIndex  = NULL;
    }
#if VCEG_AZ07_FRUC_MERGE
    if ( m_puhFRUCMgrMode ) 
    { 
      xFree(m_puhFRUCMgrMode);      
      m_puhFRUCMgrMode = NULL; 
    }
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    if ( m_peMergeType ) 
    {
      xFree(m_peMergeType ); 
      m_peMergeType      = NULL;  
    }
#endif
#if VCEG_AZ07_IMV
    if ( m_iMVFlag ) 
    { 
      delete[] m_iMVFlag;           
      m_iMVFlag = NULL; 
    }
    if ( m_piMVCandNum ) 
    { 
      delete[] m_piMVCandNum;
      m_piMVCandNum = NULL; 
    }
#endif
#if COM16_C806_OBMC
    if ( m_OBMCFlag ) 
    { 
      delete[] m_OBMCFlag;
      m_OBMCFlag = NULL; 
    }
#endif
#if VCEG_AZ06_IC
    if ( m_pbICFlag )
    { 
      delete [] m_pbICFlag;
      m_pbICFlag = NULL; 
    }
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    if ( m_puhAlfCtrlFlag ) 
    { 
      xFree(m_puhAlfCtrlFlag);      
      m_puhAlfCtrlFlag = NULL; 
    }
#else
    if ( m_puiAlfCtrlFlag ) 
    { 
      xFree(m_puiAlfCtrlFlag);      
      m_puiAlfCtrlFlag = NULL; 
    }
#endif
#endif

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
      xFree(m_ucDIMDEnabledFlag[ch]);
      m_ucDIMDEnabledFlag[ch] = NULL;
#if DIMD_NUM_INTRA_DIR_INC
      xFree(m_ucExtIntraDir[ch]);
      m_ucExtIntraDir[ch] = NULL;
#endif
      xFree(m_ucDIMDNoBTLevelFlag[ch]);
      m_ucDIMDNoBTLevelFlag[ch] = NULL;
    }
    xFree(m_puhSecIntraDir);
    m_puhSecIntraDir = NULL;
#endif

#if COM16_C806_EMT
    if ( m_puhEmtTuIdx ) 
    { 
      xFree(m_puhEmtTuIdx);
      m_puhEmtTuIdx       = NULL; 
    }

    if ( m_puhEmtCuFlag )
    { 
      xFree(m_puhEmtCuFlag);
      m_puhEmtCuFlag      = NULL; 
    }
#endif

#if COM16_C1016_AFFINE
    if ( m_affineFlag )
    {
      xFree(m_affineFlag);
      m_affineFlag        = NULL;
    }
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    if (m_affineParamFlag)
    {
      xFree(m_affineParamFlag);
      m_affineParamFlag = NULL;
    }
#endif

#if MRG_CONSIDER_COMPLEX
    if ( m_pbMergeComplexFlag )       { xFree(m_pbMergeComplexFlag);    m_pbMergeComplexFlag   = NULL; }
    if ( m_pNodeMotionInfo    )       { xFree(m_pNodeMotionInfo);       m_pNodeMotionInfo      = NULL; }
    if ( m_pNumMergeComplexCand )     { xFree(m_pNumMergeComplexCand);  m_pNumMergeComplexCand = NULL; }
    if ( m_puhMergeComplexIndex )     { xFree(m_puhMergeComplexIndex);  m_puhMergeComplexIndex = NULL; }
#endif

#if RRC_INTRA_SHARPENING
    if (m_sharpFlag)
    {
      delete[] m_sharpFlag;
      m_sharpFlag         = NULL; 
    }
#endif

#if RRC_SELECTABLE_IF
    if ( m_pSIFFlag )
    { 
      //delete[] m_pSIFFlag;
      xFree(m_pSIFFlag);
      m_pSIFFlag      = NULL; 
    }
#endif
#if MLIP_EXTEND_REF_LINE
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
        xFree(m_puhExternRef[ch]);
        m_puhExternRef[ch] = NULL;
    }
#endif
    for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
    {
      xFree(m_puhIntraDir[ch]);
      m_puhIntraDir[ch] = NULL;
    }

#if !JVET_C0024_QTBT
    if ( m_puhTrIdx )
    {
      xFree(m_puhTrIdx);
      m_puhTrIdx = NULL;
    }
#endif

    for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      if ( m_crossComponentPredictionAlpha[comp] )
      {
        xFree(m_crossComponentPredictionAlpha[comp]);
        m_crossComponentPredictionAlpha[comp] = NULL;
      }
      if ( m_puhTransformSkip[comp] )
      {
        xFree(m_puhTransformSkip[comp]);
        m_puhTransformSkip[comp] = NULL;
      }
      if ( m_puhCbf[comp] )
      {
        xFree(m_puhCbf[comp]);
        m_puhCbf[comp] = NULL;
      }
      if ( m_pcTrCoeff[comp] )
      {
        xFree(m_pcTrCoeff[comp]);
        m_pcTrCoeff[comp] = NULL;
      }
      if ( m_explicitRdpcmMode[comp] )
      {
        xFree(m_explicitRdpcmMode[comp]);
        m_explicitRdpcmMode[comp] = NULL;
      }
#if VCEG_AZ08_KLT_COMMON
      if (m_puhKLTFlag[comp])
      { 
        xFree(m_puhKLTFlag[comp]);
        m_puhKLTFlag[comp] = NULL;
      }
#endif
#if ADAPTIVE_QP_SELECTION
      if (!m_ArlCoeffIsAliasedAllocation)
      {
        if ( m_pcArlCoeff[comp] )
        {
          xFree(m_pcArlCoeff[comp]);
          m_pcArlCoeff[comp] = NULL;
        }
      }
#endif

      if ( m_pcIPCMSample[comp] )
      {
        xFree(m_pcIPCMSample[comp]);
        m_pcIPCMSample[comp] = NULL;
      }
#if JVET_C0024_ITSKIP
      if ( m_puiSkipWidth[comp] ) {xFree(m_puiSkipWidth[comp]); m_puiSkipWidth[comp] = NULL;}
      if ( m_puiSkipHeight[comp] ) {xFree(m_puiSkipHeight[comp]); m_puiSkipHeight[comp] = NULL;}
#endif
    }
    if ( m_pbIPCMFlag )
    {
      xFree(m_pbIPCMFlag );
      m_pbIPCMFlag = NULL;
    }

    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      const RefPicList rpl=RefPicList(i);
      if ( m_apiMVPIdx[rpl] )
      {
        delete[] m_apiMVPIdx[rpl];
        m_apiMVPIdx[rpl] = NULL;
      }
      if ( m_apiMVPNum[rpl] )
      {
        delete[] m_apiMVPNum[rpl];
        m_apiMVPNum[rpl] = NULL;
      }
    }

    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      const RefPicList rpl=RefPicList(i);
      m_acCUMvField[rpl].destroy();
    }
  }

#if VCEG_AZ07_FRUC_MERGE
  m_acFRUCUniLateralMVField[0].destroy();
  m_acFRUCUniLateralMVField[1].destroy();
#endif

  m_pcPic              = NULL;
  m_pcSlice            = NULL;

  m_pCtuAboveLeft      = NULL;
  m_pCtuAboveRight     = NULL;
  m_pCtuAbove          = NULL;
  m_pCtuLeft           = NULL;


  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i]  = NULL;
  }

}

Bool TComDataCU::CUIsFromSameTile            ( const TComDataCU *pCU /* Can be NULL */) const
{
  return pCU!=NULL &&
         pCU->getSlice() != NULL &&
         m_pcPic->getPicSym()->getTileIdxMap( pCU->getCtuRsAddr() ) == m_pcPic->getPicSym()->getTileIdxMap(getCtuRsAddr());
}

Bool TComDataCU::CUIsFromSameSliceAndTile    ( const TComDataCU *pCU /* Can be NULL */) const
{
  return pCU!=NULL &&
         pCU->getSlice() != NULL &&
         pCU->getSlice()->getSliceCurStartCtuTsAddr() == getSlice()->getSliceCurStartCtuTsAddr() &&
         m_pcPic->getPicSym()->getTileIdxMap( pCU->getCtuRsAddr() ) == m_pcPic->getPicSym()->getTileIdxMap(getCtuRsAddr())
         ;
}

Bool TComDataCU::CUIsFromSameSliceTileAndWavefrontRow( const TComDataCU *pCU /* Can be NULL */) const
{
  return CUIsFromSameSliceAndTile(pCU)
         && (!getSlice()->getPPS()->getEntropyCodingSyncEnabledFlag() || getPic()->getCtu(getCtuRsAddr())->getCUPelY() == getPic()->getCtu(pCU->getCtuRsAddr())->getCUPelY());
}

Bool TComDataCU::isLastSubCUOfCtu(const UInt absPartIdx)
{
  const TComSPS &sps=*(getSlice()->getSPS());

  const UInt picWidth = sps.getPicWidthInLumaSamples();
  const UInt picHeight = sps.getPicHeightInLumaSamples();
#if JVET_C0024_QTBT
  const UInt granularityWidth = sps.getCTUSize();
#else
  const UInt granularityWidth = sps.getMaxCUWidth();
#endif

  const UInt cuPosX = getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[absPartIdx] ];
  const UInt cuPosY = getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[absPartIdx] ];

  return (((cuPosX+getWidth( absPartIdx))%granularityWidth==0||(cuPosX+getWidth( absPartIdx)==picWidth ))
       && ((cuPosY+getHeight(absPartIdx))%granularityWidth==0||(cuPosY+getHeight(absPartIdx)==picHeight)));
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

// --------------------------------------------------------------------------------------------------------------------
// Initialization
// --------------------------------------------------------------------------------------------------------------------

/**
 Initialize top-level CU: create internal buffers and set initial values before encoding the CTU.
 
 \param  pcPic       picture (TComPic) class pointer
 \param  ctuRsAddr   CTU address in raster scan order
 */
Void TComDataCU::initCtu( TComPic* pcPic, UInt ctuRsAddr )
{
#if JVET_C0024_QTBT
  const UInt maxCUWidth = pcPic->getPicSym()->getSPS().getCTUSize();
  const UInt maxCUHeight= pcPic->getPicSym()->getSPS().getCTUSize();
#else
  const UInt maxCUWidth = pcPic->getPicSym()->getSPS().getMaxCUWidth();
  const UInt maxCUHeight= pcPic->getPicSym()->getSPS().getMaxCUHeight();
#endif
  m_pcPic              = pcPic;
  m_pcSlice            = pcPic->getSlice(pcPic->getCurrSliceIdx());
  m_ctuRsAddr          = ctuRsAddr;
  m_uiCUPelX           = ( ctuRsAddr % pcPic->getFrameWidthInCtus() ) * maxCUWidth;
  m_uiCUPelY           = ( ctuRsAddr / pcPic->getFrameWidthInCtus() ) * maxCUHeight;
  m_absZIdxInCtu       = 0;
  m_dTotalCost         = MAX_DOUBLE;
  m_uiTotalDistortion  = 0;
  m_uiTotalBits        = 0;
  m_uiTotalBins        = 0;
  m_uiNumPartition     = pcPic->getNumPartitionsInCtu();

#if BLOCK_BASED_BACKGROUND_REFERENCE
  m_codedQP[CHANNEL_TYPE_LUMA] = 0;
  m_codedQP[CHANNEL_TYPE_CHROMA] = 0;
  m_LCUBgFlag          = false;
#endif

#if BARC
  m_BARCFlag[0] = m_BARCFlag[1] = false;
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = false;
#endif
  memset( m_skipFlag          , false,                      m_uiNumPartition * sizeof( *m_skipFlag ) );
#if VCEG_AZ05_INTRA_MPI
  memset( m_MPIIdx            , 0,                      m_uiNumPartition * sizeof( *m_MPIIdx ) );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  memset( m_PDPCIdx           , 0,                      m_uiNumPartition * sizeof(*m_PDPCIdx) );
#endif

#if RRC_DWDIP
  memset( m_DWDIPIdx          , 0,                      m_uiNumPartition * sizeof(*m_DWDIPIdx));
#endif
#if RRC_DWDIP_MLR
  memset( m_DWDIP_GPIdx           , 0,                      m_uiNumPartition * sizeof(*m_DWDIP_GPIdx));
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
  for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ )
  {
    memset( m_ROTIdx[i]      , 0,                      m_uiNumPartition * sizeof( *m_ROTIdx[i] ) );
  }
#else
  memset( m_ROTIdx            , 0,                      m_uiNumPartition * sizeof( *m_ROTIdx ) );
#endif
#endif 

#if !JVET_C0024_QTBT
  memset( m_pePartSize        , NUMBER_OF_PART_SIZES,       m_uiNumPartition * sizeof( *m_pePartSize ) );
#endif
  memset( m_pePredMode        , NUMBER_OF_PREDICTION_MODES, m_uiNumPartition * sizeof( *m_pePredMode ) );
  memset( m_CUTransquantBypass, false,                      m_uiNumPartition * sizeof( *m_CUTransquantBypass) );
#if JVET_C0024_QTBT
  for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
  {
#if !BLOCK_BASED_BACKGROUND_REFERENCE	
#if JVET_C0024_DELTA_QP_FIX
    memset( m_phQP[ch]              , getSlice()->getSliceQp(),   m_uiNumPartition * sizeof( *m_phQP[ch] ) );
#endif
#else
		if (getSlice()->getSliceType() != I_SLICE && getSlice()->getQPList() != NULL)
		{
			memset(m_phQP[ch], getSlice()->getQPList(ctuRsAddr), m_uiNumPartition * sizeof(*m_phQP[ch]));
		}
		else
		{
			memset(m_phQP[ch], getSlice()->getSliceQp(), m_uiNumPartition * sizeof(*m_phQP[ch]));
		}
#endif
    memset( m_puhDepth[ch]          , 0,                          m_uiNumPartition * sizeof( *m_puhDepth[ch] ) );
#if JVET_C0024_CTU_256
    memset( m_puhWidth[ch]          , maxCUWidth >> MIN_CU_LOG2,  m_uiNumPartition * sizeof( *m_puhWidth[ch] ) );
    memset( m_puhHeight[ch]         , maxCUHeight>> MIN_CU_LOG2,  m_uiNumPartition * sizeof( *m_puhHeight[ch] ) );
#else
    memset( m_puhWidth[ch]          , maxCUWidth,                 m_uiNumPartition * sizeof( *m_puhWidth[ch] ) );
    memset( m_puhHeight[ch]         , maxCUHeight,                m_uiNumPartition * sizeof( *m_puhHeight[ch] ) );
#endif
#if TT_TREE_AFFINE
    memset( m_puhCUabsZIdxInCtu[ch][0] , 0,                        m_uiNumPartition * sizeof( *m_puhCUabsZIdxInCtu[ch][0]));
    memset( m_puhCUabsZIdxInCtu[ch][1] , 0,                        m_uiNumPartition * sizeof( *m_puhCUabsZIdxInCtu[ch][1]));
#endif
    memset( m_puhBTSplitMode[ch][0]    , 0,                        m_uiNumPartition * sizeof( *m_puhBTSplitMode[ch][0]));
    memset( m_puhBTSplitMode[ch][1]    , 0,                        m_uiNumPartition * sizeof( *m_puhBTSplitMode[ch][1]));
#if BT_FORCE_BF
    memset(m_puhBTSplitMode[ch][2], 0, m_uiNumPartition * sizeof(*m_puhBTSplitMode[ch][2]));
    memset(m_puhBTSplitMode[ch][3], 0, m_uiNumPartition * sizeof(*m_puhBTSplitMode[ch][3]));
#endif
  }
#else
  memset( m_puhDepth          , 0,                          m_uiNumPartition * sizeof( *m_puhDepth ) );
#endif
#if !JVET_C0024_QTBT
  memset( m_puhTrIdx          , 0,                          m_uiNumPartition * sizeof( *m_puhTrIdx ) );
#if COM16_C806_LARGE_CTU
  for( Int n = 0 ; n < m_uiNumPartition ; n++ )
  {
    m_puhWidth[n]  = maxCUWidth;
    m_puhHeight[n] = maxCUHeight;
  }
#else
  memset( m_puhWidth          , maxCUWidth,                 m_uiNumPartition * sizeof( *m_puhWidth ) );
  memset( m_puhHeight         , maxCUHeight,                m_uiNumPartition * sizeof( *m_puhHeight ) );
#endif
#endif
  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    memset( m_apiMVPIdx[rpl]  , -1,                         m_uiNumPartition * sizeof( *m_apiMVPIdx[rpl] ) );
    memset( m_apiMVPNum[rpl]  , -1,                         m_uiNumPartition * sizeof( *m_apiMVPNum[rpl] ) );
  }
#if !JVET_C0024_DELTA_QP_FIX
  memset( m_phQP              , getSlice()->getSliceQp(),   m_uiNumPartition * sizeof( *m_phQP ) );
#endif
  memset( m_ChromaQpAdj       , 0,                          m_uiNumPartition * sizeof( *m_ChromaQpAdj ) );
  for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    memset( m_crossComponentPredictionAlpha[comp] , 0,                     m_uiNumPartition * sizeof( *m_crossComponentPredictionAlpha[comp] ) );
    memset( m_puhTransformSkip[comp]              , 0,                     m_uiNumPartition * sizeof( *m_puhTransformSkip[comp]) );
    memset( m_puhCbf[comp]                        , 0,                     m_uiNumPartition * sizeof( *m_puhCbf[comp] ) );
    memset( m_explicitRdpcmMode[comp]             , NUMBER_OF_RDPCM_MODES, m_uiNumPartition * sizeof( *m_explicitRdpcmMode[comp] ) );
#if VCEG_AZ08_KLT_COMMON
    memset( m_puhKLTFlag[comp]                    , 0,                     m_uiNumPartition * sizeof( *m_puhKLTFlag[comp] ) );
#endif
  }
  memset( m_pbMergeFlag       , false,                    m_uiNumPartition * sizeof( *m_pbMergeFlag ) );
  memset( m_puhMergeIndex     , 0,                        m_uiNumPartition * sizeof( *m_puhMergeIndex ) );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  memset( m_peMergeType       , MGR_TYPE_DEFAULT_N,       m_uiNumPartition * sizeof( *m_peMergeType   ) );
#endif
#if COM16_C806_OBMC
  memset( m_OBMCFlag          , true,                     m_uiNumPartition * sizeof( *m_OBMCFlag ) );
#endif
#if VCEG_AZ07_FRUC_MERGE
  memset( m_puhFRUCMgrMode    , 0,                        m_uiNumPartition * sizeof( *m_puhFRUCMgrMode ) );
#endif
#if VCEG_AZ07_IMV
  memset( m_iMVFlag           , false,                    m_uiNumPartition * sizeof( *m_iMVFlag ) );
  memset( m_piMVCandNum       , 0,                        m_uiNumPartition * sizeof( *m_piMVCandNum ) );
#endif
#if VCEG_AZ06_IC
  memset( m_pbICFlag          , false,                    m_uiNumPartition * sizeof( *m_pbICFlag ) );
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
  memset( m_puhAlfCtrlFlag    , 0,                        m_uiNumPartition * sizeof( *m_puhAlfCtrlFlag ) );
#else
  memset( m_puiAlfCtrlFlag    , 0,                        m_uiNumPartition * sizeof( *m_puiAlfCtrlFlag ) );
#endif
#endif

#if DIMD_INTRA_PRED
  for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
  {
    memset(m_ucDIMDEnabledFlag[ch], 0, m_uiNumPartition * sizeof(*m_ucDIMDEnabledFlag[ch]));
#if DIMD_NUM_INTRA_DIR_INC
    memset(m_ucExtIntraDir[ch], ((ch == 0) ? DC_IDX : 0), m_uiNumPartition * sizeof(*(m_ucExtIntraDir[ch])));
#endif
    memset(m_ucDIMDNoBTLevelFlag[ch], 0, m_uiNumPartition * sizeof(*m_ucDIMDNoBTLevelFlag[ch]));
  }
  memset(m_puhSecIntraDir, DC_IDX, m_uiNumPartition * sizeof(*(m_puhSecIntraDir)));
#endif

#if COM16_C806_EMT
  memset( m_puhEmtTuIdx       , 0,                        m_uiNumPartition * sizeof( *m_puhEmtTuIdx ) );
  memset( m_puhEmtCuFlag      , 0,                        m_uiNumPartition * sizeof( *m_puhEmtCuFlag ) );
#endif

#if COM16_C1016_AFFINE
  memset( m_affineFlag        , 0,                        m_uiNumPartition * sizeof( *m_affineFlag ) );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  memset( m_affineParamFlag   , 0,                        m_uiNumPartition * sizeof(*m_affineParamFlag) );
#endif

#if MRG_CONSIDER_COMPLEX
  memset( m_pbMergeComplexFlag  , false, m_uiNumPartition * sizeof( *m_pbMergeComplexFlag ) );
  memset( m_pNumMergeComplexCand, 0,     m_uiNumPartition * sizeof( *m_pNumMergeComplexCand ) );
  memset( m_puhMergeComplexIndex, 0,     m_uiNumPartition * sizeof( *m_puhMergeComplexIndex ) );
#endif

#if RRC_SELECTABLE_IF
  memset(m_pSIFFlag           , 0,                        m_uiNumPartition * sizeof(*m_pSIFFlag) );
#endif

#if RRC_INTRA_SHARPENING
  memset(m_sharpFlag          , 0,                        m_uiNumPartition * sizeof(*m_sharpFlag) );
#endif
  
  for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
  {
    memset( m_puhIntraDir[ch] , ((ch==0) ? DC_IDX : 0),   m_uiNumPartition * sizeof( *(m_puhIntraDir[ch]) ) );
#if MLIP_EXTEND_REF_LINE
    memset(m_puhExternRef[ch], 0, m_uiNumPartition * sizeof(*(m_puhExternRef[ch])));
#endif
  }
  memset( m_puhInterDir       , 0,                        m_uiNumPartition * sizeof( *m_puhInterDir ) );
  memset( m_pbIPCMFlag        , false,                    m_uiNumPartition * sizeof( *m_pbIPCMFlag ) );

  const UInt numCoeffY    = maxCUWidth*maxCUHeight;
  for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    const UInt componentShift = m_pcPic->getComponentScaleX(ComponentID(comp)) + m_pcPic->getComponentScaleY(ComponentID(comp));
    memset( m_pcTrCoeff[comp], 0, sizeof(TCoeff)* numCoeffY>>componentShift );
#if ADAPTIVE_QP_SELECTION
    memset( m_pcArlCoeff[comp], 0, sizeof(TCoeff)* numCoeffY>>componentShift );
#endif
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_acCUMvField[i].clearMvField();
  }

  // Setting neighbor CU
  m_pCtuLeft        = NULL;
  m_pCtuAbove       = NULL;
  m_pCtuAboveLeft   = NULL;
  m_pCtuAboveRight  = NULL;


  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i]  = NULL;
  }

  UInt frameWidthInCtus = pcPic->getFrameWidthInCtus();
  if ( m_ctuRsAddr % frameWidthInCtus )
  {
    m_pCtuLeft = pcPic->getCtu( m_ctuRsAddr - 1 );
  }

  if ( m_ctuRsAddr / frameWidthInCtus )
  {
    m_pCtuAbove = pcPic->getCtu( m_ctuRsAddr - frameWidthInCtus );
  }

  if ( m_pCtuLeft && m_pCtuAbove )
  {
    m_pCtuAboveLeft = pcPic->getCtu( m_ctuRsAddr - frameWidthInCtus - 1 );
  }

  if ( m_pCtuAbove && ( (m_ctuRsAddr%frameWidthInCtus) < (frameWidthInCtus-1) )  )
  {
    m_pCtuAboveRight = pcPic->getCtu( m_ctuRsAddr - frameWidthInCtus + 1 );
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    if ( getSlice()->getNumRefIdx( rpl ) > 0 )
    {
      m_apcCUColocated[rpl] = getSlice()->getRefPic( rpl, 0)->getCtu( m_ctuRsAddr );
    }
  }
}

#if BARC
Void TComDataCU::initDscaleCtu( ChannelType chType, TComPic* pcPic, UInt ctuRsAddr )
{
  const UInt maxCUWidth = pcPic->getPicSym()->getSPS().getCTUSize();
  const UInt maxCUHeight= pcPic->getPicSym()->getSPS().getCTUSize();
  UInt startIdx = 0;
  UInt endIdx = 1;
  
  m_uiNumPartition     = (pcPic->getNumPartitionsInCtu()>> BARC_DScaleRatio) >> BARC_DScaleRatio;
  m_BARCFlag[chType] = true;

  if (chType == CHANNEL_TYPE_CHROMA)
  {
    startIdx = 1;
    endIdx = 2;
    Bool isLumaBARC = pcPic->getCtu(ctuRsAddr)->getBARCFlag(CHANNEL_TYPE_LUMA);
    m_BARCFlag[0] = isLumaBARC;
  }
  for (UInt ch= startIdx; ch < endIdx; ch++)
  {
    memset( m_puhWidth[ch]          , maxCUWidth>> BARC_DScaleRatio,    m_uiNumPartition * sizeof( *m_puhWidth[ch] ) );
    memset( m_puhHeight[ch]         , maxCUHeight>> BARC_DScaleRatio,   m_uiNumPartition * sizeof( *m_puhHeight[ch] ) );
    memset( m_puhDepth[ch]          , BARC_DScaleRatio,                m_uiNumPartition * sizeof( *m_puhDepth[ch] ) );
  }
}

Void TComDataCU::initDscaleCtuDec( ChannelType chType, TComPic* pcPic, UInt ctuRsAddr )
{
  const UInt maxCUWidth = pcPic->getPicSym()->getSPS().getCTUSize();
  const UInt maxCUHeight= pcPic->getPicSym()->getSPS().getCTUSize();
  UInt startIdx = 0;
  UInt endIdx = 1;

  m_uiNumPartition     = pcPic->getNumPartitionsInCtu();
  m_BARCFlag[chType] = true;
  if (chType == CHANNEL_TYPE_CHROMA)
  {
    startIdx = 1;
    endIdx = 2;
  }
  for (UInt ch= startIdx; ch < endIdx; ch++)
  {
    memset( m_puhWidth[ch]          , maxCUWidth>> BARC_DScaleRatio,    m_uiNumPartition/4 * sizeof( *m_puhWidth[ch] ) );
    memset( m_puhHeight[ch]         , maxCUHeight>> BARC_DScaleRatio,   m_uiNumPartition/4 * sizeof( *m_puhHeight[ch] ) );
  }
}

Void TComDataCU::copyCtu (TComDataCU* pCtu)
{
#if JVET_C0024_QTBT
  const UInt maxCUWidth = pCtu->getPic()->getPicSym()->getSPS().getCTUSize();
  const UInt maxCUHeight= pCtu->getPic()->getPicSym()->getSPS().getCTUSize();
#else
  const UInt maxCUWidth = pcPic->getPicSym()->getSPS().getMaxCUWidth();
  const UInt maxCUHeight= pcPic->getPicSym()->getSPS().getMaxCUHeight();
#endif
  m_pcPic              = pCtu->getPic();
  m_pcSlice            = pCtu->getPic()->getSlice(m_pcPic->getCurrSliceIdx());
  m_ctuRsAddr          = pCtu->getCtuRsAddr();
  m_uiCUPelX           = ( m_ctuRsAddr % m_pcPic->getFrameWidthInCtus() ) * maxCUWidth;
  m_uiCUPelY           = ( m_ctuRsAddr / m_pcPic->getFrameWidthInCtus() ) * maxCUHeight;
  m_absZIdxInCtu       = 0;
  m_dTotalCost         = pCtu->getTotalCost();
  m_uiTotalDistortion  = pCtu->getTotalDistortion();
  m_uiTotalBits        = pCtu->getTotalBits();
  m_uiTotalBins        = pCtu->getTotalBins();
  m_uiNumPartition     = m_pcPic->getNumPartitionsInCtu();
#if BARC
  m_BARCFlag[0]        = pCtu->getBARCFlag(CHANNEL_TYPE_LUMA);
  m_BARCFlag[1]        = pCtu->getBARCFlag(CHANNEL_TYPE_CHROMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = pCtu->getUpCNNFlag();
#endif
  for (UInt k = 0; k < m_uiNumPartition; k++)
  {
    m_skipFlag[k] = pCtu->getSkipFlag(k);
#if VCEG_AZ05_INTRA_MPI
    m_MPIIdx[k] = pCtu->getMPIIdx(k); 
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    m_PDPCIdx[k] = pCtu->getPDPCIdx(k);
#endif
#if RRC_DWDIP
    m_DWDIPIdx[k] = pCtu->getDWDIPIdx(k);
#endif
#if RRC_DWDIP_MLR
    m_DWDIP_GPIdx[k] = pCtu->getDWDIP_GPIdx(k);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ )
    {
      m_ROTIdx[i][k] = pCtu->getROTIdx(ChannelType(i),k);
    }
#endif 
    m_pePredMode[k]       = pCtu->getPredictionMode(k);
    m_CUTransquantBypass[k] = pCtu->getCUTransquantBypass(k);

    for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
    {
#if JVET_C0024_DELTA_QP_FIX
      m_phQP[ch][k]    = pCtu->getQP(ChannelType(ch))[k];
#endif
      m_puhDepth[ch][k]  = pCtu->getUpDepth(ch,k);
      m_puhWidth[ch][k]  = pCtu->getUpWidth(ch,k);
      m_puhHeight[ch][k]  = pCtu->getUpHeight(ch,k);
      m_puhBTSplitMode[ch][0][k] = pCtu->getUpBTSplitModePart(ch, 0, k);
      m_puhBTSplitMode[ch][1][k] = pCtu->getUpBTSplitModePart(ch, 1, k);
#if BARC_BTFORCE
      m_puhBTSplitMode[ch][2][k] = pCtu->getUpBTSplitModePart(ch, 2, k);
      m_puhBTSplitMode[ch][3][k] = pCtu->getUpBTSplitModePart(ch, 3, k);
#endif
    }
    m_ChromaQpAdj[k] = pCtu->getChromaQpAdj(k);

    for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      m_crossComponentPredictionAlpha[comp][k] = pCtu->getCrossComponentPredictionAlpha(k, ComponentID(comp));
      m_puhTransformSkip[comp][k] = pCtu->getTransformSkip(k, ComponentID(comp));
      m_puhCbf[comp][k] = pCtu->getCbf(k, ComponentID(comp));
      m_explicitRdpcmMode[comp][k] = pCtu->getExplicitRdpcmMode(ComponentID(comp), k);
#if VCEG_AZ08_KLT_COMMON
      m_puhKLTFlag[comp][k] = pCtu->getKLTFlag(k, ComponentID(comp));
#endif
    }

#if ALF_HM3_REFACTOR
    m_puhAlfCtrlFlag[k] = pCtu->getAlfCtrlFlag(k);
#endif

#if COM16_C806_EMT
    m_puhEmtTuIdx[k] = pCtu->getEmtTuIdx(k);
    m_puhEmtCuFlag[k] = pCtu->getEmtCuFlag(k);
#endif

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
      m_ucDIMDEnabledFlag[ch][k] = pCtu->getDIMDEnabledFlag(ChannelType(ch),k);
#if DIMD_NUM_INTRA_DIR_INC
      m_ucExtIntraDir[ch][k] = pCtu->getExtIntraDir(ChannelType(ch),k);
#endif
      m_ucDIMDNoBTLevelFlag[ch][k] = pCtu->getDIMDNoBTLevelFlag(ChannelType(ch),k);
    }
    m_puhSecIntraDir[k] = pCtu->getSecIntraDir(k);
#endif

#if COM16_C806_EMT
    m_puhEmtTuIdx[k] = pCtu->getEmtTuIdx(k);
    m_puhEmtCuFlag[k] = pCtu->getEmtCuFlag(k);
#endif

#if RRC_INTRA_SHARPENING
    m_sharpFlag[k] = pCtu->getIntraSharpFlag(k);
#endif
    m_puhIntraDir[0][k] = pCtu->getIntraDir(CHANNEL_TYPE_LUMA, k);
#if MLIP_EXTEND_REF_LINE
    m_puhExternRef[0][k] = pCtu->getExternRef(CHANNEL_TYPE_LUMA, k);
#endif
    m_pbIPCMFlag[k] = pCtu->getIPCMFlag(k);
  }

  const UInt numCoeffY = maxCUWidth*maxCUHeight;
  for (UInt ch = 0; ch < 1; ch++)
  {
    const ComponentID component = ComponentID(ch);
    const UInt componentShift = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
    memcpy(m_pcTrCoeff[ch], pCtu->getCoeff(component), sizeof(TCoeff)*(numCoeffY >> componentShift));
#if ADAPTIVE_QP_SELECTION
    memcpy(m_pcArlCoeff[ch], pCtu->getArlCoeff(component), sizeof(TCoeff)*(numCoeffY >> componentShift));
#endif
    memcpy(m_pcIPCMSample[ch], pCtu->getPCMSample(component), sizeof(Pel)*(numCoeffY >> componentShift));
  }
}
#endif

/** Initialize prediction data with enabling sub-CTU-level delta QP.
*   - set CU width and CU height according to depth
*   - set qp value according to input qp
*   - set last-coded qp value according to input last-coded qp
*
* \param  uiDepth            depth of the current CU
* \param  qp                 qp for the current CU
* \param  bTransquantBypass  true for transquant bypass
*/
#if JVET_C0024_QTBT
Void TComDataCU::initEstData( const UInt uiDepth, const Int qp, const Bool bTransquantBypass, UInt uiWidth, UInt uiHeight, Int iBTSplitMode )
#else
Void TComDataCU::initEstData( const UInt uiDepth, const Int qp, const Bool bTransquantBypass )
#endif
{
  m_dTotalCost         = MAX_DOUBLE;
  m_uiTotalDistortion  = 0;
  m_uiTotalBits        = 0;
  m_uiTotalBins        = 0;
#if JVET_C0024_QTBT
  ChannelType eCType   = getTextType();
#if BARC
  m_BARCFlag[eCType] = getPic()->getCtu(getCtuRsAddr())->getBARCFlag(eCType);
  if (!isLuma(eCType))
    m_BARCFlag[CHANNEL_TYPE_LUMA] = getPic()->getCtu(getCtuRsAddr())->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = false;
#endif
  if (uiWidth==0 || uiHeight==0)
  {
    uiHeight = getHeight(0);
    uiWidth = getWidth(0);
  }
  if (iBTSplitMode == -1)
  {
    iBTSplitMode = getBTSplitMode(0);
  }
#else
#if COM16_C806_LARGE_CTU
  const UShort uhWidth  = getSlice()->getSPS()->getMaxCUWidth()  >> uiDepth;
  const UShort uhHeight = getSlice()->getSPS()->getMaxCUHeight() >> uiDepth;
#else
  const UChar uhWidth  = getSlice()->getSPS()->getMaxCUWidth()  >> uiDepth;
  const UChar uhHeight = getSlice()->getSPS()->getMaxCUHeight() >> uiDepth;
#endif
#endif

#if JVET_C0024_QTBT
  for (UInt k=0; k<uiHeight/getPic()->getMinCUHeight(); k++)
  {
    for (UInt j=0; j<uiWidth/getPic()->getMinCUWidth(); j++)    
#else
  for (UInt ui = 0; ui < m_uiNumPartition; ui++)
#endif
  {
#if JVET_C0024_QTBT
#if TT_TREE_DS1
    UInt ui = k * getPic()->getNumPartInCtuWidth() + j;  //raster
    ui = g_auiRasterToZscan[ui];
#else
    UInt ui = k * getPic()->getNumPartInCtuWidth() + j + g_auiZscanToRaster[m_absZIdxInCtu];  //raster
    ui = g_auiRasterToZscan[ui] - m_absZIdxInCtu;
#endif
#endif
    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      const RefPicList rpl=RefPicList(i);
      m_apiMVPIdx[rpl][ui]  = -1;
      m_apiMVPNum[rpl][ui]  = -1;
    }
#if JVET_C0024_QTBT
    m_puhDepth[eCType]  [ui] = uiDepth;
#if JVET_C0024_CTU_256
    m_puhWidth[eCType]  [ui] = uiWidth >> MIN_CU_LOG2;
    m_puhHeight[eCType] [ui] = uiHeight>> MIN_CU_LOG2;
#else
    m_puhWidth[eCType]  [ui] = uiWidth;
    m_puhHeight[eCType] [ui] = uiHeight;
#endif
#if TT_TREE_AFFINE
    m_puhCUabsZIdxInCtu[eCType][0][ui] = m_absZIdxInCtu & 255;
    m_puhCUabsZIdxInCtu[eCType][1][ui] = m_absZIdxInCtu >> 8;
#endif
    m_puhBTSplitMode[eCType][0][ui] = iBTSplitMode & 255;
    m_puhBTSplitMode[eCType][1][ui] = iBTSplitMode>>8; 
#if BT_FORCE_BF
    m_puhBTSplitMode[eCType][2][ui] = iBTSplitMode >> 16;
    m_puhBTSplitMode[eCType][3][ui] = iBTSplitMode >> 24;
#endif
#else
    m_puhDepth  [ui]    = uiDepth;
    m_puhWidth  [ui]    = uhWidth;
    m_puhHeight [ui]    = uhHeight;
    m_puhTrIdx  [ui]    = 0;
#endif
    for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      m_crossComponentPredictionAlpha[comp][ui] = 0;
      m_puhTransformSkip             [comp][ui] = 0;
      m_explicitRdpcmMode            [comp][ui] = NUMBER_OF_RDPCM_MODES;
#if VCEG_AZ08_KLT_COMMON
      m_puhKLTFlag                   [comp][ui] = 0;
#endif
    }
    m_skipFlag[ui]      = false;
#if VCEG_AZ05_INTRA_MPI
    m_MPIIdx[ui]        = 0;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    m_PDPCIdx[ui]       = 0;
#endif

#if RRC_DWDIP
    m_DWDIPIdx[ui] = 0;
#endif
#if RRC_DWDIP_MLR
    m_DWDIP_GPIdx[ui] = 0;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ )
    {
      m_ROTIdx[i][ui] = 0;
    }
#else
    m_ROTIdx[ui] = 0;
#endif
#endif
#if !JVET_C0024_QTBT
    m_pePartSize[ui]    = NUMBER_OF_PART_SIZES;
#endif
    m_pePredMode[ui]    = NUMBER_OF_PREDICTION_MODES;
    m_CUTransquantBypass[ui] = bTransquantBypass;
    m_pbIPCMFlag[ui]    = 0;
#if JVET_C0024_DELTA_QP_FIX
    m_phQP[eCType][ui]  = qp;
#else
    m_phQP[ui]          = qp;
#endif
    m_ChromaQpAdj[ui]   = 0;
    m_pbMergeFlag[ui]   = 0;
    m_puhMergeIndex[ui] = 0;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    m_peMergeType[ui]   = MGR_TYPE_DEFAULT_N;
#endif
#if COM16_C806_OBMC
    m_OBMCFlag[ui]      = true;
#endif
#if VCEG_AZ07_FRUC_MERGE
    m_puhFRUCMgrMode[ui]= 0;
#endif
#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
    m_iMVFlag[ui]       = 0;
#else
    m_iMVFlag[ui]       = false;
#endif
    m_piMVCandNum[ui]   = 0;
#endif
#if VCEG_AZ06_IC
    m_pbICFlag[ui]      = false;
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    m_puhAlfCtrlFlag[ui] = 0;
#else
    m_puiAlfCtrlFlag[ui] = 0;
#endif
#endif


#if DIMD_INTRA_PRED
#if BT_FORCE_DIMD_HARM
#if BT_FORCE_QTN
    Int iMaxQtLv = g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[isLuma(getTextType()) ? getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
    const UInt uiDIMDMaxBtDepht = (getBoundaryCtuLocation() && getSlice()->isIntra()) ? iMaxQtLv : DIMD_MAX_BT_DEPTH_SIG_AI;
#else
    const UInt uiDIMDMaxBtDepht = (getBoundaryCtuLocation() && getSlice()->isIntra()) ? MAX_FQT_LV_INTRA : DIMD_MAX_BT_DEPTH_SIG_AI;
#endif
#endif
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
#if BT_FORCE_DIMD_HARM
      if (!getSlice()->isIntra() || getBTDepth(0) <= uiDIMDMaxBtDepht)
#else
      if (!getSlice()->isIntra() || getBTDepth(0) <= DIMD_MAX_BT_DEPTH_SIG_AI)
#endif
      {
        m_ucDIMDEnabledFlag[ch][ui] = 0;
        m_ucDIMDNoBTLevelFlag[ch][ui] = 0;
      }
#if DIMD_NUM_INTRA_DIR_INC
      m_ucExtIntraDir[ch][ui] = ((ch == 0) ? DC_IDX : 0);
#endif
    }
    m_puhSecIntraDir[ui] = DC_IDX;
#endif

#if COM16_C806_EMT
    m_puhEmtTuIdx [ui]  = 0;
    m_puhEmtCuFlag[ui]  = 0;
#endif

#if COM16_C1016_AFFINE
    m_affineFlag[ui]    = false;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    m_affineParamFlag[ui] = false;
#endif

#if MRG_CONSIDER_COMPLEX
    m_pbMergeComplexFlag[ui]   = 0;
    m_pNumMergeComplexCand[ui] = 0;
    m_puhMergeComplexIndex[ui] = 0;
#endif

#if RRC_INTRA_SHARPENING
    m_sharpFlag[ui]     = 0;
#endif

#if RRC_SELECTABLE_IF
    m_pSIFFlag[ui]      = m_sifToCheck;
#endif

    for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
    {
      m_puhIntraDir[ch][ui] = ((ch==0) ? DC_IDX : 0);
#if MLIP_EXTEND_REF_LINE
      m_puhExternRef[ch][ui] = 0;
#endif
    }

    m_puhInterDir[ui] = 0;
    for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      m_puhCbf[comp][ui] = 0;
    }
  }
#if JVET_C0024_QTBT
  }
#endif

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_acCUMvField[i].clearMvField();
  }

#if JVET_C0024_QTBT
  const UInt numCoeffY = uiWidth*uiHeight;
#else
  const UInt numCoeffY = uhWidth*uhHeight;
#endif

  for (UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    const ComponentID component = ComponentID(comp);
    const UInt numCoeff = numCoeffY >> (getPic()->getComponentScaleX(component) + getPic()->getComponentScaleY(component));
    memset( m_pcTrCoeff[comp],    0, numCoeff * sizeof( TCoeff ) );
#if ADAPTIVE_QP_SELECTION
    memset( m_pcArlCoeff[comp],   0, numCoeff * sizeof( TCoeff ) );
#endif
    memset( m_pcIPCMSample[comp], 0, numCoeff * sizeof( Pel) );
  }
}

#if JVET_C0024_QTBT
Void TComDataCU::initSubBT(TComDataCU* pcCU, UInt uiPartUnitIdx, UInt uiCUDepth, UInt uiBTWidth, UInt uiBTHeight, UInt uiSplitMode, Int qp )
{
#if TT_TREE_DS1 //initSubBT
  if(uiSplitMode>=3 && uiSplitMode<=6)
    assert( uiPartUnitIdx<3 );
  else
#endif
  assert( uiPartUnitIdx<2 );

  m_pcPic              = pcCU->getPic();
  m_pcSlice            = m_pcPic->getSlice(m_pcPic->getCurrSliceIdx());
  m_ctuRsAddr           = pcCU->getCtuRsAddr();

  if (uiSplitMode==1)
  {
    m_uiCUPelX           = pcCU->getCUPelX();
    m_uiCUPelY           = pcCU->getCUPelY() + uiBTHeight * uiPartUnitIdx;
  }
#if TT_ENC_OPT //initSubBT
  else if (uiSplitMode==3 || uiSplitMode==5)
  {
    m_uiCUPelX           = pcCU->getCUPelX();
    if(uiSplitMode==3)
      m_uiCUPelY         = pcCU->getCUPelY() + (uiBTHeight>>1);
    else
      m_uiCUPelY         = pcCU->getCUPelY() + uiBTHeight * (uiPartUnitIdx==0? 0 : 3);
  }
  else if (uiSplitMode==4 || uiSplitMode==6)
  {
    if(uiSplitMode==4)
      m_uiCUPelX         = pcCU->getCUPelX() + (uiBTWidth >> 1);
    else
      m_uiCUPelX         = pcCU->getCUPelX() + uiBTWidth * (uiPartUnitIdx==0? 0 : 3);
    m_uiCUPelY           = pcCU->getCUPelY();
  }
#endif
  else
  {
    assert(uiSplitMode==2);
    m_uiCUPelX           = pcCU->getCUPelX() + uiBTWidth * uiPartUnitIdx;
    m_uiCUPelY           = pcCU->getCUPelY();
  }

  UInt uiXInCTU = m_uiCUPelX % m_pcSlice->getSPS()->getCTUSize();
  UInt uiYInCTU = m_uiCUPelY % m_pcSlice->getSPS()->getCTUSize();

  UInt uiRaster = uiYInCTU/pcCU->getPic()->getMinCUHeight() * pcCU->getPic()->getNumPartInCtuWidth() + uiXInCTU/pcCU->getPic()->getMinCUWidth();
  m_absZIdxInCtu = g_auiRasterToZscan[uiRaster];

  m_dTotalCost         = MAX_DOUBLE;
  m_uiTotalDistortion  = 0;
  m_uiTotalBits        = 0;
  m_uiTotalBins        = 0;
#if TT_ENC_OPT
  if(uiSplitMode==5 || uiSplitMode==6)
    m_uiNumPartition   = pcCU->getTotalNumPart() >> 2;
  else
#endif
  m_uiNumPartition     = pcCU->getTotalNumPart() >> 1;

#if TT_FAST_SMALL_ES || TT_FAST_LARGE_ES
  setSTSplitModeSave( 0, pcCU->getSTSplitModeSave(0,0), pcCU->getSTSplitModeSave(0,1) );
  setSTSplitModeSave( 1, pcCU->getSTSplitModeSave(1,0), pcCU->getSTSplitModeSave(1,1) );
  setSTSplitModeSaveL( 0, pcCU->getSTSplitModeSave(0,2), pcCU->getSTSplitModeSave(0,3) );
  setSTSplitModeSaveL( 1, pcCU->getSTSplitModeSave(1,2), pcCU->getSTSplitModeSave(1,3) );
  m_uiSubCUIdx = uiPartUnitIdx;
#endif

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiWidthSave  = uiBTWidth;
  UInt uiHeightSave = uiBTHeight;
  UInt uiNumXDiv = (uiBTWidth  > m_pcPic->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiBTHeight > m_pcPic->getMinCUHeight())? 2 : 1;
  uiBTWidth  = (uiNumXDiv==1)? uiBTWidth  : (uiBTWidth>>1);
  uiBTHeight = (uiNumYDiv==1)? uiBTHeight : (uiBTHeight>>1);
  UInt uiLenX    = uiBTWidth  / m_pcPic->getMinCUWidth();
  UInt uiLenY    = uiBTHeight / m_pcPic->getMinCUHeight();
  UInt uiRaster2;
  UInt uiRasterSrc = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
  UInt uiRasterDst = g_auiZscanToRaster[m_absZIdxInCtu];
  assert(uiRasterSrc<=uiRasterDst);
  UInt uiDeltaY    = (uiRasterDst - uiRasterSrc) / m_pcPic->getNumPartInCtuWidth();;
  UInt uiDeltaX    = (uiRasterDst - uiRasterSrc) % m_pcPic->getNumPartInCtuWidth();;

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      uiRaster = uiDeltaX + (uiLenX * xIdx) + (uiLenY * yIdx + uiDeltaY) * m_pcPic->getNumPartInCtuWidth();
      uiRaster2= (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
#else
  uiRaster = g_auiZscanToRaster[m_absZIdxInCtu];
#endif
  UInt uiShort, uiLong;
  UInt uiStride;
  if (uiBTHeight > uiBTWidth)
  {
    uiShort = uiBTWidth;
    uiLong = uiBTHeight;
    uiStride = m_pcPic->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiBTHeight;
    uiLong = uiBTWidth;
    uiStride = 1;
  }

  UInt uiDepth = g_aucConvertToBit[m_pcSlice->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;
#if BARC
  m_BARCFlag[getTextType()] = getPic()->getCtu(getCtuRsAddr())->getBARCFlag(getTextType());
  if (!isLuma(getTextType()))
    m_BARCFlag[CHANNEL_TYPE_LUMA] =  getPic()->getCtu(getCtuRsAddr())->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = pcCU->getUpCNNFlag();
#endif
  for (UInt i=0; i<uiLong; i+=uiShort)
  {
#if TT_TREE_DS1
    UInt uiZorderDst = g_auiRasterToZscan[uiRaster2];
    UInt uiZorderSrc = g_auiRasterToZscan[uiRaster];
#else
    UInt uiZorderDst = g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu;
    UInt uiZorderSrc = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
#endif

#if JVET_C0024_DELTA_QP_FIX
    memset(m_phQP[getTextType()] + uiZorderDst, qp, uiCurrPartNumb );
#else
    memset(m_phQP + uiZorderDst, qp, uiCurrPartNumb );
#endif
    memset(m_pbMergeFlag + uiZorderDst, false, uiCurrPartNumb );
    memset(m_puhMergeIndex + uiZorderDst, 0, uiCurrPartNumb );

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    memset( m_peMergeType + uiZorderDst, 0, uiCurrPartNumb );
#endif
#if COM16_C806_OBMC
    memset( m_OBMCFlag + uiZorderDst, 1, uiCurrPartNumb  );
#endif
#if VCEG_AZ07_FRUC_MERGE
    memset( m_puhFRUCMgrMode + uiZorderDst, 0, uiCurrPartNumb );
#endif
#if VCEG_AZ07_IMV
    memset( m_iMVFlag + uiZorderDst, 0, uiCurrPartNumb  );
    memset( m_piMVCandNum + uiZorderDst, 0, uiCurrPartNumb   );
#endif
#if VCEG_AZ06_IC
    memset( m_pbICFlag + uiZorderDst, 0, uiCurrPartNumb  );
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    memset( m_puhAlfCtrlFlag + uiZorderDst, 0, sizeof(UChar) * uiCurrPartNumb );
#else
    memset( m_puiAlfCtrlFlag + uiZorderDst, 0, sizeof(UInt) * uiCurrPartNumb );
#endif
#endif

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
    {
      memset(m_ucDIMDEnabledFlag[ch] + uiZorderDst, 0, uiCurrPartNumb);
#if DIMD_NUM_INTRA_DIR_INC
      memset(m_ucExtIntraDir[ch] + uiZorderDst, ((ch == 0) ? DC_IDX : 0), uiCurrPartNumb);
#endif
      memset(m_ucDIMDNoBTLevelFlag[ch] + uiZorderDst, 0, uiCurrPartNumb);
    }
    memset(m_puhSecIntraDir + uiZorderDst, DC_IDX, uiCurrPartNumb);
#endif

#if COM16_C806_EMT
    memset( m_puhEmtTuIdx + uiZorderDst, 0, uiCurrPartNumb );
    memset( m_puhEmtCuFlag + uiZorderDst, 0, uiCurrPartNumb );
#endif
#if COM16_C1016_AFFINE
    memset( m_affineFlag + uiZorderDst, 0, uiCurrPartNumb  );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    memset( m_affineParamFlag + uiZorderDst, 0, uiCurrPartNumb );
#endif

#if MRG_CONSIDER_COMPLEX
    memset( m_pbMergeComplexFlag + uiZorderDst, false, uiCurrPartNumb );
    memset( m_pNumMergeComplexCand + uiZorderDst, 0, uiCurrPartNumb );
    memset( m_puhMergeComplexIndex + uiZorderDst, 0, uiCurrPartNumb );
#endif

#if VCEG_AZ05_INTRA_MPI
    memset( m_MPIIdx + uiZorderDst, 0, uiCurrPartNumb );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    memset( m_PDPCIdx + uiZorderDst, 0, uiCurrPartNumb );
#endif
#if RRC_DWDIP
    memset(m_DWDIPIdx + uiZorderDst, 0, uiCurrPartNumb);
#endif
#if RRC_DWDIP_MLR
    memset(m_DWDIP_GPIdx + uiZorderDst, 0, uiCurrPartNumb);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++ )
    {
      memset( m_ROTIdx[ch] + uiZorderDst, 0, uiCurrPartNumb );
    }
#else
    memset( m_ROTIdx + uiZorderDst, 0, uiCurrPartNumb );
#endif
#endif

    for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
    {
      memset( m_puhIntraDir[ch] + uiZorderDst,  ((ch==0) ? DC_IDX : 0), uiCurrPartNumb );
#if MLIP_EXTEND_REF_LINE
      memset(m_puhExternRef[ch] + uiZorderDst, 0, uiCurrPartNumb);
#endif
    }
    memset(m_puhInterDir + uiZorderDst, 0, uiCurrPartNumb );
    memset(m_skipFlag + uiZorderDst, false, uiCurrPartNumb );
    memset(m_pePredMode + uiZorderDst, NUMBER_OF_PREDICTION_MODES, uiCurrPartNumb );
    memset(m_CUTransquantBypass + uiZorderDst, 0, uiCurrPartNumb );

    for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
    {
      memset( m_crossComponentPredictionAlpha[comp], 0, uiCurrPartNumb );
      memset( m_puhTransformSkip[comp],              0, uiCurrPartNumb );
      memset( m_puhCbf[comp],                        0, uiCurrPartNumb );
      memset( m_explicitRdpcmMode[comp],             NUMBER_OF_RDPCM_MODES, uiCurrPartNumb );
    }

    memset(m_puhDepth[getTextType()] + uiZorderDst, uiCUDepth, uiCurrPartNumb );
#if TT_TREE_DS1 //data copy recover w & h
#if JVET_C0024_CTU_256
    memset(m_puhWidth[getTextType()] + uiZorderDst, uiWidthSave >> MIN_CU_LOG2, uiCurrPartNumb );
    memset(m_puhHeight[getTextType()] + uiZorderDst, uiHeightSave>> MIN_CU_LOG2, uiCurrPartNumb );
#else
    memset(m_puhWidth[getTextType()] + uiZorderDst, uiWidthSave, uiCurrPartNumb );
    memset(m_puhHeight[getTextType()] + uiZorderDst, uiHeightSave, uiCurrPartNumb );
#endif
#if TT_TREE_AFFINE
    memset(m_puhCUabsZIdxInCtu[getTextType()][0] + uiZorderDst, m_absZIdxInCtu&255, uiCurrPartNumb );
    memset(m_puhCUabsZIdxInCtu[getTextType()][1] + uiZorderDst, m_absZIdxInCtu>>8,  uiCurrPartNumb );
#endif
#else
#if JVET_C0024_CTU_256
    memset(m_puhWidth[getTextType()] + uiZorderDst, uiBTWidth >> MIN_CU_LOG2, uiCurrPartNumb );
    memset(m_puhHeight[getTextType()] + uiZorderDst, uiBTHeight>> MIN_CU_LOG2, uiCurrPartNumb );
#else
    memset(m_puhWidth[getTextType()] + uiZorderDst, uiBTWidth, uiCurrPartNumb );
    memset(m_puhHeight[getTextType()] + uiZorderDst, uiBTHeight, uiCurrPartNumb );
#endif
#endif
    memcpy(m_puhBTSplitMode[getTextType()][0] + uiZorderDst, pcCU->getBTSplitModePart(0) + uiZorderSrc, uiCurrPartNumb );
    memcpy(m_puhBTSplitMode[getTextType()][1] + uiZorderDst, pcCU->getBTSplitModePart(1) + uiZorderSrc, uiCurrPartNumb );
#if BT_FORCE_BF
    memcpy(m_puhBTSplitMode[getTextType()][2] + uiZorderDst, pcCU->getBTSplitModePart(2) + uiZorderSrc, uiCurrPartNumb);
    memcpy(m_puhBTSplitMode[getTextType()][3] + uiZorderDst, pcCU->getBTSplitModePart(3) + uiZorderSrc, uiCurrPartNumb);
#endif

    memset(m_pbIPCMFlag + uiZorderDst, 0, uiCurrPartNumb );

    for(UInt refListIdx=0; refListIdx<NUM_REF_PIC_LIST_01; refListIdx++)
    {
      const RefPicList rpl=RefPicList(refListIdx);
      memset(m_apiMVPIdx[rpl] + uiZorderDst, -1, uiCurrPartNumb );
      memset(m_apiMVPNum[rpl] + uiZorderDst, -1, uiCurrPartNumb );
    }

#if RRC_INTRA_SHARPENING
    memset(m_sharpFlag + uiZorderDst, 0, uiCurrPartNumb);
#endif

    uiRaster += uiNumPartInShort * uiStride;
#if TT_TREE_DS1
    uiRaster2 += uiNumPartInShort * uiStride;
#endif
  }
#if TT_TREE_DS1 //data copy end
  }}
  //recover uiBTWidth & uiBTHeight
  uiBTWidth = uiWidthSave;
  uiBTHeight = uiHeightSave;
#endif

  UInt numCoeffY = uiBTWidth*uiBTHeight;
  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    const UInt componentShift = m_pcPic->getComponentScaleX(ComponentID(ch)) + m_pcPic->getComponentScaleY(ComponentID(ch));
    memset( m_pcTrCoeff[ch],  0, sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
    memset( m_pcArlCoeff[ch], 0, sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
    memset( m_pcIPCMSample[ch], 0, sizeof(Pel)* (numCoeffY>>componentShift) );
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_acCUMvField[i].clearMvField();
  }

  m_pCtuLeft        = pcCU->getCtuLeft();
  m_pCtuAbove       = pcCU->getCtuAbove();
  m_pCtuAboveLeft   = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight  = pcCU->getCtuAboveRight();

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i] = pcCU->getCUColocated(RefPicList(i));
  }
#if RRC_SELECTABLE_IF
  m_sifToCheck = pcCU->m_sifToCheck;
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  UInt uiBTDepth = getBTDepth(0);
#if BT_FORCE_BF
  assert(uiBTDepth >= 1);
#else
  assert(uiBTDepth >= 1 && uiBTDepth <= 3); //now only support MaxBTD = 3
#endif
  m_uiSubCUIdxCompST = (uiPartUnitIdx<<((uiBTDepth-1)*4)) + pcCU->getSubCUIdxComp();
  deriveCURouteInfo();
#endif
}

#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
Void TComDataCU::updateLeftUpNeighbor(Int* x, Int* y, UInt uiBTDepth)
{
  if(isLeftNeighborInsideCurPart(x,y))
  {
    UInt uiWIdx = g_aucConvertToBit[x[1]-x[0]];
    UInt uiHIdx = g_aucConvertToBit[y[1]-y[0]];
    UInt uiStride = getPic()->getNumPartInCtuWidth();
    UInt uiZIdx = m_uiQTLeafZorderIdx + g_auiRasterToZscan[ (x[0]>>MIN_CU_LOG2) + (y[0]>>MIN_CU_LOG2) * uiStride ];
    m_uiCURouteInfo[1] = (uiZIdx<<8) + (uiWIdx<<5) + (uiHIdx<<2) + uiBTDepth;
    return;
  }

  if(isUpNeighborInsideCurPart(x,y))
  {
    UInt uiWIdx = g_aucConvertToBit[x[1]-x[0]];
    UInt uiHIdx = g_aucConvertToBit[y[1]-y[0]];
    UInt uiStride = getPic()->getNumPartInCtuWidth();
    UInt uiZIdx = m_uiQTLeafZorderIdx + g_auiRasterToZscan[ (x[0]>>MIN_CU_LOG2) + (y[0]>>MIN_CU_LOG2) * uiStride ];
    m_uiCURouteInfo[2] = (uiZIdx<<8) + (uiWIdx<<5) + (uiHIdx<<2) + uiBTDepth;
    return;
  }
}

Bool TComDataCU::isLeftNeighborInsideCurPart(Int* x, Int* y)
{
  if(m_iXInQTLeaf-1<x[0] || m_iXInQTLeaf-1>=x[1])
    return false;
  else if(m_iYInQTLeaf<y[0] || m_iYInQTLeaf>=y[1])
    return false;
  return true;
}

Bool TComDataCU::isUpNeighborInsideCurPart(Int* x, Int* y)
{
  if(m_iYInQTLeaf-1<y[0] || m_iYInQTLeaf-1>=y[1])
    return false;
  else if(m_iXInQTLeaf<x[0] || m_iXInQTLeaf>=x[1])
    return false;
  return true;
}

Void TComDataCU::deriveCURouteInfo()
{
  //todo
  int  i,j;
  UInt uiQTDepth = getDepth(0);
  UInt uiBTDepth = getBTDepth(0);
  //UInt uiBTSplitModeComp = getBTSplitMode(0);
  UInt uiSubCUIdxComp    = getSubCUIdxComp();
  assert(uiBTDepth>=1);
  UInt uiZorderIdx = getZorderIdxInCtu();
  UInt uiCTUSize = getSlice()->getSPS()->getCTUSize();
  UInt uiQTLeafSize = uiCTUSize>>uiQTDepth;
  UInt uiQTLeafSizeInPart = uiCTUSize >> (uiQTDepth + MIN_CU_LOG2);
  UInt uiRasterIdxCTU = g_auiZscanToRaster[uiZorderIdx];
  UInt uiXIdx = uiRasterIdxCTU % getPic()->getNumPartInCtuWidth(); //XIdx in CTU
  UInt uiYIdx = uiRasterIdxCTU / getPic()->getNumPartInCtuWidth(); //YIdx in CTU
  m_uiQTLeafZorderIdx = g_auiRasterToZscan[(uiYIdx/uiQTLeafSizeInPart)*uiQTLeafSizeInPart*getPic()->getNumPartInCtuWidth() + (uiXIdx/uiQTLeafSizeInPart)*uiQTLeafSizeInPart];
  m_iXInQTLeaf = (uiXIdx % uiQTLeafSizeInPart)<<MIN_CU_LOG2;
  m_iYInQTLeaf = (uiYIdx % uiQTLeafSizeInPart)<<MIN_CU_LOG2;
  UInt uiBTSplitMode, uiSubCUIdx;
  UInt uiWIdxL, uiHIdxL, uiWIdxU, uiHIdxU;
  uiWIdxL = uiHIdxL = uiWIdxU = uiHIdxU = g_aucConvertToBit[uiQTLeafSize];
  UInt uiSTDir;  //0: Horizontal, 1: Vertical
  UInt uiSTType; //0: BT, 1: TT
  UInt uiParentCUNumPart, uiCurrCUNumPart;
  Int xParent[2], yParent[2], xCurr[2], yCurr[2], xNext[2], yNext[2];
  Int iWidthParent, iHeightParent;
  uiParentCUNumPart = uiQTLeafSizeInPart*uiQTLeafSizeInPart;
  xParent[0] = 0; xParent[1] = uiQTLeafSize;
  yParent[0] = 0; yParent[1] = uiQTLeafSize;

  //0: CodingOrder in 4x4, 1: LeftInfo [ZIdx 24-8][WIdx 7-5][HIdx 4-2][STDepth 1-0], 2: UpperInfo  [ZIdx 24-8][WIdx 7-5][HIdx 4-2][STDepth 1-0]
  //set Left and Up as default
  m_uiCURouteInfo[0] = 0;
  m_uiCURouteInfo[1] = (uiParentCUNumPart<<8) + (uiWIdxL<<5) + (uiHIdxL<<2) + 0;
  m_uiCURouteInfo[2] = (uiParentCUNumPart<<8) + (uiWIdxU<<5) + (uiHIdxU<<2) + 0;

  //update Left and Up info based on BTSplit info
  for(i=0; i<uiBTDepth; i++)
  {
    uiBTSplitMode = getBTSplitModeForBTDepth(0, i);
    uiSubCUIdx    = (uiSubCUIdxComp >> (4*i)) & 15;
#if TT_TREE_DS1
    uiSTDir = (uiBTSplitMode>>2)==1? 0 : 1;
    uiSTType = (uiBTSplitMode>>1) & 1;
#else
    uiSTDir = uiBTSplitMode==1? 0 : 1;
    uiSTType = 0;
#endif
    iWidthParent  = xParent[1]-xParent[0];
    iHeightParent = yParent[1]-yParent[0];

    if(uiSTType==0)
    {
      //update CodingOrder
      uiCurrCUNumPart     = (uiParentCUNumPart>>1);
      m_uiCURouteInfo[0] += (uiParentCUNumPart>>1) * uiSubCUIdx;

      //check whether Left and Up is the subCU
      for(j=0; j<2; j++)
      {
        if(uiSTDir==0) //HBT
        {
          xCurr[0] = xParent[0];  xCurr[1] = xParent[1];
          yCurr[0] = yParent[0] + (iHeightParent>>1) * j;  
          yCurr[1] = yCurr[0]   + (iHeightParent>>1);
        }
        else //VBT
        {
          yCurr[0] = yParent[0];  yCurr[1] = yParent[1];
          xCurr[0] = xParent[0] + (iWidthParent>>1) * j;  
          xCurr[1] = xCurr[0]   + (iWidthParent>>1);
        }
        updateLeftUpNeighbor(xCurr, yCurr, i+1);
        if(uiSubCUIdx==j) 
        {
          xNext[0] = xCurr[0];  xNext[1] = xCurr[1];  
          yNext[0] = yCurr[0];  yNext[1] = yCurr[1];
        }
      }
    }
    else
    {
      //update CodingOrder
      uiCurrCUNumPart     = (uiParentCUNumPart>>2) * (uiSubCUIdx==1? 2 : 1);
      m_uiCURouteInfo[0] += (uiParentCUNumPart>>2) * (uiSubCUIdx==2? 3 : uiSubCUIdx);

      //check whether Left and Up is the subCU
      for(j=0; j<3; j++)
      {
        if(uiSTDir==0) //HTT
        {
          xCurr[0] = xParent[0];  xCurr[1] = xParent[1];
          yCurr[0] = yParent[0] + (iHeightParent>>2) * (j==2? 3 : j);  
          yCurr[1] = yCurr[0]   + (j==1? (iHeightParent>>1) : (iHeightParent>>2));
        }
        else //VTT
        {
          yCurr[0] = yParent[0];  yCurr[1] = yParent[1];
          xCurr[0] = xParent[0] + (iWidthParent>>2) * (j==2? 3 : j);  
          xCurr[1] = xCurr[0]   + (j==1? (iWidthParent>>1) : (iWidthParent>>2));
        }
        updateLeftUpNeighbor(xCurr, yCurr, i+1);
        if(uiSubCUIdx==j) 
        {
          xNext[0] = xCurr[0];  xNext[1] = xCurr[1];  
          yNext[0] = yCurr[0];  yNext[1] = yCurr[1];
        }
      }
    }

    //update current CU as parent CU
    uiParentCUNumPart = uiCurrCUNumPart;
    xParent[0] = xNext[0]; xParent[1] = xNext[1];
    yParent[0] = yNext[0]; yParent[1] = yNext[1];
  }
  assert(uiZorderIdx == m_uiQTLeafZorderIdx + g_auiRasterToZscan[(yParent[0]>>MIN_CU_LOG2)*getPic()->getNumPartInCtuWidth() + (xParent[0]>>MIN_CU_LOG2)] );
}
#endif

#endif

// initialize Sub partition
Void TComDataCU::initSubCU( TComDataCU* pcCU, UInt uiPartUnitIdx, UInt uiDepth, Int qp )
{
  assert( uiPartUnitIdx<4 );

  UInt uiPartOffset = ( pcCU->getTotalNumPart()>>2 )*uiPartUnitIdx;

  m_pcPic              = pcCU->getPic();
  m_pcSlice            = pcCU->getSlice();
  m_ctuRsAddr          = pcCU->getCtuRsAddr();
  m_absZIdxInCtu       = pcCU->getZorderIdxInCtu() + uiPartOffset;
#if JVET_C0024_QTBT
  ChannelType   eCType = getTextType();

  const UChar uhWidth  = getSlice()->getSPS()->getCTUSize()  >> uiDepth;
  const UChar uhHeight = getSlice()->getSPS()->getCTUSize() >> uiDepth;
#else
#if COM16_C806_LARGE_CTU
  const UShort uhWidth  = getSlice()->getSPS()->getMaxCUWidth()  >> uiDepth;
  const UShort uhHeight = getSlice()->getSPS()->getMaxCUHeight() >> uiDepth;
#else
  const UChar uhWidth  = getSlice()->getSPS()->getMaxCUWidth()  >> uiDepth;
  const UChar uhHeight = getSlice()->getSPS()->getMaxCUHeight() >> uiDepth;
#endif
#endif

  m_uiCUPelX           = pcCU->getCUPelX() + ( uhWidth )*( uiPartUnitIdx &  1 );
  m_uiCUPelY           = pcCU->getCUPelY() + ( uhHeight)*( uiPartUnitIdx >> 1 );

  m_dTotalCost         = MAX_DOUBLE;
  m_uiTotalDistortion  = 0;
  m_uiTotalBits        = 0;
  m_uiTotalBins        = 0;
  m_uiNumPartition     = pcCU->getTotalNumPart() >> 2;

  Int iSizeInUchar = sizeof( UChar  ) * m_uiNumPartition;
  Int iSizeInBool  = sizeof( Bool   ) * m_uiNumPartition;
  Int sizeInChar = sizeof( Char  ) * m_uiNumPartition;
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
  m_uiSubCUIdxCompST = 0;
#endif
#if BARC
  m_BARCFlag[eCType] = getPic()->getCtu(getCtuRsAddr())->getBARCFlag(eCType);
  if (!isLuma(eCType))
     m_BARCFlag[CHANNEL_TYPE_LUMA] = getPic()->getCtu(getCtuRsAddr())->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = pcCU->getUpCNNFlag();
#endif
#if JVET_C0024_DELTA_QP_FIX
  memset( m_phQP[eCType],      qp,  sizeInChar );
#else
  memset( m_phQP,              qp,  sizeInChar );
#endif
  memset( m_pbMergeFlag,        0, iSizeInBool  );
  memset( m_puhMergeIndex,      0, iSizeInUchar );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  memset( m_peMergeType,        0, iSizeInUchar );
#endif
#if COM16_C806_OBMC
  memset( m_OBMCFlag,           1, iSizeInBool  );
#endif
#if VCEG_AZ07_FRUC_MERGE
  memset( m_puhFRUCMgrMode,     0, iSizeInUchar );
#endif
#if VCEG_AZ07_IMV
  memset( m_iMVFlag,            0, iSizeInBool  );
  memset( m_piMVCandNum,        0, sizeInChar   );
#endif
#if VCEG_AZ06_IC
  memset( m_pbICFlag,           0, iSizeInBool  );
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
  memset( m_puhAlfCtrlFlag,     0, sizeof( UChar   ) * m_uiNumPartition );
#else
  memset( m_puiAlfCtrlFlag,     0, sizeof( UInt   ) * m_uiNumPartition );
#endif
#endif

#if DIMD_INTRA_PRED
  for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
  {
    memset(m_ucDIMDEnabledFlag[ch], 0, iSizeInUchar);
#if DIMD_NUM_INTRA_DIR_INC
    memset(m_ucExtIntraDir[ch], ((ch == 0) ? DC_IDX : 0), iSizeInUchar);
#endif
    memset(m_ucDIMDNoBTLevelFlag[ch], 0, iSizeInUchar);
  }
  memset(m_puhSecIntraDir, DC_IDX, iSizeInUchar);
#endif

#if COM16_C806_EMT
  memset( m_puhEmtTuIdx,        0, iSizeInUchar );
  memset( m_puhEmtCuFlag,       0, iSizeInUchar );
#endif

#if COM16_C1016_AFFINE
  memset( m_affineFlag,         0, iSizeInBool  );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  memset( m_affineParamFlag, 0, iSizeInBool);
#endif

#if MRG_CONSIDER_COMPLEX
  memset( m_pbMergeComplexFlag,   0, iSizeInBool  );
  memset( m_pNumMergeComplexCand, 0, iSizeInUchar );
  memset( m_puhMergeComplexIndex, 0, iSizeInUchar );
#endif

#if RRC_INTRA_SHARPENING
  memset(m_sharpFlag, 0, sizeof(*m_sharpFlag)   * m_uiNumPartition);
#endif

  for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
  {
    memset( m_puhIntraDir[ch],  ((ch==0) ? DC_IDX : 0), iSizeInUchar );
#if MLIP_EXTEND_REF_LINE
    memset(m_puhExternRef[ch], 0, iSizeInUchar);
#endif
  }

  memset( m_puhInterDir,        0, iSizeInUchar );
#if !JVET_C0024_QTBT
  memset( m_puhTrIdx,           0, iSizeInUchar );
#endif

  for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    memset( m_crossComponentPredictionAlpha[comp], 0, iSizeInUchar );
    memset( m_puhTransformSkip[comp],              0, iSizeInUchar );
    memset( m_puhCbf[comp],                        0, iSizeInUchar );
    memset( m_explicitRdpcmMode[comp],             NUMBER_OF_RDPCM_MODES, iSizeInUchar );
#if VCEG_AZ08_KLT_COMMON
    memset( m_puhKLTFlag[comp],                    0, iSizeInUchar);
#endif
  }

#if JVET_C0024_QTBT
  memset( m_puhDepth[eCType],     uiDepth, iSizeInUchar );
#if JVET_C0024_CTU_256
  memset( m_puhWidth[eCType],     uhWidth >> MIN_CU_LOG2, iSizeInUchar );
  memset( m_puhHeight[eCType],    uhHeight>> MIN_CU_LOG2, iSizeInUchar );
#else
  memset( m_puhWidth[eCType],     uhWidth, iSizeInUchar );
  memset( m_puhHeight[eCType],    uhHeight, iSizeInUchar );
#endif
#if TT_TREE_AFFINE
  memset( m_puhCUabsZIdxInCtu[eCType][0],    m_absZIdxInCtu&255, iSizeInUchar);
  memset( m_puhCUabsZIdxInCtu[eCType][1],    m_absZIdxInCtu>>8,  iSizeInUchar);
#endif
  memset( m_puhBTSplitMode[eCType][0],    0, iSizeInUchar);
  memset( m_puhBTSplitMode[eCType][1],    0, iSizeInUchar);
#if BT_FORCE_BF
  memset(m_puhBTSplitMode[eCType][2], 0, iSizeInUchar);
  memset(m_puhBTSplitMode[eCType][3], 0, iSizeInUchar);
#endif
#else
  memset( m_puhDepth,     uiDepth, iSizeInUchar );
#if COM16_C806_LARGE_CTU
  for( Int n = 0 ; n < iSizeInUchar ; n++ )
  {
    m_puhWidth[n] = uhWidth;
    m_puhHeight[n] = uhHeight;
  }
#else
  memset( m_puhWidth,          uhWidth,  iSizeInUchar );
  memset( m_puhHeight,         uhHeight, iSizeInUchar );
#endif
#endif
  memset( m_pbIPCMFlag,        0, iSizeInBool  );

  for (UInt ui = 0; ui < m_uiNumPartition; ui++)
  {
    m_skipFlag[ui]   = false;
#if VCEG_AZ05_INTRA_MPI
    m_MPIIdx[ui]     = 0;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    m_PDPCIdx[ui] = 0;
#endif
#if RRC_DWDIP
    m_DWDIPIdx[ui] = 0;
#endif
#if RRC_DWDIP_MLR
    m_DWDIP_GPIdx[ui] = 0;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
    {
      m_ROTIdx[i][ui] = 0;
    }
#else
    m_ROTIdx[ui] = 0;
#endif
#endif

#if COM16_C1016_AFFINE
    m_affineFlag[ui] = false;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    m_affineParamFlag[ui] = false;
#endif

#if MRG_CONSIDER_COMPLEX
    m_pbMergeComplexFlag[ui]   = 0;
    m_pNumMergeComplexCand[ui] = 0;
    m_puhMergeComplexIndex[ui] = 0;
#endif

#if !JVET_C0024_QTBT
    m_pePartSize[ui] = NUMBER_OF_PART_SIZES;
#endif
    m_pePredMode[ui] = NUMBER_OF_PREDICTION_MODES;
    m_CUTransquantBypass[ui] = false;
    m_ChromaQpAdj[ui] = 0;

    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      const RefPicList rpl=RefPicList(i);
      m_apiMVPIdx[rpl][ui] = -1;
      m_apiMVPNum[rpl][ui] = -1;
    }
  }

  const UInt numCoeffY    = uhWidth*uhHeight;
  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    const UInt componentShift = m_pcPic->getComponentScaleX(ComponentID(ch)) + m_pcPic->getComponentScaleY(ComponentID(ch));
    memset( m_pcTrCoeff[ch],  0, sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
    memset( m_pcArlCoeff[ch], 0, sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
    memset( m_pcIPCMSample[ch], 0, sizeof(Pel)* (numCoeffY>>componentShift) );
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_acCUMvField[i].clearMvField();
  }

  m_pCtuLeft        = pcCU->getCtuLeft();
  m_pCtuAbove       = pcCU->getCtuAbove();
  m_pCtuAboveLeft   = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight  = pcCU->getCtuAboveRight();

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i] = pcCU->getCUColocated(RefPicList(i));
  }
#if RRC_SELECTABLE_IF
  m_sifToCheck = pcCU->m_sifToCheck;
#endif
}

Void TComDataCU::setOutsideCUPart( UInt uiAbsPartIdx, UInt uiDepth )
{
  const UInt     uiNumPartition = m_uiNumPartition >> (uiDepth << 1);
  const UInt     uiSizeInUchar  = sizeof( UChar  ) * uiNumPartition;
  const TComSPS &sps            = *(getSlice()->getSPS());
#if JVET_C0024_QTBT
  UShort uhWidth  = sps.getCTUSize()  >> uiDepth;
  UShort uhHeight = sps.getCTUSize() >> uiDepth;
  ChannelType eCType = getTextType();
  memset( m_puhDepth[eCType]    + uiAbsPartIdx,     uiDepth,  uiSizeInUchar );
#if JVET_C0024_CTU_256
  memset( m_puhWidth[eCType]    + uiAbsPartIdx,     uhWidth >> MIN_CU_LOG2,  uiSizeInUchar );
  memset( m_puhHeight[eCType]    + uiAbsPartIdx,     uhHeight >> MIN_CU_LOG2,  uiSizeInUchar );
#else
  memset( m_puhWidth[eCType]    + uiAbsPartIdx,     uhWidth,  uiSizeInUchar );
  memset( m_puhHeight[eCType]    + uiAbsPartIdx,     uhHeight,  uiSizeInUchar );
#endif
#else
#if COM16_C806_LARGE_CTU
  UShort uhWidth  = sps.getMaxCUWidth()  >> uiDepth;
  UShort uhHeight = sps.getMaxCUHeight() >> uiDepth;
  memset( m_puhDepth    + uiAbsPartIdx,     uiDepth,  uiSizeInUchar );
  for( UInt n = 0 ; n < uiSizeInUchar ; n++ )
  {
    m_puhWidth[uiAbsPartIdx+n] = uhWidth;
    m_puhHeight[uiAbsPartIdx+n] = uhHeight;
  }
#else
  const UChar    uhWidth        = sps.getMaxCUWidth()  >> uiDepth;
  const UChar    uhHeight       = sps.getMaxCUHeight() >> uiDepth;
  memset( m_puhDepth    + uiAbsPartIdx,     uiDepth,  uiSizeInUchar );
  memset( m_puhWidth    + uiAbsPartIdx,     uhWidth,  uiSizeInUchar );
  memset( m_puhHeight   + uiAbsPartIdx,     uhHeight, uiSizeInUchar );
#endif
#endif
}

// --------------------------------------------------------------------------------------------------------------------
// Copy
// --------------------------------------------------------------------------------------------------------------------

Void TComDataCU::copySubCU( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiPart = uiAbsPartIdx;

  m_pcPic              = pcCU->getPic();
  m_pcSlice            = pcCU->getSlice();
  m_ctuRsAddr          = pcCU->getCtuRsAddr();
  m_absZIdxInCtu       = uiAbsPartIdx;
#if JVET_C0024_QTBT
  ChannelType eCType  = getTextType();
#endif
#if BARC
  m_BARCFlag[0] = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
  m_BARCFlag[1] = pcCU->getBARCFlag(CHANNEL_TYPE_CHROMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = pcCU->getUpCNNFlag();
#endif
  m_uiCUPelX           = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  m_uiCUPelY           = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];

  m_skipFlag           =pcCU->getSkipFlag()          + uiPart;

#if VCEG_AZ05_INTRA_MPI
  m_MPIIdx             = pcCU->getMPIIdx()           + uiPart;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  m_PDPCIdx            = pcCU->getPDPCIdx()          + uiPart;
#endif           
#if RRC_DWDIP
  m_DWDIPIdx           = pcCU->getDWDIPIdx()         + uiPart;
#endif          
#if RRC_DWDIP_MLR
  m_DWDIP_GPIdx            = pcCU->getDWDIP_GPIdx()          + uiPart;
#endif 
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
  for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ )
  {
    m_ROTIdx[i]       = pcCU->getROTIdx(ChannelType(i)) + uiPart;
  }
#else
  m_ROTIdx            = pcCU->getROTIdx()           + uiPart;
#endif
#endif
#if JVET_C0024_DELTA_QP_FIX
  m_phQP[eCType] = pcCU->getQP()          + uiPart;
#else
  m_phQP=pcCU->getQP()                    + uiPart;
#endif
  m_ChromaQpAdj = pcCU->getChromaQpAdj()  + uiPart;
#if !JVET_C0024_QTBT
  m_pePartSize = pcCU->getPartitionSize() + uiPart;
#endif
  m_pePredMode=pcCU->getPredictionMode()  + uiPart;
  m_CUTransquantBypass  = pcCU->getCUTransquantBypass()+uiPart;

  m_pbMergeFlag         = pcCU->getMergeFlag()        + uiPart;
  m_puhMergeIndex       = pcCU->getMergeIndex()       + uiPart;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  m_peMergeType         = pcCU->getMergeType()        + uiPart;
#endif
#if COM16_C806_OBMC
  m_OBMCFlag            = pcCU->getOBMCFlag()         + uiPart;
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_puhFRUCMgrMode      = pcCU->m_puhFRUCMgrMode      + uiPart;
#endif
#if VCEG_AZ07_IMV
  m_iMVFlag             = pcCU->getiMVFlag()          + uiPart;
  m_piMVCandNum         = pcCU->getiMVCandNum()       + uiPart;
#endif
#if VCEG_AZ06_IC
  m_pbICFlag            = pcCU->getICFlag()           + uiPart;
#endif
  for (UInt ch=0; ch<MAX_NUM_CHANNEL_TYPE; ch++)
  {
    m_puhIntraDir[ch]   = pcCU->getIntraDir(ChannelType(ch)) + uiPart;
#if MLIP_EXTEND_REF_LINE
    m_puhExternRef[ch] = pcCU->getExternRef(ChannelType(ch)) + uiPart;
#endif
  }

  m_puhInterDir         = pcCU->getInterDir()         + uiPart;
#if !JVET_C0024_QTBT
  m_puhTrIdx            = pcCU->getTransformIdx()     + uiPart;
#endif

#if DIMD_INTRA_PRED
  for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
  {
    m_ucDIMDEnabledFlag[ch] = pcCU->getDIMDEnabledFlag(ChannelType(ch)) + uiPart;
#if DIMD_NUM_INTRA_DIR_INC
    m_ucExtIntraDir[ch] = pcCU->getExtIntraDir(ChannelType(ch)) + uiPart;
#endif
    m_ucDIMDNoBTLevelFlag[ch] = pcCU->getDIMDNoBTLevelFlag(ChannelType(ch)) + uiPart;
  }
  m_puhSecIntraDir = pcCU->getSecIntraDir() + uiPart;
#endif

#if COM16_C806_EMT
  m_puhEmtTuIdx         = pcCU->getEmtTuIdx()         + uiPart;
  m_puhEmtCuFlag        = pcCU->getEmtCuFlag()        + uiPart;
#endif

#if COM16_C1016_AFFINE
  m_affineFlag          = pcCU->getAffineFlag()       + uiPart;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  m_affineParamFlag     = pcCU->getAffineParamFlag()  + uiPart;
#endif

#if MRG_CONSIDER_COMPLEX
  m_pbMergeComplexFlag   = pcCU->getMergeComplexFlag()    + uiPart;
  m_pNodeMotionInfo      = pcCU->getNodeMotionInfo()      + uiPart;
  m_pNumMergeComplexCand = pcCU->getNumMergeComplexCand() + uiPart;
  m_puhMergeComplexIndex = pcCU->getMergeComplexIndex()   + uiPart;
#endif

#if RRC_SELECTABLE_IF
  m_pSIFFlag            = pcCU->getSIFFlag()          + uiPart;
#endif

#if RRC_INTRA_SHARPENING
  m_sharpFlag = pcCU->getIntraSharpFlag()                  + uiPart;
#endif

  for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    m_crossComponentPredictionAlpha[comp] = pcCU->getCrossComponentPredictionAlpha(ComponentID(comp)) + uiPart;
    m_puhTransformSkip[comp]              = pcCU->getTransformSkip(ComponentID(comp))                 + uiPart;
    m_puhCbf[comp]                        = pcCU->getCbf(ComponentID(comp))                           + uiPart;
    m_explicitRdpcmMode[comp]             = pcCU->getExplicitRdpcmMode(ComponentID(comp))             + uiPart;
#if VCEG_AZ08_KLT_COMMON
    m_puhKLTFlag[comp]                    = pcCU->getKLTFlag(ComponentID(comp))                       + uiPart;
#endif

#if JVET_C0024_ITSKIP
    m_puiSkipWidth[comp]                  = pcCU->getTUSkipWidth(ComponentID(comp)) + uiPart;
    m_puiSkipHeight[comp]                 = pcCU->getTUSkipHeight(ComponentID(comp)) + uiPart;
#endif
  }

#if JVET_C0024_QTBT
  m_puhDepth[eCType]=pcCU->getDepth()                     + uiPart;
  m_puhWidth[eCType]=pcCU->getWidth()                     + uiPart;
  m_puhHeight[eCType]=pcCU->getHeight()                   + uiPart;
  m_puhBTSplitMode[eCType][0] = pcCU->getBTSplitModePart(0) + uiPart;
  m_puhBTSplitMode[eCType][1] = pcCU->getBTSplitModePart(1) + uiPart;
#if BT_FORCE_BF
  m_puhBTSplitMode[eCType][2] = pcCU->getBTSplitModePart(2) + uiPart;
  m_puhBTSplitMode[eCType][3] = pcCU->getBTSplitModePart(3) + uiPart;
#endif
#if TT_TREE_AFFINE
  m_puhCUabsZIdxInCtu[eCType][0] = pcCU->getCUAbsZIdxInCtu(0) + uiPart;
  m_puhCUabsZIdxInCtu[eCType][1] = pcCU->getCUAbsZIdxInCtu(1) + uiPart;
#endif
#else
  m_puhDepth=pcCU->getDepth()                     + uiPart;
  m_puhWidth=pcCU->getWidth()                     + uiPart;
  m_puhHeight=pcCU->getHeight()                   + uiPart;
#endif

  m_pbIPCMFlag         = pcCU->getIPCMFlag()        + uiPart;

  m_pCtuAboveLeft      = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight     = pcCU->getCtuAboveRight();
  m_pCtuAbove          = pcCU->getCtuAbove();
  m_pCtuLeft           = pcCU->getCtuLeft();

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    m_apcCUColocated[rpl] = pcCU->getCUColocated(rpl);
    m_apiMVPIdx[rpl]=pcCU->getMVPIdx(rpl)  + uiPart;
    m_apiMVPNum[rpl]=pcCU->getMVPNum(rpl)  + uiPart;
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    m_acCUMvField[rpl].linkToWithOffset( pcCU->getCUMvField(rpl), uiPart );
#if TT_TREE_DS3 //save address of CTU MV field
    // NOTE: at decoder, we shall use absolute ZIdx related to CTU not to CU
    m_acCUMvFieldCtu[rpl].linkToWithOffset( pcCU->getCUMvField(rpl), 0 );
    m_bUseAbsoluteZIdxToCtu = true;
    m_pCtuCurr = pcCU;
#endif
  }

#if JVET_C0024_QTBT 
  UInt uiCoffOffset = pcCU->getPic()->getCodedAreaInCTU();
#else
  UInt uiMaxCuWidth=pcCU->getSlice()->getSPS()->getMaxCUWidth();
  UInt uiMaxCuHeight=pcCU->getSlice()->getSPS()->getMaxCUHeight();

  UInt uiCoffOffset = uiMaxCuWidth*uiMaxCuHeight*uiAbsPartIdx/pcCU->getPic()->getNumPartitionsInCtu();
#endif

  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    const ComponentID component = ComponentID(ch);
    const UInt componentShift   = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
    const UInt offset           = uiCoffOffset >> componentShift;
    m_pcTrCoeff[ch] = pcCU->getCoeff(component) + offset;
#if ADAPTIVE_QP_SELECTION
    m_pcArlCoeff[ch] = pcCU->getArlCoeff(component) + offset;
#endif
    m_pcIPCMSample[ch] = pcCU->getPCMSample(component) + offset;
  }
}

// Copy inter prediction info from the biggest CU
Void TComDataCU::copyInterPredInfoFrom    ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefPicList )
{
  m_pcPic              = pcCU->getPic();
  m_pcSlice            = pcCU->getSlice();
  m_ctuRsAddr          = pcCU->getCtuRsAddr();
  m_absZIdxInCtu       = uiAbsPartIdx;
#if JVET_C0024_QTBT
  ChannelType eCType  = getTextType();
  assert(isLuma(eCType));
#endif
#if BARC
  m_BARCFlag[0] = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
  m_BARCFlag[1] = pcCU->getBARCFlag(CHANNEL_TYPE_CHROMA);
#endif


  Int iRastPartIdx     = g_auiZscanToRaster[uiAbsPartIdx];
  m_uiCUPelX           = pcCU->getCUPelX() + m_pcPic->getMinCUWidth ()*( iRastPartIdx % m_pcPic->getNumPartInCtuWidth() );
  m_uiCUPelY           = pcCU->getCUPelY() + m_pcPic->getMinCUHeight()*( iRastPartIdx / m_pcPic->getNumPartInCtuWidth() );

  m_pCtuAboveLeft      = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight     = pcCU->getCtuAboveRight();
  m_pCtuAbove          = pcCU->getCtuAbove();
  m_pCtuLeft           = pcCU->getCtuLeft();

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_apcCUColocated[i]  = pcCU->getCUColocated(RefPicList(i));
  }

  m_skipFlag           = pcCU->getSkipFlag ()             + uiAbsPartIdx;
#if VCEG_AZ05_INTRA_MPI
  m_MPIIdx             = pcCU->getMPIIdx ()               + uiAbsPartIdx;
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  m_PDPCIdx            = pcCU->getPDPCIdx()               + uiAbsPartIdx;
#endif

#if DIMD_INTRA_PRED
  for (UInt ch = 0; ch < MAX_NUM_CHANNEL_TYPE; ch++)
  {
    m_ucDIMDEnabledFlag[ch] = pcCU->getDIMDEnabledFlag(ChannelType(ch)) + uiAbsPartIdx;
#if DIMD_NUM_INTRA_DIR_INC
    m_ucExtIntraDir[ch] = pcCU->getExtIntraDir(ChannelType(ch)) + uiAbsPartIdx;
#endif
    m_ucDIMDNoBTLevelFlag[ch] = pcCU->getDIMDNoBTLevelFlag(ChannelType(ch)) + uiAbsPartIdx;
  }
  m_puhSecIntraDir = pcCU->getSecIntraDir() + uiAbsPartIdx;
#endif

#if RRC_DWDIP
  m_DWDIPIdx           = pcCU->getDWDIPIdx()              + uiAbsPartIdx;
#endif  
#if RRC_DWDIP_MLR
  m_DWDIP_GPIdx        = pcCU->getDWDIP_GPIdx()           + uiAbsPartIdx;
#endif

#if RRC_INTRA_SHARPENING
  m_sharpFlag          = pcCU->getIntraSharpFlag()        + uiAbsPartIdx;
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
  for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
  {
    m_ROTIdx[i]          = pcCU->getROTIdx (ChannelType(i)) + uiAbsPartIdx;
  }
#else
  m_ROTIdx             = pcCU->getROTIdx ()               + uiAbsPartIdx;
#endif
#endif
#if !JVET_C0024_QTBT
  m_pePartSize         = pcCU->getPartitionSize ()        + uiAbsPartIdx;
#endif
  m_pePredMode         = pcCU->getPredictionMode()        + uiAbsPartIdx;
  m_ChromaQpAdj        = pcCU->getChromaQpAdj()           + uiAbsPartIdx;
  m_CUTransquantBypass = pcCU->getCUTransquantBypass()    + uiAbsPartIdx;
  m_puhInterDir        = pcCU->getInterDir      ()        + uiAbsPartIdx;

#if JVET_C0024_QTBT
  m_puhDepth[eCType]           = pcCU->getDepth ()                + uiAbsPartIdx;
  m_puhWidth[eCType]           = pcCU->getWidth ()                + uiAbsPartIdx;
  m_puhHeight[eCType]          = pcCU->getHeight()                + uiAbsPartIdx;
  m_puhBTSplitMode[eCType][0]     = pcCU->getBTSplitModePart(0)   + uiAbsPartIdx;
  m_puhBTSplitMode[eCType][1]     = pcCU->getBTSplitModePart(1)   + uiAbsPartIdx;
#if BT_FORCE_BF
  m_puhBTSplitMode[eCType][2] = pcCU->getBTSplitModePart(2) + uiAbsPartIdx;
  m_puhBTSplitMode[eCType][3] = pcCU->getBTSplitModePart(3) + uiAbsPartIdx;
#endif
#if TT_TREE_AFFINE
  m_puhCUabsZIdxInCtu[eCType][0]  = pcCU->getCUAbsZIdxInCtu(0)    + uiAbsPartIdx;
  m_puhCUabsZIdxInCtu[eCType][1]  = pcCU->getCUAbsZIdxInCtu(1)    + uiAbsPartIdx;
#endif
#else
  m_puhDepth           = pcCU->getDepth ()                + uiAbsPartIdx;
  m_puhWidth           = pcCU->getWidth ()                + uiAbsPartIdx;
  m_puhHeight          = pcCU->getHeight()                + uiAbsPartIdx;
#endif

  m_pbMergeFlag        = pcCU->getMergeFlag()             + uiAbsPartIdx;
  m_puhMergeIndex      = pcCU->getMergeIndex()            + uiAbsPartIdx;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  m_peMergeType        = pcCU->getMergeType()             + uiAbsPartIdx;
#endif
#if COM16_C806_OBMC
  m_OBMCFlag           = pcCU->getOBMCFlag ()             + uiAbsPartIdx;
#endif
#if VCEG_AZ07_FRUC_MERGE
  m_puhFRUCMgrMode     = pcCU->getFRUCMgrMode()           + uiAbsPartIdx;
#endif
#if VCEG_AZ07_IMV
  m_iMVFlag            = pcCU->getiMVFlag()               + uiAbsPartIdx;
  m_piMVCandNum        = pcCU->getiMVCandNum()            + uiAbsPartIdx;
#endif
#if VCEG_AZ06_IC
  m_pbICFlag           = pcCU->getICFlag()                + uiAbsPartIdx;
#endif
#if COM16_C1016_AFFINE
  m_affineFlag         = pcCU->getAffineFlag()            + uiAbsPartIdx;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM // maybe bug fix
  m_affineParamFlag    = pcCU->getAffineParamFlag() + uiAbsPartIdx;
#endif

#if MRG_CONSIDER_COMPLEX
  m_pbMergeComplexFlag   = pcCU->getMergeComplexFlag()    + uiAbsPartIdx;
  m_pNodeMotionInfo      = pcCU->getNodeMotionInfo()      + uiAbsPartIdx;
  m_pNumMergeComplexCand = pcCU->getNumMergeComplexCand() + uiAbsPartIdx;
  m_puhMergeComplexIndex = pcCU->getMergeComplexIndex()   + uiAbsPartIdx;
#endif

#if RRC_SELECTABLE_IF
  m_pSIFFlag           = pcCU->getSIFFlag()               + uiAbsPartIdx;
#endif
  m_apiMVPIdx[eRefPicList] = pcCU->getMVPIdx(eRefPicList) + uiAbsPartIdx;
  m_apiMVPNum[eRefPicList] = pcCU->getMVPNum(eRefPicList) + uiAbsPartIdx;
  m_acCUMvField[ eRefPicList ].linkToWithOffset( pcCU->getCUMvField(eRefPicList), uiAbsPartIdx );
#if TT_TREE_DS3 //save address of CTU MV field
  m_acCUMvFieldCtu[ eRefPicList ].linkToWithOffset( pcCU->getCUMvField(eRefPicList), 0 );
  m_bUseAbsoluteZIdxToCtu = true;
  m_pCtuCurr = pcCU;
#endif
}

#if VCEG_AZ08_INTER_KLT
// Copy  CU data to this CU.
// One of quarter parts overwritten by predicted sub part.
Void TComDataCU::copySameSizeCUFrom(TComDataCU* pcCU, UInt uiPartUnitIdx, UInt uiDepth)
{
    assert(uiPartUnitIdx<4);

    m_dTotalCost = 0; 
    m_uiTotalDistortion = 0; 
    m_uiTotalBits = 0; 
    m_uiTotalBins = 0;
    UInt uiOffset = 0;

    const UInt numValidComp = pcCU->getPic()->getNumberValidComponents();
    const UInt numValidChan = pcCU->getPic()->getChromaFormat() == CHROMA_400 ? 1 : 2;

    UInt uiNumPartition = pcCU->getTotalNumPart();
    Int iSizeInUchar = sizeof(UChar)* uiNumPartition;
    Int iSizeInBool = sizeof(Bool)* uiNumPartition;
#if COM16_C806_LARGE_CTU && !JVET_C0024_QTBT
    Int iSizeInUshort = sizeof(UShort)* uiNumPartition;
#endif
#if BARC
    m_BARCFlag[0] = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
    m_BARCFlag[1] = pcCU->getBARCFlag(CHANNEL_TYPE_CHROMA);
#endif
#if BARC_UP_SAMPLE_RDO
	m_upCNNFlag = pcCU->getUpCNNFlag();
#endif
    Int sizeInChar = sizeof(Char)* uiNumPartition;
    memcpy(m_skipFlag + uiOffset, pcCU->getSkipFlag(), sizeof(*m_skipFlag)   * uiNumPartition);
#if VCEG_AZ05_INTRA_MPI
    memcpy(m_MPIIdx + uiOffset, pcCU->getMPIIdx(), sizeof(*m_MPIIdx)   * uiNumPartition);
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    memcpy(m_PDPCIdx + uiOffset, pcCU->getPDPCIdx(), sizeof(*m_PDPCIdx)   * uiNumPartition);
#endif
#if RRC_DWDIP
    memcpy(m_DWDIPIdx + uiOffset, pcCU->getDWDIPIdx(), sizeof(*m_DWDIPIdx) * uiNumPartition);
#endif
#if RRC_DWDIP_MLR
    memcpy(m_DWDIP_GPIdx + uiOffset, pcCU->getDWDIP_GPIdx(), sizeof(*m_DWDIP_GPIdx) * uiNumPartition);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for( UInt i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ ) 
    {
      memcpy(m_ROTIdx[i] + uiOffset, pcCU->getROTIdx(ChannelType(i)), sizeof(*m_ROTIdx[i])   * uiNumPartition);
    }
#else
    memcpy(m_ROTIdx + uiOffset, pcCU->getROTIdx(), sizeof(*m_ROTIdx)   * uiNumPartition);
#endif
#endif
#if JVET_C0024_DELTA_QP_FIX
    memcpy(m_phQP[getTextType()] + uiOffset, pcCU->getQP(), sizeInChar);
#else
    memcpy(m_phQP + uiOffset, pcCU->getQP(), sizeInChar);
#endif
#if !JVET_C0024_QTBT
    memcpy(m_pePartSize + uiOffset, pcCU->getPartitionSize(), sizeof(*m_pePartSize) * uiNumPartition);
#endif
    memcpy(m_pePredMode + uiOffset, pcCU->getPredictionMode(), sizeof(*m_pePredMode) * uiNumPartition);
    memcpy(m_ChromaQpAdj + uiOffset, pcCU->getChromaQpAdj(), sizeof(*m_ChromaQpAdj) * uiNumPartition);
    memcpy(m_CUTransquantBypass + uiOffset, pcCU->getCUTransquantBypass(), sizeof(*m_CUTransquantBypass) * uiNumPartition);
    memcpy(m_pbMergeFlag + uiOffset, pcCU->getMergeFlag(), iSizeInBool);
    memcpy(m_puhMergeIndex + uiOffset, pcCU->getMergeIndex(), iSizeInUchar);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    memcpy(m_peMergeType + uiOffset, pcCU->getMergeType(), iSizeInUchar);
#endif
#if COM16_C806_OBMC
    memcpy(m_OBMCFlag + uiOffset, pcCU->getOBMCFlag(), sizeof(*m_OBMCFlag)   * uiNumPartition);
#endif
#if VCEG_AZ07_FRUC_MERGE
    memcpy(m_puhFRUCMgrMode + uiOffset, pcCU->getFRUCMgrMode(), iSizeInUchar);
#endif
#if VCEG_AZ07_IMV
    memcpy(m_iMVFlag + uiOffset, pcCU->getiMVFlag(), sizeof(*m_iMVFlag)    * uiNumPartition);
    memcpy(m_piMVCandNum + uiOffset, pcCU->getiMVCandNum(), sizeof(*m_piMVCandNum)* uiNumPartition);
#endif
#if VCEG_AZ06_IC
    memcpy(m_pbICFlag + uiOffset, pcCU->getICFlag(), iSizeInBool);
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    memcpy(m_puhAlfCtrlFlag + uiOffset, pcCU->getAlfCtrlFlag(), sizeof(*m_puhAlfCtrlFlag) * uiNumPartition);
#else
    memcpy(m_puiAlfCtrlFlag + uiOffset, pcCU->getAlfCtrlFlag(), sizeof(*m_puiAlfCtrlFlag) * uiNumPartition);
#endif
#endif

    for (UInt ch = 0; ch<numValidChan; ch++)
    {
        memcpy(m_puhIntraDir[ch] + uiOffset, pcCU->getIntraDir(ChannelType(ch)), iSizeInUchar);
#if MLIP_EXTEND_REF_LINE
        memcpy(m_puhExternRef[ch] + uiOffset, pcCU->getExternRef(ChannelType(ch)), iSizeInUchar);
#endif
    }

    memcpy(m_puhInterDir + uiOffset, pcCU->getInterDir(), iSizeInUchar);
#if !JVET_C0024_QTBT
    memcpy(m_puhTrIdx + uiOffset, pcCU->getTransformIdx(), iSizeInUchar);
#endif

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < numValidChan; ch++)
    {
      memcpy(m_ucDIMDEnabledFlag[ch] + uiOffset, pcCU->getDIMDEnabledFlag(ChannelType(ch)), iSizeInUchar);
#if DIMD_NUM_INTRA_DIR_INC
      memcpy(m_ucExtIntraDir[ch] + uiOffset, pcCU->getExtIntraDir(ChannelType(ch)), iSizeInUchar);
#endif
      memcpy(m_ucDIMDNoBTLevelFlag[ch] + uiOffset, pcCU->getDIMDNoBTLevelFlag(ChannelType(ch)), iSizeInUchar);
    }
    memcpy(m_puhSecIntraDir + uiOffset, pcCU->getSecIntraDir(), iSizeInUchar);
#endif

#if COM16_C806_EMT
    memcpy(m_puhEmtTuIdx + uiOffset, pcCU->getEmtTuIdx(), iSizeInUchar);
    memcpy(m_puhEmtCuFlag + uiOffset, pcCU->getEmtCuFlag(), iSizeInUchar);
#endif

#if COM16_C1016_AFFINE
    memcpy(m_affineFlag + uiOffset, pcCU->getAffineFlag(), iSizeInBool);
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    memcpy(m_affineParamFlag + uiOffset, pcCU->getAffineParamFlag(), iSizeInBool);
#endif

#if MRG_CONSIDER_COMPLEX
    memcpy( m_pbMergeComplexFlag   + uiOffset, pcCU->getMergeComplexFlag(),        iSizeInBool  );
    memcpy( m_pNodeMotionInfo      + uiOffset, pcCU->getNodeMotionInfo(),  sizeof( *m_pNodeMotionInfo ) * uiNumPartition );
    memcpy( m_pNumMergeComplexCand + uiOffset, pcCU->getNumMergeComplexCand(), iSizeInUchar );
    memcpy( m_puhMergeComplexIndex + uiOffset, pcCU->getMergeComplexIndex(),  iSizeInUchar );
#endif

#if RRC_INTRA_SHARPENING
    memcpy(m_sharpFlag + uiOffset, pcCU->getIntraSharpFlag(), sizeof(*m_sharpFlag) * uiNumPartition);
#endif

#if RRC_SELECTABLE_IF
    memcpy(m_pSIFFlag + uiOffset, pcCU->getSIFFlag(), sizeof(*m_pSIFFlag) * uiNumPartition);
#endif
    for (UInt comp = 0; comp<numValidComp; comp++)
    {
        memcpy(m_crossComponentPredictionAlpha[comp] + uiOffset, pcCU->getCrossComponentPredictionAlpha(ComponentID(comp)), iSizeInUchar);
        memcpy(m_puhTransformSkip[comp] + uiOffset, pcCU->getTransformSkip(ComponentID(comp)), iSizeInUchar);
        memcpy(m_puhCbf[comp] + uiOffset, pcCU->getCbf(ComponentID(comp)), iSizeInUchar);
        memcpy(m_explicitRdpcmMode[comp] + uiOffset, pcCU->getExplicitRdpcmMode(ComponentID(comp)), iSizeInUchar);
#if VCEG_AZ08_KLT_COMMON
        memcpy(m_puhKLTFlag[comp] + uiOffset, pcCU->getKLTFlag(ComponentID(comp)), iSizeInUchar);
#endif
    }

#if JVET_C0024_QTBT    
    memcpy(m_puhDepth[getTextType()] + uiOffset, pcCU->getDepth(), iSizeInUchar);
    memcpy(m_puhWidth[getTextType()] + uiOffset, pcCU->getWidth(), iSizeInUchar);
    memcpy(m_puhHeight[getTextType()] + uiOffset, pcCU->getHeight(), iSizeInUchar);
#if TT_TREE_AFFINE
    memcpy(m_puhCUabsZIdxInCtu[getTextType()][0] + uiOffset, pcCU->getCUAbsZIdxInCtu(0), iSizeInUchar);
    memcpy(m_puhCUabsZIdxInCtu[getTextType()][1] + uiOffset, pcCU->getCUAbsZIdxInCtu(1), iSizeInUchar);
#endif
    memcpy(m_puhBTSplitMode[getTextType()][0] + uiOffset, pcCU->getBTSplitModePart(0), iSizeInUchar);
    memcpy(m_puhBTSplitMode[getTextType()][1] + uiOffset, pcCU->getBTSplitModePart(1), iSizeInUchar);
#if BT_FORCE_BF
    memcpy(m_puhBTSplitMode[getTextType()][2] + uiOffset, pcCU->getBTSplitModePart(2), iSizeInUchar);
    memcpy(m_puhBTSplitMode[getTextType()][3] + uiOffset, pcCU->getBTSplitModePart(3), iSizeInUchar);
#endif
#else
    memcpy(m_puhDepth + uiOffset, pcCU->getDepth(), iSizeInUchar);
#if COM16_C806_LARGE_CTU
    memcpy(m_puhWidth + uiOffset, pcCU->getWidth(), iSizeInUshort);
    memcpy(m_puhHeight + uiOffset, pcCU->getHeight(), iSizeInUshort);
#else
    memcpy(m_puhWidth + uiOffset, pcCU->getWidth(), iSizeInUchar);
    memcpy(m_puhHeight + uiOffset, pcCU->getHeight(), iSizeInUchar);
#endif
#endif

    memcpy(m_pbIPCMFlag + uiOffset, pcCU->getIPCMFlag(), iSizeInBool);

    m_pCtuAboveLeft = pcCU->getCtuAboveLeft();
    m_pCtuAboveRight = pcCU->getCtuAboveRight();
    m_pCtuAbove = pcCU->getCtuAbove();
    m_pCtuLeft = pcCU->getCtuLeft();

    for (UInt i = 0; i<NUM_REF_PIC_LIST_01; i++)
    {
        const RefPicList rpl = RefPicList(i);
        memcpy(m_apiMVPIdx[rpl] + uiOffset, pcCU->getMVPIdx(rpl), iSizeInUchar);
        memcpy(m_apiMVPNum[rpl] + uiOffset, pcCU->getMVPNum(rpl), iSizeInUchar);
        m_apcCUColocated[rpl] = pcCU->getCUColocated(rpl);
}

    for (UInt i = 0; i<NUM_REF_PIC_LIST_01; i++)
    {
        const RefPicList rpl = RefPicList(i);
#if JVET_C0024_QTBT
        m_acCUMvField[rpl].copyFromTo( pcCU->getCUMvField(rpl), uiNumPartition, 0, uiOffset );
#else
        m_acCUMvField[rpl].copyFrom(pcCU->getCUMvField(rpl), pcCU->getTotalNumPart(), uiOffset);
#endif
    }

#if JVET_C0024_QTBT
    const UInt numCoeffY = (pcCU->getSlice()->getSPS()->getCTUSize()*pcCU->getSlice()->getSPS()->getCTUSize()) >> (uiDepth << 1);
#else
    const UInt numCoeffY = (pcCU->getSlice()->getSPS()->getMaxCUWidth()*pcCU->getSlice()->getSPS()->getMaxCUHeight()) >> (uiDepth << 1);
#endif
    const UInt offsetY = uiPartUnitIdx*numCoeffY;

    for (UInt ch = 0; ch<numValidComp; ch++)
    {
        const ComponentID component = ComponentID(ch);
        const UInt componentShift = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
        const UInt offset = offsetY >> componentShift;
        memcpy(m_pcTrCoeff[ch] + offset, pcCU->getCoeff(component), sizeof(TCoeff)*(numCoeffY >> componentShift));
#if ADAPTIVE_QP_SELECTION
        memcpy(m_pcArlCoeff[ch] + offset, pcCU->getArlCoeff(component), sizeof(TCoeff)*(numCoeffY >> componentShift));
#endif
        memcpy(m_pcIPCMSample[ch] + offset, pcCU->getPCMSample(component), sizeof(Pel)*(numCoeffY >> componentShift));
    }

    m_uiTotalBins += pcCU->getTotalBins();
}

#endif

// Copy small CU to bigger CU.
// One of quarter parts overwritten by predicted sub part.
#if JVET_C0024_QTBT 
Void TComDataCU::copyPartFrom( TComDataCU* pcCU, UInt uiPartUnitIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight
#if TT_TREE_DS1
  , Bool bIsTT
#endif
  )
#else
Void TComDataCU::copyPartFrom( TComDataCU* pcCU, UInt uiPartUnitIdx, UInt uiDepth )
#endif
{
  assert( uiPartUnitIdx<4 );

  m_dTotalCost         += pcCU->getTotalCost();
  m_uiTotalDistortion  += pcCU->getTotalDistortion();
  m_uiTotalBits        += pcCU->getTotalBits();
#if JVET_C0024_QTBT
  const UInt numValidComp=pcCU->getPic()->getNumberValidComponents();
  const UInt numValidChan=pcCU->getPic()->getChromaFormat()==CHROMA_400 ? 1:2;

  UInt uiCTUSize = pcCU->getSlice()->getSPS()->getCTUSize();
  ChannelType eCType  = getTextType();
#if BARC
  m_BARCFlag[0] = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
  m_BARCFlag[1] = pcCU->getBARCFlag(CHANNEL_TYPE_CHROMA);
#endif
#if BARC_UP_SAMPLE_RDO
  m_upCNNFlag = pcCU->getUpCNNFlag();
#endif
  UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
  UInt uiShort, uiLong;
  UInt uiStride;

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiWidthSave  = uiWidth;
  UInt uiHeightSave = uiHeight;
  UInt uiNumXDiv = (uiWidth  > m_pcPic->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiHeight > m_pcPic->getMinCUHeight())? 2 : 1;
  uiWidth  = (uiNumXDiv==1)? uiWidth  : (uiWidth>>1);
  uiHeight = (uiNumYDiv==1)? uiHeight : (uiHeight>>1);
  UInt uiLenX    = uiWidth  / m_pcPic->getMinCUWidth();
  UInt uiLenY    = uiHeight / m_pcPic->getMinCUHeight();
  UInt uiRaster2;
  UInt uiRasterSrc = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
  UInt uiRasterDst = g_auiZscanToRaster[m_absZIdxInCtu];
  assert(uiRasterSrc>=uiRasterDst);
  UInt uiDeltaY    = (uiRasterSrc - uiRasterDst) / m_pcPic->getNumPartInCtuWidth();
  UInt uiDeltaX    = (uiRasterSrc - uiRasterDst) % m_pcPic->getNumPartInCtuWidth();

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      uiRaster = uiDeltaX + (uiLenX * xIdx) + (uiLenY * yIdx + uiDeltaY) * m_pcPic->getNumPartInCtuWidth();
      uiRaster2 = (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
#endif
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = m_pcPic->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  uiDepth = g_aucConvertToBit[uiCTUSize] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;

  UInt uiZorderDst, uiZorderSrc;
  for (UInt i=0; i<uiLong; i+=uiShort)
  {
#if TT_TREE_DS1
    uiZorderSrc = g_auiRasterToZscan[uiRaster2];
    uiZorderDst = g_auiRasterToZscan[uiRaster];
#else
    uiZorderSrc = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
    uiZorderDst = g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu;
#endif

    memcpy(m_skipFlag + uiZorderDst, pcCU->getSkipFlag()+uiZorderSrc, uiCurrPartNumb );
#if VCEG_AZ05_INTRA_MPI
    memcpy( m_MPIIdx     + uiZorderDst, pcCU->getMPIIdx()+uiZorderSrc, uiCurrPartNumb );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
    memcpy( m_PDPCIdx    + uiZorderDst, pcCU->getPDPCIdx()+uiZorderSrc, uiCurrPartNumb );
#endif
#if RRC_DWDIP
    memcpy( m_DWDIPIdx   + uiZorderDst, pcCU->getDWDIPIdx()+uiZorderSrc, uiCurrPartNumb);
#endif

#if RRC_DWDIP_MLR
    memcpy( m_DWDIP_GPIdx   + uiZorderDst, pcCU->getDWDIP_GPIdx()+uiZorderSrc, uiCurrPartNumb);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
    for (UInt ch=0; ch<numValidChan; ch++)
    {
      memcpy( m_ROTIdx[ch]   + uiZorderDst, pcCU->getROTIdx(ChannelType(ch))+uiZorderSrc, uiCurrPartNumb );
    }
#else
    memcpy( m_ROTIdx     + uiZorderDst, pcCU->getROTIdx()+uiZorderSrc, uiCurrPartNumb );
#endif
#endif
#if JVET_C0024_DELTA_QP_FIX
    memcpy( m_phQP[eCType] + uiZorderDst, pcCU->getQP()+uiZorderSrc, uiCurrPartNumb );
#else
    memcpy( m_phQP       + uiZorderDst, pcCU->getQP()+uiZorderSrc, uiCurrPartNumb );
#endif

    memcpy( m_pePredMode + uiZorderDst, pcCU->getPredictionMode()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_ChromaQpAdj + uiZorderDst, pcCU->getChromaQpAdj()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_CUTransquantBypass + uiZorderDst, pcCU->getCUTransquantBypass()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_pbMergeFlag         + uiZorderDst, pcCU->getMergeFlag()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_puhMergeIndex       + uiZorderDst, pcCU->getMergeIndex()+uiZorderSrc, uiCurrPartNumb );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    memcpy( m_peMergeType         + uiZorderDst, pcCU->getMergeType()+uiZorderSrc, uiCurrPartNumb );
#endif
#if COM16_C806_OBMC
    memcpy( m_OBMCFlag + uiZorderDst, pcCU->getOBMCFlag()+uiZorderSrc, uiCurrPartNumb );
#endif
#if VCEG_AZ07_FRUC_MERGE
    memcpy( m_puhFRUCMgrMode + uiZorderDst, pcCU->getFRUCMgrMode()+uiZorderSrc, uiCurrPartNumb );
#endif
#if VCEG_AZ07_IMV
    memcpy( m_iMVFlag    + uiZorderDst, pcCU->getiMVFlag()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_piMVCandNum + uiZorderDst, pcCU->getiMVCandNum()+uiZorderSrc, uiCurrPartNumb );
#endif
#if VCEG_AZ06_IC
    memcpy( m_pbICFlag + uiZorderDst, pcCU->getICFlag()+uiZorderSrc, uiCurrPartNumb );
#endif
#if ALF_HM3_REFACTOR
#if JVET_C0024_QTBT
    memcpy( m_puhAlfCtrlFlag      + uiZorderDst, pcCU->getAlfCtrlFlag()+uiZorderSrc, uiCurrPartNumb );
#else
    memcpy( m_puiAlfCtrlFlag      + uiZorderDst, pcCU->getAlfCtrlFlag()+uiZorderSrc, uiCurrPartNumb );
#endif
#endif
    for (UInt ch=0; ch<numValidChan; ch++)
    {
      memcpy( m_puhIntraDir[ch]   + uiZorderDst, pcCU->getIntraDir(ChannelType(ch))+uiZorderSrc, uiCurrPartNumb );
#if MLIP_EXTEND_REF_LINE
      memcpy(m_puhExternRef[ch] + uiZorderDst, pcCU->getExternRef(ChannelType(ch)) + uiZorderSrc, uiCurrPartNumb);
#endif
    }

    memcpy( m_puhInterDir         + uiZorderDst, pcCU->getInterDir()+uiZorderSrc, uiCurrPartNumb );

#if DIMD_INTRA_PRED
    for (UInt ch = 0; ch < numValidChan; ch++)
    {
      memcpy(m_ucDIMDEnabledFlag[ch] + uiZorderDst, pcCU->getDIMDEnabledFlag(ChannelType(ch)) + uiZorderSrc, uiCurrPartNumb);
#if DIMD_NUM_INTRA_DIR_INC
      memcpy(m_ucExtIntraDir[ch] + uiZorderDst, pcCU->getExtIntraDir(ChannelType(ch)) + uiZorderSrc, uiCurrPartNumb);
#endif
      memcpy(m_ucDIMDNoBTLevelFlag[ch] + uiZorderDst, pcCU->getDIMDNoBTLevelFlag(ChannelType(ch)) + uiZorderSrc, uiCurrPartNumb);
    }
    memcpy(m_puhSecIntraDir + uiZorderDst, pcCU->getSecIntraDir() + uiZorderSrc, uiCurrPartNumb);
#endif

#if COM16_C806_EMT
    memcpy( m_puhEmtTuIdx         + uiZorderDst, pcCU->getEmtTuIdx()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_puhEmtCuFlag        + uiZorderDst, pcCU->getEmtCuFlag()+uiZorderSrc, uiCurrPartNumb );
#endif

#if COM16_C1016_AFFINE
    memcpy( m_affineFlag          + uiZorderDst, pcCU->getAffineFlag()+uiZorderSrc, uiCurrPartNumb );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
    memcpy( m_affineParamFlag     + uiZorderDst, pcCU->getAffineParamFlag()+uiZorderSrc, uiCurrPartNumb );
#endif

#if MRG_CONSIDER_COMPLEX
    memcpy( m_pbMergeComplexFlag   + uiZorderDst, pcCU->getMergeComplexFlag()+uiZorderSrc,    uiCurrPartNumb  );
    memcpy( m_pNodeMotionInfo      + uiZorderDst, pcCU->getNodeMotionInfo()+uiZorderSrc,      sizeof( *m_pNodeMotionInfo ) * uiCurrPartNumb );
    memcpy( m_pNumMergeComplexCand + uiZorderDst, pcCU->getNumMergeComplexCand()+uiZorderSrc, sizeof( *m_pNumMergeComplexCand ) * uiCurrPartNumb );
    memcpy( m_puhMergeComplexIndex + uiZorderDst, pcCU->getMergeComplexIndex()+uiZorderSrc,   uiCurrPartNumb );
#endif

#if RRC_SELECTABLE_IF
    memcpy( m_pSIFFlag            + uiZorderDst, pcCU->getSIFFlag()+uiZorderSrc, uiCurrPartNumb );
#endif
    
#if RRC_INTRA_SHARPENING
    memcpy(m_sharpFlag            + uiZorderDst, pcCU->getIntraSharpFlag()+uiZorderSrc, sizeof(*m_sharpFlag)*uiCurrPartNumb);
#endif

    for(UInt comp=0; comp<numValidComp; comp++)
    {
      memcpy( m_crossComponentPredictionAlpha[comp] + uiZorderDst, pcCU->getCrossComponentPredictionAlpha(ComponentID(comp))+uiZorderSrc, uiCurrPartNumb );
      memcpy( m_puhTransformSkip[comp]              + uiZorderDst, pcCU->getTransformSkip(ComponentID(comp))+uiZorderSrc, uiCurrPartNumb );
      memcpy( m_puhCbf[comp]                        + uiZorderDst, pcCU->getCbf(ComponentID(comp))+uiZorderSrc, uiCurrPartNumb );
      memcpy( m_explicitRdpcmMode[comp]             + uiZorderDst, pcCU->getExplicitRdpcmMode(ComponentID(comp))+uiZorderSrc, uiCurrPartNumb );
#if VCEG_AZ08_KLT_COMMON
      memcpy( m_puhKLTFlag[comp]                    + uiZorderDst, pcCU->getKLTFlag(ComponentID(comp))+uiZorderSrc, uiCurrPartNumb );
#endif
    }

    memcpy( m_puhDepth[eCType]  + uiZorderDst, pcCU->getDepth()+uiZorderSrc, uiCurrPartNumb );
    memcpy( m_puhWidth[eCType]  + uiZorderDst, pcCU->getWidth()+uiZorderSrc, uiCurrPartNumb  );
    memcpy( m_puhHeight[eCType] + uiZorderDst, pcCU->getHeight()+uiZorderSrc, uiCurrPartNumb  );
#if TT_TREE_AFFINE
    memcpy( m_puhCUabsZIdxInCtu[eCType][0] + uiZorderDst, pcCU->getCUAbsZIdxInCtu(0) + uiZorderSrc, uiCurrPartNumb);
    memcpy( m_puhCUabsZIdxInCtu[eCType][1] + uiZorderDst, pcCU->getCUAbsZIdxInCtu(1) + uiZorderSrc, uiCurrPartNumb);
#endif
    memcpy( m_puhBTSplitMode[eCType][0] + uiZorderDst, pcCU->getBTSplitModePart(0)+uiZorderSrc, uiCurrPartNumb);
    memcpy( m_puhBTSplitMode[eCType][1] + uiZorderDst, pcCU->getBTSplitModePart(1)+uiZorderSrc, uiCurrPartNumb);
#if BT_FORCE_BF
    memcpy(m_puhBTSplitMode[eCType][2] + uiZorderDst, pcCU->getBTSplitModePart(2) + uiZorderSrc, uiCurrPartNumb);
    memcpy(m_puhBTSplitMode[eCType][3] + uiZorderDst, pcCU->getBTSplitModePart(3) + uiZorderSrc, uiCurrPartNumb);
#endif

    memcpy( m_pbIPCMFlag + uiZorderDst, pcCU->getIPCMFlag()+uiZorderSrc, uiCurrPartNumb  );

    for(UInt refListIdx=0; refListIdx<NUM_REF_PIC_LIST_01; refListIdx++)
    {
      const RefPicList rpl=RefPicList(refListIdx);
      memcpy( m_apiMVPIdx[rpl] + uiZorderDst, pcCU->getMVPIdx(rpl)+uiZorderSrc, uiCurrPartNumb  );
      memcpy( m_apiMVPNum[rpl] + uiZorderDst, pcCU->getMVPNum(rpl)+uiZorderSrc, uiCurrPartNumb  );
      m_apcCUColocated[rpl] = pcCU->getCUColocated(rpl);
    }

    for(UInt refListIdx=0; refListIdx<NUM_REF_PIC_LIST_01; refListIdx++)
    {
      const RefPicList rpl=RefPicList(refListIdx);
      m_acCUMvField[rpl].copyFromTo( pcCU->getCUMvField( rpl ), uiCurrPartNumb, uiZorderSrc, uiZorderDst );
    }
    uiRaster += uiNumPartInShort * uiStride;
#if TT_TREE_DS1
    uiRaster2 += uiNumPartInShort * uiStride;
#endif
  }
#if TT_TREE_DS1 //data copy end
  }}
  //recover uiWidth & uiHeight
  uiWidth = uiWidthSave;
  uiHeight = uiHeightSave;
#endif
  m_pCtuAboveLeft      = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight     = pcCU->getCtuAboveRight();
  m_pCtuAbove          = pcCU->getCtuAbove();
  m_pCtuLeft           = pcCU->getCtuLeft();

#else
  UInt uiOffset         = pcCU->getTotalNumPart()*uiPartUnitIdx;
  const UInt numValidComp=pcCU->getPic()->getNumberValidComponents();
  const UInt numValidChan=pcCU->getPic()->getChromaFormat()==CHROMA_400 ? 1:2;

  UInt uiNumPartition = pcCU->getTotalNumPart();
  Int iSizeInUchar  = sizeof( UChar ) * uiNumPartition;
  Int iSizeInBool   = sizeof( Bool  ) * uiNumPartition;
#if COM16_C806_LARGE_CTU
  Int iSizeInUshort = sizeof( UShort ) * uiNumPartition;
#endif

  Int sizeInChar  = sizeof( Char ) * uiNumPartition;
  memcpy( m_skipFlag   + uiOffset, pcCU->getSkipFlag(),       sizeof( *m_skipFlag )   * uiNumPartition );
#if VCEG_AZ05_INTRA_MPI
  memcpy( m_MPIIdx     + uiOffset, pcCU->getMPIIdx(),         sizeof( *m_MPIIdx )   * uiNumPartition );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  memcpy( m_PDPCIdx    + uiOffset, pcCU->getPDPCIdx(),        sizeof( *m_PDPCIdx)   * uiNumPartition);
#endif
#if RRC_DWDIP
  memcpy(m_DWDIPIdx    + uiOffset, pcCU->getDWDIPIdx(),       sizeof(*m_DWDIPIdx)   * uiNumPartition);
#endif

#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  memcpy( m_ROTIdx     + uiOffset, pcCU->getROTIdx(),         sizeof( *m_ROTIdx )   * uiNumPartition );
#endif
  memcpy( m_phQP       + uiOffset, pcCU->getQP(),             sizeInChar                        );
#if !JVET_C0024_QTBT
  memcpy( m_pePartSize + uiOffset, pcCU->getPartitionSize(),  sizeof( *m_pePartSize ) * uiNumPartition );
#endif
  memcpy( m_pePredMode + uiOffset, pcCU->getPredictionMode(), sizeof( *m_pePredMode ) * uiNumPartition );
  memcpy( m_ChromaQpAdj + uiOffset, pcCU->getChromaQpAdj(),   sizeof( *m_ChromaQpAdj ) * uiNumPartition );
  memcpy( m_CUTransquantBypass + uiOffset, pcCU->getCUTransquantBypass(), sizeof( *m_CUTransquantBypass ) * uiNumPartition );
  memcpy( m_pbMergeFlag         + uiOffset, pcCU->getMergeFlag(),         iSizeInBool  );
  memcpy( m_puhMergeIndex       + uiOffset, pcCU->getMergeIndex(),        iSizeInUchar );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  memcpy( m_peMergeType         + uiOffset, pcCU->getMergeType(),         iSizeInUchar );
#endif
#if COM16_C806_OBMC
  memcpy( m_OBMCFlag + uiOffset, pcCU->getOBMCFlag(),       sizeof( *m_OBMCFlag )   * uiNumPartition );
#endif
#if VCEG_AZ07_FRUC_MERGE
  memcpy( m_puhFRUCMgrMode + uiOffset, pcCU->getFRUCMgrMode(), iSizeInUchar );
#endif
#if VCEG_AZ07_IMV
  memcpy( m_iMVFlag    + uiOffset, pcCU->getiMVFlag(),      sizeof( *m_iMVFlag )    * uiNumPartition );
  memcpy( m_piMVCandNum + uiOffset, pcCU->getiMVCandNum(),  sizeof( *m_piMVCandNum )* uiNumPartition );
#endif
#if VCEG_AZ06_IC
  memcpy( m_pbICFlag + uiOffset, pcCU->getICFlag(), iSizeInBool );
#endif
#if ALF_HM3_REFACTOR
  memcpy( m_puiAlfCtrlFlag      + uiOffset, pcCU->getAlfCtrlFlag(),       sizeof(*m_puiAlfCtrlFlag) * uiNumPartition  );
#endif

  for (UInt ch=0; ch<numValidChan; ch++)
  {
    memcpy( m_puhIntraDir[ch]   + uiOffset, pcCU->getIntraDir(ChannelType(ch)), iSizeInUchar );
#if MLIP_EXTEND_REF_LINE
    memcpy(m_puhExternRef[ch] + uiOffset, pcCU->getExternRef(ChannelType(ch)), iSizeInUchar);
#endif
  }

  memcpy( m_puhInterDir         + uiOffset, pcCU->getInterDir(),          iSizeInUchar );
#if !JVET_C0024_QTBT
  memcpy( m_puhTrIdx            + uiOffset, pcCU->getTransformIdx(),      iSizeInUchar );
#endif

#if COM16_C806_EMT
  memcpy( m_puhEmtTuIdx         + uiOffset, pcCU->getEmtTuIdx(),          iSizeInUchar );
  memcpy( m_puhEmtCuFlag        + uiOffset, pcCU->getEmtCuFlag(),         iSizeInUchar );
#endif

#if COM16_C1016_AFFINE
  memcpy( m_affineFlag          + uiOffset, pcCU->getAffineFlag(),        iSizeInBool );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  memcpy( m_affineParamFlag + uiOffset, pcCU->getAffineParamFlag(), iSizeInBool );
#endif

#if MRG_CONSIDER_COMPLEX
  memcpy( m_pbMergeComplexFlag   + uiOffset, pcCU->getMergeComplexFlag(),    iSizeInBool  );
  memcpy( m_pNodeMotionInfo      + uiOffset, pcCU->getNodeMotionInfo(),      sizeof( *m_pNodeMotionInfo ) * uiNumPartition );
  memcpy( m_pNumMergeComplexCand + uiOffset, pcCU->getNumMergeComplexCand(), sizeof( *m_pNumMergeComplexCand ) * uiNumPartition );
  memcpy( m_puhMergeComplexIndex + uiOffset, pcCU->getMergeComplexIndex(),   iSizeInUchar );
#endif

#if RRC_SELECTABLE_IF
  memcpy(m_pSIFFlag             + uiOffset, pcCU->getSIFFlag(), sizeof(*m_pSIFFlag)*uiNumPartition);
#endif
  for(UInt comp=0; comp<numValidComp; comp++)
  {
    memcpy( m_crossComponentPredictionAlpha[comp] + uiOffset, pcCU->getCrossComponentPredictionAlpha(ComponentID(comp)), iSizeInUchar );
    memcpy( m_puhTransformSkip[comp]              + uiOffset, pcCU->getTransformSkip(ComponentID(comp))                , iSizeInUchar );
    memcpy( m_puhCbf[comp]                        + uiOffset, pcCU->getCbf(ComponentID(comp))                          , iSizeInUchar );
    memcpy( m_explicitRdpcmMode[comp]             + uiOffset, pcCU->getExplicitRdpcmMode(ComponentID(comp))            , iSizeInUchar );
#if VCEG_AZ08_KLT_COMMON
    memcpy( m_puhKLTFlag[comp]                    + uiOffset, pcCU->getKLTFlag(ComponentID(comp))                      , iSizeInUchar );
#endif
  }

#if JVET_C0024_QTBT
  memcpy( m_puhDepth[eCType]  + uiOffset, pcCU->getDepth(),  iSizeInUchar );
#else
  memcpy( m_puhDepth  + uiOffset, pcCU->getDepth(),  iSizeInUchar );
#endif
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
  memcpy( m_puhWidth[eCType]  + uiOffset, pcCU->getWidth(),  iSizeInUchar );
  memcpy( m_puhHeight[eCType] + uiOffset, pcCU->getHeight(), iSizeInUchar );
#else
  memcpy( m_puhWidth  + uiOffset, pcCU->getWidth(),  iSizeInUshort );
  memcpy( m_puhHeight + uiOffset, pcCU->getHeight(), iSizeInUshort );
#endif
#else
  memcpy( m_puhWidth  + uiOffset, pcCU->getWidth(),  iSizeInUchar );
  memcpy( m_puhHeight + uiOffset, pcCU->getHeight(), iSizeInUchar );
#endif

  memcpy( m_pbIPCMFlag + uiOffset, pcCU->getIPCMFlag(), iSizeInBool );

  m_pCtuAboveLeft      = pcCU->getCtuAboveLeft();
  m_pCtuAboveRight     = pcCU->getCtuAboveRight();
  m_pCtuAbove          = pcCU->getCtuAbove();
  m_pCtuLeft           = pcCU->getCtuLeft();

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    memcpy( m_apiMVPIdx[rpl] + uiOffset, pcCU->getMVPIdx(rpl), iSizeInUchar );
    memcpy( m_apiMVPNum[rpl] + uiOffset, pcCU->getMVPNum(rpl), iSizeInUchar );
    m_apcCUColocated[rpl] = pcCU->getCUColocated(rpl);
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    m_acCUMvField[rpl].copyFrom( pcCU->getCUMvField( rpl ), pcCU->getTotalNumPart(), uiOffset );
  }

#if RRC_INTRA_SHARPENING
  memcpy(m_sharpFlag + uiOffset, pcCU->getIntraSharpFlag(), sizeof(*m_sharpFlag)*uiNumPartition);
#endif

#endif  //JVET_C0024_QTBT

#if JVET_C0024_QTBT
  const UInt numCoeffY = uiWidth*uiHeight;
#else
  const UInt numCoeffY = (pcCU->getSlice()->getSPS()->getMaxCUWidth()*pcCU->getSlice()->getSPS()->getMaxCUHeight()) >> (uiDepth<<1);
#endif
#if TT_TREE_DS1 //correct coeff offset
  UInt offsetY;
  if(bIsTT)
  {
    if(uiPartUnitIdx==0)
      offsetY = 0;
    else if(uiPartUnitIdx==1)
      offsetY = numCoeffY>>1;
    else
      offsetY = numCoeffY*3;
  }
  else
    offsetY   = uiPartUnitIdx*numCoeffY;
#else
  const UInt offsetY   = uiPartUnitIdx*numCoeffY;
#endif

  for (UInt ch=0; ch<numValidComp; ch++)
  {
    const ComponentID component = ComponentID(ch);
    const UInt componentShift   = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
    const UInt offset           = offsetY>>componentShift;
    memcpy( m_pcTrCoeff [ch] + offset, pcCU->getCoeff(component),    sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
    memcpy( m_pcArlCoeff[ch] + offset, pcCU->getArlCoeff(component), sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
    memcpy( m_pcIPCMSample[ch] + offset, pcCU->getPCMSample(component), sizeof(Pel)*(numCoeffY>>componentShift) );
  }
  
  m_uiTotalBins += pcCU->getTotalBins();
}

// Copy current predicted part to a CU in picture.
// It is used to predict for next part
#if JVET_C0024_QTBT
Void TComDataCU::copyToPic( UChar uhDepth, UInt uiWidth, UInt uiHeight )
#else
Void TComDataCU::copyToPic( UChar uhDepth )
#endif
{
  TComDataCU* pCtu = m_pcPic->getCtu( m_ctuRsAddr );
  const UInt numValidComp=pCtu->getPic()->getNumberValidComponents();
  const UInt numValidChan=pCtu->getPic()->getChromaFormat()==CHROMA_400 ? 1:2;

  pCtu->getTotalCost()       = m_dTotalCost;
  pCtu->getTotalDistortion() = m_uiTotalDistortion;
  pCtu->getTotalBits()       = m_uiTotalBits;
#if JVET_C0024_QTBT
  ChannelType eCType  = getTextType();
#if BARC
  pCtu->setBARCFlag(eCType, m_BARCFlag[eCType]);
#endif
  UInt uiRaster = g_auiZscanToRaster[m_absZIdxInCtu];
  UInt uiShort, uiLong;
  UInt uiStride;

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiWidthSave  = uiWidth;
  UInt uiHeightSave = uiHeight;
  UInt uiNumXDiv = (uiWidth  > m_pcPic->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiHeight > m_pcPic->getMinCUHeight())? 2 : 1;
  uiWidth  = (uiNumXDiv==1)? uiWidth  : (uiWidth>>1);
  uiHeight = (uiNumYDiv==1)? uiHeight : (uiHeight>>1);
  UInt uiLenX    = uiWidth  / m_pcPic->getMinCUWidth();
  UInt uiLenY    = uiHeight / m_pcPic->getMinCUHeight();
  UInt uiRaster2;

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      uiRaster = g_auiZscanToRaster[m_absZIdxInCtu] + (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
      uiRaster2 = (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
#endif
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = m_pcPic->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }
  uhDepth = g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uhDepth << 1);
  UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uhDepth;

  UInt uiZorderSrc, uiZorderDst;
  if (isLuma(eCType))
  {
#if BARC_UP_SAMPLE_RDO
	  pCtu->setUpCNNFlag(m_upCNNFlag);
#endif
    for (UInt i=0; i<uiLong; i+=uiShort)
    {
#if TT_TREE_DS1
      uiZorderSrc = g_auiRasterToZscan[uiRaster2];
#else
      uiZorderSrc = g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu;
#endif
      uiZorderDst = g_auiRasterToZscan[uiRaster];

      memcpy( pCtu->getSkipFlag() + uiZorderDst, m_skipFlag + uiZorderSrc, uiCurrPartNumb );
#if VCEG_AZ05_INTRA_MPI
      memcpy( pCtu->getMPIIdx()   + uiZorderDst, m_MPIIdx + uiZorderSrc, uiCurrPartNumb );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
      memcpy( pCtu->getPDPCIdx() + uiZorderDst, m_PDPCIdx + uiZorderSrc, uiCurrPartNumb);
#endif
#if RRC_DWDIP
      memcpy( pCtu->getDWDIPIdx() + uiZorderDst, m_DWDIPIdx + uiZorderSrc, uiCurrPartNumb);
#endif
#if RRC_DWDIP_MLR
      memcpy( pCtu->getDWDIP_GPIdx() + uiZorderDst, m_DWDIP_GPIdx + uiZorderSrc, uiCurrPartNumb);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
      for (UInt ch=0; ch<numValidChan; ch++)
      {
        memcpy( pCtu->getROTIdx(ChannelType(ch)) + uiZorderDst, m_ROTIdx[ch] + uiZorderSrc, uiCurrPartNumb);
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      memcpy( pCtu->getQP(eCType) + uiZorderDst, m_phQP[eCType] + uiZorderSrc, uiCurrPartNumb  );
#else
      memcpy( pCtu->getQP() + uiZorderDst, m_phQP + uiZorderSrc, uiCurrPartNumb  );
#endif

      memcpy( pCtu->getPredictionMode() + uiZorderDst, m_pePredMode + uiZorderSrc, uiCurrPartNumb );
      memcpy( pCtu->getChromaQpAdj() + uiZorderDst, m_ChromaQpAdj + uiZorderSrc, uiCurrPartNumb );
      memcpy( pCtu->getCUTransquantBypass()+ uiZorderDst, m_CUTransquantBypass + uiZorderSrc, uiCurrPartNumb );
      memcpy( pCtu->getMergeFlag()         + uiZorderDst, m_pbMergeFlag + uiZorderSrc,         uiCurrPartNumb  );
      memcpy( pCtu->getMergeIndex()        + uiZorderDst, m_puhMergeIndex + uiZorderSrc,       uiCurrPartNumb );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
      memcpy( pCtu->getMergeType()         + uiZorderDst, m_peMergeType + uiZorderSrc,         uiCurrPartNumb );
#endif
#if COM16_C806_OBMC
      memcpy( pCtu->getOBMCFlag()          + uiZorderDst, m_OBMCFlag + uiZorderSrc, uiCurrPartNumb );
#endif
#if VCEG_AZ07_FRUC_MERGE
      memcpy( pCtu->getFRUCMgrMode()       + uiZorderDst, m_puhFRUCMgrMode + uiZorderSrc,      uiCurrPartNumb  );
#endif
#if VCEG_AZ07_IMV
      memcpy( pCtu->getiMVFlag()           + uiZorderDst, m_iMVFlag + uiZorderSrc, uiCurrPartNumb );
      memcpy( pCtu->getiMVCandNum()        + uiZorderDst, m_piMVCandNum + uiZorderSrc, uiCurrPartNumb );
#endif
#if VCEG_AZ06_IC
      memcpy( pCtu->getICFlag()            + uiZorderDst, m_pbICFlag + uiZorderSrc, uiCurrPartNumb );
#endif
#if ALF_HM3_REFACTOR
      memcpy( pCtu->getAlfCtrlFlag()       + uiZorderDst, m_puhAlfCtrlFlag + uiZorderSrc,      sizeof( *m_puhAlfCtrlFlag ) * uiCurrPartNumb );
#endif
      for (UInt ch=0; ch<numValidChan; ch++)
      {
        memcpy( pCtu->getIntraDir(ChannelType(ch)) + uiZorderDst, m_puhIntraDir[ch] + uiZorderSrc, uiCurrPartNumb);
#if MLIP_EXTEND_REF_LINE
        memcpy(pCtu->getExternRef(ChannelType(ch)) + uiZorderDst, m_puhExternRef[ch] + uiZorderSrc, uiCurrPartNumb);
#endif
      }

      memcpy( pCtu->getInterDir()          + uiZorderDst, m_puhInterDir + uiZorderSrc,         uiCurrPartNumb );

#if DIMD_INTRA_PRED
      for (UInt ch = 0; ch < numValidChan; ch++)
      {
        memcpy(pCtu->getDIMDEnabledFlag(ChannelType(ch)) + uiZorderDst, m_ucDIMDEnabledFlag[ch] + uiZorderSrc, uiCurrPartNumb);
#if DIMD_NUM_INTRA_DIR_INC
        memcpy(pCtu->getExtIntraDir(ChannelType(ch)) + uiZorderDst, m_ucExtIntraDir[ch] + uiZorderSrc, uiCurrPartNumb);
#endif
        memcpy(pCtu->getDIMDNoBTLevelFlag(ChannelType(ch)) + uiZorderDst, m_ucDIMDNoBTLevelFlag[ch] + uiZorderSrc, uiCurrPartNumb);
      }
      memcpy(pCtu->getSecIntraDir() + uiZorderDst, m_puhSecIntraDir + uiZorderSrc, uiCurrPartNumb);
#endif

#if COM16_C806_EMT
      memcpy( pCtu->getEmtTuIdx()          + uiZorderDst, m_puhEmtTuIdx + uiZorderSrc,         uiCurrPartNumb );
      memcpy( pCtu->getEmtCuFlag()         + uiZorderDst, m_puhEmtCuFlag + uiZorderSrc,        uiCurrPartNumb );
#endif

      for(UInt comp=0; comp<numValidComp; comp++)
      {
        memcpy( pCtu->getCrossComponentPredictionAlpha(ComponentID(comp)) + uiZorderDst, m_crossComponentPredictionAlpha[comp] + uiZorderSrc, uiCurrPartNumb );
        memcpy( pCtu->getTransformSkip(ComponentID(comp))                 + uiZorderDst, m_puhTransformSkip[comp] + uiZorderSrc,              uiCurrPartNumb );
        memcpy( pCtu->getCbf(ComponentID(comp))                           + uiZorderDst, m_puhCbf[comp] + uiZorderSrc,                        uiCurrPartNumb );
        memcpy( pCtu->getExplicitRdpcmMode(ComponentID(comp))             + uiZorderDst, m_explicitRdpcmMode[comp] + uiZorderSrc,             uiCurrPartNumb );
#if VCEG_AZ08_KLT_COMMON
        memcpy(pCtu->getKLTFlag(ComponentID(comp))                        + uiZorderDst, m_puhKLTFlag[comp] + uiZorderSrc,             uiCurrPartNumb );
#endif
      }
#if COM16_C1016_AFFINE
      memcpy( pCtu->getAffineFlag()        + uiZorderDst, m_affineFlag + uiZorderSrc,          uiCurrPartNumb );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
      memcpy( pCtu->getAffineParamFlag() + uiZorderDst, m_affineParamFlag + uiZorderSrc, uiCurrPartNumb );
#endif

#if MRG_CONSIDER_COMPLEX
      memcpy( pCtu->getMergeComplexFlag()    + uiZorderDst, m_pbMergeComplexFlag+ uiZorderSrc,   uiCurrPartNumb  );
      memcpy( pCtu->getNodeMotionInfo()      + uiZorderDst, m_pNodeMotionInfo+ uiZorderSrc,      sizeof(* m_pNodeMotionInfo)     * uiCurrPartNumb );
      memcpy( pCtu->getNumMergeComplexCand() + uiZorderDst, m_pNumMergeComplexCand+ uiZorderSrc, sizeof(* m_pNumMergeComplexCand)     * uiCurrPartNumb );
      memcpy( pCtu->getMergeComplexIndex()   + uiZorderDst, m_puhMergeComplexIndex+ uiZorderSrc, uiCurrPartNumb );
#endif

#if RRC_SELECTABLE_IF
      memcpy( pCtu->getSIFFlag()             + uiZorderDst, m_pSIFFlag  + uiZorderSrc,uiCurrPartNumb);
#endif

#if RRC_INTRA_SHARPENING
      memcpy(pCtu->getIntraSharpFlag()            + uiZorderDst, m_sharpFlag + uiZorderSrc,uiCurrPartNumb);
#endif

      memcpy( pCtu->getDepth()  + uiZorderDst, m_puhDepth[eCType] + uiZorderSrc,  uiCurrPartNumb );
      memcpy( pCtu->getWidth()  + uiZorderDst, m_puhWidth[eCType] + uiZorderSrc,  uiCurrPartNumb );
      memcpy( pCtu->getHeight() + uiZorderDst, m_puhHeight[eCType] + uiZorderSrc, uiCurrPartNumb );
#if TT_TREE_AFFINE
      memcpy( pCtu->getCUAbsZIdxInCtu(0) + uiZorderDst, m_puhCUabsZIdxInCtu[eCType][0] + uiZorderSrc, uiCurrPartNumb);
      memcpy( pCtu->getCUAbsZIdxInCtu(1) + uiZorderDst, m_puhCUabsZIdxInCtu[eCType][1] + uiZorderSrc, uiCurrPartNumb);
#endif
      memcpy( pCtu->getBTSplitModePart(0) + uiZorderDst, m_puhBTSplitMode[eCType][0] + uiZorderSrc, uiCurrPartNumb);
      memcpy( pCtu->getBTSplitModePart(1) + uiZorderDst, m_puhBTSplitMode[eCType][1] + uiZorderSrc, uiCurrPartNumb);
#if BT_FORCE_BF
      memcpy(pCtu->getBTSplitModePart(2) + uiZorderDst, m_puhBTSplitMode[eCType][2] + uiZorderSrc, uiCurrPartNumb);
      memcpy(pCtu->getBTSplitModePart(3) + uiZorderDst, m_puhBTSplitMode[eCType][3] + uiZorderSrc, uiCurrPartNumb);
#endif

      for(UInt refListIdx=0; refListIdx<NUM_REF_PIC_LIST_01; refListIdx++)
      {
        const RefPicList rpl=RefPicList(refListIdx);
        memcpy( pCtu->getMVPIdx(rpl) + uiZorderDst, m_apiMVPIdx[rpl] + uiZorderSrc, uiCurrPartNumb );
        memcpy( pCtu->getMVPNum(rpl) + uiZorderDst, m_apiMVPNum[rpl] + uiZorderSrc, uiCurrPartNumb );
      }

      for(UInt refListIdx=0; refListIdx<NUM_REF_PIC_LIST_01; refListIdx++)
      {
        const RefPicList rpl=RefPicList(refListIdx);
        pCtu->getCUMvField(rpl)->copyFromTo(&m_acCUMvField[rpl], uiCurrPartNumb, uiZorderSrc, uiZorderDst);
      }

      memcpy( pCtu->getIPCMFlag() + uiZorderDst, m_pbIPCMFlag + uiZorderSrc,         uiCurrPartNumb  );

      uiRaster += uiNumPartInShort * uiStride;
#if TT_TREE_DS1
      uiRaster2 += uiNumPartInShort * uiStride;
#endif
    }

#if TT_TREE_DS1 //data copy recover w & h
    const UInt numCoeffY    = uiWidthSave * uiHeightSave;
#else
    const UInt numCoeffY    = uiWidth * uiHeight;
#endif
    const UInt offsetY      = m_pcPic->getCodedAreaInCTU();
    for (UInt comp=0; comp<numValidComp; comp++)
    {
      const ComponentID component = ComponentID(comp);
      const UInt componentShift   = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
      memcpy( pCtu->getCoeff(component)   + (offsetY>>componentShift), m_pcTrCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
      memcpy( pCtu->getArlCoeff(component) + (offsetY>>componentShift), m_pcArlCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
      memcpy( pCtu->getPCMSample(component) + (offsetY>>componentShift), m_pcIPCMSample[component], sizeof(Pel)*(numCoeffY>>componentShift) );
    }
  }
  else
  {
    for (UInt i=0; i<uiLong; i+=uiShort)
    {
#if TT_TREE_DS1
      uiZorderSrc = g_auiRasterToZscan[uiRaster2];
#else
      uiZorderSrc = g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu;
#endif
      uiZorderDst = g_auiRasterToZscan[uiRaster];
#if BARC
      if (pCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        memcpy( pCtu->getPredictionMode() + uiZorderDst, m_pePredMode + uiZorderSrc, uiCurrPartNumb );
#endif
      for (UInt ch=1; ch<numValidChan; ch++)
      {
        memcpy( pCtu->getIntraDir(ChannelType(ch)) + uiZorderDst, m_puhIntraDir[ch] + uiZorderSrc, uiCurrPartNumb);
      }

#if DIMD_INTRA_PRED
      for (UInt ch = 1; ch < numValidChan; ch++)
      {
        memcpy(pCtu->getDIMDEnabledFlag(ChannelType(ch)) + uiZorderDst, m_ucDIMDEnabledFlag[ch] + uiZorderSrc, uiCurrPartNumb);
#if DIMD_NUM_INTRA_DIR_INC
        memcpy(pCtu->getExtIntraDir(ChannelType(ch)) + uiZorderDst, m_ucExtIntraDir[ch] + uiZorderSrc, uiCurrPartNumb);
#endif
        memcpy(pCtu->getDIMDNoBTLevelFlag(ChannelType(ch)) + uiZorderDst, m_ucDIMDNoBTLevelFlag[ch] + uiZorderSrc, uiCurrPartNumb);
      }
#endif

      for(UInt comp=1; comp<numValidComp; comp++)
      {
        memcpy( pCtu->getCrossComponentPredictionAlpha(ComponentID(comp)) + uiZorderDst, m_crossComponentPredictionAlpha[comp] + uiZorderSrc, uiCurrPartNumb );
        memcpy( pCtu->getTransformSkip(ComponentID(comp))                 + uiZorderDst, m_puhTransformSkip[comp] + uiZorderSrc,              uiCurrPartNumb );
        memcpy( pCtu->getCbf(ComponentID(comp))                           + uiZorderDst, m_puhCbf[comp] + uiZorderSrc,                        uiCurrPartNumb );
        memcpy( pCtu->getExplicitRdpcmMode(ComponentID(comp))             + uiZorderDst, m_explicitRdpcmMode[comp] + uiZorderSrc,             uiCurrPartNumb );
#if VCEG_AZ08_KLT_COMMON
        memcpy(pCtu->getKLTFlag(ComponentID(comp))                        + uiZorderDst, m_puhKLTFlag[comp] + uiZorderSrc,             uiCurrPartNumb );
#endif
      }
#if COM16_C1044_NSST || VCEG_AZ05_ROT_TR
      for (UInt ch=1; ch<numValidChan; ch++)
      {
        memcpy( pCtu->getROTIdx(ChannelType(ch)) + uiZorderDst, m_ROTIdx[ch] + uiZorderSrc, uiCurrPartNumb);
      }
#endif
#if JVET_C0024_DELTA_QP_FIX
      memcpy( pCtu->getQP(eCType) + uiZorderDst, m_phQP[eCType] + uiZorderSrc, uiCurrPartNumb  );
#endif
      memcpy( pCtu->getDepth()  + uiZorderDst, m_puhDepth[eCType] + uiZorderSrc,  uiCurrPartNumb );
      memcpy( pCtu->getWidth()  + uiZorderDst, m_puhWidth[eCType] + uiZorderSrc,  uiCurrPartNumb );
      memcpy( pCtu->getHeight() + uiZorderDst, m_puhHeight[eCType] + uiZorderSrc, uiCurrPartNumb );
#if TT_TREE_AFFINE
      memcpy( pCtu->getCUAbsZIdxInCtu(0) + uiZorderDst, m_puhCUabsZIdxInCtu[eCType][0] + uiZorderSrc, uiCurrPartNumb);
      memcpy( pCtu->getCUAbsZIdxInCtu(1) + uiZorderDst, m_puhCUabsZIdxInCtu[eCType][1] + uiZorderSrc, uiCurrPartNumb);
#endif
      memcpy( pCtu->getBTSplitModePart(0) + uiZorderDst, m_puhBTSplitMode[eCType][0] + uiZorderSrc, uiCurrPartNumb);
      memcpy( pCtu->getBTSplitModePart(1) + uiZorderDst, m_puhBTSplitMode[eCType][1] + uiZorderSrc, uiCurrPartNumb);
#if BT_FORCE_BF
      memcpy(pCtu->getBTSplitModePart(2) + uiZorderDst, m_puhBTSplitMode[eCType][2] + uiZorderSrc, uiCurrPartNumb);
      memcpy(pCtu->getBTSplitModePart(3) + uiZorderDst, m_puhBTSplitMode[eCType][3] + uiZorderSrc, uiCurrPartNumb);
#endif

      uiRaster += uiNumPartInShort * uiStride;
#if TT_TREE_DS1
      uiRaster2 += uiNumPartInShort * uiStride;
#endif
    }

#if TT_TREE_DS1 //data copy recover w & h
    const UInt numCoeffY    = uiWidthSave * uiHeightSave;
#else
    const UInt numCoeffY    = uiWidth * uiHeight;
#endif
    const UInt offsetY      = m_pcPic->getCodedAreaInCTU();
    for (UInt comp=1; comp<numValidComp; comp++)
    {
      const ComponentID component = ComponentID(comp);
      const UInt componentShift   = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
      memcpy( pCtu->getCoeff(component)   + (offsetY>>componentShift), m_pcTrCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
      memcpy( pCtu->getArlCoeff(component) + (offsetY>>componentShift), m_pcArlCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
      memcpy( pCtu->getPCMSample(component) + (offsetY>>componentShift), m_pcIPCMSample[component], sizeof(Pel)*(numCoeffY>>componentShift) );
    }
    //maybe need to add chromaQP, getCUTransquantBypass, getIPCMFlag,..., which not impact CTC, JCA
  }
#if TT_TREE_DS1 //data copy end
  }}
#endif
#else

  Int iSizeInUchar  = sizeof( UChar ) * m_uiNumPartition;
  Int iSizeInBool   = sizeof( Bool  ) * m_uiNumPartition;
  Int sizeInChar  = sizeof( Char ) * m_uiNumPartition;
#if COM16_C806_LARGE_CTU
  Int sizeInUshort = sizeof( UShort ) * m_uiNumPartition;
#endif
  memcpy( pCtu->getSkipFlag() + m_absZIdxInCtu, m_skipFlag, sizeof( *m_skipFlag ) * m_uiNumPartition );
#if VCEG_AZ05_INTRA_MPI
  memcpy( pCtu->getMPIIdx()   + m_absZIdxInCtu, m_MPIIdx, sizeof( *m_MPIIdx ) * m_uiNumPartition );
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
  memcpy( pCtu->getPDPCIdx() + m_absZIdxInCtu, m_PDPCIdx, sizeof(*m_PDPCIdx) * m_uiNumPartition);
#endif
#if RRC_DWDIP
  memcpy(pCtu->getDWDIPIdx() + m_absZIdxInCtu, m_DWDIPIdx, sizeof(*m_DWDIPIdx) * m_uiNumPartition);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
  memcpy( pCtu->getROTIdx()   + m_absZIdxInCtu, m_ROTIdx, sizeof( *m_ROTIdx ) * m_uiNumPartition );
#endif
  memcpy( pCtu->getQP() + m_absZIdxInCtu, m_phQP, sizeInChar  );

  memcpy( pCtu->getPartitionSize()  + m_absZIdxInCtu, m_pePartSize, sizeof( *m_pePartSize ) * m_uiNumPartition );
  memcpy( pCtu->getPredictionMode() + m_absZIdxInCtu, m_pePredMode, sizeof( *m_pePredMode ) * m_uiNumPartition );
  memcpy( pCtu->getChromaQpAdj() + m_absZIdxInCtu, m_ChromaQpAdj, sizeof( *m_ChromaQpAdj ) * m_uiNumPartition );
  memcpy( pCtu->getCUTransquantBypass()+ m_absZIdxInCtu, m_CUTransquantBypass, sizeof( *m_CUTransquantBypass ) * m_uiNumPartition );
  memcpy( pCtu->getMergeFlag()         + m_absZIdxInCtu, m_pbMergeFlag,         iSizeInBool  );
  memcpy( pCtu->getMergeIndex()        + m_absZIdxInCtu, m_puhMergeIndex,       iSizeInUchar );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  memcpy( pCtu->getMergeType()         + m_absZIdxInCtu, m_peMergeType,         iSizeInUchar );
#endif
#if COM16_C806_OBMC
  memcpy( pCtu->getOBMCFlag()          + m_absZIdxInCtu, m_OBMCFlag, sizeof( *m_OBMCFlag ) * m_uiNumPartition );
#endif
#if VCEG_AZ07_FRUC_MERGE
  memcpy( pCtu->getFRUCMgrMode()       + m_absZIdxInCtu, m_puhFRUCMgrMode,      iSizeInUchar  );
#endif
#if VCEG_AZ07_IMV
  memcpy( pCtu->getiMVFlag()           + m_absZIdxInCtu, m_iMVFlag, sizeof( *m_iMVFlag ) * m_uiNumPartition );
  memcpy( pCtu->getiMVCandNum()        + m_absZIdxInCtu, m_piMVCandNum, sizeof( *m_piMVCandNum ) * m_uiNumPartition );
#endif
#if VCEG_AZ06_IC
  memcpy( pCtu->getICFlag()            + m_absZIdxInCtu, m_pbICFlag, iSizeInBool );
#endif
#if ALF_HM3_REFACTOR
  memcpy( pCtu->getAlfCtrlFlag()       + m_absZIdxInCtu, m_puiAlfCtrlFlag,      sizeof( *m_puiAlfCtrlFlag ) * m_uiNumPartition );
#endif
  for (UInt ch=0; ch<numValidChan; ch++)
  {
    memcpy( pCtu->getIntraDir(ChannelType(ch)) + m_absZIdxInCtu, m_puhIntraDir[ch], iSizeInUchar);
#if MLIP_EXTEND_REF_LINE
    memcpy(pCtu->getExternRef(ChannelType(ch)) + m_absZIdxInCtu, m_puhExternRef[ch], iSizeInUchar);
#endif
  }

  memcpy( pCtu->getInterDir()          + m_absZIdxInCtu, m_puhInterDir,         iSizeInUchar );
  memcpy( pCtu->getTransformIdx()      + m_absZIdxInCtu, m_puhTrIdx,            iSizeInUchar );

#if COM16_C806_EMT
  memcpy( pCtu->getEmtTuIdx()          + m_absZIdxInCtu, m_puhEmtTuIdx,         iSizeInUchar );
  memcpy( pCtu->getEmtCuFlag()         + m_absZIdxInCtu, m_puhEmtCuFlag,        iSizeInUchar );
#endif
  
#if COM16_C1016_AFFINE
  memcpy( pCtu->getAffineFlag()        + m_absZIdxInCtu, m_affineFlag,          iSizeInBool );
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  memcpy( pCtu->getAffineParamFlag() + m_absZIdxInCtu, m_affineParamFlag, iSizeInBool );
#endif

#if MRG_CONSIDER_COMPLEX
  memcpy( pCtu->getMergeComplexFlag()    + m_absZIdxInCtu, m_pbMergeComplexFlag,   iSizeInBool  );
  memcpy( pCtu->getNodeMotionInfo()      + m_absZIdxInCtu, m_pNodeMotionInfo,      sizeof(* m_pNodeMotionInfo) * m_uiNumPartition );
  memcpy( pCtu->getNumMergeComplexCand() + m_absZIdxInCtu, m_pNumMergeComplexCand, iSizeInUchar );
  memcpy( pCtu->getMergeComplexIndex()   + m_absZIdxInCtu, m_puhMergeComplexIndex, iSizeInUchar );
#endif

#if RRC_INTRA_SHARPENING
  memcpy( pCtu->getIntraSharpFlag()           + m_absZIdxInCtu, m_sharpFlag,           sizeInChar );
#endif

#if RRC_SELECTABLE_IF
  memcpy( pCtu->getSIFFlag()             + m_absZIdxInCtu, m_pSIFFlag,        sizeof( *m_pSIFFlag ) * m_uiNumPartition );
#endif

  for(UInt comp=0; comp<numValidComp; comp++)
  {
    memcpy( pCtu->getCrossComponentPredictionAlpha(ComponentID(comp)) + m_absZIdxInCtu, m_crossComponentPredictionAlpha[comp], iSizeInUchar );
    memcpy( pCtu->getTransformSkip(ComponentID(comp))                 + m_absZIdxInCtu, m_puhTransformSkip[comp],              iSizeInUchar );
    memcpy( pCtu->getCbf(ComponentID(comp))                           + m_absZIdxInCtu, m_puhCbf[comp],                        iSizeInUchar );
    memcpy( pCtu->getExplicitRdpcmMode(ComponentID(comp))             + m_absZIdxInCtu, m_explicitRdpcmMode[comp],             iSizeInUchar );
#if VCEG_AZ08_KLT_COMMON
    memcpy(pCtu->getKLTFlag(ComponentID(comp))                        + m_absZIdxInCtu, m_puhKLTFlag[0],                       iSizeInUchar );
#endif
  }

  memcpy( pCtu->getDepth()  + m_absZIdxInCtu, m_puhDepth,  iSizeInUchar );
#if COM16_C806_LARGE_CTU
  memcpy( pCtu->getWidth()  + m_absZIdxInCtu, m_puhWidth,  sizeInUshort );
  memcpy( pCtu->getHeight() + m_absZIdxInCtu, m_puhHeight, sizeInUshort );
#else
  memcpy( pCtu->getWidth()  + m_absZIdxInCtu, m_puhWidth,  iSizeInUchar );
  memcpy( pCtu->getHeight() + m_absZIdxInCtu, m_puhHeight, iSizeInUchar );
#endif

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    memcpy( pCtu->getMVPIdx(rpl) + m_absZIdxInCtu, m_apiMVPIdx[rpl], iSizeInUchar );
    memcpy( pCtu->getMVPNum(rpl) + m_absZIdxInCtu, m_apiMVPNum[rpl], iSizeInUchar );
  }

  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    const RefPicList rpl=RefPicList(i);
    m_acCUMvField[rpl].copyTo( pCtu->getCUMvField( rpl ), m_absZIdxInCtu );
  }

  memcpy( pCtu->getIPCMFlag() + m_absZIdxInCtu, m_pbIPCMFlag,         iSizeInBool  );

  const UInt numCoeffY    = (pCtu->getSlice()->getSPS()->getMaxCUWidth()*pCtu->getSlice()->getSPS()->getMaxCUHeight())>>(uhDepth<<1);
  const UInt offsetY      = m_absZIdxInCtu*m_pcPic->getMinCUWidth()*m_pcPic->getMinCUHeight();
  for (UInt comp=0; comp<numValidComp; comp++)
  {
    const ComponentID component = ComponentID(comp);
    const UInt componentShift   = m_pcPic->getComponentScaleX(component) + m_pcPic->getComponentScaleY(component);
    memcpy( pCtu->getCoeff(component)   + (offsetY>>componentShift), m_pcTrCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#if ADAPTIVE_QP_SELECTION
    memcpy( pCtu->getArlCoeff(component) + (offsetY>>componentShift), m_pcArlCoeff[component], sizeof(TCoeff)*(numCoeffY>>componentShift) );
#endif
    memcpy( pCtu->getPCMSample(component) + (offsetY>>componentShift), m_pcIPCMSample[component], sizeof(Pel)*(numCoeffY>>componentShift) );
  }
#endif //JVET_C0024_QTBT
  pCtu->getTotalBins() = m_uiTotalBins;
}

// --------------------------------------------------------------------------------------------------------------------
// Other public functions
// --------------------------------------------------------------------------------------------------------------------

TComDataCU* TComDataCU::getPULeft( UInt& uiLPartUnitIdx,
                                   UInt uiCurrPartUnitIdx,
                                   Bool bEnforceSliceRestriction,
                                   Bool bEnforceTileRestriction )
{
  UInt uiAbsPartIdx       = g_auiZscanToRaster[uiCurrPartUnitIdx];
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[m_absZIdxInCtu];
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();

  if ( !RasterAddress::isZeroCol( uiAbsPartIdx, numPartInCtuWidth ) )
  {
    uiLPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - 1 ];
    if ( RasterAddress::isEqualCol( uiAbsPartIdx, uiAbsZorderCUIdx, numPartInCtuWidth ) )
    {
      return m_pcPic->getCtu( getCtuRsAddr() );
    }
    else
    {
#if TT_TREE_FRUC // deriveZIdxInCU
      uiLPartUnitIdx  = deriveZIdxInCU(m_absZIdxInCtu, uiLPartUnitIdx);
#else
      uiLPartUnitIdx -= m_absZIdxInCtu;
#endif
      return this;
    }
  }

  uiLPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx + numPartInCtuWidth - 1 ];
  if ( (bEnforceSliceRestriction && !CUIsFromSameSlice(m_pCtuLeft)) || (bEnforceTileRestriction && !CUIsFromSameTile(m_pCtuLeft)) )
  {
    return NULL;
  }
#if BARC_Neighbour 
  if (m_pCtuLeft && !getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()))
  {
    uiLPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] + numPartInCtuWidth / 2 - 1];
  }
  else if (m_pCtuLeft && getBARCFlag(getTextType()) && !m_pCtuLeft->getBARCFlag(getTextType()))
  {
    uiLPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + numPartInCtuWidth - 1];
  }
  else if (m_pCtuLeft && getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()))
  {
    uiLPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdx + numPartInCtuWidth / 2 - 1];
  }
#endif
  return m_pCtuLeft;
}


TComDataCU* TComDataCU::getPUAbove( UInt& uiAPartUnitIdx,
                                    UInt uiCurrPartUnitIdx,
                                    Bool bEnforceSliceRestriction,
                                    Bool planarAtCtuBoundary,
                                    Bool bEnforceTileRestriction )
{
  UInt uiAbsPartIdx       = g_auiZscanToRaster[uiCurrPartUnitIdx];
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[m_absZIdxInCtu];
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();

  if ( !RasterAddress::isZeroRow( uiAbsPartIdx, numPartInCtuWidth ) )
  {
    uiAPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - numPartInCtuWidth ];
    if ( RasterAddress::isEqualRow( uiAbsPartIdx, uiAbsZorderCUIdx, numPartInCtuWidth ) )
    {
      return m_pcPic->getCtu( getCtuRsAddr() );
    }
    else
    {
#if TT_TREE_FRUC // deriveZIdxInCU
      uiAPartUnitIdx  = deriveZIdxInCU(m_absZIdxInCtu, uiAPartUnitIdx);
#else
      uiAPartUnitIdx -= m_absZIdxInCtu;
#endif
      return this;
    }
  }

  if(planarAtCtuBoundary)
  {
    return NULL;
  }

  uiAPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx + m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth ];

  if ( (bEnforceSliceRestriction && !CUIsFromSameSlice(m_pCtuAbove)) || (bEnforceTileRestriction && !CUIsFromSameTile(m_pCtuAbove)) )
  {
    return NULL;
  }
#if BARC_Neighbour 
  if (m_pCtuAbove && !getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
  {
    uiAPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth];
  }
  else if (m_pCtuAbove && getBARCFlag(getTextType()) && !m_pCtuAbove->getBARCFlag(getTextType()))
  {
    uiAPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth];
  }
  else if (m_pCtuAbove && getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
  {
    uiAPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdx + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth];
  }
#endif
  return m_pCtuAbove;
}

TComDataCU* TComDataCU::getPUAboveLeft( UInt& uiALPartUnitIdx, UInt uiCurrPartUnitIdx, Bool bEnforceSliceRestriction )
{
  UInt uiAbsPartIdx       = g_auiZscanToRaster[uiCurrPartUnitIdx];
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[m_absZIdxInCtu];
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();

  if ( !RasterAddress::isZeroCol( uiAbsPartIdx, numPartInCtuWidth ) )
  {
    if ( !RasterAddress::isZeroRow( uiAbsPartIdx, numPartInCtuWidth ) )
    {
      uiALPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - numPartInCtuWidth - 1 ];
      if ( RasterAddress::isEqualRowOrCol( uiAbsPartIdx, uiAbsZorderCUIdx, numPartInCtuWidth ) )
      {
        return m_pcPic->getCtu( getCtuRsAddr() );
      }
      else
      {
#if TT_TREE_FRUC // deriveZIdxInCU
        uiALPartUnitIdx  = deriveZIdxInCU(m_absZIdxInCtu, uiALPartUnitIdx);
#else
        uiALPartUnitIdx -= m_absZIdxInCtu;
#endif
        return this;
      }
    }
    uiALPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx + getPic()->getNumPartitionsInCtu() - numPartInCtuWidth - 1 ];
    if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuAbove) )
    {
      return NULL;
    }
#if BARC_Neighbour 
    if (m_pCtuAbove && !getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth - 1];
    }
    else if (m_pCtuAbove && getBARCFlag(getTextType()) && !m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth - 1];
    }
    else if (m_pCtuAbove && getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdx + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth - 1];
    }
#endif
    return m_pCtuAbove;
  }

  if ( !RasterAddress::isZeroRow( uiAbsPartIdx, numPartInCtuWidth ) )
  {
    uiALPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - 1 ];
    if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuLeft) )
    {
      return NULL;
    }
#if BARC_Neighbour
    if (m_pCtuLeft && !getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()) && !RasterAddress::isZeroRow(g_auiZscanToRaster[(uiCurrPartUnitIdx >> 2)], numPartInCtuWidth))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] - numPartInCtuWidth / 2 - 1];
    }
    else if (m_pCtuLeft && getBARCFlag(getTextType()) && !m_pCtuLeft->getBARCFlag(getTextType()))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + numPartInCtuWidth - 1];
      return m_pCtuLeft;
    }
    else if (m_pCtuLeft && getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()))
    {
      uiALPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdx - numPartInCtuWidth / 2 - 1];
      return m_pCtuLeft;
    }
#endif
    return m_pCtuLeft;
  }

  uiALPartUnitIdx = g_auiRasterToZscan[ m_pcPic->getNumPartitionsInCtu() - 1 ];
  if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuAboveLeft) )
  {
    return NULL;
  }
#if BARC_Neighbour 
  if (m_pCtuAboveLeft && m_pCtuAboveLeft->getBARCFlag(getTextType()))
  {
    uiALPartUnitIdx = g_auiRasterToZscan[  m_pcPic->getNumPartitionsInCtu()/2 - 1 - numPartInCtuWidth/2 ];
  }
#endif
  return m_pCtuAboveLeft;
}

TComDataCU* TComDataCU::getPUBelowLeft(UInt& uiBLPartUnitIdx,  UInt uiCurrPartUnitIdx, UInt uiPartUnitOffset, Bool bEnforceSliceRestriction)
{
  UInt uiAbsPartIdxLB     = g_auiZscanToRaster[uiCurrPartUnitIdx];
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  UInt uiAbsZorderCUIdxLB = g_auiZscanToRaster[ m_absZIdxInCtu ] + ((getHeight(0) / m_pcPic->getMinCUHeight()) - 1)*numPartInCtuWidth;
#else
  UInt uiAbsZorderCUIdxLB = g_auiZscanToRaster[ m_absZIdxInCtu ] + ((m_puhHeight[getTextType()][0] / m_pcPic->getMinCUHeight()) - 1)*numPartInCtuWidth;
#endif
#else
  UInt uiAbsZorderCUIdxLB = g_auiZscanToRaster[ m_absZIdxInCtu ] + ((m_puhHeight[0] / m_pcPic->getMinCUHeight()) - 1)*numPartInCtuWidth;
#endif
#if BARC_Boundary
  UInt isBarc = m_pcPic->getCtu(m_ctuRsAddr)->getBARCFlag(getTextType());
  UInt lScale = isBarc;
  if( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + (g_auiRasterToPelY[uiAbsPartIdxLB])*(1+lScale) + (m_pcPic->getPicSym()->getMinCUHeight() * uiPartUnitOffset*(1+lScale))) >= m_pcSlice->getSPS()->getPicHeightInLumaSamples())
#else
  if( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + g_auiRasterToPelY[uiAbsPartIdxLB] + (m_pcPic->getPicSym()->getMinCUHeight() * uiPartUnitOffset)) >= m_pcSlice->getSPS()->getPicHeightInLumaSamples())
#endif
  {
    uiBLPartUnitIdx = MAX_UINT;
    return NULL;
  }
#if BARC_Neighbour
  if ( RasterAddress::lessThanRow( uiAbsPartIdxLB, m_pcPic->getNumPartInCtuHeight()/(1+isBarc) - uiPartUnitOffset, numPartInCtuWidth ) )
#else
  if ( RasterAddress::lessThanRow( uiAbsPartIdxLB, m_pcPic->getNumPartInCtuHeight() - uiPartUnitOffset, numPartInCtuWidth ) )
#endif
  {
    if ( !RasterAddress::isZeroCol( uiAbsPartIdxLB, numPartInCtuWidth ) )
    {
#if JVET_C0024_QTBT
      UInt uiBlkX = g_auiRasterToPelX[ uiAbsPartIdxLB + uiPartUnitOffset*numPartInCtuWidth - 1 ] >>MIN_CU_LOG2;
      UInt uiBlkY = g_auiRasterToPelY[ uiAbsPartIdxLB + uiPartUnitOffset*numPartInCtuWidth - 1 ] >>MIN_CU_LOG2;
      if (getPic()->getCodedBlkInCTU(uiBlkX, uiBlkY))
#else
      if ( uiCurrPartUnitIdx > g_auiRasterToZscan[ uiAbsPartIdxLB + uiPartUnitOffset * numPartInCtuWidth - 1 ] )
#endif
      {
        uiBLPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdxLB + uiPartUnitOffset * numPartInCtuWidth - 1 ];
        if ( RasterAddress::isEqualRowOrCol( uiAbsPartIdxLB, uiAbsZorderCUIdxLB, numPartInCtuWidth ) )
        {
          return m_pcPic->getCtu( getCtuRsAddr() );
        }
        else
        {
#if TT_TREE_FRUC // deriveZIdxInCU
          uiBLPartUnitIdx  = deriveZIdxInCU(m_absZIdxInCtu, uiBLPartUnitIdx);
#else
          uiBLPartUnitIdx -= m_absZIdxInCtu;
#endif
          return this;
        }
      }
      uiBLPartUnitIdx = MAX_UINT;
      return NULL;
    }
    uiBLPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdxLB + (1+uiPartUnitOffset) * numPartInCtuWidth - 1 ];
    if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuLeft) )
    {
      return NULL;
    }
#if BARC_Neighbour
    if (m_pCtuLeft && !getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()))
    {
      uiBLPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] + ((1 + uiPartUnitOffset) >> 1) * numPartInCtuWidth - numPartInCtuWidth / 2 - 1];
    }
    else if (m_pCtuLeft && getBARCFlag(getTextType()) && !m_pCtuLeft->getBARCFlag(getTextType()))
    {
      uiBLPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + ((1 + uiPartUnitOffset) << 1) * numPartInCtuWidth - 1];
    }
    else if (m_pCtuLeft && getBARCFlag(getTextType()) && m_pCtuLeft->getBARCFlag(getTextType()))
    {
      uiBLPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdxLB + (1 + uiPartUnitOffset) * numPartInCtuWidth - numPartInCtuWidth / 2 - 1];
    }
#endif
    return m_pCtuLeft;
  }

  uiBLPartUnitIdx = MAX_UINT;
  return NULL;
}

TComDataCU* TComDataCU::getPUAboveRight(UInt&  uiARPartUnitIdx, UInt uiCurrPartUnitIdx, UInt uiPartUnitOffset, Bool bEnforceSliceRestriction)
{
  UInt uiAbsPartIdxRT     = g_auiZscanToRaster[uiCurrPartUnitIdx];
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[ m_absZIdxInCtu ] + (getWidth(0) / m_pcPic->getMinCUWidth()) - 1;
#else
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[ m_absZIdxInCtu ] + (m_puhWidth[getTextType()][0] / m_pcPic->getMinCUWidth()) - 1;
#endif
#else
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[ m_absZIdxInCtu ] + (m_puhWidth[0] / m_pcPic->getMinCUWidth()) - 1;
#endif
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();

#if BARC_Boundary
  UInt isBarc = m_pcPic->getCtu(m_ctuRsAddr)->getBARCFlag(getTextType());
  UInt aScale = isBarc;
  if( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + (g_auiRasterToPelX[uiAbsPartIdxRT])*(1+aScale) + (m_pcPic->getPicSym()->getMinCUHeight() * uiPartUnitOffset*(1+aScale))) >= m_pcSlice->getSPS()->getPicWidthInLumaSamples() )
#else
  if( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + g_auiRasterToPelX[uiAbsPartIdxRT] + (m_pcPic->getPicSym()->getMinCUHeight() * uiPartUnitOffset)) >= m_pcSlice->getSPS()->getPicWidthInLumaSamples() )
#endif
  {
    uiARPartUnitIdx = MAX_UINT;
    return NULL;
  }
#if BARC_Neighbour
  if ( RasterAddress::lessThanCol( uiAbsPartIdxRT, numPartInCtuWidth/(1+isBarc) - uiPartUnitOffset, numPartInCtuWidth ) )
#else
  if ( RasterAddress::lessThanCol( uiAbsPartIdxRT, numPartInCtuWidth - uiPartUnitOffset, numPartInCtuWidth ) )
#endif
  {
    if ( !RasterAddress::isZeroRow( uiAbsPartIdxRT, numPartInCtuWidth ) )
    {
#if JVET_C0024_QTBT
      UInt uiBlkX = g_auiRasterToPelX[ uiAbsPartIdxRT - numPartInCtuWidth + uiPartUnitOffset ] >>MIN_CU_LOG2;
      UInt uiBlkY = g_auiRasterToPelY[ uiAbsPartIdxRT - numPartInCtuWidth + uiPartUnitOffset ] >>MIN_CU_LOG2;
      if (getPic()->getCodedBlkInCTU(uiBlkX, uiBlkY) )
#else
      if ( uiCurrPartUnitIdx > g_auiRasterToZscan[ uiAbsPartIdxRT - numPartInCtuWidth + uiPartUnitOffset ] )
#endif
      {
        uiARPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdxRT - numPartInCtuWidth + uiPartUnitOffset ];
        if ( RasterAddress::isEqualRowOrCol( uiAbsPartIdxRT, uiAbsZorderCUIdx, numPartInCtuWidth ) )
        {
          return m_pcPic->getCtu( getCtuRsAddr() );
        }
        else
        {
#if TT_TREE_FRUC // deriveZIdxInCU
          uiARPartUnitIdx  = deriveZIdxInCU(m_absZIdxInCtu, uiARPartUnitIdx);
#else
          uiARPartUnitIdx -= m_absZIdxInCtu;
#endif
          return this;
        }
      }
      uiARPartUnitIdx = MAX_UINT;
      return NULL;
    }

    uiARPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdxRT + m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth + uiPartUnitOffset ];
    if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuAbove) )
    {
      return NULL;
    }
#if BARC_Neighbour 
    if (m_pCtuAbove && !getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiARPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx >> 2] + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth + (uiPartUnitOffset >> 1)];
    }
    else if (m_pCtuAbove && getBARCFlag(getTextType()) && !m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiARPartUnitIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiCurrPartUnitIdx << 2] + m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth + (uiPartUnitOffset << 1)];
    }
    else if (m_pCtuAbove && getBARCFlag(getTextType()) && m_pCtuAbove->getBARCFlag(getTextType()))
    {
      uiARPartUnitIdx = g_auiRasterToZscan[uiAbsPartIdxRT + m_pcPic->getNumPartitionsInCtu() / 2 - numPartInCtuWidth + uiPartUnitOffset];
    }
#endif
    return m_pCtuAbove;
  }

  if ( !RasterAddress::isZeroRow( uiAbsPartIdxRT, numPartInCtuWidth ) )
  {
    uiARPartUnitIdx = MAX_UINT;
    return NULL;
  }

  uiARPartUnitIdx = g_auiRasterToZscan[ m_pcPic->getNumPartitionsInCtu() - numPartInCtuWidth + uiPartUnitOffset-1 ];
  if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuAboveRight) )
  {
    return NULL;
  }
#if BARC_Neighbour
  if (m_pCtuAboveRight && m_pCtuAboveRight->getBARCFlag(getTextType()))
  {
    uiARPartUnitIdx = g_auiRasterToZscan[ m_pcPic->getNumPartitionsInCtu()/2 - numPartInCtuWidth + uiPartUnitOffset/2-1 ];
  }
#endif
  return m_pCtuAboveRight;
}

#if DIMD_INTRA_PRED
#if DIMD_REDUCE_LINE_MEMORY_JMB
TEMPLATE_TYPE TComDataCU::deriveRefRegPosSize(Int iCurX, Int iCurY, UInt uiCurWidth, UInt uiCurHeight, Int iTemplateWidth, Int iTemplateHeight, Int& iRefX, Int& iRefY, UInt& uiRefWidth, UInt& uiRefHeight, UInt uiCTUSize)
#else
TEMPLATE_TYPE TComDataCU::deriveRefRegPosSize(Int iCurX, Int iCurY, UInt uiCurWidth, UInt uiCurHeight, Int iTemplateWidth, Int iTemplateHeight, Int& iRefX, Int& iRefY, UInt& uiRefWidth, UInt& uiRefHeight)
#endif
{
  //assert(iCurX >= 0 && iCurY >= 0);

  if (iCurX == 0 && iCurY == 0)
  {
    return NO_NEIGHBOR;
  }

  TEMPLATE_TYPE eTempType = NO_NEIGHBOR;
  iRefX = iRefY = -1;
  if (iCurX > 0 && iCurY > 0)
  {
    //assert(iCurX >= iTemplateWidth && iCurY >= iTemplateHeight);
    iRefX = iCurX - iTemplateWidth;
    iRefY = iCurY - iTemplateHeight;
    uiRefWidth = uiCurWidth + iTemplateWidth;
    uiRefHeight = uiCurHeight + iTemplateHeight;
    eTempType = LEFT_ABOVE_NEIGHBOR;
#if DIMD_REDUCE_LINE_MEMORY_JMB
    if (iCurY % uiCTUSize == 0) /* Check if the current CU location aligns with the CTU boundary */
    {
      eTempType = LEFT_NEIGHBOR;
    }
#endif
  }
  else if (iCurX == 0 && iCurY > 0)
  {
    //assert(iCurY >= iTemplateHeight);
    iRefX = iCurX;
    iRefY = iCurY - iTemplateHeight;
    uiRefWidth = uiCurWidth;
    uiRefHeight = uiCurHeight;
    eTempType = ABOVE_NEIGHBOR;
#if DIMD_REDUCE_LINE_MEMORY_JMB
    if (iCurY % uiCTUSize == 0) /* Check if the current CU location aligns with the CTU boundary */
    {
      eTempType = NO_NEIGHBOR;
    }
#endif
  }
  else if (iCurX > 0 && iCurY == 0)
  {
    //assert(iCurX >= iTemplateWidth);
    iRefX = iCurX - iTemplateWidth;
    iRefY = iCurY;
    uiRefWidth = uiCurWidth;
    uiRefHeight = uiCurHeight;
    eTempType = LEFT_NEIGHBOR;
  }
  else
  {
    assert(0);
  }
  //assert(iRefX >= 0 && iRefY >= 0);

  return eTempType;
}

TComDataCU* TComDataCU::getOrigRefReg(UInt& uiRefPartUnitIdx, Int iPelX, Int iPelY
#if BARC_DIMD
  ,Bool isBARC
#endif
)
{
  UInt uiUnitSize = getPic()->getMinCUWidth();
  UInt uiPicWidthInLCU = getPic()->getFrameWidthInCtus();
  Int  iMaxCUWidth = getSlice()->getSPS()->getCTUSize();
  Int  iMaxCUHeight = getSlice()->getSPS()->getCTUSize();
#if BARC_DIMD
  UInt uiPicWidth = getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt uiPicHeight = getSlice()->getSPS()->getPicHeightInLumaSamples();
  Int deltaX = 0;
  Int deltaY = 0;
  if (isBARC)
  {
    deltaX = uiPicWidth % iMaxCUWidth ? (uiPicWidth % iMaxCUWidth) / 2 : iMaxCUWidth / 2;
    deltaY = uiPicHeight % iMaxCUHeight ? (uiPicHeight % iMaxCUHeight) / 2 : iMaxCUHeight / 2;
  }
#endif

  if (iPelX < 0 || iPelX >= getSlice()->getSPS()->getPicWidthInLumaSamples()
#if BARC_DIMD
    - deltaX
#endif
    )
  {
    return NULL;
  }
  if (iPelY < 0 || iPelY >= getSlice()->getSPS()->getPicHeightInLumaSamples()
#if BARC_DIMD
    - deltaY
#endif
    )
  {
    return NULL;
  }

  Int iLCUHIdx = iPelX / iMaxCUWidth;
  Int iLCUVIdx = iPelY / iMaxCUHeight;
  TComDataCU* pcLCU = getPic()->getCtu(iLCUVIdx*uiPicWidthInLCU + iLCUHIdx);
  //assert(pcLCU != NULL);
  //assert(pcLCU->getZorderIdxInCtu() == 0);

#if BARC_DIMD
  if (isBARC)
  {
    UInt iCurCTUAddr = m_pcTargetCtu->getCtuRsAddr();
   // UInt iCurPixAboveRight = (iCurCTUAddr%uiPicWidthInLCU)*iMaxCUWidth + iMaxCUWidth / 2;
    UInt iCurPixBelowLeft = (iCurCTUAddr/uiPicWidthInLCU)*iMaxCUHeight + iMaxCUHeight / 2;
    if (iPelY >= iCurPixBelowLeft)
      return NULL;
  }
#endif
  Int iPelXInLCU = iPelX - iLCUHIdx*iMaxCUWidth;
  Int iPelYInLCU = iPelY - iLCUVIdx*iMaxCUHeight;
  Int iRefRasterIdx = iPelYInLCU / uiUnitSize*(getPic()->getNumPartInCtuWidth()) + (iPelXInLCU / uiUnitSize);
  uiRefPartUnitIdx = g_auiRasterToZscan[iRefRasterIdx];

  return pcLCU;
}

#if BARC_DIMD
Void TComDataCU::getBARCRefReg(UInt refCuAddr, UInt& uiBARCRefPartUnitIdx, Int iPelX, Int iPelY)
{
  UInt  iMaxCUWidth = getSlice()->getSPS()->getCTUSize();
  UInt  iMaxCUHeight = getSlice()->getSPS()->getCTUSize();
  UInt  iWidthPicLumaSamples = getSlice()->getSPS()->getPicWidthInLumaSamples();
  UInt  iHeightPicLumaSamples = getSlice()->getSPS()->getPicHeightInLumaSamples();

  UInt uiUnitSize = getPic()->getMinCUWidth();
  UInt uiPicWidthInLCU = getPic()->getFrameWidthInCtus();

  Int iLCUHIdx = iPelX / iMaxCUWidth;
  Int iLCUVIdx = iPelY / iMaxCUHeight;

  UInt deltaX = 0;
  UInt deltaY = 0;

  if ((refCuAddr % uiPicWidthInLCU) == (getCtuRsAddr()%uiPicWidthInLCU)) //Above
  {
    deltaY = iMaxCUHeight / 2;
  }
  else if (refCuAddr == (getCtuRsAddr() - 1)) //Left
  {
    deltaX = iMaxCUWidth / 2;
  }
  else //AboveLeft
  {
    deltaX = iMaxCUWidth / 2;
    deltaY = iMaxCUHeight / 2;
  }
  Int iPelXInLCU = iPelX - iLCUHIdx*iMaxCUWidth - deltaX;
  Int iPelYInLCU = iPelY - iLCUVIdx*iMaxCUHeight - deltaY;
  Int iRefRasterIdx = iPelYInLCU / uiUnitSize*(getPic()->getNumPartInCtuWidth()) + (iPelXInLCU / uiUnitSize);
  uiBARCRefPartUnitIdx = g_auiRasterToZscan[iRefRasterIdx];
}
#endif

Bool TComDataCU::isReferenceAvailable(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
  //assert(pcCU->getCtuRsAddr() != MAX_UINT);
  if (pcCU->getCtuRsAddr() < m_pcTargetCtu->getCtuRsAddr())
  {
    return true;
  }

  if (pcCU->getCtuRsAddr() > m_pcTargetCtu->getCtuRsAddr())
  {
    return false;
  }

#if TT_TREE_DIMD // deriveZIdxInCTU
  UInt uiRefZOrder = pcCU->deriveZIdxInCTU(pcCU->getZorderIdxInCtu(), uiAbsPartIdx);
#else
  UInt uiRefZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
#endif


  //UInt uiRefZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
  UInt uiBlkX = g_auiRasterToPelX[g_auiZscanToRaster[uiRefZOrder]] >> MIN_CU_LOG2;
  UInt uiBlkY = g_auiRasterToPelY[g_auiZscanToRaster[uiRefZOrder]] >> MIN_CU_LOG2;
  if (getPic()->getCodedBlkInCTU(uiBlkX, uiBlkY))
  {
    return true;
  }
  else
  {
    return false;
  }
}

Int TComDataCU::generateSecIntraDir(UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight)
{

#if TT_TREE_DIMD // deriveZIdxInCTU
  UInt partIdxLT, partIdxRT, partIdxLB;
  deriveLeftRightTopIdxGeneral(uiAbsPartIdx, 0, partIdxLT, partIdxRT);
  deriveLeftBottomIdxGeneral(uiAbsPartIdx, 0, partIdxLB);
#else
  UInt partIdxLT = m_absZIdxInCtu + uiAbsPartIdx;
  UInt partIdxRT = g_auiRasterToZscan[g_auiZscanToRaster[partIdxLT] + uiWidth / getPic()->getMinCUWidth() - 1];
  UInt partIdxLB = g_auiRasterToZscan[g_auiZscanToRaster[partIdxLT] + ((uiHeight / getPic()->getMinCUHeight()) - 1) * getPic()->getNumPartInCtuWidth()];
#endif


  UInt uiIntraDirMaxCount = 0, uiMaxIntraDir = PLANAR_IDX;
  UInt uiIntraDirCnt[NUM_INTRA_MODE], partIdx;
  memset(uiIntraDirCnt, 0, sizeof(UInt)*NUM_INTRA_MODE);

  // left
  TComDataCU *cu = getPULeft(partIdx, partIdxLB);
  if (cu && cu->isIntra(partIdx))
  {
    UChar ucTempIntraDir = cu->getSecIntraDir(partIdx);
    uiIntraDirCnt[ucTempIntraDir] ++;
    if (uiIntraDirMaxCount < uiIntraDirCnt[ucTempIntraDir])
    {
      uiIntraDirMaxCount = uiIntraDirCnt[ucTempIntraDir];
      uiMaxIntraDir = ucTempIntraDir;
    }
  }

  // above
  cu = getPUAbove(partIdx, partIdxRT);
  if (cu && cu->isIntra(partIdx))
  {
    UChar ucTempIntraDir = cu->getSecIntraDir(partIdx);
    uiIntraDirCnt[ucTempIntraDir] ++;
    if (uiIntraDirMaxCount < uiIntraDirCnt[ucTempIntraDir])
    {
      uiIntraDirMaxCount = uiIntraDirCnt[ucTempIntraDir];
      uiMaxIntraDir = ucTempIntraDir;
    }
  }

  // below left
  cu = getPUBelowLeft(partIdx, partIdxLB);
  if (cu && cu->isIntra(partIdx))
  {
    UChar ucTempIntraDir = cu->getSecIntraDir(partIdx);
    uiIntraDirCnt[ucTempIntraDir] ++;
    if (uiIntraDirMaxCount < uiIntraDirCnt[ucTempIntraDir])
    {
      uiIntraDirMaxCount = uiIntraDirCnt[ucTempIntraDir];
      uiMaxIntraDir = ucTempIntraDir;
    }
  }

  // above right
  cu = getPUAboveRight(partIdx, partIdxRT);
  if (cu && cu->isIntra(partIdx))
  {
    UChar ucTempIntraDir = cu->getSecIntraDir(partIdx);
    uiIntraDirCnt[ucTempIntraDir] ++;
    if (uiIntraDirMaxCount < uiIntraDirCnt[ucTempIntraDir])
    {
      uiIntraDirMaxCount = uiIntraDirCnt[ucTempIntraDir];
      uiMaxIntraDir = ucTempIntraDir;
    }
  }

  //above left
  cu = getPUAboveLeft(partIdx, partIdxLT);
  if (cu && cu->isIntra(partIdx))
  {
    UChar ucTempIntraDir = cu->getSecIntraDir(partIdx);
    uiIntraDirCnt[ucTempIntraDir] ++;
    if (uiIntraDirMaxCount < uiIntraDirCnt[ucTempIntraDir])
    {
      uiIntraDirMaxCount = uiIntraDirCnt[ucTempIntraDir];
      uiMaxIntraDir = ucTempIntraDir;
    }
  }
  return uiMaxIntraDir;
}
#endif


#if VCEG_AZ07_FRUC_MERGE
Bool TComDataCU::getBlockBelowRight( UInt uiAbsPartIdx, Int nCurBlkWidth , Int nCurBlkHeight , UInt & rCUAddr , UInt & rBRAbsPartIdx )
{
  assert( MIN_PU_SIZE == 4 );
  Int x = getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]] + nCurBlkWidth;
  Int y = getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]] + nCurBlkHeight;
  if( x >= getPic()->getPicYuvRec()->getWidth( COMPONENT_Y ) || y >= getPic()->getPicYuvRec()->getHeight( COMPONENT_Y ) )
    return( false );

#if JVET_C0024_QTBT
  rCUAddr = getPic()->getFrameWidthInCtus() * ( y >> ( g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] + MIN_CU_LOG2 ) ) + ( x >> ( g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] + MIN_CU_LOG2 ) );
  Int nMask = ( 1 << ( g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] + MIN_CU_LOG2 ) ) - 1;
#else
  rCUAddr = getPic()->getFrameWidthInCtus() * ( y >> ( g_aucConvertToBit[getSlice()->getSPS()->getMaxCUWidth()] + 2 ) ) + ( x >> ( g_aucConvertToBit[getSlice()->getSPS()->getMaxCUWidth()] + 2 ) );
  Int nMask = ( 1 << ( g_aucConvertToBit[getSlice()->getSPS()->getMaxCUWidth()] + 2 ) ) - 1;
#endif
  UInt uiRasterIdx = ( ( y & nMask ) >> 2 ) * getPic()->getNumPartInCtuWidth() + ( ( x & nMask ) >> 2 );
  rBRAbsPartIdx = g_auiRasterToZscan[uiRasterIdx];
  return( true );
}
#endif

/** Get left QpMinCu
*\param   uiLPartUnitIdx
*\param   uiCurrAbsIdxInCtu
*\returns TComDataCU*   point of TComDataCU of left QpMinCu
*/
TComDataCU* TComDataCU::getQpMinCuLeft( UInt& uiLPartUnitIdx, UInt uiCurrAbsIdxInCtu )
{
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();
#if JVET_C0024_DELTA_QP_FIX
  UInt absZorderQpMinCUIdx = m_uiQuPartIdx;
  UInt absRorderQpMinCUIdx = g_auiZscanToRaster[absZorderQpMinCUIdx];
#else
  const UInt maxCUDepth        = getSlice()->getSPS()->getMaxTotalCUDepth();
  const UInt maxCuDQPDepth     = getSlice()->getPPS()->getMaxCuDQPDepth();
  const UInt doubleDepthDifference = ((maxCUDepth - maxCuDQPDepth)<<1);
  UInt absZorderQpMinCUIdx = (uiCurrAbsIdxInCtu>>doubleDepthDifference)<<doubleDepthDifference;
  UInt absRorderQpMinCUIdx = g_auiZscanToRaster[absZorderQpMinCUIdx];
#endif

  // check for left CTU boundary
  if ( RasterAddress::isZeroCol(absRorderQpMinCUIdx, numPartInCtuWidth) )
  {
    return NULL;
  }

  // get index of left-CU relative to top-left corner of current quantization group
  uiLPartUnitIdx = g_auiRasterToZscan[absRorderQpMinCUIdx - 1];

  // return pointer to current CTU
  return m_pcPic->getCtu( getCtuRsAddr() );
}

/** Get Above QpMinCu
*\param   uiAPartUnitIdx
*\param   uiCurrAbsIdxInCtu
*\returns TComDataCU*   point of TComDataCU of above QpMinCu
*/
TComDataCU* TComDataCU::getQpMinCuAbove( UInt& uiAPartUnitIdx, UInt uiCurrAbsIdxInCtu )
{
  const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();
#if JVET_C0024_DELTA_QP_FIX
  UInt absZorderQpMinCUIdx = m_uiQuPartIdx;
  UInt absRorderQpMinCUIdx = g_auiZscanToRaster[absZorderQpMinCUIdx];
#else
  const UInt maxCUDepth        = getSlice()->getSPS()->getMaxTotalCUDepth();
  const UInt maxCuDQPDepth     = getSlice()->getPPS()->getMaxCuDQPDepth();
  const UInt doubleDepthDifference = ((maxCUDepth - maxCuDQPDepth)<<1);
  UInt absZorderQpMinCUIdx = (uiCurrAbsIdxInCtu>>doubleDepthDifference)<<doubleDepthDifference;
  UInt absRorderQpMinCUIdx = g_auiZscanToRaster[absZorderQpMinCUIdx];
#endif

  // check for top CTU boundary
  if ( RasterAddress::isZeroRow( absRorderQpMinCUIdx, numPartInCtuWidth) )
  {
    return NULL;
  }

  // get index of top-CU relative to top-left corner of current quantization group
  uiAPartUnitIdx = g_auiRasterToZscan[absRorderQpMinCUIdx - numPartInCtuWidth];

  // return pointer to current CTU
  return m_pcPic->getCtu( getCtuRsAddr() );
}



/** Get reference QP from left QpMinCu or latest coded QP
*\param   uiCurrAbsIdxInCtu
*\returns Char   reference QP value
*/
Char TComDataCU::getRefQP( UInt uiCurrAbsIdxInCtu )
{
  UInt lPartIdx = MAX_UINT;
  UInt aPartIdx = MAX_UINT;
#if TT_TREE_DS1 // deriveZIdxInCTU
  TComDataCU* cULeft  = getQpMinCuLeft ( lPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiCurrAbsIdxInCtu) );
  TComDataCU* cUAbove = getQpMinCuAbove( aPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiCurrAbsIdxInCtu) );
#else
  TComDataCU* cULeft  = getQpMinCuLeft ( lPartIdx, m_absZIdxInCtu + uiCurrAbsIdxInCtu );
  TComDataCU* cUAbove = getQpMinCuAbove( aPartIdx, m_absZIdxInCtu + uiCurrAbsIdxInCtu );
#endif
#if JVET_C0024_DELTA_QP_FIX
  return (((cULeft? cULeft->getQP( lPartIdx ): m_QuLastCodedQP) + (cUAbove? cUAbove->getQP( aPartIdx ): m_QuLastCodedQP) + 1) >> 1);
#else
  return (((cULeft? cULeft->getQP( lPartIdx ): getLastCodedQP( uiCurrAbsIdxInCtu )) + (cUAbove? cUAbove->getQP( aPartIdx ): getLastCodedQP( uiCurrAbsIdxInCtu )) + 1) >> 1);
#endif
}

#if JVET_C0024_DELTA_QP_FIX
Char TComDataCU::getCtuLastCodedQP()
{
  if( getPic()->getPicSym()->getCtuTsToRsAddrMap(getSlice()->getSliceCurStartCtuTsAddr()) != getCtuRsAddr() && getPic()->getPicSym()->getCtuRsToTsAddrMap(getCtuRsAddr()) > 0 )
  {
    TComDataCU* pcCtuPrev = getPic()->getCtu(getPic()->getPicSym()->getCtuTsToRsAddrMap(getPic()->getPicSym()->getCtuRsToTsAddrMap(getCtuRsAddr())-1));
    if( CUIsFromSameSliceTileAndWavefrontRow( pcCtuPrev ) )
    {
      return pcCtuPrev->getCodedQP();
    }
  }
  return getSlice()->getSliceQp();
}
#else

Int TComDataCU::getLastValidPartIdx( Int iAbsPartIdx )
{
  Int iLastValidPartIdx = iAbsPartIdx-1;
  while ( iLastValidPartIdx >= 0
       && getPredictionMode( iLastValidPartIdx ) == NUMBER_OF_PREDICTION_MODES )
  {
    UInt uiDepth = getDepth( iLastValidPartIdx );
    iLastValidPartIdx -= m_uiNumPartition>>(uiDepth<<1);
  }
  return iLastValidPartIdx;
}

Char TComDataCU::getLastCodedQP( UInt uiAbsPartIdx )
{
  UInt uiQUPartIdxMask = ~((1<<((getSlice()->getSPS()->getMaxTotalCUDepth() - getSlice()->getPPS()->getMaxCuDQPDepth())<<1))-1);
  Int iLastValidPartIdx = getLastValidPartIdx( uiAbsPartIdx&uiQUPartIdxMask ); // A idx will be invalid if it is off the right or bottom edge of the picture.
  // If this CU is in the first CTU of the slice and there is no valid part before this one, use slice QP
  if ( getPic()->getPicSym()->getCtuTsToRsAddrMap(getSlice()->getSliceCurStartCtuTsAddr()) == getCtuRsAddr() && Int(getZorderIdxInCtu())+iLastValidPartIdx<0)
  {
    return getSlice()->getSliceQp();
  }
  else if ( iLastValidPartIdx >= 0 )
  {
    // If there is a valid part within the current Sub-CU, use it
    return getQP( iLastValidPartIdx );
  }
  else
  {
    if ( getZorderIdxInCtu() > 0 )
    {
      // If this wasn't the first sub-cu within the Ctu, explore the CTU itself.
      return getPic()->getCtu( getCtuRsAddr() )->getLastCodedQP( getZorderIdxInCtu() ); // TODO - remove this recursion
    }
    else if ( getPic()->getPicSym()->getCtuRsToTsAddrMap(getCtuRsAddr()) > 0
      && CUIsFromSameSliceTileAndWavefrontRow(getPic()->getCtu(getPic()->getPicSym()->getCtuTsToRsAddrMap(getPic()->getPicSym()->getCtuRsToTsAddrMap(getCtuRsAddr())-1))) )
    {
      // If this isn't the first Ctu (how can it be due to the first 'if'?), and the previous Ctu is from the same tile, examine the previous Ctu.
      return getPic()->getCtu( getPic()->getPicSym()->getCtuTsToRsAddrMap(getPic()->getPicSym()->getCtuRsToTsAddrMap(getCtuRsAddr())-1) )->getLastCodedQP( getPic()->getNumPartitionsInCtu() );  // TODO - remove this recursion
    }
    else
    {
      // No other options available - use the slice-level QP.
      return getSlice()->getSliceQp();
    }
  }
}
#endif

/** Check whether the CU is coded in lossless coding mode.
 * \param   absPartIdx
 * \returns true if the CU is coded in lossless coding mode; false if otherwise
 */
Bool TComDataCU::isLosslessCoded(UInt absPartIdx)
{
  return (getSlice()->getPPS()->getTransquantBypassEnableFlag() && getCUTransquantBypass (absPartIdx));
}

#if JVET_E0062_MULTI_DMS
UInt TComDataCU::getDMMode  ( UInt uiAbsPartIdx, UInt uiDMIdx, UInt uiChDMMode[NUM_DM_MODES], UInt* iTotalCnt ) 
{
  UInt uiChFinalMode;
  Int  uiChMode[NUM_MOST_PROBABLE_MODES], i = 0;
  Bool bIntraSlice = getSlice()->isIntra();
  for (i = 0; i< NUM_MOST_PROBABLE_MODES; i++)
  {
    uiChMode[i] = INVALID_CHROMAMODE;  
  }
#if JVET_C0024_QTBT
  if(bIntraSlice)
  {
    if(uiChDMMode)
    {
      for(i = 0; i< NUM_DM_MODES; i++)
      {
        uiChDMMode[i] = INVALID_CHROMAMODE;
      }
    }
    const TComDataCU* pcCtu = getPic()->getCtu(getCtuRsAddr());
    UInt uiAbsZeroIdx = getZorderIdxInCtu();
#if TT_TREE_DS5 // deriveZIdxInCTU & deriveZIdxInCU
    UInt uiRasterZeroOrder =  g_auiZscanToRaster[deriveZIdxInCTU(uiAbsZeroIdx, uiAbsPartIdx)]; 
#else
    UInt uiRasterZeroOrder =  g_auiZscanToRaster[uiAbsPartIdx + uiAbsZeroIdx]; 
#endif
    UInt uiRasterOrder, uiSubBlkIdx;

    //get all candidates from luma block before pruning
    Int iDMNumFromLuma = 5;

    {
#if TT_TREE_DS5 // deriveZIdxInCTU
      //get TL
#if DIMD_PARSING_FIX_INTRA_JB
#if BARC_DIMD
      UInt uiIdxInBARCCTU = deriveZIdxInCTU(uiAbsZeroIdx, uiAbsPartIdx);
      if (pcCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiChMode[1] = bIntraSlice ? pcCtu->getSecIntraDir(uiIdxInBARCCTU) : getSecIntraDir(uiAbsPartIdx);
#else
      uiChMode[1] = bIntraSlice ? pcCtu->getSecIntraDir(deriveZIdxInCTU(uiAbsZeroIdx, uiAbsPartIdx)) : getSecIntraDir(uiAbsPartIdx);
#endif
#else
      uiChMode[1] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, deriveZIdxInCTU(uiAbsZeroIdx, uiAbsPartIdx));
#endif
      //get TR
      uiRasterOrder = uiRasterZeroOrder + (getWidth(uiAbsPartIdx) >>2 ) - 1;
      uiSubBlkIdx = deriveZIdxInCU(uiAbsZeroIdx, g_auiRasterToZscan[uiRasterOrder]);
#if DIMD_PARSING_FIX_INTRA_JB
#if BARC_DIMD
      uiIdxInBARCCTU = deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx);
      if (pcCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiChMode[2] = bIntraSlice ? pcCtu->getSecIntraDir(uiIdxInBARCCTU) : getSecIntraDir(uiSubBlkIdx);
#else
      uiChMode[2] = bIntraSlice ? pcCtu->getSecIntraDir(deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx)) : getSecIntraDir(uiSubBlkIdx);
#endif
#else
      uiChMode[2] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx));
#endif

      //get BL
      uiRasterOrder = uiRasterZeroOrder + getPic()->getNumPartInCtuWidth() *((getHeight(uiAbsPartIdx)>>2)-1);
      uiSubBlkIdx = deriveZIdxInCU(uiAbsZeroIdx, g_auiRasterToZscan[uiRasterOrder]);
#if DIMD_PARSING_FIX_INTRA_JB
#if BARC_DIMD
      uiIdxInBARCCTU = deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx);
      if (pcCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiChMode[3] = bIntraSlice ? pcCtu->getSecIntraDir(uiIdxInBARCCTU) : getSecIntraDir(uiSubBlkIdx);
#else
      uiChMode[3] = bIntraSlice ? pcCtu->getSecIntraDir(deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx)) : getSecIntraDir(uiSubBlkIdx);
#endif
#else
      uiChMode[3] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx));
#endif

      //get BR
      uiRasterOrder = uiRasterZeroOrder + getPic()->getNumPartInCtuWidth() *((getHeight(uiAbsPartIdx)>>2)-1) + (getWidth(uiAbsPartIdx)>>2)-1;

      uiSubBlkIdx = deriveZIdxInCU(uiAbsZeroIdx, g_auiRasterToZscan[uiRasterOrder]);
#if DIMD_PARSING_FIX_INTRA_JB
#if BARC_DIMD
      uiIdxInBARCCTU = deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx);
      if (pcCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiChMode[4] = bIntraSlice ? pcCtu->getSecIntraDir(uiIdxInBARCCTU) : getSecIntraDir(uiSubBlkIdx);
#else
      uiChMode[4] = bIntraSlice ? pcCtu->getSecIntraDir(deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx)) : getSecIntraDir(uiSubBlkIdx);
#endif
#else
      uiChMode[4] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, deriveZIdxInCTU(uiAbsZeroIdx, uiSubBlkIdx));
#endif
#else
#if DIMD_PARSING_FIX_INTRA_JB
      uiChMode[1] = bIntraSlice ? pcCtu->getSecIntraDir(uiAbsZeroIdx + uiAbsPartIdx) : getSecIntraDir(uiAbsPartIdx);
#else
      uiChMode[1] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsZeroIdx + uiAbsPartIdx);
#endif
      
      //get TR
      uiRasterOrder = uiRasterZeroOrder + (getWidth(uiAbsPartIdx) >>2 ) - 1;
      uiSubBlkIdx = g_auiRasterToZscan[uiRasterOrder] - uiAbsZeroIdx;
#if DIMD_PARSING_FIX_INTRA_JB
      uiChMode[2] = bIntraSlice ? pcCtu->getSecIntraDir(uiAbsZeroIdx + uiSubBlkIdx) : getSecIntraDir(uiAbsPartIdx);
#else
      uiChMode[2] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsZeroIdx + uiSubBlkIdx);
#endif

      //get BL
      uiRasterOrder = uiRasterZeroOrder + getPic()->getNumPartInCtuWidth() *((getHeight(uiAbsPartIdx)>>2)-1);
      uiSubBlkIdx = g_auiRasterToZscan[uiRasterOrder] - uiAbsZeroIdx;
#if DIMD_PARSING_FIX_INTRA_JB
      uiChMode[3] = bIntraSlice ? pcCtu->getSecIntraDir(uiAbsZeroIdx + uiSubBlkIdx) : getSecIntraDir(uiAbsPartIdx);
#else
      uiChMode[3] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsZeroIdx + uiSubBlkIdx);
#endif

      //get BR
      uiRasterOrder = uiRasterZeroOrder + getPic()->getNumPartInCtuWidth() *((getHeight(uiAbsPartIdx)>>2)-1) + (getWidth(uiAbsPartIdx)>>2)-1;

      uiSubBlkIdx = g_auiRasterToZscan[uiRasterOrder] - uiAbsZeroIdx;
#if DIMD_PARSING_FIX_INTRA_JB
      uiChMode[4] = bIntraSlice ? pcCtu->getSecIntraDir(uiAbsZeroIdx + uiSubBlkIdx) : getSecIntraDir(uiAbsPartIdx);
#else
      uiChMode[4] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsZeroIdx + uiSubBlkIdx);
#endif
#endif
      UInt uiWidth = getWidth(uiAbsPartIdx);
      UInt uiHeight = getHeight(uiAbsPartIdx);
      //get CR
      if(uiWidth==4 && uiHeight==4)
      {
        uiRasterOrder = uiRasterZeroOrder;
      }
      else if(uiHeight==4)
      {
        uiRasterOrder = uiRasterZeroOrder + (uiWidth>>2)/2 -1;
      }
      else if(uiWidth==4)
      {
        uiRasterOrder = uiRasterZeroOrder + pcCtu->getPic()->getNumPartInCtuWidth() *((uiHeight>>2)/2 - 1);
      }
      else
      {
        uiRasterOrder = uiRasterZeroOrder + pcCtu->getPic()->getNumPartInCtuWidth() *((uiHeight>>2)/2 - 1) + (uiWidth>>2)/2 -1;
      }
      uiSubBlkIdx = g_auiRasterToZscan[uiRasterOrder] - uiAbsZeroIdx;
#if DIMD_PARSING_FIX_INTRA_JB
#if BARC_DIMD
      uiIdxInBARCCTU = uiAbsZeroIdx + uiSubBlkIdx;
      if (pcCtu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiChMode[0] = bIntraSlice ? pcCtu->getSecIntraDir(uiIdxInBARCCTU) : getSecIntraDir(uiSubBlkIdx);
#else
      uiChMode[0] = bIntraSlice ? pcCtu->getSecIntraDir(uiAbsZeroIdx + uiSubBlkIdx) : getSecIntraDir(uiSubBlkIdx);
#endif
#else
      uiChMode[0] = pcCtu->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsZeroIdx + uiSubBlkIdx);
#endif
    }

    //perform pruning process
    Bool bRedundent[NUM_DM_MODES] = {false, false, false, false, false};
    for(i = 1; i < iDMNumFromLuma; i ++)
    {
      for(Int j = 0; j < i; j ++)
      {
        if(uiChMode[i] == uiChMode[j])
        {
          bRedundent[i] = true;
          break;
        }
      }
    }
    Int iAvailableCnt = 1;
    for(i = 1; i < iDMNumFromLuma; i ++)
    {
      if(bRedundent[i] == false)
      {
        uiChMode[iAvailableCnt++] = uiChMode[i];
      }
    } 
    if (iAvailableCnt < NUM_DM_MODES)
    {
      for (i = iAvailableCnt; i < iDMNumFromLuma; i++)
      {
        uiChMode[i] = INVALID_CHROMAMODE;
      }
      getIntraDirPredictor(uiAbsPartIdx, uiChMode, ComponentID(1)
#if MLIP_EXTEND_REF_LINE
          , 0
#endif
      );
    }
    if(uiChDMMode)
    {
      for(i = 0; i < NUM_DM_MODES; i ++)
      {
        uiChDMMode[i] = uiChMode[i];
      }
      (*iTotalCnt) = NUM_DM_MODES;
    } 
    uiChFinalMode = uiChMode[uiDMIdx];
  }
  else
#endif
  {
#if DIMD_PARSING_FIX_INTRA_JB  
    uiChFinalMode = getSecIntraDir(uiAbsPartIdx);
#else
    uiChFinalMode = getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
#endif
    uiChMode[0] = uiChFinalMode;
    getIntraDirPredictor(uiAbsPartIdx, uiChMode, ComponentID(1)
#if MLIP_EXTEND_REF_LINE
        , 0
#endif
    );
    if(uiChDMMode)
    {
      (*iTotalCnt) = NUM_DM_MODES;
      for(i = 0; i < NUM_DM_MODES; i ++)
      {
        uiChDMMode[i] = uiChMode[i];
      }
    }
    uiChFinalMode = uiChMode[uiDMIdx];
  }
  return uiChFinalMode;
}
#endif
/** Get allowed chroma intra modes
*   - fills uiModeList with chroma intra modes
*
*\param   [in]  uiAbsPartIdx
*\param   [out] uiModeList pointer to chroma intra modes array
*/

#if JVET_E0077_ENHANCED_LM
Int TComDataCU::getAllowedChromaDir(UInt uiAbsPartIdx, UInt uiModeList[NUM_CHROMA_MODE])
#else
Void TComDataCU::getAllowedChromaDir( UInt uiAbsPartIdx, UInt uiModeList[NUM_CHROMA_MODE] )
#endif
{  
#if JVET_E0062_MULTI_DMS
#if COM16_C806_LMCHROMA
  uiModeList[0] = LM_CHROMA_IDX;
  Int iStartIdx = 1;
  
#if JVET_E0077_ENHANCED_LM
  Int &iCurIdx = iStartIdx;

  const UInt csx = getComponentScaleX(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());
  const UInt csy = getComponentScaleY(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());

  Int iBlockSize = (getHeight(uiAbsPartIdx) >> csy) + (getWidth(uiAbsPartIdx) >> csx);

#if JVET_E0077_MMLM
  if (iBlockSize >= g_aiMMLM_MinSize[getSlice()->isIntra() ? 0 : 1])
  {
      uiModeList[iCurIdx++] = MMLM_CHROMA_IDX;
  }
#endif

#if JVET_E0077_LM_MF
  if (iBlockSize >= g_aiMFLM_MinSize[getSlice()->isIntra() ? 0 : 1])
  {
      for (Int s = 0; s < LM_FILTER_NUM; s++)
      {
          uiModeList[iCurIdx++] = LM_CHROMA_F1_IDX + s;
      }
  }
#endif
#endif

#else
  Int iStartIdx = 0;
#endif

  UInt uiAllowedDMmodes[ NUM_DM_MODES ];  
  UInt uiTotalDMmodes = 1;
  getDMMode(uiAbsPartIdx , 0, uiAllowedDMmodes, &uiTotalDMmodes); 

  for(Int iDMIdx = 0; iDMIdx < uiTotalDMmodes; iDMIdx ++)
  {
    uiModeList[iStartIdx++] = uiAllowedDMmodes[iDMIdx];
#if !JVET_E0077_ENHANCED_LM
    assert(uiModeList[iStartIdx-1]<=LM_CHROMA_IDX);
#endif
  } 
  assert(iStartIdx == NUM_CHROMA_MODE);
#if JVET_E0077_ENHANCED_LM
  return iStartIdx;
#endif
#else
  uiModeList[0] = PLANAR_IDX;
  uiModeList[1] = VER_IDX;
  uiModeList[2] = HOR_IDX;
  uiModeList[3] = DC_IDX;

#if COM16_C806_LMCHROMA
#if JVET_E0077_ENHANCED_LM
  Int iCurIdx = 4;

  const UInt csx = getComponentScaleX(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());
  const UInt csy = getComponentScaleY(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());

  Int iBlockSize = (getHeight(uiAbsPartIdx) >> csy) + (getWidth(uiAbsPartIdx) >> csx);

  uiModeList[iCurIdx++] = LM_CHROMA_IDX;
#if JVET_E0077_MMLM
  if (iBlockSize >= g_aiMMLM_MinSize[getSlice()->isIntra() ? 0 : 1])
  {
      uiModeList[iCurIdx++] = MMLM_CHROMA_IDX;
  }
#endif
#if JVET_E0077_LM_MF

  if (iBlockSize >= g_aiMFLM_MinSize[getSlice()->isIntra() ? 0 : 1])
  {
      for (Int s = 0; s < LM_FILTER_NUM; s++)
      {
          uiModeList[iCurIdx++] = LM_CHROMA_F1_IDX + s;
      }
  }
#endif

  uiModeList[iCurIdx++] = DM_CHROMA_IDX;
#else
  uiModeList[4] = LM_CHROMA_IDX;
  uiModeList[5] = DM_CHROMA_IDX;
#endif
  assert(5<NUM_CHROMA_MODE);
#else
  uiModeList[4] = DM_CHROMA_IDX;
  assert(4<NUM_CHROMA_MODE);
#endif

#if JVET_C0024_QTBT
#if DIMD_INTRA_PRED
  UInt uiLumaMode = getSlice()->isIntra() ? getPic()->getCtu(getCtuRsAddr())->getSecIntraDir(getZorderIdxInCtu() + uiAbsPartIdx) : getSecIntraDir(uiAbsPartIdx);
#else
  UInt uiLumaMode = getSlice()->isIntra() ? getPic()->getCtu(getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, getZorderIdxInCtu() + uiAbsPartIdx)
  : getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
#endif
#else
  UInt uiLumaMode = getIntraDir( CHANNEL_TYPE_LUMA, uiAbsPartIdx );
#endif

#if COM16_C806_LMCHROMA
  for( Int i = 0; i < NUM_CHROMA_MODE - 2; i++ )
#else
  for( Int i = 0; i < NUM_CHROMA_MODE - 1; i++ )
#endif
  {
    if( uiLumaMode == uiModeList[i] )
    {
#if VCEG_AZ07_INTRA_65ANG_MODES
      uiModeList[i] = VDIA_IDX; // VER+8 mode
#else
      uiModeList[i] = 34; // VER+8 mode
#endif
      break;
    }
  }

#if JVET_E0077_ENHANCED_LM
  return iCurIdx;
#endif
#endif
}



#if JVET_E0077_ENHANCED_LM

Int TComDataCU::getLMSymbolList(Int *pModeList, Int uiAbsPartIdx)
{
    const Int iNeibours = 5;
    TComDataCU* apcNeiboruingCU[iNeibours];
    UInt        auiNeiboruingPartIdx[iNeibours];

    UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;

    deriveLeftRightTopIdxGeneral(uiAbsPartIdx, 0, uiPartIdxLT, uiPartIdxRT);
    deriveLeftBottomIdxGeneral(uiAbsPartIdx, 0, uiPartIdxLB);

    apcNeiboruingCU[0] = getPULeft(auiNeiboruingPartIdx[0], uiPartIdxLB);
    apcNeiboruingCU[1] = getPUAbove(auiNeiboruingPartIdx[1], uiPartIdxRT);
    apcNeiboruingCU[2] = getPUAboveRight(auiNeiboruingPartIdx[2], uiPartIdxRT);
    apcNeiboruingCU[3] = getPUBelowLeft(auiNeiboruingPartIdx[3], uiPartIdxLB);
    apcNeiboruingCU[4] = getPUAboveLeft(auiNeiboruingPartIdx[4], m_absZIdxInCtu + uiAbsPartIdx);

    Int iCount = 0;

    for (Int i = 0; i < iNeibours; i++)
    {
        if (apcNeiboruingCU[i] && apcNeiboruingCU[i]->isIntra(auiNeiboruingPartIdx[i]))
        {
            Int iMode = apcNeiboruingCU[i]->getIntraDir(CHANNEL_TYPE_CHROMA, auiNeiboruingPartIdx[i]);
            if (!IsLMMode(iMode))
            {
                iCount++;
            }
        }
    }

    Bool bNonLMInsert = false;
    Int iIdx = 0;

    pModeList[iIdx++] = LM_CHROMA_IDX;

    if (iCount >= g_aiNonLMPosThrs[0] && !bNonLMInsert)
    {
        pModeList[iIdx++] = -1;
        bNonLMInsert = true;
    }

    const UInt csx = getComponentScaleX(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());
    const UInt csy = getComponentScaleY(COMPONENT_Cb, getSlice()->getSPS()->getChromaFormatIdc());

    Int iBlockSize = (getHeight(uiAbsPartIdx) >> csy) + (getWidth(uiAbsPartIdx) >> csx);


#if JVET_E0077_MMLM
    if (iBlockSize >= g_aiMMLM_MinSize[getSlice()->isIntra() ? 0 : 1])
    {
        pModeList[iIdx++] = MMLM_CHROMA_IDX;
    }
#endif

    if (iCount >= g_aiNonLMPosThrs[1] && !bNonLMInsert)
    {
        pModeList[iIdx++] = -1;
        bNonLMInsert = true;
    }

#if JVET_E0077_LM_MF

    if (iBlockSize >= g_aiMFLM_MinSize[getSlice()->isIntra() ? 0 : 1])
    {
        pModeList[iIdx++] = LM_CHROMA_F1_IDX;
        pModeList[iIdx++] = LM_CHROMA_F2_IDX;
        if (iCount >= g_aiNonLMPosThrs[2] && !bNonLMInsert)
        {
            pModeList[iIdx++] = -1;
            bNonLMInsert = true;
        }
        pModeList[iIdx++] = LM_CHROMA_F3_IDX;
        pModeList[iIdx++] = LM_CHROMA_F4_IDX;
    }
#endif


    if (!bNonLMInsert)
    {
        pModeList[iIdx++] = -1;
        bNonLMInsert = true;
    }

    return iIdx;
}

#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
TComDataCU* TComDataCU::getPULeftOffset( UInt& uiPartUnitIdx, 
                                        UInt uiCurrPartUnitIdx, 
                                        UInt uiPartOffset,
                                        Bool bEnforceSliceRestriction, 
                                        Bool bEnforceTileRestriction )
{
  UInt uiAbsPartIdx       = g_auiZscanToRaster[uiCurrPartUnitIdx];
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[m_absZIdxInCtu];
  UInt uiNumPartInCUWidth = m_pcPic->getNumPartInCtuWidth();
  
  uiAbsPartIdx += uiPartOffset * uiNumPartInCUWidth;

  if ( !RasterAddress::isZeroCol( uiAbsPartIdx, uiNumPartInCUWidth ) )
  {
    uiPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - 1 ];
    if ( RasterAddress::isEqualCol( uiAbsPartIdx, uiAbsZorderCUIdx, uiNumPartInCUWidth ) )
    {
      return m_pcPic->getCtu( getCtuRsAddr() );
    }
    else
    {
      uiPartUnitIdx -= m_absZIdxInCtu;
      return this;
    }
  }
  
  uiPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx + uiNumPartInCUWidth - 1 ];

  if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuLeft) )
  {
    return NULL;
  }
  return m_pCtuLeft;
}

TComDataCU* TComDataCU::getPUAboveOffset( UInt& uiPartUnitIdx, 
                                        UInt uiCurrPartUnitIdx, 
                                        UInt uiPartOffset,
                                        Bool bEnforceSliceRestriction, 
                                        Bool planarAtLCUBoundary,
                                        Bool bEnforceTileRestriction )
{
  UInt uiAbsPartIdx       = g_auiZscanToRaster[uiCurrPartUnitIdx];
  UInt uiAbsZorderCUIdx   = g_auiZscanToRaster[m_absZIdxInCtu];
  UInt uiNumPartInCUWidth = m_pcPic->getNumPartInCtuWidth();
  
  uiAbsPartIdx += uiPartOffset;

  if ( !RasterAddress::isZeroRow( uiAbsPartIdx, uiNumPartInCUWidth ) )
  {
    uiPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx - uiNumPartInCUWidth ];
    if ( RasterAddress::isEqualRow( uiAbsPartIdx, uiAbsZorderCUIdx, uiNumPartInCUWidth ) )
    {
      return m_pcPic->getCtu( getCtuRsAddr() );
    }
    else
    {
      uiPartUnitIdx -= m_absZIdxInCtu;
      return this;
    }
  }

  if(planarAtLCUBoundary)
  {
    return NULL;
  }
  
  uiPartUnitIdx = g_auiRasterToZscan[ uiAbsPartIdx + m_pcPic->getNumPartitionsInCtu() - uiNumPartInCUWidth ];

  if ( bEnforceSliceRestriction && !CUIsFromSameSliceAndTile(m_pCtuAbove) )
  {
    return NULL;
  }
  return m_pCtuAbove;
}
#endif

/** Get most probable intra modes
*\param   uiAbsPartIdx    partition index
*\param   uiIntraDirPred  pointer to the array for MPM storage
*\param   compID          colour component ID
*\param   piMode          it is set with MPM mode in case both MPM are equal. It is used to restrict RD search at encode side.
*\returns Number of MPM
*/
Void TComDataCU::getIntraDirPredictor( UInt uiAbsPartIdx, Int uiIntraDirPred[NUM_MOST_PROBABLE_MODES], const ComponentID compID
#if VCEG_AZ07_INTRA_65ANG_MODES && !JVET_C0055_INTRA_MPM
                                      , Int &iLeftAboveCase
#endif
#if MLIP_EXTEND_REF_LINE
    , Int ExternRef
#endif
                                      , Int* piMode
#if DIMD_INTRA_PRED 
  , Bool bUseRecIntraMode
#endif
)
{
#if JVET_C0055_INTRA_MPM
  // This function is not used for chroma texture type.
  // If it is used for chroma, DM mode should be converted to real mode during MPM list derivation.
#if !JVET_E0062_MULTI_DMS
  assert( isLuma(compID) );
#endif

  const ChannelType chType = toChannelType(compID);
  Bool includedMode[NUM_INTRA_MODE];
  memset( includedMode, false, sizeof(includedMode) );

  UInt modeIdx = 0, partIdx;

#if VCEG_AZ07_INTRA_65ANG_MODES
  const Int offset = 62;
#else
  const Int offset = 29;
#endif
  const Int mod = offset + 3;

  UInt width  = getWidth(uiAbsPartIdx);
  UInt height = getHeight(uiAbsPartIdx);

#if !JVET_C0024_QTBT
  if( getPartitionSize(uiAbsPartIdx) == SIZE_NxN )
  {
    width  >>= 1;
    height >>= 1;
  }
#endif

#if TT_TREE_DS5 // deriveZIdxInCTU
  UInt partIdxLT = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  UInt partIdxLT = m_absZIdxInCtu + uiAbsPartIdx;
#endif
  UInt partIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ partIdxLT ] + width / m_pcPic->getMinCUWidth() - 1 ];
  UInt partIdxLB = g_auiRasterToZscan [g_auiZscanToRaster[ partIdxLT ] + ( (height / m_pcPic->getMinCUHeight()) - 1 ) * m_pcPic->getNumPartInCtuWidth()];
#if JVET_E0062_MULTI_DMS
  Int iMaxModeNum = (compID == COMPONENT_Cb ? NUM_DM_MODES : NUM_MOST_PROBABLE_MODES);
  if(compID == COMPONENT_Cb)
  {
    Int i = 0;
    for(; i< NUM_DM_MODES; i++)
    {
      if (uiIntraDirPred[i] == INVALID_CHROMAMODE)
      {
        break;
      }
      else
      {
        includedMode[uiIntraDirPred[i]] = true;
      }
    }
#if COM16_C806_LMCHROMA
    includedMode[LM_CHROMA_IDX] = true;
#endif

#if JVET_E0077_MMLM
    includedMode[MMLM_CHROMA_IDX] = true;
#endif
#if JVET_E0077_LM_MF
    for (Int iF = 0; iF < LM_FILTER_NUM; iF++)
    {
        includedMode[LM_CHROMA_F1_IDX + iF] = true;
    }
#endif

    modeIdx = i;
  }
#endif
  // left
  TComDataCU *cu = getPULeft( partIdx, partIdxLB );
  if( cu && cu->isIntra(partIdx) )
  {
#if DIMD_INTRA_PRED
    if (bUseRecIntraMode)
    {
      uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
    }
    else
    {
#if BARC_DIMD
      UInt uiIdxInBARCCTU = partIdx;
      if (compID > 0 && cu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiIntraDirPred[modeIdx] = cu->getSecIntraDir(uiIdxInBARCCTU);
#else
      uiIntraDirPred[modeIdx] = cu->getSecIntraDir(partIdx);
#endif
    }
#else
    uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
#endif 
    if( !includedMode[uiIntraDirPred[modeIdx]] 
#if MLIP_EXTEND_REF_LINE
        && (compID == COMPONENT_Y ? (!(ExternRef && uiIntraDirPred[modeIdx] <= DC_IDX)) : 1)
#endif
        )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }   
  }
#if JVET_E0062_MULTI_DMS
  if(compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
  {
    return;
  }
#endif

  // above
  cu = getPUAbove( partIdx, partIdxRT );
  if( cu && cu->isIntra(partIdx) )
  {
#if DIMD_INTRA_PRED
    if (bUseRecIntraMode)
    {
      uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
    }
    else
    {
#if BARC_DIMD
      UInt uiIdxInBARCCTU = partIdx;
      if (compID > 0 && cu->getBARCFlag(CHANNEL_TYPE_LUMA))
        uiIdxInBARCCTU >>= 2;
      uiIntraDirPred[modeIdx] = cu->getSecIntraDir(uiIdxInBARCCTU);
#else
      uiIntraDirPred[modeIdx] = cu->getSecIntraDir(partIdx);
#endif
    }
#else
    uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
#endif
    if( !includedMode[uiIntraDirPred[modeIdx]]
#if MLIP_EXTEND_REF_LINE
        && (compID == COMPONENT_Y ? (!(ExternRef && uiIntraDirPred[modeIdx] <= DC_IDX)) : 1)
#endif
        )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }   
  }
#if JVET_E0062_MULTI_DMS
  if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
  {
    return;
  }
#endif

  if( piMode )
  {
    *piMode = Int(modeIdx > 1) + 1;
  }
#if JVET_E0062_MULTI_DMS  
  if(compID != COMPONENT_Cb)
  {
#endif  
  // PLANAR mode
#if MLIP_EXTEND_REF_LINE
      if (compID == COMPONENT_Y ? (ExternRef == 0) : 1)
      {
#endif

  uiIntraDirPred[modeIdx] = PLANAR_IDX;
  if( !includedMode[uiIntraDirPred[modeIdx]] )
  {
    includedMode[uiIntraDirPred[modeIdx]] = true;
    modeIdx++;
  }

  // DC mode
#if !VCEG_AZ07_INTRA_65ANG_MODES
  if( modeIdx < NUM_MOST_PROBABLE_MODES )
#endif
  {
    uiIntraDirPred[modeIdx] = DC_IDX;
    if( !includedMode[uiIntraDirPred[modeIdx]] )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }
  }
#if MLIP_EXTEND_REF_LINE
      }
#endif
#if JVET_E0062_MULTI_DMS  
  }
#endif 
  // below left
#if !VCEG_AZ07_INTRA_65ANG_MODES
  if( modeIdx < NUM_MOST_PROBABLE_MODES )
#endif
  {
    cu = getPUBelowLeft( partIdx, partIdxLB );
    if( cu && cu->isIntra(partIdx) )
    {
#if DIMD_INTRA_PRED
      if (bUseRecIntraMode)
      {
        uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
      }
      else
      {
#if BARC_DIMD
        UInt uiIdxInBARCCTU = partIdx;
        if (compID > 0 && cu->getBARCFlag(CHANNEL_TYPE_LUMA))
          uiIdxInBARCCTU >>= 2;
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(uiIdxInBARCCTU);
#else
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(partIdx);
#endif
      }
#else
      uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
#endif
      if( !includedMode[uiIntraDirPred[modeIdx]]
#if MLIP_EXTEND_REF_LINE
          && (compID == COMPONENT_Y ? (!(ExternRef && uiIntraDirPred[modeIdx] <= DC_IDX)) : 1)
#endif
          )
      {
        includedMode[uiIntraDirPred[modeIdx]] = true;
        modeIdx++;
      }
    }
  }
#if JVET_E0062_MULTI_DMS
  if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
  {
    return;
  }
#endif

  // above right
#if !VCEG_AZ07_INTRA_65ANG_MODES
  if( modeIdx < NUM_MOST_PROBABLE_MODES )
#endif
  {
    cu = getPUAboveRight( partIdx, partIdxRT );
    if( cu && cu->isIntra(partIdx) )
    {
#if DIMD_INTRA_PRED
      if (bUseRecIntraMode)
      {
        uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
      }
      else
      {
#if BARC_DIMD
        UInt uiIdxInBARCCTU = partIdx;
        if (compID > 0 && cu->getBARCFlag(CHANNEL_TYPE_LUMA))
          uiIdxInBARCCTU >>= 2;
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(uiIdxInBARCCTU);
#else
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(partIdx);
#endif
      }
#else
      uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
#endif
      if( !includedMode[uiIntraDirPred[modeIdx]]
#if MLIP_EXTEND_REF_LINE
          && (compID == COMPONENT_Y ? (!(ExternRef && uiIntraDirPred[modeIdx] <= DC_IDX)) : 1)
#endif
          )
      {
        includedMode[uiIntraDirPred[modeIdx]] = true;
        modeIdx++;
      }
    }
  }
#if JVET_E0062_MULTI_DMS
  if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
  {
    return;
  }
#endif

  //above left
#if JVET_E0062_MULTI_DMS  
  if( modeIdx < iMaxModeNum )
#else
  if( modeIdx < NUM_MOST_PROBABLE_MODES )
#endif
  {
    cu = getPUAboveLeft( partIdx, partIdxLT );
    if( cu && cu->isIntra(partIdx) )
    {
#if DIMD_INTRA_PRED
      if (bUseRecIntraMode)
      {
        uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
      }
      else
      {
#if BARC_DIMD
        UInt uiIdxInBARCCTU = partIdx;
        if (compID > 0 && cu->getBARCFlag(CHANNEL_TYPE_LUMA))
          uiIdxInBARCCTU >>= 2;
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(uiIdxInBARCCTU);
#else
        uiIntraDirPred[modeIdx] = cu->getSecIntraDir(partIdx);
#endif
      }
#else
      uiIntraDirPred[modeIdx] = cu->getIntraDir(chType, partIdx);
#endif
      if( !includedMode[uiIntraDirPred[modeIdx]]
#if MLIP_EXTEND_REF_LINE
          && (compID == COMPONENT_Y ? (!(ExternRef && uiIntraDirPred[modeIdx] <= DC_IDX)) : 1)
#endif
          )
      {
        includedMode[uiIntraDirPred[modeIdx]] = true;
        modeIdx++;
      }
    }
  }
#if JVET_E0062_MULTI_DMS
  if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
  {
    return;
  }
 
  if(compID == COMPONENT_Cb)
  {
    // PLANAR mode
    uiIntraDirPred[modeIdx] = PLANAR_IDX;
    if( !includedMode[uiIntraDirPred[modeIdx]] )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }

    if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
    {
      return;
    }

    // DC mode
#if !VCEG_AZ07_INTRA_65ANG_MODES
    if( modeIdx < NUM_MOST_PROBABLE_MODES )
#endif
    {
      uiIntraDirPred[modeIdx] = DC_IDX;
      if( !includedMode[uiIntraDirPred[modeIdx]] )
      {
        includedMode[uiIntraDirPred[modeIdx]] = true;
        modeIdx++;
      }
    }

    if( compID == COMPONENT_Cb && modeIdx == NUM_DM_MODES )
    {
      return;
    }
  }
#endif 
  UInt numAddedModes = modeIdx;

  // -+1 derived modes
#if JVET_E0062_MULTI_DMS  
  for( UInt idx = 0; idx < numAddedModes && modeIdx < iMaxModeNum; idx++ )
#else
  for( UInt idx = 0; idx < numAddedModes && modeIdx < NUM_MOST_PROBABLE_MODES; idx++ )
#endif
  {   
    UInt mode = uiIntraDirPred[idx];
    
    // directional mode
    if( mode <= DC_IDX )
    {
      continue;
    }

    // -1
    uiIntraDirPred[modeIdx] = ((mode + offset) % mod) + 2;

    if( !includedMode[uiIntraDirPred[modeIdx]] )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }    

#if JVET_E0062_MULTI_DMS
    if( modeIdx == iMaxModeNum )
#else
    if( modeIdx == NUM_MOST_PROBABLE_MODES )
#endif
    {
      break;
    }

    // +1
    uiIntraDirPred[modeIdx] = ((mode - 1) % mod) + 2;

    if( !includedMode[uiIntraDirPred[modeIdx]] )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }   
  }

  // default modes
  UInt defaultIntraModes[] = {PLANAR_IDX, DC_IDX, VER_IDX, HOR_IDX, 2, DIA_IDX};
#if MLIP_EXTEND_REF_LINE
  if (compID == COMPONENT_Y ? (ExternRef == 0) : 1)
#endif
  assert( modeIdx > 1 );
#if JVET_E0062_MULTI_DMS  
  for( UInt idx = 2; idx < iMaxModeNum && modeIdx < iMaxModeNum; idx++)
#else
  for( UInt idx = 2; idx < NUM_MOST_PROBABLE_MODES && modeIdx < NUM_MOST_PROBABLE_MODES; idx++ )
#endif
  {
    uiIntraDirPred[modeIdx] = defaultIntraModes[idx];
    if( !includedMode[uiIntraDirPred[modeIdx]] )
    {
      includedMode[uiIntraDirPred[modeIdx]] = true;
      modeIdx++;
    }
  }

#if MLIP_EXTEND_REF_LINE
  if ((compID == COMPONENT_Y) && ExternRef && modeIdx < NUM_MOST_PROBABLE_MODES) {
      UInt mode[] = { 10, 66, 26, 42, 30, 38 };
      for (UInt i = 0; i < NUM_MOST_PROBABLE_MODES && modeIdx < NUM_MOST_PROBABLE_MODES; i++) {
          uiIntraDirPred[modeIdx] = mode[i];
          if (!includedMode[uiIntraDirPred[modeIdx]])
          {
              includedMode[uiIntraDirPred[modeIdx]] = true;
              modeIdx++;
          }
      }
  }

#endif

#if JVET_E0062_MULTI_DMS  
  assert( modeIdx == iMaxModeNum );
#else
  assert( modeIdx == NUM_MOST_PROBABLE_MODES );
#endif
#else
  TComDataCU* pcCULeft = NULL, *pcCUAbove = NULL;
  UInt        LeftPartIdx  = MAX_UINT;
  UInt        AbovePartIdx = MAX_UINT;
  Int         iLeftIntraDir, iAboveIntraDir;
#if !VCEG_AZ07_INTRA_65ANG_MODES && !JVET_C0024_QTBT
  const TComSPS *sps=getSlice()->getSPS();
  const UInt partsPerMinCU = 1<<(2*(sps->getMaxTotalCUDepth() - sps->getLog2DiffMaxMinCodingBlockSize()));
#endif

  const ChannelType chType = toChannelType(compID);
#if !VCEG_AZ07_INTRA_65ANG_MODES && !JVET_C0024_QTBT
  const ChromaFormat chForm = getPic()->getChromaFormat();
#endif

#if VCEG_AZ07_INTRA_65ANG_MODES
#if JVET_C0024_QTBT
  UInt        uiPUWidth   = getWidth(uiAbsPartIdx) ;
  UInt        uiPUHeight  = getHeight(uiAbsPartIdx) ;
#else
  UInt        uiPUWidth   = ( getWidth(uiAbsPartIdx) >> ( getPartitionSize(uiAbsPartIdx)==SIZE_NxN ? 1 : 0 ) );
  UInt        uiPUHeight  = ( getHeight(uiAbsPartIdx) >> ( getPartitionSize(uiAbsPartIdx)==SIZE_NxN ? 1 : 0 ) );
#endif

  UInt        uiCaseIdx   = 0;

  UInt        uiPartIdxLT = m_absZIdxInCtu + uiAbsPartIdx;
  UInt        uiPartIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ uiPartIdxLT ] + uiPUWidth / m_pcPic->getMinCUWidth() - 1 ];
  UInt        uiPartIdxLB = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu + uiAbsPartIdx ] + ((uiPUHeight / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth()];;
  UChar       ucTempIntraDir;
#if JVET_C0024_QTBT
  const UInt  uiNumUnitsInPUHeight = (g_auiZscanToRaster[uiPartIdxLB] - g_auiZscanToRaster[uiPartIdxLT]) / getPic()->getNumPartInCtuWidth() + 1;
  const UInt  uiNumUnitsInPUWidth = g_auiZscanToRaster[uiPartIdxRT] - g_auiZscanToRaster[uiPartIdxLT] + 1;
#else
  const UInt  uiNumUnitsInPU = (g_auiZscanToRaster[uiPartIdxLB] - g_auiZscanToRaster[uiPartIdxLT]) / getPic()->getNumPartInCtuWidth() + 1;
#endif

  // Initializaiton
  Int uiLeftIntraDirCnt[NUM_INTRA_MODE], uiAboveIntraDirCnt[NUM_INTRA_MODE];
  memset( uiLeftIntraDirCnt,  0, sizeof(Int)*NUM_INTRA_MODE );
  memset( uiAboveIntraDirCnt, 0, sizeof(Int)*NUM_INTRA_MODE );
  iLeftIntraDir  = DC_IDX;
  iAboveIntraDir = DC_IDX;

#if !JVET_C0024_QTBT
  deriveLeftRightTopIdxGeneral( uiAbsPartIdx, 0, uiPartIdxLT, uiPartIdxRT );
#endif

  // Get intra direction from above side
  UInt uiDirMaxCount = 0;
#if JVET_C0024_QTBT
  for ( UInt uiOffset = 0; uiOffset < uiNumUnitsInPUWidth; uiOffset++ )
#else
  for ( UInt uiOffset = 0; uiOffset < uiNumUnitsInPU; uiOffset++ )
#endif
  {
    pcCUAbove = getPUAboveOffset( AbovePartIdx, uiPartIdxLT, uiOffset );

    if( pcCUAbove && pcCUAbove->isIntra( AbovePartIdx ) == MODE_INTRA )
    {
#if DIMD_INTRA_PRED
      if (bUseRecIntraMode)
      {
        ucTempIntraDir = pcCUAbove->getIntraDir(chType, AbovePartIdx);
      }
      else
      {
        ucTempIntraDir = pcCUAbove->getSecIntraDir(AbovePartIdx);
      }
#else
      ucTempIntraDir = pcCUAbove->getIntraDir(chType, AbovePartIdx);
#endif
      uiAboveIntraDirCnt[ucTempIntraDir] ++;
      if( uiDirMaxCount<uiAboveIntraDirCnt[ucTempIntraDir] )
      {
        uiDirMaxCount  = uiAboveIntraDirCnt[ucTempIntraDir];
        iAboveIntraDir = ucTempIntraDir;
      }
    }
  }

  // Get intra direction from left side
  uiDirMaxCount = 0;
#if JVET_C0024_QTBT
  for ( UInt uiOffset = 0; uiOffset < uiNumUnitsInPUHeight; uiOffset++ )
#else
  for ( UInt uiOffset = 0; uiOffset < uiNumUnitsInPU; uiOffset++ )
#endif
  {
    pcCULeft = getPULeftOffset( LeftPartIdx, uiPartIdxLT, uiOffset );

    if( pcCULeft && pcCULeft->isIntra( LeftPartIdx ) == MODE_INTRA )
    {
#if DIMD_INTRA_PRED
      assert(chType == CHANNEL_TYPE_LUMA);
      if (bUseRecIntraMode)
      {
        ucTempIntraDir = pcCULeft->getIntraDir(chType, LeftPartIdx);
      }
      else
      {
        ucTempIntraDir = pcCULeft->getSecIntraDir(LeftPartIdx);
      }
#else
      ucTempIntraDir = pcCULeft->getIntraDir(chType, LeftPartIdx);
#endif
      uiLeftIntraDirCnt[ucTempIntraDir] ++;
      if( uiDirMaxCount<uiLeftIntraDirCnt[ucTempIntraDir] )
      {
        uiDirMaxCount = uiLeftIntraDirCnt[ucTempIntraDir];
        iLeftIntraDir = ucTempIntraDir;
      }
    }
  }
#else
  // Get intra direction of left PU
  pcCULeft = getPULeft( LeftPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#if !JVET_C0024_QTBT
  if (isChroma(compID))
  {
    LeftPartIdx = getChromasCorrespondingPULumaIdx(LeftPartIdx, chForm, partsPerMinCU);
  }
#endif
  iLeftIntraDir  = pcCULeft ? ( pcCULeft->isIntra( LeftPartIdx ) ? pcCULeft->getIntraDir( chType, LeftPartIdx ) : DC_IDX ) : DC_IDX;

  // Get intra direction of above PU
  pcCUAbove = getPUAbove( AbovePartIdx, m_absZIdxInCtu + uiAbsPartIdx, true, true );
#if !JVET_C0024_QTBT
  if (isChroma(compID))
  {
    AbovePartIdx = getChromasCorrespondingPULumaIdx(AbovePartIdx, chForm, partsPerMinCU);
  }
#endif
  iAboveIntraDir = pcCUAbove ? ( pcCUAbove->isIntra( AbovePartIdx ) ? pcCUAbove->getIntraDir( chType, AbovePartIdx ) : DC_IDX ) : DC_IDX;
#endif

  if (isChroma(chType))
  {
    if (iLeftIntraDir  == DM_CHROMA_IDX)
    {
      iLeftIntraDir  = pcCULeft-> getIntraDir( CHANNEL_TYPE_LUMA, LeftPartIdx  );
    }
    if (iAboveIntraDir == DM_CHROMA_IDX)
    {
      iAboveIntraDir = pcCUAbove->getIntraDir( CHANNEL_TYPE_LUMA, AbovePartIdx );
    }
  }

#if VCEG_AZ07_INTRA_65ANG_MODES
  const UInt uiHor = HOR_IDX;
  const UInt uiVer = VER_IDX;
  const UInt uiDia = DIA_IDX;
  const Int  iOffset = 62;
  const Int  iMod = iOffset+3;
#endif

  assert (2<NUM_MOST_PROBABLE_MODES);
  if(iLeftIntraDir == iAboveIntraDir)
  {
    if( piMode )
    {
      *piMode = 1;
    }

    if (iLeftIntraDir > 1) // angular modes
    {
      uiIntraDirPred[0] = iLeftIntraDir;
#if VCEG_AZ07_INTRA_65ANG_MODES
      uiIntraDirPred[1] = PLANAR_IDX;
      uiIntraDirPred[2] = ((iLeftIntraDir - 1 ) % iMod) + 2;
      uiIntraDirPred[3] = ((iLeftIntraDir + iOffset) % iMod) + 2;
      uiIntraDirPred[4] = ((uiIntraDirPred[2] - 1 ) % iMod) + 2;
      uiIntraDirPred[5] = DC_IDX;
      uiCaseIdx         = 0;
#else
      uiIntraDirPred[1] = ((iLeftIntraDir + 29) % 32) + 2;
      uiIntraDirPred[2] = ((iLeftIntraDir - 1 ) % 32) + 2;
#endif
    }
    else //non-angular
    {
      uiIntraDirPred[0] = PLANAR_IDX;
      uiIntraDirPred[1] = DC_IDX;
#if VCEG_AZ07_INTRA_65ANG_MODES
      uiIntraDirPred[2] = uiVer; 
      uiIntraDirPred[3] = uiHor; 
      uiIntraDirPred[4] = 2;
      uiIntraDirPred[5] = uiDia;
      uiCaseIdx         = 1;
#else
      uiIntraDirPred[2] = VER_IDX;
#endif
    }
  }
  else
  {
    if( piMode )
    {
      *piMode = 2;
    }
    uiIntraDirPred[0] = iLeftIntraDir;
    uiIntraDirPred[1] = iAboveIntraDir;

    if (iLeftIntraDir && iAboveIntraDir ) //both modes are non-planar
    {
      uiIntraDirPred[2] = PLANAR_IDX;
#if VCEG_AZ07_INTRA_65ANG_MODES
      Int iMaxDir = max( iAboveIntraDir, iLeftIntraDir );
      Int iMinDir = min( iAboveIntraDir, iLeftIntraDir );
      if( iLeftIntraDir==DC_IDX || iAboveIntraDir==DC_IDX )
      {
        Int iNonDcMode = iMaxDir;
        uiIntraDirPred[3] = (( iNonDcMode + iOffset ) % iMod) + 2;
        uiIntraDirPred[4] = (( iNonDcMode - 1 ) % iMod) + 2;
        uiIntraDirPred[5] = (( uiIntraDirPred[4] -1 ) % iMod) + 2;
      }
      else
      {
        uiIntraDirPred[3] = DC_IDX;
        uiIntraDirPred[4] = (( iMaxDir - 1 ) % iMod) + 2;
        if( uiIntraDirPred[4] == iMinDir )
        {
          uiIntraDirPred[4]++;
        }
        uiIntraDirPred[5] = (( iMinDir + iOffset ) % iMod) + 2;
        if( uiIntraDirPred[5] == iMaxDir )
        {
          uiIntraDirPred[5]--;
        }
        if( uiIntraDirPred[5] == uiIntraDirPred[4] )
        {
          uiIntraDirPred[5] = iMinDir + 1;
        }
      }
      uiCaseIdx = 2;
#endif
    }
    else
    {
#if VCEG_AZ07_INTRA_65ANG_MODES
      Int iMaxDir = max( iAboveIntraDir, iLeftIntraDir );
      uiIntraDirPred[2] =  (iLeftIntraDir+iAboveIntraDir)<2? uiVer : DC_IDX;
      if( (iLeftIntraDir+iAboveIntraDir)<2 )
      {
        uiIntraDirPred[3] = uiHor;
        uiIntraDirPred[4] = 2;
        uiIntraDirPred[5] = uiDia;
      }
      else
      {
        uiIntraDirPred[3] = ((iMaxDir + iOffset ) % iMod) + 2;
        uiIntraDirPred[4] = ((iMaxDir - 1 ) % iMod) + 2;
        uiIntraDirPred[5] = ((uiIntraDirPred[4] - 1 ) % iMod) + 2;
      }
      uiCaseIdx = 3;
#else
      uiIntraDirPred[2] =  (iLeftIntraDir+iAboveIntraDir)<2? VER_IDX : DC_IDX;
#endif
    }
  }
#if VCEG_AZ07_INTRA_65ANG_MODES
  if ( iLeftAboveCase!=-1 )
  {
    iLeftAboveCase = uiCaseIdx;
  }
#endif
#endif

#if JVET_E0062_MULTI_DMS
  for (UInt i=0; i<iMaxModeNum; i++)
#else
  for (UInt i=0; i<NUM_MOST_PROBABLE_MODES; i++)
#endif
  {
#if VCEG_AZ07_INTRA_65ANG_MODES
    assert(uiIntraDirPred[i] < (NUM_INTRA_MODE - 1));

#else
    assert(uiIntraDirPred[i] < 35);
#endif
  }
}

#if JVET_C0024_QTBT
#if TT_SYNTAX_NEW
Char TComDataCU::deriveAllowedSplitType( UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, UInt uiBTDepth )
{
  bool bHBT = true, bHTT = true, bVBT = true, bVTT = true;
  UInt uiMinBTSize = getSlice()->isIntra() ? (isLuma(getTextType())?MIN_BT_SIZE:MIN_BT_SIZE_C): MIN_BT_SIZE_INTER;

  if(uiWidth==uiMinBTSize)
    bVBT = bVTT = false;
  else if(uiWidth==(uiMinBTSize<<1))
    bVTT = false;
  if(uiHeight==uiMinBTSize)
    bHBT = bHTT = false;
  else if(uiHeight==(uiMinBTSize<<1))
    bHTT = false;

#if TT_LIMIT_SYNTAX
  UInt uiParentType = (uiBTDepth==0)? 0 : getBTSplitModeForBTDepth(uiAbsPartIdx, uiBTDepth-1);
#if TT_LIMIT_SMALL
  if(uiParentType==ST_HORI_TT_S)
    bHTT = false;
  else if(uiParentType==ST_VERT_TT_S)
    bVTT = false;
#endif
#if TT_LIMIT_LARGE
  if(uiParentType==ST_HORI_TT_L)
    bHTT = false;
  else if(uiParentType==ST_VERT_TT_L)
    bVTT = false;
#endif
#endif

#if (TT_LIMIT_AR>0)
  if(uiWidth>=uiHeight*TT_LIMIT_AR)
    bHBT = bHTT = false;
  else if(uiHeight>=uiWidth*TT_LIMIT_AR)
    bVBT = bVTT = false;
#endif

#if JVET_C0024_BT_RMV_REDUNDANT
  if ( getSplitConstrain() == 1 ) 
    bHBT = false;
  else if( getSplitConstrain() == 2 )
    bVBT = false;
#endif

  return ((int)bHBT<<3) + ((int)bVBT<<2) + ((int)bHTT<<1) + (int)bVTT;
}
#endif
#if CNN_BASED_INLOOP_FILTER
UInt TComDataCU::getCNNLoopFilterCtxFlag(ComponentID compIdx, CNNLFBlkParam* &ctuParam)
{
	TComDataCU* pcTempCU;
	UInt        uiTempPartIdx;
	UInt        uiCtx = 0;

	pcTempCU = getPULeft(uiTempPartIdx, 0);
	uiCtx = (pcTempCU) ? UInt(ctuParam[pcTempCU->getCtuRsAddr()][compIdx]) : 0;
	pcTempCU = getPUAbove(uiTempPartIdx, 0);
	uiCtx += (pcTempCU) ? UInt(ctuParam[pcTempCU->getCtuRsAddr()][compIdx]) : 0;
	return uiCtx;
}
#endif
UInt TComDataCU::getCtxBTSplitFlag( UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx;
  UInt        uiTempDepth;

#if TT_SYNTAX //split ctx
  UInt uiCurrDepth = getDepth(uiAbsPartIdx)*2 + getDeltaSizeIdx(uiAbsPartIdx, uiWidth, uiHeight);
#else
  UInt uiCurrDepth = getDepth(uiAbsPartIdx)*2 + getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);
#endif

  // Get left split flag
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
#if TT_SYNTAX //split ctx
  uiTempDepth = ( pcTempCU ) ? (pcTempCU->getDepth(uiTempPartIdx)*2 + pcTempCU->getDeltaSizeIdx(uiTempPartIdx)) : 0;
#else
  uiTempDepth = ( pcTempCU ) ? (pcTempCU->getDepth(uiTempPartIdx)*2 + pcTempCU->getBTDepth(uiTempPartIdx)) : 0;
#endif

  uiCtx  = ( pcTempCU ) ? (  uiTempDepth > uiCurrDepth ? 1 : 0 ) : 0;

  // Get above split flag
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
#if TT_SYNTAX //split ctx
  uiTempDepth = ( pcTempCU ) ? (pcTempCU->getDepth(uiTempPartIdx)*2 + pcTempCU->getDeltaSizeIdx(uiTempPartIdx)) : 0;
#else
  uiTempDepth = ( pcTempCU ) ? (pcTempCU->getDepth(uiTempPartIdx)*2 + pcTempCU->getBTDepth(uiTempPartIdx)) : 0;
#endif

  uiCtx  += ( pcTempCU ) ? (  uiTempDepth > uiCurrDepth ? 1 : 0 ) : 0;

  return uiCtx;
}
#endif

UInt TComDataCU::getCtxSplitFlag( UInt uiAbsPartIdx, UInt uiDepth )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx;
  // Get left split flag
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx  = ( pcTempCU ) ? ( ( pcTempCU->getDepth( uiTempPartIdx ) > uiDepth ) ? 1 : 0 ) : 0;

  // Get above split flag
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx += ( pcTempCU ) ? ( ( pcTempCU->getDepth( uiTempPartIdx ) > uiDepth ) ? 1 : 0 ) : 0;

#if COM16_C806_LARGE_CTU
  UChar ucMinDepth = 0;
#if JVET_C0024_QTBT
  UChar ucMaxDepth = g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] 
  - g_aucConvertToBit[getSlice()->getSPS()->getMinQTSize(getSlice()->getSliceType(), getTextType())];
#else
  UChar ucMaxDepth = ( UChar )( getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() );
#endif

#if TT_TREE_DS1 // deriveZIdxInCTU
  getMaxMinCUDepth( ucMinDepth , ucMaxDepth , deriveZIdxInCTU(getZorderIdxInCtu(),uiAbsPartIdx) );
#else
  getMaxMinCUDepth( ucMinDepth , ucMaxDepth , uiAbsPartIdx + getZorderIdxInCtu() );
#endif
  if( uiDepth < ucMinDepth )
  {
    uiCtx = 3;
  }
  else if( uiDepth >= ucMaxDepth + 1 )
  {    
    uiCtx = 4;
  }
#endif

  return uiCtx;
}

UInt TComDataCU::getCtxQtCbf( TComTU &rTu, const ChannelType chType )
{
  const UInt transformDepth = rTu.GetTransformDepthRel();

  if (isChroma(chType))
  {
    return transformDepth;
  }
  else
  {
    const UInt uiCtx = ( transformDepth == 0 ? 1 : 0 );
    return uiCtx;
  }
}

#if !JVET_C0024_QTBT
UInt TComDataCU::getQuadtreeTULog2MinSizeInCU( UInt absPartIdx )
{
  UInt log2CbSize = g_aucConvertToBit[getWidth( absPartIdx )] + 2;
  PartSize  partSize  = getPartitionSize( absPartIdx );
  UInt quadtreeTUMaxDepth = isIntra( absPartIdx ) ? m_pcSlice->getSPS()->getQuadtreeTUMaxDepthIntra() : m_pcSlice->getSPS()->getQuadtreeTUMaxDepthInter();
  Int intraSplitFlag = ( isIntra( absPartIdx ) && partSize == SIZE_NxN ) ? 1 : 0;
  Int interSplitFlag = ((quadtreeTUMaxDepth == 1) && isInter( absPartIdx ) && (partSize != SIZE_2Nx2N) );

  UInt log2MinTUSizeInCU = 0;
  if (log2CbSize < (m_pcSlice->getSPS()->getQuadtreeTULog2MinSize() + quadtreeTUMaxDepth - 1 + interSplitFlag + intraSplitFlag) )
  {
    // when fully making use of signaled TUMaxDepth + inter/intraSplitFlag, resulting luma TB size is < QuadtreeTULog2MinSize
    log2MinTUSizeInCU = m_pcSlice->getSPS()->getQuadtreeTULog2MinSize();
  }
  else
  {
    // when fully making use of signaled TUMaxDepth + inter/intraSplitFlag, resulting luma TB size is still >= QuadtreeTULog2MinSize
    log2MinTUSizeInCU = log2CbSize - ( quadtreeTUMaxDepth - 1 + interSplitFlag + intraSplitFlag); // stop when trafoDepth == hierarchy_depth = splitFlag
    if ( log2MinTUSizeInCU > m_pcSlice->getSPS()->getQuadtreeTULog2MaxSize())
    {
      // when fully making use of signaled TUMaxDepth + inter/intraSplitFlag, resulting luma TB size is still > QuadtreeTULog2MaxSize
      log2MinTUSizeInCU = m_pcSlice->getSPS()->getQuadtreeTULog2MaxSize();
    }
  }
  return log2MinTUSizeInCU;
}
#endif

UInt TComDataCU::getCtxSkipFlag( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

  // Get BCBP of left PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx    = ( pcTempCU ) ? pcTempCU->isSkipped( uiTempPartIdx ) : 0;

  // Get BCBP of above PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx   += ( pcTempCU ) ? pcTempCU->isSkipped( uiTempPartIdx ) : 0;

  return uiCtx;
}

#if BARC_UP_SAMPLE_RDO
UInt TComDataCU::getCtxUpCNNFlag(UInt uiAbsPartIdx)
{
	TComDataCU* pcTempCU;
	UInt        uiTempPartIdx;
	UInt        uiCtx = 0;

	pcTempCU = getPULeft(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
	uiCtx = (pcTempCU && pcTempCU->getBARCFlag(pcTempCU->getTextType())) ? pcTempCU->getUpCNNFlag() : 0;

	pcTempCU = getPUAbove(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
	uiCtx += (pcTempCU && pcTempCU->getBARCFlag(pcTempCU->getTextType())) ? pcTempCU->getUpCNNFlag() : 0;

	return uiCtx;
}
#endif
#if BARC
UInt TComDataCU::getCtxBARCFlag( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;
  // Get BCBP of left PU
  pcTempCU = getPULeft( uiTempPartIdx, 0 );

  uiCtx    = ( pcTempCU ) ? pcTempCU->getBARCFlag( pcTempCU->getTextType() ) : 0;

  return uiCtx;
}
#endif

#if VCEG_AZ07_IMV
UInt TComDataCU::getCtxiMVFlag( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

  // Get BCBP of left PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif

#if  JVET_E0076_MULTI_PEL_MVD
  uiCtx    = ( pcTempCU ) ? (pcTempCU->getiMVFlag( uiTempPartIdx ) ? 1 : 0 ) : 0;
#else
  uiCtx    = ( pcTempCU ) ? pcTempCU->getiMVFlag( uiTempPartIdx ) : 0;
#endif

  // Get BCBP of above PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif

#if  JVET_E0076_MULTI_PEL_MVD
  uiCtx   += ( pcTempCU ) ? (pcTempCU->getiMVFlag( uiTempPartIdx ) ? 1 : 0 ) : 0;
#else
  uiCtx   += ( pcTempCU ) ? pcTempCU->getiMVFlag( uiTempPartIdx ) : 0;
#endif

  return uiCtx;
}

Bool TComDataCU::hasSubCUNonZeroMVd()
{
  Bool bNonZeroMvd = false;
#if JVET_C0024_QTBT
  Int iPartNum = 1;
#else
  Int iPartNum = ( getPartitionSize( 0 ) == SIZE_NxN ? 4 : ( getPartitionSize( 0 ) == SIZE_2Nx2N ? 1 : 2 ) );
#endif
  for( Int iPart = 0; iPart < iPartNum; iPart++ )
  {
    UInt uiPartAddr = 0;
    Int iHeight, iWidth;
    getPartIndexAndSize( iPart, uiPartAddr, iWidth, iHeight );
    for( Int iRefPicList = 0; iRefPicList < 2; iRefPicList++ )
    {
      RefPicList eRefPicList = (RefPicList) iRefPicList;
      if( getCUMvField( eRefPicList )->getRefIdx( uiPartAddr ) >= 0 )
      {
        bNonZeroMvd |= ( getCUMvField( eRefPicList )->getMvd( uiPartAddr ).getHor() != 0 );
        bNonZeroMvd |= ( getCUMvField( eRefPicList )->getMvd( uiPartAddr ).getVer() != 0 );
      }
    }
  }
  return( bNonZeroMvd );
}

Char TComDataCU::getMaxNeighboriMVCandNum( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  const Char  cDefault = 0;
  Char        cMaxiMVCandNum = 0;

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  cMaxiMVCandNum    = max( cMaxiMVCandNum , ( pcTempCU ) ? pcTempCU->getiMVCandNum( uiTempPartIdx ) : cDefault );

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  cMaxiMVCandNum    = max( cMaxiMVCandNum , ( pcTempCU ) ? pcTempCU->getiMVCandNum( uiTempPartIdx ) : cDefault );

  return cMaxiMVCandNum;
}

Bool TComDataCU::resetMVDandMV2Int( UInt uiAbsPartIdx , UInt uiPartIdx , Bool bResetMV 
#if VCEG_AZ07_FRUC_MERGE
  , TComPrediction * pPred
#endif
)
{
#if JVET_E0076_MULTI_PEL_MVD
  assert( getiMVFlag( uiAbsPartIdx ) != 0 );
#else
  assert( getiMVFlag( uiAbsPartIdx ) == true );
#endif

  if( !getMergeFlag( uiAbsPartIdx ) )
  {
    for( RefPicList eRefPicList = REF_PIC_LIST_0 ; eRefPicList <= REF_PIC_LIST_1 ; eRefPicList = ( RefPicList )( ( Int )eRefPicList + 1 ) )
    {
      TComCUMvField * pMVField = getCUMvField( eRefPicList );
      if( pMVField->getRefIdx( uiAbsPartIdx ) >= 0 )
      {
        TComMv mv = pMVField->getMv( uiAbsPartIdx );
        TComMv mvPred;

        AMVPInfo*  pcAMVPInfo = getCUMvField( eRefPicList )->getAMVPInfo();
        fillMvpCand( uiPartIdx , uiAbsPartIdx , eRefPicList , pMVField->getRefIdx( uiAbsPartIdx ) , pcAMVPInfo 
#if VCEG_AZ07_FRUC_MERGE
          , pPred
#endif
          );

        mvPred = pcAMVPInfo->m_acMvCand[getMVPIdx(eRefPicList , uiAbsPartIdx)];

        if( bResetMV )
        {
#if JVET_E0076_MULTI_PEL_MVD
          xRoundMV( mv, uiAbsPartIdx );
#else
          xRoundMV( mv );
#endif
        }

        TComMv mvDiff = mv - mvPred;
        if( getSlice()->getMvdL1ZeroFlag() && eRefPicList == REF_PIC_LIST_1 && getInterDir( uiAbsPartIdx ) == 3 )
        {
          mvDiff.setZero();
          assert( bResetMV == true || mv == mvPred );
          mv = mvPred;
        }
#if JVET_C0024_QTBT
        if( bResetMV )
          pMVField->setAllMv( mv , SIZE_2Nx2N , uiAbsPartIdx , 0 , uiPartIdx );
        pMVField->setAllMvd( mvDiff , SIZE_2Nx2N, uiAbsPartIdx , 0 , uiPartIdx );
#else
        if( bResetMV )
          pMVField->setAllMv( mv , getPartitionSize( uiAbsPartIdx ) , uiAbsPartIdx , 0 , uiPartIdx );
        pMVField->setAllMvd( mvDiff , getPartitionSize( uiAbsPartIdx ) , uiAbsPartIdx , 0 , uiPartIdx );
#endif
      }
    }
  }
#if VCEG_AZ07_FRUC_MERGE
  else if( getFRUCMgrMode( uiAbsPartIdx ) )
  {
    Bool bAvailable = pPred->deriveFRUCMV( this , getDepth( uiAbsPartIdx ) , uiAbsPartIdx , uiPartIdx );
    if( bAvailable == false )
      return( false );
  }
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  else
  {
    UChar    eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
    memset(eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, sizeof(UChar)*MRG_MAX_NUM_CANDS);
#if JVET_C0035_ATMVP_SIMPLIFICATION
    TComMvField*  pcMvFieldSP[NUM_MGR_TYPE] = {NULL};
    UChar* puhInterDirSP[NUM_MGR_TYPE] = {NULL};
    for (Int i=0; i< NUM_MGR_TYPE; i++)
    {
      pcMvFieldSP[i] = new TComMvField[getPic()->getPicSym()->getNumPartitionsInCtu()*2];
      puhInterDirSP[i] = new UChar[getPic()->getPicSym()->getNumPartitionsInCtu()]; 
    }
#else
    TComMvField*  pcMvFieldSP[2] = {NULL,NULL};
    UChar* puhInterDirSP[2] = {NULL,NULL};
    pcMvFieldSP[0] = new TComMvField[getPic()->getPicSym()->getNumPartitionsInCtu()*2]; 
    pcMvFieldSP[1] = new TComMvField[getPic()->getPicSym()->getNumPartitionsInCtu()*2]; 
    puhInterDirSP[0] = new UChar[getPic()->getPicSym()->getNumPartitionsInCtu()]; 
    puhInterDirSP[1] = new UChar[getPic()->getPicSym()->getNumPartitionsInCtu()]; 
#endif
    UInt uiMergeIndex = getMergeIndex( uiAbsPartIdx );
    UInt uiDepth = getDepth( uiAbsPartIdx );
#if JVET_C0024_QTBT
    PartSize ePartSize = SIZE_2Nx2N;
#else
    PartSize ePartSize = getPartitionSize( uiAbsPartIdx );
#endif
    TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
    UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
    Int numValidMergeCand = 0;
#if VCEG_AZ06_IC
    Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif
    getInterMergeCandidates( uiAbsPartIdx, uiPartIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand 
#if VCEG_AZ06_IC
      , abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
      , eMergeCandTypeNieghors
      , pcMvFieldSP
      , puhInterDirSP
#endif
      );

    setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiAbsPartIdx, uiPartIdx, uiDepth );
#if VCEG_AZ06_IC
#if !JVET_C0024_QTBT
    if( getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
#endif
    {
      setICFlagSubParts( getSlice()->getApplyIC() ? abICFlag[uiMergeIndex] : 0, uiAbsPartIdx, uiDepth );
    }
#endif
    TComMv cTmpMv( 0, 0 );
    for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
    {        
      if ( getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
      {
        setMVPIdxSubParts( 0, RefPicList( uiRefListIdx ), uiAbsPartIdx, uiPartIdx, uiDepth);
        setMVPNumSubParts( 0, RefPicList( uiRefListIdx ), uiAbsPartIdx, uiPartIdx, uiDepth);
        getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cTmpMv, ePartSize, uiAbsPartIdx, 0, uiPartIdx );
        getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvField( cMvFieldNeighbours[ 2*uiMergeIndex + uiRefListIdx ], ePartSize, uiAbsPartIdx, 0, uiPartIdx );
      }
    }

    setMergeTypeSubParts( eMergeCandTypeNieghors[uiMergeIndex], uiAbsPartIdx, uiPartIdx, uiDepth ); 


#if JVET_C0035_ATMVP_SIMPLIFICATION
    if( eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_SUBPU_ATMVP || eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_SUBPU_ATMVP_EXT)
#else
    if( eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_SUBPU_TMVP || eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_SUBPU_TMVP_EXT)
#endif
    {
      Int iWidth, iHeight;
      UInt uiIdx;
      getPartIndexAndSize( uiPartIdx, uiIdx, iWidth, iHeight, uiAbsPartIdx, true );

      UInt uiSPAddr;

      Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
      UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeIndex];
#else
      UInt uiSPListIndex = eMergeCandTypeNieghors[uiMergeIndex] == MGR_TYPE_SUBPU_TMVP?0:1;
#endif

      getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

      for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
      {
        getSPAbsPartIdx( uiAbsPartIdx, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
        setInterDirSP(puhInterDirSP[uiSPListIndex][iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
        getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP( this, uiSPAddr, pcMvFieldSP[uiSPListIndex][2*iPartitionIdx], iSPWidth, iSPHeight);
        getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP( this, uiSPAddr, pcMvFieldSP[uiSPListIndex][2*iPartitionIdx + 1], iSPWidth, iSPHeight);
      }
    }
    if ( ( getInterDir( uiAbsPartIdx ) == 3 ) && isBipredRestriction( uiPartIdx ) )
    {
      getCUMvField( REF_PIC_LIST_1 )->setAllMv( TComMv(0,0), ePartSize, uiAbsPartIdx, 0, uiPartIdx);
      getCUMvField( REF_PIC_LIST_1 )->setAllRefIdx( -1, ePartSize, uiAbsPartIdx, 0, uiPartIdx);
      setInterDirSubParts( 1, uiAbsPartIdx, uiPartIdx, uiDepth);
    }
#if JVET_C0035_ATMVP_SIMPLIFICATION
    for (Int i=0; i< NUM_MGR_TYPE; i++)
    {
      delete [] pcMvFieldSP[i];
      delete [] puhInterDirSP[i];
    }
#else
    delete [] pcMvFieldSP[0];
    delete [] pcMvFieldSP[1];
    delete [] puhInterDirSP[0];
    delete [] puhInterDirSP[1];
#endif
  }
#endif
  return( true );
}

Bool TComDataCU::resetMVDandMV2Int( Bool bResetMV 
#if VCEG_AZ07_FRUC_MERGE
  , TComPrediction * pPred
#endif
  )
{
#if !JVET_C0024_QTBT
  UInt uiPartOffset = ( getPic()->getNumPartitionsInCtu() >> ( getDepth(0) << 1 ) ) >> 2;
#endif
  Bool bReset = true;

#if JVET_C0024_QTBT
  UInt uiAbsPartIdx = 0;
  bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
#else
  PartSize ePart = getPartitionSize( 0 );
  UInt uiAbsPartIdx = 0;

  switch ( ePart )
  {
  case SIZE_2Nx2N:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }

  case SIZE_2NxN:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += ( uiPartOffset << 1 );
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }

  case SIZE_Nx2N:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += uiPartOffset;
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }

  case SIZE_NxN:
    {
      for ( Int iPartIdx = 0; iPartIdx < 4; iPartIdx++ )
      {
        bReset &= resetMVDandMV2Int( uiAbsPartIdx , iPartIdx , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
          , pPred
#endif
          );
        uiAbsPartIdx += uiPartOffset;
      }
      break;
    }

  case SIZE_2NxnU:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += ( uiPartOffset >> 1 );
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }

  case SIZE_2NxnD:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += ( uiPartOffset >> 1 ) + ( uiPartOffset << 1 );
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }


  case SIZE_nLx2N:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += ( uiPartOffset >> 2 );
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }

  case SIZE_nRx2N:
    {
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 0 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      uiAbsPartIdx += uiPartOffset + ( uiPartOffset >> 2 );
      bReset &= resetMVDandMV2Int( uiAbsPartIdx , 1 , bResetMV 
#if VCEG_AZ07_FRUC_MERGE
        , pPred
#endif
        );
      break;
    }
  default:
    assert( 0 );
    break;
  }
#endif

  return( bReset );
}
#endif

#if VCEG_AZ07_FRUC_MERGE
UInt TComDataCU::getCtxFRUCMgrMode( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

  // Get BCBP of left PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu+ uiAbsPartIdx );
#endif
  uiCtx    = ( pcTempCU ) ? pcTempCU->getFRUCMgrMode( uiTempPartIdx ) > 0 : 0;

  // Get BCBP of above PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx   += ( pcTempCU ) ? pcTempCU->getFRUCMgrMode( uiTempPartIdx ) > 0 : 0;

  return( uiCtx );
}

UInt TComDataCU::getCtxFRUCME( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

  // Get BCBP of left PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx    = ( pcTempCU ) ? pcTempCU->getFRUCMgrMode( uiTempPartIdx ) == FRUC_MERGE_BILATERALMV : 0;

  // Get BCBP of above PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx   += ( pcTempCU ) ? pcTempCU->getFRUCMgrMode( uiTempPartIdx ) == FRUC_MERGE_BILATERALMV : 0;

  return( uiCtx );
}

#if THIRD_FRUC_MODE
UInt TComDataCU::getCtxFRUCME_UNI(UInt uiAbsPartIdx)
{
  return((UInt)(getWidth(uiAbsPartIdx) <= THIRD_FRUC_MODE_CTX_SWITCH_LEVEL || getHeight(uiAbsPartIdx) <= THIRD_FRUC_MODE_CTX_SWITCH_LEVEL));
}
#endif
Bool TComDataCU::getMvPair( RefPicList eCurRefPicList , const TComMvField & rCurMvField , TComMvField & rMvPair )
{
  Int nTargetRefIdx = getSlice()->getRefIdx4MVPair( eCurRefPicList , rCurMvField.getRefIdx() );
  if (nTargetRefIdx < 0)
  {
      return(false);
  }
  RefPicList eTarRefPicList = ( RefPicList )( 1 - ( Int )eCurRefPicList );
  Int nCurPOC = getSlice()->getPOC();
  Int nRefPOC = getSlice()->getRefPOC( eCurRefPicList , rCurMvField.getRefIdx() );
  Int nTargetPOC = getSlice()->getRefPOC( eTarRefPicList , nTargetRefIdx );
  Int nScale = xGetDistScaleFactor( nCurPOC , nTargetPOC , nCurPOC , nRefPOC );
  rMvPair.getMv() = rCurMvField.getMv().scaleMv( nScale );
  rMvPair.setRefIdx( nTargetRefIdx );

  return( true );
}
#endif

UInt TComDataCU::getCtxInterDir( UInt uiAbsPartIdx )
{
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
  return( Clip3( 0 , 3 , 5 - ((g_aucConvertToBit[getWidth( uiAbsPartIdx )]+g_aucConvertToBit[getHeight(uiAbsPartIdx)]+1)>>1) ) ); 
#else
  return( Clip3( 0 , 3 , 4 - g_aucConvertToBit[getWidth( uiAbsPartIdx )] ) );
#endif
#else
  return getDepth( uiAbsPartIdx );
#endif
}


UChar TComDataCU::getQtRootCbf( UInt uiIdx )
{
  const UInt numberValidComponents = getPic()->getNumberValidComponents();
  return getCbf( uiIdx, COMPONENT_Y, 0 )
          || ((numberValidComponents > COMPONENT_Cb) && getCbf( uiIdx, COMPONENT_Cb, 0 ))
          || ((numberValidComponents > COMPONENT_Cr) && getCbf( uiIdx, COMPONENT_Cr, 0 ));
}

Void TComDataCU::setCbfSubParts( const UInt uiCbf[MAX_NUM_COMPONENT], UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiCbf[COMPONENT_Y], m_puhCbf[COMPONENT_Y], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  setSubPart<UChar>( uiCbf[COMPONENT_Cb], m_puhCbf[COMPONENT_Cb], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  setSubPart<UChar>( uiCbf[COMPONENT_Cr], m_puhCbf[COMPONENT_Cr], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  for(UInt comp=0; comp<MAX_NUM_COMPONENT; comp++)
  {
    memset( m_puhCbf[comp] + uiAbsPartIdx, uiCbf[comp], sizeof( UChar ) * uiCurrPartNumb );
  }
#endif
}

Void TComDataCU::setCbfSubParts( UInt uiCbf, ComponentID compID, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiCbf, m_puhCbf[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  memset( m_puhCbf[compID] + uiAbsPartIdx, uiCbf, sizeof( UChar ) * uiCurrPartNumb );
#endif
}

#if VCEG_AZ07_IMV
#if JVET_E0076_MULTI_PEL_MVD
Void TComDataCU::setiMVFlagSubParts( UChar iMV, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( iMV, m_iMVFlag, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_iMVFlag) == 1 );
  memset( m_iMVFlag + absPartIdx, iMV, m_pcPic->getNumPartitionsInCtu() >> ( 2 * depth ) );
#endif
}
#else
Void TComDataCU::setiMVFlagSubParts( Bool iMV, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<Bool>( iMV, m_iMVFlag, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_iMVFlag) == 1 );
  memset( m_iMVFlag + absPartIdx, iMV, m_pcPic->getNumPartitionsInCtu() >> ( 2 * depth ) );
#endif
}
#endif

Void TComDataCU::setiMVCandNumSubParts( Char ciMVCandNum, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<Char>( ciMVCandNum, m_piMVCandNum, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_piMVCandNum) == 1 );
  memset( m_piMVCandNum + absPartIdx, ciMVCandNum, m_pcPic->getNumPartitionsInCtu() >> ( 2 * depth ) );
#endif
}
#endif

#if COM16_C806_OBMC
Void TComDataCU::setOBMCFlagSubParts( Bool OBMC, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<Bool>( OBMC, m_OBMCFlag, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_OBMCFlag) == 1 );
  memset( m_OBMCFlag + absPartIdx, OBMC, m_pcPic->getNumPartitionsInCtu() >> ( depth << 1 ) );
#endif
}

Bool TComDataCU::isOBMCFlagCoded( UInt uiAbsPartIdx )
{
#if JVET_C0024_QTBT
  if ( isIntra( uiAbsPartIdx ) || ( getMergeFlag( uiAbsPartIdx ) ) || (getWidth( uiAbsPartIdx )*getHeight(uiAbsPartIdx) > COM16_C806_AOBMC_MAXCUSIZE*COM16_C806_AOBMC_MAXCUSIZE) )
#else
  if ( isIntra( uiAbsPartIdx ) || ( getMergeFlag( uiAbsPartIdx ) && getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N ) || getWidth( uiAbsPartIdx ) > COM16_C806_AOBMC_MAXCUSIZE )
#endif
  {
    return false;
  }
  else
  {
    return true;
  }
}
#endif

/** Sets a coded block flag for all sub-partitions of a partition
 * \param uiCbf          The value of the coded block flag to be set
 * \param compID
 * \param uiAbsPartIdx
 * \param uiPartIdx
 * \param uiDepth
 */
#if !JVET_C0024_QTBT
Void TComDataCU::setCbfSubParts ( UInt uiCbf, ComponentID compID, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
  setSubPart<UChar>( uiCbf, m_puhCbf[compID], uiAbsPartIdx, uiDepth, uiPartIdx );
}
#endif

Void TComDataCU::setCbfPartRange ( UInt uiCbf, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiCbf, m_puhCbf[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  memset((m_puhCbf[compID] + uiAbsPartIdx), uiCbf, (sizeof(UChar) * uiCoveredPartIdxes));
#endif
}

#if !JVET_C0024_QTBT
Void TComDataCU::bitwiseOrCbfPartRange( UInt uiCbf, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
  const UInt stopAbsPartIdx = uiAbsPartIdx + uiCoveredPartIdxes;

  for (UInt subPartIdx = uiAbsPartIdx; subPartIdx < stopAbsPartIdx; subPartIdx++)
  {
    m_puhCbf[compID][subPartIdx] |= uiCbf;
  }
}
#endif

Void TComDataCU::setDepthSubParts( UInt uiDepth, UInt uiAbsPartIdx )
{
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
#if JVET_C0024_QTBT
  memset( m_puhDepth[getTextType()] + uiAbsPartIdx, uiDepth, sizeof(UChar)*uiCurrPartNumb );
#else
  memset( m_puhDepth + uiAbsPartIdx, uiDepth, sizeof(UChar)*uiCurrPartNumb );
#endif
}

#if TT_TREE_AFFINE
Void TComDataCU::setCUAbsZIdxInCtuSubParts( UInt uiAbsPartIdx )
{
  setSubPart<UChar>( uiAbsPartIdx&255, m_puhCUabsZIdxInCtu[getTextType()][0], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  setSubPart<UChar>( uiAbsPartIdx>>8,  m_puhCUabsZIdxInCtu[getTextType()][1], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
}
#endif

#if !JVET_C0024_QTBT
Bool TComDataCU::isFirstAbsZorderIdxInDepth (UInt uiAbsPartIdx, UInt uiDepth)
{
  UInt uiPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  return (((m_absZIdxInCtu + uiAbsPartIdx)% uiPartNumb) == 0);
}
#endif

#if JVET_C0024_QTBT
#if TT_TREE_DS1
UInt TComDataCU::getDeltaSizeIdx(UInt uiAbsPartIdx)
{
  UInt uiBTDepth=0;
#if BT_FORCE_BF
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][3][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][2][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#else
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#endif

  UInt uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  UInt uiDeltaSizeIdx = 0;
  while(uiBTSplitModeForBTDepth != 0)
  {
    uiDeltaSizeIdx += (uiBTSplitModeForBTDepth & 1)? 2 : 1; //last bit is 1 means smaller partition in tt
    uiBTDepth++;
    uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  }
  return uiDeltaSizeIdx;
}

UInt TComDataCU::getBTDepthFromBTSplitMode(UInt uiBTSplitMode)
{
  UInt uiBTDepth=0;
  UInt uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  while(uiBTSplitModeForBTDepth != 0)
  {
    uiBTDepth++;
    uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  }
  return uiBTDepth;
}

UInt TComDataCU::getDeltaSizeIdx(UInt uiAbsPartIdx, UInt uiBTWidth, UInt uiBTHeight)
{
  UInt uiCTUSize = getSlice()->getSPS()->getCTUSize();
  UInt uiQTWidth = uiCTUSize>>getDepth(uiAbsPartIdx);
  UInt uiQTHeight = uiCTUSize>>getDepth(uiAbsPartIdx);

  UInt uiDeltaSizeIdx = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiBTWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiBTHeight];
  return uiDeltaSizeIdx;
}

UInt TComDataCU::getBTDepth(UInt uiAbsPartIdx, UInt uiBTWidth, UInt uiBTHeight)
{
  UInt uiCTUSize = getSlice()->getSPS()->getCTUSize();
  UInt uiQTWidth = uiCTUSize>>getDepth(uiAbsPartIdx);
  UInt uiQTHeight = uiCTUSize>>getDepth(uiAbsPartIdx);
  Int iDeltaSizeIdx = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiBTWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiBTHeight];

  UInt uiBTDepth=0;
#if BT_FORCE_BF
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][3][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][2][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#else
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#endif

  UInt uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  while(uiBTSplitModeForBTDepth != 0 && iDeltaSizeIdx>0)
  {
    uiBTDepth++;
    iDeltaSizeIdx -= (uiBTSplitModeForBTDepth&1)? 2:1;
    uiBTSplitModeForBTDepth = (uiBTSplitMode>>(4*uiBTDepth))&15;
  }
  assert(iDeltaSizeIdx==0);
#if !BT_FORCE_BF
  assert(uiBTDepth<=3);
#endif
  return uiBTDepth;
}
#else
UInt TComDataCU::getBTDepth(UInt uiAbsPartIdx, UInt uiBTWidth, UInt uiBTHeight)
{
  UInt uiCTUSize = getSlice()->getSPS()->getCTUSize();
  UInt uiQTWidth = uiCTUSize>>getDepth(uiAbsPartIdx);
  UInt uiQTHeight = uiCTUSize>>getDepth(uiAbsPartIdx);

  UInt uiBTDepth = g_aucConvertToBit[uiQTWidth]-g_aucConvertToBit[uiBTWidth] + g_aucConvertToBit[uiQTHeight]-g_aucConvertToBit[uiBTHeight];
  return uiBTDepth;
}
#endif

UInt TComDataCU::getBTDepth(UInt uiAbsPartIdx)
{
  UInt uiBTDepth=0;
#if BT_FORCE_BF
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][3][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][2][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#else
  UInt uiBTSplitMode = m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTSplitMode <<= 8;
  uiBTSplitMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#endif

#if TT_TREE_DS1
  while(((uiBTSplitMode>>(4*uiBTDepth))&15) != 0)
#else
  while(((uiBTSplitMode>>(2*uiBTDepth))&3) != 0)
#endif
  {
    uiBTDepth++;
  }
  return uiBTDepth;
}

UInt TComDataCU::getBTSplitModeForBTDepth (UInt uiAbsPartIdx, UInt uiBTDepth)
{
#if TT_TREE_DS1
  UInt uiBitPos = 4*uiBTDepth;
#if BT_FORCE_BF
  UInt uiBTMode = m_puhBTSplitMode[getTextType()][3][uiAbsPartIdx];
  uiBTMode <<= 8;
  uiBTMode += m_puhBTSplitMode[getTextType()][2][uiAbsPartIdx];
  uiBTMode <<= 8;
  uiBTMode += m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTMode <<= 8;
  uiBTMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#else
  assert(uiBTDepth <= 3);
  UInt uiBTMode = m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTMode <<= 8;
  uiBTMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];
#endif

  UInt uiBTSplitMode = (uiBTMode>>uiBitPos) & 15;
#else
  UInt uiBitPos = 2*uiBTDepth;
  UInt uiBTMode = m_puhBTSplitMode[getTextType()][1][uiAbsPartIdx];
  uiBTMode <<= 8;
  uiBTMode += m_puhBTSplitMode[getTextType()][0][uiAbsPartIdx];

  UInt uiBTSplitMode = (uiBTMode>>uiBitPos) & 3;
  assert(uiBTSplitMode<3);
#endif

  return uiBTSplitMode;
}

#if TT_TREE_DS1 // data copy
Void TComDataCU::setBTSplitModeSubParts( UInt uiSplitMode, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, UInt uiBTDepth )
{
  assert(uiBTDepth == getBTDepth(uiAbsPartIdx, uiWidth, uiHeight));
#if BT_FORCE_BF
  UChar mask;
  UChar* pMode;
  if (uiBTDepth <= 1)
  {
    mask = 0xff - (15 << ((uiBTDepth) << 2));
    pMode = m_puhBTSplitMode[getTextType()][0];
    assert(uiBTDepth <= 1);
  }
  else if (uiBTDepth <= 3 && uiBTDepth>1)
  {
    mask = 0xff - (15 << ((uiBTDepth - 2) << 2));
    pMode = m_puhBTSplitMode[getTextType()][1];
    uiBTDepth -= 2;
    assert(uiBTDepth <= 1);
  }
  else if (uiBTDepth <= 5 && uiBTDepth > 3)
  {
    mask = 0xff - (15 << ((uiBTDepth - 4) << 2));
    pMode = m_puhBTSplitMode[getTextType()][2];
    uiBTDepth -= 4;
    assert(uiBTDepth <= 1);
  }
  else
  {
    mask = 0xff - (15 << ((uiBTDepth - 6) << 2));
    pMode = m_puhBTSplitMode[getTextType()][3];
    uiBTDepth -= 6;
    assert(uiBTDepth <= 1);
  }
#else
  assert(uiBTDepth <= 3);
  UChar mask = 0xff - (15<<((uiBTDepth<=1?uiBTDepth:(uiBTDepth-2))<<2));
  UChar* pMode = uiBTDepth<=1 ? m_puhBTSplitMode[getTextType()][0] : (uiBTDepth -= 2, m_puhBTSplitMode[getTextType()][1]);
#endif
  pMode[uiAbsPartIdx] &= mask ;
  UInt mode,uiAbsPartIdx2,uiRaster,uiRaster2;
  
  if(uiSplitMode & 2) //if it is tt
  {
    UInt splitDir = uiSplitMode>>2;
    if(splitDir==1)
    {
      UInt uiDepth = g_aucConvertToBit[m_pcSlice->getSPS()->getCTUSize()] - g_aucConvertToBit[uiHeight>>2];
      UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;
      UInt uiStride = m_pcPic->getNumPartInCtuWidth();
      UChar ModeParent = pMode[uiAbsPartIdx];
      uiRaster = g_auiZscanToRaster[uiAbsPartIdx];
      //two small partitions
      mode = ((uiSplitMode|1)<<(uiBTDepth<<2)) | ModeParent;
      uiAbsPartIdx2 = uiAbsPartIdx;
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth, uiHeight>>2 );
      uiRaster2 = uiRaster + uiStride * uiNumPartInShort * 3;
      uiAbsPartIdx2 = g_auiRasterToZscan[uiRaster2];
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth, uiHeight>>2 );
      //one large partition
      mode = ((uiSplitMode&14)<<(uiBTDepth<<2)) | ModeParent;
      uiRaster2 = uiRaster + uiStride * uiNumPartInShort;
      uiAbsPartIdx2 = g_auiRasterToZscan[uiRaster2];
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth, uiHeight>>1 );
    }
    else if(splitDir==2)
    {
      UInt uiDepth = g_aucConvertToBit[m_pcSlice->getSPS()->getCTUSize()] - g_aucConvertToBit[uiWidth>>2];
      UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;
      UInt uiStride = 1;
      UChar ModeParent = pMode[uiAbsPartIdx];
      uiRaster = g_auiZscanToRaster[uiAbsPartIdx];
      //two small partitions
      mode = ((uiSplitMode|1)<<(uiBTDepth<<2)) | ModeParent;
      uiAbsPartIdx2 = uiAbsPartIdx;
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth>>2, uiHeight );
      uiRaster2 = uiRaster + uiStride * uiNumPartInShort * 3;
      uiAbsPartIdx2 = g_auiRasterToZscan[uiRaster2];
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth>>2, uiHeight );
      //one large partition
      mode = ((uiSplitMode&14)<<(uiBTDepth<<2)) | ModeParent;
      uiRaster2 = uiRaster + uiStride * uiNumPartInShort;
      uiAbsPartIdx2 = g_auiRasterToZscan[uiRaster2];
      setSubPart<UChar>( mode, pMode, uiAbsPartIdx2, uiWidth>>1, uiHeight );
    }
    else
      assert(0);
  }
  else
  {
    mode = (uiSplitMode<<(uiBTDepth<<2)) | pMode[uiAbsPartIdx];
    setSubPart<UChar>( mode, pMode, uiAbsPartIdx, uiWidth, uiHeight );
  }
}
#else
Void TComDataCU::setBTSplitModeSubParts( UInt uiSplitMode, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight )
{
  UInt uiBTDepth = getBTDepth(uiAbsPartIdx, uiWidth, uiHeight);

  assert(uiBTDepth<=7);
#if JVET_C0024_QTBT
  UChar mask = 0xff - (3<<((uiBTDepth<=3?uiBTDepth:(uiBTDepth-4))<<1));
#endif
  UChar* pMode = uiBTDepth<=3 ? m_puhBTSplitMode[getTextType()][0] : (uiBTDepth -= 4, m_puhBTSplitMode[getTextType()][1]);
#if JVET_C0024_QTBT
  pMode[uiAbsPartIdx] &= mask ;
#endif
  UInt mode = (uiSplitMode<<(uiBTDepth<<1)) | pMode[uiAbsPartIdx];
  assert(uiBTDepth<=3);

  setSubPart<UChar>( mode, pMode, uiAbsPartIdx, uiWidth, uiHeight );
}
#endif
#else
Void TComDataCU::setPartSizeSubParts( PartSize eMode, UInt uiAbsPartIdx, UInt uiDepth )
{
  assert( sizeof( *m_pePartSize) == 1 );
  memset( m_pePartSize + uiAbsPartIdx, eMode, m_pcPic->getNumPartitionsInCtu() >> ( 2 * uiDepth ) );
}
#endif

#if VCEG_AZ06_IC
Void TComDataCU::setICFlagSubParts( Bool bICFlag, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<Bool>( bICFlag, m_pbICFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  memset( m_pbICFlag + uiAbsPartIdx, bICFlag, (m_pcPic->getNumPartitionsInCtu() >> ( 2 * uiDepth ))*sizeof(Bool) );
#endif
}

Bool TComDataCU::isICFlagCoded( UInt uiAbsPartIdx )
{
  if ( isIntra( uiAbsPartIdx ) || !getSlice()->getApplyIC() )
  {
    return false;
  }

#if JVET_C0024_QTBT
  if( getMergeFlag( uiAbsPartIdx ) 
#else
  if( getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N && getMergeFlag( uiAbsPartIdx ) 
#endif
#if VCEG_AZ07_FRUC_MERGE
    && getFRUCMgrMode( uiAbsPartIdx ) == FRUC_MERGE_OFF 
#endif
    )
  {
    return false;
  }

#if VCEG_AZ06_IC_SPEEDUP
#if JVET_C0024_QTBT
  if ( uiWidth * uiHeight <= 32 )
  {
    return false;
  }
#else
  if( getPartitionSize( uiAbsPartIdx ) != SIZE_2Nx2N )
  {
    return false;
  }
#endif
#endif

#if COM16_C1016_AFFINE
#if JVET_C0024_QTBT
  if ( getAffineFlag( uiAbsPartIdx) )
#else
  if ( getPartitionSize( uiAbsPartIdx) == SIZE_2Nx2N && getAffineFlag( uiAbsPartIdx) )
#endif
  {
    return false;
  }
#endif

#if MRG_CONSIDER_COMPLEX
  if (getMergeComplexFlag(uiAbsPartIdx))
  {
    return false;
  }
#endif

  return true;
}
#endif

Void TComDataCU::setCUTransquantBypassSubParts( Bool flag, UInt uiAbsPartIdx, UInt uiDepth )
{
  memset( m_CUTransquantBypass + uiAbsPartIdx, flag, m_pcPic->getNumPartitionsInCtu() >> ( 2 * uiDepth ) );
}

Void TComDataCU::setSkipFlagSubParts( Bool skip, UInt absPartIdx, UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<Bool>( skip, m_skipFlag, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_skipFlag) == 1 );
  memset( m_skipFlag + absPartIdx, skip, m_pcPic->getNumPartitionsInCtu() >> ( 2 * depth ) );
#endif
}

#if VCEG_AZ05_INTRA_MPI
Void TComDataCU::setMPIIdxSubParts( Char MPIIdx, UInt absPartIdx, UInt depth  )
{
  assert( sizeof( *m_MPIIdx) == 1 );
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset(  m_MPIIdx + absPartIdx, MPIIdx, sizeof(Char)*uiCurrPartNumb );
}
#endif

#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
Void TComDataCU::setPDPCIdxSubParts(Char PDPCIdx, UInt absPartIdx, UInt depth)
{
#if JVET_C0024_QTBT
  setSubPart<Char>( PDPCIdx, m_PDPCIdx, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert(sizeof(*m_PDPCIdx) == 1);
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset(m_PDPCIdx + absPartIdx, PDPCIdx, sizeof(Char)*uiCurrPartNumb);
#endif
}
#endif

#if RRC_DWDIP
Void TComDataCU::setDWDIPIdxSubParts(Char DWDIPIdx, UInt absPartIdx, UInt depth)
{
#if JVET_C0024_QTBT
  setSubPart<Char>(DWDIPIdx, m_DWDIPIdx, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx));
#else
  assert(sizeof(*m_DWDIPIdx) == 1);
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset(m_DWDIPIdx + absPartIdx, DWDIPIdx, sizeof(Char)*uiCurrPartNumb);
#endif
}
#endif

#if RRC_DWDIP_MLR
Void TComDataCU::setDWDIP_GPIdxSubParts(Char DWDIP_GPIdx, UInt absPartIdx, UInt depth)
{
#if JVET_C0024_QTBT
  setSubPart<Char>(DWDIP_GPIdx, m_DWDIP_GPIdx, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx));
#else
  assert(sizeof(*m_DWDIP_GPIdx) == 1);
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset(m_DWDIP_GPIdx + absPartIdx, DWDIP_GPIdx, sizeof(Char)*uiCurrPartNumb);
#endif
}
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
#if JVET_C0024_QTBT
Void TComDataCU::setROTIdxSubParts( const ChannelType channelType, Char ROTIdx, UInt absPartIdx, UInt depth  )
#else
Void TComDataCU::setROTIdxSubParts( Char ROTIdx, UInt absPartIdx, UInt depth  )
#endif
{
#if JVET_C0024_QTBT
  setSubPart<Char>( ROTIdx, m_ROTIdx[channelType], absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof( *m_ROTIdx) == 1 );
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset(  m_ROTIdx + absPartIdx, ROTIdx, sizeof(Char)*uiCurrPartNumb );
#endif
}
#endif
Void TComDataCU::setPredModeSubParts( PredMode eMode, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  assert( sizeof( *m_pePredMode) == 1 );
  setSubPart<Char>( eMode, m_pePredMode, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  assert( sizeof( *m_pePredMode) == 1 );
  memset( m_pePredMode + uiAbsPartIdx, eMode, m_pcPic->getNumPartitionsInCtu() >> ( 2 * uiDepth ) );
#endif
}

Void TComDataCU::setChromaQpAdjSubParts( UChar val, Int absPartIdx, Int depth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( val, m_ChromaQpAdj, absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  assert( sizeof(*m_ChromaQpAdj) == 1 );
  memset( m_ChromaQpAdj + absPartIdx, val, m_pcPic->getNumPartitionsInCtu() >> ( 2 * depth ) );
#endif
}

#if JVET_C0024_DELTA_QP_FIX
Void TComDataCU::setQPSubCUs( Int qp, UInt absPartIdx, UInt depth, UInt uiWidth, UInt uiHeight, UInt& ruiFirstNonZeroPartIdx, Bool &foundNonZeroCbf )
#else
Void TComDataCU::setQPSubCUs( Int qp, UInt absPartIdx, UInt depth, Bool &foundNonZeroCbf )
#endif
{
  UInt currPartNumb = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  UInt currPartNumQ = currPartNumb >> 2;
  const UInt numValidComp = m_pcPic->getNumberValidComponents();

  if(!foundNonZeroCbf)
  {
    if(getDepth(absPartIdx) > depth)
    {
      for ( UInt partUnitIdx = 0; partUnitIdx < 4; partUnitIdx++ )
      {
#if JVET_C0024_DELTA_QP_FIX
        setQPSubCUs( qp, absPartIdx+partUnitIdx*currPartNumQ, depth+1, uiWidth>>1, uiHeight>>1, ruiFirstNonZeroPartIdx, foundNonZeroCbf );
#else
        setQPSubCUs( qp, absPartIdx+partUnitIdx*currPartNumQ, depth+1, foundNonZeroCbf );
#endif
      }
    }
    else
    {
#if JVET_C0024_DELTA_QP_FIX
      // check for BT split
      UInt uiBTDepth = getBTDepth(absPartIdx, uiWidth, uiHeight);

      UInt uiBTSplitMode = getBTSplitModeForBTDepth(absPartIdx, uiBTDepth);
#if TT_TREE_DS1
      if (uiBTSplitMode==ST_HORI_BT)
#else
      if (uiBTSplitMode==1)
#endif
      {
        for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
        {
          if (uiPartUnitIdx==1)
          {
            absPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] 
            + (uiHeight>>1)/getPic()->getMinCUHeight()*getPic()->getNumPartInCtuWidth()];
          }
          setQPSubCUs( qp, absPartIdx, depth, uiWidth, uiHeight>>1, ruiFirstNonZeroPartIdx, foundNonZeroCbf );
        }
        return;
      }
#if TT_TREE_DS1
      else if (uiBTSplitMode==ST_VERT_BT)
#else
      else if (uiBTSplitMode==2)
#endif
      {
        for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 2; uiPartUnitIdx++ )
        {
          if (uiPartUnitIdx==1)
          {
            absPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] 
            + (uiWidth>>1)/getPic()->getMinCUWidth()];
          }
          setQPSubCUs( qp, absPartIdx, depth, uiWidth>>1, uiHeight, ruiFirstNonZeroPartIdx, foundNonZeroCbf );
        }
        return;
      }
#if TT_ENC_OPT //set qp encoder
      else if(uiBTSplitMode==ST_HORI_TT_L || uiBTSplitMode==ST_HORI_TT_S)
      {
        UInt uiHeightSubPart, absPartIdxSubPart;
        for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
        {
          if (uiPartUnitIdx==1)
          {
            absPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + (uiHeight>>2)/getPic()->getMinCUHeight()*getPic()->getNumPartInCtuWidth()];
            uiHeightSubPart = uiHeight>>1;
          }
          else
          {
            if(uiPartUnitIdx==2)
              absPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + ((uiHeight>>2)+(uiHeight>>1))/getPic()->getMinCUHeight()*getPic()->getNumPartInCtuWidth()];
            else
              absPartIdxSubPart = absPartIdx;
            uiHeightSubPart = uiHeight>>2;
          }
          setQPSubCUs( qp, absPartIdx, depth, uiWidth, uiHeightSubPart, ruiFirstNonZeroPartIdx, foundNonZeroCbf );
        }
        return;
      }
      else if(uiBTSplitMode==ST_VERT_TT_L || uiBTSplitMode==ST_VERT_TT_S)
      {
        UInt uiWidthSubPart, absPartIdxSubPart;
        for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 3; uiPartUnitIdx++ )
        {
          if (uiPartUnitIdx==1)
          {
            absPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + (uiWidth>>2)/getPic()->getMinCUWidth()];
            uiWidthSubPart = uiHeight>>1;
          }
          else
          {
            if(uiPartUnitIdx==2)
              absPartIdxSubPart = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + ((uiWidth>>2)+(uiWidth>>1))/getPic()->getMinCUWidth()];
            else
              absPartIdxSubPart = absPartIdx;
            uiWidthSubPart = uiHeight>>2;
          }
          setQPSubCUs( qp, absPartIdx, depth, uiWidthSubPart, uiHeight, ruiFirstNonZeroPartIdx, foundNonZeroCbf );
        }
        return;
      }
#endif
      else
      {
        if( getSlice()->isIntra() )
        {
          if( getTextType() == CHANNEL_TYPE_LUMA )
          {
            if( getCbf( absPartIdx, COMPONENT_Y ) )
            {
              foundNonZeroCbf = true;
            }
          }
          else
          {
            if( ( numValidComp>COMPONENT_Cb && getCbf( absPartIdx, COMPONENT_Cb ) ) || ( numValidComp>COMPONENT_Cr && getCbf( absPartIdx, COMPONENT_Cr) ) )
            {
              foundNonZeroCbf = true;
            }
          }
        }
        else if(getCbf( absPartIdx, COMPONENT_Y ) || (numValidComp>COMPONENT_Cb && getCbf( absPartIdx, COMPONENT_Cb )) || (numValidComp>COMPONENT_Cr && getCbf( absPartIdx, COMPONENT_Cr) ) )
        {
          foundNonZeroCbf = true;
        }

        if( !foundNonZeroCbf )
        {
          setQPSubParts(qp, absPartIdx, uiWidth, uiHeight);
        }
        else
        {
          ruiFirstNonZeroPartIdx = absPartIdx;
        }
      }
#else
      if(getCbf( absPartIdx, COMPONENT_Y ) || (numValidComp>COMPONENT_Cb && getCbf( absPartIdx, COMPONENT_Cb )) || (numValidComp>COMPONENT_Cr && getCbf( absPartIdx, COMPONENT_Cr) ) )
      {
        foundNonZeroCbf = true;
      }
      else
      {
        setQPSubParts(qp, absPartIdx, depth);
      }
#endif
    }
  }
}

#if JVET_C0024_DELTA_QP_FIX
Void TComDataCU::setQPSubParts( Int qp, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight )
#else
Void TComDataCU::setQPSubParts( Int qp, UInt uiAbsPartIdx, UInt uiDepth )
#endif
{
#if JVET_C0024_QTBT
#if JVET_C0024_DELTA_QP_FIX
  setSubPart<Char>( qp, m_phQP[getTextType()], uiAbsPartIdx, uiWidth, uiHeight );
  if( !getSlice()->isIntra() )
  {
    assert( getTextType() == CHANNEL_TYPE_LUMA );
  }
#else
  setSubPart<Char>( qp, m_phQP, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#endif
#else
  const UInt numPart = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  memset(m_phQP+uiAbsPartIdx, qp, numPart);
#endif
}

Void TComDataCU::setIntraDirSubParts( const ChannelType channelType, const UInt dir, const UInt absPartIdx, const UInt depth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( dir, m_puhIntraDir[channelType], absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx) );
#else
  UInt numPart = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
  memset( m_puhIntraDir[channelType] + absPartIdx, dir,sizeof(UChar)*numPart );
#endif
}

#if MLIP_EXTEND_REF_LINE
Void TComDataCU::setExternRefSubParts(const ChannelType channelType, const UInt dir, const UInt absPartIdx, const UInt depth)
{
#if JVET_C0024_QTBT
    setSubPart<UChar>(dir, m_puhExternRef[channelType], absPartIdx, getWidth(absPartIdx), getHeight(absPartIdx));
#else
    UInt numPart = m_pcPic->getNumPartitionsInCtu() >> (depth << 1);
    memset(m_puhExternRef[channelType] + absPartIdx, dir, sizeof(UChar)*numPart);
#endif
}
#endif

template<typename T>
#if JVET_C0024_QTBT
Void TComDataCU::setSubPart( T uiParameter, T* puhBaseLCU, UInt uiCUAddr, UInt uiWidth, UInt uiHeight )
{
  assert( sizeof(T) == 1 ); // Using memset() works only for types of size 1
  assert(((uiWidth & (uiWidth-1))==0 ) && ((uiHeight & (uiHeight-1))==0 ) );

  UInt uiRaster = g_auiZscanToRaster[m_absZIdxInCtu+uiCUAddr];
  UInt uiShort, uiLong;
  UInt uiStride;

#if TT_TREE_DS1 //data copy
  //always split into four parts if possible
  UInt uiNumXDiv = (uiWidth  > m_pcPic->getMinCUWidth() )? 2 : 1;
  UInt uiNumYDiv = (uiHeight > m_pcPic->getMinCUHeight())? 2 : 1;
  uiWidth  = (uiNumXDiv==1)? uiWidth  : (uiWidth>>1);
  uiHeight = (uiNumYDiv==1)? uiHeight : (uiHeight>>1);
  UInt uiLenX    = uiWidth  / m_pcPic->getMinCUWidth();
  UInt uiLenY    = uiHeight / m_pcPic->getMinCUHeight();

  for(UInt yIdx=0; yIdx < uiNumYDiv; yIdx++)
  {
    for(UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
    {
      if(m_bUseAbsoluteZIdxToCtu==false)
        uiRaster = g_auiZscanToRaster[uiCUAddr] + (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
      else
        uiRaster = g_auiZscanToRaster[m_absZIdxInCtu+uiCUAddr] + (uiLenX * xIdx) + (uiLenY * yIdx) * m_pcPic->getNumPartInCtuWidth();
#endif
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = m_pcPic->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  UInt uiDepth = g_aucConvertToBit[m_pcSlice->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;

  for (UInt i=0; i<uiLong; i+=uiShort)
  {
#if TT_TREE_DS1 //call at both encoder & decoder
    if(m_bUseAbsoluteZIdxToCtu==false)
      memset(puhBaseLCU + g_auiRasterToZscan[uiRaster], uiParameter, uiCurrPartNumb );
    else
      memset(puhBaseLCU + g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu, uiParameter, uiCurrPartNumb );
#else
    memset(puhBaseLCU + g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu, uiParameter, uiCurrPartNumb );
#endif
    uiRaster += uiNumPartInShort * uiStride;
  }
#if TT_TREE_DS1 //data copy end
  }}
#endif
}
#else
Void TComDataCU::setSubPart( T uiParameter, T* puhBaseCtu, UInt uiCUAddr, UInt uiCUDepth, UInt uiPUIdx )
{
  assert( sizeof(T) == 1 ); // Using memset() works only for types of size 1
#if COM16_C806_OBMC
  UInt uiCurrPartNum  = (m_pcPic->getNumPartitionsInCtu() >> (2 * uiCUDepth));
  UInt uiCurrPartNumQ = uiCurrPartNum >> 2;
#else
  UInt uiCurrPartNumQ = (m_pcPic->getNumPartitionsInCtu() >> (2 * uiCUDepth)) >> 2;
#endif
  switch ( m_pePartSize[ uiCUAddr ] )
  {
    case SIZE_2Nx2N:
#if COM16_C806_OBMC
      memset( puhBaseCtu + uiCUAddr, uiParameter, uiCurrPartNum );
#else
      memset( puhBaseCtu + uiCUAddr, uiParameter, 4 * uiCurrPartNumQ );
#endif
      break;
    case SIZE_2NxN:
      memset( puhBaseCtu + uiCUAddr, uiParameter, 2 * uiCurrPartNumQ );
      break;
    case SIZE_Nx2N:
      memset( puhBaseCtu + uiCUAddr, uiParameter, uiCurrPartNumQ );
      memset( puhBaseCtu + uiCUAddr + 2 * uiCurrPartNumQ, uiParameter, uiCurrPartNumQ );
      break;
    case SIZE_NxN:
      memset( puhBaseCtu + uiCUAddr, uiParameter, uiCurrPartNumQ );
      break;
    case SIZE_2NxnU:
      if ( uiPUIdx == 0 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 1) );
        memset( puhBaseCtu + uiCUAddr + uiCurrPartNumQ, uiParameter, (uiCurrPartNumQ >> 1) );
      }
      else if ( uiPUIdx == 1 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 1) );
        memset( puhBaseCtu + uiCUAddr + uiCurrPartNumQ, uiParameter, ((uiCurrPartNumQ >> 1) + (uiCurrPartNumQ << 1)) );
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_2NxnD:
      if ( uiPUIdx == 0 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, ((uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1)) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1) + uiCurrPartNumQ, uiParameter, (uiCurrPartNumQ >> 1) );
      }
      else if ( uiPUIdx == 1 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 1) );
        memset( puhBaseCtu + uiCUAddr + uiCurrPartNumQ, uiParameter, (uiCurrPartNumQ >> 1) );
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_nLx2N:
      if ( uiPUIdx == 0 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
      }
      else if ( uiPUIdx == 1 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)) );
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_nRx2N:
      if ( uiPUIdx == 0 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)) );
        memset( puhBaseCtu + uiCUAddr + uiCurrPartNumQ + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1), uiParameter, (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1) + uiCurrPartNumQ + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
      }
      else if ( uiPUIdx == 1 )
      {
        memset( puhBaseCtu + uiCUAddr, uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1), uiParameter, (uiCurrPartNumQ >> 2) );
        memset( puhBaseCtu + uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1), uiParameter, (uiCurrPartNumQ >> 2) );
      }
      else
      {
        assert(0);
      }
      break;
    default:
      assert( 0 );
      break;
  }
}
#endif //JVET_C0024_QTBT

Void TComDataCU::setMergeFlagSubParts ( Bool bMergeFlag, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<Bool>( bMergeFlag, m_pbMergeFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart( bMergeFlag, m_pbMergeFlag, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

Void TComDataCU::setMergeIndexSubParts ( UInt uiMergeIndex, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiMergeIndex, m_puhMergeIndex, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart<UChar>( uiMergeIndex, m_puhMergeIndex, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

#if VCEG_AZ07_FRUC_MERGE
Void TComDataCU::setFRUCMgrModeSubParts ( UChar uhFRUCMgrMode, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uhFRUCMgrMode, m_puhFRUCMgrMode, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart( uhFRUCMgrMode, m_puhFRUCMgrMode, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}
#endif

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
Void TComDataCU::setMergeTypeSubParts( UChar eMergeType, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( eMergeType, m_peMergeType, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart<UChar>( eMergeType, m_peMergeType, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}
#endif

Void TComDataCU::setInterDirSubParts( UInt uiDir, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiDir, m_puhInterDir, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart<UChar>( uiDir, m_puhInterDir, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

Void TComDataCU::setMVPIdxSubParts( Int iMVPIdx, RefPicList eRefPicList, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<Char>( iMVPIdx, m_apiMVPIdx[eRefPicList], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart<Char>( iMVPIdx, m_apiMVPIdx[eRefPicList], uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

Void TComDataCU::setMVPNumSubParts( Int iMVPNum, RefPicList eRefPicList, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<Char>( iMVPNum, m_apiMVPNum[eRefPicList], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart<Char>( iMVPNum, m_apiMVPNum[eRefPicList], uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

#if DIMD_INTRA_PRED
Void TComDataCU::setDIMDEnabledFlagSubParts(ChannelType channelType, UChar uc, UInt uiAbsPartIdx, UInt uiDepth)
{
  setSubPart<UChar>(uc, m_ucDIMDEnabledFlag[channelType], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}
#if  TT_TREE_DIMD
Void TComDataCU::setDIMDEnabledFlagSubParts_TT(ChannelType channelType, UChar uc, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
  setSubPart<UChar>(uc, m_ucDIMDEnabledFlag[channelType], uiAbsPartIdx, uiWidth, uiHeight);
}
#endif
#if DIMD_NUM_INTRA_DIR_INC
Void TComDataCU::setExtIntraDirSubParts(ChannelType channelType, UInt uiDir, UInt uiAbsPartIdx, UInt uiDepth)
{
  setSubPart<UChar>(uiDir, m_ucExtIntraDir[channelType], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}
#endif
Void TComDataCU::setDIMDNoBTLevelFlagSubParts(ChannelType channelType, UChar uc, UInt uiAbsPartIdx, UInt uiDepth)
{
  setSubPart<UChar>(uc, m_ucDIMDNoBTLevelFlag[channelType], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}
#if TT_TREE_DIMD
Void TComDataCU::setDIMDNoBTLevelFlagSubParts_TT(ChannelType channelType, UChar uc, UInt uiAbsPartIdx, UInt uiDepth, UInt uiWidth, UInt uiHeight)
{
  setSubPart<UChar>(uc, m_ucDIMDNoBTLevelFlag[channelType], uiAbsPartIdx, uiWidth, uiHeight);
}
#endif

Void TComDataCU::setSecIntraDirSubParts(const UInt uiDir, const UInt uiAbsPartIdx, const UInt uiDepth)
{
  setSubPart<UChar>(uiDir, m_puhSecIntraDir, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}
#endif

#if COM16_C806_EMT
Void TComDataCU::setEmtTuIdxSubParts( UInt uiTrMode, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiTrMode, m_puhEmtTuIdx, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  
  memset( m_puhEmtTuIdx + uiAbsPartIdx, uiTrMode, sizeof(UChar)*uiCurrPartNumb );
#endif
}

Void TComDataCU::setEmtTuIdxPartsRange( UInt uiTrMode, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
  if( COMPONENT_Y==compID )
  {
#if JVET_C0024_QTBT
    setSubPart<UChar>( uiTrMode, m_puhEmtTuIdx, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
    memset( m_puhEmtTuIdx + uiAbsPartIdx, uiTrMode, (sizeof(Char) * uiCoveredPartIdxes));
#endif
  }
}

Void TComDataCU::setEmtCuFlagSubParts( UInt uiTrMode, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( uiTrMode, m_puhEmtCuFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  
  memset( m_puhEmtCuFlag + uiAbsPartIdx, uiTrMode, sizeof(UChar)*uiCurrPartNumb );
#endif
}
#endif

#if !JVET_C0024_QTBT
Void TComDataCU::setTrIdxSubParts( UInt uiTrIdx, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

  memset( m_puhTrIdx + uiAbsPartIdx, uiTrIdx, sizeof(UChar)*uiCurrPartNumb );
}
#endif

Void TComDataCU::setTransformSkipSubParts( const UInt useTransformSkip[MAX_NUM_COMPONENT], UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( useTransformSkip[COMPONENT_Y], m_puhTransformSkip[COMPONENT_Y], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  setSubPart<UChar>( useTransformSkip[COMPONENT_Cb], m_puhTransformSkip[COMPONENT_Cb], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  setSubPart<UChar>( useTransformSkip[COMPONENT_Cr], m_puhTransformSkip[COMPONENT_Cr], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

  for(UInt i=0; i<MAX_NUM_COMPONENT; i++)
  {
    memset( m_puhTransformSkip[i] + uiAbsPartIdx, useTransformSkip[i], sizeof( UChar ) * uiCurrPartNumb );
  }
#endif
}

Void TComDataCU::setTransformSkipSubParts( UInt useTransformSkip, ComponentID compID, UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( useTransformSkip, m_puhTransformSkip[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

  memset( m_puhTransformSkip[compID] + uiAbsPartIdx, useTransformSkip, sizeof( UChar ) * uiCurrPartNumb );
#endif
}

Void TComDataCU::setTransformSkipPartRange ( UInt useTransformSkip, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( useTransformSkip, m_puhTransformSkip[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  memset((m_puhTransformSkip[compID] + uiAbsPartIdx), useTransformSkip, (sizeof(UChar) * uiCoveredPartIdxes));
#endif
}

Void TComDataCU::setCrossComponentPredictionAlphaPartRange( Char alphaValue, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
#if JVET_C0024_QTBT
  setSubPart<Char>( alphaValue, m_crossComponentPredictionAlpha[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  memset((m_crossComponentPredictionAlpha[compID] + uiAbsPartIdx), alphaValue, (sizeof(Char) * uiCoveredPartIdxes));
#endif
}

Void TComDataCU::setExplicitRdpcmModePartRange ( UInt rdpcmMode, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes )
{
#if JVET_C0024_QTBT
  setSubPart<UChar>( rdpcmMode, m_explicitRdpcmMode[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  memset((m_explicitRdpcmMode[compID] + uiAbsPartIdx), rdpcmMode, (sizeof(UChar) * uiCoveredPartIdxes));
#endif
}

#if VCEG_AZ08_KLT_COMMON
Void TComDataCU::setKLTFlagSubParts(const UInt useKLT[MAX_NUM_COMPONENT], UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  for (UInt i = 0; i<MAX_NUM_COMPONENT; i++)
  {
    setSubPart<UChar>(useKLT[i], m_puhKLTFlag[i], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
  }
#else
    UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

    for (UInt i = 0; i<MAX_NUM_COMPONENT; i++)
    {
        memset(m_puhKLTFlag[i] + uiAbsPartIdx, useKLT[i], sizeof(UChar)* uiCurrPartNumb);
    }
#endif
}
Void TComDataCU::setKLTFlagSubParts(UInt useKLTY, ComponentID compID, UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  setSubPart<UChar>(useKLTY, m_puhKLTFlag[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
    UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

    memset(m_puhKLTFlag[compID] + uiAbsPartIdx, useKLTY, sizeof(UChar)* uiCurrPartNumb);
#endif
}
Void TComDataCU::setKLTPartRange(UInt useKLTY, ComponentID compID, UInt uiAbsPartIdx, UInt uiCoveredPartIdxes)
{
#if JVET_C0024_QTBT
  setSubPart<UChar>(useKLTY, m_puhKLTFlag[compID], uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
    memset((m_puhKLTFlag[compID] + uiAbsPartIdx), useKLTY, (sizeof(UChar)* uiCoveredPartIdxes));
#endif
}
#endif

Void TComDataCU::setSizeSubParts( UInt uiWidth, UInt uiHeight, UInt uiAbsPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  assert(uiWidth<=256 && uiHeight<=256); //within one byte;
  ChannelType eCType  = getTextType();
#if JVET_C0024_CTU_256
  setSubPart<UChar>(uiWidth >> MIN_CU_LOG2, m_puhWidth[eCType], uiAbsPartIdx, uiWidth , uiHeight);
  setSubPart<UChar>(uiHeight>> MIN_CU_LOG2, m_puhHeight[eCType], uiAbsPartIdx, uiWidth, uiHeight);
#else
  setSubPart<UChar>(uiWidth, m_puhWidth[eCType], uiAbsPartIdx, uiWidth, uiHeight);
  setSubPart<UChar>(uiHeight, m_puhHeight[eCType], uiAbsPartIdx, uiWidth, uiHeight);
#endif
#else
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

#if COM16_C806_LARGE_CTU
  for( UInt n = 0 ; n < uiCurrPartNumb ; n++ )
  {
    m_puhWidth[uiAbsPartIdx+n] = uiWidth;
    m_puhHeight[uiAbsPartIdx+n] = uiHeight;
  }
#else
  memset( m_puhWidth  + uiAbsPartIdx, uiWidth,  sizeof(UChar)*uiCurrPartNumb );
  memset( m_puhHeight + uiAbsPartIdx, uiHeight, sizeof(UChar)*uiCurrPartNumb );
#endif
#endif
}

UChar TComDataCU::getNumPartitions(const UInt uiAbsPartIdx)
{
#if JVET_C0024_QTBT
  return 1;
#else
  UChar iNumPart = 0;

  switch ( m_pePartSize[uiAbsPartIdx] )
  {
    case SIZE_2Nx2N:    iNumPart = 1; break;
    case SIZE_2NxN:     iNumPart = 2; break;
    case SIZE_Nx2N:     iNumPart = 2; break;
    case SIZE_NxN:      iNumPart = 4; break;
    case SIZE_2NxnU:    iNumPart = 2; break;
    case SIZE_2NxnD:    iNumPart = 2; break;
    case SIZE_nLx2N:    iNumPart = 2; break;
    case SIZE_nRx2N:    iNumPart = 2; break;
    default:            assert (0);   break;
  }

  return  iNumPart;
#endif
}

// This is for use by a leaf/sub CU object only, with no additional AbsPartIdx
Void TComDataCU::getPartIndexAndSize( UInt uiPartIdx, UInt& ruiPartAddr, Int& riWidth, Int& riHeight )
{
#if JVET_C0024_QTBT
  riWidth = getWidth(0);      riHeight = getHeight(0);      ruiPartAddr = 0;
#else
  switch ( m_pePartSize[0] )
  {
    case SIZE_2NxN:
      riWidth = getWidth(0);      riHeight = getHeight(0) >> 1; ruiPartAddr = ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 1;
      break;
    case SIZE_Nx2N:
      riWidth = getWidth(0) >> 1; riHeight = getHeight(0);      ruiPartAddr = ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 2;
      break;
    case SIZE_NxN:
      riWidth = getWidth(0) >> 1; riHeight = getHeight(0) >> 1; ruiPartAddr = ( m_uiNumPartition >> 2 ) * uiPartIdx;
      break;
    case SIZE_2NxnU:
      riWidth     = getWidth(0);
      riHeight    = ( uiPartIdx == 0 ) ?  getHeight(0) >> 2 : ( getHeight(0) >> 2 ) + ( getHeight(0) >> 1 );
      ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : m_uiNumPartition >> 3;
      break;
    case SIZE_2NxnD:
      riWidth     = getWidth(0);
      riHeight    = ( uiPartIdx == 0 ) ?  ( getHeight(0) >> 2 ) + ( getHeight(0) >> 1 ) : getHeight(0) >> 2;
      ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : (m_uiNumPartition >> 1) + (m_uiNumPartition >> 3);
      break;
    case SIZE_nLx2N:
      riWidth     = ( uiPartIdx == 0 ) ? getWidth(0) >> 2 : ( getWidth(0) >> 2 ) + ( getWidth(0) >> 1 );
      riHeight    = getHeight(0);
      ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : m_uiNumPartition >> 4;
      break;
    case SIZE_nRx2N:
      riWidth     = ( uiPartIdx == 0 ) ? ( getWidth(0) >> 2 ) + ( getWidth(0) >> 1 ) : getWidth(0) >> 2;
      riHeight    = getHeight(0);
      ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : (m_uiNumPartition >> 2) + (m_uiNumPartition >> 4);
      break;
    default:
      assert ( m_pePartSize[0] == SIZE_2Nx2N );
      riWidth = getWidth(0);      riHeight = getHeight(0);      ruiPartAddr = 0;
      break;
  }
#endif
}


Void TComDataCU::getMvField ( TComDataCU* pcCU, UInt uiAbsPartIdx, RefPicList eRefPicList, TComMvField& rcMvField )
{
  if ( pcCU == NULL )  // OUT OF BOUNDARY
  {
    TComMv  cZeroMv;
    rcMvField.setMvField( cZeroMv, NOT_VALID );
    return;
  }

  TComCUMvField*  pcCUMvField = pcCU->getCUMvField( eRefPicList );
  rcMvField.setMvField( pcCUMvField->getMv( uiAbsPartIdx ), pcCUMvField->getRefIdx( uiAbsPartIdx ) );
}

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
Void TComDataCU::getPartIndexAndSize( UInt uiPartIdx, UInt& ruiPartAddr, Int& riWidth, Int& riHeight, UInt uiAbsPartIdx, Bool bLCU)
{
#if JVET_C0024_QTBT
  UInt  uiTmpAbsPartIdx  = bLCU ? uiAbsPartIdx : 0;
  riWidth = getWidth( uiTmpAbsPartIdx );      riHeight = getHeight( uiTmpAbsPartIdx );      ruiPartAddr = 0;
#else
  UInt uiNumPartition  = bLCU ? (getWidth(uiAbsPartIdx)*getHeight(uiAbsPartIdx) >> 4) : m_uiNumPartition;
  UInt  uiTmpAbsPartIdx  = bLCU ? uiAbsPartIdx : 0;

  switch ( m_pePartSize[uiTmpAbsPartIdx] )
  {
  case SIZE_2NxN:
    riWidth = getWidth( uiTmpAbsPartIdx );      riHeight = getHeight( uiTmpAbsPartIdx ) >> 1; ruiPartAddr = ( uiPartIdx == 0 )? 0 : uiNumPartition >> 1;
    break;
  case SIZE_Nx2N:
    riWidth = getWidth( uiTmpAbsPartIdx ) >> 1; riHeight = getHeight( uiTmpAbsPartIdx );      ruiPartAddr = ( uiPartIdx == 0 )? 0 : uiNumPartition >> 2;
    break;
  case SIZE_NxN:
    riWidth = getWidth( uiTmpAbsPartIdx ) >> 1; riHeight = getHeight( uiTmpAbsPartIdx ) >> 1; ruiPartAddr = ( uiNumPartition >> 2 ) * uiPartIdx;
    break;
  case SIZE_2NxnU:
    riWidth     = getWidth( uiTmpAbsPartIdx );
    riHeight    = ( uiPartIdx == 0 ) ?  getHeight( uiTmpAbsPartIdx ) >> 2 : ( getHeight( uiTmpAbsPartIdx ) >> 2 ) + ( getHeight( uiTmpAbsPartIdx ) >> 1 );
    ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : uiNumPartition >> 3;
    break;
  case SIZE_2NxnD:
    riWidth     = getWidth( uiTmpAbsPartIdx );
    riHeight    = ( uiPartIdx == 0 ) ?  ( getHeight( uiTmpAbsPartIdx ) >> 2 ) + ( getHeight( uiTmpAbsPartIdx ) >> 1 ) : getHeight( uiTmpAbsPartIdx ) >> 2;
    ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : (uiNumPartition >> 1) + (uiNumPartition >> 3);
    break;
  case SIZE_nLx2N:
    riWidth     = ( uiPartIdx == 0 ) ? getWidth( uiTmpAbsPartIdx ) >> 2 : ( getWidth( uiTmpAbsPartIdx ) >> 2 ) + ( getWidth( uiTmpAbsPartIdx ) >> 1 );
    riHeight    = getHeight( uiTmpAbsPartIdx );
    ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : uiNumPartition >> 4;
    break;
  case SIZE_nRx2N:
    riWidth     = ( uiPartIdx == 0 ) ? ( getWidth( uiTmpAbsPartIdx ) >> 2 ) + ( getWidth( uiTmpAbsPartIdx ) >> 1 ) : getWidth( uiTmpAbsPartIdx ) >> 2;
    riHeight    = getHeight( uiTmpAbsPartIdx );
    ruiPartAddr = ( uiPartIdx == 0 ) ? 0 : (uiNumPartition >> 2) + (uiNumPartition >> 4);
    break;
  default:
    assert ( m_pePartSize[uiTmpAbsPartIdx] == SIZE_2Nx2N ); 
    riWidth = getWidth( uiTmpAbsPartIdx );      riHeight = getHeight( uiTmpAbsPartIdx );      ruiPartAddr = 0;
    break;
  }
#endif
}
#if !JVET_C0035_ATMVP_SIMPLIFICATION
Void TComDataCU::get1stTvFromSpatialNeighbor ( UInt uiAbsPartIdx, UInt uiPUIdx, Bool &bTvAva, Int &iPOC, TComMv &rcMv)
{
  TComPic *pColPic = getSlice()->getRefPic( RefPicList(getSlice()->isInterB() ? 1-getSlice()->getColFromL0Flag() : 0), getSlice()->getColRefIdx());

  iPOC             = pColPic->getPOC();

  UInt uiAbsPartAddr = m_absZIdxInCtu + uiAbsPartIdx;

  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
  deriveLeftRightTopIdxGeneral( uiAbsPartIdx, uiPUIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdxGeneral  ( uiAbsPartIdx, uiPUIdx, uiPartIdxLB );

  TComMvField tempMvField;
  bTvAva            = false;

  //left
  UInt uiLeftPartIdx = 0;
  TComDataCU* pcCULeft = 0;
  pcCULeft = getPULeft( uiLeftPartIdx, uiPartIdxLB );

  // above
  UInt uiAbovePartIdx = 0;
  TComDataCU* pcCUAbove = 0;
  pcCUAbove = getPUAbove( uiAbovePartIdx, uiPartIdxRT );

   // above right
  UInt uiAboveRightPartIdx = 0;
  TComDataCU* pcCUAboveRight = 0;
  pcCUAboveRight = getPUAboveRight( uiAboveRightPartIdx, uiPartIdxRT );

      //left bottom
  UInt uiLeftBottomPartIdx = 0;
  TComDataCU* pcCULeftBottom = 0;
  pcCULeftBottom = this->getPUBelowLeft( uiLeftBottomPartIdx, uiPartIdxLB );

  // above left
  UInt uiAboveLeftPartIdx = 0;
  TComDataCU* pcCUAboveLeft = 0;
  pcCUAboveLeft = getPUAboveLeft( uiAboveLeftPartIdx, uiAbsPartAddr );

  TComDataCU* pcNeighorPU[5]   = {pcCULeft,      pcCUAbove,        pcCUAboveRight,       pcCULeftBottom,       pcCUAboveLeft,};
  UInt        uiNeighorPartIdx[5] = {uiLeftPartIdx, uiAbovePartIdx,   uiAboveRightPartIdx,  uiLeftBottomPartIdx,  uiAboveLeftPartIdx,};
  for (UInt uiN=0; uiN< 5; uiN++) 
    for( UInt uiCurrRefListId = 0; uiCurrRefListId < (getSlice()->getSliceType() == B_SLICE ?  2 : 1 ) ; uiCurrRefListId++ )
    {
      RefPicList  eCurrRefPicList = RefPicList( RefPicList( getSlice()->isInterB() ? (getSlice()->getColFromL0Flag()? uiCurrRefListId: 1- uiCurrRefListId) : uiCurrRefListId ));

      if ( pcNeighorPU[uiN] && !pcNeighorPU[uiN]->isIntra( uiNeighorPartIdx[uiN] ) )
      {
        if ( pcNeighorPU[uiN]->getInterDir( uiNeighorPartIdx[uiN]  ) & (1<<eCurrRefPicList) )
        {
          pcNeighorPU[uiN]->getMvField( pcNeighorPU[uiN], uiNeighorPartIdx[uiN], eCurrRefPicList, tempMvField);
          iPOC = getSlice()->getRefPic( eCurrRefPicList, tempMvField.getRefIdx())->getPOC() ;
          rcMv  = tempMvField.getMv();
          bTvAva = true;
          return;
        }
      }
    }
}
#endif
#endif

Void TComDataCU::deriveLeftRightTopIdxGeneral ( UInt uiAbsPartIdx, UInt uiPartIdx, UInt& ruiPartIdxLT, UInt& ruiPartIdxRT )
{
#if TT_TREE_DS5 // deriveZIdxInCTU
  ruiPartIdxLT = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  ruiPartIdxLT = m_absZIdxInCtu + uiAbsPartIdx;
#endif
  UInt uiPUWidth = 0;

#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  uiPUWidth = getWidth(uiAbsPartIdx);  
#else
  uiPUWidth = m_puhWidth[getTextType()][uiAbsPartIdx];  
#endif
#else
  switch ( m_pePartSize[uiAbsPartIdx] )
  {
    case SIZE_2Nx2N: uiPUWidth = m_puhWidth[uiAbsPartIdx];  break;
    case SIZE_2NxN:  uiPUWidth = m_puhWidth[uiAbsPartIdx];   break;
    case SIZE_Nx2N:  uiPUWidth = m_puhWidth[uiAbsPartIdx]  >> 1;  break;
    case SIZE_NxN:   uiPUWidth = m_puhWidth[uiAbsPartIdx]  >> 1; break;
    case SIZE_2NxnU:   uiPUWidth = m_puhWidth[uiAbsPartIdx]; break;
    case SIZE_2NxnD:   uiPUWidth = m_puhWidth[uiAbsPartIdx]; break;
    case SIZE_nLx2N:
      if ( uiPartIdx == 0 )
      {
        uiPUWidth = m_puhWidth[uiAbsPartIdx]  >> 2;
      }
      else if ( uiPartIdx == 1 )
      {
        uiPUWidth = (m_puhWidth[uiAbsPartIdx]  >> 1) + (m_puhWidth[uiAbsPartIdx]  >> 2);
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_nRx2N:
      if ( uiPartIdx == 0 )
      {
        uiPUWidth = (m_puhWidth[uiAbsPartIdx]  >> 1) + (m_puhWidth[uiAbsPartIdx]  >> 2);
      }
      else if ( uiPartIdx == 1 )
      {
        uiPUWidth = m_puhWidth[uiAbsPartIdx]  >> 2;
      }
      else
      {
        assert(0);
      }
      break;
    default:
      assert (0);
      break;
  }
#endif

  ruiPartIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ ruiPartIdxLT ] + uiPUWidth / m_pcPic->getMinCUWidth() - 1 ];
}

Void TComDataCU::deriveLeftBottomIdxGeneral( UInt uiAbsPartIdx, UInt uiPartIdx, UInt& ruiPartIdxLB )
{
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  UInt uiPUHeight = getHeight(uiAbsPartIdx); 
#else
  UInt uiPUHeight = m_puhHeight[getTextType()][uiAbsPartIdx]; 
#endif
#else
  UInt uiPUHeight = 0;
  switch ( m_pePartSize[uiAbsPartIdx] )
  {
    case SIZE_2Nx2N: uiPUHeight = m_puhHeight[uiAbsPartIdx];    break;
    case SIZE_2NxN:  uiPUHeight = m_puhHeight[uiAbsPartIdx] >> 1;    break;
    case SIZE_Nx2N:  uiPUHeight = m_puhHeight[uiAbsPartIdx];  break;
    case SIZE_NxN:   uiPUHeight = m_puhHeight[uiAbsPartIdx] >> 1;    break;
    case SIZE_2NxnU:
      if ( uiPartIdx == 0 )
      {
        uiPUHeight = m_puhHeight[uiAbsPartIdx] >> 2;
      }
      else if ( uiPartIdx == 1 )
      {
        uiPUHeight = (m_puhHeight[uiAbsPartIdx] >> 1) + (m_puhHeight[uiAbsPartIdx] >> 2);
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_2NxnD:
      if ( uiPartIdx == 0 )
      {
        uiPUHeight = (m_puhHeight[uiAbsPartIdx] >> 1) + (m_puhHeight[uiAbsPartIdx] >> 2);
      }
      else if ( uiPartIdx == 1 )
      {
        uiPUHeight = m_puhHeight[uiAbsPartIdx] >> 2;
      }
      else
      {
        assert(0);
      }
      break;
    case SIZE_nLx2N: uiPUHeight = m_puhHeight[uiAbsPartIdx];  break;
    case SIZE_nRx2N: uiPUHeight = m_puhHeight[uiAbsPartIdx];  break;
    default:
      assert (0);
      break;
  }
#endif

#if TT_TREE_DS5 // deriveZIdxInCTU
  ruiPartIdxLB      = g_auiRasterToZscan [g_auiZscanToRaster[ deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) ] + ((uiPUHeight / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth()];
#else
  ruiPartIdxLB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu + uiAbsPartIdx ] + ((uiPUHeight / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth()];
#endif
}

Void TComDataCU::deriveLeftRightTopIdx ( UInt uiPartIdx, UInt& ruiPartIdxLT, UInt& ruiPartIdxRT )
{
  ruiPartIdxLT = m_absZIdxInCtu;
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  ruiPartIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ ruiPartIdxLT ] + getWidth(0) / m_pcPic->getMinCUWidth() - 1 ];
#else
  ruiPartIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ ruiPartIdxLT ] + m_puhWidth[getTextType()][0] / m_pcPic->getMinCUWidth() - 1 ];
#endif
#else
  ruiPartIdxRT = g_auiRasterToZscan [g_auiZscanToRaster[ ruiPartIdxLT ] + m_puhWidth[0] / m_pcPic->getMinCUWidth() - 1 ];

  switch ( m_pePartSize[0] )
  {
    case SIZE_2Nx2N:                                                                                                                                break;
    case SIZE_2NxN:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 1; ruiPartIdxRT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 1;
      break;
    case SIZE_Nx2N:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 2; ruiPartIdxRT -= ( uiPartIdx == 1 )? 0 : m_uiNumPartition >> 2;
      break;
    case SIZE_NxN:
      ruiPartIdxLT += ( m_uiNumPartition >> 2 ) * uiPartIdx;         ruiPartIdxRT +=  ( m_uiNumPartition >> 2 ) * ( uiPartIdx - 1 );
      break;
    case SIZE_2NxnU:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 3;
      ruiPartIdxRT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 3;
      break;
    case SIZE_2NxnD:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : ( m_uiNumPartition >> 1 ) + ( m_uiNumPartition >> 3 );
      ruiPartIdxRT += ( uiPartIdx == 0 )? 0 : ( m_uiNumPartition >> 1 ) + ( m_uiNumPartition >> 3 );
      break;
    case SIZE_nLx2N:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 4;
      ruiPartIdxRT -= ( uiPartIdx == 1 )? 0 : ( m_uiNumPartition >> 2 ) + ( m_uiNumPartition >> 4 );
      break;
    case SIZE_nRx2N:
      ruiPartIdxLT += ( uiPartIdx == 0 )? 0 : ( m_uiNumPartition >> 2 ) + ( m_uiNumPartition >> 4 );
      ruiPartIdxRT -= ( uiPartIdx == 1 )? 0 : m_uiNumPartition >> 4;
      break;
    default:
      assert (0);
      break;
  }
#endif

}

Void TComDataCU::deriveLeftBottomIdx( UInt  uiPartIdx,      UInt&      ruiPartIdxLB )
{
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  ruiPartIdxLB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( (getHeight(0) / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth()];
#else
  ruiPartIdxLB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( (m_puhHeight[getTextType()][0] / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth()];
#endif
#else
  ruiPartIdxLB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( ((m_puhHeight[0] / m_pcPic->getMinCUHeight())>>1) - 1)*m_pcPic->getNumPartInCtuWidth()];

  switch ( m_pePartSize[0] )
  {
    case SIZE_2Nx2N:
      ruiPartIdxLB += m_uiNumPartition >> 1;
      break;
    case SIZE_2NxN:
      ruiPartIdxLB += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 1;
      break;
    case SIZE_Nx2N:
      ruiPartIdxLB += ( uiPartIdx == 0 )? m_uiNumPartition >> 1 : (m_uiNumPartition >> 2)*3;
      break;
    case SIZE_NxN:
      ruiPartIdxLB += ( m_uiNumPartition >> 2 ) * uiPartIdx;
      break;
    case SIZE_2NxnU:
      ruiPartIdxLB += ( uiPartIdx == 0 ) ? -((Int)m_uiNumPartition >> 3) : m_uiNumPartition >> 1;
      break;
    case SIZE_2NxnD:
      ruiPartIdxLB += ( uiPartIdx == 0 ) ? (m_uiNumPartition >> 2) + (m_uiNumPartition >> 3): m_uiNumPartition >> 1;
      break;
    case SIZE_nLx2N:
      ruiPartIdxLB += ( uiPartIdx == 0 ) ? m_uiNumPartition >> 1 : (m_uiNumPartition >> 1) + (m_uiNumPartition >> 4);
      break;
    case SIZE_nRx2N:
      ruiPartIdxLB += ( uiPartIdx == 0 ) ? m_uiNumPartition >> 1 : (m_uiNumPartition >> 1) + (m_uiNumPartition >> 2) + (m_uiNumPartition >> 4);
      break;
    default:
      assert (0);
      break;
  }
#endif
}

/** Derive the partition index of neighbouring bottom right block
 * \param [in]  uiPartIdx     current partition index
 * \param [out] ruiPartIdxRB  partition index of neighbouring bottom right block
 */
Void TComDataCU::deriveRightBottomIdx( UInt uiPartIdx, UInt &ruiPartIdxRB )
{
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
  ruiPartIdxRB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( (getHeight(0) / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth() +  getWidth(0) / m_pcPic->getMinCUWidth() - 1];
#else
  ruiPartIdxRB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( (m_puhHeight[getTextType()][0] / m_pcPic->getMinCUHeight()) - 1)*m_pcPic->getNumPartInCtuWidth() +  m_puhWidth[getTextType()][0] / m_pcPic->getMinCUWidth() - 1];
#endif
#else
  ruiPartIdxRB      = g_auiRasterToZscan [g_auiZscanToRaster[ m_absZIdxInCtu ] + ( ((m_puhHeight[0] / m_pcPic->getMinCUHeight())>>1) - 1)*m_pcPic->getNumPartInCtuWidth() +  m_puhWidth[0] / m_pcPic->getMinCUWidth() - 1];

  switch ( m_pePartSize[0] )
  {
    case SIZE_2Nx2N:
      ruiPartIdxRB += m_uiNumPartition >> 1;
      break;
    case SIZE_2NxN:
      ruiPartIdxRB += ( uiPartIdx == 0 )? 0 : m_uiNumPartition >> 1;
      break;
    case SIZE_Nx2N:
      ruiPartIdxRB += ( uiPartIdx == 0 )? m_uiNumPartition >> 2 : (m_uiNumPartition >> 1);
      break;
    case SIZE_NxN:
      ruiPartIdxRB += ( m_uiNumPartition >> 2 ) * ( uiPartIdx - 1 );
      break;
    case SIZE_2NxnU:
      ruiPartIdxRB += ( uiPartIdx == 0 ) ? -((Int)m_uiNumPartition >> 3) : m_uiNumPartition >> 1;
      break;
    case SIZE_2NxnD:
      ruiPartIdxRB += ( uiPartIdx == 0 ) ? (m_uiNumPartition >> 2) + (m_uiNumPartition >> 3): m_uiNumPartition >> 1;
      break;
    case SIZE_nLx2N:
      ruiPartIdxRB += ( uiPartIdx == 0 ) ? (m_uiNumPartition >> 3) + (m_uiNumPartition >> 4): m_uiNumPartition >> 1;
      break;
    case SIZE_nRx2N:
      ruiPartIdxRB += ( uiPartIdx == 0 ) ? (m_uiNumPartition >> 2) + (m_uiNumPartition >> 3) + (m_uiNumPartition >> 4) : m_uiNumPartition >> 1;
      break;
    default:
      assert (0);
      break;
  }
#endif
}

#if TT_TREE_DS1 //deriveZIdxInCTU & deriveZIdxInCU
//NOTE: uiCUZIdxInCTU is the absPartIdx of CU related to CTU and uiPartZIdxInCU is the absPartIdx inside a CU
UInt TComDataCU::deriveZIdxInCTU( UInt uiCUZIdxInCTU, UInt uiPartZIdxInCU)
{
  //decoder case: mode/mv data arranged in a CTUZIdx manner
  if(m_bUseAbsoluteZIdxToCtu) //check if the cu data is pointed to ctu data with a ZIdx offset, which only performs at decoder 
    return uiCUZIdxInCTU + uiPartZIdxInCU;
  
  //encoder case: mode/mv data arranged in a CUZIdx manner
  if(uiCUZIdxInCTU==0 || uiPartZIdxInCU==0)
    return uiCUZIdxInCTU + uiPartZIdxInCU; //this is for speed-up

  UInt uiRasterInCTU = g_auiZscanToRaster[ uiCUZIdxInCTU ];
  UInt uiRasterInCU  = g_auiZscanToRaster[ uiPartZIdxInCU ];
  UInt uiStride      = m_pcPic->getNumPartInCtuWidth();
  UInt uiNewY        = (uiRasterInCTU + uiRasterInCU) / uiStride;
  UInt uiNewX        = (uiRasterInCTU + uiRasterInCU) % uiStride;
  assert(uiNewX<uiStride);
  assert(uiNewY<m_pcPic->getNumPartInCtuHeight());
  return g_auiRasterToZscan[uiNewY * uiStride + uiNewX];
}

UInt TComDataCU::deriveZIdxInCTUNoCheck( UInt uiCUZIdxInCTU, UInt uiPartZIdxInCU)
{
  //encoder case: mode/mv data arranged in a CUZIdx manner
  if(uiCUZIdxInCTU==0 || uiPartZIdxInCU==0)
    return uiCUZIdxInCTU + uiPartZIdxInCU; //this is for speed-up

  UInt uiRasterInCTU = g_auiZscanToRaster[ uiCUZIdxInCTU ];
  UInt uiRasterInCU  = g_auiZscanToRaster[ uiPartZIdxInCU ];
  UInt uiStride      = m_pcPic->getNumPartInCtuWidth();
  UInt uiNewY        = (uiRasterInCTU + uiRasterInCU) / uiStride;
  UInt uiNewX        = (uiRasterInCTU + uiRasterInCU) % uiStride;
  assert(uiNewX<uiStride);
  assert(uiNewY<m_pcPic->getNumPartInCtuHeight());
  return g_auiRasterToZscan[uiNewY * uiStride + uiNewX];
}

UInt TComDataCU::deriveZIdxInCU( UInt uiCUZIdxInCTU, UInt uiPartZIdxInCTU)
{
  //decoder case: mode/mv data arranged in a CTUZIdx manner
  if(m_bUseAbsoluteZIdxToCtu) //check if the cu data is pointed to ctu data with a ZIdx offset, which only performs at decoder
    return uiPartZIdxInCTU - uiCUZIdxInCTU;

  //encoder case: mode/mv data arranged in a CUZIdx manner
  if(uiCUZIdxInCTU==uiPartZIdxInCTU)
    return 0;
  if(uiCUZIdxInCTU==0)
    return uiPartZIdxInCTU; //this is for speed-up

  UInt uiRasterCU    = g_auiZscanToRaster[ uiCUZIdxInCTU ];
  UInt uiRasterPart  = g_auiZscanToRaster[ uiPartZIdxInCTU ];
  UInt uiStride      = m_pcPic->getNumPartInCtuWidth();
  UInt uiNewY        = (uiRasterPart - uiRasterCU) / uiStride;
  UInt uiNewX        = (uiRasterPart - uiRasterCU) % uiStride;
  assert(uiCUZIdxInCTU<=uiPartZIdxInCTU);
  assert(uiNewX>=0);
  assert(uiNewY>=0);
  return g_auiRasterToZscan[uiNewY * uiStride + uiNewX];
}

//because at the decoder side, we need both a CU-ZIdx for pixel buffer addressing and a CTU-ZIdx for cu info
//this function is used for some cases related to pixel buffer addressing
UInt TComDataCU::deriveZIdxInCUNoCheck( UInt uiCUZIdxInCTU, UInt uiPartZIdxInCTU)
{
  if(uiCUZIdxInCTU==uiPartZIdxInCTU)
    return 0;
  if(uiCUZIdxInCTU==0)
    return uiPartZIdxInCTU; //this is for speed-up

  UInt uiRasterCU    = g_auiZscanToRaster[ uiCUZIdxInCTU ];
  UInt uiRasterPart  = g_auiZscanToRaster[ uiPartZIdxInCTU ];
  UInt uiStride      = m_pcPic->getNumPartInCtuWidth();
  UInt uiNewY        = (uiRasterPart - uiRasterCU) / uiStride;
  UInt uiNewX        = (uiRasterPart - uiRasterCU) % uiStride;
  assert(uiCUZIdxInCTU<=uiPartZIdxInCTU);
  assert(uiNewX>=0);
  assert(uiNewY>=0);
  return g_auiRasterToZscan[uiNewY * uiStride + uiNewX];
}
#endif

Bool TComDataCU::hasEqualMotion( UInt uiAbsPartIdx, TComDataCU* pcCandCU, UInt uiCandAbsPartIdx )
{
  if ( getInterDir( uiAbsPartIdx ) != pcCandCU->getInterDir( uiCandAbsPartIdx ) )
  {
    return false;
  }

  for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
  {
    if ( getInterDir( uiAbsPartIdx ) & ( 1 << uiRefListIdx ) )
    {
      if ( getCUMvField( RefPicList( uiRefListIdx ) )->getMv( uiAbsPartIdx )     != pcCandCU->getCUMvField( RefPicList( uiRefListIdx ) )->getMv( uiCandAbsPartIdx ) ||
        getCUMvField( RefPicList( uiRefListIdx ) )->getRefIdx( uiAbsPartIdx ) != pcCandCU->getCUMvField( RefPicList( uiRefListIdx ) )->getRefIdx( uiCandAbsPartIdx ) )
      {
        return false;
      }
    }
  }

  return true;
}

#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
Void TComDataCU::getNeighboringMvField(TComDataCU *pcCU, UInt uiPartIdx, TComMvField *cMvField,UChar *pucInterDir)
{
   Int iRefPOCSrc,iRefPOCMirror;
   RefPicList eRefPicListSrc /*, eRefPicListMirror*/;
   UInt uiMvIdxSrc,uiMvIdxMirror;
   TComMvField cMvFieldTemp;
   if (pcCU->getInterDir(uiPartIdx)==3)
   {
     *pucInterDir=3;
     for (uiMvIdxSrc=0;uiMvIdxSrc<2;uiMvIdxSrc++)
     {
       eRefPicListSrc = (RefPicList)uiMvIdxSrc;
       pcCU->getMvField(pcCU,uiPartIdx,eRefPicListSrc, cMvFieldTemp);
       if (cMvFieldTemp.getRefIdx() ==0)
       {
         cMvField[uiMvIdxSrc] = cMvFieldTemp;
       }
       else
       {
         iRefPOCSrc    = m_pcSlice->getRefPOC( eRefPicListSrc, cMvFieldTemp.getRefIdx() );       
         iRefPOCMirror = m_pcSlice->getRefPOC( eRefPicListSrc, 0 );
         Int iScale = xGetDistScaleFactor( getSlice()->getPOC(), iRefPOCMirror, getSlice()->getPOC(), iRefPOCSrc );
         if ( iScale == 4096 )
         {
           cMvField[uiMvIdxSrc].setMvField(cMvFieldTemp.getMv(), 0);
         }
         else
         {
           cMvField[uiMvIdxSrc].setMvField(cMvFieldTemp.getMv().scaleMv( iScale ), 0);
         }
       }
     }
   }
   else 
   {
     if (pcCU->getInterDir(uiPartIdx)&1)
     {
       eRefPicListSrc    = REF_PIC_LIST_0;
       //eRefPicListMirror = REF_PIC_LIST_1;
       uiMvIdxSrc = 0;
     }
     else
     {
       eRefPicListSrc    = REF_PIC_LIST_1;
       //eRefPicListMirror = REF_PIC_LIST_0;
       uiMvIdxSrc = 1;
     }
     *pucInterDir=uiMvIdxSrc+1;
     uiMvIdxMirror = 1- uiMvIdxSrc;
     pcCU->getMvField(pcCU,uiPartIdx,eRefPicListSrc, cMvFieldTemp);
     iRefPOCSrc    = m_pcSlice->getRefPOC( eRefPicListSrc, cMvFieldTemp.getRefIdx() );       
     if (cMvFieldTemp.getRefIdx() ==0)
     {
       cMvField[uiMvIdxSrc] = cMvFieldTemp;
     }
     else
     {
       iRefPOCMirror = m_pcSlice->getRefPOC( eRefPicListSrc, 0 );
       Int iScale = xGetDistScaleFactor( getSlice()->getPOC(), iRefPOCMirror, getSlice()->getPOC(), iRefPOCSrc );
       if ( iScale == 4096 )
       {
         cMvField[uiMvIdxSrc].setMvField(cMvFieldTemp.getMv(), 0);
       }
       else
       {
         cMvField[uiMvIdxSrc].setMvField(cMvFieldTemp.getMv().scaleMv( iScale ), 0);
       }
     }
     TComMv cZeroMv;
     cZeroMv.setZero();
     cMvField[uiMvIdxMirror].setMvField(cZeroMv,-1);
  }
}

Void TComDataCU::generateMvField(TComMvField *cMvField,UChar* pucInterDir, UInt uiMvNum,TComMvField* cMvFieldMedian,UChar &ucInterDirMedian)
{
  UChar ucDisable = uiMvNum;
  TComMv cMv;
  ucInterDirMedian = 0;

  if (uiMvNum==0)
  {
    if (getSlice()->getSliceType() == P_SLICE)
    {
      ucInterDirMedian = 1;
      cMv.setZero();
      cMvFieldMedian[0].setMvField( cMv ,0);
      cMvFieldMedian[1].setMvField( cMv ,-1);
    }
    else
    {
      ucInterDirMedian = 3;
      cMv.setZero();
      cMvFieldMedian[0].setMvField( cMv ,0);
      cMvFieldMedian[1].setMvField( cMv ,0);
    }
    return;
  }
  for (UInt j=0;j<2;j++)
  {
    Int iExistMvNum =0;
    Int cMvX=0,cMvY=0;
    for (UInt i=0;i<uiMvNum;i++)
    {
      if (pucInterDir[i] & (j+1) && ucDisable!= i)
      {
        cMvX +=  cMvField[ (i<<1)+j ].getHor();
        cMvY +=  cMvField[ (i<<1)+j ].getVer();
        iExistMvNum++;
      }
    }
    if (iExistMvNum)
    {
      ucInterDirMedian |= (j+1);
      if (iExistMvNum==3)
      {
        cMv.set((Short)(cMvX*43/128), (Short)(cMvY *43/128));
      }
      else if (iExistMvNum==2)
      {
        cMv.set((Short)(cMvX/2), (Short)(cMvY/2));
      }
      else
      {
        cMv.set((Short)(cMvX), (Short)(cMvY));
      }
      cMvFieldMedian[j].setMvField(cMv,0);
    }
    else
    {
      cMv.setZero();
      cMvFieldMedian[j].setMvField(cMv,-1);
    }
  }
}




Bool TComDataCU::getInterMergeSubPURecursiveCandidate( UInt uiAbsPartIdx, UInt uiPUIdx, TComMvField* pcMvFieldNeighbours, UChar* puhInterDirNeighbours, Int& numValidMergeCand
  , UChar*          peMergeTypeNeighbors  , 
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComMvField* pcMvFieldSP[NUM_MGR_TYPE] , UChar* puhInterDirSP[NUM_MGR_TYPE] ,
#else
  TComMvField* pcMvFieldSP[2] , UChar* puhInterDirSP[2] ,
#endif
  Int iCount )
{
    Bool bAtmvpAva=false;
    Int iPartitionIdx     = 0;
    // compute the location of the current PU
    Int iCurrPosX, iCurrPosY;
    Int iWidth, iHeight;
    Int iPUWidth, iPUHeight, iNumPart, iNumPartLine;

    getPartPosition( uiPUIdx, iCurrPosX, iCurrPosY,      iWidth, iHeight);
    getSPPara      ( iWidth,  iHeight,  iNumPart, iNumPartLine, iPUWidth, iPUHeight);
#if JVET_C0035_ATMVP_SIMPLIFICATION
    UInt uiCurMRGType = MGR_TYPE_SUBPU_ATMVP_EXT;
#endif  
    UInt uiSameCount=0;
    UInt uiSameCountATMVP=0;
    for (Int i=0; i <  iHeight; i += iPUHeight)
    {
      for (Int j = 0; j <  iWidth; j += iPUWidth)
      {
          TComMvField cMvField[6],cMvFieldMedian[2];
          UChar ucInterDir[3],ucInterDirMedian=0;
          UInt uiMVCount=0;
          UInt uiSPAddr=0;
          TComDataCU* pcCULeftMedian= 0,*pcCUAboveMedian= 0;
          UInt uiLeftPartIdx_median=0, uiAbovePartIdx_median=0;
          getSPAbsPartIdx(uiAbsPartIdx, iPUWidth, iPUHeight, iPartitionIdx, iNumPartLine, uiSPAddr);
#if TT_TREE_ATMVP // deriveZIdxInCTU
          uiSPAddr  = deriveZIdxInCTU(m_absZIdxInCtu, uiSPAddr);
#else
          uiSPAddr += m_absZIdxInCtu;
#endif
          //get left
          if (iPartitionIdx%iNumPartLine ==0)
          {
              for (UInt uiCurAddrY = i/iPUHeight; uiCurAddrY <iHeight/iPUHeight ; uiCurAddrY++)
              {
                UInt uiSPAddrCur=0;
                getSPAbsPartIdx(uiAbsPartIdx, iPUWidth, iPUHeight, uiCurAddrY* iNumPartLine, iNumPartLine, uiSPAddrCur);
#if TT_TREE_ATMVP // deriveZIdxInCTU
                uiSPAddrCur  = deriveZIdxInCTU(m_absZIdxInCtu, uiSPAddrCur);
#else
                uiSPAddrCur += m_absZIdxInCtu;
#endif

                pcCULeftMedian = getPULeft( uiLeftPartIdx_median, uiSPAddrCur );
                if ( pcCULeftMedian &&!pcCULeftMedian->isIntra( uiLeftPartIdx_median ))
                {
                  getNeighboringMvField(pcCULeftMedian,uiLeftPartIdx_median,cMvField,ucInterDir);
                  uiMVCount++;
                  break;
                }
              }
          }
          else
          {
#if JVET_C0035_ATMVP_SIMPLIFICATION
            ucInterDir[0]= puhInterDirSP[uiCurMRGType][iPartitionIdx-1];
            cMvField[0]=pcMvFieldSP[uiCurMRGType][((iPartitionIdx-1)<<1)  ];
            cMvField[1]=pcMvFieldSP[uiCurMRGType][((iPartitionIdx-1)<<1)+1];
#else
            ucInterDir[0]= puhInterDirSP[1][iPartitionIdx-1];
            cMvField[0]=pcMvFieldSP[1][((iPartitionIdx-1)<<1)  ];
            cMvField[1]=pcMvFieldSP[1][((iPartitionIdx-1)<<1)+1];
#endif
            uiMVCount++;
          }
          //get above
          if (iPartitionIdx < iNumPartLine)
          {
              for (UInt uiCurAddrX = iPartitionIdx; uiCurAddrX <iNumPartLine ; uiCurAddrX++)
              {
                UInt uiSPAddrCur=0;
                getSPAbsPartIdx(uiAbsPartIdx, iPUWidth, iPUHeight, uiCurAddrX, iNumPartLine, uiSPAddrCur);
#if TT_TREE_ATMVP // deriveZIdxInCTU
                uiSPAddrCur  = deriveZIdxInCTU(m_absZIdxInCtu, uiSPAddrCur);
#else
                uiSPAddrCur += m_absZIdxInCtu;
#endif
                pcCUAboveMedian = getPUAbove( uiAbovePartIdx_median, uiSPAddrCur );
                if (pcCUAboveMedian &&!pcCUAboveMedian->isIntra( uiAbovePartIdx_median ))
                {
                  getNeighboringMvField(pcCUAboveMedian,uiAbovePartIdx_median,cMvField+(uiMVCount<<1),ucInterDir+uiMVCount);
                  uiMVCount++;
                  break;
                }
              }
          }
          else
          {
#if JVET_C0035_ATMVP_SIMPLIFICATION
            ucInterDir[uiMVCount]= puhInterDirSP[uiCurMRGType][iPartitionIdx-iNumPartLine];
            cMvField[(uiMVCount<<1)  ]=pcMvFieldSP[uiCurMRGType][((iPartitionIdx-iNumPartLine)<<1)  ];
            cMvField[(uiMVCount<<1)+1]=pcMvFieldSP[uiCurMRGType][((iPartitionIdx-iNumPartLine)<<1)+1];
#else
            ucInterDir[uiMVCount]= puhInterDirSP[1][iPartitionIdx-iNumPartLine];
            cMvField[(uiMVCount<<1)  ]=pcMvFieldSP[1][((iPartitionIdx-iNumPartLine)<<1)  ];
            cMvField[(uiMVCount<<1)+1]=pcMvFieldSP[1][((iPartitionIdx-iNumPartLine)<<1)+1];
#endif
            uiMVCount++;
          }

          {
            Bool bExistMV = false;
            ucInterDir[uiMVCount] = 0;
            if ( getSlice()->getEnableTMVPFlag())
            {
              Int iRefIdx = 0;
              TComMv cColMvTemp;
              UInt uiAbsPartIdxSP = g_auiZscanToRaster[uiSPAddr];
              UInt uiAbsPartSPAddrRB = uiSPAddr;
              //UInt uiNumPartInCUWidth = m_pcPic->getNumPartInCtuWidth();
              const UInt numPartInCtuWidth  = m_pcPic->getNumPartInCtuWidth();
              const UInt numPartInCtuHeight = m_pcPic->getNumPartInCtuHeight();
              Int uiLCUIdx = getCtuRsAddr();
              
              UInt uiPartAddrCenter = uiSPAddr;
              uiPartAddrCenter = g_auiRasterToZscan[ g_auiZscanToRaster[ uiPartAddrCenter ]
                                        + ( iPUHeight/m_pcPic->getMinCUHeight()  )/2*numPartInCtuWidth
                                        + ( iPUWidth/m_pcPic->getMinCUWidth()  )/2];

              if      ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + g_auiRasterToPelX[uiAbsPartIdxSP] + m_pcPic->getMinCUWidth() ) >= m_pcSlice->getSPS()->getPicWidthInLumaSamples() )  // image boundary check
              {
              }
              else if ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + g_auiRasterToPelY[uiAbsPartIdxSP] + m_pcPic->getMinCUHeight() ) >= m_pcSlice->getSPS()->getPicHeightInLumaSamples() )
              {
              }
              else
              {
                if ( ( uiAbsPartIdxSP % numPartInCtuWidth < numPartInCtuWidth - 1 ) &&           // is not at the last column of LCU 
                     ( uiAbsPartIdxSP / numPartInCtuWidth < numPartInCtuHeight - 1 ) ) // is not at the last row    of LCU
                {
                  uiAbsPartSPAddrRB = g_auiRasterToZscan[ uiAbsPartIdxSP + numPartInCtuWidth + 1 ];
                  uiLCUIdx = getCtuRsAddr();
                }
                else if ( uiAbsPartIdxSP % numPartInCtuWidth < numPartInCtuWidth - 1 )           // is not at the last column of LCU But is last row of LCU
                {
#if COM16_C806_GEN_MRG_IMPROVEMENT
                  uiLCUIdx = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus();
#endif
                  uiAbsPartSPAddrRB = g_auiRasterToZscan[ (uiAbsPartIdxSP + numPartInCtuWidth + 1) % m_pcPic->getNumPartitionsInCtu() ];
                }
                else if ( uiAbsPartIdxSP / numPartInCtuWidth < numPartInCtuHeight - 1 ) // is not at the last row of LCU But is last column of LCU
                {
                  uiAbsPartSPAddrRB = g_auiRasterToZscan[ uiAbsPartIdxSP + 1 ];
                  uiLCUIdx = getCtuRsAddr() + 1;
                }
                else //is the right bottom corner of LCU                       
                {
                  uiAbsPartSPAddrRB = 0;
#if COM16_C806_GEN_MRG_IMPROVEMENT
                  uiLCUIdx = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus() + 1;
#endif
                }
              }
              bExistMV = uiLCUIdx >= 0 && xGetColMVP( REF_PIC_LIST_0, uiLCUIdx, uiAbsPartSPAddrRB, cColMvTemp, iRefIdx );
              if( bExistMV == false )
              {
                bExistMV = xGetColMVP( REF_PIC_LIST_0, getCtuRsAddr(), uiPartAddrCenter, cColMvTemp, iRefIdx );
              }
              if( bExistMV )
              {
                ucInterDir[uiMVCount] |= 1;
                cMvField[ 2 * uiMVCount ].setMvField( cColMvTemp, iRefIdx );
              }
    
              if ( getSlice()->isInterB() )
              {
                bExistMV = uiLCUIdx >= 0 && xGetColMVP( REF_PIC_LIST_1, uiLCUIdx, uiAbsPartSPAddrRB, cColMvTemp, iRefIdx);
                if( bExistMV == false )
                {
                  bExistMV = xGetColMVP( REF_PIC_LIST_1, getCtuRsAddr(), uiPartAddrCenter, cColMvTemp, iRefIdx );
                }
                if( bExistMV )
                {
                  ucInterDir[uiMVCount] |= 2;
                  cMvField[ 2 * uiMVCount + 1 ].setMvField( cColMvTemp, iRefIdx );
                }
              }
            }
            if(ucInterDir[uiMVCount]>0)
            {  
              uiMVCount++;
            }
          }
          generateMvField(   cMvField, ucInterDir, uiMVCount,cMvFieldMedian,ucInterDirMedian);

#if JVET_C0035_ATMVP_SIMPLIFICATION
          puhInterDirSP[uiCurMRGType][iPartitionIdx] = ucInterDirMedian;
          pcMvFieldSP[uiCurMRGType][(iPartitionIdx<<1)  ] = cMvFieldMedian[0];
          pcMvFieldSP[uiCurMRGType][(iPartitionIdx<<1)+1] = cMvFieldMedian[1];

          if (iPartitionIdx == 0  || 
            (//uiSameCount == iPartitionIdx &&
            puhInterDirSP[uiCurMRGType][iPartitionIdx]      == puhInterDirSP[uiCurMRGType][0] && 
            pcMvFieldSP[uiCurMRGType][(iPartitionIdx<<1)  ] == pcMvFieldSP[uiCurMRGType][0] &&
            pcMvFieldSP[uiCurMRGType][(iPartitionIdx<<1)+1] == pcMvFieldSP[uiCurMRGType][1] ))
          {
            uiSameCount++;
          }
          if ( uiSameCountATMVP == iPartitionIdx && 
            puhInterDirSP[MGR_TYPE_SUBPU_ATMVP_EXT][iPartitionIdx]      == puhInterDirSP[MGR_TYPE_SUBPU_ATMVP][iPartitionIdx] &&
            pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][(iPartitionIdx<<1)  ] == pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP][((iPartitionIdx)<<1)  ] &&
            pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][(iPartitionIdx<<1)+1] == pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP][((iPartitionIdx)<<1)+1] )              
          {
            uiSameCountATMVP++;
          }
#else
          puhInterDirSP[1][iPartitionIdx] = ucInterDirMedian;
          pcMvFieldSP[1][(iPartitionIdx<<1)  ] = cMvFieldMedian[0];
          pcMvFieldSP[1][(iPartitionIdx<<1)+1] = cMvFieldMedian[1];

          if (iPartitionIdx == 0  || 
            (uiSameCount == iPartitionIdx &&
            puhInterDirSP[1][iPartitionIdx]      == puhInterDirSP[1][0] && 
            pcMvFieldSP[1][(iPartitionIdx<<1)  ] == pcMvFieldSP[1][0] &&
            pcMvFieldSP[1][(iPartitionIdx<<1)+1] == pcMvFieldSP[1][1] ))
          {
            uiSameCount++;
          }
          if ( uiSameCountATMVP == iPartitionIdx && 
            puhInterDirSP[1][iPartitionIdx]      == puhInterDirSP[0][iPartitionIdx] &&
            pcMvFieldSP[1][(iPartitionIdx<<1)  ] == pcMvFieldSP[0][((iPartitionIdx)<<1)  ] &&
            pcMvFieldSP[1][(iPartitionIdx<<1)+1] == pcMvFieldSP[0][((iPartitionIdx)<<1)+1] )              
          {
             uiSameCountATMVP++;
          }
#endif

          iPartitionIdx++;
        }
      }
      Bool bAtmvpExtAva =true;
      if (uiSameCount == iNumPart)
      {
        for (UInt uiIdx=0; uiIdx <iCount;uiIdx++)
        {
#if JVET_C0035_ATMVP_SIMPLIFICATION
          if (peMergeTypeNeighbors[uiIdx]  !=  MGR_TYPE_SUBPU_ATMVP) 
#else
          if (peMergeTypeNeighbors[uiIdx]  !=  MGR_TYPE_SUBPU_TMVP) 
#endif
          {
#if JVET_C0035_ATMVP_SIMPLIFICATION
            if (puhInterDirNeighbours[uiIdx] == puhInterDirSP[uiCurMRGType][0] &&
              pcMvFieldNeighbours[uiIdx<<1]      == pcMvFieldSP[uiCurMRGType][0] &&
              pcMvFieldNeighbours[(uiIdx<<1)+1]  == pcMvFieldSP[uiCurMRGType][1])
#else
            if (puhInterDirNeighbours[uiIdx] == puhInterDirSP[1][0] &&
                pcMvFieldNeighbours[uiIdx<<1]      == pcMvFieldSP[1][0] &&
                pcMvFieldNeighbours[(uiIdx<<1)+1]  == pcMvFieldSP[1][1])
#endif
            {
              bAtmvpExtAva = false;
              break;
            }
          }
        }
      }
      if(bAtmvpExtAva && bAtmvpAva)
      {
        if(uiSameCountATMVP == iNumPart)
        {
           bAtmvpExtAva = false;
        }
      }
      return bAtmvpExtAva;
}

#endif

//! Construct a list of merging candidates
Void TComDataCU::getInterMergeCandidates( UInt uiAbsPartIdx, UInt uiPUIdx, TComMvField* pcMvFieldNeighbours, UChar* puhInterDirNeighbours, Int& numValidMergeCand
#if VCEG_AZ06_IC
  , Bool*           pbICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  , UChar*          peMergeTypeNeighbors
#if JVET_C0035_ATMVP_SIMPLIFICATION
  , TComMvField*    pcMvFieldSP[NUM_MGR_TYPE]
  , UChar*          puhInterDirSP[NUM_MGR_TYPE]
#else
  , TComMvField*    pcMvFieldSP[2]
  , UChar*          puhInterDirSP[2]
#endif
  , UInt            uiDecCurrAbsPartIdx
  , TComDataCU*     pDecCurrCU
#endif
  , Int mrgCandIdx
  )
{
  UInt uiAbsPartAddr = m_absZIdxInCtu + uiAbsPartIdx;
  Bool abCandIsInter[ MRG_MAX_NUM_CANDS ];
  for( UInt ui = 0; ui < getSlice()->getMaxNumMergeCand(); ++ui )
  {
    abCandIsInter[ui] = false;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    pcMvFieldNeighbours[ ( ui << 1 )     ].setMvField(TComMv(0,0), NOT_VALID);
    pcMvFieldNeighbours[ ( ui << 1 ) + 1 ].setMvField(TComMv(0,0), NOT_VALID);
#else
    pcMvFieldNeighbours[ ( ui << 1 )     ].setRefIdx(NOT_VALID);
    pcMvFieldNeighbours[ ( ui << 1 ) + 1 ].setRefIdx(NOT_VALID);
#endif
  }
#if VCEG_AZ06_IC
  memset( pbICFlag, false, sizeof( Bool )*MRG_MAX_NUM_CANDS );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  Bool bICFlag = false;
#endif
#endif
  numValidMergeCand = getSlice()->getMaxNumMergeCand();
  // compute the location of the current PU
  Int xP, yP, nPSW, nPSH;
  this->getPartPosition(uiPUIdx, xP, yP, nPSW, nPSH);

  Int iCount = 0;

  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
#if JVET_C0024_QTBT
  assert(uiPUIdx==0);
#else
  PartSize cCurPS = getPartitionSize( uiAbsPartIdx );
#endif
  deriveLeftRightTopIdxGeneral( uiAbsPartIdx, uiPUIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdxGeneral( uiAbsPartIdx, uiPUIdx, uiPartIdxLB );
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  Bool bEnableATMVP = getSlice()->getSPS()->getAtmvpEnableFlag();
#if JVET_C0035_ATMVP_SIMPLIFICATION
  for (Int i3=0; i3< NUM_MGR_TYPE; i3++)
#else
  for (Int i3=0; i3< 2; i3++)
#endif
  {
    memset( puhInterDirSP[i3] , 0 , getTotalNumPart() * sizeof( UChar ) );
    for (Int i=0 ; i <getTotalNumPart(); i++)
    {
      pcMvFieldSP[i3][2*i].setMvField(TComMv(0,0), NOT_VALID);
      pcMvFieldSP[i3][2*i+1].setMvField(TComMv(0,0), NOT_VALID);
    }
  }
  memset( peMergeTypeNeighbors, MGR_TYPE_DEFAULT_N, sizeof( UChar ) * MRG_MAX_NUM_CANDS );
#if JVET_C0035_ATMVP_SIMPLIFICATION
  Bool bAtmvpAva = false;
  Int  uiAtmvpPos = -1;
#else
  Bool bAtmvpAva = false,  iAvaNBTV  =  false;
  Int  iPOCAtmvp = 0, uiAtmvpPos = -1;
  TComMv cNBTV; 
#endif 
#endif 
  //left
  UInt uiLeftPartIdx = 0;
  TComDataCU* pcCULeft = 0;
  pcCULeft = getPULeft( uiLeftPartIdx, uiPartIdxLB );

  Bool isAvailableA1 = pcCULeft &&
                       pcCULeft->isDiffMER(xP -1, yP+nPSH-1, xP, yP) &&
#if !JVET_C0024_QTBT
                       !( uiPUIdx == 1 && (cCurPS == SIZE_Nx2N || cCurPS == SIZE_nLx2N || cCurPS == SIZE_nRx2N) ) &&
#endif
                       pcCULeft->isInter( uiLeftPartIdx ) ;

  if ( isAvailableA1 )
  {
    abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
    pbICFlag[iCount] = pcCULeft->getICFlag( uiLeftPartIdx );
#endif
    // get Inter Dir
    puhInterDirNeighbours[iCount] = pcCULeft->getInterDir( uiLeftPartIdx );
    // get Mv from Left
    pcCULeft->getMvField( pcCULeft, uiLeftPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbours[iCount<<1] );
    if ( getSlice()->isInterB() )
    {
      pcCULeft->getMvField( pcCULeft, uiLeftPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbours[(iCount<<1)+1] );
    }
    if ( mrgCandIdx == iCount )
    {
      return;
    }
    iCount ++;
  }

  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }
  // above
  UInt uiAbovePartIdx = 0;
  TComDataCU* pcCUAbove = 0;
  pcCUAbove = getPUAbove( uiAbovePartIdx, uiPartIdxRT );

  Bool isAvailableB1 = pcCUAbove &&
                       pcCUAbove->isDiffMER(xP+nPSW-1, yP-1, xP, yP) &&
#if !JVET_C0024_QTBT
                       !( uiPUIdx == 1 && (cCurPS == SIZE_2NxN || cCurPS == SIZE_2NxnU || cCurPS == SIZE_2NxnD) ) &&
#endif
                       pcCUAbove->isInter( uiAbovePartIdx );
#if VCEG_AZ06_IC
  if ( isAvailableB1 && (!isAvailableA1 || !( pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCUAbove, uiAbovePartIdx ) && pcCULeft->getICFlag( uiLeftPartIdx ) == pcCUAbove->getICFlag( uiAbovePartIdx ) ) ) )
#else
  if ( isAvailableB1 && (!isAvailableA1 || !pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCUAbove, uiAbovePartIdx ) ) )
#endif
  {
    abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
    pbICFlag[iCount] = pcCUAbove->getICFlag( uiAbovePartIdx );
#endif
    // get Inter Dir
    puhInterDirNeighbours[iCount] = pcCUAbove->getInterDir( uiAbovePartIdx );
    // get Mv from Left
    pcCUAbove->getMvField( pcCUAbove, uiAbovePartIdx, REF_PIC_LIST_0, pcMvFieldNeighbours[iCount<<1] );
    if ( getSlice()->isInterB() )
    {
      pcCUAbove->getMvField( pcCUAbove, uiAbovePartIdx, REF_PIC_LIST_1, pcMvFieldNeighbours[(iCount<<1)+1] );
    }
    if ( mrgCandIdx == iCount )
    {
      return;
    }
    iCount ++;
  }
  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }

  // above right
  UInt uiAboveRightPartIdx = 0;
  TComDataCU* pcCUAboveRight = 0;
  pcCUAboveRight = getPUAboveRight( uiAboveRightPartIdx, uiPartIdxRT );

  Bool isAvailableB0 = pcCUAboveRight &&
                       pcCUAboveRight->isDiffMER(xP+nPSW, yP-1, xP, yP) &&
                       pcCUAboveRight->isInter( uiAboveRightPartIdx );
#if COM16_C806_GEN_MRG_IMPROVEMENT
#if VCEG_AZ06_IC
  if ( isAvailableB0 && ( !isAvailableB1 || !( pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCUAboveRight, uiAboveRightPartIdx ) && pcCUAbove->getICFlag( uiAbovePartIdx ) == pcCUAboveRight->getICFlag( uiAboveRightPartIdx ) ) ) 
                     && ( !isAvailableA1 || !( pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCUAboveRight, uiAboveRightPartIdx ) && pcCULeft->getICFlag( uiLeftPartIdx ) == pcCUAboveRight->getICFlag( uiAboveRightPartIdx ) ) ) )
#else
  if ( isAvailableB0 && ( !isAvailableB1 || !pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCUAboveRight, uiAboveRightPartIdx ) ) 
                     && ( !isAvailableA1 || !pcCULeft->hasEqualMotion( uiLeftPartIdx , pcCUAboveRight, uiAboveRightPartIdx ) ))
#endif
#else
#if VCEG_AZ06_IC
  if ( isAvailableB0 && ( !isAvailableB1 || !( pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCUAboveRight, uiAboveRightPartIdx ) && pcCUAbove->getICFlag( uiAbovePartIdx ) == pcCUAboveRight->getICFlag( uiAboveRightPartIdx ) ) ) )
#else
  if ( isAvailableB0 && ( !isAvailableB1 || !pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCUAboveRight, uiAboveRightPartIdx ) ) )
#endif
#endif
  {
    abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
    pbICFlag[iCount] = pcCUAboveRight->getICFlag( uiAboveRightPartIdx );
#endif
    // get Inter Dir
    puhInterDirNeighbours[iCount] = pcCUAboveRight->getInterDir( uiAboveRightPartIdx );
    // get Mv from Left
    pcCUAboveRight->getMvField( pcCUAboveRight, uiAboveRightPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbours[iCount<<1] );
    if ( getSlice()->isInterB() )
    {
      pcCUAboveRight->getMvField( pcCUAboveRight, uiAboveRightPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbours[(iCount<<1)+1] );
    }
    if ( mrgCandIdx == iCount )
    {
      return;
    }
    iCount ++;
  }
  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }

  //left bottom
  UInt uiLeftBottomPartIdx = 0;
  TComDataCU* pcCULeftBottom = 0;
  pcCULeftBottom = this->getPUBelowLeft( uiLeftBottomPartIdx, uiPartIdxLB );

  Bool isAvailableA0 = pcCULeftBottom &&
                       pcCULeftBottom->isDiffMER(xP-1, yP+nPSH, xP, yP) &&
                       pcCULeftBottom->isInter( uiLeftBottomPartIdx ) ;
#if COM16_C806_GEN_MRG_IMPROVEMENT
#if VCEG_AZ06_IC
  if ( isAvailableA0 
    && ( !isAvailableA1 || !( pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) && pcCULeft->getICFlag( uiLeftPartIdx ) == pcCULeftBottom->getICFlag( uiLeftBottomPartIdx ) ) ) 
    && ( !isAvailableB1 || !( pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCULeftBottom, uiLeftBottomPartIdx  ) && pcCUAbove->getICFlag( uiAbovePartIdx ) == pcCULeftBottom->getICFlag( uiLeftBottomPartIdx ) ) ) 
    && ( !isAvailableB0 || !( pcCUAboveRight->hasEqualMotion( uiAboveRightPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) && pcCUAboveRight->getICFlag( uiAboveRightPartIdx ) == pcCULeftBottom->getICFlag( uiLeftBottomPartIdx ) ) ))
#else
  if ( isAvailableA0 && ( !isAvailableA1 || !pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) ) && ( !isAvailableB1 || !pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCULeftBottom, uiLeftBottomPartIdx  ) ) && (!isAvailableB0 || !pcCUAboveRight->hasEqualMotion( uiAboveRightPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) ))
#endif
#else
#if VCEG_AZ06_IC
  if ( isAvailableA0 && ( !isAvailableA1 || !( pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) && pcCULeft->getICFlag( uiLeftPartIdx ) == pcCULeftBottom->getICFlag( uiLeftBottomPartIdx ) ) ) )
#else
  if ( isAvailableA0 && ( !isAvailableA1 || !pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCULeftBottom, uiLeftBottomPartIdx ) ) )
#endif
#endif
  {
    abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
    pbICFlag[iCount] = pcCULeftBottom->getICFlag( uiLeftBottomPartIdx );
#endif
    // get Inter Dir
    puhInterDirNeighbours[iCount] = pcCULeftBottom->getInterDir( uiLeftBottomPartIdx );
    // get Mv from Left
    pcCULeftBottom->getMvField( pcCULeftBottom, uiLeftBottomPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbours[iCount<<1] );
    if ( getSlice()->isInterB() )
    {
      pcCULeftBottom->getMvField( pcCULeftBottom, uiLeftBottomPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbours[(iCount<<1)+1] );
    }
    if ( mrgCandIdx == iCount )
    {
      return;
    }
    iCount ++;
  }
  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  if(bEnableATMVP && getSlice()->getEnableTMVPFlag())
  {
    Bool bMrgIdxMatchATMVPCan = (mrgCandIdx == iCount);

#if JVET_C0035_ATMVP_SIMPLIFICATION
    bAtmvpAva = getInterMergeSubPUTmvpCandidate( uiAbsPartIdx,  uiPUIdx,  iCount, pcMvFieldNeighbours, puhInterDirNeighbours,pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP], puhInterDirSP[MGR_TYPE_SUBPU_ATMVP], &pcMvFieldNeighbours[iCount<<1], &puhInterDirNeighbours[iCount], bMrgIdxMatchATMVPCan,
#if VCEG_AZ06_IC
      bICFlag,
#endif
      pDecCurrCU, uiDecCurrAbsPartIdx);
#else
    get1stTvFromSpatialNeighbor(uiAbsPartIdx,uiPUIdx,iAvaNBTV, iPOCAtmvp, cNBTV);
  
    ///////////////////////////////////////////////////////////////////////
    //////GET Motion for the ATMVP candidate based on Temporal Vector//////
    ///////////////////////////////////////////////////////////////////////
    bAtmvpAva = getInterMergeSubPUTmvpCandidate( uiPUIdx, pcMvFieldSP[0], puhInterDirSP[0], &pcMvFieldNeighbours[iCount<<1], &puhInterDirNeighbours[iCount], cNBTV, bMrgIdxMatchATMVPCan,
#if VCEG_AZ06_IC
      bICFlag,
#endif
      iPOCAtmvp, pDecCurrCU, uiDecCurrAbsPartIdx);
#endif

    if ( bAtmvpAva ) 
    {
      abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
      pbICFlag[iCount] = bICFlag;
#endif
#if JVET_C0035_ATMVP_SIMPLIFICATION
      peMergeTypeNeighbors[iCount] = MGR_TYPE_SUBPU_ATMVP;
#else
      peMergeTypeNeighbors[iCount] = MGR_TYPE_SUBPU_TMVP;
#endif
      if ( bMrgIdxMatchATMVPCan )
      {
        return;
      }
      uiAtmvpPos = iCount;
      iCount ++;
      if (iCount == getSlice()->getMaxNumMergeCand()) 
      {
        return;
      }      
    }
  }
  if(bEnableATMVP && getSlice()->getEnableTMVPFlag())
  {
    Bool bAtmvpExtAva = getInterMergeSubPURecursiveCandidate(  uiAbsPartIdx,  uiPUIdx, pcMvFieldNeighbours, puhInterDirNeighbours, numValidMergeCand, peMergeTypeNeighbors  , pcMvFieldSP , puhInterDirSP,  iCount );
    if (bAtmvpExtAva)
    {
        //store the sub-PU motion information
        Int iCurrPosX, iCurrPosY;
        Int iWidth, iHeight;
        Int iPUWidth, iPUHeight, iNumPart, iNumPartLine,iPartitionIdx;
        getPartPosition( uiPUIdx, iCurrPosX, iCurrPosY,      iWidth, iHeight);
        getSPPara      ( iWidth,  iHeight,  iNumPart, iNumPartLine, iPUWidth, iPUHeight);

#if JVET_C0035_ATMVP_SIMPLIFICATION 
        peMergeTypeNeighbors[iCount]  =  MGR_TYPE_SUBPU_ATMVP_EXT;
        puhInterDirNeighbours[iCount] = puhInterDirSP[MGR_TYPE_SUBPU_ATMVP_EXT][iNumPart-1];
        pcMvFieldNeighbours[(iCount<<1)]    =pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][ (iNumPart-1)<<1   ];
        pcMvFieldNeighbours[(iCount<<1)+1]  =pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][((iNumPart-1)<<1)+1];
#else
        peMergeTypeNeighbors[iCount]  =  MGR_TYPE_SUBPU_TMVP_EXT;
        puhInterDirNeighbours[iCount] = puhInterDirSP[1][iNumPart-1];
        pcMvFieldNeighbours[(iCount<<1)]    =pcMvFieldSP[1][ (iNumPart-1)<<1   ];
        pcMvFieldNeighbours[(iCount<<1)+1]  =pcMvFieldSP[1][((iNumPart-1)<<1)+1];
#endif
        abCandIsInter[iCount] = true;

#if VCEG_AZ06_IC
        if( getSlice()->getApplyIC() )
        {
          pbICFlag[iCount] = bAtmvpAva ? !bICFlag : false;
        }
#endif
        if ( mrgCandIdx == iCount )
        {
          for (iPartitionIdx = 0; iPartitionIdx < iNumPart; iPartitionIdx++)
          {
            UInt uiSPAddr;
            pDecCurrCU->getSPAbsPartIdx(uiDecCurrAbsPartIdx, iPUWidth, iPUHeight, iPartitionIdx, iNumPartLine, uiSPAddr);
#if JVET_C0035_ATMVP_SIMPLIFICATION 
            pDecCurrCU->setInterDirSP(puhInterDirSP[MGR_TYPE_SUBPU_ATMVP_EXT][iPartitionIdx], uiSPAddr, iPUWidth, iPUHeight);
            pDecCurrCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][2*iPartitionIdx], iPUWidth, iPUHeight);
            pDecCurrCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[MGR_TYPE_SUBPU_ATMVP_EXT][2*iPartitionIdx + 1], iPUWidth, iPUHeight);
#else
            pDecCurrCU->setInterDirSP(puhInterDirSP[1][iPartitionIdx], uiSPAddr, iPUWidth, iPUHeight);
            pDecCurrCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[1][2*iPartitionIdx], iPUWidth, iPUHeight);
            pDecCurrCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[1][2*iPartitionIdx + 1], iPUWidth, iPUHeight);
#endif
          }
          return;
        }
        iCount ++;
        if (iCount == getSlice()->getMaxNumMergeCand()) 
        {
          return;
        }
      }
  }
  // above left
  if( iCount < (bEnableATMVP? 6: 4) ) 
#else
  // above left
  if( iCount < 4 )
#endif
  {
    UInt uiAboveLeftPartIdx = 0;
    TComDataCU* pcCUAboveLeft = 0;
    pcCUAboveLeft = getPUAboveLeft( uiAboveLeftPartIdx, uiAbsPartAddr );

    Bool isAvailableB2 = pcCUAboveLeft &&
                         pcCUAboveLeft->isDiffMER(xP-1, yP-1, xP, yP) &&
                         pcCUAboveLeft->isInter( uiAboveLeftPartIdx );
#if COM16_C806_GEN_MRG_IMPROVEMENT
#if VCEG_AZ06_IC
    if ( isAvailableB2 
      && ( !isAvailableA1 || !( pcCULeft       ->hasEqualMotion(      uiLeftPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) && pcCULeft->getICFlag( uiLeftPartIdx )   == pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx ) ) )
      && ( !isAvailableB1 || !( pcCUAbove      ->hasEqualMotion(     uiAbovePartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) && pcCUAbove->getICFlag( uiAbovePartIdx ) == pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx ) ) )
      && ( !isAvailableA0 || !( pcCULeftBottom ->hasEqualMotion(uiLeftBottomPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) && pcCULeftBottom->getICFlag( uiLeftBottomPartIdx ) == pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx ) ) )
      && ( !isAvailableB0 || !( pcCUAboveRight ->hasEqualMotion(uiAboveRightPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) && pcCUAboveRight->getICFlag( uiAboveRightPartIdx ) == pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx ) ) )
      )
#else
    if ( isAvailableB2 
        && ( !isAvailableA1 || !pcCULeft       ->hasEqualMotion(      uiLeftPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) )
        && ( !isAvailableB1 || !pcCUAbove      ->hasEqualMotion(     uiAbovePartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) )
        && ( !isAvailableA0 || !pcCULeftBottom ->hasEqualMotion(uiLeftBottomPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) )
        && ( !isAvailableB0 || !pcCUAboveRight ->hasEqualMotion(uiAboveRightPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) )
       )
#endif
#else
#if VCEG_AZ06_IC
    if ( isAvailableB2 
      && ( !isAvailableB1 || !( pcCUAbove      ->hasEqualMotion(     uiAbovePartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) && pcCUAbove->getICFlag( uiAbovePartIdx ) == pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx ) ) )
      )
#else
    if ( isAvailableB2 && ( !isAvailableA1 || !pcCULeft->hasEqualMotion( uiLeftPartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) )
        && ( !isAvailableB1 || !pcCUAbove->hasEqualMotion( uiAbovePartIdx, pcCUAboveLeft, uiAboveLeftPartIdx ) ) )
#endif
#endif
    {
      abCandIsInter[iCount] = true;
#if VCEG_AZ06_IC
      pbICFlag[iCount] = pcCUAboveLeft->getICFlag( uiAboveLeftPartIdx );
#endif
      // get Inter Dir
      puhInterDirNeighbours[iCount] = pcCUAboveLeft->getInterDir( uiAboveLeftPartIdx );
      // get Mv from Left
      pcCUAboveLeft->getMvField( pcCUAboveLeft, uiAboveLeftPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbours[iCount<<1] );
      if ( getSlice()->isInterB() )
      {
        pcCUAboveLeft->getMvField( pcCUAboveLeft, uiAboveLeftPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbours[(iCount<<1)+1] );
      }
      if ( mrgCandIdx == iCount )
      {
        return;
      }
      iCount ++;
    }
  }
  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }

  if ( getSlice()->getEnableTMVPFlag() )
  {
    //>> MTK colocated-RightBottom
    UInt uiPartIdxRB;

    deriveRightBottomIdx( uiPUIdx, uiPartIdxRB );

    UInt uiAbsPartIdxTmp = g_auiZscanToRaster[uiPartIdxRB];
    const UInt numPartInCtuWidth  = m_pcPic->getNumPartInCtuWidth();
    const UInt numPartInCtuHeight = m_pcPic->getNumPartInCtuHeight();

    TComMv cColMv;
    Int iRefIdx;
    Int ctuRsAddr = -1;

    if (   ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + g_auiRasterToPelX[uiAbsPartIdxTmp] + m_pcPic->getMinCUWidth () ) < m_pcSlice->getSPS()->getPicWidthInLumaSamples () )  // image boundary check
        && ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + g_auiRasterToPelY[uiAbsPartIdxTmp] + m_pcPic->getMinCUHeight() ) < m_pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
    {
      if ( ( uiAbsPartIdxTmp % numPartInCtuWidth < numPartInCtuWidth - 1 ) &&           // is not at the last column of CTU
        ( uiAbsPartIdxTmp / numPartInCtuWidth < numPartInCtuHeight - 1 ) )              // is not at the last row    of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[ uiAbsPartIdxTmp + numPartInCtuWidth + 1 ];
        ctuRsAddr = getCtuRsAddr();
      }
      else if ( uiAbsPartIdxTmp % numPartInCtuWidth < numPartInCtuWidth - 1 )           // is not at the last column of CTU But is last row of CTU
      {
#if COM16_C806_GEN_MRG_IMPROVEMENT
        ctuRsAddr = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus();
#endif
        uiAbsPartAddr = g_auiRasterToZscan[ (uiAbsPartIdxTmp + numPartInCtuWidth + 1) % m_pcPic->getNumPartitionsInCtu() ];
      }
      else if ( uiAbsPartIdxTmp / numPartInCtuWidth < numPartInCtuHeight - 1 )          // is not at the last row of CTU But is last column of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[ uiAbsPartIdxTmp + 1 ];
        ctuRsAddr = getCtuRsAddr() + 1;
      }
      else //is the right bottom corner of CTU
      {
        uiAbsPartAddr = 0;
#if COM16_C806_GEN_MRG_IMPROVEMENT
        ctuRsAddr = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus() + 1;
#endif
      }
    }

    iRefIdx = 0;

    Bool bExistMV = false;
    UInt uiPartIdxCenter;
    Int dir = 0;
    UInt uiArrayAddr = iCount;
    xDeriveCenterIdx( uiPUIdx, uiPartIdxCenter );
#if VCEG_AZ06_IC
    Bool abTMVPICFlag[2] = { false, false }, bTMVPICFlag = false;
    bExistMV = ctuRsAddr >= 0 && xGetColMVP( REF_PIC_LIST_0, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx, &abTMVPICFlag[0] );
#else
    bExistMV = ctuRsAddr >= 0 && xGetColMVP( REF_PIC_LIST_0, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx );
#endif
    if( bExistMV == false )
    {
#if VCEG_AZ06_IC
      bExistMV = xGetColMVP( REF_PIC_LIST_0, getCtuRsAddr(), uiPartIdxCenter, cColMv, iRefIdx, &abTMVPICFlag[0] );
#else
      bExistMV = xGetColMVP( REF_PIC_LIST_0, getCtuRsAddr(), uiPartIdxCenter,  cColMv, iRefIdx );
#endif
    }
    if( bExistMV )
    {
      dir |= 1;
#if VCEG_AZ06_IC
      bTMVPICFlag |= abTMVPICFlag[0];
#endif
      pcMvFieldNeighbours[ 2 * uiArrayAddr ].setMvField( cColMv, iRefIdx );
    }

    if ( getSlice()->isInterB() )
    {
#if VCEG_AZ06_IC
      bExistMV = ctuRsAddr >= 0 && xGetColMVP( REF_PIC_LIST_1, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx, &abTMVPICFlag[1] );
#else
      bExistMV = ctuRsAddr >= 0 && xGetColMVP( REF_PIC_LIST_1, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx);
#endif
      if( bExistMV == false )
      {
#if VCEG_AZ06_IC
        bExistMV = xGetColMVP( REF_PIC_LIST_1, getCtuRsAddr(), uiPartIdxCenter, cColMv, iRefIdx, &abTMVPICFlag[1] );
#else
        bExistMV = xGetColMVP( REF_PIC_LIST_1, getCtuRsAddr(), uiPartIdxCenter, cColMv, iRefIdx );
#endif
      }
      if( bExistMV )
      {
        dir |= 2;
#if VCEG_AZ06_IC
        bTMVPICFlag |= abTMVPICFlag[1];
#endif
        pcMvFieldNeighbours[ 2 * uiArrayAddr + 1 ].setMvField( cColMv, iRefIdx );
      }
    }

    if (dir != 0)
    {
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP 
      Bool bAddT  = true;
      if (bAtmvpAva)
      {
        bAddT = false;
#if VCEG_AZ06_IC
        if ( dir!= puhInterDirNeighbours[uiAtmvpPos ] || bTMVPICFlag != pbICFlag[uiAtmvpPos] )
#else
        if ( dir!= puhInterDirNeighbours[uiAtmvpPos ] )
#endif
        {
          bAddT = true ;
        }
        else
        {
          for ( UInt uiR = 0; uiR< 2; uiR++ )
          {
            if ( dir & ( 1 << uiR ) )
            {
              if ( pcMvFieldNeighbours[ 2 * uiArrayAddr+uiR].getRefIdx()!= pcMvFieldNeighbours[2*uiAtmvpPos+uiR].getRefIdx())
              {
                bAddT = true;
                break;
              }
              else 
              {
                if (pcMvFieldNeighbours[ 2 * uiArrayAddr+uiR].getMv()  != pcMvFieldNeighbours[2*uiAtmvpPos+uiR].getMv())
                {
                  bAddT = true;
                  break;
                }
              }
            }
          }
        }
      }
#if COM16_C806_GEN_MRG_IMPROVEMENT
      Int iSpatCan = bAtmvpAva? (iCount-1): iCount;
      for(Int iCanIdx = 0; iCanIdx <iSpatCan; iCanIdx++)
      {
        if(puhInterDirNeighbours[iCanIdx]== dir && pcMvFieldNeighbours[iCanIdx<<1]==pcMvFieldNeighbours[uiArrayAddr<<1] && pcMvFieldNeighbours[1+ (iCanIdx<<1)]==pcMvFieldNeighbours[1+(uiArrayAddr<<1)]
#if VCEG_AZ06_IC
        && pbICFlag[iCanIdx] == bTMVPICFlag
#endif
        )
        {
          bAddT = false;
        }
      }
#endif
      if ( bAddT )
      {
#endif
      puhInterDirNeighbours[uiArrayAddr] = dir;
      abCandIsInter[uiArrayAddr] = true;
#if VCEG_AZ06_IC
      pbICFlag[uiArrayAddr] = bTMVPICFlag;
#endif
      if ( mrgCandIdx == iCount )
      {
        return;
      }
      iCount++;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
    }
#endif
    }
  }
  // early termination
  if (iCount == getSlice()->getMaxNumMergeCand())
  {
    return;
  }

  UInt uiArrayAddr = iCount;
  UInt uiCutoff = uiArrayAddr;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  uiCutoff      = ( uiCutoff > 4 ? 4: uiCutoff );
#endif
  if ( getSlice()->isInterB() )
  {
    static const UInt NUM_PRIORITY_LIST=12;
    static const UInt uiPriorityList0[NUM_PRIORITY_LIST] = {0 , 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3};
    static const UInt uiPriorityList1[NUM_PRIORITY_LIST] = {1 , 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2};

    for (Int idx=0; idx<uiCutoff*(uiCutoff-1) && uiArrayAddr!= getSlice()->getMaxNumMergeCand(); idx++)
    {
      assert(idx<NUM_PRIORITY_LIST);
      Int i = uiPriorityList0[idx];
      Int j = uiPriorityList1[idx];
      if (abCandIsInter[i] && abCandIsInter[j]&& (puhInterDirNeighbours[i]&0x1)&&(puhInterDirNeighbours[j]&0x2))
      {
        abCandIsInter[uiArrayAddr] = true;
        puhInterDirNeighbours[uiArrayAddr] = 3;
#if VCEG_AZ06_IC
        pbICFlag[uiArrayAddr] = pbICFlag[i] || pbICFlag[j];
#endif
        // get Mv from cand[i] and cand[j]
        pcMvFieldNeighbours[uiArrayAddr << 1].setMvField(pcMvFieldNeighbours[i<<1].getMv(), pcMvFieldNeighbours[i<<1].getRefIdx());
        pcMvFieldNeighbours[( uiArrayAddr << 1 ) + 1].setMvField(pcMvFieldNeighbours[(j<<1)+1].getMv(), pcMvFieldNeighbours[(j<<1)+1].getRefIdx());

        Int iRefPOCL0 = m_pcSlice->getRefPOC( REF_PIC_LIST_0, pcMvFieldNeighbours[(uiArrayAddr<<1)].getRefIdx() );
        Int iRefPOCL1 = m_pcSlice->getRefPOC( REF_PIC_LIST_1, pcMvFieldNeighbours[(uiArrayAddr<<1)+1].getRefIdx() );
        if (iRefPOCL0 == iRefPOCL1 && pcMvFieldNeighbours[(uiArrayAddr<<1)].getMv() == pcMvFieldNeighbours[(uiArrayAddr<<1)+1].getMv())
        {
          abCandIsInter[uiArrayAddr] = false;
        }
        else
        {
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
          if(uiArrayAddr == mrgCandIdx)
          {
            return;
          }
#endif
          uiArrayAddr++;
        }
      }
    }
  }
  // early termination
  if (uiArrayAddr == getSlice()->getMaxNumMergeCand())
  {
    return;
  }

  Int iNumRefIdx = (getSlice()->isInterB()) ? min(m_pcSlice->getNumRefIdx(REF_PIC_LIST_0), m_pcSlice->getNumRefIdx(REF_PIC_LIST_1)) : m_pcSlice->getNumRefIdx(REF_PIC_LIST_0);

  Int r = 0;
  Int refcnt = 0;
  while (uiArrayAddr < getSlice()->getMaxNumMergeCand())
  {
    abCandIsInter[uiArrayAddr] = true;
    puhInterDirNeighbours[uiArrayAddr] = 1;
    pcMvFieldNeighbours[uiArrayAddr << 1].setMvField( TComMv(0, 0), r);
#if VCEG_AZ06_IC
    pbICFlag[uiArrayAddr] = false;
#endif
    if ( getSlice()->isInterB() )
    {
      puhInterDirNeighbours[uiArrayAddr] = 3;
      pcMvFieldNeighbours[(uiArrayAddr << 1) + 1].setMvField(TComMv(0, 0), r);
    }
    uiArrayAddr++;

    if ( refcnt == iNumRefIdx - 1 )
    {
      r = 0;
    }
    else
    {
      ++r;
      ++refcnt;
    }
  }
  numValidMergeCand = uiArrayAddr;
}
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
/** Constructs the ATMVP candidate with a given temporal vector
 * \param uiPUIdx
 * \param pcMvFieldSP 
 * \param puhInterDirSP
 * \param pcMvFieldDefault
 * \param pcInterDirDefault
 * \param numValidMergeCand
 * \param cTMv
 * \param bMrgIdxMatchATMVPCan
 * \param iPocColPic
 * \param pDecCurrCU
 * \param uiDecCurrAbsPartIdx
  */
#if JVET_C0035_ATMVP_SIMPLIFICATION
Bool TComDataCU::getInterMergeSubPUTmvpCandidate ( UInt uiAbsPartIdx, UInt uiPUIdx,UInt uiCount,TComMvField* pcMvFieldNeighbours, UChar* puhInterDirNeighbours,
                                                  TComMvField* pcMvFieldSP, UChar* puhInterDirSP, TComMvField* pcMvFieldDefault, UChar* pcInterDirDefault, Bool bMrgIdxMatchATMVPCan, 
#if VCEG_AZ06_IC
                                                  Bool& rbICFlag,
#endif
                                                  TComDataCU* pDecCurrCU , UInt uiDecCurrAbsPartIdx)
#else
Bool TComDataCU::getInterMergeSubPUTmvpCandidate ( UInt uiPUIdx,  TComMvField* pcMvFieldSP, UChar* puhInterDirSP,
                                                   TComMvField* pcMvFieldDefault, UChar* pcInterDirDefault, TComMv cTMv, Bool bMrgIdxMatchATMVPCan, 
#if VCEG_AZ06_IC
                                                   Bool& rbICFlag,
#endif
                                                   Int iPocColPic, TComDataCU* pDecCurrCU, UInt uiDecCurrAbsPartIdx)
#endif
{
#if JVET_C0035_ATMVP_SIMPLIFICATION
  TComPic *pColPic = getSlice()->getRefPic( RefPicList(getSlice()->isInterB() ? 1-getSlice()->getColFromL0Flag() : 0), getSlice()->getColRefIdx());
  Int iPocColPic             = pColPic->getPOC();
  TComMv cTMv;

  RefPicList eFetchRefPicList=RefPicList(getSlice()->isInterB() ? 1-getSlice()->getColFromL0Flag() : 0);
  if (uiCount)
  {
    UInt uiN=0;
    for( UInt uiCurrRefListId = 0; uiCurrRefListId < (getSlice()->getSliceType() == B_SLICE ?  2 : 1 ) ; uiCurrRefListId++ )
    {
      RefPicList  eCurrRefPicList = RefPicList( RefPicList( getSlice()->isInterB() ? (getSlice()->getColFromL0Flag()? uiCurrRefListId: 1- uiCurrRefListId) : uiCurrRefListId ));
      if ( puhInterDirNeighbours[uiN] & (1<<eCurrRefPicList) )
      {
        iPocColPic = getSlice()->getRefPic( eCurrRefPicList, pcMvFieldNeighbours[uiN*2+eCurrRefPicList].getRefIdx())->getPOC() ;
        cTMv  = pcMvFieldNeighbours[uiN*2+eCurrRefPicList].getMv();
        eFetchRefPicList = eCurrRefPicList;
        break;
      }
    }
  }
#endif

  TComDataCU* pcTempCU;
  Int iTempCUAddr, iTempAbsPartIdx, iTempPosX, iTempPosY;
  Int iPartition     = 0;
  Int iInterDirSaved = 0;
#if VCEG_AZ06_IC
  Bool bTmpICFlag;
#endif
  ///////////////////////////////////////////////////////////////////////
  ////////          GET Initial Temporal Vector                  ////////
  ///////////////////////////////////////////////////////////////////////
  
  TComMv cTempVector = cTMv;

  // compute the location of the current PU
  Int iCurrPosX, iCurrPosY, iWidth, iHeight;
  Int iPUWidth, iPUHeight, iNumPart, iNumPartLine;

  getPartPosition( uiPUIdx, iCurrPosX, iCurrPosY,      iWidth, iHeight);
  getSPPara      (  iWidth,  iHeight,  iNumPart, iNumPartLine, iPUWidth, iPUHeight);

  Int iOffsetX = iPUWidth/2;;
  Int iOffsetY = iPUHeight/2;
  
  TComMv cColMv;
  // use coldir.
  Bool bBSlice = getSlice()->isInterB();
  UInt bColL0  = getSlice()->getColFromL0Flag();
#if JVET_C0035_ATMVP_SIMPLIFICATION
  pColPic = getPicfromPOC(iPocColPic);
#else
  TComPic *pColPic = getPicfromPOC(iPocColPic);
#endif


  Bool bATMVPAvailFlag = false;
  TComMvField cDefaultMvField[2];
  cDefaultMvField[0].getMv().set(0, 0);
  cDefaultMvField[1].getMv().set(0, 0);


  Int         iTempCenterCUAddr, iTempCenterAbsPartIdx;
  Int         iCenterPosX, iCenterPosY;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  Int nOffset = 1 << ( VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE + 1 ); // 2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE - 1
#endif
  
  Bool bInit = false;
  for( UInt uiLX = 0; uiLX < (bBSlice ? 2:1) && !bATMVPAvailFlag; uiLX++)
  {
    RefPicList eListY = RefPicList( bBSlice ? (bColL0 ? uiLX: 1- uiLX) : uiLX );
    for (Int refIdxY = (bInit ? 0 : -1); refIdxY < getSlice()->getNumRefIdx(eListY) && !bATMVPAvailFlag; refIdxY++)
    {
      if (!bInit) bInit = true;
      else
      {
        pColPic = getSlice()->getRefPic(eListY, refIdxY);
#if JVET_C0035_ATMVP_SIMPLIFICATION
        eFetchRefPicList = eListY;
#endif
      }
      Int iNewColPicPOC = pColPic->getPOC();
      if ( iNewColPicPOC!= iPocColPic)
      {
      //////////////// POC based scaling of the temporal vector /////////////
          Int iScale = xGetDistScaleFactor(getSlice()->getPOC(), iNewColPicPOC, getSlice()->getPOC(), iPocColPic);
        if ( iScale != 4096 )
          cTempVector=cTMv.scaleMv( iScale );
      }
      else
        cTempVector=cTMv;

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      iCenterPosX = iCurrPosX + ( ( iWidth /  iPUWidth ) >> 1 )  * iPUWidth + ( iPUWidth >> 1 ) +   ((cTempVector.getHor()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE)); 
      iCenterPosY = iCurrPosY + ( ( iHeight /  iPUHeight ) >> 1 )  * iPUHeight + (iPUHeight >> 1) + ((cTempVector.getVer()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
      iCenterPosX = iCurrPosX + ( ( iWidth /  iPUWidth ) >> 1 )  * iPUWidth + ( iPUWidth >> 1 ) +   ((cTempVector.getHor()+2)>>2); 
      iCenterPosY = iCurrPosY + ( ( iHeight /  iPUHeight ) >> 1 )  * iPUHeight + (iPUHeight >> 1) + ((cTempVector.getVer()+2)>>2) ;
#endif


      if(iWidth == iPUWidth && iHeight == iPUHeight)
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        iCenterPosX = iCurrPosX + (iWidth >> 1)  + ((cTempVector.getHor()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
        iCenterPosY = iCurrPosY + (iHeight >> 1) + ((cTempVector.getVer()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
        iCenterPosX = iCurrPosX + (iWidth >> 1)  + ((cTempVector.getHor()+2)>>2);
        iCenterPosY = iCurrPosY + (iHeight >> 1) + ((cTempVector.getVer()+2)>>2);
#endif
      }

      iCenterPosX = Clip3( 0, pColPic->getPicYuvRec()->getWidth (COMPONENT_Y) - 1,  iCenterPosX  );
      iCenterPosY = Clip3( 0, pColPic->getPicYuvRec()->getHeight(COMPONENT_Y) - 1,  iCenterPosY  );

      // derivation of center motion parameters from the collocated CU
      pColPic->getCUAddrAndPartIdx( iCenterPosX , iCenterPosY , iTempCenterCUAddr, iTempCenterAbsPartIdx );
      TComDataCU* pcDefaultCU    = pColPic->getCtu( iTempCenterCUAddr );
      if( pcDefaultCU->getPredictionMode( iTempCenterAbsPartIdx ) != MODE_INTRA )
      {
        for( UInt uiCurrRefListId = 0; uiCurrRefListId < (bBSlice ? 2:1) ; uiCurrRefListId++ )
        {
          RefPicList  eCurrRefPicList = RefPicList( uiCurrRefListId );
          if (deriveScaledMotionTemporalForOneDirection(pcDefaultCU,eCurrRefPicList,cColMv, iTempCenterAbsPartIdx,0 
#if VCEG_AZ06_IC
            , bTmpICFlag
#endif
            ,pColPic
#if JVET_C0035_ATMVP_SIMPLIFICATION
            ,eFetchRefPicList
#endif
            ))
          {
            cDefaultMvField[uiCurrRefListId].setMvField(cColMv,0);
            bATMVPAvailFlag = true; // keep this variable here for later algrithm tuning
#if VCEG_AZ06_IC
            rbICFlag = bTmpICFlag;
#endif
          }
        }
      }
    }
  }
  // The advanced TMVP candidate is considered as available if and only if the center block contains motion 
  if ( bATMVPAvailFlag == true )
  {   
    // perform ATMVP based on center now
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    iCurrPosX += ((cTempVector.getHor()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
    iCurrPosY += ((cTempVector.getVer()+nOffset)>>(2+VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
#else
    iCurrPosX += ((cTempVector.getHor()+2)>>2);
    iCurrPosY += ((cTempVector.getVer()+2)>>2);
#endif
    iInterDirSaved = (cDefaultMvField[0].getRefIdx() !=-1 ? 1: 0) + (cDefaultMvField[1].getRefIdx() !=-1 ? 2: 0);

    Int iPicWidth  = pColPic->getPicYuvRec()->getWidth (COMPONENT_Y) - 1;
    Int iPicHeight = pColPic->getPicYuvRec()->getHeight(COMPONENT_Y) - 1;

    for (Int i=iCurrPosY; i < iCurrPosY + iHeight; i += iPUHeight)
      for (Int j = iCurrPosX; j < iCurrPosX + iWidth; j += iPUWidth)
      {
        iTempPosX     = j + iOffsetX;
        iTempPosY     = i + iOffsetY; 

        iTempPosX = Clip3( 0, iPicWidth,   iTempPosX  );
        iTempPosY = Clip3( 0, iPicHeight,  iTempPosY  );

        pColPic->getCUAddrAndPartIdx( iTempPosX, iTempPosY, iTempCUAddr, iTempAbsPartIdx );
        pcTempCU  = pColPic->getCtu( iTempCUAddr );

        if( pcTempCU && !pcTempCU->isIntra(iTempAbsPartIdx) )
        {
          for( UInt uiCurrRefListId = 0; uiCurrRefListId < (bBSlice ? 2:1); uiCurrRefListId++ )
          {
            RefPicList  eCurrRefPicList = RefPicList( uiCurrRefListId );
            if (deriveScaledMotionTemporalForOneDirection(pcTempCU,eCurrRefPicList, cColMv, iTempAbsPartIdx, 0, 
#if VCEG_AZ06_IC
              bTmpICFlag,
#endif
              pColPic 
#if JVET_C0035_ATMVP_SIMPLIFICATION
              ,eFetchRefPicList
#endif
              ))
            {
              pcMvFieldSP[2*iPartition + uiCurrRefListId].setMvField(cColMv,0);
            }                
          }
        }
        else // intra coded, in this case, no motion vector is available for list 0 or list 1
        {
          pcMvFieldSP[2*iPartition + 0].setMvField(cDefaultMvField[0].getMv(),cDefaultMvField[0].getRefIdx());
          pcMvFieldSP[2*iPartition + 1].setMvField(cDefaultMvField[1].getMv(),cDefaultMvField[1].getRefIdx());
        }
        puhInterDirSP[iPartition] = (pcMvFieldSP[2*iPartition].getRefIdx()!=-1 ? 1: 0) + (pcMvFieldSP[2*iPartition+1].getRefIdx()!=-1 ? 2: 0);
        iPartition++;
      }  
      pcMvFieldDefault[0].setMvField(cDefaultMvField[0].getMv(),cDefaultMvField[0].getRefIdx());
      pcMvFieldDefault[1].setMvField(cDefaultMvField[1].getMv(),cDefaultMvField[1].getRefIdx());
      pcInterDirDefault[0] = iInterDirSaved ;
       
      if( bMrgIdxMatchATMVPCan) // only invoked at the decoder since here the parsing of the motion vectors are skipped
      {
        //store the sub-PU motion information
        UInt uiSPAddr;
        Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
        getSPPara( iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);
        for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
        {
          pDecCurrCU->getSPAbsPartIdx(uiDecCurrAbsPartIdx, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
          pDecCurrCU->setInterDirSP(puhInterDirSP[iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
          pDecCurrCU->getCUMvField( REF_PIC_LIST_0 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[2*iPartitionIdx], iSPWidth, iSPHeight);
          pDecCurrCU->getCUMvField( REF_PIC_LIST_1 )->setMvFieldSP(pDecCurrCU, uiSPAddr, pcMvFieldSP[2*iPartitionIdx + 1], iSPWidth, iSPHeight);
        }
      }
      return true;
  } 
  return false;
}
Void TComDataCU::getSPPara(Int iPUWidth, Int iPUHeight, Int& iNumSP, Int& iNumSPInOneLine, Int& iSPWidth, Int& iSPHeight)
{
  Int iSubPUSize = 1 << getSlice()->getSPS()->getSubPUTLog2Size();

  iNumSPInOneLine = iPUWidth/iSubPUSize;
  iNumSPInOneLine = iNumSPInOneLine < 1 ? 1: iNumSPInOneLine;
  Int iNumSPInOneColumn = iPUHeight/iSubPUSize;
  iNumSPInOneColumn = iNumSPInOneColumn < 1 ? 1: iNumSPInOneColumn;
  iNumSP = iNumSPInOneLine * iNumSPInOneColumn;

  iSPWidth = iNumSPInOneLine == 1 ? iPUWidth: iSubPUSize; 
  iSPHeight = iNumSPInOneColumn == 1 ? iPUHeight: iSubPUSize; 
}

Void TComDataCU::getSPAbsPartIdx(UInt uiBaseAbsPartIdx, Int iWidth, Int iHeight, Int iPartIdx, Int iNumPartLine, UInt& ruiPartAddr )
{
#if TT_TREE_ATMVP // deriveZIdxInCTU
  uiBaseAbsPartIdx  = deriveZIdxInCTU(m_absZIdxInCtu, uiBaseAbsPartIdx);
#else
  uiBaseAbsPartIdx += m_absZIdxInCtu;
#endif
  Int iBasePelX = g_auiRasterToPelX[g_auiZscanToRaster[uiBaseAbsPartIdx]];
  Int iBasePelY = g_auiRasterToPelY[g_auiZscanToRaster[uiBaseAbsPartIdx]];
  Int iCurrPelX = iBasePelX + iPartIdx%iNumPartLine * iWidth;
  Int iCurrPelY = iBasePelY + iPartIdx/iNumPartLine * iHeight;
  Int iCurrRaster = iCurrPelY / getPic()->getMinCUHeight() * getPic()->getNumPartInCtuWidth() + iCurrPelX/getPic()->getMinCUWidth();
  ruiPartAddr = g_auiRasterToZscan[iCurrRaster];
#if TT_TREE_ATMVP // deriveZIdxInCTU
  ruiPartAddr = deriveZIdxInCU(m_absZIdxInCtu, ruiPartAddr); 
#else
  ruiPartAddr -= m_absZIdxInCtu;  
#endif
}

Void TComDataCU::setInterDirSP( UInt uiDir, UInt uiAbsPartIdx, Int iWidth, Int iHeight )
{
#if TT_TREE_ATMVP // deriveZIdxInCTU
  uiAbsPartIdx  = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  uiAbsPartIdx += m_absZIdxInCtu;
#endif
  Int iStartPelX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  Int iStartPelY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
  Int iEndPelX = iStartPelX + iWidth;
  Int iEndPelY = iStartPelY + iHeight;

  Int iCurrRaster, uiPartAddr;

  for (Int i=iStartPelY; i<iEndPelY; i+=getPic()->getMinCUHeight())
  {
    for (Int j=iStartPelX; j < iEndPelX; j += getPic()->getMinCUWidth())
    {
      iCurrRaster = i / getPic()->getMinCUHeight() * getPic()->getNumPartInCtuWidth() + j/getPic()->getMinCUWidth();
      uiPartAddr = g_auiRasterToZscan[iCurrRaster];
#if TT_TREE_ATMVP // deriveZIdxInCU
      uiPartAddr  = deriveZIdxInCU(m_absZIdxInCtu, uiPartAddr);
#else
      uiPartAddr -= m_absZIdxInCtu;  
#endif

      m_puhInterDir[uiPartAddr] = uiDir;
    }
  }
}
Bool TComDataCU::deriveScaledMotionTemporalForOneDirection( TComDataCU* pcTempCU,RefPicList eCurrRefPicList, TComMv &cColMv, UInt uiAbsPartIdx, Int iTargetRefIdx, 
#if VCEG_AZ06_IC
  Bool& rbICFlag,
#endif
  TComPic *pColPic
#if JVET_C0035_ATMVP_SIMPLIFICATION
  ,  RefPicList eFetchRefPicList
#endif 
  )
{
  Int iColPOC, iColRefPOC, iCurrPOC, iCurrRefPOC, iScale;
#if JVET_C0035_ATMVP_SIMPLIFICATION
  Bool bAllowMirrorMV= true; 
  RefPicList eColRefPicList=  getSlice()->getCheckLDC()?eCurrRefPicList : RefPicList(1-eFetchRefPicList);
  if ( pColPic == m_pcSlice->getRefPic(RefPicList(getSlice()->isInterB() ? 1-getSlice()->getColFromL0Flag() : 0 ),getSlice()->getColRefIdx()))
  {
    eColRefPicList = eCurrRefPicList;   //67 -> disable, 64 -> enable
    bAllowMirrorMV =false;
  }
#else
  Bool bAllowMirrorMV= false; 
  if ( pColPic == m_pcSlice->getRefPic(RefPicList(getSlice()->getColFromL0Flag()),getSlice()->getColRefIdx()))
  {
    bAllowMirrorMV =true;
  }
  else
  {
    if (getSlice()->getCheckLDC())
    {
      if (getWidth(0)  <= 32 && getHeight(0)  <= 32)
      {
        bAllowMirrorMV = true;
      }
    }
    else
    {
      if (getWidth(0)  <= 8 && getHeight(0)  <= 8)
      {
        bAllowMirrorMV = true;
      }
    }
  }
  RefPicList eColRefPicList = (getSlice()->getCheckLDC() || !bAllowMirrorMV) ? eCurrRefPicList :  RefPicList(getSlice()->getColFromL0Flag());
#endif

#if VCEG_AZ07_FRUC_MERGE
  if( pcTempCU->getSlice()->getSPS()->getUseFRUCMgrMode() )
  {
    eColRefPicList = eCurrRefPicList;
  }
#endif
  // Although it might make sense to keep the unavailable motion field per direction still be unavailable, I made the MV prediction the same way as in TMVP
  // So there is an interaction between MV0 and MV1 of the corresponding blocks identified by TV.

  // Grab motion and do necessary scaling.{{
  iCurrPOC = m_pcSlice->getPOC();
  Int iColRefIdx = pcTempCU->getCUMvField(RefPicList(eColRefPicList))->getRefIdx(uiAbsPartIdx);
#if JVET_C0035_ATMVP_SIMPLIFICATION
  if (iColRefIdx < 0 && (getSlice()->getCheckLDC() || bAllowMirrorMV)
#else
  if (iColRefIdx < 0 && bAllowMirrorMV
#endif
#if VCEG_AZ07_FRUC_MERGE
    && !pcTempCU->getSlice()->getSPS()->getUseFRUCMgrMode()
#endif
    )
  {
    eColRefPicList = RefPicList(1 - eColRefPicList);
    iColRefIdx = pcTempCU->getCUMvField(RefPicList(eColRefPicList))->getRefIdx(uiAbsPartIdx);

    if (iColRefIdx < 0 )
    {
      return false;
    }
  }
  if ( iColRefIdx >=0 && iTargetRefIdx>=0 && iTargetRefIdx<m_pcSlice->getNumRefIdx(eCurrRefPicList)) 
  {
    iColPOC     = pcTempCU->getSlice()->getPOC();
    iColRefPOC  = pcTempCU->getSlice()->getRefPOC(eColRefPicList, iColRefIdx);
    //////////////////////////////////////////////////////////////
    // Set the target reference index to 0, may be changed later;//
    //////////////////////////////////////////////////////////////
    iCurrRefPOC = m_pcSlice->getRefPic(eCurrRefPicList, iTargetRefIdx)->getPOC();
    // Scale the vector.
    cColMv      = pcTempCU->getCUMvField(eColRefPicList)->getMv(uiAbsPartIdx);
    //pcMvFieldSP[2*iPartition + eCurrRefPicList].getMv();
    // Assume always short-term for now
    iScale = xGetDistScaleFactor(iCurrPOC, iCurrRefPOC, iColPOC, iColRefPOC);
    if ( iScale != 4096 )
      cColMv=cColMv.scaleMv( iScale );
#if VCEG_AZ06_IC
    rbICFlag = pcTempCU->getICFlag( uiAbsPartIdx );
#endif
    return true;
  }
  return false;
}
TComPic * TComDataCU::getPicfromPOC(Int iPocColPic)
{
  TComPic * pColPic = NULL;
  for( UInt uiCurrRefListId = 0; uiCurrRefListId < (getSlice()->getSliceType() == B_SLICE ?  2 : 1 ) ; uiCurrRefListId++ )
  {
    RefPicList  eRefPicListX = RefPicList( uiCurrRefListId );
    for( UInt uiRefIdx =0; uiRefIdx<getSlice()->getNumRefIdx(eRefPicListX); uiRefIdx++)
    {
      if( getSlice()->getRefPic(eRefPicListX,uiRefIdx)->getPOC() == iPocColPic )
      {
        pColPic   = getSlice()->getRefPic(eRefPicListX,uiRefIdx);
        return pColPic;
      }
    }
  }
  return pColPic;
}
#endif 
/** Check whether the current PU and a spatial neighboring PU are in a same ME region.
 * \param xN, yN   location of the upper-left corner pixel of a neighboring PU
 * \param xP, yP   location of the upper-left corner pixel of the current PU
 */
Bool TComDataCU::isDiffMER(Int xN, Int yN, Int xP, Int yP)
{

  UInt plevel = this->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() + 2;
  if ((xN>>plevel)!= (xP>>plevel))
  {
    return true;
  }
  if ((yN>>plevel)!= (yP>>plevel))
  {
    return true;
  }
  return false;
}

/** Calculate the location of upper-left corner pixel and size of the current PU.
 * \param partIdx       PU index within a CU
 * \param xP, yP        location of the upper-left corner pixel of the current PU
 * \param nPSW, nPSH    size of the current PU
 */
Void TComDataCU::getPartPosition( UInt partIdx, Int& xP, Int& yP, Int& nPSW, Int& nPSH)
{
  UInt col = m_uiCUPelX;
  UInt row = m_uiCUPelY;

#if JVET_C0024_QTBT
  nPSW = getWidth(0);
  nPSH = getHeight(0);
  xP   = col ;
  yP   = row ;
#else
  switch ( m_pePartSize[0] )
  {
  case SIZE_2NxN:
    nPSW = getWidth(0);
    nPSH = getHeight(0) >> 1;
    xP   = col;
    yP   = (partIdx ==0)? row: row + nPSH;
    break;
  case SIZE_Nx2N:
    nPSW = getWidth(0) >> 1;
    nPSH = getHeight(0);
    xP   = (partIdx ==0)? col: col + nPSW;
    yP   = row;
    break;
  case SIZE_NxN:
    nPSW = getWidth(0) >> 1;
    nPSH = getHeight(0) >> 1;
    xP   = col + (partIdx&0x1)*nPSW;
    yP   = row + (partIdx>>1)*nPSH;
    break;
  case SIZE_2NxnU:
    nPSW = getWidth(0);
    nPSH = ( partIdx == 0 ) ?  getHeight(0) >> 2 : ( getHeight(0) >> 2 ) + ( getHeight(0) >> 1 );
    xP   = col;
    yP   = (partIdx ==0)? row: row + getHeight(0) - nPSH;

    break;
  case SIZE_2NxnD:
    nPSW = getWidth(0);
    nPSH = ( partIdx == 0 ) ?  ( getHeight(0) >> 2 ) + ( getHeight(0) >> 1 ) : getHeight(0) >> 2;
    xP   = col;
    yP   = (partIdx ==0)? row: row + getHeight(0) - nPSH;
    break;
  case SIZE_nLx2N:
    nPSW = ( partIdx == 0 ) ? getWidth(0) >> 2 : ( getWidth(0) >> 2 ) + ( getWidth(0) >> 1 );
    nPSH = getHeight(0);
    xP   = (partIdx ==0)? col: col + getWidth(0) - nPSW;
    yP   = row;
    break;
  case SIZE_nRx2N:
    nPSW = ( partIdx == 0 ) ? ( getWidth(0) >> 2 ) + ( getWidth(0) >> 1 ) : getWidth(0) >> 2;
    nPSH = getHeight(0);
    xP   = (partIdx ==0)? col: col + getWidth(0) - nPSW;
    yP   = row;
    break;
  default:
    assert ( m_pePartSize[0] == SIZE_2Nx2N );
    nPSW = getWidth(0);
    nPSH = getHeight(0);
    xP   = col ;
    yP   = row ;

    break;
  }
#endif
}

/** Constructs a list of candidates for AMVP
 * \param uiPartIdx
 * \param uiPartAddr
 * \param eRefPicList
 * \param iRefIdx
 * \param pInfo
 */
Void TComDataCU::fillMvpCand ( UInt uiPartIdx, UInt uiPartAddr, RefPicList eRefPicList, Int iRefIdx, AMVPInfo* pInfo 
#if VCEG_AZ07_FRUC_MERGE
  , TComPrediction * pPred
#endif
  )
{
  TComMv cMvPred;
  Bool bAddedSmvp = false;

  pInfo->iN = 0;
  if (iRefIdx < 0)
  {
    return;
  }

  //-- Get Spatial MV
  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
  const UInt numPartInCtuWidth  = m_pcPic->getNumPartInCtuWidth();
  const UInt numPartInCtuHeight = m_pcPic->getNumPartInCtuHeight();
  Bool bAdded = false;

  deriveLeftRightTopIdx( uiPartIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdx( uiPartIdx, uiPartIdxLB );

  TComDataCU* tmpCU = NULL;
  UInt idx;
  tmpCU = getPUBelowLeft(idx, uiPartIdxLB);
  bAddedSmvp = (tmpCU != NULL) && (tmpCU->isInter(idx));

  if (!bAddedSmvp)
  {
    tmpCU = getPULeft(idx, uiPartIdxLB);
    bAddedSmvp = (tmpCU != NULL) && (tmpCU->isInter(idx));
  }

  // Left predictor search
  bAdded = xAddMVPCand( pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT);
  if (!bAdded)
  {
    bAdded = xAddMVPCand( pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT );
  }

  if(!bAdded)
  {
    bAdded = xAddMVPCandOrder( pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT);
    if (!bAdded)
    {
      xAddMVPCandOrder( pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT );
    }
  }

  // Above predictor search
  bAdded = xAddMVPCand( pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT);

  if (!bAdded)
  {
    bAdded = xAddMVPCand( pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE);
  }

  if(!bAdded)
  {
    xAddMVPCand( pInfo, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT);
  }

  if(!bAddedSmvp)
  {
    bAdded = xAddMVPCandOrder( pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT);
    if (!bAdded)
    {
      bAdded = xAddMVPCandOrder( pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE);
    }

    if(!bAdded)
    {
      xAddMVPCandOrder( pInfo, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT);
    }
  }

#if VCEG_AZ07_IMV
  if( getiMVFlag( uiPartAddr ) && getSlice()->getSPS()->getIMV() )
  {
    for( Int i = 0; i < pInfo->iN; i++ )
    {
#if JVET_E0076_MULTI_PEL_MVD
      xRoundMV( pInfo->m_acMvCand[i], uiPartAddr );
#else
      xRoundMV( pInfo->m_acMvCand[i] );
#endif
    }
  }
#endif

  if ( pInfo->iN == 2 )
  {
    if ( pInfo->m_acMvCand[ 0 ] == pInfo->m_acMvCand[ 1 ] )
    {
      pInfo->iN = 1;
    }
  }

  if ( getSlice()->getEnableTMVPFlag() )
  {
    // Get Temporal Motion Predictor
    Int iRefIdx_Col = iRefIdx;
    TComMv cColMv;
    UInt uiPartIdxRB;
    UInt uiAbsPartIdx;
    UInt uiAbsPartAddr;

    deriveRightBottomIdx( uiPartIdx, uiPartIdxRB );
#if TT_TREE_DS5 // deriveZIdxInCTU
    uiAbsPartAddr = deriveZIdxInCTU(m_absZIdxInCtu, uiPartAddr);
#else
    uiAbsPartAddr = m_absZIdxInCtu + uiPartAddr;
#endif

    //----  co-located RightBottom Temporal Predictor (H) ---//
    uiAbsPartIdx = g_auiZscanToRaster[uiPartIdxRB];
    Int ctuRsAddr = -1;
    if (  ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + g_auiRasterToPelX[uiAbsPartIdx] + m_pcPic->getMinCUWidth () ) < m_pcSlice->getSPS()->getPicWidthInLumaSamples () )  // image boundary check
       && ( ( m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + g_auiRasterToPelY[uiAbsPartIdx] + m_pcPic->getMinCUHeight() ) < m_pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
    {
      if ( ( uiAbsPartIdx % numPartInCtuWidth < numPartInCtuWidth - 1 ) &&  // is not at the last column of CTU
           ( uiAbsPartIdx / numPartInCtuWidth < numPartInCtuHeight - 1 ) )  // is not at the last row    of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[ uiAbsPartIdx + numPartInCtuWidth + 1 ];
        ctuRsAddr = getCtuRsAddr();
      }
      else if ( uiAbsPartIdx % numPartInCtuWidth < numPartInCtuWidth - 1 )  // is not at the last column of CTU But is last row of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[ (uiAbsPartIdx + numPartInCtuWidth + 1) % m_pcPic->getNumPartitionsInCtu() ];
      }
      else if ( uiAbsPartIdx / numPartInCtuWidth < numPartInCtuHeight - 1 ) // is not at the last row of CTU But is last column of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[ uiAbsPartIdx + 1 ];
        ctuRsAddr = getCtuRsAddr() + 1;
      }
      else //is the right bottom corner of CTU
      {
        uiAbsPartAddr = 0;
      }
    }
    if ( ctuRsAddr >= 0 && xGetColMVP( eRefPicList, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx_Col ) )
    {
      pInfo->m_acMvCand[pInfo->iN++] = cColMv;
    }
    else
    {
      UInt uiPartIdxCenter;
      xDeriveCenterIdx( uiPartIdx, uiPartIdxCenter );
      if (xGetColMVP( eRefPicList, getCtuRsAddr(), uiPartIdxCenter,  cColMv, iRefIdx_Col ))
      {
        pInfo->m_acMvCand[pInfo->iN++] = cColMv;
      }
    }
    //----  co-located RightBottom Temporal Predictor  ---//
  }

#if VCEG_AZ07_FRUC_MERGE
  if( getSlice()->getSPS()->getUseFRUCMgrMode() )
  {
    if( pPred != NULL && pPred->deriveFRUCMV( this , getDepth( uiPartAddr ) , uiPartAddr , uiPartIdx , iRefIdx , eRefPicList
#if JVET_E0060_FRUC_CAND
        , pInfo
#endif
        ) )
    {
      const TComMv & rMV = getCUMvField( eRefPicList )->getMv( uiPartAddr );
      if( pInfo->iN == 0 )
      {
        pInfo->m_acMvCand[0] = rMV;
        pInfo->iN++;
      }
      else if( pInfo->m_acMvCand[0] != rMV )
      {
        for( Int n = min( pInfo->iN , AMVP_MAX_NUM_CANDS - 1 ) ; n > 0 ; n-- )
        {
          pInfo->m_acMvCand[n] = pInfo->m_acMvCand[n-1];
        }
        pInfo->m_acMvCand[0] = rMV;
        pInfo->iN = min( pInfo->iN + 1 , AMVP_MAX_NUM_CANDS );
      }
    }
  }
#endif

  if (pInfo->iN > AMVP_MAX_NUM_CANDS)
  {
    pInfo->iN = AMVP_MAX_NUM_CANDS;
  }

  while (pInfo->iN < AMVP_MAX_NUM_CANDS)
  {
    pInfo->m_acMvCand[pInfo->iN].set(0,0);
    pInfo->iN++;
  }
  
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  for( Int i = 0 ; i < pInfo->iN ; i++ )
  {
    pInfo->m_acMvCand[i].roundMV2SignalPrecision();
  }
#endif

#if VCEG_AZ07_IMV
  if( getiMVFlag( uiPartAddr ) && getSlice()->getSPS()->getIMV() )
  {
    for( Int i = 0 ; i < pInfo->iN ; i++ )
    {
#if JVET_E0076_MULTI_PEL_MVD
      xRoundMV( pInfo->m_acMvCand[i],uiPartAddr );
#else
      xRoundMV( pInfo->m_acMvCand[i] );
#endif
    }
  }
#endif

  return ;
}


Bool TComDataCU::isBipredRestriction(UInt puIdx)
{
#if JVET_C0024_QTBT
  return false;
#endif
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
  if (getSlice()->getSPS()->getAtmvpEnableFlag())
  {
    return false;
  }
  else
  {
#endif
  Int width = 0;
  Int height = 0;
  UInt partAddr;

  getPartIndexAndSize( puIdx, partAddr, width, height );
  if ( getWidth(0) == 8 && (width < 8 || height < 8) )
  {
    return true;
  }
  return false;
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
  }
#endif
}


Void TComDataCU::clipMv    (TComMv&  rcMv)
{
  const TComSPS &sps=*(m_pcSlice->getSPS());
  Int  iMvShift = 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  iMvShift += VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
  Int iOffset = 8;
  Int iHorMax = ( sps.getPicWidthInLumaSamples() + iOffset - (Int)m_uiCUPelX - 1 ) << iMvShift;
#if JVET_C0024_QTBT
  Int iHorMin = (      -(Int)sps.getCTUSize() - iOffset - (Int)m_uiCUPelX + 1 ) << iMvShift;
#else
  Int iHorMin = (      -(Int)sps.getMaxCUWidth() - iOffset - (Int)m_uiCUPelX + 1 ) << iMvShift;
#endif

  Int iVerMax = ( sps.getPicHeightInLumaSamples() + iOffset - (Int)m_uiCUPelY - 1 ) << iMvShift;
#if JVET_C0024_QTBT
  Int iVerMin = (      -(Int)sps.getCTUSize() - iOffset - (Int)m_uiCUPelY + 1 ) << iMvShift;
#else
  Int iVerMin = (      -(Int)sps.getMaxCUHeight() - iOffset - (Int)m_uiCUPelY + 1 ) << iMvShift;
#endif

  rcMv.setHor( min (iHorMax, max (iHorMin, rcMv.getHor())) );
  rcMv.setVer( min (iVerMax, max (iVerMin, rcMv.getVer())) );
}

#if VCEG_AZ08_INTRA_KLT
Void TComDataCU::clipMvIntraConstraint(Int regionId, Int &iHorMin, Int &iHorMax, Int &iVerMin, Int &iVerMax, Int iRange, UInt uiTemplateSize, UInt uiBlkSize, Int iCurrY, Int iCurrX, Int offsetLCUY, Int offsetLCUX)
{
    Int  iMvShift = 0;
    Int iTemplateSize = uiTemplateSize;
    Int iBlkSize = uiBlkSize;
    if (regionId == 0) //above outside LCU
    {
        iHorMax = min((iCurrX + iRange) << iMvShift, (Int)((m_pcSlice->getSPS()->getPicWidthInLumaSamples() - iBlkSize) << iMvShift));
        iHorMin = max((iTemplateSize) << iMvShift, (iCurrX - offsetLCUX - iBlkSize + 1) << iMvShift);

        iVerMax = (iCurrY - iBlkSize - offsetLCUY) << iMvShift;
        iVerMin = max(((iTemplateSize) << iMvShift), ((iCurrY - iRange) << iMvShift));

        iHorMin = iHorMin - iCurrX;
        iHorMax = iHorMax - iCurrX;
        iVerMax = iVerMax - iCurrY;
        iVerMin = iVerMin - iCurrY;
    }
    else if (regionId == 1) //left outside LCU
    {
        iHorMax = (iCurrX - offsetLCUX - iBlkSize) << iMvShift;
        iHorMin = max((iTemplateSize) << iMvShift, (iCurrX - iRange) << iMvShift);

        iVerMin = max((iTemplateSize) << iMvShift, (iCurrY - iBlkSize - offsetLCUY) << iMvShift);
        iVerMax = (iCurrY) << iMvShift;

        iHorMin = iHorMin - iCurrX;
        iHorMax = iHorMax - iCurrX;
        iVerMax = iVerMax - iCurrY;
        iVerMin = iVerMin - iCurrY;
    }
    else if (regionId == 2) //left outside LCU (can reach the bottom row of LCU)
    {
        iHorMin = max((iTemplateSize) << iMvShift, (iCurrX - iRange) << iMvShift);
        iHorMax = (iCurrX - offsetLCUX - iBlkSize) << iMvShift;
        iVerMin = (iCurrY + 1) << iMvShift;
#if JVET_C0024_QTBT
        iVerMax = min(m_pcSlice->getSPS()->getPicHeightInLumaSamples() - iBlkSize, (iCurrY - offsetLCUY + m_pcSlice->getSPS()->getCTUSize() - iBlkSize) << iMvShift);
#else
        iVerMax = min(m_pcSlice->getSPS()->getPicHeightInLumaSamples() - iBlkSize, (iCurrY - offsetLCUY + m_pcSlice->getSPS()->getMaxCUHeight() - iBlkSize) << iMvShift);
#endif

        iHorMin = iHorMin - iCurrX;
        iHorMax = iHorMax - iCurrX;
        iVerMax = iVerMax - iCurrY;
        iVerMin = iVerMin - iCurrY;
    }
}
#endif

UInt TComDataCU::getIntraSizeIdx(UInt uiAbsPartIdx)
{
#if JVET_C0024_QTBT
  UInt uiShift = 0;
#else
  UInt uiShift = ( m_pePartSize[uiAbsPartIdx]==SIZE_NxN ? 1 : 0 );
#endif

#if COM16_C806_LARGE_CTU
  UShort
#else
  UChar 
#endif
#if JVET_C0024_QTBT
#if JVET_C0024_CTU_256
    uiWidth = getWidth(uiAbsPartIdx)>>uiShift;
#else
    uiWidth = m_puhWidth[getTextType()][uiAbsPartIdx]>>uiShift;
#endif
#else
    uiWidth = m_puhWidth[uiAbsPartIdx]>>uiShift;
#endif
  UInt  uiCnt = 0;
  while( uiWidth )
  {
    uiCnt++;
    uiWidth>>=1;
  }
  uiCnt-=2;
  return uiCnt > 6 ? 6 : uiCnt;
}

Void TComDataCU::clearCbf( UInt uiIdx, ComponentID compID, UInt uiNumParts )
{
  memset( &m_puhCbf[compID][uiIdx], 0, sizeof(UChar)*uiNumParts);
}

/** Set a I_PCM flag for all sub-partitions of a partition.
 * \param bIpcmFlag I_PCM flag
 * \param uiAbsPartIdx patition index
 * \param uiDepth CU depth
 * \returns Void
 */
Void TComDataCU::setIPCMFlagSubParts  (Bool bIpcmFlag, UInt uiAbsPartIdx, UInt uiDepth)
{
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

  memset(m_pbIPCMFlag + uiAbsPartIdx, bIpcmFlag, sizeof(Bool)*uiCurrPartNumb );
}

/** Test whether the block at uiPartIdx is skipped.
 * \param uiPartIdx Partition index
 * \returns true if the current the block is skipped
 */
Bool TComDataCU::isSkipped( UInt uiPartIdx )
{
  return ( getSkipFlag( uiPartIdx ) );
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

Bool TComDataCU::xAddMVPCand( AMVPInfo* pInfo, RefPicList eRefPicList, Int iRefIdx, UInt uiPartUnitIdx, MVP_DIR eDir
#if COM16_C1016_AFFINE
  , bool bAffine
#endif
  )
{
  TComDataCU* pcTmpCU = NULL;
  UInt uiIdx;
  switch( eDir )
  {
    case MD_LEFT:
    {
      pcTmpCU = getPULeft(uiIdx, uiPartUnitIdx);
      break;
    }
    case MD_ABOVE:
    {
      pcTmpCU = getPUAbove(uiIdx, uiPartUnitIdx);
      break;
    }
    case MD_ABOVE_RIGHT:
    {
      pcTmpCU = getPUAboveRight(uiIdx, uiPartUnitIdx);
      break;
    }
    case MD_BELOW_LEFT:
    {
      pcTmpCU = getPUBelowLeft(uiIdx, uiPartUnitIdx);
      break;
    }
    case MD_ABOVE_LEFT:
    {
      pcTmpCU = getPUAboveLeft(uiIdx, uiPartUnitIdx);
      break;
    }
    default:
    {
      break;
    }
  }

  if ( pcTmpCU == NULL )
  {
    return false;
  }

  if ( pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) >= 0 && m_pcSlice->getRefPic( eRefPicList, iRefIdx)->getPOC() == pcTmpCU->getSlice()->getRefPOC( eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) ))
  {
    TComMv cMvPred = pcTmpCU->getCUMvField(eRefPicList)->getMv(uiIdx);

#if COM16_C1016_AFFINE
    // Unique
    if ( bAffine )
    {
      Int i = 0;
      for ( i=0; i<pInfo->iN; i++ )
      {
        if ( pInfo->m_acMvCand[i] == cMvPred )
          break;
      }
      if ( i == pInfo->iN )
      {
        pInfo->m_acMvCand[pInfo->iN++] = cMvPred;
        return true;
      }
    }
    else
    {
      pInfo->m_acMvCand[ pInfo->iN++] = cMvPred;
      return true;
    }
#else
    pInfo->m_acMvCand[ pInfo->iN++] = cMvPred;
    return true;
#endif
  }

  RefPicList eRefPicList2nd = REF_PIC_LIST_0;
  if(       eRefPicList == REF_PIC_LIST_0 )
  {
    eRefPicList2nd = REF_PIC_LIST_1;
  }
  else if ( eRefPicList == REF_PIC_LIST_1)
  {
    eRefPicList2nd = REF_PIC_LIST_0;
  }


  Int iCurrRefPOC = m_pcSlice->getRefPic( eRefPicList, iRefIdx)->getPOC();
  Int iNeibRefPOC;


  if( pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) >= 0 )
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC( eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) );
    if( iNeibRefPOC == iCurrRefPOC ) // Same Reference Frame But Diff List//
    {
      TComMv cMvPred = pcTmpCU->getCUMvField(eRefPicList2nd)->getMv(uiIdx);

#if COM16_C1016_AFFINE
      // Unique
      if ( bAffine )
      {
        Int i = 0;
        for ( i=0; i<pInfo->iN; i++ )
        {
          if ( pInfo->m_acMvCand[i] == cMvPred )
            break;
        }
        if ( i == pInfo->iN )
        {
          pInfo->m_acMvCand[pInfo->iN++] = cMvPred;
          return true;
        }
      }
      else
      {
        pInfo->m_acMvCand[ pInfo->iN++] = cMvPred;
        return true;
      }
#else
      pInfo->m_acMvCand[ pInfo->iN++] = cMvPred;
      return true;
#endif
    }
  }
  return false;
}

/**
 * \param pInfo
 * \param eRefPicList
 * \param iRefIdx
 * \param uiPartUnitIdx
 * \param eDir
 * \returns Bool
 */
Bool TComDataCU::xAddMVPCandOrder( AMVPInfo* pInfo, RefPicList eRefPicList, Int iRefIdx, UInt uiPartUnitIdx, MVP_DIR eDir 
#if COM16_C1016_AFFINE
  , bool bAffine
#endif
  )
{
  TComDataCU* pcTmpCU = NULL;
  UInt uiIdx;
  switch( eDir )
  {
  case MD_LEFT:
    {
      pcTmpCU = getPULeft(uiIdx, uiPartUnitIdx);
      break;
    }
  case MD_ABOVE:
    {
      pcTmpCU = getPUAbove(uiIdx, uiPartUnitIdx);
      break;
    }
  case MD_ABOVE_RIGHT:
    {
      pcTmpCU = getPUAboveRight(uiIdx, uiPartUnitIdx);
      break;
    }
  case MD_BELOW_LEFT:
    {
      pcTmpCU = getPUBelowLeft(uiIdx, uiPartUnitIdx);
      break;
    }
  case MD_ABOVE_LEFT:
    {
      pcTmpCU = getPUAboveLeft(uiIdx, uiPartUnitIdx);
      break;
    }
  default:
    {
      break;
    }
  }

  if ( pcTmpCU == NULL )
  {
    return false;
  }

  RefPicList eRefPicList2nd = REF_PIC_LIST_0;
  if(       eRefPicList == REF_PIC_LIST_0 )
  {
    eRefPicList2nd = REF_PIC_LIST_1;
  }
  else if ( eRefPicList == REF_PIC_LIST_1)
  {
    eRefPicList2nd = REF_PIC_LIST_0;
  }

  Int iCurrPOC = m_pcSlice->getPOC();
  Int iCurrRefPOC = m_pcSlice->getRefPic( eRefPicList, iRefIdx)->getPOC();
  Int iNeibPOC = iCurrPOC;
  Int iNeibRefPOC;
  Bool bIsCurrRefLongTerm = m_pcSlice->getRefPic( eRefPicList, iRefIdx)->getIsLongTerm();
  Bool bIsNeibRefLongTerm = false;

  //---------------  V1 (END) ------------------//
  if( pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) >= 0)
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC( eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) );
    TComMv cMvPred = pcTmpCU->getCUMvField(eRefPicList)->getMv(uiIdx);
    TComMv rcMv;

    bIsNeibRefLongTerm = pcTmpCU->getSlice()->getRefPic( eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) )->getIsLongTerm();
    if ( bIsCurrRefLongTerm == bIsNeibRefLongTerm )
    {
      if ( bIsCurrRefLongTerm || bIsNeibRefLongTerm )
      {
        rcMv = cMvPred;
      }
      else
      {
        Int iScale = xGetDistScaleFactor( iCurrPOC, iCurrRefPOC, iNeibPOC, iNeibRefPOC );
        if ( iScale == 4096 )
        {
          rcMv = cMvPred;
        }
        else
        {
          rcMv = cMvPred.scaleMv( iScale );
        }
      }

#if COM16_C1016_AFFINE
      // Unique
      if ( bAffine )
      {
        Int i = 0;
        for ( i=0; i<pInfo->iN; i++ )
        {
          if ( pInfo->m_acMvCand[i] == rcMv )
            break;
        }
        if ( i == pInfo->iN )
        {
          pInfo->m_acMvCand[pInfo->iN++] = rcMv;
          return true;
        }
      }
      else
      {
        pInfo->m_acMvCand[ pInfo->iN++] = rcMv;
        return true;
      }
#else
      pInfo->m_acMvCand[ pInfo->iN++] = rcMv;
      return true;
#endif
    }
  }
  //---------------------- V2(END) --------------------//
  if( pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) >= 0)
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC( eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) );
    TComMv cMvPred = pcTmpCU->getCUMvField(eRefPicList2nd)->getMv(uiIdx);
    TComMv rcMv;

    bIsNeibRefLongTerm = pcTmpCU->getSlice()->getRefPic( eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) )->getIsLongTerm();
    if ( bIsCurrRefLongTerm == bIsNeibRefLongTerm )
    {
      if ( bIsCurrRefLongTerm || bIsNeibRefLongTerm )
      {
        rcMv = cMvPred;
      }
      else
      {
        Int iScale = xGetDistScaleFactor( iCurrPOC, iCurrRefPOC, iNeibPOC, iNeibRefPOC );
        if ( iScale == 4096 )
        {
          rcMv = cMvPred;
        }
        else
        {
          rcMv = cMvPred.scaleMv( iScale );
        }
      }

#if COM16_C1016_AFFINE
      // Unique
      if ( bAffine )
      {
        Int i = 0;
        for ( i=0; i<pInfo->iN; i++ )
        {
          if ( pInfo->m_acMvCand[i] == rcMv )
            break;
        }
        if ( i == pInfo->iN )
        {
          pInfo->m_acMvCand[pInfo->iN++] = rcMv;
          return true;
        }
      }
      else
      {
        pInfo->m_acMvCand[ pInfo->iN++] = rcMv;
        return true;
      }
#else
      pInfo->m_acMvCand[ pInfo->iN++] = rcMv;
      return true;
#endif
    }
  }
  //---------------------- V3(END) --------------------//
  return false;
}

#if VCEG_AZ06_IC
Bool TComDataCU::xGetColMVP( RefPicList eRefPicList, Int ctuRsAddr, Int uiPartUnitIdx, TComMv& rcMv, Int& riRefIdx, Bool* bICFlag )
#else
Bool TComDataCU::xGetColMVP( RefPicList eRefPicList, Int ctuRsAddr, Int uiPartUnitIdx, TComMv& rcMv, Int& riRefIdx )
#endif
{
  UInt uiAbsPartAddr = uiPartUnitIdx;

  RefPicList  eColRefPicList;
  Int iColPOC, iColRefPOC, iCurrPOC, iCurrRefPOC, iScale;
  TComMv cColMv;

  // use coldir.
  TComPic *pColPic = getSlice()->getRefPic( RefPicList(getSlice()->isInterB() ? 1-getSlice()->getColFromL0Flag() : 0), getSlice()->getColRefIdx());
  TComDataCU *pColCtu = pColPic->getCtu( ctuRsAddr );
#if JVET_C0024_QTBT
  if(pColCtu->getPic()==0)
#else
  if(pColCtu->getPic()==0||pColCtu->getPartitionSize(uiPartUnitIdx)==NUMBER_OF_PART_SIZES)
#endif
  {
    return false;
  }
  iCurrPOC = m_pcSlice->getPOC();
  iColPOC = pColCtu->getSlice()->getPOC();

  if (!pColCtu->isInter(uiAbsPartAddr))
  {
    return false;
  }

  eColRefPicList = getSlice()->getCheckLDC() ? eRefPicList : RefPicList(getSlice()->getColFromL0Flag());

  Int iColRefIdx = pColCtu->getCUMvField(RefPicList(eColRefPicList))->getRefIdx(uiAbsPartAddr);

  if (iColRefIdx < 0 )
  {
    eColRefPicList = RefPicList(1 - eColRefPicList);
    iColRefIdx = pColCtu->getCUMvField(RefPicList(eColRefPicList))->getRefIdx(uiAbsPartAddr);

    if (iColRefIdx < 0 )
    {
      return false;
    }
  }

  // Scale the vector.
  iColRefPOC = pColCtu->getSlice()->getRefPOC(eColRefPicList, iColRefIdx);
  cColMv = pColCtu->getCUMvField(eColRefPicList)->getMv(uiAbsPartAddr);

  iCurrRefPOC = m_pcSlice->getRefPic(eRefPicList, riRefIdx)->getPOC();

  Bool bIsCurrRefLongTerm = m_pcSlice->getRefPic(eRefPicList, riRefIdx)->getIsLongTerm();
  Bool bIsColRefLongTerm = pColCtu->getSlice()->getIsUsedAsLongTerm(eColRefPicList, iColRefIdx);

  if ( bIsCurrRefLongTerm != bIsColRefLongTerm )
  {
    return false;
  }

#if VCEG_AZ06_IC
  if( bICFlag != NULL )
  {
    *bICFlag = pColCtu->getICFlag( uiAbsPartAddr ) ;
  }
#endif
  if ( bIsCurrRefLongTerm || bIsColRefLongTerm )
  {
    rcMv = cColMv;
  }
  else
  {
    iScale = xGetDistScaleFactor(iCurrPOC, iCurrRefPOC, iColPOC, iColRefPOC);
    if ( iScale == 4096 )
    {
      rcMv = cColMv;
    }
    else
    {
      rcMv = cColMv.scaleMv( iScale );
    }
  }

  return true;
}

#if VCEG_AZ07_FRUC_MERGE
TComMv TComDataCU::scaleMV( const TComMv & rColMV , Int iCurrPOC, Int iCurrRefPOC, Int iColPOC, Int iColRefPOC )
{
  TComMv mv = rColMV;
  Int iScale = xGetDistScaleFactor(iCurrPOC, iCurrRefPOC, iColPOC, iColRefPOC);
  if ( iScale != 4096 )
  {
    mv = rColMV.scaleMv( iScale );
  }
  return( mv );
}

Bool TComDataCU::isSameMVField( RefPicList eListA , TComMvField & rMVFieldA , RefPicList eListB , TComMvField & rMVFieldB )
{
  if( rMVFieldA.getRefIdx() >= 0 && rMVFieldB.getRefIdx() >= 0 )
  {
    return( rMVFieldA.getMv() == rMVFieldB.getMv()
      && getSlice()->getRefPOC( eListA , rMVFieldA.getRefIdx() ) == getSlice()->getRefPOC( eListB , rMVFieldB.getRefIdx() ) );
  }
  else
  {
    return( false );
  }
}
#endif

Int TComDataCU::xGetDistScaleFactor(Int iCurrPOC, Int iCurrRefPOC, Int iColPOC, Int iColRefPOC)
{
  Int iDiffPocD = iColPOC - iColRefPOC;
  Int iDiffPocB = iCurrPOC - iCurrRefPOC;

  if( iDiffPocD == iDiffPocB )
  {
    return 4096;
  }
  else
  {
    Int iTDB      = Clip3( -128, 127, iDiffPocB );
    Int iTDD      = Clip3( -128, 127, iDiffPocD );
#if VCEG_AZ07_FRUC_MERGE
    Int iScale = getSlice()->getScaleFactor( iTDB , iTDD );
#else
    Int iX        = (0x4000 + abs(iTDD/2)) / iTDD;
    Int iScale    = Clip3( -4096, 4095, (iTDB * iX + 32) >> 6 );
#endif
    return iScale;
  }
}

Void TComDataCU::xDeriveCenterIdx( UInt uiPartIdx, UInt& ruiPartIdxCenter )
{
  UInt uiPartAddr;
  Int  iPartWidth;
  Int  iPartHeight;
  getPartIndexAndSize( uiPartIdx, uiPartAddr, iPartWidth, iPartHeight);

  ruiPartIdxCenter = m_absZIdxInCtu+uiPartAddr; // partition origin.
  ruiPartIdxCenter = g_auiRasterToZscan[ g_auiZscanToRaster[ ruiPartIdxCenter ]
                                        + ( iPartHeight/m_pcPic->getMinCUHeight()  )/2*m_pcPic->getNumPartInCtuWidth()
                                        + ( iPartWidth/m_pcPic->getMinCUWidth()  )/2];
}

Void TComDataCU::compressMV()
{
  Int scaleFactor = 4 * AMVP_DECIMATION_FACTOR / m_unitSize;
  if (scaleFactor > 0)
  {
    for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
    {
      m_acCUMvField[i].compress(m_pePredMode, scaleFactor);
    }
  }
}

#if DIMD_NUM_INTRA_DIR_INC || DIMD_INTRA_PRED
UInt TComDataCU::getCoefScanIdx(const UInt uiAbsPartIdx, const UInt uiWidth, const UInt uiHeight, const ComponentID compID)
#else
UInt TComDataCU::getCoefScanIdx(const UInt uiAbsPartIdx, const UInt uiWidth, const UInt uiHeight, const ComponentID compID) const
#endif
{
  //------------------------------------------------

  //this mechanism is available for intra only

  if (!isIntra(uiAbsPartIdx))
  {
    return SCAN_DIAG;
  }

  //------------------------------------------------

  //check that MDCS can be used for this TU

  const ChromaFormat format = getPic()->getChromaFormat();

  const UInt maximumWidth  = MDCS_MAXIMUM_WIDTH  >> getComponentScaleX(compID, format);
  const UInt maximumHeight = MDCS_MAXIMUM_HEIGHT >> getComponentScaleY(compID, format);

  if ((uiWidth > maximumWidth) || (uiHeight > maximumHeight))
  {
    return SCAN_DIAG;
  }

  //------------------------------------------------

  //otherwise, select the appropriate mode

  UInt uiDirMode  = getIntraDir(toChannelType(compID), uiAbsPartIdx);
#if !JVET_E0062_MULTI_DMS
  if (uiDirMode==DM_CHROMA_IDX)
  {
#if JVET_C0024_QTBT
    assert(isChroma(compID));
    uiDirMode = getSlice()->isIntra() ? getPic()->getCtu(getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, getZorderIdxInCtu()+uiAbsPartIdx): getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
#else
    const TComSPS *sps=getSlice()->getSPS();
    const UInt partsPerMinCU = 1<<(2*(sps->getMaxTotalCUDepth() - sps->getLog2DiffMaxMinCodingBlockSize()));
    uiDirMode = getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, getPic()->getChromaFormat(), partsPerMinCU));
#endif
  }
#endif

  if (isChroma(compID) && (format == CHROMA_422))
  {
    uiDirMode = g_chroma422IntraAngleMappingTable[uiDirMode];
  }

  //------------------

#if DIMD_INTRA_PRED
  if (getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx))
  {
    uiDirMode = getSecIntraDir(uiAbsPartIdx);
    assert(uiDirMode >= PLANAR_IDX && uiDirMode < NUM_INTRA_MODE - 1);
    if (abs((Int)uiDirMode - VER_IDX) <= MDCS_ANGLE_LIMIT)
    {
      return SCAN_HOR;
    }
    else if (abs((Int)uiDirMode - HOR_IDX) <= MDCS_ANGLE_LIMIT)
    {
      return SCAN_VER;
    }
    else
    {
      return SCAN_DIAG;
    }
  }
#endif

#if DIMD_NUM_INTRA_DIR_INC
  if (getDIMDEnabledFlag(toChannelType(compID), uiAbsPartIdx))
  {
    uiDirMode = getExtIntraDir(toChannelType(compID), uiAbsPartIdx);
    if (abs((Int)uiDirMode - EXT_VER_IDX) <= MDCS_EXT_ANGLE_LIMIT)
    {
      return SCAN_HOR;
    }
    else if (abs((Int)uiDirMode - EXT_HOR_IDX) <= MDCS_EXT_ANGLE_LIMIT)
    {
      return SCAN_VER;
    }
    else
    {
      return SCAN_DIAG;
    }
  }
  else
  {
    assert(uiDirMode >= PLANAR_IDX && uiDirMode < NUM_INTRA_MODE);
#endif
  if      (abs((Int)uiDirMode - VER_IDX) <= MDCS_ANGLE_LIMIT)
  {
    return SCAN_HOR;
  }
  else if (abs((Int)uiDirMode - HOR_IDX) <= MDCS_ANGLE_LIMIT)
  {
    return SCAN_VER;
  }
  else
  {
    return SCAN_DIAG;
  }
#if DIMD_NUM_INTRA_DIR_INC
  }
#endif
}

#if ALF_HM3_REFACTOR
UInt TComDataCU::getCtxAlfCtrlFlag( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

  // Get BCBP of left PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(getZorderIdxInCtu(), uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, getZorderIdxInCtu() + uiAbsPartIdx );
#endif
  uiCtx    = ( pcTempCU ) ? pcTempCU->getAlfCtrlFlag( uiTempPartIdx ) : 0;

  // Get BCBP of above PU
#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(getZorderIdxInCtu() , uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, getZorderIdxInCtu() + uiAbsPartIdx );
#endif
  uiCtx   += ( pcTempCU ) ? pcTempCU->getAlfCtrlFlag( uiTempPartIdx ) : 0;

  return uiCtx;
}

#if JVET_C0024_QTBT
Void TComDataCU::setAlfCtrlFlagSubParts         ( UChar uhFlag, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight )
{
    setSubPart<UChar>( uhFlag, m_puhAlfCtrlFlag, uiAbsPartIdx, uiWidth, uiHeight );
}
#else
Void TComDataCU::setAlfCtrlFlagSubParts         ( UInt uiFlag, UInt uiAbsPartIdx, UInt uiDepth )
{
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);

  for (UInt ui = 0; ui < uiCurrPartNumb; ui++ )
  {
    m_puiAlfCtrlFlag[uiAbsPartIdx + ui] = uiFlag;
  }
}
#endif

Void TComDataCU::createTmpAlfCtrlFlag()
{
#if JVET_C0024_QTBT
  m_puhTmpAlfCtrlFlag = (UChar* )xMalloc(UChar, m_uiNumPartition);
#else
  m_puiTmpAlfCtrlFlag = (UInt* )xMalloc(UInt, m_uiNumPartition);
#endif
}

Void TComDataCU::destroyTmpAlfCtrlFlag()
{
#if JVET_C0024_QTBT
  if(m_puhTmpAlfCtrlFlag)
  {
    xFree(m_puhTmpAlfCtrlFlag);        m_puhTmpAlfCtrlFlag = NULL;
  }
#else
  if(m_puiTmpAlfCtrlFlag)
  {
    xFree(m_puiTmpAlfCtrlFlag);        m_puiTmpAlfCtrlFlag = NULL;
  }
#endif
}

Void TComDataCU::copyAlfCtrlFlagToTmp()
{
#if JVET_C0024_QTBT
  memcpy( m_puhTmpAlfCtrlFlag, m_puhAlfCtrlFlag, sizeof(UChar)*m_uiNumPartition );
#else
  memcpy( m_puiTmpAlfCtrlFlag, m_puiAlfCtrlFlag, sizeof(UInt)*m_uiNumPartition );
#endif
}

Void TComDataCU::copyAlfCtrlFlagFromTmp()
{
#if JVET_C0024_QTBT
  memcpy( m_puhAlfCtrlFlag, m_puhTmpAlfCtrlFlag, sizeof(UChar)*m_uiNumPartition );
#else
  memcpy( m_puiAlfCtrlFlag, m_puiTmpAlfCtrlFlag, sizeof(UInt)*m_uiNumPartition );
#endif
}
#endif

#if COM16_C806_OBMC
// Function for fetching neighboring motions.
Bool TComDataCU::getNeigMotion( UInt uiAbsPartIdx, TComMvField cNeigMvField[2], Int &irNeigPredDir, Int iDir, TComMvField cCurMvField[2], Int &iCurrDir, UInt uiZeroIdx, Bool &bTobeStored )
{
  TComDataCU* pcTmpCU = NULL;
  UInt uiIdx;

#if TT_TREE_OBMC
  if( iDir == 0 ) //above
  {
    pcTmpCU = getPUAbove( uiIdx, deriveZIdxInCTU(uiZeroIdx,uiAbsPartIdx) );
  }  
  else if( iDir == 1 ) //left
  {
    pcTmpCU = getPULeft( uiIdx, deriveZIdxInCTU(uiZeroIdx,uiAbsPartIdx) );
  }
  else if( iDir == 2 ) //below
  {
    pcTmpCU = this;
    uiIdx = deriveZIdxInCU(uiZeroIdx, g_auiRasterToZscan[g_auiZscanToRaster[deriveZIdxInCTU(uiZeroIdx,uiAbsPartIdx)] + getPic()->getNumPartInCtuWidth()]);
  }
  else if( iDir == 3 ) //right
  {
    pcTmpCU = this;
    uiIdx = deriveZIdxInCU(uiZeroIdx, g_auiRasterToZscan[g_auiZscanToRaster[deriveZIdxInCTU(uiZeroIdx,uiAbsPartIdx)] + 1]);
  }
#else
  if( iDir == 0 ) //above
  {
    pcTmpCU = getPUAbove( uiIdx, uiAbsPartIdx + uiZeroIdx );
  }  
  else if( iDir == 1 ) //left
  {
    pcTmpCU = getPULeft( uiIdx, uiAbsPartIdx + uiZeroIdx );
  }
  else if( iDir == 2 ) //below
  {
    pcTmpCU = this;
    uiIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx + uiZeroIdx] + getPic()->getNumPartInCtuWidth()] - uiZeroIdx;
  }
  else if( iDir == 3 ) //right
  {
    pcTmpCU = this;
    uiIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiAbsPartIdx + uiZeroIdx] + 1] - uiZeroIdx;
  }
#endif

  if( pcTmpCU == NULL || pcTmpCU->isIntra( uiIdx ) )
  {
    return false;
  }

  irNeigPredDir = pcTmpCU->getInterDir( uiIdx );
  if(irNeigPredDir)
  { 
    if( !bTobeStored )
    {
      //backup motion information
      for(UInt iRefList = 0; iRefList < 2; iRefList ++)
      {
        TComCUMvField* pTmpMvField = getCUMvField( RefPicList(iRefList) );
        cCurMvField[iRefList].setMvField( pTmpMvField->getMv( uiAbsPartIdx ), pTmpMvField->getRefIdx( uiAbsPartIdx ) );
      }
      iCurrDir = getInterDir( uiAbsPartIdx );
      bTobeStored = true;
    }
    for(UInt iRefList = 0; iRefList < 2; iRefList ++)
    {
      TComCUMvField* pTmpMvField = pcTmpCU->getCUMvField( RefPicList(iRefList) );
      cNeigMvField[iRefList].setMvField( pTmpMvField->getMv( uiIdx ), pTmpMvField->getRefIdx( uiIdx ) );
    }

    if(irNeigPredDir != iCurrDir)
    {
      return true;
    }
    else
    {
      for(UInt iRefList = 0; iRefList < 2; iRefList ++)
      {
        if( iCurrDir & ( 1 << iRefList ) )
        {       
          if(!(cCurMvField[iRefList] == cNeigMvField[iRefList]))
          {
            return true;
          }
        }
      }
      return false;
    }
  }
  else
  {
    return false;
  }
}
#endif

#if COM16_C806_LARGE_CTU
Void TComDataCU::getMaxMinCUDepth( UChar & rucMinDepth , UChar & rucMaxDepth , UInt uiAbsPartIdx )
{
#if JVET_C0024_QTBT
  const UChar ucMaxCUDepth = g_aucConvertToBit[getSlice()->getSPS()->getCTUSize()] 
  - g_aucConvertToBit[getSlice()->getSPS()->getMinQTSize(getSlice()->getSliceType(), getTextType())];
#else
  const UChar ucMaxCUDepth = ( UChar )( getSlice()->getSPS()->getLog2DiffMaxMinCodingBlockSize() );
#endif
  rucMinDepth = ucMaxCUDepth;
  rucMaxDepth = 0;

  TComDataCU * pNeighbor;
  UInt uiNeighbor; 

  // left
  pNeighbor = getPULeft( uiNeighbor , uiAbsPartIdx );
  if( pNeighbor != NULL )
  {
    rucMinDepth = min( rucMinDepth , pNeighbor->getDepth( uiNeighbor ) );
    rucMaxDepth = max( rucMaxDepth , pNeighbor->getDepth( uiNeighbor ) );
  }
  else
  {
    rucMinDepth = 0;
    rucMaxDepth = ucMaxCUDepth;
  }

  // bottom left
  pNeighbor = getPUBelowLeft( uiNeighbor , uiAbsPartIdx );
  if( pNeighbor != NULL )
  {
    rucMinDepth = min( rucMinDepth , pNeighbor->getDepth( uiNeighbor ) );
    rucMaxDepth = max( rucMaxDepth , pNeighbor->getDepth( uiNeighbor ) );
  }
  else
  {
    rucMinDepth = 0;
    rucMaxDepth = ucMaxCUDepth;
  }

  // top
  pNeighbor = getPUAbove( uiNeighbor , uiAbsPartIdx );
  if( pNeighbor != NULL )
  {
    rucMinDepth = min( rucMinDepth , pNeighbor->getDepth( uiNeighbor ) );
    rucMaxDepth = max( rucMaxDepth , pNeighbor->getDepth( uiNeighbor ) );
  }
  else
  {
    rucMinDepth = 0;
    rucMaxDepth = ucMaxCUDepth;
  }

  // top right
  pNeighbor = getPUAboveRight( uiNeighbor , uiAbsPartIdx );
  if( pNeighbor != NULL )
  {
    rucMinDepth = min( rucMinDepth , pNeighbor->getDepth( uiNeighbor ) );
    rucMaxDepth = max( rucMaxDepth , pNeighbor->getDepth( uiNeighbor ) );
  }
  else
  {
    rucMinDepth = 0;
    rucMaxDepth = ucMaxCUDepth;
  }

  Int nDepthInc = 1;
  rucMinDepth = rucMinDepth >= nDepthInc ? rucMinDepth - nDepthInc : 0;
  rucMaxDepth = min( ( UChar )ucMaxCUDepth , ( UChar )( rucMaxDepth + 1 ) );
}
#endif

#if COM16_C1045_BIO_HARMO_IMPROV
Bool TComDataCU::isBIOLDB( UInt uiAbsPartIdx )
{
  Bool BIOLDB = false;
#if  JVET_C0027_BIO
  if (!getSlice()->getBioLDBPossible() ) return( BIOLDB );
#endif
  TComCUMvField * pCuMvField0 = getCUMvField( REF_PIC_LIST_0 );
  TComCUMvField * pCuMvField1 = getCUMvField( REF_PIC_LIST_1 );
  if( getSlice()->getCheckLDC()
    && pCuMvField0->getRefIdx( uiAbsPartIdx ) >= 0 && pCuMvField1->getRefIdx( uiAbsPartIdx ) >= 0 )
  {
    Int pocCur = getSlice()->getPOC();
    Int poc0 = getSlice()->getRefPOC( REF_PIC_LIST_0 , pCuMvField0->getRefIdx( uiAbsPartIdx ) );
    Int poc1 = getSlice()->getRefPOC( REF_PIC_LIST_1 , pCuMvField1->getRefIdx( uiAbsPartIdx ) );
    if( poc0 != poc1 && ( poc0 - pocCur ) * ( poc1 - pocCur ) > 0 )
    {
      Int dT0 = poc0 - pocCur;
      Int dT1 = poc1 - pocCur;
      Bool zeroMv0 = ( pCuMvField0->getMv( uiAbsPartIdx ).getAbsHor() + pCuMvField0->getMv( uiAbsPartIdx ).getAbsVer() ) == 0;
      Bool zeroMv1 = ( pCuMvField1->getMv( uiAbsPartIdx ).getAbsHor() + pCuMvField1->getMv( uiAbsPartIdx ).getAbsVer() ) == 0;
      if( !zeroMv0 && !zeroMv1 )
      {
        BIOLDB = dT0 * pCuMvField1->getMv( uiAbsPartIdx ).getHor() == dT1 * pCuMvField0->getMv( uiAbsPartIdx ).getHor() 
          && dT0 * pCuMvField1->getMv( uiAbsPartIdx ).getVer() == dT1 * pCuMvField0->getMv( uiAbsPartIdx ).getVer();
      }
    }
  }
  return( BIOLDB );
}
#endif

#if COM16_C1016_AFFINE
/** Test whether the current block is affine mode
 * \param uiPartIdx Block index
 * \returns Flag indicating whether the block is affine
 */
Bool TComDataCU::isAffine( UInt uiPartIdx )
{
  return ( getAffineFlag( uiPartIdx ) );
}

UInt TComDataCU::getCtxAffineFlag( UInt uiAbsPartIdx )
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx    = ( pcTempCU ) ? pcTempCU->isAffine( uiTempPartIdx ) : 0;

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove( uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove( uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx );
#endif
  uiCtx   += ( pcTempCU ) ? pcTempCU->isAffine( uiTempPartIdx ) : 0;

  return uiCtx;
}

#if JVECT_C0062_AFFINE_SIX_PARAM
UInt TComDataCU::getCtxAffineParamFlag(UInt uiAbsPartIdx)
{
  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;
  UInt        uiCtx = 0;

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPULeft(uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPULeft(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
#endif
  uiCtx = (pcTempCU) ? (pcTempCU->isAffine(uiTempPartIdx) && pcTempCU->getAffineParamFlag(uiTempPartIdx)) : 0;

#if TT_TREE_DS1 // deriveZIdxInCTU
  pcTempCU = getPUAbove(uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) );
#else
  pcTempCU = getPUAbove(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
#endif
  uiCtx += (pcTempCU) ? (pcTempCU->isAffine(uiTempPartIdx) && pcTempCU->getAffineParamFlag(uiTempPartIdx)) : 0;

  return uiCtx;
}

Void TComDataCU::setAffineParamFlagSubParts(Bool bAffineParamFlag, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  setSubPart(bAffineParamFlag, m_affineParamFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
#else
  setSubPart(bAffineParamFlag, m_affineParamFlag, uiAbsPartIdx, uiDepth, uiPartIdx);
#endif
}
#endif

Void TComDataCU::setAffineFlagSubParts( Bool bAffineFlag, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth )
{
#if JVET_C0024_QTBT
  setSubPart( bAffineFlag, m_affineFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx) );
#else
  setSubPart( bAffineFlag, m_affineFlag, uiAbsPartIdx, uiDepth, uiPartIdx );
#endif
}

Void TComDataCU::setAllAffineMvField( UInt uiAbsPartIdx, UInt uiPuIdx, TComMvField *pcMvField, RefPicList eRefPicList, UInt uiDepth
#if JVECT_C0062_AFFINE_SIX_PARAM
  , Bool b6Param
#endif
)
{
  // Set Mv
  TComMv acMv[3];
  for ( Int i=0; i<3; i++ )
  {
    acMv[i] = pcMvField[i].getMv();
  }
  setAllAffineMv( uiAbsPartIdx, uiPuIdx, acMv, eRefPicList, uiDepth
#if JVECT_C0062_AFFINE_SIX_PARAM
    , b6Param
#endif
  );

  // Set RefIdx
  assert( pcMvField[0].getRefIdx() == pcMvField[1].getRefIdx() && pcMvField[0].getRefIdx() == pcMvField[2].getRefIdx() );
#if JVET_C0024_QTBT
  getCUMvField(eRefPicList)->setAllRefIdx( pcMvField[0].getRefIdx(), SIZE_2Nx2N, uiAbsPartIdx, uiDepth, uiPuIdx );
#else
  PartSize ePartSize = getPartitionSize( uiAbsPartIdx );
  getCUMvField(eRefPicList)->setAllRefIdx( pcMvField[0].getRefIdx(), ePartSize, uiAbsPartIdx, uiDepth, uiPuIdx );
#endif
}


#if AFFINE_UNIFY_MC
Void TComDataCU::deriveAffineMCBlkSize(TComMv acMv[3], Int width, Int height, Int& blockWidth, Int& blockHeight) { // same as in xPredAffineBlk
  Int mvWx = max(abs((acMv[1] - acMv[0]).getHor()), abs((acMv[1] - acMv[0]).getVer()));
  Int mvWy = max(abs((acMv[2] - acMv[0]).getHor()), abs((acMv[2] - acMv[0]).getVer()));

  Int iMvPrecision = 4;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
  iMvPrecision -= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
#if AFFINE_MC_PRECISION_FIX
  iMvPrecision = 0;
#endif

  if (mvWx)
  {
    blockWidth = max((Int)((width >> iMvPrecision) / mvWx), 1);
    while (width % blockWidth)
    {
      blockWidth--;
    }
    blockWidth = max(AFFINE_MIN_BLOCK_SIZE, blockWidth);
  }
  if (mvWy)
  {
    blockHeight = max((Int)((height >> iMvPrecision) / mvWy), 1);
    while (height % blockHeight)
    {
      blockHeight--;
    }
    blockHeight = max(AFFINE_MIN_BLOCK_SIZE, blockHeight);
  }
}
#endif

Void TComDataCU::setAllAffineMv ( UInt uiAbsPartIdx, UInt uiPuIdx, TComMv acMv[3], RefPicList e, UInt uiDepth
#if JVECT_C0062_AFFINE_SIX_PARAM
  , Bool b6Param
#endif
)
{
  Int iWidth  = getWidth(uiAbsPartIdx);
  Int iHeight = getHeight(uiAbsPartIdx);
  Int iPartW = iWidth/getPic()->getMinCUWidth();
  Int iPartH = iHeight/getPic()->getMinCUHeight();

#if JVECT_C0062_AFFINE_SIX_PARAM
  Bool bRecShape = ((iHeight / iWidth) >= 2);
#endif

  // Get partIdx of four corner
#if TT_TREE_AFFINE // deriveZIdxInCTU
  UInt uiPartIdxLT = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  UInt uiPartIdxLT = uiAbsPartIdx + m_absZIdxInCtu;
#endif
  UInt uiPartIdxRT = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + iPartW - 1 ];
  UInt uiPartIdxLB = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + (iPartH-1)*getPic()->getNumPartInCtuWidth() ];
  UInt uiPartIdxRB = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + iPartW - 1  + (iPartH-1)*getPic()->getNumPartInCtuHeight() ];

  // Set CUMvField for reference by neighbor PU
  TComMv cMv;
  UInt uiPartIdx;
  Int iBit = 6;

#if AFFINE_UNIFY_MC // make the block size consistent with xPredAffineBlk
  // get affine sub-block width and height
  Int blockWidth = iWidth;
  Int blockHeight = iHeight;

  deriveAffineMCBlkSize(acMv, iWidth, iHeight, blockWidth, blockHeight); // same way to derive the block size as in xPredAffineBlk

  // convert to 2^(2+iBit) precision
#if JVET_C0024_QTBT
  Int bitMinW = g_aucConvertToBit[getPic()->getMinCUWidth()] + MIN_CU_LOG2;
  Int bitMinH = g_aucConvertToBit[getPic()->getMinCUHeight()] + MIN_CU_LOG2;
#else
  Int bitMinW = g_aucConvertToBit[getPic()->getMinCUWidth()] + 2;
  Int bitMinH = g_aucConvertToBit[getPic()->getMinCUHeight()] + 2;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  Int iDMvHorX, iDMvVerX, iDMvHorY, iDMvVerY;
  if (b6Param)
  {
    iDMvHorX = ((acMv[1] - acMv[0]).getHor() << (iBit + bitMinW)) / iWidth;   // deltaMvHor
    iDMvHorY = ((acMv[1] - acMv[0]).getVer() << (iBit + bitMinW)) / iWidth;
    iDMvVerX = ((acMv[2] - acMv[0]).getHor() << (iBit + bitMinH)) / iHeight;  // deltaMvVer
    iDMvVerY = ((acMv[2] - acMv[0]).getVer() << (iBit + bitMinH)) / iHeight;
  }
  else
  {
    if (bRecShape)
    {
      Int vx1 = (acMv[2].getVer() - acMv[0].getVer())*iWidth / iHeight + acMv[0].getHor();
      Int vy1 = -(acMv[2].getHor() - acMv[0].getHor())*iWidth / iHeight + acMv[0].getVer();
      acMv[1].set(vx1, vy1);
    }
    else
    {
      Int vx2 = -(acMv[1].getVer() - acMv[0].getVer()) * iHeight / iWidth + acMv[0].getHor();
      Int vy2 = (acMv[1].getHor() - acMv[0].getHor()) * iHeight / iWidth + acMv[0].getVer();
      acMv[2].set(vx2, vy2);
    }
    iDMvHorX = ((acMv[1] - acMv[0]).getHor() << (iBit + bitMinW)) / iWidth;
    iDMvHorY = ((acMv[1] - acMv[0]).getVer() << (iBit + bitMinH)) / iWidth;
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }
#else
  Int iDMvHorX = ((acMv[1] - acMv[0]).getHor() << (iBit + bitMinW)) / iWidth;
  Int iDMvHorY = ((acMv[1] - acMv[0]).getVer() << (iBit + bitMinH)) / iWidth;
  Int iDMvVerX = -iDMvHorY;
  Int iDMvVerY = iDMvHorX;
#endif

  Int iMvScaleHor = acMv[0].getHor() << iBit;
  Int iMvScaleVer = acMv[0].getVer() << iBit;
  Int iMvYHor = iMvScaleHor;
  Int iMvYVer = iMvScaleVer;

  // Calculate Mv for 4x4 Part
  Int iMvScaleTmpHor, iMvScaleTmpVer;
  for (Int hh = 0; hh < iHeight; hh += blockHeight)
  {
    Int maxH = (hh + blockHeight) >> bitMinH;

    for (Int ww = 0; ww < iWidth; ww += blockWidth)
    {
      iMvScaleTmpHor = iMvScaleHor + (iDMvHorX * blockWidth >> 3) + (iDMvVerX * blockHeight >> 3);
      iMvScaleTmpVer = iMvScaleVer + (iDMvHorY * blockWidth >> 3) + (iDMvVerY * blockHeight >> 3);

      // get the MV in hevc define precision
      Int xHevc, yHevc;
      xHevc = iMvScaleTmpHor >> iBit;
      yHevc = iMvScaleTmpVer >> iBit;
      TComMv cMvHevc(xHevc, yHevc);
      Int maxW = (ww + blockWidth) >> bitMinW;

      for (Int h = (hh >> bitMinH); h < maxH; h++)
      {
        for (Int w = (ww >> bitMinW); w < maxW; w++)
        {
#if TT_TREE_AFFINE // deriveZIdxInCU
          uiPartIdx = deriveZIdxInCU(m_absZIdxInCtu, g_auiRasterToZscan[g_auiZscanToRaster[uiPartIdxLT] + w + h*getPic()->getNumPartInCtuWidth()]);
#else
          uiPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[uiPartIdxLT] + w + h*getPic()->getNumPartInCtuWidth()] - m_absZIdxInCtu;
#endif
          getCUMvField(e)->setMv(cMvHevc, uiPartIdx);

          // save corner mv
          if (ww == 0 && hh == 0)
          {
            getCUMvField(e)->setMv(acMv[0], uiPartIdx);
          }
          if (ww == (iWidth - blockWidth) && hh == 0)
          {
            getCUMvField(e)->setMv(acMv[1], uiPartIdx);
          }
          if (ww == 0 && hh == (iHeight - blockHeight))
          {
            getCUMvField(e)->setMv(acMv[2], uiPartIdx);
          }
        }
      }
      iMvScaleHor += iDMvHorX*(blockWidth >> bitMinW);  // switch from x to x+AffineBlockSize, add deltaMvHor
      iMvScaleVer += iDMvHorY*(blockWidth >> bitMinW);
    } // hh
    iMvYHor += iDMvVerX*(blockHeight >> bitMinH);       // switch from y to y+AffineBlockSize, add deltaMvVer
    iMvYVer += iDMvVerY*(blockHeight >> bitMinH);

    iMvScaleHor = iMvYHor;
    iMvScaleVer = iMvYVer;
  } // ww
#else
  // convert to 2^(2+iBit) precision
#if JVET_C0024_QTBT
  Int bitMinW = g_aucConvertToBit[getPic()->getMinCUWidth() ] + MIN_CU_LOG2;
  Int bitMinH = g_aucConvertToBit[getPic()->getMinCUHeight()] + MIN_CU_LOG2;
#else
  Int bitMinW = g_aucConvertToBit[getPic()->getMinCUWidth() ] + 2;
  Int bitMinH = g_aucConvertToBit[getPic()->getMinCUHeight()] + 2;
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
  Int iDMvHorX, iDMvVerX, iDMvHorY, iDMvVerY;
  if (b6Param)
  {
    iDMvHorX = ((acMv[1] - acMv[0]).getHor() << (iBit + bitMinW)) / iWidth;   // deltaMvHor
    iDMvHorY = ((acMv[1] - acMv[0]).getVer() << (iBit + bitMinW)) / iWidth;
    iDMvVerX = ((acMv[2] - acMv[0]).getHor() << (iBit + bitMinH)) / iHeight;  // deltaMvVer
    iDMvVerY = ((acMv[2] - acMv[0]).getVer() << (iBit + bitMinH)) / iHeight;
  }
  else
  {
    if (bRecShape)
    {
#if !AFFINE_CLIP_REMOVE
      clipMv(acMv[0]);
      clipMv(acMv[2]);
#endif
      Int vx1 = (acMv[2].getVer() - acMv[0].getVer())*iWidth / iHeight + acMv[0].getHor();
      Int vy1 = -(acMv[2].getHor() - acMv[0].getHor())*iWidth / iHeight + acMv[0].getVer();
      acMv[1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
      clipMv(acMv[1]);
#endif
    }
    else
    {
#if !AFFINE_CLIP_REMOVE
      clipMv(acMv[0]);
      clipMv(acMv[1]);
#endif
      Int vx2 = -(acMv[1].getVer() - acMv[0].getVer()) * iHeight / iWidth + acMv[0].getHor();
      Int vy2 = (acMv[1].getHor() - acMv[0].getHor()) * iHeight / iWidth + acMv[0].getVer();
      acMv[2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
      clipMv(acMv[2]);
#endif
    }
    iDMvHorX = ((acMv[1] - acMv[0]).getHor() << (iBit + bitMinW)) / iWidth;
    iDMvHorY = ((acMv[1] - acMv[0]).getVer() << (iBit + bitMinH)) / iWidth;
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }
#else
  Int iDMvHorX = ( (acMv[1] - acMv[0]).getHor() << (iBit + bitMinW) ) / iWidth;
  Int iDMvHorY = ( (acMv[1] - acMv[0]).getVer() << (iBit + bitMinH) ) / iWidth;
  Int iDMvVerX = -iDMvHorY;
  Int iDMvVerY =  iDMvHorX;
#endif

  Int iMvScaleHor = acMv[0].getHor() << iBit;
  Int iMvScaleVer = acMv[0].getVer() << iBit;
  Int iMvYHor = iMvScaleHor;
  Int iMvYVer = iMvScaleVer;

  // Calculate Mv for 4x4 Part
  Int iMvScaleTmpHor, iMvScaleTmpVer;
  for ( Int h=0; h<iPartH; h++ )
  {
    for ( Int w=0; w<iPartW; w++ )
    {
      iMvScaleTmpHor = iMvScaleHor + ( iDMvHorX >> 1 ) + ( iDMvVerX >> 1 );
      iMvScaleTmpVer = iMvScaleVer + ( iDMvHorY >> 1 ) + ( iDMvVerY >> 1 );

      // get the MV in hevc define precision
      Int xHevc, yHevc;
      xHevc  = iMvScaleTmpHor >> iBit;
      yHevc  = iMvScaleTmpVer >> iBit;
      TComMv cMvHevc(xHevc, yHevc);
#if JVECT_C0062_AFFINE_SIX_PARAM && !AFFINE_CLIP_REMOVE //ZF
      clipMv(cMvHevc);
#endif

#if TT_TREE_AFFINE // deriveZIdxInCU
      uiPartIdx = deriveZIdxInCU(m_absZIdxInCtu, g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + w + h*getPic()->getNumPartInCtuWidth() ] );
#else
      uiPartIdx = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + w + h*getPic()->getNumPartInCtuWidth() ] - m_absZIdxInCtu;
#endif

      getCUMvField(e)->setMv(cMvHevc, uiPartIdx);

      iMvScaleHor += iDMvHorX;  // switch from x to x+AffineBlockSize, add deltaMvHor
      iMvScaleVer += iDMvHorY;
    }

    iMvYHor += iDMvVerX;        // switch from y to y+AffineBlockSize, add deltaMvVer
    iMvYVer += iDMvVerY;

    iMvScaleHor = iMvYHor;
    iMvScaleVer = iMvYVer;
  }
#endif

  // Set AffineMvField for affine motion compensation LT, RT, LB and RB
  cMv = acMv[1] + acMv[2] - acMv[0];
#if JVECT_C0062_AFFINE_SIX_PARAM && !AFFINE_CLIP_REMOVE
  clipMv(cMv);//ZF
#endif

#if !AFFINE_UNIFY_MC
#if TT_TREE_AFFINE // deriveZIdxInCU
  getCUMvField(e)->setMv( acMv[0], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxLT) );
  getCUMvField(e)->setMv( acMv[1], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxRT) );
  getCUMvField(e)->setMv( acMv[2], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxLB) );
  getCUMvField(e)->setMv( cMv,     deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxRB) );
#else
  getCUMvField(e)->setMv( acMv[0], uiPartIdxLT - m_absZIdxInCtu );
  getCUMvField(e)->setMv( acMv[1], uiPartIdxRT - m_absZIdxInCtu );
  getCUMvField(e)->setMv( acMv[2], uiPartIdxLB - m_absZIdxInCtu );
  getCUMvField(e)->setMv( cMv, uiPartIdxRB - m_absZIdxInCtu );
#endif
#endif
}

Void TComDataCU::setAllAffineMvd( UInt uiAbsPartIdx, UInt uiPuIdx, TComMv acMvd[3], RefPicList e, UInt uiDepth )
{
  Int iWidth  = getWidth(uiAbsPartIdx);
#if !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
  Int iHeight = getHeight(uiAbsPartIdx);
#endif
  Int iPartW = iWidth/getPic()->getMinCUWidth();
#if JVET_C0024_QTBT
  PartSize ePartSize = SIZE_2Nx2N;
#if JVECT_C0062_AFFINE_SIX_PARAM
  Int iPartH = iHeight / getPic()->getMinCUHeight();
#endif
#else
  Int iPartH = iHeight/getPic()->getMinCUHeight();
  PartSize ePartSize = getPartitionSize(uiAbsPartIdx);
#endif

  // Get partIdx of four corner
#if TT_TREE_AFFINE // deriveZIdxInCTU
  UInt uiPartIdxLT = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  UInt uiPartIdxLT = uiAbsPartIdx + m_absZIdxInCtu;
#endif
  UInt uiPartIdxRT = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + iPartW - 1 ];
#if !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
  UInt uiPartIdxLB = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + (iPartH-1)*getPic()->getNumPartInCtuWidth() ];
  UInt uiPartIdxRB = g_auiRasterToZscan[ g_auiZscanToRaster[uiPartIdxLT] + iPartW - 1  + (iPartH-1)*getPic()->getNumPartInCtuWidth() ];
#endif

  // Set other position
#if JVET_C0024_QTBT && !JVECT_C0062_AFFINE_SIX_PARAM
  TComMv cMv = acMvd[0] + acMvd[1];
#else
  TComMv cMv = acMvd[1] + acMvd[2];
#endif
  cMv >>= 1;
  getCUMvField(e)->setAllMvd( cMv, ePartSize, uiAbsPartIdx, uiDepth, uiPuIdx );

  // Set LT, RT, LB and RB
#if !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
  cMv = acMvd[1] + acMvd[2] - acMvd[0];
#endif

#if TT_TREE_AFFINE // deriveZIdxInCU
  getCUMvField(e)->setMvd( acMvd[0], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxLT) );
  getCUMvField(e)->setMvd( acMvd[1], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxRT) );
#if !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
  getCUMvField(e)->setMvd( acMvd[2], deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxLB) );
  getCUMvField(e)->setMvd( cMv,      deriveZIdxInCU(m_absZIdxInCtu, uiPartIdxRB) );
#endif
#else
  getCUMvField(e)->setMvd( acMvd[0], uiPartIdxLT - m_absZIdxInCtu );
  getCUMvField(e)->setMvd( acMvd[1], uiPartIdxRT - m_absZIdxInCtu );
#if !JVET_C0024_QTBT || JVECT_C0062_AFFINE_SIX_PARAM
  getCUMvField(e)->setMvd( acMvd[2], uiPartIdxLB - m_absZIdxInCtu );
  getCUMvField(e)->setMvd( cMv, uiPartIdxRB - m_absZIdxInCtu );
#endif
#endif
}

/** Constructs a list of candidates for Affine AMVP
 * \param uiPartIdx
 * \param uiPartAddr 
 * \param eRefPicList
 * \param iRefIdx
 * \param pInfo
 */
Void TComDataCU::fillAffineMvpCand ( UInt uiPartIdx, UInt uiPartAddr, RefPicList eRefPicList, Int iRefIdx, AffineAMVPInfo* pInfo
#if JVECT_C0062_AFFINE_SIX_PARAM
  , Bool b6Param
#endif
)
{
  TComMv cZeorMv(0, 0);
  pInfo->iN = 0;
  if (iRefIdx < 0)
  {
    return;
  }

  // Get Part Index in LCU
  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
  deriveLeftRightTopIdxGeneral( uiPartAddr, uiPartIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdxGeneral  ( uiPartAddr, uiPartIdx, uiPartIdxLB );

#if AFFINE_AMVP_IMPROVE // add before generate MVP
  Bool bAdded = false;
  TComDataCU* tmpCU = NULL;
  UInt idx;
  tmpCU = getPUBelowLeft(idx, uiPartIdxLB);
  Bool bAddedSmvp = (tmpCU != NULL) && (tmpCU->isAffine(idx));

  if (!bAddedSmvp)
  {
    tmpCU = getPULeft(idx, uiPartIdxLB);
    bAddedSmvp = (tmpCU != NULL) && (tmpCU->isAffine(idx));
  }

  // Left predictor search
  bAdded = xAddAffineMVPCand(pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT);
  if (!bAdded)
  {
    bAdded = xAddAffineMVPCand(pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT);
  }

  if (!bAdded)
  {
    bAdded = xAddAffineMVPCandOrder(pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT);
    if (!bAdded)
    {
      xAddAffineMVPCandOrder(pInfo, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT);
    }
  }

  // Above predictor search
  bAdded = xAddAffineMVPCand(pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT);

  if (!bAdded)
  {
    bAdded = xAddAffineMVPCand(pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE);
  }

  if (!bAdded)
  {
    xAddAffineMVPCand(pInfo, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT);
  }

  if (!bAddedSmvp)
  {
    bAdded = xAddAffineMVPCandOrder(pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT);
    if (!bAdded)
    {
      bAdded = xAddAffineMVPCandOrder(pInfo, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE);
    }

    if (!bAdded)
    {
      xAddAffineMVPCandOrder(pInfo, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT);
    }
  }

  if (pInfo->iN == 2)
  {
    if (pInfo->m_acMvCand[0][0] == pInfo->m_acMvCand[1][0] && pInfo->m_acMvCand[0][1] == pInfo->m_acMvCand[1][1])
    {
      pInfo->iN = 1;
    }
  }

  for (Int i = 0; i < pInfo->iN; i++)
  {
    Int iWidth = getWidth(uiPartAddr);
    Int iHeight = getHeight(uiPartAddr);

#if JVECT_C0062_AFFINE_SIX_PARAM
    if (b6Param)
    {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif

#if !AFFINE_CLIP_REMOVE
      clipMv(pInfo->m_acMvCand[i][0]);
      clipMv(pInfo->m_acMvCand[i][1]);
      clipMv(pInfo->m_acMvCand[i][2]);
#endif
    }
    else
    { //4 param
      if (iHeight / iWidth >= 2)
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif

#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
        Int vx1 = (pInfo->m_acMvCand[i][2].getVer() - pInfo->m_acMvCand[i][0].getVer())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getHor();
        Int vy1 = -(pInfo->m_acMvCand[i][2].getHor() - pInfo->m_acMvCand[i][0].getHor())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][1]);
#endif
      }
      else
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][1]);
#endif
        Int vx2 = -(pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
        Int vy2 = (pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
      }
    }
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
    pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][0]);
    clipMv(pInfo->m_acMvCand[i][1]);
#endif
    Int vx2 = -(pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
    Int vy2 = (pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
    pInfo->m_acMvCand[i][2].set(vx2, vy2);

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE && !JVET_C0024_QTBT
    pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][2]);
#endif
#endif
  }

  if (pInfo->iN == 2)
  {
    return;
  }
#endif


  //-------------------  V0 (START) -------------------//
  AMVPInfo cAMVPInfo0;
  cAMVPInfo0.iN = 0;

  // A->C: Above Left, Above, Left
  xAddMVPCand( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT, true );
  xAddMVPCand( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE     , true );
  xAddMVPCand( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_LEFT      , true );
  if ( cAMVPInfo0.iN < AFFINE_MAX_NUM_V0 )
  {
    xAddMVPCandOrder( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE_LEFT, true );
    if ( cAMVPInfo0.iN < AFFINE_MAX_NUM_V0 )
    {
      xAddMVPCandOrder( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_ABOVE, true );
      if ( cAMVPInfo0.iN < AFFINE_MAX_NUM_V0 )
      {
        xAddMVPCandOrder( &cAMVPInfo0, eRefPicList, iRefIdx, uiPartIdxLT, MD_LEFT, true );
      }
    }
  }

  //-------------------  V1 (START) -------------------//
  AMVPInfo cAMVPInfo1;
  cAMVPInfo1.iN = 0;

  // D->E: Above, Above Right
  xAddMVPCand( &cAMVPInfo1, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE      , true );
  xAddMVPCand( &cAMVPInfo1, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT, true );
  if ( cAMVPInfo1.iN < AFFINE_MAX_NUM_V1 )
  {
    xAddMVPCandOrder( &cAMVPInfo1, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE, true );
    if ( cAMVPInfo1.iN < AFFINE_MAX_NUM_V1 )
    {
      xAddMVPCandOrder( &cAMVPInfo1, eRefPicList, iRefIdx, uiPartIdxRT, MD_ABOVE_RIGHT, true );
    }
  }

  //-------------------  V2 (START) -------------------//
  AMVPInfo cAMVPInfo2;
  cAMVPInfo2.iN = 0;

  // F->G: Left, Below Left
  xAddMVPCand( &cAMVPInfo2, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT      , true );
  xAddMVPCand( &cAMVPInfo2, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT, true );
  if ( cAMVPInfo2.iN < AFFINE_MAX_NUM_V2 )
  {
    xAddMVPCandOrder( &cAMVPInfo2, eRefPicList, iRefIdx, uiPartIdxLB, MD_LEFT, true );
    if ( cAMVPInfo2.iN < AFFINE_MAX_NUM_V2 )
    {
      xAddMVPCandOrder( &cAMVPInfo2, eRefPicList, iRefIdx, uiPartIdxLB, MD_BELOW_LEFT, true );
    }
  }

  // Check Valid Candidates and Sort through DV
  Int   iRecord[AFFINE_MAX_NUM_COMB][3];
  Int   iDV[AFFINE_MAX_NUM_COMB];
  Int   iTempDV;
  Int   iCount = 0;
  for ( Int i=0; i<cAMVPInfo0.iN; i++ )
  {
    for ( Int j=0; j<cAMVPInfo1.iN; j++ )
    {
      for ( Int k=0; k<cAMVPInfo2.iN; k++ )
      {
#if JVECT_C0062_AFFINE_SIX_PARAM
        Bool bValid = isValidAffineCandidate(uiPartAddr, uiPartIdx, cAMVPInfo0.m_acMvCand[i], cAMVPInfo1.m_acMvCand[j], cAMVPInfo2.m_acMvCand[k], iDV[iCount], b6Param);
#else
        Bool bValid = isValidAffineCandidate( uiPartAddr, uiPartIdx, cAMVPInfo0.m_acMvCand[i], cAMVPInfo1.m_acMvCand[j], cAMVPInfo2.m_acMvCand[k], iDV[iCount] );
#endif
        if ( bValid )
        {
          // Sort
          if ( iCount==0 || iDV[iCount]>=iDV[iCount-1] )
          {
            iRecord[iCount][0] = i;
            iRecord[iCount][1] = j;
            iRecord[iCount][2] = k;
          }
          else
          {
            // save last element
            iTempDV = iDV[iCount];
            // find position and move back record
            Int m = 0;
            for ( m=iCount-1; m>=0 && iTempDV<iDV[m]; m-- )
            {
              iDV[m+1] = iDV[m];
              memcpy( iRecord[m+1], iRecord[m], sizeof(Int) * 3 );
            }
            // insert
            iDV[m+1] = iTempDV;
            iRecord[m+1][0] = i;
            iRecord[m+1][1] = j;
            iRecord[m+1][2] = k;
          }
          iCount ++;
        }
      }
    }
  }

#if AFFINE_AMVP_IMPROVE
  Int iCountBefore = pInfo->iN;
  pInfo->iN = min(iCount + pInfo->iN, AMVP_MAX_NUM_CANDS);
#else
  pInfo->iN = min(iCount, AMVP_MAX_NUM_CANDS);
#endif

  Int iWidth = getWidth(uiPartAddr);
  Int iHeight = getHeight(uiPartAddr);

#if AFFINE_AMVP_IMPROVE
  for (Int i = iCountBefore; i < pInfo->iN; i++)
  {
    pInfo->m_acMvCand[i][0] = cAMVPInfo0.m_acMvCand[iRecord[i - iCountBefore][0]];
    pInfo->m_acMvCand[i][1] = cAMVPInfo1.m_acMvCand[iRecord[i - iCountBefore][1]];
    pInfo->m_acMvCand[i][2] = cAMVPInfo2.m_acMvCand[iRecord[i - iCountBefore][2]];

#if JVECT_C0062_AFFINE_SIX_PARAM
    if (b6Param)
    {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
      clipMv(pInfo->m_acMvCand[i][0]);
      clipMv(pInfo->m_acMvCand[i][1]);
      clipMv(pInfo->m_acMvCand[i][2]);
#endif
    }
    else
    { //4 param
      if (iHeight / iWidth >= 2)
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
        Int vx1 = (pInfo->m_acMvCand[i][2].getVer() - pInfo->m_acMvCand[i][0].getVer())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getHor();
        Int vy1 = -(pInfo->m_acMvCand[i][2].getHor() - pInfo->m_acMvCand[i][0].getHor())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][1]);
#endif
      }
      else
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif

#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][1]);
#endif

        Int vx2 = -(pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
        Int vy2 = (pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
      }
    }
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
    pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif

#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][0]);
    clipMv(pInfo->m_acMvCand[i][1]);
#endif
    Int vx2 = -(pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
    Int vy2 = (pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
    pInfo->m_acMvCand[i][2].set(vx2, vy2);

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE && !JVET_C0024_QTBT  //Affine fix
    pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][2]);
#endif
#endif
  }
#else
  for ( Int i=0; i<pInfo->iN; i++ )
  {
    pInfo->m_acMvCand[i][0] = cAMVPInfo0.m_acMvCand[ iRecord[i][0] ];
    pInfo->m_acMvCand[i][1] = cAMVPInfo1.m_acMvCand[ iRecord[i][1] ];
    pInfo->m_acMvCand[i][2] = cAMVPInfo2.m_acMvCand[ iRecord[i][2] ];

#if JVECT_C0062_AFFINE_SIX_PARAM
    if (b6Param)
    {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
      pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
      pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
      clipMv(pInfo->m_acMvCand[i][0]);
      clipMv(pInfo->m_acMvCand[i][1]);
      clipMv(pInfo->m_acMvCand[i][2]);
#endif
    }
    else
    { //4 param
      if (iHeight / iWidth >= 2)
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
        Int vx1 = (pInfo->m_acMvCand[i][2].getVer() - pInfo->m_acMvCand[i][0].getVer())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getHor();
        Int vy1 = -(pInfo->m_acMvCand[i][2].getHor() - pInfo->m_acMvCand[i][0].getHor())*iWidth / iHeight + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][1]);
#endif
      }
      else
      {
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
        pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
        pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][0]);
        clipMv(pInfo->m_acMvCand[i][1]);
#endif
        Int vx2 = -(pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
        Int vy2 = (pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor()) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
        pInfo->m_acMvCand[i][2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
        clipMv(pInfo->m_acMvCand[i][2]);
#endif
      }
    }
#else
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
    pInfo->m_acMvCand[i][0].roundMV2SignalPrecision();
    pInfo->m_acMvCand[i][1].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][0]);
    clipMv(pInfo->m_acMvCand[i][1]);
#endif
    Int vx2 =  - ( pInfo->m_acMvCand[i][1].getVer() - pInfo->m_acMvCand[i][0].getVer() ) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getHor();
    Int vy2 =    ( pInfo->m_acMvCand[i][1].getHor() - pInfo->m_acMvCand[i][0].getHor() ) * iHeight / iWidth + pInfo->m_acMvCand[i][0].getVer();
    pInfo->m_acMvCand[i][2].set( vx2, vy2 );

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE && !JVET_C0024_QTBT  //Affine fix
    pInfo->m_acMvCand[i][2].roundMV2SignalPrecision();
#endif
#if !AFFINE_CLIP_REMOVE
    clipMv(pInfo->m_acMvCand[i][2]);
#endif
#endif
  }
#endif

  if ( pInfo->iN < 2 )
  {
    AMVPInfo* pAMVPInfo = getCUMvField( eRefPicList )->getAMVPInfo();
    fillMvpCand( uiPartIdx, uiPartAddr, eRefPicList, iRefIdx, pAMVPInfo );

    Int iAdd = pAMVPInfo->iN - pInfo->iN;
    for ( Int i=0; i<iAdd; i++ )
    {
#if !AFFINE_CLIP_REMOVE
      clipMv(pAMVPInfo->m_acMvCand[i]);
#endif
      pInfo->m_acMvCand[pInfo->iN][0] = pAMVPInfo->m_acMvCand[i];
      pInfo->m_acMvCand[pInfo->iN][1] = pAMVPInfo->m_acMvCand[i];
      pInfo->m_acMvCand[pInfo->iN][2] = pAMVPInfo->m_acMvCand[i];
      pInfo->iN ++;
    }
  }

  return ;
}

Bool TComDataCU::isValidAffineCandidate( UInt uiAbsPartIdx, UInt uiPuIdx, TComMv cMv0, TComMv cMv1, TComMv cMv2, Int& riDV
#if JVECT_C0062_AFFINE_SIX_PARAM
  , Bool b6Param
#endif
)
{
  TComMv zeroMv(0, 0);
  TComMv deltaHor = cMv1 - cMv0;
  TComMv deltaVer = cMv2 - cMv0;

  // same motion vector, translation model
#if JVECT_C0062_AFFINE_SIX_PARAM
  if (b6Param)
  {
    if (deltaHor == zeroMv && deltaVer == zeroMv)
      return false;
  }
  else
  {
    if (deltaHor == zeroMv)
      return false;
  }
#else
  if ( deltaHor == zeroMv )
    return false;
#endif

  // S/8, but the Mv is 4 precision, so change to S/2
  Int width = getWidth(uiAbsPartIdx);
  Int height = getHeight(uiAbsPartIdx);
  Int iDiffHor = width>>1;
  Int iDiffVer = height>>1;

  if ( deltaHor.getAbsHor() > iDiffHor || deltaHor.getAbsVer() > iDiffVer || deltaVer.getAbsHor() > iDiffHor || deltaVer.getAbsVer() > iDiffVer )
    return false;

  // Calculate DV
  riDV = abs( deltaHor.getHor() * height - deltaVer.getVer() * width ) + abs( deltaHor.getVer() * height + deltaVer.getHor() * width );
  return true;
}

Void TComDataCU::getAffineMergeCandidates( UInt uiAbsPartIdx, UInt uiPUIdx, TComMvField (*pcMvFieldNeighbours)[3], UChar* puhInterDirNeighbours, Int& numValidMergeCand,
#if JVECT_C0062_AFFINE_SIX_PARAM
  UInt& neighborAffineParam,
#endif
  Int mrgCandIdx )
{
  TComMvField affineMvField[2][3];
  UInt uiAbsPartAddr = m_absZIdxInCtu + uiAbsPartIdx;

  TComDataCU *pcTempCU = NULL;
  Bool bAvailable = false;
  UInt uiPartIdx = 0;

  for ( Int mvNum = 0; mvNum < 3; mvNum++ )
  {
    pcMvFieldNeighbours[0][mvNum].setMvField( TComMv(0,0), -1 );
    pcMvFieldNeighbours[1][mvNum].setMvField( TComMv(0,0), -1 );
  }

  // Get Part Index in LCU
  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
  deriveLeftRightTopIdxGeneral( uiAbsPartIdx, uiPUIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdxGeneral  ( uiAbsPartIdx, uiPUIdx, uiPartIdxLB );

#if AFFINE_MERGE_IMPROVE
  // left
  pcTempCU = getPULeft(uiPartIdx, uiPartIdxLB);
  bAvailable = pcTempCU && pcTempCU->isAffine(uiPartIdx);

  // above
  if (bAvailable)
  {
    Int size = pcTempCU->getWidth(uiPartIdx) + pcTempCU->getHeight(uiPartIdx);

    UInt uiPartIdxOther = 0;
    TComDataCU *pcOtherCU = getPUAbove(uiPartIdxOther, uiPartIdxRT);
    Bool bAvailableOther = pcOtherCU && pcOtherCU->isAffine(uiPartIdxOther);

    if (bAvailableOther)
    {
      Int  sizeOther = pcOtherCU->getWidth(uiPartIdxOther) + pcOtherCU->getHeight(uiPartIdxOther);
      if (sizeOther > size)
      {
        pcTempCU = pcOtherCU;
        uiPartIdx = uiPartIdxOther;
      }
    }
  }
  else
  {
    pcTempCU = getPUAbove(uiPartIdx, uiPartIdxRT);
    bAvailable = pcTempCU && pcTempCU->isAffine(uiPartIdx);
  }

  // aboveRight
  if (bAvailable)
  {
    Int size = pcTempCU->getWidth(uiPartIdx) + pcTempCU->getHeight(uiPartIdx);

    UInt uiPartIdxOther = 0;
    TComDataCU *pcOtherCU = getPUAboveRight(uiPartIdxOther, uiPartIdxRT);
    Bool bAvailableOther = pcOtherCU && pcOtherCU->isAffine(uiPartIdxOther);

    if (bAvailableOther)
    {
      Int  sizeOther = pcOtherCU->getWidth(uiPartIdxOther) + pcOtherCU->getHeight(uiPartIdxOther);
      if (sizeOther > size)
      {
        pcTempCU = pcOtherCU;
        uiPartIdx = uiPartIdxOther;
      }
    }
  }
  else
  {
    pcTempCU = getPUAboveRight(uiPartIdx, uiPartIdxRT);
    bAvailable = pcTempCU && pcTempCU->isAffine(uiPartIdx);
  }

  // leftBottom
  if (bAvailable)
  {
    Int size = pcTempCU->getWidth(uiPartIdx) + pcTempCU->getHeight(uiPartIdx);

    UInt uiPartIdxOther = 0;
    TComDataCU *pcOtherCU = getPUBelowLeft(uiPartIdxOther, uiPartIdxLB);
    Bool bAvailableOther = pcOtherCU && pcOtherCU->isAffine(uiPartIdxOther);

    if (bAvailableOther)
    {
      Int  sizeOther = pcOtherCU->getWidth(uiPartIdxOther) + pcOtherCU->getHeight(uiPartIdxOther);
      if (sizeOther > size)
      {
        pcTempCU = pcOtherCU;
        uiPartIdx = uiPartIdxOther;
      }
    }
  }
  else
  {
    pcTempCU = getPUBelowLeft(uiPartIdx, uiPartIdxLB);
    bAvailable = pcTempCU && pcTempCU->isAffine(uiPartIdx);
  }

  // aboveLeft
  if (bAvailable)
  {
    Int size = pcTempCU->getWidth(uiPartIdx) + pcTempCU->getHeight(uiPartIdx);

    UInt uiPartIdxOther = 0;
    TComDataCU *pcOtherCU = getPUAboveLeft(uiPartIdxOther, uiAbsPartAddr);
    Bool bAvailableOther = pcOtherCU && pcOtherCU->isAffine(uiPartIdxOther);

    if (bAvailableOther)
    {
      Int  sizeOther = pcOtherCU->getWidth(uiPartIdxOther) + pcOtherCU->getHeight(uiPartIdxOther);
      if (sizeOther > size)
      {
        pcTempCU = pcOtherCU;
        uiPartIdx = uiPartIdxOther;
      }
    }
  }
  else
  {
    pcTempCU = getPUAboveLeft(uiPartIdx, uiAbsPartAddr);
    bAvailable = pcTempCU && pcTempCU->isAffine(uiPartIdx);
  }
#else
  // left
  pcTempCU = getPULeft( uiPartIdx, uiPartIdxLB );
  bAvailable = pcTempCU && pcTempCU->isAffine( uiPartIdx );

  // above
  if ( !bAvailable )
  {
    pcTempCU = getPUAbove( uiPartIdx, uiPartIdxRT );
    bAvailable = pcTempCU && pcTempCU->isAffine( uiPartIdx );
  }

  // aboveRight
  if ( !bAvailable )
  {
    pcTempCU = getPUAboveRight( uiPartIdx, uiPartIdxRT );
    bAvailable = pcTempCU && pcTempCU->isAffine( uiPartIdx );
  }

  // leftBottom
  if ( !bAvailable )
  {
    pcTempCU = getPUBelowLeft( uiPartIdx, uiPartIdxLB );
    bAvailable = pcTempCU && pcTempCU->isAffine( uiPartIdx );
  }

  // aboveLeft
  if ( !bAvailable )
  {
    pcTempCU = getPUAboveLeft( uiPartIdx, uiAbsPartAddr );
    bAvailable = pcTempCU && pcTempCU->isAffine( uiPartIdx );
  }
#endif

  if ( !bAvailable )
  {
    numValidMergeCand = -1;
    return;
  }
  else
  {
    assert( pcTempCU );
    numValidMergeCand = 1;
  }

  Int width = pcTempCU->getWidth(uiPartIdx);
  Int height = pcTempCU->getHeight(uiPartIdx);
#if !JVET_C0024_QTBT
  Int depth  = pcTempCU->getDepth(uiPartIdx);
  Int numPart = pcTempCU->getPic()->getNumPartitionsInCtu() >> ( depth << 1 );
#endif
  Int iPartW = width / getPic()->getMinCUWidth();
  Int iPartH = height / getPic()->getMinCUHeight();

#if JVET_C0024_QTBT
#if TT_TREE_AFFINE //important fix
  Int aboveLeftIdx = pcTempCU->getCUAbsZIdxInCtuComp(uiPartIdx);
#else
  Int iBlkY = g_auiZscanToRaster[uiPartIdx]/(getPic()->getNumPartInCtuWidth()*iPartH);
  Int iBlkX = (g_auiZscanToRaster[uiPartIdx]%getPic()->getNumPartInCtuWidth())/iPartW;
  Int aboveLeftIdx = g_auiRasterToZscan[iBlkY * getPic()->getNumPartInCtuWidth() * iPartH + iBlkX * iPartW];
#endif
#else
  Int aboveLeftIdx = uiPartIdx - ( uiPartIdx % numPart );
#endif
  Int aboveRightIdx = g_auiRasterToZscan[ g_auiZscanToRaster[aboveLeftIdx] + iPartW - 1 ];
  Int bottomLeftIdx = g_auiRasterToZscan[ g_auiZscanToRaster[aboveLeftIdx] + ( iPartH - 1 ) * getPic()->getNumPartInCtuWidth() ];

  // list0
  affineMvField[0][0].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_0)->getMv( aboveLeftIdx ),  pcTempCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx( aboveLeftIdx ) );
  affineMvField[0][1].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_0)->getMv( aboveRightIdx ), pcTempCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx( aboveRightIdx ) );
  affineMvField[0][2].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_0)->getMv( bottomLeftIdx ), pcTempCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx( bottomLeftIdx ) );

  // list1
  affineMvField[1][0].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_1)->getMv( aboveLeftIdx ),  pcTempCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx( aboveLeftIdx ) );
  affineMvField[1][1].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_1)->getMv( aboveRightIdx ), pcTempCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx( aboveRightIdx ) );
  affineMvField[1][2].setMvField( pcTempCU->getCUMvField(REF_PIC_LIST_1)->getMv( bottomLeftIdx ), pcTempCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx( bottomLeftIdx ) );

  Int pixelOrgX    = g_auiRasterToPelX[ g_auiZscanToRaster[ aboveLeftIdx ] ] + pcTempCU->getCUPelX();
  Int pixelCurrX   = getCUPelX();
  Int pixelOrgY    = g_auiRasterToPelY[ g_auiZscanToRaster[ aboveLeftIdx ] ] + pcTempCU->getCUPelY();
  Int pixelCurrY   = getCUPelY();

  Int vx0 = Int( affineMvField[0][0].getHor() + 1.0 * ( affineMvField[0][2].getHor() - affineMvField[0][0].getHor() ) * ( pixelCurrY - pixelOrgY ) / height
    + 1.0 * ( affineMvField[0][1].getHor() - affineMvField[0][0].getHor() ) * ( pixelCurrX - pixelOrgX ) / width );
  Int vy0 = Int( affineMvField[0][0].getVer() + 1.0 * ( affineMvField[0][2].getVer() - affineMvField[0][0].getVer() ) * ( pixelCurrY - pixelOrgY ) / height
    + 1.0 * ( affineMvField[0][1].getVer() - affineMvField[0][0].getVer() ) * ( pixelCurrX - pixelOrgX ) / width );
  pcMvFieldNeighbours[0][0].setMvField( TComMv(vx0, vy0), affineMvField[0][0].getRefIdx() );
  if ( getSlice()->isInterB() )
  {
    vx0 = Int( affineMvField[1][0].getHor() + 1.0 * ( affineMvField[1][2].getHor() - affineMvField[1][0].getHor() ) * ( pixelCurrY - pixelOrgY ) / height 
      + 1.0 * ( affineMvField[1][1].getHor() - affineMvField[1][0].getHor() ) * ( pixelCurrX - pixelOrgX ) / width );
    vy0 = Int( affineMvField[1][0].getVer() + 1.0 * ( affineMvField[1][2].getVer() - affineMvField[1][0].getVer() ) * ( pixelCurrY - pixelOrgY ) / height 
      + 1.0 * ( affineMvField[1][1].getVer() - affineMvField[1][0].getVer() ) * ( pixelCurrX - pixelOrgX ) / width );
    pcMvFieldNeighbours[1][0].setMvField( TComMv(vx0, vy0), affineMvField[1][0].getRefIdx() );
  }

  Int vx1 = Int( 1.0 * ( affineMvField[0][1].getHor() - affineMvField[0][0].getHor() ) * getWidth(0) / width + pcMvFieldNeighbours[0][0].getHor() );
  Int vy1 = Int( 1.0 * ( affineMvField[0][1].getVer() - affineMvField[0][0].getVer() ) * getWidth(0) / width + pcMvFieldNeighbours[0][0].getVer() );
  pcMvFieldNeighbours[0][1].setMvField( TComMv(vx1, vy1), affineMvField[0][0].getRefIdx() );
  if ( getSlice()->isInterB() )
  {
    vx1 = Int( 1.0 * ( affineMvField[1][1].getHor() - affineMvField[1][0].getHor() ) * getWidth(0) / width + pcMvFieldNeighbours[1][0].getHor() );
    vy1 = Int( 1.0 * ( affineMvField[1][1].getVer() - affineMvField[1][0].getVer() ) * getWidth(0) / width + pcMvFieldNeighbours[1][0].getVer() );
    pcMvFieldNeighbours[1][1].setMvField( TComMv(vx1, vy1), affineMvField[1][0].getRefIdx() );
  }

  Int vx2 = Int( 1.0 * ( affineMvField[0][2].getHor() - affineMvField[0][0].getHor() ) * getHeight(0) / height + pcMvFieldNeighbours[0][0].getHor() );
  Int vy2 = Int( 1.0 * ( affineMvField[0][2].getVer() - affineMvField[0][0].getVer() ) * getHeight(0) / height + pcMvFieldNeighbours[0][0].getVer() );
  pcMvFieldNeighbours[0][2].setMvField( TComMv(vx2, vy2), affineMvField[0][0].getRefIdx() );
  if ( getSlice()->isInterB() )
  {
    vx2 = Int( 1.0 * ( affineMvField[1][2].getHor() - affineMvField[1][0].getHor() ) * getHeight(0) / height + pcMvFieldNeighbours[1][0].getHor() );
    vy2 = Int( 1.0 * ( affineMvField[1][2].getVer() - affineMvField[1][0].getVer() ) * getHeight(0) / height + pcMvFieldNeighbours[1][0].getVer() );
    pcMvFieldNeighbours[1][2].setMvField( TComMv(vx2, vy2), affineMvField[1][0].getRefIdx() );
  }

#if JVECT_C0062_AFFINE_SIX_PARAM
  neighborAffineParam = 1;
#endif
  puhInterDirNeighbours[0] = pcTempCU->getInterDir( uiPartIdx );
  return;
}

Bool TComDataCU::isAffineMrgFlagCoded( UInt uiAbsPartIdx, UInt uiPUIdx )
{
#if JVET_C0024_QTBT
  if (getWidth(uiAbsPartIdx)*getHeight(uiAbsPartIdx)< 64)
  {
    return false;
  }
#endif
#if TT_TREE_AFFINE // deriveZIdxInCTU
  UInt uiAbsPartAddr = deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx);
#else
  UInt uiAbsPartAddr = m_absZIdxInCtu + uiAbsPartIdx;
#endif

  // Get Part Index in LCU
  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
  deriveLeftRightTopIdxGeneral( uiAbsPartIdx, uiPUIdx, uiPartIdxLT, uiPartIdxRT );
  deriveLeftBottomIdxGeneral  ( uiAbsPartIdx, uiPUIdx, uiPartIdxLB );

  // left
  UInt uiLeftPartIdx = 0;
  TComDataCU* pcCULeft = 0;
  pcCULeft = getPULeft( uiLeftPartIdx, uiPartIdxLB );
  Bool isLeftAvailable = pcCULeft && pcCULeft->isAffine( uiLeftPartIdx );

  // above
  UInt uiAbovePartIdx = 0;
  TComDataCU* pcCUAbove = 0;
  pcCUAbove = getPUAbove( uiAbovePartIdx, uiPartIdxRT );
  Bool isAboveAvailable = pcCUAbove && pcCUAbove->isAffine( uiAbovePartIdx );

  // aboveRight
  UInt uiAboveRightPartIdx = 0;
  TComDataCU* pcCUAboveRight = 0;
  pcCUAboveRight = getPUAboveRight( uiAboveRightPartIdx, uiPartIdxRT );
  Bool isAboveRightAvailable = pcCUAboveRight && pcCUAboveRight->isAffine( uiAboveRightPartIdx );

  // leftBottom
  UInt uiBottomLeftPartIdx = 0;
  TComDataCU* pcCUBottomLeft = 0;
  pcCUBottomLeft = getPUBelowLeft( uiBottomLeftPartIdx, uiPartIdxLB );
  Bool isBottomLeftAvailable = pcCUBottomLeft && pcCUBottomLeft->isAffine( uiBottomLeftPartIdx );

  //aboveLeft
  UInt uiAboveLeftPartIdx = 0;
  TComDataCU* pcCUAboveLeft = 0;
  pcCUAboveLeft = getPUAboveLeft( uiAboveLeftPartIdx, uiAbsPartAddr );
  Bool isAboveLeftAvailable = pcCUAboveLeft && pcCUAboveLeft->isAffine( uiAboveLeftPartIdx );

  return ( isLeftAvailable || isAboveAvailable || isAboveRightAvailable || isBottomLeftAvailable || isAboveLeftAvailable );
}
#endif // end AFFINE

#if AFFINE_AMVP_IMPROVE
Bool TComDataCU::xAddAffineMVPCand(AffineAMVPInfo* pInfo, RefPicList eRefPicList, Int iRefIdx, UInt uiPartUnitIdx, MVP_DIR eDir)
{
  TComDataCU* pcTmpCU = NULL;
  UInt uiIdx;
  switch (eDir)
  {
  case MD_LEFT:
  {
    pcTmpCU = getPULeft(uiIdx, uiPartUnitIdx);
    break;
    }
  case MD_ABOVE:
  {
    pcTmpCU = getPUAbove(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_ABOVE_RIGHT:
  {
    pcTmpCU = getPUAboveRight(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_BELOW_LEFT:
  {
    pcTmpCU = getPUBelowLeft(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_ABOVE_LEFT:
  {
    pcTmpCU = getPUAboveLeft(uiIdx, uiPartUnitIdx);
    break;
  }
  default:
  {
    break;
  }
  }

  if (pcTmpCU == NULL || pcTmpCU->getAffineFlag(uiIdx) == false)
  {
    return false;
  }

  if (pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) >= 0 && m_pcSlice->getRefPic(eRefPicList, iRefIdx)->getPOC() == pcTmpCU->getSlice()->getRefPOC(eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx)))
  {
    TComMv cMvPred[3];
    xDerivePredAffineMv(pcTmpCU, uiIdx, eRefPicList, cMvPred);
    pInfo->m_acMvCand[pInfo->iN][0] = cMvPred[0];
    pInfo->m_acMvCand[pInfo->iN][1] = cMvPred[1];
    pInfo->m_acMvCand[pInfo->iN][2] = cMvPred[2];
    pInfo->iN++;
    return true;
  }

  RefPicList eRefPicList2nd = REF_PIC_LIST_0;
  if (eRefPicList == REF_PIC_LIST_0)
  {
    eRefPicList2nd = REF_PIC_LIST_1;
  }
  else if (eRefPicList == REF_PIC_LIST_1)
  {
    eRefPicList2nd = REF_PIC_LIST_0;
  }

  Int iCurrRefPOC = m_pcSlice->getRefPic(eRefPicList, iRefIdx)->getPOC();
  Int iNeibRefPOC;

  if (pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) >= 0)
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC(eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx));
    if (iNeibRefPOC == iCurrRefPOC) // Same Reference Frame But Diff List//
    {
      TComMv cMvPred[3];
      xDerivePredAffineMv(pcTmpCU, uiIdx, eRefPicList2nd, cMvPred);
      pInfo->m_acMvCand[pInfo->iN][0] = cMvPred[0];
      pInfo->m_acMvCand[pInfo->iN][1] = cMvPred[1];
      pInfo->m_acMvCand[pInfo->iN][2] = cMvPred[2];
      pInfo->iN++;
      return true;
    }
  }
  return false;
}

/**
* \param pInfo
* \param eRefPicList
* \param iRefIdx
* \param uiPartUnitIdx
* \param eDir
* \returns Bool
*/
Bool TComDataCU::xAddAffineMVPCandOrder(AffineAMVPInfo* pInfo, RefPicList eRefPicList, Int iRefIdx, UInt uiPartUnitIdx, MVP_DIR eDir)
{
  TComDataCU* pcTmpCU = NULL;
  UInt uiIdx;
  switch (eDir)
  {
  case MD_LEFT:
  {
    pcTmpCU = getPULeft(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_ABOVE:
  {
    pcTmpCU = getPUAbove(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_ABOVE_RIGHT:
  {
    pcTmpCU = getPUAboveRight(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_BELOW_LEFT:
  {
    pcTmpCU = getPUBelowLeft(uiIdx, uiPartUnitIdx);
    break;
  }
  case MD_ABOVE_LEFT:
  {
    pcTmpCU = getPUAboveLeft(uiIdx, uiPartUnitIdx);
    break;
  }
  default:
  {
    break;
  }
  }

  if (pcTmpCU == NULL || pcTmpCU->getAffineFlag(uiIdx) == false)
  {
    return false;
  }

  RefPicList eRefPicList2nd = REF_PIC_LIST_0;
  if (eRefPicList == REF_PIC_LIST_0)
  {
    eRefPicList2nd = REF_PIC_LIST_1;
  }
  else if (eRefPicList == REF_PIC_LIST_1)
  {
    eRefPicList2nd = REF_PIC_LIST_0;
  }

  Int iCurrPOC = m_pcSlice->getPOC();
  Int iCurrRefPOC = m_pcSlice->getRefPic(eRefPicList, iRefIdx)->getPOC();
  Int iNeibPOC = iCurrPOC;
  Int iNeibRefPOC;
  Bool bIsCurrRefLongTerm = m_pcSlice->getRefPic(eRefPicList, iRefIdx)->getIsLongTerm();
  Bool bIsNeibRefLongTerm = false;

  //---------------  V1 (END) ------------------//
  if (pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx) >= 0)
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC(eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx));
    TComMv cMvPred[3], rcMv[3];
    xDerivePredAffineMv(pcTmpCU, uiIdx, eRefPicList, cMvPred);

    bIsNeibRefLongTerm = pcTmpCU->getSlice()->getRefPic(eRefPicList, pcTmpCU->getCUMvField(eRefPicList)->getRefIdx(uiIdx))->getIsLongTerm();
    if (bIsCurrRefLongTerm == bIsNeibRefLongTerm)
    {
      if (bIsCurrRefLongTerm || bIsNeibRefLongTerm)
      {
        rcMv[0] = cMvPred[0];
        rcMv[1] = cMvPred[1];
        rcMv[2] = cMvPred[2];
      }
      else
      {
        Int iScale = xGetDistScaleFactor(iCurrPOC, iCurrRefPOC, iNeibPOC, iNeibRefPOC);
        if (iScale == 4096)
        {
          rcMv[0] = cMvPred[0];
          rcMv[1] = cMvPred[1];
          rcMv[2] = cMvPred[2];
        }
        else
        {
          rcMv[0] = cMvPred[0].scaleMv(iScale);
          rcMv[1] = cMvPred[1].scaleMv(iScale);
          rcMv[2] = cMvPred[2].scaleMv(iScale);
        }
      }

      pInfo->m_acMvCand[pInfo->iN][0] = rcMv[0];
      pInfo->m_acMvCand[pInfo->iN][1] = rcMv[1];
      pInfo->m_acMvCand[pInfo->iN][2] = rcMv[2];
      pInfo->iN++;
      return true;
    }
  }
  //---------------------- V2(END) --------------------//
  if (pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx) >= 0)
  {
    iNeibRefPOC = pcTmpCU->getSlice()->getRefPOC(eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx));
    TComMv cMvPred[3], rcMv[3];
    xDerivePredAffineMv(pcTmpCU, uiIdx, eRefPicList, cMvPred);

    bIsNeibRefLongTerm = pcTmpCU->getSlice()->getRefPic(eRefPicList2nd, pcTmpCU->getCUMvField(eRefPicList2nd)->getRefIdx(uiIdx))->getIsLongTerm();
    if (bIsCurrRefLongTerm == bIsNeibRefLongTerm)
    {
      if (bIsCurrRefLongTerm || bIsNeibRefLongTerm)
      {
        rcMv[0] = cMvPred[0];
        rcMv[1] = cMvPred[1];
        rcMv[2] = cMvPred[2];
      }
      else
      {
        Int iScale = xGetDistScaleFactor(iCurrPOC, iCurrRefPOC, iNeibPOC, iNeibRefPOC);
        if (iScale == 4096)
        {
          rcMv[0] = cMvPred[0];
          rcMv[1] = cMvPred[1];
          rcMv[2] = cMvPred[2];
        }
        else
        {
          rcMv[0] = cMvPred[0].scaleMv(iScale);
          rcMv[1] = cMvPred[1].scaleMv(iScale);
          rcMv[2] = cMvPred[2].scaleMv(iScale);
        }
      }

      pInfo->m_acMvCand[pInfo->iN][0] = rcMv[0];
      pInfo->m_acMvCand[pInfo->iN][1] = rcMv[1];
      pInfo->m_acMvCand[pInfo->iN][2] = rcMv[2];
      pInfo->iN++;
      return true;
    }
  }
  //---------------------- V3(END) --------------------//
  return false;
}

Void TComDataCU::xDerivePredAffineMv(TComDataCU* pcTempCU, UInt uiPartAddr, RefPicList eRefPicList, TComMv rcMv[3])
{
  Int width = pcTempCU->getWidth(uiPartAddr);
  Int height = pcTempCU->getHeight(uiPartAddr);

#if !JVET_C0024_QTBT
  Int depth = pcTempCU->getDepth(uiPartIdx);
  Int numPart = pcTempCU->getPic()->getNumPartitionsInCtu() >> (depth << 1);
#endif
  Int iPartW = width / getPic()->getMinCUWidth();
  Int iPartH = height / getPic()->getMinCUHeight();

#if JVET_C0024_QTBT
#if TT_TREE_AFFINE //important fix
  Int aboveLeftIdx = pcTempCU->getCUAbsZIdxInCtuComp(uiPartAddr);
#else
  Int iBlkY = g_auiZscanToRaster[uiPartAddr] / (getPic()->getNumPartInCtuWidth()*iPartH);
  Int iBlkX = (g_auiZscanToRaster[uiPartAddr] % getPic()->getNumPartInCtuWidth()) / iPartW;
  Int aboveLeftIdx = g_auiRasterToZscan[iBlkY * getPic()->getNumPartInCtuWidth() * iPartH + iBlkX * iPartW];
#endif
#else
  Int aboveLeftIdx = uiPartAddr - (uiPartAddr % numPart);
#endif
  Int aboveRightIdx = g_auiRasterToZscan[g_auiZscanToRaster[aboveLeftIdx] + iPartW - 1];
  Int bottomLeftIdx = g_auiRasterToZscan[g_auiZscanToRaster[aboveLeftIdx] + (iPartH - 1) * getPic()->getNumPartInCtuWidth()];

  TComMv acMvNeigh[3];
  acMvNeigh[0] = pcTempCU->getCUMvField(eRefPicList)->getMv(aboveLeftIdx);
  acMvNeigh[1] = pcTempCU->getCUMvField(eRefPicList)->getMv(aboveRightIdx);
  acMvNeigh[2] = pcTempCU->getCUMvField(eRefPicList)->getMv(bottomLeftIdx);

  Int pixelOrgX = g_auiRasterToPelX[g_auiZscanToRaster[aboveLeftIdx]] + pcTempCU->getCUPelX();  Int pixelCurrX = getCUPelX();
  Int pixelOrgY = g_auiRasterToPelY[g_auiZscanToRaster[aboveLeftIdx]] + pcTempCU->getCUPelY();  Int pixelCurrY = getCUPelY();

  Int vx0 = Int(acMvNeigh[0].getHor() + 1.0 * (acMvNeigh[2].getHor() - acMvNeigh[0].getHor()) * (pixelCurrY - pixelOrgY) / height
    + 1.0 * (acMvNeigh[1].getHor() - acMvNeigh[0].getHor()) * (pixelCurrX - pixelOrgX) / width);
  Int vy0 = Int(acMvNeigh[0].getVer() + 1.0 * (acMvNeigh[2].getVer() - acMvNeigh[0].getVer()) * (pixelCurrY - pixelOrgY) / height
    + 1.0 * (acMvNeigh[1].getVer() - acMvNeigh[0].getVer()) * (pixelCurrX - pixelOrgX) / width);
  rcMv[0].set(vx0, vy0);

  Int vx1 = Int(1.0 * (acMvNeigh[1].getHor() - acMvNeigh[0].getHor()) * getWidth(0) / width + rcMv[0].getHor());
  Int vy1 = Int(1.0 * (acMvNeigh[1].getVer() - acMvNeigh[0].getVer()) * getWidth(0) / width + rcMv[0].getVer());
  rcMv[1].set(vx1, vy1);

  Int vx2 = -(rcMv[1].getVer() - rcMv[0].getVer()) + rcMv[0].getHor();
  Int vy2 = (rcMv[1].getHor() - rcMv[0].getHor()) + rcMv[0].getVer();
  rcMv[2].set(vx2, vy2);
}

#endif

#if MRG_CONSIDER_COMPLEX
UInt TComDataCU::getCtxMergeComplexFlag(UInt uiAbsPartIdx)
{
  UInt uiCtx = 0;

  TComDataCU* pcTempCU;
  UInt        uiTempPartIdx;

  // Get BCBP of left PU
#if TT_TREE_AFFINE //MRG_CONSIDER_COMPLEX
  pcTempCU = getPULeft(uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx));
#else
  pcTempCU = getPULeft(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
#endif
  uiCtx = (pcTempCU) ? pcTempCU->getMergeComplexFlag(uiTempPartIdx) : 0;

  // Get BCBP of above PU
#if TT_TREE_AFFINE //MRG_CONSIDER_COMPLEX
  pcTempCU = getPUAbove(uiTempPartIdx, deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx));
#else
  pcTempCU = getPUAbove(uiTempPartIdx, m_absZIdxInCtu + uiAbsPartIdx);
#endif
  uiCtx += (pcTempCU) ? pcTempCU->getMergeComplexFlag(uiTempPartIdx) : 0;

  return uiCtx;
}

template<typename T>
#if JVET_C0024_QTBT
Void TComDataCU::setSubPartT(T uiParameter, T* puhBaseLCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight)
{
  UInt uiRaster = g_auiZscanToRaster[m_absZIdxInCtu + uiAbsPartIdx];
  UInt uiShort, uiLong;
  UInt uiStride;
  if (uiHeight > uiWidth)
  {
    uiShort = uiWidth;
    uiLong = uiHeight;
    uiStride = m_pcPic->getNumPartInCtuWidth();
  }
  else
  {
    uiShort = uiHeight;
    uiLong = uiWidth;
    uiStride = 1;
  }

  UInt uiDepth = g_aucConvertToBit[m_pcSlice->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
  UInt uiCurrPartNumb = m_pcPic->getNumPartitionsInCtu() >> (uiDepth << 1);
  UInt uiNumPartInShort = m_pcPic->getNumPartInCtuWidth() >> uiDepth;

  for (UInt i = 0; i < uiLong; i += uiShort)
  {
    for (UInt ui = 0; ui < uiCurrPartNumb; ui++)
      puhBaseLCU[g_auiRasterToZscan[uiRaster] - m_absZIdxInCtu + ui] = uiParameter;
    uiRaster += uiNumPartInShort * uiStride;
  }
}

Void TComDataCU::setMergeComplexFlagSubParts(Bool bMergeComplexFlag, UInt uiAbsPartIdx)
{
  setSubPart(bMergeComplexFlag, m_pbMergeComplexFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}

Void TComDataCU::setNodeMotionInfoSubParts(NodeInfo cNodeMotionInfo, UInt uiAbsPartIdx)
{
  setSubPartT<NodeInfo>(cNodeMotionInfo, m_pNodeMotionInfo, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}

Void TComDataCU::setNumMergeComplexCandSubParts(UChar NumMergeComplexCand, UInt uiAbsPartIdx)
{
  setSubPartT<UChar>(NumMergeComplexCand, m_pNumMergeComplexCand, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}

Void TComDataCU::setMergeComplexIndexSubParts(UInt uiMergeComplexIndex, UInt uiAbsPartIdx)
{
  setSubPart<UChar>(uiMergeComplexIndex, m_puhMergeComplexIndex, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
}
#else
Void TComDataCU::setSubPartT(T uiParameter, T* puhBaseLCU, UInt uiCUAddr, UInt uiCUDepth, UInt uiPUIdx)
{
  UInt uiCurrPartNumQ = (m_pcPic->getNumPartInCU() >> (2 * uiCUDepth)) >> 2;
  switch (m_pePartSize[uiCUAddr])
  {
  case SIZE_2Nx2N:
    for (UInt ui = 0; ui < 4 * uiCurrPartNumQ; ui++)
      puhBaseLCU[uiCUAddr + ui] = uiParameter;

    break;
  case SIZE_2NxN:
    for (UInt ui = 0; ui < 2 * uiCurrPartNumQ; ui++)
      puhBaseLCU[uiCUAddr + ui] = uiParameter;
    break;
  case SIZE_Nx2N:
    for (UInt ui = 0; ui < uiCurrPartNumQ; ui++)
      puhBaseLCU[uiCUAddr + ui] = uiParameter;
    for (UInt ui = 0; ui < uiCurrPartNumQ; ui++)
      puhBaseLCU[uiCUAddr + 2 * uiCurrPartNumQ + ui] = uiParameter;
    break;
  case SIZE_NxN:
    for (UInt ui = 0; ui < uiCurrPartNumQ; ui++)
      puhBaseLCU[uiCUAddr + ui] = uiParameter;
    break;
  case SIZE_2NxnU:
    if (uiPUIdx == 0)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + uiCurrPartNumQ + ui] = uiParameter;

    }
    else if (uiPUIdx == 1)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1) + (uiCurrPartNumQ << 1); ui++)
        puhBaseLCU[uiCUAddr + uiCurrPartNumQ + ui] = uiParameter;

    }
    else
    {
      assert(0);
    }
    break;
  case SIZE_2NxnD:
    if (uiPUIdx == 0)
    {
      for (UInt ui = 0; ui < ((uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1)); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + uiCurrPartNumQ + ui] = uiParameter;

    }
    else if (uiPUIdx == 1)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 1); ui++)
        puhBaseLCU[uiCUAddr + uiCurrPartNumQ + ui] = uiParameter;

    }
    else
    {
      assert(0);
    }
    break;
  case SIZE_nLx2N:
    if (uiPUIdx == 0)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ >> 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1) + ui] = uiParameter;

    }
    else if (uiPUIdx == 1)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ >> 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1) + ui] = uiParameter;

    }
    else
    {
      assert(0);
    }
    break;
  case SIZE_nRx2N:
    if (uiPUIdx == 0)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + uiCurrPartNumQ + (uiCurrPartNumQ >> 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ + (uiCurrPartNumQ >> 2)); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + uiCurrPartNumQ + (uiCurrPartNumQ >> 1) + ui] = uiParameter;

    }
    else if (uiPUIdx == 1)
    {
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ >> 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + ui] = uiParameter;
      for (UInt ui = 0; ui < (uiCurrPartNumQ >> 2); ui++)
        puhBaseLCU[uiCUAddr + (uiCurrPartNumQ << 1) + (uiCurrPartNumQ >> 1) + ui] = uiParameter;

    }
    else
    {
      assert(0);
    }
    break;
  default:
    assert(0);
  }
}

Void TComDataCU::setMergeComplexFlagSubParts(Bool bMergeComplexFlag, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth)
{
  setSubPart(bMergeComplexFlag, m_pbMergeComplexFlag, uiAbsPartIdx, uiDepth, uiPartIdx);
}

Void TComDataCU::setNodeMotionInfoSubParts(NodeInfo cNodeMotionInfo, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth)
{
  setSubPartT<NodeInfo>(cNodeMotionInfo, m_pNodeMotionInfo, uiAbsPartIdx, uiDepth, uiPartIdx);
}

Void TComDataCU::setNumMergeComplexCandSubParts(UChar NumMergeComplexCand, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth)
{
  setSubPartT<UChar>(NumMergeComplexCand, m_pNumMergeComplexCand, uiAbsPartIdx, uiDepth, uiPartIdx);
}

Void TComDataCU::setMergeComplexIndexSubParts(UInt uiMergeComplexIndex, UInt uiAbsPartIdx, UInt uiPartIdx, UInt uiDepth)
{
  setSubPart<UChar>(uiMergeComplexIndex, m_puhMergeComplexIndex, uiAbsPartIdx, uiDepth, uiPartIdx);
}
#endif

Void  TComDataCU::setMotionInfoField(UInt MergeComplexIndex, UChar uhInterDirComplex, PartSize ePartSize, UInt uiAbsPartIdx, Int iPUIdx)
{
  NodeInfo NodeMotionInfo = getNodeMotionInfo(uiAbsPartIdx);
  if (uhInterDirComplex == 1)
  {
    getCUMvField(REF_PIC_LIST_0)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_0, uiAbsPartIdx, 0, iPUIdx); // interprets depth relative to rpcTempCU level
    getCUMvField(REF_PIC_LIST_1)->setAllMvField(TComMvField(), ePartSize, uiAbsPartIdx, 0, iPUIdx);
  }
  if (uhInterDirComplex == 2)
  {
    getCUMvField(REF_PIC_LIST_0)->setAllMvField(TComMvField(), ePartSize, uiAbsPartIdx, 0, iPUIdx);
    getCUMvField(REF_PIC_LIST_1)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_1, uiAbsPartIdx, 0, iPUIdx);
  }
  if (uhInterDirComplex == 3)
  {
    getCUMvField(REF_PIC_LIST_0)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_0, uiAbsPartIdx, 0, iPUIdx);
    getCUMvField(REF_PIC_LIST_1)->setAllMvFieldNode(NodeMotionInfo, MergeComplexIndex, ePartSize, REF_PIC_LIST_1, uiAbsPartIdx, 0, iPUIdx);
  }
}

Void TComDataCU::preproccessMotionFieldNodes(Int* RealComplexIndex, Int MergeComplexIndex, UChar* InterDirNeighbourNodes, UChar& uhInterDirComplex, Int* numValidMergeNode, TComMvField* cMvFieldNeighbourNodes)
{
  TComSlice*  pSlice = getSlice();
  TComMvField MvFieldNbNodes[8];
  Int wTemp[4] = { 0 }; //temp weighting factor for selecting reference index
  Int NodeNum;
  if (RealComplexIndex[MergeComplexIndex] == 0)
  {
    NodeNum = 4;
  }
  else if (RealComplexIndex[MergeComplexIndex] < 5)
  {
    NodeNum = 3;
  }
  else
  {
    NodeNum = 2;
  }

  // pre process the control motion information: decide inter prediction direction
  switch (RealComplexIndex[MergeComplexIndex])
  {
  case 0:
    for (Int i = 0; i < 8; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[1] & InterDirNeighbourNodes[2] & InterDirNeighbourNodes[3];
    break;

  case 1:
    for (Int i = 0; i < 6; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[1] & InterDirNeighbourNodes[2];
    break;

  case 2:
    for (Int i = 0; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    for (Int i = 4; i < 6; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[1] & InterDirNeighbourNodes[3];
    break;

  case 3:
    for (Int i = 0; i < 2; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    for (Int i = 2; i < 6; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[2] & InterDirNeighbourNodes[3];
    break;

  case 4:
    for (Int i = 0; i < 6; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[1] & InterDirNeighbourNodes[2] & InterDirNeighbourNodes[3];
    break;

  case 5:
    for (Int i = 0; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[1];
    break;

  case 6:
    for (Int i = 0; i < 2; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    for (Int i = 2; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[2];
    break;

  case 7:
    for (Int i = 0; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[1] & InterDirNeighbourNodes[2];
    break;

  case 8:
    for (Int i = 0; i < 2; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i].getMv(), cMvFieldNeighbourNodes[i].getRefIdx());
    }
    for (Int i = 2; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 4].getMv(), cMvFieldNeighbourNodes[i + 4].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[0] & InterDirNeighbourNodes[3];
    break;

  case 9:
    for (Int i = 0; i < 2; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 2].getMv(), cMvFieldNeighbourNodes[i + 2].getRefIdx());
    }
    for (Int i = 2; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 4].getMv(), cMvFieldNeighbourNodes[i + 4].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[1] & InterDirNeighbourNodes[3];
    break;

  case 10:
    for (Int i = 0; i < 4; i++)
    {
      MvFieldNbNodes[i].setMvField(cMvFieldNeighbourNodes[i + 4].getMv(), cMvFieldNeighbourNodes[i + 4].getRefIdx());
    }
    uhInterDirComplex = InterDirNeighbourNodes[2] & InterDirNeighbourNodes[3];
    break;

  default:
    assert(0);
  }

  // pre process the control motion information: decide reference index and scale motion vector
  for (Int l = 0; l < 2; l++)
  {
    RefPicList eRefPicList = RefPicList(l);
    if (((l == 0) && (uhInterDirComplex == 2)) || ((l == 1) && (uhInterDirComplex == 1)))
    {
      continue;
    }
    Int max;
    if (numValidMergeNode[(RealComplexIndex[MergeComplexIndex] << 1) + l] > 1)
    {
      Int  countRefIdx = 1;
      for (Int i = 1; i < NodeNum; i++)
      {
        Int j;
        for (j = 0; j < i; j++)
        {
          if (MvFieldNbNodes[(i << 1) + l].getRefIdx() == MvFieldNbNodes[(j << 1) + l].getRefIdx())
            break;
        }
        if (j == i)
          countRefIdx++;
      }
      // set the target reference index
      for (Int j = 0; j < NodeNum; j++)
      {
        wTemp[j] = 1;
      }

      for (Int j = 0; j < (NodeNum - 1); j++)
      {
        if (wTemp[j])
        {
          for (Int k = j + 1;k < NodeNum;k++)
          {
            if (MvFieldNbNodes[(j << 1) + l].getRefIdx() == MvFieldNbNodes[(k << 1) + l].getRefIdx())
            {
              wTemp[j] += wTemp[k];
              wTemp[k] = 0;
            }
          }
        }
      }
      max = 0;
      for (Int j = 1; j < NodeNum; j++)
      {
        if ((wTemp[j] > wTemp[max]) || ((wTemp[j] == wTemp[max]) && (MvFieldNbNodes[(j << 1) + l].getRefIdx() < MvFieldNbNodes[(max << 1) + l].getRefIdx())))
          max = j;
      }
      //when a MVP is derived from a MV pointing to a different reference picture, the MV is scaled to the target reference picture as the final MVP
      TComMv sMV[4];
      if (countRefIdx == 1)
      {
        for (Int j = 0; j < NodeNum; j++)
        {
          sMV[j] = MvFieldNbNodes[(j << 1) + l].getMv();
        }
      }
      else
      {
        for (Int j = 0; j < NodeNum; j++)
        {
          if (pSlice->isEqualRef(eRefPicList, MvFieldNbNodes[(max << 1) + l].getRefIdx(), MvFieldNbNodes[(j << 1) + l].getRefIdx()))
          {
            sMV[j] = MvFieldNbNodes[(j << 1) + l].getMv();
          }
          else
          {
            Int iCurrPOC = pSlice->getPOC();
            Int iCurrRefPOC = pSlice->getRefPOC(eRefPicList, MvFieldNbNodes[(max << 1) + l].getRefIdx());
            Int iNeibPOC = iCurrPOC;
            Int iNeibRefPOC = pSlice->getRefPOC(eRefPicList, MvFieldNbNodes[(j << 1) + l].getRefIdx());
            Int iDiffPocD = iNeibPOC - iNeibRefPOC;
            Int iDiffPocB = iCurrPOC - iCurrRefPOC;
            Int iScale;
            if (iDiffPocD == iDiffPocB)
            {
              iScale = 1024;
            }
            else
            {
              Int iTDB = Clip3(-128, 127, iDiffPocB);
              Int iTDD = Clip3(-128, 127, iDiffPocD);
              Int iX = (0x4000 + abs(iTDD / 2)) / iTDD;
              iScale = Clip3(-1024, 1023, (iTDB * iX + 32) >> 6);
            }

            if (iScale == 1024)
            {
              sMV[j] = MvFieldNbNodes[(j << 1) + l].getMv();
            }
            else
            {
              sMV[j] = MvFieldNbNodes[(j << 1) + l].getMv().scaleMv(iScale);
            }
          }
        }
      }

      for (Int i = 0; i < NodeNum; i++)
      {
        MvFieldNbNodes[(i << 1) + l].setMvField(sMV[i], MvFieldNbNodes[(max << 1) + l].getRefIdx());
      }
    }
  }

  for (Int i = 0; i < NodeNum; i++)
  {
    cMvFieldNeighbourNodes[(i << 1)].setMvField(MvFieldNbNodes[(i << 1)].getMv(), MvFieldNbNodes[(i << 1)].getRefIdx());
    cMvFieldNeighbourNodes[(i << 1) + 1].setMvField(MvFieldNbNodes[(i << 1) + 1].getMv(), MvFieldNbNodes[(i << 1) + 1].getRefIdx());
  }
}

Void TComDataCU::computeValidNodeNum(Bool* bExistFlagTemp, TComMvField* pcMvFieldNeighbourNodes, Int* numValidMergeNode, Int ModelIndex)
{
  // for list 0
  for (UInt uiOuter = 0; uiOuter < MRG_NODE_NUM; uiOuter++)
  {
    if (bExistFlagTemp[uiOuter << 1])
    {
      for (UInt uiIter = uiOuter + 1; uiIter < MRG_NODE_NUM; uiIter++)
      {
        if (bExistFlagTemp[uiIter << 1])
        {
          UInt uiMvFieldNeighIdxCurr = uiOuter << 1;
          UInt uiMvFieldNeighIdxComp = uiIter << 1;

          if ((pcMvFieldNeighbourNodes[uiMvFieldNeighIdxCurr].getRefIdx() == pcMvFieldNeighbourNodes[uiMvFieldNeighIdxComp].getRefIdx()) && (pcMvFieldNeighbourNodes[uiMvFieldNeighIdxCurr].getMv() == pcMvFieldNeighbourNodes[uiMvFieldNeighIdxComp].getMv()))
          {
            bExistFlagTemp[uiIter << 1] = false;
          }
        }
      }
    }
  }

  // for list 1
  for (UInt uiOuter = 0; uiOuter < MRG_NODE_NUM; uiOuter++)
  {
    if (bExistFlagTemp[(uiOuter << 1) + 1])
    {
      for (UInt uiIter = uiOuter + 1; uiIter < MRG_NODE_NUM; uiIter++)
      {
        if (bExistFlagTemp[(uiIter << 1) + 1])
        {
          UInt uiMvFieldNeighIdxCurr = (uiOuter << 1) + 1;
          UInt uiMvFieldNeighIdxComp = (uiIter << 1) + 1;

          if ((pcMvFieldNeighbourNodes[uiMvFieldNeighIdxCurr].getRefIdx() == pcMvFieldNeighbourNodes[uiMvFieldNeighIdxComp].getRefIdx()) && (pcMvFieldNeighbourNodes[uiMvFieldNeighIdxCurr].getMv() == pcMvFieldNeighbourNodes[uiMvFieldNeighIdxComp].getMv()))
          {
            bExistFlagTemp[(uiIter << 1) + 1] = false;
          }
        }
      }
    }
  }

  for (Int i = 0; i < MRG_NODE_NUM; i++)
  {
    if (bExistFlagTemp[i << 1])
      numValidMergeNode[ModelIndex << 1]++;
    if (bExistFlagTemp[(i << 1) + 1])
      numValidMergeNode[(ModelIndex << 1) + 1]++;
  }
}

Void TComDataCU::getValidPredictionModels(Bool* bExistFlag, TComMvField* pcMvFieldNeighbourNodes, Int* numValidMergeNode, Int* numExistMergeNode, Int* RealComplexIndex, UChar& numValidMergeComplexCand, Int mrgComplexCandIdx)
{
  Bool bExistFlagTemp[MRG_NODE_NUM << 1];
  for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
  {
    bExistFlagTemp[i] = bExistFlag[i];
  }
  computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 0);
  Int bValidList[2] = { 0 };
  Int ComplexCount = 0;
  if (((numValidMergeNode[0] > 1) && (numExistMergeNode[0] == 4)) || ((numValidMergeNode[1] > 1) && (numExistMergeNode[1] == 4)))
  {
    // whether the index 7 is valid
    bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[4] == true));
    bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[5] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[0] = false;
      bExistFlagTemp[1] = false;
      bExistFlagTemp[6] = false;
      bExistFlagTemp[7] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 7);
    }
    if (((numValidMergeNode[14] > 1) && bValidList[0]) || ((numValidMergeNode[15] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 7;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 6 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[4] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[5] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[2] = false;
      bExistFlagTemp[3] = false;
      bExistFlagTemp[6] = false;
      bExistFlagTemp[7] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 6);
    }
    if (((numValidMergeNode[12] > 1) && bValidList[0]) || ((numValidMergeNode[13] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 6;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 1 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[2] == true) && (bExistFlag[4] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[3] == true) && (bExistFlag[5] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[6] = false;
      bExistFlagTemp[7] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 1);
    }
    if (((numValidMergeNode[2] > 1) && bValidList[0]) || ((numValidMergeNode[3] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 1;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 5 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[2] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[3] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[4] = false;
      bExistFlagTemp[5] = false;
      bExistFlagTemp[6] = false;
      bExistFlagTemp[7] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 5);
    }
    if (((numValidMergeNode[10] > 1) && bValidList[0]) || ((numValidMergeNode[11] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 5;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 9 is valid
    bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[7] == true));

    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[0] = false;
      bExistFlagTemp[1] = false;
      bExistFlagTemp[4] = false;
      bExistFlagTemp[5] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 9);
    }
    if (((numValidMergeNode[18] > 1) && bValidList[0]) || ((numValidMergeNode[19] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 9;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    // whether the index 10 is valid
    bValidList[0] = ((bExistFlag[4] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[5] == true) && (bExistFlag[7] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[0] = false;
      bExistFlagTemp[1] = false;
      bExistFlagTemp[2] = false;
      bExistFlagTemp[3] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 10);
    }
    if (((numValidMergeNode[20] > 1) && bValidList[0]) || ((numValidMergeNode[21] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 10;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 8 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[7] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[2] = false;
      bExistFlagTemp[3] = false;
      bExistFlagTemp[4] = false;
      bExistFlagTemp[5] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 8);
    }
    if (((numValidMergeNode[16] > 1) && bValidList[0]) || ((numValidMergeNode[17] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 8;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //index 0
    RealComplexIndex[ComplexCount] = 0;
    if (mrgComplexCandIdx == ComplexCount)
    {
      numValidMergeComplexCand = mrgComplexCandIdx + 1;
      return;
    }
    ComplexCount++;

    // whether the index 2 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[2] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[3] == true) && (bExistFlag[7] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[4] = false;
      bExistFlagTemp[5] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 2);
    }
    if (((numValidMergeNode[4] > 1) && bValidList[0]) || ((numValidMergeNode[5] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 2;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    // whether the index 4 is valid
    bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[4] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[5] == true) && (bExistFlag[7] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[0] = false;
      bExistFlagTemp[1] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 4);
    }
    if (((numValidMergeNode[8] > 1) && bValidList[0]) || ((numValidMergeNode[9] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 4;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    //whether the index 3 is valid
    bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[4] == true) && (bExistFlag[6] == true));
    bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[5] == true) && (bExistFlag[7] == true));
    if (bValidList[0] || bValidList[1])
    {
      for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
      {
        bExistFlagTemp[i] = bExistFlag[i];
      }
      bExistFlagTemp[2] = false;
      bExistFlagTemp[3] = false;
      computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 3);
    }
    if (((numValidMergeNode[6] > 1) && bValidList[0]) || ((numValidMergeNode[7] > 1) && bValidList[1]))
    {
      RealComplexIndex[ComplexCount] = 3;
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }
  }
  else if (((numValidMergeNode[0] > 1) && (numExistMergeNode[0] == 3)) || ((numValidMergeNode[1] > 1) && (numExistMergeNode[1] == 3)))
  {
    if ((bExistFlag[6] == false) && (bExistFlag[7] == false))
    {
      //whether the index 7 is valid
      bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[4] == true));
      bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[5] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[0] = false;
        bExistFlagTemp[1] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 7);
      }
      if (((numValidMergeNode[14] > 1) && bValidList[0]) || ((numValidMergeNode[15] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 7;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 6 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[4] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[5] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[2] = false;
        bExistFlagTemp[3] = false;

        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 6);
      }
      if (((numValidMergeNode[12] > 1) && bValidList[0]) || ((numValidMergeNode[13] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 6;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //index 1
      RealComplexIndex[ComplexCount] = 1;
      numValidMergeNode[2] = numValidMergeNode[0];
      numValidMergeNode[3] = numValidMergeNode[1];
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;

      //whether the index 5 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[2] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[3] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[4] = false;
        bExistFlagTemp[5] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 5);
      }
      if (((numValidMergeNode[10] > 1) && bValidList[0]) || ((numValidMergeNode[11] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 5;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }
    }

    if ((bExistFlag[4] == false) && (bExistFlag[5] == false))
    {
      //whether the index 5 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[2] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[3] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }

        bExistFlagTemp[6] = false;
        bExistFlagTemp[7] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 5);
      }
      if (((numValidMergeNode[10] > 1) && bValidList[0]) || ((numValidMergeNode[11] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 5;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 9 is valid
      bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[7] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[0] = false;
        bExistFlagTemp[1] = false;

        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 9);
      }
      if (((numValidMergeNode[18] > 1) && bValidList[0]) || ((numValidMergeNode[19] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 9;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 8 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[7] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[2] = false;
        bExistFlagTemp[3] = false;

        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 8);
      }
      if (((numValidMergeNode[16] > 1) && bValidList[0]) || ((numValidMergeNode[17] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 8;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //index 2
      RealComplexIndex[ComplexCount] = 2;
      numValidMergeNode[4] = numValidMergeNode[0];
      numValidMergeNode[5] = numValidMergeNode[1];
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    if ((bExistFlag[2] == false) && (bExistFlag[3] == false))
    {
      //whether the index 6 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[4] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[5] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }

        bExistFlagTemp[6] = false;
        bExistFlagTemp[7] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 6);
      }
      if (((numValidMergeNode[12] > 1) && bValidList[0]) || ((numValidMergeNode[13] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 6;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 10 is valid
      bValidList[0] = ((bExistFlag[4] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[5] == true) && (bExistFlag[7] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[0] = false;
        bExistFlagTemp[1] = false;

        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 10);
      }
      if (((numValidMergeNode[20] > 1) && bValidList[0]) || ((numValidMergeNode[21] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 10;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 8 is valid
      bValidList[0] = ((bExistFlag[0] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[1] == true) && (bExistFlag[7] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }

        bExistFlagTemp[4] = false;
        bExistFlagTemp[5] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 8);
      }
      if (((numValidMergeNode[16] > 1) && bValidList[0]) || ((numValidMergeNode[17] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 8;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //index 3
      RealComplexIndex[ComplexCount] = 3;
      numValidMergeNode[6] = numValidMergeNode[0];
      numValidMergeNode[7] = numValidMergeNode[1];
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }

    if ((bExistFlag[0] == false) && (bExistFlag[1] == false))
    {
      //whether the index 7 is valid
      bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[4] == true));
      bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[5] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[6] = false;
        bExistFlagTemp[7] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 7);
      }
      if (((numValidMergeNode[14] > 1) && bValidList[0]) || ((numValidMergeNode[15] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 7;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 9 is valid
      bValidList[0] = ((bExistFlag[2] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[3] == true) && (bExistFlag[7] == true));
      if (bValidList[0] || bValidList[1])
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[4] = false;
        bExistFlagTemp[5] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 9);
      }
      if (((numValidMergeNode[18] > 1) && bValidList[0]) || ((numValidMergeNode[19] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 9;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //whether the index 10 is valid
      bValidList[0] = ((bExistFlag[4] == true) && (bExistFlag[6] == true));
      bValidList[1] = ((bExistFlag[5] == true) && (bExistFlag[7] == true));
      {
        for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
        {
          bExistFlagTemp[i] = bExistFlag[i];
        }
        bExistFlagTemp[2] = false;
        bExistFlagTemp[3] = false;
        computeValidNodeNum(bExistFlagTemp, pcMvFieldNeighbourNodes, numValidMergeNode, 10);
      }
      if (((numValidMergeNode[20] > 1) && bValidList[0]) || ((numValidMergeNode[21] > 1) && bValidList[1]))
      {
        RealComplexIndex[ComplexCount] = 10;
        if (mrgComplexCandIdx == ComplexCount)
        {
          numValidMergeComplexCand = mrgComplexCandIdx + 1;
          return;
        }
        ComplexCount++;
      }

      //index 4
      RealComplexIndex[ComplexCount] = 4;
      numValidMergeNode[8] = numValidMergeNode[0];
      numValidMergeNode[9] = numValidMergeNode[1];
      if (mrgComplexCandIdx == ComplexCount)
      {
        numValidMergeComplexCand = mrgComplexCandIdx + 1;
        return;
      }
      ComplexCount++;
    }
  }
  else if (((numValidMergeNode[0] > 1) && (numExistMergeNode[0] == 2)) || ((numValidMergeNode[1] > 1) && (numExistMergeNode[1] == 2)))
  {
    if ((bExistFlag[6] == false) && (bExistFlag[7] == false) && (bExistFlag[4] == false) && (bExistFlag[5] == false))
    {
      numValidMergeNode[10] = numValidMergeNode[0];
      numValidMergeNode[11] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 5;
    }
    if ((bExistFlag[6] == false) && (bExistFlag[7] == false) && (bExistFlag[2] == false) && (bExistFlag[3] == false))
    {
      numValidMergeNode[12] = numValidMergeNode[0];
      numValidMergeNode[13] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 6;
    }
    if ((bExistFlag[6] == false) && (bExistFlag[7] == false) && (bExistFlag[0] == false) && (bExistFlag[1] == false))
    {
      numValidMergeNode[14] = numValidMergeNode[0];
      numValidMergeNode[15] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 7;
    }
    if ((bExistFlag[4] == false) && (bExistFlag[5] == false) && (bExistFlag[2] == false) && (bExistFlag[3] == false))
    {
      numValidMergeNode[16] = numValidMergeNode[0];
      numValidMergeNode[17] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 8;
    }
    if ((bExistFlag[4] == false) && (bExistFlag[5] == false) && (bExistFlag[0] == false) && (bExistFlag[1] == false))
    {
      numValidMergeNode[18] = numValidMergeNode[0];
      numValidMergeNode[19] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 9;
    }
    if ((bExistFlag[2] == false) && (bExistFlag[3] == false) && (bExistFlag[0] == false) && (bExistFlag[1] == false))
    {
      numValidMergeNode[20] = numValidMergeNode[0];
      numValidMergeNode[21] = numValidMergeNode[1];
      RealComplexIndex[ComplexCount++] = 10;
    }
  }
  numValidMergeComplexCand = ComplexCount;
  if (numValidMergeComplexCand > MRG_COMPLEX_NUM_CANDS_REAL)
  {
    numValidMergeComplexCand = MRG_COMPLEX_NUM_CANDS_REAL;
  }
}

Void TComDataCU::getInterMergeComplexCandidates(UInt uiAbsPartIdx, UInt uiPUIdx, UChar& numValidMergeComplexCand, TComMvField* pcMvFieldNeighbourNodes, UChar* puhInterDirNeighbourNodes, Int* numExistMergeNode, Int* numValidMergeNode, Int* RealComplexIndex, Int mrgComplexCandIdx)
{
  UInt uiAbsPartAddr = m_absZIdxInCtu + uiAbsPartIdx;
  Bool bExistFlag[MRG_NODE_NUM << 1]; //Represent motion information of four corner node exist or not, 0 ,1 LT;2,3RT;4,5LB;6,7RB

  for (UInt ui = 0; ui < MRG_NODE_NUM << 1; ++ui)
  {
    bExistFlag[ui] = false;
  }
  for (UInt ui = 0; ui < MRG_NODE_NUM; ++ui)
  {
    pcMvFieldNeighbourNodes[(ui << 1)].setRefIdx(NOT_VALID);
    pcMvFieldNeighbourNodes[(ui << 1) + 1].setRefIdx(NOT_VALID);
  }

  UInt uiPartIdxLT, uiPartIdxRT, uiPartIdxLB;
#if !JVET_C0024_QTBT
  PartSize cCurPS = getPartitionSize(uiAbsPartIdx);
#endif
  deriveLeftRightTopIdxGeneral(uiAbsPartIdx, uiPUIdx, uiPartIdxLT, uiPartIdxRT);
  deriveLeftBottomIdxGeneral(uiAbsPartIdx, uiPUIdx, uiPartIdxLB);

  //above left
  UInt uiAboveLeftPartIdx = 0;
  TComDataCU* pcCUAboveLeft = 0;
  UInt uiLeftPartIdx1 = 0;
  UInt uiAbovePartIdx1 = 0;
  TComDataCU* pcCULeft1 = 0;
  TComDataCU* pcCUAbove1 = 0;
  pcCUAboveLeft = getPUAboveLeft(uiAboveLeftPartIdx, uiAbsPartAddr);
  pcCULeft1 = getPULeft(uiLeftPartIdx1, uiAbsPartAddr);
  pcCUAbove1 = getPUAbove(uiAbovePartIdx1, uiAbsPartAddr);
  if (pcCUAboveLeft && !pcCUAboveLeft->isIntra(uiAboveLeftPartIdx))
  {
    //get motion information for node 0
    pcCUAboveLeft->getMvField(pcCUAboveLeft, uiAboveLeftPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[0 << 1]);
    if (pcMvFieldNeighbourNodes[0 << 1].getRefIdx() >= 0)
    {
      bExistFlag[0 << 1] = true;
      numExistMergeNode[0]++;
    }

    if (getSlice()->isInterB())
    {
      pcCUAboveLeft->getMvField(pcCUAboveLeft, uiAboveLeftPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(0 << 1) + 1]);
      if (pcMvFieldNeighbourNodes[(0 << 1) + 1].getRefIdx() >= 0)
      {
        bExistFlag[(0 << 1) + 1] = true;
        numExistMergeNode[1]++;
      }
    }
  }

  // left
  if (pcCULeft1 && !pcCULeft1->isIntra(uiLeftPartIdx1))
  {
    //get motion information for node 0
    if (bExistFlag[0 << 1] == false)
    {
      pcCULeft1->getMvField(pcCULeft1, uiLeftPartIdx1, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[0 << 1]);
      if (pcMvFieldNeighbourNodes[0 << 1].getRefIdx() >= 0)
      {
        bExistFlag[0 << 1] = true;
        numExistMergeNode[0]++;
      }
    }

    // get Mv from Left
    if (getSlice()->isInterB())
    {
      if (bExistFlag[(0 << 1) + 1] == false)
      {
        pcCULeft1->getMvField(pcCULeft1, uiLeftPartIdx1, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(0 << 1) + 1]);
        if (pcMvFieldNeighbourNodes[(0 << 1) + 1].getRefIdx() >= 0)
        {
          bExistFlag[(0 << 1) + 1] = true;
          numExistMergeNode[1]++;
        }
      }
    }
  }

  // above
  if (pcCUAbove1 && !pcCUAbove1->isIntra(uiAbovePartIdx1))
  {
    //get motion information for node 0
    if (bExistFlag[0 << 1] == false)
    {
      pcCUAbove1->getMvField(pcCUAbove1, uiAbovePartIdx1, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[0 << 1]);
      if (pcMvFieldNeighbourNodes[0 << 1].getRefIdx() >= 0)
      {
        bExistFlag[0 << 1] = true;
        numExistMergeNode[0]++;
      }
    }

    // get Mv from Left
    if (getSlice()->isInterB())
    {
      if (bExistFlag[(0 << 1) + 1] == false)
      {
        pcCUAbove1->getMvField(pcCUAbove1, uiAbovePartIdx1, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(0 << 1) + 1]);
        if (pcMvFieldNeighbourNodes[(0 << 1) + 1].getRefIdx() >= 0)
        {
          bExistFlag[(0 << 1) + 1] = true;
          numExistMergeNode[1]++;
        }
      }
    }
  }

  //above right
  UInt uiAboveRightPartIdx = 0;
  TComDataCU* pcCUAboveRight = 0;
  pcCUAboveRight = getPUAboveRight(uiAboveRightPartIdx, uiPartIdxRT);
  if (pcCUAboveRight && !pcCUAboveRight->isIntra(uiAboveRightPartIdx))
  {
    pcCUAboveRight->getMvField(pcCUAboveRight, uiAboveRightPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[1 << 1]);
    if (pcMvFieldNeighbourNodes[1 << 1].getRefIdx() >= 0)
    {
      bExistFlag[1 << 1] = true;
      numExistMergeNode[0]++;
    }

    if (getSlice()->isInterB())
    {
      pcCUAboveRight->getMvField(pcCUAboveRight, uiAboveRightPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(1 << 1) + 1]);
      if (pcMvFieldNeighbourNodes[(1 << 1) + 1].getRefIdx() >= 0)
      {
        bExistFlag[(1 << 1) + 1] = true;
        numExistMergeNode[1]++;
      }
    }
  }

  // above
  UInt uiAbovePartIdx = 0;
  TComDataCU* pcCUAbove = 0;
  pcCUAbove = getPUAbove(uiAbovePartIdx, uiPartIdxRT);
  if (pcCUAbove && !pcCUAbove->isIntra(uiAbovePartIdx))
  {
    if (bExistFlag[2] == false)
    {
      pcCUAbove->getMvField(pcCUAbove, uiAbovePartIdx, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[1 << 1]);
      if (pcMvFieldNeighbourNodes[1 << 1].getRefIdx() >= 0)
      {
        bExistFlag[1 << 1] = true;
        numExistMergeNode[0]++;
      }
    }

    if (getSlice()->isInterB())
    {
      if (bExistFlag[3] == false)
      {
        pcCUAbove->getMvField(pcCUAbove, uiAbovePartIdx, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(1 << 1) + 1]);
        if (pcMvFieldNeighbourNodes[(1 << 1) + 1].getRefIdx() >= 0)
        {
          bExistFlag[(1 << 1) + 1] = true;
          numExistMergeNode[1]++;
        }
      }
    }
  }

  //left bottom
  UInt uiLeftBottomPartIdx = 0;
  TComDataCU* pcCULeftBottom = 0;
  pcCULeftBottom = this->getPUBelowLeft(uiLeftBottomPartIdx, uiPartIdxLB);
  if (pcCULeftBottom && !pcCULeftBottom->isIntra(uiLeftBottomPartIdx))
  {
    pcCULeftBottom->getMvField(pcCULeftBottom, uiLeftBottomPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[2 << 1]);
    if (pcMvFieldNeighbourNodes[2 << 1].getRefIdx() >= 0)
    {
      bExistFlag[2 << 1] = true;
      numExistMergeNode[0]++;
    }

    if (getSlice()->isInterB())
    {
      pcCULeftBottom->getMvField(pcCULeftBottom, uiLeftBottomPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(2 << 1) + 1]);
      if (pcMvFieldNeighbourNodes[(2 << 1) + 1].getRefIdx() >= 0)
      {
        bExistFlag[(2 << 1) + 1] = true;
        numExistMergeNode[1]++;
      }
    }
  }

  UInt uiLeftPartIdx = 0;
  TComDataCU* pcCULeft = 0;
  pcCULeft = getPULeft(uiLeftPartIdx, uiPartIdxLB);
  if (pcCULeft && !pcCULeft->isIntra(uiLeftPartIdx))
  {
    if (bExistFlag[4] == false)
    {
      pcCULeft->getMvField(pcCULeft, uiLeftPartIdx, REF_PIC_LIST_0, pcMvFieldNeighbourNodes[2 << 1]);
      if (pcMvFieldNeighbourNodes[2 << 1].getRefIdx() >= 0)
      {
        bExistFlag[2 << 1] = true;
        numExistMergeNode[0]++;
      }
    }

    if (getSlice()->isInterB())
    {
      if (bExistFlag[5] == false)
      {
        pcCULeft->getMvField(pcCULeft, uiLeftPartIdx, REF_PIC_LIST_1, pcMvFieldNeighbourNodes[(2 << 1) + 1]);
        if (pcMvFieldNeighbourNodes[(2 << 1) + 1].getRefIdx() >= 0)
        {
          bExistFlag[(2 << 1) + 1] = true;
          numExistMergeNode[1]++;
        }
      }
    }
  }

  //get  inter directions for the first three control nodes
  for (Int i = 0; i < 3; i++)
  {
    if ((bExistFlag[i << 1] == true) && (bExistFlag[(i << 1) + 1] == true))
      puhInterDirNeighbourNodes[i] = 3;
    else if (bExistFlag[i << 1] == true)
      puhInterDirNeighbourNodes[i] = 1;
    else if (bExistFlag[(i << 1) + 1] == true)
      puhInterDirNeighbourNodes[i] = 2;
  }

  if (getSlice()->getEnableTMVPFlag())
  {
    //>> MTK colocated-RightBottom
    UInt uiPartIdxRB;

    deriveRightBottomIdx(uiPUIdx, uiPartIdxRB);

    UInt uiAbsPartIdxTmp = g_auiZscanToRaster[uiPartIdxRB];
    const UInt numPartInCtuWidth = m_pcPic->getNumPartInCtuWidth();
    const UInt numPartInCtuHeight = m_pcPic->getNumPartInCtuHeight();

    TComMv cColMv;
    Int iRefIdx;
    Int ctuRsAddr = -1;

    if (((m_pcPic->getCtu(m_ctuRsAddr)->getCUPelX() + g_auiRasterToPelX[uiAbsPartIdxTmp] + m_pcPic->getMinCUWidth()) < m_pcSlice->getSPS()->getPicWidthInLumaSamples())  // image boundary check
      && ((m_pcPic->getCtu(m_ctuRsAddr)->getCUPelY() + g_auiRasterToPelY[uiAbsPartIdxTmp] + m_pcPic->getMinCUHeight()) < m_pcSlice->getSPS()->getPicHeightInLumaSamples()))
    {
      if ((uiAbsPartIdxTmp % numPartInCtuWidth < numPartInCtuWidth - 1) &&           // is not at the last column of CTU
        (uiAbsPartIdxTmp / numPartInCtuWidth < numPartInCtuHeight - 1))              // is not at the last row    of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[uiAbsPartIdxTmp + numPartInCtuWidth + 1];
        ctuRsAddr = getCtuRsAddr();
      }
      else if (uiAbsPartIdxTmp % numPartInCtuWidth < numPartInCtuWidth - 1)           // is not at the last column of CTU But is last row of CTU
      {
#if COM16_C806_GEN_MRG_IMPROVEMENT
        ctuRsAddr = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus();
#endif
        uiAbsPartAddr = g_auiRasterToZscan[(uiAbsPartIdxTmp + numPartInCtuWidth + 1) % m_pcPic->getNumPartitionsInCtu()];
      }
      else if (uiAbsPartIdxTmp / numPartInCtuWidth < numPartInCtuHeight - 1)          // is not at the last row of CTU But is last column of CTU
      {
        uiAbsPartAddr = g_auiRasterToZscan[uiAbsPartIdxTmp + 1];
        ctuRsAddr = getCtuRsAddr() + 1;
      }
      else //is the right bottom corner of CTU
      {
        uiAbsPartAddr = 0;
#if COM16_C806_GEN_MRG_IMPROVEMENT
        ctuRsAddr = getCtuRsAddr() + m_pcPic->getFrameWidthInCtus() + 1;
#endif
      }
    }

    iRefIdx = 0;
    Bool bExistMV = false;
    bExistMV = ctuRsAddr >= 0 && xGetColMVP(REF_PIC_LIST_0, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx);

    if (bExistMV)
    {
      bExistFlag[3 << 1] = true;
      numExistMergeNode[0]++;
      // get Inter Dir
      puhInterDirNeighbourNodes[3] |= 1;
      pcMvFieldNeighbourNodes[3 << 1].setMvField(cColMv, iRefIdx);
    }

    if (getSlice()->isInterB())
    {
      bExistMV = ctuRsAddr >= 0 && xGetColMVP(REF_PIC_LIST_1, ctuRsAddr, uiAbsPartAddr, cColMv, iRefIdx);

      if (bExistMV)
      {
        bExistFlag[(3 << 1) + 1] = true;
        numExistMergeNode[1]++;
        // get Inter Dir
        puhInterDirNeighbourNodes[3] |= 2;
        pcMvFieldNeighbourNodes[(3 << 1) + 1].setMvField(cColMv, iRefIdx);
      }
    }
  }

  getValidPredictionModels(bExistFlag, pcMvFieldNeighbourNodes, numValidMergeNode, numExistMergeNode, RealComplexIndex, numValidMergeComplexCand, mrgComplexCandIdx);
}
#endif

#if RRC_INTRA_SHARPENING
Void TComDataCU::setIntraSharpFlagSubParts(Char param, UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  setSubPart(param, m_sharpFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
#else
  setSubPart(param, m_sharpFlag, uiAbsPartIdx, uiDepth, uiPartIdx);
#endif
}
Bool TComDataCU::isIntraSharpFlagCoded(UInt uiAbsPartIdx)
{
  if (!getSlice()->getSPS()->getUseIntraSharp()) 
    return false;

  Bool bCodeSharp = isIntra(uiAbsPartIdx);
#if INTRA_SHARP_PDPC_HARM
  if (getPDPCIdx(uiAbsPartIdx))
    bCodeSharp = 0;
#endif
#if RRC_DWDIP
  if (getDWDIPIdx(uiAbsPartIdx))
    bCodeSharp = 0;
#endif


#if MLIP_EXTEND_REF_LINE
  if (getExternRef(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
  {
      bCodeSharp = 0;
  }
#endif

#if INTRA_SHARP_MAX_ASPECT_RATIO
  int AspectRatio = std::max(getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx)) / std::min(getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
  if (AspectRatio > INTRA_SHARP_MAX_ASPECT_RATIO)
    bCodeSharp = 0;
#endif
#if INTRA_SHARP_MAX_SIZE
  if (std::max(getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx)) > INTRA_SHARP_MAX_SIZE)
    bCodeSharp = 0;
#endif
#if DIMD_INTRA_PRED
  if (getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
    bCodeSharp = 0;
#endif
#if MLIP_EXTEND_REF_LINE
  if (getExternRef(CHANNEL_TYPE_LUMA, uiAbsPartIdx))
      bCodeSharp = 0;
#endif

#if !JVET_C0024_QTBT
  if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N)
#endif
  {
    UInt uiOrgMode = getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
    if (uiOrgMode == DC_IDX || uiOrgMode == PLANAR_IDX)
    {
      bCodeSharp = 0;
    }
  }
  return bCodeSharp;
}
#endif

#if BT_FORCE
UInt TComDataCU::getBoundaryCtuLocation()
{
  const TComSPS &sps = *(getSlice()->getSPS());
  const UInt uiLPelX = (getCtuRsAddr() % getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  const UInt uiRPelX = uiLPelX + sps.getCTUSize() - 1;
  const UInt uiTPelY = (getCtuRsAddr() / getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  const UInt uiBPelY = uiTPelY + sps.getCTUSize() - 1;
  if (uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples())
    return 1;
  else if (uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples())
    return 2;
  else
    return 0;
}

#if BT_FORCE_INTRA_SPEED_UP
Bool TComDataCU::isFromForcedQT(UInt uiAbsPartIdx)
{
  if (!getSlice()->isIntra())
    return false;
  if (getBoundaryCtuLocation() == 1)
  {
    Int iCtuSize = getSlice()->getSPS()->getCTUSize();
    Int iRest = getSlice()->getSPS()->getPicHeightInLumaSamples() % iCtuSize;//bottom
    Int iTempSize = iCtuSize;
    Int iFromFQT = 0;
    Int iShift = 1;
#if BT_FORCE_QTN
    Int iMaxQtLv = g_aucConvertToBit[iCtuSize] - g_aucConvertToBit[isLuma(getTextType()) ? getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
    for (Int iCount = 0; iCount < iMaxQtLv; iCount++)
#else
    for (Int iCount = 0; iCount < MAX_FQT_LV_INTRA; iCount++)
#endif
    {
      iTempSize >>= iShift;
      if (iRest > iTempSize)
      {
        iFromFQT += iTempSize;
        iRest -= iTempSize;
      }
    }
    if (iFromFQT == 0)
      return false;
#if TT_TREE_DS1 // deriveZIdxInCTU
#if BARC_BTFORCE
    if (g_auiRasterToPelY[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) << (2 * getBARCFlag(getTextType()))]] < iFromFQT)
#else
    if (g_auiRasterToPelY[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx)]] < iFromFQT)
#endif
#else
    if (g_auiRasterToPelY[g_auiZscanToRaster[m_absZIdxInCtu + uiAbsPartIdx]] < iFromFQT)
#endif
      return true;
    else
      return false;
  }
  if (getBoundaryCtuLocation() == 2)
  {
    Int iCtuSize = getSlice()->getSPS()->getCTUSize();
    Int iRest = getSlice()->getSPS()->getPicWidthInLumaSamples() % iCtuSize;//right
    Int iTempSize = iCtuSize;
    Int iFromFQT = 0;
    Int iShift = 1;
#if BT_FORCE_QTN
    Int iMaxQtLv = g_aucConvertToBit[iCtuSize] - g_aucConvertToBit[isLuma(getTextType()) ? getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
    for (Int iCount = 0; iCount < iMaxQtLv; iCount++)
#else
    for (Int iCount = 0; iCount < MAX_FQT_LV_INTRA; iCount++)
#endif
    {
      iTempSize >>= iShift;
      if (iRest > iTempSize)
      {
        iFromFQT += iTempSize;
        iRest -= iTempSize;
      }
    }
    if (iFromFQT == 0)
      return false;
#if TT_TREE_DS1 // deriveZIdxInCTU
#if BARC_BTFORCE
    if (g_auiRasterToPelX[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) << (2 * getBARCFlag(getTextType()))]] < iFromFQT)
#else
    if (g_auiRasterToPelX[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx)]] < iFromFQT)
#endif
#else
    if (g_auiRasterToPelX[g_auiZscanToRaster[m_absZIdxInCtu + uiAbsPartIdx]] < iFromFQT)
#endif
      return true;
    else
      return false;
  }
  else
    return false;
}
#endif
#if BT_FORCE_TT
Bool TComDataCU::isForcedBT(UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight)
{
  if (getBoundaryCuLocation(uiAbsPartIdx, uiWidth, uiHeight) == 1)
  {
    Int iCtuSize = getSlice()->getSPS()->getCTUSize();
    Int iRest = getSlice()->getSPS()->getPicHeightInLumaSamples() % iCtuSize;//bottom
#if BARC_BTFORCE
    Int iRemainInCU = iRest - g_auiRasterToPelY[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) << (2* getBARCFlag(getTextType()))]];
    if (iRemainInCU == (uiHeight >> !getBARCFlag(getTextType())))
#else
    Int iRemainInCU = iRest - g_auiRasterToPelY[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx)]];
    if (iRemainInCU == (uiHeight >> 1))
#endif
      return true;
    else
      return false;
  }
  if (getBoundaryCuLocation(uiAbsPartIdx, uiWidth, uiHeight) == 2)
  {
    Int iCtuSize = getSlice()->getSPS()->getCTUSize();
    Int iRest = getSlice()->getSPS()->getPicWidthInLumaSamples() % iCtuSize;//right
#if BARC_BTFORCE
    Int iRemainInCU = iRest - g_auiRasterToPelX[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx) << (2 * getBARCFlag(getTextType()))]];
    if (iRemainInCU == (uiWidth >> !getBARCFlag(getTextType())))
#else
    Int iRemainInCU = iRest - g_auiRasterToPelX[g_auiZscanToRaster[deriveZIdxInCTU(m_absZIdxInCtu, uiAbsPartIdx)]];
    if (iRemainInCU == (uiWidth >> 1))
#endif
      return true;
    else
      return false;
  }
  return false;
}

UInt TComDataCU::getBoundaryCuLocation(UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight)
{
  const TComSPS &sps = *(getSlice()->getSPS());
  const UInt uiLPelX = getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  const UInt uiRPelX = uiLPelX + uiWidth - 1;
  const UInt uiTPelY = getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
  const UInt uiBPelY = uiTPelY + uiHeight - 1;

#if BARC_BTFORCE
  UInt uiPelXInCTU = uiLPelX - getPic()->getCtu(getCtuRsAddr())->getCUPelX();
  UInt uiPelYInCTU = uiTPelY - getPic()->getCtu(getCtuRsAddr())->getCUPelY();

  UInt deltaX = (uiPelXInCTU + uiWidth)*getBARCFlag(getTextType());
  UInt deltaY = (uiPelYInCTU + uiHeight)*getBARCFlag(getTextType());

  Bool bBBoundary = ((uiRPelX + deltaX) < sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) >= sps.getPicHeightInLumaSamples());
  Bool bRBoundary = ((uiRPelX + deltaX) >= sps.getPicWidthInLumaSamples()) && ((uiBPelY + deltaY) < sps.getPicHeightInLumaSamples());
#endif

#if BARC_BTFORCE
  if (bBBoundary)
#else
  if (uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples())
#endif
    return 1;
#if BARC_BTFORCE
  else if (bRBoundary)
#else
  else if (uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples())
#endif
    return 2;
  else
    return 0;
}

#endif
#if BT_FORCE_FBTD
UInt TComDataCU::getFBTDepth(UInt uiAbsPartIdx)
{
  const TComSPS &sps = *(getSlice()->getSPS());
  UInt uiLPelX = (getCtuRsAddr() % getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  UInt uiWidth = sps.getCTUSize();
  UInt uiRPelX = uiLPelX + uiWidth - 1;
  UInt uiTPelY = (getCtuRsAddr() / getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  UInt uiHeight = sps.getCTUSize();
  UInt uiBPelY = uiTPelY + uiWidth - 1;
  UInt uiFBTDepth = 0;
  UInt uiDepth = 0;
  m_uiFBTDepth = 0;
  if (uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples())
    return m_uiFBTDepth;
  xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidth, uiHeight, uiFBTDepth, uiDepth);
  return m_uiFBTDepth;
}
Void TComDataCU::xRecurseFBTD(UInt uiAbsPartIdx, UInt uiLPelX, UInt uiTPelY, UInt uiWidth, UInt uiHeight, UInt uiFBTDepth, UInt uiDepth)
{
  const TComSPS &sps = *(getSlice()->getSPS());
#if BT_FORCE_QTN
  const Int iMaxQtLv = g_aucConvertToBit[sps.getCTUSize()] - g_aucConvertToBit[isLuma(getTextType()) ? getSlice()->getMaxBTSize() : MAX_BT_SIZE_C];
#else
  Int iMaxQtLv =
#if BT_FORCE_INTRA_ON
    getSlice()->isIntra() ? MAX_FQT_LV_INTRA :
#endif
    MAX_FQT_LV;
#endif
  UInt uiLPelXOrg = (getCtuRsAddr() % getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  UInt uiTPelYOrg = (getCtuRsAddr() / getPic()->getFrameWidthInCtus())*sps.getCTUSize();
  UInt uiRPelX = uiLPelX + uiWidth - 1;
  UInt uiBPelY = uiTPelY + uiHeight - 1;
  UInt uiCurLPelX = getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  UInt uiCurTPelY = getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
  if (uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples() && (uiDepth >= iMaxQtLv)) //bottom
  {
#if BT_FORCE_TT
    Int iRemainInCU = getSlice()->getSPS()->getPicHeightInLumaSamples() % sps.getCTUSize() - (uiTPelY - uiTPelYOrg);//bottom
    if (iRemainInCU != (uiHeight >> 1))
    {
      for (UInt uiPartIdx = 0; uiPartIdx < 3; uiPartIdx++)
      {
        UInt uiHeightTmp = (uiPartIdx == 1) ? (uiHeight >> 1) :(uiHeight >> 2);
        if (uiPartIdx == 1)
        {
          uiTPelY = uiTPelY + (uiHeight >> 2);
        }
        if (uiPartIdx == 2)
        {
          uiTPelY = uiTPelY + (uiHeight >> 1);
        }
        uiBPelY = uiTPelY + uiHeightTmp - 1;
        if (uiCurTPelY >= uiTPelY && uiCurTPelY <= uiBPelY)
          xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidth, uiHeightTmp, uiFBTDepth + 1, uiDepth);
      }
    }
    else
    {
#endif
      for (UInt uiPartIdx = 0; uiPartIdx < 2; uiPartIdx++)
      {
        if (uiPartIdx == 1)
        {
          uiTPelY = uiTPelY + (uiHeight >> 1);
        }
        uiBPelY = uiTPelY + (uiHeight >> 1) - 1;
        if (uiCurTPelY >= uiTPelY && uiCurTPelY <= uiBPelY)
          xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidth, uiHeight >> 1, uiFBTDepth + 1, uiDepth);
      }
#if BT_FORCE_TT
    }
#endif
    return;
  }
  if (uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples() && (uiDepth >= iMaxQtLv)) //right
  {
#if BT_FORCE_TT
    Int iRemainInCU = getSlice()->getSPS()->getPicWidthInLumaSamples() % sps.getCTUSize() - (uiLPelX - uiLPelXOrg);//right
    if (iRemainInCU != (uiWidth >> 1))
    {
      for (UInt uiPartIdx = 0; uiPartIdx < 3; uiPartIdx++)
      {
        UInt uiWidthTmp = (uiPartIdx == 1) ? (uiWidth >> 1) : (uiWidth >> 2);
        if (uiPartIdx == 1)
        {
          uiLPelX = uiLPelX + (uiWidth >> 2);
        }
        if (uiPartIdx == 2)
        {
          uiLPelX = uiLPelX + (uiWidth >> 1);
        }
        uiRPelX = uiLPelX + uiWidthTmp - 1;
        if (uiCurLPelX >= uiLPelX && uiCurLPelX <= uiRPelX)
          xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidthTmp, uiHeight, uiFBTDepth + 1, uiDepth);
      }
    }
    else
    {
#endif
      for (UInt uiPartIdx = 0; uiPartIdx < 2; uiPartIdx++)
      {
        if (uiPartIdx == 1)
        {
          uiLPelX = uiLPelX + (uiWidth >> 1);
        }
        uiRPelX = uiLPelX + (uiWidth >> 1) - 1;
        if (uiCurLPelX >= uiLPelX && uiCurLPelX <= uiRPelX)
          xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidth >> 1, uiHeight, uiFBTDepth + 1, uiDepth);
      }
#if BT_FORCE_TT
    }
#endif
    return;
  }
  if ((uiRPelX >= sps.getPicWidthInLumaSamples() && uiBPelY >= sps.getPicHeightInLumaSamples()) || (uiDepth < iMaxQtLv)) //corner
  {
    for (UInt uiPartIdx = 0; uiPartIdx < 4; uiPartIdx++)
    {
      if (uiPartIdx == 1)
      {
        uiLPelX = uiLPelX + (uiWidth >> 1);
      }
      if (uiPartIdx == 2)
      {
        uiLPelX = uiLPelX - (uiWidth >> 1);
        uiTPelY = uiTPelY + (uiHeight >> 1);
      }
      if (uiPartIdx == 3)
      {
        uiLPelX = uiLPelX + (uiWidth >> 1);
      }
      uiRPelX = uiLPelX + (uiWidth >> 1) - 1;
      uiBPelY = uiTPelY + (uiHeight >> 1) - 1;
      if (uiCurLPelX >= uiLPelX && uiCurLPelX <= uiRPelX && uiCurTPelY >= uiTPelY && uiCurTPelY <= uiBPelY)
        xRecurseFBTD(uiAbsPartIdx, uiLPelX, uiTPelY, uiWidth >> 1, uiHeight>>1, uiFBTDepth, uiDepth+1);
    }
    return;
  }
  m_uiFBTDepth = uiFBTDepth;
}
#endif
#endif

#if RRC_SELECTABLE_IF
Int TComDataCU::getSIFModeMv(const TComMv *mv, Int& nFltTypeHor, Int& nFltTypeVer, UInt uiAbsPartIdx)
{
  nFltTypeHor = 0;
  nFltTypeVer = 0;

#if SIF_CU_LEVEL
#if RRC_SIF
#error not supported
#endif
  if (getSIFFlag(uiAbsPartIdx))
  {
#if RRC_BIF
    nFltTypeHor = 1;
    nFltTypeVer = 1;
    return SIF_BLUR;
#endif
#if RRC_MIE
    return SIF_SHARP;
#endif
  }
  else
  {
    return SIF_STD;
  }
#endif

  Int xFrac;
  Int yFrac;

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 2
  xFrac = mv->getHor() & 0x15;
  yFrac = mv->getVer() & 0x15;

#if SIF_GAUSS
  // 'cross' pattern
  Bool xFracG = xFrac >= 1 && xFrac <= 7;
  Bool yFracG = yFrac >= 1 && yFrac <= 7;
  if ((xFracG || yFracG)
#if VCEG_AZ07_FRUC_MERGE
    && getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_OFF
#endif
    )
  {
    if (xFracG)
      nFltTypeHor = 1;
    if (yFracG)
      nFltTypeVer = 1;
    return SIF_BLUR;
  }
#if SIF_GAUSS_ONLY
  else
  {
    return SIF_STD;
  }
#endif
#endif

#if SHARP_ALL_FRAC
  if (xFrac == 0 && yFrac == 0)
#else
  if (xFrac == 0 || yFrac == 0
  || (xFrac == 8 && yFrac == 8))
#endif
  {
    return SIF_STD;
  }
  assert(SIF_GAUSS_ONLY == 0);
  return SIF_SHARP;

#elif VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE == 1
  xFrac = mv->getHor() & 0x7;
  yFrac = mv->getVer() & 0x7;

#if SIF_GAUSS
  // 'cross' pattern
  Bool xFracG = xFrac >= 1 && xFrac <= 3;
  Bool yFracG = yFrac >= 1 && yFrac <= 3;
  if ((xFracG || yFracG)
#if VCEG_AZ07_FRUC_MERGE
    && getFRUCMgrMode(uiAbsPartIdx) == FRUC_MERGE_OFF
#endif
    )
  {
    if (xFracG)
      nFltTypeHor = 1;
    if (yFracG)
      nFltTypeVer = 1;
    return 1;
  }
#if SIF_GAUSS_ONLY
  else
  {
    return SIF_STD;
  }
#endif
#endif

#if SHARP_ALL_FRAC
  if (xFrac == 0 && yFrac == 0)
#else
  if (xFrac == 0 || yFrac == 0
  || (xFrac == 4 && yFrac == 4))
#endif
  {
    return SIF_STD;
  }
  assert(SIF_GAUSS_ONLY == 0);
  return SIF_SHARP;

#else
  // HEVC MV precision
  xFrac = mv->getHor() & 0x3;
  yFrac = mv->getVer() & 0x3;

#if SIF_GAUSS
  // 'cross' pattern if before xFrac==0 || yFrac == 0 check
  if (xFrac == 1 || yFrac == 1)
  {
    if (xFrac == 1)
      nFltTypeHor = 1;
    if (yFrac == 1)
      nFltTypeVer = 1;
    return SIF_BLUR;
  }
#if SIF_GAUSS_ONLY
  else
  {
    return SIF_STD;
  }
#endif
#endif

#if SHARP_ALL_FRAC
  if (xFrac == 0 && yFrac == 0)
#else
  if (xFrac == 0 || yFrac == 0 
  || (xFrac == 2 && yFrac == 2))
#endif
  {
    return SIF_STD;
  }
  assert(SIF_GAUSS_ONLY == 0);
  return SIF_SHARP;
#endif
}

Int TComDataCU::getSIFMode(const TComMv *mv, UInt uiAbsPartIdx, int check_precision, int refQP)
{
  int sifMode = SIF_STD;
  if (getSIFFlag(uiAbsPartIdx))
  {
    Int fltTypeHor, fltTypeVer;
    sifMode = getSIFModeMv(mv, fltTypeHor, fltTypeVer, uiAbsPartIdx);
  }
  return sifMode;
}

Void TComDataCU::getFilterType(const TComMv *mv, UInt uiAbsPartIdx, Int& nFltTypeHor, Int& nFltTypeVer)
{
  nFltTypeHor = 0;
  nFltTypeVer = 0;
  if (1 == getSIFMode(mv, uiAbsPartIdx, 0))
  {
    getSIFModeMv(mv, nFltTypeHor, nFltTypeVer, uiAbsPartIdx);
  }
}

Bool TComDataCU::isCodeBIF(UInt uiWidth, UInt uiDepth)
{
  if (getSlice()->isIntra())
  {
    return 0;
  }
  if (uiDepth == SIF_DEPTH)
    return 1;
  else
    return 0;
}
Bool TComDataCU::isSIFFlagCoded(UInt uiAbsPartIdx)
{
  // check if we need to decode
  if (getSlice()->isIntra()
#if SIF_CU_LEVEL
    || isIntra(uiAbsPartIdx)
#endif
#if SIF_AFFINE_HARMONIZATION
    || isAffine(uiAbsPartIdx)
#endif
#if SIF_HARMO_IMV
    || getiMVFlag(uiAbsPartIdx)
#endif
#if BIF_NO_SKIP_ONLY
    || isSkipped(uiAbsPartIdx)
#endif
    )
  {
    return false;
  }
  else
  {
    return true;
  }
}
Void TComDataCU::setSIFFlagSubParts(Char param, UInt uiAbsPartIdx, UInt uiDepth)
{
#if JVET_C0024_QTBT
  setSubPart<Char>(param, m_pSIFFlag, uiAbsPartIdx, getWidth(uiAbsPartIdx), getHeight(uiAbsPartIdx));
#else
  assert(sizeof(*m_pSIFFlag) == 1);
  memset(m_pSIFFlag + uiAbsPartIdx, param, sizeof(Char)*m_pcPic->getNumPartitionsInCtu() >> (2 * uiDepth));
#endif
}
#endif//RRC_SELECTABLE_IF
//! \}

#if BARC
static int dn_filt_taps0[4][12] = { {36, 0, -4, 0, 0, 0}, /* 3-lobe Lanczos*/
{40, 0, -12, 0, 4, 0}, /*h11*/
{26, 19, 5, -3, -4, 0 ,2},// wrong
{39, -19, -10, 4, 1}};
static int dn_filt_half_lens0[4] = {3, 5, 7,6};
static int dn_centerTap[4] = {64, 64, 64, 128};
static int dn_filt_shifts[4] = {128, 128, 128, 158};

Void  TComDataCU:: downScaleComponent2x2ZeroPhaseBoundary(const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight, Bool isAboveAvailable, Bool isLeftAvailable, ComponentID compID)
{
  const int iDestWidth = iSrcWidth/2;
  const int iDestHeight = iSrcHeight/2;

  const int HALF_NUM_DOWN_TAPS = dn_filt_half_lens0[BARC_REF_DN_FILT];
  const int* iTaps = dn_filt_taps0[BARC_REF_DN_FILT];
  const int iFilterShift = dn_filt_shifts[BARC_REF_DN_FILT];
  const int iCentreTap = dn_centerTap[BARC_REF_DN_FILT];

  const int iTempStride = BARC_MAX_CTU_BUFFER_LINE;
  Pel pTemp[BARC_MAX_CTU_BUFFER] = {0};

  int iTemp;
  Pel* tmpDest = pDest;

  // Vertical
  const Pel* pDataSrc = pSrc + iSrcStride; //select lowboundary samples

  Pel* pDataDest = pTemp + HALF_NUM_DOWN_TAPS;

  for (int j=0; j<iDestHeight; ++j) 
  {
    for (int i=0; i<iSrcWidth; ++i)
    {
      iTemp = pDataSrc[i]*iCentreTap;
      for (int k=0; k<HALF_NUM_DOWN_TAPS; ++k)
      {
        if (2*j + k + 2 >= iSrcHeight)
        {
          iTemp += iTaps[k]*(pDataSrc[i+(iSrcHeight-2-2*j)*iSrcStride] + pDataSrc[i-(k+1)*iSrcStride]);
        }
        else if (2*j - k < 0 && !isAboveAvailable)
        {
          iTemp += iTaps[k]*(pDataSrc[i+(k+1)*iSrcStride] + pDataSrc[i-(2*j)*iSrcStride]);
        }
        else
          iTemp += iTaps[k]*(pDataSrc[i+(k+1)*iSrcStride] + pDataSrc[i-(k+1)*iSrcStride]);
      }
      iTemp += (iFilterShift>>1); 
      pDataDest[i] = ClipA(iTemp/iFilterShift, compID);
    }
    if (isLeftAvailable)
    {
      for (int i = -HALF_NUM_DOWN_TAPS; i < 0; i++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
      std::fill(pDataDest - HALF_NUM_DOWN_TAPS, pDataDest, pDataDest[0] );
    std::fill(pDataDest + iSrcWidth, pDataDest+iSrcWidth+HALF_NUM_DOWN_TAPS, pDataDest[iSrcWidth-1] );
    pDataSrc += 2*iSrcStride;
    pDataDest += iTempStride;
  }

// Horizontal
  pDataSrc = pTemp + HALF_NUM_DOWN_TAPS;
  pDataDest = pDest;

  for (int j=0; j<iDestHeight; ++j)
  {
    for (int i=0,ii=1; i<iDestWidth; ++i,ii+=2) //select right samples
    {
      iTemp = pDataSrc[ii]*iCentreTap;
      for (int k=0; k<HALF_NUM_DOWN_TAPS; ++k)
      {
        iTemp += iTaps[k]*(pDataSrc[ii+k+1] + pDataSrc[ii-k-1]);
      }
      iTemp += (iFilterShift>>1);
      pDataDest[i] = ClipA(iTemp/iFilterShift, compID);
    }
    pDataSrc += iTempStride;
    pDataDest += iDestStride;
  }

  Int left, middle, right;

  pDataDest = pDest + iDestStride*(iDestHeight-1);
  pDataSrc = pSrc + iSrcStride*(iSrcHeight-1); //select lowboundary samples
  for (Int i = 0; i < iDestWidth; i++)
  {
    left = 2*i; 
    middle = 2*i+1;
    right = (2*i+2 >= iSrcWidth) ? (2*i+1) : (2*i+2);
    pDataDest[i] = ClipA((pDataSrc[left] + 2*pDataSrc[middle] + pDataSrc[right])>>2, compID);;
  }

  pDataDest = pDest + iDestWidth-1;
  pDataSrc =  pSrc + iSrcWidth-1; //select lowboundary samples
  for (Int i = 0; i < iDestHeight; i++)
  {
    left = 2*i*iSrcStride; 
    middle = (2*i+1)*iSrcStride;
    right = (2*i+2 >= iSrcHeight) ? (2*i+1)*iSrcStride : (2*i+2)*iSrcStride;
    pDataDest[i*iDestStride] = ClipA((pDataSrc[left] + 2*pDataSrc[middle] + pDataSrc[right])>>2, compID);
  }
} 

static int up_filt_taps0[4][4] = { {40, -11, 4,  -1},            /*f7*/
{36, -4, 0, 0},           /*Linear*/
{20, -5, 1, 0},
{20, -4, 0, 0} };
static int up_filt_half_lens0[4] = {4, 2, 3, 2};
static int up_filt_shifts0[4] = {6, 6, 5, 5};

/** Separable upscale 2x2 with samples at phases 0 and +1/2 with respect to low-res grid
* \returns void
*/
Void  TComDataCU:: upScaleComponent2x2( const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight, Bool isAboveAvail, Bool isLeftAvail, ComponentID compID)
{
  const int iDestHeight = iSrcHeight*2;

  UInt iFitlterSelect = compID == COMPONENT_Y ? 0 : 1;
  const int NUM_UP_TAPS = up_filt_half_lens0[BARC_REF_UP_FILT + iFitlterSelect];
  const int* iTaps = up_filt_taps0[BARC_REF_UP_FILT + iFitlterSelect];
  const int iFilterShift = up_filt_shifts0[BARC_REF_UP_FILT + iFitlterSelect];
#if LIYE_FIX_BUG
  const int iTempStride = iSrcWidth + 2 * NUM_UP_TAPS;
  Pel* pTemp = new Pel[iDestHeight*iTempStride];
#else
  const int iTempStride = BARC_MAX_CTU_BUFFER_LINE;
  Pel pTemp[BARC_MAX_CTU_BUFFER] = {0};
#endif
  // Vertical
  const Pel* pDataSrc = pSrc;

  Pel* pDataDest = pTemp + NUM_UP_TAPS;
  int iTemp;

  for (int j=0; j<iSrcHeight; ++j)
  {
    // +1/2 phase
    for (int i=0; i<iSrcWidth; ++i)
    {
      iTemp = 0;
      for (int k=0; k<NUM_UP_TAPS; ++k)
      {
        if ( j - k - 1 < 0 && !isAboveAvail )
        {
          iTemp += iTaps[k]*(pDataSrc[i-j*iSrcStride] + pDataSrc[i+k*iSrcStride]);
        }
        else if ((j +  k ) >= iSrcHeight)
          iTemp += iTaps[k]*(pDataSrc[i-(k+1)*iSrcStride] + pDataSrc[i+(iSrcHeight-j-1)*iSrcStride]);
        else
          iTemp += iTaps[k]*(pDataSrc[i-(k+1)*iSrcStride] + pDataSrc[i+k*iSrcStride]);
      }
      iTemp += 1<<(iFilterShift-1);
      pDataDest[i] = ClipA(iTemp>>iFilterShift, compID);
    }

    if (isLeftAvail)
    {
      for (int i = -NUM_UP_TAPS; i < 0; i ++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else 
    {
      std::fill(pDataDest - NUM_UP_TAPS, pDataDest, pDataDest[0] );
    }
    std::fill(pDataDest + iSrcWidth, pDataDest+iSrcWidth+NUM_UP_TAPS, pDataDest[iSrcWidth-1] ); //fill boundary samples

    pDataDest += iTempStride;

    // Zero phase
    for (int i=0; i<iSrcWidth; ++i)
    {
      pDataDest[i] = pDataSrc[i]; //fill boundary samples
    }
    if (isLeftAvail)
    {
      for (int i = -NUM_UP_TAPS; i < 0; i ++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
    {
      std::fill(pDataDest - NUM_UP_TAPS, pDataDest, pDataDest[0] );
    }
    std::fill(pDataDest + iSrcWidth, pDataDest+iSrcWidth+NUM_UP_TAPS, pDataDest[iSrcWidth-1] ); //fill boundary samples

    pDataDest += iTempStride;
    pDataSrc += iSrcStride;
  }

  pDataSrc = pTemp + NUM_UP_TAPS;
  pDataDest = pDest;

  for (int j=0; j<iDestHeight; ++j)
  {
    for (int i=0,ii=0; i<iSrcWidth; ++i, ii+=2)
    {
      pDataDest[ii+1] = pDataSrc[i];

      iTemp = 0;
      for (int k=0; k<NUM_UP_TAPS; ++k)
      {
        iTemp += iTaps[k] * (pDataSrc[i + k] + pDataSrc[i -(k+1)] );
      }
      iTemp += 1<<(iFilterShift-1);
      pDataDest[ii] = ClipA(iTemp>>iFilterShift, compID);
    }
    pDataSrc += iTempStride;
    pDataDest += iDestStride;
  }

  pDest[iSrcWidth*2-1] = pDest[(iSrcWidth*2-1)+iDestStride];
}
Void  TComDataCU::upScaleComponent2x2Post(const Pel* pSrc, Pel* pDest, const Int iSrcStride, const Int iDestStride, const Int iSrcWidth, const Int iSrcHeight, Bool isAboveAvail, Bool isLeftAvail, Bool isBottomAvail, Bool isRightAvail, ComponentID compID)
{
  const int iDestHeight = iSrcHeight * 2;
  UInt iFitlterSelect = compID == COMPONENT_Y ? 0 : 1;
  const int NUM_UP_TAPS = up_filt_half_lens0[BARC_REF_UP_FILT + iFitlterSelect];
  const int* iTaps = up_filt_taps0[BARC_REF_UP_FILT + iFitlterSelect];
  const int iFilterShift = up_filt_shifts0[BARC_REF_UP_FILT + iFitlterSelect];

  const int iTempStride = BARC_MAX_CTU_BUFFER_LINE;
  Pel pTemp[BARC_MAX_CTU_BUFFER] = {0};

  // Vertical
  const Pel* pDataSrc = pSrc;

  Pel* pDataDest = pTemp + NUM_UP_TAPS;
  int iTemp;

  for (int j = 0; j < iSrcHeight; ++j)
  {
    // +1/2 phase
    for (int i = 0; i < iSrcWidth; ++i)
    {
      iTemp = 0;
      for (int k = 0; k < NUM_UP_TAPS; ++k)
      {
        if (j - k - 1 < 0 && !isAboveAvail)
        {
          iTemp += iTaps[k] * (pDataSrc[i - j*iSrcStride] + pDataSrc[i + k*iSrcStride]);
        }
        else if ((j + k) >= iSrcHeight && !isBottomAvail)
          iTemp += iTaps[k] * (pDataSrc[i - (k + 1)*iSrcStride] + pDataSrc[i + (iSrcHeight - j - 1)*iSrcStride]);
        else
          iTemp += iTaps[k] * (pDataSrc[i - (k + 1)*iSrcStride] + pDataSrc[i + k*iSrcStride]);
      }
      iTemp += 1 << (iFilterShift - 1);
      pDataDest[i] = ClipA(iTemp >> iFilterShift, compID);
    }

    if (isLeftAvail)
    {
      for (int i = -NUM_UP_TAPS; i < 0; i++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
      std::fill(pDataDest - NUM_UP_TAPS, pDataDest, pDataDest[0]);

    if (isRightAvail)
    {
      for (int i = iSrcWidth; i < iSrcWidth + NUM_UP_TAPS; i++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
      std::fill(pDataDest + iSrcWidth, pDataDest + iSrcWidth + NUM_UP_TAPS, pDataDest[iSrcWidth - 1]); //fill boundary samples

    pDataDest += iTempStride;

    // Zero phase
    for (int i = 0; i < iSrcWidth; ++i)
    {
      pDataDest[i] = pDataSrc[i]; //fill boundary samples
    }
    if (isLeftAvail)
    {
      for (int i = -NUM_UP_TAPS; i < 0; i++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
      std::fill(pDataDest - NUM_UP_TAPS, pDataDest, pDataDest[0]);

    if (isRightAvail)
    {
      for (int i = iSrcWidth; i < iSrcWidth + NUM_UP_TAPS; i++)
      {
        pDataDest[i] = pDataSrc[i];
      }
    }
    else
      std::fill(pDataDest + iSrcWidth, pDataDest + iSrcWidth + NUM_UP_TAPS, pDataDest[iSrcWidth - 1]); //fill boundary samples

    pDataDest += iTempStride;
    pDataSrc += iSrcStride;
  }

  pDataSrc = pTemp + NUM_UP_TAPS;
  pDataDest = pDest;

  for (int j = 0; j < iDestHeight; ++j)
  {
    for (int i = 0, ii = 0; i < iSrcWidth; ++i, ii += 2)
    {
      if ((j >= iDestHeight - 2*NUM_UP_TAPS) || ( i >= iSrcWidth - NUM_UP_TAPS))
      {
        pDataDest[ii + 1] = pDataSrc[i];

        iTemp = 0;
        for (int k = 0; k < NUM_UP_TAPS; ++k)
        {
          iTemp += iTaps[k] * (pDataSrc[i + k] + pDataSrc[i - (k + 1)]);
        }
        iTemp += 1 << (iFilterShift - 1);
        pDataDest[ii] = ClipA(iTemp >> iFilterShift, compID);
      }
    }
    pDataSrc += iTempStride;
    pDataDest += iDestStride;
  }

  pDest[iSrcWidth * 2 - 1] = pDest[(iSrcWidth * 2 - 1) + iDestStride];
}
#endif
//! \}

