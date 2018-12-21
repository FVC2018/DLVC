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

/** \file     TEncSearch.cpp
\brief    encoder search class
*/

#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComRom.h"
#include "TLibCommon/TComMotionInfo.h"
#include "TEncSearch.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/Debug.h"
#if RRC_SELECTABLE_IF
#include "TLibCommon/TComSharpeningFilter.h"
#endif
#include <math.h>
#include <limits>


#if AFF_SKIP_ME
int g_FourDir = -1;
int g_FourRefIdx0 = -1;
int g_FourRefIdx1 = -1;
int g_FourRefIdxBi0 = -1;
int g_FourRefIdxBi1 = -1;
//Int g_HevcDir = -1;
//Int g_HevcRefIdx0 = -1;
//Int g_HevcRefIdx1 = -1;
//Int g_HevcRefIdxBi0 = -1;
//Int g_HevcRefIdxBi1 = -1;
#endif
//! \ingroup TLibEncoder
//! \{

#if VCEG_AZ08_KLT_COMMON
extern short **g_ppsEigenVector[USE_MORE_BLOCKSIZE_DEPTH_MAX];
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING_ENC_METHOD
namespace {
#if !!ITH_C1_SMT_RESIDUE
	template <int N>
	Int average(const Pel *r, Int i0, Int j0, Int uiHeight, Int uiWidth, UInt uiStrideRes)
	{
		Int s = 0, n = 0;
		for (int i = -N / 2 + i0; i <= N / 2 + i0; ++i)
			for (int j = -N / 2 + j0; j <= N / 2 + j0; ++j)
			{
				const int ii = (i + uiHeight) % uiHeight;
				const int jj = (j + uiWidth) % uiWidth;
				s += r[ii*uiStrideRes + jj];
				++n;
			}
		assert(n > 0);
		return (s + n / 2) / n;
	}

	void smoothResidual(const Bound prm, const Pel *piOrg, Pel *piResi, UInt uiHeight, UInt uiWidth, UInt uiStrideOrg, UInt uiStrideRes)
	{
		// find boundaries of the res

		char bmM[MAX_CU_SIZE * MAX_CU_SIZE];
		Pel r[MAX_CU_SIZE * MAX_CU_SIZE];

		for (int i = 0, k = 0; i < (int)uiHeight; ++i)
		{
			for (int j = 0; j < (int)uiWidth; ++j, k++)
			{
				if (piOrg[i*uiStrideOrg + j] <= prm.m)      bmM[k] = -1;
				else if (piOrg[i*uiStrideOrg + j] >= prm.M) bmM[k] = +1;
				else                                        bmM[k] = 0;
			}
		}

		twoDimBufCopier((UChar *)r, uiWidth * sizeof(Pel), (UChar *)piResi, uiStrideRes * sizeof(Pel), uiWidth * sizeof(Pel), uiHeight);

		const int cptmax = 4;
		int cpt = 0;
		bool cont = true;
		while (cpt < cptmax&&cont)
		{
			cont = false;
			for (int i = 0; i < (int)uiHeight; ++i)
				for (int j = 0; j < (int)uiWidth; ++j)
				{
					if (bmM[i*uiWidth + j] == -1)
					{ // we can lower the res
						Int s = average<3>(piResi, i, j, uiHeight, uiWidth, uiStrideRes);
						if (s < piResi[i * uiStrideRes + j]) r[i * uiWidth + j] = s;
						cont = true;
					}
					else if (bmM[i*uiWidth + j] == 1)
					{ // we can inc the res
						Int s = average<3>(piResi, i, j, uiHeight, uiWidth, uiStrideRes);
						if (s > piResi[i * uiStrideRes + j]) r[i * uiWidth + j] = s;
						cont = true;
					}
				}
			++cpt;

			twoDimBufCopier((UChar *)piResi, uiStrideRes * sizeof(Pel), (UChar *)r, uiWidth * sizeof(Pel), uiWidth * sizeof(Pel), uiHeight);
		}
	}

	void smoothResidual(const ComponentID compID, const Pel *piOrg, Pel *piResi, UInt uiHeight, UInt uiWidth, UInt uiStrideOrg, UInt uiStrideRes)
	{
#if !!ITH_C1_CLIPA
		Bound prm = g_ClipParam.getComponentBound(compID);
#else
		Bound prm;

		switch (compID)
		{
		case COMPONENT_Y:    prm = g_ClipParam.Y(); break;
		case COMPONENT_Cb:   prm = g_ClipParam.U(); break;
		case COMPONENT_Cr:   prm = g_ClipParam.V(); break;
		default: assert(false);
		}
#endif

		Bool activate = false;
		for (Int i = 0; i < uiHeight; ++i)
		{
			for (Int j = 0; j < uiWidth; ++j)
			{
				if (piOrg[i*uiStrideOrg + j] <= prm.m || piOrg[i*uiStrideOrg + j] >= prm.M)
				{
					activate = true;
					break;
				}
			}
			if (activate == true)
			{
				break;
			}
		}

		if (activate)
		{
			smoothResidual(prm, piOrg, piResi, uiHeight, uiWidth, uiStrideOrg, uiStrideRes);
		}
	}
#else
	template <int N>
	Int average(const std::vector<Pel> &r, Int i0, Int j0, Int uiHeight, Int uiWidth) {
		Int s = 0, n = 0;
		for (int i = -N / 2 + i0; i <= N / 2 + i0; ++i)
			for (int j = -N / 2 + j0; j <= N / 2 + j0; ++j) {
				const int ii = (i + uiHeight) % uiHeight;
				const int jj = (j + uiWidth) % uiWidth;
				s += r[ii*uiWidth + jj];
				++n;
			}
		assert(n>0);
		return (s + n / 2) / n;
	}

	void smoothResidual(const Bound prm, const std::vector<Pel> &org, std::vector<Pel> &res, UInt uiHeight, UInt uiWidth) {
		// find boundaries of the res
		static std::vector<char> bmM;
		static std::vector<Pel> r;
		r = res;
		bmM.resize(uiHeight*uiWidth); // avoir realloc

		for (int k = 0; k<(int)org.size(); ++k) {
			if (org[k] <= prm.m)      bmM[k] = -1;
			else if (org[k] >= prm.M) bmM[k] = +1;
			else                    bmM[k] = 0;
		}

		const int cptmax = 4;
		int cpt = 0;
		bool cont = true;
		while (cpt<cptmax&&cont) {
			cont = false;
			for (int i = 0; i<(int)uiHeight; ++i)
				for (int j = 0; j<(int)uiWidth; ++j) {
					if (bmM[i*uiWidth + j] == -1) { // we can lower the res
						Int s = average<3>(res, i, j, uiHeight, uiWidth);
						if (s<res[i*uiWidth + j]) r[i*uiWidth + j] = s;
						cont = true;
					}
					else if (bmM[i*uiWidth + j] == 1) { // we can inc the res
						Int s = average<3>(res, i, j, uiHeight, uiWidth);
						if (s>res[i*uiWidth + j]) r[i*uiWidth + j] = s;
						cont = true;
					}
				}
			++cpt;
			std::copy(r.begin(), r.end(), res.begin());
		}
	}

	void smoothResidual(const ComponentID compID, const Pel *piOrg, Pel *piResi, UInt uiHeight, UInt uiWidth, UInt uiStrideOrg, UInt uiStrideRes) {

		Bound prm;
		switch (compID) {
		case COMPONENT_Y:    prm = g_ClipParam.Y(); break;
		case COMPONENT_Cb:   prm = g_ClipParam.U(); break;
		case COMPONENT_Cr:   prm = g_ClipParam.V(); break;
		default: assert(false);
		}

		static std::vector<Pel> org; // avoid realloc
		org.resize(uiHeight*uiWidth);

		Bool activate = false;
		for (Int i = 0, k = 0; i<uiHeight; ++i)
			for (Int j = 0; j<uiWidth; ++j, ++k) {
				org[k] = piOrg[i*uiStrideOrg + j];
				if (org[k] <= prm.m || org[k] >= prm.M) activate = true;
			}


		if (activate) {
			static std::vector<Pel> r; // avoid realloc
			r.resize(uiHeight*uiWidth);
			for (Int i = 0, k = 0; i<uiHeight; ++i)
				for (Int j = 0; j<uiWidth; ++j, ++k) {
					r[k] = piResi[i*uiStrideRes + j];
				}
			smoothResidual(prm, org, r, uiHeight, uiWidth);
			// copy back
			for (Int i = 0, k = 0; i<uiHeight; ++i)
				for (Int j = 0; j<uiWidth; ++j, ++k) {
					piResi[i*uiStrideRes + j] = r[k];
				}
		}
	}
#endif
	void smoothResidual(const TComYuv* pcYuvOrg, TComYuv *pcYuvRes, const UInt uiTrUnitIdx, const UInt uiHeight, const UInt uiWidth) {

		for (Int comp = 0; comp<pcYuvRes->getNumberValidComponents(); ++comp) {
			const ComponentID compID = ComponentID(comp);

			const Int uiPartWidth = (uiWidth >> pcYuvRes->getComponentScaleX(compID));
			const Int uiPartHeight = (uiHeight >> pcYuvRes->getComponentScaleY(compID));

			const Pel* pOrg = pcYuvOrg->getAddr(compID, uiTrUnitIdx, uiPartWidth);
			const Int  iOrgStride = pcYuvOrg->getStride(compID);

			Pel* pRes = pcYuvRes->getAddr(compID, uiTrUnitIdx, uiPartWidth);
			const Int  iResStride = pcYuvRes->getStride(compID);

			smoothResidual(compID, pOrg, pRes, uiPartHeight, uiPartWidth, iOrgStride, iResStride);
		}
	}

} // anom
#endif
static const TComMv s_acMvRefineH[9] =
{
	TComMv(0,  0), // 0
	TComMv(0, -1), // 1
	TComMv(0,  1), // 2
	TComMv(-1,  0), // 3
	TComMv(1,  0), // 4
	TComMv(-1, -1), // 5
	TComMv(1, -1), // 6
	TComMv(-1,  1), // 7
	TComMv(1,  1)  // 8
};

static const TComMv s_acMvRefineQ[9] =
{
	TComMv(0,  0), // 0
	TComMv(0, -1), // 1
	TComMv(0,  1), // 2
	TComMv(-1, -1), // 5
	TComMv(1, -1), // 6
	TComMv(-1,  0), // 3
	TComMv(1,  0), // 4
	TComMv(-1,  1), // 7
	TComMv(1,  1)  // 8
};

static const UInt s_auiDFilter[9] =
{
	0, 1, 0,
	2, 3, 2,
	0, 1, 0
};

#if !JVET_C0024_QTBT
static Void offsetSubTUCBFs(TComTU &rTu, const ComponentID compID)
{
	TComDataCU *pcCU = rTu.getCU();
	const UInt        uiTrDepth = rTu.GetTransformDepthRel();
	const UInt        uiAbsPartIdx = rTu.GetAbsPartIdxTU(compID);
	const UInt        partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

	//move the CBFs down a level and set the parent CBF

	UChar subTUCBF[2];
	UChar combinedSubTUCBF = 0;

	for (UInt subTU = 0; subTU < 2; subTU++)
	{
		const UInt subTUAbsPartIdx = uiAbsPartIdx + (subTU * partIdxesPerSubTU);

		subTUCBF[subTU] = pcCU->getCbf(subTUAbsPartIdx, compID, uiTrDepth);
		combinedSubTUCBF |= subTUCBF[subTU];
	}

	for (UInt subTU = 0; subTU < 2; subTU++)
	{
		const UInt subTUAbsPartIdx = uiAbsPartIdx + (subTU * partIdxesPerSubTU);
		const UChar compositeCBF = (subTUCBF[subTU] << 1) | combinedSubTUCBF;

		pcCU->setCbfPartRange((compositeCBF << uiTrDepth), compID, subTUAbsPartIdx, partIdxesPerSubTU);
	}
}
#endif

TEncSearch::TEncSearch()
#if JVET_C0024_QTBT
	: m_ppcQTTempTComYuv(NULL)
#else
	: m_puhQTTempTrIdx(NULL)
	, m_pcQTTempTComYuv(NULL)
#endif
#if VCEG_AZ08_INTER_KLT
#if JVET_C0024_QTBT
	, m_ppcQTTempTComYuvRec(NULL)
#else
	, m_pcQTTempTComYuvRec(NULL)
#endif
#endif
	, m_pcEncCfg(NULL)
	, m_pcTrQuant(NULL)
	, m_pcRdCost(NULL)
	, m_pcEntropyCoder(NULL)
	, m_iSearchRange(0)
	, m_bipredSearchRange(0)
	, m_iFastSearch(0)
#if JVET_C0024_QTBT
	, m_ppppcRDSbacCoder(NULL)
#else
	, m_pppcRDSbacCoder(NULL)
#endif
	, m_pcRDGoOnSbacCoder(NULL)
	, m_pTempPel(NULL)
	, m_puiDFilter(NULL)
	, m_isInitialized(false)
#if COM16_C806_EMT
	, m_puhQTTempEmtTuIdx(NULL)
	, m_puhQTTempEmtCuFlag(NULL)
#endif
{
#if JVET_D0123_ME_CTX_LUT_BITS
	m_pcPuMeEstBitsSbac = new estPuMeBitsSbacStruct;
#endif
	for (UInt ch = 0; ch<MAX_NUM_COMPONENT; ch++)
	{
#if JVET_C0024_QTBT
		m_pppcQTTempCoeff[ch] = NULL;
#else
		m_ppcQTTempCoeff[ch] = NULL;
#endif
		m_pcQTTempCoeff[ch] = NULL;
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
		m_pppcQTTempArlCoeff[ch] = NULL;
#else
		m_ppcQTTempArlCoeff[ch] = NULL;
#endif
		m_pcQTTempArlCoeff[ch] = NULL;
#endif
		m_puhQTTempCbf[ch] = NULL;
		m_phQTTempCrossComponentPredictionAlpha[ch] = NULL;
		m_pSharedPredTransformSkip[ch] = NULL;
		m_pcQTTempTUCoeff[ch] = NULL;
#if ADAPTIVE_QP_SELECTION
		m_ppcQTTempTUArlCoeff[ch] = NULL;
#endif
		m_puhQTTempTransformSkipFlag[ch] = NULL;
#if VCEG_AZ08_KLT_COMMON
		m_puhQTTempKLTFlag[ch] = NULL;
#endif
#if COM16_C806_EMT
		m_puhQTTempExplicitRdpcmMode[ch] = NULL;
#endif
	}

	for (Int i = 0; i<MAX_NUM_REF_LIST_ADAPT_SR; i++)
	{
		memset(m_aaiAdaptSR[i], 0, MAX_IDX_ADAPT_SR * sizeof(Int));
	}
	for (Int i = 0; i<AMVP_MAX_NUM_CANDS + 1; i++)
	{
		memset(m_auiMVPIdxCost[i], 0, (AMVP_MAX_NUM_CANDS + 1) * sizeof(UInt));
	}

#if COM16_C806_LARGE_CTU
#if !JVET_C0024_QTBT
	memset(m_resiSplitBuffer, 0, sizeof(m_resiSplitBuffer));
#endif
	memset(m_resiSingleBuffer, 0, sizeof(m_resiSingleBuffer));
	memset(m_resiPUBuffer, 0, sizeof(m_resiPUBuffer));
#endif

	setWpScalingDistParam(NULL, -1, REF_PIC_LIST_X);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
	for (Int i = 0; i< NUM_MGR_TYPE; i++)
	{
		m_pMvFieldSP[i] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH * 2];
		m_phInterDirSP[i] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
		assert(m_pMvFieldSP[i] != NULL && m_phInterDirSP[i] != NULL);
	}
#else
	m_pMvFieldSP[0] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH * 2];
	m_pMvFieldSP[1] = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH * 2];
	m_phInterDirSP[0] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	m_phInterDirSP[1] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	assert(m_pMvFieldSP[0] != NULL && m_phInterDirSP[0] != NULL);
	assert(m_pMvFieldSP[1] != NULL && m_phInterDirSP[1] != NULL);
#endif
#endif

#if VCEG_AZ07_FRUC_MERGE
	m_pMvFieldFRUC = new TComMvField[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH * 2];
	m_phInterDirFRUC = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	m_phFRUCRefineDist[0] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	m_phFRUCRefineDist[1] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	m_phFRUCSBlkRefineDist[0] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	m_phFRUCSBlkRefineDist[1] = new UChar[MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH];
	assert(m_pMvFieldFRUC != NULL && m_phInterDirFRUC != NULL);
#endif

#if   JVET_E0077_LM_MF
	m_pLMMFPredSaved = new Pel*[8];//4*(Cb+Cr)
	for (Int k = 0; k < 8; k++)
	{
		m_pLMMFPredSaved[k] = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
	}
#endif
}

#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
Void  TEncSearch::cleanSaveLoadInfo(UInt uiPartIdx, Int uiWIdx, Int uiHIdx)
{
	m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx][uiPartIdx] = 0;
}

UChar TEncSearch::searchSaveLoadIdx(UInt uiPartIdx, UInt uiWIdx, UInt uiHIdx, UInt uiBTDepth, UInt* puiCURouteInfo)
{
	if (uiBTDepth == 0)
	{
		//not need to save for QTLeaf node
		m_SaveLoadTag[uiWIdx][uiHIdx][uiPartIdx] = SAVE_ENC_INFO;
		m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx][uiPartIdx]++;
		return 0;
	}

	UInt idx, i;
	UChar ucRouteIdx = MAX_NUM_ROUTE;
	UInt uiCodingOrder = puiCURouteInfo[0];
	UInt uiLeftNeighorInfo = puiCURouteInfo[1];
	UInt uiUpNeighborInfo = puiCURouteInfo[2];
	UInt uiNumRouteSaved = m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx][uiPartIdx];
#if BT_FORCE_EXT
	assert(uiNumRouteSaved < MAX_NUM_ROUTE);
#endif
	//search the matching RouteIdx
	for (i = 0; i<uiNumRouteSaved; i++)
	{
		idx = uiPartIdx + (i * MAX_PART_IN_CTU);
		if (uiCodingOrder == m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx][idx])
		{
			if (uiLeftNeighorInfo == m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx][idx])
			{
				if (uiUpNeighborInfo == m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx][idx])
					ucRouteIdx = i;
			}
		}
	}

#if FAST_ADDITIONAL_REF_ROUTE //if not find matching route in the restrict pass, search in a second pass
	//R5 match: only match coding order + two neighbor size (does not consider BTDepth)
	if (ucRouteIdx == MAX_NUM_ROUTE)
	{
		for (i = 0; i<uiNumRouteSaved; i++)
		{
			idx = uiPartIdx + (i * MAX_PART_IN_CTU);
			if (uiCodingOrder == m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx][idx])
			{
				if ((uiLeftNeighorInfo >> 2) == (m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx][idx] >> 2))
				{
					if ((uiUpNeighborInfo >> 2) == (m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx][idx] >> 2))
						ucRouteIdx = i;
				}
			}
		}
	}
#endif

	if (ucRouteIdx == MAX_NUM_ROUTE) //not match any existing route
	{
		ucRouteIdx = m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx][uiPartIdx];
		//save RouteInfo to the new node
		idx = ucRouteIdx*MAX_PART_IN_CTU;
		m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx][uiPartIdx + idx] = uiCodingOrder;
		m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx][uiPartIdx + idx] = uiLeftNeighorInfo;
		m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx][uiPartIdx + idx] = uiUpNeighborInfo;
		m_SaveLoadTag[uiWIdx][uiHIdx][uiPartIdx + idx] = SAVE_ENC_INFO;
		m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx][uiPartIdx]++;
	}

	return ucRouteIdx;
}
#endif

Void TEncSearch::destroy()
{
	assert(m_isInitialized);
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
	for (UInt uiWIdx = 0; uiWIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiWIdx++)
	{
		for (UInt uiHIdx = 0; uiHIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiHIdx++)
		{
			if (m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx])  free(m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx]);
			if (m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx])  free(m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx]);
			if (m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx])  free(m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx]);
			if (m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx])  free(m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx]);
			if (m_SaveLoadTag[uiWIdx][uiHIdx])  free(m_SaveLoadTag[uiWIdx][uiHIdx]);
#if COM16_C806_EMT
			if (m_SaveLoadEmtFlag[uiWIdx][uiHIdx])  free(m_SaveLoadEmtFlag[uiWIdx][uiHIdx]);
			if (m_SaveLoadEmtIdx[uiWIdx][uiHIdx])  free(m_SaveLoadEmtIdx[uiWIdx][uiHIdx]);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
			if (m_SaveLoadRotIdx[uiWIdx][uiHIdx])  free(m_SaveLoadRotIdx[uiWIdx][uiHIdx]);
#endif
#if COM16_C1046_PDPC_INTRA
			if (m_SaveLoadPdpcIdx[uiWIdx][uiHIdx])  free(m_SaveLoadPdpcIdx[uiWIdx][uiHIdx]);
#endif
#if VCEG_AZ07_FRUC_MERGE
			if (m_SaveLoadFrucMode[uiWIdx][uiHIdx])  free(m_SaveLoadFrucMode[uiWIdx][uiHIdx]);
#endif
#if JVET_E0076_MULTI_PEL_MVD
			if (m_SaveLoadIMVFlag[uiWIdx][uiHIdx])  free(m_SaveLoadIMVFlag[uiWIdx][uiHIdx]);
#endif
#if VCEG_AZ06_IC
			if (m_SaveLoadICFlag[uiWIdx][uiHIdx])  free(m_SaveLoadICFlag[uiWIdx][uiHIdx]);
#endif
#if COM16_C1016_AFFINE
			if (m_SaveLoadAffineFlag[uiWIdx][uiHIdx])  free(m_SaveLoadAffineFlag[uiWIdx][uiHIdx]);
#endif
			if (m_SaveLoadMergeFlag[uiWIdx][uiHIdx])  free(m_SaveLoadMergeFlag[uiWIdx][uiHIdx]);
			if (m_SaveLoadInterDir[uiWIdx][uiHIdx])  free(m_SaveLoadInterDir[uiWIdx][uiHIdx]);
			if (m_SaveLoadSplit[uiWIdx][uiHIdx])  free(m_SaveLoadSplit[uiWIdx][uiHIdx]);
		}
	}
#endif
	if (m_pTempPel)
	{
		delete[] m_pTempPel;
		m_pTempPel = NULL;
	}

#if JVET_C0024_QTBT
	const UInt uiNumLayersAllocated = g_aucConvertToBit[m_pcEncCfg->getCTUSize()] + 1;
#endif
	if (m_pcEncCfg)
	{
#if !JVET_C0024_QTBT
		const UInt uiNumLayersAllocated = m_pcEncCfg->getQuadtreeTULog2MaxSize() - m_pcEncCfg->getQuadtreeTULog2MinSize() + 1;
#endif

		for (UInt ch = 0; ch<MAX_NUM_COMPONENT; ch++)
		{
#if JVET_C0024_QTBT
			for (UInt uiWIdx = 0; uiWIdx < uiNumLayersAllocated; uiWIdx++)
			{
				for (UInt uiHIdx = 0; uiHIdx < uiNumLayersAllocated; uiHIdx++)
				{
					delete[] m_pppcQTTempCoeff[ch][uiWIdx][uiHIdx];
#if ADAPTIVE_QP_SELECTION
					delete[] m_pppcQTTempArlCoeff[ch][uiWIdx][uiHIdx];
#endif
				}
				delete[] m_pppcQTTempCoeff[ch][uiWIdx];
#if ADAPTIVE_QP_SELECTION
				delete[] m_pppcQTTempArlCoeff[ch][uiWIdx];
#endif
			}
			delete[] m_pppcQTTempCoeff[ch];
#if ADAPTIVE_QP_SELECTION
			delete[] m_pppcQTTempArlCoeff[ch];
#endif
#else
			for (UInt layer = 0; layer < uiNumLayersAllocated; layer++)
			{
				delete[] m_ppcQTTempCoeff[ch][layer];
#if ADAPTIVE_QP_SELECTION
				delete[] m_ppcQTTempArlCoeff[ch][layer];
#endif
			}
			delete[] m_ppcQTTempCoeff[ch];
#endif
			delete[] m_pcQTTempCoeff[ch];
			delete[] m_puhQTTempCbf[ch];
#if ADAPTIVE_QP_SELECTION
#if !JVET_C0024_QTBT
			delete[] m_ppcQTTempArlCoeff[ch];
#endif
			delete[] m_pcQTTempArlCoeff[ch];
#endif
		}

#if JVET_C0024_QTBT
		for (UInt uiWIdx = 0; uiWIdx < uiNumLayersAllocated; uiWIdx++)
		{
			for (UInt uiHIdx = 0; uiHIdx < uiNumLayersAllocated; uiHIdx++)
			{
				m_ppcQTTempTComYuv[uiWIdx][uiHIdx].destroy();
#if JVET_C0024_QTBT
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
				if (m_pcTrQuant->getUseInterKLT())
				{
#endif
					m_ppcQTTempTComYuvRec[uiWIdx][uiHIdx].destroy();
#if VCEG_AZ08_USE_KLT
				}
#endif
#endif
#endif
			}
			delete[] m_ppcQTTempTComYuv[uiWIdx];
#if JVET_C0024_QTBT && VCEG_AZ08_INTER_KLT
			if (m_pcTrQuant->getUseInterKLT())
			{
				delete[] m_ppcQTTempTComYuvRec[uiWIdx];
			}
#endif
		}
		delete[] m_ppcQTTempTComYuv;
#if JVET_C0024_QTBT && VCEG_AZ08_INTER_KLT
		if (m_pcTrQuant->getUseInterKLT())
		{
			delete[] m_ppcQTTempTComYuvRec;
		}
#endif
#else
		for (UInt layer = 0; layer < uiNumLayersAllocated; layer++)
		{
			m_pcQTTempTComYuv[layer].destroy();
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
			if (m_pcTrQuant->getUseInterKLT())
			{
#endif
				m_pcQTTempTComYuvRec[layer].destroy();
#if VCEG_AZ08_USE_KLT
			}
#endif
#endif
		}
#endif
	}
#if !JVET_C0024_QTBT
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
	if (m_pcTrQuant->getUseInterKLT())
	{
#endif
		delete[] m_pcQTTempTComYuvRec;
#if VCEG_AZ08_USE_KLT
	}
#endif
#endif
#endif
#if !JVET_C0024_QTBT
	delete[] m_puhQTTempTrIdx;
	delete[] m_pcQTTempTComYuv;
#endif

#if COM16_C806_EMT
	delete[] m_puhQTTempEmtTuIdx;
	delete[] m_puhQTTempEmtCuFlag;
#endif

#if COM16_C806_LARGE_CTU
	for (Int n = 0; n < NUMBER_OF_STORED_RESIDUAL_TYPES; n++)
	{
#if JVET_C0024_QTBT
		for (Int widx = 0; widx < uiNumLayersAllocated; widx++)
		{
			for (Int hidx = 0; hidx<uiNumLayersAllocated; hidx++)
			{
				if (m_resiSingleBuffer[widx][hidx][n])
				{
					delete[] m_resiSingleBuffer[widx][hidx][n];
					m_resiSingleBuffer[widx][hidx][n] = NULL;
				}
			}
		}
#else
		for (Int depth = 0; depth < MAX_CU_DEPTH; depth++)
		{
			if (m_resiSplitBuffer[depth][n])
			{
				delete[] m_resiSplitBuffer[depth][n];
				m_resiSplitBuffer[depth][n] = NULL;
			}
			if (m_resiSingleBuffer[depth][n])
			{
				delete[] m_resiSingleBuffer[depth][n];
				m_resiSingleBuffer[depth][n] = NULL;
			}
		}
#endif
		if (m_resiPUBuffer[n])
		{
			delete[] m_resiPUBuffer[n];
			m_resiPUBuffer[n] = NULL;
		}
	}
#endif

	for (UInt ch = 0; ch<MAX_NUM_COMPONENT; ch++)
	{
		delete[] m_pSharedPredTransformSkip[ch];
		delete[] m_pcQTTempTUCoeff[ch];
#if ADAPTIVE_QP_SELECTION
		delete[] m_ppcQTTempTUArlCoeff[ch];
#endif
		delete[] m_phQTTempCrossComponentPredictionAlpha[ch];
		delete[] m_puhQTTempTransformSkipFlag[ch];
#if VCEG_AZ08_KLT_COMMON
		delete[] m_puhQTTempKLTFlag[ch];
#endif
#if COM16_C806_EMT
		delete[] m_puhQTTempExplicitRdpcmMode[ch];
#endif
	}
	m_pcQTTempTransformSkipTComYuv.destroy();
#if COM16_C806_LMCHROMA
	m_pcQTTempResiTComYuv.destroy();
#endif

	m_tmpYuvPred.destroy();
	m_isInitialized = false;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
#if JVET_C0035_ATMVP_SIMPLIFICATION
	for (UInt ui = 0; ui<NUM_MGR_TYPE; ui++)
#else
	for (UInt ui = 0; ui<2; ui++)
#endif
	{
		if (m_pMvFieldSP[ui] != NULL)
		{
			delete[] m_pMvFieldSP[ui];
			m_pMvFieldSP[ui] = NULL;
		}
		if (m_phInterDirSP[ui] != NULL)
		{
			delete[] m_phInterDirSP[ui];
			m_phInterDirSP[ui] = NULL;
		}
	}
#endif
#if VCEG_AZ07_FRUC_MERGE
	if (m_pMvFieldFRUC != NULL)
	{
		delete[] m_pMvFieldFRUC;
		m_pMvFieldFRUC = NULL;
	}
	if (m_phInterDirFRUC != NULL)
	{
		delete[] m_phInterDirFRUC;
		m_phInterDirFRUC = NULL;
	}
	if (m_phFRUCRefineDist[0] != NULL)
	{
		delete[] m_phFRUCRefineDist[0];
		m_phFRUCRefineDist[0] = NULL;
	}
	if (m_phFRUCRefineDist[1] != NULL)
	{
		delete[] m_phFRUCRefineDist[1];
		m_phFRUCRefineDist[1] = NULL;
	}
	if (m_phFRUCSBlkRefineDist[0] != NULL)
	{
		delete[] m_phFRUCSBlkRefineDist[0];
		m_phFRUCSBlkRefineDist[0] = NULL;
	}
	if (m_phFRUCSBlkRefineDist[1] != NULL)
	{
		delete[] m_phFRUCSBlkRefineDist[1];
		m_phFRUCSBlkRefineDist[1] = NULL;
	}
#endif

#if COM16_C1016_AFFINE
	if (m_tmpError != NULL)
	{
		delete[] m_tmpError;
	}
	if (m_tmpDerivate[0] != NULL)
	{
		delete[] m_tmpDerivate[0];
	}

	if (m_tmpDerivate[1] != NULL)
	{
		delete[] m_tmpDerivate[1];
	}
#endif



#if JVET_E0077_LM_MF
	for (Int k = 0; k < 8; k++)
	{
		delete[]m_pLMMFPredSaved[k];
	}
	delete[]m_pLMMFPredSaved;
#endif
}

TEncSearch::~TEncSearch()
{
	if (m_isInitialized)
	{
		destroy();
	}
#if JVET_D0123_ME_CTX_LUT_BITS
	delete m_pcPuMeEstBitsSbac;
#endif
}




Void TEncSearch::init(TEncCfg*      pcEncCfg,
	TComTrQuant*  pcTrQuant,
	Int           iSearchRange,
	Int           bipredSearchRange,
	Int           iFastSearch,
	const UInt    maxCUWidth,
	const UInt    maxCUHeight,
	const UInt    maxTotalCUDepth,
	TEncEntropy*  pcEntropyCoder,
	TComRdCost*   pcRdCost,
#if JVET_C0024_QTBT
	TEncSbac****  ppppcRDSbacCoder,
#else
	TEncSbac*** pppcRDSbacCoder,
#endif
	TEncSbac*   pcRDGoOnSbacCoder
)
{
	assert(!m_isInitialized);
	m_pcEncCfg = pcEncCfg;
	m_pcTrQuant = pcTrQuant;
	m_iSearchRange = iSearchRange;
	m_bipredSearchRange = bipredSearchRange;
	m_iFastSearch = iFastSearch;
	m_pcEntropyCoder = pcEntropyCoder;
	m_pcRdCost = pcRdCost;

#if JVET_C0024_QTBT
	m_ppppcRDSbacCoder = ppppcRDSbacCoder;
#else
	m_pppcRDSbacCoder = pppcRDSbacCoder;
#endif
	m_pcRDGoOnSbacCoder = pcRDGoOnSbacCoder;

	for (UInt iDir = 0; iDir < MAX_NUM_REF_LIST_ADAPT_SR; iDir++)
	{
		for (UInt iRefIdx = 0; iRefIdx < MAX_IDX_ADAPT_SR; iRefIdx++)
		{
			m_aaiAdaptSR[iDir][iRefIdx] = iSearchRange;
		}
	}

	m_puiDFilter = s_auiDFilter + 4;

	// initialize motion cost
	for (Int iNum = 0; iNum < AMVP_MAX_NUM_CANDS + 1; iNum++)
	{
		for (Int iIdx = 0; iIdx < AMVP_MAX_NUM_CANDS; iIdx++)
		{
			if (iIdx < iNum)
			{
				m_auiMVPIdxCost[iIdx][iNum] = xGetMvpIdxBits(iIdx, iNum);
			}
			else
			{
				m_auiMVPIdxCost[iIdx][iNum] = MAX_INT;
			}
		}
	}

	const ChromaFormat cform = pcEncCfg->getChromaFormatIdc();
#if COM16_C806_LMCHROMA
	const Int bitDepth = pcEncCfg->getBitDepth(CHANNEL_TYPE_LUMA);
	initTempBuff(cform, bitDepth
#if VCEG_AZ08_INTER_KLT
		, pcEncCfg->getUseInterKLT(), pcEncCfg->getSourceWidth(), pcEncCfg->getSourceHeight(), maxCUWidth, maxCUHeight, maxTotalCUDepth
#endif
	);
#else
	initTempBuff(cform
#if VCEG_AZ08_INTER_KLT
		, pcEncCfg->getUseInterKLT(), pcEncCfg->getSourceWidth(), pcEncCfg->getSourceHeight(), maxCUWidth, maxCUHeight, maxTotalCUDepth
#endif
	);
#endif

	m_pTempPel = new Pel[maxCUWidth*maxCUHeight];

#if JVET_C0024_QTBT
	const UInt uiNumLayersToAllocate = g_aucConvertToBit[pcEncCfg->getCTUSize()] + 1;
#else
	const UInt uiNumLayersToAllocate = pcEncCfg->getQuadtreeTULog2MaxSize() - pcEncCfg->getQuadtreeTULog2MinSize() + 1;
#endif
	const UInt uiNumPartitions = 1 << (maxTotalCUDepth << 1);
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
	if (m_pcTrQuant->getUseInterKLT())
	{
#endif
#if JVET_C0024_QTBT
		m_ppcQTTempTComYuvRec = new TComYuv*[uiNumLayersToAllocate];
#else
		m_pcQTTempTComYuvRec = new TComYuv[uiNumLayersToAllocate];
#endif
#if VCEG_AZ08_USE_KLT
	}
#endif
#endif
	for (UInt ch = 0; ch<MAX_NUM_COMPONENT; ch++)
	{
		const UInt csx = ::getComponentScaleX(ComponentID(ch), cform);
		const UInt csy = ::getComponentScaleY(ComponentID(ch), cform);
#if JVET_C0024_QTBT
		m_pppcQTTempCoeff[ch] = new TCoeff**[uiNumLayersToAllocate];
#else
		m_ppcQTTempCoeff[ch] = new TCoeff*[uiNumLayersToAllocate];
#endif
		m_pcQTTempCoeff[ch] = new TCoeff[(maxCUWidth*maxCUHeight) >> (csx + csy)];
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
		m_pppcQTTempArlCoeff[ch] = new TCoeff**[uiNumLayersToAllocate];
#else
		m_ppcQTTempArlCoeff[ch] = new TCoeff*[uiNumLayersToAllocate];
#endif
		m_pcQTTempArlCoeff[ch] = new TCoeff[(maxCUWidth*maxCUHeight) >> (csx + csy)];
#endif
		m_puhQTTempCbf[ch] = new UChar[uiNumPartitions];

#if JVET_C0024_QTBT
		for (UInt uiWIdx = 0; uiWIdx<uiNumLayersToAllocate; uiWIdx++)
		{
			m_pppcQTTempCoeff[ch][uiWIdx] = new TCoeff*[uiNumLayersToAllocate];
#if ADAPTIVE_QP_SELECTION
			m_pppcQTTempArlCoeff[ch][uiWIdx] = new TCoeff*[uiNumLayersToAllocate];
#endif
			for (UInt uiHIdx = 0; uiHIdx<uiNumLayersToAllocate; uiHIdx++)
			{
				m_pppcQTTempCoeff[ch][uiWIdx][uiHIdx] = new TCoeff[(1 << (uiWIdx + uiHIdx + (MIN_CU_LOG2 << 1))) >> (csx + csy)];
#if ADAPTIVE_QP_SELECTION
				m_pppcQTTempArlCoeff[ch][uiWIdx][uiHIdx] = new TCoeff[(1 << (uiWIdx + uiHIdx + (MIN_CU_LOG2 << 1))) >> (csx + csy)];
#endif
			}
		}
#else
		for (UInt layer = 0; layer < uiNumLayersToAllocate; layer++)
		{
			m_ppcQTTempCoeff[ch][layer] = new TCoeff[(maxCUWidth*maxCUHeight) >> (csx + csy)];
#if ADAPTIVE_QP_SELECTION
			m_ppcQTTempArlCoeff[ch][layer] = new TCoeff[(maxCUWidth*maxCUHeight) >> (csx + csy)];
#endif
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
			if (m_pcTrQuant->getUseInterKLT())
			{
#endif
				if (ch == 0)
				{
					m_pcQTTempTComYuvRec[layer].create(maxCUWidth, maxCUHeight, cform);
				}
#if VCEG_AZ08_USE_KLT
			}
#endif
#endif
		}
#endif

		m_phQTTempCrossComponentPredictionAlpha[ch] = new Char[uiNumPartitions];
		m_pSharedPredTransformSkip[ch] = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];
		m_pcQTTempTUCoeff[ch] = new TCoeff[MAX_CU_SIZE*MAX_CU_SIZE];
#if ADAPTIVE_QP_SELECTION
		m_ppcQTTempTUArlCoeff[ch] = new TCoeff[MAX_CU_SIZE*MAX_CU_SIZE];
#endif
		m_puhQTTempTransformSkipFlag[ch] = new UChar[uiNumPartitions];
#if VCEG_AZ08_KLT_COMMON
		m_puhQTTempKLTFlag[ch] = new UChar[uiNumPartitions];
#endif
#if COM16_C806_EMT
		m_puhQTTempExplicitRdpcmMode[ch] = new UChar[uiNumPartitions];
#endif
	}
#if JVET_C0024_QTBT
	m_ppcQTTempTComYuv = new TComYuv*[uiNumLayersToAllocate];
#else
	m_puhQTTempTrIdx = new UChar[uiNumPartitions];
	m_pcQTTempTComYuv = new TComYuv[uiNumLayersToAllocate];
#endif
#if COM16_C806_EMT
	m_puhQTTempEmtTuIdx = new UChar[uiNumPartitions];
	m_puhQTTempEmtCuFlag = new UChar[uiNumPartitions];
#endif
#if JVET_C0024_QTBT
	for (UInt uiWIdx = 0; uiWIdx<uiNumLayersToAllocate; uiWIdx++)
	{
		m_ppcQTTempTComYuv[uiWIdx] = new TComYuv[uiNumLayersToAllocate];
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
		if (m_pcTrQuant->getUseInterKLT())
		{
#endif
			m_ppcQTTempTComYuvRec[uiWIdx] = new TComYuv[uiNumLayersToAllocate];
#if VCEG_AZ08_USE_KLT
		}
#endif
#endif
		for (UInt uiHIdx = 0; uiHIdx<uiNumLayersToAllocate; uiHIdx++)
		{
			m_ppcQTTempTComYuv[uiWIdx][uiHIdx].create(1 << (uiWIdx + MIN_CU_LOG2), 1 << (uiHIdx + MIN_CU_LOG2), pcEncCfg->getChromaFormatIdc());
#if JVET_C0024_QTBT
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
			if (m_pcTrQuant->getUseInterKLT())
			{
#endif
				m_ppcQTTempTComYuvRec[uiWIdx][uiHIdx].create(1 << (uiWIdx + MIN_CU_LOG2), 1 << (uiHIdx + MIN_CU_LOG2), cform);
#if VCEG_AZ08_USE_KLT
			}
#endif
#endif
#endif
		}
	}
#else
	for (UInt ui = 0; ui < uiNumLayersToAllocate; ++ui)
	{
		m_pcQTTempTComYuv[ui].create(maxCUWidth, maxCUHeight, pcEncCfg->getChromaFormatIdc());
	}
#endif
	m_pcQTTempTransformSkipTComYuv.create(maxCUWidth, maxCUHeight, pcEncCfg->getChromaFormatIdc());
	m_tmpYuvPred.create(MAX_CU_SIZE, MAX_CU_SIZE, pcEncCfg->getChromaFormatIdc());

#if COM16_C806_LARGE_CTU
	for (Int n = 0; n < NUMBER_OF_STORED_RESIDUAL_TYPES; n++)
	{
#if JVET_C0024_QTBT
		for (Int widx = 0; widx<uiNumLayersToAllocate; widx++)
		{
			for (Int hidx = 0; hidx<uiNumLayersToAllocate; hidx++)
			{
				m_resiSingleBuffer[widx][hidx][n] = new Pel[(size_t)1 << (widx + hidx + (MIN_CU_LOG2 << 1))];
			}
		}
#else
		for (Int depth = 0; depth < MAX_CU_DEPTH; depth++)
		{
			m_resiSplitBuffer[depth][n] = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];
			m_resiSingleBuffer[depth][n] = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];
		}
#endif
		m_resiPUBuffer[n] = new Pel[MAX_CU_SIZE*MAX_CU_SIZE];
	}
#endif
#if COM16_C806_LMCHROMA
	m_pcQTTempResiTComYuv.create(MAX_CU_SIZE, MAX_CU_SIZE, pcEncCfg->getChromaFormatIdc());
#endif

#if COM16_C1016_AFFINE
#if !ITH_C2_FXD_PT_AFF_PARAMS
	m_tmpError = new Int[MAX_CU_SIZE * MAX_CU_SIZE];
	m_tmpDerivate[0] = new Double[MAX_CU_SIZE * MAX_CU_SIZE];
	m_tmpDerivate[1] = new Double[MAX_CU_SIZE * MAX_CU_SIZE];
#else
	m_tmpDerivate[0] = new Int[MAX_CU_SIZE * MAX_CU_SIZE];
	m_tmpDerivate[1] = new Int[MAX_CU_SIZE * MAX_CU_SIZE];
	m_tmpError = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#endif
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO
	for (UInt uiWIdx = 0; uiWIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiWIdx++)
	{
		for (UInt uiHIdx = 0; uiHIdx <= MAX_CU_DEPTH - MIN_CU_LOG2; uiHIdx++)
		{
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
			m_SaveLoadNumRouteSaved[uiWIdx][uiHIdx] = (UInt*)malloc(MAX_PART_IN_CTU * sizeof(UInt));
			m_SaveLoadRouteCodingOrder[uiWIdx][uiHIdx] = (UInt*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UInt));
			m_SaveLoadRouteLeftNeigbor[uiWIdx][uiHIdx] = (UInt*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UInt));
			m_SaveLoadRouteUpNeigbor[uiWIdx][uiHIdx] = (UInt*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UInt));

			m_SaveLoadTag[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
#if COM16_C806_EMT
			m_SaveLoadEmtFlag[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
			m_SaveLoadEmtIdx[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST 
			m_SaveLoadRotIdx[uiWIdx][uiHIdx] = (Char *)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(Char));
#endif
#if COM16_C1046_PDPC_INTRA
			m_SaveLoadPdpcIdx[uiWIdx][uiHIdx] = (Char *)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(Char));
#endif
#if VCEG_AZ07_FRUC_MERGE
			m_SaveLoadFrucMode[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
#endif
#if JVET_E0076_MULTI_PEL_MVD
			m_SaveLoadIMVFlag[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
#endif
#if VCEG_AZ06_IC
			m_SaveLoadICFlag[uiWIdx][uiHIdx] = (Bool *)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(Bool));
#endif
#if COM16_C1016_AFFINE
			m_SaveLoadAffineFlag[uiWIdx][uiHIdx] = (Bool *)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(Bool));
#endif
			m_SaveLoadMergeFlag[uiWIdx][uiHIdx] = (Bool *)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(Bool));
			m_SaveLoadInterDir[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));
			m_SaveLoadSplit[uiWIdx][uiHIdx] = (UChar*)malloc(MAX_PART_IN_CTU*MAX_NUM_ROUTE * sizeof(UChar));

			for (UInt uiPartIdx = 0; uiPartIdx<MAX_PART_IN_CTU*MAX_NUM_ROUTE; uiPartIdx++) {
				//m_SaveLoadPartIdx[uiWIdx][uiHIdx] = MAX_UINT;
				m_SaveLoadTag[uiWIdx][uiHIdx][uiPartIdx] = SAVE_LOAD_INIT;
#if COM16_C806_EMT
				m_SaveLoadEmtFlag[uiWIdx][uiHIdx][uiPartIdx] = 0;
				m_SaveLoadEmtIdx[uiWIdx][uiHIdx][uiPartIdx] = 0;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
				m_SaveLoadRotIdx[uiWIdx][uiHIdx][uiPartIdx] = 0;
#endif
#if COM16_C1046_PDPC_INTRA
				m_SaveLoadPdpcIdx[uiWIdx][uiHIdx][uiPartIdx] = 0;
#endif
#if VCEG_AZ07_FRUC_MERGE
				m_SaveLoadFrucMode[uiWIdx][uiHIdx][uiPartIdx] = 0;
#endif
#if VCEG_AZ07_IMV
				m_SaveLoadIMVFlag[uiWIdx][uiHIdx][uiPartIdx] = false;
#endif
#if VCEG_AZ06_IC
				m_SaveLoadICFlag[uiWIdx][uiHIdx][uiPartIdx] = false;
#endif
#if COM16_C1016_AFFINE
				m_SaveLoadAffineFlag[uiWIdx][uiHIdx][uiPartIdx] = false;
#endif
				m_SaveLoadMergeFlag[uiWIdx][uiHIdx][uiPartIdx] = false;
				m_SaveLoadInterDir[uiWIdx][uiHIdx][uiPartIdx] = 0;
				m_SaveLoadSplit[uiWIdx][uiHIdx][uiPartIdx] = 0;
			}
#else
			m_SaveLoadPartIdx[uiWIdx][uiHIdx] = MAX_UINT;
			m_SaveLoadTag[uiWIdx][uiHIdx] = SAVE_LOAD_INIT;
#if COM16_C806_EMT
			m_SaveLoadEmtFlag[uiWIdx][uiHIdx] = 0;
			m_SaveLoadEmtIdx[uiWIdx][uiHIdx] = 0;
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST
			m_SaveLoadRotIdx[uiWIdx][uiHIdx] = 0;
#endif
#if RRC_INTRA_SHARPENING && INTRA_SHARP_D0077
			m_SaveLoadSharpIdx[uiWIdx][uiHIdx] = 0;
#endif
#if COM16_C1046_PDPC_INTRA
			m_SaveLoadPdpcIdx[uiWIdx][uiHIdx] = 0;
#endif
#if VCEG_AZ07_FRUC_MERGE
			m_SaveLoadFrucMode[uiWIdx][uiHIdx] = 0;
#endif
#if VCEG_AZ07_IMV
			m_SaveLoadIMVFlag[uiWIdx][uiHIdx] = false;
#endif
#if VCEG_AZ06_IC
			m_SaveLoadICFlag[uiWIdx][uiHIdx] = false;
#endif
#if COM16_C1016_AFFINE
			m_SaveLoadAffineFlag[uiWIdx][uiHIdx] = false;
#endif
			m_SaveLoadMergeFlag[uiWIdx][uiHIdx] = false;
			m_SaveLoadInterDir[uiWIdx][uiHIdx] = 0;
			m_SaveLoadSplit[uiWIdx][uiHIdx] = 0;
#endif
		}
	}
#endif

	m_isInitialized = true;
}

#define TZ_SEARCH_CONFIGURATION                                                                                 \
const Int  iRaster                  = 5;  /* TZ soll von aussen ?ergeben werden */                            \
const Bool bTestOtherPredictedMV    = 0;                                                                      \
const Bool bTestZeroVector          = 1;                                                                      \
const Bool bTestZeroVectorStart     = 0;                                                                      \
const Bool bTestZeroVectorStop      = 0;                                                                      \
const Bool bFirstSearchDiamond      = 1;  /* 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch */        \
const Bool bFirstSearchStop         = m_pcEncCfg->getFastMEAssumingSmootherMVEnabled();                       \
const UInt uiFirstSearchRounds      = 3;  /* first search stop X rounds after best match (must be >=1) */     \
const Bool bEnableRasterSearch      = 1;                                                                      \
const Bool bAlwaysRasterSearch      = 0;  /* ===== 1: BETTER but factor 2 slower ===== */                     \
const Bool bRasterRefinementEnable  = 0;  /* enable either raster refinement or star refinement */            \
const Bool bRasterRefinementDiamond = 0;  /* 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch */        \
const Bool bStarRefinementEnable    = 1;  /* enable either star refinement or raster refinement */            \
const Bool bStarRefinementDiamond   = 1;  /* 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch */        \
const Bool bStarRefinementStop      = 0;                                                                      \
const UInt uiStarRefinementRounds   = 2;  /* star refinement stop X rounds after best match (must be >=1) */  \


#define SEL_SEARCH_CONFIGURATION                                                                                 \
  const Bool bTestOtherPredictedMV    = 1;                                                                       \
  const Bool bTestZeroVector          = 1;                                                                       \
  const Bool bEnableRasterSearch      = 1;                                                                       \
  const Bool bAlwaysRasterSearch      = 0;  /* ===== 1: BETTER but factor 15x slower ===== */                    \
  const Bool bStarRefinementEnable    = 1;  /* enable either star refinement or raster refinement */             \
  const Bool bStarRefinementDiamond   = 1;  /* 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch */         \
  const Bool bStarRefinementStop      = 0;                                                                       \
  const UInt uiStarRefinementRounds   = 2;  /* star refinement stop X rounds after best match (must be >=1) */   \
  const UInt uiSearchRange            = m_iSearchRange;                                                          \
  const Int  uiSearchRangeInitial     = m_iSearchRange >> 2;                                                     \
  const Int  uiSearchStep             = 4;                                                                       \
  const Int  iMVDistThresh            = 8;                                                                       \

#if JVET_D0123_ME_CTX_LUT_BITS
__inline UInt TEncSearch::xGetBitsEpExGolomb(UInt uiSymbol, UInt uiCount)
{
	UInt numBins = 0;
	while (uiSymbol >= (UInt)(1 << uiCount))
	{
		numBins++;
		uiSymbol -= 1 << uiCount;
		uiCount++;
	}
	numBins++;
	numBins += uiCount;

	assert(numBins <= 32);
	return (numBins << 15);
}

__inline UInt TEncSearch::xGetMvdBitsLut(TComDataCU* pcCU, UInt uiAbsPartIdx, Int iCurVal)
{
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
	assert(iCurVal == (iCurVal >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));
	iCurVal >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif

#if VCEG_AZ07_IMV
	if (pcCU->getiMVFlag(uiAbsPartIdx) && pcCU->getSlice()->getSPS()->getIMV())
	{
		iCurVal >>= 2;

#if JVET_E0076_MULTI_PEL_MVD
		if (pcCU->getiMVFlag(uiAbsPartIdx) == 2)
		{
			iCurVal >>= MULTI_PEL_MVD_BITS;
		}
#endif
	}
#endif

	UInt uiCtxIdx = 0;
	UInt uiBits = 0;
	uiBits += m_pcPuMeEstBitsSbac->mvdBits[uiCtxIdx][(iCurVal != 0) ? 1 : 0];
	const Bool bValAbsGr0 = iCurVal != 0;
	const UInt uiValAbs = 0 > iCurVal ? -iCurVal : iCurVal;

	uiCtxIdx++;
	if (bValAbsGr0)
	{
		uiBits += m_pcPuMeEstBitsSbac->mvdBits[uiCtxIdx][uiValAbs > 1 ? 1 : 0];
	}

	if (bValAbsGr0)
	{
		if (uiValAbs > 1)
		{
			uiBits += xGetBitsEpExGolomb(uiValAbs - 2, 1);
		}
		uiBits += (UInt)EPBIT;
	}

	return uiBits;
}

__inline UInt TEncSearch::xMvdBits(TComPattern* pcPatternKey, Int iX, Int iY)
{
	TComDataCU* pcCU = pcPatternKey->m_pcCU;
	UInt uiPartAddr = pcPatternKey->m_uiPartAddr;
	Int iMvpX = pcPatternKey->m_mvPred->getHor();
	Int iMvpY = pcPatternKey->m_mvPred->getVer();
	Int iCurX = (iX << iCostScale) - iMvpX;
	Int iCurY = (iY << iCostScale) - iMvpY;

	UInt uiBits = xGetMvdBitsLut(pcCU, uiPartAddr, iCurX) + xGetMvdBitsLut(pcCU, uiPartAddr, iCurY);

	return uiBits;
}
#endif


__inline Void TEncSearch::xTZSearchHelp(TComPattern* pcPatternKey, IntTZSearchStruct& rcStruct, const Int iSearchX, const Int iSearchY, const UChar ucPointNr, const UInt uiDistance)
{
	Distortion  uiSad = 0;

	Pel*  piRefSrch;

	piRefSrch = rcStruct.piRefY + iSearchY * rcStruct.iYStride + iSearchX;
#if VCEG_AZ06_IC
	m_cDistParam.bMRFlag = pcPatternKey->getMRFlag();
#endif  
	//-- jclee for using the SAD function pointer
	m_pcRdCost->setDistParam(pcPatternKey, piRefSrch, rcStruct.iYStride, m_cDistParam);

	if (m_pcEncCfg->getFastSearch() != SELECTIVE)
	{
		// fast encoder decision: use subsampled SAD when rows > 8 for integer ME
		if (m_pcEncCfg->getUseFastEnc())
		{
			if (m_cDistParam.iRows > 8)
			{
				m_cDistParam.iSubShift = 1;
			}
		}
	}

	setDistParamComp(COMPONENT_Y);

	// distortion
	m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
	if (m_pcEncCfg->getFastSearch() == SELECTIVE)
	{
		Int isubShift = 0;
		// motion cost
#if JVET_D0123_ME_CTX_LUT_BITS
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		Distortion uiBitCost = xMvdBits(pcPatternKey, iSearchX << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, iSearchY << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
		Distortion uiBitCost = xMvdBits(pcPatternKey, iSearchX, iSearchY);
#endif
		uiBitCost = m_pcRdCost->getCostSearch(uiBitCost);
#else
		Distortion uiBitCost = m_pcRdCost->getCost(iSearchX, iSearchY
#if VCEG_AZ07_IMV
			, pcPatternKey->getImvFlag()
#endif
		);
#endif
		if (m_cDistParam.iRows > 32)
		{
			m_cDistParam.iSubShift = 4;
		}
		else if (m_cDistParam.iRows > 16)
		{
			m_cDistParam.iSubShift = 3;
		}
		else if (m_cDistParam.iRows > 8)
		{
			m_cDistParam.iSubShift = 2;
		}
		else
		{
			m_cDistParam.iSubShift = 1;
		}

		Distortion uiTempSad = m_cDistParam.DistFunc(&m_cDistParam);
		if ((uiTempSad + uiBitCost) < rcStruct.uiBestSad)
		{
			uiSad += uiTempSad >> m_cDistParam.iSubShift;
			while (m_cDistParam.iSubShift > 0)
			{
				isubShift = m_cDistParam.iSubShift - 1;
				m_cDistParam.pOrg = pcPatternKey->getROIY() + (pcPatternKey->getPatternLStride() << isubShift);
				m_cDistParam.pCur = piRefSrch + (rcStruct.iYStride << isubShift);
				uiTempSad = m_cDistParam.DistFunc(&m_cDistParam);
				uiSad += uiTempSad >> m_cDistParam.iSubShift;
				if (((uiSad << isubShift) + uiBitCost) > rcStruct.uiBestSad)
				{
					break;
				}

				m_cDistParam.iSubShift--;
			}

			if (m_cDistParam.iSubShift == 0)
			{
#if JVET_D0123_ME_CTX_LUT_BITS
				UInt uiCurSad = uiSad;
#endif
				uiSad += uiBitCost;
				if (uiSad < rcStruct.uiBestSad)
				{
#if JVET_D0123_ME_CTX_LUT_BITS
					rcStruct.uiBestSadOnly = uiCurSad;
#endif
					rcStruct.uiBestSad = uiSad;
					rcStruct.iBestX = iSearchX;
					rcStruct.iBestY = iSearchY;
					rcStruct.uiBestDistance = uiDistance;
					rcStruct.uiBestRound = 0;
					rcStruct.ucPointNr = ucPointNr;
				}
			}
		}
	}
	else
	{
		uiSad = m_cDistParam.DistFunc(&m_cDistParam);

		// motion cost
#if JVET_D0123_ME_CTX_LUT_BITS
		Distortion uiCurSad = uiSad;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		Distortion uiMvdBits = xMvdBits(pcPatternKey, iSearchX << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, iSearchY << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
		Distortion uiMvdBits = xMvdBits(pcPatternKey, iSearchX, iSearchY);
#endif
		uiSad += (UInt)m_pcRdCost->getCostSearch(uiMvdBits);
#else
		uiSad += m_pcRdCost->getCost(iSearchX, iSearchY
#if VCEG_AZ07_IMV
			, pcPatternKey->getImvFlag()
#endif
		);
#endif

		if (uiSad < rcStruct.uiBestSad)
		{
#if JVET_D0123_ME_CTX_LUT_BITS
			rcStruct.uiBestSadOnly = uiCurSad;
#endif
			rcStruct.uiBestSad = uiSad;
			rcStruct.iBestX = iSearchX;
			rcStruct.iBestY = iSearchY;
			rcStruct.uiBestDistance = uiDistance;
			rcStruct.uiBestRound = 0;
			rcStruct.ucPointNr = ucPointNr;
		}
	}
}




__inline Void TEncSearch::xTZ2PointSearch(TComPattern* pcPatternKey, IntTZSearchStruct& rcStruct, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB)
{
	Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
	Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
	Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
	Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

	// 2 point search,                   //   1 2 3
	// check only the 2 untested points  //   4 0 5
	// around the start point            //   6 7 8
	Int iStartX = rcStruct.iBestX;
	Int iStartY = rcStruct.iBestY;
	switch (rcStruct.ucPointNr)
	{
	case 1:
	{
		if ((iStartX - 1) >= iSrchRngHorLeft)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY, 0, 2);
		}
		if ((iStartY - 1) >= iSrchRngVerTop)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iStartY - 1, 0, 2);
		}
	}
	break;
	case 2:
	{
		if ((iStartY - 1) >= iSrchRngVerTop)
		{
			if ((iStartX - 1) >= iSrchRngHorLeft)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY - 1, 0, 2);
			}
			if ((iStartX + 1) <= iSrchRngHorRight)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY - 1, 0, 2);
			}
		}
	}
	break;
	case 3:
	{
		if ((iStartY - 1) >= iSrchRngVerTop)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iStartY - 1, 0, 2);
		}
		if ((iStartX + 1) <= iSrchRngHorRight)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY, 0, 2);
		}
	}
	break;
	case 4:
	{
		if ((iStartX - 1) >= iSrchRngHorLeft)
		{
			if ((iStartY + 1) <= iSrchRngVerBottom)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY + 1, 0, 2);
			}
			if ((iStartY - 1) >= iSrchRngVerTop)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY - 1, 0, 2);
			}
		}
	}
	break;
	case 5:
	{
		if ((iStartX + 1) <= iSrchRngHorRight)
		{
			if ((iStartY - 1) >= iSrchRngVerTop)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY - 1, 0, 2);
			}
			if ((iStartY + 1) <= iSrchRngVerBottom)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY + 1, 0, 2);
			}
		}
	}
	break;
	case 6:
	{
		if ((iStartX - 1) >= iSrchRngHorLeft)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY, 0, 2);
		}
		if ((iStartY + 1) <= iSrchRngVerBottom)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iStartY + 1, 0, 2);
		}
	}
	break;
	case 7:
	{
		if ((iStartY + 1) <= iSrchRngVerBottom)
		{
			if ((iStartX - 1) >= iSrchRngHorLeft)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iStartY + 1, 0, 2);
			}
			if ((iStartX + 1) <= iSrchRngHorRight)
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY + 1, 0, 2);
			}
		}
	}
	break;
	case 8:
	{
		if ((iStartX + 1) <= iSrchRngHorRight)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX + 1, iStartY, 0, 2);
		}
		if ((iStartY + 1) <= iSrchRngVerBottom)
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iStartY + 1, 0, 2);
		}
	}
	break;
	default:
	{
		assert(false);
	}
	break;
	} // switch( rcStruct.ucPointNr )
}




__inline Void TEncSearch::xTZ8PointSquareSearch(TComPattern* pcPatternKey, IntTZSearchStruct& rcStruct, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist)
{
	Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
	Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
	Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
	Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

	// 8 point search,                   //   1 2 3
	// search around the start point     //   4 0 5
	// with the required  distance       //   6 7 8
	assert(iDist != 0);
	const Int iTop = iStartY - iDist;
	const Int iBottom = iStartY + iDist;
	const Int iLeft = iStartX - iDist;
	const Int iRight = iStartX + iDist;
	rcStruct.uiBestRound += 1;

	if (iTop >= iSrchRngVerTop) // check top
	{
		if (iLeft >= iSrchRngHorLeft) // check top left
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iTop, 1, iDist);
		}
		// top middle
		xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 2, iDist);

		if (iRight <= iSrchRngHorRight) // check top right
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iRight, iTop, 3, iDist);
		}
	} // check top
	if (iLeft >= iSrchRngHorLeft) // check middle left
	{
		xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist);
	}
	if (iRight <= iSrchRngHorRight) // check middle right
	{
		xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 5, iDist);
	}
	if (iBottom <= iSrchRngVerBottom) // check bottom
	{
		if (iLeft >= iSrchRngHorLeft) // check bottom left
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iBottom, 6, iDist);
		}
		// check bottom middle
		xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist);

		if (iRight <= iSrchRngHorRight) // check bottom right
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iRight, iBottom, 8, iDist);
		}
	} // check bottom
}




__inline Void TEncSearch::xTZ8PointDiamondSearch(TComPattern* pcPatternKey, IntTZSearchStruct& rcStruct, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist)
{
	Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
	Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
	Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
	Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

	// 8 point search,                   //   1 2 3
	// search around the start point     //   4 0 5
	// with the required  distance       //   6 7 8
	assert(iDist != 0);
	const Int iTop = iStartY - iDist;
	const Int iBottom = iStartY + iDist;
	const Int iLeft = iStartX - iDist;
	const Int iRight = iStartX + iDist;
	rcStruct.uiBestRound += 1;

	if (iDist == 1) // iDist == 1
	{
		if (iTop >= iSrchRngVerTop) // check top
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 2, iDist);
		}
		if (iLeft >= iSrchRngHorLeft) // check middle left
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist);
		}
		if (iRight <= iSrchRngHorRight) // check middle right
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 5, iDist);
		}
		if (iBottom <= iSrchRngVerBottom) // check bottom
		{
			xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist);
		}
	}
	else // if (iDist != 1)
	{
		if (iDist <= 8)
		{
			const Int iTop_2 = iStartY - (iDist >> 1);
			const Int iBottom_2 = iStartY + (iDist >> 1);
			const Int iLeft_2 = iStartX - (iDist >> 1);
			const Int iRight_2 = iStartX + (iDist >> 1);

			if (iTop >= iSrchRngVerTop && iLeft >= iSrchRngHorLeft &&
				iRight <= iSrchRngHorRight && iBottom <= iSrchRngVerBottom) // check border
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 2, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iLeft_2, iTop_2, 1, iDist >> 1);
				xTZSearchHelp(pcPatternKey, rcStruct, iRight_2, iTop_2, 3, iDist >> 1);
				xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 5, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iLeft_2, iBottom_2, 6, iDist >> 1);
				xTZSearchHelp(pcPatternKey, rcStruct, iRight_2, iBottom_2, 8, iDist >> 1);
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist);
			}
			else // check border
			{
				if (iTop >= iSrchRngVerTop) // check top
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 2, iDist);
				}
				if (iTop_2 >= iSrchRngVerTop) // check half top
				{
					if (iLeft_2 >= iSrchRngHorLeft) // check half left
					{
						xTZSearchHelp(pcPatternKey, rcStruct, iLeft_2, iTop_2, 1, (iDist >> 1));
					}
					if (iRight_2 <= iSrchRngHorRight) // check half right
					{
						xTZSearchHelp(pcPatternKey, rcStruct, iRight_2, iTop_2, 3, (iDist >> 1));
					}
				} // check half top
				if (iLeft >= iSrchRngHorLeft) // check left
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist);
				}
				if (iRight <= iSrchRngHorRight) // check right
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 5, iDist);
				}
				if (iBottom_2 <= iSrchRngVerBottom) // check half bottom
				{
					if (iLeft_2 >= iSrchRngHorLeft) // check half left
					{
						xTZSearchHelp(pcPatternKey, rcStruct, iLeft_2, iBottom_2, 6, (iDist >> 1));
					}
					if (iRight_2 <= iSrchRngHorRight) // check half right
					{
						xTZSearchHelp(pcPatternKey, rcStruct, iRight_2, iBottom_2, 8, (iDist >> 1));
					}
				} // check half bottom
				if (iBottom <= iSrchRngVerBottom) // check bottom
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist);
				}
			} // check border
		}
		else // iDist > 8
		{
			if (iTop >= iSrchRngVerTop && iLeft >= iSrchRngHorLeft &&
				iRight <= iSrchRngHorRight && iBottom <= iSrchRngVerBottom) // check border
			{
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 0, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 0, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 0, iDist);
				xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 0, iDist);
				for (Int index = 1; index < 4; index++)
				{
					Int iPosYT = iTop + ((iDist >> 2) * index);
					Int iPosYB = iBottom - ((iDist >> 2) * index);
					Int iPosXL = iStartX - ((iDist >> 2) * index);
					Int iPosXR = iStartX + ((iDist >> 2) * index);
					xTZSearchHelp(pcPatternKey, rcStruct, iPosXL, iPosYT, 0, iDist);
					xTZSearchHelp(pcPatternKey, rcStruct, iPosXR, iPosYT, 0, iDist);
					xTZSearchHelp(pcPatternKey, rcStruct, iPosXL, iPosYB, 0, iDist);
					xTZSearchHelp(pcPatternKey, rcStruct, iPosXR, iPosYB, 0, iDist);
				}
			}
			else // check border
			{
				if (iTop >= iSrchRngVerTop) // check top
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iTop, 0, iDist);
				}
				if (iLeft >= iSrchRngHorLeft) // check left
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iStartY, 0, iDist);
				}
				if (iRight <= iSrchRngHorRight) // check right
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iRight, iStartY, 0, iDist);
				}
				if (iBottom <= iSrchRngVerBottom) // check bottom
				{
					xTZSearchHelp(pcPatternKey, rcStruct, iStartX, iBottom, 0, iDist);
				}
				for (Int index = 1; index < 4; index++)
				{
					Int iPosYT = iTop + ((iDist >> 2) * index);
					Int iPosYB = iBottom - ((iDist >> 2) * index);
					Int iPosXL = iStartX - ((iDist >> 2) * index);
					Int iPosXR = iStartX + ((iDist >> 2) * index);

					if (iPosYT >= iSrchRngVerTop) // check top
					{
						if (iPosXL >= iSrchRngHorLeft) // check left
						{
							xTZSearchHelp(pcPatternKey, rcStruct, iPosXL, iPosYT, 0, iDist);
						}
						if (iPosXR <= iSrchRngHorRight) // check right
						{
							xTZSearchHelp(pcPatternKey, rcStruct, iPosXR, iPosYT, 0, iDist);
						}
					} // check top
					if (iPosYB <= iSrchRngVerBottom) // check bottom
					{
						if (iPosXL >= iSrchRngHorLeft) // check left
						{
							xTZSearchHelp(pcPatternKey, rcStruct, iPosXL, iPosYB, 0, iDist);
						}
						if (iPosXR <= iSrchRngHorRight) // check right
						{
							xTZSearchHelp(pcPatternKey, rcStruct, iPosXR, iPosYB, 0, iDist);
						}
					} // check bottom
				} // for ...
			} // check border
		} // iDist <= 8
	} // iDist == 1
}





//<--

Distortion TEncSearch::xPatternRefinement(TComPattern* pcPatternKey,
	TComMv baseRefMv,
	Int iFrac, TComMv& rcMvFrac,
	Bool bAllowUseOfHadamard
#if RRC_SELECTABLE_IF
	, TComDataCU* pcCU
	, TComPic *pcRefPic
	, UInt uiPartAddr
#endif
)
{
	Distortion  uiDist;
	Distortion  uiDistBest = std::numeric_limits<Distortion>::max();
	UInt        uiDirecBest = 0;

	Pel*  piRefPos;
	Int iRefStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);

	m_pcRdCost->setDistParam(pcPatternKey, m_filteredBlock[0][0].getAddr(COMPONENT_Y), iRefStride, 1, m_cDistParam, m_pcEncCfg->getUseHADME() && bAllowUseOfHadamard);

	const TComMv* pcMvRefine = (iFrac == 2 ? s_acMvRefineH : s_acMvRefineQ);
#if RRC_SELECTABLE_IF
	Int sharpenedBlockMark[4][4];
	memset(&sharpenedBlockMark, 0, sizeof(sharpenedBlockMark));
#endif

	for (UInt i = 0; i < 9; i++)
	{
#if VCEG_AZ07_IMV
		if (pcPatternKey->getImvFlag() && (i > 0 || iFrac == 1))
		{
			continue;
		}
#endif
		TComMv cMvTest = pcMvRefine[i];
		cMvTest += baseRefMv;

		Int horVal = cMvTest.getHor() * iFrac;
		Int verVal = cMvTest.getVer() * iFrac;
		piRefPos = m_filteredBlock[verVal & 3][horVal & 3].getAddr(COMPONENT_Y);
		if (horVal == 2 && (verVal & 1) == 0)
		{
			piRefPos += 1;
		}
		if ((horVal & 1) == 0 && verVal == 2)
		{
			piRefPos += iRefStride;
		}
		cMvTest = pcMvRefine[i];
		cMvTest += rcMvFrac;
#if RRC_SELECTABLE_IF
		{
			TComMv MvSharp = cMvTest;
			if (iFrac == 2)
				MvSharp <<= 1;
#if SIF_SHIFT_MV_SEARCH
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			MvSharp <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
#endif
			Int sharpMode = pcCU->getSIFMode(&MvSharp, uiPartAddr, 0);
			if (sharpMode == SIF_SHARP && !sharpenedBlockMark[verVal & 3][horVal & 3])
			{
				TComYuv *pcYuvPred = &m_filteredBlock[verVal & 3][horVal & 3];
				sharpenedBlockMark[verVal & 3][horVal & 3] = 1;
				pcRefPic->getSharpeningFilter()->sharpBlockLumaForSearch(pcCU, uiPartAddr, 0, m_cDistParam.iCols, m_cDistParam.iRows,
					sharpMode, pcYuvPred, this, &MvSharp, false, pcRefPic->getPicYuvRec());
			}
		}
#endif

		setDistParamComp(COMPONENT_Y);
#if VCEG_AZ06_IC
		m_cDistParam.bMRFlag = pcPatternKey->getMRFlag();
#endif
		m_cDistParam.pCur = piRefPos;
		m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
		uiDist = m_cDistParam.DistFunc(&m_cDistParam);
#if JVET_D0123_ME_CTX_LUT_BITS
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		Distortion uiMvdBits = xMvdBits(pcPatternKey,
			cMvTest.getHor() << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
			cMvTest.getVer() << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
		Distortion uiMvdBits = xMvdBits(pcPatternKey, cMvTest.getHor(), cMvTest.getVer());
#endif
		uiDist += (UInt)m_pcRdCost->getCostSearch(uiMvdBits);
#else
		uiDist += m_pcRdCost->getCost(cMvTest.getHor(), cMvTest.getVer()
#if VCEG_AZ07_IMV
			, pcPatternKey->getImvFlag()
#endif
		);
#endif
		if (uiDist < uiDistBest)
		{
			uiDistBest = uiDist;
			uiDirecBest = i;
		}
	}

	rcMvFrac = pcMvRefine[uiDirecBest];

	return uiDistBest;
}



Void
TEncSearch::xEncSubdivCbfQT(TComTU      &rTu,
	Bool         bLuma,
	Bool         bChroma)
{
	TComDataCU* pcCU = rTu.getCU();
	const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
	const UInt uiTrDepth = rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
	assert(uiTrDepth == 0);
#else
	const UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
	const UInt uiSubdiv = (uiTrMode > uiTrDepth ? 1 : 0);
	const UInt uiLog2LumaTrafoSize = rTu.GetLog2LumaTrSize();

	if (pcCU->isIntra(0) && pcCU->getPartitionSize(0) == SIZE_NxN && uiTrDepth == 0)
	{
		assert(uiSubdiv);
	}
	else if (uiLog2LumaTrafoSize > pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize())
	{
		assert(uiSubdiv);
	}
	else if (uiLog2LumaTrafoSize == pcCU->getSlice()->getSPS()->getQuadtreeTULog2MinSize())
	{
		assert(!uiSubdiv);
	}
	else if (uiLog2LumaTrafoSize == pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx))
	{
		assert(!uiSubdiv);
	}
	else
	{
		assert(uiLog2LumaTrafoSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx));
		if (bLuma)
		{
#if COM16_C806_T64
			m_pcEntropyCoder->encodeTransformSubdivFlag(uiSubdiv, pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2LumaTrafoSize);
#else
			m_pcEntropyCoder->encodeTransformSubdivFlag(uiSubdiv, 5 - uiLog2LumaTrafoSize);
#endif
		}
	}
#endif

	if (bChroma)
	{
		const UInt numberValidComponents = getNumberValidComponents(rTu.GetChromaFormat());
		for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
		{
			const ComponentID compID = ComponentID(ch);
			if (rTu.ProcessingAllQuadrants(compID) && (uiTrDepth == 0 || pcCU->getCbf(uiAbsPartIdx, compID, uiTrDepth - 1)))
			{
#if JVET_C0024_QTBT
				m_pcEntropyCoder->encodeQtCbf(rTu, compID, true);
#else
				m_pcEntropyCoder->encodeQtCbf(rTu, compID, (uiSubdiv == 0));
#endif
			}
		}
	}

#if !JVET_C0024_QTBT
	if (uiSubdiv)
	{
		TComTURecurse tuRecurse(rTu, false);
#if COM16_C806_EMT
		if (0 == uiTrDepth)
		{
			m_pcEntropyCoder->encodeEmtCuFlag(pcCU, uiAbsPartIdx, pcCU->getDepth(0), true);
		}
#endif
		do
		{
			xEncSubdivCbfQT(tuRecurse, bLuma, bChroma);
		} while (tuRecurse.nextSection(rTu));
	}
	else
#endif
	{
#if COM16_C806_EMT
#if JVET_C0024_QTBT
		if (bLuma)
		{
			m_pcEntropyCoder->encodeEmtCuFlag(pcCU, uiAbsPartIdx, pcCU->getDepth(0), pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, 0) ? true : false);
		}
#else
		if (0 == uiTrDepth)
		{
			m_pcEntropyCoder->encodeEmtCuFlag(pcCU, uiAbsPartIdx, pcCU->getDepth(0), pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrMode) ? true : false);
		}
#endif
#endif
		//===== Cbfs =====
		if (bLuma)
		{
			m_pcEntropyCoder->encodeQtCbf(rTu, COMPONENT_Y, true);
		}
	}
}




Void
TEncSearch::xEncCoeffQT(TComTU &rTu,
	const ComponentID  component,
	Bool         bRealCoeff)
{
	TComDataCU* pcCU = rTu.getCU();
	const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#if JVET_C0024_QTBT
	assert(uiAbsPartIdx == 0);
#endif
	const UInt uiTrDepth = rTu.GetTransformDepthRel();

#if JVET_C0024_QTBT
	const UInt  uiTrMode = 0;
	assert(uiTrDepth == 0);
#else
	const UInt  uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
	const UInt  uiSubdiv = (uiTrMode > uiTrDepth ? 1 : 0);

	if (uiSubdiv)
	{
		TComTURecurse tuRecurseChild(rTu, false);
		do
		{
			xEncCoeffQT(tuRecurseChild, component, bRealCoeff);
		} while (tuRecurseChild.nextSection(rTu));
	}
	else
#endif
		if (rTu.ProcessComponentSection(component))
		{
			//===== coefficients =====
#if JVET_C0024_QTBT
			UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
			UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
			const UInt  uiLog2TrafoSize = rTu.GetLog2LumaTrSize();
#endif
			UInt    uiCoeffOffset = rTu.getCoefficientOffset(component);
#if JVET_C0024_QTBT
			TCoeff* pcCoeff = bRealCoeff ? pcCU->getCoeff(component) : m_pppcQTTempCoeff[component][uiWIdx][uiHIdx];
#else
			UInt    uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrafoSize;
			TCoeff* pcCoeff = bRealCoeff ? pcCU->getCoeff(component) : m_ppcQTTempCoeff[component][uiQTLayer];
#endif

			if (isChroma(component) && (pcCU->getCbf(rTu.GetAbsPartIdxTU(), COMPONENT_Y, uiTrMode) != 0) && pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
			{
				m_pcEntropyCoder->encodeCrossComponentPrediction(rTu, component);
			}

			m_pcEntropyCoder->encodeCoeffNxN(rTu, pcCoeff + uiCoeffOffset, component);
		}
}




Void
TEncSearch::xEncIntraHeader(TComDataCU*  pcCU,
	UInt         uiTrDepth,
	UInt         uiAbsPartIdx,
	Bool         bLuma,
	Bool         bChroma)
{
	if (bLuma)
	{
		// CU header
		if (uiAbsPartIdx == 0)
		{
			if (!pcCU->getSlice()->isIntra())
			{
				if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
				{
					m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
				}
				m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
				m_pcEntropyCoder->encodePredMode(pcCU, 0, true);
#if VCEG_AZ05_INTRA_MPI 
				m_pcEntropyCoder->encodeMPIIdx(pcCU, 0, true);
#endif
			}
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
			m_pcEntropyCoder->encodePDPCIdx(pcCU, 0, true);
#endif

#if DIMD_INTRA_PRED
			m_pcEntropyCoder->encodeDIMDFlag(pcCU, 0, pcCU->getDepth(0), pcCU->getWidth(0), pcCU->getHeight(0), true);
#endif

#if JVET_C0024_QTBT
			if (pcCU->isIntra(0))
#else
			m_pcEntropyCoder->encodePartSize(pcCU, 0, pcCU->getDepth(0), true);

			if (pcCU->isIntra(0) && pcCU->getPartitionSize(0) == SIZE_2Nx2N)
#endif
			{
				m_pcEntropyCoder->encodeIPCMInfo(pcCU, 0, true);

				if (pcCU->getIPCMFlag(0))
				{
					return;
				}
			}
		}
		// luma prediction mode
#if !JVET_C0024_QTBT
		if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
		{
#endif
			if (uiAbsPartIdx == 0)
			{
				m_pcEntropyCoder->encodeIntraDirModeLuma(pcCU, 0);
			}
#if !JVET_C0024_QTBT
		}
		else
		{
			UInt uiQNumParts = pcCU->getTotalNumPart() >> 2;
			if (uiTrDepth>0 && (uiAbsPartIdx%uiQNumParts) == 0)
			{
				m_pcEntropyCoder->encodeIntraDirModeLuma(pcCU, uiAbsPartIdx);
			}
		}
#endif
	}

	if (bChroma)
	{
#if !JVET_C0024_QTBT
		if (pcCU->getPartitionSize(0) == SIZE_2Nx2N || !enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()))
		{
#endif
			if (uiAbsPartIdx == 0)
			{
				m_pcEntropyCoder->encodeIntraDirModeChroma(pcCU, uiAbsPartIdx);
			}
#if !JVET_C0024_QTBT
		}
		else
		{
			UInt uiQNumParts = pcCU->getTotalNumPart() >> 2;
			assert(uiTrDepth>0);
			if ((uiAbsPartIdx%uiQNumParts) == 0)
			{
				m_pcEntropyCoder->encodeIntraDirModeChroma(pcCU, uiAbsPartIdx);
			}
		}
#endif
	}
}




UInt
TEncSearch::xGetIntraBitsQT(TComTU &rTu,
	Bool         bLuma,
	Bool         bChroma,
	Bool         bRealCoeff /* just for test */)
{
	TComDataCU* pcCU = rTu.getCU();
	const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
	const UInt uiTrDepth = rTu.GetTransformDepthRel();
	m_pcEntropyCoder->resetBits();
	xEncIntraHeader(pcCU, uiTrDepth, uiAbsPartIdx, bLuma, bChroma);
	xEncSubdivCbfQT(rTu, bLuma, bChroma);

	if (bLuma)
	{
		xEncCoeffQT(rTu, COMPONENT_Y, bRealCoeff);
	}
	if (bChroma)
	{
		xEncCoeffQT(rTu, COMPONENT_Cb, bRealCoeff);
		xEncCoeffQT(rTu, COMPONENT_Cr, bRealCoeff);
	}
	UInt   uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();

	return uiBits;
}

UInt TEncSearch::xGetIntraBitsQTChroma(TComTU &rTu,
	ComponentID compID,
	Bool         bRealCoeff /* just for test */)
{
	m_pcEntropyCoder->resetBits();
	xEncCoeffQT(rTu, compID, bRealCoeff);
	UInt   uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
	return uiBits;
}

Void TEncSearch::xIntraCodingTUBlock(TComYuv*    pcOrgYuv,
	TComYuv*    pcPredYuv,
	TComYuv*    pcResiYuv,
#if COM16_C806_LARGE_CTU
	Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
	Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
	const Bool        checkCrossCPrediction,
	Distortion& ruiDist,
	const ComponentID compID,
	TComTU&     rTu
	DEBUG_STRING_FN_DECLARE(sDebug)
	, Int        default0Save1Load2
#if COM16_C806_EMT
	, UInt*      puiSigNum
#endif
)
{
	if (!rTu.ProcessComponentSection(compID))
	{
		return;
	}
	const Bool           bIsLuma = isLuma(compID);
	const TComRectangle &rect = rTu.getRect(compID);
	TComDataCU    *pcCU = rTu.getCU();
	const UInt           uiAbsPartIdx = rTu.GetAbsPartIdxTU();
	const TComSPS       &sps = *(pcCU->getSlice()->getSPS());
#if BARC
	Bool          isBARC = pcCU->getBARCFlag(toChannelType(compID));
#endif
#if !JVET_C0024_QTBT
	const UInt           uiTrDepth = rTu.GetTransformDepthRelAdj(compID);
	const UInt           uiFullDepth = rTu.GetTransformDepthTotal();
	const UInt           uiLog2TrSize = rTu.GetLog2LumaTrSize();
#endif
	const ChromaFormat   chFmt = pcOrgYuv->getChromaFormat();
	const ChannelType    chType = toChannelType(compID);
	const Int            bitDepth = sps.getBitDepth(chType);

	const UInt           uiWidth = rect.width;
	const UInt           uiHeight = rect.height;

#if JVET_C0024_QTBT
	if (isChroma(chType))
	{
		assert(uiWidth == pcCU->getWidth(0) >> 1);
	}
	const UInt           uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
	const UInt           uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#endif
	const UInt           uiStride = pcOrgYuv->getStride(compID);
	Pel           *piOrg = pcOrgYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piPred = pcPredYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piResi = pcResiYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piReco = pcPredYuv->getAddr(compID, uiAbsPartIdx);
#if JVET_C0024_QTBT
	assert(uiAbsPartIdx == 0);
	Pel           *piRecQt = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
	const UInt           uiRecQtStride = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
	const UInt           uiQTLayer = sps.getQuadtreeTULog2MaxSize() - uiLog2TrSize;
	Pel           *piRecQt = m_pcQTTempTComYuv[uiQTLayer].getAddr(compID, uiAbsPartIdx);
	const UInt           uiRecQtStride = m_pcQTTempTComYuv[uiQTLayer].getStride(compID);
#endif
	const UInt           uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
	Pel           *piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
	UInt           uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
#if JVET_C0024_QTBT
	assert(rTu.getCoefficientOffset(compID) == 0);
	TCoeff        *pcCoeff = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx];
#else
	TCoeff        *pcCoeff = m_ppcQTTempCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
	Bool           useTransformSkip = pcCU->getTransformSkip(uiAbsPartIdx, compID);

#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
	TCoeff        *pcArlCoeff = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx];
#else
	TCoeff        *pcArlCoeff = m_ppcQTTempArlCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
#endif

	const UInt           uiChPredMode = pcCU->getIntraDir(chType, uiAbsPartIdx);
#if JVET_C0024_QTBT
#if JVET_E0062_MULTI_DMS
	const UInt uiChCodedMode = uiChPredMode;
#else
	const UInt           uiChCodedMode = (uiChPredMode == DM_CHROMA_IDX && !bIsLuma)
		? (pcCU->getSlice()->isIntra() ? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiAbsPartIdx)
			: pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx)) : uiChPredMode;
#endif
#else
	const UInt           partsPerMinCU = 1 << (2 * (sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
	const UInt           uiChCodedMode = (uiChPredMode == DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
#endif
	const UInt           uiChFinalMode = ((chFmt == CHROMA_422) && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;

	const Int            blkX = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
	const Int            blkY = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
	const Int            bufferOffset = blkX + (blkY * MAX_CU_SIZE);
	Pel  *const    encoderLumaResidual = resiLuma[RESIDUAL_ENCODER_SIDE] + bufferOffset;
	Pel  *const    reconstructedLumaResidual = resiLuma[RESIDUAL_RECONSTRUCTED] + bufferOffset;
#if JVET_E0062_MULTI_DMS
#if JVET_E0077_ENHANCED_LM
	const Bool           bUseCrossCPrediction = isChroma(compID) && !IsLMMode(uiChPredMode) && checkCrossCPrediction;
#else
	const Bool           bUseCrossCPrediction = isChroma(compID) && (uiChPredMode != LM_CHROMA_IDX) && checkCrossCPrediction;
#endif
#else
	const Bool           bUseCrossCPrediction = isChroma(compID) && (uiChPredMode == DM_CHROMA_IDX) && checkCrossCPrediction;
#endif
	const Bool           bUseReconstructedResidualForEstimate = m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate();
	Pel *const     lumaResidualForEstimate = bUseReconstructedResidualForEstimate ? reconstructedLumaResidual : encoderLumaResidual;

#if DEBUG_STRING
	const Int debugPredModeMask = DebugStringGetPredModeMask(MODE_INTRA);
#endif

#if COM16_C806_EMT
	UChar   ucTrIdx = (pcCU->getEmtCuFlag(uiAbsPartIdx) && compID == COMPONENT_Y) ? pcCU->getEmtTuIdx(uiAbsPartIdx) : (pcCU->getSlice()->getSPS()->getUseIntraEMT() ? DCT2_EMT : DCT2_HEVC);
#endif
#if RRC_INTRA_SHARPENING
	Bool doSharpening = pcCU->getIntraSharpFlag(uiAbsPartIdx);
#if !INTRA_SHARP_CHROMA
	doSharpening = (compID == COMPONENT_Y) ? doSharpening : 0;
#elif JVET_C0024_QTBT
	if (isChroma(pcCU->getTextType()))
	{
		UInt absPartIdx = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
		absPartIdx = g_auiRasterToZscan[g_auiZscanToRaster[absPartIdx] + (pcCU->getHeight(uiAbsPartIdx) / pcCU->getPic()->getMinCUHeight()) / 2 * pcCU->getPic()->getNumPartInCtuWidth() + (pcCU->getWidth(uiAbsPartIdx) / pcCU->getPic()->getMinCUWidth()) / 2];
		doSharpening = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->isIntra(absPartIdx) ? pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraSharpFlag(absPartIdx) : 0;
		if (uiChFinalMode <= DC_IDX || uiChFinalMode >= LM_CHROMA_IDX)
			doSharpening = 0;
	}
#endif
#endif

	//===== init availability pattern =====
	DEBUG_STRING_NEW(sTemp)

#if JVET_E0077_LM_MF
		if ((uiChFinalMode >= LM_CHROMA_F1_IDX) && (uiChFinalMode < LM_CHROMA_F1_IDX + LM_FILTER_NUM))
		{
			Pel* piOrgPred = m_pLMMFPredSaved[(uiChFinalMode - LM_CHROMA_F1_IDX) * 2 + compID - COMPONENT_Cb];
			Pel* piSavPred = piPred;

			twoDimBufCopier
				(
					(UChar *)piSavPred,
					uiWidth * sizeof(Pel),
					(UChar *)piOrgPred,
					uiStride * sizeof(Pel),
					uiWidth * sizeof(Pel),
					uiHeight
					);
			//for (Int i = 0; i < uiHeight; i++)
			//{
			//	memcpy(piSavPred, piOrgPred, sizeof(Pel)* uiWidth);
			//	piOrgPred += uiWidth;
			//	piSavPred += uiStride;
			//}
		}
		else
#endif
#if !DEBUG_STRING
			if (default0Save1Load2 != 2)
#endif
			{
#if !RRC_INTRA_SHARPENING
				const Bool bUseFilteredPredictions = TComPrediction::filteringIntraReferenceSamples(compID, uiChFinalMode, uiWidth, uiHeight, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag());
#else
#if DIMD_NUM_INTRA_DIR_INC
				const Bool bUseFilteredPredictions = TComPrediction::filteringIntraReferenceSamples(compID, uiChFinalMode, uiWidth, uiHeight, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag(), pcCU, uiAbsPartIdx);
#else
				const Bool bUseFilteredPredictions = doSharpening || TComPrediction::filteringIntraReferenceSamples(compID, uiChFinalMode, uiWidth, uiHeight, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag());
#endif  
#endif
#if BARC_CORE
				Bool isCTUBoundary = 0;
				if (pcCU->getBARCFlag(pcCU->getTextType()))
				{
					UInt leftBoundary = (pcCU->getCUPelX() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
					UInt aboveBoundary = (pcCU->getCUPelY() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
					isCTUBoundary = leftBoundary + aboveBoundary > 0 ? 1 : 0;
				}
#endif
				initIntraPatternChType(rTu, compID, bUseFilteredPredictions
#if RRC_INTRA_SHARPENING
					, doSharpening
#endif
#if BARC_CORE
					, isCTUBoundary
#endif
					DEBUG_STRING_PASS_INTO(sDebug));

				//===== get prediction signal =====
#if COM16_C806_LMCHROMA
				if (uiChFinalMode == LM_CHROMA_IDX
#if JVET_E0077_ENHANCED_LM
					|| IsLMMode(uiChFinalMode)
#endif
					)
				{
#if !JVET_E0077_ENHANCED_LM || !JVET_C0024_QTBT
					getLumaRecPixels(rTu, uiWidth, uiHeight);
#endif

#if JVET_E0077_ENHANCED_LM
					Int iLMType = 0;
					iLMType = uiChFinalMode;
					predLMIntraChroma(rTu, compID, piPred, uiStride, uiWidth, uiHeight, iLMType);
#else
					predLMIntraChroma(rTu, compID, piPred, uiStride, uiWidth, uiHeight);
#endif

				}
				else
				{
#endif

					predIntraAng(compID, uiChFinalMode, piOrg, uiStride, piPred, uiStride, rTu, bUseFilteredPredictions);

#if COM16_C806_LMCHROMA
					if (compID == COMPONENT_Cr && pcCU->getSlice()->getSPS()->getUseLMChroma())
					{
						addCrossColorResi(rTu, compID, piPred, uiStride, uiWidth, uiHeight, pcResiYuv->getAddr(COMPONENT_Cb, uiAbsPartIdx), pcResiYuv->getStride(COMPONENT_Cb));
					}
				}
#endif
				// save prediction
				if (default0Save1Load2 == 1)
				{
					Pel*  pPred = piPred;
					Pel*  pPredBuf = m_pSharedPredTransformSkip[compID];
					//Int k = 0;

					twoDimBufCopier
						(
							(UChar *)pPredBuf,
							uiWidth * sizeof(Pel),
							(UChar *)pPred,
							uiStride * sizeof(Pel),
							uiWidth * sizeof(Pel),
							uiHeight
							);
					//for (UInt uiY = 0; uiY < uiHeight; uiY++)
					//{
					//	for (UInt uiX = 0; uiX < uiWidth; uiX++)
					//	{
					//		pPredBuf[k++] = pPred[uiX];
					//	}
					//	pPred += uiStride;
					//}
				}
			}
#if !DEBUG_STRING
			else
			{
				// load prediction
				Pel*  pPred = piPred;
				Pel*  pPredBuf = m_pSharedPredTransformSkip[compID];
				Int k = 0;

				twoDimBufCopier
					(
						(UChar *)pPred,
						uiStride * sizeof(Pel),
						(UChar *)pPredBuf,
						uiWidth * sizeof(Pel),
						uiWidth * sizeof(Pel),
						uiHeight
						);
				//for (UInt uiY = 0; uiY < uiHeight; uiY++)
				//{
				//	for (UInt uiX = 0; uiX < uiWidth; uiX++)
				//	{
				//		pPred[uiX] = pPredBuf[k++];
				//	}
				//	pPred += uiStride;
				//}
			}
#endif

	//===== get residual signal =====
	{
		// get residual
		Pel*  pOrg = piOrg;
		Pel*  pPred = piPred;
		Pel*  pResi = piResi;

		residueComputer
			(
				piOrg,
				uiStride,
				piPred,
				uiStride,
				piResi,
				uiStride,
				uiWidth,
				uiHeight
				);
		//for (UInt uiY = 0; uiY < uiHeight; uiY++)
		//{
		//	for (UInt uiX = 0; uiX < uiWidth; uiX++)
		//	{
		//		pResi[uiX] = pOrg[uiX] - pPred[uiX];
		//	}

		//	pOrg += uiStride;
		//	pResi += uiStride;
		//	pPred += uiStride;
		//}
	}

	if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
	{
		if (bUseCrossCPrediction)
		{
			if (xCalcCrossComponentPredictionAlpha(rTu, compID, lumaResidualForEstimate, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride) == 0)
			{
				return;
			}
			TComTrQuant::crossComponentPrediction(rTu, compID, reconstructedLumaResidual, piResi, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride, uiStride, false);
		}
		else if (isLuma(compID) && !bUseReconstructedResidualForEstimate)
		{
			xStoreCrossComponentPredictionResult(encoderLumaResidual, piResi, rTu, 0, 0, MAX_CU_SIZE, uiStride);
		}
	}

	//===== transform and quantization =====
	//--- init rate estimation arrays for RDOQ ---
#if COM16_C806_EMT
	if ((useTransformSkip ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ()) && !(ucTrIdx >= 1 && ucTrIdx <= 3))
#else
	if (useTransformSkip ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ())
#endif
	{
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
		COEFF_SCAN_TYPE scanType = COEFF_SCAN_TYPE(pcCU->getCoefScanIdx(uiAbsPartIdx, uiWidth, uiHeight, compID));
#endif
		m_pcEntropyCoder->estimateBit(m_pcTrQuant->m_pcEstBitsSbac, uiWidth, uiHeight, chType
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
			, scanType
#endif
		);
	}

	//--- transform and quantization ---
	TCoeff uiAbsSum = 0;
#if !JVET_C0024_QTBT
	if (bIsLuma)
	{
		pcCU->setTrIdxSubParts(uiTrDepth, uiAbsPartIdx, uiFullDepth);
	}
#endif

	const QpParam cQP(*pcCU, compID);

#if RDOQ_CHROMA_LAMBDA
	m_pcTrQuant->selectLambda(compID);
#endif

#if COM16_C806_CR_FROM_CB_LAMBDA_ADJUSTMENT
	if (
#if JVET_E0077_ENHANCED_LM
		!IsLMMode(uiChFinalMode) &&
#endif 
		uiChFinalMode != LM_CHROMA_IDX && pcCU->getSlice()->getSPS()->getUseLMChroma())
	{
		if (compID == COMPONENT_Cb)
		{
			m_pcTrQuant->setLambda(m_pcTrQuant->getlambda() * 15 / 16);
		}
		else if (compID == COMPONENT_Cr)
		{
			m_pcTrQuant->setLambda(m_pcTrQuant->getlambda() * 16 / 15);
		}
	}
#endif

#if PVC_AQS
	m_pcTrQuant->setEstQSFactor(piPred, uiStride, uiWidth, uiHeight, pcCU->getSlice()->getSliceType(), pcCU, uiAbsPartIdx, compID);
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING_ENC_METHOD
	smoothResidual(compID, piOrg, piResi, uiHeight, uiWidth, uiStride, uiStride);
#endif
	m_pcTrQuant->transformNxN(rTu, compID, piResi, uiStride, pcCoeff,
#if ADAPTIVE_QP_SELECTION
		pcArlCoeff,
#endif
		uiAbsSum, cQP
	);

#if COM16_C806_EMT
	if (ucTrIdx != DCT2_EMT && ucTrIdx != DCT2_HEVC)
	{
		*puiSigNum = 0;
		for (UInt uiX = 0; uiX < uiHeight*uiWidth; uiX++)
		{
			if (pcCoeff[uiX])
			{
				(*puiSigNum)++;
				if (*puiSigNum>g_iEmtSigNumThr)
				{
					break;
				}
			}
		}

		if (ucTrIdx != 0 && *puiSigNum <= g_iEmtSigNumThr && !useTransformSkip)
		{
			return;
		}
	}
#endif

	//--- inverse transform ---

#if DEBUG_STRING
	if ((uiAbsSum > 0) || (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask))
#else
	if (uiAbsSum > 0)
#endif
	{
		m_pcTrQuant->invTransformNxN(rTu, compID, piResi, uiStride, pcCoeff, cQP DEBUG_STRING_PASS_INTO_OPTIONAL(&sDebug, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask)));
	}
	else
	{
		Pel* pResi = piResi;
		memset(pcCoeff, 0, sizeof(TCoeff) * uiWidth * uiHeight);
		for (UInt uiY = 0; uiY < uiHeight; uiY++)
		{
			memset(pResi, 0, sizeof(Pel) * uiWidth);
			pResi += uiStride;
		}
	}


	//===== reconstruction =====
	{
		Pel* pPred = piPred;
		Pel* pResi = piResi;
		Pel* pReco = piReco;
		Pel* pRecQt = piRecQt;
		Pel* pRecIPred = piRecIPred;

		if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
		{
			if (bUseCrossCPrediction)
			{
				TComTrQuant::crossComponentPrediction(rTu, compID, reconstructedLumaResidual, piResi, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride, uiStride, true);
			}
			else if (isLuma(compID))
			{
				xStoreCrossComponentPredictionResult(reconstructedLumaResidual, piResi, rTu, 0, 0, MAX_CU_SIZE, uiStride);
			}
		}

#if DEBUG_STRING
		std::stringstream ss(stringstream::out);
		const Bool bDebugPred = ((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
		const Bool bDebugResi = ((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
		const Bool bDebugReco = ((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));

		if (bDebugPred || bDebugResi || bDebugReco)
		{
			ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << "\n";
			for (UInt uiY = 0; uiY < uiHeight; uiY++)
			{
				ss << "###: ";
				if (bDebugPred)
				{
					ss << " - pred: ";
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						ss << pPred[uiX] << ", ";
					}
				}
				if (bDebugResi)
				{
					ss << " - resi: ";
				}
				for (UInt uiX = 0; uiX < uiWidth; uiX++)
				{
					if (bDebugResi)
					{
						ss << pResi[uiX] << ", ";
					}
#if JVET_D0033_ADAPTIVE_CLIPPING
					pReco[uiX] = Pel(ClipA<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), compID));
#else
					pReco[uiX] = Pel(ClipBD<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth));
#endif
					pRecQt[uiX] = pReco[uiX];
					pRecIPred[uiX] = pReco[uiX];
				}
				if (bDebugReco)
				{
					ss << " - reco: ";
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						ss << pReco[uiX] << ", ";
					}
				}
				pPred += uiStride;
				pResi += uiStride;
				pReco += uiStride;
				pRecQt += uiRecQtStride;
				pRecIPred += uiRecIPredStride;
				ss << "\n";
			}
			DEBUG_STRING_APPEND(sDebug, ss.str())
		}
		else
#endif
		{
			Int iMin = g_ClipParam.min(compID);
			Int iMax = g_ClipParam.max(compID);

			reconComputer
				(
					pResi,
					uiStride,
					pPred,
					uiStride,
					pReco,
					uiStride,
					iMin,
					iMax,
					uiWidth,
					uiHeight
					);

			twoDimBufCopier
				(
					(UChar *)pRecQt,
					uiRecQtStride * sizeof(Pel),
					(UChar *)pReco,
					uiStride * sizeof(Pel),
					uiWidth * sizeof(Pel),
					uiHeight
					);

			twoDimBufCopier
				(
					(UChar *)pRecIPred,
					uiRecIPredStride * sizeof(Pel),
					(UChar *)pReco,
					uiStride * sizeof(Pel),
					uiWidth * sizeof(Pel),
					uiHeight
					);
//			for (UInt uiY = 0; uiY < uiHeight; uiY++)
//			{
//				for (UInt uiX = 0; uiX < uiWidth; uiX++)
//				{
//#if JVET_D0033_ADAPTIVE_CLIPPING // encoder intrac rec
//					pReco[uiX] = Pel(ClipA<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), compID));
//#else
//					pReco[uiX] = Pel(ClipBD<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth));
//#endif
//					pRecQt[uiX] = pReco[uiX];
//					pRecIPred[uiX] = pReco[uiX];
//				}
//				pPred += uiStride;
//				pResi += uiStride;
//				pReco += uiStride;
//				pRecQt += uiRecQtStride;
//				pRecIPred += uiRecIPredStride;
//			}
		}
	}

#if JVET_F0096_BILATERAL_FILTER
	if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
	{
		if (isLuma(compID))
		{
			if (uiAbsSum && (pcCU->getQP(COMPONENT_Y) > 17))
			{
				TComBilateralFilter::instance()->bilateralFilterIntra(pcCU, uiWidth, uiHeight, piReco, uiStride, pcCU->getQP(COMPONENT_Y));
				
				twoDimBufCopier
					(
						(UChar *)piRecQt,
						uiRecQtStride * sizeof(Short),
						(UChar *)piReco,
						uiStride * sizeof(Short),
						uiWidth * sizeof(Short),
						uiHeight
						);

				twoDimBufCopier
					(
						(UChar *)piRecIPred,
						uiRecIPredStride * sizeof(Short),
						(UChar *)piReco,
						uiStride * sizeof(Short),
						uiWidth * sizeof(Short),
						uiHeight
						);
				//for (UInt uiY = 0; uiY < uiHeight; uiY++)
				//{
				//	memcpy(piRecQt + uiY * uiRecQtStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//	memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//	uiY++;
				//	memcpy(piRecQt + uiY * uiRecQtStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//	memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//}
			}
		}
	}
#endif
	//===== update distortion =====
#if WCG_LUMA_DQP_CM_SCALE
	if (m_pcEncCfg->getUseLumaDeltaQp() > 0) {
		UInt           iOrgStrideLuma = pcOrgYuv->getStride(COMPONENT_Y);
		Pel           *piOrgLuma = pcOrgYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);
		ruiDist += m_pcRdCost->getDistPart(bitDepth, piReco, uiStride, piOrg, uiStride, uiWidth, uiHeight, compID, DF_SSE_WTD, piOrgLuma, iOrgStrideLuma);  // use weighted SSE
	}
	else
#endif
		ruiDist += m_pcRdCost->getDistPart(bitDepth, piReco, uiStride, piOrg, uiStride, uiWidth, uiHeight, compID);
}

#if VCEG_AZ08_INTRA_KLT
Bool TEncSearch::xIntraCodingTUBlockTM(TComYuv*    pcOrgYuv,
	TComYuv*    pcPredYuv,
	TComYuv*    pcResiYuv,
	Distortion& ruiDist,
	const ComponentID compID,
	TComTU&     rTu
	DEBUG_STRING_FN_DECLARE(sDebug)
#if COM16_C806_EMT
	, UInt*      puiSigNum
#endif
	, Int tmpred0_tmpredklt1_ori2
)
{
	if (!rTu.ProcessComponentSection(compID))
	{
		return false;
	}
	const Bool           bIsLuma = isLuma(compID);
	const TComRectangle &rect = rTu.getRect(compID);
	TComDataCU    *pcCU = rTu.getCU();
	const UInt           uiAbsPartIdx = rTu.GetAbsPartIdxTU();
	const TComSPS       &sps = *(pcCU->getSlice()->getSPS());

#if JVET_C0024_QTBT
	UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
	UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
	assert(uiWIdx == uiHIdx);
	//const UInt           uiLog2TrSize = ((uiWIdx+uiHIdx)>>1) + MIN_CU_LOG2; 
#else
	const UInt           uiTrDepth = rTu.GetTransformDepthRelAdj(compID);
	const UInt           uiFullDepth = rTu.GetTransformDepthTotal();
	const UInt           uiLog2TrSize = rTu.GetLog2LumaTrSize();
#endif
	const ChannelType    chType = toChannelType(compID);
	const Int            bitDepth = sps.getBitDepth(chType);

	const UInt           uiWidth = rect.width;
	const UInt           uiHeight = rect.height;
	const UInt           uiStride = pcOrgYuv->getStride(compID);
	Pel           *piOrg = pcOrgYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piPred = pcPredYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piResi = pcResiYuv->getAddr(compID, uiAbsPartIdx);
	Pel           *piReco = pcPredYuv->getAddr(compID, uiAbsPartIdx);
#if JVET_C0024_QTBT
	Pel           *piRecQt = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
	const UInt           uiRecQtStride = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
	const UInt           uiQTLayer = sps.getQuadtreeTULog2MaxSize() - uiLog2TrSize;
	Pel           *piRecQt = m_pcQTTempTComYuv[uiQTLayer].getAddr(compID, uiAbsPartIdx);
	const UInt           uiRecQtStride = m_pcQTTempTComYuv[uiQTLayer].getStride(compID);
#endif
	const UInt           uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
	Pel           *piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
	UInt           uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
#if JVET_C0024_QTBT
	TCoeff        *pcCoeff = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
	TCoeff        *pcCoeff = m_ppcQTTempCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
	Bool           useTransformSkip = pcCU->getTransformSkip(uiAbsPartIdx, compID);

#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
	TCoeff        *pcArlCoeff = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
	TCoeff        *pcArlCoeff = m_ppcQTTempArlCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
#endif

#if DEBUG_STRING
	const Int debugPredModeMask = DebugStringGetPredModeMask(MODE_INTRA);
#endif

#if COM16_C806_EMT
	UChar   ucTrIdx = (pcCU->getEmtCuFlag(uiAbsPartIdx) && compID == COMPONENT_Y) ? pcCU->getEmtTuIdx(uiAbsPartIdx) : (pcCU->getSlice()->getSPS()->getUseIntraEMT() ? DCT2_EMT : DCT2_HEVC);
#endif
	//===== init availability pattern =====
	DEBUG_STRING_NEW(sTemp)

#if VCEG_AZ08_INTRA_KLT
		Bool useKLT = false;
	if (tmpred0_tmpredklt1_ori2 != 2 && bIsLuma)
	{
		UInt uiBlkSize = uiWidth;
		UInt uiTarDepth = g_aucConvertToBit[uiBlkSize];
		UInt uiTempSize = g_uiDepth2IntraTempSize[uiTarDepth];
		m_pcTrQuant->getTargetTemplate(pcCU, uiAbsPartIdx, uiBlkSize, uiTempSize);
		m_pcTrQuant->candidateSearchIntra(pcCU, uiAbsPartIdx, uiBlkSize, uiTempSize);

		Int foundCandiNum;
		Bool bSuccessful = m_pcTrQuant->generateTMPrediction(piPred, uiStride, uiBlkSize, uiTempSize, foundCandiNum);
		if (bSuccessful == false || foundCandiNum < 1)
		{
			return false;
		}
		if (1 == tmpred0_tmpredklt1_ori2 && bSuccessful)
		{
			useKLT = m_pcTrQuant->calcKLTIntra(piPred, uiStride, uiBlkSize);
		}
	}
#endif

	//===== get residual signal =====
	{
		// get residual
		Pel*  pOrg = piOrg;
		Pel*  pPred = piPred;
		Pel*  pResi = piResi;

		for (UInt uiY = 0; uiY < uiHeight; uiY++)
		{
			for (UInt uiX = 0; uiX < uiWidth; uiX++)
			{
				pResi[uiX] = pOrg[uiX] - pPred[uiX];
			}

			pOrg += uiStride;
			pResi += uiStride;
			pPred += uiStride;
		}
	}

	//===== transform and quantization =====
	//--- init rate estimation arrays for RDOQ ---
#if COM16_C806_EMT
	if ((useTransformSkip ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ()) && !(ucTrIdx >= 1 && ucTrIdx <= 3))
#else
	if (useTransformSkip ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ())
#endif
	{
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
		COEFF_SCAN_TYPE scanType = COEFF_SCAN_TYPE(pcCU->getCoefScanIdx(uiAbsPartIdx, uiWidth, uiHeight, compID));
#endif
		m_pcEntropyCoder->estimateBit(m_pcTrQuant->m_pcEstBitsSbac, uiWidth, uiHeight, chType
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
			, scanType
#endif
		);
	}

	//--- transform and quantization ---
	TCoeff uiAbsSum = 0;
#if !JVET_C0024_QTBT
	if (bIsLuma)
	{
		pcCU->setTrIdxSubParts(uiTrDepth, uiAbsPartIdx, uiFullDepth);
	}
#endif

	const QpParam cQP(*pcCU, compID);

#if RDOQ_CHROMA_LAMBDA
	m_pcTrQuant->selectLambda(compID);
#endif
	m_pcTrQuant->transformNxN(rTu, compID, piResi, uiStride, pcCoeff,
#if ADAPTIVE_QP_SELECTION
		pcArlCoeff,
#endif
		uiAbsSum, cQP
#if VCEG_AZ08_INTRA_KLT
		, useKLT
#endif
	);

#if COM16_C806_EMT
	if (ucTrIdx != DCT2_EMT && ucTrIdx != DCT2_HEVC)
	{
		*puiSigNum = 0;
		for (UInt uiX = 0; uiX < uiHeight*uiWidth; uiX++)
		{
			if (pcCoeff[uiX])
			{
				(*puiSigNum)++;
				if (*puiSigNum>g_iEmtSigNumThr)
				{
					break;
				}
			}
		}

		if (ucTrIdx != 0 && *puiSigNum <= g_iEmtSigNumThr && !useTransformSkip)
		{
			return false;
		}
	}
#endif

	//--- inverse transform ---

#if DEBUG_STRING
	if ((uiAbsSum > 0) || (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask))
#else
	if (uiAbsSum > 0)
#endif
	{
		const UInt *scan;
		if (useKLT)
		{
			TUEntropyCodingParameters codingParameters;
			getTUEntropyCodingParameters(codingParameters, rTu, compID);
			scan = codingParameters.scan;
			recoverOrderCoeff(pcCoeff, scan, uiWidth, uiHeight);
#if ADAPTIVE_QP_SELECTION
			recoverOrderCoeff(pcArlCoeff, scan, uiWidth, uiHeight);
#endif 
		}
		m_pcTrQuant->invTransformNxN(rTu, compID, piResi, uiStride, pcCoeff, cQP, useKLT DEBUG_STRING_PASS_INTO_OPTIONAL(&sDebug, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask)));
		if (useKLT)
		{
			reOrderCoeff(pcCoeff, scan, uiWidth, uiHeight);
#if ADAPTIVE_QP_SELECTION
			reOrderCoeff(pcArlCoeff, scan, uiWidth, uiHeight);
#endif 
		}
	}
	else
	{
		Pel* pResi = piResi;
		memset(pcCoeff, 0, sizeof(TCoeff)* uiWidth * uiHeight);
		for (UInt uiY = 0; uiY < uiHeight; uiY++)
		{
			memset(pResi, 0, sizeof(Pel)* uiWidth);
			pResi += uiStride;
		}
	}


	//===== reconstruction =====
	{
		Pel* pPred = piPred;
		Pel* pResi = piResi;
		Pel* pReco = piReco;
		Pel* pRecQt = piRecQt;
		Pel* pRecIPred = piRecIPred;

#if DEBUG_STRING
		std::stringstream ss(stringstream::out);
		const Bool bDebugPred = ((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
		const Bool bDebugResi = ((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
		const Bool bDebugReco = ((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));

		if (bDebugPred || bDebugResi || bDebugReco)
		{
			ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << "\n";
			for (UInt uiY = 0; uiY < uiHeight; uiY++)
			{
				ss << "###: ";
				if (bDebugPred)
				{
					ss << " - pred: ";
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						ss << pPred[uiX] << ", ";
					}
				}
				if (bDebugResi)
				{
					ss << " - resi: ";
				}
				for (UInt uiX = 0; uiX < uiWidth; uiX++)
				{
					if (bDebugResi)
					{
						ss << pResi[uiX] << ", ";
					}
#if JVET_D0033_ADAPTIVE_CLIPPING
					pReco[uiX] = Pel(ClipA<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), compID));
#else
					pReco[uiX] = Pel(ClipBD<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth));
#endif
					pRecQt[uiX] = pReco[uiX];
					pRecIPred[uiX] = pReco[uiX];
				}
				if (bDebugReco)
				{
					ss << " - reco: ";
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						ss << pReco[uiX] << ", ";
					}
				}
				pPred += uiStride;
				pResi += uiStride;
				pReco += uiStride;
				pRecQt += uiRecQtStride;
				pRecIPred += uiRecIPredStride;
				ss << "\n";
			}
			DEBUG_STRING_APPEND(sDebug, ss.str())
		}
		else
#endif
		{

			for (UInt uiY = 0; uiY < uiHeight; uiY++)
			{
				for (UInt uiX = 0; uiX < uiWidth; uiX++)
				{
#if JVET_D0033_ADAPTIVE_CLIPPING
					pReco[uiX] = Pel(ClipA<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), compID));
#else
					pReco[uiX] = Pel(ClipBD<Int>(Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth));
#endif
					pRecQt[uiX] = pReco[uiX];
					pRecIPred[uiX] = pReco[uiX];
				}
				pPred += uiStride;
				pResi += uiStride;
				pReco += uiStride;
				pRecQt += uiRecQtStride;
				pRecIPred += uiRecIPredStride;
			}
		}
	}
#if JVET_F0096_BILATERAL_FILTER
	if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
	{
		if (isLuma(compID))
		{
			if (uiAbsSum && (pcCU->getQP(COMPONENT_Y) > 17))
			{
				TComBilateralFilter::instance()->bilateralFilterIntra(pcCU, uiWidth, uiHeight, piReco, uiStride, pcCU->getQP(COMPONENT_Y));
				
				twoDimBufCopier
					(
						(UChar *)piRecQt,
						uiRecQtStride * sizeof(Short),
						(UChar *)piReco,
						uiStride * sizeof(Short),
						uiWidth * sizeof(Short),
						uiHeight
						);

				twoDimBufCopier
					(
						(UChar *)piRecIPred,
						uiRecIPredStride * sizeof(Short),
						(UChar *)piReco,
						uiStride * sizeof(Short),
						uiWidth * sizeof(Short),
						uiHeight
						);
				//for (UInt uiY = 0; uiY < uiHeight; uiY++)
				//{
				//	memcpy(piRecQt + uiY * uiRecQtStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//	memcpy(piRecIPred + uiY * uiRecIPredStride, piReco + uiY * uiStride, uiWidth * sizeof(Short));
				//}
			}
		}
	}
#endif
	//===== update distortion =====
#if WCG_LUMA_DQP_CM_SCALE
	if (m_pcEncCfg->getUseLumaDeltaQp() > 0) {
		UInt           iOrgStrideLuma = pcOrgYuv->getStride(COMPONENT_Y);
		Pel           *piOrgLuma = pcOrgYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);
		ruiDist += m_pcRdCost->getDistPart(bitDepth, piReco, uiStride, piOrg, uiStride, uiWidth, uiHeight, compID, DF_SSE_WTD, piOrgLuma, iOrgStrideLuma);  // use weighted SSE
	}
	else
#endif
		ruiDist += m_pcRdCost->getDistPart(bitDepth, piReco, uiStride, piOrg, uiStride, uiWidth, uiHeight, compID);

	return true;
}
#endif

Void
TEncSearch::xRecurIntraCodingLumaQT(TComYuv*    pcOrgYuv,
	TComYuv*    pcPredYuv,
	TComYuv*    pcResiYuv,
#if COM16_C806_LARGE_CTU
	Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
	Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
	Distortion& ruiDistY,
#if HHI_RQT_INTRA_SPEEDUP
	Bool        bCheckFirst,
#endif
	Double&     dRDCost,
	TComTU&     rTu
	DEBUG_STRING_FN_DECLARE(sDebug))
{
	TComDataCU   *pcCU = rTu.getCU();
	const UInt    uiAbsPartIdx = rTu.GetAbsPartIdxTU();
	const UInt    uiFullDepth = rTu.GetTransformDepthTotal();
	const UInt    uiTrDepth = rTu.GetTransformDepthRel();
#if BARC
	Bool    isBARC = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#if JVET_C0024_QTBT
	Bool bCheckFull = true;
	Bool bCheckSplit = false;
	assert(uiAbsPartIdx == 0 && uiTrDepth == 0);
	UInt uiWidth = pcCU->getWidth(uiAbsPartIdx);
	UInt uiHeight = pcCU->getHeight(uiAbsPartIdx);

	UInt uiWIdx = g_aucConvertToBit[uiWidth];
	UInt uiHIdx = g_aucConvertToBit[uiHeight];
#else
	const UInt    uiLog2TrSize = rTu.GetLog2LumaTrSize();
	Bool    bCheckFull = (uiLog2TrSize <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize());
	Bool    bCheckSplit = (uiLog2TrSize  >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx));
#endif

#if !COM16_C806_LARGE_CTU
	Pel     resiLumaSplit[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
	Pel     resiLumaSingle[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
#endif

	Bool    bMaintainResidual[NUMBER_OF_STORED_RESIDUAL_TYPES];
	for (UInt residualTypeIndex = 0; residualTypeIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; residualTypeIndex++)
	{
		bMaintainResidual[residualTypeIndex] = true; //assume true unless specified otherwise
	}

	bMaintainResidual[RESIDUAL_ENCODER_SIDE] = !(m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate());

#if COM16_C806_EMT
	UInt    uiSigNum;
#endif

#if HHI_RQT_INTRA_SPEEDUP
#if JVET_C0024_QTBT
	Int isIntraSlice = (pcCU->getSlice()->getSliceType() == I_SLICE);
#else
	Int maxTuSize = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize();
	Int isIntraSlice = (pcCU->getSlice()->getSliceType() == I_SLICE);
	// don't check split if TU size is less or equal to max TU size
	Bool noSplitIntraMaxTuSize = bCheckFull;
	if (m_pcEncCfg->getRDpenalty() && !isIntraSlice)
	{
		// in addition don't check split if TU size is less or equal to 16x16 TU size for non-intra slice
		noSplitIntraMaxTuSize = (uiLog2TrSize <= min(maxTuSize, 4));

		// if maximum RD-penalty don't check TU size 32x32
		if (m_pcEncCfg->getRDpenalty() == 2)
		{
			bCheckFull = (uiLog2TrSize <= min(maxTuSize, 4));
		}
	}
	if (bCheckFirst && noSplitIntraMaxTuSize)

	{
		bCheckSplit = false;
	}
#endif
#else
	Int maxTuSize = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize();
	Int isIntraSlice = (pcCU->getSlice()->getSliceType() == I_SLICE);
	// if maximum RD-penalty don't check TU size 32x32
	if ((m_pcEncCfg->getRDpenalty() == 2) && !isIntraSlice)
	{
		bCheckFull = (uiLog2TrSize <= min(maxTuSize, 4));
	}
#endif

#if COM16_C806_EMT && HHI_RQT_INTRA_SPEEDUP
	// Re-use the selected transform indexes in the previous call of xRecurIntraCodingQT
#if JVET_C0024_QTBT
	UInt    uiInitTrDepth = 0;
#else
	UInt    uiInitTrDepth = pcCU->getPartitionSize(0) == SIZE_2Nx2N ? 0 : 1;
#endif
	UChar   ucSavedEmtTrIdx = 0;
	Bool    bCheckInitTrDepth = false;
	static UInt uiInitAbsPartIdx;
	if (uiTrDepth == uiInitTrDepth)
	{
		uiInitAbsPartIdx = uiAbsPartIdx;
	}
	if (!bCheckFirst && uiTrDepth == uiInitTrDepth)
	{
		ucSavedEmtTrIdx = m_puhQTTempEmtTuIdx[uiAbsPartIdx - uiInitAbsPartIdx];
		bCheckInitTrDepth = true;
	}
#endif

#if JVET_C0024_QTBT
	assert(bCheckFull && !bCheckSplit);
#endif

	Double     dSingleCost = MAX_DOUBLE;
	Distortion uiSingleDistLuma = 0;
	UInt       uiSingleCbfLuma = 0;
	Bool       checkTransformSkip = pcCU->getSlice()->getPPS()->getUseTransformSkip();
	Int        bestModeId[MAX_NUM_COMPONENT] = { 0, 0, 0 };
#if COM16_C806_EMT
	UChar   bestTrIdx = 0;
#if FAST_EMT_v1
	UChar   nNumTrCands = pcCU->getEmtCuFlag(uiAbsPartIdx) ? 3 : 1;
#else
	UChar   nNumTrCands = pcCU->getEmtCuFlag(uiAbsPartIdx) ? 4 : 1;
#endif
	Bool    bAllIntra = (m_pcEncCfg->getIntraPeriod() == 1);
#endif

#if JVET_C0024_QTBT
	checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(pcCU->getSlice()->isIntra(), rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
#else
	checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
#endif
	checkTransformSkip &= (!pcCU->getCUTransquantBypass(0));
#if JVET_F0031_RMV_REDUNDANT_TRSKIP && COM16_C806_EMT
	checkTransformSkip &= (!pcCU->getEmtCuFlag(uiAbsPartIdx));
#endif

	assert(rTu.ProcessComponentSection(COMPONENT_Y));
#if !JVET_C0024_QTBT
	const UInt totalAdjustedDepthChan = rTu.GetTransformDepthTotalAdj(COMPONENT_Y);

	if (m_pcEncCfg->getUseTransformSkipFast())
	{
		checkTransformSkip &= (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_NxN);
	}
#endif


#if VCEG_AZ08_INTRA_KLT
#if !JVET_C0024_QTBT
	UInt  bestTMKLT = 0;
#endif
#if VCEG_AZ08_USE_KLT
	Bool checkTM = pcCU->getSlice()->getSPS()->getUseIntraKLT() && (bCheckFirst == false);
#else
	Bool checkTM = (bCheckFirst == false);
#endif
	if (checkTM)
	{
		UInt uiMaxTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
		UInt uiMinTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MIN - 1];
#if !JVET_C0024_QTBT
		UInt    uiWidth = pcCU->getWidth(0) >> uiTrDepth;
		UInt    uiHeight = pcCU->getHeight(0) >> uiTrDepth;
#endif
		Bool bCheckKLTFlag = (uiWidth == uiHeight) && (uiWidth <= uiMaxTrWidth) && (uiWidth >= uiMinTrWidth);
		checkTM &= bCheckKLTFlag;

#if VCEG_AZ05_INTRA_MPI //for speed up only
		Int iMPIidx = pcCU->getMPIIdx(uiAbsPartIdx);
		checkTM = checkTM & (iMPIidx == 0);
#endif
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
		Int iPdpcIdx = pcCU->getPDPCIdx(uiAbsPartIdx);
		checkTM = checkTM & (iPdpcIdx == 0);
#endif
#if VCEG_AZ05_ROT_TR || COM16_C1044_NSST //for speed up only
#if JVET_C0024_QTBT
		Int iRotidx = pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
#else
		Int iRotidx = pcCU->getROTIdx(uiAbsPartIdx);
#endif
		checkTM = checkTM & (iRotidx == 0);
#endif
	}
#endif

	if (bCheckFull)
	{
		if (checkTransformSkip == true)
		{
			//----- store original entropy coding status -----
#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
			m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif

			Distortion singleDistTmpLuma = 0;
			UInt       singleCbfTmpLuma = 0;
			Double     singleCostTmp = 0;
			Int        firstCheckId = 0;
#if VCEG_AZ08_INTRA_KLT
			pcCU->setKLTFlagSubParts(0, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#endif
			for (Int modeId = firstCheckId; modeId < 2; modeId++)
			{
#if COM16_C806_EMT
#if HHI_RQT_INTRA_SPEEDUP
				UChar numTrIdxCands = ((modeId == firstCheckId && !bCheckInitTrDepth) ? nNumTrCands : 1);
#else
				UChar numTrIdxCands = ((modeId == firstCheckId) ? nNumTrCands : 1);
#endif

				for (UChar ucTrIdx = 0; ucTrIdx < numTrIdxCands; ucTrIdx++)
				{
					// Skip checking other transform candidates if zero CBF is encountered
					if (ucTrIdx && !uiSingleCbfLuma && bAllIntra && m_pcEncCfg->getUseFastIntraEMT())
					{
						continue;
					}
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
					if (this->getSaveLoadTag(pcCU->getZorderIdxInCtu(), pcCU->getCURouteIdx(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && ucTrIdx && ucTrIdx != this->getSaveLoadEmtIdx(pcCU->getZorderIdxInCtu(), pcCU->getCURouteIdx(), uiWIdx, uiHIdx))
#else
					if (this->getSaveLoadTag(pcCU->getZorderIdxInCtu(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && ucTrIdx && ucTrIdx != this->getSaveLoadEmtIdx(uiWIdx, uiHIdx))
#endif
					{
						continue;
					}
#endif
#endif

					DEBUG_STRING_NEW(sModeString)
						Int  default0Save1Load2 = 0;
					singleDistTmpLuma = 0;
#if COM16_C806_EMT
					if (modeId == firstCheckId && ucTrIdx == 0)
#else
					if (modeId == firstCheckId)
#endif
					{
						default0Save1Load2 = 1;
					}
					else
					{
						default0Save1Load2 = 2;
					}

#if COM16_C806_EMT
#if HHI_RQT_INTRA_SPEEDUP
					pcCU->setEmtTuIdxSubParts(bCheckInitTrDepth ? ucSavedEmtTrIdx : ucTrIdx, uiAbsPartIdx, uiFullDepth);
#else
					pcCU->setEmtTuIdxSubParts(txIdx, uiAbsPartIdx, uiFullDepth);
#endif
#endif

#if JVET_C0024_QTBT
					pcCU->setTransformSkipSubParts(modeId, COMPONENT_Y, uiAbsPartIdx, 0);
#else
					pcCU->setTransformSkipSubParts(modeId, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
					xIntraCodingTUBlock(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
						m_resiSingleBuffer[uiWIdx][uiHIdx],
#else
						m_resiSingleBuffer[uiLog2TrSize],
#endif
#else
						resiLumaSingle,
#endif
						false, singleDistTmpLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sModeString)
						, default0Save1Load2
#if COM16_C806_EMT
						, &uiSigNum
#endif
					);

					singleCbfTmpLuma = pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrDepth);

					//----- determine rate and r-d cost -----
#if COM16_C806_EMT
					if ((modeId == 1 && singleCbfTmpLuma == 0) || (modeId == 0 && ucTrIdx && ucTrIdx != DCT2_EMT && uiSigNum <= g_iEmtSigNumThr))
#else
					if (modeId == 1 && singleCbfTmpLuma == 0)
#endif
					{
						//In order not to code TS flag when cbf is zero, the case for TS with cbf being zero is forbidden.
						singleCostTmp = MAX_DOUBLE;
					}
					else
					{
						UInt uiSingleBits = xGetIntraBitsQT(rTu, true, false, false);
						singleCostTmp = m_pcRdCost->calcRdCost(uiSingleBits, singleDistTmpLuma);
					}
					if (singleCostTmp < dSingleCost)
					{
						DEBUG_STRING_SWAP(sDebug, sModeString)
							dSingleCost = singleCostTmp;
						uiSingleDistLuma = singleDistTmpLuma;
						uiSingleCbfLuma = singleCbfTmpLuma;
#if COM16_C806_EMT
						bestTrIdx = ucTrIdx;
#endif
						bestModeId[COMPONENT_Y] = modeId;
						if (bestModeId[COMPONENT_Y] == firstCheckId)
						{
							xStoreIntraResultQT(COMPONENT_Y, rTu);
#if JVET_C0024_QTBT
							m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
							m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
						}

						if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
						{
							const Int xOffset = rTu.getRect(COMPONENT_Y).x0;
							const Int yOffset = rTu.getRect(COMPONENT_Y).y0;
							for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
							{
								if (bMaintainResidual[storedResidualIndex])
								{
									xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
										m_resiSingleBuffer[uiWIdx][uiHIdx][storedResidualIndex],
#else
										m_resiSingleBuffer[uiLog2TrSize][storedResidualIndex],
#endif
#else
										resiLumaSingle[storedResidualIndex],
#endif
										rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
								}
							}
						}
					}
					if (modeId == firstCheckId)
					{
#if JVET_C0024_QTBT
						m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
						m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif
					}
#if COM16_C806_EMT
				}
#endif
			}
#if JVET_C0024_QTBT
			pcCU->setTransformSkipSubParts(bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, 0);
#else
			pcCU->setTransformSkipSubParts(bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
#if COM16_C806_EMT
			pcCU->setEmtTuIdxSubParts(bestTrIdx, uiAbsPartIdx, uiFullDepth);
#endif

			if (bestModeId[COMPONENT_Y] == firstCheckId)
			{
				xLoadIntraResultQT(COMPONENT_Y, rTu);
				pcCU->setCbfSubParts(uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, rTu.GetTransformDepthTotalAdj(COMPONENT_Y));

#if JVET_C0024_QTBT
				m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
			}
		}
		else
		{
#if COM16_C806_EMT
#if HHI_RQT_INTRA_SPEEDUP
			UChar numTrIdxCands = (!bCheckInitTrDepth) ? nNumTrCands : 1;
#else
			UChar numTrIdxCands = nNumTrCands;
#endif

			for (UChar ucTrIdx = 0; ucTrIdx < numTrIdxCands; ucTrIdx++)
			{
				UInt   singleDistYTmp = 0;
				UInt   singleCbfYTmp = 0;
				Double singleCostTmp = 0;
				Bool   bSaveEmtResults = ucTrIdx<(numTrIdxCands - 1);

				// Skip checking other transform candidates if zero CBF is encountered
				if (ucTrIdx && !uiSingleCbfLuma && bAllIntra && m_pcEncCfg->getUseFastIntraEMT())
				{
					continue;
				}
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
				if (this->getSaveLoadTag(pcCU->getZorderIdxInCtu(), pcCU->getCURouteIdx(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && ucTrIdx && ucTrIdx != this->getSaveLoadEmtIdx(pcCU->getZorderIdxInCtu(), pcCU->getCURouteIdx(), uiWIdx, uiHIdx))
#else
				if (this->getSaveLoadTag(pcCU->getZorderIdxInCtu(), uiWIdx, uiHIdx) == LOAD_ENC_INFO && ucTrIdx && ucTrIdx != this->getSaveLoadEmtIdx(uiWIdx, uiHIdx))
#endif
				{
					continue;
				}
#endif
#endif

				//----- store original entropy coding status -----
#if VCEG_AZ08_INTRA_KLT
#if COM16_C806_EMT
				if (bCheckSplit || bSaveEmtResults || checkTM)
#else
				if (bCheckSplit)
#endif
#else
#if COM16_C806_EMT
				if (bCheckSplit || bSaveEmtResults)
#else
				if (bCheckSplit)
#endif
#endif
				{
#if JVET_C0024_QTBT
					m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
					m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif
				}

#if COM16_C806_EMT
				Int default0Save1Load2 = numTrIdxCands>1 ? (bSaveEmtResults ? 1 : 2) : 0;
#endif

				//----- code luma/chroma block with given intra prediction mode and store Cbf-----

#if COM16_C806_EMT
#if HHI_RQT_INTRA_SPEEDUP
				pcCU->setEmtTuIdxSubParts(bCheckInitTrDepth ? ucSavedEmtTrIdx : ucTrIdx, uiAbsPartIdx, uiFullDepth);
#else
				pcCU->setEmtTuIdxSubParts(trIdx, uiAbsPartIdx, uiFullDepth);
#endif
#if JVET_C0024_QTBT
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, 0);
#else
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
#if VCEG_AZ08_INTRA_KLT
				pcCU->setKLTFlagSubParts(0, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#endif
				xIntraCodingTUBlock(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
					m_resiSingleBuffer[uiWIdx][uiHIdx],
#else
					m_resiSingleBuffer[uiLog2TrSize],
#endif
#else
					resiLumaSingle,
#endif
					false, singleDistYTmp, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sDebug), default0Save1Load2, &uiSigNum);
#else
				dSingleCost = 0.0;
#if JVET_C0024_QTBT
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, 0);
#else
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
				xIntraCodingTUBlock(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
					m_resiSingleBuffer[uiWIdx][uiHIdx],
#else
					m_resiSingleBuffer[uiLog2TrSize],
#endif
#else
					resiLumaSingle,
#endif
					false, uiSingleDistLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sDebug));
#endif

#if VCEG_AZ08_INTRA_KLT
#if COM16_C806_EMT
				if (bCheckSplit || bSaveEmtResults || checkTM)
#else
				if (bCheckSplit)
#endif
#else
#if COM16_C806_EMT
				if (bCheckSplit || bSaveEmtResults)
#else
				if (bCheckSplit)
#endif
#endif
				{
#if COM16_C806_EMT
					singleCbfYTmp = pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrDepth);
#else
					uiSingleCbfLuma = pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrDepth);
#endif
				}
				//----- determine rate and r-d cost -----
#if COM16_C806_EMT
				if (ucTrIdx && ucTrIdx != DCT2_EMT && uiSigNum <= g_iEmtSigNumThr)
				{
					singleCostTmp = MAX_DOUBLE;
				}
				else
				{
#endif
					UInt uiSingleBits = xGetIntraBitsQT(rTu, true, false, false);
#if JVET_C0024_QTBT
					if (m_pcEncCfg->getRDpenalty() && (uiWIdx + uiHIdx + (1 << MIN_CU_LOG2) >= 10) && !isIntraSlice)
#else
					if (m_pcEncCfg->getRDpenalty() && (uiLog2TrSize == 5) && !isIntraSlice)
#endif
					{
						uiSingleBits = uiSingleBits * 4;
					}

#if COM16_C806_EMT
					singleCostTmp = m_pcRdCost->calcRdCost(uiSingleBits, singleDistYTmp);
				}
#else
					dSingleCost = m_pcRdCost->calcRdCost(uiSingleBits, uiSingleDistLuma);
#endif

#if COM16_C806_EMT
					if (singleCostTmp < dSingleCost)
					{
#endif
						if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
						{
							const Int xOffset = rTu.getRect(COMPONENT_Y).x0;
							const Int yOffset = rTu.getRect(COMPONENT_Y).y0;
							for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
							{
								if (bMaintainResidual[storedResidualIndex])
								{
									xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
#if JVET_C0024_QTBT
										m_resiSingleBuffer[uiWIdx][uiHIdx][storedResidualIndex],
#else
										m_resiSingleBuffer[uiLog2TrSize][storedResidualIndex],
#endif
#else
										resiLumaSingle[storedResidualIndex],
#endif
										rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
								}
							}
						}
#if COM16_C806_EMT
						dSingleCost = singleCostTmp;
						uiSingleDistLuma = singleDistYTmp;
						uiSingleCbfLuma = singleCbfYTmp;
						bestTrIdx = ucTrIdx;

						if (bSaveEmtResults && (uiSingleCbfLuma || !bAllIntra || !m_pcEncCfg->getUseFastIntraEMT()))
						{
							xStoreIntraResultQT(COMPONENT_Y, rTu);
#if JVET_C0024_QTBT
							m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
							m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
						}
					}

					if (bSaveEmtResults && (uiSingleCbfLuma || !bAllIntra || !m_pcEncCfg->getUseFastIntraEMT()))
					{
#if JVET_C0024_QTBT
						m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
						m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif
					}
			}

				pcCU->setEmtTuIdxSubParts(bestTrIdx, uiAbsPartIdx, uiFullDepth);

				if (bestTrIdx < (numTrIdxCands - 1) && (uiSingleCbfLuma || !bAllIntra || !m_pcEncCfg->getUseFastIntraEMT()))
				{
					xLoadIntraResultQT(COMPONENT_Y, rTu);
					pcCU->setCbfSubParts(uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#if JVET_C0024_QTBT
					m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
				}
#endif
		}
#if VCEG_AZ08_INTRA_KLT
			if (checkTM && uiSingleCbfLuma)
			{
				//backup the former results
				Double  dSingleCostBackUp = dSingleCost;
				UInt    uiSingleCbfYBackUp = uiSingleCbfLuma;
				UInt    uiSingleDistYBackUp = uiSingleDistLuma;

				uiSingleDistLuma = 0;
				dSingleCost = 0.0;
				uiSingleCbfLuma = 0;
#if VCEG_AZ08_INTRA_KLT
				pcCU->setKLTFlagSubParts(1, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#endif
				xStoreIntraResultQT(COMPONENT_Y, rTu);
#if JVET_C0024_QTBT
				m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
				//----- store original entropy coding status -----
				m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
				//----- store original entropy coding status -----
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif
				//----- code luma block with given intra prediction mode and store Cbf-----
				Bool bSuccessful;

#if COM16_C806_EMT
				pcCU->setEmtTuIdxSubParts(0, uiAbsPartIdx, uiFullDepth);
#if JVET_C0024_QTBT
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, 0);
#else
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
				bSuccessful = xIntraCodingTUBlockTM(pcOrgYuv, pcPredYuv, pcResiYuv,
					uiSingleDistLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sDebug), &uiSigNum, TMPRED0_TMPREDKLT1_ORI2);
#else
				dSingleCost = 0.0;
#if JVET_C0024_QTBT
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, 0);
#else
				pcCU->setTransformSkipSubParts(0, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#endif
#if JVET_C0024_QTBT
				bSuccessful =
#endif
					xIntraCodingTUBlockTM(pcOrgYuv, pcPredYuv, pcResiYuv,
						uiSingleDistLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sDebug)
#if JVET_C0024_QTBT
						, TMPRED0_TMPREDKLT1_ORI2
#else
#endif
					);
#endif
				uiSingleCbfLuma = pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrDepth);

				if (bSuccessful == true && uiSingleCbfLuma) //assume only cbf being nonzero, can use TM mode
				{
					//----- determine rate and r-d cost ----- 
					UInt uiSingleBits = xGetIntraBitsQT(rTu, true, false, false);
#if JVET_C0024_QTBT
					if (m_pcEncCfg->getRDpenalty() && (uiWIdx + uiHIdx + (MIN_CU_LOG2 << 1) >= 10) && !isIntraSlice)
#else
					if (m_pcEncCfg->getRDpenalty() && (uiLog2TrSize == 5) && !isIntraSlice)
#endif
					{
						uiSingleBits = uiSingleBits * 4;
					}
					dSingleCost = m_pcRdCost->calcRdCost(uiSingleBits, uiSingleDistLuma);
				}
				else
				{
					dSingleCost = MAX_DOUBLE;
				}

				if (dSingleCostBackUp <= dSingleCost)
				{
					dSingleCost = dSingleCostBackUp;
					uiSingleCbfLuma = uiSingleCbfYBackUp;
					uiSingleDistLuma = uiSingleDistYBackUp;
					pcCU->setKLTFlagSubParts(0, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#if COM16_C806_EMT
					pcCU->setEmtTuIdxSubParts(bestTrIdx, uiAbsPartIdx, uiFullDepth);
#endif
#if JVET_C0024_QTBT
					pcCU->setTransformSkipSubParts(bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, 0); //1.30.2016  uiFullDepth);
#else
					pcCU->setTransformSkipSubParts(bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan); //1.30.2016  uiFullDepth);
#endif
					xLoadIntraResultQT(COMPONENT_Y, rTu);
					pcCU->setCbfSubParts(uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#if JVET_C0024_QTBT
					m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
				}
				else
				{
#if !JVET_C0024_QTBT
					bestTMKLT = 1;
#endif
					bestModeId[COMPONENT_Y] = 0;
#if JVET_C0024_QTBT && COM16_C806_EMT
					bestTrIdx = 0;
#endif
#if COM16_C806_EMT
					pcCU->setEmtTuIdxSubParts(0, uiAbsPartIdx, uiFullDepth);
#endif
				}
			}
#endif
	}

#if !JVET_C0024_QTBT
		if (bCheckSplit)
		{
			//----- store full entropy coding status, load original entropy coding status -----
			if (bCheckFull)
			{
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_TEST]);
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
			}
			else
			{
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
			}
			//----- code splitted block -----
			Double     dSplitCost = 0.0;
			Distortion uiSplitDistLuma = 0;
			UInt       uiSplitCbfLuma = 0;
#if COM16_C806_EMT
			Bool    bSplitSelected = true;
#endif

			TComTURecurse tuRecurseChild(rTu, false);
			DEBUG_STRING_NEW(sSplit)
				do
				{
					DEBUG_STRING_NEW(sChild)
#if HHI_RQT_INTRA_SPEEDUP
						xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
							m_resiSplitBuffer[uiLog2TrSize],
#else
							resiLumaSplit,
#endif
							uiSplitDistLuma, bCheckFirst, dSplitCost, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild));
#else
						xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
							m_resiSplitBuffer[uiLog2TrSize],
#else
							resiLumaSplit,
#endif
							uiSplitDistLuma, dSplitCost, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild));
#endif
					DEBUG_STRING_APPEND(sSplit, sChild)
						uiSplitCbfLuma |= pcCU->getCbf(tuRecurseChild.GetAbsPartIdxTU(), COMPONENT_Y, tuRecurseChild.GetTransformDepthRel());
#if COM16_C806_EMT
					if (dSplitCost>dSingleCost && m_pcEncCfg->getUseFastIntraEMT())
					{
						bSplitSelected = false;
						break;
					}
#endif
				} while (tuRecurseChild.nextSection(rTu));

#if COM16_C806_EMT
				if (bSplitSelected)
				{
#endif
					UInt    uiPartsDiv = rTu.GetAbsPartIdxNumParts();
					{
						if (uiSplitCbfLuma)
						{
							const UInt flag = 1 << uiTrDepth;
							UChar *pBase = pcCU->getCbf(COMPONENT_Y);
							for (UInt uiOffs = 0; uiOffs < uiPartsDiv; uiOffs++)
							{
								pBase[uiAbsPartIdx + uiOffs] |= flag;
							}
						}
					}
					//----- restore context states -----
					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);

					//----- determine rate and r-d cost -----
					UInt uiSplitBits = xGetIntraBitsQT(rTu, true, false, false);
					dSplitCost = m_pcRdCost->calcRdCost(uiSplitBits, uiSplitDistLuma);
#if COM16_C806_EMT
				}
#endif

				//===== compare and set best =====
				if (dSplitCost < dSingleCost)
				{
					//--- update cost ---
					DEBUG_STRING_SWAP(sSplit, sDebug)
						ruiDistY += uiSplitDistLuma;
					dRDCost += dSplitCost;

					if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
					{
						const Int xOffset = rTu.getRect(COMPONENT_Y).x0;
						const Int yOffset = rTu.getRect(COMPONENT_Y).y0;
						for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
						{
							if (bMaintainResidual[storedResidualIndex])
							{
								xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
									m_resiSplitBuffer[uiLog2TrSize][storedResidualIndex],
#else
									resiLumaSplit[storedResidualIndex],
#endif
									rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
							}
						}
					}

					return;
				}

				//----- set entropy coding status -----
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_TEST]);

				//--- set transform index and Cbf values ---
				pcCU->setTrIdxSubParts(uiTrDepth, uiAbsPartIdx, uiFullDepth);
				const TComRectangle &tuRect = rTu.getRect(COMPONENT_Y);
				pcCU->setCbfSubParts(uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
				pcCU->setTransformSkipSubParts(bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan);
#if COM16_C806_EMT
				pcCU->setEmtTuIdxSubParts(bestTrIdx, uiAbsPartIdx, uiFullDepth);
#endif
#if VCEG_AZ08_INTRA_KLT
				pcCU->setKLTFlagSubParts(bestTMKLT, COMPONENT_Y, uiAbsPartIdx, uiFullDepth);
#endif

				//--- set reconstruction for next intra prediction blocks ---
				const UInt  uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
				const UInt  uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
				const UInt  uiWidth = tuRect.width;
				const UInt  uiHeight = tuRect.height;
				Pel*  piSrc = m_pcQTTempTComYuv[uiQTLayer].getAddr(COMPONENT_Y, uiAbsPartIdx);
				UInt  uiSrcStride = m_pcQTTempTComYuv[uiQTLayer].getStride(COMPONENT_Y);
				Pel*  piDes = pcCU->getPic()->getPicYuvRec()->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiZOrder);
				UInt  uiDesStride = pcCU->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);

				for (UInt uiY = 0; uiY < uiHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride)
				{
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						piDes[uiX] = piSrc[uiX];
					}
				}
		}
#endif

		ruiDistY += uiSingleDistLuma;
		dRDCost += dSingleCost;
}


	Void
		TEncSearch::xSetIntraResultLumaQT(TComYuv* pcRecoYuv, TComTU &rTu)
	{
		TComDataCU *pcCU = rTu.getCU();
		const UInt uiTrDepth = rTu.GetTransformDepthRel();
		const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#if JVET_C0024_QTBT
		assert(uiTrDepth == 0 && uiAbsPartIdx == 0);
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
		if (uiTrMode == uiTrDepth)
		{
			UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
			UInt uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif

			//===== copy transform coefficients =====

			const TComRectangle &tuRect = rTu.getRect(COMPONENT_Y);
			const UInt coeffOffset = rTu.getCoefficientOffset(COMPONENT_Y);
			const UInt numCoeffInBlock = tuRect.width * tuRect.height;

			if (numCoeffInBlock != 0)
			{
#if JVET_C0024_QTBT
				assert(coeffOffset == 0);
				const TCoeff* srcCoeff = m_pppcQTTempCoeff[COMPONENT_Y][uiWIdx][uiHIdx] + coeffOffset;
#else
				const TCoeff* srcCoeff = m_ppcQTTempCoeff[COMPONENT_Y][uiQTLayer] + coeffOffset;
#endif
				TCoeff* destCoeff = pcCU->getCoeff(COMPONENT_Y) + coeffOffset;
				::memcpy(destCoeff, srcCoeff, sizeof(TCoeff)*numCoeffInBlock);
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
				const TCoeff* srcArlCoeff = m_pppcQTTempArlCoeff[COMPONENT_Y][uiWIdx][uiHIdx] + coeffOffset;
#else
				const TCoeff* srcArlCoeff = m_ppcQTTempArlCoeff[COMPONENT_Y][uiQTLayer] + coeffOffset;
#endif
				TCoeff* destArlCoeff = pcCU->getArlCoeff(COMPONENT_Y) + coeffOffset;
				::memcpy(destArlCoeff, srcArlCoeff, sizeof(TCoeff) * numCoeffInBlock);
#endif
#if JVET_C0024_QTBT
				m_ppcQTTempTComYuv[uiWIdx][uiHIdx].copyPartToPartComponent(COMPONENT_Y, pcRecoYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
#else
				m_pcQTTempTComYuv[uiQTLayer].copyPartToPartComponent(COMPONENT_Y, pcRecoYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
#endif
			}

#if !JVET_C0024_QTBT
		}
		else
		{
			TComTURecurse tuRecurseChild(rTu, false);
			do
			{
				xSetIntraResultLumaQT(pcRecoYuv, tuRecurseChild);
			} while (tuRecurseChild.nextSection(rTu));
		}
#endif
	}


#if COM16_C806_LMCHROMA
	Void TEncSearch::xStoreIntraResultQT(const ComponentID compID, TComTU &rTu, TComYuv* pcResiYuv)
#else
	Void TEncSearch::xStoreIntraResultQT(const ComponentID compID, TComTU &rTu)
#endif
	{
		TComDataCU *pcCU = rTu.getCU();
		const UInt uiTrDepth = rTu.GetTransformDepthRel();
		const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#if JVET_C0024_QTBT
		assert(uiTrDepth == 0 && uiAbsPartIdx == 0);
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		const UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
		if (compID == COMPONENT_Y || uiTrMode == uiTrDepth)
		{
			assert(uiTrMode == uiTrDepth);
			const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
			const UInt uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif

			if (rTu.ProcessComponentSection(compID))
			{
				const TComRectangle &tuRect = rTu.getRect(compID);

				//===== copy transform coefficients =====
				const UInt uiNumCoeff = tuRect.width * tuRect.height;
#if JVET_C0024_QTBT
				TCoeff* pcCoeffSrc = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
				TCoeff* pcCoeffSrc = m_ppcQTTempCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
				TCoeff* pcCoeffDst = m_pcQTTempTUCoeff[compID];

				::memcpy(pcCoeffDst, pcCoeffSrc, sizeof(TCoeff) * uiNumCoeff);
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
				TCoeff* pcArlCoeffSrc = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
				TCoeff* pcArlCoeffSrc = m_ppcQTTempArlCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
				TCoeff* pcArlCoeffDst = m_ppcQTTempTUArlCoeff[compID];
				::memcpy(pcArlCoeffDst, pcArlCoeffSrc, sizeof(TCoeff) * uiNumCoeff);
#endif
				//===== copy reconstruction =====
#if JVET_C0024_QTBT
				m_ppcQTTempTComYuv[uiWIdx][uiHIdx].copyPartToPartComponent(compID, &m_pcQTTempTransformSkipTComYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
#else
				m_pcQTTempTComYuv[uiQTLayer].copyPartToPartComponent(compID, &m_pcQTTempTransformSkipTComYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
#endif
			}
#if !JVET_C0024_QTBT
		}
#endif

#if COM16_C806_LMCHROMA
		if (compID == COMPONENT_Cb && pcCU->getSlice()->getSPS()->getUseLMChroma())
		{
			const TComRectangle &tuRect = rTu.getRect(compID);
			pcResiYuv->copyPartToPartComponent(compID, &m_pcQTTempResiTComYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
		}
#endif
	}


#if COM16_C806_LMCHROMA
	Void TEncSearch::xLoadIntraResultQT(const ComponentID compID, TComTU &rTu, TComYuv* pcResiYuv)
#else
	Void TEncSearch::xLoadIntraResultQT(const ComponentID compID, TComTU &rTu)
#endif
	{
		TComDataCU *pcCU = rTu.getCU();
		const UInt uiTrDepth = rTu.GetTransformDepthRel();
		const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#if JVET_C0024_QTBT
		assert(uiTrDepth == 0 && uiAbsPartIdx == 0);
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		const UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
		if (compID == COMPONENT_Y || uiTrMode == uiTrDepth)
		{
			assert(uiTrMode == uiTrDepth);
			const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
			const UInt uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif
			const UInt uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;

			if (rTu.ProcessComponentSection(compID))
			{
				const TComRectangle &tuRect = rTu.getRect(compID);

				//===== copy transform coefficients =====
				const UInt uiNumCoeff = tuRect.width * tuRect.height;
#if JVET_C0024_QTBT
				TCoeff* pcCoeffDst = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
				TCoeff* pcCoeffDst = m_ppcQTTempCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
				TCoeff* pcCoeffSrc = m_pcQTTempTUCoeff[compID];

				::memcpy(pcCoeffDst, pcCoeffSrc, sizeof(TCoeff) * uiNumCoeff);
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
				TCoeff* pcArlCoeffDst = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#else
				TCoeff* pcArlCoeffDst = m_ppcQTTempArlCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
#endif
				TCoeff* pcArlCoeffSrc = m_ppcQTTempTUArlCoeff[compID];
				::memcpy(pcArlCoeffDst, pcArlCoeffSrc, sizeof(TCoeff) * uiNumCoeff);
#endif
				//===== copy reconstruction =====
#if JVET_C0024_QTBT
				m_pcQTTempTransformSkipTComYuv.copyPartToPartComponent(compID, &m_ppcQTTempTComYuv[uiWIdx][uiHIdx], uiAbsPartIdx, tuRect.width, tuRect.height);
#else
				m_pcQTTempTransformSkipTComYuv.copyPartToPartComponent(compID, &m_pcQTTempTComYuv[uiQTLayer], uiAbsPartIdx, tuRect.width, tuRect.height);
#endif

				Pel*    piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
				UInt    uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
#if JVET_C0024_QTBT
				Pel*    piRecQt = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
				UInt    uiRecQtStride = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
				Pel*    piRecQt = m_pcQTTempTComYuv[uiQTLayer].getAddr(compID, uiAbsPartIdx);
				UInt    uiRecQtStride = m_pcQTTempTComYuv[uiQTLayer].getStride(compID);
#endif
				UInt    uiWidth = tuRect.width;
				UInt    uiHeight = tuRect.height;
				Pel* pRecQt = piRecQt;
				Pel* pRecIPred = piRecIPred;
				for (UInt uiY = 0; uiY < uiHeight; uiY++)
				{
					for (UInt uiX = 0; uiX < uiWidth; uiX++)
					{
						pRecIPred[uiX] = pRecQt[uiX];
					}
					pRecQt += uiRecQtStride;
					pRecIPred += uiRecIPredStride;
				}
			}
#if !JVET_C0024_QTBT
		}
#endif

#if COM16_C806_LMCHROMA
		if (compID == COMPONENT_Cb && pcCU->getSlice()->getSPS()->getUseLMChroma())
		{
			const TComRectangle &tuRect = rTu.getRect(compID);
			m_pcQTTempResiTComYuv.copyPartToPartComponent(compID, pcResiYuv, uiAbsPartIdx, tuRect.width, tuRect.height);
		}
#endif
	}

	Void
		TEncSearch::xStoreCrossComponentPredictionResult(Pel    *pResiDst,
			const Pel    *pResiSrc,
			TComTU &rTu,
			const Int     xOffset,
			const Int     yOffset,
			const Int     strideDst,
			const Int     strideSrc)
	{
		const Pel *pSrc = pResiSrc + yOffset * strideSrc + xOffset;
		Pel *pDst = pResiDst + yOffset * strideDst + xOffset;

		twoDimBufCopier
			(
				(UChar *)pDst,
				strideDst * sizeof(Pel),
				(UChar *)pSrc,
				strideSrc * sizeof(Pel),
				rTu.getRect(COMPONENT_Y).width * sizeof(Pel),
				rTu.getRect(COMPONENT_Y).height
				);
		//for (Int y = 0; y < rTu.getRect(COMPONENT_Y).height; y++)
		//{
		//	::memcpy(pDst, pSrc, sizeof(Pel) * rTu.getRect(COMPONENT_Y).width);
		//	pDst += strideDst;
		//	pSrc += strideSrc;
		//}
	}

	Char
		TEncSearch::xCalcCrossComponentPredictionAlpha(TComTU &rTu,
			const ComponentID compID,
			const Pel*        piResiL,
			const Pel*        piResiC,
			const Int         width,
			const Int         height,
			const Int         strideL,
			const Int         strideC)
	{
		const Pel *pResiL = piResiL;
		const Pel *pResiC = piResiC;

		TComDataCU *pCU = rTu.getCU();
		const Int  absPartIdx = rTu.GetAbsPartIdxTU(compID);
		const Int diffBitDepth = pCU->getSlice()->getSPS()->getDifferentialLumaChromaBitDepth();

		Char alpha = 0;
		Int SSxy = 0;
		Int SSxx = 0;

		for (UInt uiY = 0; uiY < height; uiY++)
		{
			for (UInt uiX = 0; uiX < width; uiX++)
			{
				const Pel scaledResiL = rightShift(pResiL[uiX], diffBitDepth);
				SSxy += (scaledResiL * pResiC[uiX]);
				SSxx += (scaledResiL * scaledResiL);
			}

			pResiL += strideL;
			pResiC += strideC;
		}

		if (SSxx != 0)
		{
			Double dAlpha = SSxy / Double(SSxx);
			alpha = Char(Clip3<Int>(-16, 16, (Int)(dAlpha * 16)));

			static const Char alphaQuant[17] = { 0, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8 };

			alpha = (alpha < 0) ? -alphaQuant[Int(-alpha)] : alphaQuant[Int(alpha)];
		}
		pCU->setCrossComponentPredictionAlphaPartRange(alpha, compID, absPartIdx, rTu.GetAbsPartIdxNumParts(compID));

		return alpha;
	}

	Void
		TEncSearch::xRecurIntraChromaCodingQT(TComYuv*    pcOrgYuv,
			TComYuv*    pcPredYuv,
			TComYuv*    pcResiYuv,
#if COM16_C806_LARGE_CTU
			Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES],
#else
			Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
#endif
			Distortion& ruiDist,
			TComTU&     rTu
			DEBUG_STRING_FN_DECLARE(sDebug))
	{
		TComDataCU         *pcCU = rTu.getCU();
		const UInt          uiTrDepth = rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
		assert(rTu.GetAbsPartIdxTU() == 0 && uiTrDepth == 0);
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		const UInt          uiAbsPartIdx = rTu.GetAbsPartIdxTU();
#endif
		const ChromaFormat  format = rTu.GetChromaFormat();
#if !JVET_C0024_QTBT
		UInt                uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
#endif
		const UInt          numberValidComponents = getNumberValidComponents(format);

#if JVET_C0024_QTBT
		assert(uiTrDepth == 0);
#else
		if (uiTrMode == uiTrDepth)
		{
#endif
			if (!rTu.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
			{
				return;
			}

#if !JVET_C0024_QTBT
			const UInt uiFullDepth = rTu.GetTransformDepthTotal();
#endif

			Bool checkTransformSkip = pcCU->getSlice()->getPPS()->getUseTransformSkip();
#if JVET_C0024_QTBT
			checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(pcCU->getSlice()->isIntra(), rTu.getRect(COMPONENT_Cb), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
#else
			checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Cb), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
#endif

			if (m_pcEncCfg->getUseTransformSkipFast())
			{
#if JVET_C0024_QTBT
				checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(pcCU->getSlice()->isIntra(), rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
#else
				checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());

				if (checkTransformSkip)
				{
					Int nbLumaSkip = 0;
					const UInt maxAbsPartIdxSub = uiAbsPartIdx + (rTu.ProcessingAllQuadrants(COMPONENT_Cb) ? 1 : 4);
					for (UInt absPartIdxSub = uiAbsPartIdx; absPartIdxSub < maxAbsPartIdxSub; absPartIdxSub++)
					{
						nbLumaSkip += pcCU->getTransformSkip(absPartIdxSub, COMPONENT_Y);
					}
					checkTransformSkip &= (nbLumaSkip > 0);
				}
#endif
			}


			for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
			{
				const ComponentID compID = ComponentID(ch);
				DEBUG_STRING_NEW(sDebugBestMode)

					//use RDO to decide whether Cr/Cb takes TS
#if JVET_C0024_QTBT
					m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
				const Bool splitIntoSubTUs = false;
#else
					m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);

				const Bool splitIntoSubTUs = rTu.getRect(compID).width != rTu.getRect(compID).height;
#endif

				TComTURecurse TUIterator(rTu, false, (splitIntoSubTUs ? TComTU::VERTICAL_SPLIT : TComTU::DONT_SPLIT), true, compID);

#if !JVET_C0024_QTBT
				const UInt partIdxesPerSubTU = TUIterator.GetAbsPartIdxNumParts(compID);
#endif

				do
				{
					const UInt subTUAbsPartIdx = TUIterator.GetAbsPartIdxTU(compID);

					Double     dSingleCost = MAX_DOUBLE;
					Int        bestModeId = 0;
					Distortion singleDistC = 0;
					UInt       singleCbfC = 0;
					Distortion singleDistCTmp = 0;
					Double     singleCostTmp = 0;
					UInt       singleCbfCTmp = 0;
					Char       bestCrossCPredictionAlpha = 0;
					Int        bestTransformSkipMode = 0;
#if JVET_E0062_MULTI_DMS
#if JVET_E0077_ENHANCED_LM
					const Bool checkCrossComponentPrediction = (!IsLMMode(pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, subTUAbsPartIdx)))
						&& pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
						&& (pcCU->getCbf(subTUAbsPartIdx, COMPONENT_Y, uiTrDepth) != 0);
#else
					const Bool checkCrossComponentPrediction = (pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, subTUAbsPartIdx) != LM_CHROMA_IDX)
						&& pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
						&& (pcCU->getCbf(subTUAbsPartIdx, COMPONENT_Y, uiTrDepth) != 0);
#endif
#else
					const Bool checkCrossComponentPrediction = (pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, subTUAbsPartIdx) == DM_CHROMA_IDX)
						&& pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
						&& (pcCU->getCbf(subTUAbsPartIdx, COMPONENT_Y, uiTrDepth) != 0);
#endif
					const Int  crossCPredictionModesToTest = checkCrossComponentPrediction ? 2 : 1;
					const Int  transformSkipModesToTest = checkTransformSkip ? 2 : 1;
					const Int  totalModesToTest = crossCPredictionModesToTest * transformSkipModesToTest;
					Int  currModeId = 0;
					Int  default0Save1Load2 = 0;

					for (Int transformSkipModeId = 0; transformSkipModeId < transformSkipModesToTest; transformSkipModeId++)
					{
						for (Int crossCPredictionModeId = 0; crossCPredictionModeId < crossCPredictionModesToTest; crossCPredictionModeId++)
						{
#if JVET_C0024_QTBT
							pcCU->setCrossComponentPredictionAlphaPartRange(0, compID, subTUAbsPartIdx, 0);
							pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, 0);
#else
							pcCU->setCrossComponentPredictionAlphaPartRange(0, compID, subTUAbsPartIdx, partIdxesPerSubTU);
							DEBUG_STRING_NEW(sDebugMode)
								pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
							currModeId++;

							const Bool isOneMode = (totalModesToTest == 1);
							const Bool isLastMode = (currModeId == totalModesToTest); // currModeId is indexed from 1

							if (isOneMode)
							{
								default0Save1Load2 = 0;
							}
							else if (!isOneMode && (transformSkipModeId == 0) && (crossCPredictionModeId == 0))
							{
								default0Save1Load2 = 1; //save prediction on first mode
							}
							else
							{
								default0Save1Load2 = 2; //load it on subsequent modes
							}

							singleDistCTmp = 0;

							xIntraCodingTUBlock(pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, (crossCPredictionModeId != 0), singleDistCTmp, compID, TUIterator DEBUG_STRING_PASS_INTO(sDebugMode), default0Save1Load2);

							singleCbfCTmp = pcCU->getCbf(subTUAbsPartIdx, compID, uiTrDepth);

							if (((crossCPredictionModeId == 1) && (pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) == 0))
								|| ((transformSkipModeId == 1) && (singleCbfCTmp == 0))) //In order not to code TS flag when cbf is zero, the case for TS with cbf being zero is forbidden.
							{
								singleCostTmp = MAX_DOUBLE;
							}
							else if (!isOneMode)
							{
								UInt bitsTmp = xGetIntraBitsQTChroma(TUIterator, compID, false);
								singleCostTmp = m_pcRdCost->calcRdCost(bitsTmp, singleDistCTmp);
							}

							if (singleCostTmp < dSingleCost)
							{
								DEBUG_STRING_SWAP(sDebugBestMode, sDebugMode)
									dSingleCost = singleCostTmp;
								singleDistC = singleDistCTmp;
								bestCrossCPredictionAlpha = (crossCPredictionModeId != 0) ? pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) : 0;
								bestTransformSkipMode = transformSkipModeId;
								bestModeId = currModeId;
								singleCbfC = singleCbfCTmp;

								if (!isOneMode && !isLastMode)
								{
#if COM16_C806_LMCHROMA
									xStoreIntraResultQT(compID, TUIterator, pcResiYuv);
#else
									xStoreIntraResultQT(compID, TUIterator);
#endif
#if JVET_C0024_QTBT
									m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
									m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
								}
							}

							if (!isOneMode && !isLastMode)
							{
#if JVET_C0024_QTBT
								m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
								m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT]);
#endif
							}
						}
					}

					if (bestModeId < totalModesToTest)
					{
#if COM16_C806_LMCHROMA
						xLoadIntraResultQT(compID, TUIterator, pcResiYuv);
#else
						xLoadIntraResultQT(compID, TUIterator);
#endif
#if JVET_C0024_QTBT
						pcCU->setCbfPartRange(singleCbfC, compID, subTUAbsPartIdx, 0);
						m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
						pcCU->setCbfPartRange(singleCbfC << uiTrDepth, compID, subTUAbsPartIdx, partIdxesPerSubTU);

						m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiFullDepth][CI_TEMP_BEST]);
#endif
					}

					DEBUG_STRING_APPEND(sDebug, sDebugBestMode)
#if JVET_C0024_QTBT
						pcCU->setTransformSkipPartRange(bestTransformSkipMode, compID, subTUAbsPartIdx, 0);
					pcCU->setCrossComponentPredictionAlphaPartRange(bestCrossCPredictionAlpha, compID, subTUAbsPartIdx, 0);
#else
						pcCU->setTransformSkipPartRange(bestTransformSkipMode, compID, subTUAbsPartIdx, partIdxesPerSubTU);
					pcCU->setCrossComponentPredictionAlphaPartRange(bestCrossCPredictionAlpha, compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
					ruiDist += singleDistC;
				} while (TUIterator.nextSection(rTu));

#if !JVET_C0024_QTBT
				if (splitIntoSubTUs)
				{
					offsetSubTUCBFs(rTu, compID);
				}
#endif
			}
#if !JVET_C0024_QTBT
		}
		else
		{
			UInt    uiSplitCbf[MAX_NUM_COMPONENT] = { 0,0,0 };

			TComTURecurse tuRecurseChild(rTu, false);
			const UInt uiTrDepthChild = tuRecurseChild.GetTransformDepthRel();
			do
			{
				DEBUG_STRING_NEW(sChild)

					xRecurIntraChromaCodingQT(pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, ruiDist, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild));

				DEBUG_STRING_APPEND(sDebug, sChild)
					const UInt uiAbsPartIdxSub = tuRecurseChild.GetAbsPartIdxTU();

				for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
				{
					uiSplitCbf[ch] |= pcCU->getCbf(uiAbsPartIdxSub, ComponentID(ch), uiTrDepthChild);
				}
			} while (tuRecurseChild.nextSection(rTu));


			UInt uiPartsDiv = rTu.GetAbsPartIdxNumParts();
			for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
			{
				if (uiSplitCbf[ch])
				{
					const UInt flag = 1 << uiTrDepth;
					ComponentID compID = ComponentID(ch);
					UChar *pBase = pcCU->getCbf(compID);
					for (UInt uiOffs = 0; uiOffs < uiPartsDiv; uiOffs++)
					{
						pBase[uiAbsPartIdx + uiOffs] |= flag;
					}
				}
			}
		}
#endif  //#if !JVET_C0024_QTBT
	}




	Void
		TEncSearch::xSetIntraResultChromaQT(TComYuv*    pcRecoYuv, TComTU &rTu)
	{
		if (!rTu.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
		{
			return;
		}
		TComDataCU *pcCU = rTu.getCU();
		const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
		const UInt uiTrDepth = rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
		assert(uiTrDepth == 0 && uiAbsPartIdx == 0);
		const UInt    uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		const UInt    uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
		if (uiTrMode == uiTrDepth)
		{
			UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
			UInt uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif

			//===== copy transform coefficients =====
			const TComRectangle &tuRectCb = rTu.getRect(COMPONENT_Cb);
			UInt uiNumCoeffC = tuRectCb.width*tuRectCb.height;//( pcCU->getSlice()->getSPS()->getMaxCUWidth() * pcCU->getSlice()->getSPS()->getMaxCUHeight() ) >> ( uiFullDepth << 1 );
			const UInt offset = rTu.getCoefficientOffset(COMPONENT_Cb);

			const UInt numberValidComponents = getNumberValidComponents(rTu.GetChromaFormat());
			for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
			{
				const ComponentID component = ComponentID(ch);
#if JVET_C0024_QTBT
				assert(offset == 0);
				const TCoeff* src = m_pppcQTTempCoeff[component][uiWIdx][uiHIdx] + offset;//(uiNumCoeffIncC*uiAbsPartIdx);
#else
				const TCoeff* src = m_ppcQTTempCoeff[component][uiQTLayer] + offset;//(uiNumCoeffIncC*uiAbsPartIdx);
#endif
				TCoeff* dest = pcCU->getCoeff(component) + offset;//(uiNumCoeffIncC*uiAbsPartIdx);
				::memcpy(dest, src, sizeof(TCoeff)*uiNumCoeffC);
#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
				TCoeff* pcArlCoeffSrc = m_pppcQTTempArlCoeff[component][uiWIdx][uiHIdx] + offset;//( uiNumCoeffIncC * uiAbsPartIdx );
#else
				TCoeff* pcArlCoeffSrc = m_ppcQTTempArlCoeff[component][uiQTLayer] + offset;//( uiNumCoeffIncC * uiAbsPartIdx );
#endif
				TCoeff* pcArlCoeffDst = pcCU->getArlCoeff(component) + offset;//( uiNumCoeffIncC * uiAbsPartIdx );
				::memcpy(pcArlCoeffDst, pcArlCoeffSrc, sizeof(TCoeff) * uiNumCoeffC);
#endif
			}

			//===== copy reconstruction =====
#if JVET_C0024_QTBT
			m_ppcQTTempTComYuv[uiWIdx][uiHIdx].copyPartToPartComponent(COMPONENT_Cb, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height);
			m_ppcQTTempTComYuv[uiWIdx][uiHIdx].copyPartToPartComponent(COMPONENT_Cr, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height);
#else
			m_pcQTTempTComYuv[uiQTLayer].copyPartToPartComponent(COMPONENT_Cb, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height);
			m_pcQTTempTComYuv[uiQTLayer].copyPartToPartComponent(COMPONENT_Cr, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height);
		}
		else
		{
			TComTURecurse tuRecurseChild(rTu, false);
			do
			{
				xSetIntraResultChromaQT(pcRecoYuv, tuRecurseChild);
			} while (tuRecurseChild.nextSection(rTu));
		}
#endif
	}



	Void
		TEncSearch::estIntraPredLumaQT(TComDataCU* pcCU,
			TComYuv*    pcOrgYuv,
			TComYuv*    pcPredYuv,
			TComYuv*    pcResiYuv,
			TComYuv*    pcRecoYuv,
			// #if DIMD_MLIP_HARM
			//     Bool bUseDIMD,
			// #endif
#if COM16_C806_LARGE_CTU
			Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
			Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
			DEBUG_STRING_FN_DECLARE(sDebug))
	{
		const UInt         uiDepth = pcCU->getDepth(0);
#if JVET_C0024_QTBT
		const UInt    uiInitTrDepth = 0;
		const UInt    uiWidth = pcCU->getWidth(0);
		const UInt    uiHeight = pcCU->getHeight(0);
		const UInt    uiWIdx = g_aucConvertToBit[uiWidth];
		const UInt    uiHIdx = g_aucConvertToBit[uiHeight];
#if BARC
		Bool    isBARC = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#else
		const UInt         uiInitTrDepth = pcCU->getPartitionSize(0) == SIZE_2Nx2N ? 0 : 1;
		const UInt         uiNumPU = 1 << (2 * uiInitTrDepth);
		const UInt         uiQNumParts = pcCU->getTotalNumPart() >> 2;
		const UInt         uiWidthBit = pcCU->getIntraSizeIdx(0);
#endif

		const ChromaFormat chFmt = pcCU->getPic()->getChromaFormat();
		const UInt         numberValidComponents = getNumberValidComponents(chFmt);
		const TComSPS     &sps = *(pcCU->getSlice()->getSPS());
		const TComPPS     &pps = *(pcCU->getSlice()->getPPS());
		Distortion   uiOverallDistY = 0;
		UInt         CandNum;
		Double       CandCostList[FAST_UDI_MAX_RDMODE_NUM];
#if JVET_C0024_PBINTRA_FAST
		Double    CandHadList[FAST_UDI_MAX_RDMODE_NUM];
		for (Int i = 0; i<FAST_UDI_MAX_RDMODE_NUM; i++)
		{
			CandHadList[i] = MAX_DOUBLE;
		}
#endif

#if HAMONIZATION_FOR_MLIP_AND_D0127
		UInt         TempCandNum;
		Double       TempCandCostList[FAST_UDI_MAX_RDMODE_NUM];
#if JVET_C0024_PBINTRA_FAST
		Double    TempCandHadList[FAST_UDI_MAX_RDMODE_NUM];
		for (Int i = 0; i < FAST_UDI_MAX_RDMODE_NUM; i++)
		{
			TempCandHadList[i] = MAX_DOUBLE;
		}
#endif
#endif

#if !COM16_C806_LARGE_CTU
		Pel          resiLumaPU[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
#endif

		Bool    bMaintainResidual[NUMBER_OF_STORED_RESIDUAL_TYPES];
		for (UInt residualTypeIndex = 0; residualTypeIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; residualTypeIndex++)
		{
			bMaintainResidual[residualTypeIndex] = true; //assume true unless specified otherwise
		}

		bMaintainResidual[RESIDUAL_ENCODER_SIDE] = !(m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate());

		// Lambda calculation at equivalent Qp of 4 is recommended because at that Qp, the quantisation divisor is 1.
#if FULL_NBIT
		const Double sqrtLambdaForFirstPass = (m_pcEncCfg->getCostMode() == COST_MIXED_LOSSLESS_LOSSY_CODING && pcCU->getCUTransquantBypass(0)) ?
			sqrt(0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12) / 3.0)))
			: m_pcRdCost->getSqrtLambda();
#else
		const Double sqrtLambdaForFirstPass = (m_pcEncCfg->getCostMode() == COST_MIXED_LOSSLESS_LOSSY_CODING && pcCU->getCUTransquantBypass(0)) ?
			sqrt(0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12 - 6 * (sps.getBitDepth(CHANNEL_TYPE_LUMA) - 8)) / 3.0)))
			: m_pcRdCost->getSqrtLambda();
#endif

		//===== set QP and clear Cbf =====
		if (pps.getUseDQP() == true)
		{
#if JVET_C0024_DELTA_QP_FIX
			pcCU->setQPSubParts(pcCU->getQP(0), 0, pcCU->getWidth(0), pcCU->getHeight(0));
#else
			pcCU->setQPSubParts(pcCU->getQP(0), 0, uiDepth);
#endif
		}
		else
		{
#if JVET_C0024_DELTA_QP_FIX
#if BARC_AdpQP
			pcCU->setQPSubParts(pcCU->getQP(0), 0, pcCU->getWidth(0), pcCU->getHeight(0));
#else
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), 0, pcCU->getWidth(0), pcCU->getHeight(0));
#endif
#else
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), 0, uiDepth);
#endif
		}

#if JVET_D0127_REDUNDANCY_REMOVAL||HAMONIZATION_FOR_MLIP_AND_D0127
#if JVET_C0024_QTBT
		Bool NSSTFlag = (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) == 0) || (uiWidth > 64 || uiHeight > 64);
		Bool NSSTSaveFlag = (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) == 0)
#if COM16_C1046_PDPC_INTRA && !JVET_G0104_PLANAR_PDPC
			&& (pcCU->getPDPCIdx(0) == 0)
#endif
			;
#else
		Bool NSSTFlag = (pcCU->getROTIdx(0) == 0);
		Bool NSSTSaveFlag = (pcCU->getROTIdx(0) == 0)
#if COM16_C1046_PDPC_INTRA
			&& (pcCU->getPDPCIdx(0) == 0)
#endif
#if COM16_C806_EMT
			&& (pcCU->getEmtCuFlag(0) == 0)
#endif
			;
#endif
		static UInt   uiSavedRdModeListNSST[35], uiSavedNumRdModesNSST, uiSavedHadModeListNSST[35];
		static Double dSavedModeCostNSST[35], dSavedHadListNSST[FAST_UDI_MAX_RDMODE_NUM];
#if JVET_C0024_PBINTRA_FAST
		UInt uiHadModeList[67];
#endif
#endif

#if COM16_C806_EMT
#if JVET_C0024_QTBT
		static Double dBestModeCostStore; // RD cost of the best mode for each PU using DCT2
		static Double dModeCostStore[35]; // RD cost of each mode for each PU using DCT2
		static UInt   uiSavedRdModeList[35], uiSavedNumRdModes;
#if MLIP_EXTEND_REF_LINE
		static UInt   uiSavedExternRefList[FAST_UDI_MAX_RDMODE_NUM];
#endif
#else
		static Double dBestModeCostStore[4]; // RD cost of the best mode for each PU using DCT2
		static Double dModeCostStore[4][35]; // RD cost of each mode for each PU using DCT2
		static UInt   uiSavedRdModeList[4][35], uiSavedNumRdModes[4];
#endif

		// Marking EMT usage for faster EMT
		// 0: EMT not applicable for current CU (pcCU->getWidth(0) <= EMT_INTRA_MAX_CU)
		// 1: EMT can be applied for current CU, and DCT2 is being checked
		// 2: EMT is being checked for current CU. Stored results of DCT2 can be utilized for speedup
#if JVET_C0024_QTBT
		UChar ucEmtUsageFlag = ((uiWidth <= EMT_INTRA_MAX_CU && uiHeight <= EMT_INTRA_MAX_CU) ? (pcCU->getEmtCuFlag(0) == 1 ? 2 : 1) : 0);
#else
		UChar ucEmtUsageFlag = (pcCU->getWidth(0) <= EMT_INTRA_MAX_CU ? (pcCU->getEmtCuFlag(0) == 1 ? 2 : 1) : 0);
#endif
		Bool  bAllIntra = (m_pcEncCfg->getIntraPeriod() == 1);

#if JVET_C0024_QTBT
		if (uiWidth*uiHeight<64 && !bAllIntra)
#else
		if (pcCU->getPartitionSize(0) == SIZE_NxN && !bAllIntra)
#endif
		{
			ucEmtUsageFlag = 0;
		}
#if JVET_D0077_SAVE_LOAD_ENC_INFO
		UInt  uiZorderIdx = pcCU->getZorderIdxInCtu();
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
		UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, pcCU->getCURouteIdx(), uiWIdx, uiHIdx);
#else
		UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, uiWIdx, uiHIdx);
#endif
		if (saveLoadTag == LOAD_ENC_INFO /*&& getSaveLoadEmtFlag(uiWIdx, uiHIdx) == 0*/)
		{
			ucEmtUsageFlag = 0;
		}
#if JVET_D0127_REDUNDANCY_REMOVAL||HAMONIZATION_FOR_MLIP_AND_D0127
		NSSTFlag |= saveLoadTag == LOAD_ENC_INFO;
#endif
#endif
#endif

		//===== loop over partitions =====
		TComTURecurse tuRecurseCU(pcCU, 0);
		TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

#if COM16_C806_EMT
		UInt uiPU = 0;
#endif
		do
		{
			const UInt uiPartOffset = tuRecurseWithPU.GetAbsPartIdxTU();
			//  for( UInt uiPU = 0, uiPartOffset=0; uiPU < uiNumPU; uiPU++, uiPartOffset += uiQNumParts )
			//{

#if COM16_C806_EMT
#if MLIP_EXTEND_REF_LINE
			Double MPMCostList[2 * NUM_MOST_PROBABLE_MODES];
			Bool MPMCalc[2 * NUM_MOST_PROBABLE_MODES];
			for (UInt x = 0; x < 2 * NUM_MOST_PROBABLE_MODES; x++) {
				MPMCostList[x] = MAX_DOUBLE;
				MPMCalc[x] = false;
			}
			UInt MPMCandList[2 * NUM_MOST_PROBABLE_MODES];
			UInt MPMList[2 * NUM_MOST_PROBABLE_MODES];

			UInt uiExternRefList[FAST_UDI_MAX_RDMODE_NUM];
			UInt iNumOfPassesExternRef = MLIP_REF_LINE_NUM;
			if (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiPartOffset)
#if RRC_INTRA_SHARPENING
				|| pcCU->getIntraSharpFlag(0)
#endif
#if RRC_DWDIP
				|| pcCU->getDWDIPIdx(0)
#endif
#if RRC_DWDIP_MLR
				|| pcCU->getDWDIP_GPIdx(0)
#endif
#if DIMD_INTRA_PRED
				|| pcCU->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, 0)
#endif
				)
				iNumOfPassesExternRef = 1;

			Bool bSPSMLIPFlag = pcCU->getSlice()->getSPS()->getUseMLIP();
			if (!bSPSMLIPFlag)// if SPS flag is 0, then the num is 1,(not use MLIP)
			{
				iNumOfPassesExternRef = 1;
			}

#endif

			UInt uiRdModeList[FAST_UDI_MAX_RDMODE_NUM];

#if HAMONIZATION_FOR_MLIP_AND_D0127
			UInt uiTempRdModeList[FAST_UDI_MAX_RDMODE_NUM];
#endif
#if JVET_C0024_QTBT
			assert(uiPartOffset == 0);
			Int numModesForFullRD = g_aucIntraModeNumFast_UseMPM[uiWIdx][uiHIdx];
#if HAMONIZATION_FOR_MLIP_AND_D0127
			Int tempNumModesForFullRD = g_aucIntraModeNumFast_UseMPM[uiWIdx][uiHIdx];
			Int curBlkROTIdx = pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiPartOffset);
#endif

#if MLIP_EXTEND_REF_LINE
			numModesForFullRD *= iNumOfPassesExternRef;
#endif

#else
			Int numModesForFullRD = m_pcEncCfg->getFastUDIUseMPMEnabled() ? g_aucIntraModeNumFast_UseMPM[uiWidthBit] : g_aucIntraModeNumFast_NotUseMPM[uiWidthBit];
#if VCEG_AZ07_INTRA_65ANG_MODES
			numModesForFullRD -= 1;
#endif
#endif

#if MLIP_EXTEND_REF_LINE
			UInt OrinumModesForFullRD = numModesForFullRD;
#endif

			if (ucEmtUsageFlag != 2)
			{
#endif
				//===== init pattern for luma prediction =====
				DEBUG_STRING_NEW(sTemp2)

					//===== determine set of modes to be tested (using prediction signal only) =====
#if VCEG_AZ07_INTRA_65ANG_MODES
					Int numModesAvailable = 67; //total number of Intra modes
#else
					Int numModesAvailable = 35; //total number of Intra modes
#endif
#if !COM16_C806_EMT
				UInt uiRdModeList[FAST_UDI_MAX_RDMODE_NUM];
#if JVET_C0024_QTBT
				assert(uiPartOffset == 0);
				Int numModesForFullRD = g_aucIntraModeNumFast_UseMPM[uiWIdx][uiHIdx];
#else
				Int numModesForFullRD = m_pcEncCfg->getFastUDIUseMPMEnabled() ? g_aucIntraModeNumFast_UseMPM[uiWidthBit] : g_aucIntraModeNumFast_NotUseMPM[uiWidthBit];
#if VCEG_AZ07_INTRA_65ANG_MODES
				numModesForFullRD -= 1;
#endif
#endif
#endif

				// this should always be true
				assert(tuRecurseWithPU.ProcessComponentSection(COMPONENT_Y));
#if !MLIP_EXTEND_REF_LINE
#if BARC_CORE
				Bool isCTUBoundary = false;
				if (isBARC)
				{
					UInt leftBoundary = (pcCU->getCUPelX() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
					UInt aboveBoundary = (pcCU->getCUPelY() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
					isCTUBoundary = leftBoundary + aboveBoundary > 0 ? 1 : 0;
				}
#endif
				initIntraPatternChType(tuRecurseWithPU, COMPONENT_Y, true
#if RRC_INTRA_SHARPENING
					, false
#endif
#if BARC_CORE
					, isCTUBoundary
#endif 
					DEBUG_STRING_PASS_INTO(sTemp2));
#endif

				Bool doFastSearch = (numModesForFullRD != numModesAvailable);
				if (doFastSearch)
				{
					assert(numModesForFullRD < numModesAvailable);

#if JVET_D0127_REDUNDANCY_REMOVAL
					for (Int i = 0; i < numModesForFullRD + 2; i++)
#else
					for (Int i = 0; i < numModesForFullRD; i++)
#endif
					{
						CandCostList[i] = MAX_DOUBLE;
					}
					CandNum = 0;
#if HAMONIZATION_FOR_MLIP_AND_D0127
					for (Int i = 0; i < tempNumModesForFullRD + 2; i++)
					{
						TempCandCostList[i] = MAX_DOUBLE;
					}
					TempCandNum = 0;
#endif
#if VCEG_AZ07_INTRA_65ANG_MODES
					Int uiPreds[6] = { -1, -1, -1, -1, -1, -1 };
					Int iAboveLeftCase = 0, iMode = -1;

#if MLIP_EXTEND_REF_LINE
#if BARC
					Bool isCTUBoundary = false;
					if (isBARC)
					{
						UInt leftBoundary = (pcCU->getCUPelX() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
						UInt aboveBoundary = (pcCU->getCUPelY() % pcCU->getSlice()->getSPS()->getCTUSize()) == 0 ? 1 : 0;
						isCTUBoundary = leftBoundary + aboveBoundary > 0 ? 1 : 0;
					}
#endif
					for (Int ExternRef = 0; ExternRef < iNumOfPassesExternRef; ExternRef++) {

						pcCU->setExternRefSubParts(CHANNEL_TYPE_LUMA, ExternRef, uiPartOffset, uiDepth + uiInitTrDepth);
						initIntraPatternChType(tuRecurseWithPU, COMPONENT_Y, true
#if RRC_INTRA_SHARPENING
							, false
#endif
#if BARC_CORE
							, isCTUBoundary
#endif 
							DEBUG_STRING_PASS_INTO(sTemp2));
#endif


#if JVET_C0055_INTRA_MPM
						pcCU->getIntraDirPredictor(uiPartOffset, uiPreds, COMPONENT_Y,
#if MLIP_EXTEND_REF_LINE
							ExternRef,
#endif
							&iMode); // Pre-calculate the MPMs, so avoid redundant MPM calculations during the SATD loop
#else
						pcCU->getIntraDirPredictor(uiPartOffset, uiPreds, COMPONENT_Y, iAboveLeftCase, &iMode); // Pre-calculate the MPMs, so avoid redundant MPM calculations during the SATD loop
#endif

#if MLIP_EXTEND_REF_LINE
						if (ExternRef < 2) {
							for (UInt x = 0; x < NUM_MOST_PROBABLE_MODES; x++)
								MPMList[NUM_MOST_PROBABLE_MODES*ExternRef + x] = uiPreds[x];
						}
#endif

						assert(iMode >= 0);

						iMode = min(iMode + 1, 6);

						Char bSatdChecked[NUM_INTRA_MODE];
						memset(bSatdChecked, 0, NUM_INTRA_MODE * sizeof(Char));

#endif

						const TComRectangle &puRect = tuRecurseWithPU.getRect(COMPONENT_Y);
						const UInt uiAbsPartIdx = tuRecurseWithPU.GetAbsPartIdxTU();

						Pel* piOrg = pcOrgYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);
						Pel* piPred = pcPredYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);
						UInt uiStride = pcPredYuv->getStride(COMPONENT_Y);
						DistParam distParam;
						const Bool bUseHadamard = pcCU->getCUTransquantBypass(0) == 0;
						m_pcRdCost->setDistParam(distParam, sps.getBitDepth(CHANNEL_TYPE_LUMA), piOrg, uiStride, piPred, uiStride, puRect.width, puRect.height, bUseHadamard);
						distParam.bApplyWeight = false;

#if JVET_D0127_REDUNDANCY_REMOVAL||HAMONIZATION_FOR_MLIP_AND_D0127
						if (NSSTFlag) {
#endif
							for (Int modeIdx = 0; modeIdx < numModesAvailable; modeIdx++)
							{
								UInt       uiMode = modeIdx;
								Distortion uiSad = 0;

#if COM16_C1044_NSST
#if !JVET_C0024_QTBT 
								if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
#endif
								{
#if FAST_NSST_v1
									const Int iNumberOfPassesROT = (uiMode <= DC_IDX) ? 2 : 4;
#else
									const Int iNumberOfPassesROT = (uiMode <= DC_IDX) ? 3 : 4;
#endif
#if JVET_C0024_QTBT
									if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0))
#else
									if (iNumberOfPassesROT <= pcCU->getROTIdx(0))
#endif
									{
										continue;
									}
								}
#endif
#if RRC_INTRA_SHARPENING
								if (pcCU->getIntraSharpFlag(uiPartOffset)
#if !JVET_C0024_QTBT
									&& pcCU->getPartitionSize(uiPartOffset) == SIZE_2Nx2N
#endif
									)
								{
									if (uiMode == DC_IDX || uiMode == PLANAR_IDX)
									{
										continue;
									}
								}
#endif   

#if VCEG_AZ07_INTRA_65ANG_MODES
								if (uiMode>DC_IDX && (uiMode & 1))
								{
									// Skip checking extended Angular modes in the first round of SATD
									continue;
								}
								bSatdChecked[uiMode] = true;
#endif

#if MLIP_EXTEND_REF_LINE
								if (ExternRef && (uiMode == DC_IDX || uiMode == PLANAR_IDX))
									continue;
#endif
								const Bool bUseFilter = TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, uiMode, puRect.width, puRect.height, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag()
#if DIMD_NUM_INTRA_DIR_INC
									, pcCU, uiAbsPartIdx
#endif
								);

								predIntraAng(COMPONENT_Y, uiMode, piOrg, uiStride, piPred, uiStride, tuRecurseWithPU, bUseFilter
									, TComPrediction::UseDPCMForFirstPassIntraEstimation(tuRecurseWithPU, uiMode));

								// use hadamard transform here
								uiSad += distParam.DistFunc(&distParam);

								UInt   iModeBits = 0;

								// NB xModeBitsIntra will not affect the mode for chroma that may have already been pre-estimated.
								iModeBits += xModeBitsIntra(pcCU, uiMode, uiPartOffset, uiDepth, CHANNEL_TYPE_LUMA
#if VCEG_AZ07_INTRA_65ANG_MODES
									, uiPreds, iAboveLeftCase
#endif
								);

								Double cost = (Double)uiSad + (Double)iModeBits * sqrtLambdaForFirstPass;

#if DEBUG_INTRA_SEARCH_COSTS
								std::cout << "1st pass mode " << uiMode << " SAD = " << uiSad << ", mode bits = " << iModeBits << ", cost = " << cost << "\n";
#endif

#if JVET_C0024_FAST_MRG
#if HAMONIZATION_FOR_MLIP_AND_D0127

								CandNum += updateCandList(uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList
#if MLIP_EXTEND_REF_LINE
									, uiExternRefList, ExternRef
#endif
								);
								//for harmonization
								if (/*(curBlkROTIdx == 0) && */(ExternRef == 0))
								{
									TempCandNum += updateCandList(uiMode, cost, tempNumModesForFullRD + 2, uiTempRdModeList, TempCandCostList);
								}


#else

#if JVET_D0127_REDUNDANCY_REMOVAL
								CandNum += updateCandList(uiMode, cost, numModesForFullRD + 2, uiRdModeList, CandCostList);
#else
								CandNum += updateCandList(uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList
#if MLIP_EXTEND_REF_LINE
									, uiExternRefList, ExternRef
#endif
								);
#endif
#endif
#else
								CandNum += xUpdateCandList(uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList);
#endif

#if MLIP_EXTEND_REF_LINE
								if (ExternRef < 2) {
									Bool isMPM = false;
									Int MPMIdx = -1;
									for (UInt x = 0; x < NUM_MOST_PROBABLE_MODES; x++) {
										if (uiMode == uiPreds[x]) {
											isMPM = true;
											MPMIdx = x;
											MPMCalc[NUM_MOST_PROBABLE_MODES*ExternRef + x] = true;
										}
									}
									if (isMPM)
										updateCandList(NUM_MOST_PROBABLE_MODES*ExternRef + MPMIdx, cost, 2 * NUM_MOST_PROBABLE_MODES, MPMCandList, MPMCostList);
								}

#endif


#if JVET_C0024_PBINTRA_FAST

#if HAMONIZATION_FOR_MLIP_AND_D0127

								//if ((curBlkROTIdx == 0) && (ExternRef == 0))
								if (ExternRef == 0)
								{
									updateCandList(uiMode, uiSad, tempNumModesForFullRD + 2, uiHadModeList, TempCandHadList);
								}
								//else
								//{
								if (uiSad < CandHadList[0])
								{
									CandHadList[2] = CandHadList[1];
									CandHadList[1] = CandHadList[0];
									CandHadList[0] = uiSad;
								}
								else if (uiSad < CandHadList[1])
								{
									CandHadList[2] = CandHadList[1];
									CandHadList[1] = uiSad;
								}
								else if (uiSad < CandHadList[2])
								{
									CandHadList[2] = uiSad;
								}
								//}

#else

#if JVET_D0127_REDUNDANCY_REMOVAL

								updateCandList(uiMode, uiSad, numModesForFullRD + 2, uiHadModeList, CandHadList);
#else
								if (uiSad < CandHadList[0])
								{
									CandHadList[2] = CandHadList[1];
									CandHadList[1] = CandHadList[0];
									CandHadList[0] = uiSad;
								}
								else if (uiSad < CandHadList[1])
								{
									CandHadList[2] = CandHadList[1];
									CandHadList[1] = uiSad;
								}
								else if (uiSad < CandHadList[2])
								{
									CandHadList[2] = uiSad;
								}
#endif   
#endif
#endif
							}
#if JVET_D0127_REDUNDANCY_REMOVAL||HAMONIZATION_FOR_MLIP_AND_D0127
#if HAMONIZATION_FOR_MLIP_AND_D0127
							//if ((NSSTSaveFlag) && ((curBlkROTIdx == 0) && (ExternRef == 0)))
							if ((NSSTSaveFlag) && (ExternRef == 0))
							{
								uiSavedNumRdModesNSST = tempNumModesForFullRD;
								::memcpy(uiSavedRdModeListNSST, uiTempRdModeList, (tempNumModesForFullRD + 2) * sizeof(UInt));
								::memcpy(dSavedModeCostNSST, TempCandCostList, (tempNumModesForFullRD + 2) * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
								::memcpy(uiSavedHadModeListNSST, uiHadModeList, (tempNumModesForFullRD + 2) * sizeof(UInt));
								::memcpy(dSavedHadListNSST, TempCandHadList, (tempNumModesForFullRD + 2) * sizeof(Double));
#endif
							}
#else
							if (NSSTSaveFlag) {
								uiSavedNumRdModesNSST = numModesForFullRD;
								::memcpy(uiSavedRdModeListNSST, uiRdModeList, (numModesForFullRD + 2) * sizeof(UInt));
								::memcpy(dSavedModeCostNSST, CandCostList, (numModesForFullRD + 2) * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
								::memcpy(uiSavedHadModeListNSST, uiHadModeList, (numModesForFullRD + 2) * sizeof(UInt));
								::memcpy(dSavedHadListNSST, CandHadList, (numModesForFullRD + 2) * sizeof(Double));
#endif
							}
#endif
						}
						else {
#if HAMONIZATION_FOR_MLIP_AND_D0127

							//if ((curBlkROTIdx != 0) && (ExternRef == 0))
							if (ExternRef == 0)
							{
#if JVET_C0024_QTBT
#if FAST_NSST_v1
								if (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) >1){
#else
								if (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) == 3) {
#endif
#else
								if (pcCU->getROTIdx(0) == 3) {
#endif
									numModesForFullRD = uiSavedNumRdModesNSST;
									::memcpy(uiRdModeList, uiSavedRdModeListNSST, (numModesForFullRD + 2) * sizeof(UInt));
									::memcpy(CandCostList, dSavedModeCostNSST, (numModesForFullRD + 2) * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
									::memcpy(uiHadModeList, uiSavedHadModeListNSST, (numModesForFullRD + 2) * sizeof(UInt));
									::memcpy(CandHadList, dSavedHadListNSST, (numModesForFullRD + 2) * sizeof(Double));
#endif
									Int cnt = 0;
									Int i = 0;

									for (i = 0; i < numModesForFullRD; i++) {
										if (uiRdModeList[i] <= DC_IDX) {
											for (UInt j = i; j < numModesForFullRD + 1 - cnt; j++) {
												uiRdModeList[j] = uiRdModeList[j + 1];
												CandCostList[j] = CandCostList[j + 1];
											}
											cnt++;
											i--;
										}
									}

#if JVET_C0024_PBINTRA_FAST
									cnt = 0;
									for (i = 0; i < numModesForFullRD; i++) {
										if (uiHadModeList[i] <= DC_IDX) {
											for (UInt j = i; j < numModesForFullRD + 1 - cnt; j++) {
												uiHadModeList[j] = uiHadModeList[j + 1];
												CandHadList[j] = CandHadList[j + 1];
											}
											cnt++;
											i--;
										}
									}
#endif
								}
								else {
									numModesForFullRD = uiSavedNumRdModesNSST;
									::memcpy(uiRdModeList, uiSavedRdModeListNSST, numModesForFullRD * sizeof(UInt));
									::memcpy(CandCostList, dSavedModeCostNSST, numModesForFullRD * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
									::memcpy(CandHadList, dSavedHadListNSST, numModesForFullRD * sizeof(Double));
#endif
								}

								}
#else
#if JVET_C0024_QTBT
#if FAST_NSST_v1
							if (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) > 1 )
#else
							if (pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0) == 3)
#endif
							{
#else
							if (pcCU->getROTIdx(0) == 3) {

#endif

								numModesForFullRD = uiSavedNumRdModesNSST;
								::memcpy(uiRdModeList, uiSavedRdModeListNSST, (numModesForFullRD + 2) * sizeof(UInt));
								::memcpy(CandCostList, dSavedModeCostNSST, (numModesForFullRD + 2) * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
								::memcpy(uiHadModeList, uiSavedHadModeListNSST, (numModesForFullRD + 2) * sizeof(UInt));
								::memcpy(CandHadList, dSavedHadListNSST, (numModesForFullRD + 2) * sizeof(Double));
#endif
								Int cnt = 0;
								Int i = 0;

								for (i = 0; i < numModesForFullRD; i++) {
									if (uiRdModeList[i] <= DC_IDX) {
										for (UInt j = i; j < numModesForFullRD + 1 - cnt; j++) {
											uiRdModeList[j] = uiRdModeList[j + 1];
											CandCostList[j] = CandCostList[j + 1];
										}
										cnt++;
										i--;
									}
								}

#if JVET_C0024_PBINTRA_FAST
								cnt = 0;
								for (i = 0; i < numModesForFullRD; i++) {
									if (uiHadModeList[i] <= DC_IDX) {
										for (UInt j = i; j < numModesForFullRD + 1 - cnt; j++) {
											uiHadModeList[j] = uiHadModeList[j + 1];
											CandHadList[j] = CandHadList[j + 1];
										}
										cnt++;
										i--;
									}
								}
#endif
							}
							else
							{
								numModesForFullRD = uiSavedNumRdModesNSST;
								::memcpy(uiRdModeList, uiSavedRdModeListNSST, numModesForFullRD * sizeof(UInt));
								::memcpy(CandCostList, dSavedModeCostNSST, numModesForFullRD * sizeof(Double));
#if JVET_C0024_PBINTRA_FAST
								::memcpy(CandHadList, dSavedHadListNSST, numModesForFullRD * sizeof(Double));
#endif
							}

#endif
							}//if       if not  (NSSTFlag)

#endif // if JVET_D0127_REDUNDANCY_REMOVAL||HAMONIZATION_FOR_EXT_REF_LINE_AND_D0127

#if VCEG_AZ07_INTRA_65ANG_MODES
						UInt uiParentCandList[FAST_UDI_MAX_RDMODE_NUM];
						memcpy(uiParentCandList, uiRdModeList, sizeof(UInt)*numModesForFullRD);

#if MLIP_EXTEND_REF_LINE
						UInt uiParentExternRefList[FAST_UDI_MAX_RDMODE_NUM];
						memcpy(uiParentExternRefList, uiExternRefList, sizeof(UInt)*numModesForFullRD);
#if HAMONIZATION_FOR_MLIP_AND_D0127
						if (curBlkROTIdx != 0)// current NSST or ROT not use for further reference line .  
						{
							memset(uiParentExternRefList, 0, sizeof(UInt)*numModesForFullRD);
							memset(uiExternRefList, 0, sizeof(UInt)*numModesForFullRD);
						}

#endif 

#endif

						// Second round of SATD for extended Angular modes
						for (Int modeIdx = 0; modeIdx < numModesForFullRD; modeIdx++)
						{
							UInt uiParentMode = uiParentCandList[modeIdx];
#if MLIP_EXTEND_REF_LINE
							UInt uiExternRef = uiParentExternRefList[modeIdx];
							if (modeIdx > numModesForFullRD * (ExternRef + 1) / iNumOfPassesExternRef)
								break;
							if (ExternRef && uiExternRef != ExternRef)
								continue;
#endif

#if JVET_E0077_ENHANCED_LM
							if (uiParentMode>2 && uiParentMode<(NUM_INTRA_MODE - NUM_INTRA_MODE_NON_ANG - 1))
#else
							if (uiParentMode>2 && uiParentMode<(NUM_INTRA_MODE - 2))
#endif
							{
								for (Int subModeIdx = -1; subModeIdx <= 1; subModeIdx += 2)
								{
									UInt uiMode = uiParentMode + subModeIdx;

#if COM16_C1044_NSST
#if !JVET_C0024_QTBT
									if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
#endif
									{
#if FAST_NSST_v1
										const Int iNumberOfPassesROT = (uiMode <= DC_IDX) ? 2 : 4;
#else
										const Int iNumberOfPassesROT = (uiMode <= DC_IDX) ? 3 : 4;
#endif
#if JVET_C0024_QTBT
										if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0))
#else
										if (iNumberOfPassesROT <= pcCU->getROTIdx(0))
#endif
										{
											continue;
										}
									}
#endif

									if (!bSatdChecked[uiMode])
									{
										const Bool bUseFilter = TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, uiMode, puRect.width, puRect.height, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag()
#if DIMD_NUM_INTRA_DIR_INC
											, pcCU, uiAbsPartIdx
#endif
										);

										predIntraAng(COMPONENT_Y, uiMode, piOrg, uiStride, piPred, uiStride, tuRecurseWithPU, bUseFilter
											, TComPrediction::UseDPCMForFirstPassIntraEstimation(tuRecurseWithPU, uiMode));


										// use hadamard transform here
										Distortion uiSad = distParam.DistFunc(&distParam);
										UInt   iModeBits = xModeBitsIntra(pcCU, uiMode, uiPartOffset, uiDepth, CHANNEL_TYPE_LUMA, uiPreds, iAboveLeftCase);
										Double cost = (Double)uiSad + (Double)iModeBits * sqrtLambdaForFirstPass;

#if DEBUG_INTRA_SEARCH_COSTS
										std::cout << "1st pass mode for extended angular mode " << uiMode << " SAD = " << uiSad << ", mode bits = " << iModeBits << ", cost = " << cost << "\n";
#endif

#if JVET_C0024_FAST_MRG
										CandNum += updateCandList(uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList
#if MLIP_EXTEND_REF_LINE
											, uiExternRefList, ExternRef
#endif
										);
#else
										CandNum += xUpdateCandList(uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList);
#endif
#if MLIP_EXTEND_REF_LINE
										if (ExternRef < 2) {
											Bool isMPM = false;
											Int MPMIdx = -1;
											for (UInt x = 0; x < NUM_MOST_PROBABLE_MODES; x++) {
												if (uiMode == uiPreds[x]) {
													isMPM = true;
													MPMIdx = x;
													MPMCalc[NUM_MOST_PROBABLE_MODES*ExternRef + x] = true;
												}
											}
											if (isMPM)
												updateCandList(NUM_MOST_PROBABLE_MODES*ExternRef + MPMIdx, cost, 2 * NUM_MOST_PROBABLE_MODES, MPMCandList, MPMCostList);
										}
#endif
#if JVET_C0024_PBINTRA_FAST
										if (uiSad < CandHadList[0])
										{
											CandHadList[2] = CandHadList[1];
											CandHadList[1] = CandHadList[0];
											CandHadList[0] = uiSad;
										}
										else if (uiSad < CandHadList[1])
										{
											CandHadList[2] = CandHadList[1];
											CandHadList[1] = uiSad;
										}
										else if (uiSad < CandHadList[2])
										{
											CandHadList[2] = uiSad;
										}
#endif
										bSatdChecked[uiMode] = true; // Mark as checked
									}
								}
							}
						}
#endif

#if MLIP_EXTEND_REF_LINE
						if (ExternRef < 2
#if MLIP_SKIP_NO_USE_PROCESS
							&&curBlkROTIdx == 0
#endif
							)
						{
							for (UInt x = 0; x < NUM_MOST_PROBABLE_MODES; x++) {
								if (MPMCalc[NUM_MOST_PROBABLE_MODES*ExternRef + x] == false) {
									UInt uiMode = uiPreds[x];
									const Bool bUseFilter = TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, uiMode, puRect.width, puRect.height, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag()
#if DIMD_NUM_INTRA_DIR_INC
										, pcCU, uiAbsPartIdx
#endif
									);

									predIntraAng(COMPONENT_Y, uiMode, piOrg, uiStride, piPred, uiStride, tuRecurseWithPU, bUseFilter, TComPrediction::UseDPCMForFirstPassIntraEstimation(tuRecurseWithPU, uiMode));


									// use hadamard transform here
									Distortion uiSad = distParam.DistFunc(&distParam);
									UInt   iModeBits = xModeBitsIntra(pcCU, uiMode, uiPartOffset, uiDepth, CHANNEL_TYPE_LUMA, uiPreds, iAboveLeftCase);
									Double cost = (Double)uiSad + (Double)iModeBits * sqrtLambdaForFirstPass;

#if DEBUG_INTRA_SEARCH_COSTS
									std::cout << "1st pass mode for extended angular mode " << uiMode << " SAD = " << uiSad << ", mode bits = " << iModeBits << ", cost = " << cost << "\n";
#endif

									updateCandList(NUM_MOST_PROBABLE_MODES*ExternRef + x, cost, 2 * NUM_MOST_PROBABLE_MODES, MPMCandList, MPMCostList);
									MPMCalc[NUM_MOST_PROBABLE_MODES*ExternRef + x] = true;
								}
							}
						}

							}
#endif

					if (m_pcEncCfg->getFastUDIUseMPMEnabled())
					{
#if !VCEG_AZ07_INTRA_65ANG_MODES
						Int uiPreds[NUM_MOST_PROBABLE_MODES] = { -1, -1, -1 };
						Int iMode = -1;

						pcCU->getIntraDirPredictor(uiPartOffset, uiPreds, COMPONENT_Y,
#if MLIP_EXTEND_REF_LINE
							externRef,
#endif
							&iMode);
#endif
#if MLIP_EXTEND_REF_LINE
						pcCU->getIntraDirPredictor(uiPartOffset, uiPreds, COMPONENT_Y, 0);
#endif
						const Int numCand = (iMode >= 0) ? iMode : Int(NUM_MOST_PROBABLE_MODES);

						for (Int j = 0; j < numCand; j++)
						{
							Bool mostProbableModeIncluded = false;
							Int mostProbableMode = uiPreds[j];

#if COM16_C1044_NSST
#if !JVET_C0024_QTBT
							if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
#endif
							{
#if FAST_NSST_v1
								const Int iNumberOfPassesROT = (mostProbableMode <= DC_IDX) ? 2 : 4;
#else
								const Int iNumberOfPassesROT = (mostProbableMode <= DC_IDX) ? 3 : 4;
#endif
#if JVET_C0024_QTBT
								if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0))
#else
								if (iNumberOfPassesROT <= pcCU->getROTIdx(0))
#endif
								{
									continue;
								}
							}
#endif

							for (Int i = 0; i < numModesForFullRD; i++)
							{
#if MLIP_EXTEND_REF_LINE
								mostProbableModeIncluded |= (mostProbableMode == uiRdModeList[i] && uiExternRefList[i] == 0);
#else
								mostProbableModeIncluded |= (mostProbableMode == uiRdModeList[i]);
#endif
							}
							if (!mostProbableModeIncluded)
							{
#if MLIP_EXTEND_REF_LINE
								uiExternRefList[numModesForFullRD] = 0;
#endif
								uiRdModeList[numModesForFullRD++] = mostProbableMode;
							}
						}

#if DIMD_INTRA_PRED
#if MLIP_EXTEND_REF_LINE
						if (!bSPSMLIPFlag)//MX: temporally usage, not combined with DIMD.
#endif
						{
							Bool mostProbableModeIncluded = false;
							UInt uiMostProbMode = getDIMDIntraMode();
#if COM16_C1044_NSST
							const Int iNumberOfPassesROT = (uiMostProbMode <= DC_IDX) ? 3 : 4;
							if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0))
							{
								mostProbableModeIncluded = true;
							}
#endif
							for (Int i = 0; i < numModesForFullRD; i++)
							{
								mostProbableModeIncluded |= (uiMostProbMode == uiRdModeList[i]);
							}
							if (!mostProbableModeIncluded)
							{
#if MLIP_EXTEND_REF_LINE
								uiExternRefList[numModesForFullRD] = 0;
#endif
								uiRdModeList[numModesForFullRD++] = uiMostProbMode;
							}
						}
#endif
					}

#if MLIP_EXTEND_REF_LINE
					if (iNumOfPassesExternRef == 4
#if MLIP_SKIP_NO_USE_PROCESS
						&&curBlkROTIdx == 0
#endif
						)
					{

						Int passnum = (g_aucIntraModeNumFast_UseMPM[uiWIdx][uiHIdx] == 2) ? 1 : 3;

						Int pass = 0;
						Int numCand = 0;

						while (pass < passnum && numCand < 2 * NUM_MOST_PROBABLE_MODES) {
							Bool mostProbableModeIncluded = false;
							Int mostProbableMode = MPMList[MPMCandList[numCand]];

#if COM16_C1044_NSST
#if !JVET_C0024_QTBT
							if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
#endif
							{
								const Int iNumberOfPassesROT = (mostProbableMode <= DC_IDX) ? 3 : 4;
#if JVET_C0024_QTBT
								if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0))
#else
								if (iNumberOfPassesROT <= pcCU->getROTIdx(0))
#endif
								{
									continue;
								}
							}
#endif

							for (Int i = 0; i < numModesForFullRD; i++)
							{
#if MLIP_EXTEND_REF_LINE
								mostProbableModeIncluded |= (mostProbableMode == uiRdModeList[i] && (MPMCandList[numCand] > 5) == uiExternRefList[i]);
#else
								mostProbableModeIncluded |= (mostProbableMode == uiRdModeList[i]);
#endif
							}
							if (!mostProbableModeIncluded)
							{
								Bool isMPM = false;
								for (Int x = 0; x < numCand; x++) {
									if (uiRdModeList[OrinumModesForFullRD - 1 - pass] == MPMList[MPMCandList[x]] && (MPMCandList[x] > 5) == uiExternRefList[OrinumModesForFullRD - 1 - pass])
										isMPM = true;
								}
								if (!isMPM) {
									uiExternRefList[OrinumModesForFullRD - 1 - pass] = (MPMCandList[numCand] > 5);
									uiRdModeList[OrinumModesForFullRD - 1 - pass] = mostProbableMode;
								}
								pass++;
							}

							numCand++;
						}
					}
#endif

						}
				else
				{
					for (Int i = 0; i < numModesForFullRD; i++)
					{
						uiRdModeList[i] = i;
					}
				}
#if COM16_C806_EMT
				if (ucEmtUsageFlag == 1)
				{
					// Store the modes to be checked with RD
#if JVET_C0024_QTBT
					uiSavedNumRdModes = numModesForFullRD;
					::memcpy(uiSavedRdModeList, uiRdModeList, numModesForFullRD * sizeof(UInt));
#if MLIP_EXTEND_REF_LINE
					::memcpy(uiSavedExternRefList, uiExternRefList, numModesForFullRD * sizeof(UInt));
#endif
#else
					uiSavedNumRdModes[uiPU] = numModesForFullRD;
					::memcpy(uiSavedRdModeList[uiPU], uiRdModeList, numModesForFullRD * sizeof(UInt));
#endif
				}
					}
			else if (ucEmtUsageFlag == 2)
			{
				if (bAllIntra && m_pcEncCfg->getUseFastIntraEMT())
				{
#if !JVET_C0024_QTBT
					UInt uiWidth = pcCU->getWidth(0) >> uiInitTrDepth;
#endif
#if JVET_C0024_QTBT
					const Double dThrFastMode = 1.0 + 1.4 / sqrt((Double)(uiWidth*uiHeight));
#else
					double dThrFastMode;

					switch (uiWidth)
					{
					case  4: dThrFastMode = 1.47; break; // Skip checking   4x4 Intra modes using the R-D cost in the DCT2-pass 
					case  8: dThrFastMode = 1.28; break; // Skip checking   8x8 Intra modes using the R-D cost in the DCT2-pass
					case 16: dThrFastMode = 1.12; break; // Skip checking 16x16 Intra modes using the R-D cost in the DCT2-pass
					case 32: dThrFastMode = 1.06; break; // Skip checking 32x32 Intra modes using the R-D cost in the DCT2-pass
					default: dThrFastMode = 1.06; break; // Skip checking 32x32 Intra modes using the R-D cost in the DCT2-pass
					}
#endif
					numModesForFullRD = 0;

					// Skip checking the modes with much larger R-D cost than the best mode
#if JVET_C0024_QTBT
					for (Int i = 0; i < uiSavedNumRdModes; i++)
					{
						if (dModeCostStore[i] <= dThrFastMode * dBestModeCostStore)
						{
#if MLIP_EXTEND_REF_LINE
							uiExternRefList[numModesForFullRD] = uiSavedExternRefList[i];
#endif
							uiRdModeList[numModesForFullRD++] = uiSavedRdModeList[i];
						}
					}
#else
					for (Int i = 0; i < uiSavedNumRdModes[uiPU]; i++)
					{
						if (dModeCostStore[uiPU][i] <= dThrFastMode * dBestModeCostStore[uiPU])
						{
							uiRdModeList[numModesForFullRD++] = uiSavedRdModeList[uiPU][i];
						}
					}
#endif
				}
				else
				{
					// Restore the modes to be checked with RD
#if JVET_C0024_QTBT
					numModesForFullRD = uiSavedNumRdModes;
					::memcpy(uiRdModeList, uiSavedRdModeList, numModesForFullRD * sizeof(UInt));
#if MLIP_EXTEND_REF_LINE
					::memcpy(uiExternRefList, uiSavedExternRefList, numModesForFullRD * sizeof(UInt));
#endif
#else
					numModesForFullRD = uiSavedNumRdModes[uiPU];
					::memcpy(uiRdModeList, uiSavedRdModeList[uiPU], numModesForFullRD * sizeof(UInt));
#endif
				}
			}
#endif

			//===== check modes (using r-d costs) =====
#if HHI_RQT_INTRA_SPEEDUP_MOD
			UInt   uiSecondBestMode = MAX_UINT;
			Double dSecondBestPUCost = MAX_DOUBLE;
#endif
			DEBUG_STRING_NEW(sPU)
				UInt       uiBestPUMode = 0;
			Distortion uiBestPUDistY = 0;
			Double     dBestPUCost = MAX_DOUBLE;
#if MLIP_EXTEND_REF_LINE
			UInt       uiBestExternRef = 0;
#endif

#if JVET_C0024_PBINTRA_FAST
			if (pcCU->getSlice()->getSliceType() != I_SLICE
#if COM16_C806_EMT
				&& ucEmtUsageFlag != 2
#endif
				)
			{
				if (CandHadList[2] > (Double)pcCU->getInterHAD()*PBINTRA_RATIO)
				{
					numModesForFullRD = 2;
				}
				if (CandHadList[1] > (Double)pcCU->getInterHAD()*PBINTRA_RATIO)
				{
					numModesForFullRD = 1;
				}
				if (CandHadList[0] > (Double)pcCU->getInterHAD()*PBINTRA_RATIO)
				{
					pcCU->getTotalDistortion() = MAX_UINT;
					pcCU->getInterHAD() = 0;
					return;
				}
			}
#endif

#if RRC_DWDIP_MLR
			if (pcCU->getDWDIP_GPIdx(uiPartOffset))
			{
				numModesForFullRD = 1;
			}
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
			UInt max = numModesForFullRD;

			if (DebugOptionList::ForceLumaMode.isSet())
			{
				max = 0;  // we are forcing a direction, so don't bother with mode check
			}
			for (UInt uiMode = 0; uiMode < max; uiMode++)
#else
			for (UInt uiMode = 0; uiMode < numModesForFullRD; uiMode++)
#endif
			{
				// set luma prediction mode
				UInt uiOrgMode = uiRdModeList[uiMode];
#if MLIP_EXTEND_REF_LINE
				UInt ExternRef = uiExternRefList[uiMode];
				pcCU->setExternRefSubParts(CHANNEL_TYPE_LUMA, ExternRef, uiPartOffset, uiDepth + uiInitTrDepth);

#endif

#if RRC_INTRA_SHARPENING
				if (pcCU->getIntraSharpFlag(uiPartOffset)
#if !JVET_C0024_QTBT
					&& pcCU->getPartitionSize(uiPartOffset) == SIZE_2Nx2N
#endif
					)
				{
					if (uiOrgMode == DC_IDX || uiOrgMode == PLANAR_IDX)
					{
						continue;
					}
				}
#endif   

				pcCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiOrgMode, uiPartOffset, uiDepth + uiInitTrDepth);
#if DIMD_INTRA_PRED
				pcCU->setSecIntraDirSubParts(uiOrgMode, uiPartOffset, uiDepth + uiInitTrDepth);
#endif
				DEBUG_STRING_NEW(sMode)
					// set context models
#if JVET_C0024_QTBT
					m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
#endif

				// determine residual for partition
				Distortion uiPUDistY = 0;
				Double     dPUCost = 0.0;

#if HHI_RQT_INTRA_SPEEDUP
				xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
					m_resiPUBuffer,
#else
					resiLumaPU,
#endif
					uiPUDistY, true, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode));
#else
				xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
					m_resiPUBuffer,
#else
					resiLumaPU,

#endif
					uiPUDistY, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode));
#endif

#if DEBUG_INTRA_SEARCH_COSTS
				std::cout << "2nd pass [luma,chroma] mode [" << Int(pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiPartOffset)) << "," << Int(pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiPartOffset)) << "] cost = " << dPUCost << "\n";
#endif

#if COM16_C806_EMT
				if (1 == ucEmtUsageFlag && m_pcEncCfg->getUseFastIntraEMT())
				{
#if JVET_C0024_QTBT
					dModeCostStore[uiMode] = dPUCost;
#else
					dModeCostStore[uiPU][uiMode] = dPUCost;
#endif
				}
#endif

				// check r-d cost
				if (dPUCost < dBestPUCost)
				{
#if !JVET_C0024_QTBT
					if (isRSAFEnabled)
					{
						isBestRSAF = (isNonRSAF == 0);
					}
#endif

					DEBUG_STRING_SWAP(sPU, sMode)
#if HHI_RQT_INTRA_SPEEDUP_MOD
						uiSecondBestMode = uiBestPUMode;
					dSecondBestPUCost = dBestPUCost;
#endif
					uiBestPUMode = uiOrgMode;
					uiBestPUDistY = uiPUDistY;
					dBestPUCost = dPUCost;

#if MLIP_EXTEND_REF_LINE
					uiBestExternRef = ExternRef;
#endif

#if COM16_C806_EMT
					if (1 == ucEmtUsageFlag && m_pcEncCfg->getUseFastIntraEMT())
					{
#if JVET_C0024_QTBT
						dBestModeCostStore = dPUCost;
#else
						dBestModeCostStore[uiPU] = dPUCost;
#endif
					}
#endif
					xSetIntraResultLumaQT(pcRecoYuv, tuRecurseWithPU);

					if (pps.getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
					{
						const Int xOffset = tuRecurseWithPU.getRect(COMPONENT_Y).x0;
						const Int yOffset = tuRecurseWithPU.getRect(COMPONENT_Y).y0;
						for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
						{
							if (bMaintainResidual[storedResidualIndex])
							{
								xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
									m_resiPUBuffer[storedResidualIndex],
#else
									resiLumaPU[storedResidualIndex],
#endif
									tuRecurseWithPU, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
							}
						}
					}

					UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();
#if COM16_C806_EMT
					::memcpy(m_puhQTTempEmtTuIdx, pcCU->getEmtTuIdx() + uiPartOffset, uiQPartNum * sizeof(UChar));
					::memcpy(m_puhQTTempEmtCuFlag, pcCU->getEmtCuFlag() + uiPartOffset, uiQPartNum * sizeof(UChar));
#endif
#if !JVET_C0024_QTBT
					::memcpy(m_puhQTTempTrIdx, pcCU->getTransformIdx() + uiPartOffset, uiQPartNum * sizeof(UChar));
#endif
					for (UInt component = 0; component < numberValidComponents; component++)
					{
						const ComponentID compID = ComponentID(component);
						::memcpy(m_puhQTTempCbf[compID], pcCU->getCbf(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
						::memcpy(m_puhQTTempTransformSkipFlag[compID], pcCU->getTransformSkip(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
#if VCEG_AZ08_INTRA_KLT
						::memcpy(m_puhQTTempKLTFlag[compID], pcCU->getKLTFlag(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
#endif
					}
				}
#if HHI_RQT_INTRA_SPEEDUP_MOD
				else if (dPUCost < dSecondBestPUCost)
				{
					uiSecondBestMode = uiOrgMode;
					dSecondBestPUCost = dPUCost;
				}
#endif
			} // Mode loop

#if !JVET_C0024_QTBT

#if HHI_RQT_INTRA_SPEEDUP
#if HHI_RQT_INTRA_SPEEDUP_MOD
			for (UInt ui = 0; ui < 2; ++ui)
#endif
			{
#if HHI_RQT_INTRA_SPEEDUP_MOD
				UInt uiOrgMode = ui ? uiSecondBestMode : uiBestPUMode;
				if (uiOrgMode == MAX_UINT)
				{
					break;
				}
#else
				UInt uiOrgMode = uiBestPUMode;
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
				if (DebugOptionList::ForceLumaMode.isSet())
				{
					uiOrgMode = DebugOptionList::ForceLumaMode.getInt();
				}
#endif

				pcCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiOrgMode, uiPartOffset, uiDepth + uiInitTrDepth);
				DEBUG_STRING_NEW(sModeTree)

					// set context models
					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

				// determine residual for partition
				Distortion uiPUDistY = 0;
				Double     dPUCost = 0.0;

				xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
					m_resiPUBuffer,
#else
					resiLumaPU,
#endif
					uiPUDistY, false, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sModeTree));

				// check r-d cost
				if (dPUCost < dBestPUCost)
				{
					DEBUG_STRING_SWAP(sPU, sModeTree)
						uiBestPUMode = uiOrgMode;
					uiBestPUDistY = uiPUDistY;
					dBestPUCost = dPUCost;

					xSetIntraResultLumaQT(pcRecoYuv, tuRecurseWithPU);

					if (pps.getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
					{
						const Int xOffset = tuRecurseWithPU.getRect(COMPONENT_Y).x0;
						const Int yOffset = tuRecurseWithPU.getRect(COMPONENT_Y).y0;
						for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
						{
							if (bMaintainResidual[storedResidualIndex])
							{
								xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
									m_resiPUBuffer[storedResidualIndex],
#else
									resiLumaPU[storedResidualIndex],
#endif
									tuRecurseWithPU, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
							}
						}
					}

					const UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();
					::memcpy(m_puhQTTempTrIdx, pcCU->getTransformIdx() + uiPartOffset, uiQPartNum * sizeof(UChar));
#if COM16_C806_EMT 
					::memcpy(m_puhQTTempEmtTuIdx, pcCU->getEmtTuIdx() + uiPartOffset, uiQPartNum * sizeof(UChar));
					::memcpy(m_puhQTTempEmtCuFlag, pcCU->getEmtCuFlag() + uiPartOffset, uiQPartNum * sizeof(UChar));
#endif

					for (UInt component = 0; component < numberValidComponents; component++)
					{
						const ComponentID compID = ComponentID(component);
						::memcpy(m_puhQTTempCbf[compID], pcCU->getCbf(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
						::memcpy(m_puhQTTempTransformSkipFlag[compID], pcCU->getTransformSkip(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
#if VCEG_AZ08_INTRA_KLT
						::memcpy(m_puhQTTempKLTFlag[compID], pcCU->getKLTFlag(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
#endif
					}
				}
			} // Mode loop
#endif
#endif
			DEBUG_STRING_APPEND(sDebug, sPU)

				//--- update overall distortion ---
				uiOverallDistY += uiBestPUDistY;

			//--- update transform index and cbf ---
			const UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();
#if !JVET_C0024_QTBT
			::memcpy(pcCU->getTransformIdx() + uiPartOffset, m_puhQTTempTrIdx, uiQPartNum * sizeof(UChar));
#endif
#if COM16_C806_EMT
			::memcpy(pcCU->getEmtTuIdx() + uiPartOffset, m_puhQTTempEmtTuIdx, uiQPartNum * sizeof(UChar));
			::memcpy(pcCU->getEmtCuFlag() + uiPartOffset, m_puhQTTempEmtCuFlag, uiQPartNum * sizeof(UChar));
#endif
			for (UInt component = 0; component < numberValidComponents; component++)
			{
				const ComponentID compID = ComponentID(component);
				::memcpy(pcCU->getCbf(compID) + uiPartOffset, m_puhQTTempCbf[compID], uiQPartNum * sizeof(UChar));
				::memcpy(pcCU->getTransformSkip(compID) + uiPartOffset, m_puhQTTempTransformSkipFlag[compID], uiQPartNum * sizeof(UChar));
#if VCEG_AZ08_INTRA_KLT
				::memcpy(pcCU->getKLTFlag(compID) + uiPartOffset, m_puhQTTempKLTFlag[compID], uiQPartNum * sizeof(UChar));
#endif
			}

			//--- set reconstruction for next intra prediction blocks ---
			if (!tuRecurseWithPU.IsLastSection())
			{
#if JVET_C0024_QTBT
				assert(0);
#endif
				const TComRectangle &puRect = tuRecurseWithPU.getRect(COMPONENT_Y);
				const UInt  uiCompWidth = puRect.width;
				const UInt  uiCompHeight = puRect.height;

				const UInt  uiZOrder = pcCU->getZorderIdxInCtu() + uiPartOffset;
				Pel*  piDes = pcCU->getPic()->getPicYuvRec()->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiZOrder);
#if BARC
				Pel*  piDscaleDes = pcCU->getPic()->getPicYuvDScaleRef()->getAddrDscale(COMPONENT_Y, pcCU->getCtuRsAddr(), uiZOrder);
#endif
				const UInt  uiDesStride = pcCU->getPic()->getPicYuvRec()->getStride(COMPONENT_Y);
				const Pel*  piSrc = pcRecoYuv->getAddr(COMPONENT_Y, uiPartOffset);
				const UInt  uiSrcStride = pcRecoYuv->getStride(COMPONENT_Y);

				for (UInt uiY = 0; uiY < uiCompHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride
#if BARC
					, piDscaleDes += uiDesStride
#endif
					)
				{
					for (UInt uiX = 0; uiX < uiCompWidth; uiX++)
					{
						piDes[uiX] = piSrc[uiX];
#if BARC
						piDscaleDes[uiX] = piSrc[uiX];
#endif
					}
				}
			}
#if COM16_C806_EMT
			uiPU++;
#endif
			//=== update PU data ====
			pcCU->setIntraDirSubParts(CHANNEL_TYPE_LUMA, uiBestPUMode, uiPartOffset, uiDepth + uiInitTrDepth);
#if MLIP_EXTEND_REF_LINE
			pcCU->setExternRefSubParts(CHANNEL_TYPE_LUMA, uiBestExternRef, uiPartOffset, uiDepth + uiInitTrDepth);
#endif
#if DIMD_INTRA_PRED
			pcCU->setSecIntraDirSubParts(uiBestPUMode, uiPartOffset, uiDepth + uiInitTrDepth);
#endif
				} while (tuRecurseWithPU.nextSection(tuRecurseCU));

#if !JVET_C0024_QTBT
				if (uiNumPU > 1)
				{ // set Cbf for all blocks
					UInt uiCombCbfY = 0;
					UInt uiCombCbfU = 0;
					UInt uiCombCbfV = 0;
					UInt uiPartIdx = 0;
					for (UInt uiPart = 0; uiPart < 4; uiPart++, uiPartIdx += uiQNumParts)
					{
						uiCombCbfY |= pcCU->getCbf(uiPartIdx, COMPONENT_Y, 1);
						uiCombCbfU |= pcCU->getCbf(uiPartIdx, COMPONENT_Cb, 1);
						uiCombCbfV |= pcCU->getCbf(uiPartIdx, COMPONENT_Cr, 1);
					}
					for (UInt uiOffs = 0; uiOffs < 4 * uiQNumParts; uiOffs++)
					{
						pcCU->getCbf(COMPONENT_Y)[uiOffs] |= uiCombCbfY;
						pcCU->getCbf(COMPONENT_Cb)[uiOffs] |= uiCombCbfU;
						pcCU->getCbf(COMPONENT_Cr)[uiOffs] |= uiCombCbfV;
					}
				}

#endif
				//===== reset context models =====
#if JVET_C0024_QTBT
				m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
#endif

				//===== set distortion (rate and r-d costs are determined later) =====
				pcCU->getTotalDistortion() = uiOverallDistY;
			}




	Void
		TEncSearch::estIntraPredChromaQT(TComDataCU* pcCU,
			TComYuv*    pcOrgYuv,
			TComYuv*    pcPredYuv,
			TComYuv*    pcResiYuv,
			TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
			Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
			Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
			DEBUG_STRING_FN_DECLARE(sDebug))
	{
#if JVET_C0024_QTBT
		const UInt    uiInitTrDepth = 0;
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		const UInt    uiInitTrDepth = pcCU->getPartitionSize(0) != SIZE_2Nx2N && enable4ChromaPUsInIntraNxNCU(pcOrgYuv->getChromaFormat()) ? 1 : 0;
#endif

		TComTURecurse tuRecurseCU(pcCU, 0);
		TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);
#if !JVET_C0024_QTBT
		const UInt    uiQNumParts = tuRecurseWithPU.GetAbsPartIdxNumParts();
#endif
		const UInt    uiDepthCU = tuRecurseWithPU.getCUDepth();
		const UInt    numberValidComponents = pcCU->getPic()->getNumberValidComponents();
#if MLIP_EXTEND_REF_LINE
		UInt ExternRef = pcCU->getExternRef(CHANNEL_TYPE_CHROMA, 0);
		if (ExternRef > 1)
			ExternRef = 1;
#endif
		do
		{

#if JVET_E0077_LM_MF
			UInt auiSATDModeList[LM_FILTER_NUM];
			UInt auiSATDSortedcost[LM_FILTER_NUM];

			const UInt csx = getComponentScaleX(COMPONENT_Cb, pcCU->getSlice()->getSPS()->getChromaFormatIdc());
			const UInt csy = getComponentScaleY(COMPONENT_Cb, pcCU->getSlice()->getSPS()->getChromaFormatIdc());

			Int iBlockSize = (pcCU->getHeight(0) >> csy) + (pcCU->getWidth(0) >> csx);
#endif

			UInt       uiBestMode = 0;
			Distortion uiBestDist = 0;
			Double     dBestCost = MAX_DOUBLE;

			//----- init mode list -----
			if (tuRecurseWithPU.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
			{
#if JVET_E0062_MULTI_DMS
				UInt uiModeList[NUM_CHROMA_MODE];
#else
				UInt uiModeList[FAST_UDI_MAX_RDMODE_NUM];
#endif

#if !JVET_C0024_QTBT
				const UInt  uiQPartNum = uiQNumParts;
#endif
				const UInt  uiPartOffset = tuRecurseWithPU.GetAbsPartIdxTU();
				{
					UInt  uiMinMode = 0;

					UInt  uiMaxMode = NUM_CHROMA_MODE;

					//----- check chroma modes -----
#if JVET_E0077_ENHANCED_LM
					uiMaxMode =
#endif
						pcCU->getAllowedChromaDir(uiPartOffset, uiModeList);

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
					if (DebugOptionList::ForceChromaMode.isSet())
					{
						uiMinMode = DebugOptionList::ForceChromaMode.getInt();
						if (uiModeList[uiMinMode] == 34)
						{
							uiMinMode = 4; // if the fixed mode has been renumbered because DM_CHROMA covers it, use DM_CHROMA.
						}
						uiMaxMode = uiMinMode + 1;
					}
#endif

#if JVET_C0024_QTBT
					assert(uiPartOffset == 0);

					UInt uiShort, uiLong;
					UInt uiStr;

#if TT_TREE_DS1 //data copy - prepare
					//always split into four parts if possible
					UInt uiWidth = pcCU->getWidth(0);
					UInt uiHeight = pcCU->getHeight(0);
					UInt uiNumXDiv = (uiWidth  > pcCU->getPic()->getMinCUWidth()) ? 2 : 1;
					UInt uiNumYDiv = (uiHeight > pcCU->getPic()->getMinCUHeight()) ? 2 : 1;
					uiWidth = (uiNumXDiv == 1) ? uiWidth : (uiWidth >> 1);
					uiHeight = (uiNumYDiv == 1) ? uiHeight : (uiHeight >> 1);
					UInt uiLenX = uiWidth / pcCU->getPic()->getMinCUWidth();
					UInt uiLenY = uiHeight / pcCU->getPic()->getMinCUHeight();

					if (uiHeight > uiWidth)
					{
						uiShort = uiWidth;
						uiLong = uiHeight;
						uiStr = pcCU->getPic()->getNumPartInCtuWidth();
					}
					else
					{
						uiShort = uiHeight;
						uiLong = uiWidth;
						uiStr = 1;
					}
#else
					if (pcCU->getHeight(0) > pcCU->getWidth(0))
					{
						uiShort = pcCU->getWidth(0);
						uiLong = pcCU->getHeight(0);
						uiStr = pcCU->getPic()->getNumPartInCtuWidth();
					}
					else
					{
						uiShort = pcCU->getHeight(0);
						uiLong = pcCU->getWidth(0);
						uiStr = 1;
					}
#endif

					UInt uiCurrDepth = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
					UInt uiCurrPartNumb = pcCU->getPic()->getNumPartitionsInCtu() >> (uiCurrDepth << 1);
					UInt uiNumPartInShort = pcCU->getPic()->getNumPartInCtuWidth() >> uiCurrDepth;
#endif
					DEBUG_STRING_NEW(sPU)

#if JVET_E0077_ENHANCED_LM
						const TComRectangle &puRect = tuRecurseWithPU.getRect(COMPONENT_Cb);
					getLumaRecPixels(tuRecurseWithPU, puRect.width, puRect.height);
					initIntraPatternChType(tuRecurseWithPU, COMPONENT_Cb, false  DEBUG_STRING_PASS_INTO(sDebug));
					initIntraPatternChType(tuRecurseWithPU, COMPONENT_Cr, false  DEBUG_STRING_PASS_INTO(sDebug));

#endif

#if JVET_E0077_LM_MF
					const UInt uiAbsPartIdx = tuRecurseWithPU.GetAbsPartIdxTU();
					DistParam distParam;
					const Bool bUseHadamard = true;
					Int iCurLMMFIdx = 0;

					//SATD checking for LMMF candidates
					const Int iLMMFinRDNum = 1;

					if (pcCU->getSlice()->getSPS()->getUseLMChroma() && iBlockSize >= g_aiMFLM_MinSize[pcCU->getSlice()->isIntra() ? 0 : 1])
					{
						for (UInt uiMode = LM_CHROMA_F1_IDX; uiMode < LM_CHROMA_F1_IDX + LM_FILTER_NUM; uiMode++)
						{
							UInt uiSad = 0;

							Pel* piOrg = pcOrgYuv->getAddr(COMPONENT_Cb, uiAbsPartIdx);
							Pel* piPred = pcPredYuv->getAddr(COMPONENT_Cb, uiAbsPartIdx);
							UInt uiStride = pcPredYuv->getStride(COMPONENT_Cb);
							m_pcRdCost->setDistParam(distParam, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_CHROMA), piOrg, uiStride, piPred, uiStride, puRect.width, puRect.height, bUseHadamard);
							distParam.bApplyWeight = false;

							predLMIntraChroma(tuRecurseWithPU, COMPONENT_Cb, piPred, uiStride, puRect.width, puRect.height, uiMode);

							Pel* piOrgPred = piPred;
							Pel* piSavPred = m_pLMMFPredSaved[(uiMode - LM_CHROMA_F1_IDX) * 2];

							twoDimBufCopier
								(
									(UChar *)piSavPred,
									puRect.width * sizeof(Pel),
									(UChar *)piOrgPred,
									uiStride * sizeof(Pel),
									puRect.width * sizeof(Pel),
									puRect.height
									);

							//for (Int i = 0; i < puRect.height; i++)
							//{
							//	memcpy(piSavPred, piOrgPred, sizeof(Pel)* puRect.width);
							//	piOrgPred += uiStride;
							//	piSavPred += puRect.width;
							//}

							uiSad += distParam.DistFunc(&distParam);


							piOrg = pcOrgYuv->getAddr(COMPONENT_Cr, uiAbsPartIdx);
							piPred = pcPredYuv->getAddr(COMPONENT_Cr, uiAbsPartIdx);
							uiStride = pcPredYuv->getStride(COMPONENT_Cr);
							m_pcRdCost->setDistParam(distParam, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_CHROMA), piOrg, uiStride, piPred, uiStride, puRect.width, puRect.height, bUseHadamard);
							distParam.bApplyWeight = false;

							predLMIntraChroma(tuRecurseWithPU, COMPONENT_Cr, piPred, uiStride, puRect.width, puRect.height, uiMode);

							piOrgPred = piPred;
							piSavPred = m_pLMMFPredSaved[(uiMode - LM_CHROMA_F1_IDX) * 2 + 1];

							twoDimBufCopier
								(
									(UChar *)piSavPred,
									puRect.width * sizeof(Pel),
									(UChar *)piOrgPred,
									uiStride * sizeof(Pel),
									puRect.width * sizeof(Pel),
									puRect.height
									);

							//for (Int i = 0; i < puRect.height; i++)
							//{
							//	memcpy(piSavPred, piOrgPred, sizeof(Pel)* puRect.width);
							//	piOrgPred += uiStride;
							//	piSavPred += puRect.width;
							//}

							uiSad += distParam.DistFunc(&distParam);

							auiSATDSortedcost[iCurLMMFIdx] = uiSad;
							auiSATDModeList[iCurLMMFIdx] = uiMode;
							for (Int k = iCurLMMFIdx; k > 0 && auiSATDSortedcost[k] < auiSATDSortedcost[k - 1]; k--)
							{
								UInt tmp = auiSATDSortedcost[k];
								auiSATDSortedcost[k] = auiSATDSortedcost[k - 1];
								auiSATDSortedcost[k - 1] = tmp;

								tmp = auiSATDModeList[k];
								auiSATDModeList[k] = auiSATDModeList[k - 1];
								auiSATDModeList[k - 1] = tmp;
							}
							iCurLMMFIdx++;
						}
					}
#endif


#if JVET_E0077_ENHANCED_LM
					m_iCurAngMode = 0;
#endif
					for (UInt uiMode = uiMinMode; uiMode < uiMaxMode; uiMode++)
					{

#if COM16_C806_LMCHROMA
						if (!pcCU->getSlice()->getSPS()->getUseLMChroma() && uiModeList[uiMode] == LM_CHROMA_IDX)
						{
							continue;
						}
#endif
#if MLIP_EXTEND_REF_LINE
						if (ExternRef && (uiMode == 0 || uiMode == 3)) {
							continue;
						}
#endif
#if JVET_E0077_ENHANCED_LM
						if (!pcCU->getSlice()->getSPS()->getUseLMChroma() && IsLMMode(uiModeList[uiMode]))
						{
							continue;
						}
#endif

#if JVET_E0077_LM_MF
						if (uiModeList[uiMode] >= LM_CHROMA_F1_IDX &&  uiModeList[uiMode] < LM_CHROMA_F1_IDX + LM_FILTER_NUM)
						{
							Bool bNeedRD = false;
							for (Int i = 0; i < iLMMFinRDNum; i++)
							{
								if (auiSATDModeList[i] == uiModeList[uiMode])
								{
									bNeedRD = true;
									break;
								}
							}
							if (!bNeedRD)
							{
								continue;
							}
						}
#endif

#if JVET_C0024_QTBT && COM16_C1044_NSST
						if (pcCU->getSlice()->isIntra())
						{
							UInt uiIntraMode = uiModeList[uiMode];
#if JVET_E0062_MULTI_DMS && COM16_C806_LMCHROMA            
							if (uiIntraMode == LM_CHROMA_IDX)
							{
								uiIntraMode = PLANAR_IDX;
							}
#if JVET_E0077_ENHANCED_LM
							else if (IsLMMode(uiIntraMode))
							{
								uiIntraMode = PLANAR_IDX;
							}
#endif 
#else
							if (uiIntraMode == DM_CHROMA_IDX)
							{
#if DIMD_INTRA_PRED
								if (pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getDIMDEnabledFlag(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiPartOffset))
								{
									uiIntraMode = DC_IDX + 1;
								}
								else
								{
#endif
									uiIntraMode = pcCU->getPic()->getCtu(pcCU->getCtuRsAddr())->getIntraDir(CHANNEL_TYPE_LUMA, pcCU->getZorderIdxInCtu() + uiPartOffset);
#if DIMD_INTRA_PRED
								}
#endif
							}
#if COM16_C806_LMCHROMA
							else if (uiIntraMode == LM_CHROMA_IDX)
							{
								uiIntraMode = PLANAR_IDX;
							}
#endif

#if JVET_E0077_ENHANCED_LM
							else if (IsLMMode(uiIntraMode))
							{
								uiIntraMode = PLANAR_IDX;
							}
#endif 

#endif
#if FAST_NSST_v1
							const Int iNumberOfPassesROT = (uiIntraMode <= DC_IDX) ? 2 : 4;
#else
							const Int iNumberOfPassesROT = (uiIntraMode <= DC_IDX) ? 3 : 4;
#endif
							if (iNumberOfPassesROT <= pcCU->getROTIdx(CHANNEL_TYPE_CHROMA, uiPartOffset))
							{
								continue;
							}


						}
#endif

#if JVET_E0077_ENHANCED_LM
						if (!IsLMMode(uiModeList[uiMode]))
						{
							m_iCurAngMode++;

#if QC_LM_ANGULAR_PREDICTION
							if (numRDChecking > 1 && uiCPLEFlag == 1 && adModecost[m_iCurAngMode - 1] > adSortedcost[iSecondRDtime - 1])
							{
								continue;
							}
#endif

#if JVET_E0077_ENHANCED_LM && JVET_E0062_MULTI_DMS
							if (m_iCurAngMode > NUM_CHROMA_MODE)
							{
								continue;
							}
#endif
						}
#endif

						//----- restore context models -----
#if JVET_C0024_QTBT
						m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
						m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST]);
#endif

						DEBUG_STRING_NEW(sMode)
							//----- chroma coding -----
							Distortion uiDist = 0;
#if JVET_E0062_MULTI_DMS
						assert(uiModeList[uiMode] <= (LM_CHROMA_IDX + NUM_DM_MODES));
#endif
						pcCU->setIntraDirSubParts(CHANNEL_TYPE_CHROMA, uiModeList[uiMode], uiPartOffset, uiDepthCU + uiInitTrDepth);


						xRecurIntraChromaCodingQT(pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, uiDist, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode));

						if (pcCU->getSlice()->getPPS()->getUseTransformSkip())
						{
#if JVET_C0024_QTBT
							m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
							m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST]);
#endif
						}

						UInt    uiBits = xGetIntraBitsQT(tuRecurseWithPU, false, true, false);
						Double  dCost = m_pcRdCost->calcRdCost(uiBits, uiDist);


						//----- compare -----
						if (dCost < dBestCost)
						{
							DEBUG_STRING_SWAP(sPU, sMode);
							dBestCost = dCost;
							uiBestDist = uiDist;
							uiBestMode = uiModeList[uiMode];
							xSetIntraResultChromaQT(pcRecoYuv, tuRecurseWithPU);
#if JVET_C0024_QTBT            
#if TT_TREE_DS1 //data copy
							for (UInt yIdx = 0; yIdx < uiNumYDiv; yIdx++)
							{
								for (UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
								{
									UInt uiRaster = (uiLenX * xIdx) + (uiLenY * yIdx) * pcCU->getPic()->getNumPartInCtuWidth();
#else
							UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
#endif
							for (UInt i = 0; i<uiLong; i += uiShort)
							{
#if TT_TREE_DS1
								UInt uiZorder = g_auiRasterToZscan[uiRaster];
#else
								UInt uiZorder = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
#endif

								for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
								{
									const ComponentID compID = ComponentID(componentIndex);
									::memcpy(m_puhQTTempCbf[compID] + uiZorder, pcCU->getCbf(compID) + uiZorder, uiCurrPartNumb);
									::memcpy(m_puhQTTempTransformSkipFlag[compID] + uiZorder, pcCU->getTransformSkip(compID) + uiZorder, uiCurrPartNumb);
									::memcpy(m_phQTTempCrossComponentPredictionAlpha[compID] + uiZorder, pcCU->getCrossComponentPredictionAlpha(compID) + uiZorder, uiCurrPartNumb);
								}

								uiRaster += uiNumPartInShort * uiStr;
							}
#if TT_TREE_DS1 //data copy end
								}
							}
#endif
#else
							for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
							{
								const ComponentID compID = ComponentID(componentIndex);
								::memcpy(m_puhQTTempCbf[compID], pcCU->getCbf(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
								::memcpy(m_puhQTTempTransformSkipFlag[compID], pcCU->getTransformSkip(compID) + uiPartOffset, uiQPartNum * sizeof(UChar));
								::memcpy(m_phQTTempCrossComponentPredictionAlpha[compID], pcCU->getCrossComponentPredictionAlpha(compID) + uiPartOffset, uiQPartNum * sizeof(Char));
							}
#endif
						}
					}

					DEBUG_STRING_APPEND(sDebug, sPU)

						//----- set data -----
#if JVET_C0024_QTBT
#if TT_TREE_DS1 //data copy
						for (UInt yIdx = 0; yIdx < uiNumYDiv; yIdx++)
						{
							for (UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
							{
								UInt uiRaster = (uiLenX * xIdx) + (uiLenY * yIdx) * pcCU->getPic()->getNumPartInCtuWidth();
#else
						UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
#endif
					for (UInt i = 0; i<uiLong; i += uiShort)
					{
#if TT_TREE_DS1
						UInt uiZorder = g_auiRasterToZscan[uiRaster];
#else
						UInt uiZorder = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
#endif

						for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
						{
							const ComponentID compID = ComponentID(componentIndex);
							::memcpy(pcCU->getCbf(compID) + uiZorder, m_puhQTTempCbf[compID] + uiZorder, uiCurrPartNumb);
							::memcpy(pcCU->getTransformSkip(compID) + uiZorder, m_puhQTTempTransformSkipFlag[compID] + uiZorder, uiCurrPartNumb);
							::memcpy(pcCU->getCrossComponentPredictionAlpha(compID) + uiZorder, m_phQTTempCrossComponentPredictionAlpha[compID] + uiZorder, uiCurrPartNumb);
						}

						uiRaster += uiNumPartInShort * uiStr;
					}
#if TT_TREE_DS1 //data copy end
							}
						}
#endif
#else
						for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
						{
							const ComponentID compID = ComponentID(componentIndex);
							::memcpy(pcCU->getCbf(compID) + uiPartOffset, m_puhQTTempCbf[compID], uiQPartNum * sizeof(UChar));
							::memcpy(pcCU->getTransformSkip(compID) + uiPartOffset, m_puhQTTempTransformSkipFlag[compID], uiQPartNum * sizeof(UChar));
							::memcpy(pcCU->getCrossComponentPredictionAlpha(compID) + uiPartOffset, m_phQTTempCrossComponentPredictionAlpha[compID], uiQPartNum * sizeof(Char));
						}
#endif
				}

				if (!tuRecurseWithPU.IsLastSection())
				{
#if JVET_C0024_QTBT
					assert(0);
#endif
					for (UInt ch = COMPONENT_Cb; ch<numberValidComponents; ch++)
					{
						const ComponentID compID = ComponentID(ch);
						const TComRectangle &tuRect = tuRecurseWithPU.getRect(compID);
						const UInt  uiCompWidth = tuRect.width;
						const UInt  uiCompHeight = tuRect.height;
						const UInt  uiZOrder = pcCU->getZorderIdxInCtu() + tuRecurseWithPU.GetAbsPartIdxTU();
						Pel*  piDes = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
						const UInt  uiDesStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
						const Pel*  piSrc = pcRecoYuv->getAddr(compID, uiPartOffset);
						const UInt  uiSrcStride = pcRecoYuv->getStride(compID);

						for (UInt uiY = 0; uiY < uiCompHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride)
						{
							for (UInt uiX = 0; uiX < uiCompWidth; uiX++)
							{
								piDes[uiX] = piSrc[uiX];
							}
						}
					}
				}

				pcCU->setIntraDirSubParts(CHANNEL_TYPE_CHROMA, uiBestMode, uiPartOffset, uiDepthCU + uiInitTrDepth);

				pcCU->getTotalDistortion() += uiBestDist;
			}

		} while (tuRecurseWithPU.nextSection(tuRecurseCU));

		//----- restore context models -----
#if !JVET_C0024_QTBT
		if (uiInitTrDepth != 0)
		{ // set Cbf for all blocks
			UInt uiCombCbfU = 0;
			UInt uiCombCbfV = 0;
			UInt uiPartIdx = 0;
			for (UInt uiPart = 0; uiPart < 4; uiPart++, uiPartIdx += uiQNumParts)
			{
				uiCombCbfU |= pcCU->getCbf(uiPartIdx, COMPONENT_Cb, 1);
				uiCombCbfV |= pcCU->getCbf(uiPartIdx, COMPONENT_Cr, 1);
			}
			for (UInt uiOffs = 0; uiOffs < 4 * uiQNumParts; uiOffs++)
			{
				pcCU->getCbf(COMPONENT_Cb)[uiOffs] |= uiCombCbfU;
				pcCU->getCbf(COMPONENT_Cr)[uiOffs] |= uiCombCbfV;
			}
		}
#endif

#if JVET_C0024_QTBT
		m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
		m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST]);
#endif
	}

#if DIMD_INTRA_PRED
	Void  TEncSearch::estDIMDIntraPredLumaQT(TComDataCU* pcCU,
		TComYuv*    pcOrgYuv,
		TComYuv*    pcPredYuv,
		TComYuv*    pcResiYuv,
		TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
		Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
		Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
	)
	{
		const UInt    uiInitTrDepth = 0;
		const UInt    uiDepth = pcCU->getDepth(0);
		const UInt    uiWidth = pcCU->getWidth(0);
		const UInt    uiHeight = pcCU->getHeight(0);
		const UInt    uiWIdx = g_aucConvertToBit[uiWidth];
		const UInt    uiHIdx = g_aucConvertToBit[uiHeight];

		const TComPPS     &pps = *(pcCU->getSlice()->getPPS());
		Distortion        uiOverallDistY = 0;
#if BARC_DIMD
		Bool     isBARC = pcCU->getBARCFlag(CHANNEL_TYPE_LUMA);
#endif
#if !COM16_C806_LARGE_CTU
		Pel          resiLumaPU[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
#endif

		Bool    bMaintainResidual[NUMBER_OF_STORED_RESIDUAL_TYPES];
		for (UInt residualTypeIndex = 0; residualTypeIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; residualTypeIndex++)
		{
			bMaintainResidual[residualTypeIndex] = true; //assume true unless specified otherwise
		}

		bMaintainResidual[RESIDUAL_ENCODER_SIDE] = !(m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate());

		//===== set QP and clear Cbf =====
		if (pps.getUseDQP() == true)
		{
#if JVET_C0024_DELTA_QP_FIX
			pcCU->setQPSubParts(pcCU->getQP(0), 0, uiWidth, uiHeight);
#else
			pcCU->setQPSubParts(pcCU->getQP(0), 0, uiDepth);
#endif
		}
		else
		{
#if JVET_C0024_DELTA_QP_FIX
#if BARC_DIMD
			pcCU->setQPSubParts(pcCU->getQP(0), 0, uiWidth, uiHeight);
#else
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), 0, uiWidth, uiHeight);
#endif
#else
			pcCU->setQPSubParts(pcCU->getSlice()->getSliceQp(), 0, uiDepth);
#endif
		}

		//===== loop over partitions =====
		TComTURecurse tuRecurseCU(pcCU, 0);
		TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

		do
		{
			const UInt uiPartOffset = tuRecurseWithPU.GetAbsPartIdxTU();
			assert(uiPartOffset == 0);
#if !DIMD_INTRA_PRED
#if COM16_C1044_NSST
			const Int iNumberOfPassesROT = (pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiPartOffset) <= DC_IDX) ? 3 : 4;
			assert(pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiPartOffset) < iNumberOfPassesROT);
#endif
#endif

			// set context models
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);

			// determine residual for partition
			Distortion uiPUDistY = 0;
			Double     dPUCost = 0.0;

#if HHI_RQT_INTRA_SPEEDUP
			xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
				m_resiPUBuffer,
#else
				resiLumaPU,
#endif
				uiPUDistY, true, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode));
#else
			xRecurIntraCodingLumaQT(pcOrgYuv, pcPredYuv, pcResiYuv,
#if COM16_C806_LARGE_CTU
				m_resiPUBuffer,
#else
				resiLumaPU,
#endif
				uiPUDistY, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode));
#endif

			xSetIntraResultLumaQT(pcRecoYuv, tuRecurseWithPU);

			if (pps.getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
			{
				const Int xOffset = tuRecurseWithPU.getRect(COMPONENT_Y).x0;
				const Int yOffset = tuRecurseWithPU.getRect(COMPONENT_Y).y0;
				for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
				{
					if (bMaintainResidual[storedResidualIndex])
					{
						xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex],
#if COM16_C806_LARGE_CTU
							m_resiPUBuffer[storedResidualIndex],
#else
							resiLumaPU[storedResidualIndex],
#endif
							tuRecurseWithPU, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
					}
				}
			}

			//--- update overall distortion ---
			uiOverallDistY += uiPUDistY;

			assert(tuRecurseWithPU.IsLastSection());
		} while (tuRecurseWithPU.nextSection(tuRecurseCU));

		//===== reset context models =====
		m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);

		//===== set distortion (rate and r-d costs are determined later) =====
		pcCU->getTotalDistortion() = uiOverallDistY;
	}

	Void TEncSearch::estDIMDIntraPredChromaQT(TComDataCU* pcCU,
		TComYuv*    pcOrgYuv,
		TComYuv*    pcPredYuv,
		TComYuv*    pcResiYuv,
		TComYuv*    pcRecoYuv,
#if COM16_C806_LARGE_CTU
		Pel*        resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES]
#else
		Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
#endif
	)
	{
		const UInt    uiInitTrDepth = 0;
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];

		TComTURecurse tuRecurseCU(pcCU, 0);
		TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

		do
		{
			//----- init mode list -----
			if (tuRecurseWithPU.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
			{
				const UInt  uiPartOffset = tuRecurseWithPU.GetAbsPartIdxTU();
				assert(uiPartOffset == 0);
#if COM16_C1044_NSST
				UInt uiIntraMode = pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiPartOffset);
				assert(uiIntraMode == DM_CHROMA_IDX);
				assert(!pcCU->getSlice()->isIntra());
#endif

				//----- restore context models -----
				m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);

				//----- chroma coding -----
				Distortion uiDist = 0;
				xRecurIntraChromaCodingQT(pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, uiDist, tuRecurseWithPU);
				xSetIntraResultChromaQT(pcRecoYuv, tuRecurseWithPU);

				assert(tuRecurseWithPU.IsLastSection());
				pcCU->getTotalDistortion() += uiDist;
			}
		} while (tuRecurseWithPU.nextSection(tuRecurseCU));

		m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
	}
#endif


	/** Function for encoding and reconstructing luma/chroma samples of a PCM mode CU.
	* \param pcCU pointer to current CU
	* \param uiAbsPartIdx part index
	* \param pOrg pointer to original sample arrays
	* \param pPCM pointer to PCM code arrays
	* \param pPred pointer to prediction signal arrays
	* \param pResi pointer to residual signal arrays
	* \param pReco pointer to reconstructed sample arrays
	* \param uiStride stride of the original/prediction/residual sample arrays
	* \param uiWidth block width
	* \param uiHeight block height
	* \param compID texture component type
	*/
	Void TEncSearch::xEncPCM(TComDataCU* pcCU, UInt uiAbsPartIdx, Pel* pOrg, Pel* pPCM, Pel* pPred, Pel* pResi, Pel* pReco, UInt uiStride, UInt uiWidth, UInt uiHeight, const ComponentID compID)
	{
		const UInt uiReconStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
		const UInt uiPCMBitDepth = pcCU->getSlice()->getSPS()->getPCMBitDepth(toChannelType(compID));
		const Int  channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
		Pel* pRecoPic = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx);

		const Int pcmShiftRight = (channelBitDepth - Int(uiPCMBitDepth));

		assert(pcmShiftRight >= 0);

		for (UInt uiY = 0; uiY < uiHeight; uiY++)
		{
			for (UInt uiX = 0; uiX < uiWidth; uiX++)
			{
				// Reset pred and residual
				pPred[uiX] = 0;
				pResi[uiX] = 0;
				// Encode
				pPCM[uiX] = (pOrg[uiX] >> pcmShiftRight);
				// Reconstruction
				pReco[uiX] = (pPCM[uiX] << (pcmShiftRight));
				pRecoPic[uiX] = pReco[uiX];
			}
			pPred += uiStride;
			pResi += uiStride;
			pPCM += uiWidth;
			pOrg += uiStride;
			pReco += uiStride;
			pRecoPic += uiReconStride;
		}
	}


	//!  Function for PCM mode estimation.
	Void TEncSearch::IPCMSearch(TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv)
	{
		UInt        uiDepth = pcCU->getDepth(0);
		const UInt  uiDistortion = 0;
		UInt        uiBits;

		Double dCost;

		for (UInt ch = 0; ch < pcCU->getPic()->getNumberValidComponents(); ch++)
		{
			const ComponentID compID = ComponentID(ch);
			const UInt width = pcCU->getWidth(0) >> pcCU->getPic()->getComponentScaleX(compID);
			const UInt height = pcCU->getHeight(0) >> pcCU->getPic()->getComponentScaleY(compID);
			const UInt stride = pcPredYuv->getStride(compID);

			Pel * pOrig = pcOrgYuv->getAddr(compID, 0, width);
			Pel * pResi = pcResiYuv->getAddr(compID, 0, width);
			Pel * pPred = pcPredYuv->getAddr(compID, 0, width);
			Pel * pReco = pcRecoYuv->getAddr(compID, 0, width);
			Pel * pPCM = pcCU->getPCMSample(compID);

			xEncPCM(pcCU, 0, pOrig, pPCM, pPred, pResi, pReco, stride, width, height, compID);

		}

		m_pcEntropyCoder->resetBits();
		xEncIntraHeader(pcCU, uiDepth, 0, true, false);
		uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();

		dCost = m_pcRdCost->calcRdCost(uiBits, uiDistortion);

#if JVET_C0024_QTBT
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
		m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
		m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
#endif

		pcCU->getTotalBits() = uiBits;
		pcCU->getTotalCost() = dCost;
		pcCU->getTotalDistortion() = uiDistortion;

#if JVET_C0024_QTBT
		pcCU->copyToPic(uiDepth, pcCU->getWidth(0), pcCU->getHeight(0));
#else
		pcCU->copyToPic(uiDepth);
#endif
	}




	Void TEncSearch::xGetInterPredictionError(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, Distortion& ruiErr, Bool /*bHadamard*/)
	{
#if JVET_E0052_DMVR
		motionCompensation(pcCU, &m_tmpYuvPred, true, REF_PIC_LIST_X, iPartIdx);
#else
		motionCompensation(pcCU, &m_tmpYuvPred, REF_PIC_LIST_X, iPartIdx);
#endif

		UInt uiAbsPartIdx = 0;
		Int iWidth = 0;
		Int iHeight = 0;
		pcCU->getPartIndexAndSize(iPartIdx, uiAbsPartIdx, iWidth, iHeight);

		DistParam cDistParam;

		cDistParam.bApplyWeight = false;


		m_pcRdCost->setDistParam(cDistParam, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA),
			pcYuvOrg->getAddr(COMPONENT_Y, uiAbsPartIdx), pcYuvOrg->getStride(COMPONENT_Y),
			m_tmpYuvPred.getAddr(COMPONENT_Y, uiAbsPartIdx), m_tmpYuvPred.getStride(COMPONENT_Y),
			iWidth, iHeight, m_pcEncCfg->getUseHADME() && (pcCU->getCUTransquantBypass(iPartIdx) == 0));

		ruiErr = cDistParam.DistFunc(&cDistParam);
	}

	//! estimation of best merge coding
	Void TEncSearch::xMergeEstimation(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPUIdx, UInt& uiInterDir, TComMvField* pacMvField, UInt& uiMergeIndex, Distortion& ruiCost, TComMvField* cMvFieldNeighbours, UChar* uhInterDirNeighbours, Int& numValidMergeCand
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
		, UChar*          pMergeTypeNeighbor
#if JVET_C0035_ATMVP_SIMPLIFICATION
		, TComMvField*    pcMvFieldSP[NUM_MGR_TYPE]
		, UChar*          puhInterDirSP[NUM_MGR_TYPE]
#else
		, TComMvField*    pcMvFieldSP[2]
		, UChar*          puhInterDirSP[2]
#endif
#endif
	)
	{
		UInt uiAbsPartIdx = 0;
		Int iWidth = 0;
		Int iHeight = 0;
#if VCEG_AZ06_IC
		Bool abICFlag[MRG_MAX_NUM_CANDS];
#endif
		pcCU->getPartIndexAndSize(iPUIdx, uiAbsPartIdx, iWidth, iHeight);
#if JVET_C0024_QTBT
		PartSize partSize = SIZE_2Nx2N;
#else
		UInt uiDepth = pcCU->getDepth(uiAbsPartIdx);

		PartSize partSize = pcCU->getPartitionSize(0);
#endif
		if (pcCU->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() && partSize != SIZE_2Nx2N && pcCU->getWidth(0) <= 8)
		{
			if (iPUIdx == 0)
			{
#if !JVET_C0024_QTBT
				pcCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uiDepth); // temporarily set
#endif
				pcCU->getInterMergeCandidates(0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#if VCEG_AZ06_IC
					, abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
					, pMergeTypeNeighbor
					, pcMvFieldSP
					, puhInterDirSP
#endif
				);
#if !JVET_C0024_QTBT
				pcCU->setPartSizeSubParts(partSize, 0, uiDepth); // restore
#endif
			}
		}
		else
		{
			pcCU->getInterMergeCandidates(uiAbsPartIdx, iPUIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#if VCEG_AZ06_IC
				, abICFlag
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
				, pMergeTypeNeighbor
				, pcMvFieldSP
				, puhInterDirSP
#endif
			);
		}

		xRestrictBipredMergeCand(pcCU, iPUIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand);

		ruiCost = std::numeric_limits<Distortion>::max();
		for (UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand)
		{
			Distortion uiCostCand = std::numeric_limits<Distortion>::max();
			UInt       uiBitsCand = 0;

#if JVET_C0024_QTBT
			PartSize ePartSize = SIZE_2Nx2N;
#else
			PartSize ePartSize = pcCU->getPartitionSize(0);
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
			pcCU->setMergeTypeSubParts(pMergeTypeNeighbor[uiMergeCand], uiAbsPartIdx, iPUIdx, pcCU->getDepth(uiAbsPartIdx));
			if (pMergeTypeNeighbor[uiMergeCand])
			{

				UInt uiSPAddr;
				Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
#if JVET_C0035_ATMVP_SIMPLIFICATION
				UInt uiSPListIndex = pMergeTypeNeighbor[uiMergeCand];
#else
				UInt uiSPListIndex = pMergeTypeNeighbor[uiMergeCand] == MGR_TYPE_SUBPU_TMVP ? 0 : 1;
#endif
				pcCU->getSPPara(iWidth, iHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);

				for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
				{
					pcCU->getSPAbsPartIdx(uiAbsPartIdx, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
					pcCU->getCUMvField(REF_PIC_LIST_0)->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx], iSPWidth, iSPHeight);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx + 1], iSPWidth, iSPHeight);
				}
			}
			else
			{
#endif
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(cMvFieldNeighbours[0 + 2 * uiMergeCand], ePartSize, uiAbsPartIdx, 0, iPUIdx);
				pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(cMvFieldNeighbours[1 + 2 * uiMergeCand], ePartSize, uiAbsPartIdx, 0, iPUIdx);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
			}
#endif

			xGetInterPredictionError(pcCU, pcYuvOrg, iPUIdx, uiCostCand, m_pcEncCfg->getUseHADME());
#if JVET_D0123_ME_CTX_LUT_BITS
			uiBitsCand = xGetMrgIdxBits(pcCU, uiAbsPartIdx, uiMergeCand);
			uiCostCand = uiCostCand + m_pcRdCost->getCostSearch(uiBitsCand);
#else
			uiBitsCand = uiMergeCand + 1;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
			UInt uiMaxCand = m_pcEncCfg->getMaxNumMergeCand() + (m_pcEncCfg->getAtmvp() ? 2 : 0);
			if (uiMergeCand == uiMaxCand - 1)
#else
			if (uiMergeCand == m_pcEncCfg->getMaxNumMergeCand() - 1)
#endif
			{
				uiBitsCand--;
			}
			uiCostCand = uiCostCand + m_pcRdCost->getCost(uiBitsCand);
#endif
			if (uiCostCand < ruiCost)
			{
				ruiCost = uiCostCand;
				pacMvField[0] = cMvFieldNeighbours[0 + 2 * uiMergeCand];
				pacMvField[1] = cMvFieldNeighbours[1 + 2 * uiMergeCand];
				uiInterDir = uhInterDirNeighbours[uiMergeCand];
				uiMergeIndex = uiMergeCand;
			}
		}
	}

#if VCEG_AZ07_FRUC_MERGE
	Void TEncSearch::xFRUCMgrEstimation(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPUIdx, TComMvField* pacMvField,
		UChar * phInterDir, UChar ** phFRUCRefineDist, UChar ** phFRUCSBlkRefineDist,
		UInt& ruiMinCost, UChar & ruhFRUCMode)
	{
		UInt uiAbsPartIdx = 0;
		Int iWidth = 0;
		Int iHeight = 0;
		pcCU->getPartIndexAndSize(iPUIdx, uiAbsPartIdx, iWidth, iHeight);
		UInt uiDepth = pcCU->getDepth(uiAbsPartIdx);
		ruiMinCost = MAX_UINT;

		const UChar uhFRUCME[2] = { FRUC_MERGE_BILATERALMV , FRUC_MERGE_TEMPLATE };
		pcCU->setMergeFlagSubParts(true, uiAbsPartIdx, iPUIdx, uiDepth);
		for (Int nME = 0; nME < 2; nME++)
		{
			pcCU->setFRUCMgrModeSubParts(uhFRUCME[nME], uiAbsPartIdx, iPUIdx, uiDepth);
			Bool bAvailable = deriveFRUCMV(pcCU, uiDepth, uiAbsPartIdx, iPUIdx);
			if (bAvailable)
			{
				UInt uiCostCand = 0;
				xGetInterPredictionError(pcCU, pcYuvOrg, iPUIdx, uiCostCand, m_pcEncCfg->getUseHADME());
#if JVET_D0123_ME_CTX_LUT_BITS
				UInt uiBitsCand = xGetFRUCMgrModeBits(pcCU, uiAbsPartIdx, iPUIdx, uhFRUCME[nME]);
				UInt uiCost = uiCostCand + m_pcRdCost->getCostSearch(uiBitsCand);
#else
				UInt uiBitsCand = 1;
				UInt uiCost = uiCostCand + m_pcRdCost->getCost(uiBitsCand);
#endif
				if (uiCost < ruiMinCost)
				{
					ruiMinCost = uiCost;
					ruhFRUCMode = uhFRUCME[nME];
					for (Int y = 0, yRasterOffset = 0; y < iHeight; y += 4, yRasterOffset += pcCU->getPic()->getNumPartInCtuWidth())
					{
						for (Int x = 0, xRasterOffset = 0; x < iWidth; x += 4, xRasterOffset++)
						{
							UInt idx = g_auiRasterToZscan[g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiAbsPartIdx] + yRasterOffset + xRasterOffset] - pcCU->getZorderIdxInCtu();
							pacMvField[idx << 1].setMvField(pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(idx), pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(idx));
							pacMvField[(idx << 1) + 1].setMvField(pcCU->getCUMvField(REF_PIC_LIST_1)->getMv(idx), pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(idx));
							phInterDir[idx] = pcCU->getInterDir(idx);
						}
					}
				}
			}
		}
	}
#endif

	/** convert bi-pred merge candidates to uni-pred
	* \param pcCU
	* \param puIdx
	* \param mvFieldNeighbours
	* \param interDirNeighbours
	* \param numValidMergeCand
	* \returns Void
	*/
	Void TEncSearch::xRestrictBipredMergeCand(TComDataCU* pcCU, UInt puIdx, TComMvField* mvFieldNeighbours, UChar* interDirNeighbours, Int numValidMergeCand)
	{
		if (pcCU->isBipredRestriction(puIdx))
		{
			for (UInt mergeCand = 0; mergeCand < numValidMergeCand; ++mergeCand)
			{
				if (interDirNeighbours[mergeCand] == 3)
				{
					interDirNeighbours[mergeCand] = 1;
					mvFieldNeighbours[(mergeCand << 1) + 1].setMvField(TComMv(0, 0), -1);
				}
			}
		}
	}

	//! search of the best candidate for inter prediction
#if AMP_MRG
#if COM16_C806_OBMC
#if AFF_EARLY_TERMINATION
Void TEncSearch::predInterSearch(TComDataCU* pcCU, TComDataCU* pcBestCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv DEBUG_STRING_FN_DECLARE(sDebug), TComYuv *pcPredYuvWoOBMC, TComYuv *pcTmpYuv1, TComYuv *pcTmpYuv2, Bool bUseRes, Bool bUseMRG)
#else
Void TEncSearch::predInterSearch( TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv DEBUG_STRING_FN_DECLARE(sDebug), TComYuv *pcPredYuvWoOBMC, TComYuv *pcTmpYuv1, TComYuv *pcTmpYuv2, Bool bUseRes, Bool bUseMRG )
#endif
#else
	Void TEncSearch::predInterSearch(TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseRes, Bool bUseMRG)
#endif
#else
#if COM16_C806_OBMC
	Void TEncSearch::predInterSearch(TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv, TComYuv *pcPredYuvWoOBMC, TComYuv *pcTmpYuv1, TComYuv *pcTmpYuv2, Bool bUseRes)
#else
	Void TEncSearch::predInterSearch(TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv, Bool bUseRes)
#endif
#endif
	{
#if !COM16_C806_LARGE_CTU
		for (UInt i = 0; i<NUM_REF_PIC_LIST_01; i++)
		{
			m_acYuvPred[i].clear();
		}
		m_cYuvPredTemp.clear();
		pcPredYuv->clear();
#endif

		if (!bUseRes)
		{
			pcResiYuv->clear();
		}

#if !COM16_C806_LARGE_CTU
		pcRecoYuv->clear();
#endif

		TComMv       cMvSrchRngLT;
		TComMv       cMvSrchRngRB;

		TComMv       cMvZero;
		TComMv       TempMv; //kolya

		TComMv       cMv[2];
		TComMv       cMvBi[2];
		TComMv       cMvTemp[2][33];
#if COM16_C1016_AFFINE
		TComMv       cMvHevcTemp[2][33];
#endif

		Int          iNumPart = pcCU->getNumPartitions();
		Int          iNumPredDir = pcCU->getSlice()->isInterP() ? 1 : 2;

		TComMv       cMvPred[2][33];

		TComMv       cMvPredBi[2][33];
		Int          aaiMvpIdxBi[2][33];

		Int          aaiMvpIdx[2][33];
		Int          aaiMvpNum[2][33];

		AMVPInfo     aacAMVPInfo[2][33];

		Int          iRefIdx[2] = { 0,0 }; //If un-initialized, may cause SEGV in bi-directional prediction iterative stage.
		Int          iRefIdxBi[2];

		UInt         uiPartAddr;
		Int          iRoiWidth, iRoiHeight;

		UInt         uiMbBits[3] = { 1, 1, 0 };

		UInt         uiLastMode = 0;
#if COM16_C1016_AFFINE
		UInt         uiLastModeTemp = 0;
#endif
		Int          iRefStart, iRefEnd;

#if JVET_C0024_QTBT
		PartSize     ePartSize = SIZE_2Nx2N;
#else
		PartSize     ePartSize = pcCU->getPartitionSize(0);
#endif

		Int          bestBiPRefIdxL1 = 0;
		Int          bestBiPMvpL1 = 0;
		Distortion   biPDistTemp = std::numeric_limits<Distortion>::max();

#if !JVET_C0024_QTBT
		TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
		UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
		Int numValidMergeCand = 0;
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
		UChar eMergeCandTypeNieghors[MRG_MAX_NUM_CANDS];
		for (Int i = 0; i< MRG_MAX_NUM_CANDS; i++)
		{
			eMergeCandTypeNieghors[i] = MGR_TYPE_DEFAULT_N;
		}
#endif
#endif

#if JVET_D0123_ME_CTX_LUT_BITS
		loadCtxMe(pcCU);
#endif

		for (Int iPartIdx = 0; iPartIdx < iNumPart; iPartIdx++)
		{
#if COM16_C1016_AFFINE
			Distortion  uiHevcCost = std::numeric_limits<Distortion>::max();
			Distortion  uiAffineCost = std::numeric_limits<Distortion>::max();
#endif

#if COM16_C806_OBMC
			//consider OBMC in motion estimation
#if JVET_C0024_QTBT
			assert(iPartIdx == 0);
#else
			pcOrgYuv->copyFromPicYuv(pcCU->getPic()->getPicYuvOrg(), pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
#endif
#if JVET_E0052_DMVR
			subBlockOBMC(pcCU, 0, pcOrgYuv, pcTmpYuv1, pcTmpYuv2, true, true);
#else
			subBlockOBMC(pcCU, 0, pcOrgYuv, pcTmpYuv1, pcTmpYuv2, true);
#endif
#endif
			Distortion   uiCost[2] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };
			Distortion   uiCostBi = std::numeric_limits<Distortion>::max();
			Distortion   uiCostTemp;

			UInt         uiBits[3];
			UInt         uiBitsTemp;
			Distortion   bestBiPDist = std::numeric_limits<Distortion>::max();

			Distortion   uiCostTempL0[MAX_NUM_REF];
			for (Int iNumRef = 0; iNumRef < MAX_NUM_REF; iNumRef++)
			{
				uiCostTempL0[iNumRef] = std::numeric_limits<Distortion>::max();
			}
			UInt         uiBitsTempL0[MAX_NUM_REF];

			TComMv       mvValidList1;
			Int          refIdxValidList1 = 0;
			UInt         bitsValidList1 = MAX_UINT;
			Distortion   costValidList1 = std::numeric_limits<Distortion>::max();
#if JVET_D0123_ME_CTX_LUT_BITS
			//m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_PU_NEXT_BEST]);
			loadCtxMeSbacCoder(pcCU);
			//m_pcEntropyCoder->resetBits();

			m_pcEntropyCoder->estimatePuMeBit(m_pcPuMeEstBitsSbac);

			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
			xGetBlkBits(ePartSize, pcCU->getSlice()->isInterP(), iPartIdx, uiLastMode, uiMbBits, pcCU, uiPartAddr);
			xUpdateMvpIdxBits(pcCU);
#if VCEG_AZ07_IMV
			if (iPartIdx == 0)
				xGetIMVFlagBits(pcCU, uiPartAddr);
#endif
#else
			xGetBlkBits(ePartSize, pcCU->getSlice()->isInterP(), iPartIdx, uiLastMode, uiMbBits);
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO 
			UInt uiWidthIdx = g_aucConvertToBit[iRoiWidth];
			UInt uiHeightIdx = g_aucConvertToBit[iRoiHeight];
			UInt uiZorderIdx = pcCU->getZorderIdxInCtu();
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
			UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, pcCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx);
#if COM16_C1016_AFFINE
			Bool bFastSkipAffine = (saveLoadTag == LOAD_ENC_INFO && !getSaveLoadAffineFlag(uiZorderIdx, pcCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx));
#endif
			Bool bFastSkipBi = (saveLoadTag == LOAD_ENC_INFO && getSaveLoadInterDir(uiZorderIdx, pcCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx) != 3);
#else
			UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, uiWidthIdx, uiHeightIdx);
#if COM16_C1016_AFFINE
			Bool bFastSkipAffine = (saveLoadTag == LOAD_ENC_INFO && !getSaveLoadAffineFlag(uiWidthIdx, uiHeightIdx));
#endif
			Bool bFastSkipBi = (saveLoadTag == LOAD_ENC_INFO && getSaveLoadInterDir(uiWidthIdx, uiHeightIdx) != 3);
#endif
#endif

#if AMP_MRG
			Bool bTestNormalMC = true;

			if (bUseMRG && pcCU->getWidth(0) > 8 && iNumPart == 2)
			{
				bTestNormalMC = false;
			}

			if (bTestNormalMC)
			{
#endif

				//  Uni-directional prediction
				for (Int iRefList = 0; iRefList < iNumPredDir; iRefList++)
				{
					RefPicList  eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

					for (Int iRefIdxTemp = 0; iRefIdxTemp < pcCU->getSlice()->getNumRefIdx(eRefPicList); iRefIdxTemp++)
					{
						uiBitsTemp = uiMbBits[iRefList];
						if (pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1)
						{
#if JVET_D0123_ME_CTX_LUT_BITS
							uiBitsTemp += xRefFrameIdxBits(pcCU, iRefIdxTemp, (UInt)pcCU->getSlice()->getNumRefIdx(eRefPicList));
#else
							uiBitsTemp += iRefIdxTemp + 1;
							if (iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1)
							{
								uiBitsTemp--;
							}
#endif
						}
						xEstimateMvPredAMVP(pcCU, pcOrgYuv, iPartIdx, eRefPicList, iRefIdxTemp, cMvPred[iRefList][iRefIdxTemp], false, &biPDistTemp);
						aaiMvpIdx[iRefList][iRefIdxTemp] = pcCU->getMVPIdx(eRefPicList, uiPartAddr);
						aaiMvpNum[iRefList][iRefIdxTemp] = pcCU->getMVPNum(eRefPicList, uiPartAddr);

						if (pcCU->getSlice()->getMvdL1ZeroFlag() && iRefList == 1 && biPDistTemp < bestBiPDist)
						{
							bestBiPDist = biPDistTemp;
							bestBiPMvpL1 = aaiMvpIdx[iRefList][iRefIdxTemp];
							bestBiPRefIdxL1 = iRefIdxTemp;
						}

						uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];

						if (m_pcEncCfg->getFastMEForGenBLowDelayEnabled() && iRefList == 1)    // list 1
						{
							if (pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp) >= 0)
							{
								cMvTemp[1][iRefIdxTemp] = cMvTemp[0][pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp)];
								uiCostTemp = uiCostTempL0[pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp)];
								/*first subtract the bit-rate part of the cost of the other list*/
#if JVET_D0123_ME_CTX_LUT_BITS
								uiCostTemp -= m_pcRdCost->getCostSearch(uiBitsTempL0[pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp)]);
#else
								uiCostTemp -= m_pcRdCost->getCost(uiBitsTempL0[pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp)]);
#endif
								/*correct the bit-rate part of the current ref*/
								m_pcRdCost->setPredictor(cMvPred[iRefList][iRefIdxTemp]);
#if JVET_D0123_ME_CTX_LUT_BITS
								TComPattern   tmpPattern;
								TComPattern*  pcPatternKey = &tmpPattern;
								pcPatternKey->m_pcCU = pcCU;
								pcPatternKey->m_ePartSize = ePartSize;
								pcPatternKey->m_uiPartAddr = uiPartAddr;
								pcPatternKey->m_mvPred = &cMvPred[iRefList][iRefIdxTemp];

								iCostScale = 0;

								uiBitsTemp += xMvdBits(pcPatternKey, cMvTemp[1][iRefIdxTemp].getHor(), cMvTemp[1][iRefIdxTemp].getVer());
								uiCostTemp += m_pcRdCost->getCostSearch(uiBitsTemp);
#else
								uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
									cMvTemp[1][iRefIdxTemp].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, cMvTemp[1][iRefIdxTemp].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
									cMvTemp[1][iRefIdxTemp].getHor(), cMvTemp[1][iRefIdxTemp].getVer()
#endif
#if VCEG_AZ07_IMV
									, pcCU->getiMVFlag(uiPartAddr)
#endif
								);
								/*calculate the correct cost*/
								uiCostTemp += m_pcRdCost->getCost(uiBitsTemp);
#endif
							}
							else
							{
#if JVET_E0076_MULTI_PEL_MVD
								xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);
#else
								xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);
#endif
							}
						}
						else
						{
#if JVET_E0076_MULTI_PEL_MVD
							xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);
#else
							xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);
#endif
						}
						xCopyAMVPInfo(pcCU->getCUMvField(eRefPicList)->getAMVPInfo(), &aacAMVPInfo[iRefList][iRefIdxTemp]); // must always be done ( also when AMVP_MODE = AM_NONE )
						xCheckBestMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
							, uiPartAddr
#endif
						);

						if (iRefList == 0)
						{
							uiCostTempL0[iRefIdxTemp] = uiCostTemp;
							uiBitsTempL0[iRefIdxTemp] = uiBitsTemp;
						}
						if (uiCostTemp < uiCost[iRefList])
						{
							uiCost[iRefList] = uiCostTemp;
							uiBits[iRefList] = uiBitsTemp; // storing for bi-prediction

														   // set motion
							cMv[iRefList] = cMvTemp[iRefList][iRefIdxTemp];
							iRefIdx[iRefList] = iRefIdxTemp;
						}

						if (iRefList == 1 && uiCostTemp < costValidList1 && pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp) < 0)
						{
							costValidList1 = uiCostTemp;
							bitsValidList1 = uiBitsTemp;

							// set motion
							mvValidList1 = cMvTemp[iRefList][iRefIdxTemp];
							refIdxValidList1 = iRefIdxTemp;
						}
					}
				}

#if COM16_C1016_AFFINE // save regular Hevc ME result for Affine ME
				::memcpy(cMvHevcTemp, cMvTemp, sizeof(cMvTemp));
#endif

				//  Bi-directional prediction
#if COM16_C806_LARGE_CTU && !JVET_C0024_QTBT
				UChar uiMaxCUDepth = 0, uiMinCUDepth = 0;
				pcCU->getMaxMinCUDepth(uiMinCUDepth, uiMaxCUDepth, pcCU->getZorderIdxInCtu());
				Bool bCheckBi = true;
				if (m_pcEncCfg->getUseFastLCTU())
				{
					bCheckBi = pcCU->getDepth(0) < uiMaxCUDepth || iRoiWidth * iRoiHeight >= 64;
				}
				if (bCheckBi)
				{
#endif
					if ((pcCU->getSlice()->isInterB()) && (pcCU->isBipredRestriction(iPartIdx) == false)
#if JVET_D0077_SAVE_LOAD_ENC_INFO
						&& !bFastSkipBi
#endif
						)
					{

						cMvBi[0] = cMv[0];            cMvBi[1] = cMv[1];
						iRefIdxBi[0] = iRefIdx[0];    iRefIdxBi[1] = iRefIdx[1];

						::memcpy(cMvPredBi, cMvPred, sizeof(cMvPred));
						::memcpy(aaiMvpIdxBi, aaiMvpIdx, sizeof(aaiMvpIdx));

						UInt uiMotBits[2];

						if (pcCU->getSlice()->getMvdL1ZeroFlag())
						{
							xCopyAMVPInfo(&aacAMVPInfo[1][bestBiPRefIdxL1], pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo());
							pcCU->setMVPIdxSubParts(bestBiPMvpL1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
							aaiMvpIdxBi[1][bestBiPRefIdxL1] = bestBiPMvpL1;
							cMvPredBi[1][bestBiPRefIdxL1] = pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo()->m_acMvCand[bestBiPMvpL1];

							cMvBi[1] = cMvPredBi[1][bestBiPRefIdxL1];
							iRefIdxBi[1] = bestBiPRefIdxL1;
							pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv(cMvBi[1], ePartSize, uiPartAddr, 0, iPartIdx);
							pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx);
							TComYuv* pcYuvPred = &m_acYuvPred[REF_PIC_LIST_1];
#if JVET_E0052_DMVR
							motionCompensation(pcCU, pcYuvPred, true, REF_PIC_LIST_1, iPartIdx);
#else
							motionCompensation(pcCU, pcYuvPred, REF_PIC_LIST_1, iPartIdx);
#endif

							uiMotBits[0] = uiBits[0] - uiMbBits[0];
							uiMotBits[1] = uiMbBits[1];

							if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 1)
							{
#if JVET_D0123_ME_CTX_LUT_BITS
								uiMotBits[1] += xRefFrameIdxBits(pcCU, bestBiPRefIdxL1, (UInt)pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1));
#else
								uiMotBits[1] += bestBiPRefIdxL1 + 1;
								if (bestBiPRefIdxL1 == pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) - 1)
								{
									uiMotBits[1]--;
								}
#endif
							}

							uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][AMVP_MAX_NUM_CANDS];

							uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];

							cMvTemp[1][bestBiPRefIdxL1] = cMvBi[1];
						}
						else
						{
							uiMotBits[0] = uiBits[0] - uiMbBits[0];
							uiMotBits[1] = uiBits[1] - uiMbBits[1];
							uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
						}

						// 4-times iteration (default)
						Int iNumIter = 4;

						// fast encoder setting: only one iteration
						if (m_pcEncCfg->getUseFastEnc() || pcCU->getSlice()->getMvdL1ZeroFlag())
						{
							iNumIter = 1;
						}

						for (Int iIter = 0; iIter < iNumIter; iIter++)
						{
							Int         iRefList = iIter % 2;

							if (m_pcEncCfg->getUseFastEnc())
							{
								if (uiCost[0] <= uiCost[1])
								{
									iRefList = 1;
								}
								else
								{
									iRefList = 0;
								}
							}
							else if (iIter == 0)
							{
								iRefList = 0;
							}
							if (iIter == 0 && !pcCU->getSlice()->getMvdL1ZeroFlag())
							{
								pcCU->getCUMvField(RefPicList(1 - iRefList))->setAllMv(cMv[1 - iRefList], ePartSize, uiPartAddr, 0, iPartIdx);
								pcCU->getCUMvField(RefPicList(1 - iRefList))->setAllRefIdx(iRefIdx[1 - iRefList], ePartSize, uiPartAddr, 0, iPartIdx);
								TComYuv*  pcYuvPred = &m_acYuvPred[1 - iRefList];
#if JVET_E0052_DMVR
								motionCompensation(pcCU, pcYuvPred, true, RefPicList(1 - iRefList), iPartIdx);
#else
								motionCompensation(pcCU, pcYuvPred, RefPicList(1 - iRefList), iPartIdx);
#endif
							}

							RefPicList  eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

							if (pcCU->getSlice()->getMvdL1ZeroFlag())
							{
								iRefList = 0;
								eRefPicList = REF_PIC_LIST_0;
							}

							Bool bChanged = false;

							iRefStart = 0;
							iRefEnd = pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1;

							for (Int iRefIdxTemp = iRefStart; iRefIdxTemp <= iRefEnd; iRefIdxTemp++)
							{
								uiBitsTemp = uiMbBits[2] + uiMotBits[1 - iRefList];
								if (pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1)
								{
#if JVET_D0123_ME_CTX_LUT_BITS
									uiBitsTemp += xRefFrameIdxBits(pcCU, iRefIdxTemp, (UInt)pcCU->getSlice()->getNumRefIdx(eRefPicList));
#else
									uiBitsTemp += iRefIdxTemp + 1;
									if (iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1)
									{
										uiBitsTemp--;
									}
#endif
								}
								uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];
								// call ME
#if JVET_E0076_MULTI_PEL_MVD
								xCopyAMVPInfo(&aacAMVPInfo[iRefList][iRefIdxTemp], pcCU->getCUMvField(eRefPicList)->getAMVPInfo());
								xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, true);
#else
								xMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, true);
#endif

								xCopyAMVPInfo(&aacAMVPInfo[iRefList][iRefIdxTemp], pcCU->getCUMvField(eRefPicList)->getAMVPInfo());
								xCheckBestMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
									, uiPartAddr
#endif
								);

#if JVET_E0076_MULTI_PEL_MVD
								assert(cMvPredBi[iRefList][iRefIdxTemp] == aacAMVPInfo[iRefList][iRefIdxTemp].m_acMvCand[aaiMvpIdxBi[iRefList][iRefIdxTemp]]);
#endif
								if (uiCostTemp < uiCostBi)
								{
									bChanged = true;

									cMvBi[iRefList] = cMvTemp[iRefList][iRefIdxTemp];
									iRefIdxBi[iRefList] = iRefIdxTemp;

									uiCostBi = uiCostTemp;
									uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2] - uiMotBits[1 - iRefList];
									uiBits[2] = uiBitsTemp;

									if (iNumIter != 1)
									{
										//  Set motion
										pcCU->getCUMvField(eRefPicList)->setAllMv(cMvBi[iRefList], ePartSize, uiPartAddr, 0, iPartIdx);
										pcCU->getCUMvField(eRefPicList)->setAllRefIdx(iRefIdxBi[iRefList], ePartSize, uiPartAddr, 0, iPartIdx);

										TComYuv* pcYuvPred = &m_acYuvPred[iRefList];
#if JVET_E0052_DMVR
										motionCompensation(pcCU, pcYuvPred, true, eRefPicList, iPartIdx);
#else
										motionCompensation(pcCU, pcYuvPred, eRefPicList, iPartIdx);
#endif
									}
								}
							} // for loop-iRefIdxTemp

							if (!bChanged)
							{
								if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])
								{
									xCopyAMVPInfo(&aacAMVPInfo[0][iRefIdxBi[0]], pcCU->getCUMvField(REF_PIC_LIST_0)->getAMVPInfo());
									xCheckBestMVP(pcCU, REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], uiBits[2], uiCostBi
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
										, uiPartAddr
#endif
									);
									if (!pcCU->getSlice()->getMvdL1ZeroFlag())
									{
										xCopyAMVPInfo(&aacAMVPInfo[1][iRefIdxBi[1]], pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo());
										xCheckBestMVP(pcCU, REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], uiBits[2], uiCostBi
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
											, uiPartAddr
#endif
										);
									}
								}
								break;
							}
						} // for loop-iter
					} // if (B_SLICE)
#if COM16_C806_LARGE_CTU && !JVET_C0024_QTBT
				}
#endif

#if AMP_MRG
			} //end if bTestNormalMC
#endif
			  //  Clear Motion Field
			pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);

			pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if !JVET_C0024_QTBT
			UInt uiMEBits = 0;
#endif
			// Set Motion Field_
			cMv[1] = mvValidList1;
			iRefIdx[1] = refIdxValidList1;
			uiBits[1] = bitsValidList1;
			uiCost[1] = costValidList1;

#if AMP_MRG
			if (bTestNormalMC)
			{
#endif
#if COM16_C1016_AFFINE
				uiLastModeTemp = uiLastMode;
#endif

				if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])
				{
					uiLastMode = 2;
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMv(cMvBi[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iRefIdxBi[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv(cMvBi[1], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx);

					TempMv = cMvBi[0] - cMvPredBi[0][iRefIdxBi[0]];
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(TempMv, ePartSize, uiPartAddr, 0, iPartIdx);

					TempMv = cMvBi[1] - cMvPredBi[1][iRefIdxBi[1]];
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(TempMv, ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setInterDirSubParts(3, uiPartAddr, iPartIdx, pcCU->getDepth(0));

					pcCU->setMVPIdxSubParts(aaiMvpIdxBi[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(aaiMvpNum[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPIdxSubParts(aaiMvpIdxBi[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(aaiMvpNum[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if !JVET_C0024_QTBT
					uiMEBits = uiBits[2];
#endif
				}
				else if (uiCost[0] <= uiCost[1])
				{
					uiLastMode = 0;
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMv(cMv[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iRefIdx[0], ePartSize, uiPartAddr, 0, iPartIdx);

					TempMv = cMv[0] - cMvPred[0][iRefIdx[0]];
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(TempMv, ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setInterDirSubParts(1, uiPartAddr, iPartIdx, pcCU->getDepth(0));

					pcCU->setMVPIdxSubParts(aaiMvpIdx[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(aaiMvpNum[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if !JVET_C0024_QTBT
					uiMEBits = uiBits[0];
#endif
				}
				else
				{
					uiLastMode = 1;
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv(cMv[1], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdx[1], ePartSize, uiPartAddr, 0, iPartIdx);

					TempMv = cMv[1] - cMvPred[1][iRefIdx[1]];
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(TempMv, ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setInterDirSubParts(2, uiPartAddr, iPartIdx, pcCU->getDepth(0));

					pcCU->setMVPIdxSubParts(aaiMvpIdx[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(aaiMvpNum[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if !JVET_C0024_QTBT
					uiMEBits = uiBits[1];
#endif
				}
#if AMP_MRG
			} // end if bTestNormalMC
#endif

#if !JVET_C0024_QTBT
			if (pcCU->getPartitionSize(uiPartAddr) != SIZE_2Nx2N)
			{
#if JVET_D0123_ME_CTX_LUT_BITS
				UInt uiMrgFlagBits[2] = { 0, 0 };
				xGetMrgFlagBits(pcCU, uiPartAddr, uiMrgFlagBits);
#endif
				UInt uiMRGInterDir = 0;
				TComMvField cMRGMvField[2];
				UInt uiMRGIndex = 0;

				UInt uiMEInterDir = 0;
				TComMvField cMEMvField[2];

				m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(uiPartAddr));

#if AMP_MRG
				// calculate ME cost
				Distortion uiMEError = std::numeric_limits<Distortion>::max();
				Distortion uiMECost = std::numeric_limits<Distortion>::max();

				if (bTestNormalMC)
				{
					xGetInterPredictionError(pcCU, pcOrgYuv, iPartIdx, uiMEError, m_pcEncCfg->getUseHADME());
#if JVET_D0123_ME_CTX_LUT_BITS
					uiMEBits += uiMrgFlagBits[0];
					uiMECost = uiMEError + m_pcRdCost->getCostSearch(uiMEBits);
#else
					uiMECost = uiMEError + m_pcRdCost->getCost(uiMEBits);
#endif
				}
#else
				// calculate ME cost
				Distortion uiMEError = std::numeric_limits<Distortion>::max();
				xGetInterPredictionError(pcCU, pcOrgYuv, iPartIdx, uiMEError, m_pcEncCfg->getUseHADME());
#if JVET_D0123_ME_CTX_LUT_BITS
				uiMEBits += uiMrgFlagBits[0];
				Distortion uiMECost = uiMEError + m_pcRdCost->getCostSearch(uiMEBits);
#else
				Distortion uiMECost = uiMEError + m_pcRdCost->getCost(uiMEBits);
#endif
#endif
				// save ME result.
				uiMEInterDir = pcCU->getInterDir(uiPartAddr);
				pcCU->getMvField(pcCU, uiPartAddr, REF_PIC_LIST_0, cMEMvField[0]);
				pcCU->getMvField(pcCU, uiPartAddr, REF_PIC_LIST_1, cMEMvField[1]);

#if JVET_D0123_ME_CTX_LUT_BITS && VCEG_AZ07_IMV
				if (uiMEError != std::numeric_limits<Distortion>::max() && iPartIdx == iNumPart - 1)
				{
					Bool bNonZeroMvd = false;
					for (Int iPartIdxTmp = 0; iPartIdxTmp < iNumPart; iPartIdxTmp++)
					{
						UInt uiPartAddrTmp = 0;
						pcCU->getPartIndexAndSize(iPartIdxTmp, uiPartAddrTmp, iRoiWidth, iRoiHeight);
						for (UInt uiRefListIdx = 0; uiRefListIdx < 2 && (!bNonZeroMvd); uiRefListIdx++)
						{
							if (pcCU->getSlice()->getNumRefIdx(RefPicList(uiRefListIdx)) > 0)
							{
								bNonZeroMvd |= (pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd(uiPartAddrTmp).getHor() != 0);
								bNonZeroMvd |= (pcCU->getCUMvField(RefPicList(uiRefListIdx))->getMvd(uiPartAddrTmp).getVer() != 0);
							}
						}
					}
					if (bNonZeroMvd && pcCU->getSlice()->getSPS()->getIMV())
					{
						uiMEError += m_pcRdCost->getCostSearch(m_uiBitsIMVFlag[pcCU->getiMVFlag(uiPartAddr) ? 1 : 0]);
					}
				}
#endif

				// find Merge result
				Distortion uiMRGCost = std::numeric_limits<Distortion>::max();
#if VCEG_AZ07_IMV
				if (pcCU->getSlice()->getSPS()->getIMV())
				{
					pcCU->setMergeFlagSubParts(true, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
#endif
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
				memset(eMergeCandTypeNieghors, MGR_TYPE_DEFAULT_N, MRG_MAX_NUM_CANDS * sizeof(UChar));
#endif
				xMergeEstimation(pcCU, pcOrgYuv, iPartIdx, uiMRGInterDir, cMRGMvField, uiMRGIndex, uiMRGCost, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
					, eMergeCandTypeNieghors
					, m_pMvFieldSP
					, m_phInterDirSP
#endif 
				);
#if JVET_D0123_ME_CTX_LUT_BITS
				assert(uiMRGCost != std::numeric_limits<Distortion>::max() && uiMRGCost + m_pcRdCost->getCostSearch(uiMrgFlagBits[1]) > uiMRGCost);
				uiMRGCost += m_pcRdCost->getCostSearch(uiMrgFlagBits[1]);
#endif
#if VCEG_AZ07_FRUC_MERGE
				UInt uiFRUCMgrCost = MAX_UINT;
				UChar uhFRUCMode = 0;
				if (pcCU->getSlice()->getSPS()->getUseFRUCMgrMode())
				{
					xFRUCMgrEstimation(pcCU, pcOrgYuv, iPartIdx,
						m_pMvFieldFRUC, m_phInterDirFRUC, m_phFRUCRefineDist, m_phFRUCSBlkRefineDist,
						uiFRUCMgrCost, uhFRUCMode);
#if JVET_D0123_ME_CTX_LUT_BITS
					if (uiFRUCMgrCost < MAX_UINT)
					{
						uiFRUCMgrCost += m_pcRdCost->getCostSearch(uiMrgFlagBits[1]);
					}
					uiMRGCost += m_pcRdCost->getCostSearch(xGetFRUCMgrModeBits(pcCU, uiPartAddr, iPartIdx, FRUC_MERGE_OFF));
#endif
				}
#endif

				if (uiMRGCost < uiMECost
#if VCEG_AZ07_FRUC_MERGE
					|| uiFRUCMgrCost < uiMECost
#endif
					)
				{
					// set Merge result

#if VCEG_AZ07_FRUC_MERGE
					if (uiFRUCMgrCost < uiMRGCost)
					{
						pcCU->setFRUCMgrModeSubParts(uhFRUCMode, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
						pcCU->setMergeFlagSubParts(true, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
						for (Int y = 0, yRasterOffset = 0; y < iRoiHeight; y += 4, yRasterOffset += pcCU->getPic()->getNumPartInCtuWidth())
						{
							for (Int x = 0, xRasterOffset = 0; x < iRoiWidth; x += 4, xRasterOffset++)
							{
								UInt idx = g_auiRasterToZscan[g_auiZscanToRaster[pcCU->getZorderIdxInCtu() + uiPartAddr] + yRasterOffset + xRasterOffset] - pcCU->getZorderIdxInCtu();
								pcCU->getCUMvField(REF_PIC_LIST_0)->setMvFieldSP(pcCU, idx, m_pMvFieldFRUC[idx << 1], 4, 4);
								pcCU->getCUMvField(REF_PIC_LIST_1)->setMvFieldSP(pcCU, idx, m_pMvFieldFRUC[(idx << 1) + 1], 4, 4);
								pcCU->setInterDir(idx, m_phInterDirFRUC[idx]);
							}
						}
					}
					else
					{
						pcCU->setFRUCMgrModeSubParts(FRUC_MERGE_OFF, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#endif
						pcCU->setMergeFlagSubParts(true, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
						pcCU->setMergeIndexSubParts(uiMRGIndex, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
						pcCU->setMergeTypeSubParts(eMergeCandTypeNieghors[uiMRGIndex], uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
						if (eMergeCandTypeNieghors[uiMRGIndex])
						{
							UInt uiSPAddr;
							Int iNumSPInOneLine, iNumSP, iSPWidth, iSPHeight;
							pcCU->getSPPara(iRoiWidth, iRoiHeight, iNumSP, iNumSPInOneLine, iSPWidth, iSPHeight);
#if JVET_C0035_ATMVP_SIMPLIFICATION
							UInt uiSPListIndex = eMergeCandTypeNieghors[uiMRGIndex];
#else
							UInt uiSPListIndex = eMergeCandTypeNieghors[uiMRGIndex] == MGR_TYPE_SUBPU_TMVP ? 0 : 1;
#endif
							for (Int iPartitionIdx = 0; iPartitionIdx < iNumSP; iPartitionIdx++)
							{
								pcCU->getSPAbsPartIdx(uiPartAddr, iSPWidth, iSPHeight, iPartitionIdx, iNumSPInOneLine, uiSPAddr);
								pcCU->setInterDirSP(m_phInterDirSP[uiSPListIndex][iPartitionIdx], uiSPAddr, iSPWidth, iSPHeight);
								pcCU->getCUMvField(REF_PIC_LIST_0)->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx], iSPWidth, iSPHeight);
								pcCU->getCUMvField(REF_PIC_LIST_1)->setMvFieldSP(pcCU, uiSPAddr, m_pMvFieldSP[uiSPListIndex][2 * iPartitionIdx + 1], iSPWidth, iSPHeight);
							}
							if (pcCU->getInterDir(uiPartAddr) == 3 && pcCU->isBipredRestriction(iPartIdx))
							{
								pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv(TComMv(0, 0), ePartSize, uiPartAddr, 0, iPartIdx);
								pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(-1, ePartSize, uiPartAddr, 0, iPartIdx);
								pcCU->setInterDirSubParts(1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
							}
						}
						else
						{
#endif
							pcCU->setInterDirSubParts(uiMRGInterDir, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
							pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(cMRGMvField[0], ePartSize, uiPartAddr, 0, iPartIdx);
							pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(cMRGMvField[1], ePartSize, uiPartAddr, 0, iPartIdx);
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
						}
#endif
#if VCEG_AZ07_FRUC_MERGE
					}
#endif
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
				else
				{
					// set ME result
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
					pcCU->setMergeTypeSubParts(MGR_TYPE_DEFAULT_N, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#endif
#if VCEG_AZ07_FRUC_MERGE
					pcCU->setFRUCMgrModeSubParts(FRUC_MERGE_OFF, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#endif
					pcCU->setMergeFlagSubParts(false, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setInterDirSubParts(uiMEInterDir, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(cMEMvField[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(cMEMvField[1], ePartSize, uiPartAddr, 0, iPartIdx);
				}

#if COM16_C1016_AFFINE
				uiHevcCost = (uiMRGCost < uiMECost) ? uiMRGCost : uiMECost;
#endif
			}
#endif

#if COM16_C1016_AFFINE
#if !JVET_C0024_QTBT
			else    // get best regular Hevc Cost
#endif
			{
				uiHevcCost = (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1]) ? uiCostBi :
					((uiCost[0] <= uiCost[1]) ? uiCost[0] : uiCost[1]);
			}

#if JVECT_C0062_AFFINE_SIX_PARAM
			UInt uiBestAffineParam = 0;
#endif

#if JVET_C0024_QTBT
			if (pcCU->getWidth(uiPartAddr) > 8 && pcCU->getHeight(uiPartAddr) > 8 && pcCU->getSlice()->getSPS()->getUseAffine()
#else
			if (pcCU->getWidth(uiPartAddr) > 8 && ePartSize == SIZE_2Nx2N && pcCU->getSlice()->getSPS()->getUseAffine()
#endif
#if VCEG_AZ07_IMV
				&& !pcCU->getiMVFlag(uiPartAddr)
#endif
#if VCEG_AZ06_IC
				&& !pcCU->getICFlag(uiPartAddr)
#endif
#if JVET_D0077_SAVE_LOAD_ENC_INFO
				&& !bFastSkipAffine
#endif
#if SIF_AFFINE_HARMONIZATION
				&& !pcCU->getSIFFlag(uiPartAddr)
#endif
				)
			{
				// save normal hevc result
				UInt uiMRGIndex = pcCU->getMergeIndex(uiPartAddr);
				Bool bMergeFlag = pcCU->getMergeFlag(uiPartAddr);
				UInt uiInterDir = pcCU->getInterDir(uiPartAddr);

				TComMv cMvd[2];
				UInt uiMvpIdx[2], uiMvpNum[2];
				uiMvpIdx[0] = pcCU->getMVPIdx(REF_PIC_LIST_0, uiPartAddr);
				uiMvpIdx[1] = pcCU->getMVPIdx(REF_PIC_LIST_1, uiPartAddr);
				uiMvpNum[0] = pcCU->getMVPNum(REF_PIC_LIST_0, uiPartAddr);
				uiMvpNum[1] = pcCU->getMVPNum(REF_PIC_LIST_1, uiPartAddr);
				cMvd[0] = pcCU->getCUMvField(REF_PIC_LIST_0)->getMvd(uiPartAddr);
				cMvd[1] = pcCU->getCUMvField(REF_PIC_LIST_1)->getMvd(uiPartAddr);

				TComMvField cHevcMvField[2];
				pcCU->getMvField(pcCU, uiPartAddr, REF_PIC_LIST_0, cHevcMvField[0]);
				pcCU->getMvField(pcCU, uiPartAddr, REF_PIC_LIST_1, cHevcMvField[1]);

				// do affine ME & Merge
#if JVECT_C0062_AFFINE_SIX_PARAM
				TComMv cMVBest[2][3], bestMVD[2][3];
				Int bestDir = 0;
				Int iBestRefIdx[2] = { 0, 0 };
				Int bestMVPIdx[2] = { 0, 0 };
				Int bestMVPNum[2] = { 0, 0 };

#if AFFINE_INHERIT_FOUR_PARAM
			TComMv cMvAffineParam[2][33][3];
#endif
      for (UInt uiAffineParam = 0; uiAffineParam < 2; uiAffineParam++)
      {
        predAffineInterSearch(pcCU, pcOrgYuv, iPartIdx, uiLastModeTemp, uiAffineCost, cMvHevcTemp, uiAffineParam, uiBestAffineParam, cMVBest, iBestRefIdx, bestMVD, bestDir, bestMVPIdx, bestMVPNum
#if AFFINE_INHERIT_FOUR_PARAM
					, cMvAffineParam
#endif
					);
#if AFF_EARLY_TERMINATION
				if (!(pcBestCU->getMergeComplexFlag(uiPartAddr) || pcBestCU->getAffineFlag(uiPartAddr)) && (uiAffineCost > uiHevcCost))
				{
					break;
				}
#else
        if (uiAffineCost > uiHevcCost * 1.05)
        {
          break;
        }
#endif
      }
      //4,6 parameter switch flag
      uiAffineCost += m_pcRdCost->getCost(1);

				//reset to best for later MC
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
				pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);
				pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);

				pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(-1, ePartSize, uiPartAddr, 0, iPartIdx);
				pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(-1, ePartSize, uiPartAddr, 0, iPartIdx);

				//use best to set 
				if (bestDir == 3) //Bi
				{
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMVBest[0], REF_PIC_LIST_0, 0, uiBestAffineParam);
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMVBest[1], REF_PIC_LIST_1, 0, uiBestAffineParam);

					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iBestRefIdx[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iBestRefIdx[1], ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, bestMVD[0], REF_PIC_LIST_0, 0);
					pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, bestMVD[1], REF_PIC_LIST_1, 0);

					pcCU->setInterDirSubParts(3, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

					pcCU->setMVPIdxSubParts(bestMVPIdx[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(bestMVPNum[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

					pcCU->setMVPIdxSubParts(bestMVPIdx[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(bestMVPNum[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
				else if (bestDir == 1) //L0
				{
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMVBest[0], REF_PIC_LIST_0, 0, uiBestAffineParam);
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iBestRefIdx[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, bestMVD[0], REF_PIC_LIST_0, 0);
					pcCU->setInterDirSubParts(1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

					pcCU->setMVPIdxSubParts(bestMVPIdx[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(bestMVPNum[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
				else //2 L1
				{
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMVBest[1], REF_PIC_LIST_1, 0, uiBestAffineParam);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iBestRefIdx[1], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, bestMVD[1], REF_PIC_LIST_1, 0);
					pcCU->setInterDirSubParts(2, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

					pcCU->setMVPIdxSubParts(bestMVPIdx[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(bestMVPNum[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
#else
				predAffineInterSearch(pcCU, pcOrgYuv, iPartIdx, uiLastModeTemp, uiAffineCost, cMvHevcTemp);
#endif

#if JVET_D0123_ME_CTX_LUT_BITS
				UInt uiAffineFlagBits[2] = { 0, 0 };
				xGetAffineFlagBits(pcCU, uiPartAddr, iPartIdx, uiAffineFlagBits);
				assert(uiHevcCost != MAX_UINT);
				uiHevcCost += m_pcRdCost->getCostSearch(uiAffineFlagBits[0]);
				assert(uiAffineCost != MAX_UINT);
				uiAffineCost += m_pcRdCost->getCostSearch(uiAffineFlagBits[1]);
#endif
				// check best mode, normal or affine
				if (uiHevcCost <= uiAffineCost)
				{
					// set hevc me result
					pcCU->setAffineFlagSubParts(false, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMergeFlagSubParts(bMergeFlag, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMergeIndexSubParts(uiMRGIndex, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setInterDirSubParts(uiInterDir, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(cHevcMvField[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(cHevcMvField[1], ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cMvd[0], ePartSize, uiPartAddr, 0, iPartIdx);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cMvd[1], ePartSize, uiPartAddr, 0, iPartIdx);

					pcCU->setMVPIdxSubParts(uiMvpIdx[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(uiMvpNum[0], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPIdxSubParts(uiMvpIdx[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(uiMvpNum[1], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}
				else
				{
					assert(pcCU->isAffine(uiPartAddr));
					uiLastMode = uiLastModeTemp;
				}
			}
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM 
			Bool b6Param = uiBestAffineParam;
			pcCU->setAffineParamFlagSubParts(b6Param, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#endif

			//  MC
#if JVET_E0052_DMVR
			motionCompensation(pcCU, pcPredYuv, true, REF_PIC_LIST_X, iPartIdx);
#else
			motionCompensation(pcCU, pcPredYuv, REF_PIC_LIST_X, iPartIdx);
#endif

#if JVET_D0123_ME_CTX_LUT_BITS && (!JVET_C0024_QTBT) //encode motion info && update ctx model
			if (iPartIdx < iNumPart - 1)
			{
				//m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_PU_NEXT_BEST]);
				loadCtxMeSbacCoder(pcCU);
				m_pcEntropyCoder->encodePuMotionInfo(pcCU, uiPartAddr, iPartIdx);
				//m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_PU_NEXT_BEST]);
				storeCtxMeSbacCoder(pcCU);
				//m_pcRDGoOnSbacCoder->resetBits();
			}
#endif

		} //  end of for ( Int iPartIdx = 0; iPartIdx < iNumPart; iPartIdx++ )

#if COM16_C806_OBMC
		pcPredYuv->copyToPartYuv(pcPredYuvWoOBMC, 0);
		subBlockOBMC(pcCU, 0, pcPredYuv, pcTmpYuv1, pcTmpYuv2);
		//re-fetch original YUV
		pcOrgYuv->copyFromPicYuv(pcCU->getPic()->getPicYuvOrg(), pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
#endif
		setWpScalingDistParam(pcCU, -1, REF_PIC_LIST_X);

		return;
	}


	// AMVP
	Void TEncSearch::xEstimateMvPredAMVP(TComDataCU* pcCU, TComYuv* pcOrgYuv, UInt uiPartIdx, RefPicList eRefPicList, Int iRefIdx, TComMv& rcMvPred, Bool bFilled, Distortion* puiDistBiP)
	{
		AMVPInfo*  pcAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAMVPInfo();

		TComMv     cBestMv;
		Int        iBestIdx = 0;
		TComMv     cZeroMv;
		TComMv     cMvPred;
		Distortion uiBestCost = std::numeric_limits<Distortion>::max();
		UInt       uiPartAddr = 0;
		Int        iRoiWidth, iRoiHeight;
		Int        i;

		pcCU->getPartIndexAndSize(uiPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
		// Fill the MV Candidates
		if (!bFilled)
		{
			pcCU->fillMvpCand(uiPartIdx, uiPartAddr, eRefPicList, iRefIdx, pcAMVPInfo
#if VCEG_AZ07_FRUC_MERGE
				, this
#endif
			);
		}

		// initialize Mvp index & Mvp
		iBestIdx = 0;
		cBestMv = pcAMVPInfo->m_acMvCand[0];
		if (pcAMVPInfo->iN <= 1)
		{
			rcMvPred = cBestMv;

			pcCU->setMVPIdxSubParts(iBestIdx, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(pcAMVPInfo->iN, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));

			if (pcCU->getSlice()->getMvdL1ZeroFlag() && eRefPicList == REF_PIC_LIST_1)
			{
				(*puiDistBiP) = xGetTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, rcMvPred, 0, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
			}
			return;
		}

		if (bFilled)
		{
			assert(pcCU->getMVPIdx(eRefPicList, uiPartAddr) >= 0);
			rcMvPred = pcAMVPInfo->m_acMvCand[pcCU->getMVPIdx(eRefPicList, uiPartAddr)];
			return;
		}

#if !COM16_C806_LARGE_CTU
		m_cYuvPredTemp.clear();
#endif

		//-- Check Minimum Cost.
		for (i = 0; i < pcAMVPInfo->iN; i++)
		{
			Distortion uiTmpCost;
			uiTmpCost = xGetTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, pcAMVPInfo->m_acMvCand[i], i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
			if (uiBestCost > uiTmpCost)
			{
				uiBestCost = uiTmpCost;
				cBestMv = pcAMVPInfo->m_acMvCand[i];
				iBestIdx = i;
				(*puiDistBiP) = uiTmpCost;
			}
		}

#if !COM16_C806_LARGE_CTU
		m_cYuvPredTemp.clear();
#endif

		// Setting Best MVP
		rcMvPred = cBestMv;
		pcCU->setMVPIdxSubParts(iBestIdx, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
		pcCU->setMVPNumSubParts(pcAMVPInfo->iN, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
		return;
	}

#if JVET_D0123_ME_CTX_LUT_BITS
	UInt TEncSearch::xRefFrameIdxBits(TComDataCU* pcCU, Int iRefFrame, UInt uiRefNum)
	{
		UInt uiCtx = 0;
		UInt uiBits = m_pcPuMeEstBitsSbac->refIdxBits[uiCtx][(iRefFrame == 0 ? 0 : 1)];
		if (iRefFrame > 0)
		{
			//UInt uiRefNum = pcCU->getSlice()->getNumRefIdx( eRefList ) - 2;
			uiRefNum = uiRefNum - 2;
			//pCtx++;
			//iRefFrame--;
			uiCtx++;
			iRefFrame--;
			for (UInt ui = 0; ui < uiRefNum; ++ui)
			{
				const UInt uiSymbol = ui == iRefFrame ? 0 : 1;
				if (ui == 0)
				{
					//m_pcBinIf->encodeBin( uiSymbol, *pCtx );
					uiBits += m_pcPuMeEstBitsSbac->refIdxBits[uiCtx][uiSymbol];
				}
				else
				{
					//m_pcBinIf->encodeBinEP( uiSymbol );
					uiBits += (UInt)EPBIT;
				}
				if (uiSymbol == 0)
				{
					break;
				}
			}
		}
		return uiBits;
	}
#if JVET_C0024_FAST_MRG
	Void TEncSearch::getMrgCandBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiMrgIdxBits[MRG_MAX_NUM_CANDS])
	{
		loadCtxMe(pcCU);
		loadCtxMeSbacCoder(pcCU);
		m_pcEntropyCoder->estimatePuMeBit(m_pcPuMeEstBitsSbac);
		for (UInt i = 0; i < MRG_MAX_NUM_CANDS; i++)
		{
			uiMrgIdxBits[i] = xGetMrgIdxBits(pcCU, uiAbsPartIdx, i);
		}
	}
#endif
	Void TEncSearch::xUpdateMvpIdxBits(TComDataCU* pcCU)
	{
		for (Int iNum = 0; iNum < AMVP_MAX_NUM_CANDS + 1; iNum++)
		{
			for (Int iIdx = 0; iIdx < AMVP_MAX_NUM_CANDS; iIdx++)
			{
				if (iIdx < iNum)
				{
					const Int iOffset = 1;
					//so far: only works when AMVP_MAX_NUM_CANDS is equal to 2
					UInt uiSymbol = (UInt)iIdx;
					Int iMaxSymbol = iNum - 1;
					if (iMaxSymbol == 0)
					{
						m_auiMVPIdxCost[iIdx][iNum] = 0;
					}
					else
					{
						UInt uiBits = m_pcPuMeEstBitsSbac->mvpIdxBits[0][uiSymbol ? 1 : 0];
						if (uiSymbol != 0)
						{
							Bool bCodeLast = (iMaxSymbol > uiSymbol);
							while (--uiSymbol)
							{
								//m_pcBinIf->encodeBin( 1, pcSCModel[ iOffset ] );
								uiBits += m_pcPuMeEstBitsSbac->mvpIdxBits[iOffset][1];
							}
							if (bCodeLast)
							{
								//m_pcBinIf->encodeBin( 0, pcSCModel[ iOffset ] );
								uiBits += m_pcPuMeEstBitsSbac->mvpIdxBits[iOffset][0];
							}
						}
						m_auiMVPIdxCost[iIdx][iNum] = uiBits;
					}
				}
				else
				{
					m_auiMVPIdxCost[iIdx][iNum] = MAX_INT;
				}
			}
		}
	}

	Void TEncSearch::loadCtxMe(TComDataCU* pcCU)
	{
#if JVET_C0024_QTBT
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
		m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_PU_NEXT_BEST]->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
		UInt uiDepth = pcCU->getDepth(0);
		m_pppcRDSbacCoder[uiDepth][CI_PU_NEXT_BEST]->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
#endif
	}

	Void TEncSearch::loadCtxMeSbacCoder(TComDataCU* pcCU)
	{
#if JVET_C0024_QTBT
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
		m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_PU_NEXT_BEST]);
#else
		UInt uiDepth = pcCU->getDepth(0);
		m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_PU_NEXT_BEST]);
#endif
	}

	Void TEncSearch::storeCtxMeSbacCoder(TComDataCU* pcCU)
	{
#if JVET_C0024_QTBT
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
		m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_PU_NEXT_BEST]);
#else
		UInt uiDepth = pcCU->getDepth(0);
		m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_PU_NEXT_BEST]);
#endif
	}

#if COM16_C1016_AFFINE
	Void TEncSearch::xGetAffineFlagBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt iPUIdx, UInt uiAffineFlagBits[2]) {
		UInt uiCtxAffine = pcCU->getCtxAffineFlag(uiAbsPartIdx);
		for (UInt i = 0; i < 2; i++)
		{
			uiAffineFlagBits[i] = m_pcPuMeEstBitsSbac->affineFlagBits[uiCtxAffine][i];
		}
	}
#endif

	Void TEncSearch::xGetMrgFlagBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiMrgFlagBits[2])
	{
		for (UInt i = 0; i< 2; i++)
		{
			uiMrgFlagBits[i] = m_pcPuMeEstBitsSbac->mrgFlagBits[0][i];
		}
	}

	UInt TEncSearch::xGetMrgIdxBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiMergeCand)
	{
		UInt uiNumCand = pcCU->getSlice()->getMaxNumMergeCand();
		UInt uiBits = 0;
		if (uiNumCand > 1)
		{
			for (UInt ui = 0; ui < uiNumCand - 1; ++ui)
			{
				const UInt uiSymbol = ui == uiMergeCand ? 0 : 1;
#if COM16_C806_GEN_MRG_IMPROVEMENT
				//m_pcBinIf->encodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, (ui > NUM_MERGE_IDX_EXT_CTX-1? NUM_MERGE_IDX_EXT_CTX-1:ui) ) );
				uiBits += m_pcPuMeEstBitsSbac->mrgIdxBits[(ui > NUM_MERGE_IDX_EXT_CTX - 1 ? NUM_MERGE_IDX_EXT_CTX - 1 : ui)][uiSymbol];
#else
				if (ui == 0)
				{
					//m_pcBinIf->encodeBin( uiSymbol, m_cCUMergeIdxExtSCModel.get( 0, 0, 0 ) );
					uiBits += m_pcPuMeEstBitsSbac->mrgIdxBits[0][uiSymbol];
				}
				else
				{
					uiBits += (UInt)EPBIT;
				}
#endif
				if (uiSymbol == 0)
				{
					break;
				}
			}
		}
		return uiBits;
	}

#if VCEG_AZ07_FRUC_MERGE
	UInt TEncSearch::xGetFRUCMgrModeBits(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt iPUIdx, const UChar uhFRUCMode) {
		if (!pcCU->getSlice()->getSPS()->getUseFRUCMgrMode())
			return MAX_UINT;
		UInt uiBits = 0;
		//UInt uiFirstBin = pcCU->getFRUCMgrMode( uiAbsPartIdx ) != FRUC_MERGE_OFF;
		UInt uiFirstBin = uhFRUCMode != FRUC_MERGE_OFF;
		//m_pcBinIf->encodeBin( uiFirstBin, m_cCUFRUCMgrModeSCModel.get( 0 ,  0 , pcCU->getCtxFRUCMgrMode( uiAbsPartIdx ) ) );
		uiBits += m_pcPuMeEstBitsSbac->frucMrgBits[pcCU->getCtxFRUCMgrMode(uiAbsPartIdx)][uiFirstBin];
		if (uiFirstBin)
		{
			if (pcCU->getSlice()->isInterP())
			{
				//assert( pcCU->getFRUCMgrMode( uiAbsPartIdx ) == FRUC_MERGE_TEMPLATE );
				assert(uhFRUCMode == FRUC_MERGE_TEMPLATE);
			}
			else
			{
				UInt uiSecondBin = uhFRUCMode == FRUC_MERGE_BILATERALMV;
				//m_pcBinIf->encodeBin( uiSecondBin , m_cCUFRUCMESCModel.get( 0 , 0 , pcCU->getCtxFRUCME( uiAbsPartIdx ) ) );
				uiBits += m_pcPuMeEstBitsSbac->frucMeBits[pcCU->getCtxFRUCME(uiAbsPartIdx)][uiSecondBin];
			}
		}
		return uiBits;
	}
#endif

	UInt TEncSearch::getInterDirPUBits(UInt uiInterDir, UInt uiCtx, TComDataCU* pcCU, UInt uiAbsPartIdx) {
#if COM16_C806_LARGE_CTU
		assert(uiCtx < 4);  // uiCtx=4 is only for the last bin when uiInterDir < 2
#endif
		UInt uiBits = 0;
		UInt uiSymbol = 0;
#if JVET_C0024_QTBT
		if (1)
#else
#if COM16_C806_HEVC_MOTION_CONSTRAINT_REMOVAL
		if (pcCU->getSlice()->getSPS()->getAtmvpEnableFlag() || pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8)
#else
		if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8)
#endif
#endif
		{
			//m_pcBinIf->encodeBin( uiInterDir == 2 ? 1 : 0, *( pCtx + uiCtx ) );  //uiCtx
			uiSymbol = (uiInterDir == 2) ? 1 : 0;
			uiBits += m_pcPuMeEstBitsSbac->interDirBits[uiCtx][uiSymbol];
		}

		if (uiInterDir < 2)
		{
			//m_pcBinIf->encodeBin( uiInterDir, *( pCtx + 4 ) );  //4
			uiBits += m_pcPuMeEstBitsSbac->interDirBits[4][uiInterDir];
		}
		return uiBits;
	}
#endif

	UInt TEncSearch::xGetMvpIdxBits(Int iIdx, Int iNum)
	{
		assert(iIdx >= 0 && iNum >= 0 && iIdx < iNum);

		if (iNum == 1)
		{
			return 0;
		}

		UInt uiLength = 1;
		Int iTemp = iIdx;
		if (iTemp == 0)
		{
			return uiLength;
		}

		Bool bCodeLast = (iNum - 1 > iTemp);

		uiLength += (iTemp - 1);

		if (bCodeLast)
		{
			uiLength++;
		}

		return uiLength;
	}

#if VCEG_AZ07_IMV && JVET_D0123_ME_CTX_LUT_BITS
	Void TEncSearch::xGetIMVFlagBits(TComDataCU* pcCU, UInt uiAbsPartIdx)
	{
		UInt uiCtxiMV = pcCU->getCtxiMVFlag(uiAbsPartIdx);
		for (UInt i = 0; i < 2; i++)
		{
			m_uiBitsIMVFlag[i] = m_pcPuMeEstBitsSbac->iMVFlagBits[uiCtxiMV][i];
		}
	}
#endif
	Void TEncSearch::xGetBlkBits(PartSize eCUMode, Bool bPSlice, Int iPartIdx, UInt uiLastMode, UInt uiBlkBit[3]
#if JVET_D0123_ME_CTX_LUT_BITS
		, TComDataCU* pcCU, UInt uiAbsPartIdx
#endif
	)
	{
#if JVET_D0123_ME_CTX_LUT_BITS
		if (bPSlice)
		{
			memset(uiBlkBit, 0, 3 * sizeof(UInt));
			return;
		}
		UInt uiCtx = pcCU->getCtxInterDir(uiAbsPartIdx);
		for (UInt i = 0; i < 3; i++)
		{
			uiBlkBit[i] = getInterDirPUBits(i, uiCtx, pcCU, uiAbsPartIdx);
		}
#else
#if JVET_C0024_QTBT
		uiBlkBit[0] = (!bPSlice) ? 3 : 1;
		uiBlkBit[1] = 3;
		uiBlkBit[2] = 5;
#else
		if (eCUMode == SIZE_2Nx2N)
		{
			uiBlkBit[0] = (!bPSlice) ? 3 : 1;
			uiBlkBit[1] = 3;
			uiBlkBit[2] = 5;
		}
		else if ((eCUMode == SIZE_2NxN || eCUMode == SIZE_2NxnU) || eCUMode == SIZE_2NxnD)
		{
			UInt aauiMbBits[2][3][3] = { { { 0,0,3 },{ 0,0,0 },{ 0,0,0 } } ,{ { 5,7,7 },{ 7,5,7 },{ 9 - 3,9 - 3,9 - 3 } } };
			if (bPSlice)
			{
				uiBlkBit[0] = 3;
				uiBlkBit[1] = 0;
				uiBlkBit[2] = 0;
			}
			else
			{
				::memcpy(uiBlkBit, aauiMbBits[iPartIdx][uiLastMode], 3 * sizeof(UInt));
			}
		}
		else if ((eCUMode == SIZE_Nx2N || eCUMode == SIZE_nLx2N) || eCUMode == SIZE_nRx2N)
		{
			UInt aauiMbBits[2][3][3] = { { { 0,2,3 },{ 0,0,0 },{ 0,0,0 } } ,{ { 5,7,7 },{ 7 - 2,7 - 2,9 - 2 },{ 9 - 3,9 - 3,9 - 3 } } };
			if (bPSlice)
			{
				uiBlkBit[0] = 3;
				uiBlkBit[1] = 0;
				uiBlkBit[2] = 0;
			}
			else
			{
				::memcpy(uiBlkBit, aauiMbBits[iPartIdx][uiLastMode], 3 * sizeof(UInt));
			}
		}
		else if (eCUMode == SIZE_NxN)
		{
			uiBlkBit[0] = (!bPSlice) ? 3 : 1;
			uiBlkBit[1] = 3;
			uiBlkBit[2] = 5;
		}
		else
		{
			printf("Wrong!\n");
			assert(0);
		}
#endif
#endif
	}

	Void TEncSearch::xCopyAMVPInfo(AMVPInfo* pSrc, AMVPInfo* pDst)
	{
		pDst->iN = pSrc->iN;
		for (Int i = 0; i < pSrc->iN; i++)
		{
			pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
		}
	}

#if JVET_E0076_MULTI_PEL_MVD
	Void TEncSearch::intMvRefine(TComDataCU* pcCU, RefPicList eRefPicList, TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, TComMv& rcMv, TComMv& rcMvPred, Int& riMVPIdx, UInt& ruiBits, Distortion& ruiCost, Double fWeight)
	{
		m_pcRdCost->setDistParam(pcPatternKey, piRefY, iRefStride, 1, m_cDistParam, m_pcEncCfg->getUseHADME());

		assert(intMvRefineNeeded(pcCU, 0));

		//check 9 points and 2 avmp predictor candiate to find the best matched one

		AMVPInfo* pcAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAMVPInfo();

		assert(pcAMVPInfo->m_acMvCand[riMVPIdx] == rcMvPred);

		m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(0));
		m_pcRdCost->setCostScale(0);

#if JVET_D0123_ME_CTX_LUT_BITS
		iCostScale = 0;
#endif

		Distortion  uiDist, uiSATD = 0;
		Distortion  uiBestDist = std::numeric_limits<Distortion>::max();

		TComMv cBestMv = rcMv;
		TComMv cBaseMvd[2];
		Int iBestBits = 0;
		Int iBestMVPIdx = riMVPIdx;
		Int testPos[9][2] = { { 0, 0 },{ -1, -1 },{ -1, 0 },{ -1, 1 },{ 0, -1 },{ 0, 1 },{ 1, -1 },{ 1, 0 },{ 1, 1 } };

		//  pcCU->xRoundMV(rcMv, 0);

		cBaseMvd[0] = (rcMv - pcAMVPInfo->m_acMvCand[0]);
		cBaseMvd[1] = (rcMv - pcAMVPInfo->m_acMvCand[1]);

		assert((cBaseMvd[0].getHor() & 3) == 0 && (cBaseMvd[0].getVer() & 3) == 0);
		assert((cBaseMvd[1].getHor() & 3) == 0 && (cBaseMvd[1].getVer() & 3) == 0);

		Int mvOffset = 1 << (2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);

#if JVET_E0076_MULTI_PEL_MVD
		pcCU->xRoundMV(cBaseMvd[0], 0);
		pcCU->xRoundMV(cBaseMvd[1], 0);

		if (pcCU->getiMVFlag(0) == 2)
		{
			mvOffset <<= MULTI_PEL_MVD_BITS;
		}
#endif

		for (Int pos = 0; pos < 9; pos++)
		{
			TComMv cTestMv[2];
			for (Int iMVPIdx = 0; iMVPIdx < pcAMVPInfo->iN; iMVPIdx++)
			{
				cTestMv[iMVPIdx].set(testPos[pos][0] * mvOffset, testPos[pos][1] * mvOffset);
				cTestMv[iMVPIdx] += cBaseMvd[iMVPIdx];
				cTestMv[iMVPIdx] += pcAMVPInfo->m_acMvCand[iMVPIdx];

				if (iMVPIdx == 0 || cTestMv[0] != cTestMv[1])
				{
					TComMv cTempMV = cTestMv[iMVPIdx];
					pcCU->clipMv(cTempMV);

					Pel* piRefPos = piRefY + iRefStride * (cTempMV.getVer() >> (2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE)) + (cTempMV.getHor() >> (2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE));

					m_cDistParam.pCur = piRefPos;
					m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
					uiDist = uiSATD = (Distortion)(m_cDistParam.DistFunc(&m_cDistParam) * fWeight);
				}
				else
				{
					uiDist = uiSATD;
					assert(cTestMv[0] == cTestMv[1]);
				}

				Int iMvBits = m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];
#if JVET_D0123_ME_CTX_LUT_BITS
				pcPatternKey->m_mvPred = &pcAMVPInfo->m_acMvCand[iMVPIdx];
				iMvBits += xMvdBits(pcPatternKey, cTestMv[iMVPIdx].getHor(), cTestMv[iMVPIdx].getVer());
				uiDist += (UInt)m_pcRdCost->getCostSearch(iMvBits);
#else
				m_pcRdCost->setPredictor(pcAMVPInfo->m_acMvCand[iMVPIdx]);
				iMvBits += m_pcRdCost->getBits(
					cTestMv[iMVPIdx].getHor(), cTestMv[iMVPIdx].getVer()
#if VCEG_AZ07_IMV
					, pcCU->getiMVFlag(0)
#endif
				);
				uiDist += m_pcRdCost->getCost(cTestMv[iMVPIdx].getHor(), cTestMv[iMVPIdx].getVer()
#if VCEG_AZ07_IMV
					, pcPatternKey->getImvFlag()
#endif
				);
#endif

				if (uiDist < uiBestDist)
				{
					uiBestDist = uiDist;
					cBestMv = cTestMv[iMVPIdx];
					iBestMVPIdx = iMVPIdx;
					iBestBits = iMvBits;
				}
			}
		}

		rcMv = cBestMv;
		rcMvPred = pcAMVPInfo->m_acMvCand[iBestMVPIdx];
		riMVPIdx = iBestMVPIdx;

		ruiBits += iBestBits;

#if JVET_D0123_ME_CTX_LUT_BITS
		ruiCost = uiBestDist - m_pcRdCost->getCostSearch(iBestBits) + m_pcRdCost->getCostSearch(ruiBits);
#else
		ruiCost = uiBestDist - m_pcRdCost->getCost(iBestBits) + m_pcRdCost->getCost(ruiBits);
#endif
	}

#endif

	Void TEncSearch::xCheckBestMVP(TComDataCU* pcCU, RefPicList eRefPicList, TComMv cMv, TComMv& rcMvPred, Int& riMVPIdx, UInt& ruiBits, Distortion& ruiCost
#if VCEG_AZ07_IMV || JVET_D0123_ME_CTX_LUT_BITS
		, UInt uiPartAddr
#endif
	)
	{

#if JVET_E0076_MULTI_PEL_MVD
		if (intMvRefineNeeded(pcCU, 0))
		{
			return;
		}
#endif

		AMVPInfo* pcAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAMVPInfo();

		assert(pcAMVPInfo->m_acMvCand[riMVPIdx] == rcMvPred);

		if (pcAMVPInfo->iN < 2)
		{
			return;
		}

		m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(0));
		m_pcRdCost->setCostScale(0);

		Int iBestMVPIdx = riMVPIdx;

		m_pcRdCost->setPredictor(rcMvPred);
#if JVET_D0123_ME_CTX_LUT_BITS
		TComPattern   tmpPattern;
		TComPattern*  pcPatternKey = &tmpPattern;
		pcPatternKey->m_pcCU = pcCU;
		pcPatternKey->m_uiPartAddr = uiPartAddr;
		pcPatternKey->m_mvPred = &rcMvPred;
		iCostScale = 0;
		Int iOrgMvBits = xMvdBits(pcPatternKey, cMv.getHor(), cMv.getVer());
#else
		Int iOrgMvBits = m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			cMv.getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, cMv.getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
			cMv.getHor(), cMv.getVer()
#endif
#if VCEG_AZ07_IMV
			, pcCU->getiMVFlag(uiPartAddr)
#endif
		);
#endif
		iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];
		Int iBestMvBits = iOrgMvBits;

		for (Int iMVPIdx = 0; iMVPIdx < pcAMVPInfo->iN; iMVPIdx++)
		{
			if (iMVPIdx == riMVPIdx)
			{
				continue;
			}
#if JVET_D0123_ME_CTX_LUT_BITS
			pcPatternKey->m_mvPred = &pcAMVPInfo->m_acMvCand[iMVPIdx];
			Int iMvBits = xMvdBits(pcPatternKey, cMv.getHor(), cMv.getVer());
#else
			m_pcRdCost->setPredictor(pcAMVPInfo->m_acMvCand[iMVPIdx]);

			Int iMvBits = m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				cMv.getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, cMv.getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
				cMv.getHor(), cMv.getVer()
#endif
#if VCEG_AZ07_IMV
				, pcCU->getiMVFlag(uiPartAddr)
#endif
			);
#endif
			iMvBits += m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];

			if (iMvBits < iBestMvBits)
			{
				iBestMvBits = iMvBits;
				iBestMVPIdx = iMVPIdx;
			}
		}

		if (iBestMVPIdx != riMVPIdx)  //if changed
		{
			rcMvPred = pcAMVPInfo->m_acMvCand[iBestMVPIdx];

			riMVPIdx = iBestMVPIdx;
			UInt uiOrgBits = ruiBits;
			ruiBits = uiOrgBits - iOrgMvBits + iBestMvBits;
#if JVET_D0123_ME_CTX_LUT_BITS
			ruiCost = (ruiCost - m_pcRdCost->getCostSearch(uiOrgBits)) + m_pcRdCost->getCostSearch(ruiBits);
#else
			ruiCost = (ruiCost - m_pcRdCost->getCost(uiOrgBits)) + m_pcRdCost->getCost(ruiBits);
#endif
		}
	}


	Distortion TEncSearch::xGetTemplateCost(TComDataCU* pcCU,
		UInt        uiPartAddr,
		TComYuv*    pcOrgYuv,
		TComYuv*    pcTemplateCand,
		TComMv      cMvCand,
		Int         iMVPIdx,
		Int         iMVPNum,
		RefPicList  eRefPicList,
		Int         iRefIdx,
		Int         iSizeX,
		Int         iSizeY
	)
	{
		Distortion uiCost = std::numeric_limits<Distortion>::max();

		TComPicYuv* pcPicYuvRef = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdx)->getPicYuvRec();

		pcCU->clipMv(cMvCand);
#if VCEG_AZ06_IC
		Bool bICFlag = pcCU->getICFlag(uiPartAddr);
#endif
		// prediction pattern
		if (pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType() == P_SLICE
#if VCEG_AZ06_IC
			&& !bICFlag
#endif
			)
		{
			xPredInterBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, &cMvCand, iSizeX, iSizeY, pcTemplateCand, true, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA));
		}
		else
		{
			xPredInterBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, &cMvCand, iSizeX, iSizeY, pcTemplateCand, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if VCEG_AZ05_BIO
				, false
#endif
#if VCEG_AZ07_FRUC_MERGE
				, false
#endif
#if VCEG_AZ06_IC
				, bICFlag
#endif
			);
		}

#if RRC_SELECTABLE_IF
		TComPic *refPic = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdx);
		Int sharpMode = pcCU->getSIFMode(&cMvCand, uiPartAddr, 1);
		if (sharpMode == SIF_SHARP)
		{
			refPic->getSharpeningFilter()->sharpBlockLuma(pcCU, uiPartAddr, iSizeX, iSizeY, sharpMode, pcTemplateCand, this, &cMvCand, false, false, pcPicYuvRef);
		}
#endif
		if (pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType() == P_SLICE
#if VCEG_AZ06_IC
			&& !bICFlag
#endif
			)
		{
			xWeightedPredictionUni(pcCU, pcTemplateCand, uiPartAddr, iSizeX, iSizeY, eRefPicList, pcTemplateCand, iRefIdx);
		}

		// calc distortion

		uiCost = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pcTemplateCand->getAddr(COMPONENT_Y, uiPartAddr), pcTemplateCand->getStride(COMPONENT_Y), pcOrgYuv->getAddr(COMPONENT_Y, uiPartAddr), pcOrgYuv->getStride(COMPONENT_Y), iSizeX, iSizeY, COMPONENT_Y, DF_SAD);
#if JVET_D0123_ME_CTX_LUT_BITS
		uiCost = uiCost + m_pcRdCost->getCostSearch(m_auiMVPIdxCost[iMVPIdx][iMVPNum]);
#else
		uiCost = (UInt)m_pcRdCost->calcRdCost(m_auiMVPIdxCost[iMVPIdx][iMVPNum], uiCost, false, DF_SAD);
#endif
		return uiCost;
	}

	Void TEncSearch::xMotionEstimation(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, RefPicList eRefPicList, TComMv* pcMvPred, Int iRefIdxPred, TComMv& rcMv,
#if JVET_E0076_MULTI_PEL_MVD
		Int& riMVPIdx,
#endif
		UInt& ruiBits, Distortion& ruiCost, Bool bBi)
	{
		UInt          uiPartAddr;
		Int           iRoiWidth;
		Int           iRoiHeight;

		TComMv        cMvHalf, cMvQter;
		TComMv        cMvSrchRngLT;
		TComMv        cMvSrchRngRB;

		TComYuv*      pcYuv = pcYuvOrg;

		assert(eRefPicList < MAX_NUM_REF_LIST_ADAPT_SR && iRefIdxPred<Int(MAX_IDX_ADAPT_SR));
		m_iSearchRange = m_aaiAdaptSR[eRefPicList][iRefIdxPred];

		Int           iSrchRng = (bBi ? m_bipredSearchRange : m_iSearchRange);

#if FAST_ME_FOR_BG_REFERENCE
		if (pcCU->getSlice()->getRefPOC(eRefPicList, iRefIdxPred) == -1)
		{
			if (pcCU->getSlice()->getSPS()->getPicBg()->getSpliceIdx(pcCU->getCtuRsAddr()) != 0)
			{
				m_iSearchRange = 8;
			}
			else
			{
				m_iSearchRange = 1;
			}
		}
#endif

		TComPattern   tmpPattern;
		TComPattern*  pcPatternKey = &tmpPattern;

		Double        fWeight = 1.0;

		pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
#if VCEG_AZ07_IMV
		pcPatternKey->setImvFlag(pcCU->getiMVFlag(uiPartAddr));
#endif
#if VCEG_AZ06_IC
		Bool bICFlag = pcCU->getICFlag(uiPartAddr);
		pcPatternKey->setMRFlag(bICFlag);
#endif
		if (bBi)
		{
			TComYuv*  pcYuvOther = &m_acYuvPred[1 - (Int)eRefPicList];
			pcYuv = &m_cYuvPredTemp;

			pcYuvOrg->copyPartToPartYuv(pcYuv, uiPartAddr, iRoiWidth, iRoiHeight);

			pcYuv->removeHighFreq(pcYuvOther, uiPartAddr, iRoiWidth, iRoiHeight, pcCU->getSlice()->getSPS()->getBitDepths().recon, m_pcEncCfg->getClipForBiPredMeEnabled());

			fWeight = 0.5;
		}

		//  Search key pattern initialization
		pcPatternKey->initPattern(pcYuv->getAddr(COMPONENT_Y, uiPartAddr),
			iRoiWidth,
			iRoiHeight,
			pcYuv->getStride(COMPONENT_Y),
			pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA));

		Pel*        piRefY = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred)->getPicYuvRec()->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);
		Int         iRefStride = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred)->getPicYuvRec()->getStride(COMPONENT_Y);

		TComMv      cMvPred = *pcMvPred;

#if JVET_C0024_QTBT //for encoder speedup, not try intege ME for one block twice.
		Bool bResetSrchRng = false;
		if (bBi)  xSetSearchRange(pcCU, rcMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
		else
		{
			if (pcCU->getPic()->IsSetIntMv(pcCU->getZorderIdxInCtu(), pcCU->getWidth(0), pcCU->getHeight(0), eRefPicList, iRefIdxPred))
			{
				TComMv cIntMv = pcCU->getPic()->getIntMv(pcCU->getZorderIdxInCtu(), pcCU->getWidth(0), pcCU->getHeight(0), eRefPicList, iRefIdxPred);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				cIntMv <<= (2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
				cIntMv <<= 2;
#endif
				xSetSearchRange(pcCU, cIntMv, 1, cMvSrchRngLT, cMvSrchRngRB);
				bResetSrchRng = true;
			}
			else
			{
				xSetSearchRange(pcCU, cMvPred, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
			}
		}
#else
		if (bBi)
		{
			xSetSearchRange(pcCU, rcMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
		}
		else
		{
			xSetSearchRange(pcCU, cMvPred, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
		}
#endif

		m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(uiPartAddr));

#if JVET_D0123_ME_CTX_LUT_BITS
		pcPatternKey->m_pcCU = pcCU;
		pcPatternKey->m_uiPartAddr = uiPartAddr;
		pcPatternKey->m_mvPred = pcMvPred;
		iCostScale = 2;
#endif

		m_pcRdCost->setPredictor(*pcMvPred);
		m_pcRdCost->setCostScale(2);
#if VCEG_AZ06_IC
		if (bICFlag)
		{
			m_cDistParam.bApplyWeight = false;
		}
		else
#endif
			setWpScalingDistParam(pcCU, iRefIdxPred, eRefPicList);
		//  Do integer search
#if JVET_C0024_QTBT //for encoder speedup, not try intege ME for one block twice.
		if (!m_iFastSearch || bBi || bResetSrchRng)
#else
		if (!m_iFastSearch || bBi)
#endif
		{
			xPatternSearch(pcPatternKey, piRefY, iRefStride, &cMvSrchRngLT, &cMvSrchRngRB, rcMv, ruiCost);
		}
		else
		{
			rcMv = *pcMvPred;
			const TComMv *pIntegerMv2Nx2NPred = 0;
#if !JVET_C0024_QTBT  //not understand the meaning
			if (pcCU->getPartitionSize(0) != SIZE_2Nx2N || pcCU->getDepth(0) != 0)
			{
				pIntegerMv2Nx2NPred = &(m_integerMv2Nx2N[eRefPicList][iRefIdxPred]);
			}
#endif
			xPatternSearchFast(pcCU, pcPatternKey, piRefY, iRefStride, &cMvSrchRngLT, &cMvSrchRngRB, rcMv, ruiCost, pIntegerMv2Nx2NPred);
#if JVET_C0024_QTBT    
			pcCU->getPic()->setIntMv(pcCU->getZorderIdxInCtu(), pcCU->getWidth(0), pcCU->getHeight(0), eRefPicList, iRefIdxPred, rcMv);
#else
			if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
			{
				m_integerMv2Nx2N[eRefPicList][iRefIdxPred] = rcMv;
			}
#endif
		}

		m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(uiPartAddr));
		m_pcRdCost->setCostScale(1);
#if JVET_D0123_ME_CTX_LUT_BITS
		iCostScale = 1;
#endif
		const Bool bIsLosslessCoded = pcCU->getCUTransquantBypass(uiPartAddr) != 0;
#if JVET_E0076_MULTI_PEL_MVD
		if (!intMvRefineNeeded(pcCU, iPartIdx))
		{
			xPatternSearchFracDIF(bIsLosslessCoded, pcPatternKey, piRefY, iRefStride, &rcMv, cMvHalf, cMvQter, ruiCost
#if RRC_SELECTABLE_IF
				, pcCU, pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred), uiPartAddr
#endif
			);
		}
		else
		{
			cMvHalf.set(0, 0);
			cMvQter.set(0, 0);
		}
#else
		xPatternSearchFracDIF(bIsLosslessCoded, pcPatternKey, piRefY, iRefStride, &rcMv, cMvHalf, cMvQter, ruiCost
#if RRC_SELECTABLE_IF
			, pcCU, pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred), uiPartAddr
#endif
		);
#endif

		m_pcRdCost->setCostScale(0);
#if JVET_D0123_ME_CTX_LUT_BITS
		iCostScale = 0;
#endif
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		rcMv <<= (2 + VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
		rcMv += (cMvHalf <<= (VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE + 1));
		rcMv += (cMvQter <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
		rcMv <<= 2;
		rcMv += (cMvHalf <<= 1);
		rcMv += cMvQter;
#endif

#if JVET_E0076_MULTI_PEL_MVD
		if (intMvRefineNeeded(pcCU, iPartIdx))
		{
			ruiBits -= m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];

			intMvRefine(pcCU, eRefPicList, pcPatternKey, piRefY, iRefStride, rcMv, *pcMvPred, riMVPIdx, ruiBits, ruiCost, fWeight);
#if JVET_D0123_ME_CTX_LUT_BITS
			pcPatternKey->m_mvPred = pcMvPred;
#else
			m_pcRdCost->setPredictor(*pcMvPred);
#endif 
		}
#endif
#if JVET_D0123_ME_CTX_LUT_BITS
		UInt uiMvBits = xMvdBits(pcPatternKey, rcMv.getHor(), rcMv.getVer());
#else
		UInt uiMvBits = m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			rcMv.getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, rcMv.getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
			rcMv.getHor(), rcMv.getVer()
#endif
#if VCEG_AZ07_IMV
			, pcCU->getiMVFlag(uiPartAddr)
#endif
		);
#endif

		ruiBits += uiMvBits;
#if JVET_E0076_MULTI_PEL_MVD
		if (!intMvRefineNeeded(pcCU, iPartIdx))
#endif
		{
#if JVET_D0123_ME_CTX_LUT_BITS
			ruiCost = (Distortion)(floor(fWeight * ((Double)ruiCost - (Double)m_pcRdCost->getCostSearch(uiMvBits))) + (Double)m_pcRdCost->getCostSearch(ruiBits));
#else
			ruiCost = (Distortion)(floor(fWeight * ((Double)ruiCost - (Double)m_pcRdCost->getCost(uiMvBits))) + (Double)m_pcRdCost->getCost(ruiBits));
#endif
		}
	}




	Void TEncSearch::xSetSearchRange(TComDataCU* pcCU, TComMv& cMvPred, Int iSrchRng, TComMv& rcMvSrchRngLT, TComMv& rcMvSrchRngRB)
	{
		Int  iMvShift = 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		iMvShift += VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
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




	Void TEncSearch::xPatternSearch(TComPattern* pcPatternKey, Pel* piRefY, Int iRefStride, TComMv* pcMvSrchRngLT, TComMv* pcMvSrchRngRB, TComMv& rcMv, Distortion& ruiSAD)
	{
		Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
		Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
		Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
		Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

		Distortion  uiSad;
		Distortion  uiSadBest = std::numeric_limits<Distortion>::max();
		Int         iBestX = 0;
		Int         iBestY = 0;

		Pel*  piRefSrch;
#if JVET_D0123_ME_CTX_LUT_BITS
		Distortion  uiCurSad;
		Distortion  uiSadOnlyBest = std::numeric_limits<Distortion>::max();
#endif
		//-- jclee for using the SAD function pointer
		m_pcRdCost->setDistParam(pcPatternKey, piRefY, iRefStride, m_cDistParam);

		// fast encoder decision: use subsampled SAD for integer ME
		if (m_pcEncCfg->getUseFastEnc())
		{
			if (m_cDistParam.iRows > 8)
			{
				m_cDistParam.iSubShift = 1;
			}
		}

		piRefY += (iSrchRngVerTop * iRefStride);

		for (Int y = iSrchRngVerTop; y <= iSrchRngVerBottom; y++)
		{
			for (Int x = iSrchRngHorLeft; x <= iSrchRngHorRight; x++)
			{
				//  find min. distortion position
				piRefSrch = piRefY + x;
				m_cDistParam.pCur = piRefSrch;

				setDistParamComp(COMPONENT_Y);
#if VCEG_AZ06_IC
				m_cDistParam.bMRFlag = pcPatternKey->getMRFlag();
#endif
				m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
				uiSad = m_cDistParam.DistFunc(&m_cDistParam);

				// motion cost
#if JVET_D0123_ME_CTX_LUT_BITS
				uiCurSad = uiSad;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				UInt uiMvdBits = xMvdBits(pcPatternKey, x << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, y << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#else
				UInt uiMvdBits = xMvdBits(pcPatternKey, x, y);
#endif
				uiSad += (UInt)m_pcRdCost->getCostSearch(uiMvdBits);
#else
				uiSad += m_pcRdCost->getCost(x, y
#if VCEG_AZ07_IMV
					, pcPatternKey->getImvFlag()
#endif
				);
#endif
				if (uiSad < uiSadBest)
				{
#if JVET_D0123_ME_CTX_LUT_BITS
					uiSadOnlyBest = uiCurSad;
#endif
					uiSadBest = uiSad;
					iBestX = x;
					iBestY = y;
				}
			}
			piRefY += iRefStride;
		}

		rcMv.set(iBestX, iBestY);
#if JVET_D0123_ME_CTX_LUT_BITS
		ruiSAD = uiSadOnlyBest;
#else
		ruiSAD = uiSadBest - m_pcRdCost->getCost(iBestX, iBestY
#if VCEG_AZ07_IMV
			, pcPatternKey->getImvFlag()
#endif
		);
#endif
		return;
	}



	Void TEncSearch::xPatternSearchFast(TComDataCU*   pcCU,
		TComPattern*  pcPatternKey,
		Pel*          piRefY,
		Int           iRefStride,
		TComMv*       pcMvSrchRngLT,
		TComMv*       pcMvSrchRngRB,
		TComMv       &rcMv,
		Distortion   &ruiSAD,
		const TComMv* pIntegerMv2Nx2NPred)
	{
		assert(MD_LEFT < NUM_MV_PREDICTORS);
		pcCU->getMvPredLeft(m_acMvPredictors[MD_LEFT]);
		assert(MD_ABOVE < NUM_MV_PREDICTORS);
		pcCU->getMvPredAbove(m_acMvPredictors[MD_ABOVE]);
		assert(MD_ABOVE_RIGHT < NUM_MV_PREDICTORS);
		pcCU->getMvPredAboveRight(m_acMvPredictors[MD_ABOVE_RIGHT]);

		switch (m_iFastSearch)
		{
		case 1:
			xTZSearch(pcCU, pcPatternKey, piRefY, iRefStride, pcMvSrchRngLT, pcMvSrchRngRB, rcMv, ruiSAD, pIntegerMv2Nx2NPred);
			break;

		case 2:
			xTZSearchSelective(pcCU, pcPatternKey, piRefY, iRefStride, pcMvSrchRngLT, pcMvSrchRngRB, rcMv, ruiSAD, pIntegerMv2Nx2NPred);
			break;
		default:
			break;
		}
	}




	Void TEncSearch::xTZSearch(TComDataCU*  pcCU,
		TComPattern* pcPatternKey,
		Pel*         piRefY,
		Int          iRefStride,
		TComMv*      pcMvSrchRngLT,
		TComMv*      pcMvSrchRngRB,
		TComMv      &rcMv,
		Distortion  &ruiSAD,
		const TComMv* pIntegerMv2Nx2NPred)
	{
		Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
		Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
		Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
		Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

		TZ_SEARCH_CONFIGURATION

			UInt uiSearchRange = m_iSearchRange;
		pcCU->clipMv(rcMv);
		rcMv >>= 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
		rcMv >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
		// init TZSearchStruct
		IntTZSearchStruct cStruct;
		cStruct.iYStride = iRefStride;
		cStruct.piRefY = piRefY;
		cStruct.uiBestSad = MAX_UINT;
#if JVET_D0123_ME_CTX_LUT_BITS
		cStruct.uiBestSadOnly = MAX_UINT;
#endif

		// set rcMv (Median predictor) as start point and as best point
		xTZSearchHelp(pcPatternKey, cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0);

		// test whether one of PRED_A, PRED_B, PRED_C MV is better start point than Median predictor
		if (bTestOtherPredictedMV)
		{
			for (UInt index = 0; index < NUM_MV_PREDICTORS; index++)
			{
				TComMv cMv = m_acMvPredictors[index];
				pcCU->clipMv(cMv);
				cMv >>= 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				cMv >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				xTZSearchHelp(pcPatternKey, cStruct, cMv.getHor(), cMv.getVer(), 0, 0);
			}
		}

		// test whether zero Mv is better start point than Median predictor
		if (bTestZeroVector)
		{
			xTZSearchHelp(pcPatternKey, cStruct, 0, 0, 0, 0);
		}

		if (pIntegerMv2Nx2NPred != 0)
		{
			TComMv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
			integerMv2Nx2NPred <<= 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			integerMv2Nx2NPred <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->clipMv(integerMv2Nx2NPred);
			integerMv2Nx2NPred >>= 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			integerMv2Nx2NPred >>= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			xTZSearchHelp(pcPatternKey, cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);

			// reset search range
			TComMv cMvSrchRngLT;
			TComMv cMvSrchRngRB;
			Int iSrchRng = m_iSearchRange;
			TComMv currBestMv(cStruct.iBestX, cStruct.iBestY);
			currBestMv <<= 2;
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			currBestMv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			xSetSearchRange(pcCU, currBestMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
			iSrchRngHorLeft = cMvSrchRngLT.getHor();
			iSrchRngHorRight = cMvSrchRngRB.getHor();
			iSrchRngVerTop = cMvSrchRngLT.getVer();
			iSrchRngVerBottom = cMvSrchRngRB.getVer();
		}

		// start search
		Int  iDist = 0;
		Int  iStartX = cStruct.iBestX;
		Int  iStartY = cStruct.iBestY;

		// first search
		for (iDist = 1; iDist <= (Int)uiSearchRange; iDist *= 2)
		{
			if (bFirstSearchDiamond == 1)
			{
				xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
			}
			else
			{
				xTZ8PointSquareSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
			}

			if (bFirstSearchStop && (cStruct.uiBestRound >= uiFirstSearchRounds)) // stop criterion
			{
				break;
			}
		}

		// test whether zero Mv is a better start point than Median predictor
		if (bTestZeroVectorStart && ((cStruct.iBestX != 0) || (cStruct.iBestY != 0)))
		{
			xTZSearchHelp(pcPatternKey, cStruct, 0, 0, 0, 0);
			if ((cStruct.iBestX == 0) && (cStruct.iBestY == 0))
			{
				// test its neighborhood
				for (iDist = 1; iDist <= (Int)uiSearchRange; iDist *= 2)
				{
					xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, 0, 0, iDist);
					if (bTestZeroVectorStop && (cStruct.uiBestRound > 0)) // stop criterion
					{
						break;
					}
				}
			}
		}

		// calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
		if (cStruct.uiBestDistance == 1)
		{
			cStruct.uiBestDistance = 0;
			xTZ2PointSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB);
		}

		// raster search if distance is too big
		if (bEnableRasterSearch && (((Int)(cStruct.uiBestDistance) > iRaster) || bAlwaysRasterSearch))
		{
			cStruct.uiBestDistance = iRaster;
			for (iStartY = iSrchRngVerTop; iStartY <= iSrchRngVerBottom; iStartY += iRaster)
			{
				for (iStartX = iSrchRngHorLeft; iStartX <= iSrchRngHorRight; iStartX += iRaster)
				{
					xTZSearchHelp(pcPatternKey, cStruct, iStartX, iStartY, 0, iRaster);
				}
			}
		}

		// raster refinement
		if (bRasterRefinementEnable && cStruct.uiBestDistance > 0)
		{
			while (cStruct.uiBestDistance > 0)
			{
				iStartX = cStruct.iBestX;
				iStartY = cStruct.iBestY;
				if (cStruct.uiBestDistance > 1)
				{
					iDist = cStruct.uiBestDistance >>= 1;
					if (bRasterRefinementDiamond == 1)
					{
						xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
					else
					{
						xTZ8PointSquareSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
				}

				// calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
				if (cStruct.uiBestDistance == 1)
				{
					cStruct.uiBestDistance = 0;
					if (cStruct.ucPointNr != 0)
					{
						xTZ2PointSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB);
					}
				}
			}
		}

		// start refinement
		if (bStarRefinementEnable && cStruct.uiBestDistance > 0)
		{
			while (cStruct.uiBestDistance > 0)
			{
				iStartX = cStruct.iBestX;
				iStartY = cStruct.iBestY;
				cStruct.uiBestDistance = 0;
				cStruct.ucPointNr = 0;
				for (iDist = 1; iDist < (Int)uiSearchRange + 1; iDist *= 2)
				{
					if (bStarRefinementDiamond == 1)
					{
						xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
					else
					{
						xTZ8PointSquareSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
					if (bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds)) // stop criterion
					{
						break;
					}
				}

				// calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
				if (cStruct.uiBestDistance == 1)
				{
					cStruct.uiBestDistance = 0;
					if (cStruct.ucPointNr != 0)
					{
						xTZ2PointSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB);
					}
				}
			}
		}

		// write out best match
		rcMv.set(cStruct.iBestX, cStruct.iBestY);
#if JVET_D0123_ME_CTX_LUT_BITS
		ruiSAD = cStruct.uiBestSadOnly;
#else
		ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCost(cStruct.iBestX, cStruct.iBestY
#if VCEG_AZ07_IMV
			, pcPatternKey->getImvFlag()
#endif
		);
#endif
	}


	Void TEncSearch::xTZSearchSelective(TComDataCU*   pcCU,
		TComPattern*  pcPatternKey,
		Pel*          piRefY,
		Int           iRefStride,
		TComMv*       pcMvSrchRngLT,
		TComMv*       pcMvSrchRngRB,
		TComMv       &rcMv,
		Distortion   &ruiSAD,
		const TComMv* pIntegerMv2Nx2NPred)
	{
		SEL_SEARCH_CONFIGURATION

			Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
		Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
		Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
		Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();
		Int   iFirstSrchRngHorLeft = 0;
		Int   iFirstSrchRngHorRight = 0;
		Int   iFirstSrchRngVerTop = 0;
		Int   iFirstSrchRngVerBottom = 0;
		Int   iStartX = 0;
		Int   iStartY = 0;
		Int   iBestX = 0;
		Int   iBestY = 0;
		Int   iDist = 0;

		pcCU->clipMv(rcMv);
		rcMv >>= 2;
		// init TZSearchStruct
		IntTZSearchStruct cStruct;
		cStruct.iYStride = iRefStride;
		cStruct.piRefY = piRefY;
		cStruct.uiBestSad = MAX_UINT;
		cStruct.iBestX = 0;
		cStruct.iBestY = 0;
#if JVET_D0123_ME_CTX_LUT_BITS
		cStruct.uiBestSadOnly = MAX_UINT;
#endif

		// set rcMv (Median predictor) as start point and as best point
		xTZSearchHelp(pcPatternKey, cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0);

		// test whether one of PRED_A, PRED_B, PRED_C MV is better start point than Median predictor
		if (bTestOtherPredictedMV)
		{
			for (UInt index = 0; index < NUM_MV_PREDICTORS; index++)
			{
				TComMv cMv = m_acMvPredictors[index];
				pcCU->clipMv(cMv);
				cMv >>= 2;
				xTZSearchHelp(pcPatternKey, cStruct, cMv.getHor(), cMv.getVer(), 0, 0);
			}
		}

		// test whether zero Mv is better start point than Median predictor
		if (bTestZeroVector)
		{
			xTZSearchHelp(pcPatternKey, cStruct, 0, 0, 0, 0);
		}

		if (pIntegerMv2Nx2NPred != 0)
		{
			TComMv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
			integerMv2Nx2NPred <<= 2;
			pcCU->clipMv(integerMv2Nx2NPred);
			integerMv2Nx2NPred >>= 2;
			xTZSearchHelp(pcPatternKey, cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);

			// reset search range
			TComMv cMvSrchRngLT;
			TComMv cMvSrchRngRB;
			Int iSrchRng = m_iSearchRange;
			TComMv currBestMv(cStruct.iBestX, cStruct.iBestY);
			currBestMv <<= 2;
			xSetSearchRange(pcCU, currBestMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB);
			iSrchRngHorLeft = cMvSrchRngLT.getHor();
			iSrchRngHorRight = cMvSrchRngRB.getHor();
			iSrchRngVerTop = cMvSrchRngLT.getVer();
			iSrchRngVerBottom = cMvSrchRngRB.getVer();
		}

		// Initial search
		iBestX = cStruct.iBestX;
		iBestY = cStruct.iBestY;
		iFirstSrchRngHorLeft = ((iBestX - uiSearchRangeInitial) > iSrchRngHorLeft) ? (iBestX - uiSearchRangeInitial) : iSrchRngHorLeft;
		iFirstSrchRngVerTop = ((iBestY - uiSearchRangeInitial) > iSrchRngVerTop) ? (iBestY - uiSearchRangeInitial) : iSrchRngVerTop;
		iFirstSrchRngHorRight = ((iBestX + uiSearchRangeInitial) < iSrchRngHorRight) ? (iBestX + uiSearchRangeInitial) : iSrchRngHorRight;
		iFirstSrchRngVerBottom = ((iBestY + uiSearchRangeInitial) < iSrchRngVerBottom) ? (iBestY + uiSearchRangeInitial) : iSrchRngVerBottom;

		for (iStartY = iFirstSrchRngVerTop; iStartY <= iFirstSrchRngVerBottom; iStartY += uiSearchStep)
		{
			for (iStartX = iFirstSrchRngHorLeft; iStartX <= iFirstSrchRngHorRight; iStartX += uiSearchStep)
			{
				xTZSearchHelp(pcPatternKey, cStruct, iStartX, iStartY, 0, 0);
				xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, 1);
				xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, 2);
			}
		}

		Int iMaxMVDistToPred = (abs(cStruct.iBestX - iBestX) > iMVDistThresh || abs(cStruct.iBestY - iBestY) > iMVDistThresh);

		//full search with early exit if MV is distant from predictors
		if (bEnableRasterSearch && (iMaxMVDistToPred || bAlwaysRasterSearch))
		{
			for (iStartY = iSrchRngVerTop; iStartY <= iSrchRngVerBottom; iStartY += 1)
			{
				for (iStartX = iSrchRngHorLeft; iStartX <= iSrchRngHorRight; iStartX += 1)
				{
					xTZSearchHelp(pcPatternKey, cStruct, iStartX, iStartY, 0, 1);
				}
			}
		}
		//Smaller MV, refine around predictor
		else if (bStarRefinementEnable && cStruct.uiBestDistance > 0)
		{
			// start refinement
			while (cStruct.uiBestDistance > 0)
			{
				iStartX = cStruct.iBestX;
				iStartY = cStruct.iBestY;
				cStruct.uiBestDistance = 0;
				cStruct.ucPointNr = 0;
				for (iDist = 1; iDist < (Int)uiSearchRange + 1; iDist *= 2)
				{
					if (bStarRefinementDiamond == 1)
					{
						xTZ8PointDiamondSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
					else
					{
						xTZ8PointSquareSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
					}
					if (bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds)) // stop criterion
					{
						break;
					}
				}

				// calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
				if (cStruct.uiBestDistance == 1)
				{
					cStruct.uiBestDistance = 0;
					if (cStruct.ucPointNr != 0)
					{
						xTZ2PointSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB);
					}
				}
			}
		}

		// write out best match
		rcMv.set(cStruct.iBestX, cStruct.iBestY);
#if JVET_D0123_ME_CTX_LUT_BITS
		ruiSAD = cStruct.uiBestSadOnly;
#else
		ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCost(cStruct.iBestX, cStruct.iBestY);
#endif
	}


	Void TEncSearch::xPatternSearchFracDIF(
		Bool         bIsLosslessCoded,
		TComPattern* pcPatternKey,
		Pel*         piRefY,
		Int          iRefStride,
		TComMv*      pcMvInt,
		TComMv&      rcMvHalf,
		TComMv&      rcMvQter,
		Distortion&  ruiCost
#if RRC_SELECTABLE_IF
		, TComDataCU* pcCU
		, TComPic *pcRefPic
		, UInt uiPartAddr
#endif
	)
	{
		//  Reference pattern initialization (integer scale)
		TComPattern cPatternRoi;
		Int         iOffset = pcMvInt->getHor() + pcMvInt->getVer() * iRefStride;
		cPatternRoi.initPattern(piRefY + iOffset,
			pcPatternKey->getROIYWidth(),
			pcPatternKey->getROIYHeight(),
			iRefStride,
			pcPatternKey->getBitDepthY());

		//  Half-pel refinement
		xExtDIFUpSamplingH(&cPatternRoi);

		rcMvHalf = *pcMvInt;   rcMvHalf <<= 1;    // for mv-cost
		TComMv baseRefMv(0, 0);
		ruiCost = xPatternRefinement(pcPatternKey, baseRefMv, 2, rcMvHalf, !bIsLosslessCoded
#if RRC_SELECTABLE_IF
			, pcCU, pcRefPic, uiPartAddr
#endif
		);

		m_pcRdCost->setCostScale(0);
#if JVET_D0123_ME_CTX_LUT_BITS
		iCostScale = 0;
#endif
#if VCEG_AZ07_IMV
		if (pcPatternKey->getImvFlag())
		{
			rcMvQter = TComMv(0, 0);
			return;
		}
#endif

		xExtDIFUpSamplingQ(&cPatternRoi, rcMvHalf
#if SIF_GAUSS
			, pcCU, pcRefPic, uiPartAddr
#endif
		);
		baseRefMv = rcMvHalf;
		baseRefMv <<= 1;

		rcMvQter = *pcMvInt;   rcMvQter <<= 1;    // for mv-cost
		rcMvQter += rcMvHalf;  rcMvQter <<= 1;
		ruiCost = xPatternRefinement(pcPatternKey, baseRefMv, 1, rcMvQter, !bIsLosslessCoded
#if RRC_SELECTABLE_IF
			, pcCU, pcRefPic, uiPartAddr
#endif
		);
	}


	//! encode residual and calculate rate-distortion for a CU block
	Void TEncSearch::encodeResAndCalcRdInterCU(TComDataCU* pcCU, TComYuv* pcYuvOrg, TComYuv* pcYuvPred,
		TComYuv* pcYuvResi, TComYuv* pcYuvResiBest, TComYuv* pcYuvRec,
		Bool bSkipResidual
#if COM16_C806_EMT
		, Double dBestCost
#endif
		DEBUG_STRING_FN_DECLARE(sDebug))
	{
		assert(!pcCU->isIntra(0));

		const UInt cuWidthPixels = pcCU->getWidth(0);
		const UInt cuHeightPixels = pcCU->getHeight(0);
#if JVET_C0024_QTBT
		UInt uiWIdx = g_aucConvertToBit[cuWidthPixels];
		UInt uiHIdx = g_aucConvertToBit[cuHeightPixels];

#if COM16_C806_EMT
		UInt uiShort;
		UInt uiLong;
		UInt uiStr;
#if TT_TREE_DS1 //data copy - prepare
		//always split into four parts if possible
		UInt uiNumXDiv = (cuWidthPixels  > pcCU->getPic()->getMinCUWidth()) ? 2 : 1;
		UInt uiNumYDiv = (cuHeightPixels > pcCU->getPic()->getMinCUHeight()) ? 2 : 1;
		UInt uiWidth = (uiNumXDiv == 1) ? cuWidthPixels : (cuWidthPixels >> 1);
		UInt uiHeight = (uiNumYDiv == 1) ? cuHeightPixels : (cuHeightPixels >> 1);
		UInt uiLenX = uiWidth / pcCU->getPic()->getMinCUWidth();
		UInt uiLenY = uiHeight / pcCU->getPic()->getMinCUHeight();

		if (uiHeight > uiWidth)
		{
			uiShort = uiWidth;
			uiLong = uiHeight;
			uiStr = pcCU->getPic()->getNumPartInCtuWidth();
		}
		else
		{
			uiShort = uiHeight;
			uiLong = uiWidth;
			uiStr = 1;
		}
#else
		if (cuHeightPixels > cuWidthPixels)
		{
			uiShort = cuWidthPixels;
			uiLong = cuHeightPixels;
			uiStr = pcCU->getPic()->getNumPartInCtuWidth();
		}
		else
		{
			uiShort = cuHeightPixels;
			uiLong = cuWidthPixels;
			uiStr = 1;
		}
#endif
		UInt uiCurrDepth = g_aucConvertToBit[pcCU->getSlice()->getSPS()->getCTUSize()] - g_aucConvertToBit[uiShort];
		UInt uiCurrPartNumb = pcCU->getPic()->getNumPartitionsInCtu() >> (uiCurrDepth << 1);
		UInt uiNumPartInShort = pcCU->getPic()->getNumPartInCtuWidth() >> uiCurrDepth;
#endif
#endif
		const Int  numValidComponents = pcCU->getPic()->getNumberValidComponents();
		const TComSPS &sps = *(pcCU->getSlice()->getSPS());

		// The pcCU is not marked as skip-mode at this point, and its m_pcTrCoeff, m_pcArlCoeff, m_puhCbf, m_puhTrIdx will all be 0.
		// due to prior calls to TComDataCU::initEstData(  );

		if (bSkipResidual) //  No residual coding : SKIP mode
		{
			pcCU->setSkipFlagSubParts(true, 0, pcCU->getDepth(0));

			pcYuvResi->clear();

#if JVET_D0033_ADAPTIVE_CLIPPING // encoder, inter, here copy pred in rec , and skip
			if (g_ClipParam.isActive)
			{
				pcYuvPred->clipToPartYuv(pcYuvRec, 0, pcCU->getSlice()->getSPS()->getBitDepths());
			}
			else
			{
				pcYuvPred->copyToPartYuv(pcYuvRec, 0);
			}
#else
			pcYuvPred->copyToPartYuv(pcYuvRec, 0);
#endif
			Distortion distortion = 0;

			for (Int comp = 0; comp < numValidComponents; comp++)
			{
				const ComponentID compID = ComponentID(comp);
				const UInt csx = pcYuvOrg->getComponentScaleX(compID);
				const UInt csy = pcYuvOrg->getComponentScaleY(compID);
#if WCG_LUMA_DQP_CM_SCALE
				if (m_pcEncCfg->getUseLumaDeltaQp() > 0) {
					UInt           iOrgStrideLuma = pcYuvOrg->getStride(COMPONENT_Y);
					Pel           *piOrgLuma = pcYuvOrg->getAddr(COMPONENT_Y);
					distortion += m_pcRdCost->getDistPart(sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID), pcYuvRec->getStride(compID), pcYuvOrg->getAddr(compID),
						pcYuvOrg->getStride(compID), cuWidthPixels >> csx, cuHeightPixels >> csy, compID, DF_SSE_WTD, piOrgLuma, iOrgStrideLuma);  // use weighted SSE
				}
				else
#endif
					distortion += m_pcRdCost->getDistPart(sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID), pcYuvRec->getStride(compID), pcYuvOrg->getAddr(compID),
						pcYuvOrg->getStride(compID), cuWidthPixels >> csx, cuHeightPixels >> csy, compID);
			}
#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_CURR_BEST]);
#endif
			m_pcEntropyCoder->resetBits();

			if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
			{
				m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
			}

			m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL
			if (pcCU->isCodeBIF(pcCU->getWidth(0), pcCU->getDepth(0)
#if JVET_C0024_QTBT && SIF_QTBT_FIX
				+ pcCU->getBTDepth(0) // check total depth
#endif
			))
			{
				m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
			}
#endif
#if VCEG_AZ07_FRUC_MERGE
			m_pcEntropyCoder->encodeFRUCMgrMode(pcCU, 0, 0);
			if (!pcCU->getFRUCMgrMode(0))
#endif
#if MRG_CONSIDER_COMPLEX
			{
				if ((pcCU->getWidth(0) >= 16) || (pcCU->getHeight(0) >= 16))
				{
					m_pcEntropyCoder->encodeMergeComplexFlag(pcCU, 0);
					if (pcCU->getMergeComplexFlag(0))
					{
						m_pcEntropyCoder->encodeMergeComplexIndex(pcCU, 0, true);
					}
				}

				if (!pcCU->getMergeComplexFlag(0))
#endif
#if COM16_C1016_AFFINE
				{
					if (pcCU->isAffineMrgFlagCoded(0, 0))
					{
						m_pcEntropyCoder->encodeAffineFlag(pcCU, 0, 0);
					}

					if (!pcCU->isAffine(0))
					{
						m_pcEntropyCoder->encodeMergeIndex(pcCU, 0, true);
					}
				}
#else
					m_pcEntropyCoder->encodeMergeIndex(pcCU, 0, true);
#endif
#if MRG_CONSIDER_COMPLEX
			}
#endif
#if SIF_CU_LEVEL
			m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
#endif

			UInt uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
			pcCU->getTotalBits() = uiBits;
			pcCU->getTotalDistortion() = distortion;
			pcCU->getTotalCost() = m_pcRdCost->calcRdCost(uiBits, distortion);

#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
			m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_TEMP_BEST]);
#endif

#if DEBUG_STRING
			pcYuvResiBest->clear(); // Clear the residual image, if we didn't code it.
			for (UInt i = 0; i<MAX_NUM_COMPONENT + 1; i++)
			{
				sDebug += debug_reorder_data_inter_token[i];
			}
#endif

			return;
		}

		//  Residual coding.
#if JVET_C0024_QTBT
		pcYuvResi->subtract(pcYuvOrg, pcYuvPred, 0, cuWidthPixels, cuHeightPixels);
#else
		pcYuvResi->subtract(pcYuvOrg, pcYuvPred, 0, cuWidthPixels);
#endif
#if JVET_D0033_ADAPTIVE_CLIPPING_ENC_METHOD
		smoothResidual(pcYuvOrg, pcYuvResi, 0, cuHeightPixels, cuWidthPixels);
#endif

#if COM16_C806_EMT
		Bool bestIsSkip = false;
		UChar ucBestCuFlag = 0;
#if JVET_C0024_QTBT
		UChar ucEmtUsage = ((cuWidthPixels> EMT_INTER_MAX_CU || cuHeightPixels > EMT_INTER_MAX_CU) || (pcCU->getSlice()->getSPS()->getUseInterEMT() == 0)) ? 1 : 2;
#else
		UChar ucEmtUsage = ((cuWidthPixels > EMT_INTER_MAX_CU) || (pcCU->getSlice()->getSPS()->getUseInterEMT() == 0)) ? 1 : 2;
#endif

#if COM16_C806_OBMC 
		// Disable EMT when CU-level OBMC control flag is 0
		if (m_pcEncCfg->getUseFastInterEMT() && pcCU->getSlice()->getSPS()->getOBMC() && pcCU->getOBMCFlag(0) == 0)
		{
			ucEmtUsage = 1;
		}
#endif
		Bool  bZeroCu = false;
		UInt bestBits = 0;
		Double dCost = 0, bestCost = MAX_DOUBLE;

		// CU-level optimization
		for (UChar ucCuFlag = 0; ucCuFlag < ucEmtUsage; ucCuFlag++)
		{
			if (ucCuFlag && m_pcEncCfg->getUseFastInterEMT())
			{
				static const Double thEmtInterSkip = 1.1; // Skip checking EMT transforms
				if (dCost > thEmtInterSkip * dBestCost)
				{
					break;
				}
			}
			if (ucCuFlag && bZeroCu && m_pcEncCfg->getUseFastInterEMT())
			{
				break;
			}
			pcCU->setEmtCuFlagSubParts(ucCuFlag, 0, pcCU->getDepth(0));
			pcCU->setSkipFlagSubParts(false, 0, pcCU->getDepth(0));
#endif

			TComTURecurse tuLevel0(pcCU, 0);

			Double     nonZeroCost = 0;
			UInt       nonZeroBits = 0;
			Distortion nonZeroDistortion = 0;
			Distortion zeroDistortion = 0;

#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_CURR_BEST]);
#endif

#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
			xEstimateInterResidualQT(pcYuvResi, nonZeroCost, nonZeroBits, nonZeroDistortion, &zeroDistortion, tuLevel0, pcYuvPred DEBUG_STRING_PASS_INTO(sDebug));
#else
			xEstimateInterResidualQT(pcYuvResi, nonZeroCost, nonZeroBits, nonZeroDistortion, &zeroDistortion, tuLevel0 DEBUG_STRING_PASS_INTO(sDebug));
#endif
			// -------------------------------------------------------
			// set the coefficients in the pcCU, and also calculates the residual data.
			// If a block full of 0's is efficient, then just use 0's.
			// The costs at this point do not include header bits.

			m_pcEntropyCoder->resetBits();
			m_pcEntropyCoder->encodeQtRootCbfZero();
			const UInt   zeroResiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
#if WCG_LUMA_DQP_CM_SCALE
			Double zeroCost;
			if (m_pcEncCfg->getUseLumaDeltaQp() > 0)
			{
				zeroCost = (pcCU->isLosslessCoded(0)) ? (nonZeroCost + 1) : (m_pcRdCost->calcRdCost(zeroResiBits, zeroDistortion, false, DF_DEFAULT_ORI));
			}
			else {
				zeroCost = (pcCU->isLosslessCoded(0)) ? (nonZeroCost + 1) : (m_pcRdCost->calcRdCost(zeroResiBits, zeroDistortion));
			}
#else
			const Double zeroCost = (pcCU->isLosslessCoded(0)) ? (nonZeroCost + 1) : (m_pcRdCost->calcRdCost(zeroResiBits, zeroDistortion));
#endif

			if (zeroCost < nonZeroCost || !pcCU->getQtRootCbf(0))
			{
#if JVET_C0024_QTBT
				for (Int comp = 0; comp < numValidComponents; comp++)
				{
					const ComponentID component = ComponentID(comp);
					pcCU->setCbfSubParts(0, component, 0, 0);
					pcCU->setCrossComponentPredictionAlphaPartRange(0, component, 0, 0);
				}
#else
				const UInt uiQPartNum = tuLevel0.GetAbsPartIdxNumParts();
				::memset(pcCU->getTransformIdx(), 0, uiQPartNum * sizeof(UChar));
				for (Int comp = 0; comp < numValidComponents; comp++)
				{
					const ComponentID component = ComponentID(comp);
					::memset(pcCU->getCbf(component), 0, uiQPartNum * sizeof(UChar));
					::memset(pcCU->getCrossComponentPredictionAlpha(component), 0, (uiQPartNum * sizeof(Char)));
				}
#endif
				static const UInt useTS[MAX_NUM_COMPONENT] = { 0,0,0 };
				pcCU->setTransformSkipSubParts(useTS, 0, pcCU->getDepth(0));
#if DEBUG_STRING
				sDebug.clear();
				for (UInt i = 0; i<MAX_NUM_COMPONENT + 1; i++)
				{
					sDebug += debug_reorder_data_inter_token[i];
				}
#endif
			}
			else
			{
				xSetInterResidualQTData(NULL, false, tuLevel0); // Call first time to set coefficients.
			}

			// all decisions now made. Fully encode the CU, including the headers:
#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST]);
#else
			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_CURR_BEST]);
#endif

			UInt finalBits = 0;
			xAddSymbolBitsInter(pcCU, finalBits);
			// we've now encoded the pcCU, and so have a valid bit cost

#if COM16_C806_EMT
			if (m_pcEncCfg->getUseFastInterEMT())
			{
				bZeroCu = pcCU->getCbf(0, COMPONENT_Y) ? false : true;
			}

			Distortion uiDist = (zeroCost < nonZeroCost || !pcCU->getQtRootCbf(0)) ? zeroDistortion : nonZeroDistortion;
#if WCG_LUMA_DQP_CM_SCALE
			if (m_pcEncCfg->getUseLumaDeltaQp() > 0)
			{
				dCost = m_pcRdCost->calcRdCost(finalBits, uiDist, false, DF_DEFAULT_ORI);
			}
			else
#endif
				dCost = m_pcRdCost->calcRdCost(finalBits, uiDist);
			if (dCost < bestCost)
			{
				bestBits = finalBits;
				bestCost = dCost;
				ucBestCuFlag = ucCuFlag;
#if JVET_C0024_QTBT
				bestIsSkip = (pcCU->getMergeFlag(0) && !pcCU->getQtRootCbf(0));
#else
				bestIsSkip = (pcCU->getMergeFlag(0) && pcCU->getPartitionSize(0) == SIZE_2Nx2N && !pcCU->getQtRootCbf(0));
#endif
#endif

				if (!pcCU->getQtRootCbf(0))
				{
					pcYuvResiBest->clear(); // Clear the residual image, if we didn't code it.
				}
				else
				{
					xSetInterResidualQTData(pcYuvResiBest, true, tuLevel0); // else set the residual image data pcYUVResiBest from the various temp images.
				}
#if JVET_C0024_QTBT
				m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#else
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_TEMP_BEST]);
#endif

#if COM16_C806_EMT
				if (ucCuFlag < (ucEmtUsage - 1))
				{
					// Store current results to the temporary buffer
#if JVET_C0024_QTBT
#if TT_TREE_DS1 //data copy
					for (UInt yIdx = 0; yIdx < uiNumYDiv; yIdx++)
					{
						for (UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
						{
							UInt uiRaster = (uiLenX * xIdx) + (uiLenY * yIdx) * pcCU->getPic()->getNumPartInCtuWidth();
#else
					UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
#endif
					for (UInt j = 0; j<uiLong; j += uiShort)
					{
#if TT_TREE_DS1
						UInt uiZorder = g_auiRasterToZscan[uiRaster];
#else
						UInt uiZorder = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
#endif

						for (Int componentIndex = 0; componentIndex<pcCU->getPic()->getNumberValidComponents(); componentIndex++)
						{
							const ComponentID compID = ComponentID(componentIndex);
							::memcpy(m_puhQTTempCbf[compID] + uiZorder, pcCU->getCbf(compID) + uiZorder, uiCurrPartNumb);
							::memcpy(m_puhQTTempExplicitRdpcmMode[compID] + uiZorder, pcCU->getExplicitRdpcmMode(compID) + uiZorder, uiCurrPartNumb * sizeof(UChar));
							::memcpy(m_puhQTTempTransformSkipFlag[compID] + uiZorder, pcCU->getTransformSkip(compID) + uiZorder, uiCurrPartNumb * sizeof(UChar));
							::memcpy(m_phQTTempCrossComponentPredictionAlpha[compID] + uiZorder, pcCU->getCrossComponentPredictionAlpha(compID) + uiZorder, uiCurrPartNumb * sizeof(Char));
						}
						uiRaster += uiNumPartInShort * uiStr;
					}
#if TT_TREE_DS1 //data copy end
						}
					}
#endif
			const UInt numCoeffY = cuWidthPixels * cuHeightPixels;
			for (Int componentIndex = 0; componentIndex<pcCU->getPic()->getNumberValidComponents(); componentIndex++)
			{
				const ComponentID compID = ComponentID(componentIndex);
				const UInt componentShift = pcCU->getPic()->getComponentScaleX(compID) + pcCU->getPic()->getComponentScaleY(compID);
				::memcpy(m_pcQTTempCoeff[compID], pcCU->getCoeff(compID), (numCoeffY >> componentShift) * sizeof(TCoeff));
#if ADAPTIVE_QP_SELECTION
				::memcpy(m_pcQTTempArlCoeff[compID], pcCU->getArlCoeff(compID), (numCoeffY >> componentShift) * sizeof(TCoeff));
#endif
			}
#else
					const UInt uiQPartNum = pcCU->getPic()->getNumPartitionsInCtu() >> (pcCU->getDepth(0) << 1);
					::memcpy(m_puhQTTempTrIdx, pcCU->getTransformIdx(), uiQPartNum * sizeof(UChar));
					const UInt numCoeffY = cuWidthPixels * cuHeightPixels;
					for (Int i = 0; i<pcCU->getPic()->getNumberValidComponents(); i++)
					{
						const ComponentID compID = ComponentID(i);
						const UInt componentShift = pcCU->getPic()->getComponentScaleX(compID) + pcCU->getPic()->getComponentScaleY(compID);
						::memcpy(m_puhQTTempCbf[compID], pcCU->getCbf(compID), uiQPartNum * sizeof(UChar));
						::memcpy(m_pcQTTempCoeff[compID], pcCU->getCoeff(compID), (numCoeffY >> componentShift) * sizeof(TCoeff));
#if ADAPTIVE_QP_SELECTION
						::memcpy(m_pcQTTempArlCoeff[compID], pcCU->getArlCoeff(compID), (numCoeffY >> componentShift) * sizeof(TCoeff));
#endif
						::memcpy(m_puhQTTempExplicitRdpcmMode[compID], pcCU->getExplicitRdpcmMode(compID), uiQPartNum * sizeof(UChar));
						::memcpy(m_puhQTTempTransformSkipFlag[compID], pcCU->getTransformSkip(compID), uiQPartNum * sizeof(UChar));
#if VCEG_AZ08_KLT_COMMON 
						::memcpy(m_puhQTTempKLTFlag[compID], pcCU->getKLTFlag(compID), uiQPartNum * sizeof(UChar));
#endif
						::memcpy(m_phQTTempCrossComponentPredictionAlpha[compID], pcCU->getCrossComponentPredictionAlpha(compID), uiQPartNum * sizeof(Char));
					}
#endif
				}
			}
		}
#endif

#if COM16_C806_EMT
		pcCU->setEmtCuFlagSubParts(ucBestCuFlag, 0, pcCU->getDepth(0));
		pcCU->setSkipFlagSubParts(bestIsSkip, 0, pcCU->getDepth(0));
		if (ucBestCuFlag < (ucEmtUsage - 1))
		{
			// Restore the best results from the temporary buffer
#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_TEMP_BEST]);
#if TT_TREE_DS1 //data copy
			for (UInt yIdx = 0; yIdx < uiNumYDiv; yIdx++)
			{
				for (UInt xIdx = 0; xIdx < uiNumXDiv; xIdx++)
				{
					UInt uiRaster = (uiLenX * xIdx) + (uiLenY * yIdx) * pcCU->getPic()->getNumPartInCtuWidth();
#else
			UInt uiRaster = g_auiZscanToRaster[pcCU->getZorderIdxInCtu()];
#endif
			for (UInt j = 0; j<uiLong; j += uiShort)
			{
#if TT_TREE_DS1
				UInt uiZorder = g_auiRasterToZscan[uiRaster];
#else
				UInt uiZorder = g_auiRasterToZscan[uiRaster] - pcCU->getZorderIdxInCtu();
#endif

				for (Int componentIndex = 0; componentIndex<pcCU->getPic()->getNumberValidComponents(); componentIndex++)
				{
					const ComponentID compID = ComponentID(componentIndex);
					::memcpy(pcCU->getCbf(compID) + uiZorder, m_puhQTTempCbf[compID] + uiZorder, uiCurrPartNumb);
					::memcpy(pcCU->getExplicitRdpcmMode(compID) + uiZorder, m_puhQTTempExplicitRdpcmMode[compID] + uiZorder, uiCurrPartNumb * sizeof(UChar));
					::memcpy(pcCU->getTransformSkip(compID) + uiZorder, m_puhQTTempTransformSkipFlag[compID] + uiZorder, uiCurrPartNumb * sizeof(UChar));
					::memcpy(pcCU->getCrossComponentPredictionAlpha(compID) + uiZorder, m_phQTTempCrossComponentPredictionAlpha[compID] + uiZorder, uiCurrPartNumb * sizeof(Char));
				}
				uiRaster += uiNumPartInShort * uiStr;
			}
#if TT_TREE_DS1 //data copy end
				}
			}
#endif
	const UInt numCoeffY = cuWidthPixels * cuHeightPixels;
	for (Int componentIndex = 0; componentIndex<pcCU->getPic()->getNumberValidComponents(); componentIndex++)
	{
		const ComponentID compID = ComponentID(componentIndex);
		const UInt componentShift = pcCU->getPic()->getComponentScaleX(compID) + pcCU->getPic()->getComponentScaleY(compID);
		::memcpy(pcCU->getCoeff(compID), m_pcQTTempCoeff[compID], (numCoeffY >> componentShift) * sizeof(TCoeff));
#if ADAPTIVE_QP_SELECTION
		::memcpy(pcCU->getArlCoeff(compID), m_pcQTTempArlCoeff[compID], (numCoeffY >> componentShift) * sizeof(TCoeff));
#endif
	}
#else
			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_TEMP_BEST]);
			const UInt uiQPartNum = pcCU->getPic()->getNumPartitionsInCtu() >> (pcCU->getDepth(0) << 1);
			::memcpy(pcCU->getTransformIdx(), m_puhQTTempTrIdx, uiQPartNum * sizeof(UChar));
			const UInt numCoeffY = cuWidthPixels * cuHeightPixels;
			for (Int i = 0; i<pcCU->getPic()->getNumberValidComponents(); i++)
			{
				const ComponentID compID = ComponentID(i);
				const UInt componentShift = pcCU->getPic()->getComponentScaleX(compID) + pcCU->getPic()->getComponentScaleY(compID);
				::memcpy(pcCU->getCbf(compID), m_puhQTTempCbf[compID], uiQPartNum * sizeof(UChar));
				::memcpy(pcCU->getCoeff(compID), m_pcQTTempCoeff[compID], (numCoeffY >> componentShift) * sizeof(TCoeff));
#if ADAPTIVE_QP_SELECTION
				::memcpy(pcCU->getArlCoeff(compID), m_pcQTTempArlCoeff[compID], (numCoeffY >> componentShift) * sizeof(TCoeff));
#endif
				::memcpy(pcCU->getExplicitRdpcmMode(compID), m_puhQTTempExplicitRdpcmMode[compID], uiQPartNum * sizeof(UChar));
				::memcpy(pcCU->getTransformSkip(compID), m_puhQTTempTransformSkipFlag[compID], uiQPartNum * sizeof(UChar));
#if VCEG_AZ08_KLT_COMMON
				::memcpy(pcCU->getKLTFlag(compID), m_puhQTTempKLTFlag[compID], uiQPartNum * sizeof(UChar));
#endif
				::memcpy(pcCU->getCrossComponentPredictionAlpha(compID), m_phQTTempCrossComponentPredictionAlpha[compID], uiQPartNum * sizeof(Char));
			}
#endif
		}
#endif

#if JVET_C0024_QTBT
		pcYuvRec->addClip(pcYuvPred, pcYuvResiBest, 0, cuWidthPixels, cuHeightPixels, sps.getBitDepths());
#else
		pcYuvRec->addClip(pcYuvPred, pcYuvResiBest, 0, cuWidthPixels, sps.getBitDepths());
#endif

		// update with clipped distortion and cost (previously unclipped reconstruction values were used)

		Distortion finalDistortion = 0;
		for (Int comp = 0; comp<numValidComponents; comp++)
		{
			const ComponentID compID = ComponentID(comp);
#if WCG_LUMA_DQP_CM_SCALE
			if (m_pcEncCfg->getUseLumaDeltaQp() > 0) {
				UInt           iOrgStrideLuma = pcYuvOrg->getStride(COMPONENT_Y);
				Pel           *piOrgLuma = pcYuvOrg->getAddr(COMPONENT_Y);
				finalDistortion += m_pcRdCost->getDistPart(sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID), pcYuvRec->getStride(compID), pcYuvOrg->getAddr(compID),
					pcYuvOrg->getStride(compID), cuWidthPixels >> pcYuvOrg->getComponentScaleX(compID), cuHeightPixels >> pcYuvOrg->getComponentScaleY(compID), compID, DF_SSE_WTD
					, piOrgLuma, iOrgStrideLuma
				);  // use weighted SSE
			}
			else
#endif
				finalDistortion += m_pcRdCost->getDistPart(sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID), pcYuvRec->getStride(compID), pcYuvOrg->getAddr(compID), pcYuvOrg->getStride(compID), cuWidthPixels >> pcYuvOrg->getComponentScaleX(compID), cuHeightPixels >> pcYuvOrg->getComponentScaleY(compID), compID);
		}

#if COM16_C806_EMT
		pcCU->getTotalBits() = bestBits;
		pcCU->getTotalDistortion() = finalDistortion;
		pcCU->getTotalCost() = m_pcRdCost->calcRdCost(bestBits, finalDistortion);
#else
		pcCU->getTotalBits() = finalBits;
		pcCU->getTotalDistortion() = finalDistortion;
		pcCU->getTotalCost() = m_pcRdCost->calcRdCost(finalBits, finalDistortion);
#endif
	}



	Void TEncSearch::xEstimateInterResidualQT(TComYuv    *pcResi,
		Double     &rdCost,
		UInt       &ruiBits,
		Distortion &ruiDist,
		Distortion *puiZeroDist,
		TComTU     &rTu
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
		, TComYuv* pcPred
#endif
		DEBUG_STRING_FN_DECLARE(sDebug)
	)
	{
		TComDataCU *pcCU = rTu.getCU();
		const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
		const UInt uiDepth = rTu.GetTransformDepthTotal();
		const UInt uiTrMode = rTu.GetTransformDepthRel();
#if JVET_C0024_QTBT
		assert(uiAbsPartIdx == 0 && uiTrMode == 0);
		UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
		UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
		const UInt subTUDepth = uiTrMode + 1;
#endif
		const UInt numValidComp = pcCU->getPic()->getNumberValidComponents();
		DEBUG_STRING_NEW(sSingleStringComp[MAX_NUM_COMPONENT])

			assert(pcCU->getDepth(0) == pcCU->getDepth(uiAbsPartIdx));
#if JVET_C0024_QTBT
		Bool bCheckFull = true;
#else
		const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();

		UInt SplitFlag = ((pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) && pcCU->isInter(uiAbsPartIdx) && (pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N));
#if DEBUG_STRING
		const Int debugPredModeMask = DebugStringGetPredModeMask(pcCU->getPredictionMode(uiAbsPartIdx));
#endif

		Bool bCheckFull;

		if (SplitFlag && uiDepth == pcCU->getDepth(uiAbsPartIdx) && (uiLog2TrSize >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx)))
		{
			bCheckFull = false;
		}
		else
		{
			bCheckFull = (uiLog2TrSize <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize());
		}

		const Bool bCheckSplit = (uiLog2TrSize >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx));

		assert(bCheckFull || bCheckSplit);
#endif

		// code full block
		Double     dSingleCost = MAX_DOUBLE;
		UInt       uiSingleBits = 0;
		Distortion uiSingleDistComp[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };
		Distortion uiSingleDist = 0;
		TCoeff     uiAbsSum[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };
		UInt       uiBestTransformMode[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };
		//  Stores the best explicit RDPCM mode for a TU encoded without split
		UInt       bestExplicitRdpcmModeUnSplit[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 3,3 },{ 3,3 },{ 3,3 } };
		Char       bestCrossCPredictionAlpha[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };
#if COM16_C806_EMT
		UChar      bestEmtTrIdx[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };
#endif
#if VCEG_AZ08_INTER_KLT
		UInt       bestKLTMode[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = { { 0,0 },{ 0,0 },{ 0,0 } };;
#endif
#if JVET_C0024_QTBT
		m_pcRDGoOnSbacCoder->store(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
		m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
#endif

#if COM16_C806_EMT
		pcCU->setEmtTuIdxSubParts(0, uiAbsPartIdx, uiDepth); // initialization
#endif

		if (bCheckFull)
		{
			Double minCost[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/];
			Bool checkTransformSkip[MAX_NUM_COMPONENT];
#if !JVET_C0024_QTBT
			pcCU->setTrIdxSubParts(uiTrMode, uiAbsPartIdx, uiDepth);
#endif

			m_pcEntropyCoder->resetBits();

			memset(m_pTempPel, 0, sizeof(Pel) * rTu.getRect(COMPONENT_Y).width * rTu.getRect(COMPONENT_Y).height); // not necessary needed for inside of recursion (only at the beginning)

#if !JVET_C0024_QTBT
			const UInt uiQTTempAccessLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif
			TCoeff *pcCoeffCurr[MAX_NUM_COMPONENT];
#if ADAPTIVE_QP_SELECTION
			TCoeff *pcArlCoeffCurr[MAX_NUM_COMPONENT];
#endif

			for (UInt i = 0; i<numValidComp; i++)
			{
				minCost[i][0] = MAX_DOUBLE;
				minCost[i][1] = MAX_DOUBLE;
			}

#if VCEG_AZ08_INTER_KLT
			Bool checkKLT[MAX_NUM_COMPONENT] = { 0, 0, 0 };
#endif
			Pel crossCPredictedResidualBuffer[MAX_TU_SIZE * MAX_TU_SIZE];

			for (UInt i = 0; i<numValidComp; i++)
			{
				checkTransformSkip[i] = false;
				const ComponentID compID = ComponentID(i);
				const Int channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#if JVET_C0024_QTBT
				assert(rTu.getCoefficientOffset(compID) == 0);
				pcCoeffCurr[compID] = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#if ADAPTIVE_QP_SELECTION
				pcArlCoeffCurr[compID] = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx] + rTu.getCoefficientOffset(compID);
#endif
#else
				pcCoeffCurr[compID] = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + rTu.getCoefficientOffset(compID);
#if ADAPTIVE_QP_SELECTION
				pcArlCoeffCurr[compID] = m_ppcQTTempArlCoeff[compID][uiQTTempAccessLayer] + rTu.getCoefficientOffset(compID);
#endif
#endif

				if (rTu.ProcessComponentSection(compID))
				{
					const QpParam cQP(*pcCU, compID);

					checkTransformSkip[compID] = pcCU->getSlice()->getPPS()->getUseTransformSkip() &&
#if JVET_C0024_QTBT
						TUCompRectHasAssociatedTransformSkipFlag(pcCU->getSlice()->isIntra(), rTu.getRect(compID), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize()) &&
#else
						TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(compID), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize()) &&
#endif
						(!pcCU->isLosslessCoded(0));
#if JVET_F0031_RMV_REDUNDANT_TRSKIP && COM16_C806_EMT
					if (isLuma(compID))
					{
						checkTransformSkip[compID] &= (!pcCU->getEmtCuFlag(uiAbsPartIdx));
					}
#endif
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER
					UInt tuWidth = rTu.getRect(compID).width;
					UInt tuHeight = rTu.getRect(compID).height;
#endif
#if VCEG_AZ08_INTER_KLT
					UInt uiMaxTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MAX - 1];
					UInt uiMinTrWidth = g_uiDepth2Width[USE_MORE_BLOCKSIZE_DEPTH_MIN - 1];
					checkKLT[compID] = g_bEnableCheck && isLuma(compID) && ((tuWidth == tuHeight) && (tuWidth <= uiMaxTrWidth) && (tuWidth >= uiMinTrWidth));
#endif
#if JVET_C0024_QTBT
					const Bool splitIntoSubTUs = false;
#else
					const Bool splitIntoSubTUs = rTu.getRect(compID).width != rTu.getRect(compID).height;
#endif

					TComTURecurse TUIterator(rTu, false, (splitIntoSubTUs ? TComTU::VERTICAL_SPLIT : TComTU::DONT_SPLIT), true, compID);

#if JVET_C0024_QTBT
					const UInt partIdxesPerSubTU = 0;
#else
					const UInt partIdxesPerSubTU = TUIterator.GetAbsPartIdxNumParts(compID);
#endif

					do
					{
						const UInt           subTUIndex = TUIterator.GetSectionNumber();
						const UInt           subTUAbsPartIdx = TUIterator.GetAbsPartIdxTU(compID);
						const TComRectangle &tuCompRect = TUIterator.getRect(compID);
						const UInt           subTUBufferOffset = tuCompRect.width * tuCompRect.height * subTUIndex;
#if JVET_C0024_QTBT
						assert(subTUBufferOffset == 0);
#endif

						TCoeff        *currentCoefficients = pcCoeffCurr[compID] + subTUBufferOffset;
#if ADAPTIVE_QP_SELECTION
						TCoeff        *currentARLCoefficients = pcArlCoeffCurr[compID] + subTUBufferOffset;
#endif
						const Bool isCrossCPredictionAvailable = isChroma(compID)
							&& pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
							&& (pcCU->getCbf(subTUAbsPartIdx, COMPONENT_Y, uiTrMode) != 0);

						Char preCalcAlpha = 0;
#if JVET_C0024_QTBT
						assert(rTu.getRect(COMPONENT_Y).x0 == 0 && rTu.getRect(COMPONENT_Y).y0 == 0);
						const Pel *pLumaResi = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(COMPONENT_Y, rTu.getRect(COMPONENT_Y).x0, rTu.getRect(COMPONENT_Y).y0);
#else
						const Pel *pLumaResi = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(COMPONENT_Y, rTu.getRect(COMPONENT_Y).x0, rTu.getRect(COMPONENT_Y).y0);
#endif

						if (isCrossCPredictionAvailable)
						{
							const Bool bUseReconstructedResidualForEstimate = m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate();
							const Pel  *const lumaResidualForEstimate = bUseReconstructedResidualForEstimate ? pLumaResi : pcResi->getAddrPix(COMPONENT_Y, tuCompRect.x0, tuCompRect.y0);
#if JVET_C0024_QTBT
							const UInt        lumaResidualStrideForEstimate = bUseReconstructedResidualForEstimate ? m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(COMPONENT_Y) : pcResi->getStride(COMPONENT_Y);
#else
							const UInt        lumaResidualStrideForEstimate = bUseReconstructedResidualForEstimate ? m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y) : pcResi->getStride(COMPONENT_Y);
#endif

							preCalcAlpha = xCalcCrossComponentPredictionAlpha(TUIterator,
								compID,
								lumaResidualForEstimate,
								pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
								tuCompRect.width,
								tuCompRect.height,
								lumaResidualStrideForEstimate,
								pcResi->getStride(compID));
						}

						const Int transformSkipModesToTest = checkTransformSkip[compID] ? 2 : 1;
						const Int crossCPredictionModesToTest = (preCalcAlpha != 0) ? 2 : 1; // preCalcAlpha cannot be anything other than 0 if isCrossCPredictionAvailable is false
#if VCEG_AZ08_INTER_KLT
						const Int KLTModesToTest = checkKLT[compID] ? 2 : 1;
#endif
#if COM16_C806_EMT
#if JVET_C0024_QTBT
#if FAST_EMT_v1
						const Int emtTrIdxToTest = isLuma(compID) && pcCU->getEmtCuFlag(subTUAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseInterEMT() ? 3 : 1;
#else
						const Int emtTrIdxToTest = isLuma(compID) && pcCU->getEmtCuFlag(subTUAbsPartIdx) && pcCU->getSlice()->getSPS()->getUseInterEMT() ? 4 : 1;
#endif
#else
						const Int emtTrIdxToTest = isLuma(compID) && pcCU->getEmtCuFlag(subTUAbsPartIdx) ? 4 : 1;
#endif
#endif
						const Bool isOneMode = (crossCPredictionModesToTest == 1) && (transformSkipModesToTest == 1)
#if COM16_C806_EMT
							&& (emtTrIdxToTest == 1)
#endif
#if VCEG_AZ08_INTER_KLT
							&& (KLTModesToTest == 1)
#endif
							;

#if VCEG_AZ08_INTER_KLT
						for (Int KLTModeId = 0; KLTModeId < KLTModesToTest; KLTModeId++)
						{
							pcCU->setKLTPartRange(KLTModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
#if VCEG_AZ08_INTER_KLT
							for (Int transformSkipModeId = 0; transformSkipModeId < (KLTModeId ? 1 : transformSkipModesToTest); transformSkipModeId++)
#else
							for (Int transformSkipModeId = 0; transformSkipModeId < transformSkipModesToTest; transformSkipModeId++)
#endif
							{
								//pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);
#if COM16_C806_EMT
#if VCEG_AZ08_INTER_KLT
								for (Int emtTrIdx = 0; emtTrIdx < ((transformSkipModeId || KLTModeId) ? 1 : emtTrIdxToTest); emtTrIdx++)
#else
								for (Int emtTrIdx = 0; emtTrIdx < (transformSkipModeId ? 1 : emtTrIdxToTest); emtTrIdx++)
#endif
								{
									pcCU->setEmtTuIdxPartsRange(emtTrIdx, compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
#if VCEG_AZ08_INTER_KLT
									for (Int crossCPredictionModeId = 0; crossCPredictionModeId < (KLTModeId ? 1 : crossCPredictionModesToTest); crossCPredictionModeId++)
#else
									for (Int crossCPredictionModeId = 0; crossCPredictionModeId < crossCPredictionModesToTest; crossCPredictionModeId++)
#endif
									{
										const Bool isFirstMode = (transformSkipModeId == 0) && (crossCPredictionModeId == 0)
#if COM16_C806_EMT
											&& (emtTrIdx == 0)
#endif
#if VCEG_AZ08_INTER_KLT
											&& (KLTModeId == 0)
#endif
											;
										const Bool bUseCrossCPrediction = crossCPredictionModeId != 0;
#if VCEG_AZ08_INTER_KLT
										const Bool bCheckKLT = (KLTModeId != 0);
										Bool useKLT = false;
#endif

#if JVET_C0024_QTBT
										m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
										m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
#endif
										m_pcEntropyCoder->resetBits();

										pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);
										pcCU->setCrossComponentPredictionAlphaPartRange((bUseCrossCPrediction ? preCalcAlpha : 0), compID, subTUAbsPartIdx, partIdxesPerSubTU);

										if ((compID != COMPONENT_Cr) && ((transformSkipModeId == 1) ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ()))
										{
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
											COEFF_SCAN_TYPE scanType = COEFF_SCAN_TYPE(pcCU->getCoefScanIdx(uiAbsPartIdx, tuCompRect.width, tuCompRect.height, compID));
#endif
											m_pcEntropyCoder->estimateBit(m_pcTrQuant->m_pcEstBitsSbac, tuCompRect.width, tuCompRect.height, toChannelType(compID)
#if RDOQ_BIT_ESTIMATE_FIX_TICKET29
												, scanType
#endif
											);
										}

#if RDOQ_CHROMA_LAMBDA
										m_pcTrQuant->selectLambda(compID);
#endif

#if JVET_C0024_QTBT
										Pel *pcResiCurrComp = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
										UInt resiStride = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
										Pel *pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
										UInt resiStride = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
#endif

										TCoeff bestCoeffComp[MAX_TU_SIZE*MAX_TU_SIZE];
										Pel    bestResiComp[MAX_TU_SIZE*MAX_TU_SIZE];

#if ADAPTIVE_QP_SELECTION
										TCoeff bestArlCoeffComp[MAX_TU_SIZE*MAX_TU_SIZE];
#endif
										TCoeff     currAbsSum = 0;
										UInt       currCompBits = 0;
										Distortion currCompDist = 0;
										Double     currCompCost = 0;
										UInt       nonCoeffBits = 0;
										Distortion nonCoeffDist = 0;
										Double     nonCoeffCost = 0;

										if (!isOneMode && !isFirstMode)
										{
											memcpy(bestCoeffComp, currentCoefficients, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#if ADAPTIVE_QP_SELECTION
											memcpy(bestArlCoeffComp, currentARLCoefficients, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#endif
											twoDimBufCopier
												(
													(UChar *)bestResiComp,
													tuCompRect.width * sizeof(Pel),
													(UChar *)pcResiCurrComp,
													resiStride * sizeof(Pel),
													tuCompRect.width * sizeof(Pel),
													tuCompRect.height
													);
											
											//for (Int y = 0; y < tuCompRect.height; y++)
											//{
											//	memcpy(&bestResiComp[y * tuCompRect.width], (pcResiCurrComp + (y * resiStride)), (sizeof(Pel) * tuCompRect.width));
											//}
										}

										if (bUseCrossCPrediction)
										{
											TComTrQuant::crossComponentPrediction(TUIterator,
												compID,
												pLumaResi,
												pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
												crossCPredictedResidualBuffer,
												tuCompRect.width,
												tuCompRect.height,
#if JVET_C0024_QTBT
												m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(COMPONENT_Y),
#else
												m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
#endif
												pcResi->getStride(compID),
												tuCompRect.width,
												false);

#if PVC_AQS
											m_pcTrQuant->setEstQSFactor(pcPred->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0), pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, pcCU->getSlice()->getSliceType(), pcCU, uiAbsPartIdx, compID);
											//assert(absTUPartIdx==uiAbsPartIdx );
#endif

											m_pcTrQuant->transformNxN(TUIterator, compID, crossCPredictedResidualBuffer, tuCompRect.width, currentCoefficients,
#if ADAPTIVE_QP_SELECTION
												currentARLCoefficients,
#endif
												currAbsSum, cQP
											);
										}
										else
										{
#if VCEG_AZ08_INTER_KLT
											if (bCheckKLT && uiAbsSum[compID][subTUIndex])
											{
												UInt uiTarDepth = g_aucConvertToBit[tuWidth];
												UInt uiTempSize = g_uiDepth2InterTempSize[uiTarDepth];
												m_pcTrQuant->getTargetPatch(pcCU, uiAbsPartIdx, uiAbsPartIdx, pcPred, tuWidth, uiTempSize);
												m_pcTrQuant->candidateSearch(pcCU, uiAbsPartIdx, tuWidth, uiTempSize);
												useKLT = m_pcTrQuant->candidateTrain(tuWidth, uiTempSize);
											}

#if PVC_AQS
											m_pcTrQuant->setEstQSFactor(pcPred->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0), pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, pcCU->getSlice()->getSliceType(), pcCU, uiAbsPartIdx, compID);
#endif
											m_pcTrQuant->transformNxN(TUIterator, compID, pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0), pcResi->getStride(compID), currentCoefficients,
#if ADAPTIVE_QP_SELECTION
												currentARLCoefficients,
#endif
												currAbsSum, cQP
												, useKLT
											);
#else
#if PVC_AQS
											m_pcTrQuant->setEstQSFactor(pcPred->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0), pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, pcCU->getSlice()->getSliceType(), pcCU, uiAbsPartIdx, compID);
#endif
											m_pcTrQuant->transformNxN(TUIterator, compID, pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0), pcResi->getStride(compID), currentCoefficients,
#if ADAPTIVE_QP_SELECTION
												currentARLCoefficients,
#endif
												currAbsSum, cQP
											);
#endif
										}

										if (isFirstMode || (currAbsSum == 0))
										{
											if (bUseCrossCPrediction)
											{
												TComTrQuant::crossComponentPrediction(TUIterator,
													compID,
													pLumaResi,
													m_pTempPel,
#if JVET_C0024_QTBT
													m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
#else
													m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
#endif
													tuCompRect.width,
													tuCompRect.height,
#if JVET_C0024_QTBT
													m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(COMPONENT_Y),
													tuCompRect.width,
													m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID),
#else
													m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
													tuCompRect.width,
													m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
#endif
													true);

#if JVET_C0024_QTBT
												nonCoeffDist = m_pcRdCost->getDistPart(channelBitDepth, m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID), pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
#else
												nonCoeffDist = m_pcRdCost->getDistPart(channelBitDepth, m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID), pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
#endif
													pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, compID); // initialized with zero residual distortion
											}
											else
											{
												nonCoeffDist = m_pcRdCost->getDistPart(channelBitDepth, m_pTempPel, tuCompRect.width, pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, compID); // initialized with zero residual distortion
											}

											m_pcEntropyCoder->encodeQtCbfZero(TUIterator, toChannelType(compID));

											if (isCrossCPredictionAvailable)
											{
												m_pcEntropyCoder->encodeCrossComponentPrediction(TUIterator, compID);
											}

											nonCoeffBits = m_pcEntropyCoder->getNumberOfWrittenBits();
#if WCG_LUMA_DQP_CM_SCALE
											if (m_pcEncCfg->getUseLumaDeltaQp() > 0)
											{
												nonCoeffCost = m_pcRdCost->calcRdCost(nonCoeffBits, nonCoeffDist, false, DF_DEFAULT_ORI);
											}
											else
#endif
												nonCoeffCost = m_pcRdCost->calcRdCost(nonCoeffBits, nonCoeffDist);
										}

										if ((puiZeroDist != NULL) && isFirstMode)
										{
											*puiZeroDist += nonCoeffDist; // initialized with zero residual distortion
										}

										DEBUG_STRING_NEW(sSingleStringTest)

											if (currAbsSum > 0) //if non-zero coefficients are present, a residual needs to be derived for further prediction
											{
												if (isFirstMode)
												{
#if JVET_C0024_QTBT
													m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
													m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
#endif
													m_pcEntropyCoder->resetBits();
												}

												m_pcEntropyCoder->encodeQtCbf(TUIterator, compID, true);

												if (isCrossCPredictionAvailable)
												{
													m_pcEntropyCoder->encodeCrossComponentPrediction(TUIterator, compID);
												}

												m_pcEntropyCoder->encodeCoeffNxN(TUIterator, currentCoefficients, compID);
												currCompBits = m_pcEntropyCoder->getNumberOfWrittenBits();

#if JVET_C0024_QTBT
												pcResiCurrComp = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
#else
												pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
#endif

#if VCEG_AZ08_INTER_KLT
												TUEntropyCodingParameters codingParameters;
												if (useKLT)
												{
													getTUEntropyCodingParameters(codingParameters, rTu, compID);
													const UInt *scan = codingParameters.scan;
													recoverOrderCoeff(currentCoefficients, scan, tuWidth, tuHeight);
#if ADAPTIVE_QP_SELECTION
													recoverOrderCoeff(currentARLCoefficients, scan, tuWidth, tuHeight);
#endif
												}
#endif
#if JVET_C0024_QTBT
												m_pcTrQuant->invTransformNxN(TUIterator, compID, pcResiCurrComp, m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID), currentCoefficients, cQP
#if VCEG_AZ08_INTER_KLT
													, useKLT
#endif
													DEBUG_STRING_PASS_INTO_OPTIONAL(&sSingleStringTest, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask)));
#else
												m_pcTrQuant->invTransformNxN(TUIterator, compID, pcResiCurrComp, m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID), currentCoefficients, cQP
#if VCEG_AZ08_INTER_KLT
													, useKLT
#endif
													DEBUG_STRING_PASS_INTO_OPTIONAL(&sSingleStringTest, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask))
												);
#endif
#if JVET_F0096_BILATERAL_FILTER
												if (pcCU->getSlice()->getSPS()->getUseBilateralFilter())
												{
													if (isLuma(compID))
													{
														UInt minSize = std::min(tuWidth, tuHeight);
														if ((currAbsSum > 0) && (pcCU->getQP(uiAbsPartIdx) > 17) && (minSize < 16))
														{
															Pel *pcPtrPred = pcPred->getAddr(compID, uiAbsPartIdx);
															UInt uiStridePred = pcPred->getStride(compID);
#if JVET_C0024_QTBT
															Pel *pcPtrRes = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
															UInt uiStrideRes = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
															Pel *pcPtrRec = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
															UInt uiStrideRec = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
															Pel *pcPtrRes = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddr(compID, uiAbsPartIdx);
															UInt uiStrideRes = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
															Pel *pcPtrRec = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddr(compID, uiAbsPartIdx);
															UInt uiStrideRec = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
#endif
															const Int clipbd = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
															TComBilateralFilter::instance()->bilateralFilterInter(pcCU, tuWidth, tuHeight, pcPtrRes, uiStrideRes, pcPtrPred, uiStridePred, pcPtrRec, uiStrideRec, clipbd, pcCU->getQP(uiAbsPartIdx));
														}
													}
												}
#endif

#if VCEG_AZ08_INTER_KLT
												if (useKLT)
												{
													const UInt *scan = codingParameters.scan;
													reOrderCoeff(currentCoefficients, scan, tuWidth, tuHeight);
#if ADAPTIVE_QP_SELECTION
													reOrderCoeff(currentARLCoefficients, scan, tuWidth, tuHeight);
#endif
												}
#endif
												if (bUseCrossCPrediction)
												{
													TComTrQuant::crossComponentPrediction(TUIterator,
														compID,
														pLumaResi,
#if JVET_C0024_QTBT
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														tuCompRect.width,
														tuCompRect.height,
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(COMPONENT_Y),
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID),
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID),
#else
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														tuCompRect.width,
														tuCompRect.height,
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
#endif
														true);
												}

#if JVET_C0024_QTBT
												currCompDist = m_pcRdCost->getDistPart(channelBitDepth, m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID),
#else
												currCompDist = m_pcRdCost->getDistPart(channelBitDepth, m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
#endif
													pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
													pcResi->getStride(compID),
													tuCompRect.width, tuCompRect.height, compID);

#if WCG_LUMA_DQP_CM_SCALE
												currCompCost = m_pcRdCost->calcRdCost(currCompBits, currCompDist, false, DF_DEFAULT_ORI);
#else
												currCompCost = m_pcRdCost->calcRdCost(currCompBits, currCompDist);
#endif
												if (pcCU->isLosslessCoded(0))
												{
													nonCoeffCost = MAX_DOUBLE;
												}
#if VCEG_AZ08_INTER_KLT
												if (bCheckKLT && useKLT == false)
												{
													nonCoeffCost = MAX_DOUBLE;
													currCompCost = MAX_DOUBLE;
												}
#endif
											}
											else if ((transformSkipModeId == 1) && !bUseCrossCPrediction)
											{
												currCompCost = MAX_DOUBLE;
											}
#if VCEG_AZ08_INTER_KLT
											else if ((KLTModeId == 1) && !bUseCrossCPrediction)
											{
												currCompCost = MAX_DOUBLE;
											}
#endif
											else
											{
												currCompBits = nonCoeffBits;
												currCompDist = nonCoeffDist;
												currCompCost = nonCoeffCost;
											}

										// evaluate
										if ((currCompCost < minCost[compID][subTUIndex]) || ((transformSkipModeId == 1) && (currCompCost == minCost[compID][subTUIndex])))
										{
											bestExplicitRdpcmModeUnSplit[compID][subTUIndex] = pcCU->getExplicitRdpcmMode(compID, subTUAbsPartIdx);

											if (isFirstMode) //check for forced null
											{
												if ((nonCoeffCost < currCompCost) || (currAbsSum == 0))
												{
													memset(currentCoefficients, 0, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));

													currAbsSum = 0;
													currCompBits = nonCoeffBits;
													currCompDist = nonCoeffDist;
													currCompCost = nonCoeffCost;
												}
											}

#if DEBUG_STRING
											if (currAbsSum > 0)
											{
												DEBUG_STRING_SWAP(sSingleStringComp[compID], sSingleStringTest)
											}
											else
											{
												sSingleStringComp[compID].clear();
											}
#endif

											uiAbsSum[compID][subTUIndex] = currAbsSum;
											uiSingleDistComp[compID][subTUIndex] = currCompDist;
											minCost[compID][subTUIndex] = currCompCost;
											uiBestTransformMode[compID][subTUIndex] = transformSkipModeId;
											bestCrossCPredictionAlpha[compID][subTUIndex] = (crossCPredictionModeId == 1) ? pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) : 0;
#if VCEG_AZ08_INTER_KLT
											bestKLTMode[compID][subTUIndex] = KLTModeId;
#endif
#if COM16_C806_EMT
											bestEmtTrIdx[compID][subTUIndex] = emtTrIdx;
#endif
											if (uiAbsSum[compID][subTUIndex] == 0)
											{
												if (bUseCrossCPrediction)
												{
													TComTrQuant::crossComponentPrediction(TUIterator,
														compID,
														pLumaResi,
														m_pTempPel,
#if JVET_C0024_QTBT
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														tuCompRect.width,
														tuCompRect.height,
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(COMPONENT_Y),
														tuCompRect.width,
														m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID),
#else
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
														tuCompRect.width,
														tuCompRect.height,
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
														tuCompRect.width,
														m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
#endif
														true);
												}
												else
												{
#if JVET_C0024_QTBT
													pcResiCurrComp = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
													const UInt uiStride = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
#else
													pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
													const UInt uiStride = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
#endif
													for (UInt uiY = 0; uiY < tuCompRect.height; uiY++)
													{
														memset(pcResiCurrComp, 0, (sizeof(Pel) * tuCompRect.width));
														pcResiCurrComp += uiStride;
													}
												}
											}
										}
										else
										{
											// reset
											memcpy(currentCoefficients, bestCoeffComp, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#if ADAPTIVE_QP_SELECTION
											memcpy(currentARLCoefficients, bestArlCoeffComp, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#endif
											twoDimBufCopier
												(
													(UChar *)pcResiCurrComp,
													resiStride * sizeof(Pel),
													(UChar *)bestResiComp,
													tuCompRect.width * sizeof(Pel),
													tuCompRect.width * sizeof(Pel),
													tuCompRect.height
													);
											
											//for (Int y = 0; y < tuCompRect.height; y++)
											//{
											//	memcpy((pcResiCurrComp + (y * resiStride)), &bestResiComp[y * tuCompRect.width], (sizeof(Pel) * tuCompRect.width));
											//}
										}
									}
#if COM16_C806_EMT
								}
#endif
							}
#if VCEG_AZ08_INTER_KLT
						}
#endif
						pcCU->setExplicitRdpcmModePartRange(bestExplicitRdpcmModeUnSplit[compID][subTUIndex], compID, subTUAbsPartIdx, partIdxesPerSubTU);
						pcCU->setTransformSkipPartRange(uiBestTransformMode[compID][subTUIndex], compID, subTUAbsPartIdx, partIdxesPerSubTU);
						pcCU->setCbfPartRange((((uiAbsSum[compID][subTUIndex] > 0) ? 1 : 0) << uiTrMode), compID, subTUAbsPartIdx, partIdxesPerSubTU);
						pcCU->setCrossComponentPredictionAlphaPartRange(bestCrossCPredictionAlpha[compID][subTUIndex], compID, subTUAbsPartIdx, partIdxesPerSubTU);
#if COM16_C806_EMT
						pcCU->setEmtTuIdxPartsRange(bestEmtTrIdx[compID][subTUIndex], compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
#if VCEG_AZ08_INTER_KLT
						pcCU->setKLTPartRange(bestKLTMode[compID][subTUIndex], compID, subTUAbsPartIdx, partIdxesPerSubTU);
#endif
					} while (TUIterator.nextSection(rTu)); //end of sub-TU loop
				} // processing section
			} // component loop

#if JVET_C0024_QTBT
			m_pcRDGoOnSbacCoder->load(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_QT_TRAFO_ROOT]);
#else
			for (UInt ch = 0; ch < numValidComp; ch++)
			{
				const ComponentID compID = ComponentID(ch);
				if (rTu.ProcessComponentSection(compID) && (rTu.getRect(compID).width != rTu.getRect(compID).height))
				{
					offsetSubTUCBFs(rTu, compID); //the CBFs up to now have been defined for two sub-TUs - shift them down a level and replace with the parent level CBF
				}
			}

			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
#endif
			m_pcEntropyCoder->resetBits();

#if !JVET_C0024_QTBT
			if (uiLog2TrSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx))
			{
#if COM16_C806_T64
				m_pcEntropyCoder->encodeTransformSubdivFlag(0, pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize);
#else
				m_pcEntropyCoder->encodeTransformSubdivFlag(0, 5 - uiLog2TrSize);
#endif
			}
#endif

			for (UInt ch = 0; ch < numValidComp; ch++)
			{
				const UInt chOrderChange = ((ch + 1) == numValidComp) ? 0 : (ch + 1);
				const ComponentID compID = ComponentID(chOrderChange);
				if (rTu.ProcessComponentSection(compID))
				{
					m_pcEntropyCoder->encodeQtCbf(rTu, compID, true);
				}
			}

			for (UInt ch = 0; ch < numValidComp; ch++)
			{
				const ComponentID compID = ComponentID(ch);
				if (rTu.ProcessComponentSection(compID))
				{
					if (isChroma(compID) && (uiAbsSum[COMPONENT_Y][0] != 0))
					{
						m_pcEntropyCoder->encodeCrossComponentPrediction(rTu, compID);
					}

					m_pcEntropyCoder->encodeCoeffNxN(rTu, pcCoeffCurr[compID], compID);
					for (UInt subTUIndex = 0; subTUIndex < 2; subTUIndex++)
					{
						uiSingleDist += uiSingleDistComp[compID][subTUIndex];
					}
#if JVET_C0024_QTBT
					assert(uiSingleDistComp[compID][1] == 0);
#endif
				}
			}

			uiSingleBits = m_pcEntropyCoder->getNumberOfWrittenBits();

#if WCG_LUMA_DQP_CM_SCALE
			if (m_pcEncCfg->getUseLumaDeltaQp() > 0)
			{
				dSingleCost = m_pcRdCost->calcRdCost(uiSingleBits, uiSingleDist, false, DF_DEFAULT_ORI);
			}
			else
#endif
				dSingleCost = m_pcRdCost->calcRdCost(uiSingleBits, uiSingleDist);

#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
			if (pcCU->getSlice()->getSPS()->getUseInterKLT())
			{
#endif
				// The reconstruction will be updated for each TU in inter mode coding
				// Save reconstruction to (1) facilitate the following utilization in next blocks; (2) backup the possible best reconstruction.
				Bool bCodeChroma = false;
				Int numValidCompCurr = bCodeChroma ? numValidComp : 1;
				for (UInt ch = 0; ch < numValidCompCurr; ch++)
				{
					const ComponentID compID = ComponentID(ch);
#if !JVET_D0033_ADAPTIVE_CLIPPING
					const Int clipbd = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
#endif
					const TComRectangle &rect = rTu.getRect(compID);
					const UInt           uiWidth = rect.width;
					const UInt           uiHeight = rect.height;
					Pel *pcPtrPred = pcPred->getAddr(compID, uiAbsPartIdx);
					UInt uiStridePred = pcPred->getStride(compID);
#if JVET_C0024_QTBT
					Pel *pcPtrRes = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
					UInt uiStrideRes = m_ppcQTTempTComYuv[uiWIdx][uiHIdx].getStride(compID);
					Pel *pcPtrRec = m_ppcQTTempTComYuvRec[uiWIdx][uiHIdx].getAddr(compID, uiAbsPartIdx);
					UInt uiStrideRec = m_ppcQTTempTComYuvRec[uiWIdx][uiHIdx].getStride(compID);
#else
					Pel *pcPtrRes = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddr(compID, uiAbsPartIdx);
					UInt uiStrideRes = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
					Pel *pcPtrRec = m_pcQTTempTComYuvRec[uiQTTempAccessLayer].getAddr(compID, uiAbsPartIdx);
					UInt uiStrideRec = m_pcQTTempTComYuvRec[uiQTTempAccessLayer].getStride(compID);
#endif
					UInt    uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
					Pel*    piRecIPred = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
					UInt    uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);

					for (UInt uiY = 0; uiY < uiHeight; ++uiY)
					{
						for (UInt uiX = 0; uiX < uiWidth; ++uiX)
						{
							//save the temporary reconstruction to temporary buffer.
#if JVET_D0033_ADAPTIVE_CLIPPING
#if O0043_BEST_EFFORT_DECODING
							pcPtrRec[uiX] = ClipA(rightShiftEvenRounding<Pel>(pcPtrPred[uiX] + pcPtrRes[uiX], bitDepthDelta), compID);
#else
							pcPtrRec[uiX] = ClipA(pcPtrPred[uiX] + pcPtrRes[uiX], compID);
#endif
#else
#if O0043_BEST_EFFORT_DECODING
							pcPtrRec[uiX] = ClipBD(rightShiftEvenRounding<Pel>(pcPtrPred[uiX] + pcPtrRes[uiX], bitDepthDelta), clipbd);
#else
							pcPtrRec[uiX] = ClipBD(pcPtrPred[uiX] + pcPtrRes[uiX], clipbd);
#endif
#endif
							//save the temporary reconstruction to reconstruction picture.
							piRecIPred[uiX] = pcPtrRec[uiX];
						}
						pcPtrPred += uiStridePred;
						pcPtrRes += uiStrideRes;
						pcPtrRec += uiStrideRec;
						piRecIPred += uiRecIPredStride;
					}
				}
#if VCEG_AZ08_USE_KLT
			}
#endif
#endif
		} // check full


#if !JVET_C0024_QTBT
		  // code sub-blocks
		if (bCheckSplit)
		{
			if (bCheckFull)
			{
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_TEST]);
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
			}
			Distortion uiSubdivDist = 0;
			UInt       uiSubdivBits = 0;
			Double     dSubdivCost = 0.0;

			//save the non-split CBFs in case we need to restore them later

			UInt bestCBF[MAX_NUM_COMPONENT];
			UInt bestsubTUCBF[MAX_NUM_COMPONENT][2];
			for (UInt ch = 0; ch < numValidComp; ch++)
			{
				const ComponentID compID = ComponentID(ch);

				if (rTu.ProcessComponentSection(compID))
				{
					bestCBF[compID] = pcCU->getCbf(uiAbsPartIdx, compID, uiTrMode);

					const TComRectangle &tuCompRect = rTu.getRect(compID);
					if (tuCompRect.width != tuCompRect.height)
					{
						const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

						for (UInt subTU = 0; subTU < 2; subTU++)
						{
							bestsubTUCBF[compID][subTU] = pcCU->getCbf((uiAbsPartIdx + (subTU * partIdxesPerSubTU)), compID, subTUDepth);
						}
					}
				}
			}


			TComTURecurse tuRecurseChild(rTu, false);
			const UInt uiQPartNumSubdiv = tuRecurseChild.GetAbsPartIdxNumParts();

			DEBUG_STRING_NEW(sSplitString[MAX_NUM_COMPONENT])

				do
				{
					DEBUG_STRING_NEW(childString)
#if VCEG_AZ08_INTER_KLT || JVET_F0096_BILATERAL_FILTER || PVC_AQS
						xEstimateInterResidualQT(pcResi, dSubdivCost, uiSubdivBits, uiSubdivDist, bCheckFull ? NULL : puiZeroDist, tuRecurseChild, pcPred DEBUG_STRING_PASS_INTO(childString));
#else
						xEstimateInterResidualQT(pcResi, dSubdivCost, uiSubdivBits, uiSubdivDist, bCheckFull ? NULL : puiZeroDist, tuRecurseChild DEBUG_STRING_PASS_INTO(childString));
#endif
#if DEBUG_STRING
					// split the string by component and append to the relevant output (because decoder decodes in channel order, whereas this search searches by TU-order)
					std::size_t lastPos = 0;
					const std::size_t endStrng = childString.find(debug_reorder_data_inter_token[MAX_NUM_COMPONENT], lastPos);
					for (UInt ch = 0; ch < numValidComp; ch++)
					{
						if (lastPos != std::string::npos && childString.find(debug_reorder_data_inter_token[ch], lastPos) == lastPos)
						{
							lastPos += strlen(debug_reorder_data_inter_token[ch]); // skip leading string
						}
						std::size_t pos = childString.find(debug_reorder_data_inter_token[ch + 1], lastPos);
						if (pos != std::string::npos && pos>endStrng)
						{
							lastPos = endStrng;
						}
						sSplitString[ch] += childString.substr(lastPos, (pos == std::string::npos) ? std::string::npos : (pos - lastPos));
						lastPos = pos;
					}
#endif
#if COM16_C806_EMT
					if (dSubdivCost >= dSingleCost && m_pcEncCfg->getUseFastInterEMT())
					{
						break;
					}
#endif
				} while (tuRecurseChild.nextSection(rTu));

				UInt uiCbfAny = 0;
#if COM16_C806_EMT
				Bool bSplitSelected = false;
				if ((dSubdivCost < dSingleCost && m_pcEncCfg->getUseFastInterEMT()) || !m_pcEncCfg->getUseFastInterEMT())
				{
#endif
					for (UInt ch = 0; ch < numValidComp; ch++)
					{
						UInt uiYUVCbf = 0;
						for (UInt ui = 0; ui < 4; ++ui)
						{
							uiYUVCbf |= pcCU->getCbf(uiAbsPartIdx + ui * uiQPartNumSubdiv, ComponentID(ch), uiTrMode + 1);
						}
						UChar *pBase = pcCU->getCbf(ComponentID(ch));
						const UInt flags = uiYUVCbf << uiTrMode;
						for (UInt ui = 0; ui < 4 * uiQPartNumSubdiv; ++ui)
						{
							pBase[uiAbsPartIdx + ui] |= flags;
						}
						uiCbfAny |= uiYUVCbf;
					}

					m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_ROOT]);
					m_pcEntropyCoder->resetBits();

					// when compID isn't a channel, code Cbfs:
					xEncodeInterResidualQT(MAX_NUM_COMPONENT, rTu);
					for (UInt ch = 0; ch < numValidComp; ch++)
					{
						xEncodeInterResidualQT(ComponentID(ch), rTu);
					}

					uiSubdivBits = m_pcEntropyCoder->getNumberOfWrittenBits();
					dSubdivCost = m_pcRdCost->calcRdCost(uiSubdivBits, uiSubdivDist);
					if (!bCheckFull || (uiCbfAny && (dSubdivCost < dSingleCost)))
					{
						rdCost += dSubdivCost;
						ruiBits += uiSubdivBits;
						ruiDist += uiSubdivDist;
#if COM16_C806_EMT
						if (m_pcEncCfg->getUseFastInterEMT())
						{
							bSplitSelected = true;
						}
#endif
#if DEBUG_STRING
						for (UInt ch = 0; ch < numValidComp; ch++)
						{
							DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[ch])
								DEBUG_STRING_APPEND(sDebug, sSplitString[ch])
						}
#endif
					}
#if COM16_C806_EMT
				}
				if ((!bSplitSelected && m_pcEncCfg->getUseFastInterEMT()) || (!m_pcEncCfg->getUseFastInterEMT() && !(!bCheckFull || (uiCbfAny && (dSubdivCost < dSingleCost)))))
#else
					else
#endif
					{
						rdCost += dSingleCost;
						ruiBits += uiSingleBits;
						ruiDist += uiSingleDist;

						//restore state to unsplit

						pcCU->setTrIdxSubParts(uiTrMode, uiAbsPartIdx, uiDepth);

						for (UInt ch = 0; ch < numValidComp; ch++)
						{
							const ComponentID compID = ComponentID(ch);

							DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[ch])
								if (rTu.ProcessComponentSection(compID))
								{
									DEBUG_STRING_APPEND(sDebug, sSingleStringComp[compID])

										const Bool splitIntoSubTUs = rTu.getRect(compID).width != rTu.getRect(compID).height;
									const UInt numberOfSections = splitIntoSubTUs ? 2 : 1;
									const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> (splitIntoSubTUs ? 1 : 0);

									for (UInt subTUIndex = 0; subTUIndex < numberOfSections; subTUIndex++)
									{
										const UInt  uisubTUPartIdx = uiAbsPartIdx + (subTUIndex * partIdxesPerSubTU);

										if (splitIntoSubTUs)
										{
											const UChar combinedCBF = (bestsubTUCBF[compID][subTUIndex] << subTUDepth) | (bestCBF[compID] << uiTrMode);
											pcCU->setCbfPartRange(combinedCBF, compID, uisubTUPartIdx, partIdxesPerSubTU);
										}
										else
										{
											pcCU->setCbfPartRange((bestCBF[compID] << uiTrMode), compID, uisubTUPartIdx, partIdxesPerSubTU);
										}

										pcCU->setCrossComponentPredictionAlphaPartRange(bestCrossCPredictionAlpha[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
										pcCU->setTransformSkipPartRange(uiBestTransformMode[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
										pcCU->setExplicitRdpcmModePartRange(bestExplicitRdpcmModeUnSplit[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
#if COM16_C806_EMT
										pcCU->setEmtTuIdxPartsRange(bestEmtTrIdx[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
#endif
#if VCEG_AZ08_INTER_KLT
										pcCU->setKLTPartRange(bestKLTMode[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
#endif
									}
								}
						}

						m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_QT_TRAFO_TEST]);
#if VCEG_AZ08_INTER_KLT
#if VCEG_AZ08_USE_KLT
						if (pcCU->getSlice()->getSPS()->getUseInterKLT())
						{
#endif
							//Set reconstruction for the possible utility of template of next blocks
							UInt trWidth = 1 << uiLog2TrSize;
							UInt trHeight = 1 << uiLog2TrSize;
							const UInt uiQTLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
							Bool bCodeChroma = false;
							Int numValidCompCurr = bCodeChroma ? numValidComp : 1;
							for (UInt ch = 0; ch < numValidCompCurr; ch++)
							{
								const ComponentID compID = ComponentID(ch);
								UInt    uiZOrder = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
								Pel*  piSrc = m_pcQTTempTComYuvRec[uiQTLayer].getAddr(compID, uiAbsPartIdx);
								UInt  uiSrcStride = m_pcQTTempTComYuvRec[uiQTLayer].getStride(compID);
								Pel*  piDes = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), uiZOrder);
								UInt  uiDesStride = pcCU->getPic()->getPicYuvRec()->getStride(compID);
								for (UInt uiY = 0; uiY < trHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride)
								{
									for (UInt uiX = 0; uiX < trWidth; uiX++)
									{
										piDes[uiX] = piSrc[uiX];
									}
								}
							}
#if VCEG_AZ08_USE_KLT
						}
#endif
#endif
					}
		}
				else
#endif
				{
					rdCost += dSingleCost;
					ruiBits += uiSingleBits;
					ruiDist += uiSingleDist;
#if DEBUG_STRING
					for (UInt ch = 0; ch < numValidComp; ch++)
					{
						const ComponentID compID = ComponentID(ch);
						DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[compID])

							if (rTu.ProcessComponentSection(compID))
							{
								DEBUG_STRING_APPEND(sDebug, sSingleStringComp[compID])
							}
					}
#endif
				}
				DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[MAX_NUM_COMPONENT])
	}

#if !JVET_C0024_QTBT
		Void TEncSearch::xEncodeInterResidualQT(const ComponentID compID, TComTU &rTu)
		{
			TComDataCU* pcCU = rTu.getCU();
			const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
			const UInt uiCurrTrMode = rTu.GetTransformDepthRel();
			assert(pcCU->getDepth(0) == pcCU->getDepth(uiAbsPartIdx));
			const UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);

			const Bool bSubdiv = uiCurrTrMode != uiTrMode;

			const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();

			if (compID == MAX_NUM_COMPONENT)  // we are not processing a channel, instead we always recurse and code the CBFs
			{
				if (uiLog2TrSize <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() && uiLog2TrSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx))
				{
					if ((pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) && (pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N))
					{
						assert(bSubdiv); // Inferred splitting rule - see derivation and use of interSplitFlag in the specification.
					}
					else
					{
#if COM16_C806_T64
						m_pcEntropyCoder->encodeTransformSubdivFlag(bSubdiv, pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize);
#else
						m_pcEntropyCoder->encodeTransformSubdivFlag(bSubdiv, 5 - uiLog2TrSize);
#endif
					}
				}

				assert(!pcCU->isIntra(uiAbsPartIdx));

				const Bool bFirstCbfOfCU = uiCurrTrMode == 0;

				for (UInt ch = COMPONENT_Cb; ch<pcCU->getPic()->getNumberValidComponents(); ch++)
				{
					const ComponentID compIdInner = ComponentID(ch);
					if (bFirstCbfOfCU || rTu.ProcessingAllQuadrants(compIdInner))
					{
						if (bFirstCbfOfCU || pcCU->getCbf(uiAbsPartIdx, compIdInner, uiCurrTrMode - 1))
						{
							m_pcEntropyCoder->encodeQtCbf(rTu, compIdInner, !bSubdiv);
						}
					}
					else
					{
						assert(pcCU->getCbf(uiAbsPartIdx, compIdInner, uiCurrTrMode) == pcCU->getCbf(uiAbsPartIdx, compIdInner, uiCurrTrMode - 1));
					}
				}

				if (!bSubdiv)
				{
					m_pcEntropyCoder->encodeQtCbf(rTu, COMPONENT_Y, true);
				}
			}

			if (!bSubdiv)
			{
				if (compID != MAX_NUM_COMPONENT) // we have already coded the CBFs, so now we code coefficients
				{
					if (rTu.ProcessComponentSection(compID))
					{
						if (isChroma(compID) && (pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrMode) != 0))
						{
							m_pcEntropyCoder->encodeCrossComponentPrediction(rTu, compID);
						}

						if (pcCU->getCbf(uiAbsPartIdx, compID, uiTrMode) != 0)
						{
							const UInt uiQTTempAccessLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
							TCoeff *pcCoeffCurr = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + rTu.getCoefficientOffset(compID);
							m_pcEntropyCoder->encodeCoeffNxN(rTu, pcCoeffCurr, compID);
						}
					}
				}
			}
			else
			{
				if (compID == MAX_NUM_COMPONENT || pcCU->getCbf(uiAbsPartIdx, compID, uiCurrTrMode))
				{
					TComTURecurse tuRecurseChild(rTu, false);
					do
					{
						xEncodeInterResidualQT(compID, tuRecurseChild);
					} while (tuRecurseChild.nextSection(rTu));
				}
			}
		}
#endif
		Void TEncSearch::xSetInterResidualQTData(TComYuv* pcResi, Bool bSpatial, TComTU &rTu) // TODO: turn this into two functions for bSpatial=true and false.
		{
			TComDataCU* pcCU = rTu.getCU();
			const UInt uiCurrTrMode = rTu.GetTransformDepthRel();
			const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
			assert(pcCU->getDepth(0) == pcCU->getDepth(uiAbsPartIdx));
#if !JVET_C0024_QTBT
			const UInt uiTrMode = pcCU->getTransformIdx(uiAbsPartIdx);
#endif
			const TComSPS *sps = pcCU->getSlice()->getSPS();

#if JVET_C0024_QTBT
			assert(uiCurrTrMode == 0 && uiAbsPartIdx == 0);
			UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
			UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
#else
			if (uiCurrTrMode == uiTrMode)
			{
				const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
				const UInt uiQTTempAccessLayer = sps->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
#endif

				if (bSpatial)
				{
					// Data to be copied is in the spatial domain, i.e., inverse-transformed.

					for (UInt i = 0; i<pcResi->getNumberValidComponents(); i++)
					{
						const ComponentID compID = ComponentID(i);
						if (rTu.ProcessComponentSection(compID))
						{
							const TComRectangle &rectCompTU(rTu.getRect(compID));
#if JVET_C0024_QTBT
							m_ppcQTTempTComYuv[uiWIdx][uiHIdx].copyPartToPartComponentMxN(compID, pcResi, rectCompTU);
#else
							m_pcQTTempTComYuv[uiQTTempAccessLayer].copyPartToPartComponentMxN(compID, pcResi, rectCompTU);
#endif
						}
					}
				}
				else
				{
					for (UInt ch = 0; ch < getNumberValidComponents(sps->getChromaFormatIdc()); ch++)
					{
						const ComponentID compID = ComponentID(ch);
						if (rTu.ProcessComponentSection(compID))
						{
							const TComRectangle &rectCompTU(rTu.getRect(compID));
							const UInt numCoeffInBlock = rectCompTU.width * rectCompTU.height;
							const UInt offset = rTu.getCoefficientOffset(compID);
							TCoeff* dest = pcCU->getCoeff(compID) + offset;
#if JVET_C0024_QTBT
							assert(offset == 0);
							const TCoeff* src = m_pppcQTTempCoeff[compID][uiWIdx][uiHIdx] + offset;
#else
							const TCoeff* src = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + offset;
#endif
							::memcpy(dest, src, sizeof(TCoeff)*numCoeffInBlock);

#if ADAPTIVE_QP_SELECTION
#if JVET_C0024_QTBT
							TCoeff* pcArlCoeffSrc = m_pppcQTTempArlCoeff[compID][uiWIdx][uiHIdx] + offset;
#else
							TCoeff* pcArlCoeffSrc = m_ppcQTTempArlCoeff[compID][uiQTTempAccessLayer] + offset;
#endif
							TCoeff* pcArlCoeffDst = pcCU->getArlCoeff(compID) + offset;
							::memcpy(pcArlCoeffDst, pcArlCoeffSrc, sizeof(TCoeff) * numCoeffInBlock);
#endif
						}
					}
				}
#if !JVET_C0024_QTBT
			}
			else
			{

				TComTURecurse tuRecurseChild(rTu, false);
				do
				{
					xSetInterResidualQTData(pcResi, bSpatial, tuRecurseChild);
				} while (tuRecurseChild.nextSection(rTu));
			}
#endif
		}




		UInt TEncSearch::xModeBitsIntra(TComDataCU* pcCU, UInt uiMode, UInt uiPartOffset, UInt uiDepth, const ChannelType chType
#if VCEG_AZ07_INTRA_65ANG_MODES
			, Int* piModes, Int  iCase
#endif
		)
		{
			// Reload only contexts required for coding intra mode information
#if JVET_C0024_QTBT
			UInt uiWIdx = g_aucConvertToBit[pcCU->getWidth(0)];
			UInt uiHIdx = g_aucConvertToBit[pcCU->getHeight(0)];
			m_pcRDGoOnSbacCoder->loadIntraDirMode(m_ppppcRDSbacCoder[uiWIdx][uiHIdx][CI_CURR_BEST], chType);
#else
			m_pcRDGoOnSbacCoder->loadIntraDirMode(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST], chType);
#endif

			// Temporarily set the intra dir being tested, and only
			// for absPartIdx, since encodeIntraDirModeLuma/Chroma only use
			// the entry at absPartIdx.

			UChar &rIntraDirVal = pcCU->getIntraDir(chType)[uiPartOffset];
			UChar origVal = rIntraDirVal;
			rIntraDirVal = uiMode;
			//pcCU->setIntraDirSubParts ( chType, uiMode, uiPartOffset, uiDepth + uiInitTrDepth );

			m_pcEntropyCoder->resetBits();
			if (isLuma(chType))
			{
				m_pcEntropyCoder->encodeIntraDirModeLuma(pcCU, uiPartOffset
#if VCEG_AZ07_INTRA_65ANG_MODES
					, false, piModes, iCase
#endif
				);
			}
			else
			{
				m_pcEntropyCoder->encodeIntraDirModeChroma(pcCU, uiPartOffset);
			}

			rIntraDirVal = origVal; // restore

			return m_pcEntropyCoder->getNumberOfWrittenBits();
		}




#if JVET_C0024_FAST_MRG
		UInt TEncSearch::updateCandList(UInt uiMode, Double uiCost, UInt uiFastCandNum, UInt * CandModeList, Double * CandCostList
#if MLIP_EXTEND_REF_LINE
			, UInt * ExternRefList, UInt ExternRef
#endif
		)
#else
		UInt TEncSearch::xUpdateCandList(UInt uiMode, Double uiCost, UInt uiFastCandNum, UInt * CandModeList, Double * CandCostList)
#endif
		{
			UInt i;
			UInt shift = 0;

			while (shift<uiFastCandNum && uiCost<CandCostList[uiFastCandNum - 1 - shift])
			{
				shift++;
			}

			if (shift != 0)
			{
				for (i = 1; i<shift; i++)
				{
					CandModeList[uiFastCandNum - i] = CandModeList[uiFastCandNum - 1 - i];
					CandCostList[uiFastCandNum - i] = CandCostList[uiFastCandNum - 1 - i];
#if MLIP_EXTEND_REF_LINE
					if (ExternRefList)
						ExternRefList[uiFastCandNum - i] = ExternRefList[uiFastCandNum - 1 - i];
#endif
				}
				CandModeList[uiFastCandNum - shift] = uiMode;
				CandCostList[uiFastCandNum - shift] = uiCost;
#if MLIP_EXTEND_REF_LINE
				if (ExternRefList)
					ExternRefList[uiFastCandNum - shift] = ExternRef;
#endif
				return 1;
			}

			return 0;
		}





		/** add inter-prediction syntax elements for a CU block
		* \param pcCU
		* \param uiQp
		* \param uiTrMode
		* \param ruiBits
		* \returns Void
		*/
		Void  TEncSearch::xAddSymbolBitsInter(TComDataCU* pcCU, UInt& ruiBits)
		{
#if JVET_C0024_QTBT
			if (pcCU->getMergeFlag(0) && !pcCU->getQtRootCbf(0))
#else
			if (pcCU->getMergeFlag(0) && pcCU->getPartitionSize(0) == SIZE_2Nx2N && !pcCU->getQtRootCbf(0))
#endif
			{
				pcCU->setSkipFlagSubParts(true, 0, pcCU->getDepth(0));

				m_pcEntropyCoder->resetBits();
				if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
				{
					m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
				}
				m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL
				if (pcCU->isCodeBIF(pcCU->getWidth(0), pcCU->getDepth(0)
#if JVET_C0024_QTBT && SIF_QTBT_FIX
					+ pcCU->getBTDepth(0) // check total depth
#endif      
				))
				{
					m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
				}
#endif
#if VCEG_AZ07_FRUC_MERGE
				m_pcEntropyCoder->encodeFRUCMgrMode(pcCU, 0, 0);
				if (!pcCU->getFRUCMgrMode(0))
#endif
#if MRG_CONSIDER_COMPLEX
				{
					if ((pcCU->getWidth(0) >= 16) || (pcCU->getHeight(0) >= 16))
					{
						m_pcEntropyCoder->encodeMergeComplexFlag(pcCU, 0);
						if (pcCU->getMergeComplexFlag(0))
						{
							m_pcEntropyCoder->encodeMergeComplexIndex(pcCU, 0, true);
						}
					}

					if (!pcCU->getMergeComplexFlag(0))
#endif
#if COM16_C1016_AFFINE
					{
						if (pcCU->isAffineMrgFlagCoded(0, 0))
						{
							m_pcEntropyCoder->encodeAffineFlag(pcCU, 0, 0);
						}

						if (!pcCU->isAffine(0))
						{
							m_pcEntropyCoder->encodeMergeIndex(pcCU, 0, true);
						}
					}
#else
						m_pcEntropyCoder->encodeMergeIndex(pcCU, 0, true);
#endif
#if MRG_CONSIDER_COMPLEX
				}
#endif
#if SIF_CU_LEVEL
				m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
#endif

#if VCEG_AZ06_IC
				m_pcEntropyCoder->encodeICFlag(pcCU, 0);
#endif
				ruiBits += m_pcEntropyCoder->getNumberOfWrittenBits();
			}
			else
			{
				m_pcEntropyCoder->resetBits();

				if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
				{
					m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
				}

				m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
#if RRC_SELECTABLE_IF && !SIF_SPLIT_ONLY && !SIF_CU_LEVEL
				if (pcCU->isCodeBIF(pcCU->getWidth(0), pcCU->getDepth(0)
#if JVET_C0024_QTBT && SIF_QTBT_FIX
					+ pcCU->getBTDepth(0) // check total depth
#endif
				))
				{
					m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
				}
#endif
				m_pcEntropyCoder->encodePredMode(pcCU, 0, true);
#if !JVET_C0024_QTBT
				m_pcEntropyCoder->encodePartSize(pcCU, 0, pcCU->getDepth(0), true);
#endif
				m_pcEntropyCoder->encodePredInfo(pcCU, 0);
#if COM16_C806_OBMC
				m_pcEntropyCoder->encodeOBMCFlag(pcCU, 0, true);
#endif
#if SIF_CU_LEVEL
				m_pcEntropyCoder->encodeAdaptiveIfIndex(pcCU, 0);
#endif
#if VCEG_AZ06_IC
				m_pcEntropyCoder->encodeICFlag(pcCU, 0);
#endif
				Bool codeDeltaQp = false;
				Bool codeChromaQpAdj = false;
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
				Int bNonZeroCoeff = false;
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
				Int iNonZeroCoeffNonTs;
#endif
				m_pcEntropyCoder->encodeCoeff(pcCU, 0, pcCU->getDepth(0), codeDeltaQp, codeChromaQpAdj
#if VCEG_AZ05_ROT_TR  || VCEG_AZ05_INTRA_MPI || COM16_C1044_NSST || COM16_C1046_PDPC_INTRA
					, bNonZeroCoeff
#endif
#if JVET_C0045_C0053_NO_NSST_FOR_TS
					, iNonZeroCoeffNonTs
#endif
				);

				ruiBits += m_pcEntropyCoder->getNumberOfWrittenBits();
			}
		}





		/**
		* \brief Generate half-sample interpolated block
		*
		* \param pattern Reference picture ROI
		* \param biPred    Flag indicating whether block is for biprediction
		*/
		Void TEncSearch::xExtDIFUpSamplingH(TComPattern* pattern)
		{
			Int width = pattern->getROIYWidth();
			Int height = pattern->getROIYHeight();
			Int srcStride = pattern->getPatternLStride();

			Int intStride = m_filteredBlockTmp[0].getStride(COMPONENT_Y);
			Int dstStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);
			Pel *intPtr;
			Pel *dstPtr;
			Int filterSize = NTAPS_LUMA;
			Int halfFilterSize = (filterSize >> 1);
			Pel *srcPtr = pattern->getROIY() - halfFilterSize*srcStride - 1;

			const ChromaFormat chFmt = m_filteredBlock[0][0].getChromaFormat();

			m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[0].getAddr(COMPONENT_Y), intStride, width + 1, height + filterSize, 0, false, chFmt, pattern->getBitDepthY());
			m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[2].getAddr(COMPONENT_Y), intStride, width + 1, height + filterSize,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				2 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				2,
#endif
				false, chFmt, pattern->getBitDepthY());

			intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + halfFilterSize * intStride + 1;
			dstPtr = m_filteredBlock[0][0].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 0, height + 0, 0, false, true, chFmt, pattern->getBitDepthY());

			intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride + 1;
			dstPtr = m_filteredBlock[2][0].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 0, height + 1,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				2 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				2,
#endif
				false, true, chFmt, pattern->getBitDepthY());

			intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
			dstPtr = m_filteredBlock[0][2].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 1, height + 0, 0, false, true, chFmt, pattern->getBitDepthY());

			intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
			dstPtr = m_filteredBlock[2][2].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 1, height + 1,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				2 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				2,
#endif
				false, true, chFmt, pattern->getBitDepthY());
		}





		/**
		* \brief Generate quarter-sample interpolated blocks
		*
		* \param pattern    Reference picture ROI
		* \param halfPelRef Half-pel mv
		* \param biPred     Flag indicating whether block is for biprediction
		*/
		Void TEncSearch::xExtDIFUpSamplingQ(TComPattern* pattern, TComMv halfPelRef
#if SIF_GAUSS
			, TComDataCU* pcCU
			, TComPic *pcRefPic
			, UInt uiPartAddr
#endif
		)
		{
			Int width = pattern->getROIYWidth();
			Int height = pattern->getROIYHeight();
			Int srcStride = pattern->getPatternLStride();

			Pel *srcPtr;
			Int intStride = m_filteredBlockTmp[0].getStride(COMPONENT_Y);
			Int dstStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);
			Pel *intPtr;
			Pel *dstPtr;
			Int filterSize = NTAPS_LUMA;

			Int halfFilterSize = (filterSize >> 1);

			Int extHeight = (halfPelRef.getVer() == 0) ? height + filterSize : height + filterSize - 1;

			const ChromaFormat chFmt = m_filteredBlock[0][0].getChromaFormat();

#if SIF_GAUSS
			Int nFilterTypeHor = 0;
			Int nFilterTypeVer = 0;
			TComMv mv(0, 0);
#endif
			// Horizontal filter 1/4
#if SIF_GAUSS
			mv.set(1, 0);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			srcPtr = pattern->getROIY() - halfFilterSize * srcStride - 1;
			intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y);
			if (halfPelRef.getVer() > 0)
			{
				srcPtr += srcStride;
			}
			if (halfPelRef.getHor() >= 0)
			{
				srcPtr += 1;
			}
			m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				1 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				1,
#endif
				false, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeHor
#endif
			);

			// Horizontal filter 3/4
#if SIF_GAUSS
			mv.set(3, 0);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			srcPtr = pattern->getROIY() - halfFilterSize*srcStride - 1;
			intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y);
			if (halfPelRef.getVer() > 0)
			{
				srcPtr += srcStride;
			}
			if (halfPelRef.getHor() > 0)
			{
				srcPtr += 1;
			}
			m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				3 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				3,
#endif
				false, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeHor
#endif
			);

			// Generate @ 1,1
#if SIF_GAUSS
			mv.set(1, 1);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
			dstPtr = m_filteredBlock[1][1].getAddr(COMPONENT_Y);
			if (halfPelRef.getVer() == 0)
			{
				intPtr += intStride;
			}
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				1 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				1,
#endif
				false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeVer
#endif
			);

			// Generate @ 3,1
#if SIF_GAUSS
			mv.set(1, 3);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
			dstPtr = m_filteredBlock[3][1].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				3 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				3,
#endif
				false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeVer
#endif
			);

			if (halfPelRef.getVer() != 0)
			{
				// Generate @ 2,1
#if SIF_GAUSS
				mv.set(1, 2);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
				dstPtr = m_filteredBlock[2][1].getAddr(COMPONENT_Y);
				if (halfPelRef.getVer() == 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					2 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					2,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);

				// Generate @ 2,3
#if SIF_GAUSS
				mv.set(3, 2);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
				dstPtr = m_filteredBlock[2][3].getAddr(COMPONENT_Y);
				if (halfPelRef.getVer() == 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					2 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					2,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);
			}
			else
			{
				// Generate @ 0,1
#if SIF_GAUSS
				mv.set(1, 0);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
				dstPtr = m_filteredBlock[0][1].getAddr(COMPONENT_Y);
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0, false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);

				// Generate @ 0,3
#if SIF_GAUSS
				mv.set(3, 0);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
				dstPtr = m_filteredBlock[0][3].getAddr(COMPONENT_Y);
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0, false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);
			}

			if (halfPelRef.getHor() != 0)
			{
				// Generate @ 1,2
#if SIF_GAUSS
				mv.set(2, 1);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
				dstPtr = m_filteredBlock[1][2].getAddr(COMPONENT_Y);
				if (halfPelRef.getHor() > 0)
				{
					intPtr += 1;
				}
				if (halfPelRef.getVer() >= 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					1 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					1,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);

				// Generate @ 3,2
#if SIF_GAUSS
				mv.set(2, 3);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
				dstPtr = m_filteredBlock[3][2].getAddr(COMPONENT_Y);
				if (halfPelRef.getHor() > 0)
				{
					intPtr += 1;
				}
				if (halfPelRef.getVer() > 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					3 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					3,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);
			}
			else
			{
				// Generate @ 1,0
#if SIF_GAUSS
				mv.set(0, 1);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride + 1;
				dstPtr = m_filteredBlock[1][0].getAddr(COMPONENT_Y);
				if (halfPelRef.getVer() >= 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					1 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					1,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);

				// Generate @ 3,0
#if SIF_GAUSS
				mv.set(0, 3);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
				pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
				intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride + 1;
				dstPtr = m_filteredBlock[3][0].getAddr(COMPONENT_Y);
				if (halfPelRef.getVer() > 0)
				{
					intPtr += intStride;
				}
				m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					3 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
					3,
#endif
					false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
					, 0
#endif
#if SIF_GAUSS
					, nFilterTypeVer
#endif
				);
			}

			// Generate @ 1,3
#if SIF_GAUSS
			mv.set(3, 1);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
			dstPtr = m_filteredBlock[1][3].getAddr(COMPONENT_Y);
			if (halfPelRef.getVer() == 0)
			{
				intPtr += intStride;
			}
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				1 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				1,
#endif
				false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeVer
#endif
			);

			// Generate @ 3,3
#if SIF_GAUSS
			mv.set(3, 3);
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
			mv <<= VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE;
#endif
			pcCU->getFilterType(&mv, uiPartAddr, nFilterTypeHor, nFilterTypeVer);
#endif
			intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize - 1) * intStride;
			dstPtr = m_filteredBlock[3][3].getAddr(COMPONENT_Y);
			m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height,
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
				3 << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE,
#else
				3,
#endif
				false, true, chFmt, pattern->getBitDepthY()
#if VCEG_AZ07_FRUC_MERGE
				, 0
#endif
#if SIF_GAUSS
				, nFilterTypeVer
#endif
			);
		}





		//! set wp tables
		Void  TEncSearch::setWpScalingDistParam(TComDataCU* pcCU, Int iRefIdx, RefPicList eRefPicListCur)
		{
			if (iRefIdx<0)
			{
				m_cDistParam.bApplyWeight = false;
				return;
			}

			TComSlice       *pcSlice = pcCU->getSlice();
			WPScalingParam  *wp0, *wp1;

			m_cDistParam.bApplyWeight = (pcSlice->getSliceType() == P_SLICE && pcSlice->testWeightPred()) || (pcSlice->getSliceType() == B_SLICE && pcSlice->testWeightBiPred());

			if (!m_cDistParam.bApplyWeight)
			{
				return;
			}

			Int iRefIdx0 = (eRefPicListCur == REF_PIC_LIST_0) ? iRefIdx : (-1);
			Int iRefIdx1 = (eRefPicListCur == REF_PIC_LIST_1) ? iRefIdx : (-1);

			getWpScaling(pcCU, iRefIdx0, iRefIdx1, wp0, wp1);

			if (iRefIdx0 < 0)
			{
				wp0 = NULL;
			}
			if (iRefIdx1 < 0)
			{
				wp1 = NULL;
			}

			m_cDistParam.wpCur = NULL;

			if (eRefPicListCur == REF_PIC_LIST_0)
			{
				m_cDistParam.wpCur = wp0;
			}
			else
			{
				m_cDistParam.wpCur = wp1;
			}
		}

#if COM16_C1016_AFFINE
		Bool TEncSearch::xEstimateAffineAMVP(TComDataCU* pcCU, TComYuv* pcOrgYuv, UInt uiPartIdx, RefPicList eRefPicList, Int iRefIdx, TComMv acMvPred[3], Distortion* puiDistBiP
#if JVECT_C0062_AFFINE_SIX_PARAM
			, Bool b6ParamAffine
#endif
		)
		{
#if JVECT_C0062_AFFINE_SIX_PARAM
			AffineAMVPInfo* pcAffineAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo(b6ParamAffine);
#else
			AffineAMVPInfo* pcAffineAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo();
#endif
			TComMv     acBestMv[3];
			Int        iBestIdx = 0;
			Distortion uiBestCost = std::numeric_limits<Distortion>::max();
			UInt       uiPartAddr = 0;
			Int        iRoiWidth, iRoiHeight;
			Int        i;

			pcCU->getPartIndexAndSize(uiPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);

			// Fill the MV Candidates
#if JVECT_C0062_AFFINE_SIX_PARAM
			pcCU->fillAffineMvpCand(uiPartIdx, uiPartAddr, eRefPicList, iRefIdx, pcAffineAMVPInfo, b6ParamAffine);
#else
			pcCU->fillAffineMvpCand(uiPartIdx, uiPartAddr, eRefPicList, iRefIdx, pcAffineAMVPInfo);
#endif
			assert(pcAffineAMVPInfo->iN > 0);

			// initialize Mvp index & Mvp
			m_cYuvPredTemp.clear();
			iBestIdx = 0;
			for (i = 0; i < pcAffineAMVPInfo->iN; i++)
			{
#if AFFINE_CLIP_REMOVE
				TComMv acMvTemp[3];
				memcpy(acMvTemp, pcAffineAMVPInfo->m_acMvCand[i], sizeof(TComMv) * 3);
#if JVECT_C0062_AFFINE_SIX_PARAM
				Distortion uiTmpCost = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, acMvTemp, i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight, b6ParamAffine);
#else
				Distortion uiTmpCost = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, acMvTemp, i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
#endif
#else
#if JVECT_C0062_AFFINE_SIX_PARAM
				Distortion uiTmpCost = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, pcAffineAMVPInfo->m_acMvCand[i], i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight, b6ParamAffine);
#else
				Distortion uiTmpCost = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, pcAffineAMVPInfo->m_acMvCand[i], i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
#endif
#endif
				if (uiBestCost > uiTmpCost)
				{
					uiBestCost = uiTmpCost;
					memcpy(acBestMv, pcAffineAMVPInfo->m_acMvCand[i], sizeof(TComMv) * 3);
					iBestIdx = i;
					(*puiDistBiP) = uiTmpCost;
				}
			}

#if !COM16_C806_LARGE_CTU
			m_cYuvPredTemp.clear();
#endif

			// Setting Best MVP
			memcpy(acMvPred, acBestMv, sizeof(TComMv) * 3);
			pcCU->setMVPIdxSubParts(iBestIdx, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(pcAffineAMVPInfo->iN, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
			return true;
		}

		Distortion TEncSearch::xGetAffineTemplateCost(TComDataCU* pcCU, UInt uiPartAddr, TComYuv* pcOrgYuv, TComYuv* pcTemplateCand, TComMv acMvCand[3], Int iMVPIdx, Int iMVPNum, RefPicList eRefPicList, Int iRefIdx, Int iSizeX, Int iSizeY
#if JVECT_C0062_AFFINE_SIX_PARAM
			, Bool b6Param
#endif
		)
		{
			Distortion uiCost = std::numeric_limits<Distortion>::max();

			TComPicYuv* pcPicYuvRef = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdx)->getPicYuvRec();

			// prediction pattern
			if (pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType() == P_SLICE)
			{
				xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvCand, iSizeX, iSizeY, pcTemplateCand, true, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
					, b6Param
#endif
				);
			}
			else
			{
				xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvCand, iSizeX, iSizeY, pcTemplateCand, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
					, b6Param
#endif
				);
			}

			if (pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType() == P_SLICE)
			{
				xWeightedPredictionUni(pcCU, pcTemplateCand, uiPartAddr, iSizeX, iSizeY, eRefPicList, pcTemplateCand, iRefIdx);
			}

			// calc distortion
			uiCost = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pcTemplateCand->getAddr(COMPONENT_Y, uiPartAddr), pcTemplateCand->getStride(COMPONENT_Y), pcOrgYuv->getAddr(COMPONENT_Y, uiPartAddr), pcOrgYuv->getStride(COMPONENT_Y), iSizeX, iSizeY, COMPONENT_Y, DF_SAD);
#if JVET_D0123_ME_CTX_LUT_BITS
			uiCost = uiCost + (UInt)m_pcRdCost->getCostSearch(m_auiMVPIdxCost[iMVPIdx][iMVPNum]);
#else
			uiCost = (UInt)m_pcRdCost->calcRdCost(m_auiMVPIdxCost[iMVPIdx][iMVPNum], uiCost, false, DF_SAD);
#endif
			return uiCost;
		}

		Void TEncSearch::xCopyAffineAMVPInfo(AffineAMVPInfo* pSrc, AffineAMVPInfo* pDst)
		{
			pDst->iN = pSrc->iN;
			for (Int i = 0; i < pSrc->iN; i++)
			{
				memcpy(pDst->m_acMvCand[i], pSrc->m_acMvCand[i], sizeof(TComMv) * 3);
			}
		}
#if ITH_C2_FXD_PT_AFF_PARAMS
		Void TEncSearch::solveEqual(Double(*dEqualCoeff)[TEncSearch::iEqualCoeffBufDimension], Int iOrder, Double* dAffinePara)
#else
		Void TEncSearch::solveEqual(Double** dEqualCoeff, Int iOrder, Double* dAffinePara)
#endif
		{
			// row echelon
			for (Int i = 1; i < iOrder; i++)
			{
				// find column max
				Double temp = fabs(dEqualCoeff[i][i - 1]);
				Int tempIdx = i;
				for (Int j = i + 1; j < iOrder + 1; j++)
				{
					if (fabs(dEqualCoeff[j][i - 1]) > temp)
					{
						temp = fabs(dEqualCoeff[j][i - 1]);
						tempIdx = j;
					}
				}

				// swap line
				if (tempIdx != i)
				{
					for (Int j = 0; j < iOrder + 1; j++)
					{
						dEqualCoeff[0][j] = dEqualCoeff[i][j];
						dEqualCoeff[i][j] = dEqualCoeff[tempIdx][j];
						dEqualCoeff[tempIdx][j] = dEqualCoeff[0][j];
					}
				}

				// elimination first column
				for (Int j = i + 1; j < iOrder + 1; j++)
				{
					for (Int k = i; k < iOrder + 1; k++)
					{
						dEqualCoeff[j][k] = dEqualCoeff[j][k] - dEqualCoeff[i][k] * dEqualCoeff[j][i - 1] / dEqualCoeff[i][i - 1];
					}
				}
			}

			dAffinePara[iOrder - 1] = dEqualCoeff[iOrder][iOrder] / dEqualCoeff[iOrder][iOrder - 1];
			for (Int i = iOrder - 2; i >= 0; i--)
			{
				Double temp = 0;
				for (Int j = i + 1; j < iOrder; j++)
				{
					temp += dEqualCoeff[i + 1][j] * dAffinePara[j];
				}
				dAffinePara[i] = (dEqualCoeff[i + 1][iOrder] - temp) / dEqualCoeff[i + 1][i];
			}
		}

		Void TEncSearch::predAffineInterSearch(TComDataCU* pcCU, TComYuv* pcOrgYuv, Int iPartIdx, UInt& ruiLastMode, Distortion& ruiAffineCost, TComMv cHevcMv[2][33]
#if JVECT_C0062_AFFINE_SIX_PARAM
			, UInt uiAffineParam, UInt& uiBestAffineParam, TComMv cMVBest[2][3], Int iBestRefIdx[2], TComMv bestMVD[2][3], Int&bestDir, Int bestMVPIdx[2], Int bestMVPNum[2]
#if AFFINE_INHERIT_FOUR_PARAM
	, TComMv cMvAffineParam[2][33][3]
#endif
#endif
)
{
#if !JVECT_C0062_AFFINE_SIX_PARAM
			ruiAffineCost = std::numeric_limits<Distortion>::max();
#endif

#if !COM16_C806_LARGE_CTU
			for (UInt i = 0; i<NUM_REF_PIC_LIST_01; i++)
			{
				m_acYuvPred[i].clear();
			}
			m_cYuvPredTemp.clear();
#endif

			TComMv        cMvZero;
			TComMv        aacMvd[2][3];
			TComMv        aacMv[2][3];
			TComMv        cMvBi[2][3];
			TComMv        cMvTemp[2][33][3];

			Int           iNumPredDir = pcCU->getSlice()->isInterP() ? 1 : 2;

			// Mvp
			TComMv        cMvPred[2][33][3];
			TComMv        cMvPredBi[2][33][3];
			Int           aaiMvpIdxBi[2][33];
			Int           aaiMvpIdx[2][33];
			Int           aaiMvpNum[2][33];

			AffineAMVPInfo aacAffineAMVPInfo[2][33];

			Int           iRefIdx[2] = { 0,0 }; // If un-initialized, may cause SEGV in bi-directional prediction iterative stage.
			Int           iRefIdxBi[2];

			UInt          uiPartAddr;
			Int           iRoiWidth, iRoiHeight;

			UInt          uiMbBits[3] = { 1, 1, 0 };

			Int           iRefStart, iRefEnd;

#if JVET_C0024_QTBT
			PartSize      ePartSize = SIZE_2Nx2N;
#else
			PartSize      ePartSize = pcCU->getPartitionSize(0);
#endif

			Int           bestBiPRefIdxL1 = 0;
			Int           bestBiPMvpL1 = 0;
			UInt          biPDistTemp = MAX_INT;

			Distortion    uiCost[2] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };
			Distortion    uiCostBi = std::numeric_limits<Distortion>::max();
			Distortion    uiCostTemp;

			UInt          uiBits[3];
			UInt          uiBitsTemp;
			Distortion    bestBiPDist = std::numeric_limits<Distortion>::max();

			Distortion    uiCostTempL0[MAX_NUM_REF];
			for (Int iNumRef = 0; iNumRef < MAX_NUM_REF; iNumRef++)
			{
				uiCostTempL0[iNumRef] = std::numeric_limits<Distortion>::max();
			}
			Distortion    uiBitsTempL0[MAX_NUM_REF];

			TComMv        mvValidList1[4];
			Int           refIdxValidList1 = 0;
			UInt          bitsValidList1 = MAX_UINT;
			UInt          costValidList1 = MAX_UINT;
			TComMv        mvHevc[3];
#if JVET_D0123_ME_CTX_LUT_BITS
			loadCtxMeSbacCoder(pcCU);
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
			xGetBlkBits(ePartSize, pcCU->getSlice()->isInterP(), iPartIdx, ruiLastMode, uiMbBits, pcCU, uiPartAddr);
			xUpdateMvpIdxBits(pcCU);
#else
			xGetBlkBits(ePartSize, pcCU->getSlice()->isInterP(), iPartIdx, ruiLastMode, uiMbBits);

			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);
#endif

#if JVET_D0077_SAVE_LOAD_ENC_INFO
			UInt uiWidthIdx = g_aucConvertToBit[iRoiWidth];
			UInt uiHeightIdx = g_aucConvertToBit[iRoiHeight];
			UInt uiZorderIdx = pcCU->getZorderIdxInCtu();
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
			UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, pcCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx);
#else
			UChar saveLoadTag = getSaveLoadTag(uiZorderIdx, uiWidthIdx, uiHeightIdx);
#endif
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
			Bool bRecShape = (iRoiHeight / iRoiWidth) >= 2;
#endif

			// set Affine flag
			pcCU->setAffineFlagSubParts(true, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMergeFlagSubParts(false, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#if JVECT_C0062_AFFINE_SIX_PARAM
			pcCU->setAffineParamFlagSubParts(uiAffineParam == 1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
#endif

			// Uni-directional prediction
			for (Int iRefList = 0; iRefList < iNumPredDir; iRefList++)
			{
				RefPicList  eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

    for ( Int iRefIdxTemp = 0; iRefIdxTemp < pcCU->getSlice()->getNumRefIdx(eRefPicList); iRefIdxTemp++ )
    {
#if CLOSE_AFFINE_FOR_BG_REFERENCE
			if (pcCU->getSlice()->getRefPOC(eRefPicList, iRefIdxTemp) == -1)
			{
				continue;
			}
#endif

#if AFF_SKIP_ME
			if (uiAffineParam&&pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1)
			{
				if ((iRefList == 0 && iRefIdxTemp != g_FourRefIdx0) || (iRefList == 1 && iRefIdxTemp != g_FourRefIdx1))
				{
					continue;
				}
			}
#endif
      // Get RefIdx bits
      uiBitsTemp = uiMbBits[iRefList];
      if ( pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1 )
      {
#if JVET_D0123_ME_CTX_LUT_BITS
						uiBitsTemp += xRefFrameIdxBits(pcCU, iRefIdxTemp, (UInt)pcCU->getSlice()->getNumRefIdx(eRefPicList));
#else
						uiBitsTemp += iRefIdxTemp + 1;
						if (iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1)
						{
							uiBitsTemp--;
						}
#endif
					}

					// Do Affine AMVP
					xEstimateAffineAMVP(pcCU, pcOrgYuv, iPartIdx, eRefPicList, iRefIdxTemp, cMvPred[iRefList][iRefIdxTemp], &biPDistTemp
#if JVECT_C0062_AFFINE_SIX_PARAM
						, uiAffineParam == 1
#endif
					);
					aaiMvpIdx[iRefList][iRefIdxTemp] = pcCU->getMVPIdx(eRefPicList, uiPartAddr);
					aaiMvpNum[iRefList][iRefIdxTemp] = pcCU->getMVPNum(eRefPicList, uiPartAddr);

					// set hevc ME result as start search position when it is best than mvp
					for (Int i = 0; i<3; i++)
					{
						mvHevc[i] = cHevcMv[iRefList][iRefIdxTemp];
					}
					UInt uiCandCost = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, mvHevc, aaiMvpIdx[iRefList][iRefIdxTemp], AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp, iRoiWidth, iRoiHeight
#if JVECT_C0062_AFFINE_SIX_PARAM
						, uiAffineParam == 1
#endif
      );
#if AFFINE_INHERIT_FOUR_PARAM
			TComMv        mvFour[3];
			if (uiAffineParam)
			{
				for (Int i = 0; i < 3; i++)
				{
					mvFour[i] = cMvAffineParam[iRefList][iRefIdxTemp][i];
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					int iCurVal = mvFour[i].getHor();
					if (iCurVal != (iCurVal >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE))
					{
						if (mvFour[i].getHor()<0)
						{
							mvFour[i].setHor((mvFour[i].getHor() - 2) >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						}
						else
						{
							mvFour[i].setHor((mvFour[i].getHor() + 2) >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						}
					}
					iCurVal = mvFour[i].getVer();
					if (iCurVal != (iCurVal >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE))
					{
						if (mvFour[i].getVer() < 0)
						{
							mvFour[i].setVer((mvFour[i].getVer() - 2) >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						}
						else
						{
							mvFour[i].setVer((mvFour[i].getVer() + 2) >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						}
					}
#endif
				}
				UInt uiCandCost_Inherit = xGetAffineTemplateCost(pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, mvFour, aaiMvpIdx[iRefList][iRefIdxTemp], AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp, iRoiWidth, iRoiHeight
#if JVECT_C0062_AFFINE_SIX_PARAM
					, uiAffineParam == 1
#endif
					);
				if (uiCandCost_Inherit < uiCandCost)
				{
					uiCandCost = uiCandCost_Inherit;
					for (Int i = 0; i < 3; i++)
					{
						mvHevc[i] = mvFour[i];
					}
				}
			}

#endif

      if ( uiCandCost < biPDistTemp )
      {
        memcpy( cMvTemp[iRefList][iRefIdxTemp], mvHevc, sizeof(TComMv)*3 );
      }
      else
      {
        memcpy( cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], sizeof(TComMv)*3 );
      }

					// GPB list 1, save the best MvpIdx, RefIdx and Cost
					if (pcCU->getSlice()->getMvdL1ZeroFlag() && iRefList == 1 && biPDistTemp < bestBiPDist)
					{
						bestBiPDist = biPDistTemp;
						bestBiPMvpL1 = aaiMvpIdx[iRefList][iRefIdxTemp];
						bestBiPRefIdxL1 = iRefIdxTemp;
					}

					// Update bits
					uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];

					if (m_pcEncCfg->getFastMEForGenBLowDelayEnabled() && iRefList == 1)   // list 1
					{
						if (pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp) >= 0)
						{
							Int iList1ToList0Idx = pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp);
							memcpy(cMvTemp[1][iRefIdxTemp], cMvTemp[0][iList1ToList0Idx], sizeof(TComMv) * 3);
							uiCostTemp = uiCostTempL0[iList1ToList0Idx];
							/*first subtract the bit-rate part of the cost of the other list*/
#if JVET_D0123_ME_CTX_LUT_BITS
							uiCostTemp -= m_pcRdCost->getCostSearch(uiBitsTempL0[iList1ToList0Idx]);
#else
							uiCostTemp -= m_pcRdCost->getCost(uiBitsTempL0[iList1ToList0Idx]);
#endif
							/*correct the bit-rate part of the current ref*/
#if JVECT_C0062_AFFINE_SIX_PARAM
							Int numVerIdx = uiAffineParam ? 3 : 2;
							for (Int iVerIdx = 0; iVerIdx < numVerIdx; iVerIdx++)
#else
#if JVET_C0024_QTBT
							for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++) //Affine fix
#else
							for (Int iVerIdx = 0; iVerIdx<3; iVerIdx++)
#endif
#endif
							{
#if JVECT_C0062_AFFINE_SIX_PARAM
								Int iTmpIdx = (bRecShape && iVerIdx && uiAffineParam == 0) ? (iVerIdx + 1) : iVerIdx;
								m_pcRdCost->setPredictor(cMvPred[iRefList][iRefIdxTemp][iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
								TComPattern   tmpPattern;
								TComPattern*  pcPatternKey = &tmpPattern;
								pcPatternKey->m_pcCU = pcCU;
								pcPatternKey->m_uiPartAddr = uiPartAddr;
								pcPatternKey->m_mvPred = &cMvPred[iRefList][iRefIdxTemp][iTmpIdx];
								iCostScale = 0;
								uiBitsTemp += xMvdBits(pcPatternKey, cMvTemp[1][iRefIdxTemp][iTmpIdx].getHor(), cMvTemp[1][iRefIdxTemp][iTmpIdx].getVer());
#else
								uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
									cMvTemp[1][iRefIdxTemp][iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, cMvTemp[1][iRefIdxTemp][iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
									cMvTemp[1][iRefIdxTemp][iTmpIdx].getHor(), cMvTemp[1][iRefIdxTemp][iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
									, false
#endif
								);
#endif
#else
								m_pcRdCost->setPredictor(cMvPred[iRefList][iRefIdxTemp][iVerIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
								TComPattern   tmpPattern;
								TComPattern*  pcPatternKey = &tmpPattern;
								pcPatternKey->m_pcCU = pcCU;
								pcPatternKey->m_uiPartAddr = uiPartAddr;
								pcPatternKey->m_mvPred = &cMvPred[iRefList][iRefIdxTemp][iVerIdx];
								iCostScale = 0;
								uiBitsTemp += xMvdBits(pcPatternKey, cMvTemp[1][iRefIdxTemp][iVerIdx].getHor(), cMvTemp[1][iRefIdxTemp][iVerIdx].getVer());
#else
								uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
									cMvTemp[1][iRefIdxTemp][iVerIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, cMvTemp[1][iRefIdxTemp][iVerIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
									cMvTemp[1][iRefIdxTemp][iVerIdx].getHor(), cMvTemp[1][iRefIdxTemp][iVerIdx].getVer()
#endif
#if VCEG_AZ07_IMV
									, false
#endif
								);
#endif
#endif
							}
							/*calculate the correct cost*/
#if JVET_D0123_ME_CTX_LUT_BITS
							uiCostTemp += m_pcRdCost->getCostSearch(uiBitsTemp);
#else
							uiCostTemp += m_pcRdCost->getCost(uiBitsTemp);
#endif
						}
						else
						{
							xAffineMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if JVECT_C0062_AFFINE_SIX_PARAM
								, uiAffineParam
#endif
							);
						}
					}
					else
					{
						xAffineMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if JVECT_C0062_AFFINE_SIX_PARAM
							, uiAffineParam
#endif
						);
					}

					// Set best AMVP Index
#if JVECT_C0062_AFFINE_SIX_PARAM
					xCopyAffineAMVPInfo(pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo(uiAffineParam), &aacAffineAMVPInfo[iRefList][iRefIdxTemp]);
					xCheckBestAffineMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, uiAffineParam
#if JVET_D0123_ME_CTX_LUT_BITS
						, uiPartAddr
#endif
					);
#else
					xCopyAffineAMVPInfo(pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo(), &aacAffineAMVPInfo[iRefList][iRefIdxTemp]);
					xCheckBestAffineMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if JVET_D0123_ME_CTX_LUT_BITS
						, uiPartAddr
#endif
					);
#endif

					if (iRefList == 0)
					{
						uiCostTempL0[iRefIdxTemp] = uiCostTemp;
						uiBitsTempL0[iRefIdxTemp] = uiBitsTemp;
					}
					if (uiCostTemp < uiCost[iRefList])
					{
						uiCost[iRefList] = uiCostTemp;
						uiBits[iRefList] = uiBitsTemp; // storing for bi-prediction

													   // set best motion
						memcpy(aacMv[iRefList], cMvTemp[iRefList][iRefIdxTemp], sizeof(TComMv) * 3);
						iRefIdx[iRefList] = iRefIdxTemp;
					}

					if (iRefList == 1 && uiCostTemp < costValidList1 && pcCU->getSlice()->getList1IdxToList0Idx(iRefIdxTemp) < 0)
					{
						costValidList1 = uiCostTemp;
						bitsValidList1 = uiBitsTemp;

						// set motion
						memcpy(mvValidList1, cMvTemp[iRefList][iRefIdxTemp], sizeof(TComMv) * 3);
						refIdxValidList1 = iRefIdxTemp;
					}
				} // End refIdx loop
			} // end Uni-prediction

#if AFFINE_INHERIT_FOUR_PARAM // save 4 parameter affine ME result for 6 parameter Affine ME
	::memcpy(cMvAffineParam, cMvTemp, sizeof(cMvTemp));
#endif

  // Bi-directional prediction
  if ( ( pcCU->getSlice()->isInterB() ) && ( pcCU->isBipredRestriction(iPartIdx) == false ) 
#if JVET_D0077_SAVE_LOAD_ENC_INFO
#if JVET_D0077_SAVE_LOAD_ENC_INFO_EXT
				&& !(saveLoadTag == LOAD_ENC_INFO && getSaveLoadInterDir(pcCU->getZorderIdxInCtu(), pcCU->getCURouteIdx(), uiWidthIdx, uiHeightIdx) != 3)
#else
				&& !(saveLoadTag == LOAD_ENC_INFO && getSaveLoadInterDir(uiWidthIdx, uiHeightIdx) != 3)
#endif
#endif
				)
			{
				// Set as best list0 and list1
				iRefIdxBi[0] = iRefIdx[0];
				iRefIdxBi[1] = iRefIdx[1];

				::memcpy(cMvBi, aacMv, sizeof(aacMv));
				::memcpy(cMvPredBi, cMvPred, sizeof(cMvPred));
				::memcpy(aaiMvpIdxBi, aaiMvpIdx, sizeof(aaiMvpIdx));

				UInt uiMotBits[2];

				if (pcCU->getSlice()->getMvdL1ZeroFlag()) // GPB, list 1 only use Mvp
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					xCopyAffineAMVPInfo(&aacAffineAMVPInfo[1][bestBiPRefIdxL1], pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo(uiAffineParam));
#else
					xCopyAffineAMVPInfo(&aacAffineAMVPInfo[1][bestBiPRefIdxL1], pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo());
#endif
					pcCU->setMVPIdxSubParts(bestBiPMvpL1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					aaiMvpIdxBi[1][bestBiPRefIdxL1] = bestBiPMvpL1;

					// Set Mv for list1
#if JVECT_C0062_AFFINE_SIX_PARAM
					TComMv* pcMvTemp = pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo(uiAffineParam)->m_acMvCand[bestBiPMvpL1];
#else
					TComMv* pcMvTemp = pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo()->m_acMvCand[bestBiPMvpL1];
#endif
					memcpy(cMvPredBi[1][bestBiPRefIdxL1], pcMvTemp, sizeof(TComMv) * 3);
					memcpy(cMvBi[1], pcMvTemp, sizeof(TComMv) * 3);
					memcpy(cMvTemp[1][bestBiPRefIdxL1], pcMvTemp, sizeof(TComMv) * 3);
					iRefIdxBi[1] = bestBiPRefIdxL1;

					// Get list1 prediction block
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[1], REF_PIC_LIST_1, 0
#if JVECT_C0062_AFFINE_SIX_PARAM
						, uiAffineParam
#endif
					);
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx);
					TComYuv* pcYuvPred = &m_acYuvPred[1];
#if JVET_E0052_DMVR
					motionCompensation(pcCU, pcYuvPred, true, REF_PIC_LIST_1, iPartIdx);
#else
					motionCompensation(pcCU, pcYuvPred, REF_PIC_LIST_1, iPartIdx);
#endif

					// Update bits
					uiMotBits[0] = uiBits[0] - uiMbBits[0];
					uiMotBits[1] = uiMbBits[1];
					if (pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 1)
					{
#if JVET_D0123_ME_CTX_LUT_BITS
						uiMotBits[1] += xRefFrameIdxBits(pcCU, bestBiPRefIdxL1, (UInt)pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1));
#else
						uiMotBits[1] += bestBiPRefIdxL1 + 1;
						if (bestBiPRefIdxL1 == pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) - 1) uiMotBits[1]--;
#endif
					}
					uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][AMVP_MAX_NUM_CANDS];
					uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
				}
				else
				{
					uiMotBits[0] = uiBits[0] - uiMbBits[0];
					uiMotBits[1] = uiBits[1] - uiMbBits[1];
					uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
				}

				// 4-times iteration (default)
				Int iNumIter = 4;
				// fast encoder setting or GPB: only one iteration
				if (m_pcEncCfg->getUseFastEnc() || pcCU->getSlice()->getMvdL1ZeroFlag())
				{
					iNumIter = 1;
				}

				for (Int iIter = 0; iIter < iNumIter; iIter++)
				{
					// Set RefList
					Int iRefList = iIter % 2;
					if (m_pcEncCfg->getUseFastEnc()) // iNumIter = 1, reset the search list
					{
						if (uiCost[0] <= uiCost[1])
						{
							iRefList = 1;
						}
						else
						{
							iRefList = 0;
						}
					}
					else if (iIter == 0)
					{
						iRefList = 0;
					}

					// First iterate, get prediction block of opposite direction
					if (iIter == 0 && !pcCU->getSlice()->getMvdL1ZeroFlag())
					{
						pcCU->setAllAffineMv(uiPartAddr, iPartIdx, aacMv[1 - iRefList], RefPicList(1 - iRefList), 0
#if JVECT_C0062_AFFINE_SIX_PARAM
							, uiAffineParam
#endif
						);
						pcCU->getCUMvField(RefPicList(1 - iRefList))->setAllRefIdx(iRefIdx[1 - iRefList], ePartSize, uiPartAddr, 0, iPartIdx);
						TComYuv*  pcYuvPred = &m_acYuvPred[1 - iRefList];
#if JVET_E0052_DMVR
						motionCompensation(pcCU, pcYuvPred, true, RefPicList(1 - iRefList), iPartIdx);
#else
						motionCompensation(pcCU, pcYuvPred, RefPicList(1 - iRefList), iPartIdx);
#endif
					}

					RefPicList eRefPicList = (iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
					if (pcCU->getSlice()->getMvdL1ZeroFlag()) // GPB, fix List 1, search List 0
					{
						iRefList = 0;
						eRefPicList = REF_PIC_LIST_0;
					}

					Bool bChanged = false;

					iRefStart = 0;
					iRefEnd = pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1;

      for ( Int iRefIdxTemp = iRefStart; iRefIdxTemp <= iRefEnd; iRefIdxTemp++ )
      {
#if CLOSE_AFFINE_FOR_BG_REFERENCE
				if (pcCU->getSlice()->getRefPOC(eRefPicList, iRefIdxTemp) == -1)
				{
					continue;
				}
#endif

#if AFF_SKIP_ME
				if (uiAffineParam&&iRefEnd > 0)
				{
					if ((iRefList == 0 && iRefIdxTemp != g_FourRefIdx0) || (iRefList == 1 && iRefIdxTemp != g_FourRefIdx1))
					{
						continue;
					}
				}
#endif
        // update bits
        uiBitsTemp = uiMbBits[2] + uiMotBits[1-iRefList];
        if ( pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1 )
        {
#if JVET_D0123_ME_CTX_LUT_BITS
							uiBitsTemp += xRefFrameIdxBits(pcCU, iRefIdxTemp, (UInt)pcCU->getSlice()->getNumRefIdx(eRefPicList));
#else
							uiBitsTemp += iRefIdxTemp + 1;
							if (iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList) - 1) uiBitsTemp--;
#endif
						}
						uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];

						// call Affine ME
#if JVECT_C0062_AFFINE_SIX_PARAM
						xAffineMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, uiAffineParam, true);
						xCopyAffineAMVPInfo(&aacAffineAMVPInfo[iRefList][iRefIdxTemp], pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo(uiAffineParam));
						xCheckBestAffineMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, uiAffineParam
#if JVET_D0123_ME_CTX_LUT_BITS
							, uiPartAddr
#endif
						);
#else
						xAffineMotionEstimation(pcCU, pcOrgYuv, iPartIdx, eRefPicList, cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, true);
						xCopyAffineAMVPInfo(&aacAffineAMVPInfo[iRefList][iRefIdxTemp], pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo());
						xCheckBestAffineMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
#if JVET_D0123_ME_CTX_LUT_BITS
							, uiPartAddr
#endif
						);
#endif

						if (uiCostTemp < uiCostBi)
						{
							bChanged = true;
							memcpy(cMvBi[iRefList], cMvTemp[iRefList][iRefIdxTemp], sizeof(TComMv) * 3);
							iRefIdxBi[iRefList] = iRefIdxTemp;

							uiCostBi = uiCostTemp;
							uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2] - uiMotBits[1 - iRefList];
							uiBits[2] = uiBitsTemp;

							if (iNumIter != 1) // MC for next iter
							{
								//  Set motion
								pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[iRefList], eRefPicList, 0
#if JVECT_C0062_AFFINE_SIX_PARAM
									, uiAffineParam
#endif
								);
								pcCU->getCUMvField(eRefPicList)->setAllRefIdx(iRefIdxBi[iRefList], ePartSize, uiPartAddr, 0, iPartIdx);

								TComYuv* pcYuvPred = &m_acYuvPred[iRefList];
#if JVET_E0052_DMVR
								motionCompensation(pcCU, pcYuvPred, true, eRefPicList, iPartIdx);
#else
								motionCompensation(pcCU, pcYuvPred, eRefPicList, iPartIdx);
#endif
							}
						}
					} // for loop-iRefIdxTemp

					if (!bChanged)
					{
						if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])
						{
#if JVECT_C0062_AFFINE_SIX_PARAM
							if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])//ZF ToDo
							{
								xCopyAffineAMVPInfo(&aacAffineAMVPInfo[0][iRefIdxBi[0]], pcCU->getCUMvField(REF_PIC_LIST_0)->getAffineAMVPInfo(uiAffineParam));

								xCheckBestAffineMVP(pcCU, REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], uiBits[2], uiCostBi, uiAffineParam
#if JVET_D0123_ME_CTX_LUT_BITS
									, uiPartAddr
#endif
								);

								if (!pcCU->getSlice()->getMvdL1ZeroFlag())
								{
									xCopyAffineAMVPInfo(&aacAffineAMVPInfo[1][iRefIdxBi[1]], pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo(uiAffineParam));

									xCheckBestAffineMVP(pcCU, REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], uiBits[2], uiCostBi, uiAffineParam
#if JVET_D0123_ME_CTX_LUT_BITS
										, uiPartAddr
#endif
									);

								}
							}
#else
							xCopyAffineAMVPInfo(&aacAffineAMVPInfo[0][iRefIdxBi[0]], pcCU->getCUMvField(REF_PIC_LIST_0)->getAffineAMVPInfo());
							xCheckBestAffineMVP(pcCU, REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], uiBits[2], uiCostBi
#if JVET_D0123_ME_CTX_LUT_BITS
								, uiPartAddr
#endif
							);

							if (!pcCU->getSlice()->getMvdL1ZeroFlag())
							{
								xCopyAffineAMVPInfo(&aacAffineAMVPInfo[1][iRefIdxBi[1]], pcCU->getCUMvField(REF_PIC_LIST_1)->getAffineAMVPInfo());
								xCheckBestAffineMVP(pcCU, REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], uiBits[2], uiCostBi
#if JVET_D0123_ME_CTX_LUT_BITS
									, uiPartAddr
#endif
								);
							}
#endif
						}
						break;
					}
				} // for loop-iter
			} // if (B_SLICE)

			  // Clear Motion Field
			pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd(cMvZero, ePartSize, uiPartAddr, 0, iPartIdx);

			pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPIdxSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			pcCU->setMVPNumSubParts(-1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

#if JVECT_C0062_AFFINE_SIX_PARAM
			pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(-1, ePartSize, uiPartAddr, 0, iPartIdx);
			pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(-1, ePartSize, uiPartAddr, 0, iPartIdx);
#endif

			// Set Motion Field
			memcpy(aacMv[1], mvValidList1, sizeof(TComMv) * 3);
			iRefIdx[1] = refIdxValidList1;
			uiBits[1] = bitsValidList1;
			uiCost[1] = costValidList1;

			// Affine ME result set
#if JVECT_C0062_AFFINE_SIX_PARAM
#if AFF_SKIP_ME
	if (!uiAffineParam)
	{
		g_FourRefIdx0 = iRefIdx[0];
		g_FourRefIdx1 = iRefIdx[1];
	}
#endif
  if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1] && uiCostBi < ruiAffineCost) // Bi
#else
			if (uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1]) // Bi
#endif
			{
				ruiLastMode = 2;
				ruiAffineCost = uiCostBi;

#if JVECT_C0062_AFFINE_SIX_PARAM
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[0], REF_PIC_LIST_0, 0, uiAffineParam);
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[1], REF_PIC_LIST_1, 0, uiAffineParam);
#else
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[0], REF_PIC_LIST_0, 0);
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, cMvBi[1], REF_PIC_LIST_1, 0);
#endif
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iRefIdxBi[0], ePartSize, uiPartAddr, 0, iPartIdx);
				pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx);

#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int iVerIdx = 0; iVerIdx < 3; iVerIdx++)
#else
				for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++)
#endif
				{
					aacMvd[0][iVerIdx] = cMvBi[0][iVerIdx] - cMvPredBi[0][iRefIdxBi[0]][iVerIdx];
					aacMvd[1][iVerIdx] = cMvBi[1][iVerIdx] - cMvPredBi[1][iRefIdxBi[1]][iVerIdx];
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				uiBestAffineParam = uiAffineParam;
				cMVBest[0][0] = cMvBi[0][0];
				cMVBest[0][1] = cMvBi[0][1];
				cMVBest[0][2] = cMvBi[0][2];

				cMVBest[1][0] = cMvBi[1][0];
				cMVBest[1][1] = cMvBi[1][1];
				cMVBest[1][2] = cMvBi[1][2];

				iBestRefIdx[0] = iRefIdxBi[0];
				iBestRefIdx[1] = iRefIdxBi[1];

				bestMVD[0][0] = aacMvd[0][0];
				bestMVD[0][1] = aacMvd[0][1];
				bestMVD[0][2] = aacMvd[0][2];

				bestMVD[1][0] = aacMvd[1][0];
				bestMVD[1][1] = aacMvd[1][1];
				bestMVD[1][2] = aacMvd[1][2];

				bestDir = 3;
				bestMVPIdx[0] = aaiMvpIdxBi[0][iRefIdxBi[0]];
				bestMVPIdx[1] = aaiMvpIdxBi[1][iRefIdxBi[1]];
				bestMVPNum[0] = aaiMvpNum[0][iRefIdxBi[0]];
				bestMVPNum[1] = aaiMvpNum[1][iRefIdxBi[1]];
#endif

				pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, aacMvd[0], REF_PIC_LIST_0, 0);
				pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, aacMvd[1], REF_PIC_LIST_1, 0);

				pcCU->setInterDirSubParts(3, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

				pcCU->setMVPIdxSubParts(aaiMvpIdxBi[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPNumSubParts(aaiMvpNum[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPIdxSubParts(aaiMvpIdxBi[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPNumSubParts(aaiMvpNum[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			}
#if JVECT_C0062_AFFINE_SIX_PARAM
			else if (uiCost[0] <= uiCost[1] && uiCost[0] < ruiAffineCost) // List 0
#else
			else if (uiCost[0] <= uiCost[1]) // List 0
#endif
			{
				ruiLastMode = 0;
				ruiAffineCost = uiCost[0];

#if JVECT_C0062_AFFINE_SIX_PARAM
				uiBestAffineParam = uiAffineParam;
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, aacMv[0], REF_PIC_LIST_0, 0, uiAffineParam);
#else
				pcCU->setAllAffineMv(uiPartAddr, iPartIdx, aacMv[0], REF_PIC_LIST_0, 0);
#endif
				pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx(iRefIdx[0], ePartSize, uiPartAddr, 0, iPartIdx);

#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int iVerIdx = 0; iVerIdx < 3; iVerIdx++)
#else
				for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++)
#endif
				{
					aacMvd[0][iVerIdx] = aacMv[0][iVerIdx] - cMvPred[0][iRefIdx[0]][iVerIdx];
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				cMVBest[0][0] = aacMv[0][0];
				cMVBest[0][1] = aacMv[0][1];
				cMVBest[0][2] = aacMv[0][2];

				iBestRefIdx[0] = iRefIdx[0];

				bestMVD[0][0] = aacMvd[0][0];
				bestMVD[0][1] = aacMvd[0][1];
				bestMVD[0][2] = aacMvd[0][2];

				bestDir = 1;
				bestMVPIdx[0] = aaiMvpIdx[0][iRefIdx[0]];
				bestMVPNum[0] = aaiMvpNum[0][iRefIdx[0]];
#endif

				pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, aacMvd[0], REF_PIC_LIST_0, 0);

				pcCU->setInterDirSubParts(1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

				pcCU->setMVPIdxSubParts(aaiMvpIdx[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				pcCU->setMVPNumSubParts(aaiMvpNum[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
			}
			else
#if JVECT_C0062_AFFINE_SIX_PARAM
				if (uiCost[1] < ruiAffineCost)
#endif
				{
					ruiLastMode = 1;
					ruiAffineCost = uiCost[1];

#if JVECT_C0062_AFFINE_SIX_PARAM
					uiBestAffineParam = uiAffineParam;
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, aacMv[1], REF_PIC_LIST_1, 0, uiAffineParam);
#else
					pcCU->setAllAffineMv(uiPartAddr, iPartIdx, aacMv[1], REF_PIC_LIST_1, 0);
#endif
					pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx(iRefIdx[1], ePartSize, uiPartAddr, 0, iPartIdx);

#if JVECT_C0062_AFFINE_SIX_PARAM
					for (Int iVerIdx = 0; iVerIdx < 3; iVerIdx++)
#else
					for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++)
#endif
					{
						aacMvd[1][iVerIdx] = aacMv[1][iVerIdx] - cMvPred[1][iRefIdx[1]][iVerIdx];
					}

#if JVECT_C0062_AFFINE_SIX_PARAM
					cMVBest[1][0] = aacMv[1][0];
					cMVBest[1][1] = aacMv[1][1];
					cMVBest[1][2] = aacMv[1][2];

					iBestRefIdx[1] = iRefIdx[1];

					bestMVD[1][0] = aacMvd[1][0];
					bestMVD[1][1] = aacMvd[1][1];
					bestMVD[1][2] = aacMvd[1][2];

					bestDir = 2;
					bestMVPIdx[1] = aaiMvpIdx[1][iRefIdx[1]];
					bestMVPNum[1] = aaiMvpNum[1][iRefIdx[1]];
#endif

					pcCU->setAllAffineMvd(uiPartAddr, iPartIdx, aacMvd[1], REF_PIC_LIST_1, 0);

					pcCU->setInterDirSubParts(2, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

					pcCU->setMVPIdxSubParts(aaiMvpIdx[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
					pcCU->setMVPNumSubParts(aaiMvpNum[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
				}

			return;
		}

		Void TEncSearch::xCheckBestAffineMVP(TComDataCU* pcCU, RefPicList eRefPicList, TComMv acMv[3], TComMv acMvPred[3], Int& riMVPIdx, UInt& ruiBits, Distortion& ruiCost
#if JVECT_C0062_AFFINE_SIX_PARAM
			, Bool b6Param
#endif
#if JVET_D0123_ME_CTX_LUT_BITS
			, UInt uiPartAddr
#endif
		)
		{
#if JVECT_C0062_AFFINE_SIX_PARAM
			AffineAMVPInfo* pcAffineAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo(b6Param);
#else
			AffineAMVPInfo* pcAffineAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAffineAMVPInfo();
#endif

			if (pcAffineAMVPInfo->iN < 2)
			{
				return;
			}

			m_pcRdCost->getMotionCost(true, 0, pcCU->getCUTransquantBypass(0));
			m_pcRdCost->setCostScale(0);

			Int iBestMVPIdx = riMVPIdx;
#if JVECT_C0062_AFFINE_SIX_PARAM
			Int numVerIdx = (b6Param == 1) ? 3 : 2;
			Int width = pcCU->getWidth(0);
			Int height = pcCU->getHeight(0);
			Bool bRecShape = (height / width >= 2);
#endif

#if JVET_D0123_ME_CTX_LUT_BITS
			TComPattern   tmpPattern;
			TComPattern*  pcPatternKey = &tmpPattern;
			pcPatternKey->m_pcCU = pcCU;
			pcPatternKey->m_uiPartAddr = uiPartAddr;
			iCostScale = 0;
#endif

			// Get origin MV bits
			Int iOrgMvBits = 0;
#if JVECT_C0062_AFFINE_SIX_PARAM
			for (Int iVerIdx = 0; iVerIdx < numVerIdx; iVerIdx++)
#else
			for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++)
#endif
			{
#if JVECT_C0062_AFFINE_SIX_PARAM
				Int iTmpIdx = (bRecShape && iVerIdx && b6Param == 0) ? (iVerIdx + 1) : iVerIdx;
#if JVET_D0123_ME_CTX_LUT_BITS
				pcPatternKey->m_mvPred = &acMvPred[iTmpIdx];
				iOrgMvBits += xMvdBits(pcPatternKey, acMv[iTmpIdx].getHor(), acMv[iTmpIdx].getVer());
#else
				m_pcRdCost->setPredictor(acMvPred[iTmpIdx]);
				iOrgMvBits += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMv[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMv[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMv[iTmpIdx].getHor(), acMv[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
				);
#endif
#else
				m_pcRdCost->setPredictor(acMvPred[iVerIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
				pcPatternKey->m_mvPred = &acMvPred[iVerIdx];
				iOrgMvBits += xMvdBits(pcPatternKey, acMv[iVerIdx].getHor(), acMv[iVerIdx].getVer());
#else
				iOrgMvBits += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMv[iVerIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMv[iVerIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMv[iVerIdx].getHor(), acMv[iVerIdx].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
				);
#endif
#endif
			}
			iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];

			Int iBestMvBits = iOrgMvBits;
			for (Int iMVPIdx = 0; iMVPIdx < pcAffineAMVPInfo->iN; iMVPIdx++)
			{
				if (iMVPIdx == riMVPIdx)
				{
					continue;
				}

				Int iMvBits = 0;
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int iVerIdx = 0; iVerIdx < numVerIdx; iVerIdx++)
#else
				for (Int iVerIdx = 0; iVerIdx<2; iVerIdx++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (bRecShape && iVerIdx && b6Param == 0) ? (iVerIdx + 1) : iVerIdx;
					m_pcRdCost->setPredictor(pcAffineAMVPInfo->m_acMvCand[iMVPIdx][iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
					pcPatternKey->m_mvPred = &pcAffineAMVPInfo->m_acMvCand[iMVPIdx][iTmpIdx];
					iMvBits += xMvdBits(pcPatternKey, acMv[iTmpIdx].getHor(), acMv[iTmpIdx].getVer());
#else
					iMvBits += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMv[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMv[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMv[iTmpIdx].getHor(), acMv[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
					);
#endif
#else
					m_pcRdCost->setPredictor(pcAffineAMVPInfo->m_acMvCand[iMVPIdx][iVerIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
					pcPatternKey->m_mvPred = &pcAffineAMVPInfo->m_acMvCand[iMVPIdx][iVerIdx];
					iMvBits += xMvdBits(pcPatternKey, acMv[iVerIdx].getHor(), acMv[iVerIdx].getVer());
#else
					iMvBits += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMv[iVerIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMv[iVerIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMv[iVerIdx].getHor(), acMv[iVerIdx].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
					);
#endif
#endif
				}
				iMvBits += m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];

				if (iMvBits < iBestMvBits)
				{
					iBestMvBits = iMvBits;
					iBestMVPIdx = iMVPIdx;
				}
			}

			if (iBestMVPIdx != riMVPIdx)  // if changed
			{
				memcpy(acMvPred, pcAffineAMVPInfo->m_acMvCand[iBestMVPIdx], sizeof(TComMv) * 3);
				riMVPIdx = iBestMVPIdx;
				UInt uiOrgBits = ruiBits;
				ruiBits = uiOrgBits - iOrgMvBits + iBestMvBits;
#if JVET_D0123_ME_CTX_LUT_BITS
				ruiCost = (ruiCost - m_pcRdCost->getCostSearch(uiOrgBits)) + m_pcRdCost->getCostSearch(ruiBits);
#else
				ruiCost = (ruiCost - m_pcRdCost->getCost(uiOrgBits)) + m_pcRdCost->getCost(ruiBits);
#endif
			}
		}
#if ITH_C2_FXD_PT_AFF_PARAMS
#if ITH_C3_HORIZ_SOBEL_FILT==C_FXN
		static inline Void scaledHorizontalSobelFilter
			(
				Pel *const pPred,
				const Int iPredStride,
				Int *const piDerivate,
				const Int iDerivateBufStride,
				const Int iWidth,
				const Int iHeight
				)
		{
			for (Int j = 1; j < iHeight - 1; j++)
			{
				for (Int k = 1; k < iWidth - 1; k++)
				{
					Int iCenter = j * iPredStride + k;

					piDerivate[j * iDerivateBufStride + k] =
						pPred[iCenter + 1 - iPredStride] -
						pPred[iCenter - 1 - iPredStride] +
						(pPred[iCenter + 1] * 2) -
						(pPred[iCenter - 1] * 2) +
						pPred[iCenter + 1 + iPredStride] -
						pPred[iCenter - 1 + iPredStride];
				}

				piDerivate[j * iDerivateBufStride] = piDerivate[j * iDerivateBufStride + 1];
				piDerivate[j * iDerivateBufStride + iWidth - 1] = piDerivate[j * iDerivateBufStride + iWidth - 2];
			}

			piDerivate[0] = piDerivate[iDerivateBufStride + 1];
			piDerivate[iWidth - 1] = piDerivate[iDerivateBufStride + iWidth - 2];
			piDerivate[(iHeight - 1) * iDerivateBufStride] = piDerivate[(iHeight - 2) * iDerivateBufStride + 1];
			piDerivate[(iHeight - 1) * iDerivateBufStride + iWidth - 1] = piDerivate[(iHeight - 2) * iDerivateBufStride + (iWidth - 2)];

			for (Int j = 1; j < iWidth - 1; j++)
			{
				piDerivate[j] = piDerivate[iDerivateBufStride + j];
				piDerivate[(iHeight - 1) * iDerivateBufStride + j] = piDerivate[(iHeight - 2) * iDerivateBufStride + j];
			}
		}
#else
		static inline Void scaledHorizontalSobelFilter
			(
				Pel *const pPred,
				const Int iPredStride,
				Int *const piDerivate,
				const Int iDerivateBufStride,
				const Int iWidth,
				const Int iHeight
				)
		{
			__m128i mmPred[4];
			__m128i mm2xPred[2];
			__m128i mmIntermediates[4];
			__m128i mmDerivate[2];

			assert(!(iHeight % 2));
			assert(!(iWidth % 4));

			/* Derivates of the rows and columns at the boundary are done at the end of this function */
			/* The value of iCol and iRow indicate the columns and rows for which the derivates have already been ccomputed */
			for (Int iCol = 1; (iCol + 2) < iWidth; iCol += 2)
			{
				mmPred[0] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[0 * iPredStride + iCol - 1]));
				mmPred[1] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[1 * iPredStride + iCol - 1]));

				mmPred[0] = _mm_cvtepi16_epi32(mmPred[0]);
				mmPred[1] = _mm_cvtepi16_epi32(mmPred[1]);

				for (Int iRow = 1; iRow < (iHeight - 1); iRow += 2)
				{
					mmPred[2] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[(iRow + 1) * iPredStride + iCol - 1]));
					mmPred[3] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[(iRow + 2) * iPredStride + iCol - 1]));

					mmPred[2] = _mm_cvtepi16_epi32(mmPred[2]);
					mmPred[3] = _mm_cvtepi16_epi32(mmPred[3]);

					mm2xPred[0] = _mm_slli_epi32(mmPred[1], 1);
					mm2xPred[1] = _mm_slli_epi32(mmPred[2], 1);

					mmIntermediates[0] = _mm_add_epi32(mm2xPred[0], mmPred[0]);
					mmIntermediates[2] = _mm_add_epi32(mm2xPred[1], mmPred[1]);

					mmIntermediates[0] = _mm_add_epi32(mmIntermediates[0], mmPred[2]);
					mmIntermediates[2] = _mm_add_epi32(mmIntermediates[2], mmPred[3]);

					mmPred[0] = mmPred[2];
					mmPred[1] = mmPred[3];

					mmIntermediates[1] = _mm_srli_si128(mmIntermediates[0], 8);
					mmIntermediates[3] = _mm_srli_si128(mmIntermediates[2], 8);

					mmDerivate[0] = _mm_sub_epi32(mmIntermediates[1], mmIntermediates[0]);
					mmDerivate[1] = _mm_sub_epi32(mmIntermediates[3], mmIntermediates[2]);

					_mm_storel_epi64(reinterpret_cast<__m128i *> (&piDerivate[iCol + (iRow + 0) * iDerivateBufStride]), mmDerivate[0]);
					_mm_storel_epi64(reinterpret_cast<__m128i *> (&piDerivate[iCol + (iRow + 1) * iDerivateBufStride]), mmDerivate[1]);
				}
			}

			for (Int j = 1; j < (iHeight - 1); j++)
			{
				piDerivate[j * iDerivateBufStride] = piDerivate[j * iDerivateBufStride + 1];
				piDerivate[j * iDerivateBufStride + (iWidth - 1)] = piDerivate[j * iDerivateBufStride + (iWidth - 2)];
			}

			memcpy
				(
					piDerivate,
					piDerivate + iDerivateBufStride,
					iWidth * sizeof(piDerivate[0])
					);

			memcpy
				(
					piDerivate + (iHeight - 1) * iDerivateBufStride,
					piDerivate + (iHeight - 2) * iDerivateBufStride,
					iWidth * sizeof(piDerivate[0])
					);
		}
#endif

#if ITH_C3_VERT_SOBEL_FILT==C_FXN
		static inline Void scaledVerticalSobelFilter
			(
				Pel *const pPred,
				const Int iPredStride,
				Int *const piDerivate,
				const Int iDerivateBufStride,
				const Int iWidth,
				const Int iHeight
				)
		{
			for (Int k = 1; k < iWidth - 1; k++)
			{
				for (Int j = 1; j < iHeight - 1; j++)
				{
					Int iCenter = j * iPredStride + k;

					piDerivate[j * iDerivateBufStride + k] =
						pPred[iCenter + iPredStride - 1] -
						pPred[iCenter - iPredStride - 1] +
						(pPred[iCenter + iPredStride] * 2) -
						(pPred[iCenter - iPredStride] * 2) +
						pPred[iCenter + iPredStride + 1] -
						pPred[iCenter - iPredStride + 1];
				}

				piDerivate[k] = piDerivate[iDerivateBufStride + k];
				piDerivate[(iHeight - 1) * iDerivateBufStride + k] = piDerivate[(iHeight - 2) * iDerivateBufStride + k];
			}

			piDerivate[0] = piDerivate[iDerivateBufStride + 1];
			piDerivate[iWidth - 1] = piDerivate[iDerivateBufStride + iWidth - 2];
			piDerivate[(iHeight - 1) * iDerivateBufStride] = piDerivate[(iHeight - 2) * iDerivateBufStride + 1];
			piDerivate[(iHeight - 1) * iDerivateBufStride + iWidth - 1] = piDerivate[(iHeight - 2) * iDerivateBufStride + (iWidth - 2)];

			for (Int j = 1; j < iHeight - 1; j++)
			{
				piDerivate[j * iDerivateBufStride] = piDerivate[j * iDerivateBufStride + 1];
				piDerivate[j * iDerivateBufStride + iWidth - 1] = piDerivate[j * iDerivateBufStride + iWidth - 2];
			}
		}
#else
		static Void scaledVerticalSobelFilter
			(
				Pel *const pPred,
				const Int iPredStride,
				Int *const piDerivate,
				const Int iDerivateBufStride,
				const Int iWidth,
				const Int iHeight
				)
		{
			__m128i mmPred[4];
			__m128i mmIntermediates[6];
			__m128i mmDerivate[2];

			assert(!(iHeight % 2));
			assert(!(iWidth % 4));

			/* Derivates of the rows and columns at the boundary are done at the end of this function */
			/* The value of iCol and iRow indicate the columns and rows for which the derivates have already been ccomputed */
			for (Int iCol = 1; iCol < (iWidth - 1); iCol += 2)
			{
				mmPred[0] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[0 * iPredStride + iCol - 1]));
				mmPred[1] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[1 * iPredStride + iCol - 1]));

				mmPred[0] = _mm_cvtepi16_epi32(mmPred[0]);
				mmPred[1] = _mm_cvtepi16_epi32(mmPred[1]);

				for (Int iRow = 1; iRow < (iHeight - 1); iRow += 2)
				{
					mmPred[2] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[(iRow + 1) * iPredStride + iCol - 1]));
					mmPred[3] = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(&pPred[(iRow + 2) * iPredStride + iCol - 1]));

					mmPred[2] = _mm_cvtepi16_epi32(mmPred[2]);
					mmPred[3] = _mm_cvtepi16_epi32(mmPred[3]);

					mmIntermediates[0] = _mm_sub_epi32(mmPred[2], mmPred[0]);
					mmIntermediates[3] = _mm_sub_epi32(mmPred[3], mmPred[1]);

					mmPred[0] = mmPred[2];
					mmPred[1] = mmPred[3];

					mmIntermediates[1] = _mm_srli_si128(mmIntermediates[0], 4);
					mmIntermediates[4] = _mm_srli_si128(mmIntermediates[3], 4);
					mmIntermediates[2] = _mm_srli_si128(mmIntermediates[0], 8);
					mmIntermediates[5] = _mm_srli_si128(mmIntermediates[3], 8);

					mmIntermediates[1] = _mm_slli_epi32(mmIntermediates[1], 1);
					mmIntermediates[4] = _mm_slli_epi32(mmIntermediates[4], 1);

					mmIntermediates[0] = _mm_add_epi32(mmIntermediates[0], mmIntermediates[2]);
					mmIntermediates[3] = _mm_add_epi32(mmIntermediates[3], mmIntermediates[5]);

					mmDerivate[0] = _mm_add_epi32(mmIntermediates[0], mmIntermediates[1]);
					mmDerivate[1] = _mm_add_epi32(mmIntermediates[3], mmIntermediates[4]);

					_mm_storel_epi64(reinterpret_cast<__m128i *> (&piDerivate[iCol + (iRow + 0) * iDerivateBufStride]), mmDerivate[0]);
					_mm_storel_epi64(reinterpret_cast<__m128i *> (&piDerivate[iCol + (iRow + 1) * iDerivateBufStride]), mmDerivate[1]);
				}
			}

			for (Int j = 1; j < (iHeight - 1); j++)
			{
				piDerivate[j * iDerivateBufStride] = piDerivate[j * iDerivateBufStride + 1];
				piDerivate[j * iDerivateBufStride + (iWidth - 1)] = piDerivate[j * iDerivateBufStride + (iWidth - 2)];
			}

			memcpy
				(
					piDerivate,
					piDerivate + iDerivateBufStride,
					iWidth * sizeof(piDerivate[0])
					);

			memcpy
				(
					piDerivate + (iHeight - 1) * iDerivateBufStride,
					piDerivate + (iHeight - 2) * iDerivateBufStride,
					iWidth * sizeof(piDerivate[0])
					);
		}
#endif

#if ITH_C3_EQUAL_COEFF_COMPUTER==C_FXN
		Void TEncSearch::equalCoeffComputer
			(
				Pel *pResidue,
				Int iResidueStride,
				Int **ppiDerivate,
				Int iDerivateBufStride,
				Int64(*pi64EqualCoeff)[TEncSearch::iEqualCoeffBufDimension],
				Int iWidth,
				Int iHeight,
				Bool b6Param
				)
		{
			const Int iAffineParamNum = b6Param ? 6 : 4;

			for (Int j = 0; j != iHeight; j++)
			{
				for (Int k = 0; k != iWidth; k++)
				{
					Int64 i64Intermediates[2];
					Int iC[6];

					Int iIdx = j * iDerivateBufStride + k;

					if (!b6Param)
					{
						iC[0] = ppiDerivate[0][iIdx];
						iC[1] = k * ppiDerivate[0][iIdx];
						iC[1] += j * ppiDerivate[1][iIdx];
						iC[2] = ppiDerivate[1][iIdx];
						iC[3] = j * ppiDerivate[0][iIdx];
						iC[3] -= k * ppiDerivate[1][iIdx];
					}
					else
					{
						iC[0] = ppiDerivate[0][iIdx];
						iC[1] = k * ppiDerivate[0][iIdx];
						iC[2] = ppiDerivate[1][iIdx];
						iC[3] = k * ppiDerivate[1][iIdx];
						iC[4] = j * ppiDerivate[0][iIdx];
						iC[5] = j * ppiDerivate[1][iIdx];
					}

					for (Int col = 0; col < iAffineParamNum; col++)
					{
						i64Intermediates[0] = iC[col];

						for (Int row = 0; row < iAffineParamNum; row++)
						{
							i64Intermediates[1] = i64Intermediates[0] * iC[row];

							pi64EqualCoeff[col + 1][row] += i64Intermediates[1];
						}

						i64Intermediates[1] = i64Intermediates[0] * pResidue[iIdx];

						pi64EqualCoeff[col + 1][iAffineParamNum] += i64Intermediates[1] * 8;
					}
				}
			}
		}
#else
#define CALC_EQUAL_COEFF_8PXLS(x1,x2,y1,y2,tmp0,tmp1,tmp2,tmp3,inter0,inter1,inter2,inter3,loadLocation)       \
{                                                                                                              \
inter0 = _mm_mul_epi32(x1, y1);                                                                                \
inter1 = _mm_mul_epi32(tmp0, tmp2);                                                                            \
inter2 = _mm_mul_epi32(x2, y2);                                                                                \
inter3 = _mm_mul_epi32(tmp1, tmp3);                                                                            \
inter2 = _mm_add_epi64(inter0, inter2);                                                                        \
inter3 = _mm_add_epi64(inter1, inter3);                                                                        \
inter0 = _mm_loadl_epi64(loadLocation);                                                                        \
inter3 = _mm_add_epi64(inter2, inter3);                                                                        \
inter1 = _mm_srli_si128(inter3, 8);                                                                            \
inter3 = _mm_add_epi64(inter1, inter3);                                                                        \
inter3 = _mm_add_epi64(inter0, inter3);                                                                        \
}

		Void TEncSearch::equalCoeffComputer
			(
				Pel *pResidue,
				Int iResidueStride,
				Int **ppiDerivate,
				Int iDerivateBufStride,
				Int64(*pi64EqualCoeff)[TEncSearch::iEqualCoeffBufDimension],
				Int iWidth,
				Int iHeight,
				Bool b6Param
				)
		{
			__m128i mmTwo, mmFour;
			__m128i mmTmp[4];
			__m128i mmIntermediate[4];
			__m128i mmIndxK, mmIndxJ[2];
			__m128i mmResidue[2];

			// Add directly to indexes to get new index
			mmTwo = _mm_set1_epi32(2);
			mmFour = _mm_set1_epi32(4);

			if (b6Param)
			{
				__m128i mmC[12];

				//  mmC ¨C map 
				//  C for 1st row of pixels
				//  mmC[0] = iC[0][i] | iC[0][i+1] | iC[0][i+2] | iC[0][i+3]
				//  mmC[1] = iC[1][i] | iC[1][i+1] | iC[1][i+2] | iC[1][i+3]
				//  mmC[2] = iC[2][i] | iC[2][i+1] | iC[2][i+2] | iC[2][i+3]
				//  mmC[3] = iC[3][i] | iC[3][i+1] | iC[3][i+2] | iC[3][i+3]
				//  mmC[4] = iC[4][i] | iC[4][i+1] | iC[4][i+2] | iC[4][i+3]
				//  mmC[5] = iC[5][i] | iC[5][i+1] | iC[5][i+2] | iC[5][i+3]

				//  C for 2nd row of pixels
				//  mmC[6] = iC[6][i] | iC[6][i+1] | iC[6][i+2] | iC[6][i+3]
				//  mmC[7] = iC[7][i] | iC[7][i+1] | iC[7][i+2] | iC[7][i+3]
				//  mmC[8] = iC[8][i] | iC[8][i+1] | iC[8][i+2] | iC[8][i+3]
				//  mmC[9] = iC[9][i] | iC[9][i+1] | iC[9][i+2] | iC[9][i+3]
				//  mmC[10] = iC[10][i] | iC[10][i+1] | iC[10][i+2] | iC[10][i+3]
				//  mmC[11] = iC[11][i] | iC[11][i+1] | iC[11][i+2] | iC[11][i+3]

				Int iIdx1 = 0, iIdx2 = 0;
				iIdx1 = -2 * iDerivateBufStride - 4;
				iIdx2 = -iDerivateBufStride - 4;
				mmIndxJ[0] = _mm_set1_epi32(-2);
				mmIndxJ[1] = _mm_set1_epi32(-1);

				for (Int j = 0; j < iHeight; j += 2)
				{
					mmIndxJ[0] = _mm_add_epi32(mmIndxJ[0], mmTwo);
					mmIndxJ[1] = _mm_add_epi32(mmIndxJ[1], mmTwo);
					mmIndxK = _mm_set_epi32(-1, -2, -3, -4);
					iIdx1 += (iDerivateBufStride << 1);
					iIdx2 += (iDerivateBufStride << 1);

					for (Int k = 0; k < iWidth; k += 4)
					{
						iIdx1 += 4;
						iIdx2 += 4;

						mmIndxK = _mm_add_epi32(mmIndxK, mmFour);

						// 1st row
						mmC[0] = _mm_loadu_si128((const __m128i*)&ppiDerivate[0][iIdx1]);
						mmC[2] = _mm_loadu_si128((const __m128i*)&ppiDerivate[1][iIdx1]);
						// 2nd row
						mmC[6] = _mm_loadu_si128((const __m128i*)&ppiDerivate[0][iIdx2]);
						mmC[8] = _mm_loadu_si128((const __m128i*)&ppiDerivate[1][iIdx2]);

						// 1st row
						mmC[1] = _mm_mullo_epi32(mmIndxK, mmC[0]);
						mmC[3] = _mm_mullo_epi32(mmIndxK, mmC[2]);
						mmC[4] = _mm_mullo_epi32(mmIndxJ[0], mmC[0]);
						mmC[5] = _mm_mullo_epi32(mmIndxJ[0], mmC[2]);

						// 2nd row
						mmC[7] = _mm_mullo_epi32(mmIndxK, mmC[6]);
						mmC[9] = _mm_mullo_epi32(mmIndxK, mmC[8]);
						mmC[10] = _mm_mullo_epi32(mmIndxJ[1], mmC[6]);
						mmC[11] = _mm_mullo_epi32(mmIndxJ[1], mmC[8]);

						// Residue
						mmResidue[0] = _mm_loadl_epi64((const __m128i*)&pResidue[iIdx1]);
						mmResidue[1] = _mm_loadl_epi64((const __m128i*)&pResidue[iIdx2]);

						mmResidue[0] = _mm_cvtepi16_epi32(mmResidue[0]);
						mmResidue[1] = _mm_cvtepi16_epi32(mmResidue[1]);

						mmResidue[0] = _mm_slli_epi32(mmResidue[0], 3);
						mmResidue[1] = _mm_slli_epi32(mmResidue[1], 3);

						// Calculate residue coefficients first
						mmTmp[2] = _mm_srli_si128(mmResidue[0], 4);
						mmTmp[3] = _mm_srli_si128(mmResidue[1], 4);

						// 1st row
						mmTmp[0] = _mm_srli_si128(mmC[0], 4);
						mmTmp[1] = _mm_srli_si128(mmC[6], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][6], mmIntermediate[3]);

						// 2nd row
						mmTmp[0] = _mm_srli_si128(mmC[1], 4);
						mmTmp[1] = _mm_srli_si128(mmC[7], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][6], mmIntermediate[3]);

						// 3rd row
						mmTmp[0] = _mm_srli_si128(mmC[2], 4);
						mmTmp[1] = _mm_srli_si128(mmC[8], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[8], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][6], mmIntermediate[3]);

						// 4th row
						mmTmp[0] = _mm_srli_si128(mmC[3], 4);
						mmTmp[1] = _mm_srli_si128(mmC[9], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[9], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][6], mmIntermediate[3]);

						// 5th row
						mmTmp[0] = _mm_srli_si128(mmC[4], 4);
						mmTmp[1] = _mm_srli_si128(mmC[10], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[4], mmC[10], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[5][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][6], mmIntermediate[3]);

						// 6th row
						mmTmp[0] = _mm_srli_si128(mmC[5], 4);
						mmTmp[1] = _mm_srli_si128(mmC[11], 4);
						// 7th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[5], mmC[11], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[6][6]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][6], mmIntermediate[3]);

						//Start calculation of coefficient matrix
						// 1st row
						mmTmp[0] = _mm_srli_si128(mmC[0], 4);
						mmTmp[1] = _mm_srli_si128(mmC[6], 4);

						// 1st col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[0], mmC[6], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][0]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][0], mmIntermediate[3]);
						// 2nd col of row and 1st col of 2nd row
						mmTmp[2] = _mm_srli_si128(mmC[1], 4);
						mmTmp[3] = _mm_srli_si128(mmC[7], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[1], mmC[7], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][1]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][1], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][0], mmIntermediate[3]);
						// 3rd col of row and 1st col of 3rd row
						mmTmp[2] = _mm_srli_si128(mmC[2], 4);
						mmTmp[3] = _mm_srli_si128(mmC[8], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[2], mmC[8], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][2], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][0], mmIntermediate[3]);
						// 4th col of row and 1st col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[9], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[3], mmC[9], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][0], mmIntermediate[3]);
						// 5th col of row and 1st col of the 5th row
						mmTmp[2] = _mm_srli_si128(mmC[4], 4);
						mmTmp[3] = _mm_srli_si128(mmC[10], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[4], mmC[10], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][4], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][0], mmIntermediate[3]);
						// 6th col of row and 1st col of the 6th row
						mmTmp[2] = _mm_srli_si128(mmC[5], 4);
						mmTmp[3] = _mm_srli_si128(mmC[11], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[6], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][5], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][0], mmIntermediate[3]);

						// 2nd row
						mmTmp[0] = _mm_srli_si128(mmC[1], 4);
						mmTmp[1] = _mm_srli_si128(mmC[7], 4);

						// 2nd col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmC[1], mmC[7], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][1]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][1], mmIntermediate[3]);
						// 3rd col of row and 2nd col of 3rd row
						mmTmp[2] = _mm_srli_si128(mmC[2], 4);
						mmTmp[3] = _mm_srli_si128(mmC[8], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmC[2], mmC[8], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][2], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][1], mmIntermediate[3]);
						// 4th col of row and 2nd col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[9], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmC[3], mmC[9], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][1], mmIntermediate[3]);
						// 5th col of row and 1st col of the 5th row
						mmTmp[2] = _mm_srli_si128(mmC[4], 4);
						mmTmp[3] = _mm_srli_si128(mmC[10], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmC[4], mmC[10], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][4], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][1], mmIntermediate[3]);
						// 6th col of row and 1st col of the 6th row
						mmTmp[2] = _mm_srli_si128(mmC[5], 4);
						mmTmp[3] = _mm_srli_si128(mmC[11], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[7], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][5], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][1], mmIntermediate[3]);

						// 3rd row
						mmTmp[0] = _mm_srli_si128(mmC[2], 4);
						mmTmp[1] = _mm_srli_si128(mmC[8], 4);

						//3rd Col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[8], mmC[2], mmC[8], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][2], mmIntermediate[3]);
						// 4th col of row and 3rd col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[9], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[8], mmC[3], mmC[9], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][2], mmIntermediate[3]);
						// 5th col of row and 1st col of the 5th row
						mmTmp[2] = _mm_srli_si128(mmC[4], 4);
						mmTmp[3] = _mm_srli_si128(mmC[10], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[8], mmC[4], mmC[10], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][4], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][2], mmIntermediate[3]);
						// 6th col of row and 1st col of the 6th row
						mmTmp[2] = _mm_srli_si128(mmC[5], 4);
						mmTmp[3] = _mm_srli_si128(mmC[11], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[8], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][5], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][2], mmIntermediate[3]);

						// 4th row
						mmTmp[0] = _mm_srli_si128(mmC[3], 4);
						mmTmp[1] = _mm_srli_si128(mmC[9], 4);

						// 4th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[9], mmC[3], mmC[9], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][3], mmIntermediate[3]);
						// 5th col of row and 1st col of the 5th row
						mmTmp[2] = _mm_srli_si128(mmC[4], 4);
						mmTmp[3] = _mm_srli_si128(mmC[10], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[9], mmC[4], mmC[10], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][4], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][3], mmIntermediate[3]);
						// 6th col of row and 1st col of the 6th row
						mmTmp[2] = _mm_srli_si128(mmC[5], 4);
						mmTmp[3] = _mm_srli_si128(mmC[11], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[9], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][5], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][3], mmIntermediate[3]);

						// 5th row
						mmTmp[0] = _mm_srli_si128(mmC[4], 4);
						mmTmp[1] = _mm_srli_si128(mmC[10], 4);
						// 5th col of row and 1st col of the 5th row
						CALC_EQUAL_COEFF_8PXLS(mmC[4], mmC[10], mmC[4], mmC[10], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[5][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][4], mmIntermediate[3]);
						// 6th col of row and 1st col of the 6th row
						mmTmp[2] = _mm_srli_si128(mmC[5], 4);
						mmTmp[3] = _mm_srli_si128(mmC[11], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[4], mmC[10], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[5][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[5][5], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][4], mmIntermediate[3]);

						// 6th row
						mmTmp[0] = _mm_srli_si128(mmC[5], 4);
						mmTmp[1] = _mm_srli_si128(mmC[11], 4);
						// 5th col of row and 1st col of the 5th row
						CALC_EQUAL_COEFF_8PXLS(mmC[5], mmC[11], mmC[5], mmC[11], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[6][5]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[6][5], mmIntermediate[3]);
					}

					iIdx1 -= (iWidth);
					iIdx2 -= (iWidth);
				}
			}
			else
			{
				__m128i mmC[8];

				//  mmC ? map 
				//  C for 1st row of pixels
				//  mmC[0] = iC[0][i] | iC[0][i+1] | iC[0][i+2] | iC[0][i+3]
				//  mmC[1] = iC[1][i] | iC[1][i+1] | iC[1][i+2] | iC[1][i+3]
				//  mmC[2] = iC[2][i] | iC[2][i+1] | iC[2][i+2] | iC[2][i+3]
				//  mmC[3] = iC[3][i] | iC[3][i+1] | iC[3][i+2] | iC[3][i+3]
				//  C for 2nd row of pixels
				//  mmC[4] = iC[0][i] | iC[0][i+1] | iC[0][i+2] | iC[0][i+3]
				//  mmC[5] = iC[1][i] | iC[1][i+1] | iC[1][i+2] | iC[1][i+3]
				//  mmC[6] = iC[2][i] | iC[2][i+1] | iC[2][i+2] | iC[2][i+3]
				//  mmC[7] = iC[3][i] | iC[3][i+1] | iC[3][i+2] | iC[3][i+3]
				Int iIdx1 = 0, iIdx2 = 0;
				iIdx1 = -2 * iDerivateBufStride - 4;
				iIdx2 = -iDerivateBufStride - 4;
				mmIndxJ[0] = _mm_set1_epi32(-2);
				mmIndxJ[1] = _mm_set1_epi32(-1);

				for (Int j = 0; j < iHeight; j += 2)
				{
					mmIndxJ[0] = _mm_add_epi32(mmIndxJ[0], mmTwo);
					mmIndxJ[1] = _mm_add_epi32(mmIndxJ[1], mmTwo);
					mmIndxK = _mm_set_epi32(-1, -2, -3, -4);
					iIdx1 += (iDerivateBufStride << 1);
					iIdx2 += (iDerivateBufStride << 1);

					for (Int k = 0; k < iWidth; k += 4)
					{
						iIdx1 += 4;
						iIdx2 += 4;

						mmIndxK = _mm_add_epi32(mmIndxK, mmFour);

						mmC[0] = _mm_loadu_si128((const __m128i*)&ppiDerivate[0][iIdx1]);
						mmC[2] = _mm_loadu_si128((const __m128i*)&ppiDerivate[1][iIdx1]);
						mmC[4] = _mm_loadu_si128((const __m128i*)&ppiDerivate[0][iIdx2]);
						mmC[6] = _mm_loadu_si128((const __m128i*)&ppiDerivate[1][iIdx2]);

						mmC[1] = _mm_mullo_epi32(mmIndxK, mmC[0]);
						mmC[3] = _mm_mullo_epi32(mmIndxJ[0], mmC[0]);
						mmC[5] = _mm_mullo_epi32(mmIndxK, mmC[4]);
						mmC[7] = _mm_mullo_epi32(mmIndxJ[1], mmC[4]);

						mmResidue[0] = _mm_loadl_epi64((const __m128i*)&pResidue[iIdx1]);
						mmResidue[1] = _mm_loadl_epi64((const __m128i*)&pResidue[iIdx2]);

						mmTmp[0] = _mm_mullo_epi32(mmIndxJ[0], mmC[2]);
						mmTmp[1] = _mm_mullo_epi32(mmIndxK, mmC[2]);
						mmTmp[2] = _mm_mullo_epi32(mmIndxJ[1], mmC[6]);
						mmTmp[3] = _mm_mullo_epi32(mmIndxK, mmC[6]);

						mmResidue[0] = _mm_cvtepi16_epi32(mmResidue[0]);
						mmResidue[1] = _mm_cvtepi16_epi32(mmResidue[1]);

						mmC[1] = _mm_add_epi32(mmC[1], mmTmp[0]);
						mmC[3] = _mm_sub_epi32(mmC[3], mmTmp[1]);
						mmC[5] = _mm_add_epi32(mmC[5], mmTmp[2]);
						mmC[7] = _mm_sub_epi32(mmC[7], mmTmp[3]);

						mmResidue[0] = _mm_slli_epi32(mmResidue[0], 3);
						mmResidue[1] = _mm_slli_epi32(mmResidue[1], 3);

						//Start calculation of coefficient matrix
						// 1st row
						mmTmp[0] = _mm_srli_si128(mmC[0], 4);
						mmTmp[1] = _mm_srli_si128(mmC[4], 4);

						// 1st col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[4], mmC[0], mmC[4], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][0]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][0], mmIntermediate[3]);
						// 2nd col of row and 1st col of 2nd row
						mmTmp[2] = _mm_srli_si128(mmC[1], 4);
						mmTmp[3] = _mm_srli_si128(mmC[5], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[4], mmC[1], mmC[5], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][1]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][1], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][0], mmIntermediate[3]);
						// 3rd col of row and 1st col of 3rd row
						mmTmp[2] = _mm_srli_si128(mmC[2], 4);
						mmTmp[3] = _mm_srli_si128(mmC[6], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[4], mmC[2], mmC[6], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][2], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][0], mmIntermediate[3]);
						// 4th col of row and 1st col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[7], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[4], mmC[3], mmC[7], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][0], mmIntermediate[3]);
						// 5th col of row
						mmTmp[2] = _mm_srli_si128(mmResidue[0], 4);
						mmTmp[3] = _mm_srli_si128(mmResidue[1], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[0], mmC[4], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[1][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[1][4], mmIntermediate[3]);

						// 2nd row
						mmTmp[0] = _mm_srli_si128(mmC[1], 4);
						mmTmp[1] = _mm_srli_si128(mmC[5], 4);

						// 2nd col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[5], mmC[1], mmC[5], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][1]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][1], mmIntermediate[3]);
						// 3rd col of row and 2nd col of 3rd row
						mmTmp[2] = _mm_srli_si128(mmC[2], 4);
						mmTmp[3] = _mm_srli_si128(mmC[6], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[5], mmC[2], mmC[6], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][2], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][1], mmIntermediate[3]);
						// 4th col of row and 2nd col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[7], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[5], mmC[3], mmC[7], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][1], mmIntermediate[3]);
						// 5th col of row
						mmTmp[2] = _mm_srli_si128(mmResidue[0], 4);
						mmTmp[3] = _mm_srli_si128(mmResidue[1], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[1], mmC[5], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[2][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[2][4], mmIntermediate[3]);

						// 3rd row
						mmTmp[0] = _mm_srli_si128(mmC[2], 4);
						mmTmp[1] = _mm_srli_si128(mmC[6], 4);

						//3rd Col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[6], mmC[2], mmC[6], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][2]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][2], mmIntermediate[3]);
						// 4th col of row and 3rd col of 4th row
						mmTmp[2] = _mm_srli_si128(mmC[3], 4);
						mmTmp[3] = _mm_srli_si128(mmC[7], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[6], mmC[3], mmC[7], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][3], mmIntermediate[3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][2], mmIntermediate[3]);
						// 5th col of row
						mmTmp[2] = _mm_srli_si128(mmResidue[0], 4);
						mmTmp[3] = _mm_srli_si128(mmResidue[1], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[2], mmC[6], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[3][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[3][4], mmIntermediate[3]);

						// 4th row
						mmTmp[0] = _mm_srli_si128(mmC[3], 4);
						mmTmp[1] = _mm_srli_si128(mmC[7], 4);

						// 4th col of row
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[7], mmC[3], mmC[7], mmTmp[0], mmTmp[1], mmTmp[0], mmTmp[1], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][3]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][3], mmIntermediate[3]);
						// 5th col of row
						mmTmp[2] = _mm_srli_si128(mmResidue[0], 4);
						mmTmp[3] = _mm_srli_si128(mmResidue[1], 4);
						CALC_EQUAL_COEFF_8PXLS(mmC[3], mmC[7], mmResidue[0], mmResidue[1], mmTmp[0], mmTmp[1], mmTmp[2], mmTmp[3], mmIntermediate[0], mmIntermediate[1], mmIntermediate[2], mmIntermediate[3], (const __m128i*)&pi64EqualCoeff[4][4]);
						_mm_storel_epi64((__m128i*)&pi64EqualCoeff[4][4], mmIntermediate[3]);
					}

					iIdx1 -= (iWidth);
					iIdx2 -= (iWidth);
				}
			}
		}
#endif

		Void TEncSearch::xAffineMotionEstimation(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, RefPicList eRefPicList, TComMv acMvPred[3], Int iRefIdxPred, TComMv acMv[3], UInt& ruiBits, Distortion& ruiCost,
#if JVECT_C0062_AFFINE_SIX_PARAM
			Bool b6Param,
#endif  
			Bool bBi)
		{
			// Get Block information
			UInt          uiPartAddr;
			Int           iRoiWidth;
			Int           iRoiHeight;
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);

			// Get RefPic
			TComPicYuv* pcPicYuvRef = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred)->getPicYuvRec();

			// Set Origin YUV: pcYuv
			TComYuv*      pcYuv = pcYuvOrg;
			Double        fWeight = 1.0;

			// if Bi, set to ( 2 * Org - ListX )
			if (bBi)
			{
				TComYuv*  pcYuvOther = &m_acYuvPred[1 - (Int)eRefPicList];
				pcYuv = &m_cYuvPredTemp;

				pcYuvOrg->copyPartToPartYuv(pcYuv, uiPartAddr, iRoiWidth, iRoiHeight);
				pcYuv->removeHighFreq(pcYuvOther, uiPartAddr, iRoiWidth, iRoiHeight, pcCU->getSlice()->getSPS()->getBitDepths().recon, m_pcEncCfg->getClipForBiPredMeEnabled());

				fWeight = 0.5;
			}

			// pred YUV
			TComYuv* pPredYuv = &m_tmpYuvPred;

			Int iOrgStride = pcYuv->getStride(COMPONENT_Y);
			Int iPredStride = pPredYuv->getStride(COMPONENT_Y);

			// Set start Mv position, use input mv as started search mv
			TComMv acMvTemp[3];
			memcpy(acMvTemp, acMv, sizeof(TComMv) * 3);

			// Set delta mv
			TComMv acDeltaMv[3];
			for (Int i = 0; i < 3; i++)
				acDeltaMv[i].set(0, 0);

#if JVECT_C0062_AFFINE_SIX_PARAM
			Double dAffinePara[6];
			Double dDeltaMv[6];
#else
			Double dAffinePara[4];
			Double dDeltaMv[4];
#endif

			// malloc buffer
#if JVECT_C0062_AFFINE_SIX_PARAM
			Int iParaNum = b6Param ? 7 : 5;
			Int iAffineParamNum = iParaNum - 1;
#else
			Int iParaNum = 5;
#endif
			Int64 i64EqualCoeff[TEncSearch::iEqualCoeffBufDimension][TEncSearch::iEqualCoeffBufDimension];
			Double dEqualCoeff[TEncSearch::iEqualCoeffBufDimension][TEncSearch::iEqualCoeffBufDimension];

			Pel    *pError = m_tmpError;
			Int    *piDerivate[2];
			piDerivate[0] = m_tmpDerivate[0];
			piDerivate[1] = m_tmpDerivate[1];

			Distortion uiCostTemp = std::numeric_limits<Distortion>::max();;
			Distortion uiCostBest = std::numeric_limits<Distortion>::max();;
			UInt uiBitsTemp = 0;
			UInt uiBitsBest = 0;

			// do motion compensation with origin mv
#if JVECT_C0062_AFFINE_SIX_PARAM
			Bool bRecShape = (iRoiHeight / iRoiWidth >= 2);
			if (b6Param == 0)
			{
				if (bRecShape)
				{
					pcCU->clipMv(acMvTemp[0]);
					pcCU->clipMv(acMvTemp[2]);
					Int vx1 = (acMvTemp[2].getVer() - acMvTemp[0].getVer())*iRoiWidth / iRoiHeight + acMvTemp[0].getHor();
					Int vy1 = -(acMvTemp[2].getHor() - acMvTemp[0].getHor())*iRoiWidth / iRoiHeight + acMvTemp[0].getVer();
					acMvTemp[1].set(vx1, vy1);
					pcCU->clipMv(acMvTemp[1]);
				}
				else
				{
					pcCU->clipMv(acMvTemp[0]);
					pcCU->clipMv(acMvTemp[1]);
					Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
					Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
					acMvTemp[2].set(vx2, vy2);
					pcCU->clipMv(acMvTemp[2]);
				}
			}
			else
			{
				pcCU->clipMv(acMvTemp[0]);
				pcCU->clipMv(acMvTemp[1]);
				pcCU->clipMv(acMvTemp[2]);
			}
#else
			pcCU->clipMv(acMvTemp[0]);
			pcCU->clipMv(acMvTemp[1]);
			Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
			Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
			acMvTemp[2].set(vx2, vy2);
			pcCU->clipMv(acMvTemp[2]);
#endif

			xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvTemp, iRoiWidth, iRoiHeight, pPredYuv, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
				, b6Param
#endif
				);

			// get error
			uiCostBest = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pPredYuv->getAddr(COMPONENT_Y, uiPartAddr), pPredYuv->getStride(COMPONENT_Y), pcYuv->getAddr(COMPONENT_Y, uiPartAddr), pcYuv->getStride(COMPONENT_Y), iRoiWidth, iRoiHeight, COMPONENT_Y, DF_HADS);

			// get cost with mv
			m_pcRdCost->setCostScale(0);
			uiBitsBest = ruiBits;
#if JVECT_C0062_AFFINE_SIX_PARAM
			Int numVerIdx = (b6Param == 1) ? 3 : 2;
			for (Int i = 0; i < numVerIdx; i++)
#else
			for (Int i = 0; i<2; i++)
#endif
			{
#if JVECT_C0062_AFFINE_SIX_PARAM
				Int iTmpIdx = (bRecShape && i == 1 && b6Param == 0) ? (i + 1) : i;
				m_pcRdCost->setPredictor(acMvPred[iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
				TComPattern   tmpPattern;
				TComPattern*  pcPatternKey = &tmpPattern;
				pcPatternKey->m_pcCU = pcCU;
				pcPatternKey->m_uiPartAddr = uiPartAddr;
				pcPatternKey->m_mvPred = &acMvPred[iTmpIdx];
				iCostScale = 0;
				uiBitsBest += xMvdBits(pcPatternKey, acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer());
#else
				uiBitsBest += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMvTemp[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
					);
#endif
#else
				m_pcRdCost->setPredictor(acMvPred[i]);
#if JVET_D0123_ME_CTX_LUT_BITS
				TComPattern   tmpPattern;
				TComPattern*  pcPatternKey = &tmpPattern;
				pcPatternKey->m_pcCU = pcCU;
				pcPatternKey->m_uiPartAddr = uiPartAddr;
				pcPatternKey->m_mvPred = &acMvPred[i];
				iCostScale = 0;
				uiBitsBest += xMvdBits(pcPatternKey, acMvTemp[i].getHor(), acMvTemp[i].getVer());
#else
				uiBitsBest += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMvTemp[i].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[i].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMvTemp[i].getHor(), acMvTemp[i].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
					);
#endif
#endif
			}
#if JVET_D0123_ME_CTX_LUT_BITS
			uiCostBest = (UInt)(floor(fWeight * (Double)uiCostBest) + (Double)m_pcRdCost->getCostSearch(uiBitsBest));
#else
			uiCostBest = (UInt)(floor(fWeight * (Double)uiCostBest) + (Double)m_pcRdCost->getCost(uiBitsBest));
#endif

#if COM16_C1016_AFFINE
			memcpy(acMv, acMvTemp, sizeof(TComMv) * 3);
#endif

			Int iIterTime = bBi ? 5 : 7;
#if AFFINE_REDUCE_ITERATION
			if (b6Param)
			{
				iIterTime = bBi ? 3 : 5;
			}
#endif
			for (Int iter = 0; iter < iIterTime; iter++)    // iterate loop
			{
				/*********************************************************************************
				*                         use gradient to update mv
				*********************************************************************************/
				// get Error Matrix
				Pel* pOrg = pcYuv->getAddr(COMPONENT_Y, uiPartAddr);
				Pel* pPred = pPredYuv->getAddr(COMPONENT_Y, uiPartAddr);
				residueComputer
					(
						pOrg,
						iOrgStride,
						pPred,
						iPredStride,
						pError,
						iRoiWidth,
						iRoiWidth,
						iRoiHeight
						);

				// sobel x direction
				// -1 0 1
				// -2 0 2
				// -1 0 1
				pPred = pPredYuv->getAddr(COMPONENT_Y, uiPartAddr);
				scaledHorizontalSobelFilter
					(
						pPred,
						iPredStride,
						piDerivate[0],
						iRoiWidth,
						iRoiWidth,
						iRoiHeight
						);

				// sobel y direction
				// -1 -2 -1
				//  0  0  0
				//  1  2  1
				scaledVerticalSobelFilter
					(
						pPred,
						iPredStride,
						piDerivate[1],
						iRoiWidth,
						iRoiWidth,
						iRoiHeight
						);

				// solve delta x and y
				for (Int row = 0; row < iParaNum; row++)
				{
					memset
						(
							&i64EqualCoeff[row][0],
							0,
							iParaNum * sizeof(Int64)
							);
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				equalCoeffComputer
					(
						pError,
						iRoiWidth,
						piDerivate,
						iRoiWidth,
						i64EqualCoeff,
						iRoiWidth,
						iRoiHeight,
						b6Param
						);
#else
				equalCoeffComputer
					(
						pError,
						iRoiWidth,
						piDerivate,
						iRoiWidth,
						i64EqualCoeff,
						iRoiWidth,
						iRoiHeight,
						false
						);
#endif

				for (Int row = 0; row < iParaNum; row++)
				{
					LongLong2DoubleConverter
						(
							&i64EqualCoeff[row][0],
							&dEqualCoeff[row][0],
							iParaNum
							);
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				solveEqual(dEqualCoeff, iAffineParamNum, dAffinePara);
#else
				solveEqual(dEqualCoeff, 4, dAffinePara);
#endif

				// convert to delta mv
#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param)
				{
					//for MV0
					dDeltaMv[0] = dAffinePara[0];
					dDeltaMv[2] = dAffinePara[2];
					//for MV1
					dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];//ZF
					dDeltaMv[3] = dAffinePara[3] * iRoiWidth + dAffinePara[2];
					//for MV2
					dDeltaMv[4] = dAffinePara[4] * iRoiHeight + dAffinePara[0];//ZF
					dDeltaMv[5] = dAffinePara[5] * iRoiHeight + dAffinePara[2];
				}
				else
				{
					if (bRecShape)
					{
						dDeltaMv[0] = dAffinePara[0];
						dDeltaMv[2] = dAffinePara[2];

						dDeltaMv[1] = dAffinePara[3] * iRoiHeight + dAffinePara[0];
						dDeltaMv[3] = dAffinePara[1] * iRoiHeight + dAffinePara[2];
					}
					else
					{
						dDeltaMv[0] = dAffinePara[0];
						dDeltaMv[2] = dAffinePara[2];

						dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];//ZF
						dDeltaMv[3] = -dAffinePara[3] * iRoiWidth + dAffinePara[2];
					}
				}
#else
				dDeltaMv[0] = dAffinePara[0];
				dDeltaMv[2] = dAffinePara[2];

				dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];
				dDeltaMv[3] = -dAffinePara[3] * iRoiWidth + dAffinePara[2];
#endif

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param)
				{
					acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[2].setHor((Int)(dDeltaMv[4] * 4 + SIGN(dDeltaMv[4]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[2].setVer((Int)(dDeltaMv[5] * 4 + SIGN(dDeltaMv[5]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				}
				else
				{
					if (bRecShape)
					{
						acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[2].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[2].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					}
					else
					{
						acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					}
				}
#else
				acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#endif
#else
				acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5));
				acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5));
				acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5));
				acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5));
#endif

				Bool bAllZero = false;
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (bRecShape && i && b6Param == 0) ? (i + 1) : i;
					if (acDeltaMv[iTmpIdx].getHor() != 0 || acDeltaMv[iTmpIdx].getVer() != 0)
#else
					if (acDeltaMv[i].getHor() != 0 || acDeltaMv[i].getVer() != 0)
#endif
					{
						bAllZero = false;
						break;
					}
					bAllZero = true;
				}

				if (bAllZero)
					break;

				// do motion compensation with updated mv
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (bRecShape && i && b6Param == 0) ? (i + 1) : i;
					acMvTemp[iTmpIdx] += acDeltaMv[iTmpIdx];
					pcCU->clipMv(acMvTemp[iTmpIdx]);
#else
					acMvTemp[i] += acDeltaMv[i];
					pcCU->clipMv(acMvTemp[i]);
#endif
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param == 0)
				{
					if (bRecShape)
					{
						Int vx1 = (acMvTemp[2].getVer() - acMvTemp[0].getVer())*iRoiWidth / iRoiHeight + acMvTemp[0].getHor();
						Int vy1 = -(acMvTemp[2].getHor() - acMvTemp[0].getHor())*iRoiWidth / iRoiHeight + acMvTemp[0].getVer();
						acMvTemp[1].set(vx1, vy1);
						pcCU->clipMv(acMvTemp[1]);
					}
					else
					{
						Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
						Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
						acMvTemp[2].set(vx2, vy2);
						pcCU->clipMv(acMvTemp[2]);
					}
				}
#else
				vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
				vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
				acMvTemp[2].set(vx2, vy2);
				pcCU->clipMv(acMvTemp[2]);
#endif

				xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvTemp, iRoiWidth, iRoiHeight, pPredYuv, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
					, b6Param
#endif
					);

				// get error
				uiCostTemp = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pPredYuv->getAddr(COMPONENT_Y, uiPartAddr), pPredYuv->getStride(COMPONENT_Y), pcYuv->getAddr(COMPONENT_Y, uiPartAddr), pcYuv->getStride(COMPONENT_Y), iRoiWidth, iRoiHeight, COMPONENT_Y, DF_HADS);

				// get cost with mv
				m_pcRdCost->setCostScale(0);
				uiBitsTemp = ruiBits;
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (i && bRecShape && b6Param == 0) ? (i + 1) : i;
					m_pcRdCost->setPredictor(acMvPred[iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
					TComPattern   tmpPattern;
					TComPattern*  pcPatternKey = &tmpPattern;
					pcPatternKey->m_pcCU = pcCU;
					pcPatternKey->m_uiPartAddr = uiPartAddr;
					pcPatternKey->m_mvPred = &acMvPred[iTmpIdx];
					iCostScale = 0;
					uiBitsTemp += xMvdBits(pcPatternKey, acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer());
#else
					uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMvTemp[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
						);
#endif
#else
					m_pcRdCost->setPredictor(acMvPred[i]);
#if JVET_D0123_ME_CTX_LUT_BITS
					TComPattern   tmpPattern;
					TComPattern*  pcPatternKey = &tmpPattern;
					pcPatternKey->m_pcCU = pcCU;
					pcPatternKey->m_uiPartAddr = uiPartAddr;
					pcPatternKey->m_mvPred = &acMvPred[i];
					iCostScale = 0;
					uiBitsTemp += xMvdBits(pcPatternKey, acMvTemp[i].getHor(), acMvTemp[i].getVer());
#else 
					uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMvTemp[i].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[i].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMvTemp[i].getHor(), acMvTemp[i].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
						);
#endif
#endif
				}
#if JVET_D0123_ME_CTX_LUT_BITS
				uiCostTemp = (UInt)(floor(fWeight * (Double)uiCostTemp) + (Double)m_pcRdCost->getCostSearch(uiBitsTemp));
#else
				uiCostTemp = (UInt)(floor(fWeight * (Double)uiCostTemp) + (Double)m_pcRdCost->getCost(uiBitsTemp));
#endif

				// store best cost and mv
				if (uiCostTemp < uiCostBest)
				{
					uiCostBest = uiCostTemp;
					uiBitsBest = uiBitsTemp;
					memcpy(acMv, acMvTemp, sizeof(TComMv) * 3);
				}
			}

			ruiBits = uiBitsBest;
			ruiCost = uiCostBest;
		}
#else
		Void TEncSearch::xAffineMotionEstimation(TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, RefPicList eRefPicList, TComMv acMvPred[3], Int iRefIdxPred, TComMv acMv[3], UInt& ruiBits, Distortion& ruiCost,
#if JVECT_C0062_AFFINE_SIX_PARAM
			Bool b6Param,
#endif  
			Bool bBi)
		{
			// Get Block information
			UInt          uiPartAddr;
			Int           iRoiWidth;
			Int           iRoiHeight;
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight);

			// Get RefPic
			TComPicYuv* pcPicYuvRef = pcCU->getSlice()->getRefPic(eRefPicList, iRefIdxPred)->getPicYuvRec();

			// Set Origin YUV: pcYuv
			TComYuv*      pcYuv = pcYuvOrg;
			Double        fWeight = 1.0;

			// if Bi, set to ( 2 * Org - ListX )
			if (bBi)
			{
				TComYuv*  pcYuvOther = &m_acYuvPred[1 - (Int)eRefPicList];
				pcYuv = &m_cYuvPredTemp;

				pcYuvOrg->copyPartToPartYuv(pcYuv, uiPartAddr, iRoiWidth, iRoiHeight);
				pcYuv->removeHighFreq(pcYuvOther, uiPartAddr, iRoiWidth, iRoiHeight, pcCU->getSlice()->getSPS()->getBitDepths().recon, m_pcEncCfg->getClipForBiPredMeEnabled());

				fWeight = 0.5;
			}

			// pred YUV
			TComYuv* pPredYuv = &m_tmpYuvPred;

			Int iOrgStride = pcYuv->getStride(COMPONENT_Y);
			Int iPredStride = pPredYuv->getStride(COMPONENT_Y);

			// Set start Mv position, use input mv as started search mv
			TComMv acMvTemp[3];
			memcpy(acMvTemp, acMv, sizeof(TComMv) * 3);

			// Set delta mv
			TComMv acDeltaMv[3];
			for (Int i = 0; i<3; i++)
				acDeltaMv[i].set(0, 0);

#if JVECT_C0062_AFFINE_SIX_PARAM
			Double dAffinePara[6];
			Double dDeltaMv[6];
#else
			Double dAffinePara[4];
			Double dDeltaMv[4];
#endif

			// malloc buffer
#if JVECT_C0062_AFFINE_SIX_PARAM
			Int iParaNum = b6Param ? 7 : 5;
			Int iAffineParamNum = iParaNum - 1;
#else
			Int iParaNum = 5;
#endif
			Double **pdEqualCoeff;
			pdEqualCoeff = new Double *[iParaNum];
			for (Int i = 0; i < iParaNum; i++)
			{
				pdEqualCoeff[i] = new Double[iParaNum];
				for (Int j = 0; j < iParaNum; j++)
				{
					pdEqualCoeff[i][j] = 0.0;
				}
			}

			Int    *piError = m_tmpError;
			Double *pdDerivate[2];
			pdDerivate[0] = m_tmpDerivate[0];
			pdDerivate[1] = m_tmpDerivate[1];

			Distortion uiCostTemp = std::numeric_limits<Distortion>::max();;
			Distortion uiCostBest = std::numeric_limits<Distortion>::max();;
			UInt uiBitsTemp = 0;
			UInt uiBitsBest = 0;

			// do motion compensation with origin mv
#if JVECT_C0062_AFFINE_SIX_PARAM
			Bool bRecShape = (iRoiHeight / iRoiWidth >= 2);
			if (b6Param == 0)
			{
				if (bRecShape)
				{
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[0]);
					pcCU->clipMv(acMvTemp[2]);
#endif
					Int vx1 = (acMvTemp[2].getVer() - acMvTemp[0].getVer())*iRoiWidth / iRoiHeight + acMvTemp[0].getHor();
					Int vy1 = -(acMvTemp[2].getHor() - acMvTemp[0].getHor())*iRoiWidth / iRoiHeight + acMvTemp[0].getVer();
					acMvTemp[1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[1]);
#endif
				}
				else
				{
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[0]);
					pcCU->clipMv(acMvTemp[1]);
#endif
					Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
					Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
					acMvTemp[2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[2]);
#endif
				}
			}
#if !AFFINE_CLIP_REMOVE
			else
			{
				pcCU->clipMv(acMvTemp[0]);
				pcCU->clipMv(acMvTemp[1]);
				pcCU->clipMv(acMvTemp[2]);
			}
#endif
#else
#if !AFFINE_CLIP_REMOVE
			pcCU->clipMv(acMvTemp[0]);
			pcCU->clipMv(acMvTemp[1]);
#endif
			Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
			Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
			acMvTemp[2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
			pcCU->clipMv(acMvTemp[2]);
#endif
#endif

			xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvTemp, iRoiWidth, iRoiHeight, pPredYuv, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
				, b6Param
#endif
			);

			// get error
			uiCostBest = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pPredYuv->getAddr(COMPONENT_Y, uiPartAddr), pPredYuv->getStride(COMPONENT_Y), pcYuv->getAddr(COMPONENT_Y, uiPartAddr), pcYuv->getStride(COMPONENT_Y), iRoiWidth, iRoiHeight, COMPONENT_Y, DF_HADS);

			// get cost with mv
			m_pcRdCost->setCostScale(0);
			uiBitsBest = ruiBits;
#if JVECT_C0062_AFFINE_SIX_PARAM
			Int numVerIdx = (b6Param == 1) ? 3 : 2;
			for (Int i = 0; i < numVerIdx; i++)
#else
			for (Int i = 0; i<2; i++)
#endif
			{
#if JVECT_C0062_AFFINE_SIX_PARAM
				Int iTmpIdx = (bRecShape && i == 1 && b6Param == 0) ? (i + 1) : i;
				m_pcRdCost->setPredictor(acMvPred[iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
				TComPattern   tmpPattern;
				TComPattern*  pcPatternKey = &tmpPattern;
				pcPatternKey->m_pcCU = pcCU;
				pcPatternKey->m_uiPartAddr = uiPartAddr;
				pcPatternKey->m_mvPred = &acMvPred[iTmpIdx];
				iCostScale = 0;
				uiBitsBest += xMvdBits(pcPatternKey, acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer());
#else
				uiBitsBest += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMvTemp[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
				);
#endif
#else
				m_pcRdCost->setPredictor(acMvPred[i]);
#if JVET_D0123_ME_CTX_LUT_BITS
				TComPattern   tmpPattern;
				TComPattern*  pcPatternKey = &tmpPattern;
				pcPatternKey->m_pcCU = pcCU;
				pcPatternKey->m_uiPartAddr = uiPartAddr;
				pcPatternKey->m_mvPred = &acMvPred[i];
				iCostScale = 0;
				uiBitsBest += xMvdBits(pcPatternKey, acMvTemp[i].getHor(), acMvTemp[i].getVer());
#else
				uiBitsBest += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
					acMvTemp[i].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[i].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
					acMvTemp[i].getHor(), acMvTemp[i].getVer()
#endif
#if VCEG_AZ07_IMV
					, false
#endif
				);
#endif
#endif
			}
#if JVET_D0123_ME_CTX_LUT_BITS
			uiCostBest = (UInt)(floor(fWeight * (Double)uiCostBest) + (Double)m_pcRdCost->getCostSearch(uiBitsBest));
#else
			uiCostBest = (UInt)(floor(fWeight * (Double)uiCostBest) + (Double)m_pcRdCost->getCost(uiBitsBest));
#endif

#if COM16_C1016_AFFINE
			memcpy(acMv, acMvTemp, sizeof(TComMv) * 3);
#endif

  Int iIterTime = bBi ? 5 : 7;
#if AFFINE_REDUCE_ITERATION
	if (b6Param)
	{
		iIterTime = bBi ? 3 : 5;
	}
#endif

  for ( Int iter=0; iter<iIterTime; iter++ )    // iterate loop
  {
/*********************************************************************************
 *                         use gradient to update mv
 *********************************************************************************/
    // get Error Matrix
    Pel* pOrg  = pcYuv->getAddr(COMPONENT_Y, uiPartAddr);
    Pel* pPred = pPredYuv->getAddr(COMPONENT_Y, uiPartAddr);
    for ( Int j=0; j<iRoiHeight; j++ )
    {
      for ( Int i=0; i<iRoiWidth; i++ )
      {
        piError[i + j * iRoiWidth] = pOrg[i] - pPred[i];
      }
      pOrg  += iOrgStride;
      pPred += iPredStride;
    }

				// sobel x direction
				// -1 0 1
				// -2 0 2
				// -1 0 1
				pPred = pPredYuv->getAddr(COMPONENT_Y, uiPartAddr);
				for (Int j = 1; j < iRoiHeight - 1; j++)
				{
					for (Int k = 1; k < iRoiWidth - 1; k++)
					{
						Int iCenter = j*iPredStride + k;
						pdDerivate[0][j*iRoiWidth + k] = (Double)(pPred[iCenter + 1 - iPredStride] - pPred[iCenter - 1 - iPredStride]
							+ (pPred[iCenter + 1] << 1) - (pPred[iCenter - 1] << 1)
							+ pPred[iCenter + 1 + iPredStride] - pPred[iCenter - 1 + iPredStride]) / 8;
					}
					pdDerivate[0][j*iRoiWidth] = pdDerivate[0][j*iRoiWidth + 1];
					pdDerivate[0][j*iRoiWidth + iRoiWidth - 1] = pdDerivate[0][j*iRoiWidth + iRoiWidth - 2];
				}

				pdDerivate[0][0] = pdDerivate[0][iRoiWidth + 1];
				pdDerivate[0][iRoiWidth - 1] = pdDerivate[0][iRoiWidth + iRoiWidth - 2];
				pdDerivate[0][(iRoiHeight - 1)*iRoiWidth] = pdDerivate[0][(iRoiHeight - 2)*iRoiWidth + 1];
				pdDerivate[0][(iRoiHeight - 1)*iRoiWidth + iRoiWidth - 1] = pdDerivate[0][(iRoiHeight - 2)*iRoiWidth + (iRoiWidth - 2)];

				for (Int j = 1; j < iRoiWidth - 1; j++)
				{
					pdDerivate[0][j] = pdDerivate[0][iRoiWidth + j];
					pdDerivate[0][(iRoiHeight - 1)*iRoiWidth + j] = pdDerivate[0][(iRoiHeight - 2)*iRoiWidth + j];
				}

				// sobel y direction
				// -1 -2 -1
				//  0  0  0
				//  1  2  1
				pPred = pPredYuv->getAddr(COMPONENT_Y, uiPartAddr);
				for (Int k = 1; k < iRoiWidth - 1; k++)
				{
					for (Int j = 1; j < iRoiHeight - 1; j++)
					{
						Int iCenter = j*iPredStride + k;
						pdDerivate[1][j*iRoiWidth + k] = (Double)(pPred[iCenter + iPredStride - 1] - pPred[iCenter - iPredStride - 1]
							+ (pPred[iCenter + iPredStride] << 1) - (pPred[iCenter - iPredStride] << 1)
							+ pPred[iCenter + iPredStride + 1] - pPred[iCenter - iPredStride + 1]) / 8;
					}
					pdDerivate[1][k] = pdDerivate[1][iRoiWidth + k];
					pdDerivate[1][(iRoiHeight - 1) * iRoiWidth + k] = pdDerivate[1][(iRoiHeight - 2) * iRoiWidth + k];
				}

				pdDerivate[1][0] = pdDerivate[1][iRoiWidth + 1];
				pdDerivate[1][iRoiWidth - 1] = pdDerivate[1][iRoiWidth + iRoiWidth - 2];
				pdDerivate[1][(iRoiHeight - 1)*iRoiWidth] = pdDerivate[1][(iRoiHeight - 2)*iRoiWidth + 1];
				pdDerivate[1][(iRoiHeight - 1)*iRoiWidth + iRoiWidth - 1] = pdDerivate[1][(iRoiHeight - 2)*iRoiWidth + (iRoiWidth - 2)];

				for (Int j = 1; j < iRoiHeight - 1; j++)
				{
					pdDerivate[1][j*iRoiWidth] = pdDerivate[1][j*iRoiWidth + 1];
					pdDerivate[1][j*iRoiWidth + iRoiWidth - 1] = pdDerivate[1][j*iRoiWidth + iRoiWidth - 2];
				}

				// solve delta x and y
				for (Int m = 0; m != iParaNum; m++)
				{
					for (Int n = 0; n != iParaNum; n++)
					{
						pdEqualCoeff[m][n] = 0.0;
					}
				}

				for (Int j = 0; j != iRoiHeight; j++)
				{
					for (Int k = 0; k != iRoiWidth; k++)
					{
						Int iIdx = j * iRoiWidth + k;
#if JVECT_C0062_AFFINE_SIX_PARAM
						Double dC[6];
						if (b6Param)
						{
							dC[0] = pdDerivate[0][iIdx];
							dC[1] = k * pdDerivate[0][iIdx];
							dC[2] = pdDerivate[1][iIdx];
							dC[3] = k * pdDerivate[1][iIdx];
							dC[4] = j * pdDerivate[0][iIdx];
							dC[5] = j * pdDerivate[1][iIdx];
						}
						else
						{
							dC[0] = pdDerivate[0][iIdx];
							dC[1] = k * pdDerivate[0][iIdx] + j * pdDerivate[1][iIdx];
							dC[2] = pdDerivate[1][iIdx];
							dC[3] = j * pdDerivate[0][iIdx] - k * pdDerivate[1][iIdx];
						}
#else
						Double dC[4];
						dC[0] = pdDerivate[0][iIdx];
						dC[1] = k * pdDerivate[0][iIdx] + j * pdDerivate[1][iIdx];
						dC[2] = pdDerivate[1][iIdx];
						dC[3] = j * pdDerivate[0][iIdx] - k * pdDerivate[1][iIdx];
#endif

#if JVECT_C0062_AFFINE_SIX_PARAM
						for (Int col = 0; col < iAffineParamNum; col++)
						{
							for (Int row = 0; row < iAffineParamNum; row++)
							{
#else
						for (Int col = 0; col<4; col++)
						{
							for (Int row = 0; row<4; row++)
							{
#endif
								pdEqualCoeff[col + 1][row] += dC[col] * dC[row];
							}
#if JVECT_C0062_AFFINE_SIX_PARAM
							pdEqualCoeff[col + 1][iAffineParamNum] += (Double)(piError[iIdx] * dC[col]);
#else
							pdEqualCoeff[col + 1][4] += (Double)(piError[iIdx] * dC[col]);
#endif
						}
							}
						}
#if JVECT_C0062_AFFINE_SIX_PARAM
				solveEqual(pdEqualCoeff, iAffineParamNum, dAffinePara);
#else
				solveEqual(pdEqualCoeff, 4, dAffinePara);
#endif

				// convert to delta mv
#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param)
				{
					//for MV0
					dDeltaMv[0] = dAffinePara[0];
					dDeltaMv[2] = dAffinePara[2];
					//for MV1
					dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];//ZF
					dDeltaMv[3] = dAffinePara[3] * iRoiWidth + dAffinePara[2];
					//for MV2
					dDeltaMv[4] = dAffinePara[4] * iRoiHeight + dAffinePara[0];//ZF
					dDeltaMv[5] = dAffinePara[5] * iRoiHeight + dAffinePara[2];
				}
				else
				{
					if (bRecShape)
					{
						dDeltaMv[0] = dAffinePara[0];
						dDeltaMv[2] = dAffinePara[2];

						dDeltaMv[1] = dAffinePara[3] * iRoiHeight + dAffinePara[0];
						dDeltaMv[3] = dAffinePara[1] * iRoiHeight + dAffinePara[2];
					}
					else
					{
						dDeltaMv[0] = dAffinePara[0];
						dDeltaMv[2] = dAffinePara[2];

						dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];//ZF
						dDeltaMv[3] = -dAffinePara[3] * iRoiWidth + dAffinePara[2];
					}
				}
#else
				dDeltaMv[0] = dAffinePara[0];
				dDeltaMv[2] = dAffinePara[2];

				dDeltaMv[1] = dAffinePara[1] * iRoiWidth + dAffinePara[0];
				dDeltaMv[3] = -dAffinePara[3] * iRoiWidth + dAffinePara[2];
#endif

#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param)
				{
					acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[2].setHor((Int)(dDeltaMv[4] * 4 + SIGN(dDeltaMv[4]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					acDeltaMv[2].setVer((Int)(dDeltaMv[5] * 4 + SIGN(dDeltaMv[5]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				}
				else
				{
					if (bRecShape)
					{
						acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[2].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[2].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					}
					else
					{
						acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
						acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
					}
				}
#else
				acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
				acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5) << VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE);
#endif
#else
				acDeltaMv[0].setHor((Int)(dDeltaMv[0] * 4 + SIGN(dDeltaMv[0]) * 0.5));
				acDeltaMv[0].setVer((Int)(dDeltaMv[2] * 4 + SIGN(dDeltaMv[2]) * 0.5));
				acDeltaMv[1].setHor((Int)(dDeltaMv[1] * 4 + SIGN(dDeltaMv[1]) * 0.5));
				acDeltaMv[1].setVer((Int)(dDeltaMv[3] * 4 + SIGN(dDeltaMv[3]) * 0.5));
#endif

				Bool bAllZero = false;
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (bRecShape && i && b6Param == 0) ? (i + 1) : i;
					if (acDeltaMv[iTmpIdx].getHor() != 0 || acDeltaMv[iTmpIdx].getVer() != 0)
#else
					if (acDeltaMv[i].getHor() != 0 || acDeltaMv[i].getVer() != 0)
#endif
					{
						bAllZero = false;
						break;
					}
					bAllZero = true;
				}

				if (bAllZero)
					break;

				// do motion compensation with updated mv
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (bRecShape && i && b6Param == 0) ? (i + 1) : i;
					acMvTemp[iTmpIdx] += acDeltaMv[iTmpIdx];
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[iTmpIdx]);
#endif
#else
					acMvTemp[i] += acDeltaMv[i];
#if !AFFINE_CLIP_REMOVE
					pcCU->clipMv(acMvTemp[i]);
#endif
#endif
				}

#if JVECT_C0062_AFFINE_SIX_PARAM
				if (b6Param == 0)
				{
					if (bRecShape)
					{
						Int vx1 = (acMvTemp[2].getVer() - acMvTemp[0].getVer())*iRoiWidth / iRoiHeight + acMvTemp[0].getHor();
						Int vy1 = -(acMvTemp[2].getHor() - acMvTemp[0].getHor())*iRoiWidth / iRoiHeight + acMvTemp[0].getVer();
						acMvTemp[1].set(vx1, vy1);
#if !AFFINE_CLIP_REMOVE
						pcCU->clipMv(acMvTemp[1]);
#endif
					}
					else
					{
						Int vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
						Int vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
						acMvTemp[2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
						pcCU->clipMv(acMvTemp[2]);
#endif
					}
				}
#else
				vx2 = -(acMvTemp[1].getVer() - acMvTemp[0].getVer()) * iRoiHeight / iRoiWidth + acMvTemp[0].getHor();
				vy2 = (acMvTemp[1].getHor() - acMvTemp[0].getHor()) * iRoiHeight / iRoiWidth + acMvTemp[0].getVer();
				acMvTemp[2].set(vx2, vy2);
#if !AFFINE_CLIP_REMOVE
				pcCU->clipMv(acMvTemp[2]);
#endif
#endif

				xPredAffineBlk(COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, acMvTemp, iRoiWidth, iRoiHeight, pPredYuv, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA)
#if JVECT_C0062_AFFINE_SIX_PARAM
					, b6Param
#endif
				);

				// get error
				uiCostTemp = m_pcRdCost->getDistPart(pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pPredYuv->getAddr(COMPONENT_Y, uiPartAddr), pPredYuv->getStride(COMPONENT_Y), pcYuv->getAddr(COMPONENT_Y, uiPartAddr), pcYuv->getStride(COMPONENT_Y), iRoiWidth, iRoiHeight, COMPONENT_Y, DF_HADS);

				// get cost with mv
				m_pcRdCost->setCostScale(0);
				uiBitsTemp = ruiBits;
#if JVECT_C0062_AFFINE_SIX_PARAM
				for (Int i = 0; i < numVerIdx; i++)
#else
				for (Int i = 0; i<2; i++)
#endif
				{
#if JVECT_C0062_AFFINE_SIX_PARAM
					Int iTmpIdx = (i && bRecShape && b6Param == 0) ? (i + 1) : i;
					m_pcRdCost->setPredictor(acMvPred[iTmpIdx]);
#if JVET_D0123_ME_CTX_LUT_BITS
					TComPattern   tmpPattern;
					TComPattern*  pcPatternKey = &tmpPattern;
					pcPatternKey->m_pcCU = pcCU;
					pcPatternKey->m_uiPartAddr = uiPartAddr;
					pcPatternKey->m_mvPred = &acMvPred[iTmpIdx];
					iCostScale = 0;
					uiBitsTemp += xMvdBits(pcPatternKey, acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer());
#else
					uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMvTemp[iTmpIdx].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[iTmpIdx].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMvTemp[iTmpIdx].getHor(), acMvTemp[iTmpIdx].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
					);
#endif
#else
					m_pcRdCost->setPredictor(acMvPred[i]);
#if JVET_D0123_ME_CTX_LUT_BITS
					TComPattern   tmpPattern;
					TComPattern*  pcPatternKey = &tmpPattern;
					pcPatternKey->m_pcCU = pcCU;
					pcPatternKey->m_uiPartAddr = uiPartAddr;
					pcPatternKey->m_mvPred = &acMvPred[i];
					iCostScale = 0;
					uiBitsTemp += xMvdBits(pcPatternKey, acMvTemp[i].getHor(), acMvTemp[i].getVer());
#else 
					uiBitsTemp += m_pcRdCost->getBits(
#if VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
						acMvTemp[i].getHor() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE, acMvTemp[i].getVer() >> VCEG_AZ07_MV_ADD_PRECISION_BIT_FOR_STORE
#else
						acMvTemp[i].getHor(), acMvTemp[i].getVer()
#endif
#if VCEG_AZ07_IMV
						, false
#endif
					);
#endif
#endif
				}
#if JVET_D0123_ME_CTX_LUT_BITS
				uiCostTemp = (UInt)(floor(fWeight * (Double)uiCostTemp) + (Double)m_pcRdCost->getCostSearch(uiBitsTemp));
#else
				uiCostTemp = (UInt)(floor(fWeight * (Double)uiCostTemp) + (Double)m_pcRdCost->getCost(uiBitsTemp));
#endif

				// store best cost and mv
				if (uiCostTemp < uiCostBest)
				{
					uiCostBest = uiCostTemp;
					uiBitsBest = uiBitsTemp;
					memcpy(acMv, acMvTemp, sizeof(TComMv) * 3);
				}
					}

			// free buffer
			for (Int i = 0; i<iParaNum; i++)
				delete[]pdEqualCoeff[i];
			delete[]pdEqualCoeff;

			ruiBits = uiBitsBest;
			ruiCost = uiCostBest;
				}
#endif
#endif

#if MRG_CONSIDER_COMPLEX
		Int TEncSearch::getComplexIndexCodeBits(UInt numValidMergeComplexCand, UInt uiMergeComplexIndex)
		{
			return 3;
		}

		Void TEncSearch::xMergeComplexEstimation(TComDataCU* pcCU, TComYuv*  pcYuvOrg, Int iPUIdx, UInt& uiComplexInterDir, NodeInfo& NodeMotionInfoP, UInt&  uiMergeComplexIndex, UInt& ruiCost)
		{
			UInt uiAbsPartIdx = 0;
			Int iWidth = 0;
			Int iHeight = 0;
			pcCU->getPartIndexAndSize(0, uiAbsPartIdx, iWidth, iHeight);
			ruiCost = MAX_UINT;
			UChar NumCand = pcCU->getNumMergeComplexCand(uiAbsPartIdx);
			pcCU->setMergeComplexFlagSubParts(true, uiAbsPartIdx);
			NodeInfo NodeMotionInfo = pcCU->getNodeMotionInfo(uiAbsPartIdx);
			TComMvField  cMvFieldNeighbourNodes[MRG_NODE_NUM << 1];
			UChar InterDirNeighbourNodes[MRG_NODE_NUM];
			Int numValidMergeNode[MRG_COMPLEX_NUM_CANDS << 1];
			Int RealComplexIndex[MRG_COMPLEX_NUM_CANDS] = { 0 };

			for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
			{
				cMvFieldNeighbourNodes[i] = NodeMotionInfo.MvFieldNbNodes[i];
			}
			for (Int i = 0; i < MRG_NODE_NUM; i++)
			{
				InterDirNeighbourNodes[i] = NodeMotionInfo.InterDirNbNodes[i];
			}
			for (Int i = 0; i < (MRG_COMPLEX_NUM_CANDS << 1); i++)
			{
				numValidMergeNode[i] = NodeMotionInfo.numValidNodes[i];
			}
			for (Int i = 0; i < MRG_COMPLEX_NUM_CANDS; i++)
			{
				RealComplexIndex[i] = NodeMotionInfo.RealComplexIdx[i];
			}

			for (UInt uiMergeComplexCand = 0; uiMergeComplexCand < NumCand; ++uiMergeComplexCand)
			{
				UInt uiCostCand = MAX_UINT;
				UChar uhInterDirComplex = 0;
				// pre process the control motion information
				TComMvField  cMvFieldNeighbourNodesTemp[MRG_NODE_NUM << 1];
				for (Int i = 0; i < MRG_NODE_NUM; i++)
				{
					cMvFieldNeighbourNodesTemp[(i << 1)].setMvField(cMvFieldNeighbourNodes[(i << 1)].getMv(), cMvFieldNeighbourNodes[(i << 1)].getRefIdx());
					cMvFieldNeighbourNodesTemp[(i << 1) + 1].setMvField(cMvFieldNeighbourNodes[(i << 1) + 1].getMv(), cMvFieldNeighbourNodes[(i << 1) + 1].getRefIdx());
				}
				pcCU->preproccessMotionFieldNodes(RealComplexIndex, uiMergeComplexCand, InterDirNeighbourNodes, uhInterDirComplex, numValidMergeNode, cMvFieldNeighbourNodesTemp);
				for (Int i = 0; i < (MRG_NODE_NUM << 1); i++)
				{
					NodeMotionInfo.MvFieldNbNodes[i] = cMvFieldNeighbourNodesTemp[i];
				}
				pcCU->setNodeMotionInfoSubParts(NodeMotionInfo, uiAbsPartIdx);
				pcCU->setMotionInfoField(uiMergeComplexCand, uhInterDirComplex, SIZE_2Nx2N, uiAbsPartIdx, iPUIdx);
				pcCU->setMergeComplexIndexSubParts(uiMergeComplexCand, uiAbsPartIdx);
				xGetInterPredictionError(pcCU, pcYuvOrg, iPUIdx, uiCostCand, m_pcEncCfg->getUseHADME());

				if (uiCostCand < ruiCost)
				{
					ruiCost = uiCostCand;
					uiComplexInterDir = uhInterDirComplex;
					NodeMotionInfoP = NodeMotionInfo;
					uiMergeComplexIndex = uiMergeComplexCand;
				}
			}
			pcCU->setMergeComplexFlagSubParts(false, uiAbsPartIdx);
		}
#endif

		//! \}