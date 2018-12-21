
#include "trace_flags.h"

#if  TRACE_FLAGS

#include "TLibCommon/TComDataCU.h"
#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TypeDef.h"
#include "TLibCommon/TComTU.h"
#include <unordered_map>
#include <memory>
#include <stdio.h>
#include <fstream>

static std::unordered_map<std::string, std::unique_ptr<std::ofstream> >  gpTextStreams;
unsigned int gNumFrame = 0;
unsigned int gNumCu = 0;
unsigned int gNumCtu = 0;

std::string getStatFilePath(const std::string& name)
{
	std::string statPath = DebugOptionList::StatPath.getString();
	statPath = statPath.empty() ? statPath : statPath + std::string("\\");
	std::string fullPath = statPath + std::string(name);
	return fullPath;
}

std::pair<std::ostream*, bool> getTextStream(const std::string& name, std::ios_base::openmode open_mode)
{
	auto it = gpTextStreams.find(name);
	if (it == gpTextStreams.end())
	{
		std::string fullpath = getStatFilePath(name);
		std::unique_ptr<std::ofstream> stream = std::make_unique<std::ofstream>(fullpath, open_mode);
		if (!stream->is_open())
		{
			std::cout << "Unable to open text stream: " << fullpath << std::endl;
			exit(0);
		}
		auto val = gpTextStreams.insert(std::make_pair(name, std::move(stream)));
		return  std::make_pair(val.first->second.get(), true);
	}
	return std::make_pair(it->second.get(), false);
}

void DUMP_HEADER(TComDataCU* pcCU, UInt uiAbsPartIdx, std::ostream& stream)
{
	if (pcCU->isIntra(uiAbsPartIdx))
	{

		print_args_(stream, "NumFrame"
			, "NumCtu"
			, "NumCu"
			, "POC"
			, "SliceType"
            , "FrameQp"
			, "uiAbsPartIdx"
			, "CuX", "CuY"
			, "CuWidth", "CuHeight"
#if !JVET_C0024_QTBT
			, "TuX", "TuY"
			, "TuWidth", "TuHeight"
#endif

			, "Skip"
#if VCEG_AZ07_FRUC_MERGE
			, "FRUCMgrMode"
#endif
			, "Cbf"
#if VCEG_AZ05_INTRA_MPI
			, "MPIIdx"
#endif
#if COM16_C1046_PDPC_INTRA
			, "PDPCIdx"
#endif
			, "IPCMFlag"
			, "PredMode"
#if !JVET_C0024_QTBT
			, "PartitionSize"
#endif
			, "IntraDir"
#if RRC_DWDIP
			, "DWDIPIdx"
#endif
#if RRC_INTRA_SHARPENING
			, "IntraSharpFlag"
#endif
#if COM16_C1044_NSST
			, "ROTIdx"
#endif
#if COM16_C806_EMT
			, "EmtCuFlag"
			, "EmtTuIdx"
#endif
			, "TransformSkip"
		);
	}
	else
	{
		print_args_(stream, "NumFrame"
			, "NumCtu"
			, "NumCu"
			, "POC"
			, "SliceType"
            , "FrameQp"
            , "uiAbsPartIdx"
			, "CuX", "CuY"
			, "CuWidth", "CuHeight"
#if !JVET_C0024_QTBT
			, "TuX", "TuY"
			, "TuWidth", "TuHeight"
#endif
			, "Skip"
#if VCEG_AZ07_FRUC_MERGE
			, "FRUCMgrMode"
#endif
#if COM16_C1016_AFFINE
			, "AffineFlag"
#endif
			, "Cbf"
			, "MergeFlag"
			, "MergeIndex"
#if COM16_C806_OBMC
			, "OBMCFlag"
#endif
#if SIF_CU_LEVEL
			, "SIFFlag"
#endif
#if VCEG_AZ06_IC
			, "ICFlag"
#endif
			, "PredMode"
#if !JVET_C0024_QTBT
			, "PartitionSize"
#endif
			, "InterDir"

			, "RefFrmIdx_R0"
			, "RefFrmIdx_R1"
			, "MVPIdx_R0"
			, "MVPIdx_R1"
			, "CUMvField_R0_HOR"
			, "CUMvField_R0_VER"
			, "CUMvField_R1_HOR"
			, "CUMvField_R1_VER"
#if VCEG_AZ07_IMV
			, "iMVFlag"
#endif
#if COM16_C1044_NSST
			, "ROTIdx"
#endif
#if COM16_C806_EMT
			, "EmtCuFlag"
			, "EmtTuIdx"
#endif
			, "TransformSkip"
		);


	}
}

void DUMP_PU_FLAGS_INTER(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, Int uiAbsCuIdx)
{
	if (pcCU->getPredictionMode(uiAbsPartIdx) != MODE_INTER)
		return;

#if !JVET_C0024_QTBT
	PartSize ePartSize = pcCU->getPartitionSize(uiAbsPartIdx);
	UInt uiNumPU = (ePartSize == SIZE_2Nx2N ? 1 : (ePartSize == SIZE_NxN ? 4 : 2));
	UInt uiDepth = pcCU->getDepth(uiAbsPartIdx);
	UInt uiPUOffset = (g_auiPUOffset[UInt(ePartSize)] << ((pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth) << 1)) >> 4;
#else
	UInt uiNumPU = 1;
	UInt uiPUOffset = 0;
	uiAbsPartIdx = 0;
#endif
	for (UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset)
	{
#if TRACE_FLAGS
		DUMP_CU_FLAGS_INFO(name, pcCU, uiSubPartIdx, uiSubPartIdx, uiAbsPartIdx);
#endif
	}
}

void DUMP_BLOCK_FLAGS_INFO(const char* name, TComDataCU* pcCU,
#if !JVET_C0024_QTBT
	TComTU* pTU,
#endif
	UInt uiAbsPartIdx, UInt absCtuIdx, Int absCuIdx)
{
#if JVET_C0024_QTBT
	if (pcCU->getTextType() != CHANNEL_TYPE_LUMA)
		return;
#endif
	UInt cuX = pcCU->getCUPelX();
	UInt cuY = pcCU->getCUPelY();
	if (absCuIdx > 0)
	{
		cuX += g_auiRasterToPelX[g_auiZscanToRaster[absCuIdx]];
		cuY += g_auiRasterToPelY[g_auiZscanToRaster[absCuIdx]];
	}
#if !JVET_C0024_QTBT
	Int tuX = -1;
	Int tuY = -1;
	Int tuWidth = -1;
	Int tuHeight = -1;

	if (pTU)
	{
		tuX = cuX + g_auiRasterToPelX[g_auiZscanToRaster[pTU->GetRelPartIdxTU()]];
		tuY = cuY + g_auiRasterToPelY[g_auiZscanToRaster[pTU->GetRelPartIdxTU()]];
		tuWidth = pTU->getRect(COMPONENT_Y).width;
		tuHeight = pTU->getRect(COMPONENT_Y).height;
	}

#endif
	std::string sname(name);
	sname += pcCU->isIntra(uiAbsPartIdx) ? std::string("_intra.csv") : std::string("_inter.csv");

	auto stream = getTextStream(sname);
	if (stream.second)
	{
		DUMP_HEADER(pcCU, uiAbsPartIdx, *(stream.first));
	}

	if (pcCU->isIntra(uiAbsPartIdx))
	{
		print_args(*(stream.first)
			, gNumFrame
			, gNumCtu
			, gNumCu
			, pcCU->getSlice()->getPOC()
			, pcCU->getSlice()->getSliceType()
            , pcCU->getSlice()->getSliceQp()
			, absCtuIdx
			, cuX, cuY
			, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx)
#if !JVET_C0024_QTBT
			, tuX, tuY
			, tuWidth, tuHeight
#endif
			, pcCU->getSkipFlag(uiAbsPartIdx)
			///SKIP
#if VCEG_AZ07_FRUC_MERGE
			, pcCU->getFRUCMgrMode(uiAbsPartIdx)
#endif
			, pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y)
#if VCEG_AZ05_INTRA_MPI
			, pcCU->getMPIIdx(uiAbsPartIdx)
#endif
#if COM16_C1046_PDPC_INTRA
			, pcCU->getPDPCIdx(uiAbsPartIdx)
#endif
			, pcCU->getIPCMFlag(uiAbsPartIdx)

			, pcCU->getPredictionMode(uiAbsPartIdx)
#if !JVET_C0024_QTBT
			, pcCU->getPartitionSize(uiAbsPartIdx)
#endif

			, pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx)
#if RRC_DWDIP
			, pcCU->getDWDIPIdx(uiAbsPartIdx)
#endif
#if RRC_INTRA_SHARPENING
			, pcCU->getIntraSharpFlag(uiAbsPartIdx)
#endif
#if COM16_C1044_NSST
			, pcCU->getROTIdx(CHANNEL_TYPE_LUMA, uiAbsPartIdx)
#endif
#if COM16_C806_EMT
			, pcCU->getEmtCuFlag(uiAbsPartIdx)
			, pcCU->getEmtTuIdx(uiAbsPartIdx)
#endif
			, pcCU->getTransformSkip(uiAbsPartIdx, COMPONENT_Y)

		);
	}
	else
	{
		print_args(*(stream.first)
			, gNumFrame
			, gNumCtu
			, gNumCu
			, pcCU->getSlice()->getPOC()
			, pcCU->getSlice()->getSliceType()
            , pcCU->getSlice()->getSliceQp()
			, absCtuIdx
			, cuX, cuY
			, pcCU->getWidth(uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx)
#if !JVET_C0024_QTBT
			, tuX, tuY
			, tuWidth, tuHeight
#endif
			, pcCU->getSkipFlag(uiAbsPartIdx)
#if VCEG_AZ07_FRUC_MERGE
			, pcCU->getFRUCMgrMode(uiAbsPartIdx)
#endif
#if COM16_C1016_AFFINE
			, pcCU->getAffineFlag(uiAbsPartIdx)
#endif
			, pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y)
			, pcCU->getMergeFlag( uiAbsPartIdx )
			, pcCU->getMergeIndex(uiAbsPartIdx)
#if COM16_C806_OBMC
			, pcCU->getOBMCFlag(uiAbsPartIdx)
#endif
#if SIF_CU_LEVEL
			, pcCU->getSIFFlag(uiAbsPartIdx)
#endif
#if VCEG_AZ06_IC
			, pcCU->getICFlag(uiAbsPartIdx)
#endif
			, pcCU->getPredictionMode(uiAbsPartIdx)
#if !JVET_C0024_QTBT
			, pcCU->getPartitionSize(uiAbsPartIdx)
#endif
			, pcCU->getInterDir(uiAbsPartIdx)

			, pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiAbsPartIdx)
			, pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiAbsPartIdx)
			, pcCU->getMVPIdx(REF_PIC_LIST_0, uiAbsPartIdx)
			, pcCU->getMVPIdx(REF_PIC_LIST_1, uiAbsPartIdx)
			, pcCU->getCUMvField(REF_PIC_LIST_0)->getMvd(uiAbsPartIdx).getHor()
			, pcCU->getCUMvField(REF_PIC_LIST_0)->getMvd(uiAbsPartIdx).getVer()
			, pcCU->getCUMvField(REF_PIC_LIST_1)->getMvd(uiAbsPartIdx).getHor()
			, pcCU->getCUMvField(REF_PIC_LIST_1)->getMvd(uiAbsPartIdx).getVer()
#if VCEG_AZ07_IMV
			, pcCU->getiMVFlag(uiAbsPartIdx)
#endif
#if COM16_C1044_NSST
			, pcCU->getROTIdx(CHANNEL_TYPE_LUMA, 0)
#endif
#if COM16_C806_EMT
			, pcCU->getEmtCuFlag(uiAbsPartIdx)
			, pcCU->getEmtTuIdx(uiAbsPartIdx)
#endif
			, pcCU->getTransformSkip(uiAbsPartIdx, COMPONENT_Y)
		);
	}


}


void DUMP_CU_FLAGS_INFO(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absCtuIdx, Int absCuIdx)
{
	DUMP_BLOCK_FLAGS_INFO(name, pcCU,
#if !JVET_C0024_QTBT
		nullptr,
#endif
		uiAbsPartIdx, absCtuIdx, absCuIdx);
}

#if !JVET_C0024_QTBT
void DUMP_TU_FLAGS_INFO(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absCtuIdx, TComTU& rTu, Int absCuIdx)
{
	DUMP_BLOCK_FLAGS_INFO(name, pcCU,
#if !JVET_C0024_QTBT
		&rTu,
#endif
		uiAbsPartIdx, absCtuIdx, absCuIdx);
}
#endif
#endif

