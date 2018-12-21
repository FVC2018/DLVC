#pragma once

#include "trace_common.h"


#if  TRACE_FLAGS 
#include "../TLibCommon/CommonDef.h"
class TComDataCU;
class TComPicYuv;
class TComYuv;
class TComROI;
class TComTU;

extern unsigned int gNumFrame;
extern unsigned int gNumCtu;
extern unsigned int gNumCu;

void DUMP_PU_FLAGS_INTER(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, Int uiAbsCuIdx);
void DUMP_CU_FLAGS_INFO(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absCtuId, Int absCuIdx = -1);
#if !JVET_C0024_QTBT
void DUMP_TU_FLAGS_INFO(const char* name, TComDataCU* pcCU, UInt uiAbsPartIdx, UInt absCtuIdx, TComTU& rTu, Int absCuIdx = -1);
#endif
std::pair<std::ostream*, bool>  getTextStream(const std::string& name, std::ios_base::openmode open_mode = std::ostream::out | std::ostream::trunc);

template<class Head>
void print_args_(std::ostream& s, Head&& head) {
s << +std::forward<Head>(head) << std::endl;
}

template<class Head, class... Tail>
void print_args_(std::ostream& s, Head&& head, Tail&&... tail) {
	s << +(std::forward<Head>(head)) << ",";
	print_args_(s, std::forward<Tail>(tail)...);
}

template<class... Args>
void print_args(std::ostream& s, Args&&... args) {
	print_args_(s, std::forward<Args>(args)...);
}

#endif


