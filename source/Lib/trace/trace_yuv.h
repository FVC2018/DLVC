#ifndef __TRACE_YUV_H_
#define __TRACE_YUV_H_

#include "trace_common.h"

#if TRACE_YUV

class TComDataCU;
class TComPicYuv;
class TComYuv;
class TComROI;
class TComTU;

#include "../TLibCommon/CommonDef.h"
#include <stdio.h>

void TRACE_CU_PART(TComYuv *picSrc, TComDataCU* pcCU, const char *fname, float scale = 1, int norm = 0);
void TRACE_CU_PART_DIFF(TComYuv *picSrcA, TComYuv *picSrcB, TComDataCU* pcCU, const char *fname);
void TRACE_CU_PART_SUM(TComYuv *picSrcA, TComYuv *picSrcB, TComDataCU* pcCU, const char *fname, Int shift_right = 1);
TComPicYuv* TRACE_GET_PIC_YUV(const char *fname);
TComPicYuv* TRACE_GET_PIC_YUV(const char *fname, int width, int height);
void CREATE_PIC(int width, int height, const char *fname);
void TRACE_PIC(TComPicYuv *picSrc, const char *fname);
void TRACE_FLOAT(float *src, int width, int height, const char* fname, short bias, float scale);
void TRACE_INT(int *src, int width, int height, const char* fname, int shift = 0, int bias = 0, int shift_left = 0);

template <typename TPel> void TRACE_REF_SAMPLES(TPel *pSrc, Int iWidth, Int iHeight, Int iSrcStride, TComDataCU* pcCU, UInt uiAbsPartIdx, Int iVerHor, const char *fname, float scale = 1, int norm = 0, int bias = 0);

void TRACE_TEXT(int *src, int width, int height, const char* fname, int shift = 0, int bias = 0, int shift_left = 0);

void MARK_PARTIONING(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, const char *fname, int lineColor = 0, bool putChar = 1, int refList = 0);
void MARK_PARTIONING_LCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt startDepth, const char *fname, int lineColor = 0, bool putChar = 1, int refList = 0);
void TEXT_CU(TComDataCU* pcCU, UInt uiAbsPartIdx, int width, int height, const char* str, const char *fname, int drawCUborder = false, int lineColor = 0);
Bool isIdenticalMotionMy(TComDataCU* pcCU, UInt PartAddr);
void DUMP_ALL_YUV();

#else//TRACE_YUV
#define TRACE_CU_PART
#define TRACE_CU_PART_DIFF
#define TRACE_CU_PART_SUB
#define TRACE_PIC
#define TRACE_PIC_PART
#define TRACE_PIC_PART_SUB
#define TRACE_GET_PIC_YUV   0
#define TRACE_GET_PIC_YUV   0
#define TRACE_FLOAT
#define TRACE_INT
#define ZERO_CU_PART
#define CREATE_PIC
#define MARK_PARTIONING
#define DUMP_ALL_YUV

#endif //TRACE_YUV

#endif //__TRACE_YUV_H_
