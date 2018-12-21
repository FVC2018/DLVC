#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TypeDef.h"
#include "TLibCommon/TComTU.h"
#include "trace_yuv.h"

#if TRACE_YUV

#include "pel_access_extension.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>

#include "../TLibCommon/TComDataCU.h"
#include "../TLibCommon/TComList.h"
#include "../TLibCommon/TComPic.h"
#include "../TLibCommon/TComPicYuv.h"
#include "../TLibEncoder/TEncCu.h"
#include "../TLibVideoIO/TVideoIOYuv.h"
#if RRC_SELECTABLE_IF
#include "../TLibCommon/TComSharpeningFilter.h"
#endif

#include "u_vga_font.h"

#pragma warning(disable : 4996)
#pragma warning( disable : 4244 )

#define STREAMS_NUM     80
#define MAX_FNAME       255

#define GET_STREAM(fname)                 \
  FILE* stream = getStream(fname, "wb");  \
  if (!stream) return;

#define GET_STREAM_TXT(fname)             \
  FILE* stream = getStream(fname, "wt");  \
  if (!stream) return;

#define GET_YUV_CU(fname, pcCU)  TComPicYuvExt* picYuv = getPicYuv(fname, pcCU);

  
BitDepths     gBitDepthsTrace = {8,8};
TVideoIOYuv*  gpStreams[STREAMS_NUM] = {0};
char          gpStreamNames[STREAMS_NUM][MAX_FNAME] = {0};
int           gpStreamsHeader[STREAMS_NUM] = {0};
TComPicYuvExt* gpPicYuv[STREAMS_NUM] = {0};

#if SVC_EXTENSION
#define CHECK_LAYER(CU)                   \
  if (pcCU->getLayerId() != DUMP_LAYER_ID)   \
  return;
#else
#define CHECK_LAYER
#endif//SVC_EXTENSION

#define GET_STREAM_HEADER_WRITEN(fname)   \
  getStreamFlag(fname, gpStreamsHeader)

static void markPartitioning(TComDataCU* pcCU, UInt uiAbsPartIdx, TComPicYuv *pic);


TVideoIOYuv* getStream(const char *fname, const char *openKeys)
{
  int i = 0;
  int found_idx = -1;
  for (i = 0; i < STREAMS_NUM; i++)
  {
    if (!strcmp(fname, gpStreamNames[i]))
    {
      found_idx = i;
      break;
    }
  }
  if (found_idx >= 0)
  {
    return gpStreams[found_idx];
  }
  else
  {
    // find first empty slot
    for (i = 0; i < STREAMS_NUM; i++)
    {
      if (!strcmp("", gpStreamNames[i]))
      {
        found_idx = i;
        break;
      }
    }

    // open file
    if (found_idx >= 0)
    {
      strcpy(gpStreamNames[found_idx], fname);
      assert(!gpStreams[found_idx]);
      gpStreams[found_idx] = new TVideoIOYuv();
      Int       m_inputBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };         ///< bit-depth of input file
      Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };         ///< bit-depth of output file
      Int       m_MSBExtendedBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };      ///< bit-depth of input samples after MSB extension

      gpStreams[found_idx]->open((Char*)fname, true, m_inputBitDepth, m_outputBitDepth, m_MSBExtendedBitDepth);
      return gpStreams[found_idx];
    }
    else
    {
      // no space
      assert(!"No space to store log PCM streams, increase STREAMS_NUM");
    }
  }
  return 0;
}
int getStreamIdx(const char *fname)
{
  int i = 0;
  int found_idx = -1;
  for (i = 0; i < STREAMS_NUM; i++)
  {
    if (!strcmp(fname, gpStreamNames[i]))
    {
      found_idx = i;
      break;
    }
  }
  if (found_idx < 0)
  {
    // find first empty slot
    for (i = 0; i < STREAMS_NUM; i++)
    {
      if (!strcmp("", gpStreamNames[i]))
      {
        found_idx = i;
        break;
      }
    }
    // open file
    if (found_idx >= 0)
    {
      strcpy(gpStreamNames[found_idx], fname);
      assert(!gpStreams[found_idx]);
      gpStreams[found_idx] = new TVideoIOYuv();
      Int       m_inputBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };         ///< bit-depth of input file
      Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };         ///< bit-depth of output file
      Int       m_MSBExtendedBitDepth[MAX_NUM_CHANNEL_TYPE] = { 8, 8 };      ///< bit-depth of input samples after MSB extension

      gpStreams[found_idx]->open((Char*)fname, true, m_inputBitDepth, m_outputBitDepth, m_MSBExtendedBitDepth);
    }
    else
    {
      // no space
      assert(!"No space to store log PCM streams, increase STREAMS_NUM");
    }
  }
  return found_idx;
}

int *getStreamFlag(const char *fname, int pFlagArr[])
{
  int i = 0;
  int found_idx = -1;
  for (i = 0; i < STREAMS_NUM; i++)
  {
    if (!strcmp(fname, gpStreamNames[i]))
    {
      found_idx = i;
      break;
    }
  }
  if (found_idx >= 0)
    return &pFlagArr[found_idx];
  else
  {
    assert(!"file not found");
    return 0;
  }
}

TComPicYuvExt *getPicYuv(const char *fname, int width, int height)
{
  int idx = getStreamIdx(fname);
  if (idx >= 0)
  {
    if (!gpPicYuv[idx])
    {
      gpPicYuv[idx] = new TComPicYuvExt();
      assert(width != 0 && height != 0);
      gpPicYuv[idx]->create(width, height, CHROMA_420, width, height, 0, true);
      gpPicYuv[idx]->clean();
    }
    return gpPicYuv[idx];
  }
  return 0;
}

TComPicYuvExt *getPicYuv(const char *fname, TComPicYuv* picSrc)
{
  int idx = getStreamIdx(fname);
  if (idx >= 0)
  {
    if (!gpPicYuv[idx])
    {
      gpPicYuv[idx] = new TComPicYuvExt();

      gpPicYuv[idx]->create(picSrc->getWidth(COMPONENT_Y), picSrc->getHeight(COMPONENT_Y), 
        picSrc->getChromaFormat(), picSrc->getMarginX(COMPONENT_Y)-16, picSrc->getMarginY(COMPONENT_Y)-16, 0, true);
      gpPicYuv[idx]->clean();
    }
    return gpPicYuv[idx];
  }
  return 0;
}

TComPicYuvExt *getPicYuv(const char *fname, TComDataCU *pcCU)
{
  int width  = pcCU->getPic()->getPicYuvRec()->getWidth(COMPONENT_Y);
  int height = pcCU->getPic()->getPicYuvRec()->getHeight(COMPONENT_Y);
  int idx = getStreamIdx(fname);
  if (idx >= 0)
  {
    if (!gpPicYuv[idx])
    {
      gpPicYuv[idx] = new TComPicYuvExt();
      assert(width != 0 && height != 0);
      gpPicYuv[idx]->create(width, height, 
        pcCU->getSlice()->getSPS()->getChromaFormatIdc(), 
#if JVET_C0024_QTBT
        pcCU->getSlice()->getSPS()->getCTUSize(),
        pcCU->getSlice()->getSPS()->getCTUSize(),
#else
        pcCU->getSlice()->getSPS()->getMaxCUWidth(),
        pcCU->getSlice()->getSPS()->getMaxCUHeight(),
#endif
        pcCU->getSlice()->getSPS()->getMaxTotalCUDepth(), 
        true);

      gBitDepthsTrace = pcCU->getSlice()->getSPS()->getBitDepths();
      gpPicYuv[idx]->clean();
    }
    return gpPicYuv[idx];
  }
  return 0;
}

static void line_QPel(Pel* startAddr, int x0, int y0, int dx, int dy, int stride, int color, int block_width, int block_height, int norm2length = 1)
{
  if (dx == 0 && dy == 0)
    return;
  if (norm2length > 1)
  {
    // normalize to length
    float c = norm2length / (sqrt(dx*dx + dy*dy));
    dx = dx*c + 0.5;
    dy = dy*c + 0.5;
  }

  Pel *dst = startAddr;
  if (dx != 0)
  {
    float error = 0;
    float deltaerr = fabs((float)dy / dx);
    int inc = dx > 0 ? 1:-1;
    int inc_y = dy > 0 ? 1:-1;
    int y = y0; int x = x0;
    for (; x != x0+dx; x += inc)
    {
      int newpos =  x + y*stride;
      if (x <= block_width && x >= 0 && y <= block_height && y >= 0) dst[newpos] = color;
      error = error + deltaerr;
      if (error >= 0.5)
      {
        y = y + inc_y;
        error = error - 1.0;
      }
    }
    int newpos =  x + y*stride;
    if (x <= block_width && x >= 0 && y <= block_height && y >= 0)
      dst[newpos] = color;
  }
  else
  {
    int inc = dy > 0 ? 1:-1;
    for (int y = y0; y != y0+dy; y+=inc)
    {
      int newpos =  0 + y*stride;
      if (y <= block_height && y >= 0)
        dst[newpos] = color;
    }
  }
}

void AALine(int x0, int y0, int x1, int y1, Pel* startAddr, int stride, int color, int block_width, int block_height)
{
  x1 += x0;
  y1 += y0;
  Pel *finalAddr = startAddr + block_width + block_height*stride;
  int addr = (y0 * stride + x0) * 1;//4;
  int dx = x1 - x0;
  int dy = y1 - y0;
  int u, v, du, dv, uincr, vincr;
  /* By switching to (u,v), we combine all eight octants */
  if (abs(dx) > abs(dy))
  {
    /* Note: If this were actual C, these integers would be lost
    * at the closing brace.  That's not what I mean to do.  Do what
    * I mean. */
    du = abs(dx);
    dv = abs(dy);
    u = x1;
    v = y1;
    uincr = 1;//4;
    vincr = stride * 1;//4;
    if (dx < 0) uincr = -uincr;
    if (dy < 0) vincr = -vincr;
  }
  else
  {
    du = abs(dy);
    dv = abs(dx);
    u = y1;
    v = x1;
    uincr = stride * 1;//4;
    vincr = 1;//4;
    if (dy < 0) uincr = -uincr;
    if (dx < 0) vincr = -vincr;
  }
  int uend = u + 2 * du;
  int d = (2 * dv) - du;	    /* Initial value as in Bresenham's */
  int incrS = 2 * dv;	/* ?d for straight increments */
  int incrD = 2 * (dv - du);	/* ?d for diagonal increments */
  int twovdu = 0;	/* Numerator of distance; starts at 0 */
  double invD = color*1.0 / (2.0*sqrt(du*du + dv*dv));   /* Precomputed inverse denominator */
  double invD2du = 2.0 * (du*invD);   /* Precomputed constant */
  do
  {
    /* Note: this pseudocode doesn't ensure that the address is
    * valid, or that it even represents a pixel on the same side of
    * the screen as the adjacent pixel */
    if (startAddr + addr >= finalAddr)
      break;
    //DrawPixelD(addr, twovdu*invD);
    //DrawPixelD(addr + vincr, invD2du - twovdu*invD);
    //DrawPixelD(addr - vincr, invD2du + twovdu*invD);
    startAddr[addr] = twovdu*invD;
    startAddr[addr + vincr] = invD2du - twovdu*invD;
    startAddr[addr - vincr] = invD2du + twovdu*invD;

    if (d < 0)
    {
      /* choose straight (u direction) */
      twovdu = d + du;
      d = d + incrS;
    }
    else
    {
      /* choose diagonal (u+v direction) */
      twovdu = d - du;
      d = d + incrD;
      v = v + 1;
      addr = addr + vincr;
    }
    u = u + 1;
    addr = addr + uincr;
  } while (u < uend);
}

void DrawWuLine(Pel* startAddr, short X0, short Y0, short X1, short Y1,
  short BaseColor, short NumLevels, unsigned short IntensityBits, int stride, int block_width, int block_height, int norm2length = 1)
{
  if (X1 == 0 && Y1 == 0)
    return;
  if (norm2length > 1)
  {
    // normalize to length
    float c = norm2length / (sqrt(X1*X1 + Y1*Y1));
    X1 = X1*c + 0.5;
    Y1 = Y1*c + 0.5;
  }  
  X1 += X0;
  Y1 += Y0;
  unsigned short IntensityShift, ErrorAdj, ErrorAcc;
  unsigned short ErrorAccTemp, Weighting, WeightingComplementMask;
  short DeltaX, DeltaY, Temp, XDir;

  /* Make sure the line runs top to bottom */
  if (Y0 > Y1) {
    Temp = Y0; Y0 = Y1; Y1 = Temp;
    Temp = X0; X0 = X1; X1 = Temp;
  }

  if ((DeltaX = X1 - X0) >= 0) {
    XDir = 1;
  }
  else {
    XDir = -1;
    DeltaX = -DeltaX; /* make DeltaX positive */
  }

  /* Draw the initial pixel, which is always exactly intersected by
  the line and so needs no weighting */
  //DrawPixel(pDC, X0, Y0, BaseColor);
  if (DeltaX && (DeltaY = Y1 - Y0))
  startAddr[X0 + Y0*stride] = BaseColor;

  /* Special-case horizontal, vertical, and diagonal lines, which
  require no weighting because they go right through the center of
  every pixel */
  if ((DeltaY = Y1 - Y0) == 0) {
    /* Horizontal line */
    while (DeltaX-- != 0) {
      X0 += XDir;
      if (X0 >= block_width || Y0 >= block_height)
        break;
      //DrawPixel(pDC, X0, Y0, BaseColor);
      startAddr[X0 + Y0*stride] = BaseColor;
    }
    return;
  }
  if (DeltaX == 0) {
    /* Vertical line */
    do {
      Y0++;
      if (X0 >= block_width || Y0 >= block_height)
        break;
      //DrawPixel(pDC, X0, Y0, BaseColor);
      startAddr[X0 + Y0*stride] = BaseColor;
    } while (--DeltaY != 0);
    return;
  }
  if (DeltaX == DeltaY) {
    /* Diagonal line */
    do {
      X0 += XDir;
      Y0++;
      if (X0 >= block_width || Y0 >= block_height)
        break;
      //DrawPixel(pDC, X0, Y0, BaseColor);
      startAddr[X0 + Y0*stride] = BaseColor;
    } while (--DeltaY != 0);
    return;
  }
  /* Line is not horizontal, diagonal, or vertical */
  ErrorAcc = 0;  /* initialize the line error accumulator to 0 */
  /* # of bits by which to shift ErrorAcc to get intensity level */
  IntensityShift = 16 - IntensityBits;
  /* Mask used to flip all bits in an intensity weighting, producing the
  result (1 - intensity weighting) */
  WeightingComplementMask = NumLevels - 1;
  /* Is this an X-major or Y-major line? */
  if (DeltaY > DeltaX) {
    /* Y-major line; calculate 16-bit fixed-point fractional part of a
    pixel that X advances each time Y advances 1 pixel, truncating the
    result so that we won't overrun the endpoint along the X axis */
    ErrorAdj = ((unsigned long)DeltaX << 16) / (unsigned long)DeltaY;
    /* Draw all pixels other than the first and last */
    while (--DeltaY) {
      ErrorAccTemp = ErrorAcc;   /* remember currrent accumulated error */
      ErrorAcc += ErrorAdj;      /* calculate error for next pixel */
      if (ErrorAcc <= ErrorAccTemp) {
        /* The error accumulator turned over, so advance the X coord */
        X0 += XDir;
      }
      Y0++; /* Y-major, so always advance Y */
      /* The IntensityBits most significant bits of ErrorAcc give us the
      intensity weighting for this pixel, and the complement of the
      weighting for the paired pixel */
      Weighting = ErrorAcc >> IntensityShift;
      if (X0 >= block_width || Y0 >= block_height)
        break;
      //DrawPixel(pDC, X0, Y0, BaseColor + Weighting);
      //DrawPixel(pDC, X0 + XDir, Y0,
      //  BaseColor + (Weighting ^ WeightingComplementMask));
      startAddr[X0 + Y0*stride] = BaseColor + Weighting;
      startAddr[X0 + XDir, Y0*stride] = BaseColor + (Weighting ^ WeightingComplementMask);
    }
    /* Draw the final pixel, which is
    always exactly intersected by the line
    and so needs no weighting */
    //DrawPixel(pDC, X1, Y1, BaseColor);
    if (!(X1 >= block_width || Y1 >= block_height))
      startAddr[X1 + Y1*stride] = BaseColor;
    return;
  }
  /* It's an X-major line; calculate 16-bit fixed-point fractional part of a
  pixel that Y advances each time X advances 1 pixel, truncating the
  result to avoid overrunning the endpoint along the X axis */
  ErrorAdj = ((unsigned long)DeltaY << 16) / (unsigned long)DeltaX;
  /* Draw all pixels other than the first and last */
  while (--DeltaX) {
    ErrorAccTemp = ErrorAcc;   /* remember currrent accumulated error */
    ErrorAcc += ErrorAdj;      /* calculate error for next pixel */
    if (ErrorAcc <= ErrorAccTemp) {
      /* The error accumulator turned over, so advance the Y coord */
      Y0++;
    }
    X0 += XDir; /* X-major, so always advance X */
    /* The IntensityBits most significant bits of ErrorAcc give us the
    intensity weighting for this pixel, and the complement of the
    weighting for the paired pixel */
    Weighting = ErrorAcc >> IntensityShift;
    if (X0 >= block_width || Y0 >= block_height)
      break;
    //DrawPixel(pDC, X0, Y0, BaseColor + Weighting);
    //DrawPixel(pDC, X0, Y0 + 1,
    //  BaseColor + (Weighting ^ WeightingComplementMask));
    startAddr[X0 + Y0*stride] = BaseColor + Weighting;
    startAddr[X0, (Y0 + 1)*stride] = BaseColor + BaseColor + (Weighting ^ WeightingComplementMask);
  }
  /* Draw the final pixel, which is always exactly intersected by the line
  and so needs no weighting */
  //DrawPixel(pDC, X1, Y1, BaseColor);
  if (DeltaX && DeltaY)
  if (!(X1 >= block_width || Y1 >= block_height))
    startAddr[X1 + Y1*stride] = BaseColor;
}

void TRACE_CU_PART(TComYuv *picSrc, TComDataCU* pcCU, const char *fname, float scale, int norm)
{
  CHECK_LAYER(pcCU);
  GET_YUV_CU(fname, pcCU);
  TComYuvExt picSrcExt(picSrc);
  picSrcExt.copyToPicPartYuv(picYuv, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu(), scale, norm);
}

void ZERO_CU_PART(TComDataCU* pcCU, UInt uiAbsPartIdx, const char *fname)
{
  CHECK_LAYER(pcCU);
  GET_YUV_CU(fname, pcCU);
  int height = pcCU->getHeight(uiAbsPartIdx);
  int width  = pcCU->getWidth(uiAbsPartIdx);
  Pel *dstY = picYuv->getAddr(COMPONENT_Y,  pcCU->getCtuRsAddr(), uiAbsPartIdx);
  Pel *dstU = picYuv->getAddr(COMPONENT_Cb, pcCU->getCtuRsAddr(), uiAbsPartIdx);
  Pel *dstV = picYuv->getAddr(COMPONENT_Cr, pcCU->getCtuRsAddr(), uiAbsPartIdx);
  for (int y = 0; y < height; y++)
  {
    ::memset(dstY, 0, sizeof(Pel)*width);
    ::memset(dstU, 0, sizeof(Pel)*width>>1);
    ::memset(dstV, 0, sizeof(Pel)*width>>1);
  }
}

// save diff: A - B
void TRACE_CU_PART_DIFF(TComYuv *picSrcA, TComYuv *picSrcB, TComDataCU* pcCU, const char *fname)
{
  CHECK_LAYER(pcCU);
  GET_YUV_CU(fname, pcCU);
  TComYuvExt picSrcExtA(picSrcA);
  TComYuvExt picSrcExtB(picSrcB);
  picSrcExtA.copyToPicPartYuv(picYuv, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
  picSrcExtB.subToPicPartYuv(picYuv, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
}

void TRACE_CU_PART_SUM(TComYuv *picSrcA, TComYuv *picSrcB, TComDataCU* pcCU, const char *fname, Int shift_right)
{
  CHECK_LAYER(pcCU);
  GET_YUV_CU(fname, pcCU);
  TComYuvExt picSrcExtA(picSrcA);
  TComYuvExt picSrcExtB(picSrcB);
  picSrcExtA.copyToPicPartYuv(picYuv, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
  picSrcExtB.sumToPicPartYuv(picYuv, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu(), shift_right);
}

void TRACE_FLOAT(float *src, int width, int height, const char* fname, short bias, float scale)
{
  TComPicYuv* picOut = TRACE_GET_PIC_YUV(fname, width, height);
  Int stride = picOut->getStride(COMPONENT_Y);
  Pel* dst = picOut->getAddr(COMPONENT_Y);

  // luma
  for (int py = 0; py < height; py++)
  {
    for (int px = 0; px < width; px++)
    {
      short val = scale*src[px + py*width] + bias + 0.5;
      dst[px + py*stride] = val;
    }
  }
}

void TRACE_INT(int *src, int width, int height, const char* fname, int shift, int bias, int shift_left)
{
  TComPicYuv* picOut = TRACE_GET_PIC_YUV(fname, width, height);
  Int stride = picOut->getStride(COMPONENT_Y);
  Pel* dst = picOut->getAddr(COMPONENT_Y);

  // luma
  for (int py = 0; py < height; py++)
  {
    for (int px = 0; px < width; px++)
    {
      short val = ((src[px + py*width] >> shift) << shift_left) + bias;
      dst[px + py*stride] = val;
    }
  }
}

void getPicYuv(const char *fname, TComPicYuv**out)
{
  *out = 0;
  TComPicYuv *picYuv = getPicYuv(fname, 0, 0);
  *out = picYuv;
}

TComPicYuv* TRACE_GET_PIC_YUV(const char *fname)
{
  TComPicYuv *out;
  getPicYuv(fname, &out);
  return out;
}

void getPicYuv(const char *fname, TComPicYuv**out, int width, int height)
{
  *out = 0;
  TComPicYuv *picYuv = getPicYuv(fname, width, height);
  *out = picYuv;
}

TComPicYuv* TRACE_GET_PIC_YUV(const char *fname, int width, int height)
{
  TComPicYuv *out;
  getPicYuv(fname, &out, width, height);
  return out;
}

void TRACE_PIC(TComPicYuv *picSrc, const char *fname)
{
  TComPicYuv *picYuv = getPicYuv(fname, picSrc);
  picSrc->copyToPic(picYuv);
}

void markPartitioning(TComDataCU* pcCU, UInt uiAbsPartIdx, TComPicYuvExt *pic)
{
  int width,height;
  Pel *piSrc;
  if (pcCU->getZorderIdxInCtu() == 0)
  {
    // Large CU
    piSrc  = pic->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiAbsPartIdx);
    width  = pcCU->getWidth(uiAbsPartIdx);
    height = pcCU->getHeight(uiAbsPartIdx);
  }
  else
  {
    // Small CU
    piSrc = pic->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
    width  = pcCU->getWidth(0);
    height = pcCU->getHeight(0);
  }
  Int iStride = pic->getStride(COMPONENT_Y);
  memset(piSrc, 0, width*sizeof(Pel));
  for (int i = 0; i < height;i++)
  {
    piSrc[0] = 0;
    piSrc += iStride;
  }
}

#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)


static void put_char(const char ch, Pel *pFrame, int width, int height, int iStride, int &pos_x, int &pos_y, Pel luma_val = 0, 
                     bool limit_with_boundary = true, float scale = 1.f, const unsigned char *bmp_given = 0)
{
/*
  if (ch == '\n')
  {
    pos_y += gFont.Height;
    pos_x  = kStartX;
    return;
  }
*/
  int idx = ch - FONT_START_INDEX;
  const unsigned char *bmp = bmp_given ? bmp_given : &gFont.Bitmap[idx*gFont.Height];
  {
    Pel *yuv_data = pFrame + pos_y*iStride;
    for (int i = 0; i < (limit_with_boundary ? min(gFont.Height, height - pos_y) : gFont.Height); i++)
    {
      for (int k = 0; k < (limit_with_boundary ? min(gFont.Width, width - pos_x) : gFont.Width); k++)
      {
        if (bmp[i] & (1 << (gFont.Width - k - 1)))
          yuv_data[(int)((pos_x + k)*scale + 0.5)] = luma_val;        
      }
      yuv_data += (int)(iStride*scale +0.5);
    }
  }

  pos_x += gFont.Width >> 1; // slim fonts
}


static void put_str(const char *ch, Pel *pFrame, int width, int height, int iStride, int pos_x, int pos_y, Pel luma_val = 0, 
                       bool limit_with_boundary = true, float scale = 1.f, const unsigned char *bmp_given = 0)
{
  int x_start = pos_x;
  int x = pos_x;
  int y = pos_y;
  while (*ch != 0)
  {
    if (*ch == '\n')
    {
      x  = x_start;
      y += gFont.Height >> 1; // slim font
    }
    else
      put_char(*ch, pFrame, width, height, iStride, x, y, luma_val, limit_with_boundary, scale, bmp_given);

    ch++;
  }
}

static const unsigned char arrow_bmp_u[] = {
	_X______,
	X_X_____,
        ________,
        ________,
        ________,
        ________,
        ________,
        ________,
        ________,
        ________,
        ________,
        ________,
	________,
	________,
	________,
	________,
};

static const unsigned char arrow_bmp_d[] = {
  X_X_____,
  _X______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};
static const unsigned char arrow_bmp_r[] = {
  X_______,
  _X______,
  X_______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};
static const unsigned char arrow_bmp_l[] = {
  _X______,
  X_______,
  _X______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};
static const unsigned char arrow_bmp_lu[] = {
  XX______,
  X_______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};
static const unsigned char arrow_bmp_ru[] = {
  XX______,
  _X______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};

static const unsigned char arrow_bmp_ld[] = {
  X_______,
  XX______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};
static const unsigned char arrow_bmp_rd[] = {
  _X______,
  XX______,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
  ________,
};

const unsigned char* get_arrow_bmp_by_degree360(float deg)
{
  assert(deg >= 0 && deg <= 360.f);
  if (deg >= 0)
  {
    if (deg < 22.5)       return arrow_bmp_r;
    else if (deg < 67.5)  return arrow_bmp_ru;
    else if (deg < 112.5) return arrow_bmp_u;
    else if (deg < 157.5) return arrow_bmp_lu;
    else if (deg < 202.5) return arrow_bmp_l;
    else if (deg < 247.5) return arrow_bmp_ld;
    else if (deg < 292.5) return arrow_bmp_d;
    else if (deg < 337.5) return arrow_bmp_rd;
    else                  return arrow_bmp_r;
  }
  else
  {
    assert(0);
    return arrow_bmp_r;
  }
}

void markLCU(TComDataCU* pcCU, TComPicYuvExt *pic)
{
  Pel* piSrc = pic->getLumaAddr(pcCU->getCtuRsAddr(), 2, 2);
  piSrc[0] = 0;
}

void TEXT_CU(TComDataCU* pcCU, UInt uiAbsPartIdx, int width, int height, const char* str, const char *fname, int drawCUborder, int lineColor)
{
  GET_YUV_CU(fname, pcCU);

  Pel *piSrc = picYuv->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiAbsPartIdx);
  // draw borders
  if (drawCUborder)
  {
    for (int i = 0; i < width; i++)
      piSrc[i] = lineColor;
    for (int i = 0; i < height;i++)
      piSrc[i*picYuv->getStride(COMPONENT_Y)] = lineColor;
  }

  put_str(str, piSrc, width, height, picYuv->getStride(COMPONENT_Y), 2, 1, lineColor);
}


/** calculate the location of upper-left corner pixel and size of the current PU in relative to upper-left corner of CU at index uiAbsPartIdx.
* \param partIdx  PU index within a CU
* \param xP, yP   location of the upper-left corner pixel of the current PU
* \param PSW, nPSH    size of the current PU
* \returns Void
*/
Void getPartPositionMy(TComDataCU *pcCU, UInt uiAbsPartIdx, UInt partIdx, Int& xP, Int& yP, Int& nPSW, Int& nPSH)
{
  UInt col = 0;//g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt row = 0;//g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
#if !JVET_C0024_QTBT
  switch (pcCU->getPartitionSize(uiAbsPartIdx))
  {
  case SIZE_2NxN:
    nPSW = pcCU->getWidth(uiAbsPartIdx);
    nPSH = pcCU->getHeight(uiAbsPartIdx) >> 1;
    xP = col;
    yP = (partIdx == 0) ? row : row + nPSH;
    break;
  case SIZE_Nx2N:
    nPSW = pcCU->getWidth(uiAbsPartIdx) >> 1;
    nPSH = pcCU->getHeight(uiAbsPartIdx);
    xP = (partIdx == 0) ? col : col + nPSW;
    yP = row;
    break;
  case SIZE_NxN:
    nPSW = pcCU->getWidth(uiAbsPartIdx) >> 1;
    nPSH = pcCU->getHeight(uiAbsPartIdx) >> 1;
    xP = col + (partIdx & 0x1)*nPSW;
    yP = row + (partIdx >> 1)*nPSH;
    break;
  case SIZE_2NxnU:
    nPSW = pcCU->getWidth(uiAbsPartIdx);
    nPSH = (partIdx == 0) ? pcCU->getHeight(uiAbsPartIdx) >> 2 : (pcCU->getHeight(uiAbsPartIdx) >> 2) + (pcCU->getHeight(uiAbsPartIdx) >> 1);
    xP = col;
    yP = (partIdx == 0) ? row : row + pcCU->getHeight(uiAbsPartIdx) - nPSH;

    break;
  case SIZE_2NxnD:
    nPSW = pcCU->getWidth(uiAbsPartIdx);
    nPSH = (partIdx == 0) ? (pcCU->getHeight(uiAbsPartIdx) >> 2) + (pcCU->getHeight(uiAbsPartIdx) >> 1) : pcCU->getHeight(uiAbsPartIdx) >> 2;
    xP = col;
    yP = (partIdx == 0) ? row : row + pcCU->getHeight(uiAbsPartIdx) - nPSH;
    break;
  case SIZE_nLx2N:
    nPSW = (partIdx == 0) ? pcCU->getWidth(uiAbsPartIdx) >> 2 : (pcCU->getWidth(uiAbsPartIdx) >> 2) + (pcCU->getWidth(uiAbsPartIdx) >> 1);
    nPSH = pcCU->getHeight(uiAbsPartIdx);
    xP = (partIdx == 0) ? col : col + pcCU->getWidth(uiAbsPartIdx) - nPSW;
    yP = row;
    break;
  case SIZE_nRx2N:
    nPSW = (partIdx == 0) ? (pcCU->getWidth(uiAbsPartIdx) >> 2) + (pcCU->getWidth(uiAbsPartIdx) >> 1) : pcCU->getWidth(uiAbsPartIdx) >> 2;
    nPSH = pcCU->getHeight(uiAbsPartIdx);
    xP = (partIdx == 0) ? col : col + pcCU->getWidth(uiAbsPartIdx) - nPSW;
    yP = row;
    break;
  default:
    //assert ( getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N );
    nPSW = pcCU->getWidth(uiAbsPartIdx);
    nPSH = pcCU->getHeight(uiAbsPartIdx);
    xP = col;
    yP = row;

    break;
  }
#endif
}

Bool isIdenticalMotionMy(TComDataCU* pcCU, UInt PartAddr)
{
  if (pcCU->getSlice()->isInterB() && !pcCU->getSlice()->getPPS()->getWPBiPred())
  {
    if (pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(PartAddr) >= 0 && pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(PartAddr) >= 0)
    {
      Int RefPOCL0 = pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(PartAddr))->getPOC();
      Int RefPOCL1 = pcCU->getSlice()->getRefPic(REF_PIC_LIST_1, pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(PartAddr))->getPOC();
      if (RefPOCL0 == RefPOCL1 && pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(PartAddr) == pcCU->getCUMvField(REF_PIC_LIST_1)->getMv(PartAddr))
      {
        return true;
      }
    }
  }
  return false;
}

void markPU_inter(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, TComPicYuvExt *pic, bool bPutChar = false, int lineColor = 0, bool bPutCost = false, int refList = 1)
{
#if !JVET_C0024_QTBT
  PartSize ePartSize = pcCU->getPartitionSize(uiAbsPartIdx);
  UInt uiNumPU = ( ePartSize == SIZE_2Nx2N ? 1 : ( ePartSize == SIZE_NxN ? 4 : 2 ) );
  UInt uiDepth = pcCU->getDepth( uiAbsPartIdx );
  UInt uiPUOffset = (g_auiPUOffset[UInt(ePartSize)] << ((pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth) << 1)) >> 4;
#else
  UInt uiNumPU = 1;
  UInt uiPUOffset = 0;
#endif
  for ( UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset )
  {
    Int x,y, width, height;
    x = y = 0;

#if !JVET_C0024_QTBT
    getPartPositionMy(pcCU, uiAbsPartIdx, uiPartIdx, x, y, width, height);
#else
    width  = uiWidth;
    height = uiHeight;
#endif
    x += pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
    y += pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];

    TComPic *pRefPic  = 0;
    TComPic *pRefPic0 = 0;
    TComPic *pRefPic1 = 0;

    // draw borders
    Pel *piSrc = pic->getLumaAddrXY(x, y);
    for (int i = 0; i < width; i++)
      piSrc[i] = lineColor;
    for (int i = 0; i < height;i++)
      piSrc[i*pic->getStride(COMPONENT_Y)] = lineColor;

    if (pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiSubPartIdx) >= 0)
      pRefPic0 = pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiSubPartIdx));
    if (pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiSubPartIdx) >= 0)
      pRefPic1 = pcCU->getSlice()->getRefPic(REF_PIC_LIST_1, pcCU->getCUMvField(REF_PIC_LIST_1)->getRefIdx(uiSubPartIdx));

    Int biPrediction = (pcCU->getCUMvField( REF_PIC_LIST_0 )->getRefIdx( uiSubPartIdx ) >= 0) && 
                       (pcCU->getCUMvField( REF_PIC_LIST_1 )->getRefIdx( uiSubPartIdx ) >= 0);
    // check identical motion
    if (biPrediction && isIdenticalMotionMy(pcCU, uiSubPartIdx))
        biPrediction = false;

    RefPicList eRefPicList = (RefPicList)refList;

    // print mode
    if (bPutChar)
    {
      char str[255];
      Pel luma_val = lineColor;
      sprintf(str, "");

      TComMv mv = pcCU->getCUMvField(eRefPicList)->getMv(uiSubPartIdx);
      Int ref_idx = pcCU->getCUMvField(eRefPicList)->getRefIdx(uiSubPartIdx);
      Int ref_poc = -1;
      if (ref_idx >= 0)
      {
        ref_poc = pcCU->getSlice()->getRefPic(eRefPicList, ref_idx)->getPOC();
        pRefPic = pcCU->getSlice()->getRefPic(eRefPicList, ref_idx);
      }
      Int refQP = 88;
      /*
      if (ref_idx >= 0)
      {
        TComPic *refPic = pcCU->getSlice()->getRefPic(eRefPicList, ref_idx);
        refQP = refPic->getQP();
      }
      */
#if RRC_SELECTABLE_IF
      //int sharpMode = pcCU->getSIFMode(&mv, uiSubPartIdx, 1, refQP);
      int sharpMode = pcCU->getSIFFlag(uiSubPartIdx);
      if (sharpMode > 0)
      {
        luma_val = 0x88;
        if (sharpMode == 1)
          luma_val = 0xff;
        luma_val <<= gBitDepthsTrace.recon[CHANNEL_TYPE_LUMA] - 8;
        {
          sprintf(str, "%d", sharpMode);
        }
      }
#endif
#if 0// B-Pred
      if (biPrediction)
      {
        luma_val = 0x88;
        sprintf(str, "%s%d", str, ref_poc);
      }
#endif
      if (luma_val != lineColor)
      {
        // draw borders again with changed color
        Pel *piSrc = pic->getLumaAddrXY(x, y);
        for (int i = 0; i < width; i++) piSrc[i] = luma_val;
        for (int i = 0; i < height; i++) piSrc[i*pic->getStride(COMPONENT_Y)] = luma_val;
      }

#if 0//SVC_SHARPENING_BLOCK_MODE || SVC_SHARPENING_BLOCK_AFTER_ALL || SVC_SHARPENING_INTER
      TSharpnessBlockParams *par = pcCU->getSIFFlag(uiSubPartIdx);
      // get ref idx
      Int ref_idx = pcCU->getCUMvField(eRefPicList)->getRefIdx(uiPartIdx);
      if (ref_idx >= 0 && ref_idx < MAX_NUM_REF)
      {
        assert(ref_idx < MAX_NUM_REF);
        sprintf(str, "%d", ref_idx);

       if (pcCU->getSharpEnable(uiSubPartIdx))
       {
          luma_val = 0x80;
          {
            sprintf(str, "%d", par->k_idx);
          }
          assert(par->k_idx != -1);
          // draw borders again with changed color
          Pel *piSrc = pic->getLumaAddrXY(x, y);
          for (int i = 0; i < width; i++) piSrc[i] = luma_val;
          for (int i = 0; i < height;i++) piSrc[i*pic->getStride()] = luma_val;
        }
      }
#endif

#if 0// skip flag
      if (pcCU->isSkipped(uiSubPartIdx))
        sprintf(str, "%s\n", "+");

      //sprintf(str, "%s\n%d", str, pcCU->isSkipped(uiSubPartIdx));
      /*
        if (pcCU->isSkipped(uiSubPartIdx))
      {
        sprintf(str, "%s\n", "s");
      }
   
      else
      {
        if (pcCU->getROI(uiAbsPartIdx) >= 0) 
          pic->setPicPartCr(pcCU->getROI(uiAbsPartIdx), pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiAbsPartIdx, width, height);
        else
          pic->setPicPartCr(250, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiAbsPartIdx, width, height);
      }
      */
#endif

#if 0// transform skip
      sprintf(str, "%s\n%d", str, pcCU->getTransformSkip(uiSubPartIdx, TEXT_LUMA));
      //if (pcCU->getTransformSkip(uiSubPartIdx, TEXT_LUMA))
        //sprintf(str, "%s\n", "+");
#endif

      // k_idx
#if 0//SVC_SHARPENING_BLOCK_MODE || SVC_SHARPENING_BLOCK_AFTER_ALL || SVC_SHARPENING_INTER      
      if (par->is_enabled)
        sprintf(str, "%d", par->k_idx);
      else
        sprintf(str, "");
#endif

#if 0 // edge dir,abs
      //sprintf(str, "%s\n%d",str, (int)(par->edge_dir+0.5));
      if (par->is_enabled)
        sprintf(str, "%d", par->k_idx);
      else
        sprintf(str, "");
      sprintf(str, "%s,%d", str, par->mv_class);
      //sprintf(str, "%s\n%d,%d", str,par->edge_dir_idx,par->edge_abs_idx);
      //sprintf(str, "%s\n%d,%d",str, par->edge_dir_idx, (int)(par->edge_dir+0.5));
#endif

#if 0 // merge & amvp
      sprintf(str, "%s\n%d",str, pcCU->getMergeFlag(uiSubPartIdx));
      if (pcCU->getMergeFlag(uiSubPartIdx))
        sprintf(str, "%s,%d",str, pcCU->getMergeIndex(uiSubPartIdx));
      else
      {
        sprintf(str, "%s,%d",str, pcCU->getMVPIdx(eRefPicList, uiSubPartIdx));
        TComMv mvd = pcCU->getCUMvField(eRefPicList)->getMvd(uiSubPartIdx);
        sprintf(str, "%s\n%d,%d",str, mvd.getHor(),mvd.getVer());
      }
#endif

#if 1// merge
      if (pcCU->getMergeFlag(uiSubPartIdx))
      {
        //pic->setPicPartCb(0, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
        sprintf(str, "%s%s", str, "+");
      }
#endif

#if 1 && VCEG_AZ07_IMV
      if (pcCU->getiMVFlag(uiSubPartIdx))
      {
        sprintf(str, "%s\n%s", str, "imv");
      }
#endif

#if 1 && COM16_C1016_AFFINE
      if (pcCU->isAffine(uiSubPartIdx))
      {
        sprintf(str, "%s%s", str, "aff");
      }
#endif

#if 0 // split
      sprintf(str, "%s%d", str, pcCU->getDepth(uiSubPartIdx));
#endif

#if 1 && VCEG_AZ07_FRUC_MERGE
      if (pcCU->getFRUCMgrMode(uiSubPartIdx))
      {
        sprintf(str, "%s\n%s%d", str, "fr", pcCU->getFRUCMgrMode(uiSubPartIdx));
        //markPU_interFRUC(pcCU, uiSubPartIdx, width, height, eRefPicList, pic);
      }
#endif

#if 1 && COM16_C806_OBMC
      if (pcCU->getOBMCFlag(uiSubPartIdx))
      {
        sprintf(str, "%s\n%s", str, "ob");
      }
#endif

#if 1 && VCEG_AZ06_IC
      if (pcCU->getICFlag(uiSubPartIdx))
      {
        sprintf(str, "%s\n%s", str, "ic");
      }
#endif

      //sprintf(str, "%s\n", str);

#if RRC_SELECTABLE_IF // sharp enabled, possible
      {
        //pic->setPicPartCb(0, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
        Int fltTypeHor, fltTypeVer;
        Int mode = pcCU->getSIFModeMv(&mv, fltTypeHor, fltTypeVer, uiSubPartIdx);
        //Int mode = pcCU->getSIFMode(&mv, uiSubPartIdx, 0, 0);
        sprintf(str, "%s\n%d\n%d,%d", str, mode, fltTypeHor, fltTypeVer);
        //sprintf(str, "%d,%d,%d\n", mode, fltTypeHor, fltTypeVer);
        //sprintf(str, "%d\n", mode);
      }
#endif

#if 1 // mv
      sprintf(str, "%s%d,%d,%d",str, mv.getHor() & 0x3, mv.getVer() & 0x3, pcCU->getCUMvField(eRefPicList)->getRefIdx(uiSubPartIdx));
#endif


#if 0 // mv frac
      //sprintf(str, "%s\n%d,%d",str, mv.getHor() & 3, mv.getVer() & 3);
      sprintf(str, "%s %d,%d",str, mv.getHor() & 3, mv.getVer() & 3);
#endif

#if 0 && SEGMENT_CLASSES
      sprintf(str, "%s\n%d",str, pcCU->getClassId(uiSubPartIdx));
#endif

#if 1 && MOTION_MASK // motion mask
    if (gMotionMaskPic)
    {
      Float m = gMotionMaskPic->getAvgPartLuma(pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);

      //if (m > MV_MASK_THRE) pic->setPicPartCb(0, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
      //pic->setPicPartCb(m, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);

      //if (mv.m_MvMask > 0)
      {
        //pic->setPicPartCr(mv.m_MvMask, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
      }

      //assert(pcCU->getMask(uiSubPartIdx) != -1); // not initialized

      if (pcCU->getMask(uiSubPartIdx) > 0) pic->setPicPartCb(0, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);

      // AMVP
      //if (pcCU->getROI(uiSubPartIdx) >= 0) pic->setPicPartCr(pcCU->getROI(uiSubPartIdx), pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
    }
#endif

#if 0 // IF mode merge
    if (pcCU->getIFModeMerge() == 1)
    {
      pic->setPicPartCr(0, pcCU->getAddr(), pcCU->getZorderIdxInCU() + uiSubPartIdx, width, height);
    }
#endif


#if 0 // edge count
      sprintf(str, "%s\n%d",str, par->edge_cnt);
#endif

#if 0 // edge_abs
      //sprintf(str, "%s\n%.2f",str, par->edge_abs);
      sprintf(str, "%s\n%d",str, par->edge_abs_idx);
      //sprintf(str, "%s\n%d,%.2f",str, par->edge_abs_idx, par->edge_abs);
#endif

#if 0 // edge_type
      sprintf(str, "%s\n%d",str, par->edge_type);
#endif

#if 0 // edge_dir
      //sprintf(str, "%s\n%d",str, (int)(par->edge_dir+0.5));
      sprintf(str, "%s\n%d",str, par->edge_dir_idx);
      //sprintf(str, "%s\n%d,%d",str, par->edge_dir_idx, (int)(par->edge_dir+0.5));
#endif

#if BIG && RRC_SELECTABLE_IF
      if (pRefPic/* && biPrediction*/)
      {
        tFrameStat stat = {0};
        pRefPic->getSharpeningFilter()->getBlockStat(&stat, pcCU, uiSubPartIdx, width, height);

        TSharpnessBlockParams params = { 0 };
        pRefPic->getSharpeningFilter()->FillBlockParams(pcCU, uiSubPartIdx, width, height, &params, pRefPic->getPicYuvRec(), &mv);
        sprintf(str, "%s\n%.0f", str, params.direction);
        sprintf(str, "%s\n%.1f, %.1f", str, params.dx_avg, params.dy_avg);
        //line_QPel(piSrc, 2 + width / 2, 2 + height / 2, params.dx_avg, params.dy_avg, pic->getStride(COMPONENT_Y), 0x00, width, height);
        DrawWuLine(piSrc, 2 + width / 2, 2 + height / 2, -params.dy_avg, params.dx_avg, 0x00, 255, 8, pic->getStride(COMPONENT_Y), width, height, 10);

        if (0) // gradient field
        {
          for (int y = 0; y < width; y+=3)
          {
            for (int x = 0; x < height; x+=3)
            {
              int dx = pRefPic->getSharpeningFilter()->_dxPic.getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiSubPartIdx)
                [x + y*pRefPic->getSharpeningFilter()->_dxPic.getStride(COMPONENT_Y)];
              int dy = pRefPic->getSharpeningFilter()->_dyPic.getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiSubPartIdx)
                [x + y*pRefPic->getSharpeningFilter()->_dyPic.getStride(COMPONENT_Y)];
              dx /= 10;
              dy /= 10;
              //line_QPel(piSrc, x, y, dx, dy, pic->getStride(COMPONENT_Y), 0x80, width, height);
              //AALine(x, y, dx, dy, piSrc, pic->getStride(COMPONENT_Y), 0x80, width, height);
              DrawWuLine(piSrc, x, y, dx, dy, 0x00, 255, 8, pic->getStride(COMPONENT_Y), width, height);
            }
          }
         }
      }
#endif

      int start_x = 2; int start_y = 3;
      put_str(str, piSrc, width, height, pic->getStride(COMPONENT_Y), start_x, start_y, luma_val);
#if 0 // motion vector line
      int stride = pic->getStride();
      line_QPel(piSrc, width/2, height/2, mv.getHor(), mv.getVer(), pic->getStride(), 0x80, width, height);
#endif

#if 0
      if (biPrediction)
      {
        TComMv mv0 = pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(uiSubPartIdx);
        TComMv mv1 = pcCU->getCUMvField(REF_PIC_LIST_1)->getMv(uiSubPartIdx);
        int stride = pic->getStride(COMPONENT_Y);
        mv0 -= mv1;
        line_QPel(piSrc, width / 2, height / 2, mv0.getHor(), mv0.getVer(), pic->getStride(COMPONENT_Y), 0x00, width, height);
      }
#endif

#if 0//SVC_SHARPENING_BLOCK_MODE || SVC_SHARPENING_BLOCK_AFTER_ALL || SVC_SHARPENING_INTER      //direction
      start_x = -1 + width >> 1; start_y = -1 + height >> 1;
      put_char(0, piSrc, width, height, pic->getStride(), start_x, start_y, 0, true, 1, get_arrow_bmp_by_degree360(par->edge_dir));
#endif
    }
  }
}

#if !JVET_C0024_QTBT
void markPU_intra(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, TComPicYuvExt *pic, bool bPutChar = true, int lineColor = 0, bool bPutCost = false)
{
  int width,height,mode_idx,cu_idx;
  Pel *piSrc;
  PartSize partSize;
  if (pcCU->getZorderIdxInCtu() == 0)
  {
    // Large CU
    piSrc     = pic->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), uiAbsPartIdx);
    width     = pcCU->getWidth(uiAbsPartIdx);
    height    = pcCU->getHeight(uiAbsPartIdx);
    partSize  = pcCU->getPartitionSize( uiAbsPartIdx );
    mode_idx  = pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
    cu_idx    = uiAbsPartIdx;
  }
  else
  {
    // Small CU
    piSrc     = pic->getAddr(COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu());
    width     = pcCU->getWidth(0);
    height    = pcCU->getHeight(0);
    partSize  = pcCU->getPartitionSize(0);
    mode_idx  = pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiAbsPartIdx);
    cu_idx    = pcCU->getZorderIdxInCtu();
  }
  UInt partNum = partSize == SIZE_NxN ? 4 : 1;
  UInt partOffset = (pcCU->getPic()->getNumPartitionsInCtu() >> (pcCU->getDepth(uiAbsPartIdx) << 1)) >> 2;
  if (partNum > 1)
  {
    width  >>= 1;
    height >>= 1;
  }
  int bitShift = gBitDepthsTrace.recon[COMPONENT_Y] - 8;

  for (UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < partNum; uiPartIdx++, uiSubPartIdx += partOffset)
  {
    Int x, y, width2, height2;

    getPartPositionMy(pcCU, uiSubPartIdx, uiPartIdx, x, y, width2, height2);
    x += pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
    y += pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];

    Pel *piSrc = pic->getLumaAddrXY(x, y);
    Int iStride = pic->getStride(COMPONENT_Y);

    Pel *piSrcCr = pic->getAddrXY(x, y, COMPONENT_Cr);
    Pel *piSrcCb = pic->getAddrXY(x, y, COMPONENT_Cb);
    int chromaStride = pic->getStride(COMPONENT_Cr);

#if RRC_INTRA_SHARPENING
    Int sharpMode = pcCU->getIntraSharpFlag(uiSubPartIdx);
    if (sharpMode > 0)
    {
      lineColor = 0x80;
      lineColor <<= bitShift;
    }
#endif

    int lineColorTop  = lineColor;
    int lineColorLeft = lineColor;
    if (uiPartIdx == 1 || uiPartIdx == 3)
      lineColorLeft = 127;
    if (uiPartIdx == 2 || uiPartIdx == 3)
      lineColorTop = 127;
    for (int i = 0; i < width; i++)
      piSrc[i] = lineColorTop;
    for (int i = 0; i < height; i++)
      piSrc[i*iStride] = lineColorLeft;

    if (bPutChar)
    {
      int start_x = 0;
      int start_y = 2;
      char str[255] = "";
      Pel luma_val = lineColor;
      //sprintf(str, "%s","i");

#if 1 // predMode
      sprintf(str, "%s%d", str, pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiSubPartIdx));
#endif

#if COM16_C1046_PDPC_INTRA
      sprintf(str, "%s\n%d", str, pcCU->getPDPCIdx(uiSubPartIdx));
#endif

#if 0 // CBF luma
      if (!pcCU->getCbf(uiSubPartIdx, COMPONENT_Y, pcCU->getTransformIdx(uiSubPartIdx)))
      for (int y = 0; y < height >> pic->getComponentScaleY(COMPONENT_Cb); y++)
      {
        for (int x = 0; x < width >> pic->getComponentScaleX(COMPONENT_Cb); x++)
        {
          //piSrcCr[y*chromaStride + x] = Clip3(0, 255, avgX + avgY);
          piSrcCb[y*chromaStride + x] = 100;
        }
      }      
#endif
      put_str(str, piSrc, width, height, pic->getStride(COMPONENT_Y), start_x, start_y, luma_val);
    }
  }

  // TODO: don't work for encoder
/*
  UInt  uiInitTrDepth = ( pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N ? 0 : 1 );
  UInt partOffset = ( pcCU->getPic()->getNumPartInCU() >> ( pcCU->getDepth(uiAbsPartIdx) << 1 ) ) >> 2;

  PartSize mode = pcCU->getPartitionSize( uiAbsPartIdx );
  UInt partNum = mode==SIZE_NxN?4:1;

  for( UInt uiPU = 0; uiPU < partNum; uiPU++ )
  {
    UInt    width           = pcCU     ->getWidth   ( 0 ) >> uiInitTrDepth;
    UInt    height          = pcCU     ->getHeight  ( 0 ) >> uiInitTrDepth;
    int addr = uiAbsPartIdx + uiPU * partNum;
    UInt x = g_auiRasterToPelX[ g_auiZscanToRaster[addr] ];
    UInt y = g_auiRasterToPelY[ g_auiZscanToRaster[addr] ];
    Pel *piSrc = pic->getLumaAddr(pcCU->getAddr(), addr);
    Int iStride = pic->getStride();
    memset(piSrc, 0, width*sizeof(Pel));
    for (int i = 0; i < height;i++)
    {
      piSrc[0] = 0;
      piSrc += iStride;
    }
    if (bPutChar)
    {
      int start_x = 0;
      int start_y = 1;
      put_char('i', pic->getLumaAddr(pcCU->getAddr(), addr), 
        width, height, pic->getStride(), start_x, start_y, 1.f, 0);
    }
  }  
*/
}
#else
void markPU_intraQTBT(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, TComPicYuvExt *pic, bool bPutChar = true, int lineColor = 0, bool bPutCost = false)
{
#if !JVET_C0024_QTBT
  PartSize ePartSize = pcCU->getPartitionSize(uiAbsPartIdx);
  UInt uiNumPU = (ePartSize == SIZE_2Nx2N ? 1 : (ePartSize == SIZE_NxN ? 4 : 2));
  UInt uiDepth = pcCU->getDepth(uiAbsPartIdx);
  UInt uiPUOffset = (g_auiPUOffset[UInt(ePartSize)] << ((pcCU->getSlice()->getSPS()->getMaxTotalCUDepth() - uiDepth) << 1)) >> 4;
#else
  UInt uiNumPU = 1;
  UInt uiPUOffset = 0;
#endif
  int bitShift = gBitDepthsTrace.recon[COMPONENT_Y] - 8;

  for (UInt uiPartIdx = 0, uiSubPartIdx = uiAbsPartIdx; uiPartIdx < uiNumPU; uiPartIdx++, uiSubPartIdx += uiPUOffset)
  {
    Int x, y, width, height;
    x = y = 0;

#if !JVET_C0024_QTBT
    getPartPositionMy(pcCU, uiAbsPartIdx, uiPartIdx, x, y, width, height);
#else
    width = uiWidth;
    height = uiHeight;
#endif
    x += pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
    y += pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];

    Pel *piSrc = pic->getLumaAddrXY(x, y);
    Int iStride = pic->getStride(COMPONENT_Y);

    Pel *piSrcCr = pic->getAddrXY(x, y, COMPONENT_Cr);
    Pel *piSrcCb = pic->getAddrXY(x, y, COMPONENT_Cb);
    int chromaStride = pic->getStride(COMPONENT_Cr);

#if RRC_INTRA_SHARPENING
    Int sharpMode = pcCU->getIntraSharpFlag(uiSubPartIdx);
    if (sharpMode > 0)
    {
      lineColor = 0xff << bitShift;
    }
#endif

    int lineColorTop = lineColor;
    int lineColorLeft = lineColor;
    if (uiPartIdx == 1 || uiPartIdx == 3)
      lineColorLeft = 127;
    if (uiPartIdx == 2 || uiPartIdx == 3)
      lineColorTop = 127;
    for (int i = 0; i < width; i++)
      piSrc[i] = lineColorTop;
    for (int i = 0; i < height; i++)
      piSrc[i*iStride] = lineColorLeft;

    if (bPutChar)
    {
      char str[255] = "";
      Pel luma_val = lineColor;
      sprintf(str, "%s","i");

#if 1 // predMode
      sprintf(str, "%s%d", str, pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiSubPartIdx));
#endif

#if 0 && COM16_C1046_PDPC_INTRA
      sprintf(str, "%s\n%d", str, pcCU->getPDPCIdx(uiSubPartIdx));
#endif

#if 0 // CBF luma
      if (!pcCU->getCbf(uiSubPartIdx, COMPONENT_Y, pcCU->getTransformIdx(uiSubPartIdx)))
        for (int y = 0; y < height >> pic->getComponentScaleY(COMPONENT_Cb); y++)
        {
          for (int x = 0; x < width >> pic->getComponentScaleX(COMPONENT_Cb); x++)
          {
            //piSrcCr[y*chromaStride + x] = Clip3(0, 255, avgX + avgY);
            piSrcCb[y*chromaStride + x] = 100;
          }
        }
#endif
      int start_x = 2; int start_y = 3;
      put_str(str, piSrc, width, height, pic->getStride(COMPONENT_Y), start_x, start_y, luma_val);
    }
  }
}
#endif
Void xIntraRecBlkMark(TComYuv*    pcRecoYuv,
                 TComYuv*    pcPredYuv,
                 TComYuv*    pcResiYuv,
                 const ComponentID compID,
                 TComTU     &rTu
                 , TComPicYuvExt *pic, bool bPutChar = true, int lineColor = 0, bool bPutCost = false)
{
  if (!rTu.ProcessComponentSection(compID))
  {
    return;
  }
  const Bool       bIsLuma = isLuma(compID);


  TComDataCU *pcCU = rTu.getCU();
  const TComSPS &sps = *(pcCU->getSlice()->getSPS());
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();

  const TComRectangle &tuRect = rTu.getRect(compID);
  const UInt uiWidth = tuRect.width;
  const UInt uiHeight = tuRect.height;
//  const UInt uiStride = pcRecoYuv->getStride(compID);
//  Pel* piPred = pcPredYuv->getAddr(compID, uiAbsPartIdx);
  const ChromaFormat chFmt = rTu.GetChromaFormat();

  if (uiWidth != uiHeight)
  {
    //------------------------------------------------

    //split at current level if dividing into square sub-TUs

    TComTURecurse subTURecurse(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

    //recurse further
    do
    {
      xIntraRecBlkMark(pcRecoYuv, pcPredYuv, pcResiYuv, compID, subTURecurse, pic, bPutChar, lineColor, bPutCost);
    } while (subTURecurse.nextSection(rTu));

    //------------------------------------------------

    return;
  }
}
#if !JVET_C0024_QTBT
Void xIntraRecQTMark(TComYuv*    pcRecoYuv,
                    TComYuv*    pcPredYuv,
                    TComYuv*    pcResiYuv,
                    const ChannelType chType,
                    TComTU     &rTu
                    , TComPicYuvExt *pic, bool bPutChar = true, int lineColor = 0, bool bPutCost = false)
{
  UInt uiTrDepth    = rTu.GetTransformDepthRel();
  TComDataCU *pcCU  = rTu.getCU();
  UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if( uiTrMode == uiTrDepth )
  {
    if (isLuma(chType))
    {
      xIntraRecBlkMark(pcRecoYuv, pcPredYuv, pcResiYuv, COMPONENT_Y, rTu, pic, bPutChar, lineColor, bPutCost);
    }
    else
    {
      const UInt numValidComp=getNumberValidComponents(rTu.GetChromaFormat());
      for(UInt compID=COMPONENT_Cb; compID<numValidComp; compID++)
      {
        xIntraRecBlkMark(pcRecoYuv, pcPredYuv, pcResiYuv, ComponentID(compID), rTu, pic, bPutChar, lineColor, bPutCost);
      }
    }
  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xIntraRecQTMark(pcRecoYuv, pcPredYuv, pcResiYuv, chType, tuRecurseChild, pic, bPutChar, lineColor, bPutCost);
    } while (tuRecurseChild.nextSection(rTu));
  }

  Int x, y, width2, height2;
  Int width  = rTu.getRect(COMPONENT_Y).width;
  Int height = rTu.getRect(COMPONENT_Y).height;

  UInt uiSubPartIdx = uiAbsPartIdx;
  UInt uiPartIdx = 0;
  getPartPositionMy(pcCU, uiSubPartIdx, uiPartIdx, x, y, width2, height2);
  x += pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  y += pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];

  Pel *piSrc = pic->getLumaAddrXY(x, y);
  Int iStride = pic->getStride(COMPONENT_Y);

  Pel *piSrcCr = pic->getAddrXY(x, y, COMPONENT_Cr);
  Pel *piSrcCb = pic->getAddrXY(x, y, COMPONENT_Cb);
  int chromaStride = pic->getStride(COMPONENT_Cr);

  int bitShift = gBitDepthsTrace.recon[COMPONENT_Y] - 8;

#if RRC_INTRA_SHARPENING
  Int sharpMode = pcCU->getIntraSharpFlag(uiSubPartIdx);
  if (sharpMode > 0)
  {
    lineColor = 0x80;
  }
#endif
  int lineColorTop = lineColor;
  int lineColorLeft = lineColor;
  PartSize ePartSize = pcCU->getPartitionSize(uiAbsPartIdx);
  int sectionNum = rTu.GetSectionNumber();
  if (uiTrDepth == 1 && ePartSize == SIZE_NxN)
  {
    if (sectionNum == 1 || sectionNum == 3)
      lineColorLeft = 127;
    if (sectionNum == 2 || sectionNum == 3)
      lineColorTop  = 127;
  }
  else if (uiTrDepth >= 1)
  {
    if (sectionNum == 1 || sectionNum == 3)
      lineColorLeft = 180;
    if (sectionNum == 2 || sectionNum == 3)
      lineColorTop  = 180;
  }
  //PartSize ePartSize = pcCU->getPartitionSize(uiAbsPartIdx);
  for (int i = 0; i < width; i++)
    piSrc[i] = lineColorTop << bitShift;
  for (int i = 0; i < height; i++)
    piSrc[i*iStride] = lineColorLeft << bitShift;

  if (bPutChar)
  {
    int start_x = 0;
    int start_y = 2;
    char str[255] = "";

#if COM16_C806_EMT
    if (pcCU->getEmtCuFlag(uiSubPartIdx))
      sprintf(str, "%s%d", str, pcCU->getEmtCuFlag(uiSubPartIdx));
#endif
#if 0 && COM16_C1046_PDPC_INTRA
    if (pcCU->getPDPCIdx(uiSubPartIdx))
    sprintf(str, "%s%d", str, pcCU->getPDPCIdx(uiSubPartIdx));
#endif
#if 0 && (VCEG_AZ05_ROT_TR || COM16_C1044_NSST)
    if (pcCU->getROTIdx(uiSubPartIdx))
    sprintf(str, "%s\n%d", str, pcCU->getROTIdx(uiSubPartIdx));
#endif
    
    put_str(str, piSrc, width, height, pic->getStride(COMPONENT_Y), start_x, start_y, lineColor);
  }

}

void xMarkIntraQT(TComDataCU* pcCU, UInt uiDepth, UInt uiWidth, UInt uiHeight, TComPicYuvExt *pic, bool bPutChar = true, int lineColor = 0, bool bPutCost = false)
{
  const UInt numChType = 1;// pcCU->getPic()->getChromaFormat() != CHROMA_400 ? 2 : 1;
  for (UInt chType = CHANNEL_TYPE_LUMA; chType<numChType; chType++)
  {
    const ChannelType chanType = ChannelType(chType);
    const Bool NxNPUHas4Parts = ::isChroma(chanType) ? enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) : true;
    const UInt uiInitTrDepth = (pcCU->getPartitionSize(0) != SIZE_2Nx2N && NxNPUHas4Parts ? 1 : 0);

    TComTURecurse tuRecurseCU(pcCU, 0);
    TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

    do
    {
      xIntraRecQTMark(NULL, NULL, NULL, chanType, tuRecurseWithPU, pic, bPutChar, lineColor, bPutCost);
    } while (tuRecurseWithPU.nextSection(tuRecurseCU));
  }
}
#endif

void MARK_PARTIONING(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiWidth, UInt uiHeight, const char *fname, int lineColor, bool putChar, int refList)
{
  GET_YUV_CU(fname, pcCU);
  assert(picYuv->getWidth(COMPONENT_Y) && picYuv->getHeight(COMPONENT_Y));
  markLCU(pcCU, picYuv);
#if !JVET_C0024_QTBT
  switch(pcCU->getPredictionMode(uiAbsPartIdx))
#else
  if (!gIsEncoder)
    uiAbsPartIdx = 0;
  switch (pcCU->getPredictionMode(uiAbsPartIdx))
#endif
  {
//*
  case MODE_INTER:
    markPU_inter(pcCU, uiAbsPartIdx, uiWidth, uiHeight, picYuv, putChar, lineColor, 1, refList);
    break;
  case MODE_INTRA:
#if JVET_C0024_QTBT
    markPU_intraQTBT(pcCU, uiAbsPartIdx, uiWidth, uiHeight, picYuv, putChar, 0x80 << (gBitDepthsTrace.recon[COMPONENT_Y] - 8)/*lineColor*/, 1);
#else
if (gIsEncoder)
    markPU_intra(pcCU, uiAbsPartIdx, uiWidth, uiHeight, picYuv, putChar, lineColor, 1);
else
    xMarkIntraQT(pcCU, uiAbsPartIdx, uiWidth, uiHeight, picYuv, putChar, lineColor, 1);
#endif
    break;
  default:
    markPartitioning(pcCU, uiAbsPartIdx, picYuv);
    break;
  }
}

#if 0//!JVET_C0024_QTBT
Void xMarkCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, const char *fname, int lineColor, bool putChar, int refList)
{
  TComPic   *const pcPic = pcCU->getPic();
  TComSlice *const pcSlice = pcCU->getSlice();
  const TComSPS   &sps = *(pcSlice->getSPS());
  const TComPPS   &pps = *(pcSlice->getPPS());

#if !JVET_C0024_QTBT
  const UInt maxCUWidth = sps.getMaxCUWidth();
  const UInt maxCUHeight = sps.getMaxCUHeight();
#else
  const UInt maxCUWidth  = sps.getCTUSize();
  const UInt maxCUHeight = sps.getCTUSize();
#endif
  Bool bBoundary = false;
  UInt uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  const UInt uiRPelX = uiLPelX + (maxCUWidth >> uiDepth) - 1;
  UInt uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
  const UInt uiBPelY = uiTPelY + (maxCUHeight >> uiDepth) - 1;

  if ((uiRPelX < sps.getPicWidthInLumaSamples()) && (uiBPelY < sps.getPicHeightInLumaSamples()))
  {
  }
  else
  {
    bBoundary = true;
  }

  if (((uiDepth < pcCU->getDepth(uiAbsPartIdx)) && (uiDepth < sps.getLog2DiffMaxMinCodingBlockSize())) || bBoundary)
  {
    UInt uiQNumParts = (pcPic->getNumPartitionsInCtu() >> (uiDepth << 1)) >> 2;

    for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx += uiQNumParts)
    {
      uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
      uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
      if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
      {
        xMarkCU(pcCU, uiAbsPartIdx, uiDepth + 1, fname, lineColor, putChar, refList);
      }
    }
    return;
  }
  MARK_PARTIONING(pcCU, uiAbsPartIdx, fname, lineColor, putChar, refList);
}

void MARK_PARTIONING_LCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt startDepth, const char *fname, int lineColor, bool putChar, int refList)
{
  xMarkCU(pcCU, uiAbsPartIdx, startDepth, fname, lineColor, putChar, refList);
}
#endif

template <typename TPel> void TRACE_REF_SAMPLES(TPel *pSrc, Int iWidth, Int iHeight, Int iSrcStride, TComDataCU* pcCU, UInt uiAbsPartIdx, Int iVerHor, const char *fname, float scale, int norm, int bias)
{
  CHECK_LAYER(pcCU);
  GET_YUV_CU(fname, pcCU);

  ComponentID compID = COMPONENT_Y;

  Pel* pDst = picYuv->getAddr(compID, pcCU->getCtuRsAddr(), uiAbsPartIdx + pcCU->getZorderIdxInCtu());
  Int offset = 0;
  Int shiftNum = 0;
  if (norm)
  {
    shiftNum = IF_INTERNAL_PREC /*+ 1*/ - 8;
    offset = (1 << (shiftNum - 1)) +/* 2 * */IF_INTERNAL_OFFS;
  }

  iWidth = (iWidth - 1) >> 1;
  iHeight = (iHeight - 1) >> 1;

  UInt  iDstStride = picYuv->getStride(compID);
  if (iVerHor == 0)
  {
    // vertical
    for (int y = iHeight; y != 0; y--)
    {
      *pDst = *pSrc*scale + bias;
      if (norm)
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
        *pDst = ClipA((*pDst + offset) >> shiftNum, compID);
#else
        *pDst = ClipBD((*pDst + offset) >> shiftNum, 8);
#endif
      }
      pDst += iDstStride;
      pSrc += iSrcStride;
    }
  }
  else
  {
    for (int i = 0; i < iWidth; i++)
    {
      pDst[i] = (Pel)(pSrc[i] * scale + bias);
      if (norm)
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
        pDst[i] = ClipA((pDst[i] + offset) >> shiftNum, compID);
#else
        pDst[i] = ClipBD((pDst[i] + offset) >> shiftNum, 8);
#endif
      }
    }
  }
}

template void TRACE_REF_SAMPLES<Pel>(Pel *pSrc, Int iWidth, Int iHeight, Int iSrcStride, TComDataCU* pcCU, UInt uiAbsPartIdx, Int iVerHor, const char *fname, float scale, int norm, int bias);
template void TRACE_REF_SAMPLES<float>(float *pSrc, Int iWidth, Int iHeight, Int iSrcStride, TComDataCU* pcCU, UInt uiAbsPartIdx, Int iVerHor, const char *fname, float scale, int norm, int bias);
void DUMP_ALL_YUV()
{  
  for (int i = 0; i < STREAMS_NUM; i++)
  {
    if (!gpStreams[i])
      break;

    assert(gpPicYuv[i]);
//    gpPicYuv[i]->dump(gpStreamNames[i], {8,8}, true);
    gpPicYuv[i]->dump(gpStreamNames[i], gBitDepthsTrace, true);
    gpPicYuv[i]->clean();
  }
}



#endif //TRACE_YUV
