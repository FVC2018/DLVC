#include "TComSharpeningFilter.h"
#include "TComPrediction.h"
#include <math.h>
#include <fstream>

#pragma warning( disable : 4244 ) //  conversion from 'int' to 'float', possible loss of data

#if RRC_SELECTABLE_IF

#define DEF_SATURATE             (10.f)
#define DEF_FLOOR                (0.1f)
#define MAP_SHIFT                 4

#define SATURATE_FIX   (10 * 3.2 + 0.5)


static float filterDx[] = {
  -1, 0, 1,
  -1, 0, 1,
  -1, 0, 1,
};

static float filterDy[] = {
  1,  1,  1,
  0,  0,  0,
 -1, -1, -1,
};

#define NG(x) x/(273.f/25.f)
float gauss_blur_filter[] = {
  NG(1), NG(4), NG(7), NG(4), NG(1),
  NG(4), NG(16), NG(26), NG(16), NG(4),
  NG(7), NG(26), NG(41), NG(26), NG(7),
  NG(4), NG(16), NG(26), NG(16), NG(4),
  NG(1), NG(4), NG(7), NG(4), NG(1),
};

#define CHECK_BOUND(x,width)          \
  if (x < 0) x = 0;                   \
        else if(x >= width) x= width - 1;

#define COUNT(a)    (sizeof(a)/sizeof(a[0]))

//////////////////////////////////////////////////////////////////////////////////////
/////////////////////Free-standing functions definitions//////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
static void filter_kernel(float* src, float* dst, int stride, const BlockDescription& block, float* filter, int filter_half_len)
{
  src += block.get_offset();
  dst += block.get_offset();

  int filter_width = filter_half_len * 2 + 1;

  for (int py = 0; py < block.get_height(); py++)
  {
    for (int px = 0; px < block.get_width(); px++)
    {
      float res = 0;
      for (int fy = -filter_half_len; fy <= filter_half_len; fy++)
      {
        int y = py + fy;

        block.check_y_border(&y);
        for (int fx = -filter_half_len; fx <= filter_half_len; fx++)
        {
          int x = px + fx;
          block.check_x_border(&x);
          res += src[x + y*stride] * filter[fx + filter_half_len + (fy + filter_half_len)*filter_width];
        }
      }

      dst[px + py*stride] = res/ (filter_width*filter_width);
    }
  }
}

static void filter_kernel_no_check(float* src, float* dst, int stride, const BlockDescription& block, float* filter, int filter_half_len)
{
  src += block.get_offset();
  dst += block.get_offset();

  int filter_width = filter_half_len * 2 + 1;

  for (int py = 0; py < block.get_height(); py++)
  {
    for (int px = 0; px < block.get_width(); px++)
    {
      float res = 0;
      for (int fy = -filter_half_len; fy <= filter_half_len; fy++)
      {
        int y = py + fy;
        for (int fx = -filter_half_len; fx <= filter_half_len; fx++)
        {
          int x = px + fx;
          res += src[x + y*stride] * filter[fx + filter_half_len + (fy + filter_half_len)*filter_width];
        }
      }

      dst[px + py*stride] = res / (filter_width*filter_width);
    }
  }
}

static void prewitt_dx(float* src, float* dst, int stride, const BlockDescription& block)
{
  src += block.get_offset();
  dst += block.get_offset();
  for (int py = 0; py < block.get_height(); py++)
    for (int px = 0; px < block.get_width(); px++)
      dst[px + py * stride] =
      (src[px + 1 + (py - 1)*stride] - src[px - 1 + (py - 1) * stride] +
      src[px + 1 + (py)* stride] - src[px - 1 + (py)* stride] +
      src[px + 1 + (py + 1)*stride] - src[px - 1 + (py + 1) * stride]) / (3 * 3);
}
static void prewitt_dy(float* src, float* dst, int stride, const BlockDescription& block)
{
  src += block.get_offset();
  dst += block.get_offset();
  for (int py = 0; py < block.get_height(); py++)
    for (int px = 0; px < block.get_width(); px++)
      dst[px + py*stride] =
      (src[px - 1 + (py - 1) * stride] - src[px - 1 + (py + 1) * stride] +
      src[px + (py - 1) * stride] - src[px + (py + 1) * stride] +
      src[px + 1 + (py - 1) * stride] - src[px + 1 + (py + 1) * stride]) / (3 * 3);
}
static void floor_vector(float* src, float c, int width, int height, int stride)
{
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j)
      if (fabs(src[i * stride + j]) <= c)
        src[i * stride + j] = 0;
}

static void saturate(float *src, float c, int width, int height, int stride)
{
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j)
      if (src[i * stride + j] > 0 && src[i * stride + j] > c)
        src[i * stride + j] = c;
      else if (src[i * stride + j] < 0 && src[i * stride + j] < -c)
        src[i * stride + j] = -c;
}

static void der_x(float *d, float *dd, int stride, const BlockDescription& block)
{
  for (int py = 0; py < block.get_height(); py++)
  {
    for (int px = 0; px < block.get_width(); px++)
    {
      int x1 = px - 1;
      block.check_x_border(&x1);
      int x2 = px + 1;
      block.check_x_border(&x2);
      dd[px + py*stride] = d[x1 + py*stride] - d[x2 + py*stride];
    }
  }
}

static void der_y(float *d, float *dd, int stride, const BlockDescription& block)
{
  for (int py = 0; py < block.get_height(); py++)
  {
    int y1 = py - 1;
    block.check_y_border(&y1);
    int y2 = py + 1;
    block.check_y_border(&y2);
    for (int px = 0; px < block.get_width(); px++)
    {
      dd[px + py*stride] = d[px + y1*stride] - d[px + y2*stride];
    }
  }
}

static void der_x_no_check(float *d, float *dd, int stride, const BlockDescription& block)
{
  for (int py = 0; py < block.get_height(); py++)
  {
    for (int px = 0; px < block.get_width(); px++)
    {
      int x1 = px - 1;
      int x2 = px + 1;
      dd[px + py*stride] = d[x1 + py*stride] - d[x2 + py*stride];
    }
  }
}

static void der_y_no_check(float *d, float *dd, int stride, const BlockDescription& block)
{
  for (int py = 0; py < block.get_height(); py++)
  {
    int y1 = py - 1;
    int y2 = py + 1;
    for (int px = 0; px < block.get_width(); px++)
    {
      dd[px + py*stride] = d[px + y1*stride] - d[px + y2*stride];
    }
  }
}

static void displaceBlk(TComYuv* pcPicSrc, TComYuv* pcPicDst, TComYuv *pcMapX, TComYuv *pcMapY,
  int map_shift, BitDepths &bitDepths, ComponentID compId, float k, int x_start, int y_start, int block_width, int block_height, int clipOut, int margin = 0)
{
  Int widthSrc  = pcPicSrc->getWidth (compId);
  Int heightSrc = pcPicSrc->getHeight(compId);
  Int strideSrc = pcPicSrc->getStride(compId);
  Pel* src      = pcPicSrc->getAddr  (compId);

  Int widthDst  = pcPicDst->getWidth (compId);
  Int heightDst = pcPicDst->getHeight(compId);
  Int strideDst = pcPicDst->getStride(compId);
  Pel* dst      = pcPicDst->getAddr  (compId);

  Int widthMap  = pcMapX->getWidth (compId);
  Int heightMap = pcMapX->getHeight(compId);
  Int strideMap = pcMapX->getStride(compId);
  
  ComponentID mapCompId = isChroma(compId) ? COMPONENT_Cb : compId; // take map from Cb for chroma
  Pel* mapX     = pcMapX->getAddr(mapCompId);
  Pel* mapY     = pcMapY->getAddr(mapCompId);

  int bitDepth = bitDepths.recon[toChannelType(compId)];
  int maxVal = (1 << bitDepth) - 1;
  int midVal = 1 << (bitDepth - 1);

  if (x_start < 0) x_start = 0;
  if (x_start + block_width > widthSrc)
    block_width = widthSrc - x_start;
  if (y_start < 0) y_start = 0;
  if (y_start + block_height > heightSrc)
    block_height = heightSrc - y_start;

  map_shift += bitDepth - 8;
  float scale = 1.f / (1 << map_shift);
  for (int py = y_start; py < y_start + block_height; py++)
  {
    for (int px = x_start; px < x_start + block_width; px++)
    {
      float x = px;
      float y = py;
      float dx = (mapX[px + py*strideMap] - midVal) * scale;
      float dy = (mapY[px + py*strideMap] - midVal) * scale;
      x += k*dx;
      y += k*dy;

      // block borders
      if (x < x_start - margin) x = x_start - margin;
      if (x > x_start + block_width - 1 + margin)  x = x_start + block_width - 1 + margin;
      if (y < y_start - margin) y = y_start - margin;
      if (y > y_start + block_height - 1 + margin) y = y_start + block_height - 1 + margin;

      x += margin;
      y += margin;

      // bilinear interpolation
      if (1)
      {
        int base_x, base_y;
        float dx, dy, val;
        base_x = (int)x;
        base_y = (int)y;
        dx = x - base_x;
        dy = y - base_y;
        {
          val = ((1 - dx)*(1 - dy)*src[base_x +      base_y     *strideSrc] +
                 (1 - dx)*(dy)    *src[base_x +     (base_y + 1)*strideSrc] +
                     (dx)*(1 - dy)*src[base_x + 1 +  base_y     *strideSrc] +
                     (dx)*(dy)    *src[base_x + 1 + (base_y + 1)*strideSrc]);
        }
        Pel pelVal = (Pel)(val + 0.5);
        if (clipOut)
        {
          //pelVal = (pelVal < 0) ? 0 : pelVal;
          //pelVal = (pelVal > maxVal) ? maxVal : pelVal;
#if JVET_D0033_ADAPTIVE_CLIPPING // encoder intrac rec
          pelVal = Pel(ClipA<Pel>(pelVal, compId));
#else
          pelVal = Pel(ClipBD<Pel>(pelVal, bitDepth));
#endif
        }
        dst[px + py*strideDst] = pelVal;
      }
    }
  }
}

static void downsample_luma2chroma(float *src_luma, float* dst_chroma, int stride, const BlockDescription& block_location)
{
  //int width_src = block_location.get_width();
  //int height_src = block_location.get_height();
  int stride_src = stride << 1;

  for (int py = 0; py < block_location.get_height() / 2; py++)
  {
    for (int px = 0; px < block_location.get_width() / 2; px++)
    {
      float sum = 0;
      for (int fy = 2 * py; fy <= 2 * py + 1; fy++)
      {
        int fx;
        fx = 2 * px - 1;
        block_location.check_x_border(&fx);
        sum += src_luma[fx + fy*stride_src] * 1 / 8.f;
        fx = 2 * px;
        sum += src_luma[fx + fy*stride_src] * 2 / 8.f;
        fx = 2 * px + 1;
        block_location.check_x_border(&fx);
        sum += src_luma[fx + fy*stride_src] * 1 / 8.f;
      }
      dst_chroma[px + py*stride] = sum; // TODO: rounding?
    }
  }
}

static void find_abs(float *dx, float *dy, float *dxy_abs, int stride, const BlockDescription& block)
{
  dx += block.get_offset();
  dy += block.get_offset();
  dxy_abs += block.get_offset();

  for (int i = 0; i < block.get_height(); ++i)
    for (int j = 0; j < block.get_width(); ++j)
    {
      dxy_abs[i * stride + j] = abs(dx[i*stride + j]) + abs(dy[i*stride + j]);
    }
}

static void extend_yuv(const BlockDescription& block, TComYuv* d2xYuv, TComYuv* d2yYuv)
{
  //top extension
  for (int i = 0; i < block.get_top_extension(); ++i)
  {
    Pel* src_ptr_x = d2xYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + block.get_top_extension() * d2xYuv->getStride(COMPONENT_Y);
    Pel* src_ptr_y = d2yYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + block.get_top_extension() * d2yYuv->getStride(COMPONENT_Y);
    Pel* dst_ptr_x = d2xYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + i * d2xYuv->getStride(COMPONENT_Y);
    Pel* dst_ptr_y = d2yYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + i * d2xYuv->getStride(COMPONENT_Y);

    copy(src_ptr_x, src_ptr_x + block.get_width(), dst_ptr_x);
    copy(src_ptr_y, src_ptr_y + block.get_width(), dst_ptr_y);
  }

  //bottom extension
  for (int i = block.get_height(); i < block.get_height() + block.get_bottom_extension(); ++i)
  {
    Pel* src_ptr_x = d2xYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + (block.get_height() - 1) * d2xYuv->getStride(COMPONENT_Y);
    Pel* src_ptr_y = d2yYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + (block.get_height() - 1) * d2xYuv->getStride(COMPONENT_Y);
    Pel* dst_ptr_x = d2xYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + i * d2xYuv->getStride(COMPONENT_Y);
    Pel* dst_ptr_y = d2yYuv->getAddr(COMPONENT_Y) + block.get_left_extension() + i * d2xYuv->getStride(COMPONENT_Y);

    copy(src_ptr_x, src_ptr_x + block.get_width(), dst_ptr_x);
    copy(src_ptr_y, src_ptr_y + block.get_width(), dst_ptr_y);
  }

  //left extension
  for (int i = 0; i < block.get_height() + block.get_top_extension() + block.get_bottom_extension(); ++i)
  {
    Pel* dst_ptr_x    = d2xYuv->getAddr(COMPONENT_Y) + i * d2xYuv->getStride(COMPONENT_Y);
    Pel* dst_ptr_y    = d2yYuv->getAddr(COMPONENT_Y) + i * d2yYuv->getStride(COMPONENT_Y);
    Pel fill_val_it_x = d2xYuv->getAddr(COMPONENT_Y)[i * d2xYuv->getStride(COMPONENT_Y) + block.get_left_extension()];
    Pel fill_val_it_y = d2yYuv->getAddr(COMPONENT_Y)[i * d2yYuv->getStride(COMPONENT_Y) + block.get_left_extension()];

    fill(dst_ptr_x, dst_ptr_x + block.get_left_extension(), fill_val_it_x);
    fill(dst_ptr_y, dst_ptr_y + block.get_left_extension(), fill_val_it_y);
  }

  //right extension
  for (int i = 0; i < block.get_height() + block.get_top_extension() + block.get_bottom_extension(); ++i)
  {
    Pel* dst_ptr_x    = d2xYuv->getAddr(COMPONENT_Y) + i * d2xYuv->getStride(COMPONENT_Y) + block.get_width();
    Pel* dst_ptr_y    = d2yYuv->getAddr(COMPONENT_Y) + i * d2xYuv->getStride(COMPONENT_Y) + block.get_width();
    Pel fill_val_it_x = d2xYuv->getAddr(COMPONENT_Y)[i * d2xYuv->getStride(COMPONENT_Y) + block.get_width() - 1];
    Pel fill_val_it_y = d2yYuv->getAddr(COMPONENT_Y)[i * d2yYuv->getStride(COMPONENT_Y) + block.get_width() - 1];

    fill(dst_ptr_x, dst_ptr_x + block.get_right_extension(), fill_val_it_x);
    fill(dst_ptr_y, dst_ptr_y + block.get_right_extension(), fill_val_it_y);
  }
}

//////////////////////////////////////////////////////////////////////////////////////
/////////////////////class TComSharpeningFilter Definitions////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
const int TComSharpeningFilter::_interp_filter_margin(4);
const int TComSharpeningFilter::_internal_extension((int)sqrt(COUNT(filterDx)) / 2 + (int)sqrt(COUNT(gauss_blur_filter)) / 2);

TComSharpeningFilter::TComSharpeningFilter(TComPicYuv* pcBasePic, const TComSPS* sps, bool isEncoder)
  : _d2x(0)
  , _d2y(0)
  , _dx(0)
  , _dy(0)
  , _d2x_chroma(0)
  , _d2y_chroma(0)
  , _isEncoder(isEncoder)
{
  const int height = pcBasePic->getHeight(COMPONENT_Y);
  const int width  = pcBasePic->getWidth(COMPONENT_Y);
  m_chromaFormatIDC = pcBasePic->getChromaFormat();

  m_uiMaxCUWidth  = MAX_CU_SIZE;
  m_uiMaxCUHeight = MAX_CU_SIZE;
  m_uiMaxCUDepth  = sps->getMaxTotalCUDepth();
  m_bitDepths     = sps->getBitDepths();

  _picYuv.create   (MAX_CU_SIZE, MAX_CU_SIZE, m_chromaFormatIDC);
  _picYuvTmp.create(MAX_CU_SIZE, MAX_CU_SIZE, m_chromaFormatIDC);

  const int max_yuv_block_size = MAX_CU_SIZE + 2 * _interp_filter_margin + _internal_extension;
  _d2xYuv.create(max_yuv_block_size, max_yuv_block_size, m_chromaFormatIDC);
  _d2yYuv.create(max_yuv_block_size, max_yuv_block_size, m_chromaFormatIDC);

  _d2xYuvInter.create(MAX_CU_SIZE, MAX_CU_SIZE, m_chromaFormatIDC);
  _d2yYuvInter.create(MAX_CU_SIZE, MAX_CU_SIZE, m_chromaFormatIDC);

  _src = new float[width * height];
  _tmp = new float[width * height];
  _d2x = new float[width * height];
  _d2y = new float[width * height];
  _dx = new float[width * height];
  _dy = new float[width * height];
  _d2x_chroma = new float[(width / 2) * (height / 2)];
  _d2y_chroma = new float[(width / 2) * (height / 2)];
  _dxy_abs = new float[width * height];

  for (int i = 0; i < pcBasePic->getHeight(COMPONENT_Y); ++i)
    copy(pcBasePic->getAddr(COMPONENT_Y) + i * pcBasePic->getStride(COMPONENT_Y), pcBasePic->getAddr(COMPONENT_Y) + i * pcBasePic->getStride(COMPONENT_Y) + pcBasePic->getWidth(COMPONENT_Y), _src + i * width);
}

TComSharpeningFilter::~TComSharpeningFilter()
{
  _d2xPic.destroy();
  _d2yPic.destroy();
  _picYuv.destroy();
  _picYuvTmp.destroy();
  _d2xYuv.destroy();
  _d2yYuv.destroy();
  _d2xYuvInter.destroy();
  _d2yYuvInter.destroy();

  delete[] _src;
  delete[] _tmp;
  delete[] _d2x;
  delete[] _d2y;
  delete[] _dx;
  delete[] _dy;
  delete[] _d2x_chroma;
  delete[] _d2y_chroma;
  delete[] _dxy_abs;
}

void TComSharpeningFilter::clear_buffers()
{
  const int height = _d2xPic.getHeight(COMPONENT_Y);
  const int width  = _d2xPic.getWidth(COMPONENT_Y);

  fill(_d2x, _d2x + width * height, 0.f);
  fill(_d2y, _d2y + width * height, 0.f);
  fill(_dx, _dx + width * height, 0.f);
  fill(_dy, _dy + width * height, 0.f);
  fill(_dxy_abs, _dxy_abs + width * height, 0.f);
  fill(_d2x_chroma, _d2x_chroma + width * height / 4, 0.f);
  fill(_d2y_chroma, _d2y_chroma + width * height / 4, 0.f);
}

Void TComSharpeningFilter::buildBlockEdgeMap(TComPicYuv* pcBasePic, TComDataCU* cu, UInt partAddr, TComMv* mv, UInt uiWidth, UInt uiHeight)
{
  assert(0); // TODO: fix map shift for 10 bits
  const int height = pcBasePic->getHeight(COMPONENT_Y);
  const int width  = pcBasePic->getWidth(COMPONENT_Y);

  BlockDescription block(cu, partAddr, mv, uiWidth, uiHeight, width, height);
  BlockDescription block_ = block.create_extended(_interp_filter_margin);
  BlockDescription block_ex = block_.create_extended(_internal_extension, true);

  m_is_left_chroma_alignment = block.get_left_chroma_alignment();
  m_is_top_chroma_alignment = block.get_top_chroma_alignment();

  if (block_ex.is_on_frame_border())
  {
    filter_kernel(_src, _dx, width, block_ex, filterDx, sqrt((float)COUNT(filterDx)) / 2);
    filter_kernel(_src, _dy, width, block_ex, filterDy, sqrt((float)COUNT(filterDy)) / 2);
  }
  else
  {
    prewitt_dx(_src, _dx, width, block_ex);
    prewitt_dy(_src, _dy, width, block_ex);
  }

  find_abs(_dx, _dy, _dxy_abs, width, block_ex);
  saturate(_dxy_abs + block_ex.get_offset(), DEF_SATURATE, block_ex.get_width(), block_ex.get_height(), width);
  floor_vector(_dxy_abs + block_ex.get_offset(), DEF_FLOOR, block_ex.get_width(), block_ex.get_height(), width);

  if (block_ex.is_on_frame_border())
  {
    filter_kernel(_dxy_abs, _tmp, width, block_.create_extended(1, true), gauss_blur_filter, sqrt((float)COUNT(gauss_blur_filter)) / 2);
    der_x(_tmp + block_.get_offset(), _d2x + block_.get_offset(), width, block_);
    der_y(_tmp + block_.get_offset(), _d2y + block_.get_offset(), width, block_);
  }
  else
  {
    filter_kernel_no_check(_dxy_abs, _tmp, width, block_.create_extended(1, true), gauss_blur_filter, sqrt((float)COUNT(gauss_blur_filter)) / 2);
    der_x_no_check(_tmp + block_.get_offset(), _d2x + block_.get_offset(), width, block_);
    der_y_no_check(_tmp + block_.get_offset(), _d2y + block_.get_offset(), width, block_);
  }

  int i = 0;
  do
  {
    int j = 0;
    do
    {
      const int dst_index = (j + block_.get_left_extension()) + (i + block_.get_top_extension()) * _d2xYuv.getWidth(COMPONENT_Y);
      const int src_index = block_.get_offset() + j + i * width;

      _d2xYuv.getAddr(COMPONENT_Y)[dst_index] = 127 + _d2x[src_index] * (1 << MAP_SHIFT) + 0.5;
      _d2yYuv.getAddr(COMPONENT_Y)[dst_index] = 127 + _d2y[src_index] * (1 << MAP_SHIFT) + 0.5;
    } while (j++ < block_.get_width() - block_.get_right_border_displacement() - 1);
  } while (i++ < block_.get_height() - block_.get_bottom_border_displacement() - 1);

  const unsigned int chroma_offset = ((block_.get_offset() / width) / 2) * (width / 2) + (block_.get_offset() % width) / 2;

  downsample_luma2chroma(_d2x + block_.get_offset(), _d2x_chroma + chroma_offset, width >> 1, block_);
  downsample_luma2chroma(_d2y + block_.get_offset(), _d2y_chroma + chroma_offset, width >> 1, block_);

  for (int i = 0; i < block_.get_height() / 2; ++i)
    for (int j = 0; j < block_.get_width() / 2; ++j)
    {
      const int dst_index = (j + block_.get_left_extension() / 2) + (i + block_.get_top_extension() / 2) * _d2xYuv.getWidth(COMPONENT_Cb);

      _d2xYuv.getAddr(COMPONENT_Cb)[dst_index] = 127 + _d2x_chroma[chroma_offset + j + i * width / 2] * (1 << MAP_SHIFT) + 0.5;
      _d2xYuv.getAddr(COMPONENT_Cr)[dst_index] = 127 + _d2x_chroma[chroma_offset + j + i * width / 2] * (1 << MAP_SHIFT) + 0.5;

      _d2yYuv.getAddr(COMPONENT_Cb)[dst_index] = 127 + _d2y_chroma[chroma_offset + j + i * width / 2] * (1 << MAP_SHIFT) + 0.5;
      _d2yYuv.getAddr(COMPONENT_Cr)[dst_index] = 127 + _d2y_chroma[chroma_offset + j + i * width / 2] * (1 << MAP_SHIFT) + 0.5;
    }

  extend_yuv(block_, &_d2xYuv, &_d2yYuv);
}

static Void initMatrix_int(int ***m2D, int width, int height)
{
  int i;

  if (!(*m2D = (int **)calloc(height, sizeof(int *))))
    FATAL_ERROR_0("initMatrix_int: memory allocation problem\n", -1);
  if (!((*m2D)[0] = (int *)calloc(width * height, sizeof(int))))
    FATAL_ERROR_0("initMatrix_int: memory allocation problem\n", -1);

  for (i = 1; i < height; i++)
    (*m2D)[i] = (*m2D)[i - 1] + width;
}

static Void destroyMatrix_int(int **m2D)
{
  if (m2D)
  {
    if (m2D[0])
      free(m2D[0]);
    else
      FATAL_ERROR_0("destroyMatrix_int: memory free problem\n", -1);
    free(m2D);
  }
}


#define FLT_5       0 
#define FLT_6       1

#if 1
#define FLT_TYPE    FLT_5 // [1,4,7,4,1]
#define FLT_OFFSET  2
#define DER_OFFSET  1
#else
#define FLT_TYPE    FLT_6 // [1,2,3,3,2,1]
#define FLT_OFFSET  3
#define DER_OFFSET  1
#endif

#define OFFSET      (DER_OFFSET + FLT_OFFSET + DER_OFFSET)
#define OFFSET_2    (2*OFFSET)

void der_abs(Pel *srcPic, Int srcStride, Int *dstPic, Int dstStride, Int width, Int height, Int offset = 0, Int bitDepth = 8)
{
  for (Int y = 0 - offset; y < height + offset; y++)
  {
    Pel *src     = srcPic + y*srcStride;
    Pel *srcUp   = srcPic + (y + 1)*srcStride;
    Pel *srcDown = srcPic + (y - 1)*srcStride;
    Int *dst     = dstPic + y*dstStride;
    for (Int x = 0 - offset; x < width + offset; x++)
    {
      Pel ver = abs(srcDown[x] - srcUp[x]);
      Pel hor = abs(src[x + 1] - src[x - 1]);
      dst[x] = ver + hor;
    }
  }
}

void apply_limits(Int *srcPic, Int srcStride, Int *dstPic, Int dstStride, Int saturation_limit, Int width, Int height, Int offset = 0)
{
  for (Int y = 0 - offset; y < height + offset; y++)
  {
    Int *src = srcPic + y*srcStride;
    Int *dst = dstPic + y*dstStride;
    for (Int x = 0 - offset; x < width + offset; x++)
    {
      dst[x] = src[x] < saturation_limit ? src[x] : saturation_limit;
    }
  }
}

void flt_block_x(Int *srcPic, Int srcStride, Int *dstPic, Int dstStride, Int width, Int height, Int offset = 0)
{
  for (Int y = 0 - offset; y < height + offset; y++)
  {
    Int *src = srcPic + y*srcStride;
    Int *dst = dstPic + y*dstStride;
    for (Int x = 0 - offset; x < width + offset; x++)
    {
      dst[x] 
#if FLT_TYPE  == FLT_6
        = ((src[x + 0] + src[x + 5])
        + ((src[x + 1] + src[x + 4]) << 1)
        + ((src[x + 2] + src[x + 3]) * 3));
#elif FLT_TYPE == FLT_5
        = ((src[x - 2] + src[x + 2])
        + ((src[x - 1] + src[x + 1]) << 2)
        + ((src[x + 0]) << 3) - src[x + 0]);
#else
        assert(0)
#endif
    }
  }
}

void flt_block_y(Int *srcPic, Int srcStride, Int *dstPic, Int dstStride, Int width, Int height, Int offset = 0)
{
  for (Int y = 0 - offset; y < height + offset; y++)
  {
    Int *src = srcPic + y*srcStride;
    Int *dst = dstPic + y*dstStride;
    for (Int x = 0 - offset; x < width + offset; x++)
    {
      dst[x]
#if FLT_TYPE  == FLT_6
        = ((src[x + 0 * srcStride] + src[x + 5 * srcStride])
        + ((src[x + 1 * srcStride] + src[x + 4 * srcStride]) << 1)
        + ((src[x + 2 * srcStride] + src[x + 3 * srcStride]) * 3));
#elif FLT_TYPE == FLT_5
        = ((src[x - 2 * srcStride] + src[x + 2 * srcStride])
        + ((src[x - 1 * srcStride] + src[x + 1 * srcStride]) << 2)
        + ((src[x + 0 * srcStride]) << 3) - src[x + 0 * srcStride])
        ;
#else
        assert(0)
#endif
    }
  }
}

void der_xy(Int *srcPic, Int srcStride, Int *dstPicX, Int dstStrideX, Int *dstPicY, Int dstStrideY, Int width, Int height)
{
  for (Int y = 0; y < height; y++)
  {
    Int *src     = srcPic + y*srcStride;
    Int *srcUp   = srcPic + (y + 1)*srcStride;
    Int *srcDown = srcPic + (y - 1)*srcStride;
    Int *dst_x = dstPicX + y*dstStrideX;
    Int *dst_y = dstPicY + y*dstStrideY;
    for (Int x = 0; x < width; x++)
    {
      Int ver = srcDown[x] - srcUp[x];
      Int hor = src[x + 1] - src[x - 1];
      dst_x[x] = - hor;
      dst_y[x] =   ver;
    }
  }
}

void arr_toPic(Int *srcPic, Int srcStride, Pel *dstPic, Int dstStride, Int width, Int height, int shift, int bitDepth)
{
  shift += bitDepth - 8;
  int   mid = 1 << (bitDepth - 1);
  for (Int y = 0; y < height; y++)
  {
    Int *src = srcPic + y*srcStride;
    Pel *dst = dstPic + y*dstStride;
    for (Int x = 0; x < width; x++)
    {
      dst[x] = mid + (src[x] >> shift);
    }
  }
}
void downsample_luma2chroma_fx(Pel *srcPic, Int srcStride, Pel *dstPicX, Int dstStrideX, Pel *dstPicY, Int dstStrideY, Int width, Int height)
{
  for (Int y = 0; y < height >> 1; y++)
  {
    Pel *src   = srcPic +        (y<<1)*srcStride;
    Pel *srcUp = srcPic + (1 + (y << 1))*srcStride;
    Pel *dst_x = dstPicX + y*dstStrideX;
    Pel *dst_y = dstPicY + y*dstStrideY;

    for (Int x = 0; x < width >> 1; x++)
    {
      Int val = src[(x << 1) + 0] +   src[(x << 1) + 1] +
              srcUp[(x << 1) + 0] + srcUp[(x << 1) + 1];
      val >>= 2;
      dst_x[x] = val;
      dst_y[x] = val;
    }
  }
}

Void TComSharpeningFilter::buildEdgeMap(TComPicYuv* pcBasePic)
{
  int width = pcBasePic->getWidth(COMPONENT_Y);
  int height = pcBasePic->getHeight(COMPONENT_Y);
  //////////////////////////////////////////////////////////////////////////
  pcBasePic->setBorderExtension(false);
  pcBasePic->extendPicBorder();
  //////////////////////////////////////////////////////////////////////////

  Int **_tmpDer, **_tmpFltX, **_tmpFltY;
  initMatrix_int(&_tmpDer,  width + OFFSET_2, height + OFFSET_2);
  initMatrix_int(&_tmpFltX, width + OFFSET_2, height + OFFSET_2);
  initMatrix_int(&_tmpFltY, width + OFFSET_2, height + OFFSET_2);

  Int *tmpDer  = &_tmpDer [OFFSET][OFFSET];
  Int *tmpFltX = &_tmpFltX[OFFSET][OFFSET];
  Int *tmpFltY = &_tmpFltY[OFFSET][OFFSET];

  der_abs(pcBasePic->getAddr(COMPONENT_Y), pcBasePic->getStride(COMPONENT_Y), tmpDer, width + OFFSET_2, width, height, FLT_OFFSET + DER_OFFSET, m_bitDepths.recon[CHANNEL_TYPE_LUMA]);

  Int bitDepth = m_bitDepths.recon[toChannelType(COMPONENT_Y)];
  int saturateFix = SATURATE_FIX;
  saturateFix = 0;
  if (bitDepth == 8)
    saturateFix = 33;
  if (bitDepth == 10)
    saturateFix = 132;
  assert(saturateFix);
  apply_limits(tmpDer, width + OFFSET_2, tmpDer, width + OFFSET_2,
    //(Int)(SATURATE_FIX) << (bitDepth - 8), 
    saturateFix,
    width, height, FLT_OFFSET + DER_OFFSET);

  flt_block_x(tmpDer,  width + OFFSET_2, tmpFltX, width + OFFSET_2, width, height, DER_OFFSET + FLT_OFFSET);
  flt_block_y(tmpFltX, width + OFFSET_2, tmpFltY, width + OFFSET_2, width, height, DER_OFFSET);

  _d2xPic.destroy();
  _d2xPic.create(width, height, m_chromaFormatIDC, m_uiMaxCUWidth, m_uiMaxCUHeight, m_uiMaxCUDepth, true);
  _d2yPic.destroy();
  _d2yPic.create(width, height, m_chromaFormatIDC, m_uiMaxCUWidth, m_uiMaxCUHeight, m_uiMaxCUDepth, true);

  der_xy(tmpFltY, width + OFFSET_2,
    tmpFltX, width + OFFSET_2,
    tmpDer,  width + OFFSET_2,
    width, height);

  arr_toPic(tmpFltX, width + OFFSET_2,
    _d2xPic.getAddr(COMPONENT_Y), _d2xPic.getStride(COMPONENT_Y),
    width, height, 
    MAP_SHIFT + 2 - (m_bitDepths.recon[toChannelType(COMPONENT_Y)] - 8),
    m_bitDepths.recon[toChannelType(COMPONENT_Y)]);

  arr_toPic(tmpDer, width + OFFSET_2,
    _d2yPic.getAddr(COMPONENT_Y), _d2yPic.getStride(COMPONENT_Y),
    width, height, 
    MAP_SHIFT + 2 - (m_bitDepths.recon[toChannelType(COMPONENT_Y)] - 8),
    m_bitDepths.recon[toChannelType(COMPONENT_Y)]);

  downsample_luma2chroma_fx(
    _d2xPic.getAddr(COMPONENT_Y),  _d2xPic.getStride(COMPONENT_Y),
    _d2xPic.getAddr(COMPONENT_Cb), _d2xPic.getStride(COMPONENT_Cb),
    _d2xPic.getAddr(COMPONENT_Cr), _d2xPic.getStride(COMPONENT_Cr),
    width, height);
  downsample_luma2chroma_fx(
    _d2yPic.getAddr(COMPONENT_Y),  _d2yPic.getStride(COMPONENT_Y),
    _d2yPic.getAddr(COMPONENT_Cb), _d2yPic.getStride(COMPONENT_Cb),
    _d2yPic.getAddr(COMPONENT_Cr), _d2yPic.getStride(COMPONENT_Cr),
    width, height);

  _d2xPic.setBorderExtension(false);
  _d2xPic.extendPicBorder();
  _d2yPic.setBorderExtension(false);
  _d2yPic.extendPicBorder();

  destroyMatrix_int(_tmpDer);
  destroyMatrix_int(_tmpFltX);
  destroyMatrix_int(_tmpFltY);
}

Void TComSharpeningFilter::sharpBlockLumaForSearch(TComDataCU *pcCU, Int uiAbsPartIdx, Int uiAbsPartIdxBlock, UInt uiWidth, UInt uiHeight, int sharpMode,
  TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, TComPicYuv *refPicYuv)
{
  assert(SIF_GAUSS_ONLY == 0);
  // interpolate map
  TComYuv *mapIterX = &_d2xYuv;
  TComYuv *mapIterY = &_d2yYuv;
  comPred->xPredInterBlk(COMPONENT_Y, pcCU, &_d2xPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterX, false, m_bitDepths.recon[CHANNEL_TYPE_LUMA]);
  comPred->xPredInterBlk(COMPONENT_Y, pcCU, &_d2yPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterY, false, m_bitDepths.recon[CHANNEL_TYPE_LUMA]);

  // displace
  float k; int x, y;
  getSharpParams(sharpMode, pcCU, uiAbsPartIdx, k, x, y);

  dstPic->copyPartToPartComponentDst(COMPONENT_Y, &_picYuv, uiAbsPartIdxBlock, uiAbsPartIdx, uiWidth, uiHeight);
  displaceBlk(&_picYuv, &_picYuvTmp, &_d2xYuv, &_d2yYuv, MAP_SHIFT, m_bitDepths, COMPONENT_Y, k, x, y, uiWidth, uiHeight, !bi);

  // TODO: optimize copying
  _picYuvTmp.copyPartToPartComponentDst(COMPONENT_Y, dstPic, uiAbsPartIdx, uiAbsPartIdxBlock, uiWidth, uiHeight);
}

Void TComSharpeningFilter::sharpBlockLuma(TComDataCU *pcCU, Int uiAbsPartIdx, UInt uiWidth, UInt uiHeight, int sharpMode,
  TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, Bool isDecoder, TComPicYuv *refPicYuv)
{
  assert(SIF_GAUSS_ONLY == 0);
  // interpolate map
  TComYuv *mapIterX_ = &_d2xYuvInter;
  TComYuv *mapIterY_ = &_d2yYuvInter;

#if RRC_BLOCK_SHARPENING_DECODER
  if (!isDecoder)
  {
#endif
    comPred->xPredInterBlk(COMPONENT_Y, pcCU, &_d2xPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterX_, false, m_bitDepths.recon[CHANNEL_TYPE_LUMA]);
    comPred->xPredInterBlk(COMPONENT_Y, pcCU, &_d2yPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterY_, false, m_bitDepths.recon[CHANNEL_TYPE_LUMA]);
#if RRC_BLOCK_SHARPENING_DECODER
  }
  else
  {
    comPred->xPredInterLumaBlk(pcCU, &_d2xYuv, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterX_, false, m_is_left_chroma_alignment, m_is_top_chroma_alignment);
    comPred->xPredInterLumaBlk(pcCU, &_d2yYuv, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterY_, false, m_is_left_chroma_alignment, m_is_top_chroma_alignment);
  }
#endif

  // displace
  float k; int x, y;
  getSharpParams(sharpMode, pcCU, uiAbsPartIdx, k, x, y);
  dstPic->copyPartToPartComponent(COMPONENT_Y, &_picYuv, uiAbsPartIdx, uiWidth, uiHeight);
  displaceBlk(&_picYuv, dstPic, &_d2xYuvInter, &_d2yYuvInter, MAP_SHIFT, m_bitDepths, COMPONENT_Y, k, x, y, uiWidth, uiHeight, !bi);
}

Void TComSharpeningFilter::sharpBlockChroma(TComDataCU *pcCU, Int uiAbsPartIdx, UInt uiWidth, UInt uiHeight, int sharpMode,
  TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, Bool isDecoder, TComPicYuv *refPicYuv)
{
  assert(SIF_GAUSS_ONLY == 0);
#if RRC_BLOCK_SHARPENING_DECODER  
  BlockDescription block(pcCU, uiAbsPartIdx, mv, uiWidth, uiHeight, pcCU->getPic()->getPicYuvRec()->getWidth(), pcCU->getPic()->getPicYuvRec()->getHeight());
  BlockDescription block_ = block.create_extended(_interp_filter_margin, true);

  if (block_.is_on_frame_border())
  {
    return;
  }
#endif

  // interpolate map
  TComYuv *mapIterX_ = &_d2xYuvInter;
  TComYuv *mapIterY_ = &_d2yYuvInter;

#if RRC_BLOCK_SHARPENING_DECODER
  if (!isDecoder)
  {
#endif
    comPred->xPredInterBlk(COMPONENT_Cb, pcCU, &_d2xPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterX_, false, m_bitDepths.recon[CHANNEL_TYPE_CHROMA]);
    comPred->xPredInterBlk(COMPONENT_Cb, pcCU, &_d2yPic, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterY_, false, m_bitDepths.recon[CHANNEL_TYPE_CHROMA]);
#if RRC_BLOCK_SHARPENING_DECODER
  }
  else
  {
    comPred->xPredInterChromaBlk(pcCU, &_d2xYuv, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterX_, false);
    comPred->xPredInterChromaBlk(pcCU, &_d2yYuv, uiAbsPartIdx, mv, uiWidth, uiHeight, mapIterY_, false);
  }
#endif

  // displace
  float k; int x, y;
  getSharpParams(sharpMode, pcCU, uiAbsPartIdx, k, x, y);

  // TODO: it seems non need to shift size here
  dstPic->copyPartToPartComponent(COMPONENT_Cb, &_picYuv, uiAbsPartIdx, uiWidth >> dstPic->getComponentScaleX(COMPONENT_Cb), uiHeight >> dstPic->getComponentScaleY(COMPONENT_Cb));
  dstPic->copyPartToPartComponent(COMPONENT_Cr, &_picYuv, uiAbsPartIdx, uiWidth >> dstPic->getComponentScaleX(COMPONENT_Cr), uiHeight >> dstPic->getComponentScaleY(COMPONENT_Cr));

  displaceBlk(&_picYuv, dstPic, &_d2xYuvInter, &_d2yYuvInter, MAP_SHIFT, m_bitDepths, COMPONENT_Cb, k*2, x >> 1, y >> 1, uiWidth >> 1, uiHeight >> 1, !bi);
  displaceBlk(&_picYuv, dstPic, &_d2xYuvInter, &_d2yYuvInter, MAP_SHIFT, m_bitDepths, COMPONENT_Cr, k*2, x >> 1, y >> 1, uiWidth >> 1, uiHeight >> 1, !bi);
}

Void TComSharpeningFilter::getSharpParams(int sharpMode, TComDataCU *pcCU, Int uiAbsPartIdx, float &k, int &x, int &y)
{
  assert(sharpMode > 0);
  k = 0.06f;
  x = g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
  y = g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
}

//////////////////////////////////////////////////////////////////////////////////////
////////////////////////class BlockDescription Definitions////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

BlockDescription::BlockDescription(unsigned int frame_width, unsigned int frame_height)
  : m_width(frame_width)
  , m_height(frame_height)
  , m_frame_width(frame_width)
  , m_frame_height(frame_height)
  , m_pos_x(0)
  , m_pos_y(0)
  , m_offset(0)
  , m_is_left_border(true)
  , m_is_top_border(true)
  , m_is_right_border(true)
  , m_is_bottom_border(true)
  , m_left_extension(0)
  , m_right_extension(0)
  , m_top_extension(0)
  , m_bottom_extension(0)
  , m_is_left_chroma_alignment(false)
  , m_is_top_chroma_alignment(false)
  , m_bottom_border_min_dist(0)
  , m_top_border_min_dist(0)
  , m_left_border_min_dist(0)
  , m_right_border_min_dist(0)
{}

BlockDescription::BlockDescription(TComDataCU* cu, UInt partAddr, TComMv* mv, unsigned int width, unsigned int height, unsigned int frame_width, unsigned int frame_height)
  : m_width(width)
  , m_height(height)
  , m_frame_width(frame_width)
  , m_frame_height(frame_height)
  , m_is_left_border(false)
  , m_is_top_border(false)
  , m_is_right_border(false)
  , m_is_bottom_border(false)
  , m_left_extension(0)
  , m_right_extension(0)
  , m_top_extension(0)
  , m_bottom_extension(0)
  , m_is_left_chroma_alignment(false)
  , m_is_top_chroma_alignment(false)
  , m_bottom_border_min_dist(0)
  , m_top_border_min_dist(0)
  , m_left_border_min_dist(0)
  , m_right_border_min_dist(0)
{
  UInt uiMaxCUWidth  = MAX_CU_SIZE;
  UInt uiMaxCUHeight = MAX_CU_SIZE;
  UInt uiMaxCUDepth = cu->getSlice()->getSPS()->getMaxTotalCUDepth();

  const unsigned int width_in_lcu = (m_frame_width / uiMaxCUWidth) + ((m_frame_width % uiMaxCUWidth) != 0);
  const unsigned int lcu_offset = (cu->getCtuRsAddr() / width_in_lcu) * m_frame_width * uiMaxCUHeight + (cu->getCtuRsAddr() % width_in_lcu) * uiMaxCUWidth;
  const unsigned int pu_offset = 4 * (g_auiZscanToRaster[cu->getZorderIdxInCtu() + partAddr] / (uiMaxCUWidth / 4)) * m_frame_width + 4 * (g_auiZscanToRaster[cu->getZorderIdxInCtu() + partAddr] % (uiMaxCUWidth / 4));

  m_offset = lcu_offset + pu_offset;
  m_pos_y = m_offset / m_frame_width;
  m_pos_x = m_offset % m_frame_width;

#if QC_MV_STORE_PRECISION_BIT 
  const int mv_x = mv->getHor() >> QC_MV_STORE_PRECISION_BIT;
  const int mv_y = mv->getVer() >> QC_MV_STORE_PRECISION_BIT;
#else
  const int mv_x = mv->getHor() >> 2;
  const int mv_y = mv->getVer() >> 2;
#endif

  if (m_pos_x + mv_x <= 0)  //Block out of left frame border. It covers partial and whole outing cases. 
  {
    m_left_extension = -(mv_x + m_pos_x);
    m_width -= m_left_extension;
    if (m_width < 0)  //Whole block is out of frame case. It requires to compute just one column inside frame for futher margin filling.
    {
      m_left_extension += m_width;  //Excluding free space between left frame border and right block border. It doesn't necessary to store fair extension becasue whole margin is a repetition of the one real value.
      m_left_border_min_dist = -(m_width - 1);  //Computing distance between left frame border and right block border considering free space in the middle. This information will be used in futher block extensions.
      m_width = 1;  //Creating one real column inside of the frame for futhre margin filling.
    }
    m_pos_x = 0;
    m_is_left_border = true;
  }
  else if (m_pos_x + mv_x + m_width >= m_frame_width - 1) //Block out of right frame border. It covers partial and whole outing cases. 
  {
    m_right_extension = m_pos_x + mv_x + m_width - (m_frame_width - 1);
    m_width = (frame_width - 1) - (m_pos_x + mv_x);
    if (m_width < 0)  //Whole block is out of right frame border case. It requires to compute just one column inside frame for futher margin filling.
    {
      m_right_extension += m_width;  //Excluding free space between right frame border and left block border. It doesn't necessary to store fair extension becasue whole margin is a repetition of the one real value.
      m_right_border_min_dist = -m_width;
      m_width = 1;
      m_pos_x = frame_width - 1;
    }
    else
    {
      m_pos_x += mv_x;
    }
    m_is_right_border = true;
  }
  else
  {
    m_pos_x += mv_x;
  }

  if (m_pos_y + mv_y <= 0)
  {
    m_top_extension = -(mv_y + m_pos_y);
    m_height -= m_top_extension;
    if (m_height < 0)
    {
      m_top_extension += m_height;
      m_top_border_min_dist = -(m_height - 1);
      m_height = 1;
    }
    m_pos_y = 0;
    m_is_top_border = true;
  }
  else if (m_pos_y + mv_y + m_height >= m_frame_height - 1)
  {
    m_bottom_extension = m_pos_y + mv_y + m_height - (m_frame_height - 1);
    m_height = (frame_height - 1) - (m_pos_y + mv_y);
    if (m_height < 0)
    {
      m_bottom_extension -= m_height;
      m_bottom_border_min_dist = -m_height;
      m_height = 1;
      m_pos_y = frame_height - 1;
    }
    else
    {
      m_pos_y += mv_y;
    }
    m_is_bottom_border = true;
  }
  else
  {
    m_pos_y += mv_y;
  }

  //alignment for next chroma downsampling
  if (!is_on_frame_border())
  {
    if (m_pos_x % 2 != 0)
    {
      --m_pos_x;
      ++m_width;

      m_is_left_chroma_alignment = true;
    }

    if (m_pos_y % 2 != 0)
    {
      --m_pos_y;
      ++m_height;

      m_is_top_chroma_alignment = true;
    }
  }

  m_offset = m_pos_y * m_frame_width + m_pos_x;
}

BlockDescription BlockDescription::create_extended(unsigned int margin, bool is_real_extension) const
{
  BlockDescription res(*this);

  //left extension 
  if (is_real_extension || m_right_border_min_dist == 0)  //Regular case. Extening the real part of the block (inside the frame).
  {
    if (m_pos_x > margin)
    {
      res.m_pos_x -= margin;
      res.m_offset -= margin;
      res.m_width += margin;
    }
    else
    {
      res.m_offset -= m_pos_x;
      res.m_width += m_pos_x;
      res.m_left_extension += margin - m_pos_x;
      res.m_pos_x = 0;
      res.m_is_left_border = true;
    }
  }
  else  // Whole block is out of frame case.
  {
    if (m_right_border_min_dist > margin)
    {
      res.m_right_border_min_dist -= margin;
      res.m_right_extension += margin;
    }
    else
    {
      res.m_right_border_min_dist = 0;
      res.m_pos_x -= margin - m_right_border_min_dist;
      res.m_width += margin - m_right_border_min_dist;
      res.m_right_extension += m_right_border_min_dist - 1;
    }
  }

  //right extension 
  if (is_real_extension || m_left_border_min_dist == 0)
  {
    if (m_pos_x + m_width + margin < m_frame_width)
    {
      res.m_width += margin;
    }
    else
    {
      res.m_width += m_frame_width - (m_pos_x + m_width);
      res.m_right_extension += (m_pos_x + m_width + margin - m_frame_width);
      res.m_is_right_border = true;
    }
  }
  else
  {
    if (m_left_border_min_dist > margin)
    {
      res.m_left_border_min_dist -= margin;
      res.m_left_extension += margin;
    }
    else
    {
      res.m_left_border_min_dist = 0;
      res.m_width += margin - m_left_border_min_dist;
      res.m_left_extension += m_left_border_min_dist - 1;
    }
  }

  //top extension 
  if (is_real_extension || m_bottom_border_min_dist == 0)
  {
    if (m_pos_y > margin)
    {
      res.m_pos_y -= margin;
      res.m_offset -= margin * m_frame_width;
      res.m_height += margin;
    }
    else
    {
      res.m_offset -= m_pos_y * m_frame_width;
      res.m_height += m_pos_y;
      res.m_top_extension += margin - m_pos_y;
      res.m_pos_y = 0;
      res.m_is_top_border = true;
    }
  }
  else
  {
    if (m_bottom_border_min_dist > margin)
    {
      res.m_bottom_border_min_dist -= margin;
      res.m_bottom_extension += margin;
    }
    else
    {
      res.m_bottom_border_min_dist = 0;
      res.m_pos_y -= margin - m_bottom_border_min_dist;
      res.m_height += margin - m_bottom_border_min_dist;
      res.m_bottom_extension += m_bottom_border_min_dist - 1;
    }
  }

  //bottom extension 
  if (is_real_extension || m_top_border_min_dist == 0)
  {
    if (m_pos_y + m_height + margin < m_frame_height)
    {
      res.m_height += margin;
    }
    else
    {
      res.m_height += m_frame_height - (m_pos_y + m_height);
      res.m_bottom_extension += m_pos_y + m_height + margin - m_frame_height;
      res.m_is_bottom_border = true;
    }
  }
  else
  {
    if (m_top_border_min_dist > margin)
    {
      res.m_top_border_min_dist -= margin;
      res.m_top_extension += margin;
    }
    else
    {
      res.m_top_border_min_dist = 0;
      res.m_height += margin - m_top_border_min_dist;
      res.m_top_extension += m_top_border_min_dist - 1;
    }
  }

  res.m_offset = res.m_pos_y * res.m_frame_width + res.m_pos_x;
  return res;
}

void BlockDescription::check_y_border(int* y) const
{
  if (m_pos_y + *y < 0)
    *y = 0;
  else if (m_pos_y + *y >= m_frame_height)
    *y = m_frame_height - 1 - m_pos_y;
}

void BlockDescription::check_x_border(int* x) const
{
  if (m_pos_x + *x < 0)
    *x = 0;
  else if (m_pos_x + *x >= m_frame_width)
    *x = m_frame_width - 1 - m_pos_x;
}
#endif//RRC_SELECTABLE_IF
