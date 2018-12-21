
#pragma once

#include <vector>
#include "CommonDef.h"
#include <cmath>

#if RRC_NS      

#define SQR(X) ((X) * (X))
#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

struct NsParams
{
  unsigned nHW;
  unsigned kHW;
  unsigned NHW;
  unsigned pHW;
  int      bit_depth;

  unsigned width;
  unsigned height;
  unsigned width_chroma;
  unsigned height_chroma;
};

class NsImage
{
  int _width;
  int _height;
  int _width_chroma;
  int _height_chroma;
public:
  std::vector<RRC_NS_PEL_TYPE> m_Y;
  std::vector<RRC_NS_PEL_TYPE> m_U;
  std::vector<RRC_NS_PEL_TYPE> m_V;

  NsImage(int width, int height, int width_chroma, int height_chroma)
  {
    m_Y = std::vector<RRC_NS_PEL_TYPE>(width * height);
    _width = width;
    _height = height;

    m_U = std::vector<RRC_NS_PEL_TYPE>(width_chroma * height_chroma);
    m_V = std::vector<RRC_NS_PEL_TYPE>(width_chroma * height_chroma);
    _width_chroma = width_chroma;
    _height_chroma = height_chroma;
  }

  NsImage() : _width(0), _height(0), _width_chroma(0), _height_chroma(0)
  {

  }

  int width() const
  {
    return _width;
  }

  int height() const
  {
    return _height;
  }

  int width_chroma() const
  {
    return _width_chroma;
  }

  int height_chroma() const
  {
    return _height_chroma;
  }

};



int run_ns(
    const float sigma
,   NsImage &img_noisy
,   NsImage &img_denoised
,   const int bit_depth
);

#endif // RRC_NS

