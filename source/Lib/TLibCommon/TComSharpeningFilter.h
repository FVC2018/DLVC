#ifndef __TCOMSHARPENINGFILTER__
#define __TCOMSHARPENINGFILTER__

#include "TComPic.h"
#include "TComYuv.h"

class TComPrediction;

class TComSharpeningFilter
{
public:
  TComSharpeningFilter(TComPicYuv* pcBasePic, const TComSPS* sps, bool isEncoder = 0);
  ~TComSharpeningFilter();

  Void buildEdgeMap(TComPicYuv* pcBasePic);
  Void buildBlockEdgeMap(TComPicYuv* pcBasePic, TComDataCU* cu, UInt partAddr, TComMv* mv, UInt uiWidth, UInt uiHeight);
  Void getSharpParams(int sharpMode, TComDataCU *pcCU, Int uiAbsPartIdx, float &k, int &x, int &y);

  Void sharpBlockLuma(TComDataCU *pcCU, Int uiAbsPartIdx, UInt uiWidth, UInt uiHeight, int sharpMode, TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, Bool isDecoder, TComPicYuv *refPicYuv);
  Void sharpBlockChroma(TComDataCU *pcCU, Int uiAbsPartIdx, UInt uiWidth, UInt uiHeight, int sharpMode, TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, Bool isDecoder, TComPicYuv *refPicYuv);
  Void sharpBlockLumaForSearch(TComDataCU *pcCU, Int uiAbsPartIdx, Int uiAbsPartIdxBlock, UInt uiWidth, UInt uiHeight, int sharpMode, TComYuv *&dstPic, TComPrediction* comPred, TComMv *mv, Bool bi, TComPicYuv *refPicYuv);

private:
  bool _isEncoder;
  float* _src;
  float* _d2x;
  float* _d2y;
  float* _dx;
  float* _dy;
  float* _d2x_chroma;
  float* _d2y_chroma;
  float* _dxy_abs;
  float* _tmp;

  bool m_is_left_chroma_alignment;
  bool m_is_top_chroma_alignment;

  UInt m_uiMaxCUWidth;
  UInt m_uiMaxCUHeight;
  UInt m_uiMaxCUDepth;
  ChromaFormat m_chromaFormatIDC;
  BitDepths m_bitDepths;

  TComPicYuv _d2xPic;
  TComPicYuv _d2yPic;

  TComYuv _picYuv;
  TComYuv _picYuvTmp;
  TComYuv _d2xYuv;
  TComYuv _d2yYuv;
  TComYuv _d2xYuvInter;
  TComYuv _d2yYuvInter;

  static const int _interp_filter_margin;
  static const int _internal_extension;

  void float2pic(int shift);
  void clear_buffers();

};

class BlockDescription
{
  int m_pos_x;
  int m_pos_y;
  int m_width;
  int m_height;

  int m_top_border_min_dist;        //Distance between top frame border and bottom block border if whole block is oot of frame, zero otherwise.
  int m_bottom_border_min_dist;     //Distance between bottom frame border and top block border if whole block is oot of frame, zero otherwise.
  int m_left_border_min_dist;       //Distance between left frame border and right block border if whole block is oot of frame, zero otherwise.
  int m_right_border_min_dist;      //Distance between right frame border and left block border if whole block is oot of frame, zero otherwise.

  int m_frame_width;
  int m_frame_height;
  int m_offset;

  bool m_is_left_border;
  bool m_is_top_border;
  bool m_is_right_border;
  bool m_is_bottom_border;

  int m_left_extension;
  int m_right_extension;
  int m_top_extension;
  int m_bottom_extension;

  bool m_is_left_chroma_alignment;
  bool m_is_top_chroma_alignment;
public:
  BlockDescription(unsigned int frame_width, unsigned int frame_height);
  BlockDescription(TComDataCU* cu, UInt partAddr, TComMv* mv, unsigned int width, unsigned int height, unsigned int frame_width, unsigned int frame_height);
  BlockDescription create_extended(unsigned int margin, bool is_real_extension = false) const;
  void check_y_border(int* y) const;
  void check_x_border(int* x) const;

  int get_width() const { return m_width; }
  int get_height() const { return m_height; }
  int get_offset() const { return m_offset; }
  int get_pos_x() const { return m_pos_x; }
  int get_pos_y() const { return m_pos_y; }
  int get_left_extension() const { return m_left_extension; }
  int get_right_extension() const { return m_right_extension; }
  int get_top_extension() const { return m_top_extension; }
  int get_bottom_extension() const { return m_bottom_extension; }
  int get_top_border_displacement() const { return m_top_border_min_dist; }
  int get_bottom_border_displacement() const { return m_bottom_border_min_dist; }
  int get_left_border_displacement() const { return m_left_border_min_dist; }
  int get_right_border_displacement() const { return m_right_border_min_dist; }
  bool get_left_chroma_alignment() const { return m_is_left_chroma_alignment; }
  bool get_top_chroma_alignment() const { return m_is_top_chroma_alignment; }
  bool is_on_frame_border() const  { return m_is_left_border || m_is_right_border || m_is_top_border || m_is_bottom_border; }
};

#endif //__TCOMSHARPENINGFILTER__
