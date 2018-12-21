
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <immintrin.h>

#include "TComNoiseSuppressor.h"

#if RRC_NS     
#if RRC_NS_OPTIMIZE_OFF
#pragma optimize( "", off )
#endif

#if defined(_MSC_VER)
#define ALIGNED_(x) __declspec(align(x))
#define FORCE_INLINE __forceinline
#else
#if defined(__GNUC__)
#define ALIGNED_(x) __attribute__ ((aligned(x)))
#define FORCE_INLINE __attribute__((always_inline))
#endif
#endif

using namespace std;
using namespace ns;

/////////////////////////////////////////////////////////////////////////////////////
/////////////                                                           /////////////
/////////////                FORWARD DECLARATIONS                       /////////////
/////////////                                                           /////////////
/////////////////////////////////////////////////////////////////////////////////////

typedef struct
{
  unsigned Size;
  unsigned PixNum;
  unsigned Patch[RRC_NS_NHW];
} tPatchTable;

typedef struct
{
  RRC_NS_TABLE_DISTANCE_TYPE Mse;
  unsigned Patch;
}tPd;

typedef struct
{
  unsigned Size;
  tPd PatchTbl[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1) / 2];
}tTableDistance;

typedef struct
{
  unsigned Size;
  unsigned *pPxlNum;
}tLattice;

typedef struct
{
  NsParams ns_params;
  short *denominator;
  int *numerator;

  tLattice row_ind;
  tLattice column_ind;

  tPatchTable *patch_table;

  tTableDistance table_distance;

  NsImage img_sym_noisy;
}tNs;

void *NsCreate(NsParams *pNsParams);

void NsBlockMaching(tNs *pNs, const float sigma, const NsImage & img_noisy);
void NsFltLuma(tNs *pNs, const float sigma, NsImage & img_noisy);
void NsFltChroma(tNs *pNs, const float sigma, NsImage & img_noisy);

void wiener_filtering_hadamard4(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const size_t nSx_r, const unsigned kWien, float sigma, const int width);

int closest_power_of_2(size_t n);

void ind_initialize(
  tLattice &ind_set
  , const unsigned max_size
  , const unsigned N
  , const unsigned step
  );

///////////////////////////////////////////////////////////////////
/////////////////NsClassifier class implementation/////////////////
///////////////////////////////////////////////////////////////////

#define CHECK_THR (1.0)
static unsigned Check(const RRC_NS_PEL_TYPE *data, unsigned const x, unsigned const y, unsigned const stride, unsigned const size)
{
  double Mean = 0;
  double Var = 0;

  int iMean = 0;
  int iVar = 0;
  const RRC_NS_PEL_TYPE *pData = &data[y*stride + x];

  for (unsigned i = 0; i < size; ++i, pData += stride)
  {
    for (unsigned j = 0; j < size; ++j)
    {
      const int tmp = pData[j];
      iMean += tmp;
      iVar += tmp * tmp;
    }
  }

  {
    const unsigned S2 = size*size;
    Mean = (double)iMean / S2;
    Var = (double)iVar / S2;
    Var -= Mean*Mean;
  }

  if (Var > CHECK_THR)
    return 1;
  else
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////                                                           /////////////
/////////////                          HELPERS                          /////////////
/////////////                                                           /////////////
/////////////////////////////////////////////////////////////////////////////////////

void symetrize(
  const std::vector<RRC_NS_PEL_TYPE> &img
  , std::vector<RRC_NS_PEL_TYPE> &img_sym
  , const unsigned width
  , const unsigned height
  , const unsigned N
  ){
  //! Declaration
  const unsigned w = width + 2 * N;
  const unsigned h = height + 2 * N;

  if (img_sym.size() != w * h)
    img_sym.resize(w * h);

  {
    unsigned dc = 0;
    unsigned dc_2 = N * w + N;

    //! Center of the image
    for (unsigned i = 0; i < height; i++)
      for (unsigned j = 0; j < width; j++, dc++)
        img_sym[dc_2 + i * w + j] = img[dc];

    //! Top and bottom
    dc_2 = 0;
    for (unsigned j = 0; j < w; j++, dc_2++)
      for (unsigned i = 0; i < N; i++)
      {
        img_sym[dc_2 + i * w] = img_sym[dc_2 + (2 * N - i - 1) * w];
        img_sym[dc_2 + (h - i - 1) * w] = img_sym[dc_2 + (h - 2 * N + i) * w];
      }

    //! Right and left
    dc_2 = 0;
    for (unsigned i = 0; i < h; i++)
    {
      const unsigned di = dc_2 + i * w;
      for (unsigned j = 0; j < N; j++)
      {
        img_sym[di + j] = img_sym[di + 2 * N - j - 1];
        img_sym[di + w - j - 1] = img_sym[di + w - 2 * N + j];
      }
    }
  }

  return;
}


int closest_power_of_2(size_t n)
{
  unsigned r = 1;
  while (r * 2 <= n)
    r *= 2;

  return r;
}

void ind_initialize(
  tLattice &ind_set
  , const unsigned max_size
  , const unsigned N
  , const unsigned step
  )
{
  unsigned Size = 0;
  unsigned ind = N;
  while (ind < max_size - N)
  {
    Size++;
    ind += step;
  }
  if (ind - step < max_size - N - 1)
  {
    Size++;
  }
  ind_set.pPxlNum = new unsigned[Size];

  ind_set.Size = 0;
  ind = N;
  while (ind < max_size - N)
  {
    assert(ind_set.Size < Size);
    ind_set.pPxlNum[ind_set.Size++] = ind;
    ind += step;
  }
  if (ind_set.pPxlNum[ind_set.Size - 1] < max_size - N - 1)
  {
    assert(ind_set.Size < Size);
    ind_set.pPxlNum[ind_set.Size++] = (max_size - N - 1);
  }
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////                                                           /////////////
/////////////                     BLOCK MATCHING                        /////////////
/////////////                                                           /////////////
/////////////////////////////////////////////////////////////////////////////////////


#if !RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const float threshold)
{
  RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

  const RRC_NS_PEL_TYPE *pSrc = &img[k0];
  const RRC_NS_PEL_TYPE *pSrc1 = &img[k1];

  for (size_t p = 0; p < kHW && sum <= threshold; p += 1, pSrc += width, pSrc1 += width)
  {
    for (size_t q = 0; q < kHW; q += 1)
    {
#if RRC_NS_USE_SAD
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = img[k0 + p * width + q] - img[k1 + p * width + q];

      sum += (RRC_NS_TABLE_DISTANCE_TYPE)fabs((float)tmp);
#else
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = pSrc[q] - pSrc1[q];

      sum += tmp * tmp;
    }
  }
#endif //RRC_NS_USE_SAD
  return sum;
}
#endif

#if !RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance_short(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const float threshold)
{
  RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

  const RRC_NS_PEL_TYPE *pSrc = &img[k0];
  const RRC_NS_PEL_TYPE *pSrc1 = &img[k1];

  for (size_t p = 0; p < kHW / 2 && sum <= threshold; p += 1, pSrc += width, pSrc1 += width)
  {
    for (size_t q = 0; q < kHW; q += 1)
    {
#if RRC_NS_USE_SAD
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = img[k0 + p * width + q] - img[k1 + p * width + q];

      sum += (RRC_NS_TABLE_DISTANCE_TYPE)fabs((float)tmp);
#else
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = pSrc[q] - pSrc1[q];

      sum += tmp * tmp;
    }
  }
#endif //RRC_NS_USE_SAD
  return sum;
}
#endif

#if !RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance_short1(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const float threshold)
{
  RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

  const RRC_NS_PEL_TYPE *pSrc = &img[k0];
  const RRC_NS_PEL_TYPE *pSrc1 = &img[k1];

  for (size_t p = 0; p < kHW / 2; p += 1, pSrc += width, pSrc1 += width)
  {
    for (size_t q = 0; q < kHW; q += 1)
    {
#if RRC_NS_USE_SAD
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = img[k0 + p * width + q] - img[k1 + p * width + q];

      sum += (RRC_NS_TABLE_DISTANCE_TYPE)fabs((float)tmp);
#else
      const RRC_NS_TABLE_DISTANCE_TYPE tmp = pSrc[q] - pSrc1[q];

      sum += tmp * tmp;
    }
  }
#endif //RRC_NS_USE_SAD
  return sum;
}
#endif


#if RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const RRC_NS_TABLE_DISTANCE_TYPE threshold)
{

  if (kHW == 8)
  {

    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];


    for (size_t p = 0; p < 2; p += 1)
    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      {
        ALIGNED_(16) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

        if (sum > threshold)
          return sum;
      }

      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      i1 += width;
      i2 += width;

#if 0
      {
        ALIGNED_(16)) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

        if (sum > threshold)
          //break;
          return sum;
      }
#endif

      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      i1 += width;
      i2 += width;


      {
        ALIGNED_(16) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

        if (sum > threshold)
          return sum;
      }

    }

    return sum;
  }
  else
  {
    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];

    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      {
        ALIGNED_(16) int buf[2];
        _mm_storel_epi64((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1];

        if (sum > threshold)
          return sum;
      }

      mtmp1 = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      i1 += width;
      i2 += width;

      mtmp = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      mtmp1 = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      {
        ALIGNED_(16) int buf[2];
        _mm_storel_epi64((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1];
      }
    }
    return sum;
  }

}

#endif //RRC_NS_USE_SSE

#if RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance_short(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const RRC_NS_TABLE_DISTANCE_TYPE threshold)
{

  if (kHW == 8)
  {

    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];

    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      {
        ALIGNED_(16) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

        if (sum > threshold)
          return sum;
      }

      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      i1 += width;
      i2 += width;


      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      {
        ALIGNED_(16) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

      }

    }

    return sum;
  }
  else
  {
    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];

    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      {
        ALIGNED_(16) int buf[2];
        _mm_storel_epi64((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1];

        if (sum > threshold)
          return sum;
      }

      mtmp1 = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      {
        ALIGNED_(16) int buf[2];
        _mm_storel_epi64((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1];
      }
    }
    return sum;
  }

}
#endif //RRC_NS_USE_SSE

#if RRC_NS_USE_SSE
FORCE_INLINE RRC_NS_TABLE_DISTANCE_TYPE patch_distance_short1(const vector<RRC_NS_PEL_TYPE> &img, unsigned k0, unsigned k1, unsigned kHW, unsigned width, const RRC_NS_TABLE_DISTANCE_TYPE threshold)
{

  if (kHW == 8)
  {

    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];

    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;
      /*
      {
      ALIGNED_(16) int buf[4];
      _mm_store_si128((__m128i *)buf, sum_ss);
      sum = buf[0] + buf[1] + buf[2] + buf[3];

      if (sum > threshold)
      return sum;
      }
      */
      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      i1 += width;
      i2 += width;


      mtmp = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;

      mtmp1 = _mm_sub_epi16(_mm_loadu_si128((__m128i*)i1), _mm_loadu_si128((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      {
        ALIGNED_(16) int buf[4];
        _mm_store_si128((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1] + buf[2] + buf[3];

      }

    }

    return sum;
  }
  else
  {
    RRC_NS_TABLE_DISTANCE_TYPE sum = 0;

    __m128i sum_ss = { 0 };

    const short * i1 = &img[k0];
    const short * i2 = &img[k1];

    {
      __m128i mtmp, mtmp1;

      mtmp = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp, mtmp));

      i1 += width;
      i2 += width;
      /*
      {
      ALIGNED_(16) int buf[2];
      _mm_storel_epi64((__m128i *)buf, sum_ss);
      sum = buf[0] + buf[1];

      if (sum > threshold)
      return sum;
      }
      */
      mtmp1 = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)i1), _mm_loadl_epi64((__m128i*)i2));
      sum_ss = _mm_add_epi32(sum_ss, _mm_madd_epi16(mtmp1, mtmp1));

      {
        ALIGNED_(16) int buf[2];
        _mm_storel_epi64((__m128i *)buf, sum_ss);
        sum = buf[0] + buf[1];
      }
    }
    return sum;
  }

}
#endif //RRC_NS_USE_SSE


FORCE_INLINE void put(tTableDistance & table_distance, tPatchTable *patch_table, int NHW, int k_r)
{
  const unsigned nSx_r = min(closest_power_of_2(table_distance.Size + 1), NHW);
  tPatchTable *pPatchTable = &patch_table[k_r];

  pPatchTable->Patch[pPatchTable->Size++] = pPatchTable->PixNum;
  for (unsigned n = 0; n + 1< nSx_r; n++)
  {
    pPatchTable->Patch[pPatchTable->Size++] = table_distance.PatchTbl[n].Patch;
  }
}


FORCE_INLINE static void insertionSort(tTableDistance & table_distance)
{
  tPd temp;
  int item;
  for (int counter = 1; counter < table_distance.Size; counter++)
  {
    temp = table_distance.PatchTbl[counter];
    item = counter - 1;
    while (item >= 0 && table_distance.PatchTbl[item].Mse > temp.Mse)
    {
      table_distance.PatchTbl[item + 1] = table_distance.PatchTbl[item];
      table_distance.PatchTbl[item] = temp;
      item--;
    }
  }
}

FORCE_INLINE static void insertionSort(tTableDistance & table_distance, const unsigned MaxIdx)
{
  tPd temp;
  int item;
  for (unsigned counter = 1; counter < table_distance.Size; counter++)
  {
    temp = table_distance.PatchTbl[counter];
    item = min(counter - 1, MaxIdx);
    while (item >= 0 && table_distance.PatchTbl[item].Mse > temp.Mse)
    {
      table_distance.PatchTbl[item + 1] = table_distance.PatchTbl[item];
      table_distance.PatchTbl[item] = temp;
      item--;
    }
  }
}
FORCE_INLINE void find_patches_fcfs_romb_fst(unsigned k_r
  , float threshold
  , tTableDistance & table_distance
  , const NsImage &img
  , const NsParams & bp
  )
{
  const unsigned width = img.width();
  const unsigned height = img.height();
  const unsigned kHW = bp.kHW;
  const unsigned nHW = bp.nHW;
  const unsigned NHW = bp.NHW;
  const unsigned pHW = bp.pHW;

  static int Mask[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1)];
  static unsigned Size = 0;

  float cut_thres = threshold;
  float MaxVal = -1;
  unsigned Idx = 0;
  unsigned Sz = 0;


  if (Size == 0)
  {
    const int N = 1;

    for (int di = -N; di <= N; ++di)
    {
      for (int dj = -N; dj <= N; ++dj)
      {

        if (dj == 0 && di == 0)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }

    for (int di = -int(nHW); di <= int(nHW); ++di)
    {
      for (int dj = -int(nHW); dj <= int(nHW); ++dj)
      {
        if (abs(di) + abs(dj) > nHW)
        {
          continue;
        }

        if (abs(dj) <= N && abs(di) <= N)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }
  }

  for (unsigned k = 0; k < Size; k++)
  {

    float sum = (float)patch_distance(img.m_Y, k_r, k_r + Mask[k], kHW, img.width(), (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);

    if (sum > cut_thres)
      continue;

    if (table_distance.Size < (NHW - 1))
      //    if (table_distance.size() < (NHW-1+1))
      //    if (table_distance.size() < (NHW-1+2))
    {
      table_distance.PatchTbl[table_distance.Size].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
      table_distance.PatchTbl[table_distance.Size++].Patch = k_r + Mask[k];

      if (sum > MaxVal)
      {
        MaxVal = sum;
        Idx = Sz;
      }
      Sz++;
    }
    else
    {
      if (sum <= MaxVal)
      {
        table_distance.PatchTbl[Idx].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
        table_distance.PatchTbl[Idx].Patch = k_r + Mask[k];

        MaxVal = -1;
        for (unsigned i = 0; i < (NHW - 1); i++)
          //        for (unsigned i = 0; i < (NHW - 1+1); i++)
          //        for (unsigned i = 0; i < (NHW - 1+2); i++)
        {
          if (table_distance.PatchTbl[i].Mse > MaxVal)
          {
            Idx = i;
            MaxVal = (float)table_distance.PatchTbl[i].Mse;
          }
        }
      }

      cut_thres = MaxVal;
    }

  }

  insertionSort(table_distance);
}

FORCE_INLINE void find_patches_fcfs_romb(unsigned k_r
  , float threshold
  , tTableDistance & table_distance
  , const NsImage &img
  , const NsParams & bp
  )
{
  const unsigned width = img.width();
  const unsigned height = img.height();
  const unsigned kHW = bp.kHW;
  const unsigned nHW = bp.nHW;
  const unsigned NHW = bp.NHW;
  const unsigned pHW = bp.pHW;

  float cut_thres = threshold;

  static int Mask[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1)];
  static unsigned Size = 0;

  float MaxVal = -1;
  unsigned Idx = 0;
  unsigned Sz = 0;

  if (Size == 0)
  {
    int N = 1;
    for (int di = -N; di <= N; ++di)
    {
      for (int dj = -N; dj <= N; ++dj)
      {

        if (dj == 0 && di == 0)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }

    for (int di = -int(nHW); di <= int(nHW); ++di)
    {
      for (int dj = -int(nHW); dj <= int(nHW); ++dj)
      {
        if (abs(di) + abs(dj) > nHW)
        {
          continue;
        }

        if (abs(dj) <= N && abs(di) <= N)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }
  }

  for (unsigned k = 0; k < Size; k++)
  {

    float sum = (float)patch_distance(img.m_Y, k_r, k_r + Mask[k], kHW, img.width(), (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);

    if (sum > cut_thres)
      continue;

    if (table_distance.Size < (NHW - 1))
    {
      table_distance.PatchTbl[table_distance.Size].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
      table_distance.PatchTbl[table_distance.Size++].Patch = k_r + Mask[k];

      if (sum > MaxVal)
      {
        MaxVal = sum;
        Idx = Sz;
      }
      Sz++;
    }
    else
    {
      if (sum <= MaxVal)
      {
        table_distance.PatchTbl[Idx].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
        table_distance.PatchTbl[Idx].Patch = k_r + Mask[k];

        MaxVal = -1;
        for (unsigned i = 0; i < (NHW - 1); i++)
        {
          if (table_distance.PatchTbl[i].Mse > MaxVal)
          {
            Idx = i;
            MaxVal = (float)table_distance.PatchTbl[i].Mse;
          }
        }
      }

      cut_thres = MaxVal;
    }

  }

  insertionSort(table_distance);

}

FORCE_INLINE void find_patches_fcfs_romb(unsigned k_r
  , float threshold
  , tTableDistance & table_distance
  , const NsImage &img
  , const NsParams & bp
  , int *pSum
  )
{
  const unsigned width = img.width();
  const unsigned height = img.height();
  const unsigned kHW = bp.kHW;
  const unsigned kHW2 = kHW >> 1;
  const unsigned nHW = bp.nHW;
  const unsigned NHW = bp.NHW;
  const unsigned pHW = bp.pHW;

  float cut_thres = threshold;

  static int Mask[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1)];
  static unsigned Size = 0;

  float MaxVal = -1;
  unsigned Idx = 0;
  unsigned Sz = 0;


  if (Size == 0)
  {
    int N = 1;
    for (int di = -N; di <= N; ++di)
    {
      for (int dj = -N; dj <= N; ++dj)
      {

        if (dj == 0 && di == 0)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }

    for (int di = -int(nHW); di <= int(nHW); ++di)
    {
      for (int dj = -int(nHW); dj <= int(nHW); ++dj)
      {
        if (abs(di) + abs(dj) > nHW)
        {
          continue;
        }

        if (abs(dj) <= N && abs(di) <= N)
          continue;

        const int dk = di * width + dj;

        Mask[Size++] = dk;
      }
    }
  }


  for (unsigned k = 0; k < Size; k++)
  {
    float sum = (float)pSum[k];

    if (sum <= 0)
      sum = (float)patch_distance_short(img.m_Y, k_r, k_r + Mask[k], kHW, width, (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);

    pSum[k] = -1;
    if (sum > cut_thres)
      continue;

    pSum[k] = patch_distance_short1(img.m_Y, k_r + kHW2*width, k_r + Mask[k] + kHW2 * width, kHW, width, (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);
    sum += pSum[k];



    if (sum > cut_thres)
      continue;

    if (table_distance.Size < (NHW - 1))
    {
      table_distance.PatchTbl[table_distance.Size].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
      table_distance.PatchTbl[table_distance.Size++].Patch = k_r + Mask[k];

      if (sum > MaxVal)
      {
        MaxVal = sum;
        Idx = Sz;
      }
      Sz++;
    }
    else
    {
      if (sum <= MaxVal)
      {
        table_distance.PatchTbl[Idx].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
        table_distance.PatchTbl[Idx].Patch = k_r + Mask[k];

        MaxVal = -1;
        for (unsigned i = 0; i < (NHW - 1); i++)
        {
          if (table_distance.PatchTbl[i].Mse > MaxVal)
          {
            Idx = i;
            MaxVal = (float)table_distance.PatchTbl[i].Mse;
          }
        }
      }

      cut_thres = MaxVal;
    }

  }

  insertionSort(table_distance);
}

void precompute_BM_ALG2(
  tNs *pNs
  , const NsImage &img
  , float          tauMatch
  )
{
  const unsigned width = img.width();
  const unsigned height = img.height();
  const unsigned kHW = pNs->ns_params.kHW;
  const unsigned nHW = pNs->ns_params.nHW;
  const unsigned NHW = pNs->ns_params.NHW;
  const unsigned pHW = pNs->ns_params.pHW;

  const unsigned Ns = 2 * nHW + 1;
  const float threshold = tauMatch * kHW * kHW;

  int Sum[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1) / 2];

  memset(Sum, 0, sizeof(Sum));

#if RRC_NS_LIFTING

  for (unsigned ind_j = 0; ind_j < pNs->column_ind.Size; ind_j++)
  {
    for (unsigned ind_i = 0; ind_i < pNs->row_ind.Size; ind_i++)
    {
      const unsigned i_r = pNs->row_ind.pPxlNum[ind_i];
      const unsigned j_r = pNs->column_ind.pPxlNum[ind_j];
      const unsigned k_r = i_r * width + j_r;

      pNs->table_distance.Size = 0;

      if (Check(img.m_Y.data(), j_r, i_r, width, kHW))
        find_patches_fcfs_romb(k_r, threshold, pNs->table_distance, img, pNs->ns_params, Sum);
      else
        memset(Sum, 0, sizeof(Sum));

      const unsigned nSx_r = min<unsigned>(closest_power_of_2(pNs->table_distance.Size + 1), NHW);

      //! Keep a maximum of NHW similar patches
      {
        tPatchTable *pPatchTable = &pNs->patch_table[ind_i*pNs->column_ind.Size + ind_j];

        pPatchTable->Patch[pPatchTable->Size++] = pPatchTable->PixNum;
        for (unsigned n = 0; n + 1< nSx_r; n++)
        {
          pPatchTable->Patch[pPatchTable->Size++] = pNs->table_distance.PatchTbl[n].Patch;
        }
      }
    }
    memset(Sum, 0, sizeof(Sum));
  }
#else
  for (unsigned ind_i = 0; ind_i < pNs->row_ind.Size; ind_i++)
  {
    for (unsigned ind_j = 0; ind_j < pNs->column_ind.Size; ind_j++)
    {
      const unsigned i_r = pNs->row_ind.pPxlNum[ind_i];
      const unsigned j_r = pNs->column_ind.pPxlNum[ind_j];
      const unsigned k_r = i_r * width + j_r;

      pNs->table_distance.Size = 0;

      if (Check(img.m_Y.data(), j_r, i_r, width, kHW))
        find_patches_fcfs_romb(k_r, threshold, pNs->table_distance, img, pNs->NsParams);
      else
        memset(Sum, 0, sizeof(Sum));

      const unsigned nSx_r = min<unsigned>(closest_power_of_2(pNs->table_distance.Size + 1), NHW);

      //! Keep a maximum of NHW similar patches

      {
        tPatchTable *pPatchTable = &pNs->patch_table[ind_i*pNs->column_ind.Size + ind_j];

        pPatchTable->Patch[pPatchTable->Size++] = pPatchTable->PixNum;
        for (unsigned n = 0; n + 1< nSx_r; n++)
        {
          pPatchTable->Patch[pPatchTable->Size++] = pNs->table_distance.PatchTbl[n].Patch;
        }
      }

    }
    memset(Sum, 0, sizeof(Sum));
  }
#endif // RRC_NS_LIFTING

}

FORCE_INLINE void find_patches_fcfs_mask(unsigned k_r
  , float threshold
  , tTableDistance & table_distance
  , const NsImage &img
  , const NsParams & bp
  , const int *pMask
  , unsigned Size
  )
{
  const unsigned width = img.width();
  const unsigned kHW = bp.kHW;
  const unsigned kHW2 = kHW >> 1;
  const unsigned NHW = bp.NHW;

  float cut_thres = threshold;
  float MaxVal = -1;
  unsigned Idx = 0;
  unsigned Sz = 0;


  for (int k = 0; k <Size; ++k)
  {
    float sum = (float)patch_distance(img.m_Y, k_r, k_r + pMask[k], kHW, width, (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);

    if (sum > cut_thres)
      continue;

    if (table_distance.Size < (NHW - 1 + 2))
    {
      table_distance.PatchTbl[table_distance.Size].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
      table_distance.PatchTbl[table_distance.Size++].Patch = k_r + pMask[k];

      if (sum > MaxVal)
      {
        MaxVal = sum;
        Idx = Sz;
      }
      Sz++;
    }
    else
    {
      if (sum <= MaxVal)
      {
        table_distance.PatchTbl[Idx].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
        table_distance.PatchTbl[Idx].Patch = k_r + pMask[k];

        MaxVal = -1;
        for (unsigned i = 0; i < (NHW - 1 + 2); i++)
        {
          if (table_distance.PatchTbl[i].Mse > MaxVal)
          {
            Idx = i;
            MaxVal = (float)table_distance.PatchTbl[i].Mse;
          }
        }
      }
      cut_thres = MaxVal;
    }

  }

  insertionSort(table_distance, NHW - 2);
}

FORCE_INLINE void find_patches_fcfs_mask(unsigned k_r
  , float threshold
  , tTableDistance & table_distance
  , const NsImage &img
  , const NsParams & bp
  , const int *pMask
  , unsigned Size
  , int *pSum
  )
{
  const unsigned width = img.width();
  const unsigned kHW = bp.kHW;
  const unsigned kHW2 = kHW >> 1;
  const unsigned NHW = bp.NHW;

  RRC_NS_TABLE_DISTANCE_TYPE cut_thres = (RRC_NS_TABLE_DISTANCE_TYPE)threshold;
  RRC_NS_TABLE_DISTANCE_TYPE MaxVal = -1;
  unsigned Idx = 0;
  unsigned Sz = 0;

  for (int k = 0; k <Size; ++k)
  {
    RRC_NS_TABLE_DISTANCE_TYPE sum = pSum[k];
    //tPd Pd;

    if (sum <= 0)
      sum = patch_distance_short(img.m_Y, k_r, k_r + pMask[k], kHW, width, (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);

    pSum[k] = -1;
    if (sum > cut_thres)
      continue;

    pSum[k] = patch_distance_short1(img.m_Y, k_r + kHW2*width, k_r + pMask[k] + kHW2 * width, kHW, width, (const RRC_NS_TABLE_DISTANCE_TYPE)cut_thres);
    sum += pSum[k];

    if (sum > cut_thres)
      continue;

    //Pd.Mse = sum;
    //Pd.Patch = k_r + pMask[k];

    if (table_distance.Size < (NHW - 1 + 2))
    {
      table_distance.PatchTbl[table_distance.Size].Mse = (RRC_NS_TABLE_DISTANCE_TYPE)sum;
      table_distance.PatchTbl[table_distance.Size++].Patch = k_r + pMask[k];

      //table_distance.PatchTbl[table_distance.Size++] = Pd;

      if (sum > MaxVal)
      {
        MaxVal = sum;
        Idx = Sz;
      }
      Sz++;
    }
    else
    {
      if (sum <= MaxVal)
      {
        table_distance.PatchTbl[Idx].Mse = sum;
        table_distance.PatchTbl[Idx].Patch = k_r + pMask[k];

        //table_distance.PatchTbl[Idx] = Pd;

        MaxVal = -1;
        for (unsigned i = 0; i < (NHW - 1 + 2); i++)
        {
          if (table_distance.PatchTbl[i].Mse > MaxVal)
          {
            Idx = i;
            MaxVal = table_distance.PatchTbl[i].Mse;
          }
        }
      }
      cut_thres = MaxVal;
    }

  }

  insertionSort(table_distance, NHW - 2);
}


void find_patches_sb(int sbj, int sbi, tPatchTable *patch_table, const NsImage &img, const NsParams & p, const float    threshold, int SB_SIZE, const tLattice & column_ind, const tLattice & row_ind, tTableDistance & table_distance)
{
  const int Ns = 2 * p.nHW + 1;

  const int bi_min = sbi * SB_SIZE;
  const int bj_min = sbj * SB_SIZE;
  const int bi_max = min((sbi + 1) * SB_SIZE, (int)row_ind.Size);
  const int bj_max = min((sbj + 1) * SB_SIZE, (int)column_ind.Size);
  const int db_max = min(bi_max - bi_min, bj_max - bj_min);

  const unsigned width = img.width();
  int Mask[(2 * RRC_NS_nHW + 1)*(2 * RRC_NS_nHW + 1)];
#if RRC_NS_LIFTING
  int Sum[RRC_NS_BM_SB_SIZE*(RRC_NS_NHW + 5) / 2];
#endif

  memset(Mask, 0, sizeof(Mask));


  // diagonal
  for (int db = 0; db < db_max; db += 2)
  {
    table_distance.Size = 0;

    const unsigned bi = bi_min + db;
    const unsigned bj = bj_min + db;

    const int k_r = row_ind.pPxlNum[bi] * img.width() + column_ind.pPxlNum[bj];
    if (Check(img.m_Y.data(), column_ind.pPxlNum[bj], row_ind.pPxlNum[bi], img.width(), p.kHW))
      find_patches_fcfs_romb_fst(k_r, threshold, table_distance, img, p);

    for (size_t n = 0; n < table_distance.Size; ++n)
    {
      int pt = table_distance.PatchTbl[n].Patch;
      int di = pt / img.width() - (int)row_ind.pPxlNum[bi];
      int dj = pt % img.width() - (int)column_ind.pPxlNum[bj];

      Mask[(di + p.nHW) * Ns + dj + p.nHW] += 1;
    }

    put(table_distance, patch_table, p.NHW, bi*column_ind.Size + bj);
  }

  unsigned Size = 0;
  for (int k = 0; k < Ns*Ns; ++k)
  {
    if (Mask[k] == 0) continue;

    const int di = k / Ns - p.nHW;
    const int dj = k % Ns - p.nHW;

    if (di == 0 && dj == 0)
      continue;

    const int dk = di * width + dj;
#if RRC_NS_LIFTING
    Sum[Size] = -1;
#endif
    Mask[Size++] = dk;
  }

  //assert(Size < RRC_NS_BM_SB_SIZE*(RRC_NS_NHW + 5) / 2);

#if !RRC_NS_LIFTING
  for (int bi = bi_min; bi < bi_max; ++bi)
  {
    for (int bj = bj_min; bj < bj_max; ++bj)
    {
      const int k_r = row_ind.pPxlNum[bi] * img.width() + column_ind.pPxlNum[bj];

      if (patch_table[bi*column_ind.Size + bj].Size > 0) continue;

      table_distance.Size = 0;

      if (Check(img.m_Y.data(), column_ind.pPxlNum[bj], row_ind.pPxlNum[bi], img.width(), p.kHW))
        find_patches_fcfs_mask(k_r, threshold, table_distance, img, p, Mask, Size);

      put(table_distance, patch_table, p.NHW, bi*column_ind.Size + bj);
    }
  }
#else

  const unsigned SumSize = sizeof(Sum[0])*Size;

  for (int bj = bj_min; bj < bj_max; ++bj)
  {
    for (int bi = bi_min; bi < bi_max; ++bi)
    {
      const int k_r = row_ind.pPxlNum[bi] * img.width() + column_ind.pPxlNum[bj];

      if (patch_table[bi*column_ind.Size + bj].Size > 0)
      {
        memset(Sum, 0, SumSize);
        continue;
      }

      table_distance.Size = 0;

      if (Check(img.m_Y.data(), column_ind.pPxlNum[bj], row_ind.pPxlNum[bi], img.width(), p.kHW))
      {
        find_patches_fcfs_mask(k_r, threshold, table_distance, img, p, Mask, Size, Sum);
      }
      else
      {
        memset(Sum, 0, SumSize);
      }

      put(table_distance, patch_table, p.NHW, bi*column_ind.Size + bj);
    }
    memset(Sum, 0, SumSize);
  }
#endif


}


void precompute_BM_SB_2(
  tNs *pNs
  , const NsImage &img
  , float          tauMatch
  )
{
  const unsigned kHW = pNs->ns_params.kHW;
  const RRC_NS_TABLE_DISTANCE_TYPE threshold = (const RRC_NS_TABLE_DISTANCE_TYPE)(tauMatch * kHW * kHW);

  const int SB_SIZE = RRC_NS_BM_SB_SIZE; // in blocks
  const int WIDTH_IN_SB = (int)ceil(pNs->column_ind.Size / float(SB_SIZE));
  const int HEIGHT_IN_SB = (int)ceil(pNs->row_ind.Size / float(SB_SIZE));

  for (unsigned sbi = 0; sbi < HEIGHT_IN_SB; sbi++)
  {
    for (unsigned sbj = 0; sbj < WIDTH_IN_SB; sbj++)
    {
      find_patches_sb(sbj, sbi, pNs->patch_table, img, pNs->ns_params, (const float)threshold, SB_SIZE, pNs->column_ind, pNs->row_ind, pNs->table_distance);
    }
  }
}

void precompute_BM(
  tNs *pNs
  , const NsImage &img
  , float          tauMatch
  ){
#if RRC_NS_USE_SAD
  tauMatch = sqrt(tauMatch);
#endif

  TOKENPASTE2(precompute_BM_, RRC_NS_BM_ALGO)(
    pNs
    , img
    , tauMatch);
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////                                                           /////////////
/////////////                    WIENER FITERING                        /////////////
/////////////                                                           /////////////
/////////////////////////////////////////////////////////////////////////////////////

#define BIT_RND 4

#define BIT_RND8 (3+BIT_RND)
#define BIT_RND8M1 (BIT_RND8-1)

#define BIT_RND4 (2+BIT_RND)
#define BIT_RND4M1 (BIT_RND4-1)

#define BIT_RND2 (1+BIT_RND)
#define BIT_RND2M1 (BIT_RND2-1)

#define BIT_RND1 (0+BIT_RND)
#define BIT_RND1M1 (BIT_RND1-1)

#define WIENER_TABLE_SIZE (1024)
#define WIENER_TABLE_SIZE8 (WIENER_TABLE_SIZE*8-1)
#define WIENER_TABLE_SIZE4 (WIENER_TABLE_SIZE*4-1)
#define WIENER_TABLE_SIZE2 (WIENER_TABLE_SIZE*2-1)
#define WIENER_TABLE_SIZE1 (WIENER_TABLE_SIZE-1)

#define OFFSET8 (8 * WIENER_TABLE_SIZE - 1)
#define OFFSET4 (4 * WIENER_TABLE_SIZE - 1)
#define OFFSET2 (2 * WIENER_TABLE_SIZE - 1)
#define OFFSET1 (1 * WIENER_TABLE_SIZE - 1)

FORCE_INLINE void fast_wiener_filter_dct8(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const int *pTbl, const int width, const unsigned Size)
{
  const unsigned stride = width - Size;
  short *pCntOld = pCnt;

  const unsigned p0 = pPatch[0];
  const unsigned p1 = pPatch[1];
  const unsigned p2 = pPatch[2];
  const unsigned p3 = pPatch[3];
  const unsigned p4 = pPatch[4];
  const unsigned p5 = pPatch[5];
  const unsigned p6 = pPatch[6];
  const unsigned p7 = pPatch[7];


  for (unsigned r = 0; r < Size; ++r, pIn += stride, pOut += stride, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pIn++, pOut++, pCnt++)
    {
      const int u0 = pIn[p0];
      const int u1 = pIn[p1];
      const int u2 = pIn[p2];
      const int u3 = pIn[p3];

      const int u4 = pIn[p4];
      const int u5 = pIn[p5];
      const int u6 = pIn[p6];
      const int u7 = pIn[p7];

      const int s07 = u0 + u7;
      const int s16 = u1 + u6;
      const int s25 = u2 + u5;
      const int s34 = u3 + u4;
      const int d07 = u0 - u7;
      const int d16 = u1 - u6;
      const int d25 = u2 - u5;
      const int d34 = u3 - u4;

      const int a0 = s07 + s34;
      const int a1 = s16 + s25;
      const int a2 = s07 - s34;
      const int a3 = s16 - s25;
      const int a4 = d16 + d25 + (d07 + (d07 >> 1));
      const int a5 = d07 - d34 - (d25 + (d25 >> 1));
      const int a6 = d07 + d34 - (d16 + (d16 >> 1));
      const int a7 = d16 - d25 + (d34 + (d34 >> 1));

      // filtration

      const int z0 = pTbl[a0 + a1 + OFFSET8];
      const int z1 = pTbl[a4 + (a7 >> 2) + OFFSET8];
      const int z2 = pTbl[a2 + (a3 >> 1) + OFFSET8];
      const int z3 = pTbl[a5 + (a6 >> 2) + OFFSET8];
      const int z4 = pTbl[a0 - a1 + OFFSET8];
      const int z5 = pTbl[a6 - (a5 >> 2) + OFFSET8];
      const int z6 = pTbl[(a2 >> 1) - a3 + OFFSET8];
      const int z7 = pTbl[(a4 >> 2) - a7 + OFFSET8];

      // backward transform

      const int ia0 = z0 + z4;
      const int ia2 = z0 - z4;
      const int ia4 = (z2 >> 1) - z6;
      const int ia6 = (z6 >> 1) + z2;
      const int ia1 = -z3 + z5 - z7 - (z7 >> 1);
      const int ia3 = z1 + z7 - z3 - (z3 >> 1);
      const int ia5 = -z1 + z7 + z5 + (z5 >> 1);
      const int ia7 = z3 + z5 + z1 + (z1 >> 1);

      const int b0 = ia0 + ia6;
      const int b2 = ia2 + ia4;
      const int b4 = ia2 - ia4;
      const int b6 = ia0 - ia6;
      const int b1 = (ia7 >> 2) + ia1;
      const int b3 = ia3 + (ia5 >> 2);
      const int b5 = (ia3 >> 2) - ia5;
      const int b7 = ia7 - (ia1 >> 2);

      pOut[p0] += ((b0 + b7) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p1] += ((b2 + b5) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p2] += ((b4 + b3) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p3] += ((b6 + b1) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p4] += ((b6 - b1) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p5] += ((b4 - b3) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p6] += ((b2 - b5) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p7] += ((b0 - b7) + (1 << BIT_RND8M1)) >> BIT_RND8;

    }
  }

  pCnt = pCntOld;

  for (unsigned r = 0; r < Size; ++r, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pCnt++)
    {

      ++pCnt[p0];
      ++pCnt[p1];
      ++pCnt[p2];
      ++pCnt[p3];
      ++pCnt[p4];
      ++pCnt[p5];
      ++pCnt[p6];
      ++pCnt[p7];

    }
  }
}

FORCE_INLINE void fast_wiener_filter8(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const int *pTbl, const int width, const unsigned Size)
{
  const unsigned stride = width - Size;
  short *pCntOld = pCnt;

  const unsigned p0 = pPatch[0];
  const unsigned p1 = pPatch[1];
  const unsigned p2 = pPatch[2];
  const unsigned p3 = pPatch[3];
  const unsigned p4 = pPatch[4];
  const unsigned p5 = pPatch[5];
  const unsigned p6 = pPatch[6];
  const unsigned p7 = pPatch[7];


  for (unsigned r = 0; r < Size; ++r, pIn += stride, pOut += stride, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pIn++, pOut++, pCnt++)
    {
      const int u0 = pIn[p0];
      const int u1 = pIn[p1];
      const int u2 = pIn[p2];
      const int u3 = pIn[p3];

      const int u4 = pIn[p4];
      const int u5 = pIn[p5];
      const int u6 = pIn[p6];
      const int u7 = pIn[p7];

      // forward transform
      const int x0 = u0 + u4;
      const int x1 = u1 + u5;
      const int x2 = u2 + u6;
      const int x3 = u3 + u7;

      const int x4 = u0 - u4;
      const int x5 = u1 - u5;
      const int x6 = u2 - u6;
      const int x7 = u3 - u7;

      const int y0 = x0 + x2;
      const int y1 = x1 + x3;
      const int y2 = x0 - x2;
      const int y3 = x1 - x3;

      const int y4 = x4 + x6;
      const int y5 = x5 + x7;
      const int y6 = x4 - x6;
      const int y7 = x5 - x7;

      // filtration

      const int z0 = pTbl[y0 + y1 + OFFSET8];
      const int z1 = pTbl[y0 - y1 + OFFSET8];
      const int z2 = pTbl[y2 + y3 + OFFSET8];
      const int z3 = pTbl[y2 - y3 + OFFSET8];

      const int z4 = pTbl[y4 + y5 + OFFSET8];
      const int z5 = pTbl[y4 - y5 + OFFSET8];
      const int z6 = pTbl[y6 + y7 + OFFSET8];
      const int z7 = pTbl[y6 - y7 + OFFSET8];

      // backward transform

      const int ix0 = z0 + z4;
      const int ix1 = z1 + z5;
      const int ix2 = z2 + z6;
      const int ix3 = z3 + z7;

      const int ix4 = z0 - z4;
      const int ix5 = z1 - z5;
      const int ix6 = z2 - z6;
      const int ix7 = z3 - z7;

      const int iy0 = ix0 + ix2;
      const int iy1 = ix1 + ix3;
      const int iy2 = ix0 - ix2;
      const int iy3 = ix1 - ix3;

      const int iy4 = ix4 + ix6;
      const int iy5 = ix5 + ix7;
      const int iy6 = ix4 - ix6;
      const int iy7 = ix5 - ix7;


      pOut[p0] += ((iy0 + iy1) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p1] += ((iy0 - iy1) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p2] += ((iy2 + iy3) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p3] += ((iy2 - iy3) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p4] += ((iy4 + iy5) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p5] += ((iy4 - iy5) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p6] += ((iy6 + iy7) + (1 << BIT_RND8M1)) >> BIT_RND8;
      pOut[p7] += ((iy6 - iy7) + (1 << BIT_RND8M1)) >> BIT_RND8;

          }
    }

  pCnt = pCntOld;

  for (unsigned r = 0; r < Size; ++r, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pCnt++)
    {

      ++pCnt[p0];
      ++pCnt[p1];
      ++pCnt[p2];
      ++pCnt[p3];
      ++pCnt[p4];
      ++pCnt[p5];
      ++pCnt[p6];
      ++pCnt[p7];

    }
  }
}


FORCE_INLINE void fast_wiener_filter4(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const int *pTbl, const int width, const unsigned Size)
{
  const unsigned stride = width - Size;

  short *pCntOld = pCnt;

  const unsigned p0 = pPatch[0];
  const unsigned p1 = pPatch[1];
  const unsigned p2 = pPatch[2];
  const unsigned p3 = pPatch[3];

  for (unsigned r = 0; r < Size; ++r, pIn += stride, pOut += stride, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pIn++, pOut++, pCnt++)
    {
      const int x0 = pIn[p0];
      const int x1 = pIn[p1];
      const int x2 = pIn[p2];
      const int x3 = pIn[p3];

      // forward transform
      const int y0 = x0 + x2;
      const int y1 = x1 + x3;
      const int y2 = x0 - x2;
      const int y3 = x1 - x3;

      // filtration

      const int z0 = pTbl[y0 + y1 + OFFSET4];
      const int z1 = pTbl[y0 - y1 + OFFSET4];
      const int z2 = pTbl[y2 + y3 + OFFSET4];
      const int z3 = pTbl[y2 - y3 + OFFSET4];

      // backward transform

      const int iy0 = z0 + z2;
      const int iy1 = z1 + z3;
      const int iy2 = z0 - z2;
      const int iy3 = z1 - z3;


      pOut[p0] += ((iy0 + iy1) + (1 << BIT_RND4M1)) >> BIT_RND4;
      pOut[p1] += ((iy0 - iy1) + (1 << BIT_RND4M1)) >> BIT_RND4;
      pOut[p2] += ((iy2 + iy3) + (1 << BIT_RND4M1)) >> BIT_RND4;
      pOut[p3] += ((iy2 - iy3) + (1 << BIT_RND4M1)) >> BIT_RND4;
    }
  }

  pCnt = pCntOld;
  for (unsigned r = 0; r < Size; ++r, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pCnt++)
    {
      ++pCnt[p0];
      ++pCnt[p1];
      ++pCnt[p2];
      ++pCnt[p3];
    }
  }

}

FORCE_INLINE void fast_wiener_filter2(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const int *pTbl, const int width, const unsigned Size)
{

  const unsigned stride = width - Size;

  const unsigned p0 = pPatch[0];
  const unsigned p1 = pPatch[1];

  for (unsigned r = 0; r < Size; ++r, pIn += stride, pOut += stride, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pIn++, pOut++, pCnt++)
    {

      // forward transform
      const int y0 = pIn[p0];
      const int y1 = pIn[p1];

      // filtration
      const int iy0 = pTbl[y0 + y1 + OFFSET2];
      const int iy1 = pTbl[y0 - y1 + OFFSET2];

      // backward transform

      pOut[p0] += ((iy0 + iy1) + (1 << BIT_RND2M1)) >> BIT_RND2;
      ++pCnt[p0];
      pOut[p1] += ((iy0 - iy1) + (1 << BIT_RND2M1)) >> BIT_RND2;
      ++pCnt[p1];
    }
  }
}

FORCE_INLINE void fast_wiener_filter1(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const int *pTbl, const int width, const unsigned Size)
{

  const unsigned stride = width - Size;

  for (unsigned r = 0; r < Size; ++r, pIn += stride, pOut += stride, pCnt += stride)
  {
    for (unsigned c = 0; c < Size; ++c, pIn++, pOut++, pCnt++)
    {
      pOut[pPatch[0]] += (int)pIn[pPatch[0]] << 10;
      pCnt[pPatch[0]]++;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////                                                           /////////////
/////////////                CORE IMPLEMENTATION                        /////////////
/////////////                                                           /////////////
/////////////////////////////////////////////////////////////////////////////////////


int run_ns(
  const float sigma
  , NsImage &img_noisy
  , NsImage &img_denoised
  , const int bit_depth
  )
{
  static tNs *pNs = 0;

  if (pNs == 0)
  {
    NsParams NsPrm;

    NsPrm.nHW = RRC_NS_nHW; //! Half size of the search window
    NsPrm.NHW = RRC_NS_NHW; //! Must be a power of 2
    NsPrm.bit_depth = bit_depth;

    NsPrm.width = img_noisy.width();
    NsPrm.height = img_noisy.height();

    NsPrm.width_chroma = img_noisy.width_chroma();
    NsPrm.height_chroma = img_noisy.height_chroma();

    if (NsPrm.height > 480)
    {
      NsPrm.kHW = 8;
#if  RRC_NS_LIFTING
      NsPrm.pHW = 4;
#else
      NsPrm.pHW = 3;
#endif
      //NsPrm.pHW = 2;
    }
    else
    {
      NsPrm.kHW = 4;
      NsPrm.pHW = 2;
    }

    pNs = (tNs *)NsCreate(&NsPrm);
  }

  {
    //! Add boundaries and symetrize them
    const unsigned w_b = pNs->ns_params.width + 2 * pNs->ns_params.nHW;
    const unsigned w_b_c = pNs->ns_params.width_chroma + pNs->ns_params.nHW;

    symetrize(img_noisy.m_Y, pNs->img_sym_noisy.m_Y, img_noisy.width(), img_noisy.height(), pNs->ns_params.nHW);

    symetrize(img_noisy.m_U, pNs->img_sym_noisy.m_U, img_noisy.width_chroma(), img_noisy.height_chroma(), pNs->ns_params.nHW / 2);
    symetrize(img_noisy.m_V, pNs->img_sym_noisy.m_V, img_noisy.width_chroma(), img_noisy.height_chroma(), pNs->ns_params.nHW / 2);

    NsBlockMaching(pNs, sigma, pNs->img_sym_noisy);
    NsFltLuma(pNs, sigma, pNs->img_sym_noisy);
#if RRC_NS_CHROMA
    NsFltChroma(pNs, sigma, pNs->img_sym_noisy);
#endif
#if  RRC_NS_ITER
    NsFltLuma(pNs, sigma, pNs->img_sym_noisy);
#if RRC_NS_CHROMA
    NsFltChroma(pNs, sigma, pNs->img_sym_noisy);
#endif
#endif

    const unsigned dc_b = pNs->ns_params.nHW * w_b + pNs->ns_params.nHW;
    unsigned dc = 0;
    for (unsigned i = 0; i < pNs->ns_params.height; i++)
      for (unsigned j = 0; j < pNs->ns_params.width; j++, dc++)
        img_denoised.m_Y[dc] = pNs->img_sym_noisy.m_Y[dc_b + i * w_b + j];

#if RRC_NS_CHROMA
    {
      const unsigned dc_b = pNs->ns_params.nHW * w_b_c / 2 + pNs->ns_params.nHW / 2;
      unsigned dc = 0;
      for (unsigned i = 0; i < pNs->ns_params.height_chroma; i++)
      {
        for (unsigned j = 0; j < pNs->ns_params.width_chroma; j++, dc++)
        {
          img_denoised.m_U[dc] = pNs->img_sym_noisy.m_U[dc_b + i * w_b_c + j];
          img_denoised.m_V[dc] = pNs->img_sym_noisy.m_V[dc_b + i * w_b_c + j];
        }
      }
    }
#endif

  }
  return EXIT_SUCCESS;
}


void *NsCreate(NsParams *pNsParams)
{
  tNs *pNs = new tNs;
  unsigned width = pNsParams->width + 2 * pNsParams->nHW;
  unsigned height = pNsParams->height + 2 * pNsParams->nHW;

  unsigned width_chroma = pNsParams->width_chroma + pNsParams->nHW;
  unsigned height_chroma = pNsParams->height_chroma + pNsParams->nHW;

  memcpy(&pNs->ns_params, pNsParams, sizeof(*pNsParams));

  ind_initialize(pNs->row_ind, height - pNsParams->kHW + 1, pNsParams->nHW, pNsParams->pHW);
  ind_initialize(pNs->column_ind, width - pNsParams->kHW + 1, pNsParams->nHW, pNsParams->pHW);

  pNs->patch_table = new tPatchTable[pNs->row_ind.Size * pNs->column_ind.Size];
  memset(pNs->patch_table, 0, sizeof(pNs->patch_table[0]) * pNs->row_ind.Size * pNs->column_ind.Size);

  pNs->denominator = new short[width * height];
  pNs->numerator = new int[width * height];

  pNs->img_sym_noisy = NsImage(width, height, width_chroma, height_chroma);

  return pNs;
}

void NsBlockMaching(
  tNs *pNs
  , const float sigma
  , const NsImage & img_noisy
  )
{
  const unsigned width = pNs->ns_params.width + 2 * pNs->ns_params.nHW;
  const unsigned height = pNs->ns_params.height + 2 * pNs->ns_params.nHW;
  const unsigned kWien = pNs->ns_params.kHW;

  //! Parameters initialization
  float tauMatch = (sigma < 35.0f ? 400.f : 3500.f); //! threshold used to determinate similarity between patches
  tauMatch *= SQR(1 << (pNs->ns_params.bit_depth - 8));

  // Initialization
  {
    for (unsigned ind_i = 0; ind_i < pNs->row_ind.Size; ind_i++)
    {
      const unsigned i_r = pNs->row_ind.pPxlNum[ind_i];
      for (unsigned ind_j = 0; ind_j < pNs->column_ind.Size; ind_j++)
      {
        const unsigned j_r = pNs->column_ind.pPxlNum[ind_j];
        const unsigned k_r = i_r * width + j_r;
        const unsigned idx = ind_i * (unsigned)pNs->column_ind.Size + ind_j;

        pNs->patch_table[idx].Size = 0;
        pNs->patch_table[idx].PixNum = k_r;
      }
    }
  }

  // Block matching
  precompute_BM(pNs, img_noisy, tauMatch);
}

void NsFltLuma(
  tNs *pNs
  , const float sigma
  , NsImage & img_noisy
  )
{
  const unsigned width = pNs->ns_params.width + 2 * pNs->ns_params.nHW;
  const unsigned height = pNs->ns_params.height + 2 * pNs->ns_params.nHW;
  const unsigned kWien = pNs->ns_params.kHW;

  // Initialization
  {
    memset(pNs->denominator, 0, width * height * sizeof(pNs->denominator[0]));
    memset(pNs->numerator, 0, width * height * sizeof(pNs->numerator[0]));
  }

  for (unsigned ind_i = 0; ind_i < pNs->row_ind.Size; ind_i++)
  {
    for (unsigned ind_j = 0; ind_j < pNs->column_ind.Size; ind_j++)
    {
      const unsigned k_r = ind_i * (unsigned)pNs->column_ind.Size + ind_j;

      //! Number of similar patches
      const size_t nSx_r = pNs->patch_table[k_r].Size;

      wiener_filtering_hadamard4(img_noisy.m_Y.data(), pNs->numerator, pNs->denominator, pNs->patch_table[k_r].Patch, nSx_r, kWien, sigma, width);
    }
  }

  //! Final reconstruction
  const int MaxVal = (int)((1 << pNs->ns_params.bit_depth) - 1);
  for (unsigned k = 0; k < width * height; k++)
  {
    int Pxl;

    pNs->denominator[k] = max<short>(1, pNs->denominator[k]);

    Pxl = pNs->numerator[k] / pNs->denominator[k];
    Pxl += (1 << 9);
    Pxl >>= 10;
    Pxl = max(0, Pxl);
    Pxl = min(MaxVal, Pxl);
    img_noisy.m_Y[k] = (RRC_NS_PEL_TYPE)Pxl;
  }
}
void NsFltChroma(
  tNs *pNs
  , const float sigma
  , NsImage & img_noisy
  )
{
  const unsigned width = pNs->ns_params.width + 2 * pNs->ns_params.nHW;

  const unsigned width_chroma = pNs->ns_params.width_chroma + pNs->ns_params.nHW;
  const unsigned height_chroma = pNs->ns_params.height_chroma + pNs->ns_params.nHW;
  const unsigned kWien = pNs->ns_params.kHW / 2;

  short *denominatorU = pNs->denominator;
  int *numeratorU = pNs->numerator;

  short *denominatorV = &pNs->denominator[width_chroma * height_chroma];
  int *numeratorV = &pNs->numerator[width_chroma * height_chroma];

  memset(denominatorU, 0, width_chroma * height_chroma * sizeof(pNs->denominator[0]));
  memset(numeratorU, 0, width_chroma * height_chroma * sizeof(pNs->numerator[0]));

  memset(denominatorV, 0, width_chroma * height_chroma * sizeof(pNs->denominator[0]));
  memset(numeratorV, 0, width_chroma * height_chroma * sizeof(pNs->numerator[0]));


  for (unsigned ind_i = 0; ind_i < pNs->row_ind.Size; ind_i++)
  {
    for (unsigned ind_j = 0; ind_j < pNs->column_ind.Size; ind_j++)
    {
      const unsigned k_r = ind_i * (unsigned)pNs->column_ind.Size + ind_j;

      //! Number of similar patches
      const size_t nSx_r = pNs->patch_table[k_r].Size;

      unsigned PatchChroma[RRC_NS_nHW];
      unsigned row, col;

      for (unsigned i = 0; i < nSx_r; i++)
      {
        row = (pNs->patch_table[k_r].Patch[i] / width) / 2;
        col = (pNs->patch_table[k_r].Patch[i] % width) / 2;
        PatchChroma[i] = row*width_chroma + col;
      }

      wiener_filtering_hadamard4(img_noisy.m_U.data(), numeratorU, denominatorU, PatchChroma, nSx_r, kWien, sigma, width_chroma);
      wiener_filtering_hadamard4(img_noisy.m_V.data(), numeratorV, denominatorV, PatchChroma, nSx_r, kWien, sigma, width_chroma);
    }
  }

  //! Final reconstruction
  const int MaxVal = (int)((1 << pNs->ns_params.bit_depth) - 1);
  for (unsigned k = 0; k < width_chroma * height_chroma; k++)
  {
    int Pxl;

    denominatorU[k] = max<short>(1, denominatorU[k]);
    Pxl = numeratorU[k] / denominatorU[k];
    Pxl += (1 << 9);
    Pxl >>= 10;
    Pxl = max(0, Pxl);
    Pxl = min(MaxVal, Pxl);
    img_noisy.m_U[k] = (RRC_NS_PEL_TYPE)Pxl;

    denominatorV[k] = max<short>(1, denominatorV[k]);
    Pxl = numeratorV[k] / denominatorV[k];
    Pxl += (1 << 9);
    Pxl >>= 10;
    Pxl = max(0, Pxl);
    Pxl = min(MaxVal, Pxl);
    img_noisy.m_V[k] = (RRC_NS_PEL_TYPE)Pxl;

  }
}


void wiener_filtering_hadamard4(const short *pIn, int *pOut, short *pCnt, const unsigned *pPatch, const size_t nSx_r, const unsigned kWien,float sigma,const int width)
{
  const unsigned kWien_2 = kWien * kWien;
  float value;

  static float SigTbl = 0;
  static int WienerTable8[2 * WIENER_TABLE_SIZE8 + 1];
  static int WienerTable4[2 * WIENER_TABLE_SIZE4 + 1];
  static int WienerTable2[2 * WIENER_TABLE_SIZE2 + 1];
  static int WienerTable1[2 * WIENER_TABLE_SIZE1 + 1];
  const int *pWienerTable;

  if (SigTbl != sigma)
  {
    int p, m;
    float loc_sig, loc_coef;
    int *pWienerTable;
    int TblSize;

    for (m = 1; m <= 8; m *= 2)
    {
      TblSize = m * WIENER_TABLE_SIZE - 1;

      switch (m)
      {
      case 1: pWienerTable = &WienerTable1[0]; break;
      case 2: pWienerTable = &WienerTable2[0]; break;
      case 4: pWienerTable = &WienerTable4[0]; break;
      case 8: pWienerTable = &WienerTable8[0]; break;
      default: assert(0);
        break;
      }

      loc_sig = sigma * sigma * (float)m;
      loc_coef = 1.0f / (float)m;

      p = 0;
      for (int k = TblSize; k >= 0; k--)
      {
        value = (float)k;
        value *= value;
        value /= (value + loc_sig);
        value = (float)k * value;

        const int value_int = (int)(value + 0.5f) * (1 << (10 + 4));

        pWienerTable[p] = -value_int;
        pWienerTable[2 * TblSize - p] = value_int;

        p++;
      }
    }

    SigTbl = sigma;
  }

  switch (nSx_r)
  {
  case 1: pWienerTable = &WienerTable1[0];
    fast_wiener_filter1(pIn, pOut, pCnt, pPatch, pWienerTable, width, kWien);
    break;
  case 2: pWienerTable = &WienerTable2[0];
    fast_wiener_filter2(pIn, pOut, pCnt, pPatch, pWienerTable, width, kWien);
    break;
  case 4: pWienerTable = &WienerTable4[0];
    fast_wiener_filter4(pIn, pOut, pCnt, pPatch, pWienerTable, width, kWien);
    break;
  case 8: pWienerTable = &WienerTable8[0];
    fast_wiener_filter8(pIn, pOut, pCnt, pPatch, pWienerTable, width, kWien);
    break;
  default: assert(0);
    break;
  }
}

#endif
