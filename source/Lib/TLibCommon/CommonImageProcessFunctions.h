#ifndef __COMMON_IMAGE_PROCESS_FUNCTIONS__
#define __COMMON_IMAGE_PROCESS_FUNCTIONS__

#include "CommonDef.h"
#include <assert.h>

//histogram function
Void createSEHistogram(Int*& piHist, Int iNumValue);
Int* createSEHistogram(Int iNumValue);
Void deleteSEHistogram(Int*& piHist);
Void resetSEHistogram(Int* piHist, Int iNumValue, Int iResetValue=0);
Void addToSEHistogram(Int* piHist, Int iEntrance, Bool bCondition,    Int iAddValue=1);
Void showSEHistogram( Int* piHist, Int iNumValue, Bool bOutputToFile, FILE* fpLogFile, Char* pcLabelName, Int iLabelNum, int iEntranceOffset, Bool bShowZeroEntry = false);
Void showSEHistogramSlim(Int* piHist, Int iNumValue, Bool bOutputToFile, FILE* fpLogFile, Char* pcLabelName, Int iLabelNum, int iEntranceOffset, bool bShowZeroEntry = false);
Void HistogramOneArray(Int* piHist, Int* piArray, Int iNumElementInArray, Int* piMask,  Char* pcLogFileName, int iPOC);
int remove_masked_elements(int* piArrayOne, int* piMask, int iNumElement, int* piArrayOneNew);

//log file operation function
Void createLogFile(FILE*& fpLogFile, Char* pcLogFileName);
FILE* createLogFile(Char* pcLogFileName);
Void closeLogFile(FILE*& fpLogFile);
//Map file operation function
FILE* createMapFile(int width, int height, UChar*& pucMap, Char* pcFileName);
Void  closeMapFile(FILE*& fpMapFile, UChar*& pucMap);
Void DrawLine(UChar* pucMap, int x_start, int y_start, int x_end, int y_end, int width, int height, UChar LineValue, int mode);

int FindArrayMax(Int* piArray, Int iNumElementInArray);
int FindArrayMin(Int* piArray, Int iNumElementInArray);
void FindArrayMax(Int* piArray, Int iNumElementInArray, Int& max_val, Int& max_idx);
void FindArrayMin(Int* piArray, Int iNumElementInArray, Int& min_val, Int& min_idx);
Void FindArrayMaxMin(Int* piArray, Int iNumElementInArray, int& max_value, int& min_value);

void outputYUV(int width, int height, Pel* pPelY, int stride, FILE* fpFile);
void outputYUV(int width, int height, UChar* pPelY, int stride, FILE* fpFile, UChar* pPelU = NULL, UChar* pPelV = NULL);


//debuging function using CImg Lib function
double get_correlation(int* piArrayOne, int* piArrayTwo, int iNumElement);

//entropy calculation
void calcImageEntropy(double* Img_entropy, Pel* pOrg, int iWidth, int iHeight, int iStride, bool bShowEntropyImg = false);
double entropy(int* piArray, int iNumElementInArray);

int getLog2Value(int inputValue);


//template functions
template <typename T> 
double calcSimilarity(T* pBlk, T* pBlk2, int blk_idx, int blk_w, int blk_h)
{
  double simi = 0;

  int sum_diff = 0;
  int step = 0;
  T* pBlk_curr = NULL;
  T* pBlk_neig = NULL;
  if(blk_idx==4)
  {
    step = blk_w;
    pBlk_curr = pBlk  + blk_w - 1;
    pBlk_neig = pBlk2;
  }
  else if(blk_idx == 6)
  {
    step = 1;
    pBlk_curr = pBlk + blk_w * (blk_h-1);
    pBlk_neig = pBlk;
  }
  else
  {
    simi = 0;
    return simi;
  }

  //check difference
  T max_val = 0;
  T min_val = 255;
  for(int i=0;i<blk_w;i++)
  {
    sum_diff +=abs(*pBlk_curr - *pBlk_neig);
    if(*pBlk_curr > max_val)
      max_val = *pBlk_curr;
    if(*pBlk_curr < min_val)
      min_val = *pBlk_curr;
    pBlk_curr+=step;
    pBlk_neig+=step;
  }
  sum_diff = sum_diff/blk_w;

  if(sum_diff<=10.0 && max_val-min_val > 10)
    simi = 1.0 - sum_diff/10.0;

  return simi;
}

//get the average luminance of a block
//Note: this function is not accurate.
template <typename T>
int get_block_avg_weighted(T *buff, int block_w, int block_h, UChar *wMap, int norm)
{
  int sum = 0;
  int i;
  int block_size = block_w * block_h;

  for(i=0;i<block_size;i++)
    sum += (int)buff[i]*(int)wMap[i];

  return sum = (sum + (norm>>1)) / norm;
}

//get the variance luminance of a block
//Note: this function is not accurate
template <typename T>
int get_block_std_simp_weighted(T *buff, int block_w, int block_h, int average, UChar * wMap, int norm)
{
  int sum = 0;
  int i;
  int block_size = block_w*block_h;

  //not using CImg Lib
  for(i=0;i<block_size;i++)
    sum += abs( (buff[i]-average)*(int)wMap[i] );

  return sum = (sum + (norm>>1)) / norm;
}

//get the average luminance of a block
//Note: this function is not accurate.
template <typename T>
int get_block_avg(T *buff, int block_w, int block_h)
{
  int sum = 0;
  int i;
  int block_size = block_w * block_h;

  for(i=0;i<block_size;i++)
    sum += (int)buff[i];

  return ( sum+(block_size>>1) ) / block_size;
}

//get the variance luminance of a block
//Note: this function is not accurate
template <typename T>
int get_block_std_simp(T *buff, int block_w, int block_h, int average)
{
  int sum = 0;
  int i;
  int block_size = block_w*block_h;

  for(i=0;i<block_size;i++)
    sum += abs( buff[i]-average );

  return ( sum+(block_size>>1) ) / block_size;
}

template <typename T> 
void assign_block2Pic(T* pic, T val, int width, int height, int x, int y, int block_w, int block_h)
{
  int i,j,ind_p,k,posx,posy;
  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    if(posy>=height)
      break;
    k = posy * width;
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      if(posx>=width)
        break;
      ind_p = k + posx;
      pic[ind_p] = val;
    }
  }
}

template <typename T> 
void add_block2Pic(T* pic, T val, int width, int height, int x, int y, int block_w, int block_h, T min_val, T max_val, int sign)
{
  int i,j,ind_p,k,posx,posy;

  if(sign<0)
  {
    for(j=0; j<block_h; j++)
    {
      posy = y + j;
      if(posy>=height)
        break;
      k = posy * width;
      for(i=0; i<block_w; i++)
      {
        posx = x + i;
        if(posx>=width)
          break;
        ind_p = k + posx;
        pic[ind_p] = Clip3(min_val,max_val,pic[ind_p]-val);
      }
    }
    return;
  }

  //if sign > 0, i.e., val is positive
  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    if(posy>=height)
      break;
    k = posy * width;
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      if(posx>=width)
        break;
      ind_p = k + posx;
      pic[ind_p] = Clip3(min_val,max_val,pic[ind_p]+val);
    }
  }
}


//fetch a block of pixels (in one plane) from a picture, given the up-left corner of the block (x,y) in the picture and the block size
//Pel* pic,  picture buffer of one plane
//UChar* block, block buffer
//int width,  image width without padding
//int height, image height
//int x,      the x coordinate of the upper-left corner of the block
//int y,      the y --
//int block_w, width of the block
//int block_h, height of --
//int stride  padded line width
template <typename T1, typename T2>
void fetch_block(T1* pic, T2* block, int width, int height, int x, int y, int block_w, int block_h, int stride)
{
  int i,j, ind_p, ind_b, k, kk, posx, posy;

  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    posy = posy>=height?(height-1):posy;
    k    = posy * stride;
    kk   = j * block_w; 
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      posx = posx>=width?(width-1):posx;
      ind_p = k  + posx;
      ind_b = kk + i;
      block[ind_b] = (T2)pic[ind_p];
      assert(block[ind_b]>=0 && block[ind_b]<=255);
    }
  }
}

template <typename T>
void fetch_block(T* pic, T* block, int width, int height, int x, int y, int block_w, int block_h, int stride)
{
  int i,j, ind_p, ind_b, k, kk, posx, posy;

  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    posy = posy>=height?(height-1):posy;
    k    = posy * stride;
    kk   = j * block_w; 
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      posx = posx>=width?(width-1):posx;
      ind_p = k  + posx;
      ind_b = kk + i;
      block[ind_b] = pic[ind_p];
      assert(block[ind_b]>=0 && block[ind_b]<=255);
    }
  }
}

//write a block of pixels (Y component) to a picture
//UChar* pic,   picture buffer
//UChar* block, block buffer
//int width,  image width without padding
//int height, image height
//int x,      the x coordinate of the upper-left corner of the block
//int y,      the y --
//int block_w, width of the block
//int block_h, height of --
//int stride  padded line width
template <typename T>
void write_block(T* pic, T* block, int width, int height, int x, int y, int block_w, int block_h)
{
  int i,j, ind_p, ind_b, k, kk, posx, posy;

  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    posy = posy>=height?(height-1):posy;
    k    = posy * width;
    kk   = j * block_w; 
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      posx = posx>=width?(width-1):posx;
      ind_p = k  + posx;
      ind_b = kk + i;
      pic[ind_p] = block[ind_b];
    }
  }
}

//assign a block region in the picture with the same value
//UChar* pic, picture buffer
//UChar val.  value assigned to the block region in the picture
//int width,  image width without padding
//int height, image height
//int x,      the x coordinate of the upper-left corner of the block
//int y,      the y --
//int block_w, width of the block
//int block_h, height of --
//int stride  padded line width
template <typename T>
void assign_block(T* pic, T val, int width, int height, int x, int y, int block_w, int block_h)
{
  int i,j,ind_p,k,posx,posy;

  for(j=0; j<block_h; j++)
  {
    posy = y + j;
    if(posy>=height)
      break;
    k = posy * width;
    for(i=0; i<block_w; i++)
    {
      posx = x + i;
      if(posx>=width)
        break;
      ind_p = k + posx;
      pic[ind_p] = val;
    }
  }
}

template <typename T>
void composite_block(T* pDes, int strideD, T* pBlk1, int width1, int height1, int stride1, T* pBlk2, int width2, int height2, int stride2, int x, int y)
{
  int i,j;
  T* pTB  = pDes;
  T* pSB1 = pBlk1;
  T* pSB2 = pBlk2;

  for(j=0; j<height1; j++)
  {
    for(i=0; i<width1; i++)
    {
      *(pTB+i) = *(pSB1+i);
    }
    pTB  += strideD;
    pSB1 += stride1;
  }

  pTB = pDes + y*strideD + x;
  for(j=0; j<height2; j++)
  {
    for(i=0; i<width2; i++)
    {
      *(pTB+i) = *(pSB2+i);
    }
    pTB  += strideD;
    pSB2 += stride2;
  }
}

#endif //__COMMON_IMAGE_PROCESS_FUNCTIONS__