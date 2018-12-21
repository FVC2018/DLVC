#include "TLibCommon/CommonImageProcessFunctions.h"
#include "math.h"

using namespace std;

/*--------------------------------------------- output pixel memory to a YUV file ------------------------------------*/

//output the Y component of an image to a YUV 4:2:0 file
//int width,  image width without padding
//int height, image height
//int stride  padded line width
//Pel* pPelY, the pixel memory
//FILE* fpFile, the output file handle
void outputYUV(int width, int height, Pel* pPelY, int stride, FILE* fpFile)
{
  int s  = width*height;
  int hs = s/2;
  UChar* UV  = (UChar*)malloc(hs*sizeof(UChar));
  int i,j;
  Pel* pPic = pPelY;

  for(j=0;j<height;j++)
  {
    for(i=0;i<width;i++)
      UV[i] = (UChar)pPic[i];
    fwrite(UV, sizeof(UChar), width,  fpFile);
    pPic += stride;
  }

  //UV of 128
  memset(UV, 128, hs);
  fwrite(UV, sizeof(UChar), hs, fpFile);

  //free memory
  if(UV)
    free(UV);
}

//output the Y U V components (if exist) of an image to a YUV 4:2:0 file
void outputYUV(int width, int height, UChar* pPelY, int stride, FILE* fpFile, UChar* pPelU, UChar* pPelV)
{
  int s  = width*height;
  int qs = s/4;
  UChar* UV  = NULL;
  int j;
  UChar* pPic = pPelY;
  UV = (UChar*)malloc(s*sizeof(UChar));
  memset(UV, 128, s);

  if(pPelY == NULL)
    fwrite(UV, sizeof(UChar), s, fpFile);
  else
  {
    for(j=0;j<height;j++)
    {
      fwrite(pPic, sizeof(UChar), width,  fpFile);
      pPic += stride;
    }
  }

  if(pPelU == NULL)
    fwrite(UV, sizeof(UChar), qs, fpFile);
  else
    fwrite(pPelU, sizeof(UChar), qs, fpFile);

  if(pPelV == NULL)
    fwrite(UV, sizeof(UChar), qs, fpFile);
  else
    fwrite(pPelV, sizeof(UChar), qs, fpFile);

  //free memory
  if(UV)
    free(UV);
}



/*--------------------------------------------- histogram functions  ------------------------------------*/
//create histogram of an (INT-type) array recording Syntax Element values
//Int* piHist, histogram 
//Int* piArray,  the target array
//Int iNumElementInArray,  length of the array
//Int* piMask,  the mask that records whether each array element should be added to the histogram
//Char* pcLogFileName, log file for recording the histogram
//int iPOC, POC value
Void HistogramOneArray(Int* piHist, Int* piArray, Int iNumElementInArray, Int* piMask,  Char* pcLogFileName, int iPOC)
{
  Int iNumValueInArray = 0;
  Int max_value = 0;
  Int min_value = 0; 
  FindArrayMaxMin(piArray, iNumElementInArray, max_value, min_value);
  iNumValueInArray = max_value - min_value + 1;
  createSEHistogram(piHist,iNumValueInArray);
  resetSEHistogram(piHist, iNumValueInArray);

  if(piMask == NULL)
  {
    for(int i=0; i<iNumElementInArray; i++)
      addToSEHistogram(piHist, piArray[i]-min_value, true);
  }
  else
  {
    for(int i=0; i<iNumElementInArray; i++)
      addToSEHistogram(piHist, piArray[i]-min_value, piMask[i]!=0);
  }

  FILE* fpLogFile = createLogFile(pcLogFileName);

  showSEHistogramSlim(piHist, iNumValueInArray, 0, fpLogFile, "POC", iPOC, min_value, true);

  closeLogFile(fpLogFile);
  deleteSEHistogram(piHist);
}

//Find the max, min values in an array, which is used to find out the dynamic range of the array
int FindArrayMax(Int* piArray, Int iNumElementInArray)
{
  int max_value = -65535;
  int i;
  for(i=0; i< iNumElementInArray; i++)
  {
    if(piArray[i]>max_value)
      max_value = piArray[i];
  }
  return max_value;
}
int FindArrayMin(Int* piArray, Int iNumElementInArray)
{
  int min_value = 65535;
  int i;
  for(i=0; i< iNumElementInArray; i++)
  {
    if(piArray[i]<min_value)
      min_value = piArray[i];
  }
  return min_value;
}
void FindArrayMax(Int* piArray, Int iNumElementInArray, Int& max_val, Int& max_idx)
{
  max_val = -65535;
  int i;
  for(i=0; i< iNumElementInArray; i++)
  {
    if(piArray[i]>max_val)
    {
      max_val = piArray[i];
      max_idx   = i;
    }
  }
}
void FindArrayMin(Int* piArray, Int iNumElementInArray, Int& min_val, Int& min_idx)
{
  min_val = 65535;
  int i;
  for(i=0; i< iNumElementInArray; i++)
  {
    if(piArray[i]<min_val)
    {
      min_val = piArray[i];
      min_idx   = i;
    }
  }
}
Void FindArrayMaxMin(Int* piArray, Int iNumElementInArray, int& max_value, int& min_value)
{
  min_value = 65535;
  max_value = -65535;
  for(int i=0; i< iNumElementInArray; i++)
  {
    if(piArray[i]<min_value)
      min_value = piArray[i];
    if(piArray[i]>max_value)
      max_value = piArray[i];
  }
}

//create, reset and remove resources for a histogram
Void createSEHistogram(Int*& piHist, Int iNumValue)
{
  piHist  = (Int*)calloc(iNumValue, sizeof(Int));
  if(piHist == NULL)
  {
    printf("Error: Cannot allocate memory for a histogram");
    exit(-1);
  }
}
Int* createSEHistogram(Int iNumValue)
{
  return (Int*)calloc(iNumValue, sizeof(Int));
}
Void deleteSEHistogram(Int*& piHist)
{
  if(piHist)
  {
    free(piHist);
    piHist = NULL;
  }
}
Void resetSEHistogram(Int* piHist, Int iNumValue, Int iResetValue)
{
  Int i;
  for(i=0;i<iNumValue;i++)
    piHist[i] = iResetValue;
}

//add one element to the histogram counter if satisfying certain condition
Void addToSEHistogram(Int* piHist, Int iEntrance, Bool bCondition, Int iAddValue)
{
  if(bCondition)
    piHist[iEntrance] += iAddValue;
}

//show the histogram of the array in a slim mode, where only presenting the percentage in a line
Void showSEHistogramSlim(Int* piHist, Int iNumValue, Bool bOutputToFile, FILE* fpLogFile, Char* pcLabelName, Int iLabelNum, int iEntranceOffset, Bool bShowZeroEntry)
{
  int i = 0;
  int total_num = 0;
  double percet = 0;

  for(i=0;i<iNumValue;i++)
    total_num+=piHist[i];

  //print info
  if(bOutputToFile)
  {
    fprintf(fpLogFile,"\n%8s %2d :\t", pcLabelName, iLabelNum);
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(fpLogFile,"%3d\t%5.1f%%\t", i+iEntranceOffset, piHist[i]/(double)total_num*100.0);
    }
  }
  else
  {
    fprintf(stderr,"\n%8s %2d :\t", pcLabelName, iLabelNum);
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(stderr,"%3d\t%5.1f%%\t", i+iEntranceOffset, piHist[i]/(double)total_num*100.0);
    }
  }
}

//show the histogram of the array in a more detailed version, where presenting the entry index, the number, the percentage in multiple lines
Void showSEHistogram(Int* piHist, Int iNumValue, Bool bOutputToFile, FILE* fpLogFile, Char* pcLabelName, Int iLabelNum, int iEntranceOffset, Bool bShowZeroEntry)
{
  int i = 0;
  int total_num = 0;
  double percet = 0;

  for(i=0;i<iNumValue;i++)
    total_num+=piHist[i];

  //print info
  if(bOutputToFile)
  {
    fprintf(fpLogFile,"\n %8s %2d :", pcLabelName, iLabelNum);
    fprintf(fpLogFile,"\n  Value       : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(fpLogFile, "%6d\t", i+iEntranceOffset);
    }

    fprintf(fpLogFile,"\n  Numbers     : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(fpLogFile,"%6d\t", piHist[i]);
    }

    fprintf(fpLogFile,"\n  Percentage  : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(fpLogFile,"%5.1f%%\t", piHist[i]/(double)total_num*100.0);
    }
  }
  else
  {
    fprintf(stderr,"\n %8s %2d :", pcLabelName, iLabelNum);
    fprintf(stderr,"\n  Value       : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(stderr, "%6d\t", i+iEntranceOffset);
    }

    fprintf(fpLogFile,"\n  Numbers     : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(stderr,"%6d\t", piHist[i]);
    }

    fprintf(stderr,"\n  Percentage  : ");
    for(i=0;i<iNumValue;i++)
    {
      if( bShowZeroEntry || (piHist[i]!=0 && !bShowZeroEntry) )
        fprintf(stderr,"%5.1f%%\t", piHist[i]/(double)total_num*100.0);
    }
  }
}

/*--------------------------------------------- log file functions  ------------------------------------*/
Void createLogFile(FILE*& fpLogFile, Char* pcLogFileName)
{
  fpLogFile = fopen(pcLogFileName,"a");
  if(fpLogFile == NULL)
  {
    printf("\n Error Cannot open or create Debug Log File [%s]", pcLogFileName);
    exit(0);
  }
}
FILE* createLogFile(Char* pcLogFileName)
{
  FILE* fpLogFile = fopen(pcLogFileName,"a");
  if(fpLogFile == NULL)
  {
    printf("\n Error Cannot open or create Debug Log File [%s]", pcLogFileName);
    exit(0);
  }
  return fpLogFile;
}
Void closeLogFile(FILE*& fpLogFile)
{
  if(fpLogFile)
    fclose(fpLogFile);
  fpLogFile =  NULL;
}

/*--------------------------------------------- Map file functions  ------------------------------------*/
//create a map file and allocate picture memory for the file
FILE* createMapFile(int width, int height, UChar*& pucMap, Char* pcFileName)
{
  FILE* fpFile = fopen(pcFileName,"ab");
  if(pcFileName == NULL)
  {
    printf("\n Error: Cannot open or create Map File [%s]", pcFileName);
    exit(0);
  }
  pucMap = (UChar*)malloc(width*height*sizeof(UChar));
  if(pucMap == NULL)
  {
    printf("\n Error: Cannot allocate memory for Map file %s", pcFileName);
  }
  memset(pucMap, 0, width*height);

  return fpFile;
}
//close a map file and free picture memory associated with the file
Void  closeMapFile(FILE*& fpMapFile, UChar*& pucMap)
{
  if(fpMapFile)
    fclose(fpMapFile);
  fpMapFile = NULL;

  if(pucMap)
    free(pucMap);
  pucMap = NULL;
}

//function: draw a straight line from (x_start,y_start) to (x_end, y_end)
Void DrawLine(UChar* pucMap, int x_start, int y_start, int x_end, int y_end, int width, int height, UChar LineValue, int mode)
{
  int j,pos,jstep, x, y;
  double k;
  UChar BlendValue;

  if(x_start>=width)
    x_start = width - 1;
  if(x_end>=width)
    x_end   = width - 1;
  if(y_start>=height)
    y_start = height - 1;
  if(y_end>=height)
    y_end   = height - 1;

  if(mode == 0) // draw a line to replace the original pixel (i.e., a solid line)
  {
    if(x_end == x_start && y_end == y_start)
    {
      pucMap[x_end+y_end*width] = LineValue;
      return;
    }

    if( abs(y_end-y_start) >= abs(x_end-x_start))
    {
      k = (x_end-x_start)/(double)(y_end-y_start);
      if(y_end>=y_start)
        jstep = 1;
      else
        jstep = -1;
      for(j=y_start;j<=y_end;j+=jstep)
      {
        x = (int)floor(k*(j-y_start) + x_start + 0.5);
        y = j;
        pos = x + y*width;
        pucMap[pos] = LineValue;
      }
    }
    else
    {
      k = (y_end-y_start)/(double)(x_end-x_start);
      if(x_end>=x_start)
        jstep = 1;
      else
        jstep = -1;
      for(j=x_start;j<=x_end;j+=jstep)
      {
        x = j;
        y = (int)floor(k*(j-x_start) + y_start + 0.5);
        pos = x + y*width;
        pucMap[pos] = LineValue;
      }
    }
    return;
  }

  //if mode == 1, the line blended with the original pixel (i.e., a semi-transparent line)
  if(x_end == x_start && y_end == y_start)
  {
    BlendValue = pucMap[x_end+y_end*width];
    pucMap[x_end+y_end*width] = (LineValue+BlendValue)>>1;
    return;
  }

  if( abs(y_end-y_start) >= abs(x_end-x_start))
  {
    k = (x_end-x_start)/(double)(y_end-y_start);
    if(y_end>=y_start)
      jstep = 1;
    else
      jstep = -1;
    for(j=y_start;j<=y_end;j+=jstep)
    {
      x = (int)floor(k*(j-y_start) + x_start + 0.5);
      y = j;
      pos = x + y*width;
      BlendValue = pucMap[pos];
      pucMap[pos] = (LineValue + BlendValue)>>1;
    }
  }
  else
  {
    k = (y_end-y_start)/(double)(x_end-x_start);
    if(x_end>=x_start)
      jstep = 1;
    else
      jstep = -1;
    for(j=x_start;j<=x_end;j+=jstep)
    {
      x = j;
      y = (int)floor(k*(j-x_start) + y_start + 0.5);
      pos = x + y*width;
      BlendValue = pucMap[pos];
      pucMap[pos] = (LineValue + BlendValue)>>1;
    }
  }
}

int remove_masked_elements(int* piArrayOne, int* piMask, int iNumElement, int* piArrayOneNew)
{
  int idx = 0;
  for(int i=0; i<iNumElement; i++)
  {
    if(piMask[i])
      piArrayOneNew[idx++] = piArrayOne[i];
  }
  return idx;
}

//using CImgLib
#if DECODER_INFO_COLLECT
//calc the correlation of two arrays
double get_correlation(int* piArrayOne, int* piArrayTwo, int iNumElement)
{
  CImg<double> A(piArrayOne, 1, iNumElement);
  CImg<double> B(piArrayTwo, 1, iNumElement);
  double corr2 = 0;

  A -=A.mean();
  B -=B.mean();

  double val1 = A.get_mul(B).sum();
  double val2 = A.pow(2).sum();
  double val3 = B.pow(2).sum();

  if ( abs(val2 * val3) < 1e-6 )
    corr2 = 0;
  else
    corr2 = val1/sqrt(val2*val3);

  return corr2;
}
#endif

double entropy(int* piArray, int iNumElementInArray)
{
  Int max_value = 0;
  Int min_value = 0; 
  double E = 0; //entropy value
  double val = 0;
  int iNumValueInArray = 0;
  int* piHist = NULL;

  FindArrayMaxMin(piArray, iNumElementInArray, max_value, min_value);
  iNumValueInArray = max_value - min_value + 1;
  createSEHistogram(piHist,iNumValueInArray);
  resetSEHistogram(piHist, iNumValueInArray);

  for(int i=0; i<iNumElementInArray; i++)
    addToSEHistogram(piHist, piArray[i]-min_value, true, 1);

  for(int i=0;i<iNumValueInArray;i++)
  {
    if(piHist[i] == 0)
      val = 0;
    else
    {
      val = piHist[i]/double(iNumElementInArray);
      val = -val * log(val);
    }
    E += val;
  }

  deleteSEHistogram(piHist);
  return E;
}

void calcImageEntropy(double* Img_entropy, Pel* pOrg, int iWidth, int iHeight, int iStride, bool bShowEntropyImg)
{
  int es_strength = 10;
  int blk_w = 8, blk_h = 8;
  int bg_w  = 12, bg_h = 12;
  int x_start  = (bg_w-blk_w)/2;
  int x_end    = iWidth -  (bg_w-blk_w)/2 - blk_w;
  int y_start  = (bg_h-blk_h)/2;
  int y_end    = iHeight - (bg_h-blk_h)/2 - blk_h;
  int x_offset = (bg_w-blk_w)/2;
  int y_offset = (bg_h-blk_h)/2;
  int i,j,x,y;
  Pel* pBlk     = (Pel*)malloc(bg_w*bg_h*sizeof(Pel));
  Pel* pPattern = (Pel*)malloc(5*5*sizeof(Pel));
  int a = 0, b = 0;
  int max_val = 0, min_val = 0, max_idx = 0, min_idx = 0;

  int edge_th[3] = {10,30,60}; //three-level edge strength quantization
  int max_level  = 3;
  int v[4];//gradient along 4 directions
  double max_num     = 4 * max_level + 1;
  double max_entropy = log(max_num);
  double E = 0;
  int sampling_step  = 2;
  int num_primitives = (blk_w/sampling_step)*(blk_h/sampling_step);
  Int* primitive_idx = (Int*)malloc(num_primitives*sizeof(Int));

  for(y=y_start; y<=y_end; y+=blk_h)
  {
    for(x=x_start; x<=x_end; x+=blk_w)
    {
      int pos_idx = 0;
      for(int n=0;n<num_primitives;n++)
        primitive_idx[n] = 0; //reset

      fetch_block(pOrg, pBlk, iWidth, iHeight, x-x_offset, y-y_offset, bg_w, bg_h, iStride);
      for(j=y_offset; j<y_offset+blk_h; j+=sampling_step)
      {
        for(i=x_offset; i<x_offset+blk_w; i+=sampling_step)
        {
          fetch_block(pBlk,pPattern, bg_w, bg_h, i, j, 5, 5, bg_w);
          a = abs(pPattern[2*5+0]-pPattern[2*5+4]);  b = abs(pPattern[2*5+1]-pPattern[2*5+3]); //horizontal
          v[0] = max(a,b);
          a = abs(pPattern[0*5+2]-pPattern[4*5+2]);  b = abs(pPattern[1*5+2]-pPattern[3*5+2]); //vertical
          v[1] = max(a,b);
          a = abs(pPattern[0*5+4]-pPattern[4*5+0]);  b = abs(pPattern[1*5+3]-pPattern[3*5+1]); //down-left
          v[2] = max(a,b);
          a = abs(pPattern[0*5+0]-pPattern[4*5+4]);  b = abs(pPattern[1*5+1]-pPattern[3*5+3]); //down-right
          v[3] = max(a,b);
          FindArrayMax(v,4,max_val,max_idx);
          FindArrayMin(v,4,min_val,min_idx);
          for(int idx = max_level-1; idx>=0; idx--)
          {
            if(max_val > edge_th[idx])
            {
              primitive_idx[pos_idx++] = 1 + min_idx + 4 * idx;
              break;
            }
          }
        }
      }
      //calculate entropy from primitive index array
      //printf("\n"); for(int i=0;i<16;i++) printf("%d ", primitive_idx[i]);
      E = entropy(primitive_idx, num_primitives);
      E = (1.0-E/max_entropy)*0.75+0.25;
      assign_block(Img_entropy, E, iWidth, iHeight, x, y, blk_w, blk_h);
    }
  }

  if (bShowEntropyImg)
  {
    //show image
    UChar* pY = NULL;
    FILE* fpMap = createMapFile(iWidth, iHeight, pY, "entropy_map.yuv");
    UChar* pUV = (UChar*)malloc(iWidth*iHeight/2);
    memset(pUV, 128, iWidth*iHeight/2);

    for(i=0;i<iWidth*iHeight;i++)
      pY[i] = UChar(Img_entropy[i]*255);
    fwrite(pY, sizeof(UChar), iWidth*iHeight,   fpMap);
    //UV of 128
    fwrite(pUV, sizeof(UChar), iWidth*iHeight/2, fpMap);

    for(j=0;j<iHeight;j++)
    {
      for(i=0;i<iWidth;i++)
      {
        pY[j*iWidth + i] = UChar(pOrg[j*iStride + i]);
      }
    }

    fwrite(pY,  sizeof(UChar), iWidth*iHeight,   fpMap);
    fwrite(pUV, sizeof(UChar), iWidth*iHeight/2, fpMap);

    fwrite(pY,  sizeof(UChar), iWidth*iHeight,   fpMap);
    int idx = 0;
    int idx2 = 0;
    memset(pY, 128, iWidth*iHeight);
    for(j=0;j<iHeight;j+=2)
    {
      for(i=0;i<iWidth;i+=2)
      {
        idx2 = j*iWidth+i;
        pY[idx++]  = (UChar)(Img_entropy[idx2]*128+64);
      }
    }

    fwrite(pY, sizeof(UChar), iWidth*iHeight/4, fpMap);
    fwrite(pY, sizeof(UChar), iWidth*iHeight/4, fpMap);

    closeMapFile(fpMap, pY);
    free(pUV);
  }

  if(pBlk)
    free(pBlk);
  if(primitive_idx)
    free(primitive_idx);
}

//get log2 value for frequently used integer value
int getLog2Value(int inputValue)
{
  int Log2Value = 0;
  switch (inputValue)
  {
    case 2:
    {
      Log2Value = 1;      break;
    }
    case 4:
    {
      Log2Value = 2;      break;
    }
    case 8:
    {
      Log2Value = 3;      break;
    }
    case 16:
    {
      Log2Value = 4;      break;
    }
    case 32:
    {
      Log2Value = 5;      break;
    }
    case 64:
    {
      Log2Value = 6;      break;
    }
    case 128:
    {
      Log2Value = 7;      break;
    }
    default:
    {
      assert(0);
    }
  }
  return Log2Value;
}
