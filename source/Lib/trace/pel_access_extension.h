#ifndef __PEL_ACCESS_EXTENSION_H_
#define __PEL_ACCESS_EXTENSION_H_

#include "../TLibCommon/CommonDef.h"
#include "../TLibCommon/TComYuv.h"
#include "../TLibCommon/TComInterpolationFilter.h"
#include "../TLibCommon/TComPicYuv.h"

#if JVET_D0033_ADAPTIVE_CLIPPING
template <typename T> inline T ClipBD(const T x, const Int bitDepth) { return Clip3(T(0), T((1 << bitDepth) - 1), x); }
#endif

class TComYuvExt
{
  TComYuv *_comYuv;
public:
  TComYuvExt(TComYuv *comYuv) : _comYuv(comYuv) {}
  virtual ~TComYuvExt(){}
  Void    copyToPicPartYuv(TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, float scale = 1, int norm = 0);
  Void    copyToPicPartComponent(const ComponentID id, TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, float scale = 1, int norm = 0);

  Void    subToPicPartYuv(TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx);
  Void    subToPicPartComponent(const ComponentID id, TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx);

  Void    sumToPicPartYuv(TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, Int shift_right);
  Void    sumToPicPartComponent(const ComponentID id, TComPicYuv*    pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, Int shift_right);
};

class TComPicYuvExt : public TComPicYuv
{
  //TComPicYuv* _comPicYuv;
public:
  //TComPicYuvExt(TComPicYuv *comPicYuv);
  TComPicYuvExt();
  virtual ~TComPicYuvExt(){}
  Void  copyToPicPart(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);
  Void  copyToPicPartComponent(ComponentID compId, TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);

  Void  subPicPart(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);
  Void  subPicPartLuma(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);
  Void  subPicPartCb(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);
  Void  subPicPartCr(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress);

  Void clean();
  void cleanLuma();
  void set(Pel init_val);
  //*
  Pel*  getLumaAddr(Int iCuAddr, Int X, Int Y);
  Pel*  getLumaAddr(Int iCuAddr, Int uiAbsZorderIdx, Int X, Int Y);
  Pel*  getLumaAddrXY(Int X, Int Y);
  Pel*  getCrAddrXY(Int X, Int Y);
  Pel*  getAddrXY(Int X, Int Y, ComponentID compId);
};

#endif//__PEL_ACCESS_EXTENSION_H_
