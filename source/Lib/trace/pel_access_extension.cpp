#include "pel_access_extension.h"

Void TComYuvExt::copyToPicPartYuv(TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, float scale, int norm)
{
  for (Int comp = 0; comp < _comYuv->getNumberValidComponents(); comp++)
  {
    copyToPicPartComponent(ComponentID(comp), pcYuvDst, iCuAddr, uiAbsZorderIdx, scale, norm);
  }
}

Void TComYuvExt::copyToPicPartComponent(const ComponentID compID, TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, float scale, int norm)
{
  Int  y;

  Pel* pSrc = _comYuv->getAddr(compID);
  Pel* pDst = pcYuvDst->getAddr(compID, iCuAddr, uiAbsZorderIdx);
  Int offset = 0;
  Int shiftNum = 0;
  if (norm)
  {
    shiftNum = IF_INTERNAL_PREC /*+ 1*/ - 8;
    offset = (1 << (shiftNum - 1)) +/* 2 * */IF_INTERNAL_OFFS;
  }

  UInt  iSrcStride = _comYuv->getStride(compID);
  UInt  iDstStride = pcYuvDst->getStride(compID);
  const Int iWidth = _comYuv->getWidth(compID);
  const Int iHeight = _comYuv->getHeight(compID);
  for (y = iHeight; y != 0; y--)
  {
    ::memcpy(pDst, pSrc, sizeof(Pel)*iWidth);
    if (scale != 1)
    {
      for (int i = 0; i < iWidth; i++)
        pDst[i] = (Pel)(pDst[i] * scale);
    }
    if (norm)
    {
      for (int i = 0; i < iWidth; i++)
      {
#if JVET_D0033_ADAPTIVE_CLIPPING
        pDst[i] = ClipA((pDst[i] + offset) >> shiftNum, compID);
#else
        pDst[i] = ClipBD((pDst[i] + offset) >> shiftNum, 8);
#endif
      }
    }
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}

Void TComYuvExt::subToPicPartYuv(TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx)
{
  for (Int comp = 0; comp < _comYuv->getNumberValidComponents(); comp++)
  {
    subToPicPartComponent(ComponentID(comp), pcYuvDst, iCuAddr, uiAbsZorderIdx);
  }
}

Void TComYuvExt::subToPicPartComponent(const ComponentID compID, TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx)
{
  Int  y, x;

  Pel* pSrc = _comYuv->getAddr(compID);
  Pel* pDst = pcYuvDst->getAddr(compID, iCuAddr, uiAbsZorderIdx);

  UInt  iSrcStride = _comYuv->getStride(compID);
  UInt  iDstStride = pcYuvDst->getStride(compID);
  const Int iWidth = _comYuv->getWidth(compID);
  const Int iHeight = _comYuv->getHeight(compID);
  for (y = iHeight; y != 0; y--)
  {
    for (x = 0; x < iWidth; x++)
    {
      pDst[x] -= pSrc[x];
    }
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}

Void TComYuvExt::sumToPicPartYuv(TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, Int shift_right)
{
  for (Int comp = 0; comp < _comYuv->getNumberValidComponents(); comp++)
  {
    sumToPicPartComponent(ComponentID(comp), pcYuvDst, iCuAddr, uiAbsZorderIdx, shift_right);
  }
}

Void TComYuvExt::sumToPicPartComponent(const ComponentID compID, TComPicYuv* pcYuvDst, Int iCuAddr, Int uiAbsZorderIdx, Int shift_right)
{
  Int  y, x;

  Pel* pSrc = _comYuv->getAddr(compID);
  Pel* pDst = pcYuvDst->getAddr(compID, iCuAddr, uiAbsZorderIdx);

  UInt  iSrcStride = _comYuv->getStride(compID);
  UInt  iDstStride = pcYuvDst->getStride(compID);
  Int offset = 0;
  Int shiftNum = 0;
  if (shift_right)
  {
    shiftNum = IF_INTERNAL_PREC + 1 - 8;
    offset = (1 << (shiftNum - 1)) + 2 * IF_INTERNAL_OFFS;
  }
  const Int iWidth = _comYuv->getWidth(compID);
  const Int iHeight = _comYuv->getHeight(compID);
  for (y = iHeight; y != 0; y--)
  {
    for (x = 0; x < iWidth; x++)
    {
      pDst[x] += pSrc[x];
#if JVET_D0033_ADAPTIVE_CLIPPING
      pDst[x] = ClipA((pDst[x] + offset) >> shiftNum, compID);
#else
      pDst[x] = ClipBD((pDst[x] + offset) >> shiftNum, 8);
#endif
    }
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}

//////////////////////////////////////////////////////////////////////////
TComPicYuvExt::TComPicYuvExt()
{
}

Void  TComPicYuvExt::copyToPicPart(TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress)
{
  copyToPicPartComponent(COMPONENT_Y, pcPicYuvDst, iCuAddr, uiAbsZorderIdx, uiWidth, uiHeight, bSrcUseCUAdress);
  copyToPicPartComponent(COMPONENT_Cb, pcPicYuvDst, iCuAddr, uiAbsZorderIdx, uiWidth, uiHeight, bSrcUseCUAdress);
  copyToPicPartComponent(COMPONENT_Cr, pcPicYuvDst, iCuAddr, uiAbsZorderIdx, uiWidth, uiHeight, bSrcUseCUAdress);
}

Void  TComPicYuvExt::copyToPicPartComponent(ComponentID compId, TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress)
{
  assert(getWidth(compId) == pcPicYuvDst->getWidth(compId));
  assert(getHeight(compId) == pcPicYuvDst->getHeight(compId));

  Pel* pSrc = getAddr(compId);
  if (bSrcUseCUAdress)
    pSrc = getAddr(compId, iCuAddr, uiAbsZorderIdx);
  Pel* pDst = pcPicYuvDst->getAddr(compId, iCuAddr, uiAbsZorderIdx);

  UInt  iSrcStride = getStride(compId);
  UInt  iDstStride = pcPicYuvDst->getStride(compId);

  uiWidth >>= getComponentScaleX(compId);
  uiHeight >>= getComponentScaleY(compId);

  for (int y = uiHeight; y != 0; y--)
  {
    ::memcpy(pDst, pSrc, sizeof(Pel)*uiWidth);
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}

void TComPicYuvExt::clean()
{
  for (Int chan = 0; chan < getNumberValidComponents(); chan++)
  {
    const ComponentID ch = ComponentID(chan);
    ::memset(getBuf(ch), 0, sizeof(Pel) * getStride(ch) * getTotalHeight(ch));
  }
}

Pel* TComPicYuvExt::getLumaAddr(Int iCuAddr, Int X, Int Y)
{
  return getAddr(COMPONENT_Y, iCuAddr) + X + Y*getStride(COMPONENT_Y);
}

Pel* TComPicYuvExt::getLumaAddr(Int iCuAddr, Int uiAbsZorderIdx, Int X, Int Y)
{
  return getAddr(COMPONENT_Y, iCuAddr, uiAbsZorderIdx) + X + Y*getStride(COMPONENT_Y);
}

Pel* TComPicYuvExt::getLumaAddrXY(Int X, Int Y)
{
  return getAddr(COMPONENT_Y) + X + Y*getStride(COMPONENT_Y);
}

Pel* TComPicYuvExt::getAddrXY(Int X, Int Y, ComponentID compId)
{
  return getAddr(compId) + (X >> getComponentScaleX(compId)) + (Y >> getComponentScaleY(compId))*getStride(compId);
}

/*
// wrapping on already created version
TComPicYuvExt::TComPicYuvExt(TComPicYuv *comPicYuv) : _comPicYuv(comPicYuv)
{
  TComPicYuv* tmp = this;
  *tmp = *comPicYuv;
}

Void  TComPicYuvExt::copyToPicPartComponent(ComponentID compId, TComPicYuv* pcPicYuvDst, UInt iCuAddr, UInt uiAbsZorderIdx, UInt uiWidth, UInt uiHeight, Bool bSrcUseCUAdress)
{
  assert(_comPicYuv->getWidth(compId)  == pcPicYuvDst->getWidth(compId));
  assert(_comPicYuv->getHeight(compId) == pcPicYuvDst->getHeight(compId));

  Pel* pSrc = _comPicYuv->getAddr(compId);
  if (bSrcUseCUAdress)
    pSrc = _comPicYuv->getAddr(compId, iCuAddr, uiAbsZorderIdx);
  Pel* pDst = pcPicYuvDst->getAddr(compId, iCuAddr, uiAbsZorderIdx);

  UInt  iSrcStride = _comPicYuv->getStride(compId);
  UInt  iDstStride = pcPicYuvDst->getStride(compId);

  uiWidth  >>= _comPicYuv->getComponentScaleX(compId);
  uiHeight >>= _comPicYuv->getComponentScaleY(compId);

  for (int y = uiHeight; y != 0; y--)
  {
    ::memcpy(pDst, pSrc, sizeof(Pel)*uiWidth);
    pDst += iDstStride;
    pSrc += iSrcStride;
  }
}

void TComPicYuvExt::clean()
{
  for (Int chan = 0; chan < _comPicYuv->getNumberValidComponents(); chan++)
  {
    const ComponentID ch = ComponentID(chan);
    ::memset(_comPicYuv->getBuf(ch), 0, sizeof(Pel) * _comPicYuv->getStride(ch) * _comPicYuv->getTotalHeight(ch));
  }
}

Pel* TComPicYuvExt::getLumaAddr(Int iCuAddr, Int X, Int Y) 
{ 
  return _comPicYuv->getAddr(COMPONENT_Y, iCuAddr) + X + Y*_comPicYuv->getStride(COMPONENT_Y);
}

Pel* TComPicYuvExt::getLumaAddr(Int iCuAddr, Int uiAbsZorderIdx, Int X, Int Y)
{ 
  return _comPicYuv->getAddr(COMPONENT_Y, iCuAddr, uiAbsZorderIdx) + X + Y*_comPicYuv->getStride(COMPONENT_Y);
}

Pel* TComPicYuvExt::getLumaAddrXY(Int X, Int Y)
{
  return _comPicYuv->getAddr(COMPONENT_Y) + X + Y*_comPicYuv->getStride(COMPONENT_Y);
}
*/