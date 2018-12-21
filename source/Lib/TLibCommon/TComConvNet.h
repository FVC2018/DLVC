

#ifndef __TCOMCONVNET__
#define   __TCOMCONVNET__

#include "TComCaffe.h"
#include "CommonDef.h"
#if frameworkCombine
#include "TComOpencv.h"
#endif

class TComConvNet
{
private:
  string                    m_netFileName;
  string                    m_caffeModelName;
  TComCaffe                 m_pcCaffe;
  cv::Mat                   m_inputBlock;
  std::vector<Float>        m_networkOutput;
  Int                       m_iInputHeight;
  Int                       m_iInputWidth;
#if frameworkCombine
  TComOpencv                m_pcOpencv;
  Int                       m_iOutputHeight;
  Int                       m_iOutputWidth;
  string                    m_dlPlatform;
#endif

public:
  TComConvNet();
  ~TComConvNet();

  Void                create(const string& netFile, const string& modelFile, Char* query_blob_name,
#if frameworkCombine
    string platform,
#endif
    Char* mode, Int iGPUid, cv::Size inputSize
#if frameworkCombine
    , cv::Size outputSize
#endif
  );
  Void                destroy(const string& netFile, const string& modelFile, Char* mode
#if frameworkCombine
    , string platform
#endif
  );

  Int                 getInputHeight() { return m_iInputHeight; }
  Int                 getInputWidth() { return m_iInputWidth; }

  Void                setInputBlock(Pel* src, Int Stride, cv::Mat& inputMat, Int inputHeight, Int inputWidth);
  Void                setInputBlockSize(Int inputWidth, Int inputHeight);
#if frameworkCombine
  Void                setOutputBlockSize(Int outputWidth, Int outputHeight);
  string			  getPlatform() { return m_dlPlatform; }
#endif
  cv::Mat             getInputBlock() { return m_inputBlock; }
  Void                ConvNetForward(cv::Mat& inputMat);

  std::vector<Float>  getConvNetOutput();
  std::vector<Float>  getConvNetOutputCenter(Int iMargin);

#if CNN_BASED_FRAC_INTERPOLATION
  Void copyFromVectorToPointer(std::vector<float>& src, Pel* pDst, Int iStrideDst, Int iheight, Int iwidth, Bool isLast = true);
#endif
};

#endif