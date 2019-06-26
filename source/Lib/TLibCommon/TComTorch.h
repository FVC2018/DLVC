#ifndef __TCOMTORCH__
#define __TCOMTORCH__

#include <iostream>

#include "CommonDef.h"
#if USE_LIBTORCH
#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <memory>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace torch;

cv::Mat TestFrameInterp(string modelpath, cv::Mat input1, cv::Mat input2, int lqp, int rqp, bool GPUMode = false);
int TestLibTorchUtil();
#endif
#endif