#ifndef __TCOMTF__
#define __TCOMTF__

#include "CommonDef.h"
#include "TComPic.h"
#if USE_TENSORFLOW
#include "TensorFlowHeader.h"
#include <fstream>
#include <utility>
#include <vector>
#include <eigen/Dense>

#include "tensorflow/cc/ops/const_op.h"
#include "tensorflow/cc/ops/image_ops.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/graph.pb.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/graph/default_device.h"
#include "tensorflow/core/graph/graph_def_builder.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/core/threadpool.h"
#include "tensorflow/core/lib/io/path.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/lib/strings/stringprintf.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/init_main.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/util/command_line_flags.h"



// These are all common classes it's handy to reference with no namespace.
using namespace tensorflow;
using tensorflow::Flag;
using tensorflow::Tensor;
using tensorflow::Status;
using tensorflow::string;
using tensorflow::int32;

class TComTF
{

public:
	TComTF();
	virtual ~TComTF();
	Status LoadGraph(const string& graph_file_name, std::unique_ptr<tensorflow::Session>* session);
	Tensor ReadTensorFromYuv(Pel* piSrc, Int inputHeight, Int inputWidth, Int iStride);
	Void TFNetForward(const string graph_path, char* cGPUid, Pel* piSrc, Int inputHeight, Int inputWidth, Int iStride, string inputLayer, string outputLayer);
	Status Tensor2Yuv(Pel* piSrc, Int inputHeight, Int inputWidth, Int iStride, Tensor inputTensors);
	Int setenv(const char *name, const char *value, int overwrite);
};

#endif

#endif

