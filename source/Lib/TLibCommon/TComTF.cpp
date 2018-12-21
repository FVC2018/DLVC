#include "TComTF.h"
#if USE_TENSORFLOW
TComTF::TComTF() {}
TComTF::~TComTF() {}

Status TComTF::LoadGraph(const string& graph_file_name,
	std::unique_ptr<tensorflow::Session>* session) {
	tensorflow::GraphDef graph_def;
	Status load_graph_status =
		ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
	if (!load_graph_status.ok()) {
		return tensorflow::errors::NotFound("Failed to load compute graph at '",
			graph_file_name, "'");
	}
	session->reset(tensorflow::NewSession(tensorflow::SessionOptions()));
	Status session_create_status = (*session)->Create(graph_def);
	if (!session_create_status.ok()) {
		return session_create_status;
	}
	return Status::OK();
}


Tensor TComTF::ReadTensorFromYuv(Pel* piSrc, Int inputHeight, Int inputWidth, Int iStride)
{
	std::vector<float> inputYuv;
	for (int i = 0;i < inputHeight;i++)
	{
		for (int j = 0;j < inputWidth;j++)
		{
			inputYuv.push_back(piSrc[j]/ PixelRangeDouble);
		}
		piSrc += iStride;
	}
	auto mapped_X_ = Eigen::TensorMap<Eigen::Tensor<float, 4, Eigen::RowMajor>>
		(&inputYuv[0], 1, inputHeight, inputWidth,1);
	auto eigen_X_ = Eigen::Tensor<float, 4, Eigen::RowMajor>(mapped_X_);
	
	Tensor out_tensors(DT_FLOAT, TensorShape({ 1,inputHeight, inputWidth,1 }));
	out_tensors.tensor<float, 4>() = eigen_X_;
	return out_tensors;
}
Status TComTF::Tensor2Yuv(Pel* piSrc, Int inputHeight, Int inputWidth, Int iStride, Tensor inputTensors)
{
	tensorflow::TTypes<float>::Flat flatTensor = inputTensors.flat<float>();
	for (int i = 0;i < inputHeight;i++)
	{
		for (int j = 0;j < inputWidth;j++)
		{
			piSrc[j]=Clip3(0, PixelRangeInt,int(flatTensor(i*inputWidth+j)*PixelRangeDouble +0.5));
		}
		piSrc += iStride;
	}
	return Status::OK();
}
Int TComTF::setenv(const char *name, const char *value, int overwrite)
{
	int errcode = 0;
	if (!overwrite) {
		size_t envsize = 0;
		errcode = getenv_s(&envsize, NULL, 0, name);
		if (errcode || envsize) return errcode;
	}
	return _putenv_s(name, value);
}

Void TComTF::TFNetForward(const string graph_path, char* cGPUid, Pel* piSrc,Int inputHeight, Int inputWidth,Int iStride, string inputLayer, string outputLayer) {
	
	setenv("CUDA_VISIBLE_DEVICES", cGPUid, 1);

	std::unique_ptr<tensorflow::Session> session;
	
	TF_CHECK_OK(LoadGraph(graph_path, &session));
	
	Tensor inputTensor;
	inputTensor=ReadTensorFromYuv(piSrc, inputHeight, inputWidth, iStride);
	
	std::vector<Tensor> outputs;
	TF_CHECK_OK(session->Run({ { inputLayer, inputTensor } }, { outputLayer }, {}, &outputs));
	
	Tensor Y_ = outputs[0];
	TF_CHECK_OK(Tensor2Yuv(piSrc, inputHeight, inputWidth, iStride, Y_));
	session->Close();
}

#endif
