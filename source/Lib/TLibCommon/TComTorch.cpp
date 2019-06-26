#include "TComTorch.h"
#if USE_LIBTORCH
cv::Mat TestFrameInterp(string modelpath, cv::Mat input1, cv::Mat input2, int lqp, int rqp,  bool GPUMode)
{
	torch::NoGradGuard no_grad;

	std::shared_ptr<jit::script::Module> FrameInterpModule;
	if (GPUMode) FrameInterpModule = jit::load(modelpath, torch::kCUDA);
	else FrameInterpModule = jit::load(modelpath);

	float lqpf = float(lqp) / 51.0;
	float rqpf = float(rqp) / 51.0;

	Tensor tensorinput1 = torch::from_blob(input1.data, { 1,input1.rows,input1.cols,3 }, kF32);
	Tensor tensorinput2 = torch::from_blob(input2.data, { 1,input1.rows,input1.cols,3 }, kF32);
	tensorinput1 = tensorinput1.permute({ 0,3,1,2 });
	tensorinput2 = tensorinput2.permute({ 0,3,1,2 });

	int inthei = tensorinput1.size(2);
	int	intwid = tensorinput1.size(3);
	int	intpaddingwidth = 0;
	int	intpaddingheight = 0;
	if (inthei != ((inthei >> 5) << 5))
	{
		intpaddingheight = (((inthei >> 5) + 1) << 5) - inthei;
	}
	if (intwid != ((intwid >> 5) << 5))
	{
		intpaddingwidth = (((intwid >> 5) + 1) << 5) - intwid;
	}
	int	intpaddingleft = int(intpaddingwidth / 2);
	int intpaddingright = intpaddingwidth - int(intpaddingwidth / 2);
	int intpaddingtop = int(intpaddingheight / 2);
	int intpaddingbottom = intpaddingheight - int(intpaddingheight / 2);
	tensorinput1 = torch::replication_pad2d(tensorinput1, { intpaddingleft, intpaddingright, intpaddingtop, intpaddingbottom });
	tensorinput2 = torch::replication_pad2d(tensorinput2, { intpaddingleft, intpaddingright, intpaddingtop, intpaddingbottom });

	Tensor tensorinputlqp = torch::zeros({ 1,1,tensorinput1.size(2),tensorinput1.size(3) }) + lqpf;
	Tensor tensorinputrqp = torch::zeros({ 1,1,tensorinput1.size(2),tensorinput1.size(3) }) + rqpf;

	Tensor tresulttensor;
	if (GPUMode)
	{
		tresulttensor = FrameInterpModule->forward({ tensorinput1.to(torch::kCUDA),tensorinput2.to(torch::kCUDA),tensorinputlqp.to(torch::kCUDA) ,tensorinputrqp.to(torch::kCUDA) }).toTensor();
	}
	else
	{
		if (modelpath == "./LDFInterp.pt")
		{
			tresulttensor = FrameInterpModule->forward({ tensorinput1,tensorinput2,tensorinputlqp ,tensorinputrqp }).toTensor();
		}
		else
		{
			tresulttensor = FrameInterpModule->forward({ tensorinput1,tensorinput2 }).toTensor();
		}
	}


	tresulttensor = torch::replication_pad2d(tresulttensor, { 0 - intpaddingleft, 0 - intpaddingright, 0 - intpaddingtop, 0 - intpaddingbottom });

	if (GPUMode) tresulttensor = tresulttensor.to(kCPU);
	auto resulttensor = tresulttensor.detach();
	resulttensor = resulttensor.permute({ 0,2,3,1 });
	resulttensor = resulttensor.squeeze(0);
	resulttensor = resulttensor.clamp(0.0, 1.0);
	cv::Mat resimg(input1.rows, input1.cols, CV_32FC3);
	std::memcpy((void*)resimg.data, resulttensor.data_ptr(), resulttensor.numel() * sizeof(float));

	return resimg;
}

int TestLibTorchUtil()
{
	string left_filename = "left.png";
	string right_filename = "right.png";
	string default_model("./LDFInterp.pt");

	cv::Mat left_img = cv::imread(left_filename);
	cv::Mat right_img = cv::imread(right_filename);
        cv::Mat left_blob, right_blob;
        left_img.convertTo(left_blob, CV_32FC3, 1.0/255);
        right_img.convertTo(right_blob, CV_32FC3, 1.0/255);

	cv::Mat result = TestFrameInterp(default_model, left_blob, right_blob, 32, 32, false);
        cv::Mat wimg;
        result.convertTo(wimg, CV_8UC3, 255);
	printf("%.4f\n", result.at<cv::Vec3f>(5,5)[0]);
        printf("%d\n", wimg.at<cv::Vec3b>(5,5)[0]);
	cv::imwrite("output.bmp", wimg);
	return 0;
}
#endif
