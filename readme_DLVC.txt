Requirements:
Windows
Visual Studio 2015 64bit
GPU supporting CUDA8.0, don't forget to install the GPU driver.

Steps:
1. Download caffe, tensorflow and CUDA from:
https://pan.baidu.com/s/1jlJ1ux-ZVQhzj6Msxn4KgA  password: rq0l
2. Extract the above caffe and CUDA packages to this folder. As for the tensorflow packages, you should extract the corresponding one according to the running mode CPU or GPU. 
3. Use VS to open the corresponding solution, and build (Debug or Release, but only x64!)
4. Please note, after building, there will be some DLL files in the Release/Debug folder. If you want to copy the executables to another machine, do not forget these DLL files
5. Download caffe and tensorflow models, from:
https://pan.baidu.com/s/1RWMfA7OrKzKZw2upJ3E14w password: 8n4s
6. Extract the above two packages to working folder
7. When running executables, you can use the commands "-dl Caffe/Tensorflow/Opencv" "-mode CPU/GPU" and "-gpuid 0/1/2..." to select different deep learning frameworks, different running modes and specific GPU in your machine.


