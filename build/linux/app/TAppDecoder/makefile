# the SOURCE definiton lets you move your makefile to another position
CONFIG 				= CONSOLE

# set directories to your wanted values
SRC_DIR				= ../../../../source/App/TAppDecoder
INC_DIR				= ../../../../source/Lib
LIB_DIR				= ../../../../lib
BIN_DIR				= ../../../../bin

SRC_DIR1		=
SRC_DIR2		=
SRC_DIR3		=
SRC_DIR4		=

USER_INC_DIRS	= -I$(SRC_DIR) 
USER_LIB_DIRS	=

ifeq ($(HIGHBITDEPTH), 1)
HBD=HighBitDepth
else
HBD=
endif

# intermediate directory for object files
OBJ_DIR				= ./objects$(HBD)

# set executable name
PRJ_NAME			= TAppDecoder$(HBD)

# defines to set
DEFS				= -DMSYS_LINUX -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DMSYS_UNIX_LARGEFILE

# set objects
OBJS          		= 	\
					$(OBJ_DIR)/decmain.o \
					$(OBJ_DIR)/TAppDecCfg.o \
					$(OBJ_DIR)/TAppDecTop.o \

# set libs to link with
LIBS				= -ldl

DEBUG_LIBS			=
RELEASE_LIBS		=

CAFFE_LIBS          = -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64:/opt/caffe/build/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/local/lib:::::::: \
	/opt/caffe/build/lib/libcaffe.so.1.0.0 /opt/caffe/build/lib/libproto.a -lboost_system -lboost_thread -lboost_filesystem -lboost_chrono -lboost_date_time -lboost_atomic -lpthread -lglog \
	-lgflags /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.so /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.so -lsz -lz -ldl -lm -lpthread -lglog -lgflags \
	/usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.so /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.so -lsz -lz -ldl -lm -lprotobuf -lpthread -llmdb -lleveldb \
	/usr/local/cuda/lib64/libcudart.so /usr/local/cuda/lib64/libcurand.so /usr/local/cuda/lib64/libcublas.so -lcudnn /usr/local/lib/libnccl.so \
	/usr/local/lib/libopencv_highgui.so.3.4.5 /usr/local/lib/libopencv_videoio.so.3.4.5 /usr/local/lib/libopencv_imgcodecs.so.3.4.5 /usr/local/lib/libopencv_imgproc.so.3.4.5 \
	/usr/local/lib/libopencv_core.so.3.4.5 -llapack -lcblas -latlas -lboost_python /usr/local/lib/libopencv_dnn.so.3.4.5 \
	-Wl,-rpath,/opt/caffe/build/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial/lib:/usr/local/cuda/lib64:/usr/local/lib:::::::: -lglog -lopenblas

STAT_LIBS			= -lpthread
DYN_LIBS			= 
TORCH_LIBS          = -L/opt/pytorch/lib/ -ltorch -lc10 -lc10_cuda -lthnvrtc -lfoxi -lcaffe2_detectron_ops_gpu

DYN_DEBUG_LIBS		= -lTLibDecoder$(HBD)d -lTLibCommon$(HBD)d -lTLibVideoIO$(HBD)d -lTAppCommon$(HBD)d $(CAFFE_LIBS) $(TORCH_LIBS)
DYN_DEBUG_PREREQS		= $(LIB_DIR)/libTLibDecoder$(HBD)d.a $(LIB_DIR)/libTLibCommon$(HBD)d.a $(LIB_DIR)/libTLibVideoIO$(HBD)d.a $(LIB_DIR)/libTAppCommon$(HBD)d.a
STAT_DEBUG_LIBS		= -lTLibDecoder$(HBD)Staticd -lTLibCommon$(HBD)Staticd -lTLibVideoIO$(HBD)Staticd -lTAppCommon$(HBD)Staticd $(CAFFE_LIBS) $(TORCH_LIBS)
STAT_DEBUG_PREREQS		= $(LIB_DIR)/libTLibDecoder$(HBD)Staticd.a $(LIB_DIR)/libTLibCommon$(HBD)Staticd.a $(LIB_DIR)/libTLibVideoIO$(HBD)Staticd.a $(LIB_DIR)/libTAppCommon$(HBD)Staticd.a

DYN_RELEASE_LIBS	= -lTLibDecoder$(HBD) -lTLibCommon$(HBD) -lTLibVideoIO$(HBD) -lTAppCommon$(HBD) $(CAFFE_LIBS) $(TORCH_LIBS)
DYN_RELEASE_PREREQS	= $(LIB_DIR)/libTLibDecoder$(HBD).a $(LIB_DIR)/libTLibCommon$(HBD).a $(LIB_DIR)/libTLibVideoIO$(HBD).a $(LIB_DIR)/libTAppCommon$(HBD).a
STAT_RELEASE_LIBS	= -lTLibDecoder$(HBD)Static -lTLibCommon$(HBD)Static -lTLibVideoIO$(HBD)Static -lTAppCommon$(HBD)Static $(CAFFE_LIBS) $(TORCH_LIBS)
STAT_RELEASE_PREREQS	= $(LIB_DIR)/libTLibDecoder$(HBD)Static.a $(LIB_DIR)/libTLibCommon$(HBD)Static.a $(LIB_DIR)/libTLibVideoIO$(HBD)Static.a $(LIB_DIR)/libTAppCommon$(HBD)Static.a


# name of the base makefile
MAKE_FILE_NAME		= ../../common/makefile.base

# include the base makefile
include $(MAKE_FILE_NAME)
