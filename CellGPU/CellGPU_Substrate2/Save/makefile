#This makefile has crazy-looking paths and includes, reflecting a makefile that works on the different
#operating systems and computers regulating used by the developer. For your system you can definitely,
#definitely simplify this. Currently, this will auto-compile all cpp files in the main directory,
#and make objects for all cpp and cu files in the src directory

CXX := g++
CC := gcc
LINK := g++ -fPIC
#NVCC := nvcc

INCLUDES := -I. -I./inc/ -I/opt/local/include
INCLUDES += #-I/usr/lib/cuda/includes -I/usr/lib/cuda/include
INCLUDES += -I/home/diogo/CellGPULib/CGAL-4.9/include -I/home/diogo/CellGPUPackages/CGAL/include -I/home/linuxbrew/.linuxbrew/Cellar/gmp/6.1.2_2/include -I/home/linuxbrew/.linuxbrew/Cellar/mpfr/4.0.1/include -I/home/linuxbrew/.linuxbrew/Cellar/boost/1.68.0/include
INCLUDES += -I/home/diogo/CellGPUPackages/NetCDF/include -I/home/diogo/CellGPULib/netcdf-cxx4-4.3.0 -I/home/diogo/CellGPULib/netcdf-cxx-4.2 -I/home/diogo/CellGPULib/eigen-eigen-323c052e1731
INCLUDES += -I/home/user/MEGA/cenas/CellGPU_Substrate/boost_1_68_0


LIB_CUDA := #-L. -L/usr/lib/cuda/lib64 -L/usr/lib/cuda/lib -lcuda -lcudart
LIB_CGAL := -L/home/diogo/CellGPUPackages/CGAL/lib -L/home/linuxbrew/.linuxbrew/Cellar/gmp/6.1.2_2/lib -L/home/linuxbrew/.linuxbrew/Cellar/mpfr/4.0.1/lib -L/home/diogo/CellGPULib/CGAL-4.9/lib -lCGAL -lCGAL_Core -lgmp -lmpfr
LIB_NETCDF = -lnetcdf -lnetcdf_c++ -L/opt/local/lib -L/home/diogo/CellGPUPackages/NetCDF/lib

#common flags
COMMONFLAGS += -std=c++11 -DCGAL_DISABLE_ROUNDING_MATH_CHECK -O3 -D_FORCE_INLINES
#NVCCFLAGS += -arch=sm_35 $(COMMONFLAGS) -Wno-deprecated-gpu-targets #-Xptxas -fmad=false#-O0#-dlcm=ca#-G
CXXFLAGS += $(COMMONFLAGS)
CXXFLAGS += -w -frounding-math
CFLAGS += $(COMMONFLAGS) -frounding-math

CUOBJ_DIR=obj/cuobj
MODULES = databases models updaters utility analysis fit
INCLUDES += -I./inc/databases -I./inc/models -I./inc/updaters -I./inc/utility -I./inc/analysis -I./inc/fit
OBJ_DIR=obj
SRC_DIR=src
BIN_DIR=.

.SECONDARY:

PROGS := $(wildcard *.cpp)
PROG_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.main.o,$(PROGS))
PROG_MAINS := $(patsubst %.cpp,$(BIN_DIR)/%.out,$(PROGS))


CPP_FILES := $(wildcard src/*/*.cpp)
CPP_FILES += $(wildcard src/*.cpp)
#CU_FILES := $(wildcard src/*/*.cu)
#CU_FILES += $(wildcard src/*.cu)

CLASS_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_FILES))
#CU_OBJS := $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.cu.o,$(CU_FILES))


#cuda objects
#$(OBJ_DIR)/%.cu.o : $(SRC_DIR)/%.cu 
#	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(LIB_CUDA)  -o $@ -c $<

#cpp class objects
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	 $(CXX) $(CXXFLAGS) $(INCLUDES) $(LIB_NETCDF) $(LIB_CGAL) -o $@ -c $<

#program objects
$(OBJ_DIR)/%.main.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIB_NETCDF) $(LIB_CGAL) -o $@ -c $<

#Programs
%.out: $(OBJ_DIR)/%.main.o $(CLASS_OBJS)
	 $(CXX) $(CXXFLAGS) $(INCLUDES) $(LIB_CGAL) $(LIB_NETCDF) -o $@ $+

#target rules

all:build

float: CXXFLAGS += -DSCALARFLOAT
float: build

debug: CXXFLAGS += -g -DCUDATHREADSYNC -DDEBUGFLAGUP
#debug: NVCCFLAGS += -g -lineinfo -Xptxas --generate-line-info -DDEBUGFLAGUP # -G note that in debug mode noise will always be reproducible
debug: build
build: $(CLASS_OBJS) $(PROG_MAINS) $(PROGS)

clean: 
	rm -f $(PROG_OBJS) $(CLASS_OBJS) $(PROG_OBJS) $(PROG_MAINS)

print-%  : ; @echo $* = $($*)
