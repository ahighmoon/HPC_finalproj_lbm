CXX=g++
CXXFLAGS= -std=c++11 -Wall -O3 -march=native -Wno-unused-result -fopenmp
LIBS = -lm

NVCC = nvcc
NVCCFLAGS = -std=c++11 -Xcompiler -fopenmp
#CUDA_INCDIR = -I $(CUDA_HOME)/include -I $(CUDA_HOME)/samples/common/inc
#CUDA_LIBS = -lblas -L${CUDA_HOME}/lib64 
#-lcudart -lcublas

SRCS = $(wildcard *.cpp) $(wildcard *.cu)
TARGETS = $(basename $(wildcard *.cpp)) $(basename $(wildcard *.cu))
#TARGETS = $(SRCS:%=%.out)
all : $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@

%:%.cu
	$(NVCC) $(NVCCFLAGS) $< -o $@

.PHONY: all clean

clean:
	rm -f $(TARGETS)
