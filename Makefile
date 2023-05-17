CXX=g++
CXXFLAGS= -std=c++11 -Wall -O3 -march=native -Wno-unused-result -fopenmp
LIBS = -lm

NVCC = nvcc
NVCCFLAGS = -std=c++11 -Xcompiler -fopenmp

SRCS = $(wildcard *.cpp) $(wildcard *.cu)
TARGETS = $(basename $(wildcard *.cpp)) $(basename $(wildcard *.cu))
all : $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@

%:%.cu
	$(NVCC) $(NVCCFLAGS) $< -o $@

.PHONY: all clean

clean:
	rm -f $(TARGETS)
