
SIMPLEDATA := arch/sse/static_data.cc

#OMP:=-fopenmp -std=c++11 -DOMP
OMP:=-std=c++11 -O3

#CXX       := mpicxx-openmpi-devel-clang40
#CXX       := mpiicpc
#CXX       := g++
CXX       := mpicxx
CXXFLAGS  := $(OMP)

AVX512_DATA   := arch/avx512/static_data.cc
AVX2_DATA     := arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc
AVX_DATA      := arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc
SSE_DATA      := arch/sse/static_data.cc

#############################################
# Intel
#############################################
#AVX512_CXXFLAGS  := -DAVX512 -xcore-avx512 $(OMP)
#AVX2_CXXFLAGS    := -DAVX2  -march=core-avx2 -xcore-avx2 $(OMP)
#AVX_CXXFLAGS     := -DAVX1  -mavx -xavx $(OMP)
#SSE_CXXFLAGS     := -DSSE4  -msse4.2 -xsse4.2  $(OMP)

#############################################
# CLANG
#############################################
AVX512_CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 $(OMP)
AVX2_CXXFLAGS    := -DAVX2  -mavx2 -mfma $(OMP)
AVX_CXXFLAGS     := -DAVX1  -mavx $(OMP)
SSE_CXXFLAGS     := -DSSE4  -msse4.2  $(OMP)

#############################################
# G++
#############################################
#AVX512_CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 $(OMP)
#AVX2_CXXFLAGS    := -DAVX2  -mavx2 -mfma $(OMP)
#AVX_CXXFLAGS     := -DAVX1  -mavx $(OMP)
#SSE_CXXFLAGS     := -DSSE4  -msse4.2  $(OMP)

#############################################
# TriSYCL
#############################################
CL_CXX := g++-mp-7
CL_CXXFLAGS := -I/Users/ayamaguc/Grid/triSYCL-master/include/ -std=c++1z -DOPENCL -DVGPU -DGEN_SIMD_WIDTH=64
CL_LDLIBS   := -lm
CL_LDFLAGS  := 


#Generic options
GENERIC_CXXFLAGS  := -DGEN -O3 -DGEN_SIMD_WIDTH=16 $(OMP)
GENERIC_DATA      := arch/sse/static_data.cc
triSYCL_CXXFLAGS  := -I/Users/ayamaguc/Grid/triSYCL-master/include 

################################################################################
# NVCC and gpu ; 512 bit vector coalescing
################################################################################

# VOLTA
GPUARCH    := --relocatable-device-code=true -gencode arch=compute_70,code=sm_70 

# PASCAL
#GPUARCH    := --relocatable-device-code=true -gencode arch=compute_60,code=sm_60 

GPUCC      := nvcc 
GPULINK    := nvcc $(GPUARCH)
GPU_CXXFLAGS  := -x cu -DVGPU -DGEN_SIMD_WIDTH=64 -I. -O3 -ccbin g++ -std=c++11 --expt-relaxed-constexpr --expt-extended-lambda $(GPUARCH)  -Xcompiler -fno-strict-aliasing
GPU_LDFLAGS  := -link -ccbin g++
GPU_DATA      := arch/avx512/static_data.cc
################################################################################
LDLIBS    := -lm
LDFLAGS   := 

all: bench.avx512 bench.avx2 bench.avx bench.sse bench.gen bench.simple parallel_vector_add.triSYCL bench.triSYCL

bench.avx512: bench.cc $(AVX512_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(AVX512_CXXFLAGS) bench.cc $(AVX512_DATA) $(LDLIBS) $(LDFLAGS) -o bench.avx512

bench.avx2: bench.cc $(AVX2_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(AVX2_CXXFLAGS) bench.cc $(AVX2_DATA) $(LDLIBS) $(LDFLAGS) -o bench.avx2

bench.avx: bench.cc $(AVX_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(AVX_CXXFLAGS) bench.cc $(AVX_DATA) $(LDLIBS) $(LDFLAGS) -o bench.avx

bench.sse: bench.cc $(SSE_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(SSE_CXXFLAGS) bench.cc $(SSE_DATA) $(LDLIBS) $(LDFLAGS) -o bench.sse

bench.gen: bench.cc $(GENERIC_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(GENERIC_CXXFLAGS) bench.cc $(GENERIC_DATA) $(LDLIBS) $(LDFLAGS) -o bench.gen

#	nvcc -x cu -DVGPU -DGEN_SIMD_WIDTH=64 -I. -O3 -ccbin g++ -std=c++11 -Xcompiler -Wno-deprecated-gpu-targets --expt-relaxed-constexpr --expt-extended-lambda --relocatable-device-code=true -gencode arch=compute_60,code=sm_60  -Xcompiler -fno-strict-aliasing -c bench.cc -o bench.gpu.o 
bench.gpu: bench.cc $(GPU_DATA)  WilsonKernelsHand.h Makefile
	$(GPUCC) $(GPU_CXXFLAGS) -c bench.cc -o bench.gpu.o
	$(GPUCC) $(GPU_CXXFLAGS) -c $(GPU_DATA) -o data.gpu.o
	$(GPULINK) $(GPU_LDFLAGS) bench.gpu.o data.gpu.o -o bench.gpu $(LDLIBS) $(LDFLAGS)

bench.simple: bench_simple.cc $(SIMPLEDATA) dslash_simple.h Makefile
	$(CXX) $(CXXFLAGS) bench_simple.cc $(SIMPLEDATA) $(LDLIBS) $(LDFLAGS) -o bench.simple

bench.triSYCL: bench.cc $(GPU_DATA) WilsonKernelsHand.h Makefile
	$(CL_CXX) $(CL_CXXFLAGS) bench.cc -o bench.triSYCL.o 
	$(CL_CXX) $(CL_CXXFLAGS) -c $(GPU_DATA) -o data.triSYCL.o
	$(CL_CXX) $(CL_LDFLAGS) bench.triSYCL.o data.triSYCL.o -o bench.triSYCL $(CL_LDLIBS)

parallel_vector_add.triSYCL: parallel_vector_add.cpp  Makefile
	$(CL_CXX) $(CL_CXXFLAGS) parallel_vector_add.cpp $(CL_LDLIBS)  $(CL_LDFLAGS) -o parallel_vector_add.triSYCL


clean:
	rm -f  bench.avx512 bench.avx2 bench.avx bench.sse bench.gen  bench.simple  bench.triSYCL TableGenerate
	rm -rf  *.dSYM*
	rm -f  *~

