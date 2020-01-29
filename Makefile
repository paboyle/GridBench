
SIMPLEDATA := arch/sse/static_data.cc

OMP:=-std=c++11  -O3 -DUSE_SYCL
#OMP:=-std=c++11 -O3

#CXX       := mpicxx-openmpi-devel-clang40
#CXX       := mpiicpc
#CXX       := g++-7
CXX        := dpcpp
#CXX       := g++
#CXX       := mpicxx
#CXX       := clang++-mp-6.0
#CXXCL     := clang++-mp-7.0
#CXXCL     := ${HOME}/QCD/build/bin/clang++ 
CXXCL     := dpcpp
CXXFLAGSCL:= 
LDFLAGSCL:= 
CXXFLAGS  := $(OMP)

AVX512_DATA   := arch/avx512/static_data.cc
AVX2_DATA     := arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc
AVX_DATA      := arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc
SSE_DATA      := arch/sse/static_data.cc
RRII_DATA     := arch/avx512/static_data.cc

#############################################
# Intel
#############################################
#AVX512_CXXFLAGS  := -DAVX512 -xcore-avx512 $(OMP)
#AVX2_CXXFLAGS    := -DAVX2  -march=core-avx2 -xcore-avx2 $(OMP)
#AVX_CXXFLAGS     := -DAVX1  -mavx -xavx $(OMP)
#SSE_CXXFLAGS     := -DSSE4  -msse4.2 -xsse4.2  $(OMP)
#RRII_CXXFLAGS     := -DRRII  -mavx2 -mfma  $(OMP)

#############################################
# CLANG
#############################################
AVX512_CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 $(OMP)
AVX2_CXXFLAGS    := -DAVX2  -mavx2 -mfma $(OMP)
AVX_CXXFLAGS     := -DAVX1  -mavx $(OMP)
SSE_CXXFLAGS     := -DSSE4  -msse4.2  $(OMP)
RRII_CXXFLAGS     := -DRRII  -mavx2 -mfma  $(OMP) -DGEN_SIMD_WIDTH=32

#############################################
# G++
#############################################
#AVX512_CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 $(OMP)
#AVX2_CXXFLAGS    := -DAVX2  -mavx2 -mfma $(OMP)
#AVX_CXXFLAGS     := -DAVX1  -mavx $(OMP)
#SSE_CXXFLAGS     := -DSSE4  -msse4.2  $(OMP)

#Generic options
GENERIC_CXXFLAGS  := -DGEN -O3 -DGEN_SIMD_WIDTH=16 $(OMP)
GENERIC_DATA      := arch/sse/static_data.cc

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

all: bench.avx512 bench.avx2 bench.avx bench.sse bench.gen bench.simple bench.rrii bench.sycl

bench.avx512: bench.cc $(AVX512_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(AVX512_CXXFLAGS) bench.cc $(AVX512_DATA) $(LDLIBS) $(LDFLAGS) -o bench.avx512

bench.avx2: bench.cc $(AVX2_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(AVX2_CXXFLAGS) bench.cc $(AVX2_DATA) $(LDLIBS) $(LDFLAGS) -o bench.avx2

bench.rrii: bench.cc $(RRII_DATA)  WilsonKernelsHand.h Makefile
	$(CXX) $(RRII_CXXFLAGS) bench.cc $(RRII_DATA) $(LDLIBS) $(LDFLAGS) -o bench.rrii

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
	$(CXX) $(CXXFLAGS) bench_simple.cc $(SIMPLEDATA) -I/usr/local/Cellar/boost/1.68.0_1/include -o bench.simple

bench.sycl: bench_sycl.cc $(SIMPLEDATA) dslash_simple.h Makefile
	$(CXXCL) -O3 -std=c++17 $(CXXFLAGSCL) bench_sycl.cc $(SIMPLEDATA) $(LDLIBS) $(LDFLAGS) $(LDFLAGSCL) -o bench.sycl

######################
# Build a test from triSYCL distro to check compiler working
######################
#
#parallel_vector_add:
#	clang++-mp-7.0  -std=c++17 -I/Users/ayamaguc/Grid/triSYCL-master/include parallel_vector_add.cpp  -I/usr/local/Cellar/boost/1.68.0_1/include
#

clean:
	rm -f  bench.avx512 bench.avx2 bench.avx bench.sse bench.gen  bench.simple TableGenerate bench.gpu bench.sycl
	rm -rf  *.dSYM*
	rm -f  *~

