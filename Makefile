
SIMPLEDATA := arch/sse/static_data.cc

#OMP:=-fopenmp -std=c++11
OMP:=-std=c++11

CXX       := mpicxx-openmpi-devel-clang40
CXXFLAGS  := -O3 $(OMP)

#AVX512 options
AVX512_CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 $(OMP)
AVX512_DATA      := arch/avx512/static_data.cc

#AVX2 options
AVX2_CXXFLAGS  := -DAVX2 -O3 -mavx2 -mfma $(OMP)
AVX2_DATA      :=  arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc

#AVX options
AVX_CXXFLAGS  := -DAVX1 -O3 -mavx $(OMP)
AVX_DATA      := arch/avx/static_data_gauge.cc arch/avx/static_data_fermion.cc

#SSE4 options
SSE_CXXFLAGS  := -DSSE4 -O3 -msse4.2  $(OMP)
SSE_DATA      := arch/sse/static_data.cc

#Generic options
GENERIC_CXXFLAGS  := -DGEN -O3 -DGEN_SIMD_WIDTH=16 $(OMP)
GENERIC_DATA      := arch/sse/static_data.cc


LDLIBS    := -lm
LDFLAGS   := 

all: bench.avx512 bench.avx2 bench.avx bench.sse bench.gen bench.simple

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

bench.simple: bench_simple.cc $(SIMPLEDATA) dslash_simple.h Makefile
	$(CXX) $(CXXFLAGS) bench_simple.cc $(SIMPLEDATA) $(LDLIBS) $(LDFLAGS) -o bench.simple

clean:
	rm -f  bench.avx512 bench.avx2 bench.avx bench.sse bench.gen  bench.simple TableGenerate
	rm -rf  *.dSYM*
	rm -f  *~

