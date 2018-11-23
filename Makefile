
SIMPLEDATA := arch/sse/static_data.cc

CXX       := mpicxx-openmpi-devel-clang60

#AVX512 options
CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 
DATA      := arch/avx512/static_data.cc

#AVX2 options
#CXXFLAGS  := -DAVX2 -O3 
#DATA      := arch/avx/static_data.cc

#SSE4 options
#CXXFLAGS  := -DAVX2 -O3 
#DATA      := arch/avx/static_data.cc

#Generic options
#CXXFLAGS  := -DGEN -O3 
#DATA      := arch/avx/static_data.cc

LDLIBS    := -lm
LDFLAGS   := 

all: bench bench_simple

bench: bench.cpp $(DATA)  WilsonKernelsHand.h
	$(CXX) $(CXXFLAGS) bench.cpp $(DATA) $(LDLIBS) $(LDFLAGS) -o bench

bench_simple: bench_simple.cpp $(SIMPLEDATA) dslash_simple.h
	$(CXX) $(CXXFLAGS) bench_simple.cpp $(SIMPLEDATA) $(LDLIBS) $(LDFLAGS) -o bench_simple

clean:
	rm -f bench bench_simple

