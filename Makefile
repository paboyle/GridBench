
CXX       := mpicxx-openmpi-devel-clang60
CXXFLAGS  := -DAVX512 -mavx512f -mavx512pf -mavx512er -mavx512cd -O3 

#### TODO ### Implement other SIMD targets and verify
#CXXFLAGS  := -DAVX2 -mavx2 -O3 
#CXXFLAGS  := -DAVX -mavx -O3 
#CXXFLAGS  := -DSSE4 -msse4 -O3 

LDLIBS    := -lm
LDFLAGS   := 

all: bench

bench: bench.cpp static_data.o  WilsonKernelsHand.h
	$(CXX) $(CXXFLAGS) bench.cpp static_data.o $(LDLIBS) $(LDFLAGS) -o bench

static_data.o: arch/avx512/static_data.cc
	$(CXX) $(CXXFLAGS) -c  arch/avx512/static_data.cc -o static_data.o

clean:
	rm -f $(target) $(objects) $(depfile) *~

