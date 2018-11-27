
// Invoke dslash.s - test for compiler-gsnerated code
#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include <math.h>
#include <chrono>
#include <cassert>

#include "Simd.h"

#ifdef GEN
#include "arch/sse/static_data.h"
#endif

#ifdef SSE4
#include "arch/sse/static_data.h"
#endif

#if defined(AVX1) || defined (AVXFMA) || defined(AVX2) || defined(AVXFMA4)
#include "arch/avx/static_data.h"
#endif

#ifdef AVX512
#include "arch/avx512/static_data.h"
#endif

#include "WilsonKernelsHand.h"


#define  FMT std::dec
int main(int argc, char* argv[])
{

  ////////////////////////////////////////////////////////////////////
  // Option 2: copy from static arrays
  ////////////////////////////////////////////////////////////////////
  uint64_t umax   = nsite*18*8 *vComplexD::Nsimd(); 
  uint64_t fmax   = nsite*24*Ls*vComplexD::Nsimd(); 
  uint64_t nbrmax = nsite*Ls*8;
  uint64_t vol    = nsite*Ls*vComplexD::Nsimd(); 

  Vector<double> U(umax);   bcopy(U_static,&U[0],umax*sizeof(double));
  Vector<double> Psi(fmax); bzero(&Psi[0],fmax*sizeof(double));
  Vector<double> Phi(fmax);     bcopy(Phi_static,&Phi[0],fmax*sizeof(double));
  Vector<double> Psi_cpp(fmax); bcopy(Psi_cpp_static,&Psi_cpp[0],fmax*sizeof(double));
  uint64_t *nbr    = new uint64_t[nsite*Ls*8]; bcopy(nbr_static,nbr,nbrmax*sizeof(uint64_t));
  uint8_t  *prm    = new uint8_t[nsite*Ls*8]; bcopy(prm_static,prm,nbrmax*sizeof(uint8_t));

  std::cout << std::endl;
  std::cout << "Calling dslash_kernel "<<std::endl;

  typedef  std::chrono::system_clock          Clock;
  typedef  std::chrono::time_point<Clock> TimePoint;
  typedef  std::chrono::microseconds          Usecs;

  Usecs elapsed;
  double flops = 1320.0*vol;
  int nrep=300; // cache warm
  TimePoint start = Clock::now();
  for(int i=0;i<nrep;i++){
    dslash_kernel<vComplexD>((vComplexD *)&U[0],
			     (vComplexD *)&Psi[0],
			     (vComplexD *)&Phi[0],
			     &nbr[0],
			     nsite,
			     Ls,
			     &prm[0]);
  }
  elapsed = std::chrono::duration_cast<Usecs>(Clock::now()-start);   
  std::cout << std::endl;
  std::cout <<"\t"<< nrep*flops/elapsed.count()/1000. << " Gflop/s in double precision; kernel call "<<elapsed.count()/nrep <<" microseconds "<<std::endl;
  std::cout << std::endl;

  // Check results
  double err=0;
  for(uint64_t i=0; i<fmax;i++){
    err += (Psi_cpp[i]-Psi[i])*(Psi_cpp[i]-Psi[i]);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);
  
  return 0;
}
