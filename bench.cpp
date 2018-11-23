
// Invoke dslash.s - test for compiler-gsnerated code
#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include <math.h>

#include "Simd.h"

#ifdef SSE4
#warning "Including SSE static data"
#include "arch/sse4/static_data.h"
#endif

#if defined(AVX1) || defined (AVXFMA) || defined(AVX2) || defined(AVXFMA4)
#warning "Including AVX static data"
#include "arch/avx/static_data.h"
#endif

#ifdef AVX512
#warning "Including AVX512 static data"
#include "arch/avx512/static_data.h"
#endif

#include "WilsonKernelsHand.h"

#define  FMT std::dec
int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////////////
  // Option 2: copy from static arrays
  ////////////////////////////////////////////////////////////////////
  uint64_t umax   = nsite*18*8  *vComplexD::Nsimd(); 
  uint64_t fmax   = nsite*24*Ls*vComplexD::Nsimd(); 
  uint64_t nbrmax = nsite*Ls*8;

  std::cout << "umax " << umax<<std::endl;

  Vector<double> U(umax);   bcopy(U_static,&U[0],umax*sizeof(double));
  Vector<double> Psi(fmax);
  Vector<double> Phi(fmax);     bcopy(Phi_static,&Phi[0],fmax*sizeof(double));
  Vector<double> Psi_cpp(fmax); bcopy(Psi_cpp_static,&Psi_cpp[0],fmax*sizeof(double));
  uint64_t *nbr    = new uint64_t[nsite*Ls*8]; bcopy(nbr_static,nbr,nbrmax*sizeof(uint64_t));
  uint8_t  *prm    = new uint8_t[nsite*Ls*8]; bcopy(prm_static,prm,nbrmax*sizeof(uint8_t));
  //  uint64_t *lo     = new uint64_t[nsite];      bcopy(lo_static,lo,nsite*sizeof(uint64_t));

  std::cout << " calling dslash_kernel "<<std::endl;
  dslash_kernel<vComplexD>((vComplexD *)&U[0],
			   (vComplexD *)&Psi[0],
			   (vComplexD *)&Phi[0],
			   &nbr[0],
			   nsite,
			   Ls,
			   &prm[0]);
  
  double err=0;
  for(uint64_t i=0; i<fmax;i++){
    err += pow(Psi_cpp[i]-Psi[i],2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  
  return 0;
}
