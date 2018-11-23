
// Invoke dslash.s - test for compiler-gsnerated code
#include <stdio.h>
#include <vector>
#include <complex>
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include <math.h>


// 128 bit double precision
#include "AlignedAllocator.h"
#include "arch/sse/static_data.h"

typedef std::complex<double> ComplexD;

#include "dslash_simple.h"

#define  FMT std::dec
int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////////////
  // Option 2: copy from static arrays
  ////////////////////////////////////////////////////////////////////
  uint64_t umax   = nsite*18*8 ;
  uint64_t fmax   = nsite*24*Ls;
  uint64_t nbrmax = nsite*Ls*8;


  Vector<double> U(umax);   bcopy(U_static,&U[0],umax*sizeof(double));
  Vector<double> Psi(fmax);
  Vector<double> Phi(fmax);     bcopy(Phi_static,&Phi[0],fmax*sizeof(double));
  Vector<double> Psi_cpp(fmax); bcopy(Psi_cpp_static,&Psi_cpp[0],fmax*sizeof(double));
  uint64_t *nbr    = new uint64_t[nsite*Ls*8]; bcopy(nbr_static,nbr,nbrmax*sizeof(uint64_t));
  uint8_t  *prm    = new uint8_t[nsite*Ls*8]; bcopy(prm_static,prm,nbrmax*sizeof(uint8_t));

  std::cout << " calling dslash_kernel "<<std::endl;
  dslash_kernel<ComplexD>( (ComplexD *)&U[0],
			   (ComplexD *)&Psi[0],
			   (ComplexD *)&Phi[0],
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
