
// Invoke dslash.s - test for compiler-gsnerated code
#include <stdio.h>
#include <vector>
#include <complex>
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include <math.h>
#include <chrono>
#include <cassert>


// 128 bit double precision
#include "AlignedAllocator.h"
#include "arch/sse/static_data.h"
//#define RESTRICT __restrict
#define RESTRICT 

template<class Double>
class myComplex 
{
public:
  Double re;
  Double im;
  // constructors
  template<class Floating> inline myComplex(Floating _re,Floating _im) : re(_re), im(_im) {};
  inline myComplex() {};

  inline myComplex operator - () const
  {
    return myComplex(-this->re,-this->im);
  }
  inline myComplex & operator += (const myComplex & RESTRICT r) 
  {
    this->re += r.re;
    this->im += r.im;
    return *this;
  }
  inline myComplex & operator -= (const myComplex & RESTRICT r) 
  {
    this->re -= r.re;
    this->im -= r.im;
    return *this;
  }
};
template<class Double>
inline myComplex<Double> operator * (const myComplex<Double> & RESTRICT l,const myComplex<Double> &RESTRICT r)
{
  return myComplex<Double>(l.re*r.re-l.im*r.im,l.re*r.im+l.im*r.re) ;
}
template<class Double>
inline myComplex<Double> operator + (const myComplex<Double> & RESTRICT l,const myComplex<Double> &RESTRICT r)
{
  return myComplex<Double>(l.re+r.re,l.im+r.im) ;
}
template<class Double>
inline myComplex<Double> operator - (const myComplex<Double> & RESTRICT l,const myComplex<Double> &RESTRICT r)
{
  return myComplex<Double>(l.re-r.re,l.im-r.im) ;
}

//typedef std::complex<double> ComplexD;
typedef myComplex<double> ComplexD;

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
  uint64_t vol    = nsite*Ls; 

  Vector<double> U(umax);   bcopy(U_static,&U[0],umax*sizeof(double));
  Vector<double> Psi(fmax);
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
    dslash_kernel<ComplexD>((ComplexD *)&U[0],
			    (ComplexD *)&Psi[0],
			    (ComplexD *)&Phi[0],
			    &nbr[0],
			    nsite,
			    Ls,
			    &prm[0]);
  }
  elapsed = std::chrono::duration_cast<Usecs>(Clock::now()-start);   
  std::cout <<std::endl;
  std::cout <<"\t"<< nrep*flops/elapsed.count()/1000. << " Gflop/s in double precision; kernel call "<<elapsed.count()/nrep <<" microseconds "<<std::endl;
  std::cout <<std::endl;

  // Check results
  double err=0;
  for(uint64_t i=0; i<fmax;i++){
    err += pow(Psi_cpp[i]-Psi[i],2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);

  std::cout <<std::endl;
  std::cout << "Calling dslash_kernel_unroll "<<std::endl;

  start = Clock::now();
  for(int i=0;i<nrep;i++){
    dslash_kernel_unroll<ComplexD>((ComplexD *)&U[0],
				   (ComplexD *)&Psi[0],
				   (ComplexD *)&Phi[0],
				   &nbr[0],
				   nsite,
				   Ls,
				   &prm[0]);
  }
  elapsed = std::chrono::duration_cast<Usecs>(Clock::now()-start);   
  std::cout <<std::endl;
  std::cout <<"\t"<< nrep*flops/elapsed.count()/1000. << " Gflop/s in double precision; kernel call "<<elapsed.count()/nrep <<" microseconds "<<std::endl;
  std::cout <<std::endl;

  // Check results
  err=0;
  for(uint64_t i=0; i<fmax;i++){
    err += pow(Psi_cpp[i]-Psi[i],2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);

  
  return 0;
}
