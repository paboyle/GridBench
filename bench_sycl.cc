

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

#include <CL/sycl.hpp>
using namespace cl::sycl;

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
  uint64_t umax   = nsite*9*8 ;
  uint64_t fmax   = nsite*12*Ls;
  uint64_t nbrmax = nsite*Ls*8;
  uint64_t vol    = nsite*Ls; 

  Vector<ComplexD> U(umax);   bcopy(U_static,&U[0],umax*sizeof(ComplexD));
  Vector<ComplexD> Psi(fmax);
  Vector<ComplexD> Phi(fmax);     bcopy(Phi_static,&Phi[0],fmax*sizeof(ComplexD));
  Vector<ComplexD> Psi_cpp(fmax); bcopy(Psi_cpp_static,&Psi_cpp[0],fmax*sizeof(ComplexD));
  uint64_t *nbr    = new uint64_t[nsite*Ls*8]; bcopy(nbr_static,nbr,nbrmax*sizeof(uint64_t));
  uint8_t  *prm    = new uint8_t[nsite*Ls*8]; bcopy(prm_static,prm,nbrmax*sizeof(uint8_t));
  uint64_t nbr_size=nsite*Ls*8;

  std::cout << std::endl;
  std::cout << "Calling dslash_kernel "<<std::endl;

  ComplexD zero(0,0);
  for(uint64_t i=0; i<fmax;i++){
    Psi[i]=zero;
  }

  typedef  std::chrono::system_clock          Clock;
  typedef  std::chrono::time_point<Clock> TimePoint;
  typedef  std::chrono::microseconds          Usecs;

  Usecs elapsed;
  Usecs kernel;
  double flops = 1320.0*vol;
  int nrep=300; // cache warm


  TimePoint kernel_start;
  TimePoint start = Clock::now();
  for(int i=0;i<nrep;i++){
    dslash_kernel<ComplexD>(&U[0],
			    &Psi[0],
			    &Phi[0],
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
    err += pow(Psi_cpp[i].re-Psi[i].re,2)+ pow(Psi_cpp[i].im-Psi[i].im,2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);
  for(uint64_t i=0; i<fmax;i++){
    Psi[i]=zero;
  }

  std::cout <<std::endl;
  std::cout << "Calling dslash_kernel_unroll "<<std::endl;

  start = Clock::now();
  for(int i=0;i<nrep;i++){
    dslash_kernel_unroll<ComplexD>(&U[0],
				   &Psi[0],
				   &Phi[0],
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
    err += pow(Psi_cpp[i].re-Psi[i].re,2)+ pow(Psi_cpp[i].im-Psi[i].im,2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);

  for(uint64_t i=0; i<fmax;i++){
    Psi[i]=zero;
  }
  ////////////////////////////////////////
  // Create a queue to work on 
  // SYCL call
  ////////////////////////////////////////
  std::cout <<std::endl;
  std::cout << "Calling dslash_kernel_SYCL "<<std::endl;

  cl::sycl::gpu_selector selector; cl::sycl::queue q(selector);
  //  cl::sycl::queue q;

  buffer<ComplexD> U_b   { std::begin(U), std::end(U) };
  buffer<ComplexD> Phi_b { std::begin(Phi), std::end(Phi) };
  buffer<uint64_t> nbr_b { &nbr[0], &nbr[nbr_size] };
  buffer<uint8_t > prm_b { &prm[0], &prm[nbr_size] };

  // result
  cl::sycl::range<1> nnsite{nsite};

  //  int nnsite=nsite;
  start = Clock::now();
  {
    buffer<ComplexD> Psi_b { &Psi[0],  fmax };
    
    kernel_start = Clock::now();
    for(int i=0;i<nrep;i++){

      q.submit([&](handler &cgh) {

	  auto U_k   = U_b.get_access<access::mode::read>(cgh);
	  auto Phi_k = Phi_b.get_access<access::mode::read>(cgh);
	  auto nbr_k = nbr_b.get_access<access::mode::read>(cgh);
	  auto prm_k = prm_b.get_access<access::mode::read>(cgh);
	  auto Psi_k = Psi_b.get_access<access::mode::write>(cgh);
	  
	// Enqueue a parallel kernel

	  cgh.parallel_for<class vector_add>(nnsite, [=] (id<1> index) {
	    
	    int site = index[0];
	    dslash_kernel_site(site,
			       U_k,
			       Psi_k,
			       Phi_k,
			       nbr_k,
			       nsite,
			       Ls,
			       prm_k);
	});
      }); //< End of our commands for this queue
    }
    q.wait();
    kernel  = std::chrono::duration_cast<Usecs>(Clock::now()-kernel_start);   
  }    // Trigger copy back
  elapsed = std::chrono::duration_cast<Usecs>(Clock::now()-start);   

  std::cout <<std::endl;
  std::cout <<"\t"<< nrep*flops/elapsed.count()/1000. << " Gflop/s in double precision; including data motion "<<elapsed.count()/nrep <<" microseconds "<<std::endl;
  std::cout <<"\t"<< nrep*flops/elapsed.count()/1000. << " Gflop/s in double precision; kernel calls "<<kernel.count()/nrep <<" microseconds "<<std::endl;
  std::cout <<std::endl;

  err=0;
  for(uint64_t i=0; i<nsite*24;i++){
    err += pow(Psi_cpp[i].re-Psi[i].re,2)+ pow(Psi_cpp[i].im-Psi[i].im,2);
  };
  std::cout<< "normdiff "<< err<<std::endl;
  assert(err <= 1.0e-10);

  return 0;
}

