#pragma once

#include "Macros.h"

#define Xp (0)
#define Yp (1)
#define Zp (2)
#define Tp (3)
#define Xm (4)
#define Ym (5)
#define Zm (6)
#define Tm (7)



#ifdef VGPU
#inclue "WilsonKernelsHandGpu.h"

template<class Simd>
void dslash_kernel(int nrep,Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  for(int i=0;i<nrep;i++){
    dslash_kernel_gpu(Up,outp,inp,nbr,nsite,Ls,prm);
  }
}
#else
#include "WilsonKernelsHandCpu.h"

template<class Simd>
void dslash_kernel(int nrep,Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  dslash_kernel_cpu(nrep,Up,outp,inp,nbr,nsite,Ls,prm);
}
#endif
