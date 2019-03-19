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

#define LOAD_CHIMU \
  {const SiteSpinor & ref (in[offset]);	\
    Chimu_00=ref[0][0];\
    Chimu_01=ref[0][1];\
    Chimu_02=ref[0][2];\
    Chimu_10=ref[1][0];\
    Chimu_11=ref[1][1];\
    Chimu_12=ref[1][2];\
    Chimu_20=ref[2][0];\
    Chimu_21=ref[2][1];\
    Chimu_22=ref[2][2];\
    Chimu_30=ref[3][0];\
    Chimu_31=ref[3][1];\
    Chimu_32=ref[3][2];}


// To splat or not to splat depends on the implementation
#define MULT_2SPIN(A)\
  {auto & ref(U[sU][A]);			\
    U_00=ref[0][0];				\
    U_10=ref[1][0];				\
    U_20=ref[2][0];				\
    U_01=ref[0][1];				\
    U_11=ref[1][1];				\
    U_21=ref[2][1];				\
    UChi_00 = U_00*Chi_00;\
    UChi_10 = U_00*Chi_10;\
    UChi_01 = U_10*Chi_00;\
    UChi_11 = U_10*Chi_10;\
    UChi_02 = U_20*Chi_00;\
    UChi_12 = U_20*Chi_10;\
    UChi_00+= U_01*Chi_01;\
    UChi_10+= U_01*Chi_11;\
    UChi_01+= U_11*Chi_01;\
    UChi_11+= U_11*Chi_11;\
    UChi_02+= U_21*Chi_01;\
    UChi_12+= U_21*Chi_11;\
    U_00=ref[0][2];	\
    U_10=ref[1][2];	\
    U_20=ref[2][2];	\
    UChi_00+= U_00*Chi_02;\
    UChi_10+= U_00*Chi_12;\
    UChi_01+= U_10*Chi_02;\
    UChi_11+= U_10*Chi_12;\
    UChi_02+= U_20*Chi_02;\
    UChi_12+= U_20*Chi_12;}


#define PERMUTE_DIR(dir)			\
      permute##dir(Chi_00,Chi_00);\
      permute##dir(Chi_01,Chi_01);\
      permute##dir(Chi_02,Chi_02);\
      permute##dir(Chi_10,Chi_10);\
      permute##dir(Chi_11,Chi_11);\
      permute##dir(Chi_12,Chi_12);

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJ \
    Chi_00 = Chimu_00+timesI(Chimu_30);\
    Chi_01 = Chimu_01+timesI(Chimu_31);\
    Chi_02 = Chimu_02+timesI(Chimu_32);\
    Chi_10 = Chimu_10+timesI(Chimu_20);\
    Chi_11 = Chimu_11+timesI(Chimu_21);\
    Chi_12 = Chimu_12+timesI(Chimu_22);

#define YP_PROJ \
    Chi_00 = Chimu_00-Chimu_30;\
    Chi_01 = Chimu_01-Chimu_31;\
    Chi_02 = Chimu_02-Chimu_32;\
    Chi_10 = Chimu_10+Chimu_20;\
    Chi_11 = Chimu_11+Chimu_21;\
    Chi_12 = Chimu_12+Chimu_22;

#define ZP_PROJ \
  Chi_00 = Chimu_00+timesI(Chimu_20);		\
  Chi_01 = Chimu_01+timesI(Chimu_21);		\
  Chi_02 = Chimu_02+timesI(Chimu_22);		\
  Chi_10 = Chimu_10-timesI(Chimu_30);		\
  Chi_11 = Chimu_11-timesI(Chimu_31);		\
  Chi_12 = Chimu_12-timesI(Chimu_32);

#define TP_PROJ \
  Chi_00 = Chimu_00+Chimu_20;		\
  Chi_01 = Chimu_01+Chimu_21;		\
  Chi_02 = Chimu_02+Chimu_22;		\
  Chi_10 = Chimu_10+Chimu_30;		\
  Chi_11 = Chimu_11+Chimu_31;		\
  Chi_12 = Chimu_12+Chimu_32;


//      hspin(0)=fspin(0)-timesI(fspin(3));
//      hspin(1)=fspin(1)-timesI(fspin(2));
#define XM_PROJ \
    Chi_00 = Chimu_00-timesI(Chimu_30);\
    Chi_01 = Chimu_01-timesI(Chimu_31);\
    Chi_02 = Chimu_02-timesI(Chimu_32);\
    Chi_10 = Chimu_10-timesI(Chimu_20);\
    Chi_11 = Chimu_11-timesI(Chimu_21);\
    Chi_12 = Chimu_12-timesI(Chimu_22);

#define YM_PROJ \
    Chi_00 = Chimu_00+Chimu_30;\
    Chi_01 = Chimu_01+Chimu_31;\
    Chi_02 = Chimu_02+Chimu_32;\
    Chi_10 = Chimu_10-Chimu_20;\
    Chi_11 = Chimu_11-Chimu_21;\
    Chi_12 = Chimu_12-Chimu_22;

#define ZM_PROJ \
  Chi_00 = Chimu_00-timesI(Chimu_20);		\
  Chi_01 = Chimu_01-timesI(Chimu_21);		\
  Chi_02 = Chimu_02-timesI(Chimu_22);		\
  Chi_10 = Chimu_10+timesI(Chimu_30);		\
  Chi_11 = Chimu_11+timesI(Chimu_31);		\
  Chi_12 = Chimu_12+timesI(Chimu_32);

#define TM_PROJ \
  Chi_00 = Chimu_00-Chimu_20;		\
  Chi_01 = Chimu_01-Chimu_21;		\
  Chi_02 = Chimu_02-Chimu_22;		\
  Chi_10 = Chimu_10-Chimu_30;		\
  Chi_11 = Chimu_11-Chimu_31;		\
  Chi_12 = Chimu_12-Chimu_32;

//      fspin(0)=hspin(0);
//      fspin(1)=hspin(1);
//      fspin(2)=timesMinusI(hspin(1));
//      fspin(3)=timesMinusI(hspin(0));
#define XP_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesMinusI(UChi_10);\
  result_21 = timesMinusI(UChi_11);\
  result_22 = timesMinusI(UChi_12);\
  result_30 = timesMinusI(UChi_00);\
  result_31 = timesMinusI(UChi_01);\
  result_32 = timesMinusI(UChi_02);

#define XP_RECON_ACCUM\
  result_00+=UChi_00;\
  result_01+=UChi_01;\
  result_02+=UChi_02;\
  result_10+=UChi_10;\
  result_11+=UChi_11;\
  result_12+=UChi_12;\
  result_20-=timesI(UChi_10);\
  result_21-=timesI(UChi_11);\
  result_22-=timesI(UChi_12);\
  result_30-=timesI(UChi_00);\
  result_31-=timesI(UChi_01);\
  result_32-=timesI(UChi_02);

#define XM_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesI(UChi_10);\
  result_21 = timesI(UChi_11);\
  result_22 = timesI(UChi_12);\
  result_30 = timesI(UChi_00);\
  result_31 = timesI(UChi_01);\
  result_32 = timesI(UChi_02);

#define XM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_10);\
  result_21+= timesI(UChi_11);\
  result_22+= timesI(UChi_12);\
  result_30+= timesI(UChi_00);\
  result_31+= timesI(UChi_01);\
  result_32+= timesI(UChi_02);

#define YP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_10;\
  result_21+= UChi_11;\
  result_22+= UChi_12;\
  result_30-= UChi_00;\
  result_31-= UChi_01;\
  result_32-= UChi_02;

#define YM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_10;\
  result_21-= UChi_11;\
  result_22-= UChi_12;\
  result_30+= UChi_00;\
  result_31+= UChi_01;\
  result_32+= UChi_02;

#define ZP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= timesI(UChi_00);			\
  result_21-= timesI(UChi_01);			\
  result_22-= timesI(UChi_02);			\
  result_30+= timesI(UChi_10);			\
  result_31+= timesI(UChi_11);			\
  result_32+= timesI(UChi_12);

#define ZM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_00);			\
  result_21+= timesI(UChi_01);			\
  result_22+= timesI(UChi_02);			\
  result_30-= timesI(UChi_10);			\
  result_31-= timesI(UChi_11);			\
  result_32-= timesI(UChi_12);

#define TP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_00;			\
  result_21+= UChi_01;			\
  result_22+= UChi_02;			\
  result_30+= UChi_10;			\
  result_31+= UChi_11;			\
  result_32+= UChi_12;

#define TM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_00;	\
  result_21-= UChi_01;	\
  result_22-= UChi_02;	\
  result_30-= UChi_10;	\
  result_31-= UChi_11;	\
  result_32-= UChi_12;

#define HAND_STENCIL_LEG(PROJ,PERM,DIR,RECON)		\
  offset = nbr[ss*8+DIR];				\
  perm   = prm[ss*8+DIR];				\
  LOAD_CHIMU;						\
  PROJ;							\
  if ( perm) {						\
    PERMUTE_DIR(PERM);					\
  }							\
  MULT_2SPIN(DIR);					\
  RECON;					

#define HAND_RESULT(ss)				\
  {						\
    SiteSpinor & ref (out[ss]);			\
    vstream(ref[0][0],result_00);		\
    vstream(ref[0][1],result_01);		\
    vstream(ref[0][2],result_02);		\
    vstream(ref[1][0],result_10);		\
    vstream(ref[1][1],result_11);		\
    vstream(ref[1][2],result_12);		\
    vstream(ref[2][0],result_20);		\
    vstream(ref[2][1],result_21);		\
    vstream(ref[2][2],result_22);		\
    vstream(ref[3][0],result_30);		\
    vstream(ref[3][1],result_31);		\
    vstream(ref[3][2],result_32);		\
  }

#define HAND_DECLARATIONS(Simd)			\
  Simd result_00;				\
  Simd result_01;				\
  Simd result_02;				\
  Simd result_10;				\
  Simd result_11;				\
  Simd result_12;				\
  Simd result_20;				\
  Simd result_21;				\
  Simd result_22;				\
  Simd result_30;				\
  Simd result_31;				\
  Simd result_32;				\
  Simd Chi_00;					\
  Simd Chi_01;					\
  Simd Chi_02;					\
  Simd Chi_10;					\
  Simd Chi_11;					\
  Simd Chi_12;					\
  Simd UChi_00;					\
  Simd UChi_01;					\
  Simd UChi_02;					\
  Simd UChi_10;					\
  Simd UChi_11;					\
  Simd UChi_12;					\
  Simd U_00;					\
  Simd U_10;					\
  Simd U_20;					\
  Simd U_01;					\
  Simd U_11;					\
  Simd U_21;

#define ZERO_RESULT				\
  result_00=zero;				\
  result_01=zero;				\
  result_02=zero;				\
  result_10=zero;				\
  result_11=zero;				\
  result_12=zero;				\
  result_20=zero;				\
  result_21=zero;				\
  result_22=zero;				\
  result_30=zero;				\
  result_31=zero;				\
  result_32=zero;			

#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_12
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_12

template<class Simd>
void dslash_kernel_cpu(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  typedef Simd SiteSpinor[4][3];
  typedef Simd SiteHalfSpinor[2][3];
  typedef Simd SiteDoubledGaugeField[8][3][3];

  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  SiteSpinor *out = (SiteSpinor *) outp;
  SiteSpinor *in  = (SiteSpinor *) inp;
  SiteDoubledGaugeField *U  = (SiteDoubledGaugeField *) Up;

  Simd complex_i;  vsplat(complex_i, S(0.0, 1.0));

#ifdef OMP
#pragma omp parallel for
#endif
  for(uint64_t ssite=0;ssite<nsite;ssite++){

    HAND_DECLARATIONS(Simd);

    int offset,perm;
    uint64_t sU = ssite;
    uint64_t ss = sU*Ls;
    for(uint64_t s=0;s<Ls;s++){
      HAND_STENCIL_LEG(XM_PROJ,3,Xp,XM_RECON);
      HAND_STENCIL_LEG(YM_PROJ,2,Yp,YM_RECON_ACCUM);
      HAND_STENCIL_LEG(ZM_PROJ,1,Zp,ZM_RECON_ACCUM);
      HAND_STENCIL_LEG(TM_PROJ,0,Tp,TM_RECON_ACCUM);
      HAND_STENCIL_LEG(XP_PROJ,3,Xm,XP_RECON_ACCUM);
      HAND_STENCIL_LEG(YP_PROJ,2,Ym,YP_RECON_ACCUM);
      HAND_STENCIL_LEG(ZP_PROJ,1,Zm,ZP_RECON_ACCUM);
      HAND_STENCIL_LEG(TP_PROJ,0,Tm,TP_RECON_ACCUM);
      HAND_RESULT(ss);
      ss++;
    }
  }
}


#define HAND_STENCIL_LEG_GPU(PROJ,PERM,DIR,RECON)			\
  {offset = nbr[ss*8+DIR];						\
  perm   = prm[ss*8+DIR];						\
  int mask = Nsimd >> (PERM + 1);					\
  int plane= perm ? (lane ^ mask) : lane;				\
  LOAD_CHIMU_GPU;							\
  PROJ;									\
  MULT_2SPIN_GPU(DIR);							\
  RECON;}					

#define LOAD_CHIMU_GPU \
  {const SiteSpinor & ref (in[offset]);	\
    Chimu_00=ref[0][0].v[plane];\
    Chimu_01=ref[0][1].v[plane];\
    Chimu_02=ref[0][2].v[plane];\
    Chimu_10=ref[1][0].v[plane];\
    Chimu_11=ref[1][1].v[plane];\
    Chimu_12=ref[1][2].v[plane];\
    Chimu_20=ref[2][0].v[plane];\
    Chimu_21=ref[2][1].v[plane];\
    Chimu_22=ref[2][2].v[plane];\
    Chimu_30=ref[3][0].v[plane];\
    Chimu_31=ref[3][1].v[plane];\
    Chimu_32=ref[3][2].v[plane];}


// To splat or not to splat depends on the implementation
#define MULT_2SPIN_GPU(A)\
  {auto & ref(U[sU][A]);			\
    U_00=ref[0][0].v[lane];				\
    U_10=ref[1][0].v[lane];				\
    U_20=ref[2][0].v[lane];				\
    U_01=ref[0][1].v[lane];				\
    U_11=ref[1][1].v[lane];				\
    U_21=ref[2][1].v[lane];				\
    UChi_00 = U_00*Chi_00;\
    UChi_10 = U_00*Chi_10;\
    UChi_01 = U_10*Chi_00;\
    UChi_11 = U_10*Chi_10;\
    UChi_02 = U_20*Chi_00;\
    UChi_12 = U_20*Chi_10;\
    UChi_00+= U_01*Chi_01;\
    UChi_10+= U_01*Chi_11;\
    UChi_01+= U_11*Chi_01;\
    UChi_11+= U_11*Chi_11;\
    UChi_02+= U_21*Chi_01;\
    UChi_12+= U_21*Chi_11;\
    U_00=ref[0][2].v[lane];	\
    U_10=ref[1][2].v[lane];	\
    U_20=ref[2][2].v[lane];	\
    UChi_00+= U_00*Chi_02;\
    UChi_10+= U_00*Chi_12;\
    UChi_01+= U_10*Chi_02;\
    UChi_11+= U_10*Chi_12;\
    UChi_02+= U_20*Chi_02;\
    UChi_12+= U_20*Chi_12;}


#define HAND_RESULT_GPU(ss)				\
  {						\
    SiteSpinor & ref (out[ss]);			\
    ref[0][0].v[lane]=result_00;		\
    ref[0][1].v[lane]=result_01;		\
    ref[0][2].v[lane]=result_02;		\
    ref[1][0].v[lane]=result_10;		\
    ref[1][1].v[lane]=result_11;		\
    ref[1][2].v[lane]=result_12;		\
    ref[2][0].v[lane]=result_20;		\
    ref[2][1].v[lane]=result_21;		\
    ref[2][2].v[lane]=result_22;		\
    ref[3][0].v[lane]=result_30;		\
    ref[3][1].v[lane]=result_31;		\
    ref[3][2].v[lane]=result_32;		\
  }

template<class Simd>
accelerator_inline void dslash_kernel_gpu_site(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint8_t *prm,int Ls,uint64_t ssite,uint64_t s)
{
  typedef typename Simd::scalar_type scalar_type;
  typedef typename Simd::vector_type vector_type;

  // Live in memory in GPU
  typedef vector_type SiteSpinor[4][3];
  typedef vector_type SiteHalfSpinor[2][3];
  typedef vector_type SiteDoubledGaugeField[8][3][3];

  // Live in register file in GPU
  typedef scalar_type ScalarSpinor[4][3];
  typedef scalar_type ScalarHalfSpinor[2][3];
  typedef scalar_type ScalarGaugeLink[3][3];

  const int Nsimd = Simd::Nsimd();

  uint64_t lane_offset= get_my_lane_offset(Nsimd);
  uint64_t lanes      = get_my_lanes(Nsimd);

  // Memory references
  SiteSpinor *out = (SiteSpinor *) outp;
  SiteSpinor *in  = (SiteSpinor *) inp;
  SiteDoubledGaugeField *U  = (SiteDoubledGaugeField *) Up;

  scalar_type complex_i(0.0,1.0);

  // CPU loop over lanes, GPU uses CUDA threads for each lane
#ifndef __CUDA_ARCH__
  for(int lane = lane_offset;lane<lane_offset+lanes;lane++){
#else
    int lane = lane_offset; 
#endif
    HAND_DECLARATIONS(scalar_type);
    int offset,perm;
    uint64_t sU = ssite;
    uint64_t ss = sU*Ls+s;
    HAND_STENCIL_LEG_GPU(XM_PROJ,3,Xp,XM_RECON);
    HAND_STENCIL_LEG_GPU(YM_PROJ,2,Yp,YM_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(ZM_PROJ,1,Zp,ZM_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(TM_PROJ,0,Tp,TM_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(XP_PROJ,3,Xm,XP_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(YP_PROJ,2,Ym,YP_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(ZP_PROJ,1,Zm,ZP_RECON_ACCUM);
    HAND_STENCIL_LEG_GPU(TP_PROJ,0,Tm,TP_RECON_ACCUM);
    HAND_RESULT_GPU(ss);
#ifndef __CUDA_ARCH__
  }
#endif
}

const uint32_t gpu_threads = 128;

#ifdef __NVCC__
template<class Simd>
void dslash_kernel_gpu(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  const uint64_t nsimd = Simd::Nsimd();
  const uint64_t    NN = nsite*Ls*nsimd;

  // Translates to kernel call to CUDA
  accelerator_loopN( sss, NN, {
      uint64_t cur  = sss;      cur = cur / nsimd;
      uint64_t   s  = cur%Ls;   cur = cur / Ls;
      uint64_t   sU = cur;        // 4d site
      dslash_kernel_gpu_site<Simd>(Up,outp,inp,nbr,prm,Ls,sU,s);
  });
}
#else

#ifdef VGPU
#include <CL/sycl.hpp>
//
// This is where it gets complicated
//
// Example code:
/*
  // Create a queue to work on
  cl::sycl::queue q;

  // Create buffers from a & b vectors
  cl::sycl::buffer<float> A { std::begin(a), std::end(a) };
  cl::sycl::buffer<float> B { std::begin(b), std::end(b) };

  {
    // A buffer of N float using the storage of c
    cl::sycl::buffer<float> C { c, N };

    using namespace cl::sycl;

    q.submit([&](handler &cgh) {
        // In the kernel A and B are read, but C is written
        auto ka = A.get_access<access::mode::read>(cgh);
        auto kb = B.get_access<access::mode::read>(cgh);
        auto kc = C.get_access<access::mode::write>(cgh);

        // Enqueue a parallel kernel
        cgh.parallel_for<class vector_add>(N,
                                           [=] (id<1> index) {
                                             kc[index] = ka[index] + kb[index];
                                           });
      }); //< End of our commands for this queue
  } //< Buffer C goes out of scope and copies back values to c
 */
// Requirements:
//
// For each host array must create a sycl::buffer
//     We must get an accessor to the buffer with "get_access".
//     
// This is "hard" because the the lambda is capturing its arguments. 
// 
// It is easier to constrain in the Grid ET code, because we do a transformation 
// from a Lattice<object> to a "LatticeView" object. This can set up an accessor 
// in principle when the expression object is created.
//
// However, for the dslash kernel this will need to be done manually.
//

template<class Simd>
void dslash_kernel_sycl(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  using namespace cl::sycl;

  const uint64_t nsimd = Simd::Nsimd();
  const uint64_t    NN = nsite*Ls*nsimd;

  uint64_t begin=0;
  uint64_t end  =NN;
  cl::sycl::queue q;
  cl::sycl::buffer<Simd>     Up_b   { &  Up[begin],&  Up[end]};
  cl::sycl::buffer<Simd>     inp_b  { & inp[begin],& inp[end]};
  cl::sycl::buffer<uint64_t> nbr_b  { & nbr[begin],& nbr[end]};
  cl::sycl::buffer<uint8_t>  prm_b  { & prm[begin],& prm[end]};

  {

    cl::sycl::buffer<Simd>     outp_b { &outp[begin],&outp[end]};

    q.submit([&](handler &cgh) {
        // In the kernel A and B are read, but C is written
        auto Up_k   =  Up_b.get_access<access::mode::read>(cgh);
        auto inp_k  = inp_b.get_access<access::mode::read>(cgh); 
        auto nbr_k  = nbr_b.get_access<access::mode::read>(cgh); 
        auto prm_k  = prm_b.get_access<access::mode::read>(cgh); 
        auto outp_k =outp_b.get_access<access::mode::write>(cgh); 

	cgh.parallel_for<class dslash>(NN,
				       [=] (id<1> sss) {
					 uint64_t cur  = sss;      cur = cur / nsimd;
					 uint64_t   s  = cur%Ls;   cur = cur / Ls;
					 uint64_t   sU = cur;        // 4d site
					 dslash_kernel_gpu_site<Simd>(&Up_k[0],
								      &outp_k[0],
								      &inp_k[0],
								      &nbr_k[0],
								      &prm_k[0],Ls,sU,s);
				       });
      }); //< End of our commands for this queue

  } //< Buffer outp_b goes out of scope and copies back values to outp
}
#endif
#endif

#ifdef VGPU
template<class Simd>
void dslash_kernel(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
#ifdef __NVCC__
  dslash_kernel_gpu(Up,outp,inp,nbr,nsite,Ls,prm);
#else
  dslash_kernel_sycl(Up,outp,inp,nbr,nsite,Ls,prm);
#endif
}
#else
template<class Simd>
void dslash_kernel(Simd *Up,Simd *outp,Simd *inp,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  dslash_kernel_cpu(Up,outp,inp,nbr,nsite,Ls,prm);
}
#endif
