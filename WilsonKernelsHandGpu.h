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
