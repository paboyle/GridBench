#pragma once

/////////////////////////////////////////////////
// Looping constructs for OpenMP and onload
// Can we add SyCL ??
/////////////////////////////////////////////////

#ifndef VGPU 

#define accelerator_inline inline
#define accelerator
#define thread_loop( range , ... )              for range { __VA_ARGS__ ; };
#define accelerator_loopN( iterator, num, ... ) thread_loop( (int iterator = 0;iterator<num;iterator++), { __VA_ARGS__ });

#else

////////////////////
// CUDA target
////////////////////

template<typename lambda>  __global__
void LambdaApply(uint64_t base, uint64_t Num, lambda Lambda)
{
  uint64_t ss = blockIdx.x*blockDim.x + threadIdx.x;
  if ( ss < Num ) {
    Lambda(ss+base);
  }
}

#define accelerator_inline __host__ __device__ inline
#define accelerator __host__ __device__
#define thread_loop( range , ... )              for range { __VA_ARGS__ ; };

#define accelerator_loopN_debug( iterator, num, ... ) thread_loop( (int iterator = 0;iterator<num;iterator++), { __VA_ARGS__ });

#define accelerator_loopN( iterator, num, ... )			\
  typedef decltype(num) Iterator;				\
  if ( num > 0 ) {			                        \
    auto lambda = [=] accelerator (Iterator iterator) mutable { \
      __VA_ARGS__;						\
    };								\
    Iterator base = 0;						\
    Iterator num_block  = (num+gpu_threads-1)/gpu_threads;	\
    LambdaApply<<<num_block,gpu_threads>>>(base,num,lambda);	\
    cudaDeviceSynchronize();					\
    cudaError err = cudaGetLastError();				\
    if ( cudaSuccess != err ) {					\
      printf("Cuda error %s\n",cudaGetErrorString( err ));	\
      exit(0);							\
    }								\
  }

#endif


///////////////////////////////////////////////////////
// GPU each thread does one SIMD lane of work
// Host each thread must loop over SIMD lanes
// Use these inline routines to decide what to do according
// to host or device. By using a vector length Nsimd we get Nsimd level
// of read coalescing on the GPU.
// CPU loops over lanes provide compatability and debug with hopefully
// some level of compiler vectorisation in good enough compilers (not holding breath)...
///////////////////////////////////////////////////////
accelerator_inline int get_my_lanes(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return 1;
#else 
  return Nsimd;
#endif
}
accelerator_inline int get_my_lane_offset(int Nsimd) 
{
#ifdef __CUDA_ARCH__
  return ( (threadIdx.x) % Nsimd);
#else
  return 0;
#endif
}
