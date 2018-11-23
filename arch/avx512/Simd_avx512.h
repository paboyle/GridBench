#pragma once

#include <immintrin.h>

  struct Vsplat{
    //Complex float
    inline __m512 operator()(float a, float b){
      return _mm512_set_ps(b,a,b,a,b,a,b,a,b,a,b,a,b,a,b,a);
    }
    // Real float
    inline __m512 operator()(float a){
      return _mm512_set1_ps(a);
    }
    //Complex double
    inline __m512d operator()(double a, double b){
      return _mm512_set_pd(b,a,b,a,b,a,b,a);
    }
    //Real double
    inline __m512d operator()(double a){
      return _mm512_set1_pd(a);
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(__m512 a, float* F){
      _mm512_store_ps(F,a);
    }
    //Double
    inline void operator()(__m512d a, double* D){
      _mm512_store_pd(D,a);
    }
  };


  struct Vstream{
    //Float
    inline void operator()(float * a, __m512 b){
      _mm512_stream_ps(a,b);
      //      _mm512_store_ps(a,b);
    }
    //Double
    inline void operator()(double * a, __m512d b){
      _mm512_stream_pd(a,b);
      //      _mm512_store_pd(a,b);
    }

  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Complex/Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_add_ps(a,b);
    }
    //Complex/Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_add_pd(a,b);
    }
  };

  struct Sub{
    //Complex/Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_sub_ps(a,b);
    }
    //Complex/Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_sub_pd(a,b);
    }
  };

  // Note, we can beat the shuf overhead in chain with two temporaries
  // Ar Ai , Br Bi,  Ai Ar  // one shuf
  //tmpr Ar Br,  Ai Bi    // Mul/Mac/Mac
  //tmpi Br Ai,  Bi Ar    // Mul/Mac/Mac
  // add tmpi,shuf(tmpi)
  // sub tmpr,shuf(tmpi)
  // shuf(tmpr,tmpi).    // Could drop/trade for write mask

  // Gives
  //  2mul,4 mac +add+sub = 8 flop type insns
  //  3shuf + 2 (+shuf)   = 5/6 simd perm and 1/2 the load.

  struct MultComplex{
    // Complex float
    inline __m512 operator()(__m512 a, __m512 b){
      // dup, dup, perm, mul, madd
      __m512 a_real = _mm512_moveldup_ps( a ); // Ar Ar
      __m512 a_imag = _mm512_movehdup_ps( a ); // Ai Ai
      a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( b, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm512_fmaddsub_ps( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
    }
    // Complex double
    inline __m512d operator()(__m512d a, __m512d b){
      __m512d a_real = _mm512_shuffle_pd( a, a, 0x00 );
      __m512d a_imag = _mm512_shuffle_pd( a, a, 0xFF );
      a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( b, 0x55 ) ); 
      return _mm512_fmaddsub_pd( a_real, b, a_imag );
    }
  };
  
  struct Mult{

    inline void mac(__m512 &a, __m512 b, __m512 c){         
       a= _mm512_fmadd_ps( b, c, a);                         
    }
    inline void mac(__m512d &a, __m512d b, __m512d c){
      a= _mm512_fmadd_pd( b, c, a);                   
    }                                             
    // Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_mul_ps(a,b);
    }
    // Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_mul_pd(a,b);
    }
  };

/*
  struct Div{
    // Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_div_ps(a,b);
    }
    // Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_div_pd(a,b);
    }
  };
*/

  struct Conj{
    // Complex single
    inline __m512 operator()(__m512 in){
      return _mm512_mask_sub_ps(in,0xaaaa,_mm512_setzero_ps(),in); // Zero out 0+real 0-imag  
    }
    // Complex double
    inline __m512d operator()(__m512d in){
      return _mm512_mask_sub_pd(in, 0xaa,_mm512_setzero_pd(), in);
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline __m512 operator()(__m512 in, __m512 ret){
      //__m512 tmp = _mm512_mask_sub_ps(in,0xaaaa,_mm512_setzero_ps(),in); // real -imag 
      //return _mm512_shuffle_ps(tmp,tmp,_MM_SELECT_FOUR_FOUR(2,3,1,0));   // 0x4E??
      __m512 tmp = _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
      return _mm512_mask_sub_ps(tmp,0xaaaa,_mm512_setzero_ps(),tmp);
    }
    //Complex double
    inline __m512d operator()(__m512d in, __m512d ret){
      //__m512d tmp = _mm512_mask_sub_pd(in,0xaa,_mm512_setzero_pd(),in); // real -imag 
      //return _mm512_shuffle_pd(tmp,tmp,0x55);
      __m512d tmp = _mm512_shuffle_pd(in,in,0x55);
      return _mm512_mask_sub_pd(tmp,0xaa,_mm512_setzero_pd(),tmp);
    } 
  };

  struct TimesI{
    //Complex single
    inline __m512 operator()(__m512 in, __m512 ret){
      __m512 tmp = _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
      return _mm512_mask_sub_ps(tmp,0x5555,_mm512_setzero_ps(),tmp); 
    }
    //Complex double
    inline __m512d operator()(__m512d in, __m512d ret){
      __m512d tmp = _mm512_shuffle_pd(in,in,0x55);
      return _mm512_mask_sub_pd(tmp,0x55,_mm512_setzero_pd(),tmp); 
    }

  };
  
  // Gpermute utilities consider coalescing into 1 Gpermute
  struct Permute{
    
    static inline __m512 Permute0(__m512 in){
      return _mm512_shuffle_f32x4(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512 Permute1(__m512 in){
      return _mm512_shuffle_f32x4(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };
    static inline __m512 Permute2(__m512 in){
      return _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512 Permute3(__m512 in){
      return _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };

    static inline __m512d Permute0(__m512d in){
      return _mm512_shuffle_f64x2(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512d Permute1(__m512d in){
      return _mm512_shuffle_f64x2(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };
    static inline __m512d Permute2(__m512d in){
      return _mm512_shuffle_pd(in,in,0x55);
    };
    static inline __m512d Permute3(__m512d in){
      return in;
    };

  };

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef __m512i SIMD_Htype;  // Single precision type
  typedef __m512  SIMD_Ftype;  // Single precision type
  typedef __m512d SIMD_Dtype; // Double precision type

  // prefecth
  inline void v_prefetch0(int size, const char *ptr){
    for(int i=0;i<size;i+=64){ //  Define L1 linesize above
      _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
      _mm_prefetch(ptr+i+512,_MM_HINT_T0);
    }
  }
  inline void prefetch_HINT_T0(const char *ptr){
    _mm_prefetch(ptr,_MM_HINT_T0);
  }
  
  // Function name aliases
  typedef Vsplat   VsplatSIMD;
  typedef Vstore   VstoreSIMD;
//  typedef Vset     VsetSIMD;
  typedef Vstream  VstreamSIMD;


  // Arithmetic operations
  typedef Sum         SumSIMD;
  typedef Sub         SubSIMD;
  typedef Mult        MultSIMD;
//  typedef Div         DivSIMD;
  typedef MultComplex MultComplexSIMD;
//  typedef MultRealPart MultRealPartSIMD;
//  typedef MaddRealPart MaddRealPartSIMD;
  typedef Conj        ConjSIMD;
  typedef TimesMinusI TimesMinusISIMD;
  typedef TimesI      TimesISIMD;


