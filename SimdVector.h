#pragma once

#ifdef GEN
#include "arch/generic/Simd_generic.h"
#endif
#ifdef SSE4
#include "arch/sse/Simd_sse4.h"
#endif
#if defined(AVX1) || defined (AVXFMA) || defined(AVX2) || defined(AVXFMA4)
#include "arch/avx/Simd_avx.h"
#endif
#if defined AVX512
#include "arch/avx512/Simd_avx512.h"
#endif

//////////////////////////////////////
// demote a vector to real type
//////////////////////////////////////
#include <type_traits>
// type alias used to simplify the syntax of std::enable_if
template <typename T> using Invoke = typename T::type;
template <typename Condition, typename ReturnType> using EnableIf    = Invoke<std::enable_if<Condition::value, ReturnType> >;
template <typename Condition, typename ReturnType> using NotEnableIf = Invoke<std::enable_if<!Condition::value, ReturnType> >;

////////////////////////////////////////////////////////
// Check for complexity with type traits
template <typename T> struct is_complex : public std::false_type {};
template <> struct is_complex<std::complex<double> > : public std::true_type {};
template <> struct is_complex<std::complex<float> > : public std::true_type {};

template <typename T>              using IfReal    = Invoke<std::enable_if<std::is_floating_point<T>::value, int> >;
template <typename T>              using IfComplex = Invoke<std::enable_if<is_complex<T>::value, int> >;
template <typename T>              using IfNotReal    = Invoke<std::enable_if<!std::is_floating_point<T>::value, int> >;
template <typename T>              using IfNotComplex = Invoke<std::enable_if<!is_complex<T>::value, int> >;

////////////////////////////////////////////////////////
// Define the operation templates functors

template <class Out, class Input1, class Input2, class Operation>
Out binary(Input1 src_1, Input2 src_2, Operation op) {
  return op(src_1, src_2);
}

template <class Out, class Input, class Operation>
Out unary(Input src, Operation op) {
  return op(src);
}
///////////////////////////////////////////////

/*
  @brief Grid_simd class for the SIMD vector type operations
 */
template <class Scalar_type, class Vector_type>
class Simd {
 public:
  Vector_type v; // Only data payload

 public:
  typedef Vector_type vector_type;
  typedef Scalar_type scalar_type;

  static inline constexpr int Nsimd(void) {
    return sizeof(Vector_type) / sizeof(Scalar_type);
  }

  Simd &operator=(const Simd &&rhs) {
    v = rhs.v;
    return *this;
  };
  Simd &operator=(const Simd &rhs) {
    v = rhs.v;
    return *this;
  };  // faster than not declaring it and leaving to the compiler
  Simd() = default;
  Simd(const Simd &rhs) : v(rhs.v){};  // compiles in movaps
  Simd(const Simd &&rhs) : v(rhs.v){};

  /////////////////////////////
  // Constructors
  /////////////////////////////
  Simd &operator=(Zero &z) {
    vzero(*this);
    return (*this);
  }

  // Enable if complex type
  template <typename S = Scalar_type>
  Simd(const typename std::enable_if<is_complex<S>::value, S>::type a) {
    vsplat(*this, a);
  };

  Simd(const RealD a) { vsplat(*this, Scalar_type(a)); };

  ///////////////////////////////////////////////
  // mac, mult, sub, add, adj
  ///////////////////////////////////////////////

  // FIXME -- alias this to an inline MAC struct.
  friend inline void mac(Simd *__restrict__ y,
                         const Simd *__restrict__ a,
                         const Simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };

  friend inline void mult(Simd *__restrict__ y,
                          const Simd *__restrict__ l,
                          const Simd *__restrict__ r) {
    *y = (*l) * (*r);
  }

  friend inline void sub(Simd *__restrict__ y,
                         const Simd *__restrict__ l,
                         const Simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Simd *__restrict__ y,
                         const Simd *__restrict__ l,
                         const Simd *__restrict__ r) {
    *y = (*l) + (*r);
  }
  friend inline void mac(Simd *__restrict__ y,
                         const Scalar_type *__restrict__ a,
                         const Simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend inline void mult(Simd *__restrict__ y,
                          const Scalar_type *__restrict__ l,
                          const Simd *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend inline void sub(Simd *__restrict__ y,
                         const Scalar_type *__restrict__ l,
                         const Simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Simd *__restrict__ y,
                         const Scalar_type *__restrict__ l,
                         const Simd *__restrict__ r) {
    *y = (*l) + (*r);
  }

  friend inline void mac(Simd *__restrict__ y,
                         const Simd *__restrict__ a,
                         const Scalar_type *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend inline void mult(Simd *__restrict__ y,
                          const Simd *__restrict__ l,
                          const Scalar_type *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend inline void sub(Simd *__restrict__ y,
                         const Simd *__restrict__ l,
                         const Scalar_type *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Simd *__restrict__ y,
                         const Simd *__restrict__ l,
                         const Scalar_type *__restrict__ r) {
    *y = (*l) + (*r);
  }


  ///////////////////////
  // Vstore
  ///////////////////////
  friend inline void vstore(const Simd &ret, Scalar_type *a) {
    binary<void>(ret.v, (RealD *)a, VstoreSIMD());
  }

  ////////////////////////////
  // operator scalar * simd
  ////////////////////////////
  friend inline Simd operator*(const Scalar_type &a, Simd b) {
    Simd va;
    vsplat(va, a);
    return va * b;
  }
  friend inline Simd operator*(Simd b, const Scalar_type &a) {
    return a * b;
  }

  ///////////////////////
  // Unary negation
  ///////////////////////
  friend inline Simd operator-(const Simd &r) {
    Simd ret;
    vzero(ret);
    ret = ret - r;
    return ret;
  }
  // *=,+=,-= operators
  inline Simd &operator*=(const Simd &r) {
    *this = (*this) * r;
    return *this;
    // return (*this)*r; ?
  }
  inline Simd &operator+=(const Simd &r) {
    *this = *this + r;
    return *this;
  }
  inline Simd &operator-=(const Simd &r) {
    *this = *this - r;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////
  // General permute; assumes vector length is same across
  // all subtypes; may not be a good assumption, but could
  // add the vector width as a template param for BG/Q for example
  ////////////////////////////////////////////////////////////////////
  friend inline void permute0(Simd &y, Simd b) {
    y.v = Permute::Permute0(b.v);
  }
  friend inline void permute1(Simd &y, Simd b) {
    y.v = Permute::Permute1(b.v);
  }
  friend inline void permute2(Simd &y, Simd b) {
    y.v = Permute::Permute2(b.v);
  }
  friend inline void permute3(Simd &y, Simd b) {
    y.v = Permute::Permute3(b.v);
  }
  friend inline void permute(Simd &y, Simd b, int perm) {

         if(perm==3) permute3(y, b);
    else if(perm==2) permute2(y, b);
    else if(perm==1) permute1(y, b);
    else if(perm==0) permute0(y, b);
  }
  
};  // end of Simd class definition

inline void permute(ComplexD &y,ComplexD b, int perm) {  y=b; }
inline void permute(ComplexF &y,ComplexF b, int perm) {  y=b; }
inline void permute(RealD &y,RealD b, int perm) {  y=b; }
inline void permute(RealF &y,RealF b, int perm) {  y=b; }


///////////////////////
// Splat
///////////////////////

// this is only for the complex version pass real and im parts
template <class S, class V, IfComplex<S> = 0, class ABtype>
inline void vsplat(Simd<S, V> &ret, ABtype a, ABtype b) {
  ret.v = binary<V>(a, b, VsplatSIMD());
}

// overload if complex
template <class S, class V>
inline void vsplat(Simd<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), imag(c));
}


///////////////////////////////////////////////
// Initialise to 1,0,i for the correct types
///////////////////////////////////////////////
// For complex types
template <class S, class V, IfComplex<S> = 0>
inline void vone(Simd<S, V> &ret) {
  vsplat(ret, S(1.0, 0.0));
}
template <class S, class V, IfComplex<S> = 0>
inline void vzero(Simd<S, V> &ret) {
  vsplat(ret, S(0.0, 0.0));
}  // use xor?
template <class S, class V, IfComplex<S> = 0>
inline void vcomplex_i(Simd<S, V> &ret) {
  vsplat(ret, S(0.0, 1.0));
}
// if not complex overload here
template <class S, class V, IfReal<S> = 0>
inline void vone(Simd<S, V> &ret) {
  vsplat(ret, S(1.0));
}
template <class S, class V, IfReal<S> = 0>
inline void vzero(Simd<S, V> &ret) {
  vsplat(ret, S(0.0));
}

// For integral types
template <class S, class V>
inline void zeroit(Simd<S, V> &z) {
  vzero(z);
}

///////////////////////
// Vstream
///////////////////////
template <class S, class V, IfReal<S> = 0>
inline void vstream(Simd<S, V> &out, const Simd<S, V> &in) {
  binary<void>((S *)&out.v, in.v, VstreamSIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline void vstream(Simd<S, V> &out, const Simd<S, V> &in) {
  typedef typename S::value_type T;
  binary<void>((T *)&out.v, in.v, VstreamSIMD());
}

////////////////////////////////////
// Arithmetic operator overloads +,-,*
////////////////////////////////////
template <class S, class V>
inline Simd<S, V> operator+(Simd<S, V> a, Simd<S, V> b) {
  Simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SumSIMD());
  return ret;
};

template <class S, class V>
inline Simd<S, V> operator-(Simd<S, V> a, Simd<S, V> b) {
  Simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SubSIMD());
  return ret;
};

// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
inline Simd<S, V> operator*(Simd<S, V> a, Simd<S, V> b) {
  Simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultComplexSIMD());
  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
inline Simd<S, V> operator*(Simd<S, V> a, Simd<S, V> b) {
  Simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultSIMD());
  return ret;
};


///////////////////////
// Conjugate
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline Simd<S, V> conjugate(const Simd<S, V> &in) {
  Simd<S, V> ret;
  ret.v = unary<V>(in.v, ConjSIMD());
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Simd<S, V> conjugate(const Simd<S, V> &in) {
  return in;  // for real objects
}

///////////////////////
// timesMinusI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline void timesMinusI(Simd<S, V> &ret, const Simd<S, V> &in) {
  ret.v = binary<V>(in.v, ret.v, TimesMinusISIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline Simd<S, V> timesMinusI(const Simd<S, V> &in) {
  Simd<S, V> ret;
  timesMinusI(ret, in);
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Simd<S, V> timesMinusI(const Simd<S, V> &in) {
  return in;
}

///////////////////////
// timesI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline void timesI(Simd<S, V> &ret, const Simd<S, V> &in) {
  ret.v = binary<V>(in.v, ret.v, TimesISIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline Simd<S, V> timesI(const Simd<S, V> &in) {
  Simd<S, V> ret;
  timesI(ret, in);
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Simd<S, V> timesI(const Simd<S, V> &in) {
  return in;
}

///////////////////////////////
// Define available types
///////////////////////////////
typedef Simd<float , SIMD_Ftype> vRealF;
typedef Simd<double, SIMD_Dtype> vRealD;
typedef Simd<std::complex<float>, SIMD_Ftype> vComplexF;
typedef Simd<std::complex<double>, SIMD_Dtype> vComplexD;

static_assert(sizeof(SIMD_Ftype) == sizeof(SIMD_Dtype), "SIMD vector lengths incorrect");

/////////////////////////////////////////
// Some traits to recognise the types
/////////////////////////////////////////
template <typename T>
struct is_simd : public std::false_type {};
template <> struct is_simd<vRealF>     : public std::true_type {};
template <> struct is_simd<vRealD>     : public std::true_type {};
template <> struct is_simd<vComplexF>  : public std::true_type {};
template <> struct is_simd<vComplexD>  : public std::true_type {};

template <typename T> using IfSimd    = Invoke<std::enable_if<is_simd<T>::value, int> >;
template <typename T> using IfNotSimd = Invoke<std::enable_if<!is_simd<T>::value, unsigned> >;

