#include <cuda_fp16.h>

#define COALESCE_GRANULARITY ( GEN_SIMD_WIDTH )

template<class pair>
class GpuComplex {
public:
  pair z;
  typedef decltype(z.x) real;
public: 
  GpuComplex() = default;
  accelerator_inline GpuComplex(real re,real im)      { z.x=re; z.y=im; };
  accelerator_inline GpuComplex(const GpuComplex &zz) { z = zz.z;};
  template<class Float> 
  accelerator_inline GpuComplex operator= (const complex<Float> & r){
    Float re = r.real();
    Float im = r.imag();
    z.x = re; 
    z.y = im;
    return *this;
  }
  template<class Float> 
  accelerator_inline operator complex<Float>() const 
  {
    complex<Float> r(z.x,z.y);
    return r;
  }
  friend accelerator_inline  GpuComplex operator+(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x = lhs.z.x + rhs.z.x; 
    r.z.y = lhs.z.y + rhs.z.y; 
    return r; 
  }
  friend accelerator_inline GpuComplex operator-(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x = lhs.z.x - rhs.z.x; 
    r.z.y = lhs.z.y - rhs.z.y; 
    return r; 
  }
  friend accelerator_inline GpuComplex operator*(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x= lhs.z.x*rhs.z.x - lhs.z.y*rhs.z.y; // rr-ii
    r.z.y= lhs.z.x*rhs.z.y + lhs.z.y*rhs.z.x; // ri+ir
    return r;
  }
  friend accelerator_inline GpuComplex real_mult(const GpuComplex &l,const GpuComplex &r) 
  {
    GpuComplex ret;
    ret.z.x = l.z.x*r.z.x;
    ret.z.y = l.z.x*r.z.y;
    return ret;
  }
  friend std::ostream& operator<< (std::ostream& stream, const GpuComplex o){
    stream << "("<< o.z.x << ","<< o.z.y <<")";
    return stream;
  }
};

template<int _N, class _datum>
struct GpuVector {
  _datum v[_N];
  static const int N = _N;
  typedef _datum datum;
};

template<int N,class datum>
accelerator_inline GpuVector<N,datum> operator*(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]*r.v[i];
  }
  return ret;
}
template<int N,class datum>
accelerator_inline GpuVector<N,datum> operator-(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]-r.v[i];
  }
  return ret;
}
template<int N,class datum>
accelerator_inline GpuVector<N,datum> operator+(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]+r.v[i];
  }
  return ret;
}
template<int N,class datum>
accelerator_inline GpuVector<N,datum> operator/(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]/r.v[i];
  }
  return ret;
}

constexpr int NSIMD_RealH    = COALESCE_GRANULARITY / sizeof(half);
constexpr int NSIMD_ComplexH = COALESCE_GRANULARITY / sizeof(half2);
constexpr int NSIMD_RealF    = COALESCE_GRANULARITY / sizeof(float);
constexpr int NSIMD_ComplexF = COALESCE_GRANULARITY / sizeof(float2);
constexpr int NSIMD_RealD    = COALESCE_GRANULARITY / sizeof(double);
constexpr int NSIMD_ComplexD = COALESCE_GRANULARITY / sizeof(double2);
constexpr int NSIMD_Integer  = COALESCE_GRANULARITY / sizeof(Integer);

typedef GpuComplex<half2  > GpuComplexH;
typedef GpuComplex<float2 > GpuComplexF;
typedef GpuComplex<double2> GpuComplexD;

typedef GpuVector<NSIMD_RealH   , half        > GpuVectorRH;
typedef GpuVector<NSIMD_ComplexH, GpuComplexH > GpuVectorCH;
typedef GpuVector<NSIMD_RealF,    float       > GpuVectorRF;
typedef GpuVector<NSIMD_ComplexF, GpuComplexF > GpuVectorCF;
typedef GpuVector<NSIMD_RealD,    double      > GpuVectorRD;
typedef GpuVector<NSIMD_ComplexD, GpuComplexD > GpuVectorCD;
typedef GpuVector<NSIMD_Integer,  Integer     > GpuVectorI;


  struct Vsplat{
    //Complex float
    accelerator_inline GpuVectorCF operator()(float a, float b){
      GpuVectorCF ret;
      for(int i=0;i<GpuVectorCF::N;i++){
	ret.v[i] = typename GpuVectorCF::datum(a,b);
      }
      return ret;
    }
    // Real float
    accelerator_inline GpuVectorRF operator()(float a){
      GpuVectorRF ret;
      for(int i=0;i<GpuVectorRF::N;i++){
	ret.v[i] = typename GpuVectorRF::datum(a);
      }
      return ret;
    }
    //Complex double
    accelerator_inline GpuVectorCD operator()(double a, double b){
      GpuVectorCD ret;
      for(int i=0;i<GpuVectorCD::N;i++){
	ret.v[i] = typename GpuVectorCD::datum(a,b);
      }
      return ret;
    }
    //Real double
    accelerator_inline GpuVectorRD operator()(double a){
      GpuVectorRD ret; 
      for(int i=0;i<GpuVectorRD::N;i++){
	ret.v[i] = typename GpuVectorRD::datum(a);
      }
      return ret;
    }
    //Integer
    accelerator_inline GpuVectorI operator()(Integer a){
      GpuVectorI ret;
      for(int i=0;i<GpuVectorI::N;i++){
	ret.v[i] = typename GpuVectorI::datum(a);
      }
      return ret;
    }
  };

  struct Vstore{
    template<int N,class datum,class P>
    accelerator_inline void operator()(GpuVector<N,datum> a, P* Fp){
      GpuVector<N,datum> *vF = (GpuVector<N,datum> *)Fp;
      *vF = a;
    }
  };

  struct Vstream{
    template<int N,class datum, class P>
    accelerator_inline void operator()(P* F,GpuVector<N,datum> a){
      GpuVector<N,datum> *vF = (GpuVector<N,datum> *)F;
      *vF = a;
    }
  };

  struct Vset{
    // Complex float 
    accelerator_inline GpuVectorCF operator()(ComplexF *a){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i].real(),a[i].imag());
      }
      return ret;
    }
    // Complex double 
    accelerator_inline GpuVectorCD operator()(ComplexD *a){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i].real(),a[i].imag());
      }
      return ret;
    }
    // Real float 
    accelerator_inline GpuVectorRF operator()(float *a){
      typedef GpuVectorRF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(double *a){
      typedef GpuVectorRD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Integer
    accelerator_inline GpuVectorI operator()(Integer *a){
      typedef GpuVectorI vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
      }
      return ret;
    }
  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){
      return a+b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){
      return a+b;
    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a+b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a+b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){
      return a+b;
    }
  };

  struct Sub{
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){
      return a-b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){
      return a-b;
    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a-b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a-b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){
      return a-b;
    }
  };

  struct MultRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]);
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]);
      }
      return ret;
    }
  };

  struct MaddRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b,GpuVectorCF c){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]) +c.v[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b,GpuVectorCD c){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]) +c.v[i];
      }
      return ret;
    }
  };

  struct MultComplex{

    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a*b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a*b;
    }
  };

  struct Mult{
    accelerator_inline void mac(GpuVectorRF &a, GpuVectorRF b, GpuVectorRF c){
      a= a+b*c;
    }
    accelerator_inline void mac(GpuVectorRD &a, GpuVectorRD b, GpuVectorRD c){
      a= a+b*c;
    }
    // Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){
      return a*b;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){
      return a*b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a, GpuVectorI b){
      return a*b;
    }
  };

  struct Div{
    // Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){
      return a/b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){
      return a/b;
    }
    // Danger -- element wise divide fro complex, not complex div. 
    // See Grid_vector_types.h lines around 735, applied after "toReal"
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a, GpuVectorCF b){
      GpuVectorCF ret;
      for(int i=0;i< GpuVectorCF::N;i++){
	ret.v[i].z.x = a.v[i].z.x / b.v[i].z.x;
	ret.v[i].z.y = a.v[i].z.y / b.v[i].z.y;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a, GpuVectorCD b){
      GpuVectorCD ret;
      for(int i=0;i< GpuVectorCD::N;i++){
	ret.v[i].z.x = a.v[i].z.x / b.v[i].z.x;
	ret.v[i].z.y = a.v[i].z.y / b.v[i].z.y;
      }
      return ret;
    }
  };


  struct Conj{
    // Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.x;
	ret.v[i].z.y =-in.v[i].z.y;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.x;
	ret.v[i].z.y =-in.v[i].z.y;
      }
      return ret;
    }
  };

  struct TimesMinusI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.y;
	ret.v[i].z.y =-in.v[i].z.x;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.y;
	ret.v[i].z.y =-in.v[i].z.x;
      }
      return ret;
    }
  };

  struct TimesI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x =-in.v[i].z.y;
	ret.v[i].z.y = in.v[i].z.x;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x =-in.v[i].z.y;
	ret.v[i].z.y = in.v[i].z.x;
      }
      return ret;
    }
  };

  struct Permute{

    template <int n,typename vec>				       
    static accelerator_inline vec PermuteN(vec in) {   
      vec out;					
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out.v[i] = in.v[i^_mask];
      }
      return out;	
    }
    
    template <typename vec>  static accelerator_inline vec Permute0(vec in) { return PermuteN<0,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute1(vec in) { return PermuteN<1,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute2(vec in) { return PermuteN<2,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute3(vec in) { return PermuteN<3,vec>(in);  }
    
  };


//////////////////////////////////////////////
// Some Template specialization

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
//////////////////////////////////////////////////////////////////////////////////////
  typedef GpuVectorRH  SIMD_Htype; // Single precision type
  typedef GpuVectorRF  SIMD_Ftype; // Single precision type
  typedef GpuVectorRD  SIMD_Dtype; // Double precision type
  typedef GpuVectorI   SIMD_Itype; // Integer type

  typedef GpuVectorCH  SIMD_CHtype; // Single precision type
  typedef GpuVectorCF  SIMD_CFtype; // Single precision type
  typedef GpuVectorCD  SIMD_CDtype; // Double precision type

  // Function name aliases
  typedef Vsplat   VsplatSIMD;
  typedef Vstore   VstoreSIMD;
  typedef Vset     VsetSIMD;
  typedef Vstream  VstreamSIMD;

  // Arithmetic operations
  typedef Sum         SumSIMD;
  typedef Sub         SubSIMD;
  typedef Div         DivSIMD;
  typedef Mult        MultSIMD;
  typedef MultComplex MultComplexSIMD;
  typedef MultRealPart MultRealPartSIMD;
  typedef MaddRealPart MaddRealPartSIMD;
  typedef Conj        ConjSIMD;
  typedef TimesMinusI TimesMinusISIMD;
  typedef TimesI      TimesISIMD;
