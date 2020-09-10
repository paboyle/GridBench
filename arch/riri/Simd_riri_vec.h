
template<class datum> struct wordsize {  };
template<> struct wordsize<float>  { static const int bytes = sizeof(float) ; typedef float word_type; };
template<> struct wordsize<double> { static const int bytes = sizeof(double); typedef double word_type; };

//typedef half      vhalf   __attribute__ ((vector_size (GEN_SIMD_WIDTH)));

template<class _datum> struct datum2 {
  typedef _datum datum;
  datum x;
  datum y;
};

//typedef datum2<half>     half2;
typedef datum2<float>   float2;
typedef datum2<double> double2;

template<class pair>
class GpuComplex {
public:
  pair z;
  typedef decltype(z.x) real;
  typedef real value_type;
public: 
  GpuComplex() = default;
  GpuComplex(real x,real y) { z.x=x; z.y=y; };
  static const int N = sizeof(real)/(wordsize<real>::bytes);
  typedef typename wordsize<real>::word_type word_type;
  typedef GpuComplex<datum2<word_type> > scalar_type;
  accelerator_inline GpuComplex(const GpuComplex &zz) { z = zz.z;};
    
  // *=,+=,-= operators
  accelerator_inline GpuComplex &operator+=(const GpuComplex &r) {
    *this = (*this) + r;
    return *this;
  }
  accelerator_inline GpuComplex &operator-=(const GpuComplex &r) {
    *this = (*this) - r;
    return *this;
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
  friend accelerator_inline GpuComplex conj(const GpuComplex &l) 
  {
    GpuComplex ret;
    ret.z.x = l.z.x;
    ret.z.y =-l.z.y;
    return ret;
  }
  friend std::ostream& operator<< (std::ostream& stream, const GpuComplex o){
    for(int i=0;i<N;i++) {
      //      stream << i<< " ("<< o.z.x[i] << ","<< o.z.y[i] <<")";
    }
    return stream;
  }
};

template<int _N, class _datum>
struct GpuVector {
  _datum v[_N];
  static const int N = _N;
  typedef _datum datum;
  typedef typename _datum::value_type word_type;
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

//typedef std::complex<float> GpuComplexF;
//typedef std::complex<double> GpuComplexD;

typedef GpuComplex<float2> GpuComplexF;
typedef GpuComplex<double2> GpuComplexD;

constexpr int NSIMD_ComplexF = EXPAND_SIMD;
constexpr int NSIMD_ComplexD = EXPAND_SIMD;

typedef GpuVector<NSIMD_ComplexF, GpuComplexF > GpuVectorCF;
typedef GpuVector<NSIMD_ComplexD, GpuComplexD > GpuVectorCD;


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
    /*
    accelerator_inline GpuVectorRF operator()(float a){
      GpuVectorRF ret;
      for(int i=0;i<GpuVectorRF::N;i++){
	ret.v[i] = typename GpuVectorRF::datum(a);
      }
      return ret;
    }
    */
    //Complex double
    accelerator_inline GpuVectorCD operator()(double a, double b){
      GpuVectorCD ret;
      for(int i=0;i<GpuVectorCD::N;i++){
	ret.v[i] = typename GpuVectorCD::datum(a,b);
      }
      return ret;
    }
    //Real double
    /*
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
    */
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
    /*
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
    */
  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Real float
    //    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){      return a+b;    }
    //    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){      return a+b;    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){      return a+b;    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){      return a+b;    }
    //    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){      return a+b;    }
  };

  struct Sub{
    //    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){      return a-b;    }
    //    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){      return a-b;    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){      return a-b;    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){      return a-b;    }
    //    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){     return a-b;    }
  };

  struct MultRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      typedef GpuVectorCF vec;
      vec ret;
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      typedef GpuVectorCD vec;
      vec ret;
      return ret;
    }
  };

  struct MaddRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b,GpuVectorCF c){
      typedef GpuVectorCF vec;
      vec ret;
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b,GpuVectorCD c){
      typedef GpuVectorCD vec;
      vec ret;
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
    //    accelerator_inline void mac(GpuVectorRF &a, GpuVectorRF b, GpuVectorRF c){     a= a+b*c;    }
    //    accelerator_inline void mac(GpuVectorRD &a, GpuVectorRD b, GpuVectorRD c){      a= a+b*c;    }
    // Real float
    //    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){      return a*b;    }
    // Real double
    //    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){      return a*b;    }
    //    accelerator_inline GpuVectorI operator()(GpuVectorI a, GpuVectorI b){      return a*b;    }
  };

  struct Div{
    // Real float
    //    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){     return a/b;    }
    //    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){      return a/b;    }
    // Danger -- element wise divide fro complex, not complex div. 
    // See Grid_vector_types.h lines around 735, applied after "toReal"
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a, GpuVectorCF b){
      GpuVectorCF ret;
      for(int i=0;i< GpuVectorCF::N;i++){
	//	ret.v[i] = a.v[i] / b.v[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a, GpuVectorCD b){
      GpuVectorCD ret;
      for(int i=0;i< GpuVectorCD::N;i++){
	//	ret.v[i] = a.v[i] / b.v[i];
      }
      return ret;
    }
  };


  struct Conj
  {
    // Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = conj(in.v[i]);
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = conj(in.v[i]);
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
	ret.v[i] = in.v[i]*GpuComplexF(0.,-1.);
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = in.v[i]*GpuComplexD(0.,-1.);
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
	ret.v[i] = in.v[i]*GpuComplexF(0.,1.);
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = in.v[i]*GpuComplexD(0.,1.);
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
//  typedef GpuVectorRF  SIMD_Ftype; // Single precision type
//  typedef GpuVectorRD  SIMD_Dtype; // Double precision type
//  typedef GpuVectorI   SIMD_Itype; // Integer type

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

#if defined(GRID_SYCL_SIMT) 


/*Small support to allow GPU coalesced access*/

template<class vec>
typename vec::vector_type::datum coalescedRead(const vec &in, int lane){
  return in.v.v[lane];
}
template<class vec>
void coalescedWrite(vec &out,const typename vec::vector_type::datum &in,int lane){
  out.v.v[lane] = in;
}
template<int ptype,class vec>
typename vec::vector_type::datum coalescedReadPermute(const vec & __restrict__ in,int doperm,int lane)
{
  typename vec::vector_type::datum ret;
  constexpr int mask = DATA_SIMD >> (ptype + 1); // Keep the permutes fixed as SIMD expanded
  int plane= doperm ? lane ^ mask : lane;
  return in.v.v[plane];
}
#endif
