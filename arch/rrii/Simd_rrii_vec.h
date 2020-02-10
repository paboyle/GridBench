
#define COALESCE_GRANULARITY ( GEN_SIMD_WIDTH )

//typedef uint16_t half;
#include <CL/sycl.hpp>
//#include <CL/sycl/group.hpp>
//#include <CL/__spirv/spirv_vars.hpp>

// SYCL
#ifdef __SYCL_DEVICE_ONLY__
inline int SIMTlane(void) { return __spirv_BuiltInGlobalInvocationId.x ; };
#else
inline int SIMTlane(void) { return 0; };
#endif

typedef float     vfloat  __attribute__ ((vector_size (GEN_SIMD_WIDTH)));
typedef double    vdouble __attribute__ ((vector_size (GEN_SIMD_WIDTH)));
typedef Integer vinteger  __attribute__ ((vector_size (GEN_SIMD_WIDTH)));

template<class datum> struct wordsize {  };
template<> struct wordsize<float>  { static const int bytes = sizeof(float) ; typedef float word_type; };
template<> struct wordsize<double> { static const int bytes = sizeof(double); typedef double word_type; };

template<> struct wordsize<vfloat>  { static const int bytes = sizeof(float)   ; typedef float word_type; };
template<> struct wordsize<vdouble> { static const int bytes = sizeof(double)  ; typedef double word_type; };
template<> struct wordsize<vinteger>{ static const int bytes = sizeof(Integer) ; typedef Integer word_type; };

//typedef half      vhalf   __attribute__ ((vector_size (GEN_SIMD_WIDTH)));

template<class _datum> struct datum2 {
  typedef _datum datum;
  datum x;
  datum y;
};

//typedef datum2<half>     half2;
typedef datum2<float>   float2;
typedef datum2<double> double2;

//typedef datum2<vhalf>     vhalf2;
typedef datum2<vfloat>   vfloat2;
typedef datum2<vdouble> vdouble2;

template<> struct wordsize<vfloat2>  { static const int bytes = sizeof(float)   ; typedef float word_type; };
template<> struct wordsize<vdouble2> { static const int bytes = sizeof(double)  ; typedef double word_type; };

template<class pair>
class GpuComplex {
public:
  pair z;
  typedef decltype(z.x) real;
public: 
  GpuComplex() = default;
  static const int N = sizeof(real)/(wordsize<real>::bytes);
  typedef typename wordsize<real>::word_type word_type;
  accelerator_inline GpuComplex(const GpuComplex &zz) { z = zz.z;};
    
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
    for(int i=0;i<N;i++) {
      //      stream << i<< " ("<< o.z.x[i] << ","<< o.z.y[i] <<")";
    }
    return stream;
  }
};

template<class scalar>
class GpuReal {
public:
  typedef scalar real;
  real z;
  static const int N = sizeof(scalar)/wordsize<scalar>::bytes;
public: 
  GpuReal() = default;
  typedef typename wordsize<scalar>::word_type word_type;
  accelerator_inline GpuReal(const GpuReal &zz) { z = zz.z;};
  friend accelerator_inline  GpuReal operator+(const GpuReal &lhs,const GpuReal &rhs) { 
    GpuReal r ; 
    r.z = lhs.z + rhs.z ;
    return r; 
  }
  friend accelerator_inline GpuReal operator-(const GpuReal &lhs,const GpuReal &rhs) { 
    GpuReal r ; 
    r.z = lhs.z - rhs.z; 
    return r; 
  }
  friend accelerator_inline GpuReal operator*(const GpuReal &lhs,const GpuReal &rhs) { 
    GpuReal r ; 
    r.z= lhs.z*rhs.z ;
    return r;
  }
};

//typedef GpuComplex<half2  > GpuComplexH;
typedef GpuComplex<float2 > GpuComplexF;
typedef GpuComplex<double2> GpuComplexD;

//typedef GpuComplex<vhalf2  > GpuVectorCH;
typedef GpuComplex<vfloat2 > GpuVectorCF;
typedef GpuComplex<vdouble2> GpuVectorCD;

//typedef GpuReal<vhalf>    GpuVectorRH;
typedef GpuReal<vfloat>   GpuVectorRF;
typedef GpuReal<vdouble>  GpuVectorRD;
typedef GpuReal<vinteger> GpuVectorI;


struct Vsplat{

  accelerator_inline GpuVectorCF operator()(float a, float b){
      GpuVectorCF ret;
      float * x_p = (float *)&ret.z.x;
      float * y_p = (float *)&ret.z.y;
      for(int i=0;i<GpuVectorCF::N;i++){
	x_p[i] = a;
	y_p[i] = b;
      }
      return ret;
    }
    accelerator_inline GpuVectorRF operator()(float a){
      GpuVectorRF ret;
      float *z_p= (float *)&ret.z;
      for(int i=0;i<GpuVectorRF::N;i++){
	z_p[i] = a;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(double a, double b){
      GpuVectorCD ret;
      double *x_p=(double *)&ret.z.x;
      double *y_p=(double *)&ret.z.y;
      for(int i=0;i<GpuVectorCD::N;i++){
	x_p[i] = a;
	y_p[i] = b;
      }
      return ret;
    }
    accelerator_inline GpuVectorRD operator()(double a){
      GpuVectorRD ret;
      double *z_p = (double *)&ret.z;
      for(int i=0;i<GpuVectorRD::N;i++){
	z_p[i] = a;
      }
      return ret;
    }
    //Integer
    accelerator_inline GpuVectorI operator()(Integer a){
      GpuVectorI ret;
      Integer *I_p=(Integer *)&ret.z;
      for(int i=0;i<GpuVectorI::N;i++){
	I_p[i] = a;
      }
      return ret;
    }
  };

  struct Vstore{
    template<class datum,class P>
    accelerator_inline void operator()(GpuComplex<datum> a, P* Fp){
      GpuComplex<datum> *vF = (GpuComplex<datum> *)Fp;
      *vF = a;
    }
    template<class datum,class P>
    accelerator_inline void operator()(GpuReal<datum> a, P* Fp){
      GpuReal<datum> *vF = (GpuReal<datum> *)Fp;
      *vF = a;
    }
  };

  struct Vstream{
    template<class datum,class P>
      accelerator_inline void operator()( P* Fp,GpuComplex<datum> a){
      GpuComplex<datum> *vF = (GpuComplex<datum> *)Fp;
      *vF = a;
    }
    template<class datum,class P>
      accelerator_inline void operator()( P* Fp, GpuReal<datum> a){
      GpuReal<datum> *vF = (GpuReal<datum> *)Fp;
      *vF = a;
    }
  };

  struct Vset{
    // Complex float 
    accelerator_inline GpuVectorCF operator()(ComplexF *a){
      typedef GpuVectorCF vec;
      vec ret;
      float *x_p=(float *)&ret.z.x;
      float *y_p=(float *)&ret.z.y;
      for(int i=0;i<vec::N;i++){
	x_p[i] = a[i].real();
	y_p[i] = a[i].imag();
      }
      return ret;
    }
    // Complex double 
    accelerator_inline GpuVectorCD operator()(ComplexD *a){
      typedef GpuVectorCD vec;
      vec ret;
      double *x_p=(double *)&ret.z.x;
      double *y_p=(double *)&ret.z.y;
      for(int i=0;i<vec::N;i++){
	x_p[i] = a[i].real();
	y_p[i] = a[i].imag();
      }
      return ret;
    }
    // Real float 
    accelerator_inline GpuVectorRF operator()(float *a){
      typedef GpuVectorRF vec;
      vec ret;
      float *z_p = (float *)&ret.z;
      for(int i=0;i<vec::N;i++){
	z_p[i] = a[i];
      }
      return ret;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(double *a){
      typedef GpuVectorRD vec;
      vec ret;
      double *z_p = (double *)&ret.z;
      for(int i=0;i<vec::N;i++){
	z_p[i] = a[i];
      }
      return ret;
    }
    // Integer
    accelerator_inline GpuVectorI operator()(Integer *a){
      typedef GpuVectorI vec;
      vec ret;
      Integer *z_p=(Integer *)&ret.z;
      for(int i=0;i<vec::N;i++){
	z_p[i] = a[i];
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
      vec ret=real_mult(a,b);
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      typedef GpuVectorCD vec;
      vec ret=real_mult(a,b);
      return ret;
    }
  };

  struct MaddRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b,GpuVectorCF c){
      typedef GpuVectorCF vec;
      vec ret;
      ret = real_mult(a,b) +c;
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b,GpuVectorCD c){
      typedef GpuVectorCD vec;
      vec ret;
      ret = real_mult(a,b) +c;
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
      GpuVectorRF ret;
      ret.z=a.z/b.z;
      return ret;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){
      GpuVectorRD ret;
      ret.z=a.z/b.z;
      return ret;
    }
    // Danger -- element wise divide fro complex, not complex div. 
    // See Grid_vector_types.h lines around 735, applied after "toReal"
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a, GpuVectorCF b){
      GpuVectorCF ret;
      ret.z.x = a.z.x / b.z.x;
      ret.z.y = a.z.y / b.z.y;
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a, GpuVectorCD b){
      GpuVectorCD ret;
      ret.z.x = a.z.x / b.z.x;
      ret.z.y = a.z.y / b.z.y;
      return ret;
    }
  };

  struct Conj{
    // Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      ret.z.x = in.z.x;
      ret.z.y = in.z.y*(-1.0);
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      ret.z.x = in.z.x;
      ret.z.y = in.z.y*(-1.0);
      return ret;
    }
  };

  struct TimesMinusI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      ret.z.x = in.z.y;
      ret.z.y = in.z.x*(-1.0);
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      ret.z.x = in.z.y;
      ret.z.y = in.z.x*(-1.0);
      return ret;
    }
  };

  struct TimesI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      ret.z.x = in.z.y*(-1.0);
      ret.z.y = in.z.x;
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      ret.z.x = in.z.y*(-1.0);
      ret.z.y = in.z.x;
      return ret;
    }
  };

  struct Permute{

    template <int n,class pair>				       
    static accelerator_inline GpuComplex<pair> PermuteN(GpuComplex<pair> in) {
      typedef GpuComplex<pair> vec;
      typedef typename wordsize<pair>::word_type word_type;
      vec out;
      word_type *out_x_p = (word_type *)&out.z.x;
      word_type *out_y_p = (word_type *)&out.z.y;
      word_type * in_x_p = (word_type *)&in.z.x;
      word_type * in_y_p = (word_type *)&in.z.y;
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out_x_p[i] = in_x_p[i^_mask];
	out_y_p[i] = in_y_p[i^_mask];
      }
      return out;	
    }
    template <int n,class datum>				       
    static accelerator_inline GpuReal<datum> PermuteN(GpuReal<datum> in) {
      typedef GpuReal<datum> vec;
      typedef typename GpuReal<datum>::word_type word_type;
      vec out;
      word_type *out_p = (word_type *)&out.z;
      word_type * in_p = (word_type *)& in.z;
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out_p[i] = in_p[i^_mask];
      }
      return out;	
    }
    
    template <typename vec>  static accelerator_inline vec Permute0(vec in) { return PermuteN<0>(in);  }
    template <typename vec>  static accelerator_inline vec Permute1(vec in) { return PermuteN<1>(in);  }
    template <typename vec>  static accelerator_inline vec Permute2(vec in) { return PermuteN<2>(in);  }
    template <typename vec>  static accelerator_inline vec Permute3(vec in) { return PermuteN<3>(in);  }
    
  };


//////////////////////////////////////////////
// Some Template specialization

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
//////////////////////////////////////////////////////////////////////////////////////
//  typedef GpuVectorRH  SIMD_Htype; // Single precision type
  typedef GpuVectorRF  SIMD_Ftype; // Single precision type
  typedef GpuVectorRD  SIMD_Dtype; // Double precision type
  typedef GpuVectorI   SIMD_Itype; // Integer type

//  typedef GpuVectorCH  SIMD_CHtype; // Single precision type
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
