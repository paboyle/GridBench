/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_simd.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include "Simd.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////
// Extract/merge a fundamental vector type, to pointer array with offset
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Extract a fundamental vector type to scalar array 
////////////////////////////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void extract(const vsimd &y,std::vector<scalar> &extracted){

  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar *buf = (scalar *)&y;
  for(int i=0;i<Nextr;i++){
    extracted[i]=buf[i*s];
  }
};

////////////////////////////////////////////////////////////////////////
// Merge simd vector from array of scalars
////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void merge(vsimd &y,std::vector<scalar> &extracted){
  int Nextr=extracted.size();
  static const int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;
  scalar *buf = (scalar *)&y;

  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i]; // replicates value
    }
  }
};

/*
////////////////////////////////////////////////////////////////////////
// Extract to contiguous array scalar object
////////////////////////////////////////////////////////////////////////
template<class vobj> inline void extract(const vobj &vec,std::vector<typename vobj::scalar_object> &extracted)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;

  static const int Nsimd=sizeof(vector_type)/sizeof(scalar_type);
  static const int words=sizeof(vobj)/sizeof(vector_type);
  int Nextr=extracted.size();
  int s=Nsimd/Nextr;

  std::vector<scalar_type *> pointers(Nextr);
  for(int i=0;i<Nextr;i++) 
    pointers[i] =(scalar_type *)& extracted[i];

  vector_type *vp = (vector_type *)&vec;
  for(int w=0;w<words;w++){
    extract<vector_type,scalar_type>(&vp[w],pointers,w);
  }
}


template<class vsimd,class scalar>
inline void extract(vsimd * y,  std::vector<scalar *> &extracted,int offset)
{
  // FIXME: bounce off memory is painful
  static const int Nsimd=sizeof(vsimd)/sizeof(scalar);
  int Nextr=extracted.size();
  int s=Nsimd/Nextr;

  scalar*buf = (scalar *)y;
  for(int i=0;i<Nextr;i++){
    extracted[i][offset] = buf[i*s];
  }
};
////////////////////////////////////////////////////////////////////////
// Merge simd vector from array of scalars to pointer array with offset
////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void merge(vsimd  * y, std::vector<scalar *> &extracted,int offset)
{
  static const int Nsimd=sizeof(vsimd)/sizeof(scalar);

  int Nextr=extracted.size();
  int s=Nsimd/Nextr; // can have sparse occupation of simd vector if simd_layout does not fill it
                     // replicate n-fold. Use to allow Integer masks to 
                     // predicate floating point of various width assignments and maintain conformable.
  scalar *buf =(scalar *) y;
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i][offset];
    }
  }
};
*/

using namespace std;

class funcPlus {
public:
  funcPlus() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1+i2;}
  std::string name(void) const { return std::string("Plus"); }
};
class funcMinus {
public:
  funcMinus() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1-i2;}
  std::string name(void) const { return std::string("Minus"); }
};
class funcTimes {
public:
  funcTimes() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = i1*i2;}
  std::string name(void) const { return std::string("Times"); }
};
class funcConj {
public:
  funcConj() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = conjugate(i1);}
  std::string name(void) const { return std::string("Conj"); }
};
class funcTimesI {
public:
  funcTimesI() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = timesI(i1);}
  std::string name(void) const { return std::string("timesI"); }
};
class funcTimesMinusI {
public:
  funcTimesMinusI() {};
  template<class vec> void operator()(vec &rr,vec &i1,vec &i2) const { rr = timesMinusI(i1);}
  std::string name(void) const { return std::string("timesMinusI"); }
};


template<class scal, class vec,class functor > 
void Tester(const functor &func)
{
  
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  std::vector<scal> result(Nsimd);
  std::vector<scal> reference(Nsimd);

  std::vector<vec,alignedAllocator<vec> > buf(3);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];
  vec & v_result = buf[2];


  for(int i=0;i<Nsimd;i++){
    input1[i] = drand48();
    input2[i] = drand48();
    result[i] = drand48();
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);
  merge<vec,scal>(v_result,result);

  func(v_result,v_input1,v_input2);

  for(int i=0;i<Nsimd;i++) {
    func(reference[i],input1[i],input2[i]);
  }

  extract<vec,scal>(v_result,result);

  std::cout << " " << func.name() << std::endl;

  std::cout << v_input1 << std::endl;
  std::cout << v_input2 << std::endl;
  std::cout << v_result << std::endl;

  int ok=0;
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>1.0e-6){
      std::cout<< "*****" << std::endl;
      std::cout<< "["<<i<<"] "<< abs(reference[i]-result[i]) << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) {
    std::cout << " OK!" <<std::endl;
  }
  assert(ok==0);
}

class funcPermute {
public:
  int n;
  funcPermute(int _n) { n=_n;};
  template<class vec>    void operator()(vec &rr,vec &i1,vec &i2) const { permute(rr,i1,n);}
  template<class scal>   void apply(std::vector<scal> &rr,std::vector<scal> &in)  const { 
    int sz=in.size();
    int msk = sz>>(n+1);
    for(int i=0;i<sz;i++){
      rr[i] = in[ i^msk ];
    }
  }
  std::string name(void) const { return std::string("Permute"); }
};



template<class scal, class vec,class functor > 
void PermTester(const functor &func)
{
  int Nsimd = vec::Nsimd();

  std::vector<scal> input1(Nsimd);
  std::vector<scal> input2(Nsimd);
  std::vector<scal> result(Nsimd);
  std::vector<scal> reference(Nsimd);

  std::vector<vec,alignedAllocator<vec> > buf(3);
  vec & v_input1 = buf[0];
  vec & v_input2 = buf[1];
  vec & v_result = buf[2];

  for(int i=0;i<Nsimd;i++){
    input1[i] = scal(drand48());
    input2[i] = scal(drand48());
    result[i] = scal(drand48());
  }

  merge<vec,scal>(v_input1,input1);
  merge<vec,scal>(v_input2,input2);
  merge<vec,scal>(v_result,result);

  func(v_result,v_input1,v_input2);

  func.apply(reference,input1);

  extract<vec,scal>(v_result,result);
  std::cout << " " << func.name() << " " <<func.n <<std::endl;

  int ok=0;
  if (0) {
    std::cout<< "*****" << std::endl;
    for(int i=0;i<Nsimd;i++){
      std::cout<< input1[i]<<" ";
    }
    std::cout <<std::endl; 
    for(int i=0;i<Nsimd;i++){
      std::cout<< result[i]<<" ";
    }
    std::cout <<std::endl; 
    for(int i=0;i<Nsimd;i++){
      std::cout<< reference[i]<<" ";
    }
    std::cout <<std::endl; 
    std::cout<< "*****" << std::endl;
  }
  for(int i=0;i<Nsimd;i++){
    if ( abs(reference[i]-result[i])>1.0e-7){
      std::cout<< "*****" << std::endl;      
      std::cout<< "["<<i<<"] "<< abs(reference[i]-result[i]) << " " <<reference[i]<< " " << result[i]<<std::endl;
      ok++;
    }
  }
  if ( ok==0 ) {
    std::cout << " OK!" <<std::endl;
  }
  assert(ok==0);
}


int main (int argc, char ** argv)
{
  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vRealF "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  Tester<RealF,vRealF>(funcPlus());
  Tester<RealF,vRealF>(funcMinus());
  Tester<RealF,vRealF>(funcTimes());
  Tester<RealF,vRealF>(funcConj());

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vRealF permutes "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vRealF::Nsimd();i++){
    PermTester<RealF,vRealF>(funcPermute(i));
  }

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vRealD "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  Tester<RealD,vRealD>(funcPlus());
  Tester<RealD,vRealD>(funcMinus());
  Tester<RealD,vRealD>(funcTimes());
  Tester<RealD,vRealD>(funcConj());


  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vRealD permutes "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vRealD::Nsimd();i++){
    PermTester<RealD,vRealD>(funcPermute(i));
  }

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vComplexF "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  Tester<ComplexF,vComplexF>(funcTimesI());
  Tester<ComplexF,vComplexF>(funcTimesMinusI());
  Tester<ComplexF,vComplexF>(funcPlus());
  Tester<ComplexF,vComplexF>(funcMinus());
  Tester<ComplexF,vComplexF>(funcTimes());
  Tester<ComplexF,vComplexF>(funcConj());

  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vComplexF permutes "<<std::endl;
  std::cout << "==================================="<<  std::endl;

  // Log2 iteration
  for(int i=0;(1<<i)< vComplexF::Nsimd();i++){
    PermTester<ComplexF,vComplexF>(funcPermute(i));
  }


  std::cout << "==================================="<<  std::endl;
  std::cout << "Testing vComplexD "<<std::endl;
  std::cout << "==================================="<<  std::endl;


  Tester<ComplexD,vComplexD>(funcTimesI());
  Tester<ComplexD,vComplexD>(funcTimesMinusI());
  Tester<ComplexD,vComplexD>(funcPlus());
  Tester<ComplexD,vComplexD>(funcMinus());
  Tester<ComplexD,vComplexD>(funcTimes());
  Tester<ComplexD,vComplexD>(funcConj());

  std::cout << "===================================" << std::endl;
  std::cout << "Testing vComplexD permutes " << std::endl;
  std::cout << "===================================" << std::endl;

  // Log2 iteration
  for (int i = 0; (1 << i) < vComplexD::Nsimd(); i++) {
    PermTester<ComplexD, vComplexD>(funcPermute(i));
  }

}
