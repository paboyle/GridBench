const int Nc=3;

#define sitePhi(sp,co) Phi[idx*4*Nc + sp*Nc+co]

#define siteU(c1,c2)   U  [site*Nc*Nc*8 + mu*Nc*Nc + c1*Nc+c2]

#define SUN_MULT						\
  for(int c1=0;c1<Nc;c1++){					\
    siteUChi[0][c1] =						\
      siteU(c1,0)*siteChi[0][0]					\
      + siteU(c1,1)*siteChi[0][1]				\
      + siteU(c1,2)*siteChi[0][2];				\
    siteUChi[1][c1] =						\
      siteU(c1,0)*siteChi[1][0]					\
      + siteU(c1,1)*siteChi[1][1]				\
      + siteU(c1,2)*siteChi[1][2];				\
  }

template<class Complex>
void dslash_kernel(Complex * RESTRICT Up,Complex * RESTRICT Psip,Complex * RESTRICT Phip,uint64_t * RESTRICT nbr,uint64_t nsite,uint64_t Ls,uint8_t * RESTRICT prm)
{
  Complex * RESTRICT U  = (Complex *)Up;
  Complex * RESTRICT Psi= (Complex *)Psip;
  Complex * RESTRICT Phi= (Complex *)Phip;
  Complex complex_i(0.0,1.0);

#ifdef OMP
#pragma omp parallel for
#endif 
  for(uint64_t ssite=0;ssite<nsite;ssite++){
    uint64_t site = ssite;
    Complex sitePsi [4][Nc];
    Complex siteChi [2][Nc];
    Complex siteUChi [2][Nc];

    for(uint64_t s=0;s<Ls;s++){

      uint64_t site5 = site*Ls+s;
      for(uint64_t mu=0;mu<8;mu++){
	uint64_t idx = nbr[site5*8+mu]; 

	// Two spin projection
	if ( mu == 0 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) - complex_i * sitePhi(3,co);
	    siteChi[1][co] = sitePhi(1,co) - complex_i * sitePhi(2,co);
	  }
	} else if ( mu == 1 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) + sitePhi(3,co);
	    siteChi[1][co] = sitePhi(1,co) - sitePhi(2,co);
	  }
	} else if ( mu == 2 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) - complex_i * sitePhi(2,co);
	    siteChi[1][co] = sitePhi(1,co) + complex_i * sitePhi(3,co);
	  }
	} else if ( mu == 3 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) - sitePhi(2,co);
	    siteChi[1][co] = sitePhi(1,co) - sitePhi(3,co);
	  }
	} else if ( mu == 4 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) + complex_i * sitePhi(3,co);
	    siteChi[1][co] = sitePhi(1,co) + complex_i * sitePhi(2,co);
	  }
	} else if ( mu == 5 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) - sitePhi(3,co);
	    siteChi[1][co] = sitePhi(1,co) + sitePhi(2,co);
	  }
	} else if ( mu == 6 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) + complex_i * sitePhi(2,co);
	    siteChi[1][co] = sitePhi(1,co) - complex_i * sitePhi(3,co);
	  }
	} else if ( mu == 7 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi(0,co) + sitePhi(2,co);
	    siteChi[1][co] = sitePhi(1,co) + sitePhi(3,co);
	  }
	}

	// SU(N) multiply
	SUN_MULT;

	if ( mu == 0 ) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[0][co] = siteUChi[0][co];
	    sitePsi[1][co] = siteUChi[1][co];
	  }
	} else { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[0][co]+= siteUChi[0][co];
	    sitePsi[1][co]+= siteUChi[1][co];
	  }
	}

	// Four spin reconstruction
	if ( mu == 0 ) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]= complex_i * siteUChi[1][co];
	    sitePsi[3][co]= complex_i * siteUChi[0][co];
	  }
	} else if ( mu == 1) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]-= siteUChi[1][co];
	    sitePsi[3][co]+= siteUChi[0][co];
	  }
	} else if ( mu == 2) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]+= complex_i* siteUChi[0][co];
	    sitePsi[3][co]-= complex_i* siteUChi[1][co];
	  }
	} else if ( mu == 3) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]-= siteUChi[0][co];
	    sitePsi[3][co]-= siteUChi[1][co];
	  }
	} else if ( mu == 4) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]-=complex_i * siteUChi[1][co];
	    sitePsi[3][co]-=complex_i * siteUChi[0][co];
	  }
	} else if ( mu == 5) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]+= siteUChi[1][co];
	    sitePsi[3][co]-= siteUChi[0][co];
	  }
	} else if ( mu == 6) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]-= complex_i* siteUChi[0][co];
	    sitePsi[3][co]+= complex_i* siteUChi[1][co];
	  }
	} else if ( mu == 7) { 
	  for(int co=0;co<Nc;co++){
	    sitePsi[2][co]+= siteUChi[0][co];
	    sitePsi[3][co]+= siteUChi[1][co];
	  }
	}
	for(int sp=0;sp<4;sp++){
	for(int co=0;co<Nc;co++){
	  Psi[site5 *4*Nc + sp*Nc+co] = sitePsi[sp][co];
	}}

      }
    }
  }
}


// scalar version optimised on a Skylake laptop.
// Unrolled mu loop. Probably not necessary.
// Removing temporary stack arrays on U and Phi helps; encode memrefs.

template<class Complex>
void dslash_kernel_unroll(Complex * RESTRICT Up,Complex * RESTRICT Psip,Complex * RESTRICT Phip,uint64_t * RESTRICT nbr,uint64_t nsite,uint64_t Ls,uint8_t * RESTRICT prm)
{
  Complex * RESTRICT U  = (Complex *)Up;
  Complex * RESTRICT Psi= (Complex *)Psip;
  Complex * RESTRICT Phi= (Complex *)Phip;
  Complex complex_i(0.0,1.0);

#ifdef OMP
#pragma omp parallel for
#endif
  for(uint64_t ssite=0;ssite<nsite;ssite++){
    uint64_t site = ssite;
    Complex sitePsi [4][Nc];
    Complex siteChi [2][Nc];
    Complex siteUChi [2][Nc];

    for(uint64_t s=0;s<Ls;s++){

      uint64_t site5 = site*Ls+s;
      uint64_t mu,idx;

      mu=0; idx = nbr[site5*8+mu]; 
      // Two spin projection
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) - complex_i * sitePhi(3,co);
	siteChi[1][co] = sitePhi(1,co) - complex_i * sitePhi(2,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co] = siteUChi[0][co];
	sitePsi[1][co] = siteUChi[1][co];
	sitePsi[2][co] = complex_i * siteUChi[1][co];
	sitePsi[3][co] = complex_i * siteUChi[0][co];
      }

      mu=1; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) + sitePhi(3,co);
	siteChi[1][co] = sitePhi(1,co) - sitePhi(2,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]-= siteUChi[1][co];
	sitePsi[3][co]+= siteUChi[0][co];
      }

      mu=2; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) - complex_i * sitePhi(2,co);
	siteChi[1][co] = sitePhi(1,co) + complex_i * sitePhi(3,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]+= complex_i* siteUChi[0][co];
	sitePsi[3][co]-= complex_i* siteUChi[1][co];
      }

      mu=3; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) - sitePhi(2,co);
	siteChi[1][co] = sitePhi(1,co) - sitePhi(3,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]-= siteUChi[0][co];
	sitePsi[3][co]-= siteUChi[1][co];
      }

      mu=4; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) + complex_i * sitePhi(3,co);
	siteChi[1][co] = sitePhi(1,co) + complex_i * sitePhi(2,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]-=complex_i * siteUChi[1][co];
	sitePsi[3][co]-=complex_i * siteUChi[0][co];
      }

      mu=5; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) - sitePhi(3,co);
	siteChi[1][co] = sitePhi(1,co) + sitePhi(2,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]+= siteUChi[1][co];
	sitePsi[3][co]-= siteUChi[0][co];
      }

      mu=6; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) + complex_i * sitePhi(2,co);
	siteChi[1][co] = sitePhi(1,co) - complex_i * sitePhi(3,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]-= complex_i* siteUChi[0][co];
	sitePsi[3][co]+= complex_i* siteUChi[1][co];
      }

      mu=7; idx = nbr[site5*8+mu]; 
      for(int co=0;co<Nc;co++){
	siteChi[0][co] = sitePhi(0,co) + sitePhi(2,co);
	siteChi[1][co] = sitePhi(1,co) + sitePhi(3,co);
      }
      SUN_MULT;
      for(int co=0;co<Nc;co++){
	sitePsi[0][co]+= siteUChi[0][co];
	sitePsi[1][co]+= siteUChi[1][co];
	sitePsi[2][co]+= siteUChi[0][co];
	sitePsi[3][co]+= siteUChi[1][co];
      }

      for(int sp=0;sp<4;sp++){
	Psi[site5 *4*Nc + sp*Nc+0] = sitePsi[sp][0];
	Psi[site5 *4*Nc + sp*Nc+1] = sitePsi[sp][1];
	Psi[site5 *4*Nc + sp*Nc+2] = sitePsi[sp][2];
      }
    }
  }
}
