const int Nc=3;

template<class Complex>
void dslash_kernel(Complex *Up,Complex *Psip,Complex *Phip,uint64_t *nbr,uint64_t nsite,uint64_t Ls,uint8_t *prm)
{
  Complex *U  = (Complex *)Up;
  Complex *Psi= (Complex *)Psip;
  Complex *Phi= (Complex *)Phip;

  for(uint64_t ssite=0;ssite<nsite;ssite++){
    uint64_t site = ssite;
    for(uint64_t s=0;s<Ls;s++){

      Complex sitePsi [4][Nc];
      Complex siteChi [2][Nc];
      Complex siteU   [Nc][Nc];
      Complex sitePhi [4][Nc];
      Complex siteUChi [2][Nc];
      Complex complex_i(0.0,1.0);

      uint64_t site5 = site*Ls+s;
      for(uint64_t mu=0;mu<8;mu++){
	uint64_t idx = nbr[site5*8+mu]; //idx = idx+s;

	for(int sp=0;sp<4;sp++){
	for(int co=0;co<Nc;co++){
	  sitePhi[sp][co] = Phi[ idx*4*Nc + sp*Nc+ co];
	}}
	// Two spin projection
	if ( mu == 0 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] - complex_i * sitePhi[3][co];
	    siteChi[1][co] = sitePhi[1][co] - complex_i * sitePhi[2][co];
	  }
	} else if ( mu == 1 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] + sitePhi[3][co];
	    siteChi[1][co] = sitePhi[1][co] - sitePhi[2][co];
	  }
	} else if ( mu == 2 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] - complex_i * sitePhi[2][co];
	    siteChi[1][co] = sitePhi[1][co] + complex_i * sitePhi[3][co];
	  }
	} else if ( mu == 3 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] - sitePhi[2][co];
	    siteChi[1][co] = sitePhi[1][co] - sitePhi[3][co];
	  }
	} else if ( mu == 4 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] + complex_i * sitePhi[3][co];
	    siteChi[1][co] = sitePhi[1][co] + complex_i * sitePhi[2][co];
	  }
	} else if ( mu == 5 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] - sitePhi[3][co];
	    siteChi[1][co] = sitePhi[1][co] + sitePhi[2][co];
	  }
	} else if ( mu == 6 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] + complex_i * sitePhi[2][co];
	    siteChi[1][co] = sitePhi[1][co] - complex_i * sitePhi[3][co];
	  }
	} else if ( mu == 7 ) { 
	  for(int co=0;co<Nc;co++){
	    siteChi[0][co] = sitePhi[0][co] + sitePhi[2][co];
	    siteChi[1][co] = sitePhi[1][co] + sitePhi[3][co];
	  }
	}

	// SU(N) multiply
	for(int c1=0;c1<Nc;c1++){
	for(int c2=0;c2<Nc;c2++){
	  siteU[c1][c2] = U[ site*Nc*Nc*8 + mu*Nc*Nc + c1*Nc+c2];
	}}
	for(int c1=0;c1<Nc;c1++){
	  siteUChi[0][c1] = siteU[c1][0]*siteChi[0][0];
	  siteUChi[1][c1] = siteU[c1][0]*siteChi[1][0];
	  for(int c2=1;c2<Nc;c2++){
	    siteUChi[0][c1]+= siteU[c1][c2]*siteChi[0][c2];
	    siteUChi[1][c1]+= siteU[c1][c2]*siteChi[1][c2];
	  }
	}

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
