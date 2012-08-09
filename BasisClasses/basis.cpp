#include "basis.h"

// Constructor
basis::basis(int set_nmax,int set_J) : nmax_(set_nmax),J_(set_J) {
  Init();
}

// Destructor
basis::~basis() {
}

// Create the basis
void basis::Init() {

  popBasis(); // Calculate all HO states with even l and E<=nmax_+3;

}

void basis::popBasis() {

  HOstate temp(0,0,0,0,J_);

  dim=0;

  for(int nshell=0; nshell<=nmax_; nshell+=2) {
    for(int n=0; n<=nshell/2; n++) {

      temp.set_n(n);

      // Remember to require that l+s is even
      for(int l=0; l<=nshell-2*n; l+=2) {
       
	temp.set_l(l);

	for(int N=0; N <= (nshell-2*n-l)/2 ; N++) {

	  temp.set_N(N);
	  temp.set_L(nshell-2*(n+N)-l);
	  
	  if( temp.IsValid() ) {
	    fullBasis.push_back( temp );
	    dim++;
	  }
	    
	} // end N loop
      } // end l loop
    } // end n loop
  } // end nshell loop
  if ( dim!=fullBasis.size() ) cout<<"WARNING: error in size of basis"
}

void basis::printBasis() {

  int length=fullBasis.size();
  
  fullBasis[0].printHeader();
  cout<<"----------------\n";

  for(int i=0;i<length;i++) {

    fullBasis[i].printKet();
    // cout<<"\n";

  }
}

void basis::symmetrizer() {

}
