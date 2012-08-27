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
  popSymmetrizer(); // Create the matrix 1+\Pi
  calcSymBasis();

}

void basis::popBasis() {

  HOstate temp(0,0,0,0,J_);

  basisDim_ = 0;

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
	    basisDim_++;
	  }
	    
	} // end N loop
      } // end l loop
    } // end n loop
  } // end nshell loop
  if ( basisDim_ != fullBasis.size() ) cout<<"WARNING: error in size of basis"; 
}

void basis::printBasis() {

  int length=fullBasis.size();
  
  fullBasis[0].printHeader();
  cout<<"----------------\n";

  for(int i=0;i<length;i++) {

    fullBasis[i].printKet();

  }
  cout<<"\n";
}

void basis::printSymmetrizer() {

  cout<<"The symmetrizer matrix is \n"<<symmetrizer_<<"\n";

}

void basis::printStates() {

  cout<<"The basis of HO states with correct permutational symmetry contains "<<statesDim_<<" states\n";
  cout<<"----------------------------------------------------------------------------- \n\n";

  for(unsigned int i=0;i<statesDim_;i++) {
    
    for(unsigned int j=0;j<basisDim_;j++) {
      
      if(states_.col(i)(j)!=0) {
	cout<<setw(10)<<showpoint<<setprecision(6)<<states_.col(i)(j);
	fullBasis[j].printKet();	
      }

    } // End j loop
    
    cout<<"\n";

  } // End i loop
	     
}

// Sorts the states into increasing energy order
void basis::sortStates() {

  VectorXi energies(statesDim_);
  VectorXd temp(basisDim_);
  int tempIndex;
  int tempEnergy;

  // Find the energy of each basis state
  for(unsigned int i=0; i<statesDim_; i++) {

    for(unsigned int j=0; j<basisDim_; j++) {
      
      if(states_.col(i)(j)!=0) tempIndex=j;
     
    } // end j loop
    
    energies(i)=fullBasis[tempIndex].energy();
   
  } // End i loop

  // Go over sufficient times to move a column from right side to left side
  for(unsigned int counter=1; counter<=statesDim_; counter++) {
    
    for(unsigned int i=0; i<statesDim_-1; i++) {
      
      if( energies(i)>energies(i+1) ) {
	
	// Switch the columns of states_
	temp=states_.col(i);
	states_.col(i)=states_.col(i+1);
	states_.col(i+1)=temp;

	// Switch the saved energies
	tempEnergy=energies(i);
	energies(i)=energies(i+1);
	energies(i+1)=tempEnergy;

      } // End if

    } // End i loop

  } // End counter loop

}

void basis::popSymmetrizer() {

  HOstate prime,unprime;

  symmetrizer_.resize( basisDim_,basisDim_ );

  for(unsigned int i=0;i<basisDim_;i++) {
    
    prime=fullBasis[i];

    for(unsigned int j=i;j<basisDim_;j++) {

      unprime=fullBasis[j];

      symmetrizer_(i,j) = ( 2.0*TMB( unprime.N() , unprime.L() , unprime.n() , unprime.l() , 
				      prime.n() , prime.l() , prime.N() , prime.L(), 
				      prime.J(), 3) ) / 3.0;
      
      if(i==j) symmetrizer_(i,j)+=1.0/3.0; 
      else symmetrizer_(j,i)=symmetrizer_(i,j);
      

    } //end j loop
  } //end i loop

}


void basis::calcSymBasis() {

  SelfAdjointEigenSolver<MatrixXd> eigensolver(symmetrizer_);

  // double norm;

  statesDim_=0;
  
  // SHOULD ADD ERROR HANDLING
  //  if (eigensolver.info() != Success) abort();

  // cout<<"\nThe eigenvalues are: \n";
  // cout<<"--------------------------\n";
  // cout<<eigensolver.eigenvalues()<<"\n";
  // cout<<"--------------------------\n";

  for(unsigned int i=0;i<basisDim_;i++) {
    if( fabs( eigensolver.eigenvalues()(i)-1 )<=.001 ) {
      
      statesDim_++;
      states_.conservativeResize(basisDim_,statesDim_);

      states_.col(statesDim_-1)=eigensolver.eigenvectors().col(i);
      
      // Remove numerically small entries
      for(unsigned int j=0;j<basisDim_;j++) {
	if( fabs( states_.col(statesDim_-1)(j) ) < .0000001 ) states_.col(statesDim_-1)(j) = 0;
      }
      
      // // Need to normalize the state vectors;
      // norm=0;
      
      // for(unsigned int j=0;j<basisDim_;j++) {
	
      // 	if( fabs( states_.col(statesDim_-1)(j) ) < .0000001 ) states_.col(statesDim_-1)(j) = 0;
	
      // 	else norm+=states_.col(statesDim_-1)(j)*states_.col(statesDim_-1)(j);
      // }

      // states_.col(statesDim_-1)/=sqrt(norm);
      
    } // End if
  } // End i loop
  
}

