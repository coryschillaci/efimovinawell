#include "energies.h"

// Constructor
energies::energies(basis set_basis) : basis_(set_basis) {
  Init();
}

energies::energies(basis set_basis, double set_aInvMin, double set_aInvMax) : aMin_(1.0/set_aInvMin),aMax_(1.0/set_aInvMax),basis_(set_basis) {
  Init();
}

// Destructor
energies::~energies() {
}

void energies::Init() {

  // Initialize the matrices G0 and T12 to the identity 
  G0_.resize( basis_.basisDim(), 1 ); 
  G0_=MatrixXd::Identity( basis_.basisDim(),basis_.basisDim() );
  T12_.resize( basis_.basisDim(), 1 );
  
}

// Calculates the HO Green's function matrix at energy E
void energies::setG0(double E) {

  for(unsigned int i=0; i<basis_.basisDim(); i++) G0_(i,i)/= 1/(E - 2*basis_.getBasisState(i).n() - basis_.getBasisState(i).l() - 2*basis_.getBasisState(i).N()- basis_.getBasisState(i).L() - 3 );
		 
}

// Calculates the two body t-matrix for a regulated two-body contact interaction at energy E and interaction strength a
void energies::setT12(double E, double a) {

  HOstate state,primeState;

  double temp;

  for(unsigned int i=0; i<basis_.basisDim(); i++) {
    
    if(basis_.getBasisState(i).l()==0 ) {

      primeState=basis_.getBasisState(i);
      
      for (unsigned int j=0; j<basis_.basisDim(); j++ ) {
	
	state=basis_.getBasisState(j);

	// If N=N',L=L',l=l'=0 use 2-body t-matrix
	if( state.N()==primeState.N() && state.L()==primeState.L() && state.l()==0 ) {
	  
	  temp=( LogGamma(primeState.n()+1.5)+LogGamma(state.n()+1.5)-LogGamma(primeState.n()+1.0)-LogGamma(state.n()+1) )/2.0;
	  temp=exp( temp );

	  temp*=2*ROOT2/PI;
	  
	  temp*=C0( E, a,state.N(),state.L() );	  
	  
	}
	
	// Otherwise the entries are zero
	else temp=0;

	T12_(i,j)= temp;
	T12_(j,i)= temp;
	
      } // End j loop
    } // End if (l'==0)
 
    // Inserts zeros when l'!=0
    else for (unsigned int j=0; j<basis_.basisDim(); j++ ) {
	T12_(i,j)=0;
	T12_(j,i)=0;
      } // End else for

  } // End i loop
  
}

// Required by setT12, defined as in Tom's notes
double energies::C0(double E, double a, int N, int L) {

  double c0=0.0;
  int imax= (basis_.nmax()-2.0*N-L )/2.0;

  // Return zero for special case a=0
  if(a==0) return 0;

  // Otherwise calculate according to Tom's notes
  for(int i=0; i<imax; i++ ) {
    c0+=exp( LogGamma(i+1.5) - LogGamma(i+1) )/( 2.0*i+1.5-(E-2*N-L-1.5) ) ;
  }

  c0*= -2*ROOT2/PI;

  c0-= ROOT2*exp( LogGamma( .75-(E-2*N-L-1.5)/2.0 ) - LogGamma( .25-(E-2*N-L-1.5)/2.0 ) );

  c0+=1/a;

  return 1.0/c0;

}

void energies::calcPair(double a,double e, Vector2d * pair) {

  MatrixXd mm(basis_.statesDim(),basis_.statesDim());
  VectorXd vec(basis_.basisDim());
  VectorXd primeVec(basis_.basisDim());

  double best=123456789; // Start with an absurd number

  // Generate the matrix
  for(unsigned int i=0; i<basis_.statesDim(); i++) {

    primeVec=basis_.getState(i);

    for(unsigned int j=0; j<basis_.statesDim(); j++) {

      vec=basis_.getState(j);

      mm(i,j)=primeVec.transpose()*G0_*T12_*vec;

    }

  }

  // Calculate the eigenvalues
  SelfAdjointEigenSolver<MatrixXd> eigensolver(mm);

  // Find eigenvalue closes to 1
  for(unsigned int i=0; i<basis_.statesDim(); i++ ) {
    
    if( fabs(eigensolver.eigenvalues()(i)-1) < best ) best=eigensolver.eigenvalues()(i);
    
  }

  if(best==123456789) cout<<"Error, no good eigenvalues";

  (*pair)(0)=e;
  (*pair)(1)=best;
  

}

void energies::calcAllPairs(double aInvMin, double aInvMax, double aInvStep, double eMin, double eMax, double eStep) {

  numPairs_=0;

  for( double e=eMin; e<eMax; e+=eStep) {

    setG0(e);

    for( double aInv=aInvMin; aInv<aInvMax; aInv+=aInvStep) {

      setT12(e,a);

      numPairs++;

      eLambdaPairs_.conservativeResize(numPairs_);
      
      calcPair(a,e,&(eLambdaPairs_.row(numPairs_-1)));

    }

  }

}
