#ifndef basis_H
#define basis_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#include "HOstate.h"
#include "coupling.h"

using namespace std;
using namespace Eigen;

class basis {
public:
  basis(int set_nmax,int set_J);
  ~basis();

  void Init();
  
  // These functions return values of private class variables
  int nmax() {return nmax_;};
  int J() {return J_;};

  void printBasis();
  void printSymmetrizer();
  void printStates();

  void sortStates(); // Sorts the basis states in order of increasing energy

protected:
  const int nmax_;
  const int J_;
  unsigned int basisDim_;
  unsigned int statesDim_;

  vector<HOstate> fullBasis;
  MatrixXd symmetrizer_;
  MatrixXd states_;

  void popBasis(); // Finds all basis states with l even and E<=nmax
  void popSymmetrizer(); // Creates the matrix 1/3*(1+2*\Pi)
  void calcSymBasis();
  

};

#endif
