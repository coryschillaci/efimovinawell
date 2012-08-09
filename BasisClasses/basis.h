#ifndef basis_H
#define basis_H

#include <iostream>
#include <vector>
#include <Eigen>

#include "HOstate.h"

using namespace std;

class basis {
public:
  basis(int set_nmax,int set_J);
  ~basis();

  void Init();
  
  // These functions return values of private class variables
  int nmax() {return nmax_;};
  int J() {return J_;};

  void printBasis();

protected:
  const int nmax_;
  const int J_;
  unsigned int dim;

  vector<HOstate> fullBasis;
  MatrixXd symmetrizer;

  void popBasis(); // Finds all basis states with l even and E<=nmax
  void symmetrizer();
  

};

#endif
