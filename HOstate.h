/*******************************************************************************
This class creates harmonic oscillator states of the form | n l, N L, J >

-Upon construction, it verifies that the quantum numbers are physical. If they 
 are not a warning is printed.

-User visualization of state is handled via printHeader() and printKet()

-The CoM energy for a state in HO units can be obtained from the energy() 
 function

*******************************************************************************/


#ifndef HOstate_H
#define HOstate_H

#include <iostream>
#include <math.h>

using namespace std;

class HOstate {

public: 
  HOstate();
  HOstate(int set_n, int set_l,int set_N, int set_L,int set_J);
  ~HOstate();

  void Init(); // Initializer
  
  void printHeader(); // Prints | nl, NL, J>
  void printKet(); // Print a nicely formatted ket expression;

  // Return energy of the state
  int energy();

  // Print out any problems with current quantum numbers
  bool IsValid() {return valid_;}
  void IsValidVerbose();

  // These functions set values of private class variables
  void set_n(int new_n);
  void set_l(int new_l);
  void set_N(int new_N);
  void set_L(int new_L);
  void set_J(int new_J);
 
  // These functions return values of private class variables
  int n() {return n_;}
  int l() {return l_;}
  int N() {return N_;}
  int L() {return L_;}
  int J() {return J_;}
  
protected:
  int n_;
  int l_;
  int N_;
  int L_;
  int J_;
  
  bool valid_;

  // Utility that tests whether the triad (l1,l2,L) satisfies the triangular inequality |l1-l2|<L<l1+l2
  bool IsTriangular (int l1,int l2,int L);
  bool checkState();
  
};

#endif
