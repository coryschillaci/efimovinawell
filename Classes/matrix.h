#ifndef matrix_H
#define matrix_H

#include <iostream>
#include <vector>

using namespace std;

// Use c-indexing, i.e. j=0,1,2,3,...,dim_-1
// For internal purposes, the matrix is indexed as (i,j)= dim_*i+j

class matrix {

 public:
  matrix();
  matrix(unsigned int d);
  ~matrix();
  
  double get(unsigned int i,unsigned int j); // Retrieve an element. 
  void set(unsigned int i,unsigned int j,double x); // Set an element
 
  void reserve(unsigned int n); // Reserves space for an n x n matrix. 

  void print();
  

 protected:
  
  unsigned int dim_;
  
  vector<double> m_; // The matrix
  

};

#endif
