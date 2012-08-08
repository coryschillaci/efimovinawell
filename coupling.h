#ifndef COUPLING_H
#define COUPLING_H

// Simple function to return -1^n
int parity (int n);
// Simple function to test whether an integer is even
bool IsEven (int n);
// Simple function to test whether an integer is even
bool IsInteger (float x);
// Tests whether the triad (l1,l2,L) satisfies the triangular inequality |l1-l2|<L<l1+l2
bool IsTriangular (float l1,float l2,float L);

// Return the greater element of the pair (a,b)
int max(int a, int b);
float max(float a,float b);
double max(double a,double b);

// Return the lesser element of the pair (a,b)
int min(int a, int b);
float min(float a,float b);
double min(double a,double b);

// Return the absolute value of an integer as an integer
int iabs(int n);

// Computes the logarithm of the triangle coefficient 
double LogTriangle (float a,float b,float c);

// Versions of the two angular momentum coupling coefficients (3j Symbol and CG coefficient) when m1=m2=M=0 
double ThreeJ (int l1,int l2,int L);
double ClebschGordan (int l1,int l2,int L);

// The 6j symbol {j1  j2  j3
//                j4  j5  j6}
double SixJ(float j1,float j2,float j3,float j4,float j5,float j6);

// The 9j symbol {j1  j2   J12
//                j3  j4   J34
//                J13 J24  J  }
double NineJ(float j1,float j2,float J12,float j3,float j4,float J34,float J13,float J24,float J);

// The Talmi-Moshinsky Bracket
double TMB (int n,int l,int N,int L,int n1,int l1,int n2,int l2,int J,float d);

// Separate functions to calculate the first and second sums of Rashid
double sum1 (int l,int L,int l1,int l2,int J,float d,int k);
double sum2 (int l,int N,int L,int n1,int l1,int n2,int l2,float d,int k);

#endif
