#define ROOT_PI 1.77245385090552

#include <math.h>
#include <iostream>
using namespace std;

#include "gamma.h"
#include "coupling.h"

// Simple function to return -1^n
int parity (int n) {
  return ( n % 2 == 0 ) ? 1:-1;
}

// Simple function to test whether an integer is even
bool IsEven (int n) {
  return ( n % 2 == 0 ) ? 1:0;
}

// Simple function to test whether an integer is even
bool IsInteger (float x){
  return ( x- int(x) == 0 ) ? 1:0;
}

// Tests whether the triad (l1,l2,L) satisfies the triangular inequality |l1-l2|<L<l1+l2
bool IsTriangular (float l1,float l2,float L) {
  return ( fabs(l1-l2)<=L && L<=l1+l2 ) ? 1:0;
}

// Return the greater element of the pair (a,b)
int max(int a, int b) {return ( a>b ) ? a:b;}
float max(float a,float b) {return ( a>b ) ? a:b;}
double max(double a,double b) {return ( a>b ) ? a:b;}

// Return the less element of the pair (a,b)
int min(int a, int b) {return ( a<b ) ? a:b;}
float min(float a,float b) {return ( a<b ) ? a:b;}
double min(double a,double b) {return ( a<b ) ? a:b;}

// Return the absolute value of an integer as an integer
int iabs(int n) {return (n<0) ? -n:n;}

// Find the log of the triangle coefficient (a+b-c)!*(a-b+c)!*(-a+b+c)!/(a+b+c+1)!
double LogTriangle (float a,float b,float c) {
  return LogFactorial(a+b-c)+LogFactorial(a-b+c)+LogFactorial(-a+b+c)-LogFactorial(a+b+c+1);
}

// Calculate the 3j symbol when m1=m2=M=0 
double ThreeJ (int l1,int l2,int L) {
  
  int g=(L+l1+l2)/2;
  double temp;

  // Require that |l1-l2|<L<l1+l2
  if ( !IsTriangular(l1,l2,L) ) return 0;

  // This coefficient is nonzero only if g is an integer
  if ( IsEven(L+l1+l2) )  {
    temp = LogFactorial(g)-LogFactorial(g-l1)-LogFactorial(g-l2)-LogFactorial(g-L);
    temp+= 0.5*(LogFactorial(2*g-2*l1)+LogFactorial(2*g-2*l2)+LogFactorial(2*g-2*L)-LogFactorial(2*g+1));
    return parity(g)*exp( temp );
  }

  else return 0;

}

// Calculate the CG coefficient <l1 0, l2 0|L 0>
double ClebschGordan (int l1,int l2,int L) {
  return parity(l1-l2)*sqrt(2*L+1)*ThreeJ(l1,l2,L);
}


// The six j symbol {j1  j2  j3
//                   j4  j5  j6}
// TO DO: add protection for non integer/half integer arguments
double SixJ(float j1,float j2,float j3,float j4,float j5,float j6) {

  double sum=0,coeff;
  
  int imin,imax;
  
  // Find limits of summation
  imin=int(j1+j2+j3);
  imax=int(j1+j2+j4+j5);
  imin=max( imin, max( int(j1+j5+j6), max( int(j2+j4+j6) , int(j3+j4+j5) ) ) );
  imax=min( imax, min( int(j2+j3+j5+j6) , int(j1+j3+j4+j6) ) );

  // Enforce that (i,j,k) satisfies the triangular inequality and that i+j+k is an integer for the triads
  // (j1,j2,j3)    (j1,j5,j6)    (j4,j2,j6)    (j4,j5,j3)
  if( !(IsTriangular(j1,j2,j3) && IsInteger(j1+j2+j3) && IsTriangular(j1,j5,j6) && IsInteger(j1+j5+j6) && IsTriangular(j4,j2,j6) && IsInteger(j4+j2+j6) && IsTriangular(j4,j5,j3) && IsInteger(j4+j5+j3)) ) return 0;

  // This is the coefficient in front of the sum
  coeff=LogTriangle(j1,j2,j3)+LogTriangle(j1,j5,j6)+LogTriangle(j4,j2,j6)+LogTriangle(j4,j5,j3);
  coeff=exp(coeff/2);

  // Compute the sum
  for(int i=imin;i<=imax;++i) {
    sum+=parity(i)*exp( LogFactorial(i+1)-LogFactorial(i-j1-j2-j3)-LogFactorial(i-j1-j5-j6)-LogFactorial(i-j2-j4-j6)-LogFactorial(i-j3-j4-j5)-LogFactorial(j1+j2+j4+j5-i)-LogFactorial(j2+j3+j5+j6-i)-LogFactorial(j1+j3+j4+j6-i) );
  }

  return sum*coeff;

}

// The 9j symbol {j1  j2   J12
//                j3  j4   J34
//                J13 J24  J  }
double NineJ(float j1,float j2,float J12,float j3,float j4,float J34,float J13,float J24,float J) {
  
  double sum=0;

  int imin,imax;

  imin = max( fabs(J24-j3), max( fabs(j1-J),fabs(j2-J34) ) );
  imax = min( J24+j3, min( j1+J, j2+J34 ) );

  for(int i=imin;i<=imax;++i) {
    sum+=parity(2*i)*(2*i+1)*SixJ(j1,j2,J12,J34,J,i)*SixJ(j3,j4,J34,j2,i,J24)*SixJ(J13,J24,J,i,j1,j3);
  }

  return sum;
  
}

// TM Bracket calculated using Rashid's expression
// < nl, NL, J | n1 l1, n2 l2, J>_d
double TMB (int n,int l,int N,int L,int n1,int l1,int n2,int l2,int J,float d) {

  double coeff;
  double ksum=0;

  // Enforce energy conservation
  if(2.0*(n+N)+l+L!=2.0*(n1+n2)+l1+l2) return 0;
  
  // Require physical angular momentum couplings
  if( (!IsTriangular(l1,l2,J)) || (!IsTriangular(L,l,J)) ) return 0;

  // Factorials and gammas from line 3 of Eqn 40
  coeff=LogFactorial(n1)+LogFactorial(n2)+LogFactorial(n)-LogFactorial(N);
  coeff+=LogGamma(n1+l1+1.5)+LogGamma(n2+l2+1.5)-LogGamma(n+l+1.5)-LogGamma(N+L+1.5);
  //  coeff/=2.0; ERROR IN RASHID!!!!
 
  // Factorials from line 2 of Eqn 40
  coeff+=LogFactorial(l1+l2+J+1)+LogFactorial(l1+l2-J)+LogFactorial(l1-l2+J);
  coeff+=LogFactorial(l+L+J+1)+LogFactorial(l+L-J)+LogFactorial(l-L+J);
  coeff-=LogFactorial(l2-l1+J)+LogFactorial(L-l+J);
  coeff/=2.0;

  // exponentiate and multiply by (-1)^N from line 3
  coeff=parity(N)*exp(coeff);

  // Line 1 of Eqn 40
  coeff*=2*ROOT_PI*sqrt( (l1+0.5)*(l2+0.5)*(l+0.5)*(L+0.5) );
  
  for (int k=0;k<=l;k++) {
    
    ksum+=sum1(l,L,l1,l2,J,d,k)*sum2(l,N,L,n1,l1,n2,l2,d,k);
    
  }

  return coeff*ksum;

}

double sum1 (int l,int L,int l1,int l2,int J,float d,int k) {
  
  int lam2;

  int lam12min,lam12max;
  int zmin,zmax;

  double temp,total=0;

  for (int lam1=0;lam1<=l-k;lam1++) {
    
    lam2=l-k-lam1; // From Kronecker delta

    // Set the range for the lam12 index
    lam12min=max( iabs(l+l1-l2-k-2*lam1) , l-l1-l2-k );
    lam12min=max( lam12min , iabs(L-k)               );

    lam12max=min( l1+l2-l+k , L+k );

    for(int lam12=lam12min;lam12<=lam12max;lam12++) {

      // Factorials, from angular momentum algebra, require some combinations to be even
      if ( IsEven(l1+l2+l+k+lam12) && IsEven (L+k+lam12) ) {
	
	zmin=max( 0    , l-k-lam12-J           );
	zmin=max( zmin , l1-l2-J               );
	zmin=max( zmin , l1-l2-lam1+lam2+l-k-J );
	
	zmax= l-k+lam12-J ;
       	
	for(int z=zmin;z<=zmax;z++) {
	  
	  // The factorials from line 5
	  temp= -LogFactorial(lam1)-LogFactorial(lam2);
	  
	  // The factorials from line 6
	  temp+=LogFactorial(k)+LogFactorial(-l1+l2+lam1-lam2+lam12)+LogFactorial( (l1+l2-l+k+lam12)/2 );
	  temp-=LogFactorial(l1+l2-l+k+lam12+1)+LogFactorial( (l1+l2-l+k-lam12)/2.0 )+LogFactorial( (l1-l2-lam1+lam2+lam12)/2.0 )+LogFactorial( (-l1+l2+lam1-lam2+lam12)/2.0 );
	  
	  // The factorials from line 7
	  temp+=LogFactorial(L-k+lam12)+LogFactorial( (L+k+lam12)/2.0 );
	  temp-=LogFactorial(L+k+lam12+1)+LogFactorial( (L+k-lam12)/2.0 )+LogFactorial( (L-k+lam12)/2.0 )+LogFactorial( (-L+k+lam12)/2.0 );
	  
	  // The factorials from line 8
 	  temp+=LogFactorial(-l+k+lam12+J+z)+LogFactorial(-l1+l2+J+z);
	  temp-=LogFactorial(z)+LogFactorial(l-k+lam12-J-z)+LogFactorial(-l1+l2+lam1-lam2-l+k+J+z)+LogFactorial(2*J+1+z);
	  
	  // Exponentiate
	  temp=exp(temp);
	  
	  // Multipliers from lines 4 and 5
	  temp*=parity( lam1+z+(l1+l2+l-L)/2 )*pow(d,(l1-lam1+lam2+k)/2.0)*pow(2, 2.0*(k-l) )*(2*lam12+1)/pow(1+d, (l1+l2)/2.0 );
	 
	  // Add this term to the sum
	  total+=temp;

	  // cout<<"sum2 loop has total of "<<total<<".\n";
	

	}// End z loop	
      
      } // end if (IfEven)
	   
    } // End lam12 loop
    
  }// End lam1 loop
  
  return total;
  
}

double sum2 (int l,int N,int L,int n1,int l1,int n2,int l2,float d,int k) {

  int t2min,p1min,p1max;

  double temp,total=0;

  for (int t1=0;t1<=n1;++t1) {
    
    t2min=max( 0 , (L+l-l1-l2)/2 -t1 + N);
   
    for (int t2=t2min;t2<=n2;t2++) {
      
      p1max=min( k, t1 );
      p1min=max(k-t2,0);

      for (int p1=p1min;p1<=p1max;p1++) {

	// Factorials from line 9
	temp = -LogFactorial(p1)-LogFactorial(k-p1)-LogFactorial(t1-p1)-LogFactorial(t2-k+p1);

	// Factorials from line 10;
	temp+=LogFactorial( t1+t2+(l1+l2-l-L)/2 )+LogGamma( t1+t2+(l1+l2-l+L+3)/2.0 );
	temp-=LogFactorial(n1-t1)+LogFactorial(n2-t2)+LogGamma(t1+l1+1.5)+LogGamma(t2+l2+1.5)+LogFactorial(t1+t2-N+(l1+l2-l-L)/2);
       
	// Exponentiate
	temp=exp(temp);
	
	// Multiplicative factors from line 9
	temp*=parity(t1+t2+p1)*pow(d,t1-p1)/pow(1+d,t1+t2);
	
	// add to the total
	total+=temp;

      }// end p1 loop
    }// end t2 loop    
  }// end t1 loop

  return total;

} 

