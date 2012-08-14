#include <iostream>
#include <iomanip>
#include <math.h>
#include "gamma.h"
#include "coupling.h"
#include "basis.h"

#define PI 3.14159265358979
#define ROOT_PI 1.77245385090552
#define DOUBLE_INFINITY 1.7e308

using namespace std;

int main () {

  // ************************************************************************************
  // Test of log factorial function
  //float x;
  //
  // cout<<"Calculating log factorial of\n";
  // cin>>x;  
  // cin.ignore();
  // cout<<x<<"!="<< LogFactorial(x)<<"\n";


  // ************************************************************************************
  // Test of 3j Symbols/CG Coefficients
  // int l1,l2,L;

  // cout<<"Calculating CG <l1 0, l2 0|L 0>\n";
  // cout<<"Enter l1: ";
  // cin>>l1;
  // cin.ignore();
  // cout<<"\nEnter l2: ";
  // cin>>l2;
  // cin.ignore();
  // cout<<"\nEnter L: ";
  // cin>>L;
  // cin.ignore();
  // cout<<"\n< "<<l1<<" 0 ; "<<l2<<" 0 | "<<L<<" 0 >= "<<ClebschGordan(l1,l2,L)<<"\n";

  // ************************************************************************************
  // Test of 6j Symbols
  // float j1=11.5,j2=10.5,j3=10,j4=1.5,j5=10.5,j6=10;

  // cout<<"{"<<j1<<" "<<j2<<" "<<j3<<"\n "<<j4<<" "<<j5<<" "<<j6<<"} = "<<SixJ(j1,j2,j3,j4,j5,j6)<<"\n";

  // ************************************************************************************
  // Test of 9j Symbols

  // float j1=100,j2=101,J12=13,j3=100,j4=100,J34=12,J13=10,J24=10,J=17;

  // cout<<"{"<<j1<<" "<<j2<<" "<<J12<<"\n "<<j3<<" "<<j4<<" "<<J34<<"\n "<<J13<<" "<<J24<<" "<<J<<"} = "<<NineJ(j1,j2,J12,j3,j4,J34,J13,J24,J)<<"\n";


  // ************************************************************************************
  // Test of TM Brackets
  // int n=3,l=1,N=0,L=2,n1=1,l1=2,n2=2,l2=1,J=2;
  // float d=2;

  // cout<<"\n"<<"d="<<d<<"\n\n";
 
  // cout<<"< n l, N L, J | n1 l1, n2 l2, J>_d\n";
  // cout<<"< "<<n<<","<<l<<" ; "<<N<<","<<L<<" ; "<<J<<" | "<<n1<<","<<l1<<" ; "<<n2<<","<<l2<<" ; "<<J<<" >="<<TMB(n,l,N,L,n1,l1,n2,l2,J,d)<<"\n\n";
  
  // cout<<"< N L, n l, J | n1 l1, n2 l2, J>_d\n";
  // cout<<"< "<<N<<","<<L<<" ; "<<n<<","<<l<<" ; "<<J<<" | "<<n1<<","<<l1<<" ; "<<n2<<","<<l2<<" ; "<<J<<" >="<<TMB(N,L,n,l,n1,l1,n2,l2,J,d)<<"\n\n";

  // cout<<"< N L, n l, J | n2 l2, n1 l1, J>_d\n";
  // cout<<"< "<<N<<","<<L<<" ; "<<n<<","<<l<<" ; "<<J<<" | "<<n2<<","<<l2<<" ; "<<n1<<","<<l1<<" ; "<<J<<" >="<<TMB(N,L,n,l,n2,l2,n1,l1,J,d)<<"\n\n";

  // cout<<"< n l, N L, J | n2 l2, n1 l1, J>_d\n";
  // cout<<"< "<<n<<","<<l<<" ; "<<N<<","<<L<<" ; "<<J<<" | "<<n2<<","<<l2<<" ; "<<n1<<","<<l1<<" ; "<<J<<" >="<<TMB(n,l,N,L,n2,l2,n1,l1,J,d)<<"\n\n";

  

  basis NewBasis(50,0);

  // NewBasis.printBasis();

  // NewBasis.printSymmetrizer();

  NewBasis.sortStates();
  
  NewBasis.printStates();
  
}


