#include <iostream>
#include <cmath>

using namespace std; 

// f(r)
void f(double *f, const double*  const r, const double a, const double b, const double c){
  f[0]=a*(r[1]-r[0]);
  f[1]=r[0]*(b-r[2])-r[1];
  f[2]=(r[0]*r[1])-(c*r[2]);
}
   
 

int main(){

const double dt=0.001;
const double T= 100;
const int N= T/dt;

const double a=10.0;
const double b=28.0;
const double c=8/3;

double r[3];
r[0] = 1;
r[1] = 1;
r[2] = 1;

double k1[3];
double k2[3];
double k3[3];
double k4[3];
double rtemp[3];


for(int i=0;i<N;i++){

  f(k1,r,a,b,c);  // k1 = f(r_n)
  
  rtemp[0] = r[0] + dt/2 * k1[0];
  rtemp[1] = r[1] + dt/2 * k1[1];
  rtemp[2] = r[2] + dt/2 * k1[2];
  

  f(k2,rtemp,a,b,c); // k2 = f( r_n + dt/2 * k1)
    
  rtemp[0] = r[0] + dt/2 * k2[0];
  rtemp[1] = r[1] + dt/2 * k2[1];
  rtemp[2] = r[2] + dt/2 * k2[2];
  
  f(k3,rtemp,a,b,c); // k3 = f(r_n + dt/2*k2);

  rtemp[0] = r[0] + dt * k3[0];
  rtemp[1] = r[1] + dt * k3[1];
  rtemp[2] = r[2] + dt * k3[2];

  f(k4, rtemp,a,b,c); // k4 = f(r_n + dt * k3)
  
  r[0] = r[0] + dt/6 * (k1[0]+2*k2[0]+2*k3[0]+k4[0]);
  r[1] = r[1] + dt/6 * (k1[1]+2*k2[1]+2*k3[1]+k4[1]);
  r[2] = r[2] + dt/6 * (k1[2]+2*k2[2]+2*k3[2]+k4[2]);

   cout << (i+1)*dt <<'\t'<< r[0]<< "\t" << r[1] << "\t"<< r[2]<<endl;
  
}


 return 0;
}
