// Calcular funciones de Bessel usando integrales de Simpson
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const double ERR=1e-7;



double f(double t,double x,double Alpha){
  return cos(Alpha*t-x*sin(t));
}
double Simpson(double a,double b,int n,double x,double Alpha){
  n*=2; //multplicar n por 2
  double h=(b-a)/n;
  int i; double suma,t;
  for(suma=0,i=0;i<=n;i++){
    t=a+i*h;//calcular x;
    if(i==0 || i==n)
      suma+=f(t,x,Alpha);
    else if(i%2==0)
      suma+=2*f(t,x,Alpha);
    else
      suma+=4*f(t,x,Alpha);
  }
  return suma*h/3;
}

double Bessel(double Alpha,double x,int n){
  return 1.0/M_PI*Simpson(0,M_PI,n,x,Alpha);
}

double cerosPorBiseccion(double a,double b, double Alpha, double n){
  double m,fa,fm;
  //calcular fa;
  fa=Bessel(Alpha,a,n);
  while(b-a>ERR){
    //calcular m y fm;
    m=(a+b)/2; fm=Bessel(Alpha,m,n);
    if(fa*fm>0)
      {a=m; fa=fm;} //mover a hacia m
    else
      b=m;         //mover b hacia m;
  }
  return (a+b)/2;
}
int main(void){
  ofstream myfile;
  double Alpha=0;

  int n=50; double x;
  myfile.open ("Bessel_teorico.txt");
  myfile << "gamma" << std::setw(10) << "Bessel"<<"\n";
  for(x=0.1;x<15;x+=0.01)
    myfile << x<< std::setw(10) << Bessel(Alpha,x,n)<<"\n";

  myfile.close();


  myfile.open ("Bessel_teorico_raices.txt");
  myfile<<"Raices teoricas"<<"\n";
  myfile<<cerosPorBiseccion(0,3,Alpha,n)<<"\n";
  myfile<<cerosPorBiseccion(3,8,Alpha,n)<<"\n";
  myfile<<cerosPorBiseccion(8,10,Alpha,n)<<"\n";
  myfile<<cerosPorBiseccion(10,14,Alpha,n)<<"\n";
  myfile<<cerosPorBiseccion(14,15,Alpha,n)<<"\n";
  myfile.close();


  cout<<"set terminal png"<<endl;
  cout<<"set output 'Bessel_teorico.png'"<<endl;
  cout<<"set ylabel \"Bessel\""<<endl;
  cout<<"set xlabel \"Lambda\""<<endl;
  cout<<"set title \"Bessel Torico\""<<endl;
  cout<<"set key noautotitle"<<endl;
  cout<<"set xrange[0.1:15]"<<endl;
  cout<<"plot 'Bessel_teorico.txt' u 1:2 w l"<<endl;






  return 0;
}
