#ifndef THERMOLIB_H
#define THERMOLIB_H
#include "secant.hpp"
using namespace std;

class thermolib
{
   private:
      bool success;
	  int dim, i, j;
      double pressure, temperature, molevolume, mole,*molefrac;
	  double *Pc, *Tc, *omega, Z, tmp, Tr, Pr;
       void construct();
      double a(int i) {return (0.42748*pow(8.8144,2)*pow(Tc[i],2)*pow(1.0+f_omega(i)*(1.0-sqrt(temperature/Tc[i])), 2)/Pc[i]);}
	  double a_mix();
       double b(int i) {return (0.08664*8.3144*Tc[i]/Pc[i]);}
	   double b_mix();
	  double A(int i) {return (a(i)*pressure/pow(8.3144, 2)/pow(temperature, 2));}
	  double A() {return (a_mix()*pressure/pow(8.3144, 2)/pow(temperature, 2));}
      double B(int i) {return (b(i)*pressure/8.3144/temperature);}
	  double B(){return (b_mix()*pressure/8.3144/temperature);}
	  double Zv(), phiV(int), phiL(int);
	  double f_omega(int i) {return (0.48 + 1.574*omega[i] - 0.176*pow(omega[i], 2));}
	  int task;  //0=find P   1=find T   2=find v   3= find K   4=find Zv
	  secant<thermolib> *solver;

   public:
      double P();		//retruns pressure at T and v, in atm
      double T();		//returns temperature at P and v, in K
      double v();		//returns volume flow at P, T, n(), in m3/s
	  double K() {Z = Zv(); return phiL(0)/phiV(0);}  //returns the vapor-liquid equilibirum constant
      double K(int i) {return phiL(i)/phiV(i);}
      double compres_coeff(){return 1.0;};


// affectation :
  thermolib & operator = ( const thermolib & t );

  thermolib ( int d = 1 ) { dim=d; construct();}

	 void send(double pc, double tc, double w) { Pc[0] = pc*101.325; Tc[0] = tc; omega[0]=w;}
	 void send(double*, double*, double*, double*);
	 void set(double p, double t, double v, double n) {pressure=p*101.325; temperature=t; molevolume=0.001*n/v; mole=n;}
	  double f(double);
	  int get_dim() {return dim;}
      ~thermolib();
	  void reset(int);
};
#endif

