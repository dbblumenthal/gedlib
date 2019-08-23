#ifndef CHEMICAL_H
#define CHEMICAL_H

#include "thermolib.hpp"
#include <string>
using namespace std;

class chemical
{
 public:
  string name;
  string CAS;
      double M;			//molar weight, in g/mol
      bool state; 		//current state : 0=liquid 1=gas
      double P, T, m, v;	//current values : P in atm, T in K, m in kg/s, v in m3/s

      double n() {return 1000.0*m/M;} //mole flow in mol/s
      double Hvap(double t) {return dHvap*pow((Tc-t)/(Tc-Tb),0.38);} //vaporization heat at specific T (Watson correlation)
      double omega() {return ((-1.0)*log10(Psat(0.7*Tc)/Pc) -1.0);}   //Pitzer acentric factor
      double Tboil(double p) {return (Psat_param[1]/(Psat_param[0]-log(760.01*p))-Psat_param[2]);}  //boiling tempararure at specific P
      double mu(), rho(), Cp(), Cp(bool), Psat(), Psat(double);			  //T-dependant properties
      double K();												//liquid-vapor equilibrium constant
      double gamma(){return Cp(true)/(Cp(true)-8.3144);}//compressibility ratio =Cp/Cv
      double dH(double, double, double);			//enthalpy variation kJ/mol
      void find_T(), find_P(), find_v(), find_state();	//for gases only
      double Tm, Tb, Tc, Pc;											//melting, boiling and critical temp. (K); critical pressure (atm)
      double Ho;//standard formation heat in kJ/mol
      
      //private: 
      int warning, error;
      void check_error();
      double dHvap, tmp;					//vaporization heat (kJ/mol)
      double mu_param[2], Cp_param[4], Cp_liq, Psat_param[3], rho_liq; //correlations parameters
      thermolib *thermo;

      // public:
  // chemical() {};

  // copy-constr. :
  chemical ( const chemical & chem );

  chemical ( const string & chem_name );	//Contructor : initializes fields and reads CAS in the data file
  void purge() {P=T=m=v=0.0; state=false;}
  ~chemical(){delete thermo;};
};

#endif

