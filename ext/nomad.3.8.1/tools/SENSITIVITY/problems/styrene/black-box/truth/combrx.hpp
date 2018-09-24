#ifndef COMBRX_H
#define COMBRX_H
#include "chemical.hpp"
#include "defines.hpp"
using namespace std;

class combrx {

private : 
  ifstream data;
  bool stop;
  double nO2, nCO2, nH2O;
  char tmp[41];
  string CAS;
  chemical *H2O, *N2, *O2, *CO2, *COMB;
  double LFLo, UFLo, Hro, sum;	  
   
public:
  combrx( const string & cas );
  double O2_flow() { return (O2->M*nO2/1000.0); } //theoritical O2 flow, in kg/mol of COMB
  double N2_flow(){return (0.79*O2_flow()/0.21);} //theoritical N2 flow, in kg/mol of COMB
  double CO2_flow() {return (nCO2*CO2->M/1000.0);} //theoritical CO2 flow, in kg/mol of COMB
  double H2O_flow() {return (nH2O*H2O->M/1000.0);} //theoritical H2O flow, in kg/mol of COMB
  double LFL(double P, double T) {sum=LFLo + 0.03139/Hro*(T-298); if(sum<EPS) return EPS; else return sum;}		//in %vol
  double UFL(double P, double T) {sum=UFLo - 0.03139/Hro*(T-298) + 0.206*(log10(0.101325*P)+1); if(sum>1) return (1-EPS); else return sum;}	//in %vol
  double Hcomb(double T) {return (nCO2*CO2->dH(298,T,1)+ nH2O*H2O->dH(298,T,1)-nO2*O2->dH(298,T,1)-COMB->dH(298,T,1)+Hro);}			//in kJ/mol
  ~combrx(){delete H2O; delete N2; delete O2; delete CO2; delete COMB;}
};
#endif
