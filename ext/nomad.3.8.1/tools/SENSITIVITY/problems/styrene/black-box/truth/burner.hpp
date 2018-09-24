/*
This unit simulates a burner. The user must provide the air excess
and all combustion data as defined in data\\combustion.prop :
      CAS nb_moles_O2 nb_moles_CO2 nb_moles_H2O
Then, mass and energy balances are performed and flows of common
combustion pollutants are estimated.
(reference : Crowl & Louvar)

Structure in the .process file:
burner {name} {index of input stream} {index of output stream} {air excess >0 (ex.: 1.2 is a 120% excess)}

How to use:
   1- Call the constructor: burn = new burner(in, out);			burner(nb_in, chem_list)
   ->set(in, out)
   2- Set the air excess : burn->set(excess);
   3- Set the name of the unit: burn->set(name);
   4- Solve: bool=burn->solve();
*/
#ifndef BURNER_H
#define BURNER_H
#include "stream.hpp"
#include "combrx.hpp"
using namespace std;

class burner
{
private:

  string filename;
  int rem_nb;
  stream *in, *out;
  chemical *O2, *N2, *CO2, *H2O;
  combrx **rx;
  bool *can_burn, stop, OK;
  double eta, NO, NO2, N2O, CO;
  double T, LFLmix, UFLmix, composition;
  string name;
  double * m;
  double a[4], b[4], c[4], K[4];
  int i;
  double buff, Q, m_can_burn, step, num, den;
  ifstream data;
  // ofstream logf, results;
  // terminator *end;
  void fill_K_array();
	  
public:
  // burner(){};
  burner ( int , chemical ** );
  void set ( stream * s1 , stream * s2 ) { in=s1; out=s2; for(i=0;i<in->nb;i++) m[i] = in->chem[i]->m;}
  void set ( const string & n ) { name = n; }
  void set(double e) {eta = e;}
  bool solve(double * y);
  void write();
  void cost();
  double get_cost ( void );
  ~burner();
};
#endif
