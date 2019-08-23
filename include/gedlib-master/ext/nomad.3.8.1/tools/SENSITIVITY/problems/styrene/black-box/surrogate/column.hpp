/*
This unit simulates a distillation column, using
the FUG method.
(ref : Seader & Henley).

Structure in the .process file:
column {name} {pressure} {index of input stream} {index of bottoms and heads} {indexes of light key and heavy key} {x_LK and x_HK}

How to use:
   1- Call the constructor : col = new column(in, out_B, out_D);		//in is F, out_B is B, out_D is D  column(nb, chem_list)
   set(in, out_B, out_D)
   2- Set operating conditions : col->set(pressure, LK, x_LK, HK, x_HK);		// LK and HK and integer indexes, x_LK is the undesired mole fraction of LK in B, x_HK...
   3- Set the name : col->set(name);
   4- Run the model: col->solve();
*/
#ifndef COLUMN_H
#define COLUMN_H

#include "flash.hpp"
using namespace std;

class column
{
private:
  //  ofstream results, logf;
  bool OK;
  string name;
  stream *F, *B, *D, *L, *V ;
  int LK, HK, feed, i;
  double x_B, x_D, T_b, T_d, T_f, vol, money, diam;
  double Nmin, N, Rmin, Ract,  tmp, Q_condens, Q_reboil;
  double *alpha_1, *alpha_f, *alpha_N, *alpha_m;
  flash *flasher;
  void set_alpha(), first_split(), distribute(), condense(), reboil();
  double Fenske() { return log10(D->chem[LK]->n()*B->chem[HK]->n()/D->chem[HK]->n()/B->chem[LK]->n())/log10(alpha_m[LK]);}
  double Underwood() {return L->n()*(D->chem[LK]->n()/L->chem[LK]->n()-alpha_m[LK]*D->chem[HK]->n()/L->chem[HK]->n())/(D->n()*(alpha_m[LK]-1));}
  double Gilliland(){N=(Ract-Rmin)/(Ract+1); tmp=1-exp((1+54.4*N)*(N-1)/(11+117.2*N)/pow(N, 0.5)); return (tmp+Nmin)/(1-tmp);}
  int Kirkbride() {tmp=pow(B->n()*F->chem[HK]->n()*pow(x_B/x_D,2)/F->chem[LK]->n()/D->n(), 0.206); return int(N/(tmp+1));}
  
public:
  // column(){P=0.0;}
  //	  column(int, chemical*);
  //	  void set(stream*&, stream*&, stream*&);
  column(stream*, stream*, stream*);
  ~column();
  double P;
  void set(double, int, double, int, double);
  void set( const string & n ) { name = n; }
  bool solve();								
  void write();
  void cost(), water(), power();

  double get_cost ( void );
  double get_power ( void ) const { return Q_reboil/0.85-Q_condens; }

  double get_water ( void ) const { return fabs(Q_condens)/(4.185*0.85*0.25*fabs(T_d-298)); }

  int get_N ( void ) const { return (int) N; }

};
#endif
