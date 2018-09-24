/*
This unit performs the heating of a stream (mode 0) or the heat 
exchange to a stream (mode 1).
(ref : McCabe, Smith & Harriott)

Structure in the .process file:
mode 0 : heatx {name} 0 {index of input stream} {index of output stream} {T out} {efficiency}						//efficiency is a fraction between 0 and 1
mode 1 : heatx {name} 1 {index of input stream} {index of output stream} {Q} {efficiency}							//Q is the heat flow in kW
How to use:
   1- Call the constructor: heat = new heatx(mode, in, out);
   2- Set the operating conditions : heat->set(T_out, eta);	//mode 0
                                           or : heat->set(Q, eta);				//mode 1
   3- Set the name of the unit: heat->set(name);
   4- Solve: bool=heat->solve();
*/
#ifndef HEATX_H
#define HEATX_H

#include "stream.hpp"
#include "bissection.hpp"
using namespace std;

class heatx
{
private:
  int i;
  bool success, mode;
  bissection<heatx> *solver;
  // ofstream logf, results;
  double min, max;
  string name;
  stream *in, *out;	//streams of the unit
  double eta, Q, Qreal, T;
  
public:
  heatx(){}
  heatx(bool, stream*, stream*);		//defines the connectivities of this unit
  ~heatx(){delete solver;}
  void set(double, double);
  void set(const string & n) { name = n; }
  double f(double);							
  bool solve();								
  void write();
  void power(), water(), cost();
  double get_cost();
  double get_power();
  double get_water();
};
#endif
