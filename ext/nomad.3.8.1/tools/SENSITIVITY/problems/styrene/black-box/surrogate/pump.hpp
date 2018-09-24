/*
This unit takes one input stream and divides in two or more
output streams. The pressure and temparature of output streams
are the same as the input's.
(ref : McCabe, Smith & Harriott)

Structure in the .process file:
pump {name} {index of input stream} {index of output stream} {output pressure in atm} {efficiency, between 0 and 1}

How to use:
   1- Call the constructor: pump1 = new pump(in, out);
   2- Set conditions: pump1->set(P_output, efficiency);
   3- Set the name: pump1->set(name);
   4- Solve: pump1->solve();   
*/
#ifndef PUMP_H
#define PUMP_H

#include "stream.hpp"
using namespace std;

class pump
{
private:
  int i, j, n;
  double state, tmp, tmp1;
  bool success;
  string name;
  stream *in;		//pointer to input stream
  stream *out;		//pointer to output stream
  
public:
  double P, W, eta;			//output presure in atm, work in kW and efficiency
  pump(stream* s1, stream* s2) {in=s1; out=s2; success=true; W=0.0; n=0; tmp=0.0;}
  ~pump(){}
  void set(double p, double e){P = p; eta = e; state=in->quality();}
  void set(const string & n) { name = n; }
  bool solve();				  //finds the temperature and computes mass balance
  void write();
  void cost();
  double get_cost(); // calcule W aussi
  double get_power() const { return W; }

  void  power();
};
#endif
