/*
This unit takes one input stream and divides in two or more
output streams. The pressure and temparature of output streams
are the same as the input's.

Structure in the .process file:
split {name} {index of input stream}  {nb_out} {indexes of output streams and fractions of input}

How to use:
   1- Call the constructor: split1 = new split(nb_out, in, list_out);
   2- Set split fractions: split1->set(fractions);
   3- Set the name: split1->set(name);
   4- Solve: split1->solve();   
*/
#ifndef SPLIT_H
#define SPLIT_H

#include "stream.hpp"
using namespace std;

class split
{
private:

  int i, j;
  bool success;
  double tmp;
  string name;
  int nb_out;			//number of input streams
  stream *in;		//pointer to input stream
  stream **out;		//list of pointers to output streams
  double *frac;			//list of split fractions
	  
public:
  split(int, stream*, stream**);		//defines the connectivities of this unit
  ~split(){}
  void set(double* f)  {frac=f;}
  void set(const string & n) { name = n;}
  bool solve();								//finds the temperature and computes mass balance
  void write();
};

#endif
