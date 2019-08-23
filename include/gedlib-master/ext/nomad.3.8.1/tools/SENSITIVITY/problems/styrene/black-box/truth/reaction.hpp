#ifndef REACTION_H
#define REACTION_H

#include "chemical.hpp"
using namespace std;

class reaction {
private:

  int m;
  double *n, k0, E, Hr, *safe_n, *safe_a;
  // double tmp;
//   char file[41], line[31];
  chemical ** list;
//   ifstream in;
  // ofstream logf;
  // terminator *end;
	
  int find_chemical ( const string & chem_name ) const;

public:
  // reaction(){};
  reaction ( const string & , int ,  chemical ** );
  ~reaction();
  double *a;									//contains the molar coefficients
  double dHr(double);					//returns heat of rection at T, in kJ/mol
  double rate(double, double*);		//returns rate of reaction aT and C[], in mol/s.m3

//   void show_name(){cout<<name;}
};

#endif
