/** -----------------------------------------------------------
   Test: Hungarian algorithm for solving the Linear Sum Assignment Problem (LSAP)
   authors: Sebastien Bougleux
   institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.
*/

#include <fstream>
#include <iostream>
#include "lsap.h"

int main(int argc, char *argv[]) 
{
  if (argc != 2) { std::cerr << "USAGE: hungarianLSAPEtest filename" << std::endl; return 1; }
  
  //------------------------------------------------------------------
  // input : edit cost matrix as a text file
  std::ifstream ifile(argv[1]);
  if (!ifile) { std::cerr << "Problem with file " << argv[1] << std::endl; return 1; }
  
  // get the edit cost matrix size
  int nrows, ncols;
  ifile >> nrows >> ncols;
  
  // get the edit cost matrix
  double *C = new double[nrows*ncols];
  int i = 0, j;
  for (; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      ifile >> C[j*nrows+i];
  ifile.close();
  //------------------------------------------------------------------
  // outputs
  int *rho = new int[nrows], *varrho = new int[ncols];
  double *u = new double[nrows], *v = new double[ncols];
  hungarianLSAP<double,int>(C,nrows,ncols,rho,u,v,varrho);
  std::cout << "Assignment " << std::endl;
  std::cout << "varphi=[";
  for (i = 0; i < nrows; i++) std::cout << rho[i] << ",";
  std::cout << "]" << std::endl;
  std::cout << "phi=[";
  for (i = 0; i < ncols; i++) std::cout << varrho[i] << ",";
  std::cout << "]" << std::endl;
  std::cout << "Minimal cost=";
  double optC = 0;
  for (i = 0; i < nrows; i++) optC += u[i];
  for (i = 0; i < ncols; i++) optC += v[i];
  std::cout << optC << std::endl;
  delete[] rho; delete[] varrho; delete[] u; delete[] v;
}
