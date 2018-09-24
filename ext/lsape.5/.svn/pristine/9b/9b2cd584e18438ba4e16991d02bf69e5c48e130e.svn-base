/** -----------------------------------------------------------
   Test: Hungarian algorithm for solving the Linear Sum Assignment Problem with Edition (LSAPE)
   authors: Sebastien Bougleux
   institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.
*/

#include <fstream>
#include <iostream>
#include "hungarian-lsape.hh"

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
  int *rho = new int[nrows-1], *varrho = new int[ncols-1];
  double *u = new double[nrows], *v = new double[ncols];
  hungarianLSAPE<double,int>(C,nrows,ncols,rho,varrho,u,v);
  std::cout << "Assignment with edition" << std::endl;
  int *varphieps = NULL, n = nrows-1, m = ncols-1;
  int nbeps = reconstructInsertions(varrho,n,m,&varphieps);
  std::cout << "varphi=[";
  for (i = 0; i < nrows-1; i++) std::cout << rho[i] << ",";
  std::cout << "{";
  j = 0;
  for (; j < nbeps-1; j++) std::cout << varphieps[j] << ",";
  if (nbeps > 0) std::cout << varphieps[j];
  std::cout << "}]" << std::endl;
  
  delete[] rho; delete[] varrho; delete[] u; delete[] v;
  if (nbeps > 0) delete[] varphieps;  
}
