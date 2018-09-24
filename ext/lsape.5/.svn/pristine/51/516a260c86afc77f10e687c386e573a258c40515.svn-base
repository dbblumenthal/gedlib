/** -----------------------------------------------------------
 * Matlab interface for solving the LSAPE
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux and Luc Brun
 * institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072 
 *
 * ----------------------------------------------------------- 
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * ----------------------------------------------------------- 
 *
 *     [X,u,v] = hungarianLSAPE(C,forbassign)
 * [X,mincost] = hungarianLSAPE(C,forbassign)
 *         [X] = hungarianLSAPE(C,forbassign)
 * 
 * Given a cost matrix C, compute a solution X (matrix) to the LSAPE 
 * and a solution (u,v) to its dual problem
 * forbassign is an optional boolean (false by default):
 *   - true -> some value are negative in C (forbidden assignments)
 *   - false -> no forbidden assignments
 *
 * -----------------------------------------------------------
 * execute matlab file 'compile_mex' to compile this function
 * and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsape.hh>
#include <utils.hh>

template <typename DT>
void hungarian_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int nr1 = nrows-1, nc1 = ncols-1;
  int dimr2[2] = { nrows, 1 }, dimc2[2] = { 1, ncols };
  bool forbassign = false;
  if (nrhs == 2) forbassign = (bool)mxGetScalar(prhs[1]);
  //------------------------------------------------------------------
  int dimr[2] = { nr1, 1 }, dimc[2] = { 1, nc1 }, i = 0;
  int dimx[2] = { nrows, ncols };
  plhs[0] = mxCreateNumericArray(2, dimx, mxINT32_CLASS, mxREAL);
  int *X = (int*)mxGetData(plhs[0]);
  int *rho = new int[nr1];
  int *varrho = new int[nc1];
  DT *u = NULL, *v = NULL;
  if (nlhs == 3)
  {
    plhs[1] = mxCreateNumericArray(2, dimr2, mxC, mxREAL);
    u = (DT*)mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(2, dimc2, mxC, mxREAL);
    v = (DT*)mxGetPr(plhs[2]);
  }
  else { u = new DT[nrows]; v = new DT[ncols];}
  hungarianLSAPE<DT,int>(C,nrows,ncols,rho,varrho,u,v,forbassign);
  ecm2Mtx<int,int>(rho,nr1,varrho,nc1,X);
  if (nlhs == 2)
  {
    DT mincost = 0;
    for (int i = 0; i < nr1; i++) { rho[i]++; mincost += u[i]; }
    for (int j = 0; j < nc1; j++) { varrho[j]++; mincost += v[j]; }
    plhs[1] = mxCreateDoubleScalar((double)mincost);
  }
  else
  {
    for (int i = 0; i < nr1; i++) rho[i]++;
    for (int i = 0; i < nc1; i++) varrho[i]++;
  }
  if (nlhs != 3) { delete[] u; delete[] v; }
  delete[] rho; delete[] varrho;
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) 
{
  if (nrhs < 1 || nrhs > 2 || nlhs < 1 || nlhs > 3) mexErrMsgTxt("USAGE: [X] = hungarianLSAPE(C,forbassign)\nUSAGE: [X,mincost] = hungarianLSAPE(C,forbassign)\nUSAGE: [X,u,v] = hungarianLSAPE(C,forbassign)\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
    case mxINT16_CLASS: hungarian_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: hungarian_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: hungarian_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: hungarian_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: hungarian_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}

