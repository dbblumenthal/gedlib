/** -----------------------------------------------------------
 * Matlab interface for solving the LSAP with Hungarian algorithm
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux
 * institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072
 *       Date: 2 Feb. 2017
 *      Modif: March 2018
 * -----------------------------------------------------------
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
 *
 * execute matlab file 'compile_mex' to compile this function
 * and use it in matlab
 *  -----------------------------------------------------------
*/

#include <iostream>
#include <mex.h>
#include <lsap.h>

using namespace lsape;

//------------------------------------------------------------------
template <typename DT>
void allSolLSAP_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  const DT *C = (DT*)mxGetPr(prhs[0]);
  const size_t *dimrc = mxGetDimensions(prhs[0]);
  const unsigned int nrows = (unsigned int)dimrc[0];
  const unsigned int ncols = (unsigned int)dimrc[1];
  const size_t dimr[2] = { dimrc[0], 1 }, dimc[2] = { 1, dimrc[1] };
  int k = -1;
  if (nrhs == 2) k = (int)mxGetScalar(prhs[1]);

  std::list<unsigned int*> solutions;
  DT minCost = lsapSolutions<DT>(C, nrows, ncols, k, solutions);
  
  plhs[0] = mxCreateNumericMatrix(nrows, solutions.size(), mxUINT32_CLASS, mxREAL);
  unsigned int* solMtx = (unsigned int*)mxGetData(plhs[0]);
  lmatchings2mmachings(solutions,solMtx,nrows,1);
  lmatchingsFree(solutions);
  
  if (nlhs > 1) { plhs[1] = mxCreateDoubleScalar(minCost); }
  
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  if (nrhs < 1 || nrhs > 2 || nlhs < 1 || nlhs > 2) mexErrMsgTxt("USE : [sol,minCost] = lsapSolutions(C,ksol)");

  // get the type of cost values
  mxClassID tclass = mxGetClassID(prhs[0]);

  switch(tclass)
  {
    case mxINT16_CLASS: allSolLSAP_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: allSolLSAP_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: allSolLSAP_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: allSolLSAP_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: allSolLSAP_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
