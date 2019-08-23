/** -----------------------------------------------------------
  MEX file for computing an approximate solution to the LSAP with greedy algorithms
 
    [rho,varrho,cost] = lsapeGreedy(C,greedy_type)

  Given a nxm cost matrix C (integer ou floatting values), 
  with last row and last column corresponding to null elements,
  it computes an assignment with low cost

  rho is the assignment from the rows to the columns ((n-1)x1 matrix)
  varrho is the assignment from the columns to the rows (1x(m-1) matrix)
  cost is the cost of the assignment
  
  optional integer greedy_type:
  0: Basic
  1: Refined
  2: Loss (default)
  3: Basic sort
  4: Counting sort (integers only)

  supporting cost values: int16, int32, int64, single, double (default)

  This is a MEX-file for MATLAB.
  
  This file is part of LSAPE.
  LSAPE is free software: you can redistribute it and/or modify it
  under the terms of the CeCILL-C License. See README for more details.

     Copyright 2015-2017
      authors: Sebastien Bougleux
  institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
   last modif: July 6 2017
 
  execute matlab file 'compile_mex' to compile this function and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <lsape.h>

using namespace lsape;

template <typename DT>
void greedy_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  const size_t *nrc = mxGetDimensions(prhs[0]);
  const unsigned int nrows = (int)nrc[0];
  const unsigned int ncols = (int)nrc[1];
  const size_t dimr2[2] = { nrc[0]-1, 1 }, dimc2[2] = { 1, nrc[1]-1 };
  GREEDY_METHOD greedy_type = BASIC;
  if (nrhs > 1)
  {
    greedy_type = (GREEDY_METHOD)mxGetScalar(prhs[1]);
  }
  //------------------------------------------------------------------
  plhs[0] = mxCreateNumericArray(2, dimr2, mxUINT32_CLASS, mxREAL);
  unsigned int *rho = (unsigned int*)mxGetData(plhs[0]);
  plhs[1] = mxCreateNumericArray(2, dimc2, mxUINT32_CLASS, mxREAL);
  unsigned int *varrho = (unsigned int*)mxGetData(plhs[1]);
  DT cost = lsapeGreedy<DT>(C,nrows,ncols,rho,varrho,greedy_type);
  plhs[2] = mxCreateDoubleScalar((double)cost);
  for (unsigned int i = 0; i < nrows-1; i++) rho[i]++;
  for (unsigned int j = 0; j < ncols-1; j++) varrho[j]++;
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  if (nrhs < 1 || nrhs > 2 || nlhs != 3) mexErrMsgTxt("USAGE: [rho,varrho,cost] = lsapeGreedy(C,greedy_type)\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
  case mxINT16_CLASS: greedy_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxINT32_CLASS: greedy_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxINT64_CLASS: greedy_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxSINGLE_CLASS: greedy_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
  default: greedy_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
