/** -----------------------------------------------------------
  MEX file for computing an approximate solution to the LSAP with greedy algorithms
 
    [rho,varrho,cost] = greedyLSAPE(C,k,greedy_type)

  Given a nxm cost matrix C (integer ou floatting values),
  with last row and last column representing null elements,
  it computes k error-correcting matchings with low cost

  rho is the assignment from the rows to the columns ((n-1)x1 matrix)
  varrho is the assignment from the columns to the rows (1x(m-1) matrix)
  cost is the cost of the assignment
  
  optional integer greedy_type:
  0: Basic
  1: Refined
  2: Loss (default)
  3: Basic sort (std::sort)
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
#include <greedy-lsape.hh>

template <typename DT>
void greedy_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int kcheap = (int)mxGetScalar(prhs[1]);
  int dimr2[2] = { nrows-1, kcheap }, dimc2[2] = { kcheap, ncols-1 }, dims[2] = {1, kcheap };
  unsigned short init_type = 4;
  if (nrhs > 2)
  {
    init_type = (unsigned short)mxGetScalar(prhs[2]);
  }
  //------------------------------------------------------------------
  plhs[0] = mxCreateNumericArray(2, dimr2, mxINT32_CLASS, mxREAL);
  int *rho = (int*)mxGetData(plhs[0]);
  plhs[1] = mxCreateNumericArray(2, dimc2, mxINT32_CLASS, mxREAL);
  int *varrho = (int*)mxGetData(plhs[1]);
  if (init_type == 4 && (mxC == mxSINGLE_CLASS || mxC == mxDOUBLE_CLASS))
      mexErrMsgTxt("USAGE: integer values only for cost matrix with this greedy type (Counting sort)");
  else
  {
    plhs[2] = mxCreateNumericArray(2, dims, mxC, mxREAL);
    DT *costs = (DT*)mxGetData(plhs[2]);
    greedykLSAPE<DT,int>(C,nrows,ncols,rho,varrho,kcheap,costs,init_type);
    for (int k = 0; k < kcheap; k++)
    {
      for (int i = 0; i < nrows-1; i++) rho[k*(nrows-1)+i]++;
      for (int j = 0; j < ncols-1; j++) varrho[j*kcheap+k]++;
    }
  }
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  if (nrhs < 1 || nrhs > 3 || nlhs != 3) mexErrMsgTxt("USAGE: [rho,varrho,cost] = greedykLSAP(C,k,init_type)\n");
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
