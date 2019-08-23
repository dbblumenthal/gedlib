/** -----------------------------------------------------------
 * Matlab interface for solving the LSAPE
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux
 * institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, Caen, France
 *
 * ----------------------------------------------------------- 
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 *
 * -----------------------------------------------------------
 * execute matlab file 'compile_mex' to compile this function
 * and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <lsape.h>

using namespace lsape;

template <typename DT>
void lsape_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  const size_t *nrc = mxGetDimensions(prhs[0]);
  unsigned int nrows = (unsigned int)nrc[0];
  unsigned int ncols = (unsigned int)nrc[1];
  // 2nd and 3rd inputs original number of elements in each set
  LSAPE_MODEL model_type = ECBP;
  unsigned short init_type = 1;
  if (nrhs >= 2) model_type = (LSAPE_MODEL)mxGetScalar(prhs[1]);
  if (nrhs == 3) init_type = (unsigned short)mxGetScalar(prhs[2]);

  //------------------------------------------------------------------
  // outputs and solve
  const size_t dimr[2] = { nrc[0], 1 }, dimc[2] = { 1, nrc[1] };
  size_t nr1, nc1;
  if (model_type != ECBP)
  {
    plhs[0] = mxCreateNumericArray(2, dimr, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dimc, mxUINT32_CLASS, mxREAL);
    nr1 = nrc[0]; nc1 = nrc[1];
  }
  else
  {
    const size_t dimr2[2] = { nrc[0]-1, 1 }, dimc2[2] = { 1, nrc[1]-1 };
    plhs[0] = mxCreateNumericArray(2, dimr2, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dimc2, mxUINT32_CLASS, mxREAL);
    nr1 = nrc[0]-1; nc1 = nrc[1]-1;
  }

  unsigned int *rho = (unsigned int*)mxGetData(plhs[0]);
  unsigned int *varrho = (unsigned int*)mxGetData(plhs[1]);

  DT *u = NULL, *v = NULL;
  if (nlhs == 4)
  {
    plhs[2] = mxCreateNumericArray(2, dimr, mxC, mxREAL);
    u = (DT*)mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(2, dimc, mxC, mxREAL);
    v = (DT*)mxGetPr(plhs[3]);
  }
  else { u = new DT[nrows]; v = new DT[ncols];}

  lsapeSolverModel<DT>(C,nrows,ncols,rho,varrho,u,v,model_type,init_type);

  if (nlhs == 3)
  {
    DT mincost = 0;
    for (unsigned int i = 0; i < nr1; i++) { rho[i]++; mincost += u[i]; }
    for (unsigned int j = 0; j < nc1; j++) { varrho[j]++; mincost += v[j]; }
    plhs[2] = mxCreateDoubleScalar((double)mincost);
  }
  else
  {
    for (unsigned int i = 0; i < nr1; i++) rho[i]++;
    for (unsigned int j = 0; j < nc1; j++) varrho[j]++;
  }
  if (nlhs != 4) { delete[] u; delete[] v; }
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) 
{
  if (nrhs < 2 || nrhs > 3 || nlhs < 2 || nlhs > 4) mexErrMsgTxt("USAGE: [rho,varrho] = lsapeSolverModel(C,n,m,model_type,init_type)\nUSAGE: [rho,varrho,mincost] = lsapeSolverModel(C,n,m,model_type,init_type)\nUSAGE: [rho,varrho,u,v] = lsapeSolverModel(C,n,m,model_type,init_type)\n See help for more details\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
    case mxINT16_CLASS: lsape_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: lsape_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: lsape_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: lsape_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: lsape_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
