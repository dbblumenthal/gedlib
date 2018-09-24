/** -----------------------------------------------------------
 * Matlab interface for solving the LSAP with Hungarian algorithm
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux and Luc Brun
 * institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072 
 *       date: Dec 1 2015
 * last modif: Jan 30 2017
 * ----------------------------------------------------------- 
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * ----------------------------------------------------------- 
 * 
 *         [X] = hungarianRect(C,forb)
 *     [X,u,v] = hungarianLSAP(C,forb)
 * [X,mincost] = hungarianRect(C,forb)
 *      
 * Given a cost matrix C (integer ou floatting values), compute a solution X to the LSAP 
 * and a solution (u,v) to its dual problem
 * optional boolean parameter forb:
 *    - true  -> forbidden assignments are represented by negative cost values
 *    - false -> no forbidden assignments (by default)
 *
 * -----------------------------------------------------------
 *
 * execute matlab file 'compile_mex' to compile this function
 * and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsap.hh>
#include <utils.hh>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  if (nrhs < 1 || nrhs > 2 || nlhs < 1 || nlhs > 3)
    mexErrMsgTxt("Hungarian algorithm for LSAP with square or rectangular cost matrix C\nUSE: [rho] = hungarianLSAP(C,forb)\nUSE: [rho,u,v] = hungarianLSAP(C,forb)\nUSE: [rho,mincost] = hungarianLSAP(C,forb)\n forb is optional, true if C contains negative values (forbidden assignments), false else");
  //------------------------------------------------------------------
  // 1st input : edit cost matrix
  // get dimensions
  mxClassID tclass = mxGetClassID(prhs[0]);
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int dimr[2] = { nrows, 1 }, dimc[2] = { 1, ncols };
  bool forb_assign = false;
  if (nrhs == 2) forb_assign = (bool)mxGetScalar(prhs[1]);
  //------------------------------------------------------------------
  //-- 1st output: assignment function -------------------------------
  plhs[0] = mxCreateNumericArray(2, dimr, mxINT32_CLASS, mxREAL);
  int *rho = (int*)mxGetData(plhs[0]);
  //------------------------------------------------------------------
  switch(tclass)
  {
  case mxINT16_CLASS:
    {
      short *u = NULL, *v = NULL;
      if (nlhs == 3)
      {
	plhs[1] = mxCreateNumericArray(2, dimr, tclass, mxREAL);
	u = (short*)mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(2, dimc, tclass, mxREAL);
	v = (short*)mxGetPr(plhs[2]);
      }
      else { u = new short[nrows]; v = new short[ncols]; }
      const short *C = (short*)mxGetPr(prhs[0]);
      hungarianLSAP<short,int>(C,nrows,ncols,rho,u,v,forb_assign);
      //------------------------------------------------------------------
      int dim2[2] = { nrows, ncols };
      plhs[0] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
      double *X = (double*)mxGetData(plhs[0]);
      pperm2Mtx(rho,nrows,ncols,X);
      if (nlhs == 2)
      {
	short mincost = 0;
	for (int i = 0; i < nrows; i++) { rho[i]++; mincost += u[i]; }
	for (int j = 0; j < ncols; j++) { mincost += v[j]; }
	plhs[1] = mxCreateDoubleScalar((double)mincost);
      }
      else for (int i = 0; i < nrows; i++) rho[i]++;
      if (nlhs != 3) { delete[] u; delete[] v; }
    }
    break;
  case mxINT32_CLASS:
    {
      int *u = NULL, *v = NULL;
      if (nlhs == 3)
      {
	plhs[1] = mxCreateNumericArray(2, dimr, tclass, mxREAL);
	u = (int*)mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(2, dimc, tclass, mxREAL);
	v = (int*)mxGetPr(plhs[2]);
      }
      else { u = new int[nrows]; v = new int[ncols]; }
      const int *C = (int*)mxGetPr(prhs[0]);
      hungarianLSAP<int,int>(C,nrows,ncols,rho,u,v,forb_assign);
      //------------------------------------------------------------------
      int dim2[2] = { nrows, ncols };
      plhs[0] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
      double *X = (double*)mxGetData(plhs[0]);
      pperm2Mtx(rho,nrows,ncols,X);
      if (nlhs == 2)
      {
	int mincost = 0;
	for (int i = 0; i < nrows; i++) { rho[i]++; mincost += u[i]; }
	for (int j = 0; j < ncols; j++) { mincost += v[j]; }
	plhs[1] = mxCreateDoubleScalar((double)mincost);
      }
      else for (int i = 0; i < nrows; i++) rho[i]++;
      if (nlhs != 3) { delete[] u; delete[] v; }
    }
    break;
  case mxINT64_CLASS:
    {
      int64_T *u = NULL, *v = NULL;
      if (nlhs == 3)
      {
	plhs[1] = mxCreateNumericArray(2, dimr, tclass, mxREAL);
	u = (int64_T*)mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(2, dimc, tclass, mxREAL);
	v = (int64_T*)mxGetPr(plhs[2]);
      }
      else { u = new int64_T[nrows]; v = new int64_T[ncols]; }
      const int64_T *C = (int64_T*)mxGetPr(prhs[0]);
      hungarianLSAP<int64_T,int>(C,nrows,ncols,rho,u,v,forb_assign);
      //------------------------------------------------------------------
      int dim2[2] = { nrows, ncols };
      plhs[0] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
      double *X = (double*)mxGetData(plhs[0]);
      pperm2Mtx(rho,nrows,ncols,X);
      if (nlhs == 2)
      {
	int64_T mincost = 0;
	for (int i = 0; i < nrows; i++) { rho[i]++; mincost += u[i]; }
	for (int j = 0; j < ncols; j++) { mincost += v[j]; }
	plhs[1] = mxCreateDoubleScalar((double)mincost);
      }
      else for (int i = 0; i < nrows; i++) rho[i]++;
      if (nlhs != 3) { delete[] u; delete[] v; }
    }
    break;
  case mxSINGLE_CLASS:
    {
      float *u = NULL, *v = NULL;
      if (nlhs == 3)
      {
	plhs[1] = mxCreateNumericArray(2, dimr, tclass, mxREAL);
	u = (float*)mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(2, dimc, tclass, mxREAL);
	v = (float*)mxGetPr(plhs[2]);
      }
      else { u = new float[nrows]; v = new float[ncols]; }
      const float *C = (float*)mxGetPr(prhs[0]);
      hungarianLSAP<float,int>(C,nrows,ncols,rho,u,v,forb_assign);
      //------------------------------------------------------------------
      int dim2[2] = { nrows, ncols };
      plhs[0] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
      double *X = (double*)mxGetData(plhs[0]);
      pperm2Mtx(rho,nrows,ncols,X);
      if (nlhs == 2)
      {
	float mincost = 0;
	for (int i = 0; i < nrows; i++) { rho[i]++; mincost += u[i]; }
	for (int j = 0; j < ncols; j++) { mincost += v[j]; }
	plhs[1] = mxCreateDoubleScalar((double)mincost);
      }
      else for (int i = 0; i < nrows; i++) rho[i]++;
      if (nlhs != 3) { delete[] u; delete[] v; }
    }
    break;
  default:
    {
      double *u = NULL, *v = NULL;
      if (nlhs == 3)
      {
	plhs[1] = mxCreateNumericArray(2, dimr, tclass, mxREAL);
	u = mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(2, dimc, tclass, mxREAL);
	v = mxGetPr(plhs[2]);
      }
      else { u = new double[nrows]; v = new double[ncols]; }
      const double *C = mxGetPr(prhs[0]);
      hungarianLSAP<double,int>(C,nrows,ncols,rho,u,v,forb_assign);
      //------------------------------------------------------------------
      int dim2[2] = { nrows, ncols };
      plhs[0] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
      double *X = (double*)mxGetData(plhs[0]);
      pperm2Mtx(rho,nrows,ncols,X);
      if (nlhs == 2)
      {
	double mincost = 0;
	for (int i = 0; i < nrows; i++) { rho[i]++; mincost += u[i]; }
	for (int j = 0; j < ncols; j++) { mincost += v[j]; }
	plhs[1] = mxCreateDoubleScalar((double)mincost);
      }
      else for (int i = 0; i < nrows; i++) rho[i]++;
      if (nlhs != 3) { delete[] u; delete[] v; }
    }
  }
}
