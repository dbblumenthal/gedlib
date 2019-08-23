/** -----------------------------------------------------------
 * Matlab interface for solving the LSAP with Hungarian algorithm
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux and Luc Brun
 * institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072
 *       Date: 2 Feb. 2017
 *      Modif:
 * -----------------------------------------------------------
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
 *
 *         [rho] = allSolLSAP(C,forb)
 *     [rho,u,v] = allSolLSAP(C,forb)
 * [rho,mincost] = allSolLSAP(C,forb)
 *
 * Given a cost matrix C (integer ou floatting values), compute a solution rho to the LSAP
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

#include <iostream>
#include <hungarian-lsap.hh>
#include <list>
#include <mex.h>
#include <utils.hh>
#include "RandomAllPerfectMatchingsEC.h"


/**
 * \param[in] C    cost matrix
 * \param[in] rows nr of rows in C
 * \param[in] cols nr of cols in C
 * \param[in] n    size of G1
 * \param[in] m    size of G2
 * \param[out] nsol  nr of solutions
 * \param[out] solutions solutions
 * \param[in] k    max number of computed solutions
 */
template <typename DT, typename IT>
void allSolutionsFromOneLSAP(const DT *C, const IT &rows, const IT &cols, const IT &n, const IT &m,
                             IT *rho12, DT *u, DT *v, IT& nsol, IT** solutions, int k)
{
  std::cout << "r = " << rows << std::endl;
  std::cout << "c = " << cols << std::endl;
  std::cout << "n = " << n << std::endl;
  std::cout << "m = " << m << std::endl;
// construct equality directed graph
  // i.e. ({0,...,n-1} U {0,...,m-1}, rho12 U adj21)
  // i.e. find adj21
  cDigraph<IT> edg = equalityDigraph<DT,IT>(C,rows,cols,rho12,u,v);

	RandomAllPerfectMatchingsEC<IT> apm(edg, n, m);
	apm.enumPerfectMatchings(edg, k);
	std::list< IT* > match = apm.getPerfectMatchings();

  nsol = match.size()+1;
  IT* _solutions = new IT[rows*nsol];

  for(int i=0; i<rows; i++){
    _solutions[i] = rho12[i];
  }

  typename std::list< IT* >::iterator it = match.begin();
  for (int s=1; s<nsol; s++) {
    IT* tab = *it;
    for (int i=0; i<rows; i++){
      _solutions[s*rows+i] = tab[i];
    }
    it++;
  }
  *solutions = _solutions;
}

// -----------------------------------------------------------
/**
 * \brief Compute all solutions to the LSAP with Hungarian algorithm given a square cost matrix
 * \param[in] C nxm cost matrix
 * \param[in] n Number of rows of \p C (size of the 1st set)
 * \param[in] m Number of columns of \p C (size of the 2nd set)
 * \param[out] nsol Number of solutions found
 * \param[out] solutions Array of solutions
 * \param[in] forb_assign If true, forbidden assignments are marked with negative values in the cost matrix
 */
template <typename DT, typename IT>
void allSolutionsLSAP(const DT *C, const IT &rows, const IT &cols, const IT &n, const IT &m,
                      IT& nsol, IT** solutions, bool forbassign = false, int k=-1)
{
  IT *rho12 = new IT[rows];
  DT *u = new DT[rows], *v = new DT[cols];
  std::list<std::list<IT>*> *sccL = new std::list<std::list<IT>*>();

  // find a solution
  hungarianLSAP<DT,IT>(C,rows,cols,rho12,u,v,forbassign);

  // find all solutions
  allSolutionsFromOneLSAP<DT,IT>(C,rows,cols, n,m, rho12,u,v, nsol, solutions,k);

  // free memory
  delete[] rho12; delete[] u; delete[] v; delete sccL;
}


//------------------------------------------------------------------
template <typename DT>
void allSolLSAP_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  const DT *C = (DT*)mxGetPr(prhs[0]);
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  
  int n = (int)mxGetScalar(prhs[1]);
  int m = (int)mxGetScalar(prhs[2]);
  
  int dimr[2] = { nrows, 1 }, dimc[2] = { 1, ncols };
  bool forb_assign = false;
  if (nrhs >= 4) forb_assign = (bool)mxGetScalar(prhs[3]);

  int k = -1;
  if (nrhs == 5) k = (int)mxGetScalar(prhs[4]);

  int nsol;
  int* solutions = NULL;

  //if (n+m != nrows || nrows != ncols){
  //  mexErrMsgTxt("C is not (n+m)x(n+m)");
  //  return;
  //}

  allSolutionsLSAP<DT,int>(C, nrows, ncols, n, m, nsol, &solutions, forb_assign, k);

  //*
  plhs[0] = mxCreateNumericMatrix(nrows, nsol, mxINT32_CLASS, mxREAL);

  if (solutions != NULL){
    int* tab = (int*)mxGetData(plhs[0]);
    for (int s=0; s<nsol; s++) {
      for (int i=0; i<nrows; i++){
        tab[s*nrows+i] = solutions[s*nrows+i];
      }
    }

    delete[] solutions;
  }
  //*/
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  if (nrhs < 3 || nrhs > 5)
    mexErrMsgTxt("all solution LSAP : \n USE : rho = allSolLSAP(C, n, m, forb_assign, k)");

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
