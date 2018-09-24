% greedyLSAPE.m Help file for greedyLSAP MEX-file.
%  greedyLSAPE.cpp - Compute an approximate solution to (symmetric or
%  asymmetric) LSAPE with greedy algorithms, i.e. an assignment having a
%  low cost
% 
%    [rho,varrho,cost] = greedyLSAPE(C,greedy_type)
%
%  Given a nxm edit cost matrix C (integer ou floatting values),
%  with last row and last column representing null elements,
%  it computes an error-coorecting matching with low cost
%
%  rho is the assignment from the rows to the columns ((n-1)x1 matrix)
%  varrho is the assignment from the columns to the rows (1x(m-1) matrix)
%  cost is the cost of the assignment
%  
%  optional integer greedy_type:
%  0: Basic
%  1: Refined
%  2: Loss (default)
%  3: Basic sort (C++ std::sort)
%  4: Counting sort (integers only)
% 
%   supporting cost values: int16, int32, int64, single, double (default)
%
% 
%   This is a MEX-file for MATLAB.
%   
%   This file is part of LSAPE.
%   LSAPE is free software: you can redistribute it and/or modify it
%   under the terms of the CeCILL-C License. See README for more details.
% 
%      Copyright 2015-2017
%       authors: Sebastien Bougleux
%   institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
%    last modif: July 6 2017
%  
%   execute matlab file 'compile_mex.m' to compile this function and use it in matlab
