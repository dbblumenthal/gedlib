% lsapSolutions.m
%  Compute at most k solutions to an LSAP instance
% 
%    [solutions,minCost] = lsapSolutions(C,ksol)
%
%  Compute at most ksol solutions to the LSAP instance C of size (nrows x
%  ncols)
%
%  solutions is a (nrows x kcomp) matrix, each row represents an assignment
%  minCost is the optimal cost (optional), kcom<=ksol is the number of
%  solutions enumerated by the algorithm
%
%  support cost values: int16, int32, int64, single, double (default)
%
% -------------------------------------------------------------------------
% lsapSolutions.m Help file for lsapSolutions.cpp MEX-file.
% ------------------------------------------------------------------------- 
%   This file is part of LSAPE.
%   LSAPE is free software: you can redistribute it and/or modify it
%   under the terms of the CeCILL-C License. See README for more details.
% 
%      Copyright 2015
%       authors: Sebastien Bougleux
%   institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, France
%    last modif: March 2018
%  
%   execute matlab file 'compile_mex.m' to compile this function and use it in matlab
