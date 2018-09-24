function [Cext,val] = extendLSAPEinstance(C,val)
% Extend a cost matrix for error-correcting matching to its square version
%
%  [Cext,val] = extendCostMtx(C,val)
%
%  inputs:
%  C (n+1)x(m+1) cost matrix of an error-correcting bipartite graph
%  val the value used to fill 2nd and 3rd block of Cext (optional)
%  val = max(max(C)) + 1 if not given
%
%  outputs:
%  Cext (n+m)x(m+n) extended cost matrix
%  valcpt (optional) is the value val used to fill 2nd and 3rd block of
%  Cext
%  
%         | c(1,1)   ...  c(1,m-1)  | c(1,m) val   ...      val   |
%         |   .              .      |  val  c(2,m) val ...   .    | 
%         |   .              .      |   .                   val   |
%         | c(n-1,1) ... c(n-1,m-1) |  val     ...   val c(n-1,m) |
%  Cext =  -------------------------------------------------------
%         | c(n,1) val  ...    val  |                             |
%         |  val c(n,2) val ..  .   |           0_{m,n}           |
%         |   .                val  |                             |
%         |  val  ...  val c(n,m-1) |                             |
%
%
%   author: Sebastien Bougleux
%   institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC
% -----------------------------------------------------------
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% -----------------------------------------------------------
    
    if ~exist('val','var') || nargin == 1
        val = max(max(C)) + 1;
    end
    [n,m] = size(C);
    n1 = n-1;
    m1 = m-1;
    D = val .* ones(n1,n1);
    E = val .* ones(m1,m1);
    
    D(1:n:end) = C(n*m-n+1:n*m-1);
    E(1:m:end) = C(n:n:n*m-2);
    Cext = [[C(1:n1,1:m1),D];[E,zeros(m1,n1)]];

end
