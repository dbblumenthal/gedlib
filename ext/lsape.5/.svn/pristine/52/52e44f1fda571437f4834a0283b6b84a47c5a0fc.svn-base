%SinkhornKnopp algorithm with regularization (balancing)
%
%  [ X,nbit ] = sinkhornKnopp( G, lambda, tol )
%
%  G must be non-negative and square
%  
%  lambda (optional) is the regularization parameter in [0,1], 0 corresponding to no
%  regularization (basic Sinkhorn) by default
%
%  tol (optional) is the error parameter to control the convergence
%  set to 0.01/n by default
%
% -----------------------------------------------------------
% authors: Sebastien Bougleux
% institution: Normandie Univ, UNICAEN - ENSICAEN - CNRS, GREYC, Caen, France 
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 
function [ X,nbit ] = sinkhornKnopp( G, lambda, tol )

    [n, m] = size(G);
    
    if n~=m
        error('n must b equal to m');
    end
    
    if nargin < 3
        tol = 0.01/n;
    end
    
    if nargin < 2
        lambda = 0;
    end 
    
    r = ones(n,1); 
    c = r;
    d = G'*r + lambda*sum(r);
    nbit = 0;
    
    while norm(c.*d - 1,1) > tol
        c = 1./d;
        r = 1./(G*c+ lambda*sum(c));
        d = G'*r + lambda*sum(r);
        nbit = nbit+1;
    end
    
     X = diag(r)*G*diag(c);
    
end
