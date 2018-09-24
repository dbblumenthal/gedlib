% -----------------------------------------------------------
% file: test_hungarian_lsap.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux
% institution: Normandie Univ, UNICAEN - ENSICAEN - CNRS, GREYC, Caen, France 
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

clear; clc;
library_directory = '../../matlab/';
addpath(library_directory);
disp('Simple LSAP examples');
disp(' ');

% ------------------------------------------
% A simple example
% Cost matrix with random values in {1,...,n}
n = 5;  % can be modified
C = randiLSAPCosts(n,n,1,'int32');
disp('Square cost matrix with random integers');
disp(num2str(C));

% Hungarian
% rho is the optimal permutation and (u,v) is the dual solution
[rho,u,v] = hungarianLSAP(C);
mincost = sum(u)+sum(v);

disp('A solution to the LSAP instance:');
disp([num2str((1:n)'),repmat(' -> ',n,1),num2str(rho)]);
disp('Associated dual solution:');
disp(['u = ',num2str(u')]);
disp(['v = ',num2str(v)]);
disp(['Minimal cost = ',num2str(mincost)]);

P = perm2Mtx(rho);
disp('The solution as a permutation matrix:');
disp(num2str(P));

% ------------------------------------------
% Another simple example
% Rectangular cost matrix with random values in {1,...,100*min{n,m}}
n = 5;  % can be modified
m = 8;
C = randiLSAPCosts(n,m,100,'int32');
disp(' ');
disp('Rectangular cost matrix with random integers');
disp(num2str(C));

% Hungarian
% rho is the optimal permutation and (u,v) is the dual solution
[rho,varrho,u,v] = hungarianLSAP(C);
mincost = sum(u)+sum(v);

disp('A solution to the LSAP instance:');
disp([num2str((1:n)'),repmat(' -> ',n,1),num2str(rho)]);
disp('Associated dual solution:');
disp(['u = ',num2str(u')]);
disp(['v = ',num2str(v)]);
disp(['Minimal cost = ',num2str(mincost)]);

P = perm2Mtx(rho,m);
disp('The solution as a matrix:');
disp(num2str(P));
