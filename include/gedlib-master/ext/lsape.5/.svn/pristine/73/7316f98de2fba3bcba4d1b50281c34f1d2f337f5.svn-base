% -----------------------------------------------------------
% file: test_lsap_k_solutions.m
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
disp('k solutions to an LSAP instance');
disp(' ');

% ------------------------------------------
% A simple example
% Cost matrix with random values in {1,...,n}
n = 5;  % can be modified
C = randiLSAPCosts(n,n,1/2,'double');
disp('Square cost matrix with random integers');
disp(num2str(C));

% Hungarian
% rho is the optimal permutation and (u,v) is the dual solution
[sols,minCost] = lsapSolutions(C,10);

disp(['minimal cost = ',num2str(minCost)]);
disp(['The ',num2str(size(sols,2)),' solutions (ask for 10)']);
disp(num2str(sols));
