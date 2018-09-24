% -----------------------------------------------------------
% file: test_lsape_k_solutions.m
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
disp('k solutions to an LSAPE instance');
disp(' ');

% ------------------------------------------
% A simple example with random values
n = 5;
m = 8;
nbEnforcedDelIns = 3;
scaleFactor = 1;

C = randiLSAPECosts(n,m,nbEnforcedDelIns,scaleFactor,'double');
disp('Rectangular cost matrix with random integers');
disp(num2str(C));

% sols is the matrix of solutions (a solution by column)
ksol = 5;
model_type = 1; % Hungarian, model (n+1)x(m+1)
[sols,minCost] = lsapeSolutions(C,ksol,model_type);

disp(['minimal cost = ',num2str(minCost)]);
disp(['The ',num2str(size(sols,2)),' solutions (ask for ',num2str(ksol),')']);
disp(num2str(sols));

% ------------------------------------------
% the same with model (n+m)x(m+n)
model_type = 2; % EBP model

% extend cost matrix
Ce = extendLSAPEinstance(C,max(max(C))+1);

[sols,minCost] = lsapeSolutions(Ce,ksol,model_type,n,m);

disp(['minimal cost = ',num2str(minCost)]);
disp(['The same solutions with extended model']);
disp(num2str(sols));