% -----------------------------------------------------------
% file: test_hungarian_lsape.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux
% institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

clear; clc;

% edit cost matrix
n = 7; % nb rows
m = 6; % nb columns

C = [[7 3 8 7 5 7];
     [8 4 8 3 8 3];
     [3 4 3 2 9 1];
     [9 6 7 9 5 1];
     [5 3 6 6 2 2];
     [7 5 6 1 2 3];
     [5 1 2 3 1 0]];
 disp('Cost (last row = insertions, last column = removals) ='); disp(num2str(C));disp(' ');
 % ------------------------------------------

[rho,varrho,u,v] = lsapeSolver(int32(C)); % computing with int32 type

disp(['   rho = ',num2str(rho')]);
disp(['varrho = ',num2str(varrho)]);
disp(['     u = ',num2str(u')]);
disp(['     v = ',num2str(v)]);
disp(['min cost = ',num2str(sum(u)+sum(v))]);
