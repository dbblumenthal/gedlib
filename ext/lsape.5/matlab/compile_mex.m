% -----------------------------------------------------------
% file: compile_mex.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux and Évariste Daller
% institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072
% -----------------------------------------------------------
% Execute this file in matlab to compile matlab functions for LSAP and LSAPE
% -----------------------------------------------------------
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README for more
% details.
% -----------------------------------------------------------

disp('Compiling tools for LSAP and LSAPE ...');
files =  { 'hungarianLSAP.cpp', 'greedyLSAP.cpp', 'lsapSolutions.cpp', ...
           'lsapeSolver.cpp', 'lsapeGreedy.cpp', 'lsapeSolutions.cpp', ...
           'lsapeSolverModel.cpp'};
%            'greedykLSAPE.cpp', ...

str = 'mex -I../include ';

for i=1:length(files)
    eval([str files{i} ' ']);
end
