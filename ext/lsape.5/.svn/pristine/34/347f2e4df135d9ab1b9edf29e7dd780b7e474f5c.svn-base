% -----------------------------------------------------------
% file: test_hungarian_lsap_float_vs_int.m
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
disp('LSAP example');
disp(' ');

% ------------------------------------------
% Time comparison between double and int32 costs
% for worst-case time complexity
% Machol-Wien instance c_{ij} = i*j
% YOU HAVE TO WAIT SOME SECONDS

disp('double vs int32 with Machol-Wien instance, wait ... !');

bg = 100;
stp = 100;
nd = 1000;
k = 1;
nb = (nd-bg)/stp;
tpsF = zeros(1,nb); tpsI = zeros(1,nb);

for n = bg:stp:nd

    C = MacholWien(n,n);
    
    % with Hungarian algorithm
    tic;
    phi = hungarianLSAP(C);
    tpsF(k) = toc;
    
    C = int32(C);
    tic;
    phi = hungarianLSAP(C);
    tpsI(k) = toc;
    k = k+1;
    
end

figure;
plot(bg:stp:nd,tpsF,'-r');
hold on;
plot(bg:stp:nd,tpsI,'--b');
title('Double vs int32: time');
xlabel('number n');
ylabel('time (in s)');
legend('double','int32','Location','northwest');
