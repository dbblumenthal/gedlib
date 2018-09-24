% -----------------------------------------------------------
% file: test_hungarian_lsap_init.m
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
disp('Compare initializations for LSAP');
disp(' ');

% ------------------------------------------
% with vs without initialization with random costs
disp('with vs without initialization with random costs');

bg = 100;
stp = 50;
nd = 1000;
k = 1;
nb = (nd-bg)/stp;
tpsW = zeros(1,nb); tpsWO = zeros(1,nb);

for n = bg:stp:nd

    C = randiLSAPCosts(n,n,10,'int32');
    
    % with Hungarian algorithm
    tic;
    phi = hungarianLSAP(C,0);
    tpsWO(k) = toc;
    
    tic;
    phi = hungarianLSAP(C);
    tpsW(k) = toc;
    k = k+1;
    
end

figure;
plot(bg:stp:nd,tpsWO,'-r');
hold on;
plot(bg:stp:nd,tpsW,'--b');
title('With vs without initialization: time');
xlabel('number n');
ylabel('time (in s)');
legend('without initialization','with initialization','Location','northwest');
