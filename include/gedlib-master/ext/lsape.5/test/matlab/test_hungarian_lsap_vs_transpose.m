% -----------------------------------------------------------
% file: test_hungarian_lsap_vs_transpose.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux
% institution: Normandie Univ, CNRS - UNICAEN - ENSICAEN, GREYC UMR 6072 
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

clear; clc;
library_directory = '../../matlab/';
addpath(library_directory);
disp('LSAP example : C vs C^T');
disp(' ');

% ------------------------------------------
% Time comparison between C and C'
% for worst-case time complexity
% Machol-Wien instance c_{ij} = i*j
% YOU HAVE TO WAIT SOME SECONDS
n = 500;
bg = 100;
stp = 100;
nd = 1000;
k = 1;
nb = (nd-bg)/stp;
tpsC = zeros(1,nb); tpsCp = zeros(1,nb);
tpsCr = zeros(1,nb); tpsCpr = zeros(1,nb);

for m = bg:stp:nd

    C = int32(MacholWien(n,m));
    % with Hungarian algorithm
    tic;
    phi = hungarianLSAP(C);
    tpsC(k) = toc;
    
    C = C';
    tic;
    phi = hungarianLSAP(C);
    tpsCp(k) = toc;
    
    C = randiLSAPCosts(n,m,10,'int32');
    
    tic;
    phi = hungarianLSAP(C);
    tpsCr(k) = toc;
    
    C = C';
    tic;
    phi = hungarianLSAP(C);
    tpsCpr(k) = toc;
    
    k = k+1;   
    
end

figure;
plot(bg:stp:nd,tpsC,'-r');
hold on;
plot(bg:stp:nd,tpsCp,'-b');
title('Cost Matrix vs its transpose');
xlabel('number m (n=500)');
ylabel('time (in s)');
legend('C MW','C^T MW','Location','northwest');

figure;
plot(bg:stp:nd,tpsCr,'--r');
hold on;
plot(bg:stp:nd,tpsCpr,'--b');
title('Cost Matrix vs its transpose');
xlabel('number m (n=500)');
ylabel('time (in s)');
legend('C Rand', 'C^T Rand','Location','northwest');

disp('So you should consider |1st set| >= |2nd set| in any case, it is faster');
