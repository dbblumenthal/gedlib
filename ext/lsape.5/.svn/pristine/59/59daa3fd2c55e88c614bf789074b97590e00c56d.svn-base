% -----------------------------------------------------------
% file: test_hungarian_vs_greedy_lsap.m
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
disp('Greedy LSAP example');
disp(' ');

% ------------------------------------------
% optimal vs greedy
disp(' ');
disp('Optimal vs greedy with random instance, wait ...!');

bg = 100;
stp = 100;
nd = 2000;
k = 1;
nb = (nd-bg)/stp;

tpsG1 = zeros(1,nb); tpsWO = zeros(1,nb);
solG1 = zeros(1,nb); solWO = zeros(1,nb);
tpsG2 = zeros(1,nb); solG2 = zeros(1,nb);
tpsG3 = zeros(1,nb); solG3 = zeros(1,nb);
tpsG4 = zeros(1,nb); solG4 = zeros(1,nb);
tpsG5 = zeros(1,nb); solG5 = zeros(1,nb);

k = 1;

for n = bg:stp:nd

    C = randiLSAPCosts(n,n,1,'int32');
    
    % with Hungarian algorithm
    tic;
    [phi,u,v] = hungarianLSAP(C);
    tpsWO(k) = toc;
    solW0(k) = sum(u)+sum(v);
    
    tic;
    [phi,solG1(k)] = greedyLSAP(C,0); % Basic greedy algorithm
    tpsG1(k) = toc;
    
    tic;
    [phi,solG2(k)] = greedyLSAP(C,1); % Refined greedy algorithm
    tpsG2(k) = toc;
    
    tic;
    [phi,solG3(k)] = greedyLSAP(C,2); % Loss greedy algorithm
    tpsG3(k) = toc;
    
    tic;
    [phi,solG4(k)] = greedyLSAP(C,3); % Basic sort greedy algorithm
    tpsG4(k) = toc;
    
    tic;
    [phi,solG5(k)] = greedyLSAP(C,4); % Counting sort greedy algorithm
    tpsG5(k) = toc;
    
    k = k+1;
    
end

figure;
plot(bg:stp:nd,tpsWO,'-r');
hold on;
plot(bg:stp:nd,tpsG1,'--b');
plot(bg:stp:nd,tpsG2,'--k');
plot(bg:stp:nd,tpsG3,'--g');
plot(bg:stp:nd,tpsG4,'-p');
plot(bg:stp:nd,tpsG5,'--p');
title('Hungarian vs Greedy LSAP: time');
xlabel('number n');
ylabel('time (in s)');
legend('Hungarian LSAP','Basic Greedy LSAP','Refined Greedy LSAP','Loss Greedy LSAP','Basic sort Greedy LSAP','Counting sort Greedy LSAP','Location','northwest');

figure;
plot(bg:stp:nd,solW0,'-r');
hold on;
plot(bg:stp:nd,solG1,'--b');
plot(bg:stp:nd,solG2,'--k');
plot(bg:stp:nd,solG3,'--g');
plot(bg:stp:nd,solG4,'-p');
plot(bg:stp:nd,solG5,'--p');
title('Hungarian vs Greedy LSAP: costs');
xlabel('number n');
ylabel('cost');
legend('Hungarian LSAP','Basic Greedy LSAP','Refined Greedy LSAP','Loss Greedy LSAP','Basic sort Greedy LSAP','Counting sort Greedy LSAP','Location','northwest');


% ------------------------------------------
% Optimal vs approximation
disp(' ');
disp('Optimal vs greedy with random instance, unbalanced case, wait ...!');

k = 1;
n=2000;
nd = 2*n;
bg = 100;
stp = 100;
nb = (nd-bg)/stp;

tpsG1 = zeros(1,nb); tpsC = zeros(1,nb);
solG1 = zeros(1,nb); solC = zeros(1,nb);
tpsG2 = zeros(1,nb); solG2 = zeros(1,nb);
tpsG3 = zeros(1,nb); solG3 = zeros(1,nb);
tpsG4 = zeros(1,nb); solG4 = zeros(1,nb);
tpsG5 = zeros(1,nb); solG5 = zeros(1,nb);

for m = bg:stp:nd

    C = randiLSAPCosts(n,m,1,'int32');
    
    % with Hungarian algorithm
    tic;
    [phi,u,v] = hungarianLSAP(C);
    tpsC(k) = toc;
    solC(k) = sum(u)+sum(v);
    
    tic;
    [phi,solG1(k)] = greedyLSAP(C,0); % Basic greedy algorithm
    tpsG1(k) = toc;
    
    tic;
    [phi,solG2(k)] = greedyLSAP(C,1); % Refined greedy algorithm
    tpsG2(k) = toc;
    
    tic;
    [phi,solG3(k)] = greedyLSAP(C,2); % Loss greedy algorithm
    tpsG3(k) = toc;
    
    tic;
    [phi,solG4(k)] = greedyLSAP(C,3); % Basic sort greedy algorithm
    tpsG4(k) = toc;
    
    tic;
    [phi,solG5(k)] = greedyLSAP(C,4); % Counting sort greedy algorithm
    tpsG5(k) = toc;
    
    k = k+1;
    
end

figure;
plot(bg:stp:nd,tpsC,'-r');
hold on;
plot(bg:stp:nd,tpsG1,'--b');
plot(bg:stp:nd,tpsG2,'--k');
plot(bg:stp:nd,tpsG3,'--g');
plot(bg:stp:nd,tpsG4,'--y');
plot(bg:stp:nd,tpsG5,'--p');
title('Hungarian vs Greedy unbalanced LSAP: time');
xlabel('number n');
ylabel('time (in s)');
legend('Hungarian LSAP','Basic Greedy LSAP','Refined Greedy LSAP','Loss Greedy LSAP','Basic sort Greedy LSAP','Counting sort Greedy LSAP','Location','northwest');

figure;
plot(bg:stp:nd,solC,'-r');
hold on;
plot(bg:stp:nd,solG1,'--b');
plot(bg:stp:nd,solG2,'--k');
plot(bg:stp:nd,solG3,'--g');
plot(bg:stp:nd,solG4,'--y');
plot(bg:stp:nd,solG5,'--p');
title('Hungarian vs Greedy unbalanced LSAP: costs');
xlabel('number n');
ylabel('cost');
legend('Hungarian LSAP','Basic Greedy LSAP','Refined Greedy LSAP','Loss Greedy LSAP','Basic sort Greedy LSAP','Counting sort Greedy LSAP','Location','northwest');
