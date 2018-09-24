% -----------------------------------------------------------
% file: compare_lsape_models_hungarian.m
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

tools_directory   = '../../matlab/';
addpath(tools_directory);

% ---------------------------------------------------------- 
% 1st set size loop, 2nd set size is fixed
%%
n = 450;    % can be modified
beg = 50;  % can be modified
step = 100;  % can be modified
endn = 1000; % can be modified
k = 0;
nbtps = 10;
tpsEBP = []; tpsHNG = []; tpsSFBP = []; tpsFBP = []; tpsFBP0 = []; tpsRBPWC = []; tpsRBPCC = [];
mode = 1;
minRBPWC = 0;
%Cl={};

for m = beg:step:endn
    
    switch(mode)
        case 1 % Weiner instance
            C = MacholWien(n,m);
            %C = int64(C);
        case 2 % random integers
            %C = 10*ones(n,m);
            C = randCostLabel(n,m,ceil(sqrt((n-1)*(m-1))/10),5,20,40,20,0);
            [b,T] = triangularIneq(double(C));
            if ~b
                disp(['ok',num2str(sum(sum(T==0)))]);
            end
            C(n,m) = 0;
            C = int32(C);
        case 3 % cir
            C = gallery('cycol',[n,m],ceil(min(n,m)/4));
            C = C - min(min(C));
            %C = C';
        case 4 % integerdata
            C = randCost(n,m,10);
            C(n,m) = inf;
                mnC = min(C);
                mnR = min(C,[],2);
                dv = 1/10;
                for i = 1:n-1
                    C(i,m) = randi(ceil(mnR(i)/dv+1));
                end
                for j = 1:m-1
                    C(n,j) = randi(ceil(mnC(j)/dv+1));
                end
                C(n,m) = 0;
                C = int32(C);
                [b,T] = triangularIneq(double(C));
                if b == 1
                    disp('F*** !');
                    return;
                end
                
         case 5 % geom
            tresh = inf;
            X1 = randi(10*n,2,n-1);
            X2 = randi(10*n,2,m-1);
            C1 = pdist2(X1',X1','minkowski',1);
            C2 = pdist2(X2',X2','minkowski',1);
            [B1,Id1] = sort(C1,2);
            [B2,Id2] = sort(C2,2);
            B1 = B1(:,1:5);
            B2 = B2(:,1:5);
            C1 = mean(B1,2);
            C2 = mean(B2,2);
            C1 = repmat(C1,1,m-1);
            C2 = repmat(C2',n-1,1);
            C = abs(C1-C2);
            mxm = double(max(max(C))/2);
            Crem = mxm*double((1-(B1(:,2)>tresh)));
            Cins = mxm*double((1-(B2(:,2)>tresh)));
            C = [[C,Crem];[Cins',0]];
            
        otherwise
            C = Cl{k+1};
    end
    %k = k + 1;
    %Cl{k} = C;
%end
%return;
    % 
%     if n>m
%         n1 = m;
%         m1 = n;
%     else
%         n1 = n;
%         m1 = m;
%     end
%     C = n1*ones(n1,m1);
%     perm = randperm(m1-1);
%     idxp = randperm(n1-1);
%     per = floor((n1-1)*35/100);
%     idx = sub2ind(size(C),idxp(1:n1-1-per),perm(1:n1-1-per));
%     C(idx) = 1:size(idx,2);
%     C(idxp(n1-per:end),m1) = n1-per;
%     C(n1,perm(n1-per:end)) = n1-per;
%     if n>m
%         C = C';
%     end
    %C = C';
    k = k + 1;
    tps1 = 0; 
    tps1a = 0; tps1b = 0; tpsSFBPi = 0; tps3 = 0; tps4 = 0; tps5 = 0;
    %tps1 = inf; 
    %tps1a = inf; tps1b = inf; 
    %tps2 = inf; tps3 = inf; tps4 = inf; tps5 = inf;
     
    disp(['(n,m)=',num2str(n),',',num2str(m)]);
    
    for tpsk = 1:nbtps
    
        % ----------------------------------------------------------- 
        % Riesen and Bunke: BP with extended cost matrix, i.e. EBP
        tic;
        [rho,varrho,minEBP] = lsapeSolver(C,2);
        
        %tps1 = min(tps1,tps);
        tps1 = tps1 + toc;
        %sum(rho1==m)
%         prc = sum(varrho1==n)*100/m;
%         if tpsk == 1
%             if prc > 0
%                 disp(['prc=',num2str(prc)]);
%             end
%         end
        %disp(['EBP=',num2str(tps)]);
        % ----------------------------------------------------------- 
        % Serratosa 1st proposal: FBP (less general than BP), i.e. FBP0
        tic;
        [rho,varrho,minFBP0] = lsapeSolver(C,5);
        %tps1a = min(tps1a,tps);
        tps1a = tps1a + toc;
         
        % solved as described in our experiments (a bit faster and use less memory)
        tic;
        [rho,varrho,minFBP] = lsapeSolver(C,4);
        %tps1b = min(tps1b,tps);
        tps1b = tps1b + toc;
        
        % ----------------------------------------------------------- 
        % Serratosa 2nd proposal: SFBP (less general than BP)
        tic;
        [rho,varrho,minSFBP] = lsapeSolver(C,6);
        %tps2 = min(tps2,tps);
        tpsSFBPi = tpsSFBPi + toc;
         %disp(['SFBP=',num2str(minSFBP)]);
         %disp(['u=',num2str(u')]);
         %disp(num2str(rho'));
        % ----------------------------------------------------------- 
        % Our proposal: Specialization of Hungarian for LSAPE
        % more general than FBP and SFBP, as general as EBP
        tic;
        [rho,varrho,minHNG] = lsapeSolver(C,1);
        %tps3 = min(tps3,tps);
        tps3 = tps3 + toc;
        % ----------------------------------------------------------- 
        % Our proposal: BP (more general than FBP and SFBP, as general as BP)
        tic
        [rho,varrho,minRBPWC] = lsapeSolver(C,0);
        %tps4 = min(tps4,tps);
        tps4 = tps4 + toc;
        %disp(['EBP=',num2str(tps)]);
        
        % sanity check
        if minEBP ~= minFBP0 || minFBP0 ~= minFBP || minSFBP ~= minFBP || minSFBP ~= minHNG || minHNG ~= minRBPWC
            disp('error');
            disp(['EBP=',num2str(minEBP)]);
            disp(['FBP0=',num2str(minFBP0)]);
            disp(['FBP=',num2str(minFBP)]);
            disp(['SFBP=',num2str(minSFBP)]);
            disp(['HNG=',num2str(minHNG)]);
            disp(['RBPWC=',num2str(minRBPWC)]);
            
            return;
        end

    end
    
        
    tpsEBP = [tpsEBP,tps1/nbtps];
    tpsFBP0 = [tpsFBP0,tps1a/nbtps];
    tpsFBP = [tpsFBP,tps1b/nbtps];
    tpsSFBP = [tpsSFBP,tpsSFBPi/nbtps];
    tpsHNG = [tpsHNG,tps3/nbtps];
    tpsRBPWC = [tpsRBPWC,tps4/nbtps];
    %tpsARBP = [tpsARBP,tps5/nbtps];
    %break;
end

disp('done !');

%%
fctr = 1;

figure1 = figure;
axes1 = axes('Parent',figure1,'FontWeight','bold',...
    'FontSize',14);

box(axes1,'on');
hold(axes1,'all');

plot1 = plot(beg:step:endn,fctr*tpsEBP,'Parent',axes1,'LineWidth',2);
plot2 = plot(beg:step:endn,fctr*tpsSFBP,'Parent',axes1,'LineWidth',2);
plot3 = plot(beg:step:endn,fctr*tpsHNG,'Parent',axes1,'LineWidth',2);
plot4 = plot(beg:step:endn,fctr*tpsFBP0,'Parent',axes1,'LineWidth',2);
plot5 = plot(beg:step:endn,fctr*tpsFBP,'Parent',axes1,'LineWidth',2);
%plot(beg:step:endn,tpsRBP*1000,'--k');
plot6 = plot(beg:step:endn,fctr*tpsRBPWC,'Parent',axes1,'LineWidth',2);

set(plot1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'Marker','o',...
   'Color',[1 0 0],...
   'DisplayName','EBP');
set(plot2,'MarkerFaceColor',[1 1 1],'MarkerSize',7,'Marker','diamond',...
    'Color',[0 0 1],...
    'DisplayName','SFBP');
set(plot3,'MarkerEdgeColor',[0 0.800000011920929 0],'MarkerSize',10,...
    'Marker','x',...
    'Color',[0 0.800000011920929 0],...
    'DisplayName','HNG\epsilon');
set(plot4,'MarkerFaceColor',[0.749019622802734 0 0.749019622802734],...
    'MarkerEdgeColor',[0.749019622802734 0 0.749019622802734],...
    'Marker','square',...
    'Color',[0.749019622802734 0 0.749019622802734],...
    'DisplayName','FBP0');
set(plot5,'MarkerFaceColor',[1 1 1],...
   'MarkerEdgeColor',[0.749019622802734 0.749019622802734 0],...
    'MarkerSize',7,...
    'Marker','>',...
    'Color',[0.749019622802734 0.749019622802734 0],...
    'DisplayName','FBP');
set(plot6,'DisplayName','FLWC','Color',[0 0 0]);

% Create xlabel
xlabel('number m (n is fixed to 1000)','FontWeight','bold','FontSize',14);

% Create ylabel
ylabel('mean time (s)','FontWeight','bold','FontSize',14);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',12);
