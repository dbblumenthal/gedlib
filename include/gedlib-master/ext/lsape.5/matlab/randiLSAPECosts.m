function C = randiLSAPECosts( n,m,pdelins,v,classname)
% Random cost matrix of size (n+1)x(m+1) with at least pdelins removals or
% insertions
% classname is the type for C

    if nargin < 5
        classname = 'double';
    end
    
    if nargin < 4
        v = 1;
    end
    
    if nargin < 3
        pdelins = 0;
    end
    
    C = randiLSAPCosts(n+1,m+1,v,classname);
    
    if pdelins > 0
        Rn = randperm(n);
        Rm = randperm(m);
        if n<=m
            Mn = min(C(Rn(1:pdelins),1:m-1),[],2);
            Mm = min(C(1:n,Rm(1:pdelins+m-n)));
            MMn = min(Mn);
            MMm = min(Mm);
            for i=1:pdelins
                C(Rn(i),m+1) = floor(min(Mn(i),MMm)/2-0.5);
            end
            for i=1:pdelins+m-n
                C(n+1,Rm(i)) = floor(min(Mm(i),MMn)/2-0.5);
            end
        else
            Mn = min(C(Rn(1:pdelins+n-m),1:m-1),[],2);
            Mm = min(C(1:n,Rm(1:pdelins)));
            MMn = min(Mn);
            MMm = min(Mm);
            for i=1:pdelins
                C(n+1,Rm(i)) = floor(min(Mm(i),MMn)/2-0.5);
            end
            for i=pdelins+n-m
                C(Rn(i),m+1) = floor(min(Mn(i),MMm)/2-0.5);
            end
        end
    end
   
    C(end,end) = 0;
    
end

