function C = randiLSAPCosts( n, m, v, classname )
% return a nxm matrix C of random integers in {1,...,v*min(n,m)}
% with v the scaling factir (v=1 by default)
% with classname as a datatype for C ('double' by default)
% classical values for v are [1/10,1,10]
%

    if nargin < 4
        classname = 'double';
    end

    if nargin < 3
        v = 1;
    end

    mx = min(n,m);
    
    mx = int32(v*mx);
    
    if mx < 1
        mx = 1;
    end

    C = randi(mx, [n,m], classname);


end

