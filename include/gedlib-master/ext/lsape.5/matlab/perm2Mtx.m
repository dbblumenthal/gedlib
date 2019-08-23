function P = perm2Mtx( rho, m )


    n = size(rho,1);
    
    if nargin == 1 % permutation case / perfect matching

        idx = sub2ind([n;n], 1:n, rho');
        P = zeros(n,n);
        P(idx) = 1;
        
    else  % maximum matching case
        
        P = zeros(n,m);
        idx = sub2ind([n;m], 1:n, rho');
        P(idx) = 1;
        
    end

end