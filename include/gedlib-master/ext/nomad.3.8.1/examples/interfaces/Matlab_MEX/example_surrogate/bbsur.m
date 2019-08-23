function eval = bbsur(x,sur)

if (nargin==1)
    sur=false;
end

if (sur)
    eval=[9.9*(x(2)-x(1)^2); 1 - x(1)];
else
    eval=[10*(x(2)-x(1)^2); 1 - x(1)];
end
