function eval = bb(x,extra_param)

param = [0 ;0];
if (nargin==2)
    param=extra_param;
end

eval=[10*(x(2)-x(1)^2)+param(1); 1 - x(1)+param(2)];
