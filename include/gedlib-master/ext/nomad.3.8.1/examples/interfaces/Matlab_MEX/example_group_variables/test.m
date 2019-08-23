%% PROBLEM 4 [fval = -2.5]

clc
fun = @(x)   [-x(1) - x(2) - x(3);
              (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.) - 1/4;
                x(1) - x(2);
                x(1) + x(3) + x(4) - 2];          
ub = [1;10;10;5];
lb = [0;0;0;0];
x0 = [0;0;0;0];

opts = []; %not no nomadset just to keep user options
opts.display_degree = 2;
opts.bb_output_type = 'OBJ PB PB PB';
opts.model_search = 'false';
opts.param_file='param.txt';
opts.max_bb_eval = 50;

[xr,fval,ef,iter] = nomad(fun,x0,lb,ub,opts)