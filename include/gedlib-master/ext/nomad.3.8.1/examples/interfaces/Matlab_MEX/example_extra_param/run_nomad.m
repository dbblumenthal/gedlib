%% TEST FOR EXTRA BB PARAMETER
clc
x0 = [1;1];

bb_extra_param = [ 0 ; 1 ];
opts = nomadset('display_degree',2,'bb_output_type','OBJ EB');
[x,fval] = nomad(@bb,x0,[-10;-10],[10;10],opts,bb_extra_param)
