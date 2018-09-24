%% TEST SURROGATE
clc
x0 = [1;1];
opts = nomadset('display_degree',2,'bb_output_type','OBJ EB','has_sgte',1)
[x,fval] = nomad(@bbsur,x0,[-10;-10],[10;10],opts)
