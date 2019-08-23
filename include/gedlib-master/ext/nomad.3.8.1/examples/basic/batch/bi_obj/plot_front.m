fid = fopen ( 'front.txt','r');
X=fscanf(fid,'%f%f');
fclose(fid);
X=reshape(X,2,length(X)/2);
X=X';

hold on;

% 1 - x  - x^2 * sin(8 Pi x)
x=(0:0.0001:1);
fx = 1 - x.^2  - x .* sin(8 .* 3.1415926 .* x);
plot( x , fx , 'b-' );


plot( X(:,1) , X(:,2) , 'rx' );

axis([0 1 -1 9]);


hold off;