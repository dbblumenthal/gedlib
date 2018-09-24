clear all;

fid = fopen('output.txt','r');
X   = fscanf(fid,'%f%f');
fclose(fid);
X = reshape(X,3,length(X)/3);
X = X';

n = length(X);
A = [];
B = [];

for i = 1:n
    Sa=size(A);
    Sb=size(B);
    if ( X(i,3) > 0.5 )
        j = Sa(1,1);
        A(j+1,1) = X(i,1);
        A(j+1,2) = X(i,2);
    else
        j = Sb(1,1);
        B(j+1,1) = X(i,1);
        B(j+1,2) = X(i,2);
    end
end

hold on;

if ( length(B) > 0 )
    plot( B(:,1) , B(:,2) , 'r.' );
end
if ( length(A) > 0 )
    plot( A(:,1) , A(:,2) , 'b.' );
end

%title ('Title','fontsize',15);
xlabel('constraint j','fontsize',15);
ylabel('objective','fontsize',15);
%axis([0 2000 -0.7 0.02]);

hold off;