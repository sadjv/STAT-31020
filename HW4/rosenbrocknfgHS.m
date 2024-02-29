function [f,g,H] = rosenbrocknfgHS(x)

[n,m] = size(x);


f = 0;
for i=1:n-1
    f = f+100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
end

g1 = 1:n;
sz = size(g1);
g = zeros(sz);

g(1) = -400*x(1)*(x(2)-x(1)^2)+2*(x(1)-1);
g(n) = 200*(x(n)-x(n-1)^2);

for j=1:n-2
    g(j+1) = -400*x(j+1)*(x(j+2)-x(j+1)^2)+2*(x(j+1)-1)+200*(x(j+1)-x(j)^2);
end

H = sparse(n,m);  % use sparse matrix instead of dense matrix for hessian
H(1,1) = -400*x(2)+1200*x(1)^2+2;
H(1,2) = -400*x(1);
H(n,n) = 200;
H(n,n-1) = -400*x(n-1);

for j=1:n-2
    H(j+1,j) = -400*x(j);
    H(j+1,j+1) = -400*x(j+2)+1200*x(j+1)^2+2+200;
    H(j+1,j+2) = -400*x(j+1);
end