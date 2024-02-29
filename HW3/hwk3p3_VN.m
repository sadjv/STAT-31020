function x=hwk3p3_VN(x0,niter,eps)
% vanilla Newton

x = x0.';
N = 1;
%[f,g,H] = rosenbrocknfgH(x);
[f,g,H]=fentonfgH(x);
errort = [];

while norm(g,2) > eps && N <= niter
    errort(N) = norm(g,2);
    x = x - inv(H)*g.';
    [f,g,H]=fentonfgH(x);
    N = N + 1;
end

mesh = linspace(1,N-1,N-1);
plot(mesh,errort)
title('L2 norm of gradient for VN method with inital value (3,4)')