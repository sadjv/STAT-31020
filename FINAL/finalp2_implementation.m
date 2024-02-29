% implementation of the function final2
% not very meaningful, since my code is actually not convergent.
% the basic logic of my code is to use final2.m to do the algorithm
% 17.4, update parameters and check for convergence. It will call
% Linesearchp2 (or Trustregionp2) to solve the subproblem
% min grad(L(xi,lambdai,mu)) subject to lb <= x <= ub.
% Other related functions are:
% Finalp2Der.m; calculate Lx, dLx, d^2Lx
% Solverp2: linear gradient projection method for problem 2.

gamma = 0.1;
n = 60;
a = -ones(60,1);
b = ones(30,1);
[t,x,lambda]=finalp2(n,gamma,b,a);

nlist = [10,50,100];
tmean = zeros(3,1);
for i=1:3
    n = nlist(i);
    ttemp = zeros(100,1);
    for j=1:20
        [t,x,lambda]=finalp2(n,gamma,b,a);
        ttemp(j)=t;
    end
    tmean(i)=sum(ttemp)/20;
end

plot(nlist,tmean)
title('mean running time vs n for final2')
xlabel('n')
ylabel('running time in second')
