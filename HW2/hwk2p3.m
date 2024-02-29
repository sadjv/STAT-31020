function x=hwk2p3(x0, c)

ep = 10^(-5);
alphabar = 1;
rho = 0.5;
c1 = 0.5;
Q = [2 0; 0 2*c];
x = x0;
N = 1;
terror = [];
merror = [];

while norm(Q*x,2) > ep && N <= 1000
    merror(N) = norm(Q*x,2);
    terror(N) = norm(x,2);
    N = N+1;
    alpha = alphabar;
    pk = -Q*x/(norm(Q*x,2));
    while 0.5*(x+alpha*pk).'*Q*(x+alpha*pk) > 0.5*x.'*Q*x+c1*alpha*(Q*x).'*pk
        alpha = rho*alpha;
    end
    x = x+alpha*pk;
end

mesh = linspace(1,N-1,N-1);

hold on
set(gca, 'YScale', 'log')
plot(mesh,terror)
plot(mesh,merror)
legend('true error', 'convergence metric')
title(x0)
hold off
