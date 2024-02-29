function [xk] = Linesearchp2(gamma,mu,ub,lb,x0,lambda,n,a,b,omega)
% the linear search part for solving fk+gradfk*xk+xk*Bk*xk
% (nonlinear gradient projection)use backtracking to find alphak

m = length(x0);
c1 = 0.5;
rho = 0.5;
alphabar = 0.1;
[f,dL,ddL] = Finalp2Der(x0,mu,lambda,gamma,a,b,n);
% projection P(x-dL)
xk = x0;
x0 = zeros(m,1);
count = 0;
miter = 100;
alpha = alphabar;
pk = zeros(9,1);

while norm(xk(1:n)-Projp2(xk(1:n),-dL(1:n),lb,ub)) >= omega && count < miter
%while norm(xk-x0) >= omega && count < miter
    count = count + 1;
    xhat = solverp2(gamma,mu,ub,lb,xk,lambda,n,a,b);  % find x^hat
    pk = xhat-xk; % set pk
    x0 = xk;
    % find a_k using backtracking
    alpha = alphabar;
    [f1,dL1,ddL] = Finalp2Der(xk+alpha*pk,mu,lambda,gamma,a,b,n);
    while f1 > f - c1*alpha*dL.'*pk
        alpha = rho*alpha;
        f = f1;
        dL = dL1;
        [f1,dL1,ddL] = Finalp2Der(xk+alpha*pk,mu,lambda,gamma,a,b,n);
    end
    xk = xk+alpha*pk;
    [f,dL,ddL] = Finalp2Der(xk,mu,lambda,gamma,a,b,n);
    % projection P(x-dL)
end
if count==1000
    fprintf('maximun numbers of iterations reached.')
end