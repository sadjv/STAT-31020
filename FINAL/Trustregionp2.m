function [xk] = Trustregionp2(gamma,mu,ub,lb,x0,lambda,n,a,b,omega)
% the trust region part for solving nonlinear gradient projection

m = length(x0);
c1 = 0.5;
rho = 0.5;
[f,dL,ddL] = Finalp2Der(x0,mu,lambda,gamma,a,b,n);
% projection P(x-dL)
xk = x0;
x0 = zeros(m,1);
count = 0;
miter = 100;
pk = zeros(9,1);
deltahat = 1;
delta = 0.9;
xita = 1/5;


while norm(xk(1:n)-Projp2(xk(1:n),-dL(1:n),lb,ub)) >= omega && count < miter
    count = count + 1;
    % update ub and lb
    for i=1:m
        ub(i) = min(ub(i),xk(i)+delta);
        lb(i) = max(lb(i),xk(i)-delta);
    end
    xhat = solverp2(gamma,mu,ub,lb,xk,lambda,n,a,b);  % find x^hat
    pk = xhat-xk; % set pk
    xk = x0+pk
    [f1,dL1,ddL1] = Finalp2Der(xk,mu,lambda,gamma,a,b,n);
    rk = (f-f1)/(-dL.'*pk-0.5*pk.'*ddL*pk);
    if rk < 0.25
        delta = 0.25*delta;
    else
        if rk>0.75 && norm(pk)==delta
            delta = min(2*delta,deltahat);
        end
    end
    if rk > xita
        xk = x0 + pk;
    end
    x0 = xk;
    [f,dL,ddL] = Finalp2Der(xk,mu,lambda,gamma,a,b,n);
    
    % projection P(x-dL)
end
if count==1000
    fprintf('maximun numbers of iterations reached.')
end