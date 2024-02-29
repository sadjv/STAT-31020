function [t,x,lambda]=finalp2(n,gamma,b,a)

% initialization
m = n+length(b);
bs = length(b);
xk = rand(m,1);  % initial x0
lambda = 0.1*rand(bs+1,1) % inital 
omestar = 10^-3;
xitastar = 10^-3;
mu = 10;
omega = 1/mu;
xita = 1/(mu)^(0.1);
miter = 100;
[~,dl,~] = Finalp2Der(xk,mu,lambda,gamma,a,b,n);
% create lb and ub
% -1 <= xi <= 1, i=1,...,n. for slack variables si, in this problem,
% we only consider b=1 and 0.5*xk-1^2+xk^2+xk+1^2-b-sk=0, we have
% 0 <= sk <= 1.
ub = ones(m);
lb = -ones(m);
for i=1:bs
    lb(n+i) = lb(n+i) + 1;
end
count = 0;
tic
for i=1:miter
    count = count + 1;
    [xk] = Linesearchp2(gamma,mu,ub,lb,xk,lambda,n,a,b,omega)
    [~,dl,~] = Finalp2Der(xk,mu,lambda,gamma,a,b,n);
    % checking for kkt condition 1: ci(xk)=0, for all i
    done = 1;
    cxklist = zeros(bs+1,1);
    for i=1:bs-1
        cxklist(i) = abs(0.5*xk(2*i-1)^2+xk(2*i)^2+0.5*xk(2*i+1)^2-b(i)-xk(n+i));
        % x(n+i) is slack variable
    end
    cxklist(bs) = abs(0.5*xk(2*bs-1)^2+xk(2*bs)^2-b(bs)-xk(n+bs));
    cxklist(bs+1) = abs(sum(xk(1:n)));
    % equality
    if cxklist(bs+1) > xita
        done = 0;
    end
    % inequality
    for i=1:bs
        if cxklist(i) > xita % x(n+i) is slack variable
            done = 0;
        end
    end
    if done ==1
        finaldone = 1;
        for i=1:bs+1
            if cxklist(i) > xitastar
                finaldone = 0;
            end
        end
        if norm(xk-Projp2(xk,-dl,lb,ub)) < omestar && finaldone == 1
            break
        end
        % update multipliers, tignten tolerances
        lambda = lambda-mu*cxklist;
        xita = xita/(mu)^(0.9);
        omega = omega/mu;
    else
        % increase penalty parameter
        mu = 100*mu;
        xita = 1/(mu)^(0.1);
        omega = 1/(mu);
    end
end
t = toc;
if count==miter
    fprintf('maximun iteration reached')
end
x = xk;
