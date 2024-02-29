function x=hwk3p3(x0,niter,eps)
% Newton's method with Hessian modification
format long
x = x0.';
[f,g,H] = rosenbrocknfgH(x);
%[f,g,H]=fentonfgH(x);
[m,n] = size(H);
alphabar = 1;
rho = 0.5;
c1 = 0.5;
N = 1;
beta = 10^(-3);
merror = []

while norm(g,2) > eps && N <= niter
    merror(N) = norm(g,2);
    N = N+1;
    try chol(H)
        Bk = H;
    catch ME
        am = min(diag(H));
        if am > 0
            tao = 0;
        else
            tao = -am+beta;
        end
        count = 0;
        err_count = 0;
        while count == err_count
            try chol(H+tao*eye(n))
                Ek = chol(H+tao*eye(n));
            catch ME
                tao = max(2*tao,beta);
                err_count = err_count + 1;
            end
            count = count + 1;
        end
        Bk = H + Ek;
    end

    alpha = alphabar;
    pk = -inv(Bk)*g.';
    [f1,g1,H1] = rosenbrocknfgH(x+alpha*pk);
    % [f1,g1,H1]=fentonfgH(x+alpha*pk);
    while f1 > f + c1*alpha*g*pk
        alpha = rho*alpha;
        f = f1;
        g = g1;
        [f1,g1,H1] = rosenbrocknfgH(x+alpha*pk);
        % [f1,g1,H1]=fentonfgH(x+alpha*pk);
    end
    x = x+alpha*pk;
    % [f,g,H]=fentonfgH(x);
    [f,g,H] = rosenbrocknfgH(x);
end
N
% mesh = linspace(1,N-1,N-1);
% plot(mesh,merror)
% title('L2 norm of gradient for modified method with inital value [10,10,10,10,10,10,10]')



