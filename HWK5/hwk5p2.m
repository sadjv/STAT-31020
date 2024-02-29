function x=hwk5p2(x_in,tol,iter)

[n,m] = size(x_in);
[f,g,H] = rosenbrocknfgHS(x_in);
x = x_in;
i = 1;
k = 1;
error = [];


while norm(g) > tol && k<=iter
    k = k+1;
    ek = min(0.5, sqrt(norm(g)));
    z1 = 0;
    r1 = g;
    d1 = -r1;
    if d1.'*H*d1 <=0
        pk = -g;
    else
        while d1.'*H*d1 > 0
            alpha = r1.'*r1/(d1.'*H*d1);
            z2 = z1 + alpha*d1;
            r2 = r1 + alpha*H*d1;
%             if norm(r2) < ek
%                 pk = z2;
%                 break
%             end
%
            beta = r2.'*r2/(r1.'*r1);
            d1 = -r2 + beta*d1;
            r1 = r2;
            z1 = z2;
        end
        pk = z1;   
    end
    
    % backtracking
    alpha = 1;
    rho = 0.8;
    c1 = 0.5;
    [f1,g1,H1] = rosenbrocknfgHS(x+alpha*pk);
    i = i+1;
    while f1 > f+c1*alpha*g.'*pk
        alpha = rho*alpha;
        [f1,g1,H1] = rosenbrocknfgHS(x+alpha*pk);
        i = i+1;
    end
    x = x + alpha*pk;
    [f1,g1,H1] = rosenbrocknfgHS(x);
    error(k) = norm(g1)/norm(g);
    f = f1;
    g = g1;
    H = H1;
    i = i+1;
end

error