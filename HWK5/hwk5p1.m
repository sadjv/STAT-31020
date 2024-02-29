function x=hwk5p1(x_in,tol,iter)

[n,m] = size(x_in);
[f,g,H] = rosenbrocknfgHS(x_in);
x = x_in;
i = 1;
k = 1;
error = [];

while norm(g) > tol && k<=iter
    k = k+1;
    pk = -inv(H)*g;
    
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
   
    x1 = x+alpha*pk;
    s = x1-x;
    y = g1-g;
    if s.'*y >= 0.2*s.'*H*s
        theta = 1;
    else
        theta = 0.8*(s.'*H*s)/(s.'*H*s-s.'*y);
    end
    r = theta*y+(1-theta)*H*s;
    H = H - (H*s*s.'*H)/(s.'*H*s)+r*r.'/(r.'*s);
    error(k) = norm(g1)/norm(g);
    f = f1;
    g = g1;
    x = x1;
end
error