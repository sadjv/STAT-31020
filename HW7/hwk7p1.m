function [x,obj,iter]=hwk7p1(Q,A,b,c,x0)

x = x0;
iter = 1;
obj = 0.5*x.'*Q*x + x.'*c;
[m,n] = size(A);
s = ones(n,1)*2;
lambda = ones(m,1)*2;
rp = A*x - b;
rd = Q*x - A.'*lambda + c - s;
tol = 10^(-5);
s1=1;
s2=1;
objlist = [obj,obj-1];
k = 0;

% & objlist(iter) > objlist(iter+1)
while (s1 >= tol & s2 >= tol & objlist(iter) > objlist(iter+1))
    k = k+1;
    if iter == k
        iter =1;
    else
       iter = iter + 1;
    end
    % find affine delta's
    rxp = diag(x)*diag(s)*ones(n,1);
    xinv = inv(diag(x));
    J = xinv*diag(s) + Q;
    Jinv = inv(J);
    D = chol(A*Jinv*A.');
    dlamaff = (inv(D).'*inv(D))*(-rp+A*Jinv*rd+A*Jinv*xinv*rxp);
    dxaff = Jinv*(A.'*dlamaff-rd-xinv*rxp);
    dsaff = -xinv*rxp-xinv*(diag(s)*dxaff);
    
    mu = x.'*s/m;
    % find alpha^{hat}_aff
%     if x.'*s + dxxaff.'*dsaff < 0
%        alphahatx = -(x.'*s)/(dxaff.'*dsaff);
%     else
%         alphahatx = 1;   
%     end

    alphahats = 1;
    alphahatx = 1;
    for i=1:n
        if s(i)+dsaff(i) < 0
            alphahats = min(alphahats, -s(i)/dsaff(i));
        else
            alphahats = min(alphahats, 1);
        end
    end
    for i=1:n
        if x(i)+dxaff(i) < 0
            alphahatx = min(alphahatx, -x(i)/dxaff(i));
        else
            alphahatx = min(alphahatx, 1);
        end
    end
    
    alphahat = min(alphahatx,alphahats);
    

    
    % calculate u_aff
    muaff = (x+alphahat*dxaff).'*(s+alphahat*dsaff)/m;
    sigma = (muaff/mu)^3;
    
    rxp = diag(x)*diag(s)*ones(n,1) - sigma*mu*ones(n,1);
    dlam = (inv(D).'*inv(D))*(-rp+A*Jinv*rd+A*Jinv*xinv*rxp);
    dx = Jinv*(A.'*dlam-rd-xinv*rxp);
    ds = -xinv*rxp-xinv*(diag(s)*dx);
    
    tao = iter/(iter+1);
    % find alpha
    alphad = 1;
    alphap = 1;
    for i=1:n
        if s(i)+ds(i) < (1-tao)*s(i)
            alphad = min(alphad, -tao*s(i)/ds(i));
        else
            alphad = min(alphad, 1);
        end
    end
    for i=1:n
        if x(i)+dx(i) < (1-tao)*x(i)
            alphap = min(alphap, -tao*x(i)/dx(i));
        else
            alphap = min(alphap, 1);
        end
    end
    
    alpha = min(alphap,alphad);
    
    % set xk+1
    x1 = x;
    x = x + alpha*dx;
    s = s + alpha*ds;
    lambda = lambda + alpha*dlam;
    rp1 = A*x - b;
    rd1 = Q*x - A.'*lambda + c - s;
    s1 = norm(rp1-rp);
    s2 = norm(rd1-rd);
    rp = rp1;
    rd = rd1;
    obj = 0.5*x.'*Q*x + x.'*c;
    objlist(iter+1) = obj;
end

obj = objlist(iter);
x = x1;
iter = iter;
objlist = objlist(1:iter);
mesh=linspace(0,iter-1,iter);
% plot(mesh,objlist)
