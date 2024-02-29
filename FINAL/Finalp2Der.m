function [f,g,H] = Finalp2Der(x,mu,lambda,gamma,a,b,n)
% In the problem 2, use this function to find
% Lx, grad(Lx) and Hessian(Lx) at given pt xk


m = length(x); % size of x includes size of x and size of slack variables

% calculate value of f_k
f = 0;
% objective func part
for i=1:n-1
    f = f+(gamma+1)*x(i)^2-2*x(i)*x(i+1)+a(i)*x(i);
end
f = f + (gamma+1)*x(n)^2 + a(n)*x(n);
% equality constraint part
xsum = 0;
for i=1:n
    xsum=xsum+x(i);
end
f=f-lambda(n/2+1)*xsum+0.5*mu*(xsum)^2;
% inequality part
ineqt = zeros(n/2,1);
for i=1:n/2-1
    ineqt(i) = 0.5*x(2*i-1)^2+0.5*x(2*i+1)^2+x(2*i)^2-b(i)+x(n+i);
end
ineqt(n/2)=0.5*x(n-1)^2+x(n)^2-b(n/2)+x(m);
for j=1:n/2
    f=f+0.5*mu*ineqt(i)^2-lambda(i)*ineqt(i);
end


% calculate grad(f_k)
g = zeros(m,1);

g(1) = 2*(gamma+2)*x(1)-2*x(2)+a(1)-lambda(1)*x(1)+mu*ineqt(1)*x(1)- ...
lambda(n/2+1)+mu*xsum;

% for xi terms
for i=1:n/2-1
    % EVEN
    g(2*i) = 2*(gamma+2)*x(2*i)-2*x(2*i-1)-2*x(2*i+1)+a(2*i)- ...
    2*lambda(i)*x(2*i)+2*mu*ineqt(i)*x(2*i)-lambda(n/2+1)+mu*xsum;
    % ODD
    g(2*i+1) = 2*(gamma+2)*x(2*i+1)-2*x(2*i)-2*x(2*i+2)+a(2*i+1)- ...
    lambda(i)*x(2*i+1)-lambda(i+1)*x(2*i+1)+mu*ineqt(i)*x(2*i+1)+ ...
    +mu*ineqt(i+1)*x(2*i+1)-lambda(n/2+1)+mu*xsum;   
end
% for x_n
g(n) = 2*(gamma+2)*x(n)-2*x(n-1)+a(n)- ...
    2*lambda(n/2)*x(n)+2*mu*ineqt(n/2)*x(n)-lambda(n/2+1)+mu*xsum;
% for slack variables
for j=n+1:m
    g(j) = -lambda(j-n)+mu*ineqt(j-n);
end


% calculate Hessian matrix H_k
H = sparse(m,m);  % use sparse matrix instead of dense matrix for hessian
% i=1 & i=n
H(1,1)=2*(gamma+2)-lambda(1)+mu*ineqt(1);
H(1,2)=-2;
H(n,n)=2*(gamma+2)-2*lambda(n/2)+mu*ineqt(n/2);
H(n,n-1)=-2;

for i=1:n/2-1
    %even
    H(2*i,2*i)=2*(gamma+2)-2*lambda(i)+2*mu*ineqt(i);
    H(2*i,2*i-1)=-2;
    H(2*i,2*i+1)=-2;
    %odd
    H(2*i+1,2*i+1)=2*(gamma+2)-lambda(i)-lambda(i+1)+mu*ineqt(i)+mu*ineqt(i+1);
    H(2*i+1,2*i)=-2;
    H(2*i+1,2*i+2)=-2;
end
