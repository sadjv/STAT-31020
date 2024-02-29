% implementation of the function hwk5p2

mesh = linspace(30,210,10);
numfunc = zeros(18,1);
tol = 10^-4;
iter = 1000;

for i=1:10
    x_in = ones(mesh(i),1)*3;
    disp('for n=')
    mesh(i)
    disp('the number of functions evaluated is i')
    x=hwk5p2(x_in,tol,iter);
end

% for all n, the outputing x is indeed [1,1,1,1,1,...,1]

% to check superlinearity
% if norm(grad(f_{k+1}))/norm(grad(f_{k})) approaches 0 for large n
% then the algorithm converges superlinearly.
x=hwk5p2(x_in,10^-3,iter);

% error =
% 
%   Columns 1 through 8
% 
%          0    0.2923    0.4752    0.1686    1.8983    0.0537    1.8329    0.0492
% 
%   Columns 9 through 16
% 
%     4.9989    0.2200    2.4035    0.6118    1.5229    0.2823    1.9897    0.5303
% 
%   Columns 17 through 24
% 
%     1.6155    0.2196    2.6033    0.2217    1.8032    0.3428    0.9474    0.1225
% 
%   Columns 25 through 27
% 
%     0.2690    0.0046    0.0004

% it does converges superlinearly