function [Pxk]=Projp2(xk,pk,lb,ub)
m = length(xk);

Pxk = xk+pk;
for i=1:m
    if Pxk(i)<lb(i)
        Pxk(i)=lb(i);
    elseif Pxk(i)>ub(i)
        Pxk(i)=ub(i);
    end
end