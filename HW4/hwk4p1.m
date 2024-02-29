function [x]=hwk4p1(x0,niter,eps)

deltahat = 5;
delta = 4;
goma = 1/5;
tol = eps;
[f,g,H] = rosenbrocknfgHS(x0);
x = x0;
i = 1;

while norm(g)<=tol | i<=niter
    i = i+1;
    % Cauchy point
    [f,g,H] = rosenbrocknfgHS(x);
    if g*H*g.' <= 0
        tao = 1;
    else
        tao = min(norm(g)^3/(delta*g*H*g.'),1);
    end
    pc = -tao*delta/(norm(g))*g; % 1*n
    
    % dogleg
    pb = -H/g;
    
   if norm(pc) == delta
        pk = pc;
        break
    elseif norm(pc) < delta
    
        tspan = linspace(1,2,1000);
        rlist = [];
        for i=1:1000
            k = norm(pc+(tspan(i)-1)*(pb.'-pc));
            if k<=delta
                rlist(i)=tspan(i);  
            end
        end
        mpt = @(x) f+g*(pc+(x-1)*(pb.'-pc)).'+1/2*(pc+(x-1)*(pb.'-pc))*H*(pc+(x-1)*(pb.'-pc)).';
        mpteva = [];
        [m,~] = size(rlist);
        for i=1:m
            mpteva(i)=mpt(rlist(i));
        end
        [~,t] = min(mpteva);
        ta = rlist(t);
        b=mpt(ta);
        a=f+g*tao*pc.'+1/2*tao*pc*H*pc.';
        if a>b
            pk = pc;
        else
            pk = pc+(t-1)*(pb.'-pc);
        end
    
    try chol(H)  
    catch ME
       pk = pc;
    end
    
    if cond(H)>10^3
       pk = pc;
    end

    
    [f1,g1,H1] = rosenbrocknfgHS(x+pk);
    rho = (f-f1)/(-g*pc.'-0.5*pc*H*pc.');
        if rho<0.25
            delta = 0.25*delta;
        else
            if rho > 0.75 && norm(pk)==delta
                delta = min(2*delta,deltahat);
            else
                delta=delta;
                
            end
            if rho>goma
                x = x+pk;
            else
                x = x;
            end
        end
         
    end
    
    
    
    
    
end