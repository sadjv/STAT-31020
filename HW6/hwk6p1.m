function x=hwk6p1(x_in,tol,iter)


[n,m] = size(x_in);
[f,g,H] = rosenbrocknfgHS(x_in);
x = x_in;
i = 1;
k = 1;
error = [];
gamma = 0.15;


while norm(g) > tol || k<=iter
    deltahat = 2.5;
    deltak = 1.5;
    k = k+1;
    ek = min(0.5, sqrt(norm(g)));
    z1 = zeros(n,1);
    r1 = g;
    d1 = -r1;
    if d1.'*H*d1 <=0
        pk = -g*(deltak/norm(g));
        break
    else
        % first stopping criteria
        while d1.'*H*d1 > 0
            alpha = r1.'*r1/(d1.'*H*d1);
            z2 = z1 + alpha*d1;
            r2 = r1 + alpha*H*d1;
            % second stopping criteria
            if norm(z2) >= deltak
                % solve for norm(pk) = deltak
                r = roots([norm(d1)^2,2*z1.'*d1,deltak^2-norm(z1)^2]);
                r = sqrt(real(r).^2+imag(r).^2);
                h = r(1);
                pk = z1 + h*d1;
                [f2,g2,H2] = rosenbrocknfgHS(pk);
                i = i+1;
                mk = f2 + g2.'*pk + 0.5*pk.'*H2*pk;
                for k=2:length(r)
                    ptemp = z1 + r(k)*d1;
                    [f2,g2,H2] = rosenbrocknfgHS(ptemp);
                    i = i+1;
                    mtemp = f2 + g2.'*pk + 0.5*pk.'*H2*pk;
                    if mtemp > mk && r(k)>=0
                        pk = ptemp;
                        h = r(k);
                        mk = mtemp;
                    end
                end
                break
            end
            if norm(r2) < ek
                pk = z2;
                break
            end
            beta = r2.'*r2/(r1.'*r1);
            d1 = -r2 + beta*d1;
            r1 = r2;
            z1 = z2;
        end
        % update pk for one more time
        r = roots([norm(d1)^2,2*z1.'*d1,deltak^2-norm(z1)^2]);
        r = sqrt(real(r).^2 + imag(r).^2);
        h = r(1);
        pk = z1 + h*d1;
        [f2,g2,H2] = rosenbrocknfgHS(pk);
        i = i+1;
        mk = f2 + g2.'*pk + 0.5*pk.'*H2*pk;
        for k=2:length(r)
            ptemp = z1 + r(k)*d1;
            [f2,g2,H2] = rosenbrocknfgHS(ptemp);
            i = i+1;
            mtemp = f2 + g2.'*pk + 0.5*pk.'*H2*pk;
            if mtemp > mk && r(k)>=0
                pk = ptemp;
                mk = mtemp;
            end
        end
    end
    
    x = x + pk;
    [f1,g1,H1] = rosenbrocknfgHS(x);
    i = i+1;
    
    % update deltak
    rho = (f - f1)/(f-(f1+g1.'*pk+0.5*pk.'*H1*pk));
    if rho<0.25
        deltak = 0.25*deltak;
    else
        if rho > 0.75 && norm(pk)==deltak
            deltak = min(2*deltak,deltahat);
        else
            deltak=deltak;
                
        end
        if rho>gamma
            x = x+pk;
        end
    end
    
    
    [f1,g1,H1] = rosenbrocknfgHS(x);
    i = i+1;
    error(k) = norm(g1)/norm(g);
    f = f1;
    g = g1;
    H = H1;
    i = i+1;
end
