function [x] = solverp2(gamma,mu,ub,lb,xk,lambda,nl,a,b)

% linear gradient projection method similar to the function
% in problem 1 (actually, the same as the func in problem
% except G and a).

% initialization
[L,dL,ddL] = Finalp2Der(xk,mu,lambda,gamma,a,b,nl);
n = length(xk);
G = ddL;
c = dL;
% lambda = ones(n,1);
x = xk;
g = G*x + c;
miter = 10;
count = 0;
% define A and b
bsize = n-nl;
A = zeros(2*nl+bsize,nl+bsize);  % KKT1
for i=1:nl
    A(i,i)=1;
    A(nl+i,i)=-1;
end
for i=1:bsize
    A(2*nl+i,nl+i)=1;
end
b = -ones(2*nl+bsize,1);
for i=1:bsize
    b(2*nl+i)=0;
end

for o =1:miter
    count = count + 1;
    % find cauchy point
    tlist = zeros(+1,1);
    for i=1:n   % find t^{bar} list
        ltemp = (x(i)-lb(i))/g(i);
        utemp = (x(i)-ub(i))/g(i);
        if g(i) < 0
            tlist(i+1) = utemp;
        elseif g(i) > 0
            tlist(i+1) = ltemp;
        else
            tlist(i+1) = Inf;
        end
    end
    tlistS = sort(unique(tlist, 'first'));
    tlist(1) = [];


    for i=1:length(tlistS)-1
        tj = tlistS(i+1);
        tjm1 = tlistS(i);
        xt = zeros(n,1);
        for j=1:n  % find x(t_{j-1})
            if tjm1 <= tlist(j)
                xt(j) = x(j) - tjm1*g(j);
            else
                xt(j) = x(j) - tlist(j)*g(j);
            end
        end
    
        p = -g; % find p_{i}^{j-1}
        for j=1:n
            if tjm1 >= tlist(j)
                p(j) = 0;
            end
        end
    
        % calculate fj-1, fj-1', fj-1''
        fjg = c.'*p + xt.'*G*p;
        fjh = p.'*G*p;
        deltat = -fjg/fjh;
    
        % find Cauchy point
        if fjg > 0
            xc = xt;
            break
        elseif deltat >= 0 && deltat < tj-tjm1
            tjm1 = tjm1 + deltat;
            xc = zeros(n,1);
            for j=1:n  % find x(t_{j-1})
                if tjm1 <= tlist(j)
                    xc(j) = x(j) - tjm1*g(j);
                else
                    xc(j) = x(j) - tlist(j)*g(j);
                end
            end
            break
        else
            if gamma <= 0 && o==1
                xc = -0.5*ones(n,1);
                % xc = x;
            else
                xc = x;
            end     
        end
    end

    % solve the subproblem
    % find Z and Y
    xtemp = linspace(1,n,n);
    ytemp = linspace(1,n,n);
    for i=1:n
        if xc(i)==lb(i) || xc(i)==ub(i)
            ytemp(i) = 0;
        else
            xtemp(i) = 0;
        end
    end
    ytemp = unique(ytemp,'first');
    if ytemp(1) == 0
        ytemp(1) = [];
    end
    nMm = length(ytemp);
    nMmtemp = linspace(1,nMm,nMm);
    v = ones(nMm,1);
    Z = sparse(ytemp, nMmtemp, v, n, nMm);
%     if gamma==-2
%         H = eye(n,n);
%     else
%         H = diag(abs(diag(G)));
%     end
    H = eye(n,n);
    P = Z*inv(Z.'*H*Z)*Z.';

    % projected CG method
    xp = xc;
    for i=1:n
        if xp(i)>1
            xp(i)=1;
        elseif xp(i)<-1
            xp(i)=-1;
        end
    end
    r = G*xp + c;
    gp = P*r;
    d = -gp;
    k = 0;


    while  0.5*xc.'*G*xc+xc.'*c <= 0.5*xp.'*G*xp+xp.'*c
        k = k+1;
        if k >= miter
            xp = xc;
            break
        end
        if d.'*G*d==0 || norm(gp)==0
            if k==1
                gp = r;
                d = -gp;
            else
                fprintf('zero break')
                break
            end 
        end
        alpha = r.'*gp/(d.'*G*d);
        % test = 0;
%         for i=1:nMm
%             if xp(ytemp(i))<=-1
%                 test = 1;
%             elseif  xp(ytemp(i))>=1
%                 test = 1;
%             end
%         end
%         if test==1
%             fprintf('done break')
%             break
%         end
        xp = xp + alpha*d;
        r1 = r + alpha*G*d;
        gp1 = P*r1;
        if norm(gp1)==0
            gp1=r1;
        end
        beta = r1.'*gp1/(r.'*gp);
        d = -gp1 + beta*d;
        gp = gp1;
        r = r1;
    end

    % correction of xp according to l<xi<u
    for i=1:n
        if xp(i)<lb(i)
            xp(i)=lb(i);
        elseif xp(i)>ub(i)
            xp(i)=ub(1);       
        end
    end
    
    qxp = 0.5*xp.'*G*xp+xp.'*c;
    qx = 0.5*x.'*G*x+x.'*c;
    xtemp = x;
    x = xp;
    g = G*x+c;
    f=A*x-b;
    
%     use KKT to break
%     KKT1 = Gx-A.'*lambda+c=0
%     KKT2 = Ax-b >=0
%     KKT3 = (Ax-b)_i*lambda_i=0, i=1...m
%     KKT4 = lambda_i>=0, i=1...m
    if norm(xtemp-x) < eps
        testK = 1;
        for i=1:length(b)
            if f(i) < 0 % KKT2
                testK=0;
            end
        end
        if testK==1
            % use KKT1, KKT3, KKT4 to find lagragian
            break   
        end
    elseif qx < qxp && count~=1
        qxp = qx;
        x = xtemp;
        g = G*x+c;
        f=A*x-b;
        testK = 1;
        for i=1:length(b)
            if f(i) < 0 % KKT2
                testK=0;
            end
        end
        if testK==1
            % use KKT1, KKT3, KKT4 to find lagragian
            break   
        end
    end
end
