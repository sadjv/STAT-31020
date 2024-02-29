m = [5,10,15];

% linear programming
avetime = [];
for k=1:3
    temp = [];
    n = m(k)*2;
    Q = zeros(n,n);
    for i=1:10
        A = rand(n,n);
        b = rand(n,1);
        c = rand(n,1);
        x0 = ones(n,1)*2;
        tic
        [x,obj,iter]=hwk7p1(Q,A,b,c,x0);
        toc
        temp(i) = toc;
    end
    avetime(k) = sum(temp)/10;
end
plot(m*2,avetime)
title('m vs time for linear programming')

% quadratic programming
avetime = [];
for k=1:3
    temp = [];
    n = m(k)*2;
    x = ones(n,1);
    for j=1:n
        x(j) = x(j) + (j-1)/n;
    end
    Q = diag(x);
    for i=1:10
        A = rand(n,n);
        b = rand(n,1);
        c = rand(n,1);
        x0 = ones(n,1)*2;
        tic
        [x,obj,iter]=hwk7p1(Q,A,b,c,x0);
        toc
        temp(i) = toc;
    end
    avetime(k) = sum(temp)/10;
end
plot(m*2,avetime)
title('m vs time for quadratic programming')