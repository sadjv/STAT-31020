ttoc = 5;
i = 1;
n=200;
x = zeros(n,1)*2;
b = ones(n,1);
[~,~,HD] = rosenbrocknfgH(x);
[~,~,HS] = rosenbrocknfgHS(x);

tic
HD\b;
toc
t1 = toc;

tic
HS\b;
toc
t2 = toc;

tlistd = [t1];
tlists = [t2];
nlist = [n];
while max(t1,t2)<=ttoc
    i=i+1
    n = n+1000;
    x = zeros(n,1)*2;
    b = ones(n,1);
    [~,~,HD] = rosenbrocknfgH(x);
    [~,~,HS] = rosenbrocknfgHS(x);

    tic
    h=HD\b;
    toc
    t1 = toc;

    tic
    h=HS\b;
    toc
    t2 = toc;
    tlistd(i) = t1;
    tlists(i) = t2;
    nlist(i) = n;
end

hold on
plot(nlist,tlistd)
plot(nlist,tlists)
legend('a(n)','b(n)')
xlabel('n')
ylabel('time')
title('question 1 part b: time for dense vs sparse')
hold off





