% implementation of the function final1

eps = 10^-4;
nlist = [10,20,30];

tmean = zeros(3,1);
for i=1:3
    n = nlist(i);
    ttemp = zeros(100,1);
    for j=1:20
        gamma = rand(1,1)*4-2;
        a = rand(n,1)*2-1;
        tic
        [x,lambda]=finalp1(n,gamma,a,eps);
        toc
        ttemp(j)=toc;
    end
    tmean(i)=sum(ttemp)/20;
end

plot(nlist,tmean)
title('mean running time vs n for final1')
xlabel('n')
ylabel('running time in second')


