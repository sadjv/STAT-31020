x1=hwk3p3([5,5,5,5,5],1000,10^-3);

x=hwk3p3([10,10,10,10,10,10,10],1000,10^-3);



mesh = linspace(5,50,10);
timel = zeros(10,1);
for i = 1:5
    x = zeros(i*5,1) + 10;
    x = x.'
    tic
    hwk3p3(x,1000,10^-3);
    toc
    timel(i) = toc;
end

plot(mesh,timel)
title('running time vs n')