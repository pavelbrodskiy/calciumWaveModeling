close all

mean = 1e-2;
density = 1e4;
maxX = 5e-2;

x=0:maxX/density:maxX;

k = 1;
theta = mean ./ k;
gamma1=pdf('Gamma',x,k,theta);

k = 0.01;
theta = mean ./ k;
gamma2=pdf('Gamma',x,k,theta);

plot(x,gamma1,x,gamma2);
legend('k = 1.00','k = 0.01');

plot(x,gamma2);