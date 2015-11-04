function plotConcentrations( xs, CaC, CaER, IP3, IP3R, domainSize )

        subplot(1,4,1)
        plot(xs,CaC);%,xs,P);
        axis([0,domainSize,0,1])
        title('Cyto Ca2+')
        subplot(1,4,2)
        plot(xs,CaER);%,xs,C);
        axis([0,domainSize,0,60])
        title('ER Ca2+')
        subplot(1,4,3)
        plot(xs,IP3);%,xs,C);
        axis([0,domainSize,0,0.05])
        title('IP3')
        subplot(1,4,4)
        plot(xs,IP3R);%,xs,C);
        axis([0,domainSize,0,1])
        title('Active IP3R')


end

