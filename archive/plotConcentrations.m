function plotConcentrations( xs, CaC, CaER, IP3, IP3R, domainSize, p )

        subplot(2,2,1)
        plot(xs,CaC);%,xs,P);
        axis([[0,domainSize] p.CaBound])
        title('Cyto Ca2+')
        subplot(2,2,2)
        plot(xs,CaER);%,xs,C);
        axis([[0,domainSize] p.CaERBound])
        title('ER Ca2+')
        subplot(2,2,3)
        plot(xs,IP3);%,xs,C);
        axis([[0,domainSize] p.IP3Bound])
        title('IP3')
        subplot(2,2,4)
        plot(xs,IP3R);%,xs,C);
        axis([[0,domainSize] p.IP3RBound])
        title('Active IP3R')


end

