function plotOneConcentration( xs, CaC, domainSize, p )

plot(xs,CaC);
axis([[0,domainSize], p.CaBound])
title('Cyto Ca2+')

end

