opts = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','objective',...
 @(x) OptimizeSummaryStatisticsFromSignal(x),'x0',[0.3,25],'lb',[-0.5, 1],'ub',[1.5, 1000],'options',opts);
ms = MultiStart;
ms.UseParallel = true;
[x,f] = run(ms,problem,100);