% Varying gap junctions, PLCgamma kinetics
% Gap Junctions: P_IP3 from 1e-3 to 1e

clear all
close all

runNumber = 1;
tic
for i = 10:-1:1
    a = 1;
    for v8 = [ 0.005:0.001:0.012 ]
        [ frequency(i,a), amplitude(i,a), width(i,a) ] = fHofer2002 (0.08, v8, 1);
        disp(['run number ' num2str(runNumber) ' of 600 done   seconds left: ' num2str((toc/runNumber)*(600-runNumber))])
        runNumber = runNumber + 1;
        a = a + 1;
    end
    %a = 1;
%     for P_IP3 = [0, 0.1, 0.25, 0.5, 1, 2, 3.5, 5, 7.5, 10]
%         [ frequency(i,2,a), amplitude(i,2,a), width(i,2,a) ] = fHofer2002 (0.08, 1e-2, P_IP3);
%         disp(['run number ' num2str(runNumber) ' of 600 done   seconds left: ' num2str((toc/runNumber)*(600-runNumber))])
%         runNumber = runNumber + 1;
%         a = a + 1;
%     end
%     a = 1;
%     for v8 = [0, 1e-4, 8e-4, 2e-3, 1e-2, 5e-2, 0.1, 0.5, 1, 3]
%         [ frequency(i,3,a), amplitude(i,3,a), width(i,3,a) ] = fHofer2002 (0.08, v8, 1);
%         disp(['run number ' num2str(runNumber) ' of 600 done   seconds left: ' num2str((toc/runNumber)*(600-runNumber))])
%         runNumber = runNumber + 1;
%         a = a + 1;
%     end
    
end