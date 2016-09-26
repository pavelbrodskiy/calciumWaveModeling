function [ outputFlag ] = evaluateFlags( CaC, CaER, IP3, IP3R, maxIP3 )

outputFlag = 0;
if min([CaC CaER IP3 IP3R]) < 0
    outputFlag = 1;
    disp('Negative Concentration');
elseif sum(IP3>maxIP3)
    outputFlag = 2;
    disp('IP3 is exploding');
elseif sum(isnan([CaC CaER IP3 IP3R]))
    outputFlag = 3;
    disp('NaN values in simulation');
end

end

