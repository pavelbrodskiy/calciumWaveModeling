function [ xPulseCoords, yPulseCoords ] = pulseCoordinates( xSize, ySize, IP3PulseCoords, p )
xPulse = round(rand * (xSize - 1)) + 1;
yPulse = round(rand * (ySize - 1)) + 1;
xPulseCoords = round(xPulse-IP3PulseCoords/2):round(xPulse+IP3PulseCoords/2);
yPulseCoords = round(yPulse-IP3PulseCoords/2):round(yPulse+IP3PulseCoords/2);
switch p.boundCondition
    case 'per'
        xBig = sum(xPulseCoords < 1);
        yBig = sum(yPulseCoords < 1);
        xSmall = sum(xPulseCoords > xSize);
        ySmall = sum(yPulseCoords > ySize);
        
        if xBig; xPulseCoords = [xPulseCoords (xSize-xBig):(xSize)]; end
        if yBig; yPulseCoords = [yPulseCoords (ySize-yBig):(ySize)]; end
        if xSmall; xPulseCoords = [xPulseCoords 1:xSmall]; end
        if ySmall; yPulseCoords = [yPulseCoords 1:ySmall]; end
    case 'noflux'
        
    otherwise
        error('Boundary conditions unspecified');
end

xPulseCoords(xPulseCoords < 1) = [];
yPulseCoords(yPulseCoords < 1) = [];
xPulseCoords(xPulseCoords > xSize) = [];
yPulseCoords(yPulseCoords > ySize) = [];

end

