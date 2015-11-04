% Implementation of De Young and Keizer 1992 model
% Pavel Brodskiy, Zartman Lab, 10.21.15
% http://www.pnas.org/content/89/20/9895.full.pdf

%% PARAMETER DECLARATION
c0  = 2.0;      % [uM]      Total [Ca2+] in terms of cytosolic vol

c1  = 0.185;    %           (ER vol)/(cytosolic vol)
v1  = 6;        % [1/s]     Max Ca2+ channel flux
v2  = 0.11;     % [1/s]     Ca2+ leak flux constant
v3  = 0.9;      % [1/uM*s]  Max Ca2+ uptake
k3  = 0.1;      % [uM]      Activation constant for ATP-Ca2+ pump

% Receptor binding constants
a1  = 400;      % [1/uM*s]  IP3
a2  = 0.2;      % [1/uM*s]  Ca2+ (inhibition)
a3  = 400;      % [1/uM*s]  IP3
a4  = 0.2;      % [1/uM*s]  Ca2+ (inhibition)
a5  = 20;       % [1/uM*s]  Ca2+ (activation)

% Receptor dissociation constants (d = b / a)
d1  = 0.13;     % [uM]      IP3
d2  = 1.049;    % [uM]      Ca2+ (inhibition)
d3  = 943.4;    % [nM]      IP3
d4  = 144.5;    % [nM]      Ca2+ (inhibition)
d5  = 82.34;    % [nM]      Ca2+ (activation)

%% MODEL PARAMETER DECLARATION
totalTime 	= 100;          % [s]   Total simulation time
domainSize  = 100;          % [um]  Total simulation domain size
dt          = 0.1;          % [s]   Simulation timestep
dx          = 0.5;          % [um]  Simulation spatial stepsize

IP30        = 0.5;          % [uM]  Initial concentration of IP3
Cai0        = 0.5;          % [uM]  Initial concentration of Ca2+ in cytoplasm
CaER0       = 0.5;          % [uM]  Initial concentration of Ca2+ in ER

%% INITIALIZATION

% Convert to modeling units (uM, um, s)
d3      = d3 / 1e3;     % nM to uM
d4      = d4 / 1e3;     % nM to uM
d5      = d5 / 1e3;     % nM to uM

% Initialize variables
xSize   = domainSize / dx;
ySize   = domainSize / dx;
IP3     = ones(1, xSize) * IP30;
Cai     = ones(1, xSize) * Cai0;
CaER    = ones(1, xSize) * CaER0;





%% SIMULATION LOOP

for t = 0:dt:totalTime
    
    x110	= Ca.*IP3.*d2/((Ca.*IP3+IP3*d2+d1*d2+Ca*d3).*(Ca+d5)); % IP3 probability
    
    % IP3R Rates
    V1      = a4*(Cai*x000 - d4*x001);
    V2      = a4*(Cai*x010 - d4*x011);
    V3      = a5*(Cai*x000 - d5*x010);
    V4      = a5*(Cai*x001 - d5*x011);
    
    % Fluxes
    J1      = c1*(v1*x110.^3+v2).*(CaER-Cai); % Outward flux (receptor + leak)
    J2      = (v3*Cai.^2)./(Cai.^2+k3^2); % Inward flux(ATP Ca pump
    
    % Forward Euler calculations
    dCaidt  = J1 - J2;
    dCaERdt = ;
    dIP3dt  = v4 * ((Cai+(1-alpha)*k4)./(Cai+k4)) - Ir * IP3;
    dx000dt = -V1 - V3;
    dx001dt = V1 - V4;
    dx010dt = V3 - V2;
    
    
    
    
    
    
    
    
    
    
    
end







