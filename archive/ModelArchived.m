% This function contains default parameter values for the calcium wave
% simulation code. It should be called at the beginning of the code to make
% the code easier to read.
% ALL PARAMETERS SHOULD BE IN PRODUCTS OF THE FOLLOWING UNITS:
% time:             [s]
% length:           [um]
% concentration:    [uM]
% particle number:  [1e-18 moles] = [602214 particles]

function model = ModelArchived(cellSize)
%% Constants
p.MM_Ca       = 40.078;       % [g/mol]   Molar mass of Calcium atom (± 0.004 u)
p.MM_IP3      = 420.096;   	% [g/mol]   Molar mass of IP3 molecule

%% Kinetic parameters
p.v_PLC     = 1e-2;   	% [uM/s]    Rate of PLC-gamma
p.K_Ca  	= 0.2;    	% [uM]      Half-saturation constant for calcium activation of IP3R
p.K_IP3    	= 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
p.K_i    	= 0.09;    	% [uM]      Half-saturation constant for calcium inhibition of IP3R
p.k_1      	= 4e-4;     % [1/s]     Rate constant of calcium leak from ER
p.k_2     	= 0.08;     % [1/s]     Rate constant of calcium release through IP3
p.k_SERCA 	= 0.02;     % [uM]      Half-saturation constant for SERCA pump
p.v_SERCA 	= 0.08;     % [uM/s]    Maximal rate for SERCA pump
p.k_i     	= 13.3;   	% [1/s]     Rate constant of IP3R inactivation
p.k_deg  	= 0.08;     % [1/s]     Rate constant of IP3 degradation     
p.k_EC    	= 0.035;    % [uM]      Saturation constant for calcium transfer with media
p.v_EC     	= 0.51;     % [1/s]     Rate constant for calcium transfer with media
p.v_gen     = 0.01;

%% Physical parameters
p.D_IP3 	= 280;	% [uM/s^2]  Diffusion coefficient of IP3	
p.D_Ca    	= 20; 	% [uM/s^2]	Effective diffusion coefficient of calcium in cytosol
p.P_IP3   	= 1;   	% [uM/s]    Gap-junctional permeability of IP3
p.P_Ca   	= 0.01;	% [uM/s]    Effective gap-junctional permeability of calcium
p.beta    	= 20; 	%           Ratio of the effective volumes for calcium of cytoplasm and ER

%% Initialization
% Homogenize diffusivity
p.Deff_Ca = 1./(1./p.D_Ca + 1./(p.P_Ca .* cellSize));
p.Deff_IP3 = 1./(1./p.D_IP3 + 1./(p.P_IP3 .* cellSize));


%% Describe system of PDE's in m, d, c, a, f form
% No second derivatives with respect to t
m = 0;

% Describes each equation as a generation rate (du/dt = ...)
d = [1, 1, 1, 1]';
 
% Adds diffusion terms (del2(u))
c = [p.Deff_Ca ,p.Deff_IP3, 0, 0]';
 
% Adds decay terms (-ku)
a = [p.k_EC, p.k_deg, 0, p.k_i]';

% Adds flux terms
f = @(region, state)pdeFluxEquations(region, state, p);

% Create empty pde model
model = createpde(4);
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);
end

%% Flux equations
function f = pdeFluxEquations(region, state, p)
% Calculate squares of C and I to reduce computation
C2 = state.u(1,:).^2;
I2 = state.u(2,:).^2;

% Calculate flux terms to reduce compuatation
J_flux = (p.k_1+p.k_2.*state.u(4,:).*C2.*I2) ...
    .*(state.u(3,:)-state.u(1,:)) ...
    ./(p.K_Ca^2+C2) ...
    ./(p.K_IP3^2+I2);
J_SERCA = p.v_SERCA .* C2 ./ (p.k_SERCA^2 + C2);

% Solve for remaining terms
f(1,:) = J_flux - J_SERCA + p.v_EC;
f(2,:) = p.v_gen;  %gamrnd(p.gammaShape, p.theta);
f(3,:) = -p.beta*(J_flux - J_SERCA);
f(4,:) = p.k_i*(p.K_i^2./(p.K_i^2+C2));
end

