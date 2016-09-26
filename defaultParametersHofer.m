function [ p ] = defaultParameters()
% This function contains default parameter values for the calcium wave
% simulation code. It should be called at the beginning of the code to make
% the code easier to read.
% ALL PARAMETERS SHOULD BE IN PRODUCTS OF THE FOLLOWING UNITS:
% time:             [s]
% length:           [um]
% concentration:    [uM]
% particle number:  [1e-18 moles] = [602214 particles]

%% Constants
MM_Ca       = 40.078;       % [g/mol]   Molar mass of Calcium atom (± 0.004 u)
MM_IP3      = 420.096;   	% [g/mol]   Molar mass of IP3 molecule

%% Sources for parameter values
% Hofer model paper:




%% Wave simulation parameters
% Kinetic parameters
p.k_1       = 0.0004;       % [1/s]     Rate constant of calcium leak from ER
p.k_2       = 0.08;         % [1/s]     Rate constant of calcium release through IP3
p.k_3       = 0.8;        	% [1/s]     Rate constant of SERCA pump
p.k_5       = 0.5;          % [1/s]     Rate constant of calcium extrusion
p.k_6       = 4;            % [1/s]     Rate constant of IP3R inactivation
p.k_9       = 0.08;       	% [1/s]     Rate constant of IP3 degradation

p.v_40      = 0.025;        % [uM/s]    Rate of calcium leak across the plasma membrane
p.v_41      = 0.2;          % [uM/s]    Maximal rate of activation-dependent calcium influx

p.K_IP3     = 0.3;          % [uM]      Half-saturation constant for IP3 activation of IP3R
p.K_a       = 0.2;          % [uM]      Half-saturation constant for calcium activation of IP3R
p.K_i       = 0.2;          % [uM]      Half-saturation constant for calcium inhibition of IP3R
p.K_r       = 1;            % [uM]      Half-saturation constant for agonist-dependent calcium entry

p.D_IP3     = 280;          % [uM/s^2]  Diffusion coefficient of IP3
p.D_Ca      = 20;           % [uM/s^2]	Diffusion coefficient of calcium

p.P_IP3     = 1;            % [uM/s]    Gap-junctional permeability of IP3
p.P_Ca      = 0.01*p.P_IP3;	% [uM/s]    Effective gap-junctional permeability of calcium

p.beta      = 20;           % [no unit]	Ratio of the effective volumes for calcium of cytoplasm and ER

p.V_CaInd   = 1e-2;
p.V_Ca      = 8e-2;
p.K_Ca      = 3e-1;

% Postulated PLC kinetics
p.PLC = @(Ca2) p.V_CaInd + p.V_Ca.*Ca2./(p.K_Ca.^2+Ca2);

% Stochastic pulse parameters
p.pulseTimeConstant = 30;    % [s]       Average time between pulses
p.IP3Pulse          = 5;    % [uM]      Concentration of IP3 in affected coordinates
p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse

% Initial conditions
p.Ca_cyto_0         = 0;	% [uM]      Initial concentration of Ca2+ in cytoplasm
p.Ca_ER_0           = 80;	% [uM]      Initial concentration of Ca2+ in ER
p.IP3_0             = 0;	% [uM]      Initial concentration of IP3
p.IP3R_0            = 1;	% [uM]      Initial active concentration of IP3R

end

