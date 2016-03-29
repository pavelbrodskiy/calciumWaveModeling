function [ p ] = defaultParameters(varargin)
% This function contains default parameter values for the calcium wave
% simulation code. It should be called at the beginning of the code to make
% the code easier to read.
% ALL PARAMETERS SHOULD BE IN PRODUCTS OF THE FOLLOWING UNITS:
% time:             [s]
% length:           [um]
% concentration:    [uM]
% particle number:  [1e-18 moles] = [602214 particles]

%% FINAL AICHE MODEL
% Kinetic parameters
p.v_PLC         = 1e-2;   	% [uM/s]    Rate of PLC-gamma
p.gammaShape    = 1;        % [  ]      Shape constant for noise term
p.K_Ca       	= 0.2;    	% [uM]      Half-saturation constant for calcium activation of IP3R
p.K_IP3         = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
p.K_i           = 0.09;    	% [uM]      Half-saturation constant for calcium inhibition of IP3R
p.k_1           = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
p.k_2           = 0.08;     % [1/s]     Rate constant of calcium release through IP3
p.k_gam         = 0.02;     % [uM]      Half-saturation constant for SERCA pump
p.gam           = 0.08;     % [uM/s]    Maximal rate for SERCA pump
p.k_6           = 13.3;   	% [1/s]     Rate constant of IP3R inactivation
p.k_9           = 0.08;     % [1/s]     Rate constant of IP3 degradation     
p.Ca_media      = 0.035;    % [uM]      Saturation constant for calcium transfer with media
p.P_Ca_media    = 0.51;     % [1/s]     Rate constant for calcium transfer with media

% Physical parameters
p.D_IP3         = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
p.D_Ca          = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
p.P_IP3         = 1;     	% [uM/s]    Gap-junctional permeability of IP3
p.P_Ca          = 0.01;     % [uM/s]    Effective gap-junctional permeability of calcium
p.beta          = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Intial conditions
p.Ca_cyto_0   	= 0.05;    	% [uM]      Initial concentration of Ca2+ in cytoplasm
p.Ca_ER_0    	= 80;     	% [uM]      Initial concentration of Ca2+ in ER
p.IP3_0      	= 0.15;    	% [uM]      Initial concentration of IP3
p.IP3R_0      	= 1;      	% [uM]      Initial active concentration of IP3R

% Random pulse parameters
p.pulseTimeCon	= @(t) inf; % [s]       Average time between pulses as a function of time
p.IP3Pulse    	= 5;        % [uM]      Concentration of IP3 in affected coordinates
p.IP3Extent   	= 1;        % [um]      Extent of IP3 pulse

if nargin == 0
    calciumWaveSimulation();
end

end

