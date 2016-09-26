function [ p ] = defaultParameters()
% This function contains default parameter values for the calcium wave
% simulation code. It should be called at the beginning of the code to make
% the code easier to read.
% ALL PARAMETERS SHOULD BE IN PRODUCTS OF THE FOLLOWING UNITS:
% time:             [s]
% length:           [um]
% concentration:    [uM]
% particle number:  [1e-18 moles] = [602214 particles]

% % %% Constants
% % MM_Ca       = 40.078;       % [g/mol]   Molar mass of Calcium atom (± 0.004 u)
% % MM_IP3      = 420.096;   	% [g/mol]   Molar mass of IP3 molecule
% % 
% % %% Cody's paper parameters
% % % p.n         = 1;            % [   ]
% % % p.k_flux    = .5;           % [uM]      3 uM/sec affects intensity of Ca2+ flash %TUNED
% % % p.k_mu      = 0.01;         % [uM]
% % % p.b         = 0.11;         % [   ]
% % % p.k_1       = 0.70;         % [uM]
% % % p.gamma     = 1.10;         % [uM/s]
% % % p.k_gamma 	= 0.27;         % [uM]
% % % p.V_P       = 0.08;         % [1/s]     %TUNED
% % % p.k_P       = 0.5*1.00;     % [uM]
% % 
% % %% Assumed parameters
% % % p.v_PLC     = 8e-5;
% % 
% % % Emperical IP3R Fit
% % % p.v1        = 1.02;
% % % p.k1        = 0.01;
% % 
% % % Postulated PLC kinetics
% % % p.PLC = @(Ca2) p.V_CaInd + p.V_Ca.*Ca2./(p.K_Ca.^2+Ca2);
% % 
% % % Stochastic pulse parameters
% % % p.pulseTimeConstant = 10;	% [s]       Average time between pulses
% % % p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% % % p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% % 
% % %% Wave simulation parameters from Hofer paper
% % % Kinetic parameters
% % % p.k_1       = 0.0004;       % [1/s]     Rate constant of calcium leak from ER
% % % p.k_2       = 0.08;         % [1/s]     Rate constant of calcium release through IP3
% % % p.k_3       = 0.8;        	% [1/s]     Rate constant of SERCA pump
% % % p.k_5       = 0.5;          % [1/s]     Rate constant of calcium extrusion
% % % p.k_6       = 4;            % [1/s]     Rate constant of IP3R inactivation
% % % p.k_9       = 0.08;       	% [1/s]     Rate constant of IP3 degradation
% % % 
% % % p.v_40      = 0.025;        % [uM/s]    Rate of calcium leak across the plasma membrane
% % % p.v_41      = 0.2;          % [uM/s]    Maximal rate of activation-dependent calcium influx
% % % 
% % % p.K_IP3     = 0.3;          % [uM]      Half-saturation constant for IP3 activation of IP3R
% % % p.K_a       = 0.2;          % [uM]      Half-saturation constant for calcium activation of IP3R
% % % p.K_i       = 0.2;          % [uM]      Half-saturation constant for calcium inhibition of IP3R
% % % p.K_r       = 1;            % [uM]      Half-saturation constant for agonist-dependent calcium entry
% % % 
% % % 
% % % 
% % % p.beta      = 20;           % [   ]     Ratio of the effective volumes for calcium of cytoplasm and ER
% % % 
% % % p.V_CaInd   = 1e-2;
% % % p.V_Ca      = 8e-2;
% % % p.K_Ca      = 3e-1;
% % 
% % % Initial conditions
% % % p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% % % p.Ca_ER_0	= 60;           % [uM]      Initial concentration of Ca2+ in ER
% % % p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% % % p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% % 
% % %% REVISED PARAMETERS 11.6.15
% % % Stochastic pulse parameters
% % 
% % %p.k_1       = 0.0063;    	% [1/s]     Scaled, sneyd k_flux/(CaER-Ca)
% % %p.k_2       = 0.01;         % [uM]      k_mu from sneyd model
% % %p.k_3       = 0.70;         % [uM]      k_1 from sneyd
% % % p.k_4       = 1e-1;         % [uM/s]    gamma from sneyd
% % % p.k_5    	= 1e-2;         % [uM]      k_gamma from sneyd
% % % p.k_6       = 0.08;       	% [1/s]     Rate constant of IP3 degradation 
% % % p.k_7       = 4;            % [1/s]     Rate constant of IP3R inactivation in hofer (1/tau)
% % % p.k_8       = 0.2;          % [uM]      k_2 from sneyd is 0.7, K_i is 0.2 in hofer
% % %  
% % % p.beta      = 20;           % [  ]      Ratio of the effective volumes for calcium of cytoplasm and ER
% % % 
% % % p.v_PLC     = 1e-2;         % [uM/s]    Rate of IP3 regeneration
% % % %p.v_leak    = 0.0004;       % [1/s]     k_1 in Hofer
% % % 
% % % 
% % % p.b         = 0.11;          % [  ]      Fraction of open IP3R when no IP3 is present
% % % 
% % % p.D_IP3     = 280;          % [uM/s^2]  Diffusion coefficient of IP3
% % % p.D_Ca      = 20;           % [uM/s^2]	Diffusion coefficient of calcium
% % % p.P_IP3     = 1.001;            % [uM/s]    Gap-junctional permeability of IP3
% % % p.P_Ca      = 0.01*p.P_IP3;	% [uM/s]    Effective gap-junctional permeability of calcium
% % % 
% % % %p.n         = 1;
% % % 
% % %  p.k_1       = 0.0004;       % [1/s]     Rate constant of calcium leak from ER
% % %  p.k_2       = 0.08;         % [1/s]     Rate constant of calcium release through IP3
% % %  p.K_IP3     = 0.3;          % [uM]      Half-saturation constant for IP3 activation of IP3R
% % %  p.K_Ca      = 0.2;          % [uM]      Half-saturation constant for calcium activation of IP3R
% % 
% % %% HOFFER PARAMETERS
% % 
% % % Hofer Model Parameters
% % p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% % p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% % p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% % p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% % p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% % p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% % p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% % p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% % p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% % p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% % p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% % p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% % 
% % % Tunable Parameters
% % p.k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
% % p.v_7     = 0.08;         % [uM/s]    Maximal rate of PLCd
% % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % p.K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % 
% % 
% % p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% % p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% % p.alpha_0 = 1.0;          
% % p.Kappa_g = 1.0;
% % 
% % p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% % p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% % p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% % p.sigma     = 1;         % [   ]     Standard deviation of production of IP3
% % p.sigma2        =0; % DO NOT USE OR GHOSTS
% % 
% % p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% % p.Ca_ER_0	= 60;           % [uM]      Initial concentration of Ca2+ in ER
% % p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% % p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% % 
% % %% HOFFER PARAMETERS TO PLAY WITH
% % 
% % % Hofer Model Parameters
% % p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% % p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% % p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% % p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% % p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% % p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% % p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% % p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% % p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% % p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% % p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% % p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% % 
% % % Tunable Parameters
% % % p.k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
% % % p.v_7     = 1;         % [uM/s]    Maximal rate of PLCd
% % % p.v_PLC   = 0;         % [uM/s]    Rate of PLCb
% % % p.K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % p.k_3     = 2;          % [1/s]     Rate constant of SERCA pump
% % % % p.v_7     = 0.08;         % [uM/s]    Maximal rate of PLCd
% % % % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % % % p.K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % % p.k_gam   =     1e-2;    %uM
% % % % p.gam     =     p.k_gam*3;    %uM/sec
% % % % % p.v_7     = 1e0;         % [uM/s]    Maximal rate of PLCd
% % % % % p.v_PLC   = 0;         % [uM/s]    Rate of PLCb
% % % % %% p.K_Ca    = p.v_7 .* 0.1;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % % % p.k_gam   =   p.v_7 .*  .5e-2;    %uM
% % % % % p.gam     =     p.k_gam*3;    %uM/sec
% % p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% % % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % %p.K_Ca    = p.v_7 .* 0.1;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % p.k_gam   =     3e-2;    %uM
% % p.gam     =     p.k_gam*2;    %uM/sec
% % %p.k_gam   =     .27;    %uM
% % %p.gam     =     1;    %uM/sec
% % 
% % p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% % p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% % p.alpha_0 = 1.0;          
% % p.Kappa_g = 1.0;
% % 
% % p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% % %p.pulseTimeConstant = @(t) (t>0).*1e20;	% [s]       Average time between pulses
% % p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% % p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% % p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% % p.sigma2        =0; % DO NOT USE OR GHOSTS
% % 
% % p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% % p.Ca_ER_0	= 80;           % [uM]      Initial concentration of Ca2+ in ER
% % p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% % p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% % 
% % p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% % 
% % %% HOFFER PARAMETERS TO MURDER v_40, v_41, and k_5
% % 
% % % Hofer Model Parameters
% % p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% % p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% % p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% % p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% % p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% % p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% % p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% % p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% % p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% % p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% % p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% % p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% % 
% % % Tunable Parameters
% % % p.k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
% % % p.v_7     = 1;         % [uM/s]    Maximal rate of PLCd
% % % p.v_PLC   = 0;         % [uM/s]    Rate of PLCb
% % % p.K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % p.k_3     = 2;          % [1/s]     Rate constant of SERCA pump
% % % % p.v_7     = 0.08;         % [uM/s]    Maximal rate of PLCd
% % % % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % % % p.K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % % p.k_gam   =     1e-2;    %uM
% % % % p.gam     =     p.k_gam*3;    %uM/sec
% % % % % p.v_7     = 1e0;         % [uM/s]    Maximal rate of PLCd
% % % % % p.v_PLC   = 0;         % [uM/s]    Rate of PLCb
% % % % %% p.K_Ca    = p.v_7 .* 0.1;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % % % % p.k_gam   =   p.v_7 .*  .5e-2;    %uM
% % % % % p.gam     =     p.k_gam*3;    %uM/sec
% % p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% % % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % %p.K_Ca    = p.v_7 .* 0.1;      % [uM]      Half-saturation constant for calcium activation of PLCd
% % p.k_gam   =     0.27;    %uM
% % p.gam     =    1;    %uM/sec
% % %p.k_gam   =     .27;    %uM
% % %p.gam     =     1;    %uM/sec
% % 
% % p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% % p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% % % p.alpha_0 = 1.0;          
% % % p.Kappa_g = 1.0;
% % 
% % p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% % %p.pulseTimeConstant = @(t) (t>0).*1e20;	% [s]       Average time between pulses
% % p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% % p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% % p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% % p.sigma2        =0; % DO NOT USE OR GHOSTS
% % 
% % p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% % p.Ca_ER_0	= 30;           % [uM]      Initial concentration of Ca2+ in ER
% % p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% % p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% % 
% % p.v_41    = 0;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % 
% % p.k_1     = 0.00036;	% [1/s]     Rate constant of calcium leak from ER
% % % p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% % p.v_40    = 0;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.k_5     = 0;      % [1/s]     Rate constant of calcium extrusion
% % % p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% % % p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% % 
% % p.Ca_media = 0.05;
% % p.P_Ca_media = 0.5;
% % 
% % % Temporary
% % p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% % 
% % 
% % % Flux from media is around 0.5 * (0.02-0.05) = .015
% % % If media Ca is around 10, then P should be around 0.5/333= 1.5e-3
% % 
% % % p.Ca_media = 10;
% % % p.P_Ca_media = 1.5e-4;
% % % p.P_Ca_media = 0;
% 
% %% HOFFER PARAMETERS TO MURDER v_40, v_41, and k_5
% 
% % Hofer Model Parameters
% p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% 
% % Tunable Parameters
% p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% p.k_gam   =     0.27;    %uM
% p.gam     =    1;    %uM/sec
% 
% % Hold still
% p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% p.sigma2        =0; % DO NOT USE OR GHOSTS
% 
% p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% p.Ca_ER_0	= 60;           % [uM]      Initial concentration of Ca2+ in ER
% p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% 
% p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% 
% p.v_41    = 0;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% p.v_40    = 0;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.k_5     = 0;      % [1/s]     Rate constant of calcium extrusion
% 
% p.Ca_media = 0.05;
% p.P_Ca_media = 0.5;
% 
% % Flux from media is around 0.5 * (0.02-0.05) = .015
% % If media Ca is around 10, then P should be around 0.5/333= 1.5e-3
% 
% % p.Ca_media = 10;
% % p.P_Ca_media = 1.5e-4;
% % p.P_Ca_media = 0;
% 
% %% HOFFER PARAMETERS TO MURDER MEDIA CALCIUM
% 
% % % Won't affect oscillation
% % p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% % p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% % p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% % p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% % p.sigma2        =0; % DO NOT USE OR GHOSTS
% % p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% % p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% % p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% % p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% % 
% % % Might affect oscillation
% % p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% % p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% % p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% % p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% % p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% % p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% % p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% % p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% % p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% % 
% % p.k_gam   =     0.27;    %uM
% % p.gam     =    1;    %uM/sec
% % 
% % % Don't touch
% % p.v_41    = 0;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% % p.v_40    = 0;    % [uM/s]    Rate of calcium leak across the plasma membrane
% % p.k_5     = 0;      % [1/s]     Rate constant of calcium extrusion
% % p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% % 
% % p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% % p.Ca_ER_0	= 80;           % [uM]      Initial concentration of Ca2+ in ER
% % p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% % p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% % 
% % p.Ca_media = 0.05;
% % p.P_Ca_media = 0.5;
% % 
% % % Flux from media is around 0.5 * (0.02-0.05) = .015
% % % If media Ca is around 10, then P should be around 0.5/333= 1.5e-3
% % 
% % % p.Ca_media = 10;
% % % p.P_Ca_media = 1.5e-4;
% % % p.P_Ca_media = 0;



%% HOFFER PARAMETERS TO MURDER v_40, v_41, and k_5

% Working with following fluxes:

%     v_rel   = (p.k_1 + p.k_2.*IP3R.*C2.*I2./(p.K_a^2+C2)./(p.K_IP3.^2+I2)).*(CaER-CaC);
%     v_PLC  = p.v_PLC .* normrnd(1,p.sigma, [xSize, ySize]) + p.v_7*C2;
%     v_SERCA = (p.gam .* C2) ./ (p.k_gam.^2 + C2);
%     v_deg   = p.k_9.*IP3;
%     v_media    = p.v_40+p.v_41.*I2./(p.K_r.^2+I2)- p.k_5 .* CaC;
%     dIP3Rdt = p.k_6.*(p.K_i.^2./(p.K_i.^2+C2)-IP3R);

% Hofer Model Parameters
% p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% p.k_gam   = 1e-2;    %uM
% p.gam     = p.k_gam * 3;    %uM/sec
% p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% 
% p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% p.Ca_ER_0	= 30;           % [uM]      Initial concentration of Ca2+ in ER
% p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% 
% p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% p.sigma2        =0; % DO NOT USE OR GHOSTS
% 
% 
% p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% p.v_41    = 0;      % [uM/s]    Maximal rate of activation-dependent calcium influx

%% ATTEMPT TO CONSOLIDATE MEDIA TERM

% Working with following fluxes:

%     v_rel   = (p.k_1 + p.k_2.*IP3R.*C2.*I2./(p.K_a^2+C2)./(p.K_IP3.^2+I2)).*(CaER-CaC);
%     v_PLC  = p.v_PLC .* normrnd(1,p.sigma, [xSize, ySize]) + p.v_7*C2;
%     v_SERCA = (p.gam .* C2) ./ (p.k_gam.^2 + C2);
%     v_deg   = p.k_9.*IP3;
%     v_media    = p.v_40+p.v_41.*I2./(p.K_r.^2+I2)- p.k_5 .* CaC;
%     dIP3Rdt = p.k_6.*(p.K_i.^2./(p.K_i.^2+C2)-IP3R);

% Hofer Model Parameters
% p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
% p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
% p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
% p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
% p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
% p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
% p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
% p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
% p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
% p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
% p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
% p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
% p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLCb
% p.k_gam   = 1e-2;    %uM
% p.gam     = p.k_gam * 3;    %uM/sec
% p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
% p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
% p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
% 
% p.Ca_cyto_0	= 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
% p.Ca_ER_0	= 30;           % [uM]      Initial concentration of Ca2+ in ER
% p.IP3_0     = 0;            % [uM]      Initial concentration of IP3
% p.IP3R_0	= 1;            % [uM]      Initial active concentration of IP3R
% 
% p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
% p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
% p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
% p.sigma     = 0;         % [   ]     Standard deviation of production of IP3
% p.sigma2        =0; % DO NOT USE OR GHOSTS
% 
% % p.Ca_media = 0.05;
% % p.P_Ca_media = 0.5;
% 
% p.Ca_media = 0.10;
% % .015 =  p.P_Ca_media .* (p.Ca_media - CaC)
% p.P_Ca_media = .015 ./ (p.Ca_media - 0.02);
% 
% % Flux from media is around 0.5 * (0.02-0.05) = .015
% % If media Ca is around 10, then P should be around 0.5/333= 1.5e-3
% 
% % p.Ca_media = 10;
% % p.P_Ca_media = 1.5e-4;
% % p.P_Ca_media = 0;
% 
% p.v_40    = 0;    % [uM/s]    Rate of calcium leak across the plasma membrane
% p.k_5     = 0;      % [1/s]     Rate constant of calcium extrusion
% p.v_7     = 0;         % [uM/s]    Maximal rate of PLCd
% p.v_41    = 0;      % [uM/s]    Maximal rate of activation-dependent calcium influx

%% FINAL AICHE MODEL

% Hofer Model Parameters
p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
p.K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
p.beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLC-gamma
p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
p.k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
p.k_gam   = 0.03;    %uM
p.gam     = p.k_gam * 3;    %uM/sec

p.D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
p.D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
p.P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
p.P_Ca    = 0.01*p.P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium

p.Ca_cyto_0         = 0;            % [uM]      Initial concentration of Ca2+ in cytoplasm
p.Ca_ER_0           = 60;           % [uM]      Initial concentration of Ca2+ in ER
p.IP3_0             = .15;            % [uM]      Initial concentration of IP3
p.IP3R_0            = 1;            % [uM]      Initial active concentration of IP3R

% Tuneable
p.v_PLC   = 1e-2;         % [uM/s]    Rate of PLC-gamma
p.noiseMemory = 0;    
p.sigma   = 1;         % [   ]     Standard deviation of production of IP3

p.K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
p.K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
p.K_i     = 0.09;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
p.K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry

p.k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
p.v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
p.k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
p.k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
p.k_gam   = 0.02;    %uM

p.gam     = p.k_gam * 4;    %uM/sec
p.k_6     = 4 / p.K_i * 0.3;        % [1/s]     Rate constant of IP3R inactivation
p.k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     


% Control Parameters
p.pulseTimeConstant = @(t) inf;	% [s]       Average time between pulses
p.IP3Pulse          = 50;    % [uM]      Concentration of IP3 in affected coordinates
p.IP3Extent         = 1;    % [um]      Extent of IP3 pulse
p.sigma2            = 0; % DO NOT USE OR GHOSTS

end

