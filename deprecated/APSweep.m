
k_1     = 0.0004;       % [1/s]     Rate constant of calcium leak from ER
k_2     = 0.08;         % [1/s]     Rate constant of calcium release through IP3
v_40    = 0.025;        % [uM/s]    Rate of calcium leak across the plasma membrane
v_41    = 0.2;          % [uM/s]    Maximal rate of activation-dependent calcium influx
k_5     = 0.5;          % [1/s]     Rate constant of calcium extrusion
k_6     = 4;            % [1/s]     Rate constant of IP3R inactivation
k_9     = 0.08;         % [1/s]     Rate constant of IP3 degradation
K_IP3	= 0.3;          % [uM]      Half-saturation constant for IP3 activation of IP3R
K_a     = 0.2;          % [uM]      Half-saturation constant for calcium activation of IP3R
K_i     = 0.2;          % [uM]      Half-saturation constant for calcium inhibition of IP3R
K_r     = 1;            % [uM]      Half-saturation constant for agonist-dependent calcium entry
D_IP3   = 280;          % [uM/s^2]  Diffusion coefficient of IP3
D_Ca    = 20;           % [uM/s^2]	Diffusion coefficient of calcium
beta    = 20;           %           Ratio of the effective volumes for calcium of cytoplasm and ER
V_CaInd = 8e-3;
V_Ca    = 8e-2;
K_Ca    = 3e-1;
k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
noiseSigma  = 0.3;      % [1/s] Noise per second

runNumber = 100;
percentVariation = 1;
timeAtStart = now();

sweeper( k_1, 'k_1', percentVariation, runNumber, timeAtStart);
sweeper( k_2, 'k_2', percentVariation, runNumber, timeAtStart);
sweeper( v_40, 'v_40', percentVariation, runNumber, timeAtStart);
sweeper( v_41, 'v_41', percentVariation, runNumber, timeAtStart);
sweeper( k_5, 'k_5', percentVariation, runNumber, timeAtStart);
sweeper( k_6, 'k_6', percentVariation, runNumber, timeAtStart );
sweeper( k_9, 'k_9', percentVariation, runNumber, timeAtStart);
sweeper( K_IP3, 'K_IP3', percentVariation, runNumber, timeAtStart);
sweeper( K_a, 'K_a', percentVariation, runNumber, timeAtStart);
sweeper( K_i, 'K_i', percentVariation, runNumber, timeAtStart);
sweeper( K_r, 'K_r', percentVariation, runNumber, timeAtStart);
sweeper( D_IP3, 'D_IP3', percentVariation, runNumber, timeAtStart);
sweeper( D_Ca, 'D_Ca', percentVariation, runNumber, timeAtStart);
sweeper( beta, 'beta', percentVariation, runNumber, timeAtStart);
sweeper( V_CaInd, 'V_CaInd', percentVariation, runNumber, timeAtStart);
sweeper( V_Ca, 'V_Ca', percentVariation, runNumber, timeAtStart);
sweeper( K_Ca, 'K_Ca', percentVariation, runNumber, timeAtStart);
sweeper( k_3, 'k_3', percentVariation, runNumber, timeAtStart);
sweeper( P_IP3, 'P_IP3', percentVariation, runNumber, timeAtStart);
sweeper( noiseSigma, 'noiseSigma', percentVariation, runNumber, timeAtStart);


