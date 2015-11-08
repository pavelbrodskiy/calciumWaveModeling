
beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER
v_PLC   = 1e-2;         % [uM/s]    Rate of PLC-gamma
gammaShape = 1;
K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
K_IP3	  = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
K_i     = 0.09;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
k_gam   = 0.02;    %uM
gam     = k_gam * 4;    %uM/sec
k_6     = 4 / K_i * 0.3;        % [1/s]     Rate constant of IP3R inactivation
k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
Ca_media = 0.035;
P_Ca_media = 0.51;
D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
P_Ca    = 0.01*P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium


runNumber = 80;
percentVariation = 1;
timeAtStart = now();

sweeper2( D_IP3, '$PD_IP3', percentVariation, runNumber, timeAtStart);
sweeper2( P_IP3, '$PP_IP3', percentVariation, runNumber, timeAtStart);
sweeper2( P_Ca, '$PP_Ca', percentVariation, runNumber, timeAtStart);
sweeper2( D_Ca, '$PD_Ca', percentVariation, runNumber, timeAtStart);
sweeper2( k_9, '$Pk_9', percentVariation, runNumber, timeAtStart);
sweeper2( k_2, '$Pk_2', percentVariation, runNumber, timeAtStart);
sweeper2( k_gam, '$Pk_gam', percentVariation, runNumber, timeAtStart);
sweeper2( gam, '$Pgam', percentVariation, runNumber, timeAtStart);
sweeper2( k_6, '$Pk_6', percentVariation, runNumber, timeAtStart);
sweeper2( Ca_media, '$PCa_media', percentVariation, runNumber, timeAtStart);
sweeper2( P_Ca_media, '$PP_Ca_media', percentVariation, runNumber, timeAtStart);
sweeper2( beta, '$Pbeta', percentVariation, runNumber, timeAtStart);
sweeper2( v_PLC, '$Pv_PLC', percentVariation, runNumber, timeAtStart);
sweeper2( gammaShape, '$PgammaShape', percentVariation, runNumber, timeAtStart);
sweeper2( K_a, '$PK_a', percentVariation, runNumber, timeAtStart);
sweeper2( K_IP3, '$PK_IP3', percentVariation, runNumber, timeAtStart);
sweeper2( K_i, '$PK_i', percentVariation, runNumber, timeAtStart);
sweeper2( k_1, '$Pk_1', percentVariation, runNumber, timeAtStart);
