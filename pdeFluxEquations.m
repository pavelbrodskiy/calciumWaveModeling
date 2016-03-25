function hfHandle = pdeFluxEquations(p)
hfHandle = @fHandle;
% ---------------------------------------
function f = fHandle(~, state)
% Calculate squares of C and I to reduce computation
C2 = state.u(1,:).^2;
I2 = state.u(2,:).^2;

% Calculate flux terms to reduce compuatation
J_flux = (k_1+k_2.*state.u(4,:).*C2.*I2) ...
    .*(state.u(3,:)-state.u(1,:)) ...
    ./(K_Ca^2+C2) ...
    ./(K_IP3^2+I2);
J_SERCA = v_SERCA .* C2 ./ (K_SERCA^2 + C2);

% Solve for remaining terms
f(1,:) = J_flux - J_SERCA + v_EC;
f(2,:) = v_gen;
f(3,:) = -beta*(J_flux - J_SERCA);
f(4,:) = k_i*(K_i^2./(K_i^2+C2));
