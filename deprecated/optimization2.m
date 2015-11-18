function [ SSR ] = optimization2( x )

resolution = 10;
IP3 = repmat(linspace(0,.5,resolution),[resolution,1]);
Ca = repmat(linspace(0,.5,resolution),[resolution,1])';

CaER = 80;
IP3R = 1;

k_1 = 4e-4;
k_2 = 8e-2;
K_a = 0.2;
K_IP3 = 0.3;

defaultHoferFlux = (k_1+k_2*IP3R.*Ca.^2.*IP3.^2./(K_a^2+Ca.^2)./(K_IP3.^2+IP3.^2)).*(CaER-Ca);

k_flux = x(1);
k_mu = x(2);
b = x(3);
k_1 = x(4);
beta = 0.15;

defaultSneydFlux = k_flux.*IP3./(k_mu+IP3).*IP3R.*(b+(1-b).*Ca./(k_1+Ca)) + beta;

SSR = sum(sum((defaultHoferFlux-defaultSneydFlux).^2));