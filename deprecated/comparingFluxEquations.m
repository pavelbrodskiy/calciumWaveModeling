% This script compares several forms 

k_flux = 8;
k_mu = 0.01;
b = 0.11;
k_1 = 0.7;
beta = 0.15;

banana = @(x)optimization2(x);
[x,fval] = fminsearch(banana,[k_flux, k_mu, b, k_1]);

resolution = 1000;
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

defaultSneydFlux = k_flux.*IP3./(k_mu+IP3).*IP3R.*(b+(1-b).*Ca./(k_1+Ca)) + beta;

subplot(2,2,1)
imshow(defaultHoferFlux,[0,4]);

subplot(2,2,2)
imshow(defaultSneydFlux,[0,4]);

subplot(2,2,3)
imshow(IP3,[]);

subplot(2,2,4)
imshow(Ca,[]);

colormap('jet');