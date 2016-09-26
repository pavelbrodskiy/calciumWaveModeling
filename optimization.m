function [ SSR ] = optimization( x )
IP3_max = 0.1;
Ca_max = 0.2;
density = 1000;
k_flux =.5;
k_mu =0.01; 
b = 0.11;
k_1 =0.70; 
n=1;
IP3 = repmat(0:IP3_max/density:IP3_max,[density+1,1]);
Ca = repmat((0:Ca_max/density:Ca_max),[density+1,1])';
IP3R = 1;

v_in_ablation = k_flux*IP3R.*(IP3.^n./(k_mu.^n+IP3.^n).*(b+(1-b).*Ca)./(k_1+Ca));
v_new = x(1).*IP3R.*Ca.*(IP3./(IP3+x(2)));

error=(v_in_ablation-v_new).^2;
SSR = sum(error(:));