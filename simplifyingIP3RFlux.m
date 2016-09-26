close all

k_flux =.5;
k_mu =0.01; 
b = 0.11;
k_1 =0.70; 
n=1;

k1 = 1;
v1 = 1;
k2 = 1;

IP3_max = 5;
Ca_max = 0.5;

density = 10000;

IP3 = repmat(0:IP3_max/density:IP3_max,[density+1,1]);
Ca = repmat((0:Ca_max/density:Ca_max),[density+1,1])';
IP3R = 1;

banana = @(x)optimization(x);
[x,fval] = fminsearch(banana,[5, 1]);

%x = [2.2, 0.02];

v_in_ablation = k_flux*IP3R.*(IP3.^n./(k_mu.^n+IP3.^n).*(b+(1-b).*Ca)./(k_1+Ca));
v_new = x(1).*IP3R.*Ca.*(IP3./(IP3+x(2)));

subplot(1,2,1)
imshow(imresize(v_in_ablation, [1000 1000], 'nearest'), []);
subplot(1,2,2)
imshow(imresize(v_new, [1000 1000], 'nearest'), []);
%imshow(imresize(IP3, [1000 1000]), [0, IP3_max]);
%imshow(imresize(Ca, [1000 1000]), [0, Ca_max]);
colormap('jet');

% Can assume IP3 has linear effect, Ca2+ has hill effect.