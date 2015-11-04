%function [maxA, maxC, whm, shme, shmr, Anisotropy] = simulation(k_flux, b, k1, k2, tau_h, k_gam, gam, beta, k_mu, V_p, k_p, D_c, D_p, P_p, P_c, minF, maxF, minL, maxL, name)
clear all
close all

tic

k_flux  =     8.00;    % 3 uM/sec        % affects intensity of Ca2+ flash %TUNED
b       =     0.11;    %dim-less
k1      =     0.70;    %uM
k2      =     0.70;    %uM
tau_h   =     0.20;     %0.20;    %sec or 60   % affects size of flash
k_gam   =     0.27;    %uM
gam     =     1.1*1.00;    %uM/sec
beta    =     0.15;    %uM/sec
k_mu    =     0.01;    %uM
V_p     =     0.08;    %1/sec %TUNED
k_p     =     0.5*1.00;    %uM

D_c     =     20.0;    %uM^2/sec
D_p     =     300.;    %uM^2/sec

P_p     =     (1/10)*2.0;    %uM/sec   % affects dynamics of expansion /dying out %TUNED
P_c     =     (0.05*P_p);   %uM/sec %TUNED

minF = 0;
maxF = 0;
minL = 4;
maxL = 4;

%Define the spatial grid
xdir = 100; %um
ydir = xdir; %um

%Spatial grid step
dspace = 2; %um
xnodes = xdir/dspace; %grid nodes x 1px=dspace um
ynodes = ydir/dspace; %grid nodes y

u = ones(xnodes,ynodes); %for when we need it.

%Simulation time
time = 500;

%Time step
dtime = 0.01;
timesteps = time/dtime;

Css     =     sqrt((beta*k_gam^2)/(gam-beta)); %SS calcium value

if minL == maxL
    sizemat = u * maxL;
else
    sizemat = repmat((minL:((maxL-minL)/(dspace-1)):maxL),dspace,1);
end

if minF == maxF
    fmat = u * maxF;
else
    fmat = repmat((minf:((maxf-minf)/(dspace-1)):maxf),dspace,1);
end




%effective diffusivities
%x-elongation
Deffpmatxx = (1/D_p + 1./((sizemat*P_p))).^(-1); %diffusion in x direction
Deffcmatxx = (1/D_c + 1./((sizemat*P_c))).^(-1);
Deffpmatyy = (1/D_p + 1./((((1-fmat)./(1+fmat)).*sizemat)*P_p)).^(-1); % diffusion in y direction
Deffcmatyy = (1/D_c + 1./((((1-fmat)./(1+fmat)).*sizemat)*P_p)).^(-1);

%flash stimulus
stimtime                     = 0.2;        %sec
stimpulse                    = 0.01;         %uM
stimradius                   = 2.0; %um

stimMat = 0*u;
for i = 1:xnodes
    for j = 1:ynodes
        if (i-xnodes/2)^2+(j-ynodes/2)^2 <= (stimradius/dspace)^2
            stimMat(i,j) = stimpulse;
        end
    end
end

%loop parameters
n    = 0; %time loop
q    = 0; %movie frames
done = 0; %flag stop

%Initialize variables on grid
p=0.*u;
c=Css*u;
h=u;


%Define .tif stack structure
stack = zeros(xnodes, ynodes, timesteps/1000, 'uint16'); %preallocate image space for tiff stack, every 1/4 second

while ~done %time loop
    
    %initiate stimulus
    if n*dtime <= stimtime
        p = p + stimMat;
    end
    
    %calculate fluxes and compensate cell type
    mu_p    = p ./ (k_mu + p);
    J_flux  = k_flux .* mu_p .* h .* (b + ((1 - b) .* c) ./ (k1 + c));
    J_pump  = (gam .* c.^2) ./ (k_gam.^2 + c.^2);
    J_leak  = beta;
   
    %create padded matrices to incorporate neumann bcs
    pp = [[0 p(2,:) 0];[p(:,2) p p(:,end-1)];[0 p(end-1,:) 0]];
    cc = [[0 c(2,:) 0];[c(:,2) c c(:,end-1)];[0 c(end-1,:) 0]];

    %spatial derivatives approx
    pxx = (pp(2:end-1,1:end-2) + pp(2:end-1,3:end) -2*p)/dspace; 
    pyy = (pp(1:end-2,2:end-1) + pp(3:end,2:end-1) -2*p)/dspace;
    
    cxx = (cc(2:end-1,1:end-2) + cc(2:end-1,3:end) -2*c)/dspace; 
    cyy = (cc(1:end-2,2:end-1) + cc(3:end,2:end-1) -2*c)/dspace;
    
    %update p,c,h
    dpdt = Deffpmatyy .* (pyy) +Deffpmatxx.* (pxx)-((V_p * p * k_p) ./ (k_p + p));
    dcdt = Deffcmatyy .* (cyy) +Deffcmatxx.* (cxx) + J_flux - J_pump + J_leak;
    dhdt = k2^2 ./ (tau_h * (k2^2 + c^2)) - h ./ tau_h;
    
    p_new = p + dpdt*dtime;
    c_new = c + dcdt*dtime;
    h_new = h + dhdt*dtime;
    
    p = p_new; clear p_new
    c = c_new; clear c_new
    h = h_new; clear h_new
    
        if rem(n,5) == 0
            q=q+1;
            stack(:,:,q)=100000*c;
            cstack(:,:,q)=double(c);
            imshow(c,[]);
            drawnow
            n*dtime
        end
    
    n=n+1;
    done=(n > timesteps);
    %if max(c(:)) < 1.0e-4, done=1; end             % If activation extinguishes, quit early.
    if ~isempty(get(gcf,'userdata')), done=1; end   % Quit if user clicks on 'Quit' button.
end


[a1, a2] = uiputfile('*.TIF','Save Calcium Video');
imwrite(stack(:,:,1), [a1, a2])
    for k = 2:size(stack,3)
        imwrite(stack(:,:,k), name, 'writemode', 'append');
    end
[maxA, maxC, whm, shme, shmr, Anisotropy] = analyzeTemporalData(cstack, dtime*500, dspace);

toc
close all

