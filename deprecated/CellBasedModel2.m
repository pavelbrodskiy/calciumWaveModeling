clear all

%Variables
kflux = 8 ; % (microM/s)
b = 0.11;
k1 = 0.7; %(microM)
k2 = 0.7; %(microM)
Th = 0.2 ; %(s)
kg = 0.27; %(microM)
g = 1; % (MicroM/s)
B = 0.15; %(microM/s)
kmu = 0.01; % (MicroM)
Vp = 0.08; % (1/s)
kp = .5;  % (MicroM)
Dc = 20;  % (Microm^2/s)
Dp = 300;  % (Microm^2/s)

pulseLength = 0.2; %S
pulseA = 1e-4;


Fp = .2; %IP3 permeability
Fc = Fp/20;

n = 1; %(Hill coefficient)  

TMax=500;
DeltaT=0.01;

stimulus = pulseLength/DeltaT;

DeltaX = .5;
NC = 128;
L = 4;

sizeOfDomain = (L*NC/DeltaX);

C = zeros(1, sizeOfDomain);
bound = C;
P = C;
h = ones(1,sizeOfDomain);
P(int16(sizeOfDomain/2)) = 1;
P1 = P;
C1 = ones(1, sizeOfDomain)*.1;
C = C1;
h1 = h;

found = false;
for i = 12:L/DeltaX:sizeOfDomain
   bound(int16(i)) = 1; 
   if ~found && i > sizeOfDomain/2
      found = true;
      
   end
end

boundleft = logical(bound);
boundleftleft = logical(boundleft);
boundright = logical(boundleft);
boundleftleft(1) = [];
boundleftleft = [boundleftleft 0];
boundright = [0 boundright];
boundright(end) = [];
boundrightright = logical([0 boundright]);
boundrightright(end) = [];
boundleftleft = logical(boundleftleft);
boundright = logical(boundright);
doublebound = boundright | boundleft;
i = 1;
xs = DeltaX:DeltaX:(L*NC);


Dc1 = (1/Dc + 1/(Fc*L))^(-1);
Dp1 = (1/Dp + 1/(Fp*L))^(-1);

domainDim = size(P);

sigma = 0;
K_Ca    = 0.3;	% [uM]      Half-saturation constant for calcium activation of PLCd
v_7     = 0.003;	% [uM/s]    Maximal rate of PLCd

for T = 0:DeltaT:TMax
    
    if T < stimulus
        P(int16(sizeOfDomain/2)) = P(int16(sizeOfDomain/2))+pulseA;
        P1(int16(sizeOfDomain/2)) = P1(int16(sizeOfDomain/2))+pulseA;
    end
    
    
%     mu=(P.^n).*((((kmu)^n)+P.^n).^(-1));
%     Jflux= kflux*mu.*h.*((b+(1-b)*C.*((k1+C).^(-1))));
%     Jpump=g*(C.^2).*(((kg)^2 + C.^2).^(-1));
%     J1=((k2)^2)*(((k2)^2 + C.^2).^(-1)) - h;
%     
%     P = P + DeltaT * ((Dp*del2(P)).*double(~doublebound)/DeltaX^2 - (Fp*gradient(P,DeltaX).*double(boundright)) + (Fp*gradient(P,DeltaX).*double(boundleft)) - (Vp*P*kp).*((kp + P).^(-1)));
%     C = C + DeltaT * ((Dc*del2(C)).*double(~doublebound)/DeltaX^2 - (Fc*gradient(C,DeltaX).*double(boundright)) + (Fc*gradient(C,DeltaX).*double(boundleft)) + Jflux - Jpump + B);
%     h = h + DeltaT * ((1/Th)*J1);

    P1 = enforceNoFluxCondition(P1);
    C1 = enforceNoFluxCondition(C1);
    h = enforceNoFluxCondition(h);

    PLC_gamma = (v_7*C1.^2)./(K_Ca^2+C1.^2);

    mu=(P1.^n).*((((kmu)^n)+P1.^n).^(-1));
    Jflux= kflux*mu.*h1.*((b+(1-b)*C1.*((k1+C1).^(-1))));
    Jpump=g*(C1.^2).*(((kg)^2 + C1.^2).^(-1));
    J1=((k2)^2)*(((k2)^2 + C1.^2).^(-1)) - h1;
    
    P1 = P1 + DeltaT * ((Dp1*del2(P1))/DeltaX^2 - (Vp*P1*kp).*((kp + P1).^(-1)) + PLC_gamma);
    C1 = C1 + DeltaT * ((Dc1*del2(C1))/DeltaX^2 + Jflux - Jpump + B);
    h1 = h1 + DeltaT * ((1/Th)*J1);

    P1 = P1 .* normrnd(1, sigma, domainDim);
    P1(P1<0) = 0;
    
%     P(boundleftleft) = (P(boundleft) + P(boundleftleft))/2;
%     P(boundleft) = P(boundleftleft);
%     P(boundrightright) = (P(boundright) + P(boundrightright))/2;
%     P(boundright) = P(boundrightright);
%     
%     C(boundleftleft) = (C(boundleft) + C(boundleftleft))/2;
%     C(boundleft) = C(boundleftleft);
%     C(boundrightright) = (C(boundright) + C(boundrightright))/2;
%     C(boundright) = C(boundrightright);
    
    
    if mod(T,DeltaT*100) == 0
        
        subplot(1,3,1)
        plot(xs,P1);%,xs,P);
        axis([0,(L*NC),0,2])
        subplot(1,3,2)
        plot(xs,C1);%,xs,C);
        axis([0,(L*NC),0,2])
        subplot(1,3,3)
        plot(xs,h1);%,xs,C);
        axis([0,(L*NC),0,1])
        T
        
        drawnow
        
%         PStorewc(i,:) = P;
%     CStorewc(i,:) = C;
    PStorenc(i,:) = P1;
    CStorenc(i,:) = C1;
    i = i + 1;
    
    
    end
    %imshow(imresizeOfDomain(PStore,[i i]),[])
    
    
end

% figure(1);
% plot(P);
% hold on
% plot(C);
