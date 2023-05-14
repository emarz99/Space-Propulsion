%% CONFIG

N = 300; 
CDoxLim = [0.6, 0.8]; 
CDFuLim = [0.6, 0.8]; 

DinjF0 = sqrt(engineFinal.AinjF/pi)*2; 
DinjO0 = sqrt(engineFinal.AinjO/pi)*2;

DinjFLim = [DinjF0,50*(1e-6)]; 
DinjOLim = [DinjO0,50*(1e-6)];

AtLim = [1, 0.01].*engineFinal.At;
epsLim = [1, 0.01].*engineFinal.eps; 

getValUnif = @(lim, N) lim(1).*ones(1, N) + (lim(2) - lim(1)).*rand(1, N);  
getValNorm = @(lim, N) lim(1).*ones(1, N) + lim(2).*randn(1, N);    % lim(1): mu, lim(2): STD



%% MAIN
CDox = getValUnif(CDoxLim, N);
CDfu = getValUnif(CDFuLim, N);

DinjF = getValNorm(DinjFLim, N);
DinjO = getValNorm(DinjOLim, N);

AinjF = DinjF.^2*pi/4; 
AinjO = DinjO.^2*pi/4;

At = getValNorm(AtLim, N);
eps = getValNorm(epsLim, N);


data_engine = cell(1, N);  

for i = 1:N
    engineN = engineFinal; 
    
    engineN.CDox = CDox(i); 
    engineN.CDfu = CDfu(i); 
    
    engineN.AinjF = AinjF(i); 
    engineN.AinjO = AinjO(i); 

    engineN.At = At(i); 
    engineN.eps = eps(i);

    
    data_engine{i} = engineSim(engineN);

    clc; 
    fprintf('%d/%d\n', i, N)

end



%%
close all; 

T = zeros(1, N); 
Isp = zeros(1, N); 
dV = zeros(1, N); 
OF = zeros(1, N); 
pc = zeros(1, N); 
TO = zeros(1, N); 
TF = zeros(1, N); 

for i = 1:N
    de = data_engine{i}; 
    T(i) = de.Thrust; 
    Isp(i) = de.Isp; 
    dV(i) = de.deltaV; 
    OF(i) = de.OF; 
    pc(i) = de.pc; 
    TO(i) = de.TO; 
    TF(i) = de.TF; 
end

figure; 
histogram(T, 15); 
title('Thrust'); 
xlabel('Thrust [N]'); 

figure; 
histogram(Isp, 15); 
title('Specific impulse'); 
xlabel('Isp [s]'); 

figure; 
histogram(dV, 15); 
title('Delta V'); 
xlabel('deltaV [m/s]'); 

figure; 
histogram(OF, 15); 
title('Oxidizer/Fuel ratio'); 
xlabel('OF [-]'); 

figure; 
histogram(pc, 15); 
title('Chamber pressure'); 
xlabel('pc [bar]'); 

figure; 
histogram(TO, 15); 
title('Oxidizer max temperature'); 
xlabel('TOX [K]'); 

figure; 
histogram(TF, 15); 
title('Fuel max temperature'); 
xlabel('TFU [K]'); 






