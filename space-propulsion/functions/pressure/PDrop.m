function [ PDropTotal, Lenght ] = PDrop(Cooling, PCombustionChamber, MassFlowRate, FlowVelocity, Rho, Mu)
% set length in base of maximum acceptable delta p
DeltaPressurePipe = 20e3 ;

% Radius of feeding Pipes
Radius = RadiusPipe(MassFlowRate, FlowVelocity, Rho) ;

% Reynolds Number
Re = Reynolds(Rho, FlowVelocity, 2.*Radius, Mu) ;

% Darcy Friction Factor
K = 0.002 ; % Absolute Roughness Coefficient for Alluminum [mm]
Epsilon = K/(2.*Radius) ; % Relative Roughness Coefficient for Alluminum Oxydizer [mm]

%Darcy = friction(Re, Epsilon, 'iterative') ;
Darcy = 64/Re; % Laminar flow

% MAX Lenght of feeding Pipes
Lenght = (2.*DeltaPressurePipe.*2.*Radius)/(Darcy.*Rho.*(FlowVelocity.^2));

% Pressure Drop due to tank exit
PDropTank = PressureDropTank(Rho, FlowVelocity) ;

% Pressure Drop due to feed lines pipes
%PDropFeedLine = PressureDropFeedLine(Darcy, Rho, FlowVelocity, Lenght, 2.*Radius) ;  

% Pressure Drop due to Injection Head and Injectors
PDropInjectors = PCombustionChamber.*0.30 ; 

% Pressure Drop due to Cooling Jacket
if Cooling == 0                 % Liquid not used for cooling 
    PDropCoolingJacket = 0 ; 
elseif Cooling == 1             % Cooling used for fuel
    PDropCoolingJacket = PCombustionChamber.*0.25 ; 
end

PDropValves = 1.25e5;         % 0.15 for each valve, 2 valve (control + check)


PDropTotal = PDropTank + DeltaPressurePipe + PDropCoolingJacket + PDropInjectors + PDropValves;
