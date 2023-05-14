%{

Set engine fixed parameters

%}

%% constants
engine.g0 = 9.81; 
engine.T0 = 288.15; 

engine.nSub = 2;
engine.nSup = 1; 


% Injectors
engine.CDox = 0.7;        % oxidiser injection plate CD 
engine.CDfu = 0.7;        % fuel injection plate CD
engine.AinjMin = (0.0008^2)*pi/4;    % Minimum injection area [m2]
engine.nInjMax = 60;      % Number of injectors max
engine.injOFmult = nan;     % number of Ox Injector for each Fu Injector (nan mean Like on Like injection)
engine.pLossInj = 0.3;    % pressure loss through the injection plate [%]

% Pipelines
engine.pipeDiamF = 0.00635;    % Fuel pipelines diameter [m]  (1/4 in)
engine.pipeDiamO = 0.00635;    % Oxidizer pipelines diameter [m] (1/4 in)

% cooling
engine.oxChamber = 0.05;    % fraction of chamber cooled by oxidizer [-]
engine.oxConv = 0;         % fraction of convergent cooled by oxidizer [-]
engine.fuChamber = 0.95;    % fraction of chamber cooled by fuel [-]
engine.fuConv = 1;         % fraction of convergent cooled by fuel [-]
engine.Tcoat = 1300;           % Coating max wall temperature [K]
engine.kCoat = nan;          % Coating conductive coefficient [W/m K]

%tank 220 PSIA carburante  150 PSIA ox   
engine.pHelium = 235e5;   % Helium tank pressure [Pa] 
engine.TShelum = 298.15;  % Helium storage temperature [K]
engine.MMhelium = 4.00260e-3;  % Helium molar mass [kg/mol] 

%% Mission param
engine.ms = 250;         % structural mass [kg]
engine.deltaV = 2500;    % required deltaV [m/s]

engine.pc = 10e5;           % Chamber Pressure [Pa]
engine.eps = 130;           % Expansion Ratio [-]
engine.OF = 0.9347;            % Mixture Ratio (oxidizer/fuel) [-]
engine.Lstar = 0.75;        % Chamber L* [m]
engine.Thrust = 700;        % Thrust [N]  

engine.CR = 8;                 % Chamber contraction ratio [Ac/At]
engine.Twall = 900;            % Chamber wall temperature [K]
engine.TSf = 288;              % Fuel storage temperature [K]
engine.TSo = 288;              % Oxidizer storage temperature [K]
engine.thetac = 45*pi/180;     % convergent angle [rad]

% fuel data
engine.rhoF = 1008;             % Fuel density [km/m3]
engine.muF = 0.876*(1e-3);      % Fuel viscosity [Pa*s]
engine.TboilF = 460;            % Boiling temperature [K]
engine.CPf = 3430;              % Fuel CP [J/kg*K] 

% Ox data

Hvap = 38.12;   % Enthalpy of vaporization kJ/mol
R = 0.00831446261; % Gas constant [kJ/mol K]
T1 = 293;  
p1 = 1e5; 
p2 = engine.pc*(1 + engine.pLossInj); 

engine.TboilO = 1/(1/T1 + R/Hvap * log(p1/p2)); 

engine.rhoO = 1450;       % Oxidizer density [km/m3]   
engine.muO = 0.47*(1e-3);        % Oxidizer viscosity [Pa*s]
% engine.TboilO = 293;      % Oxidizer boiling temperature [K]
engine.CPo = 1548.7;      % [J/kg g]


%% Create CEA input - Toxix
input.problem.type = 'rocket'; 
input.problem.solution = 'frozen'; 
input.problem.facFlag = false; 
input.problem.acat = 0; 
input.problem.nfz = 2; 

input.fuel.name = {'N2H4'};
input.fuel.T = 298.15; 
input.fuel.wt = nan; 

input.oxidizer.name = {'N2O4'}; 
input.oxidizer.T = 298.15; 
input.oxidizer.wt = nan; 

input.OF = [engine.OF]; 
input.pc = [engine.pc]./(1e5); 

input.subar = []; 
input.supar = linspace(1,engine.eps,engine.nSup+1);
input.subar = input.supar(2:end);

input.coeffOut = 'p t mach...fz son rho gam visc cond...fz cp pran...fz o/f ivac..fz isp...fz cf...fz';

engine.input = input; 


engine.ceaPath = 'CEA/'; 

%% Nozzle Data
load('data\RAOfinal.mat');
load('data\RAOinital.mat'); 
engine.RAOfinal =  RAOfinal; 
engine.RAOinitial = RAOinitial; 


%% Optimization

lb(1) = 0.9;                  % OF 
%lb(2) = 100;                  % T

ub(1) = 1.8;                  % Ignition temperature [K]
%ub(2) = 750;                 % OF 

engine.wTemp = 10;
engine.wThrust = 5;

