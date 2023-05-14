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
engine.injOFmult = 2;     % number of Ox Injector for each Fu Injector (nan mean Like on Like injection)
engine.pLossInj = 0.3;    % pressure loss through the injection plate [%]

% Pipelines
engine.pipeDiamF = 0.00635;    % Fuel pipelines diameter [m]  (1/4 in)
engine.pipeDiamO = 0.00635;    % Oxidizer pipelines diameter [m] (1/4 in)

% cooling
engine.oxChamber = 1;      % fraction of chamber cooled by oxidizer [-]
engine.oxConv = 0.5;       % fraction of convergent cooled by oxidizer [-]
engine.fuChamber = 0;      % fraction of chamber cooled by fuel [-]
engine.fuConv = 0.5;       % fraction of convergent cooled by fuel [-]

engine.pHelium = 235e5;   % Helium tank pressure [Pa] 
engine.TShelum = 298.15;  % Helium storage temperature [K]
engine.MMhelium = 4.00260e-3;  % Helium molar mass [kg/mol] 


%% Mission param
engine.ms = 250;         % structural mass [kg]
engine.deltaV = 2500;    % required deltaV [m/s]

engine.pc = 10e5;            % Chamber Pressure [Pa]
engine.eps = 130;            % Expansion Ratio [-]
engine.OF = 5.0916;           % Mixture Ratio (oxidizer/fuel) [-]
engine.Lstar = 1.52;         % Chamber L* [m]  %For green the range is between 152-178
engine.Thrust = 1400;        % Thrust [N]  

engine.CR = 8;                 % Chamber contraction ratio [Ac/At]
engine.Twall = 900;            % Chamber wall temperature [K]
engine.Tcoat = 1300;           % Coating max wall temperature [K]
engine.kCoat = 0.5;           % Coating conductive coefficient [W/m K]
engine.TSf = 293;              % Fuel storage temperature [K]
engine.TSo = 293;              % Oxidizer storage temperature [K]
engine.thetac = 45*pi/180;     % convergent angle [rad]

% fuel data
engine.rhoF = 800;                           % Fuel density [kg/m3]
engine.muF = 0.75*(1e-3);                    % Fuel viscosity [Pa*s]
T1F= 460;                                    % Boiling temperature at 1 bar[K]
engine.CPf = 1882;                           % Fuel CP [J/kg*K] 
R=8.314; p1=10^5; p2=20*10^5; H=246000;      %J/mol
engine.TboilF = 1/(1/T1F-(R/H)*log(p2/p1));  % Clapeyron: boiling temperature at p2

% Ox data
engine.rhoO = 1437;                           % Oxidizer density [kg/m3]   
engine.muO = 1.245*(1e-3);                    % Oxidizer viscosity [Pa*s]
T1O = 150+273.15;                             % Oxidizer boiling temperature at 1 bar [K]
engine.CPo = 2619;                            %4180; % [J/kg K]
R=8.314; p1=10^5; p2=20*10^5; H=48500;        %J/mol
engine.TboilO = 1/(1/T1O-(R/H)*log(p2/p1));   % Clapeyron: boiling temperature at p2


%% Create CEA input - Toxix
input.problem.type = 'rocket'; 
input.problem.solution = 'frozen'; 
input.problem.facFlag = false; 
input.problem.acat = 0; 
input.problem.nfz = 2; 

input.fuel.name = {'RP-1'};
input.fuel.T = 298.15; 
input.fuel.wt = nan; 

input.oxidizer.name = {'H2O','O2'};         %{'H2O','O2'}; 
input.oxidizer.T = [1200,1200];              %[298.15, 298.15]; 
input.oxidizer.wt = 100*[9/17, 8/17];             %100*[9/17, 8/17]; 

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

lb(1) = 5;                  % OF 
%lb(2) = 100;                  % T

ub(1) = 9;                  % Ignition temperature [K]
%ub(2) = 750;                 % OF 

engine.wTemp = 10;
engine.wThrust = 5;

