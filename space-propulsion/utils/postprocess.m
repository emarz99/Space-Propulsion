%% -------------- CEA CODE POSTPROCESS ----------------

%% DATA
clear;clc;close all;

g = 9.81;                                       % Earth Acceleration [m/s2]
m = 250;                                        % Inert mass budget [kg]
deltaV = 2500;                                  % DeltaV [m/s]

% This thrust is decided so that the burning time is about 15min (cf. paper
% of greem propellant)

T = 1000;                                       % Thrust [N]
bt = 15*60;                                     % Burning time [s]

%% Inputs -> Legacy System : N2O4 / MMH (CH3NHNH2)

Pc = 10e5;                                      % Chamber Pressure [Pa]
eps = 200;                                      % Area Ratio [-]
MR = 1.65;                                      % Mixture Ratio [-]

%% Outputs
Isp = 315.5;                                    % Specific Impulse [s]
Temp = [0 0 0 0]';                              % Temperature [K]
P = [0 0 0 0]';                                 % Pressure [Pa]
c_star = 1200;                                  % Characteristic Velocity [m/s]
ct = 1.9;                                       % Thrust Coefficient [-]
eps_real = 0;                                   % Area Ratio [-]

%% Thermodynamic

% Tsiolkovsky
m_init = m*exp(deltaV/(Isp*g));                 % Initial Mass [kg]
m_prop = m_init-m;                              % Propellant mass [kg]
Itot = Isp*m_prop*g;                            % Total Impulse 

c = ct * c_star;                                % Characteristic Velocity [m/s]              
m_dot = T/(Isp*g);                              % Propellant mass flow rate [kg/s]
At = m_dot*c_star/Pc;                           % Throat Area [m2]
Ae = At*eps;                                    % Exit Area [m2]

% The optimal mixture ratio is MR = m_ox/m_fu
m_ox = (m_prop*MR)/(1+MR);
m_fu = m_prop - m_ox; 
% The optimal mixture ratio is MR = m_dot_ox/m_dot_fu
m_dot_ox = (m_dot*MR)/(1+MR);
m_dot_fu = m_dot - m_dot_ox; 

% Chemical properties Fuel and Oxydixer
Mfu = 46.072e-3;                                % Molar mass N2H4 [kg/mol]
dens_fu = 0.8788e3;                             % Density N2H4 298k[kg/m3]
visc_fu = 0.855e-3;                             % Viscosity at 293K [Pa.s]
hf_fu = 54.14;                                  % Enthalpy of formation 298K [kJ/mol]

Mox = 92.016e-3;                                % Molar mass N2O4 [kg/mol]
dens_ox = 1.440e3;                              % Density N2O4 298K[kg/m3]
visc_ox = 0.47e-3;                              % Viscosity at 293K [Pa.s]
hf_ox = -19.56;                                 % Enthalpy of formation 298K[kJ/mol]

% Propellant density [kg/m3]
dens_prop = (m_ox+m_fu)/((m_ox/dens_ox)+(m_fu/dens_fu));  

% Compactness
Iv = Isp * dens_prop;                           % Volumetric Specific Impulse [kg.s/m3]

% Volumes
V_ox = m_ox/dens_ox;                            % Volume of the oxydizer [m3]
V_fu = m_fu/dens_fu;                            % Volume of the fuel [m3]

%% RAO Nozzle Design

% Data
alpha = 15;
de = 2*sqrt(Ae)/sqrt(pi);
dt = 2*sqrt(At)/sqrt(pi);

% % Minimum length approach (mla)
% perc_mla = 0.75;
% L_dc = 0.5*(de-dt)/tand(alpha);
% L_rao = perc_mla * L_dc;
% theta_i = 32;
% theta_e = 10;
% alpha_rao = atand((de-dt)/(2*L_rao));
% lambda = 0.5 * (1 + cosd ((alpha_rao+theta_e)/2));
% 
% thrust_rao = lambda*T;

% Max performance approach (mpa)
perc_mpa = 1;
L_dc = 0.5*(de-dt)/tand(alpha);
L_rao = perc_mpa * L_dc;

% To compute from tables
theta_i = 27;
theta_e = 5;

alpha_rao = atand((de-dt)/(2*L_rao));
lambda = 0.5 * (1 + cosd ((alpha+theta_e)/2));

thrust_rao = lambda*T;

%% Combustion Chamber Design

% Inputs : Pc, L* CR (depends on c* from CEA .out OR fix it within a range)
L_star = 80e-2;                                 % Lstar for N2O2/Hydrazine-base [m]
CR = 3;                                         % Contraction ratio
theta = 45;                                     % Angle CC-throat [deg]
Vcc = L_star * At;                              % Combustion Chamber [m3]

Lcc = (1/CR)*(Vcc/At - At^.5*(CR^(1/3)-1)/(3*pi^.5*tand(theta))); % Length CC [m]
Aside = 2*Lcc*sqrt(At*pi*CR) + At*(CR-1)*(1/sind(theta));    % Side Area [m2]
dside = sqrt(4*Aside/pi);                       % Diameter of side area [m]

%% Injector Plate Design

% Inputs : - Area of C.C (contraction ratio) and minimum diameter of
% injection holes (additive manufacturing precision capability)

deltaPinj = 0.3*Pc;                             % Pressure drop in injector plate [Pa]
Cd = 0.7;                                       % Discharge coefficient [-]
min_diam = 0.2e-3;                              % AM minimum diameter [m]

% I found  Cd = 0.8 on a NASA experimental test paper but 0.7 is the magic number

% TO BE DEFINED
Nox = 30;                                       % Number of orifices 
Nfu = 50;                                       % Number of orifices

% This combination gives a reasonable mass flow rate per injector and a
% diameter of the orifices for ox and fu are both greater than the minimum
% diameter ensured by Additive Manufacturing : min_diam

m_dot_ox1inj = m_dot_ox/Nox;                    % Ox mass flow rate 1 injector [kg/s]
m_dot_fu1inj = m_dot_fu/Nfu;                    % Fu mass flow rate 1 injector [kg/s]

Ainj_ox = m_dot_ox1inj/...
    (Cd*sqrt(2*deltaPinj*dens_ox));             % Ox area inj [m2]
Ainj_fu = m_dot_fu1inj/...
(Cd*sqrt(2*deltaPinj*dens_fu));                 % Fu area inj [m2]

r_inj_ox = sqrt(Ainj_ox/pi);                    % Ox radius inj [m]
r_inj_fu = sqrt(Ainj_fu/pi);                    % Fu radius inj [m]

printf("The dimension of the orifices is %f")

v_inj_ox = Cd*sqrt(2*deltaPinj/dens_ox);        % Velocity inj ox [m/s] 
v_inj_fu = Cd*sqrt(2*deltaPinj/dens_fu);        % Velocity inj fu [m/s]

