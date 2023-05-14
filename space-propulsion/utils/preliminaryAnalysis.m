%% PROPULSION PROJECT (15/04/2023 - 23/04/2023)
% Design of a set of engines for inspace application

%% DATA
clear;clc;close all;

g = 9.81;                           % Earth Acceleration [m/s2]
m = 250;                            % Inert mass budget [kg]
deltaV = 2500;                      % DeltaV [m/s]

%% ------------------- TOXIC ------------------------

%% 1 -> Toxic Couple : N2H4 / N2O4
clc

% Tables values 
MR1 = 1.08;                         % Mixture ratio [-]
% In this value is embedded, the molar charachteristics of each propellant
% n = m/M
c_star1 = 1765;                     % c* coefficient [m/s]
k1 = 1.26;                          % Specific Heat Ratio [-]
Is1 = 336;                          % Specific Impulse [s]

% Specific Impulse (Sutton and Biblarz for frozen equilibrium)
% Is1 = 283;  

Tc1 = 3258;                         % Chamber temperature [K]
Pc1 = 19.25e5;                      % Chamber pressure [Pa]
Thrust1 = 450;                      % Thrust [N]

c1 = Is1*g;                         % Characteristic Velocity [m/s]
ct1 = c1/c_star1;                   % Thrust Coefficient [-]
m_dot1 = Thrust1/c1;                % Mass Flow Rate [kg/s]
At1 = m_dot1*c_star1/Pc1;           % Throat Area [m2]

% Chemical properties
Mfu1 = 32.0452e-3;                  % Molar mass N2H4 [kg/mol]
dens_fu1 = 1.005e3;                 % Density N2H4 [kg/m3]
visc_fu1 = 0.97e-3;                 % Viscosity at 293K [Pa.s]
Mox1 = 92.016e-3;                   % Molar mass N2O4 [kg/mol]
dens_ox1 = 1.440e3;                 % Density N2O4 [kg/m3]
visc_ox1 = 0.47e-3;                 % Viscosity at 293K [Pa.s]

% Tsiolkovsky
m_init = m*exp(deltaV/(Is1*g));     % Initial Mass [kg]
m_prop1 = m_init-m;                 % Propellant mass [kg]
Itot1 = Is1*m_prop1*g;              % Total Impulse 

% The optimal mixture ratio is MR = m_ox1/m_fu1
m_ox1 = (m_prop1*MR1)/(1+MR1);
m_fu1 = m_prop1 - m_ox1; 

% Propellant density [kg/m3]
dens_prop1 = (m_ox1+m_fu1)/((m_ox1/dens_ox1)+(m_fu1/dens_fu1));  

% Compactness
Iv1 = Is1 * dens_prop1;             % Volumetric Specific Impulse [kg.s/m3]

% Volumes
V_ox1 = m_ox1/dens_ox1;             % Volume of the oxydizer [m3]
V_fu1 = m_fu1/dens_fu1;             % Volume of the fuel [m3]
    
% We are in space, we can't be optimal in vacuum. Therefore, we will need
% to impose epsilon and after that get the exit area. The value must be for
% inspace systems between 40 and 400. The exit mach must be between 5 and
% 7. We suppose eps = 40 from the plot of "Data for nozzle discharging in
% VACUUM" in "Useful Tables".

% Exit Area
eps1 = 40;                          % Area ratio [-]
Ae1 = eps1*At1;                     % Exit Area [m2]

%% Pressure-fed architecture

% From plot relating tank volume and tank pressure
Ptank_ox = 25e5;                    % Oxydizer Tank Pressure [Pa]
Ptank_fu = 23e5;                    % Oxydizer Fuel Pressure [Pa]

deltaP_ip = 0.3*Pc1;                % Pressure drop injection plate [Pa]
deltaP_fl = 20000;                  % Pressure drop feed lines [Pa]
deltaP_dyn = 0.5*dens_prop1*10;     % Dynamic pressure loss after valve [Pa]

%% RAO Approximation

% Data
alpha = 15;
de1 = 2*sqrt(Ae1)/sqrt(pi);
dt1 = 2*sqrt(At1)/sqrt(pi);

% Minimum length approach (mla)
perc_mla = 0.75;
L_dc = 0.5*(de1-dt1)/tand(alpha);
L_rao = perc_mla * L_dc;
theta_i = 32;
theta_e = 10;
alpha_rao = atand((de1-dt1)/(2*L_rao));
lambda = 0.5 * (1 + cosd ((alpha_rao+theta_e)/2));

thrust_rao = lambda*Thrust1;

% Max performance approach (mpa)
perc_mpa = 1;
L_dc = 0.5*(de1-dt1)/tand(alpha);
L_rao = perc_mpa * L_dc;
theta_i = 27;
theta_e = 5;
alpha_rao = atand((de1-dt1)/(2*L_rao));
lambda = 0.5 * (1 + cosd ((alpha+theta_e)/2));

thrust_rao = lambda*Thrust1;









