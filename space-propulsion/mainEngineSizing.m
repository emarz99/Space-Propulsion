clc; 
clear; 
close all; 
configEngineSizing; 

addpath(engine.ceaPath); 
addpath(genpath("functions\")); 

normalize = @(x) x./norm(x);

wNorm = normalize([engine.wTemp, engine.wThrust]);
engine.wTemp = wNorm(1); 
engine.wThrust = wNorm(2); 


fun = @(x) engineOpt(x, engine); 

x0 = [engine.OF]'; 

fun(x0); 
%[~, par] = fun(x0); 

% x0(1) = par.TFinj;

A = [-1 ; 1 ]; 
%       0 -1 
%       0 1 ];

b = [-lb; ub]; b = b(1:end)'; 

options = optimoptions('ga', 'MaxStallGenerations', 5, 'FunctionTolerance', ...
        1e-4, 'MaxGenerations', 10, 'NonlinearConstraintAlgorithm', 'penalty',...
        'PopulationSize', 200, 'Display', 'iter');


[x, fval] = ga(fun,1,A, b,[],[],lb,ub,[],[], options);

engineFinal = engineDesign(x, engine); 

%MCanalysis; 

