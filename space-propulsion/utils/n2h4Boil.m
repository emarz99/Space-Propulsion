clc; close all; clear;

T=linspace(280,300,200);

pVapHydrazine9=@(t) exp(98.9831 - 8384.67./t - 13.269*log(t) + 0.116005*t) ./ 14.0537;

pVapHydrazine8=@(t) exp(59.5786 - 7015.15./t - 6.77418*log(t) + 0.00418762*t) ./ 14.0537;

pVapHydrazine1=@(t) exp(58.7582 - 0.707e4./t - 7.088*log(t) + 0.00457*t);

% plot(T,pVapHydrazine9(T));%,T,pVapHydrazine8(T))
% hold on
% grid on
% plot(300,pVapHydrazine9(300),"Marker",".","MarkerSize",10)
% %plot(450,pVapHydrazine8(450),"Marker","o","MarkerSize",10)
% plot(T, 13*ones(size(T)))

pT = 13;

fun9 = @(x) pVapHydrazine9(x) - pT; 
fun8 = @(x) pVapHydrazine8(x) - pT; 
fun1 = @(x) pVapHydrazine1(x) - pT; 

[x1, f1] = fzero(fun1, 1000)
[x8, f8] = fzero(fun8, 1000)
[x9, f9] = fzero(fun9, 1000)

