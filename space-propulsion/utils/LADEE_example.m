clc; close all; clear; 


%% LADEE

%N2O4+MMH

mDry=258; %[kg]
mTot=383; %[kg]
Pc=9.4; %[bar]

T=445; %[N]
Is=322; %[s]
g0=9.81; %[m/s^2]
%mProp=;
OF=1.65; %mass/mass
expRatio=375;
deltaV=Is*g0*log(mTot/mDry) %[m/s]

%% PROJECT

%N2O4+MMH

mDry=250; %[kg]
deltaV=2500; %[m/s]

Pc=20.4; %[bar]
T=2452; %[N]
Is=315.5; %[s]
OF=2.0; %mass/mass
g0=9.81; %[m/s^2]

mProp=mDry*exp(deltaV/(Is*g0))-mDry

mOx = (mProp*OF)/(1+OF)
mFu = mProp - mOx

c=Is*g0;
dm=T/c

tBurn=mProp/dm

dmOx = (dm*OF)/(1+OF)
dmFu = dm - dmOx

Pinj=1.3*Pc;

%At=





