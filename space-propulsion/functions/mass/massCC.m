function [massCChamber] = massCC(rhoMetalCC, lenght, lenghtConv, radiusThroat, radius, thickness)
% function doesn't account for cooling channels
% 3 equal combustion chamber 

volCilinder = emptyCilinderVol(3, radius, lenght, thickness) ;
volConv = emptyConeVol(radius, radiusThroat, lenghtConv, thickness) ; 

volume = volCilinder + volConv ;

massCChamber = mass(rhoMetalCC, volume) ; 

end