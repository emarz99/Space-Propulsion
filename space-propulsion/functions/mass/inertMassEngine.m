function [massTot] = inertMassEngine(thickness, rhoMetalFeedingLine, radiusTankOx, lenghtTankOx, radiusPipeOx, lenghtPipeOx, radiusTankFu, lenghtTankFu, ...
    radiusPipeFu, lenghtPipeFu, rhoMetalCC, lenghtCC, lenghtConvCC, radiusMaxCC, radiusMinCC, thicknessCC)

data.a = thickness ;



% feed line is tank and pipe
massFeedLineOx = massFeedLine(rhoMetalFeedingLine, radiusTankOx, lenghtTankOx, thickness, radiusPipeOx, lenghtPipeOx, thickness) ; 
massFeedLineFu = massFeedLine(rhoMetalFeedingLine, radiusTankFu, lenghtTankFu, thickness, radiusPipeFu, lenghtPipeFu, thickness) ;

% mass combustion chamber 
massCChamber = massCC(rhoMetalCC, lenghtCC, lenghtConvCC, radiusMaxCC, radiusMinCC, thicknessCC) ;

massTot = massFeedLineOx + massFeedLineFu + massCChamber ;

end