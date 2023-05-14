function [massTot] = massFeedLine(rhoMetalFeedingLine, radiusTank, lenghtTank, thicknessTank, radiusPipe, lenghtPipe, thicknessPipe)

volTank = emptyCilinderVol(1, radiusTank, lenghtTank, thicknessTank) ;
massTank = mass(rhoMetalFeedingLine, volTank) ;

volPipe = emptyCilinderVol(2, radiusPipe, lenghtPipe, thicknessPipe) ;
massPipe = mass(rhoMetalFeedingLine, volPipe) ;

massTot = massTank + massPipe ;

end