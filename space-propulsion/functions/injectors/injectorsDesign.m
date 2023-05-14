function [Ainj, nInj] = injectorsDesign(AfuTot, AoxTot, engine)
%
%   Ainj[Fu, Ox]
%   nInj[Fu, Ox]

    if isnan(engine.injOFmult)
        % Like on Like
        nOx = floor(AoxTot/engine.AinjMin);
        nFu = floor(AfuTot/engine.AinjMin); 
        
        % obtain even numbers
        nOx = nOx - mod(nOx, 2); 
        nFu = nFu - mod(nFu, 2); 
        
        if (nOx == 0) nOx = 2; end 
        if (nFu == 0) nFu = 2; end

    elseif AfuTot > AoxTot/engine.injOFmult
        % nOx > nFu
        nOx = floor(AoxTot/engine.AinjMin); 
        nFu = floor(nOx/engine.injOFmult);
        nOx = engine.injOFmult*nFu; 
    else
        %nFu > nOx
        nFu = floor(AfuTot/engine.AinjMin); 
        nOx = engine.injOFmult*nFu; 
    end
    
    if nFu == 0
        nFu = 1;
    end
    

    if not(isnan(engine.nInjMax)) && (nFu + nOx > engine.nInjMax)
        nFu = floor(engine.nInjMax/(1 + engine.injOFmult));
        nOx = engine.injOFmult*nFu;
    end

    AFu = AfuTot/nFu; 
    AOx = AoxTot/nOx; 
    
    Ainj = [AFu, AOx]; 
    nInj = [nFu, nOx]; 

end

