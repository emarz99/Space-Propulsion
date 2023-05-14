function [PDropTotalOx, PDropTotalFu, lMaxOx, lMaxFu] = TotalPressureDrop(CoolingOx, CoolingFu,massFlowRateOx, massFlowRateFu, engine)
 
    FlowVelocityOx = radius2velocity(massFlowRateOx, engine.rhoO, engine.pipeDiamO/2) ;
    FlowVelocityFu = radius2velocity(massFlowRateFu, engine.rhoF, engine.pipeDiamF/2) ;
    
    [PDropTotalOx, lMaxOx] = PDrop(CoolingOx, engine.pc, massFlowRateOx, FlowVelocityOx, engine.rhoO, engine.muO) ;
    [PDropTotalFu, lMaxFu] = PDrop(CoolingFu, engine.pc, massFlowRateFu, FlowVelocityFu, engine.rhoF, engine.muF) ;

end


