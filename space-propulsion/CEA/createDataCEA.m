function coeff = createDataCEA(input, ceaPath)
    
    currentPath = pwd; 
    
    cd(ceaPath)

    fIN = fopen('dataCEA.inp', 'w'); 
    
    %%% finite area combustor
    fprintf(fIN, 'problem rocket o/f=%d ', input.OF);
    if input.problem.facFlag 
        fprintf(fIN, 'fac acat=%f ', input.problem.acat); 
    end
    
    %%% equlibrium / frozen
    if isequal(input.problem.solution, 'frozen')
        fprintf(fIN, 'frozen nfz=%d\r\n',input.problem.nfz); 
    end
    
    if isequal(input.problem.solution, 'equilibrium') || (input.problem.nfz > 1)
        fprintf(fIN, 'equilibrium\r\n'); 
    end

    %%% pressure assigned
    strTemp = 'p,bar='; 
    for i = input.pc
        strTemp = strcat(strTemp,num2str(i),','); 
    end
    strTemp = strcat(strTemp, '\r\n'); 
    fprintf(fIN, strTemp); 

    %%% subar
    if not(isempty(input.subar))
        strTemp = 'subar='; 
        for i = input.subar
            strTemp = strcat(strTemp,num2str(i),','); 
        end
        strTemp = strcat(strTemp,'\r\n'); 
        fprintf(fIN, strTemp); 
    end

    %%% supar
    if not(isempty(input.supar))
        strTemp = 'supar='; 
        for i = input.supar
            strTemp = strcat(strTemp,num2str(i),','); 
        end
        strTemp = strcat(strTemp, '\r\n'); 
        fprintf(fIN, strTemp);
    end

    %%% reactants
    fprintf(fIN,'\treactants\n'); 
    for i = 1:length(input.fuel.name)
        if not(isnan(input.fuel.wt(i)))
            fprintf(fIN, strcat('fuel = ',input.fuel.name{i},' wt%%=%f t(k)=%f\r\n'), input.fuel.wt(i), input.fuel.T(i)); 
        else
            fprintf(fIN, strcat('fuel = ',input.fuel.name{i},' t(k)=%f\r\n'), input.fuel.T(i)); 
        end
    end
    
    for i = 1:length(input.oxidizer.name)
        if not(isnan(input.oxidizer.wt(i)))
            fprintf(fIN, strcat('oxid = ',input.oxidizer.name{i},' wt%%=%f t(k)=%f\r\n'), input.oxidizer.wt(i), input.oxidizer.T(i)); 
        else
            fprintf(fIN, strcat('oxid = ',input.oxidizer.name{i},' t(k)=%f\r\n'), input.oxidizer.T(i)); 
        end        
    end

    %%% plot selected data 
    strTemp = strcat('\toutput siunits transport\n\r\tplot',{' '},input.coeffOut,'\r\n\tend');
    
    fprintf(fIN,strTemp{1}); 

    fclose(fIN); 
    
    system('ceaModR.exe'); 
    coeff = CEAparser(input); 
    coeff(1, :) = coeff(1, :)*(1e5);
    coeff(7, :) = coeff(7, :)/10;
    coeff(13, :) = coeff(13, :)*(1e-4);
    
    delete("dataCEA.csv","dataCEA.out","dataCEA.plt", "dataCEA.inp");

    cd(currentPath)
end