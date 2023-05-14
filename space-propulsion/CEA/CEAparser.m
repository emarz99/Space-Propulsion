function [coeff, parsedNames] = CEAparser(input)

%   coeff [nPar, nEl]: parameter number, element numeber 'transposition of
%   the format in the .csv file'

    filename = 'dataCEA.out'; 

    dataCEA = fileread(filename);      % read data    

    %find data EQUILIBRIUM 
    resE = regexp(dataCEA,'THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM'); 
    
    %find data FROZEN
    resF = regexp(dataCEA,'THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION');
    res =[resE, resF]; 
    tok = cell(1,numel(res)); 

    %coefficients to be parsed
    parsedNames = [{'P, BAR'}, {'T, K'}, {'RHO, KG/CU M'}, {'GAMMAs'}, {'SON VEL,M/SEC'}, {'MACH NUMBER'}, {'CONDUCTIVITY  '},...
                    {'PRANDTL NUMBER'}, {'CSTAR, M/SEC'}, {'CF'}, {'Ivac, M/SEC'}, {'Isp, M/SEC'}, {'VISC,MILLIPOISE'}]; 
    nCoeff = length(parsedNames); 
    
    %total stations to be computed
    nCol = 2 + length(input.subar) + length(input.supar); 

    coeff = zeros(nCoeff, nCol); 
    
    %main loop trough the tokens with coefficients, each one contains the full set of coeffs
    for i = 1:length(tok)
        
        if i == length(tok)
            tok{i} = dataCEA(res(i):end);
        else
            tok{i} = dataCEA(res(i):res(i+1)); 
        end
        
        % if it is considering more tokens relatet with the same model (EQUILIBRIUM or FROZEN), consider
        % an offet in the position in wich saving coefficients (see later)
        if i <= length(resE)
            offset = 8*(i-1);
        elseif i <= length(resE) + length(resF)
            offset = 8*(i-1-length(resE));
        end

        % search for the line containing each coefficient's data
        for k = 1:nCoeff
            
            coeffTemp = nan(1, nCol);

            posI = regexp(tok{i}, parsedNames{k});  %beginning of the line
            posI = posI(1); 
            flag = true;
            
            % find the end of the line and save it in posF
            posF = posI; 
            while flag
                posF = posF + 1; 
                if isequal(tok{i}(posF), newline)
                    posF = posF - 1; 
                    flag = false;
                end
            end

            % string containing the data for each coefficient in the given token
            strCoeff = tok{i}(posI:posF); 
            
            nChar = length(strCoeff);           %compute the total number of char in the string with the data of one coeff
            posSpace = regexp(strCoeff, ' ');   %compute the positions of spaces in the same string
            posCoeff = ones(1, nChar);                          
            posCoeff(posSpace) = 0;             %compute the positions with data as the ones with no spaces
            
            % FAKE represent the number of fake trigger in the function that searches for numbers, it works by finding the
            % consecutive elements that are no spaces, but in the string there is also the coefficient name which is not a number    
            fake = length(regexp(parsedNames{k}, ' ')) + 1;  
            
            % this is the case of CONDUCTIVITY that in the names has two extra spaces in order to not being confuded with 
            % other descriptions in the .out file
            if k == 7
                fake = 1; % trust me bro
            end

            % In each token data regarding CHAMBER and THROAT are repeated so if it is not the first token
            % you have to skip them
            if offset > 0
                fake = fake + 2;
            end

            %%% Read the coefficients
            cont = 1;
            nData = 0; 

            % keep scrolling till the end of the string
            while cont <= nChar
               
                %flag is true if it is reading a coefficient

                if not(flag) && posCoeff(cont) == 1
                    flag = true;  %if I was not reading anything and I find a 1, it means it is a meaningfull char and I start reading
                    posI = cont;  % save initial position
                    posF = cont;  % final position, waiting for update
                end

                notInnerWhile = true; % usefull to avoid double increasing in counter
                
                while flag && cont <=nChar                  
                    
                    if posCoeff(cont) == 1 && cont ~=nChar
                        % keep final position updated wirh the last reading
                        posF = cont; 
                    else
                        % if while it is reading it reaches the end or finds a space, stops reading
                        if cont == nChar
                            posF = cont; 
                        end
                        flag = false;           % Stop reading
                        nData = nData + 1;      % Number of data found from the beginning
                        
                        % Control that the number of elements found is bigger than the "fakes" ones
                        if nData>fake           

                            % proceed to save the string considering the exponential notation
                            strVal = strCoeff(posI:posF); 
                            
                            posExp = regexp(strVal, '-');
                            if isempty(posExp)
                                posExp = regexp(strVal, '+'); 
                            end
                            
                            if not(isempty(posExp))
                                strVal = strcat(strVal(1:posExp-1), 'e',strVal(posExp:end));
                            end
                            
                            % apply shif if needed
                            col = (nData-fake)+offset;

                            % Cstar, CT, Isp, Ivac do not have data in the COMBUSTION CHAMBER, so you have to skip the first column
                            if (k>=9 && k <=12)
                                col = col + 1;  % trust me
                            end

                            if i <= length(resE) 
                                % EQUILIBRIUM
                                coeffTemp(1,col) = sscanf(strVal,'%f');
                            else 
                                % FROZEN
                                if col > 2
                                    col = col + length(input.subar); 
                                end
                                coeffTemp(1,col) = sscanf(strVal,'%f');
                            end
                        end
                    end
                    cont = cont + 1; 
                    notInnerWhile = false; 
                end
                
                if notInnerWhile
                    cont = cont + 1; 
                end
            end
            
            % remove zeros if they are exponential part of the number
            if length(coeffTemp) > nCol
                coeffTemp = coeffTemp(coeffTemp~=0); 
            end
            
            index = not(isnan(coeffTemp)); 

            coeff(k, index) = coeffTemp(index); 
        end
    end
    

end