function [ returnValue,ceq  ] = switchFunCostr( requestedCase,IN,Parameters,Skf,Skk,Skh,fvals_f,fvals_c,ISsol_f,ISsol_c,Hf,Hc)

ceq = [];

persistent fitnessValue ISsol oldInputs;

if isempty(oldInputs)
    oldInputs = []; %initialize oldInputs as cell
end

if size(oldInputs,1) == 0
    isAllreadyCalculated = 0;
else
[isAllreadyCalculated, iOldInput] =...
ismember(IN, oldInputs, 'rows'); 
end

if isempty(oldInputs) || ~isAllreadyCalculated
     
    [ftemp,istemp] = soo(IN,Parameters,Skf,Skk,Skh,fvals_f,fvals_c,ISsol_f,ISsol_c,Hf,Hc,Parameters.MM.forceCl);
    

   fitnessValue(end+1)=ftemp;
   ISsol(:,end+1) =istemp;

   
  
    oldInputs = [oldInputs;IN]; %add current inputs to persistent array
    returnValue_cell = {fitnessValue(end), ISsol(:,end)};

else
    returnValue_cell = {fitnessValue(iOldInput), ISsol(:,iOldInput)}; % return old values
end

if strcmp(requestedCase, 'FitnessValue')

    returnValue = returnValue_cell{1};
    
elseif strcmp(requestedCase, 'ISsol')
    
    if Parameters.MM.forceCl == 0    
      c = returnValue_cell{2};
      returnValue = c;
    
    elseif  Parameters.MM.forceCl == 1
      c = returnValue_cell{2};
      returnValue = c(1);
      %ceq = c(2);    %  <- GA does not solve problems with integer and equality constraints.
                      % For more help see No Equality Constraints in the documentation.

    end
      
else
    error('MyApp:SelectOutput','Requested case not supported')
end

end %function

