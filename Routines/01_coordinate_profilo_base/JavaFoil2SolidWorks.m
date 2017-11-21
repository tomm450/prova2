function [done] = JavaFoil2SolidWorks(java_out,solid_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% RICORDARSI DI SOSTITUIRE LE VIRGOLE NEL FILE DI JAVAFOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% java_out = stringa contenente nome file salvato da JavaFoil
% solid_in = stringa contenente nome file con cui andrà a disegnare
%            curva in SolidWorks

% INIZIO
z = 0;

M = load(java_out); % 

fid = fopen(solid_in,'w+');

M_mm = M*1000; % Conversione da metri a millimetri

for j = 1 : size(M,1)
    fprintf(fid,'%12.5f %12.5f %12.5f\r\n',M_mm(j,1),M_mm(j,2),z);
end

fclose(fid);

done = fprintf('JavaFoil -> SLDWK : Conversione eseguita \n');
end