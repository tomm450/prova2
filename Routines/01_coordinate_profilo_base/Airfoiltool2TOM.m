function [ b ] = Airfoiltool2TOM( file_in, flag_save, file_out)
%% ORDINA COORDINATE PROFILO SECONDO STANDARD SCRIPT
% % input 
% - file_in -> nome (compresa estensione) del file che contiene le coordinate 
%        secondo AirfoilTool.com
% - flag_save : se == 1 procedo a salvare il file di coordinate riordinate
% - file_out -> stringa contenente il nome (compresa estensione) del file in 
%               cui andr� a salvare l'uscita 
% % output
% - b -> matrice riordinata secondo standard [ [1 0 1]' , [up LE dwn]' ]

%% Inizio
if nargin == 1
    flag_save = 0;
    file_out = '';
end

s = size(file_in,1);

[~,s] = min(abs(1- file_in(1:floor(0.75*s),1)));

b = [flipud(file_in(1:s,:)); file_in(s+2:end,:)];

if flag_save == 1
    fid = fopen(file_out,'w+');
    
    for j = 1 : size(b,1)
        fprintf(fid,'%12.5f %12.5f \r\n',b(j,1),b(j,2));
    end
    
    fclose(fid);
    
    fprintf('AirfoilTools -> stndard script : Conversione eseguita e file salvato\n');
end
end

