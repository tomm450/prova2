% case to copy
CTP = [871:881];
% cartella principale output
cf = '../Cases_folder';
% new folder
cn = './backUpLight';


if exist(sprintf('%s',cn),'dir')
else
    system(sprintf('mkdir %s',cn));
end

for i = 1:size(CTP,2)
   I = CTP(i);
   if exist(sprintf('%s/%d',cn,I),'dir')
   else
       
       if exist(sprintf('%s/%d',cf,I),'dir')
       % creo cartella
       fprintf('mkdir %s/%d \n',cn,I)
       [win1,res1] = system(sprintf('mkdir %s/%d',cn,I));
%        % copio mesh
%        fprintf('cp -r %s/%d/20extrude %s/%d/ \n',cf,I,cn,I);
%        [win2,res2] = system(sprintf('cp -r %s/%d/20extrude %s/%d/',cf,I,cn,I));
       % copio log
       fprintf('cp -r %s/%d/*.txt %s/%d/ \n',cf,I,cn,I),
       [win2,res2] = system(sprintf('cp -r %s/%d/*.txt %s/%d/',cf,I,cn,I));
       % copio csv
       fprintf('cp -r %s/%d/*.csv %s/%d/ \n',cf,I,cn,I);
       [win2,res2] = system(sprintf('cp -r %s/%d/*.csv %s/%d/',cf,I,cn,I));
       % copio mat
       fprintf('cp -r %s/%d/*.mat %s/%d/ \n',cf,I,cn,I);
       [win2,res2] = system(sprintf('cp -r %s/%d/*.mat %s/%d/',cf,I,cn,I));
       
       % creo 30 simple
       system(sprintf('mkdir %s/%d/30simple',cn,I));
       
       % copio postProcess
       fprintf('cp -r %s/%d/30simple/postProcessing %s/%d/30simple/',cf,I,cn,I);
       [win2,res2] = system(sprintf('cp -r %s/%d/30simple/postProcessing %s/%d/30simple/',cf,I,cn,I));
       
       
       % copio system
       fprintf('cp -r %s/%d/30simple/system %s/%d/ \n\n',cf,I,cn,I);
       [win2,res2] = system(sprintf('cp -r %s/%d/30simple/system %s/%d/30simple/',cf,I,cn,I));
       end
   end
   
   
end





