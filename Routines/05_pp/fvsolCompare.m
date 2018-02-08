% case to copy
CTP = [500:600];
% cartella principale output
cf = './Cases_folder';
% new folder
cn = './outfSol';



%system(sprintf('mkdir %s',cn));

for i = 1:size(CTP,2)
   I = CTP(i);
   
   if exist(sprintf('%s/%d',cf,I),'dir')
   
   % copio e rinomino
   fprintf('cp -r %s/%d/30simple/system/fvSolution %s/%d \n',cf,I,cn,I);
 
   
   system(sprintf('cp -r %s/%d/30simple/system/fvSolution %s/',cf,I,cn));
   system(sprintf('mv %s/fvSolution %s/%d',cn,cn,I));
   
   end
   
end





