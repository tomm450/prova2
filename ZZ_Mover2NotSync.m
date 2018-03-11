dt = datestr(now,'yyyy_mm_dd_HH_MM_SS');

if exist('../NOT_SYNC','dir') ~= 1 && exist('./PathAdder.m','file') > 0
    system('mkdir ../NOT_SYNC');
end

system(sprintf('mkdir ../NOT_SYNC/%s',dt));

% SPARSE FILE
system(sprintf('mv ./*.txt ../NOT_SYNC/%s/',dt));
system(sprintf('mv ./*.csv ../NOT_SYNC/%s/',dt));

% OUTPUT
system(sprintf('mkdir ../NOT_SYNC/%s/Output',dt));
system(sprintf('mv ./Output/* ../NOT_SYNC/%s/Output/',dt));
segnaposto('Output','Output');

% CASES_FOLDER
system(sprintf('mkdir ../NOT_SYNC/%s/Cases_folder',dt));
system(sprintf('mv ./Cases_folder/* ../NOT_SYNC/%s/Cases_folder/',dt));
segnaposto('Cases_folder','Cases_folder');

% CALLNUMBER
fid = fopen('./callNumber.txt','w+');
fprintf(fid,'1');
fclose(fid);

function segnaposto(dir,string)

fid = fopen(sprintf('./%s/readme',dir),'w+');
fprintf(fid,'%s',string);
fclose(fid);

end

