if exist('./callNumber.txt') % sono nella cartella giusta
    % reinizializzo contatore
    fid = fopen('./callNumber.txt','w+');
    fprintf(fid,'1');
    fclose(fid);
    
    % pulisco cartelle
    [a,b] = system('cd ./Cases_folder && ls -1');
    
    b = strsplit(b,'\n');
    
    for j = 1:size(b,2)-1
        if strcmp(b{j},'readme')
            % do nothing
        else
            fprintf('rm -r ./Cases_folder/%s \n',b{j})
            system(sprintf('rm -r ./Cases_folder/%s',b{j}));
        end
    end
else
    error('Eseguire solo dalla cartella principale!')
end
