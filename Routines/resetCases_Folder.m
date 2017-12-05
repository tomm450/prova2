if exist('./callNumber.txt') % sono nella cartella giusta
   
    [a,b] = system('cd ./Cases_folder && ls -1');
    
    b = strsplit(b,'\n');
    jj = 1;
    fprintf('Da Cancellare: \n');
    for j = 1:size(b,2)-1
        if strcmp(b{j},'readme')
            % do nothing
        else
            fprintf('%s \n',sprintf('rm -r ./Cases_folder/%s',b{j}));
            cmd_c{jj} = sprintf('rm -r ./Cases_folder/%s',b{j});
            jj = jj + 1;
            %system(sprintf('rm -r ./Cases_folder/%s',b{j}));
        end
    end
    
    reply = input('Do you want proceed? Y/N [Y]:','s');
    if isempty(reply)
        reply = 'Y';
    end
    
    
    if strcmpi(reply,'y')
        
        % cancello
        fprintf('\nCancello ...\n');
        for k = 1:jj-1
            fprintf('%s \n',cmd_c{k});
            system(sprintf('%s',cmd_c{k}));
        end
        % reinizializzo contatore
         
        fid = fopen('./callNumber.txt','w+');
        fprintf(fid,'1');
        fclose(fid);
        
    end
else
    error('Eseguire solo dalla cartella principale!')
end
