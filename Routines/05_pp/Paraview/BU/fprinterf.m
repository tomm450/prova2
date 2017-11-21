n_row = 779;
fid = fopen('prova.m','w+');
for i = 1:n_row

fprintf('tail -n +%d ./p_extr.py | head -n 1 \n',i);
[a,b] = system(sprintf('tail -n +%d ./p_extr.py > ./ptemp.temp',i));
[a,b] = system(sprintf('head -n 1 ./ptemp.temp',i));
b2 = strsplit(b,'\n');
b3 = strrep(b2{1},'''','''''');
fprintf(fid,'fprintf(fid,'' %s \\n'');  \n',b3);



end

fclose all