clean

% case to restart
Ctr = 9;
IterNew = 10;


for i = 1:size(Ctr,2)
   
   % controllo se esiste directory 
   if exist(strcat('./Cases_folder/',num2str(Ctr(i))),'dir') == 0
       error('Caso %d non è stato calcolato in primo luogo \n',Ctr(i));
   end
   
   % quale solutore? 
   if     exist(strcat('./Cases_folder/',num2str(Ctr(i)),'/30simple'),'dir') == 7
       % solutore simple
       sol = '30simple';
       
   elseif exist(strcat('./Cases_folder/',num2str(Ctr(i)),'/40piso'),'dir') == 7
       % solutore piso
       sol = '40piso';
   
   else
       % qualcosa non quadra
       error('Caso %d non è presente cartella soluzione \n',Ctr(i)); 
       
   end
   
   % calcolo iniziato?
   if exist(sprintf('./Cases_folder/%d/%s/postProcessing',Ctr(i),sol),'dir') == 0
       error('Caso %d non è presente cartella postProcessing \n',Ctr(i));
   end
   
   % iterazione di arrivo?
   [winLs,Ls] = system(sprintf('cd ./Cases_folder/%d/%s/postProcessing/forceCoeffs/ && ls -1',...
       Ctr(i),sol));
   
   if winLs ~= 0
       error('Caso %d: ls forceCoeffs non a buon fine \n');
   else
       clear winLs
   end
   
   Ls = strsplit(Ls,'\n');
   ripTime = nan(1,max(size(Ls))-1);
   
   for j = 1:max(size(Ls))-1
       ripTime(j) = str2double(Ls{j});
   end
   
   [~,imaxJ] = max(ripTime);
  
   [winN,N] = system(sprintf('tail -n 1 ./Cases_folder/%d/%s/postProcessing/forceCoeffs/%s/forceCoeffs.dat',...
       Ctr(i),sol,Ls{imaxJ}));
   
   if winN ~= 0
       error('Caso %d: tail forceCoeffs non a buon fine \n');
   else
       clear winN
   end
   
   N = strsplit(N,'\t'); N = str2double(strtrim(N{1}));
   
   clear j Ls ripTime imaxJ 
   
   % numero processori?
   [winLs,Ls] = system(sprintf('cd ./Cases_folder/%d/%s && ls -1',Ctr(i),sol));
   
   if winLs ~= 0
       error('Caso %d: ls cartelle non a buon fine \n');
   else
       clear winLs
   end
   
   Ls = strsplit(Ls,'\n');
   n_proc = 0;
   
   for j = 1:max(size(Ls))
       if strncmp('processor',Ls{j},9)
           n_proc = n_proc +1 ;
       end
   end
   clear j Ls
      
   
   % quante ripartenze?
   [winLs,Ls] = system(sprintf('cd ./Cases_folder/%d/%s/system && ls -1',Ctr(i),sol));
   
   if winLs ~= 0
       error('Caso %d: ls dizionari non a buon fine \n');
   else
       clear winLs
   end
   
   Ls = strsplit(Ls,'\n');
   n_rip = 0;
   
   for j = 1:max(size(Ls))
       if strncmp('controlDict',Ls{j},11)
           n_rip = n_rip +1 ;
       end
   end
   clear j Ls
   
   % modifica controldict
   winCp = system(sprintf('cp -f ./Cases_folder/%d/%s/system/controlDict ./Cases_folder/%d/%s/system/controlDict%d',Ctr(i),sol,Ctr(i),sol,n_rip-1));
   [winIntro,Intro] = system(sprintf('head -n 18  ./Cases_folder/%d/%s/system/controlDict',Ctr(i),sol));
   [winOutro,Outro] = system(sprintf('tail -n +22 ./Cases_folder/%d/%s/system/controlDict',Ctr(i),sol));
   
   if (winCp+1)*(winIntro+1)*(winOutro+1) ~= 1
       error('Caso %d: tail e head di controlDict non a buon fine \n');
   end  
   
   clear winCp winIntro winOutro
   
   fid = fopen(sprintf('./Cases_folder/%d/%s/system/controlDict',Ctr(i),sol),'w+');
   
   fprintf(fid,'%s \n',Intro);
 
   fprintf(fid,'startTime       latestTime; \n'); 
   fprintf(fid,'stopAt          endTime;  \n'); 
   fprintf(fid,'endTime         %d; \n\n',IterNew+N);

   fprintf(fid,'%s \n',Outro);
   
   fclose(fid);
   
   %% rilancio
   tstart = tic;
   %goGoOpenFOAM(sprintf('cd ./Cases_folder/%d/%d && mpirun -n %d %sFoam -parallel >> ../3logFoam.txt',...
   fprintf('cd ./Cases_folder/%d/%s && mpirun -n %d %sFoam -parallel > ../3logFoam%d.txt',...
       Ctr(i),sol,n_proc,sol(3:end),n_rip);
   
   %%
   case_dir = sprintf('./Cases_folder/%d',Ctr(i));
   
   % checkConvergenza
   [dummy,solverReport] = system(sprintf('tail -n 19 %s/3logFoam%d.txt',case_dir,n_rip));
   solverReport
   
   crono(3) = toc(tstart)/60;
   
   fprintf('Foam in %f min \n',crono(2));
   crono
   
   % estraggo cl
   [good_cmd,reading] = system(sprintf('cat %s/3logFoam%d.txt | grep "Cl   " | tail -n 1',case_dir,n_rip));
   reading = strsplit(reading,' ');
   f = -str2num(reading{end});
   
   [good_cmd,reading] = system(sprintf('cat %s/3logFoam%d.txt | grep "Cl    = " | cut -d '' '' -f10 | tail -200',case_dir,n_rip));
   
   f_conv = -str2num(reading);
   
   f_conv = std(f_conv);
   fprintf('std(Cl(end-200:end) = %f \n',f_conv);
   cp = f_conv;
   save(sprintf('%s/ID_calcolo.mat',case_dir),'X_IN','IN','RES_struct','Parameters','f','f_conv','RES_struct');
   
   %% TAGLIA-CUCI
   %[punti_ventre_air, punti_dorso_air, punti_ventre_slat, punti_dorso_slat ] ...
   %    = cp_over_airfoil(5,my_dir);
   
   %%
   case_folder = sprintf('%s/Cases_folder/%d',pwd,Ctr(i));
   
   if strcmp(SOLVER.solver,'simple')
       system(sprintf('touch %s/30simple/%d.foam',case_folder,Ctr(i)));
       TO_BE_LOAD = sprintf('%s/30simple/%d.foam',case_folder,Ctr(i));
   else
       system(sprintf('touch %s/40piso/%d.foam',case_folder,Ctr(i)));
       TO_BE_LOAD = sprintf('%s/40piso/%d.foam',case_folder,Ctr(i));
   end
   TO_BE_OPEN = sprintf('%s/PPy.py',case_folder);
   TO_BE_SAVE = sprintf('%s/%d.csv',case_folder,Ctr(i));
   
   % scrivo script pyton
   [done] = PPywriter(TO_BE_OPEN,TO_BE_LOAD,TO_BE_SAVE);
   % eseguo
   goGoOpenFOAM(sprintf('pvbatch %s',TO_BE_OPEN));
   % modifico output in modo da poterlo direttamente importare
   system(sprintf('cd %s && tail -n +2 %d0.csv > %dc.csv',case_folder,Ctr(i),Ctr(i)));
   
   air_mat = csvread(sprintf('%s/%dc.csv',case_folder,Ctr(i)));
   
   % %% per calibrazione commentare
   % cd = load('current_design.mat');
   %
   % LE_air = [cd.xp(1,0.5*(size(cd.xp,2)+1)); cd.yp(1,0.5*(size(cd.xp,2)+1))];
   % LE_sla = [cd.xp(2,0.5*(size(cd.xp,2)+1)); cd.yp(2,0.5*(size(cd.xp,2)+1))];
   %
   % TE_air = [cd.xp(1,1); cd.yp(1,1)];
   % TE_sla = [cd.xp(2,1); cd.yp(2,1)];
   %
   % XZ = [air_mat(:,8),air_mat(:,10)];
   %
   % TGT = [LE_air,TE_air,LE_sla,TE_sla];
   %
   % for k = 1:size(TGT,2)
   %     err2 = (XZ(:,1) - TGT(1,k)).^2 + (XZ(:,2) - TGT(2,k)).^2;
   %
   %     [err2,imin] = min(abs(err2));
   %
   %     tgt(k) = imin;
   %     XZ = [XZ(1:imin-1,:);XZ(imin+1:end,:)];
   %     air_mat = [air_mat(1:imin-1,:);air_mat(imin+1:end,:)];
   % end
   %
   % punti_dorso_air = [];
   % punti_ventre_air = [];
   % punti_dorso_slat = [];
   % punti_ventre_slat = [];
   %
   % pps = 0.5*(size(cd.xp,2)+1);
   % ventre_air = [cd.xp(1,2:pps-1);
   %               cd.yp(1,2:pps-1)];
   %
   % dorso_air  = [cd.xp(1,pps+1:end);
   %               cd.yp(1,pps+1:end)];
   %
   % ventre_slat = [cd.xp(2,2:pps-1);
   %                cd.yp(2,2:pps-1)];
   %
   % dorso_slat  = [cd.xp(2,pps+1:end);
   %                cd.yp(2,pps+1:end)];
   %
   % for i = 1:size(air_mat,1)
   %
   %     [yp1] = spline(ventre_air(1,:) ,ventre_air(2,:) ,air_mat(i,8));
   %     [yp2] = spline(dorso_air(1,:)  ,dorso_air(2,:)  ,air_mat(i,8));
   %     [yp3] = spline(ventre_slat(1,:),ventre_slat(2,:),air_mat(i,8));
   %     [yp4] = spline(dorso_slat(1,:) ,dorso_slat(2,:) ,air_mat(i,8));
   %
   %     [err,imin] = min(abs([yp1;yp2;yp3;yp4] - ones(4,1)*air_mat(i,10)));
   %
   %     if imin == 1
   %         punti_ventre_air = [punti_ventre_air; air_mat(i,:)];
   %     elseif imin == 2
   %         punti_dorso_air = [punti_dorso_air; air_mat(i,:)];
   %     elseif imin == 3
   %         punti_ventre_slat = [punti_ventre_slat; air_mat(i,:)];
   %     elseif imin == 4
   %         punti_dorso_slat = [punti_dorso_slat; air_mat(i,:)];
   %     end
   %
   % end
   %
   % [~,idx] = sort(punti_ventre_air(:,8)); % sort just the first column
   % punti_ventre_air = punti_ventre_air(idx,:);   % sort the whole matrix using the sort indices
   %
   % [~,idx] = sort(punti_dorso_air(:,8)); % sort just the first column
   % punti_dorso_air = punti_dorso_air(idx,:);   % sort the whole matrix using the sort indices
   %
   % [~,idx] = sort(punti_ventre_slat(:,8)); % sort just the first column
   % punti_ventre_slat = punti_ventre_slat(idx,:);   % sort the whole matrix using the sort indices
   %
   % [~,idx] = sort(punti_dorso_slat(:,8)); % sort just the first column
   % punti_dorso_slat = punti_dorso_slat(idx,:);   % sort the whole matrix using the sort indic
   %
   % % reinterpolo sui punti noti da metodo course
   % x_ventre_air = fliplr(xc(1,1:(size(xc,2)/2)));
   % x_dorso_air  = xc(1,(size(xc,2)/2)+1:end);
   %
   % CP_va = spline(punti_ventre_air(:,8),punti_ventre_air(:,7),x_ventre_air);
   % CP_da = spline(punti_dorso_air(:,8),punti_dorso_air(:,7),x_dorso_air);
   %
   % % figure
   % % plot(punti_ventre_air(:,8),punti_ventre_air(:,7)/(0.5*BU_par.Umag^2),'b');
   % % hold on
   % % plot(punti_dorso_air(:,8),punti_dorso_air(:,7)/(0.5*BU_par.Umag^2),'r');
   % % plot(x_ventre_air,CP_va/(0.5*BU_par.Umag^2),'bo');
   % % plot(x_dorso_air,CP_da/(0.5*BU_par.Umag^2),'ro');
   %
   % cp = [fliplr(CP_va),CP_da];
   %
   % x_ventre_slat = fliplr(xc(2,1:(size(xc,2)/2)));
   % x_dorso_slat  = xc(2,(size(xc,2)/2)+1:end);
   %
   % CP_vs = spline(punti_ventre_slat(:,8),punti_ventre_slat(:,7),x_ventre_slat);
   % CP_ds = spline(punti_dorso_slat(:,8),punti_dorso_slat(:,7),x_dorso_slat);
   %
   % % plot(punti_ventre_slat(:,8),punti_ventre_slat(:,7)/(0.5*BU_par.Umag^2),'b');
   % % hold on
   % % plot(punti_dorso_slat(:,8),punti_dorso_slat(:,7)/(0.5*BU_par.Umag^2),'r');
   % % plot(x_ventre_slat,CP_vs/(0.5*BU_par.Umag^2),'bo');
   % % plot(x_dorso_slat,CP_ds/(0.5*BU_par.Umag^2),'bo');
   %
   %
   % cp = [cp,fliplr(CP_vs),CP_ds]/(0.5*BU_par.Umag^2);
   % cp = cp';
   cp = 'dummy';
   
 
 
 
   
end