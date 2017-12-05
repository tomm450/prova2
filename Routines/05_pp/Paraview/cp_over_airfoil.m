function [punti_ventre_air, punti_dorso_air, punti_ventre_slat, punti_dorso_slat ] ...
    = cp_over_airfoil(case2read,MAIN_folder)

case_folder = sprintf('%s/Cases_folder/%d',MAIN_folder,case2read);

system(sprintf('touch %s/30simple/%d.foam',case_folder,case2read));

TO_BE_LOAD = sprintf('%s/30simple/%d.foam',case_folder,case2read);
TO_BE_OPEN = sprintf('%s/PPy.py',case_folder);
TO_BE_SAVE = sprintf('%s/%d.csv',case_folder,case2read);

[done] = PPywriter(TO_BE_OPEN,TO_BE_LOAD,TO_BE_SAVE);

goGoOpenFOAM(sprintf('pvbatch %s',TO_BE_OPEN));

system(sprintf('cd %s && tail -n +2 %d0.csv > %dc.csv',case_folder,case2read,case2read));

air_mat = csvread(sprintf('%s/%dc.csv',case_folder,case2read));

% devo riordinare e reinterpolare Cp(x) sui punti noti
cd = load(['./Output/',num2str(case2read),'current_design.mat']);

LE_air = [cd.xp(1,0.5*(size(cd.xp,2)+1)); cd.yp(1,0.5*(size(cd.xp,2)+1))];
LE_sla = [cd.xp(2,0.5*(size(cd.xp,2)+1)); cd.yp(2,0.5*(size(cd.xp,2)+1))];

TE_air = [cd.xp(1,1); cd.yp(1,1)];
TE_sla = [cd.xp(2,1); cd.yp(2,1)];

XZ = [air_mat(:,8),air_mat(:,10)];

TGT = [LE_air,TE_air,LE_sla,TE_sla];

for k = 1:size(TGT,2)
    err2 = (XZ(:,1) - TGT(1,k)).^2 + (XZ(:,2) - TGT(2,k)).^2;
    
    [err2,imin] = min(abs(err2));
    
    tgt(k) = imin;
    XZ = [XZ(1:imin-1,:);XZ(imin+1:end,:)];
    air_mat = [air_mat(1:imin-1,:);air_mat(imin+1:end,:)];
end

punti_dorso_air = [];
punti_ventre_air = [];
punti_dorso_slat = [];
punti_ventre_slat = [];

pps = 0.5*(size(cd.xp,2)+1);
ventre_air = [cd.xp(1,2:pps-1);
              cd.yp(1,2:pps-1)]; 
          
dorso_air  = [cd.xp(1,pps+1:end);
              cd.yp(1,pps+1:end)];

ventre_slat = [cd.xp(2,2:pps-1);
               cd.yp(2,2:pps-1)]; 
          
dorso_slat  = [cd.xp(2,pps+1:end);
               cd.yp(2,pps+1:end)];
          
         
for i = 1:size(air_mat,1)
      
    [yp1] = spline(ventre_air(1,:) ,ventre_air(2,:) ,air_mat(i,8));
    [yp2] = spline(dorso_air(1,:)  ,dorso_air(2,:)  ,air_mat(i,8));
    [yp3] = spline(ventre_slat(1,:),ventre_slat(2,:),air_mat(i,8));
    [yp4] = spline(dorso_slat(1,:) ,dorso_slat(2,:) ,air_mat(i,8));
       
    [err,imin] = min(abs([yp1;yp2;yp3;yp4] - ones(4,1)*air_mat(i,10)));
        
    if imin == 1
        punti_ventre_air = [punti_ventre_air; air_mat(i,:)];
    elseif imin == 2
        punti_dorso_air = [punti_dorso_air; air_mat(i,:)];
    elseif imin == 3
        punti_ventre_slat = [punti_ventre_slat; air_mat(i,:)];
    elseif imin == 4
        punti_dorso_slat = [punti_dorso_slat; air_mat(i,:)];
    end
    
end

[~,idx] = sort(punti_ventre_air(:,8)); % sort just the first column
punti_ventre_air = punti_ventre_air(idx,:);   % sort the whole matrix using the sort indices

[~,idx] = sort(punti_dorso_air(:,8)); % sort just the first column
punti_dorso_air = punti_dorso_air(idx,:);   % sort the whole matrix using the sort indices

[~,idx] = sort(punti_ventre_slat(:,8)); % sort just the first column
punti_ventre_slat = punti_ventre_slat(idx,:);   % sort the whole matrix using the sort indices

[~,idx] = sort(punti_dorso_slat(:,8)); % sort just the first column
punti_dorso_slat = punti_dorso_slat(idx,:);   % sort the whole matrix using the sort indic







