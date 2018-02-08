function [f_cell] = fCl_core_OPT(X_IN,IN,RES_struct,Parameters,usefulCoor)

STL           = RES_struct.STL;
BU_par        = RES_struct.BU_par;
MESH_par      = RES_struct.MESH_par;
SOLVER        = RES_struct.SOLVER;

if isfield(SOLVER,'PLT_history') == 0
    SOLVER.PLT_history = 0;
end
    
%% funzione cannibalizzata da script...
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n');
% Aggiungo percorso
my_dir = what;
my_dir = my_dir.path;
my_dirs = genpath(my_dir);
addpath(my_dirs);

if exist('callNumber.txt','file')
    ID = load('callNumber.txt');
else
    ID = 1;
end

% creo caso
system(sprintf('cp -r ./CLEAN_CASE ./Cases_folder/%d',ID));

case_dir = sprintf('./Cases_folder/%d',ID);

system(sprintf('mkdir %s/0msh',case_dir));

save(sprintf('%s/ID_calcolo.mat',case_dir),'X_IN','IN','RES_struct','Parameters');

fid = fopen('callNumber.txt','w+');
fprintf('CASO #%d\n',ID);
fprintf(fid,'%d',ID+1);
fclose(fid);

%% calcoli preliminari

Re    = BU_par.Rho*BU_par.L*BU_par.Umag/BU_par.Nu;
Cf    = 0.025/(Re^(1/7));
Tao   = 0.5*Cf*BU_par.Rho*BU_par.Umag^2;
Ufric = sqrt(Tao./BU_par.Rho);

% y+ voluta
if BU_par.wall_function == 1
    
    yplus_tgt = 50;  % uso wall function
    MESH_par.ds1 = (yplus_tgt*BU_par.Nu)/(Ufric*BU_par.Rho);
    
else
    
    yplus_tgt = 1;
    MESH_par.ds1      = (yplus_tgt*BU_par.Nu)/(Ufric*BU_par.Rho);
    BU_par.omega_body = 60*BU_par.Nu/(0.09*(MESH_par.ds1)^2);

end

% calcolo valori inlet
%k_inlet     = 3/2*(BU_par.Umag*0.05)^2; % https://www.cfd-online.com/Forums/openfoam/67728-help-k-epsilon-values-turbulence.html
%k_inlet     = 1e-3*BU_par.Umag^2/Re;
k_inlet = (0.001*BU_par.Nu)*(10* BU_par.Umag/BU_par.L); % Quadrio

%omega_inlet = 5*BU_par.Umag/(BU_par.L);
%omega_inlet = 5*BU_par.Umag/(2*MESH_par.x_dom);
omega_inlet = (10* BU_par.Umag/BU_par.L); % Quadrio


%%
Dmean = MESH_par.l_airfoil;

[~,nlay] = min(abs(Dmean - (MESH_par.ds1*MESH_par.expRatio.^[1:50])));

MESH_par.nlay   = nlay-1;
MESH_par.th_tot = sum(MESH_par.ds1*MESH_par.expRatio.^(nlay-1));

crono = [];

np = Parameters.n_processori;

%% SCRIVO DICT

[done] = decomposeWrite(np,case_dir,SOLVER);

if strcmp(SOLVER.solver,'simple')
     system(sprintf('rm -r ./Cases_folder/%d/40piso',ID));
     [done] = controlWrite('simple',BU_par,SOLVER,case_dir);
     
elseif strcmp(SOLVER.solver,'piso')
    system(sprintf('rm -r ./Cases_folder/%d/30simple',ID));
    disp('DEROGA SU dt, ricorreggere!')
    
     %CFL_exp = 3*U_mag*deltaT/deltaX 
     SOLVERmod.deltaT        = 1*Dmean/(200*BU_par.Umag);
     SOLVERmod.startTime     = 0;
     SOLVERmod.endTime       = SOLVER.endTime*SOLVERmod.deltaT;
      
     SOLVERmod.writeInterval = floor(SOLVER.endTime/100);
     
     fprintf('deltaT = %f \nfinalT = %f \n\n',SOLVERmod.deltaT,SOLVERmod.endTime);
     
     [done] = controlWrite('piso',BU_par,SOLVERmod,case_dir);
     
else
    error('someway,somewhere,someone fuck up');
end

%% MESH

fprintf('\nMESH... \n')
tstart = tic;
if MESH_par.solver == 'gmsh'
    
    %[~] = Gmesher(X_IN,IN,STL,MESH_par,case_dir,Parameters,SOLVER,BU_par.L);
    [~] = Gmesherrd(X_IN,IN,STL,MESH_par,case_dir,Parameters,SOLVER,BU_par.L);
elseif MESH_par.solver == 'snap'
     error('Snapper.m è da aggiornare\n');        
%    [ done,crono ] = Snapper( X_IN,IN,STL,MESH_par,case_dir,Parameters );
end

crono(1) = toc(tstart)/60;

fprintf('\ncd %s/20extrude && checkMesh > ../logcheckMesh.txt \n',case_dir)
[result,~]=goGoOpenFOAM(sprintf('cd %s/20extrude && checkMesh > ../logcheckMesh.txt ',case_dir));
[result,logMesh]=system(sprintf('tail -n 20 %s/logcheckMesh.txt',case_dir));

fprintf('\n%s \n',logMesh);

fprintf('\nMesh in %f min \n\n',crono(1));


%% SCRIVO DICT
tstart = tic;
bu_type = BU_par.BU_type;  % 'freestream';

if strcmp(SOLVER.solver,'simple')
    
    [done] = ConstantDir(case_dir,SOLVER,BU_par);
    
    [done] = BC_Write('simple',BU_par,bu_type,case_dir );    
   
    [done] = BC_turbolent(k_inlet,omega_inlet,case_dir,SOLVER,BU_par);
    
    
    fprintf('decompose\ncd %s/30simple/ && decomposePar -force > ../3logdec.txt\n',case_dir);
    [done] = goGoOpenFOAM(sprintf('cd %s/30simple/ && decomposePar -force > ../3logdec.txt',case_dir));
    fprintf('simple\ncd %s/30simple/ && mpirun -n %d simpleFoam -parallel > ../3logFoam.txt \n',case_dir,Parameters.n_processori)
    
    if SOLVER.PLT_history == 0
        [done] = goGoOpenFOAM(sprintf('cd %s/30simple/ && mpirun -n %d simpleFoam -parallel > ../3logFoam.txt ',case_dir,Parameters.n_processori));
    else
        [done] = goGoOpenFOAM(sprintf('cd %s/30simple/ && gnuplot resPlot & disown',case_dir));
        [done] = goGoOpenFOAM(sprintf('cd %s/30simple/ && mpirun -n %d simpleFoam -parallel > ../3logFoam.txt ',case_dir,Parameters.n_processori));
        [done] = system('pkill gnuplot');
    end
    
    fprintf('reconstructPar \ncd %s/30simple/ && reconstructPar -latestTime > ../3logrec.txt\n',case_dir);
    [done] = goGoOpenFOAM(sprintf('cd %s/30simple/ && reconstructPar -latestTime > yp.txt',case_dir));
    system(sprintf('tail .%s/30simple/yp.txt',case_dir));
    
elseif strcmp(SOLVER.solver,'piso')
    
    [done] = transportProp(case_dir,SOLVER,BU_par); 
    
    [done] = BC_Write('piso',BU_par,bu_type,case_dir );
    
    [done] = BC_turbolent(k_inlet,omega_inlet,case_dir,SOLVER,BU_par);
    
    fprintf('decompose\ncd %s/40piso/ && decomposePar -force > ../3logdec.txt\n',case_dir);
    [done] = goGoOpenFOAM(sprintf('cd %s/40piso/ && decomposePar -force > ../3logdec.txt',case_dir));
    fprintf('piso\ncd %s/40piso/ && mpirun -n %d simpleFoam -parallel > ../3logFoam.txt \n',case_dir,Parameters.n_processori)
    
    if SOLVER.PLT_history == 0
        [done] = goGoOpenFOAM(sprintf('cd %s/40piso/ && mpirun -n %d pisoFoam -parallel > ../3logFoam.txt ',case_dir,Parameters.n_processori));
    else
        [done] = goGoOpenFOAM(sprintf('cd %s/40piso/ && gnuplot resPlot & disown',case_dir));
        [done] = goGoOpenFOAM(sprintf('cd %s/40piso/ && mpirun -n %d pisoFoam -parallel > ../3logFoam.txt ',case_dir,Parameters.n_processori));
        [done] = system('pkill gnuplot');
    end
    
    fprintf('reconstructPar \ncd %s/40piso/ && reconstructPar -latestTime > ../3logrec.txt\n',case_dir);
    [done] = goGoOpenFOAM(sprintf('cd %s/40piso/ && reconstructPar -latestTime > ../3logrec.txt',case_dir));
else
    error('someway,somewhere,someone fuck up');
end

% checkConvergenza
[dummy,solverReport] = system(sprintf('tail -n 19 %s/3logFoam.txt',case_dir));
solverReport

crono(2) = toc(tstart)/60;

fprintf('Foam in %f min \n',crono(2));
crono

% estraggo cl e cd
[good_cmd,reading] = system(sprintf('cat %s/3logFoam.txt | grep "Cl   " | tail -n 1',case_dir));
reading = strsplit(reading,' ');
CL = str2num(reading{end});

[good_cmd,reading] = system(sprintf('cat %s/3logFoam.txt | grep "Cd   " | tail -n 1',case_dir));
reading = strsplit(reading,' ');
CD = str2num(reading{end});

% oscillazioni
[good_cmd,reading] = system(sprintf('cat %s/3logFoam.txt | grep "Cl    = " | cut -d '' '' -f10 | tail -200',case_dir));
CL_conv = str2num(reading);

CL_conv = std(CL_conv);
[good_cmd,reading] = system(sprintf('cat %s/3logFoam.txt | grep "Cd    = " | cut -d '' '' -f10 | tail -200',case_dir));
CD_conv = str2num(reading);

CD_conv = std(CD_conv);
fprintf('std(end-200:end) = [%f %f] \n',CL_conv,CD_conv);


fprintf('READ SIMPLE e SCRIVO CSV \n');
[rcls_output] = readClsimple(ID,[0:4],'./Cases_folder/',1,1,strcat(case_dir,'/CSV'));

f_cell =  {CL,CD,CL_conv,CD_conv};

%% TAGLIA-CUCI
if isnan(usefulCoor) == 1
    return
end

% AGGIORNARE

xc    = usefulCoor{1}/1000;
cd.xp = usefulCoor{2}/1000; % current design
cd.yp = usefulCoor{3}/1000; % current design

case_folder = sprintf('%s/Cases_folder/%d',pwd,ID);

if strcmp(SOLVER.solver,'simple')
    system(sprintf('touch %s/30simple/%d.foam',case_folder,ID));
    TO_BE_LOAD = sprintf('%s/30simple/%d.foam',case_folder,ID);
else
    system(sprintf('touch %s/40piso/%d.foam',case_folder,ID));
    TO_BE_LOAD = sprintf('%s/40piso/%d.foam',case_folder,ID);
end
TO_BE_OPEN = sprintf('%s/PPy.py',case_folder);
TO_BE_SAVE = sprintf('%s/%d.csv',case_folder,ID);

% scrivo script pyton 
[done] = PPywriter(TO_BE_OPEN,TO_BE_LOAD,TO_BE_SAVE);
% eseguo
goGoOpenFOAM(sprintf('pvbatch %s',TO_BE_OPEN));
% modifico output in modo da poterlo direttamente importare
system(sprintf('cd %s && tail -n +2 %d0.csv > %dc.csv',case_folder,ID,ID));

air_mat = csvread(sprintf('%s/%dc.csv',case_folder,ID));

%% per calibrazione commentare
%cd = load('current_design.mat');

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
       
    [~,imin] = min(abs([yp1;yp2;yp3;yp4] - ones(4,1)*air_mat(i,10)));
        
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

%%
% reinterpolo sui punti noti da metodo course
x_ventre_air = fliplr(xc(1,1:(size(xc,2)/2)));
x_dorso_air  = xc(1,(size(xc,2)/2)+1:end);

% CP_va = spline(punti_ventre_air(:,8),punti_ventre_air(:,7),x_ventre_air);
% CP_da = spline(punti_dorso_air(:,8),punti_dorso_air(:,7),x_dorso_air);

CP_va = interp1(punti_ventre_air(:,8),punti_ventre_air(:,7),x_ventre_air,'linear','extrap');
CP_da = interp1(punti_dorso_air(:,8),punti_dorso_air(:,7),x_dorso_air,'linear','extrap');


% figure
% plot(punti_ventre_air(:,8),punti_ventre_air(:,7)/(0.5*BU_par.Umag^2),'b'); 
% hold on
% plot(punti_dorso_air(:,8),punti_dorso_air(:,7)/(0.5*BU_par.Umag^2),'r'); 
% plot(x_ventre_air,CP_va/(0.5*BU_par.Umag^2),'bo');
% plot(x_dorso_air,CP_da/(0.5*BU_par.Umag^2),'ro');

cp = [fliplr(CP_va),CP_da];

%%
x_ventre_slat = fliplr(xc(2,1:(size(xc,2)/2)));
x_dorso_slat  = xc(2,(size(xc,2)/2)+1:end);

CP_vs = interp1(punti_ventre_slat(:,8),punti_ventre_slat(:,7),x_ventre_slat,'linear','extrap');
CP_ds = interp1(punti_dorso_slat(:,8),punti_dorso_slat(:,7),x_dorso_slat,'linear','extrap');




% plot(punti_ventre_slat(:,8),punti_ventre_slat(:,7)/(0.5*BU_par.Umag^2),'b'); 
% hold on
% plot(punti_dorso_slat(:,8),punti_dorso_slat(:,7)/(0.5*BU_par.Umag^2),'r'); 
% plot(x_ventre_slat,CP_vs/(0.5*BU_par.Umag^2),'bo');
% plot(x_dorso_slat,CP_ds/(0.5*BU_par.Umag^2),'bo');


cp = [cp,fliplr(CP_vs),CP_ds]/(0.5*BU_par.Umag^2);


if sum(isnan(cp)) > 0
    save('debug.mat')
    error('isnan(cp) > 0!')
end
   

if size(cp,1) < size(cp,2)
    cp = cp';
end

% returnValue è il vettore <0
n = (size(xc,2))/2;

cp_air  = cp(1:(2*n));
cp_slat = cp(2*n+1:end);

D_me = min(cp_air - cp_air(1));   % deltaCp mainElement
D_se = min(cp_slat - cp_slat(1)); % deltaCp slatElement


inequalityC = [-14-D_me , -14-D_se];


%f_cell = {[CL;CD],cp,-[CL/CD],inequalityC,[CL_conv,CD_conv]};
f_cell = {[CL;CD],cp,-CL,inequalityC,[CL_conv,CD_conv]};

save(sprintf('%s/ID_calcolo.mat',case_dir),'X_IN','IN','RES_struct','Parameters','f_cell');

end

