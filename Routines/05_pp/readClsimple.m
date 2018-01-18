function [Cl_history,Cd_history,t_history] = readClsimple(CD,wtr,dir,PLT,writeCSV)

% LEGGO; SALVO IN CSV 

if nargin == 1
    wtr = [0:3];
    dir = 'Cases_folder/';
    PLT = 0;
    writeCSV = 0;
elseif nargin == 2
    dir = 'Cases_folder/';
    PLT = 0;
    writeCSV = 0;
elseif nargin == 3
    PLT = 0;
    writeCSV = 0;
elseif nargin == 4
    writeCSV = 0;
end
% reset plot

leg_cell = {};

%leg_cell2 = {};

n_plot = 0;

passo = 50;

CSVMAT = [];

for i = 1:size(CD,2)
    
    dirProbe = strcat(dir,num2str(CD(i)));
    
    if exist(dirProbe,'dir') > 0
        
        fprintf('Sto leggendo caso %d \n',CD(i));
        
        n_plot = n_plot +1;
        
        b = load(sprintf('./%s/ID_calcolo.mat',dirProbe));
                   
        %         whatCase = [whatCase,CD(i)];
        %         whatIdid = [whatIdid,b.RES_struct.MESH_par.wtd];
                
        %%
        for wtr_index = 1:size(wtr,2) % se non è quel che voglio leggere, ignoro
            
            if b.RES_struct.MESH_par.wtd == wtr(wtr_index) && strcmp(b.RES_struct.MESH_par.solver,'gmsh')
                
                info_vect = nan(1,53);
                info_vect(1) = CD(i);
                
                %% aggiorno vettore info secondo convenzione
                info_vect(2) = b.RES_struct.BU_par.Umag;
                info_vect(3) = b.RES_struct.BU_par.alpha;
                info_vect(4) = b.RES_struct.BU_par.p;
                
                if strcmp(b.RES_struct.BU_par.BU_type,'freestream')
                    info_vect(5) = 0;
                else % 'fixedValue'
                    info_vect(5) = 1;
                end
                info_vect(6)  = b.RES_struct.BU_par.wall_function;
                
                info_vect(7)  = b.RES_struct.MESH_par.x_dom;
                info_vect(8)  = b.RES_struct.MESH_par.expRatio;
                info_vect(9)  = b.RES_struct.MESH_par.wtd;
                info_vect(10) = b.RES_struct.MESH_par.l_dom;
                info_vect(11) = b.RES_struct.MESH_par.l_airfoil;
                info_vect(12) = b.RES_struct.MESH_par.BL;
                info_vect(13) = b.RES_struct.MESH_par.Fstruct;
                
                info_vect(14:28) = 0;
                
                for t = 1:max(size(b.RES_struct.MESH_par.ref_method))
                    
                    if strcmp(b.RES_struct.MESH_par.ref_method{t},'none')
                        info_vect(14) = 1;
                    elseif strcmp(b.RES_struct.MESH_par.ref_method{t},'clock_simple')
                        
                        info_vect(15) = 1;
                        info_vect(16) = cell2mat(b.RES_struct.MESH_par.par_method{t}(1));
                        info_vect(17) = cell2mat(b.RES_struct.MESH_par.par_method{t}(2));
                        
                    elseif strcmp(b.RES_struct.MESH_par.ref_method{t},'sublinear')
                        
                        info_vect(18) = 1;
                        info_vect(19) = cell2mat(b.RES_struct.MESH_par.par_method{t}(1));
                        info_vect(20) = cell2mat(b.RES_struct.MESH_par.par_method{t}(2));
                        info_vect(21) = cell2mat(b.RES_struct.MESH_par.par_method{t}(3));
                        info_vect(22) = cell2mat(b.RES_struct.MESH_par.par_method{t}(4));
                        info_vect(23) = cell2mat(b.RES_struct.MESH_par.par_method{t}(5));
                        
                    elseif strcmp(b.RES_struct.MESH_par.ref_method{t},'wake')
                        
                        info_vect(24) = 1;
                        info_vect(25) = cell2mat(b.RES_struct.MESH_par.par_method{t}(1));
                        info_vect(26) = cell2mat(b.RES_struct.MESH_par.par_method{t}(2));
                        info_vect(27) = cell2mat(b.RES_struct.MESH_par.par_method{t}(3));
                        info_vect(28) = cell2mat(b.RES_struct.MESH_par.par_method{t}(4));                      
                        
                    else
                        error('condizione non supportata')
                    end
                end
                
                if strcmp(b.RES_struct.SOLVER.solver,'simple')
                    info_vect(29) = 0;
                else % 'piso'
                    info_vect(29) = 1;
                end
                
                info_vect(30) = b.RES_struct.SOLVER.endTime;
                info_vect(32) = b.RES_struct.SOLVER.deltaT;
                               
                
                %%
                % vado a leggere PostProcess
                fprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat \n',dirProbe)
                [~,mat] = system(sprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat',dirProbe));
                
                % Time       	Cm           	Cd           	Cl           	Cl(f)        	Cl(r)        
                mat = str2num(mat);
                %Sdev = std(mat(end-round(0.1*size(mat,1)):end,4));
                
                [~,exetim] = system(sprintf('cat %s/3logFoam.txt | grep -a "ClockTime = " | tail -n 1 | cut -d " " -f8',dirProbe));
                
                if size(exetim) == [0 0] % calcolo fallito
                    exetim = 0;
                else
                    exetim = str2num(exetim); exetim = exetim/60;
                end
                
                info_vect(33) = exetim;
                
                % CHECKMESH
                if exist(sprintf('./%s/logcheckMesh.txt',dirProbe),'file')
                    % non devo rilanciare
                else
                    fprintf('checkMesh...\n');
                end
                
                if exist(sprintf('./%s/20extrude',dirProbe),'dir')
                    
                [~] = goGoOpenFOAM(sprintf('cd ./%s/20extrude && checkMesh > ../logcheckMesh.txt',dirProbe));
                
                elseif exist(sprintf('./%s/30simple',dirProbe),'dir')
                    
                    [~] = goGoOpenFOAM(sprintf('cd ./%s/30simple && checkMesh > ../logcheckMesh.txt',dirProbe));
                    
                elseif exist(sprintf('./%s/40piso',dirProbe),'dir')
                    
                    [~] = goGoOpenFOAM(sprintf('cd ./%s/40piso && checkMesh > ../logcheckMesh.txt',dirProbe));
                    
                else 
                    error('Where I am?');
                end
                % nface
                [~,nface] = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    faces:            ''',dirProbe)); nface = strsplit(nface,' '); nface = str2num(nface{end}); 
                % npoint
                [~,npoin] = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    points:           ''',dirProbe)); npoin = strsplit(npoin,' '); npoin = str2num(npoin{end});
                % Max Aspect Ratio
                [~,MAR]   = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    Max aspect ratio = ''',dirProbe));
                
                if size(MAR) == [0 0] % controllo fallito
                    MAR = nan;
                else
                    MAR = strsplit(MAR,' ');
                    MAR = str2num(MAR{6});
                end
                
                % Face Area (min e Max)
                [~,areaCheck]   = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    Minimum face area = ''',dirProbe)); 
                
                if size(areaCheck) == [0 0]
                    mA = nan; MA = nan;
                else
                    areaCheck = strsplit(areaCheck,' '); mA = str2num(areaCheck{6}(1:end-1)); MA = str2num(areaCheck{11}(1:end-1));
                end
                
                % Volume (min e Max)
                [~,volCheck]   = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    Min volume = ''',dirProbe)); 
                
                if size(volCheck) == [0 0]
                    mV = nan; MV = nan;
                else
                volCheck = strsplit(volCheck,' ');   mV = str2num(volCheck{5}(1:end-1)); MV = str2num(volCheck{9}(1:end-1));
                end
                
                % Volume (min e Max)
                [~,ortCheck]   = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    Mesh non-orthogonality Max: ''',dirProbe)); 
                ortCheck = strsplit(ortCheck,' '); MNonOrt = str2num(ortCheck{5}); ANonOrt = str2num(ortCheck{7});
                % Max Skew
                [~,MS] = system(sprintf('cat ./%s/logcheckMesh.txt | grep -a ''    Max skewness = ''',dirProbe)); 
                
                if size(MS) == [0 0]
                    MS = nan; 
                else
                MS = strsplit(MS,' '); MS = str2num(MS{end-1});
                end
                
                info_vect(34:43) = [nface, npoin, MAR, mA,MA,mV,MV,MNonOrt,ANonOrt,MS];
                               
                [~,meshres] = system(sprintf('tail -n 5 ./%s/logcheckMesh.txt',dirProbe));
                fprintf('%s\n',meshres);
                
                if size(mat) == [0 0]
                    % calcolo falito
                    mat = nan(1,4);
                    iterV = nan(1,4);
                    ps = 1;
                    
                    resUx = nan(1,4);
                    resUz = nan(1,4);
                    resP  = nan(1,4);
                    resO  = nan(1,4);
                    resK  = nan(1,4);
                    
                    
                else
                    
                    
                    ps = round(0.3*size(mat,1));
                    
                    iterV = ps+[0:size(mat(ps:passo:end,4),1)-1]*passo;
                    
                    %% RES
                    [a,resUx] = system(sprintf('cat ./%s/3logFoam.txt | grep ''Solving for Ux'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));     resUx = str2num(resUx);
                    [a,resUz] = system(sprintf('cat ./%s/3logFoam.txt | grep ''Solving for Uz'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));     resUz = str2num(resUz);
                    [a,resP]  = system(sprintf('cat ./%s/3logFoam.txt | grep ''Solving for p'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));      resP  = str2num(resP);
                    [a,resO]  = system(sprintf('cat ./%s/3logFoam.txt | grep ''Solving for omega'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));  resO  = str2num(resO);
                    [a,resK]  = system(sprintf('cat ./%s/3logFoam.txt | grep ''Solving for k'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));      resK  = str2num(resK);
                    
                    if size(resUx) ~= size(resP)
                       resP = resP(1:2:end-1);
                    end
                    
                    if size(resO) == [0 0] % caso laminare
                       resO = nan(size(resUx));
                       resK = nan(size(resUx));
                       
                    end
                    Sdev = std(mat(end-round(0.1*size(mat,1)):end,4));
                    info_vect(45) = Sdev;
                end
                
                
                
                % Time       	Cm           	Cd           	Cl           	Cl(f)        	Cl(r)
                info_vect(31) = mat(end,1);
                
                info_vect(44) = mat(end,4);
                
                info_vect(46) = mat(end,3);
                info_vect(47) = mat(end,2);
                
                info_vect(48:52) = [resUx(end) resUz(end) resP(end) resO(end) resK(end)];
                
                info_vect(53)    = mean(mat(ceil(0.99*size(mat,1)):end,4));
                
                plt_style = pl_style_fun(n_plot);
                
                if PLT == 1
                    
                    figure(2500); plot(iterV,mat(ps:passo:end,4),plt_style); hold on; grid on
                                    
                    
                    %                 whatCoef_ad = [mat(end,4) - Sdev; mat(end,4);mat(end,4) + Sdev];
                    %                 whatCoef = [whatCoef, whatCoef_ad];
                    %                 whatTime = [whatTime,exetim];
                    
                    figure(2501);
                    subplot(1,5,1); semilogy(iterV,resUx(ps:passo:end),plt_style); title('Ux');hold on; grid on
                    subplot(1,5,2); semilogy(iterV,resUz(ps:passo:end),plt_style); title('Uz');hold on; grid on
                    subplot(1,5,3); semilogy(resP(ps:passo:end),plt_style);  title('p') ;hold on; grid on %semilogy(iterV,resP(ps:passo:end),plt_style);  title('p') ;hold on; grid on
                    subplot(1,5,4); semilogy(iterV,resO(ps:passo:end),plt_style);  title('Omega') ;hold on; grid on
                    subplot(1,5,5); semilogy(iterV,resK(ps:passo:end),plt_style);  title('k') ;hold on; grid on
                    
                    leg_cell{end+1} = sprintf('ID = %d; l_{airfoil} = %1.4f; alpha = %2.2f; ',CD(i),b.RES_struct.MESH_par.l_airfoil,b.RES_struct.BU_par.alpha); % n_{point} = %d; n_{face} = %d; x_{dom} = %d; l_{dom} = %1.3f; l_{airfoil} = %1.4f; Fstruct = %d; Quad = %d ',CD(i),...
%                     %                     npoin,nface,b.RES_struct.MESH_par.x_dom, b.RES_struct.MESH_par.l_dom ,b.RES_struct.MESH_par.l_airfoil ,b.RES_struct.MESH_par.Fstruct,Quad);
%                                      leg_cell2{end+1} = sprintf('ID = %d',CD(i));
%                                             leg_cell{end+1} = sprintf('ID = %d (tol up)',CD(i));
%                                             leg_cell{end+1} = sprintf('ID = %d (tol dwn)',CD(i));
                                            
                   
                    %                 %         figi = 1000*i+nf(j)*10;
                    %                 %         figure(figi)
                    %                 %         plot(na(k),mat(end,4),'o'); hold on; grid on;title(sprintf('Cl @ xd = %d; nf = %d;',xd(i),nf(j)));
                    %                 Cl_history{n_plot} = mat(:,4);
                    %                 Cd_history{n_plot} = mat(:,3);
                    %                 t_history{n_plot}  = mat(:,1);
                    %                 fprintf('CL = %f; std(CL(200:end)) = %f; exetime = %2.2f min \n',mat(end,4),Sdev,exetim);
                end
                
                CSVMAT = [CSVMAT; info_vect];
            end
        end
    end
    
   
end

if writeCSV == 1
   csvwrite('report.csv',CSVMAT)
end

if PLT == 1
    figure(2500); legend(leg_cell);
    figure(2501); legend(leg_cell);
end


% %for f = 1:2
% if size(whatCase,2) > 0 
% 
%     figure(2500); legend(leg_cell);
%     figure(2501); legend(leg_cell);
%     
%     figure(30); plot(whatCase,whatIdid,'o'); grid on; title('wtd'); xlabel('Case folder');
%     figure(40); plot(whatCase,whatCoef(1,:),'rx--',whatCase,whatCoef(2,:),'bo-',whatCase,whatCoef(3,:),'rx--');
%     grid on;title('C_d'); xlabel('Case folder');
%     figure(50); plot(whatCase,whatTime,'o');
%     grid on;title('Time'); xlabel('Case folder');
%     
%     
% %end
% end
end

function plt_style = pl_style_fun(n_line)

%% differenzio plot
plotc = {'b','g','r','c','m','k'}; % 6
plots = {'o','x','+','*','s','d','v','^','<','>','p','h'}; % 12
plotl = {'-',':','-.','--'}; % 4

ic = 0;
is = 1;
il = 1;

for w = 1:n_line
    ic = ic +1;
    if ic > 6
        is = is+1;
        ic = 1;
    end
    
    if is > 12
        il = il+1;
        is = 1;
    end
end
plt_style = strcat(plotc{ic},plots{is},plotl{il});
end
