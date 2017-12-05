function [Cl_history,Cd_history,t_history] = readClsimple(CD,wtr,dir)

if nargin == 1
    wtr = [0:3];
    dir = '';
elseif nargin == 2
    dir = '';
end


close all

leg_cell = {};
leg_cell2 = {};
n_plot = 0;
passo = 5;
%CD = [4:6];

whatIdid = [];
whatCase = [];
whatCoef = [];
whatTime = [];

for i = 1:size(CD,2)
    
    
    dirProbe = strcat('./',dir,'Cases_folder/',num2str(CD(i)));
    
    if exist(dirProbe,'dir') > 0
        n_plot = n_plot +1;
        
        b = load(sprintf('./%s/ID_calcolo.mat',dirProbe));
        whatIdid = [whatIdid,b.RES_struct.MESH_par.wtd];
        whatCase = [whatCase,CD(i)];
        
        for wtr_index = 1:size(wtr,2)
            if b.RES_struct.MESH_par.wtd == wtr(wtr_index)
                
                % vado a leggere
                fprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat \n',dirProbe)
                [win,mat] = system(sprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat',dirProbe));
                
                mat = str2num(mat);
                [win,exetim] = system(sprintf('cat %s/3logsimple.txt | grep "ClockTime = " | tail -n 1 | cut -d " " -f8',dirProbe));
                exetim = str2num(exetim); exetim = exetim/60;
                           
                fprintf('checkMesh...\n');
                [winCheck] = goGoOpenFOAM(sprintf('cd ./%s/30simple && checkMesh > ../logcheckMesh.txt',dirProbe));
                [a,nface] = system(sprintf('cat ./%s/logcheckMesh.txt | grep ''    faces:            ''',dirProbe)); nface = strsplit(nface,' '); nface = str2num(nface{end});
                [a,npoin] = system(sprintf('cat ./%s/logcheckMesh.txt | grep ''    points:           ''',dirProbe)); npoin = strsplit(npoin,' '); npoin = str2num(npoin{end});
                [a,meshres] = system(sprintf('tail -n 5 ./%s/logcheckMesh.txt',dirProbe));
                fprintf('%s\n',meshres);
                
                ps = round(0.3*size(mat,1));
                plt_style = pl_style_fun(n_plot);
                iterV = ps+[0:size(mat(ps:passo:end,4),1)-1]*passo;
                
                figure(2500); plot(iterV,mat(ps:passo:end,4),plt_style); hold on; grid on
                Sdev = std(mat(end-round(0.1*size(mat,1)):end,4));
                
                
                whatCoef_ad = [mat(end,4) - Sdev; mat(end,4);mat(end,4) + Sdev];
                whatCoef = [whatCoef, whatCoef_ad];
                whatTime = [whatTime,exetim];
                %% RES
                [a,resUx] = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for Ux'' | cut -d'' '' -f9 | tr -d '',''',dirProbe)); resUx = str2num(resUx);
                [a,resUz] = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for Uz'' | cut -d'' '' -f9 | tr -d '',''',dirProbe)); resUz = str2num(resUz);
                [a,resP]  = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for p'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));  resP  = str2num(resP);
                [a,resO]  = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for omega'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));  resO  = str2num(resO);
                [a,resK]  = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for k'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));      resK  = str2num(resK);
                
                
                figure(2501);
                subplot(1,5,1); semilogy(resUx(ps:passo:end),plt_style); title('Ux');hold on; grid on
                subplot(1,5,2); semilogy(resUz(ps:passo:end),plt_style); title('Uz');hold on; grid on
                subplot(1,5,3); semilogy(resP(ps:passo:end),plt_style);  title('p') ;hold on; grid on
                subplot(1,5,4); semilogy(resO(ps:passo:end),plt_style);  title('Omega') ;hold on; grid on
                subplot(1,5,5); semilogy(resK(ps:passo:end),plt_style);  title('k') ;hold on; grid on
                
                
                if isfield('b.RES_struct.MESH_par','Fquad') == 0
                    Quad = b.RES_struct.MESH_par.Fstruct;
                else
                    Quad = b.RES_struct.MESH_par.Fquad;
                end
                
                leg_cell{end+1} = sprintf('ID = %d; n_{point} = %d; n_{face} = %d; x_{dom} = %d; l_{dom} = %1.3f; l_{airfoil} = %1.4f; Fstruct = %d; Quad = %d ',CD(i),...
                    npoin,nface,b.RES_struct.MESH_par.x_dom, b.RES_struct.MESH_par.l_dom ,b.RES_struct.MESH_par.l_airfoil ,b.RES_struct.MESH_par.Fstruct,Quad);
                %leg_cell2{end+1} = sprintf('ID = %d',CD(i));
                %         leg_cell{end+1} = sprintf('ID = %d (tol up)',CD(i));
                %         leg_cell{end+1} = sprintf('ID = %d (tol dwn)',CD(i));
                %leg_cell{end+1} = sprintf('xd = %d; nf = %d; na = %d; it = %d; struct = %d; crono = [%f %f]',xd(i),nf(j),na(k),itend(w),struct(h),crono(1),crono(2));
                
                %         figi = 1000*i+nf(j)*10;
                %         figure(figi)
                %         plot(na(k),mat(end,4),'o'); hold on; grid on;title(sprintf('Cl @ xd = %d; nf = %d;',xd(i),nf(j)));
                Cl_history{n_plot} = mat(:,4);
                Cd_history{n_plot} = mat(:,3);
                t_history{n_plot}  = mat(:,1);
                fprintf('CL = %f; std(CL(200:end)) = %f; exetime = %2.2f min \n',mat(end,4),Sdev,exetim);
            end
        end
    end
end


%for f = 1:2
if size(whatCase,2) > 0 

    figure(2500); legend(leg_cell);
    figure(2501); legend(leg_cell);
    
    figure(30); plot(whatCase,whatIdid,'o'); grid on; title('wtd'); xlabel('Case folder');
    figure(40); plot(whatCase,whatCoef(1,:),'rx--',whatCase,whatCoef(2,:),'bo-',whatCase,whatCoef(3,:),'rx--');
    grid on;title('C_d'); xlabel('Case folder');
    figure(50); plot(whatCase,whatTime,'o');
    grid on;title('Time'); xlabel('Case folder');
    
    
%end
end
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
