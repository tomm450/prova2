clc
close all
clear all


xd     = [20 ];
nf     = [20 ];
na     = [400];

itend  = 6000;%[2500 5000 6000];
struct = [0 1];

leg_cell = {};
n_plot = 0;
passo = 50;

for i = 1:size(xd,2)
    
    for j = 1:size(nf,2)
        
        for k = 1:size(na,2)
            
            for w = 1:size(itend,2)
                
                
                for h = 1:size(struct,2)
                    
                    it_s = sprintf('%d',itend(w));
                    dirProbe = strcat('x_',num2str(xd(i)),'_nfarfiel_',num2str(nf(j)),'_n_air_',num2str(na(k)),'_struct_',num2str(struct(h)),'_',it_s);
                    
                    if exist(dirProbe,'dir') > 0
                        n_plot = n_plot +1;
                        % vado a leggere
                        fprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat \n',dirProbe)
                        [win,mat] = system(sprintf('tail -n +10 ./%s/30simple/postProcessing/forceCoeffs/0/forceCoeffs.dat',dirProbe));
                        
                        mat = str2num(mat);
                        
                        b = load(sprintf('./%s/ID_calcolo.mat',dirProbe));
                        
                        if isfield(b,'crono_cell')
                            crono = b.crono_cell{end};
                        else
                            crono = b.crono;
                        end
                        plt_style = pl_style_fun(n_plot);
                        
                        figure(1); plot(mat(500:passo:end,4),plt_style); hold on; grid on
                         hold on;plot([0 max(size(mat(500:passo:end,4)))],[mat(end,4) mat(end,4)]-0.001,'r--')
                         hold on;plot([0 max(size(mat(500:passo:end,4)))],[mat(end,4) mat(end,4)]+0.001,'r--')
                        %% RES
                        [a,resUx] = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for Ux'' | cut -d'' '' -f9 | tr -d '',''',dirProbe)); resUx = str2num(resUx);
                        [a,resUz] = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for Uz'' | cut -d'' '' -f9 | tr -d '',''',dirProbe)); resUz = str2num(resUz);
                        [a,resP]  = system(sprintf('cat ./%s/3logsimple.txt | grep ''Solving for p'' | cut -d'' '' -f9 | tr -d '',''',dirProbe));  resP  = str2num(resP);
                        
                        figure(2); semilogy(resUx(500:passo:end),plt_style); title('Ux');hold on; grid on
                        figure(3); semilogy(resUz(500:passo:end),plt_style); title('Uz');hold on; grid on
                        figure(4); semilogy(resP(500:passo:end),plt_style);  title('p') ;hold on; grid on
                        
                        
                        
                        leg_cell{end+1} = sprintf('xd = %d; nf = %d; na = %d; it = %d; struct = %d; crono = [%f %f]',xd(i),nf(j),na(k),itend(w),struct(h),crono(1),crono(2));
                        
                        figi = 1000*i+nf(j)*10;
                        figure(figi)
                        plot(na(k),mat(end,4),'o'); hold on; grid on;title(sprintf('Cl @ xd = %d; nf = %d;',xd(i),nf(j)));
                        
                        
                        
                    end
                end
                
                
                
            end
            
        end
        
    end
end
for f = 1:4
    figure(f); legend(leg_cell);
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
