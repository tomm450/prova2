% caller
clear all
close all
clc

clean

my_dir = what;
my_dir = my_dir.path;
my_dirs = genpath(my_dir);
addpath(my_dirs);

diary ./Output/diary.txt

% winner geometry
IN   = [50 250 50 250 50 2 0.02 0.03];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT UTENTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Parameters.n_processori = 8;
Parameters.gmsh_cmd     = 'gmsh';             % fisso
%Parameters.gmsh_cmd     = '~/Documents/gmsh'; % portatile
% carico dati profilo
Airfoil = load('NACA64212at.txt');

% costruzione geometria
C_mesh                        = 1000; %[mm];
Parameters.Geom.d_geometry    = 500;
Parameters.Geom.toll_geometry = 1e-6;
Parameters.Dir_path.my_dir    = my_dir;

% interpolo dati grezzi e applico approssimazione elisse
% (out escono unitari)
[ up, dwn, Airfoil ] = interfoil( Airfoil,Parameters.Geom.d_geometry,1);

% aggiorno Parameters
Parameters.Airfoil.up  = C_mesh*up;
Parameters.Airfoil.dwn = C_mesh*dwn;
Parameters.Airfoil.c   = C_mesh;

Parameters.Airfoil.Mom_ref_x = 0.25*C_mesh;
Parameters.Airfoil.Mom_ref_y = 0;

% calcolo geometria slat&profilo
[GEOM,fail,log] = geometry_main(Parameters,IN);

%% BOUNDARY CONDITION
BU_par.Umag  = 45;         %m/s
BU_par.p     = 0;
%
BU_par.BU_type = 'freestream';%{'freestream','fixedValue'}; %%%%
BU_par.L     = 1;          %m
BU_par.Nu    = 0.000018375; %kg/m s
BU_par.Rho   = 1.225;       %kg/m3
BU_par.alpha = 15;
BU_par.Ux = BU_par.Umag*cos(BU_par.alpha*pi/180);
BU_par.Uz = BU_par.Umag*sin(BU_par.alpha*pi/180);

BU_par.extrusion_Thickness = 0.05; %m

BU_par.wall_function = 1;


%% HS
Parameters.HSA.npane = 150;
%% MODELLO CFD

STL.point_txt = 'NACA64212at.txt';
% mesher Gmsh o snappyHexMesh
x_dom     = 30; % semilato quadrato
expRatio  = 1.1;
n_cell_ff = 10;

% GMSH
GM_par.wtd       = 3;
GM_par.solver    = 'gmsh'; 
GM_par.x_dom     = x_dom;
% 
GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff);
GM_par.expRatio  = expRatio;
GM_par.l_airfoil = 0.005;%GM_par.l_dom/5000;%0.002;
GM_par.l_slat    = 0.005/4;%GM_par.l_dom/5000;%0.001;
GM_par.Fstruct   = 0;
GM_par.Fquad     = 1;

% GM refinement
% casi attualmente implementati
%     case {'none'}
%         method_par = {};
%         done = 1;
%     case {'clock_simple'}
%         % 24 punti su circonferenza di raggio rref di dimensione lref
%         rref = method_par{1};
%         lref = method_par{2};
%     case 'sublinear'
%         % parametri per forma ellisse
%         ellx = method_par{1};
%         elly = method_par{2};
%         
%         % raggio interpolazione lineare
%         rref  = method_par{3}; 
%         % coefficiente tc nuova lunghezza a rref sia lref*l_linear
%         lref  = method_par{4};
%         % numero circonferenze da marchiare
%         nstaz = method_par{5};

% GM_par.ref_method = {'none'};
% GM_par.par_method = {0};

GM_par.ref_method = {'clock_simple','wake'};

GM_par.par_method{1} = {1,10*GM_par.l_airfoil};
GM_par.par_method{2} = {0.5,2,2,30*GM_par.l_airfoil};

GM_par.BL = 2;
% % %
MESH_par = GM_par;

% OPENFOAM
SOLVER.solver        = 'simple';
SOLVER.T_model       = 'la'; % la, sa, ko
SOLVER.PLT_history   = 1;
SOLVER.startTime     = 0;
SOLVER.endTime       = 10000;
SOLVER.deltaT        = 1;
SOLVER.writeInterval = floor(SOLVER.endTime/50);

%
CFD.STL           = STL;
CFD.BU_par        = BU_par;
CFD.MESH_par      = MESH_par;
CFD.SOLVER        = SOLVER;

%% OPT parametri

OPT.nvars = 3;
OPT.lb = [-100 -15 -40];
OPT.ub = [   0  15 -5];

% nstart = 2;
% 
% x0M =  [linspace(OPT.lb(1),OPT.ub(1),nstart+2);...
%         linspace(OPT.lb(2),OPT.ub(2),nstart+2);...
%         linspace(OPT.lb(3),OPT.ub(3),nstart+2)];
%     
% x0M = x0M(:,2:end-1);
% OPT.x0  = zeros(nstart.^OPT.nvars,OPT.nvars);
% riempimento = 1;
% 
% for i =1:size(x0M,2)
%     for j =1:size(x0M,2)
%         for k =1:size(x0M,2)
%             OPT.x0(riempimento,:) = [x0M(1,i),x0M(2,j),x0M(3,k)];
%             riempimento = riempimento + 1;
%         end
%     end
% end
% OPT.x0(end+1,:) = OPT.lb + 0.5*[OPT.ub-OPT.lb];
% 
% clear riempimento x0M 

OPT.x0 = 0.5*(OPT.lb+OPT.ub);

OPT.Aineq = [3/2, -1,0];
OPT.bineq = -0.015*1000;

OPT.MaxFunctionEvaluations_Data = 1000;
OPT.MaxIterations_Data          = 1000;
OPT.OptimalityTolerance_Data    = 1e-4;
OPT.StepTolerance_Data          = 1e-4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ACTUAL MM
MASTER = [];
iter_f = 1;
fvals_f = {}; ISsol_f = {};
fvals_c = {}; ISsol_c = {};
Skf     = {}; Skk     = {};

max_iter = 10;

dx = 1; nonImpr = 0;
          % ripartire? da dove? con che file? da dove?
restart = {0,0,'./Output/1dump_1fmin.mat'};    


while iter_f < max_iter && dx > 1e-3 && nonImpr < 8
    
    if restart{1} == 1 && restart{1} == 1
        load(restart{3});
        restart{1} = 0;
    elseif restart{1} == 1 && restart{1} > 1
        % do nothing
    else
        
        
        [c_cell{iter_f},copt_cell{iter_f},...
            WIN_OPT{iter_f},FULL_OPT{iter_f},usefulCoor{iter_f}] = ...
            toBeMin(OPT,BU_par,Parameters,GEOM,...
            fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk);
    end
    save(sprintf('./Output/%ddump_1fmin.mat',iter_f));
    
    %OPT.x0(end+1,:) =    WIN_OPT{iter_f}{1};
    
    if restart{1} == 1 && restart{1} == 2
        load(restart{3});
        restart{1} = 0;
    else
        [f_cell{iter_f}] = fCl_core_OPT(WIN_OPT{iter_f}{1},IN,CFD,Parameters,usefulCoor{iter_f});
    end
    save(sprintf('./Output/%ddump_2cfd.mat',iter_f));
    %[CLf,CDf] = aeroCoeff( xc,dl,theta_G,cp_f,BU_par );
    
    
    %     fine{iter_f} = {xc,dl,theta_G,cp_f,CLf,CDf};
    %     fvals_f{iter_f} = cp_f;
%        
    
    cp_c    = c_cell{iter_f}{2};
    cp_copt = copt_cell{iter_f}{2};
    cp_f    = f_cell{iter_f}{2};
   
    % PLOT
    h = figure(10*iter_f);
    subplot(1,2,1)
    xc = usefulCoor{iter_f}{1};
    n = 0.5*max(size(xc));
    xc_m_plot = [xc(1,1:n) ,nan,xc(2,1:n)];
    
    cp_m_c    = [cp_c(1:n)',   nan,cp_c(2*n+1:3*n)'];
    cp_m_corr = [cp_copt(1:n)',nan,cp_copt(2*n+1:3*n)'];
    cp_m_foam = [cp_f(1:n)',   nan,cp_f(2*n+1:3*n)'];
    
    xc_c_plot = [xc(1,n:end),    nan,xc(2,n:end)];        
    cp_c_c    = [cp_c(n:2*n)',   nan,cp_c(3*n:end)'];
    cp_c_corr = [cp_copt(n:2*n)',nan,cp_copt(3*n:end)'];
    cp_c_foam = [cp_f(n:2*n)',   nan,cp_f(3*n:end)'];
    
    hold on
    plot(xc_m_plot,cp_m_c,'c>','LineWidth',2);
    plot(xc_c_plot,cp_c_c,'c>','LineWidth',2);
    
    plot(xc_m_plot,cp_m_corr,'g*','LineWidth',2);
    plot(xc_c_plot,cp_c_corr,'g+','LineWidth',2);
    
    plot(xc_m_plot,cp_m_foam,'r>','LineWidth',2);
    plot(xc_c_plot,cp_c_foam,'r>','LineWidth',2);
    
    grid on
    set(gca,'YDir','Reverse'); % reverse y axis
    %alfa=num2str(alpha*180/pi);
    %tit1 = 'COEFFICIENTE DI PRESSIONE';
    %tit2 = ['Incidenza: ',BU_par.alpha, ' gradi'];
    %title({tit1;tit2})
    xlabel('corda');
    ylabel('-C_p');
    legend('Dorso HS','Ventre HS','Dorso HS_corr','Ventre HS_corr','Dorso CFD','Ventre CFD')
    
    subplot(1,2,2)
    semilogy(abs(cp_f-cp_c)); grid on; hold on title('Pressure difference error');
    semilogy(abs(cp_f-cp_copt));
    legend('HS vs CFD','HS_Corr vs CFD');
    
    savefig(gcf,sprintf('./Output/Pdistr_%d.fig',size(fvals_f,2)+1));
    %close gcf

    % PRINT INFO
    fprintf('iter = %d; x = [%2.3f %2.3f %2.3f]\nerr_medio   = %3.6f \nerr_medio_c = %3.6f \n',...
        iter_f,WIN_OPT{iter_f}{1}(1),WIN_OPT{iter_f}{1}(2),WIN_OPT{iter_f}{1}(3),iter_f,mean(abs(cp_f-cp_c)),mean(abs(cp_f-cp_copt)));
    
    fprintf('\nHS   -> CL = %02.5f; CD = %02.5f; fit = %02.5f; ineq = [ %02.5f ; %02.5f] \n',...
        c_cell{iter_f}{1}(1),c_cell{iter_f}{1}(2), c_cell{iter_f}{3}, c_cell{iter_f}{4}(1), c_cell{iter_f}{4}(2));
    fprintf('HS_c -> CL = %02.5f; CD = %02.5f; fit = %02.5f; ineq = [ %02.5f ; %02.5f] \n',...
        copt_cell{iter_f}{1}(1),copt_cell{iter_f}{1}(2), copt_cell{iter_f}{3}, copt_cell{iter_f}{4}(1), copt_cell{iter_f}{4}(2));
    fprintf('CFD  -> CL = %02.5f; CD = %02.5f; fit = %02.5f; ineq = [ %02.5f ; %02.5f] \n',...
        f_cell{iter_f}{1}(1),f_cell{iter_f}{1}(2), f_cell{iter_f}{3}, f_cell{iter_f}{4}(1), f_cell{iter_f}{4}(2));
    
    %% correzione
    
    fvals_c{iter_f}    = c_cell{iter_f}{1};
    %cp_corr = copt_cell{iter_f}{2};
    fvals_f{iter_f}    = f_cell{iter_f}{1};
    
    
    ISsol_c{iter_f}    = c_cell{iter_f}{2};
    %cp_corr = copt_cell{iter_f}{2};
    ISsol_f{iter_f}    = f_cell{iter_f}{2};
    
     
    if iter_f == 1  % ho calcolato un fvals_c e un fvals_f, Skk e Skf saranno eye
        Skf = {eye(max(size(fvals_f{end})))};
        Skk = {eye(max(size(ISsol_f{end})))};
        % Skh = {1};
       
        DF = [];
        DC = [];
        
    else
        clear DF DC
        %for w = 1: max([size(tgt,1), size(xs,1) -1])
        %for w = 1: max([size(tgt,1),(size(fvals_f,2) -1)])
        
        for w = 1: min([OPT.nvars-1,(size(fvals_f,2) -1)])
            
            DF(:,w) = [fvals_f{end} - fvals_f{end-w}];
            DC(:,w) = [fvals_c{end} - fvals_c{end-w}];
            
            DISf(:,w) = [ISsol_f{end} - ISsol_f{end-w}(1)];
            DISc(:,w) = [ISsol_c{end} - ISsol_c{end-w}(1)];
            %
            %             DHf(:,w)  = [Hf{end} - Hf{end-w}];
            %             DHc(:,w)  = [Hc{end} - Hc{end-w}];
            
            
        end
        
        fprintf('DF(:,w) = [fvals_f{end} - fvals_f{end-w}];\nDC(:,w) = [fvals_c{end} - fvals_c{end-w}];\n')
        
        %DC
        
        %DF
        
        [Uf,Sf,Vf] = svd(DC);
        [Us,Ss,Vs] = svd(DISc);
        %         [Uh,Sh,Vh] = svd(DHc);
        
        
        % pseudoinversa
        Sf_cros = Sf';
        
        for j = 1: min([size(Sf_cros,1),size(Sf_cros,2)])
            if abs(Sf_cros(j,j)) >= 1e-10
                Sf_cros(j,j) = 1/Sf_cros(j,j);
            end
        end
        
        Ss_cros = Ss';
        for j = 1: min([size(Ss_cros,1),size(Ss_cros,2)])
            if abs(Ss_cros(j,j)) >= 1e-10
                Ss_cros(j,j) = 1/Ss_cros(j,j);
            end
        end
        %
        %         Sh_cros = Sh';
        %         for j = 1: min([size(Sh_cros,1),size(Sh_cros,2)])
        %             if abs(Sh_cros(j,j)) >= 1e-10
        %                 Sh_cros(j,j) = 1/Sh_cros(j,j);
        %             end
        %         end
               
        DCc   = Vf*Sf_cros*Uf';
        DIScc = Vs*Ss_cros*Us';
        %         DHcc = Vh*Sh_cros*Uh';
        
        Skf{end+1} = DF*DCc;
        Skk{end+1} = DISf*DIScc;
        %         Skh{end+1} = DHf*DHcc;
        
        
        %         Skf{end}
        %         Skk{end}
        %         Skh{end}
        
    end
    
    
    
    MASTER(iter_f) = f_cell{iter_f}{3}; % cifra minimizzata

    %% aggiorno
    if iter_f <= 2
        dx = 1; nonImpr = 0;
    else
        dx   = norm(WIN_OPT{end}{1}-WIN_OPT{end-1}{1});
        
    end
    
    if MASTER(end) == min(MASTER)
        nonImpr = 0;
    else
        nonImpr = nonImpr+1;
    end
    save(sprintf('./Output/%ddump.mat',iter_f));
    iter_f = iter_f+1;
end

diary off

