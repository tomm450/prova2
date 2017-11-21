% caller
%clear all
close all
clc

my_dir = what;
my_dir = my_dir.path;
my_dirs = genpath(my_dir);
addpath(my_dirs);


% winner geometry
IN   = [50 250 50 250 50 2 0.02 0.03];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT UTENTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters.n_processori = 3;

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
BU_par.Umag  = 135;         %m/s
BU_par.p     = 0;
%
BU_par.BU_type = 'freestream';%{'freestream','fixedValue'}; %%%%
BU_par.L    = 1;          %m
BU_par.Nu   = 0.000018375; %kg/m s
BU_par.Rho  = 1.225;       %kg/m3
BU_par.alpha = 15;
BU_par.Ux = BU_par.Umag*cos(BU_par.alpha*pi/180);
BU_par.Uz = BU_par.Umag*sin(BU_par.alpha*pi/180);
BU_par.extrusion_Thickness = 0.05; %m

%% HS
Parameters.HSA.npane = 100;
%% MODELLO CFD

STL.point_txt     = 'NACA64212at.txt';
% mesher Gmsh o snappyHexMesh
x_dom     = 20; % semilato quadrato
expRatio  = 1.3;
n_cell_ff = 10;

% GMSH
GM_par.wtd       = 3;
GM_par.solver    = 'gmsh'; 
GM_par.x_dom     = x_dom;
% 
GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff);
GM_par.expRatio  = expRatio;
GM_par.l_airfoil = 0.0008;%GM_par.l_dom/5000;%0.002;
GM_par.l_slat    = 0.0008;%GM_par.l_dom/5000;%0.001;
GM_par.Fstruct   = 0;
GM_par.Fquad     = 1;
% %
% GM_par.rref      = 1; %1.5
% GM_par.lref      = 0.6; %50*GM_par.l_airfoil;
% % %
MESH_par = GM_par;

% OPENFOAM
SOLVER.solver        = 'piso';
SOLVER.startTime     = 0;
SOLVER.endTime       = 25000;
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
OPT.x0 = 0.5*(OPT.lb+OPT.ub);

OPT.Aineq = [3/2, -1,0];
OPT.bineq = -0.02;

OPT.MaxFunctionEvaluations_Data = 100;
OPT.MaxIterations_Data          = 100;
OPT.OptimalityTolerance_Data    = 0.01;
OPT.StepTolerance_Data          = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ACTUAL MM

iter_f = 1;
fvals_f = {};
fvals_c = {};
Skf     = {};

max_iter = 5;

for iter_f = 1:max_iter
    
    % FMIN %%%%
    [x,fval,exitflag,output,lambda,grad,hessian,fvals_copt{iter_f}] = toBeMin(OPT,...
        BU_par,Parameters,GEOM,fvals_f,fvals_c,Skf);
    
    %OPT_res{iter_f} = {x,fval,exitflag,output,lambda,grad,hessian};
    
    copt{iter_f}     = {x,fval,exitflag,output,lambda,grad,hessian};
          
    [xc,dl,theta_G,cp_c] = HSA_core(x,BU_par,Parameters,GEOM,1,10*iter_f,1);
    [CLc,CDc] = aeroCoeff( xc,dl,theta_G,cp_c,BU_par );
    
    course{iter_f}  = {xc,dl,theta_G,cp_c,CLc,CDc};
    fvals_c{iter_f} = cp_c;
        
    [cp_f] = fCl_core(x,IN,CFD,Parameters,xc);
    [CLf,CDf] = aeroCoeff( xc,dl,theta_G,cp_f,BU_par );
    
    fine{iter_f} = {xc,dl,theta_G,cp_f,CLf,CDf};
    fvals_f{iter_f} = cp_f;
       
    h = figure(10*iter_f);
    subplot(1,2,2)
    hold on
    n = 0.5*max(size(xc));
    
    xc_m_foam = [xc(1,1:n),nan,xc(2,1:n)];
    cp_m_foam = [cp_f(1:n)',nan,cp_f(2*n+1:3*n)'];
    
    xc_c_foam = [xc(1,n:end),nan,xc(2,n:end)];        
    cp_c_foam = [cp_f(n:2*n)',nan,cp_f(3*n:end)'];
    
    plot(xc_m_foam,cp_m_foam,'c>','LineWidth',2);
    plot(xc_c_foam,cp_c_foam,'m>','LineWidth',2);
    
    grid on
    set(gca,'YDir','Reverse'); % reverse y axis
    %alfa=num2str(alpha*180/pi);
    tit1 = 'COEFFICIENTE DI PRESSIONE';
    tit2 = ['Incidenza: ',BU_par.alpha, ' gradi'];
    title({tit1;tit2})
    xlabel('corda');
    ylabel('-C_p');
    legend('Dorso HS','Ventre HS','Dorso CFD','Ventre CFD')
    
    %% correzione
    if iter_f == 1  % ho calcolato un fvals_c e un fvals_f, Skk e Skf saranno eye
        Skf = {eye(max(size(fvals_f{end})))};

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
            %             DISf(:,w) = [ISsol_f{end}(1) - ISsol_f{end-w}(1)];
            %             DISc(:,w) = [ISsol_c{end}(1) - ISsol_c{end-w}(1)];
            %
            %             DHf(:,w)  = [Hf{end} - Hf{end-w}];
            %             DHc(:,w)  = [Hc{end} - Hc{end-w}];
            
            
        end
        
        fprintf('DF(:,w) = [fvals_f{end} - fvals_f{end-w}];\nDC(:,w) = [fvals_c{end} - fvals_c{end-w}];\n')
        
        %DC
        
        %DF
        
        [Uf,Sf,Vf] = svd(DC);
        %         [Us,Ss,Vs] = svd(DISc);
        %         [Uh,Sh,Vh] = svd(DHc);
        
        
        % pseudoinversa
        Sf_cros = Sf';
        
        for j = 1: min([size(Sf_cros,1),size(Sf_cros,2)])
            if abs(Sf_cros(j,j)) >= 1e-10
                Sf_cros(j,j) = 1/Sf_cros(j,j);
            end
        end
        
        %         Ss_cros = Ss';
        %         for j = 1: min([size(Ss_cros,1),size(Ss_cros,2)])
        %             if abs(Ss_cros(j,j)) >= 1e-10
        %                 Ss_cros(j,j) = 1/Ss_cros(j,j);
        %             end
        %         end
        %
        %         Sh_cros = Sh';
        %         for j = 1: min([size(Sh_cros,1),size(Sh_cros,2)])
        %             if abs(Sh_cros(j,j)) >= 1e-10
        %                 Sh_cros(j,j) = 1/Sh_cros(j,j);
        %             end
        %         end
               
        DCc   = Vf*Sf_cros*Uf';
        %         DIScc = Vs*Ss_cros*Us';
        %         DHcc = Vh*Sh_cros*Uh';
        
        Skf{end+1} = DF*DCc;
        %         Skk{end+1} = DISf*DIScc;
        %         Skh{end+1} = DHf*DHcc;
        
        
        %         Skf{end}
        %         Skk{end}
        %         Skh{end}
        
    end
    
    
    
end
