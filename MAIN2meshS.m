% caller
%clear all
close all
clc

my_dir  = what;
my_dir  = my_dir.path;
my_dirs = genpath(my_dir);
addpath(my_dirs);

Parameters.gmsh_cmd     = '/usr/bin/gmsh';             % fisso
%Parameters.gmsh_cmd     = '/home/tom/Documents/gmsh'; % portatile
Parameters.n_processori = 8;
if exist(Parameters.gmsh_cmd,'file')
    % ok
else 
    error('percorso Gmsh non esiste')
end


%% TO BE TESTED
l_airfoil_v = [0.005];
alpha_v     = 14; %[10 12 14 16 18];
n_cell_ff_v = [200]; % divisioni su semiliato
x_dom_v     = 20;   % semilato quadrato
max_ref_v   = 6;
wf          = [1];

%%


% winner geometry
IN = [50 250 50 250 50 2 0.02 0.03];
Airfoil = load('NACA64212at.txt');

% costruzione geometria
C_mesh        = 1000; %[mm];
Parameters.Geom.d_geometry    = 500;
Parameters.Geom.toll_geometry = 1e-6;
Parameters.Dir_path.my_dir    = my_dir;

% interpolo dati grezzi e applico approssimazione elisse
% (out escono unitari)
[ up, dwn, Airfoil ] = interfoil(  Airfoil,Parameters.Geom.d_geometry,1 );

% aggiorno Parameters
Parameters.Airfoil.up  = C_mesh*up;
Parameters.Airfoil.dwn = C_mesh*dwn;
Parameters.Airfoil.c   = C_mesh;

Parameters.Airfoil.Mom_ref_x = 0.25*C_mesh;
Parameters.Airfoil.Mom_ref_y = 0;

Parameters.HSA.npane = 100;

% calcolo geometria effettiva

[GEOM,fail,log] = geometry_main(Parameters,IN);

%% BOUNDARY CONDITION

BU_par.Umag  = 45;         %m/s
BU_par.p     = 0;

BU_par.BU_type = 'freestream';%{'freestream','fixedValue'}; %%%%

BU_par.L    = 1;          %m

BU_par.Nu   = 0.000018375; %kg/m s
BU_par.Rho  = 1.225;       %kg/m3

BU_par.extrusion_Thickness = 0.05; %m

%% MODELLO CFD
STL.point_txt     = 'NACA64212at.txt';

% mesher Gmsh o snappyHexMesh 


% % % % SNAPPY
SN_par.solver    = 'snap';
SN_par.wtd = 0;
SN_par.expRatio  = 1.1;
SN_par.wakeBox = 0;


SN_par.simplegrading.flag = 0;
SN_par.LEBox = 1;
SN_par.TEBox = 1;

% % % % GMSH
% GM_par.wtd  = 0;
% GM_par.solver    = 'gmsh';
% GM_par.expRatio  = 1.1;
% GM_par.Fstruct   = 0;
% GM_par.Fquad     = 1;
% 
% % GM refinement
% % casi attualmente implementati
% %     case {'none'}
% %         method_par = {};
% %         done = 1;
% %     case {'clock_simple'}
% %         % 24 punti su circonferenza di raggio rref di dimensione lref
% %         rref = method_par{1};
% %         lref = method_par{2};
% %     case 'sublinear'
% %         % parametri per forma ellisse
% %         ellx = method_par{1};
% %         elly = method_par{2};
% %         
% %         % raggio interpolazione lineare
% %         rref  = method_par{3}; 
% %         % coefficiente tc nuova lunghezza a rref sia lref*l_linear
% %         lref  = method_par{4};
% %         % numero circonferenze da marchiare
% %         nstaz = method_par{5};
% %    case 'wake'
% %       posizione verticale centro raccordo
% %       y_c = method_par{1};
% %       altezza scia all'outlet
% %       y_w = method_par{2};
% %       raggio raccordo 
% %       r   = method_par{3};
% %       lunghezza elementi
% %       l_w = method_par{4};

% % OPENFOAM
SOLVER.solver        = 'simple';
SOLVER.startTime     = 0;
SOLVER.endTime       = 5000;
SOLVER.deltaT        = 1;
SOLVER.writeInterval = floor(SOLVER.endTime/25);

%% OPT parametri

OPT.nvars = 3;
OPT.lb = [-100 -15 -40];
OPT.ub = [   0  15 -5];
OPT.x0 = 0.5*(OPT.lb+OPT.ub);

OPT.Aineq = [3/2, -1,0];
OPT.bineq = -0.02;

OPT.MaxFunctionEvaluations_Data = 100;
OPT.MaxIterations_Data = 100;
OPT.OptimalityTolerance_Data = 0.01;
OPT.StepTolerance_Data = 0.01;

%% ACTUAL MM

iter_f = 1;
fvals_f = {};
fvals_c = {};
Skf     = {};

max_iter = 5;

% [xc,dl,theta_G,cp] = HSA_core_light([-50 2.5 -40],BU_par,xp1,yp1,xp2,yp2,1,1,0);     
% return

cp_f = {};
xc = 'dummy';

for i = 1:size(l_airfoil_v,2)
    for j = 1:size(alpha_v,2)
        for k = 1:size(n_cell_ff_v,2)
            for w = 1:size(x_dom_v,2)
                for z = 1:size(max_ref_v,2)
                    for q = 1:size(wf,2)
                        
                        
                        
                        fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
                        
                        fprintf('l_air              caso %02d/%02d \n',i,size(l_airfoil_v,2));
                        fprintf('  alpha            caso %02d/%02d \n',j,size(alpha_v,2));
                        fprintf('    n_ff           caso %02d/%02d \n',k,size(n_cell_ff_v,2));
                        fprintf('      x_dom        caso %02d/%02d \n',w,size(x_dom_v,2));
                        fprintf('        max ref    caso %02d/%02d \n',z,size(max_ref_v,2));
                        fprintf('          wallfunc caso %02d/%02d \n',q,size(wf,2));
                        
                        BU_par.wall_function = wf(q);
                        BU_par.alpha = alpha_v(j);
                        
                        BU_par.Ux    = BU_par.Umag*cos(BU_par.alpha*pi/180);
                        BU_par.Uz    = BU_par.Umag*sin(BU_par.alpha*pi/180);
                        
                        %% GMSH
                        %GM_par.x_dom     = x_dom_v(w);
                        %GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff_v(k));
                        
                        %GM_par.BL = BL_v(z);
                        
                        
                        %GM_par.l_airfoil = l_airfoil_v(i);
                        %GM_par.l_slat    = l_airfoil_v(i);
                        
                        
                        %     case {'clock_simple'}
                        %         % 24 punti su circonferenza di raggio rref di dimensione lref
                        %         rref = method_par{1};
                        %         lref = method_par{2};
                        %         nstaz = method_par{5};
                        %    case 'wake'
                        %       posizione verticale centro raccordo
                        %       y_c = method_par{1};
                        %       altezza scia all'outlet
                        %       y_w = method_par{2};
                        %       raggio raccordo
                        %       r   = method_par{3};
                        %       lunghezza elementi
                        %       l_w = method_par{4};
                        
                        %     GM_par.ref_method = {'none'};
                        %     GM_par.par_method = {0};
                        
                        %GM_par.ref_method = {'none'};%{'clock_simple'};%,'wake'};
                        %GM_par.par_method = {0};
                        %GM_par.par_method{1} = {1,10*l_airfoil_v(i)};
                        %GM_par.par_method{2} = {0.5,2,2,100*l_airfoil_v(i)};
                        
                        
                        %% MESH_par     = GM_par;
                        
                        SN_par.l_airfoil = l_airfoil_v(i);
                        SN_par.xMax    = x_dom_v(w);
                        SN_par.xMin    = -x_dom_v(w);
                        SN_par.zMax    = x_dom_v(w);
                        SN_par.zMin    = -x_dom_v(w);
                        SN_par.yMax    =  0.1;
                        SN_par.yMin    = -0.1;
                        SN_par.nx      = n_cell_ff_v(k);
                        SN_par.nz      = n_cell_ff_v(k);
                        SN_par.ny      = 2;
  %                [ 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 11, 12, 13, 14, 15, 16, 17]
  SN_par.refDist = [0.4]% , 3 , 6 , 8, 11 , 14, 20, 30, 60];
  
  %SN_par.refDist = [3,9,15,25]
  
                        SN_par.deltabox = 1;
                        SN_par.MaxrefFactor = max_ref_v(z);
                        SN_par.MinrefFactor = max_ref_v(z)-1;
                        
                        
                        MESH_par     = SN_par;
                        
                        
                        CFD.MESH_par = MESH_par;
                        
                        CFD.STL           = STL;
                        CFD.BU_par        = BU_par;
                        CFD.SOLVER        = SOLVER;
                        
                        [cp_f] = fCl_core(OPT.x0,IN,CFD,Parameters,xc);
                        
                    end
                end
                
            end
        end
    end
    
end
