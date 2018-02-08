% caller
%clear all
close all
clc

my_dir  = what;
my_dir  = my_dir.path;
my_dirs = genpath(my_dir);
addpath(my_dirs);

%Parameters.gmsh_cmd     = '/usr/bin/gmsh';             % fisso
Parameters.gmsh_cmd     = '/home/tom/Documents/gmsh'; % portatile
Parameters.n_processori = 4;
if exist(Parameters.gmsh_cmd,'file')
    % ok
else 
    error('percorso Gmsh non esiste')
end

% winner geometry
IN = [50 250 50 250 50 2 0.02 0.03];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DA TESTARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_airfoil_v  = [0.005];
n_cell_ff_v  = [10];
% grandezze in unità di corda
x_dom_v     = 10; %30:10:60; % semilato quadrato
alpha_v     = [10];
GM_par.wtd       = 4;
SOLVER.endTime   = 2000;

BL_vect = [2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

BU_par.wall_function = 1;

BU_par.extrusion_Thickness = 0.05; %m

fprintf('Re = %d \n\n',(BU_par.Rho*BU_par.L*BU_par.Umag)/(BU_par.Nu));

%% MODELLO CFD
STL.point_txt     = 'NACA64212at.txt';

% mesher Gmsh o snappyHexMesh 


% %n_cell_ff =  10; % numero celle su semilato
expRatio    = 1.1;

% % % GMSH

GM_par.solver    = 'gmsh';

GM_par.expRatio  = expRatio;


GM_par.Fstruct   = 1;
GM_par.Fquad     = 0;

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
%    case 'wake'
%       posizione verticale centro raccordo
%       y_c = method_par{1};
%       altezza scia all'outlet
%       y_w = method_par{2};
%       raggio raccordo 
%       r   = method_par{3};
%       lunghezza elementi
%       l_w = method_par{4};

% % OPENFOAM
SOLVER.solver        = 'simple';
SOLVER.T_model       = 'sa'; % la, sa, ko
SOLVER.PLT_history   = 1;
SOLVER.startTime     = 0;

SOLVER.deltaT        = 1;
SOLVER.writeInterval = floor(SOLVER.endTime/25);
SOLVER.PLT_history   = 1;
%% OPT parametri

OPT.nvars = 3;
OPT.lb = [-100 -15 -40];
OPT.ub = [   0  15 -5];
%OPT.x0 = 0.5*(OPT.lb+OPT.ub);
OPT.x0 = [-62 -6 -33]

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

for z = 1:size(BL_vect,2)
    for i = 1:size(l_airfoil_v,2)
        for j = 1:size(alpha_v,2)
            for k = 1:size(n_cell_ff_v,2)
                for w = 1:size(x_dom_v,2)
                    
                    
                    clc
                    
                    fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
                    fprintf('l_air = %f; caso %d/%d \n',l_airfoil_v(i),i,max(size(l_airfoil_v)));
                    fprintf('alpha = %f; caso %d/%d \n',alpha_v(j),j,max(size(alpha_v)));
                    fprintf('n_ff  = %f; caso %d/%d \n',n_cell_ff_v(k),k,max(size(n_cell_ff_v)));
                    fprintf('x_dom = %f; caso %d/%d \n',x_dom_v(w),w,max(size(x_dom_v)));
                    
                    
                    GM_par.x_dom     = x_dom_v(w);
                    GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff_v(k));
                    GM_par.BL = BL_vect(z); % 0 = no; 1 = native; 2 = addlayer openFoam
                    
                    xc = 'dummy';
                    
                    GM_par.l_airfoil = l_airfoil_v(i);
                    GM_par.l_slat    = l_airfoil_v(i)/4;
                    BU_par.alpha = alpha_v(j);
                    BU_par.Ux    = BU_par.Umag*cos(BU_par.alpha*pi/180);
                    BU_par.Uz    = BU_par.Umag*sin(BU_par.alpha*pi/180);
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
                    
                        GM_par.ref_method = {'none'};
                        GM_par.par_method{1} = {0};
                    
%                     GM_par.ref_method = {'clock_simple','wake'};
%                     
%                     GM_par.par_method{1} = {1,10*l_airfoil_v(i)};
%                     GM_par.par_method{2} = {0.5,2,2,30*l_airfoil_v(i)};
                    
                    
                    MESH_par     = GM_par;
                    CFD.MESH_par = MESH_par;
                    
                    CFD.STL           = STL;
                    CFD.BU_par        = BU_par;
                    CFD.SOLVER        = SOLVER;
                    
                    %[cp_f{i,j}] = fCl_core_OPT([],IN,CFD,Parameters,nan);
                    [cp_f{i,j}] = fCl_core_OPT('30p30n',IN,CFD,Parameters,nan);
                    
                    
                end
            end
        end
    end
end
