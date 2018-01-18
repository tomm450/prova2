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

% winner geometry
IN = [50 250 50 250 50 2 0.02 0.03];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT UTENTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airfoil = load('NACA64212at.txt');

% costruzione geometria
C_mesh        = 1000; %[mm];
Parameters.Geom.d_geometry    = 500;
Parameters.Geom.toll_geometry = 1e-6;
Parameters.Dir_path.my_dir    = my_dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%BU_par.alpha_v = [10 12 14 16 18];
BU_par.alpha_v = [16];
BU_par.wall_function = 1;

BU_par.extrusion_Thickness = 0.05; %m

%% MODELLO CFD
STL.point_txt     = 'NACA64212at.txt';

% mesher Gmsh o snappyHexMesh 

% grandezze in unità di corda
x_dom     = 30; % semilato quadrato
% %n_cell_ff =  10; % numero celle su semilato
expRatio  = 1.1;
n_cell_ff = 25;

% % % GMSH
GM_par.wtd       = 0;
GM_par.solver    = 'gmsh';
GM_par.x_dom     = x_dom;
% %
GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff);
GM_par.expRatio  = expRatio;
GM_par.BL        = 0; % 0 = no; 1 = native; 2 = addlayer openFoam 
GM_par.l_airfoil_v = [0.001 0.00075];

%GM_par.l_slat    = 0.0008;%GM_par.l_dom/5000;%0.001;
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
SOLVER.startTime     = 0;
SOLVER.endTime       = 5000;
SOLVER.deltaT        = 1;
SOLVER.writeInterval = floor(SOLVER.endTime/25);

%


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

%for i = 1:size(GM_par.l_airfoil_v,2)
for i = 1:size(GM_par.l_airfoil_v,2)
    for j = 1:size(BU_par.alpha_v,2)
       
    fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
        
    GM_par.BL = 2;

    fprintf('Caso %d su %d \n',i,size(GM_par.l_airfoil_v,2))    
    fprintf('Sottocaso %d su %d \n',j,size(BU_par.alpha_v,2))   
    xc = 'dummy';
    GM_par.l_airfoil = GM_par.l_airfoil_v(i);
    GM_par.l_slat    = GM_par.l_airfoil_v(i);
    BU_par.alpha = BU_par.alpha_v(j);
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

%     GM_par.ref_method = {'none'};
%     GM_par.par_method = {0};

    GM_par.ref_method = {'clock_simple'};%,'wake'};
    
     GM_par.par_method{1} = {1,10*GM_par.l_airfoil_v(i)};
     %GM_par.par_method{2} = {0.5,2,2,100*GM_par.l_airfoil_v(i)};
    
       
    MESH_par     = GM_par;
    CFD.MESH_par = MESH_par;
    
    CFD.STL           = STL;
    CFD.BU_par        = BU_par;
    CFD.SOLVER        = SOLVER;
    
    [cp_f{i,j}] = fCl_core(OPT.x0,IN,CFD,Parameters,xc);
    

    
    
    end
end
