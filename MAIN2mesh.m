%caller
%clear all
close all
clc

my_dir  = what;

% AGGIUNGO CARTELLE
my_dir  = my_dir.path;
my_dirs = genpath(strcat(my_dir,'/Routines'));
addpath(my_dirs);
my_dirs = genpath(strcat(my_dir,'/ExperimentalData'));
addpath(my_dirs);
my_dirs = genpath(strcat(my_dir,'/Airfoil_preparation'));
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
%https://turbmodels.larc.nasa.gov/naca0012_val.html%

diary log.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DA TESTARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%l_airfoil_v  = [0.2,0.015,0.01,0.009,0.007,0.005,0.004,0.0035];
l_airfoil_v  = [0.01,0.009,0.007,0.005,0.004,0.0035,0.003,0.001];
%l_airfoil_v  = [0.0035,0.0033,0.003,0.0025,0.002];

n_cell_ff_v  = 15;
T_model_c       = {'ko'};%','sa'};%, ko
% grandezze in unità di corda
x_dom_v     = [300]; %30:10:60; % semilato quadrato
alpha_v     = [15];
GM_par.wtd       = 3;

SOLVER.endTime   = 9000;
BU_par.wall_function = 1;
%BL_vect = [2];
BU_par.RE    = 6e6;
GM_par.BL = 2;%BL_vect(z); % 0 = no; 1 = native; 2 = addlayer openFoam


GM_par.Fstruct   = 1;
GM_par.Fquad     = 1;
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

BU_par.p     = 0;

BU_par.BU_type = 'freestream';%{'freestream','fixedValue'}; %%%%

BU_par.L    = 1; %m ( caso incomprimibile tenere sempre 1)

BU_par.Mu   = 0.000018375;         %kg/m s
BU_par.Rho  = 1;                   %kg/m3
BU_par.Nu   = BU_par.Mu/BU_par.Rho;%m2/s 

BU_par.Umag = BU_par.RE*BU_par.Nu/BU_par.L;

BU_par.extrusion_Thickness = 0.05; %m

%fprintf('Re = %d \n\n',(BU_par.Rho*BU_par.L*)/(BU_par.Nu));

%% MODELLO CFD
STL.point_txt     = 'NACA64212at.txt';

% mesher Gmsh o snappyHexMesh 
% %n_cell_ff =  10; % numero celle su semilato
expRatio    = 1.3;

% % % GMSH

GM_par.solver    = 'gmsh';

GM_par.expRatio  = expRatio;



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
%OPT.x0 = [-62 -6 -33]
OPT.x02 = [-16.032 -5.675 -39.993]
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
case_displayed = 0;

for z = 1:size(T_model_c,2)
    for i = 1:size(l_airfoil_v,2)
        for j = 1:size(alpha_v,2)
            for k = 1:size(n_cell_ff_v,2)
                for w = 1:size(x_dom_v,2)
                    
                    
                    
                    
                    fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
                    fprintf('lTmod = %s; caso %d/%d \n',T_model_c{z},z,max(size(T_model_c)));
                    fprintf('l_air = %f; caso %d/%d \n',l_airfoil_v(i),i,max(size(l_airfoil_v)));
                    fprintf('alpha = %f; caso %d/%d \n',alpha_v(j),j,max(size(alpha_v)));
                    fprintf('n_ff  = %f; caso %d/%d \n',n_cell_ff_v(k),k,max(size(n_cell_ff_v)));
                    fprintf('x_dom = %f; caso %d/%d \n',x_dom_v(w),w,max(size(x_dom_v)));
                    
                    
                    GM_par.x_dom     = x_dom_v(w);
                    GM_par.l_dom     = GM_par.x_dom/(2*n_cell_ff_v(k));
                    
                    SOLVER.T_model = T_model_c{z};
                    xc = 'dummy';
                    
                    GM_par.l_airfoil = l_airfoil_v(i);
                    GM_par.l_slat    = l_airfoil_v(i);%/4; %**** da caso 171
                    BU_par.alpha = alpha_v(j);
                    BU_par.Ux    = BU_par.Umag*cos(BU_par.alpha*pi/180);
                    BU_par.Uz    = BU_par.Umag*sin(BU_par.alpha*pi/180);
                    
                    
                     %GM_par.ref_method = {'none'};
                     %GM_par.par_method{1} = {0};
                    
                   GM_par.ref_method = {'wake2'};
                   GM_par.par_method{1} = {1,0,alpha_v(j),3,max([0.05,10*l_airfoil_v(i)])};
                    
                    
%                     GM_par.ref_method = {'clock_simple','wake'};
%                     GM_par.par_method{1} = {0.7,20*l_airfoil_v(i)};
%                     GM_par.par_method{2} = {0.3,1.5,2,50*l_airfoil_v(i)};
                    
%                     GM_par.ref_method = {'sublinear'};
%                     GM_par.par_method{1} = {1,0.7,10,0.1,5};
%                     
%                     GM_par.ref_method = {'sublinearwake'};
%                     GM_par.par_method{1} = {0.5,2+(2/3*x_dom_v(w)*sind(alpha_v(j))),1,20*l_airfoil_v(i),0,0.20,8};
%                     
%                     GM_par.ref_method = {'wake'};
%                     GM_par.par_method{1} = {0.5,2.6+(2/3*x_dom_v(w)*sind(alpha_v(j))),2,50*l_airfoil_v(i)};
                    
                    MESH_par     = GM_par;
                    CFD.MESH_par = MESH_par;
                    
                    CFD.STL           = STL;
                    CFD.BU_par        = BU_par;
                    CFD.SOLVER        = SOLVER;
                    
                    [cp_f{i,j}] = fCl_core_OPT(OPT.x02,IN,CFD,Parameters,nan);
                    %[cp_f{i,j}] = fCl_core_OPT('30p30n',IN,CFD,Parameters,nan);
                    %[cp_f{i,j}] = fCl_core_OPT('NACA0012',IN,CFD,Parameters,nan);
                    
                    if case_displayed == 5
                        clc
                        case_displayed = 0;
                    end
                    
                    case_displayed = case_displayed+1;
                end
            end
        end
    end
end
diary off