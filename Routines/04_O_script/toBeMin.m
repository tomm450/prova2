function [c_cell,copt_cell,WIN_OPT,FULL_OPT,usefulCoor] = toBeMin(OPT,BU_par,Parameters,GEOM,...
              fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk)
% old output 
%x,fval,exitflag,output,lambda,grad,hessian,cp_copt
          
%% Modify options setting
trial = size(OPT.x0,1);   
%nvars = OPT.nvars;
lb    = OPT.lb;
ub    = OPT.ub;
x0    = OPT.x0;

Aineq = OPT.Aineq;
bineq = OPT.bineq;

MaxFunctionEvaluations_Data = OPT.MaxFunctionEvaluations_Data;
MaxIterations_Data          = OPT.MaxIterations_Data;
OptimalityTolerance_Data    = OPT.OptimalityTolerance_Data;
StepTolerance_Data          = OPT.StepTolerance_Data;

% Start with the default options
options = optimoptions('fmincon');

options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunctionEvaluations', MaxFunctionEvaluations_Data);
options = optimoptions(options,'MaxIterations', MaxIterations_Data);
options = optimoptions(options,'OptimalityTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'FunctionTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'StepTolerance', StepTolerance_Data);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize });


%% Preparo Geometria
% interpolo punti in modo da avere 2*Parameters.HSA.npane pannelli per ogni 
% superfice
G_pane_airfoil = [flipud(GEOM.dwn_land)  ; GEOM.up_land(2:end,:)];
G_pane_slat    = [flipud(GEOM.slat_u_dwn); GEOM.slat_u_up(2:end,:)];

% inizializzo con vettore completo
xp1 = G_pane_airfoil(:,1)';
xp2 = G_pane_slat(:,1)';

[~,G_pane_air_LE]  = min(G_pane_airfoil(:,1));
[~,G_pane_slat_LE] = min(G_pane_slat(:,1));

index1 = round(linspace(1,G_pane_air_LE,Parameters.HSA.npane+1));
index1 = [index1(1:end-1),round(linspace(G_pane_air_LE,size(xp1,2),Parameters.HSA.npane+1))];

index2 = round(linspace(1,G_pane_slat_LE,Parameters.HSA.npane+1));
index2 = [index2(1:end-1),round(linspace(G_pane_slat_LE,size(xp2,2),Parameters.HSA.npane+1))];

xp1 = xp1(index1); yp1 = G_pane_airfoil(index1,2)';
xp2 = xp2(index2); yp2 = G_pane_slat(index2,2)';

xp = [xp1;xp2];%/1000;
yp = [yp1;yp2];%/1000;

%% definisco ingresso fmincon

actual_cost_function = @(x) ...
    switchFunCostr(x,BU_par,xp1,yp1,xp2,yp2,...
                   fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,'FitnessValue');
               
nonlinear_constr     = @(x) ...
    switchFunCostr(x,BU_par,xp1,yp1,xp2,yp2,...
                   fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,'Constain'); 

% % % for i = 1:trial 
% % % 
% % % MAT(:,i) = nonlinear_constr(x0(i,:))%     = @(x) ...
% % %     %switchFunCostr(x,BU_par,xp1,yp1,xp2,yp2,...
% % %     %               fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,'Constain'); 
% % % 
% % % end
tstart = tic;

% inizializzo compoenti multistart
fval    = nan(1,trial);  % vettore
exitflag= nan(1,trial); 

x       = cell(1,trial); 

output  = cell(1,trial);
lambda  = cell(1,trial);
grad    = cell(1,trial);
hessian = cell(1,trial);

for ms = 1:trial
    
    fprintf('Multistart caso %d/%d \n\n',ms,size(OPT.x0,1));
    
    [x{ms},fval(ms),exitflag(ms),output{ms},lambda{ms},grad{ms},hessian{ms}] = ...
        fmincon(actual_cost_function,x0(ms,:),Aineq,bineq,[],[],lb,ub,nonlinear_constr,options);
    optfig = gcf;
    savefig(optfig,sprintf('./Output/Optifig_%d_%d.fig',size(fvals_f,2)+1,ms));
    close gcf



end

% verifico migliore convergenza
    fprintf('exitflag = [\n');
    fprintf('%d ',exitflag);
    fprintf(']\n');

[MinTrial,iMinTrial] = min(fval+(100*(exitflag== -2)));

telapsed = toc(tstart);


x_win = [x{iMinTrial}(1),x{iMinTrial}(2),x{iMinTrial}(3)];

fprintf('Fmincon in %3.2f min\n',telapsed/60);
fprintf('Caso migliore %d\n',iMinTrial);
fprintf('Minimo trovato %3.4f \n ',MinTrial)
fprintf('X_win  = [%3.3f %3.3f %3.3f]\n ',x_win(1),x_win(2),x_win(3));


FULL_OPT = {x,fval,exitflag,output,lambda,grad,hessian,iMinTrial};

WIN_OPT  = {x{iMinTrial},MinTrial,exitflag(iMinTrial),output{iMinTrial},...
           lambda{iMinTrial},grad{iMinTrial},hessian{iMinTrial}};

       
% Course, fitness e ineq
[ fitnessC      ] =  switchFunCostr(x_win,BU_par,xp1,yp1,xp2,yp2,...
    {},{},{},{},{},{},'FitnessValue');

[ inequalityC   ] =  switchFunCostr(x_win,BU_par,xp1,yp1,xp2,yp2,...
    {},{},{},{},{},{},'Constain');
% CourseCorrected, fitness e ineq
[ fitnessCopt   ] =  switchFunCostr(x_win,BU_par,xp1,yp1,xp2,yp2,...
    fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,'FitnessValue');

[ inequalityCopt] =  switchFunCostr(x_win,BU_par,xp1,yp1,xp2,yp2,...
    fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,'Constain');

% Valutazine fvals_c e ISsol_c     
[c_cell,usefulCoor] = cost_function(x_win,BU_par,xp1,yp1,xp2,yp2);
% Valutazine fvals_copt e ISsol_copt
[copt_cell] = cost_function(x_win,BU_par,xp1,yp1,xp2,yp2,...
           fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk);

% riordino 
% ...     = {[CL;CD],       cp,           ...
c_cell    = {c_cell{1},    c_cell{2},    fitnessC,    inequalityC};
copt_cell = {copt_cell{1}, copt_cell{2}, fitnessCopt, inequalityCopt};



end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ returnValue,ceq  ] = ...
    switchFunCostr(IN,BU_par,xp1,yp1,xp2,yp2,...
    fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk,requestedCase)

% evito di lanciare due volta la stessa funzione per vincolo nonlin e f
ceq = [];

persistent fitnessValue constr oldInputs interCoord;

if isempty(oldInputs)
    oldInputs = [];     
end

% isAllreadyCalculated? vediamo
if size(oldInputs,1) == 0
    isAllreadyCalculated = 0;
else
    % ismember(A,B,'rows') for matrices A and B with the same number
    % of columns, returns a vector containing true where the rows of A are
    % also rows of B and false otherwise.

    [isAllreadyCalculated, iOldInput] =...
        ismember(IN, oldInputs, 'rows');
end

% isAllreadyCalculated = 1 ho già valutato l'individuo e l'ho salvato
% in OldInput nella posizione iOldInput
% isAllreadyCalculated = 0 non ho ancora valutato l'individuo

if isempty(oldInputs) || ~isAllreadyCalculated
    % calcolo, non correggo ancora!
    [f,usefulCoor] = cost_function(IN,BU_par,xp1,yp1,xp2,yp2);
    
    xc = usefulCoor{1};
    % salvo in colonne
    fitnessValue(:,end+1)= f{1};
    constr(:,end+1)      = f{2};
    interCoord{end+1}    = xc;

    % oldInputs è di dimensione nx3
    oldInputs = [oldInputs;IN]; %add current inputs to persistent array
    
    
    %returnValue_cell = {fitnessValue, constr};
    returnValue_cell = {f{1},f{2}};  

else
    % vado a recuperare l'elemento
     returnValue_cell = {fitnessValue(:,iOldInput), constr(:,iOldInput)}; % return old values
     xc = interCoord{iOldInput};
end


%% SELEZIONE CASO E CORREZIONE
if strcmp(requestedCase, 'FitnessValue')

     % già colonna
     COEFF = returnValue_cell{1};
    
   
    if size(fvals_f,2) >= 1

        COEFF = fvals_f{end} + Skf{end}*(COEFF-fvals_c{end});
        
    end
    
    % returnValue è la funzione da minimizzare
    %returnValue = -COEFF(1)/COEFF(2);
    returnValue = -COEFF(1);
    
    
elseif strcmp(requestedCase, 'Constain')
    
    % già colonna
    cp = returnValue_cell{2};
    
    % just for safety
    if size(cp,2) > size(cp,1)
        cp = cp';
    end
    
    
    if size(fvals_f,2) >= 1
        % SB correction
        if size(cp) ~= size(ISsol_c{end})
            if max(size(cp)) ~= max(size(Skk{end}))
                needDebug = 1;
            end
        end
        
        cp = ISsol_f{end} + Skk{end}*(cp-ISsol_c{end});
        
    end
    
    % returnValue è il vettore <0
    n = (size(xc,2))/2;
    
    cp_air  = cp(1:(2*n));
    cp_slat = cp(2*n+1:end);
    
    D_me = min(cp_air - cp_air(1));   % deltaCp mainElement
    D_se = min(cp_slat - cp_slat(1)); % deltaCp slatElement
    
    
    returnValue = [-14-D_me ; -14-D_se]; 
    
    
    
else
    error('MyApp:SelectOutput','Requested case not supported')
end

end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,usefulCoor] = cost_function(IN,BU_par,xp1,yp1,xp2,yp2,...
           fvals_f,fvals_c,Skf,ISsol_f,ISsol_c,Skk) 


if nargin == 6
   
   fvals_f = {}; ISsol_f = {};
   fvals_c = {}; ISsol_c = {};
   Skf     = {}; Skk     = {};  
      
end

% possibili output
% [xc,dl,theta_G,cp,xp,yp] = HSA_core_light(IN,BU_par,xp1,yp1,xp2,yp2,PLT,k,save_coor)  
  [xc,dl,theta_G,cp,xp,yp] = HSA_core_light(IN,BU_par,[xp1;xp2],[yp1;yp2],0,  0,0);
    
  [CL,CD] = aeroCoeff(xc,dl,theta_G,cp,BU_par );
  
  % output come colonne!
  f{1} = [CL;CD];
  
  % se riga giro come colonna
  if size(cp,1) < size(cp,2)
     cp = cp';
  end
  
  f{2} = cp; % colonna
  

  % sto calcolando le celle in uscita, necessito correzione
  if size(fvals_f,2) >= 1
      
      % corrego
       COEFF = fvals_f{end} + Skf{end}*([CL;CD]-fvals_c{end});
       cp    = ISsol_f{end} + Skk{end}*(cp-ISsol_c{end});
       
      
       if size(COEFF,2) < size(COEFF,1)
           COEFF = COEFF';
       end
       
              
       f{1} = COEFF;
       f{2} = cp; 
       
  end

  usefulCoor = {xc,xp,yp};
  
%  n = (size(xc,2))/2;
    
%   cp_air  = cp(1:(2*n));
%   cp_slat = cp(2*n+1:end);
%   
%   D_me = min(cp_air - cp_air(1));
%   D_se = min(cp_slat - cp_slat(1));  
%   f{1} = norm( [D_me; D_se] - [-14;-14] );    % f
%   f{2} = abs(  [D_me; D_se])- [14;14];        % cineq 
%   
%   if dumpInfo == 1
%       f = cp;
%       fprintf('\nMain element DCp = %f \n\n\n',D_me,D_se);
%                
%   end
  
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CL,CD] = aeroCoeff( xc,dl,theta_G,cp,BU_par )

n = size(xc,2)/2;

cp_air  = cp(1:2*n);
cp_slat = cp(2*n+1:end);
% 
% punti_ventre_air = [xc(1,1:n)',    yc(1,1:n)'];
% punti_dorso_air =  [xc(1,n+1:end)',yc(1,n+1:end)'];
% 
% punti_ventre_slat = [xc(2,1:n)',    yc(2,1:n)'];
% punti_dorso_slat  = [xc(2,n+1:end)',yc(2,n+1:end)'];
% 
% x_p_ventre_air = [xc(1,1:n)',cp(1:n)];
% x_p_dorso_air  = [xc(1,n+1:end)',cp(n+1:2*n)];
% 
% x_p_ventre_slat = [xc(2,1:n)',cp(2*n+1:3*n)];
% x_p_dorso_slat  = [xc(2,n+1:end)',cp(3*n+1:end)];


CD_ref_air = sum(  dl(1,:).*cp_air'.*sin(theta_G(1,:)) + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
CL_ref_air = sum( -dl(1,:).*cp_air'.*cos(theta_G(1,:)) - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));

CL = CL_ref_air*cosd(BU_par.alpha) - CD_ref_air*sind(BU_par.alpha);
CD = CL_ref_air*sind(BU_par.alpha) + CD_ref_air*cosd(BU_par.alpha);

CD = CD/1000;
CL = CL/1000;
end
