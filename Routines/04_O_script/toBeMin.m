function [x,fval,exitflag,output,lambda,grad,hessian,cp_copt] = toBeMin(OPT,BU_par,Parameters,GEOM,...
              fvals_f,fvals_c,Skf)
% This is an auto generated MATLAB file from Optimization Tool.
nvars = OPT.nvars;
lb    = OPT.lb;
ub    = OPT.ub;
x0    = OPT.x0;

Aineq = OPT.Aineq;
bineq = OPT.bineq;

MaxFunctionEvaluations_Data = OPT.MaxFunctionEvaluations_Data;
MaxIterations_Data          = OPT.MaxIterations_Data;
OptimalityTolerance_Data    = OPT.OptimalityTolerance_Data;
StepTolerance_Data          = OPT.StepTolerance_Data;

%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunctionEvaluations', MaxFunctionEvaluations_Data);
options = optimoptions(options,'MaxIterations', MaxIterations_Data);
options = optimoptions(options,'OptimalityTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'FunctionTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'StepTolerance', StepTolerance_Data);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize });


%% Una volta sola
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


% definisco ingresso fmincon
actual_cost_function = @(x) ...
    switchFunCostr(x,BU_par,xp1,yp1,xp2,yp2,fvals_f,fvals_c,Skf,'FitnessValue');
nonlinear_constr     = @(x) ...
    switchFunCostr(x,BU_par,xp1,yp1,xp2,yp2,fvals_f,fvals_c,Skf,'Constain'); 

tstart = tic;

[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(actual_cost_function,x0,Aineq,bineq,[],[],lb,ub,nonlinear_constr,options);

telapsed = toc(tstart);
fprintf('Fmincon in %3.2F min \n',telapsed/60);

[cp_copt] = cost_function(x,BU_par,xp1,yp1,xp2,yp2,fvals_f,fvals_c,Skf,1);

end

function [ returnValue,ceq  ] = ...
    switchFunCostr(IN,BU_par,xp1,yp1,xp2,yp2,...
    fvals_f,fvals_c,Skf,requestedCase)

% evito di lanciare due volta la stessa funzione per vincolo nonlin e f

ceq = [];

persistent fitnessValue constr oldInputs;

if isempty(oldInputs)
    oldInputs = []; %initialize oldInputs as cell
end

if size(oldInputs,1) == 0
    isAllreadyCalculated = 0;
else
    [isAllreadyCalculated, iOldInput] =...
        ismember(IN, oldInputs, 'rows');
end

if isempty(oldInputs) || ~isAllreadyCalculated
    
    [f] = cost_function(IN,BU_par,xp1,yp1,xp2,yp2,fvals_f,fvals_c,Skf,0);
    
    fitnessValue(end+1)= f{1};
    constr(:,end+1)    = f{2};
    
    oldInputs = [oldInputs;IN]; %add current inputs to persistent array
    returnValue_cell = {fitnessValue(end), constr(:,end)};
    
else
    returnValue_cell = {fitnessValue(iOldInput), constr(:,iOldInput)}; % return old values
end

if strcmp(requestedCase, 'FitnessValue')
    
    returnValue = returnValue_cell{1};
    
elseif strcmp(requestedCase, 'Constain')
    
    c = returnValue_cell{2};
    returnValue = c;
    
else
    error('MyApp:SelectOutput','Requested case not supported')
end

end 

function [f] = cost_function(IN,BU_par,xp1,yp1,xp2,yp2,fvals_f,fvals_c,Skf,dumpInfo)
  
  [xc,~,~,cp] = HSA_core_light(IN,BU_par,xp1,yp1,xp2,yp2,0,0,0);
  
  if size(fvals_f,2) >= 1
      % SB correction      
      cp = fvals_f{end} + Skf{end}*(cp-fvals_c{end});
      
  end
  
  n = (size(xc,2))/2;
    
  cp_air  = cp(1:(2*n));
  cp_slat = cp(2*n+1:end);
  
  D_me = min(cp_air - cp_air(1));
  D_se = min(cp_slat - cp_slat(1));

  f{1} = norm( [D_me; D_se] - [-14;-14] );    % f
  f{2} = abs(  [D_me; D_se])- [14;14];        % cineq 
  
  if dumpInfo == 1
      f = cp;
      fprintf('Main element DCp = %f \nSlat element DCp = %f \n\n',D_me,D_se);
  end
  
%   if dumpInfo == 2
%       f = [D_me; D_se];
%   end
  
  
end
