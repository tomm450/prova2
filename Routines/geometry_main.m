function [ GEOM,fail,log ] = geometry_main(Parameters,IN)

% VERSIONE CHE NON PREVEDE INVILUPPO,PLOT o SCRITTURE DI LOG, versione usata per ottimizzazione

%% Descrizione (DA AGGIORNARE)

% input
    
% output


%% FLAG PLOT DI DEBUG
toll = Parameters.Geom.toll_geometry; % precisione bisezione
plot_deb = 0;
log = zeros(1,5); 

%% INIZIO
% TRADUCO INPUT
%Vettore variabile inizio slat [up, dwn]
slat_start = [IN(1),IN(3)];
%Vettore variabile corda [up, dwn]
slat_x     = [IN(2),IN(4)];

% scalare estensione slat superiore
EXT = IN(5);
% coordinata y nuovo LE
y_le = IN(6);
% matrice utile a costrire curve interne up
COM_POINTS_u = [0.5 IN(7)];
% matrice utile a costrire curve interne dwn
COM_POINTS_l = [0.5 IN(8)];

% SOLUZIONI FARLOCCHE nel caso geometria non fosse realizzabile
fail = 1;
GEOM.l_u = 0; GEOM.l_l = 0;
GEOM.slat_u_ret = [ 0 0 ];
GEOM.slat_l_ret = [ 0 0 ];
GEOM.slat_u_dep = [ 0 0 ];
GEOM.contact_u_index = 0;
GEOM.slat_l_dep = [ 0 0 ];
GEOM.pivot_x     = 0;
GEOM.contact_l_index(1) = 0;
GEOM.contact_l_x(1) = 0;

GEOM.up_land = [0 0];
GEOM.dwn_land = [0 0];

GEOM.slat_u_up  = [0 0];
GEOM.slat_u_dwn = [0 0];

GEOM.dwn_shell = [ 0 0 ];
GEOM.up_shell  = [ 0 0 ];

GEOM.BODY = [GEOM.up_shell; flipud(GEOM.dwn_shell(1:end-1,:))]; 
GEOM.division = 0;
%%
% Calcolo effettivo
up  = Parameters.Airfoil.up;
dwn = Parameters.Airfoil.dwn;

seed_special = 5;
thick_toll = mean(diff(up(:,1)))/(seed_special+2);

seed = 0;

% AGGIUNGO PUNTI UTILI PARTE UP
[up1] = thickener(up,slat_start(1),thick_toll,seed_special);
[up2] = thickener(up1,(slat_start(1)+slat_x(1)),thick_toll,seed);

up = up2;
clear up1 up2 up3
% AGGIUNGO PUNTI UTILI PARTE DWN

[dwn1] = thickener(dwn,slat_start(2),thick_toll,seed_special);
[dwn2] = thickener(dwn1,(slat_start(2)+slat_x(2)),thick_toll,seed);

dwn = dwn2;
clear dwn1 dwn2

if plot_deb == 1
    figure(99)
    plot(up(:,1),up(:,2),'bx--'); 
    hold on
    plot(dwn(:,1),dwn(:,2),'cx--')
    grid on
    legend('Dorso','Ventre')
    title('Profilo Base')
end

%isolo punti 
[~,i_start_u] = min(abs(up(:,1) - slat_start(1)));
[~,i_end_u  ] = min(abs(up(:,1) - (slat_start(1)+slat_x(1))));

[~,i_start_l] = min(abs(dwn(:,1) - slat_start(2)));
[~,i_end_l  ] = min(abs(dwn(:,1) - (slat_start(2)+slat_x(2))));

slat_u_ret = up(i_start_u:i_end_u,:);    % slat retratti, curve esterna
slat_l_ret = dwn(i_start_l:i_end_l,:);

% angolo corda retratta
a_u = atan2( (slat_u_ret(end,2) - slat_u_ret(1,2)),...
    ( slat_u_ret(end,1) - slat_u_ret(1,1)));
%a_l = 2*pi - atan2( (slat_l_ret(end,2) - slat_l_ret(1,2)),...
%    ( slat_l_ret(end,1) - slat_l_ret(1,1)));


% corda slats
l_u = sqrt( (slat_u_ret(end,2) - slat_u_ret(1,2))^2 ...
    + (slat_u_ret(end,1) - slat_u_ret(1,1))^2);
GEOM.l_u = l_u;

l_l = sqrt( (slat_l_ret(end,2) - slat_l_ret(1,2))^2 ...
    + (slat_l_ret(end,1) - slat_l_ret(1,1))^2);
GEOM.l_l = l_l;


log(5) = 0;
%% INVOCO SPESSORE (1)
% curva superiore
up_old = up;   % up e dwn andranno modificato e considerati con slat estesi 
dwn_old = dwn;

[up, win,slat_u_ret,GEOM.th_up] =  curve_generator(up,COM_POINTS_u, [i_start_u i_end_u] ,1);
GEOM.slat_u_ret = slat_u_ret;
% 
 if win == 0
     GEOM.th_dwn = 0;
     return
 end
log(5) = 0.3;

[up] = thickener(up,EXT,thick_toll,seed);
[~,j] = min(abs(EXT-up(:,1)));
GEOM.contact_u_index = j;

% curva inferiore
[dwn,win,slat_l_ret,GEOM.th_dwn] = curve_generator(dwn,COM_POINTS_l,[i_start_l i_end_l],2);
GEOM.slat_l_ret = slat_l_ret;

if win == 0
    return
end
log(5) = 0.6;

if plot_deb == 1
    figure(100)
    plot(slat_u_ret(:,1),slat_u_ret(:,2),'b','LineWidth',2)
    title('Dorso')
    hold on
    plot(up(:,1),up(:,2),'r')
    grid on
    %axis equal
    
    figure(105)
    plot(slat_u_ret(:,1),slat_u_ret(:,2),'b-','LineWidth',2)
    title('Slat inferiore')
    grid on
    
    figure(200)
    plot(slat_l_ret(:,1),slat_l_ret(:,2),'b','LineWidth',2)
    title('Ventre')
    hold on
    plot(dwn(:,1),dwn(:,2),'r-')
    %axis equal
    grid on
    
    figure(205)
    plot(slat_l_ret(:,1),slat_l_ret(:,2),'b','LineWidth',2)
    title('Slat inferiore')
    grid on
end

log(5) = 1;
%%  EVENTUALI VINCOLI COSTRUTTIVI (2)
[win] = vincoli(slat_l_ret);
if win == 0
    return
end
[win] = vincoli(slat_u_ret);
if win == 0
    return
end

log(5) = 2;
%% %% FASE 1 = RUOTO SLAT SUPERIORE (3)
%slat pre_rotazioni ( il pivot della rotazione deve coincidere con origine
%                     per usare matrice R )

slat_u_prerot_x = slat_u_ret(:,1) - slat_u_ret(end,1);
slat_u_prerot_y = slat_u_ret(:,2) - slat_u_ret(end,2);
slat_u_temp = [slat_u_prerot_x,slat_u_prerot_y];

% ruoto up ( pivot su TE slat superiore)
pivot = up(j,:);

GEOM.pivot_x(1)     = pivot(1);


if pivot(:,2) < 0
    return
end

if l_u<pivot(2)
    return
end

                    
[slat_u_dep,~,iter_bu,err_bu] = Bisez_up(l_u,a_u,pivot,slat_u_temp,toll,y_le);

log(1:2) = [iter_bu,err_bu];

GEOM.slat_u_dep = slat_u_dep;

if plot_deb ==1
    figure(100)
    hold on
    plot(slat_u_dep(:,1),slat_u_dep(:,2),'g')
end

log(5) = 3;
%% FASE 2 = VERIFICO VINCONOLO COMPENETRAZIONE up (4)
%disp('fase 2')
% if plot_deb ==1
%     slat_u_dep_lc =  slat_u_dep(2:2:end,:);
%     figure(100)
%     hold on
%     plot(slat_u_dep_lc(:,1),slat_u_dep_lc(:,2),'k--')
% end

[win] = comp_verifier(up,slat_u_dep,1);
if win == 0
    %%disp('[win] = comp_verifier(up,slat_u_dep,1);')
    return
end

log(5) = 4;
%% FASE 3 = SPOSTO RIGIDAMENTE SLAT INFERIORE SU NUOVO LE_profilo (5)
%%disp('fase 3')
slat_l_temp_x = slat_l_ret(:,1) - (slat_l_ret(1,1) - slat_u_dep(1,1));
slat_l_temp_y = slat_l_ret(:,2) - (slat_l_ret(1,2) - slat_u_dep(1,2));
slat_l_temp = [ slat_l_temp_x, slat_l_temp_y];

% TE slat_l_temp
TE_l = slat_l_temp(end,:);

% pivot rotazione slat inferiore
PX = slat_u_dep(1,1); %PY = slat_u_dep(1,2);

log(5) = 5;
%% FASE 4 = VERIFICO CHE IL PROFILO SIA RAGGIUNGIBILE (6)
%%disp('fase 4')
if 0.9*abs(l_l) <= abs(PX) % assumo che LE_airfoil_or sia a 0,0
    % punto SICURAMENTE non raggiungibile
    log(5) = 6;
    return
else
    %% FASE 5.1 = PROFILO RAGGIUNGIBILE,MA NECESSITO DI PREROTAZIONE (7)
    %%disp('fase 51')
    a1 = asin(abs(TE_l(2))/l_l);
    a2 = acos(abs(PX(1))/l_l);
    
    a_low = 0.5*(a1-a2); a_up = 1.5*a1; a_mid = 0.5*(a_up+a_low);
    %else
    %alpha_s = 0;
    %end
    
    %% FASE 5.2 = PROFILO RAGGIUNGIBILE, RUOTO 
    %cerco indice intersezione
    alpha_dwn = [ a_low a_mid a_up];
    
    [slat_l_dep,~,iter_bd,err_bd] = Bisez_dwn(dwn,slat_l_temp,alpha_dwn,toll);
    
    log(3:4) = [iter_bd,err_bd];
    
    GEOM.slat_l_dep = slat_l_dep;
    
    [dwn_body] = thickener(dwn,slat_l_dep(end,1),thick_toll,seed);
    
    if abs(err_bd) > toll
        error('Bisez_dwn tolleranza non raggiunta fase 5.2')    
    end
    
    if plot_deb ==1
        figure(200)
        hold on
        plot(slat_l_dep(:,1),slat_l_dep(:,2),'g')
        
%         figure(300)
%         plot(slat_l_dep(:,1),slat_l_dep(:,2),'g')
%         hold on
%         plot(slat_u_dep(:,1),slat_u_dep(:,2),'r')
%         plot(up(:,1),up(:,2),'b')
%         plot(dwn(:,1),dwn(:,2),'b')
%         axis equal
    end
    log(5) = 7;
    %% FASE 55 = VERIFICO VINCONOLO COMPENETRAZIONE slat-profilo (8)
%     if plot_deb ==1
%         slat_l_dep_uc =  slat_l_dep(2:2:end,:);
%         figure(101)
%         hold on
%         plot(slat_l_dep_uc(:,1),slat_l_decp_uc(:,2),'k--')
%     end
    %disp('fase 55')
    [win] = comp_verifier(dwn_body,slat_l_dep,2);
    if win == 0
        %%disp('[win] = comp_verifier(dwn_body,slat_l_dep,2)')
        return
    end
    log(5) = 8;
    %% FASE 56 = VERIFICO VINCONOLO COMPENETRAZIONE slat-slat (9)
    %disp('fase 56')
    [win] = comp_verifier_slat_slat(slat_u_dep,slat_l_dep);
    if win == 0
        %%disp('[win] = comp_verifier_slat_slat(slat_u_dep,slat_l_dep);')
        return
    end
    log(5) = 9;
end

%% FASE 6 = CALCOLO INDICE RACCORDO SLAT-VENTRE (10)
% infittisco vettore dwn per calcolare in maniera pi� accurata derivare e
% curve
% ho gi� infittito
[~,j_low]  = min(abs( slat_l_dep(end,1) - dwn_body(:,1)));

GEOM.contact_l_index = j_low;
GEOM.contact_l_x(1) = dwn_body(j_low);

if j_low ==1
    return
end
log(5) = 10;
%% SONO ARRIVATO QUI GLI SPESSORI NON CREANO PROBLEMI
%ricavo profilo esterno slat
slat_l_dep_ext = slat_l_dep(1:2:end,:);
slat_l_dep_ext = [slat_l_dep_ext; slat_l_dep(end,:)];
slat_u_dep_ext = slat_u_dep(1:2:end,:);
slat_u_dep_ext = [slat_u_dep_ext; slat_u_dep(end,:)];

% CALCOLO GEOMETRIA SLAT SUPERIORE
s = slat_u_ret;
a = atan((s(end,2)-s(1,2))/(s(end,1)-s(1,1)));

s(:,1) = s(:,1) - s(end,1);
s(:,2) = s(:,2) - s(end,2);

R = [cos(a) sin(a);-sin(a) cos(a)]; 
c = R*s';
c = c';

cd = c(2:2:end,:); 
cd = [c(1,:); cd];

if cd(end,:) == s(end,:)
    % ho gia preso tutti i punti
else
    cd = [cd;c(end,:)];
end

cu = c(1:2:end,:);
if cu(end,:) == c(end,:)
    % ho gia preso tutti i punti
else
    cu = [cu;c(end,:)];
end

if cu(1,:) == cd(1,:) 
    %ok
else
    cu(1,:)
    cd(1,:)
    error('LE slat non coincide')
end

if cu(end,:) == cd(end,:) 
    %ok
else
    error('TE slat non coincide')
end


% OUTPUT GEOMETRICI
GEOM.dwn_shell = [slat_l_dep_ext;dwn_body(j_low+1:end,:)];
GEOM.up_shell  = [slat_u_dep_ext;up(j+1:end,:)];

GEOM.up_land = up;
GEOM.dwn_land = dwn_old;
GEOM.slat_u_up  = cu;
GEOM.slat_u_dwn = cd;

GEOM.BODY = [GEOM.up_shell; flipud(GEOM.dwn_shell(1:end-1,:))]; 


% caso limite 9 segmenti

division(1) = 1;
% cerco in body gli indici che dividono le varie curve
[~,i_d1] = min(abs(GEOM.up_shell(:,1) - EXT));
division(2) = i_d1; %slat_dep - profilo

[~,i_d2] = min(abs(GEOM.up_shell(:,1) - slat_start(1)));
division(3) = i_d2;  %profilo - slat neg

[~,i_d3] = min(abs(GEOM.up_shell(:,1) - (slat_start(1)+slat_x(1))));
division(4) = i_d3;  %slat neg - profilo

[~,i_d4] = max(GEOM.BODY(:,1));
division(5) = i_d4;

[~,i_d5] = min(abs(GEOM.BODY(:,2) - dwn(i_end_l,2)));
division(6) = i_d5; %slat neg - profilo

[~,i_d6] = min(abs(GEOM.BODY(:,2) - dwn(i_start_l,2)));
division(7) = i_d6; %profilo - slat neg

[~,i_d7] = min(abs(GEOM.BODY(:,2) - dwn_body(j_low,2)));
division(8) = i_d7; %slat_dep - profilo

division(9) = size(GEOM.BODY,1);

% division = GEOM.division;

delay = 0;

for w = 1: size(division,2)-1
    if w > size(division,2)-1-delay
        break
    end
    
    if division(w+1) <= division(w)
        
        division = [division(1:w),division(w+2:end)];
        delay = 1;
    end
    
end

GEOM.division = division;


%%
if plot_deb == 1
    figure(201)
    plot(slat_l_dep_ext(:,1),slat_l_dep_ext(:,2),'bo')
    hold on
    plot(slat_u_dep_ext(:,1),slat_u_dep_ext(:,2),'ro')
    plot(up(:,1),up(:,2),'gx')
    plot(dwn_body(:,1),dwn_body(:,2),'gx')
    grid on
    title('Profilo con slat estesi')
    
    figure(202)
    plot(GEOM.dwn_shell(:,1),GEOM.dwn_shell(:,2),'bo')
    hold on
    plot(GEOM.up_shell(:,1),GEOM.up_shell(:,2),'ro')
    grid on
    title('BODY')
    
end

%% sono arrivato in fondo
fail = 0;

end




