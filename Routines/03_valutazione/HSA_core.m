function [xc,dl,theta_G,cp] = HSA_core(IN,BU_par,Parameters,GEOM,PLT,k,save_coor)
% HSAfo full output
%[cp_me,xp_me,cp_se,xp_se]

if nargin == 4
    PLT = 0; k = 1; save_coor = 0;
elseif nargin == 5
    k = 1; save_coor = 0;
elseif nargin == 5
    save_coor = 0;
end

%Parameters.HSA.npane
G_pane_airfoil = [flipud(GEOM.dwn_land)  ; GEOM.up_land(2:end,:)];
G_pane_slat    = [flipud(GEOM.slat_u_dwn); GEOM.slat_u_up(2:end,:)];

% inizializzo con vettore completo
xp1 = G_pane_airfoil(:,1)';
xp2 = G_pane_slat(:,1)';

[dummy,G_pane_air_LE]  = min(G_pane_airfoil(:,1));
[dummy,G_pane_slat_LE] = min(G_pane_slat(:,1));

index1 = round(linspace(1,G_pane_air_LE,Parameters.HSA.npane+1));
index1 = [index1(1:end-1),round(linspace(G_pane_air_LE,size(xp1,2),Parameters.HSA.npane+1))];

index2 = round(linspace(1,G_pane_slat_LE,Parameters.HSA.npane+1));
index2 = [index2(1:end-1),round(linspace(G_pane_slat_LE,size(xp2,2),Parameters.HSA.npane+1))];

xp1 = xp1(index1); yp1 = G_pane_airfoil(index1,2)';
xp2 = xp2(index2); yp2 = G_pane_slat(index2,2)';


alpha_flap = -IN(3)*pi/180;
alpha = 0;

xx1 = ([xp1;yp1]'*[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)])';

%dx_flap = [IN(1),IN(2)]*[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];

xx2 = ([xp2;yp2]'*[cos(alpha+alpha_flap), sin(alpha+alpha_flap);...
                   -sin(alpha+alpha_flap), cos(alpha+alpha_flap)])';

xp1 = xx1(1,:);
xp2 = xx2(1,:)+IN(1);
yp1 = xx1(2,:);
yp2 = xx2(2,:)+IN(2);

xp = [xp1;xp2]/1000;
yp = [yp1;yp2]/1000;

if save_coor == 1
    save('current_design.mat','xp','yp');
end

U_inf = [BU_par.Ux BU_par.Uz];

[~,A,B,xc,yc,dl,N,T,x] = HS_staz_multi_2(xp,yp,U_inf);

U_inf_t = [U_inf*T(:,:,1), U_inf*T(:,:,2)];
u_pert = B*x + U_inf_t';

%% Coefficiente di pressione
cp = 1-(u_pert./norm(U_inf)).^2;

[theta_G(1,:)] = panels_inclination(xp(1,:),yp(1,:));
[theta_G(2,:)] = panels_inclination(xp(2,:),yp(2,:));

%n = size(xc,2)/2;

%cp_air  = cp(1:max(size(xc)));
%cp_slat = cp(max(size(xc))+1:end);

%punti_ventre_air = [xc(1,1:n)',    yc(1,1:n)'];
%punti_dorso_air =  [xc(1,n+1:end)',yc(1,n+1:end)'];

%punti_ventre_slat = [xc(2,1:n)',    yc(2,1:n)'];
%punti_dorso_slat  = [xc(2,n+1:end)',yc(2,n+1:end)'];

%x_p_ventre_air = [xc(1,1:n)',cp(1:n)];
%x_p_dorso_air  = [xc(1,n+1:end)',cp(n+1:2*n)];

%x_p_ventre_slat = [xc(2,1:n)',cp(2*n+1:3*n)];
%x_p_dorso_slat  = [xc(2,n+1:end)',cp(3*n+1:end)];


%CD_ref_air = sum(  dl(1,:).*cp_air'.*sin(theta_G(1,:)) + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
%CL_ref_air = sum( -dl(1,:).*cp_air'.*cos(theta_G(1,:)) - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));

%CL = CL_ref_air*cosd(BU_par.alpha) - CD_ref_air*sind(BU_par.alpha);
%CD = CL_ref_air*sind(BU_par.alpha) + CD_ref_air*cosd(BU_par.alpha);

% plot qui perchè in caso di debug non devo aspettare
if PLT == 1
    
    n = size(xc,2)/2;
    
    figure(k)
    % profilo e punti notevoli
    subplot(1,2,1)
    axis equal
    hold on
    plot(xp(1,:),yp(1,:),'Color',[0 0.5 1],'LineWidth',2);
    plot(xp(2,:),yp(2,:),'Color',[0 0.5 1],'LineWidth',2);
    title(sprintf('alpha = %2.3f deg; X_{IN} = [ %1.3f %1.3f %1.3f ]',BU_par.alpha,IN(1),IN(2),-IN(3)));
         
    hold on
    plot(xc(1,:),yc(1,:),'xr','LineWidth',2);
    plot(xc(2,:),yc(2,:),'xr','LineWidth',2);

    % distribuzione pressione
    subplot(1,2,2)
    hold on
    x_d_combo  = [xc(1,n+1:2*n),nan,xc(2,n+1:end)];
    cp_d_combo = [cp(n+1:2*n)',nan,cp(3*n+1:end)'];
    plot(x_d_combo,cp_d_combo,'ro-','LineWidth',2);
%     plot(xc(1,n+1:2*n),cp(n+1:2*n),'Color',[1 0.5 0],'LineWidth',2);
%     plot(xc(2,n+1:end),cp(3*n+1:end),'Color',[1 0.5 0],'LineWidth',2);
    %%
    
    x_v_combo  = [xc(1,1:n),nan,xc(2,1:n)];
    cp_v_combo = [cp(1:n)',nan,cp(2*n+1:3*n)'];

    plot(x_v_combo,cp_v_combo,'bo-','LineWidth',2);  

    grid on
    set(gca,'YDir','Reverse'); % reverse y axis
    %alfa=num2str(alpha*180/pi);
    tit1 = 'COEFFICIENTE DI PRESSIONE';
    tit2 = ['Incidenza: ',BU_par.alpha, ' gradi'];
    title({tit1;tit2})
    xlabel('corda');
    ylabel('-C_p');
    
    legend('Dorso HS','Ventre HS');

    %figure(33)
    %for i = 1:n

    %    quiver(xc(1,:),yc(1,:),dl(1,:).*cp_air'.*sin(theta_G(1,:)),-dl(1,:).*cp_air'.*cos(theta_G(1,:)))

    %    hold on;

    %    quiver(xc(2,:),yc(2,:),dl(2,:).*cp_slat'.*sin(theta_G(2,:)),-dl(2,:).*cp_slat'.*cos(theta_G(2,:)))

    %end

end


end

%% SUBROUTINE
function [lfun,A,B,xc,yc,dl,N,T,x] = HS_staz_multi_2(xp,yp,U_inf)
%questo � giusto
lfun = localfunctions;
if nargout == 1
    return
end

Npt = size(xp,2);
Npan = Npt-1;

%% Punti di collocazione
[xc(1,:),yc(1,:)] = collocation_points(xp(1,:),yp(1,:));
[xc(2,:),yc(2,:)] = collocation_points(xp(2,:),yp(2,:));
%% Lunghezze pannelli
[dx(1,:),dy(1,:),dl(1,:)] = panels_dimensions(xp(1,:),yp(1,:));
[dx(2,:),dy(2,:),dl(2,:)] = panels_dimensions(xp(2,:),yp(2,:));
%% Versori normali e tangenti
[N(:,:,1),T(:,:,1)] = NTversors(dx(1,:),dy(1,:),dl(1,:));
[N(:,:,2),T(:,:,2)] = NTversors(dx(2,:),dy(2,:),dl(2,:));

%% Angoli globali
[theta_G(1,:)] = panels_inclination(xp(1,:),yp(1,:));
[theta_G(2,:)] = panels_inclination(xp(2,:),yp(2,:));
%% Velocit� indotte
% Preallocazione matrici coefficienti
A1 = zeros(Npt,Npt);
A2 = A1;
C1 = A1;
C2 = C1;
B1 = zeros(Npt-1,Npt);
B2 = B1;
D1 = B1;
D2 = D1;
b1 = zeros(Npt,1);
b2 = b1;

for j = 1:Npan
    % Matrice di rotazione per portarci nel sistema pannello
    [R1] = rotation_matrix(theta_G(1,j));
    [R2] = rotation_matrix(theta_G(2,j));
    for k = 1:Npan
        % Raggi fra estremi pannello j e punti di collocazione k
        % Angoli fra estremi pannello j e punti di collocazione k
        [r1_1,r2_1,dt_1] = radius_dtheta(xc(1,k),yc(1,k),xp(1,j),...
                                   yp(1,j),xp(1,j+1),yp(1,j+1),R1,j,k);

        [r1_2,r2_2,dt_2] = radius_dtheta(xc(2,k),yc(2,k),xp(2,j),...
                                   yp(2,j),xp(2,j+1),yp(2,j+1),R2,j,k);
        % Velocit� indotte SORGENTI
        [Us1,Vs1] = source_induced_speed(r1_1,r2_1,dt_1);
        [Us2,Vs2] = source_induced_speed(r1_2,r2_2,dt_2);
        % Velocit� indotte VORTICI
        [Uv1,Vv1] = vortex_induced_speed(r1_1,r2_1,dt_1);
        [Uv2,Vv2] = vortex_induced_speed(r1_2,r2_2,dt_2);
        % Vettori velocit� GLOBALI
        S_vel1 = R1'*[Us1 Vs1]';
        S_vel2 = R2'*[Us2 Vs2]';
        V_vel1 = R1'*[Uv1 Vv1]';
        V_vel2 = R2'*[Uv2 Vv2]';

        A1(k,j) = S_vel1'*N(:,k,1);
        C1(k,j) = V_vel1'*N(:,k,1);
        B1(k,j) = S_vel1'*T(:,k,1);
        D1(k,j) = V_vel1'*T(:,k,1);

        A2(k,j) = S_vel2'*N(:,k,2);
        C2(k,j) = V_vel2'*N(:,k,2);
        B2(k,j) = S_vel2'*T(:,k,2);
        D2(k,j) = V_vel2'*T(:,k,2);


    end
    % Vettore termini noti
    b1(j,1) = -U_inf*N(:,j,1);
    b2(j,1) = -U_inf*N(:,j,2);
end

%% Somma velocit� indotte dai vortici su ogni pannello
A1(:,end) = sum(C1,2);
B1(:,end) = sum(D1,2);

A2(:,end) = sum(C2,2);
B2(:,end) = sum(D2,2);

%% Condizione di Kutta
 A1(end,:) = B1(1,:)+ B1(end,:);
 b1(end,1) = -U_inf*(T(:,1,1)+T(:,end,1));

 A2(end,:) = B2(1,:)+ B2(end,:);
 b2(end,1) = -U_inf*(T(:,1,2)+T(:,end,2));

 %matrici d'influenza incrociate
             A12 = zeros(Npt,Npt);
            A21 = A12;
            C12 = A12;
            C21 = C12;
            B12 = zeros(Npt-1,Npt);
            B21 = B12;
            D12 = B12;
            D21 = D12;

            for j = 1:Npan
                % Matrice di rotazione per portarci nel sistema pannello
                [R1] = rotation_matrix(theta_G(1,j));
                [R2] = rotation_matrix(theta_G(2,j));
                for k = 1:Npan
                    % Raggi fra estremi pannello j e punti di collocazione k
                    % Angoli fra estremi pannello j e punti di collocazione k
                    [r1_12,r2_12,dt_12] = radius_dtheta(xc(1,k),yc(1,k),xp(2,j),...
                                               yp(2,j),xp(2,j+1),yp(2,j+1),R2,j,k);

                    [r1_21,r2_21,dt_21] = radius_dtheta(xc(2,k),yc(2,k),xp(1,j),...
                                               yp(1,j),xp(1,j+1),yp(1,j+1),R1,j,k);
                    % Velocit� indotte SORGENTI
                    [Us12,Vs12] = source_induced_speed(r1_12,r2_12,dt_12);
                    [Us21,Vs21] = source_induced_speed(r1_21,r2_21,dt_21);
                    % Velocit� indotte VORTICI
                    [Uv12,Vv12] = vortex_induced_speed(r1_12,r2_12,dt_12);
                    [Uv21,Vv21] = vortex_induced_speed(r1_21,r2_21,dt_21);
                    % Vettori velocit� GLOBALI
                    S_vel12 = R1'*[Us12 Vs12]';
                    S_vel21 = R2'*[Us21 Vs21]';
                    V_vel12 = R1'*[Uv12 Vv12]';
                    V_vel21 = R2'*[Uv21 Vv21]';

                    A12(k,j) = S_vel12'*N(:,k,1);
                    C12(k,j) = V_vel12'*N(:,k,1);
                    B12(k,j) = S_vel12'*T(:,k,1);
                    D12(k,j) = V_vel12'*T(:,k,1);

                    A21(k,j) = S_vel21'*N(:,k,2);
                    C21(k,j) = V_vel21'*N(:,k,2);
                    B21(k,j) = S_vel21'*T(:,k,2);
                    D21(k,j) = V_vel21'*T(:,k,2);


                end
            end

            %% Somma velocit� indotte dai vortici su ogni pannello
            A12(:,end) = sum(C12,2);
            B12(:,end) = sum(D12,2);

            A21(:,end) = sum(C21,2);
            B21(:,end) = sum(D21,2);

            %% Condizione di Kutta
             A12(end,:) = B12(1,:)+ B12(end,:);

             A21(end,:) = B21(1,:)+ B21(end,:);




 % ASSEMBLO LE MATRICI D'INFLUENZA
 A = [A1,A12;A21,A2];
 B = [B1,B12;B21,B2];
 b = [b1;b2];

 %% Risoluzione sistema lineare
 if nargout == 8
    return
 end
 x = A\b;

end

function [xc,yc] = collocation_points(xp,yp)
    xc = 0.5*(xp(2:end)+xp(1:end-1));
    yc = 0.5*(yp(2:end)+yp(1:end-1));
end

function [dx,dy,dl] = panels_dimensions(xp,yp)
    dx = xp(2:end)-xp(1:end-1);
    dy = yp(2:end)-yp(1:end-1);
    dl = sqrt(dx.^2+dy.^2);
end

function [N,T] = NTversors(dx,dy,dl)
    N = [(-dy./dl);
         (dx./dl)];
    N=N./sqrt(N(1,:).^2+N(2,:).^2);
    T = [(dx./dl);
         (dy./dl)];
    T=T./sqrt(T(1,:).^2+T(2,:).^2);
end

function [theta_G] = panels_inclination(xp,yp)
    theta_G = atan2(yp(2:end)-yp(1:end-1),xp(2:end)-xp(1:end-1));
end

function [R] = rotation_matrix(theta)
    R = [cos(theta)  sin(theta)
        -sin(theta) cos(theta)];
end

function [r1,r2,dt] = radius_dtheta(xc,yc,xp,yp,xp2,yp2,R,varargin)
    if nargin == 7
        j=2;k=3;
    else
        j = varargin{1};
        k = varargin{2};
    end

    r1_G = [xc-xp yc-yp];
    r2_G = [xc-xp2 yc-yp2];
    r1 = R*r1_G';
    r2 = R*r2_G';

    t1 = atan2(r1(2),r1(1));
    t2 = atan2(r2(2),r2(1));
    dt = t2-t1;
    if j == k && t2-t1<1e-15
        dt=t1-t2;
    end

end

function [Us,Vs] = source_induced_speed(r1,r2,dt,varargin)
    if nargin == 3
        sigma = 1;
    else
        sigma = varargin{1};
    end

    Us = -sigma./(2*pi)*log(norm(r2)/norm(r1));
    Vs = sigma.*(dt)/(2*pi);
end

function [Uv,Vv] = vortex_induced_speed(r1,r2,dt,varargin)
    if nargin == 3
        gamma = 1;
    else
        gamma = varargin{1};
    end

    Uv = -gamma.*(dt)/(2*pi);
    Vv = -gamma./(2*pi)*log(norm(r2)/norm(r1));
end
