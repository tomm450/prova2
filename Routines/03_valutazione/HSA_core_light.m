function [xc,dl,theta_G,cp,xp,yp,AM,BM,xM,bM] = HSA_core_light(IN,BU_par,xp,yp,PLT,k,save_coor)
% HSAfo full output
%[cp_me,xp_me,cp_se,xp_se]
        
if nargin == 4
    PLT = 0; k = 1; save_coor = 0;
elseif nargin == 5
    k = 1; save_coor = 0;
elseif nargin == 5
    save_coor = 0;
end


if size(IN,2) == 0 % coordinate assolute
    %xp = [xp1;xp2];%/1000;
    %yp = [yp1;yp2];%/1000;
else
    alpha_flap = -IN(3)*pi/180;
    alpha = 0;
    xp1 = xp(1,:); 
    yp1 = yp(1,:);
    xp2 = xp(2,:); 
    yp2 = yp(2,:);
    
    xx1 = ([xp1;yp1]'*[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)])';
    
    %dx_flap = [IN(1),IN(2)]*[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
    
    xx2 = ([xp2;yp2]'*[cos(alpha+alpha_flap), sin(alpha+alpha_flap);...
        -sin(alpha+alpha_flap), cos(alpha+alpha_flap)])';
    
    xp1 = xx1(1,:);
    xp2 = xx2(1,:)+IN(1);
    yp1 = xx1(2,:);
    yp2 = xx2(2,:)+IN(2);
    xp = [xp1;xp2];%/1000;
    yp = [yp1;yp2];%/1000;
end


if save_coor == 1    
    save(['./Output/',num2str(k),'current_design.mat'],'xp','yp');
end

U_inf = [BU_par.Ux BU_par.Uz];

[~,AM,BM,xc,yc,dl,N,T,xM,bM] = HS_staz_multi_2(xp,yp,U_inf);


U_inf_t = [];

for i = 1:size(xp,1)

    U_inf_t = [U_inf_t, U_inf*T(:,:,i)];
    
    [theta_G(i,:)] = panels_inclination(xp(i,:),yp(i,:));
    
end
   



u_pert = BM*xM + U_inf_t';

%% Coefficiente di pressione
cp = 1-(u_pert./norm(U_inf)).^2;

%[theta_G(2,:)] = panels_inclination(xp(2,:),yp(2,:));

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
    subplot(2,1,1)
    axis equal
    hold on
    for i = 1:size(xp,1)
        plot(xp(i,:),yp(i,:),'Color',[0 0.5 1],'LineWidth',2);
        %plot(xp(2,:),yp(2,:),'Color',[0 0.5 1],'LineWidth',2);
        hold on
        plot(xc(i,:),yc(i,:),'xr','LineWidth',2);
        %plot(xc(2,:),yc(2,:),'xr','LineWidth',2);
    end
    
    st0 = sprintf('Alpha = %d deg;',BU_par.alpha);%
    
    if size(IN,2) >= 3
        st1 = sprintf(' X_{IN} = [ %1.3f %1.3f %1.3f ]',IN(1),IN(2),-IN(3));
    else
        st1 = '';
    end
    
    title(strcat(st0,st1));
    
    subplot(2,1,2)
    hold on
    
    cp_mat = [];
    
    for i = 1:size(xp,1)
        cp_mat = [cp_mat; cp((1+2*(i-1)*n):(2*i*n))'];
        plot(xc(i,n+1:2*n),cp_mat(i,(n+1):2*n),'ro-','LineWidth',2);
        plot(xc(i,1:n)    ,cp_mat(i,1:n),'bo-','LineWidth',2);
    end
          
    % distribuzione pressione

        
%     x_d_combo  = [xc(1,n+1:2*n),nan,xc(2,n+1:end)];
%     cp_d_combo = [cp(n+1:2*n)', nan,cp(3*n+1:end)'];
%     
%     plot(x_d_combo,cp_d_combo,'ro-','LineWidth',2);
%     
%     x_v_combo  = [xc(1,1:n),nan,xc(2,1:n)];
%     cp_v_combo = [cp(1:n)',nan,cp(2*n+1:3*n)'];
% 
%     plot(x_v_combo,cp_v_combo,'bo-','LineWidth',2);  

    grid on
    set(gca,'YDir','Reverse'); % reverse y axis
    %alfa=num2str(alpha*180/pi);
%    tit1 = 'COEFFICIENTE DI PRESSIONE';
%     tit2 = ['Incidenza: ',BU_par.alpha, ' gradi'];
%     title({tit1;tit2})
    xlabel('X coordinate');
    ylabel('-C_p');
    title('Pressure coefficient')
    %legend('Dorso HS','Ventre HS','Location','best');



end


end

%% SUBROUTINE
function [lfun,AM,BM,xc,yc,dl,N,T,x,bM] = HS_staz_multi_2(xp,yp,U_inf)
%questo � giusto

RK = 1; % 1 in teoria
        % 2 ad cazzum ( e tra l'altro non esce)
lfun = localfunctions;

if nargout == 1
    return
end

Npt = size(xp,2);
Npan = Npt-1;

for i = 1:size(xp,1)
    %% Punti di collocazione
    [xc(i,:),yc(i,:)] = collocation_points(xp(i,:),yp(i,:));
    
    %% Lunghezze pannelli
    [dx(i,:),dy(i,:),dl(i,:)] = panels_dimensions(xp(i,:),yp(i,:));
   
    %% Versori normali e tangenti
    [N(:,:,i),T(:,:,i)] = NTversors(dx(i,:),dy(i,:),dl(i,:));
       
    %% Angoli globali
    [theta_G(i,:)] = panels_inclination(xp(i,:),yp(i,:));

end

for p = 1:size(xc,1)
    for q = 1:size(xc,1)
        
            A{p,q} = zeros(Npt,Npt);
            C{p,q} = A{p,q};
            B{p,q} = zeros(Npt-1,Npt);
            D{p,q} = B{p,q};
            

            for j = 1:Npan
                % Matrice di rotazione per portarci nel sistema pannello
                [R1] = rotation_matrix(theta_G(q,j));
                
                for k = 1:Npan
                    % Raggi fra estremi pannello j e punti di collocazione k
                    % Angoli fra estremi pannello j e punti di collocazione k
                    
                    [r1{p,q},r2{p,q},dt{p,q}] = radius_dtheta(xc(p,k),yc(p,k),...
                        xp(q,j),yp(q,j),xp(q,j+1),yp(q,j+1),R1,j,k);

                    % Velocit� indotte SORGENTI
                    [Us{p,q},Vs{p,q}] = source_induced_speed(r1{p,q},r2{p,q},dt{p,q});
                    % Velocit� indotte VORTICI
                    [Uv{p,q},Vv{p,q}] = vortex_induced_speed(r1{p,q},r2{p,q},dt{p,q});
                    % Vettori velocit� GLOBALI
                    S_vel{p,q} = R1'*[Us{p,q} Vs{p,q}]';
                    V_vel{p,q} = R1'*[Uv{p,q} Vv{p,q}]';
                    
                    A{p,q}(k,j) = S_vel{p,q}'*N(:,k,p);
                    C{p,q}(k,j) = V_vel{p,q}'*N(:,k,p);
                    B{p,q}(k,j) = S_vel{p,q}'*T(:,k,p);
                    D{p,q}(k,j) = V_vel{p,q}'*T(:,k,p);


                    
                end
                
                
                if p == q
                    b{p}(j,1) = -U_inf*N(:,j,p);
                end
            end

            %% Somma velocit� indotte dai vortici su ogni pannello
            A{p,q}(:,end) = sum(C{p,q},2);
            B{p,q}(:,end) = sum(D{p,q},2);

            %% Condizione di Kutta
            A{p,q}(end,:) = B{p,q}(1,:)+B{p,q}(end,:);
            
            if p == q
                b{p}(end+1) = -U_inf*(T(:,1,p)+T(:,end,p));
            end

    end
end



 % ASSEMBLO LE MATRICI D'INFLUENZA
 AM = [];
 BM = [];
 bM = [];
 
 for i = 1:size(xc,1)
     Ar = [];
     Br = [];
     for j = 1:size(xc,1)
         Ar = [Ar,A{i,j}];
         Br = [Br,B{i,j}];     
     end
          
              
     bM = [bM;b{i}];  
              
     AM = [AM;Ar];
     BM = [BM;Br];     %
 end
 
 %% Risoluzione sistema lineare
 if nargout == 8
    return
 end

 x = AM\bM;

%  [lfun,AM,BM,xc,yc,dl,N,T,x] 
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
    
     if norm(r1) < 1e-5 || norm(r1) < 1e-5
         Uv = 0; Vv = 0;
     else
        Uv = -gamma.*(dt)/(2*pi);
        Vv = -gamma./(2*pi)*log(norm(r2)/norm(r1));
     end
end