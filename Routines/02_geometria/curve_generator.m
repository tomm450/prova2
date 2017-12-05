function [new_curve,win,slat_ret,th_per] = curve_generator(original_curve,COM_POINTS,limit_index,soprasotto)
%% versione parabola
%CURVE_GENERATOR (WORK IN PROGRESS)
% original_curve = matrice n righe 2 colonne con coordinate x e y
% COM_POINTS => vuoto tiro retta
%               size(2,1) con valori compresi tra 0 e 1 indicano le
%               cordinate normalizzate del terzo punto (parabola)
%               ...
% limit_index = vettore contenente indice inizio e fine modifica
% soprasotto  == 1 se sto facendo slat sopra
%             == 2 se sto facendo slat sotto

% new_curve = curva uscita matrice n righe 2 colonne
% win == 1 se vincoli non violati
%     == 0 se fail
% slat_ret = occhio che qui si fa particolare...
%   slat_ret = matrice kx2 -> k=1 e end sono il LE e TE (coincidono sia per curva sopra che sotto)
%                             k pari -> punti curva interna
%                             k dispari -> punti curva esterna


% isolo tronconi
A = original_curve(1:limit_index(1)-1,:);
vx = original_curve(limit_index(1):limit_index(2),1);
vy = original_curve(limit_index(1):limit_index(2),2);
C = original_curve(limit_index(2)+1:end,:);

win = 1;

if size(COM_POINTS,1) == 0 % esempio per retta
    % non fare controllo
elseif size(COM_POINTS) == [1,2]
    % parabola
    
    if COM_POINTS(1) < 0
        win = 0;
    end
    if COM_POINTS(1) > 1
        win = 0;
    end
    if COM_POINTS(2) < 0
        win = 0;
    end
    if COM_POINTS(2) > 1
        win = 0;
    end
    
end


if win == 0 % controlli non passati
    
    new_curve = [0 0];
    slat_ret  = [0 0];
    th_per = 0;
elseif win == 1 % GO
    
    % giro
    a_o = atan((vy(end)-vy(1))/(vx(end)-vx(1)));
    a_o = -a_o;
    R = [cos(a_o) -sin(a_o); sin(a_o) cos(a_o)];
            
    VX_G = R*[(vx-vx(1))';(vy-vy(1))']; % ruoto rispetto a origine
    
    vx_g = VX_G(1,:); 
    vy_g = VX_G(2,:);
     
    
    if abs(vy_g(1) - vy_g(end)) < 1e-8
        %ok
    else
        error('problemi')
    end
    
    vx_pt = COM_POINTS(1) *(vx_g(end) - vx_g(1));
    vy_pt = COM_POINTS(2) *spline(vx_g,vy_g, vx_pt );
    
    
    a_o = -a_o;
    R = [cos(a_o) -sin(a_o); sin(a_o) cos(a_o)];
    %[vx_pt;vy_pt] = 
    V_MAT = R*[vx_pt;vy_pt];
    
    
    vx_pt = V_MAT(1,:);
    vy_pt = V_MAT(2,:);
    
    vx_pt = vx_pt + vx(1);
    vy_pt = vy_pt + vy(1);
      
      
    % coordinata x per interpolazione
    bx = original_curve(limit_index(1):limit_index(2),1); % colonna
    
      
    % 3 punti
    p = [original_curve(limit_index(1),:) ; [vx_pt vy_pt]; original_curve(limit_index(2),:)];
    
    PARABOL_A = [1 p(1,1) p(1,1)^2;
                 1 p(2,1) p(2,1)^2;
                 1 p(3,1) p(3,1)^2];
         
    PARABOL_b = p(:,2);     
    
    c_vect = PARABOL_A\PARABOL_b;
    c1 = c_vect(1);
    c2 = c_vect(2);
    c3 = c_vect(3);
    
    
    PARABOLI_in = @(x) c1 + c2.*x +c3.*x.^2;
    by = PARABOLI_in(bx);
    
    % versione 0
    %by = spline(p(:,1),p(:,2),bx);
    %by = interp1(p(:,1),p(:,2),bx);
    B = [bx,by];
    
%     figure(122)
%     plot(vx_g,vy_g)
%     hold on
%     plot(vx_pt,vy_pt,'x')
    
%     figure(123)
%     plot(bx,by,'ro--'); hold on
%     plot(vx,vy,'bo--')
%     plot(p(1,1),p(1,2),'ko')
%     plot(p(2,1),p(2,2),'ko')
%     plot(p(3,1),p(3,2),'ko')
    
    
    new_curve = [A;B;C];
    
    if soprasotto == 1
        if sum(B(2:end-1,2) > vy(2:end-1)) >= 1
            win = 0;

        end
        
    elseif soprasotto == 2
        if sum(B(2:end-1,2) < vy(2:end-1)) >= 1
            win = 0;
        end
    end
    
    
    % ordino slat con spessore
    
    ext = original_curve(limit_index(1):limit_index(2),:);
    in  = B(2:end-1,:);
    
    slat_ret = zeros(2+2*size(in,1),2);
    
    slat_ret(end,:) = ext(end,:);
    
    slat_ret(2:2:end-1,:) = in;
    
    slat_ret(1:2:end-1,:) = ext(1:end-1,:);
    
    
    a_o = -a_o;
    R = [cos(a_o) -sin(a_o); sin(a_o) cos(a_o)];
    aeroref_ext = R*ext';
    aeroref_in = R*in';
    
    % vincolo spessore
     aero_ref_x   = linspace(min(aeroref_ext(1,:)),max(aeroref_ext(1,:)),100);
     aeroint_in   = spline(aeroref_in(1,:) ,aeroref_in(2,:) ,aero_ref_x);
     aeroint_ext  = spline(aeroref_ext(1,:),aeroref_ext(2,:),aero_ref_x);
    
    th     = max(abs(aeroint_ext - aeroint_in));
    c_temp = max(aero_ref_x)-min(aero_ref_x);
    th_per = th*100/c_temp;
    
%     plot(aeroref_ext(1,:),aeroref_ext(2,:),'go')
%     hold on
%     plot(aeroref_in(1,:),aeroref_in(2,:),'rx')
%     %
%     figure
%     plot(aero_ref_x,aeroint_in,'go');
%     hold on
%     plot(aero_ref_x,aeroint_ext,'rx');



    
end




