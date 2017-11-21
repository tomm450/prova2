function [ GEOM_out,ERR,cineq] = slat_position( GEOM, gap, r, alpha,PLT)

if nargin == 4
PLT =0;
end

alpha = alpha*pi/180; % <- immessa in gradi
toll = 1e-6;
ERR  = 100;
itermax = 100;

if abs(gap) >= r
    error('Come cavolo lo posiziono')
end
cineq = r - abs(gap);


% gap > 0 overlap
% r       raggio distanza
% alpha   incidenza corda slat

GEOM_out = GEOM;

% TE = [0 0]
u = GEOM.slat_u_up;
d = GEOM.slat_u_dwn;
% airfoil in condizione di atterraggio
a = GEOM.up_land;
% angolo pannelli
teta = atan(diff(a(:,2))./diff(a(:,1)));
% metà pannello

pcx = a(1:end-1,1) + 0.5*diff(a(:,1));
pcy = a(1:end-1,2) + 0.5*diff(a(:,2));

err = gap - ( pcx - r.*sin(teta));
[~,ier] = min(abs(err));
er = err(ier);

% figure
% plot(pcx,err,'x--')

if sign(er) == 1
    % errore positivo aumentando gli indici cambio segno
    ier2 = ier;
    while sign(err(ier2)) == 1
        ier2 = ier2 + 1;
    end
else
    % errore negativo aumentando gli indici cambio segno
    ier2 = ier;
    while sign(err(ier2)) ~= 1
        ier2 = ier2 - 1;
    end
end

% ora ho i due punti di partenza per una bisezione
ileft  = min([ier,ier2]);
pleft  = [ err(ileft) ,pcx(ileft), pcy(ileft) ,teta(ileft) ];
iright = max([ier,ier2]);
pright = [ err(iright),pcx(iright),pcy(iright),teta(iright)];


%err_v = []; left_v =[]; mid_v = []; right_v = [];
iter = 0;
while ERR > toll
    
    dx = abs(pright(2) - pleft(2))/4;
    
    [~,imid_ad] = thickener(a,0.5*(pleft(2)+pright(2)),dx,1);
    
    
    tmid = atan((imid_ad(end,2) - imid_ad(1,2))/(imid_ad(end,1) - imid_ad(1,1)));
    err_mid = gap - ( imid_ad(2,1) - r.*sin(tmid));
    
    pmid = [ err_mid ,imid_ad(2,1), imid_ad(2,2), tmid];
    
    
    
    %%left_v = [left_v,pleft(2)];
    
    %mid_v = [mid_v,pmid(2)];
    
    %right_v = [right_v,pright(2)];
    
    [ERR,ierr] = min(abs([pleft(1),pmid(1),pright(1)]));
    
    if sign(pleft(1)) == sign(pmid(1)) % lo zero è tra mid - right
        
        if sign(pmid(1)) == sign(pright(1))
            error('Segni uguali...')
        end
        pleft = pmid;
        
        
    else % lo zero è left - mid
        pright = pmid;
    end
    
    
    %err_v = [err_v,ERR];
    
    iter = iter + 1;
    if iter == itermax
        disp('NON CONVERGE')
        break
    end
    
    
end
x_v = [pleft(2),pleft(2),pleft(2)];
y_v = [pleft(3),pleft(3),pleft(3)];
t_v = [pleft(4),pleft(4),pleft(4)];

offset_x = x_v(ierr) - r*sin(t_v(ierr));
offset_y = y_v(ierr) + r*cos(t_v(ierr));

R = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

ur = R*u'; ur(1,:) = ur(1,:) + offset_x; ur(2,:) = ur(2,:) + offset_y;
dr = R*d'; dr(1,:) = dr(1,:) + offset_x; dr(2,:) = dr(2,:) + offset_y;

GEOM_out.slat_land_u = ur';
GEOM_out.slat_land_d = dr';


if PLT == 1
figure; 
plot(a(:,1),a(:,2),'k'); 
hold on
plot(ur(1,:),ur(2,:),'r');
plot(dr(1,:),dr(2,:),'g');
% plot(cx,cy,'b')
title(sprintf('Verifica: gap= %f; r=%f; alpha=%f',gap,r,alpha))
end




