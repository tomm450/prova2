function [ GEOM_out,cineq] = slat_position_x_y_a( GEOM, x, y, alpha,PLT)

if nargin == 4
PLT = 0;
end
alpha = alpha*pi/180; % <- immessa in gradi

GEOM_out = GEOM;

% TE = [0 0]
u = GEOM.slat_u_up;
d = GEOM.slat_u_dwn;
% airfoil in condizione di atterraggio
a = GEOM.up_land;


R = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

ur = R*u'; ur(1,:) = ur(1,:) + x; ur(2,:) = ur(2,:) + y;
dr = R*d'; dr(1,:) = dr(1,:) + x; dr(2,:) = dr(2,:) + y;

GEOM_out.slat_land_u = ur';
GEOM_out.slat_land_d = dr';

%d_min = 100;


for w = 1:size(GEOM_out.slat_land_d,1)
    
    [d_min_t(w),pmin(w)] = min(sqrt( (GEOM_out.slat_land_d(w,1) - a(:,1)).^2 +...
                      (GEOM_out.slat_land_d(w,2) - a(:,2)).^2));
    
    
    
end               

[d_min,pcrit_c] = min(d_min_t);
pcrit = pmin(pcrit_c);

dx = GEOM_out.slat_land_d(pcrit_c,1) - a(pcrit,1);
dy = GEOM_out.slat_land_d(pcrit_c,2) - a(pcrit,2);

L = (dx<=0)*(dy>=0);

cineq = 10 -L*d_min;


%check 

% figure;
if PLT == 1
plot(a(:,1),a(:,2),'k');
hold on
axis equal
plot(ur(1,:),ur(2,:),'r');
plot(dr(1,:),dr(2,:),'g');

plot(GEOM_out.slat_land_d(pcrit_c,1),GEOM_out.slat_land_d(pcrit_c,2),'kx');
plot(a(pcrit,1),a(pcrit,2),'gx')
end



