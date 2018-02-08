clean

% %% 4 profili generici
[xx(1,:),yy(1,:)] = NACA_generator('NACA4412',50,'cos');
[xx(2,:),yy(2,:)] = NACA_generator('NACA0015',50,'cos');
[xx(3,:),yy(3,:)] = NACA_generator('NACA5418',50,'cos');
[xx(4,:),yy(4,:)] = NACA_generator('NACA0012',50,'cos');
[xx(5,:),yy(5,:)] = NACA_generator('NACA0012',50,'cos');
[xx(6,:),yy(6,:)] = NACA_generator('NACA0012',50,'cos');
[xx(7,:),yy(7,:)] = NACA_generator('NACA0012',50,'cos');
close all

R = @(x) [cosd(x) sind(x)
         -sind(x) cosd(x)];
     
XX = xx; YY = yy;
%         s1   m   f1   f2  s2 f3
theta = [-30   0    10   30 -50     45  60];
dx    = [-0.25 0  1.05  1.2 -0.35  1.4  1.55];
dy    = [-0.1  0 -0.05 -0.1 -0.15 -0.25 -0.5];
K     = [0.25  1  0.15  0.2  0.15  0.3   0.2]  ;

for i = 1:size(XX,1)
    
   temp = [XX(i,:);YY(i,:)];
   
   temp = K(i)*(R(theta(i))*temp);
   
   XX(i,:) = temp(1,:) + dx(i);
   YY(i,:) = temp(2,:) + dy(i);
   
   plot(XX(i,:),YY(i,:)); hold on; axis equal
end


load('4airfoil.mat')

temp = {GEOM.slat_land_u,GEOM.slat_land_d};
GEOM.slat_land_u = {}; GEOM.slat_land_d = {};
GEOM.slat_land_u{1} = temp{1}; GEOM.slat_land_d{1} = temp{2};

% slat2 indice 5
temp = [XX(5,:)',YY(5,:)'];
temp_up = temp(50:end,:);
temp_dwn = flipud(temp(1:50,:));
GEOM.slat_land_u{2} = 1000*temp_up;
GEOM.slat_land_d{2} = 1000*temp_dwn;

save('5airfoil.mat','GEOM');



% flap3 indice 6
temp = [XX(6,:)',YY(6,:)'];
temp_up = temp(50:end,:);
temp_dwn = flipud(temp(1:50,:));
GEOM.flap_land_u{3} = 1000*temp_up;
GEOM.flap_land_d{3} = 1000*temp_dwn;

save('6airfoil.mat','GEOM');

temp = [XX(7,:)',YY(7,:)'];
temp_up = temp(50:end,:);
temp_dwn = flipud(temp(1:50,:));
GEOM.flap_land_u{4} = 1000*temp_up;
GEOM.flap_land_d{4} = 1000*temp_dwn;

save('7airfoil.mat','GEOM');
