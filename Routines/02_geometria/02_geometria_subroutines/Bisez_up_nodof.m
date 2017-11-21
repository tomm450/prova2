function [slat_u_dep,d_up,n_rot1,res_y] = Bisez_up(l_u,a_u,pivot,slat_u_temp,toll)

% confuso ma funziona
n_rot1 = 0;  
n_max  = 100;
d_up_vect = [-pi/2 0 pi/2];

res_v = (l_u.* sin(a_u + d_up_vect)) - pivot(2); 

while min(abs(res_v)) > toll %&& n_rot1 < n_max
    
    % se res_v < 0 aumento angolo                                                                %
    if sign(res_v(1)) == sign(res_v(2))
        if sign(res_v(2)) == sign(res_v(3))
            l_u
            pivot
            n_rot1
            res_v
            error('Bisezione non funziona, segno uguale in tutti i punti (Bisez_up.m)');
        else
            % segno1 = segno2 -> diverso da segno 3
            d_up_vect = [d_up_vect(2) 0.5*(d_up_vect(2)+d_up_vect(3)) d_up_vect(3)];
        end
    else
        if sign(res_v(2)) == sign(res_v(3))
            % segno 1 =~ segno 2 = segno3
            d_up_vect = [d_up_vect(1) 0.5*(d_up_vect(1)+d_up_vect(2)) d_up_vect(2)];
        else
            n_rot1
            res_v
            error('Bisezione non funziona, segno uguale in tutti i punti (Bisez_up.m)');
        end
    end
    
    if n_rot1 == n_max + 1
        break
    end
    
    for k = 1:size(d_up_vect,2)
       R = [cos(d_up_vect(k)), -sin(d_up_vect(k)) ; sin(d_up_vect(k)), cos(d_up_vect(k))]; 
       slat_u_temp2 = R * slat_u_temp(1,:)'; % ottengo LE = [LE_x; LE_y]
       res_v(k) = slat_u_temp2(2)+ pivot(2);
    end
    
    n_rot1 = n_rot1 +1;  
end



% Ruoto su minimo residuo
[res_y,i_min] = min(res_v);
d_up = d_up_vect(i_min);
R = [cos(d_up) -sin(d_up) ; sin(d_up) cos(d_up)];
    
slat_u_temp2 = R * slat_u_temp';
slat_u_dep_x = slat_u_temp2(1,:) + pivot(1); slat_u_dep_y = slat_u_temp2(2,:) + pivot(2);
slat_u_dep = [slat_u_dep_x',slat_u_dep_y'];    

if slat_u_dep_y(1) >= toll
    slat_u_dep_y
    fprintf('Residio fuori tolleranza -> Y_LE = %6.6f vs toll = %6.6f\n',slat_u_dep_y(1),toll);
    error('Qualquadra non cosa (Bisez_up.m)')
end

