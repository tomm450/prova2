function [RES,res_x,res_y,slat_dep,Up] = rot_dwn(dwn,slat_l,alpha)

LE_i = slat_l(1,:);
slat_temp_x = slat_l(:,1) - LE_i(1);
slat_temp_y = slat_l(:,2) - LE_i(2);


slat_temp = [slat_temp_x,slat_temp_y];
clear slat_temp_x slat_temp_y

    R = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)];
    
    slat_temp2 = R * slat_temp';
    slat_dep_x = slat_temp2(1,:) + LE_i(1); 
    slat_dep_y = slat_temp2(2,:) + LE_i(2);
    
    slat_dep = [slat_dep_x',slat_dep_y'];
    clear slat_dep_x slat_dep_y slat_temp2 slat_temp
    
    TE_2 = slat_dep(end,:);
         
    [res_x_bis,i_res_x_s] = min((abs(dwn(:,1) - TE_2(1))));
    
    if i_res_x_s <= 10
        lim_left = ceil(0.25*i_res_x_s);
    else
        lim_left = i_res_x_s - 10;
    end
    
    
    lim_right = i_res_x_s + 10;
    [dummy] = thickener(dwn(lim_left:lim_right,:),TE_2(1),0,0);
    %[dummy] = thickener(dwn,TE_2(1),1,0);
    
    
    [res_x,i_res_x] = min((abs(dummy(:,1) - TE_2(1))));
    res_y = dummy(i_res_x,2) - TE_2(2);
    
    %% 
    if res_y > 0 %dwn(i_res_x_s,2) > TE_2(2)
        Up = 1; % profilo sopra slat
    else 
        Up = 0;
    end
    
    res_x = abs(res_x); res_y = abs(res_y);
    RES = sqrt(res_x^2 + res_y^2);
        
end

