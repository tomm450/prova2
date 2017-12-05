function [slat_u_dep,d_up,n_rot,res] = Bisez_dwn(dwn,slat_l_temp,alpha_dwn,toll)
% alpha_s = angolo prerotazione 
n_rot = 0;
n_max  = 300;
RES = ones(1,3); 
res_x = RES; 
res_y = RES;
Up     = [ 2 2 2];
Up_old = [ 2 2 2];
res_vect = zeros(3,n_max);
alpha_dwn_mat = zeros(3,n_max);


% % %verifica condizioni iniziali
% [~,~,~,slat_u_dep1,up1] = rot_dwn(dwn,slat_l_temp,alpha_dwn(1));%,toll);
% [~,~,~,slat_u_dep2,up2] = rot_dwn(dwn,slat_l_temp,alpha_dwn(2));%,toll);
% [~,~,~,slat_u_dep3,up3] = rot_dwn(dwn,slat_l_temp,alpha_dwn(3));%,toll);
% 
% %     figure(99)
% %     hold on
% %     plot(slat_u_dep1(:,1),slat_u_dep1(:,2),'k')
% %     hold on
% %     plot(slat_u_dep2(:,1),slat_u_dep2(:,2),'g')
% %     plot(slat_u_dep3(:,1),slat_u_dep3(:,2),'b')
% %     plot(dwn(:,1),dwn(:,2))
% %     axis equal
% %     legend(num2str(alpha_dwn(1)),num2str(alpha_dwn(2)),num2str(alpha_dwn(3)));
% %     title(sprintf('%d %d %d',up1,up2,up3))

%while min(abs(RES)) >= toll
while min(abs(res_y)) >= toll
    n_rot = n_rot+1;

    for k = 1:3
        [RES(k),res_x(k),res_y(k),~,Up(k)] = rot_dwn(dwn,[slat_l_temp(1,:);slat_l_temp(end,:)],alpha_dwn(k));
        %[RES(k),res_x(k),res_y(k),~,Up(k)] = rot_dwn(dwn,slat_l_temp,alpha_dwn(k));
    end
    
    if min(abs(res_y)) < toll
        break
    end
    
    if Up(1) ~= Up(2) && Up(2) == Up(3)
        a1 = alpha_dwn(1); a2 = 0.5*(alpha_dwn(1) + alpha_dwn(2)); a3 = alpha_dwn(2);
    elseif Up(1) == Up(2) && Up(2) ~= Up(3)
        a1 = alpha_dwn(2); a2 = 0.5*(alpha_dwn(3) + alpha_dwn(2)); a3 = alpha_dwn(3);
        
    else %Up(1) == Up(2) && Up(2) == Up(3)
        fprintf('%d\n',Up)
        fprintf('%d\n',n_rot)
        error('Debug necessario Bisez_dwn.m')
    end
    
    alpha_dwn = [a1,a2,a3];
    
    Up_old = Up;
    res_vect(:,n_rot) = res_y;
    alpha_dwn_mat(:,n_rot) = [alpha_dwn(1);alpha_dwn(2);alpha_dwn(3)];
    
    if n_rot == n_max
        break
    end
end

% Ruoto verso il minimo residup
[res,i_min]=min(abs(res_y));
d_up = alpha_dwn(i_min);
[~,~,~,slat_u_dep,~] = rot_dwn(dwn,slat_l_temp,alpha_dwn(i_min));

% sdeng
if n_rot == n_max  % NON CONVERSO
    %sort_a = sort(alpha_dwn);
    %if abs(sort_a(1) - sort_a(end)) < toll^2
    % tutto bene
    %else
    format long
    alpha_dwn
    RES
    Up
    Up_old
    res_y
    res_x
    alpha_dwn
    figure(1)
    hold on
    plot(slat_u_dep(:,1),slat_u_dep(:,2),'kx')
    plot(dwn(:,1),dwn(:,2))
    axis equal
    fprintf('controllare residui,iterazione %d',n_rot)
    figure
    plot(res_vect(1,:))
    hold on
    plot(res_vect(2,:))
    plot(res_vect(3,:))
    title('Residuo y')
    figure
    plot(alpha_dwn_mat(1,:))
    hold on
    plot(alpha_dwn_mat(2,:))
    plot(alpha_dwn_mat(3,:))
    title('Alpha')
    error('non conv')
end
end
