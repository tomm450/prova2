
clean
load('./Output/8dump.mat');
close all
% ULTIMO CASO

itc = max(size(f_cell));

X_mat = nan(itc,3);

CL_cfd  = nan(3,itc); CD_cfd  = nan(3,itc);
CL_hs   = nan(1,itc); CD_hs   = nan(1,itc); 
CL_hsc = nan(1,itc);  CD_hsc = nan(1,itc);

ineq_cfd = nan(2,itc);
ineq_hs  = nan(2,itc);
ineq_hsc = nan(2,itc);

p_cfd = nan(400,itc);
p_hs  = nan(400,itc);
p_hsc = nan(400,itc);



for i = 1:max(itc)

    if sum(WIN_OPT{i}{1}>OPT.lb) == 3 &&...
            sum(WIN_OPT{i}{1}<OPT.ub) == 3 && ...
        OPT.Aineq*WIN_OPT{i}{1}'<OPT.bineq
        feas(i) = 1;
    end
    X_mat(i,:) = (WIN_OPT{i}{1}-OPT.lb)./(OPT.ub-OPT.lb);
    
    % LIFT
    CL_cfd(:,i) = [f_cell{i}{1}(1)+f_cell{i}{end}(1);...
                   f_cell{i}{1}(1);...
                   f_cell{i}{1}(1)-f_cell{i}{end}(1)];
   
    CL_hs(:,i)    = c_cell{i}{1}(1);
    CL_hsc(:,i)  = copt_cell{i}{1}(1);
    
    % DRAG
    CD_cfd(:,i) = [f_cell{i}{1}(2)+f_cell{i}{end}(2);...
                   f_cell{i}{1}(2);...
                   f_cell{i}{1}(2)-f_cell{i}{end}(2)];
   
    CD_hs(:,i)    = c_cell{i}{1}(2);
    CD_hsc(:,i)   = copt_cell{i}{1}(2);
   
    % PRESSURE
    p_cfd(:,i) = f_cell{i}{2}';
    p_hs(:,i)  = c_cell{i}{2}';
    p_hsc(:,i) = copt_cell{i}{2}';
    
    % Dp
    ineq_cfd(:,i) = f_cell{i}{4}';
    ineq_hs(:,i)  = c_cell{i}{4}';
    ineq_hsc(:,i) = copt_cell{i}{4}';
    
    
%     %% PLOT
%     cp_c    = c_cell{i}{2};
%     cp_copt = copt_cell{i}{2};
%     cp_f    = f_cell{i}{2};
%         xc = usefulCoor{i}{1};
%     n = 0.5*max(size(xc));
%     xc_m_plot = [xc(1,1:n) ,nan,xc(2,1:n)];
%     
%     cp_m_c    = [cp_c(1:n)',   nan,cp_c(2*n+1:3*n)'];
%     cp_m_corr = [cp_copt(1:n)',nan,cp_copt(2*n+1:3*n)'];
%     cp_m_foam = [cp_f(1:n)',   nan,cp_f(2*n+1:3*n)'];
%     
%     xc_c_plot = [xc(1,n:end),    nan,xc(2,n:end)];        
%     cp_c_c    = [cp_c(n:2*n)',   nan,cp_c(3*n:end)'];
%     cp_c_corr = [cp_copt(n:2*n)',nan,cp_copt(3*n:end)'];
%     cp_c_foam = [cp_f(n:2*n)',   nan,cp_f(3*n:end)'];
%     
%     
%     %% PRESSURE DISTRIBUTION
%     h = figure(1000+i);
%     
%     set(h,'Position',[0 0 1350 750]);
%     subplot(2,2,1)
%    
%     plot(xc_m_plot,cp_m_c,'ro-','LineWidth',1);
%     hold on
%     plot(xc_c_plot,cp_c_c,'r>-','LineWidth',1);
%     
%     plot(xc_m_plot,cp_m_corr,'g*-','LineWidth',1);
%     plot(xc_c_plot,cp_c_corr,'g+-','LineWidth',1);
%     
%     plot(xc_m_plot,cp_m_foam,'b^-','LineWidth',1);
%     plot(xc_c_plot,cp_c_foam,'bs-','LineWidth',1);
%     
%     grid on
%     
%     set(gca,'YDir','Reverse'); % reverse y axis
%     %     alfa=num2str(alpha*180/pi);
%     %     tit1 = 'COEFFICIENTE DI PRESSIONE';
%     %     tit2 = ['Incidenza: ',BU_par.alpha, ' gradi'];
%     %     title({tit1;tit2})
%     title(sprintf('Iteration %d',i))
%     xlabel('X coordinate');
%     ylabel('-C_p');
%     legend('Dorso HS','Ventre HS',...
%            'Dorso HS_c','Ventre HS_c',...
%            'Dorso CFD','Ventre CFD')
%     
%     subplot(2,2,2)
%        
%     semilogy(xc_m_plot,abs(cp_m_foam-cp_m_c),'mo-','LineWidth',1);
%     hold on
%     semilogy(xc_c_plot,abs(cp_c_foam-cp_c_c),'m>-','LineWidth',1);
%     
%     semilogy(xc_m_plot,abs(cp_m_foam-cp_m_corr),'c*-','LineWidth',1);
%     semilogy(xc_c_plot,abs(cp_c_foam-cp_c_corr),'c+-','LineWidth',1);
%     
%     grid on
%     xlabel('X coordinate');
%     ylabel('| C_p - C_p_{cfd} |');
%     legend('Dorso HS-CFD','Ventre HS-CFD','Dorso HS_c-CFD','Ventre HS_c-CFD');%,'Dorso CFD','Ventre CFD')
%     
%     
%     subplot(2,2,3)
%     plot(cp_f,'bo-')
%     hold on
%     plot(cp_c,'ro-')
%     grid on
%     plot(cp_copt,'go-');
%     grid on
%     set(gca,'YDir','Reverse'); % reverse y axis
%   
%        xlabel('Panel index');
%     ylabel('-C_p');
%     legend('CFD','HS','HS_c')
%     %
%     subplot(2,2,4)
%     semilogy(abs(cp_f-cp_c),'mo-'); grid on; hold on;% title('Pressure difference error');
%     semilogy(abs(cp_f-cp_copt),'co-');
%     legend('HS vs CFD','HS_c vs CFD');
%     xlabel('X coordinate');

end

x_ax = [1:itc];

%% CL
figure(10);
subplot(1,2,1)
plot(x_ax,CL_hs,'ro-',x_ax,CL_hsc,'go-')
hold on
plot(x_ax,CL_cfd(2,:),'bo-',x_ax,CL_cfd(1,:),'bo--',x_ax,CL_cfd(3,:),'bo--');
hold on
grid on
title('C_l evolution');
xlabel('OMM iteration');
ylabel('C_l');
legend('HS','HS corr','CFD')


subplot(1,2,2)
semilogy(x_ax,abs(CL_cfd(2,:)-CL_hs),'mo-',...
     x_ax,abs(CL_cfd(2,:)-CL_hsc),'co-')

grid on
title('Error of C_l evolution');
xlabel('OMM iteration');
ylabel('| C_l - C_l_{ cfd}|');
legend('|HS - CFD|','|HS_{corr}-CFD|')
%% CD
figure(11);
subplot(1,2,1)
plot(x_ax,CD_hs,'rs-',x_ax,CD_hsc,'gs-')
hold on
plot(x_ax,CD_cfd(2,:),'bs-',x_ax,CD_cfd(1,:),'bs--',x_ax,CD_cfd(3,:),'bs--');
hold on
grid on
title('C_d evolution');
xlabel('OMM iteration');
ylabel('C_d');
legend('HS','HS corr','CFD')

subplot(1,2,2)
semilogy(x_ax,abs(CD_cfd(2,:)-CD_hs),'ms-',...
     x_ax,abs(CD_cfd(2,:)-CD_hsc),'cs-')

grid on
title('Error of C_l evolution');
xlabel('OMM iteration');
ylabel('| C_d - C_d_{ cfd}|');
legend('|HS - CFD|','|HS_{corr}-CFD|')

%% Xwin

figure(20); plot(X_mat,'o-')

title('Normalized component of the X vector')
xlabel('OMM iteration');
grid on

%%
figure(30); 
plot(ineq_cfd(1,:),'bo-');
hold on
 plot(ineq_hs(1,:),'ro-');
hold on
plot(ineq_hsc(1,:),'go-');
hold on
legend('DP slat CFD','DP slat HS','DP slat HSc')
xlabel('OMM iteration');
grid on

figure(31); plot(ineq_cfd(2,:),'bo-');
hold on
 plot(ineq_hs(2,:),'ro-');
hold on
plot(ineq_hsc(2,:),'go-');
hold on
legend('DP airfoil CFD','DP airfoil HS','DP airfoil HSc')
xlabel('OMM iteration');
grid on


for i = 1:max(size(FULL_OPT))
   hi = figure(5000+i);
   set(hi,'Position',[0 0 1320 750])
    leg_cell = {};
        
    for j = 1:max(size(FULL_OPT{i}{1}))
        subplot(3,3,7)
        plot([OPT.x0(j,1) FULL_OPT{i}{1}{j}(1)],...
             [OPT.x0(j,3) FULL_OPT{i}{1}{j}(3)]);
        grid on
        hold on
         
        subplot(3,3,8)
        plot(...
            [OPT.x0(j,2) FULL_OPT{i}{1}{j}(2)],...
            [OPT.x0(j,3) FULL_OPT{i}{1}{j}(3)]);
        grid on
        hold on

              
        subplot(3,3,1:6)
        plot3([OPT.x0(j,1) FULL_OPT{i}{1}{j}(1)],...
            [OPT.x0(j,2) FULL_OPT{i}{1}{j}(2)],...
            [OPT.x0(j,3) FULL_OPT{i}{1}{j}(3)]);
        grid on
        hold on
        
        leg_cell{end+1} = num2str(j);
        
    end
    
    
        
        for j = 1:max(size(FULL_OPT{i}{1}))
            if j >9
                clr = 'ro';
            else
                clr = 'ko';
            end
            
            subplot(3,3,7)
            plot(OPT.x0(j,1),OPT.x0(j,3),clr,'LineWidth',2);
            hold on
            
            subplot(3,3,8)
            plot(OPT.x0(j,2),OPT.x0(j,3),clr,'LineWidth',2);
            hold on
            
            subplot(3,3,1:6)
            
            plot3(OPT.x0(j,1),OPT.x0(j,2),OPT.x0(j,3),clr,'LineWidth',3);
            dummy = 3;
            
        end
        
    subplot(3,3,7)
     hold on
    %    legend(leg_cell); 
        xlabel('x position');
    %    ylabel('y position');
        ylabel('slat incidence')
    
    subplot(3,3,8)
     hold on
    %    legend(leg_cell); 
    %    xlabel('x position');
        xlabel('y position');
        ylabel('slat incidence')
    
    subplot(3,3,1:6)
     hold on
        legend(leg_cell); 
        xlabel('x position');
        ylabel('y position');
        zlabel('slat incidence')
        view([-37.5 -30])
    subplot(3,3,9)
    title('Fitness function')
    bar(-FULL_OPT{i}{2})
end



