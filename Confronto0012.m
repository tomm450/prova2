clean
clc
%% DATI SPERIMENTALI
l8  = load('Ladson_80.dat');
l12 = load('Ladson_120.dat');
l18 = load('Ladson_180.dat');


%MAT  = csvread('0012full.csv');
MAT = csvread('16pm.csv',1,0);
SAVE_FLAG = 0;

alpha_tested = [15];
alpha_group  = cell(1,numel(alpha_tested));
series_group = cell(1,max(MAT(:,1)));

% divido in gruppi
for i = 1:size(MAT,1)
    
    [dum,iMin] = min(abs(MAT(i,4)-alpha_tested));
    
    alpha_group{iMin}       = [alpha_group{iMin}     ; MAT(i,:)]; 
    series_group{MAT(i,1)}  = [series_group{MAT(i,1)}; MAT(i,:)]; 
    
end

%% PLT
% cl
h1 = figure(1);
plot(l8(:,1),l8(:,2),l12(:,1),l12(:,2),l18(:,1),l18(:,2),'o-','LineWidth',3)
hold on
grid on
xlabel('\alpha'); ylabel('C_l');

% solo full converge
h2 = figure(2);
plot(l8(:,1),l8(:,2),l12(:,1),l12(:,2),l18(:,1),l18(:,2),'o-','LineWidth',3)
hold on
grid on
xlabel('\alpha'); ylabel('C_l');
 
% cd
h3 = figure(3);
plot(l8(:,1),l8(:,3),l12(:,1),l12(:,3),l18(:,1),l18(:,3),'o-','LineWidth',3)
hold on
grid on
xlabel('\alpha'); ylabel('C_d');

% solo full converge
h4 = figure(4);
plot(l8(:,1),l8(:,3),l12(:,1),l12(:,3),l18(:,1),l18(:,3),'o-','LineWidth',3)
hold on
grid on
xlabel('\alpha'); ylabel('C_d');


leg14 = {'Ladson 80','Ladson 120','Ladson 180'}; 

% plot serie
for i = 1:size(series_group,2)
   
    if isempty(series_group{i}) == 0
        
        % caratteristiche
        if isnan(series_group{i}(1,52))
            % SA
            tmod = 'SA';
        else
            tmod = 'Kw';
        end
        
        % ordino per alpha calcolati
        series_group{i} = sortrows(series_group{i},4);
        
        % legenda
        leg14{end+1} = sprintf('Serie %d; l = %d; %s',...
            series_group{i}(1,1),series_group{i}(1,12),tmod);
        
        figure(1) % GLOBALE
        plot(series_group{i}(:,4),series_group{i}(:,54),plot_style(i)); hold on
        
        figure(3) % GLOBALE
        plot(series_group{i}(:,4),series_group{i}(:,47),plot_style(i)); hold on
                
        % FIGURE
        figure(1000+i) % locale
        subplot(2,2,1);
        plot(l12(:,1),l12(:,2),'g'); hold on;
        hold on; grid on;
        plot(series_group{i}(:,4),series_group{i}(:,54),'bo-','LineWidth',2);    hold on; grid on;
        
        plot(series_group{i}(:,4),series_group{i}(:,54)+0.5*series_group{i}(:,46),'b--','LineWidth',1);
        plot(series_group{i}(:,4),series_group{i}(:,54)-0.5*series_group{i}(:,46),'b--','LineWidth',1);
        
        title(sprintf('Serie %d (l = %f; %s)',series_group{i}(1,1),series_group{i}(1,12),tmod));
        
        legend('Ladson 120','CFD','Location','best')
        xlabel('\alpha');
        ylabel('C_l');
        
        subplot(2,2,2); % cl errore
        
        c_comp = spline(l12(:,1),l12(:,2),series_group{i}(:,4));
        semilogy(series_group{i}(:,4),abs(c_comp-series_group{i}(:,54)),'bo-','LineWidth',2);
        hold on; grid on;
        title('Errore rispetto NACA')
        hold on; grid on;
        title('Errore rispetto NACA')
        
        subplot(2,2,3);
        
        plot(l12(:,1),l12(:,3),'g'); hold on;
        hold on; grid on;
        plot(series_group{i}(:,4),series_group{i}(:,47),'bo-','LineWidth',2);    hold on; grid on;
        
        %plot(series_group{i}(:,4),series_group{i}(:,54)+0.5*series_group{i}(:,46),'b--','LineWidth',1);
        %plot(series_group{i}(:,4),series_group{i}(:,54)-0.5*series_group{i}(:,46),'b--','LineWidth',1);
        
        %title(sprintf('Cd-alpha Serie %d (l = %d; %s)',series_group{i}(1,1)),series_group{i}(1,14),tmod);
        
        %legend('Ladson 120','CFD','Location','best')
        xlabel('\alpha');
        ylabel('C_l');
        
        subplot(2,2,4);
        
        c_comp = spline(l12(:,1),l12(:,3),series_group{i}(:,4));
        semilogy(series_group{i}(:,4),abs(c_comp-series_group{i}(:,47)),'bo-','LineWidth',2);
        hold on; grid on;
        title('Errore rispetto NACA')
        
        ONLYConv = [nan,nan];
        ONLYConvpol = [nan,nan];
        for k = 1:size(series_group{i},1)
            
            if series_group{i}(k,56) == 0
                subplot(2,2,1);
                hold on
                plot(series_group{i}(k,4),series_group{i}(k,54),'ro','LineWidth',2);
                
                subplot(2,2,2);
                hold on
                %           semilogy(series_group{i}(k,4),abs(Cl3M(series_group{i}(k,4)+1)-series_group{i}(k,54)),'ro','LineWidth',2);
                
                subplot(2,2,3);
                hold on
                plot(series_group{i}(k,4),series_group{i}(k,47),'ro','LineWidth',2);
                
               
            else
                
                ONLYConv    = [ONLYConv;   [series_group{i}(k,4),series_group{i}(k,54)]];
                ONLYConvpol = [ONLYConvpol;[series_group{i}(k,4),series_group{i}(k,47)]];
            end
            
            
        end
        figure(2); plot(ONLYConv(:,1),ONLYConv(:,2),plot_style(i)); hold on; grid on;
        figure(4); plot(ONLYConvpol(:,1),ONLYConvpol(:,2),plot_style(i)); hold on; grid on;
        clear ONLYConv
    end
end

figure(1); legend(leg14,'Location','bestoutside');
figure(2); title('Fully converged cases'); legend(leg14,'Location','bestoutside');
figure(3); legend(leg14,'Location','bestoutside');
figure(4); title('Fully converged cases'); legend(leg14,'Location','bestoutside');
handles=findall(0,'type','figure');



for j = 1:max(size(handles))
    
    h = handles(j);
    h.Position = [50 50 1200 750];
    set(h.Children,'fontsize',12);
    if SAVE_FLAG == 1
        saveas(h,strcat('./PNG/',num2str(h.Number),'Cal','.png'));
    end
end



%close all


%%
interesse = 12; % ordine secondo l_air

for i = 1:size(alpha_group,2)
     
    % ordino secondo interesse
    alpha_group{i} = sortrows(alpha_group{i},interesse); 
    
    text_cell = alpha_group{i}(:,1)';
    text_cell = num2str(text_cell); text_cell = strsplit(text_cell,' ');
    
    dummy_x = [alpha_group{i}(1,interesse),alpha_group{i}(end,interesse)];
        
    figure(2000+i)
    subplot(1,2,1)    
    plot(alpha_group{i}(:,interesse),alpha_group{i}(:,54),'o')
    hold on
    grid on
    plot(dummy_x,ones(1,2)*spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)),'r--')
    
    title(sprintf('A = %d deg (rosso= caso non converso)',alpha_tested(i)))
    set(gca,'xdir','reverse');
    xlabel('Airfoil element dimension');
    legend('CFD results','NACA report')
    text(1.1*alpha_group{i}(:,interesse),alpha_group{i}(:,54),text_cell,'FontSize',14);
    
    axis([ 0.8*min(alpha_group{i}(:,interesse)) 1.2*max(alpha_group{i}(:,interesse)) ...
        0.8*min([alpha_group{i}(:,54);spline(l12(:,1),l12(:,2),alpha_group{i}(1,4))]) 1.2*max([alpha_group{i}(:,54);spline(l12(:,1),l12(:,2),alpha_group{i}(1,4))])]);
    
      
    subplot(1,2,2)    
    semilogy(alpha_group{i}(:,interesse),abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)),'o')
    hold on
    grid on
    set(gca,'xdir','reverse');
    xlabel('Airfoil element dimension');
    legend('|CFD results - NACA report|')
    text(1.1*alpha_group{i}(:,interesse),abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)),text_cell,'FontSize',14);
     
    axis([ 0.8*min(alpha_group{i}(:,interesse)) 1.2*max(alpha_group{i}(:,interesse)) ...
        0.8*min(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54))) 1.2*max(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)))]);
    
    
    for k = 1:size(alpha_group{i},1)
       
       if alpha_group{i}(k,56) == 0
           subplot(1,2,1);
           hold on
           plot(alpha_group{i}(k,interesse),alpha_group{i}(k,54),'ro','LineWidth',2);
           
           subplot(1,2,2);
           hold on
           semilogy( alpha_group{i}(k,interesse),abs(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(k,54))),'ro','LineWidth',2);
           
       end
    end
    
    figure(4000+i)
    
    subplot(1,3,1)
    semilogy(alpha_group{i}(:,interesse),alpha_group{i}(:,49),'o');
    text(1.1*alpha_group{i}(:,interesse),alpha_group{i}(:,49),text_cell,'FontSize',14);
    set(gca,'xdir','reverse'); grid on
    ylabel('Res(Ux)'); title(sprintf('A = %d deg',alpha_tested(i)))
    
    subplot(1,3,2)
    semilogy(alpha_group{i}(:,interesse),alpha_group{i}(:,50),'o');
    text(1.1*alpha_group{i}(:,interesse),alpha_group{i}(:,50),text_cell,'FontSize',14);
    set(gca,'xdir','reverse'); grid on
    ylabel('Res(Uz)')
    
    subplot(1,3,3)
    semilogy(alpha_group{i}(:,interesse),alpha_group{i}(:,51),'o');
    text(1.1*alpha_group{i}(:,interesse),alpha_group{i}(:,51),text_cell,'FontSize',14);
    set(gca,'xdir','reverse'); grid on
    ylabel('Res(P)')
end
    


handles=findall(0,'type','figure');

for j = 1:max(size(handles))
    
    h = handles(j);
    h.Position = [50 50 900 750];
     set(h.Children,'fontsize',12);
    if SAVE_FLAG == 1     
     saveas(h,strcat('./PNG/',num2str(h.Number),'Cal','.png'));
end
end
%close all

%% 
interesse = 35; % ordine secondo n_face

for i = 1:size(alpha_group,2)
       
    % ordino secondo interesse
    alpha_group{i} = sortrows(alpha_group{i},interesse); 
    
    text_cell = alpha_group{i}(:,1)';
    text_cell = num2str(text_cell); text_cell = strsplit(text_cell,' ');
    
    dummy_x = [alpha_group{i}(1,interesse),alpha_group{i}(end,interesse)];
        
    figure(3000+i)
    subplot(1,2,1)    
    semilogx(alpha_group{i}(:,interesse),alpha_group{i}(:,54),'o')
    hold on
    grid on
    semilogx(dummy_x,ones(1,2)*spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)),'r--')
    title(sprintf('A = %d deg (rosso= caso non converso)',alpha_tested(i)))
    %set(gca,'xdir','reverse');
    xlabel('nFace')
    legend('CFD results','NACA report')
%     title(sprintf('Cl (rosso= caso non converso)'));

    text(1.1*alpha_group{i}(:,interesse),alpha_group{i}(:,54),text_cell,'FontSize',14);
    
    axis([ 0.8*min(alpha_group{i}(:,interesse)) 1.2*max(alpha_group{i}(:,interesse)) ...
        0.8*min([alpha_group{i}(:,54);spline(l12(:,1),l12(:,2),alpha_group{i}(1,4))]) 1.2*max([alpha_group{i}(:,54);spline(l12(:,1),l12(:,2),alpha_group{i}(1,4))])]);
   
    
    subplot(1,2,2)    
    loglog(alpha_group{i}(:,interesse),abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)),'o')
    hold on
    grid on
    %set(gca,'xdir','reverse');
    xlabel('nFace');
    legend('|CFD results - NACA report|')
    text(1.1*alpha_group{i}(:,interesse),abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)),text_cell,'FontSize',14);
     
    axis([ 0.8*min(alpha_group{i}(:,interesse)) 1.2*max(alpha_group{i}(:,interesse)) ...
        0.8*min(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54))) 1.2*max(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(:,54)))]);
    
    
    for k = 1:size(alpha_group{i},1)
       
       if alpha_group{i}(k,56) == 0
           subplot(1,2,1);
           hold on
           semilogx(alpha_group{i}(k,interesse),alpha_group{i}(k,54),'ro','LineWidth',2);
           
           subplot(1,2,2);
           hold on
           semilogx( alpha_group{i}(k,interesse),abs(abs(spline(l12(:,1),l12(:,2),alpha_group{i}(1,4)) - alpha_group{i}(k,54))),'ro','LineWidth',2);
           
       end
    end
    
end



handles=findall(0,'type','figure');

for j = 1:max(size(handles))
    
    h = handles(j);
    h.Position = [50 50 900 750];
    set(h.Children,'fontsize',12);
    if SAVE_FLAG == 1
        saveas(h,strcat('./PNG/',num2str(h.Number),'Cal','.png'));
    end
end


%close all
