function [lastpoint_new] = refineModuleGmsh(method,method_par,fid,lastpoint,l_airfoil,x_dom,l_dom)
%
% [lastpoint_new] = refineModuleGmsh(method,method_par,fid,lastpoint,l_airfoil,x_dom,l_dom)
% INPUT:
% METHOD è una cella contenente stringhe con i vari casi da applicare, se 
% la cella contiene più di un elemento le procedure vanno eseguite in
% sequenza
% METHOD_PAR è una cella della stessa dimensione di METHOD, ogni elemento
% di METHOD_PAR è a sua volta una cella con i mari parametri necessari per
% eseguire i case
% FID è il fopen del file .geo che si sta scrivendo
% LASTPOINT è l'indice dell'ultimo punto usato nel file .geo
% l_airfoil è la dimesione della cella sul profilo
% x_dom     è la distanza del farfiel
% l_dom     è la dimesione della cella sul farfield
% OUTPUT:
% LASTPOINT_NEW è l'indice dell'ultimo punto usato nel file .geo (dopo
% l'applicazione dei case richiesti)
% CASE implementati
% case 'none' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {0}
% non fa niente, par è una cella contenente in double a caso
% case 'clock_simple' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {rref,lref} 
% 24 punti su circonferenza di raggio rref di dimensione lref
% case 'sublinear' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {ellx,elly,rref,lref,nstaz}
% 24 punti su nstaz circonferenze (ellx e elly coefficienti per rendere 
% ellisse), presa la retta 
% y = l_air + (l_dom-l_air)/(x_dom-0.75)*(x-0.75)
% la dimesione sulle stazione segue la spline
% l = spline([0.75 rref x_dom],[l_air l_ref*y(rref) l_dom],x)
% case 'wake' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {y_c,y_w,r,l_w}
% posiziona punti di dimensione l_w sul "trapezio" compreso tra i punti 
% A = [K*x_dom y_c-r],B = [0 y_c-r],C = [0 y_c+r],D = [K*x_dom y_w]
% in cui tra B e C è presente un raccordo di raggio r, K è settato a 2/3
% case 'wake2' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {x_te,y_te,alpha,r,l_w}
% partendo da [x_te,y_te] si tracci una retta di pendenza tan(alpha) che 
% si estende fino a K*x_dom,su questa "scia" vengono posizionati punti con 
% "lunghezza" l_w/10; sommando agli estremi +-[0 r] si ottiene il
% romboide sul cui perimentro si mettono i punti di "lunghezza" 
% l = l_w+spline([x_te,K*x_dom],[0 Coef],x), come nel caso precedente 
% case 'sublinearwake' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par: {y_c,y_w,r,l_w,ellx,elly,rref,lref,nstaz}
% unisce i casi "sublinear" e "wake" facendo si che non pestino i piedi
%
% last update 24/02/2018
    
switch lower(method)
    case {'none'}
        
        lastpoint_new = lastpoint;
        
    case {'clock_simple'}
        
        % 24 punti su circonferenza di raggio rref di dimensione lref
        rref = method_par{1};
        lref = method_par{2};
        
        clock_angle = linspace(0,2*pi,24);
        clock_angle = clock_angle(1:end-1);
             
        for w = 1:23
                       
            lastpoint = lastpoint +1;
            fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,0.5+rref*cos(clock_angle(w)),...
                rref*sin(clock_angle(w)),lref);
            fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
            
        end
        
        lastpoint = lastpoint +1;
        
%         %fprintf(fid,'Point(%d) = { 1.001,  0.0000000,  0,  %f};\n',lastpoint,l_airfoil);
%         %fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
        
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { 1.001,  0.001000,  0,  %f};\n',lastpoint,l_airfoil);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
%         
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { 1.0001, -0.001000,  0,  %f};\n',lastpoint,l_airfoil);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
%         
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { -0.001,  0.0000000,  0,  %f};\n',lastpoint,l_airfoil/5);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
        
        lastpoint_new = lastpoint;
        
    case 'sublinear'
        
        % parametri per forma ellisse
        ellx = method_par{1};
        elly = method_par{2};
        
        % raggio interpolazione lineare
        rref  = method_par{3}; 
        % coefficiente tc nuova lunghezza a rref sia lref*l_linear
        lref  = method_par{4};
        % numero circonferenze da marchiare
        nstaz = method_par{5};
        
        clock_angle = linspace(0,2*pi,24);
        clock_angle = clock_angle(1:end-1);
        
        l_linear = spline([0.75 x_dom],[l_airfoil l_dom],rref);
        
        r_inter = linspace(0.75,x_dom,nstaz); 
        r_inter = r_inter(2:end-1);
        
        l_inter = spline([0.75 rref x_dom],[l_airfoil lref*l_linear l_dom],r_inter);
                
        spacing = figure(22); plot(r_inter,l_inter,'o',[0.75 x_dom],[l_airfoil l_dom],'--',rref,l_linear,'*');title('Size spacing');grid on;
        %savefig(spacing,strcat(case_dir,'/spacing.fig'));
        pause(1); close 22
        
        for w = 1:23
            
            for ww = 1:max(size(l_inter))
                lastpoint = lastpoint +1;
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,0.5+ellx*r_inter(ww)*cos(clock_angle(w)),...
                    elly*r_inter(ww)*sin(clock_angle(w)),l_inter(ww));
                fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
            end
            
        end
        
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { 1.001,  0.001000,  0,  %f};\n',lastpoint,l_airfoil);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
%         
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { 1.0001, -0.001000,  0,  %f};\n',lastpoint,l_airfoil);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
%         
%         lastpoint = lastpoint +1;
%         fprintf(fid,'Point(%d) = { -0.001,  0.0000000,  0,  %f};\n',lastpoint,l_airfoil/5);
%         fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);

        
        lastpoint_new = lastpoint;
        
     case 'wake'
        
%       posizione verticale centro raccordo
%       y_c = method_par{1};
%       altezza scia all'outlet
%       y_w = method_par{2};
%       raggio raccordo 
%       r   = method_par{3};
%       lunghezza elementi
%       l_w = method_par{4};
        
        y_c = method_par{1};
        y_w = method_par{2};
        r   = method_par{3};
        l_w = method_par{4};
        
        
        K = 2/3;
       
        a1 = (5*pi/180)+atan2(y_w-r,K*x_dom);
                

        clock_angle = linspace(3/2*pi,pi/2+a1,24);
        clock_angle = clock_angle(2:end-1);
        
        % segmento basso
        x_sb = linspace(K*x_dom,0,25);
        y_sb = -r*ones(size(x_sb));
        
        x_sb = x_sb(2:end); y_sb = y_sb(2:end);
        
        x_rac = r*cos(clock_angle(1:end-1));
        y_rac = r*sin(clock_angle(1:end-1));
        
        x_sa = linspace(r*cos(clock_angle(end)),K*x_dom,25);
        
        % actual y_wake 
        ayw = r*sin(clock_angle(end)) + (K*x_dom + r*sin(a1))*sin(a1);
        y_sa = linspace(r*sin(clock_angle(end)),ayw,25);
        
        x_sa = x_sa(1:end-1); y_sa = y_sa(1:end-1);
        
        PXY = [x_sb,x_rac,x_sa;...
               y_sb,y_rac,y_sa];
              
        
        for w = 1:size(PXY,2)
                       
                lastpoint = lastpoint +1;
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,PXY(1,w),...
                    y_c+PXY(2,w),l_w);
                fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
            
        end


        lastpoint_new = lastpoint;
    case 'wake2'
        n_ot = 150;
        x_te = method_par{1};
        y_te = method_par{2};
        
        alpha = method_par{3};
        r     = method_par{4}; % raggio bolla
        l_w   = method_par{5};
        
        ell_y = 0.7;
        K = 2/3;
        Coef = 20;
        
        
        
        if alpha >= 0 
            y_c = 0.5;
        else
            y_c = -0.5;
        end
        
        delta_te = ((K*x_dom)-x_te)*tand(alpha);

        clock_angle = linspace(3/2*pi,pi/2,24);
        clock_angle = clock_angle(2:end-1);
        
        % segmento basso
        x_sb = linspace(K*x_dom,x_te,n_ot);
        y_sb = linspace(delta_te-r/3,y_c-ell_y*r,n_ot);
          
        x_sb = x_sb(2:end); y_sb = y_sb(2:end);
        
        % raccordo
        x_rac = x_te+r*cos(clock_angle(1:end-1));
        y_rac = y_c+r*ell_y*sin(clock_angle(1:end-1));
        
        
        x_sa = linspace(x_te,K*x_dom,n_ot);
        y_sa = linspace(y_c+ell_y*r,delta_te+r/3,n_ot);
        
        % actual y_wake 
        x_w = linspace(x_te,K*x_dom,4*n_ot);
        y_w = linspace(y_te,delta_te+y_te,4*n_ot);
        
        x_w = x_w(2:end);
        y_w = y_w(2:end);
        
        
        x_sa = x_sa(1:end-1); y_sa = y_sa(1:end-1);
        
        PXY = [x_sb,x_rac,x_sa;...
               y_sb,y_rac,y_sa;];
              
       % figure(99); plot(PXY(1,:),PXY(2,:)); axis equal
       %    hold on; plot(x_w,y_w);
           
        for w = 1:size(PXY,2)
                       
                lastpoint = lastpoint +1;
                
                if PXY(1,w) <= x_te
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,PXY(1,w),...
                    PXY(2,w),l_w);
                else
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,PXY(1,w),...
                    PXY(2,w),l_w+spline([x_te,K*x_dom],[0 Coef],PXY(1,w)));
                    %PXY(2,w),l_w+spline([x_te,K*x_dom],[0 Coef*l_w],PXY(1,w)));
                    
                end
                fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
                
        end
        for w = 1:size(x_w,2)
            
            lastpoint = lastpoint +1;
            fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,x_w(w),...
                y_w(w),l_w/10);
            fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
            
        end
        
%             fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,x_w(1),...
%                 y_w(1),l_w/10);
%             lastpoint = lastpoint +1;
%             fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,x_w(end),...
%                 y_w(end),l_w/10);
%             
%             fprintf(fid,'Line(10000) = {%d,%d};\n',lastpoint-1,lastpoint);
%             fprintf(fid,'Line(10000) In Surface{201};\n',lastpoint);
       

        lastpoint = lastpoint +1;
        lastpoint_new = lastpoint;      
    case 'sublinearwake'
        debug = 0;
        
        y_c = method_par{1};
        y_w = method_par{2};
        r   = method_par{3};
        l_w = method_par{4};
        ellx = 1;%method_par{1};
        elly = 1;%method_par{2};
        % raggio interpolazione lineare
        rref  = method_par{5}; 
        % coefficiente tc nuova lunghezza a rref sia lref*l_linear
        lref  = method_par{6};
        % numero circonferenze da marchiare
        nstaz = method_par{7};
        
        K = 2/3;
        if r > y_w 
            error('r > y_w poco senso per alpha positiva\n')
        end
        
        a1 = (5*pi/180)+atan2(y_w-r,K*x_dom);
                
        clock_angle = linspace(3/2*pi,pi/2+a1,24);
        clock_angle = clock_angle(2:end-1);
        
        % segmento basso
        x_sb = linspace(K*x_dom,0,25);
        y_sb = -r*ones(size(x_sb));
        
        x_sb = x_sb(2:end); y_sb = y_sb(2:end);
        
        x_rac = r*cos(clock_angle(1:end-1));
        y_rac = r*sin(clock_angle(1:end-1));
        
        x_sa = linspace(r*cos(clock_angle(end)),K*x_dom,25);
        
        % actual y_wake 
        ayw = r*sin(clock_angle(end)) + (K*x_dom + r*sin(a1))*sin(a1);
        y_sa = linspace(r*sin(clock_angle(end)),ayw,25);
        
        x_sa = x_sa(1:end-1); y_sa = y_sa(1:end-1);
        
        PXY = [x_sb,x_rac,x_sa;...
               y_sb,y_rac,y_sa];
              
        
        for w = 1:size(PXY,2)
                       
                lastpoint = lastpoint +1;
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,PXY(1,w),...
                    y_c+PXY(2,w),l_w);
                fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
                if debug == 1
                    figure(123)
                    plot(PXY(1,w), y_c+PXY(2,w),'ks');
                    hold on
                end
        end
        
        
        %% sublinear 
                % parametri per forma ellisse

        
        clock_angle = linspace(a1,2*pi+a1-(5*pi/180),24);
        clock_angle = clock_angle(1:end-1);
        
        l_linear = spline([0.75 x_dom],[l_airfoil l_dom],rref);
        
        r_inter = linspace(0.75,x_dom,nstaz); 
        r_inter = r_inter(2:end-1);
        
        l_inter = spline([0.75 rref x_dom],[l_airfoil lref*l_linear l_dom],r_inter);
                
        spacing = figure(22); plot(r_inter,l_inter,'o',[0.75 x_dom],[l_airfoil l_dom],'--',rref,l_linear,'*');title('Size spacing');grid on;
        %savefig(spacing,strcat(case_dir,'/spacing.fig'));
        pause(1); close 22
        
        for w = 1:23            
            for ww = 1:max(size(l_inter))
                x_p_ref = 0.5+ellx*r_inter(ww)*cos(clock_angle(w));
                y_p_ref = y_c+elly*r_inter(ww)*sin(clock_angle(w));
                
                % do fastidio alla scia?
                if  (x_p_ref > 0 && x_p_ref < max(PXY(1,:))) == 1
                    % potrebbe
                    if interp1(x_sb,y_sb+y_c,x_p_ref) > y_p_ref || interp1(x_sa,y_sa+y_c,x_p_ref) < y_p_ref
                        
                        lastpoint = lastpoint +1;
                        fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,x_p_ref,...
                            y_p_ref,l_inter(ww));
                        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
                        
                        if debug == 1
                            figure(123)
                            plot(x_p_ref,y_p_ref,'go');
                            hold on
                        end
                        
                    else
                        if debug == 1
                            figure(123)
                            plot(x_p_ref,y_p_ref,'ro');
                            hold on
                        end
                    end
                else
                    
                    lastpoint = lastpoint +1;
                    fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,x_p_ref,...
                        y_p_ref,l_inter(ww));
                    fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
                    
                    if debug == 1
                        figure(123)
                        plot(x_p_ref,y_p_ref,'go');
                        hold on
                    end
                end
            end

        end
          lastpoint_new = lastpoint;
    otherwise
        
        error('Unknown refining method.')
end














