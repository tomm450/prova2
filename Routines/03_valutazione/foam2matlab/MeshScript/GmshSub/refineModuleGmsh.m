function [done] = refineModuleGmsh(method,method_par,fid,lastpoint,l_airfoil,x_dom)

switch lower(method)
    case {'none'}
        
        done = 1;
        
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
        
        w = 24;
        fprintf(fid,'Point(%d) = { 1.001,  0.0000000,  0,  %f};\n',lastpoint+w,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+w);
        fprintf(fid,'Point(%d) = { -0.001,  0.0000000,  0,  %f};\n',lastpoint+w+1,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+w+1);
        
        done = 1;
        
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
        savefig(spacing,strcat(case_dir,'/spacing.fig'));
        pause(1); close 22
        
        for w = 1:23
            
            for ww = 1:max(size(l_inter))
                lastpoint = lastpoint +1;
                fprintf(fid,'Point(%d) = { %f,  0.0000000,  %f,  %f};\n',lastpoint,0.5+ellx*r_inter(ww)*cos(clock_angle(w)),...
                    elly*r_inter(ww)*sin(clock_angle(w)),l_inter(ww));
                fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint);
            end
            
        end
        
        w = 24;
        fprintf(fid,'Point(%d) = { 1.001,  0.0000000,  0,  %f};\n',lastpoint+w,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+w);
        fprintf(fid,'Point(%d) = { -0.001,  0.0000000,  0,  %f};\n',lastpoint+w+1,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+w+1);
        done = 1;
        
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
        if r > y_w 
            error('r > y_w poco senso per alpha positiva\n')
        end
        
        a1 = atan2(y_w-r,K*x_dom);
                

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

        fprintf(fid,'Point(%d) = { 1.001,  0.0000000,  0,  %f};\n',lastpoint+1,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+1);
        fprintf(fid,'Point(%d) = { -0.001,  0.0000000,  0,  %f};\n',lastpoint+1+1,l_airfoil/5);
        fprintf(fid,'Point{%d} In Surface{201};\n',lastpoint+1+1);
        done = 1;
        
    otherwise
        
        error('Unknown refining method.')
end














