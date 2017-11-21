function [done] = refineModuleGmsh(method,method_par,fid,lastpoint)

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
    
    otherwise
        
        error('Unknown refining method.')
end














