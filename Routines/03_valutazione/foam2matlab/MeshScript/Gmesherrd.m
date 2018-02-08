function [telapsed] = Gmesherrd(X_IN,IN,STL,GM_par,case_dir,Parameters,SOLVER,SCALA)

% genera mesh con gmsh, aggoinge i BL con refineLayer, estrude con
% extrudeMesh e copia nella cartella simple
if nargin == 7
SCALA = 1;
end

l_airfoil = SCALA*GM_par.l_airfoil;
l_slat    = SCALA*GM_par.l_slat;
l_dom     = SCALA*GM_par.l_dom;
x_dom     = SCALA*GM_par.x_dom;
exp_ratio = GM_par.expRatio;
ds1       = SCALA*GM_par.ds1;


    % Calcolo layer
    thick_cell_try = ds1*exp_ratio.^([0:50]);
    
    if GM_par.BL == 1 % native
        
        % BL non ha limiti come nel caso 2, transizione tc
        % thick_cell(end)*exp_ratio = l_airfoil
        
        [thick_error,i_thick] = min(abs(thick_cell_try-l_airfoil/exp_ratio));
        
        GM_par.nlay = round(i_thick/2)+1;
        
        thick_sum = sum(thick_cell_try(1:i_thick));
        
    elseif GM_par.BL == 2 % foam
        
        % BL si sviluppa dividendo prima cella
        
        for j = 1:max(size(thick_cell_try))
            thick_sum_try(j) = sum(thick_cell_try(1:j));
        end
        
        [thick_error,i_thick] = min(abs(thick_sum_try-l_airfoil));
        
        GM_par.nlay = i_thick;
    
    elseif GM_par.BL == 0 % no     
        
        l_slat    = l_slat/l_airfoil*ds1;
        l_airfoil = ds1;
        
    end
    


FAKE_struct = GM_par.Fstruct;
QUAD        = GM_par.Fquad;
tstart = tic;

%% GEOM
[GEO] = GEOM_phase2(X_IN,IN,Parameters,SCALA);

%% FILE .geo

l_airfoil_str = 'l_a';
l_slat_str    = 'l_s';
%l_flap_str    = 'l_f';


l_dom_str     = 'l_d';
x_dom_str     = 'x_d';


%% SCRIVO FILE

fid = fopen(sprintf('%s/0msh/tom.geo',case_dir),'w+');
fprintf(fid,'// MADE BY TOM\n');

fprintf(fid,'%s = %f; \n',l_slat_str,l_slat);
fprintf(fid,'%s = %f; \n',l_airfoil_str,l_airfoil);
fprintf(fid,'%s = %f; \n',l_dom_str,l_dom);
fprintf(fid,'%s = %f; \n',x_dom_str,x_dom);

%% POINT
p_index = 0;
airfoil_boundary = '';
surf_string      = '';

if GM_par.wtd == 0 || GM_par.wtd == 1 || GM_par.wtd == 3 || GM_par.wtd == 4
    
    
    if GM_par.wtd == 0
        % profilo cavo
        btc = GEO.MAIN{1};
    else
        btc = GEO.MAIN{2};
    end
    
    [p_index,bl_str] = Pointer_Spliner(btc,p_index,2,10,fid,l_airfoil);
    
    airfoil_boundary = strcat(airfoil_boundary,bl_str);
    surf_string      = strcat(surf_string,num2str(2));
end

% slat
if GM_par.wtd >= 2 %|| GM_par.wtd == 3
   
    for j = 1:size(GEO.SLAT,2)
        
        [p_index,bl_str] = Pointer_Spliner(GEO.SLAT{j},p_index,(90+10*j),(90+10*j),fid,l_slat);
        
        airfoil_boundary = strcat(airfoil_boundary,',',bl_str);
        
        if max(size(surf_string)) == 0
            surf_string      = strcat(surf_string,num2str((090+10*j)));
        else
            surf_string      = strcat(surf_string,',',num2str((090+10*j)));
        end
    end
end

% flap
if GM_par.wtd == 4
    
    for j = 1:size(GEO.FLAP,2)
        
        [p_index,bl_str] = Pointer_Spliner(GEO.FLAP{j},p_index,(190+10*j),(190+10*j),fid,l_slat);
        
        airfoil_boundary = strcat(airfoil_boundary,',',bl_str);
        if max(size(surf_string)) == 0
            surf_string      = strcat(surf_string,num2str((190+10*j)));
        else
            surf_string      = strcat(surf_string,',',num2str((190+10*j)));
        end
    end
end


% farfield
fprintf(fid,'Point(%d) = { %s,  0.0000000,  %s,  %s};\n',p_index(end)+1,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { %s,  0.0000000, -%s, %s};\n', p_index(end)+2,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { -%s, 0.0000000, -%s,%s};\n',  p_index(end)+3,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { -%s, 0.0000000,  %s, %s};\n', p_index(end)+4,x_dom_str,x_dom_str,l_dom_str);

lastpoint = p_index(end)+4;


%% LINE
%//Define bounding box edges
fprintf(fid,'Line(1) = {%d, %d};\n',p_index(end)+1,p_index(end)+2);
fprintf(fid,'Line(2) = {%d, %d};\n',p_index(end)+2,p_index(end)+3);
fprintf(fid,'Line(3) = {%d, %d};\n',p_index(end)+3,p_index(end)+4);
fprintf(fid,'Line(4) = {%d, %d};\n',p_index(end)+4,p_index(end)+1);

%//Define foil spline and trailing edge
%//Define bounding box outer boundary
fprintf(fid,'Line Loop(1) = {1, 2, 3, 4};\n');

%//Define unstructured far field mesh zone
fprintf(fid,'Plane Surface(201) = {1,%s};\n',surf_string);

if FAKE_struct == 1
    fprintf(fid,' Transfinite Surface{201}={%s};\n',surf_string);
    %    // forces later meshing to contain structured triangles
    %    // e.g. Transfinite Surface{6} = {1,2,3,4};
    if QUAD  == 1
        
        fprintf(fid,'    Recombine Surface{201};\n');
        %    //combine triangles to quadrangles
        %    // e.g.Recombine Surface{6};
    end
    
end

%% REFINE MODULE

for i = 1:max(size(GM_par.ref_method))
    % casi attualmente implementati
    %     case {'none'}
    %         method_par = {};
    %         done = 1;
    %     case {'clock_simple'}
    %         % 24 punti su circonferenza di raggio rref di dimensione lref
    %         rref = method_par{1};
    %         lref = method_par{2};
    %     case 'sublinear'
    %
    %         % parametri per forma ellisse
    %         ellx = method_par{1};
    %         elly = method_par{2};
    %
    %         % raggio interpolazione lineare
    %         rref  = method_par{3};
    %         % coefficiente tc nuova lunghezza a rref sia lref*l_linear
    %         lref  = method_par{4};
    %         % numero circonferenze da marchiare
    %         nstaz = method_par{5};
    
    [lastpoint] = refineModuleGmsh(GM_par.ref_method{i},GM_par.par_method{i},...
        fid,lastpoint,l_airfoil,x_dom);
    
end

%//Extrude unstructured far field mesh
fprintf(fid,'Extrude {0, 1, 0} {\n');
fprintf(fid,'Surface {201};\n');
fprintf(fid,'Layers{1};\n');
fprintf(fid,'Recombine;\n');
fprintf(fid,'}\n');

if GM_par.wtd <= 2 % una sola superficie
    %//Define physical surfaces - numeric designations from GUI
    fprintf(fid,'Physical Surface("back") = {268};\n');
    fprintf(fid,'Physical Surface("front") = {201};\n');
    fprintf(fid,'Physical Surface("inlet") = {227, 223};\n');
    fprintf(fid,'Physical Surface("outlet") = {219, 231};\n');
    fprintf(fid,'Physical Surface("airfoil") = {251,255,259,263,267,235,239,243,247};\n');
    fprintf(fid,'Physical Volume("internal") = {1};\n');
elseif GM_par.wtd == 3
    %//Define physical surfaces - numeric designations from GUI
    fprintf(fid,'Physical Surface("back") = {313};\n');
    fprintf(fid,'Physical Surface("front") = {201};\n');
    fprintf(fid,'Physical Surface("inlet") = {236, 232};\n');
    fprintf(fid,'Physical Surface("outlet") = {228, 240};\n');
    fprintf(fid,'Physical Surface("airfoil") = {260,264,268,272,276,244,248,252,256,296,300,304,284,280,312,308,288,292};\n');
    fprintf(fid,'Physical Volume("internal") = {1};\n');
    elseif GM_par.wtd == 4
    %//Define physical surfaces - numeric designations from GUI
    fprintf(fid,'Physical Surface("back") = {365};\n');
    fprintf(fid,'Physical Surface("front") = {201};\n');
    fprintf(fid,'Physical Surface("inlet") = {252, 248};\n');
    fprintf(fid,'Physical Surface("outlet") = {244, 256};\n');
    fprintf(fid,'Physical Surface("airfoil") = {312,316,320,324,328,296,300,304,308,276,280,272,268,284,288,264,292,260,348,352,344,356,340,336,360,364,332};\n');
    fprintf(fid,'Physical Volume("internal") = {1};\n');
end

if GM_par.BL == 1
    %hwall * ratio^(dist/hwall)
    fprintf(fid,'Field[1] = BoundaryLayer;\n');
    
    %AnisoMax
    %Threshold angle for creating a mesh fan in the boundary layer
    %type: float
    %default value: 10000000000
    fprintf(fid,'Field[1].AnisoMax = 1000;\n');
    
    %EdgesList
    %Indices of curves in the geometric model for which a boundary layer
    %is needed
    %type: list
    %default value: {}
    fprintf(fid,'Field[1].EdgesList = {%s};\n',bl_str);
    
    %FanNodesList
    %Indices of vertices in the geometric model for which a fan is created
    %type: list
    %default value: {}
    
    %IntersectMetrics
    %Intersect metrics of all faces
    %type: integer
    %default value: 0
    fprintf(fid,'Field[1].IntersectMetrics = 0;\n');
    
    %NodesList
    %Indices of vertices in the geometric model for which a BL ends
    %type: list
    %default value: {}
    
    
    %Quads Generate recombined elements in the boundary layer
    %type: integer
    %default value: 0
    fprintf(fid,'Field[1].Quads = 1;\n');
    
    %hfar Element size far from the wall
    %type: float
    %default value: 1
    fprintf(fid,'Field[1].hfar = %f;\n',       thick_cell_try(i_thick));
    
    %hwall_n Mesh Size Normal to the The Wall
    %type: float
    %default value: 0.1
    fprintf(fid,'Field[1].hwall_n = %f;\n',    thick_cell_try(1));
    
    %hwall_n_nodes
    %Mesh Size Normal to the The Wall at nodes (overwrite hwall n
    %when defined)
    %type: list double
    %default value: {}
    
    %ratio
    %Size Ratio Between Two Successive Layers
    %type: float
    %default value: 1.1
    fprintf(fid,'Field[1].ratio = %f;\n',      exp_ratio);
    
    %thickness
    %Maximal thickness of the boundary layer
    %type: float
    %default value: 0.01
    fprintf(fid,'Field[1].thickness = %f;\n',  thick_sum);
    fprintf(fid,'BoundaryLayer Field = 1; \n');
    
end

fclose(fid);

%% ACTUAL MESHING
system(sprintf('touch %s/logmesh.txt',case_dir));

% mesh
%if Parameters.n_processori == 1
 [status] = goGoGmsh(sprintf('-3 -o %s/0msh/tom.msh %s/0msh/tom.geo > %s/logmesh.txt',case_dir,case_dir,case_dir),...
     Parameters.gmsh_cmd);
% else
% [status] = goGoGmsh(sprintf('-3 -o %s/0msh/tom.msh %s/0msh/tom.geo > %s/logmesh.txt',case_dir,case_dir,case_dir),...
%     Parameters.gmsh_cmd,sprintf('mpirun -n %d ',Parameters.n_processori));
% end    

% copio e converto
[status] = system(sprintf('rm %s/15dummy/*.msh >> %s/logmesh.txt',case_dir,case_dir));

[status] = system(sprintf('cp %s/0msh/tom.msh %s/15dummy >> %s/logmesh.txt',case_dir,case_dir,case_dir));
[status] = goGoOpenFOAM(sprintf('cd %s/15dummy &&  gmshToFoam ./tom.msh >> ../logmesh.txt',case_dir));

% sistemo BC
[done] = BC_cor('def',case_dir);

[status] = system(sprintf('rm  %s/20extrude/constant/polyMesh/* >> %s/logmesh.txt',case_dir,case_dir));

[status] = system(sprintf('rm -r %s/30simple/constant/polyMesh >> %s/logmesh.txt',case_dir,case_dir));
[status] = system(sprintf('rm -r %s/40piso/constant/polyMesh >> %s/logmesh.txt',case_dir,case_dir));

[status] = goGoOpenFOAM(sprintf('cd %s/20extrude && extrudeMesh >> ../logmesh.txt',case_dir));% && paraFoam');

if GM_par.BL == 2
    
    for j = 1:GM_par.nlay-1
        
        step = GM_par.nlay-j;
        thick_start = sum( ds1*exp_ratio.^([0:step]));
        thick_end   = sum( ds1*exp_ratio.^([0:step-1]));
        
        [status] = goGoOpenFOAM(sprintf('cd %s/20extrude && refineWallLayer -overwrite ''(airfoil)'' %1.3f >> ../logmesh.txt',...
            case_dir,thick_end/thick_start));
        
    end
end

[resultRenumber,logRenumber]=goGoOpenFOAM(sprintf('cd %s/20extrude && renumberMesh',case_dir));

if strcmp(SOLVER.solver,'simple')
    [~,~] = system(sprintf('cp -r %s/20extrude/1/polyMesh %s/30simple/constant/',case_dir,case_dir));
elseif strcmp(SOLVER.solver,'piso')
    [~,~] = system(sprintf('cp -r %s/20extrude/1/polyMesh %s/40piso/constant/',case_dir,case_dir));
else
    error('someway,somewhere,someone fuck up');
end

[pos_win,pos_move] = system(sprintf('cp *.pos %s/',case_dir));

telapsed = toc(tstart);

end
%% SCRIVO PUNTI DI UNA DATA GEOMETRIA
function [p_index,bl_str] = Pointer_Spliner(btc,p_index,n_loop,n_spline,fid,l_cell)

leri = [];

for k = 1:size(btc{1},1)
    
    if btc{1}(k,1) > btc{2}(1)+0.1*(btc{3}(1)-btc{2}(1)) &&...
            btc{1}(k,1) < btc{2}(1)+0.9*(btc{3}(1)-btc{2}(1))
        
        fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %s};\n',p_index(end)+k,btc{1}(k,1),btc{1}(k,2),l_cell);
        
        if size(leri,2) == 0
            leri = [leri,p_index(end)+k];
        elseif size(leri,2) == 2
            leri = [leri,p_index(end)+k];
        end
        
    elseif btc{1}(k,1) >= btc{2}(1)+0.9*(btc{3}(1)-btc{2}(1))
        
        fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %d};\n',p_index(end)+k,btc{1}(k,1),btc{1}(k,2),l_cell/3);
        if size(leri,2) == 3
            leri = [leri,p_index(end)+k];
        end
        
    elseif btc{1}(k,1) <= btc{2}(1)+0.1*(btc{3}(1)-btc{2}(1))
        
        fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %d};\n',p_index(end)+k,btc{1}(k,1),btc{1}(k,2),l_cell/3);
        if size(leri,2) == 1
            leri = [leri,p_index(end)+k];
        end
        
    end
    
end

if max(size(leri)) == 4 %
    % ok
else
    error('leri ha dimensioninon congure\n');
end

fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline,  p_index(end)+1,leri(1)-3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+1,  leri(1)-3,leri(1)+3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+2,  leri(1)+3,leri(2)-3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+3,  leri(2)-3,leri(2)+3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+4,  leri(2)+3,leri(3)-3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+5,  leri(3)-3,leri(3)+3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+6,  leri(3)+3,leri(4)-3);
fprintf(fid,'Spline(%d) = {%d:%d};\n',   n_spline+7,  leri(4)-3,leri(4)+3);
fprintf(fid,'Spline(%d) = {%d:%d,%d};\n',n_spline+8,  leri(4)+3,p_index(end)+k,p_index(end)+1);

bl_str = sprintf('%d,%d,%d,%d,%d,%d,%d,%d,%d',n_spline,...
    n_spline+1,n_spline+2,n_spline+3,n_spline+4,...
    n_spline+5,n_spline+6,n_spline+7,n_spline+8);

fprintf(fid,'Line Loop(%d) = {%s};\n',n_loop,bl_str);

p_index = [p_index,p_index(end)+k];

end

%% CORREZIONE BC DOPO CONVERSIONE
function [done] = BC_cor(input,case_dir)

if input == 'def'

fb_cor  = sprintf('%s\n%s\n','        type            empty;',...
    '        inGroups        1(empty);');

a_cor   = sprintf('%s\n%s\n','        type            wall;',...
    '        inGroups        1(wall);');

clear rewrite

[~,numBc] = system(sprintf('sed ''18!d''  %s/15dummy/constant/polyMesh/boundary',case_dir));

if str2num(numBc) == 5 % ok
else 
    str2num(numBc)
    error('ho troppi boundary')
end


[~,intro] = system(sprintf('head -n 19  %s/15dummy/constant/polyMesh/boundary',case_dir));
[~,outro] = system(sprintf('tail -n 3  %s/15dummy/constant/polyMesh/boundary',case_dir));

[win,body] = system(sprintf('tail -n +20 %s/15dummy/constant/polyMesh/boundary > temp.txt && head -n 35 temp.txt',case_dir));
[~,~] = system('rm temp.txt');
body_frag = strsplit(body,'\n'); %body_frag = body_frag{1:end-1};



% back front airfoil
fid = fopen(sprintf('%s/15dummy/constant/polyMesh/boundary',case_dir),'w+');
fprintf(fid,'%s\n',intro);
i_vect = [1,8,15,22,29];

for k = 1:5
    j = i_vect(k);
    name2comp = strtrim(body_frag{j});
    if strcmp(name2comp,'front') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},fb_cor,body_frag{j+4},body_frag{j+5},body_frag{j+6});
    elseif strcmp(name2comp,'back') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},fb_cor,body_frag{j+4},body_frag{j+5},body_frag{j+6});
    elseif strcmp(name2comp,'airfoil') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},a_cor,body_frag{j+4},body_frag{j+5},body_frag{j+6});
    else
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},body_frag{j+2},body_frag{j+3},body_frag{j+4},body_frag{j+5},body_frag{j+6});
    end
    fprintf(fid,'%s\n',rewrite{end});
end


fprintf(fid,'%s\n',outro);
fclose(fid);



done = win;

else 
    done = 0;
end


end





