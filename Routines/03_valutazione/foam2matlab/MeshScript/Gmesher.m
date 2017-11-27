function [telapsed] = Gmesher(X_IN,IN,STL,GM_par,case_dir,Parameters,SOLVER)

% genera mesh con gmsh, aggoinge i BL con refineLayer, estrude con
% extrudeMesh e copia nella cartella simple

l_airfoil = GM_par.l_airfoil;
l_slat    = GM_par.l_slat;
l_dom     = GM_par.l_dom;
x_dom     = GM_par.x_dom;
exp_ratio = GM_par.expRatio;
ds1       = GM_par.ds1;


% rref = GM_par.rref;
% lref = GM_par.lref;
%point_txt = STL.point_txt;

FAKE_struct = GM_par.Fstruct;
QUAD        = GM_par.Fquad;

% my_dir = what;
% my_dir = my_dir.path;
% my_dirs = genpath(my_dir);
% addpath(my_dirs);

tstart = tic;

%% GEOM
[GEOM] = GEOM_phase2(X_IN,IN,Parameters);

%% FILE .geo

l_airfoil_str = 'l_a';
l_slat_str    = 'l_s';
l_dom_str     = 'l_d';
x_dom_str     = 'x_d';

if GM_par.wtd == 0
    pu = Parameters.Airfoil.up(round(linspace(1,size(Parameters.Airfoil.up,1),10)),:);
    pd = Parameters.Airfoil.dwn(round(linspace(1,size(Parameters.Airfoil.dwn,1),10)),:);
    p = [flipud(pu(2:end,:));pd(1:end-1,:)]/1000;
    %p = [p(1:5:end,:)];
elseif GM_par.wtd == 1
    p = [flipud(GEOM.up_land(2:end,:));GEOM.dwn_land(1:end-1,:)]/1000;
    %p = [p(1:5:end,:)];
elseif GM_par.wtd == 2
    psu = flipud(GEOM.slat_land_u)/1000; psu = psu(round(linspace(1,size(psu,1),150)),:);
    psd = GEOM.slat_land_d/1000;         psd = psd(round(linspace(1,size(psd,1),150)),:);
    
    ps = [psu(1:end-1,:);psd(1:end-1,:)];
elseif GM_par.wtd == 3
    p = [flipud(GEOM.up_land(2:end,:));GEOM.dwn_land(1:end-1,:)]/1000;
    %p = [p(1:5:end,:)];
    
    psu = flipud(GEOM.slat_land_u)/1000; psu = psu(round(linspace(1,size(psu,1),150)),:);
    psd = GEOM.slat_land_d/1000;         psd = psd(round(linspace(1,size(psd,1),150)),:);
    
    ps = [psu(1:end-1,:);psd(1:end-1,:)];
end


fid = fopen(sprintf('%s/0msh/tom.geo',case_dir),'w+');
fprintf(fid,'// MADE BY TOM\n');

fprintf(fid,'%s = %f; \n',l_slat_str,l_slat);
fprintf(fid,'%s = %f; \n',l_airfoil_str,l_airfoil);
fprintf(fid,'%s = %f; \n',l_dom_str,l_dom);
fprintf(fid,'%s = %f; \n',x_dom_str,x_dom);

if GM_par.wtd == 0 || GM_par.wtd == 1 || GM_par.wtd == 3

  for k = 1:size(p,1)
      if p(k,1) < 0.1
        fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %f};\n',k,p(k,1),p(k,2),l_airfoil/5);
      else
        fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %s};\n',k,p(k,1),p(k,2),l_airfoil_str);
      end
  end
  pend1 = k;
else
  pend1 = 0;
end

% curva superiore slat
if GM_par.wtd == 2 || GM_par.wtd == 3
  for k = 1:size(ps,1)-1
     fprintf(fid,'Point(%d) = { %f, 0.0000000, %f, %s};\n',pend1+k,ps(k,1),ps(k,2),l_slat_str);
  end
  pend = pend1+k;
else
  pend = pend1;
end

% farfield
fprintf(fid,'Point(%d) = { %s,  0.0000000,  %s,  %s};\n',pend+1,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { %s,  0.0000000, -%s, %s};\n',pend+2,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { -%s, 0.0000000, -%s,%s};\n',pend+3,x_dom_str,x_dom_str,l_dom_str);
fprintf(fid,'Point(%d) = { -%s, 0.0000000,  %s, %s};\n',pend+4,x_dom_str,x_dom_str,l_dom_str);

lastpoint = pend+4;

%//Define bounding box edges
fprintf(fid,'Line(1) = {%d, %d};\n',pend+1,pend+2);
fprintf(fid,'Line(2) = {%d, %d};\n',pend+2,pend+3);
fprintf(fid,'Line(3) = {%d, %d};\n',pend+3,pend+4);
fprintf(fid,'Line(4) = {%d, %d};\n',pend+4,pend+1);

%//Define foil spline and trailing edge
%//Define bounding box outer boundary
fprintf(fid,'Line Loop(101) = {1, 2, 3, 4};\n');


if GM_par.wtd == 0 || GM_par.wtd == 1 
 
 fprintf(fid,'Spline(5) = {1:%d,1};\n',size(p,1));
 fprintf(fid,'Line Loop(102) = {5};\n');

 surf_string ='101, 102';

elseif GM_par.wtd == 2 

 fprintf(fid,'Spline(6) = {%d:%d, %d};\n',pend1+1,pend,pend1+1);%,pend1+1);
 fprintf(fid,'Line Loop(103) = {6};\n');
 surf_string ='101, 103';

elseif GM_par.wtd == 3 
    
 fprintf(fid,'Spline(5) = {1:%d,1};\n',size(p,1));
 fprintf(fid,'Spline(6) = {%d:%d, %d};\n',pend1+1,pend,pend1+1);%,pend1+1);
  
 fprintf(fid,'Line Loop(102) = {5};\n');
 fprintf(fid,'Line Loop(103) = {6};\n');

 surf_string ='101, 102, 103';
 
end

%//Define unstructured far field mesh zone
fprintf(fid,'Plane Surface(201) = {%s};\n',surf_string);

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

[done] = refineModuleGmsh(GM_par.ref_method,GM_par.par_method,...
    fid,lastpoint,l_airfoil,x_dom);

%//Extrude unstructured far field mesh
fprintf(fid,'Extrude {0, 1, 0} {\n');
fprintf(fid,'Surface {201};\n');
fprintf(fid,'Layers{1};\n');
fprintf(fid,'Recombine;\n');
fprintf(fid,'}\n');

if GM_par.wtd <= 2
    %//Define physical surfaces - numeric designations from GUI
    fprintf(fid,'Physical Surface("back") = {228};\n');
    fprintf(fid,'Physical Surface("front") = {201};\n');
    fprintf(fid,'Physical Surface("inlet") = {219, 215};\n');
    fprintf(fid,'Physical Surface("outlet") = {211, 223};\n');
    fprintf(fid,'Physical Surface("airfoil") = {227};\n');
    fprintf(fid,'Physical Volume("internal") = {1};\n');
elseif GM_par.wtd == 3
    %//Define physical surfaces - numeric designations from GUI
    fprintf(fid,'Physical Surface("back") = {233};\n');
    fprintf(fid,'Physical Surface("front") = {201};\n');
    fprintf(fid,'Physical Surface("inlet") = {220, 216};\n');
    fprintf(fid,'Physical Surface("outlet") = {212, 224};\n');
    fprintf(fid,'Physical Surface("airfoil") = {228, 232};\n');
    fprintf(fid,'Physical Volume("internal") = {1};\n');
end

fclose(fid);

%% ACTUAL MESHING
system(sprintf('touch %s/logmesh.txt',case_dir));

% mesh
%if Parameters.n_processori == 1
% [status] = goGoGmsh(sprintf('-3 -o %s/0msh/tom.msh %s/0msh/tom.geo > %s/logmesh.txt',case_dir,case_dir,case_dir),...
%     Parameters.gmsh_cmd);
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

% % aggiungo layer
thick_cell_try = ds1*exp_ratio.^([0:50]);
for j = 1:max(size(thick_cell_try))
    thick_sum_try(j) = sum(thick_cell_try(1:j));
end

[thick_error,i_thick] = min(abs(thick_sum_try-l_airfoil));

GM_par.nlay = i_thick;

for j = 1:GM_par.nlay-1
    
    step = GM_par.nlay-j;
    thick_start = sum( ds1*exp_ratio.^([0:step]));
    thick_end   = sum( ds1*exp_ratio.^([0:step-1]));
   
    [status] = goGoOpenFOAM(sprintf('cd %s/20extrude && refineWallLayer -overwrite ''(airfoil)'' %1.3f >> ../logmesh.txt',...
        case_dir,thick_end/thick_start));
    
end


if strcmp(SOLVER.solver,'simple')
    [~,~] = system(sprintf('cp -r %s/20extrude/constant/polyMesh %s/30simple/constant/',case_dir,case_dir));
elseif strcmp(SOLVER.solver,'piso')
    [~,~] = system(sprintf('cp -r %s/20extrude/constant/polyMesh %s/40piso/constant/',case_dir,case_dir));
else
    error('someway,somewhere,someone fuck up');
end



telapsed = toc(tstart);

end

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
[dummy,superdummy] = system('rm temp.txt');
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






