function [ done,cronos] = Snapper( X_IN,IN,STL,MESH_par,case_dir,Parameters )


[dummy,dummy] = system(sprintf('mkdir %s/10snappy',case_dir));
[dummy,dummy] = system(sprintf('cp -r %s/30simple/system %s/10snappy',case_dir,case_dir));
[dummy,dummy] = system(sprintf('cp -r %s/30simple/constant %s/10snappy',case_dir,case_dir));
[dummy,dummy] = system(sprintf('mkdir %s/10snappy/constant/triSurface',case_dir));


if MESH_par.simplegrading.flag == 1
    nucleo_costante = MESH_par.simplegrading.nucleo_costante*MESH_par.xMax;
    cella_base = MESH_par.simplegrading.cella_base;
    
    nsg_coeff = round(MESH_par.nx/2-nucleo_costante/cella_base);
    n_cost    = MESH_par.nx - nsg_coeff;
    
    
    [MESH_par.sg_coeff,nx_cor] = simple_grading_helper(nsg_coeff,cella_base,MESH_par.xMax-nucleo_costante,MESH_par.simplegrading.sg_lim);
    
    MESH_par.nx = n_cost+2*nx_cor;
    MESH_par.nz = MESH_par.nx;
    % calcolo percentuali
    % direzioni
    MESH_par.pd1 = 100*(MESH_par.xMax-nucleo_costante)/(2*MESH_par.xMax);
    MESH_par.pd2 = 100*(nucleo_costante)/(MESH_par.xMax);
    % celle
    MESH_par.pn1 = 100*(nx_cor)/MESH_par.nx;
    MESH_par.pn2 = 100*(2*nucleo_costante/cella_base)/MESH_par.nx;
end

%% scrivo SnappyHex
if MESH_par.wakeBox == 1
    MESH_par.zWake = (MESH_par.xMax-L)*tand(BU_par.alpha);
    MESH_par.xWake = MESH_par.xMax;
end

% calcolo
minDx = ((MESH_par.xMax-MESH_par.xMin)/MESH_par.nx)/...
    (2.^(MESH_par.MaxrefFactor+MESH_par.deltabox-1));
minDz = ((MESH_par.zMax-MESH_par.zMin)/MESH_par.nz)/...
    (2.^(MESH_par.MaxrefFactor+MESH_par.deltabox-1));

if MESH_par.simplegrading.flag == 0
    Dmean = 0.5*(minDx+minDz);
else
    Dmean = cella_base/(2.^(MESH_par.MaxrefFactor));
end

MESH_par.endLay = 1/MESH_par.expRatio;
MESH_par.startLay = MESH_par.ds1/Dmean;

[err,nlay] = min(abs(Dmean/MESH_par.ds1 - MESH_par.expRatio.^[1:250]));
MESH_par.imin = nlay-1;



%% STL slat
tstart = tic;

if MESH_par.wtd == 0

suffisso = '_airfoil_def';

elseif MESH_par.wtd == 1

suffisso = '_airfoil';

elseif MESH_par.wtd == 2

suffisso = '_slat';

elseif MESH_par.wtd == 3

suffisso = '';

end


[x_TE] = stl_process(X_IN,IN,STL.point_txt,STL.point_txt(1:end-4),0,nan,MESH_par.wtd,case_dir);
%system(sprintf('cp %s.stl %s',strcat('./0stl/out/',STL.nome_out,suffisso),'./10snappy/constant/triSurface'));
crono(1) = toc(tstart)/60


%% MESH
tstart = tic;
%             % uso blockMesh
[done] = blockWrite(MESH_par,case_dir);


goGoOpenFOAM(sprintf('cd %s/10snappy && blockMesh > ../1logblock.txt',case_dir));
crono(2) = toc(tstart)/60


tstart = tic;
%[done] = snappyWriteDual(STL.nome_out,suffisso,MESH_par,x_TE,{'true','true','true'});
[done] = snappyWrite(STL.point_txt(1:end-4),suffisso,MESH_par,x_TE,120,case_dir);

system(sprintf('mkdir %s/10snappy/0',case_dir));

goGoOpenFOAM(sprintf('cd %s/10snappy && decomposePar -force >../1logsnapdec.txt && foamJob -s -p renumberMesh >../1logsnapren.txt && foamJob -s -p snappyHexMesh  > ../1logbsnap.txt && reconstructParMesh -latestTime -mergeTol 1e-6 >../1logsnaprec.txt',case_dir));
%goGoOpenFOAM('./runsnappy.sh');
[dummy,snapReport] = system(sprintf('tail -n 42 %s/1logbsnap.txt',case_dir));
snapReport
crono(3) = toc(tstart)/60

[~,dirTree] = system(sprintf('ls -1 %s/10snappy/',case_dir));
dirTree = strsplit(dirTree,'\n');
%SimpledirTree = SimpledirTree{1:end-1};
Ddir = []; 
for k = 1:size(dirTree,2)-1
    if isnan(str2num(dirTree{k})) == 0
        if str2num(dirTree{k}) ~= 0
            Ddir = [Ddir,str2num(dirTree{k})];
                        
        end
    end
end



system(sprintf('rm -r %s/15dummy/constant/polyMesh',case_dir));
system(sprintf('cp -r %s/10snappy/%d/* %s/15dummy/constant/',case_dir,max(Ddir),case_dir));

[done] = BC_cor('def',case_dir);

tstart = tic;
goGoOpenFOAM(sprintf('cd %s/20extrude && extrudeMesh > ../2logextr.txt',case_dir));

% system('rm -r ./30simple/constant/polyMesh');
% system('rm -r ./40piso/constant/polyMesh');
system(sprintf('cp -r %s/20extrude/constant/polyMesh %s/30simple/constant/',case_dir,case_dir));

crono(4) = toc(tstart)/60;

done = 1;
cronos = sum(crono);
end

function result=goGoOpenFOAM(command)
result=unix(['export LD_LIBRARY_PATH=""; . /opt/openfoam5/etc/bashrc; ' command]);
end

function [sg_coeff,n_out] = simple_grading_helper(n,d_cell,tgt,sg_lim)
n_out = n;
a = zeros(1,n);
a(1) = d_cell;

a_def = a;

%tgt = 20;

simplegradK = 1:0.0001:3;

for j = 1:max(size(simplegradK))
    for w = 1:max(size(a))-1
        a(w+1) = a(w)*simplegradK(j);
    end
    
    somma(j) = sum(a);
end

errore = abs(tgt-somma);

[~,imin] = min(errore);


for w = 1:max(size(a))-1
    a_def(w+1) = a_def(w)*simplegradK(imin);
end


sg_coeff = a_def(1)/a_def(end);

n_out = n;
while sg_coeff < sg_lim
    %fprintf('correggo numero di celle\n')
    n_out = round(1.1*n_out);
    
    a = zeros(1,n_out);
    a(1) = d_cell;
    
    a_def = a;
    
    %tgt = 20;
    
    simplegradK = 1:0.0001:3;
    
    for j = 1:max(size(simplegradK))
        for w = 1:max(size(a))-1
            a(w+1) = a(w)*simplegradK(j);
        end
        
        somma(j) = sum(a);
    end
    
    errore = abs(tgt-somma);
    
    [~,imin] = min(errore);
    
    
    for w = 1:max(size(a))-1
        a_def(w+1) = a_def(w)*simplegradK(imin);
    end
    
    sg_coeff = a_def(1)/a_def(end);
    
    
end

if n_out <= n || n_out > 1.5*n
    n_out = n;
    a = zeros(1,n);
    a(1) = d_cell;
    
    a_def = a;
    
    %tgt = 20;
    
    simplegradK = 1:0.0001:3;
    
    for j = 1:max(size(simplegradK))
        for w = 1:max(size(a))-1
            a(w+1) = a(w)*simplegradK(j);
        end
        
        somma(j) = sum(a);
    end
    
    errore = abs(tgt-somma);
    
    [~,imin] = min(errore);
    
    
    for w = 1:max(size(a))-1
        a_def(w+1) = a_def(w)*simplegradK(imin);
    end
    
    
    sg_coeff = a_def(1)/a_def(end);
    
end


end


function [done] = BC_cor(input,case_dir)

if input == 'def'

fb_cor  = sprintf('%s\n%s\n','        type            empty;',...
    '        inGroups        1(empty);');

a_cor   = sprintf('%s\n%s\n','        type            wall;');

clear rewrite

[~,numBc] = system(sprintf('sed ''18!d''  %s/15dummy/constant/polyMesh/boundary',case_dir));

% if str2num(numBc) == 5 % ok
% else 
%     str2num(numBc)
%     error('ho troppi boundary')
% end


[~,intro] = system(sprintf('head -n 19  %s/15dummy/constant/polyMesh/boundary',case_dir));
[~,outro] = system(sprintf('tail -n 3  %s/15dummy/constant/polyMesh/boundary',case_dir));

[win,body] = system(sprintf('tail -n +20 %s/15dummy/constant/polyMesh/boundary > temp.txt && head -n 31 temp.txt',case_dir));
[dummy,superdummy] = system('rm temp.txt');
body_frag = strsplit(body,'\n'); %body_frag = body_frag{1:end-1};



% back front airfoil
fid = fopen(sprintf('%s/15dummy/constant/polyMesh/boundary',case_dir),'w+');
fprintf(fid,'%s\n',intro);
i_vect = [1,7,13,19,25];

for k = 1:5
    
    j = i_vect(k);
    name2comp = strtrim(body_frag{j});
    if strcmp(name2comp,'front') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},fb_cor,body_frag{j+3},body_frag{j+4},body_frag{j+5});
    elseif strcmp(name2comp,'back') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},fb_cor,body_frag{j+3},body_frag{j+4},body_frag{j+5});
    elseif strcmp(name2comp,'airfoil') == 1
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n}\n',body_frag{j},body_frag{j+1},a_cor,body_frag{j+3},body_frag{j+4},body_frag{j+5});
    else
        rewrite{k} = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n',body_frag{j},body_frag{j+1},body_frag{j+2},body_frag{j+3},body_frag{j+4},body_frag{j+5});
    end
    fprintf(fid,'%s\n',rewrite{end});
    %fprintf('%s\n',rewrite{end});
    
end


fprintf(fid,'%s\n',outro);
fclose(fid);



done = win;

else 
    done = 0;
end


end

