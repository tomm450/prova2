function [x_TE] = stl_process(X_IN,IN,point_txt,nome_out,FLAG_raccordo,x_cut,wtd,case_dir)

%nome_out = 'tronco';


%%
%Airfoil = load('NACA64212at.txt');
Airfoil = load(point_txt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT UTENTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% costruzione geometria
C_mesh        = 1000; %[mm];
Parameters.Geom.d_geometry    = 5000;
Parameters.Geom.toll_geometry = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolo dati grezzi e applico approssimazione elisse
% (out escono unitari)
[ up, dwn, Airfoil ] = interfoil(  Airfoil,Parameters.Geom.d_geometry,1 );

% aggiorno Parameters
Parameters.Airfoil.up  = C_mesh*up;
Parameters.Airfoil.dwn = C_mesh*dwn;
Parameters.Airfoil.c   = C_mesh;

Parameters.Airfoil.Mom_ref_x = 0.25*C_mesh;
Parameters.Airfoil.Mom_ref_y = 0;

%Parameters.Dir_path.main = my_dir;

[ GEOM,fail,log ] = geometry_main(Parameters,IN);

[ GEOM,ERR] = slat_position_x_y_a( GEOM, X_IN(1), X_IN(2),X_IN(3));

k     = 1000;
passo = 1;

%% RACCORDO?
x_TE = 1;
if FLAG_raccordo == 1
    [up,dwn] = raccordo( up,dwn,x_cut );
    x_TE = max(up(:,1));
end

 sy = 5;
 Y = linspace(-1,1,sy)';
 body_land_dummy = 1/k.*[GEOM.up_land(1:end-1,:);flipud(GEOM.dwn_land)];
 C = max(body_land_dummy(:,1));


if wtd == 0
   % % %% PROFILO
   body_land = [up(1:end-1,:);flipud(dwn)];
   %body_land = [body_land(1:4600,:);body_land(5401:end,:)];
   i_r = ceil(size(body_land,1)/6);
   %figure; plot(body_land(:,1),body_land(:,2),'o-')
   index =[1:(i_r-1),i_r:passo:(size(body_land,1)-i_r),(size(body_land,1)-i_r+1):size(body_land,1)];
   X = body_land(index,1);
   % %
   Y = Y*C;
   Z = repmat(body_land(index,2)',sy,1);
   h = surface(X,Y,Z);
   % % % close all
   figure(7777777);
   surface(X,Y,Z,Z)
   close 7777777
   % %
   fvc = surf2patch(h,'triangles');
   stl_writer(fvc, sprintf('%s/10snappy/constant/triSurface/%s_airfoil_def',case_dir,nome_out));

elseif wtd == 1
    body_land = 1/k.*[GEOM.up_land(1:end-1,:);flipud(GEOM.dwn_land)];
    i_r = ceil(size(body_land,1)/6);
    
    %figure; plot(body_land(:,1),body_land(:,2),'o-')
    
    index =[1:(i_r-1),i_r:passo:(size(body_land,1)-i_r),(size(body_land,1)-i_r+1):size(body_land,1)];
    X = body_land(index,1);
    
    sy = 5;
    %sy = round(1.3*max(size(X)));
    
    
    
    Y = Y*C;
    Z = repmat(body_land(index,2)',sy,1);
    figure(7777777);
    h = surface(X,Y,Z);
    
    % close all
    % surface(X,Y,Z,Z)
    
    fvc = surf2patch(h,'triangles');
    stl_writer(fvc, sprintf('%s/10snappy/constant/triSurface/%s_airfoil',case_dir,nome_out));
    close 7777777
    %system(sprintf('cp ./0stl/out/%s.stl ./10snapppy/constant/triSurface/',strcat(nome_out,'_airfoil')));
elseif wtd == 2
    %% SLAT
    body_lands = 1/k.*[GEOM.slat_land_u(1:end-1,:);flipud(GEOM.slat_land_d)];
    i_r = ceil(size(body_lands,1)/6);
    
    %figure; plot(body_land(:,1),body_land(:,2),'o-')
    
    indexs =[1:(i_r-1),i_r:passo:(size(body_lands,1)-i_r),(size(body_lands,1)-i_r+1):size(body_lands,1)];
    Xs = body_lands(indexs,1);
    
    %sy = round(1.3*max(size(X)));
    
    Ys = linspace(-1,1,sy)';
    %C = max(body_lands(:,1));
    Ys = Y*C;
    Zs = repmat(body_lands(indexs,2)',sy,1);
    
    figure(7777777);
    hs = surface(Xs,Ys,Zs,'visible','off');
    
    
    % close all
    % surface(X,Y,Z,Z)
    
    fvcs = surf2patch(hs,'triangles');
    %stl_writer(fvcs, strcat('./0stl/out/',nome_out,'_slat'))
    stl_writer(fvcs, sprintf('%s/10snappy/constant/triSurface/%s_slat',case_dir,nome_out));
    close 7777777
else
    body_land = 1/k.*[GEOM.up_land(1:end-1,:);flipud(GEOM.dwn_land)];
    i_r = ceil(size(body_land,1)/6);
    
    %figure; plot(body_land(:,1),body_land(:,2),'o-')
    
    index =[1:(i_r-1),i_r:passo:(size(body_land,1)-i_r),(size(body_land,1)-i_r+1):size(body_land,1)];
    X = body_land(index,1);
    
    sy = 5;
    %sy = round(1.3*max(size(X)));
    
    
    
    Y = Y*C;
    Z = repmat(body_land(index,2)',sy,1);
    figure(7777777);
    h = surface(X,Y,Z);
    
    % close all
    % surface(X,Y,Z,Z)
    
    fvc = surf2patch(h,'triangles');
    %stl_writer(fvc, strcat('./0stl/out/',nome_out,'_airfoil'))
    stl_writer(fvc, sprintf('%s/10snappy/constant/triSurface/%s_airfoil',case_dir,nome_out));
    close 7777777
    
        %% SLAT
    body_lands = 1/k.*[GEOM.slat_land_u(1:end-1,:);flipud(GEOM.slat_land_d)];
    i_r = ceil(size(body_lands,1)/6);
    
    %figure; plot(body_land(:,1),body_land(:,2),'o-')
    
    indexs =[1:(i_r-1),i_r:passo:(size(body_lands,1)-i_r),(size(body_lands,1)-i_r+1):size(body_lands,1)];
    Xs = body_lands(indexs,1);
    
    %sy = round(1.3*max(size(X)));
    
    Ys = linspace(-1,1,sy)';
    %C = max(body_lands(:,1));
    Ys = Y*C;
    Zs = repmat(body_lands(indexs,2)',sy,1);
    
    figure(7777777);
    hs = surface(Xs,Ys,Zs,'visible','off');
    
    
    % close all
    % surface(X,Y,Z,Z)
    
    fvcs = surf2patch(hs,'triangles');
    %stl_writer(fvcs, strcat('./0stl/out/',nome_out,'_slat'))
    stl_writer(fvc, sprintf('%s/10snappy/constant/triSurface/%s_slat',case_dir,nome_out));
    close 7777777
    %system(sprintf('cp ./0stl/out/%s.stl ./10snapppy/constant/triSurface/',strcat(nome_out,'_airfoil')));
    
    
nv1=length(fvc.vertices);

fv_combined.vertices=[fvc.vertices;fvcs.vertices];
fv_combined.faces=[fvc.faces; fvcs.faces+nv1];

%stl_writer(fv_combined, strcat('./0stl/out/',nome_out))
stl_writer(fv_combined, sprintf('%s/10snappy/constant/triSurface/%s',case_dir,nome_out));
%system(sprintf('cp ./0stl/out/%s.stl ./10snapppy/constant/triSurface/',strcat(nome_out,'_slat')));

end
%system(sprintf('cp ./0stl/out/%s.stl ./10snapppy/constant/triSurface/',nome_out));

end

%%

function stl_writer(fv, name)
%routine to creat stl file from FV structure array
% fv is face-vertex structure array
% name is the file name (character string)
% either just give the 'root' file name or give 'root.stl'
lnm = length(name);

if name((lnm-3):(lnm-1))=='stl'
    label = name(1:lnm-3);
else
    label = name;
    name = sprintf('%s.stl', name);
end

v1=fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,1),:);
v2=fv.vertices(fv.faces(:,3),:)-fv.vertices(fv.faces(:,2),:);
Norms=cross3(v1,v2);
clear v1 v2
v1(:,1:3)=fv.vertices(fv.faces(:,1),1:3);
v2(:,1:3)=fv.vertices(fv.faces(:,2),1:3);
v3(:,1:3)=fv.vertices(fv.faces(:,3),1:3);
fid = fopen(name,'w');
fprintf(fid,'solid %s\n',label);
nf = length(fv.faces); %k = (1:nf)';
for k = 1:nf
    fprintf(fid,'facet normal %5.5f %5.5f %5.5f\n outer loop\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\n endloop\n endfacet\n',...
        Norms(k,1),Norms(k,2),Norms(k,3), v1(k,1), v1(k,2),v1(k,3),v2(k,1), v2(k,2), v2(k,3),v3(k,1), v3(k,2), v3(k,3) );
end
fprintf(fid,'endsolid %s\n',label);
fclose(fid);
end
function M=cross3(r,F)
% function to calculate normalized cross product rxF/|rxF|
% handles (same-size) arrays (n by 3) for r and F
%
M = [(r(:,2).*F(:,3) - r(:,3).*F(:,2)) ...
    (r(:,3).*F(:,1) - r(:,1).*F(:,3)) ...
    (r(:,1).*F(:,2) - r(:,2).*F(:,1))];
M_mag = sqrt(sum((M.*M)')');
M(:,1) = M(:,1)./M_mag;
M(:,2) = M(:,2)./M_mag;
M(:,3) = M(:,3)./M_mag;
end

function [up2,dwn2] = raccordo( up,dwn,x_cut )
%% raccordo

[err,iup] = min(abs(x_cut-up(:,1)));
[err,idwn] = min(abs(x_cut-dwn(:,1)));

if iup == idwn
    % coordinate x uguali
else
    sdeng
end


O = [up(iup,1),0.5*(dwn(idwn,2)+up(iup,2))];
%angolo up
dup  = atan2d((up(iup,2)   - up(iup-1,2)),(up(iup,1) - up(iup-1,1)));
%angolo dwn
ddwn = atan2d((dwn(idwn,2) - dwn(idwn-1,2)),(dwn(idwn,1) - dwn(idwn-1,1)));

% angolo linea media
cusp = 0.5*(dup+ddwn);
% coordinate relative
pu = (up(iup-1:iup,:) - O)';
pd = (dwn(idwn-1:idwn,:) -O)';
raggio = pu(2,2);

R = [cosd(cusp) sind(cusp); -sind(cusp) cosd(cusp)];

pur = R*pu;
pdr = R*pd;

%
% close all
% axis equal
% plot(pu(1,:),pu(2,:),'g',pd(1,:),pd(2,:),'g',0,0,'rx')
% grid on
% hold on
% plot(pur(1,:),pur(2,:),'bo-',pdr(1,:),pdr(2,:),'bo-',0,0,'rx')
y_int = linspace(pdr(2,1),pur(2,1),100);
x_new = spline([pdr(2,:),pur(2,:)],[pdr(1,:),pur(1,:)],y_int);

% hold on
% plot(x_new,linspace(pdr(2,1),pur(2,1),100));

Rin = [cosd(-cusp) sind(-cusp); -sind(-cusp) cosd(-cusp)];
[new]= Rin*[x_new;y_int];

new(1,:) = new(1,:)+O(1);

new(2,:) = new(2,:)+O(2);

% figure
% plot(up(:,1),up(:,2),'b',dwn(:,1),dwn(:,2),'b');
% hold on
% plot(new(1,:),new(2,:),'r');

[TEnew,imax] = max(x_new);

up2  = [up(1:iup-1,:);flipud([new(1,imax:end);new(2,imax:end)]')];

dwn2 = [dwn(1:iup-1,:);[new(1,1:imax)',new(2,1:imax)']];

% up2(iup-10:iup+10,:)
% 
% dwn2(iup-10:iup+10,:)
% 
% plot(up2(:,1),up2(:,2),'o',dwn2(:,1),dwn2(:,2),'o')



end

