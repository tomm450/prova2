function [GEO] = GEOM_phase2(X_IN,IN,Parameters,SCALA) %,nome_out,FLAG_raccordo,x_cut,wtd)

[ GEOM,fail,log ] = geometry_main(Parameters,IN);

if size(X_IN,2) == 0
    % deorga
    fprintf('X_IN == [];\nCarico file ./deroga.mat contenente le cordinate assulute dei (fino a) 3 elementi\n');
    
%     load('./deroga.mat')
%     GEOM.up_land   = airfoil_up;
%     GEOM.dwn_land  = airfoil_dwn;
%     
%     GEOM.slat_land_u = slat_up;
%     GEOM.slat_land_d = slat_dwn;
%     
%     GEOM.flap_land_u = flap_up;
%     GEOM.flap_land_d = flap_dwn;
      load('./7airfoil.mat');
    
else
    
    %nome_out = 'tronco';
    
    
    %[ GEOM,ERR] = slat_position( GEOM, X_IN(1), X_IN(2),X_IN(3),1);
    
    [ GEOM,~] = slat_position_x_y_a( GEOM, X_IN(1), X_IN(2),X_IN(3));
    %GEOM.up_land  = GEOM.up_land/1000;
    %GEOM.dwn_land = GEOM.dwn_land/1000;
    
end

n_camp = 150;

%% PROFILO ISOLATO
if max(size(Parameters.Airfoil.up)) >= n_camp
    pu = Parameters.Airfoil.up(round(linspace(1,size(Parameters.Airfoil.up,1),n_camp)),:);
else
    pu = Parameters.Airfoil.up;
end

if max(size(Parameters.Airfoil.dwn)) >= n_camp
    pd = Parameters.Airfoil.dwn(round(linspace(1,size(Parameters.Airfoil.dwn,1),n_camp)),:);
else
    pd = Parameters.Airfoil.dwn;
end


p0  = SCALA*[flipud(pu(2:end,:));pd(1:end-1,:)]/1000;
LE0 = SCALA*pd(1,:)/1000;
TE0 = SCALA*pd(end,:)/1000;
c0  = norm(TE0-LE0);  

GEO.MAIN{1} = {p0,LE0,TE0,c0};  

%% PROFILO MENO SUPERFICI MOBILI

if max(size(GEOM.up_land)) >= n_camp
    pu = GEOM.up_land(round(linspace(1,size(GEOM.up_land,1),n_camp)),:)/1000;
else
    pu = GEOM.up_land/1000;
end

if max(size(GEOM.up_land)) >= n_camp
    pd = GEOM.dwn_land(round(linspace(1,size(GEOM.dwn_land,1),n_camp)),:)/1000;
else
    pd = GEOM.dwn_land/1000;
end

p1 =  SCALA*[flipud(pu(2:end,:));pd(1:end-1,:)];
LE1 = SCALA*pd(1,:);
TE1 = SCALA*pd(end,:);
c1  = norm(TE1-LE1);  

GEO.MAIN{2} = {p1,LE1,TE1,c1};

%% SLAT
if isfield(GEOM,'slat_land_u')
    if iscell(GEOM.slat_land_u)
        
        kend = size(GEOM.slat_land_u,2);
    
    else
        
        kend = 1;
        temp = {GEOM.slat_land_u,GEOM.slat_land_d};
        GEOM.slat_land_u = {temp{1}};
        GEOM.slat_land_d = {temp{2}};
        
        
        
    end
    
    for k = 1:kend
        psu = flipud(GEOM.slat_land_u{k})/1000;
        if max(size(psu)) >= n_camp
            psu = psu(round(linspace(1,size(psu,1),n_camp)),:);
        end
        
        psd = GEOM.slat_land_d{k}/1000;
        if max(size(psd)) >= n_camp
            psd = psd(round(linspace(1,size(psd,1),n_camp)),:);
        end
        
        ps  = SCALA*[psu(1:end-1,:);psd(1:end-1,:)];
        LEs = SCALA*psd(1,:);
        TEs = SCALA*psd(end,:);
        cs  = norm(TEs-LEs);
        
        GEO.SLAT{k} = {ps,LEs,TEs,cs};
    end
end

%% FLAP
if isfield(GEOM,'flap_land_u')
    if iscell(GEOM.flap_land_u)
        
        kend = size(GEOM.flap_land_u,2);
        
    else
        
        kend = 1;
        temp = {GEOM.flap_land_u,GEOM.flap_land_d};
        GEOM.flap_land_u = {temp{1}};
        GEOM.flap_land_d = {temp{2}};
    end
    
    for k = 1:kend
        pfu = flipud(GEOM.flap_land_u{k})/1000;
        if max(size(pfu)) >= n_camp
            pfu = pfu(round(linspace(1,size(pfu,1),n_camp)),:);
        end
        
        pfd = GEOM.flap_land_d{k}/1000;
        if max(size(pfd)) >= n_camp
            pfd = pfd(round(linspace(1,size(psd,1),n_camp)),:);
        end
        
        pf =  SCALA*[pfu(1:end-1,:);pfd(1:end-1,:)];
        LEf = SCALA*pfd(1,:);
        TEf = SCALA*pfd(end,:);
        cf  = norm(TEf-LEf);
        
        GEO.FLAP{k} = {pf,LEf,TEf,cf};
    end
end
    


end