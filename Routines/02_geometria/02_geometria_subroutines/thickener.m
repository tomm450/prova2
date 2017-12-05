function [v_new,v_old_add] = thickener(v_old,x_ref,dx,seed)
%%

% ho una matrice v_old composta da due colonne
% v_old(:,1) = componenti x
% v_old(:,2) = componenti y
% voglio infittire attorno a un x_ref (ed eventualmente inserire x_ref nella
% matrice) con un numero di punti pari a 2*seed+1 distanti toll



%% SEED in ogni caso sottost� all intervallo tra x_ref-1 e x_ref+1

% identifico punto più vicino al nuovo x_ref
[res_i_start_u,i_start_u] = min(abs(v_old(:,1) - x_ref));

if res_i_start_u ~= 0
    %% sistemo matrici
    if v_old(i_start_u,1) < x_ref % ho stima per difetto
        lim_left  = i_start_u;
        lim_right = i_start_u+1;
    else                          % ho stima per eccesso
        lim_left  = i_start_u-1;
        lim_right = i_start_u;
    end
    
    v_old_left = v_old(1:lim_left,:);      % matrice punti alla sinistra del punto
    v_old_right = v_old(lim_right:end,:);  % matrice punti alla destra del punto
    
    % controllo che seed sia realizzabile
    
    if x_ref -dx*seed < v_old_left(end,1)
        % metterei punti con x < x_left(end)
        % devo adattare intervallo o numero di punti
        
        % ho deciso di modificare il dx
        dx_mod(1) =  ((x_ref - v_old_left(end,1)) / (seed +1) );
        
        % ho deciso di modificare numero di punti
        % seed_mod(1) = floor((x_ref - v_old_left(end,1)) / dx );
    else
        dx_mod(1) = dx;
        %seed_mod(1) = seed;
    end
    
    if x_ref +dx*seed > v_old_right(1,1)
        % metterei punti con x > x_right(1)
        % devo adattare intervallo o numero di punti
        
        % ho deciso di modificare il dx
        dx_mod(2) = ((v_old_right(1,1) - x_ref) / (seed +1) ) ;
        
        % ho deciso di modificare numero di punti
        %seed_mod(2) = floor((v_old_right(1,1) - x_ref) / dx ) ;
    else
        dx_mod(2) = dx;
        %seed_mod(2) = seed;
    end
    
    % seed/dx sistemato non dovrebbe dare problemi
%      if min(dx_mod) == dx
%     else
%         fprintf('adattato il dx da %f a %f\n',dx,min(dx_mod));
%     end
    dx = min(dx_mod);
    %seed = min(seed_mod);
    
    %% inizio vero e proprio
    if seed == 0 % aggiungo solo punto verifica (toll diventa inutile)
        
        if i_start_u >= 41 && size(v_old,1)-i_start_u >= 41
            % posso permettermi di considerare vettore curve di 80 elementi
            % per l'interpolazione spline
            delta = 40;
            
        else
            % adatto il delta
            [delta,~] = min([i_start_u-2,size(v_old,1)-i_start_u-2]);
        end
        
        % add su asse x
        v_old_add_x  = x_ref;
        
        % interpolo
        if delta <= 0
            % overkill?
            azz = 0;
            if sum(isnan(v_old(:,1)))>= 1
                azz = 1;
            end
            if sum(isnan(v_old(:,2)))>= 1
                azz = 1;
            end
            if sum(isnan(v_old_add_x(:,1)))>= 1
                azz = 1;
            end
            
            if azz == 1
               v_old
               error('Nan')
            end
            
            
            v_old_add_y  = spline(v_old(:,1),...
                v_old(:,2),v_old_add_x');
        else
            
            v_old_add_y  = spline(v_old(i_start_u-delta:i_start_u+delta,1),...
                v_old(i_start_u-delta:i_start_u+delta,2),v_old_add_x');
            
        end
        
        %         if solo_point == 1
        %             v_new = [v_old_add_x,v_old_add_y];
        %             return
        %         end
        
        
    else
        
        
        if x_ref +dx*seed > v_old_right(1,1)
            error('riga 48 thick')
        end
        if x_ref -dx*seed < v_old_left(end,1)
            error('riga 51 thick')
        end
        
        v_old_add_x  = x_ref + dx*[-seed:seed]; % NON FIDARSI DI MATLAB,LE PARENTESI SERVONO
        
        v_old_add_y  = spline(v_old(i_start_u-seed:i_start_u+seed,1),...
            v_old(i_start_u-seed:i_start_u+seed,2),v_old_add_x');
        
        
    end
    
    v_old_add = [v_old_add_x',v_old_add_y];
    
    
    v_new = [v_old_left; v_old_add; v_old_right];
    
else
    v_new = v_old;
end