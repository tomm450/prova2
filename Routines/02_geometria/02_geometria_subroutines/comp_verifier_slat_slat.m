function [win] = comp_verifier_slat_slat(slat_up,slat_dwn)

% isolo curva interna (verificare per slat inferiore)
slat_inu =  slat_up(2:2:end,:);
slat_inu = [slat_up(1,:); slat_inu];

slat_ind =  slat_dwn(2:2:end,:);
slat_ind = [slat_dwn(1,:); slat_ind];

% figure(200)
% plot(slat_inu(:,1),slat_inu(:,2),'gx')
% hold on
% plot(slat_ind(:,1),slat_ind(:,2),'bx')


[~,imax_x] = max([slat_inu(end,1), slat_ind(end,1)]);

if imax_x == 1
    long = slat_inu;
    short = slat_ind;
    soprasotto = 1;
else
    long = slat_ind;
    short = slat_inu;
    soprasotto = 2;
end

for j = 2:size(short,2)
    
              %toll serve solo nel caso seed sia diverso da 0...guarda
              % thickener...prima o poi correggo
    [~,longP] = thickener(long,short(j,1),0,0); % genereo pt su curva long
                                                   % a coordinata
                                                   % short(j,1)
    if soprasotto == 1
        if longP(2) < short(j,2)
            win = 0;
            %fprintf('fail @ x=%f y=%f\n',short(j,1),short(j,2))
            return
        end
    end
    
    if soprasotto == 2
        if longP(2) > short(j,2)
            win = 0;
            %fprintf('fail @ x=%f y=%f\n',short(j,1),short(j,2))
            return
        end
    end
end


win = 1;
end
