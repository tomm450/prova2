function [win] = comp_verifier(curve,slat,soprasotto)

% isolo curva interna (verificare per slat inferiore)
slat_in =  slat(2:2:end,:);
slat_in = [slat(1,:); slat_in];

[err,j_piv] = min(abs(curve(:,1) - slat_in(end,1)));
if abs(err) >= 1e-3
    error('siamo fuori toll')
end
% figure(853)
% plot(slat_in(:,1),slat_in(:,2));
% hold on; title('comp very')
for j = j_piv -2 : -1 : 2
    
    toll = 0; %toll serve solo nel caso seed sia diverso da 0...guarda
              % thickener...prima o poi correggo
    %       = thickener(v_old,  x_ref,     dx,   seed)
    [~,v_old_add] = thickener(slat_in,curve(j,1),0,0);
    %plot(v_old_add(1),v_old_add(2),'o')
    
    if soprasotto == 1
        if v_old_add(2) < curve(j,2)
            win = 0;
            %fprintf('fail @ x=%f y=%f',curve(j,1),curve(j,2))
            return
        end
    elseif soprasotto == 2
        if v_old_add(2) > curve(j,2)
            win = 0;
            %fprintf('fail @ x=%f y=%f',curve(j,1),curve(j,2))
            return
        end
    end
end
win = 1;
end
