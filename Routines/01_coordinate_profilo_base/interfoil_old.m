function [up, dwn, out ] = interfoil( in,d )
%% INTERPOLO CON SPLINE PUNTI PROFILO IN MODO DA AUMENTARNE IL NUMERO
% % input 
% - in -> matrice coordinata in ingresso ([1 0 1]',[up LE dwn]')
% - d - > numero desiderato di divisioni lungo la corda (la distribuzione 
%         sarà logaritmica) [ a prescindere dall'immissione il risultato
%         sarà dispari ]
% % output
% - up  -> matrice coordinate dorso 
% - dwn -> matrice coordinate ventre
% - out -> matrice coordinate in uscita ( dimensione ~ (d*2)-1 x2 )

%% Inizio 
if max(in(:,1)) ~= 1
    error('interfoil: Dati in ingresso devono riferirsi a corda unitaria')
end

m = (size(in,1) - 1) /2;

% Può comunque essere utile, cambio distribuzione punti
% if m >= d
%     error('interfoil: Numero di punti ingresso maggiore di quello desiderato')
% end 

up_temp = flipud(in(1:m+1,:));
dwn_temp = in( (m+1):end,:);

% creo nuovo vettore x
x = logspace(0,2,ceil(d/2)); % prima metà
x = x-1;
x = x./(2*max(x));
x2 = 1-x;

x = [x,fliplr(x2(1:end-1))]';  % vettore completo
clear x2

up  = [x, spline(up_temp(:,1),up_temp(:,2),x)];
dwn = [x, spline(dwn_temp(:,1),dwn_temp(:,2),x)];

out = [flipud(up);dwn(2:end,:)];
end

