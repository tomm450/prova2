function [up, dwn, out ] = interfoil( in,d,elipse_corr )
%% INTERPOLO CON SPLINE PUNTI PROFILO IN MODO DA AUMENTARNE IL NUMERO
% % input 
% - in -> matrice coordinata in ingresso ([1 0 1]',[up LE dwn]')
% - d - > numero desiderato di divisioni lungo la corda (la distribuzione 
%         sar� logaritmica) [ a prescindere dall'immissione il risultato
%         sar� dispari ]
% - elipse_corr LE => per serie 6 il LE è un ellisse(2 in realtà se
%               profilo asimmetrico) che si raccondano a 0.005c
% % output
% - up  -> matrice coordinate dorso 
% - dwn -> matrice coordinate ventre
% - out -> matrice coordinate in uscita ( dimensione ~ (d*2)-1 x2 )

%% Inizio 
if max(in(:,1)) ~= 1
    error('interfoil: Dati in ingresso devono riferirsi a corda unitaria')
end

if nargin == 2
    elipse_corr = 0;
end

if elipse_corr == 1
    [ in  ] = ellipse_corr( in );
end

m = round((size(in,1) - 1) /2);

up_temp = in(1:m,:);
up_temp(:,1) = -up_temp(:,1); 

dwn_temp = in( (m+1):end,:);

stra = [up_temp;dwn_temp];


% creo nuovo vettore x
x = logspace(0,1.5,ceil(d/2)); % prima met�
x = x-1;
x = x./(2*max(x));
x2 = 1-x;

x = [x,fliplr(x2(1:end-1))]';  % vettore completo
xLETE = x;
s = size(xLETE,1);

x = [-flipud(x);x(2:end)];
clear x2

y = spline(stra(:,1),stra(:,2),x);

% plot(x,y,'o--',stra(:,1),stra(:,2),'x--')


up  = [xLETE, flipud(y(1:s))];
dwn = [xLETE, y(s:end)];

out = [flipud(up);dwn(2:end,:)];
end

