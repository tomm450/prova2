function [ XX,YY ] = airfoil_interpolator_cos(iLE,n_pt,XP,YP)

LE = [XP(iLE),YP(iLE)];
TE = [XP(1),YP(1)];

temp = [(XP-LE(1))'; (YP-LE(2))'];


b_ts = fliplr(temp(:,1:iLE));
t_ts = temp(:,(iLE):end);


%hold on; plot(b_ts(1,:),b_ts(2,:),'c')
%         plot(t_ts(1,:),t_ts(2,:),'m')

c = norm(TE - LE);

theta = atan2(TE(1,2)-LE(1,2),TE(1,1)-LE(1,1));

R = [cos(theta)  sin(theta)
    -sin(theta) cos(theta)];

br = R*b_ts;
tr = R*t_ts;

br = sortrows(br'); br = br';
tr = sortrows(tr'); tr = tr';

X = [-fliplr(br(1,2:end)),tr(1,:)];
Y = [ fliplr(br(2,2:end)),tr(2,:)];

np = round(n_pt/2);

%% X
%% cossspace
% x1 = cosspace(0,c,np);
% %x1 = [x1(1),x1(3:end)];
% x1 = fliplr(x1);
%% linspace
%x1 = linspace(c,0,np);
%% logspace
%       x1t = logspacing(0,c/2,round(np/2),0);
%       x2t = c - fliplr(x1t);
%       x1 = [x1t(1:end-1),x2t];

%x2 = [x1(1:end-1), -fliplr(x1)];
%COx = [ x1(1:end-1), fliplr(x1)];
%% ISO s
sx = linspace(min(X),max(X),100*np);
sy = interp1(X,Y,sx,'pchip');
%plot(sx,sy,'r.--')
dS = sqrt(diff(sx).^2 + diff(sy).^2);
dS_tgt = sum(dS)/(n_pt);
dS_tgt_v = dS_tgt*[0:n_pt-1];


CdS = [0];
for i = 1:numel(dS)
   CdS(i+1) = sum(dS(1:i)); 
end



x2 = interp1(CdS,sx,dS_tgt_v,'pchip'); 
%x2 = [min(X),x2];

%figure(1000); plot(X,Y)
y1 = interp1(X,Y,x2,'pchip');
%y1 = interp1(X,Y,x2,'spline');

COx = x2-2*x2.*(x2<0);
COy = [ y1(1:np-1), (y1(np:end))];

%       %hold on
%       figure(666)
%       plot(COx,COy,'o-'); hold on

theta = -theta;

R = [cos(theta)  sin(theta)
    -sin(theta) cos(theta)];

C = R*[COx;COy];

XX = C(1,:) +LE(1);
YY = C(2,:) +LE(2);

end

function [ x ] = cosspace( a,b,n )
%COSSPACE Summary of this function goes here
%   As linspace but x=0.5*(1-cos(pi*linspace(a,b,n)));

x=0.5*(1-cos(pi*linspace(0,1,n)));
x = a+(b-a)*x;

end

function x = logspacing(a,b,n,dOrder)

t = logspace(0,1+dOrder,n);
t = (t-1);
t = t/max(t);

x = a + (b-a)*t;

end
%end

