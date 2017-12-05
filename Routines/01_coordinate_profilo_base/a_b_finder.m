function [ upeli ] = a_b_finder( up )

rpt = 0.005;

[~,iru] = min(abs(up(:,1) - rpt));

up(:,1) = -up(:,1);

global x y m 
x = up(iru,1);
y = up(iru,2);
m = (up(iru+1,2) - up(iru,2))/(up(iru+1,1) - up(iru,1));

OPTIONS = optimoptions('fsolve','Display','off');
ab = fsolve(@a_b_eq,[rpt,y],OPTIONS);

xv = linspace(ab(1)+x,ab(1),50);
xv = xv';


a2 = ab(1)^2; b2 = ab(2).^2;

y_el = @(x,a2,b2) sqrt(b2.*(1-x.^2./a2));
yv = y_el(xv,a2,b2);
xv = -xv+max(xv);


% close all
% plot(xv,yv,'k*');
% hold on
% plot(-x,y,'rx')
% plot(-up(iru+1:end,1),up(iru+1:end,2),'o--')

upeli = [flipud([xv,yv]); [-up(iru+1:end,1),up(iru+1:end,2)] ] ;


clear x y m

end

