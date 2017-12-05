function [ F ] = a_b_eq( ab )

global x y m 
F(1) = -(ab(2).^2/ab(1)^2)*(ab(1) + x)/y - m;
F(2) = (x+ab(1))^2/ab(1)^2 + y^2/ab(2)^2 -1;

end

