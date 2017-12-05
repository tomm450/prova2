function [CL,CD] = aeroCoeff( xc,dl,theta_G,cp,BU_par )

n = size(xc,2)/2;

cp_air  = cp(1:2*n);
cp_slat = cp(2*n+1:end);
% 
% punti_ventre_air = [xc(1,1:n)',    yc(1,1:n)'];
% punti_dorso_air =  [xc(1,n+1:end)',yc(1,n+1:end)'];
% 
% punti_ventre_slat = [xc(2,1:n)',    yc(2,1:n)'];
% punti_dorso_slat  = [xc(2,n+1:end)',yc(2,n+1:end)'];
% 
% x_p_ventre_air = [xc(1,1:n)',cp(1:n)];
% x_p_dorso_air  = [xc(1,n+1:end)',cp(n+1:2*n)];
% 
% x_p_ventre_slat = [xc(2,1:n)',cp(2*n+1:3*n)];
% x_p_dorso_slat  = [xc(2,n+1:end)',cp(3*n+1:end)];


CD_ref_air = sum(  dl(1,:).*cp_air'.*sin(theta_G(1,:)) + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
CL_ref_air = sum( -dl(1,:).*cp_air'.*cos(theta_G(1,:)) - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));

CL = CL_ref_air*cosd(BU_par.alpha) - CD_ref_air*sind(BU_par.alpha);
CD = CL_ref_air*sind(BU_par.alpha) + CD_ref_air*cosd(BU_par.alpha);

end

