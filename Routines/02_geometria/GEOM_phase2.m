function [GEOM] = GEOM_phase2(X_IN,IN,Parameters) %,nome_out,FLAG_raccordo,x_cut,wtd)

%nome_out = 'tronco';
[ GEOM,fail,log ] = geometry_main(Parameters,IN);

%[ GEOM,ERR] = slat_position( GEOM, X_IN(1), X_IN(2),X_IN(3),1);

[ GEOM,~] = slat_position_x_y_a( GEOM, X_IN(1), X_IN(2),X_IN(3));