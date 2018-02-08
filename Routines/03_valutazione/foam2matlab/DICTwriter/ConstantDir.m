function [done] = ConstantDir(case_dir,SOLVER,BU_par)


%% TRANSPORT PROPERTIES
if strcmp(SOLVER.solver,'simple')
   fid = fopen(sprintf('%s/30simple/constant/transportProperties',case_dir),'w+');
else
   fid = fopen(sprintf('%s/40piso/constant/transportProperties',case_dir),'w+');
end

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
fprintf(fid,'| =========                 |                                                 | \n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n');
fprintf(fid,'|  \\    /   O peration     | Version:  5                                     | \n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 | \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/ \n');
fprintf(fid,'FoamFile \n');
fprintf(fid,'{ \n');
fprintf(fid,'    version     2.0; \n');
fprintf(fid,'    format      ascii; \n');
fprintf(fid,'    class       dictionary; \n');
fprintf(fid,'    location    "constant"; \n');
fprintf(fid,'    object      transportProperties; \n');
fprintf(fid,'} \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');

fprintf(fid,'transportModel  Newtonian; \n');

fprintf(fid,'rho             [1 -3 0 0 0 0 0] %f; \n',BU_par.Rho);

fprintf(fid,'nu              [0 2 -1 0 0 0 0] %d; \n',BU_par.Nu);

fprintf(fid,'// ************************************************************************* // \n');

fclose(fid);

%done = 0.5;


%% TURBOLENT PROPERTIES
if strcmp(SOLVER.solver,'simple')
   fid = fopen(sprintf('%s/30simple/constant/turbulenceProperties',case_dir),'w+');
else
   fid = fopen(sprintf('%s/40piso/constant/turbulenceProperties',case_dir),'w+');
end

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
fprintf(fid,'| =========                 |                                                 | \n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n');
fprintf(fid,'|  \\    /   O peration     | Version:  5                                     | \n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 | \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/ \n');

fprintf(fid,'FoamFile \n');
fprintf(fid,'{ \n');
fprintf(fid,'    version     2.0; \n');
fprintf(fid,'    format      ascii; \n');
fprintf(fid,'    class       dictionary; \n');
fprintf(fid,'    location    "constant"; \n');
fprintf(fid,'    object      turbulenceProperties; \n');
fprintf(fid,'} \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');

if strcmp(SOLVER.T_model,'la')
    
    fprintf(fid,'simulationType laminar;  \n');
    
else
    fprintf(fid,'simulationType RAS;  \n');
    
    fprintf(fid,'RAS \n');
    fprintf(fid,'{ \n');
    
    if strcmp(SOLVER.T_model,'sa')
        fprintf(fid,'   RASModel        SpalartAllmaras; \n');
    elseif strcmp(SOLVER.T_model,'ko')
        fprintf(fid,'   RASModel        kOmegaSST;  \n');
    end
    
    fprintf(fid,'   turbulence      on; \n');
        
    fprintf(fid,'   printCoeffs     on; \n');
    fprintf(fid,'} \n');
    
end
fprintf(fid,'// ************************************************************************* // \n');

fclose(fid);

done = 1;



end