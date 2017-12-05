function [done] = komega(k_inlet,omega_inlet,case_dir,SOLVER,BU_par)


if strcmp(SOLVER.solver,'simple')
   fid = fopen(sprintf('%s/30simple/0/omega',case_dir),'w+');
else
   fid = fopen(sprintf('%s/40piso/0/omega',case_dir),'w+');
end

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\  \n');
fprintf(fid,'| =========                 |                                                 |  \n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n');
fprintf(fid,'|  \\    /   O peration     | Version:  5                                     |  \n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |  \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 |  \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/  \n');
fprintf(fid,'FoamFile  \n');
fprintf(fid,'{  \n');
fprintf(fid,'    version     2.0;  \n');
fprintf(fid,'    format      ascii;  \n');
fprintf(fid,'    class       volScalarField;  \n');
fprintf(fid,'    object      omega;  \n');
fprintf(fid,'}  \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //  \n');
  
fprintf(fid,'dimensions      [0 0 -1 0 0 0 0];  \n');
fprintf(fid,'internalField   uniform %f;  \n',omega_inlet);

fprintf(fid,'boundaryField  \n');
fprintf(fid,'{  \n');

fprintf(fid,'    inlet  \n');
fprintf(fid,'    {    \n');
fprintf(fid,'        type            fixedValue;  \n');
fprintf(fid,'        value           uniform %f;  \n',omega_inlet);
fprintf(fid,'    }  \n');

fprintf(fid,'    outlet  \n');
fprintf(fid,'    {  \n');
fprintf(fid,'        type            zeroGradient;  \n');
fprintf(fid,'    }  \n');


fprintf(fid,'    airfoil  \n');
fprintf(fid,'    {  \n');
if BU_par.wall_function == 1
    fprintf(fid,'        type            omegaWallFunction;  \n');
    fprintf(fid,'        Cmu              0.09;  \n');
    fprintf(fid,'        kappa            0.41;  \n');
    fprintf(fid,'        E                9.8;  \n');
    fprintf(fid,'        beta1            0.075;  \n');
    fprintf(fid,'        value           uniform 9.8;  \n');
else
    fprintf(fid,'        type            fixedValue;  \n');
    fprintf(fid,'        value           uniform %3.9f;  \n',BU_par.omega_body);
end
fprintf(fid,'    }  \n');


fprintf(fid,'        front  \n');
fprintf(fid,'    {  \n');
fprintf(fid,'        type            empty;  \n');
fprintf(fid,'    }  \n');
fprintf(fid,'    back  \n');
fprintf(fid,'    {  \n');
fprintf(fid,'        type            empty;  \n');
fprintf(fid,'    }  \n');
fprintf(fid,'}  \n');


fclose(fid);

%%

%fid = fopen(sprintf('%s/30simple/0/k',case_dir),'w+');

if strcmp(SOLVER.solver,'simple')
   fid = fopen(sprintf('%s/30simple/0/k',case_dir),'w+');
else
   fid = fopen(sprintf('%s/40piso/0/k',case_dir),'w+');
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
fprintf(fid,'    class       volScalarField; \n');
fprintf(fid,'    object      k; \n');
fprintf(fid,'} \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');

fprintf(fid,'dimensions      [0 2 -2 0 0 0 0]; \n');
fprintf(fid,'internalField   uniform %f; \n',k_inlet);

fprintf(fid,'boundaryField \n');
fprintf(fid,'{ \n');

fprintf(fid,'    inlet \n');
fprintf(fid,'    {   \n');
fprintf(fid,'        type            fixedValue; \n');
fprintf(fid,'        value           uniform %f; \n',k_inlet);
fprintf(fid,'    } \n');

fprintf(fid,'    outlet \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            zeroGradient; \n');
fprintf(fid,'    } \n');

fprintf(fid,'    airfoil \n');
fprintf(fid,'    { \n');
if BU_par.wall_function == 1
    fprintf(fid,'        type            kqRWallFunction; \n');
    fprintf(fid,'        value           uniform 0.0096; \n');
else
    fprintf(fid,'        type            fixedValue;  \n');
    fprintf(fid,'        value           uniform 0;  \n');
end
fprintf(fid,'    } \n');

fprintf(fid,'        front \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            empty; \n');
fprintf(fid,'    } \n');
fprintf(fid,'    back \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            empty; \n');
fprintf(fid,'    } \n');
fprintf(fid,'} \n');

fclose(fid);

done = 1;

end
