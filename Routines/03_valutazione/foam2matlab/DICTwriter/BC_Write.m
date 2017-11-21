function [done] = BC_Write(solver,BU_par,BU_type,case_dir )
% BU_type = freestream o fixedValue

if strcmp(solver,'simple') == 1
  dir = '30simple';
elseif strcmp(solver,'piso') == 1
  dir = '40piso';
else
  error('solver errato')
end

if nargin == 3
    case_dir = '.';
end

%% SCRIVO U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(sprintf('%s/%s/0/U',case_dir,dir),'w+');

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
fprintf(fid,'    class       volVectorField; \n');
fprintf(fid,'    object      U; \n');
fprintf(fid,'} \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
fprintf(fid,'dimensions      [0 1 -1 0 0 0 0]; \n');
fprintf(fid,'internalField   uniform (%3.2f 0 %3.2f); \n',BU_par.Ux,BU_par.Uz);
fprintf(fid,'boundaryField \n');
fprintf(fid,'{ \n');

fprintf(fid,'   outlet \n');
fprintf(fid,'    { \n');
if strcmp(BU_type,'freestream') == 1
  fprintf(fid,'          type            freestream; \n');
  fprintf(fid,'          freestreamValue uniform (%3.2f 0 %3.2f); \n',BU_par.Ux,BU_par.Uz);
elseif strcmp(BU_type,'fixedValue') == 1
  fprintf(fid,'          type            zeroGradient; \n');
else
  error('BU_type errato')
end
fprintf(fid,'       \n');
fprintf(fid,'    } \n');

fprintf(fid,'    inlet \n');
fprintf(fid,'    { \n');
if strcmp(BU_type,'freestream') == 1
  fprintf(fid,'          type            freestream; \n');
  fprintf(fid,'          freestreamValue uniform (%3.2f 0 %3.2f); \n',BU_par.Ux,BU_par.Uz);
elseif strcmp(BU_type,'fixedValue') == 1
  fprintf(fid,'          type            fixedValue; \n');
  fprintf(fid,'          value uniform (%3.2f 0 %3.2f); \n',BU_par.Ux,BU_par.Uz);
else
  error('BU_type errato')
end
fprintf(fid,'         \n');
fprintf(fid,'    } \n');

fprintf(fid,'   front \n');
fprintf(fid,'   { \n');
fprintf(fid,'       type             empty; \n');
fprintf(fid,'    } \n');
fprintf(fid,'    back \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            empty; \n');
fprintf(fid,'    } \n');
fprintf(fid,'    airfoil \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            noSlip; \n');
fprintf(fid,'    } \n');
fprintf(fid,'    slat \n');
fprintf(fid,'    { \n');
fprintf(fid,'        type            noSlip; \n');
fprintf(fid,'    } \n');

fprintf(fid,'} \n');
fprintf(fid,'// ************************************************************************* // \n');
fclose(fid);

%% SCRIVO P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(sprintf('%s/%s/0/p',case_dir,dir),'w+');

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
fprintf(fid,'| =========                 |                                                 |\n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
fprintf(fid,'|  \\    /   O peration     | Version:  5                                     |\n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
fprintf(fid,'|    \\/     M anipulation  |                                                 |\n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
fprintf(fid,'FoamFile\n');
fprintf(fid,'{\n');
fprintf(fid,'    version     2.0;\n');
fprintf(fid,'    format      ascii;\n');
fprintf(fid,'    class       volScalarField;\n');
fprintf(fid,'    object      p;\n');
fprintf(fid,'}\n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');

fprintf(fid,'dimensions      [0 2 -2 0 0 0 0];\n');
fprintf(fid,'internalField   uniform %f;\n',BU_par.p);
fprintf(fid,'boundaryField \n');
fprintf(fid,'\n');
fprintf(fid,'{\n');
fprintf(fid,'    outlet\n');
fprintf(fid,'    {\n');
if strcmp(BU_type,'freestream') == 1
  fprintf(fid,'        type            freestreamPressure;\n');
elseif strcmp(BU_type,'fixedValue') == 1
  fprintf(fid,'        type            fixedValue;\n');
  fprintf(fid,'        value           uniform %f;',BU_par.p);
else
  error('BU_type errato')
end

fprintf(fid,'      \n');
fprintf(fid,'    }\n');

fprintf(fid,'    inlet\n');
fprintf(fid,'    {\n');
if strcmp(BU_type,'freestream') == 1
  fprintf(fid,'        type            freestreamPressure;\n');
elseif strcmp(BU_type,'fixedValue') == 1
  fprintf(fid,'        type            zeroGradient;\n');
else
  error('BU_type errato')
end
fprintf(fid,'       \n');
fprintf(fid,'    }\n');

fprintf(fid,'    front\n');
fprintf(fid,'    {\n');
fprintf(fid,'        type            empty;\n');
fprintf(fid,'    }\n');
fprintf(fid,'    back\n');
fprintf(fid,'    {\n');
fprintf(fid,'        type            empty;\n');
fprintf(fid,'    }\n');
fprintf(fid,'    airfoil\n');
fprintf(fid,'    {\n');
fprintf(fid,'          type            zeroGradient;\n');
fprintf(fid,'    }\n');
fprintf(fid,'    slat\n');
fprintf(fid,'    {\n');
fprintf(fid,'          type            zeroGradient;\n');
fprintf(fid,'    }\n');
fprintf(fid,'  }\n');

fprintf(fid,'  // ************************************************************************* //\n');

fclose(fid);

done = 1;
end
