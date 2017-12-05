function [done] = controlWrite(solver,BU_par,CONTROL,case_dir)

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
%% controlDict

fid = fopen(sprintf('%s/%s/system/controlDict',case_dir,dir),'w+');

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
fprintf(fid,'| =========                 |                                                 | \n');
fprintf(fid,'| \\      /  F ield wec        | OpenFOAM: The Open Source CFD Toolbox           | \n');
fprintf(fid,'|  \\    /   O peration     | Version:  5                                     | \n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 | \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/ \n');
fprintf(fid,'FoamFile \n');
fprintf(fid,'{ \n');
fprintf(fid,'   version     2.0; \n');
fprintf(fid,'    format      ascii; \n');
fprintf(fid,'    class       dictionary; \n');
fprintf(fid,'    location    "system"; \n');
fprintf(fid,'    object      controlDict; \n');
fprintf(fid,'} \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');

if strcmp(solver,'simple') == 1
  fprintf(fid,'application     simpleFoam; \n');
elseif strcmp(solver,'piso') == 1
  fprintf(fid,'application     pisoFoam; \n');
else
  error('solver errato')
end

fprintf(fid,'startFrom       latestTime; \n');
fprintf(fid,'startTime       %f; \n',CONTROL.startTime);
fprintf(fid,'stopAt          endTime; \n');
fprintf(fid,'endTime         %1.8f; \n',CONTROL.endTime);
fprintf(fid,'deltaT          %1.8f; \n',CONTROL.deltaT);
fprintf(fid,'writeControl    timeStep; \n');
fprintf(fid,'writeInterval   %d; \n',CONTROL.writeInterval);
fprintf(fid,'purgeWrite      0; \n');
fprintf(fid,'writeFormat     ascii; \n');
fprintf(fid,'writePrecision  6; \n');
fprintf(fid,'writeCompression off; \n');
fprintf(fid,'timeFormat      general; \n');
fprintf(fid,'timePrecision   6; \n');
fprintf(fid,'runTimeModifiable true; \n');

fprintf(fid,'adjustTimeStep on; \n');
fprintf(fid,'maxCo          0.9; \n');

fprintf(fid,'functions \n');
fprintf(fid,'{ \n');
fprintf(fid,' forceCoeffs \n');
fprintf(fid,' { \n');
fprintf(fid,' type forceCoeffs; \n');
fprintf(fid,' functionObjectLibs ("libforces.so"); \n');
fprintf(fid,' writeControl timeStep; \n');
fprintf(fid,' writeInterval 1; \n');
fprintf(fid,' enabled true; \n');
fprintf(fid,' patches \n');
fprintf(fid,'  ( \n');
fprintf(fid,'  airfoil \n');
fprintf(fid,'  ); \n');
fprintf(fid,' rho rhoInf; \n');
fprintf(fid,' log true; \n');
fprintf(fid,' rhoInf %3.4f; \n',BU_par.Rho);
fprintf(fid,' CofR (%3.4f 0 0);\n',0.25*BU_par.L);
fprintf(fid,' liftDir (%3.4f 0 %3.4f);\n',-sin(BU_par.alpha*pi/180),cos(BU_par.alpha*pi/180)); %%%%%%%%%%%%
fprintf(fid,' dragDir (%3.4f 0 %3.4f);\n',cos(BU_par.alpha*pi/180),sin(BU_par.alpha*pi/180));
fprintf(fid,' pitchAxis (0 1 0);\n');
fprintf(fid,' magUInf %3.4f;\n',BU_par.Umag);
fprintf(fid,' lRef %3.4f;\n',BU_par.L);
fprintf(fid,' Aref %3.4f;\n',BU_par.L*BU_par.extrusion_Thickness);
fprintf(fid,' } \n');
fprintf(fid,'} \n');
fprintf(fid,'// ************************************************************************* // \n');

fclose(fid);

done = 1;
end
