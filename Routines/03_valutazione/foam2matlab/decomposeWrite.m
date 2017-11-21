function [done] = decomposeWrite(nProc,case_dir,SOLVER)

if nargin == 1
     case_dir = '.';
end
%% controlDict
%fprintf('%s/30simple/system/decomposeParDict',case_dir);

if strcmp(SOLVER.solver,'simple')
fid = fopen(sprintf('%s/30simple/system/decomposeParDict',case_dir),'w+');
else
fid = fopen(sprintf('%s/40piso/system/decomposeParDict',case_dir),'w+');
end

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\\n');
fprintf(fid,'| =========                 |                                                 |\n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
fprintf(fid,'|  \\    /   O peration     | Version:  2.0.1                                 |\n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |\n');
fprintf(fid,'|    \\/     M anipulation  |                                                 |\n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
fprintf(fid,'FoamFile\n');
fprintf(fid,'{\n');
fprintf(fid,'    version     2.0;\n');
fprintf(fid,'    format      ascii;\n');
fprintf(fid,'    class       dictionary;\n');
fprintf(fid,'    location    "system";\n');
fprintf(fid,'    object      decomposeParDict;\n');
fprintf(fid,'}\n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');

fprintf(fid,'numberOfSubdomains %d;\n',nProc);

fprintf(fid,'method             scotch;\n');

fprintf(fid,'hierarchicalCoeffs\n');
fprintf(fid,'{\n');

if nProc == 1
    disp('ONLY ONE PROCESSOR,ARE YOU SERIUS?');
    fprintf(fid,'    n               ( 1 1 1 );\n');
else
    fprintf(fid,'    n               ( %d 1 %d );\n',floor(nProc/2),2);
end


fprintf(fid,'    delta           0.001;\n');
fprintf(fid,'    order	    xyz;\n');
fprintf(fid,'}\n');

fprintf(fid,'distributed     no;\n');

fprintf(fid,'roots           ( );\n');

fprintf(fid,'// ************************************************************************* //\n');


fclose(fid);
done = 1;

