function [done] = blockWrite(BLOCKMESH_par,case_dir)

fid = fopen(sprintf('%s/10snappy/system/blockMeshDict',case_dir),'w+');

fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\  \n');
fprintf(fid,'| =========                 |                                                 |  \n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n');
fprintf(fid,'|  \\    /   O peration     | Version:  2.2.0                                 |  \n');
fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |  \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 |  \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/  \n');
fprintf(fid,'FoamFile                                                                         \n');

fprintf(fid,'{\n');
fprintf(fid,'    version     2.0;           \n');
fprintf(fid,'    format      ascii;         \n');
fprintf(fid,'    class       dictionary;    \n');
fprintf(fid,'    object      blockMeshDict; \n');
fprintf(fid,'}                              \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //  \n');
fprintf(fid,'convertToMeters 1; \n');
fprintf(fid,'xMin    %f;\n',BLOCKMESH_par.xMin);
fprintf(fid,'xMax    %f;\n',BLOCKMESH_par.xMax);
fprintf(fid,'zMax    %f;\n',BLOCKMESH_par.zMax);
fprintf(fid,'zMin    %f;\n',BLOCKMESH_par.zMin);
fprintf(fid,'yMin    %f;\n',BLOCKMESH_par.yMin);
fprintf(fid,'yMax    %f;\n',BLOCKMESH_par.yMax);
fprintf(fid,'nx %d;\n',BLOCKMESH_par.nx);
fprintf(fid,'ny %d;\n',BLOCKMESH_par.ny);
fprintf(fid,'nz %d;\n',BLOCKMESH_par.nz);
fprintf(fid,'vertices \n');
fprintf(fid,'(        \n');
fprintf(fid,'    ($xMin $yMin $zMin)\n');% // 0
fprintf(fid,'    ($xMax $yMin $zMin)\n');% // 1
fprintf(fid,'    ($xMax $yMax $zMin)\n');% // 2
fprintf(fid,'    ($xMin $yMax $zMin)\n');% // 3
fprintf(fid,'    ($xMin $yMin $zMax)\n');% // 4
fprintf(fid,'    ($xMax $yMin $zMax)\n');% // 5
fprintf(fid,'    ($xMax $yMax $zMax)\n');% // 6
fprintf(fid,'    ($xMin $yMax $zMax)\n');% // 7
fprintf(fid,');\n');%
fprintf(fid,'edges\n');%
fprintf(fid,'(\n');%
fprintf(fid,');\n');%
fprintf(fid,'blocks\n');%
fprintf(fid,'(\n');%
fprintf(fid,'    hex (0 1 2 3 4 5 6 7) ($nx $ny $nz)\n');%


if BLOCKMESH_par.simplegrading.flag == 1

    fprintf(fid,'    simpleGrading\n');% ( 1 1 1 ) \n');%

    fprintf(fid,'       (\n');
    fprintf(fid,'        (\n');
    fprintf(fid,'            (%f %f %f)\n',BLOCKMESH_par.pd1,BLOCKMESH_par.pn1,BLOCKMESH_par.sg_coeff);
    fprintf(fid,'            (%f %f 1)\n' ,BLOCKMESH_par.pd2,BLOCKMESH_par.pn2);
    fprintf(fid,'            (%f %f %f)\n',BLOCKMESH_par.pd1,BLOCKMESH_par.pn1,1/BLOCKMESH_par.sg_coeff);
    fprintf(fid,'        )\n');
    fprintf(fid,'        1\n');
    fprintf(fid,'        (\n');
    fprintf(fid,'            (%f %f %f)\n',BLOCKMESH_par.pd1,BLOCKMESH_par.pn1,BLOCKMESH_par.sg_coeff);
    fprintf(fid,'            (%f %f 1)\n', BLOCKMESH_par.pd2,BLOCKMESH_par.pn2);
    fprintf(fid,'            (%f %f %f)\n',BLOCKMESH_par.pd1,BLOCKMESH_par.pn1,1/BLOCKMESH_par.sg_coeff);
    fprintf(fid,'        )\n');
    fprintf(fid,'       )\n');

else
fprintf(fid,'    simpleGrading\n');% ( 1 1 1 ) \n');%

    fprintf(fid,'       (\n');
    fprintf(fid,'        1');
    fprintf(fid,'        1\n');
    fprintf(fid,'        1\n');
    fprintf(fid,'       )\n');
   

end

fprintf(fid,');\n');%
fprintf(fid,'boundary\n');%
fprintf(fid,'(\n');%
fprintf(fid,'    outlet\n');%
fprintf(fid,'    {\n');%
fprintf(fid,'        type patch;\n');%
fprintf(fid,'        faces\n');%
fprintf(fid,'        (\n');%
fprintf(fid,'            (2 6 5 1)\n');%
fprintf(fid,'            (4 5 6 7)\n');%
fprintf(fid,'        );\n');%
fprintf(fid,'    }\n');%
fprintf(fid,'    inlet\n');%
fprintf(fid,'    {\n');%
fprintf(fid,'        type patch;\n');%
fprintf(fid,'        faces\n');%
fprintf(fid,'        (\n');%
fprintf(fid,'            (0 4 7 3)\n');%
fprintf(fid,'            (0 1 2 3)\n');%
fprintf(fid,'        );\n');%
fprintf(fid,'    }\n');%
fprintf(fid,'    front\n');%
fprintf(fid,'    {\n');%
fprintf(fid,'        type patch;\n');%
fprintf(fid,'        faces\n');%
fprintf(fid,'        (\n');%
fprintf(fid,'            (7 6 2 3)\n');%
fprintf(fid,'        );\n');%
fprintf(fid,'    }\n');%
fprintf(fid,'    back\n');%
fprintf(fid,'    {\n');%
fprintf(fid,'        type patch;\n');%
fprintf(fid,'        faces\n');%
fprintf(fid,'        (\n');%
fprintf(fid,'            (4 5 1 0)\n');%
fprintf(fid,'        );\n');%
fprintf(fid,'    }\n');%
fprintf(fid,');\n');%
fprintf(fid,'// ************************************************************************* //\n');%
fclose(fid);

done = 1;
end
