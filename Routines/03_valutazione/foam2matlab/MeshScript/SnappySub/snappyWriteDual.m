function [done] = snappyWriteDual(nome_out,suffisso,SNAPPY_par,x_TE,TF)

  fid = fopen('./10snappy/system/snappyHexMeshDict','w+');

  fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
  fprintf(fid,'| =========                 |                                                 | \n');
  fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n');
  fprintf(fid,'|  \\    /   O peration     | Version:  2.2.0                                 | \n');
  fprintf(fid,'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n');
  fprintf(fid,'|    \\/     M anipulation  |                                                 | \n');
  fprintf(fid,'\\*---------------------------------------------------------------------------*/ \n');
  fprintf(fid,'FoamFile \n');
  fprintf(fid,'{ \n');
  fprintf(fid,'    version     2.2; \n');
  fprintf(fid,'    format      ascii; \n');
  fprintf(fid,'    class       dictionary; \n');
  fprintf(fid,'    object      snappyHexMeshDict; \n');
  fprintf(fid,'} \n');
  fprintf(fid,' \n');
  fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
  fprintf(fid,' \n');
  fprintf(fid,'// Which of the steps to run \n');
  fprintf(fid,'castellatedMesh %s; \n',TF{1});
  fprintf(fid,'snap            %s; \n',TF{2});
  fprintf(fid,'addLayers       %s; \n',TF{3});
  % // Geometry. Definition of all surfaces. All surfaces are of class
  % // searchableSurface.
  % // Surfaces are used
  % // - to specify refinement for any mesh cell intersecting it
  % // - to specify refinement for any mesh cell inside/outside/near
  % // - to 'snap' the mesh boundary to the surface
  fprintf(fid,'geometry \n');
  fprintf(fid,'{ \n');
  fprintf(fid,'    %s_airfoil.stl \n',nome_out);
  fprintf(fid,'    { \n');
  fprintf(fid,'        type triSurfaceMesh; \n');
  fprintf(fid,'        name airfoil; \n');
  fprintf(fid,'        patchInfo \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            type wall; \n');
  fprintf(fid,'        } \n');
  fprintf(fid,'    } \n');
  fprintf(fid,'    %s.stl \n',strcat(nome_out,suffisso));
  fprintf(fid,'    { \n');
  fprintf(fid,'        type triSurfaceMesh; \n');
  fprintf(fid,'        name slat; \n');
  fprintf(fid,'        patchInfo \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            type wall; \n');
  fprintf(fid,'        } \n');
  fprintf(fid,'    } \n');

  if SNAPPY_par.wakeBox == 1
      fprintf(fid,'    wakebox            //USER DEFINED REGION NAME \n');
      fprintf(fid,'    { \n');
      fprintf(fid,'        type searchableCylinder; \n');
      fprintf(fid,'        point1 (0.25 0 0);    // Height \n');
      fprintf(fid,'        point2 (%2.3f 0 %2.3f);    // Vector \n',SNAPPY_par.xWake,SNAPPY_par.zWake);
      fprintf(fid,'        radius 0.25; \n');
      fprintf(fid,'} \n');
  end

  if SNAPPY_par.LEBox == 1
      fprintf(fid,'    LEBox            //USER DEFINED REGION NAME \n');
      fprintf(fid,'    { \n');
      fprintf(fid,'        type searchableCylinder; \n');
      fprintf(fid,'        point1 (0 -1 0);    // Height \n');
      fprintf(fid,'        point2 (0 1 0);    // Vector \n');
      fprintf(fid,'        radius 0.05; \n');
      fprintf(fid,'} \n');
  end

  if SNAPPY_par.TEBox == 1
      fprintf(fid,'    TEBox            //USER DEFINED REGION NAME \n');
      fprintf(fid,'    { \n');
      fprintf(fid,'        type searchableCylinder; \n');
      fprintf(fid,'        point1 (%f -1 0);    // Height \n',x_TE);
      fprintf(fid,'        point2 (%f 1 0);    // Vector \n',x_TE);
      fprintf(fid,'        radius 0.025; \n');
      fprintf(fid,'} \n');
  end
  fprintf(fid,'}; \n');
  %// Settings for the castellatedMesh generation.
  fprintf(fid,'castellatedMeshControls \n');
  fprintf(fid,'{ \n');
  %     // Refinement parameters
  %     // ~~~~~~~~~~~~~~~~~~~~~
  %     // If local number of cells is >= maxLocalCells on any processor
  %     // switches from from refinement followed by balancing
  %     // (current method) to (weighted) balancing before refinement.
  fprintf(fid,'    maxLocalCells 1000000; \n');
  %     // Overall cell limit (approximately). Refinement will stop immediately
  %     // upon reaching this number so a refinement level might not complete.
  %     // Note that this is the number of cells before removing the part which
  %     // is not 'visible' from the keepPoint. The final number of cells might
  %     // actually be a lot less.
  fprintf(fid,'    maxGlobalCells 10000000; \n');
  %
  %     // The surface refinement loop might spend lots of iterations refining just a
  %     // few cells. This setting will cause refinement to stop if <= minimumRefine
  %     // are selected for refinement. Note: it will at least do one iteration
  %     // (unless the number of cells to refine is 0)
  fprintf(fid,'    minRefinementCells 0; \n');
  %
  %     // Allow a certain level of imbalance during refining
  %     // (since balancing is quite expensive)
  %     // Expressed as fraction of perfect balance (= overall number of cells /
  %     // nProcs). 0=balance always.
  fprintf(fid,'    maxLoadUnbalance 0.10; \n');
  %
  %     // Number of buffer layers between different levels.
  %     // 1 means normal 2:1 refinement restriction, larger means slower
  %     // refinement.
  fprintf(fid,'    nCellsBetweenLevels 12; \n');
  %     // Explicit feature edge refinement
  %     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  %     // Specifies a level for any cell intersected by its edges.
  %     // This is a featureEdgeMesh, read from constant/triSurface for now.
  fprintf(fid,'    features \n');
  fprintf(fid,'    ( \n');
  fprintf(fid,'    ); \n');
  %     // Surface based refinement
  %     // ~~~~~~~~~~~~~~~~~~~~~~~~
  %     // Specifies two levels for every surface. The first is the minimum level,
  %     // every cell intersecting a surface gets refined up to the minimum level.
  %     // The second level is the maximum level. Cells that 'see' multiple
  %     // intersections where the intersections make an
  %     // angle > resolveFeatureAngle get refined up to the maximum level.
  fprintf(fid,'    refinementSurfaces \n');
  fprintf(fid,'    { \n');
  fprintf(fid,'       airfoil \n');
  fprintf(fid,'        { \n');
  %            // Surface-wise min and max refinement level
  fprintf(fid,'            level (%d %d); \n',SNAPPY_par.MinrefFactor,SNAPPY_par.MaxrefFactor);
  fprintf(fid,'        } \n');
    fprintf(fid,'       slat \n');
  fprintf(fid,'        { \n');
  %            // Surface-wise min and max refinement level
  fprintf(fid,'            level (%d %d); \n',SNAPPY_par.MinrefFactor,SNAPPY_par.MaxrefFactor);
  fprintf(fid,'        } \n');
  fprintf(fid,'    } \n');
  fprintf(fid,'    resolveFeatureAngle 50; \n');
  %     // Region-wise refinement
  %     // ~~~~~~~~~~~~~~~~~~~~~~
  %     // Specifies refinement level for cells in relation to a surface. One of
  %     // three modes
  %     // - distance. 'levels' specifies per distance to the surface the
  %     //   wanted refinement level. The distances need to be specified in
  %     //   descending order.
  %     // - inside. 'levels' is only one entry and only the level is used. All
  %     //   cells inside the surface get refined up to the level. The surface
  %     //   needs to be closed for this to be possible.
  %     // - outside. Same but cells outside.
  fprintf(fid,'    refinementRegions \n');
  fprintf(fid,'    {  \n');
  fprintf(fid,'	     airfoil \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            mode distance; \n');
  fprintf(fid,'            levels  ((%f %d) (%f %d)) ; \n',SNAPPY_par.refDist/4,SNAPPY_par.MaxrefFactor+1,SNAPPY_par.refDist,SNAPPY_par.MaxrefFactor);% // levels must be ordered nearest first
  fprintf(fid,'        } \n');
  fprintf(fid,'	     slat \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            mode distance; \n');
  fprintf(fid,'            levels  ((%f %d) (%f %d)) ; \n',SNAPPY_par.refDist/4,SNAPPY_par.MaxrefFactor+1,SNAPPY_par.refDist,SNAPPY_par.MaxrefFactor);% // levels must be ordered nearest first
  
  %fprintf(fid,'            levels  ((%f %d)) ; \n',SNAPPY_par.refDist,SNAPPY_par.MaxrefFactor);% // levels must be ordered nearest first
  fprintf(fid,'        } \n');
  if SNAPPY_par.wakeBox == 1
      fprintf(fid,'	     wakebox \n');
      fprintf(fid,'	     { \n');
      fprintf(fid,'	         mode inside; \n');
      fprintf(fid,'	         levels ((1.0 %d)); \n',SNAPPY_par.MinrefFactor-SNAPPY_par.deltabox);
      fprintf(fid,'	     } \n');
  end


  if SNAPPY_par.LEBox == 1
      fprintf(fid,'	     LEBox \n');
      fprintf(fid,'	     { \n');
      fprintf(fid,'	         mode inside; \n');
      fprintf(fid,'	         levels ((1.0 %d)); \n',SNAPPY_par.MaxrefFactor+SNAPPY_par.deltabox);
      fprintf(fid,'	     } \n');
  end


  if SNAPPY_par.TEBox == 1
      fprintf(fid,'	     TEBox \n');
      fprintf(fid,'	     { \n');
      fprintf(fid,'	         mode inside; \n');
      fprintf(fid,'	         levels ((1.0 %d)); \n',SNAPPY_par.MaxrefFactor+SNAPPY_par.deltabox);
      fprintf(fid,'	     } \n');
  end
  fprintf(fid,'    } \n');
  %     // Mesh selection
  %     // ~~~~~~~~~~~~~~
  %     // After refinement patches get added for all refinementSurfaces and
  %     // all cells intersecting the surfaces get put into these patches. The
  %     // section reachable from the locationInMesh is kept.
  %     // NOTE: This point should never be on a face, always inside a cell, even
  %     // after refinement.
  fprintf(fid,'    locationInMesh (-5.1 0 0.0); \n');
  %     // Whether any faceZones (as specified in the refinementSurfaces)
  %     // are only on the boundary of corresponding cellZones or also allow
  %     // free-standing zone faces. Not used if there are %
  %
  % no faceZones.
  fprintf(fid,'    allowFreeStandingZoneFaces false; \n');
  fprintf(fid,'} \n');
  % // Settings for the snapping.
  fprintf(fid,'snapControls \n');
  fprintf(fid,'{ \n');
  %     //- Number of patch smoothing iterations before finding correspondence
  %     //  to surface
  fprintf(fid,'     nSmoothPatch 3; \n');
  %     //- Relative distance for points to be attracted by surface feature point
  %     //  or edge. True distance is this factor times local
  %     //  maximum edge length.
  %     //    tolerance 4.0;
  fprintf(fid,'    tolerance 4.0; \n');
  %     //- Number of mesh displacement relaxation iterations.
  fprintf(fid,'    nSolveIter 50; \n');
  %     //- Maximum number of snapping relaxation iterations. Should stop
  %     //  before upon reaching a correct mesh.
  fprintf(fid,'    nRelaxIter 5; \n');
  %    // Feature snapping
  %         //- Number of feature edge snapping iterations.
  %         //  Leave out altogether to disable.
  fprintf(fid,'         nFeatureSnapIter 3; \n');
  %         //- Detect (geometric only) features by sampling the surface
  %         //  (default=false).
  fprintf(fid,'        implicitFeatureSnap false; \n');
  %         //- Use castellatedMeshControls::features (default = true)
  fprintf(fid,'        explicitFeatureSnap true; \n');
  %         //- Detect points on multiple surfaces (only for explicitFeatureSnap)
  fprintf(fid,'        multiRegionFeatureSnap true; \n');
  fprintf(fid,'} \n');
  % // Settings for the layer addition.
  fprintf(fid,'addLayersControls \n');
  fprintf(fid,'{ \n');
  %     // Are the thickness parameters below relative to the undistorted
  %     // size of the refined cell outside layer (true) or absolute sizes (false).
  fprintf(fid,'    relativeSizes true; \n');
  %     // Per final patch (so not geometry!) the layer information
  fprintf(fid,'    layers \n');
  fprintf(fid,'    { \n');
  fprintf(fid,'        airfoil \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            nSurfaceLayers %d; \n',SNAPPY_par.imin);
  fprintf(fid,'        } \n');
  fprintf(fid,'        slat \n');
  fprintf(fid,'        { \n');
  fprintf(fid,'            nSurfaceLayers %d; \n',SNAPPY_par.imin);
  fprintf(fid,'        } \n');
  fprintf(fid,'    } \n');
  %     // Expansion factor for layer mesh
  fprintf(fid,'    expansionRatio %3.2f; \n',SNAPPY_par.expRatio);
  %     //- Wanted thickness of final added cell layer. If multiple layers
  %     //  is the thickness of the layer furthest away from the wall.
  %     //  See relativeSizes parameter.
  fprintf(fid,'    finalLayerThickness %3.2f; \n',SNAPPY_par.endLay);
  %     //- Minimum thickness of cell layer. If for any reason layer
  %     //  cannot be above minThickness do not add layer.
  %     //  See relativeSizes parameter.
  fprintf(fid,'    minThickness %3.2f; \n',SNAPPY_par.startLay);
  %     //- If points get not extruded do nGrow layers of connected faces that are
  %     //  also not grown. This helps convergence of the layer addition process
  %     //  close to features.
  %     // Note: changed(corrected) w.r.t 17x! (didn't do anything in 17x)
  fprintf(fid,'    nGrow 0; \n');
  %     // Advanced settings
  %
  %     // No. of steps walking away from the surface
  %     // nMedialAxisIter 10;  // default: 2^31 (unlimited)
  %
  %     //- When not to extrude surface. 0 is flat surface, 90 is when two faces
  %     //  make straight angle.
  fprintf(fid,'    featureAngle 120;  \n');
  %     //- Maximum number of snapping relaxation iterations. Should stop
  %     //  before upon reaching a correct mesh.
  fprintf(fid,'    nRelaxIter 5; \n');
  %     // Number of smoothing iterations of surface normals
  fprintf(fid,'    nSmoothSurfaceNormals 1; \n');
  %     // Number of smoothing iterations of interior mesh movement direction
  fprintf(fid,'    nSmoothNormals 3; \n');
  %     // Smooth layer thickness over surface patches
  fprintf(fid,'    nSmoothThickness 10; \n');
  %     // Stop layer growth on highly warped cells
  fprintf(fid,'    maxFaceThicknessRatio 0.5; \n');
  %     // Reduce layer growth where ratio thickness to medial
  %     // distance is large
  fprintf(fid,'    maxThicknessToMedialRatio 0.2; \n');
  %     // Angle used to pick up medial axis points
  %     // Note: changed(corrected) w.r.t 17x! 90 degrees corresponds to 130 in 17x.
  fprintf(fid,'    minMedianAxisAngle 90; \n');
  %     // Create buffer region for new layer terminations
  fprintf(fid,'    nBufferCellsNoExtrude 0; \n');
  %     // Overall max number of layer addition iterations. The mesher will exit
  %     // if it reaches this number of iterations; possibly with an illegal
  %     // mesh.
  fprintf(fid,'    nLayerIter 50; \n');
  %     // Max number of iterations after which relaxed meshQuality controls
  %     // get used. Up to nRelaxIter it uses the settings in meshQualityControls,
  %     // after nRelaxIter it uses the values in meshQualityControls::relaxed.
  fprintf(fid,'    nRelaxedIter 20; \n');
  fprintf(fid,'} \n');
  % // Generic mesh quality settings. At any undoable phase these determine
  % // where to undo.
  fprintf(fid,'meshQualityControls \n');
  fprintf(fid,'{ \n');
  %     //- Maximum non-orthogonality allowed. Set to 180 to disable.
  fprintf(fid,'    maxNonOrtho 65; \n');%	// NB: mettere valore del checkMesh post castellated e separare 3 fasi
  %     //- Max skewness allowed. Set to <0 to disable.
  fprintf(fid,'    maxBoundarySkewness 4; \n');
  fprintf(fid,'    maxInternalSkewness 2; \n');
  %     //- Max concaveness allowed. Is angle (in degrees) below which concavity
  %     //  is allowed. 0 is straight face, <0 would be convex face.
  %     //  Set to 180 to disable.
  fprintf(fid,'    maxConcave 80; \n');
  %     //- Minimum pyramid volume. Is absolute volume of cell pyramid.
  %     //  Set to a sensible fraction of the smallest cell volume expected.
  %     //  Set to very negative number (e.g. -1E30) to disable.startLay
  fprintf(fid,'    minVol 1e-13; \n');
  %     //- Minimum quality of the tet formed by the face-centre
  %     //  and variable base point minimum decomposition triangles and
  %     //  the cell centre.  Set to very negative number (e.g. -1E30) to
  %     //  disabl \n');e.
  %     //     <0 = inside out tet,
  %     //      0 = flat tet
  %     //      1 = regular tet
  fprintf(fid,'    minTetQuality 1e-30; // ( prima era 1e-30: abilitato ); \n');
  %     //- Minimum face area. Set to <0 to disable.
  fprintf(fid,'    minArea -1; \n');
  %     //- Minimum face twist. Set to <-1 to disable. dot product of face normal
  %     //- and face centre triangles normal
  fprintf(fid,'    minTwist 0.05; \n');%	// 0.05;	// -1; disattivo
  %     //- minimum normalised cell determinant
  %     //- 1 = hex, <= 0 = folded or flattened illegal cell
  fprintf(fid,'    minDeterminant 0.001; \n');
  %     //- minFaceWeight (0 -> 0.5)
  fprintf(fid,'    minFaceWeight 0.05; \n');
  %     //- minVolRatio (0 -> 1)
  fprintf(fid,'    minVolRatio 0.01; \n');
  %     //must be >0 for Fluent compatibility
  fprintf(fid,'    minTriangleTwist -1; \n');
  %     //- if >0 : preserve single systemcells with all points on the surface if the
  %     //  resulting volume after snapping (by approximation) is larger than
  %     //  minVolCollapseRatio times old volume (i.e. not collapsed to flat cell).
  %     //  If <0 : delete always.
  %     //minVolCollapseRatio 0.5;
  %     // Advanced
  %     //- Number of error distribution iterations
  fprintf(fid,'    nSmoothScale 4; \n');
  %     //- amount to scale back displacement at error points
  fprintf(fid,'    errorReduction 0.75; \n');
  %     // Optional : some meshing phases allow usage of relaxed rules.
  %     // See e.g. addLayersControls::nRelaxedIter.
  fprintf(fid,'    relaxed \n');
  fprintf(fid,'    { \n');
  %         //- Maximum non-orthogonality allowed. Set to 180 to disable.
  %         // maxNonOrtho 180;
  fprintf(fid,'    } \n');
  fprintf(fid,'} \n');
  % // Advanced
  % // Flags for optional outputwec
  % // 0 : only write final meshes
  % // 1 : write intermediate meshes
  % // 2 : write volScalarField with cellLevel for postprocessing
  % // 4 : write current intersections as .obj files
  fprintf(fid,'debug 0; \n');
  fprintf(fid,'writeFlags \n');
  fprintf(fid,'( \n');
  fprintf(fid,'    scalarLevels    // write volScalarField with cellLevel for postprocessing \n');
  fprintf(fid,'    layerSets       // write cellSets, faceSets of faces in layer \n');
  fprintf(fid,'    layerFields     // write volScalarField for layer coverage \n');
  fprintf(fid,'); \n');
  % // Merge tolerance. Is fraction of overall bounding box of initial mesh.
  % // Note: the write tolerance needs to be higher than this.
  fprintf(fid,'mergeTolerance 1E-6; \n');
  fprintf(fid,'// ************************************************************************* // \n');
  fclose(fid);
  
  done = 1;
end
